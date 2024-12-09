! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!-----------------------------------------------------------------------
SUBROUTINE koopmans_ham_proj (delta)
  !---------------------------------------------------------------------
  ! Here the KI hamiltonian is written in terms of projectors on Wannnier 
  ! functions:
  ! \Delta_H_KI = \sum_n (1/2-P_n) * \Delta_n * |w_kn><w_kn| with 
  !   - P_n = \sum_kv f_kv <u_kv|w_kn><w_kn|u_kn> the occupation of the Wannier function n, 
  !   - \Delta_n=\alpha_n <W_n^2 | f_Hxc | W_n^2> the onsite KI correction 
  ! This correction can be applied perturbatively or not (default) on all the KS states
  ! available from the preceeding nscf calculation
  !
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : xk, ngk, igk_k
  USE control_kcw,           ONLY : num_wann, evc0, Hamlt, kcw_iverbosity, &
                                    num_wann_occ, iuwfc_wann, spin_component, nkstot_eff
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx, npw, et, nbnd, current_k
  USE units_lr,              ONLY : iuwfc
  USE wavefunctions,         ONLY : evc
  USE buffers,               ONLY : get_buffer, save_buffer
  USE io_files,              ONLY : nwordwfc
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  USE lsda_mod,              ONLY : lsda, isk, current_spin
  USE uspp,                  ONLY : nkb, vkb
  USE uspp_init,             ONLY : init_us_2
  USE noncollin_module,      ONLY : npol
  !
  IMPLICIT NONE
  !
  ! The k point index 
  INTEGER :: ik, ibnd, ik_pw, k, eig_start, eig_win
  !
  ! The on-site KI correction \alpha_n*<W_0n^2|f_Hxc|W_0n^2> (computed in dH_ki_wann.f90)
  COMPLEX(DP), INTENT(IN) :: delta (num_wann)
  ! and occupations numbers of WFs: P_n = \sum_kv f_kv <u_kv|w_kn><w_kn|u_kn> 
  REAL(DP) :: occ_mat(num_wann)
  !
  COMPLEX(DP) :: overlap
  COMPLEX(DP) :: ham(nbnd,nbnd), eigvc(nbnd,nbnd), deltah(nbnd,nbnd)
  !
  ! The new eigenalues 
  REAL(DP), ALLOCATABLE :: eigvl_aux(:)
  COMPLEX(DP) , ALLOCATABLE :: eigvc_aux(:,:), ham_aux(:,:), evc_aux(:,:)
  !
  ! The new eigenalues
  REAL(DP) :: eigvl(nbnd)
  REAL(DP) :: eigvl_pert(nbnd)
  REAL(DP) :: eigvl_ks(nbnd)
  !
  INTEGER :: i
  ! 
  REAL(DP) :: ehomo, elumo
  REAL(DP) :: ehomo_ks, elumo_ks
  REAL(DP) :: ehomo_pert, elumo_pert
  INTEGER  :: lrwannfc
  REAL(DP), EXTERNAL :: get_clock
  !
  ALLOCATE (evc_aux(npwx*npol,nbnd))
  !
  WRITE( stdout, '(/,5X, "INFO: BUILD and DIAGONALIZE the KI HAMILTONIAN")')
  WRITE( stdout, '(  5X, "INFO: Projectors scheme")')
  !
  ! The occupation matrix
  ! P_n = \sum_kv f_kv <u_kv|w_kn><w_kn|u_kn> 
  CALL occupations(occ_mat)
  WRITE(stdout, 900) get_clock('KCW')
  ! 
#ifdef DEBUG
  WRITE(stdout,'(/,5X,"Screened on-site correction:")')
  WRITE(stdout,'(5X,2(F10.6, 2x), F10.6)') (delta(iwann), occ_mat(iwann), iwann=1, num_wann)
#endif
  !
  ehomo=-1D+6
  elumo=+1D+6
  ehomo_ks=-1D+6
  elumo_ks=+1D+6
  ehomo_pert=-1D+6
  elumo_pert=+1D+6
  !
  DO ik = 1, nkstot_eff
    !
    ik_pw = ik + (spin_component-1)*(nkstot_eff)
    WRITE( stdout, 9020 ) ( xk(i,ik_pw), i = 1, 3 )
    CALL get_buffer ( evc, nwordwfc, iuwfc, ik_pw )
    npw = ngk(ik_pw)
    !
    ehomo_ks = MAX ( ehomo_ks, et(num_wann_occ  , ik_pw) )
    IF (nbnd > num_wann_occ) elumo_ks = MIN ( elumo_ks, et(num_wann_occ+1, ik_pw) )
    !
    !
    ! Build and diagonalize the projector-based KI Hamiltonian on the Hilbert space spanned
    ! by the KS states available from the preceeding nscf calculation
    ! KI contribution at k: deltah_ij = \sum_n [ (1/2-P_n)D_n <u_ki | w_kn><w_kn | u_kj>]
    !
    IF (.FALSE.) THEN
      ! In the canonicla KS basis the KS hamiltonian is already diagonal
      ! maening there is no need to re-build it. I keep this as is in case
      ! we want to use a diffferent basis
      current_k = ik_pw
      IF ( lsda ) current_spin = isk(ik_pw)
      IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik_pw), xk(1,ik_pw), vkb )
      ! FIXME: ned o modify ks_hamiltonian to pass the hammiltonian (and not use Hamlt of kcw_comm)
      CALL ks_hamiltonian (evc, ik_pw, nbnd, .false.)
      !
      ! The KS hamiltonian in the Wannier Gauge (just to check)
      ham(:,:) = Hamlt(ik,:,:) 
      CALL cdiagh( nbnd, ham, nbnd, eigvl, eigvc )
      !
      !WRITE( stdout, 9020 ) ( xk(i,ik_pw), i = 1, 3 )
      ! this shoud perfetcly matchs with the KS eigenvalues. If not, there is a problem
      WRITE( stdout, '(10x, "KS* ",8F11.4)' ) (eigvl(ibnd)*rytoev, ibnd=1,nbnd)
      !
    ELSE
      !
      ! The KS Hamiltonian in the KS basis
      ham(:,:)=CMPLX(0.D0, 0.D0, kind=DP)
      DO i = 1, nbnd 
        ham(i,i)    = et(i,ik_pw)
        eigvl_ks(i) = et(i,ik_pw)
      ENDDO
      !
    ENDIF
    !
    ! \Delta_H_KI_ij = <psi_i | \Delta_H_KI | psi_j> with
    ! \Delta_H_KI = \sum_n (1/2-P_n) * \Delta_n * |w_kn><w_kn|
    CALL dki_hamiltonian (evc, ik, nbnd, occ_mat, delta, deltah) 
    !
    !
    ! Because we have defined a uniq KI Hamiltonian, we can do a perturbative approach
    ! i.e. we keep only the diagonal part of the KI Hamiltoniana
    DO i = 1, nbnd
      eigvl_pert(i) = et(i,ik_pw) + DBLE(deltah(i,i))
    ENDDO
    ehomo_pert = MAX ( ehomo_pert, eigvl_pert(num_wann_occ ) )
    IF (nbnd > num_wann_occ) elumo_pert = MIN ( elumo_pert, eigvl_pert(num_wann_occ+1 ) )
    !
    ! Add the KI contribution to the KS Hamiltonian
    ham(:,:) = ham(:,:) + deltah(:,:) 
    ! And Diagonalize it 
    CALL cdiagh( nbnd, ham, nbnd, eigvl, eigvc )
    !
    IF (kcw_iverbosity .gt. 1) THEN
      eig_win=5
      eig_start = MAX(num_wann_occ-eig_win+1, 1) 
      WRITE(stdout,'( 12x, "KI Eigenvalues around Fermi as a function of the Hilbert space")')
      WRITE(stdout,'( 12x, "KI Eigenvaluse from", I5, " to", I5, " (i.e. +/-", I3, " around VBM),/")') &
              eig_start, num_wann_occ+eig_win, eig_win
      WRITE(stdout,'( 12x, " dim   eig1      eig2      ...")')
      DO k = eig_start, nbnd
         !
         ALLOCATE (ham_aux(k,k), eigvl_aux(k), eigvc_aux(k,k))
         !ham_aux(1:k,1:k) = ham(1:k,1:k)
         ham_aux(1:k,1:k) = ham(1:k,1:k)
         !
         CALL cdiagh( k, ham_aux, k, eigvl_aux, eigvc_aux )
         !
         !neig_max = MIN (20, nbnd)
         IF (k.le.num_wann_occ+eig_win) THEN 
            WRITE(stdout,'(12x, I3, 10F10.4)') k, eigvl_aux(eig_start:k)*rytoev   
         ELSE 
            WRITE(stdout,'(12x, I3, 10F10.4)') k, eigvl_aux(eig_start:num_wann_occ+eig_win)*rytoev 
         ENDIF
         !
         DEALLOCATE (ham_aux)
         DEALLOCATE (eigvl_aux, eigvc_aux)
         !
      ENDDO
      WRITE(stdout,*)
    ENDIF
    !
    !Overwrite et and evc
    et(1:nbnd, ik_pw) = eigvl(1:nbnd)
    ! MB
    ! This is different wrt koopmans_ham.f90:
    ! (1) the first dimension (row of A/evc) = npwx, not npw;
    ! (2) cannot use the same input matrix as output; need evc_aux
    CALL ZGEMM( 'N','N', npwx*npol, nbnd, nbnd, ONE, evc, npwx*npol, eigvc, nbnd, &
    ZERO, evc_aux, npwx*npol )
    evc(:,:) = evc_aux(:,:)
    CALL save_buffer ( evc, nwordwfc, iuwfc, ik_pw )
    !
    ehomo = MAX ( ehomo, eigvl(num_wann_occ ) )
    IF (nbnd > num_wann_occ) elumo = MIN ( elumo, eigvl(num_wann_occ+1 ) )
    !
    WRITE( stdout, '(10x, "KS  ",8F11.4)' ) (eigvl_ks(ibnd)*rytoev, ibnd=1,nbnd)
    WRITE( stdout, '(10x, "KI  ",8F11.4)' ) (eigvl   (ibnd)*rytoev, ibnd=1,nbnd)
    WRITE( stdout, '(10x, "pKI ",8F11.4)' ) (eigvl_pert(ibnd)*rytoev, ibnd=1,nbnd)
    WRITE(stdout, 901) get_clock('KCW')
    !
  ENDDO
  !
  IF ( elumo < 1d+6) THEN
    WRITE( stdout, 9042 ) ehomo_ks*rytoev, elumo_ks*rytoev
    WRITE( stdout, 9044 ) ehomo*rytoev, elumo*rytoev
    WRITE( stdout, 9046 ) ehomo_pert*rytoev, elumo_pert*rytoev
  ELSE
    WRITE( stdout, 9043 ) ehomo_ks*rytoev
    WRITE( stdout, 9045 ) ehomo*rytoev
    WRITE( stdout, 9047 ) ehomo_pert*rytoev
  END IF
  !
9043 FORMAT(/,8x, 'KS  highest occupied level (ev): ',F10.4 )
9042 FORMAT(/,8x, 'KS  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9045 FORMAT(  8x, 'KI  highest occupied level (ev): ',F10.4 )
9044 FORMAT(  8x, 'KI  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9047 FORMAT(  8x, 'pKI highest occupied level (ev): ',F10.4 )
9046 FORMAT(  8x, 'pKI highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
900 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
901 FORMAT('          total cpu time spent up to now is ',F10.1,' secs' )
  !
  DEALLOCATE (evc_aux)
  !
  RETURN
  !
  CONTAINS 
  !
  ! !----------------------------------------------------------------
  SUBROUTINE dki_hamiltonian (evc, ik, h_dim, occ_mat, delta, deltah)
    !----------------------------------------------------------------
    !
    USE buffers,               ONLY : get_buffer
    USE control_kcw,           ONLY : num_wann
    USE wvfct,                 ONLY : npwx
    USE control_flags,         ONLY : gamma_only
    USE gvect,                 ONLY : gstart
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: h_dim
    COMPLEX(DP) :: delta (num_wann)
    COMPLEX(DP), INTENT(IN) :: evc(npwx*npol,h_dim)
    INTEGER, INTENT(IN) :: ik
    REAL(DP), INTENT(IN) :: occ_mat(num_wann)
    COMPLEX(DP), INTENT(OUT) :: deltah(h_dim,h_dim)
    !
    INTEGER :: lrwannfc, ib, jb, iwann
    COMPLEX(DP) :: overlap_in, overlap_nj, overlap
    !
    lrwannfc = num_wann*npwx*npol
    CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
    !
    deltah = CMPLX(0.D0, 0.D0, kind=DP)
    DO ib = 1, nbnd
      DO jb = ib, nbnd
        !
        DO iwann = 1, num_wann
          IF ( gamma_only ) THEN 
            overlap_in = 2.D0 * SUM(DBLE(CONJG(evc(1:npw*npol,ib))*(evc0(1:npw*npol,iwann))))
            overlap_nj = 2.D0 * SUM(DBLE(CONJG(evc0(1:npw*npol,iwann))*(evc(1:npw*npol,jb))))
            IF (gstart == 2) THEN
               overlap_in = overlap_in - DBLE(CONJG(evc(1,ib))*(evc0(1,iwann)))
               overlap_nj = overlap_nj - DBLE(CONJG(evc0(1,iwann))*(evc(1,jb)))
            ENDIF
          ELSE
            overlap_in = SUM(CONJG(evc(1:npw*npol,ib))*(evc0(1:npw*npol,iwann)))
            overlap_nj = SUM(CONJG(evc0(1:npw*npol,iwann))*(evc(1:npw*npol,jb)))
          ENDIF
          CALL mp_sum (overlap_in, intra_bgrp_comm)
          CALL mp_sum (overlap_nj, intra_bgrp_comm)
          overlap = overlap_in*overlap_nj
          deltah(ib,jb) = deltah(ib,jb) + (0.5D0 - occ_mat(iwann)) * delta(iwann) * overlap
          !WRITE(*,'(3X, 2I5, 2F20.12, F20.12, 2F20.12)') ibnd, iwann, delta(iwann), occ_mat(iwann), overlap 
        ENDDO
        IF (ib /= jb) deltah(jb,ib) = CONJG(deltah(ib,jb))
        !
      ENDDO
    ENDDO
    !
    !
  END SUBROUTINE dki_hamiltonian
  !
  ! !----------------------------------------------------------------
  SUBROUTINE occupations (occ_mat)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : num_wann, spin_component
    USE wvfct,                 ONLY : nbnd
    USE lsda_mod,              ONLY : nspin
    USE klist,                 ONLY : ngk
    USE wvfct,                 ONLY : wg
    USE mp_bands,              ONLY : intra_bgrp_comm
    USE mp,                    ONLY : mp_sum
    USE control_flags,         ONLY : gamma_only
    USE gvect,                 ONLY : gstart
    !
    REAL(DP), INTENT(INOUT) :: occ_mat(num_wann)
    INTEGER :: iwann, ik, ibnd, ik_pw, spin_deg
    !
    ! The canonical occupation matrix (fermi dirac or alike)
    occ_mat = REAL(0.D0, kind=DP)
    !
    DO iwann= 1, num_wann
      !
      DO ik= 1, nkstot_eff
        ik_pw = ik + (spin_component-1)*(nkstot_eff)
        npw = ngk(ik_pw)
        lrwannfc = num_wann*npwx*npol
        CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
        CALL get_buffer ( evc, nwordwfc, iuwfc, ik_pw )
        !
        spin_deg=1
        IF(nspin == 1) spin_deg = 2
        overlap = CMPLX(0.D0, 0.D0, kind=DP)
        DO ibnd = 1, nbnd
          IF (gamma_only) THEN 
             overlap = 2.D0*SUM(DBLE(CONJG(evc(1:npw*npol,ibnd))*(evc0(1:npw*npol,iwann))))
             IF (gstart ==2 ) overlap = overlap-1.D0*(DBLE(CONJG(evc(1,ibnd))*(evc0(1,iwann))))
          ELSE
             overlap = SUM(CONJG(evc(1:npw*npol,ibnd))*(evc0(1:npw*npol,iwann)))
          ENDIF
          CALL mp_sum (overlap, intra_bgrp_comm)
          overlap = CONJG(overlap)*overlap
          occ_mat(iwann) = occ_mat(iwann) + wg(ibnd,ik)/spin_deg * REAL(overlap)
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE occupations

END SUBROUTINE koopmans_ham_proj
