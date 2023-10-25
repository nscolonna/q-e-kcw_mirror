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
SUBROUTINE koopmans_ham_proj ()
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
  USE klist,                 ONLY : nkstot, xk, ngk, igk_k
  USE control_kcw,           ONLY : num_wann, l_alpha_corr, evc0, hamlt, l_diag, &
                                    alpha_final, num_wann_occ, iuwfc_wann, spin_component
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx, npw, et, nbnd, current_k
  USE units_lr,              ONLY : iuwfc
  USE wavefunctions,         ONLY : evc
  USE buffers,               ONLY : get_buffer, save_buffer
  USE io_files,              ONLY : nwordwfc
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  USE lsda_mod,              ONLY : lsda, isk, current_spin, nspin
  USE uspp,                  ONLY : nkb, vkb
  USE uspp_init,             ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  ! The k point index 
  INTEGER :: ik, ibnd, ik_pw
  !
  ! The on-site KI correction <W_0n^2|f_Hxc|W_0n^2>, and occupations numbers of WFs
  COMPLEX(DP) :: delta (num_wann)
  REAL(DP) :: occ_mat(num_wann)
  !
  COMPLEX(DP) :: delta_k, overlap
  COMPLEX(DP) :: ham(nbnd,nbnd), eigvc(npwx,nbnd), deltah(nbnd,nbnd)
  !
  ! The new eigenalues 
  REAL(DP) :: et_ki(nbnd,nkstot)
  !
  ! The correction to the diagonal term beyond second order
  REAL(DP) :: ddH(num_wann)
  !
  ! The new eigenalues
  REAL(DP) :: eigvl(nbnd)
  !
  INTEGER :: i, iwann
  ! 
  REAL(DP) :: ehomo, elumo
  REAL(DP) :: ehomo_ks, elumo_ks
  INTEGER  :: lrwannfc
  !
  WRITE(stdout, '(/,5X, "INFO: KI Hamiltonian using Projectors")')
  !
  ! The scalar term R=0 i=j 
  delta=CMPLX(0.D0,0.D0,kind=DP)
  CALL ham_scalar (delta)
  CALL occupations(occ_mat)
  ! 
#ifdef DEBUG
  WRITE( stdout,'(/,5X," "Bare on-site correction:")')
  WRITE(stdout,'(5X,10(2F10.6, 2x))') (delta(iwann), iwann=1, num_wann)
#endif
  !
  ! ... The correction beyond 2nd order 
  IF (l_alpha_corr) CALL beyond_2nd (delta, ddH)
  !
  ! Apply alpha to get the screened on-site KI correction
  delta(1:num_wann) = alpha_final(1:num_wann) * delta(1:num_wann)
  !
#ifdef DEBUG
  WRITE( stdout,'(/,5X," Screened on-site correction:")')
  WRITE(stdout,'(5X,10(2F10.6, 2x))') (delta(iwann), iwann=1, num_wann)
#endif
  !
  IF (.NOT. l_diag) WRITE(stdout, '(/, 5x, "INFO: FULL KI")')
  IF (l_diag)       WRITE(stdout, '(/, 5x, "INFO: PERTURBATIVE KI")') 
  !
  ehomo=-1D+6
  elumo=+1D+6
  ehomo_ks=-1D+6
  elumo_ks=+1D+6
  !
  DO ik = 1, nkstot/nspin
    !
    ik_pw = ik + (spin_component-1)*(nkstot/nspin)
    CALL get_buffer ( evc, nwordwfc, iuwfc, ik_pw )
    npw = ngk(ik_pw)
    !
    ehomo_ks = MAX ( ehomo_ks, et(num_wann_occ  , ik_pw) )
    elumo_ks = MIN ( elumo_ks, et(num_wann_occ+1, ik_pw) )
    !
    WRITE( stdout, 9020 ) ( xk(i,ik_pw), i = 1, 3 )
    WRITE( stdout, '(10x, "KS  ",8F11.4)' ) (et(ibnd,ik_pw)*rytoev, ibnd=1,nbnd)
    !
    IF ( .NOT. l_diag ) THEN 
      ! Build and diagonalize the projector-based KI Hamiltonian on the Hilbert space spanned
      ! by the KS states available from the preceeding nscf calculation
      ! KI contribution at k: deltah_ij = \sum_n [ (1/2-P_n)D_n <u_ki | w_kn><w_kn | u_kj>]
      !
      IF (.FALSE.) THEN
        ! In the canonicla KS basis the KS hamiltonian is already diuagonal
        ! maening there is no need to re-build it. I keep this as is in case
        ! we want to use a diffferent basis
        current_k = ik_pw
        IF ( lsda ) current_spin = isk(ik_pw)
        IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik_pw), xk(1,ik_pw), vkb )
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
        Hamlt(ik,:,:)=CMPLX(0.D0, 0.D0, kind=DP)
        DO i = 1, nbnd 
          Hamlt(ik,i,i)=et(i,ik_pw)
        ENDDO
        !
      ENDIF
      !
      CALL dki_hamiltonian (evc, ik, nbnd, occ_mat, delta, deltah) 
      !
      ! Add to the KS Hamiltonian
      Hamlt(ik,:,:) = Hamlt(ik,:,:) + deltah(:,:) 
      !
      ham(:,:) = Hamlt(ik,:,:) 
      CALL cdiagh( nbnd, ham, nbnd, eigvl, eigvc )
      !
      et_ki(1:nbnd,ik)=eigvl(1:nbnd)
      !
    ELSE ! Perturbative, corrects only eigenvalues. 
      !
      ! Build a diagonal correction using the projector-based KI Hamiltonian 
      ! Correction at k: Dk = \sum_n [ (1/2-P_iw)D_ii <u_kv | u_kn><u_kn | uk_v>]
      !
      lrwannfc = num_wann*npwx
      CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
      !
      DO ibnd = 1, nbnd
        delta_k = CMPLx(0.D0, 0.D0, kind=DP)
        DO iwann = 1, num_wann
          overlap = SUM(CONJG(evc(1:npw,ibnd))*(evc0(1:npw,iwann)))
          CALL mp_sum (overlap, intra_bgrp_comm)
          overlap = CONJG(overlap)*overlap
          delta_k = delta_k + (0.5D0 - occ_mat(iwann)) * delta(iwann) * overlap
          !WRITE(*,'(3X, 2I5, 2F20.12, F20.12, 2F20.12)') ibnd, iwann, delta(iwann), occ_mat(iwann), overlap 
        ENDDO
        et_ki(ibnd,ik) = et(ibnd,ik_pw) + REAL(delta_k)
        !WRITE (*,'(3I5, 2X, 4F20.12)') ik, ik_pw, ibnd, et(ibnd,ik_pw)*rytoev, et_ki(ibnd,ik)*rytoev, delta_k*rytoev
      ENDDO
      !
      !
    ENDIF
    !
    ehomo = MAX ( ehomo, et_ki(num_wann_occ, ik ) )
    IF (nbnd > num_wann_occ) elumo = MIN ( elumo, et_ki(num_wann_occ+1, ik ) )
    !
    WRITE( stdout, '(10x, "KI  ",8F11.4)' ) (et_ki(ibnd,ik)*rytoev, ibnd=1,nbnd)
    !
  ENDDO
  !
  IF ( elumo < 1d+6) THEN
     WRITE( stdout, 9042 ) ehomo_ks*rytoev, elumo_ks*rytoev
     WRITE( stdout, 9044 ) ehomo*rytoev, elumo*rytoev
  ELSE
     WRITE( stdout, 9043 ) ehomo_ks*rytoev
     WRITE( stdout, 9045 ) ehomo*rytoev
  END IF
  !
  !Overwrite et
  DO ik = 1, nkstot/nspin
    ik_pw = ik + (spin_component-1)*(nkstot/nspin)
    et(1:nbnd, ik_pw) = et_ki(1:nbnd,ik)
  ENDDO
  ! 
9043 FORMAT(/,8x, 'KS       highest occupied level (ev): ',F10.4 )
9042 FORMAT(/,8x, 'KS       highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9045 FORMAT(  8x, 'KI[2nd]  highest occupied level (ev): ',F10.4 )
9044 FORMAT(  8x, 'KI[2nd]  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
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
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: h_dim
    COMPLEX(DP) :: delta (num_wann)
    COMPLEX(DP), INTENT(IN) :: evc(npwx,h_dim)
    INTEGER, INTENT(IN) :: ik
    REAL(DP), INTENT(IN) :: occ_mat(num_wann)
    COMPLEX(DP), INTENT(OUT) :: deltah(h_dim,h_dim)
    !
    INTEGER :: lrwannfc, ib, jb, iwann
    COMPLEX(DP) :: overlap_in, overlap_nj, overlap
    !
    lrwannfc = num_wann*npwx
    CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
    !
    deltah = CMPLX(0.D0, 0.D0, kind=DP)
    DO ib = 1, nbnd
      DO jb = ib, nbnd
        !
        DO iwann = 1, num_wann
          overlap_in = SUM(CONJG(evc(1:npw,ib))*(evc0(1:npw,iwann)))
          overlap_nj = SUM(CONJG(evc0(1:npw,iwann))*(evc(1:npw,jb)))
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
  ! !----------------------------------------------------------------
  SUBROUTINE occupations (occ_mat)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : num_wann, spin_component
    USE wvfct,                 ONLY : nbnd
    USE lsda_mod,              ONLY : nspin
    USE klist,                 ONLY : nkstot, ngk
    USE wvfct,                 ONLY : wg
    USE mp_bands,              ONLY : intra_bgrp_comm
    USE mp,                    ONLY : mp_sum
    !
    REAL(DP), INTENT(INOUT) :: occ_mat(num_wann)
    INTEGER :: iwann, ik, ibnd, ik_pw
    !
    ! The canonical occupation matrix (fermi dirac or alike)
    occ_mat = REAL(0.D0, kind=DP)
    !
    DO iwann= 1, num_wann
      !
      DO ik= 1, nkstot/nspin
        ik_pw = ik + (spin_component-1)*(nkstot/nspin)
        npw = ngk(ik_pw)
        lrwannfc = num_wann*npwx
        CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
        CALL get_buffer ( evc, nwordwfc, iuwfc, ik_pw )
        !
        overlap = CMPLX(0.D0, 0.D0, kind=DP)
        DO ibnd = 1, nbnd
          overlap = SUM(CONJG(evc(1:npw,ibnd))*(evc0(1:npw,iwann)))
          CALL mp_sum (overlap, intra_bgrp_comm)
          overlap = CONJG(overlap)*overlap
          occ_mat(iwann) = occ_mat(iwann) + wg(ibnd,ik) * REAL(overlap)
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE occupations

  !
  !----------------------------------------------------------------
  SUBROUTINE beyond_2nd (deltaH, ddH)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : num_wann, alpha_final, alpha_final_full
    USE control_lr,            ONLY : lrpa
    !
    IMPLICIT NONE  
    !
    COMPLEX(DP), INTENT(INOUT) :: deltaH(num_wann)
    REAL(DP), INTENT(OUT) :: ddH(num_wann)
    !
    REAL(DP) :: alpha_(num_wann), alpha_fd
    ! weight of the q points, alpha from LR, alpha after the all-order correction. 
    !
    REAL(DP) second_der(num_wann), delta 
    !
    INTEGER :: iwann
    !
    ddH = 0.D0
    !
    WRITE( stdout, '(/,5X, "INFO: Correction beyond 2nd order ...",/)')
    IF (lrpa) THEN 
      !
      WRITE(*, '(8X,"INFO: l_alpha_corr and lrpa are NOT consistent.At RPA")')
      WRITE(*, '(8X,"      level there is no contribution beyond 2nd order.")')
      WRITE(*, '(8X,"      Nothing to do here. RETURN")')
      !
      RETURN
      !
    ENDIF 
    !
    DO iwann = 1, num_wann
      second_der(iwann) = -REAL(deltaH(iwann))
    ENDDO
    !
    !DO iwann = 1, num_wann_occ
    DO iwann = 1, num_wann
      delta =0.D0 
      alpha_(iwann) = alpha_final(iwann) 
      !
      !
       ! ... Compute the difference between the parabolic extrapolation at N \pm 1 and the real 
       ! ... value of the energy in the frozen orbital approximation ...
       CALL alpha_corr (iwann, delta)
       ddH(iwann) = delta
       !deltaH(iwann,iwann) = deltaH(iwann,iwann)-ddH(iwann)
       !
       ! ... The new alpha that should be closer to the Finite-difference one ...
       ! ... Remember DeltaH is nothing but the second derivative wrt the orbital occupation ...
       alpha_fd = (alpha_final(iwann)*second_der(iwann) + delta)/ (second_der(iwann)+delta)
       IF(nkstot/nspin == 1) alpha_final_full(iwann) = alpha_fd
       !
       ! ... Since the ham in the frozen approximation is approximated to second
       ! ... order only, this is the alpha we want to use. Only the
       ! ... numerator matters.
       alpha_(iwann) = (alpha_final(iwann)*second_der(iwann) + delta)/second_der(iwann)
       !
       ! ... Write it just to compare with the FD one from CP ... 
       WRITE(stdout,'(5X, "INFO: iwann, LR-alpha, FD-alpha, alpha", i3, 3f12.8)') iwann, alpha_final(iwann),alpha_fd,  alpha_(iwann)
       !
       !WRITE(stdout,'("Nicola", i3, 6f12.8)') iwann, deltaH(iwann,iwann)
       !
       ! Re-define the corrected screening parameter. 
       alpha_final(iwann) = alpha_(iwann) 
       WRITE( stdout, '(5X,"INFO: alpha RE-DEFINED ...", i5, f12.8)') iwann, alpha_final(iwann)
      !
    ENDDO
    !
  END SUBROUTINE beyond_2nd
  !
  !----------------------------------------------------------------
  SUBROUTINE ham_scalar (delta)
    !----------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    USE control_kcw,          ONLY : num_wann, nqstot, iurho_wann, &
                                     spin_component
    USE fft_base,             ONLY : dffts
    USE cell_base,            ONLY : omega
    USE gvecs,                ONLY : ngms
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE mp,                   ONLY : mp_sum
    USE buffers,              ONLY : get_buffer
    !
    IMPLICIT NONE
    !
    ! The scalar contribution to the hamiltonian 
    COMPLEX(DP), INTENT(INOUT) :: delta (num_wann)
    !
    ! Couters for the q point, wannier index. record length for the wannier density
    INTEGER :: iq, iwann, lrrho
    !
    ! The periodic part of the wannier orbital density
    COMPLEX(DP) :: rhowann(dffts%nnr, num_wann), rhor(dffts%nnr), delta_vr(dffts%nnr,nspin), delta_vr_(dffts%nnr,nspin)
    !
    ! The self Hartree
    COMPLEX(DP) :: sh(num_wann)
    !
    ! Auxiliary variables 
    COMPLEX(DP), ALLOCATABLE  :: rhog(:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
    !
    ! The weight of each q point
    REAL(DP) :: weight(nqstot)
    !
    WRITE( stdout, '(/,5X, "INFO: KC SCALAR TERM CALCULATION ... START")')
    !
    ALLOCATE ( rhog (ngms) , delta_vg(ngms,nspin), vh_rhog(ngms), delta_vg_(ngms,nspin) )
    !
    DO iq = 1, nqstot
      !
      lrrho=num_wann*dffts%nnr
      CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
      !! Retrive the rho_wann_q(r) from buffer in REAL space
      !
      weight(iq) = 1.D0/nqstot ! No SYMM 
      !
      DO iwann = 1, num_wann  ! for each band, that is actually the perturbation
         !
         rhog(:)         = CMPLX(0.D0,0.D0,kind=DP)
         delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
         vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
         rhor(:)         = CMPLX(0.D0,0.D0,kind=DP)
         !
         rhor(:) = rhowann(:,iwann)
         !! The periodic part of the orbital desity in real space
         !
         CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
         !! The periodic part of the perturbation DeltaV_q(G)
         ! 
         sh(iwann) = sh(iwann) + 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:)                )*weight(iq)*omega
         delta(iwann) = delta(iwann) + sum (CONJG(rhog (:)) * delta_vg(:,spin_component)) &
                                     * weight(iq) * omega
         !
      ENDDO
      ! 
    ENDDO ! qpoints
    WRITE( stdout, '(/,5X, "INFO: KC SCALAR TERM CALCULATION ... END")')
    !
    DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
    !
    CALL mp_sum (delta, intra_bgrp_comm)
    CALL mp_sum (sh, intra_bgrp_comm)
   !
  END SUBROUTINE ham_scalar
  ! 
END SUBROUTINE koopmans_ham_proj
