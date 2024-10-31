!
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
SUBROUTINE koopmans_ham_uniq ( dH_wann )
  !---------------------------------------------------------------------
  !
  ! Here the KI hamiltonian is written in terms of projectors on Wannnier
  ! functions:
  ! \Delta H_KI = \sum_nm |w_n> \Delta H_nm < w_m| where \Delta H_nm = <w_n | h_m |w_m> 
  ! computed in the standard way (see dH_ki_wann.f90). 
  ! Then we build and Diagonalize the KI hamiltoinan H_KS+\Delta H_KI on the basis of the 
  ! KS orbitals from the NSCF calculation: espilon_i = Diag [ < \phi_i | H_KS + \Delta H_KI | phi_j > ] 
  ! < \phi_i | H_KS | phi_j > = \delta_ij \epsilon_i^KS 
  ! < \phi_i | \Delta H_KI | phi_j > = \sum_nm <phi_i|w_n> \Delta H_nm <w_m|phi_j> 
  !
  ! NB: In principle one can use any other basis or iterative digonalization technique.  
  
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : xk, ngk
  USE control_kcw,           ONLY : num_wann, evc0, spin_component, &
                                    num_wann_occ, iuwfc_wann, nkstot_eff, &
                                    kcw_iverbosity
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx, npw, et, nbnd
  USE units_lr,              ONLY : iuwfc
  USE wavefunctions,         ONLY : evc
  USE buffers,               ONLY : get_buffer, save_buffer
  !
  USE io_files,              ONLY : nwordwfc
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  USE noncollin_module,      ONLY : npol
  !
  IMPLICIT NONE
  !
  ! The k point index 
  INTEGER :: ik, ik_pw
  !
  ! the KI hamiltonian on the Wannier basis <w_i|dh_j|w_j> 
  COMPLEX(DP), INTENT (IN) :: dH_wann(nkstot_eff,num_wann,num_wann)
  COMPLEX(DP), ALLOCATABLE :: dH_wann_aux(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_aux(:,:)
  ! 
  ! the KI operator on the KS basis of the NSCF calculation
  COMPLEX(DP) :: deltah(nbnd,nbnd)
  ! the new hamitonain, and the new eigenvalues and eigenvectors at a given k-point
  COMPLEX(DP) :: ham(nbnd,nbnd), eigvc(nbnd,nbnd)
  !
  ! The new eigenalues 
  REAL(DP) :: eigvl(nbnd)
  REAL(DP) :: eigvl_pert(nbnd)
  REAL(DP) :: eigvl_ks(nbnd)
  !
  INTEGER :: i, ibnd
  ! 
  REAL(DP) :: ehomo, elumo
  REAL(DP) :: ehomo_ks, elumo_ks
  REAL(DP) :: ehomo_pert, elumo_pert
  REAL(DP), EXTERNAL :: get_clock
  EXTERNAL :: ZGEMM, CDIAGH
  !
  COMPLEX(DP), ALLOCATABLE :: ham_aux(:,:)
  REAL(DP), ALLOCATABLE :: eigvl_ki(:)
  COMPLEX(DP), ALLOCATABLE :: eigvc_ki(:,:)
  INTEGER i_start, i_end, k
  !
  !
  ehomo=-1D+6
  elumo=+1D+6
  ehomo_ks=-1D+6
  elumo_ks=+1D+6
  ehomo_pert=-1D+6
  elumo_pert=+1D+6
  !
  WRITE( stdout, '(/,5X, "INFO: BUILD and DIAGONALIZE the KI HAMILTONIAN")')
  WRITE( stdout, '(  5X, "INFO: Unique Hamiltonian scheme using projectors:")')
  WRITE( stdout, '(  5X, "      build and diagonalize the KI Hamiltonian in")')
  WRITE( stdout, '(  5X, "      the basis of KS orbitals")')
  !
  ALLOCATE ( dH_wann_aux(num_wann, num_wann) )
  ALLOCATE ( evc_aux(npwx*npol, nbnd) )
  !
  DO ik = 1, nkstot_eff
    !
    dH_wann_aux(:,:) = dH_wann(ik, :,:)
    !
    ! Unique Hamiltonian diagonalized on the KS basis of the NSF calculation 
    ! 
    ik_pw = ik + (spin_component-1)*(nkstot_eff)
    WRITE( stdout, 9020 ) ( xk(i,ik_pw), i = 1, 3 )
    CALL get_buffer ( evc, nwordwfc, iuwfc, ik_pw )
    npw = ngk(ik_pw)
    !
    ! The KS Hamiltonian in the KS basis
    ham(:,:)=CMPLX(0.D0, 0.D0, kind=DP)
    DO i = 1, nbnd
      ham(i,i)    = et(i,ik_pw)
      eigvl_ks(i) = et(i,ik_pw)
    ENDDO
    !
    ehomo_ks = MAX ( ehomo_ks, et(num_wann_occ  , ik_pw) )
    IF (nbnd > num_wann_occ) elumo_ks = MIN ( elumo_ks, et(num_wann_occ+1, ik_pw) )
    !
    ! The Delta H_KI_ij = \sum_nm <phi_i|w_n> \Delta H_nm <w_m|phi_j>
    CALL dki_hamiltonian (evc, ik, nbnd, dH_wann_aux(:,:), deltah)
    !
    ! Add to the KS Hamiltonian
    ham(:,:) = ham(:,:) + deltah(:,:)
    !
#ifdef DEBUG
    WRITE(stdout, '(/, "dKI Hamiltonian at k = ", i4)') ik
    DO k = 1, 10
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(deltah(k,i)),AIMAG(deltah(k,i)), i=1,10)
    ENDDO
    !
    WRITE(stdout, '(/, "KI Hamiltonian at k = ", i4)') ik
    DO k = 1, 10
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(ham(k,i)),AIMAG(ham(k,i)), i=1,10)
    ENDDO
#endif
    !
    ! Because we have defined a uniq KI Hamiltonian, we can do a perturbative approach
    ! i.e. we keep only the diagonal part of the KI Hamiltoniana
    DO i = 1, nbnd
      eigvl_pert(i) = et(i,ik_pw) + DBLE(deltah(i,i))
    ENDDO
    ehomo_pert = MAX ( ehomo_pert, eigvl_pert(num_wann_occ ) )
    IF (nbnd > num_wann_occ) elumo_pert = MIN ( elumo_pert, eigvl_pert(num_wann_occ+1 ) )
    !
    IF (kcw_iverbosity .gt. 1 ) THEN
      WRITE(stdout,'(8x, "INFO: Empty states spectrum as a function of the # of orbitals")')
      !
      DO k = 1, nbnd-num_wann_occ
         !
         i_start = num_wann_occ+1; i_end = num_wann_occ+k
         !
         ALLOCATE (ham_aux(k,k), eigvl_ki(k), eigvc_ki(k,k))
         ham_aux(1:k,1:k) = ham(i_start:i_end,i_start:i_end)
         !
         CALL cdiagh( k, ham_aux, k, eigvl_ki, eigvc_ki )
         !
         IF (k.le.10) THEN
            WRITE(stdout,'(8x, I3, 10F10.4)') k, eigvl_ki(1:k)*rytoev  ! First 10 eigenvalues
         ELSE
            WRITE(stdout,'(8x, I3, 10F10.4)') k, eigvl_ki(1:10)*rytoev  ! First 10 eigenvalues
         ENDIF
         !
         DEALLOCATE (ham_aux)
         DEALLOCATE (eigvl_ki, eigvc_ki)
         !
      ENDDO
      WRITE(stdout,*)
    ENDIF
    !
    ! Diagonalize the KI Hamiltonian
    CALL CDIAGH( nbnd, ham, nbnd, eigvl, eigvc )
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
  DEALLOCATE (dH_wann_aux)
  DEALLOCATE (evc_aux)
  !
  9043 FORMAT(/,8x, 'KS  highest occupied level (ev): ',F10.4 )
  9042 FORMAT(/,8x, 'KS  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
  9045 FORMAT(  8x, 'KI  highest occupied level (ev): ',F10.4 )
  9044 FORMAT(  8x, 'KI  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
  9047 FORMAT(  8x, 'pKI highest occupied level (ev): ',F10.4 )
  9046 FORMAT(  8x, 'pKI highest occupied, lowest unoccupied level (ev): ',2F10.4 )
  9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
  901 FORMAT('          total cpu time spent up to now is ',F10.1,' secs' )
  !
  RETURN
  CONTAINS
  !
  ! !----------------------------------------------------------------
  SUBROUTINE dki_hamiltonian (evc, ik, h_dim, delta, deltah)
    !----------------------------------------------------------------
    !
    USE buffers,               ONLY : get_buffer
    USE control_kcw,           ONLY : num_wann
    USE wvfct,                 ONLY : npwx
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: h_dim
    COMPLEX(DP) :: delta (num_wann,num_wann)
    COMPLEX(DP), INTENT(IN) :: evc(npwx*npol,h_dim)
    INTEGER, INTENT(IN) :: ik
    COMPLEX(DP), INTENT(OUT) :: deltah(h_dim,h_dim)
    !
    INTEGER :: lrwannfc, ib, jb, nwann, mwann
    COMPLEX(DP) :: overlap
    !
    COMPLEX (DP) :: overlap_mat(nbnd,num_wann)
    !
    EXTERNAL :: ZGEMM
    !
    lrwannfc = num_wann*npwx*npol
    CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
    !
    deltah = CMPLX(0.D0, 0.D0, kind=DP)
    !
    CALL ZGEMM( 'C','N', nbnd, num_wann, npwx*npol, ONE, evc, npwx*npol, evc0, npwx*npol, &
    ZERO, overlap_mat, nbnd) 
    CALL mp_sum (overlap_mat, intra_bgrp_comm)
    !     
    DO ib = 1, nbnd
      DO jb = ib, nbnd
        !
        DO nwann = 1, num_wann
          DO mwann = 1, num_wann
            ! OLD 
            !overlap_in = SUM(CONJG(evc(1:npw,ib))*(evc0(1:npw,nwann)))
            !overlap_mj = SUM(CONJG(evc0(1:npw,mwann))*(evc(1:npw,jb)))
            !CALL mp_sum (overlap_in, intra_bgrp_comm)
            !CALL mp_sum (overlap_mj, intra_bgrp_comm)
            !overlap = overlap_in*overlap_mj
            ! 
            overlap = (overlap_mat(ib,nwann )) * CONJG(overlap_mat(jb,mwann))
            deltah(ib,jb) = deltah(ib,jb) + delta(nwann,mwann) * overlap
            !WRITE(*,'(5X, 2I5, 2F20.12, 2F20.12, 2F20.12)') nwann, mwann, delta(nwann,mwann), overlap, deltah(ib,jb)
          ENDDO
        ENDDO
        !WRITE(*,'(3X, 2I5, 2F20.12)') ib, jb, deltah(ib,jb)
        IF (ib /= jb) deltah(jb,ib) = CONJG(deltah(ib,jb))
        !
      ENDDO
    ENDDO
    !
    !
  END SUBROUTINE dki_hamiltonian
  !
  !
END SUBROUTINE koopmans_ham_uniq
