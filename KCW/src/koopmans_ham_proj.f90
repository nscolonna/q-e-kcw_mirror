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
SUBROUTINE koopmans_ham_proj ( dH_wann )
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
  USE klist,                 ONLY : nkstot, xk, ngk
  USE lsda_mod,              ONLY : nspin
  USE control_kcw,           ONLY : num_wann, evc0, spin_component, &
                                    num_wann_occ, iuwfc_wann
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx, npw, et, nbnd
  USE units_lr,              ONLY : iuwfc
  USE wavefunctions,         ONLY : evc
  USE buffers,               ONLY : get_buffer, save_buffer
  !
  USE io_files,              ONLY : nwordwfc
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! The k point index 
  INTEGER :: ik, ik_pw
  !
  ! the KI hamiltonian on the Wannier basis <w_i|dh_j|w_j> 
  COMPLEX(DP), INTENT (IN) :: dH_wann(nkstot/nspin,num_wann,num_wann)
  COMPLEX(DP), ALLOCATABLE :: dH_wann_aux(:,:)
  ! 
  ! the KI operator on the KS basis of the NSCF calculation
  COMPLEX(DP) :: deltah(nbnd,nbnd)
  ! the new hamitonain, and the new eigenvalues and eigenvectors at a given k-point
  COMPLEX(DP) :: ham(nbnd,nbnd), eigvc(nbnd,nbnd), eigvc_all(nbnd,nbnd,nkstot/nspin)
  !
  REAL(DP) :: et_ki(nbnd,nkstot)
  !
  ! The new eigenalues 
  REAL(DP) :: eigvl(nbnd)
  !
  INTEGER :: i, ibnd
  ! 
  REAL(DP) :: ehomo, elumo
  REAL(DP) :: ehomo_ks, elumo_ks
  REAL(DP), EXTERNAL :: get_clock
  EXTERNAL :: ZGEMM, CDIAGH
  !
  !
  ehomo=-1D+6
  elumo=+1D+6
  ehomo_ks=-1D+6
  elumo_ks=+1D+6
  !
  WRITE( stdout, '(/,5X, "INFO: BUILD and DIAGONALIZE the KI HAMILTONIAN")')
  WRITE( stdout, '(  5X, "INFO: Projection scheme")')
  !
  ALLOCATE ( dH_wann_aux(num_wann, num_wann) )
  !
  DO ik = 1, nkstot/nspin
    !
    dH_wann_aux(:,:) = dH_wann(ik, :,:)
    !
    ! Unique Hamiltonian diagonalized on the KS basis of the NSF calculation 
    ! 
    ik_pw = ik + (spin_component-1)*(nkstot/nspin)
    WRITE( stdout, 9020 ) ( xk(i,ik_pw), i = 1, 3 )
    CALL get_buffer ( evc, nwordwfc, iuwfc, ik_pw )
    npw = ngk(ik_pw)
    !
    ! The KS Hamiltonian in the KS basis
    ham(:,:)=CMPLX(0.D0, 0.D0, kind=DP)
    DO i = 1, nbnd
      ham(i,i)=et(i,ik_pw)
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
    ! Diagonalize it 
    CALL CDIAGH( nbnd, ham, nbnd, eigvl, eigvc )
    !
    ! Store the eigenvalues and eigenvector at this k point 
    et_ki(1:nbnd,ik)=eigvl(1:nbnd)
    eigvc_all(:,:,ik) = eigvc(:,:)
    !
    ehomo = MAX ( ehomo, et_ki(num_wann_occ, ik ) )
    IF (nbnd > num_wann_occ) elumo = MIN ( elumo, et_ki(num_wann_occ+1, ik ) )
    !
    WRITE( stdout, '(10x, "KS  ",8F11.4)' ) (et(ibnd,ik_pw)*rytoev, ibnd=1,nbnd)
    WRITE( stdout, '(10x, "KI  ",8F11.4)' ) (et_ki(ibnd,ik)*rytoev, ibnd=1,nbnd)
    WRITE(stdout, 901) get_clock('KCW')
    !
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
  !Overwrite et and evc
  DO ik = 1, nkstot/nspin
    ik_pw = ik + (spin_component-1)*(nkstot/nspin)
    et(1:nbnd, ik_pw) = et_ki(1:nbnd,ik)
    !
    ! Canonical wfc at each k point (overwrite the evc from DFT)
    CALL get_buffer ( evc, nwordwfc, iuwfc, ik_pw )
    ! Retrive the ks function at k (in the Wannier Gauge)
    eigvc(:,:) = eigvc_all(:,:,ik)
    CALL ZGEMM( 'N','N', npw, nbnd, nbnd, ONE, evc, npwx, eigvc, nbnd, &
                 ZERO, evc, npwx )
    !write (*,'("NICOLA lrwfc", i20)') lrwfc, iuwfc, nbnd, SIZE(evc)
    CALL save_buffer ( evc, nwordwfc, iuwfc, ik_pw )
    !
  ENDDO
  !
  DEALLOCATE (dH_wann_aux)
  !
9043 FORMAT(/,8x, 'KS       highest occupied level (ev): ',F10.4 )
9042 FORMAT(/,8x, 'KS       highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9045 FORMAT(  8x, 'KI[2nd]  highest occupied level (ev): ',F10.4 )
9044 FORMAT(  8x, 'KI[2nd]  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
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
    COMPLEX(DP), INTENT(IN) :: evc(npwx,h_dim)
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
    lrwannfc = num_wann*npwx
    CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
    !
    deltah = CMPLX(0.D0, 0.D0, kind=DP)
    !
    CALL ZGEMM( 'C','N', nbnd, num_wann, npwx, ONE, evc, npwx, evc0, npwx, &
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
          !WRITE(*,'(3X, 2I5, 2F20.12, F20.12, 2F20.12)') ibnd, iwann, delta(iwann), occ_mat(iwann), overlap
         ENDDO
        ENDDO
        IF (ib /= jb) deltah(jb,ib) = CONJG(deltah(ib,jb))
        !
      ENDDO
    ENDDO
    !
    !
  END SUBROUTINE dki_hamiltonian
  !
  !
END SUBROUTINE koopmans_ham_proj
