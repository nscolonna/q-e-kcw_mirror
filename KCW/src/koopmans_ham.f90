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
SUBROUTINE koopmans_ham (dH_wann)
  !---------------------------------------------------------------------
  !
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : nkstot, xk
  USE lsda_mod,              ONLY : nspin
  USE control_kcw,           ONLY : num_wann, Hamlt,evc0, spin_component, &
                                    num_wann_occ, iuwfc_wann
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx, npw, et, nbnd
  USE units_lr,              ONLY : lrwfc, iuwfc
  USE wavefunctions,         ONLY : evc
  USE buffers,               ONLY : get_buffer, save_buffer
  USE control_lr,            ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  ! The k point index 
  INTEGER :: ik, ik_pw
  !
  ! the KI hamiltonian, the KI contribution, and the new eigenvectors at a given k-point
  COMPLEX(DP) :: ham(num_wann,num_wann), eigvc(num_wann,num_wann)
  !
  COMPLEX(DP), INTENT (IN) :: dH_wann(nkstot/nspin,num_wann,num_wann)
  !
  ! The new eigenalues 
  REAL(DP) :: eigvl(num_wann)
  !
  INTEGER :: i, iwann
  ! 
  REAL(DP) :: ehomo, elumo
  REAL(DP) :: ehomo_ks, elumo_ks
  INTEGER  :: lrwannfc
  REAL(DP), EXTERNAL :: get_clock
  !
  !
  ehomo=-1D+6
  elumo=+1D+6
  ehomo_ks=-1D+6
  elumo_ks=+1D+6
  !
  WRITE( stdout, '(/,5X, "INFO: BUILD and DIAGONALIZE the KI HAMILTONIAN")')
  WRITE( stdout, '(  5X, "INFO: Standar scheme: diagonalize on the Wannier basis")')
  !
  DO ik = 1, nkstot/nspin
    !
    WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
    !
    ! Diagonalize the KS hamiltonian in the Wannier Gauge (just to check)
    ham(:,:) = Hamlt(ik,:,:) 
    CALL cdiagh( num_wann, ham, num_wann, eigvl, eigvc )
    !
    WRITE( stdout, '(10x, "KS  ",8F11.4)' ) (eigvl(iwann)*rytoev, iwann=1,num_wann)
    !
    ehomo_ks = MAX ( ehomo_ks, eigvl(num_wann_occ ) )
    IF (num_wann > num_wann_occ) elumo_ks = MIN ( elumo_ks, eigvl(num_wann_occ+1 ) )
    !
    ! Add the KI contribution to the KS Hamiltonian in the wannier basis
    Hamlt(ik,:,:) = Hamlt(ik,:,:) + dH_wann(ik,:,:) 
    !
#ifdef DEBUG
    WRITE(stdout, '(/, "KI Hamiltonian at k = ", i4)') ik
    DO iwann = 1, num_wann
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(Hamlt(ik, iwann,jwann)),AIMAG(Hamlt(ik,iwann,jwann)), jwann=1,num_wann)
    ENDDO
#endif
    !
    ! Diagonalize the KI hamitlonian in the Wannier Gauge
    ham(:,:) = Hamlt(ik,:,:) 
    CALL cdiagh( num_wann, ham, num_wann, eigvl, eigvc )
    WRITE( stdout, '(10x, "KI  ",8F11.4)' ) (eigvl(iwann)*rytoev, iwann=1,num_wann)
    WRITE(stdout, 901) get_clock('KCW')
    !
    ! Canonical wfc at each k point (overwrite the evc from DFT)
    lrwannfc = num_wann*npwx
    !write (*,'("NICOLA lrwannfc", i20)') lrwannfc, iuwfc_wann
    CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
    ! Retrive the ks function at k (in the Wannier Gauge)
    CALL ZGEMM( 'N','N', npw, num_wann, num_wann, ONE, evc0, npwx, eigvc, num_wann, &
                 ZERO, evc, npwx )
    lrwfc = nbnd * npwx
    ik_pw = ik + (spin_component-1)*(nkstot/nspin)
    !write (*,'("NICOLA lrwfc", i20)') lrwfc, iuwfc, nbnd, SIZE(evc)
    CALL save_buffer ( evc, lrwfc, iuwfc, ik_pw )
    !
    nbnd = num_wann
    DO iwann = 1, nbnd
      et(iwann,ik_pw) = eigvl(iwann)
    ENDDO
    !
    ehomo = MAX ( ehomo, eigvl(num_wann_occ ) )
    IF (num_wann > num_wann_occ) elumo = MIN ( elumo, eigvl(num_wann_occ+1 ) )
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
  ! For an isolated molecule the full Hamiltonian is also available
  ! OBSOLETE: moved to dH_KI_wann_supercell.f90 
  !ik=1
  !IF (nkstot/nspin == 1) CALL full_ham ( ik )
  !! ... formats
  !
9043 FORMAT(/,8x,'KS       highest occupied level (ev): ',F10.4 )
9042 FORMAT(/,8x, 'KS       highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9045 FORMAT(  8x, 'KI[2nd]  highest occupied level (ev): ',F10.4 )
9044 FORMAT(  8x, 'KI[2nd]  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
 901 FORMAT('          total cpu time spent up to now is ',F10.1,' secs' )
  !
  RETURN
  !
END SUBROUTINE koopmans_ham
