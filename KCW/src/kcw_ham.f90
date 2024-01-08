!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
SUBROUTINE kcw_ham
  !-----------------------------------------------------------------
  !
  !!  This is one the main subroutines of the KCW code to build up the 
  !!  KC hamiltonian in Reciprocal Space. It reads the output of 
  !!  a PWscf calculation and the U matrices from W90
  !!   
  !!  Code written by Nicola Colonna and Riccardo de Gennaro (EPFL April 2019) 
  !
  USE kinds,                 ONLY : DP
  USE control_kcw,           ONLY : do_bands, write_hr, h_proj, num_wann
  USE interpolation,         ONLY : interpolate_ham, dealloc_interpolation
  !
  USE io_rho_xml,            ONLY : write_scf
  USE io_files,              ONLY : prefix, iunwfc
  USE scf,                   ONLY : rho
  USE lsda_mod,              ONLY : nspin
  USE units_lr,              ONLY : iuwfc
  USE klist,                 ONLY : nkstot
  !
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: dH_wann(:,:,:)
  !
  ! 1) Set up for the KC calculation. 
  CALL kcw_setup_ham ( )
  !
  ALLOCATE ( dH_wann(nkstot/nspin,num_wann,num_wann) )
  ! 2) compute the KI correction on the Wannier basis
  CALL dH_ki_wann ( dH_wann )
  !
  IF ( h_proj ) THEN 
     ! Use Projectors to build a unique Hamiltonian and diagonalize it 
     ! in the space of the KS orbital from the NSCF calculation
     CALL koopmans_ham_proj ( dH_wann )
     ! This is an alternative formulation based on a DFT+U like hamiltonian 
     ! see koopmans_ham_proj_alternative.f90 for details
     !CALL koopmans_ham_proj_alternative ( )
  ELSE 
    ! Standard procedure using MLWFs
    CALL koopmans_ham ( dH_wann )
    !
    ! 3) If do_bands=TRUE interpolate H(k) and prints bands
    IF ( do_bands ) CALL interpolate_ham( )
    !
    ! 4) If write_hr=TRUE write H(R) to file
    IF ( write_hr ) CALL write_hr_to_file( )
    !
    IF (do_bands) CALL dealloc_interpolation( )
    !
  ENDIF
  !
  DEALLOCATE (dH_wann) 
  !
  ! WRITE data file
  iunwfc = iuwfc
  prefix = TRIM(prefix)//"_kcw"
  CALL write_scf(rho, nspin)
  !CALL punch('config-only')
  CALL punch ('all')
  !
  CALL clean_pw( .TRUE. )
  CALL close_kcw ( ) 
  !
END SUBROUTINE kcw_ham
