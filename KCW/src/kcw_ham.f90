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
  !!
  !!  Non-collinear code written in 2022-24 by 
  !!  Antimo Marrazzo (SISSA, UniTS) and Nicola Colonna (PSI)
  !
  USE kinds,                 ONLY : DP
  USE control_kcw,           ONLY : do_bands, write_hr, h_uniq, num_wann, h_proj, corr_pc, corr_sc, nkstot_eff
  USE interpolation,         ONLY : interpolate_ham, dealloc_interpolation
  !
  USE io_rho_xml,            ONLY : write_scf
  USE io_files,              ONLY : prefix, iunwfc
  USE scf,                   ONLY : rho
  USE units_lr,              ONLY : iuwfc
  USE noncollin_module,      ONLY : noncolin, nspin_mag
  !
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: dH_wann(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dH_wann_proj(:)
  !
  ! 1) Set up for the KC calculation. 
  CALL kcw_setup_ham( )
  !
  ALLOCATE ( dH_wann(nkstot_eff,num_wann,num_wann) )
  ALLOCATE ( dH_wann_proj(num_wann) )
  !
  ! 2) compute the KI correction on the Wannier basis
  IF (corr_pc ) THEN 
    CALL dH_ki_quadratic ( dH_wann, dH_wann_proj )
  ELSEIF( corr_sc) THEN 
    ! This should never happen (see check in kcw_readin). Just for extra safety
    IF (noncolin) CALL errore( 'kcw_ham', 'which_odd="ki | pkipz" not implemented for noncollinear calculations',1 )
    CALL dH_ki_full ( 1, dH_wann)
    !!! Generalization for kpoint sampling. Eventually will replace dH_ki_full
    !CALL dH_ki_full_qpoints (1, dH_wann)
  ELSE
    ! This should never happen (see check in kcw_readin). Just for extra safety
    CALL errore( 'kcw_ham', 'Something wrong with which_odd',1 )
  ENDIF 
  !
  ! 3) Build UP the Koopmans Hamiltonian (depending on the scheme) 
  IF (h_proj) THEN 
    ! This is an alternative formulation based on a DFT+U like hamiltonian 
    ! see koopmans_ham_proj.f90 for details
    CALL koopmans_ham_proj ( dH_wann_proj )
    !
  ELSEIF ( h_uniq ) THEN 
    ! Use Projectors to build a unique Hamiltonian and diagonalize it 
    ! in the space of the KS orbital from the NSCF calculation
    CALL koopmans_ham_uniq ( dH_wann )
    !
  ELSE 
    ! Standard procedure: build and diagonalize the Hamiltonian in the 
    ! space spanned by the MLWFs
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
  DEALLOCATE (dH_wann, dH_wann_proj) 
  !
  ! WRITE data file
  iunwfc = iuwfc
  prefix = TRIM(prefix)//"_kcw"
  ! Append an extra postfix to save the temporary files of different Hamiltonian schemes on an different dirs
  IF ( h_uniq ) prefix = TRIM(prefix)//"_uniq-H"
  IF ( h_proj ) prefix = TRIM(prefix)//"_proj-H"
  CALL write_scf(rho, nspin_mag)
  !CALL punch('config-only')
  CALL punch ('all')
  !
  CALL clean_pw( .TRUE. )
  CALL close_kcw ( ) 
  !
END SUBROUTINE kcw_ham
