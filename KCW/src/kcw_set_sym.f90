!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
#define ONE (0.D0,1.D0)
!#define DEBUG
!-----------------------------------------------------------------------
SUBROUTINE kcw_set_symm(nr1, nr2, nr3, nr1x, nr2x, nr3x)
  !-----------------------------------------------------------------------
  ! based on exx_set_symm in exx_base.f90
  !
  USE symm_base,        ONLY : s, ft, nsym
  USE control_kcw,      ONLY : rir
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x 
  !
  ! ... local variables
  !
  INTEGER :: ikq, isym, i,j,k, ri,rj,rk, ir, nxxs
  INTEGER, allocatable :: ftau(:,:), s_scaled(:,:,:)
  !
  nxxs = nr1x*nr2x*nr3x
  !
  ALLOCATE( rir(nxxs,nsym) )
  !
  rir = 0
  ALLOCATE ( ftau(3,nsym), s_scaled(3,3,nsym) )
  CALL scale_sym_ops (nsym, s, ft, nr1, nr2, nr3, s_scaled, ftau)

  DO isym = 1, nsym
     DO k = 1, nr3
        DO j = 1, nr2
           DO i = 1, nr1
              CALL rotate_grid_point( s_scaled(1,1,isym), ftau(1,isym), &
                   i, j, k, nr1, nr2, nr3, ri, rj, rk )
              ir = i + (j-1)*nr1x + (k-1)*nr1x*nr2x
              rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE ( s_scaled, ftau )
  !
END SUBROUTINE kcw_set_symm

