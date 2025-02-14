!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kcw_R_points
  !---------------------------------------------------------------------
  !
  !! This routine generates the R-points grid. Every R point
  !! corresponds to the position of primitive cell in a virtual
  !! supercell. R=0 is the origin and it corresponds to the 
  !! real primitive cell which all virtual cells are generated from.
  !
  USE cell_base,            ONLY : at
  USE control_kcw,          ONLY : Rvect, mp1, mp2 ,mp3, irvect, nkstot_eff
  USE control_kcw,          ONLY : Rvect_shifted, irvect_shifted, get_coulomb
  USE klist,                ONLY : nkstot
  USE lsda_mod,             ONLY : nspin
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, k, icell, num_R, mptot
  ! Number of unit cells ( = number of q points)
  !
  num_R = nkstot_eff
  mptot=mp1*mp2*mp3
  IF (num_R .ne. mptot) &
     CALL errore('kcw_R_points', ' Mismatch between num of kpoints and MP grid from input', num_R)
  !
  ALLOCATE (Rvect(3,num_R))
  ALLOCATE (iRvect(3,num_R))
  IF( get_coulomb ) THEN 
    ALLOCATE (Rvect_shifted(3,num_R))
    ALLOCATE (iRvect_shifted(3,num_R))
  END IF
  !
  WRITE(stdout,'(/,5X, "INFO: total number of primitive cell", i5)') num_R
  !
  IF ( nkstot == 1 ) THEN
    !
    Rvect(:,1) = 0.0d0
    irvect(:,1) = (/0,0,0/)
    !
    IF ( get_coulomb ) THEN
      Rvect_shifted(:,1) = (/0,0,0/)
    END IF
    !
  ELSE
    !
    ! "at" are in units of alat
    !
    icell = 0
    !
    DO i = 1, mp1
      DO j = 1, mp2
        DO k = 1, mp3
           !
           icell = icell + 1
           !
           Rvect(:,icell) = DBLE(i-1) * at(:,1) + &
                            DBLE(j-1) * at(:,2) + &
                            DBLE(k-1) * at(:,3)
           irvect(:,icell) = (/i-1,j-1,k-1/)
           IF( get_coulomb ) THEN
             iRvect_shifted(:, icell)  = (/i-1-mp1/2,j-1-mp2/2,k-1-mp3/2/)
             Rvect_shifted(:, icell) =  Rvect_shifted(1, icell) * at(:,1) + &
                                         Rvect_shifted(2, icell) * at(:,2) + &
                                         Rvect_shifted(3, icell) * at(:,3) 
           END IF
           !
        ENDDO
      ENDDO
    ENDDO
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE kcw_R_points

