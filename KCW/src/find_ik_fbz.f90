!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!-----------------------------------------------------------------------
SUBROUTINE find_ik_fbz (ik, ik_fbz)
  !-----------------------------------------------------------------------
  !
  !! Given the index of the kpoint in the IBZ (ik), this routine 
  !! finds the corresponding index of the kpoint in the full BZ
  !! Needed to "select" the correct U matrix read from W90
  !! 
  !
  USE kinds,                 ONLY : DP
  USE control_kcw,           ONLY: xk_fbz, mp1, mp2, mp3
  USE klist,                 ONLY :xk, nkstot, nks
  USE io_global,             ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: ik
  INTEGER, INTENT(OUT) :: ik_fbz
  INTEGER :: nkstot_fbz
  REAL(DP) :: dk(3) 
  LOGICAL :: found 
  !
  ! 
  nkstot_fbz = mp1*mp2*mp3 
  found = .FALSE.
  !
  !WRITE(*,'("NICOLA", i5, 3F12.4)') ik, xk(:,ik)
  DO ik_fbz = 1, nkstot_fbz 
     !WRITE(*,'("NICOLA FBZ", i5, 3F12.4)') ik_fbz, xk_fbz(:,ik_fbz)
     dk(:)=xk(:,ik)-xk_fbz(:,ik_fbz)
     IF ( all( abs( dk ) < 1d-5 ) ) THEN 
        !WRITE(stdout,'("Match found", 2I5, 3x, 6F7.4 )') ik, ik_fbz, xk(:,ik),xk_fbz(:,ik_fbz)
        found = .TRUE. 
        RETURN 
     ENDIF
  ENDDO 
  IF ( .NOT. found) CALL errore('find_ik_fbz', 'No match found',1)
  !
  RETURN
  !
END subroutine find_ik_fbz
