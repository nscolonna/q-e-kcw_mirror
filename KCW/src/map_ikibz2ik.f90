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
SUBROUTINE  map_ikibz2ik (xk_ibz, nkstot_ibz, ik_ibz2ik, isk_ibz)
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
  USE lsda_mod,              ONLY : lsda, isk
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (OUT) :: ik_ibz2ik(nkstot_ibz)
  INTEGER, INTENT (IN)  :: nkstot_ibz
  INTEGER, INTENT (IN)  :: isk_ibz(nkstot_ibz)
  REAL(DP), INTENT (IN) :: xk_ibz(3,nkstot_ibz)
  INTEGER  :: ik_ibz, ik
  REAL(DP) :: dk(3) 
  LOGICAL :: found(nkstot_ibz)
  !
  ! 
  found = .FALSE.
  !
  DO ik_ibz = 1, nkstot_ibz
  !WRITE(*,'("NICOLA", i5, 3F12.4)') ik_ibz, xk_ibz(:,ik_ibz)
    DO ik = 1, nkstot
      IF (lsda .AND. isk(ik) /= isk_ibz(ik_ibz)) CYCLE 
      !WRITE(*,'(3x, "NICOLA FBZ", i5, 3F12.4)') ik, xk(:,ik)
      dk(:) = xk_ibz(:,ik_ibz)-xk(:,ik)
      IF ( all( abs( dk ) < 1d-5 ) .AND. .NOT. found(ik_ibz) ) THEN 
         ik_ibz2ik(ik_ibz) = ik
         found(ik_ibz) = .TRUE. 
         CYCLE
      ENDIF
    ENDDO
    write(*,*) ik_ibz, ik_ibz2ik(ik_ibz)
  ENDDO 
  IF ( ANY (.NOT. found)) CALL errore('map ikibz2ik', 'No match found',1)
  !
  RETURN
  !
END subroutine map_ikibz2ik
