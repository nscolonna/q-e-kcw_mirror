!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE kcw_kpoint_grid( nrot, time_reversal, skip_equivalence, s, t_rev, &
                        bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk_, wk )
  !-----------------------------------------------------------------------
  !!  Automatic generation of a uniform grid of k-points.
  !
  USE kinds, ONLY: DP
  USE klist,                 ONLY : xk
  USE cell_base,             ONLY : at
  USE io_global,             ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nrot
  !! number of bravais lattice symmetries
  INTEGER, INTENT(IN) :: npk
  !! max number of k-points
  INTEGER, INTENT(IN) :: k1
  !! the offset from the origin, direction 1
  INTEGER, INTENT(IN) :: k2
  !! the offset from the origin, direction 2
  INTEGER, INTENT(IN) :: k3
  !! the offset from the origin, direction 3
  INTEGER, INTENT(IN) :: nk1
  !! the special-point grid, direction 1
  INTEGER, INTENT(IN) :: nk2
  !! the special-point grid, direction 2
  INTEGER, INTENT(IN) :: nk3
  !! the special-point grid, direction 3
  INTEGER, INTENT(IN) :: t_rev(48)
  !! time reversal flag, for noncolinear magnetism
  INTEGER, INTENT(IN) :: s(3,3,48)
  !! symmetry matrices, in crystal axis
  LOGICAL, INTENT(IN) :: time_reversal
  !! if .TRUE. the system has time reversal symmetry
  LOGICAL, INTENT(IN) :: skip_equivalence
  !! if .TRUE. skip check of k-points equivalence
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! bg(:,i) are the reciprocal lattice vectors, b_i,
  !! in tpiba=2pi/alat units: b_i(:) = bg(:,i)/tpiba
  INTEGER,  INTENT(out) :: nks
  !! number of k points
  REAL(DP), INTENT(out) :: xk_(3,npk)
  !! coordinates of k points
  REAL(DP), INTENT(out) :: wk(npk)
  !! weight of k points
  !
  ! ... local variables
  !
  REAL(DP) :: xkr(3), fact, xx, yy, zz
  REAL(DP), ALLOCATABLE :: xkg(:,:), wkk(:)
  INTEGER :: nkr, i, j, k, ns, n, nk, nkk
  INTEGER, ALLOCATABLE :: equiv(:)
  LOGICAL :: in_the_list
  REAL(DP), PARAMETER :: eps=1.0d-5
  !
  nkr=nk1*nk2*nk3
  ALLOCATE (xkg( 3,nkr),wkk(nkr))
  xkg=xk
  CALL cryst_to_cart(nkr, xkg, at, -1)
  ALLOCATE (equiv( nkr))
  !
  !DO i=1,nk1
  !   DO j=1,nk2
  !      DO k=1,nk3
  !         !  this is nothing but consecutive ordering
  !         n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
  !         !  xkg are the components of the complete grid in crystal axis
  !         xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
  !         xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
  !         xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
  !      ENDDO
  !   ENDDO
  !ENDDO

  !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
  !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)

  DO nk=1,nkr
     equiv(nk)=nk
  ENDDO

  IF ( skip_equivalence ) THEN
    CALL infomsg('kcw_kpoint_grid', 'ATTENTION: skip check of k-points equivalence')
    wkk = 1.d0
  ELSE
    DO nk=1,nkr
    !WRITE(*,*) "NICOLA ik=", nk
    !WRITE(*,'(a,3F12.4)') "NICOLA xk(ik)=", xkg(:,nk)
    !  check if this k-point has already been found equivalent to another
      IF (equiv(nk) == nk) THEN
        wkk(nk)   = 1.0d0
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
        DO ns=1,nrot
           !WRITE(*,*) "NICOLA is=", ns
           DO i=1,3
              xkr(i) = s(i,1,ns) * xkg(1,nk) &
                     + s(i,2,ns) * xkg(2,nk) &
                     + s(i,3,ns) * xkg(3,nk)
              xkr(i) = xkr(i) - nint( xkr(i) )
           ENDDO
           IF(t_rev(ns)==1) xkr = -xkr
           xx = xkr(1)*nk1 - 0.5d0*k1
           yy = xkr(2)*nk2 - 0.5d0*k2
           zz = xkr(3)*nk3 - 0.5d0*k3
           in_the_list = abs(xx-nint(xx))<=eps .and. &
                         abs(yy-nint(yy))<=eps .and. &
                         abs(zz-nint(zz))<=eps
           IF (in_the_list) THEN
              i = mod ( nint ( xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
              j = mod ( nint ( xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
              k = mod ( nint ( xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
              IF (n>nk .and. equiv(n)==n) THEN
                 equiv(n) = nk
                 wkk(nk)=wkk(nk)+1.0d0
                 !WRITE(*,'(a,3i5, 1F8.4)') "NICOLA", n, equiv(n), nk, wkk(nk)
              ELSE
                 IF (equiv(n)/=nk .or. n<nk ) CALL errore('kcw_kpoint_grid', &
                    'something wrong in the checking algorithm',1)
              ENDIF
           ENDIF
           IF ( time_reversal ) THEN
              xx =-xkr(1)*nk1 - 0.5d0*k1
              yy =-xkr(2)*nk2 - 0.5d0*k2
              zz =-xkr(3)*nk3 - 0.5d0*k3
              in_the_list=abs(xx-nint(xx))<=eps.and.abs(yy-nint(yy))<=eps &
                                                 .and. abs(zz-nint(zz))<=eps
              IF (in_the_list) THEN
                 i = mod ( nint (-xkr(1)*nk1 - 0.5d0 * k1 + 2*nk1), nk1 ) + 1
                 j = mod ( nint (-xkr(2)*nk2 - 0.5d0 * k2 + 2*nk2), nk2 ) + 1
                 k = mod ( nint (-xkr(3)*nk3 - 0.5d0 * k3 + 2*nk3), nk3 ) + 1
                 n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                 IF (n>nk .and. equiv(n)==n) THEN
                    equiv(n) = nk
                    wkk(nk)=wkk(nk)+1.0d0
                 ELSE
                    IF (equiv(n)/=nk.or.n<nk) CALL errore('kpoint_grid', &
                    'something wrong in the checking algorithm',2)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  !  count irreducible points and order them
  nks=0
  fact=0.0d0
  DO nk=1,nkr
     IF (equiv(nk)==nk) THEN
        nks=nks+1
        IF (nks>npk) CALL errore('kpoint_grid','too many k-points',1)
        wk(nks) = wkk(nk)
        fact    = fact+wk(nks)
        !  bring back into to the first BZ
        DO i=1,3
           xk_(i,nks) = xkg(i,nk)!-nint(xkg(i,nk))
        ENDDO
     ENDIF
  ENDDO
  !  go to cartesian axis (in units 2pi/a0)
  CALL cryst_to_cart(nks,xk_,bg,1)
  !  normalize weights to one
  DO nk=1,nks
     wk(nk) = wk(nk)/fact
  ENDDO
  !
  ! Summary of the search 
  WRITE(stdout, '(3x,"ik_IBZ   ik        equivalent ik in the full BZ")') 
  DO nk=1,nkr
     IF (equiv(nk)==nk) THEN
        WRITE (stdout,'(i5, 4x, i5, a)', advance='no') nks, nk, " <-- "
        DO nkk = 1, nkr
          IF (equiv(nkk)==nk) THEN
             WRITE (stdout,'(i5)', advance='no') nkk
          ENDIF
        ENDDO
     write(stdout,*)
     ENDIF
     FLUSH(stdout)
  ENDDO
  !
  WRITE(stdout,'(/, 3x, "nkstot_IBZ=", i5, " nsym", I5)') nks, nrot
  WRITE(stdout,'(3x, "xk (",i5, ")", 3F10.4,3x, F10.6)') (nk, xk_(1:3,nk), wk(nk), nk=1,nks) 
  !
  DEALLOCATE(equiv)
  DEALLOCATE(xkg,wkk)

  RETURN
END SUBROUTINE kcw_kpoint_grid
