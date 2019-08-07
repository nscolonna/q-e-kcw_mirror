!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------=!
   MODULE fft_smallbox
!=----------------------------------------------------------------------=!

!! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
       USE iso_c_binding
       USE fft_param
       IMPLICIT NONE
       SAVE

       PRIVATE
       PUBLIC :: cft_b, cft_b_omp_init, cft_b_omp

       INTERFACE cft_b
          MODULE PROCEDURE cft_b_dp, cft_b_sp
       END INTERFACE
       INTERFACE cft_b_omp
          MODULE PROCEDURE cft_b_omp_dp, cft_b_omp_sp
       END INTERFACE


! ...   Local Parameter

        !   Workspace that is statically allocated is defined here
        !   in order to avoid multiple copies of the same workspace
        !   lwork:   Dimension of the work space array (if any)

        INTEGER   :: cft_b_dims( 3 )
!$omp threadprivate (cft_b_dims)
        TYPE(C_PTR) :: cft_b_bw_planz = C_NULL_PTR
!$omp threadprivate (cft_b_bw_planz)
        TYPE(C_PTR) :: cft_b_bw_planx = C_NULL_PTR
!$omp threadprivate (cft_b_bw_planx)
        TYPE(C_PTR) :: cft_b_bw_plany = C_NULL_PTR
!$omp threadprivate (cft_b_bw_plany)
        LOGICAL   :: use_single_precision
!$omp threadprivate (use_single_precision)

!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!

!
!=----------------------------------------------------------------------=!
!
!
!
!         3D parallel FFT on sub-grids
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_b_dp ( f, nx, ny, nz, ldx, ldy, ldz, imin2, imax2, imin3, imax3, sgn )

!     driver routine for 3d complex fft's on box grid, parallel case
!     fft along z for all xy values 
!     fft along y is done only for the local  z values: i.e. z-planes with imin3 <= nz <= imax3.
!     fft along x is done only for the local yz values: i.e. z-planes with imin3 <= nz <= imax3 
!                                                        and y-planes with imin2 <= ny <= imax2.
!     implemented for FFTW, only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
      USE fftw_interfaces
      implicit none
      integer nx,ny,nz,ldx,ldy,ldz,imin2,imax2,imin3,imax3,sgn
      complex(dp) :: f(:)

      integer isign, naux, ibid, k
      integer nplanes
      real(DP) :: tscale

      integer :: ip, i, first_index, how_many_y
      integer, save :: icurrent = 1
      integer, save :: dims( 3, ndims ) = -1

      TYPE(C_PTR), save :: bw_planz(  ndims ) = C_NULL_PTR
      TYPE(C_PTR), save :: bw_planx(  ndims ) = C_NULL_PTR
      TYPE(C_PTR), save :: bw_plany(  ndims ) = C_NULL_PTR

      isign = -sgn
      tscale = 1.0_DP

      if ( isign > 0 ) then
         call fftx_error__('cft_b','not implemented',isign)
      end if
!
! 2d fft on xy planes - only needed planes are transformed
! note that all others are left in an unusable state
!
      nplanes = imax3 - imin3 + 1

      !
      !   Here initialize table only if necessary
      !

      ip = -1
      DO i = 1, ndims

        !   first check if there is already a table initialized
        !   for this combination of parameters

        IF ( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. ( nz == dims(3,i) ) ) THEN
           ip = i
           EXIT
        END IF

      END DO

      IF( ip == -1 ) THEN

        !   no table exist for these parameters
        !   initialize a new one

        if ( C_ASSOCIATED(bw_planz(icurrent)) ) &
             call DESTROY_PLAN_1D( bw_planz(icurrent) )
        call CREATE_PLAN_1D( bw_planz(icurrent), nz, 1 )

        if ( C_ASSOCIATED(bw_planx(icurrent)) ) &
             call DESTROY_PLAN_1D( bw_planx(icurrent) )
        call CREATE_PLAN_1D( bw_planx(icurrent), nx, 1 )

        if ( C_ASSOCIATED(bw_plany(icurrent)) ) &
             call DESTROY_PLAN_1D( bw_plany(icurrent) )
        call CREATE_PLAN_1D( bw_plany(icurrent), ny, 1 )

        dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
        ip = icurrent
        icurrent = MOD( icurrent, ndims ) + 1

      END IF

      !
      !  fft along Z
      !
      call FFTW_INPLACE_DRV_1D( bw_planz(ip), ldx*ldy, f(1), ldx*ldy, 1 )
     
      do k = imin3, imax3
      !
      !  fft along Y
      !
        first_index = (k-1)*ldx*ldy + 1
        call FFTW_INPLACE_DRV_1D( bw_plany(ip), nx, f(first_index), ldx, 1 )
      !
      !  fft along X
      !
        first_index = first_index + (imin2-1)*ldx ; how_many_y = imax2 + 1 - imin2
        call FFTW_INPLACE_DRV_1D( bw_planx(ip), how_many_y, f(first_index), 1, ldx )

      end do   

      RETURN
   END SUBROUTINE cft_b_dp

!
!=----------------------------------------------------------------------=!
!
!
!
!   3D parallel FFT on sub-grids, to be called inside OpenMP region
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_b_omp_init ( nx, ny, nz, single_precision )

!     driver routine for 3d complex fft's on box grid, init subroutine
!
      USE fftw_interfaces
      implicit none
      integer, INTENT(IN) :: nx,ny,nz
      LOGICAL, OPTIONAL, INTENT(IN) :: single_precision
      LOGICAL :: sp
      !
      !   Here initialize table 
      !

!$omp parallel default(none) private(sp) shared(nx, ny, nz, single_precision)

      sp = .FALSE.
      IF( PRESENT( single_precision ) ) THEN
         sp = single_precision
      END IF
      use_single_precision = sp

      IF( .NOT. C_ASSOCIATED(cft_b_bw_planz) ) THEN
         IF(sp) THEN
            CALL FLOAT_CREATE_PLAN_1D( cft_b_bw_planz, nz, 1 )
         ELSE
            CALL CREATE_PLAN_1D( cft_b_bw_planz, nz, 1 )
         ENDIF
         cft_b_dims(3) = nz
      END IF
      IF( .NOT. C_ASSOCIATED(cft_b_bw_planx) ) THEN
         IF(sp) THEN
            CALL FLOAT_CREATE_PLAN_1D( cft_b_bw_planx, nx, 1 )
         ELSE
            CALL CREATE_PLAN_1D( cft_b_bw_planx, nx, 1 )
         ENDIF
         cft_b_dims(1) = nx
      END IF
      IF( .NOT. C_ASSOCIATED(cft_b_bw_plany) ) THEN
         IF(sp) THEN
            CALL FLOAT_CREATE_PLAN_1D( cft_b_bw_plany, ny, 1 )
         ELSE
            CALL CREATE_PLAN_1D( cft_b_bw_plany, ny, 1 )
         ENDIF
         cft_b_dims(2) = ny
      END IF

!$omp end parallel

     RETURN
   END SUBROUTINE cft_b_omp_init


   SUBROUTINE cft_b_omp_dp ( f, nx, ny, nz, ldx, ldy, ldz, imin2, imax2, imin3, imax3, sgn )

!     driver routine for 3d complex fft's on box grid, parallel (MPI+OpenMP) case
!     fft along z for all xy values 
!     fft along y is done only for the local  z values: i.e. z-planes with imin3 <= nz <= imax3.
!     fft along x is done only for the local yz values: i.e. z-planes with imin3 <= nz <= imax3 
!                                                        and y-planes with imin2 <= ny <= imax2.
!     implemented ONLY for internal fftw, and only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
!     This driver is meant for calls inside parallel OpenMP sections
!
      USE fftw_interfaces
      implicit none
      integer, INTENT(IN) :: nx,ny,nz,ldx,ldy,ldz,imin2,imax2,imin3,imax3,sgn
      complex(dp) :: f(:)

      INTEGER, SAVE :: k, first_index, how_many_y
!$omp threadprivate (k,first_index,how_many_y)

      if ( use_single_precision ) then
         CALL fftx_error__('cft_b_omp','double precision driver with single precision initialization',1)
      end if
      if ( -sgn > 0 ) then
         CALL fftx_error__('cft_b_omp','forward transform not implemented',1)
      end if

      IF ( .NOT. C_ASSOCIATED(cft_b_bw_planz) .or. &
           .NOT. C_ASSOCIATED(cft_b_bw_planx) .or. &
           .NOT. C_ASSOCIATED(cft_b_bw_plany) ) THEN
         CALL fftx_error__('cft_b_omp','plan not initialized',1)
      END IF

      !  consistency check

      IF ( ( nx /= cft_b_dims(1) ) .or. ( ny /= cft_b_dims(2) ) .or. ( nz /= cft_b_dims(3) ) ) THEN
         CALL fftx_error__('cft_b_omp', 'dimensions are inconsistent with the existing plan',1) 
      END IF

      !
      !  fft along Z
      !
      call FFTW_INPLACE_DRV_1D( cft_b_bw_planz, ldx*ldy, f(1), ldx*ldy, 1 )

      do k = imin3, imax3
      !
      !  fft along Y
      !
        first_index = (k-1)*ldx*ldy + 1
        call FFTW_INPLACE_DRV_1D( cft_b_bw_plany, nx, f(first_index), ldx, 1 )
      !
      !  fft along X
      !
        first_index = first_index + (imin2-1)*ldx ; how_many_y = imax2 + 1 - imin2
        call FFTW_INPLACE_DRV_1D( cft_b_bw_planx, how_many_y, f(first_index), 1, ldx )
      end do   

     RETURN
   END SUBROUTINE cft_b_omp_dp
!
!=----------------------------------------------------------------------=!
!
!
!
!         NOW Single precision drivers
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_b_sp ( f, nx, ny, nz, ldx, ldy, ldz, imin2, imax2, imin3, imax3, sgn )

!     driver routine for 3d complex fft's on box grid, parallel case
!     fft along z for all xy values 
!     fft along y is done only for the local  z values: i.e. z-planes with imin3 <= nz <= imax3.
!     fft along x is done only for the local yz values: i.e. z-planes with imin3 <= nz <= imax3 
!                                                        and y-planes with imin2 <= ny <= imax2.
!     implemented for FFTW, only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
      USE fftw_interfaces
      implicit none
      integer nx,ny,nz,ldx,ldy,ldz,imin2,imax2,imin3,imax3,sgn
      complex(sp) :: f(:)

      integer isign, naux, ibid, k
      integer nplanes
      real(sp) :: tscale

      integer :: ip, i, first_index, how_many_y
      integer, save :: icurrent = 1
      integer, save :: dims( 3, ndims ) = -1

      TYPE(C_PTR), save :: bw_planz(  ndims ) = C_NULL_PTR
      TYPE(C_PTR), save :: bw_planx(  ndims ) = C_NULL_PTR
      TYPE(C_PTR), save :: bw_plany(  ndims ) = C_NULL_PTR

      isign = -sgn
      tscale = 1.0_DP

      if ( isign > 0 ) then
         call fftx_error__('cft_b','not implemented',isign)
      end if
!
! 2d fft on xy planes - only needed planes are transformed
! note that all others are left in an unusable state
!
      nplanes = imax3 - imin3 + 1

      !
      !   Here initialize table only if necessary
      !

      ip = -1
      DO i = 1, ndims

        !   first check if there is already a table initialized
        !   for this combination of parameters

        IF ( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. ( nz == dims(3,i) ) ) THEN
           ip = i
           EXIT
        END IF

      END DO

      IF( ip == -1 ) THEN

        !   no table exist for these parameters
        !   initialize a new one

        if ( C_ASSOCIATED(bw_planz(icurrent)) ) &
             call FLOAT_DESTROY_PLAN_1D( bw_planz(icurrent) )
        call FLOAT_CREATE_PLAN_1D( bw_planz(icurrent), nz, 1 )

        if ( C_ASSOCIATED(bw_planx(icurrent)) ) &
             call FLOAT_DESTROY_PLAN_1D( bw_planx(icurrent) )
        call FLOAT_CREATE_PLAN_1D( bw_planx(icurrent), nx, 1 )

        if ( C_ASSOCIATED(bw_plany(icurrent)) ) &
             call FLOAT_DESTROY_PLAN_1D( bw_plany(icurrent) )
        call FLOAT_CREATE_PLAN_1D( bw_plany(icurrent), ny, 1 )

        dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
        ip = icurrent
        icurrent = MOD( icurrent, ndims ) + 1

      END IF

      !
      !  fft along Z
      !
      call FLOAT_FFTW_INPLACE_DRV_1D( bw_planz(ip), ldx*ldy, f(1), ldx*ldy, 1 )
     
      do k = imin3, imax3
      !
      !  fft along Y
      !
        first_index = (k-1)*ldx*ldy + 1
        call FLOAT_FFTW_INPLACE_DRV_1D( bw_plany(ip), nx, f(first_index), ldx, 1 )
      !
      !  fft along X
      !
        first_index = first_index + (imin2-1)*ldx ; how_many_y = imax2 + 1 - imin2
        call FLOAT_FFTW_INPLACE_DRV_1D( bw_planx(ip), how_many_y, f(first_index), 1, ldx )

      end do   

      RETURN
   END SUBROUTINE cft_b_sp


   SUBROUTINE cft_b_omp_sp ( f, nx, ny, nz, ldx, ldy, ldz, imin2, imax2, imin3, imax3, sgn )

!     driver routine for 3d complex fft's on box grid, parallel (MPI+OpenMP) case
!     fft along z for all xy values 
!     fft along y is done only for the local  z values: i.e. z-planes with imin3 <= nz <= imax3.
!     fft along x is done only for the local yz values: i.e. z-planes with imin3 <= nz <= imax3 
!                                                        and y-planes with imin2 <= ny <= imax2.
!     implemented ONLY for internal fftw, and only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
!     This driver is meant for calls inside parallel OpenMP sections
!
      USE fftw_interfaces
      implicit none
      integer, INTENT(IN) :: nx,ny,nz,ldx,ldy,ldz,imin2,imax2,imin3,imax3,sgn
      complex(sp) :: f(:)

      INTEGER, SAVE :: k, first_index, how_many_y
!$omp threadprivate (k,first_index,how_many_y)

      if ( .NOT. use_single_precision ) then
         CALL fftx_error__('cft_b_omp','single precision driver with double precision initialization',1)
      end if
      if ( -sgn > 0 ) then
         CALL fftx_error__('cft_b_omp','forward transform not implemented',1)
      end if

      IF ( .NOT. C_ASSOCIATED(cft_b_bw_planz) .or. &
           .NOT. C_ASSOCIATED(cft_b_bw_planx) .or. &
           .NOT. C_ASSOCIATED(cft_b_bw_plany) ) THEN
         CALL fftx_error__('cft_b_omp','plan not initialized',1)
      END IF

      !  consistency check

      IF ( ( nx /= cft_b_dims(1) ) .or. ( ny /= cft_b_dims(2) ) .or. ( nz /= cft_b_dims(3) ) ) THEN
         CALL fftx_error__('cft_b_omp', 'dimensions are inconsistent with the existing plan',1) 
      END IF

      !
      !  fft along Z
      !
      call FLOAT_FFTW_INPLACE_DRV_1D( cft_b_bw_planz, ldx*ldy, f(1), ldx*ldy, 1 )

      do k = imin3, imax3
      !
      !  fft along Y
      !
        first_index = (k-1)*ldx*ldy + 1
        call FLOAT_FFTW_INPLACE_DRV_1D( cft_b_bw_plany, nx, f(first_index), ldx, 1 )
      !
      !  fft along X
      !
        first_index = first_index + (imin2-1)*ldx ; how_many_y = imax2 + 1 - imin2
        call FLOAT_FFTW_INPLACE_DRV_1D( cft_b_bw_planx, how_many_y, f(first_index), 1, ldx )
      end do   

     RETURN
   END SUBROUTINE cft_b_omp_sp



!=----------------------------------------------------------------------=!
   END MODULE fft_smallbox
!=----------------------------------------------------------------------=!

