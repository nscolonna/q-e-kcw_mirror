!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE laxlib_rdiaghg( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
  !!----------------------------------------------------------------------------
  !!
  !! Called by diaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem.
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! real matrices version.
  !! On output both matrix are unchanged.
  !!
  !! LAPACK version - uses both DSYGV and DSYGVX
  !
  USE laxlib_parallel_include
  !
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: m
  !! number of eigenstates to be calculated
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  REAL(DP), INTENT(INOUT) :: h(ldh,n)
  !! matrix to be diagonalized
  REAL(DP), INTENT(INOUT) :: s(ldh,n)
  !! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
  !! eigenvalues
  REAL(DP), INTENT(OUT) :: v(ldh,m)
  !! eigenvectors (column-wise)
  INTEGER,  INTENT(IN)  :: me_bgrp
  !! index of the processor within a band group
  INTEGER,  INTENT(IN)  :: root_bgrp
  !! index of the root processor within a band group
  INTEGER,  INTENT(IN)  :: intra_bgrp_comm
  !! intra band group communicator
  !
  INTEGER               :: lwork, nb, mm, info, i, j
    ! mm = number of calculated eigenvectors
  REAL(DP)              :: abstol
  REAL(DP), PARAMETER   :: one = 1_DP
  REAL(DP), PARAMETER   :: zero = 0_DP
  INTEGER,  ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP), ALLOCATABLE :: work(:), sdiag(:), hdiag(:)
  LOGICAL               :: all_eigenvalues
  INTEGER,  EXTERNAL    :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  CALL start_clock( 'rdiaghg' )
  !
  ! ... only the first processor diagonalize the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
     !
     ! ... save the diagonal of input S (it will be overwritten)
     !
     ALLOCATE( sdiag( n ) )
     DO i = 1, n
        sdiag(i) = s(i,i)
     END DO
     !
     all_eigenvalues = ( m == n )
     !
     ! ... check for optimal block size
     !
     nb = ILAENV( 1, 'DSYTRD', 'U', n, -1, -1, -1 )
     !
     IF ( nb < 5 .OR. nb >= n ) THEN
        !
        lwork = 8*n
        !
     ELSE
        !
        lwork = ( nb + 3 )*n
        !
     END IF
     !
     ALLOCATE( work( lwork ) )
     !
     IF ( all_eigenvalues ) THEN
        !
        ! ... calculate all eigenvalues
        !
        !$omp parallel do
        do i =1, n
           v(1:ldh,i) = h(1:ldh,i)
        end do
        !$omp end parallel do
        !
        CALL DSYGV( 1, 'V', 'U', n, v, ldh, s, ldh, e, work, lwork, info )
        !
     ELSE
        !
        ! ... calculate only m lowest eigenvalues
        !
        ALLOCATE( iwork( 5*n ) )
        ALLOCATE( ifail( n ) )
        !
        ! ... save the diagonal of input H (it will be overwritten)
        !
        ALLOCATE( hdiag( n ) )
        DO i = 1, n
           hdiag(i) = h(i,i)
        END DO
        !
        abstol = 0.D0
       ! abstol = 2.D0*DLAMCH( 'S' )
        !
        CALL DSYGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, iwork, ifail, info )
        !
        DEALLOCATE( ifail )
        DEALLOCATE( iwork )
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
        !$omp parallel do
        DO i = 1, n
           h(i,i) = hdiag(i)
           DO j = i + 1, n
              h(i,j) = h(j,i)
           END DO
           DO j = n + 1, ldh
              h(j,i) = 0.0_DP
           END DO
        END DO
        !$omp end parallel do
        !
        DEALLOCATE( hdiag )
        !
     END IF
     !
     DEALLOCATE( work )
     !
     IF ( info > n ) THEN
        CALL lax_error__( 'rdiaghg', 'S matrix not positive definite', ABS( info ) )
     ELSE IF ( info > 0 ) THEN
        CALL lax_error__( 'rdiaghg', 'eigenvectors failed to converge', ABS( info ) )
     ELSE IF ( info < 0 ) THEN
        CALL lax_error__( 'rdiaghg', 'incorrect call to DSYGV*', ABS( info ) )
     END IF

     ! ... restore input S matrix from saved diagonal and lower triangle
     !
     !$omp parallel do
     DO i = 1, n
        s(i,i) = sdiag(i)
        DO j = i + 1, n
           s(i,j) = s(j,i)
        END DO
        DO j = n + 1, ldh
           s(j,i) = 0.0_DP
        END DO
     END DO
     !$omp end parallel do
     !
     DEALLOCATE( sdiag )
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#if defined __MPI
  CALL MPI_BCAST( e, SIZE(e), MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'rdiaghg', 'error broadcasting array e', ABS( info ))
  CALL MPI_BCAST( v, SIZE(v), MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'rdiaghg', 'error broadcasting array v', ABS( info ))
#endif
  !
  CALL stop_clock( 'rdiaghg' )
  !
  RETURN
  !
END SUBROUTINE laxlib_rdiaghg

!----------------------------------------------------------------------------
#if defined(__OPENMP_GPU)
SUBROUTINE laxlib_rdiaghg_gpu( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
  !----------------------------------------------------------------------------
  ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - uses both DSYGV and DSYGVX
  !
  USE laxlib_parallel_include
  USE onemkl_lapack_omp_offload
  !USE dmr
#define __USE_GLOBAL_BUFFER
#if defined(__USE_GLOBAL_BUFFER)
  USE device_fbuff_m, ONLY : dev=>dev_buf
#else
  USE omp_lib
#endif
  !
  IMPLICIT NONE
  INCLUDE 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculated
    ! leading dimension of h, as declared in the calling pgm unit
  REAL(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
    ! matrix to be diagonalized
    ! overlap matrix
  !
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  REAL(DP), INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
  INTEGER,  INTENT(IN)  :: me_bgrp, root_bgrp, intra_bgrp_comm
  !
  INTEGER               :: lwork, info, i, j, nb, mm, omp_device
  REAL(DP)              :: abstol
  REAL(DP), PARAMETER   :: one = 1_DP
  REAL(DP), PARAMETER   :: zero = 0_DP
  INTEGER,  POINTER, CONTIGUOUS :: ifail(:)
  !
  INTEGER,  POINTER, CONTIGUOUS :: iwork(:)
  REAL(DP), POINTER, CONTIGUOUS :: work(:)
  !
  REAL(DP)                      :: dummy(2)
  !
  ! Temp arrays to save H and S.
  INTEGER                       :: h_meig, itype=1
  REAL(DP), POINTER, CONTIGUOUS :: h_bkp(:,:), s_bkp(:,:)
  character*1 :: jobz = 'V'
  character*1 :: ranged = 'I'
  character*1 :: uplo = 'U'
  !
  !
  CALL start_clock_gpu( 'rdiaghg_gpu' )
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
#if !defined(__USE_GLOBAL_BUFFER)
     !omp_device = omp_get_default_device()
     !call omp_target_alloc_f(fptr_dev=ifail, dimensions=[n],     omp_dev=omp_device)
     !call omp_target_alloc_f(fptr_dev=h_bkp, dimensions=[n,n],   omp_dev=omp_device)
     !call omp_target_alloc_f(fptr_dev=s_bkp, dimensions=[n,n],   omp_dev=omp_device)
     !call omp_target_alloc_f(fptr_dev=iwork, dimensions=[3+5*n], omp_dev=omp_device)
     !$omp allocate allocator(omp_target_device_mem_alloc)
     allocate(ifail(n))
     allocate(h_bkp(n,n))
     allocate(s_bkp(n,n))
     allocate(iwork(3+5*n))
#else
     CALL dev%lock_buffer( ifail,  n, info )
     IF( info /= 0 ) CALL lax_error__( ' rdiaghg_gpu ', ' cannot allocate ifail ', ABS( info ) )
     CALL dev%lock_buffer( h_bkp,  (/ n, n /), info )
     IF( info /= 0 ) CALL lax_error__( ' rdiaghg_gpu ', ' cannot allocate h_bkp ', ABS( info ) )
     CALL dev%lock_buffer( s_bkp,  (/ n, n /), info )
     IF( info /= 0 ) CALL lax_error__( ' rdiaghg_gpu ', ' cannot allocate s_bkp ', ABS( info ) )
     CALL dev%lock_buffer( iwork, 3+5*n, info )
     IF( info /= 0 ) CALL lax_error__( ' rdiaghg_gpu ', ' cannot allocate iwork ', ABS( info ) )
#endif
     !
     !$omp target teams distribute parallel do collapse(2) is_device_ptr(h_bkp, s_bkp)
     DO j=1,n
        DO i=1,n
           h_bkp(i,j) = h(i,j)
           s_bkp(i,j) = s(i,j)
        ENDDO
     ENDDO
     !$omp end target teams distribute parallel do
     !
     ! Workspace query
     !
     !$omp target data map(from:h_meig, dummy, info)
     !$omp target variant dispatch use_device_ptr(h, s, h_meig, e, v, dummy, iwork, ifail, info)
     call dsygvx(itype=itype, jobz=jobz, range=ranged, uplo=uplo, n=n, a=h, lda=ldh, b=s, ldb=ldh, vl=0.0D0, vu=0.0D0, il=1, iu=m, abstol=0.0D0, &
        m=h_meig, w=e, z=v, ldz=ldh, work=dummy, lwork=-1, iwork=iwork, ifail=ifail, info=info)
     !$omp end target variant dispatch
     !$omp end target data
     lwork = nint(dummy(1))
#if ! defined(__USE_GLOBAL_BUFFER)
     !call omp_target_alloc_f(fptr_dev=work,  dimensions=[lwork], omp_dev=omp_device)
     !$omp allocate allocator(omp_target_device_mem_alloc)
     allocate(work(lwork))
#else
     CALL dev%lock_buffer( work, lwork, info )
     IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate work ', ABS( info ) )
#endif
     !$omp target data map(from:h_meig, info)
     !$omp target variant dispatch use_device_ptr(h, s, h_meig, e, v, work, iwork, ifail, info)
     call dsygvx(itype=itype, jobz=jobz, range=ranged, uplo=uplo, n=n, a=h, lda=ldh, b=s, ldb=ldh, vl=0.0D0, vu=0.0D0, il=1, iu=m, abstol=0.0D0, &
        m=h_meig, w=e, z=v, ldz=ldh, work=work, lwork=lwork, iwork=iwork, ifail=ifail, info=info)
     !$omp end target variant dispatch
     !$omp end target data
     !
     IF ( info > n ) THEN
        CALL lax_error__( 'rdiaghg', 'S matrix not positive definite', ABS( info ) )
     ELSE IF ( info > 0 ) THEN
        CALL lax_error__( 'rdiaghg', 'eigenvectors failed to converge', ABS( info ) )
     ELSE IF ( info < 0 ) THEN
        CALL lax_error__( 'rdiaghg', 'incorrect call to DSYGV*', ABS( info ) )
     END IF
     !
     !$omp target teams distribute parallel do collapse(2) is_device_ptr(h_bkp, s_bkp)
     DO j=1,n
        DO i=1,n
           h(i,j) = h_bkp(i,j)
           s(i,j) = s_bkp(i,j)
        ENDDO
     ENDDO
!$omp end target teams distribute parallel do
#if !defined(__USE_GLOBAL_BUFFER)
     !call omp_target_free_f(fptr_dev=ifail)
     !call omp_target_free_f(fptr_dev=h_bkp)
     !call omp_target_free_f(fptr_dev=s_bkp)
     !call omp_target_free_f(fptr_dev=work)
     !call omp_target_free_f(fptr_dev=iwork)
     deallocate(ifail)
     deallocate(h_bkp)
     deallocate(s_bkp)
     deallocate(work)
     deallocate(iwork)
#else
     CALL dev%release_buffer(ifail, info)
     CALL dev%release_buffer(h_bkp, info)
     CALL dev%release_buffer(s_bkp, info)
     CALL dev%release_buffer(work,  info)
     CALL dev%release_buffer(iwork, info)
#endif
  !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#if defined __MPI
#if defined __GPU_MPI
!  info = cudaDeviceSynchronize()
!  IF ( info /= 0 ) &
!        CALL lax_error__( 'rdiaghg', 'error synchronizing device (first)', ABS( info ) )
!  CALL MPI_BCAST( e_d, n, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
!  IF ( info /= 0 ) &
!        CALL lax_error__( 'rdiaghg', 'error broadcasting array e_d', ABS( info ) )
!  CALL MPI_BCAST( v_d, ldh*m, MPI_DOUBLE_COMPLEX, root_bgrp, intra_bgrp_comm, info )
!  IF ( info /= 0 ) &
!        CALL lax_error__( 'rdiaghg', 'error broadcasting array v_d', ABS( info ) )
!  info = cudaDeviceSynchronize() ! this is probably redundant...
!  IF ( info /= 0 ) &
!        CALL lax_error__( 'rdiaghg', 'error synchronizing device (second)', ABS( info ) )
#else
!$omp target update from(e, v)
  CALL MPI_BCAST( e, SIZE(e), MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'rdiaghg', 'error broadcasting array e', ABS( info ))
  CALL MPI_BCAST( v, SIZE(v), MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'rdiaghg', 'error broadcasting array v', ABS( info ))
! Are this copies necessary???
!$omp target update to(e)
!$omp target update to(v)
#endif
#endif
  !
  CALL stop_clock_gpu( 'rdiaghg_gpu' )
  !
  RETURN
  !
END SUBROUTINE laxlib_rdiaghg_gpu
#else
SUBROUTINE laxlib_rdiaghg_gpu( n, m, h_d, s_d, ldh, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm )
  !----------------------------------------------------------------------------
  !!
  !! Called by diaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! real matrices version.
  !! On output both matrix are unchanged.
  !!
  !! GPU VERSION.
  !
#if defined(_OPENMP)
  USE omp_lib
#endif
  USE laxlib_parallel_include
#if defined(__CUDA)
  USE cudafor
  USE cusolverdn
#endif
  !
  ! NB: the flag below can be used to decouple LAXlib from devXlib.
  !     This will make devXlib an optional dependency of LAXlib when
  !     the library will be decoupled from QuantumESPRESSO.
#define __USE_GLOBAL_BUFFER
#if defined(__USE_GLOBAL_BUFFER) && defined(__CUDA)
  USE device_fbuff_m,        ONLY : dev=>dev_buf, pin=>pin_buf
#define VARTYPE POINTER
#else
#define VARTYPE ALLOCATABLE
#endif
  !
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: m
  !! number of eigenstates to be calculated
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  REAL(DP), INTENT(INOUT) :: h_d(ldh,n)
  !! matrix to be diagonalized, allocated on the device
  REAL(DP), INTENT(INOUT) :: s_d(ldh,n)
  !! overlap matrix, allocated on the device
  REAL(DP), INTENT(OUT) :: e_d(n)
  !! eigenvalues, allocated on the device
  REAL(DP), INTENT(OUT) :: v_d(ldh, n)
  !! eigenvectors (column-wise), allocated on the device
  INTEGER,  INTENT(IN)  :: me_bgrp
  !! index of the processor within a band group
  INTEGER,  INTENT(IN)  :: root_bgrp
  !! index of the root processor within a band group
  INTEGER,  INTENT(IN)  :: intra_bgrp_comm
  !! intra band group communicator
#if defined(__CUDA)
    ATTRIBUTES(DEVICE) :: h_d, s_d, e_d, v_d
#endif
  !
  INTEGER               :: lwork, nb, mm, info, i, j
    ! mm = number of calculated eigenvectors
  REAL(DP)              :: abstol
  REAL(DP), PARAMETER   :: one = 1_DP
  REAL(DP), PARAMETER   :: zero = 0_DP
  INTEGER,  ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP), ALLOCATABLE :: work(:), sdiag(:), hdiag(:)
#if defined(__CUDA)
  ATTRIBUTES( PINNED )          :: work, iwork
#endif
  REAL(DP), ALLOCATABLE :: v_h(:,:)
  REAL(DP), ALLOCATABLE :: e_h(:)
#if defined(__CUDA)
  ATTRIBUTES( PINNED )  :: v_h, e_h
#endif
  !
  INTEGER               :: lwork_d, liwork
  REAL(DP), VARTYPE     :: work_d(:)
#if defined(__CUDA)
  ATTRIBUTES( DEVICE )  :: work_d
#endif
  !
  ! Temp arrays to save H and S.
  REAL(DP), ALLOCATABLE :: h_diag_d(:), s_diag_d(:)
#if defined(__CUDA)
  ATTRIBUTES( DEVICE )  :: h_diag_d, s_diag_d
  !
  INTEGER                      :: devInfo_d, h_meig
  ATTRIBUTES( DEVICE )         :: devInfo_d
  TYPE(cusolverDnHandle), SAVE :: cuSolverHandle
  LOGICAL, SAVE                :: cuSolverInitialized = .FALSE.
  REAL(DP), VARTYPE            :: h_bkp_d(:,:), s_bkp_d(:,:)
  ATTRIBUTES( DEVICE )         :: h_bkp_d, s_bkp_d
#endif
#undef VARTYPE
  !
  CALL start_clock_gpu( 'rdiaghg' )
  !
  ! ... only the first processor diagonalize the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
     !
#if defined(__CUDA)
     !
#if ! defined(__USE_GLOBAL_BUFFER)
      ALLOCATE(h_bkp_d(n,n), s_bkp_d(n,n), STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' rdiaghg_gpu ', ' cannot allocate h_bkp_d or s_bkp_d ', ABS( info ) )
#else
      CALL dev%lock_buffer( h_bkp_d,  (/ n, n /), info )
      IF( info /= 0 ) CALL lax_error__( ' rdiaghg_gpu ', ' cannot allocate h_bkp_d ', ABS( info ) )
      CALL dev%lock_buffer( s_bkp_d,  (/ n, n /), info )
      IF( info /= 0 ) CALL lax_error__( ' rdiaghg_gpu ', ' cannot allocate s_bkp_d ', ABS( info ) )
#endif

!$cuf kernel do(2)
      DO j=1,n
         DO i=1,n
            h_bkp_d(i,j) = h_d(i,j)
            s_bkp_d(i,j) = s_d(i,j)
         ENDDO
      ENDDO

#if defined(_OPENMP)
      IF (omp_get_num_threads() > 1) CALL lax_error__( ' rdiaghg_gpu ', 'rdiaghg_gpu is not thread-safe',  ABS( info ) )
#endif
      IF ( .NOT. cuSolverInitialized ) THEN
         info = cusolverDnCreate(cuSolverHandle)
         IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' rdiaghg_gpu ', ' cusolverDnCreate failed ', ABS( info ) )
         cuSolverInitialized = .TRUE.
      ENDIF
      info = cusolverDnDsygvdx_bufferSize(cuSolverHandle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, &
                                                         CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                                               n, h_d, ldh, s_d, ldh, 0.D0, 0.D0, 1, m, h_meig, e_d, lwork_d)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' rdiaghg_gpu ', ' cusolverDnDsygvdx_bufferSize failed ', ABS( info ) )

#if ! defined(__USE_GLOBAL_BUFFER)
      ALLOCATE(work_d(1*lwork_d), STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' rdiaghg_gpu ', ' cannot allocate work_d ', ABS( info ) )
#else
      CALL dev%lock_buffer( work_d,  lwork_d, info )
      IF( info /= 0 ) CALL lax_error__( ' rdiaghg_gpu ', ' allocate work_d ', ABS( info ) )
#endif
      info = cusolverDnDsygvdx(cuSolverHandle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, &
                                               CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                                               n, h_d, ldh, s_d, ldh, 0.D0, 0.D0, 1, m, h_meig,&
                                               e_d, work_d, lwork_d, devInfo_d)
    IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' rdiaghg_gpu ', ' cusolverDnDsygvdx failed ', ABS( info ) )
!$cuf kernel do(2)
      DO j=1,n
         DO i=1,n
            IF(j <= m) v_d(i,j) = h_d(i,j)
            h_d(i,j) = h_bkp_d(i,j)
            s_d(i,j) = s_bkp_d(i,j)
         ENDDO
      ENDDO
      !
      !
      ! Do not destroy the handle to save the (re)creation time on each call.
      !
      ! info = cusolverDnDestroy(cuSolverHandle)
      ! IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' rdiaghg_gpu ', ' cusolverDnDestroy failed ', ABS( info ) )
      !
#if ! defined(__USE_GLOBAL_BUFFER)
      DEALLOCATE(work_d)
      DEALLOCATE(h_bkp_d, s_bkp_d)
#else
      CALL dev%release_buffer( work_d,  info )
      CALL dev%release_buffer( h_bkp_d, info )
      CALL dev%release_buffer( s_bkp_d, info )
#endif

#else
     CALL lax_error__( 'cdiaghg', 'Called GPU eigensolver without GPU support', 1 )
#endif
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#if defined __MPI
#if defined __GPU_MPI
  info = cudaDeviceSynchronize()
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error synchronizing device (first)', ABS( info ) )
  CALL MPI_BCAST( e_d(1), n, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'rdiaghg', 'error broadcasting array e_d', ABS( info ))
  CALL MPI_BCAST( v_d(1,1), ldh*m, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'rdiaghg', 'error broadcasting array v_d', ABS( info ))
  info = cudaDeviceSynchronize() ! this is probably redundant...
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error synchronizing device (second)', ABS( info ) )
#else
  ALLOCATE(e_h(n), v_h(ldh,m))
  e_h(1:n) = e_d(1:n)
  v_h(1:ldh, 1:m) = v_d(1:ldh, 1:m)
  CALL MPI_BCAST( e_h, n, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array e_d', ABS( info ) )
  CALL MPI_BCAST( v_h, ldh*m, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array v_d', ABS( info ) )
  e_d(1:n) = e_h(1:n)
  v_d(1:ldh, 1:m) = v_h(1:ldh, 1:m)
  DEALLOCATE(e_h, v_h)
#endif
#endif
  !
  CALL stop_clock_gpu( 'rdiaghg' )
  !
  RETURN
  !
END SUBROUTINE laxlib_rdiaghg_gpu
#endif
!
!----------------------------------------------------------------------------
SUBROUTINE laxlib_prdiaghg( n, h, s, ldh, e, v, idesc )
  !----------------------------------------------------------------------------
  !
  !! Called by pdiaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem.
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! real matrices version.
  !! On output both matrix are unchanged.
  !!
  !! Parallel version with full data distribution
  !!
  !
  USE laxlib_parallel_include
  USE laxlib_descriptor, ONLY : la_descriptor, laxlib_intarray_to_desc
  USE laxlib_processors_grid, ONLY : ortho_parent_comm
#if defined __SCALAPACK
  USE laxlib_processors_grid, ONLY : ortho_cntx, np_ortho, me_ortho, ortho_comm
  USE dspev_module,      ONLY : pdsyevd_drv
#endif
  !
  IMPLICIT NONE
  !
  include 'laxlib_kinds.fh'
  include 'laxlib_param.fh'
  include 'laxlib_low.fh'
  include 'laxlib_mid.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  REAL(DP), INTENT(INOUT) :: h(ldh,ldh)
  !! matrix to be diagonalized
  REAL(DP), INTENT(INOUT) :: s(ldh,ldh)
  !! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
  !! eigenvalues
  REAL(DP), INTENT(OUT) :: v(ldh,ldh)
  !! eigenvectors (column-wise)
  INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
  !! laxlib descriptor
  INTEGER, PARAMETER    :: root = 0
  INTEGER               :: nx, info
    ! local block size
  REAL(DP), PARAMETER   :: one = 1_DP
  REAL(DP), PARAMETER   :: zero = 0_DP
  REAL(DP), ALLOCATABLE :: hh(:,:)
  REAL(DP), ALLOCATABLE :: ss(:,:)
  TYPE(la_descriptor) :: desc
#if defined(__SCALAPACK)
  INTEGER     :: desch( 16 )
#endif
  INTEGER               :: i
  !
  CALL start_clock( 'rdiaghg' )
  !
  CALL laxlib_intarray_to_desc(desc,idesc)
  !
  IF( desc%active_node > 0 ) THEN
     !
     nx   = desc%nrcx
     !
     IF( nx /= ldh ) &
        CALL lax_error__(" prdiaghg ", " inconsistent leading dimension ", ldh )
     !
     ALLOCATE( hh( nx, nx ) )
     ALLOCATE( ss( nx, nx ) )
     !
     !$omp parallel do
     do i=1,nx
        hh(1:nx,i) = h(1:nx,i)
        ss(1:nx,i) = s(1:nx,i)
     end do
     !$omp end parallel do
     !
  END IF
  !
  CALL start_clock( 'rdiaghg:choldc' )
  !
  ! ... Cholesky decomposition of s ( L is stored in s )
  !
  IF( desc%active_node > 0 ) THEN
     !
#if defined(__SCALAPACK)
     CALL descinit( desch, n, n, desc%nrcx, desc%nrcx, 0, 0, ortho_cntx, SIZE( hh, 1 ) , info )

     IF( info /= 0 ) CALL lax_error__( ' rdiaghg ', ' descinit ', ABS( info ) )
#endif
     !
#if defined(__SCALAPACK)
     CALL PDPOTRF( 'L', n, ss, 1, 1, desch, info )
     IF( info /= 0 ) CALL lax_error__( ' rdiaghg ', ' problems computing cholesky ', ABS( info ) )
#else
     CALL laxlib_pdpotrf( ss, nx, n, idesc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'rdiaghg:choldc' )
  !
  ! ... L is inverted ( s = L^-1 )
  !
  CALL start_clock( 'rdiaghg:inversion' )
  !
  IF( desc%active_node > 0 ) THEN
     !
#if defined(__SCALAPACK)
     !
     CALL sqr_setmat( 'U', n, zero, ss, size(ss,1), idesc )

     CALL PDTRTRI( 'L', 'N', n, ss, 1, 1, desch, info )
     !
     IF( info /= 0 ) CALL lax_error__( ' rdiaghg ', ' problems computing inverse ', ABS( info ) )
#else
     CALL laxlib_pdtrtri ( ss, nx, n, idesc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'rdiaghg:inversion' )
  !
  ! ... v = L^-1*H
  !
  CALL start_clock( 'rdiaghg:paragemm' )
  !
  IF( desc%active_node > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, idesc )
     !
  END IF
  !
  ! ... h = ( L^-1*H )*(L^-1)^T
  !
  IF( desc%active_node > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'T', n, ONE, v, nx, ss, nx, ZERO, hh, nx, idesc )
     !
  END IF
  !
  CALL stop_clock( 'rdiaghg:paragemm' )
  !
  IF ( desc%active_node > 0 ) THEN
     !
     !  Compute local dimension of the cyclically distributed matrix
     !
#if defined(__SCALAPACK)
     CALL pdsyevd_drv( .true., n, desc%nrcx, hh, SIZE(hh,1), e, ortho_cntx, ortho_comm )
#else
     CALL laxlib_pdsyevd( .true., n, idesc, hh, SIZE(hh,1), e )
#endif
     !
  END IF
  !
  ! ... v = (L^T)^-1 v
  !
  CALL start_clock( 'rdiaghg:paragemm' )
  !
  IF ( desc%active_node > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'T', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, idesc )
     !
     DEALLOCATE( ss )
     DEALLOCATE( hh )
     !
  END IF
  !
#if defined __MPI
  CALL MPI_BCAST( e, SIZE(e), MPI_DOUBLE_PRECISION, root, ortho_parent_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'prdiaghg', 'error broadcasting array e', ABS( info ))
#endif
  !
  CALL stop_clock( 'rdiaghg:paragemm' )
  !
  CALL stop_clock( 'rdiaghg' )
  !
  RETURN
  !
END SUBROUTINE laxlib_prdiaghg
