#if defined(__CUDA) || defined(__OPENMP_GPU)
program test_diaghg_gpu

    USE laxlib_parallel_include
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : mp_world_start, mp_world_end, mpime, &
                              root, nproc, world_comm
    USE mp_bands_util, ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
    USE tester
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    !
    CALL test%init()

#if defined(__MPI)
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    !
    me_bgrp = mpime; root_bgrp=root; intra_bgrp_comm=world_comm
    !
    CALL real_1(test)
    !
    CALL complex_1(test)
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
  CONTAINS
  !
  SUBROUTINE real_1(test)
    USE LAXlib
#if !defined(__OPENMP_GPU)
    USE cudafor
#endif
    implicit none
    !
    TYPE(tester_t) :: test
    ! variables on host
    real(DP) :: h(2,2)
    real(DP) :: s(2,2)
    real(DP) :: e(2)
    real(DP) :: v(2,2)
#if !defined(__OPENMP_GPU)
    ! variables on device
    real(DP), device :: h_d(2,2)
    real(DP), device :: s_d(2,2)
    real(DP), device :: e_d(2)
    real(DP), device :: v_d(2,2)
#endif

    h      = 0.d0
    h(1,1) = 1.d0
    h(2,2) = 1.d0
    !
    s      = 0.d0
    s(1,1) = 1.d0
    s(2,2) = 1.d0
    !
    v = 0.d0
    e = 0.d0
    !
#if defined(__OPENMP_GPU)
    !$omp target data map(tofrom:h, s, e, v)
    CALL diaghg(  2, 2, h, s, 2, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
    !$omp end target data
#else
    h_d = h
    !
    s_d = s
    !
    CALL diaghg(  2, 2, h_d, s_d, 2, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
    !
    v = v_d
    e = e_d
    h = h_d
    s = s_d
#endif
    !
    CALL test%assert_close( e, [1.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(v, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(s, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    !
#if defined(__OPENMP_GPU)
    v = 0.d0
    e = 0.d0
    !
    !$omp target data map(tofrom:h, s, e, v)
    CALL diaghg(  2, 2, h, s, 2, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .false., .true. )
    !$omp end target data
#else
    v_d = 0.d0
    e_d = 0.d0
    !
    CALL diaghg(  2, 2, h_d, s_d, 2, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, .true. )
    !
    v = v_d
    e = e_d
    h = h_d
    s = s_d
#endif
    !
    CALL test%assert_close( e, [1.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(v, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(s, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    !
  END SUBROUTINE real_1
  !
  SUBROUTINE complex_1(test)
    USE LAXlib
#if !defined(__OPENMP_GPU)
    USE cudafor
#endif
    implicit none
    !
    TYPE(tester_t) :: test
    ! variables on host
    complex(DP) :: h(2,2)
    complex(DP) :: s(2,2)
    complex(DP) :: v(2,2)
    real(DP)    :: e(2)
#if !defined(__OPENMP_GPU)
    ! variables on device
    complex(DP), device :: h_d(2,2)
    complex(DP), device :: s_d(2,2)
    complex(DP), device :: v_d(2,2)
    real(DP),    device :: e_d(2)
#endif

    complex(DP) :: s_save(2,2)
    complex(DP) :: h_save(2,2)

    !
    h = 0.d0
    h(1,1) = (1.d0,  0.d0)
    h(1,2) = (0.d0, -2.d0)
    h(2,1) = ( 0.d0, 2.d0)
    h(2,2) = ( 5.d0, 0.d0)
    s = 0.d0
    s(1,1) = (1.d0, 0.d0)
    s(2,2) = (1.d0, 0.d0)
    v = (0.d0, 0.d0)
    e = 0.d0
    !
    ! save for later comparison
    h_save = h
    s_save = s
    !
#if defined(__OPENMP_GPU)
    !$omp target data map(tofrom:h, s, e, v)
    CALL diaghg(  2, 2, h, s, 2, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
    !$omp end target data
#else
    ! Update device
    h_d = h
    s_d = s
    v_d = v; e_d = e
    !
    CALL diaghg(  2, 2, h_d, s_d, 2, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
    v = v_d
    e = e_d
    h = h_d
    s = s_d
#endif
    !
    ! 0.1715728752538099, 5.82842712474619
    CALL test%assert_close( e, [0.1715728752538099d0,  5.82842712474619d0] )
    CALL test%assert_close( v(:,1), [( 0.d0, -0.9238795325112867d0), (-0.3826834323650898d0, 0.d0)] )
    CALL test%assert_close( v(:,2), [( 0.d0, -0.3826834323650898d0), ( 0.9238795325112867d0, 0.d0)] )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]))
    CALL test%assert_close( RESHAPE(s, [4]), RESHAPE(s_save, [4]))
    !
    v= (0.d0, 0.d0)
    e= 0.d0
    !
#if defined(__OPENMP_GPU)
    !$omp target data map(tofrom:h, s, e, v)
    CALL diaghg(  2, 2, h, s, 2, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .false., .true. )
    !$omp end target data
#else
    ! Update device
    h_d = h
    s_d = s
    v_d = v; e_d = e
    CALL diaghg(  2, 2, h_d, s_d, 2, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, .true. )
    v = v_d
    e = e_d
    h = h_d
    s = s_d
#endif
    !
    CALL test%assert_close( e, [0.1715728752538099d0,  5.82842712474619d0] )
    CALL test%assert_close( v(:,1), [( 0.d0, -0.9238795325112867d0), (-0.3826834323650898d0, 0.d0)] )
    CALL test%assert_close( v(:,2), [( 0.d0, -0.3826834323650898d0), ( 0.9238795325112867d0, 0.d0)] )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]))
    CALL test%assert_close( RESHAPE(s, [4]), RESHAPE(s_save, [4]))
    !
  END SUBROUTINE complex_1
end program test_diaghg_gpu
#else
program test_diaghg_gpu
end program test_diaghg_gpu
#endif
