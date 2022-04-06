#if defined(MKL_ILP64)

include "mkl_blas_omp_offload_ilp64.f90"
include "mkl_lapack_omp_offload_ilp64.f90"

module onemkl_blas_omp_offload
  use onemkl_blas_omp_offload_ilp64
endmodule onemkl_blas_omp_offload

module onemkl_lapack_omp_offload
  use onemkl_lapack_omp_offload_ilp64
endmodule onemkl_lapack_omp_offload

#else

include "mkl_blas_omp_offload_lp64.f90"
include "mkl_lapack_omp_offload_lp64.f90"

module onemkl_blas_omp_offload
  use onemkl_blas_omp_offload_lp64
endmodule onemkl_blas_omp_offload

module onemkl_lapack_omp_offload
  use onemkl_lapack_omp_offload_lp64
endmodule onemkl_lapack_omp_offload

#endif
