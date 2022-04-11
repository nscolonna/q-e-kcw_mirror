!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE cegterg_gpu( h_psi_gpu, s_psi_gpu, uspp, g_psi_gpu, &
                    npw, npwx, nvec, nvecx, npol, evc_d, ethr, &
                    e_d, btype, notcnv, lrot, dav_iter, nhpsi )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  !
#if defined(__CUDA)
  use cudafor
  use cublas
#elif defined(__OPENMP_GPU)
  use omp_lib
  use onemkl_blas_omp_offload
#endif
  USE LAXlib,        ONLY : diaghg
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
                            nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum, mp_gather, mp_bcast, mp_size,&
#if defined(__OPENMP_GPU)
                            mp_sum_mapped, mp_bcast_mapped, &
#endif
                            mp_type_create_column_section, mp_type_free
#if defined(__OPENMP_GPU)
  USE device_fbuff_m,      ONLY : gbuf => pin_buf
#endif
  USE device_memcpy_m, ONLY : dev_memcpy, dev_memset
  !
  IMPLICIT NONE
  !
  !include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! umber of spin polarizations
  COMPLEX(DP), INTENT(INOUT) :: evc_d(npwx*npol,nvec)
    !  evc contains the  refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence :
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e_d(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
#if defined(__CUDA)
  attributes(DEVICE) :: evc_d, e_d
#endif
  INTEGER, INTENT(OUT) :: nhpsi
    ! total number of indivitual hpsi
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! adapted npw and npwx
    ! do-loop counters
  INTEGER :: n_start, n_end, my_n
  INTEGER :: column_section_type
    ! defines a column section for communication
  INTEGER :: ierr
  COMPLEX(DP), ALLOCATABLE :: hc_d(:,:), sc_d(:,:), vc_d(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
#if !defined(__OPENMP_GPU)
  REAL(DP), ALLOCATABLE :: ew_d(:), ew_host(:)
  REAL(DP), ALLOCATABLE :: e_host(:)
#else
  REAL(DP), ALLOCATABLE :: ew(:)
  INTEGER :: omp_host, omp_device
#endif
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi_d(:,:), hpsi_d(:,:), spsi_d(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr
    ! threshold for empty bands
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
    ! receive counts and memory offsets
  COMPLEX(DP), POINTER  :: pinned_buffer(:,:)
    ! auxiliary variable for performing MPI operation and overcome CUDAFortran limitations
    ! For OpenMP offload pinned memory is not yet available
    ! auxiliary variables for performing dot product
  INTEGER :: i,j,k, ipol
    !
#if defined(__CUDA)
  attributes(DEVICE) :: hc_d, sc_d, vc_d, ew_d, psi_d, hpsi_d, spsi_d
#endif
  !
  REAL(DP), EXTERNAL :: KSddot
  !
  EXTERNAL  h_psi_gpu,    s_psi_gpu,    g_psi_gpu
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  nhpsi = 0
  CALL start_clock( 'cegterg' ); !write(*,*) 'start cegterg' ; FLUSH(6)
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  !$omp allocate allocator(omp_target_device_mem_alloc)
  ALLOCATE(  psi_d( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate psi ', ABS(ierr) )
  !$omp allocate allocator(omp_target_device_mem_alloc)
  ALLOCATE( hpsi_d( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !$omp allocate allocator(omp_target_device_mem_alloc)
     ALLOCATE( spsi_d( npwx*npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  !$omp allocate allocator(omp_target_device_mem_alloc)
  ALLOCATE( sc_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate sc_d ', ABS(ierr) )
  !$omp allocate allocator(omp_target_device_mem_alloc)
  ALLOCATE( hc_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hc_d ', ABS(ierr) )
  !$omp allocate allocator(omp_target_device_mem_alloc)
  ALLOCATE( vc_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc_d ', ABS(ierr) )
#if defined(__OPENMP_GPU)
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew ', ABS(ierr) )
  !$omp target enter data map(alloc:ew)
#else
  ALLOCATE( ew_d( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew_d ', ABS(ierr) )
  ALLOCATE( ew_host( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew_host ', ABS(ierr) )
  ALLOCATE( e_host( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate e_host ', ABS(ierr) )
#endif
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate conv ', ABS(ierr) )
  ALLOCATE( recv_counts(mp_size(inter_bgrp_comm)), displs(mp_size(inter_bgrp_comm)) )
  !
  ! This buffer is used to perform MPI calls with non-contiguous slices.
  ! In order to limit the number of allocated buffers, a rather large,
  ! but hopefully 'repetitive' size is selected (as of today buffers are
  ! selected according to the leading dimension(s) )
  !
#if defined(__OPENMP_GPU)
  allocate(pinned_buffer(nvecx, nvecx), stat=ierr)
#else
  CALL gbuf%lock_buffer(pinned_buffer, (/nvecx, nvecx/), ierr)
#endif
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
#if defined(__OPENMP_GPU)
  !$omp target teams loop collapse(2)
  do j=1,nvec
     do i=1, npwx*npol
        psi_d(i,j) = evc_d(i,j)
     enddo
  enddo
  !$omp end target teams loop
#else
  CALL dev_memcpy(psi_d, evc_d, (/ 1 , npwx*npol /), 1, &
                                (/ 1 , nvec /), 1)
#endif

  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_gpu( npwx, npw, nvec, psi_d, hpsi_d ) ; nhpsi = nhpsi + nvec
  !
  ! ... spsi contains s times the basis vectors
  !
  IF ( uspp ) CALL s_psi_gpu( npwx, npw, nvec, psi_d, spsi_d )
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced
  ! ... space vc contains the eigenvectors of hc
  !
  CALL start_clock( 'cegterg:init' )
  !
  CALL divide_all(inter_bgrp_comm,nbase,n_start,n_end,recv_counts,displs)
  CALL mp_type_create_column_section(sc_d(1,1), 0, nbase, nvecx, column_section_type)
  my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
  !
  if (n_start .le. n_end) then
     !$omp target variant dispatch use_device_ptr(psi_d, hpsi_d, hc_d)
     CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi_d, kdmx, hpsi_d(1,n_start), &
                 kdmx, ZERO, hc_d(1,n_start), nvecx )
     !$omp end target variant dispatch
  endif
  !
  if (n_start .le. n_end) then
     !
     !pinned_buffer(1:nbase, n_start:n_end) = hc_d( 1:nbase, n_start:n_end )
     !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, hc_d( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
#if defined(__OPENMP_GPU)
     !$omp target teams distribute parallel do collapse(2) map(from:pinned_buffer)
     do j=n_start,n_end
        do i=1,nbase
           pinned_buffer(i,j) = hc_d(i,j)
        enddo
     enddo
     !$omp end target teams distribute parallel do
#else
     CALL dev_memcpy( pinned_buffer, hc_d, (/ 1, nbase /), 1, (/ n_start, n_end /), 1 )
#endif
     CALL mp_sum( pinned_buffer(1:nbase, n_start:n_end), intra_bgrp_comm )
     !hc_d( 1:nbase, n_start:n_end ) = pinned_buffer(1:nbase, n_start:n_end)
     !ierr = cudaMemcpy2D( hc_d(1, n_start) , nvecx, pinned_buffer( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
#if defined(__OPENMP_GPU)
     !$omp target teams distribute parallel do collapse(2) map(to:pinned_buffer)
     do j=n_start,n_end
        do i=1,nbase
           hc_d(i,j) = pinned_buffer(i,j)
        enddo
     enddo
     !$omp end target teams distribute parallel do
#else
     CALL dev_memcpy( hc_d, pinned_buffer, (/ 1, nbase /), 1, (/ n_start, n_end /), 1 )
#endif
     !
  end if
  !$omp dispatch
  CALL mp_gather( hc_d, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
  !
  IF ( uspp ) THEN
     !
     if (n_start .le. n_end) then
        !$omp target variant dispatch use_device_ptr(psi_d, spsi_d, sc_d)
        CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi_d, kdmx, spsi_d(1,n_start), kdmx, &
                    ZERO, sc_d(1,n_start), nvecx )
        !$omp end target variant dispatch
     endif
     !
  ELSE
     !
     if (n_start .le. n_end) then
        !$omp target variant dispatch use_device_ptr(psi_d, sc_d)
        CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi_d, kdmx, psi_d(1,n_start), kdmx, &
                    ZERO, sc_d(1,n_start), nvecx )
        !$omp end target variant dispatch
     endif
     !
  END IF
  !
  if ((n_start .le. n_end) .and. (mp_size(intra_bgrp_comm) > 1 )) then
     !pinned_buffer(1:nbase, n_start:n_end) = sc_d( 1:nbase, n_start:n_end )
     !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, sc_d( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
#if defined(__OPENMP_GPU)
     !$omp target teams distribute parallel do collapse(2) map(from:pinned_buffer)
     do j=n_start,n_end
        do i=1,nbase
           pinned_buffer(i,j) = sc_d(i,j)
        enddo
     enddo
     !$omp end target teams distribute parallel do
#else
     CALL dev_memcpy( pinned_buffer, sc_d, (/ 1, nbase /), 1, (/ n_start, n_end /), 1 )
#endif
     CALL mp_sum( pinned_buffer( 1:nbase, n_start:n_end ), intra_bgrp_comm )
     !sc_d( 1:nbase, n_start:n_end ) = pinned_buffer(1:nbase, n_start:n_end)
     !ierr = cudaMemcpy2D( sc_d(1, n_start) , nvecx, pinned_buffer( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
#if defined(__OPENMP_GPU)
     !$omp target teams distribute parallel do collapse(2) map(to:pinned_buffer)
     do j=n_start,n_end
        do i=1,nbase
           sc_d(i,j) = pinned_buffer(i,j)
        enddo
     enddo
     !$omp end target teams distribute parallel do
#else
     CALL dev_memcpy( sc_d, pinned_buffer, (/ 1, nbase /), 1, (/ n_start, n_end /), 1 )
#endif
  end if
  !$omp dispatch
  CALL mp_gather( sc_d, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
  !
  CALL mp_type_free( column_section_type )
  !
!$cuf kernel do
!$omp target teams distribute parallel do
  DO n = 1, nbase
     !
     ! ... the diagonal of hc and sc must be strictly real
     !
     hc_d(n,n) = CMPLX( REAL( hc_d(n,n) ), 0.D0 ,kind=DP)
     sc_d(n,n) = CMPLX( REAL( sc_d(n,n) ), 0.D0 ,kind=DP)
     !
     DO m = n + 1, nbase
        !
        hc_d(n,m) = CONJG( hc_d(m,n) )
        sc_d(n,m) = CONJG( sc_d(m,n) )
        !
     END DO
     !
  END DO
  !
  CALL stop_clock( 'cegterg:init' )
  !
  IF ( lrot ) THEN
     !
     CALL dev_memset(vc_d, ZERO, (/1, nbase/), 1, (/1, nbase/), 1)
     !
!$cuf kernel do(1) <<<*,*>>>
!$omp target teams distribute parallel do
     DO n = 1, nbase
        !
        e_d(n) = REAL( hc_d(n,n) )
        !
        vc_d(n,n) = ONE
        !
     END DO
     !
#if defined(__OPENMP_GPU)
     CALL mp_bcast_mapped( e_d, root_bgrp_id, inter_bgrp_comm )
#else
     CALL mp_bcast( e_d, root_bgrp_id, inter_bgrp_comm )
#endif
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
#if defined(__OPENMP_GPU)
        CALL diaghg( nbase, nvec, hc_d, sc_d, nvecx, ew, vc_d, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
#else
        CALL diaghg( nbase, nvec, hc_d, sc_d, nvecx, ew_d, vc_d, me_bgrp, root_bgrp, intra_bgrp_comm )
#endif
     END IF
     IF( nbgrp > 1 ) THEN
        !$omp dispatch
        CALL mp_bcast( vc_d, root_bgrp_id, inter_bgrp_comm )
#if defined(__OPENMP_GPU)
        CALL mp_bcast_mapped( ew, root_bgrp_id, inter_bgrp_comm )
#else
        CALL mp_bcast( ew_d, root_bgrp_id, inter_bgrp_comm )
#endif
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !
#if defined(__OPENMP_GPU)
     !$omp target teams loop
     do i=1, nvec
        e_d(i) = ew(i)
     enddo
     !$omp end target teams loop
#else
     CALL dev_memcpy (e_d, ew_d, (/ 1, nvec /), 1 )
#endif
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !
     CALL start_clock( 'cegterg:update' )
     !
     !  ======== FROM HERE =====
     !np = 0
     !
     !DO n = 1, nvec
     !   !
     !   IF ( .NOT. conv(n) ) THEN
     !      !
     !      ! ... this root not yet converged ...
     !      !
     !      np = np + 1
     !      !
     !      ! ... reorder eigenvectors so that coefficients for unconverged
     !      ! ... roots come first. This allows to use quick matrix-matrix
     !      ! ... multiplications to set a new basis vector (see below)
     !      !
     !      IF ( np /= n ) vc_d(:,np) = vc_d(:,n)
     !      !
     !      ! ... for use in g_psi
     !      !
     !      ew_d(nbase+np) = e_d(n)
     !      !
     !   END IF
     !   !
     !END DO
     ! ========= TO HERE, REPLACED BY =======

#if defined(__OPENMP_GPU)
     CALL reorder_evals_cevecs(nbase, nvec, nvecx, conv, e_d, ew, vc_d)
#else
     CALL reorder_evals_cevecs(nbase, nvec, nvecx, conv, e_d, ew_d, vc_d)
#endif
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
     IF ( uspp ) THEN
        !
        if (n_start .le. n_end) then
           !$omp target variant dispatch use_device_ptr(spsi_d, vc_d, psi_d)
           CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, spsi_d(1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                       ZERO, psi_d(1,nb1), kdmx )
           !$omp end target variant dispatch
         endif
        !
     ELSE
        !
        if (n_start .le. n_end) then
           !$omp target variant dispatch use_device_ptr(psi_d, vc_d, psi_d)
           CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, psi_d(1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                       ZERO, psi_d(1,nb1), kdmx )
           !$omp end target variant dispatch
         endif
        !
     END IF
! NB: must not call mp_sum over inter_bgrp_comm here because it is done later to the full correction
     !
     !$cuf kernel do(3) <<<*,*>>>
     !$omp target teams distribute parallel do collapse(3)
     DO np=1,notcnv
        DO ipol = 1, npol
           DO k=1,npwx
#if defined(__OPENMP_GPU)
             psi_d(k + (ipol-1)*npwx,nbase+np) = - ew  (nbase+np)*psi_d(k + (ipol-1)*npwx,nbase+np)
#else
             psi_d(k + (ipol-1)*npwx,nbase+np) = - ew_d(nbase+np)*psi_d(k + (ipol-1)*npwx,nbase+np)
#endif
           END DO
        END DO
     END DO
     !
     if (n_start .le. n_end) then
        !$omp target variant dispatch use_device_ptr(hpsi_d, vc_d, psi_d)
        CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, hpsi_d(1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                    ONE, psi_d(1,nb1), kdmx )
        !$omp end target variant dispatch
     endif
     !
     !$omp dispatch
     CALL mp_sum( psi_d(:,nb1:nbase+notcnv), inter_bgrp_comm )
     !
     ! clean up garbage if there is any
     IF (npw < npwx) then
        CALL dev_memset(psi_d, ZERO, [npw+1,npwx], 1, [nb1, nbase+notcnv])
     endif
     IF (npol == 2)  then
        CALL dev_memset(psi_d, ZERO, [npwx+npw+1,2*npwx], 1, [nb1, nbase+notcnv])
     endif
     !
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration
     !
#if defined(__OPENMP_GPU)
     CALL g_psi_gpu( npwx, npw, notcnv, npol, psi_d(1,nb1), ew(nb1) )
#else
     CALL g_psi_gpu( npwx, npw, notcnv, npol, psi_d(1,nb1), ew_d(nb1) )
#endif
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     !!! == OPTIMIZE HERE ==
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
#if defined(__OPENMP_GPU)
           ew(n) = KSDdot( 2*npw, psi_d(1,nbn), 1, psi_d(1,nbn), 1 )
#else
           ew_host(n) = KSDdot( 2*npw, psi_d(1,nbn), 1, psi_d(1,nbn), 1 )
#endif
           !
        ELSE
           !
#if defined(__OPENMP_GPU)
           ew(n) = KSDdot( 2*npw, psi_d(1,nbn), 1, psi_d(1,nbn), 1 ) + &
                   KSDdot( 2*npw, psi_d(npwx+1,nbn), 1, psi_d(npwx+1,nbn), 1 )
#else
           ew_host(n) = KSDdot( 2*npw, psi_d(1,nbn), 1, psi_d(1,nbn), 1 ) + &
                        KSDdot( 2*npw, psi_d(npwx+1,nbn), 1, psi_d(npwx+1,nbn), 1 )
#endif
           !
        END IF
        !
     END DO
     !
#if defined(__OPENMP_GPU)
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !$omp target update to(ew(1:notcnv))
#else
     CALL mp_sum( ew_host( 1:notcnv ), intra_bgrp_comm )
     ew_d(1:notcnv) = ew_host(1:notcnv)
#endif
     !
!$cuf kernel do(3) <<<*,*>>>
!$omp target teams distribute parallel do collapse(3)
     DO i = 1,notcnv
        DO ipol = 1,npol
           DO k=1,npw
#if defined(__OPENMP_GPU)
             psi_d(k + (ipol-1)*npwx,nbase+i) = psi_d(k+(ipol-1)*npwx,nbase+i)/SQRT( ew(i) )
#else
             psi_d(k + (ipol-1)*npwx,nbase+i) = psi_d(k+(ipol-1)*npwx,nbase+i)/SQRT( ew_d(i) )
#endif
           END DO
        END DO
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi_gpu( npwx, npw, notcnv, psi_d(:,nb1), hpsi_d(:,nb1) ) ; nhpsi = nhpsi + notcnv
     !
     IF ( uspp ) CALL s_psi_gpu( npwx, npw, notcnv, psi_d(1,nb1), spsi_d(1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:overlap' )
     !
     CALL divide_all(inter_bgrp_comm,nbase+notcnv,n_start,n_end,recv_counts,displs)
     CALL mp_type_create_column_section(sc_d(1,1), nbase, notcnv, nvecx, column_section_type)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     !
     !$omp target variant dispatch use_device_ptr(hpsi_d, psi_d, hc_d)
     CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, hpsi_d(1,nb1), kdmx, psi_d(1,n_start), kdmx, &
                 ZERO, hc_d(nb1,n_start), nvecx )
     !$omp end target variant dispatch
     !
     if ((n_start .le. n_end) .and. (mp_size(intra_bgrp_comm) > 1 )) then
        !pinned_buffer(nb1:nbase+notcnv, n_start:n_end) = hc_d( nb1:nbase+notcnv, n_start:n_end )
        !ierr = cudaMemcpy2D( pinned_buffer(nb1, n_start) , nvecx, hc_d( nb1, n_start ), nvecx, notcnv, n_end-n_start+1 )
#if defined(__OPENMP_GPU)
        !$omp target teams distribute parallel do collapse(2) map(from:pinned_buffer)
        do j=n_start, n_end
           do i=nb1, nbase+notcnv
              pinned_buffer(i,j) = hc_d(i,j)
           enddo
        enddo
        !$omp end target teams distribute parallel do
#else
        CALL dev_memcpy( pinned_buffer, hc_d, (/ nb1, nbase + notcnv /), 1, (/ n_start, n_end /), 1 )
#endif
        CALL mp_sum( pinned_buffer( nb1:nbase+notcnv, n_start:n_end ), intra_bgrp_comm )
        !hc_d( nb1:nbase+notcnv, n_start:n_end ) = pinned_buffer(nb1:nbase+notcnv, n_start:)
        !ierr = cudaMemcpy2D(  hc_d( nb1, n_start ), nvecx, pinned_buffer(nb1,n_start), nvecx, notcnv, n_end-n_start+1 )
#if defined(__OPENMP_GPU)
        !$omp target teams distribute parallel do collapse(2) map(to:pinned_buffer)
        do j=n_start, n_end
           do i=nb1, nbase+notcnv
              hc_d(i,j) = pinned_buffer(i,j)
           enddo
        enddo
        !$omp end target teams distribute parallel do
#else
        CALL dev_memcpy( hc_d, pinned_buffer, (/ nb1, nbase + notcnv /), 1, (/ n_start, n_end /), 1 )
#endif
     end if
     !$omp dispatch
     CALL mp_gather( hc_d, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
     !
     CALL divide(inter_bgrp_comm,nbase+notcnv,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     IF ( uspp ) THEN
        !
        !$omp target variant dispatch use_device_ptr(spsi_d, psi_d, sc_d)
        CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, spsi_d(1,nb1), kdmx, psi_d(1,n_start), kdmx, &
                    ZERO, sc_d(nb1,n_start), nvecx )
        !$omp end target variant dispatch
        !
     ELSE
        !
        !$omp target variant dispatch use_device_ptr(psi_d, sc_d)
        CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, psi_d(1,nb1), kdmx, psi_d(1,n_start), kdmx, &
                    ZERO, sc_d(nb1,n_start), nvecx )
        !$omp end target variant dispatch
        !
     END IF
     !
     if ( (n_start .le. n_end) .and. (mp_size(intra_bgrp_comm) > 1 ) ) then
        !pinned_buffer( nb1:nbase+notcnv, n_start:n_end ) = sc_d( nb1:nbase+notcnv, n_start:n_end )
        !ierr = cudaMemcpy2D( pinned_buffer(nb1, n_start) , nvecx, sc_d( nb1, n_start ), nvecx, notcnv, n_end-n_start+1 )
#if defined(__OPENMP_GPU)
        !$omp target teams distribute parallel do collapse(2) map(from:pinned_buffer)
        do j=n_start, n_end
           do i=nb1, nbase+notcnv
              pinned_buffer(i,j) = sc_d(i,j)
           enddo
        enddo
        !$omp end target teams distribute parallel do
#else
        CALL dev_memcpy( pinned_buffer, sc_d, (/ nb1, nbase + notcnv /), 1, (/ n_start, n_end /), 1 )
#endif
        CALL mp_sum( pinned_buffer( nb1:nbase+notcnv, n_start:n_end ), intra_bgrp_comm )
        !sc_d( nb1:nbase+notcnv, n_start:n_end ) = pinned_buffer( nb1:nbase+notcnv, n_start:n_end )
        !ierr = cudaMemcpy2D(  sc_d( nb1, n_start ), nvecx, pinned_buffer(nb1,n_start), nvecx, notcnv, n_end-n_start+1 )
#if defined(__OPENMP_GPU)
        !$omp target teams distribute parallel do collapse(2) map(to:pinned_buffer)
        do j=n_start, n_end
           do i=nb1, nbase+notcnv
              sc_d(i,j) = pinned_buffer(i,j)
           enddo
        enddo
        !$omp end target teams distribute parallel do
#else
        CALL dev_memcpy( sc_d, pinned_buffer, (/ nb1, nbase + notcnv /), 1, (/ n_start, n_end /), 1 )
#endif
     end if
     !$omp dispatch
     CALL mp_gather( sc_d, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
     !
     CALL mp_type_free( column_section_type )
     !
     CALL stop_clock( 'cegterg:overlap' )
     !
     nbase = nbase + notcnv
     !
!$cuf kernel do(1) <<<*,*>>>
!$omp target teams distribute parallel do
     DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real
        !
        IF( n>=nb1 ) THEN
           hc_d(n,n) = CMPLX( REAL( hc_d(n,n) ), 0.D0 ,kind=DP)
           sc_d(n,n) = CMPLX( REAL( sc_d(n,n) ), 0.D0 ,kind=DP)
        ENDIF
        !
        DO m = MAX(n+1,nb1), nbase
           !
           hc_d(n,m) = CONJG( hc_d(m,n) )
           sc_d(n,m) = CONJG( sc_d(m,n) )
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
#if defined(__OPENMP_GPU)
        CALL diaghg( nbase, nvec, hc_d, sc_d, nvecx, ew, vc_d, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
#else
        CALL diaghg( nbase, nvec, hc_d, sc_d, nvecx, ew_d, vc_d, me_bgrp, root_bgrp, intra_bgrp_comm )
#endif
     END IF
     IF( nbgrp > 1 ) THEN
        !$omp dispatch is_device_ptr(vc_d)
        CALL mp_bcast( vc_d, root_bgrp_id, inter_bgrp_comm )
#if defined(__OPENMP_GPU)
        CALL mp_bcast_mapped( ew, root_bgrp_id, inter_bgrp_comm )
#else
        CALL mp_bcast( ew_d, root_bgrp_id, inter_bgrp_comm )
#endif
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !
     ! ... test for convergence (on the CPU)
     !
#if defined(__OPENMP_GPU)
     !$omp target update from(ew(1:nvec))
     !$omp target update from(e_d(1:nvec))
#else
     ew_host(1:nvec) = ew_d(1:nvec)
     e_host(1:nvec) = e_d(1:nvec)
#endif
     WHERE( btype(1:nvec) == 1 )
#if defined(__OPENMP_GPU)
        conv(1:nvec) = ( ( ABS( ew     (1:nvec) - e_d   (1:nvec) ) < ethr ) )
#else
        conv(1:nvec) = ( ( ABS( ew_host(1:nvec) - e_host(1:nvec) ) < ethr ) )
#endif
     ELSEWHERE
#if defined(__OPENMP_GPU)
        conv(1:nvec) = ( ( ABS( ew     (1:nvec) - e_d   (1:nvec) ) < empty_ethr ) )
#else
        conv(1:nvec) = ( ( ABS( ew_host(1:nvec) - e_host(1:nvec) ) < empty_ethr ) )
#endif
     END WHERE
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
#if defined(__OPENMP_GPU)
     !$omp target teams loop
     do i=1, nvec
        e_d(i) = ew(i)
     enddo
     !$omp end target teams loop
#else
     CALL dev_memcpy (e_d, ew_d, (/ 1, nvec /) )
#endif
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. &
          nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'cegterg:last' )
        !
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
        !$omp target variant dispatch use_device_ptr(psi_d, vc_d, evc_d)
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, psi_d(1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                    ZERO, evc_d, kdmx )
        !$omp end target variant dispatch
#if defined(__OPENMP_GPU)
        CALL mp_sum_mapped( evc_d, inter_bgrp_comm )
#else
        CALL mp_sum( evc_d, inter_bgrp_comm )
#endif
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
#if defined(__OPENMP_GPU)
        !$omp target teams loop collapse(2)
        do j=1,nvec
           do i=1, npwx*npol
              psi_d(i,j) = evc_d(i,j)
           enddo
        enddo
        !$omp end target teams loop
#else
        CALL dev_memcpy(psi_d, evc_d, (/ 1, npwx*npol /), 1, &
                                      (/ 1, nvec /), 1)
#endif
        !
        IF ( uspp ) THEN
           !
           !$omp target variant dispatch use_device_ptr(spsi_d, vc_d, psi_d)
           CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, spsi_d(1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                       ZERO, psi_d(1,nvec+1), kdmx)
           !$omp end target variant dispatch
#if defined(__OPENMP_GPU)
           !$omp target teams loop collapse(2)
           do j=1,nvec
              do i=1, npwx*npol
                 spsi_d(i,j) = psi_d(i,j+nvec)
              enddo
           enddo
           !$omp end target teams loop
#else
           CALL dev_memcpy(spsi_d, psi_d(:,nvec+1:), &
                                        (/1, npwx*npol/), 1, &
                                        (/1, nvec/), 1)
#endif
           !$omp dispatch is_device_ptr(spsi_d)
           CALL mp_sum( spsi_d(:,1:nvec), inter_bgrp_comm )
           !
        END IF
        !
        !$omp target variant dispatch use_device_ptr(hpsi_d, vc_d, psi_d)
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, hpsi_d(1,n_start), kdmx, vc_d(n_start,1), nvecx, &
                    ZERO, psi_d(1,nvec+1), kdmx )
        !$omp end target variant dispatch
#if defined(__OPENMP_GPU)
        !$omp target teams loop collapse(2)
        do j=1,nvec
           do i=1, npwx*npol
              hpsi_d(i,j) = psi_d(i,j+nvec)
           enddo
        enddo
        !$omp end target teams loop
#else
        CALL dev_memcpy(hpsi_d, psi_d(:,nvec+1:), &
                                        (/1, npwx*npol/), 1, &
                                        (/1, nvec/), 1)
#endif
        !$omp dispatch is_device_ptr(hpsi_d)
        CALL mp_sum( hpsi_d(:,1:nvec), inter_bgrp_comm )
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        ! These variables are set to ZERO in the CUF Kernel below
        !hc_d(1:nbase,1:nbase) = ZERO
        !sc_d(1:nbase,1:nbase) = ZERO
        !vc_d(1:nbase,1:nbase) = ZERO
        !
        !$cuf kernel do(2) <<<*,*>>>
        !$omp target teams distribute parallel do collapse(2)
        DO n = 1, nbase
           DO j = 1, nbase
              !
              IF ( j == n ) THEN
                 hc_d(j,n) = CMPLX( e_d(n), 0.0_DP ,kind=DP)
                 !
                 sc_d(j,n) = ONE
                 vc_d(j,n) = ONE
              ELSE
                 hc_d(j,n) = ZERO; sc_d(j,n) = ZERO; vc_d(j,n) = ZERO
              END IF
           END DO
           !
        END DO
        !
        CALL stop_clock( 'cegterg:last' )
        !
     END IF
     !
  END DO iterate
  !
#if defined(__OPENMP_GPU)
  deallocate(pinned_buffer)
#else
  CALL gbuf%release_buffer(pinned_buffer, ierr)
#endif
  DEALLOCATE( recv_counts )
  DEALLOCATE( displs )
  DEALLOCATE( conv )
#if defined(__OPENMP_GPU)
  !$omp target exit data map(delete:ew)
  deallocate(ew)
#else
  DEALLOCATE( e_host, ew_host, ew_d )
#endif
  !
  IF ( uspp ) DEALLOCATE( spsi_d )
  !
  DEALLOCATE( psi_d )
  DEALLOCATE( hpsi_d )
  DEALLOCATE( vc_d )
  DEALLOCATE( hc_d )
  DEALLOCATE( sc_d )
  !
  !
  CALL stop_clock( 'cegterg' ); !write(*,*) 'stop cegterg' ; FLUSH(6)
  !call print_clock( 'cegterg' )
  !call print_clock( 'cegterg:init' )
  !call print_clock( 'cegterg:diag' )
  !call print_clock( 'cegterg:update' )
  !call print_clock( 'cegterg:overlap' )
  !call print_clock( 'cegterg:last' )
  !
  RETURN
  !
END SUBROUTINE cegterg_gpu

SUBROUTINE reorder_evals_cevecs(nbase, nvec, nvecx, conv, e_d, ew_d, v_d)
   USE util_param,    ONLY : DP
   USE device_fbuff_m,  ONLY : buffer => dev_buf
   implicit none
   INTEGER, INTENT(IN) :: nbase, nvec, nvecx
   LOGICAL, INTENT(IN) :: conv(nvec)
   REAL(DP)            :: e_d(nvecx), ew_d(nvecx)
   COMPLEX(DP)         :: v_d(nvecx,nvecx)
#if defined(__CUDA)
   attributes(DEVICE)  :: e_d, ew_d, v_d
#endif
   !
   INTEGER :: j, k, n, np, info
   INTEGER, ALLOCATABLE :: conv_idx(:)
   INTEGER, POINTER     :: conv_idx_d(:)
   COMPLEX(DP), POINTER :: vtmp_d(:,:)
#if defined(__CUDA)
   attributes(DEVICE)   :: conv_idx_d, vtmp_d
#endif
   !
   np = 0
   ALLOCATE(conv_idx(nvec))
   DO n = 1, nvec
      conv_idx(n) = -1
      IF ( .NOT. conv(n) ) THEN
         np = np + 1
         conv_idx(n) = np
      END IF
   END DO

   CALL buffer%lock_buffer(conv_idx_d, nvec, info)
   CALL buffer%lock_buffer(vtmp_d, (/nvecx, nvecx/), info)
   IF( info /= 0 ) &
     CALL errore( ' reorder_evals_cevecs ',' cannot allocate vtmp_d ', ABS(info) )

#if defined(__OPENMP_GPU)
   !$omp target enter data map(to:conv_idx)
#else
   CALL buffer%lock_buffer(conv_idx_d, nvec, info)
   conv_idx_d(1:nvec) = conv_idx(1:nvec)
#endif

!$cuf kernel do(2) <<<*,*>>>
!$omp target teams distribute parallel do collapse(2) is_device_ptr(vtmp_d,v_d)
   DO j=1,nvec
      DO k=1,nvecx
         vtmp_d(k,j) = v_d(k,j)
      END DO
   END DO
!$omp end target teams distribute parallel do

!$cuf kernel do(2) <<<*,*>>>
!$omp target teams distribute parallel do collapse(2) is_device_ptr(vtmp_d,v_d)
   DO j=1,nvec
      DO k=1,nvecx
#if defined(__OPENMP_GPU)
         IF(conv_idx(j) /= -1) THEN
           v_d(k,conv_idx(j)) = vtmp_d(k,j)
           IF(k==1) ew_d(nbase+conv_idx(j)) = e_d(j)
         END IF
#else
         IF(conv_idx_d(j) /= -1) THEN
           v_d(k,conv_idx_d(j)) = vtmp_d(k,j)
           IF(k==1) ew_d(nbase+conv_idx_d(j)) = e_d(j)
         END IF
#endif
      END DO
   END DO
   !
#if defined(__OPENMP_GPU)
   !$omp target exit data map(delete:conv_idx)
#else
   CALL buffer%release_buffer(conv_idx_d, info)
#endif
   CALL buffer%release_buffer(vtmp_d, info)
   !
   DEALLOCATE(conv_idx)
END SUBROUTINE reorder_evals_cevecs

!
!  Wrapper for subroutine with distributed matrixes (written by Carlo Cavazzoni)
!
!----------------------------------------------------------------------------
SUBROUTINE pcegterg_gpu(h_psi_gpu, s_psi_gpu, uspp, g_psi_gpu, &
                    npw, npwx, nvec, nvecx, npol, evc_d, ethr, &
                    e_d, btype, notcnv, lrot, dav_iter , nhpsi )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  !
  USE util_param,       ONLY : DP, stdout
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier, &
                               mp_size, mp_type_free, mp_allgather
#if defined(__OPENMP_GPU)
  USE omp_lib
#endif
  USE device_memcpy_m,  ONLY : dev_memcpy, dev_memset
  USE device_fbuff_m,   ONLY : buffer => dev_buf

  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! number of spin polarizations
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
    ! chunking parameters
  COMPLEX(DP), INTENT(INOUT) :: evc_d(npwx*npol,nvec)
#if defined(__CUDA)
   attributes(DEVICE)   :: evc_d
#endif
    !  evc   contains the  refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : S|psi> not needed
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e_d(nvec)
#if defined(__CUDA)
   attributes(DEVICE)   :: e_d
#endif
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  INTEGER, INTENT(OUT) :: nhpsi
    ! total number of indivitual hpsi
  !
  ! ... LOCAL variables
  !
#if !defined(__OPENMP_GPU)
  COMPLEX(DP), ALLOCATABLE :: evc(:,:)
  REAL(DP), ALLOCATABLE :: e(:)
  REAL(DP), POINTER :: ew_d(:)
  COMPLEX(DP), POINTER :: psi_d(:,:), hpsi_d(:,:), spsi_d(:,:)
#endif

  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, ipol, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
  INTEGER :: i, j, k, ierr
  REAL(DP), ALLOCATABLE :: ew(:)
#if defined(__CUDA)
  attributes(DEVICE) :: ew_d
#endif
  COMPLEX(DP), ALLOCATABLE :: hl(:,:), sl(:,:), vl(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
#if defined(__CUDA)
  attributes(DEVICE) ::  psi_d, hpsi_d, spsi_d
#endif
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr
    ! threshold for empty bands
  INTEGER :: idesc(LAX_DESC_SIZE), idesc_old(LAX_DESC_SIZE)
  INTEGER, ALLOCATABLE :: irc_ip( : )
  INTEGER, ALLOCATABLE :: nrc_ip( : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: notcnv_ip( : )
  INTEGER, ALLOCATABLE :: ic_notcnv( : )
  !
  INTEGER :: np_ortho(2), ortho_parent_comm
  LOGICAL :: do_distr_diag_inside_bgrp
  !
  REAL(DP), EXTERNAL :: ddot
  !
  EXTERNAL  h_psi_gpu, s_psi_gpu, g_psi_gpu
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  nhpsi = 0
  CALL start_clock( 'cegterg' )
  !
  CALL laxlib_getval( np_ortho = np_ortho, ortho_parent_comm = ortho_parent_comm, &
    do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'pcegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  ! compute the number of chuncks
  numblock  = (npw+blocksize-1)/blocksize

#if !defined(__OPENMP_GPU)
  ALLOCATE(  evc( npwx*npol, nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate evc (host) ', ABS(ierr) )
  !
  ALLOCATE(  e( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate e (host) ', ABS(ierr) )
#endif
  !
  ALLOCATE(  psi( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate psi ', ABS(ierr) )
  !
  ALLOCATE( hpsi( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx*npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ! ... Initialize the matrix descriptor
  !
  ALLOCATE( ic_notcnv( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ic_notcnv ', ABS(ierr) )
  !
  ALLOCATE( notcnv_ip( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate notcnv_ip ', ABS(ierr) )
  !
  ALLOCATE( irc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate irc_ip ', ABS(ierr) )
  !
  ALLOCATE( nrc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate nrc_ip ', ABS(ierr) )
  !
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate rank_ip ', ABS(ierr) )
  !
  CALL desc_init( nvec, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
  !
  IF( la_proc ) THEN
     !
     ! only procs involved in the diagonalization need to allocate local
     ! matrix block.
     !
     ALLOCATE( vl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
     !
  ELSE
     !
     ALLOCATE( vl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ew ', ABS(ierr) )
  !
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  IF ( uspp ) spsi = ZERO
  !
  hpsi = ZERO
  psi  = ZERO
#if defined(__OPENMP_GPU)
  !$omp target enter data map(alloc:psi, hpsi, spsi, ew)
#else
  CALL buffer%lock_buffer(psi_d, (/npwx*npol, nvecx/), ierr)
  CALL buffer%lock_buffer(hpsi_d, (/npwx*npol, nvecx/), ierr)
  CALL buffer%lock_buffer(spsi_d, (/npwx*npol, nvecx/), ierr)
  CALL buffer%lock_buffer(ew_d, nvecx, ierr)
#endif

#if defined(__OPENMP_GPU)
  !$omp target update from(evc_d)
  psi(:,1:nvec) = evc_d(:,1:nvec)  !With OpenMP, evc_d is present on host (already mapped)
  !$omp target teams loop collapse(2)
  do j=1, nvec
     do i=1, npwx*npol
        psi(i,j) = evc_d(i,j)
     enddo
  enddo
  !$omp end target teams loop
#else
  evc(:,1:nvec) = evc_d(:,1:nvec)
  psi(:,1:nvec) = evc(:,1:nvec)
  CALL dev_memcpy(psi_d, evc_d, (/1, npwx*npol /), 1 , (/ 1, nvec /) )
#endif
  !
  ! ... hpsi contains h times the basis vectors
  !
#if defined(__OPENMP_GPU)
  associate(psi_d => psi, spsi_d => spsi, evc => evc_d, e => e_d, ew_d => ew)
#endif
  CALL h_psi_gpu( npwx, npw, nvec, psi_d, hpsi ) ; nhpsi = nhpsi + nvec
#if defined(__OPENMP_GPU)
  !$omp target update from(hpsi)
#else
  hpsi(1:npwx*npol, 1:nvec) = hpsi_d(1:npwx*npol, 1:nvec)
#endif
  !
  IF ( uspp ) CALL s_psi_gpu( npwx, npw, nvec, psi_d, spsi_d )
  IF ( uspp ) THEN
#if defined(__OPENMP_GPU)
     !$omp target update from(spsi)
#else
     spsi(1:npwx*npol, 1:nvec) = spsi_d(1:npwx*npol, 1:nvec)
#endif
  ENDIF
  !
  ! ... hl contains the projection of the hamiltonian onto the reduced
  ! ... space, vl contains the eigenvectors of hl. Remember hl, vl and sl
  ! ... are all distributed across processors, global replicated matrixes
  ! ... here are never allocated
  !
  CALL start_clock( 'cegterg:init' )

  CALL compute_distmat( hl, psi, hpsi )
  !
  IF ( uspp ) THEN
     !
     CALL compute_distmat( sl, psi, spsi )
     !
  ELSE
     !
     CALL compute_distmat( sl, psi, psi )
     !
  END IF
  CALL stop_clock( 'cegterg:init' )
  !
  IF ( lrot ) THEN
     !
     CALL set_e_from_h(e, hl)
#if defined(__OPENMP_GPU)
     !$omp target update to(e_d)
#else
     e_d = e
#endif
     !
     CALL set_to_identity( vl, idesc )
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !     Calling block parallel algorithm
     !
     CALL start_clock( 'cegterg:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
     END IF
     CALL stop_clock( 'cegterg:diag' )
     !
     e(1:nvec) = ew(1:nvec)
#if defined(__OPENMP_GPU)
     !$omp target update to(e)
#else
     e_d(1:nvec) = ew(1:nvec)
#endif
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !
     CALL start_clock( 'cegterg:update' )
     !
     CALL reorder_v(ew, e)
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL hpsi_dot_v()
     !
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration
     !
#if defined(__OPENMP_GPU)
     !$omp target update to(ew, psi(1:npwx*npol, nb1:nb1+notcnv))
#else
     ew_d = ew
     psi_d(1:npwx*npol, nb1:nb1+notcnv) = psi(1:npwx*npol, nb1:nb1+notcnv)
#endif
     CALL g_psi_gpu( npwx, npw, notcnv, npol, psi_d(1,nb1), ew_d(nb1) )
#if defined(__OPENMP_GPU)
     !$omp target update from(psi(1:npwx*npol, nb1:nb1+notcnv))
#else
     psi(1:npwx*npol, nb1:nb1+notcnv) = psi_d(1:npwx*npol, nb1:nb1+notcnv)
#endif
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
           ew(n) = ddot( 2*npw, psi(1,nbn), 1, psi(1,nbn), 1 )
           !
        ELSE
           !
           ew(n) = ddot( 2*npw, psi(1,nbn), 1, psi(1,nbn), 1 ) + &
                   ddot( 2*npw, psi(npwx+1,nbn), 1, psi(npwx+1,nbn), 1 )
           !
        END IF
        !
     END DO
     !
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !
     !$omp parallel do collapse(3)
     DO n = 1, notcnv
        DO ipol = 1, npol
           DO m = 1, numblock
              psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw, m*blocksize)+(ipol-1)*npwx,nbase+n) = &
              psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw, m*blocksize)+(ipol-1)*npwx,nbase+n) / &
                    SQRT( ew(n) )
           END DO
        END DO
     END DO
     !$omp end parallel do
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
#if defined(__OPENMP_GPU)
     !$omp target update to(psi(1:npwx*npol, nb1:nb1+notcnv))
#else
     psi_d(1:npwx*npol, nb1:nb1+notcnv) = psi(1:npwx*npol, nb1:nb1+notcnv)
#endif
#if defined(__OPENMP_GPU)
     CALL h_psi_gpu( npwx, npw, notcnv, psi_d(1,nb1), hpsi(1,nb1) ) ; nhpsi = nhpsi + notcnv
     !$omp target update from(hpsi(1:npwx*npol, nb1:nb1+notcnv))
#else
     CALL h_psi_gpu( npwx, npw, notcnv, psi_d(1,nb1), hpsi_d(1,nb1) ) ; nhpsi = nhpsi + notcnv
     hpsi(1:npwx*npol, nb1:nb1+notcnv) = hpsi_d(1:npwx*npol, nb1:nb1+notcnv)
#endif
     !
     IF ( uspp ) CALL s_psi_gpu( npwx, npw, notcnv, psi_d(1,nb1), spsi_d(1,nb1) )
     IF ( uspp ) THEN
#if defined(__OPENMP_GPU)
        !$omp target update from(spsi)
#else
        spsi = spsi_d
#endif
     ENDIF
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:overlap' )
     !
     ! we need to save the old descriptor in order to redistribute matrices
     !
     idesc_old = idesc
     !
     ! ... RE-Initialize the matrix descriptor
     !
     CALL desc_init( nbase+notcnv, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
     !
     IF( la_proc ) THEN

        !  redistribute hl and sl (see dsqmred), since the dimension of the subspace has changed
        !
        vl = hl
        DEALLOCATE( hl )
        ALLOCATE( hl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )

        CALL laxlib_zsqmred( nbase, vl, idesc_old(LAX_DESC_NRCX), idesc_old, nbase+notcnv, hl, nx, idesc )

        vl = sl
        DEALLOCATE( sl )
        ALLOCATE( sl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )

        CALL laxlib_zsqmred( nbase, vl, idesc_old(LAX_DESC_NRCX), idesc_old, nbase+notcnv, sl, nx, idesc )

        DEALLOCATE( vl )
        ALLOCATE( vl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )

     END IF
     !
     !
     CALL update_distmat( hl, psi, hpsi )
     !
     IF ( uspp ) THEN
        !
        CALL update_distmat( sl, psi, spsi )
        !
     ELSE
        !
        CALL update_distmat( sl, psi, psi )
        !
     END IF
     !
     CALL stop_clock( 'cegterg:overlap' )
     !
     nbase = nbase + notcnv
     !
     ! ... diagonalize the reduced hamiltonian
     !     Call block parallel algorithm
     !
     CALL start_clock( 'cegterg:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
     END IF
     CALL stop_clock( 'cegterg:diag' )
     !
     ! ... test for convergence
     !
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
#if defined(__OPENMP_GPU)
     !$omp target update to(e_d)
#else
     e_d(1:nvec) = e(1:nvec)
#endif
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'cegterg:last' )
        !
        CALL refresh_evc()
#if defined(__OPENMP_GPU)
        !$omp target update to(evc_d)
#else
        evc_d = evc
#endif
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        CALL threaded_memcpy(psi, evc, nvec*npol*npwx*2)
        !
        IF ( uspp ) THEN
           !
           CALL refresh_spsi()
           !
        END IF
        !
        CALL refresh_hpsi()
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        CALL desc_init( nvec, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
        !
        IF( la_proc ) THEN
           !
           ! note that nx has been changed by desc_init
           ! we need to re-alloc with the new size.
           !
           DEALLOCATE( vl, hl, sl )
           ALLOCATE( vl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
           ALLOCATE( hl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
           ALLOCATE( sl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
           !
        END IF
        !
        CALL set_h_from_e(hl, e)
        !
        CALL set_to_identity( vl, idesc )
        CALL set_to_identity( sl, idesc )
        !
        CALL stop_clock( 'cegterg:last' )
        !
     END IF
     !
  END DO iterate
  !
#if defined(__OPENMP_GPU)
  endassociate
  !$omp target exit data map(delete:psi, hpsi, spsi, ew)
#else
  CALL buffer%release_buffer(psi_d, ierr)
  CALL buffer%release_buffer(hpsi_d, ierr)
  CALL buffer%release_buffer(spsi_d, ierr)
  CALL buffer%release_buffer(ew_d, ierr)
  DEALLOCATE( evc )
  DEALLOCATE( e )
#endif
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
  !
  DEALLOCATE( vl, hl, sl )
  !
  DEALLOCATE( rank_ip )
  DEALLOCATE( ic_notcnv )
  DEALLOCATE( irc_ip )
  DEALLOCATE( nrc_ip )
  DEALLOCATE( notcnv_ip )
  DEALLOCATE( conv )
  DEALLOCATE( ew )

  CALL stop_clock( 'cegterg' )
  !call print_clock( 'cegterg' )
  !call print_clock( 'cegterg:init' )
  !call print_clock( 'cegterg:diag' )
  !call print_clock( 'cegterg:update' )
  !call print_clock( 'cegterg:overlap' )
  !call print_clock( 'cegterg:last' )
  !
  RETURN
  !
  !
CONTAINS
  !
  SUBROUTINE set_to_identity( distmat, idesc )
     INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
     COMPLEX(DP), INTENT(OUT) :: distmat(:,:)
     INTEGER :: i
     distmat = ( 0_DP , 0_DP )
     IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) .AND. idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
        DO i = 1, idesc(LAX_DESC_NC)
           distmat( i, i ) = ( 1_DP , 0_DP )
        END DO
     END IF
     RETURN
  END SUBROUTINE set_to_identity
  !
  !
  SUBROUTINE reorder_v(ew, e)
     !
     REAL(DP), INTENT(OUT) :: ew(:)
     REAL(DP), INTENT(IN)  :: e(:)
     INTEGER :: ipc
     INTEGER :: nc, ic
     INTEGER :: nl, npl
     !
     np = 0
     !
     notcnv_ip = 0
     !
     n = 0
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        npl = 0
        !
        IF( ic <= nvec ) THEN
           !
           DO nl = 1, min( nvec - ic + 1, nc )
              !
              n  = n  + 1
              !
              IF ( .NOT. conv(n) ) THEN
                 !
                 ! ... this root not yet converged ...
                 !
                 np  = np  + 1
                 npl = npl + 1
                 IF( npl == 1 ) ic_notcnv( ipc ) = np
                 !
                 ! ... reorder eigenvectors so that coefficients for unconverged
                 ! ... roots come first. This allows to use quick matrix-matrix
                 ! ... multiplications to set a new basis vector (see below)
                 !
                 notcnv_ip( ipc ) = notcnv_ip( ipc ) + 1
                 !
                 IF ( npl /= nl ) THEN
                    IF( la_proc .AND. idesc(LAX_DESC_MYC) == ipc-1 ) THEN
                       vl( :, npl) = vl( :, nl )
                    END IF
                 END IF
                 !
                 ! ... for use in g_psi
                 !
                 ew(nbase+np) = e(n)
                 !
              END IF
              !
           END DO
           !
        END IF
        !
     END DO
     !
  END SUBROUTINE reorder_v
  !
  !
  SUBROUTINE hpsi_dot_v()
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, ir, ic, notcl, root, np, ipol, ib
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP), ALLOCATABLE :: ptmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     ALLOCATE( ptmp( npwx*npol, nx ) )

     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        IF( notcnv_ip( ipc ) > 0 ) THEN

           notcl = notcnv_ip( ipc )
           ic    = ic_notcnv( ipc )

           beta = ZERO

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 vtmp(:,1:notcl) = vl(:,1:notcl)
              END IF

              CALL mp_bcast( vtmp(:,1:notcl), root, ortho_parent_comm )
              !
              IF ( uspp ) THEN
                 !
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    spsi(1, ir), kdmx, vtmp, nx, beta, psi(1,nb1+ic-1), kdmx )
                 !
              ELSE
                 !
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    psi(1, ir), kdmx, vtmp, nx, beta, psi(1,nb1+ic-1), kdmx )
                 !
              END IF
              !
              CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                      hpsi(1, ir), kdmx, vtmp, nx, beta, ptmp, kdmx )

              beta = ONE

           END DO

           !$omp parallel do collapse(3)
           DO np = 1, notcl
              DO ipol = 1, npol
                 DO ib = 1, numblock
                    !
                    psi( (ib-1)*blocksize+(ipol-1)*npwx+1: &
                         MIN(npw, ib*blocksize)+(ipol-1)*npwx,nbase+np+ic-1) = &
                    ptmp((ib-1)*blocksize+(ipol-1)*npwx+1: &
                         MIN(npw, ib*blocksize)+(ipol-1)*npwx,np) - &
                    ew(nbase+np+ic-1) * psi((ib-1)*blocksize+(ipol-1)*npwx+1:&
                       MIN(npw, ib*blocksize)+(ipol-1)*npwx,nbase+np+ic-1)
                    !
                 END DO
              END DO
           END DO
           !$omp end parallel do
           !
           ! clean up garbage if there is any
           IF (npw < npwx) psi(npw+1:npwx,nbase+ic:nbase+notcl+ic-1) = ZERO
           IF (npol == 2)  psi(npwx+npw+1:2*npwx,nbase+ic:nbase+notcl+ic-1) = ZERO
           !
        END IF
        !
     END DO

     DEALLOCATE( vtmp )
     DEALLOCATE( ptmp )

     RETURN
  END SUBROUTINE hpsi_dot_v
  !
  !
  SUBROUTINE refresh_evc( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
#if defined(__OPENMP_GPU)
                          psi(1,ir), kdmx, vl, nx, beta, evc_d(1,ic), kdmx )
#else
                          psi(1,ir), kdmx, vl, nx, beta, evc(1,ic), kdmx )
#endif
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
#if defined(__OPENMP_GPU)
                          psi(1,ir), kdmx, vtmp, nx, beta, evc_d(1,ic), kdmx )
#else
                          psi(1,ir), kdmx, vtmp, nx, beta, evc(1,ic), kdmx )
#endif
              END IF
              !

              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_evc
  !
  !
  SUBROUTINE refresh_spsi( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO
           !
           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi(1,ir), kdmx, vl, nx, beta, psi(1,nvec+ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi(1,ir), kdmx, vtmp, nx, beta, psi(1,nvec+ic), kdmx )
              END IF
              !
              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     CALL threaded_memcpy(spsi, psi(1,nvec+1), nvec*npol*npwx*2)
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_spsi
  !
  !
  !
  SUBROUTINE refresh_hpsi( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO
           !
           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi(1,ir), kdmx, vl, nx, beta, psi(1,nvec+ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi(1,ir), kdmx, vtmp, nx, beta, psi(1,nvec+ic), kdmx )
              END IF
              !
              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )
     !
     CALL threaded_memcpy(hpsi, psi(1,nvec+1), nvec*npol*npwx*2)
     !
     RETURN
  END SUBROUTINE refresh_hpsi
  !
  !
  SUBROUTINE compute_distmat( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = ZERO
     !
     !  Only upper triangle is computed, then the matrix is hermitianized
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        DO ipr = 1, ipc ! idesc(LAX_DESC_NPR) ! ipc ! use symmetry for the loop on row procs
           !
           nr = nrc_ip( ipr )
           ir = irc_ip( ipr )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE , &
                       v(1,ir), kdmx, w(1,ic), kdmx, ZERO, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     !  The matrix is hermitianized using upper triangle
     !
     CALL laxlib_zsqmher( nbase, dm, nx, idesc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  !
  SUBROUTINE update_distmat( dm, v, w )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, icc, ii
     COMPLEX(DP) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )

     ALLOCATE( vtmp( nx, nx ) )
     !
     vtmp = ZERO
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic+nc-1 >= nb1 ) THEN
           !
           nc = MIN( nc, ic+nc-1 - nb1 + 1 )
           IF( ic >= nb1 ) THEN
              ii = ic
              icc = 1
           ELSE
              ii = nb1
              icc = nb1-ic+1
           END IF
           !
           ! icc to nc is the local index of the unconverged bands
           ! ii is the global index of the first unconverged bands
           !
           DO ipr = 1, ipc ! idesc(LAX_DESC_NPR) use symmetry
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE, v(1, ir), &
                          kdmx, w(1,ii), kdmx, ZERO, vtmp, nx )
              IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) vtmp = vtmp/nbgrp
              !
              IF(  (idesc(LAX_DESC_ACTIVE_NODE) > 0) .AND. &
                   (ipr-1 == idesc(LAX_DESC_MYR)) .AND. (ipc-1 == idesc(LAX_DESC_MYC)) ) THEN
                 CALL mp_root_sum( vtmp(:,1:nc), dm(:,icc:icc+nc-1), root, ortho_parent_comm )
              ELSE
                 CALL mp_root_sum( vtmp(:,1:nc), dm, root, ortho_parent_comm )
              END IF

           END DO
           !
        END IF
        !
     END DO
     !
     CALL laxlib_zsqmher( nbase+notcnv, dm, nx, idesc )
     !
     DEALLOCATE( vtmp )
     RETURN
  END SUBROUTINE update_distmat
  !
  !
  SUBROUTINE set_e_from_h(e, h)
     IMPLICIT NONE
     REAL(DP), INTENT(OUT) :: e(:)
     COMPLEX(DP), INTENT(IN) :: h(:,:)
     INTEGER :: nc, ic, i
     e(1:nbase) = 0_DP
     IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) .AND. la_proc ) THEN
        nc = idesc(LAX_DESC_NC)
        ic = idesc(LAX_DESC_IC)
        DO i = 1, nc
           e( i + ic - 1 ) = REAL( hl( i, i ) )
        END DO
     END IF
     CALL mp_sum( e(1:nbase), ortho_parent_comm )
     RETURN
  END SUBROUTINE set_e_from_h
  !
  SUBROUTINE set_h_from_e(h, e)
     IMPLICIT NONE
     COMPLEX(DP), INTENT(OUT) :: h(:,:)
     REAL(DP), INTENT(IN) :: e(:)
     INTEGER :: nc, ic, i
     IF( la_proc ) THEN
        h = ZERO
        IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           nc = idesc(LAX_DESC_NC)
           ic = idesc(LAX_DESC_IC)
           DO i = 1, nc
              h(i,i) = CMPLX( e( i + ic - 1 ), 0_DP ,kind=DP)
           END DO
        END IF
     END IF
     RETURN
  END SUBROUTINE set_h_from_e
  !
END SUBROUTINE pcegterg_gpu
