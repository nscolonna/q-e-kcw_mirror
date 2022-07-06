!
! Copyright (C) 2003-2013 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_gamma( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - Gamma point.
  !
  USE parallel_include
  USE kinds,                   ONLY : DP
  USE mp_bands,                ONLY : me_bgrp
  USE fft_base,                ONLY : dffts, dfftp
  USE fft_interfaces,          ONLY : fwfft, invfft
  USE wavefunctions,           ONLY : psic, psic_omp
  USE fft_helper_subroutines,  ONLY : fftx_ntgrp, tg_get_nnr, &
                                      tg_get_group_nr3, tg_get_recip_inc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,m)
  !! Hamiltonian dot psi
  REAL(DP), INTENT(IN) :: v(dffts%nnr)
  !! the total pot. in real space (smooth grid) for current spin
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, incr
  INTEGER :: right_nnr, right_nr3, right_inc
  COMPLEX(DP) :: fp, fm
  INTEGER, ALLOCATABLE :: dffts_nl(:), dffts_nlm(:)
  !
  !Variables for task groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, v_siz_p, idx, ioff
  !
  CALL start_clock ('vloc_psi')
  incr = 2
  !
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     !
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz =  dffts%nnr_tg
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( v_siz ) )
     !
     CALL tg_gather( dffts, v, tg_v )
     CALL stop_clock ('vloc_psi:tg_gather')
     !
     incr = 2 * fftx_ntgrp(dffts)
     !
  ELSE
     !
     ALLOCATE(dffts_nl(1:dffts%ngm))
     ALLOCATE(dffts_nlm(1:dffts%ngm))
     !
     dffts_nl = dffts%nl
     dffts_nlm = dffts%nlm
     v_siz = dffts%nnr
     v_siz_p = dfftp%nnr
#if defined(__OPENMP_GPU)
     !$omp target enter data map(to: dffts_nl, dffts_nlm, v, psi, hpsi)
#endif
     !
  ENDIF
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  DO ibnd = 1, m, incr
     !
     IF( use_tg ) THEN
        !
        CALL tg_get_nnr( dffts, right_nnr )
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, 2*fftx_ntgrp(dffts), 2
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 tg_psic(dffts%nl (j)+ioff) =        psi(j,idx+ibnd-1) + &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd)
                 tg_psic(dffts%nlm(j)+ioff) = conjg( psi(j,idx+ibnd-1) - &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd) )
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 tg_psic(dffts%nl (j)+ioff) =        psi(j,idx+ibnd-1)
                 tg_psic(dffts%nlm(j)+ioff) = conjg( psi(j,idx+ibnd-1) )
              ENDDO
           ENDIF

           ioff = ioff + right_nnr

        ENDDO
        !
     ELSE
        !
        !$omp target teams distribute parallel do
        DO j = 1, v_siz_p
          psic_omp(j) = (0.d0, 0.d0)
        END DO
        IF (ibnd < m) THEN
           ! two ffts at the same time
           !$omp target teams distribute parallel do
           DO j = 1, n
              psic_omp(dffts_nl (j))=      psi(j,ibnd) + (0.0d0,1.d0)*psi(j,ibnd+1)
              psic_omp(dffts_nlm(j))=conjg(psi(j,ibnd) - (0.0d0,1.d0)*psi(j,ibnd+1))
           ENDDO
        ELSE
           !$omp target teams distribute parallel do
           DO j = 1, n
              psic_omp (dffts_nl (j)) =       psi(j, ibnd)
              psic_omp (dffts_nlm(j)) = conjg(psi(j, ibnd))
           ENDDO
        ENDIF
        !
     ENDIF
     !
     !   fft to real space
     !   product with the potential v on the smooth grid
     !   back to reciprocal space
     !
     IF( use_tg ) THEN
        !
        CALL invfft ('tgWave', tg_psic, dffts )
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
        DO j = 1, dffts%nr1x * dffts%nr2x * right_nr3
           tg_psic (j) = tg_psic (j) * tg_v(j)
        ENDDO
        !
        CALL fwfft ('tgWave', tg_psic, dffts )
        !
     ELSE
        !
#if defined(__USE_DISPATCH)
        !$omp dispatch
#endif
        CALL invfft ('Wave', psic_omp, dffts)
        !
!$omp target teams distribute parallel do
        DO j = 1, v_siz
           psic_omp (j) = psic_omp (j) * v(j)
        ENDDO
        !
#if defined(__USE_DISPATCH)
        !$omp dispatch
#endif
        CALL fwfft ('Wave', psic_omp, dffts)
        !
     ENDIF
     !
     !   addition to the total product
     !
     IF( use_tg ) THEN
        !
        ioff   = 0
        !
        CALL tg_get_recip_inc( dffts, right_inc )
        !
        DO idx = 1, 2*fftx_ntgrp(dffts), 2
           !
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 fp= ( tg_psic( dffts%nl(j) + ioff ) +  &
                       tg_psic( dffts%nlm(j) + ioff ) ) * 0.5d0
                 fm= ( tg_psic( dffts%nl(j) + ioff ) -  &
                       tg_psic( dffts%nlm(j) + ioff ) ) * 0.5d0
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                        cmplx( dble(fp), aimag(fm),kind=DP)
                 hpsi (j, ibnd+idx  ) = hpsi (j, ibnd+idx  ) + &
                                        cmplx(aimag(fp),- dble(fm),kind=DP)
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                         tg_psic( dffts%nl(j) + ioff )
              ENDDO
           ENDIF
           !
           ioff = ioff + right_inc
           !
        ENDDO
        !
     ELSE
        IF (ibnd < m) THEN
           ! two ffts at the same time
           !$omp target teams distribute parallel do
           DO j = 1, n
              fp = (psic_omp (dffts_nl(j)) + psic_omp (dffts_nlm(j)))*0.5d0
              fm = (psic_omp (dffts_nl(j)) - psic_omp (dffts_nlm(j)))*0.5d0
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + &
                                 cmplx( dble(fp), aimag(fm),kind=DP)
              hpsi (j, ibnd+1) = hpsi (j, ibnd+1) + &
                                 cmplx(aimag(fp),- dble(fm),kind=DP)
           ENDDO
        ELSE
           !$omp target teams distribute parallel do
           DO j = 1, n
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic_omp (dffts_nl(j))
           ENDDO
        ENDIF
     ENDIF
     !
  ENDDO
  !
  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
     !
   ELSE
     !
#if defined(__OPENMP_GPU)
     !$omp target exit data map(delete:dffts_nl,dffts_nlm,v,psi) map(from:hpsi)
#endif
     DEALLOCATE( dffts_nl )
     DEALLOCATE( dffts_nlm )
     !
  ENDIF
  CALL stop_clock ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - k-points:
  !
  !! * fft to real space;
  !! * product with the potential v on the smooth grid;
  !! * back to reciprocal space;
  !! * addition to the hpsi.
  !
  USE parallel_include
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE klist,                  ONLY : igk_k
  USE mp_bands,               ONLY : me_bgrp
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : fwfft, invfft
  USE fft_helper_subroutines, ONLY : fftx_ntgrp, tg_get_nnr, &
                                     tg_get_group_nr3, tg_get_recip_inc
  USE wavefunctions,          ONLY : psic, psic_omp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,m)
  !! Hamiltonian dot psi
  REAL(DP), INTENT(IN) :: v(dffts%nnr)
  !! the total pot. in real space (smooth grid) for current spin
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, incr
  INTEGER :: i, right_nnr, right_nr3, right_inc
  INTEGER, ALLOCATABLE :: dffts_nl(:)
  !
  ! chunking parameters
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
  !
  ! Task Groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx
  !
  CALL start_clock ('vloc_psi')
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     !
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz =  dffts%nnr_tg
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( v_siz ) )
     !
     CALL tg_gather( dffts, v, tg_v )
     CALL stop_clock ('vloc_psi:tg_gather')
     !
  ELSE
     !     
     ALLOCATE(dffts_nl(1:dffts%ngm))
     dffts_nl = dffts%nl
     v_siz=dffts%nnr
     !
#if defined(__OPENMP_GPU)
     !$omp target enter data map(to:v,igk_k,dffts_nl,psi,hpsi)
#endif
     !
  ENDIF
  !
  IF( use_tg ) THEN

     CALL tg_get_nnr( dffts, right_nnr )

     ! compute the number of chuncks
     numblock  = (n+blocksize-1)/blocksize

     DO ibnd = 1, m, fftx_ntgrp(dffts)
        !
!$omp parallel
        CALL threaded_barrier_memset(tg_psic, 0.D0, fftx_ntgrp(dffts)*right_nnr*2)
        !$omp do collapse(2)
        DO idx = 0, MIN(fftx_ntgrp(dffts)-1, m-ibnd)
           DO j = 1, numblock
              tg_psic(dffts%nl (igk_k((j-1)*blocksize+1:MIN(j*blocksize, n),current_k))+right_nnr*idx) = &
                 psi((j-1)*blocksize+1:MIN(j*blocksize, n),idx+ibnd)
           ENDDO
        ENDDO
        !$omp end do nowait
!$omp end parallel
        !
        CALL  invfft ('tgWave', tg_psic, dffts )
        !write (6,*) 'wfc R ' 
        !write (6,99) (tg_psic(i), i=1,400)
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
!$omp parallel do
        DO j = 1, dffts%nr1x*dffts%nr2x* right_nr3
           tg_psic (j) = tg_psic (j) * tg_v(j)
        ENDDO
!$omp end parallel do
        !write (6,*) 'v psi R ' 
        !write (6,99) (tg_psic(i), i=1,400)
        !
        CALL fwfft ('tgWave',  tg_psic, dffts )
        !
        !   addition to the total product
        !
        CALL tg_get_recip_inc( dffts, right_inc )
        !
!$omp parallel do collapse(2)
        DO idx = 0, MIN(fftx_ntgrp(dffts)-1, m-ibnd)
           DO j = 1, numblock
              hpsi ((j-1)*blocksize+1:MIN(j*blocksize, n), ibnd+idx) = &
                 hpsi ((j-1)*blocksize+1:MIN(j*blocksize, n), ibnd+idx) + &
                 tg_psic( dffts%nl(igk_k((j-1)*blocksize+1:MIN(j*blocksize, n),current_k)) + right_inc*idx )
           ENDDO
        ENDDO
!$omp end parallel do
        !
     ENDDO
  ELSE
     DO ibnd = 1, m
        !
#if defined(__OPENMP_GPU)
        !$omp target teams distribute parallel do
        do i=1,v_siz
           psic_omp(i)=(0.d0, 0.d0)
        enddo
        !$omp target teams distribute parallel do
        DO j = 1, n
          psic_omp (dffts_nl (igk_k(j, current_k))) = psi(j, ibnd)
        END DO
        !
#else
!$omp parallel
        CALL threaded_barrier_memset(psic_omp, 0.D0, dffts%nnr*2)
        !$omp do
        DO j = 1, n
           psic_omp (dffts_nl (igk_k(j,current_k))) = psi(j, ibnd)
        ENDDO
        !$omp end do nowait
!$omp end parallel

#endif

        !write (6,*) 'wfc G ', ibnd
        !write (6,99) (psic(i), i=1,400)
        !
#if defined(__USE_DISPATCH)
        !$omp dispatch
#endif
        CALL invfft ('Wave', psic_omp, dffts)
        !write (6,*) 'wfc R ' 
        !write (6,99) (psic(i), i=1,400)
        !
        !$omp target teams distribute parallel do 
        DO j = 1, v_siz
           psic_omp (j) = psic_omp (j) * v(j)
        ENDDO
        !write (6,*) 'v psi R ' 
        !write (6,99) (psic(i), i=1,400)
        !
#if defined(__USE_DISPATCH)
        !$omp dispatch
#endif
        CALL fwfft ('Wave', psic_omp, dffts)
!!$omp target update from(psic)
        !
        !   addition to the total product
        !
        !$omp target teams distribute parallel do 
        DO j = 1, n
           hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic_omp (dffts_nl(igk_k(j,current_k)))
        ENDDO
        !write (6,*) 'v psi G ', ibnd
        !write (6,99) (psic(i), i=1,400)
        !
     ENDDO
  ENDIF
  !
  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
     !
  ELSE
     !
#if defined(__OPENMP_GPU)
     !$omp target exit data map(delete:dffts_nl,v,igk_k,psi) map(from:hpsi)
#endif
     DEALLOCATE(dffts_nl)
     !
  ENDIF
  CALL stop_clock ('vloc_psi')
  !
99 format ( 20 ('(',2f12.9,')') )

  RETURN
END SUBROUTINE vloc_psi_k
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_nc( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - noncollinear.
  !
  USE parallel_include
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE klist,                  ONLY : igk_k
  USE mp_bands,               ONLY : me_bgrp
  USE fft_base,               ONLY : dffts, dfftp
  USE fft_interfaces,         ONLY : fwfft, invfft
  USE lsda_mod,               ONLY : nspin
  USE noncollin_module,       ONLY : npol, domag
  USE wavefunctions,          ONLY : psic_nc, psic_nc_omp
  USE fft_helper_subroutines, ONLY : fftx_ntgrp, tg_get_nnr, &
                                     tg_get_group_nr3, tg_get_recip_inc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  REAL(DP), INTENT(IN) :: v(dfftp%nnr,4) ! beware dimensions!
  !! the total pot. in real space (smooth grid)
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,npol,m)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j,ipol, incr, is
  COMPLEX(DP) :: sup, sdwn
  INTEGER, ALLOCATABLE :: dffts_nl(:)
  !
  ! Variables for task groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v(:,:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:,:)
  INTEGER :: v_siz, v_siz_p, idx, ioff
  INTEGER :: right_nnr, right_nr3, right_inc
  !
  CALL start_clock ('vloc_psi')
  !
  incr = 1
  !
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz = dffts%nnr_tg
     IF (domag) THEN
        ALLOCATE( tg_v( v_siz, 4 ) )
        DO is=1,nspin
           CALL tg_gather( dffts, v(:,is), tg_v(:,is) )
        ENDDO
     ELSE
        ALLOCATE( tg_v( v_siz, 1 ) )
        CALL tg_gather( dffts, v(:,1), tg_v(:,1) )
     ENDIF
     ALLOCATE( tg_psic( v_siz, npol ) )
     CALL stop_clock ('vloc_psi:tg_gather')

     incr = fftx_ntgrp(dffts)
     !
  ELSE
     !
     ALLOCATE(dffts_nl(1:dffts%ngm))
     !
     dffts_nl = dffts%nl
     v_siz = dffts%nnr
     v_siz_p = dfftp%nnr
     !
#if defined(__OPENMP_GPU)
     !$omp target enter data map(to:v,igk_k,dffts_nl,psi,hpsi)
#endif
     !
  ENDIF
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  DO ibnd = 1, m, incr

     IF( use_tg ) THEN
        !
        CALL tg_get_nnr( dffts, right_nnr )
        !
        DO ipol = 1, npol
           !
           tg_psic(:,ipol) = ( 0.D0, 0.D0 )
           ioff   = 0
           !
           DO idx = 1, fftx_ntgrp(dffts)
              !
              IF( idx + ibnd - 1 <= m ) THEN
                 DO j = 1, n
                    tg_psic( dffts%nl( igk_k(j,current_k) ) + ioff, ipol ) = &
                       psi( j +(ipol-1)*lda, idx+ibnd-1 )
                 ENDDO
              ENDIF

              ioff = ioff + right_nnr

           ENDDO
           !
           CALL invfft ('tgWave', tg_psic(:,ipol), dffts )
           !
        ENDDO
        !
     ELSE
#if defined(__OPENMP_GPU)
        !$omp target teams distribute parallel do collapse(2)
        DO ipol=1, npol
          DO j = 1, v_siz_p
              psic_nc_omp(j,ipol) = (0.d0,0.d0)
          END DO
        END DO
#else
        psic_nc_omp = (0.d0,0.d0)
#endif
        DO ipol=1,npol
        !$omp target teams distribute parallel do
           DO j = 1, n
              psic_nc_omp(dffts_nl(igk_k(j,current_k)),ipol) = psi(j+(ipol-1)*lda,ibnd)
           ENDDO
#if defined(__OMP_DISPATCH)
           !$omp dispatch
#endif
           CALL invfft ('Wave', psic_nc_omp(:,ipol), dffts)
        ENDDO
     ENDIF

     !
     !   product with the potential v = (vltot+vr) on the smooth grid
     !
     IF( use_tg ) THEN
        CALL tg_get_group_nr3( dffts, right_nr3 )
        IF (domag) THEN
           DO j=1, dffts%nr1x*dffts%nr2x*right_nr3
              sup = tg_psic(j,1) * (tg_v(j,1)+tg_v(j,4)) + &
                    tg_psic(j,2) * (tg_v(j,2)-(0.d0,1.d0)*tg_v(j,3))
              sdwn = tg_psic(j,2) * (tg_v(j,1)-tg_v(j,4)) + &
                     tg_psic(j,1) * (tg_v(j,2)+(0.d0,1.d0)*tg_v(j,3))
              tg_psic(j,1)=sup
              tg_psic(j,2)=sdwn
           ENDDO
        ELSE
           DO j=1, dffts%nr1x*dffts%nr2x*right_nr3
              tg_psic(j,:) = tg_psic(j,:) * tg_v(j,1)
           ENDDO
        ENDIF
     ELSE
        IF (domag) THEN
           !$omp target teams distribute parallel do
           DO j=1, v_siz
              sup = psic_nc_omp(j,1) * (v(j,1)+v(j,4)) + &
                    psic_nc_omp(j,2) * (v(j,2)-(0.d0,1.d0)*v(j,3))
              sdwn = psic_nc_omp(j,2) * (v(j,1)-v(j,4)) + &
                     psic_nc_omp(j,1) * (v(j,2)+(0.d0,1.d0)*v(j,3))
              psic_nc_omp(j,1)=sup
              psic_nc_omp(j,2)=sdwn
           ENDDO
        ELSE
#if defined(__OPENMP_GPU)
           !$omp target teams distribute parallel do collapse(2)
           DO ipol=1, npol
              DO j=1, v_siz
                 psic_nc_omp(j,ipol) = psic_nc_omp(j,ipol) * v(j,1)
              END DO
           ENDDO
#else
           DO j=1, v_siz
                 psic_nc_omp(j,:) = psic_nc_omp(j,:) * v(j,1)
           END DO
#endif
        ENDIF
     ENDIF
     !
     !   back to reciprocal space
     !
     IF( use_tg ) THEN
        !
        DO ipol = 1, npol

           CALL fwfft ('tgWave', tg_psic(:,ipol), dffts )
           !
           ioff   = 0
           !
           CALL tg_get_recip_inc( dffts, right_inc )
           !
           DO idx = 1, fftx_ntgrp(dffts)
              !
              IF( idx + ibnd - 1 <= m ) THEN
                 DO j = 1, n
                    hpsi (j, ipol, ibnd+idx-1) = hpsi (j, ipol, ibnd+idx-1) + &
                                 tg_psic( dffts%nl(igk_k(j,current_k)) + ioff, ipol )
                 ENDDO
              ENDIF
              !
              ioff = ioff + right_inc
              !
           ENDDO

        ENDDO
        !
     ELSE
        !
        DO ipol=1,npol
#if defined(__OPENMP_GPU)
        !$omp dispatch
#endif
           CALL fwfft ('Wave', psic_nc_omp(:,ipol), dffts)
        ENDDO
        !
        !   addition to the total product
        !
        !$omp target teams distribute parallel do collapse(2)
        DO ipol=1,npol
           DO j = 1, n
              hpsi(j,ipol,ibnd) = hpsi(j,ipol,ibnd) + &
                                  psic_nc_omp(dffts_nl(igk_k(j,current_k)),ipol)
           ENDDO
        ENDDO

     ENDIF

  ENDDO

  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_v )
     DEALLOCATE( tg_psic )
     !
  ELSE
     !
#if defined(__OPENMP_GPU)
     !$omp target exit data map(delete:dffts_nl,igk_k,v,psi) map(from:hpsi)
#endif
     DEALLOCATE( dffts_nl )
     !
  ENDIF
  CALL stop_clock ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_nc
