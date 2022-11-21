SUBROUTINE rotate_evc(isym, phase, psi_r)
   !-----------------------------------------------------------------------
   ! g psi_k = e^{ik1.rS} u_k(rS) = e^{iSk1.r} u_k(rS)
   !         = e^{ik2.r} [ e^{-iGr} u_k(rS) ]
   ! k2 = s.k1 + G
   ! u_S(k)(r) = psi_k(S^-1(r))
   ! S=s(:,:,isym) 
   !
   USE kinds,           ONLY : DP
   USE wvfct,           ONLY : nbnd, npwx
   USE wavefunctions,   ONLY : evc, psic, psic_nc
   USE fft_base,        ONLY : dffts, dfftp
   USE fft_interfaces,  ONLY : fwfft, invfft
   USE cell_base,       ONLY : bg
   USE constants,       ONLY : tpi
   USE gvect,           ONLY : g, ngm
   USE mp,              ONLY : mp_sum
   USE mp_pools,        ONLY : intra_pool_comm
   USE fft_interfaces,  ONLY : invfft
   USE scatter_mod,     ONLY : gather_grid, scatter_grid
   USE control_kcw,     ONLY : rir
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN):: isym
   COMPLEX(DP), INTENT(INOUT):: psi_r(dffts%nnr)
   COMPLEX(DP), INTENT(IN)   :: phase(dffts%nnr)
   !
   !INTEGER:: ig, igk_local, igk_global, npw1, npw2, n, nxxs, ipol, istart, isym0
   INTEGER:: nxxs
   REAL(DP)                 :: g_vect_cart(3), srt(3,3)
   COMPLEX(DP), ALLOCATABLE :: psic_all(:), temppsic_all(:)
   !
   nxxs = dffts%nr1x *dffts%nr2x *dffts%nr3x
   ALLOCATE( psic_all(nxxs), temppsic_all(nxxs))
   !
   IF (isym > 1) THEN
#if defined(__MPI)
     ! gather among all the CPUs
     CALL gather_grid(dffts, psi_r, temppsic_all)
     ! apply rotation
     psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
     ! scatter back a piece to each CPU
     CALL scatter_grid(dffts, psic_all, psi_r)
#else
     psi_r(1:nxxs) = psi_r(rir(1:nxxs,isym))
#endif
   ENDIF
   !
   ! Apply phase factor exp[-iG.r]
   psi_r = psi_r * phase 
   !
   DEALLOCATE( psic_all, temppsic_all)
   !   
END SUBROUTINE rotate_evc
