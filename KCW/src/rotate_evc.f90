SUBROUTINE rotate_evc(isym, gvector, psi_r)
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
   USE control_kcw,     ONLY : rir, rvect
   IMPLICIT NONE
   !input -> wf on regular grid
   !output -> wf rotated 
   !
   INTEGER, INTENT(IN)       :: isym
   REAL(DP), INTENT(IN)      :: gvector(3)
   COMPLEX(DP), INTENT(INOUT):: psi_r(dffts%nnr)
   !
   INTEGER                   :: ir
   COMPLEX(DP)               :: imag = (0.D0,1.D0)
   INTEGER                   :: nxxs ! number of points in r grid
   COMPLEX(DP), ALLOCATABLE :: psic_all(:), temppsic_all(:)
   COMPLEX(DP), ALLOCATABLE :: phase(:)
   !
   nxxs = dffts%nr1x *dffts%nr2x *dffts%nr3x
   WRITE(*,*) dffts%nnr
   ALLOCATE( phase(dffts%nnr))
   !
   !
   IF (isym > 1) THEN
#if defined(__MPI)
     ALLOCATE( psic_all(nxxs), temppsic_all(nxxs))
     ! gather among all the CPUs
     CALL gather_grid(dffts, psi_r, temppsic_all)
     ! apply rotation
     psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
     ! scatter back a piece to each CPU
     CALL scatter_grid(dffts, psic_all, psi_r)
     DEALLOCATE( psic_all, temppsic_all)
#else
     psi_r(1:nxxs) = psi_r(rir(1:nxxs,isym))
#endif
   ENDIF
   !
   ! Apply phase factor exp[-iG.r]
   CALL calculate_phase (gvector, phase)
   psi_r(:) = psi_r(:) * phase(:)

   DEALLOCATE ( phase )
   !
   !   
END SUBROUTINE rotate_evc
