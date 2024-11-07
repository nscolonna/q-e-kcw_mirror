!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
#define ONE (0.D0,1.D0)
!#define DEBUG
!-----------------------------------------------------------------------
SUBROUTINE rho_of_q_gamma (rhowann, ngk_all, igk_k_all)
  !-----------------------------------------------------------------------
  !
  !! This sunroutine compute the periodic part of the Wannier density 
  !! rho(r) = \sum_q exp[iqr]rho_q(r) 
  !! rho_q(r) = \sum_k u^*_{k,n}(r) u_{k+q,n}(r)
  !! The k+q is rewritten as p+G with p in the 1BZ, then u_{k+q,n}(r) = u_p(r)exp{-iGr}
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE klist,                ONLY : nkstot, xk, igk_k, ngk, nks
  USE mp,                   ONLY : mp_sum
  USE control_kcw,          ONLY : evc0, iuwfc_wann, iuwfc_wann_allk, spin_component, num_wann
  USE buffers,              ONLY : get_buffer, save_buffer
  USE wvfct,                ONLY : npwx !, wg
  USE noncollin_module,     ONLY : npol,nspin_mag
  USE control_kcw,          ONLY : map_ikq, shift_1bz, nrho
  USE cell_base,            ONLY : at
  USE mp_pools,             ONLY : inter_pool_comm
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE wavefunctions,        ONLY : psic
  USE fft_wave,             ONLY : wave_g2r
  ! USE mp_world,             ONLY : mpime
  ! 
#ifdef DEBUG
  USE gvect,                ONLY : g, ngm
#endif
  !
  IMPLICIT NONE
  ! 
  INTEGER :: ik, ikq, npw
  !! Counter for the k/q points in the BZ, total number of q points and number of pw 
  !
  INTEGER :: iband, lrwfc, nkstot_eff
  !! Band counter
  !
  COMPLEX(DP), INTENT(OUT) :: rhowann(dffts%nnr, num_wann,nrho)
  !! The periodic part of the wannier orbital density
  !
  COMPLEX(DP) ::  evc0_kq(npwx*npol, num_wann)
  !! Auxiliary vector to store the wfc at k+q
  !
  REAL(DP) :: g_vect(3), xk_(3)
  !! G vector that shift the k+q inside the 1BZ
  !
  COMPLEX(DP) :: evc_k_g (npwx*npol), evc_k_r (dffts%nnr,npol), phase(dffts%nnr)
  !! Auxiliary wfc in reciprocal and real space, the phase associated to the hift k+q-> k'
  !
  COMPLEX(DP) :: evc_kq_g (npwx*npol), evc_kq_r (dffts%nnr,npol)
  !! Auxiliary wfc in reciprocal and real space
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !! The global index of k-points
  !
  INTEGER :: global_ik, ik_eff, ip, ipp
  !
  INTEGER, INTENT(IN) :: &
       igk_k_all(npwx,nkstot),&    ! index of G corresponding to a given index of k+G
       ngk_all(nkstot)             ! number of plane waves for each k point
  !
  INTEGER :: ebnd
  REAL(DP) :: w1, w2
  !
#ifdef DEBUG
  INTEGER :: ig_save, ig, ip, ir
  COMPLEX (DP ) :: pippo
  REAL(DP) :: pippo_real, xq(3)
#endif
  !  
  CALL start_clock ( 'rho_of_q' )
  IF (nspin == 4) THEN
    nkstot_eff = nkstot
  ELSE
    nkstot_eff = nkstot/nspin
  ENDIF
  DO ik = 1, nks
    ! CHECK: Need to understand/think more about pool parallelization
    ! what happen if k+q is outside the pool??
    ! Some problem in EXX (have a look to /PW/src/exx.f90 exx_grid_init)
    ! Solved for now with WFs Broadcast (See bcast_wfc.f90) 
    !
    IF ( lsda ) current_spin = isk(ik)
    IF ( lsda .AND. current_spin /= spin_component) CYCLE
    !
    xk_ = xk(:,ik)
    CALL cryst_to_cart(1, xk_, at, -1)
    !
#ifdef DEBUG
     WRITE(*,'(10x, "ik = ", i5, 3x, "xk = ", 3f12.6, " [Cryst]")') ik, (xk_(ip), ip=1,3) 
#endif
    !
    global_ik = global_kpoint_index (nkstot,ik)
    global_ik = global_ik - (spin_component -1)*nkstot_eff
    lrwfc = num_wann*npwx*npol
    CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik )
    !! ... Retrive the ks function (in the Wannier Gauge)
    !
    npw = ngk(ik)
    !
    DO iband = 1, num_wann, 2
       !
       ebnd = iband
       IF ( iband < num_wann ) ebnd = ebnd + 1
       !
       CALL wave_g2r( evc0(1:npw,iband:ebnd), psic, dffts )
       !! ... The wfc in R-space at k
       !
       w1 = 1.D0 
       w2 = w1
       !
       rhowann(:,iband,1) =  w1 * DBLE ( psic(:) )**2 
       rhowann(:,ebnd,1)  =  w2 * AIMAG( psic(:) )**2 
       !
    ENDDO ! bands
    ! 
  ENDDO ! kpoints
  !
  CALL mp_sum( rhowann, inter_pool_comm )
  ! ... Sum up the results over k point in different pools
  !
  CALL stop_clock ('rho_of_q')
  !
END subroutine
