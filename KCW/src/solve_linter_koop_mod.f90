! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"

MODULE solve_linter_koop_mod

CONTAINS
!
!-----------------------------------------------------------------------
subroutine solve_linter_koop ( spin_ref, i_ref, delta_vr, drhog_scf, delta_vg, drhor_scf)
  !-----------------------------------------------------------------------
  !
  !!    Driver routine for the solution of the linear system which defines the 
  !!    change of the wavefunction due to an external perturbation. It work 
  !!    exactly as solve_linter in PH but for a genenral perturbation dv delta_vr.
  !!    In genereal, it performs the following tasks:
  !!     a) computes the bare potential term Delta V | psi >, add the scf contribution
  !!     b) applies P_c^+ (orthogonalization to valence states)
  !!     c) calls cgstabsolve_all to solve the linear system
  !!     d) computes Delta rho
  !
  ! ## FIXME delta_vg passed only for debug reason (to test symmetries) Needs to be removed once the
  !          problem will be solved
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout
  USE wavefunctions,         ONLY : evc, psic
  USE klist,                 ONLY : lgauss, igk_k, ngk
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE fft_base,              ONLY : dffts, dfftp
  USE fft_interfaces,        ONLY : fwfft, invfft, fft_interpolate
  USE gvect,                 ONLY : gstart
  USE gvecs,                 ONLY : doublegrid, ngms
  USE becmod,                ONLY : calbec
  USE wvfct,                 ONLY : npw, npwx, nbnd
  USE uspp_param,            ONLY : nhm
  USE control_lr,            ONLY : lgamma, nbnd_occ
  USE units_lr,              ONLY : iuwfc, lrwfc
  USE buffers,               ONLY : save_buffer, get_buffer
  USE eqv,                   ONLY : dvpsi
  USE qpoint,                ONLY : npwq, nksq, ikks, ikqs
  USE uspp,                  ONLY : okvan  
  USE uspp_init,             ONLY : init_us_2
  USE mp,                    ONLY : mp_sum
  USE mp_global,             ONLY : intra_pool_comm, inter_pool_comm 
  USE noncollin_module,      ONLY : npol
  USE control_kcw
  USE dv_of_drho_lr
  USE mp_pools,              ONLY : inter_pool_comm, intra_pool_comm
  USE response_kernels,      ONLY : sternheimer_kernel
  !
  !USE cell_base,            ONLY : omega
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: npe = 1
  !
  COMPLEX(dp), INTENT(in) ::delta_vr (dffts%nnr, nspin)
  COMPLEX(dp), INTENT(out) ::  drhog_scf (ngms, nspin)
  COMPLEX(dp), INTENT(out), OPTIONAL :: drhor_scf(dffts%nnr,nspin)
  !
  ! input: the imaginary frequency
  REAL(DP) :: averlt
  !
  COMPLEX(DP), ALLOCATABLE :: drhoscf (:,:,:)
  REAL(DP) :: thresh, dr2
  ! thresh: convergence threshold
  ! dr2   : self-consistency error
  REAL(DP) :: dos_ef
  ! dos_ef: DOS at Fermi energy (in calculation for a metal)
  COMPLEX(DP), ALLOCATABLE, TARGET :: dvscfin(:,:,:)
  ! change of the scf potential 
  COMPLEX(DP), POINTER :: dvscfins (:,:,:)
  ! change of the scf potential (smooth part only)
  COMPLEX(DP), ALLOCATABLE :: drhoscfh (:,:,:), dvscfout (:,:,:)
  ! change of rho / scf potential (output)
  ! change of scf potential (output)
  COMPLEX(DP), ALLOCATABLE :: ldos (:,:), ldoss (:,:),&
       dbecsum (:,:,:,:), aux(:), aux1 (:,:), aux2(:,:)
  COMPLEX(DP), ALLOCATABLE :: drhoc(:)
  ! drhoc: response core charge density
  REAL(DP), ALLOCATABLE :: becsum1(:,:,:)
  !
  LOGICAL :: lmetq0        ! true if xq=(0,0,0) in a metal
  LOGICAL :: all_conv

  INTEGER :: kter,       & ! counter on iterations
             iter0,      & ! starting iteration
             ibnd,       & ! counter on bands
             iter,       & ! counter on iterations
             lintercall, & ! average number of calls to cgsolve_all
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             ltaver,     & ! average counter
             ig,         & ! counter on G vectors
             ir,         & ! counter on mesh points
             is,         & ! counter on spin polarizations
             ipert,      & ! counter on perturbations
             spin_ref,   & ! the spin of the reference orbital
             i_ref         ! the orbital we want to keep fix
           


  REAL(DP) :: tcpu, get_clock ! timing variables

  CHARACTER(LEN=256) :: flmixdpot = 'mixd'
  LOGICAL :: convt 
  !
  !!## DEBUG
  COMPLEX(DP) :: delta_vg(ngms,nspin)
  !!## DEBUG 
  !
  CALL start_clock ('solve_linter')
  !
  ALLOCATE (dvscfin (dfftp%nnr, nspin , npe))    
  IF (doublegrid) THEN
     ALLOCATE (dvscfins (dffts%nnr , nspin , npe))    
  ELSE
     dvscfins => dvscfin
  ENDIF
  ALLOCATE (drhoscf  (dffts%nnr, nspin, 1) ) !! NsC
  ALLOCATE (drhoscfh (dfftp%nnr, nspin , npe))    
  ALLOCATE (dvscfout (dfftp%nnr, nspin , npe))    
  ALLOCATE (dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin , npe))    
  ALLOCATE (aux ( dffts%nnr ))    
  ALLOCATE (aux1 ( dffts%nnr, npol ))    
  ALLOCATE (aux2(npwx*npol, nbnd))
  ALLOCATE (drhoc(dfftp%nnr))
  !
  IF (kcw_at_ks .AND. fix_orb) WRITE (stdout, '("FREEZING ORBITAL #", i4, 3x , "spin", i4)') i_ref, spin_ref
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  lmetq0 = lgauss .AND. lgamma
  IF (lmetq0) THEN
     ALLOCATE ( ldos ( dfftp%nnr, nspin) )    
     ALLOCATE ( ldoss( dffts%nnr, nspin) )    
     ALLOCATE (becsum1 ( (nhm * (nhm + 1))/2 , nat , nspin))
     call localdos ( ldos , ldoss , becsum1, dos_ef )
  ENDIF
  ! 
  DO ik = 1, nksq 
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     !
     IF (lsda) current_spin = isk (ikk)
     !
     ! read unperturbed wavefunctions psi(k) and psi(k+q)
     !
     IF (nksq .GT. 1) call get_buffer (evc, lrwfc, iuwfc, ikk)
     ! 
     ! ... Compute dv_bare*psi and set dvscf to zero
     !
     dvpsi(:,:) = (0.D0, 0.D0)
     DO ibnd = 1, nbnd_occ (ik)
        aux(:) = (0.d0, 0.d0)
        DO ig = 1, npw
           aux (dffts%nl(igk_k(ig,ikk)))=evc(ig,ibnd)
        ENDDO
        CALL invfft ('Wave', aux, dffts)
        DO ir = 1, dffts%nnr
            aux(ir)=aux(ir)*delta_vr(ir,current_spin) 
        ENDDO
        !
        CALL fwfft ('Wave', aux, dffts)
        DO ig = 1, npwq
           dvpsi(ig,ibnd)=aux(dffts%nl(igk_k(ig,ikq)))
        ENDDO
        !
     ENDDO
     !
     IF (okvan) THEN
        call errore('solve_linter_koop', 'USPP not implemented yet', 1)
     ENDIF
     !
     CALL save_buffer (dvpsi, lrdvwfc, iudvwfc, ik)
     !
  ENDDO
  !
  !   The outside loop is over the iterations
  !
  dr2=0.d0
  iter0 = 0 ! We do not have a restart procedure yet: NsC
  !
  DO kter = 1, niter
     iter = kter + iter0
     !
     ltaver = 0
     !
     lintercall = 0
     drhoscf(:,:,:) = (0.d0, 0.d0)
     dvscfout(:,:,:)    = (0.d0, 0.d0)
     dbecsum(:,:,:,:) = (0.d0, 0.d0)
     !
     IF (iter == 1 ) THEN 
        thresh = 1.d-6
     ELSE 
        thresh = min (1.d-2 * sqrt (dr2), 1.d-6)
     ENDIF
     !
     CALL sternheimer_kernel(iter==1, .FALSE., 1, lrdvwfc, iudvwfc, &
         thresh, dvscfins, all_conv, averlt, drhoscf, dbecsum, exclude_hubbard=.TRUE.)
     !
     IF ((.NOT. all_conv) .AND. (iter == 1)) THEN
        WRITE(stdout, '(6x, "sternheimer_kernel not converged. Try to increase thresh_init.")')
     ENDIF
     !
#ifdef __MPI
     !
     !  The calculation of dbecsum is distributed across processors (see addusdbec)
     !  Sum over processors the contributions coming from each slice of bands
     !
     call mp_sum (dbecsum, intra_pool_comm)
     ! 
#endif
     !
     IF (doublegrid) THEN
        DO is = 1, nspin
           call fft_interpolate (dffts, drhoscf(:,is,1), dfftp, drhoscfh(:,is,1))
        ENDDO
     ELSE
        CALL zcopy (npe*nspin*dfftp%nnr, drhoscf, 1, drhoscfh, 1)
     ENDIF
     !
     ! if q=0, make sure that charge conservation is guaranteed
     !
     IF ( lgamma ) THEN
        psic(:) = drhoscfh(:, nspin, npe)
        CALL fwfft ('Rho', psic, dfftp)
        !CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
        IF ( gstart==2) psic(dfftp%nl(1)) = (0.d0, 0.d0)
        CALL invfft ('Rho', psic, dfftp)
        !CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
        drhoscfh(:, nspin, npe) = psic(:)
     ENDIF
     !
     !    Now we compute for all perturbations the total charge and potential
     !
     !call addusddens (drhoscfh, dbecsum, irr, imode0, npe, 0)
     !
#ifdef __MPI
     !
     !   Reduce the delta rho across pools
     !
     call mp_sum (drhoscf, inter_pool_comm)
     call mp_sum (drhoscfh, inter_pool_comm)
     !
#endif
     !
     !   ... save them on disk and
     !   compute the corresponding change in scf potential
     !
     DO ipert = 1, npe
        !
        CALL zcopy (dfftp%nnr*nspin,drhoscfh(1,1,ipert),1,dvscfout(1,1,ipert),1)
        ! NB: always call with imode=0 to avoid call to addcore in dv_of_drho for 
        !     nlcc pseudo. The call is not needed since we are not moving atoms!!
        !
        CALL dv_of_drho (dvscfout(1,1,ipert), .false.)
     ENDDO
     !
     !
     ! ... On output in dvscfin we have the mixed potential
     !
     CALL mix_potential (2*npe*dfftp%nnr*nspin, dvscfout, dvscfin, &
                         alpha_mix(kter), dr2, npe*tr2/npol, iter, &
                         nmix, flmixdpot, convt)
     !WRITE(mpime+1000, '(1i5,es10.3,1l1,1i5)') my_pool_id, dr2, convt, iter
     !
     ! check that convergent have been reached on ALL processors in this image
     CALL check_all_convt(convt)

     IF (doublegrid) THEN
        DO ipert = 1, npe
           DO is = 1, nspin
              !call cinterpolate (dvscfin(1,is,ipert), dvscfins(1,is,ipert), -1)
              call fft_interpolate (dffts, drhoscf(:,is,ipert), dfftp, drhoscfh(:,is,ipert))
           ENDDO
        ENDDO
     ENDIF
     !
     tcpu = get_clock ('KCW')

     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / npe
     WRITE( stdout, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
          &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (kter) , dr2
     !
     !    Here we save the information for recovering the run from this poin
     ! 
     FLUSH( stdout )
     !
     IF (convt) GOTO 155
     !
  ENDDO  ! loop over iteration
  !
155 iter0=0
  !
  ipert =1  
  ! The density variation in G-space
  !
  DO is = 1, nspin
     !aux(:) = drhoscf(:,is,ipert)
     aux(:) = drhoscfh(:,is,ipert)
     CALL fwfft ('Rho', aux, dffts)
     drhog_scf(:,is) = aux(dffts%nl(:))
  ENDDO
  !
  IF(present(drhor_scf)) drhor_scf = drhoscf(:,:,ipert)
  !
  ! The induced density in G space
  !
  IF (lmetq0) DEALLOCATE (ldoss)
  IF (lmetq0) DEALLOCATE (ldos)
  DEALLOCATE (aux)
  DEALLOCATE (aux1)
  DEALLOCATE (aux2)
  DEALLOCATE (drhoc)
  DEALLOCATE (dbecsum)
  DEALLOCATE (drhoscf )
  DEALLOCATE (dvscfout)
  DEALLOCATE (drhoscfh)
  IF (doublegrid) DEALLOCATE (dvscfins)
  DEALLOCATE (dvscfin)
  !
  CALL stop_clock ('solve_linter')
  !
  RETURN
  !
END SUBROUTINE solve_linter_koop

!------------------------------------------------------------------
SUBROUTINE check_all_convt( convt )
  !---------------------------------------------------------------
  !! Work out how many processes have converged.
  !
  USE mp,        ONLY : mp_sum
  USE mp_images, ONLY : nproc_image, intra_image_comm
  !
  IMPLICIT NONE
  !
  LOGICAL,INTENT(in) :: convt
  INTEGER            :: tot_conv
  !
  IF(nproc_image==1) RETURN
  !
  tot_conv = 0
  IF(convt) tot_conv = 1
  CALL mp_sum(tot_conv, intra_image_comm)
  !
  IF ((tot_conv > 0) .and. (tot_conv < nproc_image)) THEN
    CALL errore('check_all_convt', 'Only some processors converged: '&
               &' either something is wrong with solve_linter, or a different'&
               &' parallelism scheme should be used.', 1)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE

END MODULE
