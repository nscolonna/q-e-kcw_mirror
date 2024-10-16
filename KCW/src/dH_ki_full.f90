!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
!#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!-----------------------------------------------------------------------
SUBROUTINE dH_ki_full (ik, dH_wann)
  !-----------------------------------------------------------------------
  !
  ! This routine compute the full KI or the pKIPZ correction to the spectrum
  ! KI potential:   See Eq. A8-A11 in PRB 90, 075135 (2014).
  ! KIPZ potential: See Eq. A15 in PRB 90, 075135 (2014).
  ! NB: Only the insulating case -i.e. fi=1 or fi=0
  ! 
  ! KI potenital:
  ! First the diagonal corrections are computed (i==j):

  !                     | -E_H[n_i] + E_xc[\rho] - E_xc[\rho-n_i] -\int v_xc[\rho](r) n_i(r)  if f_i=1 (OCC) (Eq. A8)
  ! <i | v^{KI}_i |i> = | 
  !                     | -E_H[n_i] - E_xc[\rho] + E_xc[\rho+n_i] -\int v_xc[\rho](r) n_i(r)  if f_i=0 (EMP) (Eq. A8 and A10)
  !
  ! Then possibly the off-diagonal (i /= j) (ONLY apply to EMPTY states) 
  ! (The scalar term C does not contribute as <i|C|j> = 0 if i/=j):
  !
  ! <i | v^{KI}_i |j> = \int \phi*_i(r) [v_H[n_i](r) + v_xc[\rho+n_i](r) - v_xc[\rho](r) ] \phi_j(r)  if (f_i=f_j=0) (Eq. A10)
  !
  ! pKIPZ potenital adds to KI potential:
  !
  ! First the diagonal corrections are computed (i==j):
  ! <i | Dv^{PZ}_i |i> = -E_{Hxc}[n_i]  (Eq. A15)
  !
  ! Then possibly the off-diagonal (i /= j) 
  ! (The scalar term C does not contribute as <i|C|j> = 0 if i/=j):
  ! <i | Dv^{PZ}_i |j> = \int \phi*_i(r) [ -v_{Hxc}[n_i](r)] \phi_j(r)
  !
  ! NB: The Hamiltonian is not Hermitian in general because we do not enforce the pederson condition 
  !     whien using KS or MLWFs as variational orbitals. We force it to be hermitina in two ways:
  !     1) Just keep the upper part of the Hamiltonian and build the lower part by CC (this is not actually done 
  !        because in any case the diagonalization (cdiagh in ../PW/src/cdiagh.f90) implicitely assume H to be hermitian
  !        and only consider the upper half of H).
  !     2) We define an hermitian operator as follow 
  !        H_{ij} = 0.5*(<i|hj|j> + <i|h_i|j>)  = 
  !               = 0.5*(  h^R_ij  + h^L_ij  )  = 
  !               = 0.5*(<i|hj|j> + <j|h_i|i>*) =
  !               = 0.5*(h^R_{ij} + [h^R_ji]*)  = 
  !     --> H = 0.5*(h^R+CONGJ(TRANSPOSE(h^R))
  !     where h^R(h^L) is the Right(Left) hamiltonian: h^R_{ij}= <i | h_j | j> (h^L_{ij}= <i | h_i | j>)
  !     H is Hermitian by construction H_ij = [H_ji]*
  !
  !     IF the pederson condition would have been satisfied (i.e.  <i |v_j |j > = <i |v_i|j>), then H == h^L == h^R.
  !     This is the case for the KI hamiltonian for the occupied manifold (as the potential is simply a scalar and the 
  !     Pederson condition is trivially satisfied), or the case of the on_site_only apporximation (as we only keep the 
  !     diagonal term, i.e.  <j|h_i|i> = <j|h_j|i> ==0 for each i/=j.)
  ! 
  !
  USE wavefunctions,         ONLY : psic
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE wvfct,                 ONLY : npw, nbnd, npwx
  USE fft_base,              ONLY : dffts, dfftp
  USE lsda_mod,              ONLY : lsda, current_spin,nspin
  USE klist,                 ONLY : ngk, init_igk, igk_k, nkstot
  USE gvect,                 ONLY : ngm
  USE gvecs,                 ONLY : ngms
  USE buffers,               ONLY : get_buffer
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE control_kcw,           ONLY : kcw_at_ks, homo_only, alpha_final, &
                                    num_wann_occ, iuwfc_wann, kcw_iverbosity, &
                                    qp_symm, evc0, kipz_corr, num_wann, &
                                    spin_component, l_alpha_corr, on_site_only
  USE control_lr,            ONLY : lrpa
  USE mp,                    ONLY : mp_sum
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE cell_base,             ONLY : omega
  USE scf,                   ONLY : rho, rho_core, rhog_core, create_scf_type, scf_type
  USE constants,             ONLY : rytoev
  !
  USE exx,                   ONLY : vexx, exxinit, exxenergy2
  USE input_parameters,      ONLY : exxdiv_treatment
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT (OUT) :: dH_wann(nkstot/nspin,num_wann,num_wann)
  !
  INTEGER, INTENT(IN) :: ik
  ! the k-point index in hte orignal mesh of k-points
  !
  INTEGER :: ibnd, ir, k, dim_ham, nspin_aux, lrwfc
  !
  REAL(DP) :: n_r(dfftp%nnr), num1, num2, sh, aux_r(dfftp%nnr) 
  ! ... orbital density in rela space
  !
  COMPLEX(DP) n_g(ngm), n_g_aux(ngm,nspin), aux_g(ngm)
  ! ... orbital density in G space
  !
  REAL(DP) :: ehart, v(dfftp%nnr,nspin), vxc_minus1(dfftp%nnr,2), vxc(dfftp%nnr,nspin), charge, w1
  !  
  TYPE (scf_type) :: rho_minus1
  !
  REAL(DP) ::  vtxc, etxc, vtxc_minus1, etxc_minus1, etmp, delta_eig(nbnd)
  !
  COMPLEX(DP) , ALLOCATABLE :: psic_1(:) 
  COMPLEX(DP) , ALLOCATABLE :: ham (:,:), ham_left(:,:), ham_right(:,:), vpsi_r(:), vpsi_g(:)
  !
  COMPLEX(DP) :: delta_vr(dffts%nnr,nspin), delta_vr_(dffts%nnr,nspin)
  COMPLEX(DP), ALLOCATABLE  :: delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
  COMPLEX(DP) :: deltah_scal
  REAL(DP) :: ddH
  COMPLEX(DP) :: dhkipz
  INTEGER :: segno
  !
  !
  ALLOCATE (psic_1( dfftp%nnr), vpsi_r(dffts%nnr), vpsi_g(npwx))
  ALLOCATE ( delta_vg(ngms,nspin), vh_rhog(ngms), delta_vg_(ngms,nspin) )
  !
  !
  IF (kipz_corr) THEN 
     WRITE(stdout,'(/,5x,"INFO: pKIPZ calcualtion ... ",/ )')
  ELSE
     WRITE(stdout,'(/,5x,"INFO: KI calcualtion: Full Hamiltonian ... ",/ )')
  ENDIF
  IF (on_site_only) WRITE(stdout, '(/,5X, "INFO: Skipping off-diag: only R=0 and i=j")') 
  !
  w1 = 1.D0 / omega
  !
  nspin_aux=nspin
  nspin=2
  CALL create_scf_type (rho_minus1)
  nspin=nspin_aux
  !
  lrwfc = num_wann*npwx
  CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik )
  ! Retrive the variational orbital at kpoint k (includes spin)  
  IF (kcw_iverbosity .gt. 1 ) WRITE(stdout,'(8X, "INFO: u_k(g) RETRIEVED"/)')
  !
  CALL compute_map_ikq_single (ik)
  ! find tha map k+q --> k'+G and store the res 
  !
  dim_ham = num_wann
  ALLOCATE ( ham (dim_ham,dim_ham) )
  ALLOCATE ( ham_left (dim_ham,dim_ham) )
  ALLOCATE ( ham_right (dim_ham,dim_ham) )
  ham (:,:) = (0.D0,0.D0)
  ham_left (:,:) = (0.D0,0.D0)
  ham_right (:,:) = (0.D0,0.D0)
  !
  orb_loop: DO ibnd = 1, num_wann
     !
     vpsi_g(:) = (0.D0,0.D0)
     ! keep track of occ ...
     segno = -1 
     ! ... vs empty states
     IF (ibnd .gt. num_wann_occ) segno = +1
     !
     IF (kcw_at_ks .AND. homo_only .AND. (ibnd .ne. num_wann_occ) ) CYCLE orb_loop ! only homo orbilal (for fast debug) 
     ! 
     IF ( lsda ) current_spin = spin_component
     ! 
     ! ... Compute the orbital density ...
     !
     npw = ngk(ik)
     psic(:) = ( 0.D0, 0.D0 )
     psic(dffts%nl(igk_k(1:npw,ik))) = evc0(1:npw,ibnd)
     CALL invfft ('Wave', psic, dffts)
     !
     ! Store the result needed below
     psic_1(:) = psic(:) 
     !
     ! ... orbital density in real space ...
     num1 = 0.D0
     num2 = 0.D0
     n_r(:) = 0.0
     DO ir = 1, dffts%nnr
       n_r(ir) = n_r(ir) + ( DBLE( psic(ir) )**2 + AIMAG( psic(ir) )**2 )*w1
#ifdef DEBUG
       num1 = num1 + DBLE( DBLE( psic(ir) )**2 + AIMAG( psic(ir) )**2  )*w1
       num2 = num2 + rho%of_r(ir,1)
#endif
     ENDDO 
     !
#ifdef DEBUG
     CALL mp_sum (num1, intra_bgrp_comm) 
     CALL mp_sum (num2, intra_bgrp_comm) 
     WRITE(stdout,'(8x, "orbital charge", 2F18.12)') num1/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
     WRITE(stdout,'(8x, "spin-up charge", 2F18.12)') num2/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
#endif
     !
     ! ... orbital density in reciprocal space ...
     !
     CALL bare_pot ( CMPLX(n_r*omega, 0.0, kind=DP), n_g, vh_rhog, delta_vr, delta_vg, 1, delta_vr_, delta_vg_ )
     !! The periodic part of the perturbation DeltaV_q(G)
     !
     deltah_scal = - 0.5D0 * sum (CONJG(n_g (:)) * delta_vg(:,spin_component)) * omega
     !sh = 0.5D0 * sum (CONJG(n_g (:)) * vh_rhog(:)                )*omega
     !WRITE(stdout,'(8x, "self_hatree", 2i5, 1F15.8)') ibnd, current_spin, sh
     !
     IF (l_alpha_corr) CALL beyond_2nd (deltah_scal, ddH, ibnd)
     !
     IF (alpha_final(ibnd) .gt. 1.02 ) THEN
        WRITE(stdout,'("WARNING: alpha for orbital", i5, i3, "  bigger than 1.02.", F15.8, "Set it to 1.00",/)') ibnd, &
            ik, alpha_final(ibnd)
        alpha_final(ibnd) = 1.D0
     ENDIF
     !
     IF (alpha_final(ibnd) .lt. 0.00 ) THEN
        WRITE(stdout,'(8x, "WARNING: alpha for orbital", i5, i3, "  smaller than 0.00.", F15.8, "Set it to 1.00",/)') ibnd, &
            ik, alpha_final(ibnd)
        alpha_final(ibnd) = 1.D0
     ENDIF
     !
     psic(:) = (0.D0, 0.D0)
     psic(:) =  CMPLX(n_r(:),0.D0,kind=dp)
     CALL fwfft ('Rho', psic, dfftp)
     n_g(:) = psic(dfftp%nl(:))
     !
     ! ... Compute Int[e^2/|r-r'|*n_i(r')] ...
     !
     sh = 0.D0
     v(:,:)=0.D0
     n_g_aux(:,:) = (0.D0, 0.D0)
     n_g_aux(:,1) = n_g(:)
     CALL v_h( n_g_aux, ehart, charge, v )
     sh = ehart
     !
!     WRITE(stdout,'("v_hatree", 2i5, 3F15.8)') ibnd, current_spin ,REAL(v(1:3,1))
     !WRITE(stdout,'(8x, "self_hatree", 2i5, 1F15.8)') ibnd, current_spin, sh
     !
#ifdef DEBUG
     WRITE(stdout,'(8x, "orbital=", i3, 2x, "Self-Hartree", F15.10, 3x, "Ry",/)') ibnd, sh
     WRITE(stdout,'(8x,"orbital charge from v_h",F15.12,/)') charge
#endif
     !
     ! .. Add the xc contribution ...
     !
     etxc = 0.D0; vtxc = 0.D0; vxc(:,:) = 0.D0
     CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )
     !
     etmp =0.D0
     etmp = sum ( vxc(1:dfftp%nnr,current_spin) * n_r(1:dfftp%nnr) )
     etmp = etmp/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
     CALL mp_sum (etmp, intra_bgrp_comm) 
     !
     IF (nspin == 1 ) THEN
        rho_minus1%of_r(:,1) = rho%of_r(:,1)
        rho_minus1%of_r(:,2) = 0.D0
        rho_minus1%of_g(:,1) = rho%of_g(:,1)
        rho_minus1%of_g(:,2) = 0.D0
     ELSE
        rho_minus1%of_r(:,:) = rho%of_r(:,:)
        rho_minus1%of_g(:,:) = rho%of_g(:,:)
     ENDIF
     !
     delta_eig(ibnd) = 0.D0
     !
     rho_minus1%of_r(:,1) = rho_minus1%of_r(:,1) + segno * n_r(:)  ! denisty rho \pm rho_i in real space
     rho_minus1%of_g(:,1) = rho_minus1%of_g(:,1) + segno * n_g(:)  ! denisty rho \pm rho_i in reciprocal scape
     ! The magnetization depends on which channel we remove the orbital from
     ! we reduce by n_i(r) if current_spin=1, we increase by n_i(r) if current_spin=2
     ! This is taken care by the factor (3-2*current_spin)
     rho_minus1%of_r(:,2) = rho_minus1%of_r(:,2) + segno * (3-2*spin_component)*n_r(:)  ! magnetization m-m_i in real space
     rho_minus1%of_g(:,2) = rho_minus1%of_g(:,2) + segno * (3-2*spin_component)*n_g(:)  ! magnetization m-m_i in reciprocal scape
     ! 
     etxc_minus1 = 0.D0; vtxc_minus1 = 0.D0; vxc_minus1(:,:) = 0.D0
     nspin_aux=nspin; nspin=2
     CALL v_xc( rho_minus1, rho_core, rhog_core, etxc_minus1, vtxc_minus1, vxc_minus1 )
     nspin=nspin_aux
     !
     delta_eig(ibnd) = (sh - etxc + etxc_minus1 - segno * etmp)
     delta_eig(ibnd) = segno * delta_eig(ibnd)
     !
     WRITE(stdout, '(8x, "KI corr const term, sh[n_i], Exc[n], Exc[n-n_i], int{v_xc[n] n_i} ", 4F14.8)') sh, &
         etxc, etxc_minus1, etmp 
     IF (lrpa) delta_eig(ibnd) = segno * sh  !! hartree only for debug
     WRITE(stdout,'(8X, "Delta KI", 2F15.8)')  delta_eig(ibnd), delta_eig(ibnd)*alpha_final(ibnd)
     !
     ham_right(ibnd,ibnd) = delta_eig(ibnd) * alpha_final(ibnd)
     !WRITE(stdout,*) ibnd, ibnd, REAL(ham_right(ibnd,ibnd)), AIMAG(ham_right(ibnd,ibnd))
     !
     ! Eventually the KI off-diagonal contribution. Only apply to empty manifold. 
     IF ( .NOT. on_site_only .AND. ibnd .GT. num_wann_occ) THEN
        !
        vpsi_r(:) = (0.D0, 0.D0)
        DO ir = 1, dffts%nnr
           vpsi_r (ir) = CMPLX( ( v(ir,current_spin) + vxc_minus1(ir,current_spin) - &
                   vxc(ir,current_spin) ),0.D0, kind=DP) * psic_1(ir)
           IF (lrpa) vpsi_r (ir) = CMPLX( v(ir,current_spin), 0.D0, kind=DP) * psic_1(ir)
        ENDDO
        !
        ! 1) GO to G-space and store the ki gradient. v_i|phi_i> 
        CALL fwfft ('Wave', vpsi_r, dffts)
        vpsi_g(:) = vpsi_r(dffts%nl(igk_k(:,ik)))
        !
        DO k = num_wann_occ+1, num_wann
           IF (k == ibnd) CYCLE
           ham_right(k,ibnd) = SUM ( CONJG(evc0(:,k)) * vpsi_g(:) * alpha_final(ibnd) ) 
           CALL mp_sum (ham_right(k,ibnd), intra_bgrp_comm)
           !WRITE(stdout,*) k, ibnd, REAL(ham_right(k,ibnd)), AIMAG(ham_right(k,ibnd))
        ENDDO 
        !
     ENDIF
     !
     ! pKIPZ
     IF (kipz_corr) THEN 
        !
        ! First the diagonal correction 
        !
        ! Use rho_minus as a workspace 
        rho_minus1%of_r(:,:) = 0.D0
        rho_minus1%of_g(:,:) = (0.D0, 0.D0)
        rho_minus1%of_r(:,1) = n_r(:)  ! orbital denisty rho_i in real space
        rho_minus1%of_r(:,2) = (3-2*spin_component)*n_r(:)  ! orbital magnetization m_i in real space
        rho_minus1%of_g(:,1) = n_g(:)  ! orbital denisty rho_i in reciprocal scape
        rho_minus1%of_g(:,2) = (3-2*spin_component)*n_g(:)  ! orbital magnetization m_i in reciprocal scape
        !
        etxc_minus1 = 0.D0; vtxc_minus1 = 0.D0; vxc_minus1(:,:) = 0.D0
        nspin_aux=nspin; nspin=2
        aux_r =0.D0 ; aux_g = (0.D0, 0.D0)
        CALL v_xc( rho_minus1, aux_r, aux_g, etxc_minus1, vtxc_minus1, vxc_minus1 )
        nspin=nspin_aux
        !
        WRITE(stdout , '(8x, "PZ corr const term, sh[n_i], Exc[n_i], int{v_xc[n_i] n_i}", 3F15.8)') &
            sh, etxc_minus1, vtxc_minus1
        !
        WRITE(stdout,'(8X, "Delta PZ", 2F15.8)')  -(sh+etxc_minus1), -(sh+etxc_minus1)*alpha_final(ibnd)
        !
        delta_eig(ibnd) = delta_eig(ibnd) -(sh+etxc_minus1)
        ham_right(ibnd,ibnd) = ham_right(ibnd,ibnd) -(sh+etxc_minus1) *alpha_final(ibnd)
        !
        !
        ! Eventually the off-diagonal contribution
        IF ( .NOT. on_site_only) THEN 
           !
           vpsi_r(:) = (0.D0, 0.D0)
           DO ir = 1, dffts%nnr
              vpsi_r (ir) = CMPLX( ( - v(ir,current_spin) - vxc_minus1(ir,current_spin) ),0.D0, kind=DP) * psic_1(ir)
           ENDDO 
           !
           ! 1) GO to G-space and store the ki gradient 
           CALL fwfft ('Wave', vpsi_r, dffts)
           vpsi_g(:) = vpsi_r(dffts%nl(igk_k(:,ik)))
           !
           DO k = 1, num_wann
              IF (k == ibnd ) CYCLE
              IF (ibnd .LE. num_wann_occ .AND. k .GT. num_wann_occ ) CYCLE  ! NO occ-empty matrix elements
              IF (ibnd .GT. num_wann_occ .AND. k .LE. num_wann_occ ) CYCLE  ! NO empty-occ matrix elements
              dhkipz = SUM ( CONJG(evc0(:,k)) * vpsi_g(:) * alpha_final(ibnd) )
              CALL mp_sum (dhkipz, intra_bgrp_comm)
              ham_right(k,ibnd) = ham_right(k,ibnd) + dhkipz
              !WRITE(*,*) k, ibnd, REAL(ham_right(k,ibnd)), AIMAG(ham_right(k,ibnd))
           ENDDO
           !
        ENDIF 
        !
     ENDIF
     !
     WRITE(stdout,'(8x, "orbital", i3, 3x, "spin", i3, 5x, "uKC_diag", F15.8 ," Ry", 3x, "rKC_diag", F15.8, " Ry", 3x, &
         &"alpha=", F15.8, 3x,/ )') ibnd, current_spin, delta_eig(ibnd), delta_eig(ibnd)*alpha_final(ibnd), &
         alpha_final(ibnd)
     !
  ENDDO orb_loop
  !
  ! h^R_{ij} = <\phi_i | H_j | \phi_j> 
  ! h^L_{ij} = <\phi_i | H_i | \phi_j> = [<\phi_j | H_i | \phi_i>]^* = [h^R_{ji}]^*
  IF (qp_symm) ham_right = 0.5D0*(ham_right + CONJG(TRANSPOSE(ham_right)))
  !
  ! Store the res in the global variable
  dH_wann(ik,:,:) = ham_right(:,:)
  !
  WRITE(stdout,'(5x,"INFO: KI calcualtion: Full Hamiltonian ... DONE",/ )')
  ! 
  DEALLOCATE (ham) 
  DEALLOCATE (ham_left) 
  DEALLOCATE (ham_right) 
  !
  DEALLOCATE (psic_1, vpsi_r, vpsi_g) 
  DEALLOCATE ( delta_vg, vh_rhog, delta_vg_)
  !
  CONTAINS
  !
  !----------------------------------------------------------------
  SUBROUTINE beyond_2nd (dH_wann, ddH, iwann)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : alpha_final, alpha_final, group_alpha, l_do_alpha
    USE control_lr,            ONLY : lrpa
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: dH_wann
    REAL(DP), INTENT(OUT) :: ddH
    !
    REAL(DP) second_der, delta, alpha_
    !
    REAL(DP), EXTERNAL :: get_clock
    !
    INTEGER, INTENT(IN) :: iwann
    !
    ddH = 0.D0
    !
    WRITE( stdout, '(/,8X, "INFO: Correction beyond 2nd order ...",/)')
    IF (lrpa) THEN
      !
      WRITE(*, '(8X,"INFO: l_alpha_corr and lrpa are NOT consistent.At RPA")')
      WRITE(*, '(8X,"      level there is no contribution beyond 2nd order.")')
      WRITE(*, '(8X,"      Nothing to do here. RETURN")')
      !
      RETURN
      !
    ENDIF
    !
    second_der = -REAL(dH_wann)
    !
    !DO iwann = 1, num_wann
      !
      delta =0.D0
      alpha_ = 0.0
      !
      ! Only if this is one of the unique orbitals
      IF ( l_do_alpha (iwann)) THEN
        !
        !
        ! ... Compute the difference between the parabolic extrapolation at N \pm 1 and the real
        ! ... value of the energy in the frozen orbital approximation ...
        CALL alpha_corr (iwann, delta)
        ddH = delta
        !dH_wann(iwann,iwann) = dH_wann(iwann,iwann)-ddH(iwann)
        !
        ! ... The new alpha that should be closer to the Finite-difference one ...
        ! ... Remember DeltaH is nothing but the second derivative wrt the orbital occupation ...
        alpha_ = (alpha_final(iwann)*second_der + delta)/ (second_der + delta)
        !
      ELSE
        !
        alpha_ = (alpha_final(group_alpha(iwann))*second_der + delta)/ (second_der + delta)
        !
      ENDIF
      !
      ! ... Write it just to compare with the FD one from CP ...
      IF (l_do_alpha (iwann)) THEN
       WRITE(stdout,'(8X, "INFO: iwann , LR-alpha, alpha", i3, 2f12.8)') &
               iwann, alpha_final(iwann), alpha_
      ELSE
       WRITE(stdout,'(8X, "INFO: iwann*, LR-alpha, alpha", i3, 2f12.8)') &
               iwann, alpha_final(iwann), alpha_
      ENDIF
      !
      ! Re-define the corrected screening parameter. 
      alpha_final(iwann) = alpha_ 
      WRITE( stdout, '(8X,"INFO: alpha RE-DEFINED ...", i5, f12.8)') iwann, alpha_final(iwann)
      WRITE(stdout, 900) get_clock('KCW')
      !
    !ENDDO
900 FORMAT('        total cpu time spent up to now is ',F10.1,' secs', / )
    !
  END SUBROUTINE beyond_2nd
  !
END subroutine dH_ki_full
