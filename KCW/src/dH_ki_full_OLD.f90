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
  ! This routine compute the KI correction to the spectrum
  ! Nov 2024
  ! OBSOLETE: Refactorized in dH_ki_full.f90 
  !           Kept here for reference. Will disappear soon.
  !
  USE wavefunctions,         ONLY : psic
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE wvfct,                 ONLY : npw, nbnd, npwx
  USE fft_base,              ONLY : dffts, dfftp
  USE lsda_mod,              ONLY : lsda, current_spin,nspin
  USE klist,                 ONLY : nks, ngk, init_igk, igk_k, nkstot
  USE gvect,                 ONLY : ngm
  USE gvecs,                 ONLY : ngms
  USE buffers,               ONLY : get_buffer
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE control_kcw,           ONLY : kcw_at_ks, homo_only, alpha_final, &
                                    num_wann_occ, iuwfc_wann, kcw_iverbosity, &
                                    qp_symm, evc0, kipz_corr, num_wann, &
                                    spin_component, l_alpha_corr
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
  INTEGER :: ibnd, ir, k, dim_ham, nspin_aux, n_orb, lrwfc
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
  COMPLEX(DP) :: etmp2
  REAL(DP) :: etmp1
  !
  COMPLEX(DP) , ALLOCATABLE :: psic_1(:) , eigvc_ki(:,:)
  COMPLEX(DP) , ALLOCATABLE :: ham (:,:), ham_up(:,:), ham_dw(:,:), vpsi(:), vpsi_r(:), ham_aux(:,:), v_ki(:,:)
  REAL(DP), ALLOCATABLE :: eigvl_ki(:), et_aux(:,:)
  !
  LOGICAL :: off_diag = .TRUE.
  !
  COMPLEX(DP) :: delta_vr(dffts%nnr,nspin), delta_vr_(dffts%nnr,nspin)
  COMPLEX(DP), ALLOCATABLE  :: delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
  COMPLEX(DP) :: deltah_scal
  REAL(DP) :: ddH
  !
  !
  ALLOCATE (psic_1( dfftp%nnr), vpsi_r(dffts%nnr), vpsi(npwx), v_ki(npwx,nbnd))
  ALLOCATE (et_aux(nbnd,nks))
  ALLOCATE ( delta_vg(ngms,nspin), vh_rhog(ngms), delta_vg_(ngms,nspin) )
  !
  !
  WRITE(stdout,'(/,5x,"INFO: KC calcualtion: Full Hamiltonian ... ",/ )')
  !
  ! ... Loop over k_point: actually it's a loop over the spin ik=1 ==> spin_up; ik=2 ==> spin_dw ...
  !
  w1 = 1.D0 / omega
  !
  !alpha_final(:,:)=1.0  ! Just for debug
  !
  nspin_aux=nspin
  nspin=2
  CALL create_scf_type (rho_minus1)
  nspin=nspin_aux
  !
  lrwfc = num_wann*npwx
  CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik )
  ! Retrive the ks function at k 
  IF (kcw_iverbosity .gt. 1 ) WRITE(stdout,'(8X, "INFO: u_k(g) RETRIEVED"/)')
  !
  CALL compute_map_ikq_single (ik)
  ! find tha map k+q --> k'+G and store the res 
  !
  dim_ham = num_wann
  ALLOCATE ( ham (dim_ham,dim_ham) )
  ALLOCATE ( ham_up (dim_ham,dim_ham) )
  ALLOCATE ( ham_dw (dim_ham,dim_ham) )
  ham (:,:) = (0.D0,0.D0)
  ham_up (:,:) = (0.D0,0.D0)
  ham_dw (:,:) = (0.D0,0.D0)
  !
  !! ... KS Hamiltonian ....
  !ik_eff = ik + (spin_component -1)*nkstot/nspin
  !CALL ks_hamiltonian (evc0, ik_eff, dim_ham, .true.)
  !
  v_ki(:,:) = (0.D0,0.D0)
  !GOTO 101
  !
  !
  n_orb = num_wann
  orb_loop: DO ibnd = 1, n_orb
     !
     IF (kcw_at_ks .AND. homo_only .AND. (ibnd .ne. num_wann) ) CYCLE orb_loop ! only homo orbilal (for fast debug) 
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
     etmp=0.D0
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
     etmp1 = 0.D0
     delta_eig(ibnd) = 0.D0
     !
     IF ( ibnd .LE. num_wann_occ ) THEN
        !
        rho_minus1%of_r(:,1) = rho_minus1%of_r(:,1) - n_r(:)  ! denisty rho-rho_i in real space
        rho_minus1%of_g(:,1) = rho_minus1%of_g(:,1) - n_g(:)  ! denisty rho-rho_i in reciprocal scape
        ! The magnetization depends on which channel we remove the orbital from
        ! we reduce by n_i(r) if current_spin=1, we increase by n_i(r) if current_spin=2
        ! This is taken care by the factor (3-2*current_spin)
        rho_minus1%of_r(:,2) = rho_minus1%of_r(:,2) - (3-2*spin_component)*n_r(:)  ! magnetization m-m_i in real space
        rho_minus1%of_g(:,2) = rho_minus1%of_g(:,2) - (3-2*spin_component)*n_g(:)  ! magnetization m-m_i in reciprocal scape
        ! 
        etxc_minus1 = 0.D0; vtxc_minus1 = 0.D0; vxc_minus1(:,:) = 0.D0
        nspin_aux=nspin; nspin=2
        CALL v_xc( rho_minus1, rho_core, rhog_core, etxc_minus1, vtxc_minus1, vxc_minus1 )
        nspin=nspin_aux
        ! 
        !
        delta_eig(ibnd) = (-sh+(etxc-etxc_minus1-etmp)) 
        WRITE(stdout, '(8x, "KI corr const term, sh[n_i], Exc[n], Exc[n-n_i], int{v_xc[n] n_i} ", 4F14.8)') sh, &
            etxc, etxc_minus1, etmp 
        IF (lrpa) delta_eig(ibnd) = (-sh)  !! hartree only for debug
        !
        vpsi_r(:) = (0.D0, 0.D0)
        DO ir = 1, dffts%nnr
           vpsi_r (ir) = CMPLX( ( delta_eig(ibnd)),0.D0, kind=DP) * psic_1(ir)
        ENDDO
        !
        CALL fwfft ('Wave', vpsi_r, dffts)
        v_ki(:,ibnd) = vpsi_r(dffts%nl(igk_k(:,ik))) 
        !
!        WRITE(stdout,'("evc_occ", i5, 6F15.8)') ibnd, evc0(1:3,ibnd)
!        WRITE(stdout,'("vpsi", i5, 6F15.8)') ibnd, v_ki(1:3,ibnd)
        !
     ELSE
        !
        rho_minus1%of_r(:,1) = rho_minus1%of_r(:,1) + n_r(:)  ! denisty rho+rho_i in real space
        rho_minus1%of_g(:,1) = rho_minus1%of_g(:,1) + n_g(:)  ! denisty rho+rho_i in reciprocal scape
        ! The magnetization depends on which channel we remove the orbital from
        ! we reduce by n_i(r) if current_spin=1, we increase by n_i(r) if current_spin=2
        ! This is taken care by the factor (3-2*current_spin)
        rho_minus1%of_r(:,2) = rho_minus1%of_r(:,2) + (3-2*spin_component)*n_r(:)  ! magnetization m+rho_i in real space
        rho_minus1%of_g(:,2) = rho_minus1%of_g(:,2) + (3-2*spin_component)*n_g(:)  ! magnetization m+m_i in reciprocal scape
        ! 
        etxc_minus1 = 0.D0; vtxc_minus1 = 0.D0; vxc_minus1(:,:) = 0.D0
        nspin_aux=nspin; nspin=2
        CALL v_xc( rho_minus1, rho_core, rhog_core, etxc_minus1, vtxc_minus1, vxc_minus1 )
        nspin=nspin_aux
        ! 
        etmp1 = sum ( vxc_minus1(1:dfftp%nnr,current_spin) * n_r(1:dfftp%nnr) )
        etmp1= etmp1/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
        CALL mp_sum (etmp1, intra_bgrp_comm)
        !
        ! Scalar-term correction for Diagonal elements only
        delta_eig(ibnd) = -sh + (etxc_minus1 - etxc - etmp1)
        ! This is to compare with KCP implementation.
        WRITE(stdout, '(8x, "KC corr const term, sh[n_i], Exc[n], Exc[n+n_i], int{v_xc[n] n_i} ", 4F14.8)') sh, &
            etxc, etxc_minus1, etmp 
        IF (lrpa) delta_eig(ibnd) = (-sh)  !! hartree only for debug
!        delta_eig(ibnd) = -sh
!        WRITE(stdout,'("const_term_empty_hartree", i5, 3F15.8)') ibnd, -sh
!        WRITE(stdout,'("const_term_empty_xc", i5, 3F15.8)') ibnd, etxc_minus1, etxc, etmp1
!        WRITE(stdout,'("Vhartree", i5, 6F15.8)') ibnd, v(1:3,current_spin) 
        ! 
        ! Non-scalar component of the KI functional for both diagonal and off-diagonal elements
        ! This is the KI gradient (the non-scalar part)
        vpsi_r(:) = (0.D0, 0.D0)
        etmp2 = (0.D0, 0.D0)
        DO ir = 1, dffts%nnr
           vpsi_r (ir) = CMPLX( ( v(ir,current_spin) + vxc_minus1(ir,current_spin) - vxc(ir,current_spin) + &
                                  delta_eig(ibnd)),0.D0, kind=DP) * psic_1(ir)
           IF(lrpa) vpsi_r (ir) = CMPLX( ( v(ir,current_spin) + delta_eig(ibnd)),0.D0, kind=DP) * psic_1(ir)
           etmp2 = etmp2 + CONJG(psic_1(ir))*vpsi_r(ir)
        ENDDO 
!!        WRITE(*,'("NICOLA i, vi_r", i5, 10F16.8)') ibnd, v(1:5,spin_component)
        etmp2=etmp2/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
        CALL mp_sum (etmp2, intra_bgrp_comm)
        delta_eig(ibnd) = REAL(etmp2, kind=DP)
        !
        ! 1) GO to G-space and store the ki gradient 
        CALL fwfft ('Wave', vpsi_r, dffts)
        v_ki(:,ibnd) = vpsi_r(dffts%nl(igk_k(:,ik)))
        !
!        WRITE(stdout,'("evc_empty", i5, 6F15.8)') ibnd, evc0(1:3,ibnd)
!        WRITE(stdout,'("vpsi", i5, 6F15.8)') ibnd, v_ki(1:3,ibnd)
        !
     ENDIF
     !
     IF (kipz_corr) THEN !! Add the PZ part of the KIPZ potential (this is in the spirit of Perturbation theory) 
        ! 
        WRITE(stdout,'(5x, "ADDING KIPZ CORRECTION ....")') 
        !
        ! Use rho_minus as a workspace 
        rho_minus1%of_r(:,:) = 0.D0
        rho_minus1%of_g(:,:) = (0.D0, 0.D0) 
        rho_minus1%of_r(:,1) = n_r(:)  ! orbital denisty rho_i in real space
        rho_minus1%of_r(:,2) = n_r(:)  ! orbital denisty rho_i in real space
        rho_minus1%of_g(:,1) = n_g(:)  ! orbital denisty rho_i in reciprocal scape
        rho_minus1%of_g(:,2) = n_g(:)  ! orbital denisty rho_i in reciprocal scape
        !
        etxc_minus1 = 0.D0; vtxc_minus1 = 0.D0; vxc_minus1(:,:) = 0.D0
        nspin_aux=nspin; nspin=2
        aux_r =0.D0 ; aux_g = (0.D0, 0.D0)
        CALL v_xc( rho_minus1, aux_r, aux_g, etxc_minus1, vtxc_minus1, vxc_minus1 )
        nspin=nspin_aux
        !
        ! The constant term 
        etmp1 =0.D0
        etmp1 = sum ( vxc_minus1(1:dfftp%nnr,current_spin) * n_r(1:dfftp%nnr) )
        etmp1= etmp1/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
        CALL mp_sum (etmp1, intra_bgrp_comm)
        ! This is to compare with KCP implementation.
        WRITE(stdout , '(8x, "PZ corr const term, sh[n_i], Exc[n_i], int{v_xc[n_i] n_i}, int{v_xc[n_i] n_i}", 4F15.8)') &
            sh, etxc_minus1, etmp1, vtxc_minus1
        etmp1 = + sh - etxc_minus1 + etmp1
        !
        vpsi_r(:) = (0.D0, 0.D0)
        etmp2 = (0.D0, 0.D0)
        DO ir = 1, dffts%nnr
           vpsi_r (ir) = CMPLX( ( etmp1 - v(ir,current_spin) - vxc_minus1(ir,current_spin) ),0.D0, kind=DP) * psic_1(ir)
           etmp2 = etmp2 + CONJG(psic_1(ir))*vpsi_r(ir)
!           vpsi_r (ir) = CMPLX( ( v(ir,current_spin) + delta_eig(ibnd)),0.D0) * psic_1(ir)
        ENDDO
        ! The diagonal term. It's the shift of the KS eigenvalue in the canonical representation
        etmp2=etmp2/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
        CALL mp_sum (etmp2, intra_bgrp_comm)
        delta_eig(ibnd) = delta_eig(ibnd) + REAL(etmp2)
        WRITE(stdout,*) "8x, Delta KIPZ",  REAL(etmp2)
        !
        ! 1) GO to G-space and store the kipz gradient 
        CALL fwfft ('Wave', vpsi_r, dffts)
        v_ki(:,ibnd) = v_ki(:,ibnd) + vpsi_r(dffts%nl(igk_k(:,ik)))
        !
        WRITE(stdout,'(5x, ".... DONE ")') 
        !
     ENDIF
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
     v_ki(:,ibnd) = v_ki(:,ibnd) * alpha_final(ibnd)
     !
     WRITE(stdout,'(8x, "orbital", i3, 3x, "spin", i3, 5x, "uKC_diag", F15.8 ," Ry", 3x, "rKC_diag", F15.8, " Ry", 3x, &
         &"alpha=", F15.8, 3x,/ )') ibnd, current_spin, delta_eig(ibnd), delta_eig(ibnd)*alpha_final(ibnd), &
         alpha_final(ibnd)
     !
  ENDDO orb_loop
  !
!101  CONTINUE
  !
  ! ##### Build up the KI Hamiltonian 
  !
  ham_up(:,:)= (0.D0,0D0)
  ham_dw(:,:)= (0.D0,0D0)
  DO ibnd = 1, n_orb
     ! 
     DO k = ibnd, n_orb
        !
        etmp2 = 0.D0
        DO ir = 1, npw
!           etmp2 = etmp2 + CONJG(evc0(ir,ibnd)) * (v_ki(ir,k)+hpsi(ir,k))
           etmp2 = etmp2 + CONJG(evc0(ir,ibnd)) * (v_ki(ir,k))
        ENDDO
        CALL mp_sum (etmp2, intra_bgrp_comm)
!        WRITE(*,'("Nic: KI    matrix ele", 2i5, 3x, 2F20.15)') ibnd, k , etmp2
!        WRITE(*,'("Nic: KS    matrix ele", 2i5, 3x, 2F20.15)') ibnd, k , Hamlt(ik,ibnd,k)
!        WRITE(*,'("Nic: KS+KI matrix ele", 2i5, 3x, 2F20.15)') ibnd, k , Hamlt(ik,ibnd,k)+etmp2
        !
        ham_up(ibnd,k) = etmp2
        ! kill occ-empty matrix elements
        IF ( (ibnd .le. num_wann_occ)  .AND. (k .gt. num_wann_occ) ) ham_up (ibnd, k) = (0.D0,0D0)  
        ham_up(k,ibnd) = CONJG(ham_up(ibnd,k))
        !
     ENDDO
     !
     IF (qp_symm) THEN 
        !
        ! ... The lower part of the Hamiltonian ...
        ! ... In apost processing fashion the KI hamiltonian is not hermitian ...
        ! ... Here we build an Hermitian operator in the spirit of scQP GW ...
        DO k = 1, ibnd
           !
        etmp2 = 0.D0
           DO ir = 1, npw
!              etmp2 = etmp2 + CONJG(evc0(ir,ibnd)) * (v_ki(ir,k)+hpsi(ir,k))
              etmp2 = etmp2 + CONJG(evc0(ir,ibnd)) * (v_ki(ir,k))
           ENDDO
           CALL mp_sum (etmp2, intra_bgrp_comm)
           !
           ham_dw(ibnd,k) = etmp2
           IF (ibnd .gt. num_wann_occ .AND. k .le. num_wann_occ ) ham_dw (ibnd, k) = (0.D0,0D0)  ! kill occ-empty matrix elements
           ham_dw(k,ibnd) = CONJG(ham_dw(ibnd,k))
           !
        ENDDO
        !
     ENDIF
     !
     ! Shift of the KS eigenvalues (only diagonal KI)
!     delta_eig(ibnd) = DBLE( ham_up(ibnd,ibnd) )-et(ibnd,ik)
     delta_eig(ibnd) = DBLE( ham_up(ibnd,ibnd) )
!     WRITE(stdout,'("NICOLA", 3F15.8,/)')  DBLE( ham_up(ibnd,ibnd) ), DBLE(ham_dw(ibnd,ibnd))
     !
  ENDDO
  !
!!!!!! DEBUG
  IF ( .NOT. off_diag )THEN 
    ham_up = (0.D0,0.D0)
    ham_up = (0.D0,0.D0)
    DO ibnd = 1, n_orb; ham_up(ibnd,ibnd)= delta_eig(ibnd); enddo
  ENDIF
!!!!!!!!!!
  !
#ifdef DEBUG
  WRITE(stdout,'(5x, "###  KC HAMILTONIAN EMPTY ###")')
  DO k = num_wann_occ+1, num_wann_occ+4
     WRITE(stdout,'(5x, 16F20.15)') ham_up(num_wann_occ+1:num_wann_occ+4,k)
  ENDDO
  WRITE(stdout,*) 
#endif
  !
  !OPEN(987,FILE='hamiltonian_emp.dat',FORM='formatted',status='UNKNOWN')
  !DO k = num_wann_occ+1, num_wann
  !   WRITE(987,'(2E18.10)') ham_up(num_wann_occ+1:num_wann, k)
  !ENDDO
  !CLOSE(987)
  !
  ham=(0.D0,0.D0)
  !
  ! The KS hamiltonian in the Wannier Gauge (just to check)
  ham(:,:) = ham_up(:,:)
  IF (qp_symm)  ham(:,:) = 0.5D0*(ham_up(:,:)+ham_dw(:,:))
!  WRITE(stdout,'("NICOLA", 2F15.8)')  ham(2,3), ham(3,2) 
  !
  ! Store the res in the global variable
  dH_wann(ik,:,:) = ham(:,:)
  !
  WRITE(stdout,'(5x,"INFO: KC calcualtion: Full Hamiltonian ... DONE",/ )')
  ! 
  DEALLOCATE (ham) 
  DEALLOCATE (ham_up) 
  DEALLOCATE (ham_dw) 
  IF (ALLOCATED(eigvl_ki)) DEALLOCATE (eigvl_ki)
  IF (ALLOCATED(eigvc_ki)) DEALLOCATE (eigvc_ki)
  IF (ALLOCATED(ham_aux)) DEALLOCATE (ham_aux)
  !
  DEALLOCATE (psic_1, vpsi_r, vpsi, v_ki) 
  DEALLOCATE (et_aux)
  DEALLOCATE ( delta_vg, vh_rhog, delta_vg_)
  !
  CONTAINS
  !
  !----------------------------------------------------------------
  SUBROUTINE beyond_2nd (dH_wann, ddH, iwann)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : num_wann, alpha_final, alpha_final, group_alpha, l_do_alpha
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
