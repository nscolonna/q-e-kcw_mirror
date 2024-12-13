!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!-----------------------------------------------------------------------
SUBROUTINE dH_ki_quadratic (dH_wann, dH_wann_proj)
  !---------------------------------------------------------------------
  !
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE lsda_mod,              ONLY : nspin
  USE control_kcw,           ONLY : num_wann, nqstot, l_alpha_corr, qp_symm, &
                                    alpha_final, num_wann_occ, on_site_only, h_proj, nkstot_eff
  USE constants,             ONLY : rytoev
  USE buffers,               ONLY : get_buffer, save_buffer
  !
  IMPLICIT NONE
  !
  ! The k point index 
  INTEGER :: ik
  !
  ! The scalar part (independent on k) <rho_0i|v_0i|rho_0i>delta_ij
  COMPLEX(DP) :: deltah_scal (num_wann, num_wann)
  !
  ! The Hamiltonian due to the real contribution to the potential \int f_Hxc(r,r') rho_0i(r')
  COMPLEX(DP) deltah_real (num_wann, num_wann)
  COMPLEX, ALLOCATABLE :: ham_right(:,:)
#ifdef DEBUG
  COMPLEX(DP) ham (num_wann, num_wann)
#endif
  !
  COMPLEX(DP), INTENT (OUT) :: dH_wann(nkstot_eff,num_wann,num_wann)
  COMPLEX(DP), INTENT (OUT) :: dH_wann_proj(num_wann)
  !
  ! The correction to the diagonal term beyond second order
  REAL(DP) :: ddH(num_wann)
  !
  INTEGER :: iwann, jwann
  ! 
  REAL(DP), EXTERNAL :: get_clock
  !
  ALLOCATE (ham_right (num_wann, num_wann) )
  dH_wann = CMPLX(0.D0,0.D0, kind=DP) 
  !
  WRITE(stdout,'(/,5x,"INFO: qKI calculation: quadratic Hamiltonian ... ",/ )')
  IF (on_site_only) WRITE(stdout, '(/,5X, "INFO: Skipping off-diag: only R=0 and i=j")') 
  !
  ! The scalar term R=0 i=j 
  deltah_scal=CMPLX(0.D0,0.D0,kind=DP)
  CALL ham_scalar (deltah_scal)
  WRITE(stdout, 900) get_clock('KCW')
  !
  ! This is Just <w_n | f_Hxc | w_n> (needed for the the projection scheme a-la DFT+U)
  DO iwann = 1, num_wann 
    dH_wann_proj(iwann) = -2.D0*deltah_scal(iwann,iwann)
  ENDDO
  ! 
#ifdef DEBUG
  WRITE( stdout,'(/,5X," Scalar term Hamiltonian:")')
  DO iwann=1, num_wann
    WRITE(stdout,'(5X,10(2F10.6, 2x))') (deltah_scal(iwann,jwann), jwann=1,num_wann)
  ENDDO
  ! DEBUG
#endif
  !
  ! ... The correction beyond 2nd order 
  IF (l_alpha_corr) CALL beyond_2nd (deltah_scal, ddH)
  !
#ifdef DEBUG
  WRITE( stdout,'(/,5X," Scalar term Hamiltonian plus correction:")')
  ham = deltah_scal
  DO iwann=1, num_wann
    ham(iwann,iwann) = ham(iwann,iwann)-ddH(iwann)
    WRITE(stdout,'(5X,10(2F10.6, 2x))') (ham(iwann,jwann), jwann=1,num_wann)
  ENDDO
  ! DEBUG
#endif
  !
  ! Apply the screening coefficient
  dH_wann_proj(1:num_wann) = dH_wann_proj(1:num_wann)*alpha_final(1:num_wann)
  ! If formulation a-la DFT+U nothing else to do (only on-site term for now)
  IF (h_proj) GOTO 999
  !
  DO ik = 1, nkstot_eff
    !
    deltah_real = CMPLX(0.D0, 0.D0, kind = DP)
    !
    IF (.NOT. on_site_only) THEN 
       ! General routine for empty state hamiltonian at k
       CALL koopmans_ham_real_k ( ik, deltah_real )
       WRITE(stdout, 900) get_clock('KCW')
    ELSE
      ! SKIP off-diagonal elements in REAL SPACE (i.e. R/=0 i/=j) 
      ! If only R=0 and i=j, then the "real" contribution for empty states
      ! is <wi| {\int f_Hxc wi^2} |wi> (which is simply -2.D0 times
      ! the scalar controbution -0.5*<wi^2|f_Hxc|wi^2>) 
      DO iwann = num_wann_occ+1, num_wann
        deltah_real(iwann,iwann) = -2.D0*deltah_scal(iwann,iwann)
      ENDDO
    ENDIF
    !
    !
#ifdef DEBUG
    WRITE( stdout,'(/,5X," Real term Hamiltonian ik =:", i5)') ik
    DO iwann=1, num_wann
      WRITE(stdout,'(5X,10(2F10.6, 2x))') (deltah_real(iwann,jwann), jwann=1,num_wann)
    ENDDO
#endif
    !
    ! The KI contribution at k
    ham_right = deltah_scal + deltah_real
    !
    ! Apply the screening coefficient
!    DO iwann = 1, num_wann; dH_wann(:,iwann) = alpha_final (iwann)* dH_wann(:,iwann) ; ENDDO  
    DO iwann = 1, num_wann
      DO jwann = iwann , num_wann
        ham_right (iwann,jwann) = ham_right(iwann,jwann)*alpha_final(jwann)
      ENDDO
    ENDDO
    !
    ! h^R_{ij} = <\phi_i | H_j | \phi_j> 
    ! h^L_{ij} = <\phi_i | H_i | \phi_j> = [<\phi_j | H_i | \phi_i>]^* = [h^R_{ji}]^*
    IF (qp_symm) ham_right = 0.5D0*(ham_right + CONJG(TRANSPOSE(ham_right)))
    !
    ! Store the res in the global variable
    dH_wann(ik,:,:) = ham_right(:,:)
    !
#ifdef DEBUG
    WRITE(stdout, '("qKI Contribution to the Hamiltonian at k = ", i4)') ik
    DO iwann = 1, num_wann
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(dH_wann(ik, iwann,jwann)),AIMAG(dH_wann(ik, iwann,jwann)), jwann=1,num_wann)
    ENDDO
#endif
    !
  ENDDO
  !
900 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
  !
  DEALLOCATE ( ham_right ) 
999 CONTINUE
  RETURN
  !
  CONTAINS 
  !
  !----------------------------------------------------------------
  SUBROUTINE beyond_2nd (dH_wann, ddH)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : num_wann, alpha_final, group_alpha, l_do_alpha
    USE control_lr,            ONLY : lrpa
    !
    IMPLICIT NONE  
    !
    COMPLEX(DP), INTENT(INOUT) :: dH_wann(num_wann,num_wann)
    REAL(DP), INTENT(OUT) :: ddH(num_wann)
    !
    REAL(DP) :: alpha_
    ! alpha after the all-order correction. 
    !
    REAL(DP) second_der, delta 
    !
    INTEGER :: iwann
    !
    !REAL(DP), EXTERNAL :: get_clock
    !
    ddH = 0.D0
    !
    WRITE( stdout, '(/,5X, "INFO: Correction beyond 2nd order ...",/)')
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
    !DO iwann = 1, num_wann_occ
    DO iwann = 1, num_wann
      !
      delta  = 0.D0 
      alpha_ = 0.D0
      second_der = -REAL(dH_wann(iwann,iwann))
      !
      ! Only if this is one of the unique orbitals
      IF ( l_do_alpha (iwann)) THEN 
        !
        !
        ! ... Compute the difference between the parabolic extrapolation at N \pm 1 and the real 
        ! ... value of the energy in the frozen orbital approximation ...
        CALL alpha_corr (iwann, delta)
        ddH(iwann) = delta
        !
        ! ... Since the ham in the frozen approximation is approximated to second
        ! ... order only, this is the alpha we want to use. Only the
        ! ... numerator matters.
        alpha_ = (alpha_final(iwann)*second_der + delta)/second_der
      ELSE 
        !
        alpha_ = (alpha_final(group_alpha(iwann))*second_der + delta)/second_der
        !
      ENDIF
      !
      ! ... Write it just to compare with the FD one from CP ... 
      IF (l_do_alpha (iwann)) THEN
       WRITE(stdout,'(5X, "INFO: iwann , LR-alpha, alpha", i3, 3f12.8)') &
               iwann, alpha_final(iwann),  alpha_
      ELSE
       WRITE(stdout,'(5X, "INFO: iwann*, LR-alpha, alpha", i3, 3f12.8)') &
               iwann, alpha_final(iwann), alpha_
      ENDIF
      !
      !WRITE(stdout,'("Nicola", i3, 6f12.8)') iwann, dH_wann(iwann,iwann)
      !
      ! Re-define the corrected screening parameter. 
      alpha_final(iwann) = alpha_ 
      WRITE( stdout, '(5X,"INFO: alpha RE-DEFINED ...", i5, f12.8)') iwann, alpha_final(iwann)
      WRITE(stdout, 900) get_clock('KCW')
      !
    ENDDO
900 FORMAT('     total cpu time spent up to now is ',F10.1,' secs' )
    !
  END SUBROUTINE beyond_2nd
  !
  !----------------------------------------------------------------
  SUBROUTINE ham_scalar (deltah_scal)
    !----------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    USE control_kcw,          ONLY : num_wann, nqstot, iurho_wann, &
                                     spin_component, nrho
    USE fft_base,             ONLY : dffts
    USE cell_base,            ONLY : omega
    USE gvecs,                ONLY : ngms
    USE gvect,                ONLY : gstart
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE mp,                   ONLY : mp_sum
    USE buffers,              ONLY : get_buffer
    USE noncollin_module,     ONLY : nspin_mag
    USE control_flags,        ONLY : gamma_only
    !
    IMPLICIT NONE
    !
    ! The scalar contribution to the hamiltonian 
    COMPLEX(DP), INTENT(INOUT) :: deltah_scal (num_wann, num_wann)
    !
    ! Couters for the q point, wannier index. record length for the wannier density
    INTEGER :: iq, iwann, lrrho, ip
    !
    ! The periodic part of the wannier orbital density
    COMPLEX(DP) :: rhowann(dffts%nnr, num_wann, nrho), rhor(dffts%nnr, nrho)
    COMPLEX(DP) :: delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag)
    !
    ! The self Hartree
    COMPLEX(DP) :: sh(num_wann)
    !
    ! Auxiliary variables 
    COMPLEX(DP), ALLOCATABLE  :: rhog(:,:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
    !
    ! The weight of each q point
    REAL(DP) :: weight(nqstot)
    !
    WRITE( stdout, '(/,5X, "INFO: KC SCALAR TERM CALCULATION ... START")')
    !
    ALLOCATE ( rhog (ngms,nrho) , delta_vg(ngms,nspin_mag), vh_rhog(ngms), delta_vg_(ngms,nspin_mag) )
    !
    DO iq = 1, nqstot
      !
      lrrho=num_wann*dffts%nnr*nrho
      CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
      !! Retrive the rho_wann_q(r) from buffer in REAL space
      !
      weight(iq) = 1.D0/nqstot ! No SYMM 
      !
      DO iwann = 1, num_wann  ! for each band, that is actually the perturbation
         !
         rhog(:,:)       = CMPLX(0.D0,0.D0,kind=DP)
         delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
         vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
         rhor(:,:)       = CMPLX(0.D0,0.D0,kind=DP)
         !
         rhor(:,:) = rhowann(:,iwann,:)
         !! The periodic part of the orbital desity in real space
         !
         CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
         !! The periodic part of the perturbation DeltaV_q(G)
         ! 
         IF (gamma_only) THEN
           sh(iwann) = sh(iwann) + DBLE(sum (CONJG(rhog (:,1)) * vh_rhog(:))) * weight(iq)*omega
           IF (gstart == 2) sh(iwann) = sh(iwann) - 0.5D0*DBLE(CONJG(rhog (1,1)) * vh_rhog(1)) *weight(iq)*omega
         ELSE
           sh(iwann) = sh(iwann) + 0.5D0 * sum (CONJG(rhog (:,1)) * vh_rhog(:))*weight(iq)*omega
         ENDIF
         !
         IF (nspin_mag ==2 ) THEN
           !
           IF (gamma_only) THEN
              deltah_scal(iwann, iwann) = deltah_scal(iwann,iwann) &
                      - DBLE( sum (CONJG(rhog (:,1)) * delta_vg(:,spin_component))) * weight(iq) * omega
              IF (gstart == 2) deltah_scal(iwann, iwann) = deltah_scal(iwann,iwann) &
                      + 0.5D0*DBLE(CONJG(rhog (1,1)) * delta_vg(1,spin_component)) * weight(iq) * omega
           ELSE
              deltah_scal(iwann, iwann) = deltah_scal(iwann,iwann) - 0.5D0 * sum (CONJG(rhog (:,1)) * delta_vg(:,spin_component)) &
                                       * weight(iq) * omega
           ENDIF
           !
         ELSE 
            DO ip = 1, nspin_mag
              deltah_scal(iwann, iwann) = deltah_scal(iwann,iwann) - 0.5D0 * sum (CONJG(rhog (:,ip)) * delta_vg(:,ip)) &
                                          * weight(iq) * omega
            ENDDO
         ENDIF
         !
      ENDDO
      ! 
    ENDDO ! qpoints
    WRITE( stdout, '(/,5X, "INFO: KC SCALAR TERM CALCULATION ... END")')
    !
    DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
    !
    CALL mp_sum (deltah_scal, intra_bgrp_comm)
    CALL mp_sum (sh, intra_bgrp_comm)
   !
  END SUBROUTINE ham_scalar
  !
  !-----------------------------------------------------------------------
  SUBROUTINE koopmans_ham_real_k (ik, dH_wann)
    !---------------------------------------------------------------------
    !
    ! This routine compute the KI real term correction to second order to the KS
    ! Hamiltonian at a given k point ik given from input
    !
    USE kinds,                ONLY : DP
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE fft_wave,             ONLY : invfft_wave
    USE klist,                ONLY : igk_k, ngk
    USE mp,                   ONLY : mp_sum
    USE control_kcw,          ONLY : spin_component, num_wann, x_q, &
                                     num_wann_occ, evc0, iuwfc_wann, &
                                     map_ikq, shift_1bz, nrho, &
                                     map_ikq_minus, shift_1bz_minus
    USE buffers,              ONLY : get_buffer, save_buffer
    USE io_global,            ONLY : stdout
    USE control_kcw,          ONLY : iurho_wann
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE dv_of_drho_lr,        ONLY : dv_of_drho
    USE control_lr,           ONLY : lgamma
    USE lsda_mod,             ONLY : nspin
    USE gvecs,                ONLY : ngms
    USE solve_linter_koop_mod 
    USE qpoint,               ONLY : xq
    USE wvfct,                ONLY : npwx 
    USE cell_base,            ONLY : omega
    USE constants,            ONLY : rytoev
    USE noncollin_module,     ONLY : nspin_mag, npol
    USE control_flags,        ONLY: gamma_only
    USE gvect,                ONLY : gstart

    !
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(IN) :: ik
    ! the k-point index in hte orignal mesh of k-points
    !
    COMPLEX(DP), INTENT (INOUT) :: dH_wann(num_wann, num_wann)
    ! The KI real term contribution to the Hamiltonian
    !
    COMPLEX(DP) :: sh
    ! The self-Hartree of the Wannier function
    !
    INTEGER :: iq, nqs
    ! Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
    !
    INTEGER :: iwann, jwann, lrrho, lrwfc
    ! Band counters, leght of the rho record
    !
    COMPLEX(DP) :: rhowann(dffts%nnr, num_wann, nrho), rhor(dffts%nnr, nrho), delta_vr(dffts%nnr,nspin_mag), &
                   delta_vr_(dffts%nnr,nspin_mag)
    ! The periodic part of the wannier orbital density in r space
    ! The perturbig potential in real space
    ! The perturbig potential in real space (without the g=0 contribution) 
    !
    COMPLEX(DP), ALLOCATABLE  :: rhog(:,:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
    ! The periodic parte of the wannier density in G-space 
    ! The perturbig potential in G space
    ! The hartree potential due to the wannier in G-space 
    ! The perturbig potential in G space without the g=0 contribution 
    !
    REAL(DP) :: weight(nqstot)
    ! weight of the q points, alpha from LR, alpha after the all-order correction. 
    !
    INTEGER :: ikq, npw_k, npw_kq
    ! Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
    INTEGER :: ikq_m, npw_kq_m
    !
    COMPLEX(DP) ::  evc0_kq(npwx*npol, num_wann)
    ! Auxiliary vector to store the wfc at k+q
    !
    COMPLEX(DP) :: rho_r_nm(dffts%nnr, nrho), rho_g_nm(ngms, nrho), aux(dffts%nnr)
    REAL(DP) :: g_vect(3), g_vect_m(3)
    ! G vector that shift the k+q inside the 1BZ
    !
    COMPLEX(DP) :: evc_k_g (npwx*npol), evc_k_r (dffts%nnr,npol), phase(dffts%nnr)
    ! Auxiliary wfc in reciprocal and real space, the phase associated to the hift k+q-> k'
    !
    COMPLEX(DP) :: evc_kq_g (npwx*npol), evc_kq_r (dffts%nnr, npol)
    ! Auxiliary wfc in reciprocal and real space
    !
    LOGICAL :: off_diag = .TRUE., debug_nc = .true.
    ! compute Off-diagonal elements. NsC: not sure off_diag=.false. here makes sense: DO NOT CHANGE!!!!
    !
    INTEGER :: is
    !
    ! If there are no empty wannier functions, RETURN 
    IF (num_wann == num_wann_occ) RETURN 
    !
    IF (nspin_mag==2 .AND. debug_nc) &
      WRITE (stdout, '(/,5X,"INFO: debug_nc = ", L2, " Note: the k-q formula will be used.")') debug_nc
    !
    WRITE( stdout, '(/,/,5X,"INFO: qKI  HAMILTONIAN CALCULATION ik= ", i4, " ...", /)') ik
    !
    nqs = nqstot
    !
    !
    dH_wann  = CMPLX(0.D0,0.D0,kind=DP)
    sh       = CMPLX(0.D0,0.D0,kind=DP)
    rho_r_nm = CMPLX(0.D0,0.D0,kind=DP)
    !
    lrwfc = num_wann*npwx*npol
    CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik )
    ! Retrive the ks function at k (in the Wannier Gauge)
    ! IF (kcw_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: u_k(g) RETRIEVED"/)') 
    !
    CALL compute_map_ikq_single (ik, .true.)
    ! find tha map k+q --> k'+G and store the res 
    ! find also the map k-q --> k'+G and store the res
    !
    DO iq = 1, nqs
      !! Sum over the BZ 
      !
      xq(1:3)  = x_q(1:3,iq)
      !
      lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
      !
      lrrho=num_wann*dffts%nnr*nrho
      CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
      ! Retrive the rho_wann_q(r) from buffer in REAL space
      !IF (kcw_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: rhowan_q(r) RETRIEVED"/)') 
      !
      ALLOCATE ( rhog (ngms, nrho) , delta_vg(ngms,nspin_mag), vh_rhog(ngms), delta_vg_(ngms,nspin_mag) )
      !
      weight(iq) = 1.D0/nqs ! No SYMM 
      !
      !WRITE(stdout, '("weight =", i5, f12.8)') iq, weight(iq)
      !
      IF (nspin==4 .or. debug_nc ) THEN 
        ikq_m = map_ikq_minus(iq)
        g_vect_m(:) = shift_1bz_minus(:,iq)
        phase(:) = 0.D0
        CALL calculate_phase(g_vect_m, phase) 
        ! Calculate the phase associated to the k-q-> ikq map: exp[ -i(G_bar * r) ]
      ELSE
        ikq = map_ikq(iq)
        g_vect(:) = shift_1bz(:,iq)
        ! The index ikq in the 1BZ corresponding at the k+q, and the vector G_bar defining the shift 
        ! xk(:,ikq)+G_bar = xk(:,k+q)
        ! see compute_map_ikq_single.f90
        phase(:) = 0.D0
        CALL calculate_phase(g_vect, phase) 
        ! Calculate the phase associated to the k+q-> ikq map: exp[ -i(G_bar * r) ]
      ENDIF 
      !
      evc0_kq = CMPLX(0.D0,0.D0,kind=DP)
      lrwfc = num_wann * npwx * npol
      IF (nspin==4 .or. debug_nc ) THEN
        CALL get_buffer ( evc0_kq, lrwfc, iuwfc_wann, ikq_m )
        !Retrive the ks function at k-q (in the Wannier Gauge):
      ELSE
        CALL get_buffer ( evc0_kq, lrwfc, iuwfc_wann, ikq )
        ! Retrive the ks function at k+q (in the Wannier Gauge): 
      ENDIF 
      !
      !IF (kcw_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: u_kq(g) RETRIEVED")') 
      !
      DO iwann = num_wann_occ+1, num_wann
!      DO iwann = 1, num_wann
         !
         npw_k = ngk(ik)
         evc_k_g(:)   =  evc0(:,iwann)
         evc_k_r(:,:) = CMPLX(0.D0,0.D0,kind=DP)
         !
         IF (gamma_only) THEN
           ! NOTA: non collinear and Gamma_trick not compatible --> npol will always be 1 here
           evc_k_r(dffts%nl(igk_k(1:npw_k,ik)),1)  = evc_k_g(1:npw_k)
           evc_k_r(dffts%nlm(igk_k(1:npw_k,ik)),1)  = CONJG(evc_k_g(1:npw_k))
           CALL invfft ('Wave', evc_k_r(:,1), dffts)
         ELSE
           CALL invfft_wave (npwx, npw_k, igk_k (1,ik), evc_k_g , evc_k_r )
         ENDIF
         !! The wfc R=0 n=iwann in R-space at k
         !
         !DO jwann = iwann+1, num_wann 
         DO jwann = iwann, num_wann 
            !
            IF (.NOT. off_diag .AND. jwann /= iwann) CYCLE 
            !
            rhog(:,:)       = CMPLX(0.D0,0.D0,kind=DP)
            delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
            vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
            rhor(:,:)       = CMPLX(0.D0,0.D0,kind=DP)
            !
            rhor(:,:) = rhowann(:,jwann,:)
            ! The periodic part of the orbital desity R=0, n=iwann in real space
            !
            CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
            ! The periodic part of the perturbation DeltaV_q(G)
            !
            evc_kq_g = evc0_kq(:,jwann)
            evc_kq_r = CMPLX(0.D0,0.D0,kind=DP)
            !
            IF(nspin==4 .or. debug_nc ) THEN
              npw_kq_m = ngk(ikq_m)
              IF (gamma_only) THEN
                evc_kq_r(dffts%nl(igk_k(1:npw_kq_m,ik)),1)  = evc_kq_g(1:npw_kq_m)
                evc_kq_r(dffts%nlm(igk_k(1:npw_kq_m,ik)),1)  = CONJG(evc_kq_g(1:npw_kq_m))
                CALL invfft ('Wave', evc_kq_r(:,1), dffts)
              ELSE
                CALL invfft_wave (npwx, npw_kq_m, igk_k (1,ikq_m), evc_kq_g , evc_kq_r )
              ENDIF
            ELSE
              npw_kq = ngk(ikq)
              IF (gamma_only) THEN
                evc_kq_r(dffts%nl(igk_k(1:npw_kq,ik)),1)  = evc_kq_g(1:npw_kq)
                evc_kq_r(dffts%nlm(igk_k(1:npw_kq,ik)),1)  = CONJG(evc_kq_g(1:npw_kq))
                CALL invfft ('Wave', evc_kq_r(:,1), dffts)
              ELSE
                CALL invfft_wave (npwx, npw_kq, igk_k (1,ikq), evc_kq_g , evc_kq_r )
              ENDIF
            END IF
            !
            ! The wfc in R-space at k' <-- k+q where k' = (k+q)-G_bar
            ! evc_k+q(r) = sum_G exp[iG r] c_(k+q+G) = sum_G exp[iG r] c_k'+G_bar+G 
            !            = exp[-iG_bar r] sum_G' exp[iG'r] c_k'+G' = exp[-iG_bar r] *evc_k'(r)
            !
            IF (nspin==4) THEN
              !rho_r_nm(:,1) = rho_r_nm(:,1) + conjg(evc_k_r(:,2))*evc_kq_r(:,2)*phase(:)/nqs 
              !rho_r_nm(:,1) = CONJG(rho_r_nm(:,1))
              !rho_r_nm(:,2) = (conjg(evc_k_r(:,1))*evc_kq_r(:,2)+conjg(evc_k_r(:,2))*evc_kq_r(:,1))*phase(:)/nqs 
              !rho_r_nm(:,3) = CMPLX(0.D0,1.D0, kind=DP)*(conjg(evc_k_r(:,2))*evc_kq_r(:,1)-&
              !                                                          conjg(evc_k_r(:,1))*evc_kq_r(:,2))*phase(:)/nqs 
              !rho_r_nm(:,4) = (conjg(evc_k_r(:,1))*evc_kq_r(:,1)-conjg(evc_k_r(:,2))*evc_kq_r(:,2))*phase(:)/nqs
              !
              !phase = conjg(phase)
              ! \sum_{s1,s2} [u^{s1}_{k-q}(r)]^* \sigma_i^{s1,s2} [u^{}s2_{k}(r)] = \rho^i_{k-q,k}(r)
              ! i is the pauli matrices index i=0,x,y,z
              rho_r_nm(:,1) = ( conjg(evc_kq_r(:,1))*evc_k_r(:,1)*phase(:) + conjg(evc_kq_r(:,2))*evc_k_r(:,2)*phase(:) )/nqs 
              rho_r_nm(:,2) = ( conjg(evc_kq_r(:,1))*evc_k_r(:,2)*phase(:) + conjg(evc_kq_r(:,2))*evc_k_r(:,1)*phase(:) )/nqs 
              rho_r_nm(:,3) = (-CMPLX(0.D0,1.D0, kind=DP) * conjg(evc_kq_r(:,1))*evc_k_r(:,2)*phase(:) & 
                               +CMPLX(0.D0,1.D0, kind=DP) * conjg(evc_kq_r(:,2))*evc_k_r(:,1)*phase(:) )/nqs
              rho_r_nm(:,4) = ( conjg(evc_kq_r(:,1))*evc_k_r(:,1)*phase(:) - conjg(evc_kq_r(:,2))*evc_k_r(:,2)*phase(:) )/nqs
            ELSE IF (nspin==2 .and. debug_nc ) THEN
              !phase = conjg(phase)
              rho_r_nm(:,1) = conjg(evc_kq_r(:,1))*evc_k_r(:,1)*phase(:)/nqs 
              ! rho_{k-q,k}^{ji}
            ELSE IF (nspin==2 .and. .not. debug_nc) THEN
              rho_r_nm(:,1) = conjg(evc_k_r(:,1))*evc_kq_r(:,1)*phase(:)/nqs 
              ! rho_{k,k+q}^{ij}
            ENDIF
            !
            ! generalized density in r-space 
  !          IF (jwann == iwann) THEN
  !            WRITE(*,'("NICOLA evc i", 6F20.15)') evc_k_r(1:3)
  !            WRITE(*,'("NICOLA evc j", 6F20.15)') evc_kq_r(1:3)
  !            WRITE(*,'("NICOLA rho_ij", 6F20.15)') rho_r_nm(1:3)
  !            WRITE(*,'("NICOLA rho_ii", 6F20.15)') rhor(1:3)
  !          ENDIF
            !WRITE(*,'("NICOLA R", 2i5, 2F20.15)'), iwann, jwann, SUM( delta_vr(:,spin_component)*rho_r_nm(:) )/( dffts%nr1*dffts%nr2*dffts%nr3 )
            !
            DO is =1, nrho
              aux(:) = rho_r_nm(:,is)/omega
              CALL fwfft ('Rho', aux, dffts) 
              rho_g_nm(:,is) = aux(dffts%nl(:))
            ENDDO
            ! generalized density in G-spage
  !          IF (jwann == iwann) THEN
  !            WRITE(*,'("NICOLA rho_ij", 6F20.15)') rho_g_nm(1:3)
  !            WRITE(*,'("NICOLA rho_ii", 6F20.15)') rhog(1:3)
  !            WRITE(*,'("NICOLA vi_g", 6F20.15)') delta_vg(1:3,spin_component)
  !            WRITE(*,'("NICOLA i, q r ", i5, 10F16.8)') iwann, SUM (CONJG(rhog (:)) * delta_vg(:,spin_component))*weight(iq)*omega
  !            WRITE(*,'("NICOLA i, q r ", i5, 10F16.8)') iwann, SUM (CONJG(rho_g_nm (:)) * delta_vg(:,spin_component))*weight(iq)*omega
  !          ENDIF
            !
            IF (nspin==2 .and. .not. debug_nc) THEN
              !
              IF (gamma_only) THEN 
                 !
                 dH_wann(iwann, jwann) = dH_wann(iwann,jwann) + &
                    2.D0*DBLE(SUM((rho_g_nm(:,1))*CONJG(delta_vg(:,spin_component))))*weight(iq)*omega
                 IF (gstart == 2) dH_wann(iwann, jwann) = dH_wann(iwann,jwann) - &
                                     1.D0*DBLE((rho_g_nm(1,1))*CONJG(delta_vg(1,spin_component)))*weight(iq)*omega
                 !
              ELSE
                !
                dH_wann(iwann, jwann) = dH_wann(iwann,jwann) + SUM((rho_g_nm(:,1))*CONJG(delta_vg(:,spin_component)))& 
                        *weight(iq)*omega
                !
              ENDIF
              !
            ELSE IF (nspin==2 .and. debug_nc ) THEN 
              !
              IF (gamma_only) THEN
                 !
                 dH_wann(iwann, jwann) = dH_wann(iwann, jwann) + &
                    2.D0*DBLE(SUM(delta_vg(:,spin_component)*conjg(rho_g_nm(:,1))))*weight(iq)*omega
                 IF (gstart == 2) dH_wann(iwann, jwann) = dH_wann(iwann, jwann) - &
                                     1.D0*DBLE(delta_vg(1,spin_component)*conjg(rho_g_nm(1,1)))*weight(iq)*omega
                 !
              ELSE
                 !
                 dH_wann(iwann, jwann) = dH_wann(iwann, jwann) + SUM(delta_vg(:,spin_component)*conjg(rho_g_nm(:,1)))& 
                         *weight(iq)*omega
                 !
              ENDIF
              !
            ELSE
              DO is = 1, nspin_mag
                dH_wann(iwann, jwann) = dH_wann(iwann, jwann) + SUM(CONJG(rho_g_nm(:,is))*(delta_vg(:,is)))*weight(iq)*omega
              ENDDO
            ENDIF 
            !dH_wann(jwann, iwann) = dH_wann(jwann,iwann) + SUM(rho_g_nm(:)*CONJG(delta_vg(:,spin_component)))*weight(iq)*omega
            !WRITE(*,'("NICOLA G", 2i5, 2F20.15)') iwann, jwann, SUM (CONJG(rho_g_nm (:)) * delta_vg(:,spin_component))*weight(iq)*omega
            !WRITE(*,'("NICOLA G", 2i5, 2F20.15)') iwann, jwann, SUM (rho_g_nm (:) * CONJG(delta_vg(:,spin_component)))*weight(iq)*omega
            !
         ENDDO ! jwann
         ! 
         !
      ENDDO ! iwann
      !
      DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
      !
      !    
    ENDDO ! qpoints
    !
    dH_wann = nqstot*dH_wann
    CALL mp_sum (dH_wann, intra_bgrp_comm)
    CALL mp_sum (sh, intra_bgrp_comm)
    ! Sum over different processes (G vectors) 
    !
    RETURN 
    !
  END subroutine 
  !
END SUBROUTINE dH_ki_quadratic
