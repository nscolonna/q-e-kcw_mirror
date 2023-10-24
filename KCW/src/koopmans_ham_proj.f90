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
SUBROUTINE koopmans_ham_proj ()
  !---------------------------------------------------------------------
  ! Here the KI hamiltonian is written in terms of projectors on Wannnier 
  ! functions:
  ! \Delta_H_KI = \sum_n (1/2-P_n) * \Delta_n * |w_kn><w_kn| with 
  !   - P_n = \sum_kv f_kv <u_kv|w_kn><w_kn|u_kn> the occupation of the Wannier function n, 
  !   - \Delta_n=\alpha_n <W_n^2 | f_Hxc | W_n^2> the onsite KI correction 
  ! This correction is applied (for the moment) perturbatively on all the KS states
  ! available from the preceeding nscf calculation.
  !
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : nkstot, xk, ngk
  USE lsda_mod,              ONLY : nspin
  USE control_kcw,           ONLY : num_wann, Hamlt, nqstot, l_alpha_corr, evc0, &
                                    alpha_final, num_wann_occ, iuwfc_wann, spin_component
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx, npw, et, nbnd
  USE units_lr,              ONLY : lrwfc, iuwfc
  USE wavefunctions,         ONLY : evc
  USE buffers,               ONLY : get_buffer, save_buffer
  USE io_files,              ONLY : nwordwfc
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! The k point index 
  INTEGER :: ik, ibnd, ik_pw
  !
  ! The scalar part (independent on k) <rho_0i|v_0i|rho_0i>delta_ij, and occupations numbers of WFs
  COMPLEX(DP) :: deltah_scal (num_wann)
  REAL(DP) :: occ_mat(num_wann)
  !
  ! the KI hamiltonian, the KI contribution, and the new eigenvectors at a given k-point
  COMPLEX(DP) :: ham(num_wann,num_wann), deltaH(num_wann,num_wann), eigvc(npwx,num_wann)
  !
  COMPLEX(DP) :: delta_k, overlap
  !
  ! The new eigenalues 
  REAL(DP) :: et_ki(nbnd,nkstot)
  !
  ! The correction to the diagonal term beyond second order
  REAL(DP) :: ddH(num_wann)
  !
  INTEGER :: i, iwann, jwann
  ! 
  REAL(DP) :: ehomo, elumo
  REAL(DP) :: ehomo_ks, elumo_ks
  INTEGER  :: lrwannfc
  !
  WRITE(stdout, '(/,5X, "INFO: KI Hamiltonian using Projections")')
  !
  ! The scalar term R=0 i=j 
  deltah_scal=CMPLX(0.D0,0.D0,kind=DP)
  CALL ham_scalar (deltah_scal)
  CALL occupations(occ_mat)
  ! 
#ifdef DEBUG
  WRITE( stdout,'(/,5X," "Bare on-site correction:")')
  WRITE(stdout,'(5X,10(2F10.6, 2x))') (deltah_scal(iwann), iwann=1, num_wann)
#endif
  !
  ! ... The correction beyond 2nd order 
  IF (l_alpha_corr) CALL beyond_2nd (deltah_scal, ddH)
  !
  ! Apply alpha
  deltah_scal = alpha_final * deltah_scal
  !
#ifdef DEBUG
  WRITE( stdout,'(/,5X," Screened on-site correction:")')
  WRITE(stdout,'(5X,10(2F10.6, 2x))') (deltah_scal(iwann), iwann=1, num_wann)
#endif
  !
  ehomo=-1D+6
  elumo=+1D+6
  ehomo_ks=-1D+6
  elumo_ks=+1D+6
  !
  DO ik = 1, nkstot/nspin
    !
    ik_pw = ik + (spin_component-1)*(nkstot/nspin)
    !WRITE(*,*) ik_pw
    npw = ngk(ik_pw)
    lrwannfc = num_wann*npwx
    CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
    CALL get_buffer ( evc, nwordwfc, iuwfc, ik_pw )
    !
    ehomo_ks = MAX ( ehomo_ks, et(num_wann_occ  , ik_pw) )
    elumo_ks = MIN ( elumo_ks, et(num_wann_occ+1, ik_pw) )
    !
    WRITE( stdout, 9020 ) ( xk(i,ik_pw), i = 1, 3 )
    WRITE( stdout, '(10x, "KS  ",8F11.4)' ) (et(ibnd,ik_pw)*rytoev, ibnd=1,nbnd)
    !
    ! Correction at k: Dk = \sum_n [ (1/2-P_iw)D_ii <u_kv | u_kn><u_kn | uk_v>]
    DO ibnd = 1, nbnd
      delta_k = CMPLx(0.D0, 0.D0, kind=DP)
      DO iwann = 1, num_wann
        overlap = SUM(CONJG(evc(1:npw,ibnd))*(evc0(1:npw,iwann)))
        CALL mp_sum (overlap, intra_bgrp_comm)
        overlap = CONJG(overlap)*overlap
        delta_k = delta_k + (0.5D0 - occ_mat(iwann)) * deltah_scal(iwann) * overlap
        !WRITE(*,'(3X, 2I5, 2F20.12, F20.12, 2F20.12)') ibnd, iwann, deltah_scal(iwann), occ_mat(iwann), overlap 
      ENDDO
      et_ki(ibnd,ik) = et(ibnd,ik_pw) + REAL(delta_k)
      !WRITE (*,'(3I5, 2X, 4F20.12)') ik, ik_pw, ibnd, et(ibnd,ik_pw)*rytoev, et_ki(ibnd,ik)*rytoev, delta_k*rytoev
    ENDDO
    !
    ehomo = MAX ( ehomo, et_ki(num_wann_occ, ik ) )
    IF (num_wann > num_wann_occ) elumo = MIN ( elumo, et_ki(num_wann_occ+1, ik ) )
    !
    WRITE( stdout, '(10x, "KI  ",8F11.4)' ) (et_ki(ibnd,ik)*rytoev, ibnd=1,nbnd)
    !
  ENDDO
  !
  IF ( elumo < 1d+6) THEN
     WRITE( stdout, 9042 ) ehomo_ks*rytoev, elumo_ks*rytoev
     WRITE( stdout, 9044 ) ehomo*rytoev, elumo*rytoev
  ELSE
     WRITE( stdout, 9043 ) ehomo_ks*rytoev
     WRITE( stdout, 9045 ) ehomo*rytoev
  END IF
  !
  !Overwrite et
  DO ik = 1, nkstot/nspin
    ik_pw = ik + (spin_component-1)*(nkstot/nspin)
    et(1:nbnd, ik_pw) = et_ki(1:nbnd,ik)
  ENDDO
  ! 
9043 FORMAT(/,8x,'KS       highest occupied level (ev): ',F10.4 )
9042 FORMAT(/,8x, 'KS       highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9045 FORMAT(  8x, 'KI[2nd]  highest occupied level (ev): ',F10.4 )
9044 FORMAT(  8x, 'KI[2nd]  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
  !
  RETURN
  !
  CONTAINS 
  !
  ! !----------------------------------------------------------------
  SUBROUTINE occupations (occ_mat)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : num_wann, spin_component
    USE wvfct,                 ONLY : nbnd
    USE lsda_mod,              ONLY : nspin
    USE klist,                 ONLY : nkstot, ngk, wk
    USE wvfct,                 ONLY : wg
    USE mp_bands,              ONLY : intra_bgrp_comm
    USE mp,                    ONLY : mp_sum
    !
    REAL(DP), INTENT(INOUT) :: occ_mat(num_wann)
    INTEGER :: iwann, ik, ibnd, ik_pw
    COMPLEX(DP) :: f
    !
    ! The canonical occupation matrix (fermi dirac or alike)
    f = CMPLX(0.D0, 0.D0, kind=DP)
    occ_mat = CMPLX(0.D0,0.D0, kind=DP)
    !
    DO iwann= 1, num_wann
      !
      DO ik= 1, nkstot/nspin
        ik_pw = ik + (spin_component-1)*(nkstot/nspin)
        npw = ngk(ik_pw)
        lrwannfc = num_wann*npwx
        CALL get_buffer ( evc0, lrwannfc, iuwfc_wann, ik )
        CALL get_buffer ( evc, nwordwfc, iuwfc, ik_pw )
        !
        overlap = CMPLX(0.D0, 0.D0, kind=DP)
        DO ibnd = 1, nbnd
          overlap = SUM(CONJG(evc(1:npw,ibnd))*(evc0(1:npw,iwann)))
          CALL mp_sum (overlap, intra_bgrp_comm)
          overlap = CONJG(overlap)*overlap
          f=CMPLX(wg(ibnd,ik)/wk(ik), 0.D0, kind = DP)
          occ_mat(iwann) = occ_mat(iwann) + f * REAL(overlap)
        ENDDO
      ENDDO
    ENDDO
    !
    occ_mat = occ_mat/(nkstot/nspin)
    !
  END SUBROUTINE occupations

  !
  !----------------------------------------------------------------
  SUBROUTINE beyond_2nd (deltaH, ddH)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : num_wann, alpha_final, alpha_final_full
    USE control_lr,            ONLY : lrpa
    !
    IMPLICIT NONE  
    !
    COMPLEX(DP), INTENT(INOUT) :: deltaH(num_wann)
    REAL(DP), INTENT(OUT) :: ddH(num_wann)
    !
    REAL(DP) :: alpha_(num_wann), alpha_fd
    ! weight of the q points, alpha from LR, alpha after the all-order correction. 
    !
    REAL(DP) second_der(num_wann), delta 
    !
    INTEGER :: iwann
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
    DO iwann = 1, num_wann
      second_der(iwann) = -REAL(deltaH(iwann))
    ENDDO
    !
    !DO iwann = 1, num_wann_occ
    DO iwann = 1, num_wann
      delta =0.D0 
      alpha_(iwann) = alpha_final(iwann) 
      !
      !
       ! ... Compute the difference between the parabolic extrapolation at N \pm 1 and the real 
       ! ... value of the energy in the frozen orbital approximation ...
       CALL alpha_corr (iwann, delta)
       ddH(iwann) = delta
       !deltaH(iwann,iwann) = deltaH(iwann,iwann)-ddH(iwann)
       !
       ! ... The new alpha that should be closer to the Finite-difference one ...
       ! ... Remember DeltaH is nothing but the second derivative wrt the orbital occupation ...
       alpha_fd = (alpha_final(iwann)*second_der(iwann) + delta)/ (second_der(iwann)+delta)
       IF(nkstot/nspin == 1) alpha_final_full(iwann) = alpha_fd
       !
       ! ... Since the ham in the frozen approximation is approximated to second
       ! ... order only, this is the alpha we want to use. Only the
       ! ... numerator matters.
       alpha_(iwann) = (alpha_final(iwann)*second_der(iwann) + delta)/second_der(iwann)
       !
       ! ... Write it just to compare with the FD one from CP ... 
       WRITE(stdout,'(5X, "INFO: iwann, LR-alpha, FD-alpha, alpha", i3, 3f12.8)') iwann, alpha_final(iwann),alpha_fd,  alpha_(iwann)
       !
       !WRITE(stdout,'("Nicola", i3, 6f12.8)') iwann, deltaH(iwann,iwann)
       !
       ! Re-define the corrected screening parameter. 
       alpha_final(iwann) = alpha_(iwann) 
       WRITE( stdout, '(5X,"INFO: alpha RE-DEFINED ...", i5, f12.8)') iwann, alpha_final(iwann)
      !
    ENDDO
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
                                     spin_component
    USE fft_base,             ONLY : dffts
    USE cell_base,            ONLY : omega
    USE gvecs,                ONLY : ngms
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE mp,                   ONLY : mp_sum
    USE buffers,              ONLY : get_buffer
    !
    IMPLICIT NONE
    !
    ! The scalar contribution to the hamiltonian 
    COMPLEX(DP), INTENT(INOUT) :: deltah_scal (num_wann)
    !
    ! Couters for the q point, wannier index. record length for the wannier density
    INTEGER :: iq, iwann, lrrho
    !
    ! The periodic part of the wannier orbital density
    COMPLEX(DP) :: rhowann(dffts%nnr, num_wann), rhor(dffts%nnr), delta_vr(dffts%nnr,nspin), delta_vr_(dffts%nnr,nspin)
    !
    ! The self Hartree
    COMPLEX(DP) :: sh(num_wann)
    !
    ! Auxiliary variables 
    COMPLEX(DP), ALLOCATABLE  :: rhog(:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
    !
    ! The weight of each q point
    REAL(DP) :: weight(nqstot)
    !
    WRITE( stdout, '(/,5X, "INFO: KC SCALAR TERM CALCULATION ... START")')
    !
    ALLOCATE ( rhog (ngms) , delta_vg(ngms,nspin), vh_rhog(ngms), delta_vg_(ngms,nspin) )
    !
    DO iq = 1, nqstot
      !
      lrrho=num_wann*dffts%nnr
      CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
      !! Retrive the rho_wann_q(r) from buffer in REAL space
      !
      weight(iq) = 1.D0/nqstot ! No SYMM 
      !
      DO iwann = 1, num_wann  ! for each band, that is actually the perturbation
         !
         rhog(:)         = CMPLX(0.D0,0.D0,kind=DP)
         delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
         vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
         rhor(:)         = CMPLX(0.D0,0.D0,kind=DP)
         !
         rhor(:) = rhowann(:,iwann)
         !! The periodic part of the orbital desity in real space
         !
         CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
         !! The periodic part of the perturbation DeltaV_q(G)
         ! 
         sh(iwann) = sh(iwann) + 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:)                )*weight(iq)*omega
         deltah_scal(iwann) = deltah_scal(iwann) + sum (CONJG(rhog (:)) * delta_vg(:,spin_component)) &
                                     * weight(iq) * omega
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
END SUBROUTINE koopmans_ham_proj
