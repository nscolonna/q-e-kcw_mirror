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
SUBROUTINE dH_ki_full_qpoints (ik, dH_wann)
  !-----------------------------------------------------------------------
  !
  ! This routine compute the full KI or the pKIPZ correction to the spectrum
  ! KI potential:   See Eq. A8-A11 in PRB 90, 075135 (2014).
  ! KIPZ potential: See Eq. A15 in PRB 90, 075135 (2014).
  ! NB: Only the insulating case -i.e. fi=1 or fi=0 for each i- is considered
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
  ! <i | v^{KI}_j |j> = \int \phi*_i(r) [v_H[n_j](r) + v_xc[\rho+n_j](r) - v_xc[\rho](r) ] \phi_j(r)  if (f_i=f_j=0) (Eq. A10)
  !
  ! pKIPZ potenital adds to KI potential:
  !
  ! First the diagonal corrections are computed (i==j):
  ! <i | Dv^{PZ}_i |i> = -E_{Hxc}[n_i]  (Eq. A15)
  !
  ! Then possibly the off-diagonal (i /= j) 
  ! (The scalar term C does not contribute as <i|C|j> = 0 if i/=j):
  ! <i | Dv^{PZ}_j |j> = \int \phi*_i(r) [ -v_{Hxc}[n_j](r)] \phi_j(r)
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
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE wvfct,                 ONLY : nbnd, npwx
  USE fft_base,              ONLY : dffts
  USE lsda_mod,              ONLY : lsda, current_spin,nspin
  USE klist,                 ONLY : init_igk, nkstot
  USE gvecs,                 ONLY : ngms
  USE buffers,               ONLY : get_buffer
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE control_kcw,           ONLY : kcw_at_ks, homo_only, alpha_final, &
                                    num_wann_occ, &
                                    qp_symm, kipz_corr, num_wann, &
                                    spin_component, l_alpha_corr, on_site_only
  USE control_lr,            ONLY : lrpa
  USE mp,                    ONLY : mp_sum
  USE scf,                   ONLY : create_scf_type, scf_type
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
  INTEGER :: ibnd, dim_ham
  !
  REAL(DP) :: sh
  COMPLEX(DP) :: shxc 
  !
  REAL(DP) ::  etxc, etxc_minus1, etmp, delta_eig(nbnd)
  !
  COMPLEX(DP) , ALLOCATABLE :: ham_right(:,:), vpsi_r(:), vpsi_g(:)
  !
  COMPLEX(DP), ALLOCATABLE  :: vh_rhog(:), delta_vg_(:,:)
  REAL(DP) :: ddH
  INTEGER :: segno
  REAL(DP) :: krnl
  LOGICAL :: is_emp
  !
  ALLOCATE (vpsi_r(dffts%nnr), vpsi_g(npwx))
  ALLOCATE ( vh_rhog(ngms), delta_vg_(ngms,nspin) )
  !
  IF (kipz_corr) THEN 
     WRITE(stdout,'(/,5x,"INFO: pKIPZ calcualtion ... ",/ )')
  ELSE
     WRITE(stdout,'(/,5x,"INFO: KI calcualtion: Full Hamiltonian ... ",/ )')
  ENDIF
  IF (on_site_only) WRITE(stdout, '(/,5X, "INFO: Skipping off-diag: only R=0 and i=j")') 
  IF ( .NOT. on_site_only) THEN
     CALL infomsg('dH_ki_full_qpoints','WARNING: off-diagonal matrix element not implemented yet with k point sampling.')
  ENDIF 
  !
  dim_ham = num_wann
  ALLOCATE ( ham_right (dim_ham,dim_ham) )
  ham_right (:,:) = (0.D0,0.D0)
  !
  orb_loop: DO ibnd = 1, num_wann
     !
     ! keep track of occ ...
     segno = -1 
     ! ... vs empty states
     IF (ibnd .gt. num_wann_occ) segno = +1
     !
     IF (kcw_at_ks .AND. homo_only .AND. (ibnd .ne. num_wann_occ) ) CYCLE orb_loop ! only homo orbilal (for fast debug) 
     ! 
     IF ( lsda ) current_spin = spin_component
     ! 
     shxc = 0.0
     CALL self_hxc(ibnd, shxc)
     shxc=-shxc
     IF (l_alpha_corr) CALL beyond_2nd (shxc, ddH, ibnd)
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
     sh = 0.D0
     CALL self_hartree (ibnd, sh)
     krnl   = 0.d0
     is_emp = .false.
     IF (ibnd .gt. num_wann_occ) is_emp = .true.
     CALL xc_energy_n    ( ibnd, etxc, etmp, krnl ) 
     CALL xc_energy_npm1 ( ibnd, etxc_minus1, is_emp )
     !
     delta_eig(ibnd) = (sh - etxc + etxc_minus1 - segno * etmp)
     delta_eig(ibnd) = segno * delta_eig(ibnd)
     !
     IF (segno == -1) THEN
       WRITE(stdout, '(8x, "KI corr const term, sh[n_i], Exc[n], Exc[n-n_i], int{v_xc[n] n_i} ", 4F14.8)') sh, &
           etxc, etxc_minus1, etmp 
     ELSE
       WRITE(stdout, '(8x, "KI corr const term, sh[n_i], Exc[n], Exc[n+n_i], int{v_xc[n] n_i} ", 4F14.8)') sh, &
           etxc, etxc_minus1, etmp 
     ENDIF
     IF (lrpa) delta_eig(ibnd) = segno * sh  !! hartree only for debug
     WRITE(stdout,'(8X, "Delta KI", 2F15.8)')  delta_eig(ibnd), delta_eig(ibnd)*alpha_final(ibnd)
     !
     ham_right(ibnd,ibnd) = delta_eig(ibnd) * alpha_final(ibnd)
     !WRITE(stdout,*) ibnd, ibnd, REAL(ham_right(ibnd,ibnd)), AIMAG(ham_right(ibnd,ibnd))
     !
     ! pKIPZ
     IF (kipz_corr) THEN 
        !
        ! First the diagonal correction 
        !
        ! The xc energy of rho_i
        CALL xc_energy_iwann ( ibnd, etxc_minus1) 
        !
        WRITE(stdout , '(8x, "PZ corr const term, sh[n_i], Exc[n_i]}", 3F15.8)') &
            sh, etxc_minus1
        !
        WRITE(stdout,'(8X, "Delta PZ", 2F15.8)')  -(sh+etxc_minus1), -(sh+etxc_minus1)*alpha_final(ibnd)
        !
        delta_eig(ibnd) = delta_eig(ibnd) -(sh+etxc_minus1)
        ham_right(ibnd,ibnd) = ham_right(ibnd,ibnd) -(sh+etxc_minus1) *alpha_final(ibnd)
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
  IF (kipz_corr) THEN 
          WRITE(stdout,'(5x,"INFO: pKIPZ calcualtion ... DONE",/ )')
  ELSE
          WRITE(stdout,'(5x,"INFO: KI calcualtion: Full Hamiltonian ... DONE",/ )')
  ENDIF
  ! 
  DEALLOCATE (ham_right) 
  !
  DEALLOCATE (vpsi_r, vpsi_g) 
  DEALLOCATE ( vh_rhog, delta_vg_)
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
END subroutine dH_ki_full_qpoints
