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
SUBROUTINE bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
  !---------------------------------------------------------------------
  !
  !! This routine compute the Hxc potential due to the (peridoc part of the) wannier
  !! charge density. V^{0n}_{Hxc}(r) = \int dr' f_Hxc(r,r') w^{0n}(r')
  !  
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvecs,                ONLY : ngms
  USE cell_base,            ONLY : tpiba2, omega
  USE control_kcw,          ONLY : spin_component, kcw_iverbosity, x_q, nrho
  USE gvect,                ONLY : g
  USE qpoint,               ONLY : xq
  USE constants,            ONLY : e2, fpi
  USE control_lr,           ONLY : lrpa
  USE martyna_tuckerman,    ONLY : wg_corr_h, do_comp_mt
  USE io_global,            ONLY : stdout
  USE coulomb,              ONLY : g2_convolution
  USE noncollin_module,     ONLY : domag, nspin_mag
  USE xc_lib,               ONLY : xclib_dft_is
  USE control_flags,      ONLY : gamma_only

  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT (OUT) :: rhog(ngms,nrho)
  ! ... periodic part of wannier density in G-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vg(ngms,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vg_(ngms,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space without g=0 contribution
  !
  COMPLEX(DP), INTENT (OUT) :: vh_rhog(ngms)
  ! ... Hartree perturbing potential [f_hxc(r,r') x wann(r')] in G-space
  ! 
  COMPLEX(DP), INTENT (OUT) :: delta_vr(dffts%nnr,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vr_(dffts%nnr,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space without g=0 contribution
  !
  COMPLEX(DP), INTENT (IN)  :: rhor(dffts%nnr,nrho)
  ! ... periodic part of wannier density in r-space
  !
  INTEGER, INTENT (IN)      :: iq
  ! ... q-point index
  !
  COMPLEX(DP)               :: aux(dffts%nnr), aux_(dffts%nnr)
  ! ... auxiliary vectors 
  !
  COMPLEX(DP), ALLOCATABLE  :: vaux(:)
  ! ... auxiliary vector
  !
  INTEGER                   :: ig, is, ip
  ! ... counters 
  !
  REAL(DP)                  :: qg2, eh_corr
  ! ... |q+G|^2, g=0 correction to the hartree energy 
  !
  LOGICAL                   :: lgamma
  !
  REAL(DP)                  :: xkq(3), xk(3)
  ! ... coordinate of k and k+q 
  !
  REAL(DP)                  :: fac(ngms)
  ! ... Coulomb kernel possibly with the special treatment of the q+g+0 component 
  !
  COMPLEX(DP)               :: vh_rhog_g0eq0(ngms)
  ! ... The hartree potential with th q+g=0 component set to zero
  !
  COMPLEX(DP), ALLOCATABLE :: rhor_(:,:)
  !
  !! The periodic part of the orbital density in g space  
  DO ip=1,nrho !<---CONSIDER TO SUBSTITUTE WITH nspin_mag
      aux(:) = rhor(:,ip)/omega
      CALL fwfft ('Rho', aux, dffts)  
      rhog(:,ip) = aux(dffts%nl(:))
  END DO
  !
  delta_vr = ZERO
  delta_vr_ = ZERO
  aux      = ZERO
  !
  ! .. First the xc contribution
  !
  IF (.NOT. lrpa) THEN 
    ALLOCATE ( rhor_(dffts%nnr,nspin_mag) )
    rhor_ = CMPLX(0.D0,0.D0,kind=DP)
    IF (nspin_mag == 4) THEN
      rhor_=rhor/omega
    ELSE
      rhor_(:,spin_component) = rhor(:,1)/omega
    ENDIF
    CALL dv_of_drho_xc(delta_vr, rhor_)
    DEALLOCATE (rhor_)
  ENDIF 
  !
  delta_vr_ = delta_vr 
  !
  xk(:) = 0.D0
  xkq(:) = -x_q(:,iq)
  ! ... auxiliary coordinate to pass to g2_convoltion 
  ! 
  ! ... The Hartree contribution first 
  !
  CALL g2_convolution(ngms, g, xk, xkq, fac)
  ! ... the hartree kernel (eventually within the 
  ! ... Gygi-Baldereschi sheme, see setup_coulomb) 
  !
  DO ig = 1, ngms
    !
    qg2 = SUM ( (g(:,ig)+x_q(:,iq))**2 )
    !
    vh_rhog_g0eq0(ig) =  e2 * fpi * rhog(ig,1) / (tpiba2 * qg2)
    IF (qg2 .lt. 1e-8) vh_rhog_g0eq0(ig) = (0.D0, 0.D0)
    ! ... set to zero the q+g=0 component
    !
    vh_rhog(ig) =  rhog(ig,1) * cmplx(fac(ig), 0.d0, KIND=DP)
    ! ... the hartree potential possibly with the special treatment of the q+g=0 component  
    !
  ENDDO
  !
  ! ... eventually add MT corrections
  !
  lgamma = (xq(1)==0.D0.AND.xq(2)==0.D0.AND.xq(3)==0.D0)
  IF (do_comp_mt .AND. lgamma) then
     ! ... this make sense only for a GAMMA only calculation 
     ! ... (do_compt_mt = .false if more than one k-point (see kcw_readin.f90) 
     !
     IF (kcw_iverbosity .gt. 1 ) WRITE(stdout,'(5x, " ADDING Martyna-Tuckerman correction" ,/)')
     !
     ALLOCATE( vaux( ngms ) )
     CALL wg_corr_h (omega, ngms, rhog, vaux, eh_corr)
     vh_rhog(1:ngms) = vh_rhog(1:ngms) +  vaux(1:ngms)
     vh_rhog_g0eq0(1:ngms) = vh_rhog_g0eq0(1:ngms) +  vaux(1:ngms)
     DEALLOCATE( vaux )
  ENDIF
  !
  ! ... Go to r-space 
  !
  aux=(0.d0,0.d0)
  aux_=(0.d0,0.d0)
  aux(dffts%nl(:))  = vh_rhog(:)                    ! G-space components of the Hartree potential
  aux_(dffts%nl(:)) = vh_rhog_g0eq0(:)
  !
  IF (gamma_only) THEN
    aux(dffts%nlm(:))  = CONJG(vh_rhog(:))
    aux_(dffts%nlm(:)) = CONJG(vh_rhog_g0eq0(:))
  ENDIF 
  !
  CALL invfft ('Rho', aux, dffts)
  CALL invfft ('Rho', aux_, dffts)
  !
  IF (nspin_mag==4 .and. domag) THEN    ! Perturbing potential due to Hartree
    delta_vr(:,1) = delta_vr(:,1)   + aux(:)
    delta_vr_(:,1) = delta_vr_(:,1) + aux_(:)
  ELSEIF (nspin_mag==4 .and. .NOT. domag) THEN
    delta_vr(:,1) = delta_vr(:,1) + aux(:)
    delta_vr_(:,1) = delta_vr_(:,1) + aux_(:)
  ELSE
    DO is = 1, nspin_mag
      delta_vr(:,is)  = delta_vr(:,is)   + aux(:)
      delta_vr_(:,is) = delta_vr_(:,is) + aux_(:)
    END DO
  END IF
  !
  DO is = 1, nspin_mag
    !
    aux(:) = delta_vr(:,is)
    aux_(:) = delta_vr_(:,is) 
    !
    CALL fwfft ('Rho', aux, dffts)
    CALL fwfft ('Rho', aux_, dffts)
    !
    delta_vg(:,is)  = aux(dffts%nl(:))
    delta_vg_(:,is) = aux_(dffts%nl(:))
    !
  ENDDO
  !
  !
END subroutine bare_pot

!-----------------------------------------------------------------------
subroutine dv_of_drho_xc (dv, drho)
  !-----------------------------------------------------------------------
  !
  !  This routine computes the change of the self consistent potential
  !  (Hartree and XC) due to the perturbation.
  !  Note: gamma_only is disregarded for PHonon calculations,
  !  TDDFPT purposes only.
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : e2, fpi
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft, invfft
  USE gvect,             ONLY : g
  USE noncollin_module,  ONLY : nspin_lsda, nspin_mag, nspin_gga
  USE funct,             ONLY : dft_is_nonlocc
  USE xc_lib,            ONLY : xclib_dft_is
  USE scf,               ONLY : rho, rho_core
  USE uspp,              ONLY : nlcc_any
  USE Coul_cut_2D_ph,    ONLY : cutoff_dv_of_drho
  USE qpoint,            ONLY : xq
  USE gc_lr,             ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE eqv,               ONLY : dmuxc

  IMPLICIT NONE
  COMPLEX(DP), INTENT(INOUT) :: dv(dfftp%nnr, nspin_mag)
  ! output: response XC potential
  COMPLEX(DP), INTENT(IN) :: drho(dfftp%nnr, nspin_mag)
  ! input:  response charge density

  INTEGER :: ir, is, is1
  ! counter on r vectors
  ! counter on spin polarizations
  REAL(DP) :: fac
  ! fac: 1 / nspin_lsda
  COMPLEX(DP), ALLOCATABLE :: drhotot(:, :)
  ! Total charge density. Sum of electronic (drho) and nlcc (drhoc) terms.
  !
  ALLOCATE(drhotot(dfftp%nnr, nspin_mag))
  !
  drhotot = drho
  !
  do is = 1, nspin_mag
     do is1 = 1, nspin_mag
        do ir = 1, dfftp%nnr
           dv(ir,is) = dv(ir,is) + dmuxc(ir,is,is1) * drhotot(ir,is1)
        enddo
     enddo
  enddo
  !
  ! Add gradient correction to the response XC potential.
  ! NB: If nlcc=.true. we need to add here its contribution.
  ! grho contains already the core charge
  !
  if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) + rho_core (:)
  !
  if ( xclib_dft_is('gradient') ) call dgradcorr(dfftp, rho%of_r, grho, dvxc_rr, &
                                dvxc_sr, dvxc_ss, dvxc_s, xq, drhotot, &
                                nspin_mag, nspin_gga, g, dv)
  !
  if ( dft_is_nonlocc() )  call dnonloccorr(rho%of_r, drhotot, xq, dv)
  !
  if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) - rho_core (:)
  !
end subroutine dv_of_drho_xc

