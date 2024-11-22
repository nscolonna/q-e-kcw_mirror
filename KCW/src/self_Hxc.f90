!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
SUBROUTINE self_hxc (iwann, shxc)
  !----------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_kcw,          ONLY : num_wann, nqstot, iurho_wann, nrho, spin_component
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
  INTEGER, INTENT(IN) :: iwann
  !
  ! Couters for the q point, wannier index. record length for the wannier density
  INTEGER :: iq, lrrho
  !
  ! The periodic part of the wannier orbital density
  COMPLEX(DP) :: rhowann(dffts%nnr, num_wann, nrho), rhor(dffts%nnr, nrho)
  COMPLEX(DP) :: delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag)
  !
  ! The self Hartree
  COMPLEX(DP) :: shxc
  !
  ! Auxiliary variables 
  COMPLEX(DP), ALLOCATABLE  :: rhog(:,:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
  !
  ! The weight of each q point
  REAL(DP) :: weight(nqstot)
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
    ! FIXME: Is the factor 0.5 correct? 
    IF (gamma_only) THEN 
      shxc = shxc + DBLE ( sum (CONJG(rhog (:,1)) * delta_vg(:,spin_component) ) ) * weight(iq) * omega
      IF (gstart == 2)  shxc = shxc - 0.5D0 * DBLE ( CONJG(rhog (1,1)) * delta_vg(1,spin_component) ) * weight(iq) * omega
    ELSE
      shxc = shxc + 0.5D0 * sum (CONJG(rhog (:,1)) * delta_vg(:,spin_component)) * weight(iq) * omega
    ENDIF
    !
    ! 
  ENDDO ! qpoints
  !
  DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
  !
  CALL mp_sum (shxc, intra_bgrp_comm)
 !
END SUBROUTINE self_hxc

