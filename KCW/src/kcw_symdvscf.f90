!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
SUBROUTINE kcw_symdvscf (dvtosym)
  !---------------------------------------------------------------------
  !
  ! Symmetrize the self-consistent change in the density due to the perturnbation
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE fft_base,         ONLY : dfftp
  USE cell_base,        ONLY : at
  USE symm_base,        ONLY : s, ft
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag
  USE lr_symm_base,     ONLY : minus_q, irotmq, nsymq, gi, gimq
  
  implicit none

  complex(DP) :: dvtosym (dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, nspin_mag)
  ! the potential to be symmetrized
  integer :: ftau(3,48), s_scaled(3,3,48)
  integer :: is, ri, rj, rk, i, j, k, isym, irot
  !  counters
  real(DP) :: gf(3), n(3)
  !  temp variables
  complex(DP), allocatable :: dvsym (:,:,:)
  ! the symmetrized potential
  complex(DP) ::  aux2, term(3,48), phase(48)
  ! auxiliary space
  ! the multiplication factor
  ! the phase factor
  !
  WRITE(*,*) "NICOLA", nsymq, minus_q
  if (nsymq == 1 .and. (.not.minus_q) ) return
  !
  call start_clock ('kcw_symdvscf')
  !
  !WRITE(*,*) "NICOLA drho(1:3) in =", dvtosym(1:3,1,1,1)
  allocate ( dvsym(dfftp%nr1x, dfftp%nr2x, dfftp%nr3x) )
  !
  n(1) = tpi / DBLE(dfftp%nr1)
  n(2) = tpi / DBLE(dfftp%nr2)
  n(3) = tpi / DBLE(dfftp%nr3)
  !
  CALL scale_sym_ops( nsymq, s, ft, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
       s_scaled, ftau )
  !
  ! Symmetrize with -q if present (Sq = -q + G)
  !
  IF (minus_q) THEN
     !
     ! Compute the phase factor exp(iG*r)
     ! where G = Sq + q
     !
     gf(:) =  gimq(1) * at(1,:) * n(:) + &
              gimq(2) * at(2,:) * n(:) + &
              gimq(3) * at(3,:) * n(:)
     term(:,1) = CMPLX(cos(gf(:)), sin(gf(:)), kind=DP)
     !
     do is = 1, nspin_lsda
        !
        phase(1) = (1.d0, 0.d0)
        !
        do k = 1, dfftp%nr3
           !
           do j = 1, dfftp%nr2
              !
              do i = 1, dfftp%nr1
                 !
                 ! Rotation and fractional translation: S^-1 * r - ftau
                 !
                 CALL rotate_grid_point( s_scaled(1,1,irotmq), ftau(1,irotmq),&
                      i,j,k,dfftp%nr1,dfftp%nr2,dfftp%nr3,ri,rj,rk)
                 !
                 aux2 = (0.d0, 0.d0)
                 aux2 = aux2 + dvtosym(ri,rj,rk,is) * phase(1) 
                 !
                 dvsym(i,j,k) = ( dvtosym(i,j,k,is) + CONJG(aux2) ) * 0.5d0
                 !
                 phase(1) = phase(1) * term(1,1)
                 !
              enddo
              !
              phase(1) = phase(1) * term(2,1)
              !
           enddo
           !
           phase(1) = phase(1) * term(3,1)
           !
        enddo
        !
        dvtosym(:,:,:,is) = dvsym(:,:,:)
        !
     enddo
     !
  ENDIF
  !
  ! Here we symmetrize with respect to the small group of q (Sq = q + G)
  !
  ! Compute the phase factor exp(iG*r)
  ! where G = Sq - q
  !
  DO isym = 1, nsymq
     gf(:) = gi(1,isym) * at(1,:) * n(:) + &
             gi(2,isym) * at(2,:) * n(:) + &
             gi(3,isym) * at(3,:) * n(:)
     term(:,isym) = CMPLX( cos(gf(:)), sin(gf(:)), kind=DP)
  ENDDO
  !
  DO is = 1, nspin_lsda
     !
     dvsym(:,:,:) = (0.d0, 0.d0)
     !
     do isym = 1, nsymq
        phase(isym) = (1.d0, 0.d0)
     enddo
     !
     do k = 1, dfftp%nr3
        !
        do j = 1, dfftp%nr2
           !
           do i = 1, dfftp%nr1
              !
              do isym = 1, nsymq
                 !
                 irot = isym
                 !
                 ! Rotation and fractional translation: S^-1 * r - ftau
                 !
                 CALL rotate_grid_point( s_scaled(1,1,irot), ftau(1,irot),&
                      i,j,k,dfftp%nr1,dfftp%nr2,dfftp%nr3,ri,rj,rk)
                 !
                 ! Calculate drho(S^-1 * r - ftau) * exp(iG*(r-tau_pert))
                 !
                 dvsym(i,j,k) = dvsym(i,j,k) + dvtosym(ri,rj,rk,is) * phase(isym)
                 !
              enddo
              !
              do isym = 1, nsymq
                 phase(isym) = phase(isym) * term(1,isym)
              enddo
              !
           enddo
           !
           do isym = 1, nsymq
              phase(isym) = phase(isym) * term(2,isym)
           enddo
           !
        enddo
        !
        do isym = 1, nsymq
           phase(isym) = phase(isym) * term(3,isym)
        enddo
        !
     enddo
     !
     ! Final normalization on the number of the symmetry operations
     ! in the small group of q.
     !
     dvtosym(:,:,:,is) = dvsym(:,:,:) / DBLE(nsymq)
     !
  ENDDO
  !
  !WRITE(*,*) "NICOLA drho(1:3) OUT =", dvtosym(1:3,1,1,1)
  deallocate (dvsym)
  !
  call stop_clock ('kcw_symdvscf')
  !
  RETURN
  !
END SUBROUTINE kcw_symdvscf
