!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gcx( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud, &
                   run_on_gpu_ )
  !-------------------------------------------------------------------------
  !! Wrapper to gpu or non gpu version of \(\texttt{xc_gcx}\).
  !
  USE kind_l,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: ns
  !! spin dimension for input
  REAL(DP), INTENT(IN) :: rho(:,:)
  !! Charge density
  REAL(DP), INTENT(IN) :: grho(:,:,:)
  !! gradient
  REAL(DP), INTENT(OUT) :: ex(:)
  !! exchange energy
  REAL(DP), INTENT(OUT) :: ec(:)
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1x(:,:)
  !! exchange potential (density part)
  REAL(DP), INTENT(OUT) :: v2x(:,:)
  !! exchange potential (gradient part)
  REAL(DP), INTENT(OUT) :: v1c(:,:)
  !! correlation potential (density part)
  REAL(DP), INTENT(OUT) :: v2c(:,:)
  !! correlation potential (gradient part)
  REAL(DP), INTENT(OUT), OPTIONAL :: v2c_ud(:)
  !! correlation potential, cross term
  LOGICAL, INTENT(IN), OPTIONAL :: run_on_gpu_
  !! whether you wish to run on gpu in case use_gpu is true
  !
  LOGICAL :: run_on_gpu
  REAL(DP), ALLOCATABLE :: v2c_dummy(:)
  !
  run_on_gpu = .FALSE.
  !
  IF ( PRESENT(run_on_gpu_) ) run_on_gpu = run_on_gpu_
  !
  IF (ns==2 .AND. .NOT. PRESENT(v2c_ud)) CALL xclib_infomsg( 'xc_gcx', 'WARNING: cross &
                &term v2c_ud not found xc_gcx (gga) call with polarized case' )
  !
  IF ( run_on_gpu ) THEN
    !
    !$acc data deviceptr( rho(length,ns), grho(3,length,ns), ex(length), ec(length),      &
    !$acc&                v1x(length,ns), v2x(length,ns), v1c(length,ns), v2c(length,ns), &
    !$acc&                v2c_ud(length) )
    IF (PRESENT(v2c_ud)) THEN
      CALL xc_gcx_gpu( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
    ELSE
      ALLOCATE( v2c_dummy(length) )
      CALL xc_gcx_gpu( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_dummy )
      DEALLOCATE( v2c_dummy )
    ENDIF
    !$acc end data
    !
  ELSE
    !
    IF (PRESENT(v2c_ud)) THEN
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
    ELSE
      ALLOCATE( v2c_dummy(length) )
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_dummy )
      DEALLOCATE( v2c_dummy )
    ENDIF
    !
  ENDIF  
  !
  RETURN
  !
END SUBROUTINE
!
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
  !-------------------------------------------------------------------------
  !! Wrapper routine. Calls xc_gga-driver internal routines or the external
  !! ones from libxc, depending on the input choice.
  !
  !! NOTE: differently from 'xc_lda_drivers', here the input rho is in (up,down)
  !!       form (in the LSDA case).
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_setting_params,   ONLY: xc_func, xc_info, libxc_flags
#endif
  !
  USE kind_l,               ONLY: DP
  USE dft_setting_params,   ONLY: igcx, igcc, is_libxc, rho_threshold_gga, &
                                  grho_threshold_gga, rho_threshold_lda
  USE qe_drivers_gga
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: ns
  !! spin dimension for input
  REAL(DP), INTENT(IN) :: rho(length,ns)
  !! Charge density
  REAL(DP), INTENT(IN) :: grho(3,length,ns)
  !! gradient
  REAL(DP), INTENT(OUT) :: ex(length)
  !! exchange energy
  REAL(DP), INTENT(OUT) :: ec(length)
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1x(length,ns)
  !! exchange potential (density part)
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  !! exchange potential (gradient part)
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  !! correlation potential (density part)
  REAL(DP), INTENT(OUT) :: v2c(length,ns)
  !! correlation potential (gradient part)
  REAL(DP), INTENT(OUT) :: v2c_ud(length)
  !! correlation potential, cross term
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  REAL(DP), ALLOCATABLE :: rho_lxc(:), sigma(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_rho(:), vx_sigma(:)
  REAL(DP), ALLOCATABLE :: vc_rho(:), vc_sigma(:)
  !
  INTEGER :: fkind_x, np
  REAL(DP), PARAMETER :: pi34 = 0.6203504908994_DP
  !
  LOGICAL :: POLARIZED
  INTEGER :: pol_unpol
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  REAL(DP), ALLOCATABLE :: rh(:), zeta(:)
  REAL(DP), ALLOCATABLE :: grho2(:,:), grho_ud(:)
  !
  INTEGER :: k, is
  REAL(DP) :: sgn(2)
  REAL(DP) :: rho_up, rho_dw, grho_up, grho_dw
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  ex = 0.0_DP ;  v1x = 0.0_DP ;  v2x = 0.0_DP
  ec = 0.0_DP ;  v1c = 0.0_DP ;  v2c = 0.0_DP
  v2c_ud = 0.0_DP
  !
#if defined(__LIBXC)
  !
  fkind_x = -1
  lengthxc = length
  !
  POLARIZED = .FALSE.
  IF (ns == 2) THEN
     POLARIZED = .TRUE.
  ENDIF
  !
  pol_unpol = 1
  np = 1
  IF ( ns == 2 ) THEN
     pol_unpol = 2
     np = 3
  ENDIF
  !
  ALLOCATE( rho_lxc(length*ns) )
  ALLOCATE( sigma(length*np) )
  !
  ALLOCATE( ex_lxc(length)    , ec_lxc(length)      )
  ALLOCATE( vx_rho(length*ns) , vx_sigma(length*np) )
  ALLOCATE( vc_rho(length*ns) , vc_sigma(length*np) )
  !
  !
  IF ( ns == 1 ) THEN
    !
    DO k = 1, length
      rho_lxc(k) = ABS( rho(k,1) )
      IF ( rho_lxc(k) > rho_threshold_gga ) &
         sigma(k) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
    ENDDO
    !
  ELSE
    !
    DO k = 1, length
       rho_lxc(2*k-1) = rho(k,1)
       rho_lxc(2*k)   = rho(k,2)
       !
       sigma(3*k-2) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
       sigma(3*k-1) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                      grho(3,k,1) * grho(3,k,2)
       sigma(3*k)   = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
    ENDDO
    !
  ENDIF
  !
  IF ( ns==1 .AND. ANY(.NOT.is_libxc(3:4)) ) THEN
     !
     !$acc data copyin( rho_lxc, sigma ), copyout( ex, ec, v1x, v2x, v1c, v2c )
     !$acc host_data use_device( rho_lxc, sigma, ex, ec, v1x, v2x, v1c, v2c )
     CALL gcxc( length, rho_lxc, sigma, ex, ec, v1x(:,1), v2x(:,1), &
                v1c(:,1), v2c(:,1) )
     !$acc end host_data
     !$acc end data    
     !
     DO k = 1, length
        sgn(1) = SIGN(1._DP, rho(k,1))
        ex(k) = ex(k) * sgn(1)
        ec(k) = ec(k) * sgn(1)
     ENDDO
     !
  ENDIF
  !
  ! ---- GGA CORRELATION
  !
  IF ( is_libxc(4) ) THEN  !lda part of LYP not present in libxc (still so? - check)
    !
    fkind_x = xc_f03_func_info_get_kind( xc_info(4) )
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(4), small )!rho_threshold_gga )
    IF (libxc_flags(4,0)==1) THEN
      CALL xc_f03_gga_exc_vxc( xc_func(4), lengthxc, rho_lxc(1), sigma(1), ec_lxc(1), vc_rho(1), vc_sigma(1) )
    ELSE
      CALL xc_f03_gga_vxc( xc_func(4), lengthxc, rho_lxc(1), sigma(1), vc_rho(1), vc_sigma(1) )
      ec_lxc = 0.d0
    ENDIF
    !
    IF (.NOT. POLARIZED) THEN
      DO k = 1, length
        IF ( rho_lxc(k) <= rho_threshold_lda ) CYCLE
        ec(k) = ec_lxc(k) * rho_lxc(k) * SIGN(1.0_DP, rho(k,1))
        v1c(k,1) = vc_rho(k)
        IF ( rho_lxc(k) <= rho_threshold_gga .OR. &
             SQRT(ABS(sigma(k))) <= grho_threshold_gga) CYCLE
        v2c(k,1) = vc_sigma(k)*2.d0
      ENDDO
    ELSE
      DO k = 1, length
        rho_up = rho_lxc(2*k-1)
        rho_dw = rho_lxc(2*k)
        grho_up = SQRT(ABS(sigma(3*k-2)))
        grho_dw = SQRT(ABS(sigma(3*k)))
        IF ( rho_up <= rho_threshold_lda .OR. rho_dw <= rho_threshold_lda ) CYCLE
        ec(k) = ec_lxc(k) * (rho_up+rho_dw)
        v1c(k,1) = vc_rho(2*k-1)
        v1c(k,2) = vc_rho(2*k)
        IF ( rho_up <= rho_threshold_gga .OR. rho_dw <= rho_threshold_gga .OR. &
             grho_up<=grho_threshold_gga .OR. grho_dw<=grho_threshold_gga ) CYCLE
        v2c(k,1) = vc_sigma(3*k-2)*2.d0
        v2c_ud(k)= vc_sigma(3*k-1)
        v2c(k,2) = vc_sigma(3*k)*2.d0
      ENDDO
    ENDIF
    !  
  ELSEIF ( (.NOT.is_libxc(4)) .AND. fkind_x/=XC_EXCHANGE_CORRELATION ) THEN
    !
    IF ( ns /= 1 ) THEN
       !
       ALLOCATE( grho2(length,ns) )
       !
       DO is = 1, 2
          grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
       ENDDO
       !
       IF (igcc==3 .OR. igcc==7 .OR. igcc==13 ) THEN
          !
          ALLOCATE( grho_ud(length) )
          !
          grho_ud = grho(1,:,1) * grho(1,:,2) + grho(2,:,1) * grho(2,:,2) + &
                    grho(3,:,1) * grho(3,:,2)
          !
          !$acc data copyin( rho, grho2, grho_ud ), copyout( ec, v1c, v2c, v2c_ud )
          !$acc host_data use_device( rho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
          CALL gcc_spin_more( length, rho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
          !$acc end host_data
          !$acc end data
          !
          DEALLOCATE( grho_ud )
          !
       ELSE
          !
          ALLOCATE( rh(length), zeta(length) )
          !
          rh = rho(:,1) + rho(:,2)
          !
          zeta = 2.0_DP ! trash value, gcc-routines get rid of it when present
          WHERE ( rh > rho_threshold_gga ) zeta = ( rho(:,1) - rho(:,2) ) / rh(:)
          !
          grho2(:,1) = ( grho(1,:,1) + grho(1,:,2) )**2 + &
                       ( grho(2,:,1) + grho(2,:,2) )**2 + &
                       ( grho(3,:,1) + grho(3,:,2) )**2
          !
          IF ( igcc/=14 ) THEN
            !$acc data copyin( rh, zeta, grho2 ), copyout( ec, v1c, v2c )
            !$acc host_data use_device( rh, zeta, grho2, ec, v1c, v2c )          
            CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
            !$acc end host_data
            !$acc end data            
          ENDIF  
          !
          v2c(:,2) = v2c(:,1)
          IF ( ns==2 ) v2c_ud(:) = v2c(:,1)
          !
          DEALLOCATE( rh, zeta )
          !
       ENDIF
       !   
       DEALLOCATE( grho2 )
       !
    ENDIF
    !
  ENDIF  
  !
  ! --- GGA EXCHANGE
  !
  IF ( is_libxc(3) ) THEN
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(3), grho_threshold_gga )
    IF (libxc_flags(3,0)==1) THEN
      CALL xc_f03_gga_exc_vxc( xc_func(3), lengthxc, rho_lxc(1), sigma(1), ex_lxc(1), vx_rho(1), vx_sigma(1) )
    ELSE
      CALL xc_f03_gga_vxc( xc_func(3), lengthxc, rho_lxc(1), sigma(1), vx_rho(1), vx_sigma(1) )
      ex_lxc = 0.d0
    ENDIF
    !
    IF (.NOT. POLARIZED) THEN
      DO k = 1, length
        IF ( rho_lxc(k) <= rho_threshold_lda ) CYCLE
        ex(k) = ex_lxc(k) * rho_lxc(k) * SIGN(1.0_DP, rho(k,1))
        v1x(k,1) = vx_rho(k)
        IF ( rho_lxc(k) <= rho_threshold_gga .OR. &
             SQRT(ABS(sigma(k))) <= grho_threshold_gga) CYCLE
        v2x(k,1) = vx_sigma(k)*2.d0
      ENDDO
    ELSE
      DO k = 1, length
        rho_up = rho_lxc(2*k-1)
        rho_dw = rho_lxc(2*k)
        grho_up = SQRT(ABS(sigma(3*k-2)))
        grho_dw = SQRT(ABS(sigma(3*k)))
        IF ( rho_up <= rho_threshold_lda .OR. rho_dw <= rho_threshold_lda ) CYCLE
        ex(k) = ex_lxc(k) * (rho_up+rho_dw)
        v1x(k,1) = vx_rho(2*k-1)
        v1x(k,2) = vx_rho(2*k)
        IF ( rho_up <= rho_threshold_gga .OR. rho_dw <= rho_threshold_gga .OR. &
             grho_up<=grho_threshold_gga .OR. grho_dw<=grho_threshold_gga ) CYCLE
        v2x(k,1) = vx_sigma(3*k-2)*2.d0
        v2x(k,2) = vx_sigma(3*k)*2.d0
      ENDDO
    ENDIF
    !
  ELSE
    !
    IF ( ns /= 1 ) THEN
       !
       ALLOCATE( grho2(length,ns) )
       !
       DO is = 1, 2
          grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
       ENDDO
       !
       !$acc data copyin( rho, grho2 ), copyout( ex, v1x, v2x )
       !$acc host_data use_device( rho, grho2, ex, v1x, v2x )
       CALL gcx_spin( length, rho, grho2, ex, v1x, v2x )
       !$acc end host_data
       !$acc end data
       !
    ENDIF
    !
  ENDIF
  !
  DEALLOCATE( rho_lxc, sigma )
  DEALLOCATE( ex_lxc , ec_lxc   )
  DEALLOCATE( vx_rho , vx_sigma )
  DEALLOCATE( vc_rho , vc_sigma )
  !
#else
  !
  ALLOCATE( grho2(length,ns) )
  grho2 = 0.0_DP
  !
  IF ( ns == 1 ) THEN
     !
     DO k = 1, length
        IF ( ABS(rho(k,1)) > rho_threshold_gga ) &
          grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
     ENDDO
     !
     !$acc data copyin( rho, grho2 ), copyout( ex, ec, v1x, v2x, v1c, v2c )
     !$acc host_data use_device( rho, grho2, ex, ec, v1x, v2x, v1c, v2c )
     CALL gcxc( length, rho(:,1), grho2(:,1), ex, ec, v1x(:,1), v2x(:,1), &
                v1c(:,1), v2c(:,1) )
     !$acc end host_data
     !$acc end data
     !
     DO k = 1, length
        sgn(1) = SIGN(1._DP, rho(k,1))
        ex(k) = ex(k) * sgn(1)
        ec(k) = ec(k) * sgn(1)
     ENDDO
     !
  ELSE
     !
     DO is = 1, 2
        grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
     ENDDO
     !
     !$acc data copyin( rho, grho2 ), copyout( ex, v1x, v2x )
     !$acc host_data use_device( rho, grho2, ex, v1x, v2x )
     CALL gcx_spin( length, rho, grho2, ex, v1x, v2x )
     !$acc end host_data
     !$acc end data
     !
     IF (igcc==3 .OR. igcc==7 .OR. igcc==13 ) THEN
        !
        ALLOCATE( grho_ud(length) )
        !
        grho_ud = grho(1,:,1) * grho(1,:,2) + grho(2,:,1) * grho(2,:,2) + &
                  grho(3,:,1) * grho(3,:,2)
        !
        !$acc data copyin( rho, grho2, grho_ud ), copyout( ec, v1c, v2c, v2c_ud )
        !$acc host_data use_device( rho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
        CALL gcc_spin_more( length, rho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
        !$acc end host_data
        !$acc end data
        !
        DEALLOCATE( grho_ud )
        !
     ELSE
        !
        ALLOCATE( rh(length), zeta(length) )
        !
        rh = rho(:,1) + rho(:,2)
        !
        zeta = 2.0_DP ! trash value, gcc-routines get rid of it when present
        WHERE ( rh > rho_threshold_gga ) zeta = ( rho(:,1) - rho(:,2) ) / rh(:)
        !
        grho2(:,1) = ( grho(1,:,1) + grho(1,:,2) )**2 + &
                     ( grho(2,:,1) + grho(2,:,2) )**2 + &
                     ( grho(3,:,1) + grho(3,:,2) )**2
        !
        IF ( igcc/=14 ) THEN
          !$acc data copyin( rh, zeta, grho2 ), copyout( ec, v1c, v2c )
          !$acc host_data use_device( rh, zeta, grho2, ec, v1c, v2c )
          CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
          !$acc end host_data
          !$acc end data
        ENDIF
        !
        v2c(:,2)  = v2c(:,1)
        IF ( ns==2 ) v2c_ud(:) = v2c(:,1)
        !
        DEALLOCATE( rh, zeta )
        !
     ENDIF
     !   
  ENDIF
  !
  DEALLOCATE( grho2 )
  !
#endif
  !
  !
  RETURN
  !
END SUBROUTINE xc_gcx_
!
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gcx_gpu( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
  !-------------------------------------------------------------------------
  !! GGA wrapper routine - gpu double.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_setting_params,   ONLY: xc_func, xc_info, libxc_flags
#endif
  !
  USE kind_l,               ONLY: DP
  USE dft_setting_params,   ONLY: igcx, igcc, is_libxc, rho_threshold_gga, &
                                  grho_threshold_gga, rho_threshold_lda
  USE qe_drivers_gga
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  INTEGER,  INTENT(IN) :: ns
  REAL(DP), INTENT(IN) :: rho(length,ns)
  REAL(DP), INTENT(IN) :: grho(3,length,ns)
  REAL(DP), INTENT(OUT) :: ex(length)
  REAL(DP), INTENT(OUT) :: ec(length)
  REAL(DP), INTENT(OUT) :: v1x(length,ns)
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  REAL(DP), INTENT(OUT) :: v2c(length,ns)
  REAL(DP), INTENT(OUT) :: v2c_ud(length)
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  REAL(DP), ALLOCATABLE :: rho_lxc(:), sigma(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_rho(:), vx_sigma(:)
  REAL(DP), ALLOCATABLE :: vc_rho(:), vc_sigma(:)
  !
  INTEGER :: fkind_x, np
  REAL(DP) :: rs, rtot, zet, vc_2(2), arho_k
  REAL(DP), PARAMETER :: pi34 = 0.6203504908994_DP
  !
  LOGICAL :: POLARIZED
  INTEGER :: ildax, ildac, pol_unpol
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  REAL(DP), ALLOCATABLE :: rh(:), zeta(:)
  REAL(DP), ALLOCATABLE :: grho2(:,:), grho_ud(:)
  !
  INTEGER :: k, is
  REAL(DP) :: rho_up, rho_dw, grho_up, grho_dw, sgn1
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  !$acc data deviceptr( rho(length,ns), grho(3,length,ns), ex(length), ec(length),      &
  !$acc&                v1x(length,ns), v2x(length,ns), v1c(length,ns), v2c(length,ns), &
  !$acc&                v2c_ud(length) )
  !
  !ex = 0.0_DP ;  v1x = 0.0_DP ;  v2x = 0.0_DP
  !ec = 0.0_DP ;  v1c = 0.0_DP ;  v2c = 0.0_DP
  !IF ( PRESENT(v2c_ud) ) v2c_ud = 0.0_DP
  !
  !
#if defined(__LIBXC)
  !
  fkind_x = -1
  lengthxc = length
  !
  POLARIZED = .FALSE.
  IF (ns == 2) THEN
     POLARIZED = .TRUE.
  ENDIF
  !
  pol_unpol = 1
  np = 1
  IF ( ns == 2 ) THEN
     pol_unpol = 2
     np = 3
  ENDIF
  !
  ALLOCATE( rho_lxc(length*ns) )
  ALLOCATE( sigma(length*np) )
  !
  ALLOCATE( ex_lxc(length),    ec_lxc(length)      )
  ALLOCATE( vx_rho(length*ns), vx_sigma(length*np) )
  ALLOCATE( vc_rho(length*ns), vc_sigma(length*np) )
  !$acc data copyout( rho_lxc, sigma )
  !$acc host_data use_device( rho_lxc, sigma )
  !
  IF ( ns == 1 ) THEN
    !
    !$acc parallel loop
    DO k = 1, length
      arho_k = ABS(rho(k,1))
      rho_lxc(k) = arho_k
      IF ( arho_k > rho_threshold_gga ) &
         sigma(k) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
    ENDDO
    !
  ELSE
    !
    !$acc parallel loop
    DO k = 1, length
       rho_lxc(2*k-1) = rho(k,1)
       rho_lxc(2*k)   = rho(k,2)
       !
       sigma(3*k-2) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
       sigma(3*k-1) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                      grho(3,k,1) * grho(3,k,2)
       sigma(3*k)   = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
    ENDDO
    !
  ENDIF
  !
  IF ( ns==1 .AND. ANY(.NOT.is_libxc(3:4)) ) THEN
     !
     CALL gcxc( length, rho_lxc, sigma, ex, ec, v1x(:,1), v2x(:,1), v1c(:,1), v2c(:,1) )
     !
     !$acc parallel loop
     DO k = 1, length
        sgn1 = SIGN(1._DP, rho(k,1))
        ex(k) = ex(k) * sgn1
        ec(k) = ec(k) * sgn1
     ENDDO
     !
  ENDIF
  !$acc end host_data
  !$acc end data
  !
  ! ---- GGA CORRELATION
  !
  IF ( is_libxc(4) ) THEN  !lda part of LYP not present in libxc (still so? - check)
    !
    fkind_x = xc_f03_func_info_get_kind( xc_info(4) )
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(4), small )!rho_threshold_gga )
    IF (libxc_flags(4,0)==1) THEN
      CALL xc_f03_gga_exc_vxc( xc_func(4), lengthxc, rho_lxc(1), sigma(1), ec_lxc(1), vc_rho(1), vc_sigma(1) )
    ELSE
      CALL xc_f03_gga_vxc( xc_func(4), lengthxc, rho_lxc(1), sigma(1), vc_rho(1), vc_sigma(1) )
      ec_lxc = 0.d0
    ENDIF
    !
    !$acc data copyin( rho_lxc, sigma, ec_lxc, vc_rho, vc_sigma )
    IF (.NOT. POLARIZED) THEN
      !$acc parallel loop
      DO k = 1, length
        IF ( rho_lxc(k) <= rho_threshold_lda ) THEN
          ec(k)=0.d0 ;  v1c(k,1)=0.d0 ;  v2c(k,1)=0.d0
          CYCLE
        ENDIF  
        ec(k) = ec_lxc(k) * rho_lxc(k) * SIGN(1.0_DP, rho(k,1))
        v1c(k,1) = vc_rho(k)
        IF ( rho_lxc(k) <= rho_threshold_gga .OR. &
             SQRT(ABS(sigma(k))) <= grho_threshold_gga) THEN
           v2c(k,1) = 0.d0
           CYCLE
        ENDIF
        v2c(k,1) = vc_sigma(k)*2.d0
      ENDDO
    ELSE
      !$acc parallel loop
      DO k = 1, length
        rho_up = rho_lxc(2*k-1)
        rho_dw = rho_lxc(2*k)
        grho_up = SQRT(ABS(sigma(3*k-2)))
        grho_dw = SQRT(ABS(sigma(3*k)))
        IF ( rho_up <= rho_threshold_lda .OR. rho_dw <= rho_threshold_lda ) THEN
          ec(k) = 0.d0    ;  v1c(k,1) = 0.d0 ;  v1c(k,2) = 0.d0
          v2c(k,1) = 0.d0 ;  v2c_ud(k)= 0.d0 ;  v2c(k,2) = 0.d0
          CYCLE
        ENDIF
        ec(k) = ec_lxc(k) * (rho_up+rho_dw)
        v1c(k,1) = vc_rho(2*k-1)
        v1c(k,2) = vc_rho(2*k)
        IF ( rho_up <= rho_threshold_gga .OR. rho_dw <= rho_threshold_gga .OR. &
             grho_up<=grho_threshold_gga .OR. grho_dw<=grho_threshold_gga ) THEN
          v2c(k,1) = 0.d0 ; v2c_ud(k)= 0.d0 ; v2c(k,2) = 0.d0
          CYCLE
        ENDIF
        v2c(k,1) = vc_sigma(3*k-2)*2.d0
        v2c_ud(k)= vc_sigma(3*k-1)
        v2c(k,2) = vc_sigma(3*k)*2.d0
      ENDDO
    ENDIF
    !$acc end data
    !
  ELSEIF ( (.NOT.is_libxc(4)) .AND. fkind_x/=XC_EXCHANGE_CORRELATION ) THEN
    !
    ALLOCATE( grho2(length,ns) )
    !
    IF ( ns /= 1 ) THEN
       !
       IF (igcc==3 .OR. igcc==7 .OR. igcc==13 ) THEN
          !
          ALLOCATE( grho_ud(length) )
          !$acc data create( grho2(length,ns), grho_ud(length) )
          !$acc host_data use_device( grho2, grho_ud )
          !$acc parallel loop
          DO k = 1, length
            grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
            grho2(k,2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
            grho_ud(k) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                         grho(3,k,1) * grho(3,k,2)
          ENDDO
          CALL gcc_spin_more( length, rho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
          !$acc end host_data
          !$acc end data
          DEALLOCATE( grho_ud )
          !
       ELSE
          !
          ALLOCATE( rh(length), zeta(length) )
          !$acc data create( rh(length), zeta(length), grho2(length,ns))
          !$acc host_data use_device( rh, zeta, grho2 )
          !$acc parallel loop
          DO k = 1, length
            rh(k) = rho(k,1) + rho(k,2)
            IF ( rh(k) > rho_threshold_gga ) THEN
              zeta(k) = ( rho(k,1) - rho(k,2) ) / rh(k)
            ELSE
              zeta(k) = 2.0_DP ! trash value, gcc-routines get rid of it when present
            ENDIF
            grho2(k,1) = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                         ( grho(2,k,1) + grho(2,k,2) )**2 + &
                         ( grho(3,k,1) + grho(3,k,2) )**2
            grho2(k,2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
          ENDDO
          !
          CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
          !
          !$acc parallel loop
          DO k = 1, length
            v2c(k,2) = v2c(k,1)
            IF ( ns==2 ) v2c_ud(k) = v2c(k,1)
          ENDDO
          !$acc end host_data
          !$acc end data
          DEALLOCATE( rh, zeta )
          !
       ENDIF
       !
    ENDIF
    !
    DEALLOCATE( grho2 )
    !
  ENDIF
  !
  ! --- GGA EXCHANGE
  !
  IF ( is_libxc(3) ) THEN
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(3), grho_threshold_gga )
    IF (libxc_flags(3,0)==1) THEN
      CALL xc_f03_gga_exc_vxc( xc_func(3), lengthxc, rho_lxc(1), sigma(1), ex_lxc(1), vx_rho(1), vx_sigma(1) )
    ELSE
      CALL xc_f03_gga_vxc( xc_func(3), lengthxc, rho_lxc(1), sigma(1), vx_rho(1), vx_sigma(1) )
      ex_lxc = 0.d0
    ENDIF
    !
    !$acc data copyin( rho_lxc, sigma, ex_lxc, vx_rho, vx_sigma )
    IF (.NOT. POLARIZED) THEN
      !$acc parallel loop
      DO k = 1, length
        IF ( rho_lxc(k) <= rho_threshold_lda ) THEN
          ex(k) = 0.d0 ;  v1x(k,1) = 0.d0 ;  v2x(k,1) = 0.d0
          CYCLE
        ENDIF
        ex(k) = ex_lxc(k) * rho_lxc(k) * SIGN(1.0_DP, rho(k,1))
        v1x(k,1) = vx_rho(k)
        IF ( rho_lxc(k) <= rho_threshold_gga .OR. &
             SQRT(ABS(sigma(k))) <= grho_threshold_gga) THEN
          v2x(k,1) = 0.d0
          CYCLE
        ENDIF
        v2x(k,1) = vx_sigma(k)*2.d0
      ENDDO
    ELSE
      !$acc parallel loop
      DO k = 1, length
        rho_up = rho_lxc(2*k-1)
        rho_dw = rho_lxc(2*k)
        grho_up = SQRT(ABS(sigma(3*k-2)))
        grho_dw = SQRT(ABS(sigma(3*k)))
        IF ( rho_up <= rho_threshold_lda .OR. rho_dw <= rho_threshold_lda ) THEN
          ex(k) = 0.d0 
          v1x(k,1) = 0.d0 ;  v1x(k,2) = 0.d0
          v2x(k,1) = 0.d0 ;  v2x(k,2) = 0.d0
          CYCLE
        ENDIF
        ex(k) = ex_lxc(k) * (rho_up+rho_dw)
        v1x(k,1) = vx_rho(2*k-1)
        v1x(k,2) = vx_rho(2*k)
        IF ( rho_up <= rho_threshold_gga .OR. rho_dw <= rho_threshold_gga .OR. &
             grho_up<=grho_threshold_gga .OR. grho_dw<=grho_threshold_gga ) THEN
          v2x(k,1) = 0.d0 ;  v2x(k,2) = 0.d0
          CYCLE
        ENDIF     
        v2x(k,1) = vx_sigma(3*k-2)*2.d0
        v2x(k,2) = vx_sigma(3*k)*2.d0
      ENDDO
    ENDIF
    !$acc end data
    !
  ELSE
    !
    IF ( ns /= 1 ) THEN
      !
      ALLOCATE( grho2(length,ns) )
      !$acc data create( grho2(length,ns) )
      !$acc host_data use_device( grho2 )
      !$acc parallel loop collapse(2)
      DO is = 1, ns
        DO k = 1, length
          grho2(k,is) = grho(1,k,is)**2 + grho(2,k,is)**2 + grho(3,k,is)**2
        ENDDO
      ENDDO
      CALL gcx_spin( length, rho, grho2, ex, v1x, v2x )
      !$acc end host_data
      !$acc end data
      DEALLOCATE( grho2 )
      !
    ENDIF
    !
  ENDIF
  !
  DEALLOCATE( rho_lxc, sigma )
  DEALLOCATE( ex_lxc , ec_lxc   )
  DEALLOCATE( vx_rho , vx_sigma )
  DEALLOCATE( vc_rho , vc_sigma )
  !
#else
  !
  ALLOCATE( grho2(length,ns) )
  !$acc data create( grho2(length,ns) )
  !
  IF ( ns == 1 ) THEN
     !
     ALLOCATE( rh(length) )
     !$acc data create( rh(length) )
     !$acc host_data use_device( rh, grho2 )
     !$acc parallel loop
     DO k = 1, length
        rh(k) = ABS(rho(k,1))
        IF ( rh(k) > rho_threshold_gga ) &
             grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
     ENDDO
     !
     CALL gcxc( length, rh, grho2(:,1), ex, ec, v1x(:,1), v2x(:,1), v1c(:,1), v2c(:,1) )
     !
     !$acc parallel loop
     DO k = 1, length
        sgn1 = SIGN(1._DP, rho(k,1))
        ex(k) = ex(k) * sgn1
        ec(k) = ec(k) * sgn1
     ENDDO
     !
     !$acc end host_data
     !$acc end data
     DEALLOCATE( rh )
     !
  ELSE
     !
     !$acc parallel loop collapse(2)
     DO is = 1, ns
       DO k = 1, length
         grho2(k,is) = grho(1,k,is)**2 + grho(2,k,is)**2 + grho(3,k,is)**2
       ENDDO
     ENDDO
     !
     CALL gcx_spin( length, rho, grho2, ex, v1x, v2x )
     !
     IF (igcc==3 .OR. igcc==7 .OR. igcc==13 ) THEN
        !
        ALLOCATE( grho_ud(length) )
        !$acc data create( grho_ud(length) )
        !$acc host_data use_device( grho_ud, grho2 )
        !$acc parallel loop
        DO k = 1, length
          grho_ud(k) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                       grho(3,k,1) * grho(3,k,2)
        ENDDO
        !
        CALL gcc_spin_more( length, rho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
        !
        !$acc end host_data
        !$acc end data
        DEALLOCATE( grho_ud )
        !
     ELSE
        !
        ALLOCATE( rh(length), zeta(length) )
        !$acc data create( rh(length), zeta(length) )
        !$acc host_data use_device( rh, zeta, grho2 )
        !$acc parallel loop
        DO k = 1, length
          rh(k) = rho(k,1) + rho(k,2)
          IF ( rh(k) > rho_threshold_gga ) THEN
            zeta(k) = ( rho(k,1) - rho(k,2) ) / rh(k)
          ELSE
            zeta(k) = 2.0_DP ! trash value, gcc-routines get rid of it when present
          ENDIF
          grho2(k,1) = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                       ( grho(2,k,1) + grho(2,k,2) )**2 + &
                       ( grho(3,k,1) + grho(3,k,2) )**2
        ENDDO
        !
        CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
        !
        !$acc parallel loop
        DO k = 1, length
          v2c(k,2)  = v2c(k,1)
          IF ( ns==2 ) v2c_ud(k) = v2c(k,1)
        ENDDO
        !
        !$acc end host_data
        !$acc end data
        DEALLOCATE( rh, zeta )
        !
     ENDIF
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( grho2 )
  !
#endif
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE xc_gcx_gpu
