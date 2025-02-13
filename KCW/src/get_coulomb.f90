
SUBROUTINE write_coulomb()
USE kinds,                ONLY : DP
USE control_kcw,          ONLY : Vcoulomb, Wcoulomb
USE control_kcw,          ONLY : spin_component, get_coulomb, rvect_shifted
USE control_kcw,          ONLY : num_wann, nqstot
!
IMPLICIT NONE
!
INTEGER              ::    iun_coulomb
CHARACTER(len=1024)  ::    filename  
INTEGER              ::    iwann, jwann, ir, is
INTEGER              ::    spin_index

!
!
  filename = 'barecoulomb.txt'
  iun_coulomb = 237 
  OPEN (iun_coulomb, file = filename)

  DO ir = 1, nqstot
    WRITE(iun_coulomb, *) rvect_shifted(:, ir) 
    DO iwann=1, num_wann
      DO jwann=1, num_wann
        DO is = 1, 2 !one for spin component, other for non spin component
          WRITE(iun_coulomb, *) iwann, jwann, spin_index(is), Vcoulomb(ir, iwann, jwann, is) 
        END DO!is
      END DO!jwann
    END DO!iwann
  END DO!ir
  WRITE (filename, *) 'barecoulomb.txt'
  CLOSE(iun_coulomb)
  !
  iun_coulomb = 238 
  filename = 'screencoulomb.txt'
  OPEN (iun_coulomb, file = filename)
  !
  !
  DO ir = 1, nqstot
    WRITE(iun_coulomb, *) rvect_shifted(:, ir) 
    DO iwann=1, num_wann
      DO jwann=1, num_wann
        DO is = 1, 2 !one for spin component, other for non spin component
          WRITE(iun_coulomb, *) iwann, jwann, spin_index(is), Wcoulomb(ir, iwann, jwann, is) 
        END DO!is
      END DO!jwann
    END DO!iwann
  END DO!ir
!
END SUBROUTINE write_coulomb

SUBROUTINE coulomb_me( iwann, iq, drhog_scf, rhowann, weight_q)

!this function calculates the Coulomb interaction matrix element:
! <rho_{R,iwann}|V_{Hxc}|rho_{0,jwann}>
!
USE kinds,                ONLY : DP
USE control_kcw,          ONLY : nqstot, tmp_dir_save, num_wann, nrho
USE control_kcw,          ONLY : iurho_wann, x_q, rvect_shifted
USE control_kcw,          ONLY : Vcoulomb, Wcoulomb
USE control_kcw,          ONLY : spin_component, get_coulomb
USE buffers,              ONLY : get_buffer 
USE io_files,             ONLY : tmp_dir
USE fft_base,             ONLY : dffts
USE klist,                ONLY : nkstot
USE gvecs,                ONLY : ngms
USE noncollin_module,     ONLY : nspin_mag
USE constants,            ONLY : tpi
USE lsda_mod,             ONLY : nspin
USE cell_base,            ONLY : omega
!
IMPLICIT NONE
COMPLEX(DP)                  :: IMAG = (0.D0,1.D0)
INTEGER, INTENT(IN)          :: iwann
! (fixed) index of wannier function
INTEGER, INTENT(IN)          :: iq
! (fixed) index of q point
! quantity V_xc|rho_i> in G space
COMPLEX(DP), INTENT(IN)      :: drhog_scf(ngms,nspin_mag)
!LR wannier density for fixed q and iwann
COMPLEX(DP), INTENT(IN)      :: rhowann(dffts%nnr, num_wann, nrho)
INTEGER                      :: jwann
! index of wannier functions to loop over
REAL(DP), INTENT(IN)         :: weight_q 
INTEGER                      :: ip
! index of nrho to loop over
INTEGER                      :: ir
!
!quantities needed for computing bare_pot for jwann
!
COMPLEX(DP)                  :: rhog(ngms,nrho)
! ... periodic part of wannier density in G-space
COMPLEX(DP)                  :: delta_vg(ngms,nspin_mag)
! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space
COMPLEX(DP)                  :: delta_vg_(ngms,nspin_mag)
! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space without g=0 contribution
COMPLEX(DP)                  :: vh_rhog(ngms)
! ... Hartree perturbing potential [f_hxc(r,r') x wann(r')] in G-space
COMPLEX(DP)                  :: delta_vr(dffts%nnr,nspin_mag)
! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space
COMPLEX(DP)                  :: delta_vr_(dffts%nnr,nspin_mag)
! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space without g=0 contribution
COMPLEX(DP)                  :: rhor(dffts%nnr,nrho)
! ... periodic part of wannier density in r-space
INTEGER                      :: lrrho
INTEGER                      :: is, is_, is1
INTEGER                      :: spin_index
!
! we are already inside loops over iq and iwann
!
COMPLEX(DP)                  :: pi_q_unrelax(2), pi_q_relax(2), pi_q_unrelax_(2)
!
IF (nrho==4) THEN
    CALL errore('output_coulomb', 'non-collinear not implemented &
      for coulomb matrix elements.', 1)
END IF
DO jwann = 1, num_wann
  !here rhowann is already filled with the wannier density in real space 
  !for fixed iq
  !
  ! get rhog for jwann by fourier transforming rhor 
  !WARNING! WON'T WORK FOR NON_COLLINEAR
  DO ip = 1, nrho
    rhor(:,ip) = rhowann(:,jwann,ip)
    !
    DO ir = 1, nkstot/nspin
      WRITE(*,*) x_q(:,iq)
      WRITE(*,*) rvect_shifted(:,ir)
      rhor(:, ip) = rhor(:, ip) * EXP( IMAG*tpi*DOT_PRODUCT(x_q(:,iq),rvect_shifted(:,ir)) )
      CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_) 
      !
      ! for now we only work with the density, i.e. ip = 1
      ! we evaluate 
      !   <jwann, R| Vxc | iwann, 0> = <deltaVg | rhog(iwann)>
      DO is = 1, 2 !one for spin component, other for non spin component
          pi_q_unrelax (spin_index(is)) = weight_q * omega * SUM( CONJG(delta_vg (:,is)) * rhog(:,1) )  
          pi_q_unrelax_(spin_index(is)) = weight_q * omega * SUM( CONJG(delta_vg_(:,is)) * rhog(:,1) )  
      END DO  
      Vcoulomb(iwann, jwann, ir, :) = Vcoulomb(iwann, jwann, ir, :) + pi_q_unrelax(:)
      !
      ! screened potential
      ! < jwann, R| W | 0, iwann>
      pi_q_relax(:) = (0.D0, 0.D0)
      DO is =1, 2
        DO is1 = 1, 2
          pi_q_relax(spin_index(is)) = pi_q_relax(spin_index(is)) + &
                                       sum (CONJG(drhog_scf (:,is1)) * delta_vg(:,is_(is1)))*weight_q*omega
        END DO 
      END DO
      Wcoulomb(iwann, jwann, ir, :) = Wcoulomb(iwann, jwann, ir, :) + pi_q_unrelax_(:) + pi_q_relax(:)
    END DO !ir
  END DO !ip
END DO !jwann
END SUBROUTINE



FUNCTION spin_index(is) 
  !! spin index for up-up / down-down -> 1 
  !! spin index for up-down /down-up  -> 2
  USE control_kcw,     ONLY : spin_component
  IMPLICIT NONE
  INTEGER    :: is 
  INTEGER    :: spin_index 

  IF( is .eq. spin_component ) THEN 
    spin_index = 1
  ELSE 
    spin_index = 2 
  END IF
END FUNCTION

FUNCTION is_(is)        
  USE control_kcw,     ONLY : spin_component
  IMPLICIT NONE
  INTEGER    :: is 
  INTEGER    :: is_
  IF ( is .eq. spin_component ) THEN 
    is_ = is
  ELSE 
    IF ( is .eq. 1 ) THEN 
      is_ = 2
    ELSE 
      is_ = 1
    END IF
  END IF
END FUNCTION