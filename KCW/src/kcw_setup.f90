!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
!-----------------------------------------------------------------------
subroutine kcw_setup
  !-----------------------------------------------------------------------
  !
  !!  This subroutine prepares several variables which are needed in the
  !!  KC calculation:
  !!  * computes the total local potential (external+scf) on the smooth
  !!    grid to be used in h_psi and similia
  !!  * set the inverse of every matrix invs
  !!  * for metals sets the occupied bands
  !!  * Open buffer for the KS and eventualy the minimizing WFCs
  !!  * Open buffer for the KS states in the WANNIER gauge
  !!  * Read the U matrix from Wannier
  !!  * Rotate the KS state to the localized gauge
  !!  * Compute the periodic part of the wannier orbital and store in the buffer 
  !!    a separate directory is generated for each q 
  !
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ityp
  USE io_files,          ONLY : tmp_dir
  USE lsda_mod,          ONLY : nspin, starting_magnetization, lsda, isk
  USE scf,               ONLY : v, vrs, vltot,  kedtau
  USE fft_base,          ONLY : dfftp, dffts
  USE gvecs,             ONLY : doublegrid, ngms
  USE gvect,             ONLY : ig_l2g, mill
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux!, nspin_mag, npol
  USE wvfct,             ONLY : nbnd
  !USE funct,             ONLY : dft_is_gradient
  USE xc_lib,            ONLY : xclib_dft_is
  !
  USE units_lr,          ONLY : iuwfc
  USE wvfct,             ONLY : npwx
  USE control_flags,     ONLY : io_level, gamma_only
  USE io_files,          ONLY : prefix
  USE buffers,           ONLY : open_buffer, save_buffer, close_buffer, get_buffer
  USE control_kcw,       ONLY : evc0, iuwfc_wann, iuwfc_wann_allk, kcw_iverbosity, lgamma_iq, &
                                spin_component, isq, read_unitary_matrix, x_q, tmp_dir_save, & 
                                num_wann, num_wann_occ, occ_mat, tmp_dir_kcw, tmp_dir_kcwq, &
                                io_sp, io_real_space, nqstot !, wq, nqstot
  USE io_global,         ONLY : stdout
  USE klist,             ONLY : nkstot, xk
  USE cell_base,         ONLY : at, omega, bg, tpiba
  USE fft_base,          ONLY : dffts
  !
  USE control_lr,       ONLY : nbnd_occ, lgamma
  USE scf,              ONLY : rho
  USE io_global,        ONLY : ionode, ionode_id
  USE mp_images,        ONLY : intra_image_comm
  USE mp,               ONLY : mp_bcast
  USE io_files,         ONLY : create_directory
  USE io_rho_xml,       ONLY : write_scf
  !
  USE io_kcw,           ONLY : write_rhowann, write_rhowann_g
  !
  USE mp,               ONLY : mp_sum
  USE control_lr,       ONLY : lrpa
  USE coulomb,          ONLY : setup_coulomb
  !
  USE mp_pools,         ONLY : my_pool_id
  USE mp_bands,         ONLY : my_bgrp_id, root_bgrp_id, &
                               root_bgrp, intra_bgrp_comm, &
                               inter_bgrp_comm
  USE symm_base,        ONLY : s, time_reversal, sr, ft
  USE control_kcw,      ONLY : fbz2ibz, wq_ibz, irr_bz
  USE io_kcw,                ONLY : read_rhowann, read_rhowann_g
  USE fft_interfaces,        ONLY : invfft, fwfft

  !
  implicit none
  !
  integer :: na, i, ik, iq_ibz
  ! counters
  !
  INTEGER   :: lrwfc, iun_qlist, iun_qlist_ibz
  LOGICAL   :: exst, exst_mem
  INTEGER :: iq, nqs
  REAL(DP) :: xq(3)
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:,:), rhowann_aux(:)
  CHARACTER (LEN=256) :: filename, file_base
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  !
  INTEGER :: &
       igk_k_all(npwx,nkstot),&    ! index of G corresponding to a given index of k+G
       ngk_all(nkstot)             ! number of plane waves for each k point
  !
  ! Auxiliary variables for SH calculation
  COMPLEX(DP), ALLOCATABLE  :: rhog(:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:), sh(:)
  ! The periodic part of the wannier orbital density
  COMPLEX(DP), ALLOCATABLE  :: rhowann_g(:,:)
  COMPLEX(DP) :: rhor(dffts%nnr), delta_vr(dffts%nnr,nspin), delta_vr_(dffts%nnr,nspin)
  !
  ! The weight of each q point
  REAL(DP), ALLOCATABLE :: weight(:)
  LOGICAL :: lrpa_save
  !
  ALLOCATE ( rhog (ngms) , delta_vg(ngms,nspin), vh_rhog(ngms), delta_vg_(ngms,nspin) )
  !
  call start_clock ('kcw_setup')
  !
  ! ... Computes the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  !  ... If necessary calculate the local magnetization. This information is
  !      needed in find_sym
  !
  IF (.not.ALLOCATED(m_loc)) ALLOCATE( m_loc( 3, nat ) )
  IF (noncolin.and.domag) THEN
     DO na = 1, nat
        !
        m_loc(1,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
        m_loc(2,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
        m_loc(3,na) = starting_magnetization(ityp(na)) * &
                      COS( angle1(ityp(na)) )
     END DO
     ux=0.0_DP
     IF ( xclib_dft_is('gradient') ) CALL compute_ux(m_loc,ux,nat)
     !if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
  ENDIF
  !
  ! ... Computes the number of occupied bands for each k point
  !
  call setup_nbnd_occ ( ) 
  !
  ! ... Open buffers for the KS and eventualy the minimizing WFCs 
  ! 
  iuwfc = 20; lrwfc = nbnd * npwx; io_level = 1
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  !
  IF (.NOT.exst.AND..NOT.exst_mem) &
       CALL errore ('kcw_setup', 'file '//trim(prefix)//'.wfc not found', 1)
  !
  if (kcw_iverbosity .gt. 1) &
       WRITE(stdout,'(/,5X, "INFO: Buffer for KS wfcs, OPENED")')
  !
  ! ... READ the U matrix from Wannier and set the toal number of WFs
  !
  IF (read_unitary_matrix) THEN 
    !
    CALL read_wannier ( )
    if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Unitary matrix, READ from file")')
    !
  ELSE
    !
    num_wann = nbnd
    num_wann_occ = nbnd_occ(1) ! Assuming insulating
    !
  ENDIF
  ! 
  !  ... Allocate relevant quantities ...
  !
  ALLOCATE (rhowann ( dffts%nnr, nkstot/nspin, num_wann), rhowann_aux(dffts%nnr) )
  ALLOCATE ( rhowann_g (ngms, num_wann) )
  ALLOCATE ( evc0(npwx, num_wann) )
  ALLOCATE ( occ_mat (num_wann, num_wann, nkstot) )
  ALLOCATE ( sh(num_wann) )
  sh(:) = CMPLX(0.D0,0.D0,kind=DP)
  occ_mat = 0.D0
  !
  ! ... Open a new buffer to store the KS states in the WANNIER gauge
  !
  iuwfc_wann = 21
  io_level = 1
  lrwfc = num_wann * npwx 
  CALL open_buffer ( iuwfc_wann, 'wfc_wann', lrwfc, io_level, exst )
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WFs, OPENED")')
  !
  ! ... Open an other buffer for the KS states in the WANNIER gauge which contains
  !     all the k points (not just the one in this pool). This is needed for each k-point 
  !     to have access to all the other k-points. MEMORY INTENSE
  !
  iuwfc_wann_allk = 210
  io_level = 1
  lrwfc = num_wann * npwx
  CALL open_buffer ( iuwfc_wann_allk, 'wfc_wann_allk', lrwfc, io_level, exst )
  !
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WFs ALL-k, OPENED")')
  !
  ! ... Rotate the KS state to the localized gauge nd save on a buffer
  !
  CALL rotate_ks () 
  !
  ! ... pass all the WFs to all the pool (needed to have pool parallelization)
  !
  CALL bcast_wfc ( igk_k_all, ngk_all )
  !
  DEALLOCATE ( nbnd_occ )  ! otherwise allocate_ph complains: FIXME
  ! 8)
  !CALL compute_map_ikq ()
  !
  ! ... Compute the periodic part of the wannier orbital and store it on file 
  !
  WRITE(stdout,'(/)')
  WRITE( stdout, '(5X,"INFO: PREPARING THE KCW CALCULATION ...")')
  !
  !     Write the list of q points on a file and store the q coordinates
  !
  iun_qlist = 127
  OPEN (iun_qlist, file = TRIM(tmp_dir_kcw)//'qlist.txt')
  !
  nqs = nkstot/nspin
  nqstot = nqs ! FIXME: replace nqs (local) with nqstot (global)
  ALLOCATE (x_q (3, nqs), isq(nqs) )
  ALLOCATE (weight(nqs) )
  iq=1
  IF (ionode) THEN 
     WRITE(iun_qlist,'(i5)') nkstot/nspin 
     DO ik = 1, nkstot
       IF (lsda .AND. isk(ik) /= spin_component) CYCLE
       WRITE(iun_qlist, '(3f12.8)') xk(:,ik)
       x_q(:,iq) = xk(:,ik)  
       isq(iq) = isk(ik) 
       iq = iq + 1
     ENDDO 
  ENDIF
  CALL mp_bcast (x_q, ionode_id, intra_image_comm)
  CALL mp_bcast (isq, ionode_id, intra_image_comm)
  !
  ALLOCATE ( lgamma_iq(nqs) )
  lgamma_iq(:) = .FALSE.
  !
  WRITE( stdout, '(/, 5X,"INFO: Compute Wannier-orbital Densities ...")')
  !
  call setup_coulomb()
  !
  rhowann(:,:,:)=ZERO
  !! Initialize the periodic part of the wannier orbtal density at this q
  !
  DO iq = 1, nqs
    !! For each q in the mesh 
    !
    xq = x_q(:,iq)
    !
    ! IF (ionode) WRITE(iun_qlist,'(3f12.8)') xq
    !
    lgamma_iq(iq)=(x_q(1,iq)==0.D0.AND.x_q(2,iq)==0.D0.AND.x_q(3,iq)==0.D0)
    CALL cryst_to_cart(1, xq, at, -1)
    WRITE( stdout, '(/,8X, 78("="))')
    WRITE( stdout, '(  8X, "iq = ", i5)') iq
    WRITE( stdout, '(  8X, "The Wannier density at  q = ",3F12.7, "  [Cart ]")') x_q(:,iq)
    WRITE( stdout, '(  8X, "The Wannier density at  q = ",3F12.7, "  [Cryst]")') xq(:)
    WRITE( stdout, '(  8X, 78("="),/)')
    !
    CALL compute_map_ikq_single (iq)
    ! The map to identify which k point in the 1BZ corresponds to k+q and the G vector that produce the mapping
    ! The results are stored in the global variable map_ikq and shift_1bz (used inside rho_of_q) 
    ! can (should) be moved inside rho_of_q ( )
    !
    !
    CALL rho_of_q (rhowann(:,iq,:), ngk_all, igk_k_all)
    ! Compute the peridic part rho_q(r) of the wannier density rho(r)
    ! rho(r)   = \sum_q exp[iqr]rho_q(r)
    !
    WRITE( stdout, '(8X,"INFO: rho_q(r) DONE ",/)')
    !
    ! Compute the Self Hartree
    weight(iq) = 1.D0/nqs ! No SYMM 
    lrpa_save=lrpa
    lrpa = .true.
    DO i = 1, num_wann
      !
      rhog(:)         = CMPLX(0.D0,0.D0,kind=DP)
      delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
      vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
      rhor(:)         = CMPLX(0.D0,0.D0,kind=DP)
      !
      rhor(:) = rhowann(:,iq,i) 
      !! The periodic part of the orbital desity in real space
      !
      CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
      rhowann_g(:,i) = rhog
      !! The periodic part of the perturbation DeltaV_q(G)
      ! 
      sh(i) = sh(i) + 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:) )*weight(iq)*omega
      !WRITE(1005,*) "iwann=", i, "q=", iq, "SH=", REAL(.5*sum (CONJG(rhog (:)) * vh_rhog(:) )*weight(iq)*omega)
      !
    ENDDO
    lrpa=lrpa_save
    !
    ! ... each q /= gamma is saved on a different directory
    lgamma = lgamma_iq(iq)
    !
    tmp_dir_kcwq= TRIM (tmp_dir_kcw) //'q' &
                  & // TRIM(int_to_char(iq))//'/'
    filename=TRIM(tmp_dir_kcwq)//'charge-density.dat'
    IF (ionode) inquire (file =TRIM(filename), exist = exst)
    !
    CALL mp_bcast( exst, ionode_id, intra_image_comm )
    !
    CALL create_directory( tmp_dir_kcwq )
    tmp_dir=tmp_dir_kcwq
    CALL write_scf( rho, nspin )
    ! write the periodic part of the wannier orbital density on file
    DO i = 1, num_wann
      !
      IF ( .NOT. io_real_space) THEN 
        !
        rhog = rhowann_g(:,i)
        file_base=TRIM(tmp_dir_kcwq)//'rhowann_g_iwann_'//TRIM(int_to_char(i))
        IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
             CALL write_rhowann_g( file_base, &
             root_bgrp, intra_bgrp_comm, &
             bg(:,1)*tpiba, bg(:,2)*tpiba, bg(:,3)*tpiba, &
             gamma_only, mill, ig_l2g, rhog(:), .FALSE. )
        !
      ELSE  
        !
        rhowann_aux (:) = rhowann(:,iq, i)
        file_base=TRIM(tmp_dir_kcwq)//'rhowann_iwann_'//TRIM(int_to_char(i))
        CALL write_rhowann( file_base, rhowann_aux, dffts, ionode, inter_bgrp_comm )
        !
      ENDIF
    ENDDO
    tmp_dir=tmp_dir_save
    !
  ENDDO
  !
  ! Print on output the self-Hatree
  CALL mp_sum (sh, intra_bgrp_comm)
  !
  OPEN (128, file = TRIM(tmp_dir_kcw)//'sh.txt')
  WRITE(stdout,'(5X, "INFO: Orbital Self-Hartree (SH)")') 
  DO i = 1, num_wann
    WRITE(stdout,'(5X, "orb ", 1i5, 5X, "SH ", 1F10.6)') i, REAL(sh(i))
    WRITE(128,*) sh(i)
  ENDDO
  CLOSE (128)
  !
  WRITE( stdout, '(/,5X,"INFO: PREPARING THE KCW CALCULATION ... DONE")')
  WRITE(stdout,'(/)')
  !
  ! Find symmetries respected by each wannier function
  !
  IF(irr_bz) THEN
    !
    CALL symmetries_of_wannier_function(rhowann)
    !
    ! Define IBZ for every wannier function
    !
    CALL kcw_kpoint_grid()
    !
    ! Compute self - Hartree with the IBZ and the weights  
    !
    sh=0 
    DO iq = 1, nqs
      DO i = 1, num_wann
        !
        ! check if the q point is in the IBZ
        !
        IF( fbz2ibz(iq, i) .eq. -1 ) CYCLE
        iq_ibz = fbz2ibz(iq, i)
        ! 
        lrpa_save=lrpa
        lrpa = .true.
        !
        rhog(:)         = CMPLX(0.D0,0.D0,kind=DP)
        delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
        vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
        rhor(:)         = CMPLX(0.D0,0.D0,kind=DP)
        !
        rhor(:) = rhowann(:,iq, i) 
        !! The periodic part of the orbital desity in real space
        !
        CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
        rhowann_g(:,i) = rhog
        !! The periodic part of the perturbation DeltaV_q(G)
        ! 
        sh(i) = sh(i) + 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:) )*wq_ibz(iq_ibz, i)*omega
      END DO!iwann
    END DO!iq

    ! Print on output the self-Hatree
    CALL mp_sum (sh, intra_bgrp_comm)
    !
    WRITE(stdout,'(/, 5X, "INFO: Orbital Self-Hartree (SH) with Symmetries")') 
    DO i = 1, num_wann
      WRITE(stdout,'(5X, "orb ", 1i5, 5X, "SH ", 1F10.6)') i, REAL(sh(i))
    END DO
    !
    !create a file with all the info about IBZ for each wannier function
    !
    IF (ionode) THEN 
      CALL write_qlist_ibz()
      DO i = 1, num_wann
        CALL write_symmetry_op(i)
      END DO 
    ENDIF
      !
    !CALL compute_dmn()
  END IF!if for symmetries
  !
  CALL close_buffer  ( iuwfc, 'KEEP' )
  !
  CALL stop_clock ('kcw_setup')
  !
  DEALLOCATE (rhowann, rhowann_aux )
  DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
  !
  ! 
  RETURN
  !
END SUBROUTINE kcw_setup
