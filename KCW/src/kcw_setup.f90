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
  USE parameters,        ONLY : npk
  USE ions_base,         ONLY : nat, ityp
  USE io_files,          ONLY : tmp_dir
  USE lsda_mod,          ONLY : nspin, starting_magnetization, lsda, isk
  USE scf,               ONLY : v, vrs, vltot,  kedtau
  USE fft_base,          ONLY : dfftp, dffts
  USE gvecs,             ONLY : doublegrid, ngms
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux!, nspin_mag, npol
  USE wvfct,             ONLY : nbnd
  !USE funct,             ONLY : dft_is_gradient
  USE xc_lib,            ONLY : xclib_dft_is
  !
  USE units_lr,          ONLY : iuwfc
  USE wvfct,             ONLY : npwx
  USE control_flags,     ONLY : io_level
  USE io_files,          ONLY : prefix
  USE buffers,           ONLY : open_buffer, save_buffer, close_buffer, get_buffer
  USE control_kcw,       ONLY : evc0, iuwfc_wann, iuwfc_wann_allk, kcw_iverbosity, lgamma_iq, &
                                spin_component, isq, read_unitary_matrix, x_q, tmp_dir_save, & 
                                num_wann, num_wann_occ, occ_mat, tmp_dir_kcw, tmp_dir_kcwq, &
                                irr_bz, mp1, mp2, mp3, seedname!, wq, nqstot
  USE io_global,         ONLY : stdout
  USE klist,             ONLY : nkstot, xk, wk, nks
  USE cell_base,         ONLY : at, omega, bg
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
  USE mp_bands,         ONLY : inter_bgrp_comm, intra_bgrp_comm
  USE io_kcw,           ONLY : write_rhowann
  !
  USE mp,               ONLY : mp_sum
  USE mp_world,         ONLY :  mpime
  USE control_lr,       ONLY : lrpa
  !
  USE symm_base,        ONLY : s, t_rev, nrot, nsym, time_reversal, sr
  USE start_k,          ONLY : nk1, nk2, nk3, k1, k2, k3
  !
  USE lr_symm_base,     ONLY : nsymq
  USE qpoint,           ONLY : xq
  !
  implicit none

  integer :: na, i, ik, iwann, isym
  ! counters
  !
  INTEGER   :: lrwfc, iun_qlist
  LOGICAL   :: exst, exst_mem
  INTEGER :: iq, nqs, iq_eff
  REAL(DP) :: xq_(3)
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:), rhowann_aux(:), rhowann_k(:,:,:)
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
  COMPLEX(DP) :: rhor(dffts%nnr), delta_vr(dffts%nnr,nspin), delta_vr_(dffts%nnr,nspin)
  !
  ! The weight of each q point
  REAL(DP), ALLOCATABLE :: weight(:)
  LOGICAL :: lrpa_save
  !
  INTEGER  :: t_rev_eff(48)
  LOGICAL  :: skip_equivalence
  !
  integer :: isk_ibz(npk), nkstot_ibz
  real(DP) :: xk_ibz(3,npk), wk_ibz(npk)
  !
  INTEGER                  :: nsym_w(num_wann_occ), nsym_aux
  REAL(DP)                 :: sr_w (3,3,48,num_wann_occ) ! rotation in cart. coord
  REAL(DP)                 ::  s_w (3,3,48,num_wann_occ) ! rotation in Crystal coord
  REAL(DP)                 :: tvec_w (3,48,num_wann_occ) ! fractional translation in cartesian coord.
  INTEGER                  :: wsym2sym(48, num_wann_occ)
  INTEGER                  :: iun_lg
  CHARACTER (len=60)       :: header
  INTEGER, ALLOCATABLE     :: ik_ibz2ik(:)
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
  IF (irr_bz) THEN 
     wsym2sym  = 0
     sr_w =0.D0
     tvec_w = 0.D0
     nsym_w = 0
     !
     filename=TRIM(seedname)//'.lg'
     OPEN (101, file=filename, form='formatted')
     write(*,*) filename
     READ (101,*) header
     WRITE(*,*) header 
     READ (101,*)
     DO i = 1, num_wann_occ
        WRITE(*,*) "iwann =", i
        READ(101,*) iwann, nsym_w(i)
        WRITE(*,*) iwann, nsym_w(i)
        !nsym_w(i)=nsym_aux
        READ(101, "(2x, 10i4)") (wsym2sym(isym, i), isym=1,nsym_w(i))
        write(stdout,"(2x, 10i4)") wsym2sym(1:nsym_w(i), i)
        DO isym = 1, nsym_w(i)
           WRITE(*,*) isym, wsym2sym(isym,i)
           READ(101,"(3f15.7)") sr_w(:,:,isym, i), tvec_w(:,isym, i) !reading rotation matrix and translation vector in Cartesian coordinates.
           WRITE(*,"(3f15.7)")  sr_w(:,:,isym, i), tvec_w(:,isym, i) 
           s_w(:,:,isym,i) = s(:,:,wsym2sym(isym,i))
           sr_w(:,:,isym,i) = sr(:,:,wsym2sym(isym,i)) ! actually we know the map and we can select the rotation belonging to the
                                                       ! little group among all the symmetry of the system
           ! here NEED to add a check sr_w(from file) == sr_w (from index)
           WRITE(*,"(3f15.7)")  sr_w(:,:,isym, i), tvec_w(:,isym, i) 
           CALL cryst_to_cart(1, tvec_w(1:3,isym,i), at, -1)        ! fractional translation in crystal coord.
           ! here NEED to can add a check sr_w(from file) == sr_w (from index)
           WRITE(*,"(3f15.7)")  s_w(:,:,isym, i), tvec_w(:,isym, i) ! rotation and fractional translation in cryst. coord.
           !
           READ(101,*)
           WRITE(*,*)
        ENDDO
     ENDDO
     !
     t_rev_eff=0
     skip_equivalence=.false.
     nrot = nsym_w(1)
     s(:,:,1:nsym_w(1)) = s_w(:,:,1:nsym_w(1),1)
     CALL kcw_kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev_eff, &
                      bg, mp1*mp2*mp3, 0, 0, 0, mp1, mp2, mp3, nkstot_ibz, xk_ibz, wk_ibz)
     WRITE(*,'("NICOLA nkstot, nsym", 2I5)') nkstot_ibz, nrot
     CALL cryst_to_cart(nkstot_ibz, xk_ibz, at, -1)
     WRITE(*,'("NICOLA xk", 3F10.4)') (xk_ibz(1:3,ik), ik=1,nkstot_ibz) 
     CALL cryst_to_cart(nkstot_ibz, xk_ibz, bg, 1)
     IF (lsda) CALL set_kup_and_kdw( xk_ibz, wk_ibz, isk_ibz, nkstot_ibz, npk )
     ALLOCATE (ik_ibz2ik (nkstot_ibz) )
     CALL map_ikibz2ik (xk_ibz, nkstot_ibz, ik_ibz2ik, isk_ibz)
     WRITE(*,*) ik_ibz2ik
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
  ALLOCATE (rhowann ( dffts%nnr, num_wann), rhowann_aux(dffts%nnr) )
  ALLOCATE (rhowann_k ( dffts%nnr, nkstot/nspin, num_wann) )
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
  !DEALLOCATE ( nbnd_occ )  ! otherwise allocate_ph complains: FIXME
  !
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
  ALLOCATE (x_q (3, nqs), isq(nqs) )
  ALLOCATE (weight(nqs) )
  iq=1
  IF (ionode) THEN 
     WRITE(iun_qlist,'(i5)') nqs
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
  IF (irr_bz) nqs = nkstot_ibz/nspin
  DO iq = 1, nqs
    !! For each q in the mesh 
    !
    iq_eff = iq 
    IF (irr_bz) iq_eff = ik_ibz2ik(iq)
    xq_ = x_q(:,iq_eff)
    xq  = x_q(:,iq_eff)
    !
    ! IF (ionode) WRITE(iun_qlist,'(3f12.8)') xq_
    !
    lgamma_iq(iq)=(x_q(1,iq)==0.D0.AND.x_q(2,iq)==0.D0.AND.x_q(3,iq)==0.D0)
    CALL cryst_to_cart(1, xq_, at, -1)
    WRITE( stdout, '(/,8X, 78("="))')
    WRITE( stdout, '(  8X, "iq = ", 2i5)') iq, iq_eff
    WRITE( stdout, '(  8X, "The Wannier density at  q = ",3F12.7, "  [Cart ]")') xq(:)
    WRITE( stdout, '(  8X, "The Wannier density at  q = ",3F12.7, "  [Cryst]")') xq_(:)
    WRITE( stdout, '(  8X, 78("="),/)')
    !
    !IF (irr_bz) THEN
      !! This is to find the small group of q and the corresponding IBZ(q)
      !CALL setup_nscf ( .FALSE., xq, .TRUE. )
      !CALL cryst_to_cart(nkstot, xk, at, -1)
      !WRITE(*,'("NICOLA nkstot, nsymq", 2I5)') nkstot, nsymq 
      !WRITE(*,'("NICOLA xk", 3F10.4)') (xk(1:3,ik), ik=1,nkstot) 
      !CALL cryst_to_cart(nkstot, xk, bg, 1)
    !ENDIF
    !
    ! CALL compute_map_ikq_single (iq)
    ! The map to identify which k point in the 1BZ corresponds to k+q and the G vector that produce the mapping
    ! The results are stored in the global variable map_ikq and shift_1bz (used inside rho_of_q) 
    ! OBSOLETE: now moved inside rho_of_q. The search is done o the fly 
    !
    rhowann(:,:)=ZERO
    !! Initialize the periodic part of the wannier orbtal density at this q
    !
    CALL rho_of_q (rhowann, ngk_all, igk_k_all, iq_eff, rhowann_k)
    ! Compute the peridic part rho_q(r) of the wannier density rho(r)
    ! rho(r)   = \sum_q exp[iqr]rho_q(r)
    !
    WRITE( stdout, '(/, 8X,"INFO: rho_q(r) DONE ",/)')
    !
    ! Compute the Self Hartree
    weight(iq) = 1./nqs
    IF (irr_bz) weight(iq) = wk_ibz(iq) 
    WRITE(*,*) "NICOLA", weight(iq)
    lrpa_save=lrpa
    lrpa = .true.
    DO i = 1, num_wann
      !
      rhog(:)         = CMPLX(0.D0,0.D0,kind=DP)
      delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
      vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
      rhor(:)         = CMPLX(0.D0,0.D0,kind=DP)
      !
      rhor(:) = rhowann(:,i) 
      !! The periodic part of the orbital desity in real space
      !
      CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq_eff, delta_vr_, delta_vg_ )
      !! The periodic part of the perturbation DeltaV_q(G)
      DO ik = 1, nkstot/nspin
         !WRITE(*,'(3x, " ik = ", i5, " iwann = ", i5, " rhokq(1:3) =", 6F12.8)') ik, i,  rhowann_k(1:3,ik,i)
         WRITE(*,'(3x, " ik = ", i5, " iwann = ", i5, " <rhokq|DV> =", F12.8)') ik, i, &
                  0.5D0 * REAL(sum (CONJG(rhowann_k (:,ik,i)) * delta_vr(:,spin_component))/( dffts%nnr))*weight(iq)
      ENDDO
      WRITE(*,'(" iq = ", i5, " iwann = ", i5, " <rho|DV> =", F12.8)') iq, i, &
               0.5D0 * REAL(sum (CONJG(rhor (:)) * delta_vr(:,spin_component))/dffts%nnr)*weight(iq)
      WRITE(*,'(" iq = ", i5, " iwann = ", i5, " SH =", F12.8)') iq, i, &
               0.5D0 * REAL(sum (CONJG(rhog (:)) * vh_rhog(:))*weight(iq)*omega)
      ! 
      sh(i) = sh(i) + 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:) )*weight(iq)*omega
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
      rhowann_aux (:) = rhowann(:,i) 
      file_base=TRIM(tmp_dir_kcwq)//'rhowann_iwann_'//TRIM(int_to_char(i))
      CALL write_rhowann( file_base, rhowann_aux, dffts, ionode, inter_bgrp_comm )
    ENDDO
    tmp_dir=tmp_dir_save
    !
  ENDDO
  !
  ! Print on output the self-Hatree
  CALL mp_sum (sh, intra_bgrp_comm)
  
  WRITE(stdout,'(5X, "INFO: Orbital Self-Hartree (SH)")') 
  DO i = 1, num_wann
    WRITE(stdout,'(5X, "orb ", 1i5, 5X, "SH ", 1F10.6)') i, REAL(sh(i))
  ENDDO
  !
  WRITE( stdout, '(/,5X,"INFO: PREPARING THE KCW CALCULATION ... DONE")')
  WRITE(stdout,'(/)')
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
