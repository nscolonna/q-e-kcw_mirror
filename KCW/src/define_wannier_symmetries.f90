SUBROUTINE define_wannier_symmetries ()
    !
    ! This subroutine checks which symmetries are respected by each wannier
    ! function that has already been stored by rotate_ks()
    ! We try to understand for all wannier functions if the following identity
    ! is respected
    !          evc0_{R.k}(r) ?= evc0_k(R^-1.r-f)
    !
    USE symm_base,            ONLY : s, sr, ft, nsym    
    USE cell_base,            ONLY : bg

    USE wvfct,                ONLY : npwx
    USE klist,                ONLY : nkstot, nks
    USE lsda_mod,             ONLY : nspin
    USE fft_base,             ONLY : dffts
    USE buffers,              ONLY : get_buffer
    USE fft_interfaces,       ONLY : invfft

    USE control_kcw,          ONLY : irr_bz, num_wann, evc0, iuwfc_wann, spin_component
    USE io_global,            ONLY : stdout
    USE kcw_symm
    !
    IMPLICIT NONE 
    !
    INTEGER,  EXTERNAL       :: find_free_unit
    INTEGER                  :: iun_lg, ipol, i, ieff
    CHARACTER (len=60)       :: header
    LOGICAL                  :: exst
    CHARACTER (LEN=256)      :: filename
    INTEGER                  :: ierr
    !
    !Counters for do loops
    INTEGER                  :: ik, isym, iwann
    INTEGER                  :: ik_eff 
    INTEGER                  :: ikrot
    COMPLEX(DP)              :: phase !Will contain the phase to bring the k point back to FBZ
    INTEGER                  :: lrwfc !length of buffer of wannier function
    !
    !
    COMPLEX(DP), ALLOCATABLE ::  Revc0(:,:)
    REAL(DP), ALLOCATABLE    ::  gvector(:)
    ! 
    !IF ( .NOT. irr_bz) RETURN 
    !
    WRITE(stdout, *) "Checking which symmetries are satisfied by each wannier function..."
    !
    !
    wsym2sym  = 0
    sr_w =0.D0
    ft_w = 0.D0
    ftcart_w = 0.00
    nsym_w = 0
    !
    ALLOCATE( Revc0(npwx, num_wann) )
    ALLOCATE( gvector(3) )
    !
    !
    !
    !
    DO iwann = 1, num_wann
      WRITE(stdout, *) "iwann  ", iwann, "/", num_wann
      DO ik = 1, nkstot
        !
        WRITE(stdout, *) "ik   ", ik, "/", nkstot
        !
        !read wf (already rotated to wannier gauge) at kpoint ik
        !
        lrwfc = num_wann * npwx 
        !ik_eff = ik-(spin_component-1)*nkstot/nspin
        WRITE(stdout, *) "Reading evc0 from", iuwfc_wann, "with record", ik
        CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik )
        WRITE(stdout, *) "Read evc0 from", iuwfc_wann, "with record", ik
        !
        !Go to r space from G space
        !
        CALL invfft ('Wave', evc0(:,iwann), dffts)
        !
        !check which symmetries leave wannier function ik unchanged
        !
        DO isym = 1, nsym
          !
          ! Rotate k point with isym to get irot and gvector
          !
          CALL find_index_1bz_smart_ibz(ik, 1, gvector, ikrot, isym) !second argument -> iq=1 means (0,0,0)          
          !
          ! Get wf at rotated point and go to real space
          ! 
          WRITE(stdout, *) "Reading Revc0 from", iuwfc_wann, "with record", ikrot
          CALL get_buffer ( Revc0, lrwfc, iuwfc_wann, ikrot )
          WRITE(stdout, *) "Read Revc0 from", iuwfc_wann, "with record", ikrot
          CALL invfft ('Wave', Revc0(:,iwann), dffts)
          !
          ! rotate wannier function iwann in real space
          !
          CALL rotate_evc(isym, gvector, evc0)
          !
          ! Check if evc_k(R^-1.r-f) = evc_(R.k)(r)
          !
          IF ( ANY( ABS( evc0-Revc0 ) > 1.d-08 ) ) THEN
            EXIT
          END IF 
          !
          ! Check if we went through all k points and in this case
          ! add the symmetry to the ones respected by the wannier
          ! function.
          !        
          IF ( ik == nkstot ) THEN
            nsym_w(iwann) = nsym_w(iwann) + 1
            isym_w(nsym_w(iwann), iwann) = isym 
          END IF
        END DO !isym
      END DO !ik
    END DO ! iwann 
END SUBROUTINE define_wannier_symmetries




