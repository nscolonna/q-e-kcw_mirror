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
subroutine read_wannier ()
  !-----------------------------------------------------------------------
  !
  !! This routine read the relevant infos from a previous Wannier90 run
  !! The unitary matrix U are stored in unimatrx, unimatrx, unimatrx_opt.
  !! U(i,j) = <KS_i|wann_j> 
  !! Two cases are addressed:
  !! 1) A standard wannierization (unique manifold)
  !! 2) A separate wannierization where occ and emp states are treated independently. 
  !! In both cases I define a unique U whose dimension is num_wann x num_wann 
  !! and eventually a U_opt whose dimension in num_bands x num_wann
  !
  USE control_kcw,          ONLY : l_unique_manifold
  !
  IMPLICIT NONE 
  !
  IF (l_unique_manifold) THEN 
    !
    ! ... standard wannierization (occ + emp all together) 
    CALL read_wannier_unique_manifold ( )
    !
  ELSE
   !
   ! ... Separate wannierization (No occ-emp mixing)
   CALL read_wannier_two_manifold ( )
   CALL read_wannier_symmetry ()
   !
  ENDIF
  !
END subroutine read_wannier
  !
  !-----------------------------------------------------------------------
  subroutine read_wannier_unique_manifold ()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_kcw,          ONLY : unimatrx, seedname, has_disentangle, &
                                   unimatrx_opt, num_wann, kcw_iverbosity, &
                                   mp1, mp2 ,mp3, xk_fbz
  USE mp_global,            ONLY : intra_image_comm
  USE mp,                   ONLY : mp_bcast
  USE io_global,            ONLY : ionode, ionode_id
  USE cell_base,            ONLY : bg
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : nbnd
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE 
  !
  ! Local Variable
  !
  INTEGER :: ik
  ! ... the kpoint index
  !
  INTEGER :: num_wann_file, num_kpoint_file, num_ks_file
  ! ... the manifold space
  !
  REAL (DP) check
  ! ... The unitary matrix elemnts 
  !
  CHARACTER (len=256) :: u_file, u_dis_file, dum
  ! ... the name of the file containing the U matrix
  !
  INTEGER i, j, nkstot_fbz
  ! 
  INTEGER ierr
  !
  nkstot_fbz = mp1*mp2*mp3
  !! Here and in W90 We treat always one spin component at the time.
  !! While in PWscf nkstot contain k point for spin up (the first nkstot/nspin 
  !! and spin down (the last nkstot.nspin). See also rotate_ks
  !
  ALLOCATE ( xk_fbz (3,nkstot_fbz) )
  num_wann = 0
  !
  u_file=TRIM( seedname ) // '_u.mat'
  !
  ! ... The U matrix for empty states 
  !
  IF (has_disentangle) THEN
     !
     u_dis_file=TRIM( seedname ) // '_u_dis.mat'
     !! The Optimal subspace Matrix 
     !
     IF (ionode) THEN
        !
        OPEN (UNIT = 1002, FILE = u_dis_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
        ! 
        IF (ierr /= 0 ) call errore('rotate_orbitals', 'Error while reading Optimal unitary matrix', abs (ierr) )
        !
        READ (1002,*) dum
        READ (1002,*) num_kpoint_file, num_wann_file, num_ks_file
        !
        IF (num_kpoint_file /= nkstot_fbz) &
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U Optimal matrix', nkstot_fbz)
        !
        IF (num_wann_file /= num_wann .AND. num_wann /= 0) &
              CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U Optimal matrix', 1)
        !
        IF (num_ks_file /=  nbnd ) &
             CALL errore ('read_wannier', 'Mismatch between num KS state from PW and Wann90', 1)
        !
     ENDIF
     !
  ENDIF
  !
  IF (ionode) THEN 
     !
     OPEN (UNIT = 1001, FILE =u_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
     ! 
     IF (ierr /= 0 ) call errore('rotate_orbitals', 'reading Empty states unitary matrix', abs (ierr) )
     !
     READ (1001,*) dum 
     READ (1001,*) num_kpoint_file, num_wann_file, num_wann_file
     !
     IF (num_kpoint_file /= nkstot_fbz) &
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U Empty matrix', nkstot_fbz)
        !
     IF (num_wann_file /= num_wann .AND. num_wann /= 0) &
           CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U Empty matrix', 1)
     ! 
     num_wann = num_wann_file
     !! Store the numebr of wannier in a global variable
     !
  ENDIF
  !
  CALL mp_bcast( num_wann, ionode_id, intra_image_comm )
  CALL mp_bcast( num_ks_file, ionode_id, intra_image_comm )
  !
  ALLOCATE ( unimatrx( num_wann, num_wann, nkstot_fbz) )
  IF ( .NOT. has_disentangle) num_ks_file = num_wann
  ALLOCATE (unimatrx_opt(num_ks_file,num_wann,nkstot_fbz))
  !
  unimatrx = CMPLX(0.D0, 0.D0, kind=DP)
  unimatrx_opt = CMPLX(0.D0, 0.D0, kind=DP)
  !
  IF (ionode) THEN
    !
    IF (has_disentangle) THEN 
      !
      DO ik = 1, nkstot_fbz
        ! 
        READ (1002, *) xk_fbz(1,ik), xk_fbz(2,ik), xk_fbz(3,ik)
        READ (1002,'(f15.10,sp,f15.10)') ((unimatrx_opt(i,j,ik),i=1,num_ks_file),j=1,num_wann)
        !
      ENDDO 
      !
      ! ... transform the kpoints read from U file in cartesian coordinates 
      CALL cryst_to_cart(nkstot_fbz, xk_fbz, bg, 1)
      !
     ! check = 0.D0
     ! DO ik=1, nkstot_fbz
     !   check=check+(sum(xk_fbz(:,ik)-xk(:,ik)))
     ! ENDDO
     ! IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
      !
    ELSE
      !
      ! ... If no disentangle than the U-opt is a square matrix of dimension num_wan x num_wann 
      ! ... set it to the identity (in apply_u_mat we anyway apply both Ui_opt and U)
      !
      DO ik = 1, nkstot_fbz
        DO i=1, num_wann; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind = DP); ENDDO
      ENDDO
    ENDIF
    !
  ENDIF
  !
  CALL mp_bcast( unimatrx_opt, ionode_id, intra_image_comm )
  if ( kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Optimal Matrix READ")') 
  !
  IF (ionode) THEN   
    !
    DO ik = 1, nkstot_fbz
      !
      READ (1001, *) xk_fbz(1,ik), xk_fbz(2,ik), xk_fbz(3,ik)
      READ (1001,'(f15.10,sp,f15.10)') ((unimatrx(i,j,ik),i=1,num_wann),j=1,num_wann)
      !
    ENDDO
    !
    ! ... transform the kpoints read from U file in cartesian coordinates 
    CALL cryst_to_cart(nkstot_fbz, xk_fbz, bg, 1)
    !
    !check = 0.D0
    !DO ik=1, nkstot_eff
    !  check=check+(sum(xk_fbz(:,ik)-xk(:,ik)))
    !ENDDO
    !IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
    !
  ENDIF
  !
  CALL mp_bcast( unimatrx, ionode_id, intra_image_comm )
  !
  CLOSE (1001)
  CLOSE (1002)
  !
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: total number of Wannier functions", i5)') num_wann
  CALL mp_bcast( xk_fbz, ionode_id, intra_image_comm )
  !
#if defined (DEBUG)
  ! WRITE
  DO ik =1, nkstot_fbz
    WRITE (*, *) xk_fbz(1,ik), xk_fbz(2,ik), xk_fbz(3,ik)
    DO i=1,num_ks_file
       WRITE (*,'(10(f8.4,sp,f8.4))') (unimatrx_opt(i,j,ik),j=1,num_wann)
    ENDDO
    WRITE(*,*)
  ENDDO
  !
  DO ik =1, nkstot_fbz
    WRITE (*, *) xk_fbz(1,ik), xk_fbz(2,ik), xk_fbz(3,ik)
    DO i = 1, num_wann
       WRITE (*,'(10(f8.4,sp,f8.4))') (unimatrx(i,j,ik),j=1,num_wann)
    ENDDO
    WRITE(*,*)
  ENDDO
  !
#endif
  !
  RETURN
  !
END subroutine read_wannier_unique_manifold 
  !
  !
  !-----------------------------------------------------------------------
  subroutine read_wannier_two_manifold ()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_kcw,          ONLY : unimatrx, have_empty, num_wann_occ, seedname, &
                                   has_disentangle, have_empty, num_wann_emp, & 
                                   unimatrx_opt, num_wann, kcw_iverbosity, &
                                   mp1, mp2, mp3, xk_fbz
  USE control_lr,           ONLY : nbnd_occ
  USE mp_global,            ONLY : intra_image_comm
  USE mp,                   ONLY : mp_bcast
  USE io_global,            ONLY : ionode, ionode_id
  USE cell_base,            ONLY : bg
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : nbnd
  USE io_global,            ONLY : stdout
  !
  ! Local Variable
  !
  IMPLICIT NONE
  !
  INTEGER :: ik 
  ! ... the kpoint index
  !
  INTEGER :: num_wann_file, num_kpoint_file, num_ks_file
  ! ... the manifold space
  !
  REAL (DP) U_re, U_im, check
  ! ... The unitary matrix elemnts 
  !
  CHARACTER (len=256) :: u_file, u_dis_file, dum
  ! ... the name of the file containing the U matrix
  !
  INTEGER i, j, norb_occ, nkstot_fbz, norb_emp
  ! 
  INTEGER ierr
  !
  COMPLEX(DP), ALLOCATABLE :: unimatrx_occ(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: unimatrx_emp(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: unimatrx_emp_opt(:,:,:)
  INTEGER ieff, jeff
  !
  nkstot_fbz = mp1*mp2*mp3
  !! Here and in W90 We treat always one spin component at the time.
  !! While in PWscf nkstot contain k point for spin up (the first nkstot/nspin 
  !! and spin down (the last nkstot.nspin). See also rotate_ks
  !
  ALLOCATE ( xk_fbz (3, nkstot_fbz) )
  u_file=TRIM( seedname ) // '_u.mat'
  !
  IF (ionode) THEN
     !
     OPEN (UNIT = 1001, FILE = u_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
     ! 
     IF (ierr /= 0 ) call errore('rotate_orbitals', 'Error while reading unitary matrix', abs (ierr) )
     !
     READ (1001,*) dum
     READ (1001,*) num_kpoint_file, num_wann_file, num_wann_file
     !
     IF (num_kpoint_file /= nkstot_fbz) & 
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U matrix', nkstot_fbz)
     !
     IF (num_wann_file /= num_wann_occ .AND. num_wann_occ /= 0) &
              CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U matrix', 1)
     !
     num_wann_occ = num_wann_file
     !! Store the number of occupied wannier in a  global variable
     !
     IF (num_wann_occ /= nbnd_occ(1)) & 
          CALL errore ('read_wannier', 'Mismatch between  num occ bands and num wann', 1)
     !
  ENDIF
  !
  CALL mp_bcast( num_wann_occ, ionode_id, intra_image_comm )
  norb_occ = num_wann_occ ! Store the total number of occupied states
  ALLOCATE (unimatrx_occ(num_wann_occ, num_wann_occ, nkstot_fbz))
  !  
  IF ( ionode ) THEN
    !
    check = 0.D0 
    DO ik = 1, nkstot_fbz
    !
    !#### In the case occupied and empty manifolds are treated separately
    !#### we assume an insulating systems and NO disentanglement for occupied manifold
      ! 
      READ (1001, *) xk_fbz(1,ik), xk_fbz(2,ik), xk_fbz(3,ik)
      !
      DO i = 1, num_wann_occ
        DO j = 1, num_wann_occ
          READ (1001,*) U_re, U_im
          unimatrx_occ(j,i,ik)=CMPLX(U_re, U_im, kind=DP)
        ENDDO
      ENDDO
      !
    ENDDO 
    !
    ! ... transform the kpoints read from U file in cartesian coordinates 
    CALL cryst_to_cart(nkstot_fbz, xk_fbz, bg, 1)
    !
    !check = 0.D0
    !DO ik=1, nkstot_fbz
    !  check=check+(sum(xk_fbz(:,ik)-xk(:,ik)))
    !ENDDO 
    !IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
    !
  ENDIF
  !
  CALL mp_bcast( unimatrx_occ, ionode_id, intra_image_comm )
  CALL mp_bcast( xk_fbz, ionode_id, intra_image_comm )
  !
  CLOSE (1001) 
  !
  num_wann = num_wann_occ
  if (.NOT. have_empty .AND. kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: total number of Wannier functions", i5)') num_wann
  !
  IF (.NOT. have_empty) THEN
    !
    ! Store the unitary matrix in a global vabiable
    ALLOCATE (unimatrx (num_wann_occ, num_wann_occ, nkstot_fbz)) 
    unimatrx = unimatrx_occ
    !
    ! The optimal subspace Matrix
    ALLOCATE (unimatrx_opt (num_wann_occ, num_wann_occ, nkstot_fbz)) 
    unimatrx_opt=CMPLX(0.D0,0.D0, kind=DP)
    !
    ! The optimal matrix is simply the identity
    DO i = 1, num_wann_occ; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind=DP); ENDDO
    DEALLOCATE (unimatrx_occ)
    !
    RETURN  !nothing else to do
    !
  ENDIF
  !
  !####################
  !   EMPT YSTATES 
  !###################
  !
  ! ... read unitary matrix empty states ...
  !
  u_file=TRIM( seedname ) // '_emp_u.mat'
  !
  ! ... The U matrix for empty states 
  IF (has_disentangle) THEN
     !
     u_dis_file=TRIM( seedname ) // '_emp_u_dis.mat'
     !! The Optimal subspace Matrix 
     !
     IF (ionode) THEN
        !
        OPEN (UNIT = 1002, FILE = u_dis_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
        ! 
        IF (ierr /= 0 ) call errore('rotate_orbitals', 'Error while reading Optimal unitary matrix', abs (ierr) )
        !
        READ (1002,*) dum
        READ (1002,*) num_kpoint_file, num_wann_file, num_ks_file
        !
        IF (num_kpoint_file /= nkstot_fbz) &
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U Optimal matrix', nkstot_fbz)
        !
        IF (num_wann_file /= num_wann_emp .AND. num_wann_emp /= 0) &
              CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U Optimal matrix', 1)
        !
        IF (num_ks_file /= ( nbnd - num_wann_occ ) ) &
             CALL errore ('read_wannier', 'Mismatch between num empty bands and num empty wann', 1)
        !
     ENDIF
     !
  ENDIF
  !
  IF (ionode) THEN 
     !
     OPEN (UNIT = 1001, FILE =u_file, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
     ! 
     IF (ierr /= 0 ) call errore('rotate_orbitals', 'reading Empty states unitary matrix', abs (ierr) )
     !
     READ (1001,*) dum 
     READ (1001,*) num_kpoint_file, num_wann_file, num_wann_file
     !
     IF (num_kpoint_file /= nkstot_fbz) &
              CALL errore('read_wannier', 'Mismatch in num_kpoints input vs U Empty matrix', nkstot_fbz)
        !
     IF (num_wann_file /= num_wann_emp .AND. num_wann_emp /= 0) &
           CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U Empty matrix', 1)
     ! 
     num_wann_emp = num_wann_file
     !! Store the numebr of empty wannier in a global variable
     !
  ENDIF
  !
  CALL mp_bcast( num_wann_emp, ionode_id, intra_image_comm )
  CALL mp_bcast( num_ks_file, ionode_id, intra_image_comm )
  norb_emp = num_wann_emp ! store the total numebr of empty variational orbitals
  !
  ALLOCATE ( unimatrx_emp( num_wann_emp, num_wann_emp, nkstot_fbz) )
  IF (has_disentangle) ALLOCATE (unimatrx_emp_opt(num_ks_file,num_wann_emp,nkstot_fbz))
  !
  IF (ionode) THEN
    !
    IF (has_disentangle) THEN 
      !
      DO ik = 1, nkstot_fbz
        ! 
        READ (1002, *) xk_fbz(1,ik), xk_fbz(2,ik), xk_fbz(3,ik)
        READ (1002,'(f15.10,sp,f15.10)') ((unimatrx_emp_opt(i,j,ik),i=1,num_ks_file),j=1,num_wann_emp)
        !
      ENDDO 
      !
      ! ... transform the kpoints read from U file in cartesian coordinates 
      CALL cryst_to_cart(nkstot_fbz, xk_fbz, bg, 1)
      !
      !check = 0.D0
      !DO ik=1, nkstot_fbz
      !  check=check+(sum(xk_fbz(:,ik)-xk(:,ik)))
      !ENDDO
      !IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
      !
    ENDIF
    !
  ENDIF
  !
  CALL mp_bcast( unimatrx_emp_opt, ionode_id, intra_image_comm )
  if ( kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Optimal Matrix READ")') 
  !
  IF (ionode) THEN   
    !
    DO ik = 1, nkstot_fbz
      !
      READ (1001, *) xk_fbz(1,ik), xk_fbz(2,ik), xk_fbz(3,ik)
      READ (1001,'(f15.10,sp,f15.10)') ((unimatrx_emp(i,j,ik),i=1,num_wann_emp),j=1,num_wann_emp)
      !
    ENDDO
    !
    ! ... transform the kpoints read from U file in cartesian coordinates 
    CALL cryst_to_cart(nkstot_fbz, xk_fbz, bg, 1)
    !
    !check = 0.D0
    !DO ik=1, nkstot_eff
    !  check=check+(sum(xk_fbz(:,ik)-xk(:,ik)))
    !ENDDO
    !IF (check .gt. 1.D-06) CALL errore('read_wannier','Mismatch between Kpoints',1)
    !
  ENDIF
  !
  CALL mp_bcast( unimatrx_emp, ionode_id, intra_image_comm )
  !
  CLOSE (1001)
  CLOSE (1002)
  !
  num_wann = num_wann_occ + num_wann_emp
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: total number of Wannier functions", i5)') num_wann
  !
  ! Store the result in a unique matrix
  ALLOCATE (unimatrx (num_wann, num_wann, nkstot_fbz))
  ALLOCATE (unimatrx_opt (nbnd, num_wann, nkstot_fbz))
  !
  unimatrx=CMPLX(0.D0,0.D0, kind=DP)
  unimatrx_opt=CMPLX(0.D0,0.D0, kind=DP)
  !
  ! Build the unique disentangle matrix (identity for the occupied state and 
  ! unimatrx_emp_opt for empty states
  !
  ! 1)occ
  DO i = 1, num_wann_occ; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind=DP); ENDDO
  !
  ! 2)emp
  IF (.not. has_disentangle) THEN  
    !
    ! ...In this case the optimal matrix is the ideintity (nbnd=num_wann)
    DO i = 1, num_wann; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind=DP); ENDDO
    !
  ELSE
    !
    DO i = 1,num_ks_file
      !
      ieff=i+num_wann_occ
      DO j = 1, num_wann_emp
        jeff=j+num_wann_occ
        unimatrx_opt(ieff,jeff,:) = unimatrx_emp_opt(i,j,:)
      ENDDO
      !
    ENDDO
    !
  ENDIF 
  !
  ! Build the unique unitary matrix
  !
  ! 1)occ
  unimatrx(1:num_wann_occ, 1:num_wann_occ, :) = unimatrx_occ(1:num_wann_occ, 1:num_wann_occ, :)
  !
  ! 2)emp
  DO i = 1, num_wann_emp
    !
    ieff=i+num_wann_occ
    DO j = 1, num_wann_emp
      ! 
      jeff=j+num_wann_occ
      unimatrx(ieff,jeff,:) = unimatrx_emp(i,j,:)
      !
    ENDDO
    !
  ENDDO
  !
  DEALLOCATE (unimatrx_emp)
  IF (has_disentangle) DEALLOCATE (unimatrx_emp_opt )
  !

#if defined (DEBUG)
  ! WRITE
  DO ik =1, nkstot_fbz
    WRITE (*, *) xk_fbz(1,ik), xk_fbz(2,ik), xk_fbz(3,ik)
    DO i=1,nbnd
       WRITE (*,'(10(f8.4,sp,f8.4))') (unimatrx_opt(j,i,ik),j=1,num_wann)
    ENDDO
    WRITE(*,*)
  ENDDO
  !
  DO ik =1, nkstot_fbz
    WRITE (*, *) xk_fbz(1,ik), xk_fbz(2,ik), xk_fbz(3,ik)
    DO i = 1, num_wann
       WRITE (*,'(10(f8.4,sp,f8.4))') (unimatrx(j,i,ik),j=1,num_wann)
    ENDDO
    WRITE(*,*)
  ENDDO
  !
#endif
  !
  RETURN
  !
END subroutine read_wannier_two_manifold
!
!
!-----------------------------------------------------------------------
SUBROUTINE read_wannier_symmetry ()
  !
  USE symm_base,            ONLY : s, sr, ft
  USE cell_base,            ONLY : bg
  USE control_kcw,          ONLY : irr_bz, num_wann, num_wann_occ, &
                                   num_wann_emp, seedname, have_empty
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
  INTEGER                  :: isym, iwann, ierr
  !
  IF ( .NOT. irr_bz) RETURN 
  !
  ALLOCATE ( nsym_w(num_wann) )
  ALLOCATE ( sr_w(3,3,48,num_wann) )
  ALLOCATE ( s_w(3,3,48,num_wann) )
  ALLOCATE ( ft_w(3,48,num_wann) )
  ALLOCATE ( ftcart_w(3,48,num_wann) )
  ALLOCATE ( wsym2sym(48,num_wann) )
  !
  wsym2sym  = 0
  sr_w =0.D0
  ft_w = 0.D0
  ftcart_w = 0.00
  nsym_w = 0
  !
  ! #################
  ! # OCC. MANIFOLD #
  ! #################
  filename=TRIM(seedname)//'.lg'
  INQUIRE( file=filename, exist=exst )
  IF ( .NOT. exst) THEN
     CALL infomsg('kcw_setup','WARNING: file with Wannier symmetry subgroup occ. manifold &
             NOT FOUND. Going to Use the full BZ')
     irr_bz = .false.
     GOTO 101
  ENDIF
  !
  iun_lg=find_free_unit()
  OPEN (UNIT = iun_lg, FILE = filename, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
  IF (ierr /= 0 ) call errore('kcw_setup', 'Error while reading Wannier symmetry subgroup', abs (ierr) )
  WRITE(stdout,'(/,5x, a,2x, a)') "Reading Wanniey symmetry subgroup from:", TRIM(filename)
  READ (iun_lg,*) header
  READ (iun_lg,*)
  !
  DO i = 1, num_wann_occ
     READ(iun_lg,*) iwann, nsym_w(i)
     write(stdout,'(/, 5x, a, 2x, i5, 2x, a, i4)') "Symmetry subgroup of iwann =", i
     write(stdout,"(5x,a,i5,/)") "Number of symmetry operations = ", nsym_w(i)
     READ(iun_lg, "(2x, 10i4)") (wsym2sym(isym, i), isym=1,nsym_w(i))
     !write(stdout,"(5x, a)") "Map isym_wann --> isym"
     !write(stdout,"(5x, 10i4)")  wsym2sym(1:nsym_w(i), i)
     DO isym = 1, nsym_w(i)
        READ(iun_lg,"(3f15.7)") sr_w(:,:,isym, i), ftcart_w(:,isym, i) !reading rotation matrix and translation vector in Cartesian coordinates.
        ft_w(:,isym,i) = ft(:,wsym2sym(isym,i))
        ftcart_w(:,:,:) = ft_w(:,:,:)
        CALL cryst_to_cart(1, ftcart_w(1:3,isym,i), bg, 1)        ! fractional translation in crystal coord.
        s_w(:,:,isym,i) = s(:,:,wsym2sym(isym,i))
        sr_w(:,:,isym,i) = sr(:,:,wsym2sym(isym,i)) ! actually we know the map and we can select the rotation belonging to the
                                                    ! little group among all the symmetry of the system
        ! here NEED to add a check sr_w(from file) == sr_w (from index)
        ! Cartisian coordinates
        write(stdout,"(2x,a, i4, a, i4, a/)") "isym_w = ", isym, " ( --> isym =", wsym2sym(isym, i), ")"
        WRITE( stdout, '(1x,"cart. ",3x,"s_w(",i2,") = (",3f11.7, &
             &        " )    f =( ",f10.7," )")') &
             isym, (sr_w(1,ipol,isym,i),ipol=1,3), ftcart_w(1,isym,i)
        WRITE( stdout, '(19x," (",3f11.7, " )       ( ",f10.7," )")') &
             (sr_w(2,ipol,isym,i),ipol=1,3), ftcart_w(2,isym,i)
        WRITE( stdout, '(19x," (",3f11.7, " )       ( ",f10.7," )"/)') &
             (sr_w(3,ipol,isym,i),ipol=1,3), ftcart_w(3,isym,i)
        !
        ! Crystal coordinates
        WRITE( stdout, '(1x,"cryst.",3x,"s_w(",i2,") = (",3(i6,5x), &
             &        " )    f =( ",f10.7," )")') &
             isym, (s_w(1,ipol,isym,i),ipol=1,3), ft_w(1,isym,i)
        WRITE( stdout, '(19x," (",3(i6,5x), " )       ( ",f10.7," )")') &
             (s_w(2,ipol,isym,i),ipol=1,3), ft_w(2,isym,i)
        WRITE( stdout, '(19x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
             (s_w(3,ipol,isym,i),ipol=1,3), ft_w(3,isym,i)
        !
        READ(iun_lg,*)
        !WRITE(*,*)
     ENDDO
  ENDDO
  CLOSE (iun_lg)
101 CONTINUE
  !
  IF (.NOT. have_empty) RETURN ! nothing else to do
  !
  filename=TRIM(seedname)//'_emp.lg'
  INQUIRE( file=filename, exist=exst )
  IF ( .NOT. exst) THEN
     CALL infomsg('kcw_setup','WARNING: file with Wannier symmetry subgroup emp. manifold &
             NOT FOUND. Going to Use the full BZ')
     irr_bz = .false.
     GOTO 102
  ENDIF
  !
  iun_lg=find_free_unit()
  OPEN (UNIT = iun_lg, FILE = filename, FORM = 'formatted', STATUS = 'old', IOSTAT=ierr )
  IF (ierr /= 0 ) call errore('read_wannier_symmetry', 'Error while reading Wannier symmetry subgroup', abs (ierr) )
  WRITE(stdout,'(/,5x, a,2x, a)') "Reading Wanniey symmetry subgroup from:", TRIM(filename)
  READ (iun_lg,*) header
  READ (iun_lg,*)
  !
  DO i = 1, num_wann_emp
     ieff=num_wann_occ+i
     READ(iun_lg,*) iwann, nsym_w(ieff)
     write(stdout,'(/, 5x, a, 2x, i5, 2x, a, i4)') "Symmetry subgroup of iwann =", ieff
     write(stdout,"(5x,a,i5,/)") "Number of symmetry operations = ", nsym_w(ieff)
     READ(iun_lg, "(2x, 10i4)") (wsym2sym(isym, ieff), isym=1,nsym_w(ieff))
     !write(stdout,"(5x, a)") "Map isym_wann --> isym"
     !write(stdout,"(5x, 10i4)")  wsym2sym(1:nsym_w(ieff), ieff)
     DO isym = 1, nsym_w(ieff)
        READ(iun_lg,"(3f15.7)") sr_w(:,:,isym, ieff), ftcart_w(:,isym, ieff) !reading rotation matrix and translation vector in Cartesian coordinates.
        ft_w(:,isym,ieff) = ft(:,wsym2sym(isym,ieff))
        ftcart_w(:,:,:) = ft_w(:,:,:)
        CALL cryst_to_cart(1, ftcart_w(1:3,isym,ieff), bg, 1)        ! fractional translation in crystal coord.
        s_w(:,:,isym,ieff) = s(:,:,wsym2sym(isym,ieff))
        sr_w(:,:,isym,ieff) = sr(:,:,wsym2sym(isym,ieff)) ! actually we know the map and we can select the rotation belonging to the
                                                    ! little group among all the symmetry of the system
        ! here NEED to add a check sr_w(from file) == sr_w (from index)
        ! Cartisian coordinates
        write(stdout,"(2x,a, i4, a, i4, a/)") "isym_w = ", isym, " ( --> isym =", wsym2sym(isym, ieff), ")"
        WRITE( stdout, '(1x,"cart. ",3x,"s_w(",i2,") = (",3f11.7, &
             &        " )    f =( ",f10.7," )")') &
             isym, (sr_w(1,ipol,isym,ieff),ipol=1,3), ftcart_w(1,isym,ieff)
        WRITE( stdout, '(19x," (",3f11.7, " )       ( ",f10.7," )")') &
             (sr_w(2,ipol,isym,ieff),ipol=1,3), ftcart_w(2,isym,ieff)
        WRITE( stdout, '(19x," (",3f11.7, " )       ( ",f10.7," )"/)') &
             (sr_w(3,ipol,isym,ieff),ipol=1,3), ftcart_w(3,isym,ieff)
        !
        ! Crystal coordinates
        WRITE( stdout, '(1x,"cryst.",3x,"s_w(",i2,") = (",3(i6,5x), &
             &        " )    f =( ",f10.7," )")') &
             isym, (s_w(1,ipol,isym,ieff),ipol=1,3), ft_w(1,isym,ieff)
        WRITE( stdout, '(19x," (",3(i6,5x), " )       ( ",f10.7," )")') &
             (s_w(2,ipol,isym,ieff),ipol=1,3), ft_w(2,isym,ieff)
        WRITE( stdout, '(19x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
             (s_w(3,ipol,isym,ieff),ipol=1,3), ft_w(3,isym,ieff)
        !
        READ(iun_lg,*)
        !WRITE(*,*)
     ENDDO
  ENDDO
  CLOSE (iun_lg)
102 CONTINUE
  !
  RETURN
  !
END SUBROUTINE
