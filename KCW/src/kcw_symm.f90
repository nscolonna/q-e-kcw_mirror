MODULE kcw_symm
  USE kinds,             ONLY : DP
  !
  INTEGER, ALLOCATABLE      :: nsym_w(:)        ! dimension: num_wann. nsym_w(i) is the number of symmetries
                                                ! respected by the i-th wannier function.
  REAL(DP), ALLOCATABLE     :: sr_w (:,:,:,:)   ! rotation in cart. coord (3,3,48,num_Wann)
  INTEGER(DP), ALLOCATABLE  :: s_w (:,:,:,:)    ! rotation in Crystal coord (3,3,48,num_Wann)
  REAL(DP), ALLOCATABLE     :: ft_w (:,:,:)     ! fractional translation in cryst. coord. (3,48,num_wann)
  REAL(DP), ALLOCATABLE     :: ftcart_w (:,:,:) ! fractional translation in cartesian coord. (3,48,num_wann)
  INTEGER, ALLOCATABLE      :: wsym2sym(:, :)   ! table that map the index of Wannier symmetry to the 
                                                ! global index (all the symmetries of the system) (48,num_wann) DEPRECATED?
  INTEGER, ALLOCATABLE      :: isym_w(:,:)        ! index of i-th simmetry for each wannier function 
END MODULE
