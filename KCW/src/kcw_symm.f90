MODULE kcw_symm
  USE kinds,             ONLY : DP
  !
  INTEGER, ALLOCATABLE      :: nsym_w(:)
  REAL(DP), ALLOCATABLE     :: sr_w (:,:,:,:)   ! rotation in cart. coord (3,3,48,num_Wann)
  INTEGER(DP), ALLOCATABLE  :: s_w (:,:,:,:)    ! rotation in Crystal coord (3,3,48,num_Wann)
  REAL(DP), ALLOCATABLE     :: ft_w (:,:,:)     ! fractional translation in cryst. coord. (3,48,num_wann)
  REAL(DP), ALLOCATABLE     :: ftcart_w (:,:,:) ! fractional translation in cartesian coord. (3,48,num_wann)
  INTEGER, ALLOCATABLE      :: wsym2sym(:, :)   ! table that map the index of Wannier symmetry to the 
                                                ! global index (all the symmetry of the system) (48,num_wann)
END MODULE
