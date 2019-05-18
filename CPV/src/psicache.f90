MODULE psicache
  USE kinds
  IMPLICIT NONE
  SAVE

  COMPLEX(DP), ALLOCATABLE, TARGET :: psi_save(:,:)
  LOGICAL :: use_psicache = .FALSE.
END MODULE
