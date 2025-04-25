!------------------------------------------------------------------------------
!
! jsbalone_parctl.inc
! 
!------------------------------------------------------------------------------
!
! variable name   & description & default \\
!
!------------------------------------------------------------------------------
NAMELIST /jsbalone_parctl/      &
  nproca,              & ! number of processors for set A
  nprocb,              & ! number of processors for set B
  nprocio,             & ! number of processors for I/O server
  npedim                 ! Working dimension for blocks in each domain.
                         ! Default (-1): each domain is processed in one call
!------------------------------------------------------------------------------
