!------------------------------------------------------------------------------
!
! runctl.inc
! 
!------------------------------------------------------------------------------
!
! variable name   & description & default \\
!
!------------------------------------------------------------------------------
NAMELIST /parctl/      &
  nproca,              & ! number of processors for set A
  nprocb,              & ! number of processors for set B
  nprocio,             & ! number of processors for I/O server
  iomode,              & ! I/O mode for parallel I/O
  network_logger,      & ! hostname for network logging
  db_host,             & ! hostname of db server of timing database 
  lyaxt_transposition    ! if .true. use YAXT library for transpositions
!------------------------------------------------------------------------------
