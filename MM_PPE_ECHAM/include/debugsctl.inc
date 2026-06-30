! debugsctl.inc
!    sebastian.rast@zmaw.de, Hamburg October 20th, 2008
!    namelist for debug stream. Detailed documentation: debugs.pdf
!    -----------------------------------------------------------------------
NAMELIST /debugsctl/ &
  nddf,             &    ! number of additional 3-D fields on full levels
  nddfh,            &    ! number of additional 3-D fields on level interfaces 
                         ! (half levels)
  nzdf,             &    ! number of additional 2-D fields
  putdebug_stream        ! time interval of output of debug stream 
                         ! (default: 6-hourly output)
