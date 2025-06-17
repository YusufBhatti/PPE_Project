!------------------------------------------------------------------------------
!
! runctl.inc
! 
!------------------------------------------------------------------------------
!
! variable name   & description \\
!
!------------------------------------------------------------------------------
NAMELIST /runctl/    &
  nproma,            &! working dimension for grid-point computations
  out_expname,       &! prefix of file names with experiment name
  out_datapath,      &! name of data path with trailing '/'
  out_filetype,      &! output file type: (GRIB=1, NetCDF=2, NetCDF4=6)         
  out_ztype,         &! output compression type: (SZIP=1, ZIP=2)
  default_output,    &! disable output for standard streams [true]
  rerun_filetype,    &! rerun file type: (NetCDF=2, NetCDF2=4)
  trac_filetype,     &! output file type: (GRIB=1, NetCDF=2)
  iadvec,            &! selects advection scheme: 
                      ! 0 - none
                      ! 1 - Semi Lagrangian Transport
                      ! 2 - SPITFIRE
                      ! 3 - FFSL
  lresume,           &! switch resume model from restart files
  ltimer,            &! switch use timer
  ltctest,           &! switch to test time control
  lnwp,              &! switch NWP mode on/off
  lmidatm,           &! switch middle atmosphere model version on/off
  lmlo,              &! switch mixed layer ocean on/off
  lmeltpond,         &! switch calculation of meltponds on/off
  lcouple,           &! switch ocean coupling on/off
  lcouple_co2,       &! switch to couple co2 in addition
  ldailysst,         &! switch daily SST and SIC on/off
  lamip,             &! switch AMIP sst on/off
  ldiagamip,         &! switch AMIP diagnostics on/off
  lport,             &! switch on random perturbation for port test
  ltdiag,            &! switch tendency diagnostic on/off -> mo_diag_tendency
  lnmi,              &! switch NMI on/off
  lnudge,            &! switch nudging on/off
  l_orbvsop87,       &! use orbit functions from VSOP87 (T) or PCMDI (AMIP) (F)
  ly360,             &! use 360 days year
  lfractional_mask,  &! use fractional land sea mask
  dt_start,          &! external defined start date
  dt_stop,           &! external defined stop date
  dt_resume,         &! external defined restart date
  delta_time,        &! external defined time step
  no_cycles,         &! stop after number of rerun cycles
  no_steps,          &! stop after number of steps
  no_days,           &! stop after number of days
  putocean,          &! receive ocean data coupling intervals
  getocean,          &! send atmospheric data coupling intervals
  putrerun,          &! store-rerun-files interval and stops current job
  putcheckpoint,     &! store-rerun-files interval 
  nsub,              &! number of subjobs
  trigfiles,         &! compose new filenames
  trigjob,           &! subjob submit interval
  putdata,           &! put model data to file
  subflag,           &! bind output stream subjobs
  l_volc,            &! volcanic aerosols on/off
  ldebug,            &! switch mass fixer diagnostics on/off
  ldebugev,          &! switch debugging events
  ldebugs,           &! switch for debug stream on/off
  ldebugio,          &! switch debug I/O
  ldebugmem,         &! switch debug memory
  lindependent_read, &! switch to read initial/restart input on each MPI rank separatly 
  lcollective_write, &! switch to write restart in parallel on each MPI rank
  ndiahdf,           &! I/O unit for hdiff diagnostics
  lcolumn,           &! Switch for Single Column Model
  earth_angular_velocity,&! value for solid earth angular velocity (default=7.29212e-5_wp)
  lrce,              &! if .true.: - turn on radiative-convective eq. radiation calculation,
                      !            - use constant ocean surface albedo (0.07)
                      !            - ignore dynamical pbl height in pbl calculation.
                      ! specify also: ldiur & isolrad (radctl.inc), apsurf (dynctl.inc),
                      !               earth_angular_velocity = 0.0
  lmlo_ice,          &! switch to turn on ice for mixed layer ocean
  rmlo_depth,        &! depth of the mixed layer ocean
  lprint_m0,         &! print cpu time for every timestep (default = F)
!------------------------------------------------------------------------------
!>>SF: switchable internal reruns
  lforcererun          ! flag to switch on/off internal reruns
!<<SF 

#ifdef __prism
  NAMELIST /runctl/  &
  lcouple_parallel,  &! switch for all/one PE communication with oasis
  ldebugcpl           ! switch debug coupling
#endif
!------------------------------------------------------------------------------
