!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
! ==================================================================================================================================
!
!                                               DRIVER OF THE STANDALONE-JSBACH
!
! ==================================================================================================================================

PROGRAM jsbalone_driver

  USE mo_kind,             ONLY : dp
  USE mo_jsbach,           ONLY : options_type, nml_unit, debug
  USE mo_mpi,              ONLY : p_io, p_parallel, p_start, p_stop, p_bcast, p_parallel_io, &
                                  p_global_comm, p_init_communicators
  USE mo_time_control,     ONLY : get_time_step, time_set, time_reset, l_putrerun, lstop, l_trigfiles, lstart
  USE mo_timer,            ONLY : init_timer, cleanup_timer
  USE mo_jsbalone_forcing, ONLY : update_forcing, read_forcing, finish_forcing, forcing_type
  USE mo_zenith,           ONLY : compute_orbit_and_solar, cos_zenith
  USE mo_output,           ONLY : open_output_streams, close_output_streams, out_streams
  USE mo_io,               ONLY : write_streams
  USE mo_jsbach_grid,      ONLY : grid_type, domain_type
  USE mo_exception,        ONLY : message, message_text
  USE mo_util_string,      ONLY : int2string
  USE mo_jsbach_interface, ONLY : jsbach_inter_1d, jsbach_init, stepon_jsbach, stepoff_jsbach, theOptions
  USE mo_jsbach_constants, ONLY : Tmelt
  USE mo_kind,             ONLY : dp
  USE mo_jsbalone,         ONLY : update_driving, drive_type
  USE mo_machine,          ONLY : machine_setup
  USE mo_namelist,         ONLY : POSITIONED, open_nml, position_nml
  USE mo_param_switches,   ONLY : lsurf
  USE mo_jsbach_version,   ONLY : jsbach_init_version, jsbach_label_run
  USE mo_input,            ONLY : InputClose
  USE mo_hydrology,        ONLY : cleanup_hydrology
  USE mo_test,             ONLY : read_interface_variables

  implicit none

  ! Local structures to get grid and domain description back from the interface (jsbach_init)
  TYPE(grid_type)   :: grid
  TYPE(domain_type) :: domain

  ! Local structure to get JSBACH options from interface (jsbach_init)
  TYPE(options_type) :: options

  ! Local structure to hold forcing fields
  TYPE(forcing_type) :: forcing

  ! Local structure for driving field (derived vars from forcing fields)
  TYPE(drive_type)      :: driving

  INTEGER, EXTERNAL :: util_cputime !  External functions
  REAL(dp), EXTERNAL:: util_walltime
  REAL(dp):: zutime, zstime, zrtime, zwtime

  INTEGER           :: istep,status
  INTEGER           :: read_status, f_unit
  INTEGER           :: ntile

  !! Variables of namelist jsbalone_parctl
  INTEGER :: nproca  ! number of processors for jsbach
  INTEGER :: nprocb  ! has to be one
  INTEGER :: nprocio ! number of processors for I/O server
  INTEGER :: npedim  ! Working dimension for blocks in each domain.
                     ! Default (-1): each domain is processed in one call

  INTEGER :: buffer(4) ! buffer to broadcast the 4 above variables

#ifdef STANDALONE
  EXTERNAL jsbalone_iniphy
#endif

  INCLUDE 'jsbalone_parctl.inc'

  !=================================================================================================
  ! Start MPI
  CALL p_start
  !
  !=================================================================================================
  CALL jsbach_init_version
  !
  !=================================================================================================
  ! Initialize wallclock timer
!$OMP PARALLEL
!$OMP MASTER
  zwtime = util_walltime(0)
!$OMP END MASTER
!$OMP END PARALLEL
  !=================================================================================================
  ! read processor decomposition

  IF (p_parallel_io) THEN

     ! define default values
     nproca = 1
     nprocb = 1
     nprocio = 0
     npedim = -1

     ! read the namelist
     nml_unit = open_nml ('namelist.jsbach')
     f_unit = position_nml ('JSBALONE_PARCTL', nml_unit, status=read_status)
     SELECT CASE (read_status)
     CASE (POSITIONED)
        READ (f_unit, jsbalone_parctl)
     END SELECT
  ENDIF

  IF (p_parallel) THEN
     buffer(:) = (/ nproca, nprocb, nprocio, npedim /)
     CALL p_bcast (buffer, p_io, comm=p_global_comm)
     nproca  = buffer(1)
     nprocb  = buffer(2)
     nprocio = buffer(3)
     npedim  = buffer(4)
  END IF
  CALL p_init_communicators(nproca, nprocb, nprocio)

  CALL machine_setup

  CALL init_timer

  CALL jsbach_init(grid, domain, options, ntile &
#ifdef STANDALONE
                        , nproca, nprocb, npedim, forcing, driving &
#endif
                            )
  IF (debug) CALL message('jsbalone_driver','jsbach_init done')
  CALL jsbach_label_run
#ifdef STANDALONE
  CALL jsbalone_iniphy
#endif

  lsurf = .true.

  ! Main loop

  IF (debug) CALL message('jsbalone_driver','beginning of time step loop')
  main_loop:DO

     ! evaluate events at begin of time step

     IF (debug) CALL message('jsbalone_driver','call time_set')
     CALL time_set

     IF (l_trigfiles) THEN
        IF (debug) CALL message('jsbalone_driver','l_trigfiles .true.')
        CALL close_output_streams
        CALL open_output_streams
     ENDIF

     IF (.NOT. lstart) THEN 
        IF (debug) CALL message('jsbalone_driver','call stepon_jsbach')
        CALL stepon_jsbach
     END IF

     ! compute declination and cosine of zenith angle
     IF (debug) CALL message('jsbalone_driver','call compute_orbit_and_solar')
     CALL compute_orbit_and_solar(domain)

     IF (theOptions%ReadInterfaceVars) THEN
        IF (debug) CALL message('jsbalone_driver','call read_interface_variables')
        CALL read_interface_variables(grid, domain, forcing, driving, cos_zenith)
     ELSE
        ! reads forcing date from file
        IF (debug) CALL message('jsbalone_driver','call read_forcing')
        CALL read_forcing (grid, domain, forcing)

        ! generate forcing data for particular time step 
        IF (debug) CALL message('jsbalone_driver','call update_forcing')
        CALL update_forcing(grid, domain, driving%evap_act2pot, forcing)

        IF (debug) CALL message('jsbalone_driver','call update_driving')
        CALL update_driving( domain%nland, &
                             domain%nland, &
                             driving, &
                             options%HeightWind, &
                             options%HeightHumidity, &
                             cos_zenith(1:domain%nland), &
                             forcing%wind_speed(1:domain%nland), &
                             forcing%air_temp(1:domain%nland) + Tmelt, &
                             forcing%spec_humidity(1:domain%nland), &
                             forcing%precip_rain(1:domain%nland), &
                             forcing%precip_snow(1:domain%nland), &
                             forcing%rad_lw_down(1:domain%nland), &
                             forcing%rad_UV_down(1:domain%nland), &
                             forcing%rad_PAR_down(1:domain%nland), &
                             forcing%rad_NIR_down(1:domain%nland), &
                             forcing%air_pressure(1:domain%nland), &
                             forcing%CO2_concentr(1:domain%nland))

        !vg: just to make it more obvious ...
        forcing%wind_speed10(1:domain%nland) = forcing%wind_speed(1:domain%nland)

     END IF

     ! call JSBACH-interface to land surface

     IF (debug) CALL message('jsbalone_driver','call jsbach_inter_1d')
     call jsbach_inter_1d(domain%nland, &
                          domain%nland, &
                          wind = forcing%wind_speed(1:domain%nland), &
                          wind10 = forcing%wind_speed10(1:domain%nland), &
                          temp_air = forcing%air_temp(1:domain%nland) + Tmelt, &
                          qair = forcing%spec_humidity(1:domain%nland), &
                          precip_rain = forcing%precip_rain(1:domain%nland), &
                          precip_snow = forcing%precip_snow(1:domain%nland), &
                          lwdown = forcing%rad_lw_down(1:domain%nland), &
                          sw_vis_net = driving%vis_net(1:domain%nland), &
                          sw_nir_net = driving%nir_net(1:domain%nland), &
                          sw_par_frac_diffuse = forcing%frac_PAR_diffuse(1:domain%nland), &
                          sw_par_down = forcing%rad_PAR_down(1:domain%nland), &
                          czenith = cos_zenith(1:domain%nland), &
                          pressure = forcing%air_pressure(1:domain%nland), &
                          CO2_concentration = forcing%CO2_concentr(1:domain%nland), &
                          evap_act = driving%evap_act(1:domain%nland), &                ! output
                          evap_pot = driving%evap_pot(1:domain%nland), &                ! output
                          temp_soil_new = driving%soil_temp(1:domain%nland), &          ! output
                          etacoef = driving%etacoef(1:domain%nland), &
                          etbcoef = driving%etbcoef(1:domain%nland), &
                          eqacoef = driving%eqacoef(1:domain%nland), &
                          eqbcoef = driving%eqbcoef(1:domain%nland), &
                          cair = driving%cair(1:domain%nland), &             ! output
                          csat = driving%csat(1:domain%nland), &             ! output
                          z0h = driving%z0h(1:domain%nland), &               ! output
                          z0m = driving%z0m(1:domain%nland), &               ! output
                          albedo_vis = driving%albedo_vis, &                 ! output
                          albedo_nir = driving%albedo_nir, &                 ! output
                          echam_zchl = driving%zchl(1:domain%nland),&
                          cdrag = driving%cdrag(1:domain%nland), &
                          veg_height = driving%veg_height(1:domain%nland))

     IF (debug) CALL message('jsbalone_driver','call out_streams')
     CALL out_streams

     IF (debug) CALL message('jsbalone_driver','call stepoff_jsbach')
     CALL stepoff_jsbach

     IF (l_putrerun) THEN
        IF (debug) CALL message('jsbalone_driver','call write_streams')
        CALL write_streams
     ENDIF

     IF (debug) CALL message('jsbalone_driver','call time_reset')
     CALL time_reset 

     ! check for last time step (then lstop=.true.) and eventually exit the main loop
     IF (lstop) THEN 
        istep = get_time_step() ! get the number of the time step
        CALL message('jsbach_driver','Step '//int2string(istep)//'completed.')
        EXIT main_loop   
     END IF

  END DO main_loop

  CALL InputClose()

  IF (lstop) THEN
     CALL close_output_streams
  END IF

  ! Cleanup timer
  CALL cleanup_timer

  ! Deallocate memory of forcing fields

  CALL finish_forcing(options, forcing)

  IF (theOptions%withHD) CALL cleanup_hydrology

!$OMP PARALLEL
!$OMP MASTER
  status = util_cputime(zutime, zstime)
!$OMP END MASTER
!$OMP END PARALLEL
  IF (status == -1) THEN
     CALL message('stepon','Cannot determine used CPU time')
  ELSE
!$OMP PARALLEL
!$OMP MASTER
     zwtime = util_walltime(0)
!$OMP END MASTER
!$OMP END PARALLEL
     zrtime = (zutime+zstime)/zwtime
     CALL message ('', '')
     WRITE (message_text,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
     CALL message('',message_text)
     WRITE (message_text,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
     CALL message('',message_text)
     WRITE (message_text,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
     CALL message('',message_text)
     WRITE (message_text,'(a,f10.2,a)') ' Ratio            : ', 100*zrtime, ' %'
     CALL message('',message_text)
     CALL message ('', '')
  END IF

  ! Stop MPI
  CALL p_stop

END PROGRAM jsbalone_driver
