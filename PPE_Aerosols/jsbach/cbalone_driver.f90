!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                            !!
!!                          D R I V E R                       !!
!!                                                            !!
!!   for running the carbon and nitrogen pool model offline   !!
!!                                                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CHR: 09-06-04: adding  Nitrogen
!! Tea Thum, FMI     Oct 2012: first implementation of YASSO soil carbon model (J. Liski, EFI)
!! Daniel Goll, MPI  May 2013: implementation of YASSO soil carbon model

program cbalone_driver

  USE mo_mpi,                 ONLY: p_start, p_stop, p_init_communicators
  USE mo_machine,             ONLY: machine_setup
  USE mo_time_control,        ONLY: l_trigfiles, dt_start, dt_stop, delta_time, init_manager, init_times, ec_manager_init, &
                                    time_reset, lstart, lresume
  USE mo_jsbach_interface,    ONLY: get_dates
  USE mo_jsbach_lctlib,       ONLY: init_lctlib, lctlib_type
  USE mo_jsbach_grid,         ONLY: grid_type, domain_type, init_grid, init_domain, nidx, kstart, kend
  USE mo_jsbach_version,      ONLY : jsbach_init_version, jsbach_label_run
  USE mo_cbal_cpools,         ONLY: update_Cpools, N_process
  USE mo_kind,                ONLY: dp 
  USE mo_cbalone_io,          ONLY: getDailyData, writeSingleTimeStep
  USE mo_cbal_bethy,          ONLY: cbalance_type, nbalance_type,    &
                                    read_nitrogen_deposition_fluxes, &
                                    read_nitrogen_fertilizer_fluxes  !dk: Ndep originally other function in cbalone_memory
  USE mo_land_surface,        ONLY: land_surface_type
  USE mo_soil,                ONLY: soil_param_type
  USE mo_cbalone_memory,      ONLY: cbal_offline_type, Nbal_offline_type, &
                                    initSurface, initSoil, initializePools, initializeYassoPools, initCbalance, &
                                    initializeNPools, initNbalance,  &
                                    dynveg_clim_type
  USE mo_cbal_landcover_change,    ONLY: landuse_transitions_type, &
                                         read_landcover_fractions, read_landuse_transitions, read_harvest, &
                                         do_landcover_change, do_landuse_transitions, &
                                         landcover_change, landuse_transitions
  USE mo_cbalone_landcover_change, ONLY: initLandCoverChange
  USE mo_cbalone_dynveg,           ONLY: init_offline_dynveg, update_offline_climbuf
  USE mo_cbal_parameters,     ONLY: frac_wood_2_atmos, frac_green_2_atmos, frac_reserve_2_atmos, frac_mobile_2_atmos, &
                                    cn_green, cn_woods, cn_litter_green, cn_litter_wood, cn_slow
  USE mo_data_scaling,        ONLY: data_scaling_type, initDataScaling, scaleDailyData
  USE mo_cbal_parameters,     ONLY: config_cbal_parameters
  USE mo_util_string,         ONLY: tolower
  USE mo_netcdf,              ONLY: nf_max_name
  USE mo_jsbach_constants,    ONLY: molarMassCO2_kg
  USE mo_jsbach,              ONLY: nml_unit, new_day, new_month, new_year, debug, debug_Cconservation, missing_value
  USE mo_dynveg,              ONLY: dynveg, dynveg_options, config_dynveg, update_dynveg
  USE mo_climbuf,             ONLY: climbuf_type
  USE mo_namelist,            ONLY: open_nml, position_nml, POSITIONED
  USE mo_disturbance,         ONLY: disturbance_type, disturbance, config_disturbance, update_disturbance
  USE mo_disturbance_thonicke,ONLY: fire_TH_diag
  USE mo_exception,           ONLY: finish
  USE mo_input,               ONLY: InputUpdate, InputClose
  USE mo_input_calendar,      ONLY: INPUT_CALENDAR_JULIAN, INPUT_CALENDAR_GREGORIAN
  USE mo_input_interface,     ONLY: InputInitialize, INPUT_MODEL_CBALANCE
#ifdef NAG
  USE f90_unix_io,            ONLY: flush
#endif

  implicit none

  !! === PARAMETERS ====================================================================================
  
  ! --- Standardnames for input files ----

  character(nf_max_name),parameter :: iniFileVegetation =  "jsbach.nc" !! File with information on grid and initial landcover
  character(nf_max_name),parameter :: lctlib_file = "lctlib.def" !! File name of the Land cover library
  character(len=*),parameter :: outStemName =        "Cbalone"   !! Output file names are constructed as:
                                                                 !!     outStemName.xxxxxx.yyyymm.nc
                                                                 !!           where xxxxxx is the repetion number

  ! --- Elements of filename for driver data (NPP,LAI,alpha,temperature)
  ! Note: It is assumed that the filenames have the form  experiment_yyyymm.specificName.nc
  !       The infornmation on experiment and specific name are taken from run.def

  character(len=128) :: experiment
  character(len=128) :: specificName

  ! --- start and end of integration (read in from namelist)

  integer   :: climate_yearstart   !! first year from external (climate) data
  integer   :: climate_yearend     !! last year  from external (climate) data
  integer   :: run_year_first !! First absolute year of run 
  integer   :: run_year_last  !! Last absolute year of run

  character(len=5) :: out_interval !! output interval

  character(len=1024) :: driver_data_path !! path where the driver data (LAI,NPP,alpha,temperature) are found

  ! --- Control of initialization of carbon pools (options read in from run.def)

  logical             :: read_cpools         !! true: Cpools are initialised from 'cpool_file_name' 
  logical             :: read_npools         !! true: Npools are initialised from 'npool_file_name' 
  logical             :: read_ndepo          !! true: Nitrogen deposition fields are initialised from 'ndepo_file_name' 
  logical             :: read_nfert

  character(len=1024) :: cpool_file_name     !! file name, including path, for initial data of the carbon pools
  character(len=1024) :: npool_file_name     !! file name, including path, for initial data of the nitrogen pools 
  character(len=1024) :: ndepo_file_name     !! file name, including path, for nitrogen deposition fields
  character(len=1024) :: nfert_file_name
 

  ! --- Submodel controls

  logical             :: use_dynveg           !! true if dynamic vegetation module is used
  logical             :: use_disturbance      !! true if disturbances should be calculated independent of dynveg
  logical             :: with_nitrogen        !! true if besides carbon als nitrogen is cycled  
  logical             :: with_yasso           !! true if yasso model is used ! DSG 2.6.2013 
  logical             :: use_external_landcover_maps  !! true to simulate landuse change by reading external maps
  logical             :: with_landuse_transitions     !! true to calculate cmip5 landuse transitions
  character(len=128)  :: input_verbose

  ! --- structures -----------

  TYPE(lctlib_type)             :: lctlib
  TYPE(grid_type)               :: grid
  TYPE(domain_type)             :: domain
  TYPE(cbal_offline_type)       :: cbalance_diag
  TYPE(cbalance_type)           :: cbalance
  TYPE(nbalance_type)           :: nbalance
  TYPE(Nbal_offline_type)       :: nbalance_offline
  TYPE(soil_param_type)         :: soil_param
  TYPE(land_surface_type)       :: surface
  TYPE(data_scaling_type)       :: data_scaling
  TYPE(dynveg_clim_type)        :: dynveg_clim
  TYPE(climbuf_type)            :: climbuf

  ! --- allocatables

  real(dp),allocatable,dimension(:,:)  ::  specificLeafArea_C
  real(dp),allocatable,dimension(:,:)  ::  reserveC2leafC
  real(dp),allocatable,dimension(:,:)  ::  Max_LAI !! maximum value of leaf area index
  real(dp),allocatable,dimension(:,:)  ::  veg_fract_correction
  real(dp),allocatable,dimension(:,:)  ::  frac_npp_2_woodPool
  real(dp),allocatable,dimension(:,:)  ::  frac_npp_2_reservePool
  real(dp),allocatable,dimension(:,:)  ::  frac_npp_2_exudates
  real(dp),allocatable,dimension(:,:)  ::  frac_green_2_herbivory
  real(dp),allocatable,dimension(:,:)  ::  tau_Cpool_litter_green
  real(dp),allocatable,dimension(:,:)  ::  tau_Cpool_litter_wood
  real(dp),allocatable,dimension(:,:)  ::  tau_Cpool_woods
  real(dp),allocatable,dimension(:,:)  ::  LAI_shed_constant
  real(dp),allocatable,dimension(:,:)  ::  frac_C_litter_green2atmos
  real(dp),allocatable,dimension(:,:)  ::  Max_C_content_woods
  real(dp),allocatable,dimension(:,:)  ::  coverFractBox
  real(dp),allocatable,dimension(:,:)  ::  coverFract_interpolated  !! The landcover fractions interpolated throughout the year
  real(dp),allocatable,dimension(:)    ::  carbon_2_atmos, carbon_2_litterGreenPools
  real(dp),allocatable,dimension(:)    ::  CO2_emission_landcover_change, CO2_emission_harvest, CO2_flux_dynveg
  real(dp),allocatable,dimension(:)    ::  carbon_2_litterWoodPool_ag, carbon_2_litterWoodPool_bg
  real(dp),allocatable,dimension(:,:,:) ::  LeafLit_coef    !! Yasso parameters from lctlib
  real(dp),allocatable,dimension(:,:,:) ::  WoodLit_coef    !!         ''
  real(dp),allocatable,dimension(:,:)   ::  WoodLitterSize  !!         ''
  real(dp),allocatable,dimension(:)    ::  nitrogen_2_atmos, nitrogen_2_litterGreenPools
  real(dp),allocatable,dimension(:)    ::  nitrogen_2_litterWoodPool_ag,nitrogen_2_litterWoodPool_bg,nitrogen_2_SMINNpool
  real(dp),allocatable,dimension(:,:)  ::  litter_flux   !! dummy for subroutine call; not used in cbalone
  LOGICAL, allocatable                 ::  glacier_1d(:) !! auxiliary array (workaround for memory leak on aix) 

  ! --- other locals

  integer             ::  itile
  integer             ::  day,nday
  integer             ::  climate_year,month,number_of_days_in_year
  character(len=128)  ::  outFile
  character(len=128)  ::  inFile
  integer             ::  run_year !! counts years including repeats
  integer             ::  outputInterval
  character(len=1024) ::  CO2_file,climatology_diff_file,ndeposition_diff_file
  integer             ::  ref_year_past,ref_year_recent !! Only used when CO2-scaling is switched on
  LOGICAL             ::  input_scaling, input_scaling_ndep
  CHARACTER(len=16)   ::  lcc_forcing_type
  INTEGER             ::  lcc_scheme
  REAL(dp)            ::  step
  INTEGER             ::  read_status, f_unit   ! error status reading the namelist
  INTEGER             ::  IO_timestep, istep
  INTEGER             ::  apu


  ! --- constants ----

  integer,parameter   :: sec_per_day         = 86400
  integer,parameter   :: INTERVAL_DAY=1,INTERVAL_MONTH=2,INTERVAL_YEAR=3

  INCLUDE 'cbalone_ctl.inc'

  !! === GO =====================================================================

  CALL jsbach_init_version

  !! --- initialization of MPI
  call p_start
  CALL p_init_communicators(1, 1, 0) !! nproca=1, nprocb=1, nprocio=0
  CALL machine_setup

  !! --- read namelist cbalone_ctl
  driver_data_path = ''
  experiment = ''
  specificName = ''
  climate_yearstart = -9999
  climate_yearend = -9999
  run_year_first = -9999
  run_year_last = -9999
  out_interval = 'YEAR'
  input_scaling = .FALSE.
  input_scaling_ndep = .FALSE.
  co2_file = 'co2.nc'
  climatology_diff_file = 'climdiff.nc'
  ndeposition_diff_file = 'ndepdiff.nc'
  ref_year_past = -9999
  ref_year_recent = -9999
  use_dynveg = .FALSE.
  use_disturbance = .FALSE.
  with_nitrogen = .FALSE.
  with_yasso = .FALSE.
  lcc_forcing_type = 'none'
  lcc_scheme = 1
  use_external_landcover_maps = .false.
  with_landuse_transitions = .false.
  debug = .FALSE.
  debug_Cconservation = .FALSE.
  use_external_landcover_maps = .false.
  with_landuse_transitions   = .false.     !dk initialisation was missing
  !!
  read_cpools=.FALSE.      !! false: Carbon pools are initialized with zero
                           !! true: Carbon pools initialized from an external file
  cpool_file_name='Cpools.nc'
  !!
  read_npools=.FALSE.      !! false: Nitrogen pools are initialized with default values
                           !! true: Nitrogen pools initialized from an external file
  npool_file_name='Npools.nc'
  !!
  read_ndepo =.FALSE.      !! false: Default initialization is used for the Nitrogen deposition field
                           !! true: Nitrogen deposition field is initialized from external file
  ndepo_file_name ='Ndepo.nc'
  !!
  read_nfert =.FALSE.      !! false: Default initialization is used for the Nitrogen deposition field
                           !! true: Nitrogen deposition field is initialized from external file
  nfert_file_name ='Nfert.nc'
  input_verbose = ''
  !!
  !!
  nml_unit =  open_nml ('namelist.cbalone')
  f_unit = position_nml ('CBALONE_CTL', nml_unit, status=read_status)
  SELECT CASE (read_status)
  CASE (POSITIONED)
     READ (f_unit, cbalone_ctl)
  END SELECT

  IF (TRIM(driver_data_path) == '') &
       STOP 'No entry for DRIVER_DATA_PATH found in namelist.'
  IF (TRIM(experiment) == '') &
       STOP 'No entry for keyword EXPERIMENT found in namelist.'
  IF (TRIM(specificName) == '') &
       STOP 'No entry for keyword SPECIFICNAME found in namelist.'
  IF (climate_yearstart > climate_yearend) &
       STOP 'Namelist paramter climate_yearstart > climate_yearend.'
  IF (run_year_first > run_year_last) &
       STOP 'Namelist paramter run_year_first > run_year_last.'

  out_interval=toLower(out_interval)
  select case(TRIM(out_interval))
     case('day')
        outputInterval=INTERVAL_DAY
     case('month')
        outputInterval=INTERVAL_MONTH
     case('year')
        outputInterval=INTERVAL_YEAR
     case default
        stop 'No valid value for keyword OUT_INTERVAL found in namelist.'
  end select

  write(*,*) " === Cbalone info =============="
  write(*,*) " Climate data: First year: ",climate_yearstart
  write(*,*) " Climate data: Last  year: ",climate_yearend
  write(*,*) " First year of run : ",run_year_first
  write(*,*) " Last  year of run : ",run_year_last
  write(*,*) " Output interval   : ",trim(out_interval)
  write(*,*) " ==============================="

  IF (read_cpools) THEN
     WRITE(*,*) "Init file for carbon pools: ",TRIM(cpool_file_name)
  END IF

  IF (with_yasso) THEN
     WRITE(*,*) "The Yasso soil C model is active."
  END IF

  IF (with_nitrogen) THEN
     IF (with_yasso) THEN
        STOP "Sorry -- current implimentation does not allow yasso soil model together with nitrogen cycling."
     END IF
     IF (lcc_scheme /= 1) CALL finish('cbalone','Currently nitrogen cycle cannot be used togeter with anthropogenic pools')
     WRITE(*,*) "Nitrogen cycling is active."
     IF (read_npools) THEN
        WRITE(*,*) "Init file for nitrogen pools: ",TRIM(npool_file_name)
     END IF
     IF (read_ndepo) THEN
        WRITE(*,*) "Input file for nitrogen deposition fields: ",TRIM(ndepo_file_name)
     END IF
     IF (read_nfert) THEN
        WRITE(*,*) "Input file for nitrogen fertilizer fields: ",TRIM(nfert_file_name)
     END IF
  END IF

  IF (use_dynveg .AND. .NOT. use_disturbance) THEN
     WRITE(*,*) 'Dynamic vegetation does not make sense without disturbances.'
     WRITE(*,*) '---------'
     WRITE(*,*) 'WARNING: Disturbance calculations are switched on'
     WRITE(*,*) '---------'
     use_disturbance = .TRUE.
  ENDIF
  IF (use_dynveg)      WRITE(*,*) "Dynamic vegetation module is active."
  IF (use_disturbance) WRITE(*,*) "Disturbance module is active."

  IF(use_dynveg .and. with_nitrogen) then
     WRITE(*,*) "This run is with dynamic vegetation and nitrogen."
  END IF

  IF (TRIM(lcc_forcing_type) == 'maps') THEN
     use_external_landcover_maps = .TRUE.
     WRITE(*,*) "This is a run with anthropogenic landcover change."
     IF (use_dynveg) THEN
        STOP "Sorry -- current implimentation does not allow dynamical vegetation together with landcover change."
     END IF
     IF (with_yasso) THEN
        STOP "Sorry -- current implimentation does not allow yasso soil model together with landcover change."
     END IF
  ELSE IF (TRIM(lcc_forcing_type) == 'transitions') THEN
     with_landuse_transitions = .TRUE.
     WRITE(*,*) "This is a run with anthropogenic landuse transitions"
  END IF
 
  IF (lcc_scheme==1)   WRITE(*,*) "This is a run with lcc carbon allocation following JSBACH standard"
  IF (lcc_scheme==2)   WRITE(*,*) &
       "This is a run with lcc carbon allocation following the Grand Slam Protocol (Houghton et al., 1983)"

  data_scaling%do_input_scaling = input_scaling
  if (data_scaling%do_input_scaling) then
     write(*,*) "This is a run with climate input scaling"
     write(*,*) "CO2 development taken from file: ",trim(CO2_file)
     write(*,*) "Climatological differences of driver data are taken from file: ",TRIM(climatology_diff_file)
     write(*,*) "For input data scaling reference years are ",ref_year_past," and ",ref_year_recent," are used."
  else
     write(*,*) "No input scaling performed in this run"
  end if

  IF (with_nitrogen) THEN
     IF (data_scaling%do_input_scaling) THEN
        write(*,*) "This is a run with nitrogen input scaling"
        write(*,*) "CO2 development taken from file: ",trim(CO2_file)
        write(*,*) "Ndeposition differences of driver data are taken from file: ",TRIM(ndeposition_diff_file)
        write(*,*) "For input data scaling reference years are ",ref_year_past," and ",ref_year_recent," are used."
     ELSE
        write(*,*) "No ndeposition input scaling performed in this run"
     END IF
  END IF

  !! --- initialize the echam time manager
  dt_start = (/run_year_first,1,1,0,0,0/)
  dt_stop  = (/run_year_last,12,31,0,0,0/)
  delta_time = 86400                                ! time step in seconds (one day)
  lresume=.FALSE.                                   ! initialized run
  CALL get_dates(IO_timestep, istep)                ! get timestep
  CALL ec_manager_init(IO_timestep, istep)          ! time manager initialization
  CALL init_manager
  CALL init_times                                   ! initialize all dates
  new_day = .TRUE.                                  ! the internal timestep is a day
  l_trigfiles = .FALSE.

  !! --- initialize input module
  climate_year = climate_yearend-climate_yearstart+1
  day = INPUT_CALENDAR_GREGORIAN
  IF (climate_year < 400 .AND. mod(climate_year,4)==0) day = INPUT_CALENDAR_JULIAN
  CALL InputInitialize(INPUT_MODEL_CBALANCE,Debug=input_verbose,Calendar=day)
  
  !! --- initialization of lctlib
  CALL init_lctlib(lctlib_file, lctlib)
  WRITE(*,*) "LCT-lib initialized"

  !! --- initialization of the grid
  CALL init_grid(grid, .TRUE., lresume, iniFileVegetation) !! standalone=.TRUE.
  WRITE(*,*) "Grid initialized."

  CALL jsbalone_init_decomp(grid%nlon, grid%nlat, grid%mask, 1, 1, -1) !! nproca=1, nprocb=1, nproma=-1

  !! --- initialization of the domain
  CALL init_domain(grid, domain, cbalone_flag=.TRUE.)
  write(*,*) "Domain initialized."

  ! --- initialization of surface specific values
  call initSurface(iniFileVegetation, grid, lctlib, surface)
  write(*,*) "Land surface initialized."

  ! --- initialization of soil specific values
  call initSoil(iniFileVegetation, grid, soil_param)
  write(*,*) "Soil initialized."

  ! --- init memory for cbalance structure

  CALL initCbalance(grid, surface%ntiles, cbalance, cbalance_diag,with_landuse_transitions,lcc_scheme,with_yasso)
  write(*,*) "Cbalance initialized."

  ! --- init memory for nbalance structure

  IF (with_nitrogen) THEN
     call initNbalance(grid, surface%ntiles, nbalance, nbalance_offline)
     write(*,*) "Nbalance initialized."
  END IF

  ! --- initialization of cbal_parameters
  CALL config_cbal_parameters(nml_unit)

  ! --- initialization of memory for lancover change emissions
  IF (use_external_landcover_maps .OR. with_landuse_transitions) THEN

     CALL initLandCoverChange(grid, surface%ntiles, with_landuse_transitions, with_nitrogen, landcover_change, &
          landuse_transitions, lcc_scheme)

     !! --- printout parameters 
     
     WRITE(*,*) "=== Parameters Landcover Change ============"
     WRITE(*,*) "   frac_wood_2_atmos=", frac_wood_2_atmos
     WRITE(*,*) "  frac_green_2_atmos=", frac_green_2_atmos
     WRITE(*,*) "frac_reserve_2_atmos=",frac_reserve_2_atmos
     IF(with_nitrogen) then
        WRITE(*,*) " frac_mobile_2_atmos=", frac_mobile_2_atmos
     END IF
     WRITE(*,*) "============================================"
  END IF

  ! --- initialization of Cpools
  IF (read_cpools) THEN
     CALL initializePools(grid, cbalance, cbalance_diag, surface, TRIM(cpool_file_name))
  ELSE
     CALL initializePools(grid, cbalance, cbalance_diag, surface)
  END IF

  ! --- initialization of Yasso Cpools 
  IF (with_yasso) THEN
    IF (read_cpools) THEN
       CALL initializeYassoPools(grid, cbalance, surface, TRIM(cpool_file_name))
    ELSE
       CALL initializeYassoPools(grid, cbalance, surface)
    END IF
  END IF

  ! --- initialization of Npools
  IF (with_nitrogen) THEN
     IF (read_Npools) THEN
        CALL initializeNPools(grid, surface%ntiles, nbalance, TRIM(npool_file_name))
     ELSE
        CALL initializeNPools(grid, surface%ntiles, nbalance, cbalance=cbalance)
     END IF
  END IF 

  ! --- initialise memory for dynamic vegetation
  IF (use_disturbance) THEN
    CALL config_disturbance(surface, grid, domain, lctlib, nml_unit, with_yasso, with_nitrogen, cbalone=.true.)
    IF (.NOT. read_cpools) disturbance%fuel(:) = missing_value ! Needed to avoid problems on coldstart (non-initialized Cpools)
  ENDIF
  IF (use_dynveg) THEN
     CALL config_dynveg()
     IF (use_dynveg .and. read_cpools .AND..NOT. dynveg_options%read_fpc) THEN
        STOP "If Cpools are read, dynveg needs matching cover fractions! set read_fpc=.true." 
     END IF
  END IF
  IF (use_dynveg .OR. use_disturbance) &
     CALL init_offline_dynveg(grid, lctlib, surface, dynveg_clim, climbuf, use_dynveg, use_disturbance)

  ! --- initialization of memory for input data scaling
  IF(data_scaling%do_input_scaling) THEN
     CALL initDataScaling(CO2_file, climatology_diff_file, grid,surface%ntiles, ref_year_past, ref_year_recent, &
          run_year_first,run_year_last,data_scaling)
  END IF

  !! Allocate local memory

  allocate(frac_npp_2_woodPool(1:grid%nland,surface%ntiles))
  allocate(frac_npp_2_reservePool(1:grid%nland,surface%ntiles))
  allocate(frac_npp_2_exudates(1:grid%nland,surface%ntiles))
  allocate(frac_green_2_herbivory(1:grid%nland,surface%ntiles))
  allocate(tau_Cpool_litter_green(1:grid%nland,surface%ntiles))
  allocate(tau_Cpool_litter_wood(1:grid%nland,surface%ntiles))
  allocate(tau_Cpool_woods      (1:grid%nland,surface%ntiles))
  allocate(LAI_shed_constant(1:grid%nland,surface%ntiles))
  allocate(frac_C_litter_green2atmos(1:grid%nland,surface%ntiles))
  allocate(Max_C_content_woods(1:grid%nland,surface%ntiles))
  allocate(specificLeafArea_C(1:grid%nland,surface%ntiles))
  allocate(reserveC2leafC(1:grid%nland,surface%ntiles))
  allocate(Max_LAI(1:grid%nland,surface%ntiles))
  allocate(veg_fract_correction(1:grid%nland,surface%ntiles))
  allocate(coverFractBox(1:grid%nland,surface%ntiles))
  allocate(carbon_2_atmos(1:grid%nland))
  allocate(CO2_emission_landcover_change(1:grid%nland))
  CO2_emission_landcover_change(:)=0.0_dp
  allocate(CO2_emission_harvest(1:grid%nland))
  CO2_emission_harvest(:)=0.0_dp
  allocate(CO2_flux_dynveg(1:grid%nland))
  CO2_flux_dynveg(:)=0.0_dp
  allocate(carbon_2_litterGreenPools(1:grid%nland))
  allocate(carbon_2_litterWoodPool_ag(1:grid%nland))
  allocate(carbon_2_litterWoodPool_bg(1:grid%nland))
  allocate(LeafLit_coef(1:grid%nland,surface%ntiles,5))
  allocate(WoodLit_coef(1:grid%nland,surface%ntiles,5))
  allocate(WoodLitterSize(1:grid%nland,surface%ntiles))

  IF (with_nitrogen) THEN
     allocate(nitrogen_2_atmos(1:grid%nland))
     nitrogen_2_atmos(:)=0.0_dp
     allocate(nitrogen_2_litterGreenPools(1:grid%nland))
     allocate(nitrogen_2_litterWoodPool_ag(1:grid%nland))
     allocate(nitrogen_2_litterWoodPool_bg(1:grid%nland))
     allocate(nitrogen_2_SMINNpool(1:grid%nland))
  END IF
  allocate(coverFract_interpolated(1:grid%nland,surface%ntiles))
  allocate(litter_flux(1:grid%nland,surface%ntiles))
  litter_flux(:,:)=0.0_dp 
  allocate(glacier_1d(1:grid%nland))

  DO itile=1,surface%ntiles
    frac_npp_2_woodPool(1:grid%nland,itile)       = lctLib%frac_npp_2_woodPool(surface%cover_type(1:grid%nland,itile))
    frac_npp_2_reservePool(1:grid%nland,itile)    = lctLib%frac_npp_2_reservePool(surface%cover_type(1:grid%nland,itile))
    frac_npp_2_exudates(1:grid%nland,itile)       = lctLib%frac_npp_2_exudates(surface%cover_type(1:grid%nland,itile))
    frac_green_2_herbivory(1:grid%nland,itile)    = lctLib%frac_green_2_herbivory(surface%cover_type(1:grid%nland,itile))
    tau_Cpool_litter_green(1:grid%nland,itile)    = lctLib%tau_Cpool_litter_leaf(surface%cover_type(1:grid%nland,itile))
    tau_Cpool_litter_wood(1:grid%nland,itile)     = lctLib%tau_Cpool_litter_wood(surface%cover_type(1:grid%nland,itile))
    tau_Cpool_woods      (1:grid%nland,itile)     = lctLib%tau_Cpool_woods      (surface%cover_type(1:grid%nland,itile))
    LAI_shed_constant(1:grid%nland,itile)         = lctLib%LAI_shed_constant(surface%cover_type(1:grid%nland,itile))
    frac_C_litter_green2atmos(1:grid%nland,itile) = lctLib%frac_C_litter_green2atmos(surface%cover_type(1:grid%nland,itile))
    Max_C_content_woods(1:grid%nland,itile)       = lctLib%Max_C_content_woods(surface%cover_type(1:grid%nland,itile))
    specificLeafArea_C(1:grid%nland,itile)        = lctLib%specificLeafArea_C(surface%cover_type(1:grid%nland,itile))
    reserveC2leafC(1:grid%nland,itile)            = lctLib%reserveC2leafC(surface%cover_type(1:grid%nland,itile))
    Max_LAI(1:grid%nland,itile)                   = lctLib%MaxLAI(surface%cover_type(1:grid%nland,itile))
    WoodLitterSize(1:grid%nland,itile)            = lctLib%WoodLitterSize(surface%cover_type(1:grid%nland,itile))
    DO apu=1,5
      LeafLit_coef(1:grid%nland,itile,apu)        = lctLib%LeafLit_coef(surface%cover_type(1:grid%nland,itile),apu)
      WoodLit_coef(1:grid%nland,itile,apu)        = lctLib%WoodLit_coef(surface%cover_type(1:grid%nland,itile),apu)
    END DO
    !DSG: check to see if the yasso parameters are correct
    !   the sum over the last dimension must be 1.0

    IF ( SUM(LeafLit_coef(1,itile,:)) /= 1.0 ) THEN
         WRITE(*,*)  "Leaf litter coefficient has invalid values; must be 1, but is:" 
         WRITE(*,*)  SUM(LeafLit_coef(1,itile,:)), "for tile=", itile 
    END IF        
    IF ( SUM(WoodLit_coef(1,itile,:)) /= 1.0 ) THEN
         WRITE(*,*)  "Wood litter coefficient has invalid values; must be 1, but is:" 
         WRITE(*,*)  SUM(WoodLit_coef(1,itile,:)), "for tile=", itile 
    END IF        
    
    veg_fract_correction(1:grid%nland,itile)      = 1.0_dp - exp(-lctLib%MaxLAI(surface%cover_type(1:grid%nland,itile)) &
                                                           / lctLib%ClumpinessFactor(surface%cover_type(1:grid%nland,itile)))
  END DO

  CALL nullifyFields(cbalance_diag, disturbance, fire_TH_diag, landcover_change,&
       landuse_transitions, lcc_scheme)

  ! --- read in land cover fractions in case of anthropogenic landcover change

  IF (use_external_landcover_maps) THEN !! Note: Here the landcover map of the year before the first year is needed!
     !! Overwrite cover fractions read in from jsbach init file
     CALL read_landcover_fractions(run_year_first -1, grid, domain, surface)
     surface%cover_fract = landcover_change%LCC_coverFract_target
  end if

  CALL jsbach_label_run

  ! --- start of repeat loop

  run_year=run_year_first - 1
  step = -1                     !! initialize timestep counter

  RepeatLoop: do

     ! --- loop over all years of the data

     do climate_year=climate_yearstart,climate_yearend

        run_year = run_year+1 !! update year counter
        new_year = .true.

        write(*,*) "run_year/climate_year=",run_year,"/",climate_year

        !! Determine name of output file

        if(outputInterval /= INTERVAL_DAY) then
           write (outfile,'(a,I6.6,a,I4.4,a)')trim(outStemName)//".", run_year,".", climate_year, ".nc"
        end if

        !! read in fertilizer for new year
        if (With_Nitrogen .AND. read_nfert) then
           call read_nitrogen_fertilizer_fluxes(run_year, Grid,Domain, Nbalance,     &
                                                surface%is_crop(:,1:surface%ntiles), &
                                                surface%ntiles)
        end if

        !! read in deposition for new year
        if (with_Nitrogen .and. read_Ndepo) then
           call read_nitrogen_deposition_fluxes(run_year, grid, domain, nbalance)
        end if

        !! Read in new landcover fractions for new year

        if (use_external_landcover_maps) then
           CALL read_landcover_fractions(run_year, grid, domain, surface)
        end if

        !! Read in new landuse transition matrix for new year

        if (with_landuse_transitions) then
           CALL read_landuse_transitions(run_year, grid, domain)
           CALL read_harvest(run_year, grid, domain)
        end if

        !! loop over all months

        number_of_days_in_year = 0
        do month=1,12

           new_month = .true.

           write(*,*) "run_year,month=",run_year,month

           !! Construct name of output file
           if(outputInterval == INTERVAL_DAY) then
              write(outFile,'(a,I6.6,I2.2,a,I4.4,a)') trim(outStemName)//".", run_year, month, ".", climate_year, ".nc"
           end if

           ! --- get data for a full month
           write(inFile,'(a,I4.4,I2.2,a,a)') trim(experiment), climate_year, month, trim(specificName)//".nc"
           call getDailyData(trim(driver_data_path)//trim(inFile),grid, use_dynveg,use_disturbance, with_nitrogen, &
                                  with_yasso, cbalance_diag, nbalance_offline, dynveg_clim, surface%ntiles, nday)
           ! --- scale the data according to CO2 concentration

           IF (data_scaling%do_input_scaling) CALL scaleDailyData(nday, data_scaling, run_year, month, cbalance_diag)

           !! --- start day loop

           do day=1,nday

              CALL InputUpdate()

              step = step + 1  !! update timestep counter

              !! determine fractions in box to compute total C and N for mass conservation check

              coverFractBox(:,:) = surface%cover_fract(:,:) * veg_fract_correction(:,:)      &
                                   * spread(surface%veg_ratio_max(:),DIM=2,NCOPIES=surface%ntiles)

              !! remember total carbon and nitrogen for testing C and N conservation

              IF (with_yasso) THEN
                 cbalance_diag%box_test_Ccons(:) = SUM(coverFractBox(:,:)                                                       &
                                                   * (   cbalance%Cpool_green(:,:)                                              &
                                                       + cbalance%Cpool_woods(:,:)                                              &
                                                       + cbalance%YCpool_acid_ag1(:,:)                                           &
                                                       + cbalance%YCpool_water_ag1(:,:)                                          &
                                                       + cbalance%YCpool_ethanol_ag1(:,:)                                        &
                                                       + cbalance%YCpool_nonsoluble_ag1(:,:)                                     &
                                                       + cbalance%YCpool_acid_bg1(:,:)                                           &
                                                       + cbalance%YCpool_water_bg1(:,:)                                          &
                                                       + cbalance%YCpool_ethanol_bg1(:,:)                                        &
                                                       + cbalance%YCpool_nonsoluble_bg1(:,:)                                     &
                                                       + cbalance%YCpool_humus_1(:,:)                                             &
                                                       + cbalance%YCpool_acid_ag2(:,:)                                           &
                                                       + cbalance%YCpool_water_ag2(:,:)                                          &
                                                       + cbalance%YCpool_ethanol_ag2(:,:)                                        &
                                                       + cbalance%YCpool_nonsoluble_ag2(:,:)                                     &
                                                       + cbalance%YCpool_acid_bg2(:,:)                                           &
                                                       + cbalance%YCpool_water_bg2(:,:)                                          &
                                                       + cbalance%YCpool_ethanol_bg2(:,:)                                        &
                                                       + cbalance%YCpool_nonsoluble_bg2(:,:)                                     &
                                                       + cbalance%YCpool_humus_2(:,:)                                             &
                                                       + cbalance%Cpool_reserve(:,:) )                                          &
                                                      ,DIM=2)

              ELSE

                 cbalance_diag%box_test_Ccons(:) = sum(coverFractBox(:,:)                                                       &
                                                   * (  cbalance%Cpool_green(:,:)        + cbalance%Cpool_woods(:,:)            &
                                                      + cbalance%Cpool_litter_green_ag(:,:)                                     &
                                                      + cbalance%Cpool_litter_wood_ag(:,:) + cbalance%Cpool_litter_wood_bg(:,:) &
                                                      + cbalance%Cpool_litter_green_bg(:,:) + cbalance%Cpool_slow(:,:)          &
                                                      + cbalance%Cpool_reserve(:,:) )                                           &
                                                     ,DIM=2)
              END IF

              IF (use_external_landcover_maps .OR. with_landuse_transitions) THEN
                 IF (lcc_scheme == 2) THEN
                    cbalance_diag%box_test_Ccons(:) = cbalance_diag%box_test_Ccons(:) &
                         + (  cbalance%Cpool_onSite(:)       &
                            + cbalance%Cpool_paper(:)        &
                            + cbalance%Cpool_construction(:) &
                           ) * surface%veg_ratio_max(:)
                 ENDIF
                 IF (with_landuse_transitions .AND. (lcc_scheme == 2)) THEN
                    cbalance_diag%box_test_Ccons(:) = cbalance_diag%box_test_Ccons(:) &
                         + (  cbalance%Cpool_paper_harvest(:)        &
                            + cbalance%Cpool_construction_harvest(:) &
                           ) * surface%veg_ratio_max(:)
                 ENDIF
              ENDIF
                        
              IF (with_nitrogen) THEN
                 nbalance_offline%box_test_Ncons(:) = sum(coverFractBox(:,:)                    &
                      * (  nbalance%Npool_green(:,:)        + nbalance%Npool_woods(:,:)         &
                      + nbalance%Npool_litter_green_ag(:,:)                                     &
                      + nbalance%Npool_litter_wood_ag(:,:) + nbalance%Npool_litter_wood_bg(:,:) &
                      + nbalance%Npool_litter_green_bg(:,:) + nbalance%Npool_slow(:,:)          &
                      + nbalance%Npool_mobile(:,:)       + nbalance%SMINN_pool(:,:)  )          &
                      ,DIM=2)
              END IF

              !
              ! --- update dynamic vegetation
              !
              IF (use_dynveg) THEN
                 CALL update_offline_climbuf(day, dynveg_clim, climbuf, dynveg)
                 CALL update_dynveg(lstart, lresume, nidx, kstart, kend, surface%ntiles,                       &
                      use_disturbance, .FALSE., with_yasso, with_nitrogen, lctlib, surface, cbalance, climbuf, &
                      specificLeafArea_C(:,:), cbalance_diag%LAI_yDayMean(:,:,day), veg_fract_correction(:,:), &
                      CO2_flux_dynveg, LeafLit_coef(:,:,:),WoodLit_coef(:,:,:),                                &
                      nbalance=nbalance,cbalone_flag=.TRUE.)
              ELSE
                IF (use_disturbance .and. new_day) then
                  CALL update_offline_climbuf(day, dynveg_clim, climbuf, dynveg)
                  CALL update_disturbance(grid%nland,1,grid%nland,climbuf,Surface,Cbalance,with_yasso,with_nitrogen,Nbalance, &
                                          LctLib, veg_fract_correction, CO2_flux_dynveg, LeafLit_coef,WoodLit_Coef)
                ENDIF
              END IF

              ! --- do landcover change
             
              IF (use_external_landcover_maps) THEN
                 CALL do_landcover_change(surface%ntiles, lctlib, surface, cbalance, nbalance, with_nitrogen, &
                      with_yasso, LeafLit_coef, WoodLit_coef, &
                      surface%veg_ratio_max, veg_fract_correction, surface%cover_fract, &
                      CO2_emission_landcover_change, lcc=landcover_change, &
                      lcc_scheme=lcc_scheme)
              ELSE IF (with_landuse_transitions) THEN
                 CALL do_landuse_transitions(lctlib,surface,surface%ntiles, surface%is_naturalVeg, &
                      surface%is_C4vegetation, surface%is_grass, &
                      surface%is_crop, surface%is_pasture, surface%is_vegetation, surface%is_glacier, surface%veg_ratio_max,   &
                      veg_fract_correction, surface%cover_fract_pot, surface%cover_fract, cbalance, &
                      CO2_emission_landcover_change, CO2_emission_harvest, &
                      landcover_change, with_yasso, &
                      LeafLit_coef, WoodLit_coef, &
                      nbalance, with_nitrogen, lcc_scheme=lcc_scheme)
              END IF

              !! vegetation cover fractions may have changed, so redetermine fractions in box for later use here

              coverFractBox(:,:) = surface%cover_fract(:,:) * veg_fract_correction(:,:)      &
                                   * spread(surface%veg_ratio_max(:),DIM=2,NCOPIES=surface%ntiles)

              !! === update C and N pools =================

              IF (day > 1) THEN
                cbalance_diag%LAI_previousDayMean(1:grid%nland,1:surface%ntiles) = &
                  cbalance_diag%LAI_yDaymean(1:grid%nland,1:surface%ntiles,day-1)
              END IF

              cbalance%NPP_flux_correction(:,:) = 0._dp     
              cbalance%soil_respiration(:,:) = 0._dp
              cbalance%soil_respiration_pot(:,:) = 0._dp
              cbalance%excess_NPP(:,:)    = 0._dp        
              cbalance%root_exudates(:,:) = 0._dp          
              cbalance%Cflux_herbivory(:,:) = 0._dp   
              cbalance%Cflux_herbivory_LG(:,:)    = 0._dp         
              cbalance%Cflux_herbivory_2_atm(:,:) = 0._dp     

              !! Prepare for test of Carbon conservation across update Cpools by remembering the sum of the Carbon pools at the
              !! beginning of the day
              IF (.NOT.with_yasso) THEN
                   cbalance_diag%testCconserv(:,:) = cbalance%Cpool_green(:,:)  + cbalance%Cpool_woods(:,:)                   &
                                                   + cbalance%Cpool_litter_green_ag(:,:)                                      &
                                                   + cbalance%Cpool_litter_wood_ag(:,:) + cbalance%Cpool_litter_wood_bg(:,:)  &
                                                   + cbalance%Cpool_litter_green_bg(:,:) + cbalance%Cpool_slow(:,:)           &
                                                   + cbalance%Cpool_reserve(:,:)
              ELSE
                   cbalance_diag%testCconserv(:,:) = cbalance%Cpool_green          (:,:)         &
                                                   + cbalance%Cpool_woods          (:,:)         &
                                                   + cbalance%Cpool_reserve        (:,:)         &
                                                       + cbalance%YCpool_acid_ag1(:,:)                                           &
                                                       + cbalance%YCpool_water_ag1(:,:)                                          &
                                                       + cbalance%YCpool_ethanol_ag1(:,:)                                        &
                                                       + cbalance%YCpool_nonsoluble_ag1(:,:)                                     &
                                                       + cbalance%YCpool_acid_bg1(:,:)                                           &
                                                       + cbalance%YCpool_water_bg1(:,:)                                          &
                                                       + cbalance%YCpool_ethanol_bg1(:,:)                                        &
                                                       + cbalance%YCpool_nonsoluble_bg1(:,:)                                     &
                                                       + cbalance%YCpool_humus_1(:,:)                                             &
                                                       + cbalance%YCpool_acid_ag2(:,:)                                           &
                                                       + cbalance%YCpool_water_ag2(:,:)                                          &
                                                       + cbalance%YCpool_ethanol_ag2(:,:)                                        &
                                                       + cbalance%YCpool_nonsoluble_ag2(:,:)                                     &
                                                       + cbalance%YCpool_acid_bg2(:,:)                                           &
                                                       + cbalance%YCpool_water_bg2(:,:)                                          &
                                                       + cbalance%YCpool_ethanol_bg2(:,:)                                        &
                                                       + cbalance%YCpool_nonsoluble_bg2(:,:)                                     &
                                                       + cbalance%YCpool_humus_2(:,:)                                           
              END IF

              IF (with_nitrogen) THEN                                                   ! Note: reserve pool is free of nitrogen
                 nbalance_offline%testNconserv(:,:) = nbalance%Npool_green(:,:) + nbalance%Npool_woods(:,:)         &
                                         + nbalance%Npool_mobile(:,:)                                               &
                                         + nbalance%Npool_litter_green_ag(:,:)                                      &
                                         + nbalance%Npool_litter_wood_ag(:,:) + nbalance%Npool_litter_wood_bg(:,:)  &
                                         + nbalance%Npool_litter_green_bg(:,:) + nbalance%Npool_slow(:,:)           &
                                         + nbalance%SMINN_pool(:,:)     
              END IF

              IF (with_nitrogen) THEN
                 CALL N_process(cbalance%NPP_act_yDayMean(1:grid%nland,1:surface%ntiles),  &
                      nbalance%nfix_to_sminn(1:grid%nland,1:surface%ntiles),               &
                      nbalance%ndep_to_sminn(1:grid%nland,1:surface%ntiles),               &
                      SPREAD(nbalance%ndep_forc(1:grid%nland),DIM=2,NCOPIES=surface%ntiles), &
                      nbalance%nfert_forc_2d(1:grid%nland,1:surface%ntiles),               &
                      nbalance%nfert_to_sminn(1:grid%nland,1:surface%ntiles),              &
                      nbalance%sminn_to_denit(1:grid%nland,1:surface%ntiles),              &
                      nbalance%sminn_leach(1:grid%nland,1:surface%ntiles),                 &
                      nbalance%SMINN_pool(1:grid%nland,1:surface%ntiles),                  &
                      nbalance%NetEcosyst_N_flux(1:grid%nland,1:surface%ntiles),           &
                      nbalance_offline%Runoff_yDayMean(1:grid%nland,1:surface%ntiles,day), &  
                      cbalance_diag%alpha_yDayMean(1:grid%nland,1:surface%ntiles,day),     &
                      SPREAD(soil_param%maxMoisture(1:grid%nland),DIM=2,NCOPIES=surface%ntiles), &
                      nbalance%NPP_run_mean(1:grid%nland,1:surface%ntiles),                &
                      nbalance%N2O_emissions_depfix(1:grid%nland,1:surface%ntiles),        &
                      surface%is_vegetation(1:grid%nland,1:surface%ntiles),                &
                      nbalance%N2O_emissions_nfert(1:grid%nland,1:surface%ntiles)          &
                      )  

                 CALL update_Cpools(cbalance_diag%LAI_yDaymean(1:grid%nland,1:surface%ntiles,day),         &
                                    cbalance_diag%LAI_previousDayMean(1:grid%nland,1:surface%ntiles),      &
                                    cbalance_diag%NPP_yDayMean(1:grid%nland,1:surface%ntiles,day),         &
                                    cbalance_diag%topSoilTemp_yDayMean(1:grid%nland,1:surface%ntiles,day), &
                                    cbalance_diag%alpha_yDayMean(1:grid%nland,1:surface%ntiles,day),       &
                                    frac_npp_2_woodPool(1:grid%nland,1:surface%ntiles),               &
                                    frac_npp_2_reservePool(1:grid%nland,1:surface%ntiles),            &
                                    frac_npp_2_exudates(1:grid%nland,1:surface%ntiles),               &
                                    frac_green_2_herbivory(1:grid%nland,1:surface%ntiles),            &
                                    tau_Cpool_litter_green(1:grid%nland,1:surface%ntiles),            &
                                    tau_Cpool_litter_wood(1:grid%nland,1:surface%ntiles),             &
                                    tau_Cpool_woods(1:grid%nland,1:surface%ntiles),                   &
                                    LAI_shed_constant(1:grid%nland,1:surface%ntiles),                 &
                                    frac_C_litter_green2atmos(1:grid%nland,1:surface%ntiles),         &
                                    Max_C_content_woods(1:grid%nland,1:surface%ntiles),               &
                                    specificLeafArea_C(1:grid%nland,1:surface%ntiles),                &
                                    reserveC2leafC(1:grid%nland,1:surface%ntiles),                    &
                                    Max_LAI(1:grid%nland,1:surface%ntiles),                           &
                                    surface%is_vegetation(1:grid%nland,1:surface%ntiles),             &
                                    cbalance%Cpool_green(1:grid%nland,1:surface%ntiles),              &
                                    cbalance%Cpool_woods(1:grid%nland,1:surface%ntiles),              &
                                    cbalance%Cpool_reserve(1:grid%nland,1:surface%ntiles),            &
                                    cbalance%Cpool_litter_green_ag(1:grid%nland,1:surface%ntiles),    &
                                    cbalance%Cpool_litter_green_bg(1:grid%nland,1:surface%ntiles),    &
                                    cbalance%Cpool_litter_wood_ag(1:grid%nland,1:surface%ntiles),     &
                                    cbalance%Cpool_litter_wood_bg(1:grid%nland,1:surface%ntiles),     &
                                    cbalance%Cpool_slow(1:grid%nland,1:surface%ntiles),               &
                                    cbalance%soil_respiration(1:grid%nland,1:surface%ntiles),         &
                                    cbalance%soil_respiration_pot(1:grid%nland,1:surface%ntiles),     &
                                    cbalance%NPP_flux_correction(1:grid%nland,1:surface%ntiles),      &
                                    cbalance%excess_NPP(1:grid%nland,1:surface%ntiles),               &
                                    cbalance%root_exudates(1:grid%nland,1:surface%ntiles),            &
                                    litter_flux(1:grid%nland,1:surface%ntiles),                       &
                                    cbalance%Cflux_herbivory(1:grid%nland,1:surface%ntiles),          &
                                    cbalance%Cflux_herbivory_LG(1:grid%nland,1:surface%ntiles),       &
                                    cbalance%Cflux_herbivory_2_atm(1:grid%nland,1:surface%ntiles),    &
                                    cbalance%NPP_act_yDayMean(1:grid%nland,1:surface%ntiles),         & 
  frac_litter_wood_new      =       cbalance%frac_litter_wood_new(1:grid%nland,1:surface%ntiles),     &
  redFact_Nlimit            =       nbalance%redFact_Nlimitation(1:grid%nland,1:surface%ntiles),      &
  Npool_green               =       nbalance%Npool_green(1:grid%nland,1:surface%ntiles),              &
  Npool_woods               =       nbalance%Npool_woods(1:grid%nland,1:surface%ntiles),              &
  Npool_mobile              =       nbalance%Npool_mobile(1:grid%nland,1:surface%ntiles),             &
  Npool_litter_green_ag     =       nbalance%Npool_litter_green_ag(1:grid%nland,1:surface%ntiles),    &
  Npool_litter_green_bg     =       nbalance%Npool_litter_green_bg(1:grid%nland,1:surface%ntiles),    &
  Npool_litter_wood_ag      =       nbalance%Npool_litter_wood_ag(1:grid%nland,1:surface%ntiles),     &
  Npool_litter_wood_bg      =       nbalance%Npool_litter_wood_bg(1:grid%nland,1:surface%ntiles),     &
  Npool_slow                =       nbalance%Npool_slow(1:grid%nland,1:surface%ntiles),               &
  SMINN_pool                =       nbalance%SMINN_pool(1:grid%nland,1:surface%ntiles),               &
  minNflux_litter_green_ag  =       nbalance%minNflux_litter_green_ag(1:grid%nland,1:surface%ntiles), &
  minNflux_litter_green_bg  =       nbalance%minNflux_litter_green_bg(1:grid%nland,1:surface%ntiles), &
  minNflux_litter_wood_ag   =       nbalance%minNflux_litter_wood_ag(1:grid%nland,1:surface%ntiles),  &
  minNflux_litter_wood_bg   =       nbalance%minNflux_litter_wood_bg(1:grid%nland,1:surface%ntiles),  &
  minNflux_slow             =       nbalance%minNflux_slow(1:grid%nland,1:surface%ntiles),            &
  Nplant_demand             =       nbalance%Nplant_demand(1:grid%nland,1:surface%ntiles),            &
  Nsoil_demand              =       nbalance%Nsoil_demand(1:grid%nland,1:surface%ntiles),             &
  Ntotal_demand             =       nbalance%Ntotal_demand(1:grid%nland,1:surface%ntiles),            &
  SMINN_herbivory           =       nbalance%SMINN_herbivory(1:grid%nland,1:surface%ntiles),          &
  N2O_emissions_mineraliz   =       nbalance%N2O_emissions_mineraliz(1:grid%nland,1:surface%ntiles),  &
  N2O_emissions_slow        =       nbalance%N2O_emissions_slow(1:grid%nland,1:surface%ntiles),       &
  N2O_emissions_grazing     =       nbalance%N2O_emissions_grazing(1:grid%nland,1:surface%ntiles) ) 
 
              ELSE !! Carbon-only mode
                 IF (with_yasso) THEN
                    CALL update_Cpools(cbalance_diag%LAI_yDaymean(1:grid%nland,1:surface%ntiles,day),         &
                                       cbalance_diag%LAI_previousDayMean(1:grid%nland,1:surface%ntiles),      &
                                       cbalance_diag%NPP_yDayMean(1:grid%nland,1:surface%ntiles,day),         &
                                       cbalance_diag%topSoilTemp_yDayMean(1:grid%nland,1:surface%ntiles,day), &
                                       cbalance_diag%alpha_yDayMean(1:grid%nland,1:surface%ntiles,day),       &
                                       frac_npp_2_woodPool(1:grid%nland,1:surface%ntiles),                    &
                                       frac_npp_2_reservePool(1:grid%nland,1:surface%ntiles),                 &
                                       frac_npp_2_exudates(1:grid%nland,1:surface%ntiles),                    &
                                       frac_green_2_herbivory(1:grid%nland,1:surface%ntiles),                 &
                                       tau_Cpool_litter_green(1:grid%nland,1:surface%ntiles),                 &
                                       tau_Cpool_litter_wood(1:grid%nland,1:surface%ntiles),                  &
                                       tau_Cpool_woods(1:grid%nland,1:surface%ntiles),                        &
                                       LAI_shed_constant(1:grid%nland,1:surface%ntiles),                      &
                                       frac_C_litter_green2atmos(1:grid%nland,1:surface%ntiles),              &
                                       Max_C_content_woods(1:grid%nland,1:surface%ntiles),                    &
                                       specificLeafArea_C(1:grid%nland,1:surface%ntiles),                     &
                                       reserveC2leafC(1:grid%nland,1:surface%ntiles),                         &
                                       Max_LAI(1:grid%nland,1:surface%ntiles),                                &
                                       surface%is_vegetation(1:grid%nland,1:surface%ntiles),                  &
                                       cbalance%Cpool_green(1:grid%nland,1:surface%ntiles),                   &
                                       cbalance%Cpool_woods(1:grid%nland,1:surface%ntiles),                   &
                                       cbalance%Cpool_reserve(1:grid%nland,1:surface%ntiles),                 &
                                       cbalance%Cpool_litter_green_ag(1:grid%nland,1:surface%ntiles),         &
                                       cbalance%Cpool_litter_green_bg(1:grid%nland,1:surface%ntiles),         &
                                       cbalance%Cpool_litter_wood_ag(1:grid%nland,1:surface%ntiles),          &
                                       cbalance%Cpool_litter_wood_bg(1:grid%nland,1:surface%ntiles),          &
                                       cbalance%Cpool_slow(1:grid%nland,1:surface%ntiles),                    &
                                       cbalance%soil_respiration(1:grid%nland,1:surface%ntiles),              &
                                       cbalance%soil_respiration_pot(1:grid%nland,1:surface%ntiles),          &
                                       cbalance%NPP_flux_correction(1:grid%nland,1:surface%ntiles),           &
                                       cbalance%excess_NPP(1:grid%nland,1:surface%ntiles),                    &
                                       cbalance%root_exudates(1:grid%nland,1:surface%ntiles),                 &
                                       litter_flux(1:grid%nland,1:surface%ntiles),                            &
                                       cbalance%Cflux_herbivory(1:grid%nland,1:surface%ntiles),               &
                                       cbalance%Cflux_herbivory_LG(1:grid%nland,1:surface%ntiles),            &
                                       cbalance%Cflux_herbivory_2_atm(1:grid%nland,1:surface%ntiles),         &
                                       cbalance%NPP_act_yDayMean(1:grid%nland,1:surface%ntiles),              &
     frac_litter_wood_new    =         cbalance%frac_litter_wood_new(1:grid%nland,1:surface%ntiles),          &
     temp2_30d               =         cbalance_diag%pseudo_temp_yDay(1:grid%nland,1:surface%ntiles,day),     &
     precip_30d              =         cbalance_diag%pseudo_precip_yDay(1:grid%nland,1:surface%ntiles,day),   &
     YCpool_acid_ag1         =         cbalance%YCpool_acid_ag1(1:grid%nland,1:surface%ntiles),               &
     YCpool_water_ag1        =         cbalance%YCpool_water_ag1(1:grid%nland,1:surface%ntiles),              &
     YCpool_ethanol_ag1      =         cbalance%YCpool_ethanol_ag1(1:grid%nland,1:surface%ntiles),            &
     YCpool_nonsoluble_ag1   =         cbalance%YCpool_nonsoluble_ag1(1:grid%nland,1:surface%ntiles),         &
     YCpool_acid_bg1         =         cbalance%YCpool_acid_bg1(1:grid%nland,1:surface%ntiles),               &
     YCpool_water_bg1        =         cbalance%YCpool_water_bg1(1:grid%nland,1:surface%ntiles),              &
     YCpool_ethanol_bg1      =         cbalance%YCpool_ethanol_bg1(1:grid%nland,1:surface%ntiles),            &
     YCpool_nonsoluble_bg1   =         cbalance%YCpool_nonsoluble_bg1(1:grid%nland,1:surface%ntiles),         &
     YCpool_humus_1          =         cbalance%YCpool_humus_1(1:grid%nland,1:surface%ntiles),                &
     YCpool_acid_ag2         =         cbalance%YCpool_acid_ag2(1:grid%nland,1:surface%ntiles),               &
     YCpool_water_ag2        =         cbalance%YCpool_water_ag2(1:grid%nland,1:surface%ntiles),              &
     YCpool_ethanol_ag2      =         cbalance%YCpool_ethanol_ag2(1:grid%nland,1:surface%ntiles),            &
     YCpool_nonsoluble_ag2   =         cbalance%YCpool_nonsoluble_ag2(1:grid%nland,1:surface%ntiles),         &
     YCpool_acid_bg2         =         cbalance%YCpool_acid_bg2(1:grid%nland,1:surface%ntiles),               &
     YCpool_water_bg2        =         cbalance%YCpool_water_bg2(1:grid%nland,1:surface%ntiles),              &
     YCpool_ethanol_bg2      =         cbalance%YCpool_ethanol_bg2(1:grid%nland,1:surface%ntiles),            &
     YCpool_nonsoluble_bg2   =         cbalance%YCpool_nonsoluble_bg2(1:grid%nland,1:surface%ntiles),         &
     YCpool_humus_2          =         cbalance%YCpool_humus_2(1:grid%nland,1:surface%ntiles),                &
     LeafLit_coef1           =         LeafLit_coef(1:grid%nland,1:surface%ntiles,1),                         &
     LeafLit_coef2           =         LeafLit_coef(1:grid%nland,1:surface%ntiles,2),                         &
     LeafLit_coef3           =         LeafLit_coef(1:grid%nland,1:surface%ntiles,3),                         &
     LeafLit_coef4           =         LeafLit_coef(1:grid%nland,1:surface%ntiles,4),                         &
     LeafLit_coef5           =         LeafLit_coef(1:grid%nland,1:surface%ntiles,5),                         &
     WoodLit_coef1           =         WoodLit_coef(1:grid%nland,1:surface%ntiles,1),                         &
     WoodLit_coef2           =         WoodLit_coef(1:grid%nland,1:surface%ntiles,2),                         &
     WoodLit_coef3           =         WoodLit_coef(1:grid%nland,1:surface%ntiles,3),                         &
     WoodLit_coef4           =         WoodLit_coef(1:grid%nland,1:surface%ntiles,4),                         &
     WoodLit_coef5           =         WoodLit_coef(1:grid%nland,1:surface%ntiles,5),                         &
     WoodLitterSize          =         WoodLitterSize(1:grid%nland,1:surface%ntiles)                          &
                                    )
                 ELSE
                    CALL update_Cpools(cbalance_diag%LAI_yDaymean(1:grid%nland,1:surface%ntiles,day),         &
                                    cbalance_diag%LAI_previousDayMean(1:grid%nland,1:surface%ntiles),         &
                                    cbalance_diag%NPP_yDayMean(1:grid%nland,1:surface%ntiles,day),            &
                                    cbalance_diag%topSoilTemp_yDayMean(1:grid%nland,1:surface%ntiles,day),    &
                                    cbalance_diag%alpha_yDayMean(1:grid%nland,1:surface%ntiles,day),          &
                                    frac_npp_2_woodPool(1:grid%nland,1:surface%ntiles),                       &
                                    frac_npp_2_reservePool(1:grid%nland,1:surface%ntiles),                    &
                                    frac_npp_2_exudates(1:grid%nland,1:surface%ntiles),                       &
                                    frac_green_2_herbivory(1:grid%nland,1:surface%ntiles),                    &
                                    tau_Cpool_litter_green(1:grid%nland,1:surface%ntiles),                    &
                                    tau_Cpool_litter_wood(1:grid%nland,1:surface%ntiles),                     &
                                    tau_Cpool_woods(1:grid%nland,1:surface%ntiles),                           &
                                    LAI_shed_constant(1:grid%nland,1:surface%ntiles),                         &
                                    frac_C_litter_green2atmos(1:grid%nland,1:surface%ntiles),                 &
                                    Max_C_content_woods(1:grid%nland,1:surface%ntiles),                       &
                                    specificLeafArea_C(1:grid%nland,1:surface%ntiles),                        &
                                    reserveC2leafC(1:grid%nland,1:surface%ntiles),                            &
                                    Max_LAI(1:grid%nland,1:surface%ntiles),                                   &
                                    surface%is_vegetation(1:grid%nland,1:surface%ntiles),                     &
                                    cbalance%Cpool_green(1:grid%nland,1:surface%ntiles),                      &
                                    cbalance%Cpool_woods(1:grid%nland,1:surface%ntiles),                      &
                                    cbalance%Cpool_reserve(1:grid%nland,1:surface%ntiles),                    &
                                    cbalance%Cpool_litter_green_ag(1:grid%nland,1:surface%ntiles),            &
                                    cbalance%Cpool_litter_green_bg(1:grid%nland,1:surface%ntiles),            &
                                    cbalance%Cpool_litter_wood_ag(1:grid%nland,1:surface%ntiles),             &
                                    cbalance%Cpool_litter_wood_bg(1:grid%nland,1:surface%ntiles),             &
                                    cbalance%Cpool_slow(1:grid%nland,1:surface%ntiles),                       &
                                    cbalance%soil_respiration(1:grid%nland,1:surface%ntiles),                 &
                                    cbalance%soil_respiration_pot(1:grid%nland,1:surface%ntiles),             &
                                    cbalance%NPP_flux_correction(1:grid%nland,1:surface%ntiles),              &
                                    cbalance%excess_NPP(1:grid%nland,1:surface%ntiles),                       &
                                    cbalance%root_exudates(1:grid%nland,1:surface%ntiles),                    &
                                    litter_flux(1:grid%nland,1:surface%ntiles),                               &
                                    cbalance%Cflux_herbivory(1:grid%nland,1:surface%ntiles),                  &
                                    cbalance%Cflux_herbivory_LG(1:grid%nland,1:surface%ntiles),               &
                                    cbalance%Cflux_herbivory_2_atm(1:grid%nland,1:surface%ntiles),            &
                                    cbalance%NPP_act_yDayMean(1:grid%nland,1:surface%ntiles),                 &
          frac_litter_wood_new    = cbalance%frac_litter_wood_new(1:grid%nland,1:surface%ntiles)              &
                                    )
                 END IF
              END IF

              !! Update Carbon and Nitrogen conservation test 
              !!    (should be zero if everything is correct, except for green litter pool)
              IF (.NOT. with_yasso) THEN
                 cbalance_diag%testCconserv(:,:) =                                           &
                     (cbalance_diag%testCconserv(:,:)                                        &  !! Sum of pools at beginning of day 
                     + cbalance%NPP_act_yDayMean(:,:)* 86400._dp                             &  !! ... plus NPP added 
                     + cbalance%soil_respiration(:,:)* 86400._dp                             &  !! ... minus Carbon lost to atmosph.
                     + cbalance%Cflux_herbivory_2_atm(:,:)* 86400._dp                        &  !!              (times s/day)
                     -   (  cbalance%Cpool_green(:,:)       + cbalance%Cpool_woods(:,:)      &  !! ... minus the pools at end of day
                          + cbalance%Cpool_litter_green_ag(:,:)                              &  !! ... should be very close to zero!
                          + cbalance%Cpool_litter_wood_ag(:,:) + cbalance%Cpool_litter_wood_bg(:,:) &
                          + cbalance%Cpool_litter_green_bg(:,:) + cbalance%Cpool_slow(:,:)          &
                          + cbalance%Cpool_reserve(:,:) )) &
                     / MAX(1._dp, cbalance_diag%testCconserv(:,:))       !! Get relative error: With dynveg the (non area weighted)
                                                                         !! C-pools span several orders of magnitudes, and so do the
                                                                         !! errors.  
              ELSE
                 cbalance_diag%testCconserv(:,:) =                                           &
                     (cbalance_diag%testCconserv(:,:)                                        &  !! Sum of pools at beginning of day 
                     + cbalance%NPP_act_yDayMean(:,:)* 86400._dp                             &  !! ... plus NPP added 
                     + cbalance%soil_respiration(:,:)* 86400._dp                             &  !! ... minus Carbon lost to atmosph.
                     + cbalance%Cflux_herbivory_2_atm(:,:)* 86400._dp                        &  !!              (times s/day)
                     -   (   cbalance%Cpool_green(:,:)                                       &
                           + cbalance%Cpool_woods(:,:)                                       &
                           + cbalance%Cpool_reserve(:,:)                                     &
                           + cbalance%YCpool_acid_ag1(:,:)                                   &
                           + cbalance%YCpool_water_ag1(:,:)                                  &
                           + cbalance%YCpool_ethanol_ag1(:,:)                                &
                           + cbalance%YCpool_nonsoluble_ag1(:,:)                             &
                           + cbalance%YCpool_acid_bg1(:,:)                                   &
                           + cbalance%YCpool_water_bg1(:,:)                                  &
                           + cbalance%YCpool_ethanol_bg1(:,:)                                &
                           + cbalance%YCpool_nonsoluble_bg1(:,:)                             &
                           + cbalance%YCpool_humus_1(:,:)                                    &
                           + cbalance%YCpool_acid_ag2(:,:)                                   &
                           + cbalance%YCpool_water_ag2(:,:)                                  &
                           + cbalance%YCpool_ethanol_ag2(:,:)                                &
                           + cbalance%YCpool_nonsoluble_ag2(:,:)                             &
                           + cbalance%YCpool_acid_bg2(:,:)                                   &
                           + cbalance%YCpool_water_bg2(:,:)                                  &
                           + cbalance%YCpool_ethanol_bg2(:,:)                                &
                           + cbalance%YCpool_nonsoluble_bg2(:,:)                             &
                           + cbalance%YCpool_humus_2(:,:)         ))                         &                                  






                     / MAX(1._dp, cbalance_diag%testCconserv(:,:))       !! Get relative error: With dynveg the (non area weighted)
                                                                         !! C-pools span several orders of magnitudes, and so do the
                                                                         !! errors.
              END IF

              IF (with_nitrogen) THEN
                 nbalance_offline%testNconserv(:,:) =                                        &
                      nbalance_offline%testNconserv(:,:)                                     &
                      -   (  nbalance%Npool_green(:,:)          + nbalance%Npool_woods(:,:)  & !! Note: reserve pool contains no N!
                           + nbalance%Npool_litter_green_ag(:,:)                             &
                           + nbalance%Npool_litter_wood_ag(:,:) + nbalance%Npool_litter_wood_bg(:,:) &
                           + nbalance%Npool_mobile(:,:)                                              &
                           + nbalance%Npool_litter_green_bg(:,:) + nbalance%Npool_slow(:,:)          &
                           + nbalance%SMINN_pool(:,:) - nbalance%NetEcosyst_N_flux(:,:)* 86400._dp )

                 !! Filling arrays for for testing presevation C/N-ratios (should be zero except for green litter pools)

                 nbalance_offline%test_Ngreen(:,:)       = nbalance%Npool_green(:,:)*cn_green - cbalance%Cpool_green(:,:)
                 nbalance_offline%test_Nwoods(:,:)       = nbalance%Npool_woods(:,:)*cn_woods - cbalance%Cpool_woods(:,:)
                 nbalance_offline%test_Nlitter_green_ag(:,:)= &
                      nbalance%Npool_litter_green_ag(:,:)*cn_litter_green - cbalance%Cpool_litter_green_ag(:,:)
                 nbalance_offline%test_Nlitter_green_bg(:,:)= &
                      nbalance%Npool_litter_green_bg(:,:)*cn_litter_green - cbalance%Cpool_litter_green_bg(:,:)
                 nbalance_offline%test_Nlitter_wood_ag(:,:) = &
                      nbalance%Npool_litter_wood_ag(:,:)*cn_litter_wood - cbalance%Cpool_litter_wood_ag(:,:)
                 nbalance_offline%test_Nlitter_wood_bg(:,:) = &
                      nbalance%Npool_litter_wood_bg(:,:)*cn_litter_wood - cbalance%Cpool_litter_wood_bg(:,:)
                 nbalance_offline%test_Nslow(:,:)        = nbalance%Npool_slow(:,:)*cn_slow - cbalance%Cpool_slow(:,:)
              END IF

              !! remember LAI of second last day
              
              IF (day == nday) THEN
                cbalance_diag%LAI_previousDayMean(1:grid%nland,1:surface%ntiles) = &
                  cbalance_diag%LAI_yDaymean(1:grid%nland,1:surface%ntiles,day)                 
              ENDIF

              !! further diagnostics

              cbalance_diag%avg_Cpool_green(:,:)           = cbalance_diag%avg_Cpool_green(:,:) &
                                                            + coverFractBox(:,:)*cbalance%Cpool_green(:,:)
              cbalance_diag%avg_Cpool_woods(:,:)           = cbalance_diag%avg_Cpool_woods(:,:) &
                                                            + coverFractBox(:,:)*cbalance%Cpool_woods(:,:)
              cbalance_diag%avg_Cpool_reserve(:,:)         = cbalance_diag%avg_Cpool_reserve(:,:) &
                                                            + coverFractBox(:,:)*cbalance%Cpool_reserve(:,:)
              cbalance_diag%avg_Cpool_litter_green_ag(:,:) = cbalance_diag%avg_Cpool_litter_green_ag(:,:) &
                                                            + coverFractBox(:,:)*cbalance%Cpool_litter_green_ag(:,:)
              cbalance_diag%avg_Cpool_litter_green_bg(:,:) = cbalance_diag%avg_Cpool_litter_green_bg(:,:) &
                                                            + coverFractBox(:,:)*cbalance%Cpool_litter_green_bg(:,:)
              cbalance_diag%avg_Cpool_litter_wood_ag(:,:)  = cbalance_diag%avg_Cpool_litter_wood_ag(:,:) &
                                                            + coverFractBox(:,:)*cbalance%Cpool_litter_wood_ag(:,:)
              cbalance_diag%avg_Cpool_litter_wood_bg(:,:)  = cbalance_diag%avg_Cpool_litter_wood_bg(:,:) &
                                                            + coverFractBox(:,:)*cbalance%Cpool_litter_wood_bg(:,:)
              cbalance_diag%avg_Cpool_slow(:,:)            = cbalance_diag%avg_Cpool_slow(:,:) &
                                                            + coverFractBox(:,:)*cbalance%Cpool_slow(:,:)
              IF (with_yasso) THEN
                 cbalance_diag%avg_YCpool_acid_ag1(:,:)        = cbalance_diag%avg_YCpool_acid_ag1(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_acid_ag1(:,:)
                 cbalance_diag%avg_YCpool_water_ag1(:,:)       = cbalance_diag%avg_YCpool_water_ag1(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_water_ag1(:,:)
                 cbalance_diag%avg_YCpool_ethanol_ag1(:,:)     = cbalance_diag%avg_YCpool_ethanol_ag1(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_ethanol_ag1(:,:)
                 cbalance_diag%avg_YCpool_nonsoluble_ag1(:,:)  = cbalance_diag%avg_YCpool_nonsoluble_ag1(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_nonsoluble_ag1(:,:)
                 cbalance_diag%avg_YCpool_acid_bg1(:,:)        = cbalance_diag%avg_YCpool_acid_bg1(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_acid_bg1(:,:)
                 cbalance_diag%avg_YCpool_water_bg1(:,:)       = cbalance_diag%avg_YCpool_water_bg1(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_water_bg1(:,:)
                 cbalance_diag%avg_YCpool_ethanol_bg1(:,:)     = cbalance_diag%avg_YCpool_ethanol_bg1(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_ethanol_bg1(:,:)
                 cbalance_diag%avg_YCpool_nonsoluble_bg1(:,:)  = cbalance_diag%avg_YCpool_nonsoluble_bg1(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_nonsoluble_bg1(:,:)
                 cbalance_diag%avg_YCpool_humus_1(:,:)         =cbalance_diag%avg_YCpool_humus_1(:,:) & 
                                                               + coverFractBox(:,:)*cbalance%YCpool_humus_1(:,:)
                 cbalance_diag%avg_YCpool_acid_ag2(:,:)        = cbalance_diag%avg_YCpool_acid_ag2(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_acid_ag2(:,:)
                 cbalance_diag%avg_YCpool_water_ag2(:,:)       = cbalance_diag%avg_YCpool_water_ag2(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_water_ag2(:,:)
                 cbalance_diag%avg_YCpool_ethanol_ag2(:,:)     = cbalance_diag%avg_YCpool_ethanol_ag2(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_ethanol_ag2(:,:)
                 cbalance_diag%avg_YCpool_nonsoluble_ag2(:,:)  = cbalance_diag%avg_YCpool_nonsoluble_ag2(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_nonsoluble_ag2(:,:)
                 cbalance_diag%avg_YCpool_acid_bg2(:,:)        = cbalance_diag%avg_YCpool_acid_bg2(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_acid_bg2(:,:)
                 cbalance_diag%avg_YCpool_water_bg2(:,:)       = cbalance_diag%avg_YCpool_water_bg2(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_water_bg2(:,:)
                 cbalance_diag%avg_YCpool_ethanol_bg2(:,:)     = cbalance_diag%avg_YCpool_ethanol_bg2(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_ethanol_bg2(:,:)
                 cbalance_diag%avg_YCpool_nonsoluble_bg2(:,:)  = cbalance_diag%avg_YCpool_nonsoluble_bg2(:,:) &
                                                               + coverFractBox(:,:)*cbalance%YCpool_nonsoluble_bg2(:,:)
                 cbalance_diag%avg_YCpool_humus_2(:,:)         =cbalance_diag%avg_YCpool_humus_2(:,:) & 
                                                               + coverFractBox(:,:)*cbalance%YCpool_humus_2(:,:)
              END IF
              IF (use_external_landcover_maps .OR. with_landuse_transitions) THEN
                IF (lcc_scheme==2) THEN
                  cbalance_diag%avg_Cpool_onSite      (:) = &
                  cbalance_diag%avg_Cpool_onSite      (:) + surface%veg_ratio_max(:) * cbalance%Cpool_onSite(:)
                  cbalance_diag%avg_Cpool_paper       (:) = &
                  cbalance_diag%avg_Cpool_paper       (:) + surface%veg_ratio_max(:) * cbalance%Cpool_paper       (:)
                  cbalance_diag%avg_Cpool_construction(:) = &
                  cbalance_diag%avg_Cpool_construction(:) + surface%veg_ratio_max(:) * cbalance%Cpool_construction(:)
                  IF (with_landuse_transitions) THEN
                    cbalance_diag%avg_Cpool_paper_harvest       (:) = &
                    cbalance_diag%avg_Cpool_paper_harvest       (:) + surface%veg_ratio_max(:) * cbalance%Cpool_paper_harvest(:)
                    cbalance_diag%avg_Cpool_construction_harvest(:) = &
                    cbalance_diag%avg_Cpool_construction_harvest(:)+surface%veg_ratio_max(:)*cbalance%Cpool_construction_harvest(:)
                  ENDIF
                ENDIF
              ENDIF

#ifndef __PGI
              IF (with_yasso) THEN
                  cbalance_diag%box_Cpools_total(:)    = sum(coverFractBox(:,:)                                                   &
                                                    * (  cbalance%Cpool_green(:,:)                                                &
                                                       + cbalance%Cpool_woods(:,:)                                                &
                                                       + cbalance%YCpool_acid_ag1(:,:)                                            &
                                                       + cbalance%YCpool_water_ag1(:,:)                                           &
                                                       + cbalance%YCpool_ethanol_ag1(:,:)                                         &
                                                       + cbalance%YCpool_nonsoluble_ag1(:,:)                                      &
                                                       + cbalance%YCpool_acid_bg1(:,:)                                            &
                                                       + cbalance%YCpool_water_bg1(:,:)                                           &
                                                       + cbalance%YCpool_ethanol_bg1(:,:)                                         &
                                                       + cbalance%YCpool_nonsoluble_bg1(:,:)                                      &
                                                       + cbalance%YCpool_humus_1(:,:)                                             &
                                                       + cbalance%YCpool_acid_ag2(:,:)                                            &
                                                       + cbalance%YCpool_water_ag2(:,:)                                           &
                                                       + cbalance%YCpool_ethanol_ag2(:,:)                                         &
                                                       + cbalance%YCpool_nonsoluble_ag2(:,:)                                      &
                                                       + cbalance%YCpool_acid_bg2(:,:)                                            &
                                                       + cbalance%YCpool_water_bg2(:,:)                                           &
                                                       + cbalance%YCpool_ethanol_bg2(:,:)                                         &
                                                       + cbalance%YCpool_nonsoluble_bg2(:,:)                                      &
                                                       + cbalance%YCpool_humus_2(:,:)                                             &
                                                       + cbalance%Cpool_reserve(:,:)),DIM=2)
              ELSE
                  cbalance_diag%box_Cpools_total(:)    = sum(coverFractBox(:,:)                                                   &
                                                    * (  cbalance%Cpool_green          (:,:) + cbalance%Cpool_woods(:,:)          &
                                                       + cbalance%Cpool_litter_green_ag(:,:)                                      &
                                                       + cbalance%Cpool_litter_wood_ag (:,:) + cbalance%Cpool_litter_wood_bg(:,:) &
                                                       + cbalance%Cpool_litter_green_bg(:,:) + cbalance%Cpool_slow(:,:)           &
                                                       + cbalance%Cpool_reserve        (:,:)),DIM=2)
              END IF                                     
              
              IF (with_nitrogen) THEN
                 nbalance_offline%box_Npools_total(:) = sum(coverFractBox(:,:)                                                &
                                                * (  nbalance%Npool_green          (:,:) + nbalance%Npool_woods(:,:)          &
                                                   + nbalance%Npool_litter_green_ag(:,:)                                      &
                                                   + nbalance%Npool_litter_wood_ag (:,:) + nbalance%Npool_litter_wood_bg(:,:) &
                                                   + nbalance%Npool_litter_green_bg(:,:) + nbalance%Npool_slow(:,:)           &
                                                   + nbalance%Npool_mobile         (:,:) + nbalance%SMINN_pool(:,:)           &
                                                  ),DIM=2)
              END IF
#else
              !
              ! Different notation to allow compilation with the PGI compiler (pgf95 6.1-1)
              !
              IF (with_yasso) THEN
                 cbalance_diag%box_Cpools_total(:) = 0._dp
                 DO itile=1,surface%ntiles
                    cbalance_diag%box_Cpools_total(:) = cbalance_diag%box_Cpools_total(:) + coverFractBox(:,itile)                 &
                                               * (  cbalance%Cpool_green(:,itile)                                                  &
                                                       + cbalance%Cpool_woods(:,itile)                                             &
                                                       + cbalance%Cpool_reserve(:,itile)                                           &
                                                       + cbalance%YCpool_acid_ag1(:,itile)                                         &
                                                       + cbalance%YCpool_water_ag1(:,itile)                                        &
                                                       + cbalance%YCpool_ethanol_ag1(:,itile)                                      &
                                                       + cbalance%YCpool_nonsoluble_ag1(:,itile)                                   &
                                                       + cbalance%YCpool_acid_bg1(:,itile)                                         &
                                                       + cbalance%YCpool_water_bg1(:,itile)                                        &
                                                       + cbalance%YCpool_ethanol_bg1(:,itile)                                      &
                                                       + cbalance%YCpool_nonsoluble_bg1(:,itile)                                   &
                                                       + cbalance%YCpool_humus_1(:,itile)                                          &
                                                       + cbalance%YCpool_acid_ag2(:,itile)                                         &
                                                       + cbalance%YCpool_water_ag2(:,itile)                                        &
                                                       + cbalance%YCpool_ethanol_ag2(:,itile)                                      &
                                                       + cbalance%YCpool_nonsoluble_ag2(:,itile)                                   &
                                                       + cbalance%YCpool_acid_bg2(:,itile)                                         &
                                                       + cbalance%YCpool_water_bg2(:,itile)                                        &
                                                       + cbalance%YCpool_ethanol_bg2(:,itile)                                      &
                                                       + cbalance%YCpool_nonsoluble_bg2(:,itile)                                   &
                                                       + cbalance%YCpool_humus_2(:,itile)                                          &
                                                 )
                 END DO
              ELSE
                 cbalance_diag%box_Cpools_total(:) = 0._dp
                 DO itile=1,surface%ntiles
                    cbalance_diag%box_Cpools_total(:) = cbalance_diag%box_Cpools_total(:) + coverFractBox(:,itile) &
                                               * (  cbalance%Cpool_green(:,itile) + cbalance%Cpool_woods(:,itile) &
                                                  + cbalance%Cpool_litter_green_ag(:,itile) &
                                                  + cbalance%Cpool_litter_wood_ag(:,itile) + cbalance%Cpool_litter_wood_bg(:,itile)&
                                                  + cbalance%Cpool_litter_green_bg(:,itile) + cbalance%Cpool_slow(:,itile) &
                                                  + cbalance%Cpool_reserve(:,itile) &
                                                 )
                 END DO
              END IF

              IF (with_nitrogen) THEN
                 nbalance_offline%box_Npools_total(:) = 0._dp
                 DO itile=1,surface%ntiles
                    nbalance_offline%box_Npools_total(:) = nbalance_offline%box_Npools_total(:)+ coverFractBox(:,itile)          &
                                            * ( nbalance%Npool_green(:,itile)        + nbalance%Npool_woods(:,itile)             &
                                               + nbalance%Npool_litter_green_ag(:,itile)                                         &
                                               + nbalance%Npool_litter_wood_ag(:,itile) + nbalance%Npool_litter_wood_bg(:,itile) &
                                               + nbalance%Npool_litter_green_bg(:,itile)         + nbalance%Npool_slow(:,itile)  &
                                               + nbalance%Npool_mobile(:,itile)       + nbalance%SMINN_pool(:,itile)             &
                                              )
                 END DO
              END IF
#endif
              IF ((use_external_landcover_maps .OR. with_landuse_transitions) .AND. (lcc_scheme==2)) THEN 
                cbalance_diag%box_Cpools_total(:) = cbalance_diag%box_Cpools_total(:) + surface%veg_ratio_max(:) &
                                                    * (   cbalance%Cpool_onsite      (:)                         &
                                                        + cbalance%Cpool_paper       (:)                         &
                                                        + cbalance%Cpool_construction(:)                         &
                                                      )
              ENDIF
              IF (with_landuse_transitions .AND. (lcc_scheme==2)) THEN
                cbalance_diag%box_Cpools_total(:) = cbalance_diag%box_Cpools_total(:) + surface%veg_ratio_max(:) &
                                                    * (   cbalance%Cpool_paper_harvest       (:)                 &
                                                        + cbalance%Cpool_construction_harvest(:)                 &
                                                      )
              ENDIF

              !! Testarrays for C and N conservation including landcover change:

              glacier_1d(:) = ANY(surface%is_glacier(:,:),DIM=2) ! direct use of the ANY construct in the WHERE below leads to a
                                                                 ! memory leak on blizzard (~10MB per year in T63)
              WHERE (.NOT. glacier_1d(:))
                 carbon_2_atmos(:) = (CO2_emission_landcover_change(:) + CO2_emission_harvest(:) + CO2_flux_dynveg(:)) &
                                    / surface%veg_ratio_max(:) / molarMassCO2_kg * 86400._dp
              ELSEWHERE
                 carbon_2_atmos(:) = 0._dp
              END WHERE

              cbalance_diag%box_test_Ccons(:) = cbalance_diag%box_test_Ccons(:)                                           &
                                                - surface%veg_ratio_max(:)*carbon_2_atmos(:)                              &
                                                + sum(coverFractBox(:,:)*cbalance%NPP_act_yDayMean(:,:)* 86400._dp,DIM=2) &
                                                + sum(coverFractBox(:,:)*cbalance%soil_respiration(:,:)* 86400._dp,DIM=2) &
                                                + sum(coverFractBox(:,:)*cbalance%Cflux_herbivory_2_atm(:,:)* 86400._dp,DIM=2) &
                                                - cbalance_diag%box_Cpools_total(:)
              IF (with_nitrogen) THEN
                 nbalance_offline%box_test_Ncons(:) = nbalance_offline%box_test_Ncons(:) & ! leaching and fixation are missing by !c
                                                - surface%veg_ratio_max(:)*nitrogen_2_atmos(:) &
                                                - nbalance_offline%box_Npools_total(:)          &
                                                + sum(coverFractBox(:,:)*nbalance%NetEcosyst_N_flux(:,:)* 86400._dp,DIM=2)     &
                                                - sum(coverFractBox(:,:)*nbalance%N2O_emissions_mineraliz(:,:)* 86400._dp,DIM=2)
                                                                                          !! leaching, denit,fix, depo are included 

              END IF

              cbalance_diag%avg_soil_respiration(:,:)      = cbalance_diag%avg_soil_respiration(:,:)                            &
                                                            + coverFractBox(:,:) * cbalance%soil_respiration(:,:)
              cbalance_diag%avg_soil_respiration_pot(:,:)  = cbalance_diag%avg_soil_respiration_pot(:,:)                        &
                                                             + coverFractBox(:,:) * cbalance%soil_respiration_pot(:,:)
              cbalance_diag%avg_NPP_yDayMean(:,:)          = cbalance_diag%avg_NPP_yDayMean(:,:)                                &
                                                            + coverFractBox(:,:) * cbalance_diag%NPP_yDayMean(:,:,day)
              cbalance_diag%avg_NPP_flux_correction(:,:)   = cbalance_diag%avg_NPP_flux_correction(:,:)                         &
                                                            + coverFractBox(:,:) * cbalance%NPP_flux_correction(:,:)
              cbalance_diag%avg_excess_NPP(:,:)            = cbalance_diag%avg_excess_NPP(:,:)                                  &
                                                            + coverFractBox(:,:) * cbalance%excess_NPP(:,:)
              cbalance_diag%avg_root_exudates(:,:)         = cbalance_diag%avg_root_exudates(:,:)                               &
                                                            + coverFractBox(:,:) * cbalance%root_exudates(:,:)
              cbalance_diag%avg_Cflux_herbivory(:,:)       = cbalance_diag%avg_Cflux_herbivory(:,:)                             &
                                                            + coverFractBox(:,:) * cbalance%Cflux_herbivory(:,:)
              cbalance_diag%avg_Cflux_herbivory_LG(:,:)    = cbalance_diag%avg_Cflux_herbivory_LG(:,:)                          &
                                                            + coverFractBox(:,:) * cbalance%Cflux_herbivory_LG(:,:)
              cbalance_diag%avg_Cflux_herbivory_2_atm(:,:) = cbalance_diag%avg_Cflux_herbivory_2_atm(:,:)                       &
                                                            + coverFractBox(:,:) * cbalance%Cflux_herbivory_2_atm(:,:)
              cbalance_diag%box_CO2_flux_2_atm(:)          = cbalance_diag%box_CO2_flux_2_atm(:) + carbon_2_atmos(:)* &
                                                             surface%veg_ratio_max(:) * molarMassCO2_kg / 86400._dp

              cbalance_diag%avg_box_NEP(:) = cbalance_diag%avg_box_NEP(:)                                                        &
                   + sum(coverFractBox(:,:)                                                                                      &
                         * (cbalance_diag%NPP_yDayMean(:,:,day) + cbalance%NPP_flux_correction(:,:) &
                            + cbalance%soil_respiration(:,:) + cbalance%Cflux_herbivory_2_atm(:,:)),DIM=2) * molarMassCO2_kg

              cbalance_diag%avg_NPP_act_yDayMean(:,:)        = cbalance_diag%avg_NPP_act_yDayMean(:,:) &
                                                        + coverFractBox(:,:) * cbalance%NPP_act_yDayMean(:,:)

              IF (with_nitrogen) THEN
                 nbalance_offline%avg_Npool_green(:,:) = &
                      nbalance_offline%avg_Npool_green(:,:) + coverFractBox(:,:)*nbalance%Npool_green(:,:)
                 nbalance_offline%avg_Npool_woods(:,:) = &
                      nbalance_offline%avg_Npool_woods(:,:) + coverFractBox(:,:)*nbalance%Npool_woods(:,:)
                 nbalance_offline%avg_Npool_mobile(:,:) = &
                      nbalance_offline%avg_Npool_mobile(:,:) + coverFractBox(:,:)*nbalance%Npool_mobile(:,:)
                 nbalance_offline%avg_Npool_litter_green_ag(:,:) = &
                      nbalance_offline%avg_Npool_litter_green_ag(:,:) + coverFractBox(:,:)*nbalance%Npool_litter_green_ag(:,:)
                 nbalance_offline%avg_Npool_litter_green_bg(:,:) = &
                      nbalance_offline%avg_Npool_litter_green_bg(:,:) + coverFractBox(:,:)*nbalance%Npool_litter_green_bg(:,:)
                 nbalance_offline%avg_Npool_litter_wood_ag(:,:) = &
                      nbalance_offline%avg_Npool_litter_wood_ag(:,:) + coverFractBox(:,:)*nbalance%Npool_litter_wood_ag(:,:)
                 nbalance_offline%avg_Npool_litter_wood_bg(:,:) = &
                      nbalance_offline%avg_Npool_litter_wood_bg(:,:) + coverFractBox(:,:)*nbalance%Npool_litter_wood_bg(:,:)
                 nbalance_offline%avg_Npool_slow(:,:) = &
                      nbalance_offline%avg_Npool_slow(:,:) + coverFractBox(:,:)*nbalance%Npool_slow(:,:)
                 nbalance_offline%avg_SMINN_pool(:,:) = &
                      nbalance_offline%avg_SMINN_pool(:,:) + coverFractBox(:,:)*nbalance%SMINN_pool(:,:)


                 nbalance_offline%avg_minNflux_litter_green_ag(:,:)  = nbalance_offline%avg_minNflux_litter_green_ag(:,:) + &
                                                           coverFractBox(:,:)*nbalance%minNflux_litter_green_ag(:,:)
                 nbalance_offline%avg_minNflux_litter_green_bg(:,:) = nbalance_offline%avg_minNflux_litter_green_bg(:,:) + &
                                                           coverFractBox(:,:)*nbalance%minNflux_litter_green_bg(:,:)
                 nbalance_offline%avg_minNflux_litter_wood_ag(:,:)= nbalance_offline%avg_minNflux_litter_wood_ag(:,:) + &
                                                            coverFractBox(:,:)*nbalance%minNflux_litter_wood_ag(:,:)
                 nbalance_offline%avg_minNflux_litter_wood_bg(:,:) = nbalance_offline%avg_minNflux_litter_wood_bg(:,:) + &
                                                             coverFractBox(:,:)*nbalance%minNflux_litter_wood_bg(:,:)
                 nbalance_offline%avg_minNflux_slow(:,:) = nbalance_offline%avg_minNflux_slow(:,:)+ &
                                                   coverFractBox(:,:)*nbalance%minNflux_slow(:,:)
                 nbalance_offline%avg_Nplant_demand(:,:) = &
                      nbalance_offline%avg_Nplant_demand(:,:) + coverFractBox(:,:)*nbalance%Nplant_demand(:,:)
                 nbalance_offline%avg_Nsoil_demand(:,:) = &
                      nbalance_offline%avg_Nsoil_demand(:,:) + coverFractBox(:,:)*nbalance%Nsoil_demand(:,:)
                 nbalance_offline%avg_Ntotal_demand(:,:) = &
                      nbalance_offline%avg_Ntotal_demand(:,:) + coverFractBox(:,:)*nbalance%Ntotal_demand(:,:)

                 nbalance_offline%avg_SMINN_herbivory(:,:) = &
                      nbalance_offline%avg_SMINN_herbivory(:,:) + coverFractBox(:,:)*nbalance%SMINN_herbivory(:,:)
                 nbalance_offline%avg_nfix_to_sminn(:,:)           = nbalance_offline%avg_nfix_to_sminn(:,:) & 
                                                             + coverFractBox(:,:)*nbalance%nfix_to_sminn(:,:)
                 nbalance_offline%avg_sminn_to_denit(:,:)          = nbalance_offline%avg_sminn_to_denit(:,:)  &
                                                             + coverFractBox(:,:)*nbalance%sminn_to_denit(:,:) 
                 nbalance_offline%avg_sminn_leach(:,:)             = nbalance_offline%avg_sminn_leach(:,:) &
                                                             + coverFractBox(:,:)*nbalance%sminn_leach(:,:)
                 nbalance_offline%avg_N2O_emissions_ExternalN(:,:) = nbalance_offline%avg_N2O_emissions_ExternalN(:,:) &
                                                             + coverFractBox(:,:)*nbalance%N2O_emissions_depfix(:,:)
                 nbalance_offline%avg_N2O_emissions_mineraliz(:,:) = nbalance_offline%avg_N2O_emissions_mineraliz(:,:) &
                                                             + coverFractBox(:,:)*nbalance%N2O_emissions_mineraliz(:,:)
                 nbalance_offline%avg_NetEcosyst_N_flux(:,:)       = nbalance_offline%avg_NetEcosyst_N_flux(:,:)  &
                                                             + coverFractBox(:,:)*nbalance%NetEcosyst_N_flux(:,:)
                 nbalance_offline%avg_nfert_to_sminn(:,:)    = nbalance_offline%avg_nfert_to_sminn(:,:)  &
                                                             + veg_fract_correction(:,:)                  &
                                                             * spread(surface%veg_ratio_max(:),DIM=2,NCOPIES=surface%ntiles)  &
                                                             * nbalance%nfert_to_sminn(:,:)   !! without coverFract(:,:)
                 nbalance_offline%avg_N2O_emissions_nfert(:,:) = nbalance_offline%avg_N2O_emissions_nfert(:,:) &
                                                             + veg_fract_correction(:,:)                       &
                                                             * spread(surface%veg_ratio_max(:),DIM=2,NCOPIES=surface%ntiles)  &
                                                             * nbalance%N2O_emissions_nfert(:,:)   !! without coverFractBOX(:,:)
                 nbalance_offline%avg_N2O_emissions_grazing(:,:) = nbalance_offline%avg_N2O_emissions_grazing(:,:) &
                                                             + coverFractBox(:,:)*nbalance%N2O_emissions_grazing(:,:)
                 nbalance%sum_N_pools(:,:) = nbalance%Npool_green(:,:)          &  
                                           + nbalance%Npool_woods(:,:)          &
                                           + nbalance%Npool_mobile(:,:)         &
                                           + nbalance%Npool_litter_green_ag(:,:)&
                                           + nbalance%Npool_litter_wood_ag(:,:) &
                                           + nbalance%Npool_litter_wood_bg(:,:) &
                                           + nbalance%Npool_litter_green_bg(:,:)&
                                           + nbalance%Npool_slow(:,:)           &
                                           + nbalance%SMINN_pool(:,:)         
              END IF

              number_of_days_in_year = number_of_days_in_year + 1

              IF (outputInterval == INTERVAL_DAY) THEN
!                 step = REAL(run_year) * 10000._dp + REAL(month) * 100._dp + REAL(day)  ! time step for absolute time
                 CALL averageFields(cbalance_diag, disturbance,fire_TH_diag, 1,landcover_change, &
                      landuse_transitions, lcc_scheme)
                 CALL writeSingleTimeStep(trim(outFile), grid, surface, cbalance, cbalance_diag, nbalance, nbalance_offline, &
                      dynveg, landcover_change, landuse_transitions, with_nitrogen, with_yasso, use_dynveg,&
                      use_disturbance, &
                      dynveg_options%dynveg_feedback, use_external_landcover_maps, with_landuse_transitions, lcc_scheme, &
                      run_year_first, step, finish = (day==nday))
                 CALL nullifyFields(cbalance_diag, disturbance, fire_TH_diag, landcover_change,&
                      landuse_transitions,lcc_scheme)
              END IF

              new_year = .false.
              new_month = .false.
              CALL time_reset  !! update echam time step counter (for echam routines called by cbalone)

           END DO ! end of day loop

           IF (outputInterval == INTERVAL_MONTH) THEN

              CALL averageFields(cbalance_diag, disturbance,fire_TH_diag, nday,landcover_change,&
                   landuse_transitions, lcc_scheme)
!              step = REAL(run_year) * 10000._dp + REAL(month) * 100._dp  ! time step for absolute time
              CALL writeSingleTimeStep(trim(outFile), grid, surface, cbalance, cbalance_diag, nbalance, nbalance_offline, &
                   dynveg, landcover_change, landuse_transitions, with_nitrogen, with_yasso, use_dynveg, &
                   use_disturbance, &
                   dynveg_options%dynveg_feedback, use_external_landcover_maps, with_landuse_transitions, lcc_scheme,&
                   run_year_first, &
                   step, finish = (month==12))
              CALL nullifyFields(cbalance_diag, disturbance, fire_TH_diag, landcover_change,&
                   landuse_transitions, lcc_scheme)

              ! write a restart file with data of Dec. 31st
              !    contains all output variables; all box_ and avg_ variables should be ignored 
              IF (run_year == run_year_last .AND. month == 12) THEN
                 CALL writeSingleTimeStep('restart_'//trim(outFile), grid, surface, cbalance, cbalance_diag, nbalance, &
                      nbalance_offline, dynveg, landcover_change, landuse_transitions, with_nitrogen,&
                      with_yasso, use_dynveg, &
                      use_disturbance, dynveg_options%dynveg_feedback, use_external_landcover_maps, &
                      with_landuse_transitions, lcc_scheme, run_year_first, step, restart=.true.)
              END IF

           END IF

        END DO ! end month loop

        IF (outputInterval == INTERVAL_YEAR) THEN

          CALL averageFields(cbalance_diag, disturbance, fire_TH_diag, number_of_days_in_year,landcover_change,&
               landuse_transitions, lcc_scheme)
!          step = REAL(run_year) * 10000._dp  ! time step for absolute time
          CALL writeSingleTimeStep(trim(outFile), grid, surface, cbalance, cbalance_diag, nbalance, nbalance_offline, dynveg, &
               landcover_change, landuse_transitions, with_nitrogen, with_yasso, use_dynveg, use_disturbance, &
               dynveg_options%dynveg_feedback, use_external_landcover_maps, with_landuse_transitions, lcc_scheme,  &
               run_year_first, step, finish = .true.)
          CALL nullifyFields(cbalance_diag, disturbance, fire_TH_diag, landcover_change,&
               landuse_transitions, lcc_scheme)

           ! write a restart file with data of Dec. 31st
           !    contains all output variables; all box_ and avg_ variables should be ignored 
           IF (run_year == run_year_last) THEN
              CALL writeSingleTimeStep('restart_'//trim(outFile), grid, surface, cbalance, cbalance_diag, nbalance, &
                   nbalance_offline, dynveg, landcover_change, landuse_transitions, with_nitrogen, with_yasso, use_dynveg, &
                   use_disturbance, dynveg_options%dynveg_feedback, use_external_landcover_maps, &
                   with_landuse_transitions, lcc_scheme, run_year_first, step, restart=.true.)
           END IF
       END IF

       !! Check for end of run

       IF (run_year == run_year_last) THEN
          EXIT RepeatLoop
       END IF

       !! Force writing to stdout

       CALL flush(6)

     END DO               ! end loop over all climate data years
     
  END DO RepeatLoop                            ! end repeat loop

  CALL InputClose()

  !! deallocate local memory

  DEALLOCATE(frac_npp_2_woodPool, frac_npp_2_reservePool, frac_npp_2_exudates, frac_green_2_herbivory, frac_C_litter_green2atmos)
  DEALLOCATE(tau_Cpool_litter_green, tau_Cpool_litter_wood, tau_Cpool_woods)
  DEALLOCATE(LAI_shed_constant, Max_C_content_woods)
  DEALLOCATE(specificLeafArea_C, reserveC2leafC, veg_fract_correction, coverFractBox)
  DEALLOCATE(CO2_emission_landcover_change, CO2_emission_harvest, CO2_flux_dynveg)
  DEALLOCATE(carbon_2_atmos, carbon_2_litterGreenPools, carbon_2_litterWoodPool_ag, carbon_2_litterWoodPool_bg)
  DEALLOCATE(coverFract_interpolated, litter_flux, glacier_1d)
  DEALLOCATE(LeafLit_coef, WoodLit_coef, WoodLitterSize)

  IF (with_nitrogen) THEN
     DEALLOCATE(Nitrogen_2_atmos)
     DEALLOCATE(nitrogen_2_atmos, nitrogen_2_litterGreenPools, nitrogen_2_litterWoodPool_ag)
     DEALLOCATE(nitrogen_2_litterWoodPool_bg, nitrogen_2_SMINNpool)
  END IF

  !! --- finalization of MPI
  CALL p_stop

  WRITE(*,*) "Normal End!"
  
CONTAINS

  SUBROUTINE averageFields(cbalance_diag, disturbance, fire_TH_diag, noOfTimeSteps,landcover_change, &
         landuse_transitions, lcc_scheme)
    USE mo_cbalone_memory,        ONLY: cbal_offline_type
    USE mo_cbal_landcover_change, ONLY: landcover_change_type
    USE mo_disturbance,           ONLY: dist_opts
    USE mo_disturbance_thonicke,  ONLY: fire_TH_diag_type
    TYPE(cbal_offline_type),intent(inout)     :: cbalance_diag
    TYPE(disturbance_type),intent(inout)      :: disturbance
    TYPE(fire_TH_diag_type), intent(inout)    :: fire_TH_diag
    integer,intent(in) :: noOfTimeSteps
    TYPE(landcover_change_type),   intent(inout) :: landcover_change 
    TYPE(landuse_transitions_type),intent(inout) :: landuse_transitions

    INTEGER,                       intent(in) :: lcc_scheme

    INTEGER, PARAMETER :: scdsPerDay = 86400

    REAL(dp)           :: OutPer

    OutPer = REAL(noOfTimeSteps,dp) * REAL(scdsPerDay,dp)

    IF (ASSOCIATED(disturbance%box_CO2_flux_2_atmos)) &
      disturbance%box_CO2_flux_2_atmos(:) = disturbance%box_CO2_flux_2_atmos(:) / OutPer
    IF (ASSOCIATED(disturbance%box_burned_frac_avg)) &
      disturbance%box_burned_frac_avg(:,:) = disturbance%box_burned_frac_avg(:,:) / OutPer
    IF (            ASSOCIATED(disturbance%box_burned_frac_diag_avg)                                  &
        .AND. .NOT. ASSOCIATED(disturbance%box_burned_frac_diag_avg,disturbance%box_burned_frac_avg)) &
      disturbance%box_burned_frac_diag_avg(:,:) = disturbance%box_burned_frac_diag_avg(:,:) / OutPer
    IF (ASSOCIATED(disturbance%box_damaged_frac_avg)) &
      disturbance%box_damaged_frac_avg(:,:) = disturbance%box_damaged_frac_avg(:,:) / OutPer

    IF (dist_opts%fire_algorithm==3) THEN
      IF (dist_opts%ldiag) THEN
        fire_TH_diag%avg_FDI               (:)   = fire_TH_diag%avg_FDI                (:)   / OutPer
        fire_TH_diag%avg_ROS               (:)   = fire_TH_diag%avg_ROS                (:)   / OutPer
        fire_TH_diag%avg_fuel_1hr_total    (:)   = fire_TH_diag%avg_fuel_1hr_total     (:)   / OutPer
        fire_TH_diag%avg_fuel_10hr_total   (:)   = fire_TH_diag%avg_fuel_10hr_total    (:)   / OutPer
        fire_TH_diag%avg_fuel_100hr_total  (:)   = fire_TH_diag%avg_fuel_100hr_total   (:)   / OutPer
        fire_TH_diag%avg_numfire           (:)   = fire_TH_diag%avg_numfire            (:)   / OutPer
        fire_TH_diag%avg_HumanIgnition     (:)   = fire_TH_diag%avg_HumanIgnition      (:)   / OutPer 
        fire_TH_diag%avg_LightningIgnition (:)   = fire_TH_diag%avg_LightningIgnition  (:)   / OutPer 
        fire_TH_diag%avg_FireDuration      (:)   = fire_TH_diag%avg_FireDuration       (:)   / OutPer 
        fire_TH_diag%avg_fire_intensity    (:)   = fire_TH_diag%avg_fire_intensity     (:)   / OutPer 
        fire_TH_diag%avg_postfire_mortality(:,:) = fire_TH_diag%avg_postfire_mortality (:,:) / OutPer 
        fire_TH_diag%avg_vegetation_height (:,:) = fire_TH_diag%avg_vegetation_height  (:,:) / OutPer
        fire_TH_diag%avg_carbon_2_atmos    (:,:) = fire_TH_diag%avg_carbon_2_atmos     (:,:) / OutPer
      ENDIF
      IF (ASSOCIATED(fire_TH_diag%avg_FRP_gridcell)) &
        fire_TH_diag%avg_FRP_gridcell      (:)   = fire_TH_diag%avg_FRP_gridcell     (:)   / OutPer
      IF (ASSOCIATED(fire_TH_diag%avg_numfires_gridcell)) &
        fire_TH_diag%avg_numfires_gridcell (:)   = fire_TH_diag%avg_numfires_gridcell(:)   / OutPer
    ENDIF
    
    cbalance_diag%avg_Cpool_green(:,:) = cbalance_diag%avg_Cpool_green(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cpool_woods(:,:) = cbalance_diag%avg_Cpool_woods(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cpool_reserve(:,:) = cbalance_diag%avg_Cpool_reserve(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cpool_litter_green_ag(:,:) = cbalance_diag%avg_Cpool_litter_green_ag(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cpool_litter_green_bg(:,:) = cbalance_diag%avg_Cpool_litter_green_bg(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cpool_litter_wood_ag(:,:)  = cbalance_diag%avg_Cpool_litter_wood_ag(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cpool_litter_wood_bg(:,:)  = cbalance_diag%avg_Cpool_litter_wood_bg(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cpool_slow(:,:) = cbalance_diag%avg_Cpool_slow(:,:)/real(noOfTimeSteps)

    cbalance_diag%avg_soil_respiration(:,:)    = cbalance_diag%avg_soil_respiration(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_soil_respiration_pot(:,:) = cbalance_diag%avg_soil_respiration_pot(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_NPP_yDayMean(:,:)        = cbalance_diag%avg_NPP_yDayMean(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_NPP_flux_correction(:,:) = cbalance_diag%avg_NPP_flux_correction(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_excess_NPP(:,:)          = cbalance_diag%avg_excess_NPP(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_root_exudates(:,:)       = cbalance_diag%avg_root_exudates(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cflux_herbivory(:,:)     = cbalance_diag%avg_Cflux_herbivory(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cflux_herbivory_LG(:,:)  = cbalance_diag%avg_Cflux_herbivory_LG(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_Cflux_herbivory_2_atm(:,:)= cbalance_diag%avg_Cflux_herbivory_2_atm(:,:)/real(noOfTimeSteps)
    cbalance_diag%avg_box_NEP(:)               = cbalance_diag%avg_box_NEP(:)/real(noOfTimeSteps)
    cbalance_diag%avg_NPP_act_yDayMean(:,:)    = cbalance_diag%avg_NPP_act_yDayMean(:,:)/real(noOfTimeSteps)

    IF (with_yasso) THEN
       cbalance_diag%avg_YCpool_acid_ag1(:,:)       = cbalance_diag%avg_YCpool_acid_ag1(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_water_ag1(:,:)      = cbalance_diag%avg_YCpool_water_ag1(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_ethanol_ag1(:,:)    = cbalance_diag%avg_YCpool_ethanol_ag1(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_nonsoluble_ag1(:,:) = cbalance_diag%avg_YCpool_nonsoluble_ag1(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_acid_bg1(:,:)       = cbalance_diag%avg_YCpool_acid_bg1(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_water_bg1(:,:)      = cbalance_diag%avg_YCpool_water_bg1(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_ethanol_bg1(:,:)    = cbalance_diag%avg_YCpool_ethanol_bg1(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_nonsoluble_bg1(:,:) = cbalance_diag%avg_YCpool_nonsoluble_bg1(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_humus_1(:,:)        = cbalance_diag%avg_YCpool_humus_1(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_acid_ag2(:,:)       = cbalance_diag%avg_YCpool_acid_ag2(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_water_ag2(:,:)      = cbalance_diag%avg_YCpool_water_ag2(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_ethanol_ag2(:,:)    = cbalance_diag%avg_YCpool_ethanol_ag2(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_nonsoluble_ag2(:,:) = cbalance_diag%avg_YCpool_nonsoluble_ag2(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_acid_bg2(:,:)       = cbalance_diag%avg_YCpool_acid_bg2(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_water_bg2(:,:)      = cbalance_diag%avg_YCpool_water_bg2(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_ethanol_bg2(:,:)    = cbalance_diag%avg_YCpool_ethanol_bg2(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_nonsoluble_bg2(:,:) = cbalance_diag%avg_YCpool_nonsoluble_bg2(:,:)/real(noOfTimeSteps)
       cbalance_diag%avg_YCpool_humus_2(:,:)         = cbalance_diag%avg_YCpool_humus_2(:,:)/real(noOfTimeSteps)

    END IF

    IF (with_nitrogen) THEN
       nbalance_offline%avg_Npool_green(:,:) = nbalance_offline%avg_Npool_green(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_Npool_woods(:,:) = nbalance_offline%avg_Npool_woods(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_Npool_mobile(:,:)= nbalance_offline%avg_Npool_mobile(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_Npool_litter_green_ag(:,:) = nbalance_offline%avg_Npool_litter_green_ag(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_Npool_litter_green_bg(:,:) = nbalance_offline%avg_Npool_litter_green_bg(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_Npool_litter_wood_ag(:,:)  = nbalance_offline%avg_Npool_litter_wood_ag(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_Npool_litter_wood_bg(:,:)  = nbalance_offline%avg_Npool_litter_wood_bg(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_Npool_slow(:,:) = nbalance_offline%avg_Npool_slow(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_SMINN_pool(:,:) = nbalance_offline%avg_SMINN_pool(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_minNflux_litter_green_ag(:,:) = nbalance_offline%avg_minNflux_litter_green_ag(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_minNflux_litter_green_bg(:,:) = nbalance_offline%avg_minNflux_litter_green_bg(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_minNflux_litter_wood_ag(:,:)  = nbalance_offline%avg_minNflux_litter_wood_ag(:,:)  /real(noOfTimeSteps)
       nbalance_offline%avg_minNflux_litter_wood_bg(:,:)  = nbalance_offline%avg_minNflux_litter_wood_bg(:,:)  /real(noOfTimeSteps)
       nbalance_offline%avg_minNflux_slow(:,:)   = nbalance_offline%avg_minNflux_slow(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_Nplant_demand(:,:)   = nbalance_offline%avg_Nplant_demand(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_Nsoil_demand(:,:)    = nbalance_offline%avg_Nsoil_demand(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_Ntotal_demand(:,:)   = nbalance_offline%avg_Ntotal_demand(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_SMINN_herbivory(:,:)= nbalance_offline%avg_SMINN_herbivory(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_nfix_to_sminn(:,:)  = nbalance_offline%avg_nfix_to_sminn(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_sminn_to_denit(:,:) = nbalance_offline%avg_sminn_to_denit(:,:)/real(noOfTimeSteps)
       nbalance_offline%avg_sminn_leach(:,:)    = nbalance_offline%avg_sminn_leach(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_N2O_emissions_ExternalN(:,:)  = nbalance_offline%avg_N2O_emissions_ExternalN(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_N2O_emissions_mineraliz(:,:)  = nbalance_offline%avg_N2O_emissions_mineraliz(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_NetEcosyst_N_flux(:,:) = nbalance_offline%avg_NetEcosyst_N_flux(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_nfert_to_sminn(:,:)    = nbalance_offline%avg_nfert_to_sminn(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_N2O_emissions_nfert(:,:) = nbalance_offline%avg_N2O_emissions_nfert(:,:) /real(noOfTimeSteps)
       nbalance_offline%avg_N2O_emissions_grazing(:,:)  = nbalance_offline%avg_N2O_emissions_grazing(:,:) /real(noOfTimeSteps)
    END IF

    IF (use_external_landcover_maps .OR. with_landuse_transitions) THEN

       landcover_change%LCC_flux_box_C2atmos(:) = landcover_change%LCC_flux_box_C2atmos(:) / real(noOfTimeSteps)

       IF (with_landuse_transitions) THEN
          landuse_transitions%Box_flux_harvest       (:) = landuse_transitions%Box_flux_harvest       (:) / OutPer
          landuse_transitions%Box_flux_harvest_2atmos(:) = landuse_transitions%Box_flux_harvest_2atmos(:) / OutPer
       END IF

       IF (lcc_scheme==1) THEN 
          landcover_change%LCC_flux_box_C2litterGreenPools(:) = &
          landcover_change%LCC_flux_box_C2litterGreenPools(:) / real(noOfTimeSteps)
          landcover_change%LCC_flux_box_C2litterWoodPool  (:) = &
          landcover_change%LCC_flux_box_C2litterWoodPool  (:) / real(noOfTimeSteps)
       ENDIF

       IF (lcc_scheme == 2) THEN

          !atmos
          landcover_change%boxC_onSite_2_atmos       (:) = landcover_change%boxC_onSite_2_atmos       (:) / OutPer
          landcover_change%boxC_paper_2_atmos        (:) = landcover_change%boxC_paper_2_atmos        (:) / OutPer
          landcover_change%boxC_construction_2_atmos (:) = landcover_change%boxC_construction_2_atmos (:) / OutPer
          
          !transfer
          landcover_change%boxC_2_onSite              (:) = landcover_change%boxC_2_onSite              (:) / OutPer
          landcover_change%boxC_2_paper               (:) = landcover_change%boxC_2_paper               (:) / OutPer
          landcover_change%boxC_2_construction        (:) = landcover_change%boxC_2_construction        (:) / OutPer

          IF (with_landuse_transitions) THEN
             landcover_change%boxC_paper_harv_2_atmos       (:) = landcover_change%boxC_paper_harv_2_atmos       (:) / OutPer
             landcover_change%boxC_construction_harv_2_atmos(:) = landcover_change%boxC_construction_harv_2_atmos(:) / OutPer
             landcover_change%boxC_2_paper_harv             (:) = landcover_change%boxC_2_paper_harv             (:) / OutPer
             landcover_change%boxC_2_construction_harv      (:) = landcover_change%boxC_2_construction_harv      (:) / OutPer
          ENDIF

          ! Pools 
          cbalance_diag%avg_Cpool_onSite      (:) = cbalance_diag%avg_Cpool_onSite      (:) / real(noOfTimeSteps)
          cbalance_diag%avg_Cpool_paper       (:) = cbalance_diag%avg_Cpool_paper       (:) / real(noOfTimeSteps)
          cbalance_diag%avg_Cpool_construction(:) = cbalance_diag%avg_Cpool_construction(:) / real(noOfTimeSteps)

          IF (with_landuse_transitions) THEN
             cbalance_diag%avg_Cpool_paper_harvest       (:) = cbalance_diag%avg_Cpool_paper_harvest       (:)/real(noOfTimeSteps)
             cbalance_diag%avg_Cpool_construction_harvest(:) = cbalance_diag%avg_Cpool_construction_harvest(:)/real(noOfTimeSteps)
          ENDIF
       ENDIF
          
    ENDIF ! landcover_maps or landuse_transitions

  END SUBROUTINE averageFields

  SUBROUTINE nullifyFields(cbalance_diag, disturbance, fire_TH_diag, landcover_change, &
       landuse_transitions, lcc_scheme)
    USE mo_cbalone_memory,        ONLY: cbal_offline_type
    USE mo_cbal_landcover_change, ONLY: landcover_change_type
    USE mo_disturbance,           ONLY: dist_opts
    USE mo_disturbance_thonicke,  ONLY: fire_TH_diag_type

    TYPE(cbal_offline_type), INTENT(inout)       :: cbalance_diag
    TYPE(disturbance_type), INTENT(inout)        :: disturbance
    TYPE(fire_TH_diag_type), intent(inout)       :: fire_TH_diag
    TYPE(landcover_change_type),intent(inout)    :: landcover_change
    TYPE(landuse_transitions_type),intent(inout) :: landuse_transitions

    INTEGER,                       INTENT(in)    :: lcc_scheme

    cbalance_diag%avg_Cpool_green(:,:)   = 0._dp
    cbalance_diag%avg_Cpool_woods(:,:)   = 0._dp
    cbalance_diag%avg_Cpool_reserve(:,:) = 0._dp
    cbalance_diag%avg_Cpool_litter_green_ag(:,:)   = 0._dp
    cbalance_diag%avg_Cpool_litter_green_bg(:,:)   = 0._dp
    cbalance_diag%avg_Cpool_litter_wood_ag(:,:)    = 0._dp
    cbalance_diag%avg_Cpool_litter_wood_bg(:,:)    = 0._dp
    cbalance_diag%avg_Cpool_slow(:,:)    = 0._dp

    cbalance_diag%avg_soil_respiration(:,:)   = 0._dp
    cbalance_diag%avg_soil_respiration_pot(:,:) = 0._dp
    cbalance_diag%avg_NPP_yDayMean(:,:)       = 0._dp
    cbalance_diag%avg_NPP_flux_correction(:,:)= 0._dp
    cbalance_diag%avg_excess_NPP(:,:)         = 0._dp
    cbalance_diag%avg_root_exudates(:,:)      = 0._dp
    cbalance_diag%avg_Cflux_herbivory(:,:)    = 0._dp
    cbalance_diag%avg_Cflux_herbivory_LG(:,:) = 0._dp
    cbalance_diag%avg_Cflux_herbivory_2_atm(:,:) = 0._dp
    cbalance_diag%avg_box_NEP(:)              = 0._dp
    cbalance_diag%avg_NPP_act_yDayMean(:,:)   = 0._dp

    IF (with_yasso) THEN
       cbalance_diag%avg_YCpool_acid_ag1      (:,:) =0._dp
       cbalance_diag%avg_YCpool_water_ag1     (:,:) =0._dp
       cbalance_diag%avg_YCpool_ethanol_ag1   (:,:) =0._dp
       cbalance_diag%avg_YCpool_nonsoluble_ag1(:,:) =0._dp
       cbalance_diag%avg_YCpool_acid_bg1      (:,:) =0._dp
       cbalance_diag%avg_YCpool_water_bg1     (:,:) =0._dp
       cbalance_diag%avg_YCpool_ethanol_bg1   (:,:) =0._dp
       cbalance_diag%avg_YCpool_nonsoluble_bg1(:,:) =0._dp
       cbalance_diag%avg_YCpool_humus_1       (:,:) =0._dp
       cbalance_diag%avg_YCpool_acid_ag2      (:,:) =0._dp
       cbalance_diag%avg_YCpool_water_ag2     (:,:) =0._dp
       cbalance_diag%avg_YCpool_ethanol_ag2   (:,:) =0._dp
       cbalance_diag%avg_YCpool_nonsoluble_ag2(:,:) =0._dp
       cbalance_diag%avg_YCpool_acid_bg2      (:,:) =0._dp
       cbalance_diag%avg_YCpool_water_bg2     (:,:) =0._dp
       cbalance_diag%avg_YCpool_ethanol_bg2   (:,:) =0._dp
       cbalance_diag%avg_YCpool_nonsoluble_bg2(:,:) =0._dp
       cbalance_diag%avg_YCpool_humus_2       (:,:) =0._dp
    END IF

    IF (with_nitrogen) THEN
       nbalance_offline%avg_Npool_green(:,:)    = 0._dp
       nbalance_offline%avg_Npool_woods(:,:)    = 0._dp
       nbalance_offline%avg_Npool_mobile(:,:)   = 0._dp
       nbalance_offline%avg_Npool_litter_green_ag(:,:) = 0._dp
       nbalance_offline%avg_Npool_litter_green_bg(:,:) = 0._dp
       nbalance_offline%avg_Npool_litter_wood_ag(:,:)  = 0._dp
       nbalance_offline%avg_Npool_litter_wood_bg(:,:)  = 0._dp
       nbalance_offline%avg_Npool_slow(:,:)     = 0._dp
       nbalance_offline%avg_SMINN_pool(:,:)     = 0._dp
       nbalance_offline%avg_minNflux_litter_green_ag(:,:) = 0._dp
       nbalance_offline%avg_minNflux_litter_green_bg(:,:) = 0._dp
       nbalance_offline%avg_minNflux_litter_wood_ag(:,:)  = 0._dp
       nbalance_offline%avg_minNflux_litter_wood_bg(:,:)  = 0._dp
       nbalance_offline%avg_minNflux_slow(:,:)  = 0._dp
       nbalance_offline%avg_Nplant_demand(:,:)  = 0._dp
       nbalance_offline%avg_Nsoil_demand(:,:)   = 0._dp
       nbalance_offline%avg_Ntotal_demand(:,:)  = 0._dp
       nbalance_offline%avg_SMINN_herbivory(:,:)= 0._dp
       nbalance_offline%avg_nfix_to_sminn(:,:)  = 0._dp
       nbalance_offline%avg_sminn_to_denit(:,:) = 0._dp
       nbalance_offline%avg_sminn_leach(:,:)    = 0._dp
       nbalance_offline%avg_N2O_emissions_ExternalN(:,:)   = 0._dp
       nbalance_offline%avg_N2O_emissions_mineraliz(:,:)   = 0._dp
       nbalance_offline%avg_NetEcosyst_N_flux(:,:) = 0._dp
       nbalance_offline%avg_nfert_to_sminn(:,:)    = 0._dp
       nbalance_offline%avg_N2O_emissions_nfert(:,:) = 0._dp
       nbalance_offline%avg_N2O_emissions_grazing(:,:) = 0._dp           
    END IF

    IF (use_external_landcover_maps .OR. with_landuse_transitions) THEN
        landcover_change%LCC_flux_box_C2atmos(:) = 0.0_dp
       IF (lcc_scheme == 1) THEN
          landcover_change%LCC_flux_box_C2litterGreenPools(:) = 0.0_dp
          landcover_change%LCC_flux_box_C2litterWoodPool(:)   = 0.0_dp
          
          IF(with_nitrogen) THEN 
             landcover_change%LCC_flux_box_N2atmos(:)            = 0.0_dp
             landcover_change%LCC_flux_box_N2litterGreenPools(:) = 0.0_dp
             landcover_change%LCC_flux_box_N2litterWoodPool(:)   = 0.0_dp
             landcover_change%LCC_flux_box_N2SMINNpool(:)        = 0.0_dp
          END IF
       ENDIF

       IF (lcc_scheme==2) THEN
          ! fluxes
          landcover_change%boxC_onSite_2_atmos              (:) = 0.0_dp
          landcover_change%boxC_paper_2_atmos               (:) = 0.0_dp
          landcover_change%boxC_construction_2_atmos        (:) = 0.0_dp
          
          !transfer
          landcover_change%boxC_2_onSite              (:) = 0.0_dp
          landcover_change%boxC_2_paper               (:) = 0.0_dp
          landcover_change%boxC_2_construction        (:) = 0.0_dp

          IF ( with_landuse_transitions ) THEN
             landcover_change%boxC_paper_harv_2_atmos       (:) = 0.0_dp
             landcover_change%boxC_construction_harv_2_atmos(:) = 0.0_dp
             landcover_change%boxC_2_paper_harv             (:) = 0.0_dp
             landcover_change%boxC_2_construction_harv      (:) = 0.0_dp
          ENDIF

          !pools 
          cbalance_diag%avg_Cpool_onSite              (:) = 0.0_dp
          cbalance_diag%avg_Cpool_paper               (:) = 0.0_dp
          cbalance_diag%avg_Cpool_construction        (:) = 0.0_dp
          IF (with_landuse_transitions) THEN
             cbalance_diag%avg_Cpool_paper_harvest       (:) = 0.0_dp
             cbalance_diag%avg_Cpool_construction_harvest(:) = 0.0_dp
          ENDIF
       ENDIF
    END IF
    
    IF (with_landuse_transitions) THEN
       landuse_transitions%Box_flux_harvest(:) = 0.0_dp
       landuse_transitions%Box_flux_harvest_2atmos(:) = 0.0_dp
    END IF

    IF (use_disturbance) THEN
      disturbance%box_burned_frac_avg     (:,:) = 0.0_dp
      disturbance%box_burned_frac_diag_avg(:,:) = 0.0_dp
      disturbance%box_damaged_frac_avg    (:,:) = 0.0_dp
      disturbance%box_CO2_flux_2_atmos    (:)   = 0.0_dp

      IF (dist_opts%fire_algorithm == 3) THEN
        IF (dist_opts%ldiag) THEN
          fire_TH_diag%avg_FDI(:) = 0.0_dp
          fire_TH_diag%avg_ROS(:) = 0.0_dp
          fire_TH_diag%avg_fuel_1hr_total(:)   = 0.0_dp
          fire_TH_diag%avg_fuel_10hr_total(:)  = 0.0_dp
          fire_TH_diag%avg_fuel_100hr_total(:) = 0.0_dp
          fire_TH_diag%avg_numfire(:)          = 0.0_dp
          fire_TH_diag%avg_HumanIgnition       = 0.0_dp
          fire_TH_diag%avg_LightningIgnition   = 0.0_dp
          fire_TH_diag%avg_FireDuration        = 0.0_dp
          fire_TH_diag%avg_fire_intensity      = 0.0_dp
          fire_TH_diag%avg_postfire_mortality  = 0.0_dp
          fire_TH_diag%avg_vegetation_height   = 0.0_dp
          fire_TH_diag%avg_carbon_2_atmos      = 0.0_dp
        ENDIF
        IF (ASSOCIATED(fire_TH_diag%avg_FRP_gridcell     )) fire_TH_diag%avg_FRP_gridcell      = 0.0_dp
        IF (ASSOCIATED(fire_TH_diag%avg_numfires_gridcell)) fire_TH_diag%avg_numfires_gridcell = 0.0_dp
      ENDIF
    ENDIF

  END SUBROUTINE nullifyFields

  !! --- getYearDay() -----------------------------
  !!
  !! Returns the number of a day in a year, when counts starts with 1 at first of january
  !!
  integer function getYearDay(year,month,day)
    integer,intent(in) :: year
    integer,intent(in) :: month
    integer,intent(in) :: day

    logical :: leapYear

    integer,parameter :: lastDayOfPrevMonth_normal(1:12)   =(/0,31,59,90,120,151,181,212,243,273,304,334/)
    integer,parameter :: lastDayOfPrevMonth_leapYear(1:12) =(/0,31,60,91,121,152,182,213,244,274,305,335/)

    leapYear = (MOD(year,4)==0 .AND. MOD(year,100)/=0) .OR. MOD(year,400)==0 
    if(leapYear) then
       getYearDay = lastDayOfPrevMonth_leapYear(month) + day
    else
       getYearDay = lastDayOfPrevMonth_normal(month) + day
    end if

  end function getYearDay

  !! --- get_GregorianYearLength() -----------------------------
  !!
  !! Returns the number days in a year according to the (proleptic) gregorian calender 
  !!
  integer function get_GregorianYearLength(year)
    integer,intent(in) :: year

    if ( (MOD(year,4)==0 .AND. MOD(year,100)/=0) .OR. MOD(year,400)==0 ) then
       get_gregorianYearLength = 366
    else
       get_gregorianYearLength = 365
    end if

  end function get_gregorianYearLength

end program cbalone_driver
