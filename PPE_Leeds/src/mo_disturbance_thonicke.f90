!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!! The algorithms implemented in this file are used and distributed with permission from 
!! Kirsten Thonicke.
!!
MODULE mo_disturbance_thonicke
!
! Calculation of daily fire disturbance - Thonicke et al (2010)(SPITFIRE) (TH2010 in this file)
!
!
! History:
! 2011-10 GL,  MPI-M: started implementing spitfire into JSBACH
! 2014-02 StW, MPI-M: Included in cosmos-landveg-fire
USE mo_kind,          ONLY: dp
USE mo_jsbach_lctlib, ONLY: lctlib_type
USE mo_exception,     ONLY: finish, message, message_text
USE mo_jsbach_grid,   ONLY: grid_type, domain_type
USE mo_land_surface,  ONLY: land_surface_type
USE mo_input,         ONLY: InputFileAdd, InputVarAdd, InputVarGroupStart, InputVarGroupEnd, InputVarIsUpdated, &
                            INPUT_FILE_REPEAT_ALL, INPUT_FILE_REPEAT_NONE, INPUT_VAR_VALID_ACTUAL, INPUT_VAR_VALID_MIDPOINT, &
                            InputGetData, input_file_list, input_var_list
IMPLICIT NONE

! Default namelist parameter values
REAL(dp), PARAMETER :: def_SurfArea2Vol(3)          = (/ 66.0_dp,3.58_dp,0.98_dp /)
REAL(dp), PARAMETER :: def_a_nd              =     0.251_dp
REAL(dp), PARAMETER :: def_human_para        =     0.65_dp
REAL(dp), PARAMETER :: def_ign_para          =     1.0_dp
REAL(dp), PARAMETER :: def_wind_limit        =     0.01_dp
REAL(dp), PARAMETER :: def_wind_max          =   100._dp
REAL(dp), PARAMETER :: def_wind_slope        =     0.5_dp
REAL(dp), PARAMETER :: def_moisture_scaling  = 13000.0_dp
REAL(dp), PARAMETER :: def_Lethal2ResTime    = 5.0_dp
LOGICAL,  PARAMETER :: def_llight_ground     = .FALSE.
LOGICAL,  PARAMETER :: def_lwind_speed_limit = .TRUE.
LOGICAL,  PARAMETER :: def_lcalc_frp         = .FALSE.
LOGICAL,  PARAMETER :: def_lduration_popd    = .FALSE.

! Other parameters
REAL(dp), PARAMETER :: MINER_TOT             =     0.055_dp
REAL(dp), PARAMETER :: MINER_DAMP            =     0.41739_dp
REAL(dp), PARAMETER :: fbd10hr_2_1hr         =     1.2_dp
REAL(dp), PARAMETER :: fbd100hr_2_1hr        =     1.4_dp     ! from fuel bulk densities to 1hr fuel type equivs
REAL(dp), PARAMETER :: fbd_C3_livegrass      =     4.0_dp     ! fuel density of livegrass [kg/m3]
REAL(dp), PARAMETER :: fbd_C4_livegrass      =     4.0_dp     ! fuel density of livegrass [kg/m3]
REAL(dp), PARAMETER :: SurfArea2Vol_1hr      =    66.0_dp     ! surface area to volume ratio
REAL(dp), PARAMETER :: SurfArea2Vol_10hr     =     3.58_dp    ! surface area to volume ratio
REAL(dp), PARAMETER :: SurfArea2Vol_100hr    =     0.98_dp    ! surface area to volume ratio 
REAL(dp), PARAMETER :: moistfactor_livegrass =     0.2_dp     
REAL(dp), PARAMETER :: frac_1hr_wood_new     =     0.0_dp     ! frac of    1hr fuel in new carbon in woody litter pools
REAL(dp), PARAMETER :: frac_10hr_wood_new    =     0.12_dp    ! frac of   10hr fuel in new carbon in woody litter pools
REAL(dp), PARAMETER :: frac_100hr_wood_new   =     0.21_dp    ! frac of  100hr fuel in new carbon in woody litter pools
REAL(dp), PARAMETER :: frac_1000hr_wood_new  =     0.67_dp    ! frac of 1000hr fuel in new carbon in woody litter pools
REAL(dp), PARAMETER :: frac_green_active     =     0.5_dp     ! fraction of woody pft green pool which is active with fire
REAL(dp), PARAMETER :: HeatContentFuel       = 18000.0_dp     ! heat content of fuel in [kJ kg-1]
REAL(dp), PARAMETER :: hmax                  =    24.19_dp    ! maximum vegetation height
REAL(dp), PARAMETER :: hl                    =     0.19_dp    ! parameter to compute vegetation height
REAL(dp), PARAMETER :: ft_per_m              =     3.281_dp   ! Conversion from m -> feet
REAL(dp), PARAMETER :: part_dens             =   513.0_dp
REAL(dp), PARAMETER :: alpha_1hr             =     0.001_dp
REAL(dp), PARAMETER :: alpha_10hr            =     0.00005424_dp
REAL(dp), PARAMETER :: alpha_100hr           =     0.00001485_dp
REAL(dp), PARAMETER :: alpha_1000hr          =     0.000001_dp
REAL(dp), PARAMETER :: allom2                =    40._dp
REAL(dp), PARAMETER :: allom3                =     0.67_dp

! fire thonicke option type
TYPE fire_TH_options_type
  REAL(dp)          :: a_nd
  REAL(dp)          :: human_para
  REAL(dp)          :: ign_para
  REAL(dp)          :: wind_limit
  REAL(dp)          :: SurfArea2Vol(3)
  REAL(dp)          :: moisture_scaling
  REAL(dp)          :: wind_max
  REAL(dp)          :: wind_slope
  REAL(dp)          :: Lethal2ResTime
  LOGICAL           :: llight_ground
  LOGICAL           :: lwind_speed_limit
  LOGICAL           :: read_fuel_frac
  LOGICAL           :: lduration_popd
  LOGICAL           :: lcalc_frp
  LOGICAL           :: ldiag
END TYPE fire_TH_options_type

! Type holding module state variables 
TYPE fire_TH_state_type
  ! used in fuel_consumption and relocate_carbon_fire_thonicke  
  REAL(dp), POINTER :: NI_acc(:)
  REAL(dp), POINTER :: postfire_mortality(:,:)
  REAL(dp), POINTER :: mortality_crownkill(:,:)
  REAL(dp), POINTER :: burned_Cpool_litter_wood_ag(:,:)
  REAL(dp), POINTER :: burned_Cpool_litter_green_ag(:,:)
  REAL(dp), POINTER :: burned_Cpool_wood(:,:)
  REAL(dp), POINTER :: burned_Cpool_green(:,:)
  REAL(dp), POINTER :: burned_Cpool_reserve(:,:)
  ! Pools for yasso
  REAL(dp), POINTER :: burned_YC_acid_ag1      (:,:)
  REAL(dp), POINTER :: burned_YC_acid_ag2      (:,:)
  REAL(dp), POINTER :: burned_YC_water_ag1     (:,:)
  REAL(dp), POINTER :: burned_YC_water_ag2     (:,:)
  REAL(dp), POINTER :: burned_YC_ethanol_ag1   (:,:)
  REAL(dp), POINTER :: burned_YC_ethanol_ag2   (:,:)
  REAL(dp), POINTER :: burned_YC_nonsoluble_ag1(:,:)
  REAL(dp), POINTER :: burned_YC_nonsoluble_ag2(:,:)
  ! fractions of different fuel types needed for spitfire
  REAL(dp), POINTER :: frac_1hr_wood(:,:)
  REAL(dp), POINTER :: frac_10hr_wood(:,:)
  REAL(dp), POINTER :: frac_100hr_wood(:,:)
  REAL(dp), POINTER :: frac_1000hr_wood(:,:)
END TYPE fire_TH_state_type

! Type for diagnostic module variables
TYPE fire_TH_diag_type
  REAL(dp), POINTER :: avg_FDI(:)
  REAL(dp), POINTER :: avg_ROS(:)
  REAL(dp), POINTER :: avg_fuel_1hr_total(:)
  REAL(dp), POINTER :: avg_fuel_10hr_total(:)
  REAL(dp), POINTER :: avg_fuel_100hr_total(:)
  REAL(dp), POINTER :: avg_numfire(:)
  REAL(dp), POINTER :: avg_HumanIgnition(:)
  REAL(dp), POINTER :: avg_LightningIgnition(:)
  REAL(dp), POINTER :: avg_FireDuration(:)
  REAL(dp), POINTER :: avg_fire_intensity(:)
  REAL(dp), POINTER :: avg_postfire_mortality(:,:)
  REAL(dp), POINTER :: avg_vegetation_height(:,:)
  REAL(dp), POINTER :: avg_carbon_2_atmos(:,:)
  REAL(dp), POINTER :: IR_ROS(:)
  REAL(dp), POINTER :: pop(:)
  REAL(dp), POINTER :: lightning(:)
  REAL(dp), POINTER :: lat(:)
  REAL(dp), POINTER :: a_nd_array(:)
  REAL(dp), POINTER :: cellarea(:)
  REAL(dp), POINTER :: avg_FRP_gridcell(:)
  REAL(dp), POINTER :: avg_numfires_gridcell(:)
END TYPE fire_TH_diag_type

! Module variables
TYPE (fire_TH_diag_type),       PUBLIC  :: fire_TH_diag
TYPE (fire_TH_options_type),    PUBLIC  :: fire_TH_options 
TYPE (fire_TH_state_type),      PUBLIC  :: state_FTH
TYPE (fire_TH_state_type),      PRIVATE :: init_state_FTH
LOGICAL,                        PRIVATE :: module_initialized = .FALSE. ! Have config and init  been called?
TYPE (input_var_list), POINTER, PRIVATE :: lightning_var

!function declarations
PUBLIC  config_fire_thonicke
PUBLIC  burned_frac_thonicke
PUBLIC  init_fire_thonicke
PUBLIC  init_fire_thonicke_memory
PUBLIC  relocate_carbon_fire_thonicke
PUBLIC  FRP_nfire_per_gridcell

CONTAINS

  SUBROUTINE config_fire_thonicke(nml_unit,ldiag)
  USE mo_namelist,         ONLY: position_nml, POSITIONED, LENGTH_ERROR, READ_ERROR
  USE mo_io_units,         ONLY: nout
  USE mo_mpi,              ONLY: p_parallel_io, p_io, p_parallel, p_bcast

! Input parameters
  INTEGER, INTENT(in) :: nml_unit
  LOGICAL, INTENT(in) :: ldiag

! Local namelist parameters
  REAL(dp) :: a_nd
  REAL(dp) :: human_para
  REAL(dp) :: ign_para
  REAL(dp) :: wind_limit
  REAL(dp) :: wind_max
  REAL(dp) :: wind_slope
  REAL(dp) :: Lethal2ResTime
  REAL(dp) :: SurfArea2Vol(3)
  REAL(dp) :: moisture_scaling
  LOGICAL  :: llight_ground
  LOGICAL  :: lwind_speed_limit
  LOGICAL  :: read_fuel_frac
  LOGICAL  :: lduration_popd
  LOGICAL  :: lcalc_frp

!other locals
  INTEGER  :: read_status
  INTEGER  :: f_unit
  REAL(dp) :: rbuf(11)
  LOGICAL  :: lbuf(5)

  INCLUDE 'fire_thonicke_ctl.inc'

    IF (module_initialized) CALL finish('config_fire_thonicke','Attempted to initialize module twice')
    IF (p_parallel_io) THEN
      ! Map default values
      a_nd               = def_a_nd
      human_para         = def_human_para
      ign_para           = def_ign_para
      wind_limit         = def_wind_limit
      SurfArea2Vol       = def_SurfArea2Vol
      moisture_scaling   = def_moisture_scaling
      wind_max           = def_wind_max
      wind_slope         = def_wind_slope
      Lethal2ResTime     = def_Lethal2ResTime
      llight_ground      = def_llight_ground
      lwind_speed_limit  = def_lwind_speed_limit
      lcalc_frp          = def_lcalc_frp
      lduration_popd     = def_lduration_popd
      read_fuel_frac     = .FALSE.
      f_unit = position_nml ('FIRE_THONICKE_CTL', nml_unit, status=read_status)
      SELECT CASE (read_status)
        CASE (POSITIONED)
          READ (f_unit, fire_thonicke_ctl)
          CALL message('config_fire_thonicke','Namelist FIRE_THONICKE_CTL: ')
          WRITE(nout, fire_thonicke_ctl)
        CASE (LENGTH_ERROR)
          CALL finish ('config_fire_thonicke','Length error in namelist fire_thonicke_ctl')
        CASE (READ_ERROR)
          CALL finish ('config_fire_thonicke','Error reading namelist fire_thonicke_ctl ')
      END SELECT
    ENDIF

    ! Distribute namelist parameters
    rbuf = (/a_nd,human_para,ign_para,wind_limit,SurfArea2Vol(1),SurfArea2Vol(2),& 
             SurfArea2Vol(3),moisture_scaling, wind_max, wind_slope, Lethal2ResTime/)
    lbuf = (/llight_ground,lwind_speed_limit,read_fuel_frac,lduration_popd,lcalc_frp/)

    IF (p_parallel) THEN
      CALL p_bcast(rbuf,p_io)
      CALL p_bcast(lbuf,p_io)
    ENDIF
 
    fire_TH_options%llight_ground     = lbuf(1)
    fire_TH_options%lwind_speed_limit = lbuf(2)
    fire_TH_options%read_fuel_frac    = lbuf(3)
    fire_TH_options%lduration_popd    = lbuf(4)
    fire_TH_options%lcalc_frp         = lbuf(5)
    fire_TH_options%a_nd              = rbuf(1)
    fire_TH_options%human_para        = rbuf(2)
    fire_TH_options%ign_para          = rbuf(3)
    fire_TH_options%wind_limit        = rbuf(4)
    fire_TH_options%SurfArea2Vol      = rbuf(5:7)
    fire_TH_options%moisture_scaling  = rbuf(8)
    fire_TH_options%wind_max          = rbuf(9)
    fire_TH_options%wind_slope        = rbuf(10)
    fire_TH_options%Lethal2ResTime    = rbuf(11)
    fire_TH_options%ldiag             = ldiag

  END SUBROUTINE config_fire_thonicke

  SUBROUTINE init_fire_thonicke(grid,domain,surface,lcbalone,ldiag,with_yasso,nml_unt)
  ! input parameters
  TYPE(grid_type),             INTENT(in) :: grid
  TYPE(domain_type),           INTENT(in) :: domain
  TYPE(land_surface_type),     INTENT(in) :: surface
  LOGICAL,                     INTENT(in) :: lcbalone
  LOGICAL,                     INTENT(in) :: ldiag       ! Produce diagnostic output?
  LOGICAL,                     INTENT(in) :: with_yasso
  INTEGER,                     INTENT(in) :: nml_unt

  ! parameters
  CHARACTER (len=10), PARAMETER :: dim1n = 'landpoints'
  CHARACTER (len=17), PARAMETER :: dim2n = 'landpoints,ntiles'

  ! local variables
  TYPE(input_file_list), POINTER :: IFile
  INTEGER                        :: st

    CALL config_fire_thonicke(nml_unt,ldiag)

    NULLIFY(IFile,fire_TH_diag%pop,fire_TH_diag%lightning,fire_TH_diag%lat, &
            fire_TH_diag%a_nd_array, fire_TH_diag%cellarea)
    IFile => InputFileAdd('population_density.nc',cycle_act=INPUT_FILE_REPEAT_NONE)
    CALL InputVarAdd('landpoints',fire_TH_diag%pop,'population_density',IFile,dt_update=-1, & 
                     valid_time=INPUT_VAR_VALID_ACTUAL)
    fire_TH_diag%lat => domain%lat
    IFile => InputFileAdd('lightning.nc',cycle_act=INPUT_FILE_REPEAT_ALL)
    CALL InputVarAdd('landpoints',fire_TH_diag%lightning,'lightning_frq',IFile,dt_update=86400, &
                     valid_time=INPUT_VAR_VALID_MIDPOINT,var_ref=lightning_var) ! Update every day
    
    IFile => InputFileAdd('a_nd_file.nc',dt_file=0)
    CALL InputVarAdd('landpoints',fire_TH_diag%a_nd_array,'a_nd',IFile,def_val=fire_TH_options%a_nd,dt_update=0)

    IF (fire_TH_options%lcalc_frp) THEN
      IFILE => InputFileAdd('cell_area.nc',dt_file=0)
      CALL InputVarAdd('landpoints',fire_TH_diag%cellarea,'cell_area',IFile,dt_update=0)
    ENDIF

    ! Read all initial data
    CALL InputGetData()

    ! Allocate memory for cbalance (in JSBACH the stream functions in init_fire_arora_memory takes care of the allocation 
    IF (lcbalone) THEN
      ! Allocate of state variables
      ALLOCATE(state_FTH%NI_acc                      (grid%nland),                &
               state_FTH%burned_Cpool_litter_wood_ag (grid%nland,surface%ntiles), &
               state_FTH%burned_Cpool_litter_green_ag(grid%nland,surface%ntiles), &
               state_FTH%burned_Cpool_wood           (grid%nland,surface%ntiles), &
               state_FTH%burned_Cpool_green          (grid%nland,surface%ntiles), &
               state_FTH%burned_Cpool_reserve        (grid%nland,surface%ntiles), &
               state_FTH%postfire_mortality          (grid%nland,surface%ntiles), &
               state_FTH%mortality_crownkill         (grid%nland,surface%ntiles), &
               state_FTH%frac_1hr_wood               (grid%nland,surface%ntiles), &
               state_FTH%frac_10hr_wood              (grid%nland,surface%ntiles), &
               state_FTH%frac_100hr_wood             (grid%nland,surface%ntiles), &
               state_FTH%frac_1000hr_wood            (grid%nland,surface%ntiles), &
               stat=st)
      IF (st /= 0) CALL finish('init_fire_thonicke','Allocation of module variables ')
      IF (with_yasso) THEN
        ALLOCATE(state_FTH%burned_YC_acid_ag1      (grid%nland,surface%ntiles), &
                 state_FTH%burned_YC_acid_ag2      (grid%nland,surface%ntiles), &
                 state_FTH%burned_YC_water_ag1     (grid%nland,surface%ntiles), &
                 state_FTH%burned_YC_water_ag2     (grid%nland,surface%ntiles), &
                 state_FTH%burned_YC_ethanol_ag1   (grid%nland,surface%ntiles), &
                 state_FTH%burned_YC_ethanol_ag2   (grid%nland,surface%ntiles), &
                 state_FTH%burned_YC_nonsoluble_ag1(grid%nland,surface%ntiles), &
                 state_FTH%burned_YC_nonsoluble_ag2(grid%nland,surface%ntiles),stat = st)
        IF (st /= 0) CALL finish('init_fire_thonicke','Allocation of module variables ')
      ENDIF

      state_FTH%NI_acc                      (:)   = 0._dp
      state_FTH%burned_Cpool_wood           (:,:) = 0._dp
      state_FTH%burned_Cpool_green          (:,:) = 0._dp
      state_FTH%burned_Cpool_reserve        (:,:) = 0._dp
      state_FTH%postfire_mortality          (:,:) = 1._dp
      state_FTH%mortality_crownkill         (:,:) = 0._dp
      state_FTH%frac_1hr_wood               (:,:) = frac_1hr_wood_new
      state_FTH%frac_10hr_wood              (:,:) = frac_10hr_wood_new
      state_FTH%frac_100hr_wood             (:,:) = frac_100hr_wood_new
      state_FTH%frac_1000hr_wood            (:,:) = frac_1000hr_wood_new
      state_FTH%burned_Cpool_litter_wood_ag (:,:) = 0._dp
      state_FTH%burned_Cpool_litter_green_ag(:,:) = 0._dp
      IF (with_yasso) THEN
        state_FTH%burned_YC_acid_ag1       = 0._dp 
        state_FTH%burned_YC_acid_ag2       = 0._dp 
        state_FTH%burned_YC_water_ag1      = 0._dp 
        state_FTH%burned_YC_water_ag2      = 0._dp 
        state_FTH%burned_YC_ethanol_ag1    = 0._dp 
        state_FTH%burned_YC_ethanol_ag2    = 0._dp 
        state_FTH%burned_YC_nonsoluble_ag1 = 0._dp
        state_FTH%burned_YC_nonsoluble_ag2 = 0._dp
      ENDIF

      ALLOCATE(fire_TH_diag%IR_ROS(grid%nland),stat=st)                                                       
      IF (st /= 0) CALL finish('init_fire_thonicke','Allocation of fire model variables ')
      fire_TH_diag%IR_ROS                 (:)   = 0._dp
      ! Allocation of diagnostic variables if requested
      IF (ldiag) THEN
        ALLOCATE(fire_TH_diag%avg_FDI                (grid%nland),                &
                 fire_TH_diag%avg_ROS                (grid%nland),                &
                 fire_TH_diag%avg_fuel_1hr_total     (grid%nland),                &
                 fire_TH_diag%avg_fuel_10hr_total    (grid%nland),                &
                 fire_TH_diag%avg_fuel_100hr_total   (grid%nland),                &
                 fire_TH_diag%avg_numfire            (grid%nland),                &
                 fire_TH_diag%avg_HumanIgnition      (grid%nland),                &  
                 fire_TH_diag%avg_LightningIgnition  (grid%nland),                &
                 fire_TH_diag%avg_FireDuration       (grid%nland),                &
                 fire_TH_diag%avg_fire_intensity     (grid%nland),                & 
                 fire_TH_diag%avg_postfire_mortality (grid%nland,surface%ntiles), &
                 fire_TH_diag%avg_vegetation_height  (grid%nland,surface%ntiles), &
                 fire_TH_diag%avg_carbon_2_atmos     (grid%nland,surface%ntiles), &
                 stat=st)                                                       
        IF (st /= 0) CALL finish('init_fire_thonicke','Allocation of fire model variables ')
        fire_TH_diag%avg_FDI                (:)   = 0._dp
        fire_TH_diag%avg_ROS                (:)   = 0._dp
        fire_TH_diag%avg_fuel_1hr_total     (:)   = 0._dp
        fire_TH_diag%avg_fuel_10hr_total    (:)   = 0._dp
        fire_TH_diag%avg_fuel_100hr_total   (:)   = 0._dp
        fire_TH_diag%avg_numfire            (:)   = 0._dp
        fire_TH_diag%avg_HumanIgnition      (:)   = 0._dp
        fire_TH_diag%avg_LightningIgnition  (:)   = 0._dp
        fire_TH_diag%avg_FireDuration       (:)   = 0._dp
        fire_TH_diag%avg_fire_intensity     (:)   = 0._dp
        fire_TH_diag%avg_postfire_mortality (:,:) = 0._dp
        fire_TH_diag%avg_vegetation_height  (:,:) = 0._dp
        fire_TH_diag%avg_carbon_2_atmos     (:,:) = 0._dp
      ENDIF
    ENDIF 
    IF (fire_TH_options%lcalc_frp) THEN
      ALLOCATE(fire_TH_diag%avg_FRP_gridcell     (grid%nland), &
               fire_TH_diag%avg_numfires_gridcell(grid%nland), &
               stat=st)
      fire_TH_diag%avg_FRP_gridcell     (:)   = 0._dp
      fire_TH_diag%avg_numfires_gridcell(:)   = 0._dp
    ELSE
      NULLIFY(fire_TH_diag%avg_FRP_gridcell,fire_TH_diag%avg_numfires_gridcell)
    ENDIF

    ! Get initial values from fuel_frac_file if requested
    NULLIFY(init_state_FTH%frac_1hr_wood,init_state_FTH%frac_10hr_wood, &
            init_state_FTH%frac_100hr_wood,init_state_FTH%frac_1000hr_wood,init_state_FTH%NI_acc) 
    IF (fire_TH_options%read_fuel_frac) THEN
      IFile => InputFileAdd('Cpools.nc', dt_file=0)
      CALL InputVarGroupStart('fuel_fractions',dt_update=0,IFile=IFile)
      CALL InputVarAdd(dim2n,init_state_FTH%frac_1hr_wood,   'frac_1hr_wood',   def_val=frac_1hr_wood_new   )
      CALL InputVarAdd(dim2n,init_state_FTH%frac_10hr_wood,  'frac_10hr_wood',  def_val=frac_10hr_wood_new  )
      CALL InputVarAdd(dim2n,init_state_FTH%frac_100hr_wood, 'frac_100hr_wood', def_val=frac_100hr_wood_new )
      CALL InputVarAdd(dim2n,init_state_FTH%frac_1000hr_wood,'frac_1000hr_wood',def_val=frac_1000hr_wood_new)
      CALL InputVarAdd(dim1n,init_state_FTH%NI_acc,          'NI_acc'          ,def_val=0._dp)
      CALL InputVarGroupEnd
      CALL InputGetData()
    ENDIF

    module_initialized = .TRUE.

  END SUBROUTINE init_fire_thonicke

  SUBROUTINE init_fire_thonicke_memory(IO_disturbance,IO_veg,dim1p,dim1,dim1n,dim2p,dim2,dim2n,with_yasso)
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: add=>add_stream_element
  USE mo_netcdf,      ONLY: max_dim_name
  USE mo_jsbach,      ONLY: mv=>missing_value
  TYPE(t_stream),             POINTER :: IO_disturbance, IO_veg
  INTEGER,                 INTENT(IN) :: dim1p(:), dim1(:),dim2p(:),dim2(:)
  CHARACTER(max_dim_name), INTENT(IN) :: dim1n(:), dim2n(:)
  LOGICAL,                 INTENT(IN) :: with_yasso

    ! --- Define variables as stream elements
    ! code numbers of this routine are 60-74 using the GRIB disturb_table
    ! code numbers of this routine are 43, 46-49 and 64-73 using the GRIB veg_table
    IF (fire_TH_options%ldiag) THEN
      CALL add(IO_disturbance,'avg_FDI',fire_TH_diag%avg_FDI,longname='average fire danger index',lpost=.TRUE.,laccu=.TRUE.,      &
               units='[]',ldims=dim1p,gdims=dim1,dimnames=dim1n,code=60,lmiss=.TRUE.,missval=mv,contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_ROS',fire_TH_diag%avg_ROS,longname='average rate of spread',lmiss=.TRUE.,missval=mv,           &
               units='m min-1',ldims=dim1p,gdims=dim1,dimnames=dim1n,code=61,lpost=.TRUE.,laccu=.TRUE.,contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_fuel_1hr_total',fire_TH_diag%avg_fuel_1hr_total,longname='average total 1hr fuel',code=62,     &
               units='g(biomass) m-2',ldims=dim1p,gdims=dim1,dimnames=dim1n,lpost=.TRUE.,laccu=.TRUE.,lmiss=.TRUE.,missval=mv,    &
               contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_fuel_10hr_total',fire_TH_diag%avg_fuel_10hr_total,longname='average total 10hr fuel',code=63,  &
               units='g(biomass) m-2',ldims=dim1p,gdims=dim1,dimnames=dim1n,lpost=.TRUE.,laccu=.TRUE.,lmiss=.TRUE.,missval=mv,    &
               contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_fuel_100hr_total',fire_TH_diag%avg_fuel_100hr_total,longname='average total 100hr fuel',       &
               units='g(biomass) m-2',ldims=dim1p,gdims=dim1,dimnames=dim1n,code=64,lpost=.TRUE.,laccu=.TRUE.,                    &
               lmiss=.TRUE.,missval=mv,contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_numfire',fire_TH_diag%avg_numfire,longname='average number of fires',lpost=.TRUE.,laccu=.TRUE.,&
               units='ha-1 day-1',ldims=dim1p, gdims=dim1,dimnames=dim1n,code=65,lmiss=.TRUE.,missval=mv,contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_HumanIgnition',fire_TH_diag%avg_HumanIgnition,longname='average number of human caused fires', &
               units='km-2 day-1',ldims=dim1p,gdims=dim1,dimnames=dim1n,code=66,lpost=.TRUE.,laccu=.TRUE.,lmiss=.TRUE.,missval=mv &
               ,contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_LightningIgnition',fire_TH_diag%avg_LightningIgnition,lmiss=.TRUE., missval=mv,                &
               longname='average number of lightning caused fires',units='km-2 day-1',ldims=dim1p,gdims=dim1,dimnames=dim1n,      &
               code=67,lpost=.TRUE.,laccu=.TRUE.,contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_FireDuration',fire_TH_diag%avg_FireDuration,longname='average fire duration',code=68,          &
               units='min',ldims=dim1p,gdims=dim1,dimnames=dim1n,lpost=.TRUE.,laccu=.TRUE.,lmiss=.TRUE.,missval=mv,               &
               contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_fire_intensity',fire_TH_diag%avg_fire_intensity,longname='average fire intensity',             &
               units='kW m-1',ldims=dim1p,gdims=dim1,dimnames=dim1n,code=69,lpost=.TRUE.,laccu=.TRUE.,lmiss=.TRUE.,missval=mv,    &
               contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_postfire_mortality',fire_TH_diag%avg_postfire_mortality,lmiss=.TRUE.,missval=mv,               &
               longname='average probability of postfire mortality of trees',units='[]',ldims=dim2p,gdims=dim2,dimnames=dim2n,    &
               code=70,lpost=.TRUE.,laccu=.TRUE.,contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_vegetation_height',fire_TH_diag%avg_vegetation_height,longname='average vegetation height',    &
               units='[m]',ldims=dim2p,gdims=dim2,dimnames=dim2n,code=71,lpost=.TRUE.,laccu=.TRUE.,lmiss=.TRUE.,missval=mv,       &
               contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_carbon_2_atmos',fire_TH_diag%avg_carbon_2_atmos,&
               longname='carbon emitted by fires',&
              units='mol Carbon day-1 m-2 vegetated area',ldims=dim2p,gdims=dim2,dimnames=dim2n,code=74,&
              lpost=.TRUE.,laccu=.TRUE.,lmiss=.TRUE.,missval=mv ,    &
            contnorest=.TRUE.)
    ENDIF
    IF (fire_TH_options%lcalc_frp) THEN
      CALL add(IO_disturbance,'avg_FRP_gridcell',fire_TH_diag%avg_FRP_gridcell,longname='average total fire radiative power',     &
               units='W',ldims=dim1p,gdims=dim1,dimnames=dim1n,lpost=.TRUE.,laccu=.TRUE.,lmiss=.TRUE.,missval=mv,code=72,         &
               contnorest=.TRUE.)
      CALL add(IO_disturbance,'avg_numfires_gridcell',fire_TH_diag%avg_numfires_gridcell,                                         &
               longname='average number of fires per gridcell',                                                                   &
               units='',ldims=dim1p,gdims=dim1,dimnames=dim1n,lpost=.TRUE.,laccu=.TRUE.,lmiss=.TRUE.,missval=mv,code=73,          &
               contnorest=.TRUE.)
    ENDIF

    CALL add(IO_veg,'NI_acc',state_FTH%NI_acc, longname='average Nesterov index', lrerun=.TRUE.,lmiss=.TRUE.,missval=mv,          &
             units='[]',ldims=dim1p, gdims=dim1, dimnames=dim1n, code=69, lpost=.TRUE.,laccu=.FALSE., contnorest=.true.)
    IF (with_yasso) THEN
      CALL add(IO_veg,'burned_YC_acid_ag_green',state_FTH%burned_YC_acid_ag1,units='[mol(C) m-2(canopy)]',                        &
               longname='Yasso C-pool of above ground burned green acid-soluble litter',lmiss=.TRUE.,missval=mv,code=64,          &
               lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p, gdims=dim2, dimnames=dim2n, contnorest=.TRUE.)
      CALL add(IO_veg,'burned_YC_acid_ag_wood' ,state_FTH%burned_YC_acid_ag2,units='[mol(C) m-2(canopy)]',                        &
               longname='Yasso C-pool of above ground burned woody acid-soluble litter',lmiss=.TRUE.,missval=mv,code=43,          &
               lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p, gdims=dim2, dimnames=dim2n, contnorest=.TRUE.)
      CALL add(IO_veg,'burned_YC_water_ag_green',state_FTH%burned_YC_water_ag1,units='[mol(C) m-2(canopy)]',                      &
               longname='Yasso C-pool of above ground burned green water-soluble litter',lmiss=.TRUE.,missval=mv,code=46,         &
               lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p, gdims=dim2, dimnames=dim2n, contnorest=.TRUE.)
      CALL add(IO_veg,'burned_YC_water_ag_wood' ,state_FTH%burned_YC_water_ag2,units='[mol(C) m-2(canopy)]',                      &
               longname='Yasso C-pool of above ground burned woody water-soluble litter',lmiss=.TRUE.,missval=mv,code=47,         &
               lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p, gdims=dim2, dimnames=dim2n, contnorest=.TRUE.)
      CALL add(IO_veg,'burned_YC_ethanol_ag_green',state_FTH%burned_YC_ethanol_ag1,units='[mol(C) m-2(canopy)]',                  &
               longname='Yasso C-pool of above ground burned green ethanol-soluble litter',lmiss=.TRUE.,missval=mv,code=48,       &
               lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p, gdims=dim2, dimnames=dim2n, contnorest=.TRUE.)
      CALL add(IO_veg,'burned_YC_ethanol_ag_wood' ,state_FTH%burned_YC_ethanol_ag2,units='[mol(C) m-2(canopy)]',                  &
               longname='Yasso C-pool of above ground burned woody ethanol-soluble litter',lmiss=.TRUE.,missval=mv,code=49,       &
               lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p, gdims=dim2, dimnames=dim2n, contnorest=.TRUE.)
      CALL add(IO_veg,'burned_YC_nonsoluble_ag_green',state_FTH%burned_YC_nonsoluble_ag1,units='[mol(C) m-2(canopy)]',            &
               longname='Yasso C-pool of above ground burned green nonsoluble litter',lmiss=.TRUE.,missval=mv,code=70,            &
               lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p, gdims=dim2, dimnames=dim2n, contnorest=.TRUE.)
      CALL add(IO_veg,'burned_YC_nonsoluble_ag_wood' ,state_FTH%burned_YC_nonsoluble_ag2,units='[mol(C) m-2(canopy)]',            &
               longname='Yasso C-pool of above ground burned woody nonsoluble litter',lmiss=.TRUE.,missval=mv,code=71,            &
               lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p, gdims=dim2, dimnames=dim2n, contnorest=.TRUE.)
    ENDIF
    CALL add(IO_veg,'burned_Cpool_litter_wood_ag',state_FTH%burned_Cpool_litter_wood_ag,                                        &
             longname='Cpool of above ground burned woody litter',lmiss=.TRUE.,missval=mv,units='[mol(C) m-2(canopy)]',code=72, &
             lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p, gdims=dim2, dimnames=dim2n, contnorest=.TRUE.)
    CALL add(IO_veg,'burned_Cpool_litter_green_ag',state_FTH%burned_Cpool_litter_green_ag,                                      &
             longname='Cpool of above ground burned green litter', lmiss=.TRUE.,missval=mv,units='[mol(C) m-2(canopy)]',code=73,&
             lpost=.TRUE.,laccu=.FALSE.,ldims=dim2p,gdims=dim2,dimnames=dim2n, contnorest=.TRUE.)
    CALL add(IO_veg,'burned_Cpool_green',state_FTH%burned_Cpool_green,longname='burned green Cpool',missval=mv,                   &
             units='[g(?) m-2(canopy)]',code=76,ldims=dim2p, gdims=dim2, dimnames=dim2n, lpost=.FALSE.,lmiss=.TRUE.,              &
             laccu=.FALSE., contnorest=.TRUE.)
    CALL add(IO_veg,'postfire_mortality',state_FTH%postfire_mortality,units='[]',ldims=dim2p,gdims=dim2,dimnames=dim2n,           &
             longname='mortality due to crown kill and cambial kill after fire',lpost=.FALSE.,code=77,laccu=.FALSE.,              &
             contnorest=.TRUE.)

    ! "cpool" variables in the veg stream
    CALL add(IO_veg,'frac_1hr_wood',state_FTH%frac_1hr_wood, code=65, units='', longname='fraction of 1hr fuel in wood pool',     &
             contnorest=.TRUE.,lmiss=.TRUE.,missval=mv,ldims=dim2p,gdims=dim2,dimnames=dim2n,lpost=.TRUE.,laccu=.FALSE.,          &
             lrerun=.TRUE.)
    CALL add(IO_veg,'frac_10hr_wood',state_FTH%frac_10hr_wood, code=66, units='', longname='fraction of 10hr fuel in wood pool',  &
             contnorest=.TRUE.,lmiss=.TRUE.,missval=mv,ldims=dim2p,gdims=dim2,dimnames=dim2n,lpost=.TRUE.,laccu=.FALSE.,          &
             lrerun=.TRUE.)
    CALL add(IO_veg,'frac_100hr_wood',state_FTH%frac_100hr_wood, code=67,units='',longname='fraction of 100hr fuel in wood pool', &
             contnorest=.TRUE., lmiss=.TRUE.,missval=mv,ldims=dim2p,gdims=dim2,dimnames=dim2n,lpost=.TRUE.,laccu=.FALSE.,         &
             lrerun=.TRUE.)
    CALL add(IO_veg,'frac_1000hr_wood',state_FTH%frac_1000hr_wood,code=68,longname='fraction of 1000hr fuel in wood pool',        &
             contnorest=.TRUE.,lmiss=.TRUE.,missval=mv,ldims=dim2p,gdims=dim2,dimnames=dim2n,units='',lpost=.TRUE.,    &
             laccu=.FALSE.,lrerun=.TRUE.)

    ! module variables without output
    CALL add(IO_veg,'IR_ROS',fire_TH_diag%IR_ROS,longname='reaction intensity rate of spread',                            &
             units='[kJ m-2 min-1]',ldims=dim1p,gdims=dim1,dimnames=dim1n,lpost=.FALSE.,laccu=.TRUE., contnorest=.TRUE.)
    CALL add(IO_veg,'mortality_crownkill',state_FTH%mortality_crownkill, longname='fraction of trees killed by crown fires',      &
             lrerun=.FALSE.,lmiss=.TRUE.,missval=mv,units='[]',ldims=dim2p,gdims=dim2,dimnames=dim2n,lpost=.FALSE.,laccu=.FALSE.)
    CALL add(IO_veg,'burned_Cpool_reserve',state_FTH%burned_Cpool_reserve,longname='reserve pool carbon burned',                  &
             lrerun=.FALSE.,lmiss=.TRUE.,missval=mv,units='[]',ldims=dim2p,gdims=dim2,dimnames=dim2n,lpost=.FALSE.,laccu=.FALSE.)
    CALL add(IO_veg,'burned_Cpool_wood',state_FTH%burned_Cpool_wood,longname='wood pool carbon burned',                           &
             lrerun=.FALSE.,lmiss=.TRUE.,missval=mv,units='[]',ldims=dim2p,gdims=dim2,dimnames=dim2n,lpost=.FALSE.,laccu=.FALSE.)

    ! Initial values
    state_FTH%frac_1hr_wood   (:,:) = frac_1hr_wood_new
    state_FTH%frac_10hr_wood  (:,:) = frac_10hr_wood_new
    state_FTH%frac_100hr_wood (:,:) = frac_100hr_wood_new
    state_FTH%frac_1000hr_wood(:,:) = frac_1000hr_wood_new
    state_FTH%NI_acc          (:)   = 0._dp
    state_FTH%burned_Cpool_litter_wood_ag (:,:) = 0._dp
    state_FTH%burned_Cpool_litter_green_ag(:,:) = 0._dp
    state_FTH%burned_Cpool_wood           (:,:) = 0._dp
    state_FTH%burned_Cpool_green          (:,:) = 0._dp
    state_FTH%burned_Cpool_reserve        (:,:) = 0._dp

    IF (fire_TH_options%ldiag) THEN
      fire_TH_diag%avg_FDI               (:)   = 0._dp
      fire_TH_diag%avg_ROS               (:)   = 0._dp
      fire_TH_diag%avg_fuel_1hr_total    (:)   = 0._dp
      fire_TH_diag%avg_fuel_10hr_total   (:)   = 0._dp
      fire_TH_diag%avg_fuel_100hr_total  (:)   = 0._dp
      fire_TH_diag%avg_numfire           (:)   = 0._dp
      fire_TH_diag%avg_HumanIgnition     (:)   = 0._dp
      fire_TH_diag%avg_LightningIgnition (:)   = 0._dp
      fire_TH_diag%avg_FireDuration      (:)   = 0._dp
      fire_TH_diag%avg_fire_intensity    (:)   = 0._dp
      fire_TH_diag%avg_postfire_mortality(:,:) = 0._dp
      fire_TH_diag%avg_vegetation_height (:,:) = 0._dp
    ENDIF
    if (fire_TH_options%lcalc_frp) THEN
      fire_TH_diag%avg_FRP_gridcell     (:) = 0._dp
      fire_TH_diag%avg_numfires_gridcell(:) = 0._dp
    ENDIF

  END SUBROUTINE init_fire_thonicke_memory

! calculation of burned fraction due to wildfires after Thonicke et al. 2010 BG (TH2010)
  SUBROUTINE burned_frac_thonicke(lctlib,nidx,kidx0,kidx1,ntiles,with_yasso,                 &
                                  used_pft,is_present,is_glacier,                            &
                                  veg_fract_correction,veg_ratio_max,cover_fract,cover_type, &
                                  Cpool_litter_green_ag,Cpool_litter_wood_ag,                &
                                  Cpool_green,Cpool_woods,Cpool_reserve,                     &
                                  frac_litter_wood_new,                                      &
                                  WindSpeed,SoilMoisture,temp_max,temp_min,precip,           &
                                  YCpool_acid_ag1, YCpool_acid_ag2, YCpool_water_ag1,        &
                                  YCpool_water_ag2,YCpool_ethanol_ag1,YCpool_ethanol_ag2,    &
                                  YCpool_nonsoluble_ag1, YCpool_nonsoluble_ag2,              &
                                  fuel,burned_frac_diag,mortality)
  USE mo_jsbach_constants, ONLY: pi, molarMassC_kg, C_2_biomass_ratio, secperday
  USE mo_cbal_cpools,      ONLY: frac_green_aboveGround
  USE mo_land_surface,     ONLY: fract_small
  USE mo_time_control,     ONLY: lstart, lresume

! Input parameters
  TYPE (lctlib_type), INTENT(in)    :: lctlib
  INTEGER,            INTENT(in)    :: nidx, kidx0, kidx1         ! vector length, indeces for global field of NI_acc  
  INTEGER,            INTENT(in)    :: ntiles                     ! #tiles (different land cover type mosaics per gridcell)
  LOGICAL,            INTENT(in)    :: with_yasso
  LOGICAL,            INTENT(in)    :: used_pft             (:)   ! PFT's which are allowed to burn
  LOGICAL,            INTENT(in)    :: is_present           (:,:) ! determines, if tile is handled anyway
  LOGICAL,            INTENT(in)    :: is_glacier           (:,:) ! determines, if tile is handled anyway
  INTEGER,            INTENT(in)    :: cover_type           (:,:)
  REAL(dp),           INTENT(in)    :: veg_fract_correction (:,:)
  REAL(dp),           INTENT(in)    :: veg_ratio_max        (:)   ! fractional cover of vegetated area in a grid box
  REAL(dp),           INTENT(in)    :: cover_fract          (:,:) ! fractional cover of all cover types on all tiles
  REAL(dp),           INTENT(in)    :: Cpool_litter_green_ag(:,:) ! above ground green litter carbon pool [mol(C)/m2(canopy)]
  REAL(dp),           INTENT(in)    :: Cpool_litter_wood_ag (:,:) ! Value of wood litter carbon pool [mol(C)/m^2(canopy)]
  REAL(dp),           INTENT(in)    :: Cpool_green          (:,:) ! live green carbon pool
  REAL(dp),           INTENT(in)    :: Cpool_woods          (:,:) ! live wood carbon pool
  REAL(dp),           INTENT(in)    :: Cpool_reserve        (:,:) ! live reserve carbon pool
  REAL(dp),           INTENT(in)    :: frac_litter_wood_new (:,:) ! frac of wood litter pool that is new since the last call
  REAL(dp),           INTENT(in)    :: WindSpeed            (:)   ! wind speed [m/s]
  REAL(dp),           INTENT(in)    :: SoilMoisture         (:)   ! top layer soil moisture
  REAL(dp),           INTENT(in)    :: temp_min(:),temp_max (:)   ! min and max temperatures
  REAL(dp),           INTENT(in)    :: precip               (:)   ! precipitation
! Yasso pools
  REAL(dp),           INTENT(inout) :: YCpool_acid_ag1      (:,:), YCpool_acid_ag2      (:,:)
  REAL(dp),           INTENT(inout) :: YCpool_water_ag1     (:,:), YCpool_water_ag2     (:,:)
  REAL(dp),           INTENT(inout) :: YCpool_ethanol_ag1   (:,:), YCpool_ethanol_ag2   (:,:)
  REAL(dp),           INTENT(inout) :: YCpool_nonsoluble_ag1(:,:), YCpool_nonsoluble_ag2(:,:)
!Output
  REAL(dp),           INTENT(inout) :: fuel                 (:)   ! total aboveground litter  [g(C)/m^2(grid box)]
  REAL(dp),           INTENT(inout) :: burned_frac_diag     (:,:) ! burned fraction per tile
  REAL(dp),           INTENT(inout) :: mortality            (:,:) ! postfire mortality per tile

!local variables
  REAL(dp) :: fuel_tiled(nidx,ntiles)                           ! fuel per tile, used to compute fuel,[g(C)/m2(canopy)]
  REAL(dp) :: moist_tiled(nidx,ntiles)                          ! moisture factor per tile
  REAL(dp) :: Cpool_litter_green_ag_g(nidx,ntiles)              ! above ground green litter carbon pool [g(C)/m2(canopy)]
  REAL(dp) :: Cpool_litter_wood_ag_g(nidx,ntiles)               ! Value of wood litter carbon pool [g(C)/m^2(canopy)]
  REAL(dp) :: Cpool_green_g(nidx,ntiles)                        ! live green carbon pool in [g(C)/m^2(canopy)]
  REAL(dp) :: fuel_1hr(nidx,ntiles), fuel_10hr(nidx,ntiles)     ! different fuel classes  [g(C,biomass)/m2(gridcell)]
  REAL(dp) :: fuel_100hr(nidx,ntiles),fuel_1000hr(nidx,ntiles)  ! gC -> g biomass: 1/C_2_biomass_ratio
  REAL(dp) :: moistfactor(nidx),Length2BreadthRatio(nidx)       ! moisture factor,len to breadth ratio of burned ellipse
  REAL(dp) :: fuel_1hr_total(nidx), fuel_10hr_total(nidx)       ! fuel in [g(biomass)/m2(gridcell)]
  REAL(dp) :: fuel_100hr_total(nidx)                            ! accumulated over used pfts
  REAL(dp) :: canopy2activeveg(nidx,ntiles)                     ! conversion from [g(C)/m2(canopy)] to [g(C)/m2] active vegetation
  REAL(dp) :: Fuel1To100hrSum(nidx)                             ! sum over 1hr,10hr,100hr fuel_totals, [g(biomass)/m2(gridcell)]
  REAL(dp) :: net_fuel(nidx)                                    ! kg(biomass)/m2(gridcell), reduced for mineral content->MINER_TOT
  REAL(dp) :: livegrass(nidx),LiveBiomassTiled(nidx,ntiles)            ! live grass in [g(biomass?)/m2],live biomass 
  REAL(dp) :: sum_grass_fract(nidx),sum_wood_fract(nidx)        ! summed fpc of grasses and woody species
  REAL(dp) :: WindSpeed_m_min(nidx),WindSpeed_ft_min(nidx)      ! wind speed in [m/min] and [ft/min]
  REAL(dp) :: ratio_C3_livegrass(nidx),ratio_C4_livegrass(nidx) ! ratio of C3,C4 grass /tot grass biomass
  REAL(dp) :: ratio_fbd, dens_fuel_ave(nidx)                    ! fuel bulk dens weighted per fuel class,averaged over tiles
  REAL(dp) :: ratio_dead_fuel(nidx), ratio_live_fuel(nidx)      ! fraction of dead and live fuel
  REAL(dp) :: char_dens_fuel_ave(nidx), char_net_fuel(nidx)     ! average fuel density (incl. livegrass),net fuel incl. live grass
                                                                !   and reduced for mineral compounds [kg(biomass)]
  REAL(dp) :: SurfArea2Vol(nidx), char_moistfactor(nidx)        ! surface area to volume ratio [cm-1], moisture of extinction 
                                                                !   averaged over available fuels
  REAL(dp) :: RelFuelMoisture(nidx)                                         ! relative fuel moisture eq. 6 in TH2010
  REAL(dp) :: leaf_moisture(nidx)                               ! relative leaf moisture of grasses
  REAL(dp) :: FDI(nidx), numfire(nidx)                          ! Fire Danger Index, number of fires per [ha-1]
  REAL(dp) :: ROS(nidx), gamma(nidx)                            ! rate of spread[m min-1] and gamma for postfire mortality
  REAL(dp) :: ROS_b(nidx), FireDuration(nidx)                     ! backward ROS[m min-1] and fire duration [min]
  REAL(dp) :: LengthOfMajAxis(nidx)                                ! Length of major axis of the burned ellipse [m]
                                                                ! ROS*FireDuration
  REAL(dp) :: active_fractions(nidx,ntiles)                     ! fractions of fire active pfts(only crops excluded, sum up to one)
  REAL(dp) :: fuel_consum(nidx) ! fuel_consumption of 1,10,100hr fuels per m2 vegetated area
  REAL(dp) :: fuel_consum1000perAreaBurned(nidx),fuel_consum_tiled(nidx,ntiles) 
  REAL(dp) :: com_fac(nidx,ntiles), com_fac_lg(nidx,ntiles)     ! averaged combustion factor
  REAL(dp) :: com_fac_1hr(nidx,ntiles),   com_fac_10hr(nidx,ntiles)
  REAL(dp) :: com_fac_100hr(nidx,ntiles), com_fac_1000hr(nidx,ntiles) ! combustion factors computed in fuel_consumption
  REAL(dp) :: vegetation_height(nidx, ntiles)
  REAL(dp) :: scorch_height(nidx, ntiles)                       ! scorch height in [m]
  REAL(dp) :: mortality_camb(nidx,ntiles)
  REAL(dp) :: I_surface(nidx)                                   ! fire surface intensity
  LOGICAL  :: used_pft_array(nidx,ntiles)                       ! used_pft extended to array
  LOGICAL  :: present_and_used(nidx,ntiles), woody_pft(nidx, ntiles) ! combination of is_present and used_pft
  REAL(dp) :: dens_livegrass_ave(nidx)                          ! fuel density of livegrass weighted average
  REAL(dp) :: human_i(nidx), lightning_i(nidx)                  ! human ignitions,lightning ignitions,[km-2] sum up to pot. fires
  REAL(dp) :: FRP_gridcell(nidx),numfires_gridcell(nidx)                    ! fire radiative power and number of fires per gridcell
  INTEGER  :: i,itile, ki, n                                    ! ki: index equivalent to i but shifted for kidx0

    IF (fire_TH_options%read_fuel_frac .AND. (lstart .OR. lresume)) THEN
      state_FTH%frac_1hr_wood        (kidx0:kidx1,:) = init_state_FTH%frac_1hr_wood        (kidx0:kidx1,:)
      state_FTH%frac_10hr_wood       (kidx0:kidx1,:) = init_state_FTH%frac_10hr_wood       (kidx0:kidx1,:)
      state_FTH%frac_100hr_wood      (kidx0:kidx1,:) = init_state_FTH%frac_100hr_wood      (kidx0:kidx1,:)
      state_FTH%frac_1000hr_wood     (kidx0:kidx1,:) = init_state_FTH%frac_1000hr_wood     (kidx0:kidx1,:)
      state_FTH%NI_acc               (kidx0:kidx1)   = init_state_FTH%NI_acc               (kidx0:kidx1)
    ENDIF

    ! calculate cloud-to-ground flash density from total flash density
    IF (InputVarIsUpdated(Lightning_var) .AND. .NOT. fire_TH_options%llight_ground) THEN
      fire_TH_diag%lightning(kidx0:kidx1) = 0.1_dp*(1._dp+fire_TH_diag%lat(kidx0:kidx1)*fire_TH_diag%lat(kidx0:kidx1)/900._dp) &
                                          * fire_TH_diag%lightning(kidx0:kidx1)
    ENDIF

    ! Cpools not yet initialized (first step): no sensible fire can be calculated
    IF (ANY((fuel <= -1.e36_dp) .OR. (fuel >= 1.e36_dp))) THEN
      fuel       (:)   = 0._dp
      burned_frac_diag(:,:) = 0._dp
      mortality  (:,:) = 1._dp
      RETURN
    ENDIF

    ! initializations
    sum_grass_fract   (:)   = 0._dp
    sum_wood_fract    (:)   = 0._dp
    dens_fuel_ave     (:)   = 0._dp
    ratio_dead_fuel   (:)   = 0._dp
    ratio_live_fuel   (:)   = 0._dp
    char_dens_fuel_ave(:)   = 0._dp
    livegrass         (:)   = 0._dp
    Fuel1To100hrSum   (:)   = 0._dp

    ! Masks and conversions:
    !    conversion from C m-2 canopy to m-2 active vegetation
    DO i=1,ntiles
      used_pft_array(:,i) = used_pft(cover_type(:,i)) .AND. (COUNT(is_glacier(:,:),DIM=2) == 0)
      woody_pft     (:,i) = lctlib%woody_pft(cover_type(:,i))
    ENDDO
    active_fractions    = cover_fract
    WHERE (.NOT. used_pft_array(:,:))
      active_fractions(:,:) = 0._dp
    ENDWHERE
    active_fractions(:,:) = active_fractions(:,:) / MAX(1e-10_dp,SPREAD(SUM(active_fractions(:,:),2),DIM=2,NCOPIES=ntiles))
    canopy2activeveg(:,:) = veg_fract_correction(:,:) * active_fractions(:,:)

    !    combines used_pft_array to not include crops and is_present, 
    present_and_used = is_present(1:nidx,1:ntiles) .AND. used_pft_array(1:nidx,1:ntiles)

    ! convert carbon pools from mol to g
    IF (with_yasso) THEN
      Cpool_litter_green_ag_g(:,:) = (  YCpool_acid_ag1   (:,:) + YCpool_water_ag1     (:,:) & 
                                      + YCpool_ethanol_ag1(:,:) + YCpool_nonsoluble_ag1(:,:) &
                                     ) * molarMassC_kg * 1000._dp
      Cpool_litter_wood_ag_g (:,:) = (  YCpool_acid_ag2   (:,:) + YCpool_water_ag2     (:,:) & 
                                      + YCpool_ethanol_ag2(:,:) + YCpool_nonsoluble_ag2(:,:) &
                                     ) * molarMassC_kg * 1000._dp
    ELSE
      Cpool_litter_green_ag_g(:,:) = Cpool_litter_green_ag(:,:) * molarMassC_kg * 1000._dp
      Cpool_litter_wood_ag_g (:,:) = Cpool_litter_wood_ag (:,:) * molarMassC_kg * 1000._dp
    ENDIF
    Cpool_green_g          (:,:) = Cpool_green          (:,:) * molarMassC_kg * 1000._dp

    ! Divide carbon pools into fuel classes
    !   total grass green litter pool and a fraction of woody pft green litter pool is 1hr fuel
    !   for grasses also live green pool adds to 1 hr fuel, crops will be filtered in average tiles
    !   10hour fuels: 7.5 % of hard and sapwood, 100h: 21% and 1000h: 67% of hard and sapwood
    !   and initialize woody_pft

    ! update fuel fractions
    state_FTH%frac_1hr_wood   (kidx0:kidx1,:) = state_FTH%frac_1hr_wood   (kidx0:kidx1,:) * (1._dp-frac_litter_wood_new (:,:)) &
                                              + frac_1hr_wood_new                         *        frac_litter_wood_new (:,:)
    state_FTH%frac_10hr_wood  (kidx0:kidx1,:) = state_FTH%frac_10hr_wood  (kidx0:kidx1,:) * (1._dp-frac_litter_wood_new (:,:)) &
                                              + frac_10hr_wood_new                        *        frac_litter_wood_new (:,:)
    state_FTH%frac_100hr_wood (kidx0:kidx1,:) = state_FTH%frac_100hr_wood (kidx0:kidx1,:) * (1._dp-frac_litter_wood_new (:,:)) &
                                              + frac_100hr_wood_new                       *        frac_litter_wood_new (:,:)
    state_FTH%frac_1000hr_wood(kidx0:kidx1,:) = state_FTH%frac_1000hr_wood(kidx0:kidx1,:) * (1._dp-frac_litter_wood_new (:,:)) &
                                              + frac_1000hr_wood_new                      *        frac_litter_wood_new (:,:)

    ! Consistency checks
    n = 0
    DO itile=1,ntiles
      DO i=1,nidx
        ki=i+kidx0-1
        IF (Cpool_litter_green_ag(i,itile) < 0._dp) THEN
          WRITE(message_text,*) 'Error: Cpool_litter_green_ag < 0. Point (',ki,itile,')=',Cpool_litter_green_ag(i,itile)
          CALL message('burned_frac_thonicke',TRIM(message_text))
          n = n + 1
        ENDIF
        IF ((state_FTH%frac_1hr_wood(ki,itile) < 0._dp) .OR. (state_FTH%frac_1hr_wood(ki,itile) > 1._dp)) THEN
          WRITE(message_text,*) 'Error: frac_1hr_wood not in range 0-1. Point (',ki,itile,')=',state_FTH%frac_1hr_wood(ki,itile)
          CALL message('burned_frac_thonicke',TRIM(message_text))
          n = n + 1
        ENDIF
        IF ((state_FTH%frac_10hr_wood(ki,itile) < 0._dp) .OR. (state_FTH%frac_10hr_wood(ki,itile) > 1._dp)) THEN
          WRITE(message_text,*) 'Error: frac_10hr_wood not in range 0-1. Point (',ki,itile,')=',state_FTH%frac_10hr_wood(ki,itile)
          CALL message('burned_frac_thonicke',TRIM(message_text))
          n = n + 1
        ENDIF
        IF ((state_FTH%frac_100hr_wood(ki,itile) < 0._dp) .OR. (state_FTH%frac_100hr_wood(ki,itile) > 1._dp)) THEN
          WRITE(message_text,*) 'Error: frac_100hr_wood not in range 0-1. Point (',ki,itile,')=',state_FTH%frac_100hr_wood(ki,itile)
          CALL message('burned_frac_thonicke',TRIM(message_text))
          n = n + 1
        ENDIF
        IF ((state_FTH%frac_1000hr_wood(ki,itile) < 0._dp) .OR. (state_FTH%frac_1000hr_wood(ki,itile) > 1._dp)) THEN
          WRITE(message_text,*) &
            'Error: frac_1000hr_wood not in range 0-1. Point (',ki,itile,')=',state_FTH%frac_1000hr_wood(ki,itile)
          CALL message('burned_frac_thonicke',TRIM(message_text))
          n = n + 1
        ENDIF
        IF (n > 100) CALL finish('burned_frac_thonicke','Inconsistent fuel fractions')
      ENDDO
    ENDDO

    ! Calculate fuel load in g(biomass)/m2(vegetation)
    DO itile=1,ntiles
      DO i=1,nidx
        ki=i+kidx0-1
        IF (woody_pft(i,itile)) THEN
          fuel_1hr   (i,itile) = ( Cpool_litter_green_ag_g(i,itile)*frac_green_active + & 
                                   state_FTH%frac_1hr_wood   (ki,itile) * Cpool_litter_wood_ag_g (i,itile)             &
                                  )* canopy2activeveg(i,itile) / C_2_biomass_ratio
          fuel_10hr  (i,itile) =    state_FTH%frac_10hr_wood  (ki,itile) * Cpool_litter_wood_ag_g (i,itile)             &
                                  * canopy2activeveg(i,itile) / C_2_biomass_ratio
          fuel_100hr (i,itile) =    state_FTH%frac_100hr_wood (ki,itile) * Cpool_litter_wood_ag_g (i,itile)             &
                                  * canopy2activeveg(i,itile) / C_2_biomass_ratio
          fuel_1000hr(i,itile) =    state_FTH%frac_1000hr_wood(ki,itile) * Cpool_litter_wood_ag_g (i,itile)             &
                                  * canopy2activeveg(i,itile) / C_2_biomass_ratio
        ELSE
          IF (used_pft_array(i,itile)) THEN
            fuel_1hr (i,itile) = Cpool_litter_green_ag_g(i,itile) * canopy2activeveg(i,itile) / C_2_biomass_ratio
          ELSE
            fuel_1hr (i,itile) = 0._dp
          ENDIF
          fuel_10hr  (i,itile) = 0._dp
          fuel_100hr (i,itile) = 0._dp
          fuel_1000hr(i,itile) = 0._dp
        ENDIF
      ENDDO
    ENDDO

    ! Fuel aggregation
    fuel_1hr_total  (:) = SUM(fuel_1hr  (:,:),DIM=2,MASK=present_and_used(:,:))
    fuel_10hr_total (:) = SUM(fuel_10hr (:,:),DIM=2,MASK=present_and_used(:,:))
    fuel_100hr_total(:) = SUM(fuel_100hr(:,:),DIM=2,MASK=present_and_used(:,:))
    Fuel1To100hrSum (:) = fuel_1hr_total(:) + fuel_10hr_total(:) + fuel_100hr_total(:) ! 1000 hr not included here.

    ! live biomass per tile, needed only for grasses 
    LiveBiomassTiled(:,:) = Cpool_green_g(:,:) * frac_green_aboveGround * canopy2activeveg(:,:) / C_2_biomass_ratio

    ! amount of live grass, only above ground
    livegrass(:) = SUM(LiveBiomassTiled(:,:),DIM=2,MASK=present_and_used(:,:) .AND. (.NOT. woody_pft(:,:)))

    ! net fuel in kg(biomass)/m2
    net_fuel(:) = (1.0_dp - MINER_TOT) * (Fuel1To100hrSum(:) / 1000.0_dp)

    ! compute above ground litter converted to g/m2 gridcell
    fuel_tiled(:,:) = fuel_1hr(:,:) + fuel_10hr(:,:) + fuel_100hr(:,:) + fuel_1000hr(:,:)

    ! sum available fuel over tiles
    fuel(:) = SUM(fuel_tiled(:,:),DIM=2,MASK=present_and_used(:,:)) 

    ! calculate litter moisture weighting factor (moisture of extinction me) ; here 1000hr fuel is included
    !   pft dependend moisture of extinction weighted with the tile specific contribution to the available fuel
    DO i=1,nidx
      moist_tiled(i,:) = 0.0_dp
      IF (fuel(i) > 0.0_dp) moist_tiled(i,:) = fuel_tiled(i,:) / fuel(i) * lctlib%moist_extinction(cover_type(i,:))
    ENDDO
    moistfactor(:) = SUM(moist_tiled(:,:),DIM=2,MASK=present_and_used(:,:))
 
    ! Calculation for fire perimeter, calculate tree and grass coverage. One could also account for veg ratio max here  
    sum_grass_fract(:) = SUM(cover_fract(:,:),DIM=2,MASK=(.NOT. woody_pft(:,:)))
    sum_wood_fract (:) = SUM(cover_fract(:,:),DIM=2,MASK=(      woody_pft(:,:)))

    ! wind speed in m/min needed instead of m/s and reduction with respect to grass_fpc (0.6) and wood_fpc (0.4) (TH2010 p.1997)
    WindSpeed_m_min (:) = (  (sum_wood_fract (:) * 0.4_dp) &
                           + (sum_grass_fract(:) * 0.6_dp) &
                          ) * WindSpeed(:) * 60._dp
    WindSpeed_ft_min(:) = WindSpeed_m_min(:) * ft_per_m ! wind speed in ft/min for input into Rothermel's formula

    ! length to breadth ratio of burned ellipse according to presence of grasses and forests eq. 12 and 13 in TH2010
    ! in paper units of wind in eq 12,13 is m/min in code factor 0.06 included --> km/h, in original reference km/h is correct
    WHERE (WindSpeed_m_min(:) < 16.67) 
      Length2BreadthRatio(:) = 1.0_dp
    ELSEWHERE
      Length2BreadthRatio(:) = MAX(MIN(sum_wood_fract(:)                                                                     &
                                       * (1.0_dp+(8.729_dp*((1.0_dp-(EXP(-0.03_dp*0.06_dp*WindSpeed_m_min(:))))**2.155_dp))) &
                                       + sum_grass_fract(:) * (1.1_dp+((0.06_dp*WindSpeed_m_min(:))**0.0464_dp))             &
                                      ,8._dp),0.001_dp)
    ENDWHERE

    DO itile=1,ntiles
      DO i=1,nidx
        ki=i+kidx0-1
        IF (is_present(i,itile) .AND. used_pft_array(i,itile)) THEN
          ! transform fuel bulk densities of 10 and 100 hr fuel types to 1 hr equivalents
          IF (Fuel1To100hrSum(i) > fract_small) THEN
            ratio_fbd = (  fuel_1hr  (i,itile)                  &
                         + fuel_10hr (i,itile) * fbd10hr_2_1hr  &
                         + fuel_100hr(i,itile) * fbd100hr_2_1hr &
                        ) / Fuel1To100hrSum(i)
            dens_fuel_ave (i)       = dens_fuel_ave(i) + ratio_fbd * lctlib%fuel_dens(cover_type(i,itile)) ! * cover_fract(i,itile)
            IF ((ratio_fbd > 1.4_dp) .OR. (dens_fuel_ave(i) > 1.4_dp*25._dp)) THEN
              WRITE(message_text,*) 'for point ',ki,' values: ',ratio_fbd, dens_fuel_ave(i)
              CALL message('burned_frac_thonicke','Error: ratio_fbd or dens_fuel_ave out of range '//TRIM(message_text))
            ENDIF
          ELSE
            dens_fuel_ave(i) = 0._dp
          ENDIF

          ! and compute the ratio of C4 grass biomass to total grass biomass
          IF (lctlib%C4flag(cover_type(i,itile)) .AND. (livegrass(i) > fract_small) .AND. used_pft_array(i,itile)) THEN
            ratio_C4_livegrass(i) = (Cpool_green_g(i,itile) / C_2_biomass_ratio & 
                                  * canopy2activeveg(i,itile) * frac_green_aboveGround) / livegrass(i)
          ELSE
            ratio_C4_livegrass(i) = 0.0_dp
          ENDIF
        ELSE 
          ratio_C4_livegrass(i) = 0.0_dp
        ENDIF
      ENDDO
    ENDDO

    ! ratio of C3 grass biomass to toal grass biomass
    WHERE(livegrass(:) > fract_small)
      ratio_C3_livegrass(:) = 1._dp - ratio_C4_livegrass(:)
    ELSEWHERE
      ratio_C3_livegrass(:) = 0.0_dp
    ENDWHERE

    ! leaf moisture content of grasses (RelFuelMoisture_lg in spitfire) eq. B2
    leaf_moisture(:) = MAX(0._dp,((10._dp/9._dp) * SoilMoisture(:) - (1._dp/9._dp)))

    !      influence of livegrass on FBD
    !    Kirsten in lpj code: presence of livegrass additionally lowers FBD, when trees mixed 
    !    with grasses fbd for C3 and C4 = 4
    dens_livegrass_ave(:) = fbd_C3_livegrass * ratio_C3_livegrass(:) &
                          + fbd_C4_livegrass * ratio_C4_livegrass(:)

    DO i=1,nidx
      IF ((Fuel1To100hrSum(i) > fract_small) .OR. livegrass(i) > fract_small) THEN 
        ratio_dead_fuel   (i) = Fuel1To100hrSum(i) / (Fuel1To100hrSum(i) + livegrass(i))
        ratio_live_fuel   (i) = livegrass    (i) / (Fuel1To100hrSum(i) + livegrass(i))
        char_dens_fuel_ave(i) = dens_fuel_ave(i) * ratio_dead_fuel(i) + dens_livegrass_ave(i) * ratio_live_fuel(i)
      ENDIF
    ENDDO

    ! combining dead net fuel and net fuel of livegrass
    char_net_fuel   (:) = net_fuel(:) + (1.0 - MINER_TOT) * livegrass(:) / 1000.0 ! in kg biomass
    char_moistfactor(:) = moistfactor(:) * ratio_dead_fuel(:) + moistfactor_livegrass * ratio_live_fuel(:)

    WHERE (Fuel1To100hrSum(:) > fract_small)
      SurfArea2Vol(:) = (  fuel_1hr_total  (:) * fire_TH_options%SurfArea2Vol(1) &
                         + fuel_10hr_total (:) * fire_TH_options%SurfArea2Vol(2) &
                         + fuel_100hr_total(:) * fire_TH_options%SurfArea2Vol(3) &
                         ) / Fuel1To100hrSum(:)
    ELSEWHERE
      SurfArea2Vol(:) = 0.00001_dp
    END WHERE

    ! end preparation of variables
 
    ! Calculate fire dange index
    CALL fire_danger_index(nidx,fuel_1hr_total,fuel_10hr_total,fuel_100hr_total,leaf_moisture,char_moistfactor,  &
                           Fuel1To100hrSum,ratio_live_fuel,ratio_dead_fuel,fire_TH_options%SurfArea2Vol,&
                           temp_max,temp_min,precip,RelFuelMoisture, &
                           FDI,fire_TH_options%moisture_scaling,state_FTH%NI_acc(kidx0:kidx1))

    ! Calculate number of fires
    ! lightning was converted with a factor 0.2 for cloud to ground flashes and 0.003 for the efficiency to start a fire
    ! updated: now cloud to ground ratio depends on the latitude
    ! here its the density of fires instead of the absolute number to avoid the usage of the grid cell area
    human_i(:) = human_ign(nidx,fire_TH_diag%pop(kidx0:kidx1),fire_TH_diag%a_nd_array(kidx0:kidx1),fire_TH_options%human_para)

    ! only a fraction of cloud to ground flashes ignites a fire
    lightning_i(:) = fire_TH_diag%lightning(kidx0:kidx1) * 0.04_dp
    WHERE (net_fuel(:) > 0.001)
      numfire(:) = FDI(:) * (human_i(:) + lightning_i(:)) * 0.01_dp * fire_TH_options%ign_para
      ! probability of human ignitions +probability of lightning ignition per ha
    ELSEWHERE  !not enough fuel
      numfire(:) = 0.0_dp
    ENDWHERE

    ! rate of spread, ROS, in m per minute
    CALL rate_of_spread(nidx,WindSpeed_ft_min,SurfArea2Vol,RelFuelMoisture,char_dens_fuel_ave,char_net_fuel,&
                        char_moistfactor,          &
                        ROS,gamma,fire_TH_diag%IR_ROS(kidx0:kidx1),fire_TH_options%lwind_speed_limit,&
                        fire_TH_options%wind_limit, &
                        fire_TH_options%wind_max,fire_TH_options%wind_slope)
 
    ! fire duration as a function of Fire danger index and optinal population density
    IF (fire_TH_options%lduration_popd) THEN
      FireDuration(:) = 241._dp / (1.0_dp + (240.0_dp * exp(-11.06_dp * FDI(:)))) &
                      * (-log10(MIN(MAX(fire_TH_diag%pop(kidx0:kidx1),0.01_dp),100._dp)) + 4_dp) * 0.5_dp
    ELSE
      FireDuration(:) = 241._dp / (1.0_dp + (240._dp * exp(-11.06_dp * FDI(:))))
    ENDIF

    ros_b(:) = ROS(:) * EXP(-0.012_dp * WindSpeed_m_min(:)) ! windspeed in m/min and averaged for grass and trees
    WHERE (ros_b(:) < 0.05)
      ros_b(:) = 0.0_dp 
    ENDWHERE

    ! HIRSCH (1996)
    ! backward rate of spread
    WHERE (net_fuel(:) > fract_small)
      LengthOfMajAxis(:) = (ROS_b(:) + ROS) * FireDuration(:) ! in min
    ELSEWHERE
      numfire        (:) = 0._dp
      LengthOfMajAxis(:) = 0._dp
      FireDuration   (:) = 0._dp
      ros_b          (:) = 0._dp
      ROS            (:) = 0._dp
    ENDWHERE

    ! burned fraction
    DO i=1,nidx
      IF ((Length2BreadthRatio(i) > fract_small) .AND. (ROS(i) > fract_small)) THEN
        burned_frac_diag(i,:) = MIN(1._dp,numfire(i) * ((pi / (4._dp * Length2BreadthRatio(i)))                           &
                                          * ((LengthOfMajAxis(i))**2._dp)) / 10000._dp) ! burned area per ha [m2]/1ha[m2]
      ELSE
        burned_frac_diag(i,:) = 0._dp
      ENDIF
    ENDDO
    WHERE (.NOT. used_pft_array(:,:))
      burned_frac_diag(:,:) = 0._dp
    ENDWHERE

    ! calculate fuel consumption
    CALL  fuel_consumption(nidx,ntiles,used_pft_array,woody_pft,                             &
                           state_FTH%NI_acc(kidx0:kidx1),burned_frac_diag,RelFuelMoisture,   &
                           fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,fuel_1hr_total,         &
                           leaf_moisture,livegrass,char_moistfactor,com_fac,fuel_consum,     &
                           fuel_consum_tiled, fuel_consum1000perAreaBurned,                  &
                           com_fac_1hr,com_fac_10hr,com_fac_100hr,com_fac_1000hr,com_fac_lg)

    ! Fire intensity eq. 15 in TH2010
    WHERE(ANY(burned_frac_diag(:,:) > fract_small,DIM=2) .AND. (fuel_consum(:) > fract_small))
      I_surface(:) = HeatContentFuel * (fuel_consum(:) / 1000._dp / MAXVAL(burned_frac_diag(:,:),DIM=2)) * (ros / 60._dp)
    ELSEWHERE
      I_surface(:) = 0._dp
    ENDWHERE

    ! Apply lower intensity threshold for fires to be accepted
    DO itile=1,ntiles
      WHERE (I_surface(:) < 50._dp)
        burned_frac_diag(:,itile) = 0._dp 
      ENDWHERE
    ENDDO
    WHERE (I_surface(:) < 50._dp)
      I_surface(:) = 0._dp
      numfire  (:) = 0._dp
      fuel_consum(:) = 0._dp
    ENDWHERE
    IF (fire_TH_options%lcalc_frp) THEN
      CALL FRP_nfire_per_gridcell(nidx,ntiles,cover_fract,veg_ratio_max,numfire,burned_frac_diag,                   & 
                                  fire_TH_diag%IR_ROS(kidx0:kidx1),FireDuration,fire_TH_diag%cellarea(kidx0:kidx1), &
                                  FRP_gridcell,numfires_gridcell)    
    ENDIF

    ! vegetation height (Shevlikova et al. 2009, GBC), h_max=24.19 m, hl=0.19 m2/kgC
    vegetation_height(:,:) = 0._dp
    WHERE (woody_pft(:,:))
      vegetation_height(:,:) = hmax * (1.0_dp - EXP(-hl * veg_fract_correction(:,:) &
                                * (Cpool_green_g(:,:) * 1.e-3_dp + (Cpool_woods(:,:) + Cpool_reserve(:,:)) * molarMassC_kg)))
    ENDWHERE

    ! tree mortality after fire
    CALL fire_mortality(lctlib,nidx,ntiles,woody_pft,cover_type,I_surface,&
                        com_fac,fuel_consum1000perAreaBurned,gamma,vegetation_height, &
                        fire_TH_options%Lethal2ResTime,state_FTH%mortality_crownkill(kidx0:kidx1,:), &
                        mortality_camb,scorch_height,        &
                        state_FTH%postfire_mortality (kidx0:kidx1,:))
    ! mortality returned to mo_disturbance to correct burned_frac_diag to burned_frac for dynveg
    mortality(:,:) = state_FTH%postfire_mortality(kidx0:kidx1,:)
    IF (with_yasso) THEN
      CALL pool_consumption_yasso(nidx,ntiles,woody_pft,used_pft_array,burned_frac_diag,           &
                            com_fac_1hr,com_fac_10hr,com_fac_100hr,com_fac_1000hr,com_fac_lg, &
                            Cpool_green,                                                      &
                            YCpool_acid_ag1,       YCpool_acid_ag2,       YCpool_water_ag1,   &
                            YCpool_water_ag2,      YCpool_ethanol_ag1,    YCpool_ethanol_ag2, &
                            YCpool_nonsoluble_ag1, YCpool_nonsoluble_ag2,                     &
                            state_FTH%frac_1hr_wood               (kidx0:kidx1,:),            &
                            state_FTH%frac_10hr_wood              (kidx0:kidx1,:),            &
                            state_FTH%frac_100hr_wood             (kidx0:kidx1,:),            &
                            state_FTH%frac_1000hr_wood            (kidx0:kidx1,:),            &
                            state_FTH%burned_YC_acid_ag1          (kidx0:kidx1,:),            &
                            state_FTH%burned_YC_acid_ag2          (kidx0:kidx1,:),            &
                            state_FTH%burned_YC_water_ag1         (kidx0:kidx1,:),            &
                            state_FTH%burned_YC_water_ag2         (kidx0:kidx1,:),            &
                            state_FTH%burned_YC_ethanol_ag1       (kidx0:kidx1,:),            &
                            state_FTH%burned_YC_ethanol_ag2       (kidx0:kidx1,:),            &
                            state_FTH%burned_YC_nonsoluble_ag1    (kidx0:kidx1,:),            &
                            state_FTH%burned_YC_nonsoluble_ag2    (kidx0:kidx1,:),            &
                            state_FTH%burned_Cpool_litter_wood_ag (kidx0:kidx1,:),            &
                            state_FTH%burned_Cpool_litter_green_ag(kidx0:kidx1,:),            &
                            state_FTH%burned_Cpool_green          (kidx0:kidx1,:))
    ELSE
      CALL pool_consumption(nidx,ntiles,woody_pft,used_pft_array,burned_frac_diag,                 &
                            com_fac_1hr,com_fac_10hr,com_fac_100hr,com_fac_1000hr,com_fac_lg, &
                            Cpool_litter_wood_ag, Cpool_litter_green_ag, Cpool_green,         &
                            state_FTH%frac_1hr_wood               (kidx0:kidx1,:),            &
                            state_FTH%frac_10hr_wood              (kidx0:kidx1,:),            &
                            state_FTH%frac_100hr_wood             (kidx0:kidx1,:),            &
                            state_FTH%frac_1000hr_wood            (kidx0:kidx1,:),            &
                            state_FTH%burned_Cpool_litter_wood_ag (kidx0:kidx1,:),            &
                            state_FTH%burned_Cpool_litter_green_ag(kidx0:kidx1,:),            &
                            state_FTH%burned_Cpool_green          (kidx0:kidx1,:))
    ENDIF
    ! These don't have any "memory" across subroutines
    state_FTH%burned_Cpool_reserve(kidx0:kidx1,:) = 0._dp
    state_FTH%burned_Cpool_wood   (kidx0:kidx1,:) = 0._dp

    ! Cumulating diagnostics for output
    IF (fire_TH_options%ldiag) THEN
      fire_TH_diag%avg_fuel_1hr_total     (kidx0:kidx1)   = &
      fire_TH_diag%avg_fuel_1hr_total     (kidx0:kidx1)   + secperday * fuel_1hr_total
      fire_TH_diag%avg_fuel_10hr_total    (kidx0:kidx1)   = &
      fire_TH_diag%avg_fuel_10hr_total    (kidx0:kidx1)   + secperday * fuel_10hr_total
      fire_TH_diag%avg_fuel_100hr_total   (kidx0:kidx1)   = &
      fire_TH_diag%avg_fuel_100hr_total   (kidx0:kidx1)   + secperday * fuel_100hr_total
      fire_TH_diag%avg_FDI                (kidx0:kidx1)   = &
      fire_TH_diag%avg_FDI                (kidx0:kidx1)   + secperday * FDI
      fire_TH_diag%avg_HumanIgnition      (kidx0:kidx1)   = &
      fire_TH_diag%avg_HumanIgnition      (kidx0:kidx1)   + secperday * human_i
      fire_TH_diag%avg_LightningIgnition  (kidx0:kidx1)   = &
      fire_TH_diag%avg_LightningIgnition  (kidx0:kidx1)   + secperday * lightning_i
      fire_TH_diag%avg_ROS                (kidx0:kidx1)   = &
      fire_TH_diag%avg_ROS                (kidx0:kidx1)   + secperday * ROS
      fire_TH_diag%avg_FireDuration       (kidx0:kidx1)   = &
      fire_TH_diag%avg_FireDuration       (kidx0:kidx1)   + secperday * FireDuration 
      fire_TH_diag%avg_numfire            (kidx0:kidx1)   = &
      fire_TH_diag%avg_numfire            (kidx0:kidx1)   + secperday * numfire / 10000._dp
      fire_TH_diag%avg_fire_intensity     (kidx0:kidx1)   = &
      fire_TH_diag%avg_fire_intensity     (kidx0:kidx1)   + secperday * I_surface
      fire_TH_diag%avg_postfire_mortality (kidx0:kidx1,:) = &
      fire_TH_diag%avg_postfire_mortality (kidx0:kidx1,:) + secperday * state_FTH%postfire_mortality(kidx0:kidx1,:)
      fire_TH_diag%avg_vegetation_height  (kidx0:kidx1,:) = &
      fire_TH_diag%avg_vegetation_height  (kidx0:kidx1,:) + secperday * vegetation_height
    ENDIF
    IF (fire_TH_options%lcalc_frp) THEN
      fire_TH_diag%avg_FRP_gridcell     (kidx0:kidx1) = &
      fire_TH_diag%avg_FRP_gridcell     (kidx0:kidx1) + secperday * FRP_gridcell
      fire_TH_diag%avg_numfires_gridcell(kidx0:kidx1) = &
      fire_TH_diag%avg_numfires_gridcell(kidx0:kidx1) + secperday * numfires_gridcell
    ENDIF

  END SUBROUTINE burned_frac_thonicke

  SUBROUTINE fire_danger_index(nidx,fuel_1hr_total,fuel_10hr_total,fuel_100hr_total,leaf_moisture,&
                               char_moistfactor,Fuel1To100hrSum,  &
                               ratio_live_fuel,ratio_dead_fuel,SurfArea2Vol,temp_max,temp_min,precip,&
                               RelFuelMoisture,FDI,moisture_scaling,NI_acc)
  USE mo_jsbach_constants, ONLY: secperday, tmelt
  USE mo_land_surface,     ONLY: fract_small
  INTEGER , INTENT(in)    :: nidx
  REAL(dp), INTENT(in)    :: fuel_1hr_total(:), fuel_10hr_total(:)
  REAL(dp), INTENT(in)    :: fuel_100hr_total(:), leaf_moisture(:)
  REAL(dp), INTENT(in)    :: ratio_live_fuel(:), ratio_dead_fuel(:)
  REAL(dp), INTENT(in)    :: Fuel1To100hrSum(:),SurfArea2Vol(3)
  REAL(dp), INTENT(in)    :: temp_max(:), temp_min(:),precip(:)   ! minimum, maximum temperature, precipitation
  REAL(dp), INTENT(in)    :: char_moistfactor(:)
  REAL(dp), INTENT(in)    :: moisture_scaling
  REAL(dp), INTENT(inout) :: NI_acc(:)
  REAL(dp), INTENT(out)   :: FDI(:)
  REAL(dp), INTENT(out)   :: RelFuelMoisture(:)

  REAL(dp) :: d_NI(nidx),alpha_fuel(nidx),char_alpha_fuel(nidx)
  REAL(dp) :: alpha_livegrass(nidx),precip_mm(nidx)
  REAL(dp) :: alpha_fuel_1hr, alpha_fuel_10hr, alpha_fuel_100hr 
  LOGICAL  :: only10 = .false.

    ! convert precipitation from kg/m2s to mm/day and calc moisture scaling
    precip_mm(:)     = precip(:)   * secperday
    alpha_fuel_1hr   = SurfArea2Vol(1) / moisture_scaling
    alpha_fuel_10hr  = SurfArea2Vol(2) / moisture_scaling
    alpha_fuel_100hr = SurfArea2Vol(3) / moisture_scaling

    ! compute Nesterov index (eq. 5 in TH2010)
    WHERE ((precip_mm(:) < 3._dp) .AND.((temp_min(:)-4._dp-tmelt) >= 0._dp)) 
      d_NI(:) = (temp_max(:) - tmelt) * ((temp_max(:) - tmelt) - (temp_min(:) - tmelt - 4._dp))
    ELSEWHERE
      d_NI(:) = 0._dp
    ENDWHERE

    ! accumulates Nesterov index (d_NI) using the global variable NI_acc for cells without rain
    WHERE (d_NI(:) > fract_small) 
      NI_acc(:) = NI_acc(:) + d_NI   
    ELSEWHERE
      NI_acc(:) = 0._dp
    ENDWHERE

    ! sum over fuel classes in eq. 6 TH2010 + influence of livegrass
    IF (only10) THEN
      WHERE (Fuel1To100hrSum(:) > fract_small)
        alpha_fuel(:) = (   alpha_fuel_1hr  * fuel_1hr_total (:) &
                          + alpha_fuel_10hr * fuel_10hr_total(:) &
                        ) / (fuel_1hr_total + fuel_10hr_total(:))
      ELSEWHERE
        alpha_fuel(:) = 0._dp
      ENDWHERE
    ELSE
      WHERE (Fuel1To100hrSum(:) > fract_small)    
        alpha_fuel(:) = (  alpha_fuel_1hr   * fuel_1hr_total  (:) &
                         + alpha_fuel_10hr  * fuel_10hr_total (:) &
                         + alpha_fuel_100hr * fuel_100hr_total(:) &
                        ) / Fuel1To100hrSum(:)
      ELSEWHERE
        alpha_fuel(:) = 0._dp
      ENDWHERE
    ENDIF 

    WHERE ((NI_acc(:) > fract_small) .AND. (leaf_moisture(:) > fract_small)) 
      alpha_livegrass(:) = -LOG(leaf_moisture(:)) / NI_acc(:)
    ELSEWHERE
      alpha_livegrass(:) = 0.0_dp
    ENDWHERE
    char_alpha_fuel(:) = alpha_fuel(:) * ratio_dead_fuel(:) + alpha_livegrass(:) * ratio_live_fuel(:)

    ! weighted average of the relative moisture content eq. 6 in TH2010
    RelFuelMoisture(:) = EXP(-char_alpha_fuel(:) * NI_acc(:))  ! used for ROS calcs. 1000hr therefore ignored for alpha.

    ! calculate fire danger index d_fdi eq. 8 in TH2010
    WHERE ((d_NI(:) <= fract_small) .OR. (char_moistfactor(:) == 0._dp) .OR. (Fuel1To100hrSum(:) <= fract_small)) 
       FDI(:) = 0.0_dp
    ELSEWHERE
       FDI(:) = MAX(0.0_dp,(1.0_dp-((1.0_dp/char_moistfactor(:)) * (EXP(-char_alpha_fuel(:)*NI_acc(:))))))
    ENDWHERE

  END SUBROUTINE fire_danger_index

  SUBROUTINE rate_of_spread(nidx,WindSpeed_ft_min,SurfArea2Vol,RelFuelMoisture,&
                            char_dens_fuel_ave,char_net_fuel,char_moistfactor, &
                            ROS,gamma,IR,lwind_speed_limit,wind_limit,wind_max,wind_slope)
  USE mo_land_surface,     ONLY: fract_small
  INTEGER , INTENT(in)    :: nidx
  REAL(dp), INTENT(in)    :: WindSpeed_ft_min(:)
  REAL(dp), INTENT(in)    :: SurfArea2Vol(:), RelFuelMoisture(:)
  REAL(dp), INTENT(in)    :: char_net_fuel(:)
  REAL(dp), INTENT(in)    :: char_moistfactor(:)
  REAL(dp), INTENT(in)    :: char_dens_fuel_ave(:)
  REAL(dp), INTENT(inout) :: gamma(:), ROS(:)      ! ROS:rate of spread, gamma: used for FireResidenceTime for mortality
  REAL(dp), INTENT(inout) :: IR(:)                 ! reaction intensity
  LOGICAL,  INTENT(in)    :: lwind_speed_limit
  REAL(dp), INTENT(in)    :: wind_limit
  REAL(dp), INTENT(in)    :: wind_max, wind_slope

! locals
  REAL(dp) :: beta(nidx), beta_op(nidx), bet(nidx)
  REAL(dp) :: q_ig(nidx)               ! heat of pre-ignition
  REAL(dp) :: eps(nidx)                ! effective heating number
  REAL(dp) :: b(nidx), c(nidx), e(nidx), a(nidx)
  REAL(dp) :: phi_wind(nidx)
  REAL(dp) :: xi(nidx)
  REAL(dp) :: dummy(nidx), gamma_max(nidx),gamma_aptr(nidx)
  REAL(dp) :: mw_weight(nidx), moist_damp(nidx)

    ! moisture dampening coefficient, see tab. A1
    WHERE (char_moistfactor(:) > 0._dp)
      mw_weight(:) = RelFuelMoisture(:) / char_moistfactor(:)
    ELSEWHERE
      mw_weight(:) = 0._dp
    ENDWHERE

    WHERE ((SurfArea2Vol(:) > 0.0001_dp) .AND. (char_dens_fuel_ave(:) > fract_small))

      ! influence of wind speed (eq. A5-A9)
      beta   (:) = char_dens_fuel_ave(:) / part_dens
      beta_op(:) = 0.200395_dp * (SurfArea2Vol(:)**(-0.8189_dp))
      bet    (:) = beta(:) / beta_op(:)
      
      b(:) = 0.15988_dp * (SurfArea2Vol(:)**0.54_dp)
      c(:) = 7.47_dp    * (EXP(-0.8711_dp  * (SurfArea2Vol(:)**0.55_dp)))
      e(:) = 0.715_dp   * (EXP(-0.01094_dp *  SurfArea2Vol(:)))

      ! heat of pre-ignition eq. A4 in  TH2010
      q_ig(:) = 581.0_dp + 2594.0_dp * RelFuelMoisture(:)

      ! effective heating number eq. A3 in  TH2010
      eps(:) = EXP(-4.528_dp/SurfArea2Vol(:))
    ELSEWHERE
      ROS  (:) = 0._dp
      gamma(:) = 0._dp
      bet  (:) = 0._dp
      q_ig (:) = 0._dp
      eps  (:) = 0._dp
    ENDWHERE

    WHERE (bet(:) > 0._dp)
      phi_wind(:) = c(:) * (WindSpeed_ft_min(:)**b(:)) * (bet(:)**(-e(:)))
    ELSEWHERE
      phi_wind(:) = 0._dp
    ENDWHERE

    WHERE ((SurfArea2Vol(:) > 0.0001_dp) .AND. (char_dens_fuel_ave(:) > fract_small))

      ! propagating flux ratio xi Eq. A2 in  TH2010
      xi(:) = EXP((0.792_dp + 3.7597_dp * (SurfArea2Vol(:)**0.5_dp)) * (beta(:) + 0.1_dp)) / &
                  (192._dp + 7.9095_dp * SurfArea2Vol(:))
 
      ! optimum reaction velocity table A1 in TH2010
      a         (:) = 8.9033_dp * (SurfArea2Vol(:)**(-0.7913_dp))
      dummy     (:) = EXP(a(:) * (1._dp - bet(:)))
      gamma_max (:) = 1._dp / (0.0591_dp + 2.926_dp * (SurfArea2Vol(:)**(-1.5_dp))) ! maximum reaction velocity
      gamma_aptr(:) = gamma_max(:) * (bet(:)**a(:)) * dummy(:)               ! optimum reaction velocity

      moist_damp(:) = MAX(0._dp,(  1._dp - (2.59_dp * mw_weight(:)) &
                                 + (5.11_dp * (mw_weight(:)**2._dp)) - (3.52_dp * (mw_weight(:)**3._dp))))

      ! reaction intensity IR eq. A1 in TH2010
      IR(:) = gamma_aptr(:) * char_net_fuel(:) * HeatContentFuel * moist_damp(:) * MINER_DAMP

      ! for use in postfire mortality (computation of FireResidenceTime)
      gamma(:) = gamma_aptr(:) * moist_damp(:) * MINER_DAMP * char_net_fuel(:)

    ENDWHERE

    ! wind speed limitation
    IF (lwind_speed_limit) THEN
      WHERE ((ir(:) > 0._dp) .AND. (WindSpeed_ft_min(:) > wind_max))
      ! if (WindSpeed_ft_min(i)**3.0_dp/(ir(i)*0.08_dp) .gt. wind_limit) then
      !     limited_wind_speed=wind_limit*0.08*ir(i)**(1.0_dp/3.0_dp)
      !    phi_wind(i)=c(i)*(limited_wind_speed**b(i))*(bet(i)**(-e(i)))
        phi_wind(:) = c(:)*(MAX(wind_max*(1._dp + wind_slope) - wind_slope*WindSpeed_ft_min(:),0._dp)**b(:))*(bet(:)**(-e(:)))
      ENDWHERE
    ENDIF

    ! rate of spread eq. 9in TH2010
    WHERE ((char_dens_fuel_ave(:) <= 0._dp) .OR. (eps(:) <= 0._dp) .OR. (q_ig(:) <= 0._dp))
      ROS(:)=0._dp
    ELSEWHERE
      ROS(:)=(IR(:) * xi(:) * (1._dp + phi_wind(:))) / (char_dens_fuel_ave(:) * eps(:) * q_ig(:))
    ENDWHERE

  END SUBROUTINE rate_of_spread

  SUBROUTINE fuel_consumption(nidx, ntiles, used_pft_array, woody_pft, NI_acc, &
                              burned_frac_diag, RelFuelMoisture, fuel_1hr, fuel_10hr, &
                              fuel_100hr, fuel_1000hr, fuel_1hr_total, &
                              leaf_moisture, livegrass, char_moistfactor, &
                              com_fac, fuel_consum, fuel_consum_tiled, fuel_consum1000perAreaBurned, & 
                              com_fac_1hr, com_fac_10hr, com_fac_100hr, com_fac_1000hr, com_fac_lg)
!  USE mo_cbal_cpools, ONLY: frac_green_aboveGround
  INTEGER,  INTENT(in)    :: nidx, ntiles
  LOGICAL,  INTENT(in)    :: used_pft_array(:,:)
  LOGICAL,  INTENT(in)    :: woody_pft(:,:)
  REAL(dp), INTENT(in)    :: NI_acc(:)
  REAL(dp), INTENT(in)    :: burned_frac_diag(:,:)              ! fraction burned per pft, grid cell
  REAL(dp), INTENT(in)    :: RelFuelMoisture(:)
  REAL(dp), INTENT(in)    :: fuel_1hr(:,:), fuel_10hr(:,:), fuel_100hr(:,:), fuel_1000hr(:,:)
  REAL(dp), INTENT(in)    :: fuel_1hr_total(:)
  REAL(dp), INTENT(in)    :: leaf_moisture(:)
  REAL(dp), INTENT(in)    :: livegrass(:)
  REAL(dp), INTENT(in)    :: char_moistfactor(:) ! extinction moisture, weighted with relative amount of fuel per tile
  REAL(dp), INTENT(inout) :: com_fac(:,:)
  REAL(dp), INTENT(inout) :: fuel_consum(:)                ! fuel consumption per vegetated area
  REAL(dp), INTENT(inout) :: fuel_consum1000perAreaBurned(:)              ! fuel consumption per area burned
  REAL(dp), INTENT(inout) :: fuel_consum_tiled(:,:)
  REAL(dp), INTENT(out)   :: com_fac_lg(:,:), com_fac_1hr(:,:),com_fac_10hr(:,:)
  REAL(dp), INTENT(out)   :: com_fac_100hr(:,:), com_fac_1000hr(:,:)             ! com_fac : combustion factor of different fuels

! local variables:
  REAL(dp) :: fc_1hr(nidx,ntiles), fc_10hr(nidx,ntiles) 
  REAL(dp) :: fc_100hr(nidx,ntiles), fc_1000hr(nidx,ntiles) 
  REAL(dp) :: moist_lg_d1hr(nidx), moist_1hr(nidx), moist_10_100hr(nidx)
  REAL(dp) :: RelFuelMoisture_1hr(nidx)
  INTEGER  :: i, itiles

    ! initialize variables
    fuel_consum1000perAreaBurned(:)  = 0._dp
    fuel_consum   (:)   = 0._dp
    com_fac       (:,:) = 0._dp
    com_fac_lg    (:,:) = 0._dp
    com_fac_1hr   (:,:) = 0._dp
    com_fac_10hr  (:,:) = 0._dp
    com_fac_100hr (:,:) = 0._dp
    com_fac_1000hr(:,:) = 0._dp
    fc_1hr        (:,:) = 0._dp
    fc_10hr       (:,:) = 0._dp
    fc_100hr      (:,:) = 0._dp
    fc_1000hr     (:,:) = 0._dp

    RelFuelMoisture_1hr   (:) = EXP(-alpha_1hr * NI_acc(:))
    WHERE(fuel_1hr_total(:) > 0._dp)
! GL 3/2015 modification of moisture averaging, now different from SPITFIRE documentation but more reasonable 
      moist_lg_d1hr(:) = (RelFuelMoisture_1hr * fuel_1hr_total + leaf_moisture * livegrass) &
                       / (fuel_1hr_total + livegrass)
    ELSEWHERE
      moist_lg_d1hr(:) = 0._dp
    ENDWHERE
    WHERE(char_moistfactor(:) > 0._dp)
      moist_1hr     (:) = MIN(moist_lg_d1hr  (:) / char_moistfactor(:),1.0_dp)
      moist_10_100hr(:) = MIN(RelFuelMoisture(:) / char_moistfactor(:),1.0_dp)
    ELSEWHERE
      moist_1hr     (:) = 1._dp
      moist_10_100hr(:) = 1._dp
    ENDWHERE

    ! Combustion factors 
    DO itiles=1,ntiles
      DO i=1,nidx
        IF (used_pft_array(i,itiles)) THEN
          IF (.NOT. woody_pft(i,itiles)) THEN
            ! green pool for grasses, potentially combusted --> livegrass in spitfire
            IF (leaf_moisture(i) <= 0.18_dp) THEN
              com_fac_lg(i,itiles) = 1._dp * (1.0_dp - MINER_TOT)
            ELSEIF ((leaf_moisture(i) > 0.18_dp) .AND. (leaf_moisture(i) <= 0.73_dp)) THEN
              com_fac_lg(i,itiles) = (1.11_dp -  0.62_dp * leaf_moisture(i)) * (1._dp - MINER_TOT)
            ELSEIF ((leaf_moisture(i) > 0.73_dp) .AND. (leaf_moisture(i) <= 1._dp)) THEN
              com_fac_lg(i,itiles) = (2.45_dp - 2.45_dp * leaf_moisture(i)) * (1._dp - MINER_TOT)
            ENDIF
          ENDIF ! not woody
          IF (moist_1hr(i) <= 0.18_dp) THEN
            com_fac_1hr(i,itiles) = 1.0_dp * (1.0_dp - MINER_TOT)
          ELSEIF ((moist_1hr(i) > 0.18_dp) .AND. (moist_1hr(i) <= 0.73_dp)) THEN
            com_fac_1hr(i,itiles) = (1.11_dp - 0.62_dp * moist_1hr(i)) * (1.0_dp - MINER_TOT)
          ELSEIF ((moist_1hr(i) > 0.73_dp) .AND. (moist_1hr(i) <= 1._dp)) THEN
            com_fac_1hr(i,itiles) = (2.45_dp - 2.45_dp * moist_1hr(i)) * (1.0_dp - MINER_TOT)
          ENDIF
          ! 10hr fuel consumption
          IF (moist_10_100hr(i) <= 0.13_dp) THEN
            com_fac_10hr(i,itiles) = 1._dp * (1._dp - MINER_TOT)
          ELSEIF ((moist_10_100hr(i) > 0.13_dp) .AND. (moist_10_100hr(i) <= 0.51_dp)) THEN
            com_fac_10hr(i,itiles) = (1.09_dp - 0.72_dp * moist_10_100hr(i)) * (1._dp - MINER_TOT)
          ELSEIF ((moist_10_100hr(i) > 0.51_dp) .AND. (moist_10_100hr(i) <= 1._dp)) THEN
            com_fac_10hr(i,itiles) = (1.47_dp - 1.47_dp * moist_10_100hr(i)) * (1._dp - MINER_TOT)
          ENDIF
          ! 100hr fuel consumption
          com_fac_100hr(i,itiles) = MIN(0.45_dp,(0.98_dp - 0.85_dp * moist_10_100hr(i))) * (1._dp - MINER_TOT)
          ! 1000hr fuel consumption
          com_fac_1000hr(i,itiles) = MIN(0.45_dp,(0.8_dp - 0.8_dp  * moist_10_100hr(i))) * (1._dp - MINER_TOT)
        ENDIF
      ENDDO
    ENDDO

    ! Tiled fuel consumption
    WHERE (used_pft_array(:,:))
      fc_1hr   (:,:) = com_fac_1hr   (:,:) * fuel_1hr   (:,:) 
      fc_10hr  (:,:) = com_fac_10hr  (:,:) * fuel_10hr  (:,:) 
      fc_100hr (:,:) = com_fac_100hr (:,:) * fuel_100hr (:,:)
      fc_1000hr(:,:) = com_fac_1000hr(:,:) * fuel_1000hr(:,:) 
    ENDWHERE

    ! Aggregation
    com_fac(:,:) = (com_fac_1hr(:,:) + com_fac_10hr(:,:) + com_fac_100hr(:,:)) / 3._dp
    fuel_consum_tiled(:,:) = fc_1hr(:,:) + fc_10hr(:,:) + fc_100hr(:,:) 
    fuel_consum1000perAreaBurned(:) = sum(fc_1hr + fc_10hr + fc_100hr + fc_1000hr, dim=2)
    fuel_consum(:) = sum((fc_1hr + fc_10hr + fc_100hr)*burned_frac_diag, dim=2)
        
    ! Consistency check
    IF (ANY((fc_1hr(:,:) > fuel_1hr(:,:)) .OR. (fc_10hr(:,:) > fuel_10hr(:,:)) .OR. (fc_100hr(:,:) > fuel_100hr(:,:)))) THEN
      DO itiles=1,ntiles
        DO i=1,nidx
          IF (used_pft_array(i,itiles)) THEN
            IF ((fc_1hr  (i,itiles) > fuel_1hr  (i,itiles)) .OR. &
                (fc_10hr (i,itiles) > fuel_10hr (i,itiles)) .OR. &
                (fc_100hr(i,itiles) > fuel_100hr(i,itiles))) THEN
              WRITE(message_text,*) 'Fuel consumed higher than available fuel for block point (',i,itiles,')'
              CALL message('fuel_consumption',TRIM(message_text))
              WRITE(message_text,*) 'FCs. ' , fc_1hr(i,itiles), fc_10hr(i,itiles), fc_100hr(i,itiles)
              CALL message('fuel_consumption',TRIM(message_text))
              WRITE(message_text,*) 'Fuels:', fuel_1hr(i,itiles),fuel_10hr(i,itiles),fuel_100hr(i,itiles)
              CALL message('fuel_consumption',TRIM(message_text))
              WRITE(message_text,*) 'com_facs:', com_fac_1hr(i,itiles),moist_1hr(i),com_fac_10hr(i,itiles),com_fac_100hr(i,itiles)
              CALL message('fuel_consumption',TRIM(message_text))

            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE fuel_consumption

! FUNCTION HUMAN IGNITION
! Calculation of the numbers of ignition caused by humans.
  FUNCTION human_ign(nidx,popden,a_nd,param)
  INTEGER,  INTENT(in) :: nidx
  REAL(dp), INTENT(in) :: param
  REAL(dp), INTENT(in) :: popden(:)
  REAL(dp), INTENT(in) :: a_nd(:)
  REAL(dp)             :: human_ign(nidx)

    ! using equation from c code: parameter=6.8, now in namelist default= 5
    human_ign(:) = param * (EXP(-0.5_dp * (popden(:)**0.5_dp))) * a_nd(:) / 100._dp * popden(:)
    
  END FUNCTION human_ign

! transformation of combusted fuels to carbon pools
  SUBROUTINE pool_consumption(nidx,ntiles,woody_pft,used_pft_array,burned_frac_diag,                             &
                              com_fac_1hr,com_fac_10hr,com_fac_100hr,com_fac_1000hr,com_fac_lg,           &
                              Cpool_litter_wood_ag, Cpool_litter_green_ag,Cpool_green,                    &
                              frac_1hr_wood,frac_10hr_wood,frac_100hr_wood,frac_1000hr_wood,  &
                              burned_Cpool_litter_wood_ag,burned_Cpool_litter_green_ag,burned_Cpool_green)
  USE mo_land_surface, ONLY: fract_small
  USE mo_cbal_cpools,  ONLY: frac_green_aboveGround
  INTEGER,  INTENT(in)    :: nidx, ntiles
  REAL(dp), INTENT(in)    :: com_fac_1hr(:,:),com_fac_10hr(:,:),com_fac_100hr(:,:),com_fac_1000hr(:,:),com_fac_lg(:,:)
  LOGICAL,  INTENT(in)    :: woody_pft(:,:),used_pft_array(:,:)
  REAL(dp), INTENT(in)    :: burned_frac_diag(:,:)
  REAL(dp), INTENT(in)    :: Cpool_green(:,:)
  REAL(dp), INTENT(in)    :: Cpool_litter_wood_ag(:,:), Cpool_litter_green_ag(:,:)
  REAL(dp), INTENT(inout) :: frac_1hr_wood(:,:),frac_10hr_wood(:,:)
  REAL(dp), INTENT(inout) :: frac_100hr_wood(:,:),frac_1000hr_wood(:,:)
  REAL(dp), INTENT(inout) :: burned_Cpool_litter_wood_ag(:,:), burned_Cpool_litter_green_ag(:,:), burned_Cpool_green(:,:)

! local variables 
  REAL(dp)    :: burned_litter_green_leaves(nidx, ntiles)
  REAL(dp)    :: burned_litter_green_1hr(nidx, ntiles) 
  REAL(dp)    :: burned_litter_wood_1hr(nidx,ntiles),burned_litter_wood_10hr(nidx,ntiles)
  REAL(dp)    :: burned_litter_wood_100hr(nidx,ntiles),burned_litter_wood_1000hr(nidx,ntiles)
  INTEGER(dp) :: i, itile

    ! initialize burned pools
    burned_Cpool_green          (:,:) = 0._dp
    burned_litter_green_1hr     (:,:) = 0._dp
    burned_litter_green_leaves  (:,:) = 0._dp
    burned_litter_wood_1hr      (:,:) = 0._dp
    burned_litter_wood_10hr     (:,:) = 0._dp
    burned_litter_wood_100hr    (:,:) = 0._dp
    burned_litter_wood_1000hr   (:,:) = 0._dp

    ! compute combusted litter in pools
    WHERE (woody_pft(:,:))
      burned_litter_green_leaves= Cpool_litter_green_ag(:,:)* frac_green_active  * com_fac_1hr(:,:) * burned_frac_diag(:,:)
      burned_litter_wood_1hr    = Cpool_litter_wood_ag(:,:) * frac_1hr_wood(:,:) * com_fac_1hr(:,:) * burned_frac_diag(:,:)
      burned_litter_wood_10hr   = Cpool_litter_wood_ag(:,:) * frac_10hr_wood(:,:)* com_fac_10hr(:,:) * burned_frac_diag(:,:)
      burned_litter_wood_100hr  = Cpool_litter_wood_ag(:,:) * frac_100hr_wood(:,:)*com_fac_100hr (:,:)*burned_frac_diag(:,:)
      burned_litter_wood_1000hr = Cpool_litter_wood_ag(:,:) * frac_1000hr_wood(:,:)*com_fac_1000hr(:,:)*burned_frac_diag(:,:)
    ELSEWHERE (used_pft_array(:,:))
      ! combustion of live grass
      burned_Cpool_green(:,:) = frac_green_aboveGround * Cpool_green(:,:) * com_fac_lg (:,:) * burned_frac_diag(:,:)
      burned_litter_green_1hr(:,:) =           Cpool_litter_green_ag(:,:) * com_fac_1hr(:,:) * burned_frac_diag(:,:)
    ENDWHERE

    ! sum fuels for burned litter
    burned_Cpool_litter_wood_ag (:,:) = burned_litter_wood_1hr    (:,:) + burned_litter_wood_10hr   (:,:) &
                                      + burned_litter_wood_100hr  (:,:) + burned_litter_wood_1000hr (:,:)

    burned_Cpool_litter_green_ag(:,:) = burned_litter_green_leaves(:,:) + burned_litter_green_1hr(:,:)

    ! Consistency check
    DO itile=1,ntiles
      DO i=1,nidx
        IF (burned_Cpool_litter_green_ag(i,itile) > Cpool_litter_green_ag(i,itile)) THEN
          WRITE(message_text,*) 'pool_consumption, burned_Cpool_litter_green larger than Cpool at local point (',i,itile,')'
          CALL message('pool_consumption',TRIM(message_text))
          WRITE(message_text,*) burned_Cpool_litter_green_ag(i,itile), Cpool_litter_green_ag(i,itile)
          CALL message('pool_consumption',TRIM(message_text))
          WRITE(message_text,*) burned_litter_green_leaves(i,itile), burned_litter_green_1hr(i,itile),com_fac_1hr(i,itile) 
          CALL message('pool_consumption',TRIM(message_text))
          WRITE(message_text,*) frac_1hr_wood(i,itile),burned_frac_diag(i,itile)
          CALL message('pool_consumption',TRIM(message_text))
        ENDIF
      ENDDO
    ENDDO

    ! --- update fractions
    WHERE (woody_pft(:,:) .AND.  ((Cpool_litter_wood_ag(:,:) - burned_Cpool_litter_wood_ag(:,:)) > fract_small))
      frac_1hr_wood   (:,:) = (Cpool_litter_wood_ag(:,:) * frac_1hr_wood   (:,:) - burned_litter_wood_1hr   (:,:)) &
                            / (Cpool_litter_wood_ag(:,:) - burned_Cpool_litter_wood_ag(:,:))
      frac_10hr_wood  (:,:) = (Cpool_litter_wood_ag(:,:) * frac_10hr_wood  (:,:) - burned_litter_wood_10hr  (:,:)) &
                            / (Cpool_litter_wood_ag(:,:) - burned_Cpool_litter_wood_ag(:,:))
      frac_100hr_wood (:,:) = (Cpool_litter_wood_ag(:,:) * frac_100hr_wood (:,:) - burned_litter_wood_100hr (:,:)) &
                            / (Cpool_litter_wood_ag(:,:) - burned_Cpool_litter_wood_ag(:,:))
      frac_1000hr_wood(:,:) = (Cpool_litter_wood_ag(:,:) * frac_1000hr_wood(:,:) - burned_litter_wood_1000hr(:,:)) &
                            / (Cpool_litter_wood_ag(:,:) - burned_Cpool_litter_wood_ag(:,:))
    ELSEWHERE
      frac_1hr_wood   (:,:) = frac_1hr_wood_new
      frac_10hr_wood  (:,:) = frac_10hr_wood_new
      frac_100hr_wood (:,:) = frac_100hr_wood_new
      frac_1000hr_wood(:,:) = frac_1000hr_wood_new
    ENDWHERE

  END SUBROUTINE pool_consumption

  SUBROUTINE pool_consumption_yasso(nidx,ntiles,woody_pft,used_pft_array,burned_frac_diag,                       &
                              com_fac_1hr,com_fac_10hr,com_fac_100hr,com_fac_1000hr,com_fac_lg,           &
                              Cpool_green,                                                                &
                              YCpool_acid_ag1,       YCpool_acid_ag2,       YCpool_water_ag1,             &
                              YCpool_water_ag2,      YCpool_ethanol_ag1,    YCpool_ethanol_ag2,           &
                              YCpool_nonsoluble_ag1, YCpool_nonsoluble_ag2,                               &
                              frac_1hr_wood,frac_10hr_wood,frac_100hr_wood,frac_1000hr_wood,              &
                              burned_acid_ag1, burned_acid_ag2, burned_water_ag1, burned_water_ag2,       &
                              burned_ethanol_ag1, burned_ethanol_ag2, burned_nonsoluble_ag1,              &
                              burned_nonsoluble_ag2, burned_litter_wood, burned_litter_green,             &
                              burned_Cpool_green)
  USE mo_land_surface, ONLY:  fract_small
  USE mo_cbal_cpools,  ONLY: frac_green_aboveGround
  INTEGER,  INTENT(in)    :: nidx, ntiles
  REAL(dp), INTENT(in)    :: com_fac_1hr(:,:),com_fac_10hr(:,:),com_fac_100hr(:,:),com_fac_1000hr(:,:),com_fac_lg(:,:)
  LOGICAL,  INTENT(in)    :: woody_pft(:,:),used_pft_array(:,:)
  REAL(dp), INTENT(in)    :: burned_frac_diag(:,:)
  REAL(dp), INTENT(in)    :: Cpool_green(:,:)
  REAL(dp), INTENT(inout) :: YCpool_acid_ag1      (:,:), YCpool_acid_ag2       (:,:)
  REAL(dp), INTENT(inout) :: YCpool_water_ag1     (:,:), YCpool_water_ag2      (:,:)
  REAL(dp), INTENT(inout) :: YCpool_ethanol_ag1   (:,:), YCpool_ethanol_ag2    (:,:)
  REAL(dp), INTENT(inout) :: YCpool_nonsoluble_ag1(:,:), YCpool_nonsoluble_ag2 (:,:)
  REAL(dp), INTENT(inout) :: frac_1hr_wood(:,:),frac_10hr_wood(:,:)
  REAL(dp), INTENT(inout) :: frac_100hr_wood(:,:),frac_1000hr_wood(:,:)
  REAL(dp), INTENT(inout) :: burned_acid_ag1(:,:), burned_acid_ag2(:,:), burned_water_ag1(:,:), burned_water_ag2(:,:)
  REAL(dp), INTENT(inout) :: burned_ethanol_ag1(:,:),burned_ethanol_ag2(:,:),burned_nonsoluble_ag1(:,:),burned_nonsoluble_ag2(:,:)
  REAL(dp), INTENT(inout) :: burned_litter_green(:,:), burned_litter_wood(:,:), burned_Cpool_green(:,:)

! local variables 
  REAL(dp)    :: burned_acid_leaves          (nidx,ntiles), burned_water_leaves         (nidx,ntiles)
  REAL(dp)    :: burned_ethanol_leaves       (nidx,ntiles), burned_nonsoluble_leaves    (nidx,ntiles)
  REAL(dp)    :: burned_acid_ag1_1hr         (nidx,ntiles), burned_acid_ag2_1hr         (nidx,ntiles)
  REAL(dp)    :: burned_water_ag1_1hr        (nidx,ntiles), burned_water_ag2_1hr        (nidx,ntiles)
  REAL(dp)    :: burned_ethanol_ag1_1hr      (nidx,ntiles), burned_ethanol_ag2_1hr      (nidx,ntiles)
  REAL(dp)    :: burned_nonsoluble_ag1_1hr   (nidx,ntiles), burned_nonsoluble_ag2_1hr   (nidx,ntiles)
  REAL(dp)    :: burned_acid_ag2_10hr        (nidx,ntiles)
  REAL(dp)    :: burned_water_ag2_10hr       (nidx,ntiles)
  REAL(dp)    :: burned_ethanol_ag2_10hr     (nidx,ntiles)
  REAL(dp)    :: burned_nonsoluble_ag2_10hr  (nidx,ntiles)
  REAL(dp)    :: burned_acid_ag2_100hr       (nidx,ntiles)
  REAL(dp)    :: burned_water_ag2_100hr      (nidx,ntiles)
  REAL(dp)    :: burned_ethanol_ag2_100hr    (nidx,ntiles)
  REAL(dp)    :: burned_nonsoluble_ag2_100hr (nidx,ntiles)
  REAL(dp)    :: burned_acid_ag2_1000hr      (nidx,ntiles)
  REAL(dp)    :: burned_water_ag2_1000hr     (nidx,ntiles)
  REAL(dp)    :: burned_ethanol_ag2_1000hr   (nidx,ntiles)
  REAL(dp)    :: burned_nonsoluble_ag2_1000hr(nidx,ntiles)
  REAL(dp)    :: wood_fac_1hr                (nidx,ntiles), wood_fac_10hr               (nidx,ntiles)
  REAL(dp)    :: wood_fac_100hr              (nidx,ntiles), wood_fac_1000hr             (nidx,ntiles)
  REAL(dp)    :: litter_wood                 (nidx,ntiles) !, litter_green                (nidx,ntiles)
  REAL(dp)    :: burned_wood                 (nidx,ntiles) !, burned_green                (nidx,ntiles)
!  REAL(dp)    :: burned_leaves               (nidx,ntiles)
  REAL(dp)    :: burned_wood_1hr             (nidx,ntiles), burned_wood_10hr            (nidx,ntiles)
  REAL(dp)    :: burned_wood_100hr           (nidx,ntiles), burned_wood_1000hr          (nidx,ntiles)

    ! Initializations
    burned_Cpool_green          (:,:) = 0._dp
    burned_acid_ag1_1hr         (:,:) = 0._dp
    burned_water_ag1_1hr        (:,:) = 0._dp
    burned_ethanol_ag1_1hr      (:,:) = 0._dp
    burned_nonsoluble_ag1_1hr   (:,:) = 0._dp
    burned_acid_leaves          (:,:) = 0._dp
    burned_water_leaves         (:,:) = 0._dp
    burned_ethanol_leaves       (:,:) = 0._dp
    burned_nonsoluble_leaves    (:,:) = 0._dp
    burned_acid_ag2_1hr         (:,:) = 0._dp
    burned_water_ag2_1hr        (:,:) = 0._dp
    burned_ethanol_ag2_1hr      (:,:) = 0._dp
    burned_nonsoluble_ag2_1hr   (:,:) = 0._dp
    burned_acid_ag2_10hr        (:,:) = 0._dp
    burned_water_ag2_10hr       (:,:) = 0._dp
    burned_ethanol_ag2_10hr     (:,:) = 0._dp
    burned_nonsoluble_ag2_10hr  (:,:) = 0._dp
    burned_acid_ag2_100hr       (:,:) = 0._dp
    burned_water_ag2_100hr      (:,:) = 0._dp
    burned_ethanol_ag2_100hr    (:,:) = 0._dp
    burned_nonsoluble_ag2_100hr (:,:) = 0._dp
    burned_acid_ag2_1000hr      (:,:) = 0._dp
    burned_water_ag2_1000hr     (:,:) = 0._dp
    burned_ethanol_ag2_1000hr   (:,:) = 0._dp
    burned_nonsoluble_ag2_1000hr(:,:) = 0._dp

    wood_fac_1hr   (:,:) = com_fac_1hr   (:,:) * frac_1hr_wood   (:,:) * burned_frac_diag(:,:)
    wood_fac_10hr  (:,:) = com_fac_10hr  (:,:) * frac_10hr_wood  (:,:) * burned_frac_diag(:,:)
    wood_fac_100hr (:,:) = com_fac_100hr (:,:) * frac_100hr_wood (:,:) * burned_frac_diag(:,:)
    wood_fac_1000hr(:,:) = com_fac_1000hr(:,:) * frac_1000hr_wood(:,:) * burned_frac_diag(:,:)

    ! compute combusted litter in pools
    WHERE (woody_pft(:,:))
      burned_acid_leaves          (:,:) = YCpool_acid_ag1      (:,:) * com_fac_1hr    (:,:) * burned_frac_diag(:,:)
      burned_water_leaves         (:,:) = YCpool_water_ag1     (:,:) * com_fac_1hr    (:,:) * burned_frac_diag(:,:)
      burned_ethanol_leaves       (:,:) = YCpool_ethanol_ag1   (:,:) * com_fac_1hr    (:,:) * burned_frac_diag(:,:)
      burned_nonsoluble_leaves    (:,:) = YCpool_nonsoluble_ag1(:,:) * com_fac_1hr    (:,:) * burned_frac_diag(:,:)
      burned_acid_ag2_1hr         (:,:) = YCpool_acid_ag2      (:,:) * wood_fac_1hr   (:,:)
      burned_water_ag2_1hr        (:,:) = YCpool_water_ag2     (:,:) * wood_fac_1hr   (:,:)
      burned_ethanol_ag2_1hr      (:,:) = YCpool_ethanol_ag2   (:,:) * wood_fac_1hr   (:,:)
      burned_nonsoluble_ag2_1hr   (:,:) = YCpool_nonsoluble_ag2(:,:) * wood_fac_1hr   (:,:)
      burned_acid_ag2_10hr        (:,:) = YCpool_acid_ag2      (:,:) * wood_fac_10hr  (:,:)
      burned_water_ag2_10hr       (:,:) = YCpool_water_ag2     (:,:) * wood_fac_10hr  (:,:)
      burned_ethanol_ag2_10hr     (:,:) = YCpool_ethanol_ag2   (:,:) * wood_fac_10hr  (:,:)
      burned_nonsoluble_ag2_10hr  (:,:) = YCpool_nonsoluble_ag2(:,:) * wood_fac_10hr  (:,:)
      burned_acid_ag2_100hr       (:,:) = YCpool_acid_ag2      (:,:) * wood_fac_100hr (:,:)
      burned_water_ag2_100hr      (:,:) = YCpool_water_ag2     (:,:) * wood_fac_100hr (:,:)
      burned_ethanol_ag2_100hr    (:,:) = YCpool_ethanol_ag2   (:,:) * wood_fac_100hr (:,:)
      burned_nonsoluble_ag2_100hr (:,:) = YCpool_nonsoluble_ag2(:,:) * wood_fac_100hr (:,:)
      burned_acid_ag2_1000hr      (:,:) = YCpool_acid_ag2      (:,:) * wood_fac_1000hr(:,:)
      burned_water_ag2_1000hr     (:,:) = YCpool_water_ag2     (:,:) * wood_fac_1000hr(:,:)
      burned_ethanol_ag2_1000hr   (:,:) = YCpool_ethanol_ag2   (:,:) * wood_fac_1000hr(:,:)
      burned_nonsoluble_ag2_1000hr(:,:) = YCpool_nonsoluble_ag2(:,:) * wood_fac_1000hr(:,:)
    ELSEWHERE (used_pft_array(:,:))
      ! combustion of live grass
      burned_Cpool_green          (:,:) = frac_green_aboveGround * Cpool_green(:,:) * com_fac_lg (:,:) * burned_frac_diag(:,:)
      burned_acid_ag1_1hr         (:,:) = YCpool_acid_ag1                     (:,:) * com_fac_1hr(:,:) * burned_frac_diag(:,:)
      burned_water_ag1_1hr        (:,:) = YCpool_water_ag1                    (:,:) * com_fac_1hr(:,:) * burned_frac_diag(:,:)
      burned_ethanol_ag1_1hr      (:,:) = YCpool_ethanol_ag1                  (:,:) * com_fac_1hr(:,:) * burned_frac_diag(:,:)
      burned_nonsoluble_ag1_1hr   (:,:) = YCpool_nonsoluble_ag1               (:,:) * com_fac_1hr(:,:) * burned_frac_diag(:,:)
    ENDWHERE

    ! sum fuels for litter according to yasso classes
    ! ... burned green
    burned_acid_ag1      (:,:) = burned_acid_ag1_1hr        (:,:) + burned_acid_leaves          (:,:)
    burned_water_ag1     (:,:) = burned_water_ag1_1hr       (:,:) + burned_water_leaves         (:,:)
    burned_ethanol_ag1   (:,:) = burned_ethanol_ag1_1hr     (:,:) + burned_ethanol_leaves       (:,:)
    burned_nonsoluble_ag1(:,:) = burned_nonsoluble_ag1_1hr  (:,:) + burned_nonsoluble_leaves    (:,:)

    burned_litter_green  (:,:) = burned_acid_ag1   (:,:) + burned_water_ag1     (:,:) &
                               + burned_ethanol_ag1(:,:) + burned_nonsoluble_ag1(:,:)

    ! ... burned woody
    burned_acid_ag2      (:,:) = burned_acid_ag2_1hr        (:,:) + burned_acid_ag2_10hr        (:,:) & 
                               + burned_acid_ag2_100hr      (:,:) + burned_acid_ag2_1000hr      (:,:)
    burned_water_ag2     (:,:) = burned_water_ag2_1hr       (:,:) + burned_water_ag2_10hr       (:,:) & 
                               + burned_water_ag2_100hr     (:,:) + burned_water_ag2_1000hr     (:,:)
    burned_ethanol_ag2   (:,:) = burned_ethanol_ag2_1hr     (:,:) + burned_ethanol_ag2_10hr     (:,:) & 
                               + burned_ethanol_ag2_100hr   (:,:) + burned_ethanol_ag2_1000hr   (:,:)
    burned_nonsoluble_ag2(:,:) = burned_nonsoluble_ag2_1hr  (:,:) + burned_nonsoluble_ag2_10hr  (:,:) & 
                               + burned_nonsoluble_ag2_100hr(:,:) + burned_nonsoluble_ag2_1000hr(:,:)

    burned_litter_wood   (:,:) = burned_acid_ag2   (:,:) + burned_water_ag2     (:,:) &
                               + burned_ethanol_ag2(:,:) + burned_nonsoluble_ag2(:,:)

    ! ... and total
!    litter_green         (:,:) = YCpool_acid_ag1            (:,:) + YCpool_water_ag1            (:,:) &
!                               + YCpool_ethanol_ag1         (:,:) + YCpool_nonsoluble_ag1       (:,:)
    litter_wood          (:,:) = YCpool_acid_ag2            (:,:) + YCpool_water_ag2            (:,:) &
                               + YCpool_ethanol_ag2         (:,:) + YCpool_nonsoluble_ag2       (:,:)
!    burned_green         (:,:) = burned_acid_ag1            (:,:) + burned_water_ag1            (:,:) &
!                               + burned_ethanol_ag1         (:,:) + burned_nonsoluble_ag1       (:,:)
    burned_wood          (:,:) = burned_acid_ag2            (:,:) + burned_water_ag2            (:,:) &
                               + burned_ethanol_ag2         (:,:) + burned_nonsoluble_ag2       (:,:)
    burned_wood_1hr      (:,:) = burned_acid_ag2_1hr        (:,:) + burned_water_ag2_1hr        (:,:) &
                               + burned_ethanol_ag2_1hr     (:,:) + burned_nonsoluble_ag2_1hr   (:,:)
    burned_wood_10hr     (:,:) = burned_acid_ag2_10hr       (:,:) + burned_water_ag2_10hr       (:,:) &
                               + burned_ethanol_ag2_10hr    (:,:) + burned_nonsoluble_ag2_10hr  (:,:)
    burned_wood_100hr    (:,:) = burned_acid_ag2_100hr      (:,:) + burned_water_ag2_100hr      (:,:) &
                               + burned_ethanol_ag2_100hr   (:,:) + burned_nonsoluble_ag2_100hr (:,:)
    burned_wood_1000hr   (:,:) = burned_acid_ag2_1000hr     (:,:) + burned_water_ag2_1000hr     (:,:) &
                               + burned_ethanol_ag2_1000hr  (:,:) + burned_nonsoluble_ag2_1000hr(:,:)
!    burned_leaves        (:,:) = burned_acid_leaves         (:,:) + burned_water_leaves         (:,:) &
!                               + burned_ethanol_leaves      (:,:) + burned_nonsoluble_leaves    (:,:)

    WHERE (woody_pft  (:,:) .AND. ((litter_wood(:,:) - burned_wood     (:,:)) > fract_small))
      frac_1hr_wood   (:,:) =      (litter_wood(:,:) * frac_1hr_wood   (:,:) - burned_wood_1hr   (:,:)) &
                            /      (litter_wood(:,:) - burned_wood     (:,:))
      frac_10hr_wood  (:,:) =      (litter_wood(:,:) * frac_10hr_wood  (:,:) - burned_wood_10hr  (:,:)) &
                            /      (litter_wood(:,:) - burned_wood     (:,:))
      frac_100hr_wood (:,:) =      (litter_wood(:,:) * frac_100hr_wood (:,:) - burned_wood_100hr (:,:)) &
                            /      (litter_wood(:,:) - burned_wood     (:,:))
      frac_1000hr_wood(:,:) =      (litter_wood(:,:) * frac_1000hr_wood(:,:) - burned_wood_1000hr(:,:)) &
                            /      (litter_wood(:,:) - burned_wood     (:,:))
    ELSEWHERE
      frac_1hr_wood   (:,:) = frac_1hr_wood_new
      frac_10hr_wood  (:,:) = frac_10hr_wood_new
      frac_100hr_wood (:,:) = frac_100hr_wood_new
      frac_1000hr_wood(:,:) = frac_1000hr_wood_new
    ENDWHERE

  END SUBROUTINE pool_consumption_yasso

  SUBROUTINE fire_mortality(lctlib,nidx,ntiles,woody_pft,cover_type,I_surface,& 
                            com_fac,fuel_consum,gamma,vegetation_height,Lethal2ResTime, &
                            mortality_crownkill,mortality_camb,scorch_height,postfire_mortality)
  TYPE (lctlib_type), INTENT(in)    :: lctlib
  INTEGER,            INTENT(in)    :: nidx, ntiles
  LOGICAL,            INTENT(in)    :: woody_pft(:,:)
  INTEGER,            INTENT(in)    :: cover_type(:,:)
  REAL(dp),           INTENT(in)    :: I_surface(:)             ! fire intensity
  REAL(dp),           INTENT(in)    :: com_fac(:,:)             ! fraction of combusted fuels averaged over 1-100hr fuels
  REAL(dp),           INTENT(in)    :: fuel_consum(:)           
  REAL(dp),           INTENT(in)    :: gamma(:)                 ! gamma:max reaction velocity*mineral dampening*moisture dampening
  REAL(dp),           INTENT(in)    :: vegetation_height(:,:)
  REAL(dp),           INTENT(in)    :: Lethal2ResTime
  REAL(dp),           INTENT(inout) :: mortality_crownkill(:,:)
  REAL(dp),           INTENT(inout) :: mortality_camb(:,:)      ! fraction killed by cambial kill [0,1]
  REAL(dp),           INTENT(inout) :: scorch_height(:,:)       ! [m]
  REAL(dp),           INTENT(inout) :: postfire_mortality(:,:)

  ! local variables
!  REAL(dp) :: crownkill(nidx,ntiles)                ! probability of crownkill
  REAL(dp) :: FireResidenceTime(nidx,ntiles)                    ! residence time of fire
  REAL(dp) :: crownkill_class(5)
  REAL(dp) :: mortality_crownkill_class(5)
  REAL(dp) :: mortality_camb_class(5)
  REAL(dp) :: tau_c,bark_thickness,DBH,stem_length  ! tau_c: critical time of fire residence
  REAL(dp) :: crown_length
  REAL(dp) :: height                                
  INTEGER  :: i,itile,class

    ! initializations
    scorch_height      (:,:) = 0._dp
    FireResidenceTime              (:,:) = 0._dp
!    crownkill          (:,:) = 0._dp
    mortality_crownkill(:,:) = 1._dp
    mortality_camb     (:,:) = 1._dp

    ! scorch height, FireResidenceTime residence time of fire
    DO itile=1,ntiles
      DO i=1,nidx
        IF (woody_pft(i,itile) .AND. (vegetation_height(i,itile) > 0._dp) .AND. (gamma(i) > 0._dp)) then
          scorch_height(i,itile) = lctlib%flame_length_f(cover_type(i,itile)) * (I_surface(i)**0.667_dp)
          ! FireResidenceTime: residence time of fire
!          FireResidenceTime(i,itile) = 2._dp/gamma(i)*com_fac(i,itile) ! gamma: max reaction velocity*mineral damp*moisture damp
           FireResidenceTime(i,itile) = Lethal2ResTime * fuel_consum(i) / 1000._dp / gamma(i)
          ! crownkill
          ! Version using crownkill classes
          crownkill_class          (1:5) = 0._dp
          mortality_crownkill_class(1:5) = 0._dp
          mortality_camb_class     (1:5) = 0._dp
          DO class=1,5
            height       = vegetation_height(i,itile) * (REAL(class,dp) * 0.25_dp + 0.25_dp)  
            crown_length = lctlib%crown_length(cover_type(i,itile)) * height
            stem_length  = height-crown_length
            ! crown kill [0, 1]
            IF (scorch_height(i,itile) < stem_length) THEN
              crownkill_class(class) = 0._dp
            ELSEIF ((scorch_height(i,itile) >= stem_length) .AND. (scorch_height(i,itile) < height)) THEN
              crownkill_class(class) = (scorch_height(i,itile) - stem_length) / crown_length
            ELSE
              crownkill_class(class) = 1._dp
            ENDIF
            mortality_crownkill_class(class) = lctlib%RCK      (cover_type(i,itile))  &
                    * crownkill_class(class)**(lctlib%mort_prob(cover_type(i,itile)))
            ! mortality cambial damage
            ! diameter at breast height
            DBH            = (height / allom2)**(1._dp / allom3) * 100._dp
            ! bark_thickness
            bark_thickness = lctlib%bark_par1(cover_type(i,itile)) * DBH + lctlib%bark_par2(cover_type(i,itile))
            tau_c          = 2.9_dp * bark_thickness**2._dp ! critical time for cambial damage
            ! eq. 19 in TH2010
            IF     (FireResidenceTime(i,itile) / tau_c >= 2._dp) THEN
              mortality_camb_class(class) = 1._dp
            ELSEIF (FireResidenceTime(i,itile) / tau_c >  0.22_dp) THEN
              mortality_camb_class(class) = (0.563_dp * FireResidenceTime(i,itile) / tau_c ) - 0.125_dp
            ELSE
              mortality_camb_class(class) = 0._dp
            ENDIF
          ENDDO ! class
!          ! Version without classes
!          height       = vegetation_height(i,itile)
!          crown_length = lctlib%crown_length(cover_type(i,itile)) * height
!          stem_length  = height - crown_length
!          IF (scorch_height(i,itile) < stem_length) THEN
!            crownkill(i,itile) = 0._dp
!          ELSEIF ((scorch_height(i,itile) >= stem_length) .AND. (scorch_height(i,itile) < height)) THEN
!            crownkill(i,itile) = (scorch_height(i,itile) - stem_length) / crown_length
!          ELSE
!            crownkill(i,itile) = 1._dp
!          ENDIF
!          mortality_crownkill(i,itile) = lctlib%RCK(cover_type(i,itile)) *      &
!                     crownkill(i,itile)**(lctlib%mort_prob(cover_type(i,itile)))
!          ! mortality cambial damage
!          ! diameter at breast height
!          DBH            = (height / allom2)**(1.0_dp / allom3) * 100._dp
!          ! bark_thickness
!          bark_thickness = lctlib%bark_par1(cover_type(i,itile)) * DBH + lctlib%bark_par2(cover_type(i,itile))
!          tau_c          = 2.9_dp * bark_thickness**2.0_dp ! critical time for cambial damage
!          ! eq. 19 in TH2010
!          IF     (FireResidenceTime(i,itile) / tau_c >= 2._dp) THEN
!            mortality_camb(i,itile)=1.0_dp
!          ELSEIF (FireResidenceTime(i,itile) / tau_c >  0.22_dp) THEN
!            mortality_camb(i,itile) = (0.563_dp * FireResidenceTime(i,itile) / tau_c) - 0.125_dp
!          ELSE
!            mortality_camb(i,itile) = 0._dp
!          ENDIF

!          crownkill          (i,itile) = SUM(crownkill_class)           / 5._dp
          mortality_crownkill(i,itile) = SUM(mortality_crownkill_class) / 5._dp
          mortality_camb     (i,itile) = SUM(mortality_camb_class)      / 5._dp
        ENDIF
      ENDDO ! nidx
    ENDDO ! ntiles

    postfire_mortality(:,:) = mortality_crownkill(:,:) + mortality_camb(:,:) &
                            - mortality_crownkill(:,:) * mortality_camb(:,:)
!     postfire_mortality(:,:) = 1.0_dp
    WHERE(.NOT.(woody_pft))
      postfire_mortality = 1._dp
    ENDWHERE
  END SUBROUTINE fire_mortality

  SUBROUTINE relocate_carbon_fire_thonicke(lctlib,surface,nidx,kidx0,ntiles,with_yasso, used_pft, &
                                           cf_burned,cover_fract,veg_fract_correction,            &
                                           Cpool_green, Cpool_woods, Cpool_reserve,               &
                                           Cpool_litter_green_ag,    Cpool_litter_green_bg,       &
                                           Cpool_litter_wood_ag,     Cpool_litter_wood_bg,        &
                                           YCpool_acid_ag1,          YCpool_acid_ag2,             &
                                           YCpool_water_ag1,         YCpool_water_ag2,            &
                                           YCpool_ethanol_ag1,       YCpool_ethanol_ag2,          &
                                           YCpool_nonsoluble_ag1,    YCpool_nonsoluble_ag2,       &
                                           YCpool_acid_bg1,          YCpool_acid_bg2,             &
                                           YCpool_water_bg1,         YCpool_water_bg2,            &
                                           YCpool_ethanol_bg1,       YCpool_ethanol_bg2,          &
                                           YCpool_nonsoluble_bg1,    YCpool_nonsoluble_bg2,       &
                                           YCpool_humus_1,           YCpool_humus_2,              &
                                           LeafLit_coef,             WoodLit_coef,                &
                                           carbon_2_GreenLitterPools,carbon_2_WoodLitterPools,    &
                                           carbon_2_atmos)
  USE mo_jsbach_constants,ONLY: secperday
  USE mo_land_surface,  ONLY: land_surface_type, fract_small
  USE mo_cbal_cpools,   ONLY: frac_wood_aboveGround, frac_green_aboveGround
  TYPE (lctlib_type),       INTENT(in)    :: lctlib
  TYPE (land_surface_type), INTENT(in)    :: surface            ! Access to cover_type
  INTEGER,                  INTENT(in)    :: nidx, kidx0        ! Vector length
  LOGICAL,                  INTENT(in)    :: used_pft(:)        ! PFTs which can burn
  INTEGER,                  INTENT(in)    :: ntiles             ! Number of tiles
  LOGICAL,                  INTENT(in)    :: with_yasso
  REAL(dp),                 INTENT(in)    :: cf_burned(:,:)     ! Frac of the veg area burned since the last call of this routine
  REAL(dp),                 INTENT(in)    :: cover_fract(:,:)   ! Cover fractions
  REAL(dp),                 INTENT(in)    :: veg_fract_correction(:,:)  ! Sparseness of vegetation)
  REAL(dp),                 INTENT(inout) :: Cpool_green(:,:)           ! green carbon pool [mol(C)/m2(canopy)]
  REAL(dp),                 INTENT(inout) :: Cpool_woods(:,:)           ! wood carbon pool [mol(C)/m2(canopy)]
  REAL(dp),                 INTENT(inout) :: Cpool_reserve(:,:)         ! reserve carbon pool [mol(C)/m2(canopy)]
  REAL(dp),                 INTENT(inout) :: Cpool_litter_green_ag(:,:) ! above ground green litter carbon pool [mol(C)/m2(canopy)]
  REAL(dp),                 INTENT(inout) :: Cpool_litter_green_bg(:,:) ! below ground green litter carbon pool [mol(C)/m2(canopy)]
  REAL(dp),                 INTENT(inout) :: Cpool_litter_wood_ag(:,:)  ! wood litter carbon pool [mol(C)/m2(canopy)]
  REAL(dp),                 INTENT(inout) :: Cpool_litter_wood_bg(:,:)  ! wood litter carbon pool [mol(C)/m2(canopy)]
  REAL(dp),                 INTENT(inout) :: YCpool_acid_ag1      (:,:), YCpool_acid_ag2       (:,:)
  REAL(dp),                 INTENT(inout) :: YCpool_water_ag1     (:,:), YCpool_water_ag2      (:,:)
  REAL(dp),                 INTENT(inout) :: YCpool_ethanol_ag1   (:,:), YCpool_ethanol_ag2    (:,:)
  REAL(dp),                 INTENT(inout) :: YCpool_nonsoluble_ag1(:,:), YCpool_nonsoluble_ag2 (:,:)
  REAL(dp),                 INTENT(inout) :: YCpool_acid_bg1      (:,:), YCpool_acid_bg2       (:,:)
  REAL(dp),                 INTENT(inout) :: YCpool_water_bg1     (:,:), YCpool_water_bg2      (:,:)
  REAL(dp),                 INTENT(inout) :: YCpool_ethanol_bg1   (:,:), YCpool_ethanol_bg2    (:,:)
  REAL(dp),                 INTENT(inout) :: YCpool_nonsoluble_bg1(:,:), YCpool_nonsoluble_bg2 (:,:)
  REAL(dp),                 INTENT(inout) :: YCpool_humus_1       (:,:), YCpool_humus_2        (:,:)
  REAL(dp),                 INTENT(in)    :: LeafLit_coef       (:,:,:), WoodLit_coef        (:,:,:)
  REAL(dp),                 INTENT(inout) :: carbon_2_GreenLitterPools(:) ! carbon to green litter pools[mol(C)/m2(vegetated area)]
  REAL(dp),                 INTENT(inout) :: carbon_2_WoodLitterPools(:)  ! carbon to wood  litter pools[mol(C)/m2(vegetated area)]
  REAL(dp),                 INTENT(out)   :: carbon_2_atmos(:)            ! C to atmosphere [mol(C)/m2(vegetated area)/timestep]
    
  ! local variables
  REAL(dp) :: cpools_total_bf(nidx,ntiles)                       
  REAL(dp) :: cpools_total_af(nidx,ntiles)                       
  REAL(dp) :: carbon_2_GreenLitterPools_tiled(nidx,ntiles)
  REAL(dp) :: carbon_2_WoodLitterPools_tiled(nidx,ntiles)
  REAL(dp) :: carbon_2_atmos_tiled(nidx,ntiles)
  REAL(dp) :: green_2_litter_green_ag(nidx,ntiles),green_2_litter_green_bg(nidx,ntiles)
  REAL(dp) :: reserve_2_litter_green_ag(nidx,ntiles),reserve_2_litter_green_bg(nidx,ntiles)
  REAL(dp) :: wood_2_litter_wood_ag(nidx,ntiles),wood_2_litter_wood_bg(nidx,ntiles)
  REAL(dp) :: wood_ck_2_litter_wood_ag, wood_camb_2_litter_wood_ag
!  REAL(dp) :: green_ck_2_litter_green_ag, 
  REAL(dp) :: green_camb_2_litter_green_ag
!  REAL(dp) :: leaves_2_green_litter_pool, 
  REAL(dp) :: f1hr_2_wood_litter_pool, f10hr_2_wood_litter_pool
  REAL(dp) :: f100hr_2_wood_litter_pool, f1000hr_2_wood_litter_pool
!  REAL(dp) :: leaves_old_pool
  REAL(dp) :: mortality_nocrown,cc_ck_wood, cc_ck_green, cc_green
  REAL(dp) :: f1hr_old_pool, f10hr_old_pool, f100hr_old_pool, f1000hr_old_pool
  REAL(dp) :: cpool_litter_wood
  REAL(dp) :: sum_fractions(nidx,ntiles), YCrem_frac(nidx,ntiles)
  REAL(dp) :: scale_fac(nidx,ntiles) , sum_scale_fac
  INTEGER  :: i, ki, itile, ct, kidx1

    ! initialisation
    kidx1 = kidx0 + nidx - 1
    carbon_2_GreenLitterPools_tiled(:,:) = 0._dp
    carbon_2_WoodLitterPools_tiled (:,:) = 0._dp
    green_2_litter_green_ag        (:,:) = 0._dp
    reserve_2_litter_green_ag      (:,:) = 0._dp
    wood_2_litter_wood_ag          (:,:) = 0._dp

    ! Total carbon which must be preserved
    IF (with_yasso) THEN
      cpools_total_bf(:,:) = Cpool_green(:,:) + Cpool_woods(:,:) + Cpool_reserve(:,:) &
                           + YCpool_acid_ag1      (:,:) + YCpool_acid_ag2       (:,:) &
                           + YCpool_water_ag1     (:,:) + YCpool_water_ag2      (:,:) &
                           + YCpool_ethanol_ag1   (:,:) + YCpool_ethanol_ag2    (:,:) &
                           + YCpool_nonsoluble_ag1(:,:) + YCpool_nonsoluble_ag2 (:,:) &
                           + YCpool_acid_bg1      (:,:) + YCpool_acid_bg2       (:,:) &
                           + YCpool_water_bg1     (:,:) + YCpool_water_bg2      (:,:) &
                           + YCpool_ethanol_bg1   (:,:) + YCpool_ethanol_bg2    (:,:) &
                           + YCpool_nonsoluble_bg1(:,:) + YCpool_nonsoluble_bg2 (:,:) &
                           + YCpool_humus_1       (:,:) + YCpool_humus_2        (:,:)
    ELSE
      cpools_total_bf(:,:) = Cpool_green(:,:) + Cpool_woods(:,:) + Cpool_reserve(:,:) &
                           + Cpool_litter_green_ag(:,:) +   Cpool_litter_wood_ag(:,:) &
                           + Cpool_litter_green_bg(:,:) +   Cpool_litter_wood_bg(:,:)
    ENDIF

    ! below ground carbon, same for woody and grass types, postfire_mortality=1 for grasses (initialisation)
    ! also the update of below ground litter pools is the same: litter below ground is not burned
    green_2_litter_green_bg  (:,:) = (1._dp-frac_green_aboveGround)*cf_burned(:,:)*state_FTH%postfire_mortality(kidx0:kidx1,:) &
                                   * Cpool_green  (:,:)    
    reserve_2_litter_green_bg(:,:) = (1._dp-frac_green_aboveGround)*cf_burned(:,:)*state_FTH%postfire_mortality(kidx0:kidx1,:) &
                                   * Cpool_reserve(:,:)
    wood_2_litter_wood_bg    (:,:) = (1._dp-frac_wood_aboveGround )*cf_burned(:,:)*state_FTH%postfire_mortality(kidx0:kidx1,:) &
                                   * Cpool_woods  (:,:)

    ! combustion coefficient for living biomass of killed trees (see TH2010 appendix B2)
    ! cc for wood pool when crown kill occurs: above ground, 100% of twigs and small branches (1hr), 5% big branches (10hr)
    cc_ck_wood  = frac_wood_aboveGround  * (frac_1hr_wood_new + 0.05_dp * frac_10hr_wood_new)
    ! combustion coeffiecient for green pool (leaves and sap wood):
    ! all leaves, 1hr, and 5% of 10hr fuels above ground are burned
    cc_ck_green = 1.0_dp

    DO itile=1,ntiles
      DO i=1,nidx
        IF (cf_burned(i,itile) > 0.5_dp * fract_small) THEN
          ki = i + kidx0 - 1
          ct = surface%cover_type(ki,itile)
          mortality_nocrown = state_FTH%postfire_mortality(ki,itile) - state_FTH%mortality_crownkill(ki,itile)

          IF (lctlib%woody_pft(ct) .AND. surface%is_naturalVeg(ki,itile)) THEN

            ! add biomass burned of woody species:crown_kill: all leaves, all 1, and 5% 10 hour woods
            state_FTH%burned_Cpool_green  (ki,itile) = state_FTH%mortality_crownkill(ki,itile) * cf_burned(i,itile) &
                                                     * Cpool_green  (i,itile) * cc_ck_green
            state_FTH%burned_Cpool_wood   (ki,itile) = state_FTH%mortality_crownkill(ki,itile) * cf_burned(i,itile) &
                                                     * Cpool_woods  (i,itile) * cc_ck_wood
            state_FTH%burned_Cpool_reserve(ki,itile) = state_FTH%mortality_crownkill(ki,itile) * cf_burned(i,itile) &
                                                     * Cpool_reserve(i,itile) * frac_green_aboveGround

            ! carbon affected by fire and transferred into the wood litter pool - diagnostic needed for C conversion test in JSBACH
            wood_ck_2_litter_wood_ag       = state_FTH%mortality_crownkill(ki,itile) * Cpool_woods(i,itile)              &
                                           * frac_wood_aboveGround * cf_burned(i,itile)                              &
                                           * (0.95_dp * frac_10hr_wood_new + frac_100hr_wood_new + frac_1000hr_wood_new)
            wood_camb_2_litter_wood_ag     = mortality_nocrown*Cpool_woods(i,itile)*cf_burned(i,itile)*frac_wood_aboveGround
            wood_2_litter_wood_ag(i,itile) = wood_ck_2_litter_wood_ag + wood_camb_2_litter_wood_ag

            ! carbon transferred into the green litter pool above ground
            green_camb_2_litter_green_ag       = mortality_nocrown * Cpool_green(i,itile)        &
                                               * frac_green_aboveGround * cf_burned(i,itile)
            reserve_2_litter_green_ag(i,itile) = mortality_nocrown * Cpool_reserve(i,itile)      &
                                               * frac_green_aboveGround * cf_burned(i,itile)
            green_2_litter_green_ag  (i,itile) = green_camb_2_litter_green_ag ! crownkill kills all leaves 


            ! available wood litter
            IF (with_yasso) THEN
              cpool_litter_wood = YCpool_acid_ag2   (i,itile) + YCpool_water_ag2     (i,itile) &
                                + YCpool_ethanol_ag2(i,itile) + YCpool_nonsoluble_ag2(i,itile)
            ELSE
              cpool_litter_wood = Cpool_litter_wood_ag(i,itile)
            ENDIF

            ! new 1hr fuel to wood litter pool
            f1hr_2_wood_litter_pool     = frac_1hr_wood_new * wood_camb_2_litter_wood_ag
            f1hr_old_pool               = state_FTH%frac_1hr_wood(ki,itile) * cpool_litter_wood

            ! new 10hr fuel to wood litter pool
            f10hr_2_wood_litter_pool    = frac_10hr_wood_new        * wood_camb_2_litter_wood_ag                            &
                                        + state_FTH%mortality_crownkill(ki,itile) * Cpool_woods(i,itile)                    &
                                          * frac_wood_aboveGround   * cf_burned(i,itile) * frac_10hr_wood_new * 0.95_dp
            f10hr_old_pool              = state_FTH%frac_10hr_wood(ki,itile)   * cpool_litter_wood

            ! new 100hr fuel to wood litter pool
            f100hr_2_wood_litter_pool   = frac_100hr_wood_new       * wood_camb_2_litter_wood_ag                   &
                                        + state_FTH%mortality_crownkill(ki,itile) * Cpool_woods(i,itile)           &
                                          * frac_wood_aboveGround   * cf_burned(i,itile) * frac_100hr_wood_new
            f100hr_old_pool             = state_FTH%frac_100hr_wood(ki,itile)  * cpool_litter_wood

            ! new 1000hr fuel to wood litter pool
            f1000hr_2_wood_litter_pool  = frac_1000hr_wood_new      * wood_camb_2_litter_wood_ag                    &
                                        + state_FTH%mortality_crownkill(ki,itile) * Cpool_woods(i,itile)            &
                                          * frac_wood_aboveGround   * cf_burned(i,itile) * frac_1000hr_wood_new
            f1000hr_old_pool            = state_FTH%frac_1000hr_wood(ki,itile) * cpool_litter_wood

            ! Update fractions
            IF ((wood_2_litter_wood_ag(i,itile) + cpool_litter_wood) > 0._dp) THEN
              state_FTH%frac_1hr_wood   (ki,itile) = (f1hr_2_wood_litter_pool        + f1hr_old_pool)    &
                                                   / (wood_2_litter_wood_ag(i,itile) + cpool_litter_wood)
              state_FTH%frac_10hr_wood  (ki,itile) = (f10hr_2_wood_litter_pool       + f10hr_old_pool)   &
                                                   / (wood_2_litter_wood_ag(i,itile) + cpool_litter_wood)
              state_FTH%frac_100hr_wood (ki,itile) = (f100hr_2_wood_litter_pool      + f100hr_old_pool)  &
                                                   / (wood_2_litter_wood_ag(i,itile) + cpool_litter_wood)
              state_FTH%frac_1000hr_wood(ki,itile) = (f1000hr_2_wood_litter_pool     + f1000hr_old_pool) &
                                                   / (wood_2_litter_wood_ag(i,itile) + cpool_litter_wood)
            ELSE
              state_FTH%frac_1hr_wood   (ki,itile) = frac_1hr_wood_new
              state_FTH%frac_10hr_wood  (ki,itile) = frac_10hr_wood_new
              state_FTH%frac_100hr_wood (ki,itile) = frac_100hr_wood_new
              state_FTH%frac_1000hr_wood(ki,itile) = frac_1000hr_wood_new
            ENDIF

          END IF ! Natural woody pft

          IF (.NOT. lctlib%woody_pft(ct) .AND. used_pft(ct)) THEN
            green_2_litter_green_ag(i,itile) = Cpool_green(i,itile) * frac_green_aboveGround * cf_burned(i,itile) &
                                             - state_FTH%burned_Cpool_green(ki,itile)
            ! combustion completeness of grasses, burned_Cpool_green is computed in fuel_consumtion
            IF (state_FTH%burned_Cpool_green(ki,itile) > fract_small) THEN
              cc_green=state_FTH%burned_Cpool_green(ki,itile)/(Cpool_green(i,itile)*frac_green_aboveGround*cf_burned(i,itile))
              reserve_2_litter_green_ag     (i,itile) =      Cpool_reserve(i,itile)*frac_green_aboveGround*cf_burned(i,itile) &
                                                      * (1._dp-cc_green)
              state_FTH%burned_Cpool_reserve(i,itile) = Cpool_reserve(i,itile)*frac_green_aboveGround*cc_green*cf_burned(i,itile)
            ELSE
              reserve_2_litter_green_ag(i,itile) = Cpool_reserve(i,itile) * frac_green_aboveGround * cf_burned(i,itile)
              state_FTH%burned_Cpool_reserve(i,itile) = 0._dp
            ENDIF
          ENDIF ! not woody and not crop
        ENDIF ! PFT burning
      ENDDO
    ENDDO 

    carbon_2_GreenLitterPools_tiled(:,:) = green_2_litter_green_ag  (:,:) + green_2_litter_green_bg  (:,:)  &
                                         + reserve_2_litter_green_ag(:,:) + reserve_2_litter_green_bg(:,:)
    carbon_2_WoodLitterPools_tiled(:,:)  = wood_2_litter_wood_ag    (:,:) + wood_2_litter_wood_bg    (:,:)

    ! litter consumed and live biomass burned is transferred to atmosphere
    carbon_2_atmos_tiled(:,:) = state_FTH%burned_Cpool_litter_wood_ag (kidx0:kidx1,:) &
                              + state_FTH%burned_Cpool_litter_green_ag(kidx0:kidx1,:) &
                              + state_FTH%burned_Cpool_green          (kidx0:kidx1,:) &
                              + state_FTH%burned_Cpool_wood           (kidx0:kidx1,:) &
                              + state_FTH%burned_Cpool_reserve        (kidx0:kidx1,:)

    Cpool_reserve(:,:) = Cpool_reserve(:,:) - (reserve_2_litter_green_bg(:,:) + reserve_2_litter_green_ag(:,:) &
                                               + state_FTH%burned_Cpool_reserve(kidx0:kidx1,:))
    Cpool_woods  (:,:) = Cpool_woods  (:,:) - (carbon_2_WoodLitterPools_tiled(:,:) + state_FTH%burned_Cpool_wood(kidx0:kidx1,:))
    Cpool_green  (:,:) =  Cpool_green (:,:) &
                       - (green_2_litter_green_ag(:,:)+green_2_litter_green_bg(:,:)+state_FTH%burned_Cpool_green(kidx0:kidx1,:))
    IF (with_yasso) THEN
      ! Aboveground green pools
      YCrem_frac(:,:) = (1._dp - state_FTH%burned_Cpool_litter_green_ag(kidx0:kidx1,:) &
        / MAX(1.e-30_dp,YCpool_acid_ag1(:,:) + YCpool_water_ag1(:,:) + YCpool_ethanol_ag1(:,:) + YCpool_nonsoluble_ag1(:,:)))
      YCpool_acid_ag1      (:,:) = YCpool_acid_ag1      (:,:) * YCrem_frac(:,:) + LeafLit_coef(:,:,1) &
                                 * (green_2_litter_green_ag(:,:) + reserve_2_litter_green_ag(:,:))
      YCpool_water_ag1     (:,:) = YCpool_water_ag1     (:,:) * YCrem_frac(:,:) + LeafLit_coef(:,:,2) &
                                 * (green_2_litter_green_ag(:,:)+reserve_2_litter_green_ag(:,:))
      YCpool_ethanol_ag1   (:,:) = YCpool_ethanol_ag1   (:,:) * YCrem_frac(:,:) + LeafLit_coef(:,:,3) &
                                 * (green_2_litter_green_ag(:,:)+reserve_2_litter_green_ag(:,:))
      YCpool_nonsoluble_ag1(:,:) = YCpool_nonsoluble_ag1(:,:) * YCrem_frac(:,:) + LeafLit_coef(:,:,4) &
                                 * (green_2_litter_green_ag(:,:)+reserve_2_litter_green_ag(:,:))
      ! Aboveground wood pools
      YCrem_frac(:,:) = (1._dp - state_FTH%burned_Cpool_litter_wood_ag(kidx0:kidx1,:) &
        / MAX(1.e-30_dp, YCpool_acid_ag2(:,:) + YCpool_water_ag2(:,:) + YCpool_ethanol_ag2(:,:) + YCpool_nonsoluble_ag2(:,:)))
      YCpool_acid_ag2      (:,:) = YCpool_acid_ag2      (:,:) * YCrem_frac(:,:) &
                                 + WoodLit_coef(:,:,1) * wood_2_litter_wood_ag(:,:)
      YCpool_water_ag2     (:,:) = YCpool_water_ag2     (:,:) * YCrem_frac(:,:) &
                                 + WoodLit_coef(:,:,2) * wood_2_litter_wood_ag(:,:)
      YCpool_ethanol_ag2   (:,:) = YCpool_ethanol_ag2   (:,:) * YCrem_frac(:,:) &
                                 + WoodLit_coef(:,:,3) * wood_2_litter_wood_ag(:,:)
      YCpool_nonsoluble_ag2(:,:) = YCpool_nonsoluble_ag2(:,:) * YCrem_frac(:,:) &
                                 + WoodLit_coef(:,:,4) * wood_2_litter_wood_ag(:,:)
      ! Belowground green pools
      YCpool_acid_bg1      (:,:) = YCpool_acid_bg1      (:,:) + LeafLit_coef(:,:,1) &
                                 * (green_2_litter_green_bg(:,:) + reserve_2_litter_green_bg(:,:))
      YCpool_water_bg1     (:,:) = YCpool_water_bg1     (:,:) + LeafLit_coef(:,:,2) &
                                 * (green_2_litter_green_bg(:,:) + reserve_2_litter_green_bg(:,:))
      YCpool_ethanol_bg1   (:,:) = YCpool_ethanol_bg1   (:,:) + LeafLit_coef(:,:,3) &
                                 * (green_2_litter_green_bg(:,:) + reserve_2_litter_green_bg(:,:))
      YCpool_nonsoluble_bg1(:,:) = YCpool_nonsoluble_bg1(:,:) + LeafLit_coef(:,:,4) &
                                 * (green_2_litter_green_bg(:,:)+reserve_2_litter_green_bg(:,:))
      YCpool_humus_1       (:,:) = YCpool_humus_1       (:,:) + LeafLit_coef(:,:,5) &
                                 * (green_2_litter_green_bg(:,:)+reserve_2_litter_green_bg(:,:))
      ! Belowground wood pools
      YCpool_acid_bg2      (:,:) = YCpool_acid_bg2      (:,:) + WoodLit_coef(:,:,1) * wood_2_litter_wood_bg(:,:)
      YCpool_water_bg2     (:,:) = YCpool_water_bg2     (:,:) + WoodLit_coef(:,:,2) * wood_2_litter_wood_bg(:,:)
      YCpool_ethanol_bg2   (:,:) = YCpool_ethanol_bg2   (:,:) + WoodLit_coef(:,:,3) * wood_2_litter_wood_bg(:,:)
      YCpool_nonsoluble_bg2(:,:) = YCpool_nonsoluble_bg2(:,:) + WoodLit_coef(:,:,4) * wood_2_litter_wood_bg(:,:)
      YCpool_humus_2       (:,:) = YCpool_humus_2       (:,:) + WoodLit_coef(:,:,5) * wood_2_litter_wood_bg(:,:)
    ELSE
      Cpool_litter_green_ag(:,:) = Cpool_litter_green_ag(:,:) + green_2_litter_green_ag(:,:) + reserve_2_litter_green_ag(:,:) &
                                                              - state_FTH%burned_Cpool_litter_green_ag(kidx0:kidx1,:)
      Cpool_litter_green_bg(:,:) = Cpool_litter_green_bg(:,:) + green_2_litter_green_bg(:,:) + reserve_2_litter_green_bg(:,:)
      Cpool_litter_wood_ag (:,:) = Cpool_litter_wood_ag (:,:) + wood_2_litter_wood_ag  (:,:)                         &
                                                              - state_FTH%burned_Cpool_litter_wood_ag(kidx0:kidx1,:)
      Cpool_litter_wood_bg (:,:) = Cpool_litter_wood_bg (:,:) + wood_2_litter_wood_bg  (:,:)
    ENDIF

    ! bring to total grid box 
    carbon_2_atmos_tiled           (:,:) = carbon_2_atmos_tiled           (:,:) * veg_fract_correction(:,:) * cover_fract(:,:)
    carbon_2_GreenLitterPools_tiled(:,:) = carbon_2_GreenLitterPools_tiled(:,:) * veg_fract_correction(:,:) * cover_fract(:,:)
    carbon_2_WoodLitterPools_tiled (:,:) = carbon_2_WoodLitterPools_tiled (:,:) * veg_fract_correction(:,:) * cover_fract(:,:)
    
    carbon_2_atmos            = SUM(carbon_2_atmos_tiled,           DIM=2)
    carbon_2_GreenLitterPools = SUM(carbon_2_GreenLitterPools_tiled,DIM=2)
    carbon_2_WoodLitterPools  = SUM(carbon_2_WoodLitterPools_tiled, DIM=2) 

    ! normalize fractions
    sum_fractions   (:,:) = state_FTH%frac_1hr_wood  (kidx0:kidx1,:) + state_FTH%frac_10hr_wood  (kidx0:kidx1,:) &
                          + state_FTH%frac_100hr_wood(kidx0:kidx1,:) + state_FTH%frac_1000hr_wood(kidx0:kidx1,:)
    state_FTH%frac_1hr_wood   (kidx0:kidx1,:) = state_FTH%frac_1hr_wood   (kidx0:kidx1,:) / sum_fractions(:,:)
    state_FTH%frac_10hr_wood  (kidx0:kidx1,:) = state_FTH%frac_10hr_wood  (kidx0:kidx1,:) / sum_fractions(:,:)
    state_FTH%frac_100hr_wood (kidx0:kidx1,:) = state_FTH%frac_100hr_wood (kidx0:kidx1,:) / sum_fractions(:,:)
    state_FTH%frac_1000hr_wood(kidx0:kidx1,:) = state_FTH%frac_1000hr_wood(kidx0:kidx1,:) / sum_fractions(:,:)

    ! Diagnose sum of carbon pools after relocate (without slow pool) for C mass balance check 
    IF (with_yasso) THEN
      cpools_total_af(:,:) = Cpool_green(:,:) + Cpool_woods(:,:) + Cpool_reserve(:,:) &
                           + YCpool_acid_ag1      (:,:) + YCpool_acid_ag2       (:,:) &
                           + YCpool_water_ag1     (:,:) + YCpool_water_ag2      (:,:) &
                           + YCpool_ethanol_ag1   (:,:) + YCpool_ethanol_ag2    (:,:) &
                           + YCpool_nonsoluble_ag1(:,:) + YCpool_nonsoluble_ag2 (:,:) &
                           + YCpool_acid_bg1      (:,:) + YCpool_acid_bg2       (:,:) &
                           + YCpool_water_bg1     (:,:) + YCpool_water_bg2      (:,:) &
                           + YCpool_ethanol_bg1   (:,:) + YCpool_ethanol_bg2    (:,:) &
                           + YCpool_nonsoluble_bg1(:,:) + YCpool_nonsoluble_bg2 (:,:) &
                           + YCpool_humus_1       (:,:) + YCpool_humus_2        (:,:)
    ELSE
      cpools_total_af(:,:) = Cpool_green(:,:) + Cpool_woods(:,:) + Cpool_reserve(:,:) &
                           + Cpool_litter_green_ag(:,:) +   Cpool_litter_wood_ag(:,:) &
                           + Cpool_litter_green_bg(:,:) +   Cpool_litter_wood_bg(:,:)
    ENDIF
    cpools_total_af(:,:) = cpools_total_af(:,:) * veg_fract_correction(:,:) * cover_fract(:,:)
    cpools_total_bf(:,:) = cpools_total_bf(:,:) * veg_fract_correction(:,:) * cover_fract(:,:)

    WHERE(cpools_total_bf(:,:) > 1.0_dp .AND. cover_fract(:,:) > 0.5_dp*fract_small)
      scale_fac(:,:) = (carbon_2_atmos_tiled(:,:) + cpools_total_af(:,:) - cpools_total_bf(:,:)) / cpools_total_bf(:,:)
    ELSEWHERE
      scale_fac(:,:) = 0._dp
    ENDWHERE
    sum_scale_fac = SUM(SUM(scale_fac,DIM=1))
    IF(abs(sum_scale_fac) > 1.e-13_dp) THEN
       WRITE(message_text,*) 'relocate_carbon_fire_thonicke, carbon conservation in relocation rather low (',sum_scale_fac,')'
       CALL message('pool_consumption',TRIM(message_text))
    END IF

    IF (fire_TH_options%ldiag) &
      fire_TH_diag%avg_carbon_2_atmos(kidx0:kidx1,:) = fire_TH_diag%avg_carbon_2_atmos(kidx0:kidx1,:) &
                                                     + secperday * carbon_2_atmos_tiled

  END SUBROUTINE relocate_carbon_fire_thonicke

  SUBROUTINE FRP_nfire_per_gridcell(nidx, ntiles, cover_fract,veg_ratio_max,numfire,burned_frac_diag,&
                                    IR_ROS, FireDuration,cellarea,FRP_gridcell,numfires_gridcell)
  USE mo_land_surface,     ONLY: fract_small
  INTEGER,                  INTENT(in)    :: nidx,ntiles
  REAL(dp),                 INTENT(in)    :: cover_fract(:,:)   ! Cover fractions
  REAL(dp),                 INTENT(in)    :: veg_ratio_max(:)
  REAL(dp),                 INTENT(in)    :: numfire(:)         ! [per ha]
  REAL(dp),                 INTENT(in)    :: burned_frac_diag(:,:)
  REAL(dp),                 INTENT(in)    :: IR_ROS(:)       ! kJ m-2 min-1
  REAL(dp),                 INTENT(in)    :: FireDuration(:) ! min
  REAL(dp),                 INTENT(in)    :: cellarea(:)     ! gridcell area [m2]
  REAL(dp),                 INTENT(out)   :: FRP_gridcell(:),numfires_gridcell(:)! FRP [W], nfire [-]
!local
  REAL(dp)                                :: box_burned_frac_diag(nidx,ntiles)
  REAL(dp),                 PARAMETER     :: Heat2FRP_Ratio = 0.15_dp     !Conversion Factor total heat release -> FRP
  REAL(dp),                 PARAMETER     :: SecPerMin      = 60.0_dp
  REAL(dp),                 PARAMETER     :: kJ2J           = 1000.0_dp
  REAL(dp),                 PARAMETER     :: HoursPerDay    = 24.0_dp
  REAL(dp),                 PARAMETER     :: sqm2ha         = 10000._dp
   
    FRP_gridcell      = 0.0_dp
    numfires_gridcell = 0.0_dp

    box_burned_frac_diag = burned_frac_diag * cover_fract * SPREAD(veg_ratio_max,DIM=2,ncopies=ntiles)

    WHERE(FireDuration > fract_small)
      ! FRP = IR_ROS/60*1000 * BF * area / fireduration/60 * fireduration/60/24 * 0.15
      ! factor fireduration/60/24 converts mean FRP of the gridcell at the time of fire to a 
      ! mean frp per day
      ! same for the number of fires, mean number of fires burning over one day
      FRP_gridcell = IR_ROS   / SecPerMin * kJ2J * SUM(box_burned_frac_diag,DIM=2) & 
                   * cellarea / SecPerMin / SecPerMin / HoursPerDay * Heat2FRP_Ratio 
      numfires_gridcell = numfire * cellarea / sqm2ha * FireDuration / SecPerMin / HoursPerDay 
    ENDWHERE

  END SUBROUTINE FRP_nfire_per_gridcell

END MODULE mo_disturbance_thonicke
!References:
!Thonicke, K., Spessa, A., Prentice, I. C., Harrison, S. P., Dong, L., and Carmona-Moreno, C.: The influence of vegetation,
!fire spread and fire behaviour on biomass burning and trace gas emissions: results from a process-based model, 
!Biogeosciences, 7, 1991-2011, doi:10.5194/bg-7-1991-2010, 2010. 
