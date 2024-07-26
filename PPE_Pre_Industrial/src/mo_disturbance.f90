!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_disturbance
!
! Driver for calculating vegetation disturbances (fire and windbreak)
! using different fire/windbreak models
!
! History:
! 2011-04-29 - StW - Original code supporting jsbach fire-algorithm only
! 2011-12-09 - StW - Accumulated variables changed to averages

USE mo_kind,                 ONLY: dp
USE mo_mpi,                  ONLY: p_parallel_io, p_io, p_parallel, p_bcast
USE mo_exception,            ONLY: finish, message
USE mo_linked_list,          ONLY: t_stream
USE mo_land_surface,         ONLY: land_surface_type
USE mo_jsbach_grid,          ONLY: grid_type, domain_type
USE mo_climbuf,              ONLY: climbuf_type
USE mo_cbal_bethy,           ONLY: cbalance_type, nbalance_type
USE mo_jsbach_lctlib,        ONLY: lctlib_type
USE mo_utils,                ONLY: average_tiles
USE mo_cbal_cpools,          ONLY: relocate_carbon_fire, relocate_carbon_damage
USE mo_disturbance_jsbach,   ONLY: init_fire_jsbach, init_windbreak_jsbach, burned_frac_jsbach, broken_woody_frac_jsbach
USE mo_disturbance_thonicke, ONLY: burned_frac_thonicke,init_fire_thonicke,init_fire_thonicke_memory,relocate_carbon_fire_thonicke

IMPLICIT NONE

TYPE(t_stream),  POINTER :: IO_disturbance, IO_dist_diag   !! Memory streams for disturbance model state

! Submodel identifiers
INTEGER , PARAMETER, PUBLIC  :: DIST_FIRE_WOOD        =  1 !! Calculation wrt. burning woody pfts
INTEGER , PARAMETER, PUBLIC  :: DIST_FIRE_GRASS       =  2 !! Calculation wrt. burning grass pfts
INTEGER , PARAMETER, PUBLIC  :: DIST_FIRE             =  3 !! Calculation wrt. any burning pfts
INTEGER , PARAMETER, PUBLIC  :: DIST_WINDBREAK_WOOD   =  4 !! Calculation wrt. windbread of woody types
INTEGER , PARAMETER, PUBLIC  :: DIST_WINDBREAK_GRASS  =  8 !! Calculation wrt. windbreak of grass types 
                                                           !! (ok: 8 is unlikely to ever be used - maybe for extreme 
                                                           !! precipitation events?!)
INTEGER , PARAMETER, PUBLIC  :: DIST_WINDBREAK        = 12 !! Calculation wrt. windbreak of any pft types 

! ----------------------------------------------------------
! Default parameters for the various algorithms

! JSBACH wind break algorithm parameters
LOGICAL , PARAMETER, PRIVATE :: def_ldiag                  = .FALSE.  !! output all variables?
REAL(dp), PARAMETER, PRIVATE :: def_fire_frac_wood_2_atmos =  0.2_dp  !! fraction of wood carbon emitted to the atm by fire

! ---------------------------------------------------------------------
! Type to hold options for the different fire/windbreak algorithm types
TYPE disturbance_option_type                   !! See additional explanations above.
  INTEGER           :: fire_algorithm          !! Select fire algorithm: 0 : none, 1 : jsbach (default), 
                                               !! 2 : Arora and Boer (2005) (not yet), 3 : Thonicke et al.(2001)
  INTEGER           :: fuel_algorithm          !! Select fuel algorithm: 0 : none, 1 : jsbach (default), 
                                               !! 2 : Arora and Boer (2005), 3 : Thonicke et al.(2001), 4 : Read/GFED
  INTEGER           :: windbreak_algorithm     !! Select windbreak algorithm: 0 : none, 1 : jsbach (default)
  LOGICAL           :: lburn_pasture
  LOGICAL           :: ldiag
  REAL(dp)          :: fire_frac_wood_2_atmos
END TYPE disturbance_option_type

! Type holding the necessary state variables for disturbances
TYPE disturbance_type
  REAL(dp), POINTER   :: damaged_frac(:,:)
  REAL(dp), POINTER   :: burned_frac(:,:)
  REAL(dp), POINTER   :: burned_frac_diag(:,:)
  REAL(dp), POINTER   :: carbon_2_GreenLitterPools(:)
  REAL(dp), POINTER   :: carbon_2_WoodLitterPools(:)
  REAL(dp), POINTER   :: carbon_2_atmos(:)
  REAL(dp), POINTER   :: fuel(:)
  REAL(dp), POINTER   :: box_burned_frac_diag_avg(:,:) ! Temporal average of fraction of box affected by fire (cmp. with obs.)
  REAL(dp), POINTER   :: box_burned_frac_avg(:,:)      ! Temporal average of fraction of box freed due to fire(s)
  REAL(dp), POINTER   :: box_damaged_frac_avg(:,:)
  REAL(dp), POINTER   :: box_CO2_flux_2_atmos(:)
  ! Variables needed with "with_nitrogen"
  REAL(dp), POINTER   :: nitrogen_2_GreenLitterPools(:)
  REAL(dp), POINTER   :: nitrogen_2_WoodLitterPools(:)
  REAL(dp), POINTER   :: nitrogen_2_atmos(:) 
  REAL(dp), POINTER   :: nitrogen_2_sminn(:)   
END TYPE disturbance_type

! Global variables
TYPE (disturbance_option_type), PUBLIC  :: dist_opts   ! Actual set of disturbance parameters
TYPE (disturbance_type),        PUBLIC  :: disturbance ! State variables of disturbances
LOGICAL,                        PRIVATE :: module_initialized = .FALSE. ! Has config_disturbance been called?

! Subroutine declarations
PUBLIC config_disturbance
PUBLIC disturbed_frac
PUBLIC relocate_disturbed_carbon
PUBLIC update_disturbance

CONTAINS

! --- config_disturbance ----------------------------------------------
!
! Reads in parameters from namelist file
!
! History:
! 2011-04-29 - StW - Original code

  SUBROUTINE config_disturbance(surface, grid, domain, lctlib, nml_unit, with_yasso, with_nitrogen, cbalone)
  USE mo_namelist,         ONLY: position_nml, POSITIONED, LENGTH_ERROR, READ_ERROR
  USE mo_io_units,         ONLY: nout
  USE mo_input_strings,    ONLY: GetStringIndex

  ! Input parameters
  TYPE(land_surface_type) , INTENT(in) :: surface
  TYPE(grid_type)         , INTENT(in) :: grid
  TYPE(domain_type)       , INTENT(in) :: domain
  TYPE(lctlib_type)       , INTENT(in) :: lctlib
  INTEGER                 , INTENT(in) :: nml_unit
  LOGICAL                 , INTENT(in) :: with_yasso
  LOGICAL                 , INTENT(in) :: with_nitrogen
  LOGICAL, OPTIONAL       , INTENT(in) :: cbalone

  ! Parameters for determining named parameters
  CHARACTER (len=*), PARAMETER :: fire_names  =  'none jsbach arora thonicke read '
  CHARACTER (len=*), PARAMETER :: fuel_names  =  'none jsbach arora thonicke read '
  CHARACTER (len=*), PARAMETER :: wind_names  =  'none jsbach '

  ! Local namelist variables
  CHARACTER(30) :: fire_name
  CHARACTER(30) :: fuel_name
  CHARACTER(30) :: windbreak_name
  LOGICAL       :: lburn_pasture
  LOGICAL       :: ldiag
  REAL(dp)      :: fire_frac_wood_2_atmos 

  ! Other locals
  INTEGER  :: read_status
  INTEGER  :: ibuf(3)
  INTEGER  :: f_unit
  LOGICAL  :: lcbalone

  INCLUDE 'disturbance_ctl.inc'

    lcbalone = .FALSE.
    IF (PRESENT(cbalone)) lcbalone = cbalone 

    IF (p_parallel_io) THEN
      ! Map default values
      dist_opts%fire_algorithm      =  0
      dist_opts%fuel_algorithm      = -1
      dist_opts%windbreak_algorithm =  0
      fire_name              = 'jsbach'
      fuel_name              = ''
      windbreak_name         = 'jsbach'
      lburn_pasture          = .FALSE.
      ldiag                  = def_ldiag
      fire_frac_wood_2_atmos = def_fire_frac_wood_2_atmos

      ! Read namelist if available
      f_unit = position_nml ('DISTURBANCE_CTL', nml_unit, status=read_status)
      SELECT CASE (read_status)
        CASE (POSITIONED)
          READ (f_unit, disturbance_ctl)
          CALL message('config_disturbance','Namelist DISTURBANCE_CTL: ')
        CASE (LENGTH_ERROR)
          CALL finish ('config_disturbance','Length error in namelist disturbance_ctl')
        CASE (READ_ERROR)
          CALL finish ('config_disturbance','Error reading namelist disturbance_ctl ')
      END SELECT

      ! Process named parameter equivalent to enumerated algorithms
      IF (LEN_TRIM(     fire_name) > 0) &
        dist_opts%fire_algorithm      = GetStringIndex(fire_name,     fire_names,ErrMsg='disturbance_ctl%fire_name'     )-1
      IF (LEN_TRIM(     fuel_name) > 0) &
        dist_opts%fuel_algorithm      = GetStringIndex(fuel_name,     fuel_names,ErrMsg='disturbance_ctl%fuel_name'     )-1
      IF (LEN_TRIM(windbreak_name) > 0) &
        dist_opts%windbreak_algorithm = GetStringIndex(windbreak_name,wind_names,ErrMsg='disturbance_clt%windbreak_name')-1
      IF (dist_opts%fuel_algorithm < 0) dist_opts%fuel_algorithm = dist_opts%fire_algorithm

      ! Final output of namelist
      WRITE(nout, disturbance_ctl)

      ! Tests for sainity of provided parameters
      ! ... individually
      IF (dist_opts%fire_algorithm       < 0  .OR. dist_opts%fire_algorithm      > 4) &
        CALL finish('config_disturbance','Unknown fire algorithm!')
      IF (dist_opts%fuel_algorithm       < 0  .OR. dist_opts%fuel_algorithm      > 4) & 
        CALL finish('config_disturbance','Unknown fuel algorithm!')
      IF (dist_opts%windbreak_algorithm  < 0  .OR. dist_opts%windbreak_algorithm > 1) &
        CALL finish('config_disturbance','Unknown windbreak algorithm!')

      ! ... and combined
      IF ((dist_opts%fire_algorithm==2) .OR. (dist_opts%fuel_algorithm==2)) &
        CALL finish('config_disturbance','Sorry: Fire/fuel algorithm Arora not yet implemented in this version')
!      IF ((dist_opts%fire_algorithm==3) .OR. (dist_opts%fuel_algorithm==3)) &
!        CALL finish('config_disturbance','Sorry: Fire/fuel algorithm Thonicke not yet implemented in this version')
      IF ((dist_opts%fire_algorithm==4) .OR. (dist_opts%fuel_algorithm==4)) &
        CALL finish('config_disturbance','Sorry: Fire/fuel algorithm Read not yet implemented in this version')

      ! Distribute namelist parameters
      ibuf = (/dist_opts%fire_algorithm,dist_opts%fuel_algorithm,dist_opts%windbreak_algorithm/)
    ENDIF

    IF (p_parallel) THEN
      CALL p_bcast(ibuf,p_io)
      CALL p_bcast(fire_frac_wood_2_atmos,p_io)
      CALL p_bcast(lburn_pasture,p_io)
      CALL p_bcast(ldiag,p_io)
    ENDIF
     
    dist_opts%fire_algorithm         = ibuf( 1)
    dist_opts%fuel_algorithm         = ibuf( 2)
    dist_opts%windbreak_algorithm    = ibuf( 3)
    dist_opts%fire_frac_wood_2_atmos = fire_frac_wood_2_atmos
    dist_opts%lburn_pasture          = lburn_pasture
    dist_opts%ldiag                  = ldiag

    ! --- Initialize relevant modules ---
    ! Fire
    SELECT CASE (dist_opts%fire_algorithm)
      CASE (0) !! No fire
      CASE (1) !! jsbach fire
        CALL init_fire_jsbach(lctlib,nml_unit)
      CASE (2) !! Arora & Boer fire
      CASE (3) !! Thonicke fire
         CALL init_fire_thonicke(grid,domain,surface,lcbalone,dist_opts%ldiag,with_yasso,nml_unit)
      CASE (4) !! Read burned_frac
    END SELECT
    ! Wind break    
    SELECT CASE (dist_opts%windbreak_algorithm)
      CASE (0) !! No wind break 
      CASE (1) !! jsbach wind break
        CALL init_windbreak_jsbach(nml_unit)
    END SELECT

    ! Allocate and initialize memory for the disturbance type
    IF (lcbalone) THEN ! For jsbach the allocation is done through the stream functions in init_disturbance
      ALLOCATE(disturbance%burned_frac              (grid%nland,surface%ntiles), &
               disturbance%damaged_frac             (grid%nland,surface%ntiles), &
               disturbance%carbon_2_GreenLitterPools(grid%nland),                &
               disturbance%carbon_2_WoodLitterPools (grid%nland),                &
               disturbance%carbon_2_atmos           (grid%nland),                &
               disturbance%box_burned_frac_avg      (grid%nland,surface%ntiles), &
               disturbance%box_damaged_frac_avg     (grid%nland,surface%ntiles), &
               disturbance%fuel                     (grid%nland),                &
               disturbance%box_CO2_flux_2_atmos     (grid%nland))

      disturbance%burned_frac                (:,:) = 0._dp
      disturbance%damaged_frac               (:,:) = 0._dp
      disturbance%carbon_2_GreenLitterPools  (:)   = 0._dp
      disturbance%carbon_2_WoodLitterPools   (:)   = 0._dp
      disturbance%carbon_2_atmos             (:)   = 0._dp
      disturbance%box_burned_frac_avg        (:,:) = 0._dp
      disturbance%box_damaged_frac_avg       (:,:) = 0._dp
      disturbance%box_CO2_flux_2_atmos       (:)   = 0._dp
      disturbance%fuel                       (:)   = 0._dp

      IF (with_nitrogen) THEN
        ALLOCATE(disturbance%nitrogen_2_GreenLitterPools(grid%nland), &
                 disturbance%nitrogen_2_WoodLitterPools (grid%nland), &
                 disturbance%nitrogen_2_atmos           (grid%nland), &
                 disturbance%nitrogen_2_sminn           (grid%nland))
        disturbance%nitrogen_2_GreenLitterPools(:)   = 0._dp
        disturbance%nitrogen_2_WoodLitterPools (:)   = 0._dp
        disturbance%nitrogen_2_atmos           (:)   = 0._dp
        disturbance%nitrogen_2_sminn           (:)   = 0._dp
      ENDIF

      IF (dist_opts%fire_algorithm == 3) THEN ! Thonicke uses sep. diag frac
        ALLOCATE(disturbance%burned_frac_diag        (grid%nland,surface%ntiles), &
                 disturbance%box_burned_frac_diag_avg(grid%nland,surface%ntiles))
        disturbance%burned_frac_diag        (:,:) = 0._dp
        disturbance%box_burned_frac_diag_avg(:,:) = 0._dp
      ELSE
        disturbance%burned_frac_diag         => disturbance%burned_frac
        disturbance%box_burned_frac_diag_avg => disturbance%box_burned_frac_avg
      ENDIF
    ENDIF

  END SUBROUTINE config_disturbance

  ! --- init_disturbance --------------------------------------------------------------------------------------------------------
  !
  ! Initialises this disturbance module. In particular the structure "disturbance" is initialised,
  ! this means, that the initialisation of disturbance is done here.
  
  SUBROUTINE init_disturbance(grid,domain,lctlib,surface,lrestart,with_yasso,with_nitrogen, &
                              nml_unit,fileformat,fileztype,stream,dist_stream)
  USE mo_linked_list,            ONLY: LAND, TILES
  USE mo_memory_base,            ONLY: new_stream,default_stream_setting, &
                                       add =>add_stream_element
  USE mo_output,                 ONLY: veg_table, dist_table
  USE mo_netcdf,                 ONLY: max_dim_name
  USE mo_jsbach,                 ONLY: missing_value, veg_putdata
  TYPE(grid_type),   INTENT(in)       :: grid
  TYPE(domain_type), INTENT(in)       :: domain
  TYPE(lctlib_type), INTENT(in)       :: lctlib
  TYPE(land_surface_type), INTENT(in) :: surface
  LOGICAL,           INTENT(in)       :: lrestart      ! true for restarted runs (throughout the whole run)
  INTEGER,           INTENT(in)       :: nml_unit    
  INTEGER,           INTENT(in)       :: fileformat    ! output file format
  INTEGER,           INTENT(in)       :: fileztype     ! output file compression
  LOGICAL,           INTENT(in)       :: with_yasso
  LOGICAL,           INTENT(in)       :: with_nitrogen
  TYPE(t_stream),   POINTER, OPTIONAL :: stream, dist_stream

  ! local variables
  INTEGER                 :: g_nland, l_nland
  INTEGER                 :: dim1p(1), dim1(1)
  INTEGER                 :: dim2p(2), dim2(2)
  CHARACTER(max_dim_name) :: dim1n(1), dim2n(2)
    
    ! Check initialization
    IF (module_initialized) CALL message("init_disturbance","Warning: Re-initializing module.")

    ! Read and distribute disturbance configuration
    CALL config_disturbance(surface, grid, domain, lctlib, nml_unit, with_yasso, with_nitrogen)

    ! Stream definition
    IF (PRESENT(stream)) THEN 
      IF (.NOT. ASSOCIATED(stream)) THEN
        ! Add new stream
        CALL new_stream(stream, 'disturbance', filetype=fileformat, ztype=fileztype)
        ! Set default stream options
        CALL default_stream_setting(stream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
      ENDIF
      IO_disturbance => stream
    ELSE
      ! Add new stream
      CALL new_stream(IO_disturbance, 'disturbance', filetype=fileformat, ztype=fileztype)
      ! Set default stream options
      CALL default_stream_setting(IO_disturbance, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    ! Add new stream specially for disturbances if extra variables are needed
    IF (dist_opts%ldiag .AND. ((dist_opts%fire_algorithm==2) .OR. (dist_opts%fire_algorithm==3) &
                                                             .OR. (dist_opts%fire_algorithm==4))) THEN 
      CALL new_stream(IO_dist_diag, 'disturbance', filetype=fileformat, ztype=fileztype, &
                      lpost=.TRUE., lrerun=.TRUE., lcontnorest=.TRUE.,interval=veg_putdata)
      CALL default_stream_setting(IO_dist_diag, lpost=.TRUE., lrerun=.TRUE., contnorest=.TRUE., &
                                  repr=LAND, table=dist_table, leveltype=TILES)
      IF (PRESENT(dist_stream)) dist_stream => IO_dist_diag
    ENDIF

    g_nland = grid%nland
    l_nland = domain%nland

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim2p = (/ l_nland, surface%ntiles /)
    dim2  = (/ g_nland, surface%ntiles /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'tiles'

    ! --- Define variables as stream elements
    ! code numbers of this routine range from 50 to 63 using the GRIB veg_table
    ! note that more variables may be added to stream in the various disturbance modules
    CALL add(IO_disturbance,'fuel',disturbance%fuel,longname='fuel used for burned area calculations', units='mol(C)/m^2(veg)', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=51, lpost=.TRUE., lrerun = .TRUE.,lmiss=.TRUE.,missval=missing_value)
    CALL add(IO_disturbance,'burned_frac'     ,disturbance%burned_frac, longname='burned area fraction',   units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=52, lpost=.true., lmiss=.TRUE., missval=missing_value, &
             contnorest=.true.)
    CALL add(IO_disturbance,'box_burned_frac_avg',disturbance%box_burned_frac_avg, lpost=.true., &
             longname='Burned area fraction averaged over output period', units='m2m-2(grid box)', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=53, laccu=.true., lrerun=.TRUE., lmiss=.TRUE., &
             missval=missing_value, contnorest=.true.)
    CALL add(IO_disturbance,'damaged_frac'    ,disturbance%damaged_frac,longname='damaged area fraction',  units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=54, lpost=.true., lmiss=.TRUE., missval=missing_value, &
             contnorest=.true.)
    CALL add(IO_disturbance,'box_damaged_frac_avg',disturbance%box_damaged_frac_avg, lpost=.true., &
             longname='Damaged area fraction averaged over output period', units='m2m-2(grid box)', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=55, laccu=.true., lrerun=.FALSE., lmiss=.TRUE., &
             missval=missing_value, contnorest=.true.)
    CALL add(IO_disturbance,'box_fire_CO2_flux_2_atmos',disturbance%box_CO2_flux_2_atmos, contnorest=.TRUE., &
             lpost=.true.,longname='Average CO2 flux to atmosphere',   units='kg(CO2) m-2(grid box) s-1', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=56, laccu=.true., lrerun=.TRUE., lmiss=.TRUE., &
             missval=missing_value)
    CALL add(IO_disturbance,'dist_carbon_2_GreenLitterPools', disturbance%carbon_2_GreenLitterPools, lpost=.false., &
             longname='C transferred to green litter (ag+bg) by vegetation disturbances', &
             units='mol(C) m-2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
             code=57, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value, contnorest=.true.)
    CALL add(IO_disturbance,'dist_carbon_2_WoodLitterPools', disturbance%carbon_2_WoodLitterPools, lpost=.false., &
             longname='C transferred to wood litter (ag+bg) by vegetation disturbances', &
             units='mol(C) m-2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
             code=58, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value, contnorest=.true.)
    CALL add(IO_disturbance,'dist_carbon_2_atmos', disturbance%carbon_2_atmos, &
             longname='CO2 emitted to the atmosphere by fire', &
             units='mol(C) m-2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
             code=59, lpost=.false., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    IF (with_nitrogen) THEN
      CALL add(IO_disturbance,'dist_nitrogen_2_GreenLitterPools', disturbance%nitrogen_2_GreenLitterPools, &
               longname='N transferred to green litter (ag+bg) by vegetation disturbances', &
               units='mol(N) m-2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
               code=60, lpost=.false., lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
      CALL add(IO_disturbance,'dist_nitrogen_2_WoodLitterPools', disturbance%nitrogen_2_WoodLitterPools, &
               longname='N transferred to wood litter (ag+bg) by vegetation disturbances', &
               units='mol(N) m-2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
               code=61, lpost=.false., lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
      CALL add(IO_disturbance,'dist_nitrogen_2_atmos', disturbance%nitrogen_2_atmos, &
               longname='N compounds emitted to the atmosphere by fire', &
               units='mol(N) m-2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
               code=62, lpost=.false., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
      CALL add(IO_disturbance,'dist_nitrogen_2_sminn', disturbance%nitrogen_2_sminn, &
               longname='N compounds transferred to mineral N pool by fire', &
               units='mol(N) m-2(grid box)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
               code=63, lpost=.false., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    ENDIF
    IF (dist_opts%fire_algorithm==3) THEN ! Thonicke needs extra variables
      CALL add(IO_disturbance,'burned_frac_diag',disturbance%burned_frac_diag,ldims=dim2p, gdims=dim2, dimnames=dim2n, &
               longname='burned area fraction including survivors',units='', lpost=.false., lmiss=.TRUE., &
               missval=missing_value, contnorest=.true.)
      CALL add(IO_disturbance,'box_burned_frac_diag_avg',disturbance%box_burned_frac_diag_avg, lpost=.true., &
               longname='Burned area fraction including survivors', units='m2m-2(grid box)', contnorest=.true., &
               ldims=dim2p, gdims=dim2, dimnames=dim2n, code=50, laccu=.true., lmiss=.TRUE., missval=missing_value)
    ELSE
      disturbance%burned_frac_diag         => disturbance%burned_frac
      disturbance%box_burned_frac_diag_avg => disturbance%box_burned_frac_avg
    ENDIF

    ! -- Add variables local to the individual disturbance algorithms. Make sure that no grib-code for any fire algorithm overlap
    ! Fire
    SELECT CASE (dist_opts%fire_algorithm)
      CASE (0) !! No fire
      CASE (1) !! jsbach fire
      CASE (2) !! Arora & Boer fire
      CASE (3) !! Thonicke fire
        CALL init_fire_thonicke_memory(IO_dist_diag,IO_disturbance,dim1p,dim1,dim1n,dim2p,dim2,dim2n,with_yasso)
      CASE (4) !! Read burned_frac
    END SELECT
    ! Wind break
    SELECT CASE (dist_opts%windbreak_algorithm)
      CASE (0) !! No wind break
      CASE (1) !! jsbach wind break
    END SELECT 

    disturbance%box_burned_frac_avg     (:,:) = 0._dp
    disturbance%box_burned_frac_diag_avg(:,:) = 0._dp
    disturbance%box_damaged_frac_avg    (:,:) = 0._dp
    disturbance%carbon_2_GreenLitterPools (:) = 0._dp
    disturbance%carbon_2_WoodLitterPools  (:) = 0._dp
    disturbance%carbon_2_atmos            (:) = 0._dp
    disturbance%box_CO2_flux_2_atmos      (:) = 0._dp
    IF (with_nitrogen) THEN
      disturbance%nitrogen_2_GreenLitterPools(:) = 0._dp
      disturbance%nitrogen_2_WoodLitterPools (:) = 0._dp
      disturbance%nitrogen_2_atmos           (:) = 0._dp 
      disturbance%nitrogen_2_sminn           (:) = 0._dp
    ENDIF

    ! Mark module as initialized
    module_initialized = .TRUE.

    CALL message('init_disturbance','disturbances initialised')

    IF (lrestart) RETURN

    ! -- Initialisation at beginning of an experiment
    disturbance%burned_frac     (:,:) = 0._dp
    disturbance%burned_frac_diag(:,:) = 0._dp
    disturbance%damaged_frac    (:,:) = 0._dp

  END SUBROUTINE init_disturbance

! ---- disturbed_fpc ----------------------------------------------------
!
! Driver for calculation of damaged surface area (presently by fire and wind breaks)
!
! History:
! 2011-04-29 - StW - Original code
  SUBROUTINE disturbed_frac(lctlib, nidx,kidx0,kidx1,ntiles,dist_type,pft_mask,with_yasso,act_fpc, veg_fract_correction, &
                            surface, cbalance, climbuf, dist_frac, dist_frac_diag, fuel)
  USE mo_jsbach_lctlib,    ONLY: lctlib_type 

  ! Input parameters
  TYPE(lctlib_type),INTENT(in) :: lctlib                   !! PFT-specific constants
  INTEGER                 , INTENT(in)    :: nidx          !! Number of points to calculate
  INTEGER                 , INTENT(in)    :: kidx0,kidx1   !! Start and ending of vector
  INTEGER                 , INTENT(in)    :: ntiles        !! Number of tiles
  INTEGER                 , INTENT(in)    :: dist_type     !! Sub-disturbance type (one of the DIST_ constants)
  LOGICAL                 , INTENT(in)    :: pft_mask(:)   !! True for pfts to which the disturbance should be applied
  LOGICAL                 , INTENT(in)    :: with_yasso
  REAL(dp)                , INTENT(in)    :: act_fpc(:,:)
  REAL(dp)                , INTENT(in)    :: veg_fract_correction(:,:)
  TYPE(land_surface_type) , INTENT(in)    :: surface
  TYPE(cbalance_type)     , INTENT(inout) :: cbalance
  TYPE(climbuf_type)      , INTENT(in)    :: climbuf
  ! Input/output parameters
  REAL(dp)                , INTENT(inout) :: dist_frac     (:,:) !! Fraction damaged by wind or fire, excluding survivors
  REAL(dp), OPTIONAL      , INTENT(inout) :: dist_frac_diag(:,:) !! Fraction damaged by wind or fire, including survivors
  ! Output parameter
  REAL(dp), OPTIONAL      , INTENT(inout) :: fuel(:)       !! Amount of fuel on which the fire was based 

  ! Locals
  REAL(dp) :: mortality(nidx,surface%ntiles)

    dist_frac(:,:) = 0._dp
    SELECT CASE (dist_type)
      CASE (DIST_FIRE)
        SELECT CASE (dist_opts%fire_algorithm)
          CASE (0) !! No fire algorithm
          CASE (1) !! jsbach algorithm
             CALL burned_frac_jsbach(nidx,ntiles,with_yasso,pft_mask(:),&
                         surface%is_present            (kidx0:kidx1,:), &
                         surface%cover_type            (kidx0:kidx1,:), &
                         surface%cover_fract           (kidx0:kidx1,:), &
                         veg_fract_correction                         , & 
                         surface%veg_ratio_max         (kidx0:kidx1  ), &
                         climbuf%rel_hum_air           (kidx0:kidx1  ), &
                         fuel, dist_frac,                               &
                         cbalance%cpool_litter_green_ag(kidx0:kidx1,:), &
                         cbalance%cpool_litter_wood_ag (kidx0:kidx1,:), &
                         cbalance%YCpool_acid_ag1      (kidx0:kidx1,:), &
                         cbalance%YCpool_water_ag1     (kidx0:kidx1,:), &
                         cbalance%YCpool_ethanol_ag1   (kidx0:kidx1,:), &
                         cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:), &
                         cbalance%YCpool_acid_ag2      (kidx0:kidx1,:), &
                         cbalance%YCpool_water_ag2     (kidx0:kidx1,:), &
                         cbalance%YCpool_ethanol_ag2   (kidx0:kidx1,:), &
                         cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:))
          CASE (2) !! Arora & Boer algorithm
          CASE (3) !! Thonicke algorithm
             CALL burned_frac_thonicke(lctlib,nidx,kidx0,kidx1,ntiles,with_yasso,pft_mask            (     :       ), &
                       surface%is_present            (kidx0:kidx1,:), surface%is_glacier             (kidx0:kidx1,:), &
                       veg_fract_correction,                                                                          &
                       surface%veg_ratio_max         (kidx0:kidx1  ), surface%cover_fract            (kidx0:kidx1,:), &
                       surface%cover_type            (kidx0:kidx1,:), cbalance%Cpool_litter_green_ag (kidx0:kidx1,:), &
                       cbalance%Cpool_litter_wood_ag (kidx0:kidx1,:), cbalance%Cpool_green           (kidx0:kidx1,:), &
                       cbalance%Cpool_woods          (kidx0:kidx1,:), cbalance%Cpool_reserve         (kidx0:kidx1,:), &
                       cbalance%frac_litter_wood_new  (kidx0:kidx1,:), &
                       climbuf%prev_day_mean_wind10  (kidx0:kidx1  ), climbuf%prev_day_mean_vol_moist(kidx0:kidx1  ), &
                       climbuf%prev_day_temp_max     (kidx0:kidx1  ), climbuf%prev_day_temp_min      (kidx0:kidx1  ), &
                       climbuf%prev_day_precip_mean  (kidx0:kidx1  ), & !climbuf%dew_point_temp         (kidx0:kidx1  ), &
                       cbalance%YCpool_acid_ag1      (kidx0:kidx1,:), cbalance%YCpool_acid_ag2       (kidx0:kidx1,:), &
                       cbalance%YCpool_water_ag1     (kidx0:kidx1,:), cbalance%YCpool_water_ag2      (kidx0:kidx1,:), &
                       cbalance%YCpool_ethanol_ag1   (kidx0:kidx1,:), cbalance%YCpool_ethanol_ag2    (kidx0:kidx1,:), &
                       cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:), cbalance%YCpool_nonsoluble_ag2 (kidx0:kidx1,:), &
                       fuel, dist_frac, mortality)
          CASE (4) !! Read burned_frac
          CASE DEFAULT
            CALL finish('disturbed_frac','Unknown fire algorithm')
        END SELECT
        ! Correct burned fpc for survivors dependent on fuel_consumption algorithm when dynveg = TRUE
        SELECT CASE (dist_opts%fuel_algorithm)
          CASE (0,1) ! No/jsbach fuel consumption. With and without survivors - data "copy" by initial pointer assignments
          CASE (2)   ! Arora & Boer fuel consumption - survivors considered 
          CASE (3)!! Thonicke spitfire
            dist_frac_diag(:,:) = dist_frac(:,:)
            dist_frac(:,:) = dist_frac(:,:) * mortality(:,:)
          CASE (4)   ! Arora & Boer fuel consumption - survivors considered 
        END SELECT

      CASE (DIST_WINDBREAK)
        SELECT CASE (dist_opts%windbreak_algorithm)
          CASE (0) !! No windbreak algorithm
          CASE (1) !! jsbach algorithm
            CALL broken_woody_frac_jsbach(nidx,ntiles,pft_mask,act_fpc,surface%cover_type(kidx0:kidx1,:), &
                                          climbuf%prev_day_max_wind10(kidx0:kidx1),climbuf%max_wind10(kidx0:kidx1),dist_frac)
          ! Insert other wind break algorithms here
          CASE DEFAULT
            CALL finish('disturbed_frac','Unknown windbreak algorithm')
        END SELECT
      CASE DEFAULT
        CALL finish('disturbed_fpc','Unknown disturbance type')
    END SELECT

  END SUBROUTINE disturbed_frac

! ---- relocate_disturbed_carbon ----------------------------------------------------
!
! Driver for calculation carbon relocation caused by disturbances
!
! History:
! 2011-04-29 - StW - Original code
  SUBROUTINE relocate_disturbed_carbon(lctlib,surface,nidx,kidx0,ntiles,with_yasso,dist_type,   &
                                       dist_frac, cover_fract, veg_fract_correction,            &
                                       Cpool_green, Cpool_woods, Cpool_reserve,                 &
                                       carbon_2_GreenLitterPools,carbon_2_WoodLitterPools,      &
                                       Cpool_litter_green_ag, Cpool_litter_green_bg,            &
                                       Cpool_litter_wood_ag,  Cpool_litter_wood_bg,             &
                                       YCpool_acid_ag1,       YCpool_acid_ag2,                  &
                                       YCpool_water_ag1,      YCpool_water_ag2,                 &
                                       YCpool_ethanol_ag1,    YCpool_ethanol_ag2,               &
                                       YCpool_nonsoluble_ag1, YCpool_nonsoluble_ag2,            &
                                       YCpool_acid_bg1,       YCpool_acid_bg2,                  &
                                       YCpool_water_bg1,      YCpool_water_bg2,                 &
                                       YCpool_ethanol_bg1,    YCpool_ethanol_bg2,               &
                                       YCpool_nonsoluble_bg1, YCpool_nonsoluble_bg2,            &
                                       YCpool_humus_1,        YCpool_humus_2,                   &
                                       LeafLit_coef,                                            &
                                       WoodLit_coef,                                            &
                                       Npool_green, Npool_woods,                                &
                                       Npool_mobile, SMINN_pool,                                &
                                       Npool_litter_green_ag, Npool_litter_green_bg,            &
                                       Npool_litter_wood_ag, Npool_litter_wood_bg,              &
                                       nitrogen_2_GreenLitterPools, nitrogen_2_WoodLitterPools, & 
                                       nitrogen_2_atmos, nitrogen_2_sminn, carbon_2_atmos)
  USE mo_jsbach_lctlib,    ONLY: lctlib_type 
  USE mo_land_surface,     ONLY: land_surface_type

  ! Input parameters
  TYPE (lctlib_type)      ,   INTENT(in) :: lctlib
  TYPE (land_surface_type),   INTENT(in) :: surface
  INTEGER                 ,   INTENT(in) :: nidx, kidx0 !! Number of points to calculate
  INTEGER                 ,   INTENT(in) :: ntiles      !! Number of tiles
  LOGICAL                 ,   INTENT(in) :: with_yasso
  INTEGER                 ,   INTENT(in) :: dist_type   !! Sub-disturbance type (one of the DIST_ constants)
  REAL(dp), DIMENSION(:,:),   INTENT(in) :: dist_frac 
  REAL(dp), DIMENSION(:,:),   INTENT(in) :: cover_fract !!
  REAL(dp), DIMENSION(:,:),   INTENT(in) :: veg_fract_correction
  ! Input/output parameters
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: Cpool_green 
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: Cpool_woods 
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: Cpool_reserve
  REAL(dp), DIMENSION(  :),   INTENT(inout)           :: carbon_2_GreenLitterPools 
  REAL(dp), DIMENSION(  :),   INTENT(inout)           :: carbon_2_WoodLitterPools 
  REAL(dp), DIMENSION(  :),   INTENT(inout), OPTIONAL :: carbon_2_atmos 
  ! JSBACH litter pools
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: Cpool_litter_green_ag
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: Cpool_litter_green_bg
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: Cpool_litter_wood_ag
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: Cpool_litter_wood_bg
  ! Yasso litter pools
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: YCpool_acid_ag1, YCpool_acid_ag2
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: YCpool_acid_bg1, YCpool_acid_bg2
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: YCpool_water_ag1, YCpool_water_ag2
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: YCpool_water_bg1, YCpool_water_bg2
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: YCpool_ethanol_ag1,YCpool_ethanol_ag2
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: YCpool_ethanol_bg1, YCpool_ethanol_bg2
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: YCpool_nonsoluble_ag1,YCpool_nonsoluble_ag2
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: YCpool_nonsoluble_bg1,YCpool_nonsoluble_bg2
  REAL(dp), DIMENSION(:,:),   INTENT(inout)           :: YCpool_humus_1, YCpool_humus_2
  REAL(dp), DIMENSION(:,:,:), INTENT(in)              :: LeafLit_coef 
  REAL(dp), DIMENSION(:,:,:), INTENT(in)              :: WoodLit_coef 
  ! Nitrogen variables
  REAL(dp), DIMENSION(:,:),   INTENT(inout), OPTIONAL :: Npool_green 
  REAL(dp), DIMENSION(:,:),   INTENT(inout), OPTIONAL :: Npool_woods 
  REAL(dp), DIMENSION(:,:),   INTENT(inout), OPTIONAL :: Npool_mobile
  REAL(dp), DIMENSION(:,:),   INTENT(inout), OPTIONAL :: SMINN_pool     
  REAL(dp), DIMENSION(:,:),   INTENT(inout), OPTIONAL :: Npool_litter_green_ag
  REAL(dp), DIMENSION(:,:),   INTENT(inout), OPTIONAL :: Npool_litter_green_bg
  REAL(dp), DIMENSION(:,:),   INTENT(inout), OPTIONAL :: Npool_litter_wood_ag
  REAL(dp), DIMENSION(:,:),   INTENT(inout), OPTIONAL :: Npool_litter_wood_bg
  REAL(dp), DIMENSION(  :),   INTENT(inout), OPTIONAL :: nitrogen_2_GreenLitterPools 
  REAL(dp), DIMENSION(  :),   INTENT(inout), OPTIONAL :: nitrogen_2_WoodLitterPools 
  REAL(dp), DIMENSION(  :),   INTENT(inout), OPTIONAL :: nitrogen_2_atmos   
  REAL(dp), DIMENSION(  :),   INTENT(inout), OPTIONAL :: nitrogen_2_sminn           

    IF (IAND(dist_type,DIST_FIRE) /= 0) THEN 
      SELECT CASE (dist_opts%fuel_algorithm)
        CASE (0) !! No fuel consumption
        CASE (1) !!  jsbach fuel consumption
          IF (.NOT. PRESENT(carbon_2_atmos)) CALL finish('relocate_disturbed_carbon', &
            'Input parameter carbon_2_atmos is needed for relocation by fire')
          CALL relocate_carbon_fire(nidx,ntiles,with_yasso,                               &
                               dist_frac(:,:),                                            &
                               cover_fract(:,:), veg_fract_correction(:,:),               &
                               dist_opts%fire_frac_wood_2_atmos,Cpool_green(:,:),         &
                               Cpool_woods(:,:), Cpool_reserve(:,:),                      &
                               carbon_2_GreenLitterPools(:),carbon_2_WoodLitterPools(:),  &
                               carbon_2_atmos       (:),                                  &
                               Cpool_litter_green_ag(:,:), Cpool_litter_green_bg(:,:),    &
                               Cpool_litter_wood_ag (:,:), Cpool_litter_wood_bg (:,:),    &
                               YCpool_acid_ag1      (:,:), YCpool_acid_ag2      (:,:),    &
                               YCpool_water_ag1     (:,:), YCpool_water_ag2     (:,:),    &
                               YCpool_ethanol_ag1   (:,:), YCpool_ethanol_ag2   (:,:),    &
                               YCpool_nonsoluble_ag1(:,:), YCpool_nonsoluble_ag2(:,:),    &
                               YCpool_acid_bg1      (:,:), YCpool_acid_bg2      (:,:),    &
                               YCpool_water_bg1     (:,:), YCpool_water_bg2     (:,:),    &
                               YCpool_ethanol_bg1   (:,:), YCpool_ethanol_bg2   (:,:),    &
                               YCpool_nonsoluble_bg1(:,:), YCpool_nonsoluble_bg2(:,:),    &
                               YCpool_humus_1       (:,:), YCpool_humus_2       (:,:),    &
                               LeafLit_coef       (:,:,:), WoodLit_coef       (:,:,:),    &
                               ! Nitrogen parameters below may not be present and can thus not take indices through the interface
                               Npool_green                 = Npool_green,                 &
                               Npool_woods                 = Npool_woods,                 &
                               Npool_mobile                = Npool_mobile,                &
                               SMINN_pool                  = SMINN_pool,                  &
                               Npool_litter_green_ag       = Npool_litter_green_ag,       &
                               Npool_litter_green_bg       = Npool_litter_green_bg,       &
                               Npool_litter_wood_ag        = Npool_litter_wood_ag,        &
                               Npool_litter_wood_bg        = Npool_litter_wood_bg,        &
                               nitrogen_2_GreenLitterPools = nitrogen_2_GreenLitterPools, &
                               nitrogen_2_WoodLitterPools  = nitrogen_2_WoodLitterPools,  &
                               nitrogen_2_atmos            = nitrogen_2_atmos,            &
                               nitrogen_2_sminn            = nitrogen_2_sminn)
        CASE (3)
          !! relocation of carbon according to the spitfire algorithm
          CALL relocate_carbon_fire_thonicke(lctlib,surface, nidx,kidx0,ntiles,with_yasso,        &
                                             lctlib%dynamic_pft(:) .OR. (dist_opts%lburn_pasture .AND. lctlib%pasture_pft(:)), &
                                             dist_frac(:,:),                                      &
                                             cover_fract(:,:), veg_fract_correction,              &
                                             Cpool_green, Cpool_woods,Cpool_reserve,              &
                                             Cpool_litter_green_ag,Cpool_litter_green_bg,         &
                                             Cpool_litter_wood_ag, Cpool_litter_wood_bg,          &
                                             YCpool_acid_ag1,      YCpool_acid_ag2,               &
                                             YCpool_water_ag1,     YCpool_water_ag2,              &
                                             YCpool_ethanol_ag1,   YCpool_ethanol_ag2,            &
                                             YCpool_nonsoluble_ag1,YCpool_nonsoluble_ag2,         &
                                             YCpool_acid_bg1,      YCpool_acid_bg2,               &
                                             YCpool_water_bg1,     YCpool_water_bg2,              &
                                             YCpool_ethanol_bg1,   YCpool_ethanol_bg2,            &
                                             YCpool_nonsoluble_bg1,YCpool_nonsoluble_bg2,         &
                                             YCpool_humus_1,       YCpool_humus_2,                &
                                             LeafLit_coef,         WoodLit_coef,                  &
                                             carbon_2_GreenLitterPools,carbon_2_WoodLitterPools,  &
                                             carbon_2_atmos)
        CASE (4) !! Read fuel consumption
        CASE DEFAULT
          CALL finish('relocate_disturbed_carbon','Unknown fire algorithm')
      END SELECT
    ELSEIF (IAND(dist_type,DIST_WINDBREAK) /= 0) THEN
      SELECT CASE (dist_opts%windbreak_algorithm)
        CASE (0) !! No windbreak algorithm
        CASE (1) !! jsbach windbreak algorithm
           CALL relocate_carbon_damage(nidx, ntiles, with_yasso, dist_frac(:,:),                  &
                                       cover_fract(:,:), veg_fract_correction(:,:),               &
                                       Cpool_green(:,:), Cpool_woods(:,:), Cpool_reserve(:,:),    &
                                       carbon_2_GreenLitterPools(:), carbon_2_WoodLitterPools(:), &
                                       Cpool_litter_green_ag(:,:), Cpool_litter_green_bg(:,:),    &
                                       Cpool_litter_wood_ag(:,:), Cpool_litter_wood_bg(:,:),      &
                                       YCpool_acid_ag1,       YCpool_acid_ag2,                    &
                                       YCpool_water_ag1,      YCpool_water_ag2,                   &
                                       YCpool_ethanol_ag1,    YCpool_ethanol_ag2,                 &
                                       YCpool_nonsoluble_ag1, YCpool_nonsoluble_ag2,              &
                                       YCpool_acid_bg1,       YCpool_acid_bg2,                    &
                                       YCpool_water_bg1,      YCpool_water_bg2,                   &
                                       YCpool_ethanol_bg1,    YCpool_ethanol_bg2,                 &
                                       YCpool_nonsoluble_bg1, YCpool_nonsoluble_bg2,              &
                                       YCpool_humus_1,        YCpool_humus_2,                     &
                                       LeafLit_coef,          WoodLit_coef,                       &
                                       Npool_green                 =  Npool_green,                &
                                       Npool_woods                 =  Npool_woods,                &
                                       Npool_mobile                =  Npool_mobile,               &
                                       SMINN_pool                  =  SMINN_pool,                 &
                                       Npool_litter_green_ag       =  Npool_litter_green_ag,      &
                                       Npool_litter_green_bg       =  Npool_litter_green_bg,      &
                                       Npool_litter_wood_ag        =  Npool_litter_wood_ag,       &
                                       Npool_litter_wood_bg        =  Npool_litter_wood_bg,       &
                                       nitrogen_2_GreenLitterPools = nitrogen_2_GreenLitterPools, &
                                       nitrogen_2_WoodLitterPools  = nitrogen_2_WoodLitterPools,  &
                                       nitrogen_2_sminn            = nitrogen_2_sminn)
        ! Insert other windbreak algorithms here
        CASE DEFAULT
          CALL finish('relocate_disturbed_carbon','Unknown windbreak algorithm')
      END SELECT
    ENDIF

  END SUBROUTINE relocate_disturbed_carbon

! ---- calc_disturbance -----------------------------------------------
!
! Calculate disturbances (wind breaks and fires) and their associated
! carbon tranfers.
! Should exclusively be called, wenn dynveg is NOT in use.
!
! HISTORY:
! 2011-05-03 - StW - Original code extracted from mo_dynveg

  SUBROUTINE update_disturbance(nidx, kidx0, kidx1, climbuf, surface, cbalance, with_yasso, with_nitrogen, nbalance, &
                                lctlib, veg_fract_correction, CO2_emission, LeafLit_coef, WoodLit_coef)

  USE mo_jsbach_constants,    ONLY: molarMassCO2_kg
  USE mo_jsbach,              ONLY: new_day

  INTEGER                 , INTENT(IN) :: nidx,kidx0,kidx1
  TYPE(land_surface_type) , INTENT(IN) :: surface
  TYPE(climbuf_type)      , INTENT(IN) :: climbuf
  TYPE(cbalance_type)     , INTENT(INOUT) :: cbalance 
  LOGICAL                 , INTENT(IN) :: with_yasso
  LOGICAL                 , INTENT(IN) :: with_nitrogen
  TYPE(nbalance_type)     , INTENT(INOUT) :: nbalance
  TYPE(lctlib_type)       , INTENT(IN) :: lctlib
  REAL(dp), DIMENSION(:,:), INTENT(IN) :: veg_fract_correction

  REAL(dp), DIMENSION(:)  , INTENT(OUT) :: CO2_emission
  REAL(dp),                 INTENT(IN)  :: LeafLit_coef(:,:,:)
  REAL(dp),                 INTENT(IN)  :: WoodLit_coef(:,:,:)

    IF (new_day) THEN
      disturbance%burned_frac (kidx0:kidx1,:) = 0._dp
      disturbance%damaged_frac(kidx0:kidx1,:) = 0._dp

      ! Calculate burned and damaged fractions
      CALL disturbed_frac(lctlib,nidx,kidx0,kidx1,surface%ntiles,DIST_FIRE,                                 &
                          lctlib%dynamic_pft(:) .OR. (dist_opts%lburn_pasture .AND. lctlib%pasture_pft(:)), &
                          with_yasso,surface%cover_fract_pot(kidx0:kidx1,:),veg_fract_correction(1:nidx,:), &
                          surface,cbalance,climbuf,disturbance%burned_frac(kidx0:kidx1,:),                  &
                          disturbance%burned_frac_diag(kidx0:kidx1,:), disturbance%fuel(kidx0:kidx1))

      CALL disturbed_frac(lctlib,nidx,kidx0,kidx1,surface%ntiles,DIST_WINDBREAK,                 &
                          lctlib%woody_pft(:) .AND. lctlib%dynamic_pft(:), with_yasso,           &
                          surface%cover_fract_pot(kidx0:kidx1,:),veg_fract_correction(1:nidx,:), &
                          surface,cbalance,climbuf,disturbance%damaged_frac(kidx0:kidx1,:))

      ! Cumulate disturbed areas for diagnostics
      disturbance%box_burned_frac_avg     (kidx0:kidx1,:) =                                                    & 
      disturbance%box_burned_frac_avg     (kidx0:kidx1,:) + disturbance%burned_frac       (kidx0:kidx1,:)      &
                                                            * surface%cover_fract         (kidx0:kidx1,:)      &
                                                            * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2, &
                                                                     ncopies=surface%ntiles)                   &
                                                            * 86400._dp
      IF (.NOT. ASSOCIATED(disturbance%box_burned_frac_diag_avg,disturbance%box_burned_frac_avg))                &
        disturbance%box_burned_frac_diag_avg(kidx0:kidx1,:) =                                                    &
        disturbance%box_burned_frac_diag_avg(kidx0:kidx1,:) + disturbance%burned_frac_diag  (kidx0:kidx1,:)      &
                                                              * surface%cover_fract         (kidx0:kidx1,:)      &
                                                              * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2, &
                                                                       ncopies=surface%ntiles)                   &
                                                              * 86400._dp
      disturbance%box_damaged_frac_avg    (kidx0:kidx1,:) =                                                    &
      disturbance%box_damaged_frac_avg    (kidx0:kidx1,:) + disturbance%damaged_frac      (kidx0:kidx1,:)      &
                                                            * surface%cover_fract         (kidx0:kidx1,:)      &
                                                            * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2, &
                                                                     ncopies=surface%ntiles)                   &
                                                            * 86400._dp

      ! Calculate carbon relocation connected with damages and fires
      IF (with_nitrogen) THEN
        CALL relocate_disturbed_carbon(lctlib,surface,nidx,kidx0,surface%ntiles,with_yasso,    &
                                       DIST_WINDBREAK,                                         &
                                       disturbance%damaged_frac(kidx0:kidx1,:),                &
                                       surface%cover_fract(kidx0:kidx1,:),                     &
                                       veg_fract_correction(1:nidx,:),                         &
                                       cbalance%cpool_green(kidx0:kidx1,:),                    &
                                       cbalance%cpool_woods(kidx0:kidx1,:),                    &
                                       cbalance%cpool_reserve(kidx0:kidx1,:),                  &
                                       disturbance%carbon_2_GreenLitterPools(kidx0:kidx1),     &
                                       disturbance%carbon_2_WoodLitterPools(kidx0:kidx1),      &
                                       cbalance%cpool_litter_green_ag(kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_green_bg(kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_wood_ag(kidx0:kidx1,:),           &
                                       cbalance%cpool_litter_wood_bg(kidx0:kidx1,:),           &
                                       cbalance%YCpool_acid_ag1      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_ag2      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_ag1     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_ag2     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_ag1   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_ag2   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_bg1      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_bg2      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_bg1     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_bg2     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_bg1   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_bg2   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:),          &
                                       cbalance%YCpool_humus_1       (kidx0:kidx1,:),          &
                                       cbalance%YCpool_humus_2       (kidx0:kidx1,:),          &
                                       LeafLit_coef                  (1:nidx,:,:),             &
                                       WoodLit_coef                  (1:nidx,:,:),             &
             npool_green            =  nbalance%npool_green(kidx0:kidx1,:),                    &
             npool_woods            =  nbalance%npool_woods(kidx0:kidx1,:),                    &
             npool_mobile           =  nbalance%npool_mobile(kidx0:kidx1,:),                   &
             sminn_pool             =  nbalance%sminn_pool(kidx0:kidx1,:),                     &
             npool_litter_green_ag  =  nbalance%npool_litter_green_ag(kidx0:kidx1,:),          &
             npool_litter_green_bg  =  nbalance%npool_litter_green_bg(kidx0:kidx1,:),          &
             npool_litter_wood_ag   =  nbalance%npool_litter_wood_ag(kidx0:kidx1,:),           &
             npool_litter_wood_bg   =  nbalance%npool_litter_wood_bg(kidx0:kidx1,:),           &
             nitrogen_2_GreenLitterPools=disturbance%nitrogen_2_GreenLitterPools(kidx0:kidx1), &
             nitrogen_2_WoodLitterPools=disturbance%nitrogen_2_WoodLitterPools(kidx0:kidx1))

        CALL relocate_disturbed_carbon(lctlib,surface,nidx,kidx0,surface%ntiles,with_yasso,    &
                                       DIST_FIRE,                                              &
                                       disturbance%burned_frac(kidx0:kidx1,:),                 &
                                       surface%cover_fract(kidx0:kidx1,:),                     &
                                       veg_fract_correction(1:nidx,:),                         &
                                       cbalance%cpool_green(kidx0:kidx1,:),                    &
                                       cbalance%cpool_woods(kidx0:kidx1,:),                    &
                                       cbalance%cpool_reserve(kidx0:kidx1,:),                  &
                                       disturbance%carbon_2_GreenLitterPools(kidx0:kidx1),     &
                                       disturbance%carbon_2_WoodLitterPools(kidx0:kidx1),      &
                                       cbalance%cpool_litter_green_ag(kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_green_bg(kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_wood_ag (kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_wood_bg (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_ag1      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_ag2      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_ag1     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_ag2     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_ag1   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_ag2   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_bg1      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_bg2      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_bg1     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_bg2     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_bg1   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_bg2   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:),          &
                                       cbalance%YCpool_humus_1       (kidx0:kidx1,:),          &
                                       cbalance%YCpool_humus_2       (kidx0:kidx1,:),          &
                                       LeafLit_coef                  (1:nidx,:,:),             &
                                       WoodLit_coef                  (1:nidx,:,:),             &
             carbon_2_atmos         =  disturbance%carbon_2_atmos(kidx0:kidx1),                &
             npool_green            =  nbalance%npool_green(kidx0:kidx1,:),                    &
             npool_woods            =  nbalance%npool_woods(kidx0:kidx1,:),                    &
             npool_mobile           =  nbalance%npool_mobile(kidx0:kidx1,:),                   &
             sminn_pool             =  nbalance%sminn_pool(kidx0:kidx1,:),                     &
             npool_litter_green_ag  =  nbalance%npool_litter_green_ag(kidx0:kidx1,:),          &
             npool_litter_green_bg  =  nbalance%npool_litter_green_bg(kidx0:kidx1,:),          &
             npool_litter_wood_ag   =  nbalance%npool_litter_wood_ag(kidx0:kidx1,:),           &
             npool_litter_wood_bg   =  nbalance%npool_litter_wood_bg(kidx0:kidx1,:),           &
             nitrogen_2_GreenLitterPools=disturbance%nitrogen_2_GreenLitterPools(kidx0:kidx1), &
             nitrogen_2_WoodLitterPools=disturbance%nitrogen_2_WoodLitterPools(kidx0:kidx1),   &
             nitrogen_2_atmos       =   disturbance%nitrogen_2_atmos(kidx0:kidx1))
      ELSE
        CALL relocate_disturbed_carbon(lctlib,surface,nidx,kidx0,surface%ntiles,with_yasso,    &
                                       DIST_WINDBREAK,                                         &
                                       disturbance%damaged_frac(kidx0:kidx1,:),                &
                                       surface%cover_fract(kidx0:kidx1,:),                     &
                                       veg_fract_correction(1:nidx,:),                         &
                                       cbalance%cpool_green(kidx0:kidx1,:),                    &
                                       cbalance%cpool_woods(kidx0:kidx1,:),                    &
                                       cbalance%cpool_reserve(kidx0:kidx1,:),                  &
                                       disturbance%carbon_2_GreenLitterPools(kidx0:kidx1),     &
                                       disturbance%carbon_2_WoodLitterPools(kidx0:kidx1),      &
                                       cbalance%cpool_litter_green_ag(kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_green_bg(kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_wood_ag (kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_wood_bg (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_ag1      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_ag2      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_ag1     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_ag2     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_ag1   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_ag2   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_bg1      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_bg2      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_bg1     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_bg2     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_bg1   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_bg2   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:),          &
                                       cbalance%YCpool_humus_1       (kidx0:kidx1,:),          &
                                       cbalance%YCpool_humus_2       (kidx0:kidx1,:),          &
                                       LeafLit_coef                  (1:nidx,:,:),             &
                                       WoodLit_coef                  (1:nidx,:,:))

        CALL relocate_disturbed_carbon(lctlib,surface,nidx,kidx0,surface%ntiles,with_yasso,    &
                                       DIST_FIRE,                                              &
                                       disturbance%burned_frac(kidx0:kidx1,:),                 &
                                       surface%cover_fract(kidx0:kidx1,:),                     &
                                       veg_fract_correction(1:nidx,:),                         &
                                       cbalance%cpool_green(kidx0:kidx1,:),                    &
                                       cbalance%cpool_woods(kidx0:kidx1,:),                    &
                                       cbalance%cpool_reserve(kidx0:kidx1,:),                  &
                                       disturbance%carbon_2_GreenLitterPools(kidx0:kidx1),     &
                                       disturbance%carbon_2_WoodLitterPools(kidx0:kidx1),      &
                                       cbalance%cpool_litter_green_ag(kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_green_bg(kidx0:kidx1,:),          &
                                       cbalance%cpool_litter_wood_ag(kidx0:kidx1,:),           &
                                       cbalance%cpool_litter_wood_bg(kidx0:kidx1,:),           &
                                       cbalance%YCpool_acid_ag1      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_ag2      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_ag1     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_ag2     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_ag1   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_ag2   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_bg1      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_acid_bg2      (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_bg1     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_water_bg2     (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_bg1   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_ethanol_bg2   (kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:),          &
                                       cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:),          &
                                       cbalance%YCpool_humus_1       (kidx0:kidx1,:),          &
                                       cbalance%YCpool_humus_2       (kidx0:kidx1,:),          &
                                       LeafLit_coef                  (1:nidx,:,:),             &
                                       WoodLit_coef                  (1:nidx,:,:),             &
                      carbon_2_atmos = disturbance%carbon_2_atmos    (kidx0:kidx1))
      END IF

      ! Scale fluxes from vegetated area to whole grid box
      disturbance%carbon_2_GreenLitterPools(kidx0:kidx1) = disturbance%carbon_2_GreenLitterPools(kidx0:kidx1) * &
                                                                   surface%veg_ratio_max(kidx0:kidx1)
      disturbance%carbon_2_WoodLitterPools(kidx0:kidx1)  = disturbance%carbon_2_WoodLitterPools(kidx0:kidx1)  * &
                                                                   surface%veg_ratio_max(kidx0:kidx1)
      disturbance%carbon_2_atmos(kidx0:kidx1) =  disturbance%carbon_2_atmos(kidx0:kidx1) * surface%veg_ratio_max(kidx0:kidx1)

      ! Cumulate CO2 emissions for flux output
      disturbance%box_CO2_flux_2_atmos(kidx0:kidx1) = disturbance%box_CO2_flux_2_atmos(kidx0:kidx1) + &
                                                      disturbance%carbon_2_atmos(kidx0:kidx1) * molarMassCO2_kg
    ENDIF

    ! Calculate CO2 flux to the atmosphere (updated once a day)
    CO2_emission(:) = disturbance%carbon_2_atmos(kidx0:kidx1) * molarMassCO2_kg/86400._dp

  END SUBROUTINE update_disturbance

END MODULE mo_disturbance
