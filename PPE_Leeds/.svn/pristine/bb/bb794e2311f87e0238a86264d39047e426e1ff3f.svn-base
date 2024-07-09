!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_dynveg
!
! Computation of pft fractions per gridbox - scheme adopted from lpj
!

  USE mo_jsbach_grid,      ONLY: grid_type, domain_type
  USE mo_kind,             ONLY: dp
  USE mo_exception,        ONLY: finish, message, message_text
  USE mo_linked_list,      ONLY: t_stream
  USE mo_mpi,              ONLY: p_parallel_io, p_io, p_parallel, p_bcast
  USE mo_netCDF,           ONLY: file_info, nf_max_name
  USE mo_io,               ONLY: IO_READ
  USE mo_land_surface,     ONLY: scale_cover_fract, fract_small, land_surface_type
  USE mo_jsbach,           ONLY: debug_Cconservation
  USE mo_disturbance,      ONLY: disturbance, disturbed_frac, DIST_FIRE, DIST_WINDBREAK,relocate_disturbed_carbon
  USE mo_cbal_bethy,       ONLY: cbalance_type, nbalance_type
  USE mo_climbuf,          ONLY: climbuf_type

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART ======================================================================================================

  TYPE dynveg_type  !! contains the state variables
     REAL(dp),POINTER :: bio_exist(:,:) 
     REAL(dp),POINTER :: litter_moisture(:)
     REAL(dp),POINTER :: act_fpc(:,:)
     REAL(dp),POINTER :: pot_fpc(:,:) 
     REAL(dp),POINTER :: bare_fpc(:)
     REAL(dp),POINTER :: desert_fpc(:)
     REAL(dp),POINTER :: dynveg_testCconserv_1(:)
     REAL(dp),POINTER :: dynveg_testCconserv_2(:)
     REAL(dp),POINTER :: max_green_bio(:,:)      !! maximum value of green biomass within a year
     REAL(dp),POINTER :: sum_green_bio_memory(:) !! vegetated fraction calculated from green biomass
     REAL(dp),POINTER :: Dummy_null(:,:)         !! Permanently filled with 0 for replacement in fpc_daily if no disturbance
  END TYPE dynveg_type
  TYPE(dynveg_type), SAVE :: dynveg
  PUBLIC :: dynveg

  TYPE dynveg_options_type
     LOGICAL                :: dynveg_all         !! competition between all PFTs, not just between woody PFTs
     LOGICAL                :: dynveg_feedback    !! activate climate feedback of the dynamic vegetation
     LOGICAL                :: read_fpc           !! Read fractional plant cover from file
     CHARACTER(nf_max_name) :: fpc_file_name
     REAL(dp)               :: accelerate_dynveg  !! factor to accelerate vegetation and desert dynamics
  END TYPE dynveg_options_type
  TYPE(dynveg_options_type), SAVE :: dynveg_options
  PUBLIC :: dynveg_options

  PUBLIC :: dynveg_type
  PUBLIC :: dynveg_options_type
  PUBLIC :: update_dynveg 
  PUBLIC :: init_dynveg
  PUBLIC :: config_dynveg
  PUBLIC :: fpc_daily
  PUBLIC :: calc_veg_ratio_max
  PUBLIC :: fpc_to_cover_fract_pot
  PUBLIC :: cover_fract_pot_to_cover_fract
  PUBLIC :: scale_fpc

!  PUBLIC :: frac_wood_2_atmos_fire, 
  PUBLIC :: act_fpc_min

  REAL(dp), PARAMETER :: act_fpc_min                = 0.005_dp   !! minimum actual FPC to be considered to calculate potential FPC

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ==============================================================================
  PRIVATE

  REAL(dp), PARAMETER :: sum_npp_min           = 1.e-12_dp  !! minimum total npp for tree cover, kg C/m2
  REAL(dp), PARAMETER :: tree_fpc_max          = 1.0_dp     !! maximum FPC of woody PFTs
  REAL(dp), PARAMETER :: grass_mortality       = 0.01_dp / 365._dp !! grass mortality (not caused by fire, but by other unresolved 
                                                                   !! processes)
  REAL(dp), PARAMETER :: npp_nonlinearity      = 1.5_dp     !! parameter controlling the non-linearity in the dynamic equation with
                                                            !!    respect to NPP
  REAL(dp), PARAMETER :: tau_Cpool_woods_to_tau_pft = 1._dp !! factor to translate life-time of wood in time scale of establishment
                                                            !!     and mortality
  REAL(dp), PARAMETER :: desert_extend         = 0.65_dp    !! parameter controlling the extend of desert (the lower the value the
                                                            !!    more desert)
  REAL(dp), PARAMETER :: desert_margin         = 2.0_dp     !! parameter controlling the transition from vegetated land to desert
                                                            !! (the higher the value the sharper the transition)
  REAL(dp), PARAMETER :: tau_desert            = 50._dp     !! time constant by which the extension of deserts is adapted
  REAL(dp), PARAMETER :: life_to_estab         = 20._dp     !! ratio of life-time of plants to their establishment time scale
                                                            !! (only used, when competiton of woody types and grass is active)

  LOGICAL, SAVE      :: module_initialized = .false. !! Signifies whether module has been initialized

  TYPE(t_stream),  POINTER :: IO_dynveg    !! Memory stream for dynveg model state

  TYPE(file_info),   SAVE  :: fpc_file               !! Input file for FPC
  REAL(dp), POINTER, SAVE  :: init_act_fpc(:,:)      !! initial values for act_fpc if read from file
  REAL(dp), POINTER, SAVE  :: init_pot_fpc(:,:)      !! initial values for pot_fpc if read from file
  REAL(dp), POINTER, SAVE  :: init_bare_fpc(:)       !! initial values for bare_fpc if read from file
  REAL(dp), POINTER, SAVE  :: init_desert_fpc(:)     !! initial values for desert_fpc if read from file
  REAL(dp), POINTER, SAVE  :: init_max_green_bio(:,:)      !! initial values for max_green_bio if read from file
  REAL(dp), POINTER, SAVE  :: init_sum_green_bio_memory(:) !! initial values for sum_green_bio_memory if read from file
  REAL(dp), POINTER, SAVE  :: init_cover_fract(:,:)        !! initial values for cover_fract if read from file
  REAL(dp), POINTER, SAVE  :: init_cover_fract_pot(:,:)    !! initial values for cover_fract_pot if read from file

! !VERSION CONTROL:
  CHARACTER(len=*), PARAMETER :: version = '$ID$'


! !GLOBAL VARIABLES:

! initialization of bioclimatic limits from LPJ lookup table
! initial values for FPC: from original vegetation map of JSBACH
 
!      6  flammability threshold
!      8  fire resistance nindex
!
!     BIOCLIMATIC LIMITS
!
!     28 minimum coldest monthly mean temperature
!     29 maximum coldest monthly mean temperature
!     30 minimum growing degree days (at or above 5 deg C)
!     31 upper limit of temperature of the warmest month 
!     32 lower limit of growth efficiency (g/m2)
!     
!     PARAMETERS ADDED LATER
!
!     33 GDD base
!     34 20-year average min warmest - coldest month temperature range
!     10 wind damage threshold (m/s)

! Plant Functional Types (PFT) currently regarded in this scheme (and correspondence to LPJ PFTs):

! JSBACH: 1. tropical broadleaved evergreen forest       - LPJ type 1
! JSBACH: 2. tropical deciduous broadleaved forest       - LPJ type 2
! JSBACH: 3. temper./boreal evergreen forest             - LPJ types 4 (extended cold tolerance), 3 and 6
! JSBACH: 4. temper./boreal deciduous forest             - LPJ types 5, 7  and 8
! JSBACH: 5. raingreen shrubs                            - LPJ types 1, 3, 4
! JSBACH: 6. cold shrubs                                 - LPJ type 8 (tundra)
! JSBACH: 7. C3 perennial grass	                         - LPJ type 9
! JSBACH: 8. C4 perennial grass                          - LPJ type 10 
! JSBACH: 9. crops
! JSBACH: 10. pasture

CONTAINS

  ! --- config_dynveg ------------------------------------------------------------------
  !
  ! Reads in parameters from lctlib-file

  SUBROUTINE config_dynveg()

    USE mo_namelist,         ONLY: position_nml, POSITIONED
    USE mo_jsbach,           ONLY: nml_unit
    USE mo_io_units,         ONLY: nout

    !! --- locals
    INTEGER                :: read_status, f_unit

    !! Namelist Parameters
    CHARACTER(NF_MAX_NAME) :: fpc_file_name      !! file name for initial PFTs
    LOGICAL                :: dynveg_all         !! competition between all PFTs, not just between woody PFTs
    LOGICAL                :: dynveg_feedback    !! activate climate feedback of the dynamic vegetation
    LOGICAL                :: read_fpc           !! Read fractional plant cover from file
    REAL(dp)               :: accelerate_dynveg  !! factor by which vegetation and desert dynamics are accelerated, 
                                                 !!    should be 1 for production runs

    INCLUDE 'dynveg_ctl.inc'


    !! Read namelist dynveg_ctl
    IF (p_parallel_io) THEN

       ! define default values
       read_fpc = .FALSE. 
       fpc_file_name = 'fpc.nc'
       dynveg_all = .FALSE.
       dynveg_feedback = .TRUE.
       accelerate_dynveg = 1._dp
       
       f_unit = position_nml ('DYNVEG_CTL', nml_unit, status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (f_unit, dynveg_ctl)
          CALL message('config_dynveg', 'Namelist DYNVEG_CTL: ')
          WRITE(nout, dynveg_ctl)
       END SELECT
    ENDIF
    IF (p_parallel_io) THEN

       dynveg_options%read_fpc = read_fpc
       dynveg_options%fpc_file_name = fpc_file_name
       fpc_file%file_name = fpc_file_name

       IF (read_fpc) THEN
          CALL message('config_dynveg','Fractional plant cover read from external file')
          CALL message('config_dynveg','File containing initial values for fpc: '//TRIM(fpc_file_name))
       ENDIF

       dynveg_options%dynveg_all = dynveg_all
       IF (dynveg_all) THEN
          CALL message('config_dynveg', 'Competition between woody types and grass active')
       ENDIF

       dynveg_options%dynveg_feedback = dynveg_feedback
       IF (dynveg_feedback) THEN
          CALL message('config_dynveg', 'Dynamic vegetation feedback active')
       ENDIF

       dynveg_options%accelerate_dynveg = accelerate_dynveg
       IF (dynveg_options%accelerate_dynveg <= 0._dp) THEN
          CALL finish('config_dynveg', 'DYNVEG_CTL parameter accelerate_dynveg has to be > 0.')
       ENDIF
       IF (dynveg_options%accelerate_dynveg /= 1._dp) THEN
          WRITE (message_text, *) 'Run with accelerated vegetation dynamics. Factor: ', accelerate_dynveg
          CALL message('config_dynveg', message_text)
       ENDIF

    ENDIF
    IF (p_parallel) THEN
       CALL p_bcast(dynveg_options%read_fpc, p_io)
       CALL p_bcast(dynveg_options%dynveg_all, p_io)
       CALL p_bcast(dynveg_options%dynveg_feedback, p_io)
       CALL p_bcast(dynveg_options%accelerate_dynveg, p_io)
    ENDIF

  END SUBROUTINE config_dynveg

  ! --- init_dynveg --------------------------------------------------------------------------------------------------------
  !
  ! Initialises this dynveg module. In particular the structure "dynveg" is initialised,
  ! this means, that the initialisation of dynveg is done here.
  
  SUBROUTINE init_dynveg(run_disturbance,grid, domain, ntiles, lrestart, fileformat, fileztype, stream)

    USE mo_linked_list,            ONLY: LAND, TILES
    USE mo_memory_base,            ONLY: new_stream,default_stream_setting, &
                                         add =>add_stream_element
    USE mo_output,                 ONLY: veg_table
    USE mo_netcdf,                 ONLY: max_dim_name, io_inq_varid, io_get_var_double
    USE mo_io,                     ONLY: io_open, io_close
    USE mo_tr_scatter,             ONLY: scatter_field
    USE mo_temp,                   ONLY: zreal2d, zreal3d, zzreal2d, zreal2d_ptr
    USE mo_jsbach,                 ONLY: missing_value

    LOGICAL,           INTENT(in)        :: run_disturbance
    TYPE(grid_type),   INTENT(in)        :: grid
    TYPE(domain_type), INTENT(in)        :: domain
    LOGICAL,           INTENT(in)        :: lrestart      ! true for restarted runs (throughout the whole run)
    INTEGER,           INTENT(in)        :: ntiles        ! maximum number of living pfts in a grid cell
    INTEGER,           INTENT(in)        :: fileformat    ! output file format
    INTEGER,           INTENT(in)        :: fileztype     ! output file compression
    TYPE(t_stream), POINTER, OPTIONAL    :: stream

    ! local variables

    INTEGER                     :: IO_file_id, IO_var_id, i
    INTEGER                     :: g_nland, l_nland
    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim2p(2), dim2(2)
    CHARACTER(max_dim_name) :: dim1n(1), dim2n(2)
    
    !! Check initialization

    IF (module_initialized)  CALL finish("init_dynveg","Attempt to initialize twice.")

    !! Read namelist for dynamical vegetation

    CALL config_dynveg() ! For disturbance-only-case, Read_fpc is the only relevant option

    !! Stream definition
    
    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'dynveg', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_dynveg => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_dynveg, 'dynveg', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(IO_dynveg, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    g_nland = grid%nland
    l_nland = domain%nland

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim2p = (/ l_nland, ntiles /)
    dim2  = (/ g_nland, ntiles /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'tiles'

    ! code numbers of this routine are 31, 34-36, 38-40, 44-45 and 255 using the GRIB veg_table
    CALL add(IO_dynveg,'pot_fpc'         ,dynveg%pot_fpc,     longname='potential plant cover fraction', units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=38, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'bio_exist'       ,dynveg%bio_exist,   longname='bio exist',              units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=39, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'act_fpc'         ,dynveg%act_fpc,      longname='fractional plant cover', units='', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=31, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'max_green_bio', dynveg%max_green_bio, longname='maximum amount of green biomass in the year', &
             units='mol(C) m-2(canopy)', &
             ldims=dim2p, gdims=dim2, dimnames=dim2n, code=36, lmiss=.TRUE., missval=missing_value,   lpost=.false.)
    CALL add(IO_dynveg,'sum_green_bio_memory', dynveg%sum_green_bio_memory, &
             longname='vegetated fraction calculated from green biomass', units='mol(C) m-2(canopy)', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=40, lmiss=.TRUE., missval=missing_value,   lpost=.true.)
    CALL add(IO_dynveg,'desert_fpc'      ,dynveg%desert_fpc,  longname='desert fraction',        units='', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=34, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_dynveg,'bare_fpc'        ,dynveg%bare_fpc,    longname='bare soil fraction',     units='', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=35, lmiss=.TRUE., missval=missing_value)

    IF (debug_Cconservation) THEN
       CALL add(IO_dynveg,'dynveg_testCconserv_2', dynveg%dynveg_testCconserv_2, &
             longname='Second test for carbon conservation in dynveg (should be zero)', &
             units='mol(C)/m^2(vegetated area)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
             code=45, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_dynveg,'dynveg_testCconserv_1', dynveg%dynveg_testCconserv_1, &
             longname='First test for carbon conservation in dynveg (should be zero)', &
             units='mol(C)/m^2(vegetated area)', ldims=dim1p, gdims=dim1, dimnames=dim1n, &
             code=44, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    END IF

    if (.not. run_disturbance) then
       CALL add(IO_dynveg,'zero', dynveg%Dummy_null, longname='Permanent null array for various purposes', &
             units='', ldims=dim2p, gdims=dim2, dimnames=dim2n, &
             code=255, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       dynveg%Dummy_null(:,:) = 0._dp
    endif

   IF (dynveg_options%read_fpc) THEN

       ALLOCATE(init_act_fpc(l_nland,ntiles))
       ALLOCATE(init_pot_fpc(l_nland,ntiles))
       ALLOCATE(init_bare_fpc(l_nland))
       ALLOCATE(init_desert_fpc(l_nland))
       ALLOCATE(init_max_green_bio(l_nland,ntiles))
       ALLOCATE(init_sum_green_bio_memory(l_nland))
       ALLOCATE(init_cover_fract(l_nland,ntiles))
       ALLOCATE(init_cover_fract_pot(l_nland,ntiles))
       ALLOCATE(zzreal2d(domain%ndim,domain%nblocks))

       IF (p_parallel_io) THEN
          fpc_file%opened = .FALSE.
          CALL IO_open(TRIM(fpc_file%file_name), fpc_file, IO_READ)
          IO_file_id = fpc_file%file_id

          ALLOCATE(zreal3d(grid%nlon,grid%nlat,ntiles))
          ALLOCATE(zreal2d(grid%nlon,grid%nlat))
       ENDIF

       !! --- Get initial values for act_fpc (and allocate further memory)

       IF (p_parallel_io) THEN 
          CALL IO_inq_varid(IO_file_id, 'act_fpc', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zzreal2d)
          init_act_fpc(:,i) = PACK(zzreal2d, MASK=domain%mask)
       END DO

       !! --- Get initial values for pot_fpc

       IF (p_parallel_io) THEN 
          CALL IO_inq_varid(IO_file_id, 'pot_fpc', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zzreal2d)
          init_pot_fpc(:,i) = PACK(zzreal2d, MASK=domain%mask)
       END DO

       !! --- Get bare_fpc

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'bare_fpc', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       IF (p_parallel_io) zreal2d_ptr => zreal2d(:,:)
       CALL scatter_field(zreal2d_ptr, zzreal2d)
       init_bare_fpc(:) = PACK(zzreal2d, MASK=domain%mask)

       !! --- Get desert_fpc

       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'desert_fpc', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       IF (p_parallel_io) zreal2d_ptr => zreal2d(:,:)
       CALL scatter_field(zreal2d_ptr, zzreal2d)
       init_desert_fpc(:) = PACK(zzreal2d, MASK=domain%mask)

       !! --- Get initial values for max_green_bio

       IF (p_parallel_io) THEN 
          CALL IO_inq_varid(IO_file_id, 'max_green_bio', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zzreal2d)
          init_max_green_bio(:,i) = PACK(zzreal2d, MASK=domain%mask)
       END DO

       !! --- Get initial values for sum_green_bio_memory

       IF (p_parallel_io) THEN 
          CALL IO_inq_varid(IO_file_id, 'sum_green_bio_memory', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       IF (p_parallel_io) zreal2d_ptr => zreal2d(:,:)
       CALL scatter_field(zreal2d_ptr, zzreal2d)
       init_sum_green_bio_memory(:) = PACK(zzreal2d, MASK=domain%mask)

       !! --- Get initial values for cover_fract

       IF (p_parallel_io) THEN 
          CALL IO_inq_varid(IO_file_id, 'cover_fract', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zzreal2d)
          init_cover_fract(:,i) = PACK(zzreal2d, MASK=domain%mask)
       END DO

       !! --- Get initial values for cover_fract_pot

       IF (p_parallel_io) THEN 
          CALL IO_inq_varid(IO_file_id, 'cover_fract_pot', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zzreal2d)
          init_cover_fract_pot(:,i) = PACK(zzreal2d, MASK=domain%mask)
       END DO

       !! --- Finish

       IF (p_parallel_io) THEN
          CALL IO_close(fpc_file)
          DEALLOCATE(zreal3d)
          DEALLOCATE(zreal2d)
       ENDIF
       DEALLOCATE(zzreal2d)

    ENDIF

    IF (debug_Cconservation) THEN
       dynveg%dynveg_testCconserv_1 = 0._dp
       dynveg%dynveg_testCconserv_2 = 0._dp
    ENDIF

    IF (lrestart) RETURN

    ! -- Initialisation at beginning of an experiment

    dynveg%bio_exist = 1._dp
    dynveg%pot_fpc = 0._dp
    dynveg%max_green_bio = 0._dp
    dynveg%sum_green_bio_memory = 0._dp

    CALL message('init_dynveg','fpc initialised')

  END SUBROUTINE init_dynveg

  ! --- update_dynveg------------------------------------------------------------------------------------------------------


  SUBROUTINE update_dynveg(lstart, lresume, nidx, kidx0, kidx1, ntiles, &
                           use_disturbance, read_cover_fract, with_yasso, with_nitrogen, &
                           lctlib, surf, cbal, climbuf, &
                           sla, lai, veg_fract_correction, co2_emission,  &
                           LeafLit_coef, WoodLit_coef, cbalone_flag, nbalance)

! !DESCRIPTION
! 
! the module simulates vegetation dynamics based on averaged annual npp values
!
!------------------------------------------------------------------------------
    USE mo_time_control,     ONLY: get_time_step, get_year_day, previous_date, delta_time
    USE mo_cbal_cpools,      ONLY: relocate_CarbonAndNitrogen, relocate_carbon_desert
    USE mo_jsbach_constants, ONLY: molarMassCO2_kg
    USE mo_jsbach,           ONLY: debug, new_year, new_day
    USE mo_jsbach_lctlib,    ONLY: lctlib_type
    USE mo_climbuf,          ONLY: climbuf_type

    ! input parameters

    LOGICAL,         INTENT(in) :: lstart                    ! .TRUE. at the beginning of an experiment
    LOGICAL,         INTENT(in) :: lresume                   ! .TRUE. at the beginning of restarted runs
    INTEGER,         INTENT(in) :: nidx                      ! Vector lenght
    INTEGER,         INTENT(in) :: kidx0, kidx1              ! first and last index in the global array
    INTEGER,         INTENT(in) :: ntiles                    ! number of tiles
    LOGICAL,         INTENT(in) :: use_disturbance           ! Calculate disturbances through fire and windbreak
    LOGICAL,         INTENT(in) :: read_cover_fract          ! read cover fractions from initial rather than restart file
    LOGICAL,         INTENT(in) :: with_yasso
    LOGICAL,         INTENT(in) :: with_nitrogen             ! run includes nitrogen cycle
    TYPE(lctlib_type),  INTENT(in) :: lctlib          ! parameters of the land cover type library
    TYPE(land_surface_type), INTENT(inout) :: surf    ! surface parameters 
    TYPE(cbalance_type)    , INTENT(inout) :: cbal    ! c-pool parameters
    TYPE(climbuf_type), INTENT(in) :: climbuf         ! climate parameters
    REAL(dp), INTENT(in) :: sla(:,:)                  ! specific leaf area
    REAL(dp), INTENT(in) :: lai(:,:)                  ! leaf area index
    REAL(dp), INTENT(in) :: veg_fract_correction(:,:) ! corrects vegetated area 1-exp(-LAI_max/2)
    REAL(dp), INTENT(out) :: co2_emission(:)         ! co2-flux to the atmosphere due to fires [kg(CO2)/m^2(gridbox)s]
    LOGICAL, OPTIONAL, INTENT(in) :: cbalone_flag     ! indicates runs with the carbon offline model  
    REAL(dp), OPTIONAL, INTENT(in) :: LeafLit_coef(:,:,:)    ! yasso parameter
    REAL(dp), OPTIONAL, INTENT(in) :: WoodLit_coef(:,:,:)    ! yasso parameter
    TYPE(nbalance_type)    , OPTIONAL, INTENT(inout) :: nbalance   ! n-pool parameters

    ! local variables

    INTEGER    :: istep                               ! current time step (since initialization)
    INTEGER    :: time_steps_per_day                  ! number of time steps per day
    INTEGER    :: time_steps_in_last_year             ! number of time steps within the previous year
    INTEGER    :: i, itile
    LOGICAL    :: glacier(nidx)                       ! two dimensional logical glacier mask
!    REAL(dp)   :: cpool_litter(nidx)                  ! biomass of litter [mol(C)/m^2(grid box)] 
    REAL(dp)   :: veg_ratio_max_old(nidx)             ! fractional cover of vegetated area of last year
    REAL(dp)   :: cover_fract_old(nidx,ntiles)        ! JSBACH cover fractions entering the routine
    REAL(dp)   :: cover_fract_pot_previous_year(nidx,ntiles) ! last years potential cover fractions of natural vegetation
    REAL(dp)   :: act_fpc_previous_day(nidx,ntiles)   ! act_fpc from previous day [] 
    REAL(dp)   :: bare_fpc_previous_day(nidx)         ! bare_fpc from previous day [] 
    LOGICAL    :: cbalone      = .FALSE.              ! flag to indicate runs with the cbalance offline model

    IF (PRESENT(cbalone_flag)) cbalone = cbalone_flag

    !-- define a 1 dimensional glacier mask
    glacier(:) = ANY(surf%is_glacier(kidx0:kidx1,:),DIM=2)

    !-- initialise FPC and memory of vegetated fraction (i.e. desert fraction)

    IF ((lstart .OR. lresume) .AND. .NOT. cbalone) THEN

       IF (dynveg_options%read_fpc) THEN

          dynveg%act_fpc(kidx0:kidx1,:)  = init_act_fpc(kidx0:kidx1,:)
          dynveg%pot_fpc(kidx0:kidx1,:)  = init_pot_fpc(kidx0:kidx1,:)
          dynveg%bare_fpc(kidx0:kidx1)   = init_bare_fpc(kidx0:kidx1)
          dynveg%desert_fpc(kidx0:kidx1) = init_desert_fpc(kidx0:kidx1)
          dynveg%max_green_bio(kidx0:kidx1,:) = init_max_green_bio(kidx0:kidx1,:)
          dynveg%sum_green_bio_memory(kidx0:kidx1) = init_sum_green_bio_memory(kidx0:kidx1)
          surf%cover_fract(kidx0:kidx1,:) = init_cover_fract(kidx0:kidx1,:)
          surf%cover_fract_pot(kidx0:kidx1,:) = init_cover_fract_pot(kidx0:kidx1,:)

          !-- make sure, sum of fpc plus bare fraction is one in each grid cell
          CALL scale_fpc(nidx, ntiles, glacier(:), lctlib%dynamic_pft(:), surf%cover_type(kidx0:kidx1,:), &
                         dynveg%act_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1))

          !-- make sure, sum of cover fractions is one in each grid cell
          CALL scale_cover_fract (nidx, ntiles, surf%is_present(kidx0:kidx1,:), surf%is_glacier(kidx0:kidx1,:), &
                                                surf%is_naturalveg(kidx0:kidx1,:), surf%cover_fract(kidx0:kidx1,:))

          IF (debug) CALL message('update_dynveg','Initial values of FPC read from file' )

       ELSE IF (.NOT. lresume .OR. read_cover_fract) THEN

          !-- calculate initial fpc from cover_fract of the jsbach initial data file by using the map of potential natural vegetation
          dynveg%act_fpc(kidx0:kidx1,:) = surf%cover_fract_pot(kidx0:kidx1,:)
          dynveg%bare_fpc(kidx0:kidx1) = 0._dp
          dynveg%desert_fpc(kidx0:kidx1) =  MIN(1._dp - 2._dp * EPSILON(1._dp),MAX(0._dp,1._dp - surf%veg_ratio_max(kidx0:kidx1)))
          !-- make sure, sum of fpc plus bare fraction is one in each grid cell
          CALL scale_fpc(nidx, ntiles, glacier(:), lctlib%dynamic_pft(:), surf%cover_type(kidx0:kidx1,:), &
                         dynveg%act_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1))
          dynveg%pot_fpc(kidx0:kidx1,:) = MAX(act_fpc_min,dynveg%act_fpc(kidx0:kidx1,:))
          IF (debug) CALL message('update_dynveg','FPC initialized' )

          dynveg%sum_green_bio_memory(kidx0:kidx1) = MIN(1._dp,MAX(0._dp,surf%veg_ratio_max(kidx0:kidx1)))

          WHERE (glacier(:))
             surf%rock_fract(kidx0:kidx1) = 1._dp
          END WHERE

       END IF

       IF (dynveg_options%read_fpc .and. dynveg_options%dynveg_feedback) THEN

          !-- for consitency, calcuate maximum vegetated fraction and cover fractions from FPCs

          cover_fract_pot_previous_year(:,:) = surf%cover_fract_pot(kidx0:kidx1,:)
          CALL calc_veg_ratio_max(dynveg%desert_fpc(kidx0:kidx1), surf%rock_fract(kidx0:kidx1), &
                                  surf%is_glacier(kidx0:kidx1,:), surf%veg_ratio_max(kidx0:kidx1))
          CALL cover_fract_pot_to_cover_fract(nidx, ntiles, lctlib%nlct, surf%is_present(kidx0:kidx1,:), &
               surf%is_glacier(kidx0:kidx1,:), surf%is_naturalveg(kidx0:kidx1,:), &
               lctlib%woody_pft(:), lctlib%dynamic_pft(:), lctlib%pasture_pft(:), &
               surf%is_crop(kidx0:kidx1,:), surf%cover_type(kidx0:kidx1,:), &
               surf%cover_fract_pot(kidx0:kidx1,:), cover_fract_pot_previous_year(:,:), surf%cover_fract(kidx0:kidx1,:))
       ENDIF

    ENDIF

    !-- annual calculations

    IF (new_year) THEN

       !-- find out time step and whether new year just started
       istep=get_time_step()
       time_steps_per_day = NINT(86400._dp / delta_time)
       time_steps_in_last_year = NINT(get_year_day(previous_date)) * time_steps_per_day

       IF (istep >= time_steps_in_last_year - time_steps_per_day .OR. cbalone) THEN  ! not before one year was completed

          ! -- calculation of bioclimatic limits
          IF (.NOT. cbalone) THEN
             CALL bioclim_limits(nidx, ntiles, lctlib, climbuf%min_mmtemp20(kidx0:kidx1), climbuf%max_mmtemp20(kidx0:kidx1), &
                                 climbuf%prev_year_gdd(kidx0:kidx1,:), surf%cover_type(kidx0:kidx1,:), &
                                 dynveg%bio_exist(kidx0:kidx1,:))
          END IF
    
          ! -- potential FPC
          CALL potential_tree_fpc(nidx, ntiles, dynveg_options%dynveg_all, &
                                  lctlib%woody_pft(:), lctlib%dynamic_pft(:), &
                                  climbuf%ave_npp5(kidx0:kidx1,:), dynveg%bio_exist(kidx0:kidx1,:), &
                                  dynveg%act_fpc(kidx0:kidx1,:), surf%cover_type(kidx0:kidx1,:), dynveg%pot_fpc(kidx0:kidx1,:))

          ! -- calculation of desert extend
          CALL desert_fraction(nidx, ntiles, lctlib%nlct, dynveg_options%accelerate_dynveg, &
                               lctlib%woody_pft(:), lctlib%dynamic_pft(:), &
                               dynveg%act_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1), sla(:,:), &
                               dynveg%max_green_bio(kidx0:kidx1,:), surf%cover_type(kidx0:kidx1,:), &
                               dynveg%sum_green_bio_memory(kidx0:kidx1), dynveg%desert_fpc(kidx0:kidx1))
       ENDIF

       IF (dynveg_options%dynveg_feedback) THEN

          ! -- update veg_ratio_max from the desert fraction

          veg_ratio_max_old(:) = surf%veg_ratio_max(kidx0:kidx1)
          CALL calc_veg_ratio_max(dynveg%desert_fpc(kidx0:kidx1), surf%rock_fract(kidx0:kidx1), &
                                  surf%is_glacier(kidx0:kidx1,:), surf%veg_ratio_max(kidx0:kidx1))

          ! -- conversion of FPCs to cover fractions of potential natural vegetation

          cover_fract_pot_previous_year(:,:) = surf%cover_fract_pot(kidx0:kidx1,:)
          CALL fpc_to_cover_fract_pot(nidx, ntiles, dynveg%act_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1), &
                                      surf%is_present(kidx0:kidx1,:), surf%is_glacier(kidx0:kidx1,:), &
                                      surf%is_naturalveg(kidx0:kidx1,:), lctlib%woody_pft(:), &
                                      lctlib%dynamic_pft(:), surf%cover_type(kidx0:kidx1,:), surf%cover_fract_pot(kidx0:kidx1,:))

          ! -- conversion of potential natural vegetation to actual cover fractions including agricultural areas
          cover_fract_old(:,:) = surf%cover_fract(kidx0:kidx1,:)
          CALL cover_fract_pot_to_cover_fract(nidx, ntiles, lctlib%nlct, surf%is_present(kidx0:kidx1,:), &
                                              surf%is_glacier(kidx0:kidx1,:), surf%is_naturalveg(kidx0:kidx1,:), &
                                              lctlib%woody_pft(:), lctlib%dynamic_pft(:), lctlib%pasture_pft(:),   &
                                              surf%is_crop(kidx0:kidx1,:), surf%cover_type(kidx0:kidx1,:), &
                                              surf%cover_fract_pot(kidx0:kidx1,:),  &
                                              cover_fract_pot_previous_year(:,:), surf%cover_fract(kidx0:kidx1,:))

          IF (debug_Cconservation) THEN

             ! -- initialize first carbon conservation test
             IF (.NOT.with_yasso) THEN
                dynveg%dynveg_testCconserv_1(kidx0:kidx1) =  veg_ratio_max_old(:) *                  &
                     SUM( cover_fract_old(:,:) * veg_fract_correction(:,:) *                         &
                          (  cbal%Cpool_green(kidx0:kidx1,:) + cbal%Cpool_woods(kidx0:kidx1,:) + cbal%Cpool_reserve(kidx0:kidx1,:) &
                           + cbal%Cpool_litter_green_ag(kidx0:kidx1,:)        &
                           + cbal%Cpool_litter_green_bg(kidx0:kidx1,:)        &
                           + cbal%Cpool_litter_wood_ag(kidx0:kidx1,:)         &
                           + cbal%Cpool_litter_wood_bg(kidx0:kidx1,:)         &
                           + cbal%Cpool_slow(kidx0:kidx1,:)                   &
                          ),DIM=2 )
             ELSE
                dynveg%dynveg_testCconserv_1(kidx0:kidx1) =  veg_ratio_max_old(:) *                  &
                     SUM( cover_fract_old(:,:) * veg_fract_correction(:,:) *                         &
                          (  cbal%Cpool_green(kidx0:kidx1,:) + cbal%Cpool_woods(kidx0:kidx1,:) + cbal%Cpool_reserve(kidx0:kidx1,:) &
                           + cbal%YCpool_acid_ag1(kidx0:kidx1,:)       &
                           + cbal%YCpool_water_ag1(kidx0:kidx1,:)      &
                           + cbal%YCpool_ethanol_ag1(kidx0:kidx1,:)    &
                           + cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:) &
                           + cbal%YCpool_acid_bg1(kidx0:kidx1,:)       &
                           + cbal%YCpool_water_bg1(kidx0:kidx1,:)      &
                           + cbal%YCpool_ethanol_bg1(kidx0:kidx1,:)    &
                           + cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:) &
                           + cbal%YCpool_humus_1(kidx0:kidx1,:)        &
                           + cbal%YCpool_acid_ag2(kidx0:kidx1,:)       &
                           + cbal%YCpool_water_ag2(kidx0:kidx1,:)      &
                           + cbal%YCpool_ethanol_ag2(kidx0:kidx1,:)    &
                           + cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:) &
                           + cbal%YCpool_acid_bg2(kidx0:kidx1,:)       &
                           + cbal%YCpool_water_bg2(kidx0:kidx1,:)      &
                           + cbal%YCpool_ethanol_bg2(kidx0:kidx1,:)    &
                           + cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:) &
                           + cbal%YCpool_humus_2(kidx0:kidx1,:)        &
                          ),DIM=2 )
             END IF
          END IF

          ! -- relocate carbon accouting for the change in the cover fraction of each PFT
          IF (.NOT. with_yasso) THEN
             ! -- dynveg with N
             IF (with_nitrogen) THEN ! with N, no yasso

                CALL relocate_CarbonAndNitrogen(lctlib, surf, cover_fract_old(:,:), surf%cover_fract(kidx0:kidx1,:), &
                                                veg_fract_correction(:,:),                                           &
                                                cbal%Cpool_green(kidx0:kidx1,:),                                     &
                                                cbal%Cpool_woods(kidx0:kidx1,:),                                     &
                                                cbal%Cpool_reserve(kidx0:kidx1,:),                                   &
                        Cpool_litter_green_ag = cbal%Cpool_litter_green_ag(kidx0:kidx1,:),                           &
                        Cpool_litter_green_bg = cbal%Cpool_litter_green_bg(kidx0:kidx1,:),                           &
                         Cpool_litter_wood_ag = cbal%Cpool_litter_wood_ag(kidx0:kidx1,:),                            &
                         Cpool_litter_wood_bg = cbal%Cpool_litter_wood_bg(kidx0:kidx1,:),                            &
                                   Cpool_slow = cbal%Cpool_slow(kidx0:kidx1,:),                                      &
                                  Npool_green = nbalance%Npool_green(kidx0:kidx1,:),                                 &
                                  Npool_woods = nbalance%Npool_woods(kidx0:kidx1,:),                                 &
                                 Npool_mobile = nbalance%Npool_mobile(kidx0:kidx1,:),                                &
                        Npool_litter_green_ag = nbalance%Npool_litter_green_ag(kidx0:kidx1,:),                       &
                        Npool_litter_green_bg = nbalance%Npool_litter_green_bg(kidx0:kidx1,:),                       &
                         Npool_litter_wood_ag = nbalance%Npool_litter_wood_ag(kidx0:kidx1,:),                        &
                         Npool_litter_wood_bg = nbalance%Npool_litter_wood_bg(kidx0:kidx1,:),                        &
                                   Npool_slow = nbalance%Npool_slow(kidx0:kidx1,:),                                  &
                                   SMINN_pool = nbalance%SMINN_pool(kidx0:kidx1,:),                                  &
                                   lcc_scheme = 1 )

                CALL relocate_carbon_desert(nidx, ntiles, surf%is_present(kidx0:kidx1,:), surf%is_glacier(kidx0:kidx1,:), &
                                            surf%veg_ratio_max(kidx0:kidx1), veg_ratio_max_old(:), lai(:,:), sla(:,:), &
                                            cbal%Cpool_green(kidx0:kidx1,:), cbal%Cpool_woods(kidx0:kidx1,:),          &
                                            cbal%Cpool_reserve(kidx0:kidx1,:),                                         &
                    Cpool_litter_green_ag = cbal%Cpool_litter_green_ag(kidx0:kidx1,:),                                 &
                    Cpool_litter_green_bg = cbal%Cpool_litter_green_bg(kidx0:kidx1,:),                                 &
                     Cpool_litter_wood_ag = cbal%Cpool_litter_wood_ag(kidx0:kidx1,:),                                  &
                     Cpool_litter_wood_bg = cbal%Cpool_litter_wood_bg(kidx0:kidx1,:),                                  &
                               Cpool_slow = cbal%Cpool_slow(kidx0:kidx1,:),                                            &
                              Npool_green = nbalance%Npool_green(kidx0:kidx1,:),                                       &
                              Npool_woods = nbalance%Npool_woods(kidx0:kidx1,:),                                       &
                             Npool_mobile = nbalance%Npool_mobile(kidx0:kidx1,:),                                      &
                    Npool_litter_green_ag = nbalance%Npool_litter_green_ag(kidx0:kidx1,:),                             &
                    Npool_litter_green_bg = nbalance%Npool_litter_green_bg(kidx0:kidx1,:),                             &
                     Npool_litter_wood_ag = nbalance%Npool_litter_wood_ag(kidx0:kidx1,:),                              &
                     Npool_litter_wood_bg = nbalance%Npool_litter_wood_bg(kidx0:kidx1,:),                              &
                               Npool_slow = nbalance%Npool_slow(kidx0:kidx1,:),                                        &
                               SMINN_pool = nbalance%SMINN_pool(kidx0:kidx1,:))

             ELSE ! no N, no yasso

                CALL relocate_CarbonAndNitrogen(lctlib, surf, cover_fract_old(:,:), surf%cover_fract(kidx0:kidx1,:), &
                                                veg_fract_correction(:,:),                                           &
                                                cbal%Cpool_green(kidx0:kidx1,:),                                     & 
                                                cbal%Cpool_woods(kidx0:kidx1,:),                                     &
                                                cbal%Cpool_reserve(kidx0:kidx1,:),                                   &
                        Cpool_litter_green_ag = cbal%Cpool_litter_green_ag(kidx0:kidx1,:),                           &
                        Cpool_litter_green_bg = cbal%Cpool_litter_green_bg(kidx0:kidx1,:),                           &
                         Cpool_litter_wood_ag = cbal%Cpool_litter_wood_ag(kidx0:kidx1,:),                            &
                         Cpool_litter_wood_bg = cbal%Cpool_litter_wood_bg(kidx0:kidx1,:),                            &
                                   Cpool_slow = cbal%Cpool_slow(kidx0:kidx1,:),                                      &
                                   lcc_scheme = 1)
   
                ! -- scaling cpools to account for change in vegetated fraction in order to conserve the total mass of carbon
                CALL relocate_carbon_desert(nidx, ntiles, surf%is_present(kidx0:kidx1,:), surf%is_glacier(kidx0:kidx1,:), &
                                            surf%veg_ratio_max(kidx0:kidx1), veg_ratio_max_old(:), lai(:,:), sla(:,:),    &
                                            cbal%Cpool_green(kidx0:kidx1,:), cbal%Cpool_woods(kidx0:kidx1,:),             &
                                            cbal%Cpool_reserve(kidx0:kidx1,:),                                            &
                    Cpool_litter_green_ag = cbal%Cpool_litter_green_ag(kidx0:kidx1,:),                                    &
                    Cpool_litter_green_bg = cbal%Cpool_litter_green_bg(kidx0:kidx1,:),                                    &
                     Cpool_litter_wood_ag = cbal%Cpool_litter_wood_ag(kidx0:kidx1,:),                                     &
                     Cpool_litter_wood_bg = cbal%Cpool_litter_wood_bg(kidx0:kidx1,:),                                     &
                               Cpool_slow = cbal%Cpool_slow(kidx0:kidx1,:))
             END IF ! nitrogen
          ELSE ! no N, yasso

             ! -- relocate carbon accouting for the change in the cover fraction of each PFT
             CALL relocate_CarbonAndNitrogen(lctlib, surf, cover_fract_old(:,:),        &
                                             surf%cover_fract(kidx0:kidx1,:),           &
                                             veg_fract_correction(:,:),                 &
                                             cbal%Cpool_green(kidx0:kidx1,:),           & 
                                             cbal%Cpool_woods(kidx0:kidx1,:),           &
                                             cbal%Cpool_reserve(kidx0:kidx1,:),         &
                     YCpool_acid_ag1       = cbal%YCpool_acid_ag1(kidx0:kidx1,:),       &
                     YCpool_water_ag1      = cbal%YCpool_water_ag1(kidx0:kidx1,:),      &
                     YCpool_ethanol_ag1    = cbal%YCpool_ethanol_ag1(kidx0:kidx1,:),    &
                     YCpool_nonsoluble_ag1 = cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:), &
                     YCpool_acid_bg1       = cbal%YCpool_acid_bg1(kidx0:kidx1,:),       &
                     YCpool_water_bg1      = cbal%YCpool_water_bg1(kidx0:kidx1,:),      &
                     YCpool_ethanol_bg1    = cbal%YCpool_ethanol_bg1(kidx0:kidx1,:),    &
                     YCpool_nonsoluble_bg1 = cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:), &
                     YCpool_humus_1        = cbal%YCpool_humus_1(kidx0:kidx1,:),        &
                     YCpool_acid_ag2       = cbal%YCpool_acid_ag2(kidx0:kidx1,:),       &
                     YCpool_water_ag2      = cbal%YCpool_water_ag2(kidx0:kidx1,:),      &
                     YCpool_ethanol_ag2    = cbal%YCpool_ethanol_ag2(kidx0:kidx1,:),    &
                     YCpool_nonsoluble_ag2 = cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:), &
                     YCpool_acid_bg2       = cbal%YCpool_acid_bg2(kidx0:kidx1,:),       &
                     YCpool_water_bg2      = cbal%YCpool_water_bg2(kidx0:kidx1,:),      &
                     YCpool_ethanol_bg2    = cbal%YCpool_ethanol_bg2(kidx0:kidx1,:),    &
                     YCpool_nonsoluble_bg2 = cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:), &
                     YCpool_humus_2        = cbal%YCpool_humus_2(kidx0:kidx1,:),        &
                     LeafLit_coef          = LeafLit_coef(1:nidx,:,:),                  &
                     lcc_scheme            = 1)
   
             ! -- scaling cpools to account for change in vegetated fraction in order to conserve the total mass of carbon
             CALL relocate_carbon_desert(nidx, ntiles, surf%is_present(kidx0:kidx1,:), surf%is_glacier(kidx0:kidx1,:), &
                                         surf%veg_ratio_max(kidx0:kidx1), veg_ratio_max_old(:), lai(:,:), sla(:,:),    &
                                         cbal%Cpool_green(kidx0:kidx1,:), cbal%Cpool_woods(kidx0:kidx1,:),             &
                                         cbal%Cpool_reserve(kidx0:kidx1,:),                                            &
                 YCpool_acid_ag1       = cbal%YCpool_acid_ag1(kidx0:kidx1,:),                                          &
                 YCpool_water_ag1      = cbal%YCpool_water_ag1(kidx0:kidx1,:),                                         &
                 YCpool_ethanol_ag1    = cbal%YCpool_ethanol_ag1(kidx0:kidx1,:),                                       &
                 YCpool_nonsoluble_ag1 = cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:),                                    &
                 YCpool_acid_bg1       = cbal%YCpool_acid_bg1(kidx0:kidx1,:),                                          &
                 YCpool_water_bg1      = cbal%YCpool_water_bg1(kidx0:kidx1,:),                                         &
                 YCpool_ethanol_bg1    = cbal%YCpool_ethanol_bg1(kidx0:kidx1,:),                                       &
                 YCpool_nonsoluble_bg1 = cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:),                                    &
                 YCpool_humus_1        = cbal%YCpool_humus_1(kidx0:kidx1,:),                                           &
                 YCpool_acid_ag2       = cbal%YCpool_acid_ag2(kidx0:kidx1,:),                                          &
                 YCpool_water_ag2      = cbal%YCpool_water_ag2(kidx0:kidx1,:),                                         &
                 YCpool_ethanol_ag2    = cbal%YCpool_ethanol_ag2(kidx0:kidx1,:),                                       &
                 YCpool_nonsoluble_ag2 = cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:),                                    &
                 YCpool_acid_bg2       = cbal%YCpool_acid_bg2(kidx0:kidx1,:),                                          &
                 YCpool_water_bg2      = cbal%YCpool_water_bg2(kidx0:kidx1,:),                                         &
                 YCpool_ethanol_bg2    = cbal%YCpool_ethanol_bg2(kidx0:kidx1,:),                                       &
                 YCpool_nonsoluble_bg2 = cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:),                                    &
                 YCpool_humus_2        = cbal%YCpool_humus_2(kidx0:kidx1,:),                                           &
                 LeafLitcoef           = LeafLit_coef(1:nidx,:,:))
          END IF ! yasso

          IF (debug_Cconservation) THEN

             ! -- finish carbon conservation test
             IF (.NOT. with_yasso) THEN
                dynveg%dynveg_testCconserv_1(kidx0:kidx1) =  dynveg%dynveg_testCconserv_1(kidx0:kidx1) &
                   - surf%veg_ratio_max(kidx0:kidx1) * SUM( surf%cover_fract(kidx0:kidx1,:) * veg_fract_correction(:,:) *  &
                          (  cbal%Cpool_green(kidx0:kidx1,:) + cbal%Cpool_woods(kidx0:kidx1,:) + cbal%Cpool_reserve(kidx0:kidx1,:) &
                           + cbal%Cpool_litter_green_ag(kidx0:kidx1,:) + cbal%Cpool_litter_green_bg(kidx0:kidx1,:)         &
                           + cbal%Cpool_litter_wood_ag(kidx0:kidx1,:) + cbal%Cpool_litter_wood_bg(kidx0:kidx1,:)           &
                           + cbal%Cpool_slow(kidx0:kidx1,:)                                                      &
                          ),DIM=2 ) 
             ELSE 
                dynveg%dynveg_testCconserv_1(kidx0:kidx1) =  dynveg%dynveg_testCconserv_1(kidx0:kidx1) &
                   - surf%veg_ratio_max(kidx0:kidx1) * SUM( surf%cover_fract(kidx0:kidx1,:) * veg_fract_correction(:,:) *  &
                          (  cbal%Cpool_green(kidx0:kidx1,:) + cbal%Cpool_woods(kidx0:kidx1,:) + cbal%Cpool_reserve(kidx0:kidx1,:) &
                           + cbal%YCpool_acid_ag1(kidx0:kidx1,:)       &
                           + cbal%YCpool_water_ag1(kidx0:kidx1,:)      &
                           + cbal%YCpool_ethanol_ag1(kidx0:kidx1,:)    &
                           + cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:) &
                           + cbal%YCpool_acid_bg1(kidx0:kidx1,:)       &
                           + cbal%YCpool_water_bg1(kidx0:kidx1,:)      &
                           + cbal%YCpool_ethanol_bg1(kidx0:kidx1,:)    &
                           + cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:) &
                           + cbal%YCpool_humus_1(kidx0:kidx1,:)        &
                           + cbal%YCpool_acid_ag2(kidx0:kidx1,:)       &
                           + cbal%YCpool_water_ag2(kidx0:kidx1,:)      &
                           + cbal%YCpool_ethanol_ag2(kidx0:kidx1,:)    &
                           + cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:) &
                           + cbal%YCpool_acid_bg2(kidx0:kidx1,:)       &
                           + cbal%YCpool_water_bg2(kidx0:kidx1,:)      &
                           + cbal%YCpool_ethanol_bg2(kidx0:kidx1,:)    &
                           + cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:) &
                           + cbal%YCpool_humus_2(kidx0:kidx1,:)        &
                          ),DIM=2 ) 
             ENDIF
          END IF

       ENDIF

    ENDIF

    ! -- daily calculations

    IF (new_day) THEN
       ! -- determine maximum of green biomass for the current year

       IF (new_year) dynveg%max_green_bio(kidx0:kidx1,:) = 0._dp
       DO i = kidx0,kidx1
          DO itile = 1,ntiles
             IF (lctlib%dynamic_pft(surf%cover_type(i,itile))) &
                dynveg%max_green_bio(i,itile) = MAX(dynveg%max_green_bio(i,itile),cbal%cpool_green(i,itile))
          END DO
       END DO

       ! -- Save FPC-fractions from past time step (day) for consistent carbon relocation

       bare_fpc_previous_day(1:nidx)  = dynveg%bare_fpc(kidx0:kidx1)
       act_fpc_previous_day(1:nidx,:) = dynveg%act_fpc(kidx0:kidx1,:)

       ! -- calculation of FPC (competition of vegetation types and disturbances (fire, wind break))
       IF (use_disturbance) THEN
          CALL fpc_daily(lctlib, nidx, kidx0, kidx1, ntiles, with_yasso, dynveg_options%dynveg_all, &
                         dynveg_options%accelerate_dynveg,                                          &
                         lctlib%woody_pft(:), lctlib%dynamic_pft(:), lctlib%tau_Cpool_woods(:),     &
                         climbuf%ave_npp5(kidx0:kidx1,:), dynveg%bio_exist(kidx0:kidx1,:),          &
                         dynveg%pot_fpc(kidx0:kidx1,:), veg_fract_correction, surf, cbal, climbuf,  &
                         act_fpc_previous_day(1:nidx,:), dynveg%act_fpc(kidx0:kidx1,:),             &
                         dynveg%bare_fpc(kidx0:kidx1),                                              &
                         disturbance%burned_frac(kidx0:kidx1,:), disturbance%burned_frac_diag(kidx0:kidx1,:), &
                         disturbance%damaged_frac(kidx0:kidx1,:),disturbance%fuel(kidx0:kidx1))

          ! Cumulate disturbed areas for diagnostics
          IF (.NOT. ASSOCIATED(disturbance%box_burned_frac_diag_avg, disturbance%box_burned_frac_avg))             &
               disturbance%box_burned_frac_diag_avg(kidx0:kidx1,:) =                                               &
               disturbance%box_burned_frac_diag_avg(kidx0:kidx1,:) + disturbance%burned_frac_diag(kidx0:kidx1,:)   &
                                                                 * surf%cover_fract          (kidx0:kidx1,:)       &
                                                                 * SPREAD(surf%veg_ratio_max (kidx0:kidx1)         &
                                                                          * (1.0_dp-bare_fpc_previous_day(1:nidx)) &
                                                                          ,DIM=2,ncopies=ntiles)                   &
                                                                 * 86400._dp
          disturbance%box_burned_frac_avg      (kidx0:kidx1,:) =                                                   &
          disturbance%box_burned_frac_avg      (kidx0:kidx1,:) + disturbance%burned_frac    (kidx0:kidx1,:)        &
                                                                 * surf%cover_fract         (kidx0:kidx1,:)        &
                                                                 * SPREAD(surf%veg_ratio_max(kidx0:kidx1)          &
                                                                          * (1.0_dp-bare_fpc_previous_day(1:nidx)) &
                                                                          ,DIM=2,ncopies=ntiles)                   &
                                                                 * 86400._dp
          disturbance%box_damaged_frac_avg     (kidx0:kidx1,:) =                                                   &
          disturbance%box_damaged_frac_avg     (kidx0:kidx1,:) + disturbance%damaged_frac   (kidx0:kidx1,:)        &
                                                                 * surf%cover_fract         (kidx0:kidx1,:)        &
                                                                 * SPREAD(surf%veg_ratio_max(kidx0:kidx1)          &
                                                                          * (1.0_dp-bare_fpc_previous_day(1:nidx)) &
                                                                          ,DIM=2,ncopies=ntiles)                   &
                                                                 * 86400._dp
       ELSE ! if (use_disturbance) else
          CALL fpc_daily(lctlib, nidx, kidx0, kidx1, ntiles, with_yasso, dynveg_options%dynveg_all, &
                         dynveg_options%accelerate_dynveg,                                          &
                         lctlib%woody_pft(:), lctlib%dynamic_pft(:), lctlib%tau_Cpool_woods(:),     &
                         climbuf%ave_npp5(kidx0:kidx1,:), dynveg%bio_exist(kidx0:kidx1,:),          &
                         dynveg%pot_fpc(kidx0:kidx1,:), veg_fract_correction, surf, cbal, climbuf,  &
                         act_fpc_previous_day(1:nidx,:), dynveg%act_fpc(kidx0:kidx1,:),             &
                         dynveg%bare_fpc(kidx0:kidx1))
       ENDIF

       ! -- rescaling of actual fpc and bare fpc
 
       CALL scale_fpc(nidx, ntiles, glacier(:), lctlib%dynamic_pft(:), surf%cover_type(kidx0:kidx1,:), &
                      dynveg%act_fpc(kidx0:kidx1,:), dynveg%bare_fpc(kidx0:kidx1))

       IF (use_disturbance) THEN
         IF (debug_Cconservation) THEN

           ! -- initialize second carbon conservation test
           IF (.NOT.with_yasso) THEN
               dynveg%dynveg_testCconserv_2(kidx0:kidx1) =                                          &
                    SUM( surf%cover_fract(kidx0:kidx1,:) * veg_fract_correction(:,:) *              &
                         (  cbal%Cpool_green(kidx0:kidx1,:) + cbal%Cpool_woods(kidx0:kidx1,:)       & 
                          + cbal%Cpool_reserve(kidx0:kidx1,:) &
                          + cbal%Cpool_litter_green_ag(kidx0:kidx1,:) + cbal%Cpool_litter_green_bg(kidx0:kidx1,:)       &
                          + cbal%Cpool_litter_wood_ag(kidx0:kidx1,:) + cbal%Cpool_litter_wood_bg(kidx0:kidx1,:)         &
                          + cbal%Cpool_slow(kidx0:kidx1,:) &
                         ),DIM=2 ) 
            ELSE
               dynveg%dynveg_testCconserv_2(kidx0:kidx1) =                                          &
                    SUM( surf%cover_fract(kidx0:kidx1,:) * veg_fract_correction(:,:) *              &
                         (  cbal%Cpool_green(kidx0:kidx1,:) + cbal%Cpool_woods(kidx0:kidx1,:)       & 
                          + cbal%Cpool_reserve(kidx0:kidx1,:) &
                          + cbal%YCpool_acid_ag1(kidx0:kidx1,:)       &
                          + cbal%YCpool_water_ag1(kidx0:kidx1,:)      &
                          + cbal%YCpool_ethanol_ag1(kidx0:kidx1,:)    &
                          + cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:) &
                          + cbal%YCpool_acid_bg1(kidx0:kidx1,:)       &
                          + cbal%YCpool_water_bg1(kidx0:kidx1,:)      &
                          + cbal%YCpool_ethanol_bg1(kidx0:kidx1,:)    &
                          + cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:) &
                          + cbal%YCpool_humus_1(kidx0:kidx1,:)        &
                          + cbal%YCpool_acid_ag2(kidx0:kidx1,:)       &
                          + cbal%YCpool_water_ag2(kidx0:kidx1,:)      &
                          + cbal%YCpool_ethanol_ag2(kidx0:kidx1,:)    &
                          + cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:) &
                          + cbal%YCpool_acid_bg2(kidx0:kidx1,:)       &
                          + cbal%YCpool_water_bg2(kidx0:kidx1,:)      &
                          + cbal%YCpool_ethanol_bg2(kidx0:kidx1,:)    &
                          + cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:) &
                          + cbal%YCpool_humus_2(kidx0:kidx1,:)        &
                         ),DIM=2 ) 
            END IF
         END IF
         ! -- changes in carbon pools caused by wind break (to be called before the routine relocate_carbon_fire)
         IF (with_nitrogen) THEN ! with N, no yasso
            CALL relocate_disturbed_carbon(lctlib,surf,nidx,kidx0,ntiles,with_yasso,DIST_WINDBREAK,                  &
                                           disturbance%damaged_frac(kidx0:kidx1,:),                                  &
                                           surf%cover_fract(kidx0:kidx1,:),veg_fract_correction(:,:),                &
                                           cbal%Cpool_green(kidx0:kidx1,:), cbal%Cpool_woods(kidx0:kidx1,:),         &
                                           cbal%Cpool_reserve(kidx0:kidx1,:),                                        &
                                           disturbance%carbon_2_GreenLitterPools(kidx0:kidx1),                       &
                                           disturbance%carbon_2_WoodLitterPools(kidx0:kidx1),                        &
                                           cbal%Cpool_litter_green_ag(kidx0:kidx1,:),                                &
                                           cbal%Cpool_litter_green_bg(kidx0:kidx1,:),                                &
                                           cbal%Cpool_litter_wood_ag(kidx0:kidx1,:),                                 &
                                           cbal%Cpool_litter_wood_bg(kidx0:kidx1,:),                                 &
                                           cbal%YCpool_acid_ag1(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_acid_ag2(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_water_ag1(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_water_ag2(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_ethanol_ag1(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_ethanol_ag2(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:),                                &
                                           cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:),                                &
                                           cbal%YCpool_acid_bg1(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_acid_bg2(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_water_bg1(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_water_bg2(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_ethanol_bg1(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_ethanol_bg2(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:),                                &
                                           cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:),                                &
                                           cbal%YCpool_humus_1(kidx0:kidx1,:),                                       & 
                                           cbal%YCpool_humus_2(kidx0:kidx1,:),                                       &
                                           LeafLit_coef(1:nidx,:,:),                                                 &
                                           WoodLit_coef(1:nidx,:,:),                                                 &
                             Npool_green = nbalance%Npool_green(kidx0:kidx1,:),                                      &
                             Npool_woods = nbalance%Npool_woods(kidx0:kidx1,:),                                      &
                            Npool_mobile = nbalance%Npool_mobile(kidx0:kidx1,:),                                     &
                              sminn_pool = nbalance%sminn_pool(kidx0:kidx1,:),                                       &
                   Npool_litter_green_ag = nbalance%Npool_litter_green_ag(kidx0:kidx1,:),                            &
                   Npool_litter_green_bg = nbalance%Npool_litter_green_bg(kidx0:kidx1,:),                            &
                    Npool_litter_wood_ag = nbalance%Npool_litter_wood_ag(kidx0:kidx1,:),                             &
                    Npool_litter_wood_bg = nbalance%Npool_litter_wood_bg(kidx0:kidx1,:),                             &
             nitrogen_2_GreenLitterPools = disturbance%nitrogen_2_GreenLitterPools(kidx0:kidx1),                     &
              nitrogen_2_WoodLitterPools = disturbance%nitrogen_2_WoodLitterPools(kidx0:kidx1),                      &
                        nitrogen_2_sminn = disturbance%nitrogen_2_sminn(kidx0:kidx1))

            ! -- changes in carbon pools caused by fire
            CALL relocate_disturbed_carbon(lctlib,surf,nidx,kidx0,ntiles,with_yasso,DIST_FIRE,                       &
                                           disturbance%burned_frac(kidx0:kidx1,:),                                   &
                                           surf%cover_fract(kidx0:kidx1,:),veg_fract_correction(:,:),                &
                                           cbal%Cpool_green(kidx0:kidx1,:), cbal%Cpool_woods(kidx0:kidx1,:),         &
                                           cbal%Cpool_reserve(kidx0:kidx1,:),                                        &
                                           disturbance%carbon_2_GreenLitterPools(kidx0:kidx1),                       &
                                           disturbance%carbon_2_WoodLitterPools(kidx0:kidx1),                        &
                                           cbal%Cpool_litter_green_ag(kidx0:kidx1,:),                                &
                                           cbal%Cpool_litter_green_bg(kidx0:kidx1,:),                                &
                                           cbal%Cpool_litter_wood_ag(kidx0:kidx1,:),                                 &
                                           cbal%Cpool_litter_wood_bg(kidx0:kidx1,:),                                 &
                                           cbal%YCpool_acid_ag1(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_acid_ag2(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_water_ag1(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_water_ag2(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_ethanol_ag1(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_ethanol_ag2(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:),                                &
                                           cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:),                                &
                                           cbal%YCpool_acid_bg1(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_acid_bg2(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_water_bg1(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_water_bg2(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_ethanol_bg1(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_ethanol_bg2(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:),                                &
                                           cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:),                                &
                                           cbal%YCpool_humus_1(kidx0:kidx1,:),                                       & 
                                           cbal%YCpool_humus_2(kidx0:kidx1,:),                                       &
                                           LeafLit_coef(1:nidx,:,:),                                                 &
                                           WoodLit_coef(1:nidx,:,:),                                                 &
                          carbon_2_atmos = disturbance%carbon_2_atmos(kidx0:kidx1),                                  &
                             Npool_green = nbalance%Npool_green(kidx0:kidx1,:),                                      &
                             Npool_woods = nbalance%Npool_woods(kidx0:kidx1,:),                                      &
                            Npool_mobile = nbalance%Npool_mobile(kidx0:kidx1,:),                                     &
                              SMINN_pool = nbalance%SMINN_pool(kidx0:kidx1,:),                                       &
                   Npool_litter_green_ag = nbalance%Npool_litter_green_ag(kidx0:kidx1,:),                            &
                   Npool_litter_green_bg = nbalance%Npool_litter_green_bg(kidx0:kidx1,:),                            &
                    Npool_litter_wood_ag = nbalance%Npool_litter_wood_ag(kidx0:kidx1,:),                             &
                    Npool_litter_wood_bg = nbalance%Npool_litter_wood_bg(kidx0:kidx1,:),                             &
             nitrogen_2_GreenLitterPools = disturbance%nitrogen_2_GreenLitterPools(kidx0:kidx1),                     &
              nitrogen_2_WoodLitterPools = disturbance%nitrogen_2_WoodLitterPools(kidx0:kidx1),                      &
                        nitrogen_2_atmos = disturbance%nitrogen_2_atmos(kidx0:kidx1),                                &
                        nitrogen_2_sminn = disturbance%nitrogen_2_sminn(kidx0:kidx1))                                
         ELSE ! no nitrogen

            ! -- changes in carbon pools caused by wind break (to be called before the routine relocate_carbon_fire)
            CALL relocate_disturbed_carbon(lctlib,surf,nidx,kidx0,ntiles,with_yasso,DIST_WINDBREAK,                  &
                                           disturbance%damaged_frac(kidx0:kidx1,:),                                  &
                                           surf%cover_fract(kidx0:kidx1,:),veg_fract_correction(:,:),                &
                                           cbal%Cpool_green(kidx0:kidx1,:), cbal%Cpool_woods(kidx0:kidx1,:),         &
                                           cbal%Cpool_reserve(kidx0:kidx1,:),                                        &
                                           disturbance%carbon_2_GreenLitterPools(kidx0:kidx1),                       &
                                           disturbance%carbon_2_WoodLitterPools(kidx0:kidx1),                        &
                                           cbal%Cpool_litter_green_ag(kidx0:kidx1,:),                                &
                                           cbal%Cpool_litter_green_bg(kidx0:kidx1,:),                                &
                                           cbal%Cpool_litter_wood_ag(kidx0:kidx1,:),                                 &
                                           cbal%Cpool_litter_wood_bg(kidx0:kidx1,:),                                 &
                                           cbal%YCpool_acid_ag1(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_acid_ag2(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_water_ag1(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_water_ag2(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_ethanol_ag1(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_ethanol_ag2(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:),                                &
                                           cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:),                                &
                                           cbal%YCpool_acid_bg1(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_acid_bg2(kidx0:kidx1,:),                                      &
                                           cbal%YCpool_water_bg1(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_water_bg2(kidx0:kidx1,:),                                     &
                                           cbal%YCpool_ethanol_bg1(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_ethanol_bg2(kidx0:kidx1,:),                                   &
                                           cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:),                                &
                                           cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:),                                &
                                           cbal%YCpool_humus_1(kidx0:kidx1,:),                                       & 
                                           cbal%YCpool_humus_2(kidx0:kidx1,:),                                       &
                                           LeafLit_coef(1:nidx,:,:),                                                 &
                                           WoodLit_coef(1:nidx,:,:))

            ! -- changes in carbon pools caused by fire
            CALL relocate_disturbed_carbon(lctlib, surf, nidx, kidx0, ntiles, with_yasso, DIST_FIRE,           &
                                           disturbance%burned_frac(kidx0:kidx1,:),                             &
                                           surf%cover_fract(kidx0:kidx1,:), veg_fract_correction(:,:),         &
                                           cbal%Cpool_green(kidx0:kidx1,:),                                    &
                                           cbal%Cpool_woods(kidx0:kidx1,:), cbal%Cpool_reserve(kidx0:kidx1,:), &
                                           disturbance%carbon_2_GreenLitterPools(kidx0:kidx1),                 &
                                           disturbance%carbon_2_WoodLitterPools(kidx0:kidx1),                  &
                                           cbal%Cpool_litter_green_ag(kidx0:kidx1,:),                          &
                                           cbal%Cpool_litter_green_bg(kidx0:kidx1,:),                          &
                                           cbal%Cpool_litter_wood_ag(kidx0:kidx1,:),                           &
                                           cbal%Cpool_litter_wood_bg(kidx0:kidx1,:),                           &
                  LeafLit_coef           =  LeafLit_coef(1:nidx,:,:),                                          &
                  WoodLit_coef           =  WoodLit_coef(1:nidx,:,:),                                          &
                  YCpool_acid_ag1        =  cbal%YCpool_acid_ag1(kidx0:kidx1,:),                               &
                  YCpool_water_ag1       =  cbal%YCpool_water_ag1(kidx0:kidx1,:),                              &
                  YCpool_ethanol_ag1     =  cbal%YCpool_ethanol_ag1(kidx0:kidx1,:),                            &
                  YCpool_nonsoluble_ag1  =  cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:),                         &
                  YCpool_acid_bg1        =  cbal%YCpool_acid_bg1(kidx0:kidx1,:),                               &
                  YCpool_water_bg1       =  cbal%YCpool_water_bg1(kidx0:kidx1,:),                              &
                  YCpool_ethanol_bg1     =  cbal%YCpool_ethanol_bg1(kidx0:kidx1,:),                            &
                  YCpool_nonsoluble_bg1  =  cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:),                         &
                  YCpool_humus_1         =  cbal%YCpool_humus_1(kidx0:kidx1,:),                                & 
                  YCpool_acid_ag2        =  cbal%YCpool_acid_ag2(kidx0:kidx1,:),                               &
                  YCpool_water_ag2       =  cbal%YCpool_water_ag2(kidx0:kidx1,:),                              &
                  YCpool_ethanol_ag2     =  cbal%YCpool_ethanol_ag2(kidx0:kidx1,:),                            &
                  YCpool_nonsoluble_ag2  =  cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:),                         &
                  YCpool_acid_bg2        =  cbal%YCpool_acid_bg2(kidx0:kidx1,:),                               &
                  YCpool_water_bg2       =  cbal%YCpool_water_bg2(kidx0:kidx1,:),                              &
                  YCpool_ethanol_bg2     =  cbal%YCpool_ethanol_bg2(kidx0:kidx1,:),                            &
                  YCpool_nonsoluble_bg2  =  cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:),                         &
                  YCpool_humus_2         =  cbal%YCpool_humus_2(kidx0:kidx1,:),                                &
                  carbon_2_atmos         =  disturbance%carbon_2_atmos(kidx0:kidx1))
         END IF  ! Nitrogen1 

         IF (debug_Cconservation) THEN
            ! -- finish third carbon conservation test
            IF (.NOT. with_yasso) THEN
               dynveg%dynveg_testCconserv_2(kidx0:kidx1) = dynveg%dynveg_testCconserv_2(kidx0:kidx1) &
                  - disturbance%carbon_2_atmos(kidx0:kidx1)                                          &
                  - SUM( surf%cover_fract(kidx0:kidx1,:) * veg_fract_correction(:,:) *               &
                         (  cbal%Cpool_green(kidx0:kidx1,:) + cbal%Cpool_woods(kidx0:kidx1,:)        &
                          + cbal%Cpool_reserve(kidx0:kidx1,:)  &
                          + cbal%Cpool_litter_green_ag(kidx0:kidx1,:) + cbal%Cpool_litter_green_bg(kidx0:kidx1,:) &  
                          + cbal%Cpool_litter_wood_ag(kidx0:kidx1,:) + cbal%Cpool_litter_wood_bg(kidx0:kidx1,:)   &
                          + cbal%Cpool_slow(kidx0:kidx1,:) &
                         ),DIM=2 )
            ELSE
               dynveg%dynveg_testCconserv_2(kidx0:kidx1) = dynveg%dynveg_testCconserv_2(kidx0:kidx1) &
                  - disturbance%carbon_2_atmos(kidx0:kidx1)                                          &
                  - SUM( surf%cover_fract(kidx0:kidx1,:) * veg_fract_correction(:,:) *               &
                         (  cbal%Cpool_green(kidx0:kidx1,:) + cbal%Cpool_woods(kidx0:kidx1,:)        &
                          + cbal%Cpool_reserve(kidx0:kidx1,:)                                        &
                          + cbal%YCpool_acid_ag1(kidx0:kidx1,:)                                       &
                          + cbal%YCpool_water_ag1(kidx0:kidx1,:)                                      &
                          + cbal%YCpool_ethanol_ag1(kidx0:kidx1,:)                                    &
                          + cbal%YCpool_nonsoluble_ag1(kidx0:kidx1,:)                                 &
                          + cbal%YCpool_acid_bg1(kidx0:kidx1,:)                                       &
                          + cbal%YCpool_water_bg1(kidx0:kidx1,:)                                      &
                          + cbal%YCpool_ethanol_bg1(kidx0:kidx1,:)                                    &
                          + cbal%YCpool_nonsoluble_bg1(kidx0:kidx1,:)                                 &
                          + cbal%YCpool_humus_1(kidx0:kidx1,:)                                         &
                          + cbal%YCpool_acid_ag2(kidx0:kidx1,:)                                       &
                          + cbal%YCpool_water_ag2(kidx0:kidx1,:)                                      &
                          + cbal%YCpool_ethanol_ag2(kidx0:kidx1,:)                                    &
                          + cbal%YCpool_nonsoluble_ag2(kidx0:kidx1,:)                                 &
                          + cbal%YCpool_acid_bg2(kidx0:kidx1,:)                                       &
                          + cbal%YCpool_water_bg2(kidx0:kidx1,:)                                      &
                          + cbal%YCpool_ethanol_bg2(kidx0:kidx1,:)                                    &
                          + cbal%YCpool_nonsoluble_bg2(kidx0:kidx1,:)                                 &
                          + cbal%YCpool_humus_2(kidx0:kidx1,:)                                         &
                         ),DIM=2 )
            END IF
         END IF

         ! -- scale fluxes from vegetated area to whole grid box
         disturbance%carbon_2_GreenLitterPools(kidx0:kidx1) = disturbance%carbon_2_GreenLitterPools(kidx0:kidx1) * &
                                                              surf%veg_ratio_max(kidx0:kidx1)
         disturbance%carbon_2_WoodLitterPools(kidx0:kidx1)  = disturbance%carbon_2_WoodLitterPools(kidx0:kidx1)  * &
                                                              surf%veg_ratio_max(kidx0:kidx1)
         disturbance%carbon_2_atmos(kidx0:kidx1) = disturbance%carbon_2_atmos(kidx0:kidx1) * surf%veg_ratio_max(kidx0:kidx1)

         IF (with_nitrogen) THEN
           disturbance%nitrogen_2_GreenLitterPools(kidx0:kidx1) = disturbance%nitrogen_2_GreenLitterPools(kidx0:kidx1) * &
                                                              surf%veg_ratio_max(kidx0:kidx1)
           disturbance%nitrogen_2_WoodLitterPools(kidx0:kidx1)  = disturbance%nitrogen_2_WoodLitterPools(kidx0:kidx1)  * &
                                                              surf%veg_ratio_max(kidx0:kidx1)
           disturbance%nitrogen_2_atmos(kidx0:kidx1) = disturbance%nitrogen_2_atmos(kidx0:kidx1) * surf%veg_ratio_max(kidx0:kidx1)
           disturbance%nitrogen_2_sminn(kidx0:kidx1) = disturbance%nitrogen_2_sminn(kidx0:kidx1) * surf%veg_ratio_max(kidx0:kidx1) 
         END IF

         ! Cumulate CO2 emission for flux output
         disturbance%box_CO2_flux_2_atmos(kidx0:kidx1) = disturbance%box_CO2_flux_2_atmos(kidx0:kidx1) + &
                                                         disturbance%carbon_2_atmos(kidx0:kidx1) * molarMassCO2_kg
       ENDIF

    ENDIF

    !-- Calculate CO2 flux to the atmosphere (updated once a day)
    IF (use_disturbance) THEN
       CO2_emission(:) = disturbance%carbon_2_atmos(kidx0:kidx1) * molarMassCO2_kg/86400._dp
    ELSE
       CO2_emission(:) = 0._dp
    ENDIF

  END SUBROUTINE update_dynveg

!------------------------------------------------------------------------------
!
! !IROUTINE potential_tree_fpc
!
! !SUBROUTINE INTERFACE:

  SUBROUTINE potential_tree_fpc(nidx, ntiles, dynveg_all, woody_pft, dynamic_pft, npp_ave, &
                                bio_exist, act_fpc, cover_type, pot_fpc)

! !DESCRIPTION:
!
! subroutine calculates potential FPC (in absence of disturbances) based on NPP for each PFT
!
! called once a year!
!------------------------------------------------------------------------------

! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: dynveg_all      ! flag to activate competition between woody types and grass
    LOGICAL,  INTENT(in)  :: woody_pft(:)
    LOGICAL,  INTENT(in)  :: dynamic_pft(:)
    REAL(dp), INTENT(in)  :: bio_exist(:,:)
    REAL(dp), INTENT(in)  :: npp_ave(:,:)    ! NPP averaged
    REAL(dp), INTENT(in)  :: act_fpc(:,:)
    INTEGER,  INTENT(in)  :: cover_type(:,:) ! cover types

! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(out) :: pot_fpc(:,:)

! !LOCAL VARIABLES:
!
    INTEGER           :: i, itile
    REAL(dp)          :: sum_npp(nidx)
      
!------------------------------------------------------------------------------
!
!-- summing up NPP
!

    sum_npp(:) = 0._dp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (dynveg_all) THEN
             IF (dynamic_pft(cover_type(i,itile)) .AND. npp_ave(i,itile) > sum_npp_min .AND. bio_exist(i,itile) > 0.5_dp) &
                 sum_npp(i) = sum_npp(i) + npp_ave(i,itile) ** npp_nonlinearity * MAX(act_fpc_min,act_fpc(i,itile))
          ELSE
             IF (woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile)) .AND. &
              npp_ave(i,itile) > sum_npp_min .AND. bio_exist(i,itile) > 0.5_dp) &
              sum_npp(i) = sum_npp(i) + npp_ave(i,itile) ** npp_nonlinearity * MAX(act_fpc_min,act_fpc(i,itile))
          END IF
       END DO
    END DO
!
!-- determine potential FPC
!
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (dynveg_all) THEN
             IF (dynamic_pft(cover_type(i,itile)) .AND. npp_ave(i,itile) > sum_npp_min .AND. &
                 sum_npp(i) > sum_npp_min .AND. bio_exist(i,itile) > 0.5_dp) THEN
                pot_fpc(i,itile) = npp_ave(i,itile) ** npp_nonlinearity / sum_npp(i) * MAX(act_fpc_min,act_fpc(i,itile))
             ELSE
                pot_fpc(i,itile) = 0._dp
             END IF
          ELSE
             IF (woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile)) .AND. npp_ave(i,itile) > sum_npp_min .AND. &
                 bio_exist(i,itile) > 0.5_dp) THEN
                pot_fpc(i,itile) = npp_ave(i,itile) ** npp_nonlinearity / sum_npp(i) * &
                                 tree_fpc_max * MAX(act_fpc_min,act_fpc(i,itile))
             ELSE
                pot_fpc(i,itile) = 0._dp
             END IF
          END IF
       END DO
    END DO

  END SUBROUTINE potential_tree_fpc

!------------------------------------------------------------------------------
!
! !IROUTINE desert_fraction
!
! !SUBROUTINE INTERFACE:

  SUBROUTINE desert_fraction(nidx, ntiles, nlct, accelerate_dynveg,  &
                             woody_pft, dynamic_pft,                 &
                             act_fpc, bare_fpc, sla,                 &
                             max_green_bio, cover_type,              &
                             sum_green_bio_memory, desert_fpc)
! !DESCRIPTION:
!
! Calculation of desert fraction from green biomass
!
!------------------------------------------------------------------------------


! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    INTEGER,  INTENT(in)  :: nlct            ! number of land cover types
    REAL(dp), INTENT(in)  :: accelerate_dynveg  ! acceleration factor for vegetation and desert dynamics
    LOGICAL,  INTENT(in)  :: woody_pft(:)    ! flag to label a woody PFT
    LOGICAL,  INTENT(in)  :: dynamic_pft(:)  ! flag to label PFTs that take part in the natural vegetation dynamics
    REAL(dp), INTENT(in)  :: act_fpc(:,:)    ! actual FPC
    REAL(dp), INTENT(in)  :: bare_fpc(:)     ! bare FPC
    REAL(dp), INTENT(in)  :: sla(:,:)        ! specific leaf area
    REAL(dp), INTENT(in)  :: max_green_bio(:,:) ! maximum value of green biomass within a year
    INTEGER,  INTENT(in)  :: cover_type(:,:) ! cover types

! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(inout) :: sum_green_bio_memory(:) ! vegetated fraction calculated from green biomass

! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(out) :: desert_fpc(:) ! desert FPC

! !LOCAL VARIABLES:
!
    INTEGER           :: i, itile
    REAL(dp)          :: sum_green_bio(nidx)
    REAL(dp)          :: sum_act_fpc(nidx)
    REAL(dp)          :: sum_grass_fpc(nidx)
    INTEGER           :: n_grass_pft
      
!------------------------------------------------------------------------------

    sum_green_bio(:) = 0._dp
    sum_act_fpc(:) = 0._dp
    sum_grass_fpc(:) = 0._dp
    n_grass_pft = 0

    DO i = 1,nlct
       IF (.NOT. woody_pft(i) .AND. dynamic_pft(i)) n_grass_pft = n_grass_pft + 1
    END DO

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (.NOT. woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) THEN
             sum_grass_fpc(i) = sum_grass_fpc(i) + act_fpc(i,itile)
          END IF
       END DO
    END DO

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) THEN
             sum_green_bio(i) = sum_green_bio(i) + MAX(0._dp, act_fpc(i,itile) * (1._dp -          &
               exp(-desert_extend * ((max_green_bio(i,itile) * sla(i,itile)) ** desert_margin) /   &
               (3._dp ** (desert_margin - 1._dp)))))
             sum_act_fpc(i) = sum_act_fpc(i) + act_fpc(i,itile)
          ELSE IF (.NOT. woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) THEN
             IF (sum_grass_fpc(i) > EPSILON(1._dp)) THEN
                sum_green_bio(i) = sum_green_bio(i) +                                              &
                 MAX(0._dp, act_fpc(i,itile) * (1._dp + bare_fpc(i) / sum_grass_fpc(i)) * (1._dp - &
                 exp(-desert_extend * ((max_green_bio(i,itile) * sla(i,itile)) ** desert_margin) / &
                (3._dp ** (desert_margin - 1._dp)))))
                sum_act_fpc(i) = sum_act_fpc(i) + act_fpc(i,itile) * (1._dp + bare_fpc(i) / sum_grass_fpc(i))
             ELSE
                sum_green_bio(i) = sum_green_bio(i) +                                              &
                 MAX(0._dp, (bare_fpc(i) / REAL(n_grass_pft)) * (1._dp -                           &
                 exp(-desert_extend * ((max_green_bio(i,itile) * sla(i,itile)) ** desert_margin) / &
                 (3._dp ** (desert_margin - 1._dp)))))
                sum_act_fpc(i) = sum_act_fpc(i) + (bare_fpc(i) / REAL(n_grass_pft))
             END IF
          END IF
       END DO
    END DO

    WHERE (sum_act_fpc(:) > EPSILON(1._dp)) 
       sum_green_bio(:) = sum_green_bio(:) / sum_act_fpc(:) 
    END WHERE
    sum_green_bio_memory(:) = (sum_green_bio_memory(:) * (tau_desert/accelerate_dynveg - 1._dp) + sum_green_bio(:)) &
                               / (tau_desert/accelerate_dynveg)

    desert_fpc(:) = MIN(1._dp - 2._dp * EPSILON(1._dp),MAX(0._dp,1._dp - sum_green_bio_memory(:)))

  END SUBROUTINE desert_fraction

!------------------------------------------------------------------------------
!
! !IROUTINE: fpc_daily
!
! !SUBROUTINE INTERFACE:
!
  SUBROUTINE fpc_daily(lctlib, nidx, kidx0, kidx1, ntiles, with_yasso, dynveg_all, accelerate_dynveg, &
                       woody_pft, dynamic_pft, tau_Cpool_woods,                                       &
                       npp_ave, bio_exist, pot_fpc, veg_fract_correction, surf, cbal, climbuf,        &
                       act_fpc_previous_day, act_fpc, bare_fpc,                                       &
                       burned_frac, burned_frac_diag, damaged_frac, fuel)
!
! !DESCRIPTION:
!
! Calculation of daily FPC for trees, shrubs, and herbaceous plants
!
! called: each day
!------------------------------------------------------------------------------
!
    USE mo_jsbach_lctlib,    ONLY: lctlib_type
    USE mo_disturbance,      ONLY: dist_opts
! !INPUT PARAMETERS:
!
    TYPE(lctlib_type),INTENT(in) :: lctlib         !! PFT-specific constants
    INTEGER,         INTENT(in) :: nidx            ! vector length
    INTEGER,         INTENT(in) :: kidx0, kidx1    ! vector limits
    INTEGER,         INTENT(in) :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,         INTENT(in) :: with_yasso
    LOGICAL,         INTENT(in) :: dynveg_all      ! flag to activate competition between woody types and grass
    REAL(dp),        INTENT(in) :: accelerate_dynveg  ! factor to accellerate vegetation and desert dynamics
    LOGICAL,         INTENT(in) :: woody_pft(:)
    LOGICAL,         INTENT(in) :: dynamic_pft(:)
    REAL(dp),        INTENT(in) :: tau_Cpool_woods(:) ! lifetime of Cpool_woods [days]
    REAL(dp),        INTENT(in) :: npp_ave(:,:)
    REAL(dp),        INTENT(in) :: bio_exist(:,:)
    REAL(dp),        INTENT(in) :: pot_fpc(:,:)
    REAL(dp),        INTENT(in) :: veg_fract_correction(:,:)
    TYPE(land_surface_type), INTENT(in)    :: surf
    TYPE(cbalance_type)    , INTENT(inout) :: cbal
    TYPE(climbuf_type)     , INTENT(in)    :: climbuf
    REAL(dp),        INTENT(in) :: act_fpc_previous_day(:,:)

! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(inout) :: act_fpc(:,:)
    REAL(dp), INTENT(inout) :: bare_fpc(:)

! !OUTPUT PARAMETERS:
!
    REAL(dp), target, intent(out), optional :: burned_frac(:,:)
    REAL(dp), target, intent(out), optional :: burned_frac_diag(:,:)
    REAL(dp), target, intent(out), optional :: damaged_frac(:,:)
    REAL(dp), target, intent(inout), optional :: fuel(:)

! !LOCAL VARIABLES:
!  
    INTEGER     :: i, j, itile
    REAL(dp)    :: total_act_fpc(nidx)
    REAL(dp)    :: grass_sum(nidx)
    REAL(dp)    :: non_woody_fpc(nidx)
    REAL(dp)    :: woody_estab_fpc(nidx)
    REAL(dp)    :: theta, veg_inc, tau_pft
    REAL(dp), pointer :: burned_local(:,:)
    REAL(dp), pointer :: damaged_local(:,:)
    integer     :: ct
  

!------------------------------------------------------------------------------
! Initialisations
! ---------------

    IF (PRESENT(burned_frac)) THEN
      burned_local     => burned_frac
      damaged_local    => damaged_frac
      burned_local (:,:) = 0._dp
      damaged_local(:,:) = 0._dp
    ELSE
      burned_local     => dynveg%Dummy_null ! Array (nland,ntiles) of permanent zeros preallocated as
      damaged_local    => dynveg%Dummy_null ! replacement for real disturbance variables in case of no disturbance
    END IF

    non_woody_fpc(:)=bare_fpc(:)

    DO itile = 1,ntiles
       DO i = 1,nidx
          j = i + kidx0 - 1
          IF (.NOT. dynamic_pft(surf%cover_type(j,itile))) act_fpc(i,itile) = 0._dp
          IF (.NOT. woody_pft(surf%cover_type(j,itile))) non_woody_fpc(i) = non_woody_fpc(i) + act_fpc(i,itile)
       END DO
    END DO

    !-- adjust grass fractions to changes in woody type fpc
    IF (dynveg_all) THEN

      !-- reset temporary arrays, only used when dynveg_all = .true.
      woody_estab_fpc(:)=0._dp

      ! fraction to convert
      DO itile= 1,ntiles
        DO i = 1,nidx
          ct = surf%cover_type(i+kidx0-1,itile)
          tau_pft = tau_Cpool_woods(ct) * tau_Cpool_woods_to_tau_pft
          IF (dynamic_pft(ct) .AND. woody_pft(ct)) THEN
            woody_estab_fpc(i) = woody_estab_fpc(i) + accelerate_dynveg * pot_fpc(i,itile) * &
                                            non_woody_fpc(i) / tau_pft
          ENDIF
        ENDDO
      ENDDO

      ! do the conversion
      DO itile= 1,ntiles
        DO i = 1,nidx
          ct = surf%cover_type(i+kidx0-1,itile)
          IF (dynamic_pft(ct) .AND. .NOT. woody_pft(ct) .AND. non_woody_fpc(i) > EPSILON(1._dp)) THEN
                 act_fpc(i,itile) = MAX(0._dp,act_fpc(i,itile) * (1._dp - woody_estab_fpc(i) / non_woody_fpc(i)))
          ENDIF
        ENDDO
      ENDDO

    ELSE ! dynveg_all

      !-- reset temporary arrays, only used when dynveg_all = .false.
      grass_sum(:)     = 0._dp

      !-- temporary sums for grass types, dynveg_all = .false.
      DO itile = 1,ntiles
        DO i = 1,nidx
          ct=surf%cover_type(i+kidx0-1,itile)
          IF (.NOT. woody_pft(ct) .AND. dynamic_pft(ct) .AND. bio_exist(i,itile) > 0.5_dp .AND. &
              npp_ave(i,itile) > EPSILON(1._dp)) THEN
            grass_sum(i) = grass_sum(i) + npp_ave(i,itile) * MAX(act_fpc_min,act_fpc(i,itile))
          END IF
        END DO
      END DO

    ENDIF ! .not. dynveg_all

    !-- area wind break of woody types (trees and shrubs) and burned
    IF (PRESENT(burned_frac)) THEN
      CALL disturbed_frac(lctlib, nidx,kidx0,kidx1,ntiles,DIST_FIRE, &
                          dynamic_pft(:) .OR. (dist_opts%lburn_pasture .AND. lctlib%pasture_pft(:)), with_yasso, &
                          act_fpc,veg_fract_correction,surf,cbal,climbuf,burned_frac,burned_frac_diag,fuel)
      CALL disturbed_frac(lctlib, nidx,kidx0,kidx1,ntiles,DIST_WINDBREAK,woody_pft(:) .AND. dynamic_pft(:), with_yasso, &
                          act_fpc,veg_fract_correction,surf,cbal,climbuf,damaged_frac) 
    END IF

    !-- for dynveg_all, woody and grass type loops can be sensibly joined
    IF (dynveg_all) THEN
      DO itile = 1,ntiles
        DO i = 1,nidx
          ct = surf%cover_type(i+kidx0-1,itile)
          tau_pft = tau_Cpool_woods(ct) * tau_Cpool_woods_to_tau_pft
          IF (dynamic_pft(ct)) THEN
            IF (woody_pft(ct)) THEN
              ! woody types
              act_fpc(i,itile) = MAX(0._dp,act_fpc(i,itile) + accelerate_dynveg * &
                ((pot_fpc(i,itile) * non_woody_fpc(i) - act_fpc(i,itile) / life_to_estab) / tau_pft &
                 - (burned_local(i,itile) + damaged_local(i,itile)) * act_fpc_previous_day(i,itile)))
            ELSE
              ! grass types
              IF (bio_exist(i,itile) > 0.5_dp .AND. npp_ave(i,itile) > EPSILON(1._dp)) THEN
                act_fpc(i,itile) = MAX(0._dp,act_fpc(i,itile) + accelerate_dynveg * &
                  ((pot_fpc(i,itile) * bare_fpc(i)) / tau_pft - burned_local(i,itile) * act_fpc_previous_day(i,itile)))
              ELSE
                act_fpc(i,itile) = MAX(0._dp,act_fpc(i,itile) + accelerate_dynveg * &
                  ((pot_fpc(i,itile) * bare_fpc(i) - act_fpc(i,itile) / life_to_estab) / tau_pft &
                  - burned_local(i,itile) * act_fpc_previous_day(i,itile)))
              ENDIF
            ENDIF ! .not. woody_pft
          ENDIF ! dynamic_pft
        ENDDO
      ENDDO

    ELSE ! dynveg_all

      !-- dynamic equation for act_fpc, woody types, dynveg_all = .false.
      DO i = 1,nidx
        veg_inc = 0._dp
        DO itile = 1,ntiles
          ct = surf%cover_type(i+kidx0-1,itile)
          tau_pft = tau_Cpool_woods(ct) * tau_Cpool_woods_to_tau_pft
          IF (dynamic_pft(ct) .AND. woody_pft(ct)) THEN
            veg_inc = veg_inc + MAX(fract_small,act_fpc(i,itile) + accelerate_dynveg * &
            ((pot_fpc(i,itile) - act_fpc(i,itile)) / tau_pft &
            - (burned_local(i,itile) + damaged_local(i,itile)) * act_fpc_previous_day(i,itile))) &
            - act_fpc(i,itile)
          ELSEIF(dynamic_pft(ct) .AND. .NOT. woody_pft(ct)) THEN
            IF (bio_exist(i,itile) > 0.5_dp .AND. npp_ave(i,itile) > EPSILON(1._dp)) THEN
              veg_inc = veg_inc + MAX(fract_small,act_fpc(i,itile) - accelerate_dynveg * &
              (burned_local(i,itile) * act_fpc_previous_day(i,itile) + &
              grass_mortality * act_fpc(i,itile) - &
              bare_fpc(i) / grass_sum(i) * npp_ave(i,itile) * MAX(act_fpc_min,act_fpc(i,itile)) / tau_pft)) &
              - act_fpc(i,itile)
            ELSE
              veg_inc = veg_inc + MAX(fract_small,act_fpc(i,itile) - accelerate_dynveg * &
              (burned_local(i,itile) * act_fpc_previous_day(i,itile) + &
              act_fpc(i,itile) / tau_pft)) &
              - act_fpc(i,itile)
            ENDIF
          ENDIF
        END DO
        theta = 1._dp
        IF (bare_fpc(i) < (veg_inc + fract_small)) theta = 0._dp
        DO itile = 1,ntiles
          ct = surf%cover_type(i+kidx0-1,itile)
          tau_pft = tau_Cpool_woods(ct) * tau_Cpool_woods_to_tau_pft
          IF (dynamic_pft(ct) .AND. woody_pft(ct)) THEN
            act_fpc(i,itile) = MAX(fract_small,act_fpc(i,itile) + accelerate_dynveg * &
              ((theta * pot_fpc(i,itile) - act_fpc(i,itile)) / tau_pft                &
              - (burned_local(i,itile) + damaged_local(i,itile)) * act_fpc_previous_day(i,itile)))
          ENDIF
        ENDDO
      ENDDO

      !-- dynamic equation for act_fpc, grass types, dynveg_all = .false.
      DO i = 1,nidx
        DO itile = 1,ntiles
          ct = surf%cover_type(i+kidx0-1,itile)
          tau_pft = tau_Cpool_woods(ct) * tau_Cpool_woods_to_tau_pft
          IF (dynamic_pft(ct) .AND. .NOT. woody_pft(ct)) THEN
            IF (bio_exist(i,itile) > 0.5_dp .AND. npp_ave(i,itile) > EPSILON(1._dp)) THEN
              act_fpc(i,itile) = MAX(fract_small,act_fpc(i,itile) - accelerate_dynveg * &
                 (burned_local(i,itile) * act_fpc_previous_day(i,itile) + &
                 grass_mortality * act_fpc(i,itile) - &
                 bare_fpc(i) / grass_sum(i) * npp_ave(i,itile) * MAX(act_fpc_min,act_fpc(i,itile)) / tau_pft))
            ELSE
              act_fpc(i,itile) = MAX(fract_small,act_fpc(i,itile) - accelerate_dynveg * &
                 (burned_local(i,itile) * act_fpc_previous_day(i,itile) + &
                 act_fpc(i,itile) / tau_pft))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF ! .not. dynveg_all

    !-- calculate total area fraction covered by dynamic vegetation

    total_act_fpc(:)=0._dp

    DO itile = 1,ntiles
       DO i = 1,nidx
          IF (dynamic_pft(surf%cover_type(i+kidx0-1,itile))) total_act_fpc(i) = total_act_fpc(i) + act_fpc(i,itile)
       END DO
    END DO

    !-- calculation of bare soil fraction
    
    bare_fpc(:) = MAX(0._dp,1._dp - total_act_fpc(:))

  END SUBROUTINE fpc_daily

!------------------------------------------------------------------------------
!
! !IROUTINE: bioclim_limits
!
! !SUBROUTINE INTERFACE:
!
  SUBROUTINE bioclim_limits (nidx, ntiles, lctlib, min_mmtemp20, max_mmtemp20, prev_year_gdd, cover_type, bio_exist)
!
! !DESCRIPTION:
!
! calculation of bioclimatic limits for each PFT based on updated
! LPJ lookup-table
!
! Limits based on 20-year running averages of coldest-month mean
! temperature and growing degree days (5 degree base, except larches
! (2 degrees)), and minimal temperature range required (larches).
! For SURVIVAL, coldest month temperature and GDD should be
! at least as high as PFT-specific limits.
! For REGENERATION, PFT must be able to survive AND coldest month
! temperature should be no higher than a PFT-specific limit.
!------------------------------------------------------------------------------
    USE mo_jsbach_lctlib, ONLY: lctlib_type
 
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx             ! vector length
    INTEGER,  INTENT(in)  :: ntiles           ! number of tiles (different land cover type mosaics per gridcell)
    TYPE(lctlib_type), INTENT(in) :: lctlib   ! parameters of the land cover type library
    REAL(dp), INTENT(in)  :: min_mmtemp20(:)  ! Minimum monthly mean temp. (20yr climatology)
    REAL(dp), INTENT(in)  :: max_mmtemp20(:)  ! Maximum monthly mean temp. (20yr climatology)
    REAL(dp), INTENT(in)  :: prev_year_gdd(:,:) ! GDD of previous year
    INTEGER,  INTENT(in)  :: cover_type(:,:)  ! land cover types

! ! OUTPUT PARAMETERS: 
!
    REAL(dp), INTENT(out) :: bio_exist(:,:)

! !LOCAL VARIABLES:
!
    INTEGER   :: i, itile
    REAL(dp)  :: tcmin         !PFT-specific minimum coldest-month temperature limit
    REAL(dp)  :: tcmax         !PFT-specific maximum coldest-month temperature limit
    REAL(dp)  :: gddmin        !PFT-specific minimum GDD
    REAL(dp)  :: twmax         !PFT-specific maximum warmest-month temperature limit
    REAL(dp)  :: min_temprange !PFT-specific minimum difference of 20-year average warmest 
                                  !                    minus coldest month temperature
!------------------------------------------------------------------------------

    bio_exist(:,:) = 1._dp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (lctlib%dynamic_pft(cover_type(i,itile))) THEN
             tcmin=lctlib%bclimit_min_cold_mmtemp(cover_type(i,itile))
             tcmax=lctlib%bclimit_max_cold_mmtemp(cover_type(i,itile))
             twmax=lctlib%bclimit_max_warm_mmtemp(cover_type(i,itile))
             min_temprange=lctlib%bclimit_min_temprange(cover_type(i,itile))
             gddmin=lctlib%bclimit_min_gdd(cover_type(i,itile))

             bio_exist(i,itile) = 1._dp

             IF (min_mmtemp20(i) < tcmin .OR. min_mmtemp20(i) > tcmax              &
               .OR. prev_year_gdd(i,itile) < gddmin .OR. max_mmtemp20(i) > twmax   &
               .OR. (max_mmtemp20(i) - min_mmtemp20(i)) < min_temprange) bio_exist(i,itile) = 0._dp
          END IF
       END DO
    END DO
      
  END SUBROUTINE bioclim_limits

!------------------------------------------------------------------------------
!
! !IROUTINE: cover_fract_pot_to_cover_fract
!
! !SUBROUTINE INTERFACE:
!
  SUBROUTINE cover_fract_pot_to_cover_fract (nidx, ntiles, nlct, is_present, is_glacier, is_naturalveg, &
                      woody_pft, dynamic_pft, pasture_pft, is_crop, cover_type, &
                      cover_fract_pot, cover_fract_pot_previous_year, cover_fract)
!
! !DESCRIPTION:
!
! convert the fractional plant cover of potential natural vegetation to plant cover fractions including agricultural areas.
! Pasture is first established on grasslands, whereas crops are established at the expense of all natural cover types.
!
! In runs with land use transitions, the transitions from natural forests or grass lands to crops and pastures are predefined.
! It is not wanted, that the dynamic vegetation interferes with these transitions by e.g. re-establishing forests where they
! had just been removed. Thus the changes in forest fraction calculated by the dynamic vegetation are not instantanously 
! applied. This results in an imbalance between actual cover_fractions and the ones the dynamic vegetation calculates for
! instantanous establishement of crops and pasture (cover_fract_inst).
! The final cover_fractions are based on cover_fract_inst with the exception, that the forest fraction is based on the forest
! fraction after land use transitions (cover_fract entering the routine), and only the natural change in cover_fractions in 
! the direction to reduce the inbalance is taken into account (see below).
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx                 ! vector length
    INTEGER,  INTENT(in)  :: ntiles               ! number of tiles (different land cover type mosaics per gridcell)
    INTEGER,  INTENT(in)  :: nlct                 ! number of land cover types (lctlib)
    LOGICAL,  INTENT(in)  :: is_present(:,:)      ! logical mask for cells treated by jsbach
    LOGICAL,  INTENT(in)  :: is_glacier(:,:)      ! logical glacier mask on tiles
    LOGICAL,  INTENT(in)  :: is_naturalveg(:,:)   ! PFTs that are natural vegetation
    LOGICAL,  INTENT(in)  :: woody_pft(:)         ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: dynamic_pft(:)       ! PFTs taking part in competition
    LOGICAL,  INTENT(in)  :: pasture_pft(:)       ! PFTs that are pasture
    LOGICAL,  INTENT(in)  :: is_crop(:,:)         ! logical mask for crops
    INTEGER,  INTENT(in)  :: cover_type(:,:)      ! land cover types
    REAL(dp), INTENT(in)  :: cover_fract_pot(:,:) ! cover fractions if there was only natural vegetation
    REAL(dp), INTENT(in)  :: cover_fract_pot_previous_year(:,:) ! last years cover_fract_pot
!
! !IN- and OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(inout) :: cover_fract(:,:)  ! cover fraction within jsbach


! !LOCAL VARIABLES:
!
    INTEGER   :: i, itile
    REAL(dp)  :: nwoody, woody(nlct)            ! number of woody PFTs (as REAL) and helper array
    REAL(dp)  :: nexcluded, dyn(nlct)           ! number of non-dynamic PFTs (crops, pastures) and helper array
    REAL(dp)  :: cover_fract_inst(nidx,ntiles)  ! cover fractions that would result from an instantanous implementation of 
                                                ! land use change and vegetation dynamics
    REAL(dp)  :: excluded_fract(nidx)           ! fraction of the vegetated part of the grid box not considered for vegetation
                                                ! dynamics (e.g. agricultural areas)
    REAL(dp)  :: sum_woody_fract(nidx)          ! fraction of the vegetated part of the grid box covered by woody plants
    REAL(dp)  :: sum_grass_fract(nidx)          ! fraction of the vegetated part of the grid box covered by grass
    REAL(dp)  :: sum_woody_fract_inst(nidx)     ! fraction of the vegetated part of the grid box covered by woody plants
                                                ! based on cover_fract_inst
    REAL(dp)  :: sum_woody_fract_pot(nidx)      ! woody fraction of the vegetated part of the grid box considering natural
                                                ! vegetation, only
    REAL(dp)  :: sum_grass_fract_pot(nidx)      ! grass fraction of the vegetated part of the grid box considering natural
                                                ! vegetation, only
    REAL(dp)  :: sum_woody_fract_old(nidx)      ! fraction of the vegetated part of the grid box covered by woody plants
                                                ! (based on act_fpc_previous_year) that would have resulted from an
                                                ! instantaneous establishment of crop land and pasture
    REAL(dp)  :: sum_woody_fract_act(nidx)      ! fraction of the vegetated part of the grid box covered by woody plants
    REAL(dp)  :: sum_grass_fract_act(nidx)      ! fraction of the vegetated part of the grid box covered by grasses
    REAL(dp)  :: delta_woody(nidx)              ! deviation of the fraction of the vegetated part of the grid box covered 
                                                ! by woody plants from the fraction that would result from an instantaneous
                                                !  establishment of crop land and pasture
    LOGICAL   :: glacier(nidx)           ! 2 dim. glacier mask
!------------------------------------------------------------------------------

    !-- define a 2 dimensional glacier mask

    glacier(:) = ANY(is_glacier(:,:),DIM=2)

    !-- find out number of woody types and number of exculded types

    woody(:) = 0._dp
    WHERE (dynamic_pft .AND. woody_pft) woody(:) = 1._dp
    nwoody = SUM(woody(:))
    dyn(:) = 0._dp
    WHERE (dynamic_pft) dyn(:) = 1._dp
    nexcluded = REAL(ntiles) - SUM(dyn(:))

    !-- Sum up cover fractions and of woody types, grasses, pasture etc. The fractions (cover_fract) were last updated by 
    !   last years land use change, and do not yet include current vegetation dynamics.

    sum_woody_fract(:) = 0._dp
    sum_grass_fract(:) = 0._dp
    excluded_fract(:) = 0._dp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) THEN
             sum_woody_fract(i) = sum_woody_fract(i) + cover_fract(i,itile)
          ELSE IF (.NOT. woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) THEN
             sum_grass_fract(i) = sum_grass_fract(i) + cover_fract(i,itile)
          ELSE IF (.NOT. dynamic_pft(cover_type(i,itile))) THEN
             excluded_fract(i) = excluded_fract(i) + cover_fract(i,itile)
          END IF
       END DO
    END DO

    !-- Sum up woody and grass fractions of the potential natural vegetation based on current vegetation dynamics

    sum_woody_fract_pot(:) = 0._dp
    sum_grass_fract_pot(:) = 0._dp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) THEN
             sum_woody_fract_pot(i) = sum_woody_fract_pot(i) + cover_fract_pot(i,itile)
          ELSE IF (.NOT. woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) THEN
             sum_grass_fract_pot(i) = sum_grass_fract_pot(i) + cover_fract_pot(i,itile)
          END IF
       END DO
    END DO

    !-- Calculate the woody fraction resulting from last years potential vegetation assuming an instantanous 
    !   establishment of crops and pastures

    CALL calc_cover_fract_inst(nidx, ntiles, nlct, cover_fract_pot_previous_year(:,:), is_present(:,:), is_glacier(:,:),   &
                               is_naturalveg(:,:), woody_pft(:), dynamic_pft(:), pasture_pft(:), is_crop(:,:), &
                               cover_type(:,:), cover_fract(:,:), cover_fract_inst(:,:))

    sum_woody_fract_old(:) = 0._dp
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) &
               sum_woody_fract_old(i) = sum_woody_fract_old(i) + cover_fract_inst(i,itile)
       END DO
    END DO

    !-- Calculate the woody fraction resulting from current potential vegetation assuming an instantanous 
    !   establishment of crops and pastures

    CALL calc_cover_fract_inst(nidx, ntiles, nlct, cover_fract_pot(:,:), is_present(:,:), is_glacier(:,:),   &
                               is_naturalveg(:,:), woody_pft(:), dynamic_pft(:), pasture_pft(:), is_crop(:,:), &
                               cover_type(:,:), cover_fract(:,:), cover_fract_inst(:,:))

    sum_woody_fract_inst(:) = 0._dp
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) THEN
             sum_woody_fract_inst(i) = sum_woody_fract_inst(i) + cover_fract_inst(i,itile)
          END IF
       END DO
    END DO

    !-- Calculation of the final woody fraction. A difference between sum_woody_fract and sum_woody_fract_old indicates
    !   an imbalance of the state of the dynamic vegetation and the cover fractions calculated from landuse transitions.
    !   On the other hand, the difference between sum_woody_fract_inst and sum_woody_fract_old show the trend of the woody
    !   fraction caused by vegetation dynamics.
    !   The idea of this routine is, to allow shifts in forest fractions only if they reduce the imbalance between dynveg
    !   and landuse change.

    DO i = 1,nidx
       IF (sum_woody_fract(i) - sum_woody_fract_old(i) > 0._dp) THEN

          ! -- woody fraction of cover_fract entering the routine (last updated by last years land use change) is greater
          !    than the the woody fraction resulting from last years potential vegetation and instantanous establishment of
          !    crops and pastures. Thus there is an inbalance, the woody fraction of cover_fract is too high.
 
          IF (sum_woody_fract_inst(i) - sum_woody_fract_old(i) > 0._dp) THEN

             ! -- vegetation dynamics even enlarge woody fraction
             !    ==> no change in actual woody fraction, unless the enlargement by vegetation dynamics exceeds the initial
             !        woody surplus

             delta_woody(i) = MAX(0._dp, sum_woody_fract_inst(i) - sum_woody_fract(i))

          ELSE
             ! -- vegetation dynamics shrink woody fraction

             delta_woody(i) = sum_woody_fract_inst(i) - sum_woody_fract_old(i)

          END IF
       ELSE
          ! -- The woody fraction of cover_fract is too small

          IF (sum_woody_fract_inst(i) - sum_woody_fract_old(i) > 0._dp) THEN

             ! -- vegetation dynamics enlarge woody fraction

             delta_woody(i) = sum_woody_fract_inst(i) - sum_woody_fract_old(i)

          ELSE
             ! -- vegetation dynamics shrink woody fraction
             !    ==> no change in actual woody fraction, unless the reduction by vegetation dynamics exceeds the initial
             !        woody deficit

             delta_woody(i) = MIN(0._dp, sum_woody_fract_inst(i) - sum_woody_fract(i))

          END IF
       END IF

       ! -- in some cases delta_woody can not be used to calculate the actual woody fraction. If there is no grass left
       !    after instantanous establishement of agricultural areas, the woody fractions of cover_fract_inst and
       !    cover_fract_inst_old are identical and delta_woody is zero, even if vegetation dynamics change considerably.
       !    The following formulation is needed to assure that sum_woody_fract_act is smaller than sum_woody_fract_pot and
       !    sum_grass_fract_act is smaller than sum_grass_fract_pot in all cases.
 
       sum_woody_fract_act(i) = MIN(sum_woody_fract_pot(i),&
                                MAX(nwoody*fract_small, sum_woody_fract_pot(i)-excluded_fract(i)+nexcluded*fract_small, &
                                    sum_woody_fract(i) + delta_woody(i)))
        
    END DO

    sum_grass_fract_act(:) = 1._dp - sum_woody_fract_act(:) - excluded_fract(:)


    !-- Use actual woody fraction to determine all actual fractions

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (.NOT. glacier(i) .AND. dynamic_pft(cover_type(i,itile))) THEN
             IF (woody_pft(cover_type(i,itile))) THEN
                cover_fract(i,itile) = cover_fract_pot(i,itile) * sum_woody_fract_act(i) / sum_woody_fract_pot(i)
             ELSE
                cover_fract(i,itile) = cover_fract_pot(i,itile) * sum_grass_fract_act(i) / sum_grass_fract_pot(i)
             END IF
          END IF
       END DO
    END DO

    CALL scale_cover_fract (nidx, ntiles, is_present(:,:), is_glacier(:,:), is_naturalveg(:,:), cover_fract(:,:))

    IF (ANY(SUM(cover_fract(:,:),DIM=2) > 1._dp + ntiles*EPSILON(1._dp)) .OR. &
        ANY(SUM(cover_fract(:,:),DIM=2) < 1._dp - ntiles*EPSILON(1._dp))) THEN
       WRITE (message_text,*) 'sum of cover_fract /= 1: ', &
            MINVAL(SUM(cover_fract(:,:),DIM=2)), MAXVAL(SUM(cover_fract(:,:),DIM=2)), &
            MAXLOC(SUM(cover_fract(:,:),DIM=2))
       CALL finish ('cover_fract_pot_to_cover_fract', message_text)
    END IF
    DO itile = 1,ntiles
       IF (ANY(.NOT. glacier(:) .AND. is_present(:,itile) .AND. cover_fract(:,itile) < fract_small)) THEN
          WRITE (message_text,*) 'cover_fract too small: ', MINVAL(cover_fract(:,itile)), &
               MINLOC(cover_fract(:,itile),MASK=.NOT. glacier(:))
          CALL finish ('cover_fract_pot_to_cover_fract', message_text)
       END IF
    END DO

  END SUBROUTINE cover_fract_pot_to_cover_fract

!------------------------------------------------------------------------------

  SUBROUTINE calc_cover_fract_inst(nidx, ntiles, nlct, cover_fract_pot, is_present, is_glacier, is_naturalveg, &
                                   woody_pft, dynamic_pft, pasture_pft, is_crop, &
                                   cover_type, cover_fract, cover_fract_inst)
!
! !DESCRIPTION:
!
! convert fractional plant cover calculated within dynveg to jsbach cover fractions
! that would result from an instantaneous establishment of crop land and pasture
! (pasture is preferrentially established on grasslands)
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx              ! vector length
    INTEGER,  INTENT(in)  :: ntiles            ! number of tiles (different land cover type mosaics per gridcell)
    INTEGER,  INTENT(in)  :: nlct              ! number of land cover types (from lctlib)
    REAL(dp), INTENT(in)  :: cover_fract_pot(:,:)  ! potential cover fractions (without land use)
    LOGICAL,  INTENT(in)  :: is_present(:,:)   ! logical mask for cells treated by jsbach
    LOGICAL,  INTENT(in)  :: is_glacier(:,:)   ! logical glacier mask on tiles
    LOGICAL,  INTENT(in)  :: is_naturalveg(:,:) ! logical mask of natural vegetation
    LOGICAL,  INTENT(in)  :: woody_pft(:)      ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: dynamic_pft(:)    ! PFTs taking part in competition
    LOGICAL,  INTENT(in)  :: pasture_pft(:)    ! PFTs that are pasture
    LOGICAL,  INTENT(in)  :: is_crop(:,:)      ! logical crop mask
    INTEGER,  INTENT(in)  :: cover_type(:,:)   ! cover types
    REAL(dp), INTENT(in)  :: cover_fract(:,:)  ! cover fraction within jsbach
!
! !OUTPUT PARAMETERS:
! 
    REAL(dp), INTENT(out) :: cover_fract_inst(:,:) ! cover fraction within jsbach resulting from an
                                                   ! instantaneous establishment of crop land and pasture
!
! !LOCAL VARIABLES:
!
    INTEGER   :: i, itile
    REAL(dp)  :: ngrass, grass(nlct)     ! number of grass types (as real) and helper array
    REAL(dp)  :: sum_grass_fract(nidx)   ! potential fraction of grass if there was no land use
    REAL(dp)  :: sum_pasture_fract(nidx) ! fraction of the vegetated part of the grid box covered with pasture
    REAL(dp)  :: sum_crop_fract(nidx)    ! fraction of the vegetated part of the grid box covered with crops
    LOGICAL   :: glacier(nidx)           ! 2 dim. glacier mask


    !-- define a 2 dimensional glacier mask

    glacier(:) = ANY(is_glacier(:,:),DIM=2)

    !-- find out number of grass types

    grass(:) = 0._dp
    WHERE (dynamic_pft .AND. .NOT. woody_pft) grass(:) = 1._dp
    ngrass = SUM(grass(:))

    !-- sum up FPC of grasses, pastures and crops

    sum_grass_fract(:) = 0._dp
    sum_pasture_fract(:) = 0._dp
    sum_crop_fract(:) = 0._dp
    DO i = 1,nidx
      DO itile = 1,ntiles
         IF (.NOT. woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) THEN
            sum_grass_fract(i) = sum_grass_fract(i) + cover_fract_pot(i,itile)
         END IF
         IF (pasture_pft(cover_type(i,itile))) sum_pasture_fract(i) = sum_pasture_fract(i) + cover_fract(i,itile)
         IF (is_crop(i,itile)) sum_crop_fract(i) = sum_crop_fract(i) + cover_fract(i,itile)
      END DO
    END DO

    !-- establishment of pastures (at the expense of grass land)

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (.NOT. glacier(i) .AND. dynamic_pft(cover_type(i,itile))) THEN
             IF (sum_pasture_fract(i) <= sum_grass_fract(i)) THEN
 
               ! There are enough grasslands

                IF (woody_pft(cover_type(i,itile))) THEN
                   cover_fract_inst(i,itile) = cover_fract_pot(i,itile)
                ELSE
                   cover_fract_inst(i,itile) = cover_fract_pot(i,itile) &
                                      * (sum_grass_fract(i) - sum_pasture_fract(i) + ngrass * fract_small) / sum_grass_fract(i)
                END IF
             ELSE

                ! There are not enough grasslands for all pastures, so that woody types are also reduced

                IF (woody_pft(cover_type(i,itile))) THEN
                   cover_fract_inst(i,itile) = cover_fract_pot(i,itile) &
                                              * (1._dp - sum_pasture_fract(i) - 3 * fract_small) &
                                              / (1._dp - sum_grass_fract(i) - 3 * fract_small)
                ELSE
                   cover_fract_inst(i,itile) = fract_small
                END IF
             END IF
          ELSE IF (.NOT. glacier(i) .AND. pasture_pft(cover_type(i,itile))) THEN
             cover_fract_inst(i,itile) = cover_fract(i,itile)
          ELSE IF (.NOT. glacier(i)) THEN
             cover_fract_inst(i,itile) = fract_small
          ELSE  ! glacier
             cover_fract_inst(i,itile) = cover_fract(i,itile)
          END IF
       END DO
    END DO

    !-- establishment of crops (at the expense of all natural pfts)

    DO itile = 1,ntiles
       WHERE (is_crop(:,itile))
          cover_fract_inst(:,itile) = cover_fract(:,itile)
       ELSEWHERE (.NOT. pasture_pft(cover_type(:,itile)))
          cover_fract_inst(:,itile) = cover_fract_inst(:,itile) * (1._dp - sum_pasture_fract(:) - sum_crop_fract(:)) &
                                     / (1._dp - sum_pasture_fract(:) - fract_small)
       END WHERE
    END DO

    CALL scale_cover_fract (nidx, ntiles, is_present(:,:), is_glacier(:,:), is_naturalveg(:,:), cover_fract_inst(:,:))

    IF (ANY(SUM(cover_fract_inst(:,:),DIM=2) > 1._dp + ntiles*EPSILON(1._dp)) .OR. &
        ANY(SUM(cover_fract_inst(:,:),DIM=2) < 1._dp - ntiles*EPSILON(1._dp))) THEN
       WRITE (message_text,*) 'sum of cover_fract_inst /= 1: ', &
            MINVAL(SUM(cover_fract_inst(:,:),DIM=2)), MAXVAL(SUM(cover_fract_inst(:,:),DIM=2)), &
            MAXLOC(SUM(cover_fract_inst(:,:),DIM=2))
       CALL finish ('calc_cover_fract_inst', message_text)
    END IF
    DO itile = 1,ntiles
       IF (ANY(.NOT. glacier(:) .AND. is_present(:,itile) .AND. cover_fract_inst(:,itile) < fract_small)) THEN
          WRITE (message_text,*) 'cover_fract_inst too small: ', MINVAL(cover_fract_inst(:,itile)), &
               MINLOC(cover_fract_inst(:,itile),MASK=.NOT. glacier(:))
          CALL finish ('calc_cover_fract_inst', message_text)
       END IF
    END DO

  END SUBROUTINE calc_cover_fract_inst

!------------------------------------------------------------------------------

  SUBROUTINE fpc_to_cover_fract_pot (nidx, ntiles, act_fpc, bare_fpc, is_present, is_glacier, is_naturalveg, &
                                     woody_pft, dynamic_pft, cover_type, cover_fract_pot)
!
! !DESCRIPTION:
!
! Convert fractional plant cover (act_fpc) calculated within dynveg to jsbach
! cover fractions reflecting natural vegetation, only (cover_fract_pot).
! The only difference between act_fpc and cover_fract_pot is caused by bare_fpc,
! which is added to the grass cover fractions.
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx          ! vector length
    INTEGER,  INTENT(in)  :: ntiles        ! number of tiles (different land cover type mosaics per gridcell)
    REAL(dp), INTENT(in)  :: act_fpc(:,:)  ! actual fpc in dynveg
    REAL(dp), INTENT(in)  :: bare_fpc(:)   ! dynveg bare soil fraction
    LOGICAL,  INTENT(in)  :: is_present(:,:) ! logical mask for surface handled by jsbach
    LOGICAL,  INTENT(in)  :: is_glacier(:,:) ! logical glacier mask on tiles
    LOGICAL,  INTENT(in)  :: is_naturalveg(:,:) ! logical mask of natural vegetation
    LOGICAL,  INTENT(in)  :: woody_pft(:)  ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: dynamic_pft(:)! PFTs taking part in competition
    INTEGER,  INTENT(in)  :: cover_type(:,:) ! cover types
!
! !OUTPUT PARAMETERS: 
!
    REAL(dp), INTENT(out) :: cover_fract_pot(:,:) ! cover fractions if there was only natural vegetation

! !LOCAL VARIABLES:
!
    INTEGER   :: i, itile
    REAL(dp)  :: sum_grass_fpc(nidx)     ! part of the vegetated area covered by grass
    LOGICAL   :: glacier(nidx)           ! 2 dim. glacier mask
!------------------------------------------------------------------------------

!-- Initialization

    cover_fract_pot(:,:) = 0._dp

!-- define a 2 dimensional glacier mask
    glacier(:) = ANY(is_glacier(:,:),DIM=2)

!-- Find out fpc of grass

    sum_grass_fpc(:) = 0._dp

    DO i = 1,nidx
      DO itile = 1,ntiles
         IF (.NOT. woody_pft(cover_type(i,itile)) .AND. dynamic_pft(cover_type(i,itile))) &
            sum_grass_fpc(i) = sum_grass_fpc(i) + act_fpc(i,itile)
      END DO
    END DO

!-- calculate new cover_fractions

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (.NOT. glacier(i) .AND. dynamic_pft(cover_type(i,itile))) THEN
             IF (woody_pft(cover_type(i,itile))) THEN
                cover_fract_pot(i,itile) = act_fpc(i,itile)
             ELSE
                cover_fract_pot(i,itile) = act_fpc(i,itile) * (1._dp + bare_fpc(i) / MAX(sum_grass_fpc(i),EPSILON(1._dp)))
             END IF
          END IF
       END DO
    END DO

    CALL scale_cover_fract (nidx, ntiles, is_present(:,:), is_glacier(:,:), is_naturalveg(:,:), cover_fract_pot(:,:))

    IF (ANY(SUM(cover_fract_pot(:,:),DIM=2) > 1._dp + ntiles*EPSILON(1._dp)) .OR. &
        ANY(SUM(cover_fract_pot(:,:),DIM=2) < 1._dp - ntiles*EPSILON(1._dp))) THEN
       WRITE (message_text,*) 'sum of cover_fract_pot /= 1: ', &
            MINVAL(SUM(cover_fract_pot(:,:),DIM=2)), MAXVAL(SUM(cover_fract_pot(:,:),DIM=2)), &
            MAXLOC(SUM(cover_fract_pot(:,:),DIM=2))
       CALL finish ('fpc_to_cover_fract_pot', message_text)
    END IF
    DO itile = 1,ntiles
       IF (ANY(.NOT. glacier(:) .AND. is_present(:,itile) .AND. cover_fract_pot(:,itile) < fract_small)) THEN
          WRITE (message_text,*) 'cover_fract_pot too small: ', MINVAL(cover_fract_pot(:,itile)), itile
          CALL finish ('fpc_to_cover_fract_pot', message_text)
       END IF
    END DO

  END SUBROUTINE fpc_to_cover_fract_pot

!------------------------------------------------------------------------------

  SUBROUTINE calc_veg_ratio_max (desert_fpc, rock_fract, is_glacier, veg_ratio_max)
!
! !DESCRIPTION:
!
! Calculate the maximum vegetated fraction of the grid cells from the dynveg 
! desert fraction.
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(in)  :: desert_fpc(:)   ! desert fraction
    REAL(dp), INTENT(in)  :: rock_fract(:)   ! grid cell fraction not vegetated nor desert
    LOGICAL,  INTENT(in)  :: is_glacier(:,:) ! logical glacier mask on tiles
!
! !OUTPUT PARAMETERS: 
!
    REAL(dp), INTENT(out) :: veg_ratio_max(:)

!------------------------------------------------------------------------------

!-- calculate veg_ratio_max

    WHERE (.NOT. ANY(is_glacier(:,:),DIM=2))
       veg_ratio_max(:) = MAX(fract_small, 1._dp - (desert_fpc(:) + rock_fract(:)))
    ELSEWHERE
       veg_ratio_max(:) = 0._dp
    END WHERE

  END SUBROUTINE calc_veg_ratio_max
!------------------------------------------------------------------------------

  SUBROUTINE scale_fpc (nidx, ntiles, glacier, dynamic_pft, cover_type, act_fpc, bare_fpc)
!
! !DESCRIPTION:
!
! Rescaling of fractional plant coverage to assure that
!  - the sum of act_fpc and bare soil is exactly one
!  - all tiles have a minimum vegetated fraction
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx           ! vector length
    INTEGER,  INTENT(in)  :: ntiles         ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: glacier(:)     ! logical glacier mask
    LOGICAL,  INTENT(in)  :: dynamic_pft(:) ! PFTs taking part in competition
    INTEGER,  INTENT(in)  :: cover_type(:,:)! land cover types
!
! !IN- and OUTPUT PARAMETERS:
! 
    REAL(dp), INTENT(inout) :: act_fpc(:,:) ! actual fpc in dynveg
    REAL(dp), INTENT(inout) :: bare_fpc(:)  ! fraction of bare ground

! !LOCAL VARIABLES:
!
    LOGICAL   :: just_return
    INTEGER   :: i, itile       ! indices
    INTEGER   :: nsparce(nidx)  ! number of PFTs with less then fract_small vegetation
    REAL(dp)  :: sum_fpc(nidx)  ! sum of the different cover fractions
    REAL(dp)  :: excess         ! extra fraction, that needs to be redistibuted
    REAL(dp)  :: scalable       ! sum of fraction of tiles that can be scaled
    REAL(dp)  :: rescale        ! factor to rescale fraction of tiles
!------------------------------------------------------------------------------

    ! glacier
    WHERE (glacier(:))
       bare_fpc(:) = 1._dp
    END WHERE
    DO itile = 1,ntiles
       WHERE (glacier(:))
          act_fpc(:,itile) = 0._dp
       END WHERE
    END DO

    ! test, if scaling is necessary
    just_return = .true.
    sum_fpc(:) = bare_fpc(:)
    DO i = 1,nidx
       DO itile = 1,ntiles
          sum_fpc(i) = sum_fpc(i) + act_fpc(i,itile)
          IF (.NOT. glacier(i) .AND. dynamic_pft(cover_type(i,itile)) .AND. &
              act_fpc(i,itile) < fract_small) just_return = .false.
       END DO
       IF (ABS(1._dp - sum_fpc(i)) > 3._dp * EPSILON(1._dp)) just_return = .false.
    END DO
    
    IF (just_return) RETURN
    
    CALL message('update_dynveg','Warning: act_fpc is rescaled!' )

    ! sum act_fpc and bare_fpc
    sum_fpc(:) = bare_fpc(:)
    nsparce(:) = 0
    
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (dynamic_pft(cover_type(i,itile))) THEN
             IF (act_fpc(i,itile) > fract_small) THEN
                sum_fpc(i) = sum_fpc(i) + act_fpc(i,itile)
             ELSE
                nsparce(i) = nsparce(i) + 1
                act_fpc(i,itile) = fract_small
             END IF
          ELSE
             act_fpc(i,itile) = 0._dp
          END IF
       END DO
    END DO

    sum_fpc(:) = sum_fpc(:) + REAL(nsparce(:)) * fract_small

    ! scaling of bare_fpc
    WHERE (sum_fpc(:) <= 1._dp - EPSILON(1._dp))
       bare_fpc(:) = bare_fpc(:) + (1._dp - sum_fpc(:))
    ELSEWHERE (sum_fpc(:) >= 1._dp + EPSILON(1._dp))
       bare_fpc(:) = bare_fpc(:) / sum_fpc(:)
    END WHERE

    ! scaling of act_fpc
    DO i = 1,nidx
       excess = 0._dp
       DO itile = 1,ntiles
          IF (dynamic_pft(cover_type(i,itile))) THEN
             IF (sum_fpc(i) >= 1._dp + EPSILON(1._dp) .AND. act_fpc(i,itile) > fract_small * sum_fpc(i)) THEN
                act_fpc(i,itile) = act_fpc(i,itile) / sum_fpc(i)
             ELSE IF (sum_fpc(i) >= 1._dp + EPSILON(1._dp)) THEN
                excess = excess + fract_small - act_fpc(i,itile) / sum_fpc(i)
                act_fpc(i,itile) = fract_small
             END IF
          END IF
       END DO
       scalable = bare_fpc(i)
       DO itile = 1,ntiles
          IF (dynamic_pft(cover_type(i,itile))) THEN
             IF (act_fpc(i,itile) > 2._dp * fract_small) scalable = scalable + act_fpc(i,itile)
          END IF
       END DO
       rescale = (scalable - excess) / scalable
       bare_fpc(i) = bare_fpc(i) * rescale
       DO itile = 1,ntiles
          IF (dynamic_pft(cover_type(i,itile))) THEN
             IF (act_fpc(i,itile) > 2._dp * fract_small) act_fpc(i,itile) = act_fpc(i,itile) * rescale
          END IF
       END DO
    END DO

  END SUBROUTINE scale_fpc

END MODULE mo_dynveg
