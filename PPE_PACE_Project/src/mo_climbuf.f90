!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Module to calculate the climate variables
! used by mo_dynveg
!
!-----------------------------------------------------------------------------------------------------------------------------------

MODULE mo_climbuf

!-----------------------------------------------------------------------------------------------------------------------------------

  USE mo_jsbach_grid,  ONLY: grid_type, domain_type
  USE mo_kind,         ONLY: dp
  USE mo_jsbach,       ONLY: debug, missing_value
  USE mo_exception,    ONLY: finish, message, message_text
  USE mo_netCDF,       ONLY: FILE_INFO, NF_MAX_NAME
  USE mo_io,           ONLY: IO_open, IO_READ
  USE mo_mpi,          ONLY: p_parallel_io, p_parallel, p_io, p_bcast

  IMPLICIT NONE

  ! --- public subroutines ---

  PUBLIC :: config_climbuf           ! Reads the climbuf namelist
  PUBLIC :: init_climbuf             ! Allocates and initializes memory
  PUBLIC :: update_climbuf           ! Collects weather data for climatology
  PUBLIC :: climbuf_type
  PUBLIC :: climbuf_options_type

  TYPE climbuf_type                  ! contains the state variables
     REAL(dp), POINTER :: min_mmtemp20(:)
     REAL(dp), POINTER :: max_mmtemp20(:)
     REAL(dp), POINTER :: prev_year_gdd(:,:)
     REAL(dp), POINTER :: gdd_upper_tlim(:,:)
     REAL(dp), POINTER :: prev_year_npp(:,:)
     REAL(dp), POINTER :: ave_npp5(:,:)
     REAL(dp), POINTER :: prev_year_precip(:)
     REAL(dp), POINTER :: prev_year_soiltemp(:,:)
     REAL(dp), POINTER :: prev_year_soilmoist(:,:)
     REAL(dp), POINTER :: max_wind10(:)
     REAL(dp), POINTER :: prev_day_max_wind10(:)
     REAL(dp), POINTER :: rel_hum_air(:)          !! smoothend relative humidity
     REAL(dp), POINTER :: rel_hum_air_yDay(:)     !! smoothend relative humidity at new_day (needed as cbalone forcing)
     REAL(dp), POINTER :: prev_day_mean_wind10(:)      ! ToDo: Is this really the windspeed we need?
     REAL(dp), POINTER :: prev_day_mean_vol_moist(:)       ! top layer volumetric_moisture daily average
     REAL(dp), POINTER :: curr_day_temp_max(:)             ! Maximum air temp of current day until current time step
     REAL(dp), POINTER :: curr_day_temp_min(:)             ! Minimum air temp of current day until current time step
     REAL(dp), POINTER :: prev_day_temp_max(:)             ! Maximum air temp of previous day
     REAL(dp), POINTER :: prev_day_temp_min(:)             ! Minimum air temp of previous day
!     REAL(dp), POINTER :: dew_point_temp(:)
     REAL(dp), POINTER :: prev_day_precip_mean(:)
  END TYPE climbuf_type

  TYPE(climbuf_type), PUBLIC :: climbuf

  TYPE climbuf_options_type
     LOGICAL                :: init_running_means !! initialize running means of the climate buffer
     LOGICAL                :: read_climbuf       !! Read long term climate variables from file
     CHARACTER(nf_max_name) :: climbuf_file_name
  END TYPE climbuf_options_type
  TYPE(climbuf_options_type), PUBLIC :: climbuf_options


  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE                             ! Make ALL following objects private

  REAL(dp), PARAMETER :: persist_rel_hum = 0.95_dp   ! factor concerning the smoothing of relative air humidity (~14 days)
  REAL(dp), PARAMETER :: persist_wind10  = 0.9995_dp ! factor concerning the smoothing of daily maximum wind speed
  LOGICAL,  PARAMETER :: test_output     = .FALSE.   ! Output all fields of climbuf_type in stream

  ! --- private fields (state variables and fields for intra-module communication)
  INTEGER, SAVE       :: time_steps_per_day   ! number of time steps per day

  REAL(dp), POINTER :: temp_sum_month(:)
  REAL(dp), POINTER :: temp_sum_day(:)
  REAL(dp), POINTER :: surftemp_sum(:,:)
  REAL(dp), POINTER :: wind_sum_day(:)
  REAL(dp), POINTER :: precip_sum_day(:)
  REAL(dp), POINTER :: rel_soilmoist_sum(:,:)
  REAL(dp), POINTER :: max_wind10_act(:)
  REAL(dp), POINTER :: min_mmtemp_of_yr(:)
  REAL(dp), POINTER :: max_mmtemp_of_yr(:)
  REAL(dp), POINTER :: ann_npp_sum(:,:)
  REAL(dp), POINTER :: ann_precip_sum(:)
  REAL(dp), POINTER :: ann_gdd_sum(:,:)
  REAL(dp), POINTER :: gdd_upper_tlim_sum(:,:)
  REAL(dp), POINTER :: vol_moist_sum_day(:)

CONTAINS

  ! --- config_climbuf ------------------------------------------------------------------
  !
  ! Reads the namelist

  SUBROUTINE config_climbuf(climbuf_options)

    USE mo_namelist,         ONLY: position_nml, POSITIONED
    USE mo_jsbach,           ONLY: nml_unit
    USE mo_io_units,         ONLY: nout

    TYPE(climbuf_options_type), INTENT(out) :: climbuf_options

    !! --- locals

    INTEGER                :: read_status, f_unit

    !! Namelist Parameters
    CHARACTER(NF_MAX_NAME) :: climbuf_file_name  !! file name for initial climate data
    LOGICAL                :: init_running_means !! initialize running means of the climate buffer
    LOGICAL                :: read_climbuf       !! Read long term climate variables from file

    INCLUDE 'climbuf_ctl.inc'


    !! Read namelist climbuf_ctl
    IF (p_parallel_io) THEN

       ! define default values
       init_running_means = .FALSE.
       read_climbuf = .FALSE.
       climbuf_file_name = 'climbuf.nc'

       f_unit = position_nml ('CLIMBUF_CTL', nml_unit, status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (f_unit, climbuf_ctl)
          CALL message('config_climbuf', 'Namelist CLIMBUF_CTL: ')
          WRITE(nout, climbuf_ctl)
       END SELECT
    ENDIF
    IF (p_parallel_io) THEN

       climbuf_options%init_running_means = init_running_means
       IF (init_running_means) THEN
          CALL message('config_climbuf', 'Climate buffer running means are initialised')
       ENDIF

       climbuf_options%read_climbuf = read_climbuf
       climbuf_options%climbuf_file_name = climbuf_file_name
       IF (read_climbuf) THEN
          CALL message('config_climbuf','Climate buffer read from external file')
          CALL message('config_climbuf','File containing initial values for climbuf: '//TRIM(climbuf_file_name))
       ENDIF
    ENDIF
    IF (p_parallel) THEN
       CALL p_bcast(climbuf_options%read_climbuf, p_io)
       CALL p_bcast(climbuf_options%init_running_means, p_io)
    ENDIF

  END SUBROUTINE config_climbuf

  ! --- init_climbuf() -------------------------------------------------------------------------------------------------------------
  !
  ! This routine initializes the climate buffer module. It has to be called before the first time step.

  SUBROUTINE init_climbuf (grid, domain, ntiles, isRestart, fileformat, fileztype, climbuf, stream)

    USE mo_time_control,           ONLY: delta_time !! time step in seconds
    USE mo_memory_base,            ONLY: new_stream,default_stream_setting, &
                                         add =>add_stream_element
    USE mo_linked_list,            ONLY: t_stream, LAND, TILES
    USE mo_netCDF,                 ONLY: max_dim_name, io_inq_varid, io_get_var_double
    USE mo_output,                 ONLY: veg_table
    USE mo_temp,                   ONLY: zreal2d, zreal2d_ptr, zzreal2d, zzreal3d
    USE mo_tr_scatter,             ONLY: scatter_field


    ! INTENT(in)
    !-----------
    TYPE(grid_type),         INTENT(in)           :: grid
    TYPE(domain_type),       INTENT(in)           :: domain
    INTEGER,                 INTENT(in)           :: ntiles
    LOGICAL,                 INTENT(in)           :: isRestart
    INTEGER,                 INTENT(in)           :: fileformat      ! output file format
    INTEGER,                 INTENT(in)           :: fileztype       ! output file compression
    TYPE(t_stream), POINTER,                OPTIONAL :: stream

    ! INTENT(out)
    !------------
    TYPE(climbuf_type),         INTENT(out)       :: climbuf


    ! local variables
    !----------------
    TYPE(t_stream),POINTER      :: climbufStream

    INTEGER :: g_nland, l_nland
    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim3p(2), dim3(2)
    CHARACTER(LEN=max_dim_name) :: dim3n(2), dim1n(1)

    TYPE(FILE_INFO)       :: climbuf_file    !! Input file for climate buffer
    INTEGER :: IO_file_id, IO_var_id, i
    LOGICAL :: lAvail

    IF (debug) CALL message('init_climbuf','Start initialization of climate buffer')

    g_nland        = grid%nland
    l_nland        = domain%nland

    ! --- read the climbuf namelist

    CALL config_climbuf(climbuf_options)

    ! --- compute the number of time steps per day

    IF (MOD(86400,INT(delta_time)) .NE. 0) THEN ! For computing the mean day temperature (see below) it is assumed that each day ..
                                                ! is computed with a fixed number of time steps. Therefore the program is  ..
                                                ! stopped when day length is not an integer multiple of the time step.
       CALL finish("init_climbuf","ERROR: Day length is not an integer multiple of the time step!")
    ELSE
       time_steps_per_day = 86400/INT(delta_time)
    END IF

    ! --- Open new stream for climate data

    IF (PRESENT(stream)) THEN
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'climbuf', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
        climbufStream => stream
    ELSE
       ! Add new stream
       CALL new_stream(climbufStream, 'climbuf', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(climbufStream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    ! --- Define state variables as stream elements

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim3p = (/ l_nland, ntiles /)
    dim3  = (/ g_nland, ntiles /)
    dim3n(1) = 'landpoint'
    dim3n(2) = 'tiles'

    ! code numbers of this routine are 1-24, 27, 29, 30, 32, 33, 37, 42, 180 and 198 using the GRIB veg_table
    CALL add(climbufStream,'temp_sum_day',         temp_sum_day, dim1p, dim1, dimnames=dim1n, laccu=.false., code=1, &
             lmiss=.TRUE., missval=missing_value,  reset=0._dp, lpost=test_output)
    CALL add(climbufStream,'vol_moist_sum_day',      vol_moist_sum_day, dim1p, dim1, dimnames=dim1n, laccu=.false., &
             lmiss=.TRUE., missval=missing_value,  reset=0._dp, lpost=test_output, contnorest=.true.)
    CALL add(climbufStream,'temp_sum_month',     temp_sum_month, dim1p, dim1, dimnames=dim1n, laccu=.false., code=2, &
             lmiss=.TRUE., missval=missing_value,  reset=0._dp, lpost=test_output)
    CALL add(climbufStream,'min_mmtemp_of_yr', min_mmtemp_of_yr, dim1p, dim1, dimnames=dim1n, laccu=.false., code=3, &
             lmiss=.TRUE., missval=missing_value,  lpost=test_output)
    CALL add(climbufStream,'max_mmtemp_of_yr', max_mmtemp_of_yr, dim1p, dim1, dimnames=dim1n, laccu=.false., code=4, &
             lmiss=.TRUE., missval=missing_value,  lpost=test_output)
    CALL add(climbufStream,'min_mmtemp20', climbuf%min_mmtemp20, dim1p, dim1, dimnames=dim1n, laccu=.false., code=5, &
             lmiss=.TRUE., missval=missing_value,  reset=0._dp,  lpost=test_output)
    CALL add(climbufStream,'max_mmtemp20', climbuf%max_mmtemp20, dim1p, dim1, dimnames=dim1n, laccu=.false., code=6, &
             lmiss=.TRUE., missval=missing_value,  reset=0._dp,  lpost=test_output)
    CALL add(climbufStream,'surftemp_sum',         surftemp_sum, dim3p, dim3, dimnames=dim3n, laccu=.false., code=7, &
             lmiss=.TRUE., missval=missing_value,  reset=0._dp, lpost=test_output)
    CALL add(climbufStream,'prev_year_soiltemp',  &
                                     climbuf%prev_year_soiltemp, dim3p, dim3, dimnames=dim3n, laccu=.false., code=8, &
             lmiss=.TRUE., missval=missing_value, lpost=test_output)
    CALL add(climbufStream,'rel_soilmoist_sum',rel_soilmoist_sum,dim3p, dim3, dimnames=dim3n, laccu=.false., code=9, &
             lmiss=.TRUE., missval=missing_value, reset=0._dp,  lpost=test_output)
    CALL add(climbufStream,'prev_year_soilmoist', &
                                    climbuf%prev_year_soilmoist, dim3p, dim3, dimnames=dim3n, laccu=.false., code=10, &
             lmiss=.TRUE., missval=missing_value, lpost=test_output)
    CALL add(climbufStream,'ann_precip_sum',     ann_precip_sum, dim1p, dim1, dimnames=dim1n, laccu=.false., code=11, &
             lmiss=.TRUE., missval=missing_value, reset=0._dp, lpost=test_output)
    CALL add(climbufStream,'prev_year_precip',    &
                                       climbuf%prev_year_precip, dim1p, dim1, dimnames=dim1n, laccu=.false., code=12, &
             lmiss=.TRUE., missval=missing_value, lpost=test_output)
    CALL add(climbufStream,'ann_gdd_sum',           ann_gdd_sum, dim3p, dim3, dimnames=dim3n, laccu=.false., code=13, &
             lmiss=.TRUE., missval=missing_value, reset=0.0_dp, lpost=test_output)
    CALL add(climbufStream,'prev_year_gdd',climbuf%prev_year_gdd,dim3p, dim3, dimnames=dim3n, laccu=.false., code=14, &
             lmiss=.TRUE., missval=missing_value, lpost=test_output)
    CALL add(climbufStream,'gdd_upper_tlim_sum',gdd_upper_tlim_sum,dim3p,dim3,dimnames=dim3n, laccu=.false., code=15, &
             lmiss=.TRUE., missval=missing_value, reset=0._dp, lpost=test_output)
    CALL add(climbufStream,'gdd_upper_tlim',climbuf%gdd_upper_tlim,dim3p,dim3,dimnames=dim3n, laccu=.false., code=16, &
             lmiss=.TRUE., missval=missing_value, lpost=test_output)
    CALL add(climbufStream,'ann_npp_sum',           ann_npp_sum, dim3p, dim3, dimnames=dim3n, laccu=.false., code=17, &
             lmiss=.TRUE., missval=missing_value, reset=0._dp, lpost=test_output)
    CALL add(climbufStream,'prev_year_npp',climbuf%prev_year_npp,dim3p, dim3, dimnames=dim3n, laccu=.false., code=18, &
             longname='NPP of the previous year', units='mol(CO2) m-2(canopy)',                                       &
             lmiss=.TRUE., missval=missing_value, lpost=.true.)
    CALL add(climbufStream,'ave_npp5',         climbuf%ave_npp5, dim3p, dim3, dimnames=dim3n, laccu=.false., code=19, &
             longname='five year running mean of annual NPP', units='mol(CO2) m-2(canopy)',                           &
             lmiss=.TRUE., missval=missing_value, lpost=.true.)
    CALL add(climbufStream,'max_wind10_act',     max_wind10_act, dim1p, dim1, dimnames=dim1n, laccu=.false., code=22, &
             lmiss=.TRUE., missval=missing_value, reset=0._dp, lpost=test_output)
    CALL add(climbufStream,'max_wind10',     climbuf%max_wind10, dim1p, dim1, dimnames=dim1n, laccu=.false., code=23, &
             longname='typical maximum daily 10m wind speed of the previous years', units='m/s',                      &
             lmiss=.TRUE., missval=missing_value, lpost=.true.)
    CALL add(climbufStream,'prev_day_max_wind10',                                                                     &
                                    climbuf%prev_day_max_wind10, dim1p, dim1, dimnames=dim1n, laccu=.false., code=20, &
             longname='maximum 10m wind speed of the previous day', units='m/s',                                      &
             lmiss=.TRUE., missval=missing_value, lpost=.true.)
    CALL add(climbufStream,'rel_hum_air',   climbuf%rel_hum_air, dim1p, dim1, dimnames=dim1n, laccu=.false., code=21, &
             lmiss=.TRUE., missval=missing_value, lpost=test_output)
    CALL add(climbufStream,'rel_hum_air_yDay',                                                                        &
                                       climbuf%rel_hum_air_yDay, dim1p, dim1, dimnames=dim1n, laccu=.false., code=24, &
             longname='smoothed (~14 day) relative humidity of the air', units='%',                                   &
             lmiss=.TRUE., missval=missing_value, lpost=.true., contnorest=.TRUE.)
    CALL add(climbufStream,'prev_day_mean_wind10', &
                                    climbuf%prev_day_mean_wind10,dim1p, dim1, dimnames=dim1n, laccu=.false., code=42, &
             lmiss=.TRUE., missval=missing_value,lpost=.true., contnorest=.true.)
    CALL add(climbufStream,'wind_sum_day',         wind_sum_day, dim1p, dim1, dimnames=dim1n, laccu=.false., code=27, &
             lmiss=.TRUE., missval=missing_value,  reset=0._dp, lpost=.false., contnorest=.true.)
    CALL add(climbufStream,'prev_day_temp_max',climbuf%prev_day_temp_max,dim1p, dim1, dimnames=dim1n,        code=29, &
             laccu=.false.,lmiss=.TRUE., missval=missing_value, lpost=.true., contnorest=.true.)
    CALL add(climbufStream,'prev_day_temp_min',climbuf%prev_day_temp_min,dim1p,dim1,dimnames=dim1n,laccu=.false.,code=30, &
             lmiss=.TRUE., missval=missing_value, lpost=.true., contnorest=.true.)
    CALL add(climbufStream,'curr_day_temp_max',climbuf%curr_day_temp_max,dim1p,dim1,dimnames=dim1n,laccu=.false.,code=32, &
             lmiss=.TRUE., missval=missing_value, lpost=.true., contnorest=.true.)
    CALL add(climbufStream,'curr_day_temp_min',climbuf%curr_day_temp_min,dim1p,dim1,dimnames=dim1n,laccu=.false.,code=33, &
             lmiss=.TRUE., missval=missing_value, lpost=.true., contnorest=.true.)
    CALL add(climbufStream,'prev_day_precip_mean',climbuf%prev_day_precip_mean,dim1p,dim1,dimnames=dim1n,        code=37, &
             laccu=.false.,lpost=.true.,lmiss=.TRUE., missval=missing_value, reset=0._dp, contnorest=.true.)
    CALL add(climbufStream,'prev_day_mean_vol_moist',climbuf%prev_day_mean_vol_moist, dim1p,dim1,dimnames=dim1n,code=180, &
             laccu=.false., lmiss=.TRUE., missval=missing_value, lpost=.true., contnorest=.true.)
    CALL add(climbufStream,'precip_sum_day',precip_sum_day,dim1p,dim1,dimnames=dim1n,                           code=198, &
             laccu=.false.,lpost=.false.,lmiss=.TRUE.,missval=missing_value, contnorest=.true.)

    ! -- Initialisation needed if neither restart nor climate buffer file available, or not all of the variables are 
    !         found in the restart file and option contnorest is true. 
    !         The arrays will be overwritten from restart or climate buffer file respectively.
    !          
    climbuf%prev_year_soiltemp      =   0._dp
    climbuf%prev_year_soilmoist     =   0._dp
    climbuf%prev_year_precip        =   0._dp
    climbuf%prev_year_gdd           =   0._dp
    climbuf%gdd_upper_tlim          =   0._dp
    climbuf%prev_year_npp           =   0._dp
    climbuf%ave_npp5                =   0._dp
    climbuf%prev_day_max_wind10     =   0._dp
    climbuf%prev_day_precip_mean    =   0._dp
    climbuf%max_wind10              =  15._dp
    climbuf%rel_hum_air             =  50._dp
    climbuf%rel_hum_air_yDay        =  50._dp
    climbuf%prev_day_mean_wind10    =   0._dp
    climbuf%prev_day_temp_min       = 400._dp
    climbuf%prev_day_temp_max       = 150._dp
    climbuf%curr_day_temp_min       = 400._dp
    climbuf%curr_day_temp_max       = 150._dp
    climbuf%prev_day_mean_vol_moist =  0.5_dp
    climbuf%min_mmtemp20            =   0._dp
    climbuf%max_mmtemp20            =   0._dp

    ! -- Initialisation at beginning of an experiment - will be overwritten in restarted runs
    min_mmtemp_of_yr     =  1000._dp
    max_mmtemp_of_yr     = -1000._dp
    temp_sum_day         =     0._dp
    temp_sum_month       =     0._dp
    surftemp_sum         =     0._dp
    rel_soilmoist_sum    =     0._dp
    ann_precip_sum       =     0._dp
    ann_gdd_sum          =     0._dp
    gdd_upper_tlim_sum   =     0._dp
    ann_npp_sum          =     0._dp
    max_wind10_act       =     0._dp

    ! If this is a restart run we can exit now since the model state variables are read from restart file
    IF (isRestart) THEN
       IF (debug) THEN
          CALL message('init_climbuf()','State variables read from restart file')
       END IF
       RETURN
    ENDIF

    ! --- initializations -------------------------------
    IF (climbuf_options%read_climbuf) THEN

!     Initialize reading
      ALLOCATE(zreal2d(domain%ndim,domain%nblocks), &
               zzreal3d(grid%nlon,grid%nlat,ntiles), &
               zzreal2d(grid%nlon,grid%nlat))
      NULLIFY(zreal2d_ptr)
      IF (p_parallel_io) THEN
         climbuf_file%opened = .FALSE.
         CALL IO_open(TRIM(climbuf_options%climbuf_file_name), climbuf_file, IO_READ)
         IO_file_id = climbuf_file%file_id
         IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
      END IF

!     min_mmtemp_of_yr
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'min_mmtemp_of_yr', IO_var_id,lAvail)
         IF (.NOT. lAvail) THEN
           CALL message('InitClimbuf','min_mmtemp_of_yr not found in climbuf file - please manually set "read_climbuf"')
           CALL message('InitClimbuf','to ".false." in the namelist definitions of your run-script or re-create climbuf')
           CALL message('InitClimbuf','file (which is typically the cpool file) and re-run the model.')
           CALL finish ('InitClimBuf','Necessary data not found. Terminating!')
         ENDIF
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      min_mmtemp_of_yr(:) = PACK(zreal2d, MASK=domain%mask)

!     max_mmtemp_of_yr
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'max_mmtemp_of_yr', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      max_mmtemp_of_yr(:) = PACK(zreal2d, MASK=domain%mask)

!     min_mmtemp20
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'min_mmtemp20', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%min_mmtemp20(:) = PACK(zreal2d, MASK=domain%mask)

!     max_mmtemp20
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'max_mmtemp20', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%max_mmtemp20(:) = PACK(zreal2d, MASK=domain%mask)

!     prev_year_precip 
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_year_precip', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%prev_year_precip(:) = PACK(zreal2d, MASK=domain%mask)

!     rel_hum_air 
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'rel_hum_air', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%rel_hum_air(:) = PACK(zreal2d, MASK=domain%mask)

!     max_wind10 
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'max_wind10', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%max_wind10(:) = PACK(zreal2d, MASK=domain%mask)

!     prev_day_mean_wind10 
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_day_mean_wind10', IO_var_id,lAvail)
         IF (lAvail) THEN
           CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
         ELSE
           write(message_text,*) 'Warning: prev_day_mean_wind10 not found in ', &
                                  TRIM(climbuf_options%climbuf_file_name),'. Starting with 0.'
           CALL message('init_climbuf',message_text)
           zzreal2d(:,:) = 0._dp
         ENDIF
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      min_mmtemp_of_yr(:) = PACK(zreal2d, MASK=domain%mask)

!     prev_day_max_wind10 
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_day_max_wind10', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%prev_day_max_wind10(:) = PACK(zreal2d, MASK=domain%mask)

!     ave_npp5 
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'ave_npp5', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal3d)
      END IF
      DO i=1,ntiles
         IF (p_parallel_io) zreal2d_ptr => zzreal3d(:,:,i)
         CALL scatter_field(zreal2d_ptr, zreal2d)
         climbuf%ave_npp5(:,i) = PACK(zreal2d, MASK=domain%mask)
      ENDDO

!     prev_year_npp 
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_year_npp', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal3d)
      END IF
      DO i=1,ntiles
         IF (p_parallel_io) zreal2d_ptr => zzreal3d(:,:,i)
         CALL scatter_field(zreal2d_ptr, zreal2d)
         climbuf%prev_year_npp(:,i) = PACK(zreal2d, MASK=domain%mask)
      ENDDO

!     prev_year_gdd 
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_year_gdd', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal3d)
      END IF
      DO i=1,ntiles
         IF (p_parallel_io) zreal2d_ptr => zzreal3d(:,:,i)
         CALL scatter_field(zreal2d_ptr, zreal2d)
         climbuf%prev_year_gdd(:,i) = PACK(zreal2d, MASK=domain%mask)
      ENDDO

!     gdd_upper_tlim 
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'gdd_upper_tlim', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal3d)
      END IF
      DO i=1,ntiles
         IF (p_parallel_io) zreal2d_ptr => zzreal3d(:,:,i)
         CALL scatter_field(zreal2d_ptr, zreal2d)
         climbuf%gdd_upper_tlim(:,i) = PACK(zreal2d, MASK=domain%mask)
      ENDDO

!     prev_year_soiltemp - Does this also make sense for ntsoil > 1 ?
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_year_soiltemp', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal3d)
      END IF
      DO i=1,ntiles
         IF (p_parallel_io) zreal2d_ptr => zzreal3d(:,:,i)
         CALL scatter_field(zreal2d_ptr, zreal2d)
         climbuf%prev_year_soiltemp(:,i) = PACK(zreal2d, MASK=domain%mask)
      ENDDO

!     prev_year_soilmoist - Does this also make sense for nsoil > 1 ?
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_year_soilmoist', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal3d)
      END IF
      DO i=1,ntiles
         IF (p_parallel_io) zreal2d_ptr => zzreal3d(:,:,i)
         CALL scatter_field(zreal2d_ptr, zreal2d)
         climbuf%prev_year_soilmoist(:,i) = PACK(zreal2d, MASK=domain%mask)
      ENDDO

!     prev_day_precip_mean
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_day_precip_mean', IO_var_id, lstat=lAvail)
         IF (.NOT. lAvail) CALL IO_inq_varid(IO_file_id, 'precip', IO_var_id) ! For compatibility with older files
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%prev_day_precip_mean(:) = PACK(zreal2d, MASK=domain%mask)

!     prev_day_temp_max
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_day_temp_max', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%prev_day_temp_max(:) = PACK(zreal2d, MASK=domain%mask)

!     prev_day_temp_min
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'prev_day_temp_min', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%prev_day_temp_min(:) = PACK(zreal2d, MASK=domain%mask)

!     curr_day_temp_max
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'curr_day_temp_max', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%curr_day_temp_max(:) = PACK(zreal2d, MASK=domain%mask)

!     curr_day_temp_min
      IF (p_parallel_io) THEN
         CALL IO_inq_varid(IO_file_id, 'curr_day_temp_min', IO_var_id)
         CALL IO_get_var_double(IO_file_id, IO_var_id, zzreal2d)
      END IF
      CALL scatter_field(zreal2d_ptr, zreal2d)
      climbuf%curr_day_temp_min(:) = PACK(zreal2d, MASK=domain%mask)

      DEALLOCATE(zzreal3d, zreal2d, zzreal2d)

      IF (debug) THEN
         CALL message('init_climbuf()','Read climbuf_fields: Initialization of ClimateBuffer finished.')
      END IF

    ENDIF

    CALL message('init_climbuf','climate buffer initialised')

  END SUBROUTINE init_climbuf

!-----------------------------------------------------------------------------------------------------------------------------------
! Update the climate buffer
!
  SUBROUTINE update_climbuf(kidx, kidx0, kidx1, ntiles, nsoil, gdd_base, upper_tlim, &
       temp_air, temp_surf, wind10, precip, rel_soil_water, &
       npp_rate, rel_hum_air, init_running_means, climbuf, volumetric_moisture)

!-----------------------------------------------------------------------------------------------------------------------------------
   USE mo_exception,              ONLY: message
   USE mo_time_control,           ONLY: get_year_day, get_date_components, previous_date, &
                                        lstart, get_time_step, delta_time
   USE mo_jsbach,                 ONLY: new_day, new_month, new_year
   USE mo_jsbach_constants,       ONLY: tmelt

! INTENT(IN)
!------------
   INTEGER,                    INTENT(in)  ::  kidx, kidx0, kidx1
   INTEGER,                    INTENT(in)  ::  ntiles                     ! number of tiles
   INTEGER,                    INTENT(in)  ::  nsoil                      ! number of soil layers for the hydrology)
   REAL(dp), DIMENSION(:),     INTENT(in)  ::  gdd_base                   ! base temperature to calculate growing degree days
   REAL(dp), DIMENSION(:),     INTENT(in)  ::  upper_tlim                 ! base temperature to calculate gdd_upper_tlim
   REAL(dp), DIMENSION(:),     INTENT(in)  ::  temp_air                   ! air temperature of the lowest level [degC]
   REAL(dp), DIMENSION(:,:),   INTENT(in)  ::  temp_surf                  ! surface temperature [degC]
   REAL(dp), DIMENSION(:),     INTENT(in)  ::  wind10                     ! 10m wind speed [m/s]
   REAL(dp), DIMENSION(:),     INTENT(in)  ::  precip                     ! precipitation [m/s]
   REAL(dp), DIMENSION(:,:),   INTENT(in)  ::  rel_soil_water             ! relative soil water content
   REAL(dp), DIMENSION(:,:),   INTENT(in)  ::  npp_rate                   ! NPP rate [mol(C)/m^2(canopy) s]
   REAL(dp), DIMENSION(:),     INTENT(in)  ::  rel_hum_air                ! relative air humidity
   REAL(dp), DIMENSION(:),     INTENT(in)  ::  volumetric_moisture        ! volumetric soil moisture of the top soil layer

! INTENT(OUT)
!------------
   LOGICAL,                    INTENT(out) ::  init_running_means         ! initialise long term climate variables
   TYPE(climbuf_type),       INTENT(inout) ::  climbuf                    ! buffer for climate variables

! LOCAL VARIABLES
!----------------
   INTEGER :: day_of_month_at_prev_ts
   INTEGER :: istep, time_steps_in_last_year
   LOGICAL :: first_day
   REAL(dp):: month_mean_temp(kidx), prev_day_mean_temp(kidx)
   REAL(dp):: prev_day_mean_vol_moist(kidx)
   REAL(dp):: average_filling(kidx, ntiles)
   INTEGER :: pft
   REAL(dp):: persist

!-----------------------------------------------------------------------------------------------------------------------------------
   IF (debug) CALL message ('update climbuf','starting routine')

   ! --- New day, new month, new year ?

   istep=get_time_step()
   first_day = ( istep < time_steps_per_day)
   time_steps_in_last_year = INT(get_year_day(previous_date)) * time_steps_per_day

   IF (debug .AND. new_day)   CALL message('update climbuf','First time step of new day')
   IF (debug .AND. new_month) CALL message('update climbuf','First time step of new month')
   IF (debug .AND. new_year)  CALL message('update climbuf','First time step of new year')

   !   -------------------
   ! -- Climate Variables
   !   -------------------

   ! --- Variables based on air temperature
   !
   IF (.NOT. new_day .OR. lstart) THEN
      temp_sum_day(kidx0:kidx1) = temp_sum_day(kidx0:kidx1) + temp_air(:)
      vol_moist_sum_day(kidx0:kidx1) = vol_moist_sum_day(kidx0:kidx1) + volumetric_moisture(:)
      WHERE (temp_air(:)+tmelt > climbuf%curr_day_temp_max(kidx0:kidx1))
        climbuf%curr_day_temp_max(kidx0:kidx1) = temp_air(:)+tmelt
      ENDWHERE
      WHERE (temp_air(:)+tmelt < climbuf%curr_day_temp_min(kidx0:kidx1))
        climbuf%curr_day_temp_min(kidx0:kidx1) = temp_air(:)+tmelt
      ENDWHERE
      IF (lstart .AND. climbuf%prev_day_temp_max(kidx0)<climbuf%prev_day_temp_min(kidx0)) THEN
         climbuf%prev_day_temp_max(kidx0:kidx1) = temp_air(:)+tmelt
         climbuf%prev_day_temp_min(kidx0:kidx1) = temp_air(:)+tmelt
      ENDIF
   ELSE
      IF (first_day) THEN
         prev_day_mean_temp(:) = temp_sum_day(kidx0:kidx1) / istep
         prev_day_mean_vol_moist(:)=vol_moist_sum_day(kidx0:kidx1)/istep
      ELSE
         prev_day_mean_temp(:) = temp_sum_day(kidx0:kidx1) / time_steps_per_day
         prev_day_mean_vol_moist(:)=vol_moist_sum_day(kidx0:kidx1) / time_steps_per_day
         climbuf%prev_day_temp_max(kidx0:kidx1) = climbuf%curr_day_temp_max(kidx0:kidx1)
         climbuf%prev_day_temp_min(kidx0:kidx1) = climbuf%curr_day_temp_min(kidx0:kidx1)
      ENDIF
      climbuf%curr_day_temp_max(kidx0:kidx1) = temp_air(:)+tmelt
      climbuf%curr_day_temp_min(kidx0:kidx1) = temp_air(:)+tmelt
      temp_sum_day(kidx0:kidx1) = temp_air(:)
      vol_moist_sum_day(kidx0:kidx1) = volumetric_moisture(:)
      climbuf%prev_day_mean_vol_moist(kidx0:kidx1) = prev_day_mean_vol_moist(:)

   ! -- Growing Degree Days

      DO pft = 1,ntiles
         ann_gdd_sum(kidx0:kidx1,pft) = &
              ann_gdd_sum(kidx0:kidx1,pft) + MAX(prev_day_mean_temp(:) - gdd_base(pft), 0._dp)
         gdd_upper_tlim_sum(kidx0:kidx1,pft) = &
              gdd_upper_tlim_sum(kidx0:kidx1,pft) + MAX(prev_day_mean_temp(:) - upper_tlim(pft), 0._dp)
      END DO

      IF (new_year) THEN

      ! --- provide previous year GDD and reset for current year
         DO pft = 1, ntiles
            climbuf%prev_year_gdd(kidx0:kidx1,pft) = ann_gdd_sum(kidx0:kidx1,pft)
            ann_gdd_sum(kidx0:kidx1,pft) =  0._dp
            climbuf%gdd_upper_tlim(kidx0:kidx1,pft) = gdd_upper_tlim_sum(kidx0:kidx1,pft)
            gdd_upper_tlim_sum(kidx0:kidx1,pft) = 0._dp
         END DO
      ENDIF

   !
   ! --- Minimum/Maximum monthly mean temperature (20 year climatology)
   !
      temp_sum_month(kidx0:kidx1) = temp_sum_month(kidx0:kidx1) + prev_day_mean_temp(:)

      IF (new_month) THEN

      ! --- provide previous month mean temperature and reset for current month
         CALL get_date_components(previous_date, DAY=day_of_month_at_prev_ts)
         month_mean_temp(:) = temp_sum_month(kidx0:kidx1) / day_of_month_at_prev_ts
         temp_sum_month(kidx0:kidx1) = 0._dp

      ! --- calculate coldest/warmest month of the year so far
         min_mmtemp_of_yr(kidx0:kidx1) = min(min_mmtemp_of_yr(kidx0:kidx1), month_mean_temp(:))
         max_mmtemp_of_yr(kidx0:kidx1) = max(max_mmtemp_of_yr(kidx0:kidx1), month_mean_temp(:))

         IF (new_year) THEN

         ! --- bulid kind of 20yr running mean of coldest and warmest monthly temperature
            IF (climbuf_options%init_running_means) THEN
               IF (debug) CALL message ('update climbuf','initialisation of mmtemp20')
               climbuf%min_mmtemp20(kidx0:kidx1) = min_mmtemp_of_yr(kidx0:kidx1)
               climbuf%max_mmtemp20(kidx0:kidx1) = max_mmtemp_of_yr(kidx0:kidx1)
            ELSE
               climbuf%min_mmtemp20(kidx0:kidx1) = (climbuf%min_mmtemp20(kidx0:kidx1)*19._dp + min_mmtemp_of_yr(kidx0:kidx1)) &
                                                    /20._dp
               climbuf%max_mmtemp20(kidx0:kidx1) = (climbuf%max_mmtemp20(kidx0:kidx1)*19._dp + max_mmtemp_of_yr(kidx0:kidx1)) &
                                                    /20._dp
            ENDIF

         ! --- reset for current year
            min_mmtemp_of_yr(kidx0:kidx1) =  1000._dp
            max_mmtemp_of_yr(kidx0:kidx1) = -1000._dp

         END IF

      END IF
   END IF

   IF (lstart .AND. climbuf_options%read_climbuf .AND. new_year) THEN
      climbuf%min_mmtemp20(kidx0:kidx1) = (climbuf%min_mmtemp20(kidx0:kidx1)*19._dp + min_mmtemp_of_yr(kidx0:kidx1))/20._dp
      climbuf%max_mmtemp20(kidx0:kidx1) = (climbuf%max_mmtemp20(kidx0:kidx1)*19._dp + max_mmtemp_of_yr(kidx0:kidx1))/20._dp
   ! --- reset for current year
      min_mmtemp_of_yr(kidx0:kidx1) =  1000._dp
      max_mmtemp_of_yr(kidx0:kidx1) = -1000._dp
   ENDIF

   ! --- Maximum daily wind speed and relative air humidity smoothed in time

   persist = persist_rel_hum ** (1._dp/REAL(time_steps_per_day,dp))
   climbuf%rel_hum_air(kidx0:kidx1) = climbuf%rel_hum_air(kidx0:kidx1) * persist + &
                                      MIN(rel_hum_air(1:kidx),100._dp) * (1._dp - persist)

   IF (.NOT. new_day .OR. lstart) THEN
      WHERE (wind10(:) > max_wind10_act(kidx0:kidx1))
         max_wind10_act(kidx0:kidx1) = wind10(:)
      END WHERE
   ELSE
      climbuf%prev_day_max_wind10(kidx0:kidx1) = max_wind10_act(kidx0:kidx1)
      climbuf%max_wind10(kidx0:kidx1) = climbuf%max_wind10(kidx0:kidx1) * persist_wind10 + &
                                        max_wind10_act(kidx0:kidx1) * (1._dp - persist_wind10)
      max_wind10_act(kidx0:kidx1) = wind10(:)
      climbuf%rel_hum_air_yDay(kidx0:kidx1) = climbuf%rel_hum_air(kidx0:kidx1)
   END IF

   ! --- Mean wind speed of previous day
   IF (.NOT. new_day .OR. lstart) THEN
      wind_sum_day(kidx0:kidx1) = wind_sum_day(kidx0:kidx1) + wind10(:)
   ELSE
      IF (first_day) THEN
         climbuf%prev_day_mean_wind10(kidx0:kidx1) = wind_sum_day(kidx0:kidx1) / istep
      ELSE
         climbuf%prev_day_mean_wind10(kidx0:kidx1) = wind_sum_day(kidx0:kidx1) / time_steps_per_day
      ENDIF
      wind_sum_day(kidx0:kidx1) = wind10(:) 
   ENDIF

   ! Average soil moisture integrated column
   average_filling(:,:) = rel_soil_water(:,:)

   IF (.NOT. new_year .OR. lstart) THEN
      ann_npp_sum(kidx0:kidx1,:)  = ann_npp_sum(kidx0:kidx1,:) + NPP_Rate(1:kidx,:) * delta_time
      ann_precip_sum(kidx0:kidx1) = ann_precip_sum(kidx0:kidx1) + precip(1:kidx) * delta_time
      surftemp_sum(kidx0:kidx1,:)  = surftemp_sum(kidx0:kidx1,:) + temp_surf(:,:)
      rel_soilmoist_sum(kidx0:kidx1,:) = rel_soilmoist_sum(kidx0:kidx1,:) + average_filling(:,:)
   ELSE

      ! --- provide previous year NPP and reset for current year

      climbuf%prev_year_npp(kidx0:kidx1,:) = ann_npp_sum(kidx0:kidx1,:)
      ann_npp_sum(kidx0:kidx1,:) = NPP_Rate(1:kidx,:) * delta_time

      ! --- bulid kind of 5yr running mean of NPP
      IF (climbuf_options%init_running_means) THEN
         climbuf%ave_npp5(kidx0:kidx1,:) = climbuf%prev_year_npp(kidx0:kidx1,:)
      ELSE
         climbuf%ave_npp5(kidx0:kidx1,:) = (climbuf%ave_npp5(kidx0:kidx1,:)*4._dp + climbuf%prev_year_npp(kidx0:kidx1,:)) / 5._dp
      ENDIF

      ! --- provide previous year precipitation and reset for current year

      climbuf%prev_year_precip(kidx0:kidx1) = ann_precip_sum(kidx0:kidx1)
      ann_precip_sum(kidx0:kidx1) = precip(:) * delta_time

      ! --- provide previous year soil temperature and reset for current year

      climbuf%prev_year_soiltemp(kidx0:kidx1,:) = surftemp_sum(kidx0:kidx1,:) / time_steps_in_last_year
      surftemp_sum(kidx0:kidx1,:) = temp_surf(:,:)

      ! --- provide previous year soil moisture and reset for current year

      climbuf%prev_year_soilmoist(kidx0:kidx1,:) = rel_soilmoist_sum(kidx0:kidx1,:) / time_steps_in_last_year
      rel_soilmoist_sum(kidx0:kidx1,:) = average_filling(:,:)

   END IF

   ! Precipitation
   IF (.NOT. new_day .OR. lstart) THEN
     precip_sum_day(kidx0:kidx1) = precip_sum_day(kidx0:kidx1) + precip(:)
   ELSE
     IF (first_day) THEN
       climbuf%prev_day_precip_mean(kidx0:kidx1) = precip_sum_day(kidx0:kidx1) / istep
     ELSE
       climbuf%prev_day_precip_mean(kidx0:kidx1) = precip_sum_day(kidx0:kidx1) / time_steps_per_day
     ENDIF
     precip_sum_day(kidx0:kidx1) = precip(:)
   ENDIF

   init_running_means = climbuf_options%init_running_means
   IF (debug) CALL message ('update climbuf','end of routine')

  END SUBROUTINE update_climbuf

END MODULE mo_climbuf
