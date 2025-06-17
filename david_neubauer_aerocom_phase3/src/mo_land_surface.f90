!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_land_surface

  USE mo_jsbach,        ONLY: debug
  USE mo_jsbach_grid,   ONLY: grid_type, domain_type, kstart, kend, nidx
  USE mo_jsbach_lctlib, ONLY: lctlib_type

  USE mo_netcdf,        ONLY: FILE_INFO
  USE mo_linked_list,   ONLY: t_stream
  USE mo_mpi,           ONLY: p_parallel, p_parallel_io, p_bcast, p_io
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, message_text, finish
  USE mo_util_string,   ONLY: int2string

  IMPLICIT NONE

  PUBLIC ::  init_land_surface
  PUBLIC ::  init_albedo
  PUBLIC ::  init_albedo_land
  PUBLIC ::  config_albedo
  PUBLIC ::  update_land_surface_fast
  PUBLIC ::  update_albedo
  PUBLIC ::  update_albedo_diag
  PUBLIC ::  land_surface_diagnostics
  PUBLIC ::  scale_cover_fract
  PUBLIC ::  define_masks

  TYPE land_surface_type
     INTEGER              :: ntiles
     INTEGER, POINTER, DIMENSION(:,:)     :: & !! (nland, ntiles)
          cover_type                           !! Index into LctLibrary
     REAL(dp),    POINTER, DIMENSION(:,:) :: & !! (nland, ntiles)
          cover_type_real,                   & !! cover_type converted to REAL for netCDF files
          cover_fract,                       & !! Fraction of coverage for each land cover type
          cover_fract_pot,                   & !! Fraction of coverage for potential natural vegetation 
          albedo,                            & !! Surface albedo
          albedo_vis,                        & !! Surface albedo in the visible range
          albedo_nir,                        & !! Surface albedo in the NIR range
          veg_ratio,                         & !! This is the actual fraction of a tile covered by a vegetation canopy, i.e. this
                                               !! value includes the actual leaf area index (i.e. foliar projective cover)
          box_veg_ratio                        !! in contrast to veg_ratio box_veg_ratio takes into account actual cover fractions
     REAL(dp),    POINTER, DIMENSION(:,:) :: & !! (nland, month)
          veg_fract
     REAL(dp), POINTER, DIMENSION(:)      :: & !! (nland)
          veg_ratio_max,                     & !! maximal fraction of grid box covered by leaves, assuming infinite LAI
          rock_fract,                        & !! fraction of grid box not suitable for plants
          elevation,                         & !! Mean land grid cell elevation (m)
          oro_std_dev,                       & !! Standard deviation of orography
          forest_fract,                      & !! Forest fraction (TEMPORARY!)
          lake_fract,                        & !! Lake fraction (only used within HD model)
          land_fract,                        & !! Fraction of land in grid cell
          area,                              & !! Area of land fraction of grid cell (m^2)
          swdown_acc,                        & !! accumulated solar downward radiation
          swdown_reflect_acc                   !! accumulated reflected solar radiation
     LOGICAL, POINTER, DIMENSION(:,:) ::     & !! (nland,  ntiles)
          is_bare_soil,     &     !! Tile is bare soil (but not glacier)?
          is_vegetation,    &     !! Tile is vegetation (includes all natural and anthropogenically used types)
          is_C4vegetation,  &     !! Tile is vegetation with C4 photosynthesis
          is_naturalVeg,    &     !! Tile is natural vegetation
          is_forest,        &     !! Tile is forest
          is_grass,         &     !! Tile is covered with grassland          
          is_pasture,       &     !! Tile is covered with pasture
          is_crop,          &     !! Tile is covered with crops
          is_glacier,       &     !! Tile is glacier?
          is_lake,          &     !! Tile is lake?
          is_present              !! Is tile present, i.e. should it be handled by the land model?
  END TYPE land_surface_type
  PUBLIC ::  land_surface_type

  TYPE land_surface_diag_type
     REAL(dp), POINTER, DIMENSION(:)     :: & !! (nland)
          albedo
  END TYPE land_surface_diag_type

  TYPE albedo_options_type
     LOGICAL :: UseAlbedoCanopy
     LOGICAL :: UseAlbedoSoil
     LOGICAL :: UseAlbedoSoilConst
     CHARACTER(len=6) :: UseSOC
     LOGICAL :: UseLitter
  END TYPE albedo_options_type
  PUBLIC :: albedo_options_type
  TYPE(albedo_options_type), SAVE :: albedo_options

  TYPE albedo_params_type
     REAL(dp) :: AlbedoAgeWeight          !! 0: ECHAM5 scheme for snow albedo is used
                                          !! 1: snow age scheme is used
                                          !! 0 < albedo_age_weight< 1: snow albedo is calcualted by linearly weighting
                                          !! the snow albedo resulting from both schemes
  END TYPE albedo_params_type
  PUBLIC :: albedo_params_type
  TYPE(albedo_params_type), SAVE :: albedo_params

  !! The following parameter fract_small is the smallest allowed value for cover_fract and veg_ratio_max
  REAL(dp), PARAMETER :: fract_small = 1.e-10_dp       !! very small fraction (e.g. minimum value of cover_fract)

  !! Instead, the parameter box_fract_small below is used as an estimate of the smallest box-cover-fraction of PFTs -- this  
  !! should not be confused with fract_small, the smallest value for cover_fract, which refers to canopy-area instead of   
  !! box-area. The relation between these two is obtained from
  !!         box-cover-fraction=veg_ratio_max * (1-exp(-a*LAI_max)) * cover_fract. 
  !! Since the clumping term is even for a small LAI_max=0.5 larger than 0.1, and because veg_ratio_max and cover_fract are
  !! prepared in the initial files such that they are larger than fract_small, box-cover-fraction > 0.1*fract_small*fract_small. 
  !! This explains the choice of box_fract_small as
  real(dp),parameter :: box_fract_small = 0.1*fract_small*fract_small !! Estimate of smallest value for box-cover-fractions of PFTs

  PUBLIC :: fract_small,box_fract_small

  REAL(dp), PARAMETER :: initial_albedo = 0.2_dp       ! initial value for albedo in the whole solar range
  REAL(dp), PARAMETER :: initial_albedo_vis = 0.1_dp   ! initial value for albedo in the visible range
  REAL(dp), PARAMETER :: initial_albedo_nir = 0.3_dp   ! initial value for albedo in the near infrared
  REAL(dp), PARAMETER :: minimum_rad = 1.0e-09_dp      ! minimum amount of radiation for albedo calculations [W/m2]


  REAL(dp), POINTER, SAVE :: init_cover_fract(:,:)
  REAL(dp), POINTER, SAVE :: init_cover_fract_pot(:,:)
  REAL(dp), POINTER, SAVE :: init_veg_ratio_max(:)
  PUBLIC :: init_cover_fract, init_cover_fract_pot, init_veg_ratio_max

  PRIVATE

  INTEGER, SAVE :: nlct = -1        !! Number of land cover types
  INTEGER, SAVE :: ntiles = -1      !! Maximum number of land cover types actually used (<= nclt)

  TYPE(t_stream), POINTER :: IO_land_surface
  TYPE(t_stream), POINTER :: IO_diag
  TYPE(FILE_INFO), SAVE   :: land_surface_file

  TYPE(land_surface_diag_type), SAVE    :: land_surface_diag

CONTAINS

  SUBROUTINE land_surface_init_io(ntiles, nlct_help, IO_file_name)

  ! land_surface_init_io is called from init_land_surface

    USE mo_netCDF, ONLY: add_dim, IO_inq_dimid, IO_inq_dimlen, io_get_att_int, &
                         NF_MAX_NAME, NF_GLOBAL
    USE mo_io,  ONLY: IO_open, IO_READ, IO_close
    USE mo_linked_list, ONLY : TILES

    INTEGER, INTENT(in)  :: ntiles
    INTEGER, INTENT(in)  :: nlct_help
    CHARACTER(NF_MAX_NAME), INTENT(in) :: IO_file_name

    INTEGER :: IO_file_id, IO_dim_id
    INTEGER :: ntiles_help
    REAL(dp)  :: tile_values(ntiles)
    INTEGER :: i
    LOGICAL :: exist

    IF (p_parallel_io) THEN


       ! --- get number of landcover types actually used in grid boxes

       land_surface_file%opened = .FALSE.
       CALL IO_open(TRIM(IO_file_name), land_surface_file, IO_READ)
       IO_file_id = land_surface_file%file_id

       CALL IO_get_att_int(IO_file_id, NF_GLOBAL, 'nlct', nlct, exist)
       IF (exist) THEN
          CALL IO_get_att_int(IO_file_id, NF_GLOBAL, 'nlct', nlct)
       ELSE
          CALL IO_inq_dimid(IO_file_id, 'lct', IO_dim_id)
          CALL IO_inq_dimlen(IO_file_id, IO_dim_id, nlct)
       END IF
       CALL IO_inq_dimid(IO_file_id, 'ntiles', IO_dim_id)
       CALL IO_inq_dimlen(IO_file_id, IO_dim_id, ntiles_help)

       CALL IO_close(land_surface_file)

       IF (nlct_help   /= nlct)   CALL finish('land_surface_init_io', 'nlct from IO_file differs from the one from LctLibrary')
       IF (ntiles_help /= ntiles) CALL finish('land_surface_init_io', 'ntiles from IO_file differs from the one in run.def')

    END IF

    IF (p_parallel) THEN
       CALL p_bcast(nlct, p_io)
    END IF

    IF (debug) CALL message('land_surface_init_io','Adding dimensions')
    CALL add_dim("lct",nlct,"available land cover types")

    DO i=1, ntiles
      tile_values(i)=REAL(i,dp)
    END DO
    CALL add_dim ("tiles", ntiles, longname="land surface tile", units="", value=tile_values(:), levtyp=70, indx=TILES)

  END SUBROUTINE land_surface_init_io

  !=================================================================================================
  SUBROUTINE land_surface_init_memory(g_nland, l_nland, ntiles_help, land_surface, useDynveg, withHD, &
                                      isStandalone, fileformat, fileztype, diag_stream, stream)

  ! land_surface_init_memory is called from init_land_surface

    USE mo_jsbach,      ONLY : missing_value, lpost_echam
    USE mo_linked_list, ONLY : LAND, TILES
    USE mo_memory_base, ONLY : new_stream, default_stream_setting, &
                               add =>add_stream_element
    USE mo_netCDF,      ONLY : max_dim_name
    USE mo_output,      ONLY : land_table

    INTEGER,                 INTENT(in)    :: g_nland, l_nland, ntiles_help
    LOGICAL,                 INTENT(in)    :: useDynveg
    LOGICAL,                 INTENT(in)    :: withHD
    TYPE(land_surface_type), INTENT(inout) :: land_surface
    LOGICAL,                 INTENT(in)    :: isStandalone
    INTEGER,                 INTENT(in)    :: fileformat     ! output file format
    INTEGER,                 INTENT(in)    :: fileztype      ! output file compression
    TYPE(t_stream),          POINTER       :: diag_stream
    TYPE(t_stream), POINTER, OPTIONAL :: stream

    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim2p(2), dim2(2)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim2n(2)

    IF (debug) CALL message('land_surface_init_memory','Entering ...')

    IF (ASSOCIATED(diag_stream)) THEN
       IO_diag => diag_stream
    ELSE
       CALL finish('land_surface_init_memory', 'Diagnostic stream not present')
    END IF

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'land surface', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_land_surface => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_land_surface, 'land surface', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(IO_land_surface, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF


    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim2p = (/ l_nland, ntiles_help /)
    dim2  = (/ g_nland, ntiles_help /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'tiles'

    ! code numbers of this routine range from 10 to 29, using the GRIB land_table

    CALL add(IO_land_surface, 'cover_fract_pot',   land_surface%cover_fract_pot,longname='Potential Natural Land Cover Fraction',&
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=10, lrerun=.TRUE., &
             lpost=useDynveg,lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'cover_type',        land_surface%cover_type_real, longname='Land Cover', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=11, lrerun=.FALSE., lpost=.FALSE.)
    CALL add(IO_land_surface, 'cover_fract',       land_surface%cover_fract,     longname='Land Cover Fraction', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=12, lrerun=.TRUE.,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'albedo' ,           land_surface%albedo,          longname='Land Surface Albedo', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=13,  lpost=.FALSE.)
    CALL add(IO_land_surface, 'albedo_vis',        land_surface%albedo_vis,      longname='Surface Albedo in the Visible Range',  &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=14,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'albedo_nir' ,       land_surface%albedo_nir,      longname='Surface Albedo in the NIR', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=15,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'elevation',         land_surface%elevation,       longname='Surface Altitude', &
             units='m',    ldims=dim1p, gdims=dim1, dimnames=dim1n, code=16,  lpost=.FALSE.)
    CALL add(IO_land_surface, 'orography_std_dev', land_surface%oro_std_dev,     longname='Standard Deviation of the Orography', &
             units='m',    ldims=dim1p, gdims=dim1, dimnames=dim1n, code=17,  lpost=.FALSE.)
    CALL add(IO_land_surface, 'forest_fract',      land_surface%forest_fract,    longname='Forest Fraction', &
             units='',     ldims=dim1p, gdims=dim1, dimnames=dim1n, code=18,  lpost=.FALSE.)
    CALL add(IO_land_surface, 'lake_fract',       land_surface%lake_fract,    longname='Lake Fraction', &
             units='',     ldims=dim1p, gdims=dim1, dimnames=dim1n, code=26,  lrerun=withHD, lpost=.FALSE.)
    CALL add(IO_land_surface, 'veg_ratio',         land_surface%veg_ratio,       longname='vegetated fraction of grid box', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=19,  lpost=.FALSE.)
    CALL add(IO_land_surface, 'veg_ratio_max',     land_surface%veg_ratio_max,   longname='Maximum Vegetation Fraction', &
             units='',     ldims=dim1p, gdims=dim1, dimnames=dim1n, code=20, lpost=useDynveg, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'box_veg_ratio',     land_surface%box_veg_ratio,   longname='vegetated fraction of grid box', &
             units='',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=24,  lmiss=.TRUE., missval=missing_value, lrerun=.FALSE.)
    CALL add(IO_land_surface, 'land_fract',     land_surface%land_fract,      longname='Land Fraction', &
             units='',  ldims=dim1p, gdims=dim1, dimnames=dim1n, code=25, lpost=.NOT. isStandalone, contnorest=.TRUE.)
    IF (useDynveg) THEN
       CALL add(IO_land_surface, 'rock_fract',     land_surface%rock_fract,      longname='Rock Fraction', &
                units='',  ldims=dim1p, gdims=dim1, dimnames=dim1n, code=23, lpost=.FALSE., contnorest=.TRUE.)
    ENDIF
    CALL add(IO_land_surface, 'swdown_acc',        land_surface%swdown_acc,      longname='Surface Downwelling Solar Radiation', &
             units='W m-2',ldims=dim1p, gdims=dim1, dimnames=dim1n, code=21, laccu=.true., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_land_surface, 'swdown_reflect_acc',land_surface%swdown_reflect_acc,longname='Surface Upwelling Solar Radiation', &
             units='W m-2',ldims=dim1p, gdims=dim1, dimnames=dim1n, code=22, laccu=.true., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag,         'albedo',            land_surface_diag%albedo,     longname='Land Surface Albedo', &
             units='',     ldims=dim1p, gdims=dim1, dimnames=dim1n, code=13, laccu=.FALSE., lpost=lpost_echam, &
             lmiss=.TRUE., missval=missing_value)

    if (debug) CALL message('land_surface_init_memory','    exit')

    ALLOCATE(land_surface%cover_type(l_nland,ntiles_help))
    ALLOCATE(land_surface%is_bare_soil(l_nland, ntiles_help),    &
             land_surface%is_vegetation(l_nland, ntiles_help),   &
             land_surface%is_C4vegetation(l_nland, ntiles_help), &
             land_surface%is_naturalVeg(l_nland, ntiles_help),   &
             land_surface%is_forest(l_nland, ntiles_help),       &
             land_surface%is_grass(l_nland, ntiles_help),        &
             land_surface%is_pasture(l_nland, ntiles_help),      &
             land_surface%is_crop(l_nland, ntiles_help),         &
             land_surface%is_lake(l_nland, ntiles_help),         &
             land_surface%is_glacier(l_nland, ntiles_help),      &
             land_surface%is_present(l_nland, ntiles_help) )

  END SUBROUTINE land_surface_init_memory

  !=================================================================================================
  SUBROUTINE init_land_surface(grid, domain, lctlib, land_surface, surf_file,  &
       useDynveg, useLanduseTransitions, isStandalone, isRestart, &
       read_cover_fract, withHD, fileformat, fileztype, IO_diag_stream, IO_stream)

    USE mo_tr_scatter,    ONLY: scatter_field
    USE mo_io,            ONLY: IO_open, IO_READ
    USE mo_netcdf,        ONLY: io_inq_dimid, io_inq_dimlen, io_inq_varid, io_get_var_double, nf_max_name
    Use mo_temp                          ! Provides temporary arrays

    TYPE(grid_type),   INTENT(in)        :: grid
    TYPE(domain_type), INTENT(in)        :: domain
    TYPE(lctlib_type), INTENT(in)        :: lctlib
    TYPE(land_surface_type), INTENT(inout) :: land_surface
    CHARACTER(nf_max_name), INTENT(in)   :: surf_file
    LOGICAL,           INTENT(in)        :: useDynveg
    LOGICAL,           INTENT(in)        :: useLanduseTransitions
    LOGICAL,           INTENT(in)        :: isStandalone
    LOGICAL,           INTENT(in)        :: isRestart
    LOGICAL,           INTENT(in)        :: read_cover_fract ! read cover fractions from initial
                                                             ! instead of restart file
    LOGICAL,           INTENT(in)        :: withHD           ! HD model is active
    INTEGER,           INTENT(in)        :: fileformat       ! output file format
    INTEGER,           INTENT(in)        :: fileztype        ! output file compression
    TYPE(t_stream),    POINTER           :: IO_diag_stream
    TYPE(t_stream),    POINTER, OPTIONAL :: IO_stream

    TYPE(FILE_INFO) :: IO_file
    INTEGER  :: IO_file_id, IO_var_id, IO_dim_id
    INTEGER  :: i, znlon, znlat
    integer  :: status

    if(debug) call message("init_land_surface","Start initialization of mo_land_surface")

    ntiles = land_surface%ntiles  ! land_surface%ntiles is set in calling jsbach_init

    call land_surface_init_io(ntiles, lctlib%nlct, surf_file)  ! Sets nlct  and test whether lctlib%nlct is consistent with IO file

    ! --- generate land surface stream

    CALL land_surface_init_memory(grid%nland, domain%nland, ntiles, land_surface, useDynveg, withHD, isStandalone, &
         fileformat, fileztype, IO_diag_stream, stream=IO_stream)

    IF (.NOT. ASSOCIATED(IO_land_surface)) &
         CALL finish('init_land_surface','No memory stream for land surface')

    CALL init_albedo(albedo_params, albedo_options)

    IF (p_parallel_io) THEN

       ! Open ini file
       CALL message('init_land_surface','Reading land surface fields from '//TRIM(land_surface_file%file_name))
       IO_file%opened = .FALSE.
       CALL IO_open(TRIM(land_surface_file%file_name), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

       ! Check resolution
       CALL IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       CALL IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       CALL IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       CALL IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)

       IF (znlon /= grid%nlon .OR. znlat /= grid%nlat) THEN
          CALL finish('init_land_surface', 'Unexpected resolution:'//int2string(znlon)//' '//int2string(znlat))
       ENDIF

    ENDIF

    ! Temporary storage for local domain fields
    ALLOCATE(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('init_land_surface','Allocation failure (1)')

    ! Land surface cover types
    IF (p_parallel_io) THEN
       ALLOCATE(zreal3d(grid%nlon,grid%nlat,ntiles),STAT=status)
       if(status .ne. 0) call finish('init_land_surface','Allocation failure (2)')
       CALL IO_inq_varid(IO_file_id, 'cover_type', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
    END IF
    NULLIFY(zreal2d_ptr)
    DO i=1,ntiles
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_field(zreal2d_ptr, zreal2d)
       land_surface%cover_type_real(:,i) = PACK(zreal2d, MASK=domain%mask)
    END DO

    ! Convert cover_type from REAL to INTEGER
    land_surface%cover_type = NINT(land_surface%cover_type_real)

!! -- Read in landcover fractions from the initial file. 
!!    In restarted runs they are read from restart file, unless read_cover_fract is true (namelist jsbach_ctl)

    IF (.NOT. isRestart .OR. read_cover_fract ) THEN

       ! Land surface cover fractions
       ALLOCATE(init_cover_fract(domain%nland,ntiles))
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'cover_fract', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_cover_fract(:,i) = PACK(zreal2d, MASK=domain%mask)
       ENDDO

       ! Read veg_ratio_max
       NULLIFY(zreal2d_ptr)
       ALLOCATE(init_veg_ratio_max(domain%nland))
       IF (p_parallel_io) THEN
          ALLOCATE(zreal2d_ptr(grid%nlon,grid%nlat))
          CALL IO_inq_varid(IO_file_id, 'veg_ratio_max', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
       ENDIF
       CALL scatter_field(zreal2d_ptr, zreal2d)
       init_veg_ratio_max(:) = PACK(zreal2d, MASK=domain%mask)

    END IF

    !! === Read in map of potential vegetation in case landuse transitions shall be computed

    IF (useLanduseTransitions .OR. useDynveg) THEN
       ALLOCATE(init_cover_fract_pot(domain%nland,ntiles))
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'natural_veg', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       ENDIF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_cover_fract_pot(:,i) = PACK(zreal2d, MASK=domain%mask)
       ENDDO
    END IF

    IF (p_parallel_io) DEALLOCATE(zreal3d)

    ! ECHAM5 compatibility: read monthly vegetation ratio and use this as fraction x of first tile (vegetation), and
    ! use 1-x as fraction of second tile (bare soil). Only the third tile is kept from previous read (glacier)
    ALLOCATE(land_surface%veg_fract(domain%nland, 0:13))
    IF (p_parallel_io) THEN
       ALLOCATE(zreal3d(grid%nlon,grid%nlat,12), STAT=status)
       if(status .ne. 0) call finish('init_land_surface','Allocation failure')
       CALL IO_inq_varid(IO_file_id, 'veg_fract', IO_var_id)
       CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
    END IF
    NULLIFY(zreal2d_ptr)
    DO i=1,12
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_field(zreal2d_ptr, zreal2d)
       land_surface%veg_fract(:,i) = PACK(zreal2d, MASK=domain%mask)
    END DO
    land_surface%veg_fract(:,0) = land_surface%veg_fract(:,12)
    land_surface%veg_fract(:,13) = land_surface%veg_fract(:,1)

    ! If this is a restart run the model state variables are read in jsbach_init (Standalone) or the GCM
    ! by calling io_read_streams. We can therefore exit now:
    IF (isRestart) THEN
       NULLIFY(zreal2d_ptr)
       DEALLOCATE(zreal2d)
       IF (p_parallel_io) DEALLOCATE(zreal3d)
       RETURN
    ENDIF
    
    ! If this is not a restart run we continue and get the land surface model state from ini file

    land_surface%albedo = initial_albedo
    land_surface%albedo_vis = initial_albedo_vis
    land_surface%albedo_nir = initial_albedo_nir
    land_surface_diag%albedo = 0.0_dp
    land_surface%swdown_acc = 0.0_dp
    land_surface%swdown_reflect_acc = 0.0_dp

    ! Elevation
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) THEN
       ALLOCATE(zreal2d_ptr(grid%nlon,grid%nlat))
       CALL io_inq_varid (IO_file_id, 'elevation', IO_var_id)
       CALL io_get_var_double (IO_file_id, IO_var_id, zreal2d_ptr)
    ENDIF
    CALL scatter_field(zreal2d_ptr, zreal2d)
    land_surface%elevation = PACK(zreal2d, MASK=domain%mask)
    
    ! Standard deviation of orography
    IF (p_parallel_io) THEN
       CALL io_inq_varid(IO_file_id, 'orography_std_dev', IO_var_id)
       CALL io_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
    END IF
    CALL scatter_field(zreal2d_ptr, zreal2d)
    land_surface%oro_std_dev = PACK(zreal2d, MASK=domain%mask)

    ! Forest fraction (temporary until JSBACH is verified against ECHAM5)
    IF (p_parallel_io) THEN
       CALL io_inq_varid(IO_file_id, 'forest_fract', IO_var_id)
       CALL io_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
    END IF
    CALL scatter_field(zreal2d_ptr, zreal2d)
    land_surface%forest_fract = PACK(zreal2d, MASK=domain%mask)

    ! Lake fraction
    IF (withHD) THEN
       IF (p_parallel_io) THEN
          CALL io_inq_varid(IO_file_id, 'lake', IO_var_id)
          CALL io_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
       END IF
       CALL scatter_field(zreal2d_ptr, zreal2d)
       land_surface%lake_fract = PACK(zreal2d, MASK=domain%mask)
    END IF

    ! Grid cell area
    !ALLOCATE(grid%area(grid%nland))
    !IF (p_parallel_io) THEN
    !   CALL io_inq_varid (IO_file_id, 'area', IO_var_id)
    !   CALL io_get_var_double (IO_file_id, IO_var_id, zreal2d_ptr)
    !ENDIF
    !CALL scatter_field(zreal2d_ptr, zreal2d)
    !grid%area = PACK(zreal2d, MASK=grid%mask_land)

    IF (p_parallel_io) THEN
       DEALLOCATE(zreal2d_ptr)
       DEALLOCATE(zreal3d)
    END IF
    NULLIFY(zreal2d_ptr)
    DEALLOCATE(zreal2d)

    land_surface%land_fract(:) = 1._dp

    IF (useDynveg) land_surface%rock_fract(:) = 0._dp

  END SUBROUTINE init_land_surface

  !----------------------------------------------------------------------------
  SUBROUTINE define_masks(land_surface, lctlib)

    ! Routine to define various masks from land cover distribution.
    ! 
    ! As the cover_fractions are not yet available in init_land_surface,
    ! these masks have to be defined separately.
    !

    TYPE(land_surface_type), INTENT(inout) :: land_surface
    TYPE(lctlib_type),       INTENT(in)    :: lctlib

    INTEGER :: ilct

    land_surface%is_bare_soil   = .FALSE.
    land_surface%is_vegetation  = .FALSE.
    land_surface%is_C4vegetation  = .FALSE.
    land_surface%is_naturalVeg  = .FALSE.
    land_surface%is_forest      = .FALSE.
    land_surface%is_grass       = .FALSE.
    land_surface%is_pasture     = .FALSE.
    land_surface%is_crop        = .FALSE.
    land_surface%is_lake        = .FALSE.
    land_surface%is_glacier     = .FALSE.
    land_surface%is_present     = .FALSE.

    DO ilct=1,nlct
       WHERE (land_surface%cover_type == ilct .AND. land_surface%cover_fract > 0._dp)
          land_surface%is_vegetation   = lctlib%NaturalVegFlag(ilct) .OR. lctlib%CropFlag(ilct) .OR. lctlib%PastureFlag(ilct)
          land_surface%is_C4vegetation = lctlib%C4flag(ilct) 
          land_surface%is_naturalVeg   = lctlib%NaturalVegFlag(ilct)
          land_surface%is_forest       = lctlib%ForestFlag    (ilct)
          land_surface%is_grass        = lctlib%GrassFlag     (ilct)   
          land_surface%is_pasture      = lctlib%PastureFlag   (ilct)
          land_surface%is_crop         = lctlib%CropFlag      (ilct)
          land_surface%is_lake         = lctlib%LakeFlag      (ilct)
          land_surface%is_glacier      = lctlib%GlacierFlag   (ilct)
          land_surface%is_bare_soil    = lctlib%BareSoilFlag  (ilct)
       END WHERE
    END DO

    ! Logical mask for tiles that should be processed by the land model
    land_surface%is_present    = land_surface%is_vegetation .OR. land_surface%is_bare_soil .OR. &
                                 land_surface%is_lake       .OR. land_surface%is_glacier

    ! In the current model version no fractional glacier cells are allowed.  
    IF (ANY(land_surface%is_glacier(:,:) .AND. &
         land_surface%cover_fract(:,:) > fract_small + EPSILON(1._dp) .AND. &
         land_surface%cover_fract(:,:) < 1._dp - EPSILON(1._dp))) THEN
       CALL finish('init_land_surface','Fractional glacier cells are not allowed in the current model version')
    END IF

  END SUBROUTINE define_masks

  !=================================================================================================

  SUBROUTINE init_albedo(albedo_params, albedo_options)

    TYPE(albedo_params_type), INTENT(inout) :: albedo_params
    TYPE(albedo_options_type), INTENT(inout) :: albedo_options

    !! Read namelist and lctlib-file for albedo

    CALL config_albedo(albedo_params, albedo_options)
  END SUBROUTINE init_albedo

  !=================================================================================================

  SUBROUTINE config_albedo(albedo_params, albedo_options)

    USE mo_namelist,         ONLY: position_nml, POSITIONED
    USE mo_jsbach,           ONLY: nml_unit
    USE mo_io_units,         ONLY: nout

    TYPE(albedo_params_type), INTENT(inout) :: albedo_params
    TYPE(albedo_options_type), INTENT(inout) :: albedo_options

    !! Locals
    INTEGER :: read_status, f_unit

    !! Namelist Parameters
    LOGICAL :: use_albedocanopy   !! true: use map of canopy albedo
                                  !! false: use PFT specific albedo values
    LOGICAL :: use_albedosoil     !! true: calculate albedo of soil surface denpending on soil carbon and litter
                                  !! false: albedo of the soil surface as read from the jsbach.nc file is used
    LOGICAL :: use_albedosoilconst!! true: the base albedo of the soil surface (without soil carbon and leaf litter)
                                  !!       is set to a global constant
                                  !! false: the base albedo of the soil surface as read from the jsbach.nc file is used
    CHARACTER(len=6)  :: use_soc  !! linear: albedo of the soil is linearly reduced by soil carbon
                                  !! log: logarithmic dependence of the albedo of the soil surface on soil carbon
    LOGICAL :: use_litter         !! true: albedo of the soil surface depends on leaf litter
    REAL(dp) :: albedo_age_weight !! 0: ECHAM5 scheme for snow albedo is used
                                  !! 1: snow age scheme is used
                                  !! 0 < albedo_age_weight< 1: snow albedo is calcualted by linearly weighting
                                  !! the snow albedo resulting from both schemes
    INCLUDE 'albedo_ctl.inc'

    !! Read namelist albedo_ctl

    IF (p_parallel_io) THEN

       ! define default values
       use_albedocanopy = .FALSE.
       use_albedosoil = .FALSE.
       use_albedosoilconst = .FALSE.
       use_soc = 'linear'
       use_litter = .TRUE.
       albedo_age_weight = 0.5_dp

       f_unit = position_nml ('ALBEDO_CTL', nml_unit, status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (f_unit, albedo_ctl)
          CALL message('config_albedo', 'Namelist ALBEDO_CTL: ')
          WRITE(nout, albedo_ctl)
       END SELECT
    END IF

    IF (p_parallel_io) THEN
       WRITE (message_text,*) 'use_albedocanopy: ', use_albedocanopy
       CALL message('config_albedo', message_text)
       albedo_options%UseAlbedoCanopy = use_albedocanopy

       WRITE (message_text,*) 'use_albedosoil: ', use_albedosoil
       CALL message('config_albedo', message_text)
       albedo_options%UseAlbedoSoil = use_albedosoil

       WRITE (message_text,*) 'use_albedosoilconst: ', use_albedosoilconst
       CALL message('config_albedo', message_text)
       albedo_options%UseAlbedoSoilConst = use_albedosoilconst

       WRITE (message_text,*) 'use_soc: ', TRIM(use_soc)
       CALL message('config_albedo', message_text)
       albedo_options%UseSOC = TRIM(use_soc)

       WRITE (message_text,*) 'use_litter: ', use_litter
       CALL message('config_albedo', message_text)
       albedo_options%UseLitter = use_litter

       WRITE (message_text,*) 'albedo_age_weight: ', albedo_age_weight
       CALL message('config_albedo', message_text)
       albedo_params%AlbedoAgeWeight = albedo_age_weight
    END IF

    IF (p_parallel) THEN
       CALL p_bcast(albedo_options%UseAlbedoCanopy, p_io)
       CALL p_bcast(albedo_options%UseAlbedoSoil, p_io)
       CALL p_bcast(albedo_options%UseAlbedoSoilConst, p_io)
       CALL p_bcast(albedo_options%UseSOC, p_io)
       CALL p_bcast(albedo_options%UseLitter, p_io)
       CALL p_bcast(albedo_params%AlbedoAgeWeight, p_io)
    END IF
 
  END SUBROUTINE config_albedo

  !=================================================================================================

  SUBROUTINE update_land_surface_fast(kidx, lctlib, land_surface, &
          surface_temperature, snow_fract, background_albedo, canopy_snow_fract, lai )

! This routine is coded for ECHAM5 compatibility!

    USE mo_jsbach_constants, ONLY: tmelt

    INTEGER, INTENT(in) :: kidx
    TYPE(lctlib_type),       INTENT(in)    :: lctlib
    TYPE(land_surface_type), INTENT(inout) :: land_surface
    REAL(dp), INTENT(in), DIMENSION(kidx, ntiles) ::   &
         surface_temperature, &
         snow_fract,          &                    ! Fraction of snow covered ground
         background_albedo,   &                    ! background albedo
         canopy_snow_fract,   &                    ! Fraction of snow covered canopy
         lai

    ! Local variables
    REAL(dp)    :: min_temp_snow_albedo             ! Temperature threshold below which maximum snow albedo is used
    REAL(dp)    :: snow_albedo_ground(kidx,ntiles)  ! Temperature dependend snow albedo over ground
    REAL(dp)    :: sky_view_fract(kidx,ntiles)      ! Fraction of bare ground below canopy for albedo calculation
    REAL(dp)    :: forest_fract(kidx,ntiles)        ! Factor for albedo from canopy (uses forest fraction and sky view factor)
    INTEGER :: kidx0, kidx1, i, j, ilct, itile
    REAL(dp), PARAMETER :: SkyViewFactor = 1.0_dp        !! Constant in albedo calculation
    REAL(dp), PARAMETER :: AlbedoCanopySnow = 0.20_dp    !! Albedo of snow covered canopy
!!$    INTEGER :: ctype(kidx)

    kidx0 = kstart
    kidx1 = kend

    min_temp_snow_albedo = tmelt - 5.0_dp
    snow_albedo_ground = 0._dp

    do itile=1,ntiles
      do i=1,kidx
        j=kidx0+i-1
        ilct = land_surface%cover_type(j,itile)
        ! Temperature dependend snow albedo over ground
        IF (surface_temperature(i,itile) >= tmelt) THEN
          snow_albedo_ground(i,itile) = lctlib%AlbedoSnowMin(ilct)
        ELSE IF (surface_temperature(i,itile) < min_temp_snow_albedo) THEN
          snow_albedo_ground(i,itile) = lctlib%AlbedoSnowMax(ilct)
        ELSE
         snow_albedo_ground(i,itile) = lctlib%AlbedoSnowMin(ilct) + &
                     (tmelt - surface_temperature(i,itile)) * (lctlib%AlbedoSnowMax(ilct) - lctlib%AlbedoSnowMin(ilct)) / &
                     (tmelt - min_temp_snow_albedo)
        END IF
      END DO
    ENDDO

    sky_view_fract(:,:) = EXP(-SkyViewFactor * MAX(lai(:,:),2.0_dp)) ! the value 2.0 should be replaced by StemArea

    ! calculate fraction for which albedo is computed from canopy
    forest_fract(:,:) = SPREAD(land_surface%forest_fract(kidx0:kidx1), DIM=2, NCOPIES=ntiles) * (1._dp - sky_view_fract(:,:))

    WHERE (land_surface%is_glacier(kidx0:kidx1,:))
       land_surface%albedo(kidx0:kidx1,:) = snow_albedo_ground
    ELSEWHERE
       ! albedo = weighted mean of albedo of ground below canopy and albedo of canopy
       land_surface%albedo(kidx0:kidx1,:) = &
            MAX(((1._dp - forest_fract) * (snow_fract * snow_albedo_ground + (1._dp - snow_fract) * background_albedo) + &
            forest_fract * (canopy_snow_fract * AlbedoCanopySnow + (1._dp - canopy_snow_fract) * background_albedo)),    &
            background_albedo)
    END WHERE

  END SUBROUTINE update_land_surface_fast

  !=================================================================================================

  ! albedo interface
  ! Thomas Raddatz, MPI for Meteorology, 2010
  SUBROUTINE update_albedo(lctlib,                                     &
          nidx, ntiles, cover_type, cover_fract, veg_fract_correction, &
          is_glacier, is_forest,                                       &
          veg_ratio_max, radiation_net_vis, radiation_net_nir,         &
          snow_age, surface_temperature, snow_fract,                   &
          Cpool_litter_green_ag, Cpool_slow,                           &
          background_albedo_soil_vis_base,                             &
          background_albedo_soil_nir_base,                             &
          background_albedo_veg_vis, background_albedo_veg_nir,        &
          lai, canopy_snow_fract,                                      &
          albedo_vis, albedo_nir)

    TYPE(lctlib_type),       INTENT(in)    :: lctlib
    INTEGER, INTENT(in)  ::  nidx                   ! domain
    INTEGER, INTENT(in)  ::  ntiles                 ! number of tiles
    INTEGER, INTENT(in), DIMENSION(nidx,ntiles) ::          &
         cover_type                                 ! id number of plant functional type (PFT)
    REAL(dp), INTENT(in), DIMENSION(nidx,ntiles)::          &
         cover_fract,                  &            ! cover fractions of the tiles
         veg_fract_correction                       !
    LOGICAL, INTENT(in), DIMENSION(nidx,ntiles) ::          &
         is_glacier,                   &            ! glacier flag
         is_forest                                  ! forest flag
    REAL(dp), INTENT(in), DIMENSION(nidx)  ::                   &
         veg_ratio_max,                &            ! maximal fraction of the grid boy covered by vegetation
         radiation_net_vis,            &            ! net solar radiation in the visible range [W/m2]
         radiation_net_nir,            &            ! net solar radiation in the NIR range [W/m2]
         snow_age                                   ! non-dimensional age of snow
    REAL(dp), INTENT(in), DIMENSION(nidx,ntiles)  ::            &
         surface_temperature,          &            !
         snow_fract,                   &            ! fraction of snow covered ground
         Cpool_litter_green_ag,        &            ! leaf litter carbon [mol(C)/m2(canopy)]
         Cpool_slow,                   &            ! soil carbon [mol(C)/m2(canopy)]
         background_albedo_soil_vis_base, &         ! background albedo soil visible (without litter and soil carbon)
         background_albedo_soil_nir_base, &         ! background albedo soil NIR (without litter and soil carbon)
         background_albedo_veg_vis,    &            ! background albedo vegetation visible
         background_albedo_veg_nir,    &            ! background albedo vegetation NIR
         lai,                          &            ! leaf area index
         canopy_snow_fract                          ! fraction of snow covered canopy (forest)
    REAL(dp), INTENT(inout), DIMENSION(nidx,ntiles)  ::         &
         albedo_vis,                   &            ! albedo of the visible range
         albedo_nir                                 ! albedo of the NIR range

    ! local variables
    INTEGER   :: i, itile
    INTEGER   :: well_mixed_tiles                         ! number of "well mixed" tiles
    REAL(dp)  :: cover_fract_mixed_sum                    ! sum of cover fractions of tiles that are "well mixed"
    REAL(dp)  :: cover_fract_mixed(nidx,ntiles)           ! cover fractions with respect to all surfaces that are "well mixed"
                                                          ! (currently only glaciers are not "well mixed", in future this should also
                                                          ! be lakes etc.)
    REAL(dp)  :: boxC_litter_green_ag_mixed(nidx,ntiles)  ! leaf litter carbon [mol(C)/m2(well mixed area)]
    REAL(dp)  :: boxC_slow_mixed(nidx,ntiles)             ! soil carbon [mol(C)/m2(well mixed area)]
    REAL(dp)  :: background_albedo_soil_vis(nidx,ntiles)  ! background albedo soil visible (including litter and soil carbon)
    REAL(dp)  :: background_albedo_soil_nir(nidx,ntiles)  ! background albedo soil NIR (including litter and soil carbon)

    !Calculate the albedo of the soil surface (or take it from the init_file)
    IF (albedo_options%UseAlbedoSoil) THEN
       ! calculate cover_fract_mixed
       DO i = 1,nidx
          well_mixed_tiles = 0
          cover_fract_mixed_sum = 0.0_dp
          DO itile = 1,ntiles
             IF (.NOT. is_glacier(i,itile)) THEN
                well_mixed_tiles = well_mixed_tiles + 1
                cover_fract_mixed_sum = cover_fract_mixed_sum + cover_fract(i,itile)
             END IF
          END DO
          DO itile = 1,ntiles
             IF (.NOT. is_glacier(i,itile) .AND. cover_fract_mixed_sum > EPSILON(1.0_dp)) THEN
                cover_fract_mixed(i,itile) = cover_fract(i,itile) / cover_fract_mixed_sum
             ELSE IF (.NOT. is_glacier(i,itile)) THEN
                cover_fract_mixed(i,itile) = 1.0_dp / REAL(well_mixed_tiles)
             ELSE
                cover_fract_mixed(i,itile) = 0.0_dp
             END IF
          END DO          
       END DO
       ! calculate soil carbon and leaf litter carbon [mol(C)/m2(well mixed area)] 
       DO i = 1,nidx
          DO itile = 1,ntiles
             IF (.NOT. is_glacier(i,itile)) THEN
                boxC_litter_green_ag_mixed(i,itile) = Cpool_litter_green_ag(i,itile) * cover_fract_mixed(i,itile) * &
                                                       veg_ratio_max(i) * veg_fract_correction(i,itile)
                boxC_slow_mixed(i,itile) = Cpool_slow(i,itile) * cover_fract_mixed(i,itile) * &
                                                       veg_ratio_max(i) * veg_fract_correction(i,itile)
             ELSE
                boxC_litter_green_ag_mixed(i,itile) = 0.0_dp
                boxC_slow_mixed(i,itile) = 0.0_dp
             END IF
          END DO
       END DO
       ! calculate albedo of the soil surface depending on soil carbon and leaf litter carbon
       CALL update_soil_albedo(lctlib,          & 
          nidx, ntiles, is_glacier, cover_type, &
          boxC_litter_green_ag_mixed,           &
          boxC_slow_mixed,                      &
          background_albedo_soil_vis_base,      &
          background_albedo_soil_nir_base,      &
          background_albedo_soil_vis,           &
          background_albedo_soil_nir)
    ELSE
       background_albedo_soil_vis(:,:) = background_albedo_soil_vis_base(:,:)
       background_albedo_soil_nir(:,:) = background_albedo_soil_nir_base(:,:)
    END IF

    ! calculate albedo_vis and albedo_nir depending on snow age and temperature
    CALL update_albedo_snowage_temp(lctlib,                       &
       nidx, ntiles, cover_type, is_glacier, is_forest,           &
       veg_ratio_max, radiation_net_vis, radiation_net_nir,       &
       snow_age,                                                  &
       surface_temperature, snow_fract,                           &
       background_albedo_soil_vis, background_albedo_soil_nir,    &
       background_albedo_veg_vis, background_albedo_veg_nir,      &
       lai, canopy_snow_fract,                                    &
       albedo_vis, albedo_nir)

  END SUBROUTINE update_albedo

  !=================================================================================================

  ! background soil albedo routine
  ! Freja Vambourg, MPI for Meteorology, 2009
  ! Thomas Raddatz, MPI for Meteorology, 2013
  SUBROUTINE update_soil_albedo(lctlib,                  &
            nidx, ntiles, is_glacier, cover_type,        &
            boxC_litter_green_ag, boxC_slow,             &
            albedo_soil_vis_base, albedo_soil_nir_base,  &
            albedo_soil_vis, albedo_soil_nir)
   
    TYPE(lctlib_type),       INTENT(in)    :: lctlib
    INTEGER, INTENT(in)  ::  nidx                   ! domain
    INTEGER, INTENT(in)  ::  ntiles                 ! number of tiles
    LOGICAL, INTENT(in), DIMENSION(nidx,ntiles)  ::             &
         is_glacier                                 ! glacier flag
    INTEGER, INTENT(in), DIMENSION(nidx,ntiles) ::              &
         cover_type                                 ! id number of plant functional type (PFT)
    REAL(dp), INTENT(in), DIMENSION(nidx,ntiles)  ::            &
         boxC_litter_green_ag,   &                  ! leaf litter carbon [mol(C)/m2(glacier free grid box)]
         boxC_slow,              &                  ! soil carbon [mol(C)/m2(glacier free grid box)]
         albedo_soil_vis_base,   &                  ! base albedo soil visible
         albedo_soil_nir_base                       ! base albedo soil NIR
    REAL(dp), INTENT(out), DIMENSION(nidx,ntiles)  ::           &
         albedo_soil_vis,        &                  ! background albedo soil visible 
         albedo_soil_nir                            ! background albedo soil NIR

    ! parameters for spatially constant albedo of the soil surface (UseAlbedoSoilConst)
    REAL(dp), PARAMETER :: AlbedoSoilConstVis = 0.2_dp ! 
    REAL(dp), PARAMETER :: AlbedoSoilConstNir = 0.2_dp ! 
    ! parameters for albedo dependence on soil carbon (UseSOC)
    REAL(dp), PARAMETER :: AlbedoSoilSOCVis = 3.e-4_dp ! 
    REAL(dp), PARAMETER :: AlbedoSoilSOCNir = 3.e-4_dp ! 
    REAL(dp), PARAMETER :: AlbedoSoilClimit = 500._dp  ! 

    ! local variables
    INTEGER  ::  i, itile
    REAL(dp)  :: albedo_desert_dry_vis(nidx,ntiles) ! desert soil albedo in a dry state 
    REAL(dp)  :: albedo_desert_dry_nir(nidx,ntiles) ! desert soil albedo in a dry state 
    REAL(dp)  :: albedo_litter_vis(nidx)            ! albedo of the leaf litter layer (visible)
    REAL(dp)  :: albedo_litter_nir(nidx)            ! albedo of the leaf litter layer (NIR)
    REAL(dp)  :: soil_carbon_sum(nidx)              ! soil carbon of all PFT [mol(C)/m2(glacier free grid box)]
    REAL(dp)  :: litter_area_sum(nidx)
    REAL(dp)  :: litter_area, albedo_litter_vis_sum, albedo_litter_nir_sum, litter_view_fract

    !!define the desert_soil_albedo
    albedo_desert_dry_vis(:,:) = albedo_soil_vis_base(:,:)
    albedo_desert_dry_nir(:,:) = albedo_soil_nir_base(:,:)

    !!option: all background soil albedo values are set to the same value
    IF (albedo_options%UseAlbedoSoilConst) THEN
       WHERE (.NOT. is_glacier(:,:))
          albedo_desert_dry_vis(:,:) = AlbedoSoilConstVis
          albedo_desert_dry_nir(:,:) = AlbedoSoilConstNir
       END WHERE
    END IF

    !!option: albedo depends on soil carbon content ...
    soil_carbon_sum(:) = 0.0_dp
    soil_carbon_sum(:) = SUM(boxC_slow(:,:), DIM=2)
    DO i = 1,nidx
       IF(albedo_options%UseSOC .eq. 'linear') THEN
          WHERE (is_glacier(i,:))
             albedo_soil_vis(i,:) = albedo_desert_dry_vis(i,:)
             albedo_soil_nir(i,:) = albedo_desert_dry_nir(i,:)
          ELSEWHERE
             albedo_soil_vis(i,:) = albedo_desert_dry_vis(i,:) -          &
                                    AlbedoSoilSOCVis * MIN(soil_carbon_sum(i),AlbedoSoilClimit)
             albedo_soil_nir(i,:) = albedo_desert_dry_nir(i,:) -          &
                                    AlbedoSoilSOCNir * MIN(soil_carbon_sum(i),AlbedoSoilClimit)
          ENDWHERE
       ELSE IF(albedo_options%UseSOC .eq. 'log') THEN
          WHERE (is_glacier(i,:))
             albedo_soil_vis(i,:) = albedo_desert_dry_vis(i,:)
             albedo_soil_nir(i,:) = albedo_desert_dry_nir(i,:)
          ELSEWHERE
             albedo_soil_vis(i,:) = MAX(6._dp * albedo_desert_dry_vis(i,:) / &
                                    (6._dp + log(1._dp + soil_carbon_sum(i))),0.1_dp)
             albedo_soil_nir(i,:) = MAX(6._dp * albedo_desert_dry_nir(i,:) / &
                                    (6._dp + log(1._dp + soil_carbon_sum(i))),0.2_dp)
          ENDWHERE
       ELSE
          albedo_soil_vis(i,:) = albedo_desert_dry_vis(i,:)
          albedo_soil_nir(i,:) = albedo_desert_dry_nir(i,:)
       END IF ! UseSOC
    END DO

    !! calculate the dependence of albedo on soil moisture (to be implemented later)
 
    ! option: albedo depends on litter
    IF(albedo_options%UseLitter) THEN
       DO i = 1,nidx
          litter_area_sum(i) = 0._dp
          albedo_litter_vis_sum = 0._dp
          albedo_litter_nir_sum = 0._dp
          DO itile = 1,ntiles
             IF (.NOT. is_glacier(i,itile)) THEN
                litter_area = boxC_litter_green_ag(i,itile) * lctlib%specificLeafArea_C(cover_type(i,itile))
                litter_area_sum(i) = litter_area_sum(i) + litter_area
                albedo_litter_vis_sum = albedo_litter_vis_sum + &
                                        litter_area * lctlib%AlbedoLitterVIS(cover_type(i,itile))
                albedo_litter_nir_sum = albedo_litter_nir_sum + &
                                        litter_area * lctlib%AlbedoLitterNIR(cover_type(i,itile))
             END IF
          END DO
          IF (litter_area_sum(i) > EPSILON(1._dp)) THEN
             albedo_litter_vis(i) = albedo_litter_vis_sum / litter_area_sum(i)
             albedo_litter_nir(i) = albedo_litter_nir_sum / litter_area_sum(i)
          ELSE
             albedo_litter_vis(i) = 0.1_dp
             albedo_litter_nir(i) = 0.2_dp
          END IF
       END DO

       DO i = 1,nidx
          !calculate the litter_view_fraction, i.e. what fraction of the box is covered by litter 
          !this determines how much soil is seen through litter 
          litter_view_fract = 1.0_dp - EXP(-0.5_dp * litter_area_sum(i))
          DO itile = 1,ntiles
             IF (.NOT. is_glacier(i,itile)) THEN
                albedo_soil_vis(i,itile) = (1.0_dp - litter_view_fract) * albedo_soil_vis(i,itile) + &
                                           litter_view_fract * albedo_litter_vis(i)
                albedo_soil_nir(i,itile) = (1.0_dp - litter_view_fract) * albedo_soil_nir(i,itile) + &
                                           litter_view_fract * albedo_litter_nir(i)
             END IF
          END DO
       END DO
    END IF ! UseLitter

  END SUBROUTINE update_soil_albedo


  !=================================================================================================
  ! give initial albedo values
  !   replaces former call of update_albedo during the inquiring time step.
  !
  SUBROUTINE init_albedo_land (albedo_vis, albedo_nir, albedo)

    REAL(dp), INTENT(out) ::  albedo_vis(:)        ! albedo of the visible range
    REAL(dp), INTENT(out) ::  albedo_nir(:)        ! albedo of the NIR range
    REAL(dp), INTENT(out) ::  albedo(:)            ! albedo of the whole solar range

    albedo(:)     = initial_albedo
    albedo_vis(:) = initial_albedo_vis
    albedo_nir(:) = initial_albedo_nir


  END SUBROUTINE init_albedo_land
  !=================================================================================================
  !
  ! main albedo routine
  ! Calculation of the albedo of the visible range (albedo_vis) and the albedo of the NIR range 
  ! (albedo_nir).
  ! The albedo of snow depends on surface temperature (ECHAM5 scheme) or snow age (BATS) or a
  ! combination of these two schemes (with the weighting factor AlbedoAgeWeight).
  ! Thomas Raddatz, MPI for Meteorology, 2013
  !
  SUBROUTINE update_albedo_snowage_temp(lctlib,                      &
          nidx, ntiles, cover_type, is_glacier, is_forest,           &
          veg_ratio_max, radiation_net_vis, radiation_net_nir,       &
          snow_age, surface_temperature, snow_fract,                 &
          background_albedo_soil_vis, background_albedo_soil_nir,    &
          background_albedo_veg_vis, background_albedo_veg_nir,      &
          lai, canopy_snow_fract,                                    &
          albedo_vis, albedo_nir)

    USE mo_jsbach_constants, ONLY: tmelt
#if defined (__SX__) && defined (_OPENMP)
    USE omp_lib,          ONLY: omp_get_thread_num, omp_get_num_threads
#endif
    USE mo_zenith, ONLY: cos_zenith

    TYPE(lctlib_type),       INTENT(in)    :: lctlib
    INTEGER, INTENT(in)  ::  nidx                   ! domain
    INTEGER, INTENT(in)  ::  ntiles                 ! number of tiles
    INTEGER, INTENT(in), DIMENSION(nidx,ntiles) ::          &
         cover_type                                 ! number of plant functional type (PFT)
    LOGICAL, INTENT(in), DIMENSION(nidx,ntiles) ::          &
         is_glacier,                   &            ! glacier flag
         is_forest                                  ! forest flag
    REAL(dp), INTENT(in), DIMENSION(nidx)  ::                   &
         veg_ratio_max,                &            ! maximal fraction of the grid boy covered by vegetation
         radiation_net_vis,            &            ! net solar radiation in the visible+UV range [W/m2]
         radiation_net_nir,            &            ! net solar radiation in the NIR range [W/m2]
         snow_age                                   ! non-dimensional age of snow
    REAL(dp), INTENT(in), DIMENSION(nidx,ntiles)  ::            &
         surface_temperature,          &            !
         snow_fract,                   &            ! fraction of snow covered ground
         background_albedo_soil_vis,   &            ! albedo of the soil surface in visible+UV range
         background_albedo_soil_nir,   &            ! albedo of the soil surface in NIR range
         background_albedo_veg_vis,    &            ! albedo of canopy covered surface in visible+UV range
         background_albedo_veg_nir,    &            ! albedo of canopy covered surface NIR range
         lai,                          &            ! leaf area index
         canopy_snow_fract                          ! fraction of snow covered canopy (forest)
    REAL(dp), INTENT(inout), DIMENSION(nidx,ntiles)  ::         &
         albedo_vis,                   &            ! albedo of the visible range
         albedo_nir                                 ! albedo of the NIR range
 
    ! parameters snow masking
    REAL(dp), PARAMETER :: SkyViewFactor = 0.5_dp          !! constant in calculating the sky view fraction depending on
                                                           !! lai and stem area

    ! parameters for glaciers (values from ECHAM5)
    REAL(dp), PARAMETER :: AlbedoGlacierVisMin = 0.78_dp   !! albedo of glacier in the visible+UV range at the melting point
    REAL(dp), PARAMETER :: AlbedoGlacierVisMax = 0.9_dp    !! albedo of glacier in the visible+UV range at hard frost 
    REAL(dp), PARAMETER :: AlbedoGlacierNirMin = 0.4_dp    !! albedo of glacier in the NIR range at at the melting point
    REAL(dp), PARAMETER :: AlbedoGlacierNirMax = 0.75_dp   !! albedo of glacier in the NIR range at hard frost
    REAL(dp), PARAMETER :: TempAlbedoGlacierMax = 5.0_dp   !! maximum glacier albedo at this temperature below melting point of H2O

    ! parameters of snow age scheme
    REAL(dp), PARAMETER :: AlbedoCanopySnow_age = 0.20_dp  !! albedo of snow covered canopy for snow age scheme (BATS)
    REAL(dp), PARAMETER :: AlbedoSnowVisMax = 0.90_dp      !! maximum albedo of fresh snow in the visible+UV range
    REAL(dp), PARAMETER :: AlbedoSnowNirMax = 0.60_dp      !! maximum albedo of fresh snow in the NIR range
    REAL(dp), PARAMETER :: AlbedoSnowVisAge = 0.15_dp      !! maximal rel. reduction of snow albedo by aging in the visible+UV range
    REAL(dp), PARAMETER :: AlbedoSnowNirAge = 0.5_dp       !! maximal rel. reduction of snow albedo by aging in the NIR range
    REAL(dp), PARAMETER :: AlbedoSnowAngle = 0.4_dp        !! maximal rel. reduction of absorption of solar radiation at large solar
                                                           !!    zenith angle
    REAL(dp), PARAMETER :: ZenithAngleFactor = 2._dp       !! factor in solar zenith angle dependence of snow albedo
                                                           !! (the increase of snow albedo is the higher this factor

    ! parameters of snow temperature (ECHAM5) scheme
    REAL(dp), PARAMETER :: AlbedoCanopySnow_temp = 0.20_dp !! Albedo of snow covered canopy for ECHAM5 scheme
    REAL(dp), PARAMETER :: TempAlbedoSnowMax = 5.0_dp      !! Maximum snow albedo at this temperature below melting point of H2O

    ! Local variables
    REAL(dp)  ::  background_albedo_vis                ! albedo without snow in the visible range
    REAL(dp)  ::  background_albedo_nir                ! albedo without snow in the NIR range
    REAL(dp)  ::  background_albedo_canopy_vis         ! albedo without snow of the canopy in the visible range
    REAL(dp)  ::  background_albedo_canopy_nir         ! albedo without snow of the canopy in the NIR range
    REAL(dp)  ::  snow_albedo_age_vis(nidx,ntiles)     ! albedo of snow (covering the soil) snow age scheme, visible range
    REAL(dp)  ::  snow_albedo_age_nir(nidx,ntiles)     ! albedo of snow (covering the soil) snow age scheme, NIR range
    REAL(dp)  ::  snow_albedo_temp_vis(nidx,ntiles)    ! albedo of snow (covering the soil) snow temperature scheme, visible range
    REAL(dp)  ::  snow_albedo_temp_nir(nidx,ntiles)    ! albedo of snow (covering the soil) snow temperature scheme, NIR range
    REAL(dp)  ::  snow_albedo_vis(nidx,ntiles)         ! albedo of snow (covering the soil), visible range
    REAL(dp)  ::  snow_albedo_nir(nidx,ntiles)         ! albedo of snow (covering the soil), NIR range
    REAL(dp)  ::  snow_age_factor                      ! snow aging factor
    REAL(dp)  ::  snow_albedo_angle_factor             ! function of solar zenith angle (for albedo of snow on land)
    REAL(dp)  ::  sky_view_fract                       ! fraction of bare ground below canopy for albedo calculation
    REAL(dp)  ::  sky_view_fract_stem                  ! fraction added to snow covered canopy due to stem area
    REAL(dp)  ::  z_zenith                             ! zenith angle
    REAL(dp)  ::  AlbedoCanopySnow                     ! albedo of snow covered canopy
    LOGICAL  ::  l_rad(nidx)                           ! flag to indicate if it is day or night
    INTEGER  ::  i,itile

#if defined (__SX__) && defined (_OPENMP)
    INTEGER  ::  tid, nt

    IF (debug) THEN
       tid = omp_get_thread_num()
       nt = omp_get_num_threads()
       CALL message('update_albedo_snowage_temp', &
                    'OpenMP thread #'//int2string(tid)//' of '//int2string(nt)//' nidx: '//int2string(nidx))
    END IF
#endif

    snow_albedo_age_vis(:,:) = 0.0_dp
    snow_albedo_age_nir(:,:) = 0.0_dp
    snow_albedo_temp_vis(:,:) = 0.0_dp
    snow_albedo_temp_nir(:,:) = 0.0_dp
    snow_albedo_vis(:,:) = 0.0_dp
    snow_albedo_nir(:,:) = 0.0_dp

    ! Weight AlbedoCanopySnow
    AlbedoCanopySnow = albedo_params%AlbedoAgeWeight * AlbedoCanopySnow_age + &
                       (1._dp - albedo_params%AlbedoAgeWeight) * AlbedoCanopySnow_temp

    ! mask where albedo is to be calculated
    l_rad(:) = .FALSE.
    DO i = 1,nidx
       IF (radiation_net_vis(i) + radiation_net_nir(i) > minimum_rad) l_rad(i) = .TRUE.
    END DO

    ! calculate albedo of snow (snow age scheme)
    DO itile = 1,ntiles
       DO i = 1,nidx
          z_zenith = cos_zenith(kstart+i-1)
          IF (snow_fract(i,itile) > EPSILON(1._dp)) THEN
             snow_age_factor = snow_age(i) / (1._dp + snow_age(i))
             IF (z_zenith < 0.5_dp) THEN
                snow_albedo_angle_factor = ((1.0_dp + ZenithAngleFactor)/ &
                                           (1.0_dp + 2.0_dp * ZenithAngleFactor * z_zenith) - 1.0_dp) / &
                                           ZenithAngleFactor
             ELSE
                snow_albedo_angle_factor = 0._dp
             END IF
             snow_albedo_age_vis(i,itile) = AlbedoSnowVisMax * (1._dp - AlbedoSnowVisAge * snow_age_factor)
             snow_albedo_age_vis(i,itile) = snow_albedo_age_vis(i,itile) + AlbedoSnowAngle * snow_albedo_angle_factor * &
                                             (1._dp - snow_albedo_age_vis(i,itile))
             snow_albedo_age_nir(i,itile) = AlbedoSnowNirMax * (1._dp - AlbedoSnowNirAge * snow_age_factor)
             snow_albedo_age_nir(i,itile) = snow_albedo_age_nir(i,itile) + AlbedoSnowAngle * snow_albedo_angle_factor * &
                                             (1._dp - snow_albedo_age_nir(i,itile))
          ELSE
             snow_albedo_age_vis(i,itile) = 0._dp
             snow_albedo_age_nir(i,itile) = 0._dp
          END IF
       END DO
    END DO

    ! calculate albedo of snow (ECHAM5 scheme)
    DO itile = 1,ntiles
       DO i = 1,nidx
          IF (surface_temperature(i,itile) >= tmelt) THEN
             snow_albedo_temp_vis(i,itile) = lctlib%AlbedoSnowVisMin(cover_type(i,itile))
             snow_albedo_temp_nir(i,itile) = lctlib%AlbedoSnowNirMin(cover_type(i,itile))
          ELSE IF (surface_temperature(i,itile) < tmelt - TempAlbedoSnowMax) THEN
             snow_albedo_temp_vis(i,itile) = lctlib%AlbedoSnowVisMax(cover_type(i,itile))
             snow_albedo_temp_nir(i,itile) = lctlib%AlbedoSnowNirMax(cover_type(i,itile))
          ELSE
             snow_albedo_temp_vis(i,itile) = lctlib%AlbedoSnowVisMin(cover_type(i,itile)) +             &
                  (tmelt - surface_temperature(i,itile)) * (lctlib%AlbedoSnowVisMax(cover_type(i,itile)) - &
                  lctlib%AlbedoSnowVisMin(cover_type(i,itile))) / TempAlbedoSnowMax
             snow_albedo_temp_nir(i,itile) = lctlib%AlbedoSnowNirMin(cover_type(i,itile)) +             &
                  (tmelt - surface_temperature(i,itile)) * (lctlib%AlbedoSnowNirMax(cover_type(i,itile)) - &
                  lctlib%AlbedoSnowNirMin(cover_type(i,itile))) / TempAlbedoSnowMax
          END IF
       END DO
    END DO

    ! weight snow_albedo_vis and snow_albedo_nir
    DO itile = 1,ntiles
       DO i = 1,nidx
          snow_albedo_vis(i,itile) = albedo_params%AlbedoAgeWeight * snow_albedo_age_vis(i,itile) + &
                                    (1._dp - albedo_params%AlbedoAgeWeight) * snow_albedo_temp_vis(i,itile)
          snow_albedo_nir(i,itile) = albedo_params%AlbedoAgeWeight * snow_albedo_age_nir(i,itile) + &
                                    (1._dp - albedo_params%AlbedoAgeWeight) * snow_albedo_temp_nir(i,itile)
       END DO
    END DO

    ! calculate new glacier albedo (visible and NIR) only at grid points with solar radiation or at model initialisation
    DO itile = 1,ntiles
       DO i = 1,nidx
          IF (l_rad(i) .AND. is_glacier(i,itile)) THEN
             IF (surface_temperature(i,itile) >= tmelt) THEN
                albedo_vis(i,itile) = AlbedoGlacierVisMin
                albedo_nir(i,itile) = AlbedoGlacierNirMin          
             ELSE IF (surface_temperature(i,itile) < tmelt - TempAlbedoGlacierMax) THEN
                albedo_vis(i,itile) = AlbedoGlacierVisMax
                albedo_nir(i,itile) = AlbedoGlacierNirMax
             ELSE
                albedo_vis(i,itile) = AlbedoGlacierVisMin + &
                     (tmelt - surface_temperature(i,itile)) * (AlbedoGlacierVisMax - AlbedoGlacierVisMin) / &
                     TempAlbedoGlacierMax
                albedo_nir(i,itile) = AlbedoGlacierNirMin + &
                     (tmelt - surface_temperature(i,itile)) * (AlbedoGlacierNirMax - AlbedoGlacierNirMin) / &
                     TempAlbedoGlacierMax
             END IF
          END IF
       END DO
    END DO

    ! calculate new albedo (visible and NIR) only at grid points with solar radiation or at model initialisation
    DO itile = 1,ntiles
       DO i = 1,nidx
          ! albedo of the canopy (vegetation)
          IF (albedo_options%useAlbedoCanopy) THEN
             background_albedo_canopy_vis = background_albedo_veg_vis(i,itile)
             background_albedo_canopy_nir = background_albedo_veg_nir(i,itile)
          ELSE
             background_albedo_canopy_vis = lctlib%AlbedoCanopyVIS(cover_type(i,itile))
             background_albedo_canopy_nir = lctlib%AlbedoCanopyNIR(cover_type(i,itile))
          END IF
          ! fraction for which albedo is computed from canopy
          sky_view_fract = 1.0_dp - (veg_ratio_max(i) * &
                          (1.0_dp - EXP(-SkyViewFactor * lai(i,itile))))
          ! fraction which is added to the (snow covered) canopy fraction due to stem area
          sky_view_fract_stem = sky_view_fract - (1.0_dp - (veg_ratio_max(i) * &
                                (1.0_dp - EXP(-SkyViewFactor * (lai(i,itile) + &
                                lctlib%StemArea(cover_type(i,itile)) *         &
                                (2.0_dp - cos_zenith(kstart+i-1)))))))

          IF (l_rad(i) .AND. .NOT. is_glacier(i,itile)) THEN
             background_albedo_vis = (1.0_dp - sky_view_fract) * background_albedo_canopy_vis + &
                                     sky_view_fract * background_albedo_soil_vis(i,itile)
             background_albedo_nir = (1.0_dp - sky_view_fract) * background_albedo_canopy_nir + &
                                     sky_view_fract * background_albedo_soil_nir(i,itile)
             ! albedo of forests = weighted mean of albedo of ground below canopy and albedo of canopy
             IF (is_forest(i,itile)) THEN
                albedo_vis(i,itile) =                                                                                 &
                     MAX((sky_view_fract - sky_view_fract_stem) * (snow_fract(i,itile) * snow_albedo_vis(i,itile) +   &
                     (1.0_dp - snow_fract(i,itile)) * background_albedo_soil_vis(i,itile)) +                          &
                     sky_view_fract_stem  * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                          &
                     (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_soil_vis(i,itile)) +                   &
                     (1.0_dp - sky_view_fract) * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                     &
                     (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_canopy_vis),                           &
                     background_albedo_vis)
                albedo_nir(i,itile) =                                                                                 &
                     MAX((sky_view_fract - sky_view_fract_stem) * (snow_fract(i,itile) * snow_albedo_nir(i,itile) +   &
                     (1.0_dp - snow_fract(i,itile)) * background_albedo_soil_nir(i,itile)) +                          &
                     sky_view_fract_stem * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                           &
                     (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_soil_nir(i,itile)) +                   &
                     (1.0_dp - sky_view_fract) * (canopy_snow_fract(i,itile) * AlbedoCanopySnow +                     &
                     (1.0_dp - canopy_snow_fract(i,itile)) * background_albedo_canopy_nir),                           &
                     background_albedo_nir)
             ELSE
                albedo_vis(i,itile) =                                              &
                     MAX(snow_fract(i,itile) * snow_albedo_vis(i,itile) +             &
                     (1.0_dp - snow_fract(i,itile)) * background_albedo_vis,          &
                     background_albedo_vis)
                albedo_nir(i,itile) =                                              &
                     MAX(snow_fract(i,itile) * snow_albedo_nir(i,itile) +             &
                     (1.0_dp - snow_fract(i,itile)) * background_albedo_nir,          &
                     background_albedo_nir)
             END IF
          END IF
       END DO
    END DO

  END SUBROUTINE update_albedo_snowage_temp

  !--------------------------------------------------------------------------------------------------
  ! Calculation of the diagnostic albedo of the whole solar range from albedo_vis and albedo_nir
  ! of the previous time step.
  ! The routine corresponds to update_albedo_snowage_temp for (l_trigrad=true).
  !--------------------------------------------------------------------------------------------------
  SUBROUTINE update_albedo_diag(nidx, ntiles, radiation_net_vis, radiation_net_nir, &
       albedo_vis, albedo_nir, albedo)

    INTEGER,  INTENT(in)    ::  nidx                         ! domain
    INTEGER,  INTENT(in)    ::  ntiles                       ! number of tiles
    REAL(dp), INTENT(in)    ::  radiation_net_vis(nidx)      ! net solar radiation in the visible range [W/m2]
    REAL(dp), INTENT(in)    ::  radiation_net_nir(nidx)      ! net solar radiation in the NIR range [W/m2]
    REAL(dp), INTENT(in)    ::  albedo_vis(nidx,ntiles)      ! albedo of the visible range
    REAL(dp), INTENT(in)    ::  albedo_nir(nidx,ntiles)      ! albedo of the NIR range
    REAL(dp), INTENT(inout) ::  albedo(nidx,ntiles)          ! albedo of the whole solar range

    ! Local variables
    REAL(dp)  ::  fraction_down_vis(nidx,ntiles)       ! fraction of solar downward radiation in the visible range
    LOGICAL   ::  l_rad(nidx)                          ! flag to indicate if it is day or night
    INTEGER   ::  i,itile


    ! mask where albedo is to be calculated
    l_rad(:) = .FALSE.
    DO i = 1,nidx
       IF (radiation_net_vis(i) + radiation_net_nir(i) > minimum_rad) l_rad(i) = .TRUE.
    END DO

    ! diagnose albedo (whole solar spectral range) at grid points with solar radiation
    DO itile = 1,ntiles
       DO i = 1,nidx
          IF (l_rad(i)) THEN
             fraction_down_vis(i,itile) = (radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile))) /    &
                                            ((radiation_net_vis(i) / (1.0_dp - albedo_vis(i,itile)))   &
                                           + (radiation_net_nir(i) / (1.0_dp - albedo_nir(i,itile))))
             albedo(i,itile) = fraction_down_vis(i,itile) * albedo_vis(i,itile) +  &
                               (1.0_dp - fraction_down_vis(i,itile)) * albedo_nir(i,itile)
          END IF
       END DO
    END DO

  END SUBROUTINE update_albedo_diag
!----------------------------------------------------------------------------------------------------------------

  SUBROUTINE land_surface_diagnostics(surface)

    USE mo_utils,        ONLY: average_tiles

    TYPE(land_surface_type), INTENT(in) :: surface

    !Local variables
    LOGICAL  :: mask(nidx,surface%ntiles)
    REAL(dp)     :: fract(nidx,surface%ntiles)
    INTEGER  :: ntiles
    INTEGER  :: kidx0, kidx1

    ntiles = surface%ntiles
    kidx0   = kstart
    kidx1   = kend

    ! Compute grid box averages
    mask  = surface%is_present(kidx0:kidx1,1:ntiles)
    fract = surface%cover_fract(kidx0:kidx1,1:ntiles)
    CALL average_tiles(surface%albedo(kidx0:kidx1,1:ntiles), mask, fract, land_surface_diag%albedo(kidx0:kidx1))

  END SUBROUTINE land_surface_diagnostics

!------------------------------------------------------------------------------
SUBROUTINE scale_cover_fract (ndim, ntiles, is_present, is_glacier, is_naturalveg, cover_fract)
!
! !DESCRIPTION:
!
! Rescaling of cover fractions to assure that
!  - the sum of cover fractions is one
!  - all non-glacier grid cells have at least a minimum vegetated fraction
!  - all glacier grid cells have a cover fraction of one on the glacier tile
!    and zero on the other tiles
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: ndim               ! vector length
    INTEGER,  INTENT(in)  :: ntiles             ! number of tiles
    LOGICAL,  INTENT(in)  :: is_present(:,:)    ! mask for land surface treated by jsbach 
    LOGICAL,  INTENT(in)  :: is_glacier(:,:)    ! logical glacier mask on tiles
    LOGICAL,  INTENT(in)  :: is_naturalveg(:,:) ! logical mask for natural vegetation
!
! !IN- and OUTPUT PARAMETERS:
! 
    REAL(dp), INTENT(inout) :: cover_fract(:,:) ! vegetated fraction

! !LOCAL VARIABLES:
!
    INTEGER   :: i, iter
    INTEGER   :: niter             ! number of iterations needed
    INTEGER   :: nsparce(ndim)     ! number of PFTs with a vegetated fraction of less then fract_small
    REAL(dp)  :: sum_fract(ndim)   ! sum of all cover fractions
    REAL(dp)  :: excluded_fract(ndim)   ! sum of all cover fractions

!------------------------------------------------------------------------------

    ! Make sure, crop and pasture have a cover fraction of at least fract_small

    WHERE (is_present(:,:) .AND. .NOT. is_naturalveg(:,:) .AND. .NOT. is_glacier(:,:))
       cover_fract(:,:) = MAX(fract_small,cover_fract(:,:))
    END WHERE

    ! The more tiles we have, the more iterations are needed. For binary identical results
    ! with and without an extra glacier tile, the number of iterations must not change.
    niter = ntiles
    IF (ANY(ALL(.NOT. is_glacier(:,:),DIM=2) .AND. ANY(.NOT. is_present(:,:),DIM=2))) niter = ntiles -1

    DO iter = 1, niter

       sum_fract(:) = 0._dp
       excluded_fract(:) = 0._dp
       nsparce(:) = 0

       DO i = 1,ntiles
          WHERE (cover_fract(:,i) > fract_small .AND. is_naturalveg(:,i))
             sum_fract(:) = sum_fract(:) + cover_fract(:,i)
          ELSEWHERE (is_naturalveg(:,i))
             nsparce(:) = nsparce(:) + 1
          ELSEWHERE
             excluded_fract(:) = excluded_fract(:) + cover_fract(:,i)
          END WHERE
       END DO
       DO i = 1,ntiles
          WHERE (cover_fract(:,i) > fract_small .AND. is_naturalveg(:,i))
             cover_fract(:,i) = cover_fract(:,i) * (1._dp - excluded_fract(:) - REAL(nsparce,dp)*fract_small) / sum_fract(:)
          ELSEWHERE (is_naturalveg(:,i))
             cover_fract(:,i) = fract_small
          ELSEWHERE (is_present(:,i))
             cover_fract(:,i) = MAX(fract_small,cover_fract(:,i))
          ELSEWHERE
             cover_fract(:,i) = 0._dp
          END WHERE
       END DO

    END DO

    IF (ANY(cover_fract(:,:) < fract_small .AND. is_present(:,:) .AND. .NOT. is_glacier(:,:))) THEN
       WRITE(message_text,*) 'cover_fract still smaller ', fract_small, ' after ', ntiles, ' iterations:', &
            MINVAL(MERGE(cover_fract(:,:), 1._dp, cover_fract > 0._dp .AND. .NOT. is_glacier(:,:)))
       CALL message ('scale_cover_fract', message_text)

       ! sometimes cover_fract remains slightly smaller than fract_small for numerical reasons
       WHERE (cover_fract(:,:) < fract_small .AND. cover_fract(:,:) >= fract_small - REAL(ntiles,dp)*EPSILON(1._dp) &
            .AND. is_present(:,:) .AND. .NOT. is_glacier(:,:))
          cover_fract(:,:) = fract_small
       END WHERE
       IF (ANY(cover_fract(:,:) < fract_small .AND. is_present(:,:) .AND. .NOT. is_glacier(:,:))) THEN
          WRITE(message_text,*) 'cover_fract still smaller ', fract_small, ' after ', ntiles, ' iterations:', &
               MINVAL(MERGE(cover_fract(:,:), 1._dp, cover_fract > 0._dp .AND. .NOT. is_glacier(:,:)))
          CALL finish ('scale_cover_fract', message_text)
       END IF
    END IF

    DO i = 1, ntiles
       WHERE (ANY(is_glacier(:,:), DIM=2) .AND. is_glacier(:,i))
          cover_fract(:,i) = 1._dp
       ELSEWHERE (ANY(is_glacier(:,:), DIM=2) .AND. .NOT. is_glacier(:,i))
          cover_fract(:,i) = 0._dp
       END WHERE
    END DO

  END SUBROUTINE scale_cover_fract

END MODULE mo_land_surface
