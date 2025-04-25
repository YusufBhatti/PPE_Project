!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
Module mo_jsbach_veg
  !
  ! This module serves to
  ! (i) provide and initialize the infrastructure needed for the communication between the various submodels, and
  ! (ii) to call the initialization routines of the submodels.
  !
  ! I (CHR) would prefer to put the whole code of this module into mo_jsbach_interface, because it is actually
  ! an initialization of the interface. This would allow to make the fields (streams) initialized here PRIVATE instead of PUBLIC ---
  ! I am a real adict of data encapsulation!!
  !
  USE mo_netCDF, ONLY: FILE_INFO, NF_MAX_NAME
  USE mo_jsbach, ONLY: debug
  USE mo_jsbach_grid, ONLY: kstart, kend
  USE mo_linked_list, ONLY: t_stream
  USE mo_mpi, ONLY: p_parallel, p_parallel_io, p_bcast, p_io
  USE mo_kind, ONLY: dp 
  USE mo_exception, ONLY : message, finish

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART =======================================================================================================

  ! --- public subroutines

  PUBLIC :: vegetation_type, init_vegetation, veg_diagnostics

  ! --- public parameters and variables

  ! --- public fields used for communication between submodels in the jsbach-interface AND ONLY IN THE JSBACH-INTERFACE!!!!!

  TYPE vegetation_type
     INTEGER              :: ntiles
     INTEGER              :: nroot_zones
     REAL(dp),    POINTER, DIMENSION(:,:) :: & !! (nland, ntiles)
          canopy_conductance_bethy,          &       !! Canopy resistance
          canopy_conductance_limited,        &       !! Canopy resistance
          lai,                               &       !! Leaf Area Index
          lai_max,                           &       !! Projected annual maximum LAI for each phenology type ...
          veg_fract_correction,              &       !! Correction factor for cover fraction 1-exp(-LAI_max/2)
          snow_depth_canopy,                 &       !! Snow depth on canopy [m]
          snow_fract_canopy,                 &       !! Fraction of snow covered canopy
          veg_height                                 !! Vegetation Height (from lctlib file)
     REAL(dp),  POINTER, DIMENSION(:,:,:) :: & !! (nland, ntiles, 0:13)
          lai_clim                                   !! Climatological LAI
  END TYPE vegetation_type

  TYPE vegetation_diag_type
     REAL(dp), POINTER, DIMENSION(:) :: &
          canopy_conductance_bethy,     &
          canopy_conductance_limited,   &
          lai,                          &
          lai_max,                      &
          snow_depth_canopy,            &
          snow_fract_canopy
  END TYPE vegetation_diag_type


  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE 

  ! === Private declarations =======================================================================================================

  INTEGER, SAVE :: nroot_zones = -1 !! Number of root zones
 
  TYPE(t_stream),  POINTER, SAVE  :: IO_veg
  TYPE(t_stream),  POINTER, SAVE  :: IO_diag      !! Memory stream for diagnostic output
  TYPE(FILE_INFO), SAVE       :: veg_file
  TYPE(vegetation_diag_type), SAVE    :: veg_diag

Contains 

  !
  !=================================================================================================
  SUBROUTINE veg_init_io(IO_file_name)

    USE mo_netCDF,      ONLY: add_dim
    USE mo_io,          ONLY: IO_open, IO_READ, IO_close
    USE mo_linked_list, ONLY: ROOTZONES

    CHARACTER(NF_MAX_NAME), INTENT(in) :: IO_file_name
    INTEGER :: i    
    REAL(dp), ALLOCATABLE  :: root_values(:)


    IF (p_parallel_io) THEN
       ! Get number of root zones
       veg_file%opened = .FALSE.
       CALL IO_open(TRIM(IO_file_name), veg_file, IO_READ)
       CALL IO_close(veg_file)
       nroot_zones = 1 !! ONLY FOR TESTING!!!

    END IF

    IF (p_parallel) THEN
       CALL p_bcast(nroot_zones, p_io)
    ENDIF

    ALLOCATE (root_values(nroot_zones))
    DO i=1,nroot_zones 
       root_values(i) = REAL(i,dp)
    END DO
    CALL add_dim ("root_zone", nroot_zones, longname="root zone", units="", value=root_values(:), levtyp=72, indx=ROOTZONES)
    DEALLOCATE (root_values)

  END SUBROUTINE veg_init_io
  !
  !=================================================================================================
  SUBROUTINE veg_init_memory(g_nland, l_nland, ntiles, veg, fileformat, fileztype, diag_stream, stream)

    USE mo_jsbach,      ONLY : missing_value
    USE mo_linked_list, ONLY : LAND, TILES
    USE mo_memory_base, ONLY : new_stream, default_stream_setting, &
                               add =>add_stream_element
    USE mo_netCDF,      ONLY : max_dim_name
    USE mo_output,      ONLY : land_table

    INTEGER,               INTENT(in)    :: g_nland, l_nland, ntiles
    TYPE(vegetation_type), INTENT(inout) :: veg
    INTEGER,               INTENT(in)    :: fileformat            ! output file format
    INTEGER,               INTENT(in)    :: fileztype             ! output file compression
    TYPE(t_stream), POINTER              :: diag_stream
    TYPE(t_stream), POINTER, OPTIONAL    :: stream

    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim2p(2), dim2(2)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim2n(2)

    if (debug) CALL message('veg_init_memory','Entering ...')

    IF (ASSOCIATED(diag_stream)) THEN
       IO_diag => diag_stream
    ELSE
       CALL finish('veg_init_memory', 'Diagnostic stream not present')
    END IF

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'veg', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_veg => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_veg, 'veg', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(IO_veg, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim2p = (/ l_nland, ntiles /)
    dim2  = (/ g_nland, ntiles /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'tiles'

    ! code numbers of this routine range from 100 to 119, using the GRIB land_table

    CALL add(IO_veg, 'canopy_cond_bethy',  veg%canopy_conductance_bethy,       longname='Canopy Conductance BETHY', &
             units='m/s', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=105, lpost=.FALSE.)
    CALL add(IO_veg, 'canopy_cond_limited',veg%canopy_conductance_limited,     longname='Water Limited Canopy Conductance', &
             units='m/s', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=106, lpost=.FALSE.)
    CALL add(IO_veg, 'lai',                veg%lai,                            longname='Leaf Area Index', &
             units='',    ldims=dim2p, gdims=dim2, dimnames=dim2n, code=107, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_veg, 'lai_max',            veg%lai_max,                        longname='Maximum Annual Leaf Area Index', &
             units='',    ldims=dim2p, gdims=dim2, dimnames=dim2n, code=108, lpost=.FALSE.)
    CALL add(IO_veg, 'veg_fract_correction', veg%veg_fract_correction, &
             longname='Correction factor for cover fraction 1-exp(-LAI_max/2)', &
             units='',    ldims=dim2p, gdims=dim2, dimnames=dim2n, code=116, lrerun=.FALSE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_veg, 'snow_depth_canopy',  veg%snow_depth_canopy,              longname='Snow Depth on Canopy', &
             units='m',   ldims=dim2p, gdims=dim2, dimnames=dim2n, code=109, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_veg, 'snow_fract_canopy',  veg%snow_fract_canopy,              longname='Snow Fraction on Canopy', &
             units='',    ldims=dim2p, gdims=dim2, dimnames=dim2n, code=110, lmiss=.TRUE., missval=missing_value)
    CALL add(IO_veg, 'veg_height',         veg%veg_height,                     longname='Vegetation Height', &
             units='m',   ldims=dim2p, gdims=dim2, dimnames=dim2n, code=111, lpost=.FALSE.)

    CALL add(IO_diag, 'canopy_cond_bethy',       veg_diag%canopy_conductance_bethy,  longname='Canopy Conductance BETHY', &
         units='m/s',     ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=105, lpost=.FALSE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'canopy_cond_limited',     veg_diag%canopy_conductance_limited, longname='Water Limited Canopy Conductance', &
         units='m/s',     ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=106, lpost=.FALSE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'lai',                     veg_diag%lai,                       longname='Leaf Area Index', &
         units='',        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=107, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'lai_max',                 veg_diag%lai_max,                   longname='Maximum Annual Leaf Area Index', &
         units='',        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=108, lpost=.FALSE., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'snow_depth_canopy',       veg_diag%snow_depth_canopy,         longname='Snow Depth on Canopy', &
         units='m',       ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=109, &
         lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'snow_fract_canopy',       veg_diag%snow_fract_canopy,         longname='Snow Fraction on Canopy', &
         units='',        ldims=dim1p, gdims=dim1, dimnames=dim1n, laccu=.FALSE., code=110, &
         lmiss=.TRUE., missval=missing_value)

    if (debug) CALL message('veg_init_memory','    exit')

  END SUBROUTINE veg_init_memory
  !
  !=================================================================================================
  SUBROUTINE init_vegetation(grid, domain, vegetation, IO_file_name, &
       isRestart, fileformat, fileztype, IO_diag_stream, IO_stream)

    USE mo_jsbach_grid,   ONLY: grid_type, domain_type
    USE mo_tr_scatter,    ONLY: scatter_field
    USE mo_netcdf,        ONLY: nf_max_name, io_inq_varid, io_get_var_double
    USE mo_io,            ONLY: IO_open, IO_READ, IO_close
    Use mo_temp                          ! Provides temporary arrays

    TYPE(grid_type),         INTENT(in)    :: grid
    TYPE(domain_type),       INTENT(in)    :: domain
    TYPE(vegetation_type),   INTENT(inout) :: vegetation
    CHARACTER(nf_max_name),  INTENT(in)    :: IO_file_name
    LOGICAL,                 INTENT(in)    :: isRestart
    INTEGER,                 INTENT(in)    :: fileformat   ! output file format
    INTEGER,                 INTENT(in)    :: fileztype    ! output file compression
    TYPE(t_stream), POINTER           :: IO_diag_stream
    TYPE(t_stream), POINTER, OPTIONAL :: IO_stream

    TYPE(FILE_INFO) :: IO_file
    INTEGER :: ntiles
    INTEGER :: IO_file_id, IO_var_id
    INTEGER :: i
    integer :: status

    if(debug) call message("init_vegetation","Start initialization of mo_jsbach_veg")

    ntiles = vegetation%ntiles

    call veg_init_io(IO_file_name)  ! Sets nroot_zones
    vegetation%nroot_zones = nroot_zones

    ! --- generate vegetation stream

    call veg_init_memory(grid%nland, domain%nland, ntiles, vegetation, fileformat, fileztype, IO_diag_stream, stream=IO_stream)

    IF (.NOT. ASSOCIATED(IO_veg)) &
         CALL finish('init_vegetation','No memory stream for vegetation')

    IF (p_parallel_io) THEN

       ! Open ini file
       call message('init_vegetation','Reading vegetation fields from '//TRIM(veg_file%file_name))
       IO_file%opened = .FALSE.
       CALL IO_open(TRIM(veg_file%file_name), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

    ENDIF

    ! Temporary storage for local domain fields
    ALLOCATE(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('init_vegetation','Allocation failure (1)')

    !! Read in global LAI field for each calendar month
    ! Attention: The LAI read in here does not distinguish LCTs. It is used for all LCTs,
    !            limitted with the LCT specific lai_max.

    ALLOCATE(vegetation%lai_clim(domain%nland,ntiles,0:13))
    IF (p_parallel_io) THEN
       ALLOCATE(zreal3d(grid%nlon,grid%nlat,12),STAT=status)
       if(status .ne. 0) call finish('init_vegetation','Allocation failure')
       CALL IO_inq_varid(IO_file_id, 'lai_clim', IO_var_id)
       CALL io_get_var_double(IO_file_id, IO_var_id, zreal3d)
    ENDIF
    
    NULLIFY(zreal2d_ptr)
    DO i=1,12
       IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       CALL scatter_field(zreal2d_ptr, zreal2d)
       ! set clim. LAI to total LAI
       vegetation%lai_clim(:,:,i) = MAX(SPREAD(PACK(zreal2d, MASK=domain%mask),DIM=2, NCOPIES=ntiles),0.001_dp)
    END DO

    vegetation%lai_clim(:,:,0) = vegetation%lai_clim(:,:,12)
    vegetation%lai_clim(:,:,13) = vegetation%lai_clim(:,:,1)

    IF (p_parallel_io) DEALLOCATE(zreal3d)
    DEALLOCATE(zreal2d)

    ! If this is a restart run the model state variables are read in jsbach_init (Standalone) or the GCM
    ! by calling io_read_streams. We can therefore exit now:
    IF (isRestart) THEN
       IF (p_parallel_io) CALL IO_close(IO_file)
       RETURN
    ENDIF

    ! ----------------------------------------------------------------------------------------------------------------------
    ! If this is not a restart run we continue and get the vegetation model state from ini file
    !

    IF (p_parallel_io) CALL IO_close(IO_file)

    vegetation%snow_fract_canopy = 0._dp
    vegetation%snow_depth_canopy = 0._dp
    vegetation%lai = 0._dp

    veg_diag%canopy_conductance_bethy = 0._dp
    veg_diag%canopy_conductance_limited = 0._dp
    veg_diag%lai = 0._dp
    veg_diag%lai_max = 0._dp
    veg_diag%snow_depth_canopy = 0._dp
    veg_diag%snow_fract_canopy = 0._dp

    if(debug) call message("init_vegetation","Initialization of mo_jsbach_veg finished.")

  END SUBROUTINE init_vegetation

  !
  !=================================================================================================
  SUBROUTINE veg_diagnostics(surface, veg)

    USE mo_land_surface, ONLY: land_surface_type
    USE mo_utils,        ONLY: average_tiles

    TYPE(land_surface_type), INTENT(in) :: surface
    TYPE(vegetation_type),   INTENT(inout) :: veg

    !Local variables
    INTEGER  :: ntiles
    INTEGER  :: kidx0, kidx1

    ntiles  = veg%ntiles
    kidx0   = kstart
    kidx1   = kend

    CALL average_tiles(veg%canopy_conductance_bethy(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%canopy_conductance_bethy(kidx0:kidx1))
    CALL average_tiles(veg%canopy_conductance_limited(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%canopy_conductance_limited(kidx0:kidx1))
    CALL average_tiles(veg%lai(kidx0:kidx1,1:ntiles) * veg%veg_fract_correction(kidx0:kidx1,:), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%lai(kidx0:kidx1))
    veg_diag%lai(kidx0:kidx1) = veg_diag%lai(kidx0:kidx1) * surface%veg_ratio_max(kidx0:kidx1)
    CALL average_tiles(veg%lai_max(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%lai_max(kidx0:kidx1))
    CALL average_tiles(veg%snow_depth_canopy(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%snow_depth_canopy(kidx0:kidx1))
    CALL average_tiles(veg%snow_fract_canopy(kidx0:kidx1,1:ntiles), &
         surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
         veg_diag%snow_fract_canopy(kidx0:kidx1))

  END SUBROUTINE veg_diagnostics

End module mo_jsbach_veg

!Local Variables:
!mode: f90
!fill-column: 100
!End
