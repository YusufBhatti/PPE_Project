!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
Module mo_jsbach_grid

  ! 
  ! Description: 
  !   Definitions describing JSBACH grid, i.e. the grid used by the topography,
  !   vegetation and soil files given to JSBACH. The grid over which the LSS is
  !   actually run can be a subset of this grid, but it must be consistent with
  !   the larger grid (resolution, Gaussian or not, etc.)
  ! 
  ! Current Code Owner: jsbach_admin
  ! 
  ! History: 
  !  
  ! Version   Date        Comment 
  ! -------   ----        ------- 
  ! 0.1       2001/06/28  Original code. Reiner Schnur
  !
  ! CHR 28.02.03: in init_grid(): (1) added reading of Markos grid file 
  !                              (2) added deallocation of memory
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  ! Modules used: 
  ! 
  USE mo_jsbach, ONLY :  debug
  USE mo_mpi, ONLY: p_parallel, p_parallel_io, p_bcast, p_send, p_send, p_recv, p_io, p_pe, p_nprocs
  Use mo_kind,   Only : dp
  USE mo_exception,     ONLY : finish, message, message_text
  USE mo_util_string, ONLY: int2string
  USE mo_io_units, ONLY : nout
  USE mo_netCDF, ONLY: file_info, NF_NOERR, nf_inq_varid, nf_get_var_double
  USE mo_io, ONLY: IO_READ, IO_open, IO_close
  USE mo_linked_list, ONLY: t_stream, SURFACE
  USE mo_decomposition, ONLY: global_decomposition, local_decomposition

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: grid_type, domain_type
  PUBLIC :: init_grid, init_domain, update_comm_to_echam5mods
  PUBLIC :: kstart, kend, nidx, kindex

  ! Global (i.e. public) Declarations: 
  ! Global Type Definitions: 
  ! Grid for local decomposition (i.e. specific to each PE) 
  TYPE domain_type
     INTEGER           :: nlon           !! Number of longitudes
     INTEGER           :: nlat           !! Number of latitudes
     INTEGER           :: nland          !! Number of land grid cells in grid
     LOGICAL           :: lreg           !! Regular grid (nlon,nlat) or irregular (ndim,nblocks)?
     INTEGER           :: ndim           !! Maximum number of grid boxes in each block (nproma in ECHAM5)
     INTEGER           :: ndimz          !! Actual number of grid boxes in last block (npromz in ECHAM5)
     INTEGER           :: nblocks        !! Number of blocks
     LOGICAL, POINTER  :: mask(:,:)      !! Logical 2-d land mask (ndim,nblocks)
                                         !! Note: if lreg=.TRUE. then ndim=nlon and nblocks=nlat
     REAL(dp), POINTER :: mask_real(:,:) !! Real representation of mask (for netcdf output)
     INTEGER           :: kblock         !! Index of block on this PE that is processed for current call to LSS
                                         !! Note: if kblock=nblocks this is the last call to the LSS for the
                                         !!       current PE and time step
     LOGICAL           :: LastBlock = .FALSE. !! TRUE iff this is the last call to LSS for current PE and timestep
     REAL(dp), POINTER :: lon(:),lat(:)  !! Longitudes/latitudes of PE's land boxes
     REAL(dp), POINTER :: coslon(:), coslat(:), sinlon(:), sinlat(:)
     REAL(dp), POINTER :: elev(:)        !! Mean land grid cell elevation (m) (this is copied from land_surface_type,
                                         !! for convenience
  END TYPE domain_type

  TYPE grid_type
     INTEGER :: nlon, nlat                !! Number of global longitudes/latitudes
     INTEGER :: nland                     !! Number of global land points
     REAL(dp), POINTER :: lon(:), lat(:)  !! Global longitudes and latitudes
     LOGICAL,  POINTER :: mask(:,:)       !! Global land mask
     REAL(dp), POINTER :: mask_real(:,:)  !! Real representation of mask (for netcdf output)
     INTEGER,  POINTER :: kpoints(:)      !! global grid indeces of the land grid cells
     REAL(dp), POINTER :: landindex(:)    !! Indices of land grid cells
     !TYPE(file_info) :: IO_file
     !TYPE(t_stream)    :: IO_stream
  END TYPE grid_type

  TYPE(file_info) :: grid_file
  TYPE(t_stream), POINTER :: IO_grid            !! Grid stream (may be equal to IO_jsbach)

  LOGICAL, SAVE         :: grid_initialized   = .FALSE.
  LOGICAL, SAVE, PUBLIC :: domain_initialized = .FALSE.

  INTEGER               :: nidx           !! Number of land boxes for current call to LSS
  INTEGER               :: kstart         !! Start point of kindex
  INTEGER               :: kend           !! End point of kindex
  INTEGER, ALLOCATABLE  :: kindex(:)      !! Index of land boxes for current call to LSS into packed vector
                                          !!   of all land boxes of this PE's domain
!$OMP THREADPRIVATE(kstart,kend,nidx,kindex)

CONTAINS

  SUBROUTINE init_grid(grid, isStandalone, isRestart, IO_file_name)

    ! Initialize global grid (nlon, nlat, nland and mask) so that the decomposition can be defined

    USE mo_util_string, ONLY: int2string,real2string
    USE mo_netCDF,      ONLY: add_dim, NF_MAX_NAME, io_inq_dimid, io_inq_dimlen, io_inq_varid, io_get_var_double
    USE mo_temp,        ONLY: zreal2d, zint2d, zint1d
    USE mo_namelist,    ONLY: position_nml, POSITIONED
    USE mo_jsbach,      ONLY: nml_unit
    USE mo_gaussgrid,   ONLY: philon, philat

#ifdef STANDALONE
    ! To have correct values in mo_gaussgrid::inigau() when writing restart files
    USE mo_control, ONLY: echam_nlon => nlon, echam_nlat => ngl, echam_nhgl => nhgl
    USE mo_control, ONLY: nvclev
#endif
    TYPE(grid_type), INTENT(out) :: grid
    LOGICAL, INTENT(in) :: isStandalone
    LOGICAL, INTENT(in) :: isRestart
    CHARACTER(NF_MAX_NAME), INTENT(in) :: IO_file_name

    INTEGER                :: IO_file_id, IO_dim_id, IO_var_id

    INTEGER znlon, znlat                       !! number of longitudes/latitudes in global grid
    REAL(dp) :: lon_west, lon_east,      &     !! Edges of region over which the
                lat_south, lat_north           !! LSS should be run (only in standalone mode)
    REAL(dp), ALLOCATABLE :: zlon(:), zlat(:)  !! global longitudes, latitudes
    LOGICAL(dp), ALLOCATABLE :: l_zlon(:), l_zlat(:)
    INTEGER :: status, i, f_unit
    LOGICAL :: set_lon_zero                    !! whether to ignore lon from input file and set lon=0
                                               !! assuming that forcing is local time
    REAL(dp) :: time_offset,shift               !! adjust local time to timezone time (only needed if lon=0
                                               !! does not work (hours)

    INCLUDE 'jsbgrid_ctl.inc'


    IF (debug) CALL message('init_grid','BEGIN')

    IF (p_parallel_io) THEN

       IF (isStandalone) THEN

          !! Read namelist jsbgrid_ctl

          ! define default values
          lon_east = 360._dp
          lon_west = 0._dp
          lat_north = 90._dp
          lat_south = -90._dp
          set_lon_zero = .FALSE.
          time_offset = 0.

          f_unit = position_nml ('JSBGRID_CTL', nml_unit, status=status)
          SELECT CASE (status)
          CASE (POSITIONED)
             READ (f_unit, jsbgrid_ctl)
             CALL message('init_grid', 'Namelist JSBGRID_CTL: ')
             WRITE(nout, jsbgrid_ctl)
          END SELECT

          IF (.NOT. isRestart) THEN
             ! Get region over which to run LSS from namelist options
             ! Note that these define the boundaries of the region, not the centers of grid boxes
             
             WRITE(*,*) 'Running model for region:',lon_west, lon_east, lat_south, lat_north
          ENDIF

       END IF ! isStandalone

       grid_file       %opened = .FALSE.
       CALL io_open(IO_file_name, grid_file, IO_READ)
       IF (debug) WRITE(nout,*) 'Grid file successfully opened'
       IO_file_id = grid_file%file_id

       ! Get global grid dimensions from grid file
       IF (debug) WRITE(nout,*) 'Getting longitude and latitude dimensions'
       CALL io_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       CALL io_inq_dimlen (IO_file_id, IO_dim_id, znlon)
       CALL io_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       CALL io_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       IF (debug) THEN
          WRITE(nout,*) 'Got longitude and latitude dimensions'
          WRITE(nout,*) '  Dimensions in file:', znlon, znlat
       ENDIF
    END IF ! p_parallel_io


    IF (debug) CALL message('init_grid','Getting longitudes and latitudes')
    IF (isStandalone) THEN
       IF (p_parallel_io) THEN

          ! Get global lat/lon from grid file (in stand-alone runs, only)

          ALLOCATE(zlon(znlon))
          ALLOCATE(zlat(znlat))
          CALL io_inq_varid  (IO_file_id, 'lat', IO_var_id)
          CALL io_get_var_double (IO_file_id, IO_var_id, zlat)
          CALL io_inq_varid  (IO_file_id, 'lon', IO_var_id)
          CALL io_get_var_double (IO_file_id, IO_var_id, zlon)

          IF (.NOT. isRestart) THEN
             grid%nlon = COUNT(zlon >= lon_west .AND. zlon <= lon_east)
             grid%nlat = COUNT(zlat >= lat_south .AND. zlat <= lat_north)
          ELSE
             grid%nlon = znlon
             grid%nlat = znlat
          ENDIF
       END IF
       IF (p_parallel) THEN
          CALL p_bcast(grid%nlon, p_io)
          CALL p_bcast(grid%nlat, p_io)
       END IF
       ALLOCATE(grid%lon(grid%nlon), grid%lat(grid%nlat))
       IF (p_parallel_io) THEN
          ! In standalone mode it is possible to select a subregion of the global grid
          ! Note: In parallel mode, the decomposition has to be set up for the subregion.
          !       The option values defining the subregion are only used here to select
          !       the subregion from the global fields read from the files.
          IF (.NOT. isRestart) THEN
             ALLOCATE(l_zlon(znlon))
             l_zlon = zlon>=lon_west .AND. zlon<=lon_east
             grid%lon = PACK(zlon, MASK=l_zlon)
             DEALLOCATE(l_zlon)
             ALLOCATE(l_zlat(znlat))
             l_zlat = zlat>=lat_south .AND. zlat<=lat_north
             grid%lat = PACK(zlat, MASK=l_zlat)
             DEALLOCATE(l_zlat)
          ELSE
             grid%lon = zlon
             grid%lat = zlat
          END IF

          IF ( set_lon_zero .OR. time_offset /= 0._dp ) THEN
             CALL message('init_grid', 'Shifting longitudes.')
             IF( set_lon_zero ) THEN 
               grid%lon(:)=0._dp
               CALL message('init_grid', 'Setting all longitudes to zero degrees to force local time.')
               IF ( grid%nlon > 1 ) CALL finish('init_grid', &
                 'ERROR: setting ALL longitudes to ZERO although nlon = '// &
                 TRIM(int2string(grid%nlon))//'! It is not wise to continue...')
             ENDIF
             IF (time_offset /= 0._dp ) THEN
               shift = time_offset * 360._dp/24._dp ! convert hours to degrees
               grid%lon(:)=grid%lon(:) + shift 
               CALL message('init_grid', 'shifting longitudes by '//TRIM(real2string(shift))//' degrees.')
             ENDIF
          ENDIF

          DEALLOCATE(zlon, zlat)
       END IF
    ELSE
       IF (p_parallel_io) THEN
          grid%nlon = znlon
          grid%nlat = znlat
       END IF
    END IF
    WRITE(message_text,*) 'nlon: ', grid%nlon, ' nlat: ', grid%nlat
    CALL message('init_grid', message_text)

    IF (p_parallel) THEN
       CALL p_bcast(grid%nlon, p_io)
       CALL p_bcast(grid%nlat, p_io)
    ENDIF

    IF (.NOT. isStandalone) ALLOCATE(grid%lon(grid%nlon), grid%lat(grid%nlat))
    IF (p_parallel_io) THEN
       IF (.NOT. isStandalone) THEN

          ! in coupled runs longitudes and latitudes of the atmosphere model are used
          !   to assure restart reproducibility starting from the first run

          grid%lon = philon
          grid%lat = philat
       ENDIF
       lon_west  = MINVAL(grid%lon) ; lon_east  = MAXVAL(grid%lon)
       lat_south = MINVAL(grid%lat) ; lat_north = MAXVAL(grid%lat)
    ENDIF

    IF (p_parallel) THEN
       CALL p_bcast(lon_east, p_io)
       CALL p_bcast(lon_west, p_io)
       CALL p_bcast(lat_north, p_io)
       CALL p_bcast(lat_south, p_io)
       CALL p_bcast(grid%lon, p_io)
       CALL p_bcast(grid%lat, p_io)
    ENDIF
    IF (debug .AND. p_parallel_io) WRITE(nout,*) 'Got longitudes and latitudes'

#ifdef STANDALONE
    !! To have correct values of longitudes and latitudes in mo_gaussgrid::inigau() when writing 
    !! restart files in mo_io::write_streams() (problem arises because io_init(), where latitudes and longitudes
    !! are set, is not called in standalone mode)
    echam_nlat = grid%nlat
    echam_nlon = grid%nlon
    echam_nhgl = grid%nlat / 2
    nvclev = 0 ! needed in mo_grib. It has to be 0, otherwise vct is referenced without being allocated
#endif

    IF (debug .AND. p_parallel_io) WRITE(nout,*) 'Starting to get land-sea mask'
    ALLOCATE(grid%mask(grid%nlon,grid%nlat))
    IF (p_parallel_io) THEN
       ! Allocate temporary global 2-d fields
       ALLOCATE(zreal2d(grid%nlon,grid%nlat), zint2d(grid%nlon,grid%nlat))
       ! Get land-sea mask 
       ! Look for fraction of land 'slm', 'landseamask' or for 'gridnum'
       status = nf_inq_varid(IO_file_id, 'slm', IO_var_id)
       IF (status == NF_NOERR) THEN
          status = nf_get_var_double(IO_file_id, IO_var_id, zreal2d)
          IF (status /= NF_NOERR) CALL finish('init_grid', 'Error reading slm')
          grid%mask = zreal2d > 0.5_dp
       ELSE 
          status = nf_inq_varid(IO_file_id, 'landseamask', IO_var_id)
          IF (status == NF_NOERR) THEN
             status = nf_get_var_double(IO_file_id, IO_var_id, zreal2d)
             IF (status /= NF_NOERR) CALL finish('init_grid', 'Error reading landseamask')
             grid%mask = zreal2d > 0.5_dp
          ELSE
             status = nf_inq_varid(IO_file_id, 'gridnum', IO_var_id)
             IF (status /= NF_NOERR) &
                  CALL finish('init_grid', 'No land mask found in input file')
             status = nf_get_var_double(IO_file_id, IO_var_id, zint2d)
             IF (status /= NF_NOERR) CALL finish('init_grid', 'Error reading gridnum')
             grid%mask = zint2d > 0
          END IF
       END IF
       DEALLOCATE(zint2d, zreal2d)
    ENDIF
    IF (p_parallel) CALL p_bcast(grid%mask, p_io)

    grid%nland = COUNT(grid%mask)
    IF (debug .AND. p_parallel_io) WRITE(nout,*) 'Got land-sea mask'
    IF (p_parallel_io) WRITE(nout,*) 'nland:', grid%nland

    ! Global grid indices of land points
    ALLOCATE(grid%kpoints(grid%nland), zint2d(grid%nlon,grid%nlat), zint1d(grid%nlon*grid%nlat))
    zint1d = (/ (i, i=1,grid%nlon*grid%nlat) /)
    zint2d = RESHAPE(zint1d, (/grid%nlon,grid%nlat/))
    grid%kpoints = PACK(zint2d, MASK=grid%mask)
    DEALLOCATE(zint1d, zint2d)

    ! Index of land points
    ALLOCATE(grid%landindex(grid%nland))
    grid%landindex = (/ (i, i=1,grid%nland) /)

    IF (p_parallel_io) THEN
       CALL IO_close(grid_file)
    ENDIF

    ! Define dimensions
    CALL message('init_grid', 'Adding dimensions')
    CALL add_dim ("lon",     grid%nlon, "longitude", "degrees_east" )
    CALL add_dim ("lat",     grid%nlat,  "latitude" , "degrees_north")
    IF (grid%nland.eq.0) grid%nland = 1  ! Fix to allow for aqua planet experiments with no land.
    CALL add_dim ("landpoint", grid%nland)
    CALL add_dim ("surface",1, levtyp = 1, single = .TRUE., indx = SURFACE)

    grid_initialized = .TRUE.

    IF (debug) CALL message('init_grid','END')

  END SUBROUTINE init_grid

  SUBROUTINE init_domain(grid, domain, fileformat, fileztype, IO_stream, cbalone_flag)

    USE mo_jsbach_constants, ONLY: pi
    USE mo_tr_scatter, ONLY: scatter_field
    USE mo_temp,      ONLY: zreal2d, zreal2d_ptr, zint2d

    TYPE(grid_type),   INTENT(inout)  :: grid
    TYPE(domain_type), INTENT(inout)  :: domain
    INTEGER, OPTIONAL, INTENT(in)     :: fileformat  ! output file format
    INTEGER, OPTIONAL, INTENT(in)     :: fileztype   ! output file compression
    TYPE(t_stream), POINTER, OPTIONAL :: IO_stream
    LOGICAL, OPTIONAL, INTENT(in)     :: cbalone_flag

    ! Local variables
    TYPE(file_info)                :: IO_file
    INTEGER                        :: IO_file_id, IO_var_id
    REAL(dp), ALLOCATABLE, TARGET  :: lon2d(:,:), lat2d(:,:)
    INTEGER                        :: status, i, ibuf0, ibufn
    LOGICAL                        :: cbalone = .FALSE. !! flag for cbalance offline runs

    IF (PRESENT(cbalone_flag)) cbalone = cbalone_flag
!!$    IF (debug) CALL message('init_domain','BEGIN on PE '//int2string(p_pe))
    IF (debug) print*, 'init_domain ','BEGIN on PE ',p_pe

    IF (.NOT. grid_initialized) &
         CALL finish('init_domain', 'Call init_grid first')

    ! Copy domain grid information from local decomposition into grid structure
    domain%nlon    = local_decomposition%nglon
    domain%nlat    = local_decomposition%nglat
    domain%lreg    = local_decomposition%lreg
    domain%ndim    = local_decomposition%nproma
    domain%ndimz   = local_decomposition%npromz
    domain%nblocks = local_decomposition%ngpblks

    IF (p_parallel_io) THEN

       ! Open grid ini file or restart file
       IO_file%opened = .FALSE.
       CALL io_open(grid_file%file_name, IO_file, IO_READ)
       IO_file_id = IO_file%file_id

       ! Make sure that global grid dimensions from grid file are the same as those 
       ! defined in global decomposition
       IF (grid%nlon /= local_decomposition%nlon .OR. grid%nlat /= local_decomposition%nlat) THEN
          CALL message('  ','Global grid dimensions inconsistent with grid in ini/restart file')
          CALL message('  ','            nlon, nlat = '//int2string(grid%nlon)//'  '//int2string(grid%nlat))
          CALL message('  ','    from decomposition = '//int2string(local_decomposition%nlon)//'  '// &
                                                         int2string(local_decomposition%nlat))
          CALL finish ('JSBACH - init_domain:','exit')
       ENDIF

       ! Allocate temporary global 2-d fields
       ALLOCATE(zreal2d_ptr(grid%nlon,grid%nlat))
       ! Get land-sea mask (fraction of land)
       ! Look for fraction of land 'slm', or for 'gridnum'
       status = nf_inq_varid(IO_file_id, 'slm', IO_var_id)
       IF (status == NF_NOERR) THEN
          status = nf_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
          IF (status /= NF_NOERR) CALL finish('init_domain', 'Error reading slm')
       ELSE
          status = nf_inq_varid (IO_file_id, 'landseamask', IO_var_id)
          IF (status == NF_NOERR) THEN
             status = nf_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
             IF (status /= NF_NOERR) CALL finish('init_domain', 'Error reading landseamask')
          ELSE
             ALLOCATE(zint2d(grid%nlon,grid%nlat))
             status = nf_inq_varid (IO_file_id, 'gridnum', IO_var_id)
             IF (status /= NF_NOERR) &
                  CALL finish('init_domain', 'No land mask found in input file')
             status = nf_get_var_double(IO_file_id, IO_var_id, zint2d)
             IF (status /= NF_NOERR) CALL finish('init_domain', 'Error reading gridnum')
             WHERE(zint2d>0)
                zreal2d_ptr = 1._dp
             ELSEWHERE
                zreal2d_ptr = 0._dp
             END WHERE
             DEALLOCATE(zint2d)
          END IF
       END IF
    ENDIF

    IF (local_decomposition%npts > 0 .AND. grid%nland /= local_decomposition%npts) THEN
       CALL finish('init_domain', 'Number of global land points inconsistent with mask')
    ENDIF
    IF (local_decomposition%npts < 0) THEN
       local_decomposition%npts = grid%nland
       IF (debug) PRINT*,'init_domain - PE',p_pe,': Setting local_decomposition%npts to',grid%nland
       IF (p_parallel_io) THEN
          IF (debug) PRINT*,'init_domain - Setting global_decomposition%npts on all PEs to',grid%nland
          global_decomposition(:)%npts = grid%nland
       END IF
       IF (p_parallel) CALL p_bcast(global_decomposition(1:p_nprocs)%npts, p_io)
    END IF

    ! Temporary storage for local fields on domain
    ALLOCATE(zreal2d(domain%ndim, domain%nblocks))

    CALL scatter_field(zreal2d_ptr, zreal2d)
    ALLOCATE(domain%mask(domain%ndim,domain%nblocks))
    domain%mask = zreal2d > 0.5_dp
    domain%mask(domain%ndimz+1:domain%ndim,domain%nblocks) = .FALSE.

    IF (local_decomposition%ngpts > 0) THEN
       domain%nland = local_decomposition%ngpts
       IF (domain%nland /= COUNT(domain%mask)) &
            CALL finish('init_domain', 'Number of land points in decomposition inconsistent with mask')
    ELSE
       domain%nland = COUNT(domain%mask)
       local_decomposition%ngpts = domain%nland
       DO i=1,p_nprocs
          IF (global_decomposition(i)%pe == p_pe) THEN
             ibufn = domain%nland
          END IF
       END DO
       IF (p_pe == p_io) THEN
          global_decomposition(1)%ngpts = domain%nland
          DO i=2,p_nprocs
             CALL p_recv(ibuf0, global_decomposition(i)%pe, 1234)
             global_decomposition(i)%ngpts = ibuf0
          END DO
       ELSE
          CALL p_send(ibufn, p_io, 1234)
       END IF
    ENDIF
    IF (debug) THEN
       PRINT*,'init_domain (PE',p_pe,') domain%nland = ',domain%nland
       PRINT*,'init_domain (PE',p_pe,') local_decomposition%ngpts = ',local_decomposition%ngpts
       IF (p_parallel_io) PRINT*,'init_domain: global_decomposition%ngpts = ',global_decomposition%ngpts
    END IF

    IF (p_parallel_io) DEALLOCATE(zreal2d_ptr)
    !
    CALL update_comm_to_echam5mods(grid, domain)

    !------------------------------------------------------------------------------------------------------------
    ! Now that the domain dimensions and mask are defined the streams can be initialized and filled below
    IF (.NOT. cbalone) THEN
       CALL grid_init_memory(grid, domain, fileformat, fileztype, IO_stream)
       domain%mask_real = MERGE(1._dp,0._dp,domain%mask)
!!$    grid%mask_real = MERGE(1., 0., grid%mask)
    END IF
    !
    !------------------------------------------------------------------------------------------------------------
    ! Get local packed variables for each processor

    ! Fraction of land
    !domain%land_fract = PACK(zreal2d, MASK=domain%mask)

    ! Longitudes and Latitdes
    ALLOCATE(domain%lon(domain%nland), domain%lat(domain%nland))
    ALLOCATE(domain%coslon(domain%nland), domain%sinlon(domain%nland))
    ALLOCATE(domain%coslat(domain%nland), domain%sinlat(domain%nland))

    ALLOCATE(lon2d(grid%nlon,grid%nlat), lat2d(grid%nlon,grid%nlat))
    lon2d = SPREAD(grid%lon, DIM=2, NCOPIES=grid%nlat)
    lat2d = SPREAD(grid%lat, DIM=1, NCOPIES=grid%nlon)
    NULLIFY(zreal2d_ptr)
    zreal2d_ptr => lon2d
    CALL scatter_field(zreal2d_ptr, zreal2d)
    domain%lon = PACK(zreal2d, MASK=domain%mask)
    NULLIFY(zreal2d_ptr)
    zreal2d_ptr => lat2d
    CALL scatter_field(zreal2d_ptr, zreal2d)
    domain%lat = PACK(zreal2d, MASK=domain%mask)
    NULLIFY(zreal2d_ptr)
    DEALLOCATE(lon2d, lat2d)
    domain%coslon = COS(pi * domain%lon / 180._dp) ; domain%sinlon = SIN(pi * domain%lon / 180._dp)
    domain%coslat = COS(pi * domain%lat / 180._dp) ; domain%sinlat = SIN(pi * domain%lat / 180._dp)

    DEALLOCATE(zreal2d)
    IF (p_parallel_io .AND. .NOT. cbalone) THEN
       ! Free memory

       ! Close files
       CALL IO_close(IO_file)
    ENDIF

    IF (cbalone) THEN
       kstart = 1
       kend = domain%nland
       nidx = domain%nland
    END IF

    domain_initialized = .TRUE.

!!$    IF (debug) CALL message('init_domain','END on PE '//int2string(p_pe))
    IF (debug) print*,'init_domain ','END on PE ',p_pe

    RETURN
    
  END SUBROUTINE init_domain

  SUBROUTINE grid_init_memory(grid, domain, fileformat, fileztype, stream)

    ! Called from within init_grid after the decomposed grid has been computed. We need grid%nland from
    ! the decomposition to define the streams, on the other hand the streams have to be defined before 
    ! we can continue to read the initial fields

    USE mo_linked_list, ONLY : LAND, GAUSSIAN
    USE mo_memory_base, ONLY : new_stream, default_stream_setting, &
                               add =>add_stream_element
    USE mo_netCDF,      ONLY : max_dim_name
    USE mo_output,      ONLY : land_table

    TYPE(grid_type)                   :: grid        ! Global grid structure
    TYPE(domain_type)                 :: domain      ! PE domain structure
    INTEGER, INTENT(in)               :: fileformat  ! output file format
    INTEGER, INTENT(in)               :: fileztype   ! output file compression
    TYPE(t_stream), POINTER, OPTIONAL :: stream

    INTEGER                     :: dim2p(2), dim2(2)
    CHARACTER(LEN=max_dim_name) :: dim2n(2)

    IF (debug) CALL message('grid_init_memory','Begin')

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'grid', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE.)
       ENDIF
       IO_grid => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_grid, 'grid', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(IO_grid, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE.)
    ENDIF

    dim2p = (/ domain%nlon, domain%nlat /)
    dim2  = (/ grid%nlon, grid%nlat /)
    dim2n = (/ 'lon', 'lat' /)

    ! code numbers of this routine range from 1 to 9, using the GRIB land_table

    CALL add(IO_grid, 'landseamask', domain%mask_real, longname='Land Sea Mask', units='', &
            ldims=dim2p, gdims=dim2, dimnames=dim2n, repr=GAUSSIAN, code=1, lpost=.FALSE.)

    IF (debug) CALL message('grid_init_memory','END')

  END SUBROUTINE grid_init_memory

  SUBROUTINE update_comm_to_echam5mods(grid, domain)

    USE mo_jsbach_comm_to_echam5mods, ONLY: nlon, nlat, nland, lon, lat, kpoints, mask, &
                                     domain_nlon, domain_nlat, domain_nland, domain_mask

    TYPE(grid_type), INTENT(in) :: grid
    TYPE(domain_type), INTENT(in) :: domain

    nlon    = grid%nlon
    nlat    = grid%nlat
    nland   = grid%nland
    lon     => grid%lon
    lat     => grid%lat
    kpoints => grid%kpoints
    mask    => grid%mask

    domain_nlon  = domain%nlon
    domain_nlat  = domain%nlat
    domain_nland = domain%nland
    domain_mask  => domain%mask

  END SUBROUTINE update_comm_to_echam5mods



END MODULE mo_jsbach_grid


!Local Variables:
!mode: f90
!fill-column: 100
!End:
