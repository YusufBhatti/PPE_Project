!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_hydrology

  !
  ! Authors:
  !
  ! S. Hagemann, MPI, October 1999, original source
  ! U. Schlese , MPI, January 2001, cleanup and introduction
  ! L. Kornblueh, MPI, April 2002, cleanup, parallelization,
  !                                and packed in one module
  ! K. Ketelsen, NEC, February 2003, optimization
  ! T. Jahns, DKRZ, November 2009, Speedup
  ! V. Gayler, MPI, Mai 2011, cleanup, correction for lakes,
  !                           new interpolation routines,
  !                           assure water conservation
  ! S. Hagemann, MPI, Feb. 2015,  routing put in subroutine
  !            Target box for routing can be read from input file (pre-requisite for ICON-HD)
  !            Correction of bug in HD offline version
  !            Remapping alternative implemented (no remapping)
  !

  USE mo_kind,          ONLY: dp
  USE mo_jsbach_grid,   ONLY: grid_type, domain_type
  USE mo_exception,     ONLY: message, message_text, finish
  USE mo_hd_highres_io, ONLY: hd_highres_open, hd_highres_write, hd_highres_close
  USE mo_gaussgrid,     ONLY: gridarea, philat, philon
  USE mo_io,            ONLY: io_open, io_close, io_read, io_write
  USE mo_netcdf,        ONLY: nf_global, nf_double, file_info, io_enddef, &
                              io_inq_dimid, io_inq_varid, io_inq_dimlen, &
                              io_def_dim, io_def_var, io_put_att_text, &
                              io_put_att_int, io_get_var_double, io_put_var_double, &
                              nf_get_att_int, nf_noerr
  USE mo_mpi,           ONLY: p_io, p_parallel_io, p_parallel, p_bcast
  USE mo_time_control,  ONLY: lstart, get_time_step, puthd, gethd, l_puthd, l_gethd, &
                              delta_time, ev_puthd, get_interval_seconds, &
                              get_time_step
  USE mo_tr_gather,     ONLY: gather_field
  USE mo_tr_scatter,    ONLY: scatter_field
  USE mo_array_utils,   ONLY: dec_monotonic_closest_midpoint, &
                              inc_monotonic_closest_midpoint

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_hydrology, cleanup_hydrology
  PUBLIC :: hydrology_model, hydrology_collect_land, hydrology_collect_lake, hydrology_restart

  INTEGER ::  nlon, ngl  ! number of Gaussian grid longitudes, latitudes

  TYPE cart_idx_2d
    INTEGER :: ilon, ilat
  END TYPE cart_idx_2d

  TYPE cart_coord_2d
    REAL(dp) :: longitude, latitude
    INTEGER  :: dir
  END TYPE cart_coord_2d

  TYPE cart_xidx_2d
    INTEGER :: ilon, ilat, extlen
    REAL(dp), ALLOCATABLE :: amod(:), akdiv(:)
  END TYPE cart_xidx_2d

  ! HD variables defined on land only (formerly defined in mo_memory_g3b)  
  REAL(dp), POINTER,  PUBLIC :: aros(:)
  REAL(dp), POINTER,  PUBLIC :: adrain(:)
  REAL(dp), POINTER,  PUBLIC :: apmecal(:)

  ! HD variables on the Gaussian grid (formerly defined in mo_memory_g3b)
  REAL(dp), POINTER,  PUBLIC :: disch(:,:)
  REAL(dp), POINTER,  PUBLIC :: awfre(:,:)

  ! HD grid parameters
  INTEGER,  PARAMETER, PUBLIC :: nl = 720    ! number of longitudes
  INTEGER,  PARAMETER, PUBLIC :: nb = 360    ! number of latitudes
  REAL(dp), PARAMETER :: fullcirc = 360.0_dp
  REAL(dp), PARAMETER :: hd_scal_lon = fullcirc/nl
  REAL(dp), PARAMETER :: hd_scal_lat = 0.5_dp*fullcirc/nb
  REAL(dp), PARAMETER :: florg = -180.0_dp   ! north-west corner longitude of gridbox(1,1)
  REAL(dp), PARAMETER :: fborg =   90.0_dp   ! north-west corner latitude of gridbox(1,1)
  REAL(dp), PARAMETER :: fscal =    0.5_dp   ! resolution

  ! corresponding coordinates on echam grid
  REAL(dp), PARAMETER :: oclorg = 0.0_dp, ocborg = 90.0_dp

  ! HD model time steps
  INTEGER         :: hd_calling_interval     ! calling interval in seconds
  INTEGER         :: riverflow_timestep      ! sub time step used for riverflow
  INTEGER, PUBLIC :: hd_steps_per_day        ! number of hd_model calls per day
  INTEGER         :: riverflow_steps_per_day ! number of riverflow time steps per day
  REAL(dp)        :: div_riverflow_timestep  ! 1/riverflow_timestep

  REAL(dp), ALLOCATABLE :: alf_k(:,:)    ! retention constant k, overflow
  REAL(dp), ALLOCATABLE :: alf_n(:,:)    ! number of reservoirs n, overflow
  REAL(dp), ALLOCATABLE :: arf_k(:,:)    ! retention constant k, riverflow
  REAL(dp), ALLOCATABLE :: arf_n(:,:)    ! number of reservoirs n, riverflow
  REAL(dp), ALLOCATABLE :: agf_k(:,:)    ! retention constant k, baseflow
  REAL(dp), ALLOCATABLE :: agf_n(:,:)    ! number of reservoirs n, baseflow
  INTEGER,  ALLOCATABLE :: fdir(:,:)     ! river direction
  REAL(dp), ALLOCATABLE :: hd_lsm(:,:)   ! hd model land mask
  INTEGER,  ALLOCATABLE :: filnew(:,:)   ! Longitude index of flow destination according to fdir
  INTEGER,  ALLOCATABLE :: fibnew(:,:)   ! Latitude index of flow destination according to fdir

  REAL(dp), ALLOCATABLE :: finfl(:,:)    ! inflow data
  REAL(dp), ALLOCATABLE :: fgmem(:,:,:)  ! intermediate linear baseflow reservoir
  REAL(dp), ALLOCATABLE :: frfmem(:,:,:) ! intermediate reservoirs, inflow cascade
  REAL(dp), ALLOCATABLE :: flfmem(:,:,:) ! intermediate reservoir, linear overflow

  REAL(dp), ALLOCATABLE :: hd_area(:)    ! grid cell area [m2]
  REAL(dp)              :: lon_hd(nl,nb) ! longitudes of the HD model grid
  REAL(dp)              :: lat_hd(nl,nb) ! latitudes of the HD model grid
  INTEGER, PARAMETER    :: nmemrf = 5    ! number of riverflow reservoir cascades
  LOGICAL               :: calculate_lonlat_hd = .TRUE.  ! calculate longitudes and latitudes
                                                         !   if not in hdpara file
  REAL(dp), PUBLIC      :: friv(nl,nb)            ! river flow
  REAL(dp), PUBLIC      :: water_to_ocean(nl,nb)  ! water going to the ocean directly

  ! Arrays for offline model if forcing on HD grid
  REAL(dp), ALLOCATABLE, PUBLIC :: runoff_s(:,:)         ! surface water runoff [m/s]
  REAL(dp), ALLOCATABLE, PUBLIC :: runoff_dr(:,:)        ! Drainage [m/s]

  TYPE(cart_idx_2d),  ALLOCATABLE :: oclook_cache(:,:)   ! closest ECHAM grid ocean cell
  TYPE(cart_idx_2d),  ALLOCATABLE :: intpol_mapping(:,:) ! no longer used (replaced by scrip remapping)
  TYPE(cart_xidx_2d), ALLOCATABLE :: arf_n_kas(:)        !
  TYPE(cart_xidx_2d), ALLOCATABLE :: alf_n_kas(:)        !
  TYPE(cart_xidx_2d), ALLOCATABLE :: agf_n_kas(:)        !

  ! Arrays on the Gaussian grid
  REAL(dp), ALLOCATABLE, TARGET :: gl_aros(:,:)    ! runoff
  REAL(dp), ALLOCATABLE, TARGET :: gl_adrain(:,:)  ! drainage
  REAL(dp), ALLOCATABLE, TARGET :: gl_apmecal(:,:) ! P-E of glaciers
  REAL(dp), ALLOCATABLE, TARGET :: gl_disch(:,:)   ! discharge
  REAL(dp), ALLOCATABLE, TARGET :: gl_awfre(:,:)   ! P-E over ocean and lakes
  REAL(dp), ALLOCATABLE         :: no_ocean_mask(:,:)    ! [1,0] land ocean mask: lakes count to land
  REAL(dp), ALLOCATABLE         :: land_fract(:,:)       ! land fraction
  REAL(dp), ALLOCATABLE         :: lake_fract(:,:)       ! lake fraction
  REAL(dp), ALLOCATABLE         :: ocean_fract(:,:)      ! ocean fraction

  ! Arrays for Scrip remapping
  INTEGER,  ALLOCATABLE :: src_address(:)   ! address indices of the source grid
  INTEGER,  ALLOCATABLE :: dst_address(:)   ! address indices of the destination grid
  REAL(dp), ALLOCATABLE :: remap_matrix(:,:)
  INTEGER               :: num_links        ! number of overlaps between source and destination grid
  INTEGER               :: num_weights      ! number of weights (1 for first order conservative remapping)
  INTEGER, PUBLIC       :: do_remapping = 1 ! remapping to the HD model grid (0: no remapping, 1: remapping)

  ! parameters of the hydrology namelist
  LOGICAL, PUBLIC :: ldebughd                     ! true for debugging
  LOGICAL, PUBLIC :: diag_water_budget            ! true to print water budget diagnostics
  INTEGER, PUBLIC :: nhd_diag                     ! index for hd diagnostics
  LOGICAL, PUBLIC :: lhd_highres = .FALSE.        ! true for additional output on the HD grid 
  LOGICAL, PUBLIC :: lbase = .TRUE.               ! baseflow ON or OFF
  LOGICAL, PUBLIC :: locean = .TRUE.              ! closure of water budget for ocean coupling
  REAL(dp)        :: fllog1, fblog1               ! grid cell for diagnostics (with nhd_diag=99)
  REAL(dp)        :: fllog2, fblog2               ! grid cell for diagnostics (with nhd_diag=99)
  LOGICAL, PUBLIC :: lhd_rout = .FALSE.           ! true: routing via index arrays, false: via direction arrays 

  INTEGER, SAVE   :: isolog_unit                  ! unit for outflow diagnostics
  LOGICAL, PUBLIC :: hydrology_offline = .FALSE.  ! true for HD offline simulations

CONTAINS

  SUBROUTINE init_hydrology (grid, domain, jsbach_offline, fileformat, fileztype, stream)

    USE mo_linked_list, ONLY: T_stream

    TYPE(grid_type),   INTENT(in)    :: grid
    TYPE(domain_type), INTENT(in)    :: domain
    LOGICAL, OPTIONAL, INTENT(in)    :: jsbach_offline
    INTEGER, OPTIONAL, INTENT(in)    :: fileformat       ! output file format
    INTEGER, OPTIONAL, INTENT(in)    :: fileztype        ! output file compression
    TYPE(t_stream), POINTER, OPTIONAL, INTENT(inout) :: stream

    REAL(dp), ALLOCATABLE       :: slm(:,:)           ! integer land sea mask
    REAL(dp), ALLOCATABLE       :: slf(:,:)           ! fractional land sea mask
    REAL(dp), ALLOCATABLE       :: lake(:,:)          ! fractional lake mask


    ! read HD namelist
    CALL config_hydrology

    ! initialize HD model memory
    IF (PRESENT(stream)) THEN
       CALL hd_init_memory(grid, domain, jsbach_offline, slm, slf, lake, fileformat, fileztype, stream)
    ELSE
       CALL hd_init_memory(grid, domain, jsbach_offline, slm, slf, lake)
    END IF

    CALL set_riverflow_timestep

    IF (p_parallel_io) THEN

     ! Read parameter fields and restart file for the HD Model
      CALL read_hydrology(slm, slf, lake)

      ! open file for 'isolog' timeseries
      IF (nhd_diag > 0) CALL hd_open_timeseries(nhd_diag)

    END IF

    CALL hydrology_slm_invariants(slm, slf, lake)
    DEALLOCATE(slm, slf, lake)

    IF (lhd_highres .AND. p_parallel_io) CALL hd_highres_open

  END SUBROUTINE init_hydrology

  SUBROUTINE config_hydrology

    !----------------------------------------------------------------------------------
    ! parameters of namelist HYDROLOGY_CTL
    !
    !          ldebughd    additional output for debugging  
    ! diag_water_budget    switch for additional water budget diagnostics  
    !       lhd_highres    switch for outflow diagnostic on HD model grid
    !             lbase    switch for baseflow calculations
    !            locean    closure of water budget for ocean coupling
    !          nhd_diag    region number for outflow diagnostic (formerly isolog)
    !                         0   none
    !                         1   Bothnian Bay/Sea
    !                         2   Torneaelven
    !                         4   St.Lawrence
    !                         5   Paraguay
    !                         6   Oder
    !                         7   Elbe
    !                         8   Oranje
    !                         9   Amudarya
    !                        10   Lena
    !                        99   user defined (fblog1, fllog1, fblog2, fllog2)
    !            fllog1    user defined grid cells for diagnostics (with nhd_diag=99)
    !            fblog1        fllog1, fblog1: longitude, latitude of grid cell 1
    !            fllog2        fllog2, fblog2: longitude, latitude of grid cell 2
    !            fblog2 
    !          lhd_rout    switch for routing scheme
    !                          true: via direction arrays;  false: via index arrays
    !----------------------------------------------------------------------------------
    USE mo_namelist,         ONLY: open_nml, position_nml, POSITIONED
    USE mo_io_units,         ONLY: nout
    USE mo_time_control,     ONLY: p_bcast_event, ev_gethd, ev_puthd, echam_ev_init
 
    ! local variables

    INTEGER                :: read_status, inml, iunit

    INCLUDE 'hydrology_ctl.inc'

    ! set default values of the namelist parmeters

    ldebughd = .FALSE.          ! additional output for debugging
    diag_water_budget = .FALSE. ! prints to diagnose the water budget
    lbase = .TRUE.              ! base flow calculation swiched on
    locean = .TRUE.             ! close water budget for ocean coupling
    nhd_diag = 0
    lhd_highres = .FALSE.
    fllog1 = 0.0_dp
    fblog1 = 0.0_dp
    fllog2 = 0.0_dp
    fblog2 = 0.0_dp
    lhd_rout = .FALSE.          ! routing following direction arrays

    ! read namelist hydrology_ctl

    IF (p_parallel_io) THEN
       inml = open_nml('namelist.jsbach')
       iunit = position_nml ('HYDROLOGY_CTL', inml, status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (iunit, hydrology_ctl)
          CALL message('config_hydrology', 'Namelist HYDROLOGY_CTL: ')
          WRITE(nout, hydrology_ctl)
       END SELECT
    END IF
    IF (p_parallel) THEN
       CALL p_bcast(ldebughd, p_io)
       CALL p_bcast(diag_water_budget, p_io)
       CALL p_bcast(locean, p_io)
       CALL p_bcast(lhd_highres, p_io)
       CALL p_bcast_event(puthd, p_io)
       CALL p_bcast_event(gethd, p_io)
       CALL p_bcast(lhd_rout, p_io)
    ENDIF

    CALL echam_ev_init(ev_gethd, gethd, 'couple get-from-hd', 'present')
    CALL echam_ev_init(ev_puthd, puthd, 'couple put-to-hd',   'next')

  END SUBROUTINE config_hydrology

  SUBROUTINE hd_init_memory(grid, domain, jsbach_offline, slm, slf, lake, fileformat, fileztype, stream)

    USE mo_output,      ONLY: land_table
    USE mo_netcdf,      ONLY: max_dim_name
    USE mo_linked_list, ONLY: GAUSSIAN, LAND, T_stream
    USE mo_memory_base, ONLY: new_stream, default_stream_setting, &
                              add =>add_stream_element

    TYPE(grid_type),         INTENT(in)  :: grid
    TYPE(domain_type),       INTENT(in)  :: domain
    LOGICAL,                 INTENT(in)  :: jsbach_offline     ! true for jsbach offline runs
    REAL(dp), ALLOCATABLE,   INTENT(out) :: slm(:,:)           ! integer land sea mask
    REAL(dp), ALLOCATABLE,   INTENT(out) :: slf(:,:)           ! fractional land sea mask
    REAL(dp), ALLOCATABLE,   INTENT(out) :: lake(:,:)          ! fractional lake mask
    INTEGER,  OPTIONAL,      INTENT(in)  :: fileformat         ! output file format
    INTEGER,  OPTIONAL,      INTENT(in)  :: fileztype          ! output file compression
    TYPE(t_stream), POINTER, OPTIONAL, INTENT(inout) :: stream

    TYPE(t_stream), POINTER     :: IO_hd
    INTEGER                     :: dim1p(1), dim1(1)
    CHARACTER(LEN=max_dim_name) :: dim1n(1)


    ! initialize memory of HD variables on the Gaussian grid (formerly in mo_memory_g3b)
    IF (PRESENT(stream)) THEN
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'hdjsbach', filetype=fileformat, ztype=fileztype)
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE.)
       ENDIF
       IO_hd => stream

       dim1p = (/ domain%nland /)
       dim1  = (/ grid%nland /)
       dim1n = (/ 'landpoint' /)
       CALL add (IO_hd, 'aros',    aros,    ldims=dim1p, gdims=dim1, dimnames=dim1n, lpost=.FALSE.)
       CALL add (IO_hd, 'adrain',  adrain,  ldims=dim1p, gdims=dim1, dimnames=dim1n, lpost=.FALSE.)
       CALL add (IO_hd, 'apmecal', apmecal, ldims=dim1p, gdims=dim1, dimnames=dim1n, lpost=.FALSE.)

       ! Variables defined on the Gaussian grid. The most intresting cells are water cells, however these cells 
       ! are masked out in jsbach offline runs. Thus no output is written in the offline case.
       CALL add (IO_hd, 'disch', disch, code=218, longname='river discharge', units='m/s', bits=24, laccu=.FALSE., &
                 repr=GAUSSIAN, lpost=(.NOT. jsbach_offline))
       CALL add (IO_hd, 'awfre', awfre, code=242, longname='P-E over ocean and lakes',  units='kg/m**2s', laccu=.FALSE., &
                 repr=GAUSSIAN, lpost=(.NOT. jsbach_offline))
    END IF

    ! initialize Gaussian grid dimensions
    nlon = grid%nlon
    ngl  = grid%nlat

    ! Initialize memory for the HD model variables on HD grid
    ! corresponds to offline routine 'hdini.f' by S. Hagemann

    IF (p_parallel_io) THEN

      ALLOCATE (alf_k(nl,nb))          ; alf_k(:,:)    = 0.0_dp
      ALLOCATE (alf_n(nl,nb))          ; alf_n(:,:)    = 0.0_dp
      ALLOCATE (arf_k(nl,nb))          ; arf_k(:,:)    = 0.0_dp
      ALLOCATE (arf_n(nl,nb))          ; arf_n(:,:)    = 0.0_dp
      ALLOCATE (agf_k(nl,nb))          ; agf_k(:,:)    = 0.0_dp
      ALLOCATE (agf_n(nl,nb))          ; agf_n(:,:)    = 0.0_dp
      ALLOCATE (fdir(nl,nb))           ; fdir(:,:)     = 0
      ALLOCATE (hd_lsm(nl,nb))         ; hd_lsm(:,:)   = 0.0_dp
      ALLOCATE (finfl(nl,nb))          ; finfl(:,:)    = 0.0_dp
      ALLOCATE (fgmem(nl,nb,1))        ; fgmem(:,:,:)  = 0.0_dp
      ALLOCATE (frfmem(nl,nb,nmemrf))  ; frfmem(:,:,:) = 0.0_dp
      ALLOCATE (flfmem(nl,nb,1))       ; flfmem(:,:,:) = 0.0_dp
      ALLOCATE (hd_area(nb))           ; hd_area(:)    = 0.0_dp
      ALLOCATE (oclook_cache(nl,nb))   ; oclook_cache  = cart_idx_2d(-1, -1)
      ALLOCATE (intpol_mapping(nl, nb)) ; intpol_mapping = cart_idx_2d(-1, -1)
      IF (do_remapping == 0) THEN
         ALLOCATE (runoff_s(nl,nb))     ; runoff_s(:,:) = 0.0_dp
         ALLOCATE (runoff_dr(nl,nb))    ; runoff_dr(:,:) = 0.0_dp
      ENDIF
      IF (lhd_rout) THEN
         ALLOCATE (filnew(nl,nb))      ; filnew(:,:)   = 0
         ALLOCATE (fibnew(nl,nb))      ; fibnew(:,:)   = 0
      ENDIF
    END IF

    ALLOCATE (gl_aros(nlon,ngl))    ; gl_aros(:,:)    = 0.0_dp
    ALLOCATE (gl_adrain(nlon,ngl))  ; gl_adrain(:,:)  = 0.0_dp
    ALLOCATE (gl_apmecal(nlon,ngl)) ; gl_apmecal(:,:) = 0.0_dp
    ALLOCATE (gl_disch(nlon,ngl))   ; gl_disch(:,:)   = 0.0_dp
    ALLOCATE (gl_awfre(nlon,ngl))   ; gl_awfre(:,:)   = 0.0_dp
    ALLOCATE (slf(nlon,ngl))        ; slf(:,:)        = 0.0_dp
    ALLOCATE (slm(nlon,ngl))        ; slm(:,:)        = 0.0_dp
    ALLOCATE (lake(nlon,ngl))       ; lake(:,:)       = 0.0_dp
    ALLOCATE (no_ocean_mask(nlon,ngl)) ; no_ocean_mask(:,:) = 0.0_dp
    ALLOCATE (land_fract(nlon,ngl))    ; land_fract(:,:)     = 0.0_dp
    ALLOCATE (lake_fract(nlon,ngl))    ; lake_fract(:,:)     = 0.0_dp
    ALLOCATE (ocean_fract(nlon,ngl))   ; ocean_fract(:,:)    = 0.0_dp

  END SUBROUTINE hd_init_memory

  SUBROUTINE read_hydrology (slm, slf, lake)
    !
    ! reads parameter fields and restart file for the HD model
    !
    !   hdrestart.nc: restart file with reservoir cascade arrays, ...
    !      hdpara.nc: parameter file with land sea mask, runoff directions, ...
    !      jsbach.nc: lake mask from jsbach initial file

    REAL(dp),       INTENT(out)    :: slm(:,:)       ! integer land sea mask
    REAL(dp),       INTENT(out)    :: slf(:,:)       ! fractional land sea mask
    REAL(dp),       INTENT(out)    :: lake(:,:)      ! fractional lake mask

    TYPE (FILE_INFO)  :: maskfile, parafile, hdfile

    INTEGER nvarid, fileid, i, status
    INTEGER hd_steps_per_day_restart    ! hd model timestep of the run the restart file originates from
    INTEGER riverflow_timestep_restart  ! riverflow timestep of the run the restart file originates from

    CHARACTER(len=80) :: fname
    CHARACTER(len= 7) :: varname
    REAL(dp) :: hd_dp_read(nl,nb), factor
    LOGICAL :: lex

    ! read masks from jsbach initial file
    !    these arrays cannot be provided from jsbach, as fractions of non-land
    !    cells are not available from the packed jsbach arrays

    fname = 'jsbach.nc'
    INQUIRE (file=fname, exist=lex)
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(fname),'>'
      CALL message('read_hydrology', message_text)
      CALL finish ('read_hydrology', 'run terminated.')
    ENDIF

    maskfile%opened = .FALSE.
    CALL IO_open (fname, maskfile, IO_READ)
    WRITE (message_text,*) 'Reading masks from file ', TRIM(fname)
    CALL message('read_hydrology', message_text)
    fileID = maskfile%file_id
    CALL IO_inq_varid (fileID, 'slm', nvarid, lex)
    IF (lex) CALL IO_get_var_double (fileID, nvarid, slm(:,:))
    CALL IO_inq_varid (fileID, 'slf', nvarid, lex)
    IF (lex) CALL IO_get_var_double (fileID, nvarid, slf(:,:))
    CALL IO_inq_varid (fileID, 'lake', nvarid, lex)
    IF (lex) CALL IO_get_var_double (fileID, nvarid, lake(:,:))
    CALL IO_close(maskfile)
   

    ! Read parameter: Land sea mask, RDF, ...

    fname = 'hdpara.nc'
    INQUIRE (file=fname, exist=lex)
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(fname),'>'
      CALL message('read_hydrology', message_text)
      CALL finish ('read_hydrology', 'run terminated.')
    ENDIF

    parafile%opened = .FALSE.
    CALL IO_open (fname, parafile, IO_READ)
    CALL message('', '')
    WRITE (message_text,*) 'Reading hdpara from file ', TRIM(fname)
    CALL message('read_hydrology', message_text)

    fileID = parafile%file_id

    CALL IO_inq_varid (fileID, 'lon', nvarid, lex)
    IF (lex) THEN
       CALL IO_get_var_double (fileID, nvarid, lon_hd(:,1))
       lon_hd(:,:) = SPREAD(lon_hd(:,1), NCOPIES=nb, DIM=2)
       CALL IO_inq_varid (fileID, 'lat', nvarid)
       CALL IO_get_var_double (fileID, nvarid, lat_hd(1,:))
       lat_hd(:,:) = SPREAD(lat_hd(1,:), NCOPIES=nl, DIM=1)
       calculate_lonlat_hd = .FALSE.
    END IF

    CALL IO_inq_varid (fileID, 'FLAG', nvarid)
    CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
    hd_lsm = REAL(NINT(hd_dp_read), dp)
    CALL IO_inq_varid (fileID, 'FDIR', nvarid)
    CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
    ! convert back to integer directions pointlessly stored as double
    fdir = MIN(NINT(hd_dp_read), 9)
    CALL IO_inq_varid (fileID, 'ALF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, alf_k)
    CALL IO_inq_varid (fileID, 'ALF_N', nvarid)
    CALL IO_get_var_double (fileID, nvarid, alf_n)
    CALL IO_inq_varid (fileID, 'ARF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, arf_k)
    CALL IO_inq_varid (fileID, 'ARF_N', nvarid)
    CALL IO_get_var_double (fileID, nvarid, arf_n)
    CALL IO_inq_varid (fileID, 'AGF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, agf_k)
    IF (lhd_rout) THEN
       CALL IO_inq_varid (fileID, 'FILNEW', nvarid)
       CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
       filnew = NINT(hd_dp_read)
       CALL IO_inq_varid (fileID, 'FIBNEW', nvarid)
       CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
       fibnew = NINT(hd_dp_read)
    ENDIF
    CALL IO_close(parafile)


    ! Read restart information: Reservoirs and inflow

    IF (lstart) THEN
      fname = 'hdstart.nc'
    ELSE
      fname = 'hdrestart.nc'
    END IF

    hdfile%opened = .FALSE.
    INQUIRE (file=fname, exist=lex)
    IF ( .NOT. lex ) THEN
      WRITE (message_text,*) 'Could not open file <', TRIM(fname), '>'
      CALL message('read_hydrology', message_text)
      CALL finish ('read_hydrology', 'run terminated.')
    ENDIF

    CALL IO_open (fname, hdfile, IO_READ)
    WRITE (message_text,*) 'Reading hdrestart from file ', TRIM(fname)
    CALL message('read_hydrology', message_text)
    CALL message('', '')

    fileID = hdfile%file_id

    CALL IO_inq_varid (fileID, 'FLFMEM', nvarid)
    CALL IO_get_var_double (fileID, nvarid, flfmem)

    varname = 'FRFMEM'
    DO i=1, nmemrf
      WRITE(varname(7:7), '(i1)') i
      CALL IO_inq_varid (fileID, varname, nvarid)
      CALL IO_get_var_double (fileID, nvarid, frfmem(:,:,I))
    ENDDO

    CALL IO_inq_varid (fileID, 'FGMEM', nvarid)
    CALL IO_get_var_double (fileID, nvarid, fgmem)
    fgmem(:,:,1) = fgmem(:,:,1) * hd_lsm(:,:) ! third array dimension is 1 

    CALL IO_inq_varid (fileID, 'FINFL', nvarid)
    CALL IO_get_var_double (fileID, nvarid, finfl)

    ! In principle the intermediate reservoirs have the unit [m3]. However,
    ! routine kasglob was designed for daily calls, and the reservoirs are
    ! treated as volume flows with unit [m3/day]. The state of the reservoir
    ! restart variables depend on the hd model time step, the inflow on the
    ! riverflow timestep. If the timesteps are changeing within a simulation
    ! the restart variables need to be adapted.

    status = NF_GET_ATT_INT (fileID, NF_GLOBAL,'hd_steps_per_day', hd_steps_per_day_restart)
    IF (status /= NF_NOERR)  hd_steps_per_day_restart = 1    ! old restart files
    IF (hd_steps_per_day_restart /= hd_steps_per_day) THEN
       factor = REAL(hd_steps_per_day,dp) / REAL(hd_steps_per_day_restart,dp)
       flfmem = flfmem * factor
       frfmem = frfmem * factor
       fgmem = fgmem * factor
    END IF

    status = NF_GET_ATT_INT (fileID, NF_GLOBAL,'riverflow_timestep', riverflow_timestep_restart)
    IF (status /= NF_NOERR) THEN
       ! handling of old restart files without riverflow_timestep attribute
       riverflow_timestep_restart = 4
       finfl = finfl / riverflow_timestep_restart
    END IF
    IF (riverflow_timestep_restart /= riverflow_timestep) THEN
        finfl = finfl * REAL(riverflow_timestep_restart,dp) * div_riverflow_timestep
    ENDIF

    CALL IO_close(hdfile)

  END SUBROUTINE read_hydrology

  SUBROUTINE hydrology_restart

    !
    ! **** Routine that writes the restart file for the HD model
    !
    ! ***** Version 1.0 - Oktober 1999
    !            Programmed and developed by Stefan Hagemann, MPI
    !
    !            Remark: Input data of Runoff and Drainage should have the
    !                       unit m/s.
    !
    ! ***** Version 1.1 - January 2001
    !       ECHAM5-Version
    !
    ! S.Legutke MPI M&D, Jan 2002, deallocate variables at end of
    !                              rerun cycle
    !
    ! ****** list of variables
    !
    !  frfmem(nl, nb, nmemrf) = Intermediate content of reservoir cascade
    !                           for the inflows per Gridbox (=5)
    !
    !  flfmem(nl, nb) = Intermediate content of linear reservoir for
    !                           Overland Flow
    !
    !   fgmem = Array of linear baseflow reservoir (Intermediate content)
    !           At Initialization it has the unit [m^3/s]
    !               (daily time step inherently implemented)
    !   finfl = Inflow data array for each gridbox for time step nstep
    !

    TYPE (FILE_INFO)  :: restartfile

    INTEGER :: nvarid, fileID, i, dims(2), xdimid, ydimid, xvarid, yvarid
    INTEGER :: istep

    CHARACTER(len=80) :: fname, string
    CHARACTER(len=7)  :: varname

    REAL(dp) :: lons(nl), lats(nb)

    IF (p_parallel_io) THEN

      fname = 'hdrestart.nc'

      istep = get_time_step()

      !    Open HD model restart file

      restartfile%opened = .FALSE.
      CALL IO_open (fname, restartfile, IO_WRITE)
      WRITE (message_text,*) 'Writing hdrestart to file ', TRIM(fname)
      CALL message('hydrology_restart', message_text)

      fileID = restartfile%file_id

      CALL IO_put_att_int (fileID, NF_GLOBAL, 'istep', istep)
      CALL IO_put_att_int (fileID, NF_GLOBAL, 'hd_steps_per_day', hd_steps_per_day)
      CALL IO_put_att_int (fileID, NF_GLOBAL, 'riverflow_timestep', riverflow_timestep)
      CALL IO_def_dim (fileID, 'lon', nl, xdimid)
      CALL IO_def_dim (fileID, 'lat', nb, ydimid)

      dims(1) = xdimid
      CALL IO_def_var (fileID, 'lon', NF_DOUBLE, 1, dims, xvarid)
      dims(1) = ydimid
      CALL IO_def_var (fileID, 'lat', NF_DOUBLE, 1, dims, yvarid)

      string = 'degrees_east'
      CALL IO_put_att_text (fileID, xvarid, 'units', string)
      string = 'Longitude'
      CALL IO_put_att_text (fileID, xvarid, 'long_name', string)
      string = 'degrees_north'
      CALL IO_put_att_text (fileID, yvarid, 'units', string)
      string = 'Latitude'
      CALL IO_put_att_text (fileID, yvarid, 'long_name', string)

      dims(1) = xdimid
      dims(2) = ydimid

      CALL IO_def_var (fileID, 'FLFMEM', NF_DOUBLE, 2, dims, nvarid)
      string = 'Linear overlandflow reservoir'
      CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
      string = 'm**3'
      CALL IO_put_att_text (fileID, nvarid, 'units', string)
      CALL IO_put_att_int (fileID, nvarid, 'code', 710)

      varname = 'FRFMEM'
      DO i=1, nmemrf
        WRITE(varname(7:7), '(i1)') i
        CALL IO_def_var (fileID, varname, NF_DOUBLE, 2, dims, nvarid)
        string = 'Inflow reservoir cascade'
        CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
        string = 'm**3'
        CALL IO_put_att_text (fileID, nvarid, 'units', string)
        CALL IO_put_att_int (fileID, nvarid, 'code', 710+i)
      ENDDO

      CALL IO_def_var (fileID, 'FGMEM', NF_DOUBLE, 2, dims, nvarid)
      string = 'Linear baseflow reservoir'
      CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
      string = 'm**3'
      CALL IO_put_att_text (fileID, nvarid, 'units', string)
      CALL IO_put_att_int (fileID, nvarid, 'code', 716)
      CALL IO_def_var (fileID, 'FINFL', NF_DOUBLE, 2, dims, nvarid)
      string = 'Inflow for each gridbox'
      CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
      string = 'm**3'
      CALL IO_put_att_text (fileID, nvarid, 'units', string)
      CALL IO_put_att_int (fileID, nvarid, 'code', 717)

      CALL IO_enddef (fileID)

      DO i = 1, nl
        lons(i) = -180.0_dp + 180.0_dp/nl + (i-1) * fullcirc /nl
      END DO
      DO i = 1, nb
        lats(i) = 90.0_dp - 90.0_dp/nb - (i-1)*180.0_dp/nb
      END DO

      CALL IO_put_var_double (fileID, xvarid, lons)
      CALL IO_put_var_double (fileID, yvarid, lats)

      CALL IO_inq_varid (fileID, 'FLFMEM', nvarid)
      CALL IO_put_var_double (fileID, nvarid, flfmem)

      DO i=1, nmemrf
        WRITE(varname(7:7), '(i1)') i
        CALL IO_inq_varid (fileID, varname, nvarid)
        CALL IO_put_var_double (fileID, nvarid, frfmem(:,:,I))
      ENDDO

      CALL IO_inq_varid (fileID, 'FGMEM', nvarid)
      CALL IO_put_var_double (fileID, nvarid, fgmem)

      CALL IO_inq_varid (fileID, 'FINFL', nvarid)
      CALL IO_put_var_double (fileID, nvarid, finfl)

      CALL IO_close (restartfile)

      ! close output streams

      IF (nhd_diag > 0) CALL hd_close_timeseries
      IF (lhd_highres) CALL hd_highres_close

    END IF

  END SUBROUTINE hydrology_restart

  SUBROUTINE cleanup_hydrology

    CALL cleanup_hd_slm_invariants

    IF (p_parallel_io) THEN

      DEALLOCATE  (alf_k)
      DEALLOCATE  (alf_n)
      DEALLOCATE  (arf_k)
      DEALLOCATE  (arf_n)
      DEALLOCATE  (agf_k)
      DEALLOCATE  (agf_n)
      DEALLOCATE  (fdir)
      DEALLOCATE  (hd_lsm)
      DEALLOCATE  (finfl)
      DEALLOCATE  (fgmem)
      DEALLOCATE  (frfmem)
      DEALLOCATE  (flfmem)
      DEALLOCATE  (hd_area)
      DEALLOCATE  (oclook_cache)
      DEALLOCATE  (intpol_mapping)
      IF (do_remapping == 0) THEN
         DEALLOCATE  (runoff_s)
         DEALLOCATE  (runoff_dr)
      ENDIF
      IF (lhd_rout) THEN
         DEALLOCATE  (filnew)
         DEALLOCATE  (fibnew)
      ENDIF
    END IF

    DEALLOCATE (no_ocean_mask)
    DEALLOCATE (land_fract)
    DEALLOCATE (lake_fract)
    DEALLOCATE (ocean_fract)

  END SUBROUTINE cleanup_hydrology

  SUBROUTINE hydrology_model(mask_domain, jsbach_offline, &
       aros_offline, adrain_offline, disch_offline, awfre_offline, &
       apmecal_offline)

    ! HD Model - Constants and Switches
    !
    ! **** Global/Regional Discharge Simulation as Subroutine for ECHAM5
    !
    !
    ! ***** Version 1.0 - November 1999
    !   Programmed and Developed by Stefan Hagemann, MPI
    !   Program code is based on Offline-Version of the HD model
    !   which is also regionally applicable. (regsim.f)
    !
    !   Anmerkung: Input data of Runoff und Drainage should have the unit m/s.
    !
    ! **** Remarks: Changes with regard to offline version
    !   Runoff array is now passed into routine instead of reading it
    !   in kasglob via echread. In echread, now named hdech and called
    !   before kasglob, only the transformation of the runoff array
    !   to the resolution of 0.5 degree is done if necessary.
    !   Parameter/Variables luinp, area are deleted from kasglob.
    !
    !   Since the input array to be transformed is only passed to echread
    !   but not read in ECHREAD itself, in ECHREAD only
    !   the transformation of the input array to 0.5 degree is done.
    !   Therefor, new calling parameter and routine names are set/given.
    !   old: CALL echread(luinp, ihead, tocode, istep, ique)
    !   new: CALL hdech(code_t42, tocode_0.5grad, ique)
    !
    !
    ! ***** River Direction File (RDF) format:
    !
    !                    7  8  9
    !                     \ | /
    !                      \|/
    !                    4--5--6
    !                      /|\
    !                     / | \
    !                    1  2  3
    !
    !       Remark: Direction 5 = Discharge Trap
    !               Direction -1 = Ocean Point
    !
    ! ****** List of variables
    !
    !  lbase = Baseflow ON or OFF
    ! locean = Closure of Water budget for ocean coupling.
    !
    ! isolog = Logfile output into Iso file (ASCII file , two columns)
    !      0 = no, 1 = Bothnian Bay/Sea, 2 = Torneaelven, 3 = Global...
    !      4 = St.Lawrence, 5 = Paraguay 6 = Odra
    !lhd_que = Log-Output switch  (.false. No Log-Output to STDOUT)
    !
    !  istep = Chosen  time step for Reading of Input
    !     nl = Number of Longitudes
    !     nb = Number of Latitudes
    !
    !     riverflow_timestep = internal sub-timestep for riverflow computation
    !
    ! **** Global Arrays:
    !
    !  finp = local input data array for time step istep
    ! fdata = local output data array for time step istep
    ! finfl = Inflow data array for each gridbox for time step istep
    !  fdir = River direction array
    !  hd_lsm = Land mask array
    ! alf_k = Array of retention constants k  - Overland flow [day]
    ! alf_n = Array of number of reservoirs n - Overland flow
    ! arf_k = Array of retention constants k  - Riverflow [day]
    ! arf_n = Array of number of reservoirs n - Riverflow
    ! agf_k = Array of retention constants k  - Baseflow [day]
    ! agf_n = Array of number of reservoirs n - Baseflow
    !
    !  frfmem(nl, nb, nmemrf) = Intermediate array of reservoir cascade for
    !                           the inflows per Gridbox (new: = nmemrf = 5)
    !
    !  flfmem(nl, nb, nmemlf) = Intermediate array of reservoir for
    !                           Surface Runoffs per Gridbox (new := nmemlf = 1)
    !
    ! fgmem = Array of linear baseflow reservoir (intermediate content)
    !         At initialization it has the unit [m^3/s]
    !  friv = Array of mean riverflow = Mean Inflow per Gridbox
    !
    ! hd_area(jb) = Array of gridbox arreas, Unit = [m^2]
    !
    !
    ! ***** Parameter and arrays of ECHAM grid
    !
    !  nlon = Longitudes of atmosphere grid
    !  ngl  = Latitudes of atmosphere grid
    !
    !  oclorg = Longitudinal origin of global atmosphere grid
    !  ocborg = Latitudinal origin of global atmosphere grid
    !  ocscal = resolution  = Latitudinal width of atmosphere gridbox in degree
    !
    !  aros   = runoff array
    !  adrain = drainage array
    !
    !  slm   = land sea mask on the Gaussian grid
    !  disch = discharge to the ocean (on the Gaussian grid)

    LOGICAL,  INTENT(in)              :: mask_domain(:,:)
    LOGICAL,  INTENT(in)              :: jsbach_offline      ! true for jsbach offline runs
    REAL(dp), INTENT(inout), OPTIONAL :: aros_offline(:,:)   ! INTENT(in)
    REAL(dp), INTENT(inout), OPTIONAL :: adrain_offline(:,:) ! INTENT(in)
    REAL(dp), INTENT(out),   OPTIONAL :: disch_offline(:,:)
    REAL(dp), INTENT(inout), OPTIONAL :: awfre_offline(:,:)
    REAL(dp), INTENT(in),    OPTIONAL :: apmecal_offline(:,:)

    REAL(dp) :: fdata(nl, nb)
    REAL(dp) :: finp(nl, nb)
    REAL(dp) :: outflow(nl,nb)
    REAL(dp) :: water_budget           ! budget of all water within the HD model
    REAL(dp) :: conservation_test      ! water budget difference at beginning and end of routine

    REAL(dp), ALLOCATABLE :: hlp2d(:,:)

    !  Parameter and switches

    INTEGER :: jl, jg, isub
    INTEGER :: dim(2)

    REAL(dp), POINTER :: gl(:,:)


    IF (hydrology_offline) THEN

       gl_aros = aros_offline
       gl_adrain = adrain_offline
       gl_apmecal = apmecal_offline
       gl_awfre = awfre_offline

    ELSE

       ! gather data from different processes

       dim(:) = SHAPE(mask_domain)
       ALLOCATE (hlp2d(dim(1),dim(2)))

       hlp2d = UNPACK(aros, FIELD=0._dp, MASK=mask_domain)
       gl => gl_aros
       CALL gather_field (gl, hlp2d)

       hlp2d = UNPACK(adrain, FIELD=0._dp, MASK=mask_domain)
       gl => gl_adrain
       CALL gather_field (gl, hlp2d)

       hlp2d = UNPACK(apmecal, FIELD=0._dp, MASK=mask_domain)
       gl => gl_apmecal
       CALL gather_field (gl, hlp2d)

       IF (jsbach_offline) THEN
          gl_awfre = 0._dp
          gl_disch = 0._dp
       ELSE
          gl => gl_awfre
          CALL gather_field (gl, awfre)
          gl => gl_disch
          CALL gather_field (gl, disch)
       END IF

    END IF

    ! from now on only work on IO node

    IF (p_parallel_io) THEN

      ! initializations

      water_to_ocean(:,:) = 0._dp   ! water that is not handled by the HD model
      outflow(:,:) = 0._dp

      ! initialization of water conservation test

      water_budget = &
             ! runoff: JSBACH land grid cells (slm)
             SUM(gl_aros*SPREAD(gridarea,1,nlon)*land_fract) &
             ! P-E over glaciers: JSBACH land grid cells (slm)
           + SUM(gl_apmecal*SPREAD(gridarea,1,nlon)*land_fract) &
             ! drainage: JSBACH land grid cells (slm)
           + SUM(gl_adrain*SPREAD(gridarea,1,nlon)*land_fract) &
             ! fresh water flux (awfre): defined on grid cells with at least a 
             ! small water fraction (SLF<1; collect). It comprises rain, snow
             ! and evaporation over water (evapw). Summed here is the part not
             ! falling on land, i.e. not treated by jsbach and contributing to
             ! runoff and drainage. The values of awfre over lakes will be
             ! added to the runoff. awfre is zero in atmoshere-only runs.
           + SUM(gl_awfre*SPREAD(gridarea,1,nlon)*(1._dp-land_fract)) &
             ! water in the overland flow reservoirs
           + SUM(flfmem) &
             ! water in the baseflow reservoir
           + SUM(fgmem) &
             ! water in the riverflow reservoirs
           + SUM(frfmem) &
             ! inflow data from the hd restart file
           + SUM(finfl)
      conservation_test=water_budget

      IF (ldebughd) THEN
         CALL message ('hydrology_model','------------------------------------------------------')
         CALL message ('hydrology_model','water budget at start of routine hydrology_model:' )
         WRITE (message_text,*) 'gl_aros: ', SUM(gl_aros*SPREAD(gridarea,1,nlon)*land_fract), 'm3'
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) 'gl_apmecal: ', SUM(gl_apmecal*SPREAD(gridarea,1,nlon)*land_fract), 'm3'
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) 'gl_drain: ', SUM(gl_adrain*SPREAD(gridarea,1,nlon)*land_fract), ' m3'
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) 'gl_awfre: ', SUM(gl_awfre*SPREAD(gridarea,1,nlon)*(1._dp-land_fract)), ' m3'
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) 'flfmem: ', SUM(flfmem), ' m3'
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) 'fgmem: ', SUM(fgmem), ' m3'
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) 'frfmem: ', SUM(frfmem), ' m3'
         CALL message ('hydrology_model',message_text)
         CALL message ('hydrology_model','------------------------------------------------------')
      END IF

      !! Uwe Mikolajewicz, 2009/10/27
      !! Put P-E on glaciers into surface runoff field
      !! disable de facto the old glacier calving by setting input field to 0!!!
      !! requires ability of the HD model to transport negative runoff,
      !! which is given in ECHAM5.
      
      gl_aros(:,:) = gl_aros(:,:) + gl_apmecal(:,:)
      gl_apmecal(:,:) = 0._dp


      ! P-E on lakes is added to the surface runoff to be transported by the HD model
      !     Note: gl_aros (from jsbach) is only calculated on land, not over lakes.

      gl_aros(:,:) = gl_aros(:,:)*(1._dp-lake_fract) + gl_awfre(:,:)*(lake_fract)
      gl_awfre(:,:) = gl_awfre(:,:)*(1._dp-lake_fract)

      ! ----------
      !  1 Runoff
      ! ----------

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 1: ', SUM(gl_aros*(1._dp-ocean_fract)*SPREAD(gridarea,1,nlon))
         CALL message ('hydrology_model',message_text)
      END IF

      ! Interpolation of the runoff from the ECHAM to the HD model grid

      IF (do_remapping /= 0) THEN
         CALL hd_remap(nlon, ngl, gl_aros*(1._dp-ocean_fract), nl, nb, locean, SPREAD(gridarea,1,nlon), &
                     SPREAD(hd_area,1,nl), finp)
      ELSE
         finp(:,:) = runoff_s(:,:)
      ENDIF

      !  Attention: Runoff in m/s --> Trafo with  AREA to m^3/s

      DO jl = 1,nl
         finp(jl,:) = finp(jl,:) * hd_area(:)
      END DO

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 1: ', SUM(finp)
         CALL message ('hydrology_model',message_text)
      END IF

      fdata(:,:) = 0.0_dp

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 2: ', SUM(flfmem) + SUM(finp)
         CALL message ('hydrology_model',message_text)
      END IF

      ! kasglob handles HD land points with positive reservoir numbers, only.
      ! Land points without outflow (fdir=5) have reservoir number zero.
      ! Besides, some HD ocean points have runoff values, due to the missmatch
      ! of the ECHAM and the HD land sea masks. This water will be given to
      ! the ocean directly to close the water balance.

      CALL kasglob(finp, fdata, alf_k, alf_n, hd_steps_per_day, flfmem, &
           alf_n_kas)
      WHERE (fdir == -1 .OR. fdir == 0 .OR. fdir == 5)
          water_to_ocean(:,:) = water_to_ocean(:,:) + finp(:,:)
      END WHERE     

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 2: ', &
              SUM(flfmem) + SUM(fdata) + SUM(water_to_ocean)
         CALL message ('hydrology_model',message_text)
      END IF

      ! ------------
      !  2 Drainage
      ! ------------

      IF (lbase) THEN

         IF (ldebughd) THEN
            WRITE (message_text,*) 'conservation test 3: ', &
                 SUM(gl_adrain*land_fract*SPREAD(gridarea,1,nlon)) + SUM(water_to_ocean)
            CALL message ('hydrology_model',message_text)
         END IF

         ! Interpolation of the drainage from the ECHAM to the HD model grid 

         IF (do_remapping /= 0) THEN
            CALL hd_remap(nlon, ngl, gl_adrain*land_fract, nl, nb, locean, SPREAD(gridarea,1,nlon), &
                       SPREAD(hd_area,1,nl), finp)
         ELSE
            finp(:,:) = runoff_dr(:,:)
         ENDIF

         ! *** Attention: Drainage in m/s --> Trafo with AREA to m^3/s

         DO jl = 1, nl
            finp(jl,:) = finp(jl,:) * hd_area(:)
         ENDDO

         IF (ldebughd) THEN
            WRITE (message_text,*) 'conservation test 3: ', &
                 SUM(finp) + SUM(water_to_ocean)
            CALL message ('hydrology_model',message_text)
            WRITE (message_text,*) 'conservation test 4: ', &
                 SUM(fgmem) + SUM(finp) + SUM(water_to_ocean)
            CALL message ('hydrology_model',message_text)
         END IF

         ! Some HD ocean points have drainage values, due to the missmatch of the ECHAM
         ! and the HD land sea masks. This water is ignored by kasglob. It is given to
         ! the ocean directly.

         CALL kasglob(finp, outflow, agf_k, agf_n, hd_steps_per_day, fgmem, agf_n_kas)
         WHERE (hd_lsm < 0.5)
            water_to_ocean(:,:) = water_to_ocean(:,:) + finp(:,:)
         END WHERE

         IF (ldebughd) THEN
            WRITE (message_text,*) 'conservation test 4: ', &
                 SUM(fgmem) + SUM(outflow) + SUM(water_to_ocean)
            CALL message ('hydrology_model',message_text)
         END IF

      END IF

      ! -------------
      !  3 Riverflow
      ! -------------

      ! initialization of riverflow

      friv(:,:) = 0.0_dp

      ! input for routing: fdata (overlandflow) + outflow (outflow from drainage)

      finp(:,:) = fdata(:,:) + outflow(:,:)

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 5: ', &
              SUM(finfl) + SUM(frfmem) + SUM(finp) + SUM(water_to_ocean)
         CALL message ('hydrology_model',message_text)
      END IF

      ! computation of riverflow in internal sub-time steps

      DO isub = 1, riverflow_timestep

         ! computing riverflow with input finfl from preceeding sub-time step

         ! kasglob handles HD land points with positive reservoir numbers, only.
         ! Land points without outflow (fdir=5) have reservoir number zero.
         ! finfl has non-zero values in regions with fdir=5.

         CALL kasglob(finfl, fdata, arf_k, arf_n, riverflow_steps_per_day, frfmem, &
              arf_n_kas)

         WHERE (fdir == -1 .OR. fdir == 0 .OR. fdir == 5)
            water_to_ocean(:,:) = water_to_ocean(:,:) + finfl(:,:)
         END WHERE

         ! calculate routing
         IF (lhd_rout) THEN
           call routing_via_index(finp, fdata)
         ELSE
           call routing(finp, fdata)
         ENDIF  

         WHERE (fdir == -1 .OR. fdir == 0 .OR. fdir == 5) ! non land point
            water_to_ocean(:,:)  = water_to_ocean(:,:) + finfl(:,:)
            finfl(:,:) = 0._dp
         END WHERE

         friv(:,:) = friv(:,:) + finfl(:,:)

      ENDDO   ! loop over sub-time steps

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 5: ', &
              SUM(frfmem) + SUM(finfl) + SUM(water_to_ocean)
         CALL message ('hydrology_model',message_text)
      END IF

      ! HD outflow diagnostics
      ! friv is defined on HD land (without internal drainage) and water_to_ocean is defined on hd 
      ! ocean (with internal drainage). The outflow diagnostics take into account both of these arrays.

      IF (nhd_diag /= 0 .OR. lhd_highres) CALL hydrology_diags(nhd_diag, water_to_ocean + friv)

      ! Conversion from the HD to the ECHAM Grid

      IF (locean) THEN
         IF (ldebughd) THEN
            WRITE (message_text,*) 'conservation test 6: ', SUM(water_to_ocean)
            CALL message ('hydrology_model',message_text)
         END IF
         CALL hydrology_to_ocean(nlon, ngl, water_to_ocean, fdir, gl_disch)
         IF (ldebughd) THEN
            WRITE (message_text,*) 'conservation test 6: ', SUM(gl_disch)
            CALL message ('hydrology_model',message_text)
         END IF
      END IF

      ! Convert discharge from m**3/s to m/s for ocean model

      DO jg = 1, ngl
        DO jl = 1, nlon
          gl_disch(jl, jg) = gl_disch(jl, jg)/gridarea(jg)
        END DO
      END DO

      water_budget = &
             ! Rain, snow and evaporation over water from ECHAM/JSBACH. P-E from
             ! lakes was given to the runoff.
             SUM(gl_awfre*SPREAD(gridarea,1,nlon)*(ocean_fract)) &
             ! discharge to the ocean
           + SUM(gl_disch*SPREAD(gridarea,1,nlon)) &
             ! water in the overland flow reservoirs
           + SUM(flfmem) &
             ! water in the baseflow reservoir
           + SUM(fgmem) &
             ! water in the riverflow reservoirs
           + SUM(frfmem) &
             ! inflow data
           + SUM(finfl)

      conservation_test = conservation_test - water_budget
      IF (diag_water_budget) THEN
         WRITE (message_text,*) 'Water budget change: ', conservation_test,' m3/s'
      END IF
      CALL message ('hydrology_model',message_text)
      IF (locean .AND. ABS(conservation_test/water_budget) > 1.E-7_dp) THEN
         WRITE (message_text,*) 'Water conservation problem: budget change: ', &
              conservation_test,' m3/s'
         CALL finish ('hydrology_model',message_text)
      END IF
      IF (diag_water_budget) THEN
         WRITE (message_text,*) '  global discharge: ', &
              SUM(gl_disch*SPREAD(gridarea,1,nlon)), ' m3/s'
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) '       precip-evap (on ocean): ', &
         SUM((gl_awfre)*SPREAD(gridarea,1,nlon)*(ocean_fract)), ' m3/s'
         CALL message ('hydrology_model',message_text)
      END IF

    END IF

    IF (.NOT. hydrology_offline) THEN
       gl => gl_aros
       CALL scatter_field (gl, hlp2d)
       aros = PACK(hlp2d, MASK=mask_domain)

       gl => gl_adrain
       CALL scatter_field (gl, hlp2d)
       adrain = PACK(hlp2d, MASK=mask_domain)

       IF (.NOT. jsbach_offline) THEN
          gl => gl_disch
          CALL scatter_field (gl, disch)

          gl => gl_awfre
          CALL scatter_field (gl, awfre)
       END IF
    ELSE
       disch_offline =  gl_disch
    END IF

  END SUBROUTINE hydrology_model

  SUBROUTINE hydrology_diags(isolog, hd_out)
    INTEGER,  INTENT(in) :: isolog
    REAL(dp), INTENT(in) :: hd_out(nl,nb)

    REAL(dp) :: f1, f2, fb, fl
    INTEGER :: jl, jb
    INTEGER :: istep

    !  Filling  F1, F2 at chosen coordiantes with
    !  OUTFLOW per Gridbox --> FINP, not HD_OUT or FINFL

    istep = get_time_step()

    IF (isolog == 1) THEN

      ! *** Log Output for Inflow into Gulf of Bothnia
      ! *** Since INFLOW ==> FINFL bzw. HD_OUT!
      ! *** Bothnian Bay: B=65.5 ,L=21.5 .... (glob = (404, 50))

      fl = 21.5_dp
      fb = 65.5_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 =  hd_out(jl,  jb  ) + hd_out(jl+1,jb  ) + hd_out(jl+2,jb  ) &
          + hd_out(jl+3,jb  ) + hd_out(jl+4,jb  ) + hd_out(jl+5,jb  ) &
          + hd_out(jl+6,jb  ) + hd_out(jl-1,jb+1) + hd_out(jl+5,jb+1) &
          + hd_out(jl-1,jb+2) + hd_out(jl+1,jb+2) + hd_out(jl+2,jb+2) &
          + hd_out(jl+3,jb+2) + hd_out(jl+4,jb+2) + hd_out(jl-2,jb+3) &
          + hd_out(jl-1,jb+4) + hd_out(jl,  jb+4)

      ! *** Bothnian Sea: B=63.5 ,L=19.0 ....

      f2 =  hd_out(jl-5,jb+4 ) + hd_out(jl-4,jb+4 ) + hd_out(jl-7,jb+5 ) &
          + hd_out(jl-8,jb+6 ) + hd_out(jl-2,jb+6 ) + hd_out(jl-8,jb+7 ) &
          + hd_out(jl-2,jb+7 ) + hd_out(jl-1,jb+7 ) + hd_out(jl-9,jb+8 ) &
          + hd_out(jl-1,jb+8 ) + hd_out(jl-8,jb+9 ) + hd_out(jl-1,jb+9 ) &
          + hd_out(jl-6,jb+10) + hd_out(jl-1,jb+10) + hd_out(jl,  jb+10) &
          + hd_out(jl,  jb+11) + hd_out(jl+1,jb+11) + hd_out(jl+2,jb+11)

    ELSE IF (isolog == 2 .OR. isolog == 3) THEN

      ! *** Torneaelven-Outflow = (22.0 E, 65.5 N), (22.5 E, 65.5 N)
      ! ***                       (23.5 E, 65.5 N)
      ! *** regional System:

      fl = 22.0_dp
      fb = 65.5_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 = REAL(istep,dp)
      f2 = hd_out(jl,jb) + hd_out(jl+1,jb) + hd_out(jl+3,jb)

    ELSE IF (isolog == 4) THEN

      ! *** St.Lawrence-Outflow = (-71.5 W, 47.0 N)
      ! *** regional System: Measurement station at (-75.5 W, 45.5 N)

      fl = -71.5_dp
      fb =  47.0_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 = hd_out(jl,jb)
      f2 = hd_out(jl-8,jb+3)

    ELSE IF (isolog == 5) THEN

      ! *** Paraguay-Outflow = (-59 W, -27 N)

      fl = -59.0_dp
      fb = -27.0_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 = REAL(istep,dp)
      f2 = hd_out(jl,jb)

    ELSE IF (isolog == 6) THEN

      ! *** Oder-Outflow = (14.0 E, 54.5 N), Hohensaaten-Finow (14 E, 53W)

      fl = 14.0_dp
      fb = 54.5_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 = hd_out(jl,jb)
      f2 = hd_out(jl+1,jb+2)+hd_out(jl+2,jb+2)

   ELSE IF (isolog == 7) THEN

      ! *** Elbe-Outflow = (8.5 E, 54.5 N), Neu-Darchau (10.5 E, 53.5 N)

      fl = 8.5_dp
      fb = 54.5_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 = hd_out(jl,jb)
      f2 = hd_out(jl+4,jb+2)

   ELSE IF (isolog == 8) THEN

      ! *** Oranje-Outflow = (-28.5 S, 16.0 E), Congo (-6.0 S, 12.0 E)

      fl = 16.0_dp
      fb = -28.5_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 = hd_out(jl,jb)

      fl = 12.0_dp
      fb = -6.0_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f2 = hd_out(jl,jb)

   ELSE IF (isolog == 9) THEN

      ! *** Amudarya-Outflow (47) = (43 N, 59 E), Syrdarya (49) = (46 N, 62 E)

      fl = 59.0_dp
      fb = 43.0_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 = hd_out(jl,jb)

      fl = 62.0_dp
      fb = 46.0_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f2 = hd_out(jl,jb)

   ELSE IF (isolog == 10) THEN

      ! *** Lena-Outflow (40) = (72 N, 127 E), Ob (46) = (67 N, 71.5 E)

      fl = 127.0_dp
      fb = 72.0_dp
      jl = NINT((fl-florg)/fscal + 1.+dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 = hd_out(jl,jb)

      fl = 71.5_dp
      fb = 67.0_dp
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f2 = hd_out(jl,jb)

   ELSE IF (isolog == 99) THEN

      ! *** user defined outflow coordinates (namelist hydrology_ctl)

      fl = fllog1
      fb = fblog1
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f1 = hd_out(jl,jb)

      fl = fllog2
      fb = fblog2
      jl = NINT((fl-florg)/fscal + 1._dp)
      jb = NINT(1._dp + (fborg-fb)/fscal)
      f2 = hd_out(jl,jb)

   ENDIF

   IF (isolog > 0) CALL hd_write_timeseries (f1, f2)

   IF (lhd_highres) CALL hd_highres_write (hd_out)

  END SUBROUTINE hydrology_diags

  SUBROUTINE hydrology_collect_land (kidx0, kidx1, nidx, pros_hd, pdrain_hd, palac)

    !  Collects runoff and drainage from jsbach (land grid cells only)
    !
    !  hydrology_collect is called from jsbach_inter_1d
    !

    INTEGER,  INTENT(in) :: kidx0, kidx1, nidx
    REAL(dp), INTENT(in) :: pros_hd(nidx),   & ! surface runoff [m]
                            pdrain_hd(nidx), & ! drainage [m]
                            palac(nidx)        ! P-E on glaciers [m]

    REAL(dp) ::  zrmean

    ! set accumulated runoff variables zero after HD/coupling time step

    IF (l_gethd) THEN
       aros(kidx0:kidx1)    = 0.0_dp
       adrain(kidx0:kidx1)  = 0.0_dp
       apmecal(kidx0:kidx1) = 0.0_dp
    END IF

    ! accumulate variables for HD-model

    aros(kidx0:kidx1)    = aros(kidx0:kidx1) + pros_hd(:)
    adrain(kidx0:kidx1)  = adrain(kidx0:kidx1) + pdrain_hd(:)
    apmecal(kidx0:kidx1) = apmecal(kidx0:kidx1) + palac(:)

    ! convert from [m/s] to [m]

    IF (l_puthd) THEN
       zrmean = get_interval_seconds(ev_puthd)
       IF (zrmean > 0.0_dp) zrmean = 1.0_dp/zrmean
       aros(kidx0:kidx1)    = aros(kidx0:kidx1)    * zrmean
       adrain(kidx0:kidx1)  = adrain(kidx0:kidx1)  * zrmean
       apmecal(kidx0:kidx1) = apmecal(kidx0:kidx1) * zrmean
    END IF

  END SUBROUTINE hydrology_collect_land

  SUBROUTINE hydrology_collect_lake (kdim, kblock, slf, ice_fract, &
       alake, rain_convective, &
       snow_convective, rain_largescale, snow_largescale, evapw)
    !  Collects P-E over lake and ocean
    !    this formerly happened in echam routine collect.f90
    !
    !  hydrology_collect_lake is called from echam routine physc
    !
    USE mo_jsbach_constants, ONLY: rhoh2o

    INTEGER,   INTENT(in)  :: kdim, kblock
    REAL(dp),  INTENT(in)  :: slf(kdim)             ! land fraction
    REAL(dp),  INTENT(in)  :: ice_fract(kdim)       ! lake/sea ice fraction
    REAL(dp),  INTENT(in)  :: alake(kdim)           ! lake fraction
    REAL(dp),  INTENT(in)  :: rain_convective(kdim) ! convective rain fall [kg/m2s]
    REAL(dp),  INTENT(in)  :: snow_convective(kdim) ! convective snow fall [kg/m2s]
    REAL(dp),  INTENT(in)  :: rain_largescale(kdim) ! large scale rain fall [kg/m2s]
    REAL(dp),  INTENT(in)  :: snow_largescale(kdim) ! large scale snow fall [kg/m2s]
    REAL(dp),  INTENT(in)  :: evapw(kdim)           ! evaporation over lake and ocean [kg/m2s]

    REAL(dp)   zrcouple

    ! accumulate P-E over lakes and ocean
    WHERE (slf(:) < 1._dp)
       awfre(1:kdim,kblock) = awfre(1:kdim,kblock) + (rain_convective(:)+rain_largescale(:)) * delta_time      &
                        + (snow_convective(:)+snow_largescale(:)+evapw(:)) * (1._dp-ice_fract(:)) * delta_time
    END WHERE

    ! add snow over frozen lakes (as otherwise the equivalent water is missing in the coupled atmosph/ocean system)
    WHERE (slf(:) < 1._dp .AND. alake(:) > 0.5_dp)
       awfre(1:kdim,kblock) = awfre(1:kdim,kblock) + (snow_convective(:) + snow_largescale(:)) * &
                              ice_fract(:) * delta_time
    END WHERE

    ! average accumulated values and convert flux from [kg/m2s] to [m/s]
    zrcouple = 1.0_dp/REAL(hd_calling_interval,dp)
    IF (l_puthd) THEN
       awfre(1:kdim,kblock) = awfre(1:kdim,kblock) * zrcouple / rhoh2o
    END IF

  END SUBROUTINE hydrology_collect_lake

  SUBROUTINE hydrology_echam (field_in, field_out, lhd_que)

    !*************************************************************************
    !
    ! **** This program interpolates data from Gaussian grids to a half
    !    degree grid
    !
    !  Programmierung und Entwicklung: Uwe Schulzweida (echamto30min)
    !  Modified to Subroutine by Stefan Hagemann -- September 1995
    !
    ! ***** Version 1.1 -- Dezember 1995
    !          Instead of longitude centred coordinate array trcode, now
    !          the 0.5 degree coordinate field field_out is given back to the calling
    !          routine which has a perfect boundary with the Northpole/dateline
    !          --> Origin has centre coordinate 89.75 N, -179.75 W
    !
    ! ***** Version 2.0 -- November 1999
    !   Since the input array to be transformed is only passed to ECHREAD
    !   but not read in ECHREAD itself, in ECHREAD only
    !   the transformation of the input array to 0.5 degree is done.
    !   Therefor, new calling parameter and routine names are set/given.
    !          alt: SUBROUTINE echread(luinp, ihead, field_out, istep, lhd_que)
    !          neu: SUBROUTINE hydrology_model(field_in, field_out, lhd_que)
    !
    ! ****** List of variables
    !
    !  field_out = Interpolated, transposed Array
    !  lhd_que   = Log-output switch  ( 0 = No Log-Output )

    REAL(dp), INTENT(in) :: field_in(nlon,ngl)
    REAL(dp), INTENT(out) :: field_out(nl,nb)
    LOGICAL, INTENT(in) :: lhd_que

    REAL(dp) :: acode(nlon + 1,ngl + 2), axlon(nlon + 1), axlat(ngl + 2)
    REAL(dp) :: xr(nl), yr(nb)

    INTEGER :: jlat, jlon


    CALL intpol_coord_axis_setup(axlon, axlat, xr, yr)

    IF (lhd_que) THEN
      WRITE(message_text,*) philat(1), philat(2), philat(ngl)
      CALL message('hydrology_echam', message_text)
      WRITE(message_text,*) axlat(1), axlat(2), axlat(ngl + 2)
      CALL message('hydrology_echam', message_text)

      WRITE(message_text,*) philon(1), philon(2), philon(nlon)
      CALL message('hydrology_echam', message_text)
      WRITE(message_text,*) axlon(1), axlon(2), axlon(nlon + 1)
      CALL message('hydrology_echam', message_text)

      WRITE(message_text,*) xr(1), xr(2), xr(nl)
      CALL message('hydrology_echam', message_text)

      WRITE(message_text,*) yr(1), yr(2), yr(nb)
      CALL message('hydrology_echam', message_text)
    END IF

    ! generate a copy of field_in with cyclic extension of longitudes and extra
    ! latitudes in the north and south

    DO jlat = 1, ngl
      DO jlon = 1, nlon
        acode(jlon,jlat+1) = field_in(jlon,jlat)
      ENDDO
    ENDDO

    DO jlat = 2, ngl + 1
      acode(nlon + 1,jlat) = acode(1,jlat)
    ENDDO

    DO jlon = 1, nlon + 1
      acode(jlon,1)      = acode(jlon,2)
      acode(jlon,ngl + 2) = acode(jlon,ngl + 2-1)
    ENDDO

    ! interpolation to output grid
    ! (includes transformation on bounded coordinates)
    CALL intpol_with_mapping(nlon + 1, ngl + 2, acode, axlon, axlat, &
         nl, nb, field_out, xr, yr, intpol_mapping)

  END SUBROUTINE hydrology_echam

  SUBROUTINE intpol_coord_axis_setup(axlon, axlat, xr, yr)
    REAL(dp), INTENT(inout) :: axlon(nlon + 1), axlat(ngl + 2), xr(nl), yr(nb)

    INTEGER :: j

    ! definition of the echam gp data grid
    axlat(1)       =  90.0_dp
    axlat(2:ngl+1) =  philat(1:ngl)
    axlat(ngl + 2) = -90.0_dp

    axlon(1:nlon)  = philon(:)
    axlon(nlon + 1) = fullcirc

    ! definition of hd data grid
    DO j = 1, nl
      xr(j) = 0.5_dp*hd_scal_lon + (j-1)/REAL(nl,dp)*fullcirc
    END DO

    DO j = 1, nb
      ! perhaps use this formulation? (tj, 20091009)
      ! yr(j) = -(0.5*hd_scal_lat + fullcirc * hd_scal_lat * ( (j-1) &
      !          / REAL(nb,dp) - 0.5))
      yr(j) = -(0.5_dp*hd_scal_lat - 0.5_dp*hd_scal_lat * fullcirc &
           &    + hd_scal_lat * fullcirc * REAL(j-1, dp) / REAL(nb,dp))
    ENDDO
  END SUBROUTINE intpol_coord_axis_setup


  !> interpolate field src to dest
  SUBROUTINE intpol (src_i_size, src_j_size, src, src_x, src_y, &
       dest_i_size, dest_j_size, dest, dest_x, dest_y)

    INTEGER, INTENT(in) :: src_i_size, src_j_size, dest_i_size, dest_j_size

    REAL(dp), INTENT(in)  :: dest_x(dest_i_size), dest_y(dest_j_size)
    REAL(dp), INTENT(in)  :: src_x(src_i_size), src_y(src_j_size)
    REAL(dp), INTENT(in)  :: src(src_i_size,src_j_size)
    REAL(dp), INTENT(out) :: dest(dest_i_size,dest_j_size)

    INTEGER :: irun, js, jd, is, id, idt

    irun = 0
    DO js = 2, src_j_size
      DO jd = 1, dest_j_size
        ! if (dest_y(jd) < src_y(js-1) .or. dest_y(jd) > src_y(js)) cycle
        IF (dest_y(jd) < MIN(src_y(js - 1), src_y(js)) .OR.   &
             dest_y(jd) > MAX(src_y(js - 1), src_y(js))) CYCLE
        irun = irun+1
        DO is = 2, src_i_size
          DO id = 1, dest_i_size
            IF(dest_x(id) < src_x(is-1) .OR. dest_x(id) > src_x(is)) CYCLE
            idt = MOD(id - 1 + dest_i_size/2, dest_i_size) + 1
            dest(idt,jd) = src(is - 1, js - 1) &
                 * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js)) &
                 / ((src_x(is - 1) - src_x(is)) * (src_y(js - 1) - src_y(js))) &
                 + src(is, js - 1) &
                 * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js)) &
                 / ((src_x(is) - src_x(is - 1)) * (src_y(js - 1) - src_y(js))) &
                 + src(is - 1, js) &
                 * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js - 1)) &
                 / ((src_x(is - 1) - src_x(is)) * (src_y(js) - src_y(js - 1))) &
                 + src(is, js) &
                 * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js - 1)) &
                 / ((src_x(is) - src_x(is - 1)) * (src_y(js) - src_y(js - 1)))
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE intpol

  !> find midpoint is,js for every index pair from field
  !> src(src_i_size, src_j_size) to dest(dest_i_size, dest_j_size)
  !> every element cart_idx_2d(is,js) of mapping(id,jd) later defines
  !> that dest(id,jd) will be computed from
  !> src(is,js), src(is - 1, js), src(is, js - 1), src(is - 1, js - 1)
  SUBROUTINE intpol_compute_mapping(mapping, &
       src_i_size, src_j_size, src_x, src_y, &
       dest_i_size, dest_j_size, dest_x, dest_y)

    INTEGER, INTENT(in) :: src_i_size, src_j_size, dest_i_size, dest_j_size

    REAL(dp), INTENT(in)  :: dest_x(dest_i_size),dest_y(dest_j_size)
    REAL(dp), INTENT(in)  :: src_x(src_i_size),src_y(src_j_size)
    TYPE(cart_idx_2d), INTENT(out) :: mapping(dest_i_size,dest_j_size)

    INTEGER :: js, jd, is, id

    mapping = cart_idx_2d(-1, -1)
    DO js = 2, src_j_size
      DO jd = 1, dest_j_size
        ! if (dest_y(jd) < src_y(js-1) .or. dest_y(jd) > src_y(js)) cycle
        IF (dest_y(jd) < MIN(src_y(js - 1), src_y(js)) .OR.   &
            dest_y(jd) > MAX(src_y(js - 1), src_y(js))) CYCLE
        DO is = 2, src_i_size
          DO id = 1, dest_i_size
            IF(dest_x(id) < src_x(is-1) .OR. dest_x(id) > src_x(is)) CYCLE
            mapping(id, jd) = cart_idx_2d(is, js)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE intpol_compute_mapping

  !> interpolate field src to dest
  SUBROUTINE intpol_with_mapping(src_i_size, src_j_size, src, src_x, src_y, &
       dest_i_size, dest_j_size, dest, dest_x, dest_y, mapping)

    INTEGER, INTENT(in) :: src_i_size, src_j_size, dest_i_size, dest_j_size

    REAL(dp), INTENT(in)  :: dest_x(dest_i_size),dest_y(dest_j_size)
    REAL(dp), INTENT(in)  :: src_x(src_i_size),src_y(src_j_size)
    REAL(dp), INTENT(in)  :: src(src_i_size,src_j_size)
    REAL(dp), INTENT(out) :: dest(dest_i_size,dest_j_size)
    TYPE(cart_idx_2d), INTENT(in) :: mapping(dest_i_size, dest_j_size)
    INTEGER :: js, jd, is, id, idt

    dest(:,:) = 0._dp

    DO jd = 1, dest_j_size
      DO id = 1, dest_i_size
        IF (mapping(id,jd)%ilat /= -1) THEN
          is = mapping(id,jd)%ilon
          js = mapping(id,jd)%ilat
          idt = id + dest_i_size/2 - MERGE(dest_i_size, 0, id > dest_i_size/2)
          dest(idt,jd) = src(is - 1, js - 1) &
               * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js)) &
               / ((src_x(is - 1) - src_x(is)) * (src_y(js - 1) - src_y(js))) &
               + src(is, js - 1) &
               * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js)) &
               / ((src_x(is) - src_x(is - 1)) * (src_y(js - 1) - src_y(js))) &
               + src(is - 1, js) &
               * (dest_x(id) - src_x(is)) * (dest_y(jd) - src_y(js - 1)) &
               / ((src_x(is - 1) - src_x(is)) * (src_y(js) - src_y(js - 1))) &
               + src(is, js) &
               * (dest_x(id) - src_x(is - 1)) * (dest_y(jd) - src_y(js - 1)) &
               / ((src_x(is) - src_x(is - 1)) * (src_y(js) - src_y(js - 1)))
        END IF
      END DO
    END DO
  END SUBROUTINE intpol_with_mapping

  FUNCTION reassign_runoff(nlon, nlat, jlon, jlat, foclsm, fatmos, fdat) &
       RESULT(reassigned)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nlon, nlat, jlon, jlat
    REAL(dp), INTENT(in) :: foclsm(nlon, nlat), fatmos(nlon, nlat)
    REAL(dp), INTENT(inout) :: fdat
    LOGICAL :: reassigned
    !
    ! ndd = If no land point is found as direct neighbour,
    !       it is searched in NWSE direction until the maximum distance of
    !       NDD Boxes is reached.

    INTEGER, PARAMETER :: ndd = 3

    REAL(dp) :: x1, inverted_neighbour_weight
    INTEGER :: idd

    reassigned = .FALSE.
    ! HD Land but OA Water
    ! Considered neighbour gridboxes in OA grid
    ! N,S,W,E-Directions
    x1 = 0.0_dp
    inverted_neighbour_weight = 0.0_dp
    IF (jlon /= 1) THEN
      IF (foclsm(jlon-1,jlat) > 0.5_dp) THEN
        x1 = fatmos(jlon-1,jlat)
        inverted_neighbour_weight = 1.0_dp
      ENDIF
    ELSE
      IF (foclsm(nlon,jlat) > 0.5_dp) THEN
        x1 = fatmos(nlon,jlat)
        inverted_neighbour_weight = 1.0_dp
      ENDIF
    ENDIF
    IF (jlon /= nlon) THEN
      IF (foclsm(jlon+1,jlat) > 0.5_dp) THEN
        x1 = x1+fatmos(jlon+1,jlat)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ELSE
      IF (foclsm(1,jlat) > 0.5_dp) THEN
        x1 = x1+fatmos(1,jlat)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ENDIF
    IF (jlat /= 1) THEN
      IF (foclsm(jlon,jlat-1) > 0.5_dp) THEN
        x1 = x1+fatmos(jlon,jlat-1)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ENDIF
    IF (jlat /= nlat) THEN
      IF (foclsm(jlon,jlat+1) > 0.5_dp) THEN
        x1 = x1+fatmos(jlon,jlat+1)
        inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
      ENDIF
    ENDIF
    ! Land point found?
    IF (inverted_neighbour_weight > 0.5_dp) THEN
      fdat = x1/inverted_neighbour_weight
      reassigned = .TRUE.
      RETURN
    END IF

    IF (jlon /= 1) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon-1,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon-1,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon-1,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon-1,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(nlon,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(nlon,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(nlon,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(nlon,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ENDIF
    IF (jlon /= nlon) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon+1,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+1,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon+1,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+1,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(1,jlat-1) > 0.5_dp) THEN
          x1 = x1+fatmos(1,jlat-1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(1,jlat+1) > 0.5_dp) THEN
          x1 = x1+fatmos(1,jlat+1)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
    ENDIF
    IF (inverted_neighbour_weight > 0.5_dp) THEN
      ! Second next points in OA grid in N,S,W,E-Directions
      fdat = x1/inverted_neighbour_weight
      reassigned = .TRUE.
      RETURN
    END IF
    extended_surround_loop: DO idd = 2, ndd
      ! HD Land but OA Water
      x1 = 0.0_dp
      inverted_neighbour_weight = 0.0_dp
      IF (jlon-idd >= 1) THEN
        IF (foclsm(jlon-idd,jlat) > 0.5_dp) THEN
          x1 = fatmos(jlon-idd,jlat)
          inverted_neighbour_weight = 1.0_dp
        ENDIF
      ELSE
        IF (foclsm(nlon+jlon-idd,jlat) > 0.5_dp) THEN
          x1 = fatmos(nlon+jlon-idd,jlat)
          inverted_neighbour_weight = 1.0_dp
        ENDIF
      ENDIF
      IF (jlon+idd <= nlon) THEN
        IF (foclsm(jlon+idd,jlat) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+idd,jlat)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ELSE
        IF (foclsm(jlon+idd-nlon,jlat) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon+idd-nlon,jlat)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat-idd >= 1) THEN
        IF (foclsm(jlon,jlat-idd) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon,jlat-idd)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      IF (jlat+idd <= nlat) THEN
        IF (foclsm(jlon,jlat+idd) > 0.5_dp) THEN
          x1 = x1+fatmos(jlon,jlat+idd)
          inverted_neighbour_weight = inverted_neighbour_weight + 1.0_dp
        ENDIF
      ENDIF
      ! End of Do (IDD) -Loop for Land Point found
      IF (inverted_neighbour_weight > 0.5_dp) THEN
        fdat = x1/inverted_neighbour_weight
        reassigned = .TRUE.
        EXIT extended_surround_loop
      END IF
    END DO extended_surround_loop
    ! end of HD land, OA water
  END FUNCTION reassign_runoff

  SUBROUTINE kasglob(finp, ymod, a_k, a_n, steps_per_day, fmem, a_n_kas)
    !
    ! ***** Global Flow Simulation with the conceptual model reservoir cascade
    !   Program was partailly written following the routine lfsim (in gate.for)
    !   and funkas (in modfunct.for).
    !
    ! ***** Programmed and developed by Stefan Hagemann
    !
    ! ***** Version 2.2 - November 1995
    !   Finer resolution of system function for Riverflow
    !   Re-arranging of IN/OUTput-arrays and Passing within a
    !      single reservoir field fmem
    !
    ! ***** Version 3.0 - November 1995
    !   Computation of outflow via Differential Equation
    !   of lin. reservoir cascade, comprising of nn reservoirs with
    !   nn = INT(a_n) = INT(n)
    !
    ! ***** Version 3.1 - Februar 1996
    !   Implementation of possible computation of riverflow with mm
    !   sub (internal) time steps.
    !
    ! ***** Version 5.0 - Oktober 1999
    !   Runoff-Input-data are passed to kasglob instead of reading it within
    !   kasglob itself.
    !   Calling parameters/variables luinp and area are deleted.
    !

    INTEGER,  INTENT(in)           :: steps_per_day  ! number of routine calls per day
    REAL(dp), INTENT(in)           :: a_k(nl, nb)    ! array of k-parameter [day]
    REAL(dp), INTENT(in)           :: a_n(nl, nb)    ! array of n-parameter
    REAL(dp), INTENT(in)           :: finp(nl, nb)   ! input overlandflow/riverflow array
    REAL(dp), INTENT(inout)        :: ymod(nl, nb)   ! simulated overlandflow/riverflow
    REAL(dp), INTENT(inout)        :: fmem(:,:,:)    ! intermediate content of reservoir cascade

    TYPE(cart_xidx_2d), INTENT(in) :: a_n_kas(:)     ! kasglob list

    REAL(dp) :: akdiv, fdum, amod, divmm, fmd_sum

    INTEGER :: j, jl, jb
    INTEGER :: nn
    INTEGER :: nx, i, extlen, extelem

    divmm = 1._dp/steps_per_day

    nx = SIZE(a_n_kas)
    ! **** Computing modeled value at each grid point

    DO i = 1, nx
      jb = a_n_kas(i)%ilat
      jl = a_n_kas(i)%ilon
      extlen = a_n_kas(i)%extlen
      DO extelem = 1, extlen
        nn = NINT(a_n(jl,jb))
        amod = a_n_kas(i)%amod(extelem)
#if 0
        ! Test for the amod values defined in create_kasglob_list
        IF (ABS(amod) > ABS(a_k(jl,jb) * a_n(jl,jb)/AINT(a_n(jl,jb)))+EPSILON(1._dp)) THEN
          WRITE (message_text,*) 'amod problem: ', jl, jb, &
               amod, a_k(jl,jb) * a_n(jl,jb)/AINT(a_n(jl,jb)), extelem
          CALL message('hd: kasglob', message_text)
        END IF
#endif

        ! *** Dt=1 day ==> AKDIV = 1./ (AMOD+1.)
        akdiv = a_n_kas(i)%akdiv(extelem)*divmm
        fdum = finp(jl,jb)

        ! *** Nash-Cascade
        ! *** It is [AMOD] = day ==> AMOD(sec) = AMOD(day) * 1 day
        ! *** Remember: In principle,it is: FDUM = FINP * 1 day
        ! ***           ==> FMEM = x * 1 day
        ! ***           ==> FDUM = x * 1 day * 1 / AMOD(sec)
        ! ***                    = x * 1 day / (AMOD(day) * 1 day)
        ! ***                    = x / AMOD(day)
        ! ***           ==> FMEM = x * 1 day - FDUM * 1 day
        ! ***                    = (x - FDUM) * 1 day
        ! *** Outflow FDUM is computed correctly, Intermediate reservoir unit is
        ! *** a volume flow instead of a volume. This is to avoid
        ! *** back and forth multiplication with factor 1 day = 86400 sec

        DO j = 1, nn
          fmd_sum = fmem(jl,jb,j) + fdum
          fdum = fmd_sum * akdiv
          fmem(jl,jb,j) = fmd_sum - fdum
        END DO
        ymod(jl,jb) = fdum
        jl = jl + 1
      ENDDO
    ENDDO

  END SUBROUTINE kasglob

  SUBROUTINE hydrology_to_ocean(nlon, nlat, friv, fdir, disch)

    !
    ! ******* This programs distributes the river discharge from the 0.5
    !         degree inflow points to the considered ocean gridbox
    !
    !  Programmed and developed by Stefan Hagemann, MPI
    !
    ! ***** Version 1.0 -- Oktober 1999
    !
    ! ***** Version 2.0 -- January 2001
    !     ECHAM5- Version incl. Gaussian latitudes
    !
    ! ****** List of Variables
    !
    !  friv = Inflow array on HD model grid
    !  fdir = River direction file that defines river mouthes (destinations) as 0
    !         on HD model grid
    !  xidb = Summation array of inflows, for which no inflowbox into the
    !         ocean was found, e.g. Kaspian Sea and Interior
    !         Drainage Basins
    !  any_ocinflow = Inflow-Point found on ocean grid .TRUE./.FALSE.
    !
    !  nlon = Longitudes of global ocean grid
    !  nlat = Latitudes of global ocean grid
    !  oclorg = Longitudinal origin of global ocean grid
    !  ocscal = Scale/Resolution = Width of Ocean Gridbox in degree
    !  philat = Gaussian latitude of global ocean grid (centre coordinates)
    !
    !   disch = Inflow array on Ocean grid
    ! lhd_que = Log-Output switch (.FALSE. = No Log-output to STDOUT)
    !

    INTEGER,  INTENT(in)  :: nlon, nlat
    REAL(dp), INTENT(in)  :: friv(nl,nb)
    INTEGER,  INTENT(in)  :: fdir(nl,nb)
    REAL(dp), INTENT(out) :: disch(nlon,nlat)
#if 0
    REAL(dp) :: xjlat
#endif
    REAL(dp) :: xidb
    INTEGER :: jl,jb
    TYPE(cart_idx_2d) :: dest
    LOGICAL :: any_ocinflow

    xidb = 0.0_dp
    disch(:,:) = 0.0_dp
    any_ocinflow = .FALSE.

    ! ******* Loop over all inflow points

    DO jb = 1, nb
       DO jl = 1, nl

          ! HD ocean cell
          dest = oclook_cache(jl, jb)
          IF (dest%ilon /= -1) THEN
             disch(dest%ilon,dest%ilat) = disch(dest%ilon,dest%ilat) + friv(jl,jb)
             any_ocinflow = .TRUE.
          END IF

          ! internal drainage cells
          IF (fdir(jl,jb) == 5) THEN
             xidb = xidb + friv(jl,jb)
          END IF

      ENDDO
    ENDDO

    ! Distributing the water in XIDB to all Ocean Inflow Points
    ! Applying a weight to treat arid and humid regions differently

    IF (any_ocinflow) THEN
      disch(:,:) = disch(:,:)+disch(:,:)/SUM(disch(:,:)) * xidb
    ELSE
      WRITE(message_text,*) 'error no inflow points on ocean grid found'
      CALL message('hydrology_to_ocean', message_text)

    ENDIF

  END SUBROUTINE hydrology_to_ocean

  SUBROUTINE hydrology_slm_invariants(slm, slf, lake)

    USE mo_physical_constants, ONLY: earth_radius
    USE mo_math_constants,     ONLY: pi

    REAL(dp),       INTENT(in)    :: slm(:,:)       ! integer land sea mask
    REAL(dp),       INTENT(in)    :: slf(:,:)       ! fractional land sea mask
    REAL(dp),       INTENT(in)    :: lake(:,:)      ! fractional lake mask

    INTEGER  :: il, jb, i
    REAL(dp) :: ra, rb, rh, rd

    IF (p_parallel_io) THEN

       ! setup area in m^2 of the HD model internal grid

       ra = 2.0_dp*pi*earth_radius*earth_radius/REAL(nl,dp)
       rb = pi/REAL(nb,dp)
       rd = 0.5_dp*pi

       DO i = 1, nb
          rh = SIN(-rd+(i-1)*rb)-SIN(-rd+i*rb)
          hd_area(i) = ABS(rh)*ra
       END DO

       ! Gaussian latitudes in degrees

       IF (ldebughd) THEN
          WRITE(message_text,*) 'philon= ', philon(1),philon(2),philon(3),' ... ',philon(nlon) 
          CALL message('hydrology_model', message_text)
          WRITE(message_text,*) 'philat= ', philat(1),philat(2),philat(3),' ... ',philat(ngl) 
          CALL message('hydrology_model', message_text)
       END IF

       ! Definition of ocean land masks, in contrast to slm/slf with lakes represented as land

#ifdef ECHAM_FRACTIONAL
       no_ocean_mask(:,:)  = MERGE (1._dp, 0._dp, (slf(:,:) + lake(:,:) == 1.0_dp ))
       land_fract(:,:)     = slf(:,:)
       lake_fract(:,:)     = lake(:,:)
       ocean_fract(:,:)    = 1._dp - (slf(:,:) + lake(:,:))
#else
       no_ocean_mask(:,:)  = MERGE(slm(:,:), 1._dp, lake(:,:) < 0.5_dp)
       land_fract(:,:)     = slm(:,:)
       lake_fract(:,:)     = MERGE(1._dp, 0._dp, (lake(:,:) > 0.5_dp ))
       ocean_fract(:,:)    = 1._dp - no_ocean_mask(:,:)
#endif

       ! define array with number of reservoirs for the base flow.
       ! There is just 1 baseflow reservoir (linear) at each grid point. The array 
       ! is needed to be able to use routine kasglob.

       WHERE (hd_lsm > 0.5_dp)
          agf_n = 1._dp
       ELSEWHERE
          agf_n = 0._dp
       END WHERE

       ! check consistency of HD input data

       IF (ANY((hd_lsm > 0.5_dp .AND. (fdir == -1 .OR. fdir == 0)) &
            .OR. (hd_lsm < 0.5_dp .AND. (fdir /= -1 .AND. fdir /= 0 .AND. fdir /= 5)))) THEN
          CALL finish('hydrology_slm_invariants', &
               'hd slm and runoff directions do not match')
       END IF

       IF (ANY((hd_lsm > 0.5_dp .AND. arf_n < 0.5_dp .AND. fdir /= 5) &
            .OR. (hd_lsm < 0.5_dp .AND. arf_n > 0.5_dp))) THEN
          CALL finish('hydrology_slm_invariants', &
               'hd slm and riverflow reservoir numbers do not match')
       END IF

       IF (ANY((hd_lsm > 0.5_dp .AND. alf_n < 0.5_dp .AND. fdir /= 5) &
            .OR. (hd_lsm < 0.5_dp .AND. alf_n > 0.5_dp))) THEN
          CALL finish('hydrology_slm_invariants', &
               'hd slm and overlandflow reservoir numbers do not match')
       END IF

       IF (ANY((hd_lsm > 0.5_dp .AND. agf_n <= 0._dp) &
            .OR. (hd_lsm < 0.5_dp .AND. agf_n /= 0._dp))) THEN
          CALL finish('hydrology_slm_invariants', &
               'hd slm and baseflow reservoir numbers do not match')
       END IF

       IF (ANY((fdir == 5._dp) .AND. agf_k > 0._dp )) THEN
          CALL finish('hydrology_slm_invariants', &
               'baseflow over interior darinage basin')
       END IF

       ! grid cell centers

       IF (calculate_lonlat_hd) THEN
          DO il = 1, nl
             lon_hd(il,:) = florg + REAL(il,dp)*fscal - 0.5_dp*fscal
          END DO
          DO jb = 1, nb
             lat_hd(:,jb) = fborg - REAL(jb,dp)*fscal + 0.5_dp*fscal
          END DO
       END IF

       IF (do_remapping /= 0) THEN
          CALL fill_oclook_caches(no_ocean_mask, fdir)
          CALL read_remap_matrix
       ENDIF

       CALL create_kasglob_list(alf_n, alf_k, hd_steps_per_day, alf_n_kas)
       CALL create_kasglob_list(agf_n, agf_k, hd_steps_per_day, agf_n_kas)
       CALL create_kasglob_list(arf_n, arf_k, riverflow_steps_per_day, arf_n_kas)
    END IF
  END SUBROUTINE hydrology_slm_invariants

  SUBROUTINE cleanup_hd_slm_invariants
    IF (p_parallel_io) THEN
      CALL cleanup_kasglob_list(alf_n_kas)
      CALL cleanup_kasglob_list(agf_n_kas)
      CALL cleanup_kasglob_list(arf_n_kas)
    END IF
  END SUBROUTINE cleanup_hd_slm_invariants

  SUBROUTINE set_riverflow_timestep
    !*************************************************************************
    !
    ! **** Routine that defines the internal time step used for riverflow
    !      calculations
    !
    USE mo_time_control, ONLY: ev_puthd, get_interval_seconds_next

    IF (p_parallel_io) THEN

       IF (.NOT. hydrology_offline) THEN
          hd_calling_interval = get_interval_seconds_next(ev_puthd)
       ELSE
          hd_calling_interval = delta_time
       END IF

       ! the riverflow time step should not be greater than 6 hours (= 21.600 s).
       ! The variable riverflow_timestep defines the number of internal timesteps
       ! per HD model time steps.

       IF (hd_calling_interval <= 21600) THEN
          riverflow_timestep = 1                  ! no internal riverflow time steps needed
       ELSE IF (hd_calling_interval <= 43200) THEN
          riverflow_timestep = 2                  ! two riverflow time steps per HD time step
       ELSE IF (hd_calling_interval <= 64800) THEN
          riverflow_timestep = 3
       ELSE IF (hd_calling_interval <= 86400) THEN
          riverflow_timestep = 4
       ELSE
          WRITE (message_text,*) 'The hydrology model should be called at least once a day. '&
               //'hd_calling_interval = ', hd_calling_interval, ' seconds.'
          CALL finish ('mo_hydrology: set_riverflow_timestep', message_text)
       END IF

       hd_steps_per_day = NINT(86400._dp/hd_calling_interval)
       riverflow_steps_per_day = hd_steps_per_day * riverflow_timestep
       div_riverflow_timestep = 1._dp/REAL(riverflow_timestep,dp)
    END IF

    IF (p_parallel) THEN
       CALL p_bcast(hd_calling_interval, p_io)
       CALL p_bcast(hd_steps_per_day, p_io)
       CALL p_bcast(riverflow_timestep, p_io)
       CALL p_bcast(riverflow_steps_per_day, p_io)
    END IF

    IF (ldebughd) THEN
       WRITE (message_text,*) 'hd_steps_per_day = ', hd_steps_per_day
       CALL message ('set_riverflow_timestep',message_text)
       WRITE (message_text,*) 'riverflow_steps_per_day = ', riverflow_steps_per_day
       CALL message ('set_riverflow_timestep',message_text)
    END IF


  END SUBROUTINE set_riverflow_timestep

  SUBROUTINE fill_oclook_caches(slm, fdir)

    !*************************************************************************
    !
    ! **** For each HD grid cell the routine finds the closest ocean grid cell
    !      on the ECHAM grid.
    !
    !      First, the Echam grid cell is found, in which the HD grid cell falls.
    !      Then the closest ocean grid cell is seached within a range of ndd
    !      Echam grid cells.

    REAL(dp), INTENT(in) :: slm(nlon, ngl)    ! ocean land mask on ECHAM grid
    INTEGER,  INTENT(in) :: fdir(nl, nb)      ! runoff directions on HD grid

    REAL(dp) :: fb, fl                        ! latitudes/longitudes on HD grid
    INTEGER  :: jb, jl                        ! indices of HD grid latitudes/longirudes
    INTEGER  :: jlon, jlat                    ! indices of the ECHAM grid cell, in which HD coordinates (jl,jb) fall.
    TYPE(cart_coord_2d) :: coord              ! REAL(dp) indices of HD cells relative to the ECHAM indices. 

    DO jb = 1, nb
       fb = lat_hd(1,jb)
       coord%latitude = fb
       jlat = dec_monotonic_closest_midpoint(philat, fb, aub=90._dp, alb=-90._dp)
       DO jl = 1, nl
          fl = MOD(lon_hd(jl,jb) + fullcirc, fullcirc)
          jlon = inc_monotonic_closest_midpoint(philon, fl, alb=0._dp, aub=360._dp)
          coord%longitude = fl
          coord%dir = fdir(jl,jb)
          oclook_cache(jl,jb) = oclook(slm, jlon, jlat, coord)
       END DO
    END DO
  END SUBROUTINE fill_oclook_caches

  SUBROUTINE create_kasglob_list(a_n, a_k, steps_per_day, a_n_kas)
    !
    REAL(dp), INTENT(in) :: a_n(:, :), a_k(:, :)
    INTEGER, INTENT(in) :: steps_per_day
    TYPE(cart_xidx_2d), ALLOCATABLE, INTENT(inout) :: a_n_kas(:)
    REAL(dp) :: divmm
    INTEGER :: size_i, size_j, i, j, num_extents, extent, extlen
    INTEGER :: jl, extelem

    divmm = 1._dp/REAL(steps_per_day,dp)

    size_i = SIZE(a_n, 1)   ! number of longitudes of the HD-grid
    size_j = SIZE(a_n, 2)   ! number of latitudes of the HD-grid
    num_extents = 0
    IF (ALLOCATED(a_n_kas)) THEN
      CALL cleanup_kasglob_list(a_n_kas)
    END IF

    ! Each grid row is devided into pieces with positive reservoir numbers.
    ! The number of these pieces is counted (num_extents).  
    DO j = 1, size_j
      i = 1
      DO WHILE(i <= size_i)
        IF (a_n(i, j) > 0.5_dp) THEN
          num_extents = num_extents + 1
          i = i + 1
          DO WHILE(i <= size_i)
            IF (a_n(i, j) <= 0.5_dp) EXIT
            i = i + 1
          END DO
        ELSE
          i = i + 1
        END IF
      END DO
    END DO

    ! Each stripe with positive reservoir numbers is filled with data:
    !   ilat: y-index of the first cell of the stripe
    !   ilon: x-index of the first cell of the stripe
    !   extlen: number of grid cells in the stripe
    !   amod: amod for each grid cell in the stripe
    !   akdiv: value for each grid cell in the stripe (depending on time step)

    ALLOCATE(a_n_kas(num_extents))
    extent = 0
    DO j = 1, size_j
      i = 1
      DO WHILE(i <= size_i)
        IF (a_n(i, j) > 0.5_dp) THEN
          extent = extent + 1
          extlen = 0
          a_n_kas(extent)%ilat = j
          a_n_kas(extent)%ilon = i
          DO WHILE(i <= size_i)
            IF (a_n(i, j) <= 0.5_dp) EXIT
            i = i + 1
            extlen = extlen + 1
          END DO
          a_n_kas(extent)%extlen = extlen
          ALLOCATE(a_n_kas(extent)%amod(extlen), a_n_kas(extent)%akdiv(extlen))
          DO extelem = 1, extlen
            jl = a_n_kas(extent)%ilon + extelem - 1
            a_n_kas(extent)%amod(extelem) = a_k(jl,j) * a_n(jl,j) &
                 / AINT(a_n(jl,j))
            a_n_kas(extent)%akdiv(extelem) = 1.0_dp &
                 / (a_n_kas(extent)%amod(extelem) + divmm)
          END DO
        ELSE
          i = i + 1
        END IF
      END DO
    END DO
  END SUBROUTINE create_kasglob_list

  SUBROUTINE cleanup_kasglob_list(a_n_kas)
    TYPE(cart_xidx_2d), ALLOCATABLE, INTENT(inout) :: a_n_kas(:)
    INTEGER :: i, n
    n = SIZE(a_n_kas)
    DO i = 1, n
      DEALLOCATE(a_n_kas(i)%amod, a_n_kas(i)%akdiv)
    END DO
    DEALLOCATE(a_n_kas)
  END SUBROUTINE cleanup_kasglob_list

  SUBROUTINE hd_remap(nlon_src, nlat_src, src_array, nlon_dst, nlat_dst, global_corr, &
                      src_area, dst_area, dst_array)

    ! do SCRIP remapping: first order conservative remapping 

    INTEGER,  INTENT(in)  :: nlon_src, nlat_src           ! source grid dimensions
    INTEGER,  INTENT(in)  :: nlon_dst, nlat_dst           ! target grid dimensions
    REAL(dp), INTENT(in)  :: src_array(nlon_src,nlat_src) ! field on source grid
    REAL(dp), INTENT(in)  :: src_area(nlon_src,nlat_src)  ! grid cell area on source grid
    REAL(dp), INTENT(in)  :: dst_area(nlon_dst,nlat_dst)  ! grid cell area on target grid
    LOGICAL,  INTENT(in)  :: global_corr                  ! switch for global conservation
    REAL(dp), INTENT(out) :: dst_array(nlon_dst,nlat_dst) ! field on target grid

    REAL(dp) :: src(nlon_src*nlat_src)  ! source field as one dimensional array
    REAL(dp) :: dst(nlon_dst*nlat_dst)  ! destination field as one dimensional array
    REAL(dp) :: integral_src
    REAL(dp) :: integral_dst

    INTEGER :: n

    src = RESHAPE(src_array,(/nlon_src*nlat_src/))
    dst = 0._dp
    DO n = 1, num_links
       dst(dst_address(n)) = dst(dst_address(n)) + src(src_address(n)) * remap_matrix(1,n)
    END DO

    dst_array = RESHAPE(dst,(/nlon_dst,nlat_dst/))

    ! assure global conservation
    IF (global_corr) THEN
       integral_src = SUM(src_array * src_area)
       integral_dst = SUM(dst_array * dst_area)
       dst_array = dst_array * integral_src/integral_dst
    END IF

  END SUBROUTINE hd_remap

  SUBROUTINE read_remap_matrix
    !
    ! read matrix for remapping of data on the ECHAM grid to the HD model grid
    ! The matrix was generated offline with the SCRIP library (cdo gencon).

    TYPE (FILE_INFO)  :: fileinfo
    INTEGER           :: dimid, varid, fileid
    CHARACTER(len=80) :: rmpfile
    LOGICAL           :: lex
    REAL(dp), ALLOCATABLE :: array_dp(:)    !! double precision dummy array

    ! File names

    rmpfile = 'rmp_hd.nc'

    ! Read parameter: Land sea mask, RDF, ...

    INQUIRE (file=rmpfile, exist=lex)
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(rmpfile),'>'
      CALL finish ('read_remap_matrix', message_text)
    ENDIF

    fileinfo%opened = .FALSE.
    CALL IO_open (rmpfile, fileinfo, IO_READ)
    WRITE (message_text,*) 'Reading remap matrix from file ', TRIM(rmpfile)
    CALL message('read_remap_matrix', message_text)

    fileID = fileinfo%file_id

    CALL IO_inq_dimid (fileID, 'num_links', dimid)
    CALL IO_inq_dimlen (fileID, dimid, num_links)
    CALL IO_inq_dimid (fileID, 'num_wgts', dimid)
    CALL IO_inq_dimlen (fileID, dimid, num_weights)

    ALLOCATE (src_address(num_links))    !! address indices of the source grid
    ALLOCATE (dst_address(num_links))    !! address indices of the destination grid
    ALLOCATE (remap_matrix(num_weights,num_links))
    ALLOCATE (array_dp(num_links))       !! double precision dummy array

    CALL IO_inq_varid (fileID, 'src_address', varid)
    CALL IO_get_var_double (fileID, varid, array_dp)
    src_address = NINT(array_dp)
    CALL IO_inq_varid (fileID, 'dst_address', varid)
    CALL IO_get_var_double (fileID, varid, array_dp)
    dst_address = NINT(array_dp)
    CALL IO_inq_varid (fileID, 'remap_matrix', varid)
    CALL IO_get_var_double (fileID, varid, remap_matrix)

    CALL IO_close(fileinfo)
    DEALLOCATE (array_dp)
  
  END SUBROUTINE read_remap_matrix

  SUBROUTINE cleanup_remap_matrix

    DEALLOCATE(src_address, dst_address, remap_matrix)

  END SUBROUTINE cleanup_remap_matrix

  FUNCTION oclook(olm, jlon, jlat, coord) RESULT(idx)

    !*************************************************************************
    !
    ! **** Routine that looks for the next ocean gridbox on the ECHAM grid.
    !      The HD grid cell center is represented by the coord structure.

    ! Input

    REAL(dp),            INTENT(in) :: olm(:,:)  ! ocean land mask on ECHAM grid  
    INTEGER,             INTENT(in) :: jlon      ! index of ECHAM grid longitude
    INTEGER,             INTENT(in) :: jlat      ! index of ECHAM grid latitude
    TYPE(cart_coord_2d), INTENT(in) :: coord     ! indices, coordinates and mask
                                                 ! value of the HD grid cell
    ! Result

    TYPE(cart_idx_2d) :: idx                     ! nearest ocean grid cell found

    ! Local parameters

    REAL(dp), PARAMETER :: deg2rad = 1.74532925199432957692e-2  ! Degree to rad: 2pi/360

    REAL(dp) :: dx                           ! distance between ECHAM and HD grid cell centers
    REAL(dp) :: dxmin                        ! distance to the closet ECHAM ocean cell
    INTEGER  :: ndd                          ! maximum search distance (in grid cells)
    INTEGER  :: nlon, nlat                   ! number of ECHAM longitudes, latitudes
    INTEGER  :: ii, i1, i2, j1, j2, i, j, jj ! longitude, latitude indices
    REAL(dp) :: lon1, lat1, lon2, lat2       ! spherical coordinates   
    REAL(dp) :: x1, y1, z1, x2, y2, z2       ! x,y,z-coordinates   

    !
    nlon = SIZE(olm, 1)   ! number of ECHAM longitudes 
    nlat = SIZE(olm, 2)   ! number of ECHAM latiudes
 
    ! Nothing needs to be done for HD non-ocean cells

    IF (coord%dir /= -1 .AND. coord%dir /= 0) THEN  ! land or internal drainage
       idx%ilon = -1
       idx%ilat = -1
       RETURN
    END IF

    ! Check whether central ECHAM cell is an ocean grid cell

    IF (olm(jlon,jlat) < 0.5_dp) THEN  ! ocean cell
       idx%ilon = jlon
       idx%ilat = jlat
       RETURN
    END IF

    !
    !  Check the surrounding box with a diameter of ndd grid cells
    !
    ndd = INT(nlat/16) ! maximum distance: T31 3, T63 6 grid cells
    
    ! cell indicees for the search
    i1 = jlon - ndd
    i2 = jlon + ndd
    j1 = jlat - ndd
    j2 = jlat + ndd

    ! initializations
    dxmin = HUGE(dp)
    idx%ilon = -1
    idx%ilat = -1

    lon1 = coord%longitude * deg2rad  ! HD longitude [rad]
    lat1 = coord%latitude  * deg2rad  ! HD latitude [rad]

    ! find the nearest ocean grid cell
    DO ii = i1,i2
       DO jj = j1,j2
          i = MOD(nlon-1 + ii, nlon) + 1
          j = jj
          IF (jj < 1) THEN                   ! search radius south of S-pole (or north of N-pole):
            j = ABS(jj) + 1                        !  continue search on the opposite hemisphere,
            i = MOD(i-1 + NINT(nlon/2.), nlon) + 1 !  shifted by 180 deg.
          ELSE IF (jj > nlat) THEN
            j = nlat - (jj - nlat) + 1
            i = MOD(i-1 + NINT(nlon/2.), nlon) + 1
          END IF
          IF (olm(i,j) < 0.5_dp) THEN  ! ocean cell

             lon2 = philon(i) * deg2rad     ! ECHAM longitude [rad] 
             lat2 = philat(j) * deg2rad     ! ECHAM latitude [rad]
             !
             ! Transformation to x,y,z-coordinates
             !
             x1 = cos(lat1)*cos(lon1)
             y1 = cos(lat1)*sin(lon1)
             z1 = sin(lat1)

             x2 = cos(lat2)*cos(lon2)
             y2 = cos(lat2)*sin(lon2)
             z2 = sin(lat2)
             !
             ! Calculation of the distance
             ! 
             ! direct distance:
             dx = SQRT((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

             ! distance along the surface:
             dx = 2*ASIN(dx/2)

             IF (dx < dxmin) THEN
                dxmin = dx
                idx%ilon = i
                idx%ilat = j
             ENDIF
          END IF
       END DO
    END DO

    ! no ocean grid cell found within the search area?
    IF (dxmin == HUGE(dp)) THEN
       CALL message('', '')
       WRITE (message_text,*) 'no ocean cell found for HD grid cell with lon ',  &
            coord%longitude, ', lat ', coord%latitude, ' and dir ', coord%dir
       CALL message('oclook', message_text)
       WRITE (message_text,*) 'corresponding echam cell:', jlon, jlat, ' with olm: ', olm(jlon,jlat)
       CALL message('oclook', message_text)
       WRITE (message_text,*) 'search radius ndd =', ndd, ' needs to be increased.'
       CALL finish('oclook',  message_text)
    END IF
  END FUNCTION oclook

  !***********************************************************************
  ! Open outflow timeseries file
  !***********************************************************************
  SUBROUTINE hd_open_timeseries (isolog)

    USE mo_filename, ONLY: find_next_free_unit

    INTEGER,      INTENT(in)  :: isolog

    CHARACTER(80)             :: filename

    IF (p_parallel_io) THEN
       isolog_unit = find_next_free_unit (80, 100)
       WRITE (filename,'(a,i2.2,a)') 'hd_outflow_', isolog, '.log'
       OPEN (isolog_unit, file=TRIM(filename))
    ENDIF

  END SUBROUTINE hd_open_timeseries

  !***********************************************************************
  ! Open outflow timeseries file
  !***********************************************************************
  SUBROUTINE hd_close_timeseries

    IF (p_parallel_io) THEN
       CLOSE (isolog_unit)
    ENDIF

  END SUBROUTINE hd_close_timeseries

  !***********************************************************************
  ! Write outflow timeseries of specific grid cells (compare subroutine 
  ! hydrology_diags)
  !***********************************************************************
  SUBROUTINE hd_write_timeseries (f1, f2)

    USE mo_time_control,  ONLY: current_date, get_date_components

    REAL(dp),      INTENT(in) :: f1, f2
    INTEGER                   :: yr, mo, dy, hr, mn, se

    IF (p_parallel_io) THEN
       CALL get_date_components(current_date, yr, mo, dy, hr, mn, se)
       WRITE(isolog_unit,'(i6.4,i2.2,i2.2,a,i2.2,a,i2.2,a,i2.2,f14.4,f14.4)') &
            yr, mo, dy, ' ', hr, ':', mn, ':', se, F1, F2
    ENDIF

  END SUBROUTINE hd_write_timeseries

  !***********************************************************************
  ! Lateral Routing of discharge 
  !    Routine 1: as in original model at every riverflow time step
  !***********************************************************************
  SUBROUTINE routing(finp, fdata)
  !
  !        *** Conduct routing in a subroutine for easier exchange of methods
  !
  !  finp = local input data array for time step istep
  ! fdata = local output data array for time step istep
  ! finfl = Inflow data array for each gridbox for time step istep
  !  fdir = River direction array
  !
  ! **** Indices
  !
  !    jl = Longitudinal index
  !    jb = Latitudinal index
  !    il = relative change in longitude for the routing
  !    ib = relative change in latitude for the routing
  ! jlnew = jl+il
  ! jbnew = jb+ib

    REAL(dp), INTENT(in)  :: fdata(:,:)
    REAL(dp), INTENT(in)  :: finp(:,:)
    INTEGER :: jl, il, jlnew, jb, ib, jbnew, idir
!  
    ! re-initialize finfl
    finfl(:,:) = 0.0_dp

    ! ---------
    !  routing of outflow to finfl ==> new inflow per gridbox
    ! ---------

    DO jb = 1, nb
       DO jl = 1, nl
          idir = fdir(jl, jb) ! from HD parameter input file
          IF (idir > 0) THEN  ! internal land

             ! il, ib = relative direction coordinates [-1,0,1]

             ib = 1 - (idir - 1)/3
             il = MOD(idir - 1, 3) - 1

             jlnew = MOD(jl + il - 1 + nl, nl) + 1
             jbnew = jb + ib

          ELSE                ! ocean and coast
             jlnew = jl
             jbnew = jb
          END IF

          ! inflow of the new grid cell:  inflow from other cells 
          !                             + outflow from drainage and overlandflow (finp)
          !                             + outflow from riverflow calculations (fdata)

          finfl(jlnew,jbnew) = finfl(jlnew,jbnew) + finp(jl,jb)*div_riverflow_timestep + fdata(jl,jb)

       ENDDO
    ENDDO

  END SUBROUTINE routing

  !***********************************************************************
  ! Lateral Routing of discharge 
  !    Routine 2: routing via index arrays that are read in from HD parameter file
  !***********************************************************************
  SUBROUTINE routing_via_index(finp, fdata)
  !
  !        *** Conduct routing in a subroutine for easier exchange of methods
  !
  !  finp = local input data array for time step istep
  ! fdata = local output data array for time step istep
  ! finfl = Inflow data array for each gridbox for time step istep
  !  fdir = River direction array

    REAL(dp), INTENT(in)  :: fdata(:,:)
    REAL(dp), INTENT(in)  :: finp(:,:)
    INTEGER :: jl, jb, jlnew, jbnew
!  
    ! re-initialize finfl
    finfl(:,:) = 0.0_dp

    ! ---------
    !  routing of outflow to finfl ==> new inflow per gridbox
    ! ---------

    DO jb = 1, nb
       DO jl = 1, nl
         jlnew = filnew(jl, jb)
         jbnew = fibnew(jl, jb)

          ! inflow of the new grid cell:  inflow from other cells 
          !                             + outflow from drainage and overlandflow (finp)
          !                             + outflow from riverflow calculations (fdata)

          finfl(jlnew,jbnew) = finfl(jlnew,jbnew) + finp(jl,jb)*div_riverflow_timestep + fdata(jl,jb)
        ENDDO
     ENDDO

  END SUBROUTINE routing_via_index


END MODULE mo_hydrology
