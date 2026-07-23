!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_bethy
  !------------------------------------------------------------------------------
  ! WOLFGANG KNORR/GEORG HOFFMANN 050398
  ! MATTHIAS O'CUNTZ since 11/98
  ! Kalle Schnitzler MPI Hamburg, May 2000
  ! Reiner Schnur MPI Hamburg, 7/2003
  ! Christian Reick, MPI Hamburg, 8/2003
  !------------------------------------------------------------------------------

  USE mo_jsbach,          ONLY : debug
  USE mo_jsbach_grid,     ONLY : grid_type, domain_type, kstart, kend, nidx
  USE mo_netCDF,          ONLY : CANOPY
  USE mo_linked_list,     ONLY : t_stream
  USE mo_mpi,             ONLY : p_parallel_io, p_parallel, p_io, p_bcast
  USE mo_exception,       ONLY : message, finish
  USE mo_util_string,     ONLY : int2string
  USE mo_kind,            ONLY : dp
  USE mo_jsbach_lctlib,   ONLY : lctlib_type
  USE mo_bethy_photosyn,  ONLY : photosyn

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART =======================================================================================================

  ! --- public subroutines ---

  PUBLIC :: init_bethy, update_bethy, bethy_diagnostics

  TYPE canopy_type  !! Separate information for each canopy layer
     INTEGER                             :: ncanopy    ! Number of canopy layers
     INTEGER                             :: ntiles     ! Number of tiles
     REAL(dp), POINTER, DIMENSION(:,:,:) ::   & ! (nland, ncanopy, ntiles)
          lai, &                                    ! LAI per canopy layer
          canopy_conductance, &                     ! Canopy conductance per layer and per leaf area
          gross_assimilation,       &               ! GROSS PHOTOSYNTHESIS [MOL(CO2) / M^2(leaf area) S]
          dark_respiration,         &               ! DARK RESPIRATION OF LEAF [MOL(CO2) / M^2(leaf area) S]
          max_carbox_rate,          &
          max_e_transport_rate,     &
          carbox_rate,              &
          e_transport_rate,         &
          faPAR                         ! Absorbed PAR per canopy layer [fraction of irradiance]
  END TYPE canopy_type

  TYPE bethy_type
     INTEGER                           :: ntiles  ! Number of tiles
     TYPE(canopy_type)                 :: Canopy
     REAL(dp), POINTER, DIMENSION(:,:) :: & ! (nland, ntiles)
          CO2_conc_leaf,                  & ! CO2 concentration inside leaf [MOL(CO2)/MOL(AIR)]
          CO2_conc_leaf_acc,              & ! CO2 concentration inside leaf [MOL(CO2)/MOL(AIR)] (mean value)
          CO2_conc_leaf_unlimited_acc,    & ! CO2 concentration inside leaf for no water stress [MOL(CO2)/MOL(AIR)] (mean value)
          canopy_conductance,             & ! Canopy conductance for water vapour [m/s] (mean value)
          canopy_conductance_unlimited,   & ! Unlimited canopy conductance for water vapour [m/s] (mean value)
          gross_assimilation_acc,         & ! GROSS PHOTOSYNTHESIS [MOL(CO2) / M^2(ground) S] (mean value)
          dark_respiration_acc,           & ! DARK RESPIRATION OF LEAF [MOL(CO2) / M^2(ground) S] (mean value)
          max_carbox_rate,                & ! Maximum carboxylation rate (=VCmax) [MOL(CO2)/M^2 S]
          max_e_transport_rate,           & ! Maximum electron transport rate (=Jmax) (only C3 plants)[MOL(CO2)/M^2 S]
          carbox_rate,                    & ! Actual carboxylation rate (=JC) [MOL(CO2)/M^2 S]
          e_transport_rate,               & ! Actual electron transport rate (=JE) [MOL(CO2)/M^2 S]
          net_assimilation,               & ! Canopy net carbon assimilation [MOL(CO2)/m^2 s] (mean value)
          gross_assimilation,             & ! GROSS PHOTOSYNTHESIS at each time step [MOL(CO2) / M^2(ground) S]
          dark_respiration,               & ! DARK RESPIRATION OF LEAF at each time step [MOL(CO2) / M^2(ground) S]
          faPAR,                          & ! fraction of PAR absorbed by the canopy
          apar_acc,                       & ! PAR absorbed by the canopy [MOL PHOTONS/M^2 S] (mean value)
          par_acc                           ! incoming PAR [MOL PHOTONS/M^2 S] (mean value)
     REAL(dp), POINTER, DIMENSION(:) ::   &
          CO2_flux_net_acc,               & ! net CO2 flux to the atmosphere
          CO2_flux_herbivory_acc,         & ! CO2 flux to the atmosphere from grazing
          CO2_emission_landcover_change_acc, & ! CO2 flux to the atmosphere from land cover changes
          CO2_emission_harvest_acc,        & ! CO2 flux to the atmosphere from harvest
          CO2_flux_dynveg_acc                ! CO2 flux to the atmosphere from vegetation dynamics (disturbance?)
  END TYPE bethy_type

  TYPE bethy_diag_type
     REAL(dp), POINTER, DIMENSION(:)   :: & ! (nland)
          net_assimilation,               & ! Canopy net carbon assimilation [MOL(CO2)/m^2 s] (mean value)
          gross_assimilation,             & ! GROSS PHOTOSYNTHESIS [MOL(CO2) / M^2(ground) S] (mean value)
          dark_respiration,               & ! DARK RESPIRATION OF LEAF [MOL(CO2) / M^2(ground) S] (mean value)
          faPAR,                          & ! fraction of PAR absorbed by the canopy
          apar_acc,                       & ! PAR absorbed by the canopy [MOL PHOTONS/M^2 S] (mean value)
          par_acc                           ! incoming PAR [MOL PHOTONS/M^2 S] (mean value)
  END TYPE bethy_diag_type

  PUBLIC :: bethy_type, canopy_type

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE ! Make ALL following objects private

  ! === Private declarations =======================================================================================================

  INTEGER, SAVE             ::  ncanopy = -1                 ! Number of canopy layers

  LOGICAL, SAVE :: module_configured = .FALSE. 

  ! --- Stream variables ----------------------

  TYPE(t_stream), SAVE, POINTER :: IO_bethy                 !! Bethy stream
  TYPE(t_stream), SAVE, POINTER :: IO_diag                  !! Bethy diagnostic stream

  TYPE(bethy_diag_type), SAVE    :: bethy_diag

  REAL(dp), SAVE, POINTER :: canopy_boundaries_lai(:)       !! Canopy layers (ncanopy+1)
                                                            !! .. [mol(CO2)/mol(air) = 10^6 ppm] (nland x nlct)

CONTAINS

  SUBROUTINE config_bethy

    USE mo_jsbach,           ONLY: nml_unit
    USE mo_namelist,         ONLY: position_nml, POSITIONED
    USE mo_io_units,         ONLY: nout

    INTEGER :: read_status, f_unit

    INCLUDE 'bethy_ctl.inc'

    !! Read namelist bethy_ctl

    IF (p_parallel_io) THEN

       ! define default values
       ncanopy = 3

       f_unit = position_nml ('BETHY_CTL', nml_unit, status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (f_unit, bethy_ctl)
          CALL message('config_bethy', 'Namelist BETHY_CTL: ')
          WRITE(nout, bethy_ctl)
       END SELECT
    ENDIF
    IF (p_parallel_io) THEN
       IF(debug) THEN
          CALL message('config_bethy', 'number of canopy layers: '//int2string(ncanopy))
       END IF
    END IF

    IF (p_parallel) THEN
       CALL p_bcast(ncanopy, p_io)
    ENDIF

    module_configured = .TRUE.

  END SUBROUTINE config_bethy
  !
  !=================================================================================================
  SUBROUTINE bethy_init_io

    USE mo_netCDF, ONLY: add_dim

    REAL(dp), ALLOCATABLE :: canopy_values(:)
    INTEGER :: i

    ALLOCATE (canopy_values(ncanopy))
    DO i=1,ncanopy
       canopy_values(i) = REAL(i,dp)
    END DO

    IF (debug) CALL message('bethy_init_io','Adding dimensions')
    !CALL add_dim("canopy", ncanopy, longname='layers in canopy', units = 'cm', value=(/1._dp/), levtyp=73, indx=CANOPY)
    CALL add_dim("canopy_layer", ncanopy, longname='layers in canopy', units="", value=canopy_values(:), levtyp=73, indx=CANOPY)

    DEALLOCATE(canopy_values)

    IF (debug) CALL message('bethy_init_io','Exit')

  END SUBROUTINE bethy_init_io
  !
  !=================================================================================================
  SUBROUTINE bethy_init_memory(g_nland, l_nland, ntiles, bethy, fileformat, fileztype, diag_stream, stream)

    USE mo_jsbach,      ONLY : missing_value
    USE mo_linked_list, ONLY : NETCDF, LAND, TILES
    USE mo_memory_base, ONLY : new_stream, default_stream_setting, &
                               add =>add_stream_element
    USE mo_netCDF,      ONLY : max_dim_name
    USE mo_output,      ONLY : land_table

    INTEGER,          INTENT(in)     :: g_nland, l_nland, ntiles
    TYPE(bethy_type), INTENT(inout)  :: bethy
    INTEGER,          INTENT(in)     :: fileformat    ! output file format
    INTEGER,          INTENT(in)     :: fileztype     ! output file compression
    TYPE(t_stream), POINTER           :: diag_stream
    TYPE(t_stream), POINTER, OPTIONAL :: stream

    INTEGER                     :: dim1p(1), dim1(1), dim2p(2), dim2(2), dim3p(3), dim3(3)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim2n(2), dim3n(3)

    INTEGER :: high_accuracy

    IF (debug) CALL message('bethy_init_memory','Allocating memory for bethy streams')

    IF (ASSOCIATED(diag_stream)) THEN
       IO_diag => diag_stream
    ELSE
       CALL finish('bethy_init_memory', 'Diagnostic stream not present')
    END IF

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'bethy', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_bethy => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_bethy, 'bethy', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(IO_bethy, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    IF (stream%filetype==NETCDF) THEN
      high_accuracy=64
    ELSE
      high_accuracy=24
    END IF

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n(1) = 'landpoint'

    dim2p = (/ l_nland, ntiles /)
    dim2  = (/ g_nland, ntiles /)
    dim2n(1) = 'landpoint' ; dim2n(2) = 'tiles'

    dim3p = (/ l_nland, ncanopy, ntiles /)
    dim3  = (/ g_nland, ncanopy, ntiles /)
    dim3n(1) = 'landpoint' ; dim3n(2) = 'canopy_layer'; dim3n(3) = 'tiles'

    ! code numbers of this routine range from 120 to 149, using the GRIB land_table

    CALL add(IO_bethy, 'canopy_conductance',         bethy%canopy_conductance,      longname='Canopy Conductance', &
             units='m/s',                      ldims=dim2p, gdims=dim2, dimnames=dim2n, code=120, laccu=.TRUE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'canopy_conductance_unlimited',bethy%canopy_conductance_unlimited,longname='Unlimited Canopy Conductance', &
             units='m/s',                      ldims=dim2p, gdims=dim2, dimnames=dim2n, code=145, laccu=.TRUE., lpost=.FALSE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'CO2_conc_leaf',              bethy%CO2_conc_leaf,           longname='CO2 Concentration inside Leaf', &
             units='mol(CO2)/mol(air)',        ldims=dim2p, gdims=dim2, dimnames=dim2n, code=123,               lpost=.FALSE.)
    CALL add(IO_bethy, 'CO2_conc_leaf_acc',          bethy%CO2_conc_leaf_acc, &
             longname='CO2 Concentration inside Leaf (avg.)', &
             units='mol(CO2)/mol(air)',        ldims=dim2p, gdims=dim2, dimnames=dim2n, code=146, laccu=.true., lpost=.FALSE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'CO2_conc_leaf_unlimited_acc',bethy%CO2_conc_leaf_unlimited_acc, &
             longname='CO2 Concentration inside Leaf for no water stress (avg.)', &
             units='mol(CO2)/mol(air)',        ldims=dim2p, gdims=dim2, dimnames=dim2n, code=147, laccu=.true., lpost=.FALSE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'net_assimilation_acc',       bethy%net_assimilation,        longname='Canopy Net Carbon Assimilation', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=124, laccu=.TRUE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'fapar',                      bethy%faPAR,                   longname='Fraction of PAR absorbed by Canopy', &
             units='',                         ldims=dim2p, gdims=dim2, dimnames=dim2n, code=125, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_bethy, 'apar_acc',                   bethy%apar_acc,                longname='PAR absorbed by Canopy (avg.)', &
             units='mol(PHOTONS) m-2(canopy) s-1',ldims=dim2p, gdims=dim2, dimnames=dim2n, code=148, laccu=.TRUE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'par_acc',             bethy%par_acc, longname='incoming photosynthetically active radiation (PAR) (avg.)', &
             units='mol(PHOTONS) m-2 s-1',     ldims=dim2p, gdims=dim2, dimnames=dim2n, code=149, lpost=.FALSE., laccu=.TRUE., &
             lmiss=.TRUE., missval=missing_value, contnorest=.TRUE.)
    CALL add(IO_bethy, 'gross_assimilation_acc',     bethy%gross_assimilation_acc,  longname='Gross Assimilation (avg.)', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=126, laccu=.TRUE., &
             lmiss=.TRUE., missval=missing_value, bits=high_accuracy)
    CALL add(IO_bethy, 'dark_respiration_acc',       bethy%dark_respiration_acc,    longname='Dark Respiration of Leaf (avg.)', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=127, laccu=.TRUE., lpost=.FALSE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'max_carbox_rate'   ,         bethy%max_carbox_rate,         longname='Maximum Carboxylation Rate (avg.)', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=129, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_bethy, 'max_e_transport_rate',       bethy%max_e_transport_rate,    longname='Maximum Electron Transport Rate', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=130, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_bethy, 'carbox_rate',                bethy%carbox_rate,             longname='Carboxylation Rate (avg.)', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=131, laccu=.TRUE., lpost=.FALSE.)
    CALL add(IO_bethy, 'e_transport_rate',           bethy%e_transport_rate,        longname='Electron Transport Rate (avg.)', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=132, laccu=.TRUE., lpost=.FALSE., &
             bits=high_accuracy)
    CALL add(IO_bethy, 'gross_assimilation_inst',    bethy%gross_assimilation,      longname='Gross Assimilation (inst.)', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=133, bits=high_accuracy, lpost=.FALSE.)
    CALL add(IO_bethy, 'dark_respiration_inst',      bethy%dark_respiration,        longname='Dark Respiration of Leaf (inst.)', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=134,               lpost=.FALSE.)

    CALL add(IO_bethy, 'CO2_flux_net',               bethy%CO2_flux_net_acc, &
             longname='Net CO2 flux to atmosphere (acc.)',  units='mol(CO2) m-2(grid box) s-1', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=160, laccu=.TRUE., lpost=.TRUE., lrerun=.TRUE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'CO2_flux_herbivory',         bethy%CO2_flux_herbivory_acc, &
             longname='Grazing CO2 flux to atmosphere (acc.)',  units='mol(CO2) m-2(grid box) s-1', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=161, laccu=.TRUE., lpost=.TRUE., lrerun=.TRUE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'CO2_emission_landcover_change', bethy%CO2_emission_landcover_change_acc,  &
             longname='Land-cover change CO2 flux to the atmosphere (acc.)',  units='mol(CO2) m-2(grid box) s-1', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=162, laccu=.TRUE., lpost=.TRUE., lrerun=.TRUE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'CO2_emission_harvest',       bethy%CO2_emission_harvest_acc, &
             longname='Harvest CO2 flux to atmosphere (acc.)',  units='mol(CO2) m-2(grid box) s-1', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=163, laccu=.TRUE., lpost=.TRUE., lrerun=.TRUE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_bethy, 'CO2_flux_dynveg',            bethy%CO2_flux_dynveg_acc, &
             longname='Disturbance CO2 flux to atmosphere (acc.)',  units='mol(CO2) m-2(grid box) s-1', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=164, laccu=.TRUE., lpost=.TRUE., lrerun=.TRUE., &
             lmiss=.TRUE., missval=missing_value)

    CALL add(IO_bethy, 'lai_layer',                  bethy%Canopy%lai,                  dim3p, dim3, dimnames=dim3n, code=135, &
             lpost=.FALSE., leveltype=CANOPY)
    CALL add(IO_bethy, 'canopy_conductance_layer',   bethy%Canopy%canopy_conductance,   dim3p, dim3, dimnames=dim3n, code=136, &
             lpost=.FALSE., leveltype=CANOPY)
    CALL add(IO_bethy, 'gross_assimiliation_layer',  bethy%Canopy%gross_assimilation,   dim3p, dim3, dimnames=dim3n, code=137, &
             bits=high_accuracy, lpost=.FALSE., leveltype=CANOPY)
    CALL add(IO_bethy, 'dark_respiration_layer',     bethy%Canopy%dark_respiration,     dim3p, dim3, dimnames=dim3n, code=138, &
             lpost=.FALSE., leveltype=CANOPY)
    CALL add(IO_bethy, 'max_carbox_rate_layer',      bethy%Canopy%max_carbox_rate,      dim3p, dim3, dimnames=dim3n, code=140, &
             lpost=.FALSE., leveltype=CANOPY)
    CALL add(IO_bethy, 'max_e_transport_rate_layer', bethy%Canopy%max_e_transport_rate, dim3p, dim3, dimnames=dim3n, code=141, &
             lpost=.FALSE., leveltype=CANOPY)
    CALL add(IO_bethy, 'carbox_rate_layer',          bethy%Canopy%carbox_rate,          dim3p, dim3, dimnames=dim3n, code=142, &
             lpost=.FALSE., leveltype=CANOPY)
    CALL add(IO_bethy, 'e_transport_rate_layer',     bethy%Canopy%e_transport_rate,     dim3p, dim3, dimnames=dim3n, code=143, &
             lpost=.FALSE., leveltype=CANOPY)
    CALL add(IO_bethy, 'faPAR_layer',                bethy%Canopy%faPAR,                dim3p, dim3, dimnames=dim3n, code=144, &
             lpost=.FALSE., leveltype=CANOPY)

    CALL add(IO_diag, 'net_assimilation',      bethy_diag%net_assimilation,      longname='Canopy Net Carbon Assimilation (avg.)', &
             units='mol(CO2) m-2(grid box) s-1', ldims=dim1p,gdims=dim1, dimnames=dim1n,   code=124, &
             laccu=.TRUE.,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'fapar',                 bethy_diag%faPAR,             longname='Fraction of PAR Absorbed by Canopy (avg.)', &
             units='',                         ldims=dim1p, gdims=dim1, dimnames=dim1n,    code=125, lpost=.FALSE., &
             laccu=.TRUE.,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'apar_acc',              bethy_diag%apar_acc,              longname='PAR Absorbed by Canopy (avg.)', &
             units='mol(PHOTONS) m-2(grid box) s-1', ldims=dim1p, gdims=dim1, dimnames=dim1n,code=148, lpost=.TRUE., &
             laccu=.TRUE.,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'par_acc',         bethy_diag%par_acc, longname='incoming photosynthetically active radiation (PAR) (avg.)', &
             units='mol(PHOTONS) m-2 s-1',       ldims=dim1p, gdims=dim1, dimnames=dim1n,  code=149, lpost=.TRUE., &
             laccu=.TRUE.,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'gross_assimilation',    bethy_diag%gross_assimilation,    longname='Gross Assimilation (avg.)', &
             units='mol(CO2) m-2(grid box) s-1', ldims=dim1p, gdims=dim1, dimnames=dim1n,  code=126, &
             laccu=.TRUE.,  lmiss=.TRUE., missval=missing_value)
    CALL add(IO_diag, 'dark_respiration',      bethy_diag%dark_respiration,      longname='Dark Respiration of Leaf (avg.)', &
             units='mol(CO2) m-2(grid box) s-1', ldims=dim1p, gdims=dim1, dimnames=dim1n,  code=127, lpost=.FALSE., &
             laccu=.TRUE.,  lmiss=.TRUE., missval=missing_value)

  END SUBROUTINE bethy_init_memory

  !=========================================================================================================

  SUBROUTINE init_bethy(grid, domain, bethy, fileformat, fileztype, IO_diag_stream, IO_stream)

    TYPE(grid_type),   INTENT(in)     :: grid
    TYPE(domain_type), INTENT(in)     :: domain
    TYPE(bethy_type),  INTENT(inout)  :: bethy
    INTEGER,           INTENT(in)     :: fileformat
    INTEGER,           INTENT(in)     :: fileztype
    TYPE(t_stream), POINTER           :: IO_diag_stream
    TYPE(t_stream), POINTER, OPTIONAL :: IO_stream

    INTEGER :: ntiles
    INTEGER :: i

    CALL message('init_bethy','Start initialization of BETHY')
    IF (debug) CALL message('init_bethy','Start initialization of BETHY')

    IF (.NOT. module_configured) CALL config_bethy

    ! Save number of tiles
    ntiles = bethy%ntiles
    bethy%Canopy%ntiles = ntiles

    CALL bethy_init_io

    ! Save number of canopy layers in structure; ncanopy is read from run.def in bethy_config
    bethy%Canopy%ncanopy = ncanopy

    ! initialization of memory and generation of streams:
    CALL bethy_init_memory(grid%nland, domain%nland, ntiles, bethy, fileformat, fileztype, IO_diag_stream, stream=IO_stream)

    !------------------------------------------------------------------------------

    ! compute borders of the canopy layers in terms of LAI 

    ALLOCATE(canopy_boundaries_lai(0:ncanopy))
    DO i=0,ncanopy
       canopy_boundaries_lai(i) = REAL(i,dp) / REAL(ncanopy,dp)
    END DO

    bethy%e_transport_rate = 0._dp
    bethy%canopy_conductance = 0._dp
    bethy%canopy_conductance_unlimited = 0._dp
    bethy%CO2_conc_leaf_acc = 0._dp
    bethy%CO2_conc_leaf_unlimited_acc = 0._dp
    bethy%net_assimilation = 0._dp
    bethy%apar_acc = 0._dp
    bethy%par_acc = 0._dp
    bethy%gross_assimilation_acc = 0._dp
    bethy%dark_respiration_acc = 0._dp
    bethy%CO2_flux_net_acc = 0._dp
    bethy%CO2_flux_herbivory_acc = 0._dp
    bethy%CO2_emission_landcover_change_acc = 0._dp
    bethy%CO2_emission_harvest_acc = 0._dp
    bethy%CO2_flux_dynveg_acc = 0._dp

    bethy_diag%net_assimilation = 0._dp
    bethy_diag%faPAR = 0._dp
    bethy_diag%apar_acc = 0._dp
    bethy_diag%par_acc = 0._dp
    bethy_diag%gross_assimilation = 0._dp
    bethy_diag%dark_respiration = 0._dp

    IF (debug) CALL message('init_bethy','Initialization of BETHY finished')

  END SUBROUTINE init_bethy

  SUBROUTINE bethy_diagnostics(surface, bethy, veg_fract_correction)
    
    USE mo_time_control, ONLY: delta_time
    USE mo_land_surface, ONLY: land_surface_type
    USE mo_utils,        ONLY: average_tiles

    TYPE(land_surface_type), INTENT(in) :: surface
    REAL(dp),                INTENT(in) :: veg_fract_correction(:,:)   ! factor representing clumpyness of vegetation
    TYPE(bethy_type), INTENT(inout) :: bethy

    !LOCAL VARIABLES
    INTEGER  :: icanopy, ncanopy
    INTEGER  :: itile, ntiles
    INTEGER  :: kidx0, kidx1
    REAL(dp) :: veg_fract(nidx,bethy%ntiles)

    ncanopy = bethy%Canopy%ncanopy
    ntiles  = bethy%ntiles
    kidx0 = kstart
    kidx1 = kend

    ! Integrate canopy variables over canopy layers
    ! Note that these accumulated variables should have their stream definition with laccu=.TRUE.
    ! At each output interval, the accumulated values are then divided by the number of seconds in the output interval

    DO itile=1,ntiles
    ! These two were already integrated in bethy from instantaneous variables
    bethy%gross_assimilation_acc(kidx0:kidx1,itile)   = bethy%gross_assimilation_acc(kidx0:kidx1,itile)   &
         + bethy%gross_assimilation(kidx0:kidx1,itile) * delta_time
    bethy%dark_respiration_acc(kidx0:kidx1,itile)     = bethy%dark_respiration_acc(kidx0:kidx1,itile)       &
         + bethy%dark_respiration(kidx0:kidx1,itile) * delta_time
    ! Net assimilation = gross - dark
    bethy%net_assimilation(kidx0:kidx1,itile) = bethy%gross_assimilation_acc(kidx0:kidx1,itile) &
         - bethy%dark_respiration_acc(kidx0:kidx1,itile)

    DO icanopy = 1,ncanopy
       bethy%max_carbox_rate(kidx0:kidx1,itile)      = bethy%max_carbox_rate (kidx0:kidx1,itile)       &
            + bethy%Canopy%max_carbox_rate(kidx0:kidx1,icanopy,itile) &
            * bethy%Canopy%lai(kidx0:kidx1,icanopy,itile)* delta_time
       bethy%max_e_transport_rate(kidx0:kidx1,itile) = bethy%max_e_transport_rate (kidx0:kidx1,itile)  &
            + bethy%Canopy%max_e_transport_rate(kidx0:kidx1,icanopy,itile) &
            * bethy%Canopy%lai(kidx0:kidx1,icanopy,itile)* delta_time
       bethy%carbox_rate (kidx0:kidx1,itile)         = bethy%carbox_rate (kidx0:kidx1,itile)           &
            + bethy%Canopy%carbox_rate(kidx0:kidx1,icanopy,itile) &
            * bethy%Canopy%lai(kidx0:kidx1,icanopy,itile)* delta_time
       bethy%e_transport_rate(kidx0:kidx1,itile)     = bethy%e_transport_rate (kidx0:kidx1,itile)      &
            + bethy%Canopy%e_transport_rate(kidx0:kidx1,icanopy,itile) &
            * bethy%Canopy%lai(kidx0:kidx1,icanopy,itile)* delta_time
       bethy%faPAR(kidx0:kidx1,itile)                = bethy%faPAR (kidx0:kidx1,itile)                 &
            + bethy%Canopy%faPAR(kidx0:kidx1,icanopy,itile) * delta_time
    END DO

    END DO

    veg_fract(:,:) = SPREAD(surface%veg_ratio_max(kidx0:kidx1), NCOPIES=ntiles, DIM=2) * veg_fract_correction(:,:)
    CALL average_tiles(bethy%gross_assimilation_acc(kidx0:kidx1,1:ntiles) * veg_fract(:,:), &
             surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
             bethy_diag%gross_assimilation(kidx0:kidx1))
    CALL average_tiles(bethy%dark_respiration_acc(kidx0:kidx1,1:ntiles) * veg_fract(:,:), &
             surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
             bethy_diag%dark_respiration(kidx0:kidx1))
    CALL average_tiles(bethy%faPAR(kidx0:kidx1,1:ntiles) * veg_fract(:,:), &
             surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
             bethy_diag%faPAR(kidx0:kidx1))
    CALL average_tiles(bethy%apar_acc(kidx0:kidx1,1:ntiles) * veg_fract(:,:), &
             surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
             bethy_diag%apar_acc(kidx0:kidx1))
    CALL average_tiles(bethy%par_acc(kidx0:kidx1,1:ntiles), &
             surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
             bethy_diag%par_acc(kidx0:kidx1))
    CALL average_tiles(bethy%net_assimilation(kidx0:kidx1,1:ntiles) * veg_fract(:,:), &
             surface%is_present(kidx0:kidx1,1:ntiles), surface%cover_fract(kidx0:kidx1,1:ntiles), &
             bethy_diag%net_assimilation(kidx0:kidx1))

  END SUBROUTINE bethy_diagnostics
  
  ! --- bethy() --------------------------------------------------------------------------------------------------------------------
  !
  ! This is the main routine of BETHY
  !
  
SUBROUTINE update_bethy(kidx, domain, mask, bethy, cover_type, lctlib, &
     UseAlbedo, waterLimitationFlag, lai, &
     par, frac_par_direct, pressure, &
     canopy_temp, soil_albedo, CO2_concentration_air, &
     canopy_conductance, canopy_conductance_limited &
     )

    USE mo_time_control,     ONLY: delta_time
    USE mo_jsbach_constants, ONLY: molarMassCO2_kg, molarMassDryAir_kg, Tmelt
    USE mo_bethy_constants,  ONLY: EPar, SoilReflectivityParMin,FCI1C3,FCI1C4
    USE mo_bethy_fapar,      ONLY: faparl
    USE mo_zenith,           ONLY: cos_zenith

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kidx
    TYPE(domain_type), INTENT(in) :: domain
    LOGICAL, INTENT(in)           :: mask(:,:)
    TYPE(bethy_type),INTENT(inout) :: bethy
    INTEGER, INTENT(in)           :: cover_type(1:kidx,1:bethy%ntiles)  ! land cover type of <ntiles> tiles
    TYPE(lctlib_type), INTENT(in) :: lctlib
    LOGICAL,     INTENT(in)    :: UseAlbedo                      ! Indicates which albedo scheme is used ...
                                                                 ! ... (interactive albedo (true) or ECHAM scheme (false)
    LOGICAL,     INTENT(in)    :: waterLimitationFlag            ! Chooses between computation of water limited (true) and ...
                                                                 ! ... water unlimited (flase) photosynthesis
    REAL(dp),    INTENT(in)    :: lai                  (1:kidx,1:bethy%ntiles)  ! Leaf Area Index [m^2/m^2]
    REAL(dp),    INTENT(in)    :: par                  (1:kidx)  ! Photosynth. active radiation (400-700nm) [W/m^2] 
    REAL(dp),    INTENT(in)    :: frac_par_direct      (1:kidx)  ! fraction of direct radiation in par  
    REAL(dp),    INTENT(in)    :: pressure             (1:kidx)  ! Surface pressure
    REAL(dp),    INTENT(in)    :: canopy_temp          (1:kidx)  ! Temperature of canopy/vegetation
    REAL(dp),    INTENT(in)    :: soil_albedo          (1:kidx,1:bethy%ntiles)  ! Soil reflectance
    REAL(dp),    INTENT(in)    :: CO2_concentration_air(1:kidx)  ! CO2 concentration of ambient air 
                                                                 ! (mass mixing ratio) [kg(CO2)/kg(air)]
    REAL(dp),    INTENT(out)   :: canopy_conductance   (1:kidx,1:bethy%ntiles)  ! Canopy conductance no limitation
    REAL(dp),    INTENT(in)    :: canopy_conductance_limited(1:kidx,1:bethy%ntiles)  ! Canopy conductance water limit

    ! Local variables per LCT (used inside the LCT loop)
    REAL(dp):: nitrogen_scaling_factors (kidx,ncanopy)
    REAL(dp):: apar_layer               (kidx,ncanopy)   ! Absorbed PAR per canopy layer ...
                                                         ! ... [mol(absorbed photons)/(m^2(leaf area) s)]
    REAL(dp):: pardown_mol              (kidx)           ! Downward PAR flux in mol (photons)/(m^2 s)
    REAL(dp):: soil_reflectivity_par    (kidx)           ! Soil reflectivity in the PAR region
    INTEGER :: lct                      (kidx)
    REAL(dp):: CO2_concentration_mol    (kidx)
    INTEGER :: ntiles, itile
    INTEGER :: kidx0, kidx1, icanopy

    kidx0 = kstart
    kidx1 = kend

    ntiles = SIZE(cover_type, DIM=2)

    ! Initialize local variables
    soil_reflectivity_par(:) = 0._dp
    pardown_mol(:) = 0._dp

    ! Expresse radiation in mol(photons) / (m^2 s)
    pardown_mol(:) = par(:) / EPar

    ! Accumulate incoming par in mol(photons) / (m^2 s) (only once per timestep!)
    IF (.NOT. waterLimitationFlag) THEN
       bethy%par_acc(kidx0:kidx1,1:ntiles) = bethy%par_acc(kidx0:kidx1,1:ntiles) + &
            SPREAD(par(1:nidx), NCOPIES=ntiles, DIM=2) / Epar * delta_time
    END IF

    ! Loop over natural vegetation landcover types
    tile_loop: DO itile=1,ntiles
       lct(:) = cover_type(:,itile)  ! Vector of land cover types for this tile

       !if(.not. lctlib%naturalVegFlag(lct(:))) cycle !! Continue with next landcover type when not natural vegetation

       ! Initilize local variables
       nitrogen_scaling_factors = 1._dp
       apar_layer = 0._dp

       ! Compute soil reflectivity in PAR region
       IF (UseAlbedo) THEN
         ! soil reflectivity is set to soil albedo of the visible range
         soil_reflectivity_par(:) = soil_albedo(:,itile)
       ELSE
         ! soil reflectivity is derived from background albedo of the whole solar spectrum (compare eq. (122) in Knorr)
         soil_reflectivity_par(:) = MAX(0.92_dp * soil_albedo(:,itile) - 0.015_dp, SoilReflectivityParMin)
       END IF

       ! calculate nitrogen scaling factors 
       CALL calc_nitrogen_scaling_factors(lai(1:nidx,itile),                                     &
                                          canopy_boundaries_lai(0:ncanopy),                      &
                                          nitrogen_scaling_factors(1:nidx,1:ncanopy))

       DO icanopy=1,ncanopy
          WHERE(.NOT. mask(:,itile) .OR. .NOT. lctlib%NitrogenScalingFlag(lct(:)))
             nitrogen_scaling_factors(:,icanopy) = 1._dp ! to some vegetation classes nitrogen scaling is not applied
          END WHERE
       END DO

       ! Compute the absorbed PAR per leaf area index and per canopy layer with two flux approximation, normalized to
       ! incoming radiation and calculate LAI per canopy layer
       CALL faparl(&
            nidx,                            &  ! input: number of grid points to be handled in this call
            mask(:,itile),                   &
            ncanopy,                         &  ! input: number of canopy layers
            lai(1:nidx,itile),               &  ! input: leaf area index
            soil_reflectivity_par(1:nidx),   &  ! input:
            cos_zenith(kidx0:kidx1),              &  ! input:
            frac_par_direct(1:nidx),         &  ! input:
            canopy_boundaries_lai(0:ncanopy),&  ! input:
            bethy%Canopy%lai(kidx0:kidx1,1:ncanopy,itile),  & ! output: LAI per canopy layer
            bethy%Canopy%faPAR(kidx0:kidx1,1:ncanopy,itile) & ! output: normalized absorbed PAR per leaf area and canopy layer
            )

       ! Compute absorbed PAR per leaf area in canopy layer [units: (absorbed photons) / (m^2(leaf area) s)] from 
       ! par and fraction of absorbed PAR (Epar is needed to convert radiation intensity from W/m^2 to mol/(m^2 s))
       DO icanopy=1,ncanopy
             apar_layer(1:nidx,icanopy) = (par(1:nidx) / Epar) * bethy%Canopy%faPAR(kidx0:kidx1,icanopy,itile) / &
                  (MAX(bethy%Canopy%lai(kidx0:kidx1,icanopy,itile),1.e-10_dp))
             apar_layer(1:nidx,icanopy) = MERGE(apar_layer(1:nidx,icanopy),0._dp,mask(1:nidx,itile))
             IF (waterLimitationFlag) &  ! Only accumulate once per timestep
                  bethy%apar_acc(kidx0:kidx1,itile) = bethy%apar_acc(kidx0:kidx1,itile) + &
                                   (par(1:nidx) / Epar) * bethy%Canopy%faPAR(kidx0:kidx1,icanopy,itile) * delta_time
       END DO

       ! Convert CO2 mass mixing ratio [kg/kg] to volume mixing ratio [mol/mol]
       CO2_concentration_mol = CO2_concentration_air * molarMassDryAir_kg / molarMassCO2_kg

       ! Calculate stressed or unstressed photosynthesis, i.e. with or without water limitation

       IF (.NOT. waterLimitationFlag) THEN ! calculation WITHOUT water limitation
                                           ! input  to photosynthesis routine: leaf internal CO2-concentration
                                           ! output of photosynthesis routine: stomata conductance of canopy layer
 
          ! --- estimate leaf internal CO2-concentration

          WHERE(lctlib%C4flag(lct))               ! landcover type is C4 plants
             bethy%CO2_conc_leaf(kidx0:kidx1,itile) = FCI1C4*CO2_concentration_mol(1:nidx)
          ELSEWHERE                               ! landcover type is C3 plants
             bethy%CO2_conc_leaf(kidx0:kidx1,itile) = FCI1C3*CO2_concentration_mol(1:nidx)
          END WHERE

          ! Accumulate CO2 concentration inside the leaves for unstresssed photosynthesis
          bethy%CO2_conc_leaf_unlimited_acc(kidx0:kidx1,itile) = bethy%CO2_conc_leaf_unlimited_acc(kidx0:kidx1,itile) + &
                                                                 bethy%CO2_conc_leaf(kidx0:kidx1,itile) * delta_time

          DO icanopy=1,ncanopy ! go through all canopy layers 
             CALL photosyn(kidx, &
                  mask(1:nidx,itile) .AND.     &
                   ( lctlib%NaturalVegFlag(lct(1:nidx)) .OR. lctlib%CropFlag(lct(1:nidx)) .OR. lctlib%PastureFlag(lct(1:nidx)) ),&
                  waterLimitationFlag,                                  & ! .false. in this case
                  lctlib%CarboxRate               (lct(1:nidx)),        & ! parameters from lctlib
                  lctlib%ETransport               (lct(1:nidx)),        &
                  lctlib%C4flag                   (lct(1:nidx)),        & 
                  apar_layer                      (1:nidx,icanopy),     & ! input: absorbed par per LAI [mol phot./(m^2 s)]
                  (canopy_temp                    (1:nidx) - Tmelt),    & 
                  pressure                        (1:nidx),             &
                  nitrogen_scaling_factors        (1:nidx,icanopy),     &
                  pardown_mol                     (1:nidx),             & ! input: par downward [mol phot./(m^2 s)]
                  CO2_concentration_mol           (1:nidx),             & ! atm_co2_conc for photosynthesis, co2 from atmosphere
                  bethy%Canopy%canopy_conductance (kidx0:kidx1,icanopy,itile), & ! output (needed for second call of photosynthesis)
                  bethy%CO2_conc_leaf             (kidx0:kidx1,itile)   & ! input (identical for all layers)
                  )

          END DO

       ELSE ! With water limitation
          ! Compute unlimited conductance for the whole canopy
          canopy_conductance(:,itile) = MAX(1.e-20_dp,SUM(bethy%Canopy%canopy_conductance(kidx0:kidx1,:,itile) * &
               bethy%Canopy%lai(kidx0:kidx1,:,itile),DIM=2))
          ! Scale unlimited conductance per leaf area with the ratio of limited to unlimited conductance for
          ! whole canopy, ie. compute new limited conductance per leaf area and canopy layer
          DO icanopy=1,ncanopy
             bethy%Canopy%canopy_conductance(kidx0:kidx1,icanopy,itile) = &
                  bethy%Canopy%canopy_conductance(kidx0:kidx1,icanopy,itile) * canopy_conductance_limited(:,itile) / &
                  canopy_conductance(:,itile)
             bethy%Canopy%canopy_conductance(kidx0:kidx1,icanopy,itile) = &
                  MERGE(bethy%Canopy%canopy_conductance(kidx0:kidx1,icanopy,itile),0._dp,mask(1:nidx,itile))
          END DO
          DO icanopy=1,ncanopy ! go through all canopy layers 
             CALL photosyn(kidx, &
                  mask(1:nidx,itile) .AND. & 
                   ( lctlib%NaturalVegFlag(lct(1:nidx)) .OR. lctlib%CropFlag(lct(1:nidx)) .OR. lctlib%PastureFlag(lct(1:nidx)) ),&
                  waterLimitationFlag,                               & ! .true. in this case
                  lctlib%CarboxRate(lct(1:nidx)),                    &
                  lctlib%ETransport(lct(1:nidx)),                    &
                  lctlib%C4flag(lct(1:nidx)),                        & ! parameters from lctlib
                  apar_layer              (1:nidx,icanopy),          & ! input: absorbed par per LAI [mol phot./(m^2 s)]
                  (canopy_temp            (1:nidx) - Tmelt),         &
                  pressure                (1:nidx),                  &
                  nitrogen_scaling_factors(1:nidx,icanopy),          &
                  pardown_mol             (1:nidx),                  & ! input: par downward [mol phot./(m^2 s)]
                  CO2_concentration_mol            (1:nidx),         & ! atm_co2_conc for photosynthesis, co2 from atmosphere
                  bethy%Canopy%canopy_conductance  (kidx0:kidx1,icanopy,itile), & ! input from water-unlimited case
                  bethy%CO2_conc_leaf              (kidx0:kidx1,itile),         & ! output
                  bethy%Canopy%gross_assimilation  (kidx0:kidx1,icanopy,itile), & ! output
                  bethy%Canopy%dark_respiration    (kidx0:kidx1,icanopy,itile), & ! output 
                  bethy%Canopy%max_carbox_rate     (kidx0:kidx1,icanopy,itile), & ! output
                  bethy%Canopy%max_e_transport_rate(kidx0:kidx1,icanopy,itile), & ! output
                  bethy%Canopy%carbox_rate         (kidx0:kidx1,icanopy,itile), & ! output
                  bethy%Canopy%e_transport_rate    (kidx0:kidx1,icanopy,itile)  & ! output
                  )

             ! Accumulate CO2 concentration inside the leaves for stressed photosynthesis (all canopy layers weighted equally)
             bethy%CO2_conc_leaf_acc(kidx0:kidx1,itile) = bethy%CO2_conc_leaf_acc(kidx0:kidx1,itile) + &
                                                          bethy%CO2_conc_leaf(kidx0:kidx1,itile) * delta_time / REAL(ncanopy,dp)
          END DO

       END IF

    END DO tile_loop

    ! Compute canopy conductance, gross and dark respiration by integrating over canopy layers
    canopy_conductance(1:nidx,1:ntiles) = &
         MAX(1.e-20_dp, SUM(bethy%Canopy%canopy_conductance(kidx0:kidx1,:,1:ntiles) * bethy%Canopy%lai(kidx0:kidx1,:,1:ntiles),  &
                         DIM=2))

#ifndef __PGI
    bethy%gross_assimilation(kidx0:kidx1,1:ntiles) = &
         SUM(bethy%Canopy%gross_assimilation(kidx0:kidx1,:,1:ntiles) * bethy%Canopy%lai(kidx0:kidx1,:,1:ntiles), DIM=2)
    bethy%dark_respiration(kidx0:kidx1,1:ntiles) = &
         SUM(bethy%Canopy%dark_respiration(kidx0:kidx1,:,1:ntiles) * bethy%Canopy%lai(kidx0:kidx1,:,1:ntiles), DIM=2)
#else
    !
    ! Different notation to allow compilation with the PGI compiler (pgf95 6.1-1)
    !
    bethy%gross_assimilation(kidx0:kidx1,1:ntiles) = 0._dp
    DO icanopy=1,ncanopy
       bethy%gross_assimilation(kidx0:kidx1,1:ntiles) = bethy%gross_assimilation(kidx0:kidx1,1:ntiles) &
            + bethy%Canopy%gross_assimilation(kidx0:kidx1,icanopy,1:ntiles) * bethy%Canopy%lai(kidx0:kidx1,icanopy,1:ntiles)
    END DO
    bethy%dark_respiration(kidx0:kidx1,1:ntiles) = 0._dp
    DO icanopy=1,ncanopy
       bethy%dark_respiration(kidx0:kidx1,1:ntiles) =  bethy%dark_respiration(kidx0:kidx1,1:ntiles) &
            + bethy%Canopy%dark_respiration(kidx0:kidx1,icanopy,1:ntiles) * bethy%Canopy%lai(kidx0:kidx1,icanopy,1:ntiles)
    END DO
#endif

    DO itile=1,ntiles
       WHERE (.NOT. mask(:,itile))
          canopy_conductance(:,itile) = 0._dp
          bethy%gross_assimilation(kidx0:kidx1,itile) = 0._dp
          bethy%dark_respiration(kidx0:kidx1,itile) = 0._dp
       END WHERE
    END DO

    ! Accumulate unlimited or limited canopy conductance (laccu=.TRUE. for these variables)
    IF (.NOT. waterLimitationFlag) THEN
       bethy%canopy_conductance_unlimited(kidx0:kidx1,1:ntiles) = bethy%canopy_conductance_unlimited(kidx0:kidx1,1:ntiles) + &
                                                                  canopy_conductance(1:nidx,1:ntiles) * delta_time
    ELSE
       bethy%canopy_conductance(kidx0:kidx1,1:ntiles) = bethy%canopy_conductance(kidx0:kidx1,1:ntiles) + &
                                                        canopy_conductance(1:nidx,1:ntiles) * delta_time
    END IF

    RETURN

CONTAINS

    SUBROUTINE calc_nitrogen_scaling_factors(zlai,layer_bounds, factors)

      ! Calculates the nitrogen scaling factors for each grid point and canopy layer for a single PFT
      ! Eqs. (106), (107) in Knorr (1997)
      ! Note that this routine should only be applied to those vegetation classes for which nitrogen scaling shall be applied!!

      USE mo_bethy_constants,  ONLY: LaiLimit
      USE mo_orbit,            ONLY: inquire_declination
      USE mo_jsbach,           ONLY: interface_test

      REAL(dp), INTENT(in)  :: zlai(nidx)              ! LAI
      REAL(dp), INTENT(in)  :: layer_bounds(0:ncanopy) ! Borders of canopy layers
      REAL(dp), INTENT(out) :: factors(nidx,ncanopy)   ! Nitrogen scaling factors

      REAL(dp) :: declination                          ! Solar Declination
      REAL(dp) :: k12(nidx), cos_zenith_noon(nidx)
      INTEGER  :: i

      factors = 1._dp

      CALL inquire_declination(declination)

      ! Cosine of zenith angle at local noon: 

      cos_zenith_noon(:) = cos(declination) * domain%coslat(kidx0:kidx1) + sin(declination) * domain%sinlat(kidx0:kidx1)

      ! Extinction factor
      
      k12(:) = 0.5_dp / MAX(cos_zenith_noon(:), 1.E-3_dp) ! K = 1 / 2mue , Knorr (119c)

      !vg: in the current echam version inquire_declination at radiation time steps returns the declination at radiation time
      !    and not the current declination. This leads to slightly different result for offline and online jsbach simulations
      !    (https://code.zmaw.de/issues/4803). For that reason we use a k12 independand of the declination for the interface test. 
      IF (interface_test) k12(:)=1._dp

      ! Condition: LAI>LaiLimit
      DO i=1,ncanopy
!CDIR NOASSOC 
         WHERE(zlai(:) >= LaiLimit)
            factors(:,i) = MAX(1.e-10_dp,EXP( -k12(:) * layer_bounds(i-1) * zlai(:))) !! see Eqs. (107), (108) in Knorr
         END WHERE
      END DO

      RETURN

    END SUBROUTINE calc_nitrogen_scaling_factors

  END SUBROUTINE update_bethy

END MODULE mo_bethy
