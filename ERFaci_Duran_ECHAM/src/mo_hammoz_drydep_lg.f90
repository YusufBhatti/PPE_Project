!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_drydep_lg.f90
!!
!! \brief
!! This module handles all input terms required for the calculations of the surface
!! trace gas and aerosol dry deposition.
!!
!! \author M. Schultz        (FZ Juelich)
!! \author Hans-Stefan Bauer (MPIfM)
!! \author Laurens Ganzeveld (MPIfM)
!! \author Andreas Rhodin    (MPIfM)
!! \author Philip Stier      (MPIfM)
!! \author Grazia Frontoso   (C2SM)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# Martin Schultz (FZ Juelich), Hans-Stefan Bauer (MPIfM) - original code (2000-07)
!!   -# Laurens Ganzeveld (MPIfM), Andreas Rhodin (MPIfM) - revision (2001-10)
!!   -# Philip Stier (MPIfM) - revision (2002-2006)
!!   -# Martin Schultz (FZ Juelich) - harmonized tracer indices for HAMMOZ
!!                                    allocatable arrays (jptrac->ntrac)
!!                                    merge with moz_drydep (2009-01)
!!   -# Martin Schultz (FZ Juelich) - further restructuring (2009-11)
!!   -# Grazia Frontoso (C2SM) - usage of the input variables  defined over land, water, ice
!!                               to account for the non-linearity in the drydep
!!                               calculations for gridboxes containing both water
!!                               and sea ice (see #78 on redmine)
!!   -# Martin Schultz (FZ Juelich) 2014/05/19 - bug fix vdice, re-ordering of vd calculations, removed MIN(vd)
!!
!! \limitations
!! None
!!
!! \details
!! Most of the input parameters are vegetation and soil data derived from satellite
!! data a high-resolution geograhical databases. For more details see the
!! routine where the actual reading occurs. The data are monthly mean
!! values, whenever there is an annual cycle in the data. In this 
!! module, there is also the initialisation of the surface resistances
!! used for the dry deposition calculations.
!!
!! \bibliographic_references
!! None
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_hammoz_drydep_lg
  USE mo_kind,          ONLY: dp

  IMPLICIT NONE
  !----------------
  ! Public entities
  !----------------
  PRIVATE

  ! -- subroutines
  PUBLIC :: drydep_lg_init     ! initialisation 
  PUBLIC :: init_drydep_lg_stream
  PUBLIC :: drydep_lg_calcra   ! calculate aerodynamc resistance
  PUBLIC :: drydep_lg_vdbl     ! calculate dry deposition velocities for gases


  ! -- variables
  PUBLIC :: diff, diffrb, rmes, rcut, rsoil, rws, rwater, rsnow    ! needed for ham_vdaer

  ! mz_LG-20020115 declaration of the resistances, being used in the dry
  !     deposition routines, for details see vdbl.
  !     One scalar value per tracer
 
  REAL(dp), ALLOCATABLE :: diff(:),       &  
                           diffrb(:),     &
                           rmes(:),       &    ! mesophyilic resistance
                           rcut(:),       &    ! cuticular resistance
                           rsoil(:),      &    ! soil resistance
                           rws(:),        &    ! wet skin reservoir resistance
                           rwater(:),     &    ! water surface resistance
                           rsnow(:)            ! snow surface resistance

  INTEGER         :: idt_so2, idt_h2so4, idt_o3, idt_hno3, &
                     idt_no,  idt_no2,                   &
                     idt_o3s                  !SF for moz

  INTEGER         :: ibc_lai, ibc_hc, ibc_z0m, ibc_forest, &   ! boundary condition indices
                     ibc_soilph(7)

  LOGICAL         :: lhas_sulfur               ! flag to indicate if sulphate calculations are necessary

!++mgs 20140519: additional pointers for detailed diagnostics
  REAL(dp), POINTER, PUBLIC  :: dpzlai(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpzrmes(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpzhc(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpz0m(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprahcan(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpzrco_leaf(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dppustveg(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dppustslsn(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dppustarw(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dppustari(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprleaf(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprstom(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprsveg(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprbslsn(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprws(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprbveg(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dppcvs(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dppcvw(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dppvgrat(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpprahveg(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpprahslsn(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpvdveg(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpvdsoil(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpvdws(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpvdsn(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpvsn(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpvdland(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpvdwat(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpvdice(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpprahwat(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprbw(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dpprahice(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprbi(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprc0x(:,:) => NULL()
  REAL(dp), POINTER, PUBLIC  :: dprsnowhno3so2(:,:) => NULL()
!--mgs

  CONTAINS

 
 
  !! ------------------------------------------------------------------------------------------
  !! initialize Laurens Ganzeveld's drydep scheme
  !!

  SUBROUTINE drydep_lg_init

    USE mo_tracdef,            ONLY: ntrac
    USE mo_tracer,             ONLY: get_tracer
    USE mo_boundary_condition, ONLY: bc_nml, bc_define,          &
                                     BC_REPLACE, BC_BOTTOM
    USE mo_external_field_processor, ONLY: EF_FILE, EF_MODULE,   &
                                           EF_LONLAT,            &
                                           EF_IGNOREYEAR, EF_CONSTANT, &
                                           EF_NOINTER
    USE mo_ham,                ONLY: nlai_drydep_ef_type !gf #244
    USE mo_exception,          ONLY: message, em_error

    INTEGER                     :: j, ierr
    CHARACTER(len=7), PARAMETER :: cnum = '1234567'
    TYPE(bc_nml)                :: bc_lai, bc_hc, bc_z0m, bc_forest, bc_soilph

    ! allocate arrays
    IF (.NOT. ALLOCATED(diff))    ALLOCATE (diff(ntrac))
    IF (.NOT. ALLOCATED(diffrb))  ALLOCATE (diffrb(ntrac))
    IF (.NOT. ALLOCATED(rmes))    ALLOCATE (rmes(ntrac))
    IF (.NOT. ALLOCATED(rcut))    ALLOCATE (rcut(ntrac))
    IF (.NOT. ALLOCATED(rsoil))   ALLOCATE (rsoil(ntrac))
    IF (.NOT. ALLOCATED(rws))     ALLOCATE (rws(ntrac))
    IF (.NOT. ALLOCATED(rwater))  ALLOCATE (rwater(ntrac))
    IF (.NOT. ALLOCATED(rsnow))   ALLOCATE (rsnow(ntrac))

    ! define boundary conditions for dry deposition scheme
    ! Note: for coupling with JSBACH replace EF_FILE with EF_MODULE and insert
    ! a bc_set call into JSBACH. Make this namelist controllable.
    ! explanation: %bc_... = how to apply
    !              %ef_... = how to get values


!>>gf - modify and re-organize the part below #244 and #226

    ! 1) -- LAI
    bc_lai%bc_domain = BC_BOTTOM
    bc_lai%bc_mode = BC_REPLACE
    bc_lai%ef_type = nlai_drydep_ef_type
    bc_lai%ef_template = 'surface_properties.nc'
    bc_lai%ef_varname = 'lai'
    bc_lai%ef_geometry = EF_LONLAT
    bc_lai%ef_timedef = EF_IGNOREYEAR
    !! additional options (just for illustration) - see mo_boundary_condition
    !! bc_lai%ef_timeoffset = 0._dp
    !! bc_lai%ef_timeindex = 1
    !! bc_lai%ef_value = 0._dp
    bc_lai%ef_factor = 1._dp
    bc_lai%ef_interpolate = EF_NOINTER   ! none
    bc_lai%ef_actual_unit = 'm2 m-2'


    ! 2) -- canopy height
    bc_hc%ef_varname = 'hc'
    bc_hc%bc_domain = BC_BOTTOM
    bc_hc%bc_mode = BC_REPLACE
    bc_hc%ef_type = EF_FILE
    bc_hc%ef_template = 'surface_properties.nc'
    bc_hc%ef_geometry = EF_LONLAT
    bc_hc%ef_timedef = EF_IGNOREYEAR
    bc_hc%ef_factor = 1._dp
    bc_hc%ef_interpolate = EF_NOINTER
    bc_hc%ef_actual_unit = 'm2 m-2'   


    ! 3) -- roughness length
    bc_z0m%ef_varname = 'z0m'
    bc_z0m%bc_domain = BC_BOTTOM
    bc_z0m%bc_mode = BC_REPLACE
    bc_z0m%ef_type = EF_FILE
    bc_z0m%ef_template = 'surface_properties.nc'
    bc_z0m%ef_geometry = EF_LONLAT
    bc_z0m%ef_timedef = EF_IGNOREYEAR
    bc_z0m%ef_factor = 1._dp
    bc_z0m%ef_interpolate = EF_NOINTER
    bc_z0m%ef_actual_unit = 'm2 m-2'  


    ! 4) -- forest fraction
    bc_forest%ef_varname = 'forest'
    bc_forest%bc_domain = BC_BOTTOM
    bc_forest%bc_mode = BC_REPLACE
    bc_forest%ef_type = EF_FILE
    bc_forest%ef_template = 'surface_properties.nc'
    bc_forest%ef_geometry = EF_LONLAT
    bc_forest%ef_timedef = EF_IGNOREYEAR
    bc_forest%ef_factor = 1._dp
    bc_forest%ef_interpolate = EF_NOINTER
    bc_forest%ef_actual_unit = 'm2 m-2'  

    ibc_lai    = bc_define('leaf area index',  bc_lai, 2, .TRUE.)
    ibc_hc     = bc_define('canopy height',    bc_hc, 2, .TRUE.)
    ibc_z0m    = bc_define('roughness length', bc_z0m, 2, .TRUE.)
    ibc_forest = bc_define('forest fraction',  bc_forest, 2, .TRUE.)

    ! 5) -- soil pH : variable names are F1..F7

    bc_soilph%bc_domain = BC_BOTTOM
    bc_soilph%bc_mode = BC_REPLACE
    bc_soilph%ef_type = EF_FILE
    bc_soilph%ef_template = 'xtsoil.nc'
    bc_soilph%ef_geometry = EF_LONLAT
    bc_soilph%ef_timedef = EF_CONSTANT
    bc_soilph%ef_factor = 1._dp
    bc_soilph%ef_interpolate = EF_NOINTER
    bc_soilph%ef_actual_unit = 'm2 m-2'  
    bc_soilph%ef_timeindex = 1

    DO j = 1, 7
      bc_soilph%ef_varname = 'F' // cnum(j:j)
      ibc_soilph(j) = bc_define('soil pH for soil type '//cnum(j:j), bc_soilph,  &
                                2, .TRUE.)
    END DO

!<<gf

    ! locate tracers for which resistances are directly calculated
    CALL get_tracer('SO2',idx=idt_so2, ierr=ierr)
    CALL get_tracer('H2SO4',idx=idt_h2so4, ierr=ierr)
    CALL get_tracer('O3',idx=idt_o3, ierr=ierr)
    CALL get_tracer('HNO3',idx=idt_hno3, ierr=ierr)
    CALL get_tracer('NO',idx=idt_no, ierr=ierr)
    CALL get_tracer('NO2',idx=idt_no2, ierr=ierr)
    CALL get_tracer('O3S',idx=idt_o3s, ierr=ierr)

    lhas_sulfur = (idt_so2 > 0 .OR. idt_h2so4 > 0)

!>>SF #57: additional security regarding to NO and NO2 mesophyllic resistance calculation
    IF (idt_no > 0 .AND. idt_o3 == 0) &
     CALL message('drydep_lg_init', 'NO defined as a tracer, but not O3', level=em_error)
    IF (idt_no2 > 0 .AND. idt_o3 == 0) &
     CALL message('drydep_lg_init', 'NO2 defined as a tracer, but not O3', level=em_error)
!<<SF #57

    ! mz_LG_20020115 call to subroutine in which the surface resistances,
    !     used for the dry deposition calculations, are initialized. 

    CALL calc_rs
  END SUBROUTINE drydep_lg_init

  SUBROUTINE init_drydep_lg_stream (sdrydep, trdetail, idt_ddep_detail)
  USE mo_tracdef,       ONLY: ln, trlist
  USE mo_tracer,        ONLY: get_tracer
  USE mo_exception,     ONLY: message, em_warn, finish
  USE mo_memory_base,   ONLY: t_stream, default_stream_setting, add_stream_element
  USE mo_decomposition, ONLY: ldc => local_decomposition

  TYPE (t_stream), POINTER      :: sdrydep
  CHARACTER(len=ln), INTENT(IN) :: trdetail
  INTEGER, INTENT(OUT)          :: idt_ddep_detail

  INTEGER                       :: ierr

 !++mgs 20140519 : detailed diagnostics
 ! identify tracer for detailed diagnostics
  CALL get_tracer(TRIM(trdetail), idx=idt_ddep_detail, ierr=ierr)
  IF (ierr /= 0) THEN
    IF (idt_o3 /= 0) THEN
      idt_ddep_detail=idt_o3
      CALL message('init_drydep_stream', &
      'Tracer '//TRIM(trdetail)//' for detailed diagnostics not defined! Using O3 instead!', &
      level=em_warn)
    ELSE
      IF (idt_so2 /= 0) THEN
        idt_ddep_detail=idt_so2
        CALL message('init_drydep_stream', &
        'Tracer '//TRIM(trdetail)//' for detailed diagnostics not defined! Using SO2 instead!', &
        level=em_warn)
      ELSE
        CALL finish('init_drydep_lg_stream','None of the tracers '//TRIM(trdetail)//', O3, SO2 '// &
                    'for detailed diagnostics are defined!')
      ENDIF
    ENDIF
  ENDIF
  ! set stream to instantaneous
  CALL default_stream_setting (sdrydep, laccu=.FALSE.)
  ! add stream elements
  CALL add_stream_element(sdrydep, 'lai', dpzlai,     &
                          longname='leaf area index', &
                          units='1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rmes', dpzrmes,     &
                          longname='bulk mesophyllic resistance', &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'zhc', dpzhc,     &
                          longname='canopy height', &
                          units='m', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'z0m', dpz0m,     &
                          longname='roughness length', &
                          units='m', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rahcan', dprahcan,     &
                          longname='aerodynamic canopy resistance', &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rco_leaf', dpzrco_leaf,     &
                          longname='leaf stomatal resistance', &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'ustveg', dppustveg,     &
                          longname='ustar vegetation', &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'ustslsn', dppustslsn,     &
                          longname='ustar soil and snow', &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'ustarw', dppustarw,     &
                          longname='ustar water', &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'ustari', dppustari,     &
                          longname='ustar ice', &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'snow_cover', dppcvs,     &
                          longname='fractional snow cover', &
                          units='1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'wetskin_cover', dppcvw,     &
                          longname='fractional wet skin cover', &
                          units='1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'vegetation_ratio', dppvgrat,     &
                          longname='fraction of vegetated land (old ECHAM variable)', &
                          units='1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rah_veg', dpprahveg,     &
                          longname='aerodynamic resistance in vegetation', &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rah_slsn', dpprahslsn,     &
                          longname='aerodynamic resistance for bare soil and snow', &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rah_ice', dpprahice,     &
                          longname='aerodynamic resistance over ice', &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rc0x', dprc0x,     &
                          longname='bulk stomatal resistance', &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'r_snow_hno3so2', dprsnowhno3so2,     &
                          longname='snow resistance hno3so2', &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'r_leaf', dprleaf,     &
                          longname='leaf resistance of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'r_stom', dprstom,     &
                          longname='stomatal resistance of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                        units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rs_veg', dprsveg,     &
                          longname='vegetation resistance Rs of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rah_water', dpprahwat,     &
                          longname='aerodynamic resistance over water of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rb_soil_snow', dprbslsn,     &
                          longname='surface resistance for bare soil and snow of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rb_veg', dprbveg,     &
                          longname='vegetation surface resistance of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rb_water', dprbw,     &
                          longname='water surface resistance of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'rb_ice', dprbi,     &
                          longname='ice surface resistance of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='s m-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'vd_veg', dpvdveg,     &
                          longname='deposition velocity in vegetation of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'vd_soil', dpvdsoil,     &
                          longname='deposition velocity for soil of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'vd_wetskin', dpvdws,     &
                          longname='deposition velocity for wet skin of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'vd_snow', dpvdsn,     &
                          longname='deposition velocity for snow of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'vd_land', dpvdland,     &
                          longname='total deposition velocity over land of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'vd_wat', dpvdwat,     &
                          longname='total deposition velocity over water of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='m s-1', lpost=.TRUE.)
  CALL add_stream_element(sdrydep, 'vd_ice', dpvdice,     &
                          longname='total deposition velocity over ice of '//TRIM(trlist%ti(idt_ddep_detail)%fullname), &
                          units='m s-1', lpost=.TRUE.)
  END SUBROUTINE init_drydep_lg_stream


  !=============================================================================

  SUBROUTINE calc_rs
  !-----------------------------------------------------------------------------
  ! subroutine calc_rs, to calculate the values of the uptake resistances 
  ! required to calculate the trace gas dry deposition velocity. this routine 
  ! is based on an approach by wesely, 1989, in which the uptake resistances of
  ! trace gases, for which the dry deposition velocities have not been observed,
  ! are estimated based on the henry coefficient and a reactivity coefficient and
  ! the uptake resistances of so2 and o3, of the "big leaf" dry deposition scheme
  ! by ganzeveld and j. lelieveld j. geophys. res., 100, 20,999-21,012,1995,
  ! ganzeveld et al.,j. geophys. res., 103, 5679-5694, 1998 and ganzeveld et al, 
  ! submitted to j. geophys. res., 2001. for more information of the wesely 
  ! approach see atmospheric environment vol 23, no 6, 1293-1304.
  !
  ! the program needs as input data the molecular mass of the defined trace 
  ! gases, the henry coefficient [m atm-1] and an estimated reactivity 
  ! coefficient which has 3 distinct values: 0 for non-reactive species, 
  ! (e.g, so2, acetaldehyde), 0.1 for moderately reactive species (e.g., pan),
  ! and 1 for reactive species (e.g., o3, hno3). these values are defined in
  ! the module mo_moz_init.
  !-----------------------------------------------------------------------------
  !   Author:
  !   -------
  !   Laurens Ganzeveld, MPI Mainz                                   2001
  !
  !   Modifications:
  !   --------------
  !   Laurens Ganzeveld, MPI Mainz and
  !   Philip Stier,      MPI Hamburg (implementation in ECHAM5) 2001-2002

  USE mo_exception,    ONLY: message_text, message, em_param, em_warn
  USE mo_util_string,  ONLY: separator
  USE mo_tracdef,      ONLY: ntrac, trlist, GAS
  USE mo_species,      ONLY: speclist

  ! local variables
  INTEGER  :: jt, ispec

  REAL(dp) :: diffrb_so2, rsoil_so2, rwater_so2, rws_so2,  &
              rsnow_so2,  rmes_so2,  rcut_so2,   diff_so2, &
              diffrb_o3,  rsoil_o3,  rwater_o3,  rws_o3,   &
              rsnow_o3,   rmes_o3,   rcut_o3,    diff_o3

  LOGICAL :: lo_derived

  !--- attribute specific parameters in the following order:
  !
  !    - soil resistance
  !    - sea water resistance, which is generally similar to the wet skin
  !    - wet skin reservoir resistance
  !    - snow resistance
  !    - mesophyll resistance
  !    - cuticle resistance
  !    - diffusivity coefficient, to correct stomatal resistance for 
  !      differences in diffusivity between water vapour and the 
  !      specific trace gas (sqrt(molmass trace gas)/sqrt(molmass h2o))

  !--- Values for SO2 and O3 are required in any case to estimate the 
  !    resistances of the other species considered in the deposition scheme.
  !    => Define resistances for SO2 and O3 as defauilt in case they were
  !       not defined as tracers.


  diffrb_so2=1.6_dp
  rsoil_so2=250._dp
  rwater_so2=1._dp
  rws_so2=100._dp
  rsnow_so2=1._dp
  rmes_so2=1._dp
  rcut_so2=1.e5_dp
  diff_so2=1.9_dp

  IF (idt_so2 > 0) THEN
     diffrb(idt_so2)=diffrb_so2
     rsoil(idt_so2)=rsoil_so2
     rwater(idt_so2)=rwater_so2
     rws(idt_so2)=rws_so2
     rsnow(idt_so2)=rsnow_so2
     rmes(idt_so2)=rmes_so2
     rcut(idt_so2)=rcut_so2
     diff(idt_so2)=diff_so2
  END IF

  diffrb_o3=1.2_dp
  rsoil_o3=400._dp
  rwater_o3=2000._dp
  rws_o3=2000._dp
  rsnow_o3=2000._dp
  rmes_o3=1._dp
  rcut_o3=1.e5_dp
  diff_o3=1.6_dp

  IF (idt_o3 > 0) THEN
     diffrb(idt_o3)=diffrb_o3
     rsoil(idt_o3)=rsoil_o3
     rwater(idt_o3)=rwater_o3
     rws(idt_o3)=rws_o3
     rsnow(idt_o3)=rsnow_o3
     rmes(idt_o3)=rmes_o3
     rcut(idt_o3)=rcut_o3
     diff(idt_o3)=diff_o3
  END IF

  IF (idt_h2so4 > 0) THEN
     diffrb(idt_h2so4)=1.8_dp
     rsoil(idt_h2so4)=1.e5_dp
     rwater(idt_h2so4)=1.e5_dp
     rws(idt_h2so4)=1.e5_dp
     rsnow(idt_h2so4)=1.e5_dp
     rmes(idt_h2so4)=1.e5_dp
     rcut(idt_h2so4)=1.e5_dp
     diff(idt_h2so4)=2.7_dp
  END IF

  IF (idt_hno3 > 0) THEN
     diffrb(idt_hno3)=1.4_dp
     rsoil(idt_hno3)=1._dp
     rwater(idt_hno3)=1._dp
     rws(idt_hno3)=1._dp
     rsnow(idt_hno3)=1._dp
     rmes(idt_hno3)=1._dp
     rcut(idt_hno3)=1._dp
     diff(idt_hno3)=1.9_dp
  END IF

  IF (idt_no > 0) THEN
     diffrb(idt_no)=1.1_dp
     rsoil(idt_no)=1.e5_dp
     rwater(idt_no)=1.e5_dp
     rws(idt_no)=1.e5_dp
     rsnow(idt_no)=1.e5_dp
     rmes(idt_no)=500._dp
     rcut(idt_no)=1.e5_dp
     diff(idt_no)=1.3_dp
  END IF

  IF (idt_no2 > 0) THEN
     diffrb(idt_no2)=1.2_dp
     rsoil(idt_no2)=600._dp
!!   rwater(idt_no2)=1.e5_dp
!!   rws(idt_no2)=1.e5_dp
!!   rsnow(idt_no2)=1.e5_dp
     rwater(idt_no2)=2000._dp    ! new value of www.emep.int/acid/ladm.html (Tsyro), old: 1.e5
     rws(idt_no2)=2000._dp       ! wet skin resistance (is wet urban or desert value of lit. above, old: 1.e5)
     rsnow(idt_no2)=2000._dp     ! snow value of lit above

     rmes(idt_no2)=1._dp
     rcut(idt_no2)=1.e5_dp
     diff(idt_no2)=1.6_dp
  END IF

  !-- header for reporting on dry deposition parameters
  CALL message('', separator)
  CALL message('calcrs', 'Derived parameters for dry deposition:', level=em_param)
  write(message_text,'(A20,7A13)') 'tracer', 'diffrb', 'rmes', 'rcut', 'rsoil',   &
                                   'rwater', 'rsnow', 'rws'
  CALL message('', message_text, level=em_param)

  DO jt=1, ntrac
     ispec = trlist%ti(jt)%spid
     lo_derived=trlist%ti(jt)%basename/='SO2'        .AND. &
!==DT                trlist%ti(jt)%fullname/='SO4_gas'    .AND. &     ! HAM 
                trlist%ti(jt)%basename/='H2SO4'      .AND. &     ! MOZ
                trlist%ti(jt)%basename/='O3'         .AND. &
                trlist%ti(jt)%basename/='HNO3'       .AND. &
                trlist%ti(jt)%basename/='NO'         .AND. &
                trlist%ti(jt)%basename/='NO2'

!! ### where do we compute resistances for aerosols??? ###
     IF ((trlist%ti(jt)%ndrydep==2) .AND. (trlist%ti(jt)%nphase==GAS)) THEN

        IF (lo_derived) THEN

          ! calculation of term which is used to correct the stomatal resistance
          ! for differences in the diffusitivy (see also equation 4).

          diff(jt)=SQRT(speclist(ispec)%moleweight/18._dp)

          ! calculation of the term to correct for differences in diffusivity 
          ! between water vapor and the trace gas. it is calculated from: 
          ! diff bl=(v/dx)**2/3, with v=0.189 sm-1 and dx= dh2o/sqrt(mh2o/mx), 
          ! with dh2o=0.212

          diffrb(jt)=(0.189_dp/(0.212_dp/diff(jt)))**(2._dp/3._dp)

          ! calculation of rmx, the mesophyll resistance
!### WORKAROUND
IF(speclist(ispec)%henry(1)==0._dp .AND. speclist(ispec)%dryreac==0._dp) THEN
  WRITE(message_text,'(2a)') 'henry and dryreac == 0 for ',speclist(ispec)%shortname
  CALL message('calc_rs', message_text, level=em_warn)
  speclist(ispec)%henry(1) = 1.e-15_dp
END IF

          rmes(jt)=1._dp/(speclist(ispec)%henry(1)/3000._dp+100._dp*speclist(ispec)%dryreac)

          ! calculation of rlux, the cuticular resistance, equation 7 of wesely's
          ! paper

          rcut(jt)=1._dp/(1.e-5_dp*speclist(ispec)%henry(1)+speclist(ispec)%dryreac)*rcut_o3

          ! calculation of rgsx, the soil resistance, equation 9 of wesely's
          ! paper

          rsoil(jt)=1._dp/(speclist(ispec)%henry(1)/(1.e5_dp*rsoil_so2)+ &
               speclist(ispec)%dryreac/rsoil_o3)

          ! the snow resistance is similar as the soil resistance

          rsnow(jt)=rsoil(jt)

          ! calculation of rlux-wet, the wet skin resistance, equation 14 of 
          ! wesely's paper
!++sschr 20140519: fix for methane!
!!        rws(jt)=1._dp/(1._dp/(3._dp*rws_so2)+1.e-7_dp*speclist(ispec)%henry(1)+ &
!!             speclist(ispec)%dryreac/rws_o3) 
          rws(jt)=1._dp/(1._dp/(3._dp*15000._dp)+1.e-7_dp*speclist(ispec)%henry(1)+ &
               speclist(ispec)%dryreac/rws_o3) 
!--sschr

          ! calculation of sea uptake resistance, using equation 9 of wesely's
          ! paper

          rwater(jt)=1._dp/(speclist(ispec)%henry(1)/(1.e5_dp*rwater_so2)+ &
               speclist(ispec)%dryreac/rwater_o3)

        END IF

      ! output derived information 
      write(message_text,'(A20,7G13.4)') trlist%ti(jt)%fullname, diffrb(jt), rmes(jt),  &
                                        rcut(jt), rsoil(jt), rwater(jt), rsnow(jt), rws(jt)
      CALL message('', message_text, level=em_param)
    END IF

  ENDDO

  CALL message('', separator)

  END SUBROUTINE calc_rs


  !! ------------------------------------------------------------------------------------------
  !! LG- Subroutine in which the aerodynamic resistance (Ra) is being 
  !!     calculated This resistance is used in the calculation of the 
  !!     dry deposition velocity. Ra (and also the friction velocity u*
  !!     are being calculated for each surface cover fraction to consider
  !!     pronounced differences in the surface roughness   
  !!   
  !!     Laurens Ganzeveld, 1998 (see for references the papers mentioned
  !!     in the dry deposition model vdbl), modified for implementation
  !!     in echam5, October, 2001


  SUBROUTINE drydep_lg_calcra (kproma,  kbdim,   klev,    krow,                    &
!>>gf change in argument list (see #78)
                               pepdu2,  pkap,    pum1,    pvm1,    pgeom1,         &
                               pril,    priw,    prii,                             &
                               ptvir1,  ptvl,    ptvw,    ptvi,    ptsm1m,         &
                               loland,                                             &
                               pcdnl,   pcdnw,   pcdni,   pcfml,                   &
                               pcfmw,   pcfmi,   pcfncl,  pcfncw,  pcfnci,         &
!<<gf
                               paz0w,   paz0i,   paz0l,                            &
                               prahwat, prahice, prahveg, prahslsn,                &
!>>gf change in argument list (see #78)
                               pustarl, pustarw, pustari, pustveg, pustslsn        )
!<<gf

  USE mo_physical_constants, ONLY: grav
  USE mo_boundary_condition, ONLY: bc_apply

  ! -- parameters
  INTEGER, INTENT(in)   :: kproma, kbdim, klev, krow

  LOGICAL, INTENT(in)   :: loland(kbdim)

  REAL(dp), INTENT(in)  :: pcfml(kbdim),      & !
                           pcfmw(kbdim),      & !
                           pcfmi(kbdim),      & !

                           pcdnl(kbdim),      & !
                           pcdnw(kbdim),      & !
                           pcdni(kbdim),      & !

                           pcfncl(kbdim),     & !
                           pcfncw(kbdim),     & !
                           pcfnci(kbdim),     & !
                           pum1(kbdim,klev),  & !
                           pvm1(kbdim,klev),  & !
                           pgeom1(kbdim,klev),& !
                           pril(kbdim),       & ! 
                           priw(kbdim),       & !
                           prii(kbdim),       & !
                           ptvir1(kbdim,klev),& !
                           ptvl(kbdim),       & ! virtual potential temp. over land
                           ptvw(kbdim),       & ! virtual potential temp. over water
                           ptvi(kbdim),       & ! virtual potential temp. over ice
                           ptsm1m(kbdim),     & !    unused (???)
                           paz0l(kbdim),      & !    unused (???)
                           paz0w(kbdim),      & !
                           paz0i(kbdim)

  REAL(dp), INTENT(out) :: prahwat(kbdim),    & !
                           prahice(kbdim),    & !
                           prahveg(kbdim),    & !
                           prahslsn(kbdim),   & !
                           pustarl(kbdim),    & !
                           pustarw(kbdim),    & !
                           pustari(kbdim),    & !
                           pustveg(kbdim),    & !
                           pustslsn(kbdim)  

!>>gf see #75 changed zkmh. The previous value of 0.74, was derived from a von Karman constant
!             of 0.34, whereas we use 0.40. Therefore zkmkh should be equal to 0.95.
  REAL(dp), PARAMETER :: zkmkh=0.95_dp
  REAL(dp), PARAMETER :: z0slsn = 1.E-4_dp
!<<gf
!>>gf see #78
  REAL(dp), PARAMETER :: cfncmin = 1.E-6_dp
  REAL(dp), PARAMETER :: ustarmin = 1.E-5_dp
!<<gf

  ! -- local variables
  INTEGER             :: jl

  REAL(dp) ::                                           &
       cmveg, cmslsn,cml, cmw, cmi,              &
       zcfncl(kbdim), zcfncw(kbdim), zcfnci(kbdim),     &
       zcdnveg, zcdnslsn
  REAL(dp)            :: z0m(kbdim)

  REAL(dp)            :: pepdu2,pkap
  REAL(dp) :: zsurf, zmoninw, zmonini, &
              zoverlw,zoverli,zxzsurf,   &
              zws
  REAl(dp) :: zmon_pre, zmoninveg, zmoninslsn,       &
              zoverlveg, zoverlslsn, zxzsurfveg,     & 
              zxzsurfslsn
  REAl(dp) :: zpsihw, zpsihi, zpsihveg, zpsihslsn
  REAL(dp) :: z0h, zxz0h

  ! -- code
  ! -- obtain roughness length
  CALL bc_apply(ibc_z0m, kproma, krow, z0m)

!>>gf-safety
  zcfncl(1:kproma) = MAX(pcfncl(1:kproma), cfncmin)
  zcfncw(1:kproma) = MAX(pcfncw(1:kproma), cfncmin)
  zcfnci(1:kproma) = MAX(pcfnci(1:kproma), cfncmin)
!<<gf-safety

  DO jl=1,kproma

     !---surface wind speed:
     zws = SQRT(MAX(pepdu2, pum1(jl,klev)**2+pvm1(jl,klev)**2))  !gf #75

     ! LG- calculation of drag coefficient and u*. A change with
     !     the previous calculations in echam4 is that the exchange 
     !     coefficients, calculated in echam5 over land, ice and water,
     !     is being used to explicitly calculate the aerodynamic 
     !     resistances over these surface

!>>gf see #78 and #75
     cml=pcdnl(jl)*pcfml(jl)/zcfncl(jl)
     cmw=pcdnw(jl)*pcfmw(jl)/zcfncw(jl)
     cmi=pcdni(jl)*pcfmi(jl)/zcfnci(jl)

     pustarl(jl)=SQRT(cml)*zws
     pustarw(jl)=SQRT(cmw)*zws
     pustari(jl)=SQRT(cmi)*zws

     !---vegetated land surface
     zcdnveg = (pkap/LOG(1._dp+pgeom1(jl,klev)/ &
               (grav*MAX(0.02_dp,z0m(jl)))))**2

     cmveg = zcdnveg*pcfml(jl)/zcfncl(jl)

     pustveg(jl) = SQRT(cmveg)*zws

     !---bare soil and snow cover
     zcdnslsn = (pkap/LOG(1._dp+pgeom1(jl,klev)/ &
                (grav*z0slsn)))**2

     cmslsn = zcdnslsn*pcfml(jl)/zcfncl(jl)

     pustslsn(jl) = SQRT(cmslsn)*zws

!>>gf-safety
     pustarl(jl)  = MAX(pustarl(jl), ustarmin)
     pustarw(jl)  = MAX(pustarw(jl), ustarmin)
     pustari(jl)  = MAX(pustari(jl), ustarmin)
     pustveg(jl)  = MAX(pustveg(jl), ustarmin)
     pustslsn(jl) = MAX(pustslsn(jl),ustarmin)
!<<gf-safety

!<<gf
     ! LG- Computation of stability correction term, The stability
     !     correction functions are taken from Stull (page 383-385) 
     !     (08-11-98) and are slightly different from those by Williams 
     !     and Hicks et al., which were originally being used in the dry 
     !     deposition scheme. 

     zsurf=(pgeom1(jl,klev)/grav)

!>>gf see #78 and #75

     !--- Land:
     IF (pril(jl).GE.0._dp) THEN

        ! LG- calculating the Monin-Obukhov length directly applying the 
        !     formula given by Stull, 9.7.5k, page 386

        zmon_pre = (((ptvir1(jl,klev)+ptvl(jl))/2._dp)* zws)/ &
             (pkap*grav*(ptvir1(jl,klev)-ptvl(jl)))
        zmoninveg = pustveg(jl) * zmon_pre
        zoverlveg = zsurf/zmoninveg

        zpsihveg = zkmkh*LOG(zsurf / MAX(0.02_dp,z0m(jl))) + 7.8_dp*zoverlveg

        zmoninslsn = pustslsn(jl) * zmon_pre
        zoverlslsn = zsurf/zmoninslsn

        zpsihslsn = zkmkh*LOG(zsurf / MAX(0.02_dp,z0m(jl))) + 7.8_dp*zoverlslsn
     ELSE
        zmoninveg = zsurf/pril(jl)
        zoverlveg = zsurf / zmoninveg

        zxzsurfveg=SQRT(1._dp-11.6_dp*zoverlveg)

        z0h = MIN(1._dp,z0m(jl))
        z0h = MAX(z0h, 1.e-5_dp)
        zxz0h=SQRT(1._dp-11.6_dp*z0h/zmoninveg)

        zpsihveg = zkmkh * ( LOG((zxzsurfveg-1._dp)/(zxzsurfveg+1._dp)) - LOG((zxz0h-1._dp)/(zxz0h+1._dp)) )

        zmoninslsn = zsurf/pril(jl)
        zoverlslsn = zsurf / zmoninslsn

        zxzsurfslsn=SQRT(1._dp-11.6_dp*zoverlslsn)
        zxz0h=SQRT(1._dp-11.6_dp*z0h/zmoninslsn)

        zpsihslsn = zkmkh * ( LOG((zxzsurfslsn-1._dp)/(zxzsurfslsn+1._dp)) - LOG((zxz0h-1._dp)/(zxz0h+1._dp)) )
     ENDIF

     prahveg(jl)=MAX(1._dp,(zpsihveg/(pustveg(jl)*pkap))) 

     prahslsn(jl)=MAX(1._dp,(zpsihslsn/(pustslsn(jl)*pkap)))

     !--- Water:
     IF (priw(jl).GE.0._dp) THEN

        zmoninw = (pustarw(jl)*((ptvir1(jl,klev)+ptvw(jl))/2._dp)* zws)/ &
             (pkap*grav*(ptvir1(jl,klev)-ptvw(jl)))

        zoverlw=zsurf/zmoninw
        zpsihw = zkmkh*LOG(zsurf/paz0w(jl)) + 7.8_dp*zoverlw
     ELSE
        zmoninw=zsurf/priw(jl)
        zoverlw = priw(jl)  !---Stull 9.7.5j
        zxzsurf=SQRT(1._dp-11.6_dp*zoverlw)
         
        z0h = MIN(1._dp,paz0w(jl))
        z0h = MAX(z0h, 1.e-5_dp)
        zxz0h=SQRT(1._dp-11.6_dp*z0h/zmoninw)

        zpsihw = zkmkh * ( LOG((zxzsurf-1._dp)/(zxzsurf+1._dp)) - LOG((zxz0h-1._dp)/(zxz0h+1._dp)) )

     ENDIF

     prahwat(jl)=MAX(1._dp,(zpsihw/(pustarw(jl)*pkap)))

     !--- Ice:
     IF (prii(jl).GE.0._dp) THEN

        zmonini = (pustari(jl)*((ptvir1(jl,klev)+ptvi(jl))/2._dp)* zws)/ &
             (pkap*grav*(ptvir1(jl,klev)-ptvi(jl)))

        zoverli=zsurf/zmonini
        zpsihi = zkmkh*LOG(zsurf/paz0i(jl)) + 7.8_dp*zoverli
     ELSE

        zmonini=zsurf/prii(jl)
        zoverli=zsurf/zmonini

        zxzsurf=SQRT(1._dp-11.6_dp*zoverli)

        z0h = MIN(1._dp,paz0i(jl))
        z0h = MAX(z0h, 1.e-5_dp)
        zxz0h=SQRT(1._dp-11.6_dp*z0h/zmonini)

        zpsihi = zkmkh * ( LOG((zxzsurf-1._dp)/(zxzsurf+1._dp)) - LOG((zxz0h-1._dp)/(zxz0h+1._dp)) )

     ENDIF
!<<gf

     prahice(jl)=MAX(1._dp,(zpsihi/(pustari(jl)*pkap)) )

     ! Computation of relative humidity out T and Tdew (2m) (Monteith)

     ! LG- assigning the temperature, originally the 2 m temperature was
     !     used in the dry deposition scheme, however, since the 2 m
     !     temperature is calculated at the end of this routine and these
     !     particular calculations have been moved to this location,
     !     the surface temperature which is available here is used instead
     !     of the parameter PTEMP2M.

!     ztsurf=MAX(1.,ptsm1m(jl))
!     zes=0.611*EXP(MIN(10.,(19.59*(ztsurf-tmelt)/ztsurf)))
!     ze=0.611*EXP(MIN(10.,(19.59*(1.-tmelt/ &
!          MIN(ztsurf,(tmelt-zfrac*zcvm4)/(1.-zfrac))))))
!     rh(jl)=ze/zes

     ! LG- end

  ENDDO


  END SUBROUTINE drydep_lg_calcra


  !! ------------------------------------------------------------------------------------------
  !! drydep_lg_vdbl calculates deposition velocities for gases using the big leaf approach
  !! This is a merger from the former moz_vdbl and ham_vdbl
  !! This routine uses the ECHAM surface cover characterization and surface parameters
  !! See for more information about the code the papers by
  !!     Ganzeveld and Lelieveld, JGR 100, 1995,
  !!     Ganzeveld et al., JGR 103, 1998
  !!     Ganzeveld et al., 2001

  SUBROUTINE drydep_lg_vdbl (kproma,   kbdim,     klev,    krow,    laland,  &
                             psrfl,    ptsm,      pum1,    pvm1,    prh,     &
                             pfrl,     pfrw,      pfri,    pcvbs,   pcvs,    &
                             pcvw,     pvgrat,    prahwat, prahice, prahveg, &
                             prahslsn, pws,       pwsmx,                     &
                             pustarw,  pustari,   pustveg, pustslsn,         &
                             pvd,      pvdstom,   idt_ddep_detail            )

  USE mo_tracdef,            ONLY: trlist, ntrac, GAS
  USE mo_boundary_condition, ONLY: bc_apply
  USE mo_exception,          ONLY: finish
  USE mo_physical_constants, ONLY: tmelt
  USE mo_submodel_streams,   ONLY: drydep_ldetail

  !--- parameters
  INTEGER, INTENT(in)      :: kproma, kbdim, klev, krow
  LOGICAL, INTENT(in)      :: laland(kbdim)
  REAL(dp), INTENT(in)     :: psrfl(kbdim),        & ! surface solar flux
                              ptsm(kbdim),         & ! surface temperature
                              pum1(kbdim,klev),    & ! wind vector component
                              pvm1(kbdim,klev),    & ! wind vector component!
                              prh(kbdim),          & ! relative humidity (lowest model layer)
                              pfrl(kbdim),         & ! land fraction
                              pfrw(kbdim),         & ! water fraction
                              pfri(kbdim),         & ! ice fraction
                              pcvbs(kbdim),        & ! bare soil fraction  (### not used)
                              pcvs(kbdim),         & ! snow cover fraction
                              pcvw(kbdim),         & ! wet skin fraction
                              pvgrat(kbdim),       & ! vegetation ratio
                              prahwat(kbdim),      & ! aerodynamic resistance (?) over water
                              prahice(kbdim),      & ! dto. over ice
                              prahveg(kbdim),      & ! dto. over vegetation
                              prahslsn(kbdim),     & ! dto. over soil and snow
                              pws(kbdim),          & ! soil moisture
                              pwsmx(kbdim),        & ! field capacity of soil
                              pustarw(kbdim),      & ! friction velocity over water
                              pustari(kbdim),      & ! friction velocity over ice
                              pustveg(kbdim),      & ! ustar for vegetation
                              pustslsn(kbdim)        ! ustar for (bare) soil and snow
  REAL(dp), INTENT(inout)  :: pvd(kbdim, ntrac),   & ! deposition velocities
                              pvdstom(kbdim, ntrac)  ! dto. for stomatal uptake
  INTEGER, INTENT(in)      :: idt_ddep_detail

  !--- local variables
  REAL(dp) :: zrmes(kbdim)    ! locally dependent rmes
  REAL(dp) :: fsic(kbdim)
  REAL(dp) :: rbw(kbdim,ntrac),                       & !              
              rbi(kbdim,ntrac),                       & !
              rbveg(kbdim,ntrac),rbslsn(kbdim,ntrac), & !
              rahcan(kbdim),                          & !
              rc0x(kbdim,ntrac),                      & !
              rsveg(kbdim,ntrac), rleaf(kbdim,ntrac), & !
              rstom(kbdim,ntrac),                     & !   stomatal resistance ++mgs
              rsoilso2(kbdim), rsnowhno3so2(kbdim),   & !
              stheta(kbdim)
  REAL(dp) :: vdveg(kbdim,ntrac), vdsoil(kbdim,ntrac),   &
              vdstom(kbdim,ntrac),                       &  ! stomatal drydep velocity
              vdland(kbdim,ntrac), vdwat(kbdim,ntrac),   &
              vdice(kbdim,ntrac),                        & !gf see #78
              vdws(kbdim,ntrac), vdsn(kbdim,ntrac),      &
              vdso4slsn(kbdim), vdso4wat(kbdim), vdso4ice(kbdim) !gf added vdso4ice, see #78

  REAL(dp) :: vdvfac(kbdim)
  REAL(dp) :: zlai(kbdim), zhc(kbdim), z0m(kbdim)
  REAL(dp) :: zfws(kbdim),               & ! soil moisture attenuation factor
              zrco_leaf(kbdim)             ! leaf stomatal resistance

  INTEGER  :: jl, jt
  LOGICAL, SAVE :: lfirst = .TRUE.


  ! Error check
  IF (lfirst .AND. idt_h2so4 > 0 .AND. idt_so2 == 0) &
     CALL finish('drydep_lg_vdbl', 'Tracer H2SO4 defined, but SO2 not defined (idt_so2==0)!')    

  ! obtain lai, hc, z0m and soil pH from boundary conditions
  CALL bc_apply(ibc_lai, kproma, krow, zlai)
  CALL bc_apply(ibc_hc,  kproma, krow, zhc)
  CALL bc_apply(ibc_z0m, kproma, krow, z0m)

  DO jl=1,kproma

    ! mz_LG_20020115 from now on using LAI instead of the parameter AREA (05-2000),
    !     for calculation of within canopy aerodynamic resistance

    rahcan(jl)=MAX(1._dp,14._dp*zlai(jl)*zhc(jl)/pustveg(jl))

    ! mz_LG_20020115 calculation of the HNO3 and SO2 snow resistance as a function of the
    !     snow temperature

    rsnowhno3so2(jl)=MIN(MAX(10._dp,10._dp**(-0.09_dp*(ptsm(jl)-tmelt) &
                             +2.4_dp)),1.e5_dp)

  ENDDO

  !--- calculate soil moisture attenuation factor and ...
  CALL vd_prep(kbdim, kproma, krow,              &
               pws,   pwsmx,  psrfl,             &
               zfws,  zrco_leaf                   )

  !--- get soil resistance for SO2 and vd for sulphate
  IF (lhas_sulfur) CALL vd_prep_sulfur(kbdim, kproma, krow,                      &
                                       pustarw, pustari, pustslsn, ptsm,         &
                                       prh, psrfl, prahveg,                      &
                                       pum1(:,klev), pvm1(:,klev), &
                                       rsoilso2, vdso4wat, vdso4ice,    &
                                       vdso4slsn,                       &
                                       stheta                           )

!++jsr, 20040901 compare LG's new vdbl routine
!       must compute this outside jt loop, because the rc0x values of ozone
!       are used in the calculation of the resistances for NO,NO2
  IF (idt_o3 > 0) THEN
    rc0x(1:kproma,idt_O3)=diff(idt_O3)*zrco_leaf(1:kproma)/ &
         MAX(1.e-5_dp,zfws(1:kproma)) ! correcting for soil moisture
  ENDIF
!--jsr


  ! Loop for tracers, for calculation of tracer specific resistances

  DO jt=1,ntrac

    ! limit calculation of vd to gas-phase tracers fo rwhich Ganzeveld scheme has 
    ! been selected (ndrydep==2)

    IF (trlist%ti(jt)%ndrydep /= 2 .OR. trlist%ti(jt)%nphase /= GAS) CYCLE

    ! initialize return variables
    pvd(1:kproma, jt)     = 0._dp
    pvdstom(1:kproma, jt) = 0._dp

    ! mz_LG_20020115 calculation of quasi-laminar boundary layer resistances

    rbveg(1:kproma,jt)=(2._dp/(pustveg(1:kproma)*0.40_dp))*diffrb(jt)
    rbslsn(1:kproma,jt)=(1._dp/(pustslsn(1:kproma)*0.40_dp))*diffrb(jt)
!>>gf see #78
    rbw(1:kproma,jt)=(1._dp/(pustarw(1:kproma)*0.40_dp))*diffrb(jt)
    rbi(1:kproma,jt)=(1._dp/(pustari(1:kproma)*0.40_dp))*diffrb(jt)
!<<gf

    ! mz_LG_20020115 Stomatal resistance computation for each species, by
    !     multiplying the leaf stomatal resistance with the
    !     diffusivity term

    IF (jt /= idt_O3) THEN
       rc0x(1:kproma,jt)=diff(jt)*zrco_leaf(1:kproma)
       ! correcting for soil moisture
       rc0x(1:kproma,jt)=rc0x(1:kproma,jt)/MAX(1.e-5_dp,zfws(1:kproma)) 
    END IF

    ! mz_LG_20020115 definition of the NO and NO2 mesophyllic resistance as a function
    !     of the ozone mesophyllic resistances (see Ganzeveld and Lelieveld, '95)

    zrmes(1:kproma) = rmes(jt)
    IF (idt_o3 > 0 .AND. jt == idt_no) THEN         !SF bugfix #57
      zrmes(1:kproma)=5._dp*rc0x(1:kproma,idt_o3)
    END IF
    IF (idt_o3 > 0 .AND. jt == idt_no2) THEN        !SF bugfix #57
      zrmes(1:kproma)=0.5_dp*rc0x(1:kproma,idt_o3)
    END IF

    ! mz_LG_20020115  update March 2000, RLEAF instead of RVEG, and RSVEG (surface
    !      resistance for vegetated fraction) instead of RSVEG

    rleaf(1:kproma,jt)=(1._dp/((1._dp/rcut(jt))+(1._dp/(rc0x(1:kproma,jt)+zrmes(1:kproma)))))
    !++mgs added stomatal resistance for diagnostics
    rstom(1:kproma,jt)=(rc0x(1:kproma,jt)+zrmes(1:kproma))/MAX(1.e-5_dp,zlai(1:kproma))
    rsveg(1:kproma,jt)=(1._dp/((1._dp/(rahcan(1:kproma)+ &
                        (1._dp/(pustveg(1:kproma)*0.40_dp))*diffrb(jt)+rsoil(jt)))+ &
                        (1._dp/(rleaf(1:kproma,jt)/MAX(1.e-5_dp,zlai(1:kproma))))))

    ! mz_LG_20020115 end
    !--- compute vd dependent on surface type
    !    special treatment for HNO3, SO2, H2SO4

    IF (jt==idt_hno3) THEN
      rsveg(1:kproma,jt)=MAX(10._dp,rsveg(1:kproma,jt))
      vdveg(1:kproma,jt)=(1._dp/(prahveg(1:kproma)+rbveg(1:kproma,jt)+rsveg(1:kproma,jt)))
      vdsoil(1:kproma,jt)=(1._dp/(prahslsn(1:kproma)+rbslsn(1:kproma,jt)+rsoil(jt)))
      vdws(1:kproma,jt)=(1._dp/(prahveg(1:kproma)+rbveg(1:kproma,jt)+rws(jt)))
      vdwat(1:kproma,jt)=1._dp/(prahwat(1:kproma)+rbw(1:kproma,jt)+rwater(jt))
      vdsn(1:kproma,jt)= (1._dp/(prahslsn(1:kproma)+rbslsn(1:kproma,jt)+rsnowhno3so2(1:kproma)))
      vdice(1:kproma,jt)= (1._dp/(prahice(1:kproma)+rbi(1:kproma,jt)+rsnowhno3so2(1:kproma)))

    ELSE IF (jt==idt_so2) THEN
      vdveg(1:kproma,jt)=(1._dp/(prahveg(1:kproma)+rbveg(1:kproma,jt)+rsveg(1:kproma,jt)))
      vdsoil(1:kproma,jt)=(1._dp/(prahslsn(1:kproma)+rbslsn(1:kproma,jt) &
                                  +rsoilso2(1:kproma)))
      vdws(1:kproma,jt)=(1._dp/(prahveg(1:kproma)+rbveg(1:kproma,jt)+rws(jt)))
      vdwat(1:kproma,jt)=1._dp/(prahwat(1:kproma)+rbw(1:kproma,jt)+rwater(jt))
      vdsn(1:kproma,jt)= (1._dp/(prahslsn(1:kproma)+rbslsn(1:kproma,jt)+rsnowhno3so2(1:kproma)))
      vdice(1:kproma,jt)= (1._dp/(prahice(1:kproma)+rbi(1:kproma,jt)+rsnowhno3so2(1:kproma)))

    ELSE IF (jt==idt_h2so4) THEN
      IF (trlist%ti(jt)%nphase==GAS) THEN
      !--- Assume for gaseous H2SO4 the same dry deposition velocity as for SO2:
        vdveg(1:kproma,jt)=(1._dp/(prahveg(1:kproma)+         &
                                   rbveg(1:kproma,idt_so2)+   &
                                   rsveg(1:kproma,idt_so2)))
        vdsoil(1:kproma,jt)=(1._dp/(prahslsn(1:kproma)+       &
                                    rbslsn(1:kproma,idt_so2)+ &
                                    rsoilso2(1:kproma)))
        vdws(1:kproma,jt)=(1._dp/(prahveg(1:kproma)+          &
                                  rbveg(1:kproma,idt_so2)+    &
                                  rws(idt_so2)))
        vdwat(1:kproma,jt)=1._dp/(prahwat(1:kproma)+rbw(1:kproma,idt_so2)+rwater(idt_so2))
        vdsn(1:kproma,jt)= (1._dp/(prahslsn(1:kproma)+        &
                                   rbslsn(1:kproma,idt_so2)+  &
                                   rsnowhno3so2(1:kproma)))
        vdice(1:kproma,jt)= (1._dp/(prahice(1:kproma)+        &
                                   rbi(1:kproma,idt_so2)+     &
                                   rsnowhno3so2(1:kproma)))
      ELSE    ! aerosol phase
!++mgs 20140519: I doubt that this part of the code will ever be executed: the drydep interface calls this routine only for 
!                gas-phase tracers as far as I understand...
      !--- Use dry deposition velocity for total sulfate for the case of
      !    no partitioning between gas and aerosol phase:
        rsveg(1:kproma,jt)=MAX(10._dp,rsveg(1:kproma,jt))
        !>>SF #458 (replacing WHERE statements)
        vdvfac(1:kproma) = MERGE( &
                                 0.01_dp, &
                                 0.002_dp, &
                                 stheta(1:kproma) > 0.175_dp .AND. psrfl(1:kproma) > 250._dp)
        !<<SF #458 (replacing WHERE statements)

        vdveg(1:kproma,jt)=1._dp/(prahveg(1:kproma)+1._dp/(vdvfac(1:kproma)*pustveg(1:kproma)))
        vdsoil(1:kproma,jt)=vdso4slsn(1:kproma)/100._dp
        vdws(1:kproma,jt)=vdveg(1:kproma,jt)
        vdwat(1:kproma,jt)=vdso4wat(1:kproma)/100._dp
        vdsn(1:kproma,jt)=vdso4slsn(1:kproma)/100._dp
        vdice(1:kproma,jt)=vdso4ice(1:kproma)/100._dp
      END IF

    ELSE    ! all other tracers

      vdveg(1:kproma,jt)=(1._dp/(prahveg(1:kproma)+rbveg(1:kproma,jt)+rsveg(1:kproma,jt)))
      vdsoil(1:kproma,jt)=(1._dp/(prahslsn(1:kproma)+rbslsn(1:kproma,jt)+rsoil(jt)))
      vdws(1:kproma,jt)=(1._dp/(prahveg(1:kproma)+rbveg(1:kproma,jt)+rws(jt)))
 
      ! mz_LG_20020115 for the snow uptake resistance, a minimum of the
      !     resistance of HNO3 and SO2 is being applied to avoid
      !     having very large deposition velocities over ice/snow
      !     covered surfaces. Despite the lack of observations
      !     to support this assumption, it is not very likely that
      !     for other species than HNO3 and SO2 the uptake by ice/
      !     snow surfaces will be larger
 
      vdwat(1:kproma,jt)= 1._dp/(prahwat(1:kproma)+             &
                            rbw(1:kproma,jt)+rwater(jt))
      vdsn(1:kproma,jt)=(1._dp/(prahslsn(1:kproma)+rbslsn(1:kproma,jt)+ &
                                MAX(rsnowhno3so2(1:kproma),rsnow(jt))))
      vdice(1:kproma,jt)=(1._dp/(prahice(1:kproma)+rbi(1:kproma,jt)+    &
                                MAX(rsnowhno3so2(1:kproma),rsnow(jt))))
    END IF

    ! mz_LG_20020115 calculation of the over-land dry deposition velocity, considering
    !     the surface cover fractions

    vdland(1:kproma,jt)=pcvs(1:kproma)*vdsn(1:kproma,jt)+ &
         (1._dp-pcvs(1:kproma))*(1._dp-pcvw(1:kproma))*pvgrat(1:kproma)*vdveg(1:kproma,jt)+ &
         (1._dp-pcvs(1:kproma))*(1._dp-pvgrat(1:kproma))*(1._dp-pcvw(1:kproma))*vdsoil(1:kproma,jt)+ &
         (1._dp-pcvs(1:kproma))*pcvw(1:kproma)*vdws(1:kproma,jt)

    !++mgs: compute stomatal dry dep velocity 
    ! (perhaps limit to ozone only...)
    vdstom(1:kproma,jt)= (1._dp-pcvs(1:kproma))*(1._dp-pcvw(1:kproma))*pvgrat(1:kproma) * &
                         (1._dp/(prahveg(1:kproma)+rbveg(1:kproma,jt)+rstom(1:kproma,jt)))


    fsic(1:kproma) = MERGE(0._dp, pfri(1:kproma), laland(1:kproma))

    ! mz_LG_20020115 finally, calculation of the grid deposition velocity from the
    !     over-land and ocean dry deposition velocities, the deposition
    !     velocities in vdrydep are given in [m s-1]!

    pvd(1:kproma,jt)=vdland(1:kproma,jt)*pfrl(1:kproma)+ &
                     vdwat(1:kproma,jt)*pfrw(1:kproma)+  &
                     vdice(1:kproma,jt)*fsic(1:kproma)

    ! Set pvd to minimum value
!!++mgs 20140519: in the new scheme where also "slow" tracers deposit, this may not be what we want.
!!  pvd(1:kproma,jt) = MAX(1.e-5_dp, pvd(1:kproma,jt))
!!--mgs


!++mgs 20140519
    IF (drydep_ldetail) THEN
      dpzlai(1:kproma,krow) = zlai(1:kproma)
      dpzrmes(1:kproma,krow) = zrmes(1:kproma)
      dpzhc(1:kproma,krow) = zhc(1:kproma)
      dpz0m(1:kproma,krow) = z0m(1:kproma)
      dprahcan(1:kproma,krow) = rahcan(1:kproma)
      dpzrco_leaf(1:kproma,krow) = zrco_leaf(1:kproma)
      dppustveg(1:kproma,krow) = pustveg(1:kproma)
      dppustslsn(1:kproma,krow) = pustslsn(1:kproma)
      dppustarw(1:kproma,krow) = pustarw(1:kproma)
      dppustari(1:kproma,krow) = pustari(1:kproma)
      dppcvs(1:kproma,krow) = pcvs(1:kproma)
      dppcvw(1:kproma,krow) = pcvw(1:kproma)
      dppvgrat(1:kproma,krow) = pvgrat(1:kproma)
      dpprahveg(1:kproma,krow) = prahveg(1:kproma)
      dpprahslsn(1:kproma,krow) = prahslsn(1:kproma)
      dpprahwat(1:kproma,krow) = prahwat(1:kproma)
      dpprahice(1:kproma,krow) = prahice(1:kproma)
      dprsnowhno3so2(1:kproma,krow) = rsnowhno3so2(1:kproma)
      dprleaf(1:kproma, krow) = rleaf(1:kproma,jt)
      dprstom(1:kproma, krow) = rstom(1:kproma,jt)
      dprsveg(1:kproma, krow) = rsveg(1:kproma,jt)
      dprbslsn(1:kproma, krow) = rbslsn(1:kproma,jt)
      dprbveg(1:kproma, krow) = rbveg(1:kproma,jt)
      dpvdveg(1:kproma, krow) = vdveg(1:kproma,jt)
      dpvdsoil(1:kproma, krow) = vdsoil(1:kproma,jt)
      dpvdws(1:kproma, krow) = vdws(1:kproma,jt)
      dpvdsn(1:kproma, krow) = vdsn(1:kproma,jt)
      dpvdland(1:kproma, krow) =  vdland(1:kproma,jt)
      dpvdwat(1:kproma, krow) = vdwat(1:kproma,jt)
      dpvdice(1:kproma, krow) = vdice(1:kproma,jt)
      dprbw(1:kproma, krow) = rbw(1:kproma,jt)
      dprbi(1:kproma, krow) = rbi(1:kproma,jt)
      dprc0x(1:kproma, krow) = rc0x(1:kproma,jt)
    END IF
!--mgs
  END DO      ! jt=1, ntrac

  lfirst = .FALSE.

  END SUBROUTINE drydep_lg_vdbl


  !! ------------------------------------------------------------------------------------------
  !! vd_prep: calculate soil moisture attenuation factor and stomatal resistance
  !! (formerly in vdiff)

  SUBROUTINE vd_prep(kbdim, kproma, krow,              &
                     pws,   pwsmx,  psrfl,             &
                     pfws,  prco_leaf                   )

  USE mo_hammoz_vegetation, ONLY: cva, cvb, cvc, cvbc,     &
                                  cvk, cvkc, cvabc, cvrad 


  !--- arguments
  INTEGER, INTENT(in)   :: kbdim, kproma, krow

  REAL(dp), INTENT(in)  :: pws(kbdim)       ! soil moisture
  REAL(dp), INTENT(in)  :: pwsmx(kbdim)     ! max. moisture content
  REAL(dp), INTENT(in)  :: psrfl(kbdim)     ! surface solar flux 
  REAL(dp), INTENT(out) :: pfws(kbdim)      ! soil moisture attenutation factor
  REAL(dp), INTENT(out) :: prco_leaf(kbdim) ! stomatal resistance

  !--- local variables
  INTEGER       :: jl
  REAL(dp)      :: zwcrit, zwpwp, zqwevap
  REAL(dp)      :: zsrfl, zabcs
  REAL(dp), PARAMETER :: zplmin = 0.35_dp,   &
                         zplmax = 0.75_dp,   &
                         zepsr  = 1.e-10_dp

  !--- 1) Computations
  DO jl=1, kproma

     zwcrit  = zplmax*MAX(1.0E-13_dp,pwsmx(jl))
     zwpwp   = zplmin*MAX(1.0E-13_dp,pwsmx(jl))
     zqwevap = 1._dp/(zwcrit-zwpwp)
     pfws(jl)= MAX(0._dp,MIN(1._dp,(pws(jl)-zwpwp)*zqwevap))

     ! calculation of stomatal resistance
     zsrfl=MAX(zepsr,psrfl(jl)*cvrad)
     zabcs=(cva+cvbc)/(cvc*zsrfl)
     prco_leaf(jl)=1._dp/((cvb*(LOG((zabcs*EXP(cvk)+1._dp)/(zabcs+1._dp))) &
                        /cvabc-(LOG((zabcs+EXP(-cvk))/(zabcs+1._dp))))/cvkc)


  END DO


  END SUBROUTINE vd_prep
 
  !! ------------------------------------------------------------------------------------------
  !! vd_prep_sulfur: calculate specific parameters for sulphate and SO2 deposition
  !! (dependence on soil pH ...)

  SUBROUTINE vd_prep_sulfur(kbdim, kproma, krow,                    &
!>>gf change in argument list (see #78)
                            pustarw, pustari, pustslsn, ptsm, prh,  &
!<<gf
                            psrfl, prahveg, pusfc, pvsfc,           &
!>>gf change in argument list (see #78)
                            rsoilso2, vdso4wat, vdso4ice, vdso4slsn,&
!<<gf
                            stheta                                  )

  USE mo_boundary_condition, ONLY: bc_apply

  !--- arguments
  INTEGER, INTENT(in)   :: kbdim, kproma, krow
  REAL(dp),INTENT(in)   :: pustarw(kbdim),    &
                           pustari(kbdim),    &
                           pustslsn(kbdim),   &
                           ptsm(kbdim),       &     ! temperature (soil?)
                           prh(kbdim),        &     ! relative humidity
                           psrfl(kbdim),      &
                           prahveg(kbdim),    &
                           pusfc(kbdim),      &     ! wind vector u
                           pvsfc(kbdim)             ! wind vector v
  REAL(dp),INTENT(out)  :: rsoilso2(kbdim),   &
                           vdso4wat(kbdim),   &
                           vdso4ice(kbdim),   &     !gf see #78
                           vdso4slsn(kbdim),  &
                           stheta(kbdim)

  !--- local variables
  REAL(dp) :: zsoilph(kbdim,7)

  INTEGER  :: i, jl

  !--- retrieve soilpH values from boundary condition manager
  DO i=1,7 
    CALL bc_apply(ibc_soilph(i), kproma, krow, zsoilph(1:kproma,i))
  END DO

  DO jl=1, kproma
    ! mz_LG_20020115 calculation of standard deviation of wind direction, used
    !     for the calculation of the H2SO4-- deposition velocity, this should
    !     be replaced once by using the Richardson number as criteria to
    !     select between two stability classes relevant to the H2SO4-- deposition
    !     calculations.

    IF (psrfl(jl).GT.10) THEN
      stheta(jl)=MIN(0.5_dp,SQRT(MAX(1.e-5_dp,1._dp/(prahveg(jl)/9._dp+ &
                     SQRT(MAX(0.1_dp,pusfc(jl)**2+pvsfc(jl)**2))))))
    ELSE
      stheta(jl)=MIN(0.5_dp,SQRT(MAX(1.e-5_dp,1._dp/(prahveg(jl)/4._dp+ &
                     SQRT(MAX(0.1_dp,pusfc(jl)**2+pvsfc(jl)**2))))))
    END IF

    ! -- rural continental Vd H2SO4, parameterized deposition velocities
    !    as a function of the friction velocity. The parameterization is
    !    developped using the observed rural continental H2SO4 mass size
    !    distribution [Mehlmann, see Ganzeveld et al., 1998 for reference]

    vdso4slsn(jl)=MAX(1.e-5_dp,0.08_dp-0.57_dp*pustslsn(jl)+ &
                      1.66_dp*pustslsn(jl)**2-0.36_dp*pustslsn(jl)**3)

    ! -- Marine Vd H2SO4 for marine sulfate distribution

!>>gf see #78
    vdso4wat(jl)=MAX(1e-5_dp,0.12_dp-0.24_dp*pustarw(jl)+ &
                     5.25_dp*pustarw(jl)**2-1.84_dp*pustarw(jl)**3)

    vdso4ice(jl)=MAX(1.e-5_dp,0.08_dp-0.57_dp*pustari(jl)+ &
                      1.66_dp*pustari(jl)**2-0.36_dp*pustari(jl)**3)
!<<gf

    ! -- Marine Vd H2SO4 for rural continental H2SO4 mass size distribution (???)
    !      vdso4wat(jl)=MAX(1e-5_dp,0.07_dp+0.1_dp*pustar(jl)+
    !                   *    2.1_dp*pustar(jl)**2-0.20_dp*pustar(jl)**3)

    ! mz_LG_20020115 calculation of soil resistance of SO2 as a function of the RH
    !     and the soil pH. There are five different soil pH classes and
    !     each class has a assigned surface resistance, which is subsequently
    !     corrected for the RH

    !   the soil pH classes are:
    !   SOILPH(JL,1,JR) - soil pH <5.5
    !   SOILPH(JL,2,JR) - soil 5.5 <pH <7.3
    !   SOILPH(JL,3,JR) - soil 7.3< pH <8.5
    !   SOILPH(JL,4,JR) - soil 8.5 <pH
    !   SOILPH(JL,5,JR) - soil 4 < pH <8.5

    rsoilso2(jl)=MAX(50._dp, &
                     zsoilph(jl,1)*115._dp + &
                     zsoilph(jl,2)*65._dp  + &
                     zsoilph(jl,3)*25._dp  + &
                     zsoilph(jl,4)*25._dp  + &
                     zsoilph(jl,5)*70._dp  + &
                     MIN(MAX(0._dp,1000._dp*EXP(MIN(10._dp,269._dp-ptsm(jl)))), &
                     1.e5_dp))

    IF (prh(jl) < 0.6_dp.AND.prh(jl) > 0.4_dp) &  ! semi-arid regions
          rsoilso2(jl)=MAX(50._dp,(rsoilso2(jl)*3.41_dp-85._dp)+ &
                       MIN(MAX(0._dp,1000._dp*EXP(MIN(10._dp,269._dp-ptsm(jl)))), &
                       1.e5_dp))
    IF (prh(jl) <= 0.4_dp) &   ! arid regions
          rsoilso2(jl)=MAX(50._dp,MIN(1000._dp,(rsoilso2(jl)*3.41_dp-85._dp+ &
                       ((0.4_dp-prh(jl))/0.4_dp)*1.e5_dp)+ &
                       MIN(MAX(0._dp,1000._dp*EXP(MIN(10._dp,269._dp-ptsm(jl)))), &
                       1.e5_dp)))

  END DO

  END SUBROUTINE vd_prep_sulfur


END MODULE mo_hammoz_drydep_lg
