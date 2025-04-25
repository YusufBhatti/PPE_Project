!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz
!!
!! \brief
!!  This module contains the mozctl namelist variables, diagnostic streams and
!!  memory management (e.g. boundary conditions) for MOZART. Variables from the 
!!  former mo_moz_constants are also incldued here.
!!
!! \author Martin G. Schultz (FZ-Juelich)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - MGS: original version (2009-08-27)
!!
!! \limitations
!!  none
!!
!! \details
!!
!! \bibliographic_references
!!  none
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Note on boundary conditions for moz:
!   bc_aerosol: aerosol climatology from MOZART
!   bc_ch4:     methane concentration field as lower bc (separate from general lbc, because
!               input file is lat-time and relaxation approach is used.)
!   uv-albedo:  provide filename in namelist; rest done in mo_moz_uvalbedo
!   lbc:        tropospheric lower bc: provide BC structure for filename and settings, use
!               lbc_species for names of species
!   ubc:        stratospheric upper bc: provide BC structure for filename and settings, use
!               ubc_species for names of species
!


MODULE mo_moz

  USE mo_kind,               ONLY: dp
  USE mo_moz_mods,           ONLY: plon, plat, pcnstm1, grpcnt    ! preprocessor generated!
  USE mo_boundary_condition, ONLY: bc_nml
  USE mo_tracdef,            ONLY: jptrac
#ifdef __PGI  
  USE mo_decomposition,      ONLY : ldc => local_decomposition
#endif


  IMPLICIT NONE

  PRIVATE

  SAVE

  ! module routines
  PUBLIC :: setmoz

  !------------------------------------------------------------------------------------------------------
  ! mozctl namelist
  !-- flags controlling MOZART processes
  LOGICAL, PUBLIC            :: lchemsolv,       & ! activate chemistry solver ("class solution algorithms")
                                lphotolysis,     & ! activate photolysis calculation
                                lfastj,          & ! use fastJ-x algorithm for photolysis rates
                                lfastjaero,      & ! use fastJ-x with aerosols (on/off)
                                lfastjcloud,     & ! use fastJ-x with clouds (on/off)
                                lstrathet,       & ! compute heterogeneous stratospheric chemistry
                                ltrophet           ! compute heterogeneous tropospheric chemistry
           
  ! ### NOTE: ndrydep and nwetdep choices not fully implemented yet! 
  INTEGER, PUBLIC            :: ndrydep            ! choice of dry deposition scheme:
                                                   ! 0=OFF, 1=fixed vd, 2=interactive (Ganzeveld)
  INTEGER, PUBLIC            :: nwetdep            ! choice of wet deposition scheme:
                                                   ! 0=OFF, 1=ON    (????)

  !-- species lists for boundary conditions
  CHARACTER(len=32), PUBLIC   :: lbc_species(pcnstm1)      ! lower boundary conditions
  CHARACTER(len=32), PUBLIC   :: ubc_species(pcnstm1)      ! upper boundary conditions

  !-- species lists for output and diagnostics
  CHARACTER(len=32), PUBLIC   :: out_species(pcnstm1)     ! species which shall appear in _tracer stream
  CHARACTER(len=32), PUBLIC   :: burden_species(pcnstm1)  ! species with burden diagnostics
  CHARACTER(len=32), PUBLIC   :: budget_species(pcnstm1)  ! species with budget diagnostics ### TO BE IMPLEMENTED!

  !-- variable names for photolysis output
  CHARACTER(len=32), PUBLIC  :: photovars(jptrac)         

  !-- file name for UV albedo data (green and white)
  CHARACTER(len=256), PUBLIC :: uvalbedo_file

  !-- name list fine control of boundary conditions
  TYPE(bc_nml), PUBLIC       :: bc_sad,                   &
                                bc_ch4,                   &
                                bc_lbc,                   &
                                bc_ubc                     


  !------------------------------------------------------------------------------------------------------
  ! other flags and parameters (no namelist control)
  LOGICAL, PUBLIC            :: lhammoniagases          ! indicates that simulation has O2, CO2, N2O and CH4

  ! MOZART species names (mo_moz_subs)
  CHARACTER(len=24), PUBLIC  :: tracnam(pcnstm1)        ! MOZART species names
  CHARACTER(len=24), PUBLIC  :: natsnam(max(1,grpcnt))  ! names of non-advected trace species (grouping)

  ! tracer indices 
  INTEGER, PUBLIC            :: idt_o3,                   & ! ozone
                                e2m(jptrac) = -1            ! mapping ECHAM -> MOZART

  ! code management stuff
  INTEGER, PUBLIC            :: nlbc                        ! number of lower boundary conditions for tracers
  INTEGER, PUBLIC            :: nubc                        ! number of upper boundary conditions for tracers
  INTEGER, PUBLIC            :: noutspec                    ! number of tracers in output
  INTEGER, PUBLIC            :: nburdenspec                 ! number of tracers in burden diagnostics
  INTEGER, PUBLIC            :: nbudgetspec                 ! number of tracers in budget diagnostics
  INTEGER, PUBLIC            :: nphotovars                  ! number of variables in photolysis output

  ! boundary condition indices
  INTEGER, PUBLIC            :: ibc_sad,                  & ! index for SAD boundary condition
                                ibc_ch4                     ! index for methane boundary condition

  !------------------------------------------------------------------------------------------------------
  ! mo_moz_constants
!++mgs: changed to ECHAM5 value (see mo_constants.f90)
! REAL(dp), PARAMETER, PUBLIC ::  gravit = 9.80616_dp       ! m/s
  REAL(dp), PARAMETER, PUBLIC ::  gravit = 9.80665_dp       ! m/s
!--mgs
  REAL(dp), PARAMETER, PUBLIC ::  rgrav  = 1._dp/gravit
  REAL(dp), PARAMETER, PUBLIC ::  dayspy = 365._dp          ! days per year
  REAL(dp), PARAMETER, PUBLIC ::  rearth = 6.37122e6_dp     ! radius earth (m)
  REAL(dp), PARAMETER, PUBLIC ::  aomega = .7292e-4_dp      ! angular velocity of Earth rotation (1/s)

  REAL(dp), PUBLIC ::  pi                                   ! radians
  REAL(dp), PUBLIC ::  twopi                                ! 2*pi (radians)
  REAL(dp), PUBLIC ::  pid2                                 ! pi/2 (radians)
  REAL(dp), PUBLIC ::  r2d                                  ! radians to degrees
  REAL(dp), PUBLIC ::  d2r                                  ! degrees to radians
  REAL(dp), PUBLIC ::  lat25 = 0._dp                        ! 25 latitude (radians)
  REAL(dp), PUBLIC ::  lat45 = 0._dp                        ! 45 latitude (radians)
  REAL(dp), PUBLIC ::  lat59 = 0._dp                        ! 59 latitude (radians)
  REAL(dp), PUBLIC ::  lat60 = 0._dp                        ! 60 latitude (radians)
  REAL(dp), PUBLIC ::  lat70 = 0._dp                        ! 70 latitude (radians)

  !------------------------------------------------------------------------------------------------------
  ! solver diagnostics (mo_moz_imp_sol)
  PUBLIC :: t_pdiag,                  & ! type structure for diagnostics pointer
            pdiags                      ! diagnostic field

  !------------------------------------------------------------------------------------------------------
  ! Variable declarations

  type t_pdiag                             ! routine specific print diagnostics
     logical :: adv
     logical :: physlic
     logical :: imp_slv
     logical :: negtrc
  end type t_pdiag

  type(t_pdiag) :: pdiags
  
  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! setmoz: set defaults and parse mozctl namelist; set boundary conditions
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Externals:
!! <ol>
!! <li>none
!! </ol>
!!
!! @par Notes
!! 
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE setmoz (lchemistry, bc_aerosol)
  USE mo_namelist,           ONLY: position_nml, POSITIONED, open_nml
  USE mo_exception,          ONLY: message, em_param, em_error, message_text
  USE mo_submodel,           ONLY: print_value, lham
  USE mo_util_string,        ONLY: separator
  USE mo_mpi,                ONLY: p_parallel, p_parallel_io, p_bcast, p_io, p_pe
#ifndef __PGI  
  USE mo_decomposition,      ONLY: ldc => local_decomposition
#endif
  USE mo_moz_mods,           ONLY: plon, plonl, plat, pplon, plev, plevp, plevm, plnplv, nodes
  USE mo_boundary_condition, ONLY: bc_define, p_bcast_bc

!++mgs #249: added lchemistry
  LOGICAL, INTENT(in)           :: lchemistry
!--mgs
  TYPE(bc_nml), INTENT(inout)   :: bc_aerosol

  INTEGER                       :: ii, ierr, idum, inml, iunit, emlevel

  CHARACTER(LEN=64)             :: text


  INCLUDE 'mozctl.inc'

  !------------------------------------------------------------
  ! 0) initialize
  ! mo_moz_constants values
  !
  pi     = 4._dp*ATAN( 1._dp )
  twopi  = 2._dp*pi
  pid2   = .5_dp*pi
  d2r    = pi/180._dp
  r2d    = 180._dp/pi

  lat25 = d2r * 25._dp
  lat45 = d2r * 45._dp
  lat59 = d2r * 59._dp
  lat60 = d2r * 60._dp
  lat70 = d2r * 70._dp

  ! MOZART grid definition 
  !
  plon   = ldc%nlon    
  plonl  = ldc%nproma  
  plat   = ldc%nlat    
  pplon  = 1
  plev   = ldc%nlev
  plevp  = plev + 1
  plevm  = plev - 1
  plnplv = plonl*plev
  nodes  = -1      ! just for safety - shouldn't be used.

  ! counters and indices
  !
  nlbc         = -1
  nubc         = -1
  noutspec     = -1
  nburdenspec  = -1
  nbudgetspec  = -1
  nphotovars   = -1
  ibc_sad      = -1
  ibc_ch4      = -1

  ! set tracer indices to 0
  idt_o3 = 0

  ! mozctl parameters
  !
  ! switches and process selection
!++mgs #249: improved flag settings
  lchemsolv    = lchemistry
  lphotolysis  = lchemistry
!--mgs
  lfastj       = .false.     ! ## perhaps later = .true. ???
  lfastjaero   = .false.     ! ## perhaps later = .true. ???
  lfastjcloud  = .false.     ! ## perhaps later = .true. ???
  lstrathet    = .true.
  ltrophet     = .true.
            
  ndrydep       = 2      ! interactive scheme as default
  nwetdep       = 1      ! on per default

  ! boundary conditions, output and diagnostics
  budget_species(:)    = ''
  out_species(:)       = ''
  out_species(1)       = 'default'
  burden_species(:)    = ''
  burden_species(1)    = 'default'
  ubc_species(:)       = ''
  lbc_species(:)       = ''
  lbc_species(1)       = 'default'     ! (see mo_moz_lbc_ubc)
  photovars(:)         = ''
  photovars(1)         = 'default'     ! (see mo_moz_photo)

  uvalbedo_file  = 'moz_uvalbedo.%T0.nc'    ! %T0 shall be replaced by hor. resolution

  ! chemical solver diagnostics (imp_sol)
  !
  pdiags%negtrc  = .TRUE. 
  pdiags%imp_slv = .FALSE.
  pdiags%adv     = .FALSE.
  pdiags%physlic = .FALSE.

  lhammoniagases = .FALSE.   ! will be diagnosed in moz_initialize

  ! default settings for boundary conditions
  !
  ! bc_nml structure for stratospheric aerosol density climatology (default: EF_FILE)
  bc_sad%ef_type = 2          ! ef_type: EF_INACTIVE=0, EF_VALUE=1, EF_FILE=2, EF_MODULE=3
  bc_sad%ef_template = 'moz_sad_sulf.%T0.nc'      ! %T0 shall be replaced with hor resolution
  bc_sad%ef_varname = "sad_sage"                  ! variable name in file
  bc_sad%ef_geometry = 3      ! ef_geometry: EF_3D=1 , EF_LONLAT=2 , EF_LATLEV=3 , EF_LEV=4 , EF_LAT=5 , EF_SINGLE=6
  bc_sad%ef_timedef = 1       ! ef_timedef: EF_TIMERESOLVED=1 , EF_IGNOREYEAR=2 , EF_CONSTANT=3
  bc_sad%ef_interpolate = 0   ! none    ### may want to do time interp. later! ###
  bc_sad%ef_actual_unit = 'cm2/cm3'
  bc_sad%bc_domain = 0        ! apply everywhere
  bc_sad%bc_mode = 1          ! bc_mode: BC_REPLACE=1, BC_ADD=2, BC_RELAX=3, BC_SPECIAL=4

  ! bc_nml structure for generic lower boundary condition (default: EF_FILE)
  bc_lbc%ef_type = 2          ! ef_type: EF_INACTIVE=0, EF_VALUE=1, EF_FILE=2, EF_MODULE=3
  bc_lbc%ef_template = 'moz_lbc.%T0.nc'    ! %T0%L0 shall be replaced by hor. resolution
  bc_lbc%ef_varname = "*"     ! variable name in file (to be set in loop)
  bc_lbc%ef_geometry = 0      ! ef_geometry: EF_3D=1 , EF_LONLAT=2 , EF_LATLEV=3 , EF_LEV=4 , EF_LAT=5 , EF_SINGLE=6
  bc_lbc%ef_timedef = 1       ! ef_timedef: EF_TIMERESOLVED=1 , EF_IGNOREYEAR=2 , EF_CONSTANT=3
  bc_lbc%ef_interpolate = 0   ! none
  bc_lbc%ef_actual_unit = 'mole mole-1'
  bc_lbc%bc_domain = 1        ! bc_domain: BC_EVERYWHERE=0, BC_BOTTOM=1, BC_TOP=2, BC_LEVEL=3, BC_ALTITUDE=4, BC_PRESSURE=5
  bc_lbc%bc_mode = 1          ! bc_mode: BC_REPLACE=1, BC_ADD=2, BC_RELAX=3, BC_SPECIAL=4

  ! bc_nml structure for generic upper boundary condition (default: OFF)
  bc_ubc%ef_type = 0          ! ef_type: EF_INACTIVE=0, EF_VALUE=1, EF_FILE=2, EF_MODULE=3
  bc_ubc%ef_template = 'moz_ubc.%T0%L0.nc'    ! %T0%L0 shall be replaced by hor and vert. resolution
  bc_ubc%ef_varname = "*"     ! variable name in file (to be set in loop)
  bc_ubc%ef_geometry = 3      ! ef_geometry: EF_3D=1 , EF_LONLAT=2 , EF_LATLEV=3 , EF_LEV=4 , EF_LAT=5 , EF_SINGLE=6
  bc_ubc%ef_timedef = 2       ! ef_timedef: EF_TIMERESOLVED=1 , EF_IGNOREYEAR=2 , EF_CONSTANT=3
  bc_ubc%ef_interpolate = 0   ! none
  bc_ubc%ef_actual_unit = 'VMR'
  bc_ubc%bc_domain = 3        ! bc_domain: BC_EVERYWHERE=0, BC_BOTTOM=1, BC_TOP=2, BC_LEVEL=3, BC_ALTITUDE=4, BC_PRESSURE=5
  bc_ubc%bc_mode = 3          ! bc_mode: BC_REPLACE=1, BC_ADD=2, BC_RELAX=3, BC_SPECIAL=4
  bc_ubc%bc_minlev = 2
  bc_ubc%bc_maxlev = 10
  bc_ubc%bc_relaxtime = 10._dp * 86400._dp ! relaxation time

  ! template for boundary condition (default except variable names)
  bc_aerosol%ef_type = 2          ! ef_type: EF_INACTIVE=0, EF_VALUE=1, EF_FILE=2, EF_MODULE=3
  bc_aerosol%ef_template = 'ham_aerosol_climatology.%T0%L0.nc'      ! %T0%L0 shall be replaced by hor and vert. resolution
  bc_aerosol%ef_varname = "*******"                  ! variable name in file
  bc_aerosol%ef_geometry = 1      ! ef_geometry: EF_3D=1 , EF_LONLAT=2 , EF_LATLEV=3 , EF_LEV=4 , EF_LAT=5 , EF_SINGLE=6
  bc_aerosol%ef_timedef = 2       ! ef_timedef: EF_TIMERESOLVED=1 , EF_IGNOREYEAR=2 , EF_CONSTANT=3
  bc_aerosol%ef_interpolate = 0   ! none
  bc_aerosol%ef_actual_unit = '*******'
  bc_aerosol%bc_domain = 0        ! apply everywhere
  bc_aerosol%bc_mode = 1          ! bc_mode: BC_REPLACE=1, BC_ADD=2, BC_RELAX=3, BC_SPECIAL=4

!++mgs #249: added deactivation of boundary conditions if lchemistry==false
  IF (.NOT. lchemistry) THEN
    bc_aerosol%ef_type = 0     ! inactive
    bc_sad%ef_type = 0         ! inactive
    ! don't touch lbc and ubc, because these can be run in a pure transport simulation!
  END IF
!--mgs

  !------------------------------------------------------------
  ! 1) read mozctl namelist
  !
  IF (p_parallel_io) THEN
    inml = open_nml('namelist.echam')
    iunit =  position_nml ('MOZCTL', inml, status=ierr)
    SELECT CASE (ierr)
    CASE (POSITIONED)
      READ(iunit, mozctl)
    END SELECT
  END IF

!++mgs #249: overwrite lphotolysis, lchemsolv, lstrathet, and ltrophet if lchemistry==false
  IF (.NOT. lchemistry) THEN
    lchemsolv = .false.
    lphotolysis = .false.
    lstrathet = .false.
    ltrophet  = .false.
  END IF
!--mgs

  IF (p_parallel) THEN
    CALL p_bcast (lchemsolv, p_io)
    CALL p_bcast (lphotolysis, p_io)
    CALL p_bcast (lfastj, p_io)
    CALL p_bcast (lfastjaero, p_io)
    CALL p_bcast (lfastjcloud, p_io)
    CALL p_bcast (lstrathet, p_io)
    CALL p_bcast (ltrophet, p_io)
    CALL p_bcast (ndrydep, p_io)
    CALL p_bcast (nwetdep, p_io)
    CALL p_bcast (lbc_species, p_io)
    CALL p_bcast (ubc_species, p_io)
    CALL p_bcast (out_species, p_io)
    CALL p_bcast (burden_species, p_io)
    CALL p_bcast (budget_species, p_io)
    CALL p_bcast (photovars, p_io)
    CALL p_bcast (uvalbedo_file, p_io)
    CALL p_bcast_bc (bc_sad, p_io)
    CALL p_bcast_bc (bc_ch4, p_io)
    CALL p_bcast_bc (bc_lbc, p_io)
    CALL p_bcast_bc (bc_ubc, p_io)
    CALL p_bcast_bc (bc_aerosol, p_io)
  END IF
  ! switch SAD off if lstrathet == false
  IF (.NOT. lstrathet) bc_sad%ef_type = 0

  ! determine number of tracer name entries
  DO ii = 1, pcnstm1
    IF ( TRIM(budget_species(ii) ) == '' .AND. nbudgetspec < 0 ) nbudgetspec = ii-1
    IF ( TRIM(out_species(ii) ) == '' .AND. noutspec < 0 ) noutspec = ii-1
    IF ( TRIM(burden_species(ii) ) == '' .AND. nburdenspec < 0 ) nburdenspec = ii-1
    IF ( TRIM(ubc_species(ii) ) == '' .AND. nubc < 0 ) nubc = ii-1
    if ( TRIM(lbc_species(ii) ) == '' .AND. nlbc < 0 ) nlbc = ii-1
  END DO
  ! determine number of variables in (photo) output
  ! see mo_moz_photo. List will be expanded automatically to include all 
  ! variables if the first entry is "all"
  DO ii = 1, jptrac
    if ( TRIM(photovars(ii) ) == '' .AND. nphotovars < 0 ) nphotovars = ii-1
  END DO

!>>SF
  !--- Security:
  IF (ltrophet) THEN
     !SF note: the following will disapear when het chem is made independent of HAM
     !         and adapted for all aerosol microphysics schemes
     IF (.NOT. lham) THEN
        CALL message('setmoz','Het. tropos. chemistry requested but HAM is .FALSE.', level=em_error)
     ENDIF
  ENDIF
!<<SF

  !--- report settings
  CALL message('', separator)
  CALL message('setmoz', 'Parameter settings for MOZ module:')
  CALL print_value('Chemistry solver active (lchemsolv)', lchemsolv)
  CALL print_value('Photolysis calculation active (lphotolysis)', lphotolysis)
!++mgs #249: added lphotolysis condition
  IF (lphotolysis) THEN
     IF (lfastj) THEN
        CALL message('', 'Using fastJ-x photolysis scheme', level=em_param)
        CALL print_value('aerosols in FastJ-x (lfastjaero)', lfastjaero)
        CALL print_value('clouds in FastJ-x (lfastjcloud)', lfastjcloud)
      ELSE
        CALL message('', 'Using WACCM photolysis scheme', level=em_param)
     END IF
  END IF
!--mgs
  CALL print_value('Heterogeneous reactions in stratosphere active (lstrathet)', lstrathet)
  SELECT CASE (ndrydep)
    CASE (0)
      text="OFF"
      emlevel=em_param
    CASE (1)
      text="prescribed velocities (not implemented)"
      emlevel=em_error
    CASE (2)
      text="Ganzeveld scheme"
      emlevel=em_param
    CASE DEFAULT
      text="not defined!"
      emlevel=em_error
  END SELECT
  CALL print_value('Heterogeneous reactions in troposphere active (ltrophet)', ltrophet)
  CALL message ('','Scheme selection for dry deposition (ndrydep) : '//TRIM(text), level=emlevel)
  SELECT CASE (nwetdep)
    CASE (0)
      text="OFF"
      emlevel=em_param
    CASE (1)
      text="HAMMOZ wetdep scheme"
      emlevel=em_param
    CASE DEFAULT
      text="not defined!"
      emlevel=em_error
  END SELECT
  CALL message ('','Scheme selection for wet deposition (nwetdep) : '//TRIM(text), level=emlevel)
  CALL message('', separator)
!### double-check default values at a later time. Idea: set bc_XXX%ef_type from 0 to 2 in namelist and GO! ###
  ! bc_nml structure for sulfate climatology (default: OFF)
!
! 2) define boundary conditions
!

  ibc_sad     = bc_define('Stratospheric aerosol densities', bc_sad, 3, .TRUE.)
  ibc_ch4     = bc_define('CH4 lower boundary', bc_ch4, 3, .TRUE.)

END SUBROUTINE setmoz

END MODULE mo_moz
