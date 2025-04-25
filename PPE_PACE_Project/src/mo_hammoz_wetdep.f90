!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_wetdep.f90
!!
!! \brief
!! module to interface ECHAM submodules with wet deposition module(s)
!!
!! \author M. Schultz   (FZ Juelich)
!! \author S. Ferrachat (ETH-Zuerich)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2009-10-02)
!!   -# S. Ferrachat (ETH-Zuerich) - revision (2009-12-16)
!!
!! \limitations
!! None
!!
!! \details
!! This module initializes the scheme based on the namelist parameters
!! in submodeldiagctl and creates a stream for variable pointers and 
!! diagnostic quantities used in the wet deposition scheme. It also
!! provides a generic interface to the actual wet deposition routine(s).
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
MODULE mo_hammoz_wetdep

  USE mo_kind,             ONLY: dp
  USE mo_submodel_diag,    ONLY: t_diag_list
  USE mo_ham,              ONLY: nclass

  IMPLICIT NONE

  PRIVATE

  ! public variables  (see declaration below)

  ! helper flags for use in the computational routines
  PUBLIC                       :: lwetdepdetail

  ! subprograms
  PUBLIC                       :: init_wetdep_stream, &
                                  wetdep_interface

  ! wetdep_stream
  INTEGER, PARAMETER           :: nwetdepvars=24
  CHARACTER(LEN=32)            :: wetdepvars(1:nwetdepvars)= &
                                (/'precipform       ', &   ! precip. formation rate in kg m-2 s-1 (3D)
                                  'precipevap       ', &   ! precip. evaporation rate in kg m-2 s-1 (3D)
                                  'uparfrac         ', &   ! updraft grid box fraction (3D)
                                  'wdep             ', &   ! total wetdep flux
                                  'wdep_conv        ', &   ! wet dep. in convective clouds
                                  'wdep_strat       ', &   ! wet dep. in stratiform clouds
                                  'wdep_incl        ', &   ! total wet dep. in clouds
                                  'wdep_blcl        ', &   ! wet dep. below (convective(?) clouds
                                  'wdep_incl_swn    ', &   ! .. in cloud stratiform warm nucleation
                                  'wdep_incl_swi    ', &   ! .. in cloud stratiform warm impaction
                                  'wdep_incl_smn    ', &   ! .. in cloud stratiform mixed nucleation
                                  'wdep_incl_smi    ', &   ! .. in cloud stratiform mixed impaction
                                  'wdep_incl_scn    ', &   ! .. in cloud stratiform cold nucleation
                                  'wdep_incl_sci    ', &   ! .. in cloud stratiform cold impaction
                                  'wdep_incl_cwn    ', &   ! .. in cloud convective warm nucleation
                                  'wdep_incl_cwi    ', &   ! .. in cloud convective warm impaction
                                  'wdep_incl_cmn    ', &   ! .. in cloud convective mixed nucleation
                                  'wdep_incl_cmi    ', &   ! .. in cloud convective mixed impaction
                                  'wdep_incl_ccn    ', &   ! .. in cloud convective cold nucleation
                                  'wdep_incl_cci    ', &   ! .. in cloud convective cold impaction
                                  'wdep_blcl_sr     ', &   ! .. below cloud stratiform rain
                                  'wdep_blcl_ss     ', &   ! .. below cloud stratiform snow
                                  'wdep_blcl_cr     ', &   ! .. below cloud convective rain
                                  'wdep_blcl_cs     '   /) ! .. below cloud convective snow


  ! variable pointers and diagnostic lists
  REAL(dp), POINTER,  PUBLIC   :: precipform(:,:,:),   precipevap(:,:,:),      &
                                  uparfrac(:,:,:)
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep         ! wet deposition flux
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_conv    ! .. in convective precip.
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_strat   ! .. in stratiform clouds
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl    ! .. total in clouds
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_blcl    ! .. below clouds
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_blcl_sr ! detailed diagnostics
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_blcl_ss
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_blcl_cr
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_blcl_cs
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_swn
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_swi
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_smn
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_smi
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_scn
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_sci
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_cwn
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_cwi
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_cmn
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_cmi
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_ccn
  TYPE (t_diag_list), PUBLIC, TARGET   :: wdep_incl_cci

  LOGICAL :: lwetdepdetail  ! produce detailed wetdep diagnostics

  CONTAINS

  SUBROUTINE init_wetdep_stream

    USE mo_string_utls,         ONLY: st1_in_st2_proof
    USE mo_util_string,         ONLY: tolower
    USE mo_exception,           ONLY: finish
    USE mo_memory_base,         ONLY: t_stream, new_stream,   &
                                      default_stream_setting, &
                                      add_stream_reference,   &
                                      add_stream_element,     &
                                      AUTO
    USE mo_ham_m7_trac,         ONLY: ham_get_class_flag
    USE mo_tracer,              ONLY: validate_traclist
    USE mo_tracdef,             ONLY: ln, ntrac, trlist, GAS, AEROSOL
    USE mo_species,             ONLY: nspec, speclist
    USE mo_submodel_streams,    ONLY: wetdep_lpost, wetdep_tinterval, wetdepnam,   &
                                      wetdep_gastrac, wetdep_keytype
    USE mo_submodel_diag,       ONLY: new_diag_list, new_diag,                     &
                                      BYTRACER, BYSPECIES, BYMODE, BYNUMMODE
    USE mo_submodel,            ONLY: lham !SF, see #228

    ! local variables
    INTEGER, PARAMETER             :: ndefault = 1
    CHARACTER(LEN=32)              :: defnam(1:ndefault)   = &   ! default output
                                (/ 'wdep             ' /)   ! total wetdep flux
 
    INTEGER, PARAMETER             :: nall     = 8
    CHARACTER(LEN=32)              :: allnam(1:nall)   = &       ! output for ALL
                                (/'precipform       ', &   ! precip. formation rate in kg m-2 s-1 (3D)
                                  'precipevap       ', &   ! precip. evaporation rate in kg m-2 s-1 (3D)
                                  'uparfrac         ', &   ! updraft grid box fraction (3D)
                                  'wdep             ', &   ! total wetdep flux
                                  'wdep_conv        ', &   ! wet dep. in convective clouds
                                  'wdep_strat       ', &   ! wet dep. in stratiform clouds
                                  'wdep_incl        ', &   ! total wet dep. in clouds
                                  'wdep_blcl        '  /)  ! wet dep. below (convective(?) clouds

    CHARACTER(len=ln)              :: defaultgas(3)        = &   ! default gas-phase tracers for diagnost
                                (/ 'SO2     ',               &
                                   'H2SO4   ',               &
                                   'HNO3    '          /)

    LOGICAL                        :: tracflag(ntrac), specflag(nspec), modflag(MAX(nclass,1))
    CHARACTER(LEN=ln)              :: tracname(ntrac), specname(nspec), modname(MAX(nclass,1)), &
                                      modnumname(MAX(nclass,1)) !SF #299
    CHARACTER(LEN=32)              :: cdiagname
    CHARACTER(LEN=64)              :: clongname
    TYPE (t_stream), POINTER       :: swetdep
    INTEGER                        :: ierr, jt, idiag
    LOGICAL                        :: lpost
    TYPE (t_diag_list), POINTER    :: ptrwdep         ! to avoid code repetition

    !++mgs: default values and namelist read are done in init_submodel_streams !

    !-- initialize output control variables. Will be set automatically!
    lwetdepdetail = .FALSE.

    !-- handle ALL, DETAIL and DEFAULT options for wetdep output variables
    !-- Note: ALL means "all you would normally want", use DETAIL to get everything
    IF (TRIM(tolower(wetdepnam(1))) == 'detail')  THEN
       wetdepnam(1:nwetdepvars) = wetdepvars(:)
       lwetdepdetail = .TRUE. !SF #567
    ENDIF
    IF (TRIM(tolower(wetdepnam(1))) == 'all')     wetdepnam(1:nall) = allnam(:)
    IF (TRIM(tolower(wetdepnam(1))) == 'default') wetdepnam(1:ndefault) = defnam(:)

    !-- check that all variable names from namelist are valid
    IF (.NOT. st1_in_st2_proof( wetdepnam, wetdepvars, ierr=ierr) ) THEN
      IF (ierr > 0) CALL finish ( 'ini_wetdep_stream', 'variable '// &
                                  wetdepnam(ierr)//' does not exist in wetdep stream' )
    END IF

    !-- find out which gas-phase tracers shall be included in diagnostics
    CALL validate_traclist(wetdep_gastrac, defaultgas, nphase=GAS,              &
                           lwetdep=.true.)                     !>>dod bugfix <<<dod

    !-- define the flags and names for the diagnostic lists. We need one set of flags and
    !   names for each key_type (BYTRACER, BYSPECIES, BYMODE)
    !   gas-phase tracers will always be defined BYTRACER, for aerosol tracers one of the
    !   following lists will be empty.
    !   Note: wdep uses BYTRACER, BYSPECIES or BYMODE
    tracflag(:) = .FALSE.
    DO jt = 1,ntrac
      tracname(jt) = trlist%ti(jt)%fullname
      IF (trlist%ti(jt)%nphase == GAS) THEN
        tracflag(jt) = st1_in_st2_proof(trlist%ti(jt)%fullname, wetdep_gastrac)
      ELSE     ! aerosol tracer
        IF (wetdep_keytype == BYTRACER .AND. nclass > 0) THEN
          tracflag(jt) = trlist%ti(jt)%nwetdep > 0
        END IF
      END IF
    END DO
    specflag(:) = .FALSE.
    DO jt = 1,nspec
      specname(jt) = speclist(jt)%shortname
      IF (wetdep_keytype == BYSPECIES .AND.             &
          IAND(speclist(jt)%nphase, AEROSOL) /= 0 .AND. & !>>dod soa: removed check of trtype <<dod
          nclass > 0) THEN
        specflag(jt) = speclist(jt)%lwetdep
      END IF
    END DO
    modflag(:) = .FALSE.
    modname(:) = ''

    !SF #228, adding a condition to check that HAM is active:
    !SF #299, adding a condition to check if BYMODE is relevant:
    IF (lham .AND. nclass > 0 .AND. (wetdep_keytype == BYMODE)) &
       CALL ham_get_class_flag(nclass, modflag, modname, modnumname, lwetdep=.true.)

    !-- open new diagnostic stream
    CALL new_stream (swetdep,'wetdep',lpost=wetdep_lpost,lrerun=.FALSE., &
         interval=wetdep_tinterval)
    CALL default_stream_setting (swetdep, lrerun = .FALSE., &
         contnorest = .TRUE., table = 199, &
         laccu = .false., code = AUTO)
   
    !-- add standard ECHAM variables
    IF (wetdep_lpost) THEN
      CALL add_stream_reference (swetdep, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
      CALL add_stream_reference (swetdep, 'lsp'     ,'sp'    ,lpost=.TRUE.)
      CALL add_stream_reference (swetdep, 'aps'     ,'g3b'   ,lpost=.TRUE.)
      CALL add_stream_reference (swetdep, 'gboxarea','geoloc',lpost=.TRUE.)
    END IF

    !-- add physical variables to stream
    lpost = st1_in_st2_proof( 'precipform', wetdepnam)
    CALL add_stream_element (swetdep, 'precipform', precipform, &
         longname = 'precipitation formation rate', &
         units = 'kg m-2 s-1', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'precipevap', wetdepnam)
    CALL add_stream_element (swetdep, 'precipevap', precipevap, &
         longname = 'precipitation evaporation rate', &
         units = 'kg m-2 s-1', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'uparfrac', wetdepnam)
    CALL add_stream_element (swetdep, 'uparfrac', uparfrac, &
         longname = 'fraction of grid cell covered by convective clouds', &
         units = '1', lpost = lpost)
   
    !-- pure diagnostic quantities (averaged)

    CALL default_stream_setting (swetdep, lrerun=.FALSE., laccu=.TRUE.)

    ! add total wetdep flux and detailed diagnostics
    DO idiag = 1,21
      SELECT CASE (idiag)
      CASE (1) 
        ptrwdep => wdep
        cdiagname = 'wdep'
        clongname = 'wet deposition flux total'
      CASE (2) 
        ptrwdep => wdep_conv
        cdiagname = 'wdep_conv'
        clongname = 'wet deposition flux in convective clouds'
      CASE (3)
        ptrwdep => wdep_strat
        cdiagname = 'wdep_strat'
        clongname = 'wet deposition flux in stratiform clouds'
      CASE (4)
        ptrwdep => wdep_incl
        cdiagname = 'wdep_incl'
        clongname = 'wet deposition flux in-cloud'
      CASE (5)
        ptrwdep => wdep_blcl
        cdiagname = 'wdep_blcl'
        clongname = 'wet deposition flux below-cloud'
      CASE (6)
        ptrwdep => wdep_blcl_sr
        cdiagname = 'wdep_blcl_sr'
        clongname = 'wet deposition flux below-cloud stratiform rain'
      CASE (7)
        ptrwdep => wdep_blcl_ss
        cdiagname = 'wdep_blcl_ss'
        clongname = 'wet deposition flux below-cloud stratiform snow'
      CASE (8)
        ptrwdep => wdep_blcl_cr
        cdiagname = 'wdep_blcl_cr'
        clongname = 'wet deposition flux below-cloud convective rain'
      CASE (9)
        ptrwdep => wdep_blcl_cs
        cdiagname = 'wdep_blcl_cs'
        clongname = 'wet deposition flux below-cloud convective snow'
      CASE (10)
        ptrwdep => wdep_incl_swn
        cdiagname = 'wdep_incl_swn'
        clongname = 'wet deposition flux in-cloud stratiform warm nucleation'
      CASE (11)
        ptrwdep => wdep_incl_swi
        cdiagname = 'wdep_incl_swi'
        clongname = 'wet deposition flux in-cloud stratiform warm impaction'
      CASE (12)
        ptrwdep => wdep_incl_smn
        cdiagname = 'wdep_incl_smn'
        clongname = 'wet deposition flux in-cloud stratiform mixed nucleation'
      CASE (13)
        ptrwdep => wdep_incl_smi
        cdiagname = 'wdep_incl_smi'
        clongname = 'wet deposition flux in-cloud stratiform mixed impaction'
      CASE (14)
        ptrwdep => wdep_incl_scn
        cdiagname = 'wdep_incl_scn'
        clongname = 'wet deposition flux in-cloud stratiform cold nucleation'
      CASE (15)
        ptrwdep => wdep_incl_sci
        cdiagname = 'wdep_incl_sci'
        clongname = 'wet deposition flux in-cloud stratiform cold impaction'
      CASE (16)
        ptrwdep => wdep_incl_cwn
        cdiagname = 'wdep_incl_cwn'
        clongname = 'wet deposition flux in-cloud convective warm nucleation'
      CASE (17)
        ptrwdep => wdep_incl_cwi
        cdiagname = 'wdep_incl_cwi'
        clongname = 'wet deposition flux in-cloud convective warm impaction'
      CASE (18)
        ptrwdep => wdep_incl_cmn
        cdiagname = 'wdep_incl_cmn'
        clongname = 'wet deposition flux in-cloud convective mixed nucleation'
      CASE (19)
        ptrwdep => wdep_incl_cmi
        cdiagname = 'wdep_incl_cmi'
        clongname = 'wet deposition flux in-cloud convective mixed impaction'
      CASE (20)
        ptrwdep => wdep_incl_ccn
        cdiagname = 'wdep_incl_ccn'
        clongname = 'wet deposition flux in-cloud convective cold nucleation'
      CASE (21)
        ptrwdep => wdep_incl_cci
        cdiagname = 'wdep_incl_cci'
        clongname = 'wet deposition flux in-cloud convective cold impaction'
      END SELECT
      lpost = st1_in_st2_proof( cdiagname, wetdepnam) .AND. wetdep_lpost
      CALL new_diag_list (ptrwdep, swetdep, diagname=cdiagname, tsubmname='',    &
                          longname=TRIM(clongname),&
                          units='kg m-2 s-1', ndims=2,                    &
                          nmaxkey=(/ntrac, nspec, nclass, nclass, 0 /), lpost=lpost )
      ! add diagnostic elements only when output is activated
      IF (lpost) THEN
        IF (ANY(tracflag)) CALL new_diag(ptrwdep, ntrac, tracflag, tracname, BYTRACER)
        IF (ANY(specflag)) CALL new_diag(ptrwdep, nspec, specflag, specname, BYSPECIES)
        IF (ANY(modflag))  THEN
            CALL new_diag(ptrwdep, nclass , modflag,  modname,    BYMODE)    
            CALL new_diag(ptrwdep, nclass , modflag,  modnumname, BYNUMMODE) !SF #299: add mode number
        END IF
      END IF
    END DO

  END SUBROUTINE init_wetdep_stream
  !! ---------------------------------------------------------------------------------------
  !!SF prep_wetdep: preparation routine to wet deposition calculations 
  !!                (collect hydrological variables)

  SUBROUTINE prep_wetdep_hydro(kproma, kbdim,    klev, ktop, krow,      lstrat, &
                               pdpg,   pmratepr, pmrateps,   pmsnowacl,         &
                               pmlwc,  pmiwc,    paclc,                         &
                               pice,   peffwat,  peffice,    pmfu,      prhou,  &
                               pfrain, pfsnow,   pfevapr,    pfsubls,   prevap, &
                               pclceff  )

  USE mo_time_control,     ONLY: time_step_len

  INTEGER,  INTENT(in)    :: kproma                     ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                      ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: klev                       ! numer of levels
  INTEGER,  INTENT(in)    :: ktop                       ! top layer index
  INTEGER,  INTENT(in)    :: krow                       ! geographic block number
  LOGICAL,  INTENT(in)    :: lstrat                     ! indicates stratiform clouds (call from cloud)
  REAL(dp), INTENT(in)    :: pdpg(kbdim,klev)           ! dp/g
  REAL(dp), INTENT(in)    :: pmratepr(kbdim,klev)       ! rain formation in cloudy part
  REAL(dp), INTENT(in)    :: pmrateps(kbdim,klev)       ! snow formation in cloudy part
  REAL(dp), INTENT(in)    :: pmsnowacl(kbdim,klev)      ! accretion rate of snow with cloud drop. in cl [kg/kg]
  REAL(dp), INTENT(in)    :: pmfu(kbdim,klev)           ! convective flux
  REAL(dp), INTENT(in)    :: prhou(kbdim,klev)          ! air density
  REAL(dp), INTENT(in)    :: pfrain(kbdim,klev)         ! rain rate
  REAL(dp), INTENT(in)    :: pfsnow(kbdim,klev)         ! snow rate
  REAL(dp), INTENT(in)    :: pfevapr(kbdim,klev)        ! evaporation of rain
  REAL(dp), INTENT(in)    :: pfsubls(kbdim,klev)        ! sublimation of snow
  REAL(dp), INTENT(inout) :: pmlwc(kbdim,klev)          ! liquid water content before rain
  REAL(dp), INTENT(inout) :: pmiwc(kbdim,klev)          ! ice water content before snow
  REAL(dp), INTENT(inout) :: paclc(kbdim,klev)          ! cloud cover
  REAL(dp), INTENT(out)   :: pice(kbdim,klev)           ! ice mass fraction over total water
  REAL(dp), INTENT(out)   :: peffwat(kbdim,klev)        ! autoconversion rate (liq water)
  REAL(dp), INTENT(out)   :: peffice(kbdim,klev)        ! autoconversion rate (ice)
  REAL(dp), INTENT(out)   :: prevap(kbdim,klev)         !
  REAL(dp), INTENT(out)   :: pclceff(kbdim,klev)        !

  !--- Local arrays:
  LOGICAL  :: lo_2d(1:kbdim,klev)
  REAL(dp) :: zilwc(kbdim,klev), zfprec(kbdim,klev)

  !--- Constants:
  REAL(dp), PARAMETER :: zmin  = 1.e-10_dp, &
                         zwu   = 2.0_dp        ! Assumed updraft velocity in convective clouds [m s-1]:

  !--- 1.1) Calculate ice mass-fraction of the total water:

  pmiwc(1:kproma,:) = MAX(pmiwc(1:kproma,:),0._dp)
  pmlwc(1:kproma,:) = MAX(pmlwc(1:kproma,:),0._dp)

  zilwc(1:kproma,:) = pmiwc(1:kproma,:)+pmlwc(1:kproma,:)

  lo_2d(1:kproma,:) = (zilwc(1:kproma,:) > zmin)

  paclc(1:kproma,:) = MERGE(paclc(1:kproma,:)                  , 0._dp, lo_2d(1:kproma,:))
  pice(1:kproma,:)  = MERGE(pmiwc(1:kproma,:)/zilwc(1:kproma,:), 0._dp, lo_2d(1:kproma,:))

  !--- 1.2) Calculate autoconversion rate:

  lo_2d(1:kproma,:) = (pmiwc(1:kproma,:) > zmin)

  peffice(1:kproma,:) = MERGE(pmrateps(1:kproma,:)/pmiwc(1:kproma,:), 0._dp, lo_2d(1:kproma,:))
  peffice(1:kproma,:) = MAX(0._dp,MIN(1._dp,peffice(1:kproma,:)))

  lo_2d(1:kproma,:) = (pmlwc(1:kproma,:) > zmin)

  peffwat(1:kproma,:) = MERGE( (pmratepr(1:kproma,:)+pmsnowacl(1:kproma,:))/pmlwc(1:kproma,:) , &
                               0._dp, lo_2d(1:kproma,:))
  peffwat(1:kproma,:) = MAX(0._dp,MIN(1._dp,peffwat(1:kproma,:)))

  !--- 1.3) Calculate the effective grid-box fraction
  !         affected by precipitation (zclceff):

  !--- Precipitation formation (store for diagnostics [kg m-2 s-1]):
  precipform(1:kproma,:,krow) = (pmratepr(1:kproma,:)+pmrateps(1:kproma,:)+pmsnowacl(1:kproma,:)) &
                                   * pdpg(1:kproma,:) / time_step_len

  !--- 1.3.1) Stratiform clouds:
  !--- 1.3.2) Convective clouds:

  IF (.not. lstrat) THEN

    !--- Assume updraft area fraction as precipitiating
    !    fraction of grid box
    !
    !    Estimate updraft area from prescribed updraft
    !    velocities

    pclceff(1:kproma,:)       = pmfu(1:kproma,:) / (zwu*prhou(1:kproma,:))
    uparfrac(1:kproma,:,krow) = pclceff(1:kproma,:)

  END IF

  !--- 1.4) Re-evaporation

  zfprec(1:kproma,:) = pfrain(1:kproma,:) + pfsnow(1:kproma,:)
                               
  lo_2d(1:kproma,:) = (zfprec(1:kproma,:) > zmin)

  prevap(1:kproma,:) = MERGE((pfevapr(1:kproma,:)+pfsubls(1:kproma,:))/zfprec(1:kproma,:), &
                             0._dp, lo_2d(1:kproma,:))

  prevap(1:kproma,:) = MAX(0._dp,MIN(1._dp,prevap(1:kproma,:)))
  
  !--- Store re-evaporation + sublimation for diagnostics:
  precipevap(1:kproma,:,krow) = pfevapr(1:kproma,:)+pfsubls(1:kproma,:)
  
  END SUBROUTINE prep_wetdep_hydro
  
  !! ---------------------------------------------------------------------------------------
  !! SF: get_lfrac: gas tracer liquid fraction calculation

  SUBROUTINE get_lfrac(kproma, kbdim, krow, klev, ktop, &
                     ptemp,  pmlwc, prho, phenry, plfrac)

  USE mo_physical_constants, ONLY: argas,p0sl_bg

  INTEGER,  INTENT(in)    :: kproma                  ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                   ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: krow                    ! for diagnostic purposes
  INTEGER,  INTENT(in)    :: klev                    ! number of levels
  INTEGER,  INTENT(in)    :: ktop                    ! top layer index
  REAL(dp), INTENT(IN)    :: ptemp(kbdim,klev)       ! temperature
  REAL(dp), INTENT(IN)    :: pmlwc(kbdim,klev)       ! liquid water content
  REAL(dp), INTENT(IN)    :: prho(kbdim,klev)        ! density of air  
  REAL(dp), INTENT(IN)    :: phenry(2)               ! Henry's law vars (constant, and activation energy)
  REAL(dp), INTENT(OUT)   :: plfrac(kbdim,klev)      ! liquid tracer fraction
  
  !--- Gas constant in atm M-1 K-1 (Seinfeld & Pandis 2ed. p290)
  REAL(dp), PARAMETER :: zrgas = argas / p0sl_bg * 1.e03_dp,     &
                         zq298 = 1._dp / 298._dp
  
  !--- Local arrays:
  REAL(dp) :: zp(kbdim,klev), zq(kbdim,klev)

  zq(1:kproma,:) = 1._dp / ptemp(1:kproma,:) - zq298
  zp(1:kproma,:) = zrgas*ptemp(1:kproma,:)*pmlwc(1:kproma,:)*prho(1:kproma,:)*1.e-03_dp*phenry(1) &
                 * exp(phenry(2)*zq(1:kproma,:))
  
  plfrac(1:kproma,:) = zp(1:kproma,:) / (1._dp + zp(1:kproma,:))

  ! make sure vector is filled
  IF (kproma.lt.kbdim) plfrac(kproma+1:kbdim,:)=0._dp
  
  END SUBROUTINE get_lfrac
   
  !! ---------------------------------------------------------------------------------------
  !! gas_setscav: rotuine to set the wetdep flags to handle:
  !!  - in-cloud and/or below-cloud scav
  !!  - water and/or ice scav (resp. rain and/or snow)

  SUBROUTINE gas_setscav(kt,                       &
                         kscavICtype, kscavBCtype, &
                         kscavICphase, kscavBCphase)

  ! kscavICtype  = 0 no in-cloud scavenging
  !                1 prescribed in-cloud scavenging params
  !
  ! kscavBCtype  = 0 no below-cloud scavenging
  !                1 prescribed below-cloud scavenging params
  !
  ! kscavICphase = 0 no in-cloud scavenging
  !                1 in-cloud water-only scavenging
  !                2 in-cloud ice-only scavenging
  !                3 in-cloud water+ice scavenging
  !
  ! kscavBCphase = 0 no below-cloud scavenging
  !                1 below-cloud rain-only scavenging
  !                2 below-cloud snow-only scavenging
  !                3 below-cloud rain+snow scavenging

  USE mo_tracdef, ONLY: trlist
  USE mo_species, ONLY: speclist

  INTEGER, INTENT(in)  :: kt
  INTEGER, INTENT(out) :: kscavICtype, kscavBCtype, &
                          kscavICphase, kscavBCphase

  ! local variables
  INTEGER :: ispec

  !--- Normal setting:

  SELECT CASE (trlist%ti(kt)%nwetdep)
     CASE(0) !wet dep off
        kscavICtype  = 0
        kscavBCtype  = 0
        kscavICphase = 0
        kscavBCphase = 0
     CASE(1:) !wet dep on
        kscavICtype  = 1
        kscavBCtype  = 0
        !SF kscavICphase = 1 !SFnote: a value of 1 reproduces the echam5.5-ham2 gas scavenging settings
        kscavICphase = 3  !SFnote: a value of 3 reproduces the echam5-hammoz gas scavenging settings   
        kscavBCphase = 3
  END SELECT

  !--- Tracer- or species-specific setting:

  ispec = trlist%ti(kt)%spid

  ! allow below-cloud scavenging for HNO3:
  IF (speclist(ispec)%shortname == 'HNO3') kscavBCtype = 1

  ! allow below-cloud scavenging for H2SO4
  !SF #508: now allows H2SO4 below-cloud scav *also* in case of pure HAM
  IF (speclist(ispec)%shortname == 'H2SO4') kscavBCtype  = 1

  END SUBROUTINE gas_setscav

  !! ---------------------------------------------------------------------------------------
  !! gas_wetdep: routine to handle gas wet deposition

  SUBROUTINE gas_wetdep(kproma, kbdim,    klev, ktop, &
                        kt,                           &
                        kscavICtype, kscavBCtype,     &
                        kscavICphase, kscavBCphase,   &
                        lstrat,                       &
                        pxtm1,                        &
                        paclc,                        &
                        pmfu,                         &
                        pdpg,                         &
                        peffwat,                      &
                        peffice,                      &
                        prevap,                       &
                        pice,                         &
                        plfrac,                       &
                        pclc,                         &
                        pfrain,                       &
                        pfsnow,                       &
                        pxtte,                        &
                        pxtp10,                       &
                        pxtp1c,                       &
                        pdepint,                      &
                        pdepintbc,                    &
                        pdepintbcr,                   &
                        pdepintbcs,                   &
                        pdepintic,                    &
                        pmfuxt                        )

  USE mo_time_control,  ONLY: time_step_len
  USE mo_tracdef,       ONLY: ntrac

  INTEGER, INTENT(in) :: kproma        ! geographic block number of locations
  INTEGER, INTENT(in) :: kbdim         ! geographic block maximum number of locations
  INTEGER, INTENT(in) :: klev          ! numer of levels
  INTEGER, INTENT(in) :: ktop          ! top layer index
  INTEGER, INTENT(in) :: kt            ! tracer index
  INTEGER, INTENT(in) :: kscavICtype   ! indicates in-cloud scavenging scheme
  INTEGER, INTENT(in) :: kscavBCtype   ! indicates below-cloud scavenging scheme
  INTEGER, INTENT(in) :: kscavICphase  ! indicates in-cloud scavenging by water and/or ice
  INTEGER, INTENT(in) :: kscavBCphase  ! indicates below-cloud scavenging by water and/or ice

  LOGICAL, INTENT(in) :: lstrat        ! flag for stratiform or convective case

  REAL(dp), INTENT(in)   :: pxtm1(kbdim,klev,ntrac), & ! tracer mixing ratio
                            paclc(kbdim,klev),       & ! cloud cover
                            pmfu(kbdim,klev),        & ! convective flux
                            pdpg(kbdim,klev),        & ! grid box thickness
                            peffwat(kbdim,klev),     & !
                            peffice(kbdim,klev),     & !
                            prevap(kbdim,klev),      & ! evaporation rate
                            pice(kbdim,klev),        & ! ice fraction
                            plfrac(kbdim,klev),      & ! liquid fraction of corresponding tracer
                            pclc(kbdim,klev),        & ! fraction of grid covered by precip
                            pfrain(kbdim,klev),      & ! rain flux
                            pfsnow(kbdim,klev)         ! snow flux

  REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ntrac),  & ! tracer tendency
                             pxtp10(kbdim,klev,ntrac), & ! cloud-free mixing ratio
                             pxtp1c(kbdim,klev,ntrac), & ! cloudy mixing ratio
                             pdepint(kbdim),           & ! global scavenged mr
                             pdepintbc(kbdim),         & ! below-cloud scavenged mr
                             pdepintbcr(kbdim),        & ! below-cloud scavenged by rain mr
                             pdepintbcs(kbdim),        & ! below-cloud scavenged by snow mr
                             pdepintic(kbdim),         & ! in-cloud scavenged mr
                             pmfuxt(kbdim,klev,ntrac)    ! updraft mmr

  !--- local variables
  INTEGER  :: jk

  LOGICAL :: ll_cloud_cov(kbdim,klev), ll_prcp(kbdim,klev), ll1(1:kproma,klev) !SF #458

  REAL(dp) :: ztmst
  REAL(dp) :: zdxtevapic(kbdim,klev), zdxtevapbc(kbdim,klev),     &
              zdxtwat(kbdim,klev),  zdxtice(kbdim,klev),          &
              zxtwat(kbdim,klev), zxtice(kbdim,klev),             &
              zxtp10(kbdim,klev), zmf(kbdim,klev),                &
              zdep(kbdim,klev),                                   &
              zxtp1(kbdim,klev), zxtte(kbdim,klev),               &
              zxtfrac_col(kbdim,klev),                            &
              zxtfrac_colr(kbdim,klev), zxtfrac_cols(kbdim,klev), &
              zscavcoefbcr(kbdim,klev), zscavcoefbcs(kbdim,klev), &
              zcoeffr(kbdim,klev), zcoeffs(kbdim,klev),           &
              zcoeff(kbdim,klev),                                 &
              zdxtcol(kbdim,klev),                                &
              zdxtcolr(kbdim,klev), zdxtcols(kbdim,klev)

  REAL(dp) :: ztmp1(kbdim,klev) !SF #458 dummy variable

  !--- Constants:
  REAL(dp), PARAMETER :: zmin =1.e-10_dp

  ztmst = time_step_len

  zdxtevapic(1:kproma,:) = 0._dp
  zdxtevapbc(1:kproma,:) = 0._dp

  zxtfrac_col(1:kproma,:)  = 0._dp
  zxtfrac_colr(1:kproma,:) = 0._dp
  zxtfrac_cols(1:kproma,:) = 0._dp

  zscavcoefbcr(1:kproma,:) = 0._dp
  zscavcoefbcs(1:kproma,:) = 0._dp

  zdxtwat(1:kproma,:) = 0._dp
  zdxtice(1:kproma,:) = 0._dp

  IF (lstrat) THEN !stratiform case

     !--- Weight mixing ratios with cloud fraction:

     pxtp1c(1:kproma,:,kt) = pxtp1c(1:kproma,:,kt)*paclc(1:kproma,:)
     pxtp10(1:kproma,:,kt) = pxtp10(1:kproma,:,kt)*(1._dp-paclc(1:kproma,:))
     zxtp10(1:kproma,:)    = pxtp10(1:kproma,:,kt)
     zmf(1:kproma,:)       = pdpg(1:kproma,:) / ztmst
     !SF note: zxtp10 and zmf are needed in order to transparently
     !         handle the strat and the conv cases without writing too many separate,
     !         but very similar, equations
  ELSE
     zxtp10(1:kproma,:) = 0._dp
     zmf(1:kproma,:)    = pmfu(1:kproma,:)
  ENDIF

  !--- Associate tracer masses in the cloud fraction to water/ice phase
  !    to the respective mass fractions:
  zxtwat(1:kproma,:) = pxtp1c(1:kproma,:,kt)*(1._dp-pice(1:kproma,:))
  zxtice(1:kproma,:) = pxtp1c(1:kproma,:,kt)*pice(1:kproma,:)

  !--- Setting logical masks necessary for MERGE statements (!SF #458)
  ll_cloud_cov(1:kproma,:) = (paclc(1:kproma,:) > zmin)
  ll_prcp(1:kproma,:)      = (pclc(1:kproma,:)  > zmin)

  !--- 1/ Process:

  !--- 1.1/ In-cloud scavenging

  IF(kscavICtype > 0) THEN

     !--- 1.1.1/ Phase-specific calculations:

    IF (IAND(kscavICphase,1) /= 0) THEN !water scavenging on (kscavICphase==1 .or. 3)

        !--- Change in in-cloud (strat) or updraft (conv) tracer concentration:
        !>>SF #458 (replacing where statements)
        zdxtwat(1:kproma,:) = MERGE( &
                                   zxtwat(1:kproma,:)*plfrac(1:kproma,:)*peffwat(1:kproma,:), &
                                   0._dp, &
                                   ll_cloud_cov(1:kproma,:))
        !<<SF #458

        zxtwat(1:kproma,:) = zxtwat(1:kproma,:) - zdxtwat(1:kproma,:)

    ENDIF !water scavenging (kscavICphase)

    IF (IAND(kscavICphase,2) /= 0) THEN !ice scavenging on (kscavICphase==2 .or. 3)

        !--- Change in in-cloud (strat) or updraft (conv) tracer concentration:
        !>>SF #458 (replacing where statements)
        zdxtice(1:kproma,:) = MERGE( &
                                   zxtice(1:kproma,:)*plfrac(1:kproma,:)*peffice(1:kproma,:), &
                                   0._dp, &
                                   ll_cloud_cov(1:kproma,:))
        !<<SF #458 (replacing where statements)
      
        zxtice(1:kproma,:) = zxtice(1:kproma,:) - zdxtice(1:kproma,:)

    ENDIF !ice scavenging (kscavICphase)
    !--- 1.1.2/ Put everything together:

    pxtp1c(1:kproma,:,kt) = zxtwat(1:kproma,:) + zxtice(1:kproma,:)

    !--- Local deposition mass-flux [grid-box mean kg m-2 s-1]:
    zdep(1:kproma,:) = (zdxtwat(1:kproma,:) + zdxtice(1:kproma,:))*zmf(1:kproma,:)

    DO jk=ktop,klev
       !--- Integrated deposition mass flux:
       pdepintic(1:kproma) = pdepintic(1:kproma) + zdep(1:kproma,jk)

       !--- Re-evaporation:
       zdxtevapic(1:kproma,jk) = pdepintic(1:kproma)*prevap(1:kproma,jk)

       !--- Reduce integrated deposition mass flux by re-evap
       pdepintic(1:kproma) = pdepintic(1:kproma) - zdxtevapic(1:kproma,jk)
    ENDDO

     IF (.NOT. lstrat) THEN !conv case

        zxtp10(1:kproma,:) = -zdep(1:kproma,:)/pdpg(1:kproma,:)*ztmst
!SF note: previously, the tendency - instead of the ambient value - was set here.
!         tendency is indeed the only relevent quantity for later,
!         but this has been done to establish more symetry with regard to the stratiform case scavenging,
!         so that things can be generalized more easily (cf my note up where I set zxtp10 in the strat

        !--- Updraft mass flux:

        pmfuxt(1:kproma,:,kt) = pxtp1c(1:kproma,:,kt)*pmfu(1:kproma,:)

     ENDIF !end conv case

  ENDIF !in-cloud scavenging (kscavICtype)

  !--- 1.2/ Below-cloud scavenging
  IF(kscavBCtype > 0) THEN

     !--- 1.2.1/ Phase-specific calculations:

     IF (IAND(kscavBCphase,1) /= 0) THEN !rain scavenging on (kscavBCphase==1 .or. 3)
        !SFtemporary dummy values for now (in waiting to define two separate rain and snow coeffs)
        zscavcoefbcr(1:kproma,:) = 0._dp
     ENDIF

     IF (IAND(kscavBCphase,2) /= 0) THEN !snow scavenging on (kscavBCphase==2 .or. 3)
        !SFtemporary dummy values for now (in waiting to define two separate rain and snow coeffs)
        zscavcoefbcs(1:kproma,:) = 0._dp
     ENDIF
     !--- 1.2.2/ Put everything together:

        !--- Calculate fraction of below cloud scavenged tracer:
        !>>SF #458 (replacing where statements)
        ll1(1:kproma,:) = .NOT. ll_cloud_cov(1:kproma,:) .AND. ll_prcp(1:kproma,:)

        zcoeffr (1:kproma,:) = MERGE( &
                                    -ztmst*zscavcoefbcr(1:kproma,:), &
                                    0._dp, &
                                    ll1(1:kproma,:) )

        zcoeffs (1:kproma,:) = MERGE( &
                                    -ztmst*zscavcoefbcs(1:kproma,:), &
                                    0._dp, &
                                    ll1(1:kproma,:) )

        ztmp1(1:kproma,:) = MERGE(pclc(1:kproma,:), 1._dp, ll1(1:kproma,:)) !SF 1._dp is a dummy value
        !>>SFtemporary (in waiting to define two separate rain and snow coeffs):
        ztmp1(1:kproma,:) = -ztmst * 1.626e-2_dp &
                          * (MAX(0._dp,pfrain(1:kproma,:)+pfsnow(1:kproma,:)) &
                          / ztmp1(1:kproma,:))**0.6169_dp
        !SF The final expression should be: ztmp1(1:kproma,:) = zcoeffr(1:kproma,:)+zcoeffs(1:kproma,:)
        !<<SFtemporary

        zcoeff(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:))

        zxtfrac_colr(1:kproma,:) = 1._dp - EXP(zcoeffr(1:kproma,:))
        zxtfrac_cols(1:kproma,:) = 1._dp - EXP(zcoeffs(1:kproma,:))
        zxtfrac_col(1:kproma,:)  = 1._dp - EXP(zcoeff(1:kproma,:))
        
        zxtfrac_colr(1:kproma,:) = MAX(0._dp, MIN(1._dp, zxtfrac_colr(1:kproma,:) ) )
        zxtfrac_cols(1:kproma,:) = MAX(0._dp, MIN(1._dp, zxtfrac_cols(1:kproma,:) ) )
        zxtfrac_col(1:kproma,:)  = MAX(0._dp, MIN(1._dp, zxtfrac_col(1:kproma,:) ) )

        ztmp1(1:kproma,:) = pxtp10(1:kproma,:,kt)*pclc(1:kproma,:)*zxtfrac_col(1:kproma,:)
        zdxtcol(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll_prcp(1:kproma,:))

        zxtp10(1:kproma,:) = zxtp10(1:kproma,:) - zdxtcol(1:kproma,:)

        zdxtcol(1:kproma,:) = zdxtcol(1:kproma,:)*pdpg(1:kproma,:)/ztmst

        ztmp1(1:kproma,:) = pxtp10(1:kproma,:,kt) * pclc(1:kproma,:) &
                          * zxtfrac_colr(1:kproma,:) * pdpg(1:kproma,:) / ztmst

        zdxtcolr(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll_prcp(1:kproma,:))

        ztmp1(1:kproma,:) = pxtp10(1:kproma,:,kt) * pclc(1:kproma,:) &
                          * zxtfrac_cols(1:kproma,:) * pdpg(1:kproma,:) / ztmst

        zdxtcols(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll_prcp(1:kproma,:))
        !<<SF #458 (replacing where statements)

        DO jk=ktop,klev

           pdepintbc(1:kproma)  = pdepintbc(1:kproma)  + zdxtcol(1:kproma,jk) *pdpg(1:kproma,jk)/ztmst
           pdepintbcr(1:kproma) = pdepintbcr(1:kproma) + zdxtcolr(1:kproma,jk)*pdpg(1:kproma,jk)/ztmst
           pdepintbcs(1:kproma) = pdepintbcs(1:kproma) + zdxtcols(1:kproma,jk)*pdpg(1:kproma,jk)/ztmst

!temporary, don't take re-evaporation into account
!SF note: inherited from former code. Is it still relevant to dismiss re-evaporation ?

!!$           !--- Re-evaporation:
!!$           zdxtevapbc(1:kproma,jk)  = pdepintbc(1:kproma) *prevap(1:kproma,jk)
!!$           zdxtevapbcr(1:kproma,jk) = pdepintbcr(1:kproma)*prevap(1:kproma,jk)
!!$           zdxtevapbcs(1:kproma,jk) = pdepintbcs(1:kproma)*prevap(1:kproma,jk)
!!$
!!$           !--- Reduce integrated deposition mass flux by re-evap
!!$           pdepintbc(1:kproma)  = pdepintbc(1:kproma)  - zdxtevapbc(1:kproma,jk)
!!$           pdepintbcr(1:kproma) = pdepintbcr(1:kproma) - zdxtevapbcr(1:kproma,jk)
!!$           pdepintbcs(1:kproma) = pdepintbcs(1:kproma) - zdxtevapbcs(1:kproma,jk)
!end temporary
        ENDDO !jk

  ENDIF !end below-cloud scavenging (kscavBCtype)

  !--- put together all contributions to integrated deposition mass flux
  pdepint(1:kproma) = pdepintic(1:kproma) + pdepintbc(1:kproma)

  !--- Total tendency + update updraft flux in conv case:
  IF (lstrat) THEN !SF stratiform case

     pxtp10(1:kproma,:,kt) = zxtp10(1:kproma,:)                              &
                           + (zdxtevapic(1:kproma,:)+zdxtevapbc(1:kproma,:)) &
                             / pdpg(1:kproma,:) * ztmst

     zxtp1(1:kproma,:) = pxtm1(1:kproma,:,kt)+pxtte(1:kproma,:,kt)*ztmst
     zxtte(1:kproma,:) = (pxtp10(1:kproma,:,kt)+pxtp1c(1:kproma,:,kt)-zxtp1(1:kproma,:)) / ztmst

  ELSE !SF conv case

     zxtte(1:kproma,:) = zxtp10(1:kproma,:) / ztmst                                           &
                       + (zdxtevapic(1:kproma,:) + zdxtevapbc(1:kproma,:)) / pdpg(1:kproma,:)

  ENDIF

  pxtte(1:kproma,:,kt)=pxtte(1:kproma,:,kt)+zxtte(1:kproma,:)

  END SUBROUTINE gas_wetdep

  !! ---------------------------------------------------------------------------------------
  !! wetdep_interface: generic interface routine to wet deposition

  SUBROUTINE wetdep_interface(kproma, kbdim, klev, ktop, krow, lstrat, &
                              pdpg, pmratepr, pmrateps, pmsnowacl,     &
                              pmlwc, pmiwc,                            &
                              ptm1, pxtm1, plfrac_so2,                 &
                              pxtte, pxtp10, pxtp1c,                   &
                              pfrain, pfsnow, pfevapr, pfsubls,        &
                              pmfu, pmfuxt,                            &
                              paclc, pclc, prhou, pxtbound)

  USE mo_time_control,  ONLY: delta_time
  USE mo_submodel,      ONLY: lham
  USE mo_tracdef,       ONLY: ln, ntrac, trlist, AEROSOLMASS, AEROSOLNUMBER, GAS
  USE mo_species,       ONLY: speclist
  USE mo_ham_wetdep,    ONLY: ham_wetdep, ham_setscav, prep_ham_mode_init
  USE mo_submodel_diag, ONLY: get_diag_pointer

  INTEGER,  INTENT(in)    :: kproma                      ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                       ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: klev                        ! numer of levels
  INTEGER,  INTENT(in)    :: ktop                        ! top layer index
  INTEGER,  INTENT(in)    :: krow                        ! geographic block number
  LOGICAL,  INTENT(in)    :: lstrat                      ! indicates stratiform clouds (call from cloud)
  REAL(dp), INTENT(in)    :: ptm1     (kbdim,klev)       ! air temperature (t-dt)
  REAL(dp), INTENT(in)    :: pclc     (kbdim,klev)       ! fraction of grid box covered by precip
  REAL(dp), INTENT(in)    :: pfrain   (kbdim,klev)       ! rain flux before evaporation [kg/m2/s]
  REAL(dp), INTENT(in)    :: pfsnow   (kbdim,klev)       ! snow flux before evaporation [kg/m2/s]
  REAL(dp), INTENT(in)    :: pfevapr  (kbdim,klev)       ! evaporation of rain [kg/m2/s]
  REAL(dp), INTENT(in)    :: pfsubls  (kbdim,klev)       ! sublimation of snow [kg/m2/s]
  REAL(dp), INTENT(in)    :: pmsnowacl(kbdim,klev)       ! accretion rate of snow with cloud droplets 
                                                         ! in cloudy part [kg/kg]
  REAL(dp), INTENT(in)    :: pdpg     (kbdim,klev)       ! dp/g
  REAL(dp), INTENT(in)    :: pmfu     (kbdim,klev)       ! see cuflx
  REAL(dp), INTENT(in)    :: prhou    (kbdim,klev)       ! air density
  REAL(dp), INTENT(in)    :: paclc    (kbdim,klev)       ! cloud fraction
  REAL(dp), INTENT(in)    :: pmratepr (kbdim,klev)       ! rain formation rate in cloudy part
  REAL(dp), INTENT(in)    :: pmrateps (kbdim,klev)       ! ice  formation rate in cloudy part
  REAL(dp), INTENT(in)    :: pxtm1    (kbdim,klev,ntrac) ! tracer mass/number mixing ratio (t-dt)
  REAL(dp), INTENT(in)    :: plfrac_so2 (kbdim,klev)     ! liquid fraction of SO2 (for HAM)
  REAL(dp), INTENT(inout) :: pmlwc    (kbdim,klev)       ! cloud liquid water
  REAL(dp), INTENT(inout) :: pmiwc    (kbdim,klev)       ! cloud ice
  REAL(dp), INTENT(inout) :: pxtte    (kbdim,klev,ntrac) ! tracer mass/number mixing ratio tendency
  REAL(dp), INTENT(inout) :: pxtp1c   (kbdim,klev,ntrac) ! ambient  tracer mass mixing ratio (t+dt)
  REAL(dp), INTENT(inout) :: pxtp10   (kbdim,klev,ntrac) ! in-cloud tracer mass mixing ratio (t+dt)
  REAL(dp), INTENT(inout) :: pmfuxt   (kbdim,klev,ntrac) ! updraft mass flux
  REAL(dp), INTENT(inout) :: pxtbound (kbdim,ntrac)      ! conv massfix boundary condition

  !--- local variables
  INTEGER  :: jt, ispec, ierr,          &
              jscavICtype, jscavBCtype, &
              jscavICphase, jscavBCphase

  CHARACTER(len=ln) :: basename, modulename

  REAL(dp) :: zaclc(kbdim,klev), zrevap(kbdim,klev),        &
              zice(kbdim,klev),  zeffwat(kbdim,klev),       &
              zeffice(kbdim,klev),                          &
              zdepint(kbdim), zdepintbc(kbdim),             &
              zdepintbcr(kbdim), zdepintbcs(kbdim),         &
              zdepintic(kbdim), zdepintic_nucw(kbdim),      &
              zdepintic_nucm(kbdim), zdepintic_nucc(kbdim), &
              zdepintic_impw(kbdim), zdepintic_impm(kbdim), &
              zdepintic_impc(kbdim),                        &
              zhenry(2), zlfrac(kbdim,klev),                &
              zclceff(kbdim,klev), zclc(kbdim,klev),        &
              zxtp1c_sav(kbdim,klev,ntrac) !SF this variable is somehow a duplicate of pxtp1c, but untouched
                                           !   upon current scavenging, in contrast to pxtp1c which is 
                                           !   an INOUT of ham_wetdep.
                                           !   zxtp1c_sav is necessary for nucleation ice scavenging
                                           !   because this calculation presents some mode inter-dependencies.
                                           !   Taking a pristine pxtp1c allows then to avoid any risk of
                                           !   getting spurious tracer-loop order dependence in the final
                                           !   result.

  REAL(dp), POINTER    :: fld2d(:,:)         ! pointer for diagnostics

  !--- Prepare wet deposition calculations:

  !----- Hydrological variables:
  zaclc(1:kproma,:) = paclc(1:kproma,:)  !SF zaclc is modified from paclc in prep_wetdep_hydro
                                         ! for wetdep calcs
                                         ! that's why I don't directly use paclc here

  CALL prep_wetdep_hydro(kproma, kbdim,    klev, ktop, krow, lstrat, &
                         pdpg,   pmratepr, pmrateps,   pmsnowacl,    &
                         pmlwc,  pmiwc,    zaclc,                    &
                         zice,   zeffwat,  zeffice, pmfu,    prhou,  &
                         pfrain, pfsnow,   pfevapr, pfsubls, zrevap, &
                         zclceff)

  !----- Mode-wise initializations (HAM-specific, but not M7-specific):
  IF (lham) CALL prep_ham_mode_init(kproma, kbdim, klev)

  !--- Back up the current in-cloud mixing ratios (multiplied by cloud cover)
  DO jt=1,ntrac
     zxtp1c_sav(1:kproma,:,jt) = pxtp1c(1:kproma,:,jt) * paclc(1:kproma,:)
  ENDDO
!SFNote the above to refactor. probably should use zaclc also....

  !--- tracer loop
  DO jt=1, ntrac

    !--- initialize diagnostics-related variables
    zdepint(1:kproma)        = 0._dp
    zdepintbc(1:kproma)      = 0._dp
    zdepintbcr(1:kproma)     = 0._dp
    zdepintbcs(1:kproma)     = 0._dp
    zdepintic(1:kproma)      = 0._dp
    zdepintic_nucw(1:kproma) = 0._dp
    zdepintic_nucm(1:kproma) = 0._dp
    zdepintic_nucc(1:kproma) = 0._dp
    zdepintic_impw(1:kproma) = 0._dp
    zdepintic_impm(1:kproma) = 0._dp
    zdepintic_impc(1:kproma) = 0._dp

    IF(lstrat) THEN
       zclc(1:kproma,:) = pclc(1:kproma,:)
    ELSE
       zclc(1:kproma,:) = zclceff(1:kproma,:)
    ENDIF

    IF (trlist%ti(jt)%nwetdep == 0) CYCLE    ! do nothing if tracer doesn't get wet deposited

       IF (trlist%ti(jt)%nphase == AEROSOLMASS .OR. trlist%ti(jt)%nphase == AEROSOLNUMBER) THEN
       !--- aerosol wet deposition

        IF (lham) THEN
            !--- set scavenging type and phase flags according to the nwetdep scheme:
            !    note: this can be used as a hook to set some tracer- or species-specific behaviours
            CALL ham_setscav(jt,                       &
                             jscavICtype, jscavBCtype, &
                             jscavICphase, jscavBCphase)

            !SF: avoid unnecessary calculations and potentially spurious effects
            !    in case a peculiar scavenging setup has been defined via the previous 
            !    call to ham_setscav
            IF ( (jscavICtype+jscavBCtype) == 0 ) CYCLE

            !--- process scavenging:
            CALL ham_wetdep(kproma, kbdim, klev, krow, ktop,                & 
                            jt,                                             &
                            jscavICtype, jscavBCtype,                       &
                            jscavICphase, jscavBCphase,                     &
                            lstrat,                                         &
                            ptm1, pxtm1, pxtte, pxtp10, pxtp1c, zxtp1c_sav, &
                            pfrain, pfsnow, zaclc, pmfu, pmfuxt,            &
                            prhou,                                          &
                            pdpg, zice, zeffice, zeffwat, zclc,             &
                            zrevap, zdepint, zdepintbc, zdepintbcr,         &
                            zdepintbcs, zdepintic, zdepintic_nucw,          &
                            zdepintic_nucm, zdepintic_nucc,                 &
                            zdepintic_impw, zdepintic_impm,                 &
                            zdepintic_impc)                          

        ENDIF

     ELSEIF (trlist%ti(jt)%nphase == GAS) THEN
     !--- gas wet deposition

        !--- set scavenging type and phase flags according to the nwetdep scheme:
        !    note: this can be used as a hook to set some tracer- or species-specific behaviours 
        CALL gas_setscav(jt,                       &
                         jscavICtype, jscavBCtype, &
                         jscavICphase, jscavBCphase)

       !SF: avoid unnecessary calculations and potentially spurious effects
       !    in case a peculiar scavenging setup has been defined via the previous 
       !    call to ham_setscav
       IF ( (jscavICtype+jscavBCtype) == 0 ) CYCLE

        !--- computes the liquid fraction:
        ispec = trlist%ti(jt)%spid
        zhenry(:) = speclist(ispec)%henry(:)    ! Henry's law constant and activation energy

        CALL get_lfrac(kproma, kbdim, krow, klev, ktop, &
                       ptm1, pmlwc, prhou, zhenry, zlfrac)

        !--- overrides the liquid fraction value for HAM SO2 and HAM SO4_gas
        IF (lham) THEN

           basename   = trlist%ti(jt)%basename
           modulename = trlist%ti(jt)%modulename

           IF (modulename == "HAM") THEN
              IF (basename == "SO2") THEN
                 zlfrac(1:kproma,:) = plfrac_so2(1:kproma,:)
              ELSEIF (basename == "H2SO4") THEN
                 zlfrac(1:kproma,:) = 1._dp
              ENDIF
           ENDIF
        ENDIF

        !--- process scavenging:
        CALL gas_wetdep(kproma, kbdim, klev, ktop,                  &
                        jt,                                         &
                        jscavICtype, jscavBCtype,                   &
                        jscavICphase, jscavBCphase,                 &
                        lstrat,                                     &
                        pxtm1, zaclc, pmfu, pdpg, zeffwat, zeffice, &
                        zrevap, zice, zlfrac, zclc, pfrain, pfsnow, &
                        pxtte, pxtp10, pxtp1c, zdepint, zdepintbc,  &
                        zdepintbcr, zdepintbcs, zdepintic, pmfuxt   )

     ENDIF !end aerosol or gas scavenging

     !--- update the boundary condition (instantaneous wet deposition) for xt_conv_massfix
     pxtbound(1:kproma,jt) = zdepint(1:kproma)

     !--- diagnostics
     CALL get_diag_pointer(wdep, fld2d, jt, ierr=ierr)
     IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow)  &
                                            + zdepint(1:kproma)*delta_time

     CALL get_diag_pointer(wdep_incl, fld2d, jt, ierr=ierr)
     IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                           + zdepintic(1:kproma)*delta_time

     CALL get_diag_pointer(wdep_blcl, fld2d, jt, ierr=ierr)
     IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                           + zdepintbc(1:kproma)*delta_time

     IF (lstrat) THEN

        !stratiform clouds
        CALL get_diag_pointer(wdep_strat, fld2d, jt, ierr=ierr)
        IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                            + zdepint(1:kproma)*delta_time

        IF (lwetdepdetail) THEN !detailed diag

           CALL get_diag_pointer(wdep_incl_swn, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_nucw(1:kproma)*delta_time

           CALL get_diag_pointer(wdep_incl_smn, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_nucm(1:kproma)*delta_time

           CALL get_diag_pointer(wdep_incl_scn, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_nucc(1:kproma)*delta_time

           CALL get_diag_pointer(wdep_incl_swi, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_impw(1:kproma)*delta_time

           CALL get_diag_pointer(wdep_incl_smi, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_impm(1:kproma)*delta_time

           CALL get_diag_pointer(wdep_incl_sci, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_impc(1:kproma)*delta_time

           CALL get_diag_pointer(wdep_blcl_sr, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintbcr(1:kproma)*delta_time
        
           CALL get_diag_pointer(wdep_blcl_ss, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintbcs(1:kproma)*delta_time
        
        ENDIF

     ELSE

        !convective clouds             
        CALL get_diag_pointer(wdep_conv, fld2d, jt, ierr=ierr)
        IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) & 
                                            + zdepint(1:kproma)*delta_time
        
        IF (lwetdepdetail) THEN !detailed diag
        
           CALL get_diag_pointer(wdep_incl_cwn, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_nucw(1:kproma)*delta_time
   
           CALL get_diag_pointer(wdep_incl_cmn, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_nucm(1:kproma)*delta_time
   
           CALL get_diag_pointer(wdep_incl_ccn, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_nucc(1:kproma)*delta_time
           
           CALL get_diag_pointer(wdep_incl_cwi, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_impw(1:kproma)*delta_time

           CALL get_diag_pointer(wdep_incl_cmi, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_impm(1:kproma)*delta_time

           CALL get_diag_pointer(wdep_incl_cci, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintic_impc(1:kproma)*delta_time

           CALL get_diag_pointer(wdep_blcl_cr, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintbcr(1:kproma)*delta_time
                                               
           CALL get_diag_pointer(wdep_blcl_cs, fld2d, jt, ierr=ierr)
           IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow) &
                                               + zdepintbcs(1:kproma)*delta_time
              
         ENDIF
     ENDIF                                     

  END DO !end loop over tracers

  END SUBROUTINE wetdep_interface

END MODULE mo_hammoz_wetdep
