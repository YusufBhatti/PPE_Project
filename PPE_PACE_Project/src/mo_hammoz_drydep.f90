!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_drydep.f90
!!
!! \brief
!! Module to interface ECHAM submodules with dry deposition module(s).
!!
!! \author M. Schultz (FZ Juelich)
!! \author Grazia Frontoso (C2SM)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2009-10-26)
!!   -# M. Schultz (FZ Juelich) - improved diag routines (2010-04-16)
!!   -# Grazia Frontoso (C2SM)  - usage of the input variables defined over land, water, ice to account
!!                                for the non-linearity in the drydep calculations for gridboxes
!!                                containing both water and sea ice (2012-02-01) 
!!
!! \limitations
!! All diag_lists must be defined in order to avoid problems with
!! get_diag_pointer in the actual sedi_interface routine. Lists can be empty.
!!
!! \details
!! Currently there is only one unified interactive drydep scheme for
!! aerosols (HAM) and gas-phase species (MOZ).
!! This module initializes the scheme based on the namelist parameters
!! in submodeldiagctl and creates a stream for variable pointers and 
!! diagnostic quantities used in the dry deposition scheme. It also
!! provides a generic interface to the actual dry deposition routine(s).
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
MODULE mo_hammoz_drydep

  USE mo_kind,             ONLY: dp
  USE mo_submodel_diag,    ONLY: t_diag_list

  IMPLICIT NONE

  PRIVATE

  ! public variables  (see declaration below)

  ! subprograms
  PUBLIC                       :: init_drydep_stream, &
                                  drydep_interface, drydep_init

  ! drydep_stream
  INTEGER, PARAMETER           :: ndrydepvars=2
  CHARACTER(LEN=32)            :: drydepvars(1:ndrydepvars)= &
                                (/'ddep             ', &   ! total drydep flux
                                  'vddep            ' /)   ! deposition velocity

  ! variable pointers and diagnostic lists
!!  TYPE (t_diag_list), PUBLIC   :: ddepinst     ! inst. dry deposition flux     ### needed ??? ###
  TYPE (t_diag_list), PUBLIC   :: ddep         ! dry deposition flux
  TYPE (t_diag_list), PUBLIC   :: vddep        ! dry deposition velocity
 
  INTEGER :: idt_ddep_detail

  CONTAINS

!>>gf #244
  SUBROUTINE drydep_init

    USE mo_hammoz_drydep_lg,    ONLY: drydep_lg_init
    USE mo_tracdef,             ONLY: ntrac, trlist

    !--- initialize the interactive drydep scheme (ndrydep==2)
    IF (ANY(trlist%ti(1:ntrac)%ndrydep == 2)) CALL drydep_lg_init

  END SUBROUTINE drydep_init
!<<gf

  SUBROUTINE init_drydep_stream

    USE mo_string_utls,         ONLY: st1_in_st2_proof
    USE mo_util_string,         ONLY: tolower
    USE mo_exception,           ONLY: finish
    USE mo_memory_base,         ONLY: t_stream, new_stream, &
                                      default_stream_setting, &
                                      add_stream_reference, &
                                      AUTO, SURFACE
    USE mo_ham_m7_trac,         ONLY: ham_get_class_flag
    USE mo_tracer,              ONLY: validate_traclist
    USE mo_tracdef,             ONLY: ln, ntrac, trlist, GAS, AEROSOL
    USE mo_species,             ONLY: nspec, speclist 
    USE mo_ham,                 ONLY: nclass
    USE mo_submodel_streams,    ONLY: drydep_lpost, drydep_tinterval, drydepnam,  &
                                      drydep_gastrac, drydep_keytype, drydep_ldetail, &    ! ++mgs 20140519
                                      drydep_trac_detail
    USE mo_submodel_diag,       ONLY: new_diag_list, new_diag,  &
                                      BYTRACER, BYSPECIES, BYNUMMODE, BYMODE !SF #299 added BYMODE
    USE mo_hammoz_drydep_lg,    ONLY: init_drydep_lg_stream
    USE mo_submodel,            ONLY: lham !SF, see #228

    ! local variables
    INTEGER, PARAMETER             :: ndefault = 2
    CHARACTER(LEN=32)              :: defnam(1:ndefault)   = &   ! default diagnostics
                                (/ 'ddep             ', &        ! total drydep flux
                                   'vddep            ' /)        ! dry deposition velocity
    CHARACTER(len=ln)              :: defaultgas(5)        = &   ! default gas-phase tracers for diagnostics
                                (/ 'SO2     ',               &
                                   'H2SO4   ',               &
                                   'HNO3    ',               &
                                   'O3      ',               &
                                   'NO2     '           /)
    LOGICAL                        :: tracflag(ntrac), specflag(nspec), modflag(MAX(nclass,1))
    CHARACTER(LEN=ln)              :: tracname(ntrac), specname(nspec), modname(MAX(nclass,1)), &
                                      modnumname(MAX(nclass,1)) !SF #299
    TYPE (t_stream), POINTER       :: sdrydep
    INTEGER                        :: ierr, jt
    LOGICAL                        :: lpost

    !++mgs: default values and namelist read are done in init_submodel_streams !

    !-- handle ALL, DETAIL and DEFAULT options for drydep output variables
    !-- Note: ALL and DETAIL are identical for drydep output
    IF (TRIM(tolower(drydepnam(1))) == 'detail')  drydepnam(1:ndrydepvars) = drydepvars(:)
    IF (TRIM(tolower(drydepnam(1))) == 'all')     drydepnam(1:ndrydepvars) = drydepvars(:)
    IF (TRIM(tolower(drydepnam(1))) == 'default') drydepnam(1:ndefault) = defnam(:)

    !-- check that all diagnostic names from namelist are valid
    IF (.NOT. st1_in_st2_proof( drydepnam, drydepvars, ierr=ierr) ) THEN
      IF (ierr > 0) CALL finish ( 'ini_drydep_stream', 'variable '// &
                                  drydepnam(ierr)//' does not exist in drydep stream' )
    END IF

    !-- find out which gas-phase tracers shall be included in diagnostics
    CALL validate_traclist(drydep_gastrac, defaultgas, nphase=GAS,              &
                           ldrydep=.true.)                   !>>dod SOA: removed ltran <<dod

    !-- define the flags and names for the diagnostic lists. We need one set of flags and
    !   names for each key_type (BYTRACER, BYSPECIES, BYMODE)
    !   gas-phase tracers will always be defined BYTRACER, for aerosol tracers one of the
    !   following lists will be empty.
    !   Note: vddep uses BYTRACER or BYMODE, ddep uses BYTRACER or BYSPECIES
!!++mgs 2015-02-10 : split loop so that gas-phase tracers are always output last
    tracflag(:) = .FALSE.
    DO jt = 1,ntrac
      tracname(jt) = trlist%ti(jt)%fullname
      IF (trlist%ti(jt)%nphase /= GAS) THEN
        IF (drydep_keytype == BYTRACER .AND. nclass > 0) THEN
          tracflag(jt) = trlist%ti(jt)%ndrydep > 0
        END IF
      END IF
    END DO
    specflag(:) = .FALSE.
    DO jt = 1,nspec
      specname(jt) = speclist(jt)%shortname
      IF (drydep_keytype == BYSPECIES .AND.                     &
          IAND(speclist(jt)%nphase, AEROSOL) /= 0 .AND.         &   !>>dod SOA removed check of trtype <<dod
          nclass > 0) THEN
        specflag(jt) = speclist(jt)%ldrydep
      END IF
    END DO
    DO jt = 1,ntrac
      tracname(jt) = trlist%ti(jt)%fullname
      IF (trlist%ti(jt)%nphase == GAS) THEN
        tracflag(jt) = st1_in_st2_proof(trlist%ti(jt)%fullname, drydep_gastrac)
      END IF
    END DO
!--mgs
    modflag(:) = .FALSE.
    modname(:) = ''
    !SF #228, adding a condition to check that HAM is active:
    !SF #299, adding a condition to check if BYMODE is relevant:
    IF (lham .AND. nclass > 0 .AND. (drydep_keytype == BYMODE)) &
       CALL ham_get_class_flag(nclass, modflag, modname, modnumname, ldrydep=.true.)

    !-- open new diagnostic stream
    CALL new_stream (sdrydep,'drydep',lpost=drydep_lpost,lrerun=.FALSE., &
                     interval=drydep_tinterval)
    CALL default_stream_setting (sdrydep, lrerun = .FALSE., &
                     leveltype=SURFACE, &                          !++mgs 20140519: added this for safety
                     contnorest = .TRUE., table = 199, &
                     laccu = .false., code = AUTO)
   
    !-- add standard ECHAM variables
    IF (drydep_lpost) THEN
      CALL add_stream_reference (sdrydep, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
      CALL add_stream_reference (sdrydep, 'lsp'     ,'sp'    ,lpost=.TRUE.)
      CALL add_stream_reference (sdrydep, 'aps'     ,'g3b'   ,lpost=.TRUE.)
      CALL add_stream_reference (sdrydep, 'gboxarea','geoloc',lpost=.TRUE.)
    END IF

    !-- add instantaneous drydep rates (always by tracer)
    !-- these are used in the calculation and must always be present in the stream 
    !-- it is never output 
!!    CALL default_stream_setting (drydep,  lrerun=.TRUE.) 
!!     CALL new_diag_list (ddepinst, drydep, diagname='ddep_inst', tsubmname='',    &
!!                        longname='dry deposition mass flux', units='kg m-2 s-1', &
!!                        ndims=2, nmaxkey=(/ntrac, 0, 0, 0 /) )
!!    DO jt = 1,ntrac
!!      IF (trlist%ti(jt)%ndrydep > 0) THEN
!!        cunit = 'kg m-2 s-1'
!!        CALL new_diag(ddepinst, 'ddepinst_'//TRIM(trlist%ti(jt)%fullname), &
!!                      'dry deposition mass flux of '//TRIM(trlist%ti(jt)%fullname),  &
!!                      cunit, BYTRACER, jt, lpost=.false.)
!!      END IF
!!    END DO

    !-- instantaneous diagnostic quantities

    CALL default_stream_setting (sdrydep,  lrerun=.FALSE., laccu=.FALSE.) 

    !-- drydep velocities
    !-- these are used in the calculation and must always be present in the stream 
    !-- drydep_gastrac only controls output
    lpost = st1_in_st2_proof( 'vddep', drydepnam) .AND. drydep_lpost
    CALL new_diag_list (vddep, sdrydep, diagname='vddep', tsubmname='',    &
                        longname='dry deposition velocity', units='m s-1', &
                        ndims=2, nmaxkey=(/ntrac, 0, nclass, nclass, 0 /), lpost=lpost )
    CALL new_diag(vddep, ntrac, tracflag, tracname, BYTRACER)
    IF (ANY(modflag)) THEN !SF #299 added mode mass diags and fix for mode number name
      CALL new_diag(vddep, nclass, modflag, modname, BYMODE)
      CALL new_diag(vddep, nclass, modflag, modnumname, BYNUMMODE)
    END IF

    !-- averaged diagnostic quantities

    CALL default_stream_setting (sdrydep, lrerun=.FALSE., laccu=.TRUE.)

    !-- total drydep flux
    lpost = st1_in_st2_proof( 'ddep', drydepnam) .AND. drydep_lpost
    CALL new_diag_list (ddep, sdrydep, diagname='ddep', tsubmname='',    &
                        longname='accumulated dry deposition flux',&
                        units='kg m-2 s-1', ndims=2,                    &
                        nmaxkey=(/ntrac, nspec, nclass, nclass, 0 /), lpost=lpost )
    ! add diagnostic elements only when output is activated
    IF (lpost) THEN
      CALL new_diag(ddep, ntrac, tracflag, tracname, BYTRACER)
      CALL new_diag(ddep, nspec, specflag, specname, BYSPECIES)
      IF (ANY(modflag)) THEN !SF #299 added mode mass diags and fix for mode number name
        CALL new_diag(ddep, nclass, modflag, modname, BYMODE)
        CALL new_diag(ddep, nclass, modflag, modnumname, BYNUMMODE)
      END IF
    END IF

 !++mgs 20140519 : detailed diagnostics
    IF (drydep_ldetail) CALL init_drydep_lg_stream(sdrydep, drydep_trac_detail, idt_ddep_detail)

  END SUBROUTINE init_drydep_stream

  !! ---------------------------------------------------------------------------------------
  !! drydep_interface: generic interface routine to dry deposition

  SUBROUTINE drydep_interface(kbdim, kproma,  klev,     krow,           &
!>>gf modified argument list (see#78)
                     pqsfc, pqsssfc, ptsfc, pcfml, pcfmw, pcfmi,        &
                     pcfncl, pcfncw, pcfnci,                            &
                     pepdu2, pkap, pum1, pvm1, pgeom1, pril, priw,      &
                     prii,                                              &
                     ptvir1, ptvl, ptvw, ptvi, paz0,                    &
                     ptslm1, loland,                                    &
                     pfrl,   pfrw,  pfri,     pcvs,   pcvw,     pvgrat, &
                     psrfl,                           pu10,     pv10,   &
                     pxtems, pxtm1, pdensair, paphp1, pforest,  ptsi,   &
                     paz0l, paz0w, paz0i, pcdnl, pcdnw, pcdni           )
!<<gf

  USE mo_exception,        ONLY: finish
  USE mo_time_control,     ONLY: delta_time, time_step_len
  USE mo_physical_constants, ONLY: grav
  USE mo_memory_g3b,       ONLY: wsmx, ws
  USE mo_tracdef,          ONLY: trlist, ntrac
  USE mo_submodel,         ONLY: lham
  USE mo_submodel_diag,    ONLY: get_diag_pointer
  USE mo_hammoz_drydep_lg, ONLY: drydep_lg_calcra, drydep_lg_vdbl
  USE mo_ham_drydep,       ONLY: ham_vdaer, ham_vd_presc
!! USE mo_moz_diag,          ONLY: moz_drydep_diag   ### to be written

  IMPLICIT NONE
  
  
  INTEGER,  INTENT(in)     :: kproma                      ! geographic block number of locations
  INTEGER,  INTENT(in)     :: kbdim                       ! geographic block maximum number of locations 
  INTEGER,  INTENT(in)     :: klev                        ! numer of levels
  INTEGER,  INTENT(in)     :: krow                        ! geographic block number

  REAL(dp), INTENT(in)     :: pqsfc    (kbdim)            ! humidity (lowest level)
  REAL(dp), INTENT(in)     :: pqsssfc  (kbdim)            ! saturation humidity (lowest level)
  REAL(dp), INTENT(in)     :: ptsfc    (kbdim)            ! temperature (lowest level)
!>> gf see #78
  REAL(dp), INTENT(in)     :: pcfml    (kbdim)            ! stability dependend transfer coeff. for momentum over land
  REAL(dp), INTENT(in)     :: pcfmw    (kbdim)            ! stability dependend transfer coeff. for momentum over water
  REAL(dp), INTENT(in)     :: pcfmi    (kbdim)            ! stability dependend transfer coeff. for momentum over ice
  REAL(dp), INTENT(in)     :: pcfncl   (kbdim)            ! function of heat transfer coeff. over land
  REAL(dp), INTENT(in)     :: pcfncw   (kbdim)            ! function of heat transfer coeff. over water
  REAL(dp), INTENT(in)     :: pcfnci   (kbdim)            ! function of heat transfer coeff. over ice
!<<gf
  REAL(dp), INTENT(in)     :: pepdu2                      ! constant
  REAL(dp), INTENT(in)     :: pkap                        ! constant
  REAL(dp), INTENT(in)     :: pum1     (kbdim,klev)       ! u-wind (t-dt)
  REAL(dp), INTENT(in)     :: pvm1     (kbdim,klev)       ! v-wind (t-dt)
  REAL(dp), INTENT(in)     :: pgeom1   (kbdim,klev)       ! geopertential (t-dt)
!>>gf see #78
  REAL(dp), INTENT(in)     :: pril     (kbdim)            ! moist richardson number ocer land
  REAL(dp), INTENT(in)     :: priw     (kbdim)            ! moist richardson number over water
  REAL(dp), INTENT(in)     :: prii     (kbdim)            ! moist richardson number over ice
!<<gf
  REAL(dp), INTENT(in)     :: ptvir1   (kbdim,klev)       ! see vdiff
!>>gf see #78
  REAL(dp), INTENT(in)     :: ptvl     (kbdim)            ! virtual potential temp. over land
  REAL(dp), INTENT(in)     :: ptvw     (kbdim)            ! virtual potential temp. over ocean
  REAL(dp), INTENT(in)     :: ptvi     (kbdim)            ! virtual potential temp. over ice
!<<gf
  REAL(dp), INTENT(in)     :: paz0     (kbdim)            ! roughness length
  REAL(dp), INTENT(in)     :: ptslm1   (kbdim)            ! surface temperature
  LOGICAL,  INTENT(in)     :: loland   (kbdim)            ! land mask
  REAL(dp), INTENT(in)     :: pfrl     (kbdim)            ! land fraction
  REAL(dp), INTENT(in)     :: pfrw     (kbdim)            ! water fraction
  REAL(dp), INTENT(in)     :: pfri     (kbdim)            ! ice fraction
  REAL(dp), INTENT(in)     :: pcvs     (kbdim)            ! snow cover fraction
  REAL(dp), INTENT(in)     :: pcvw     (kbdim)            ! wet skin fraction
  REAL(dp), INTENT(in)     :: pvgrat   (kbdim)            ! vegetation ratio
  REAL(dp), INTENT(in)     :: psrfl    (kbdim)            ! surface solar flux
  REAL(dp), INTENT(in)     :: pu10     (kbdim)            ! 10m w-wind
  REAL(dp), INTENT(in)     :: pv10     (kbdim)            ! 10m v-wind
  REAL(dp), INTENT(inout)  :: pxtems   (kbdim,ntrac)      ! surface emissions
  REAL(dp), INTENT(inout)  :: pxtm1    (kbdim,klev,ntrac) ! tracer mass/number mixing ratio (t-dt)
  REAL(dp), INTENT(in)     :: pdensair (kbdim)            ! air density
  REAL(dp), INTENT(in)     :: paphp1   (kbdim,klev+1)     ! air pressure at layer interface (t+dt)
  REAL(dp), INTENT(in)     :: pforest  (kbdim)            ! forest fraction
  REAL(dp), INTENT(in)     :: ptsi     (kbdim)            ! surface temperature over ice
  REAL(dp), INTENT(in)     :: paz0l    (kbdim)            ! roughness length over land
  REAL(dp), INTENT(in)     :: paz0w    (kbdim)            ! roughness length over water
  REAL(dp), INTENT(in)     :: paz0i    (kbdim)            ! roughness length over ice
!>>gf see #78
  REAL(dp), INTENT(in)     :: pcdnl    (kbdim)            ! see mo_surface_land
  REAL(dp), INTENT(in)     :: pcdnw    (kbdim)            ! see mo_surface_ocean
  REAL(dp), INTENT(in)     :: pcdni    (kbdim)            ! see mo_surface_ice
!<<gf

  !--- Local variables

  INTEGER :: jl, jt, ierr

  REAL(dp), PARAMETER  :: zephum=5.e-2_dp

  REAL(dp):: zdz
  REAL(dp):: zrahwat(kbdim),   zrahice(kbdim),   zrahveg(kbdim),   zrahslsn(kbdim),   &
!>>gf see #78
             zustveg(kbdim),   zustslsn(kbdim),                                       &
             zustarl(kbdim),   zustarw(kbdim),   zustari(kbdim),                      & 
             zvgrat(kbdim),    zcvbs(kbdim),     zrh(kbdim)
!<<gf
  REAL(dp):: zws(kbdim), zwsmx(kbdim)

  REAL(dp):: zvd(kbdim,ntrac), zvdstom(kbdim,ntrac)      ! dry deposition velocity and stomatal ..
  REAL(dp):: zdrydepflux(kbdim)                          ! dry deposition mass flux

  !--- Diagnostic stream:
  REAL(dp), POINTER :: fld2d(:,:) 

  !--- 0) Initialisations: -------------------------------------------------------------

  zvd(1:kproma,:)       = 0._dp
  zvdstom(1:kproma,:)   = 0._dp
  zdrydepflux(1:kproma) = 0._dp

  zws   (1:kproma)      = ws   (1:kproma,krow) !soil wetness
  zwsmx (1:kproma)      = wsmx (1:kproma,krow) !field capacity of soil

  !--- 1) Calculate relative humidity and other parameters needed below
  zrh(1:kproma)=MIN(1._dp,MAX(zephum,pqsfc(1:kproma)/pqsssfc(1:kproma)))

  DO jl=1,kproma
    !--- Calculate bare soil fraction:
    !    It is calculated as residual term over land.
    IF (loland(jl)) THEN
      zvgrat(jl)=pvgrat(jl)                                           ! vegetation fraction
      zcvbs(jl)=(1._dp-pcvs(jl))*(1._dp-pvgrat(jl))*(1._dp-pcvw(jl))  ! bare soil fraction
    ELSE
      zvgrat(jl)=0._dp
      zcvbs(jl)=0._dp
    ENDIF
  ENDDO

  !--- 2) Calculate dry deposition velocities ------------------------------------------
  ! This is done sequentially for 
  !    2.1: prescribed velocities
  !    2.2: gas-phase and aerosol tracers using the Ganzeveld scheme
  ! Each element of trlist has a ndrydep flag which determines the scheme to be used.
  ! Ndrydep==0 means no dry deposition for this tracer.
  ! Each of the routines called in 1.1, 1.2 contain their own tracer loop and
  ! modify zdv only for the tracers they are responsible for.

  IF ( ANY(trlist%ti(:)%ndrydep>2) ) THEN
    CALL finish('drydep_interface', 'Dry deposition with ndrydep > 2 not implemented')
  END IF

  !--- 2.1) Prescibed dry deposition velocities:

  IF( ANY(trlist%ti(:)%ndrydep==1) ) THEN
    CALL ham_vd_presc(kproma, kbdim,  klev,    krow,   loland,   &
                      paphp1, pcvs,   pforest, pfri,   ptsi,     &
                      pcvw,   ptslm1, zws,     zwsmx,  pdensair, &
                      zvd                                         )
  END IF

  !--- 2.2) Explicitly calculated dry deposition velocities (Ganzeveld scheme):

  IF ( ANY(trlist%ti(:)%ndrydep==2) ) THEN
    !--- Calculate the aerodynamic resistance:
    CALL drydep_lg_calcra (kproma,   kbdim,    klev,    krow,                     &
!>>gf modified argument list (see #78)
                           pepdu2, pkap, pum1, pvm1, pgeom1, pril, priw, prii,    &
                           ptvir1, ptvl, ptvw, ptvi, ptslm1, loland,              &
                           pcdnl, pcdnw, pcdni, pcfml, pcfmw, pcfmi,              &
                           pcfncl, pcfncw, pcfnci,                                &
                           paz0w,    paz0i,   paz0l,                              &
                           zrahwat,  zrahice, zrahveg,  zrahslsn,                 &
                           zustarl, zustarw, zustari, zustveg, zustslsn           )
!<<gf

    !--- Calculate the dry deposition velocity for gas-phase species:
    CALL drydep_lg_vdbl (kproma,   kbdim,     klev,    krow,    loland,  &
!>>gf modified argument list (see#78) 
                         psrfl,    ptslm1,    pum1,    pvm1,    zrh,     &
                         pfrl,     pfrw,      pfri,    zcvbs,   pcvs,    &
                         pcvw,     zvgrat,    zrahwat, zrahice, zrahveg, &
                         zrahslsn, zws,       zwsmx,                     &
                         zustarw,  zustari,   zustveg, zustslsn,         &
                         zvd,      zvdstom,   idt_ddep_detail            )
!<<gf

    !--- Calculate the dry deposition velocity for aerosols:
    !### argument ordering should be similar to drydep_lg_vdbl
    IF (lham) THEN
       CALL ham_vdaer  (kproma,   kbdim,    klev,    krow,     loland,  zvgrat,  &
                        pcvs,     pcvw,     pfri,    zcvbs,    pfrw,    pum1,    &  
!>>gf modified argument list (see#78) 
                        pvm1,     zustarl,  zustarw, zustari,  zustveg, zustslsn,&
                        pu10,     pv10,                                          &
                        paz0w,    paz0i,    zrahwat, zrahice,  zrahveg, zrahslsn,&
                        ptslm1,   zrh,      pxtm1,    zvd                        )
!<<gf
    END IF
  END IF

  !--- prevent overflow and store dry deposition velocities for diagnostics
  DO jt=1, ntrac
    DO jl=1, kproma
      ! Security check, limit the deposition velocity to 2x the 
      ! vertical grid velocity of the lowest model layer.
      zdz=(paphp1(jl,klev+1)-paphp1(jl,klev))/(pdensair(jl)*grav)
      zvd(jl,jt) = MIN(zvd(jl,jt) , 2._dp*(zdz/time_step_len))
    END DO

    CALL get_diag_pointer(vddep, fld2d, jt, ierr=ierr)

    IF (ierr == 0) fld2d(1:kproma,krow) = zvd(1:kproma, jt)
!### also store zvdstom (only ozone ??)


  !--- 3) Change emissions tendencies due to dry deposition: --------------
  !--- and store dry depositon mass flux in diagnostics

    IF (trlist%ti(jt)%ndrydep > 0) THEN
      !--- Calculate the tracer flux to the surface that is 
      !    equivalent to the deposition velocity:
      zdrydepflux(1:kproma)=pxtm1(1:kproma,klev,jt)*pdensair(1:kproma)*zvd(1:kproma,jt)

      !--- Reduce emission flux:
      pxtems(1:kproma,jt)=pxtems(1:kproma,jt)-zdrydepflux(1:kproma)
    ELSE ! trlist%ti(jt)%ndrydep <= 0
      zdrydepflux(1:kproma)=0._dp
    END IF
    ! get diagnostics pointer
    CALL get_diag_pointer(ddep, fld2d, jt, ierr=ierr)
    IF (ierr == 0) fld2d(1:kproma,krow)=fld2d(1:kproma,krow)+zdrydepflux(1:kproma)*delta_time
!### add stomatal deposition flux ... ###

    ! special diagnostics for MOZ
    !! CALL moz_drydep_diag()  ....    ### to be done in mo_moz_diag
  END DO

  END SUBROUTINE drydep_interface


END MODULE mo_hammoz_drydep
