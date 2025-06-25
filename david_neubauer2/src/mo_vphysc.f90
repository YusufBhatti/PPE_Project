!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_vphysc
  !-----------------------------------------------------------------------
  ! Module that defines a new stream vphysc containing important variables
  ! of the physc subprogram.
  ! 
  ! contains following subprograms:
  ! 
  ! init_vphysc_stream
  ! set_vphysc_var
  !
  ! is called from init_submodel_streams in mo_submodeldiagctl.f90
  !
  ! Authors:
  ! J.S. Rast, MPI-M, June 2003, original source
  ! J. Kazil, MPI-M, October 2008, added PBL top levepl
  ! L. Kornblueh, MPI-M, November 2008, disentangle
  ! J.S. Rast and M. Schultz, October 2009, extended to disentangle aerosol code
  !                       and added ALL and DEFAULT options
  ! J.S. Rast and M. Schultz, January 2010, changed vphysc variables into structure
  !-----------------------------------------------------------------------

  USE mo_kind,             ONLY: wp

  IMPLICIT NONE

  PRIVATE

  ! variables
  PUBLIC                       :: vphysc, t_vphysc

  ! subprograms
  PUBLIC                       :: init_vphysc_stream, &
                                  set_vphysc_var

  ! type definition
  TYPE t_vphysc
    REAL(wp), POINTER          :: geom1(:,:,:)       ! geopotential height at mid points [m]
    REAL(wp), POINTER          :: geohm1(:,:,:)      ! geopotential height at interfaces [m]
    REAL(wp), POINTER          :: aphm1(:,:,:)       ! air pressure at interfaces
    REAL(wp), POINTER          :: grmassm1(:,:,:)    ! grid box mass at t-dt [kg]
    REAL(wp), POINTER          :: grvolm1(:,:,:)     ! grid box volume at t-dt [m3]
    REAL(wp), POINTER          :: grheightm1(:,:,:)  ! grid box height [m]
    REAL(wp), POINTER          :: rhoam1(:,:,:)      ! dry air mass density [kg m-3]
    REAL(wp), POINTER          :: trpwmo(:,:)        ! tropopause model level index
    REAL(wp), POINTER          :: trpwmop1(:,:)      ! one level below the tropopause
    REAL(wp), POINTER          :: pbl(:,:)           ! PBL height model level
    REAL(wp), POINTER          :: velo10m(:,:)       ! 10-m wind speed [m s-1]
    REAL(wp), POINTER          :: tsw(:,:)           ! sea water temperature
    REAL(wp), POINTER          :: smelt(:,:)         ! snow melt rate [ m s-1 ???]
    REAL(wp), POINTER          :: precip(:,:)        ! total precipitation rate [mm]
    REAL(wp), POINTER          :: precipinsoil(:,:)  ! accumulated precipitation minus drying (see dust emissions)
    REAL(wp), POINTER          :: precipconv(:,:)    ! convective precipitation at surface [mm]
    REAL(wp), POINTER          :: precipstrat(:,:)   ! convective precipitation at surface [mm]
    REAL(wp), POINTER          :: sw_flux_surf(:,:)  ! shortwave radiation flux at surface [W m-2]
    REAL(wp), POINTER          :: sw_flux_toa(:,:)   ! shortwave radiation flux at top of atmosphere [W m-2]
    REAL(wp), POINTER          :: cdnc(:,:,:)        ! cloud droplet number concentration [m-3]
    REAL(wp), POINTER          :: clw(:,:,:)         ! cloud liquid water [kg kg-1]
    REAL(wp), POINTER          :: aclc(:,:,:)        ! cloud fraction [1]
  END TYPE t_vphysc

  ! vphysc_stream
  ! vphyscvars is used to check namelist input for valid names
  INTEGER, PARAMETER           :: nvphyscvars=22 ! s.stadtler
  CHARACTER(LEN=32)            :: vphyscvars(1:nvphyscvars)= &
                                (/'geom1            ', &
                                  'geohm1           ', &
                                  'aphm1            ', &
                                  'grmassm1         ', &  
                                  'grvolm1          ', &  
                                  'grheightm1       ', &  
                                  'rhoam1           ', &  
                                  'trpwmo           ', &  
                                  'trpwmop1         ', &  
                                  'pbl              ', &  
                                  'velo10m          ', &  
                                  'tsw              ', &  
                                  'smelt            ', &  
                                  'precip           ', &  
                                  'precipinsoil     ', &  
                                  'precipconv       ', &  
                                  'precipstrat      ', &
                                  'sw_flux_surf     ', &
                                  'sw_flux_toa      ', &
                                  'cdnc             ', &
                                  'clw              ', &
                                  'aclc             '  /)


  ! variable pointers
  TYPE(t_vphysc)               :: vphysc


  CONTAINS

  SUBROUTINE init_vphysc_stream

    USE mo_control,             ONLY: nlev
    USE mo_submodel_streams,    ONLY: vphysc_lpost, vphysc_tinterval, vphyscnam
    USE mo_string_utls,         ONLY: st1_in_st2_proof
    USE mo_util_string,         ONLY: tolower
    USE mo_exception,           ONLY: finish
    USE mo_memory_base,         ONLY: t_stream, new_stream, &
                                      default_stream_setting, &
                                      add_stream_element, &
                                      AUTO, BELOWSUR
    ! local variables
    INTEGER, PARAMETER             :: ndefault = 9
    CHARACTER(LEN=32)              :: defnam(1:ndefault)   = &   ! default output
                                (/'geom1            ', &
                                  'geohm1           ', &
                                  'aphm1            ', &
                                  'grmassm1         ', &
                                  'trpwmo           ', &
                                  'pbl              ', &
                                  'precipconv       ', &
                                  'precipstrat      ', &
                                  'velo10m          ' /)
 
    TYPE (t_stream), POINTER       :: vphysc_stream
    INTEGER                        :: ierr
    LOGICAL                        :: lpost

    !-- handle ALL and DEFAULT options
    IF (TRIM(tolower(vphyscnam(1))) == 'all')     vphyscnam(1:nvphyscvars) = vphyscvars(:)
    IF (TRIM(tolower(vphyscnam(1))) == 'default') vphyscnam(1:ndefault) = defnam(:)

    !-- check that all variable names from namelist are valid
    IF (.NOT. st1_in_st2_proof( vphyscnam, vphyscvars, ierr=ierr) ) THEN
      IF (ierr > 0) CALL finish ( 'ini_vphysc_stream', 'variable '// &
                                  vphyscnam(ierr)//' does not exist in vphysc stream' )
    END IF

    !-- open new stream
    !SF #383: registering vphysc to rerun storage
    CALL new_stream (vphysc_stream,'vphysc',lpost=vphysc_lpost,lrerun=.TRUE., &
         interval=vphysc_tinterval)
    CALL default_stream_setting (vphysc_stream, lrerun = .TRUE., &
         contnorest = .TRUE., table = 199, &
         laccu = .false., code = AUTO)

    !-- add individual variables to stream
    lpost = st1_in_st2_proof( 'geom1', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'geom1', vphysc%geom1, &
         longname = 'geopotential height', &
         units = 'm3 s-2', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'geohm1', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'geohm1', vphysc%geohm1, &
         klev=nlev+1,                                                &
         longname = 'geopotential height at interfaces', &
         units = 'm3 s-2', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'aphm1', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'aphm1', vphysc%aphm1, &
         klev=nlev+1,                                                &
         longname = 'atmospheric pressure at interfaces', &
         units = 'Pa', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'grmassm1', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'grmassm1', vphysc%grmassm1, &
         longname = 'air mass in grid box', &
         units = 'kg', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'grvolm1', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'grvolm1', vphysc%grvolm1, &
         longname = 'volume of grid box', &
         units = 'm3', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'grheightm1', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'grheightm1', vphysc%grheightm1, &
         longname = 'height of grid box', &
         units = 'm', lpost=lpost)
    
    lpost = st1_in_st2_proof( 'rhoam1', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'rhoam1', vphysc%rhoam1, &
         longname = 'air density', &
         units = 'kg m-3', lpost=lpost)
    
    lpost = st1_in_st2_proof( 'trpwmo', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'trpwmo', vphysc%trpwmo, &
         longname = 'upper tropopause level index', &
         units = 'levels', lpost=lpost)
    
    lpost = st1_in_st2_proof( 'trpwmop1', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'trpwmop1', vphysc%trpwmop1, &
         longname = 'lower tropopause level index', &
         units = 'levels', lpost=lpost)
    
    lpost = st1_in_st2_proof( 'pbl', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'pbl', vphysc%pbl, &
         longname = 'planetary boundary layer top level', &
         units = 'levels', lpost=lpost)

    lpost = st1_in_st2_proof( 'velo10m', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'velo10m', vphysc%velo10m, &
         longname = 'wind velocity at 10m height', &
         units = 'm s-1', lpost=lpost)

    lpost = st1_in_st2_proof( 'tsw', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'tsw', vphysc%tsw, &
         longname = 'sea water temperature', &
         units = 'K', lpost=lpost)

    lpost = st1_in_st2_proof( 'smelt', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'smelt', vphysc%smelt, &
         longname = '', &
         units = 'm s-1', lpost=lpost)
   
    lpost = st1_in_st2_proof( 'precip', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'precip', vphysc%precip, &
         longname = 'precipitation at surface', &
         units = 'mm', lpost=lpost)
    
    lpost = st1_in_st2_proof( 'precipinsoil', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'precipinsoil', vphysc%precipinsoil, &
         longname = 'precipitation penetrated into soil layers', &
         units = 'mm', leveltype=BELOWSUR, lpost=lpost)
    
    lpost = st1_in_st2_proof( 'precipconv', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'precipconv', vphysc%precipconv, &
         longname = 'surface convective precipitation flux (rain+snow)', &
         units = 'mm', lpost=lpost)

    lpost = st1_in_st2_proof( 'precipstrat', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'precipstrat', vphysc%precipstrat, &
         longname = 'surface stratiform precipitation flux (rain+snow)', &
         units = 'mm', lpost=lpost)

    lpost = st1_in_st2_proof( 'sw_flux_surf', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'sw_flux_surf', vphysc%sw_flux_surf, &
         longname = 'shortwave radiation flux at surface', &
         units = 'W m-2', lpost=lpost)

    lpost = st1_in_st2_proof( 'sw_flux_toa', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'sw_flux_toa', vphysc%sw_flux_toa, &
         longname = 'shortwave radiation flux at top of atmosphere', &
         units = 'W m-2', lpost=lpost)

    lpost = st1_in_st2_proof( 'cdnc', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'cdnc', vphysc%cdnc, &
         longname = 'cloud droplet number concentration', &
         units = 'm-3', lpost=lpost)
    
    lpost = st1_in_st2_proof( 'clw', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'clw', vphysc%clw, &
         longname = 'cloud liquid water', &
         units = 'kg kg-1', lpost=lpost)

    lpost = st1_in_st2_proof( 'aclc', vphyscnam)
    CALL add_stream_element (vphysc_stream, 'aclc', vphysc%aclc, &
         longname = 'cloud fraction', &
         units = '1', lpost=lpost)
  END SUBROUTINE init_vphysc_stream

  !! set_vphysc_var: set the value of a vphysc stream component. 
  !! use klev = -1 if you pass a 2-D variable
  !! It is recommended to use set_vphysc_var if you set several variables at
  !! once or if some computations are necessary to obtain the vphysc variable
  !! from other variables. For simple settings, use vphysc%VARNAME(1:kproma,krow) = zVARNAME(1:kproma).

  SUBROUTINE set_vphysc_var(kproma,        klev,           krow,      &
                            paphm1,        papm1,          ptvm1,     &
                            ppbl,          ktrpwmo,        ktrpwmop1, &
                            prflconv,      psflconv,                  &
                            prflstrat,     psflstrat,      pcdnc,     &
                            pclw,          paclc        )! s.stadtler

    USE mo_physical_constants,    ONLY: rgrav, rd, rhoh2o
    USE mo_geoloc,       ONLY: gboxarea
    USE mo_time_control, ONLY: delta_time

    INTEGER,  INTENT(in) :: kproma, klev, krow
    REAL(wp), INTENT(in), OPTIONAL :: paphm1(:,:), papm1(:,:), ptvm1(:,:), ppbl(:) 
    REAL(wp), INTENT(in), OPTIONAL :: prflconv(:), psflconv(:), prflstrat(:), psflstrat(:)
    REAL(wp), INTENT(in), OPTIONAL :: pcdnc(:,:), pclw(:,:), paclc(:,:) ! s.stadtler
    INTEGER,  INTENT(in), OPTIONAL :: ktrpwmo(:), ktrpwmop1(:)

    LOGICAL, SAVE :: linit = .TRUE.

    INTEGER :: jlev

    IF (linit) THEN
      IF (PRESENT(paphm1)) vphysc%grmassm1(:,:,krow) = 0.0_wp
      IF (PRESENT(papm1) .AND. PRESENT(ptvm1)) THEN
         vphysc%rhoam1(:,:,krow) = 0.0_wp
         IF (PRESENT(paphm1)) THEN
            vphysc%grvolm1(:,:,krow)    = 0.0_wp
            vphysc%grheightm1(:,:,krow) = 0.0_wp
         END IF
      END IF
      IF (PRESENT(ppbl))       vphysc%pbl(:,krow) = 0.0_wp
      IF (PRESENT(ktrpwmo))    vphysc%trpwmo(:,krow) = 0.0_wp
      IF (PRESENT(ktrpwmop1))  vphysc%trpwmop1(:,krow) = 0.0_wp
      IF (PRESENT(prflconv)  .AND. PRESENT(psflconv))  vphysc%precipconv(:,krow) = 0.0_wp 
      IF (PRESENT(prflstrat) .AND. PRESENT(psflstrat)) vphysc%precipstrat(:,krow) = 0.0_wp 
      IF (PRESENT(pcdnc)) vphysc%cdnc(:,:,krow) = 0.0_wp 
      IF (PRESENT(pclw)) vphysc%clw(:,:,krow) = 0.0_wp 
      IF (PRESENT(paclc)) vphysc%aclc(:,:,krow) = 0.0_wp 
      linit = .false.
    ENDIF

    !-- grid box mass
    IF (PRESENT(paphm1)) THEN
       DO jlev = 1, klev
          vphysc%grmassm1(1:kproma,jlev,krow) = (paphm1(1:kproma,jlev+1)-paphm1(1:kproma,jlev)) &
                                         *gboxarea(1:kproma)*rgrav
       END DO
    END IF
    !-- grid box volume and height
    IF (PRESENT(papm1) .AND. PRESENT(ptvm1)) THEN
       DO jlev = 1, klev
          vphysc%rhoam1(1:kproma,jlev,krow) = papm1(1:kproma,jlev)/(rd * ptvm1(1:kproma,jlev))
       END DO
       IF (PRESENT(paphm1)) THEN
          DO jlev = 1, klev
             vphysc%grvolm1(1:kproma,jlev,krow) = &
                       vphysc%grmassm1(1:kproma,jlev,krow)/vphysc%rhoam1(1:kproma,jlev,krow)
             vphysc%grheightm1(1:kproma,jlev,krow) = vphysc%grvolm1(1:kproma,jlev,krow)/gboxarea(1:kproma)
          END DO
       END IF
    END IF
    !-- PBL height
    IF (PRESENT(ppbl)) vphysc%pbl(1:kproma,krow) = ppbl(1:kproma)
    !-- tropopause height
    IF (PRESENT(ktrpwmo)) vphysc%trpwmo(1:kproma,krow) = REAL(ktrpwmo(1:kproma),wp)
    IF (PRESENT(ktrpwmop1)) vphysc%trpwmop1(1:kproma,krow) = &
                            REAL(ktrpwmop1(1:kproma),wp)
    !-- precipitation
    ! convective precip comes first; compute total when stratiform precip is given
    IF (PRESENT(prflconv) .AND. PRESENT(psflconv)) THEN
       vphysc%precipconv(1:kproma,krow)= &
                      ((prflconv(1:kproma)+psflconv(1:kproma))/rhoh2o)*delta_time*1.0E+3_wp
       vphysc%precip(1:kproma,krow)=vphysc%precipconv(1:kproma,krow)
    END IF
    IF (PRESENT(prflstrat) .AND. PRESENT(psflstrat)) THEN
      vphysc%precipstrat(1:kproma,krow)= &
                      ((prflstrat(1:kproma)+psflstrat(1:kproma))/rhoh2o)*delta_time*1.0E+3_wp
      vphysc%precip(1:kproma,krow)=vphysc%precip(1:kproma,krow)+vphysc%precipstrat(1:kproma,krow)
    END IF
    
    ! cloud droplet number concentration cdnc
    IF (PRESENT(pcdnc)) THEN
      DO jlev =1, klev
        vphysc%cdnc(1:kproma,jlev,krow) = pcdnc(1:kproma,jlev) ! s.stadtler 
      END DO
    END IF
    IF (PRESENT(pclw)) THEN
      DO jlev =1, klev
        vphysc%clw(1:kproma,jlev,krow) = pclw(1:kproma,jlev) ! s.stadtler 
      END DO
    END IF
    IF (PRESENT(paclc)) THEN
      DO jlev =1, klev
        vphysc%aclc(1:kproma,jlev,krow) = paclc(1:kproma,jlev) ! s.stadtler 
      END DO
    END IF

  END SUBROUTINE set_vphysc_var

END MODULE mo_vphysc
