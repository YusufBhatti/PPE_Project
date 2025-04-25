!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_sedimentation.f90
!!
!! \brief
!! module to interface ECHAM submodules with sedimentation module(s)
!!
!! \author M. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2009-10-26)
!!   -# M. Schultz (FZ Juelich) - improved diag routines (2010-04-16)
!!
!! \limitations
!! All diag_lists must be defined in order to avoid problems with
!! get_diag_pointer in the actual sedi_interface routine. Lists can be empty.
!! Currently there is only one unified interactive sedimentation scheme for
!! aerosols (HAM).
!!
!! \details
!! This module initializes the scheme based on the namelist parameters
!! in submodeldiagctl and creates a stream for variable pointers and 
!! diagnostic quantities used in the sedimentation scheme. It also
!! provides a generic interface to the actual sedimentation routine(s).
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
MODULE mo_hammoz_sedimentation

  USE mo_kind,             ONLY: dp
  USE mo_submodel_diag,    ONLY: t_diag_list

  IMPLICIT NONE

  PRIVATE

  ! public variables  (see declaration below)

  ! subprograms
  PUBLIC                       :: init_sedi_stream, &
                                  sedi_interface

  ! sedi_stream
  INTEGER, PARAMETER           :: nsedivars=2
  CHARACTER(LEN=32)            :: sedivars(1:nsedivars)= &
                                (/'sed             ', &   ! total sedimentation flux
                                  'vsedi           '  /)  ! sedimentation velocity

  ! variable pointers and diagnostic lists
  TYPE (t_diag_list), PUBLIC   :: sed         ! sedimentation flux
  TYPE (t_diag_list), PUBLIC   :: vsedi       ! sedimentation velocity


  CONTAINS

  SUBROUTINE init_sedi_stream

    USE mo_string_utls,         ONLY: st1_in_st2_proof
    USE mo_util_string,         ONLY: tolower
    USE mo_exception,           ONLY: finish
    USE mo_memory_base,         ONLY: t_stream, new_stream, &
                                      default_stream_setting, &
                                      add_stream_reference, &
                                      AUTO
    USE mo_ham_m7_trac,         ONLY: ham_get_class_flag
    USE mo_tracdef,             ONLY: ln, ntrac, trlist, AEROSOL
    USE mo_species,             ONLY: nspec, speclist
    USE mo_ham,                 ONLY: nclass
    USE mo_submodel_streams,    ONLY: sedi_lpost, sedi_tinterval, sedinam, sedi_keytype
    USE mo_submodel_diag,       ONLY: new_diag_list, new_diag,   &
                                      BYTRACER, BYSPECIES, BYNUMMODE, BYMODE !SF #299 added BYMODE
    USE mo_submodel,            ONLY: lham !SF, see #228


    ! local variables
    INTEGER, PARAMETER             :: ndefault = 2
    CHARACTER(LEN=32)              :: defnam(1:ndefault)   = &   ! default output
                                (/ 'sed             ', &    ! total sedimentation flux
                                   'vsedi           '  /)   ! sedimentation velocity
 
    LOGICAL                        :: tracflag(ntrac), specflag(nspec), modflag(MAX(nclass,1))
    CHARACTER(LEN=ln)              :: tracname(ntrac), specname(nspec), modname(MAX(nclass,1)), &
                                      modnumname(MAX(nclass,1)) !SF #299
    TYPE (t_stream), POINTER       :: ssedi
    INTEGER                        :: ierr, jt
    LOGICAL                        :: lpost

    !++mgs: default values and namelist read are done in init_submodel_streams !

    !-- handle ALL, DETAIL and DEFAULT options for sedi output variables
    !-- Note: ALL and DETAIL are identical for sedi output
    IF (TRIM(tolower(sedinam(1))) == 'detail')  sedinam(1:nsedivars) = sedivars(:)
    IF (TRIM(tolower(sedinam(1))) == 'all')     sedinam(1:nsedivars) = sedivars(:)
    IF (TRIM(tolower(sedinam(1))) == 'default') sedinam(1:ndefault) = defnam(:)

    !-- check that all variable names from namelist are valid
    IF (.NOT. st1_in_st2_proof( sedinam, sedivars, ierr=ierr) ) THEN
      IF (ierr > 0) CALL finish ( 'ini_sedi_stream', 'variable '// &
                                  sedinam(ierr)//' does not exist in sedi stream' )
    END IF

    !-- define the flags and names for the diagnostic lists. We need one set of flags and
    !   names for each key_type (BYTRACER, BYSPECIES, BYMODE)
    !   gas-phase tracers will always be defined BYTRACER, for aerosol tracers one of the
    !   following lists will be empty.
    tracflag(:) = .FALSE.
    DO jt = 1,ntrac
      tracname(jt) = trlist%ti(jt)%fullname
      IF (IAND(trlist%ti(jt)%nphase,AEROSOL) /= 0 .AND.       &  !>>dod diagnostics bugfix <<dod
          sedi_keytype == BYTRACER .AND.              &
          nclass > 0) THEN
        tracflag(jt) = trlist%ti(jt)%nsedi > 0
      END IF
    END DO
    specflag(:) = .FALSE.
    DO jt = 1,nspec
      specname(jt) = speclist(jt)%shortname
      IF (sedi_keytype == BYSPECIES .AND.                       &
          IAND(speclist(jt)%nphase, AEROSOL) /= 0 .AND.         &
          nclass > 0) THEN
        specflag(jt) = .TRUE.
      END IF
    END DO
    modflag(:) = .false.
    modname(:) = ''
    !SF #228, adding a condition to check that HAM is active:
    !SF #299, adding a condition to check if BYMODE is relevant:
    IF (lham .AND. nclass > 0 .AND. (sedi_keytype == BYMODE)) &
       CALL ham_get_class_flag(nclass, modflag, modname, modnumname, lsedi=.TRUE.) ! get all modes

    !-- open new stream
    CALL new_stream (ssedi,'sedi',lpost=sedi_lpost,lrerun=.FALSE., &
         interval=sedi_tinterval)
    CALL default_stream_setting (ssedi, lrerun = .FALSE., &
         contnorest = .TRUE., table = 199, &
         laccu = .false., code = AUTO)
   
    !-- add standard ECHAM variables
    IF (sedi_lpost) THEN
      CALL add_stream_reference (ssedi, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
      CALL add_stream_reference (ssedi, 'lsp'     ,'sp'    ,lpost=.TRUE.)
      CALL add_stream_reference (ssedi, 'aps'     ,'g3b'   ,lpost=.TRUE.)
      CALL add_stream_reference (ssedi, 'gboxarea','geoloc',lpost=.TRUE.)
    END IF

    !-- instantaneous diagnostic quantities

    CALL default_stream_setting (ssedi,  lrerun=.FALSE., laccu=.FALSE.)

    !-- sedimentation velocities by tracer or by aerosol mode
    lpost = st1_in_st2_proof( 'vsedi', sedinam) .AND. sedi_lpost
    CALL new_diag_list (vsedi, ssedi, diagname='vsedi', tsubmname='',     &
                        longname='sedimentation velocity', units='m s-1', &
                        ndims=2, nmaxkey=(/ntrac, 0, nclass, nclass, 0 /), lpost=lpost )
    CALL new_diag(vsedi, ntrac, tracflag, tracname, BYTRACER)
    IF (ANY(modflag)) THEN !SF #299 added mode mass diags and fix for mode number name
      CALL new_diag(vsedi, nclass, modflag, modname, BYMODE)
      CALL new_diag(vsedi, nclass, modflag, modnumname, BYNUMMODE)
    END IF

    !-- average diagnostic quantities

    CALL default_stream_setting (ssedi, lrerun=.FALSE., laccu=.TRUE.)

    !-- total sedi flux
    lpost = st1_in_st2_proof( 'sed', sedinam) .AND. sedi_lpost
    CALL new_diag_list (sed, ssedi, diagname='sed', tsubmname='',    &
                        longname='accumulated sedimentation flux',&
                        units='kg m-2 s-1', ndims=2,                    &
                        nmaxkey=(/ntrac, nspec, nclass, nclass, 0 /), lpost=lpost )
    ! add diagnostic elements only when output is activated
    IF (lpost) THEN
      CALL new_diag(sed, ntrac, tracflag, tracname, BYTRACER)
      CALL new_diag(sed, nspec, specflag, specname, BYSPECIES)
      IF (ANY(modflag)) THEN !SF #299 added mode mass diags and fix for mode number name
        CALL new_diag(sed, nclass, modflag, modname, BYMODE)
        CALL new_diag(sed, nclass, modflag, modnumname, BYNUMMODE)
      END IF
    END IF

  END SUBROUTINE init_sedi_stream

  !! ---------------------------------------------------------------------------------------
  !! sedi_interface: generic interface routine to sedimentation
  !! currently, the only sedimentation scheme implemented is that of HAM

  SUBROUTINE sedi_interface(kbdim, kproma, klev, krow,      & 
                            pt,    pq,                      &
                            pap,   paph,                    &
                            pxtm1,  pxtte               )


  USE mo_tracdef,              ONLY: ntrac, trlist
  USE mo_time_control,         ONLY: time_step_len, delta_time
  USE mo_submodel_diag,        ONLY: get_diag_pointer
  USE mo_ham_sedimentation,    ONLY: ham_prep_sedi, ham_sedimentation

  !--- parameters
  INTEGER,  INTENT(in)    :: kbdim, kproma, klev, krow
  REAL(dp), INTENT(in)    :: pt(kbdim, klev),         & ! temperature
                             pq(kbdim, klev),         & ! specific humidity 
                             pap(kbdim, klev),        & ! full level pressure
                             paph(kbdim, klev+1),     & ! half level pressure
                             pxtm1(kbdim,klev,ntrac)    ! tracer mass/number mixing ratio
  REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ntrac)    ! tracer tendency

  !--- local variables
  INTEGER       :: jt, ierr
  REAL(dp)      :: ztempc(kbdim, klev),   &  ! temp. above melting
                   zvis(kbdim, klev),     &  ! air viscosity
                   zlair(kbdim, klev),    &  ! mean free path
                   zrho(kbdim, klev),     &  ! air density
                   zdpg(kbdim, klev),     &  ! layer thickness (pressure)
                   zdz(kbdim, klev),      &  ! layer thickness (length)
                   zxtp1(kbdim, klev),    &  ! updated tracer(jt) 
                   zxtte(kbdim, klev),    &  ! tracer(jt) tendency
                   zvsedi(kbdim, klev),   &  ! sedimentation velocity
                   zsediflux(kbdim, klev)    ! sedimentation flux
  REAL(dp), POINTER    :: fld2d(:,:)         ! pointer for diagnostics

  !--- calculate tracer independent physical variables
  ! note: IF (ANY(trlist%ti(:)%nsedi)/=0) THEN  not needed -- this is checked upon calling
  CALL ham_prep_sedi(kproma, kbdim, klev, &
                     pt,     pq,          &
                     pap,    paph,        &
                     ztempc, zvis,        &
                     zlair,  zrho,        &
                     zdpg,   zdz          )

  !--- tracer loop
  DO jt=1, ntrac
    IF (trlist%ti(jt)%nsedi==0) CYCLE    ! do nothing if tracer doesn't sediment

    !--- update tracer concentration
    zxtp1(1:kproma,:) = pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt) * time_step_len
    zxtte(1:kproma,:) = pxtte(1:kproma,:,jt)

    !--- call HAM sedimentation routine    
    CALL ham_sedimentation(kproma, kbdim, klev, krow, &
                           jt, zvis, zlair, zrho,     &
                           zdpg, zdz,                 & 
                           zxtp1, zxtte,              &
                           zvsedi, zsediflux)

    !--- store diagnostics
    CALL get_diag_pointer(sed, fld2d, jt, ierr=ierr)
    IF (ierr == 0) fld2d(1:kproma,krow) = fld2d(1:kproma,krow)      &
                                          + zsediflux(1:kproma,klev)*delta_time
    CALL get_diag_pointer(vsedi, fld2d, jt, ierr=ierr)
    IF (ierr == 0) fld2d(1:kproma,krow) = zvsedi(1:kproma,klev)

    pxtte(1:kproma,:,jt)=zxtte(1:kproma,:)  !csld write tendency out

  END DO

  END SUBROUTINE sedi_interface


END MODULE mo_hammoz_sedimentation
