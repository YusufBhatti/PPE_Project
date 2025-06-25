!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file defines and sets the control flags for the HAM, MOZ and 
!! HAMMONIA submodels including the settings that control the HAM and 
!! MOZ coupling and the switches for MEGAN, AIRCRAFT (emissions) and 
!! LIGHTNING.
!!
!! These switches are controlled by the *submodelctl* namelist.
!! This namelist defines general submodel switches which are needed in 
!! the interface layer (*mo_submodel_interface* or other parts of the 
!! standard ECHAM code). Also included are switches which define the 
!! coupling between various submodels (for example in HAMMOZ, which coupled 
!! aerosol and gas-phase chemical processes). Other submodel-specific switches 
!! should be defined and maintained in extra namelists which carry the name 
!! of the submodel itself (e.g. mozctl or hamctl).
!!
!! Change this file to attach different or additional submodels to ECHAM.
!!
!!
!! @author 
!! <ol>
!! <li>M. Schultz (FZ-Juelich)
!! <li>S. Rast    (MPI-Met)
!! <li>K. Zhang   (MPI-Met)
!! </ol>
!!
!! $Id: 1423$
!!
!! @par Revision History
!! <ol>
!! <li>M. Schultz   (FZ-Juelich) 
!!     -  original idea and code structure (2009-05-xx) 
!! <li>S. Rast      (MPI-Met)    
!!     -  original idea and code structure (2009-06-xx) 
!! <li>K. Zhang     (MPI-Met)    
!!     -  restucture and new style, implementation in ECHAM6 (2009-07-xx)
!! <li>L. Kornblueh (MPI-Met)
!!     -  remove spitfire (2012-02-xx)
!! </ol>
!!
!! @par This module is used by
!! ...
!! 
!! @par Notes
!! mgs: new module - contains module variables for submodels and submodel 
!! stuff from old mo_tracer
!! linterh2o: info from Hauke: even if linterh2o=true, output of H2O as tracer may differ from 
!! ECHAM's sh output. Reason unclear.
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  MODULE mo_submodel

  USE mo_tracdef,        ONLY: ln     ! length of caracter string (module name)

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: submlist                       ! submodel list  (formerly in mo_tracdef)
  PUBLIC :: t_submlist                     ! submodel list element data type  (formerly in mo_tracdef)
  PUBLIC :: nsubm                          ! number of submodel list entries  (formerly in mo_tracdef)
  PUBLIC :: new_submodel                   ! define a new submodel
  PUBLIC :: setsubmodel                    ! read namelist 
  PUBLIC :: query_submodel                 ! check if a submodel was defined
  PUBLIC :: starttracdef                   ! flag first tracer index of a submodel
  PUBLIC :: endtracdef                     ! flag last tracer index of a submodel
  
  PUBLIC :: print_status                   ! utility routines
  PUBLIC :: print_value                    ! utility routines
  
  PUBLIC :: lanysubmodel                   ! general submodel flag for fast testing
  PUBLIC :: lchemistry                     ! general process switches
  PUBLIC :: lemissions                     ! switch for emission diagnostics
  PUBLIC :: ldrydep                        ! switch for dry diagnostics
  PUBLIC :: lwetdep                        ! switch for wet diagnostics
  PUBLIC :: lsedimentation                 ! switch for sedimentation
  PUBLIC :: laero_micro                    ! switch for microphysics
  PUBLIC :: lburden                        ! switch for burden diagnostics
  
! central switches for submodels
  PUBLIC :: lmethox                        ! switch for upper atmospheric H2O production from methane
  PUBLIC :: lco2                           ! switch on/off CO2 transport submodel
  PUBLIC :: ltransdiag                     ! switch for turning on atmospheric energy transport diagnostic
                                           ! submodel switches:
  PUBLIC :: lhammoz                        ! HAMMOZ gas-phase and aerosol chemistry (lmoz+lham)
  PUBLIC :: lham                           ! HAM aerosol module
  PUBLIC :: lmoz                           ! MOZ chemistry module
  PUBLIC :: lhammonia                      ! HAMMOZ with upper atmosphere extensions
  PUBLIC :: llght                          ! (MOZ) lightning module (can run independently)
  PUBLIC :: lbioemi_stdalone               ! Biogenic emissions as a standalone submodel
  PUBLIC :: lxt                            ! simple tracer module [yet to be completed]
  PUBLIC :: losat                          ! satellite simulator switches
  PUBLIC :: loisccp                        ! ISCCP simulator switches
  PUBLIC :: lccnclim                       ! switch for CCN climatology (as submodel)
  PUBLIC :: lflighttrack                   ! switch for flight-track simulator
  
  PUBLIC :: lhmzphoto                      ! HAMMOZ coupling switches
  PUBLIC :: lhmzoxi                        !
  PUBLIC :: lhmzhet                        ! 
  PUBLIC :: lhmzhetwet                     !
  PUBLIC :: lhmzhetdust                    !
  
  PUBLIC :: lchemfeedback                  ! Switches for chemistry feedbacks on ECHAM physics
  PUBLIC :: lchemrad                       !
  PUBLIC :: linterh2o                      ! feedback water content from MOZ to ECHAM and vice versa
  PUBLIC :: lchemheat                      ! HAMMONIA coupling switches
  PUBLIC :: linteram                       !
  PUBLIC :: lintercp                       !
  PUBLIC :: laoa                           ! age-of-air submodel switch

  PUBLIC :: emi_basepath
  PUBLIC :: emi_scenario !sschr See #411 (HAMMOZ)

  PUBLIC :: id_xt
  PUBLIC :: id_ham 
  PUBLIC :: id_bioemi
  PUBLIC :: id_moz 
  PUBLIC :: id_hammonia
  PUBLIC :: id_lightning
  PUBLIC :: id_isccp
  PUBLIC :: id_sat
  PUBLIC :: id_hrates
  PUBLIC :: id_ccnclim
  PUBLIC :: id_flighttrack
  PUBLIC :: id_aoa
  
  !
  ! interfaces
  !                                
  INTERFACE print_value            ! report on a parameter value  
    MODULE PROCEDURE print_lvalue  ! logical
    MODULE PROCEDURE print_ivalue  ! integer
    MODULE PROCEDURE print_rvalue  ! real
  END INTERFACE
  !
  ! Type declarations
  !
  INTEGER, PARAMETER :: ns = 20 ! max number of submodels

  TYPE t_submlist
    CHARACTER(len=ln) :: modulename ! name of sub-model
    INTEGER           :: idtfirst   ! id of first tracer defined for this submodel
    INTEGER           :: idtlast    ! id of last tracer defined for this submodel
  END TYPE t_submlist

  
  !
  ! module variables
  !

  TYPE(t_submlist) ,SAVE         :: submlist (ns) ! submodel list
  INTEGER          ,SAVE         :: nsubm = 0     ! number of submodels defined

! submodel master switch: default false, set to true if any submodel registers
  LOGICAL :: lanysubmodel = .FALSE.

! general chemical process interface switches for debugging purposes:
! Default is true - will be turned off if no submodel is defined (*mo_submodel_interface*)
! or if switch is turned off in SUBMODELCTL namelist
  LOGICAL :: lemissions     = .TRUE. 
  LOGICAL :: lchemistry     = .TRUE.
  LOGICAL :: ldrydep        = .TRUE.
  LOGICAL :: lwetdep        = .TRUE. 
  LOGICAL :: lsedimentation = .TRUE.
  LOGICAL :: laero_micro    = .TRUE.
!!$  LOGICAL :: lburden        = .TRUE.    ! .true. for burden diagnostics 
!!$                                    ! (defaults to false if neither HAM nor MOZ are active)
  LOGICAL :: lburden        = .FALSE. ! has to be set by namelist  
  
! submodel switches (default: all submodels turned off)
  LOGICAL :: lxt               = .FALSE. ! .true. to activate simple generic tracer submodel 
  LOGICAL :: lmethox           = .FALSE. ! .true. for upper atmospheric H2O production from methane
  LOGICAL :: lco2              = .FALSE. ! .true. for interactive transport CO2 subm.
  LOGICAL :: ltransdiag        = .FALSE. ! .true. for atmospheric energy transport diagnostic
  LOGICAL :: lhammoz           = .FALSE. ! .true. to turn on HAM and MOZ and activate coupling
  LOGICAL :: lham              = .FALSE. ! .true. for aerosol module HAM
  LOGICAL :: lmoz              = .FALSE. ! .true. for gas-phase chemistry module MOZ
  LOGICAL :: llght             = .FALSE. ! .true. for enabling lightning emissions
  LOGICAL :: lbioemi_stdalone  = .FALSE. ! .true. for enabling biogenic emissions as
                                         ! a standalone submodel
  ! Diagnostics submodel switches
  ! Satellite and ISCCP cloud diagnostic processors
  LOGICAL :: losat          = .FALSE. ! .true. for satellite simulator
  LOGICAL :: loisccp        = .FALSE. ! .true. for ISCCP diagnostics processor
  LOGICAL :: lccnclim       = .FALSE. ! .true. for CCN climatology
  LOGICAL :: lflighttrack   = .FALSE. ! .true. for flight-track simulator

  
! HAMMOZ coupling switches (default: true -- they are automatically turned off if one module is inactive)
! see also Liao et al, JGR, 2005, table 1.
  LOGICAL :: lhmzphoto  = .TRUE.  ! switch HAMMOZ photolysis coupling off
  LOGICAL :: lhmzoxi    = .TRUE.  ! switch HAMMOZ oxidant coupling off
  LOGICAL :: lhmzhet    = .TRUE.  ! switch on/off all the HAMMOZ heterogenous reactions
  LOGICAL :: lhmzhetwet = .TRUE.  ! switch on/off het. reactions on SU,SS,BC,OC and wet aerosol
  LOGICAL :: lhmzhetdust= .TRUE.  ! switch on/off het. reactions on mineral dust
  
! chemistry feedback switches 
  LOGICAL :: lchemfeedback = .FALSE. ! combined switch for linterh2o, lchemheat, linteram, lintercp + radiation
  LOGICAL :: lchemrad   = .FALSE. ! switch chemistry feedback on radiation
  LOGICAL :: linterh2o  = .FALSE. ! switch return of water vapour from chemistry to ECHAM physics
! HAMMONIA switches
  LOGICAL :: lhammonia  = .FALSE. ! switch HAMMONIA (not possible in this version)
  LOGICAL :: lchemheat  = .FALSE. ! switch chemical heating (HAMMONIA only)
  LOGICAL :: linteram   = .FALSE. ! switch calculation of air mass from chemistry (HAMMONIA only)
  LOGICAL :: lintercp   = .FALSE. ! switch calculation of cp from chemistry (HAMMONIA only)
  LOGICAL :: laoa       = .FALSE. ! .true. for enabling age-of-air submodel
  
  CHARACTER(LEN=256) :: emi_basepath = ''
  CHARACTER(LEN=8)   :: emi_scenario = ''

!!mgs&jsr!! questionable if we need these id values. Might be removed later...  
  INTEGER :: id_methox
  INTEGER :: id_co2
  INTEGER :: id_transdiag
  INTEGER :: id_xt
  INTEGER :: id_ham 
  INTEGER :: id_bioemi
  INTEGER :: id_moz 
  INTEGER :: id_hammonia
  INTEGER :: id_lightning
  INTEGER :: id_isccp
  INTEGER :: id_sat
  INTEGER :: id_hrates
  INTEGER :: id_ccnclim
  INTEGER :: id_flighttrack
  
  INTEGER :: id_aoa
  
  CONTAINS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Check status of a given submodel. Return id of submodel if requested.
!! Function result is logical value indicating proper initialisation.
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Externals
!!
!! @par Notes
!! not yet used
!!  
!! @par Responsible coder
!! M.Schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  LOGICAL FUNCTION query_submodel (modulename, id)

  
  USE mo_util_string,    ONLY : toupper

  CHARACTER(len=*) ,INTENT(in)            :: modulename   ! name of submodel
  INTEGER          ,INTENT(out) ,OPTIONAL :: id           ! index of submodel, 0 if undefined
 
  !-- local variables --
 
  INTEGER           :: i
    
  IF (PRESENT(id)) id = 0
    
  query_submodel = .false.
  DO i=1, nsubm
    IF (submlist(i)% modulename == toupper(modulename)) THEN  
       query_submodel = .true.
       RETURN
    END IF
  END DO

  END FUNCTION query_submodel

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Request a new submodel with name 'name'.
!! Name is always stored in uppercase characters
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! setsubmodel
!!
!! @par Externals
!!
!! @par Notes
!!
!! @par Responsible coder
!! M.Schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE new_submodel(name,id)

  USE mo_util_string,    ONLY: toupper
  USE mo_exception,      ONLY: message_text, finish

  CHARACTER (len=*) ,INTENT(in)  :: name
  INTEGER           ,INTENT(out) ,OPTIONAL :: id           ! index of submodel, 0 if undefined

    IF (PRESENT(id)) id = 0
    nsubm = nsubm + 1
    IF (nsubm > ns) THEN
      WRITE(message_text,*) 'Submodel list full. ns = ', ns
      CALL finish ('new_submodel', message_text)
    END IF
    IF (PRESENT(id)) id = nsubm
    submlist (nsubm)% modulename = toupper(name)
    lanysubmodel = .true.
    
  END SUBROUTINE new_submodel


  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! It reads the submodel.ctl namelist to set the submodel control switches.
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! initialize
!!
!! @par Externals: 
!!
!! @par Notes 
!!
!! @par Responsible coder
!! M.Schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE setsubmodel

  
  USE mo_mpi,            ONLY: p_io, p_parallel, p_parallel_io, p_bcast
  USE mo_exception,      ONLY: message, finish, em_info, em_param,         &
                               message_text, em_warn, em_error
  USE mo_namelist,       ONLY: open_nml, position_nml, POSITIONED
  USE mo_util_string,    ONLY: separator   ! format string (----)
  USE mo_advection,      ONLY: iadvec,           &
                               no_advection,     & ! for diagnostic printout
                               semi_lagrangian,  &
                               tpcore
! local variables

  INTEGER :: i, inml, iunit, ierr

  INCLUDE 'submodelctl.inc'


  CALL message('', ' ') 
  CALL message('', '|') 
  CALL message('', '| setsubmodel ') 
  CALL message('', '|') 
  CALL message('', ' ') 
  
  ! initialize submodel list
  
  submlist (:)% modulename = ''
  submlist (:)% idtfirst = 0
  submlist (:)% idtlast  = 0

  IF (p_parallel_io) THEN

    ! read namelist
  
    inml = open_nml('namelist.echam')
    iunit = position_nml ('SUBMODELCTL', inml, status=ierr)
     
    SELECT CASE (ierr)
    CASE (POSITIONED)
      READ (iunit, submodelctl)
    END SELECT

     ! evaluate HAMMONIA switch: exit model if activated. This code does not run HAMMONIA
     
     IF ( lhammonia ) THEN
        CALL finish('setsubmodel', &
                    'LHAMMONIA=.true. --> Abort program execution. This code version ' // &
                    'does not contain HAMMONIA routines.')
     END IF
     
     ! overwrite HAMMONIA flags for safety
     
     lchemheat = .false.
     linteram = .false.
     lintercp = .false.

  ENDIF
  
  IF (p_parallel) THEN

     ! submodel switches
     CALL p_bcast (lmethox, p_io)
     CALL p_bcast (lco2, p_io)
     CALL p_bcast (ltransdiag, p_io)
     CALL p_bcast (lxt, p_io)
     CALL p_bcast (lham, p_io)
     CALL p_bcast (lmoz, p_io)
     CALL p_bcast (lhammoz, p_io)
     CALL p_bcast (lhammonia, p_io)
     CALL p_bcast (llght, p_io)
     CALL p_bcast (lbioemi_stdalone, p_io)
     CALL p_bcast (losat, p_io)
     CALL p_bcast (loisccp, p_io)
     CALL p_bcast (lccnclim, p_io)
     CALL p_bcast (lflighttrack, p_io)

     ! generic process switches
     CALL p_bcast (lchemistry, p_io)
     CALL p_bcast (ldrydep, p_io)
     CALL p_bcast (lwetdep, p_io)
     CALL p_bcast (lemissions, p_io)
     CALL p_bcast (lsedimentation, p_io)
     CALL p_bcast (laero_micro, p_io)

     ! coupling and diagnostic switches
     CALL p_bcast (lhmzphoto, p_io)
     CALL p_bcast (lhmzoxi, p_io)
     CALL p_bcast (lhmzhet, p_io)
     CALL p_bcast (lchemfeedback, p_io)
     CALL p_bcast (lchemrad, p_io)
     CALL p_bcast (linterh2o, p_io)
     CALL p_bcast (lchemheat, p_io)
     CALL p_bcast (linteram, p_io)
     CALL p_bcast (lintercp, p_io)
     CALL p_bcast (lburden, p_io)

     CALL p_bcast (emi_basepath, p_io)
     CALL p_bcast (emi_scenario, p_io)
     
     ! age-of-air tracer switch
     CALL p_bcast (laoa, p_io)
     
  END IF

!sschr: this is the right place for the following commands!

  ! evaluate HAMMOZ switches
  ! hammoz always means moz and ham are on. By default all hammoz couplings are also on,
  ! but you can switch them off in the namelist. The coupling is always off if lhammoz is not 
  ! set.
  
  IF ( lhammoz ) THEN
     lham = .true.
     lmoz = .true.
  ELSE 
     lhmzphoto = .false.
     lhmzoxi   = .false.
     lhmzhet   = .false.
  END IF
  
  ! set individual switches if lchemfeedback is true
  IF ( lchemfeedback .AND. lmoz ) THEN
    lchemrad = .TRUE.
    linterh2o = .TRUE.
    IF ( lhammonia ) THEN
      lchemheat = .TRUE.
      linteram  = .TRUE.
      lintercp  = .TRUE.
    ENDIF
  ENDIF

  ! global flag for burden diagnostics
  
!!$  IF (.NOT. lmoz .AND. .NOT. lham) THEN
!!$     lburden = .false.
!!$  END IF

  ! make sure that chemheat, linteram and lintercp are only set if HAMMONIA is active
  IF ( .NOT. lhammonia ) THEN
     IF ( lchemheat ) CALL message('setsubmodel',     &
             'Switching off lchemheat now, because LHAMMONIA = false.', &
             level=em_warn)
     IF ( linteram ) CALL message('setsubmodel',     &
             'Switching off linteram now, because LHAMMONIA = false.', &
             level=em_warn)
     IF ( lintercp ) CALL message('setsubmodel',     &
             'Switching off lintercp now, because LHAMMONIA = false.', &
             level=em_warn)
  ENDIF

  ! evaluate MOZ and chemistry switches
  
  IF ( .NOT. lmoz .OR. .NOT. lchemistry ) THEN
     IF ( lchemfeedback ) CALL message('setsubmodel',     &
             'Switching off lchemfeedback now, because LMOZ = false. or lchemistry = .false.', &
             level=em_warn)
     IF ( lchemrad ) CALL message('setsubmodel',     &
             'Switching off lchemrad now, because LMOZ = false. or lchemistry = .false.', &
             level=em_warn)
     IF ( linterh2o ) CALL message('setsubmodel',     &
             'Switching off linterh2o now, because LMOZ = false. or lchemistry = .false.', &
             level=em_warn)
     IF ( lchemheat ) CALL message('setsubmodel',     &
             'Switching off lchemheat now, because LMOZ = false. or lchemistry = .false.', &
             level=em_warn)
     IF ( linteram ) CALL message('setsubmodel',     &
             'Switching off linteram now, because LMOZ = false. or lchemistry = .false.', &
             level=em_warn)
     IF ( lintercp ) CALL message('setsubmodel',     &
             'Switching off lintercp now, because LMOZ = false. or lchemistry = .false.', &
             level=em_warn)
     lchemfeedback = .false.
     lchemrad      = .false.
     linterh2o     = .false.
     lchemheat     = .false.
     linteram      = .false.
     lintercp      = .false.
  END IF

  ! make sure that ECHAM methox is not called if MOZART is active and linterh2o=true

  IF ( linterh2o .AND. lmethox ) THEN
    CALL message('setsubmodel',     &
                 'Switching off lmethox, because stratospheric H2O is handled by MOZ and linterh2o=true.', &
                 level=em_warn)
    lmethox = .FALSE.
  ENDIF

  ! evaluate CCNCLIM and HAM switches

  IF ( lccnclim .AND. lham ) THEN
    CALL message('setsubmodel', 'De-activating CCNCLIM, because LHAM=.TRUE.', level=em_warn)
    lccnclim = .false.
  END IF

  ! force chemical water vapour if radiative feedback is activated (???)   !!baustelle!!
!!IF ( linterchem ) THEN
!!   linterh2o = .true.
!!END IF 
  
  ! individual het coupling processes
  ! presently these are not independently controlled by namelist switches
  
  IF ( lhmzhet ) THEN
     lhmzhetwet   = .true.
     lhmzhetdust  = .true.
  ELSE 
     lhmzhetwet   = .false.
     lhmzhetdust  = .false.
  END IF   
  
!!++mgs: potentially need to turn off lhmz switches if lchemistry=false...
  
  ! turn general chemistry process switches off if no submodel has been activated
!++mgs  
  IF (.NOT. (lham .OR. lmoz .OR. llght .OR. lbioemi_stdalone .OR. lco2) ) THEN
!--mgs
     lchemistry     = .FALSE.
     lemissions     = .FALSE.
     ldrydep        = .FALSE.
     lwetdep        = .FALSE.
     lsedimentation = .FALSE.
     laero_micro    = .FALSE.
     CALL message('setsubmodel', &
                  'This run is without any chemical processes (lchemistry=F, lemissions=F'// &
                  ', ldrydep=F, lwetdep=F, lsedimentation=F, laero_micro=F)', level=em_info )
  END IF

!>>SF check for biogenic emissions inconsistency:
!     Biogenic emissions can be both handled by the standalone submodel and by 
!     HAM or MOZ at the same time
  IF ((lham .OR. lmoz) .AND. lbioemi_stdalone) THEN
     CALL message('setsubmodel', &
                  'Biogenic emissions cannot be run as standalone if HAM or MOZ are active! '// &
                  'Try and turn off lbioemi_stdalone or lham and/or lmoz', level=em_error)
  ENDIF
!<<SF

  ! Register submodels 
  
  IF (lmethox)          CALL new_submodel('METHOX',           id_methox   )
  IF (lco2)             CALL new_submodel('CO2',              id_co2      )
  IF (ltransdiag)       CALL new_submodel('TRANSDIAG',        id_transdiag)
  IF (lxt)              CALL new_submodel('XT',               id_xt       )
  IF (lham)             CALL new_submodel('HAM',              id_ham      )
  IF (lmoz)             CALL new_submodel('MOZ',              id_moz      )
  IF (lhammonia)        CALL new_submodel('HAMMONIA',         id_hammonia )
  IF (llght)            CALL new_submodel('LIGHTNING',        id_lightning)
  IF (lbioemi_stdalone) CALL new_submodel('BIOEMI_STANDALONE',id_bioemi   ) 
  IF (losat)            CALL new_submodel('SAT',              id_sat      )
  IF (loisccp)          CALL new_submodel('ISCCP',            id_isccp    )
  IF (lccnclim)         CALL new_submodel('CCNCLIM',          id_ccnclim  )
  IF (lflighttrack)     CALL new_submodel('FLIGHTTRACK',    id_flighttrack)
  IF (laoa)             CALL new_submodel('AOA',              id_aoa      )

  ! report submodel status
  
  IF (p_parallel_io) THEN
     CALL message('', separator)
     IF ( nsubm > 0 ) THEN
       WRITE (message_text,*) nsubm, ' submodels registered:', &
                              (' '//TRIM(submlist(i)%modulename),i=1,nsubm)
     ELSE
       WRITE (message_text,*) 'No submodels registered (=> lanysubmodel = .false.)'
     END IF
     CALL message('', message_text, level=em_param)
     CALL message('', 'Submodel switches processed.', level=em_param)
     CALL print_status('METHOX module', lmethox)
     CALL print_status('TRANSDIAG module', ltransdiag)
     CALL print_status('HAM aerosol module', lham)
     CALL print_status('MOZ chemistry module', lmoz)
     
     IF ( lham .AND. lmoz ) THEN
        CALL print_status('HAMMOZ photolysis coupling', lhmzphoto)
        CALL print_status('HAMMOZ oxidant coupling', lhmzoxi)
        CALL print_status('HAMMOZ het. chemistry coupling', lhmzhet)
     END IF
     
     IF ( lmoz ) THEN
        CALL print_status('Chemistry feedback on radiation', lchemrad)
        CALL print_status('Use of water vapour from chemistry', linterh2o)
        IF (lhammonia) THEN
          CALL print_status('Chemical heating', lchemheat)
          CALL print_status('Use of air mass from chemistry', linteram)
          CALL print_status('Use of specific heat from chemistry', lintercp)
        END IF
        CALL print_status('Lightning NOx emissions', llght)
     END IF
    
     IF ( lbioemi_stdalone ) & 
       CALL print_status('Biogenic NMVOC emissions as standalone submodel', lbioemi_stdalone)
     CALL print_status('Satellite processor diagnostics', losat)
     CALL print_status('ISCCP cloud diagnostics', loisccp)
     CALL print_status('Burden diagnostics', lburden)
     CALL print_status('CCN climatology', lccnclim)
     CALL print_status('Flight-track simulator', lflighttrack)
     
     SELECT CASE ( iadvec )
        CASE (no_advection) 
           CALL message('', 'Run without tracer advection!', level=em_warn)
        CASE (semi_lagrangian) 
           CALL message('', 'Run with semi-lagrangian tracer advection! (are you sure?)', &
                        level=em_warn)
        CASE (tpcore)
           CALL message('', 'Run with Lin&Rood (tpcore) tracer advection.', level=em_param)
     END SELECT
    
!++mgs 
     IF (lham .OR. lmoz .OR. llght .OR. lbioemi_stdalone .OR. lco2) THEN
!--mgs
       IF ( .NOT. lchemistry )     CALL message('', 'Run with lchemistry = false !', level=em_warn)
       IF ( .NOT. lemissions )     CALL message('', 'Run with lemissions = false !', level=em_warn)
       IF ( .NOT. ldrydep )        CALL message('', 'Run with ldrydep = false !', level=em_warn)
       IF ( .NOT. lwetdep )        CALL message('', 'Run with lwetdep = false !', level=em_warn)
       IF ( .NOT. lsedimentation ) CALL message('', 'Run with lsedimentation = false !', level=em_warn)
       IF ( .NOT. laero_micro )    CALL message('', 'Run with laero_micro = false !', level=em_warn)
     END IF
     IF (lmoz .AND. .NOT. llght) THEN
       CALL message('', 'Run with MOZ chemistry but without lightning NO!', level=em_warn)
     END IF
            
     CALL print_status('Age-of-air submodel', laoa)

     CALL message('', separator)
     
  END IF

END SUBROUTINE setsubmodel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Flag first tracer defined in a given submodel.
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Externals
!!
!! @par Notes
!! not yet used
!!  
!! @par Responsible coder
!! M.Schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE starttracdef (idm)
  ! helper routine to manage tracer ids in submodel list
    USE mo_tracdef,        ONLY: ntrac
    INTEGER, INTENT(in)     :: idm

    submlist (idm) % idtfirst = ntrac+1
  END SUBROUTINE starttracdef


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Flag last tracer defined in a given submodel.
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Externals
!!
!! @par Notes
!! not yet used
!!  
!! @par Responsible coder
!! M.Schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE endtracdef (idm)
  ! helper routine to manage tracer ids in submodel list
    USE mo_tracdef,        ONLY: ntrac
    USE mo_exception,      ONLY: message, message_text, em_info, em_error

    INTEGER, INTENT(in)     :: idm

    submlist (idm) % idtlast  = ntrac

    IF (submlist (idm) % idtfirst > submlist (idm) % idtlast ) THEN
       CALL message('init_subm','No tracers defined for submodel '//submlist (idm) % modulename, &
                    level=em_error)
    ELSE
       WRITE (message_text,*) submlist (idm) % idtlast - submlist (idm) % idtfirst + 1, &
                              ' tracers defined for submodel '//submlist (idm) % modulename
       CALL message('init_subm', message_text, level=em_info)
    END IF

  END SUBROUTINE endtracdef



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Print status of a submodel or parameterisation (active or not active).
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! setsubmodel
!!
!! @par Externals
!!
!! @par Notes
!! not yet used
!!  
!! @par Responsible coder
!! M.Schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE print_status (mstring, flag)

  USE mo_exception,      ONLY: message, message_text, em_param

  IMPLICIT NONE

  CHARACTER(len=*), intent(in)   :: mstring
  LOGICAL, intent(in)            :: flag

  IF ( flag ) THEN
     write(message_text,'(a60,1x,":",a)') mstring,'active'
  ELSE
     write(message_text,'(a60,1x,":",a)') mstring,'*not* active'
  END IF
  CALL message('', message_text, level=em_param)

END SUBROUTINE print_status



SUBROUTINE print_lvalue (mstring, lvalue)

!------------------------------------------------------------------------------
! Report the value of a logical, integer or real variable
! Convenience routine interfaced by print_value(mstring, value)

  USE mo_exception,      ONLY: message, message_text, em_param

  IMPLICIT NONE

  CHARACTER(len=*), intent(in)   :: mstring
  LOGICAL, intent(in)            :: lvalue

  IF (lvalue) THEN
    write(message_text,'(a60,1x,": ",a)') mstring,'TRUE'
  ELSE
    write(message_text,'(a60,1x,": ",a)') mstring,'FALSE'
  END IF
  CALL message('', message_text, level=em_param)

END SUBROUTINE print_lvalue



SUBROUTINE print_ivalue (mstring, ivalue)

  USE mo_exception,      ONLY: message, message_text, em_param

  IMPLICIT NONE

  CHARACTER(len=*), intent(in)   :: mstring
  INTEGER, intent(in)            :: ivalue

  write(message_text,'(a60,1x,":",i10)') mstring, ivalue
  CALL message('', message_text, level=em_param)

END SUBROUTINE print_ivalue



SUBROUTINE print_rvalue (mstring, rvalue)

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: message, message_text, em_param

  IMPLICIT NONE

  CHARACTER(len=*), intent(in)   :: mstring
  REAL(wp), intent(in)           :: rvalue

  write(message_text,'(a60,1x,":",g12.5)') mstring, rvalue
  CALL message('', message_text, level=em_param)

END SUBROUTINE print_rvalue



END MODULE mo_submodel
