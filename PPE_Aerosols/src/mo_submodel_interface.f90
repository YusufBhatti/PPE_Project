!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! #define HAMMOZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file is supposed to provide all necessary interfaces of submodels
!! to the main ECHAM6 model. It controls submodel initialisation, including
!! tracer and stream definitions and it contains standard interfaces to allow
!! chemistry and aerosol processes interacting with the physical climate model.
!! The coding objective is to use mo_submodel_interface as "single entry point"
!! into the details of each individual submodel. The core ECHAM code should
!! only make use of flags and switches in mo_submodel and *not* contain any 
!! references to submodel-specific switches or values. 
!!
!! Note: this is the only file which may contain preprocessor definitions for submodels!!
!!
!! Change this file to attach extra submodels to ECHAM.
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
!! <li>M. Schultz   (FZ-Juelich) -  original idea and code structure - (2009-05-xx) 
!! <li>S. Rast      (MPI-Met)    -  original idea and code structure - (2009-06-xx) 
!! <li>K. Zhang     (MPI-Met)    -  restucture and new style, implementation in ECHAM6  - (2009-07-16)
!! <li>K. Zhang     (MPI-Met)    -  doxygen support - (2009-07-24)
!! <li>M. Schultz   (FZ-Juelich) -  Zuerich code structure and cleanup - (2009-08-25) 
!! <li>M. Schultz&S. Rast   (FZ-Juelich, MPI-Met) -  merge to include HAMMOZ revision 0004 - (2010-01-27)
!! <li>M. Schultz   (FZ-Juelich) -  merge for HAMMOZ rev0008 - (2010-04-14)
!!                                  
!! </ol>
!!
!! @par This module is used by
!! physc
!! vdiff
!! cloud
!! and to_be_added
!! 
!! @par Responsible coder
!! sebastian.rast@zmaw.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE mo_submodel_interface

  IMPLICIT NONE

  PRIVATE

! public interfaces                 ! 
                                    ! === initialisation and administration part 
                                    ! 
  PUBLIC :: init_subm               ! submodel registration, general settings, before tracer definition
                                    !
  PUBLIC :: init_subm_memory        ! initialize memory for submodels
  PUBLIC :: free_subm_memory        ! free memory for submodels
                                    ! 
                                    ! === chemistry interfaces 
                                    ! 
  PUBLIC :: stepon_subm             ! read monthly fields
  PUBLIC :: physc_subm_1            ! called from physc before radiation. Prepare stuff.
                                    ! used specifically for SOA gas-phase/particle-phase equilibrium
                                    ! 
  PUBLIC :: radiation_subm_1        ! initialize submodel fields needed for radiation calculation
                                    ! 
  PUBLIC :: radiation_subm_2        ! diagnose radiation after radiation calculation
!
  PUBLIC :: radheat_subm            ! diagnose heating of radiation here
                                    !
  PUBLIC :: physc_subm_2            ! optional starting point for chemistry calculations
                                    ! called after vdiff but before cuflx and cloud.
                                    ! MOZART2 chemistry used to be called from here --
                                    ! now all chemistry calls in physc_interface_4 (needs testing!)
                                    !
  PUBLIC :: scan1_subm
  PUBLIC :: cuflx_subm              ! called from cuflx.: first interface to wetdep
                                    ! 
  PUBLIC :: vdiff_subm              ! called from vdiff - use for chemical flux form boundary conditions
                                    ! such as surface emissions and dry deposition
                                    !
  PUBLIC :: cloud_subm_1            ! called from cloud_cdnc_icnc: preliminary calculations
                                    !
  PUBLIC :: cloud_subm_2            ! called from cloud / cloud_cdnc, etc.: aerosol wet chemistry
                                    ! 
  PUBLIC :: cloud_subm_3            ! called from cloud / cloud_cdnc, etc.: exports for AeroCom diags

  PUBLIC :: physc_subm_3            ! called from physc. 
                                    ! aerosol and chemistry interface:
                                    ! m7 or salsa, MOZART chemistry, sedimentation, etc.
                                    ! WILL INCLUDE DIAGNOSTICS ###
                                    !
  PUBLIC :: physc_subm_4            ! called from physc after the ocean coupling.
                                    ! CO2 tendency limitation and CO2 diagnostics

!++mgs: presently unused: see code in hammoz_rev0003 if needed

  LOGICAL, PARAMETER  :: luse_p1_vars = .TRUE.    ! use updated physical variables in submodel 
                                                  ! routines. This affects pressure, temperature,
                                                  ! humidity and cloud water/ice.
                                                  ! Submodel routines will only receive tendency
                                                  ! variables if they modify them (often only xtte)
!>>>mgs
  INTEGER          :: idt_NO                      ! tracer index for NO (lightning module)
!<<<mgs

  CONTAINS
                                                             
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!
  !!!   subroutine interfaces to the host atmospheric model 
  !!!
  !!!   the following interface routines are arranged in the order 
  !!!   in which they are called
  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Initialization of submodels: 
!! 
!! <ol>
!! <li>register submodels 
!! <li>set namelist variables for submodels
!! <li>set output control (yes/no simplified output)
!! <li>initialize module variables
!! <li>register species
!! <li>register modes
!! </ol>
!! 
!! All generic submodel switches (mo_submodel) were processed and the active 
!! submodels have been registered. Use this routine to read submodel specific 
!! namelists, initialize submodel tracers and perform other initialisations. 
!! Please avoid inflation of subroutine calls from this routine.
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
!! <ol>
!! <li>init_ham
!! <li>init_megan
!! <li>init_moz
!! <li>isccp_initialize
!! </ol>
!!
!! @par Notes
!! mo_xt shall replace old tracer functionality (mo_sub_nml): to be done 
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE init_subm
  
  USE mo_exception,        ONLY: message, em_warn, em_info
  USE mo_species,          ONLY: init_splist, printspec
!>>>mgs / ++mgs 20140313
  USE mo_tracer,           ONLY: get_tracer, new_tracer
  USE mo_tracdef,          ONLY: ON, GAS
  USE mo_kind,             ONLY : dp
!<<<mgs / --mgs 20140313
  USE mo_submodel,         ONLY: starttracdef,     &
                                 endtracdef,       &
                                 lmethox,          &
                                 ltransdiag,       &
                                 lbioemi_stdalone, &
                                 lxt,              &
                                 lham,             &
                                 lmoz,             &
                                 lhammoz,          &
                                 lhmzoxi,          &
                                 llght,            &
                                 lccnclim,         &
                                 lchemistry,       & !mgs #249
                                 ldrydep,          & !gf #244
                                 lemissions,       &
                                 lflighttrack,     &
                                 laoa,             &
                                 laerocom_diag,    &
                                 id_ham,           &
                                 id_moz,           &
                                 id_ccnclim
  USE mo_sub_nml,          ONLY: request_tracer_nml, set_tracer_nml
  USE mo_aoa,              ONLY: init_aoa
  USE mo_exception,        ONLY: message, em_error
#ifdef HAMMOZ
  USE mo_ham_init,         ONLY: start_ham, ham_define_tracer, ham_initialize
  USE mo_moz_init,         ONLY: start_moz, moz_define_tracer, moz_initialize
!>>>mgs
  USE mo_moz_lightning,    ONLY: start_lightning
!<<<mgs
  USE mo_hammoz,           ONLY: hammoz_initialize
  USE mo_ham,              ONLY: print_aerocomp_info
  USE mo_hammoz_aerocom_diags, ONLY: init_aerocom
#endif

  USE mo_hyb,              ONLY: cetah
!!  USE mo_isccpsimulator,   ONLY: isccp_initialize
  USE mo_methox,           ONLY: init_methox
  !! rs USE mo_co2,              ONLY: init_co2
  USE mo_ccnclim,          ONLY: init_ccnclim_submodel, & !SF CCN climatology
                                 ccnclim_define_tracer
  USE mo_activ,            ONLY: activ_initialize         !SF
  USE mo_param_switches,   ONLY: lcdnc_progn              !SF

#ifdef HAMMOZ
  USE mo_hammoz_drydep,    ONLY: drydep_init !gf #244

  USE mo_hammoz_emissions, ONLY: init_emissions
  !>>dod is not, and will not be (for the foreseeable future) tested outside of HAMMOZ. So moved inside ifdef.
  USE mo_hammoz_emi_biogenic, ONLY: start_biogenic_emissions
  !<<dod
  !>>dod split of mo_timer (redmine #51)
  USE mo_hammoz_timer,     ONLY: init_hammoz_timers
  !<<dod
#endif
  USE mo_control,          ONLY: ltimer, &
                                 lforcererun !SF (#141): disable internal rerun in case of HAMMOZ
  USE mo_flighttrack,      ONLY: flighttrack_initialize, flighttrack_init_stream_elements
!>>SF
  USE mo_cloud_utils, ONLY: init_cloud_micro_2m
!<<SF

!!mgs!!   USE mo_xt,               ONLY: idm_xt, setxt, xt_define_tracers, xt_init    ! to be completed


!>>ps 
  USE mo_hammoz_perturbations, ONLY: init_hammoz_perturbations, init_hammoz_emi_perturbations
!<<ps

  INTEGER:: ierr

  ! 0) --- Preparations: initialize species list

  CALL init_splist

  ! 1) --- Start various submodels

#ifdef HAMMOZ
  !-- disable internal reruns (SF, #141)
  IF (lforcererun) THEN
     CALL message('init_subm', "lforcererun set to FALSE, as internal reruns are incompatible with HAMMOZ", &
                 level=em_warn)
     CALL message('','')
     lforcererun = .FALSE.
  ENDIF

  !>>ps
  !--- Initialise hammoz perturbed physics setup - before emissions and before ham:
  IF (lham) CALL init_hammoz_perturbations
  !<<ps
 
  ! -- initialize aerosol module and its species
  IF (lham) CALL start_ham              !ham aerosol module

  ! -- initialize chemistry module and its species
  IF (lmoz) CALL start_moz(lchemistry, lhmzoxi)   !moz chemistry module !mgs #249
#else
  ! -- give warnings in the non-HAMMOZ version of ECHAM if ham or moz are activated
  IF (lham) THEN
    CALL message('init_subm', "LHAM set to false, because this ECHAM version doesn't support HAMMOZ!", &
                 level=em_warn)
    lham = .FALSE.
  END IF
  IF (lmoz) THEN
    CALL message('init_subm', "LMOZ set to false, because this ECHAM version doesn't support HAMMOZ!", &
                 level=em_warn)
    lmoz = .FALSE.
  END IF
#endif

!!rs  init_co2 is now call in initialize (in all cases, also for lco2=false)
  !if (lco2) CALL init_co2(co2mmr)
!!rs


  ! Test prognostic cdnc switch (cannot be done in setsubmodel, because setphys is called afterwards!)
  IF (lcdnc_progn .AND. .NOT. (lham .OR. lccnclim)) THEN
    CALL message('init_subm', 'Cannot run with lcdnc_progn=true and neither lham nor lccnclim active!', &
                 level=em_error)
    lcdnc_progn = .FALSE.
  END IF

!>>SF Initialize boundary conditions necessary for 2-moment cloud microphysics
  IF (lcdnc_progn) THEN
     CALL init_cloud_micro_2m
  ENDIF
!<<SF

  ! -- initialize CCN climatology
  IF (lccnclim) CALL init_ccnclim_submodel !SF CCN climatology module

  IF (lmethox) CALL init_methox                !upper atmospheric H2O production

  !!mgs!!IF (lxt)    CALL setxt                !simple tracer module

  ! 2) --- Define submodel species and tracers
  !!mgs!!! Simple tracer module (former XT_ functionality)
  !!mgs!!IF (lxt) THEN
  !!mgs!!  CALL starttracdef(idm_xt)
  !!mgs!!  CALL xt_define_tracers
  !!mgs!!  CALL endtracdef(idm_xt)
  !!mgs!!END IF

  ! -- CCN climatology  module
  IF (lccnclim) THEN
    CALL starttracdef(id_ccnclim)
    CALL ccnclim_define_tracer
    CALL endtracdef(id_ccnclim)
  END IF

#ifdef HAMMOZ
  ! -- HAM aerosol module
  IF (lham) THEN
    CALL starttracdef(id_ham)
    CALL ham_define_tracer
    CALL endtracdef(id_ham)
  END IF

  ! -- MOZ chemistry module
  IF (lmoz) THEN
    CALL starttracdef(id_moz)
    CALL moz_define_tracer
    CALL endtracdef(id_moz)
  END IF

  ! 3) --- Submodel-specific initialisation routines

  !!mgs!!    !
  !!mgs!!    ! Simple tracer module (former XT_ functionality)
  !!mgs!!    IF (lxt) CALL xt_init

 !--- For AeroCom diags
  IF (laerocom_diag) CALL init_aerocom

  ! -- HAM aerosol module
  IF (lham) THEN
    CALL ham_initialize
  END IF

  ! -- MOZ chemistry module
  IF (lmoz) CALL moz_initialize

!>>gf #244
  !-- initialize dry deposition scheme
  IF (ldrydep) CALL drydep_init
!<<gf

  ! -- initialize biogenic emissions module
  !>>dod is not, and will not be (for the foreseeable future) tested outside of HAMMOZ. So moved inside ifdef.
  IF (lbioemi_stdalone) CALL start_biogenic_emissions
  !<<dod

  ! -- Lightning module (can be run independent of MOZ)
!>>>mgs
  IF (llght) THEN
    CALL get_tracer('NO',idx=idt_NO,ierr=ierr)
!++mgs 20140313
    IF (idt_NO <= 0) CALL new_tracer('NO', 'LGHT', units='kg kg-1', moleweight=30.0_dp, nphase=GAS, nwrite=ON, idx=idt_NO)
!--mgs 20140313
    CALL start_lightning
  END IF
!<<<mgs

  !
  !--- 4) perform cross-submodel initialisations
  !
  IF (lhammoz) CALL hammoz_initialize

  ! --- report on tracer, species and mode definitions
  IF (lham .OR. lmoz .OR. lbioemi_stdalone) THEN
    CALL printspec
    CALL print_tracer_species_mode_info
    IF (lham) CALL print_aerocomp_info
  END IF

 !--- parse emission matrix and prepare use of emissions
  IF (lemissions) CALL init_emissions

#endif
 
  IF (lcdnc_progn) CALL activ_initialize

  IF (laoa) CALL init_aoa !initialization of age of air tracers

  !>>dod
#ifdef HAMMOZ
  !---start submodel timers
  IF (ltimer) CALL init_hammoz_timers
#endif
  !<<dod

  IF (lflighttrack) THEN
     !---Initialize flight-track simulator
     CALL flighttrack_initialize
  ENDIF
 
  !>>ps
  !--- Initialise hammoz perturbed physics setup - after emissions:
#ifdef HAMMOZ
  if (lham) CALL init_hammoz_emi_perturbations
#endif
  !<<ps

  END SUBROUTINE init_subm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Routines are called from here to initialize diagnostiv streams and allocate buffers
!!
!! @author see module info of mo_submodel_interface
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info of mo_submodel_interface
!!
!! @par This subroutine is called by
!! init_memory in mo_memory_streams
!!
!! @par Externals:
!! <ol>
!! </ol>
!!
!! @par Notes
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  SUBROUTINE init_subm_memory

    USE mo_submodel,         ONLY: lco2, lham, lmoz, llght,  &
                                   lemissions, lflighttrack, lwetdep, ldrydep, &
                                   lsedimentation, lanysubmodel, ltransdiag, &
                                   laoa, laerocom_diag
!--mgs
    USE mo_tracdef,          ONLY: trlist

    !   general process streams and diagnostics
    USE mo_submodel_streams, ONLY: init_submodel_streams
    USE mo_vphysc,           ONLY: init_vphysc_stream
#ifdef HAMMOZ
    USE mo_hammoz_wetdep,    ONLY: init_wetdep_stream
    USE mo_hammoz_drydep,    ONLY: init_drydep_stream
    USE mo_hammoz_sedimentation,   ONLY: init_sedi_stream
    USE mo_hammoz_emissions, ONLY: init_emi_stream
    USE mo_hammoz_aerocom_diags, ONLY: init_aerocom_streams
#endif
    USE mo_param_switches,   ONLY: lcdnc_progn 
    USE mo_conv,             ONLY: construct_stream_conv
    USE mo_activ,            ONLY: construct_activ_stream 
!!    USE mo_localtime,      ONLY: mo_localtime_init

    !   specific submodel streams
    USE mo_co2,              ONLY: init_co2_memory, reference_co2
    USE mo_radiation_parameters, ONLY: co2mmr
    USE mo_transdiag,        ONLY: construct_transdiag_stream
#ifdef HAMMOZ
    USE mo_ham_init,         ONLY: ham_init_memory
    USE mo_moz_diag,         ONLY: init_moz_streams
    USE mo_moz_init,         ONLY: init_moz_ic
    USE mo_moz_lightning,    ONLY: init_lightning_stream
    !SF #171 USE mo_hammoz_global_diag, ONLY: init_hammoz_gdiag
!>>NAJS: high frequency output
    USE mo_ham_stream_hifreq,ONLY: construct_stream_hifreq
!<<NAJS 
#endif
!!    USE mo_isccpsimulator,   ONLY: construct_stream_isccp
    USE mo_flighttrack,      ONLY: flighttrack_init_output, flighttrack_init_stream_elements
    USE mo_aoa,              ONLY: get_pointer2trac

    IMPLICIT NONE

    ! Initialize memory for CO2 module
    ! This is done even if interactive CO2 tracer is not used
    CALL init_co2_memory(co2mmr)
    
    IF (.NOT. lanysubmodel) RETURN

    IF (lco2) CALL reference_co2  ! assign a pointer to co2 tracer (only for interactive CO2 transport)

    ! initialization of generic submodel variables and diagnostics:
    ! * vphysc
    ! * wetdep    IF (.NOT. lanysubmodel) RETURN
    ! * drydep
    ! * sedi
    CALL init_submodel_streams      ! namelist read etc.
    CALL init_vphysc_stream

#ifdef HAMMOZ
    !++mgs #249: acknowledge general lwetdep, ldrydep, and lsedimentation flags! (fix from Sabine)
    IF (lwetdep .AND. ANY(trlist%ti(:)%nwetdep > 0)) CALL init_wetdep_stream
    IF (ldrydep .AND. ANY(trlist%ti(:)%ndrydep > 0)) CALL init_drydep_stream
    IF (lsedimentation .AND. ANY(trlist%ti(:)%nsedi   > 0)) CALL init_sedi_stream
    !--mgs
  ! -- initialize the emissions diagnostics
  ! NB: this is done here and not in alloc_mem, because the initialisation of the dust module
  !     makes use of this stream reference already (called from ham_initialize)
  IF (lemissions) CALL init_emi_stream
#endif

  IF (lcdnc_progn) THEN 
    CALL construct_activ_stream
    CALL construct_stream_conv
  ENDIF

#ifdef HAMMOZ
    !! allocate memory for HAM aerosol model
    IF (lham) CALL ham_init_memory

    !! allocate memory for MOZ chemistry module
    IF (lmoz) THEN
      CALL init_moz_streams     ! MOZ diagnostic output
      CALL init_moz_ic          ! MOZ initial conditions pre-set
    END IF
   
    IF (laerocom_diag) CALL init_aerocom_streams
#endif

    IF (ltransdiag) CALL construct_transdiag_stream

#ifdef HAMMOZ
    IF (llght)      CALL init_lightning_stream

!!mgs!!!--- initialisation of generic submodel streams
!!mgs!!
!!mgs!!!IF (llocaltimed) THEN
!!mgs!!!   CALL mo_localtime_init
!!mgs!!!END IF

    !! initialize global diagnostics for HAMMOZ components
    !SF #171 CALL init_hammoz_gdiag

!>>NAJS: high frequency output
    CALL construct_stream_hifreq
!<<NAJS 
#endif

  IF (lflighttrack) THEN
      !---Initialize flight-track simulator
      CALL flighttrack_init_stream_elements
      CALL flighttrack_init_output
  ENDIF
  IF (laoa) CALL get_pointer2trac  ! get pointers to age of air tracers

  END SUBROUTINE init_subm_memory
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Routines are called from here to deallocate memory
!!
!! @author see module info of mo_submodel_interface
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info of mo_submodel_interface
!!
!! @par This subroutine is called by
!! free_memory in mo_memory_streams
!!
!! @par Externals:
!! <ol>
!! </ol>
!!
!! @par Notes
!! each init_memory routine must have free_memory counterpart; init_streams doesn't need this.
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!### may need to add various routines to avoid memory leaks!!!

  SUBROUTINE free_subm_memory

    USE mo_submodel,       ONLY: lflighttrack, lham, llght, lmoz
    USE mo_co2,            ONLY: cleanup_co2_memory
#ifdef HAMMOZ
    USE mo_ham_init,       ONLY: ham_free_memory
#endif
    USE mo_flighttrack,    ONLY: flighttrack_finalize_output
    IMPLICIT NONE
    IF (lflighttrack) CALL flighttrack_finalize_output
    CALL cleanup_co2_memory

#ifdef HAMMOZ
    IF (lham)   CALL ham_free_memory
!    IF (lmoz)   CALL moz_free_memory
#endif

    !!mgs!!IF (llght) CALL lightning_free_memory
    !!mgs!!!--- deallocate memory of generic submodel components
    !!mgs!!! ...
    !!mgs!!! CALL drydep_free_memory
    !!mgs!!! CALL wetdep_free_memory

  END SUBROUTINE free_subm_memory

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Interface for reading submodel boundary conditions
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! stepon before first scans over gaussian latitudes
!!
!! @par Notes
!! Most emissions now make use of the boundary condition scheme in mo_boundary_condition
!! Their reading is handled automatically via the call to bc_list_read in stepon
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE stepon_subm (current_date, next_date)

  USE mo_time_conversion,     ONLY: time_days
  USE mo_time_control,        ONLY: get_date_components, lstart, lresume
!>>SW: added lmoz for feature #415 in HAMMOZ readmine
  USE mo_submodel,            ONLY: lco2, lham, lccnclim, lmoz
!<<SW
#ifdef HAMMOZ
  USE mo_ham_tools,           ONLY: calc_daylength
  USE mo_ham,                 ONLY: idsec_dust, idsec_biogenic
  USE mo_ham_dust,            ONLY: bgc_dust_read_monthly, ndurough
!>>SW: see feature #415 in HAMMOZ readmine
!>>SW: see feature #469 in HAMMOZ readmine
  USE mo_moz_dailyetf,        ONLY: read_dailyetf, set_dailyetf
  USE mo_radiation_parameters,ONLY: isolrad, lyr_perp
!<<SW
#endif
  USE mo_co2,                 ONLY: read_co2_emission
  USE mo_ccnclim,             ONLY: read_ccnclim    !++mgs ### could this use a boundary condition?

  IMPLICIT NONE

  TYPE(time_days)     :: current_date
  TYPE(time_days)     :: next_date

  INTEGER             :: iyear, imonth, iday
  INTEGER, SAVE       :: imonthm1 = -1           ! previous month
#ifdef HAMMOZ
!>>SW see feature #415 in HAMMOZ readmine
  INTEGER             :: inextyear, inextmonth, inextday
  LOGICAL             :: lnewyear, lnewday
!<<SW
#endif

  !--- obtain current date
  CALL get_date_components (current_date, year=iyear, month=imonth, day=iday)

#ifdef HAMMOZ
  !--- calculate daylength for HAM
  IF (lham) CALL calc_daylength
!>>SW see feature #415 and feature #469 in HAMMOZ readmine
  ! read in daily etf-flux for photolysis if spectrally resolved solar
  ! radiation is used (isolrad == 1)
  IF (lmoz .and. isolrad == 1 .and. .not. lyr_perp) THEN

     CALL get_date_components (next_date, month=inextmonth, year=inextyear, day=inextday)

     lnewyear=iyear/=inextyear
     lnewday=iday/=inextday

     IF (lresume .OR. lstart) THEN
        IF (lnewyear) THEN
           CALL read_dailyetf(inextyear)
        ELSE
           CALL read_dailyetf(iyear)
        END IF
     ELSE IF (lnewyear) THEN
        CALL read_dailyetf(inextyear)
     END IF
     IF (lnewday .OR. lresume .OR. lstart) CALL set_dailyetf

   ENDIF
!<<SW
#endif

  !--- read monthly fields (### should be replaced by BC scheme)
  IF (imonth /= imonthm1) THEN
#ifdef HAMMOZ
    !-- BGC dust emissions
!sschr: bgc_dust_read_monthly must be called in order to get fpar_field
!       (even if dust is processed via bc scheme)
    IF (lham .AND. idsec_dust > 0) CALL bgc_dust_read_monthly(imonth, ndurough)

#endif

    IF (lccnclim) CALL read_ccnclim(iyear,imonth)

    !-- JSBACH CO2 emissions
!!mgs!!    IF (lco2) CALL read_co2_emission  
  END IF

  IF (lco2) CALL read_co2_emission       ! all date calculations done internally

  !--- save month
  imonthm1 = imonth

  END SUBROUTINE stepon_subm

#ifdef HAMMOZ
  SUBROUTINE print_tracer_species_mode_info

  USE mo_mpi,         ONLY: p_parallel_io
  USE mo_tracdef,     ONLY: ntrac, trlist
  USE mo_species,     ONLY: nspec, speclist
  USE mo_ham,         ONLY: naerocomp, aerocomp, &
                            subm_ngasspec, subm_gasspec, subm_naerospec, subm_aerospec
  USE mo_exception,   ONLY: message_text, message, em_param, em_info
  USE mo_util_string, ONLY: separator

  INTEGER     :: jt

     CALL message('',separator)
     WRITE(message_text,'(a)') 'Tracer species and size class information'
     CALL message('',message_text,level=em_info)
     CALL message('','',level=em_param)
     WRITE(message_text,'(a)') ':  id  fullname                   spid  class  ntran' 
     CALL message('',message_text,level=em_param)
     DO jt=1,ntrac
        WRITE(message_text,'(a,i4,2x,a,3(1x,i6))') ':',jt, trlist%ti(jt)%fullname, &
                                                   trlist%ti(jt)%spid, trlist%ti(jt)%mode, trlist%ti(jt)%ntran 
        CALL message('',message_text,level=em_param)
     ENDDO

     CALL message('','',level=em_param)
     CALL message('',separator)

  END SUBROUTINE print_tracer_species_mode_info
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! First interface to submodels that is called in physc before WMO tropopause 
!! and radiation calculation.
!! Historically this was call_chem0 in call_submodels. 
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! physc
!!
!! @par Externals:
!! <ol>
!! <li>soa_equi0
!! </ol>
!!
!! @par Notes
!! this subroutine calls routines that must be executed before the radiation 
!! calculation
!! pressure has not been updated yet to t+dt
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE physc_subm_1 (kproma, kbdim, klev, &
                           klevp1, ktrac, krow, &
                           papm1,  paphm1,      &
                           ptm1,   ptte,        &
                           pxtm1,  pxtte,       &  
                           pqm1,   pqte          ) 


  USE mo_kind,          ONLY : wp
  USE mo_control,       ONLY : ltimer
  USE mo_timer,         ONLY : timer_start, timer_stop, timer_physc_subm_1
  USE mo_time_control,  ONLY : time_step_len
  USE mo_submodel,      ONLY : lco2
  USE mo_co2,           ONLY : co2_mbalance
#ifdef HAMMOZ
  USE mo_ham,           ONLY : nsoa
  USE mo_ham_soa_processes, ONLY : soa_equi0
#endif

  IMPLICIT NONE

! arguments
! pxtm1 must be intent(inout) for soa
  
  INTEGER,  INTENT(in)    :: kproma                 ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                  ! geographic block maximum number of locations 
  INTEGER,  INTENT(in)    :: klev                   ! number of levels
  INTEGER,  INTENT(in)    :: klevp1                 ! number of levels + 1 
  INTEGER,  INTENT(in)    :: ktrac                  ! number of tracers
  INTEGER,  INTENT(in)    :: krow                   ! geographic block number

  REAL(wp), INTENT(in)    :: ptm1  (kbdim,klev)     ! air temperature (t-dt)   [K]
  REAL(wp), INTENT(in)    :: ptte  (kbdim,klev)     ! air temperature tendency [K]
  REAL(wp), INTENT(in)    :: pqm1  (kbdim,klev)     ! specific humidity (t-dt) [kg/kg]
  REAL(wp), INTENT(in)    :: pqte  (kbdim,klev)     ! specific humidity tendency (t-dt) [kg/kg]
  REAL(wp), INTENT(in)    :: papm1 (kbdim,klev)     ! air pressure at layer center    (t-dt) [Pa]
  REAL(wp), INTENT(in)    :: paphm1(kbdim,klevp1)   ! air pressure at layer interface (t-dt) [Pa]
  REAL(wp), INTENT(inout) :: pxtm1 (kbdim,klev,ktrac) ! tracer mass/number mixing ratio (t-dt)
  REAL(wp), INTENT(inout) :: pxtte (kbdim,klev,ktrac) ! tracer mass/number mixing ratio tendency

! local variables
  REAL(wp) :: zt(kbdim,klev)
  REAL(wp) :: zq(kbdim,klev)
  REAL(wp) :: zxt(kbdim,klev,ktrac)
  INTEGER  :: jt

  IF (ltimer) CALL timer_start(timer_physc_subm_1)

  !-- update physical variables
  IF (luse_p1_vars) THEN
    zt(1:kproma,:) = ptm1(1:kproma,:) + ptte(1:kproma,:)*time_step_len
    zq(1:kproma,:) = pqm1(1:kproma,:) + pqte(1:kproma,:)*time_step_len
  ELSE
    zt(1:kproma,:) = ptm1(1:kproma,:)
    zq(1:kproma,:) = pqm1(1:kproma,:)
  END IF

  CALL co2_mbalance (kproma,          kbdim,             klev,             &
                     klevp1,          krow,            paphm1)

#ifdef HAMMOZ
  ! set diagnostic tracers for SOA and calculate gas/aerosol equilibrium
  IF (nsoa == 1) THEN 
     CALL soa_equi0(kproma, kbdim, klev, &
                    krow,   ktrac,       &
                    papm1,               &
                    zt,     zq,          &
                    pxtm1,  pxtte     )
  END IF
#endif
  IF (ltimer) CALL timer_stop(timer_physc_subm_1)

  END SUBROUTINE physc_subm_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Radiative heating diagnostics
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! radheat
!!
!! @par Externals:
!! <ol>
!! </ol>
!!
!! @par Notes
!!
!! @par Responsible coder
!! sebastian.rast@zmaw.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE radheat_subm(kproma       ,kbdim           ,klev             ,&
                          klevp1       ,krow            ,pconvfact        ,&
                          pflxs        ,pflxt)

  USE mo_kind,        ONLY: wp

!! This subroutine is for diagnostics of heat rates
  IMPLICIT NONE
! arguments
  
  INTEGER, INTENT(in)  :: kproma     ! actual length of grid point block
  INTEGER, INTENT(in)  :: kbdim      ! maximum length of grod point block
  INTEGER, INTENT(in)  :: klev       ! number of levels
  INTEGER, INTENT(in)  :: klevp1     ! number of levels + 1
  INTEGER, INTENT(in)  :: krow       ! geographic block number
  REAL(wp), INTENT(in) :: pconvfact(kbdim,klev) ! conversion factor
  REAL(wp), INTENT(in) :: pflxs(kbdim,klevp1), pflxt(kbdim,klevp1) ! heat fluxes for solar and thermal wave lengths

  END SUBROUTINE radheat_subm

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! This subroutine is one of the main "portals" to the chemistry submodules. 
!! It is called from physc after vdiff and radheat, before gwspectrum, 
!! cucall, etc.
!! In the past, call_chem1 (equivalent to this interface routine) was used to 
!! drive the MOZART2 chemistry, while HAMMONIA/MOZART3 was called via call_chem2
!! (which is now physc_interface_3). 
!! This should now be controlled by a namelist switch (at least until we know 
!! what works best)
!! Note that HAM was always called from call_chem2 and radionucl_sink always 
!! from call_chem1.
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! physc
!!
!! @par Externals:
!! <ol>
!! <li>to_be_added
!! </ol>
!!
!! @par Notes
!! not yet used
!! try to keep arguments similar to physc_subm_2 in order to maintain flexibility for
!! chemistry calls
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE physc_subm_2                               &
            (kproma, kbdim, klev,  klevp1, ktrac, krow, &
             itrpwmo, itrpwmop1,                        &
             paphm1, papm1, paphp1, papp1,              &
             ptm1,   ptte,  ptsurf,                     &
             pqm1,   pqte,  pxlm1, pxlte, pxim1, pxite, &
             pxtm1,  pxtte,                             &
             paclc,  ppbl,                              &
             loland, loglac                            )   

  
  USE mo_kind,              ONLY: wp
  USE mo_control,           ONLY: ltimer
  USE mo_time_control,      ONLY: time_step_len
  USE mo_timer,             ONLY: timer_start,   &
                                  timer_stop,    &
                                  timer_physc_subm_2
  USE mo_physical_constants,ONLY: vtmpc1
  USE mo_param_switches,    ONLY: ncd_activ
  USE mo_submodel,          ONLY: lmoz, llght, lemissions
  USE mo_sub_echam,         ONLY: radionucl_sink
  USE mo_vphysc,            ONLY: set_vphysc_var
#ifdef HAMMOZ
  USE mo_ham,               ONLY: nsoa
! USE mo_moz,               ONLY: ichem
! USE mo_moz_aircraft,      ONLY: moz_aircraft_emissions
! USE mo_lightning,         ONLY: lightning_emissions
  USE mo_ham_soa_processes, ONLY: soa2prod, soa_progn, soa_equi0
  ! >> thk: volatility basis set (VBS)
  USE mo_ham_vbs_partition, ONLY:  vbs_gas_phase_chem
  ! << thk

#endif


  IMPLICIT NONE

  
! arguments
  
  INTEGER, INTENT(in)  :: kproma                      ! geographic block number of locations
  INTEGER, INTENT(in)  :: kbdim                       ! geographic block maximum number of locations 
  INTEGER, INTENT(in)  :: klev                        ! numer of levels
  INTEGER, INTENT(in)  :: klevp1                      ! numer of levels + 1
  INTEGER, INTENT(in)  :: ktrac                       ! numer of tracers
  INTEGER,  INTENT(in) :: krow                        ! geographic block number
  INTEGER,  INTENT(in) :: itrpwmo  (kbdim)            ! tropopause on model levels
  INTEGER,  INTENT(in) :: itrpwmop1(kbdim)            ! tropopause on model levels + 1 

  REAL(wp), INTENT(in) :: paphm1   (kbdim,klev+1)     ! air pressure at layer interface (t-dt) [Pa]
  REAL(wp), INTENT(in) :: papm1    (kbdim,klev)       ! air pressure at layer center    (t-dt) [Pa]
  REAL(wp), INTENT(in) :: paphp1   (kbdim,klev+1)     ! air pressure at layer interface (t+dt) [Pa]
  REAL(wp), INTENT(in) :: papp1    (kbdim,klev)       ! air pressure at layer center    (t+dt) [Pa]
  REAL(wp), INTENT(in) :: ptm1     (kbdim,klev)       ! air temperature (t-dt)   [K]
  REAL(wp), INTENT(in) :: ptte     (kbdim,klev)       ! air temperature tendency [K/s]
  REAL(wp), INTENT(in) :: ptsurf   (kbdim)            ! surface temperature [K]
  REAL(wp), INTENT(in) :: pqm1     (kbdim,klev)       ! specific humidity             [kg/kg]
  REAL(wp), INTENT(in) :: pqte     (kbdim,klev)       ! specific humidity tendency    [kg/kg/s]
  REAL(wp), INTENT(in) :: pxlm1    (kbdim,klev)       ! cloud liquid content (t-dt)   [kg/kg]
  REAL(wp), INTENT(in) :: pxlte    (kbdim,klev)       ! cloud liquid content tendency [kg/kg/s]
  REAL(wp), INTENT(in) :: pxim1    (kbdim,klev)       ! cloud ice    content (t-dt)   [kg/kg]
  REAL(wp), INTENT(in) :: pxite    (kbdim,klev)       ! cloud ice    content tendency [kg/kg/s]
  REAL(wp), INTENT(in) :: paclc    (kbdim,klev)       ! cloud fraction [0~1] 
  REAL(wp), INTENT(in) :: ppbl     (kbdim)            ! level of PBL upper boundary
  REAL(wp), INTENT(inout) :: pxtm1 (kbdim,klev,ktrac) ! tracer mass/number mixing ratio (t-dt)   [kg/kg]/[#/kg]
  REAL(wp), INTENT(inout) :: pxtte (kbdim,klev,ktrac) ! tracer mass/number mixing ratio tendency [kg/kg]/[#/kg]
  LOGICAL,  INTENT(in) :: loland   (kbdim)            ! land mask 
  LOGICAL,  INTENT(in) :: loglac   (kbdim)            ! glacier mask 
  
! local variables
  REAL(wp) :: zaph(kbdim, klevp1), zap(kbdim, klev)
  REAL(wp) :: zt(kbdim,klev)
  REAL(wp) :: ztv(kbdim,klev)         ! virtual temperature
  REAL(wp) :: zq(kbdim,klev)
  REAL(wp) :: zxl(kbdim,klev), zxi(kbdim,klev)
  REAL(wp) :: zxt(kbdim,klev,ktrac)
  INTEGER  :: jt

  IF (ltimer) call timer_start(timer_physc_subm_2)

  !-- update physical variables
  IF (luse_p1_vars) THEN
    zaph(1:kproma,:) = paphp1(1:kproma,:)
    zap(1:kproma,:) = papp1(1:kproma,:)
    zt(1:kproma,:) = ptm1(1:kproma,:) + ptte(1:kproma,:)*time_step_len
    zq(1:kproma,:) = pqm1(1:kproma,:) + pqte(1:kproma,:)*time_step_len
    zxl(1:kproma,:) = pxlm1(1:kproma,:) + pxlte(1:kproma,:)*time_step_len
    zxi(1:kproma,:) = pxim1(1:kproma,:) + pxite(1:kproma,:)*time_step_len
  ELSE
    zaph(1:kproma,:) = paphm1(1:kproma,:)
    zap(1:kproma,:) = papm1(1:kproma,:)
    zt(1:kproma,:) = ptm1(1:kproma,:)
    zq(1:kproma,:) = pqm1(1:kproma,:)
    zxl(1:kproma,:) = pxlm1(1:kproma,:)
    zxi(1:kproma,:) = pxim1(1:kproma,:)
  END IF

  !-- compute virtual temperature (pure diagnostics at present)  
  ztv(1:kproma,:) = zt(1:kproma,:)*(1._wp+vtmpc1*zq(1:kproma,:)-(zxl(1:kproma,:)+zxi(1:kproma,:)))

  !-- save physical variables for diagnostic output
  CALL set_vphysc_var &
                      (kproma,                   klev,                &
                       krow,                     paphm1=zaph,         &
                       papm1=zap,                ptvm1=ztv,           &
                       ppbl=ppbl,                ktrpwmo=itrpwmo,     &
                       pqm1=pqm1,                                     &
                       ktrpwmop1=itrpwmop1)

#ifdef HAMMOZ 
  !--- secondary organic aerosol processes
  SELECT CASE (nsoa)
     CASE (1) !D. O'Donnell's SOA implementation
        CALL soa2prod(kproma, kbdim, klev, klevp1, krow,  ktrac, &
                      zaph,   zap,   zt,   pxtm1, pxtte)
        CALL soa_progn(kproma, kbdim, klev, krow,  ktrac, pxtm1, pxtte)
        CALL soa_equi0(kproma, kbdim, klev, krow,  ktrac, &
                       zap,  zt,    zq,   pxtm1, pxtte)

     CASE (2)
        !>>thk #512: VBS scheme
        !gas phase oxidation
        CALL vbs_gas_phase_chem(&
             kproma, kbdim, klev, klevp1, krow,  ktrac,  zaph,   zap,   zt,   pxtm1, pxtte &
                               )
         !<<thk #512
  END SELECT
#endif

  IF (ktrac > 0) THEN
    !
    ! commonly used processes calculated by ECHAM:
    !
    CALL radionucl_sink (kproma, kbdim, klev, pxtm1, pxtte)
    
  ENDIF

  IF (ltimer) call timer_stop(timer_physc_subm_2)
  
  END SUBROUTINE physc_subm_2

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! This is the interface routine for submodel surface flux calculations
!! (deposition and emissions)
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! vdiff
!!
!! @par Externals:
!! <ol>
!! <li>call_chem_bcond
!! <li>xt_drydep
!! </ol>
!!
!! @par Notes
!! 
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!! Grazia Frontoso, C2SM      (2012-02-01)
!! Usage of the input variables defined over land, water, ice to account
!! for the non-linearity in the drydep calculations for gridboxes 
!! containing both water and sea ice.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE vdiff_subm      (kproma, kbdim,  klev,   klevp1,        &
                              ktrac,  krow,                          &
                              ptm1,   pum1,   pvm1,   pqm1,          &
                              papm1,  paphm1, paphp1, pgeom1, ptslm1,&
                              pxtm1,  pseaice,pforest,               &
                              pfrl,   pfrw,   pfri,   pcvs,   pcvw,  &
                              pvgrat, ptsw,   ptsi,                  &
                              pu10,   pv10,                          &
                              paz0,   paz0l,  paz0w,  paz0i,         &
#ifdef HAMMOZ
!SF gf #78
                              pcfml,  pcfmw,  pcfmi,                 &
                              pcfncl, pcfncw, pcfnci, pepdu2, pkap,  &
                              pril,   priw,   prii,   ptvir1,        & 
                              ptvl, ptvw, ptvi,                      &
                              psrfl,  pcdnl,  pcdnw,  pcdni,         &
                              pqss,   pvlt,                          &
#else 
                              pcfm,   pcfnc,  pepdu2, pkap,          &       !!mgs!! pcfm !
                              pri,    ptvir1, ptvl,                  &
! TS                              psrfl,  pcdnl,  pqss,   pvlt,          &
                              psrfl,  pcdn,  pqss,   pvlt,           &
#endif
                              loland,                                &
                              pxtte,  pxtems,                        &
                              pxlm1,  pxim1                          &
#ifdef HAMMOZ
!SF gf #161
                             ,ihpbl                                  & 
#endif
                             )
  USE mo_kind,               ONLY: wp
  USE mo_control,            ONLY: ltimer
  USE mo_timer,              ONLY: timer_start,       &
                                   timer_stop,        &
                                   timer_vdiff_subm 
#ifdef HAMMOZ
  USE mo_hammoz_timer,       ONLY: timer_hammoz_emissions,   &
                                   timer_hammoz_drydep
#endif
  USE mo_physical_constants, ONLY: vtmpc1, rd
  USE mo_submodel,           ONLY: lham,       &
                                   losat,      &
                                   lmoz,       &
                                   lco2,       &
                                   ltransdiag, &
                                   lemissions, &
                                   ldrydep,    &
                                   lccnclim        !++mgs ### ?? why here??
  USE mo_tracdef,            ONLY: trlist
  USE mo_mpi,                ONLY: p_parallel_io
  USE mo_time_control,       ONLY: lstart
  USE mo_co2,                ONLY: co2_exchange
  USE mo_transdiag,          ONLY: transdiag
!  USE mo_meganv2,            ONLY: biogenic_emissions
#ifdef HAMMOZ
  USE mo_hammoz_drydep,      ONLY: drydep_interface
  USE mo_hammoz_emissions,   ONLY: emi_interface
#endif
  USE mo_ccnclim,            ONLY: ccn_3d     !++mgs ### ?? why here??
  USE mo_tracdef,            ONLY: ntrac

!+++sschr: for use of updated tracers (luse_p1_vars)
  USE mo_time_control,ONLY: time_step_len
!---sschr


! arguments

  INTEGER,  INTENT(in) :: kproma                      ! geographic block number of locations
  INTEGER,  INTENT(in) :: kbdim                       ! geographic block maximum number of locations
  INTEGER,  INTENT(in) :: klev                        ! number of levels
  INTEGER,  INTENT(in) :: klevp1                      ! number of levels + 1
  INTEGER,  INTENT(in) :: ktrac                       ! number of tracers
  INTEGER,  INTENT(in) :: krow                        ! geographic block number

  LOGICAL,  INTENT(in) :: loland   (kbdim)            ! logical land mask (including glaciers)
!>>SF gf #78
#ifdef HAMMOZ
  REAL(wp), INTENT(in) :: pcfml    (kbdim)            ! stability dependend momentum transfer coeff. over land
  REAL(wp), INTENT(in) :: pcfmw    (kbdim)            ! stability dependend momentum transfer coeff. water
  REAL(wp), INTENT(in) :: pcfmi    (kbdim)            ! stability dependend momentum transfer coeff. ice

  REAL(wp), INTENT(in) :: pcfncl   (kbdim)            ! function of heat transfer coeff. over land
  REAL(wp), INTENT(in) :: pcfncw   (kbdim)            ! function of heat transfer coeff. over water
  REAL(wp), INTENT(in) :: pcfnci   (kbdim)            ! function of heat transfer coeff. over ice
  REAL(wp), INTENT(in) :: pril     (kbdim)            ! moist richardson number over land
  REAL(wp), INTENT(in) :: priw     (kbdim)            ! moist richardson number over water
  REAL(wp), INTENT(in) :: prii     (kbdim)            ! moist richardson number over ice
  REAL(wp), INTENT(in) :: ptvl     (kbdim)            ! virtual potential temp. over land
  REAL(wp), INTENT(in) :: ptvw     (kbdim)            ! virtual potential temp. over water
  REAL(wp), INTENT(in) :: ptvi     (kbdim)            ! virtual potential temp. over ice
  REAL(wp), INTENT(in) :: pcdnl     (kbdim)           ! see mo_surface_land
  REAL(wp), INTENT(in) :: pcdnw     (kbdim)           ! see mo_surface_ocean
  REAL(wp), INTENT(in) :: pcdni     (kbdim)           ! see mo_surface_ice
#else
  REAL(wp), INTENT(in) :: pcfm     (kbdim,klev)       ! stability dependend momentum transfer coeff.
  REAL(wp), INTENT(in) :: pcfnc    (kbdim)            ! function of heat transfer coeff.
  REAL(wp), INTENT(in) :: pri      (kbdim)            ! moist richardson number
  REAL(wp), INTENT(in) :: ptvl     (kbdim)            ! see mo_surface_land
  REAL(wp), INTENT(in) :: pcdn     (kbdim)            ! see mo_surface_land
#endif
!<<SF gf
  REAL(wp), INTENT(in) :: pepdu2                      ! constant
  REAL(wp), INTENT(in) :: pkap                        ! constant
  REAL(wp), INTENT(in) :: pum1     (kbdim,klev)       ! u-wind (t-dt)
  REAL(wp), INTENT(in) :: pvm1     (kbdim,klev)       ! v-wind (t-dt)
  REAL(wp), INTENT(in) :: ptm1     (kbdim,klev)       ! temperature (t-dt)
  REAL(wp), INTENT(in) :: pqm1     (kbdim,klev)       ! specific humidity (t-dt)
  REAL(wp), INTENT(in) :: papm1    (kbdim,klev)       ! air pressure at layer center (t-dt)
  REAL(wp), INTENT(in) :: paphp1   (kbdim,klev+1)     ! air pressure at layer interface (t+dt)
  REAL(wp), INTENT(in) :: paphm1   (kbdim,klev+1)     ! air pressure at layer interface (t-dt)  (### unus
  REAL(wp), INTENT(in) :: pgeom1   (kbdim,klev)       ! geopotential (t-dt)
  REAL(wp), INTENT(in) :: ptvir1   (kbdim,klev)       ! see vdiff
  REAL(wp), INTENT(in) :: paz0     (kbdim)            ! roughness length
  REAL(wp), INTENT(in) :: ptslm1   (kbdim)            ! surface temperature (t-dt)
  REAL(wp), INTENT(in) :: pfrl     (kbdim)            ! land fraction
  REAL(wp), INTENT(in) :: pfrw     (kbdim)            ! water fraction
  REAL(wp), INTENT(in) :: pfri     (kbdim)            ! ice fraction
  REAL(wp), INTENT(in) :: pcvs     (kbdim)            ! snow cover fraction
  REAL(wp), INTENT(in) :: pcvw     (kbdim)            ! wet skin fraction
  REAL(wp), INTENT(in) :: pvgrat   (kbdim)            ! vegetation ratio
  REAL(wp), INTENT(in) :: psrfl    (kbdim)            ! surface solar flux
  REAL(wp), INTENT(in) :: pu10     (kbdim)            ! 10m w-wind
  REAL(wp), INTENT(in) :: pv10     (kbdim)            ! 10m v-wind
  REAL(wp), INTENT(in) :: pforest  (kbdim)            ! forest fraction
  REAL(wp), INTENT(in) :: pseaice  (kbdim)            ! sea ice fraction
  REAL(wp), INTENT(in) :: ptsi     (kbdim)            ! surface temperature over ice
  REAL(wp), INTENT(in) :: paz0l    (kbdim)            ! roughness length over land
  REAL(wp), INTENT(in) :: paz0w    (kbdim)            ! roughness length over water
  REAL(wp), INTENT(in) :: paz0i    (kbdim)            ! roughness length over ice
  REAL(wp), INTENT(in) :: ptsw     (kbdim)            ! surface temperature over water
  REAL(wp), INTENT(in) :: pvlt     (kbdim)            ! see vdiff    (### unused ###)
  REAL(wp), INTENT(in) :: pqss     (kbdim,klev)       ! see vdiff

  REAL(wp), INTENT(in) :: pxlm1    (kbdim,klev)       ! cloud liquid content (t-dt)   [kg/kg]
  REAL(wp), INTENT(in) :: pxim1    (kbdim,klev)       ! cloud ice    content (t-dt)   [kg/kg]

!>>SF gf #161
#ifdef HAMMOZ
  INTEGER,  INTENT(in) :: ihpbl(kbdim)                ! level of PBL top
#endif
!<<SF gf #161

  REAL(wp), INTENT(inout) :: pxtems   (kbdim,ktrac)      ! surface emissions
  REAL(wp), INTENT(inout) :: pxtm1    (kbdim,klev,ktrac) ! tracer mass/number mixing ratio (t-dt)
!NB: pxtm1 should definitively be INTENT(in) !!!
  REAL(wp), INTENT(inout) :: pxtte    (kbdim,klev,ktrac) ! tracer mass/number mixing ratio tendency

! local variables

  REAL(wp) :: zdensair (kbdim)            ! air density in surface layer
  REAL(wp) :: zxt(kbdim, klev, ktrac)
  INTEGER  :: jt

!++mgs
  REAL(wp) :: zqsfc(kbdim), zqsssfc(kbdim), ztsfc(kbdim)     ! lowest model level values
                                                             ! for humidity, saturation and temperature

  IF (ltimer) CALL timer_start(timer_vdiff_subm)

  ! preparations
  zqsfc(:)   = 0._wp
  zqsssfc(:) = 0._wp
  ztsfc(:)   = 0._wp
  ztsfc(1:kproma)   = ptm1(1:kproma,klev)
  zqsfc(1:kproma)   = pqm1(1:kproma,klev)
  zqsssfc(1:kproma) = pqss(1:kproma,klev)
  zdensair(1:kproma)= papm1(1:kproma,klev)/(ptm1(1:kproma,klev)*rd*(1._wp+vtmpc1*pqm1(1:kproma,klev)))
  IF (luse_p1_vars) THEN
    DO jt = 1,ntrac
      zxt(1:kproma,:,jt)  = pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt) * time_step_len
    END DO
!   ztsfc(1:kproma) = ptm1(1:kproma,klev) + ptte(1:kproma,klev)*time_step_len
!   zqsfc(1:kproma) = pqm1(1:kproma,klev) + pqte(1:kproma,klev)*time_step_len
!?? zqsssfc(1:kproma) = pqss(1:kproma,klev)  
!   zdensair(1:kproma)= papp1(1:kproma,klev)/(ztsfc(1:kproma,klev)*rd*(1._dp+vtmpc1*zqsfc(1:kproma,klev)))
  ELSE
    DO jt = 1,ntrac
      zxt(1:kproma,:,jt)  = pxtm1(1:kproma,:,jt)
    END DO
!   ztsfc(1:kproma)   = ptm1(1:kproma,klev)
!   zqsfc(1:kproma)   = pqm1(1:kproma,klev)
!   zqsssfc(1:kproma) = pqss(1:kproma,klev)
!   zdensair(1:kproma)= papm1(1:kproma,klev)/(ptm1(1:kproma,klev)*rd*(1._dp+vtmpc1*pqm1(1:kproma,klev)))
  ENDIF

  !-------------------------------------
  ! submodel interface: emissions
  !-------------------------------------

  IF (lemissions) THEN
#ifdef HAMMOZ    
    IF (ltimer) CALL timer_start(timer_hammoz_emissions)
    !--- general emissions controlled by emi_matrix
    CALL emi_interface(kproma, kbdim, ktrac, klev, klevp1, krow, paphp1, pgeom1, loland,  &
!gf #161               ptm1, pxtems, pxtte)
                       ptm1, ihpbl, pxtems, pxtte)

    IF (ltimer) CALL timer_stop(timer_hammoz_emissions)
#endif

  END IF  !   lemissions

  !--- CO2 module - this puts all CO2 fluxes from ocean and land into atmosphere
  IF (lco2) CALL co2_exchange(kproma, kbdim, ktrac, pfrl, pxtems, krow)

  !-------------------------------------
  ! submodel interface: dry depostion
  !-------------------------------------

#ifdef HAMMOZ
  IF (ldrydep .AND. trlist%anydrydep>0) THEN

    IF (ltimer) CALL timer_start(timer_hammoz_drydep)
    CALL drydep_interface(kbdim, kproma,  klev,     krow,                    &
!>>gf change in argument list (see #78)
                          zqsfc,  zqsssfc,ztsfc,    pcfml,    pcfmw,         &
                          pcfmi,  pcfncl,   pcfncw, pcfnci,                  &
                          pepdu2, pkap,  pum1,     pvm1,   pgeom1,           &
                          pril,   priw,  prii,                               &
                          ptvir1, ptvl, ptvw, ptvi, paz0,                    &
                          ptslm1,   loland,                                  &
!<<gf
                          pfrl,   pfrw,  pfri,     pcvs,   pcvw,     pvgrat, &
                          psrfl,                           pu10,     pv10,   &
                          pxtems, zxt, zdensair, paphp1, pforest,  ptsi,     &
!>>gf change in argument list (see #78)
                          paz0l, paz0w, paz0i, pcdnl, pcdnw, pcdni           )
!<<gf
    IF (ltimer) CALL timer_stop(timer_hammoz_drydep)

  END IF
#endif

  IF (lccnclim) CALL ccn_3d(kproma, kbdim, klev, krow, papm1)

  IF (ltransdiag) CALL transdiag(kproma,   kbdim,   klev,   klevp1,   krow,  &
                                 pum1,     pvm1,    ptm1,   pqm1,            &
                                 pxlm1,    pxim1,                            &
                                 paphm1,   pgeom1                            )

  IF (ltimer) CALL timer_stop(timer_vdiff_subm)

  END SUBROUTINE vdiff_subm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! interface between cloud physics and aerosol calculations
!!
!! @author Sylvaine Ferrachat, ETHZ
!!
!! $Id: 1423$
!!
!! @par Revision History
!! <ol>
!! <li>
!! </ol>
!!
!! @par This subroutine is called by
!! cloud_cdnc_icnc
!!
!! @par Externals:
!! <ol>
!! </ol>
!!
!! @par Notes
!! This routine provides an interface to decouple active aerosol module calculations (HAM, etc...) from
!! echam. The non-active aerosols case uses a CCN climatology to feed into cloud_cdnc_icnc
!!
!!
!! @par Responsible coder
!! sylvaine.ferrachat@env.ethz.ch
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE cloud_subm_1(                                               &
             kproma,     kbdim,      klev,       krow, ktdia,            &
             ptkem1,   pwcape, pvervel, prho,                            &
             prho_rcp, pxtm1, pxtte, ptm1, papm1, pqm1, pesw,            &
             pcdncact,                                                   &
             prwetki, prwetai, prwetci,                                  &
             pfracdusol, pfracduai, pfracduci, pfracbcsol, pfracbcinsol, &
             pascs, papnx, paprx, papsigx, ld_het                        )

   USE mo_kind,           ONLY: wp
   USE mo_exception,      ONLY: em_error
   USE mo_submodel,       ONLY: lham, lccnclim
   USE mo_tracdef,        ONLY: ntrac
   USE mo_param_switches, ONLY: ncd_activ
   USE mo_activ,          ONLY: activ_updraft, activ_lin_leaitch, nfrzmod, nw
   USE mo_ccnclim,        ONLY: ccnclim_avail_activ_lin_leaitch, ccnclim_IN_setup
#ifdef HAMMOZ
   USE mo_conv,           ONLY: cdncact_cv
   USE mo_ham,            ONLY: nclass, nccndiag, sizeclass !>>dod<<
   USE mo_ham_activ,      ONLY: ham_avail_activ_lin_leaitch,           &
                                ham_activ_diag_lin_leaitch,            &
                                ham_activ_koehler_ab,                  &
                                ham_activ_abdulrazzak_ghan,            &
                                ham_activ_diag_abdulrazzak_ghan_strat, &
                                ham_activ_diag_abdulrazzak_ghan_conv
   USE mo_ham_freezing,   ONLY: ham_IN_setup
   USE mo_ham_streams,    ONLY: a, b, sc, rdry
   ! --> thk: A-R&G for SALSA
   USE mo_ham_salsa_cloud, ONLY: salsa_abdul_razzak_ghan
   USE mo_ham,             ONLY: nham_subm, HAM_M7, HAM_SALSA
   ! <-- thk
#endif
   !>>dod activation timer
   USE mo_control,        ONLY: ltimer
   USE mo_timer,          ONLY: timer_start, timer_stop
   USE mo_hammoz_timer,   ONLY: timer_ham_activation
   !<<dod
   !>>DN
   USE mo_hammoz_aerocom_diags, ONLY: lHEaci, lHEaci_activ, lMMPPE
   USE mo_hammoz_aerocom_HEaci, ONLY: ccn_inst, ccn1_inst, ccn3_inst
   USE mo_activ,                ONLY: swat
   !<<DN

   INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow, ktdia
   REAL(wp), INTENT(IN) :: ptkem1(kbdim,klev)
   REAL(wp), INTENT(IN) :: pwcape(kbdim)
   REAL(wp), INTENT(IN) :: pvervel(kbdim,klev)
   REAL(wp), INTENT(IN) :: prho(kbdim,klev)
   REAL(wp), INTENT(IN) :: prho_rcp(kbdim,klev)
   REAL(wp), INTENT(IN) :: pxtm1(kbdim,klev,ntrac)
   REAL(wp), INTENT(IN) :: pxtte(kbdim,klev,ntrac)
   REAL(wp), INTENT(IN) :: ptm1(kbdim,klev)
   REAL(wp), INTENT(IN) :: papm1(kbdim,klev)
   REAL(wp), INTENT(IN) :: pqm1(kbdim,klev)
   REAL(wp), INTENT(IN) :: pesw(kbdim,klev)

   !-- activation:
   REAL(wp), INTENT(OUT) :: pcdncact(kbdim,klev)
   !-- mixed-phase freezing (despite their names, these quantities ARE NOT dependent on HAM):
   REAL(wp), INTENT(out) :: prwetki(kbdim,klev)  ! wet radius, aitken insoluble mode
   REAL(wp), INTENT(out) :: prwetai(kbdim,klev)  ! wet radius, accumulation insoluble mode
   REAL(wp), INTENT(out) :: prwetci(kbdim,klev)  ! wet radius, coarse insoluble mode
   REAL(wp), INTENT(out) :: pfracdusol(kbdim,klev)   ! total fraction of dust particules in all soluble mo
   REAL(wp), INTENT(out) :: pfracduai(kbdim,klev)    ! fraction of dust particules in the accum. soluble m
   REAL(wp), INTENT(out) :: pfracduci(kbdim,klev)    ! fraction of dust particules in the coarse soluble m
   REAL(wp), INTENT(out) :: pfracbcsol(kbdim,klev)   ! total fraction of BC particules in all soluble mode
   REAL(wp), INTENT(out) :: pfracbcinsol(kbdim,klev) ! total fraction of BC particules in all insoluble mo
   !-- cirrus freezing:
   REAL(wp), INTENT(out) :: pascs(kbdim,klev)           ! soluble aerosol number conc.
   REAL(wp), INTENT(out) :: papnx(kbdim,klev,nfrzmod)   ! aerosol number available for freezing [1/cm3]
   REAL(wp), INTENT(out) :: paprx(kbdim,klev,nfrzmod)   ! radius of aerosols avail. for freezing  [cm]
   REAL(wp), INTENT(out) :: papsigx(kbdim,klev,nfrzmod) ! std. dev. of aerosols available for freezing
   LOGICAL, INTENT(out)  :: ld_het  !switch to set heterogeneous freezing below 235K (cirrus scheme)

   REAL(wp), ALLOCATABLE :: zw(:,:,:),   & ! mean or bins of updraft velocity [m s-1]
                            zwpdf(:,:,:)   ! updraft velocity pdf over bins
#ifdef HAMMOZ
   REAL(wp) :: znact(kbdim,klev,nclass), & ! number of activated particles per mode [m-3]
               zfracn(kbdim,klev,nclass),& ! fraction of activated particles per mode
               zsc(kbdim,klev,nclass),   & ! critical supersaturation [% 0-1]
               za(kbdim,klev,nclass),    & ! curvature parameter A of the Koehler equation
               zb(kbdim,klev,nclass),    & ! hygroscopicity parameter B of the Koehler equation
               zrdry(kbdim,klev,nclass)    ! dry radius for each classe 
                                           !@@@ currently re-stored, check if avoidable!
   REAL(wp), ALLOCATABLE :: zrc(:,:,:,:), & ! critical radius of activation per mode and w bin [m]
                            zsmax(:,:,:)    ! maximum supersaturation per w bin [% 0-1]
   INTEGER :: jclass
   !>>DN
   REAL(wp) :: swattemp(kbdim,klev)
   !<<DN
#endif

   !--- Cloud droplets (activation)

     ALLOCATE(zw(kbdim,klev,nw))
     ALLOCATE(zwpdf(kbdim,klev,nw))
     !--- Updraft velocity calculation:
     CALL activ_updraft(kproma,   kbdim,  klev,    krow, &
                        ptkem1,   pwcape, pvervel, prho, &
                        zw, zwpdf                        )

#ifdef HAMMOZ
     !--- Koehler coefficients: used for Abdul-Razzak & Ghan activation
     IF (nham_subm == HAM_M7 .AND. (ncd_activ == 2 .OR. nccndiag > 0)) THEN
         CALL ham_activ_koehler_ab(kproma, kbdim, klev, krow, ktdia, &
                                   pxtm1,  ptm1,  za,   zb    )

         !--- Write the values into the A and B output streams (which is
         !--- how they get to ham_ccn. See also #586)
         DO jclass=1,nclass
           IF (sizeclass(jclass)%lactivation) THEN
             a(jclass)%ptr(1:kproma,:,krow) = za(1:kproma,:,jclass)
             b(jclass)%ptr(1:kproma,:,krow) = zb(1:kproma,:,jclass)
           END IF
         END DO
     END IF
#endif

     !--- Activation:
!!++mgs: double-check logic for lcdnc_progn and ncd_activ !! ncd_activ=0 ??  ###
     IF (ncd_activ == 1) THEN ! Lin & Leaitch scheme

       !--- Computes the available number of particules for activation:
       IF (lccnclim) THEN ! CCN climatology
         CALL ccnclim_avail_activ_lin_leaitch(kproma, kbdim, klev, krow)
#ifdef HAMMOZ
       ELSEIF (lham) THEN ! HAM aerosols
         CALL ham_avail_activ_lin_leaitch(kproma, kbdim, klev, krow, prho, pxtm1)
#endif
       ENDIF

       !--- Computes the real number of activated particules
       !--- Lin & Leaitch does not support updraft PDF (nw > 1); this
       !--- is enforced in setphys.
       CALL activ_lin_leaitch(kproma, kbdim, klev, krow,       &
                              zw(:,:,1), pcdncact)

#ifdef HAMMOZ
       !--- Do some HAM-specific diagnostics
       IF (lham) CALL ham_activ_diag_lin_leaitch(kproma, kbdim, klev, krow, prho, pxtm1, pcdncact)
#endif

     ELSE IF (ncd_activ == 2) THEN !  Abdul-Razzak and Ghan scheme
                                   ! (only possible with active aerosols, security check done in setphys)

#ifdef HAMMOZ
        IF (ltimer) CALL timer_start(timer_ham_activation)
         ALLOCATE(zrc(kbdim,klev,nclass,nw))
         ALLOCATE(zsmax(kbdim,klev,nw))

         SELECT CASE(nham_subm)

            CASE(HAM_M7)

               DO jclass=1, nclass
                  zrdry(1:kproma,:,jclass)=rdry(jclass)%ptr(1:kproma,:,krow)
               END DO

               !>>DN
               IF(lHEaci.AND.lHEaci_activ)THEN
                  swattemp(1:kproma,:)  = swat(1:kproma,:,krow)
                  swat(1:kproma,:,krow) = 0.001_wp
                  CALL ham_activ_abdulrazzak_ghan(kproma, kbdim, klev, krow, ktdia, &
                                               ccn1_inst(:,:,krow), pesw, prho,             &
                                               pxtm1, ptm1, papm1, pqm1,         &
                                               zw, zwpdf, za, zb, zrdry,         &
                                               znact, zfracn, zsc, zrc, zsmax)
                  swat(1:kproma,:,krow) = 0.003_wp
                  CALL ham_activ_abdulrazzak_ghan(kproma, kbdim, klev, krow, ktdia, &
                                               ccn3_inst(:,:,krow), pesw, prho,             &
                                               pxtm1, ptm1, papm1, pqm1,         &
                                               zw, zwpdf, za, zb, zrdry,         &
                                               znact, zfracn, zsc, zrc, zsmax)
                  swat(1:kproma,:,krow)  = swattemp(1:kproma,:)
                  lHEaci_activ = .FALSE.
               END IF
               !<<DN
               
               CALL ham_activ_abdulrazzak_ghan(kproma, kbdim, klev, krow, ktdia, &
                                               pcdncact, pesw, prho,             &
                                               pxtm1, ptm1, papm1, pqm1,         &
                                               zw, zwpdf, za, zb, zrdry,         &
                                               znact, zfracn, zsc, zrc, zsmax)
               !>>DN
               IF(lHEaci) ccn_inst(1:kproma,:,krow) = pcdncact(1:kproma,:)             
               !<<DN
   
               CALL ham_activ_diag_abdulrazzak_ghan_strat(kproma, kbdim, klev,   &
                                                          krow, znact, zfracn,   &
                                                          zrc, zsmax)
   
               ! Convective activation uses the stratiform values. This is not
               ! really correct, but keeps the results identical to those before
               ! factoring out the koehler_ab and diag routines.
               cdncact_cv(1:kproma,:,krow) = pcdncact(1:kproma,:)
   
               CALL ham_activ_diag_abdulrazzak_ghan_conv(kproma, kbdim, klev,   &
                                                         krow, znact, zrc, zsmax)
   
            CASE(HAM_SALSA)
               !>> thk #511: AR&G scheme for SALSA

               ! for now we decided to not use diagnostics routines
               ! in order to cut down on output
               CALL salsa_abdul_razzak_ghan(&
                    kproma,   kbdim, klev,  krow, ktdia, &
                    pcdncact, pesw,  prho,               &
                    pxtm1,    ptm1,  papm1, pqm1,        &
                    zw,       zwpdf,                     &
                    zsc,      zsmax                      &
                    )
   
               ! like for M7:
               ! Convective activation uses the stratiform values. This is not
               ! really correct, but keeps the results identical to those before
               ! factoring out the koehler_ab and diag routines.
               cdncact_cv(1:kproma,:,krow) = pcdncact(1:kproma,:)
               !<< thk #511

         END SELECT

         DO jclass=1,nclass
            !>>dod #377
            IF (sizeclass(jclass)%lactivation) THEN
               sc(jclass)%ptr(1:kproma,:,krow) = zsc(1:kproma,:,jclass)
            END IF
            !<<dod
         END DO

         DEALLOCATE(zsmax)
         DEALLOCATE(zrc)

         IF (ltimer) CALL timer_stop(timer_ham_activation)
         
#endif

     ENDIF

     DEALLOCATE(zwpdf)
     DEALLOCATE(zw)

   !--- Ice crystals (freezing: mixed-phase and cirrus)

     IF (lccnclim) THEN ! CCN climatology

        CALL ccnclim_IN_setup(kproma, kbdim, klev, krow,                                  &
                              prho, prho_rcp,                                             &
                              prwetki, prwetai, prwetci,                                  &
                              pfracdusol, pfracduai, pfracduci, pfracbcsol, pfracbcinsol, &
                              pascs, papnx, paprx, papsigx, ld_het                        )

     ELSE IF (lham) THEN ! HAM aerosols
#ifdef HAMMOZ
       CALL ham_IN_setup(kproma, kbdim, klev, krow,                                  &
                         prho, pxtm1, pxtte, pcdncact,                               &
                         prwetki, prwetai, prwetci,                                  &
                         pfracdusol, pfracduai, pfracduci, pfracbcsol, pfracbcinsol, &
                         pascs, papnx, paprx, papsigx, ld_het                        )
#endif
     ENDIF

  END SUBROUTINE cloud_subm_1

         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! interface between cloud physics and sulfur chemistry
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! <ol>
!! <li> S. Ferrachat, ETH Zurich. Replacement of old xt_chemistry call by new ham_wet_chemistry
!!                                and wetdep_interface calls (2009.08; 2009.12)
!! </ol>
!!
!! @par This subroutine is called by
!! cloud
!!
!! @par Externals:
!! <ol>
!! </ol>
!!
!! @par Notes
!! The wet deposition interface is also called from cuflx_subm.
!! 
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  SUBROUTINE cloud_subm_2(                                 &
             kproma,     kbdim,      klev,       ktop,     &
             krow,                                         &
             pmlwc,      pmiwc,      pmratepr,   pmrateps, &
             pfrain,     pfsnow,     pfevapr,    pfsubls,  &
             pmsnowacl,  paclc,      ptm1,       ptte,     &
             pxtm1,      pxtte,      paphp1,     papp1,    &
             prhop1,     pclcpre)

  USE mo_kind,              ONLY: wp
  USE mo_physical_constants,ONLY: grav
  USE mo_control,           ONLY: ltimer
  USE mo_timer,             ONLY: timer_start,          &
                                  timer_stop,           &
                                  timer_cloud_subm 
  !>>dod split mo_timer (redmine #51) 
  USE mo_hammoz_timer,      ONLY: timer_ham_wetchem,    &
                                  timer_hammoz_wetdep
  !<<dod
  USE mo_time_control,      ONLY: time_step_len
  USE mo_submodel,          ONLY: lchemistry,           &
                                  lham,                 &
                                  lhammoz,              &
                                  lwetdep
  USE mo_tracdef,           ONLY: ntrac, trlist
  USE mo_tracer_processes,  ONLY: xt_borrow   
#ifdef HAMMOZ
  USE mo_hammoz,            ONLY: hammoz_set_oxidants
  USE mo_ham_chemistry,     ONLY: ham_wet_chemistry
  USE mo_hammoz_wetdep,     ONLY: wetdep_interface
#endif

! arguments

  INTEGER,  INTENT(in)    :: kproma                      ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                       ! maximum number of locations in block
  INTEGER,  INTENT(in)    :: klev                        ! number of levels
  INTEGER,  INTENT(in)    :: ktop                        ! top level index
  INTEGER,  INTENT(in)    :: krow                        ! geographic block number
  REAL(wp), INTENT(in)    :: pclcpre  (kbdim,klev)       ! fraction of grid box covered by precip
  REAL(wp), INTENT(in)    :: pfrain   (kbdim,klev)       ! rain flux before evaporation [kg/m2/s]
  REAL(wp), INTENT(in)    :: pfsnow   (kbdim,klev)       ! snow flux before evaporation [kg/m2/s]
  REAL(wp), INTENT(in)    :: pfevapr  (kbdim,klev)       ! evaporation of rain [kg/m2/s]
  REAL(wp), INTENT(in)    :: pfsubls  (kbdim,klev)       ! sublimation of snow [kg/m2/s]
  REAL(wp), INTENT(in)    :: pmsnowacl(kbdim,klev)       ! accretion rate of snow with cloud droplets in
  REAL(wp), INTENT(in)    :: ptm1     (kbdim,klev)       ! air temperature (t-dt)   [K]
  REAL(wp), INTENT(in)    :: ptte     (kbdim,klev)       ! air temperature tendency [K]
  REAL(wp), INTENT(in)    :: prhop1   (kbdim,klev)       ! air density (t-dt)
  REAL(wp), INTENT(in)    :: papp1    (kbdim,klev)       ! air pressure at layer center    (t+dt) [Pa]
  REAL(wp), INTENT(in)    :: paphp1   (kbdim,klev+1)     ! air pressure at layer interface (t+dt) [Pa]
  REAL(wp), INTENT(inout) :: paclc    (kbdim,klev)       ! cloud fraction [0~1]
  REAL(wp), INTENT(inout) :: pmlwc    (kbdim,klev)       ! cloud liquid content before rain [kg/kg]
  REAL(wp), INTENT(inout) :: pmiwc    (kbdim,klev)       ! cloud ice    content before rain [kg/kg]
  REAL(wp), INTENT(inout) :: pmratepr (kbdim,klev)       ! rain formation rate in cloudy part
  REAL(wp), INTENT(inout) :: pmrateps (kbdim,klev)       ! ice  formation rate in cloudy part
  REAL(wp), INTENT(in)    :: pxtm1    (kbdim,klev,ntrac) ! tracer mass/number mixing ratio (t-dt)   [kg/k
  REAL(wp), INTENT(inout) :: pxtte    (kbdim,klev,ntrac) ! tracer mass/number mixing ratio tendency [kg/k

! local variables

  REAL(wp) :: zdum2d (kbdim,klev)
  REAL(wp) :: zdum3d (kbdim,klev,ntrac)
  REAL(wp) :: zxtp1    (kbdim,klev,ntrac) ! updated tracer mass/number mixing ratio
  REAL(wp) :: zxtp10   (kbdim,klev,ntrac) ! ambient tracer mass/number mixing ratio
  REAL(wp) :: zxtp1c   (kbdim,klev,ntrac) ! in-cloud tracer mass/number mixing ratio
  REAL(wp) :: zlfrac_so2   (kbdim,klev)   ! liquid tracer fraction (SO2) -ham specific-
  REAL(wp) :: zdpg(kbdim,klev)
  REAL(wp) :: zdummy(kbdim,ntrac)         ! placeholder for pxtbound, which is only necessary in the conv. 
                                          ! case

  LOGICAL, PARAMETER :: lstrat = .TRUE.   !SF this switch is used to keep track in wetdep_interface
                                          !   whether the call comes from the stratiform routine or
                                          !   the convective one - that's why it's hardcoded here
  INTEGER  :: jk, jl, jt

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!baustelle!! Check following condition:
  !!if MOZART also has het. SO2 chemistry, we run into trouble here
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  IF (ltimer) CALL timer_start(timer_cloud_subm)

  zlfrac_so2(1:kproma,:) = 0._wp                   ! liquid fraction of SO2 for HAM
  zdpg(:,:) = 0._wp

  !-- (re-) compute geopotential height
  DO jk=1,klev
    DO jl=1,kproma
      zdpg(jl,jk)=(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
    END DO
  END DO

  !--- Mass conserving correction of negative tracer values (!csld: necessary for part of the tracers, see #330):
  CALL xt_borrow(kproma, kbdim,  klev,  klev+1, ntrac, &
                 papp1,  paphp1,                       &
                 pxtm1,  pxtte                         )

  !-- initialise in-cloud and interstitial mixing ratios
  !   set both equal to tracer mixing ratio as starting point
  !   ham_wet_chemistry will re-compute these values if lham=true
  DO jt = 1,ntrac
    zxtp1(1:kproma,:,jt)  = pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt) * time_step_len
    zxtp1c(1:kproma,:,jt) = zxtp1(1:kproma,:,jt)
    zxtp10(1:kproma,:,jt) = zxtp1(1:kproma,:,jt)
  END DO

#ifdef HAMMOZ
  !-- heterogeneous (wet) chemistry calculation
  IF ( lham .AND. lchemistry ) THEN  !SF aero wet chemistry

       IF (lhammoz) CALL hammoz_set_oxidants(kproma, kbdim, klev, krow, ntrac, pxtm1, pxtte)

       IF (ltimer) CALL timer_start(timer_ham_wetchem)
       CALL ham_wet_chemistry(kproma,     kbdim,      klev,                     &
                              ktop,       krow,       pmlwc,      pmiwc,        &
                              paclc,      ptm1,       ptte,                     &
                              pxtm1,      pxtte,      zxtp1,                    &
                              zxtp10,     zxtp1c,                               &
                              paphp1,     papp1,      prhop1,     zlfrac_so2,   &
                              zdpg  )
       IF (ltimer) CALL timer_stop(timer_ham_wetchem)

  END IF

  !-- interface to wet deposition routine (also from cuflx_subm)
  IF ( lwetdep .AND. ANY(trlist%ti(:)%nwetdep > 0) ) THEN

      IF (ltimer) CALL timer_start(timer_hammoz_wetdep)
      CALL wetdep_interface(kproma, kbdim,    klev, ktop, krow,      lstrat, &
                            zdpg,   pmratepr, pmrateps,   pmsnowacl,         &
                            pmlwc,  pmiwc,                                   &
                            ptm1, pxtm1, zlfrac_so2, pxtte, zxtp10, zxtp1c,  &
                            pfrain, pfsnow, pfevapr, pfsubls,                &
                            zdum2d, zdum3d,                                  &
                            paclc,  pclcpre, prhop1, zdummy)
      IF (ltimer) CALL timer_stop(timer_hammoz_wetdep)

  END IF
#endif

  IF (ltimer) CALL timer_stop(timer_cloud_subm)

  END SUBROUTINE cloud_subm_2

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! interface between cloud physics and (so far) AerChemMIP
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! <ol>
!! <li> S. Ferrachat, ETH Zurich. Replacement of old xt_chemistry call by new ham_wet_chemistry
!!                                and wetdep_interface calls (2009.08; 2009.12)
!! </ol>
!!
!! @par This subroutine is called by
!! cloud
!!
!! @par Externals:
!! <ol>
!! </ol>
!!
!! @par Notes
!! The wet deposition interface is also called from cuflx_subm.
!! 
!!
!! @par Responsible coder
!! sylvaine.ferrachat@env.ethz.ch
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  SUBROUTINE cloud_subm_3(                                 &
             kproma,     kbdim,      klev,       krow,     &
             kcld_top, ld_liqcld,                          &
             pclcov, pcdnc, prho, paclc, prleff, prsfl, pssfl,&
             paut, pacc, pdpg, pqccol, pfrain, pfevapr, paclci)             

    USE mo_kind, ONLY: wp
#ifdef HAMMOZ
    !>>SF 
    USE mo_submodel,               ONLY: laerocom_diag
    USE mo_hammoz_aerocom_diags,   ONLY: lHEaci,lMMPPE
    USE mo_hammoz_aerocom_HEaci,   ONLY: nstr_cld_top, nstr_cld_liq_top, cdnc3d,&
         f3d, clt_inst, cdr3d, sprecip_inst, accret3d, auto3d, dpg, riming3d, &
         qcsedten3d, qcevap3d
    USE mo_hammoz_aerocom_mmppe,   ONLY: srain_inst
    USE mo_cloud_utils,    ONLY: clc_min
    !<<SF
#endif

    INTEGER, INTENT(in) :: kproma     !< geographic block number of locations
    INTEGER, INTENT(in) :: kbdim      !< maximum number of locations in block
    INTEGER, INTENT(in) :: klev       !< number of levels
    INTEGER, INTENT(in) :: krow       !< geographic block number
    INTEGER, INTENT(in) :: kcld_top(kbdim,klev)   !< level index of stratiform cloud tops
    LOGICAL, INTENT(in) :: ld_liqcld(kbdim,klev)  !< flag for presence of droplets in cloud
    REAL(wp), INTENT(in) :: pclcov(kbdim) !< total cloud cover [fraction]
    REAL(wp), INTENT(in) :: pcdnc(kbdim,klev) !< grid-mean CDNC [m-3]
    REAL(wp), INTENT(in) :: prho(kbdim,klev) !< Air density [kg/m3]
    REAL(wp), INTENT(in) :: paclc(kbdim,klev) !< cloud cover  (now diagnosed in cover)
    REAL(wp), INTENT(in) :: prleff(kbdim,klev) !< liquid effective radius [um]
    REAL(wp), INTENT(in) :: prsfl(kbdim) !< surface rain flux (kg/m2/s)
    REAL(wp), INTENT(in) :: pssfl(kbdim) !< surface snow flux (kg/m2/s)
    REAL(wp), INTENT(in) :: paut(kbdim,klev) !< Rain formation rate for mass (autoconversion) (kg/kg/s)
    REAL(wp), INTENT(in) :: pacc(kbdim,klev) !< Rain formation rate for mass (accretion) (kg/kg/s)
    REAL(wp), INTENT(in) :: pdpg(kbdim,klev)     !< delta p over g [kg/m2]
    REAL(wp), INTENT(in) :: pqccol(kbdim,klev) !< accretion (riming) of droplets by snow - mass
    REAL(wp), INTENT(in) :: pfrain   (kbdim,klev)       !< rain flux before evaporation [kg/m2/s]
    REAL(wp), INTENT(in) :: pfevapr  (kbdim,klev)       !< evaporation of rain [kg/m2/s]
    REAL(wp), INTENT(in) :: paclci(kbdim,klev)      !< cloud cover used for ice clouds [0,1]

    INTEGER :: icld_liq_top(kbdim, klev)
  
#ifdef HAMMOZ
    !>>SF
    IF (laerocom_diag) THEN
       IF (lHEaci) THEN
          !-- Set level of highest stratiform cloud top
          nstr_cld_top(1:kproma,krow) = &
               MINLOC(kcld_top(1:kproma,:),2, &
               (kcld_top(1:kproma,:) > 0 ))

          !-- Set lev index of all liq cloud tops
          icld_liq_top(1:kproma,:) = MERGE( &
               kcld_top(1:kproma,:), 0, &
               ld_liqcld(1:kproma,:)) 
          
          !-- Set level of highest stratiform liq cloud top
          nstr_cld_liq_top(1:kproma,krow) = MINLOC( &
               icld_liq_top(1:kproma,:),2, &
               (icld_liq_top(1:kproma,:) > 0))

          !-- Set cdnc
           cdnc3d(1:kproma,:,krow) = pcdnc(1:kproma,:)*prho(1:kproma,:)
          
          !-- Set cloud cover
          clt_inst(1:kproma,krow) = pclcov(1:kproma)
          f3d(1:kproma,:,krow) = paclc(1:kproma,:)
          cdr3d(1:kproma,:,krow) = MERGE(prleff(1:kproma,:)*0.000001_wp,0._wp,&
               (paclc(1:kproma,:)>clc_min.AND.ld_liqcld(1:kproma,:)))

          !-- Set sprecip
          sprecip_inst(1:kproma,krow) = prsfl(1:kproma)+pssfl(1:kproma)

          !-- Set accr, auto, riming
          accret3d(1:kproma,:,krow) = pacc(1:kproma,:)
          auto3d(1:kproma,:,krow  ) = paut(1:kproma,:)
          riming3d(1:kproma,:,krow) = pqccol(1:kproma,:)

          !-- Set dpg
          dpg(1:kproma,:,krow) = pdpg(1:kproma,:)

          !-- Set qcsedten, qcevap
          qcsedten3d(1:kproma,:,krow) = pfrain(1:kproma,:)*paclci(1:kproma,:)/pdpg(1:kproma,:)
          qcevap3d(1:kproma,:,krow)   = pfevapr(1:kproma,:)*paclci(1:kproma,:)/pdpg(1:kproma,:)
          
       ENDIF
       IF (lMMPPE) THEN
          !-- Set srain
          srain_inst(1:kproma,krow) = prsfl(1:kproma)
       END IF
    ENDIF
    !<<SF 
#endif
  END SUBROUTINE cloud_subm_3

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! This subroutine is the "portal" to the chemistry submodules. 
!! It is called from physc after cloud and before surf.
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! cloud
!!
!! @par Externals:
!! <ol>
!! <li>ham_subm_interface
!! <li>xt_aero_prop
!! <li>xt_sedimentation
!! <li>timer_start
!! <li>timer_stop
!! </ol>
!!
!! @par Notes
!! 
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!>>>mgs: replaced loland and loglac with pfrl, pfrw and pfri <<<
 
  SUBROUTINE physc_subm_3                               &
            (kproma, kbdim, klev,  klevp1, ktrac, krow, &
             paphm1, papm1, paphp1, papp1,              &     ! pressure variables
             ptm1,   ptte,  ptsurf,                     &     ! temperature variables
             pqm1,   pqte,                              &     ! gas-pase water (humidity)
             pxlm1,  pxlte, pxim1, pxite,               &     ! cloud water and ice
             pxtm1,  pxtte,                             &     ! tracer mass mixing ratios
             pgeom1, pgeohm1,                           &     ! geopotential height
             paclc,                                     &     ! cloud cover
             ppbl,   pvervel,                           &     ! PBL height and vertical velocity
             pfrl,   pfrw,  pfri,                       &
             pfrg                                        )    ! land, water, ice and glacier fractions



  USE mo_kind,            ONLY: wp
  USE mo_time_control,    ONLY: time_step_len
  USE mo_control,         ONLY: ltimer
  USE mo_physical_constants, ONLY: vtmpc1
  USE mo_timer,           ONLY: timer_start,          &
                                timer_stop,           &
                                timer_physc_subm_3 
  !>>dod split mo_timer (redmine #51)                                
  USE mo_hammoz_timer,    ONLY: timer_ham_bulk,       & !sf
                                timer_ham_m7,         &
                                timer_ham_salsa,      & !sf
                                timer_hammoz_sedimentation,  &
                                timer_ham_gaschem,    &
                                timer_moz_chem,       &
                                timer_hammoz_burden
  USE mo_tracdef,         ONLY: ntrac, trlist
  USE mo_submodel,        ONLY: lmethox,              &
                                lham,                 &
                                lmoz,                 &
                                llght,                &
                                lhmzoxi,              &
                                lchemistry,           &
                                lsedimentation,       &
                                laero_micro,          &
                                lburden,              &
                                losat,                &
                                loisccp
  USE mo_vphysc,          ONLY: vphysc
  USE mo_methox,          ONLY: methox
!!  USE mo_isccpsimulator,  ONLY: call_isccpsimulator
  USE mo_tracer_processes,   ONLY: xt_burden
#ifdef HAMMOZ
  USE mo_hammoz,          ONLY: hammoz_set_aerosol_bc
  USE mo_hammoz_sedimentation, ONLY: sedi_interface
  USE mo_ham,             ONLY: nsoa,                 &
                                nham_subm,            &
                                nccndiag,             &
                                HAM_BULK,             &
                                HAM_M7,               &
                                HAM_SALSA
  USE mo_ham_subm,        ONLY: ham_subm_interface
  USE mo_ham_chemistry,   ONLY: ham_gas_chemistry
  USE mo_moz_chemistry,   ONLY: moz_chemistry        !++mgs rev0008
!!  USE mo_moz,             ONLY: ichem
!!  USE mo_moz_aircraft,    ONLY: moz_aircraft_emissions
!!  USE mo_lightning,       ONLY: lightning_emissions
  USE mo_ham_soa_processes,   ONLY: soa_progn, soa_equi0
! USE mo_moz_test_zen,      ONLY: test_zen      !++mgs : testing of zenith angle calculation MOZ vs ECHAM
! USE mo_moz,               ONLY: setmoz        !++mgs needed for test_zen - d2r
  USE mo_ham_ccn,         ONLY: ham_ccn
#endif


! arguments
! pxtm1 must be intent(inout) for soa

  INTEGER,  INTENT(in)    :: kproma                      ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                       ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: klev                        ! number of levels
  INTEGER,  INTENT(in)    :: klevp1                      ! number of levels + 1
  INTEGER,  INTENT(in)    :: ktrac                       ! number of tracers
  INTEGER,  INTENT(in)    :: krow                        ! geographic block number
  REAL(wp), INTENT(in)    :: paphm1   (kbdim,klevp1)     ! air pressure at layer interface (t-dt) [Pa]
  REAL(wp), INTENT(in)    :: papm1    (kbdim,klev)       ! air pressure at layer center    (t-dt) [Pa]
  REAL(wp), INTENT(in)    :: paphp1   (kbdim,klevp1)     ! air pressure at layer interface (t+dt) [Pa]
  REAL(wp), INTENT(in)    :: papp1    (kbdim,klev)       ! air pressure at layer center    (t+dt) [Pa]
  REAL(wp), INTENT(in)    :: ptm1     (kbdim,klev)       ! air temperature (t-dt)   [K]
  REAL(wp), INTENT(inout) :: ptte     (kbdim,klev)       ! air temperature tendency [K]
  REAL(wp), INTENT(in)    :: ptsurf   (kbdim)            ! surface temperature [K]
  REAL(wp), INTENT(in)    :: pqm1     (kbdim,klev)       ! specific humidity             [kg/kg]
  REAL(wp), INTENT(inout) :: pqte     (kbdim,klev)       ! specific humidity tendency    [kg/kg]
  REAL(wp), INTENT(in)    :: pxlm1    (kbdim,klev)       ! cloud liquid content (t-dt)   [kg/kg]
  REAL(wp), INTENT(inout) :: pxlte    (kbdim,klev)       ! cloud liquid content tendency [kg/kg]
  REAL(wp), INTENT(in)    :: pxim1    (kbdim,klev)       ! cloud ice    content (t-dt)   [kg/kg]
  REAL(wp), INTENT(inout) :: pxite    (kbdim,klev)       ! cloud ice    content tendency [kg/kg]
  REAL(wp), INTENT(inout) :: pxtm1    (kbdim,klev,ktrac) ! tracer mass/number mixing ratio (t-dt)   [kg/kg]/[#/kg]
  REAL(wp), INTENT(inout) :: pxtte    (kbdim,klev,ktrac) ! tracer mass/number mixing ratio tendency [kg/kg]/[#/kg]
  REAL(wp), INTENT(in)    :: pgeom1   (kbdim,klev)       ! geopotential height at full levels (t-dt) [m]
  REAL(wp), INTENT(in)    :: pgeohm1  (kbdim,klevp1)     ! geopotential height at interfaces (t-dt) [m]
  REAL(wp), INTENT(in)    :: paclc    (kbdim,klev)       ! cloud fraction [0~1]
  REAL(wp), INTENT(in)    :: ppbl     (kbdim)            ! PBL top level
  REAL(wp), INTENT(in)    :: pvervel  (kbdim,klev)       ! lage scale vertical velocity [Pa s-1]
  REAL(wp), INTENT(in)    :: pfrl     (kbdim)            ! land fraction
  REAL(wp), INTENT(in)    :: pfrw     (kbdim)            ! water fraction
  REAL(wp), INTENT(in)    :: pfri     (kbdim)            ! (sea) ice fraction
  REAL(wp), INTENT(in)    :: pfrg     (kbdim)            ! (land) ice fraction

! local variables
  REAL(wp) :: zaph(kbdim, klevp1), zap(kbdim, klev)
  REAL(wp) :: zt(kbdim, klev)
  REAL(wp) :: zq(kbdim, klev)
  REAL(wp) :: zxl(kbdim, klev), zxi(kbdim, klev)
  REAL(wp) :: zxt(kbdim, klev, ktrac)
#ifdef HAMMOZ
  REAL(wp) :: ztvm1(kbdim,klev)         ! virtual temperature at t-dt
#endif

  REAL(wp) :: zgrvolm1 (kbdim,klev)        ! grid box volume

  INTEGER  :: jt
  INTEGER  :: itimer !SF dummy var to be used for timing any aerosol microphysics scheme

  IF (ltimer) CALL timer_start(timer_physc_subm_3)

  !--- Copy or update physical variables
  IF (luse_p1_vars) THEN
    zaph(1:kproma,:) = paphp1(1:kproma,:)
    zap(1:kproma,:)  = papp1(1:kproma,:)
    zt(1:kproma,:)   = ptm1(1:kproma,:) + ptte(1:kproma,:)*time_step_len
    zq(1:kproma,:)   = pqm1(1:kproma,:) + pqte(1:kproma,:)*time_step_len
    zxl(1:kproma,:)  = pxlm1(1:kproma,:) + pxlte(1:kproma,:)*time_step_len
    zxi(1:kproma,:)  = pxim1(1:kproma,:) + pxite(1:kproma,:)*time_step_len
  ELSE
    zaph(1:kproma,:) = paphm1(1:kproma,:)
    zap(1:kproma,:)  = papm1(1:kproma,:)
    zt(1:kproma,:)   = ptm1(1:kproma,:)
    zq(1:kproma,:)   = pqm1(1:kproma,:)
    zxl(1:kproma,:)  = pxlm1(1:kproma,:)
    zxi(1:kproma,:)  = pxim1(1:kproma,:)
  END IF

  !--- compute virtual temperature (not needed (?))
!!  ztv(1:nproma,:) = zt(1:nproma,:)*(1._wp+vtmpc1*zq(1:nproma,:)-(zxl(1:nproma,:)+zxi(1:nproma,:)))

  !--- Upper atmospheric H2O production from methane oxidation and destruction from HO photolysis
  IF (lmethox) THEN
    CALL methox(kproma, kbdim, klev, zq, pqte, papm1)
  ENDIF

#ifdef HAMMOZ
  !--- Chemistry solvers  !gf: now done before the aerosol microphysics

  ! ... interface to HAM gas-phase chemistry
  IF (lchemistry .AND. lham .AND. .NOT. lhmzoxi ) THEN

    IF (ltimer) CALL timer_start(timer_ham_gaschem)

    CALL ham_gas_chemistry(kbdim, klev, kproma, krow, &
                           zaph,  zap,  zt,     zq,   &
                           pxtm1, pxtte               )

    IF (ltimer) CALL timer_stop(timer_ham_gaschem)

  END IF

!!mgs!!     ! ... interface to MOZART chemistry including independent submodules
!!mgs!!     IF (lemissions .AND. llght .and. ichem >= 10) THEN
!!mgs!!        CALL lightning_emissions(kproma, kbdim, klev, ktrac, pxtte, krow)
!!mgs!!     END IF

  IF (lchemistry .AND. lmoz) THEN

    IF (ltimer) CALL timer_start(timer_moz_chem)

    IF (luse_p1_vars) THEN
      DO jt = 1, ktrac
        zxt(1:kproma,:,jt) = pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt)*time_step_len
      END DO
    ELSE
      DO jt = 1, ktrac
        zxt(1:kproma,:,jt) = pxtm1(1:kproma,:,jt)
      END DO
    END IF

    CALL hammoz_set_aerosol_bc(kproma, kbdim, klev, krow, ntrac, zxt)

!++mgs 20130304: divide ptsurf by time_step_len
    CALL moz_chemistry (kproma, kbdim, klev, klevp1, ktrac, krow,  &
                        zaph,   zap,   zt,   ptte,                 &
                        zq,     pqte,  zxl,  pxlte,  zxi,   pxite, &
                        zxt,    pxtte,                             &
                        pgeom1, pgeohm1,     ptsurf/time_step_len, paclc,        & 
                        pfrl,   pfrw,  pfri, pfrg                    )
!--mgs
    IF (ltimer) CALL timer_stop(timer_moz_chem)

  END IF

  !--- aerosol microphysics

  IF (lham .AND. laero_micro) THEN

    zgrvolm1(:,:) = vphysc%grvolm1(:,:,krow)

    SELECT CASE(nham_subm)
          
        CASE(HAM_BULK)
            itimer = timer_ham_bulk
        CASE(HAM_M7)
            itimer = timer_ham_m7
        CASE(HAM_SALSA)
            itimer = timer_ham_salsa
        END SELECT

    IF(ltimer) CALL timer_start(itimer)
    
       CALL ham_subm_interface(kproma,  kbdim, klev, krow, ktrac, & 
                               zap,     zaph,  zt,   zq,          & 
                               pxtm1,   pxtte,                    & 
                               paclc,   zgrvolm1,    ppbl         )
       
    IF (ltimer) CALL timer_stop(itimer)

  END IF

  !--- aerosol sedimentation

  IF (lsedimentation .AND. trlist%anysedi>0) THEN

    IF (ltimer) CALL timer_start(timer_hammoz_sedimentation)

    IF (ANY(trlist%ti(:)%nsedi > 0)) THEN

      CALL sedi_interface(kbdim, kproma, klev, krow, &
                          zt,    zq,     zap,  zaph, &
                          pxtm1, pxtte               )
    END IF

    IF (ltimer) CALL timer_stop(timer_hammoz_sedimentation)

  END IF


  !--- extended aerosol diagnostics:

!!mgs!!  IF (lham .AND. laerocom_diag) CALL aero_diag(kproma, kbdim, klev, krow, &
!!mgs!!                                               pxtm1,  pxtte, &
!!mgs!!                                               pxlm1,  pxlte, &
!!mgs!!                                               pxim1,  pxite, &
!!mgs!!                                               ptm1,   ptte,  &
!!mgs!!                                               papp1,  paphp1  )


!    !--- Extended MOZ diagnostics
!    IF (lmoz) THEN
!    ! save total tracer tendency
!       IF (llocaltimed) CALL buffer_localtime(krow)
!       CALL moz_diagn(krow)
!    END IF
! soa_progn
! soa_equi0: set diagnostic tracers for SOA and calculate gas/aerosol equilibrium

  IF (nsoa == 1) THEN
      CALL soa_progn(kproma, kbdim, klev, krow,  ktrac, pxtm1, pxtte)
      CALL soa_equi0(kproma, kbdim, klev, krow,  ktrac, &
                     zap,    zt,    zq,   pxtm1, pxtte)
  END IF

  ! -- calculate CCN diagnostics (see also #586)
  ! --   -- should this use t+dt values like other calls in here? -ZK
  IF (nccndiag > 0) THEN

     !--- compute virtual temperature at t-dt
     ztvm1(1:kproma,:) = ptm1(1:kproma,:)*(1._wp+vtmpc1*pqm1(1:kproma,:)-(pxlm1(1:kproma,:)+pxim1(1:kproma,:)))

     CALL ham_ccn(kproma, kbdim, klev,  klevp1, krow, ktrac, &
                  pxtm1,  ztvm1, papm1, paphm1               )
  END IF
#endif

  ! -- calculate tracer burdens
  IF (lburden) THEN
    IF (ltimer) CALL timer_start(timer_hammoz_burden)
    CALL xt_burden(kproma, kbdim, klev, klevp1, krow,  &
                   zap,    zaph,                       &
                   pxtm1,  pxtte                      )
    IF (ltimer) CALL timer_stop(timer_hammoz_burden)
  END IF

  IF (ltimer) CALL timer_stop(timer_physc_subm_3)

  END SUBROUTINE physc_subm_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! interface after hydrology_collect in physc. Needed for CO2 diagnostics and coupling
!!
!! @author see module info of mo_submodel_interface
!!
  SUBROUTINE physc_subm_4 (kproma,                kbdim,              klev,            &
                           klevp1,                ktrac,              krow,            &
                           paphm1,                pfrl,               pfrw,            &
                           pfri,                  loland,             pxtm1,           &
                           pxtte)

    USE mo_kind,       ONLY : wp
    USE mo_control,    ONLY : ltimer
    USE mo_timer,      ONLY : timer_start, timer_stop, timer_physc_subm_4
    USE mo_submodel,   ONLY : lco2, laoa
    USE mo_co2,        ONLY : co2_te_check, diag_co2
    USE mo_aoa,        ONLY : bcond_aoa, tracer_reset  ! age of air tracers

    INTEGER, INTENT(in)   :: kproma, kbdim, klev, klevp1, ktrac, krow
    REAL(wp), INTENT(in)  :: paphm1(kbdim,klevp1), & ! pressure at half levels t-dt
                             pfrl(kproma),         & ! fractional land coverage
                             pfrw(kproma),         & ! fractional water (ocean+lake) cov.
                             pfri(kproma),         & ! fractional ice coverage
                             pxtm1(kbdim,klev,ktrac) ! tracer mixing ratios at t-dt
    REAL(wp), INTENT(inout)::pxtte(kbdim,klev,ktrac) ! tracer tendencies
    LOGICAL, INTENT(in)   :: loland(kproma)          ! logical land mask

    IF (ltimer) CALL timer_start(timer_physc_subm_4)
    IF (lco2) THEN
       CALL co2_te_check(kproma,                 kbdim,               klev,            &
                         klevp1,                 ktrac,               krow,            &
                         paphm1,                 loland,              pxtm1,           &
                         pxtte)
    END IF
    CALL diag_co2(kproma,                                         krow,             &
                  pfrl,                      pfrw,                pfri)

    IF (laoa) THEN
       CALL tracer_reset(kproma, kbdim, klev, krow, pxtte)
       CALL bcond_aoa(kproma, kbdim, klev, krow, pxtte, pxtm1)
    END IF 

    IF (ltimer) CALL timer_stop(timer_physc_subm_4)

  END SUBROUTINE physc_subm_4
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! interface before the radiation calculations
!! modified from call_init_rad_submodels in call_submodels
!!
!! @author see module info of mo_submodel_interface
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info of mo_submodel_interface
!!
!! @par This subroutine is called by
!! radiation (rrtm_interface)
!!
!! @par Externals:
!! <ol>
!! <li>ham_rad_cache
!! <li>ham_rad
!! </ol>
!!
!! @par Notes
!! new radiation submodels must be introduced from here
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE radiation_subm_1(kproma           ,kbdim            ,klev         ,krow  ,&
                              ktrac            ,kaero            ,kpband       ,kb_sw ,&
                              aer_tau_sw_vr    ,aer_piz_sw_vr    ,aer_cg_sw_vr        ,&
                              aer_tau_lw_vr                                           ,&
                              ppd_hl           ,pxtm1                                  )



  USE mo_kind,           ONLY: wp
  USE mo_control,        ONLY: ltimer
  USE mo_timer,          ONLY: timer_start, timer_stop, &
                               timer_radiation_subm_1
  USE mo_submodel,       ONLY: lham
#ifdef HAMMOZ
  USE mo_ham_rad,       ONLY: ham_rad_cache,          &
                               ham_rad_cache_cleanup,  &
                               ham_rad
#endif

  IMPLICIT NONE

  ! arguments

  INTEGER, INTENT(in) :: kproma                     ! geographic block number of locations
  INTEGER, INTENT(in) :: kbdim                      ! geographic block maximum number of locations
  INTEGER, INTENT(in) :: klev                       ! numer of levels
  INTEGER, INTENT(in) :: ktrac                      ! numer of tracers
  INTEGER, INTENT(in) :: krow                       ! geographic block number
  INTEGER, INTENT(in) :: kaero                      ! switch for aerosol radiation coupling
  INTEGER, INTENT(in) :: kpband                     ! number of LW bands
  INTEGER, INTENT(in) :: kb_sw                      ! number of SW bands
  REAL(wp), INTENT(in):: ppd_hl(kbdim,klev)         ! pressure diff between half levels [Pa]
  REAL(wp), INTENT(in):: pxtm1(kbdim,klev,ktrac)    ! tracer mass/number mixing ratio (t-dt)   [kg/kg]/[

  REAL(wp), INTENT(inout) :: aer_tau_lw_vr(kbdim,klev,kpband),& !< LW optical thickness of aerosols
                             aer_tau_sw_vr(kbdim,klev,kb_sw), & !< aerosol optical thickness
                             aer_cg_sw_vr(kbdim,klev,kb_sw),  & !< aerosol asymmetry factor
                             aer_piz_sw_vr(kbdim,klev,kb_sw)    !< aerosol single scattering albedo


  IF (ltimer) CALL timer_start(timer_radiation_subm_1)

#ifdef HAMMOZ
  ! a) initialize HAM radiation submodule:
  IF (kaero==1) THEN !explicit HAM aerosol-radiation calculation
    ! for HAM-M7
    IF (lham) THEN
      CALL ham_rad_cache(kbdim, klev) !allocate memory
      CALL ham_rad(kproma, kbdim, klev, krow, kpband, kb_sw,                 &
                    pxtm1,         ppd_hl,                                    &
                    aer_tau_sw_vr, aer_piz_sw_vr, aer_cg_sw_vr, aer_tau_lw_vr )
    END IF
  END IF
#endif

  IF (ltimer) CALL timer_stop(timer_radiation_subm_1)


  END SUBROUTINE radiation_subm_1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Interface after the radiation calculations. 
!! Modified from call_finalize_rad_submodels in call_submodels. 
!!
!! @author see see module info of mo_submodel_interface
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info of mo_submodel_interface
!!
!! @par This subroutine is called by
!! radiation
!!
!! @par Externals:
!! <ol>
!! <li>ham_rad_diag
!! <li>ham_rad_cache_cleanup
!! </ol>
!!
!! @par Notes
!! 
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE radiation_subm_2(kproma, kbdim, krow, klev, &
                              ktrac, kaero,              &
                              pxtm1,                     &
!>>NAJS: lidar backscatter
                              ppd_hl,pp_fl,tk_fl,xm_vap)
!<<NAJS

  USE mo_kind,           ONLY: wp
  USE mo_control,        ONLY: ltimer
  USE mo_timer,          ONLY: timer_start, timer_stop, &
                               timer_radiation_subm_2
  USE mo_submodel,       ONLY: lham
#ifdef HAMMOZ
  USE mo_ham_rad,       ONLY: ham_rad_diag, ham_rad_cache_cleanup
#endif

  IMPLICIT NONE

  INTEGER,  INTENT(in) :: kproma                   ! geographic block number of locations
  INTEGER,  INTENT(in) :: kbdim                    ! geographic block maximum number of locations
  INTEGER,  INTENT(in) :: klev                     ! numer of levels
  INTEGER,  INTENT(in) :: ktrac                    ! numer of tracers
  INTEGER,  INTENT(in) :: krow                     ! geographic block number
  INTEGER,  INTENT(in) :: kaero                    ! geographic block number
  REAL(wp), INTENT(in) :: pxtm1 (kbdim,klev,ktrac) ! tracer mass/number mixing ratio (t-dt)   [kg/kg]/[#/k
!>>NAJS: lidar backscatter
  REAL(wp), INTENT(IN) :: ppd_hl(kbdim,klev)       ! level pressure thickness in Pa
  REAL(wp), INTENT(IN) :: pp_fl(kbdim,klev)        ! full level pressure in Pa
  REAL(wp), INTENT(IN) :: tk_fl(kbdim,klev)        ! full level temperature in K
  REAL(wp), INTENT(IN) :: xm_vap(kbdim,klev)       ! full level specific humidity in g/g
!<<NAJS


  IF (ltimer) CALL timer_start(timer_radiation_subm_2)

#ifdef HAMMOZ
  !--- 1) HAM radiation submodule:
  IF (kaero==1) THEN
    !--- Diagnostics:
    IF (lham) THEN
      CALL ham_rad_diag(kproma, kbdim, klev, krow, &
                        pxtm1,                     & 
!>>NAJS: lidar backscatter
                        ppd_hl,pp_fl,tk_fl,xm_vap)
!<<NAJS

      ! free memory for HAM aerosol optical properties:
      CALL ham_rad_cache_cleanup
    END IF
  END IF
#endif

  IF (ltimer) CALL timer_stop(timer_radiation_subm_2)

  END SUBROUTINE radiation_subm_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Interface of scan1
!! 
!!
!! @author see module info of mo_submodel_interface
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info of mo_submodel_interface
!!
!! @par This subroutine is called by
!! scan1
!!
!! @par Externals:
!! <ol>
!! </ol>
!!
!! @par Notes
!! 
!!
!! @par Responsible coder
!! sebastian.rast@mpimet.mpg.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scan1_subm

    USE mo_submodel,      ONLY : laoa
    USE mo_flighttrack,   ONLY: flighttrack_diag
    USE mo_submodel,      ONLY: lflighttrack
    USE mo_aoa,           ONLY : tf_reset

    IMPLICIT NONE
    IF (lflighttrack) CALL flighttrack_diag
    IF (laoa) CALL tf_reset

  END SUBROUTINE scan1_subm

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! interface between convective cloud physics and sulfur chemistry
!!
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! <ol>
!! <li> S. Ferrachat, ETH Zurich. Wetdep_interface call (2009.08; 2009.12)
!! </ol>
!!
!! @par This subroutine is called by
!! cuflx
!!
!! @par Externals:
!! <ol>
!! </ol>
!!
!! @par Notes
!! The wet deposition interface is also called from cloud_subm.
!!
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Interface for wet chemical processes called from cuflx
  !! calculates liquid fraction of SO2 (HAM) and calls the wet deposition interface
!>>>mgs: modified call for lightning routine

  SUBROUTINE cuflx_subm(kbdim, kproma,   klev,     ktop,      krow,      &
#ifdef HAMMOZ /* SF */
!SF mgs for lightning
                        ktype, kcbot,    kctop,                          & ! conv. cloud extent
#endif
                        pxtenh, pxtu, prhou,                             & ! tracers
                        pmfu,   pmfuxt,                                  & ! convective fluxes and corresp. mmr
                        pmlwc,  pmiwc,     pmratepr, pmrateps,           & ! cloud properties
                        pfrain, pfsnow,    pfevapr,  pfsubls,            & !   "       "
                        paclc,  pmsnowacl,                               & !   "       "           !--mgs: lstrat remo
                        ptu,    pdpg,                                    & ! thermodynamic quantities
#ifdef HAMMOZ /* SF */
!SF mgs for lightning
                        pten,   ptenh,     paphp1,                       & !
#endif
                        pxtte                                            &
#ifdef HAMMOZ /* SF */
!SF for wetdep
                       ,pxtbound                                         &
#endif
                       )
!<<<mgs


  USE mo_kind,               ONLY: wp
  USE mo_control,            ONLY: ltimer
  USE mo_timer,              ONLY: timer_start,            &
                                   timer_stop,             &
                                   timer_cuflx_subm 
  !>>dod split mo_timer (redmine #51)
  USE mo_hammoz_timer,       ONLY: timer_hammoz_wetdep
  !<<dod
  USE mo_tracdef,            ONLY: ntrac, trlist
  USE mo_submodel,           ONLY: lwetdep, llght, lchemistry, lham     !>>>mgs<<<
#ifdef HAMMOZ
  USE mo_hammoz_wetdep,      ONLY: wetdep_interface
  USE mo_ham_wetdep,         ONLY: ham_conv_lfraq_so2
  USE mo_moz_lightning,      ONLY: moz_lightning       !>>>mgs<<<
  USE mo_memory_g3b,         ONLY: slm                 !>>>mgs<<<
  USE mo_exception,          ONLY: message, em_debug
#endif

  !--- Arguments:

  INTEGER, INTENT(in)     :: kbdim, kproma, klev, ktop, & ! dimensions
                             krow
!>>SF mgs for lightning
#ifdef HAMMOZ
INTEGER, INTENT(in)       :: ktype(kbdim),              & ! convective cloud type
                             kcbot(kbdim), kctop(kbdim)   ! cloud extent in levels
#endif
!<< SF mgs

  REAL(wp), INTENT(in)    :: pdpg(kbdim,klev),          & ! geopotential height
                             pmratepr(kbdim,klev),      & !
                             pmrateps(kbdim,klev),      &
                             pmsnowacl(kbdim,klev),     &
                             ptu(kbdim,klev),           & ! temperature at beginning of time step
                             pfrain(kbdim,klev),        &
                             pfsnow(kbdim,klev),        &
                             pfevapr(kbdim,klev),       &
                             pfsubls(kbdim,klev),       &
                             pmfu(kbdim,klev),          & !
                             paclc(kbdim,klev),         & ! 3D cloud cover
                             prhou(kbdim,klev)            ! air density
!>>SF mgs for lightning
#ifdef HAMMOZ
  REAL(wp), INTENT(inout) :: pten(kbdim,klev),          &
                             ptenh(kbdim,klev),         &
                             paphp1(kbdim,klev+1)         ! half level pressures
#endif
!<<SF mgs for lightning

  REAL(wp), INTENT(inout) :: pxtte(kbdim,klev,ntrac),   & ! tracer tendency
                             pmlwc(kbdim,klev),         & ! liquid water content in ???
                             pmiwc(kbdim,klev),         & ! ice water content in ???
                             pxtenh(kbdim,klev,ntrac),  & ! ambient tracer mmr at t=t+dt
                                                          ! units: [kg(tracer) kg-1(air)] in/out
                                                          !        [kg(tracer) m-2]       inside
                             pxtu(kbdim,klev,ntrac),    & ! in-cloud tracer mmr at t=t+dt
                                                          ! units as pxtenh
                             pmfuxt(kbdim,klev,ntrac)    
!>>SF
#ifdef HAMMOZ
  REAL(wp), INTENT(inout) :: pxtbound(kbdim,ntrac)        !SF boundary condition for xt_conv_massfix
#endif
!<<SF

  !--- Local variables
  LOGICAL, PARAMETER  :: lstrat = .FALSE.    !mgs: hardcoded here for call to wetdep_interface

  ! for scavenging
  REAL(wp)            :: zdummy(kbdim,klev,ntrac), & ! for wetdep_interface call (pxtm1)
                         zdum2d(kbdim,klev),       & ! dto. (pclc)
                         zlfrac_so2(kbdim,klev)      ! liquid fraction of SO2 (for HAM)

!>>>mgs
  ! for lightning
  REAL(wp)            :: znoems_3d(kbdim, klev)      ! lightning NO emissions (MOZ)
!<<<mgs


  IF (ltimer) CALL timer_start(timer_cuflx_subm)

  ! 1)-- calculate liquid fraction of SO2 for HAM
  zlfrac_so2(1:kproma,:) = 0._wp

#ifdef HAMMOZ
  IF (lham .AND. lchemistry) THEN
    IF (ltimer) CALL timer_start(timer_hammoz_wetdep)
    CALL ham_conv_lfraq_so2(kproma, kbdim, klev, &
                            ptu, pxtu, prhou,    &
                            pmlwc, zlfrac_so2    )
    IF (ltimer) CALL timer_stop(timer_hammoz_wetdep)

  END IF

  ! 2)-- call wet deposition routine
  IF (lwetdep .AND. ANY(trlist%ti(:)%nwetdep > 0)) THEN
    IF (ltimer) CALL timer_start(timer_hammoz_wetdep)

    zdummy(1:kproma,:,:) = 0._wp
    zdum2d(1:kproma,:)   = 0._wp

    CALL wetdep_interface(kproma, kbdim,    klev,     ktop,      krow, lstrat,  &
                          pdpg,   pmratepr, pmrateps, pmsnowacl,                &
                          pmlwc,  pmiwc,                                        &
                          ptu,    zdummy,   zlfrac_so2, pxtte,    pxtenh, pxtu, &
                          pfrain, pfsnow,   pfevapr,  pfsubls,                  &
                          pmfu,   pmfuxt,                                       &
                          paclc,  zdum2d,   prhou, pxtbound)
    IF (ltimer) CALL timer_stop(timer_hammoz_wetdep)
  END IF
  IF (llght) THEN
!!  CALL message('cuflx_subm','Calling moz_lightning from cuflx_subm', level=em_debug)
    CALL moz_lightning(kproma,      kbdim,      klev,        krow,              &
                       ktype(:),    kcbot(:),   kctop(:),                       &
                       slm(:,krow), pten(:,:),  ptenh(:,:),  paphp1(:,:),       &
                       pmfu(:,:),   znoems_3d(:,:)                          )
    IF (idt_NO > 0) pxtte(1:kproma,:,idt_NO) = pxtte(1:kproma,:,idt_NO) + znoems_3d(1:kproma,:)
  END IF
#endif

  IF (ltimer) CALL timer_stop(timer_cuflx_subm)

  END SUBROUTINE cuflx_subm

END MODULE mo_submodel_interface

