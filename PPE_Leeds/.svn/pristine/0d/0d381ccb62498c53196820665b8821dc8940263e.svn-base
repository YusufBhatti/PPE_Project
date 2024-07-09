#ifdef __xlC__
@PROCESS STRICT
#endif

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE stepon 
  !
  ! Description:
  !
  ! Controls the time step.
  !
  ! Method:
  !
  ! This subroutine controls the structure of the scanning
  ! over the latitude lines and of the computations in spectral
  ! space. It also increments the time step and check for the
  ! completion of the run.
  !
  ! *stepon* is called from *control*.
  !
  ! Externals:
  ! *scan1*     1st scans over gaussian latitudes.
  ! *lti*       inverse Legendre transform.
  ! *hdiff*     horizontal diffusion.
  ! *scctp*     computations in spectral space for
  !             temperature and surface pressure equations.
  ! *sccd*      computations in spectral space for
  !             divergence equation.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, March 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, August 1998, tendency diagnostics, nudging and nmi
  ! I. Kirchner, MPI, January 1999, add nmi
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! U. Schlese, DKRZ, October 1999, ECHAM5-modifications
  ! I. Kirchner, MPI, October 2000, date/time control
  ! S. Legutke, MPI,M&D, Jan 2002, coupling interface
  ! I. Kirchner, MPI, Aug 2002, nudging revision
  ! L. Kornblueh, MPI, Apr 2003, time control changes
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! M. Schulz, FZJ, July 2010,  moved pre_radiation call from scan1 so 
  !                 that submodels can use mo_radiation_parameters and 
  !                 activated bc_list_read and added timer
  ! D. Klocke, MPI, Nov 2010, add lnwp for NWP restarts
  !
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, message_text, finish, number_of_errors
  USE mo_control,       ONLY: lmidatm, lcouple, lcolumn, lprint_m0,      &
                              lnudge, nlev,                              &
                              ltimer, ltctest, lnwp,                     &
                              lforcererun !SF: switchable internal reruns in case of HAMMOZ
  USE mo_memory_sp,     ONLY: stp, sd, svo
  USE mo_hyb,           ONLY: aktlrd, altrcp, rpr
  USE mo_hdiff,         ONLY: hdiff
  USE mo_upper_sponge,  ONLY: uspnge

  USE mo_nudging,       ONLY: Nudging
  USE mo_nudging_init,  ONLY: NudgingInit, NDG_INI_MEM
  USE mo_nudging_sst,   ONLY: NudgingReadSST
  USE mo_nudging_buffer,ONLY: nio_index, NdgSetCounter, NdgRemoveCounter, &
                              NdgCorrBuffer, NdgCleanBuffer

  USE mo_nmi,           ONLY: NMI_Make, NMI_MAKE_NMI, lnmi_run
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_couple,        ONLY: couple_get_o2a, couple_put_a2o, couple_calendar

  USE mo_time_control,  ONLY: time_set, time_reset,                          &
                              lstart, lresume, lstop, lbreak,                &
                              l2nd_day, l_trigrad, get_time_step,            &
                              nsub, l_trigjob, l_putdata, l_putrerun,        &
                              l_putcheckpoint,                               &
                              l_trigfiles, subflag, delta_time

  USE mo_output,        ONLY: open_output_streams, close_output_streams, &
                              out_streams
  USE mo_io,            ONLY: write_streams

  USE mo_timer,         ONLY: timer_stop, print_timer, timer_total,      &
                              timer_start, timer_output, timer_restart,  &
                              timer_bclistread, timer_loop,              &
                              timer_time_set, timer_trigfiles,           &
                              timer_bcond, timer_prerad, timer_subm,     &
                              timer_jsbach, timer_scan1, timer_sccd,     &
                              timer_scctp, timer_uspnge, timer_hdiff,    &
                              timer_inv_legendre, timer_time_reset,      &
                              timer_couple_get, timer_couple_put,        &
                              timer_hd

  USE mo_co2,           ONLY: co2_flux_correction
  
  USE mo_time_control,  ONLY: current_date, next_date, prev_radiation_date, &
                              radiation_date, get_date_components
  USE mo_sst,           ONLY: readsst, readice, readdailysst, readdailyice
  USE mo_control,       ONLY: lamip, ldailysst, ltimer, ltdiag

  USE mo_boundary_condition, ONLY: bc_list_read
  USE mo_aero_kinne,         ONLY: read_aero_kinne
  USE mo_aero_volc,          ONLY: read_aero_volc
  USE mo_aero_volc_tab,      ONLY: read_aero_prop_ham, read_aero_prop_crow
  USE mo_o3clim,             ONLY: read_o3clim_4
  USE mo_radiation_parameters,          ONLY: iaero,io3, jpsw => nb_sw
  USE mo_radiation,          ONLY: pre_radiation
  USE mo_submodel,           ONLY: lanysubmodel
  USE mo_submodel_interface, ONLY: stepon_subm
  USE mo_jsbach_interface,   ONLY: stepon_jsbach, stepoff_jsbach
  USE mo_input,              ONLY: InputClose 
  USE mo_memory_cfdiag,      ONLY: locfdiag, calc_avps
  USE mo_mvstream,           ONLY: mvstream_accumulate
  USE mo_util_logging,       ONLY: send_network_log
  USE mo_diag_tendency_new,  ONLY: tdiag_vars, set_tendency
  USE mo_call_trans,    ONLY: spectral_to_legendre 

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: idt, jlev, jsu, istep, nb_sw
  REAL(dp):: ztw0, ztw1

  !  External functions 
  REAL(dp), EXTERNAL :: util_walltime

  LOGICAL :: lnewyear, lnewday
!>>SF: switchable internal rerun in case of HAMMOZ
  LOGICAL :: lexit
!<<SF 
  INTEGER :: icurrentyear, icurrentmonth, icurrentday
  INTEGER :: inextyear, inextmonth, inextday

  !  External subroutines 
  EXTERNAL :: helmo, scan1, lti, sccd, scctp, subjob

  !  Executable statements

  IF (number_of_errors > 0) THEN
    WRITE (message_text, *) number_of_errors, ' errors encountered. Abort run now.'
    CALL finish('stepon', message_text)
  END IF

  IF (ltimer) THEN
    CALL timer_start(timer_total)
  ENDIF

  integration_loop: DO

    IF (ltimer) CALL timer_start(timer_loop)

     ztw0 = util_walltime(1)

     IF (ltimer) CALL timer_start(timer_time_set)     
     CALL time_set
     IF (ltimer) CALL timer_stop(timer_time_set)     

     IF (.NOT. ltctest) THEN     ! time control testing
 
     !-- run NMI part 1
 
     IF (lnmi_run) CALL NMI_Make(NMI_MAKE_NMI)

     IF (lnudge) THEN
        CALL NudgingInit(NDG_INI_MEM)
        CALL NudgingReadSST
     END IF

     IF (l_trigfiles) THEN
       IF (ltimer) CALL timer_start(timer_trigfiles)
       CALL close_output_streams
       IF (.NOT.lcolumn) THEN
         DO jsu = 1, nsub  ! filter subjobs connected to output streams
           IF(l_trigjob(jsu).AND.subflag(jsu)) CALL subjob(jsu)
         END DO
       ENDIF
       CALL open_output_streams
       IF (ltimer) CALL timer_stop(timer_trigfiles)
     ENDIF


     IF (ltimer) CALL timer_start(timer_bcond)

     IF (lamip) THEN
        CALL get_date_components (current_date, month=icurrentmonth, year=icurrentyear,  &
                              day=icurrentday                                            )
        CALL get_date_components (next_date,  month=inextmonth, year=inextyear,          &
                              day=inextday                                               )
        lnewyear=icurrentyear/=inextyear
        IF (lnewyear) THEN 
            CALL readsst
            CALL readice
        ENDIF
     ENDIF

     IF (ldailysst) THEN
        CALL get_date_components (current_date, month=icurrentmonth, year=icurrentyear,  &
                              day=icurrentday                                            )
        CALL get_date_components (next_date,  month=inextmonth, year=inextyear,          &
                              day=inextday                                               )
        lnewday=icurrentday/=inextday
        IF (lnewday) THEN
            CALL readdailysst
            CALL readdailyice
        ENDIF
     ENDIF

     IF (ltimer) CALL timer_start(timer_bclistread)
     CALL bc_list_read
     IF (ltimer) CALL timer_stop(timer_bclistread)

     ! read boundary transient boundary conditions for radiation

     IF (l_trigrad) THEN
       nb_sw = jpsw
       IF (iaero==3) THEN
         CALL read_aero_kinne(prev_radiation_date, radiation_date, nb_sw)
       END IF
       IF (iaero==5) THEN
         CALL read_aero_kinne(prev_radiation_date, radiation_date, nb_sw)
         CALL read_aero_volc(prev_radiation_date, radiation_date, nb_sw)
       END IF
       IF (iaero==6) THEN
         CALL read_aero_kinne(prev_radiation_date, radiation_date, nb_sw)
         CALL read_aero_volc(prev_radiation_date, radiation_date, nb_sw)
         CALL read_aero_prop_ham(prev_radiation_date, radiation_date)
       END IF
       IF (iaero==7) THEN
         CALL read_aero_kinne(prev_radiation_date, radiation_date, nb_sw)
         CALL read_aero_prop_crow(prev_radiation_date, radiation_date)
       END IF
       IF (io3==4) THEN
          CALL read_o3clim_4(prev_radiation_date, radiation_date)
       END IF
     END IF
     IF (ltimer) CALL timer_stop(timer_bcond)

     IF (ltimer) CALL timer_start(timer_prerad)
     CALL pre_radiation
     IF (ltimer) CALL timer_stop(timer_prerad)

     IF (lanysubmodel) THEN
       IF (ltimer) CALL timer_start(timer_subm)
       CALL stepon_subm (current_date, next_date)
       IF (ltimer) CALL timer_stop(timer_subm)
     END IF
     
     IF (.NOT. lstart .AND. .NOT. lresume) THEN
       IF (ltimer) CALL timer_start(timer_jsbach)
       CALL stepon_jsbach
       IF (ltimer) CALL timer_stop(timer_jsbach)
     ENDIF
     
     ! --- If exchange-data read event, read data 

     IF (lcouple) THEN
       IF (ltimer) CALL timer_start(timer_couple_get)       
       CALL couple_get_o2a
       IF (ltimer) CALL timer_stop(timer_couple_get)       
     ENDIF

     IF (ltimer) CALL timer_start(timer_scan1)
     CALL scan1
     IF (ltimer) CALL timer_stop(timer_scan1)

     !-- 1.2 Completion of divergence calculation

     IF (ltimer) CALL timer_start(timer_sccd)
     CALL sccd
     IF (ltimer) CALL timer_stop(timer_sccd)

     !-- 1.3 Completion of temperature and surface pressure equations.

     IF (ltimer) CALL timer_start(timer_scctp)
     CALL scctp
     IF (ltimer) CALL timer_stop(timer_scctp)

     !-- 1.4 Upper sponge layer 

     IF (lmidatm) THEN
       IF (ltimer) CALL timer_start(timer_uspnge)
       CALL uspnge
       IF (ltimer) CALL timer_stop(timer_uspnge)
     ENDIF

     !-- 1.5 Horizontal diffusion

     IF (ltdiag) THEN
       IF (ASSOCIATED(tdiag_vars%dsddt_hdiff)) THEN
         CALL set_tendency(tdiag_vars%dsddt_hdiff,sd,'sub')
       END IF
       IF (ASSOCIATED(tdiag_vars%dsvodt_hdiff)) THEN
         CALL set_tendency(tdiag_vars%dsvodt_hdiff,svo,'sub')
       END IF
       IF (ASSOCIATED(tdiag_vars%dstdt_hdiff)) THEN
         CALL set_tendency(tdiag_vars%dstdt_hdiff,stp(1:nlev,:,:),'sub')
       END IF
     END IF

     IF (ltimer) CALL timer_start(timer_hdiff)
     CALL hdiff
     IF (ltimer) CALL timer_stop(timer_hdiff)

     IF (ltdiag) THEN
       IF (ASSOCIATED(tdiag_vars%dsddt_hdiff)) THEN
         CALL set_tendency(tdiag_vars%dsddt_hdiff,sd,'add')
       END IF
       IF (ASSOCIATED(tdiag_vars%dsvodt_hdiff)) THEN
         CALL set_tendency(tdiag_vars%dsvodt_hdiff,svo,'add')
       END IF
       IF (ASSOCIATED(tdiag_vars%dstdt_hdiff)) THEN
         CALL set_tendency(tdiag_vars%dstdt_hdiff,stp(1:nlev,:,:),'add')
       END IF
     END IF

     ! -- Call nudging after the horizontal diffusion

     IF (lnudge) THEN
       CALL Nudging

     ELSE IF (lnmi_run) THEN
       !-- run NMI part 2
       CALL NMI_Make(NMI_MAKE_NMI)
       
     END IF

     ! -- Postprocessing of spectral data

     IF (lnudge) THEN
       IF(l_putdata(nio_index)) THEN
         CALL NdgCorrBuffer
       END IF
     END IF

!cms++     
     IF ( locfdiag ) THEN
        CALL calc_avps
     END IF   
!cms--

     ! --- jsbach routines that need to be updated at the end of the time step
     CALL stepoff_jsbach

     IF (lcouple) THEN
       IF (ltimer) CALL timer_start(timer_couple_put)       
       CALL couple_put_a2o
       IF (ltimer) CALL timer_stop(timer_couple_put)       
     ENDIF

     ! Store current values into accumulation buffers.
     CALL mvstream_accumulate

     IF (ltimer) CALL timer_start(timer_output)
     CALL out_streams
     IF (ltimer) CALL timer_stop(timer_output)

     IF (lnudge) THEN
       IF(l_putdata(nio_index)) CALL NdgCleanBuffer
     END IF

     ! --- Update coupling time

     IF (lcouple) THEN
        CALL  couple_calendar(delta_time)
     END IF

     ! --- Compute flux correction for CO2 tracer
     CALL co2_flux_correction

     !-- 2. Continuation of run

     !-- 2.0 Inverse Legendre transforms

     IF (ltimer) CALL timer_start(timer_inv_legendre)
     CALL spectral_to_legendre
     CALL lti
     IF (ltimer) CALL timer_stop(timer_inv_legendre)

     !-- 2.1 Recompute matrix *cn* (used to compute divergence
     !       in sccd) after the first time-step.
     
     IF (lstart .OR. lnwp) THEN
        idt = 2
        CALL helmo(idt)
        
        !-- 2.2 Multiply by 2. arrays *aktlrd* and *altrcp* (used
        !       by conteq) and *rpr* after the first time step.

        DO jlev = 1, nlev
           aktlrd(jlev) = aktlrd(jlev)*2._dp
           altrcp(jlev) = altrcp(jlev)*2._dp
        END DO
        rpr = rpr*2._dp

     END IF

     END IF      ! end if for ltctest - time control test

     ! -- store final parts of rerun files

     IF ((l_putrerun .OR. l_putcheckpoint) .AND..NOT.lcolumn) THEN
       IF (lnudge) CALL NdgSetCounter
       IF (ltimer) CALL timer_start(timer_restart)
       CALL write_streams
       IF (ltimer) CALL timer_stop(timer_restart)
       IF (lnudge) CALL NdgRemoveCounter
     END IF

     ! process subjobs not connected to output streams
     DO jsu = 1, nsub
       IF(l_trigjob(jsu).AND. .NOT.subflag(jsu)) CALL subjob(jsu)
     END DO

     IF (ltimer) CALL timer_start(timer_time_reset)
     CALL time_reset
     IF (ltimer) CALL timer_stop(timer_time_reset)

     istep = get_time_step()

     IF (ldc%nsnm0 > 0) THEN   
       IF (ldc%snn0(1) == 0) THEN   
         WRITE (message_text,'(a,i4,a,f7.2,a)') &
              ' PE',ldc%pe, ' - global mean surface temperature: ', stp(nlev,1,1), ' K'
         CALL send_network_log(message_text)
       ENDIF
     ENDIF
     
     IF (lprint_m0 .OR. l2nd_day) THEN
       IF (ldc%nsnm0 > 0) THEN   
         IF (ldc%snn0(1)==0) THEN   
           WRITE (message_text,'(a,i4,a,i10,f20.13," K")') &
                ' PE',ldc%pe,' stepon: ', istep, stp(nlev,1,1)
           CALL message('',message_text)
         END IF
       END IF

       ztw1 = util_walltime(1)

       WRITE (message_text,'(a,i4,a,i10,f10.3," s")') &
            ' PE',ldc%pe,' stepon: ', istep, (ztw1-ztw0)
       CALL message('',message_text)
     END IF

     IF (ltimer) CALL timer_stop(timer_loop)

     !-- 2.3   Test for end of run

     IF (number_of_errors > 0) THEN
       WRITE (message_text, *) number_of_errors, ' errors encountered. Abort run now.'
       CALL finish('stepon', message_text)
     END IF

!>>SF: switchable internal rerun 
    IF (lforcererun) THEN
       lexit = l_putrerun .OR. lbreak .OR. lstop !default for ECHAM
    ELSE
       lexit = (l_putrerun .AND. lbreak) .OR. lstop !compulsory for HAMMOZ
    ENDIF
!<<SF 

     !SF IF (l_putrerun .OR. lbreak .OR. lstop) THEN
     IF (lexit) THEN
       WRITE (message_text,'(a,i7,a)') 'Step ', istep, ' completed.'
       CALL message('',message_text)
       EXIT integration_loop
     END IF

  END DO integration_loop

  CALL InputClose()

  IF (.NOT.lcolumn) THEN
    ! filter subjobs connected to output streams
    DO jsu = 1, nsub
      IF(l_trigjob(jsu).AND.subflag(jsu)) CALL subjob(jsu)
    END DO
  END IF

  IF (ltimer) THEN
    CALL timer_stop(timer_total)
    CALL message('stepon','Stop Integration loop timer ...')
#ifdef DEBUG
    CALL print_timer()
#else
    CALL print_timer(short=.TRUE.)
#endif
  END IF

END SUBROUTINE stepon
