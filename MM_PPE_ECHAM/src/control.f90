!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE control

  ! Description:
  !
  ! Control routine for the model.
  !
  ! Method:
  !
  ! This subroutine controls the running of the model.
  !
  ! *control* is called from the main program (*master*).
  !
  ! Externals:
  !   *initialize*    called to initialize modules.
  !   *stepon*        controls time integration
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! R.G and M.J, ECMWF, December 1982, changed
  ! U. Schlese, MPI, August 1989, new structure
  ! U. Schlese, DKRZ, September 1994, interface *drive* added
  ! U. Schlese, DKRZ, January 1995, reading of optional files added
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ and M. Esch, MPI, July 1999, modifications for ECHAM5
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! U. Schlese, DKRZ, November 1999, "drive" removed, "stepon" called directly
  ! U. Schlese, DKRZ, December 1999, modifications for coupling
  ! M.A. Giorgetta, MPI, May 2000, ozone initialization removed, see setrad
  ! S. Legutke, MPI M&D, July 00, modifications for coupling interface
  ! I. Kirchner, MPI, December 2000, date/time control/nudging
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! I. Kirchner, MPI, Aug 2002, tendency diagnostics revision
  ! A. Rhodin, DWD, June 2002, call subroutines cleanup_...
  ! M. Esch, MPI, September 2002, modifications for mixed layer ocean
  ! L. Kornblueh, MPI, October 2003, added setup for AMIP2 global diagnostics
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! U. Schlese, MPI, November 2008, setting of "iaero" removed, "lso4" removed
  ! D. Klocke, MPI, Nov 2010, added switch for NWP restarts
  ! L. Kornblueh, MPI, Februrary 2012, remove spitfire
  ! S.K. Cheedela, MPI, March 2010, initalize and cleanup single column
  !
  ! for more details see file AUTHORS
  !

  USE mo_machine,         ONLY: machine_setup
  USE mo_control,         ONLY: lcolumn, lcouple, lnudge,               &
                                lnmi, lmlo, ldiagamip, ldailysst, lnwp, &
                                ltimer, nn, nlev,                       &
                                nproca, nprocb, nprocio
  USE mo_time_control,    ONLY: lstart, lfirst_cycle, lresume,          &
                                print_events, lstop, putrerun,          &
                                delta_time, get_cycle
  USE mo_sst,             ONLY: readsst, readice, readflux,             &
                                readdailysst, readdailyice, cleanup_sst
  USE mo_clim,            ONLY: readtslclim, readvltclim, readvgratclim
  USE mo_legendre,        ONLY: inileg, cleanup_legendre
  USE mo_nmi,             ONLY: NMI_Init, NMI_Close
  USE mo_nudging_init,    ONLY: NudgingInit, NDG_CLEAN_MEM, NDG_CLOSE
  USE mo_jsbach_interface, ONLY: stepon_jsbach
  USE mo_diag_amip2,      ONLY: init_amip2_diag
  USE mo_couple,          ONLY: couple_init
  USE mo_column,          ONLY: init_column, resetcolumn
  USE mo_output,          ONLY: init_output, cleanup_output, &
                                open_output_streams, close_output_streams
  USE mo_linked_list,     ONLY: print_stream
  USE mo_memory_base,     ONLY: ostreams, nstreams
  USE mo_memory_streams,  ONLY: init_memory, free_memory
  USE mo_geoloc,          ONLY: init_geoloc, cleanup_geoloc
  USE mo_advection,       ONLY: iadvec, tpcore, semi_lagrangian
  USE mo_semi_lagrangian, ONLY: init_semi_lagrangian,                  &
                                cleanup_semi_lagrangian
  USE mo_tpcore,          ONLY: init_tpcore, cleanup_tpcore
  USE mo_timer,           ONLY: init_timer, cleanup_timer
  USE mo_radiation_parameters,       ONLY: iaero, io3
  USE mo_param_switches,  ONLY: iconv
  USE mo_decomposition,   ONLY: ldc=>local_decomposition,              &
                                cleanup_decomposition
  USE mo_exception,       ONLY: message_text, message, finish
  USE mo_scan_buffer,     ONLY: cleanup_scanbuffer
  USE mo_o3clim,          ONLY: read_o3clim_4, cleanup_o3clim
  USE mo_clim,            ONLY: cleanup_clim

  USE mo_tmp_buffer,      ONLY: cleanup_tmp_buffer
  USE mo_gaussgrid,       ONLY: cleanup_gaussgrid
#ifdef FFT991
  USE mo_fft991,          ONLY: cleanup_fft991
#elif MKL_DFT
  USE mo_mkl_dft,         ONLY: cleanup_mkl_dft
#elif ESSL_DFT
  USE mo_essl_dft,        ONLY: cleanup_essl_dft
#else
  USE mo_fft992,          ONLY: cleanup_fft992
#endif
  USE mo_call_trans,      ONLY: spectral_to_legendre
  USE mo_hdiff,           ONLY: ldiahdf, init_hdiff_diag, finalize_hdiff_diag
  USE m_alloc_mods,       ONLY: dealloc_mods
  USE mo_netcdf,          ONLY: cleanup_netcdf
  USE mo_io,              ONLY: cleanup_io
  USE mo_greenhouse_gases,ONLY: cleanup_greenhouse_gases
  USE mo_solar_irradiance,ONLY: cleanup_solar_irradiance
  USE mo_so4,             ONLY: cleanup_so4
  USE mo_aero_kinne,      ONLY: cleanup_aero_kinne
  USE mo_aero_volc,       ONLY: cleanup_aero_volc
  USE mo_aero_volc_tab,   ONLY: cleanup_aero_volc_tab_ham
  USE mo_aero_volc_tab,   ONLY: cleanup_aero_volc_tab_crow
  USE mo_station_diag,    ONLY: cleanup_station_diag, lostation
  USE mo_surface,         ONLY: init_surface
  USE mo_truncation,      ONLY: cleanup_truncation
  USE mo_version,         ONLY: init_version, label_run
  USE mo_filename,        ONLY: out_expname, out_datapath
  USE mo_util_db_timings, ONLY: set_experiment_name, set_job_name, &
                                set_job_data
  USE mo_tracdef,         ONLY: ntrac

  USE mo_station_diag,    ONLY: lostation, init_station_diag
  USE mo_cosp_offline,    ONLY: locospoffl, init_cosp_offline

#ifdef _OPENMP
  USE omp_lib,            ONLY: omp_get_num_threads
#endif
  IMPLICIT NONE

  INTEGER :: i, nthreads, ncycle

  !  External subroutines
  EXTERNAL :: stepon, initialize, inhysi, ioinitial, lti

  !  Executable statements

  WRITE (message_text,'(a)') &
       '======================================================================'
  CALL message('',message_text)
  WRITE (message_text,'(a)') &
       'Start initialization in control'
  CALL message('',message_text)
  WRITE (message_text,'(a)') &
       '======================================================================'
  CALL message('',message_text)


  ! Section declaring default selection of multiple available schemes:
  !
  ! - advection scheme
  iadvec = tpcore
  ! - convection
  iconv = 1

  !--  Print machine specific values

  CALL machine_setup

  ! set model version strings as far as available

  CALL init_version

  ! initialize generic model timers (poor mans profiling)

  CALL init_timer

  ! initialize modules and parallel decomposition

  CALL initialize

  ! initialize time independent surface parameters

  CALL init_geoloc(io3,iaero)

  ! error handling with respect to nproma

  IF ( .NOT. ldc% lreg ) THEN
    IF ( iadvec .EQ. semi_lagrangian .OR. lnudge ) THEN
      WRITE(message_text,*) 'For running'
      CALL message('control', TRIM(message_text))
      WRITE(message_text,*) '  - semi_lagrangian advection (iadvec=1)'
      CALL message('control', TRIM(message_text))
      WRITE(message_text,*) '  - nudging (lnudge=T)'
      CALL message('control', TRIM(message_text))
      CALL finish('control','NPROMA is not implemented. Delete it from namelist runctl.')
    END IF
  END IF

  ! initialize memory

  CALL init_memory

  ! preset values needed in the advection scheme

  SELECT CASE (iadvec)
  CASE (semi_lagrangian)
    CALL init_semi_lagrangian
  CASE (tpcore)
    CALL init_tpcore
  END SELECT


  ! compute Legendre polynomials and parameters needed for the Legendre transforms.

  CALL inileg

  ! read restart/initial conditions

  IF (lresume .AND. .NOT. lnwp) THEN
     CALL iorestart
  ELSE IF (lstart .AND. .NOT. lnwp) THEN
     CALL ioinitial
  ELSE IF (lnwp) THEN
    CALL ionwp
  END IF

  !> read boundary, some initial fields for Single Column

  IF(lcolumn) CALL init_column ! pending fields with JSBach

  CALL stepon_jsbach ! Needs to be called the first time after initial or restart
                     ! files are read, and before init_surface!

  CALL init_surface

  CALL inhysi

  IF (lstart .OR. lnwp) THEN
    CALL spectral_to_legendre
    CALL lti
  END IF

  IF (lnmi .AND. lfirst_cycle) CALL NMI_Init(lnudge)

  IF (ldiagamip) CALL init_amip2_diag

  ! read optional sst- and seaice-file if not in coupled mode

  IF (.NOT. lcouple) THEN
    IF ( ldailysst ) THEN
      CALL readdailysst
      CALL readdailyice
    ELSE
      CALL readsst
      CALL readice
    END IF
  END IF

  ! read new ozone climatalogy

  IF (io3 == 4) THEN
    CALL read_o3clim_4
  ENDIF

  ! read optional flux correction if in mixed layer mode

  IF (lmlo) CALL readflux

  ! read  climate land surface temperatures

  CALL readtslclim

  ! read climate leaf-area index and climate vegetation ratio

  CALL readvltclim

  CALL readvgratclim

  ! initialize coupled run:
  !  - initialize  data exchange with coupler
  !  - check control parameters for consistency with coupler

  IF (lcouple .AND. lfirst_cycle) THEN
    CALL couple_init
  END IF

  ! define CDI resources here in STAGE_DEFINITION
  CALL init_output
  CALL open_output_streams

#ifdef HAVE_CDIPIO
  ! change from I/O STAGE_DEFINITION to STAGE_TIMELOOP
  IF (nprocio > 0 ) THEN
    CALL pioEndDef
  ENDIF
#endif

  ! label run

  CALL label_run

  ! current cycle number

  ncycle = get_cycle()

  ! setup of database writing of timers

  IF (ltimer) THEN
#ifdef _OPENMP
!$OMP PARALLEL
    nthreads = omp_get_num_threads()
!$OMP END PARALLEL
#else
    nthreads = 1
#endif

    CALL set_experiment_name(out_expname, out_datapath)
    CALL set_job_name(ncycle)
    CALL set_job_data(nproca, nprocb, nprocio, nthreads,    &
                      putrerun%unit, putrerun%counter,      &
                      NINT(delta_time), nn, nlev, ncycle)
  ENDIF

  ! print status of streams

  IF (ncycle == 1) THEN 
    CALL message('','')
    DO i = 1, nstreams
      CALL print_stream (ostreams (i))
    ENDDO
  ENDIF
  ! prints events for model run (model event handling)

  CALL print_events

  ! start hdiff diagnostics, if requested

  IF (ldiahdf) CALL init_hdiff_diag

  IF (lostation) CALL init_station_diag

  IF (locospoffl) CALL init_cosp_offline 

  WRITE (message_text,'(a)') &
       '======================================================================'
  CALL message('',message_text)
  WRITE (message_text,'(a)') &
       'Finished initialization in control'
  CALL message('',message_text)
  WRITE (message_text,'(a)') &
       '======================================================================'
  CALL message('',message_text)

  ! start time integration

  CALL stepon

  WRITE (message_text,'(a)') &
       '======================================================================'
  CALL message('',message_text)
  WRITE (message_text,'(a)') &
       'Start cleanup in control'
  CALL message('',message_text)
  WRITE (message_text,'(a)') &
       '======================================================================'
  CALL message('',message_text)

#ifdef HAVE_CDIPIO
  ! changes from I/O STAGE_DEFINITION to STAGE_CLEANUP
  IF (nprocio > 0 ) THEN
    CALL pioEndTimestepping
  ENDIF
#endif
  ! this call comes from stepon, just after integration_loop
  CALL close_output_streams

  ! stop hdiff diagnostics, if requested

  IF (ldiahdf) CALL finalize_hdiff_diag

  ! clean up

  CALL resetcolumn
  IF (iadvec == tpcore) CALL cleanup_tpcore
  IF (iadvec == semi_lagrangian) CALL cleanup_semi_lagrangian

  CALL cleanup_timer

  CALL free_memory
  CALL cleanup_scanbuffer
  CALL cleanup_greenhouse_gases
  CALL cleanup_o3clim
  CALL cleanup_clim
  CALL cleanup_sst

  ! cleanup of CDI resources
  CALL cleanup_output

  CALL cleanup_legendre
  CALL cleanup_solar_irradiance
  CALL cleanup_tmp_buffer
  CALL cleanup_gaussgrid
#ifdef FFT991
  CALL cleanup_fft991
#elif MKL_DFT
  CALL cleanup_mkl_dft
#elif ESSL_DFT
  CALL cleanup_ESSL_dft
#else
  CALL cleanup_fft992
#endif
  CALL cleanup_decomposition
  CALL dealloc_mods
  CALL cleanup_netcdf
  CALL cleanup_io
  CALL cleanup_geoloc
  CALL cleanup_truncation
  
  SELECT CASE (iaero)
  CASE (3)
    CALL cleanup_aero_kinne
  CASE (4)
    CALL cleanup_so4
  CASE (5) 
    CALL cleanup_aero_kinne
    CALL cleanup_aero_volc 
  CASE (6)
    CALL cleanup_aero_kinne
    CALL cleanup_aero_volc 
    CALL cleanup_aero_volc_tab_ham
  CASE (7)
    CALL cleanup_aero_kinne
    CALL cleanup_aero_volc_tab_crow
  END SELECT

  IF (lostation) THEN
    CALL cleanup_station_diag
  END IF

! set number of tracers to zero because the tracer memory is deallocated
  ntrac=0

  IF (lnudge) THEN
     CALL NudgingInit(NDG_CLEAN_MEM)
     IF (lstop) CALL NudgingInit(NDG_CLOSE)
  END IF

  IF (lnmi .AND. lstop) CALL NMI_Close

  WRITE (message_text,'(a)') &
       '======================================================================'
  CALL message('',message_text)
  WRITE (message_text,'(a)') &
       'Finished cleanup in control'
  CALL message('',message_text)
  WRITE (message_text,'(a)') &
       '======================================================================'
  CALL message('',message_text)

END SUBROUTINE control
