#ifdef __xlC__
@PROCESS STRICT
#endif

PROGRAM master

  !----------------------------------------------------------------------------
  !
  ! Copyright 2014 by Max Planck Institute for Meteorology
  !
  ! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
  ! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
  ! file COPYING in the root of the source tree for this code.
  ! Where software is supplied by third parties, it is indicated in the headers of the routines.
  !
  ! The ECHAM6 web page for developers (login required) is at:
  !
  ! https://code.zmaw.de/projects/echam
  !
  !----------------------------------------------------------------------------
  !
  ! Call the control subroutine (*control*).
  !
  ! Externals:
  !
  ! *control*   called to control the run.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 2000, date/time control
  ! S. Legutke, MPI M&D, July 2001, redirect stdout for coupling
  ! L. Kornblueh, MPI, April 2006, OpenMP explicit environment setting
  ! L. Kornblueh, MPI, July 2008, Profiling
  ! L. Kornblueh, MPI, July 2009, parallel I/O

  USE mo_kind,           ONLY: dp, sp
  USE mo_mpi,            ONLY: p_start, p_stop, p_bcast,                       &
                               p_global_comm, p_bcast_int_1d,                  &
                               p_all_comm, p_init_model_communicator,          &
#ifdef _PROFILE
                               p_pe,                                           &
#endif
                               p_parallel, p_parallel_io, p_io, p_comm_null
#ifdef _OPENMP
  USE omp_lib,           ONLY: omp_set_dynamic, omp_set_num_threads
#endif
#ifdef _PROFILE
  USE mo_profile,        ONLY: trace_init, trace_start, trace_stop, &
                               trace_finalize
#endif
  USE mo_exception,      ONLY: finish, message, message_text, &
                               open_log, close_log
  USE mo_time_control,   ONLY: lbreak, lstop, lresume_cmd
  USE mo_control,        ONLY: nproca, nprocb, nprocio, lyaxt_transposition
  USE mo_ensemble,       ONLY: init_ensemble, finalize_ensemble
  USE mo_namelist,       ONLY: open_nml, position_nml, positioned, close_nml
  USE mo_util_logging,   ONLY: open_network_log, close_network_log
  USE mo_util_db_timings,ONLY: set_db_connection
  USE mo_util_vcs,       ONLY: util_repository_url, &
       &                       util_branch_name,    &
       &                       util_revision_key

  USE mo_echam_yaxt,     ONLY: yaxt_initialize, yaxt_finalize

#ifdef __prism
  USE mo_couple_wrap,    ONLY: oasis_ok, oasis_terminate
  USE mo_couple,         ONLY: couple_end
  USE mo_timer,          ONLY: timer_stop, timer_start, timer_couple_end
#endif
  

  IMPLICIT NONE

#ifdef HAVE_CDIPIO
  INCLUDE 'cdipio.inc'
#endif

  !  External functions 
  REAL(dp), EXTERNAL :: util_walltime
  INTEGER,  EXTERNAL :: util_cputime
  EXTERNAL           :: util_set_malloc_alignment   

  !  External subroutines 
  EXTERNAL :: control
#ifdef __XT3__
  ! routine to set io buffers size on Catamount
  EXTERNAL :: util_base_iobuf
#endif

  REAL(dp) :: zwinit, zwtime, zutime, zstime, zrtime
  REAL(sp) :: partInFlate

  LOGICAL :: lbuffer(3)
  INTEGER :: ibuffer(4), iret, inml, iunit, iomode

  CHARACTER(len=132) :: network_logger = ''
  CHARACTER(len=132) :: db_host        = ''

  CHARACTER(len=256) :: repository = ''
  CHARACTER(len=256) :: branch     = ''
  CHARACTER(len=256) :: revision   = ''
  
  INTEGER :: nlen

  INTEGER :: i
  CHARACTER(len=32) :: arg

  ! namelist definition for parctl (parallel decomposition)
  INCLUDE 'parctl.inc'

  ! For parallel I/O preset numbers

  nprocio = 0
  iomode  = 0

  ! Initialize wallclock timer
  ! zwinit is an unused variable, but necessary to initialize wallclock timer

  zwinit = util_walltime(0)
  
  ! this aligns allocates on L1 data cache lines on AIX

  CALL util_set_malloc_alignment()

#ifdef __XT3__
  ! set buffer size for stderr and stdout on Catamount
  CALL util_base_iobuf
#endif

  ! start MPI, do not use any OpenMP before p_start

  CALL p_start

  ! Print version

  nlen = 256
  call util_repository_url(repository, nlen)
  nlen = 256
  call util_branch_name(branch, nlen)
  nlen = 256
  call util_revision_key(revision, nlen)
    
  CALL message ('','')
  CALL message ('','')
  message_text = '==========================================================='
  CALL message ('',TRIM(message_text))        
  CALL message ('','')
  message_text = '  ECHAM - Release 6.3.02 '
  CALL message ('',TRIM(message_text))        
  message_text = '  Copyright by Max-Planck-Institute for Meteorology, 2015'
  CALL message ('',TRIM(message_text))        
  message_text = '  Read echam6.f90 and mpi-m_sla_201202.pdf before using ECHAM6'
  CALL message ('',TRIM(message_text))        
  CALL message('','')
  WRITE(message_text,'(a,a)') 'Repository: ', TRIM(repository)
  CALL message('',message_text)
  WRITE(message_text,'(a,a)') 'Branch    : ', TRIM(branch)
  CALL message('',message_text)
  WRITE(message_text,'(a,a)') 'Revision  : ', TRIM(revision)
  CALL message('',message_text)
#ifdef HAMMOZ
  CALL message ('','')
  message_text = '  HAMMOZ - Release candidate for: HAM2.3-MOZ1.0 '
  CALL message ('',TRIM(message_text))        
  message_text = '  Copyright by the ECHAM-HAMMOZ consortium, 2017'
  CALL message ('',TRIM(message_text))        
  message_text = '  Project info: https://redmine.hammoz.ethz.ch/projects/hammoz'
  CALL message ('',TRIM(message_text))        
  message_text = '  Read doc/ECHAM-HAMMOZ_namelists.pdf for details on the HAMMOZ namelists.'
  CALL message ('',TRIM(message_text))        
#endif
   
  CALL message ('','')
  message_text = '==========================================================='
  CALL message ('',TRIM(message_text))        
  CALL message ('','')
  CALL message ('','')
   
  ! read command line argument

  DO i = 1, COMMAND_ARGUMENT_COUNT()
    CALL GET_COMMAND_ARGUMENT(i, arg)
    SELECT CASE (arg)
    CASE ('-i', '--init')
      lresume_cmd = .TRUE.
    CASE DEFAULT
      CALL finish('','Unrecognized command-line option: '//TRIM(arg))
    END SELECT
  END DO

  IF (p_parallel) THEN
    CALL p_bcast (lresume_cmd, p_io, comm=p_global_comm)
  END IF

  ! start profiling

#ifdef _PROFILE
  CALL trace_init ('echam6', p_pe)
  CALL trace_start ('all', 0)
#endif
  
  ! check, if an ensemble is setup
  CALL init_ensemble(finish, p_bcast_int_1d, p_io)

  ! read processor decomposition

  IF (p_parallel_io) THEN
    inml = open_nml ('namelist.echam')
    iunit = position_nml ('PARCTL', inml, status=iret)
    SELECT CASE (iret)
    CASE (POSITIONED)
      READ (iunit, parctl)
    END SELECT
    CALL close_nml(inml)
  ENDIF

  IF (p_parallel) THEN
    ibuffer(:) = (/ nproca, nprocb, nprocio, iomode /)
    CALL p_bcast (ibuffer, p_io, comm=p_global_comm)
    nproca  = ibuffer(1)
    nprocb  = ibuffer(2)
    nprocio = ibuffer(3)
    iomode  = ibuffer(4)
    CALL p_bcast (network_logger, p_io, comm=p_global_comm)
    CALL p_bcast (db_host,        p_io, comm=p_global_comm)
    CALL p_bcast (lyaxt_transposition, p_io, comm=p_global_comm)
  END IF

  IF (nprocio > 0 .OR. lyaxt_transposition) CALL yaxt_initialize(p_global_comm)

  partInFlate = 1.1_sp   ! ceil (nlon * nlat/(nproca * nprocb))/
                         ! floor(nlon * nlat/(nproca * nprocb))
  ! create model PEs communicator p_all_comm
  CALL p_init_model_communicator(p_global_comm, nproca, nprocb, nprocio, &
                                 iomode, PartInFlate)
#ifdef HAVE_CDIPIO
  ! finish OASIS, YAXT and MPI on IO-servers
  IF (p_all_comm == p_comm_null) THEN
#ifdef __prism
    CALL oasis_terminate (iret)
    IF (iret /= oasis_ok) THEN
      CALL finish('master: oasis_terminate failed!')
    ENDIF
#endif
    CALL yaxt_finalize()
    CALL p_stop
    STOP
  ENDIF
#endif  /* HAVE_CDIPIO */

  ! start other services

  IF (network_logger /= '') THEN
    CALL open_network_log(TRIM(network_logger), 19876)
  ENDIF

  IF (db_host /= '') THEN
    CALL set_db_connection(db_host, 'echam-timings', 'echam')
  ENDIF

#if defined (__prism) 
    !--  Redirect standard output file to atmout if coupled run.
    CALL open_log('atmout')
    WRITE(message_text,*)'Atmosphere standard output is assigned to file atmout.'
#else
    CALL open_log('echam6.log')
    WRITE(message_text,*)'Atmosphere standard output is assigned to file echam6.log.'
#endif
    
    CALL message ('',TRIM(message_text))          
    CALL message ('','')
    
    ! Loop over rerun cycles
    
    rerun_cycles: DO 
      
      ! do the model run
      
      CALL control
      
      ! model run ended, maybe another cycle will come 
      
      IF (lbreak .OR. lstop) EXIT
      
      CALL message('','Start next rerun cycle.')
      
    END DO rerun_cycles

#ifdef HAVE_CDIPIO
  IF (nprocio > 0 ) THEN
    ! Finalize parallel I/O with CDI
    CALL pioFinalize()
  ENDIF
#endif

#ifdef __prism
    CALL timer_start(timer_couple_end)
    CALL couple_end
    CALL timer_stop(timer_couple_end)
#endif
    
    ! stop profiling
    
#ifdef _PROFILE
    CALL trace_stop ('all', 0)
    CALL trace_finalize (p_pe)
#endif
  
    ! wallclock timer output

    iret = util_cputime(zutime, zstime)
    IF (iret == -1) THEN
      CALL message('','Cannot determine used CPU time')
    ELSE
      zwtime = util_walltime(0)
      zrtime = (zutime+zstime)/zwtime
      CALL message ('', '')
      WRITE (message_text,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
      CALL message('',message_text)
      WRITE (message_text,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
      CALL message('',message_text)
      WRITE (message_text,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
      CALL message('',message_text)
      WRITE (message_text,'(a,f10.2,a)') ' Ratio            : ', 100*zrtime, ' %'
      CALL message('',message_text)     
      CALL message ('', '')
    END IF

  CALL close_log

  CALL close_network_log

  CALL finalize_ensemble

  ! broadcast communicator is model PEs communicator
  CALL p_bcast(lstop, p_io, comm=p_all_comm)


  IF (lyaxt_transposition .OR. nprocio > 0) CALL yaxt_finalize()

  ! stop MPI
  
  CALL p_stop                      
  
  ! bail out
  
  IF (lstop) THEN
    ! POE synchronizes only if exit return value < 128 and > 1 or 0 ...
    CALL finish('','Experiment finished.', 127)
  ELSE
    CALL message('','Experiment checkpointed.')
  END IF

END PROGRAM master
