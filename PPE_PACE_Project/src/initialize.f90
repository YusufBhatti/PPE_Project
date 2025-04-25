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
SUBROUTINE initialize

  ! Description:
  !
  ! Set up constants in various modules.
  !
  ! Method:
  !
  ! This subroutine initializes all the variables and arrays
  ! in modules.
  !
  ! *initialize* is called from *control*
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! R.G AND M.J, ECMWF, December 1982, changed
  ! U. Schlese, DKRZ, in 1994, and 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_io,                 ONLY: IO_init 
  USE mo_time_control,       ONLY: lfirst_cycle, init_manager, init_events, init_times
  USE mo_nudging_init,       ONLY: NudgingInit, NDG_INI_IO, NDG_INI_STREAM
  USE mo_greenhouse_gases,   ONLY: init_ghg
  USE mo_radiation_parameters,  ONLY: ighg, co2mmr
  USE mo_radiation,          ONLY: setup_radiation
  USE mo_column,             ONLY: setup_column
  USE m_alloc_mods,          ONLY: alloc_mods ! module subroutine
  USE mo_control,            ONLY: lcolumn, ldebugs, nprocio, &
                                   lyaxt_transposition, nlev
  USE mo_debugs,             ONLY: init_debugs
  USE mo_jsbach_interface,   ONLY: jsbach_init
  USE mo_co2,                ONLY: init_co2
  USE mo_cosp_simulator,     ONLY: cosp_initialize
  USE mo_station_diag,       ONLY: station_diag_nml
  USE mo_memory_cfdiag,      ONLY: setup_cfdiag
  USE mo_cosp_offline,       ONLY: cosp_offline_nml
!++mgs
  USE mo_submodel,           ONLY: setsubmodel, lanysubmodel
!--mgs
  USE mo_tracer,             ONLY: init_trlist, finish_tracer_definition
  USE mo_submodel_interface, ONLY: init_subm
!++mgs
!!  USE mo_species,            ONLY: init_splist 
!--mgs
  USE mo_exception,          ONLY: message_text, finish, number_of_errors
  USE mo_echam_yaxt,         ONLY: setup_yaxt_decomposition, generate_yaxt_redist, &
       & add_yaxt_gp_nlevs
  USE mo_mpi,                    ONLY: p_all_comm

  IMPLICIT NONE

  !  External subroutines 
  EXTERNAL :: inictl, setdyn, setphys, setgws

  !  Executable statements 

  !-- 1. Set control variables

  !-- 1.1 Set general control variables and time stepping
  !--     Set I/O units and buffer indices

  CALL inictl

  !-- 1.2 Initialize netCDF IO

  CALL IO_init

  !-- 1.3 Initialize column model

  CALL setup_column (lcolumn)

  CALL alloc_mods

  !-- 1.4 Initialize time manager

  IF (lfirst_cycle) CALL init_manager

  !-- 1.5 Initialize nudging if selected

  CALL NudgingInit(NDG_INI_IO)

  !-- 1.6 Derive all time management dates and times 

  CALL init_times

  !-- 2. Compute decomposition 

  CALL init_decomposition

  CALL NudgingInit(NDG_INI_STREAM)

  ! initialize column model. jsr: the following line is not 
  ! commented in e.g. in rev 1362

!  CALL setcolumn

  !-- 3. Preset, modify and derive values needed in the
  !      dynamics and the initialisation and call helmo the first time.

  CALL setdyn

  ! read submodel name list and register submodels
  
  CALL setsubmodel
  
  !-- 4. Preset, modify and derive values needed in the physics.

  CALL setphys

  !-- 5. Preset, modify and derive values needed in the radiation.

  CALL setup_radiation

  !-- 6. Preset, modify and derive values needed in the gwspectrum param.

  CALL setgws

  !-- 7. Prepare greenhouse gas scenario

  IF(ighg .NE. 0) CALL init_ghg(ighg)

  !-- 8. Final event evaluation

  CALL init_events

  !-- 9. initialize debug stream
  IF (ldebugs) THEN
     CALL init_debugs
  END IF

  !-- Initialize special diagnostics (CNam,mns)
  CALL cosp_initialize
  CALL station_diag_nml
  CALL setup_cfdiag
  CALL cosp_offline_nml

  ! first interface to initialize submodels
  
  CALL init_trlist
!++mgs
  IF (lanysubmodel) CALL init_subm  ! This does not call init_co2
!--mgs

  CALL init_co2(co2mmr)
 
  CALL finish_tracer_definition     ! Needs to be after init_co2!

  ! jsbach init
  CALL jsbach_init !after init_co2

  ! terminate model run if any errors were detected during the initialisation
  IF (number_of_errors > 0) THEN
    WRITE (message_text, *) number_of_errors, ' error(s) encountered during initialisation. Run aborted.'
    CALL finish('initialize', message_text)
  END IF

  IF (lyaxt_transposition .OR. nprocio > 0) THEN
    CALL add_yaxt_gp_nlevs((/1, nlev, nlev+1/))
    CALL setup_yaxt_decomposition
  ENDIF
  IF (lyaxt_transposition) THEN
    CALL generate_yaxt_redist(p_all_comm)
  ENDIF

END SUBROUTINE initialize

