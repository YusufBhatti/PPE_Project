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
SUBROUTINE iorestart

  ! Description:
  !
  ! Reads netCDF history files for a resumed run.
  !
  ! Method:
  !
  ! *iorestart* positions data sets at the beginning of a rerun,
  ! writing data description records, and setting up necessary work
  ! files.
  !
  ! Information is written to the data description records of
  ! appropriate files, and work files are written if necessary.
  !
  !
  ! Authors:
  !
  ! L. Kornblueh, MPI, May 1999, f90 rewrite
  ! U. Schulzweida, MPI, May 1999, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_tracer_processes, ONLY: xt_initialize
  USE mo_exception,        ONLY: message
  USE mo_io,               ONLY: IO_read_streams
  USE mo_memory_gl,        ONLY: xt
  USE mo_memory_g1a,       ONLY: xtm1
  USE mo_control,          ONLY: ltimer
  USE mo_timer,            ONLY: timer_start, timer_stop, timer_restart
  USE mo_jsbach_interface, ONLY: jsbach_restart

  IMPLICIT NONE


  !  Executable statements 

  ! Restart from history files

  ! Read all restart files

  IF (ltimer) CALL timer_start(timer_restart)
  CALL IO_read_streams
  IF (ltimer) CALL timer_stop(timer_restart)

  CALL xt_initialize (xt, xtm1)

  CALL jsbach_restart

  CALL message('','')

END SUBROUTINE iorestart
