!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_port_test

  ! L. Kornblueh, MPI, April 2003, implemented port test

  IMPLICIT NONE

  PRIVATE

  LOGICAL :: lport = .FALSE.

  PUBLIC :: lport
  PUBLIC :: init_port_test

CONTAINS

  SUBROUTINE init_port_test

    USE mo_memory_base, ONLY: new_stream, add_stream_reference
    USE mo_linked_list, ONLY: t_stream, NETCDF
    USE mo_time_event,  ONLY: io_time_event

    TYPE(t_stream), POINTER :: port

    ! This routine is called from subroutine 'init_memory', module
    ! 'mo_memory_streams'. Routines are called from here to allocate memory
    ! and to define output streams.

    IF (lport) THEN
    

      
      CALL new_stream (port, 'port', filetype=NETCDF, lpost=.TRUE., &
           lrerun=.FALSE., interval=io_time_event(1,'steps','first',0)) 
      
      CALL add_stream_reference (port, 'tm1', 'g1a', lpost=.TRUE., kprec=64) 
      CALL add_stream_reference (port, 'alpsm1', 'g1a', lpost=.TRUE., kprec=64) 
    END IF

  END SUBROUTINE init_port_test

END MODULE mo_port_test
