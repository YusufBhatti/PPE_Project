MODULE mo_ham_stream_hifreq

  ! *mo_ham_stream_hifreq* used for high frequency output of limited number of parameters
  !  for CCN/AOD study
  !
  ! Author:
  ! -------
  ! Philip Stier, University of Oxford, 2011-2014
  !

  USE mo_linked_list,   ONLY: t_stream
  USE mo_memory_base,   ONLY: new_stream, add_stream_element,                    &
                              default_stream_setting, add_stream_reference, AUTO
  USE mo_time_event,    ONLY: io_time_event
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_stream_hifreq ! construct the diag stream

  !--- Declarations for stream hifreq: --------------------------------------------------------------------------

  TYPE (t_stream), PUBLIC, POINTER :: hifreq

!----------------------------------------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE construct_stream_hifreq

    IMPLICIT NONE

    !--- Create new stream and add stream references for diagnostics:

    CALL new_stream (hifreq ,'hifreq', interval=io_time_event(3,'hours','first',0))

    CALL add_stream_reference (hifreq, 'aps'               ,'g3b'    ,lpost=.TRUE.)    
    CALL add_stream_reference (hifreq, 'gboxarea'          ,'geoloc' ,lpost=.TRUE.)

    CALL add_stream_reference (hifreq, 'TAU_2D_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'ABS_2D_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'ANG_550nm_865nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_DRY_2D_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'ABS_DRY_2D_550nm','rad',lpost=.TRUE.)

  END SUBROUTINE construct_stream_hifreq

END MODULE mo_ham_stream_hifreq
