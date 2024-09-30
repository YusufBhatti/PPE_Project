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

    CALL new_stream (hifreq ,'hifreq', interval=io_time_event(3,'hours','first',0), lrerun=.FALSE.)

  !  CALL add_stream_reference (hifreq, 'aps'               ,'g3b'    ,lpost=.TRUE.)    
  !  CALL add_stream_reference (hifreq, 'gboxarea'          ,'geoloc' ,lpost=.TRUE.)

    CALL add_stream_reference (hifreq, 'TAU_2D_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'ABS_2D_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'ANG_550nm_865nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_DRY_2D_550nm','rad',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'ABS_DRY_2D_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_2D_865nm','rad',lpost=.TRUE.)

    CALL add_stream_reference (hifreq, 'REFFL_CT','activ',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'CDNC_BURDEN','activ',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'CDNC_CT','activ',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'ICNC_BURDEN','activ',lpost=.TRUE.)

    CALL add_stream_reference (hifreq, 'CN_BURDEN','ham',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'CCN_0.200','ham',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'CCN_0.500','ham',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'CCN_1.500','ham',lpost=.TRUE.)

    CALL add_stream_reference (hifreq, 'TAU_MODE_KS_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_MODE_CS_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_MODE_AS_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_MODE_KI_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_MODE_AI_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_MODE_CI_550nm','rad',lpost=.TRUE.)


    CALL add_stream_reference (hifreq, 'burden_DMS','burden',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'burden_SO2','burden',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'burden_H2SO4','burden',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'burden_SO4','burden',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'burden_BC','burden',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'burden_OC','burden',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'burden_SS','burden',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'burden_DU','burden',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'burden_WAT','burden',lpost=.TRUE.)



  END SUBROUTINE construct_stream_hifreq

END MODULE mo_ham_stream_hifreq
