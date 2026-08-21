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

    CALL new_stream (hifreq ,'hifreq', interval=io_time_event(6,'hours','first',0), lrerun=.FALSE.)


    CALL add_stream_reference (hifreq, 'TAU_2D_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'OMEGA_2D_440nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'ABS_2D_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'ANG_550nm_865nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'ANG_440nm_670nm' ,'rad',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'TAU_3D_355nm','rad',lpost=.TRUE.)
!
!    CALL add_stream_reference (hifreq, 'REFFL_CT','activ',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'CDNC_BURDEN','activ',lpost=.TRUE.)
!!    CALL add_stream_reference (hifreq, 'CDNC_CT','activ',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'CDNC_INCL_CT','activ',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'ICNC_BURDEN','activ',lpost=.TRUE.)
!
!    CALL add_stream_reference (hifreq, 'CN_BURDEN','ham',lpost=.TRUE.) ! zrmin minimum radius changed to 1.0e-7 for PACE comparison
!    CALL add_stream_reference (hifreq, 'CCN_0.300','ham',lpost=.TRUE.)
!!    CALL add_stream_reference (hifreq, 'CCN_0.500','ham',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'CCN_1.500','ham',lpost=.TRUE.)
!
    CALL add_stream_reference (hifreq, 'TAU_2D_MODE_KS_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_2D_MODE_CS_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_2D_MODE_AS_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_2D_MODE_KI_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_2D_MODE_AI_550nm','rad',lpost=.TRUE.)
    CALL add_stream_reference (hifreq, 'TAU_2D_MODE_CI_550nm','rad',lpost=.TRUE.)
!

    CALL add_stream_reference (hifreq, 'conccnmodeNS','achemon',lpost=.TRUE.) ! Nuc Sol 
    CALL add_stream_reference (hifreq, 'conccnmodeKS','achemon',lpost=.TRUE.) ! Aitken Sol 
    CALL add_stream_reference (hifreq, 'conccnmodeAS','achemon',lpost=.TRUE.) !  Accum Sol
    CALL add_stream_reference (hifreq, 'conccnmodeCS','achemon',lpost=.TRUE.) ! 
!
!    CALL add_stream_reference (hifreq, 'ccn1_inst','acheaci',lpost=.TRUE.) ! CCN1 
!    CALL add_stream_reference (hifreq, 'ccn3_inst','acheaci',lpost=.TRUE.) ! CCN3
    CALL add_stream_reference (hifreq, 'dpg','acheaci',lpost=.TRUE.) ! dpg to intigrate ccn with
    CALL add_stream_reference (hifreq, 't3d','acheaci',lpost=.TRUE.) ! t3d to intigrate ccn with
    CALL add_stream_reference (hifreq, 'aps','acheaci',lpost=.TRUE.) ! aps to intigrate ccn with
    CALL add_stream_reference (hifreq, 'hyam','acheaci',lpost=.TRUE.) ! aps to intigrate ccn with


    CALL add_stream_reference (hifreq, 'rhoam1','vphysc',lpost=.TRUE.) ! aps to intigrate ccn with
    CALL add_stream_reference (hifreq, 'rhoam1_moist','vphysc',lpost=.TRUE.) ! aps to intigrate ccn with
!!    CALL add_stream_reference (hifreq, 'burden_DMS','burden',lpost=.TRUE.)
!!    CALL add_stream_reference (hifreq, 'burden_SO2','burden',lpost=.TRUE.)
!!    CALL add_stream_reference (hifreq, 'burden_H2SO4','burden',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'burden_SO4','burden',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'burden_BC','burden',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'burden_OC','burden',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'burden_SS','burden',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'burden_DU','burden',lpost=.TRUE.)
!    CALL add_stream_reference (hifreq, 'burden_WAT','burden',lpost=.TRUE.)
!
!    CALL add_stream_reference (hifreq, 'xlvi','g3b',lpost=.TRUE.) ! LWP
!    CALL add_stream_reference (hifreq, 'aps','g3b',lpost=.TRUE.) !Suface PressureP
!!    CALL add_stream_reference (hifreq, 'water_volume_fraction','ham',lpost=.TRUE.) ! water volume fraction
!    CALL add_stream_reference (hifreq, 'lcc','achemon',lpost=.TRUE.) ! liquid Cloud Fraction 2D
!    CALL add_stream_reference (hifreq, 'iccl','achemon',lpost=.TRUE.) ! ice Cloud Fraction 3D
!    CALL add_stream_reference (hifreq, 'aclcov','g3b',lpost=.TRUE.) ! Cloud Fraction 2D

  END SUBROUTINE construct_stream_hifreq

END MODULE mo_ham_stream_hifreq
