!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_HEpro.f90
!!
!! \brief
!! Module for AeroCom Holuhraun Process diagnostics
!!
!! \author D. Neubauer (ETH Zurich)
!!  adapted from AerChemMIP diagnostics from S. Ferrachat
!!
!! D. Neubauer, david.neubauer@env.ethz.ch
!!
!! \revision_history
!!   -# D. Neubauer (ETH Zurich) - original code (2018-08-14)
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_hammoz_aerocom_HEpro

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  
  IMPLICIT NONE  

  PRIVATE

  PUBLIC :: construct_HEpro_stream
  PUBLIC :: update_HEpro_diags

  TYPE (t_stream), PUBLIC, POINTER :: achepro
  
  !SFNote: the variables in this group are output as zonal means.
  !        The zonal mean calculations can't be handled by echam,
  !        and will be done at p-proc.
  !        This means that in this stream the zonal mean aspect
  !        is simply neglected.

  CONTAINS

  !------------------------------------------------

  SUBROUTINE construct_HEpro_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE 
    USE mo_memory_base,       ONLY: new_stream, add_stream_element, AUTO,  &
                                    default_stream_setting, add_stream_reference
    USE mo_time_event,        ONLY: io_time_event

    TYPE(io_time_event) :: put_interval

    !-- set output interval
    put_interval%counter      = 6
    put_interval%unit         = 'hours'
    put_interval%adjustment   = 'last'
    put_interval%offset       = 0
    
    !-- Create new stream:
    CALL new_stream (achepro ,'achepro', &
                     filetype = AEROCOM_FILETYPE, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)
  
    !-- Add standard fields for post-processing:
    CALL default_stream_setting (achepro, &
                                 table = 199, &
                                 code = AUTO )

    !-- Basic vars (may sometimes be necessary)
    CALL add_stream_reference (achepro, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (achepro, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (achepro, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (achepro, 'gboxarea','geoloc',lpost=.TRUE.)
  
    !-- Diagnostics table'

  END SUBROUTINE construct_HEpro_stream

  SUBROUTINE update_HEpro_diags(kproma, kbdim, klev, krow)

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow

  END SUBROUTINE update_HEpro_diags

END MODULE mo_hammoz_aerocom_HEpro
