!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_MMPPEmon.f90
!!
!! \brief
!! Module for AeroCom Multi Model Perturbed Physics Ensemble (MMPPE) monthly diagnostics
!!
!! \author D. Neubauer (ETH Zurich)
!!  adapted from AerChemMIP diagnostics from S. Ferrachat
!!
!! D. Neubauer, david.neubauer@env.ethz.ch
!!
!! \revision_history
!!   -# D. Neubauer (ETH Zurich) - original code (2019-03-12)
!!   -# A. Arifi - monthly version of the MMPPE diagnostics (2026-08-21)
!!
!! \limitations
!! None
!!
!! \details
!! Monthly counterpart of the MMPPE stream (mo_hammoz_aerocom_MMPPE): the
!! diagnostics are meant to be accumulated over the output interval
!! (laccu = .TRUE.), so that the post-processed fields are monthly means rather
!! than instantaneous snapshots.
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
MODULE mo_hammoz_aerocom_MMPPEmon

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_MMPPEmon_stream
  PUBLIC :: update_MMPPEmon_diags

  TYPE (t_stream), PUBLIC, POINTER :: acmmppemon

  CONTAINS

  !------------------------------------------------

  SUBROUTINE construct_MMPPEmon_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE
    USE mo_memory_base,       ONLY: new_stream, add_stream_element, AUTO,  &
                                    default_stream_setting, add_stream_reference
    USE mo_time_event,        ONLY: io_time_event

    TYPE(io_time_event) :: put_interval

    !-- set output interval
    put_interval%counter      = 1
    put_interval%unit         = 'months'
    put_interval%adjustment   = 'last'
    put_interval%offset       = 0

    !-- Create new stream:
    CALL new_stream (acmmppemon ,'acmmppemon', &
                     filetype = AEROCOM_FILETYPE, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)

    !-- Add standard fields for post-processing:
    CALL default_stream_setting (acmmppemon, &
                                 table = 199, &
                                 code = AUTO )

    !-- Basic vars (may sometimes be necessary)
    CALL add_stream_reference (acmmppemon, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (acmmppemon, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (acmmppemon, 'aps'     ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (acmmppemon, 'gboxarea','geoloc',lpost=.TRUE.)

    !-- Diagnostics table'

  END SUBROUTINE construct_MMPPEmon_stream

  SUBROUTINE update_MMPPEmon_diags(kproma, kbdim, klev, krow)

    USE mo_time_control,   ONLY: delta_time

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow

  END SUBROUTINE update_MMPPEmon_diags

END MODULE mo_hammoz_aerocom_MMPPEmon
