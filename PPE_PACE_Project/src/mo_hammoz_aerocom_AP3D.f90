!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_AP3D.f90
!!
!! \brief
!! Module for AeroCom Phase 3 experiments daily (AP3D) diagnostics
!!
!! \author D. Neubauer (ETH Zurich)
!!  adapted from AerChemMIP diagnostics from S. Ferrachat
!!
!! D. Neubauer, david.neubauer@env.ethz.ch
!!
!! \revision_history
!!   -# D. Neubauer (ETH Zurich) - original code (2019-03-12)
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
MODULE mo_hammoz_aerocom_AP3D

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  
  IMPLICIT NONE  

  PRIVATE

  PUBLIC :: construct_AP3D_stream
  PUBLIC :: update_AP3D_diags

  TYPE (t_stream), PUBLIC, POINTER :: acp3d
 
  REAL(dp), PUBLIC, POINTER :: rsuscs(:,:)  => NULL()
  REAL(dp), PUBLIC, POINTER :: rsnscs(:,:)  => NULL()
  REAL(dp), PUBLIC, POINTER :: rlus(:,:)  => NULL()
  REAL(dp), PUBLIC, POINTER :: rlds(:,:)  => NULL()
  REAL(dp), PUBLIC, POINTER :: ec550aer_int(:,:)  => NULL()
  REAL(dp), PUBLIC, POINTER :: abs550aer_int(:,:)  => NULL()
  REAL(dp), PUBLIC, POINTER :: ec550aer(:,:)  => NULL()
  REAL(dp), PUBLIC, POINTER :: abs550aer(:,:)  => NULL()

CONTAINS

  SUBROUTINE construct_AP3D_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE 
    USE mo_memory_base,         ONLY: new_stream, add_stream_element, AUTO,  &
                                    default_stream_setting, add_stream_reference
    USE mo_time_event,          ONLY: io_time_event
    USE mo_ham,                 ONLY: nclass
    USE mo_exception,           ONLY: message, em_error

    TYPE(io_time_event) :: put_interval
    INTEGER :: ipm
    
    !-- set output interval
    put_interval%counter      = 1
    put_interval%unit         = 'days'
    put_interval%adjustment   = 'last'
    put_interval%offset       = 0
    
    !-- Create new stream:
    CALL new_stream (acp3d ,'acp3d', &
                     filetype = AEROCOM_FILETYPE, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)
  
    !-- Add standard fields for post-processing:
    CALL default_stream_setting (acp3d, &
                                 table = 199, &
                                 code = AUTO )

    !-- Basic vars (may sometimes be necessary)
    CALL add_stream_reference (acp3d, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (acp3d, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (acp3d, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (acp3d, 'gboxarea','geoloc',lpost=.TRUE.)

    CALL add_stream_element (acp3d, 'rsuscs', rsuscs, &
        longname = 'Surface Upwelling Clear-Sky Shortwave Radiation', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3d, 'rsnscs', rsnscs, &
        longname = 'Surface Net Clear-Sky Shortwave Radiation', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3d, 'rlus', rlus, &
        longname = 'Surface Upwelling Longwave Radiation', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acp3d, 'rlds', rlds, &
        longname = 'Surface Downwelling Longwave Radiation', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acp3d, 'ec550aer', ec550aer, &
        longname = 'Surface ambient aerosol extinction at 550 nm', &
        units = 'm-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acp3d, 'abs550aer', abs550aer, &
        longname = 'Surface ambient aerosol absorption at 550 nm', &
        units = 'm-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acp3d, 'ec550aer_int', ec550aer_int, &
        longname = 'Surface ambient aerosol integrated extinction at 550 nm', &
        units = 'm-1', &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )
    
    CALL add_stream_element (acp3d, 'abs550aer_int', abs550aer_int, &
        longname = 'Surface ambient aerosol integrated absorption at 550 nm', &
        units = 'm-1', &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )
        
  END SUBROUTINE construct_AP3D_stream

  SUBROUTINE update_AP3D_diags(kproma, kbdim, klev, krow)

    USE mo_time_control, ONLY: delta_time
    USE mo_memory_cfdiag,ONLY: irlu, irld
    USE mo_vphysc,       ONLY: vphysc

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow
    
    !-- rlus,rlds
    rlus(1:kproma,krow) = rlus(1:kproma,krow) + delta_time*irlu(1:kproma,klev+1,krow)
    rlds(1:kproma,krow) = rlds(1:kproma,krow) + delta_time*irld(1:kproma,klev+1,krow)
          
    !-- ec550aer, abs550aer, ec550dryaer, abs550dryaer
    ec550aer(1:kproma,krow)  = ec550aer(1:kproma,krow)  + &
         ec550aer_int(1:kproma,krow)/vphysc%grheightm1(1:kproma,klev,krow)*delta_time
    abs550aer(1:kproma,krow) = abs550aer(1:kproma,krow) + &
         abs550aer_int(1:kproma,krow)/vphysc%grheightm1(1:kproma,klev,krow)*delta_time

  END SUBROUTINE update_AP3D_diags

END MODULE mo_hammoz_aerocom_AP3D
