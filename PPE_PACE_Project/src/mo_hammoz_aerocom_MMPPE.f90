!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_MMPPE.f90
!!
!! \brief
!! Module for AeroCom Multi Model Perturbed Physics Ensemble (MMPPE) diagnostics
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
MODULE mo_hammoz_aerocom_MMPPE

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  
  IMPLICIT NONE  

  PRIVATE

  PUBLIC :: construct_MMPPE_stream
  PUBLIC :: update_MMPPE_diags

  TYPE (t_stream), PUBLIC, POINTER :: acmmppe
  
  REAL(dp), PUBLIC, POINTER :: n3(:,:,:)
  REAL(dp), PUBLIC, POINTER :: n50(:,:,:)
  REAL(dp), PUBLIC, POINTER :: fliq3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: fliq2d(:,:)
  REAL(dp), PUBLIC, POINTER :: srain_inst(:,:)

  CONTAINS

  !------------------------------------------------

  SUBROUTINE construct_MMPPE_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE 
    USE mo_memory_base,       ONLY: new_stream, add_stream_element, AUTO,  &
                                    default_stream_setting, add_stream_reference
    USE mo_time_event,        ONLY: io_time_event

    TYPE(io_time_event) :: put_interval

    !-- set output interval
    put_interval%counter      = 1
    put_interval%unit         = 'days'
    put_interval%adjustment   = 'last'
    put_interval%offset       = 0
    
    !-- Create new stream:
    CALL new_stream (acmmppe ,'acmmppe', &
                     filetype = AEROCOM_FILETYPE, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)
  
    !-- Add standard fields for post-processing:
    CALL default_stream_setting (acmmppe, &
                                 table = 199, &
                                 code = AUTO )

    !-- Basic vars (may sometimes be necessary)
    CALL add_stream_reference (acmmppe, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (acmmppe, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (acmmppe, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (acmmppe, 'gboxarea','geoloc',lpost=.TRUE.)
  
    CALL add_stream_element (acmmppe,'N3',n3, &
         longname='number_concentration_for_larger_3_nm', &
         units='cm-3',   &
         laccu = .FALSE.,&
         lpost = .TRUE., &
         lrerun = .TRUE. )

    CALL add_stream_element (acmmppe,'N50',n50, &
         longname='number_concentration_for_larger_50_nm', &
         units='cm-3',   &
         laccu = .FALSE.,&
         lpost = .TRUE., &
         lrerun = .TRUE. )

    CALL add_stream_element (acmmppe, 'fliq3d', fliq3d, &
         longname = 'liquid_phase_cloud_fraction', &
         units = '1', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )
    
    CALL add_stream_element (acmmppe, 'fliq2d', fliq2d, &
         longname = 'liquid_phase_cloud_area_fraction', &
         units = '1', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )

    CALL add_stream_element (acmmppe, 'srain', srain_inst, &
         longname = 'stratiform_rain_rate', &
         units = 'km m-2 s-1', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )
    
   !-- Diagnostics table'

  END SUBROUTINE construct_MMPPE_stream

  SUBROUTINE update_MMPPE_diags(kproma, kbdim, klev, krow)

    USE mo_ham,            ONLY: nclass, sizeclass
    USE mo_ham_streams,    ONLY: rdry
    USE mo_ham_tools,      ONLY: ham_m7_logtail
    USE mo_vphysc,         ONLY: vphysc
    USE mo_memory_g1a,     ONLY: xtm1
    USE mo_hammoz_aerocom_HEaci, ONLY: f3d, phase3d

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow

    INTEGER :: jclass,it

    REAL(dp) :: zfracn3(kbdim,klev,nclass),zfracn50(kbdim,klev,nclass)
    REAL(dp) :: zrdry(kbdim,klev,nclass),rcrit(kbdim,klev)
    REAL(dp) :: cfracn(nclass)
    
    cfracn(:) = (/1.0_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp/)

!-- Aerosol number    
    n3(1:kproma,:,krow) = 0._dp
    n50(1:kproma,:,krow) = 0._dp
    DO jclass=1, nclass
       it = sizeclass(jclass)%idt_no
       zrdry(1:kproma,:,jclass)=rdry(jclass)%ptr(1:kproma,:,krow)
       rcrit(1:kproma,:)=1.5e-9_dp!=3E-9 dry diameter
       CALL ham_m7_logtail(kproma,    kbdim,  klev,   krow, jclass, &
                           .TRUE.,   zrdry(:,:,jclass),            &
                           rcrit, zfracn3(:,:,jclass))
       rcrit(1:kproma,:)=25.e-9_dp!=50E-9 dry diameter
       CALL ham_m7_logtail(kproma,    kbdim,  klev,   krow, jclass, &
                           .TRUE.,   zrdry(:,:,jclass),            &
                           rcrit, zfracn50(:,:,jclass))
       n3(1:kproma,:,krow) = n3(1:kproma,:,krow)                   &
                           + xtm1(1:kproma,:,it,krow)*vphysc%rhoam1(1:kproma,:,krow) &
                             *zfracn3(1:kproma,:,jclass)*cfracn(jclass)*1.0e-6_dp
       n50(1:kproma,:,krow) = n50(1:kproma,:,krow)                   &
                           + xtm1(1:kproma,:,it,krow)*vphysc%rhoam1(1:kproma,:,krow) &
                             *zfracn50(1:kproma,:,jclass)*cfracn(jclass)*1.0e-6_dp
    END DO
    
!-- fliq3d
    fliq3d(1:kproma,:,krow) = f3d(1:kproma,:,krow) * phase3d(1:kproma,:,krow)
    
  END SUBROUTINE update_MMPPE_diags

END MODULE mo_hammoz_aerocom_MMPPE
