!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_HEaci.f90
!!
!! \brief
!! Module for AeroCom Holuhraun ACI diagnostics
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
MODULE mo_hammoz_aerocom_HEaci

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  
  IMPLICIT NONE  

  PRIVATE

  PUBLIC :: construct_HEaci_stream
  PUBLIC :: init_HEaci
  PUBLIC :: update_HEaci_diags

  TYPE (t_stream), PUBLIC, POINTER :: acheaci
  
  INTEGER, PUBLIC, ALLOCATABLE :: nstr_cld_top(:,:)
  INTEGER, PUBLIC, ALLOCATABLE :: nstr_cld_liq_top(:,:)

  REAL(dp), PUBLIC, POINTER :: emiso2(:,:)
  REAL(dp), PUBLIC, POINTER :: emiso2_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn01bl(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn03bl(:,:)
  REAL(dp), PUBLIC, POINTER :: od550aer(:,:)
  REAL(dp), PUBLIC, POINTER :: od440aer(:,:)
  REAL(dp), PUBLIC, POINTER :: od865aer(:,:)
  REAL(dp), PUBLIC, POINTER :: aerindex(:,:) 
  REAL(dp), PUBLIC, POINTER :: aerindex_inst(:,:) 
  REAL(dp), PUBLIC, POINTER :: precip_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: sprecip_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: ttop_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: tmf_cld_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: cdnc_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: cdr_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: lwp(:,:)
  REAL(dp), PUBLIC, POINTER :: cod_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: clt_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: accretn(:,:)
  REAL(dp), PUBLIC, POINTER :: autoconv(:,:)
  REAL(dp), PUBLIC, POINTER :: rsut_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: rsutcs_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: rsutnoa_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: rsutcsnoa_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: rlutnoa_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: rlutcsnoa_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: rlut_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: rlutcs_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: rsut_tmp(:,:)
  REAL(dp), PUBLIC, POINTER :: rsutcs_tmp(:,:)
  REAL(dp), PUBLIC, POINTER :: rsutnoa_tmp(:,:)
  REAL(dp), PUBLIC, POINTER :: rsutcsnoa_tmp(:,:)
  REAL(dp), PUBLIC, POINTER :: f3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: fice3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cod3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cod3dswl(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cod3dswi(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cdnc3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cdr3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: phase3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: lts(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn1_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn3_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn5_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn10_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: z_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: accret3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: auto3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: dpg(:,:,:)
  REAL(dp), PUBLIC, POINTER :: riming3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: qcsedten3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: qcevap3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: od550aer3d(:,:,:)
  REAL(dp), PUBLIC, POINTER :: wet3Dso2_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: wet3Dso4_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: dry3Dso2_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: dry3Dso4_inst(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cth_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: sconcso2(:,:)
  REAL(dp), PUBLIC, POINTER :: sconcso4(:,:)
  REAL(dp), PUBLIC, POINTER :: rsdt_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: rsdt_tmp(:,:)
  REAL(dp), PUBLIC, POINTER :: lcc_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: eis(:,:)
  REAL(dp), PUBLIC, POINTER :: eis_inst(:,:)
  REAL(dp), PUBLIC, POINTER :: t3d(:,:,:)

  CONTAINS

  !------------------------------------------------
  SUBROUTINE init_HEaci

    USE mo_decomposition, ONLY: dc => local_decomposition
  
    ALLOCATE(nstr_cld_top(dc%nproma, dc%ngpblks))
    ALLOCATE(nstr_cld_liq_top(dc%nproma, dc%ngpblks))
  
  END SUBROUTINE init_HEaci

  !------------------------------------------------

  !------------------------------------------------

  SUBROUTINE construct_HEaci_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE 
    USE mo_memory_base,       ONLY: new_stream, add_stream_element, AUTO,  &
                                    default_stream_setting, add_stream_reference
    USE mo_time_event,        ONLY: io_time_event

    TYPE(io_time_event) :: put_interval

    !-- set output interval
    put_interval%counter      = 3
    put_interval%unit         = 'hours'
    put_interval%adjustment   = 'last'
    put_interval%offset       = 0
    
    !-- Create new stream:
    CALL new_stream (acheaci ,'acheaci', &
                     filetype = AEROCOM_FILETYPE, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)
  
    !-- Add standard fields for post-processing:
    CALL default_stream_setting (acheaci, &
                                 table = 199, &
                                 code = AUTO )

    !-- Basic vars (may sometimes be necessary)
    CALL add_stream_reference (acheaci, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (acheaci, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (acheaci, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (acheaci, 'gboxarea','geoloc',lpost=.TRUE.)
  
    !-- Diagnostics table'
    !-- Diagnostics table'
    CALL add_stream_element (acheaci, 'lts', lts, &
        longname = 'lower tropospheric stability', &
        units = 'K', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'eis', eis, &
        longname = 'estimated inversion strength', &
        units = 'K', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acheaci, 'emiso2', emiso2, &
        longname = 'total emission of SO2', &
        units = 'kg m-2 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'sconcso2', sconcso2, &
        longname = 'Surface concentration of SO2', &
        units = 'kg m-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
   
    CALL add_stream_element (acheaci, 'sconcso4', sconcso4, &
        longname = 'Surface concentration of SO4', &
        units = 'kg m-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
   
    CALL add_stream_element (acheaci, 'emiso2_inst', emiso2_inst, &
        longname = 'total emission of SO2 (inst.)', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acheaci, 'ccn01bl', ccn01bl, &
        longname = 'cloud_condensation_nuclei_0.1_pbl', &
        units = 'm-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'ccn03bl', ccn03bl, &
        longname = 'cloud_condensation_nuclei_0.3_pbl', &
        units = 'm-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'od550aer', od550aer, &
        longname = 'atmosphere_optical_thickness_due_to_aerosol', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acheaci, 'od440aer', od440aer, &
        longname = 'atmosphere_optical_thickness_due_to_aerosol', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acheaci, 'od865aer', od865aer, &
        longname = 'atmosphere_optical_thickness_due_to_aerosol', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (acheaci, 'aerindex', aerindex, &
        longname = 'aerosol_index', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
     
     CALL add_stream_element (acheaci, 'aerindex_inst', aerindex_inst, &
        longname = 'aerosol_index (inst.)', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

     CALL add_stream_element (acheaci, 'precip', precip_inst, &
        longname = 'total_precipitation_rate', &
        units = 'km m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'sprecip', sprecip_inst, &
        longname = 'stratiform_precipitation_rate', &
        units = 'km m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'ttop', ttop_inst, &
        longname = 'air_temperature_at_cloud_top', &
        units = 'K', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    ! auxiliary var to correct for intermittent cloud presence
    ! in the 'ttop' calc.
    CALL add_stream_element (acheaci, 'tmf_cld_inst', tmf_cld_inst, &
        longname = 'time fraction of cloud presence', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acheaci, 'cdnc', cdnc_inst, &
        longname = 'liquid_cloud_droplet_number_concentration', &
        units = 'm-3', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'cth', cth_inst, &
        longname = 'liquid_cloud_top_height', &
        units = 'm', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acheaci, 'cdr', cdr_inst, &
        longname = 'liquid_cloud-top_droplet_effective_radius', &
        units = 'm', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'lwp', lwp, &
        longname = 'atmosphere_cloud_liquid_path', &
        units = 'kg m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
   
    CALL add_stream_element (acheaci, 't3d', t3d, &
        longname = 'temperature', &
        units = 'K', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
   
    CALL add_stream_element (acheaci, 'cod', cod_inst, &
        longname = 'cloud_optical_depth', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
        
    CALL add_stream_element (acheaci, 'clt', clt_inst, &
        longname = 'cloud_area_fraction', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'lcc', lcc_inst, &
        longname = 'liquid_cloud_area_fraction', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (acheaci, 'accretn', accretn, &
        longname = 'column_accretion_rate', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'autoconv', autoconv, &
        longname = 'column_autoconversion_rate', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'rsut', rsut_inst, &
        longname = 'toa_upward_shortwave_flux', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (acheaci, 'rsdt', rsdt_inst, &
        longname = 'TOA Incident Shortwave Radiation', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acheaci, 'rsutcs', rsutcs_inst, &
        longname = 'toa_upward_shortwave_flux_assuming_clear_sky', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     !Set lradforcing to true
     CALL add_stream_element (acheaci, 'rsutnoa', rsutnoa_inst, &
        longname = 'toa_upward_shortwave_flux_no_aerosol', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (acheaci, 'rsutcsnoa', rsutcsnoa_inst, &
        longname = 'toa_upward_shortwave_flux_clear_sky_no_aerosol', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (acheaci, 'rlutnoa', rlutnoa_inst, &
        longname = 'toa_upward_longwave_flux_no_aerosol', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (acheaci, 'rlutcsnoa', rlutcsnoa_inst, &
        longname = 'toa_upward_longwave_flux_clear_sky_no_aerosol', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (acheaci, 'rlut', rlut_inst, &
        longname = 'toa_upward_longwave_flux', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (acheaci, 'rlutcs', rlutcs_inst, &
        longname = 'toa_upward_longwave_flux_assuming_clear_sky', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acheaci, 'sut', rsut_tmp, &
        longname = 'toa_upward_shortwave_flux', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )
    
    CALL add_stream_element (acheaci, 'sdt', rsdt_tmp, &
        longname = 'TOA Incident Shortwave Radiation', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acheaci, 'sutcs', rsutcs_tmp, &
        longname = 'toa_upward_shortwave_flux_assuming_clear_sky', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )

     !Set lradforcing to true
     CALL add_stream_element (acheaci, 'sutnoa', rsutnoa_tmp, &
        longname = 'toa_upward_shortwave_flux_no_aerosol', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )
    
     CALL add_stream_element (acheaci, 'sutcsnoa', rsutcsnoa_tmp, &
        longname = 'toa_upward_shortwave_flux_clear_sky_no_aerosol', &
        units = 'W m-2', &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )
    
     CALL add_stream_element (acheaci, 'f3d', f3d, &
        longname = 'cloud fraction', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (acheaci, 'fice3d', fice3d, &
        longname = 'ice+mixed phase cloud fraction', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (acheaci, 'cod3d', cod3d, &
        longname = 'cloud optical thickness', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (acheaci, 'cod3dswl', cod3dswl, &
        longname = 'cloud optical thickness-water', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
     
     CALL add_stream_element (acheaci, 'cod3dswi', cod3dswi, &
        longname = 'cloud optical thickness-ice', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

     CALL add_stream_element (acheaci, 'cdnc3d', cdnc3d, &
        longname = 'in-cloud droplet number concentration', &
        units = 'm-3', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (acheaci, 'cdr3d', cdr3d, &
        longname = 'in-cloud  droplet effective radius', &
        units = 'm', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (acheaci, 'phase3d', phase3d, &
        longname = 'cloud thermodynamic phase', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (acheaci, 'ccn_inst', ccn_inst, &
        longname = 'cloud_condensation_nuclei (inst.)', &
        units = 'm-3', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
     
     CALL add_stream_element (acheaci, 'ccn1_inst', ccn1_inst, &
        longname = 'cloud_condensation_nuclei_0.1 (inst.)', &
        units = 'm-3', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
    
     CALL add_stream_element (acheaci, 'ccn3_inst', ccn3_inst, &
        longname = 'cloud_condensation_nuclei_0.3 (inst.)', &
        units = 'm-3', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
     
     CALL add_stream_element (acheaci, 'ccn5_inst', ccn5_inst, &
        longname = 'cloud_condensation_nuclei_0.5 (inst.)', &
        units = 'm-3', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
     
     CALL add_stream_element (acheaci, 'ccn10_inst', ccn10_inst, &
        longname = 'cloud_condensation_nuclei_1.0 (inst.)', &
        units = 'm-3', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
     
     CALL add_stream_element (acheaci, 'z_inst', z_inst, &
        longname = 'altitude (inst.)', &
        units = 'm', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

     CALL add_stream_element (acheaci, 'accret3d', accret3d, &
        longname = 'Accretion of cloud liquid by rain(3D)', &
        units = 'kg kg-1 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
     
     CALL add_stream_element (acheaci, 'auto3d', auto3d, &
        longname = 'Autoconversion of cloud liquid(3D)', &
        units = 'kg kg-1 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
     
     CALL add_stream_element (acheaci, 'dpg', dpg, &
        longname = 'delta p over g', &
        units = 'kg m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
    
     CALL add_stream_element (acheaci, 'riming3d', riming3d, &
        longname = 'collection of liquid cloud by ice(3D)', &
        units = 'kg kg-1 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

     CALL add_stream_element (acheaci, 'qcsedten3d', qcsedten3d, &
        longname = 'rate of droplet settling(3D)', &
        units = 'kg kg-1 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

     CALL add_stream_element (acheaci, 'qcevap3d', qcevap3d, &
        longname = 'Rate of evaporation of falling cloud liquid(3D)', &
        units = 'kg kg-1 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

     CALL add_stream_element (acheaci, 'od550aer3d', od550aer3d, &
        longname = 'atmosphere_optical_thickness_due_to_aerosol(3D)', &
        units = '1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

     CALL add_stream_element (acheaci, 'wet3Dso2_inst', wet3Dso2_inst, &
        longname = 'wet deposition of SO2(inst.)', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
    
     CALL add_stream_element (acheaci, 'wet3Dso4_inst', wet3Dso4_inst, &
        longname = 'wet deposition of SO4(inst.)', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
    
     CALL add_stream_element (acheaci, 'dry3Dso2_inst', dry3Dso2_inst, &
        longname = 'dry deposition of SO2(inst.)', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
    
     CALL add_stream_element (acheaci, 'dry3Dso4_inst', dry3Dso4_inst, &
        longname = 'dry deposition of SO4(inst.)', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )
    
    CALL add_stream_element (acheaci, 'eis_inst', eis_inst, &
        longname = 'estimated inversion strength', &
        units = 'K', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

  END SUBROUTINE construct_HEaci_stream

  SUBROUTINE update_HEaci_diags(kproma, kbdim, klev, krow, pi0)

    USE mo_vphysc,        ONLY: vphysc
    USE mo_memory_g3b,    ONLY: tpot, geosp, tsurf_na, relhum
    USE mo_memory_g1a,    ONLY: tm1, xlm1, xim1, xtm1
    USE mo_time_control,  ONLY: delta_time
    USE mo_ham_streams,   ONLY: tau_2d
    USE mo_ham_rad_data,  ONLY: Nwv_sw
    USE mo_cosp_simulator, ONLY: cisccp_cldtau3d
    USE mo_physical_constants, ONLY: grav
    USE mo_tracdef,       ONLY: ntrac, trlist, AEROSOLMASS, GAS
    USE mo_species,       ONLY: speclist
    USE mo_physical_constants, ONLY : cpd, grav, rd, alv, rv, vtmpc1

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(IN )  :: pi0(kbdim)
    INTEGER :: jl,jk,jt
    INTEGER :: icld_top(kbdim) !< highest cloud top index (strat or conv)
    INTEGER :: itmp(kbdim), zindarr(kbdim), zind850arr(kbdim)
    REAL(dp) :: zeps
    REAL(dp) :: ztemp(kbdim,klev), gamma_m_850(kbdim), t850(kbdim), qsat(kbdim), esat(kbdim), lcl(kbdim), z700(kbdim)
    LOGICAL :: ll1(kbdim), lll1(kbdim,klev)
    
    zeps=EPSILON(1.0_dp)
    
    !-- t3d
    t3d(1:kproma,:,krow) = tm1(1:kproma,:,krow)

    !-- lts
    lts(1:kproma,:,krow) = lts(1:kproma,:,krow) + tpot(1:kproma,:,krow)*delta_time

    !-- eis
    zind850arr(1:kproma) = MINLOC(ABS(vphysc%apm1(1:kproma,:,krow)-85000._dp),2,(.true.))
    DO jl=1,kproma
       t850(jl) = tm1(jl,zind850arr(jl),krow)
    END DO
    esat(1:kproma) = 610.78_dp*EXP(17.27_dp*(t850(1:kproma)-273.15_dp)/(t850(1:kproma)-273.15_dp+237.3_dp)) !Tetens equation
    qsat(1:kproma) = esat(1:kproma)*rd/rv / (85000._dp - vtmpc1 * esat(1:kproma)*rd/rv)
    qsat(1:kproma) = MERGE(qsat(1:kproma),0._dp,qsat(1:kproma)>=0._dp)
    gamma_m_850(1:kproma) = grav/cpd*(1._dp-(1._dp+alv*qsat(1:kproma)/rd/t850(1:kproma))/&
         (1._dp+alv*alv*qsat(1:kproma)/cpd/rv/t850(1:kproma)/t850(1:kproma)))
    z700(1:kproma) = rd*t850(1:kproma)/grav*log(vphysc%aphm1(1:kproma,klev+1,krow)/70000._dp)
    z700(1:kproma) = MERGE(z700(1:kproma),0._dp,z700(1:kproma)>=0._dp)
    lcl(1:kproma) = cpd/grav*(t850(1:kproma)-55._dp-1._dp/(1._dp/&
         (t850(1:kproma)-55._dp)-log(MAX(relhum(1:kproma,klev,krow),0.01_dp))/2840._dp)) ! Eq. (22) of Bolton (1980) for lifting condensation level
    lcl(1:kproma) = MERGE(lcl(1:kproma),0._dp,lcl(1:kproma)>=0._dp)
    eis_inst(1:kproma,krow) = gamma_m_850(1:kproma)*(z700(1:kproma)-lcl(1:kproma)) ! Naud et al. (2016)
    eis_inst(1:kproma,krow) = MERGE(eis_inst(1:kproma,krow),0._dp,eis_inst(1:kproma,krow)>=0._dp)
    eis(1:kproma,krow) = eis(1:kproma,krow) + eis_inst(1:kproma,krow)*delta_time
 
    !-- emiso2
    emiso2(1:kproma,krow) = emiso2(1:kproma,krow) + emiso2_inst(1:kproma,krow)*delta_time

    !-- sconcso2, sconcso4
    DO jt=1, ntrac
       IF (trlist%ti(jt)%spid/=0)THEN
          IF (trlist%ti(jt)%nphase == GAS .AND. speclist(trlist%ti(jt)%spid)%shortname == 'SO2') &
               sconcso2(1:kproma,krow) = sconcso2(1:kproma,krow) + xtm1(1:kproma,klev,jt,krow)*delta_time
          IF (trlist%ti(jt)%nphase == AEROSOLMASS .AND. speclist(trlist%ti(jt)%spid)%shortname == 'SO4') &
               sconcso4(1:kproma,krow) = sconcso4(1:kproma,krow) + xtm1(1:kproma,klev,jt,krow)*delta_time
       END IF
    END DO

    ! ccn PBL at 0.1% and 0.3% supersaturation
    DO jk=1,klev
       z_inst(1:kproma,jk,krow) = (vphysc%geom1(1:kproma,jk,krow)+geosp(1:kproma,krow)) / grav
    END DO
    ztemp(1:kproma,:) = z_inst(1:kproma,:,krow)
    zindarr(1:kproma) = MINLOC(ABS(ztemp(1:kproma,:)-1000._dp),2,(.true.))
    DO jl=1,kproma
       ccn01bl(jl,krow) = ccn1_inst(jl,zindarr(jl),krow)
       ccn03bl(jl,krow) = ccn3_inst(jl,zindarr(jl),krow)
    END DO
    
    !-- od550aer
    od550aer(1:kproma,krow) = od550aer(1:kproma,krow) &
         + tau_2d(Nwv_sw+1)%ptr(1:kproma,krow) * delta_time

    !-- od440aer
    od440aer(1:kproma,krow) = od440aer(1:kproma,krow) &
         + tau_2d(Nwv_sw+3)%ptr(1:kproma,krow) * delta_time

    !-- od865aer
    od865aer(1:kproma,krow) = od865aer(1:kproma,krow) &
         + tau_2d(Nwv_sw+2)%ptr(1:kproma,krow) * delta_time
    
    !-- aerindex
    aerindex(1:kproma,krow) = aerindex(1:kproma,krow) + &
         aerindex_inst(1:kproma,krow)*delta_time

    !-- ttop
    icld_top(1:kproma) = nstr_cld_top(1:kproma,krow)

    ll1(1:kproma) = (icld_top(1:kproma) >= 1) &
        .AND. (icld_top(1:kproma) <= klev)

    !--- Set dummy index values where icld_top is out of range
    !    in order to circumvent a techn. pb with some compilers
    itmp(1:kproma) = MERGE(icld_top(1:kproma),1,ll1(1:kproma))

    DO jk=1,kproma
        ttop_inst(jk,krow) =MERGE(tm1(jk,itmp(jk),krow)*f3d(jk,itmp(jk),krow), 0._dp, ll1(jk))
    ENDDO
    ! Records cloud presence over the accumulation interval
    ! in order to correct ttop at p-proc (ttop = ttop / tmf_cld)
    tmf_cld_inst(1:kproma,krow) = MERGE(1.0_dp, 0._dp, ll1(1:kproma))

    !-- cdnc, cdr, cth
    ll1(1:kproma) = (nstr_cld_liq_top(1:kproma,krow) >= 1) &
         .AND. (nstr_cld_liq_top(1:kproma,krow) <= klev) 
 
    !--- Set dummy index values where nstr_cld_liq_top is out of range
    !    in order to circumvent a techn. pb with some compilers
    itmp(1:kproma) = MERGE(nstr_cld_liq_top(1:kproma,krow),1,ll1(1:kproma))
 
    DO jk=1,kproma
       cdnc_inst(jk,krow) = MERGE(cdnc3d(jk,itmp(jk),krow)*f3d(jk,itmp(jk),krow),0._dp,ll1(jk))
       cdr_inst(jk,krow) = MERGE(cdr3d(jk,itmp(jk),krow)*f3d(jk,itmp(jk),krow),0._dp,ll1(jk))
       cth_inst(jk,krow) = MERGE(z_inst(jk,itmp(jk),krow)*f3d(jk,itmp(jk),krow),0._dp,ll1(jk))
    END DO

    !-- cod
    cod3d(1:kproma,:,krow)=cisccp_cldtau3d(1:kproma,:,krow)*f3d(1:kproma,:,krow)
    cod_inst(1:kproma,krow)     = 0._dp
    DO jk=1,klev
       cod_inst(1:kproma,krow) = cod_inst(1:kproma,krow) + cod3d(1:kproma,jk,krow)
    END DO

    !-- phase3d
    phase3d(1:kproma,:,krow) = xlm1(1:kproma,:,krow) / &
         MAX(xlm1(1:kproma,:,krow)+xim1(1:kproma,:,krow), zeps)
    phase3d(1:kproma,:,krow) = MERGE(MIN(MAX(phase3d(1:kproma,:,krow),0._dp),1._dp),0._dp,&
         f3d(1:kproma,:,krow)>0._dp)
    
    !-- fice3d
    fice3d(1:kproma,:,krow) = f3d(1:kproma,:,krow) * (1._dp - phase3d(1:kproma,:,krow))
    
    !-- precip
    precip_inst(1:kproma,krow) = vphysc%precip_na(1:kproma,krow)
    
    !-- accr, auto, lwp
    accretn(1:kproma,krow)     = 0._dp
    autoconv(1:kproma,krow)     = 0._dp
    lwp(1:kproma,krow)     = 0._dp
    DO jk=1,klev
       accretn(1:kproma,krow)  = accretn(1:kproma,krow)  + accret3d(1:kproma,jk,krow)*&
            dpg(1:kproma,jk,krow)
       autoconv(1:kproma,krow) = autoconv(1:kproma,krow) + auto3d(1:kproma,jk,krow)*&
            dpg(1:kproma,jk,krow)
       lwp(1:kproma,krow)      = lwp(1:kproma,krow)      + xlm1(1:kproma,jk,krow)*&
            dpg(1:kproma,jk,krow)
    END DO
 
    !-- rsut, rsutcs, rsutnoa, rsutcsnoa, rsdt
    rsut_inst(1:kproma,krow)      = pi0(1:kproma) * rsut_tmp(1:kproma,krow) 
    rsutcs_inst(1:kproma,krow)    = pi0(1:kproma) * rsutcs_tmp(1:kproma,krow) 
    rsutnoa_inst(1:kproma,krow)   = pi0(1:kproma) * rsutnoa_tmp(1:kproma,krow) 
    rsutcsnoa_inst(1:kproma,krow) = pi0(1:kproma) * rsutcsnoa_tmp(1:kproma,krow)
    rsdt_inst(1:kproma,krow)      = pi0(1:kproma) * rsdt_tmp(1:kproma,krow)
    
  END SUBROUTINE update_HEaci_diags

END MODULE mo_hammoz_aerocom_HEaci
