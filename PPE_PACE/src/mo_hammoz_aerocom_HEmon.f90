!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_HEmon.f90
!!
!! \brief
!! Module for AeroCom Holuhraun Monthly diagnostics
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
MODULE mo_hammoz_aerocom_HEmon

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  USE mo_submodel_diag, ONLY: vmem3d

  IMPLICIT NONE  

  PRIVATE

  PUBLIC :: construct_HEmon_stream
  PUBLIC :: update_HEmon_diags

  TYPE (t_stream), PUBLIC, POINTER :: achemon

  REAL(dp), PUBLIC, POINTER :: emiso2(:,:)
  REAL(dp), PUBLIC, POINTER :: od550aer(:,:)
  REAL(dp), PUBLIC, POINTER :: precip(:,:)
  REAL(dp), PUBLIC, POINTER :: sprecip(:,:)
  REAL(dp), PUBLIC, POINTER :: ttop(:,:)
  REAL(dp), PUBLIC, POINTER :: tmf_cld(:,:)
  REAL(dp), PUBLIC, POINTER :: cdnc(:,:)
  REAL(dp), PUBLIC, POINTER :: cdr(:,:)
  REAL(dp), PUBLIC, POINTER :: cod(:,:)
  REAL(dp), PUBLIC, POINTER :: codliq(:,:)
  REAL(dp), PUBLIC, POINTER :: codice(:,:)
  REAL(dp), PUBLIC, POINTER :: clt(:,:)
  REAL(dp), PUBLIC, POINTER :: lcc(:,:)
  REAL(dp), PUBLIC, POINTER :: rsut(:,:)
  REAL(dp), PUBLIC, POINTER :: rsutcs(:,:)
  REAL(dp), PUBLIC, POINTER :: rsutnoa(:,:)
  REAL(dp), PUBLIC, POINTER :: rsutcsnoa(:,:)
  REAL(dp), PUBLIC, POINTER :: rlut(:,:)
  REAL(dp), PUBLIC, POINTER :: rlutcs(:,:)
  REAL(dp), PUBLIC, POINTER :: rlutnoa(:,:)
  REAL(dp), PUBLIC, POINTER :: rlutcsnoa(:,:)
  
  REAL(dp), PUBLIC, POINTER :: rho(:,:,:)
  REAL(dp), PUBLIC, POINTER :: t(:,:,:)
  REAL(dp), PUBLIC, POINTER :: pfull(:,:,:)
  REAL(dp), PUBLIC, POINTER :: phalf(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hus(:,:,:)
  REAL(dp), PUBLIC, POINTER :: airmass(:,:,:)
  REAL(dp), PUBLIC, POINTER :: z(:,:,:)
  REAL(dp), PUBLIC, POINTER :: laythick(:,:,:)
  REAL(dp), PUBLIC, POINTER :: mmrso2(:,:,:)
  REAL(dp), PUBLIC, POINTER :: mmrso4(:,:,:)
  REAL(dp), PUBLIC, POINTER :: dry3Dso2(:,:,:)
  REAL(dp), PUBLIC, POINTER :: dry3Dso4(:,:,:)
  REAL(dp), PUBLIC, POINTER :: wet3Dso2(:,:,:)
  REAL(dp), PUBLIC, POINTER :: wet3Dso4(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ec5503Daer(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn1(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn3(:,:,:)
  REAL(dp), PUBLIC, POINTER :: lwc(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rel(:,:,:)
  REAL(dp), PUBLIC, POINTER :: lccl(:,:,:)
  REAL(dp), PUBLIC, POINTER :: mccl(:,:,:)
  REAL(dp), PUBLIC, POINTER :: iwc(:,:,:)
  REAL(dp), PUBLIC, POINTER :: iccl(:,:,:)
  REAL(dp), PUBLIC, POINTER :: accretl(:,:,:)
  REAL(dp), PUBLIC, POINTER :: autocl(:,:,:)
  REAL(dp), PUBLIC, POINTER :: qcsedten(:,:,:)
  REAL(dp), PUBLIC, POINTER :: qcevap(:,:,:)
  REAL(dp), PUBLIC, POINTER :: riming(:,:,:)
  REAL(dp), PUBLIC, POINTER :: lts(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cth(:,:)
  REAL(dp), PUBLIC, POINTER :: eis(:,:)
  REAL(dp), PUBLIC, POINTER :: rh(:,:,:)

  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: conccnmode(:)        
  
  CONTAINS

  !------------------------------------------------

  SUBROUTINE construct_HEmon_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE 
    USE mo_memory_base,       ONLY: new_stream, add_stream_element, AUTO,  &
                                    default_stream_setting, add_stream_reference
    USE mo_time_event,        ONLY: io_time_event
    USE mo_ham,           ONLY: nclass, sizeclass
    USE mo_control,           ONLY: nlev

    TYPE(io_time_event) :: put_interval

    INTEGER :: jclass
    
    IF (.NOT. ALLOCATED(conccnmode)) ALLOCATE(conccnmode(nclass))

    !-- set output interval
    put_interval%counter      = 1
    put_interval%unit         = 'months'
    put_interval%adjustment   = 'last'
    put_interval%offset       = 0
    
    !-- Create new stream:
    CALL new_stream (achemon ,'achemon', &
                     filetype = AEROCOM_FILETYPE, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)
  
    !-- Add standard fields for post-processing:
    CALL default_stream_setting (achemon, &
                                 table = 199, &
                                 code = AUTO )

    !-- Basic vars (may sometimes be necessary)
    CALL add_stream_reference (achemon, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (achemon, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (achemon, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (achemon, 'gboxarea','geoloc',lpost=.TRUE.)
  
    !-- Diagnostics table'
    CALL add_stream_element (achemon, 'lts', lts, &
        longname = 'lower tropospheric stability', &
        units = 'K', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    
    CALL add_stream_element (achemon, 'eis', eis, &
        longname = 'estimated inversion strength', &
        units = 'K', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (achemon, 'cth', cth, &
        longname = 'liquid_cloud_top_height', &
        units = 'm', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    !emi-stream must be set to monthly output
    CALL add_stream_reference (achemon, 'emi_SO2',     'emi' , &
         ref_name = 'emiso2', &
         ref_longname = 'emission of SO2', &
         lpost=.TRUE.)
    
    CALL add_stream_element (achemon, 'od550aer', od550aer, &
        longname = 'atmosphere_optical_thickness_due_to_aerosol', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (achemon, 'precip', precip, &
        longname = 'total_precipitation_rate', &
        units = 'km m-2 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (achemon, 'sprecip', sprecip, &
        longname = 'stratiform_precipitation_rate', &
        units = 'km m-2 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (achemon, 'ttop', ttop, &
        longname = 'air_temperature_at_cloud_top', &
        units = 'K', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    ! auxiliary var to correct for intermittent cloud presence
    ! in the 'ttop' calc.
    CALL add_stream_element (achemon, 'tmf_cld', tmf_cld, &
        longname = 'time fraction of cloud presence', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (achemon, 'cdnc', cdnc, &
        longname = 'liquid_cloud_droplet_number_concentration', &
        units = 'm-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (achemon, 'cdr', cdr, &
        longname = 'liquid_cloud-top_droplet_effective_radius', &
        units = 'm', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    !make sure the echam stream is output at monthly frequency
    CALL add_stream_reference (achemon, 'xlvi', 'g3b', &
        ref_name = 'lwp', &
        ref_longname = 'atmosphere_cloud_liquid_path', &
        ref_units = 'kg m-2', &
        lpost = .TRUE. )

    !make sure the echam stream is output at monthly frequency
    CALL add_stream_reference (achemon, 'xivi', 'g3b', &
        ref_name = 'iwp', &
        ref_longname = 'atmosphere_cloud_ice_path', &
        ref_units = 'kg m-2', &
        lpost = .TRUE. )

    CALL add_stream_element (achemon, 'cod', cod, &
        longname = 'cloud_optical_depth', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (achemon, 'codliq', codliq, &
        longname = 'cloud_optical_depth_due_to_liq', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (achemon, 'codice', codice, &
        longname = 'cloud_optical_depth_due_to_ice', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (achemon, 'clt', clt, &
        longname = 'cloud_area_fraction', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'lcc', lcc, &
        longname = 'liquid_cloud_area_fraction', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rsut', rsut, &
        longname = 'toa_upward_shortwave_flux', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rsutcs', rsutcs, &
        longname = 'toa_upward_shortwave_flux_assuming_clear_sky', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rsutnoa', rsutnoa, &
        longname = 'toa_upward_shortwave_flux_no_aerosol', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rsutcsnoa', rsutcsnoa, &
        longname = 'toa_upward_shortwave_flux_clear_sky_no_aerosol', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rlut', rlut, &
        longname = 'toa_upward_longwave_flux', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rlutcs', rlutcs, &
        longname = 'toa_upward_longwave_flux_assuming_clear_sky', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rlutnoa', rlutnoa, &
        longname = 'toa_upward_longwave_flux_no_aerosol', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rlutcsnoa', rlutcsnoa, &
        longname = 'toa_upward_longwave_flux_clear_sky_no_aerosol', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rho', rho, &
        longname = 'air density', &
        units = 'kg m-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 't', t, &
        longname = 'absolute temperature', &
        units = 'K', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'pfull', pfull, &
        longname = 'Air_pressure', &
        units = 'Pa', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'phalf', phalf, &
        longname = 'air_pressure_at_interfaces', &
        klev = nlev+1, &
        units = 'Pa', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'hus', hus, &
        longname = 'specific_humidity', &
        units = 'kg kg-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'rh', rh, &
        longname = 'Relative Humidity', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (achemon, 'airmass', airmass, &
        longname = 'atmosphere_mass_content_of_air', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'z', z, &
        longname = 'altitude', &
        units = 'm', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
      
     CALL add_stream_element (achemon, 'laythick', laythick, &
        longname = 'layer geometrical thickness', &
        units = 'm', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     !ham stream must be set to monthly output
     DO jclass = 1, nclass
        CALL add_stream_reference (achemon, 'rdry_'//TRIM(sizeclass(jclass)%shortname), 'ham',&
             ref_name = 'ddrymode'//TRIM(sizeclass(jclass)%shortname), &
             ref_longname = 'dry diameter of mode '//TRIM(sizeclass(jclass)%shortname), &
             lpost=.TRUE.)
        
        CALL add_stream_element (achemon, 'conccnmode'//TRIM(sizeclass(jclass)%shortname), &
             conccnmode(jclass)%ptr, &
             longname = 'number concentration of mode '//TRIM(sizeclass(jclass)%shortname), &
             units = 'm-3', &
             laccu = .FALSE., &
             lpost = .TRUE., &
             lrerun = .TRUE. )
     END DO
     
     CALL add_stream_element (achemon, 'mmrso2', mmrso2, &
        longname = 'Mass Mixing Ratio of SO2', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'mmrso4', mmrso4, &
        longname = 'Mass Mixing Ratio of SO4', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     !ham stream must be set to monthly output
     CALL add_stream_reference (achemon, 'D_PROD_H2SO4', 'ham', &
         ref_name = 'chegpso4', &
         ref_longname = 'gas phase production so4', &
         ref_units = 'kg m-2 s-1', &
         lpost=.TRUE.)

     !ham stream must be set to monthly output
     CALL add_stream_reference (achemon, 'D_PROD_SO4_LIQ', 'ham', &
         ref_name = 'cheaqpso4', &
         ref_longname = 'aqu phase production so4', &
         ref_units = 'kg m-2 s-1', &
         lpost=.TRUE.)
       
     CALL add_stream_element (achemon, 'dry3Dso2', dry3Dso2, &
        longname = 'dry deposition of SO2', &
        units = 'kg m-2 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'dry3Dso4', dry3Dso4, &
        longname = 'dry deposition of SO4', &
        units = 'kg m-2 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'wet3Dso2', wet3Dso2, &
        longname = 'wet deposition of SO2', &
        units = 'kg m-2 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'wet3Dso4', wet3Dso4, &
        longname = 'wet deposition of SO4', &
        units = 'kg m-2 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'ec5503Daer', ec5503Daer, &
        longname = 'Aerosol Extinction@550nm', &
        units = 'm-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'ccn', ccn, &
        longname = 'cloud_condensation_nuclei', &
        units = 'm-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
     
     CALL add_stream_element (achemon, 'ccn1', ccn1, &
        longname = 'cloud_condensation_nuclei_0.1', &
        units = 'm-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'ccn3', ccn3, &
        longname = 'cloud_condensation_nuclei_0.3', &
        units = 'm-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
        
     CALL add_stream_element (achemon, 'lwc', lwc, &
        longname = 'cloud_liquid_water_content', &
        units = 'kg m-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
     
     CALL add_stream_element (achemon, 'rel', rel, &
        longname = 'droplet_effective_radius', &
        units = 'm', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
     
     CALL add_stream_element (achemon, 'lccl', lccl, &
        longname = 'liquid_cloud_fraction', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
     
     CALL add_stream_element (achemon, 'mccl', mccl, &
        longname = 'mixed_phase_cloud_fraction', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
     
     CALL add_stream_element (achemon, 'iwc', iwc, &
        longname = 'cloud_ice_water_content', &
        units = 'kg m-3', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
     CALL add_stream_element (achemon, 'iccl', iccl, &
        longname = 'ice_cloud_fraction', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
     
     CALL add_stream_element (achemon, 'accretl', accretl, &
        longname = 'Accretion of cloud liquid by rain', &
        units = 'kg kg-1 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
     
     CALL add_stream_element (achemon, 'autocl', autocl, &
        longname = 'Autoconversion of cloud liquid', &
        units = 'kg kg-1 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
     
     CALL add_stream_element (achemon, 'qcsedten', qcsedten, &
        longname = 'rate of droplet settling', &
        units = 'kg kg-1 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (achemon, 'qcevap', qcevap, &
        longname = 'Rate of evaporation of falling cloud liquid', &
        units = 'kg kg-1 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     CALL add_stream_element (achemon, 'riming', riming, &
        longname = 'collection of liquid cloud by ice', &
        units = 'kg kg-1 s-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

     !Tracer stream has to be accumulated monthly
     
    END SUBROUTINE construct_HEmon_stream

  SUBROUTINE update_HEmon_diags(kproma, kbdim, klev, krow)

    USE mo_vphysc,        ONLY: vphysc
    USE mo_memory_g3b,    ONLY: tpot, geosp, relhum
    USE mo_memory_g1a,    ONLY: tm1, qm1, xtm1
    USE mo_memory_gl ,    ONLY: xl, xi
    USE mo_time_control,  ONLY: delta_time
    USE mo_ham_streams,   ONLY: tau_2d
    USE mo_ham_rad_data,  ONLY: Nwv_sw
    USE mo_hammoz_aerocom_HEaci, ONLY: ttop_inst, cdnc_inst, cdr_inst, cod_inst, ccn_inst,&
         ccn1_inst, ccn3_inst, z_inst, cod3dswl, cod3dswi, precip_inst, sprecip_inst,&
         auto3d, accret3d, rsut_inst, rsutcs_inst, rsutnoa_inst, rsutcsnoa_inst, rlut_inst,&
         rlutcs_inst, clt_inst, f3d, phase3d, cdr3d, riming3d, qcsedten3d, qcevap3d,&
         od550aer3d, wet3Dso2_inst, wet3Dso4_inst, tmf_cld_inst, dry3Dso2_inst, dry3Dso4_inst,&
         cth_inst, lcc_inst, eis_inst, rlutnoa_inst, rlutcsnoa_inst
    USE mo_physical_constants,ONLY: tmelt, grav
    USE mo_echam_cloud_params, ONLY: cthomi
    USE mo_ham,          ONLY: nclass, sizeclass
    USE mo_hammoz_aerocom_MMPPE, ONLY: fliq2d
    USE mo_param_switches, ONLY: lMMPPE
    
    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow
    INTEGER :: jk
    REAL(dp) :: codliq_inst(kbdim),codice_inst(kbdim), tcc(kbdim)
    REAL(dp) :: zacltot(kbdim)     ! total cloud cover seen from above down to the current
    REAL(dp) :: zdeltacc(kbdim)    ! utility var 
    REAL(dp) :: ztmp1_1d(kbdim), ztmp2_1d(kbdim), ztmp3_1d(kbdim), lcc_tmp(kbdim)
    LOGICAL :: ll_vis(kbdim,klev)
    LOGICAL :: ll_liq(kbdim,klev)
    LOGICAL :: ll_ice(kbdim,klev)
    LOGICAL :: ll_mix(kbdim,klev)
    INTEGER :: jclass
    REAL(dp), POINTER :: conccnmode_p(:,:,:)

    !-- lts
    lts(1:kproma,:,krow) = lts(1:kproma,:,krow) + tpot(1:kproma,:,krow)*delta_time
    
    !-- eis
    eis(1:kproma,krow) = eis(1:kproma,krow) + eis_inst(1:kproma,krow)*delta_time

    !-- cth
    cth(1:kproma,krow) = cth(1:kproma,krow) + cth_inst(1:kproma,krow)*delta_time

    !-- od550aer
    od550aer(1:kproma,krow) = od550aer(1:kproma,krow) &
         + tau_2d(Nwv_sw+1)%ptr(1:kproma,krow) * delta_time

    !-- clt
    clt(1:kproma,krow) = clt(1:kproma,krow) + clt_inst(1:kproma,krow)*delta_time

    !-- lcc
    lcc_tmp(1:kproma) = 0._dp
    tcc(1:kproma)     = 1._dp
    zacltot(1:kproma) = 0._dp
    ll_vis(1:kproma,:) = (f3d(1:kproma,:,krow) > 0.001_dp)
    DO jk = 2,klev
       ztmp1_1d(1:kproma) = MAX(f3d(1:kproma,jk,krow), f3d(1:kproma,jk-1,krow))
       ztmp2_1d(1:kproma) = MIN(f3d(1:kproma,jk-1,krow), 1._dp - 0.001_dp)

       ztmp3_1d(1:kproma) = tcc(1:kproma) * (1._dp - ztmp1_1d(1:kproma)) &
                                               / (1._dp - ztmp2_1d(1:kproma))

       zacltot(1:kproma) = MERGE(ztmp3_1d(1:kproma), zacltot(1:kproma), ll_vis(1:kproma,jk))
 
       zdeltacc(1:kproma) = tcc(1:kproma) - zacltot(1:kproma) ! additional cl frac
                                                                   ! seen from above this level
                                                                   ! because at this stage tcc = zacltot
                                                                   ! one level above current:
                      ! delta_cc = 1-tcc(jk) - (1-tcc(jk-1)) = tcc(jk-1)-tcc(jk) 

 
  
       ztmp3_1d(1:kproma) = lcc_tmp(1:kproma) + phase3d(1:kproma,jk,krow) * zdeltacc(1:kproma)
       lcc_tmp(1:kproma) = MERGE(ztmp3_1d(1:kproma), lcc_tmp(1:kproma), ll_vis(1:kproma,jk))

       !-- final:
       tcc(1:kproma) = MERGE(zacltot(1:kproma), tcc(1:kproma), ll_vis(1:kproma,jk))
    ENDDO
    IF (lMMPPE) fliq2d(1:kproma,krow) = lcc_tmp(1:kproma)
    lcc_inst(1:kproma,krow) = lcc_tmp(1:kproma)
    lcc(1:kproma,krow) = lcc(1:kproma,krow) + lcc_tmp(1:kproma)*delta_time
        
    !--lccl, mccl, iccl
    ll_liq(1:kproma,:) = (tm1(1:kproma,:,krow) >= tmelt)
    ll_mix(1:kproma,:) = (tm1(1:kproma,:,krow) < tmelt) .AND. (tm1(1:kproma,:,krow) >= cthomi)
    ll_ice(1:kproma,:) = .NOT. ll_liq(1:kproma,:) .AND. .NOT. ll_mix(1:kproma,:)
    lccl(1:kproma,:,krow) = lccl(1:kproma,:,krow) + delta_time *&
         MERGE(f3d(1:kproma,:,krow), 0._dp, ll_vis(1:kproma,:) .AND. ll_liq(1:kproma,:))
    mccl(1:kproma,:,krow) = mccl(1:kproma,:,krow) + delta_time *&
         MERGE(f3d(1:kproma,:,krow), 0._dp, ll_vis(1:kproma,:) .AND. ll_mix(1:kproma,:))
    iccl(1:kproma,:,krow) = iccl(1:kproma,:,krow) + delta_time *&
         MERGE(f3d(1:kproma,:,krow), 0._dp, ll_vis(1:kproma,:) .AND. ll_ice(1:kproma,:))

    !-- ccn
    ccn(1:kproma,:,krow) = ccn(1:kproma,:,krow) + ccn_inst(1:kproma,:,krow)*delta_time

    !-- ccn1
    ccn1(1:kproma,:,krow) = ccn1(1:kproma,:,krow) + ccn1_inst(1:kproma,:,krow)*delta_time

    !-- ccn3
    ccn3(1:kproma,:,krow) = ccn3(1:kproma,:,krow) + ccn3_inst(1:kproma,:,krow)*delta_time

    !-- z
    z(1:kproma,:,krow) = z(1:kproma,:,krow) + z_inst(1:kproma,:,krow)*delta_time

    !-- ttop
    ttop(1:kproma,krow) = ttop(1:kproma,krow) + ttop_inst(1:kproma,krow)*delta_time
    
    !-- tmf_cld
    tmf_cld(1:kproma,krow) = tmf_cld(1:kproma,krow) + tmf_cld_inst(1:kproma,krow)*delta_time

    !-- cdnc
    cdnc(1:kproma,krow) = cdnc(1:kproma,krow) + cdnc_inst(1:kproma,krow)*delta_time

    !-- cdr
    cdr(1:kproma,krow) = cdr(1:kproma,krow) + cdr_inst(1:kproma,krow)*delta_time
    
    !-- cod
    cod(1:kproma,krow) = cod(1:kproma,krow) + cod_inst(1:kproma,krow)*delta_time    
    codliq_inst(1:kproma)     = 0._dp
    codice_inst(1:kproma)     = 0._dp
    DO jk=1,klev
       codliq_inst(1:kproma) = codliq_inst(1:kproma) + cod3dswl(1:kproma,jk,krow)*&
            f3d(1:kproma,jk,krow)
       codice_inst(1:kproma) = codice_inst(1:kproma) + cod3dswi(1:kproma,jk,krow)*&
            f3d(1:kproma,jk,krow)
    END DO
    codliq(1:kproma,krow) = codliq(1:kproma,krow) + codliq_inst(1:kproma)*delta_time    
    codice(1:kproma,krow) = codice(1:kproma,krow) + codice_inst(1:kproma)*delta_time    

    !-- precip
    precip(1:kproma,krow)  = precip(1:kproma,krow)  + precip_inst(1:kproma,krow)*delta_time
    sprecip(1:kproma,krow) = sprecip(1:kproma,krow) + sprecip_inst(1:kproma,krow)*delta_time
    
    !-- accr, auto, riming
    accretl(1:kproma,:,krow) = accretl(1:kproma,:,krow) + accret3d(1:kproma,:,krow)*delta_time
    autocl(1:kproma,:,krow)  = autocl(1:kproma,:,krow)  + auto3d(1:kproma,:,krow)*delta_time
    riming(1:kproma,:,krow)  = riming(1:kproma,:,krow)  + riming3d(1:kproma,:,krow)*delta_time

    !-- rsut etc.
    rsut(1:kproma,krow)       = rsut(1:kproma,krow)       + rsut_inst(1:kproma,krow)*delta_time
    rsutcs(1:kproma,krow)     = rsutcs(1:kproma,krow)     + rsutcs_inst(1:kproma,krow)*delta_time
    rsutnoa(1:kproma,krow)    = rsutnoa(1:kproma,krow)    + rsutnoa_inst(1:kproma,krow)*delta_time
    rsutcsnoa(1:kproma,krow)  = rsutcsnoa(1:kproma,krow)  + rsutcsnoa_inst(1:kproma,krow)*delta_time
    rlut(1:kproma,krow)       = rlut(1:kproma,krow)       + rlut_inst(1:kproma,krow)*delta_time
    rlutcs(1:kproma,krow)     = rlutcs(1:kproma,krow)     + rlutcs_inst(1:kproma,krow)*delta_time
    rlutnoa(1:kproma,krow)    = rlutnoa(1:kproma,krow)    + rlutnoa_inst(1:kproma,krow)*delta_time
    rlutcsnoa(1:kproma,krow)  = rlutcsnoa(1:kproma,krow)  + rlutcsnoa_inst(1:kproma,krow)*delta_time

    !-- rho
    rho(1:kproma,:,krow) = rho(1:kproma,:,krow) + vphysc%rhoam1_moist(1:kproma,:,krow)*delta_time

    !-- t
    t(1:kproma,:,krow) = t(1:kproma,:,krow) + tm1(1:kproma,:,krow)*delta_time
 
    !-- pfull
    pfull(1:kproma,:,krow) = pfull(1:kproma,:,krow) + vphysc%apm1(1:kproma,:,krow)*delta_time
    
    !-- phalf
    phalf(1:kproma,:,krow) = phalf(1:kproma,:,krow) + vphysc%aphm1(1:kproma,:,krow)*delta_time
    
    !-- hus
    hus(1:kproma,:,krow) = hus(1:kproma,:,krow) + qm1(1:kproma,:,krow)*delta_time

    !-- rh
    rh(1:kproma,:,krow) = rh(1:kproma,:,krow) + relhum(1:kproma,:,krow)*delta_time

    !-- airmass
    airmass(1:kproma,:,krow) =  airmass(1:kproma,:,krow) + vphysc%rhoam1_moist(1:kproma,:,krow) &
         * vphysc%grheightm1(1:kproma,:,krow)*delta_time
    
    !-- laythick
    laythick(1:kproma,:,krow) = laythick(1:kproma,:,krow) + vphysc%grheightm1(1:kproma,:,krow)&
         *delta_time

     !-- conccnmodeXX
    DO jclass = 1, nclass
       conccnmode_p     => conccnmode(jclass)%ptr
       conccnmode_p(1:kproma,:,krow) = xtm1(1:kproma,:,sizeclass(jclass)%idt_no,krow) * &
            vphysc%rhoam1(1:kproma,:,krow)
    END DO
    
    !-- lwc
    lwc(1:kproma,:,krow) = lwc(1:kproma,:,krow) + xl(1:kproma,:,krow)*delta_time
    
    !-- iwc
    iwc(1:kproma,:,krow) = iwc(1:kproma,:,krow) + xi(1:kproma,:,krow)*delta_time
    
    !-- rel
    rel(1:kproma,:,krow) = rel(1:kproma,:,krow) + &
         cdr3d(1:kproma,:,krow)*f3d(1:kproma,:,krow)*delta_time

    !-- qcsedten, qcevap
    qcsedten(1:kproma,:,krow) = qcsedten(1:kproma,:,krow) + qcsedten3d(1:kproma,:,krow)*delta_time
    qcevap(1:kproma,:,krow)   = qcevap(1:kproma,:,krow)   + qcevap3d(1:kproma,:,krow)*delta_time
    
    !-- ec5503Daer
    ec5503Daer(1:kproma,:,krow) = ec5503Daer(1:kproma,:,krow) + &
         od550aer3d(1:kproma,:,krow)/vphysc%grheightm1(1:kproma,:,krow)*delta_time

    !-- wet deposition
    wet3Dso2(1:kproma,:,krow) = wet3Dso2(1:kproma,:,krow) + wet3Dso2_inst(1:kproma,:,krow)*delta_time
    wet3Dso4(1:kproma,:,krow) = wet3Dso4(1:kproma,:,krow) + wet3Dso4_inst(1:kproma,:,krow)*delta_time

    !-- dry deposition
    dry3Dso2(1:kproma,:,krow) = dry3Dso2(1:kproma,:,krow) + dry3Dso2_inst(1:kproma,:,krow)*delta_time
    dry3Dso4(1:kproma,:,krow) = dry3Dso4(1:kproma,:,krow) + dry3Dso4_inst(1:kproma,:,krow)*delta_time

  END SUBROUTINE update_HEmon_diags

END MODULE mo_hammoz_aerocom_HEmon
