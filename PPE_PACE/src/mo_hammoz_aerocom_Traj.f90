!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_Traj.f90
!!
!! \brief
!! Module for AeroCom TRAJ_DIAG_PACKAGE diagnostics
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
MODULE mo_hammoz_aerocom_Traj

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  USE mo_submodel_diag, ONLY: vmem3d
   
  IMPLICIT NONE  

  PRIVATE

  PUBLIC :: construct_Traj_stream
  PUBLIC :: construct_Trajg_stream
  PUBLIC :: update_Traj_diags

  TYPE (t_stream), PUBLIC, POINTER :: actraj
  TYPE (t_stream), PUBLIC, POINTER :: actrajg
  
  REAL(dp), PUBLIC, POINTER :: ps(:,:)
  REAL(dp), PUBLIC, POINTER :: zmla(:,:)
  REAL(dp), PUBLIC, POINTER :: ua10m(:,:)
  REAL(dp), PUBLIC, POINTER :: va10m(:,:)
  REAL(dp), PUBLIC, POINTER :: t2m(:,:)
  REAL(dp), PUBLIC, POINTER :: precip(:,:)
  REAL(dp), PUBLIC, POINTER :: hfss(:,:)
  REAL(dp), PUBLIC, POINTER :: emissa(:,:)
  REAL(dp), PUBLIC, POINTER :: emissc(:,:)
  REAL(dp), PUBLIC, POINTER :: emidust(:,:)
  REAL(dp), PUBLIC, POINTER :: emidms(:,:)
  REAL(dp), PUBLIC, POINTER :: emibc(:,:)
  
  REAL(dp), PUBLIC, POINTER :: zgeo(:,:,:)
  REAL(dp), PUBLIC, POINTER :: t(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ua(:,:,:)
  REAL(dp), PUBLIC, POINTER :: va(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hus(:,:,:)
  REAL(dp), PUBLIC, POINTER :: omega(:,:,:)
  REAL(dp), PUBLIC, POINTER :: hur(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rho(:,:,:)
  REAL(dp), PUBLIC, POINTER :: plev(:,:,:)
  REAL(dp), PUBLIC, POINTER :: airmass(:,:,:)
  REAL(dp), PUBLIC, POINTER :: zh(:,:,:)
  REAL(dp), PUBLIC, POINTER :: laythick(:,:,:)
  REAL(dp), PUBLIC, POINTER :: so4nucl(:,:,:)

  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: conccnmode(:)        

  CONTAINS

  !------------------------------------------------

  SUBROUTINE construct_Trajg_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE_GRIB 
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
    CALL new_stream (actrajg ,'actrajg', &
                     filetype = AEROCOM_FILETYPE_GRIB, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)

    !-- Add standard fields for post-processing:   
    CALL default_stream_setting (actrajg, &
                                 table = 199, &
                                 code = AUTO )
 
    !-- Diagnostics table'
    CALL add_stream_element (actrajg, 'ps', ps, &
        longname = 'surface_air_pressure', &
        units = 'Pa', &
        code = 134 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (actrajg, 'zmla', zmla, &
        longname = 'atmosphere_boundary_layer_thickness', &
        units = 'm', &
        code = 159 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (actrajg, 'ua10m', ua10m, &
        longname = 'eastward_wind', &
        units = 'm s-1', &
        code = 165 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actrajg, 'va10m', va10m, &
        longname = 'northward_wind', &
        units = 'm s-1', &
        code = 166 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actrajg, 't2m', t2m, &
        longname = 'air_temperature', &
        units = 'K', &
        code = 167 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actrajg, 'precip', precip, &
        longname = 'lwe_thickness_of_precipitation_amount', &
        units = 'm', &
        code = 228 , &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actrajg, 'hfss', hfss, &
        longname = 'surface_upward_sensible_heat_flux', &
        units = 'J m-2', &
        code = 146 , &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (actrajg, 'zgeo', zgeo, &
        longname = 'geopotential', &
        units = 'm2 s-2', &
        code = 129 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actrajg, 't', t, &
        longname = 'air_temperature', &
        units = 'K', &
        code = 130 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (actrajg, 'ua', ua, &
        longname = 'eastward_wind', &
        units = 'm s-1', &
        code = 131 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (actrajg, 'va', va, &
        longname = 'northward_wind', &
        units = 'm s-1', &
        code = 132 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (actrajg, 'hus', hus, &
        longname = 'specific_humidity', &
        units = 'kg kg-1', &
        code = 133 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
   CALL add_stream_element (actrajg, 'omega', omega, &
        longname = 'lagrangian_tendency_of_air_pressure', &
        units = 'Pa s-1', &
        code = 135 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (actrajg, 'hur', hur, &
        longname = 'relative_humidity', &
        units = '%', &
        code = 157 , &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
  END SUBROUTINE construct_Trajg_stream

  SUBROUTINE construct_Traj_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE
    USE mo_memory_base,       ONLY: new_stream, add_stream_element, AUTO,  &
                                    default_stream_setting, add_stream_reference
    USE mo_time_event,        ONLY: io_time_event
    USE mo_ham,           ONLY: nclass, sizeclass

    TYPE(io_time_event) :: put_interval

    INTEGER :: jclass
    
    IF (.NOT. ALLOCATED(conccnmode)) ALLOCATE(conccnmode(nclass))
    
    !-- set output interval
    put_interval%counter      = 3
    put_interval%unit         = 'hours'
    put_interval%adjustment   = 'last'
    put_interval%offset       = 0
    
    !-- Create new stream:
    CALL new_stream (actraj ,'actraj', &
                     filetype = AEROCOM_FILETYPE, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)
  
    !-- Add standard fields for post-processing:
    CALL default_stream_setting (actraj, &
                                 table = 199, &
                                 code = AUTO )
    
    !-- Diagnostics table'
!>> Only needed for development experiment
    CALL add_stream_reference (actraj, 'ps',     'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'zmla',   'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'ua10m',  'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'va10m',  'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 't2m',    'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'precip', 'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'hfss',   'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'zgeo',   'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 't',      'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'ua',     'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'va',     'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'hus',    'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'omega',  'actrajg' ,lpost=.TRUE.)
    CALL add_stream_reference (actraj, 'hur',    'actrajg' ,lpost=.TRUE.)
!<<
    
    CALL add_stream_element (actraj, 'rho', rho, &
        longname = 'air_density', &
        units = 'kg m-3', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actraj, 'plev', plev, &
        longname = 'air_pressure', &
        units = 'Pa', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actraj, 'airmass', airmass, &
        longname = 'atmosphere_mass_of_air_per_unit_aera', &
        units = 'kg m-2', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actraj, 'zh', zh, &
        longname = 'height', &
        units = 'm', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actraj, 'laythick', laythick, &
        longname = 'cell_thickness', &
        units = 'm', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    DO jclass = 1, nclass
       CALL add_stream_reference (actraj, 'rdry_'//TRIM(sizeclass(jclass)%shortname), 'ham',&
            ref_name = 'ddrymode'//TRIM(sizeclass(jclass)%shortname), &
            ref_longname = 'dry diameter of mode '//TRIM(sizeclass(jclass)%shortname), &
            lpost=.TRUE.)

       CALL add_stream_element (actraj, 'conccnmode'//TRIM(sizeclass(jclass)%shortname), &
            conccnmode(jclass)%ptr, &
            longname = 'number concentration of mode '//TRIM(sizeclass(jclass)%shortname), &
            units = 'm-3', &
            laccu = .FALSE., &
            lpost = .TRUE., &
            lrerun = .TRUE. )
    END DO

    CALL add_stream_element (actraj, 'emissa', emissa, &
        longname = 'emission rate of seasalt dry aerosol particles due to emission into accumulation mode', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actraj, 'emissc', emissc, &
        longname = 'emission rate of seasalt dry aerosol particles due to emission into coarse mode', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actraj, 'emidust', emidust, &
        longname = 'tendency_of_atmosphere_mass_content_of_dust_dry_aerosol_particles_due_to_emission', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actraj, 'emidms', emidms, &
        longname = 'tendency_of_atmosphere_mass_content_of_dimethyl_sulfide_due_to_emission', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    CALL add_stream_element (actraj, 'emibc', emibc, &
        longname = 'tendency_of_atmosphere_mass_content_of_elemental_carbon_dry_aerosol_particles_due_to_emission', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (actraj, 'so4nucl', so4nucl, &
        longname = 'nucleation of H2SO4 to nucleation mode', &
        units = 'kg m-2 s-1', &
        laccu = .FALSE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )
    
    !Tracer stream has to be set to 3 hourly output
    
  END SUBROUTINE construct_Traj_stream
  
  SUBROUTINE update_Traj_diags(kproma, kbdim, klev, krow)
    
    USE mo_vphysc,       ONLY: vphysc
    USE mo_physical_constants, ONLY: grav, rhoh2o
    USE mo_memory_g3b,   ONLY: u10, v10, temp2, aps, ahfs_na
    USE mo_ham_streams,  ONLY: relhum
    USE mo_time_control, ONLY: delta_time
    USE mo_memory_g2a,   ONLY: um1, vm1
    USE mo_memory_g1a,   ONLY: tm1, qm1, xtm1
    USE mo_activ,        ONLY: w_large
    USE mo_scan_buffer,  ONLY: vervel
    USE mo_ham,          ONLY: nclass, sizeclass

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow
    INTEGER :: jk, jl
    INTEGER :: ipbl
    INTEGER :: jclass
    REAL(dp), POINTER :: conccnmode_p(:,:,:)

    !-- ps
    ps(1:kproma,krow) = aps(1:kproma,krow)
    
    !-- zmla
    DO jk=1,kproma
       ipbl = INT(vphysc %pbl(jk,krow))
       zmla(jk,krow) = vphysc%geom1(jk,ipbl,krow) / grav
    ENDDO
     
    !-- ua10m
    ua10m(1:kproma,krow) = u10(1:kproma,krow)

    !-- va10m
    va10m(1:kproma,krow) = v10(1:kproma,krow)

    !-- t2m
    t2m(1:kproma,krow) =temp2(1:kproma,krow)
    
    !-- precip
    precip(1:kproma,krow) = precip(1:kproma,krow) + vphysc%precip_na(1:kproma,krow)/rhoh2o*delta_time &
         * delta_time

    !-- hfss
    hfss(1:kproma,krow) = hfss(1:kproma,krow) + ahfs_na(1:kproma,krow)*delta_time &
         * delta_time
    
    !-- zgeo
    zgeo(1:kproma,:,krow) = vphysc%geom1(1:kproma,:,krow)

    !-- t
    t(1:kproma,:,krow) = tm1(1:kproma,:,krow)
    
    !-- ua
    ua(1:kproma,:,krow) = um1(1:kproma,:,krow)

    !-- va
    va(1:kproma,:,krow) = vm1(1:kproma,:,krow)

    !-- hus
    hus(1:kproma,:,krow) = qm1(1:kproma,:,krow)

     !-- omega
    omega(1:kproma,:,krow) = vervel(1:kproma,:,krow)
    
    !-- hur
    hur(1:kproma,:,krow) = 100.0_dp*relhum(1:kproma,:,krow)
    
    !-- rho
    rho(1:kproma,:,krow) = vphysc%rhoam1_moist(1:kproma,:,krow)
    
    !-- plev
    plev(1:kproma,:,krow) = vphysc%apm1(1:kproma,:,krow)

    !-- airmass
    airmass(1:kproma,:,krow) =  vphysc%rhoam1_moist(1:kproma,:,krow) &
         * vphysc%grheightm1(1:kproma,:,krow)
    
    !-- zh
    zh(1:kproma,:,krow) = vphysc%geom1(1:kproma,:,krow) / grav

    !-- laythick
    laythick(1:kproma,:,krow) = vphysc%grheightm1(1:kproma,:,krow)

    !-- conccnmodeXX
    DO jclass = 1, nclass
       conccnmode_p     => conccnmode(jclass)%ptr
       conccnmode_p(1:kproma,:,krow) = xtm1(1:kproma,:,sizeclass(jclass)%idt_no,krow) * &
            vphysc%rhoam1(1:kproma,:,krow)
    END DO
       
  END SUBROUTINE update_Traj_diags

END MODULE mo_hammoz_aerocom_Traj
