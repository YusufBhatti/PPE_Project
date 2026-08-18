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
  REAL(dp), PUBLIC, POINTER :: ccn01(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn03(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn05(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ccn01vi(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn03vi(:,:)
  REAL(dp), PUBLIC, POINTER :: ccn05vi(:,:)

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
         units = 'kg m-2 s-1', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )

    CALL add_stream_element (acmmppe, 'ccn01', ccn01, &
         longname = 'ccn number concentration at SS=0.1%', &
         units = 'm-3', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )

    CALL add_stream_element (acmmppe, 'ccn03', ccn03, &
         longname = 'ccn number concentration at SS=0.3%', &
         units = 'm-3', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )

    CALL add_stream_element (acmmppe, 'ccn05', ccn05, &
         longname = 'ccn number concentration at SS=0.5%', &
         units = 'm-3', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )

    CALL add_stream_element (acmmppe, 'ccn01vi', ccn01vi, &
         longname = 'vertically integrated CCN number concentration at S=0.1%', &
         units = 'm-2', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )

    CALL add_stream_element (acmmppe, 'ccn03vi', ccn03vi, &
         longname = 'vertically integrated CCN number concentration at S=0.3%', &
         units = 'm-2', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )

    CALL add_stream_element (acmmppe, 'ccn05vi', ccn05vi, &
         longname = 'vertically integrated CCN number concentration at S=0.5%', &
         units = 'm-2', &
         laccu = .FALSE., &
         lpost = .TRUE., &
         lrerun = .TRUE. )
    
   !-- Diagnostics table'

  END SUBROUTINE construct_MMPPE_stream

  SUBROUTINE update_MMPPE_diags(kproma, kbdim, klev, krow)

    USE mo_ham,            ONLY: nclass, sizeclass
    USE mo_ham_streams,    ONLY: rdry, a, b
    USE mo_ham_tools,      ONLY: ham_m7_logtail
    USE mo_vphysc,         ONLY: vphysc
    USE mo_memory_g1a,     ONLY: xtm1
    USE mo_hammoz_aerocom_HEaci, ONLY: f3d, phase3d
    USE mo_physical_constants, ONLY: grav

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow

    INTEGER :: jclass,it

    REAL(dp) :: zfracn3(kbdim,klev,nclass),zfracn50(kbdim,klev,nclass)
    REAL(dp) :: zrdry(kbdim,klev,nclass),rcrit(kbdim,klev)
    REAL(dp) :: cfracn(nclass)

    ! -- CCN at fixed supersaturations (SS=0.1/0.3/0.5%), same method as mo_ham_ccn
    INTEGER,  PARAMETER :: nsat_mmppe = 3
    REAL(dp), PARAMETER :: zsat_mmppe(nsat_mmppe) = (/0.001_dp, 0.003_dp, 0.005_dp/)
    REAL(dp), PARAMETER :: zeps_mmppe = EPSILON(1._dp)
    INTEGER  :: jsat_mmppe
    REAL(dp) :: ztmp_mmppe
    REAL(dp) :: zra_mmppe(kbdim,klev), zfracn_mmppe(kbdim,klev)
    REAL(dp) :: zccn_mmppe(kbdim,klev,nsat_mmppe)
    REAL(dp) :: zdpg_mmppe(kbdim,klev)
    
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

!-- ccn01, ccn03, ccn05 (SS=0.1/0.3/0.5%), following mo_ham_ccn's ham_ccn method
    zccn_mmppe(1:kproma,:,:) = 0._dp
    DO jsat_mmppe = 1, nsat_mmppe
       ztmp_mmppe = (2._dp/zsat_mmppe(jsat_mmppe))**(2._dp/3._dp)
       DO jclass=1, nclass
          IF (.NOT. sizeclass(jclass)%lactivation) CYCLE
          zra_mmppe(1:kproma,:) = 1.0_dp ! large value => zero activated fraction by default
          WHERE (b(jclass)%ptr(1:kproma,:,krow) > zeps_mmppe)
             zra_mmppe(1:kproma,:) = a(jclass)%ptr(1:kproma,:,krow)                    &
                  / (3._dp * b(jclass)%ptr(1:kproma,:,krow)**(1._dp/3._dp)) * ztmp_mmppe
          END WHERE
          CALL ham_m7_logtail(kproma,    kbdim,  klev,   krow, jclass, &
                              .TRUE.,   rdry(jclass)%ptr(:,:,krow),   &
                              zra_mmppe, zfracn_mmppe)
          zccn_mmppe(1:kproma,:,jsat_mmppe) = zccn_mmppe(1:kproma,:,jsat_mmppe)         &
               + zfracn_mmppe(1:kproma,:)                                              &
                 * xtm1(1:kproma,:,sizeclass(jclass)%idt_no,krow)                      &
                 * vphysc%rhoam1(1:kproma,:,krow)
       END DO
    END DO

    ccn01(1:kproma,:,krow) = zccn_mmppe(1:kproma,:,1)
    ccn03(1:kproma,:,krow) = zccn_mmppe(1:kproma,:,2)
    ccn05(1:kproma,:,krow) = zccn_mmppe(1:kproma,:,3)

    !-- vertical integration (dp/g weighting, same as burden diagnostics)
    zdpg_mmppe(1:kproma,1) = 2._dp*(vphysc%aphm1(1:kproma,2,krow)-vphysc%apm1(1:kproma,1,krow))/grav
    zdpg_mmppe(1:kproma,2:klev) = (vphysc%aphm1(1:kproma,3:klev+1,krow)-vphysc%aphm1(1:kproma,2:klev,krow))/grav

    ccn01vi(1:kproma,krow) = SUM(zccn_mmppe(1:kproma,:,1)/vphysc%rhoam1(1:kproma,:,krow)*zdpg_mmppe(1:kproma,:), DIM=2)
    ccn03vi(1:kproma,krow) = SUM(zccn_mmppe(1:kproma,:,2)/vphysc%rhoam1(1:kproma,:,krow)*zdpg_mmppe(1:kproma,:), DIM=2)
    ccn05vi(1:kproma,krow) = SUM(zccn_mmppe(1:kproma,:,3)/vphysc%rhoam1(1:kproma,:,krow)*zdpg_mmppe(1:kproma,:), DIM=2)
    
  END SUBROUTINE update_MMPPE_diags

END MODULE mo_hammoz_aerocom_MMPPE
