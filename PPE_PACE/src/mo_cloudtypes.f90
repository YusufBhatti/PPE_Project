! #define _kmeans_types
MODULE mo_cloudtypes

  USE mo_kind,                    ONLY: dp
  USE mo_memory_base,             ONLY: new_stream, add_stream_element, &
                                        default_stream_setting, AUTO, t_stream, &
                                        add_stream_reference
  USE mo_time_event,              ONLY: io_time_event
  USE mo_time_control,            ONLY: p_bcast_event
  USE mo_linked_list,             ONLY: HYBRID, SURFACE
  USE mo_netcdf,                  ONLY: HYBRID_H
  USE mo_histogram,               ONLY: add_multiple_stream_elements, t_globptr, t_surfptr, t_hist, &
                                        LINSCALE, LOGSCALE, prcp_spec, mxt_spec, frac_spec, &
                                        ten_spec, dp_spec, double_histogram
  USE mo_echam_cloud_params,      ONLY: cqtmin, cthomi
  USE mo_cloud_utils,             ONLY: cloud_type_helper, get_cloud_bounds, fact_tke, &
                                        get_util_var, clc_min, xsec, eps, effective_liquid_radius, &
                                        effective_ice_radius
  USE mo_memory_base,             ONLY: HYBRID_H
  USE mo_p3_fields,               ONLY: lfull_diag


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: init_clouds, cloudtype_diagnostic, nclusters, nclusters_matus, compute_cloud_labels, &
            aggregate_labels_matus
  PUBLIC :: csemter_ct, cstrsol_ct, isccp_ct


! Cloud essentials
  REAL(dp), PUBLIC, POINTER :: cscld(:,:,:), cslf(:,:,:), cshetf(:,:,:), csliqf(:,:,:), csvvel(:,:,:), &
                               csxi(:,:,:), csxl(:,:,:), cst(:,:,:), csp(:,:,:), csq(:,:,:), &
                               csttop(:,:,:), csptop(:,:,:), csdp(:,:,:), csfr(:,:,:), csqhet(:,:,:), &
                               csqwtot(:,:,:), csmass(:,:,:), csfmask(:,:,:), &
                               csprcl(:,:,:), csprci(:,:,:), csprct(:,:,:), cslabel(:,:,:), &
                               csemi(:,:,:), cstrans(:,:,:), csliqof(:,:,:), cslabel_matus(:,:,:), &
                               cscdnc(:,:,:), csicnc(:,:,:), cseffl(:,:,:), cseffi(:,:,:)

  REAL(dp), PUBLIC, POINTER :: csmfll(:,:), csfrain_ls(:,:), csfrain_tot(:,:), &
                               csfmstaed_liq(:,:), csfmstaed_ice(:,:), csfmstaed_mxp(:,:), &
                               cscrelw(:,:), cscresw(:,:), csdayl(:,:), csemter(:,:), cstrsol(:,:), &
                               csclcov(:,:), isccp_freq(:,:), isccp_sunlit(:,:)

  TYPE(t_globptr), POINTER  :: csflabel(:), csemter_ct(:), cstrsol_ct(:)
  TYPE(t_surfptr), POINTER  :: csflc(:), csflc_prcp(:), csflc_prcpr(:), csflc_prcps(:), csmlabel(:), &
                               csfllabel(:), csfllabel_rain(:), clabelprcp_tot(:), clabelprcp_rain(:), &
                               clabelprcp_snow(:), cscrelw_ct(:), cscresw_ct(:), clabeltemp(:), &
                               clabelliqf(:), clabelliqof(:), cscovlabel(:), clabelttemp(:), &
                               clabellf(:), clabelsosif(:), clabelsocmf(:), cscovlabel_matus(:), &
                               clabeldp(:), clabelaxi(:), clabelaxl(:), clabelaxtot(:), clabelattemp(:), &
                               isccp_ct(:)

  INTEGER, parameter :: npredictors = 5
#ifdef _kmeans_types
  INTEGER, parameter :: nclusters   = 5
#else
  INTEGER, parameter :: nclusters   = 7
#endif
  INTEGER, parameter :: nclusters_matus = 4
  REAL(dp) :: centers(nclusters, npredictors)

  ! histogram spec for labels
  TYPE(t_hist) :: label_spec = t_hist(nclusters, 0.5_dp, REAL(nclusters,dp)+0.5_dp, LINSCALE)

  CONTAINS

  SUBROUTINE init_clouds

    USE mo_mpi,                   ONLY: p_parallel_io, p_barrier, p_bcast, p_io
    USE mo_namelist,              ONLY: open_nml, position_nml, POSITIONED
    USE mo_exception,             ONLY: finish

    CHARACTER(LEN=32) :: ichar
    INTEGER           :: i, j, ierr, inml, iunit, ilength, nlength
    TYPE (t_stream), POINTER      :: clouds, iclouds

    !--- Initialize centers

    ! predictors:     lf     hetf   dp     ttop   qwtot
    centers(1,:) = (/ 0.074, 0.041, 0.664, 0.259, 0.674 /)
    centers(2,:) = (/ 0.943, 0.071, 0.225, 0.662, 0.686 /)
    centers(3,:) = (/ 0.017, 0.014, 0.177, 0.279, 0.453 /)
    centers(4,:) = (/ 0.920, 0.754, 0.164, 0.622, 0.670 /)
    centers(5,:) = (/ 0.112, 0.575, 0.179, 0.517, 0.634 /)

    CALL new_stream (clouds,'clouds',lrerun=.false.)
    ! CALL new_stream (iclouds,'iclouds',lrerun=.false., interval=io_time_event(1,'days','last',0) )
    CALL new_stream (iclouds,'iclouds',lrerun=.false.)

    !--- Add standard fields for post-processing:

    CALL add_stream_reference (clouds, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (clouds, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (clouds, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (clouds, 'gboxarea','geoloc',lpost=.TRUE.)
    CALL add_stream_reference (clouds, 'slm','g3b',lpost=.TRUE.)

    CALL add_stream_reference (iclouds, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (iclouds, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (iclouds, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (iclouds, 'gboxarea','geoloc',lpost=.TRUE.)
    CALL add_stream_reference (iclouds, 'slm','g3b',lpost=.TRUE.)

    CALL default_stream_setting (clouds, lpost     = .TRUE. , &
                                        lrerun    = .FALSE. , &
                                        leveltype = HYBRID , &
                                        table     = 199,     &
                                        code      = AUTO     )

    CALL default_stream_setting (clouds, laccu=.FALSE.)

    CALL default_stream_setting (iclouds, lpost     = .TRUE. , &
                                        lrerun    = .FALSE. , &
                                        leveltype = HYBRID , &
                                        table     = 199,     &
                                        code      = AUTO     )
    CALL default_stream_setting (iclouds, laccu=.FALSE.)

    ! CALL default_stream_setting (iclouds, units     = '',          &
    !                                      lrerun    = .FALSE. ,    &
    !                                      lpost     = .TRUE. ,    &
    !                                      laccu     = .FALSE. ,   &
    !                                      reset     = 1.e-50_dp,  &
    !                                      leveltype = SURFACE ,   &
    !                                      contnorest = .TRUE. )

  
  ! cloud parameters
    CALL add_stream_element (clouds, 'cscld', cscld,                   &
                             longname='instantaneous cloud cover', units='1')
    CALL add_stream_element (clouds, 'cslf', cslf,                   &
                             longname='instantaneous liquid fraction', units='1')
    CALL add_stream_element (clouds, 'cshetf', cshetf,                   &
                             longname='instantaneous het. formed ice fraction', units='1')
    CALL add_stream_element (clouds, 'csttop', csttop,                   &
                             longname='instantaneous cloud top temperature', units='K')
    CALL add_stream_element (clouds, 'csdp', csdp,                   &
                             longname='instantaneous cloud pressure difference', units='Pa')
    CALL add_stream_element (clouds, 'csqwtot', csqwtot,                   &
                             longname='instantaneous total cloud water', units='1')
    CALL add_stream_element (clouds, 'cst', cst,                   &
                             longname='instantaneous temperature', units='K')
    CALL add_stream_element (clouds, 'csliqf', csliqf,                   &
                             longname='instantaneous ice formed through liquid fraction', units='1')
    CALL add_stream_element (clouds, 'csliqof', csliqof,                   &
                             longname='instantaneous liquid origin ice fraction', units='1')

  ! observables ONLY OUTPUT IF lfull_diag == .TRUE.
    CALL default_stream_setting (clouds, lpost=lfull_diag)
    CALL add_stream_element (clouds, 'csvvel', csvvel,                   &
                             longname='instantaneous large scale vert. velocity', units='m s-1')
    CALL add_stream_element (clouds, 'csxi', csxi,                   &
                             longname='instantaneous cloud ice', units='1')
    CALL add_stream_element (clouds, 'csxl', csxl,                   &
                             longname='instantaneous cloud water', units='1')
    CALL add_stream_element (clouds, 'csp', csp,                   &
                             longname='instantaneous pressure', units='Pa')
    CALL add_stream_element (clouds, 'csq', csq,                   &
                             longname='instantaneous humidity', units='1')
    CALL add_stream_element (clouds, 'csptop', csptop,                   &
                             longname='instantaneous cloud top pressure', units='Pa')
    CALL add_stream_element (clouds, 'csfr', csfr,                   &
                             longname='instantaneous rime fraction', units='1')
    CALL add_stream_element (clouds, 'csqhet', csqhet,                   &
                             longname='instantaneous het. formed ice', units='1')
    CALL add_stream_element (clouds, 'csprcl', csprcl,                   &
                             longname='instantaneous surface precip liquid', units='mm h-1')
    CALL add_stream_element (clouds, 'csprci', csprci,                   &
                             longname='instantaneous surface precip ice', units='mm h-1')
    CALL add_stream_element (clouds, 'csprct', csprct,                   &
                             longname='instantaneous surface precip total', units='mm h-1')
    CALL add_stream_element (clouds, 'csmass', csmass,                   &
                             longname='instantaneous air mass (dp/g)', units='kg m-2')
    CALL add_stream_element (clouds, 'csdayl', csdayl,                   &
                             longname='instantaneous day light flag', units='1')
    CALL add_stream_element (clouds, 'cscdnc', cscdnc,                   &
                             longname='instantaneous cloud droplet number', units='kg-1')
    CALL add_stream_element (clouds, 'csicnc', csicnc,                   &
                             longname='instantaneous ice crystal number', units='kg-1')
    CALL add_stream_element (clouds, 'cseffl', cseffl,                   &
                             longname='instantaneous cloud droplet radius', units='m')
    CALL add_stream_element (clouds, 'cseffi', cseffi,                   &
                             longname='instantaneous ice crystal radius', units='m')

    CALL default_stream_setting (clouds, lpost=.TRUE.)
  ! cloud parameters
    CALL add_stream_element (clouds, 'cslabel', cslabel,                   &
                             longname='instantaneous labels', units='1')
    CALL add_stream_element (clouds, 'cslabel_matus', cslabel_matus,       &
                             longname='instantaneous labels acc. matus', units='1')
    CALL add_stream_element (clouds, 'csmfll', csmfll,                   &
                             longname='instantaneous most frequent low label', units='1')

  ! cloud radiative effect per cloud type
    CALL add_multiple_stream_elements(clouds, 'csemter_ct', 2*nclusters+nclusters_matus, csemter_ct, &
                                      lpost=.FALSE., leveltype=HYBRID_H)
    CALL add_multiple_stream_elements(clouds, 'cstrsol_ct', 2*nclusters+nclusters_matus, cstrsol_ct, &
                                      lpost=.FALSE., leveltype=HYBRID_H)

! --------------------------------ACCUMULATED STREAMS---------------------------
    CALL default_stream_setting (clouds, laccu=.TRUE.)

  ! frequency diagnostics
    CALL add_stream_element (clouds, 'csfmask', csfmask,                   &
                             longname='classification frequency', units='1')
    CALL add_stream_element (clouds, 'csfrain_ls', csfrain_ls,                   &
                             longname='large-scale rain frequency', units='1')
    CALL add_stream_element (clouds, 'csfrain_tot', csfrain_tot,                   &
                             longname='rain frequency (ls + conv)', units='1')
    CALL add_stream_element (clouds, 'csclcov', csclcov,                   &
                             longname='maximum random overlap cover', units='1')

  ! label frequency
    CALL add_multiple_stream_elements(clouds, 'csflabel', nclusters, csflabel)
    CALL add_multiple_stream_elements(clouds, 'csmlabel', nclusters, csmlabel)
    CALL add_multiple_stream_elements(clouds, 'cscovlabel', nclusters, cscovlabel)
    CALL add_multiple_stream_elements(clouds, 'cscovlabel_matus', nclusters_matus, cscovlabel_matus)
    CALL add_multiple_stream_elements(clouds, 'csfllabel', nclusters, csfllabel)
    CALL add_multiple_stream_elements(clouds, 'csfllabel_rain', nclusters, csfllabel_rain)

  ! precip diagnostics
    CALL add_multiple_stream_elements(clouds, 'csflc', nclusters, csflc)
    CALL add_multiple_stream_elements(clouds, 'csflc_prcp', nclusters, csflc_prcp)
    CALL add_multiple_stream_elements(clouds, 'csflc_prcpr', nclusters, csflc_prcpr)
    CALL add_multiple_stream_elements(clouds, 'csflc_prcps', nclusters, csflc_prcps)

  ! radiation diagnostics
    CALL add_multiple_stream_elements(clouds, 'cscrelw_ct', 2*nclusters+nclusters_matus, cscrelw_ct)
    CALL add_multiple_stream_elements(clouds, 'cscresw_ct', 2*nclusters+nclusters_matus, cscresw_ct)
    CALL add_stream_element (clouds, 'cscrelw', cscrelw,                   &
                             longname='longwave cloud effect', units='W m-2')
    CALL add_stream_element (clouds, 'cscresw', cscresw,                   &
                             longname='shortwave cloud effect', units='W m-2')

  ! muelmenstaed diagnostics
    CALL add_stream_element (clouds, 'csfmstaed_liq', csfmstaed_liq,                   &
                             longname='liquid origin rain frequency', units='1')
    CALL add_stream_element (clouds, 'csfmstaed_ice', csfmstaed_ice,                   &
                             longname='ice origin rain frequency', units='1')
    CALL add_stream_element (clouds, 'csfmstaed_mxp', csfmstaed_mxp,                   &
                             longname='mixed-phase origin rain frequency', units='1')

! --------------------------------PRECIPITATION HISTOGRAMS----------------------
    ! label - precip histograms
    CALL add_multiple_stream_elements(clouds, 'clabelprcp_tot', label_spec%nbin*prcp_spec%nbin, clabelprcp_tot)
    CALL add_multiple_stream_elements(clouds, 'clabelprcp_rain', label_spec%nbin*prcp_spec%nbin, clabelprcp_rain)
    CALL add_multiple_stream_elements(clouds, 'clabelprcp_snow', label_spec%nbin*prcp_spec%nbin, clabelprcp_snow)

    ! label - temperature histograms
    CALL add_multiple_stream_elements(clouds, 'clabeltemp', label_spec%nbin*mxt_spec%nbin, clabeltemp)
    ! label - top temperature histograms
    CALL add_multiple_stream_elements(clouds, 'clabelttemp', label_spec%nbin*mxt_spec%nbin, clabelttemp)
    ! label - liquid fraction
    CALL add_multiple_stream_elements(clouds, 'clabellf', label_spec%nbin*frac_spec%nbin, clabellf)
    ! label - frozen liquid fraction
    CALL add_multiple_stream_elements(clouds, 'clabelliqf', label_spec%nbin*frac_spec%nbin, clabelliqf)
    ! label - liquid origin ice fraction
    CALL add_multiple_stream_elements(clouds, 'clabelliqof', label_spec%nbin*frac_spec%nbin, clabelliqof)
    ! label - source / sink ratio
    CALL add_multiple_stream_elements(clouds, 'clabelsosif', label_spec%nbin*frac_spec%nbin, clabelsosif)
    ! label - source / cloud mass ratio
    CALL add_multiple_stream_elements(clouds, 'clabelsocmf', label_spec%nbin*ten_spec%nbin, clabelsocmf)
    ! label - cloud pressure difference
    CALL add_multiple_stream_elements(clouds, 'clabeldp', label_spec%nbin*dp_spec%nbin, clabeldp)

    ! --------------------------------SIMPLE AVERAGES PER LABEL----------------------
    ! total, ice and liquid water
    CALL add_multiple_stream_elements(clouds, 'clabelaxtot', nclusters, clabelaxtot)
    CALL add_multiple_stream_elements(clouds, 'clabelaxl', nclusters, clabelaxl)
    CALL add_multiple_stream_elements(clouds, 'clabelaxi', nclusters, clabelaxi)
    CALL add_multiple_stream_elements(clouds, 'clabelattemp', nclusters, clabelattemp)

    ! --------------------------------CRAZY ISCCP HISTOGRAMS FOR EACH CLOUD TYPE ----------------------
    CALL add_multiple_stream_elements(iclouds, 'isccp_ct', (nclusters+1)*7*7, isccp_ct)
    CALL add_stream_element (iclouds, 'isccp_freq', isccp_freq,                   &
                             longname='Frequency of cosp output', units='h')
    CALL add_stream_element (iclouds, 'isccp_sunlit', isccp_sunlit,                   &
                             longname='Frequency of sunlit grid points', units='h')

  END SUBROUTINE init_clouds

SUBROUTINE cloudtype_diagnostic(&
               !--IN
               kproma, kbdim, klev, klevp1, ktrac, ktdia, krow, &
               paphm1, papm1, pgeo, paclc, &
               pt, pq, pxl, pxi, pxt, pssfl, prsfl, pssfc, prsfc, ictop, &
               ptkem1, pvervel, ptvm1, pemter, pemtef0, ptrsol, ptrsof0, &
               pemter_ct, ptrsol_ct, &
               pi0, rdayl_x)

  USE mo_p3_fields, ONLY: idt_qihet, idt_qirim, idt_qiliq, idt_qioliq, &
                          idt_qsrc, idt_qprc, idt_birim
  USE mo_activ,     ONLY: idt_icnc, idt_cdnc
  USE mo_physical_constants, ONLY: rgrav, rd
  USE mo_memory_g3b, ONLY: trad0_na, traf0_na
  USE mo_time_control, ONLY: delta_time, l_trigrad

  INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevp1, ktrac, ktdia, krow

  REAL(dp), INTENT(IN) :: paphm1(kbdim,klevp1), papm1(kbdim,klev)
  REAL(dp), INTENT(IN) :: paclc(kbdim,klev)   !< cloud cover (1)
  REAL(dp), INTENT(IN) :: pgeo(kbdim,klev)    !< geopotential height at full levels
  REAL(dp), INTENT(IN) :: pt(kbdim,klev)      !< temperature (K)
  REAL(dp), INTENT(IN) :: pq(kbdim,klev)      !< water vapor mmr (1)
  REAL(dp), INTENT(IN) :: pxl(kbdim,klev)     !< liquid water mmr (1)
  REAL(dp), INTENT(IN) :: pxi(kbdim,klev)     !< ice water mmr (1)
  REAL(dp), INTENT(IN) :: pxt(kbdim,klev, ktrac) !< tracers

  REAL(dp), INTENT(IN) :: ptkem1(kbdim,klev)  !< turbulent kinetic energy
  REAL(dp), INTENT(IN) :: pvervel(kbdim,klev) !< large-scale vertical velocity (m s-1)
  REAL(dp), INTENT(IN) :: ptvm1(kbdim,klev)   !< potential temperature (K)

  REAL(dp), INTENT(IN) :: pssfl(kbdim)        !< large-scale snow surface flux (kg m-2 s-1)
  REAL(dp), INTENT(IN) :: prsfl(kbdim)        !< large-scale rain surface flux (kg m-2 s-1)
  REAL(dp), INTENT(IN) :: pssfc(kbdim)        !< convective snow surface flux (kg m-2 s-1)
  REAL(dp), INTENT(IN) :: prsfc(kbdim)        !< convective rain surface flux (kg m-2 s-1)

  REAL(dp), INTENT(IN) :: pemter(kbdim,klevp1)                !< cloudy LW flux (W/m2)
  REAL(dp), INTENT(IN) :: pemtef0(kbdim,klevp1)               !< clear sky LW flux (W/m2)
  REAL(dp), INTENT(IN) :: ptrsol(kbdim,klevp1)                !< cloudy SW flux (W/m2)
  REAL(dp), INTENT(IN) :: ptrsof0(kbdim,klevp1)                !< clear sky SW flux (W/m2)
  !< LW flux per cloud type (W/m2)
  REAL(dp), INTENT(IN) :: pemter_ct(kbdim,klevp1,2*nclusters+nclusters_matus)
  !< SW flux per cloud type (W/m2)
  REAL(dp), INTENT(IN) :: ptrsol_ct(kbdim,klevp1,2*nclusters+nclusters_matus)
  REAL(dp), INTENT(IN) :: pi0(kbdim)                         !< solar zenith angle (?)
  REAL(dp), INTENT(IN) :: rdayl_x(kbdim)                     !< day time flag
  INTEGER,  INTENT(IN) :: ictop(kbdim)                       !< convective cloud top level

  REAL(dp) :: zlf(kbdim,klev)                 !< liquid fraction (1)
  REAL(dp) :: zhetf(kbdim,klev)               !< ice formed through liquid fraction (1)
  REAL(dp) :: zliqf(kbdim,klev)               !< heterogeneously formed ice fraction (1)
  REAL(dp) :: zliqof(kbdim,klev)              !< liquid origin ice fraction (1)
  REAL(dp) :: zsosif(kbdim,klev)              !< source to sink ratio (1)
  REAL(dp) :: zsocmf(kbdim,klev)              !< source to cloud mass ratio (1)
  REAL(dp) :: zfr(kbdim,klev)                 !< rimed ice fraction (1)
  REAL(dp) :: zrho_rcp(kbdim,klev)            !< inverse air density (m3 kg-1)
  REAL(dp) :: zvervx(kbdim,klev)              !< updraft velocity (cm s-1)
  REAL(dp) :: zclcov(kbdim)                   !< maximum-random overlap cover
  REAL(dp) :: zclmass(kbdim)                  !< total cloud mass (sum of cloudy*dpg)
  REAL(dp) :: zaclc(kbdim,klev)               !< effective cloud cover (no empty clouds)
  LOGICAL  :: ll_cc(kbdim,klev)               !< cloud flag

  REAL(dp) :: zcdnc(kbdim,klev)     !< cloud droplet number (kg-1)
  REAL(dp) :: zicnc(kbdim,klev)     !< ice crystal number (kg-1)
  REAL(dp) :: zbirim(kbdim,klev)    !< rime volume (m-3)
  REAL(dp) :: zqirim(kbdim,klev)    !< rimed ice (kg kg-1)

  REAL(dp) :: zeffi(kbdim,klev)     !< effective liquid radius (m)
  REAL(dp) :: zeffl(kbdim,klev)     !< effective ice radius (m)

  LOGICAL, DIMENSION(kbdim,klev) :: ll1_2d
  LOGICAL, DIMENSION(kbdim) :: ll1, ll2, ll3
  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1_2d, ztmp2_2d, ztmp3_2d
  REAL(dp), DIMENSION(kbdim) :: ztmp1, ztmp2, ztmp3
  INTEGER :: jk, jl, nc

! ----------------------------GRIDBOX DESCRIPTION---------------------------------
  REAL(dp) :: zgeoh(kbdim,klevp1)  !< Geopotential height at half levels
  REAL(dp) :: zdz(kbdim,klev)      !< layer thickness [m]
  REAL(dp) :: zdp(kbdim,klev)      !< pressure difference of the layer [Pa]
  REAL(dp) :: zdpg(kbdim,klev)     !< delta p over g [kg/m2]
  REAL(dp) :: zaaa(kbdim,klev)     !< Air density correction needed for the ice crystal fall velocity
  REAL(dp) :: zviscos(kbdim,klev)  !< Dynamic viscosity of water in air

! --------------------------------CLOUD TYPES-----------------------------------
  REAL(dp) :: zdp_cld(kbdim,klev) ! cloud pressure difference (Pa)
  REAL(dp) :: zt_top(kbdim,klev)  ! cloud top temperature (K)
  REAL(dp) :: zp_top(kbdim,klev)  ! cloud top pressure (Pa)
  INTEGER  :: zalltop(kbdim,klev) ! cloud top index
  INTEGER  :: zallbas(kbdim,klev) ! cloud base index
  INTEGER  :: itop(kbdim,klev),         & !< flag for cloud tops
              ibas(kbdim,klev),         & !< flag for cloud bases
              icl_minusbas(kbdim,klev), & !< flag for all cloud levels excepted their base
              icl_minustop(kbdim,klev), & !<  flag for all cloud levels excepted their top (useless for now)
              iclnb(kbdim)                !< number of clouds per column

! --------------------------------CLOUD TYPES-----------------------------------
  REAL(dp) :: predictors(kbdim,klev,npredictors) ! predictors for cloud classification
  INTEGER  :: labels(kbdim,klev)                 ! cloud labels
  INTEGER  :: labels_matus(kbdim,klev)           ! cloud labels aggregates as in Matus&Ecuyer (2017)
  REAL(dp) :: low_mask(kbdim,klev)               ! 1 until first cloud top is detected
  REAL(dp) :: mask(kbdim,klev)                   ! grid-boxes to be considered for labelling
  REAL(dp) :: lmask(kbdim,klev)                  ! multi-purpose label mask
  REAL(dp) :: mask_surf(kbdim)                   ! multi-purpose mask
  REAL(dp) :: rain_mask(kbdim)                   ! mask rainy surface gridboxes
  REAL(dp) :: cnt_low_labels(kbdim,nclusters) ! count the number of boxes for label_i
  REAL(dp) :: cnt_labels(kbdim,nclusters)        ! count the number of boxes for label_i
  REAL(dp) :: most_frequent_low_label(kbdim)     ! most frequent low label (column)
  REAL(dp) :: most_frequent_label(kbdim)         ! most frequent label (column)
  REAL(dp) :: label_flag(kbdim,nclusters)        ! flag for each label (0 or 1)
  REAL(dp) :: ice_flag(kbdim)                    ! flag for ice labels (0 or 1)
  REAL(dp) :: liq_flag(kbdim)                    ! flag for liquid labels (0 or 1)
  REAL(dp) :: mxp_flag(kbdim)                    ! flag for mixed-phase labels (0 or 1)
  REAL(dp) :: cliq_flag(kbdim)                   ! flag for convective liquid origin (0 or 1)
  REAL(dp) :: cice_flag(kbdim)                   ! flag for convective ice origin (0 or 1)

! --------------------------------HISTOGRAMS------------------------------------------
  REAL(dp) :: zlabelprcp_tot(label_spec%nbin*prcp_spec%nbin,kbdim)   ! precipitation label histogram (tot)
  REAL(dp) :: zlabelprcp_rain(label_spec%nbin*prcp_spec%nbin,kbdim)  ! precipitation label histogram (rain)
  REAL(dp) :: zlabelprcp_snow(label_spec%nbin*prcp_spec%nbin,kbdim)  ! precipitation label histogram (snow)
  REAL(dp) :: zlabeltemp(label_spec%nbin*mxt_spec%nbin,kbdim,klev)   ! label temperature histogram
  REAL(dp) :: zlabelttemp(label_spec%nbin*mxt_spec%nbin,kbdim,klev)   ! label top-temperature histogram
  REAL(dp) :: zlabellf(label_spec%nbin*frac_spec%nbin,kbdim,klev)    ! label lf histogram
  REAL(dp) :: zlabelliqf(label_spec%nbin*frac_spec%nbin,kbdim,klev)  ! label liqf histogram
  REAL(dp) :: zlabelliqof(label_spec%nbin*frac_spec%nbin,kbdim,klev) ! label liqof histogram
  REAL(dp) :: zlabelsosif(label_spec%nbin*frac_spec%nbin,kbdim,klev) ! label sosif histogram
  REAL(dp) :: zlabelsocmf(label_spec%nbin*ten_spec%nbin,kbdim,klev)  ! label socmf histogram
  REAL(dp) :: zlabeldp(label_spec%nbin*dp_spec%nbin,kbdim,klev)  ! label dp histogram


!-------------------------------------------------------------------------------------
!                                      CLOUD TYPES
!-------------------------------------------------------------------------------------
ll1_2d(1:kproma,:) = (pxl(1:kproma,:)+pxi(1:kproma,:)) > cqtmin
ztmp1_2d(1:kproma,:) = MERGE(paclc(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))

! compute effective cloud cover
zaclc(1:kproma,:) = MERGE(paclc(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))

! cloud flag
ll_cc(1:kproma,:) = ztmp1_2d(1:kproma,:) > clc_min

CALL get_cloud_bounds( &
        !-- IN
        kproma, kbdim, ktdia, klev, ztmp1_2d, &
        !- OUT
        itop, ibas, icl_minustop, icl_minusbas, iclnb)

CALL cloud_type_helper(&
        !--IN
        kproma, kbdim, ktdia, klev, klevp1, &
        ztmp1_2d, itop, ibas, pt, paphm1, &
        !--OUT
        zdp_cld, zt_top, zp_top, zalltop, zallbas)

!-------------------------------------------------------------------------------------
!                                  UTILITY VARIABLES
!-------------------------------------------------------------------------------------

! compute total cover assuming maximum-random overlap
zclcov(1:kproma) = 1.0_dp-zaclc(1:kproma,1)
   
DO jk = 2,klev
   ztmp1(1:kproma) = MAX(zaclc(1:kproma,jk), zaclc(1:kproma,jk-1))
   ztmp2(1:kproma) = MIN(zaclc(1:kproma,jk-1), xsec)
   
   zclcov(1:kproma) = zclcov(1:kproma) * (1._dp - ztmp1(1:kproma))  &
                      / (1._dp - ztmp2(1:kproma))
END DO
zclcov(1:kproma)  = 1.0_dp-zclcov(1:kproma)
csclcov(1:kproma,krow) = csclcov(1:kproma,krow) + delta_time*zclcov(1:kproma)

! vertical velocity
zrho_rcp(1:kproma,:) = rd*ptvm1(1:kproma,:)/papm1(1:kproma,:)
ztmp1_2d(1:kproma,:) = 100._dp*fact_tke*SQRT(ptkem1(1:kproma,:))
zvervx(1:kproma,:)   = -100._dp*rgrav*pvervel(1:kproma,:)*zrho_rcp(1:kproma,:) + ztmp1_2d(1:kproma,:)

!--- Get several utility variables:
CALL get_util_var( &
        !-- IN
        kproma, kbdim, ktdia, klev, klevp1, &
        paphm1(:,:), pgeo(:,:), papm1(:,:), pt(:,:), &
        !-- OUT
        zgeoh(:,:), zdp(:,:), zdpg(:,:), &
        zdz(:,:), zaaa(:,:), zviscos(:,:) )

! total column-integrated cloud mass
ztmp1_2d(1:kproma,:) = MERGE(1._dp, 0._dp, ll_cc(1:kproma,:))
zclmass(1:kproma) = SUM(ztmp1_2d(1:kproma,:)*zdpg(1:kproma,:), 2)

! liquid fraction
ztmp1_2d(1:kproma,:) = pxi(1:kproma,:) + pxl(1:kproma,:)
zlf(1:kproma,:) = MERGE(pxl(1:kproma,:)/ztmp1_2d(1:kproma,:), 0._dp, &
                        ztmp1_2d(1:kproma,:) > eps)
zlf(1:kproma,:) = MAX(MIN(zlf(1:kproma,:), 1._dp), 0._dp)

! heterogeneously formed ice fraction
zhetf(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qihet)/pxi(1:kproma,:), 0._dp, &
                          pxi(1:kproma,:) > eps)
zhetf(1:kproma,:) = MAX(MIN(zhetf(1:kproma,:), 1._dp), 0._dp)

! ice formed through liquid
zliqf(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qiliq)/pxi(1:kproma,:), 0._dp, &
                          pxi(1:kproma,:) > eps)
zliqf(1:kproma,:) = MAX(MIN(zliqf(1:kproma,:), 1._dp), 0._dp)

! liquid origin ice
zliqof(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qioliq)/pxi(1:kproma,:), 0._dp, &
                           pxi(1:kproma,:) > eps)
zliqof(1:kproma,:) = MAX(MIN(zliqof(1:kproma,:), 1._dp), 0._dp)

! rimed ice fraction
zfr(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qirim)/pxi(1:kproma,:), 0._dp, &
                        pxi(1:kproma,:) > eps)
zfr(1:kproma,:) = MAX(MIN(zfr(1:kproma,:), 1._dp), 0._dp)

! fraction of source and sink terms
ll1_2d(1:kproma,:) = pxt(1:kproma,:,idt_qsrc) > eps
zsosif(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qprc)/pxt(1:kproma,:,idt_qsrc), 0._dp, &
                           pxt(1:kproma,:,idt_qsrc) > eps)
zsosif(1:kproma,:) = MAX(0._dp, MIN(1._dp, zsosif(1:kproma,:)))

! fraction of source and cloud mass
ztmp1_2d(1:kproma,:) = pxi(1:kproma,:) + pxl(1:kproma,:)
zsocmf(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qsrc)/ztmp1_2d(1:kproma,:), 0._dp, &
                           ztmp1_2d(1:kproma,:) > eps)
zsocmf(1:kproma,:) = MAX(ten_spec%vmin, MIN(ten_spec%vmax, zsocmf(1:kproma,:)))

zcdnc(1:kproma,:) = pxt(1:kproma,:,idt_cdnc)
zicnc(1:kproma,:) = pxt(1:kproma,:,idt_icnc)
zbirim(1:kproma,:) = pxt(1:kproma,:,idt_birim)
zqirim(1:kproma,:) = pxt(1:kproma,:,idt_qirim)

! effective hydrometeor radii
zeffl(1:kproma,:) = effective_liquid_radius(kbdim, kproma, klev, zcdnc, pxl)
DO jk=1,klev
   zeffi(1:kproma,jk) = 1.e6_dp*effective_ice_radius(kbdim, kproma, pxi(:,jk), zqirim(:,jk), &
                                                     zicnc(:,jk), zbirim(:,jk))
END DO

!-------------------------------------------------------------------------------------
!                                  INSTANTANEOUS CLOUD
!-------------------------------------------------------------------------------------

cscld(1:kproma,:,krow) = zaclc(1:kproma,:)
cshetf(1:kproma,:,krow) = zhetf(1:kproma,:)
csliqf(1:kproma,:,krow) = zliqf(1:kproma,:)
csliqof(1:kproma,:,krow) = zliqof(1:kproma,:)
csfr(1:kproma,:,krow) = zfr(1:kproma,:)
cslf(1:kproma,:,krow) = zlf(1:kproma,:)
csvvel(1:kproma,:,krow) = zvervx(1:kproma,:)
csxi(1:kproma,:,krow) = pxi(1:kproma,:)
csxl(1:kproma,:,krow) = pxl(1:kproma,:)
cst(1:kproma,:,krow) = pt(1:kproma,:)
csp(1:kproma,:,krow) = papm1(1:kproma,:)
csq(1:kproma,:,krow) = pq(1:kproma,:)
csttop(1:kproma,:,krow) = zt_top(1:kproma,:)
csptop(1:kproma,:,krow) = zp_top(1:kproma,:)
csdp(1:kproma,:,krow) = zdp_cld(1:kproma,:)
csqhet(1:kproma,:,krow) = pxt(1:kproma,:,idt_qihet)
csqwtot(1:kproma,:,krow) = pxi(1:kproma,:) + pxl(1:kproma,:)
csmass(1:kproma,:,krow) = zdpg(1:kproma,:)
csicnc(1:kproma,:,krow) = zicnc(1:kproma,:)
cscdnc(1:kproma,:,krow) = zcdnc(1:kproma,:)
cseffi(1:kproma,:,krow) = zeffi(1:kproma,:)
cseffl(1:kproma,:,krow) = zeffl(1:kproma,:)

! precipitation in mm/h
DO jk=ktdia,klev
   csprcl(1:kproma,jk,krow) = prsfl(1:kproma)*3600._dp
   csprci(1:kproma,jk,krow) = pssfl(1:kproma)*3600._dp
   csprct(1:kproma,jk,krow) = (prsfl(1:kproma) + pssfl(1:kproma))*3600._dp
END DO

! accumulated total radiation fluxes
cscrelw(1:kproma,krow) = cscrelw(1:kproma,krow) + &
                         delta_time*(pemter(1:kproma,1) - pemtef0(1:kproma,1))
cscresw(1:kproma,krow) = cscresw(1:kproma,krow) + &
                         delta_time*(ptrsol(1:kproma,1) - ptrsof0(1:kproma,1))*pi0(1:kproma)
csdayl(1:kproma,krow)  = rdayl_x(1:kproma)

! accumulate the radiative effect per cloud type
DO nc=1,2*nclusters+nclusters_matus
   ! cloud radiative effect: all-sky - (all-sky - type) (1:nclusters)
   IF(nc < nclusters + 1) THEN
      cscrelw_ct(nc)%v(1:kproma,krow) = cscrelw_ct(nc)%v(1:kproma,krow) + &
                                        delta_time*(pemter(1:kproma,1) - pemter_ct(1:kproma,1,nc))
      cscresw_ct(nc)%v(1:kproma,krow) = cscresw_ct(nc)%v(1:kproma,krow) + &
                                        delta_time*(ptrsol(1:kproma,1) - ptrsol_ct(1:kproma,1,nc))*pi0(1:kproma)
   ! cloud radiative effect: type - clear-sky (nclusters:2*nclusters)
   ! cloud radiative effect: type - clear-sky for aggregated (Matus & L'Ecuyer 2017) (2*nclusters,-1) 
   ELSE
      cscrelw_ct(nc)%v(1:kproma,krow) = cscrelw_ct(nc)%v(1:kproma,krow) + &
                                        delta_time*(pemter_ct(1:kproma,1,nc) - pemtef0(1:kproma,1))
      cscresw_ct(nc)%v(1:kproma,krow) = cscresw_ct(nc)%v(1:kproma,krow) + &
                                        delta_time*(ptrsol_ct(1:kproma,1,nc) - ptrsof0(1:kproma,1))*pi0(1:kproma)
   END IF
END DO !nc

!-------------------------------------------------------------------------------------
!                                   CLOUD CLASSIFIER   
!-------------------------------------------------------------------------------------

#ifdef _kmeans_types
   predictors(:,:,1) = remap(kproma, kbdim, klev, 0._dp, 1._dp, .FALSE., zlf)             ! liquid fraction
   predictors(:,:,2) = remap(kproma, kbdim, klev, 0._dp, 1._dp, .FALSE., zhetf)           ! het. formed ice fraction
   predictors(:,:,3) = remap(kproma, kbdim, klev, 0._dp, 90000._dp, .FALSE., zdp_cld)     ! cloud thickness
   predictors(:,:,4) = remap(kproma, kbdim, klev, 180._dp, 330._dp, .FALSE., zt_top)      ! cloud top temperature
   predictors(:,:,5) = remap(kproma, kbdim, klev, 1.e-10_dp, 1.e-2_dp, .TRUE., pxi + pxl) ! total condensate

   labels = classify(kproma, kbdim, klev, predictors)
#else
   labels = classify_phys(kproma, kbdim, klev, zlf, zhetf, zliqof, zdp_cld, zt_top)
#endif

! store frequency of classification for label specific time of occurence
mask(1:kproma,:) = MERGE(1._dp, 0._dp, ll_cc(1:kproma,:))
mask(1:kproma,:) = MERGE(mask(1:kproma,:), 0._dp, pt(1:kproma,:) > 180._dp)
mask(1:kproma,:) = MERGE(mask(1:kproma,:), 0._dp, pt(1:kproma,:) < 330._dp)
csfmask(1:kproma,:,krow) = csfmask(1:kproma,:,krow) + delta_time*mask(1:kproma,:)

labels(1:kproma,:) = labels(1:kproma,:)*mask(1:kproma,:)
cslabel(1:kproma,:,krow) = labels(1:kproma,:)

CALL aggregate_labels_matus(&
        !--IN
        kproma, kbdim, klev, labels, paphm1, &
        !--OUT
        labels_matus)

cslabel_matus(1:kproma,:,krow) = labels_matus(1:kproma,:)


!-------------------------------------------------------------------------------------
!                                     CLOUD COVER
!-------------------------------------------------------------------------------------

! store accumulated label properties
DO nc=1,nclusters
   ztmp1_2d(1:kproma,:) = MERGE(1._dp, 0._dp, labels(1:kproma,:) == nc)
   ztmp1(1:kproma) = SUM(zdpg(1:kproma,:)*ztmp1_2d(1:kproma,:), 2)

   ! label frequency (zonal mean plots)
   csflabel(nc)%v(1:kproma,:,krow) = csflabel(nc)%v(1:kproma,:,krow) + delta_time*ztmp1_2d(1:kproma,:)

   ! label mass (pie plots)
   csmlabel(nc)%v(1:kproma,krow) = csmlabel(nc)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)

   ! total label cover   
   ! compute relative weight of given label
   ztmp2(1:kproma) = MERGE(ztmp1(1:kproma)/zclmass(1:kproma), 0._dp, zclmass(1:kproma) > cqtmin)

   ! weight clouds by their contribution to the total 3D cloud volume in the column
   ztmp1(1:kproma) = zclcov(1:kproma)*ztmp2(1:kproma)
   
   cscovlabel(nc)%v(1:kproma,krow) = cscovlabel(nc)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
END DO !nc

! as above, but only for ice, liquid, mixed-phase and multi-layer (see. Matus&Ecuyer 2017)
DO nc=1,nclusters_matus
   ! check if there are any clouds of the given type in the column. Note that the matus labels
   ! are defined column wide
   ll1(1:kproma) = ANY(labels_matus(1:kproma,:) == nc, 2)
   ztmp1(1:kproma) = MERGE(zclcov(1:kproma), 0._dp, ll1(1:kproma))

   cscovlabel_matus(nc)%v(1:kproma,krow) = cscovlabel_matus(nc)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
END DO !nc


!-------------------------------------------------------------------------------------
!                               PRECIPITATION DIAGNOSTIC
!-------------------------------------------------------------------------------------

! lowest-clouds mask
low_mask(1:kproma,:) = 0._dp
ll1(1:kproma) = .TRUE.
DO jk=klev,1,-1
   low_mask(1:kproma,jk) = MERGE(1._dp, 0._dp, ll1(1:kproma))
   ! keep true until first top is detected. From then on, keep false
   ll1(1:kproma) = (.NOT. itop(1:kproma,jk) == jk) .AND. ll1(1:kproma)
END DO !jk
! apply general mask
low_mask(1:kproma,:) = low_mask(1:kproma,:)*mask(1:kproma,:)

DO nc=1,nclusters
   lmask(1:kproma,:) = MERGE(1._dp, 0._dp, labels(1:kproma,:) == nc)
   cnt_low_labels(1:kproma,nc) = SUM(lmask(1:kproma,:)*low_mask(1:kproma,:), 2)
   cnt_labels(1:kproma,nc) = SUM(lmask(1:kproma,:), 2)
END DO !nc

! discard columns without a lowest cloud
ll1(1:kproma) = SUM(low_mask(1:kproma,:), 2) > 0.5_dp
most_frequent_low_label(1:kproma) = MERGE(MAXLOC(cnt_low_labels(1:kproma,:), 2), 0, ll1(1:kproma))
! discard columns without a cloud
ll1(1:kproma) = SUM(mask(1:kproma,:), 2) > 0.5_dp
most_frequent_label(1:kproma) = MERGE(MAXLOC(cnt_labels(1:kproma,:), 2), 0, ll1(1:kproma))

csmfll(1:kproma,krow) = most_frequent_low_label(1:kproma)

rain_mask(1:kproma) = MERGE(1._dp, 0._dp, prsfl(1:kproma)*3600._dp > 1.e-3_dp)
csfrain_ls(1:kproma,krow) = csfrain_ls(1:kproma,krow) + delta_time*rain_mask(1:kproma)

! precipitation weighted label frequency
DO nc=1,nclusters
   mask_surf(1:kproma) = MERGE(1._dp, 0._dp , most_frequent_low_label(1:kproma) == nc)
   csflc(nc)%v(1:kproma,krow) = csflc(nc)%v(1:kproma,krow) + delta_time*mask_surf(1:kproma)
   csflc_prcp(nc)%v(1:kproma,krow) = csflc_prcp(nc)%v(1:kproma,krow) &
                                   + delta_time*mask_surf(1:kproma)*(prsfl(1:kproma) + pssfl(1:kproma))
   csflc_prcpr(nc)%v(1:kproma,krow) = csflc_prcpr(nc)%v(1:kproma,krow) &
                                    + delta_time*mask_surf(1:kproma)*prsfl(1:kproma)
   csflc_prcps(nc)%v(1:kproma,krow) = csflc_prcps(nc)%v(1:kproma,krow) &
                                    + delta_time*mask_surf(1:kproma)*pssfl(1:kproma)
   csfllabel(nc)%v(1:kproma,krow) = csfllabel(nc)%v(1:kproma,krow) &
                                    + delta_time*mask_surf(1:kproma)
   csfllabel_rain(nc)%v(1:kproma,krow) = csfllabel_rain(nc)%v(1:kproma,krow) &
                                       + delta_time*mask_surf(1:kproma)*rain_mask(1:kproma)
END DO !nc

!-------------------------------------------------------------------------------------
!                             MUELMENSTAEDT DIAGNOSTIC
!-------------------------------------------------------------------------------------

! compute flags for each label
DO nc=1,nclusters
   label_flag(1:kproma,nc) = MERGE(1._dp, 0._dp, most_frequent_low_label(1:kproma) == nc)
END DO !nc

#ifdef _kmeans_types
   ! ice labels are 1 (thick), 3 (cirrus)
   ice_flag(1:kproma) = MIN(1._dp, label_flag(1:kproma,1) + label_flag(1:kproma,3))
   ! liquid label is 2
   liq_flag(1:kproma) = label_flag(1:kproma,2)
   ! mixed-phase labels are 4 (liquid dominated) and 5 (ice dominated)
   mxp_flag(1:kproma) = MIN(1._dp, label_flag(1:kproma,4) + label_flag(1:kproma,5))
#else
   ! ice labels are 1 (thick), 4 (in-situ cirrus), 5 (liquid origin cirrus)
   ice_flag(1:kproma) = MIN(1._dp, label_flag(1:kproma,1) + label_flag(1:kproma,4) &
                                 + label_flag(1:kproma,5))
   ! liquid labels are 2 (low) and 3 (high) liquid clouds
   liq_flag(1:kproma) = MIN(1._dp, label_flag(1:kproma,2) + label_flag(1:kproma,3))
   ! mixed-phase labels are 6 (liquid dominated) and 7 (ice dominated)
   mxp_flag(1:kproma) = MIN(1._dp, label_flag(1:kproma,6) + label_flag(1:kproma,7))
#endif

! estimate convective phase
DO jl=1,kbdim
   ! temperature at convective cloud top
   ztmp1(jl) = pt(jl,ictop(jl))
   ! attribute to 
   ! ice         if ttop < cthomi, 
   ! liquid      if ttop > cthomi
   ! we neglect a potential mixed-phase contribution due to the lack of appropriate
   ! microphysics in the convection module
   cice_flag(jl) = MERGE(1._dp, 0._dp, ztmp1(jl) < cthomi)
   cliq_flag(jl) = MERGE(1._dp, 0._dp, ztmp1(jl) > cthomi)
END DO ! jl

! compute the appropriate rain mask
ztmp1(1:kproma) = MERGE(1._dp, 0._dp, prsfl(1:kproma)*3600._dp > 1.e-3_dp) ! large-scale
ztmp2(1:kproma) = MERGE(1._dp, 0._dp, prsfc(1:kproma)*3600._dp > 1.e-3_dp) ! convective
rain_mask(1:kproma) = MIN(1._dp, ztmp1(1:kproma) + ztmp2(1:kproma)) ! logical or

! accumulated rain time for frequency
csfrain_tot(1:kproma,krow) = csfrain_tot(1:kproma,krow) + delta_time*rain_mask(1:kproma)

! there might be convective and large-scale rain at the same time
! in that case, weight the associated labels with 0.5
ztmp3(1:kproma) = ztmp1(1:kproma) + ztmp2(1:kproma)
ztmp3(1:kproma) = MERGE(1._dp/ztmp3(1:kproma), 0._dp, ztmp3(1:kproma) > 0.5_dp) ! 0.5: any number 1 > x > 0

! write to stream. note that the generic case is assumed to be that there is
! *EITHER* large-scale or convective rain where the lines below simplify to
! csfmstaed_X = csfmstaed_X + delta_time*X_flag
csfmstaed_liq(1:kproma,krow) = csfmstaed_liq(1:kproma,krow) + delta_time*          &
                               ztmp3(1:kproma)*(ztmp1(1:kproma)*liq_flag(1:kproma) & ! large-scale contribution
                                              + ztmp2(1:kproma)*cliq_flag(1:kproma)) ! convective contribution

csfmstaed_ice(1:kproma,krow) = csfmstaed_ice(1:kproma,krow) + delta_time*          &
                               ztmp3(1:kproma)*(ztmp1(1:kproma)*ice_flag(1:kproma) & ! large-scale contribution
                                              + ztmp2(1:kproma)*cice_flag(1:kproma)) ! convective contribution

csfmstaed_mxp(1:kproma,krow) = csfmstaed_mxp(1:kproma,krow) + delta_time*          &
                               ztmp3(1:kproma)*(ztmp1(1:kproma)*mxp_flag(1:kproma))  ! large-scale contribution
                                                                                     ! no convective mixed-phase


!-------------------------------------------------------------------------------------
!                                 LABEL HISTOGRAMS
!-------------------------------------------------------------------------------------

CALL double_histogram(&
        !--IN
        kproma, kbdim, label_spec, prcp_spec, &
        most_frequent_label, (prsfl + pssfl)*3600._dp, &
        !--OUT
        zlabelprcp_tot)
CALL double_histogram(&
        !--IN
        kproma, kbdim, label_spec, prcp_spec, &
        most_frequent_label, prsfl*3600._dp, &
        !--OUT
        zlabelprcp_rain)
CALL double_histogram(&
        !--IN
        kproma, kbdim, label_spec, prcp_spec, &
        most_frequent_label, pssfl*3600._dp, &
        !--OUT
        zlabelprcp_snow)
CALL double_histogram(&
        !--IN
        kproma, kbdim, klev, label_spec, mxt_spec, &
        REAL(labels, dp), pt, &
        !--OUT
        zlabeltemp)
CALL double_histogram(&
        !--IN
        kproma, kbdim, klev, label_spec, mxt_spec, &
        REAL(labels, dp), zt_top, &
        !--OUT
        zlabelttemp)
CALL double_histogram(&
        !--IN
        kproma, kbdim, klev, label_spec, frac_spec, &
        REAL(labels, dp), zlf, &
        !--OUT
        zlabellf)
CALL double_histogram(&
        !--IN
        kproma, kbdim, klev, label_spec, frac_spec, &
        REAL(labels, dp), zliqf, &
        !--OUT
        zlabelliqf)
CALL double_histogram(&
        !--IN
        kproma, kbdim, klev, label_spec, frac_spec, &
        REAL(labels, dp), zliqof, &
        !--OUT
        zlabelliqof)
CALL double_histogram(&
        !--IN
        kproma, kbdim, klev, label_spec, frac_spec, &
        REAL(labels, dp), zsosif, &
        !--OUT
        zlabelsosif)
CALL double_histogram(&
        !--IN
        kproma, kbdim, klev, label_spec, ten_spec, &
        REAL(labels, dp), zsocmf, &
        !--OUT
        zlabelsocmf)
CALL double_histogram(&
        !--IN
        kproma, kbdim, klev, label_spec, dp_spec, &
        REAL(labels, dp), zdp_cld, &
        !--OUT
        zlabeldp)

! NOTE: if there is no cloud, the histogram counter is 0, so we dont need the cloud
!       flag as in the histograms in mo_histogram
DO jk=1,label_spec%nbin*prcp_spec%nbin
   ! count precipitation bin per label
   clabelprcp_tot(jk)%v(1:kproma,krow) = clabelprcp_tot(jk)%v(1:kproma,krow) &
                                       + delta_time*zlabelprcp_tot(jk,1:kproma)
   clabelprcp_rain(jk)%v(1:kproma,krow) = clabelprcp_rain(jk)%v(1:kproma,krow) &
                                        + delta_time*zlabelprcp_rain(jk,1:kproma)
   clabelprcp_snow(jk)%v(1:kproma,krow) = clabelprcp_snow(jk)%v(1:kproma,krow) &
                                        + delta_time*zlabelprcp_snow(jk,1:kproma)
END DO !jk

DO jk=1,label_spec%nbin*mxt_spec%nbin
   ! temperature histogram per label
   ztmp1(1:kproma) = SUM(zlabeltemp(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
   clabeltemp(jk)%v(1:kproma,krow) = clabeltemp(jk)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
   ! top temperature histogram per label
   ztmp1(1:kproma) = SUM(zlabelttemp(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
   clabelttemp(jk)%v(1:kproma,krow) = clabelttemp(jk)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
END DO !jk

DO jk=1,label_spec%nbin*frac_spec%nbin
   ! liquid fraction origin per cloud type
   ztmp1(1:kproma) = SUM(zlabellf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
   clabellf(jk)%v(1:kproma,krow) = clabellf(jk)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
   ! frozen liquid fraction origin per cloud type
   ztmp1(1:kproma) = SUM(zlabelliqf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
   clabelliqf(jk)%v(1:kproma,krow) = clabelliqf(jk)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
   ! liquid origin ice fraction per cloud type
   ztmp1(1:kproma) = SUM(zlabelliqof(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
   clabelliqof(jk)%v(1:kproma,krow) = clabelliqof(jk)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
   ! source to sink ratio per cloud type
   ztmp1(1:kproma) = SUM(zlabelsosif(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
   clabelsosif(jk)%v(1:kproma,krow) = clabelsosif(jk)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
END DO !jk

DO jk=1,label_spec%nbin*ten_spec%nbin
   ! source / cloud mass ratio
   ztmp1(1:kproma) = SUM(zlabelsocmf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
   clabelsocmf(jk)%v(1:kproma,krow) = clabelsocmf(jk)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
END DO

DO jk=1,label_spec%nbin*dp_spec%nbin
   ! cloud pressure difference
   ztmp1(1:kproma) = SUM(zlabeldp(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
   clabeldp(jk)%v(1:kproma,krow) = clabeldp(jk)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
END DO

!-------------------------------------------------------------------------------------
!                                 SIMPLE AVERAGES PER TYPE
!-------------------------------------------------------------------------------------

! average values (weighted by air mass, i.e. need to be divided by total label mass: csmlabel)
DO nc=1,nclusters
   ! liquid, ice and total water contents
   ! get label mask
   ztmp1_2d(1:kproma,:) = MERGE(1._dp, 0._dp, labels(1:kproma,:) == nc)

   ! total cloud water
   ztmp1(1:kproma) = SUM((pxi(1:kproma,:)+pxl(1:kproma,:))*zdpg(1:kproma,:)*ztmp1_2d(1:kproma,:), 2)
   clabelaxtot(nc)%v(1:kproma,krow) = clabelaxtot(nc)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)

   ! liquid water
   ztmp1(1:kproma) = SUM(pxl(1:kproma,:)*zdpg(1:kproma,:)*ztmp1_2d(1:kproma,:), 2)
   clabelaxl(nc)%v(1:kproma,krow) = clabelaxl(nc)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)

   ! ice water
   ztmp1(1:kproma) = SUM(pxi(1:kproma,:)*zdpg(1:kproma,:)*ztmp1_2d(1:kproma,:), 2)
   clabelaxi(nc)%v(1:kproma,krow) = clabelaxi(nc)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)

   ! top temperature
   ztmp1(1:kproma) = SUM(zt_top(1:kproma,:)*zdpg(1:kproma,:)*ztmp1_2d(1:kproma,:), 2)
   clabelattemp(nc)%v(1:kproma,krow) = clabelattemp(nc)%v(1:kproma,krow) + delta_time*ztmp1(1:kproma)
END DO  


END SUBROUTINE cloudtype_diagnostic

! helper routine to calculate cloud labels with a minimal interface
SUBROUTINE compute_cloud_labels(&
              !--IN
              kproma, kbdim, klev, krow, &
              paclc, paphm1, &
              !--OUT
              labels)
              
  USE mo_memory_g1a,        ONLY: xlm1, xim1, tm1, xtm1
  USE mo_p3_fields,         ONLY: idt_qihet, idt_qioliq

  INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
  REAL(dp), INTENT(IN) :: paclc(kbdim,klev)    !< cloud cover
  REAL(dp), INTENT(IN) :: paphm1(kbdim,klev+1) !< pressure at level interfaces
  INTEGER, INTENT(OUT) :: labels(kbdim,klev)   ! cloud labels

  REAL(dp) :: zxi(kbdim,klev)       !< cloud ice
  REAL(dp) :: zxl(kbdim,klev)       !< cloud liquid
  REAL(dp) :: zt(kbdim,klev)        !< temperature (K)
  REAL(dp) :: zxihet(kbdim,klev)    !< tracers
  REAL(dp) :: zqioliq(kbdim,klev)   !< tracers

  REAL(dp) :: zlf(kbdim,klev)       !< liquid fraction
  REAL(dp) :: zhetf(kbdim,klev)     !< het. formed ice fraction
  REAL(dp) :: zliqof(kbdim,klev)    !< liquid origin ice fraction (1)

  REAL(dp) :: predictors(kbdim,klev,npredictors) ! predictors for cloud classification

  INTEGER  :: ktdia = 1             !< quick hack, is it ever not 1?
  INTEGER  :: klevp1                !< quick hack, is it ever not klev+1?

  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1_2d
  LOGICAL,  DIMENSION(kbdim,klev) :: ll1_2d
  LOGICAL,  DIMENSION(kbdim,klev) :: ll_cc

! --------------------------------CLOUD TYPES-----------------------------------
  REAL(dp) :: zdp_cld(kbdim,klev) ! cloud pressure difference (Pa)
  REAL(dp) :: zt_top(kbdim,klev)  ! cloud top temperature (K)
  REAL(dp) :: zp_top(kbdim,klev)  ! cloud top pressure (Pa)
  INTEGER  :: zalltop(kbdim,klev) ! cloud top index
  INTEGER  :: zallbas(kbdim,klev) ! cloud base index
  INTEGER  :: itop(kbdim,klev),         & !< flag for cloud tops
              ibas(kbdim,klev),         & !< flag for cloud bases
              icl_minusbas(kbdim,klev), & !< flag for all cloud levels excepted their base
              icl_minustop(kbdim,klev), & !< flag for all cloud levels excepted their top (useless for now)
              iclnb(kbdim)                !< number of clouds per column

  klevp1 = klev + 1

  ! copy data from global memory, do not touch global data after this
  zxi(1:kproma,:) = xim1(1:kproma,:,krow)
  zxl(1:kproma,:) = xlm1(1:kproma,:,krow)
  zt(1:kproma,:) = tm1(1:kproma,:,krow)
  zxihet(1:kproma,:) = xtm1(1:kproma,:,idt_qihet,krow)
  zqioliq(1:kproma,:) = xtm1(1:kproma,:,idt_qioliq,krow)

  !-------------------------------------------------------------------------------------
  !                                      CLOUD TYPES
  !-------------------------------------------------------------------------------------
  ll1_2d(1:kproma,:) = (zxl(1:kproma,:)+zxi(1:kproma,:)) > cqtmin
  ztmp1_2d(1:kproma,:) = MERGE(paclc(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))

  ! cloud flag
  ll_cc(1:kproma,:) = ztmp1_2d(1:kproma,:) > clc_min
  
  CALL get_cloud_bounds( &
          !-- IN
          kproma, kbdim, ktdia, klev, ztmp1_2d, &
          !- OUT
          itop, ibas, icl_minustop, icl_minusbas, iclnb)

  CALL cloud_type_helper(&
          !--IN
          kproma, kbdim, ktdia, klev, klevp1, &
          ztmp1_2d, itop, ibas, zt, paphm1, &
          !--OUT
          zdp_cld, zt_top, zp_top, zalltop, zallbas)

  !-------------------------------------------------------------------------------------
  !                                  UTILITY VARIABLES
  !-------------------------------------------------------------------------------------

  ! liquid fraction
  ztmp1_2d(1:kproma,:) = zxi(1:kproma,:) + zxl(1:kproma,:)
  zlf(1:kproma,:) = MERGE(zxl(1:kproma,:)/ztmp1_2d(1:kproma,:), 0._dp, &
                          ztmp1_2d(1:kproma,:) > cqtmin)
  zlf(1:kproma,:) = MAX(MIN(zlf(1:kproma,:), 1._dp), 0._dp)

  ! heterogeneously formed ice fraction
  zhetf(1:kproma,:) = MERGE(zxihet(1:kproma,:)/zxi(1:kproma,:), 0._dp, &
                            zxi(1:kproma,:) > cqtmin)
  zhetf(1:kproma,:) = MAX(MIN(zhetf(1:kproma,:), 1._dp), 0._dp)

  ! liquid origin ice
  zliqof(1:kproma,:) = MERGE(zqioliq(1:kproma,:)/zxi(1:kproma,:), 0._dp, &
                             zxi(1:kproma,:) > cqtmin)
  zliqof(1:kproma,:) = MAX(MIN(zliqof(1:kproma,:), 1._dp), 0._dp)  

  !-------------------------------------------------------------------------------------
  !                                   CLOUD CLASSIFIER   
  !-------------------------------------------------------------------------------------

  ! remap the predictors onto [0,1] for cluster classification
#ifdef _kmeans_types
  predictors(:,:,1) = remap(kproma, kbdim, klev, 0._dp, 1._dp, .FALSE., zlf)             ! liquid fraction
  predictors(:,:,2) = remap(kproma, kbdim, klev, 0._dp, 1._dp, .FALSE., zhetf)           ! het. formed ice fraction
  predictors(:,:,3) = remap(kproma, kbdim, klev, 0._dp, 90000._dp, .FALSE., zdp_cld)     ! cloud thickness
  predictors(:,:,4) = remap(kproma, kbdim, klev, 180._dp, 330._dp, .FALSE., zt_top)      ! cloud top temperature
  predictors(:,:,5) = remap(kproma, kbdim, klev, 1.e-10_dp, 1.e-2_dp, .TRUE., zxi + zxl) ! total condensate

  labels = classify(kproma, kbdim, klev, predictors) ! use cluster classification
#else
  labels = classify_phys(kproma, kbdim, klev, zlf, zhetf, zliqof, zdp_cld, zt_top) ! use masked classification
#endif

  ! set label range
  labels(1:kproma,:) = MERGE(labels(1:kproma,:), 0, ll_cc(1:kproma,:))
  labels(1:kproma,:) = MERGE(labels(1:kproma,:), 0, zt(1:kproma,:) > 180._dp)
  labels(1:kproma,:) = MERGE(labels(1:kproma,:), 0, zt(1:kproma,:) < 330._dp)

END SUBROUTINE compute_cloud_labels

! aggregate labels along columns to mimick the analysis in
! Matus&Ecuyer 2017.
SUBROUTINE aggregate_labels_matus(&
              !--IN
              kproma, kbdim, klev, labels, paph, &
              !--OUT
              labels_matus)

  USE mo_physical_constants, ONLY: grav
  
  INTEGER, INTENT(IN) :: kproma, klev, kbdim
  INTEGER, INTENT(IN) :: labels(kbdim,klev)
  REAL(dp), INTENT(IN) :: paph(kbdim,klev+1)

  INTEGER, INTENT(OUT) :: labels_matus(kbdim,klev) ! 


  REAL(dp), DIMENSION(kbdim,klev) :: zdp, zdpg
  REAL(dp), DIMENSION(kbdim) :: zclmass

  LOGICAL, DIMENSION(kbdim,klev) :: ll_cc
  LOGICAL, DIMENSION(kbdim,klev) :: ll1_2d, ll2_2d, ll3_2d
  LOGICAL, DIMENSION(kbdim) :: ll1, ll2, ll3
  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1_2d, ztmp2_2d, ztmp3_2d
  REAL(dp), DIMENSION(kbdim) :: ztmp1, ztmp2, ztmp3

  INTEGER  :: itop(kbdim,klev),         & !< flag for cloud tops
              ibas(kbdim,klev),         & !< flag for cloud bases
              icl_minusbas(kbdim,klev), & !< flag for all cloud levels excepted their base
              icl_minustop(kbdim,klev), & !< flag for all cloud levels excepted their top
              iclnb(kbdim)                !< number of clouds per column

  ! cloud flag:
  ll_cc(1:kproma,:) = labels(1:kproma,:) > 0._dp

  ztmp1_2d(1:kproma,:) = MERGE(1._dp, 0._dp, ll_cc(1:kproma,:))
  CALL get_cloud_bounds( &
          !-- IN
          kproma, kbdim, 1, klev, ztmp1_2d, &
          !- OUT
          itop, ibas, icl_minustop, icl_minusbas, iclnb)
  
  !Pressure differences:
  zdp(1:kproma,1:klev)  = paph(1:kproma,1+1:klev+1) - paph(1:kproma,1:klev)
  zdpg(1:kproma,1:klev) = zdp(1:kproma,1:klev) / grav

  ! total column-integrated cloud mass
  ztmp1_2d(1:kproma,:) = MERGE(1._dp, 0._dp, ll_cc(1:kproma,:))
  zclmass(1:kproma) = SUM(ztmp1_2d(1:kproma,:)*zdpg(1:kproma,:), 2)

  ! 4 (in-situ cirrus), 5 (liquid origin cirrus)
  ! Note that thick clouds are probably counted as multi-layered clouds in Matus
  ztmp1_2d(1:kproma,:) = MERGE(1._dp, 0._dp, labels(1:kproma,:) == 4 .OR. labels(1:kproma,:) == 5)
  ! compute relative mass of ice clouds
  ztmp1(1:kproma) = SUM(zdpg(1:kproma,:)*ztmp1_2d(1:kproma,:), 2)
  ztmp1(1:kproma) = MERGE(ztmp1(1:kproma)/zclmass(1:kproma), 0._dp, zclmass(1:kproma) > cqtmin)

  ! liquid labels are 2 (low) and 3 (high) liquid clouds
  ztmp2_2d(1:kproma,:) = MERGE(1._dp, 0._dp, labels(1:kproma,:) == 2 .OR. labels(1:kproma,:) == 3)
  ! compute relative mass of liquid clouds
  ztmp2(1:kproma) = SUM(zdpg(1:kproma,:)*ztmp2_2d(1:kproma,:), 2)
  ztmp2(1:kproma) = MERGE(ztmp2(1:kproma)/zclmass(1:kproma), 0._dp, zclmass(1:kproma) > cqtmin)

  ! mixed-phase labels are 6 (liquid dominated) and 7 (ice dominated)
  ztmp3_2d(1:kproma,:) = MERGE(1._dp, 0._dp, labels(1:kproma,:) == 6 .OR. labels(1:kproma,:) == 7)
  ! compute relative mass of mixed-phase clouds
  ztmp3(1:kproma) = SUM(zdpg(1:kproma,:)*ztmp3_2d(1:kproma,:), 2)
  ztmp3(1:kproma) = MERGE(ztmp3(1:kproma)/zclmass(1:kproma), 0._dp, zclmass(1:kproma) > cqtmin)

  ! cases:
  ! 1) there is only 1 cloud -> take the dominant phase (rel. weight > 0.5)
  ! 2) there are multiple clouds -> if rel. weight > 0.9, take this label
  !                                 else: multi-layer cloud
!>>DN bugfix
!  ll1(1:kproma) = iclnb == 1 .AND. ztmp1(1:kproma) > 0.5_dp
!  ll1(1:kproma) = ll1(1:kproma) .OR. (iclnb > 1  .AND. ztmp1(1:kproma) > 0.9_dp)
  ll1(1:kproma) = (iclnb(1:kproma) == 1 .AND. ztmp1(1:kproma) > 0.5_dp)
  ll1(1:kproma) = (ll1(1:kproma) .OR. (iclnb(1:kproma) > 1  .AND. ztmp1(1:kproma) > 0.9_dp))
!<<DN bugfix
  ll1_2d(1:kproma,:) = SPREAD(ll1(1:kproma), 2, klev)

!>>DN bugfix
!  ll2(1:kproma) = iclnb == 1 .AND. ztmp2(1:kproma) > 0.5_dp
!  ll2(1:kproma) = ll2(1:kproma) .OR. (iclnb > 1  .AND. ztmp2(1:kproma) > 0.9_dp)
  ll2(1:kproma) = (iclnb(1:kproma) == 1 .AND. ztmp2(1:kproma) > 0.5_dp)
  ll2(1:kproma) = (ll2(1:kproma) .OR. (iclnb(1:kproma) > 1  .AND. ztmp2(1:kproma) > 0.9_dp))
!<<DN bugfix
  ll2_2d(1:kproma,:) = SPREAD(ll2(1:kproma), 2, klev)

!>>DN bugfix
!  ll3(1:kproma) = iclnb == 1 .AND. ztmp3(1:kproma) > 0.5_dp
!  ll3(1:kproma) = ll3(1:kproma) .OR. (iclnb > 1  .AND. ztmp3(1:kproma) > 0.9_dp)
  ll3(1:kproma) = (iclnb(1:kproma) == 1 .AND. ztmp3(1:kproma) > 0.5_dp)
  ll3(1:kproma) = (ll3(1:kproma) .OR. (iclnb(1:kproma) > 1  .AND. ztmp3(1:kproma) > 0.9_dp))
!<<DN bugfix
    ll3_2d(1:kproma,:) = SPREAD(ll3(1:kproma), 2, klev)

  ! initialize labels to zero
  labels_matus(1:kproma,:) = 0

  ! aggregate by phase
  labels_matus(1:kproma,:) = MERGE(1, labels_matus(1:kproma,:), ll1_2d(1:kproma,:)) ! ice clouds
  labels_matus(1:kproma,:) = MERGE(2, labels_matus(1:kproma,:), ll2_2d(1:kproma,:)) ! liquid clouds
  labels_matus(1:kproma,:) = MERGE(3, labels_matus(1:kproma,:), ll3_2d(1:kproma,:)) ! mixed-phase clouds
  labels_matus(1:kproma,:) = MERGE(labels_matus(1:kproma,:), 4, &                   ! multi-layered clouds
                                   ll1_2d(1:kproma,:) .OR. ll2_2d(1:kproma,:) .OR. ll3_2d(1:kproma,:))

  ! apply same restrictions as for original labels
  labels_matus(1:kproma,:) = MERGE(labels_matus(1:kproma,:), 0, ll_cc(1:kproma,:))

END SUBROUTINE aggregate_labels_matus

! truncate the field to [fmin,fmax] and map to [0,1]
FUNCTION remap(kproma, kbdim, klev, fmin, fmax, llog, field) RESULT(field_r)

  INTEGER :: kproma, kbdim, klev
  LOGICAL :: llog
  REAL(dp) :: fmin, fmax
  REAL(dp) :: field(kbdim,klev)

  REAL(dp) :: field_r(kbdim,klev)
  REAL(dp) :: fmin_r, fmax_r

  field_r(1:kproma,:) = MIN(MAX(field(1:kproma,:), fmin), fmax)

  ! compute log10 if necessary
  IF(llog) THEN
     field_r(1:kproma,:) = LOG10(field_r(1:kproma,:))
     fmin_r = LOG10(fmin)
     fmax_r = LOG10(fmax)
  ELSE
     field_r(1:kproma,:) = field(1:kproma,:)
     fmin_r = fmin
     fmax_r = fmax
  ENDIF

  ! remap THE FIELD TO [0,1]
  field_r(1:kproma,:) = (field_r(1:kproma,:) - fmin_r)/(fmax_r - fmin_r)

END FUNCTION remap

! takes predictors on [0,1]
FUNCTION classify(kproma, kbdim, klev, predictors) RESULT(labels)

  INTEGER :: kproma, kbdim, klev
  REAL(dp) :: predictors(kbdim,klev,npredictors)
  INTEGER :: labels(kbdim,klev)

  REAL(dp) :: distance(kbdim,klev,nclusters)
  INTEGER :: jl, jk, nclu
  REAL(dp) :: dist

  DO nclu=1,nclusters
     DO jl=1,kproma
        DO jk=1,klev
           ! compute Euclidian distance to centers in predictor-space
           distance(jl,jk,nclu) = SUM((centers(nclu,:) - predictors(jl,jk,:))**2)
        END DO !jk
     END DO !jl
  END DO !nclu
  distance(1:kproma,:,:) = SQRT(distance(1:kproma,:,:))

  ! Get the center which has minimal distance to the points
  ! Minloc returns the index to the smalles element along dimension (-> clusters)
  labels(1:kproma,:) = MINLOC(distance(1:kproma,:,:), 3)

END FUNCTION classify

! takes predictors in physical space as it is not based on distances
FUNCTION classify_phys(kproma, kbdim, klev, plf, phetf, pliqof, pdp, pttop) RESULT(labels)
  USE mo_echam_cloud_params, ONLY: cthomi
  USE mo_physical_constants, ONLY: tmelt

  INTEGER :: kproma, kbdim, klev

  REAL(dp) :: plf(kbdim,klev)           ! liquid fraction
  REAL(dp) :: phetf(kbdim,klev)         ! heterogeneously formed ice fraction
  REAL(dp) :: pliqof(kbdim,klev)        ! liquid origin ice fraction
  REAL(dp) :: pdp(kbdim,klev)           ! cloud thickness (delta p) (Pa)
  REAL(dp) :: pttop(kbdim,klev)         ! cloud top temperature (K)
  INTEGER  :: labels(kbdim,klev)

  LOGICAL :: thick_mask(kbdim,klev)    ! mask for thick clouds
  LOGICAL :: lliq_mask(kbdim,klev)      ! mask for liquid clouds
  LOGICAL :: hliq_mask(kbdim,klev)      ! mask for liquid clouds
  LOGICAL :: lcirrus_mask(kbdim,klev)  ! mask for cirrus+altostratus clouds (liquid origin)
  LOGICAL :: icirrus_mask(kbdim,klev)  ! mask for cirrus+altostratus clouds (in-situ)
  LOGICAL :: mxliq_mask(kbdim,klev)    ! mask for mixed-phase (liquid dominated) clouds
  LOGICAL :: mxice_mask(kbdim,klev)    ! mask for mixed-phase (ice dominated) clouds

  REAL(dp) :: dpthresh = 50000

  ! ------------ THICK CLOUDS ------------
  thick_mask(1:kproma,:)  = pdp(1:kproma,:) > dpthresh  .AND. &
                            phetf(1:kproma,:) < 0.5_dp

  ! ------------ LOW TOP LIQUID CLOUDS -----------
  lliq_mask(1:kproma,:)    = plf(1:kproma,:) > 0.5_dp   .AND. &
                            phetf(1:kproma,:) < 0.5_dp .AND. &
                            pdp(1:kproma,:) < dpthresh .AND. &
                            pttop(1:kproma,:) > tmelt

  ! ------------ HIGH TOP LIQUID CLOUDS -----------
  hliq_mask(1:kproma,:)    = plf(1:kproma,:) > 0.5_dp   .AND. &
                            phetf(1:kproma,:) < 0.5_dp .AND. &
                            pdp(1:kproma,:) < dpthresh .AND. &
                            pttop(1:kproma,:) < tmelt

  ! ------------ IN-SITU CIRRUS CLOUDS -----------
  icirrus_mask(1:kproma,:) = plf(1:kproma,:) < 0.5_dp   .AND. &
                            pdp(1:kproma,:) < dpthresh    .AND. &
                            phetf(1:kproma,:) < 0.5_dp  .AND. &
                            pliqof(1:kproma,:) < 0.5_dp

  ! ------------ LIQUID ORIGIN CIRRUS CLOUDS -----------
  lcirrus_mask(1:kproma,:) = plf(1:kproma,:) < 0.5_dp   .AND. &
                            pdp(1:kproma,:) < dpthresh    .AND. &
                            phetf(1:kproma,:) < 0.5_dp  .AND. &
                            pliqof(1:kproma,:) > 0.5_dp

  ! ------------ MX-LIQ CLOUDS -----------
  mxliq_mask(1:kproma,:)  = plf(1:kproma,:) > 0.5_dp   .AND. &
                            phetf(1:kproma,:) > 0.5_dp

  ! ------------ MX-ICE CLOUDS -----------
  mxice_mask(1:kproma,:)  = plf(1:kproma,:) < 0.5_dp   .AND. &
                            phetf(1:kproma,:) > 0.5_dp

  ! -------------- UN CLOUDS -------------
  ! the above parameter combination is a complete set
  labels(1:kproma,:) = 0

  ! store labels note that they are pairwise exclusive
  labels(1:kproma,:) = MERGE(1, labels(1:kproma,:), thick_mask(1:kproma,:))
  labels(1:kproma,:) = MERGE(2, labels(1:kproma,:), lliq_mask(1:kproma,:))
  labels(1:kproma,:) = MERGE(3, labels(1:kproma,:), hliq_mask(1:kproma,:))
  labels(1:kproma,:) = MERGE(4, labels(1:kproma,:), icirrus_mask(1:kproma,:))
  labels(1:kproma,:) = MERGE(5, labels(1:kproma,:), lcirrus_mask(1:kproma,:))
  labels(1:kproma,:) = MERGE(6, labels(1:kproma,:), mxliq_mask(1:kproma,:))
  labels(1:kproma,:) = MERGE(7, labels(1:kproma,:), mxice_mask(1:kproma,:))

END FUNCTION classify_phys


END MODULE mo_cloudtypes
