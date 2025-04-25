#define _read_fall_velocity read_my_fall_velocity
#define _secondary_ice_properties_p3 my_secondary_ice_properties_p3
MODULE mo_cloud_micro_p3

  USE mo_kind,               ONLY : dp
  USE mo_math_constants,     ONLY : pi
  USE mo_control,            ONLY : lcolumn, ltimer
  USE mo_physical_constants, ONLY : cpd, cpv, grav, rgrav, rd, alv, als, rv,   &
                                    vtmpc1, vtmpc2, rhoh2o, ak, tmelt, p0sl_bg, con0_h
  USE mo_time_control,   ONLY : delta_time, time_step_len
  USE mo_param_switches, ONLY : icover, nauto, lvdiff, &   !++mgs
                                ncd_activ, nic_cirrus, lrad, &
#ifdef HAMMOZ
                                lorocirrus, &
#endif
                                lsecprod, lconv
  USE mo_echam_cloud_params, ONLY : cqtmin, cvtfall, crhosno, cn0s, ccwmin      &
                                  , cthomi,  clmax, clmin, jbmin, jbmax, lonacc &
                                  , ccraut, ceffmin, ceffmax, crhoi, ccsaut, csecfrl
  USE mo_cloud_utils,    ONLY: epsec, xsec, qsec, eps, mi, &
                               ri_vol_mean_1, ri_vol_mean_2, &
                               alfased_1, alfased_2, alfased_3, &
                               betased_1, betased_2, betased_3, &
                               icemin, rcd_vol_max, &
                               cdi, mw0, mi0, mi0_rcp, ka, kb, &
                               alpha, xmw, fall, rhoice, conv_effr2mvr, clc_min, icemax, &
                               dw0, exm1_1, exp_1, exm1_2, &
                               exp_2, pirho_rcp, cap, cons4, cons5, &
                               fact_PK, pow_PK, & !SF
                               get_util_var, get_cloud_bounds, fact_coll_eff, fact_tke, &
                               get_precip_fraction, &
                               effective_ice_radius, effective_liquid_radius, &
                               riming_variables, &
                               eii, eci, rho_rim, rho_frz
                               
  USE mo_echam_convect_tables, ONLY: lookuperror, jptlucu1    &
                                   , jptlucu2, tlucua, tlucub, tlucuaw, fdeltat

  USE mo_p3_fields,      ONLY : p3isalloc, idt_qirim, idt_birim, idt_qihet, idt_qiliq, idt_qioliq, &
                                idt_qsrc, idt_qprc, &
                                idt_nihet, idt_nihom, idt_ninuc, idt_nidet, &
                                access_lookup_table, calculate_lookup_table_indices, &
                                calculate_lookup_table_indices_1d, &
                                calculate_my_lookup_table_indices, access_my_lookup_table, &
                                get_ni_limits, &
                                ! CONTROL
                                iprog, lpiggy, l2moment, ldisableupdraftcond, lconstnc, &
                                ! ACTIVATION
                                lprescribeaerosols, lprescribeaerosols_every,                  &
                                lprescribe_numbers, lprescribe_numbers_every, &
                                ! SEDIMENTATION
                                lsnow_sed, lfalling_ice, lrain_sed, vt_rdc,                     &
                                ! CLOUD COVER
                                ccclpwr, lallow_si,                                            &
                                ! PROCESS RATES
                                lact, lctrl_frz_below_238K, lctrl_het_mxphase_frz, levap_rain, &
                                lmelting, lriming, lself_collection, ladjustment,              &
                                inumberadjustment, lcover_adjustment, lmass_transport,         &
                                lsubsat_cnd_dep, lcirrus, lwbf, isublimation, lconvice,        &
                                ihomfrz,                                                       &
                                ! PRESCRIBE
                                nconstnc, nximult, nxlmult, nmicro, nmicro_max, nsedi, csubw,  &
                                iintscheme, ldiabheat, microprct,                              &
                                ! MASS MINIMUM
                                xismall, xlsmall,                                              &
                                ! TUNING
                                rcmax, adjsupsat, ccnislf, ccqccol, ccftau, ccrsnow

  USE mo_preaero,        ONLY : prescribe_aerosols_mpace
  USE mo_timer,          ONLY : timer_start, timer_stop, timer_sed, timer_micro
  USE mo_debugs

  IMPLICIT NONE

  PUBLIC :: cloud_micro_interface_p3

  INTEGER :: ibc_cvcbot, ibc_wcape, ibc_tconv, ibc_detr_cond !< indices for vars passed with the 
                                                             !< boundary condition scheme

  LOGICAL :: ll_het !SF now set by cloud_subm_1

  REAL(dp) ::  zcons1, zcons2, zcons3      , & !< various constants
               zdtime, zdtime_rcp, ztmst, ztmst_rcp,ztmstdt, zdt !< various delta t's and reciprocal

  INTERFACE set_lookup_index
     MODULE PROCEDURE set_lookup_index_1d
     MODULE PROCEDURE set_lookup_index_2d
  END INTERFACE set_lookup_index

  INTERFACE update_saturation_values
     MODULE PROCEDURE update_saturation_values_1d
     MODULE PROCEDURE update_saturation_values_2d
  END INTERFACE update_saturation_values

  INTERFACE sat_spec_hum
     MODULE PROCEDURE sat_spec_hum_1d
     MODULE PROCEDURE sat_spec_hum_2d
  END INTERFACE sat_spec_hum

  INTERFACE secondary_prognostics
     MODULE PROCEDURE secondary_prognostics_1d
     MODULE PROCEDURE secondary_prognostics_2d
  END INTERFACE secondary_prognostics

  INTERFACE subtimestep
     MODULE PROCEDURE subtimestep_1d
     MODULE PROCEDURE subtimestep_2d
  END INTERFACE subtimestep

  INTERFACE threshold_vert_vel
     MODULE PROCEDURE threshold_vert_vel_1d
     MODULE PROCEDURE threshold_vert_vel_2d
  END INTERFACE threshold_vert_vel

  INTERFACE effective_2_volmean_radius_param_Schuman_2011
     MODULE PROCEDURE effective_2_volmean_radius_param_Schuman_2011_1d
     MODULE PROCEDURE effective_2_volmean_radius_param_Schuman_2011_2d
  END INTERFACE effective_2_volmean_radius_param_Schuman_2011

  INTERFACE secondary_ice_properties_2m
     MODULE PROCEDURE secondary_ice_properties_2m_1d
     MODULE PROCEDURE secondary_ice_properties_2m_2d
  END INTERFACE secondary_ice_properties_2m

  INTERFACE conservation_reduction
     MODULE PROCEDURE conservation_reduction_1d
     MODULE PROCEDURE conservation_reduction_2d
  END INTERFACE conservation_reduction

  INTERFACE ice_switch
     MODULE PROCEDURE ice_switch_1d
     MODULE PROCEDURE ice_switch_2d
  END INTERFACE ice_switch

  INTERFACE deposition_l
     MODULE PROCEDURE deposition_l_1d
     MODULE PROCEDURE deposition_l_2d
  END INTERFACE deposition_l

  CONTAINS

SUBROUTINE cloud_micro_interface_p3( &
              !-- IN:
              kproma, kbdim, klev, klevp1, ktrac, ktdia, krow, knvb, &
              paphm1, papm1, papp1, pqm1, ptm1, ptvm1, pxlm1, &
              pxim1, &
              pvervel, pgeo, pxtm1, paphp1, ptkem1, &
              !-- INOUT:
              pqtec, paclc, paclcac, &
              paclcov, paprl, pqvi, pxlvi, pxivi, pqte, ptte, &
              pxlte, pxite, pxtte, paprs, &
              !-- OUT:
              pacdnc, picnc, prelhum, pssfl, prsfl)

USE mo_vphysc,             ONLY : set_vphysc_var      !++mgs
USE mo_activ,              ONLY : swat,            &
                                  qnuc, qaut, qacc, qfre, qmel,            &
                                  cdnc_acc, lwc_acc, cloud_time, cloud_cover_duplic,&
                                  cdnc_burden_acc, reffl_acc, burden_time, &
                                  icnc_burden_acc, reffi_acc, icnc_acc,    &
                                  cliwc_time, burdic_time, iwc_acc,        &
                                  cdnc_burden, icnc_burden, cdnc, icnc,    &
                                  sice, reffl_ct, reffl_time, cdnc_ct,     &
                                  reffi_tovs, reffi_time, iwp_tovs,        &
                                  idt_cdnc, idt_icnc, nfrzmod, reffl, reffi

!>>> Remos 2m rates
USE mo_remos,              ONLY : dqsnow, dqsflx, dsacl, dsacln, dsaci,    &
                                  dsacin, dsaut, dsautn, driv_2m, &
                                  drieff_2m, dqssub, dqsmlt, dvtim_2m,     &
                                  dvtin_2m, dqissub, dnissub, dqismlt,     &
                                  dnismlt, dqisflx, dnisflx, dqifal,       &
                                  dnifal, dvtimfal, dvtinfal,              &
!<<< Remos 2m rates
!>>> Remos process rates
                                  dqiflx, dniflx, dqrimflx, dbrimflx, dfrl_het, &
                                  dfrl_hom, dfrln_hom, dfrln_het, ddep_wbf, &
                                  dqimlt, dnimlt, drpr, dqrflx, dqrevp, dxite_cv, &
                                  ddep_co, ddep_aj, dni_aj, dxlte_cv, dcnd_co, &
                                  dcnd_aj, dnc_aj, dxisub_cv, dxlevp_cv, dxisub_tp, &
                                  dxlevp_tp, dnislf, dnifrz_ci, dximlt_cv,          &
                                  dnccol, dqccol, dncnuc, dcdncact, dicnc_cv, &
                                  dcdnc_cv, ddep_ci1, ddep_ci2, dxisub_cc, dxlevp_cc,  &
                                  dcdnc_wbf, dqisub, dni_lkp, dqimlt_rain,  &
                                  drprn, dnisub, dnisub_cc, dnisub_tp,      &
                                  dncevp_tp, dncevp_cc, dnsub_co, dnevp_co, &
                                  ddeltaqf, ddeltaqs, dqcdif, dnimlt_rain,  &
                                  dqisten, dqristen, dnisten, dbgsten,      &
                                  dqhsten, dslf, daut, dacc, ddep_a,        &
                                  dxlte_ls, dxite_ls, dicncte_ls, dcdncte_ls, &
!<<< Remos process rates
!>>> Remos ice properties         
                                  dvtim, dvtin, drcm, drim, driv, drieff,   &
                                  dqirim, dbirim, dqihet, dqiliq, dqioliq,  &
                                  dnihet, dnihom, dninuc, dnidet,           &
                                  dfr, drhop, drhoice,                      &
                                  daclci, dfr_het, dfr_liq, dfr_oliq,       &
                                  dfr_sosi, dfr_socm,                       &
                                  dfr_nihet, dfr_nihom, dfr_ninuc, dfr_nidet,&
                                  dicnc, dcdnc, dxib, dxlb,                 &
                                  dicncb, dcdncb,                           &
                                  dxi_cloud, dxi_snow, dxivi_snow, dxivi_cloud, &
                                  dice_mu, dice_lam,                        &
!<<< Remos ice properties
!>>> Remos liquid properties
                                  dqr,                                      &
!<<< Remos liquid properties
!>>> Remos budgets
                                  dmete, dncdiff, dnidiff, dnisurf,          &
                                  dnidiff2, dncdiff2, dnictrl, dncctrl,      &
!<<< Remos budgets
!>>> Remos scavenging
                                  dfrain, dfevapr, dmratepr, dfsnow, dfsubls, &
                                  dmrateps, dmlwc, dmiwc, dmsnowacl, &
!<<< Remos scagenging
!>>> Remos dynamics
                                  dtke, dupdraft, dupdraftmax, drhoair,     &
                                  dcfl, ddpg, ddz, dsedtm, dsedtn,          &
                                  dnmicro, dnsedi,                          &
!<<< Remos cloud types
                                  ddp_cld, dt_top, diclnb, dclcol_time
!>>> Remos cloud types

USE mo_conv,               ONLY : cdncact_cv,     &
                                  conv_time, twc_conv
USE mo_cirrus,             ONLY : xfrzmstr
USE mo_cirrus_4_globalModelInterface, ONLY: CirrusModelInterface
USE mo_submodel_interface, ONLY: cloud_subm_1, cloud_subm_2, cloud_subm_3
USE mo_memory_g2a,         ONLY: vm1, um1
USE mo_memory_g3b,         ONLY: orostd,oromea,orogam,orothe, &
                                 oropic,oroval,orosig
#ifdef HAMMOZ
USE mo_orocirrus,          ONLY: orocirrus_w, orocirrus_cc
#endif

USE mo_boundary_condition, ONLY: bc_find, bc_apply
USE mo_time_control,       ONLY: lstart, lresume

  IMPLICIT NONE

!--- Arguments

! Input:
  INTEGER, INTENT(in)  :: kproma       !< block size
  INTEGER, INTENT(in)  :: kbdim        !< max. block size
  INTEGER, INTENT(in)  :: klev         !< number of vertical levels
  INTEGER, INTENT(in)  :: klevp1       !< number of vertical levels+1
  INTEGER, INTENT(in)  :: ktrac        !< number of tracers
  INTEGER, INTENT(in)  :: ktdia        !< highest vertical level for "diagnostics" (use??)-careful with this if /=1!
  INTEGER, INTENT(in)  :: krow         !< block number
  INTEGER, INTENT(in)  :: knvb(kbdim)  !< ???

  REAL(dp), INTENT(in)    :: paphm1  (kbdim,klevp1)  !< pressure at half levels         (t-1)
  REAL(dp), INTENT(in)    :: papm1   (kbdim,klev)    !< pressure at full levels         (t-1)
  REAL(dp), INTENT(in)    :: papp1   (kbdim,klev)    !< pressure at full levels         (t-1)
  REAL(dp), INTENT(in)    :: pqm1    (kbdim,klev)    !< specific humidity               (t-1)
  REAL(dp), INTENT(in)    :: ptm1    (kbdim,klev)    !< temperature                     (t-1)
  REAL(dp), INTENT(in)    :: ptvm1   (kbdim,klev)    !< virtual temperature             (t-1)
  REAL(dp), INTENT(in)    :: pxlm1   (kbdim,klev)    !< cloud liquid water              (t-1)
  REAL(dp), INTENT(in)    :: pxim1   (kbdim,klev)    !< cloud ice                       (t-1)
  REAL(dp), INTENT(in)    :: pvervel (kbdim,klev)    !< large scale vertical velocity [Pa s-1]
  REAL(dp), INTENT(in)    :: pgeo    (kbdim,klev)    !< geopotential height at full levels
! RD: CHANGED TO 'inout' TO ALLOW FOR AEROSOL PRESCRIPTION
  REAL(dp), INTENT(inout)    :: pxtm1   (kbdim,klev,ktrac) !< tracer mmr                    (t-1)
  REAL(dp), INTENT(in)    :: paphp1  (kbdim,klevp1)  !< pressure at half levels         (t)
  REAL(dp), INTENT(in)    :: ptkem1  (kbdim,klev)    !< turbulent kinetic energy         (t-1)

! Input / output:
  REAL(dp), INTENT(inout) :: pqtec   (kbdim,klev)    !< large-scale vertical velocity 
  REAL(dp), INTENT(inout) :: paclc   (kbdim,klev)    !< cloud cover  (now diagnosed in cover)
  REAL(dp), INTENT(inout) :: paclcac (kbdim,klev)    !< cloud cover, accumulated
  REAL(dp), INTENT(inout) :: paclcov (kbdim)         !< total cloud cover
  REAL(dp), INTENT(inout) :: paprl   (kbdim)         !< total stratiform precip. (rain+snow), accumulated
  REAL(dp), INTENT(inout) :: pqvi    (kbdim)         !< vertically integrated spec. humidity, accumulated
  REAL(dp), INTENT(inout) :: pxlvi   (kbdim)         !< vertically integrated cloud liquid water, accumul.
  REAL(dp), INTENT(inout) :: pxivi   (kbdim)         !< vertically integrated cloud ice, accumulated
  REAL(dp), INTENT(inout) :: pqte    (kbdim,klev)    !< tendency of specific humidity
  REAL(dp), INTENT(inout) :: ptte    (kbdim,klev)    !< tendency of temperature
  REAL(dp), INTENT(inout) :: pxlte   (kbdim,klev)    !< tendency of cloud liquid water
  REAL(dp), INTENT(inout) :: pxite   (kbdim,klev)    !< tendency of cloud ice
  REAL(dp), INTENT(inout) :: pxtte   (kbdim,klev,ktrac) !< tracer tendency (cdnc, icnc)
  REAL(dp), INTENT(inout) :: paprs   (kbdim)         !< snowfall accumulated

! Output:
  REAL(dp), INTENT(out) :: pacdnc  (kbdim,klev)    !< cloud droplet number concentration (specified)
  REAL(dp), INTENT(out) :: picnc   (kbdim,klev)    !< Ice crystal number concentration (ICNC) [1/m3]
  REAL(dp), INTENT(out) :: prelhum (kbdim,klev)    !< relative humidity
  REAL(dp), INTENT(out) :: pssfl   (kbdim)         !< surface snow flux (kg/m2/s)
  REAL(dp), INTENT(out) :: prsfl   (kbdim)         !< surface rain flux (kg/m2/s)

! -------------------------------------------------------------------------------- !
!                               LOCAL VARIABLES                                    !
! -------------------------------------------------------------------------------- !

! ----------------------------GRIDBOX DESCRIPTION---------------------------------
  REAL(dp):: zgeoh(kbdim,klevp1)  !< Geopotential height at half levels
  REAL(dp):: zdz(kbdim,klev)      !< layer thickness [m]
  REAL(dp):: zdp(kbdim,klev)      !< pressure difference of the layer [Pa]
  REAL(dp):: zdpg(kbdim,klev)     !< delta p over g [kg/m2]

  REAL(dp):: zaclc_tm1(kbdim,klev)       !< cloud cover at t-1 [1]
  REAL(dp):: zclcov(kbdim)               !< Total cloud cover (overlap considered)
  REAL(dp):: zclcpre(kbdim)              !< fraction of grid box covered by precip (column-wise)
  REAL(dp):: zclcfi(kbdim)               !< fraction of grid box covered by sedimenting ice (column-wise)
  REAL(dp):: zauloc(kbdim)               !< Part of the grid box allowed to participation in accretion with
                                         !< newly formed condensate
  LOGICAL :: ll_precip(kbdim)            !< traces precipitation
  LOGICAL :: ll_prcp_warm(kbdim)         !< traces warm precipitation
  INTEGER :: itop(kbdim,klev),         & !< flag for cloud tops
             ibas(kbdim,klev),         & !< flag for cloud bases
             icl_minusbas(kbdim,klev), & !< flag for all cloud levels excepted their base
             icl_minustop(kbdim,klev), & !<  flag for all cloud levels excepted their top (useless for now)
             iclbas(kbdim,klev),       & !< conv. cloud base
             itm1_look(kbdim,klev),    & !< index for temperature lookup table
             itp1_look(kbdim,klev),    & !< itm1_look + 1
             iclnb(kbdim)                !< number of clouds per column

! -----------------------------FLUID DYNAMICS-------------------------------------
  REAL(dp):: zaaa(kbdim,klev)     !< Air density correction needed for the ice crystal fall velocity
  REAL(dp):: zkair(kbdim,klev)    !< Thermal conductivity of air ???
  REAL(dp):: zviscos(kbdim,klev)  !< Dynamic viscosity of water in air
  REAL(dp):: zrho(kbdim,klev)     !< Air density [kg/m3]
  REAL(dp):: zrho_rcp(kbdim,klev) !< Inverse air density
  REAL(dp):: zrho_corr(kbdim,klev)!< Air density correstion factor (rho0/rho where rho0 air dens. at p0=1000hPa)
  REAL(dp):: zapref                      !< reference pressure [Pa]
  REAL(dp):: zlsdcp(kbdim,klev)      !< latent heat of sublimation divided by
                                     !< the specific heat at constant pressure
  REAL(dp):: zlvdcp(kbdim,klev)      !< latent heat of vaporization divided by
                                     !< the specific heat at constant pressure
  REAL(dp):: zdv(kbdim,klev)         !< ???
  REAL(dp):: zeta(kbdim,klev)        !< Variable needed for the Bergeron-Findeisen process

! ----------------------------WATER  VARIABLES-------------------------------------
  REAL(dp):: zesw(kbdim,klev)        !< Saturation vapor pressure w.r.t. water [Pa]
  REAL(dp):: zqsw(kbdim,klev)        !< Saturation specific humidity w.r.t. water (t-1) [kg/kg]
  REAL(dp):: zsupw(kbdim,klev)       !< Supersaturation w.r.t. water (0-1)
  REAL(dp):: zesi(kbdim,klev)        !< Saturation vapor pressure w.r.t. ice [Pa]
  REAL(dp):: zqsi(kbdim,klev)        !< Saturation specific humidity w.r.t. ice [kg/kg]
  REAL(dp):: zsupi(kbdim,klev)       !< Supersaturation w.r.t. water (0-1)
  REAL(dp):: zqswp1(kbdim,klev)      !< Saturation specific humidity w.r.t. water [kg/kg] at T + 0.001 K
  REAL(dp):: zqsip1(kbdim,klev)      !< Saturation specific humidity w.r.t. ice [kg/kg] at T + 0.001 K
  REAL(dp):: zqsw0(kbdim,klev)       !< Saturation specific humidity at tmelt
  REAL(dp):: zqvi(kbdim)             !< Vertically integrated specific humidity [kg/m2]
  REAL(dp):: zxlvi(kbdim)            !< Vertically integrated cloud water [kg/m2]
  REAL(dp):: zxivi(kbdim)            !< Vertically integrated cloud ice [kg/m2]
  REAL(dp):: zxirimvi(kbdim)         !< Vertically integrated cloud rimed ice [kg/m2]
! local water vars
  REAL(dp):: zcdnc(kbdim,klev)       !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp):: zxl(kbdim,klev)         !< local cloud liquid water [kg/kg]
  REAL(dp):: zxlb(kbdim,klev)        !< in-cloud liquid water mass [kg/kg]
  REAL(dp):: zcdncb(kbdim,klev)      !< in-cloud liquid water number [1/kg]
  REAL(dp):: zxi(kbdim,klev)         !< local cloud ice water [kg/kg]
  REAL(dp):: zqirim(kbdim,klev)      !< rimed ice mass [kg/kg]
  REAL(dp):: zbirim(kbdim,klev)      !< rimed ice volume [kg/kg]
  REAL(dp):: zqihet(kbdim,klev)      !< heterogeneously formed ice [kg/kg]
  REAL(dp):: zqiliq(kbdim,klev)      !< ice formed by frozen liquid [kg/kg]
  REAL(dp):: zqioliq(kbdim,klev)     !< ice initiated by frozen liquid [kg/kg]
  REAL(dp):: zqsrc(kbdim,klev)       !< accumulated water mass source [kg/kg]
  REAL(dp):: zqprc(kbdim,klev)       !< accumulated water mass precip sink [kg/kg]
  REAL(dp):: znihet(kbdim,klev)      !< number of heterogeneously nucleated crystals [1/kg]
  REAL(dp):: znidet(kbdim,klev)      !< number of detrained crystals [1/kg]
  REAL(dp):: znihom(kbdim,klev)      !< number of homogeneously frozen droplets [1/kg]
  REAL(dp):: zninuc(kbdim,klev)      !< number of nucleated crystals (cirrus scheme) [1/kg]
  REAL(dp):: zicnc(kbdim,klev)       !< ice number [1/kg]
  REAL(dp):: zxib(kbdim,klev)        !< in-cloud ice water mass [kg/kg]
  REAL(dp):: zicncb(kbdim,klev)        !< in-cloud ice water number [1/kg]
  REAL(dp):: zqirimb(kbdim,klev)     !< in-cloud value of qirim [kg/kg]
  REAL(dp):: zbirimb(kbdim,klev)     !< in-cloud value of birim [m3/kg]
  REAL(dp):: zqihetb(kbdim,klev)     !< in-cloud value of qihet [m3/kg]
  REAL(dp):: zqiliqb(kbdim,klev)     !< in-cloud value of qiliq [m3/kg]
  REAL(dp):: zqioliqb(kbdim,klev)    !< in-cloud value of qioliq [m3/kg]
  REAL(dp):: zq(kbdim,klev)          !< local humidity [kg/kg]
  REAL(dp):: zt(kbdim,klev)          !< local temperature [kg/kg]
  REAL(dp):: ztv(kbdim,klev)         !< local virtual temperature [kg/kg]
  REAL(dp):: zqr(kbdim)              !< rain water mixing ratio [kg/kg]
! ice cover
  REAL(dp):: zaclci(kbdim,klev)      !< cloud cover used for ice clouds [0,1]
  REAL(dp):: zxi_cloud(kbdim,klev)   !< cloud portion of ice
  REAL(dp):: zxi_snow(kbdim,klev)    !< snow portion of ice
! initial guess for nmicro calculation
  REAL(dp):: zxi_guess(kbdim,klev)   !< initial guess for zxi with large-scale tendencies
  REAL(dp):: zicnc_guess(kbdim,klev) !< initial guess for zicnc with large-scale tendencies
  REAL(dp):: zqirim_guess(kbdim,klev)!< initial guess for zqirim with large-scale tendencies
  REAL(dp):: zbirim_guess(kbdim,klev)!< initial guess for zbirim with large-scale tendencies
  REAL(dp):: zrhop_guess(kbdim,klev) !< initial guess for zrhop with large-scale tendencies
  REAL(dp):: zf_guess(kbdim,klev)    !< initial guess for zrimfrac with large-scale tendencies
  REAL(dp):: zvm_guess(kbdim,klev)   !< initial guess for zvtim with large-scale tendencies
  REAL(dp):: zvn_guess(kbdim,klev)   !< initial guess for zvtin with large-scale tendencies
! liquid fraction diagnostics for precip
  REAL(dp):: zxi_avg(kbdim,klev)     !< ice water content averaged over nmicro
  REAL(dp):: zxl_avg(kbdim,klev)     !< liquid water content averaged over nmicro

! ----------------------------LOCAL TENDENCIES-------------------------------------
  REAL(dp):: zxite(kbdim,klev)       !< ice water tendency
  REAL(dp):: zqirimte(kbdim,klev)    !< rimed ice tendency
  REAL(dp):: zbirimte(kbdim,klev)    !< rimed ice volume tendency
  REAL(dp):: zqihette(kbdim,klev)    !< heterogeneously formed ice tendency
  REAL(dp):: zqiliqte(kbdim,klev)    !< ice formed by frozen liquid tendency
  REAL(dp):: zqioliqte(kbdim,klev)   !< ice initiated by frozen liquid tendency
  REAL(dp):: zqsrcte(kbdim,klev)     !< water mass source tendency
  REAL(dp):: zqprcte(kbdim,klev)     !< water mass precip sink tendency
  REAL(dp):: znihette(kbdim,klev)    !< heterogeneously frozen droplet number tendency
  REAL(dp):: znidette(kbdim,klev)    !< detrained ice crystal number tendency
  REAL(dp):: znihomte(kbdim,klev)    !< homogeneously frozen droplet number tendency
  REAL(dp):: zninucte(kbdim,klev)    !< nucleated ice number tendency (cirrus-scheme)
  REAL(dp):: zxlte(kbdim,klev)       !< liquid water tendency
  REAL(dp):: zqte(kbdim,klev)        !< water vapor tendency
  REAL(dp):: ztte(kbdim,klev)        !< temperature tendency
  REAL(dp):: zcdncte(kbdim,klev)  !< cdnc tracer tendency
  REAL(dp):: zicncte(kbdim,klev)  !< icnc tracer tendency

  ! ----------------------------LOCAL LARGE-SCALE TENDENCIES-----------------------
  REAL(dp):: zxite_ls(kbdim,klev)       !< ice water tendency
  REAL(dp):: zqirimte_ls(kbdim,klev)    !< rimed ice tendency
  REAL(dp):: zbirimte_ls(kbdim,klev)    !< rimed ice volume tendency
  REAL(dp):: zqihette_ls(kbdim,klev)    !< heterogeneously formed ice tendency
  REAL(dp):: zqiliqte_ls(kbdim,klev)    !< ice formed by frozen liquid tendency
  REAL(dp):: zqioliqte_ls(kbdim,klev)   !< ice initiated by frozen liquid tendency
  REAL(dp):: zqsrcte_ls(kbdim,klev)     !< accumulated water mass source tendency, large scale
  REAL(dp):: zqprcte_ls(kbdim,klev)     !< accumulated water mass precip sink tendency, large scale
  REAL(dp):: znihette_ls(kbdim,klev)    !< heterogeneously frozen droplet number tendency
  REAL(dp):: znidette_ls(kbdim,klev)    !< detrained ice crystal number tendency
  REAL(dp):: zninucte_ls(kbdim,klev)    !< nucleated ice tendency (cirrus-scheme)
  REAL(dp):: znihomte_ls(kbdim,klev)    !< homogeneously frozen droplet number tendency
  REAL(dp):: zxlte_ls(kbdim,klev)       !< liquid water tendency
  REAL(dp):: zqte_ls(kbdim,klev)        !< water vapor tendency
  REAL(dp):: zqte_0(kbdim,klev)         !< water vapor tendency input used for cloud form.
  REAL(dp):: ztte_ls(kbdim,klev)        !< temperature tendency
  REAL(dp):: zcdncte_ls(kbdim,klev)  !< cdnc tracer tendency
  REAL(dp):: zicncte_ls(kbdim,klev)  !< icnc tracer tendency

! ----------------------------SEDIMENTATION----------------------------------------
  REAL(dp):: zqiflx(kbdim)           !< total ice mass flux out of gridbox
  REAL(dp):: zqiflx_2d(kbdim,klev)   !< total ice mass flux out of gridbox
  REAL(dp):: zniflx(kbdim)           !< ice number flux out of gridbox
  REAL(dp):: zqirimflx(kbdim)        !< rimed ice mass flux out of gridbox
  REAL(dp):: zbirimflx(kbdim)        !< riming volume flux out of gridbox

  REAL(dp):: zqisten(kbdim,klev)      !< total ice mass sedimentation tendency
  REAL(dp):: znisten(kbdim,klev)     !< ice number sedimentation tendency
  REAL(dp):: znissrc(kbdim,klev)     !< ice number sedimentation tendency
  REAL(dp):: zqristen(kbdim,klev)     !< rimed ice mass sedimentation tendency
  REAL(dp):: zbgsten(kbdim,klev)     !< riming volume sedimentation tendency
  REAL(dp):: zqhsten(kbdim,klev)     !< het. formed ice sedimentation tendency
  REAL(dp):: zqlsten(kbdim,klev)     !< ice through liquid sedimentation tendency
  REAL(dp):: zqlosten(kbdim,klev)    !< liquid origin ice sedimentation tendency
  REAL(dp):: znihetsten(kbdim,klev)    !< het. ice number sed. tendency
  REAL(dp):: znidetsten(kbdim,klev)    !< det. ice number sed. tendency
  REAL(dp):: znihomsten(kbdim,klev)    !< hom. ice number sed. tendency
  REAL(dp):: zninucsten(kbdim,klev)    !< nuc. ice number sed. tendency

  REAL(dp):: zqrflx(kbdim)           !< rain flux (kg/m2/s)
  REAL(dp):: zqrflx_bfevp(kbdim)     !< rain flux before evaporation for scavenging (kg/m2/s)
  REAL(dp):: zrpr(kbdim,klev)        !< Rain formation rate for mass   [1/m3]
  REAL(dp):: zacc(kbdim,klev)        !< Rain formation rate for mass (accretion)   [1/m3]
  REAL(dp):: zaut(kbdim,klev)        !< Rain formation rate for mass (autoconversion)   [1/m3]
  REAL(dp):: zrprn(kbdim,klev)       !< Rain formation rate for number [1/m3]
  INTEGER :: zprcp_top(kbdim,klev)   !< contains index of highest precip level for each index
  INTEGER :: zprcp_bas(kbdim,klev)   !< contains index of lowerst precip level for each index
  INTEGER :: zice_top(kbdim,klev)    !< contains index of highest ice level for each index
  INTEGER :: zice_bas(kbdim,klev)    !< contains index of lowerst ice level for each index

  REAL(dp):: zmltflx(kbdim)          !< flux of molten ice in rain

  REAL(dp):: zvmpot(kbdim,klev)   !< Potential velocity -> velocity at cloud top (mass weighted)
  REAL(dp):: zvnpot(kbdim,klev)   !< Potential velocity -> velocity at cloud top (number weighted)
  REAL(dp):: dumm1(kbdim)            !< helper varliables to calculate potential velocities
  REAL(dp):: dumn1(kbdim)
  REAL(dp):: dumm2(kbdim)
  REAL(dp):: dumn2(kbdim)

  REAL(dp):: zsedtm(kbdim,klev)       !< time of sedimentation mass (substep-time)
  REAL(dp):: zsedtn(kbdim,klev)       !< time of sedimentation number (substep-time)

! ----------------------------HUMIDITY LOOKUPTABLE--------------------------------
  REAL(dp):: zlucua(kbdim,klev)      !< Temporary variable needed for the calculation of the sat. vapor pressure (t-1)
  REAL(dp):: zlucuaw(kbdim,klev)     !< Temporary variable needed for the calculation of the sat. vapor pressure (t-1)
  REAL(dp):: zlucuawp1(kbdim,klev)   !< Temporary variable needed for the calculation of the sat. vapor pressure (t)
  REAL(dp):: zlucuap1(kbdim,klev)    !< Temporary variable needed for the calculation of the sat. vapor pressure (t

! ----------------------------CONVECTION VARIABLES--------------------------------
  REAL(dp) :: zcvcbot (kbdim)             !< conv. cloud base index
  REAL(dp) :: zwcape  (kbdim)             !< CAPE contrib. to conv. vertical velocity [m s- 1]
  REAL(dp) :: ztconv  (kbdim,klev)        !< temperature as it was in the conv scheme [K]
  REAL(dp) :: zvervx(kbdim,klev)          !< Updraft velocity [cm/s]
  REAL(dp) :: zxtec(kbdim,klev)           !< tendency for detr. conv. cloud liq water or cloud ice (t)



! -----------------------------AEROSOLS-------------------------------------------
! activation:
  REAL(dp) :: zcdncact(kbdim,klev)     !< Number of newly activated cloud droplets
! mixed-phase freezing (despite their names, these quantities ARE NOT dependent on HAM):
  REAL(dp) :: zrwetki(kbdim,klev)  ! wet radius, aitken insoluble mode
  REAL(dp) :: zrwetai(kbdim,klev)  ! wet radius, accumulation insoluble mode
  REAL(dp) :: zrwetci(kbdim,klev)  ! wet radius, coarse insoluble mode
  REAL(dp) :: zfracdusol(kbdim,klev)   ! total fraction of dust particules in all soluble mo
  REAL(dp) :: zfracduai(kbdim,klev)    ! fraction of dust particules in the accum. soluble m
  REAL(dp) :: zfracduci(kbdim,klev)    ! fraction of dust particules in the coarse soluble m
  REAL(dp) :: zfracbcsol(kbdim,klev)   ! total fraction of BC particules in all soluble mode
  REAL(dp) :: zfracbcinsol(kbdim,klev) ! total fraction of BC particules in all insoluble mo
! cirrus freezing:
  REAL(dp) :: zascs(kbdim,klev)           !< soluble aerosol number conc.
  REAL(dp) :: zapnx(kbdim,klev,nfrzmod)   !< aerosol number available for freezing per mode [1/m3]
  REAL(dp) :: zap(kbdim,klev)             !< aerosol number available for freezing total [1/m3]
  REAL(dp) :: zaprx(kbdim,klev,nfrzmod)   !< radius of aerosols avail. for freezing  [cm]
  REAL(dp) :: zapsigx(kbdim,klev,nfrzmod) !< std. dev. of aerosols available for freezing
  REAL(dp) :: zncnuc_bas(kbdim,klev)      !< Nucleated CDNC at base of stratiform clouds
  REAL(dp) :: zcdnc_bas(kbdim,klev)       !< CDNC at base of stratiform clouds

! -----------------------------PARTICLE PROPERTIES--------------------------------
  REAL(dp) :: zrim(kbdim,klev)            !< Mass-weighted mean ice radius
  REAL(dp) :: zrieff(kbdim,klev)          !< effective ice crystal radius
  REAL(dp) :: zrin(kbdim,klev)            !< Number-weighted effective ice crystal radius
  REAL(dp) :: zriv(kbdim,klev)            !< volume mean ice radius (from zrieff)
  REAL(dp) :: zvtim(kbdim,klev)           !< mass-weighted terminal velocity of ice
  REAL(dp) :: zvtin(kbdim,klev)           !< number-weighted terminal velocity of ice
  REAL(dp) :: zri(kbdim, klev)            !< size of newly nucleated IC from cirrus scheme [m]
  REAL(dp) :: zrcm(kbdim,klev)            !< mass-weighted mean radius cloud droplets

! -----------------------------2M - PARTICLE PROPERTIES----------------------------
  REAL(dp) :: zrieff_2m(kbdim,klev)       !< effective ice crystal radius for 2M scheme
  REAL(dp) :: zriv_2m(kbdim,klev)         !< volume mean ice crystal radius for 2M scheme
  REAL(dp) :: zvtim_2m(kbdim,klev)        !< mass-weighted terminal velocity of ice
  REAL(dp) :: zvtin_2m(kbdim,klev)        !< number-weighted terminal velocity of ice

! ----------------------------------DYNAMICS--------------------------------------
  REAL(dp) :: zvervmax(kbdim,klev)        !< Threshold vertical velocity
  REAL(dp), DIMENSION(kbdim) :: zb2, zgtp, zvth, zfuchs, zfre, zre !< for sublimation


! ---------------------------------UTILITIES--------------------------------------
  REAL(dp):: zastbstw(kbdim,klev)    !< Thermodynamic term needed for water nucleation
  REAL(dp):: zastbsti(kbdim,klev)    !< Thermodynamic term needed for ice nucleation
  REAL(dp):: zmix(kbdim,klev)        !< fraction to determine supersaturation in mixed
                                     !< phase supersaturation for cnd/dep to be consistent
                                     !< with the cover scheme
  REAL(dp):: zmixp1(kbdim,klev)      !< zmix at t = t + fdeltat (the spacing in lkp table for
                                     !< clausius-clapeyron
  REAL(dp):: zcc_flag(kbdim,klev)    !< cloudy levels (xi+xl > 0, clc > clc_min)
  REAL(dp):: zcci_flag(kbdim,klev)   !< icy levels (xi > 0, clc > clc_min)
  REAL(dp):: zccl_flag(kbdim,klev)   !< liquid levels (xl > 0, clc > clc_min)
  REAL(dp):: zccil_flag(kbdim,klev)  !< ice-liquid levels (xl > 0, xi > 0, clc > clc_min)

! ----------------------------TEMPORARY VARIABLES---------------------------------
  LOGICAL, DIMENSION(kbdim,klev) :: ll1_2d, ll2_2d, ll3_2d          !< 2d logical variable
  LOGICAL, DIMENSION(kbdim,klev) :: llt_2d                          !< 2d logical variable
  LOGICAL, DIMENSION(kbdim,klev) :: lo2_2d                          !< 2d logical variable
  LOGICAL, DIMENSION(kbdim)      :: ll1, ll2, ll3                   !< 1d logical variable
  LOGICAL, DIMENSION(kbdim,klev) :: ll_cc, ll_cci                   !< cloud cover flag
  LOGICAL, DIMENSION(kbdim,klev) :: ll_cv                           !< conv. ice form. possible
  LOGICAL, DIMENSION(kbdim,klev) :: ll_ice                          !< ice formation possible
  REAL(dp), DIMENSION(kbdim)     :: ztmp1, ztmp2, ztmp3, ztmp4, zsrcs, zsnks
  REAL(dp), DIMENSION(kbdim,klev):: ztmp1_2d, ztmp2_2d, ztmp3_2d, ztmp4_2d, zsrcs_2d, zsnks_2d  
  REAL(dp) :: zcd2ic(kbdim) ! amount of cloud droplets erroneously present
                            ! at temp < cthomi that should be credited to ice crystal number concentration
  REAL(dp) :: zic2cd(kbdim) ! amount of ice crystals erroneously present
                            ! at temp > tmelt that should be credited to cloud droplet number concentration
  REAL(dp) :: zcd2unphys(kbdim)      !< amount of cloud droplets erroneously present at temp < cthomi that
                                     !< are unphysical
  REAL(dp) :: zcd2ic_2d(kbdim,klev)       !< same as zcd2ic, 2D version
  REAL(dp) :: zcd2unphys_2d(kbdim,klev)   !< same as zcd2unphys, 2D version
  INTEGER :: jl, jk, jkk, jfrzmod, iqidx, ixidx, ierr, jkkk
  REAL(dp) :: zratio(kbdim)          !< ratio used for process rate adjustments
  REAL(dp) :: zrid(kbdim,klev)       !< Schumann size assumed for particles from conv. updraft

! ----------------------------------CONSTANTS-----------------------------------
  REAL(dp), PARAMETER :: zdummy = 1._dp !SF this is a dummy var to handle special non-physical cases 
                                        !   in a MERGE statement. This var is solely needed to prevent
                                        !   spurious calculations like div by zero for example
  LOGICAL, PARAMETER :: nosize = .true. !< .true. ---> aerosol size effects are ignored for homogeneous freezing
  INTEGER  :: zntmst                    !< number of subtimesteps in column_processes loop
  REAL(dp) :: zstmst                    !< subtimestep length
  REAL(dp) :: zstmst_rcp                !< subtimestep length reciprocal
  REAL(dp) :: zntmst_rcp                !< number of subtimesteps reciprocal
  INTEGER  :: nstep                     !< substep running variable
  INTEGER  :: znsedi                    !< substep running variable


! -----------------------------PROCESS RATES------------------------------------
! Nucleation:
  REAL(dp):: zncnuc(kbdim,klev)      !< Nucleation rate of CDNC [1/m3]
  REAL(dp):: zcdnc_cv(kbdim,klev)    !< Nucleation rate of CDNC by convection [1/m3] (immediately updated)
  REAL(dp):: zicnc_cv(kbdim,klev)    !< Nucleation rate of ICNC by convection [1/m3] (immediately updated)
  REAL(dp):: znifrz_ci(kbdim,klev)   !< ICNC by freezing in cirrus scheme (xfrzmstr)
  REAL(dp):: zdep_ci1(kbdim,klev)    !< Deposition rate from cirrus scheme directly
  REAL(dp):: zdep_ci2(kbdim)         !< Deposition rate in cirrus regime by exlicit deposition
! entrainment
  REAL(dp):: zxite_cv(kbdim,klev)    !< convective ice tendency
  REAL(dp):: zxlte_cv(kbdim,klev)    !< convective liquid water tendency
  REAL(dp):: zxisub_cv(kbdim,klev)   !< sublimated convective ice mass in clear-sky part
  REAL(dp):: zxlevp_cv(kbdim,klev)   !< evaporated convective liquid water mass in clear-sky part
  REAL(dp):: zximlt_cv(kbdim,klev)   !< ice that needs to be melted to merge convective and strat. schemes
! transport
  REAL(dp):: zxisub_tp(kbdim,klev)   !< sublimated transported ice mass in clear-sky part
  REAL(dp):: znisub_tp(kbdim,klev)   !< sublimated transported ice number in clear-sky part
  REAL(dp):: zqirimsub_tp(kbdim,klev)!< sublimated transported rimed ice in clear-sky part
  REAL(dp):: zbirimsub_tp(kbdim,klev)!< sublimated transported rimed volume in clear-sky part
  REAL(dp):: zqihetsub_tp(kbdim,klev)!< sublimated transported heterogeneously formed ice in clear-sky part
  REAL(dp):: zqiliqsub_tp(kbdim,klev)!< sublimated transported frozen liquid mass in in clear-sky part
  REAL(dp):: zqioliqsub_tp(kbdim,klev)!< sublimated transported frozen liquid origin in clear-sky part
  REAL(dp):: zqsrcsub_tp(kbdim,klev) !< transport sublimation: acc. source term
  REAL(dp):: zqprcsub_tp(kbdim,klev) !< transport sublimation: acc. prcp. term
  REAL(dp):: znihetsub_tp(kbdim,klev)!< sub. transport of nihet in clear-sky part
  REAL(dp):: znidetsub_tp(kbdim,klev)!< sub. transport of nidet in clear-sky part
  REAL(dp):: znihomsub_tp(kbdim,klev)!< sub. transport of nihom in clear-sky part
  REAL(dp):: zninucsub_tp(kbdim,klev)!< sub. transport of ninuc in clear-sky part
  REAL(dp):: zxlevp_tp(kbdim,klev)   !< evaporated transported liquid water mass in clear-sky part
  REAL(dp):: zncevp_tp(kbdim,klev)   !< evaporated transported liquid water number in clear-sky part
! dissipation
  REAL(dp):: zxisub_cc(kbdim)        !< sublimated ice mass in all-clear-sky condition
  REAL(dp):: znisub_cc(kbdim)        !< sublimated ice number in all-clear-sky condition
  REAL(dp):: zqirimsub_cc(kbdim)     !< sublimated rimed ice in all-clear-sky condition
  REAL(dp):: zbirimsub_cc(kbdim)     !< sublimated rimed volume in all-clear-sky condition
  REAL(dp):: zqihetsub_cc(kbdim)     !< sublimated heterogeneously formed ice in all-clear-sky condition
  REAL(dp):: zqiliqsub_cc(kbdim)     !< sublimated frozen liquid mass in all-clear-sky condition
  REAL(dp):: zqioliqsub_cc(kbdim)    !< sublimated frozen liquid origin in all-clear-sky condition
  REAL(dp):: zqsrcsub_cc(kbdim)      !< clear-sky sub: acc. source term
  REAL(dp):: zqprcsub_cc(kbdim)      !< clear-sky sub: acc. precip sink term
  REAL(dp):: znihetsub_cc(kbdim)     !< sub. of nihet in all-clear-sky condition
  REAL(dp):: znidetsub_cc(kbdim)     !< sub. of nidet in all-clear-sky condition
  REAL(dp):: znihomsub_cc(kbdim)     !< sub. of nihom in all-clear-sky condition
  REAL(dp):: zninucsub_cc(kbdim)     !< sub. of ninuc in all-clear-sky condition
  REAL(dp):: zxlevp_cc(kbdim)        !< evaporated transported liquid water mass in all-clear-sky condition
  REAL(dp):: zncevp_cc(kbdim)        !< evaporated transported liquid water number in all-clear-sky condition
! humidity convergence (transport + entrainment)
  REAL(dp):: zdep_a(kbdim)           !< deposition capability
  REAL(dp):: zcnd_co(kbdim)          !< condensation from convergence
  REAL(dp):: znevp_co(kbdim)         !< evaporation number from divergence
  REAL(dp):: zdep_co(kbdim)          !< deposition from convergence
  REAL(dp):: znsub_co(kbdim)         !< sublimation number from divergence
! lookup
  REAL(dp):: zni_lkp(kbdim,klev)     !< change due to lookup over-/unterflow
! adjustment
  REAL(dp):: zcnd_aj(kbdim)          !< condensation mass from adjustment
  REAL(dp):: zdep_aj(kbdim)          !< deposition mass from adjustment
  REAL(dp):: zni_aj(kbdim,klev)      !< deposition number from adjustment
  REAL(dp):: znc_aj(kbdim)           !< condensation number from adjustment
! melting
  REAL(dp):: zqimlt(kbdim)           !< mass melting of ice
  REAL(dp):: zqimlt_rain(kbdim)      !< mass melting of ice that is converted to rain directly
  REAL(dp):: znimlt(kbdim)           !< number melting of ice
  REAL(dp):: znimlt_rain(kbdim)           !< number melting of ice that is converted to rain
! evaporation
  REAL(dp):: zqrevp(kbdim)           !< evaporation of rain
! freezing
  REAL(dp):: zfrl_hom(kbdim)         !< homogeneous freezing rate mass
  REAL(dp):: zfrln_hom(kbdim)        !< homogeneous freezing rate number
  REAL(dp):: zfrl_het(kbdim)         !< heterogeneous freezing rate mass
  REAL(dp):: zfrln_het(kbdim)        !< heterogeneous freezing rate mass
! WBF process
  REAL(dp):: zdep_wbf(kbdim)         !< deposition due to WBF process
  REAL(dp):: zcdnc_wbf(kbdim)        !< cdnc loss due to WBF process
! P3 rates
  REAL(dp):: znislf(kbdim)           !< self-collection of ice
  REAL(dp):: zqccol(kbdim)           !< ice collection of cloud droplets, mass
  REAL(dp):: znccol(kbdim)           !< ice collection of cloud droplets, number
! sublimation
  REAL(dp):: zqisub(kbdim)           !< sublimation of ice mass
  REAL(dp):: znisub(kbdim)           !< sublimation of ice number

! ----------------------------SATURATION ADJUSTMENT-----------------------------
  REAL(dp):: zlc(kbdim,klev)         !< Latent heat of vaporation/sublimation over specific heat depending on temp.
  REAL(dp):: zal(kbdim)              !< Latent heat of vaporation/sublimation depending on temp.
  REAL(dp):: zqcdif(kbdim)           !< difference of timestep change of q and qs
  REAL(dp):: zqs(kbdim,klev)         !< sat. specific humidity depending on t 
  REAL(dp):: zqsp1(kbdim,klev)       !< sat. specific humidity depending on t+0.001K
  REAL(dp):: zqp1(kbdim,klev)        !< temporary updated version of q
  REAL(dp):: ztp1(kbdim,klev)        !< temporary updated version of t

! -----------------------------P3 VALUES----------------------------------------
  REAL(dp):: zrimfrac(kbdim,klev)    !< mass contribution of frozen liquid water to ice (without WBF)
  REAL(dp):: zhetfrac(kbdim,klev)    !< heterogeneous origin ice (includes subsequent growth)
  REAL(dp):: zliqfrac(kbdim,klev)    !< mass contribution of liquid water to ice (with WBF)
  REAL(dp):: zoliqfrac(kbdim,klev)   !< liquid origin ice (includes subsequent growth)
  REAL(dp):: zsosifrac(kbdim,klev)   !< source sink ratio
  REAL(dp):: zsocmfrac(kbdim,klev)   !< source cloud mass ratio
  REAL(dp):: znihetfrac(kbdim,klev)  !< number contribution of nihet
  REAL(dp):: znidetfrac(kbdim,klev)  !< number contribution of nidet
  REAL(dp):: znihomfrac(kbdim,klev)  !< number contribution of nihom
  REAL(dp):: zninucfrac(kbdim,klev)  !< number contribution of ninuc
  REAL(dp):: zrhop(kbdim,klev)       !< riming density
  REAL(dp):: zrhoice(kbdim,klev)     !< ice density
  INTEGER :: p3idx
  REAL(dp), dimension(kbdim,klev) :: zlkp_upper, zlkp_lower, zlkp_col, zlkp_slf, &
                                     zlkp_depx1, zlkp_depx2
  REAL(dp):: rhofacr(kbdim,klev)     !< pressure correction factor ???
  REAL(dp):: rhofaci(kbdim,klev)     !< pressure correction factor ???
  REAL(dp):: rhosui, rhosur          !< ???
!>>DN bugfix: arithmetic exception
!  REAL(dp):: rho_frz                 !< ice density due to freezing of cloud droplets
!  REAL(dp):: rho_rim                 !< ice density due to riming
!<<DN bugfix
  REAL(dp):: p3_mu(kbdim,klev)       !< viscosity of air from P3 code
  REAL(dp):: p3_dv(kbdim,klev)       !< vapor diffusivity of air from P3 code
  REAL(dp):: p3_sc(kbdim,klev)       !< ? of air from P3 code
  REAL(dp):: p3_kap(kbdim,klev)      !< ? of air from P3 code
  REAL(dp):: zice_mu(kbdim,klev)     !< parameter mu of ice PSD
  REAL(dp):: zice_lam(kbdim,klev)    !< parameter lambda of ice PSD

! ----------------------------2 CATEGORY ICE------------------------------------
  REAL(dp):: zqisflx(kbdim)        !< ice crystal mass flux
  REAL(dp):: znisflx(kbdim)        !< ice crystal number flux
  REAL(dp):: zqifal(kbdim,klev)      !< ice mass tendency due to sedimentation
  REAL(dp):: znifal(kbdim,klev)      !< ice number tendency due to sedimentation
  REAL(dp):: zvtimfal(kbdim)      !< fall-speed of sedimenting ice (mass-weighted)
  REAL(dp):: zvtinfal(kbdim)      !< fall-speed of sedimenting ice (number-weighted)
  REAL(dp):: zqsflx(kbdim)           !< snow flux
  REAL(dp):: zqsnow(kbdim)           !< snow mass
  REAL(dp):: zqssub(kbdim)           !< sublimation of snow
  REAL(dp):: zqissub(kbdim)           !< sublimation of falling ice mass
  REAL(dp):: znissub(kbdim)           !< sublimation of falling ice number
  REAL(dp):: zsacl(kbdim,klev)       !< accretion (riming) of droplets by snow - mass
  REAL(dp):: zsacln(kbdim)      !< accretion (riming) of droplets by snow - number
  REAL(dp):: zsaci(kbdim,klev)       !< accretion of ice crystals by snow - mass
  REAL(dp):: zsacin(kbdim,klev)       !< accretion of ice crystals by snow - number
  REAL(dp):: zsaut(kbdim,klev)       !< autoconversion of ice to snow - mass
  REAL(dp):: zsautn(kbdim,klev)       !< autoconversion of ice to snow - number
  REAL(dp):: zqsmlt(kbdim)       !< melting of snow
  REAL(dp):: zqismlt(kbdim)       !< melting of falling ice mass
  REAL(dp):: znismlt(kbdim)       !< melting of falling ice number

  INTEGER, save :: icnt = 0

! --------------------------------IN-CLOUD SCAVENGING---------------------------
  REAL(dp) :: zfrain   (kbdim,klev)       !< rain flux before evaporation [kg/m2/s]
  REAL(dp) :: zfsnow   (kbdim,klev)       !< snow flux before evaporation [kg/m2/s]
  REAL(dp) :: zfevapr  (kbdim,klev)       !< evaporation of rain [kg/m2/s]
  REAL(dp) :: zfsubls  (kbdim,klev)       !< sublimation of snow [kg/m2/s]
  REAL(dp) :: zmsnowacl(kbdim,klev)       !< accretion rate of snow with cloud droplets [kg/kg]
  REAL(dp) :: zmlwc    (kbdim,klev)       !< cloud liquid content before rain [kg/kg]
  REAL(dp) :: zmiwc    (kbdim,klev)       !< cloud ice    content before rain [kg/kg]
  REAL(dp) :: zmratepr (kbdim,klev)       !< rain formation rate in cloudy part [kg/kg]
  REAL(dp) :: zmrateps (kbdim,klev)       !< ice  formation rate in cloudy part [kg/kg]

! --------------------------------DIAGNOSTICS-----------------------------------
  REAL(dp):: ztte_0(kbdim,klev)      !< Temperature tendency before micro
  REAL(dp):: tmax                    !< Maximal temperature at beginning of routine
  REAL(dp):: tmin                    !< Minimal temperature at beginning of routine
  REAL(dp):: zdnc_0(kbdim,klev)      !< CDNC before micro
  REAL(dp):: zdni_0(kbdim,klev)      !< ICNC before micro
  REAL(dp):: zdnc2_0(kbdim,klev)      !< CDNC before micro
  REAL(dp):: zdni2_0(kbdim,klev)      !< ICNC before micro
  REAL(dp):: znisurf(kbdim)           !< ICNC hitting the floor
  REAL(dp):: zreffct(kbdim)          !< cloud top effective radius [m]
  REAL(dp):: zrleff(kbdim,klev)       !< liquid effective radius [um]

! ------------------------------------------------------------------------------
!>>DN
  REAL(dp):: zacc_stm(kbdim,klev)        !< Rain formation rate for mass (accretion)   [kg/kg]
  REAL(dp):: zaut_stm(kbdim,klev)        !< Rain formation rate for mass (autoconversion)   [kg/kg]
  REAL(dp):: zqccol_stm(kbdim,klev)      !< ice collection of cloud droplets, mass
!<<DN

  zaclc_tm1(1:kproma,:) = cloud_cover_duplic(1:kproma,:,krow)

!--- diagnostics
  ztte_0(1:kproma,:) = ptte(1:kproma,:)
  zdni_0(1:kproma,:) = pxtte(1:kproma,:,idt_icnc)
  zdnc_0(1:kproma,:) = pxtte(1:kproma,:,idt_cdnc)
  zdni2_0(1:kproma,:) = pxtm1(1:kproma,:,idt_icnc)
  zdnc2_0(1:kproma,:) = MERGE(nconstnc, pxtm1(1:kproma,:,idt_cdnc), lconstnc)

  pacdnc(1:kbdim,klev)    = 0._dp
  zicnc(1:kbdim,klev)     = 0._dp
  prelhum(1:kbdim,1:klev) = 0._dp
  pssfl(1:kbdim)          = 0._dp
  prsfl(1:kbdim)          = 0._dp
  znisurf(1:kbdim)        = 0._dp
  zsedtm                  = 0._dp
  zsedtn                  = 0._dp

  zclcpre(:)       = 1._dp
  zclcfi(:)        = 0._dp
  zreffct(:)       = 0._dp

! initialize P3 sedimentation
  zqiflx(:)   = 0._dp
  zniflx(:)   = 0._dp
  zqirimflx(:) = 0._dp
  zbirimflx(:) = 0._dp
  zni_lkp(:,:) = 0._dp

! initialize tendencies
  zdep_ci1(:,:) = 0._dp
  znifrz_ci(:,:) = 0._dp
  zqirimsub_tp(:,:) = 0._dp
  zbirimsub_tp(:,:) = 0._dp
  zqihetsub_tp(:,:) = 0._dp
  zqiliqsub_tp(:,:) = 0._dp
  zqioliqsub_tp(:,:) = 0._dp
  zqsrcsub_tp(:,:) = 0._dp
  zqprcsub_tp(:,:) = 0._dp
  znihetsub_tp(:,:) = 0._dp
  znidetsub_tp(:,:) = 0._dp
  znihomsub_tp(:,:) = 0._dp
  zninucsub_tp(:,:) = 0._dp
  zncnuc(:,:) = 0._dp
  zqte_0 = pqte

! P3 tendencies
  zqirimte(:,:) = 0._dp
  zbirimte(:,:) = 0._dp
  zqihette(:,:) = 0._dp
  zqiliqte(:,:) = 0._dp
  zqioliqte(:,:) = 0._dp
  zqsrcte(:,:) = 0._dp
  zqprcte(:,:) = 0._dp
  znihette(:,:) = 0._dp
  znidette(:,:) = 0._dp
  zninucte(:,:) = 0._dp
  znihomte(:,:) = 0._dp

! cloud cover
  ll_cc(1:kproma,:) = (paclc(1:kproma,:) > clc_min)

  CALL get_precip_fraction(kproma, kbdim, ktdia, klev, paclc, zaclci)

  ll_cci(1:kproma,:) = zaclci(1:kproma,:) > clc_min

  ! Initialize time variables
  zdtime     = delta_time
  zdtime_rcp = 1._dp / zdtime
  ztmst      = time_step_len
  ztmst_rcp  = 1._dp / time_step_len
  ztmstdt    = ztmst * zdtime
  zdt        = ztmst_rcp * zdtime

!SF note: the following constant (zcons1) can't be declared as a param in mo_cloud_utils like many others
!         because vtmpc2 must not be not a parameter (not clear why)
  zcons1     = cpd*vtmpc2
  zcons2     = ztmst_rcp * rgrav
  zcons3     = 1._dp / ( pi*crhosno*cn0s*cvtfall**(1._dp/1.16_dp) )**0.25_dp
!SF note: zcons3 can't be declared as a param in mo_cloud_utils like many others
!         because cvtfall is evaluated at runtime only (tuning param)

!--- Get several utility variables:
  CALL get_util_var( &
          !-- IN
          kproma, kbdim, ktdia, klev, klevp1, &
          paphm1(:,:), pgeo(:,:), papm1(:,:), ptm1(:,:), &
          !-- OUT
          zgeoh(:,:), zdp(:,:), zdpg(:,:), &
          zdz(:,:), zaaa(:,:), zviscos(:,:) )

!>>SF ToDo: move this into get_util_var
  !-- thermodynamic utilities:
  ztmp1_2d(1:kproma,:) = MAX(pqm1(1:kproma,:),0.0_dp)
  ztmp1_2d(1:kproma,:) = 1._dp/(cpd+zcons1*ztmp1_2d(1:kproma,:))
  zlvdcp(1:kproma,:)   = alv*ztmp1_2d(1:kproma,:)
  zlsdcp(1:kproma,:)   = als*ztmp1_2d(1:kproma,:)
!<<SF ToDo

  !--- Retrieve some quantities from boundary condition scheme 
  IF (lstart .OR. lresume) THEN
     IF (lconv) THEN !csld #455
       CALL bc_find('Convective cloud base index', ibc_cvcbot, ierr=ierr)
       CALL bc_find('CAPE contrib. to conv. vertical velocity', ibc_wcape, ierr=ierr)
       CALL bc_find('Temperature in convective scheme', ibc_tconv, ierr=ierr) !SF #368
       CALL bc_find('Detrained condensate', ibc_detr_cond, ierr=ierr) !SF #518
     ENDIF !csld #455
  ENDIF

  !>>csld #455 Security for cases when lconv == .FALSE.
  IF (lconv) THEN
     CALL bc_apply(ibc_cvcbot, kproma, krow, zcvcbot)
     CALL bc_apply(ibc_wcape,  kproma, krow, zwcape)  
     CALL bc_apply(ibc_tconv,  kproma, krow, ztconv) !SF #368
     CALL bc_apply(ibc_detr_cond, kproma, krow, zxtec) !SF #518
  ELSE
     zcvcbot(1:kproma) = 0._dp
     zwcape(1:kproma)  = 0._dp
     !SF note that ztconv does not need to be initialized in this case
     zxtec(1:kproma,:) = 0._dp
  ENDIF
  !<<csld #455

 !------------------------------------------------------------------
!       1.   Top boundary conditions, air density

!       1.2   Air density
  zrho(1:kproma,:)     = papm1(1:kproma,:)/(rd*ptvm1(1:kproma,:))
  zrho_rcp(1:kproma,:) = 1._dp / zrho(1:kproma,:)          
  zrho_corr(1:kproma,:)= 1.3_dp * zrho_rcp(1:kproma,:) ! RDnote: differs from MM code
  rhosur = 100000._dp/(rd*273.15_dp)
  rhosui = 60000._dp/(rd*253.15_dp)
  rhofacr(1:kproma,:) = (rhosur*zrho_rcp(1:kproma,:))**0.54
  rhofaci(1:kproma,:) = (rhosui*zrho_rcp(1:kproma,:))**0.54

  p3_mu(1:kproma,:) = 1.496e-6_dp*ptm1(1:kproma,:)**1.5_dp/(ptm1(1:kproma,:)+120._dp)
  p3_dv(1:kproma,:) = 8.794e-5_dp*ptm1(1:kproma,:)**1.81_dp/papm1(1:kproma,:)
  p3_sc(1:kproma,:) = p3_mu(1:kproma,:)/(zrho(1:kproma,:)*p3_dv(1:kproma,:))
  p3_kap(1:kproma,:)= 1.414e3_dp*p3_mu(1:kproma,:)

!SF Detrained liq water/ice calc.: 

!SF ToDo: convert this loop to array syntax
  DO 122 jk = 1,klev  

     ll1(1:kproma) = (zxtec(1:kproma,jk) > 0.0_dp)

     ztmp1(1:kproma)            = twc_conv(1:kproma,jk,krow) + ztmstdt*zxtec(1:kproma,jk)
     twc_conv(1:kproma,jk,krow) = MERGE(ztmp1(1:kproma), twc_conv(1:kproma,jk,krow), ll1(1:kproma))

     ztmp1(1:kproma)             = conv_time(1:kproma,jk,krow) + zdtime
     conv_time(1:kproma,jk,krow) = MERGE(ztmp1(1:kproma), conv_time(1:kproma,jk,krow), ll1(1:kproma))
 
     pqtec(1:kproma,jk)  = MAX(pqtec(1:kproma,jk),0.0_dp)

122 END DO

!       1.3   Calculate cloud base and cloud top

  CALL get_cloud_bounds( &
          !-- IN
          kproma, kbdim, ktdia, klev, paclc, &
          !- OUT
          itop, ibas, icl_minustop, icl_minusbas, iclnb)

  IF(lcolumn .and. lprescribeaerosols .and. (icnt < 3 .or. lprescribeaerosols_every)) THEN
    CALL prescribe_aerosols_mpace(&
                                  !--IN
                                  kproma, kbdim, klev, ktrac, 10,       & ! int = aerosol case
                                  !--OUT
                                  pxtm1, pxtte)
  END IF
  IF(lcolumn .and. lprescribe_numbers .and. (icnt < 3 .or. lprescribe_numbers_every)) THEN
     pxtm1(1:kproma,:,idt_cdnc) = MERGE(pxtm1(1:kproma,:,idt_cdnc) + nxlmult/600.0*ztmst, &
                                        pxtm1(1:kproma,:,idt_cdnc), &
                                        ll_cc(1:kproma,:))

  ENDIF

!-------------------------------------------------------------------------------------
!                                INITIALIZE LOCAL VARS
!-------------------------------------------------------------------------------------

! grid-mean variables: ice
  zxi(1:kproma,:) = pxim1(1:kproma,:)
  zicnc(1:kproma,:) = pxtm1(1:kproma,:,idt_icnc)
  zqirim(1:kproma,:) = pxtm1(1:kproma,:,idt_qirim)
  zbirim(1:kproma,:) = pxtm1(1:kproma,:,idt_birim)
  zqihet(1:kproma,:) = pxtm1(1:kproma,:,idt_qihet)
  zqiliq(1:kproma,:) = pxtm1(1:kproma,:,idt_qiliq)
  zqioliq(1:kproma,:) = pxtm1(1:kproma,:,idt_qioliq)
  zqsrc(1:kproma,:) = pxtm1(1:kproma,:,idt_qsrc)
  zqprc(1:kproma,:) = pxtm1(1:kproma,:,idt_qprc)
  znihet(1:kproma,:) = pxtm1(1:kproma,:,idt_nihet)
  znidet(1:kproma,:) = pxtm1(1:kproma,:,idt_nidet)
  znihom(1:kproma,:) = pxtm1(1:kproma,:,idt_nihom)
  zninuc(1:kproma,:) = pxtm1(1:kproma,:,idt_ninuc)

! grid-mean variables: liquid
  zxl(1:kproma,:) = pxlm1(1:kproma,:)
  ! zcdnc(1:kproma,:) = MAX(1._dp, MERGE(nconstnc, pxtm1(1:kproma,:,idt_cdnc), lconstnc))
  zcdnc(1:kproma,:) = MERGE(nconstnc, pxtm1(1:kproma,:,idt_cdnc), lconstnc)

  zq(1:kproma,:)  = pqm1(1:kproma,:)
  zt(1:kproma,:)  = ptm1(1:kproma,:)
  ztv(1:kproma,:) = ptvm1(1:kproma,:)*(1._dp+vtmpc1*zq(1:kproma,:)-(zxl(1:kproma,:)+zxi(1:kproma,:)))

  tmax = MAXVAL(zt(1:kproma,:))
  tmin = MINVAL(zt(1:kproma,:))

! in-cloud variables: ice
  ztmp1_2d(1:kproma,:) = 1._dp/MAX(zaclci(1:kproma,:),clc_min)
  zxib(1:kproma,:) = zxi(1:kproma,:)*ztmp1_2d(1:kproma,:)
  zicncb(1:kproma,:) = zicnc(1:kproma,:)*ztmp1_2d(1:kproma,:)
  zqirimb(1:kproma,:) = zqirim(1:kproma,:)*ztmp1_2d(1:kproma,:)
  zbirimb(1:kproma,:) = zbirim(1:kproma,:)*ztmp1_2d(1:kproma,:)
  zqihetb(1:kproma,:) = zqihet(1:kproma,:)*ztmp1_2d(1:kproma,:)
  zqiliqb(1:kproma,:) = zqiliq(1:kproma,:)*ztmp1_2d(1:kproma,:)
  zqioliqb(1:kproma,:) = zqioliq(1:kproma,:)*ztmp1_2d(1:kproma,:)

! in-cloud variables: liquid
  ztmp1_2d(1:kproma,:) = 1._dp/MAX(paclc(1:kproma,:),clc_min)
  zxlb(1:kproma,:) = zxl(1:kproma,:)*ztmp1_2d(1:kproma,:)
  zcdncb(1:kproma,:) = zcdnc(1:kproma,:)*ztmp1_2d(1:kproma,:)

  CALL update_saturation_values(&
          !--IN
          'update_saturation_values: init', tmin, tmax, &
          kproma, kbdim, ktdia, klev, zt, papm1, zq,    &
          !--OUT
          zesw, zqsw, zsupw, zesi, zqsi, zsupi,      &
          zqswp1, zqsip1, zqsw0)

  tmax = MAXVAL(zt(1:kproma,:))
  tmin = MINVAL(zt(1:kproma,:))

! note that these are the same for in-cloud and grid mean (fractions -> b cancels)
  CALL riming_variables(&
          !--IN
          kbdim, klev, kproma, &
          zxi, zqirim, zbirim, &
          !--OUT
          zrimfrac, zrhop)

  CALL secondary_prognostics(&
          !--IN
          kbdim, klev, kproma, &
          zxi, zqihet, zqiliq, zqioliq, zqsrc, zqprc, &
          zicnc, znihet, znihom, zninuc, znidet, &
          !--OUT
          zhetfrac, zliqfrac, zoliqfrac, zsosifrac, &
          znihetfrac, znihomfrac, zninucfrac, znidetfrac)

  ! use gridbox mean b/c the cloud fraction is accounted for when calculating
  ! the riming and self-collection rates. This keeps old and new lookup tables
  ! consistent.
  CALL _secondary_ice_properties_p3(&
          !--IN
          kbdim, klev, kproma, &
          zrhop, zrimfrac, zxi, zicnc, &
          rhofaci, &
          !--OUT
          zrim, zrieff, zriv, &
          zvtin, zvtim, zrhoice, &
          zlkp_col, zlkp_slf, &
          zlkp_depx1, zlkp_depx2, &
          zlkp_upper, zlkp_lower)

  CALL secondary_ice_properties_2m(&
          !--IN
          kbdim, klev, kproma, &
          zrho, zxi, zicnc, &
          zaaa, &
          !--OUT
          zriv_2m, zrieff_2m, zvtim_2m, zvtin_2m)


!-------------------------------------------------------------------------------------
!                           INITIALIZE UTILITY VARIABLES
!-------------------------------------------------------------------------------------

! output stream
  IF(.not. lcolumn) THEN
     swat(1:kproma,:,krow) = MAX(0._dp, zsupw(1:kproma,:))
     sice(1:kproma,:,krow) = MAX(0._dp, zsupi(1:kproma,:))
  ELSE
     swat(1:kproma,:,krow) = zsupw(1:kproma,:)
     sice(1:kproma,:,krow) = zsupi(1:kproma,:)
  ENDIF

  prelhum(1:kproma,:) = zq(1:kproma,:)/zqsw(1:kproma,:)

  CALL vertical_velocity(&
              !--IN
              kproma, kbdim, klev, ktdia, krow, &
              ptkem1(:,:), paphm1(:,:), papm1(:,:), pgeo(:,:), &
              zt(:,:), zq(:,:), zrho_rcp(:,:), pvervel(:,:), zqsi(:,:), &
              !--OUT
              zvervx(:,:))

!Set some more utility variables:
!SF ToDo: put that in get_util_var
  zastbstw(1:kproma,:) = &
                         !zast for water:
                         alv * (alv/(rv*zt(1:kproma,:)) - 1.0_dp) / (ka*zt(1:kproma,:)) &
                         !zbst for water:
                       + rv*zt(1:kproma,:)/(2.21_dp/papm1(1:kproma,:)*zesw(1:kproma,:))

  zkair(1:kproma,:) = 4.1867e-3_dp * (5.69_dp + 0.017_dp*(zt(1:kproma,:)-tmelt)) ! eq. 13-18a Pruppacher & Klett 

  zastbsti(1:kproma,:) = &
                         !zast for ice:
                         als * (als/(rv*zt(1:kproma,:)) - 1.0_dp) / (zkair(1:kproma,:)*zt(1:kproma,:)) &
                         !zbst for ice:
                       + rv*zt(1:kproma,:)/(2.21_dp/papm1(1:kproma,:)*zesi(1:kproma,:))

  ztmp1_2d(1:kproma,:) = MAX(zq(1:kproma,:),0.0_dp)
  ztmp1_2d(1:kproma,:) = 1._dp/(cpd+zcons1*ztmp1_2d(1:kproma,:))
  zlvdcp(1:kproma,:)   = alv*ztmp1_2d(1:kproma,:)
  zlsdcp(1:kproma,:)   = als*ztmp1_2d(1:kproma,:)

  zdv(1:kproma,:) = 2.21_dp/papm1(1:kproma,:)

  ztmp2_2d(1:kproma,:) = 1._dp/MAX(zq(1:kproma,:),eps) &
                    + zlsdcp(1:kproma,:)*alv/(rv*zt(1:kproma,:)**2) !SF prot. against divide by zero 
  ztmp3_2d(1:kproma,:) = grav*(zlvdcp(1:kproma,:)*rd/rv/zt(1:kproma,:)-1._dp)/(rd*zt(1:kproma,:))
  ztmp4_2d(1:kproma,:) = 1._dp/(crhoi*als**2/(zkair(1:kproma,:)*rv*zt(1:kproma,:)**2) & 
                          + crhoi*rv*zt(1:kproma,:)/(zesi(1:kproma,:)*zdv(1:kproma,:)))

  zeta(1:kproma,:) = ztmp2_2d(1:kproma,:)/ztmp3_2d(1:kproma,:)*ztmp4_2d(1:kproma,:)*4._dp*pi*crhoi*cap/zrho(1:kproma,:)

  !--- Interface to aerosol calculations (avoids HAM dependencies in cloud_micro_interface):
  !    (includes activation)

  CALL cloud_subm_1( &
          !-- IN
          kproma, kbdim, klev, krow, ktdia, &
          ptkem1, zwcape, pvervel, zrho, & 
          zrho_rcp, pxtm1, pxtte, zt, papm1, zq, zesw, &
          !-- OUT
          zcdncact, zrwetki, zrwetai, zrwetci, &
          zfracdusol, zfracduai, zfracduci, zfracbcsol, zfracbcinsol, &
          zascs, zapnx, zaprx, zapsigx, ll_het )

!-------------------------------------------------------------------------------------
!                                   MASS TRANSPORT
!-------------------------------------------------------------------------------------

  IF(lmass_transport) THEN
     ! CLEAR-SKY EVAPORATION
     ! 4 moments ice
     zxisub_tp(1:kproma,:) = MAX(0._dp, pxite(1:kproma,:))*(1-zaclci(1:kproma,:))
     znisub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_icnc))*(1-zaclci(1:kproma,:))
     zqirimsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_qirim))*(1-zaclci(1:kproma,:))
     zbirimsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_birim))*(1-zaclci(1:kproma,:))
     zqihetsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_qihet))*(1-zaclci(1:kproma,:))
     zqiliqsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_qiliq))*(1-zaclci(1:kproma,:))
     zqioliqsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_qioliq))*(1-zaclci(1:kproma,:))
     zqsrcsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_qsrc))*(1-zaclci(1:kproma,:))
     zqprcsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_qsrc))*(1-zaclci(1:kproma,:))

     znihetsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_nihet))*(1-zaclci(1:kproma,:))
     znihomsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_nihom))*(1-zaclci(1:kproma,:))
     zninucsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_ninuc))*(1-zaclci(1:kproma,:))
     znidetsub_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_nidet))*(1-zaclci(1:kproma,:))

     ! 2 moments liquid
     zxlevp_tp(1:kproma,:) = MAX(0._dp, pxlte(1:kproma,:))*(1-paclc(1:kproma,:))
     zncevp_tp(1:kproma,:) = MAX(0._dp, pxtte(1:kproma,:,idt_cdnc))*(1-paclc(1:kproma,:))
  ELSE
     zxisub_tp(1:kproma,:) = 0._dp
     znisub_tp(1:kproma,:) = 0._dp
     zqirimsub_tp(1:kproma,:) = 0._dp
     zbirimsub_tp(1:kproma,:) = 0._dp
     zqihetsub_tp(1:kproma,:) = 0._dp
     zqiliqsub_tp(1:kproma,:) = 0._dp
     zqioliqsub_tp(1:kproma,:) = 0._dp
     zqsrcsub_tp(1:kproma,:) = 0._dp
     zqprcsub_tp(1:kproma,:) = 0._dp

     znihetsub_tp(1:kproma,:) = 0._dp
     znihomsub_tp(1:kproma,:) = 0._dp
     zninucsub_tp(1:kproma,:) = 0._dp
     znidetsub_tp(1:kproma,:) = 0._dp

     zxlevp_tp(1:kproma,:) = 0._dp
     zncevp_tp(1:kproma,:) = 0._dp
  ENDIF !lmass_transport

!-------------------------------------------------------------------------------------
!                                NUMBER AND MASS CONVECTION
!-------------------------------------------------------------------------------------

  !RDnote: unsure wether to use zicncb or zicnc. zrim is the same for in-cloud and grid-mean.
  zvervmax(:,:) = threshold_vert_vel(kbdim, kproma, klev, zesw(:,:), &
                                     zesi(:,:), zicnc(:,:), zrho(:,:), &
                                     zrieff(:,:), zeta(:,:))

  !-- Final criterion used for the W-B-F process
  lo2_2d(:,:) = ice_switch(kbdim, kproma, klev, csecfrl, zt, zxi, zvervx, zvervmax)


! MASS
! ------------------------------

  ll1_2d(1:kproma,:) = zt(1:kproma,:) < tmelt
  ll2_2d(1:kproma,:) = zt(1:kproma,:) < cthomi

! Convection has ice for zt < 0 and liquid above
  zxite_cv(1:kproma,:) = MERGE(zxtec(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))
  zxlte_cv(1:kproma,:) = MERGE(0._dp, zxtec(1:kproma,:), ll1_2d(1:kproma,:))

! evaporate clear-sky part
  zxisub_cv(1:kproma,:) = zxite_cv(1:kproma,:)*(1-paclc(1:kproma,:))
  zxlevp_cv(1:kproma,:) = zxlte_cv(1:kproma,:)*(1-paclc(1:kproma,:))

! evaporate all if paclc is below clc_min
  zxisub_cv(1:kproma,:) = MERGE(zxisub_cv(1:kproma,:), zxite_cv(1:kproma,:), ll_cc(1:kproma,:))
  zxlevp_cv(1:kproma,:) = MERGE(zxlevp_cv(1:kproma,:), zxlte_cv(1:kproma,:), ll_cc(1:kproma,:))

  IF(lconvice) THEN
     ! melt ice that would not have formed according to updraft scheme in stratiform cloud
     zximlt_cv(1:kproma,:) = MERGE(0._dp, zxite_cv(1:kproma,:)-zxisub_cv(1:kproma,:), lo2_2d(1:kproma,:))
  ELSE
     ! melt all ice in mixed-phase regime
     zximlt_cv(1:kproma,:) = MERGE(0._dp, zxite_cv(1:kproma,:)-zxisub_cv(1:kproma,:), ll2_2d(1:kproma,:))
  END IF ! lconvice

! NUMBER
! ---------------------------------------

  DO jk=ktdia,klev
     iclbas(1:kproma,jk) = NINT(zcvcbot(1:kproma)) 
     ll1_2d(1:kproma,jk) = (jk == iclbas(1:kproma,jk))
  ENDDO

  ll2_2d(1:kproma,:) = (iclbas(1:kproma,:) > 0) 

!>>DN 2013-07-05
  ztmp1_2d(1:kproma,:) = MERGE(cdncact_cv(1:kproma,:,krow), 0._dp, ll1_2d(1:kproma,:))
!<<DN 2013-07-05

!> RD: convection allows for liquid cloud formation
  ll3_2d(1:kproma,:) = ll2_2d(1:kproma,:)              .AND. &
                       (zxlte_cv(1:kproma,:) > 0._dp)     .AND. &
                       ( (zt(1:kproma,:) > tmelt)    .OR.  &
                         ( (zt(1:kproma,:) > cthomi) .AND. &
                           (zt(1:kproma,:) < tmelt)  .AND. &
                           (.NOT.lo2_2d(1:kproma,:))         &
                         ) &
                        ) !DN #286: modified condition for consistency with stratiform ice / liquid water split

!>>DN 2013-07-05 #271
  ztmp2_2d(1:kproma,:) = cdncact_cv(1:kproma,:,krow)
!<<DN 2013-07-05

  DO jk=ktdia,klev
     DO jl=1,kproma
        jkk = iclbas(jl,jk)
        jkk = MAX(1, jkk) !SF prevents cases where iclbas=0
        ztmp3_2d(jl,jk) = ztmp1_2d(jl,jkk) 
     ENDDO
  ENDDO

  ztmp1_2d(1:kproma,:) = MIN(ztmp2_2d(1:kproma,:),ztmp3_2d(1:kproma,:))
  ztmp1_2d(1:kproma,:) = MAX(0._dp,ztmp1_2d(1:kproma,:))
  zcdnc_cv(1:kproma,:)  = MERGE(ztmp1_2d(1:kproma,:), 0._dp, ll3_2d(1:kproma,:))

  ll1_2d(1:kproma,:) = ll_cc(1:kproma,:) .AND. (zxl(1:kproma,:) >= epsec)

!> RD: cloud .OR. convection allows for liquid cloud formation
  ll2_2d(1:kproma,:) = ll1_2d(1:kproma,:) .OR. ll3_2d(1:kproma,:)

  ztmp1_2d(1:kproma,:) = MERGE(zxlte_cv(1:kproma,:)*ztmst, 0._dp, ll3_2d(1:kproma,:))    !> detr. water
  ztmp2_2d(1:kproma,:) = zxlb(1:kproma,:) + ztmp1_2d(1:kproma,:)                      !> total water
  ztmp2_2d(1:kproma,:) = MERGE(ztmp2_2d(1:kproma,:), zdummy, ll2_2d(1:kproma,:))

  !SF: weighs both stratiform and detrained contribs by their respective liquid water content:
  ztmp3_2d(1:kproma,:) = ( zcdncb(1:kproma,:)*zxlb(1:kproma,:)  &
                      + zcdnc_cv(1:kproma,:)*ztmp1_2d(1:kproma,:) ) &
                      / ztmp2_2d(1:kproma,:)

! RD: change due to entrainment
  zcdnc_cv(1:kproma,:) = ztmp3_2d(1:kproma,:)*paclc(1:kproma,:) - zcdnc(1:kproma,:)        !RD: compute change due to cv

  zcdnc_cv(1:kproma,:) = MAX(0._dp, zcdnc_cv(1:kproma,:))                !DN this ensures that the above weighting
  zcdnc_cv(1:kproma,:) = zcdnc_cv(1:kproma,:)*ztmst_rcp                  !will never yield an unphysical
                                                                         !decrease #271

  ! no numbers in clear-sky case
  zcdnc_cv(1:kproma,:) = MERGE(zcdnc_cv(1:kproma,:), 0._dp, ll_cc(1:kproma,:))

! ice

  ztmp1_2d(1:kproma,:) = zt(1:kproma,:)-tmelt
  ztmp1_2d(1:kproma,:) = MIN(0._dp,ztmp1_2d(1:kproma,:))

  ztmp2_2d(1:kproma,:) = 0.015_dp*ztmp1_2d(1:kproma,:)
  ztmp2_2d(1:kproma,:) = EXP(ztmp2_2d(1:kproma,:))
  ztmp2_2d(1:kproma,:) = 23.2_dp*ztmp2_2d(1:kproma,:)
  ztmp2_2d(1:kproma,:) = MAX(ztmp2_2d(1:kproma,:),1.0_dp)*1e-6
     
!       Simple param of r/re approximating the Schumann et al. 2011 data

  zrid(:,:) = effective_2_volmean_radius_param_Schuman_2011( &
                         kbdim, kproma, klev, ztmp2_2d(:,:))

!> RD: cloud .OR. convection allows for ice cloud formation
  ll_cv(1:kproma,:) = (zxite_cv(1:kproma,:) > 0._dp)     .AND. &
                     ( (zt(1:kproma,:)  < cthomi)  .OR.  &
                       ( (zt(1:kproma,:)  < tmelt) .AND. &
                          lo2_2d(1:kproma,:) .AND. lconvice &
                       ) &
                     ) !DN #286: modified condition for consistency with stratiform ice / liquid water split

!UL changed to consistenly use plates, 1.5.2012 [based on Lohmann et al, ERL 2008, expression (1)]
!>>SF #176 now this is done consistenly between cl. micro and radiation, with a general formula
  ztmp1_2d(1:kproma,:) = conv_effr2mvr*(0.5e-2_dp)**pow_PK*1000._dp/fact_PK*ztmst &
                       * zrho(1:kproma,:)*zxite_cv(1:kproma,:) &
                       / ( MAX(paclc(1:kproma,:),clc_min) * zrid(1:kproma,:)**pow_PK)
!<<SF
  ztmp1_2d(1:kproma,:) = MERGE(MAX(ztmp1_2d(1:kproma,:), 0._dp), 0._dp, ll1_2d(1:kproma,:))

!>>DN 2013-07-05 #271: removed the former IWC weigthing
  zicnc_cv(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:),  0._dp, ll_cv(1:kproma,:))
  zicnc_cv(1:kproma,:) = MAX(zicnc_cv(1:kproma,:), 0._dp) !SF actually this could probably be rearranged more
  zicnc_cv(1:kproma,:) = zicnc_cv(1:kproma,:)*ztmst_rcp    !smartly with the above... (will be done when
                                                           !cleaning up all min values, see #241)

  ! no numbers in clear-sky case
  zicnc_cv(1:kproma,:) = MERGE(zicnc_cv(1:kproma,:), 0._dp, ll_cc(1:kproma,:))
!<<DN 2013-07-05 #271: removed the former IWC weigthing 

!-------------------------------------------------------------------------------------
!                           CIRRUS SCHEME NUMBER AND MASS
!-------------------------------------------------------------------------------------

  IF(lcirrus) THEN
     ll1_2d(1:kproma,:) = (zsupi(1:kproma,:) > 0._dp) .AND. &
                          (zt(1:kproma,:)      < cthomi)

   !SF modified the cond. statement
     IF ( nic_cirrus == 1 ) THEN ! Use ICNC scheme after Lohmann, JAS, 2002

        ztmp2_2d(1:kproma,:) = 0.75_dp*zrho(1:kproma,:)*zsupi(1:kproma,:)*zqsi(1:kproma,:) &
                             / (pi*rhoice*zriv(1:kproma,:)**3)-zicnc(1:kproma,:)
        ztmp3_2d(1:kproma,:) = zascs(1:kproma,:)-zicnc(1:kproma,:)
        ztmp1_2d(1:kproma,:) = MIN(ztmp2_2d(1:kproma,:),ztmp3_2d(1:kproma,:))
        ztmp1_2d(1:kproma,:) = MAX(ztmp1_2d(1:kproma,:),0._dp)

        znifrz_ci(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))

     ELSE IF ( nic_cirrus == 2 ) THEN !--- Use Bernd's cirrus scheme 

   !
   !--- Kaercher and Lohmann, JGR, 2002b; Lohmann et al., JGR, 2004)
   !--- This implies that supersaturation with respect to ice is allowed, thus the depositional
   !--- growth equation needs to be solved
   !

   !--- the freezing rate is limited by the number of aerosols available in each mode. 
        DO jk = klev, ktdia, -1 !SFnote: is it really necessary to go upward?

           ! note that zapnx is in 1/m3, other than whats stated in the submodel interface
           zapnx(1:kproma,jk,1) = 1.e-6_dp*(zapnx(1:kproma,jk,1) - zrho(1:kproma,jk)*zicnc(1:kproma,jk)) ![1/cm3]
   !>>SF: extended the following statement to all freezing modes for security, 
   !      even though nfrzmod is so far limited to 1 anyway
           zapnx(1:kproma,jk,:) = MAX(zapnx(1:kproma,jk,:),1.e-6_dp)
   !<<SF

           zap(1:kproma,jk) = 0._dp !total number of aerosols available (over all modes)

           DO jfrzmod=1,nfrzmod
              zap(1:kproma,jk) = zap(1:kproma,jk) + zapnx(1:kproma,jk,jfrzmod)
           ENDDO 

   !----------- the freezing parameterization returns the number of newly formed ice crystals
   !----------- (znicex) and their size (zri)

               ztmp1(:) = 0._dp
               ztmp1(1:kproma) = MAX(0._dp, zsupi(1:kproma,jk))
               CALL xfrzmstr( &
                       !-- IN
                       ll_het, nosize, ztmst, &
                       klev, kbdim, kproma, jk, nfrzmod, &
                       ztmp1(:), zvervx(:,jk), &
                       zapnx(:,jk,:), zaprx(:,jk,:), zapsigx(:,jk,:), &
                       zt, tmelt, eps, papm1, cthomi, &
                       !-- OUT
                       zri, znifrz_ci)

            ENDDO !jk

        znifrz_ci(1:kproma,:) = MIN(1.e6_dp*zap(1:kproma,:),znifrz_ci(1:kproma,:))
        znifrz_ci(1:kproma,:) = MAX(znifrz_ci(1:kproma,:), 0._dp)

        znifrz_ci(1:kproma,:) = MERGE(znifrz_ci(1:kproma,:), 0._dp,ll1_2d(1:kproma,:))
        znifrz_ci(1:kproma,:) = znifrz_ci(1:kproma,:)*ztmst_rcp*zrho_rcp(1:kproma,:)
           
        ! Schumann parameterization for T-dependant size
        ztmp1_2d(1:kproma,:) = zt(1:kproma,:)-tmelt
        ztmp1_2d(1:kproma,:) = MIN(0._dp,ztmp1_2d(1:kproma,:))

        ztmp2_2d(1:kproma,:) = 0.015_dp*ztmp1_2d(1:kproma,:)
        ztmp2_2d(1:kproma,:) = EXP(ztmp2_2d(1:kproma,:))
        ztmp2_2d(1:kproma,:) = 23.2_dp*ztmp2_2d(1:kproma,:)
        ztmp2_2d(1:kproma,:) = MAX(ztmp2_2d(1:kproma,:),1.0_dp)
!>>DN bugfix
!        zrid(1:kproma,:) = effective_2_volmean_radius_param_Schuman_2011( &
!                               kbdim, kproma, klev, ztmp2_2d(1:kproma,:))
        zrid(:,:) = effective_2_volmean_radius_param_Schuman_2011( &
                               kbdim, kproma, klev, ztmp2_2d(:,:))
!<<DN bugfix

        ! handle the case where zri is too small
        ll1_2d(1:kproma,:) = (zri(1:kproma,:) >= epsec) 
        zri(1:kproma,:) = MERGE(zri(1:kproma,:), zrid(1:kproma,:), ll1_2d(1:kproma,:))
        zri(1:kproma,:) = MAX(zri(1:kproma,:), 1.e-6_dp)

        ! get mass of newly formed particles of size zri
        zdep_ci1(1:kproma,:) = 0.0121_dp*zri(1:kproma,:)**1.9*znifrz_ci(1:kproma,:)

        ! amount missing to reach ice saturation
        ztmp1_2d(1:kproma,:) = (zq(1:kproma,:)-zqsi(1:kproma,:))/  &
                               (1+als**2*zqsi(1:kproma,:)/(cpd*rv*zt(1:kproma,:)**2))

        ! limit deposition to available super-saturation
        zdep_ci1(1:kproma,:) = MERGE(zdep_ci1(1:kproma,:), 0._dp, zt(1:kproma,:) < cthomi)
        zdep_ci1(1:kproma,:) = MAX(MIN(zdep_ci1(1:kproma,:), ztmp1_2d(1:kproma,:)*ztmst_rcp), 0._dp)

        ! no cirrus when there is no cloud
        zdep_ci1(1:kproma,:) = MERGE(zdep_ci1(1:kproma,:), 0._dp, ll_cc(1:kproma,:))
        znifrz_ci(1:kproma,:) = MERGE(znifrz_ci(1:kproma,:), 0._dp, ll_cc(1:kproma,:))

   ! RDskipped: ifdef HAMMOZ, orocirrus_cc

     ELSE IF (nic_cirrus == 3) THEN
        DO jl=1,kproma
           DO jk=ktdia,klev
               CALL CirrusModelInterface( jl,jk,krow, & ! multi processing stuff
                       ! input ambient conditions
                       ztmst, pgeo(jl,jk), & ! timestep, output timestep, geopotential
                       papm1(jl,jk)*0.01, zt(jl,jk), zvervx(jl,jk), & ! pres, temp, updraft
                       !-100.0*rgrav*pvervel(jl,jk)*zrho_rcp(jl,jk), & ! large scale updraft
                       zq(jl,jk), zqsi(jl,jk), & ! for sat calculation and available vapor for deposition

                       ! input aerosols
                       zapnx(jl,jk,1), & ! homogeneous aerosol concentration in m-3
                       zaprx(jl,jk,1)*1e-2, & ! homogeneous aerosol radius in (cm -> m)
                       ! dust concentrations and radii are taken from the aerosol module in the CirrusModelInterface

                       ! input existing ice
                       zxi(jl,jk), zrho(jl,jk), zicnc(jl,jk), & ! stratiform cloudIce and IC in kg/kg and m-3
                       zicnc_cv(jl,jk), zrid(jl,jk), & ! detrainedConc in m-3, detrainedRad in m

                       pxtm1, pxtte, & ! tracer stream and tendency for aerosol concentrations

                       ! input for deposition on ice crystals calculation
                       paclc(jl,jk), zaaa(jl,jk), & ! airDensity, cloudIce, cloudCover, correctionForICfallVelocity
                       zastbsti(jl,jk), zdv(jl,jk), zviscos(jl,jk), & ! themodynTermForIceNucleation, diffusionCoeff, viscosity

                       ! output
                       znifrz_ci(jl,jk), & ! new ice crystals from nucleation
                       zdep_ci1(jl,jk)  & ! deposited water vapor on new ice crystals
                       )
            END DO  !jl
         END DO !jk
         ! convert kg/kg to kg/kg/s
         zdep_ci1(1:kproma,:) = zdep_ci1(1:kproma,:)*ztmst_rcp

         ! convert 1/m3 to 1/kg/s
         znifrz_ci(1:kproma,:) = znifrz_ci(1:kproma,:)*ztmst_rcp*zrho_rcp(1:kproma,:)

         ! amount missing to reach ice saturation
         ztmp1_2d(1:kproma,:) = (zq(1:kproma,:)-zqsi(1:kproma,:))/  &
                                (1+als**2*zqsi(1:kproma,:)/(cpd*rv*zt(1:kproma,:)**2))

         ! add up contribution from pre-existing and newly formed crystals
         zdep_ci1(1:kproma,:) = MERGE(zdep_ci1(1:kproma,:), 0._dp, zt(1:kproma,:) < cthomi)
         zdep_ci1(1:kproma,:) = MAX(MIN(zdep_ci1(1:kproma,:), ztmp1_2d(1:kproma,:)*ztmst_rcp), 0._dp)

         ! no cirrus when there is no cloud
         zdep_ci1(1:kproma,:) = MERGE(zdep_ci1(1:kproma,:), 0._dp, ll_cc(1:kproma,:))
         znifrz_ci(1:kproma,:) = MERGE(znifrz_ci(1:kproma,:), 0._dp, ll_cc(1:kproma,:))

     ENDIF   !which cirrus scheme (nic_cirrus)
  ELSE
     zdep_ci1(1:kproma,:) = 0._dp
     znifrz_ci(1:kproma,:) = 0._dp
  ENDIF !lcirrus

!-------------------------------------------------------------------------------------
!                                NUCLEATION CDNC/ICNC
!-------------------------------------------------------------------------------------

  IF(lact) THEN
     CALL aerosol_activation(&
             !--IN
             kbdim, kproma, klev, ktdia, &
             ibas(:,:), icl_minusbas(:,:), paclc(:,:), &
             zt(:,:), zxl(:,:), zcdnc(:,:), &
             zcdncact(:,:), &
             !--OUT
             zncnuc(:,:))

     ! to mitigate dt-dependance, we need a 'timescale' tau on which aerosol_activation happens.
     ! since the standard basically uses tau=dt, we use 10' here. (== dt in standard)
     ztmp1_2d(1:kproma,:) = 1._dp/(600._dp)
     zncnuc(1:kproma,:) = zncnuc(1:kproma,:)*ztmp1_2d(1:kproma,:)
  ELSE
     zncnuc(1:kproma,:) = 0._dp
  END IF !lact

!-------------------------------------------------------------------------------------
!                                LARGE-SCALE TENDENCIES
!-------------------------------------------------------------------------------------

  zxite_ls(1:kproma,:)      = pxite(1:kproma,:) - zxisub_tp(1:kproma,:) &
                            + zxite_cv(1:kproma,:) - zximlt_cv(1:kproma,:) - zxisub_cv(1:kproma,:)
  zxlte_ls(1:kproma,:)      = pxlte(1:kproma,:) - zxlevp_tp(1:kproma,:) &
                            + zxlte_cv(1:kproma,:) + zximlt_cv(1:kproma,:) - zxlevp_cv(1:kproma,:)

  zqte_ls(1:kproma,:)       = pqte(1:kproma,:) + zxisub_cv(1:kproma,:) + zxlevp_cv(1:kproma,:)  &
                            + zxisub_tp(1:kproma,:) + zxlevp_tp(1:kproma,:)
  ztte_ls(1:kproma,:)       = ptte(1:kproma,:) &
                            - zlvdcp(1:kproma,:)*(zxlevp_cv(1:kproma,:) + zxlevp_tp(1:kproma,:)) &
                            - zlsdcp(1:kproma,:)*(zxisub_cv(1:kproma,:) + zxisub_tp(1:kproma,:)) &
                            - (zlsdcp(1:kproma,:) - zlvdcp(1:kproma,:))&
                               *(zximlt_cv(1:kproma,:))

  zicncte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_icnc) &
                         + zicnc_cv(1:kproma,:) - znisub_tp(1:kproma,:)
  
  zcdncte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_cdnc) &
                         + zcdnc_cv(1:kproma,:) - zncevp_tp(1:kproma,:)

  zqihette_ls(1:kproma,:) = pxtte(1:kproma,:,idt_qihet)  &                                ! large scale tendency...
                          - zqihetsub_tp(1:kproma,:) &                                    ! modified by cover.
                          + MERGE(zxite_cv(1:kproma,:) - zximlt_cv(1:kproma,:) - zxisub_cv(1:kproma,:), &
                                  0._dp, zt(1:kproma,:) > cthomi)                         ! conv below cirrus

  zqiliqte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_qiliq) &
                          - zqiliqsub_tp(1:kproma,:) &                                    ! modified by cover.
                          + zxite_cv(1:kproma,:) - zximlt_cv(1:kproma,:) &                ! convective ice is ass.
                          - zxisub_cv(1:kproma,:)                                         ! to have formed by liq

  ! convective mass counts as liquid origin
  zqioliqte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_qioliq) &
                           - zqioliqsub_tp(1:kproma,:) &
                           + zxite_cv(1:kproma,:) - zximlt_cv(1:kproma,:) &               ! convective ice is ass.
                           - zxisub_cv(1:kproma,:)                                        ! to have formed by liq


  ! accumulate source terms
  zqsrcte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_qsrc) &
                         - zqsrcsub_tp(1:kproma,:) &
                         + zxite_cv(1:kproma,:) + zxlte_cv(1:kproma,:) &
                         - zxisub_cv(1:kproma,:) - zxlevp_cv(1:kproma,:)

  ! accumulate precip sink terms
  zqprcte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_qprc) - zqprcsub_tp(1:kproma,:)

  ! heterogeneous ice number tendency
  znihette_ls(1:kproma,:) = pxtte(1:kproma,:,idt_nihet) - znihetsub_tp(1:kproma,:)

  ! homogeneous ice number tendency
  znihomte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_nihom) - znihomsub_tp(1:kproma,:)

  ! nucleated ice number tendency
  zninucte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_ninuc) - zninucsub_tp(1:kproma,:)

  ! detrained ice number tendency
  znidette_ls(1:kproma,:) = pxtte(1:kproma,:,idt_nidet) - znidetsub_tp(1:kproma,:) &
                          + zicnc_cv(1:kproma,:)

  
  IF(iprog == 1) THEN
     zqirimte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_qirim) - zqirimsub_tp(1:kproma,:) &   ! modify by cover.
                             + zrimfrac(1:kproma,:) &
                               *(zxite_cv(1:kproma,:) - zximlt_cv(1:kproma,:) - zxisub_cv(1:kproma,:))

     ztmp1_2d(1:kproma,:) = MERGE(1._dp/zrhop(1:kproma,:), 0._dp, zrhop(1:kproma,:) > cqtmin)
     zbirimte_ls(1:kproma,:) = pxtte(1:kproma,:,idt_birim) - zbirimsub_tp(1:kproma,:) &   ! modified by cover.
                             + (zxite_cv(1:kproma,:) - zximlt_cv(1:kproma,:) - zxisub_cv(1:kproma,:)) & ! phase changes
                               *ztmp1_2d(1:kproma,:)*zrimfrac(1:kproma,:)
                          
  ELSE IF (iprog == 2) THEN
     zqirimte_ls(1:kproma,:) = 0._dp
     zbirimte_ls(1:kproma,:) = 0._dp
  ENDIF !iprog

  ! Those need to be added to local vars, therefore they appear in the z-tendencies.
  ! They will be reset in the subtimestep call. This could be arranged more intuitively.
  pxite(1:kproma,:) = 0._dp
  pxlte(1:kproma,:) = 0._dp
  ptte(1:kproma,:) = 0._dp
  pqte(1:kproma,:) = 0._dp
  pxtte(1:kproma,:,idt_icnc) = 0._dp
  pxtte(1:kproma,:,idt_cdnc) = 0._dp
  pxtte(1:kproma,:,idt_qirim) = 0._dp
  pxtte(1:kproma,:,idt_birim) = 0._dp
  pxtte(1:kproma,:,idt_qihet) = 0._dp
  pxtte(1:kproma,:,idt_qiliq) = 0._dp
  pxtte(1:kproma,:,idt_qioliq) = 0._dp
  pxtte(1:kproma,:,idt_qsrc) = 0._dp
  pxtte(1:kproma,:,idt_qprc) = 0._dp

  pxtte(1:kproma,:,idt_nihet) = 0._dp
  pxtte(1:kproma,:,idt_nihom) = 0._dp
  pxtte(1:kproma,:,idt_ninuc) = 0._dp
  pxtte(1:kproma,:,idt_nidet) = 0._dp

!-------------------------------------------------------------------------------------
!                        COMPUTE NUMBER OF MICROPHYSICS TIMESTEPS
!-------------------------------------------------------------------------------------
! See also Dietlicher 2018, GMD, Section 4 (Fig 3)

  ! diagnose liquid droplet size
  ztmp1_2d(1:kproma,:) = zxl(1:kproma,:)/MAX(zcdnc(1:kproma,:), cqtmin)
  ztmp1_2d(1:kproma,:) = MAX(ztmp1_2d(1:kproma,:), 0._dp)
  ! mass-weighted mean cloud droplet radius:
  zrcm(1:kproma,:)     = (ztmp1_2d(1:kproma,:)/(1.333*pi*rhoh2o))**0.333

  zxi_guess(1:kproma,:) = zxi(1:kproma,:) + ztmst*zxite_ls(1:kproma,:)
  zicnc_guess(1:kproma,:) = zicnc(1:kproma,:) + ztmst*zicncte_ls(1:kproma,:)
  zbirim_guess(1:kproma,:) = zbirim(1:kproma,:) + ztmst*zbirimte_ls(1:kproma,:)
  zqirim_guess(1:kproma,:) = zqirim(1:kproma,:) + ztmst*zqirimte_ls(1:kproma,:)

  ! calculate riming variables from initial guess
  CALL riming_variables(&
          !--IN
          kbdim, klev, kproma, &
          zxi_guess, zqirim_guess, zbirim_guess, &
          !--OUT
          zf_guess, zrhop_guess)

  ! RD: could be replaced by _read_fall_velocity?
  CALL read_my_fall_velocity_2d(&
          !--IN
          kbdim, kproma, klev, &
          zxi_guess, zicnc_guess, zf_guess, zrhop_guess, zbirim_guess, &
          zrho_rcp, rhofaci, vt_rdc, &
          !--OUT
          zvn_guess, zvm_guess, ztmp1_2d)

  ! diagnose CFL-stability
  dcfl(1:kproma,:,krow) = dcfl(1:kproma,:,krow) + zdtime*zvm_guess(1:kproma,:)*ztmst/zdz(1:kproma,:)

  
  ! note that we use the regular timestep here. This may be changed.
  zvmpot(1:kproma,:) = zvm_guess(1:kproma,:)
  zvnpot(1:kproma,:) = zvn_guess(1:kproma,:)
  ll1_2d(1:kproma,:) = zxi_guess(1:kproma,:) > 1e-6_dp
  dumm1(1:kproma) = 0._dp
  dumn1(1:kproma) = 0._dp
  dumm2(1:kproma) = 0._dp
  dumn2(1:kproma) = 0._dp
  ztmp1(1:kproma) = ztmst
  DO jk=1,klev
     ! reset zvmpot to 0 if there is a significant cloud
     dumm1(1:kproma) = MERGE(0._dp, dumm1(1:kproma), ll1_2d(1:kproma,jk))
     dumn1(1:kproma) = MERGE(0._dp, dumn1(1:kproma), ll1_2d(1:kproma,jk))

     ! ... and set zvmpot to given value
     dumm1(1:kproma) = MAX(zvm_guess(1:kproma,jk), dumm1(1:kproma))
     dumn1(1:kproma) = MAX(zvn_guess(1:kproma,jk), dumn1(1:kproma))

     ! use that zvmpot if there is a cloud, use zvmpot from above if there is no cloud  
     dumm2(1:kproma) = MERGE(dumm1(1:kproma), dumm2(1:kproma), ll1_2d(1:kproma,jk))
     dumn2(1:kproma) = MERGE(dumn1(1:kproma), dumn2(1:kproma), ll1_2d(1:kproma,jk))

     ! handle sedimentation time
     ! subtract time spent in level
     ztmp1(1:kproma) = MERGE(ztmp1(1:kproma) - zdz(1:kproma,jk)/dumm2(1:kproma), 0._dp, &
                             dumm2(1:kproma) > 0._dp)
     ! reset if there is a cloud
     ztmp1(1:kproma) = MERGE(ztmst, ztmp1(1:kproma), ll1_2d(1:kproma,jk))

     ! only set potential fall velocity if level is reachable
     ! here we assume that the lowest cloud level has the fastest ice
     ! this should usually be the case due to gravitational size selection
     zvmpot(1:kproma,jk) = MERGE(dumm2(1:kproma), 0._dp, ztmp1(1:kproma) > 0)
     zvnpot(1:kproma,jk) = MERGE(dumn2(1:kproma), 0._dp, ztmp1(1:kproma) > 0)
  END DO !jk

  ! calculate time spent per level
  ztmp2_2d(1:kproma,:) = MERGE(zdz(1:kproma,:)/zvmpot(1:kproma,:), 0._dp, &
                               zvmpot(1:kproma,:) > cqtmin)
  ! accumulate residence time from bottom up
  ztmp3(1:kproma) = 0._dp
  DO jk=klev,1,-1
     ! add up time spent for levels with ice
     ztmp3(1:kproma) = MERGE(ztmp3(1:kproma) + ztmp2_2d(1:kproma,jk), 0._dp, &
                             zvmpot(1:kproma,jk) > cqtmin)
     ztmp3_2d(1:kproma,jk) = ztmp3(1:kproma)
  END DO ! jk

  ! calculate CFL number up to microprct% of time-step to avoid large nmicro for an
  ! insignificant part of sedimentation
  ztmp1_2d(1:kproma,:) = MERGE(1._dp, 0._dp, &
                               ztmp3_2d(1:kproma,:) > ztmst*microprct .OR. nmicro_max == -1)

  zntmst     = CEILING(MAXVAL(zvmpot(1:kproma,1:klev-1)*ztmst/zdz(1:kproma,1:klev-1)*ztmp1_2d(1:kproma,1:klev-1)))
  zntmst     = MAX(1, zntmst)
  zntmst     = MERGE(nmicro_max, zntmst, nmicro_max > 0 .AND. zntmst > nmicro_max)
  
  zntmst     =  MERGE(nmicro, zntmst, nmicro > 0)
  zstmst     = ztmst/zntmst
  zstmst_rcp = 1._dp/zstmst
  zntmst_rcp = 1._dp/zntmst


!--- End included for CDNC/IC scheme -----------------------------------

!-------------------------------------------------------------------------------------
!                             SEDIMENTATION LOOP
!-------------------------------------------------------------------------------------

!in-cloud scavenging: accumulate over all sub-steps
zfsnow(:,:) = 0._dp
zfrain(:,:) = 0._dp
zfsubls(:,:) = 0._dp
zfevapr(:,:) = 0._dp
zmsnowacl(:,:) = 0._dp
zmratepr(:,:) = 0._dp
zmrateps(:,:) = 0._dp
!>>DN
zacc_stm(:,:) = 0._dp
zaut_stm(:,:) = 0._dp
zqccol_stm(:,:) = 0._dp
!<<DN

! set to zero diagnostics
zxi_avg(:,:) = 0._dp
zxl_avg(:,:) = 0._dp

IF(ltimer) CALL timer_start(timer_micro)

substep_column: DO nstep=1,zntmst

! set to zero ice fluxes
zqiflx(:) = 0._dp
zniflx(:) = 0._dp
zqirimflx(:) = 0._dp
zbirimflx(:) = 0._dp

zqisten(:,:) = 0._dp
zqristen(:,:) = 0._dp
zbgsten(:,:) = 0._dp
znisten(:,:) = 0._dp
zqhsten(:,:) = 0._dp
zqlsten(:,:) = 0._dp
zqlosten(:,:) = 0._dp

znihetsten(:,:) = 0._dp
znihomsten(:,:) = 0._dp
zninucsten(:,:) = 0._dp
znidetsten(:,:) = 0._dp

! set to zero snow flux
zqisflx(:) = 0._dp
znisflx(:) = 0._dp
zvtimfal(:) = 0._dp
zvtinfal(:) = 0._dp
zqsflx(:) = 0._dp
zqsnow(:) = 0._dp
zsacl(:,:) = 0._dp
zsacln(:) = 0._dp
zsaci(:,:) = 0._dp
zsacin(:,:) = 0._dp
zsaut(:,:) = 0._dp
zsautn(:,:) = 0._dp

!set to zero rain flux
zqrflx(:)  = 0._dp
zqrflx_bfevp(:) = 0._dp
zmltflx(:) = 0._dp
zrpr(:,:)  = 0._dp
zrprn(:,:) = 0._dp
zaut(:,:)  = 0._dp
zacc(:,:)  = 0._dp
znissrc(:,:)  = 0._dp

! compute average cloud water content for liquid fraction diagnostics
zxi_avg(1:kproma,:) = zxi(1:kproma,:)
zxl_avg(1:kproma,:) = zxl(1:kproma,:)

column_processes:  DO jk=ktdia,klev

!     ------------------------------------------------------------------
!     Set zero level-tendencies
!

    zqimlt      = 0._dp
    zqimlt_rain = 0._dp
    znimlt      = 0._dp
    znimlt_rain = 0._dp
    znislf      = 0._dp
    zqccol      = 0._dp
    znccol      = 0._dp
    zqisub      = 0._dp
    znisub      = 0._dp
    zqrevp      = 0._dp
    zfrl_hom    = 0._dp
    zfrln_hom   = 0._dp
    zfrl_het    = 0._dp
    zfrln_het   = 0._dp

!     ------------------------------------------------------------------
!     Microphysics
!

!............................................................
! melting ice

    IF(lmelting) THEN
       ll1(1:kproma) = (zxi(1:kproma,jk) .ge. epsec) .and. &
                       (zt(1:kproma,jk) .ge. tmelt)

       ! Use formulation from Straka Eq. (11.6), neglecting collisions.
       ! ventilation term (f_v)
       ztmp1(1:kproma) = 0.65_dp*zlkp_depx1(1:kproma,jk) + 0.44_dp*zlkp_depx2(1:kproma,jk)* &
                          p3_sc(1:kproma,jk)**(1._dp/3._dp)
       ! lookup table values are for normalized size distribution
       ztmp1(1:kproma) = zicnc(1:kproma,jk)*ztmp1(1:kproma)
       ! Temperature term
       ztmp2(1:kproma) = (rhofaci(1:kproma,jk)*zrho(1:kproma,jk)/p3_mu(1:kproma,jk))**0.5* &
                         (zt(1:kproma,jk)-tmelt)*2*pi*p3_kap(1:kproma,jk)
       ! Humidity term
       ztmp2(1:kproma) = zrho(1:kproma,jk)*p3_dv(1:kproma,jk)*alv*(zqsw0(1:kproma,jk)-zq(1:kproma,jk))

       ! Combine terms
       zqimlt(1:kproma) = ztmp1(1:kproma)*(ztmp2(1:kproma)+ztmp3(1:kproma))/(als - alv)
       zqimlt(1:kproma) = MIN(MAX(zqimlt(1:kproma), 0._dp), zxi(1:kproma,jk)/zstmst)
       zqimlt(1:kproma) = MERGE(zqimlt(1:kproma), 0._dp, ll1(1:kproma))

       ! infer number melting
       znimlt(1:kproma) = MERGE(zqimlt(1:kproma)/zxi(1:kproma,jk)*zicnc(1:kproma,jk), 0._dp, ll1(1:kproma))

       ! If melting is occuring where b==0, we assume it was formed by falling snow.
       zqimlt_rain(1:kproma) = MERGE(0._dp, zqimlt(1:kproma), ll_cc(1:kproma,jk))
       znimlt_rain(1:kproma) = MERGE(0._dp, znimlt(1:kproma), ll_cc(1:kproma,jk))
       zqimlt(1:kproma) = MERGE(zqimlt(1:kproma), 0._dp, ll_cc(1:kproma,jk))
       znimlt(1:kproma) = MERGE(znimlt(1:kproma), 0._dp, ll_cc(1:kproma,jk))
    ELSE
       zqimlt(1:kproma) = 0._dp
       znimlt(1:kproma) = 0._dp
       zqimlt_rain(1:kproma) = 0._dp
       znimlt_rain(1:kproma) = 0._dp
    ENDIF !lmelting

    IF(lmelting .AND. (lpiggy .OR. iprog == 2)) THEN
       ll1(1:kproma) = (zqsflx(1:kproma) .ge. epsec) .and. &
                       (zt(1:kproma,jk) .gt. tmelt)
       ll2(1:kproma) = (zqisflx(1:kproma) .ge. epsec) .and. &
                       (zt(1:kproma,jk) .gt. tmelt)
       ll3(1:kproma) = (znisflx(1:kproma) .ge. epsec) .and. &
                       (zt(1:kproma,jk) .gt. tmelt)
       zqsmlt(1:kproma) = MERGE(zqsflx(1:kproma)/zdpg(1:kproma,jk), 0._dp, ll1(1:kproma))
       zqismlt(1:kproma) = MERGE(zqisflx(1:kproma)/zdpg(1:kproma,jk), 0._dp, ll2(1:kproma))
       znismlt(1:kproma) = MERGE(znisflx(1:kproma)/zdpg(1:kproma,jk), 0._dp, ll3(1:kproma))
    ELSE
       zqsmlt(1:kproma) = 0._dp
       zqismlt(1:kproma) = 0._dp
       znismlt(1:kproma) = 0._dp
    ENDIF !lmelting

!............................................................
! self-collection of ice

    IF(lself_collection) THEN
       ll1(1:kproma) = zxi(1:kproma,jk) .ge. epsec .and. &
                       zicnc(1:kproma,jk) .ge. epsec

       ! inverse of cover as zlkp is calculated using grid-mean values
       ztmp1(1:kproma) = MERGE(1._dp/zaclci(1:kproma,jk), 0._dp, zaclci(1:kproma,jk) > clc_min)

       ! gridbox-mean value for selfcollection
       znislf(1:kproma) = MERGE(zlkp_slf(1:kproma,jk)*rhofaci(1:kproma,jk)*eii*zrho(1:kproma,jk), 0._dp, ll1(1:kproma))
       ! as the process rate scales with ni**2, we have to account for the sub-grid-ness
       znislf(1:kproma) = znislf(1:kproma)*ztmp1(1:kproma)**2*ccnislf
    ELSE
       znislf(1:kproma) = 0._dp
    ENDIF !lself_collection

!............................................................
! riming with cloud
    
    IF(lriming .AND. iprog == 1) THEN
       ll1(1:kproma) = zxi(1:kproma,jk) > cqtmin .AND. &
                       zxl(1:kproma,jk) > cqtmin .AND. &
                       ll_cc(1:kproma,jk)

       ! inverse of cover as zlkp is calculated using grid-mean values
       ztmp1(1:kproma) = MERGE(1._dp/zaclci(1:kproma,jk), 0._dp, zaclci(1:kproma,jk) > clc_min)

       ! riming rate from lookup table
       zqccol(1:kproma)  = rhofaci(1:kproma,jk)*zlkp_col(1:kproma,jk)*zxl(1:kproma,jk)*eci*zrho(1:kproma,jk)
       znccol(1:kproma)  = rhofaci(1:kproma,jk)*zlkp_col(1:kproma,jk)*zcdnc(1:kproma,jk)*eci*zrho(1:kproma,jk)

       ! as the above rates scale with ni and nc, we have to account for the sub-grid-ness
       zqccol(1:kproma) = zqccol(1:kproma)*ztmp1(1:kproma)**2*ccqccol
       znccol(1:kproma) = znccol(1:kproma)*ztmp1(1:kproma)**2*ccqccol

       ! limit the rate to gridbox-mean values to avoid negative tendencies
       zqccol(1:kproma) = MIN(zxl(1:kproma,jk)*zstmst_rcp, zqccol(1:kproma))
       znccol(1:kproma) = MIN(zcdnc(1:kproma,jk)*zstmst_rcp, znccol(1:kproma))

       zqccol(1:kproma) = MERGE(zqccol(1:kproma), 0._dp, ll1(1:kproma))
       znccol(1:kproma) = MERGE(znccol(1:kproma), 0._dp, ll1(1:kproma))
    ELSE
       zqccol(1:kproma) = 0._dp
       znccol(1:kproma) = 0._dp
    ENDIF !lriming

    IF(lriming .AND. (lpiggy .OR. iprog == 2)) THEN
       CALL snow_cloud_coll( &
              !--IN
              kbdim, kproma, zstmst, zstmst_rcp, &
              zrho(:,jk), zrho_rcp(:,jk), ll_cc(:,jk), paclc(:,jk), zaclci(:,jk), &
              zt(:,jk), zqsnow, zxib(:,jk), zxlb(:,jk), zicncb(:,jk), zcdncb(:,jk), &
              zviscos, &
              !--OUT
              zsacl(:,jk), zsacln)
    ELSE
       zsacl(1:kproma,jk) = 0._dp
       zsacln(1:kproma) = 0._dp
    ENDIF !lriming

!............................................................
! Warm precip formation

    ! rain mass mixing ratio
    ztmp1(1:kproma) = ( zqrflx(1:kproma)/(12.45_dp*MAX(zclcpre(1:kproma),eps) &
                      * SQRT(zrho_corr(1:kproma,jk))) )**(8._dp/9._dp) 
                      !SF eq. 10.70 of Roeckner et al 2003 (ie MPI report 349),
                      ! 12.45 comes from:
                      ! a_10*(n_0r)**(-1/8) = 90.8*(8.e6)**(-1/8)

    zqr(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll_precip(1:kproma))

    ll1(1:kproma) = ll_cc(1:kproma,jk) &
                    .AND. (zxlb(1:kproma,jk) > cqtmin)

    ztmp1(1:kproma) = MIN(paclc(1:kproma,jk), zclcpre(1:kproma))

    IF(lrain_sed) THEN
       CALL precip_formation_warm_2( &
               !--IN
               kbdim, kproma, zstmst, zstmst_rcp, &
               zaclci(:,jk), zrho(:,jk), zrho_rcp(:,jk),       &
               zcdncb(:,jk), zxlb(:,jk), zqr(:), zrcm(:,jk),     &
               !--OUT
               zrpr(:,jk), zrprn(:,jk), zaut(:,jk), zacc(:,jk))

       zrpr(1:kproma,jk) = MIN(zrpr(1:kproma,jk), zstmst_rcp*zxl(1:kproma,jk))
       zrprn(1:kproma,jk) = MIN(zrprn(1:kproma,jk), zstmst_rcp*zcdnc(1:kproma,jk))
       !>>DN
       zaut(1:kproma,jk) = MIN(zaut(1:kproma,jk), zstmst_rcp*zxl(1:kproma,jk))       
       zacc(1:kproma,jk) = MIN(zacc(1:kproma,jk), zstmst_rcp*zxl(1:kproma,jk))       
       !<<DN

       ! only allow in cloudy area
       zrpr(1:kproma,jk) = MERGE(zrpr(1:kproma,jk), 0._dp, ll_cc(1:kproma,jk))
       zrprn(1:kproma,jk) = MERGE(zrprn(1:kproma,jk), 0._dp, ll_cc(1:kproma,jk))
    ELSE
       zrpr(1:kproma,jk) = 0._dp
       zaut(1:kproma,jk) = 0._dp
       zacc(1:kproma,jk) = 0._dp
       zrprn(1:kproma,jk) = 0._dp
    END IF !rain_sed

!............................................................
! cold precip formation

   IF(lsnow_sed .AND. (lpiggy .OR. iprog == 2)) THEN
      CALL precip_formation_cold( &
              !--IN
              kbdim, kproma, zstmst, zstmst_rcp, &
              zrho(:,jk), zrho_rcp(:,jk), ll_cc(:,jk), paclc(:,jk), zaclci(:,jk), &
              zt(:,jk), zqsnow, zxib(:,jk), zxlb(:,jk), zicncb(:,jk), zcdncb(:,jk), &
              zrim(:,jk), &
              !--OUT
              zsaci(:,jk), zsacin(:,jk), zsaut(:,jk), zsautn(:,jk))
   ELSE 
      zsaci(1:kproma,jk) = 0._dp
      zsacin(1:kproma,jk) = 0._dp
      zsaut(1:kproma,jk) = 0._dp
      zsautn(1:kproma,jk) = 0._dp
   ENDIF !lsnow_sed

!............................................................
! condensation, deposition and WBF

  IF(lsubsat_cnd_dep) THEN
     ztmp1(1:kproma) = MERGE(zqsw(1:kproma,jk)/zqsi(1:kproma,jk)-1._dp, zsupi(1:kproma,jk), &
                             zxl(1:kproma,jk) > xlsmall)
     CALL deposition_l(&
             !--IN
             kproma, kbdim,  &
             zxi(:,jk), zt(:,jk), ztmp1(:), zesi(:,jk), papm1(:,jk),  & 
             zastbsti(:,jk), zlsdcp(:,jk), zvtim(:,jk), &
             zrim(:,jk), zicnc(:,jk), zrimfrac(:,jk), zrho(:,jk), zviscos(:,jk), &
             zlkp_depx1(:,jk), zlkp_depx2(:,jk), rhofaci(:,jk), p3_sc(:,jk), &
             !--OUT
             zdep_a(:))

     zdep_a(1:kproma) = MERGE(zdep_a(1:kproma), 0._dp, zxi(1:kproma,jk) > cqtmin)

     CALL cloud_formation(&
              !--IN
              kbdim, kproma, zstmst, zstmst_rcp, &
              zxl(:,jk), zxi(:,jk), zcdnc(:,jk), zicnc(:,jk), zt(:,jk), zq(:,jk), &
              zlvdcp(:,jk), zlsdcp(:,jk), zqte_0(:,jk), ztte_ls(:,jk), zdep_a(:), &
              zqsi(:,jk), zqsip1(:,jk), zqsw(:,jk), zqswp1(:,jk), zsupi(:,jk), paclc(:,jk), &
              !--OUT
              zdep_co, zcnd_co, znsub_co, znevp_co, &
              zdep_wbf, zcdnc_wbf, zqcdif)
     
     ! turn off condensation/evaporation below cloud threshold (clc_min)
     ! this must be done to be consistent with the clear-sky evaporation.
     ! Otherwise, condensation is allowed (which is the case as soon as paclc>0) but
     ! clear-sky evaporation/sublimation is still active (which is the case for paclc<clc_min)
     zdep_co(1:kproma) = MERGE(zdep_co(1:kproma), 0._dp, ll_cc(1:kproma,jk))
     zcnd_co(1:kproma) = MERGE(zcnd_co(1:kproma), 0._dp, ll_cc(1:kproma,jk))

     znsub_co(1:kproma) = MERGE(znsub_co(1:kproma), 0._dp, ll_cc(1:kproma,jk))
     znevp_co(1:kproma) = MERGE(znevp_co(1:kproma), 0._dp, ll_cc(1:kproma,jk))

     zdep_wbf(1:kproma) = MERGE(zdep_wbf(1:kproma), 0._dp, ll_cc(1:kproma,jk))
     zcdnc_wbf(1:kproma) = MERGE(zcdnc_wbf(1:kproma), 0._dp, ll_cc(1:kproma,jk))

     IF(.NOT. lwbf) THEN
        zdep_wbf(1:kproma) = 0._dp
     END IF !lwbf

     ! only allow aerosol activation to form new droplets if there is condensation
     zncnuc(1:kproma,jk) = MERGE(zncnuc(1:kproma,jk), 0._dp, zcnd_co(1:kproma) > cqtmin)

     ! deposition in cirrus regime (this could be rearranged to use with sublimation)
     CALL deposition_l(&
        !--IN
        kproma, kbdim,  &
        zxi(:,jk), zt(:,jk), zsupi(:,jk), zesi(:,jk), papm1(:,jk),  & 
        zastbsti(:,jk), zlsdcp(:,jk), zvtim(:,jk), &
        zrim(:,jk), zicnc(:,jk), zrimfrac(:,jk), zrho(:,jk), zviscos(:,jk), &
        zlkp_depx1(:,jk), zlkp_depx2(:,jk), rhofaci(:,jk), p3_sc(:,jk), &
        !--OUT
        zdep_ci2(:))

     ! deposition to reach ice saturation
     ztmp1(1:kproma) =(zq(1:kproma,jk)-zqsi(1:kproma,jk))/  &
                      (1+alv**2*zqsi(1:kproma,jk)/(cpd*rv*zt(1:kproma,jk)**2))

     ! only active in cirrus regime where cloud formation is allowed
     zdep_ci2(1:kproma) = MERGE(zdep_ci2(1:kproma), 0._dp, zt(1:kproma,jk) < cthomi .AND. ll_cc(1:kproma,jk))
     ! sublimation is handled separately (again, this might be rearranged)
     zdep_ci2(1:kproma) = MIN(zdep_ci2(1:kproma)+zdep_ci1(1:kproma,jk), ztmp1(1:kproma)*zstmst_rcp)
     zdep_ci2(1:kproma) = MAX(0._dp, zdep_ci2(1:kproma))

  ELSE
     zcnd_co(1:kproma) = 0._dp
     zdep_co(1:kproma) = 0._dp
     znsub_co(1:kproma) = 0._dp
     znevp_co(1:kproma) = 0._dp
     zdep_wbf(1:kproma) = 0._dp
     zdep_a(1:kproma) = 0._dp
     zcdnc_wbf(1:kproma) = 0._dp
     zqcdif(1:kproma) = 0._dp
     zdep_ci2(1:kproma) = 0._dp
     znc_aj(1:kproma) = 0._dp
  ENDIF !lsubsat_cnd_dep

!............................................................
! cloud droplet number adjustment

  IF(inumberadjustment > 0) THEN
     ! compute growth explicitely in case of adjustment
     ! F_k
     ztmp1(1:kproma) = (zlsdcp(1:kproma,jk)/(rv*zt(1:kproma,jk))-1)*zlsdcp(1:kproma,jk) &
                       /(con0_h*zt(1:kproma,jk))
     ! F_d
     ztmp2(1:kproma) = 2.11*1e-5*(zt(1:kproma,jk)/tmelt)**1.94 & ! thermal conductivity of
                     *(101325/papm1(1:kproma,jk))                ! water vapor in air
     ztmp2(1:kproma) = rv*zt(1:kproma,jk)/(ztmp2(1:kproma)*zesi(1:kproma,jk))

     ! condensation (Lohmann, Lueoend, Mahrt 2016 Eq. 7.22)
     ! assuming in-cloud supersaturation of 0.1% (as specified by adjsupsat)
     ztmp3(1:kproma) = zcdnc(1:kproma,jk)*zrcm(1:kproma,jk)*4*pi*adjsupsat &
                      /(ztmp1(1:kproma)+ztmp2(1:kproma))

     ! set minimum size for newly formed cloud droplets subtracting possible condensation
     znc_aj(1:kproma) = MAX(0._dp,zcnd_co(1:kproma)-ztmp3(1:kproma))/(1.333*pi*rhoh2o*rcmax**3._dp)
     ! only adjust if nucleation is too weak
     znc_aj(1:kproma) = MAX(0._dp, znc_aj(1:kproma) - zncnuc(1:kproma,jk))

     ! no adjustment if there is no cloud
     znc_aj(1:kproma) = MERGE(znc_aj(1:kproma), 0._dp, ll_cc(1:kproma,jk))
  ELSE
     znc_aj(1:kproma) = 0._dp
  ENDIF !inumberadjustment

!............................................................
! saturation adjustment

! simple equation from Straka applied so far
  IF(ladjustment) THEN
     ! case 1: liquid adjustment
     ll1(1:kproma) = zt(1:kproma,jk) > cthomi

     ! case 2: adjustment (to liquid) in cirrus (this is a sanity constraint)
     ll2(1:kproma) = zt(1:kproma,jk) < cthomi .AND. zsupw(1:kproma,jk) > 0_dp .AND. &
                     (zicnc(1:kproma,jk) > 1._dp .OR. znifrz_ci(1:kproma,jk) > 0._dp .OR. &
                      zsupw(1:kproma,jk) > 1._dp) ! if there are absolutely no aerosols...

     ! amount missing to reach liquid water saturation
     ztmp1(1:kproma) =(zqsw(1:kproma,jk)-zq(1:kproma,jk))/  &
                      (1+alv**2*zqsw(1:kproma,jk)/(cpd*rv*zt(1:kproma,jk)**2))
     
     ! previous change in humidity
     ztmp2(1:kproma) = zstmst*(zqte_ls(1:kproma,jk) + &
                               zdep_co(1:kproma) + zcnd_co(1:kproma) + zdep_ci2(1:kproma))

     ! if previous change < amount missing to saturation: do nothing
     ! if previous change > amount missing to saturation: adjust
     ztmp3(1:kproma) = ztmp2(1:kproma) - ztmp1(1:kproma)
     ztmp3(1:kproma) = MAX(0._dp, ztmp3(1:kproma))

     ! apply a constant condensation/deposition rate but limit to supersaturation
     ztmp3(1:kproma) = MIN(1e-8_dp, ztmp3(1:kproma)*zstmst_rcp)

     ! 1: produce water, 2: produce ice
     zcnd_aj(1:kproma) = MERGE(ztmp3(1:kproma), 0._dp, ll1(1:kproma))
     zdep_aj(1:kproma) = MERGE(ztmp3(1:kproma), 0._dp, ll2(1:kproma))

     ! do not allow condensation if there is no cloud
     zcnd_aj(1:kproma) = MERGE(zcnd_aj(1:kproma), 0._dp, ll_cc(1:kproma,jk))
     zdep_aj(1:kproma) = MERGE(zdep_aj(1:kproma), 0._dp, ll_cc(1:kproma,jk))
  ELSE
     zcnd_aj(1:kproma) = 0._dp
     zdep_aj(1:kproma) = 0._dp
  ENDIF !ladjustment

!............................................................
! sublimation

  IF(iprog == 1) THEN
     IF(isublimation == 1) THEN
        CALL deposition_l(&
                !--IN
                kproma, kbdim,  &
                zxi(:,jk), zt(:,jk), zsupi(:,jk), zesi(:,jk), papm1(:,jk),  & 
                zastbsti(:,jk), zlsdcp(:,jk), zvtim(:,jk), &
                zrim(:,jk), zicnc(:,jk), zrimfrac(:,jk), zrho(:,jk), zviscos(:,jk), &
                zlkp_depx1(:,jk), zlkp_depx2(:,jk), rhofaci(:,jk), p3_sc(:,jk), &
                !--OUT
                zqisub(:))

        ! sublimation to reach ice saturation
        ztmp1(1:kproma) =(zqsi(1:kproma,jk)-zq(1:kproma,jk))/  &
                         (1+alv**2*zqsi(1:kproma,jk)/(cpd*rv*zt(1:kproma,jk)**2))
        ! we are only interested in subsaturation regions where zqsi > zq
        ztmp1(1:kproma) = MAX(0._dp, ztmp1(1:kproma))

        ! zqisub is negative, we are only interested in sublimation part.
        zqisub(1:kproma) = MAX(0._dp, -zqisub(1:kproma))
        ! limit sublimation to ice saturation
        zqisub(1:kproma) = MIN(zqisub(1:kproma), ztmp1(1:kproma)*zstmst_rcp)

        ll1(1:kproma)    = zxi(1:kproma,jk) > cqtmin
        znisub(1:kproma) = MERGE(zicnc(1:kproma,jk)/zxi(1:kproma,jk)*zqisub(1:kproma), 0._dp, ll1(1:kproma))
     ELSE IF(isublimation == 2) THEN
        ztmp1(1:kproma) = MIN(0._dp, zsupi(1:kproma,jk))
        CALL deposition(&
                !--IN
                kproma, kbdim,  &
                zxi(:,jk), zt(:,jk), ztmp1(:),      & 
                zastbsti(:,jk), paclc(:,jk), zrho(:,jk), zaaa(:,jk), zviscos(:,jk), &
                zdv(:,jk), zvtim(:,jk), zrim(:,jk), zicncb(:,jk),                    &
                !--OUT
                zqisub(:))

        ! '-' because we look at SUBlimation
        zqisub(1:kproma) = MIN(-zqisub(1:kproma),zxi(1:kproma,jk)*zstmst_rcp)
        zqisub(1:kproma) = MAX(0._dp, zqisub(1:kproma))

        ll1(1:kproma)    = zxi(1:kproma,jk) > cqtmin
        znisub(1:kproma) = MERGE(zicnc(1:kproma,jk)/zxi(1:kproma,jk)*zqisub(1:kproma), 0._dp, ll1(1:kproma))
     ENDIF
  ENDIF !iprog == 1
  
!............................................................
! cover adjustment
  
  IF(lcover_adjustment) THEN
     ! liquid cover adjustment
     zxlevp_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zxl(1:kproma,jk), ll_cc(1:kproma,jk))
     zncevp_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zcdnc(1:kproma,jk), ll_cc(1:kproma,jk))

     ! ice cover adjustment
     zxisub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zxi(1:kproma,jk), ll_cci(1:kproma,jk))
     znisub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zicnc(1:kproma,jk), ll_cci(1:kproma,jk))
     zqirimsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zqirim(1:kproma,jk), ll_cci(1:kproma,jk))
     zbirimsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zbirim(1:kproma,jk), ll_cci(1:kproma,jk))
     
     ! diagnostic tracers
     ll1(1:kproma) = ll_cci(1:kproma,jk) .AND. zxi(1:kproma,jk) > cqtmin
     zqihetsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zqihet(1:kproma,jk), ll1(1:kproma))
     ll1(1:kproma) = ll_cci(1:kproma,jk) .AND. zxi(1:kproma,jk) > cqtmin
     zqiliqsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zqiliq(1:kproma,jk), ll1(1:kproma))
     ll1(1:kproma) = ll_cci(1:kproma,jk) .AND. zxi(1:kproma,jk) > cqtmin
     zqioliqsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zqioliq(1:kproma,jk), ll1(1:kproma))

     ll1(1:kproma) = ll_cc(1:kproma,jk) .AND. zxi(1:kproma,jk) + zxl(1:kproma,jk) > cqtmin
     zqsrcsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zqsrc(1:kproma,jk), ll1(1:kproma))
     ll1(1:kproma) = ll_cc(1:kproma,jk) .AND. zxi(1:kproma,jk) + zxl(1:kproma,jk) > cqtmin
     zqprcsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zqprc(1:kproma,jk), ll1(1:kproma))

     ll1(1:kproma) = ll_cci(1:kproma,jk) .AND. zicnc(1:kproma,jk) > cqtmin
     znihetsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*znihet(1:kproma,jk), ll1(1:kproma))
     ll1(1:kproma) = ll_cci(1:kproma,jk) .AND. zicnc(1:kproma,jk) > cqtmin
     znihomsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*znihom(1:kproma,jk), ll1(1:kproma))
     ll1(1:kproma) = ll_cci(1:kproma,jk) .AND. zicnc(1:kproma,jk) > cqtmin
     zninucsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*zninuc(1:kproma,jk), ll1(1:kproma))
     ll1(1:kproma) = ll_cci(1:kproma,jk) .AND. zicnc(1:kproma,jk) > cqtmin
     znidetsub_cc(1:kproma) = MERGE(0._dp, zstmst_rcp*znidet(1:kproma,jk), ll1(1:kproma))
  ELSE
     zxlevp_cc(1:kproma) = 0._dp
     zncevp_cc(1:kproma) = 0._dp
     
     zxisub_cc(1:kproma) = 0._dp
     znisub_cc(1:kproma) = 0._dp
     zqirimsub_cc(1:kproma) = 0._dp
     zbirimsub_cc(1:kproma) = 0._dp
     zqihetsub_cc(1:kproma) = 0._dp
     zqiliqsub_cc(1:kproma) = 0._dp
     zqioliqsub_cc(1:kproma) = 0._dp

     zqsrcsub_cc(1:kproma) = 0._dp
     zqprcsub_cc(1:kproma) = 0._dp

     znihetsub_cc(1:kproma) = 0._dp
     znihomsub_cc(1:kproma) = 0._dp
     zninucsub_cc(1:kproma) = 0._dp
     znidetsub_cc(1:kproma) = 0._dp
  ENDIF !lcover_adjustment

!............................................................
! evaporation of rain

! precip logicals
    zclcpre(1:kproma) = MERGE(zaclci(1:kproma,jk), 0._dp, zqrflx(1:kproma) > 0._dp)
    zauloc(1:kproma)   = 3./5000._dp*zdz(1:kproma,jk)
    zauloc(1:kproma)   = MAX(MIN(zauloc(1:kproma), clmax), clmin)
    ll1(1:kproma) = (knvb(1:kproma) >= jbmin) .AND. &
                    (knvb(1:kproma) <= jbmax) .AND. &
                    (pvervel(1:kproma,jk) > 0._dp)
    ll2(1:kproma) = (jk == knvb(1:kproma)  ) .OR. &
                    (jk == knvb(1:kproma)+1)
    ll3(1:kproma) = ll1(1:kproma) .AND. ll2(1:kproma) .AND. lonacc
    zauloc(1:kproma) = MERGE(0._dp, zauloc(1:kproma), ll3(1:kproma))

    ll_precip(1:kproma) = (zclcpre(1:kproma) > 0._dp)
    ll_prcp_warm(1:kproma) = ll_cc(1:kproma,jk) &
                       .AND. (zxlb(1:kproma,jk) > cqtmin)! &

    ztmp1(1:kproma) = MERGE(zclcpre(1:kproma), 1._dp, ll_precip(1:kproma))

    IF(levap_rain) THEN
       zqrevp(1:kproma) = -870._dp * MIN(0._dp, zsupw(1:kproma,jk)) &
                         * (zqrflx(1:kproma)/ztmp1(1:kproma))**0.61_dp &
                         / (SQRT(zrho(1:kproma,jk))*zastbstw(1:kproma,jk))

       zqrevp(1:kproma) = MIN(zqrflx(1:kproma)/zdpg(1:kproma,jk), zqrevp(1:kproma))
!>>DN bugfix: limit evaporation of rain to water supersaturation
       zqrevp(1:kproma) = MIN(MAX(zqsw(1:kproma,jk)-pqm1(1:kproma,jk),0._dp), zqrevp(1:kproma))
!<<DN bugfix
    ELSE
       zqrevp(1:kproma) = 0._dp
    ENDIF !levap_rain

!............................................................
! sublimation of snow and falling ice

    IF(isublimation > 0 .AND. (lpiggy .OR. iprog == 2)) THEN
       ! snow
       CALL sublimation_snow(&
               !--IN
               kbdim, kproma, zdpg(:,jk), &
               zq(:,jk), zt(:,jk), zaclci(:,jk), ll_cci(:,jk), zqsflx, &
               zsupi(:,jk), zqsi(:,jk), zrho_rcp(:,jk), zlsdcp(:,jk), &
               !--OUT
               zqssub)
       zqssub(1:kproma) = zqssub(1:kproma)

       ! falling ice
       CALL sublimation_snow(&
               !--IN
               kbdim, kproma, zdpg(:,jk), &
               zq(:,jk), zt(:,jk), zaclci(:,jk), ll_cci(:,jk), zqisflx, &
               zsupi(:,jk), zqsi(:,jk), zrho_rcp(:,jk), zlsdcp(:,jk), &
               !--OUT
               zqissub)
       zqissub(1:kproma) = zqissub(1:kproma)

       ! virtual falling ice mmr
       ll1(1:kproma) = zvtim_2m(1:kproma,jk) > eps
       ztmp1(1:kproma) = zqisflx(1:kproma)*zrho_rcp(1:kproma,jk)/MAX(zvtim_2m(1:kproma,jk), eps)
       ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll1(1:kproma))

       ! virtual falling ice nmr
       ll1(1:kproma) = zvtin_2m(1:kproma,jk) > eps
       ztmp2(1:kproma) = znisflx(1:kproma)*zrho_rcp(1:kproma,jk)/MAX(zvtin_2m(1:kproma,jk), eps)
       ztmp2(1:kproma) = MERGE(ztmp2(1:kproma), 0._dp, ll1(1:kproma))

       ! number fraction
       ll1(1:kproma) = ztmp1(1:kproma) > 0._dp
       ztmp3(1:kproma) = ztmp2(1:kproma)/MAX(ztmp1(1:kproma), eps)
       ztmp3(1:kproma) = MERGE(ztmp2(1:kproma), 0._dp, ll1(1:kproma))

       znissub(1:kproma) = zqissub(1:kproma)*ztmp2(1:kproma)
    ELSE
       zqssub(1:kproma) = 0._dp
       zqissub(1:kproma) = 0._dp
       znissub(1:kproma) = 0._dp
    ENDIF !isublimation

!............................................................
! freezing homogeneous and heterogeneous

! homogeneous freezing logical
    ll1(1:kproma) = (zt(1:kproma,jk) <= cthomi)
! heterogeneous freezing logical
    ll2(1:kproma) = (zxlb(1:kproma,jk) > cqtmin)     &
                     .AND. (zt(1:kproma,jk) < tmelt )  &
                     .AND. (zt(1:kproma,jk) > cthomi ) & 
                     .AND. ll_cc(1:kproma,jk)

    IF(ihomfrz == 1) THEN
       IF (ANY(ll1(1:kproma))) THEN
          CALL freezing_below_238K( &
                  !-- IN
                  kbdim, kproma, zstmst, zstmst_rcp, &
                  ll1(:), paclc(:,jk), &
                  zxl(:,jk), zcdnc(:,jk),   &
                  !-- INOUT
                  zfrl_hom(:), zfrln_hom(:))
         ! no homogeneous freezing outside of clouds
          zfrl_hom(1:kproma) = MERGE(zfrl_hom(1:kproma), 0._dp, ll_cc(1:kproma,jk))
          zfrln_hom(1:kproma) = MERGE(zfrln_hom(1:kproma), 0._dp, ll_cc(1:kproma,jk))
       ENDIF
    ELSE IF(ihomfrz == 2) THEN
       IF (ANY(ll1(1:kproma))) THEN
          CALL homogeneous_freezing( &
                  !-- IN
                  kbdim, kproma, zstmst, zstmst_rcp, &
                  zt(:,jk), paclc(:,jk), zxl(:,jk), zcdnc(:,jk), zrcm(:,jk),     &
!>>DN bugfix
                  ll1(:),&
!<<DN bugfix
                  !-- OUT
                  zfrl_hom, zfrln_hom)
         ! no homogeneous freezing outside of clouds
          zfrl_hom(1:kproma) = MERGE(zfrl_hom(1:kproma), 0._dp, ll_cc(1:kproma,jk))
          zfrln_hom(1:kproma) = MERGE(zfrln_hom(1:kproma), 0._dp, ll_cc(1:kproma,jk))
       ENDIF
    ELSE 
       zfrl_hom(1:kproma) = 0._dp
       zfrln_hom(1:kproma) = 0._dp
    ENDIF !ihomfrz

    IF(lctrl_het_mxphase_frz) THEN
       IF (ANY(ll2(1:kproma))) THEN
          CALL het_mxphase_freezing( &
              !-- IN
              kbdim, kproma, zstmst, zstmst_rcp, &
              ll2(:), papp1(:,jk), ptkem1(:,jk), &
              pvervel(:,jk), paclc(:,jk), zfracbcsol(:,jk), &
              zfracbcinsol(:,jk), zfracdusol(:,jk), zfracduai(:,jk), &
              zfracduci(:,jk),  &
              zrho(:,jk), zrho_rcp(:,jk), zrwetki(:,jk), &
              zrwetai(:,jk), zrwetci(:,jk), zt(:,jk), &
              zicnc(:,jk), zcdnc(:,jk), &
              zxi(:,jk), zxl(:,jk),   &
              !-- OUT
              zfrl_het(:), zfrln_het(:))
         ! no heterogeneous freezing below -35 deg C
         zfrl_het(1:kproma) = MERGE(0._dp, zfrl_het(1:kproma), zt(1:kproma,jk) < cthomi)
         zfrln_het(1:kproma) = MERGE(0._dp, zfrln_het(1:kproma), zt(1:kproma,jk) < cthomi)

         ! no heterogeneous freezing outside of clouds
         zfrl_het(1:kproma) = MERGE(zfrl_het(1:kproma), 0._dp, ll_cc(1:kproma,jk))
         zfrln_het(1:kproma) = MERGE(zfrln_het(1:kproma), 0._dp, ll_cc(1:kproma,jk))
       ENDIF
    ELSE
       zfrl_het(1:kproma) = 0._dp
       zfrln_het(1:kproma) = 0._dp
    ENDIF !lctrl_het_mxphase_frz

!............................................................
! Korolev updraft
    zvervmax(:,jk) = threshold_vert_vel(kbdim, kproma, zesw(:,jk), &
                               zesi(:,jk), zicnc(:,jk), zrho(:,jk), &
                               zrieff(:,jk), zeta(:,jk))
    

!.................................................................
! conservation of water
! vapor
    zsnks(1:kproma) = (MAX(0._dp, -zqte_ls(1:kproma,jk)) + &
                       MAX(0._dp, zcnd_co(1:kproma)) + &
                       MAX(0._dp, zdep_co(1:kproma)) + &
                       zdep_ci2(1:kproma) + zdep_aj(1:kproma) + zcnd_aj(1:kproma))*zstmst

    ! sources are excluded because they are subject to change below (chicken and egg...)
    zsrcs(1:kproma) = zq(1:kproma,jk)

    ztmp3(:) = conservation_reduction(kbdim, kproma, zsrcs, zsnks)

    ! note: *_co are w.r.t. solid and liquid phase. for vapor we need to switch sign
    zcnd_co(1:kproma) = MERGE(zcnd_co(1:kproma), zcnd_co(1:kproma)*ztmp3(1:kproma), &
                              zcnd_co(1:kproma) < 0._dp)
    zdep_co(1:kproma) = MERGE(zdep_co(1:kproma), zdep_co(1:kproma)*ztmp3(1:kproma), &
                              zdep_co(1:kproma) < 0._dp)
    zqte_ls(1:kproma,jk) = MERGE(zqte_ls(1:kproma,jk), zqte_ls(1:kproma,jk)*ztmp3(1:kproma), &
                                 zqte_ls(1:kproma,jk) > 0._dp)
    zdep_ci2(1:kproma) = zdep_ci2(1:kproma)*ztmp3(1:kproma)
    zdep_aj(1:kproma) = zdep_aj(1:kproma)*ztmp3(1:kproma)
    zcnd_aj(1:kproma) = zcnd_aj(1:kproma)*ztmp3(1:kproma)
                       

    ! liquid
    zsnks(1:kproma) = (zfrl_hom(1:kproma) + zfrl_het(1:kproma) + zrpr(1:kproma,jk) &
                       + MAX(0._dp, -zcnd_co(1:kproma)) &
                       + MAX(0._dp, -zxlte_ls(1:kproma,jk)) &
                       + zdep_wbf(1:kproma) + zxlevp_cc(1:kproma))*zstmst
    ! sources are excluded because they are subject to change below (chicken and egg...)
    zsrcs(1:kproma) = zxl(1:kproma,jk)

    IF(iprog == 1) THEN
       zsnks(1:kproma) = zsnks(1:kproma) + zstmst*zqccol(1:kproma)
    ELSE IF(iprog == 2) THEN
       zsnks(1:kproma) = zsnks(1:kproma) + zstmst*zsacl(1:kproma,jk)
    ENDIF !iprog

    ztmp3(:) = conservation_reduction(kbdim, kproma, zsrcs, zsnks)

    zcnd_co(1:kproma) = MERGE(zcnd_co(1:kproma), zcnd_co(1:kproma)*ztmp3(1:kproma), &
                              zcnd_co(1:kproma) > 0._dp)
    zxlte_ls(1:kproma,jk) = MERGE(zxlte_ls(1:kproma,jk), zxlte_ls(1:kproma,jk)*ztmp3(1:kproma), &
                                  zxlte_ls(1:kproma,jk) > 0._dp)
    zdep_wbf(1:kproma) = zdep_wbf(1:kproma)*ztmp3(1:kproma)
    zfrl_hom(1:kproma) = zfrl_hom(1:kproma)*ztmp3(1:kproma)
    zfrl_het(1:kproma) = zfrl_het(1:kproma)*ztmp3(1:kproma)
    zrpr(1:kproma,jk) = zrpr(1:kproma,jk)*ztmp3(1:kproma)
    zxlevp_cc(1:kproma) = zxlevp_cc(1:kproma)*ztmp3(1:kproma)

    IF (iprog == 1) THEN
       zqccol(1:kproma) = zqccol(1:kproma)*ztmp3(1:kproma)
    ELSE IF (iprog ==2) THEN
       zsacl(1:kproma,jk) = zsacl(1:kproma,jk)*ztmp3(1:kproma)
    ENDIF !iprog

    ! liquid number. This needs to be done bc Nc can either vanish or become Ni
    zsnks(1:kproma) = (zfrln_het(1:kproma) + zfrln_hom(1:kproma) + zrprn(1:kproma,jk) &
                       + MAX(0._dp, -zcdncte_ls(1:kproma,jk)) &
                       + MAX(0._dp, -zncnuc(1:kproma,jk)) &
                       + znevp_co(1:kproma) + zcdnc_wbf(1:kproma) &
                       + zncevp_cc(1:kproma))*zstmst
    ! sources are excluded because they are subject to change below (chicken and egg...)
    zsrcs(1:kproma) = zcdnc(1:kproma,jk)

    IF(iprog == 1) THEN
       zsnks(1:kproma) = zsnks(1:kproma) + zstmst*znccol(1:kproma)
    ELSE IF(iprog == 2) THEN
       zsnks(1:kproma) = zsnks(1:kproma) + zstmst*zsacln(1:kproma)
    ENDIF !iprog

    ztmp3(:) = conservation_reduction(kbdim, kproma, zsrcs, zsnks)

    zcdncte_ls(1:kproma,jk) = MERGE(zcdncte_ls(1:kproma,jk), zcdncte_ls(1:kproma,jk)*ztmp3(1:kproma), &
                                    zcdncte_ls(1:kproma,jk) > 0._dp)
    zncnuc(1:kproma,jk) = MERGE(zncnuc(1:kproma,jk), zncnuc(1:kproma,jk)*ztmp3(1:kproma), &
                                zncnuc(1:kproma,jk) > 0._dp)
    znevp_co(1:kproma) = znevp_co(1:kproma)*ztmp3(1:kproma)
    zcdnc_wbf(1:kproma) = zcdnc_wbf(1:kproma)*ztmp3(1:kproma)
    zfrln_het(1:kproma) = zfrln_het(1:kproma)*ztmp3(1:kproma)
    zfrln_hom(1:kproma) = zfrln_hom(1:kproma)*ztmp3(1:kproma)
    zncevp_cc(1:kproma) = zncevp_cc(1:kproma)*ztmp3(1:kproma)
    zrprn(1:kproma,jk) = zrprn(1:kproma,jk)*ztmp3(1:kproma)

    IF (iprog == 1) THEN
       znccol(1:kproma) = znccol(1:kproma)*ztmp3(1:kproma)
    ELSE IF (iprog ==2) THEN
       zsacln(1:kproma) = zsacln(1:kproma)*ztmp3(1:kproma)
    ENDIF !iprog
    

    ! ice
    zsnks(1:kproma) = (zqimlt(1:kproma)+zqimlt_rain(1:kproma) &
                      + MAX(0._dp, -zdep_co(1:kproma)) &
                      + MAX(0._dp, -zxite_ls(1:kproma,jk)) &
                      + zxisub_cc(1:kproma))*zstmst
    zsrcs(1:kproma) = zxi(1:kproma,jk)

    IF(iprog == 1) THEN
       zsnks(1:kproma) = zsnks(1:kproma) + zstmst*zqisub(1:kproma)
    ELSE IF(iprog == 2) THEN
       zsrcs(1:kproma) = zsrcs(1:kproma) + zstmst*zsacl(1:kproma,jk)
    ENDIF !iprog

    ztmp3(:) = conservation_reduction(kbdim, kproma, zsrcs, zsnks)
    
    zdep_co(1:kproma) = MERGE(zdep_co(1:kproma), zdep_co(1:kproma)*ztmp3(1:kproma), &
                              zdep_co(1:kproma) > 0._dp)
    zxite_ls(1:kproma,jk) = MERGE(zxite_ls(1:kproma,jk), zxite_ls(1:kproma,jk)*ztmp3(1:kproma), &
                                  zxite_ls(1:kproma,jk) > 0._dp)
    zqimlt(1:kproma) = zqimlt(1:kproma)*ztmp3(1:kproma)
    zqimlt_rain(1:kproma) = zqimlt_rain(1:kproma)*ztmp3(1:kproma)
    zxisub_cc(1:kproma) = zxisub_cc(1:kproma)*ztmp3(1:kproma)

    IF(iprog == 1) THEN
       zqisub(1:kproma) = zqisub(1:kproma)*ztmp3(1:kproma)
    ENDIF !iprog

    ! ice number
    zsnks(1:kproma) = (znimlt(1:kproma) + znimlt_rain(1:kproma) &
                       + MAX(0._dp, -zicncte_ls(1:kproma,jk)) &
                       + znsub_co(1:kproma) + znisub_cc(1:kproma))*zstmst
    ! sources are excluded because they are subject to change below (chicken and egg...)
    zsrcs(1:kproma) = zicnc(1:kproma,jk)

    IF(iprog == 1) THEN
       zsnks(1:kproma) = zsnks(1:kproma) + zstmst*znisub(1:kproma)
       zsnks(1:kproma) = zsnks(1:kproma) + zstmst*znislf(1:kproma)
    ELSE IF(iprog == 2) THEN
       zsnks(1:kproma) = zsnks(1:kproma) + zstmst*zsacin(1:kproma,jk)
       zsnks(1:kproma) = zsnks(1:kproma) + zstmst*zsautn(1:kproma,jk)
    ENDIF !iprog

    ztmp3(:) = conservation_reduction(kbdim, kproma, zsrcs, zsnks)

    zicncte_ls(1:kproma,jk) = MERGE(zicncte_ls(1:kproma,jk), zicncte_ls(1:kproma,jk)*ztmp3(1:kproma), &
                                    zicncte_ls(1:kproma,jk) > 0._dp)
    znsub_co(1:kproma) = znsub_co(1:kproma)*ztmp3(1:kproma)
    znimlt(1:kproma) = znimlt(1:kproma)*ztmp3(1:kproma)
    znimlt_rain(1:kproma) = znimlt_rain(1:kproma)*ztmp3(1:kproma)
    znisub_cc(1:kproma) = znisub_cc(1:kproma)*ztmp3(1:kproma)

    IF(iprog == 1) THEN
       znisub(1:kproma) = znisub(1:kproma)*ztmp3(1:kproma)
       znislf(1:kproma) = znislf(1:kproma)*ztmp3(1:kproma)
    ELSE IF(iprog == 2) THEN
       zsautn(1:kproma,jk) = zsautn(1:kproma,jk)*ztmp3(1:kproma)
       zsacin(1:kproma,jk) = zsacin(1:kproma,jk)*ztmp3(1:kproma)
    ENDIF !iprog

! snow category
    IF(iprog == 2) THEN
       ! snow
       zsnks(1:kproma) = zqsmlt(1:kproma) + zqssub(1:kproma)
       ! sources are excluded because they are subject to change below (chicken and egg...)
       zsrcs(1:kproma) = zqsflx(1:kproma)/zdpg(1:kproma,jk)
       
       ztmp3(:) = conservation_reduction(kbdim, kproma, zsrcs, zsnks)
    
       zqsmlt(1:kproma) = zqsmlt(1:kproma)*ztmp3(1:kproma)
       zqssub(1:kproma) = zqssub(1:kproma)*ztmp3(1:kproma)

       ! falling ice
       zsnks(1:kproma) = zqismlt(1:kproma) + zqissub(1:kproma)
       zsrcs(1:kproma) = zqisflx(1:kproma)/zdpg(1:kproma,jk)
       
       ztmp3(:) = conservation_reduction(kbdim, kproma, zsrcs, zsnks)
    
       zqismlt(1:kproma) = zqismlt(1:kproma)*ztmp3(1:kproma)
       zqissub(1:kproma) = zqissub(1:kproma)*ztmp3(1:kproma)
       znismlt(1:kproma) = znismlt(1:kproma)*ztmp3(1:kproma)
       znissub(1:kproma) = znissub(1:kproma)*ztmp3(1:kproma)
    ENDIF !iprog

!.................................................................
! in-cloud scavenging: ice microphysics

    ztmp1(1:kproma) = MERGE(1._dp/zaclci(1:kproma,jk), 0._dp, zaclci(1:kproma,jk) > clc_min)
    ! zfsubls is the fraction of sublimating snow. 
    ! Assume that all sublimation and riming is due to 'snow'.
    zfsubls(1:kproma,jk) = zfsubls(1:kproma,jk) &
                         + zdpg(1:kproma,jk)*zqisub(1:kproma)*zntmst_rcp
    zmsnowacl(1:kproma,jk) = zmsnowacl(1:kproma,jk) &
                           + zqccol(1:kproma)*zstmst


!.................................................................
! update

! 4 MOMENT UPDATE
    IF (iprog == 1) THEN 
       zxite(1:kproma,jk) = zxite_ls(1:kproma,jk) - zxisub_cc(1:kproma) &
                            - zqimlt(1:kproma) - zqimlt_rain(1:kproma) &
                            + zfrl_het(1:kproma) + zfrl_hom(1:kproma) &
                            + zqccol(1:kproma) - zqisub(1:kproma) &
                            + zdep_co(1:kproma) + zdep_wbf(1:kproma) &
                            + zdep_ci2(1:kproma) + zdep_aj(1:kproma)

       zxlte(1:kproma,jk) = zxlte_ls(1:kproma,jk) - zxlevp_cc(1:kproma) &
                            - zfrl_hom(1:kproma) - zrpr(1:kproma,jk) &
                            - zfrl_het(1:kproma) + zqimlt(1:kproma) - zqccol(1:kproma) &
                            + zcnd_co(1:kproma) + zcnd_aj(1:kproma) - zdep_wbf(1:kproma)

       zicncte(1:kproma,jk) = zicncte_ls(1:kproma,jk) - znisub_cc(1:kproma) &
                              - znimlt(1:kproma) - znimlt_rain(1:kproma)  &
                              + zfrln_hom(1:kproma) + zfrln_het(1:kproma) &
                              - znislf(1:kproma) - znisub(1:kproma) &
                              - znsub_co(1:kproma) + znifrz_ci(1:kproma,jk)

       zcdncte(1:kproma,jk) = zcdncte_ls(1:kproma,jk) - zncevp_cc(1:kproma) &
                              + znimlt(1:kproma) - znccol(1:kproma) - zrprn(1:kproma,jk)&
                              - zfrln_het(1:kproma) - zfrln_hom(1:kproma) &
                              - znevp_co(1:kproma) - zcdnc_wbf(1:kproma) + zncnuc(1:kproma,jk) &
                              + znc_aj(1:kproma)

       ztte(1:kproma,jk) = ztte_ls(1:kproma,jk) &
                           - (zlsdcp(1:kproma,jk)-zlvdcp(1:kproma,jk)) &
                              *(zqimlt(1:kproma) + zqimlt_rain(1:kproma) - zdep_wbf(1:kproma) &
                                - zfrl_hom(1:kproma) - zfrl_het(1:kproma) - zqccol(1:kproma)) &
                           - zlvdcp(1:kproma,jk)*(zqrevp(1:kproma) - zcnd_co(1:kproma) &
                                                - zcnd_aj(1:kproma) + zxlevp_cc(1:kproma)) &
                           - zlsdcp(1:kproma,jk)*(zqisub(1:kproma) - zdep_co(1:kproma) &
                                                - zdep_aj(1:kproma) + zxisub_cc(1:kproma) &
                                                - zdep_ci2(1:kproma))

       zqte(1:kproma,jk) = zqte_ls(1:kproma,jk) & 
                           + zqrevp(1:kproma) + zqisub(1:kproma) - zcnd_co(1:kproma) - zdep_co(1:kproma) &
                           - zdep_aj(1:kproma) - zcnd_aj(1:kproma) + zxisub_cc(1:kproma) + zxlevp_cc(1:kproma) &
                           - zdep_ci2(1:kproma)

       zqirimte(1:kproma,jk) = zqirimte_ls(1:kproma,jk) - zqirimsub_cc(1:kproma) &
                               + zfrl_het(1:kproma) + zfrl_hom(1:kproma) + zqccol(1:kproma) &
                               + (-zqimlt(1:kproma) - zqimlt_rain(1:kproma) &
                                 - zqisub(1:kproma) + MIN(0._dp, zdep_co(1:kproma))) &
                                 *zrimfrac(1:kproma,jk)

       ztmp1(1:kproma) = MERGE(1._dp/zrhop(1:kproma,jk), 0._dp, zrhop(1:kproma,jk) > cqtmin)

       zbirimte(1:kproma,jk) = zbirimte_ls(1:kproma,jk) - zbirimsub_cc(1:kproma) &
                               + zqccol(1:kproma)/rho_rim &
                               + (zfrl_het(1:kproma) + zfrl_hom(1:kproma))/rho_frz &
                               + (-zqimlt(1:kproma) - zqimlt_rain(1:kproma) &
                                - zqisub(1:kproma) + MIN(0._dp, zdep_co(1:kproma))) &
                                *ztmp1(1:kproma)*zrimfrac(1:kproma,jk)

       zqihette(1:kproma,jk) = zqihette_ls(1:kproma,jk) - zqihetsub_cc(1:kproma) & 
                               + zhetfrac(1:kproma,jk)*(-zqimlt(1:kproma) - zqimlt_rain(1:kproma) &
                               + zqccol(1:kproma) - zqisub(1:kproma) &
                               + zdep_wbf(1:kproma) + zdep_aj(1:kproma) + zdep_co(1:kproma)) &
                               + zfrl_het(1:kproma)

       ! discard formation directly from vapor (co+, aj, ci)
       zqiliqte(1:kproma,jk) = zqiliqte_ls(1:kproma,jk) - zqiliqsub_cc(1:kproma) &
                               + zfrl_het(1:kproma) + zfrl_hom(1:kproma) + zqccol(1:kproma) + zdep_wbf(1:kproma) &
                               + (-zqimlt(1:kproma) - zqimlt_rain(1:kproma) &
                                 - zqisub(1:kproma) + MIN(0._dp, zdep_co(1:kproma))) &
                                 *zliqfrac(1:kproma,jk)

       ! take into account growth by vapor aswell, also in cirrus regime
       ! substract zdep_ci1 as this includes freezing of solution droplets
       zqioliqte(1:kproma,jk) = zqioliqte_ls(1:kproma,jk) - zqioliqsub_cc(1:kproma) &
                              + zfrl_het(1:kproma) + zfrl_hom(1:kproma) &
                              + (zqccol(1:kproma) + zdep_wbf(1:kproma) &
                                 - zqimlt(1:kproma) - zqimlt_rain(1:kproma) &
                                 - zqisub(1:kproma) + zdep_co(1:kproma) + zdep_aj(1:kproma) &
                                 + MAX(0._dp, zdep_ci2(1:kproma) - zdep_ci1(1:kproma,jk))) &
                                 *zoliqfrac(1:kproma,jk)

       ! accumulated source terms
       zqsrcte(1:kproma,jk) = zqsrcte_ls(1:kproma,jk) - zqsrcsub_cc(1:kproma) &
                            + MAX(0._dp, zdep_co(1:kproma)) + zdep_aj(1:kproma)&
                            + MAX(0._dp, zcnd_co(1:kproma)) + zcnd_aj(1:kproma) &
                            + MAX(0._dp, zdep_ci2(1:kproma))

       ! accumulated sink terms
       zqprcte(1:kproma,jk) = zqprcte_ls(1:kproma,jk) - zqprcsub_cc(1:kproma) &
                            + zqimlt_rain(1:kproma) + zrpr(1:kproma,jk)

       ! heterogeneously formed ice crystal number
       znihette(1:kproma,jk) = znihette_ls(1:kproma,jk) - znihetsub_cc(1:kproma) &
                             + zfrln_het(1:kproma) &
                             + znihetfrac(1:kproma,jk)*(-znimlt(1:kproma) - znimlt_rain(1:kproma) &
                             - znisub(1:kproma) - znislf(1:kproma) - znsub_co(1:kproma))

       ! homogeneously formed ice crystal number
       znihomte(1:kproma,jk) = znihomte_ls(1:kproma,jk) - znihomsub_cc(1:kproma) &
                             + zfrln_hom(1:kproma) &
                             + znihomfrac(1:kproma,jk)*(-znimlt(1:kproma) - znimlt_rain(1:kproma) &
                             - znisub(1:kproma) - znislf(1:kproma) - znsub_co(1:kproma))

       ! nucleated ice crystal number
       zninucte(1:kproma,jk) = zninucte_ls(1:kproma,jk) - zninucsub_cc(1:kproma) &
                             + znifrz_ci(1:kproma,jk) &
                             + zninucfrac(1:kproma,jk)*(-znimlt(1:kproma) - znimlt_rain(1:kproma) &
                             - znisub(1:kproma) - znislf(1:kproma) - znsub_co(1:kproma))

       ! detrained ice crystal number
       znidette(1:kproma,jk) = znidette_ls(1:kproma,jk) - znidetsub_cc(1:kproma) &
                             + znidetfrac(1:kproma,jk)*(-znimlt(1:kproma) - znimlt_rain(1:kproma) &
                             - znisub(1:kproma) - znislf(1:kproma) - znsub_co(1:kproma))


       ! RAIN
       ! Diagnostics: compute the fraction of rain
       ztmp1(1:kproma) = MERGE(zmltflx(1:kproma)/zqrflx(1:kproma), 0._dp, zqrflx(1:kproma) > cqtmin)
       ztmp1(1:kproma) = MIN(MAX(0._dp, ztmp1(1:kproma)), 1._dp)

       ! update the rain flux
       zqrflx(1:kproma)  = zqrflx(1:kproma) + zdpg(1:kproma,jk) &
                                *(-zqrevp(1:kproma) + zqimlt_rain(1:kproma) + zrpr(1:kproma,jk))
       zqrflx(1:kproma)  = MAX(0._dp, zqrflx(1:kproma))
       
       ! Diagnostics: store the rain flux that formed through the ice phase
       zmltflx(1:kproma) = zmltflx(1:kproma) + zdpg(1:kproma,jk) &
                           *(-ztmp1(1:kproma)*zqrevp(1:kproma) + zqimlt_rain(1:kproma))
       zmltflx(1:kproma) = MAX(0._dp, MIN(zqrflx(1:kproma), zmltflx(1:kproma)))

! 2 MOMENT UPDATE
    ELSE IF(iprog == 2) THEN
       zxite(1:kproma,jk) = zxite_ls(1:kproma,jk) &
                            - zqimlt(1:kproma) - zqimlt_rain(1:kproma)                    &
                            + zfrl_het(1:kproma) + zfrl_hom(1:kproma)  &
                            - zsaci(1:kproma,jk) - zsaut(1:kproma,jk) &
                            + zdep_co(1:kproma) + zdep_aj(1:kproma) + zdep_wbf(1:kproma)

       zxlte(1:kproma,jk) = zxlte_ls(1:kproma,jk) &
                            - zfrl_hom(1:kproma) - zfrl_het(1:kproma) + zqimlt(1:kproma) &
                            - zsacl(1:kproma,jk) - zrpr(1:kproma,jk) &
                            + zcnd_co(1:kproma) + zcnd_aj(1:kproma) - zdep_wbf(1:kproma)

       zicncte(1:kproma,jk) = zicncte_ls(1:kproma,jk) &
                              - znimlt(1:kproma) - znsub_co(1:kproma) &
                              + zfrln_hom(1:kproma) + zfrln_het(1:kproma) &
                              - zsacin(1:kproma,jk) - zsautn(1:kproma,jk)

       zcdncte(1:kproma,jk) = zcdncte_ls(1:kproma,jk) &
                              + znimlt(1:kproma) - zsacln(1:kproma) &
                              - zfrln_het(1:kproma) - zfrln_hom(1:kproma) - zrprn(1:kproma,jk) &
                              - znevp_co(1:kproma) - zcdnc_wbf(1:kproma)

       ztte(1:kproma,jk) = ztte_ls(1:kproma,jk) &
                           - (zlsdcp(1:kproma,jk)-zlvdcp(1:kproma,jk)) &
                              *(zqimlt(1:kproma) + zqimlt_rain(1:kproma) + zqsmlt(1:kproma)  &
                                + zqismlt(1:kproma)  - zdep_wbf(1:kproma) &
                                - zfrl_hom(1:kproma) - zfrl_het(1:kproma) - zsacl(1:kproma,jk)) &
                           - zlvdcp(1:kproma,jk)*(zqrevp(1:kproma) - zcnd_co(1:kproma) - zcnd_aj(1:kproma)) &
                           - zlsdcp(1:kproma,jk)*(zqssub(1:kproma) + zqissub(1:kproma) &
                                                - zdep_co(1:kproma) - zdep_aj(1:kproma))

       zqte(1:kproma,jk) = zqte_ls(1:kproma,jk) &
                           + zqrevp(1:kproma) + zqssub(1:kproma) + zqissub(1:kproma) &
                           - zcnd_co(1:kproma) - zdep_co(1:kproma) &
                           - zdep_aj(1:kproma) - zcnd_aj(1:kproma)

       zqirimte(1:kproma,jk) = 0._dp
       zbirimte(1:kproma,jk) = 0._dp
       zqihette(1:kproma,jk) = 0._dp
       zqiliqte(1:kproma,jk) = 0._dp
       zqioliqte(1:kproma,jk) = 0._dp
       zqsrcte(1:kproma,jk) = 0._dp
       zqprcte(1:kproma,jk) = 0._dp

       znihette(1:kproma,jk) = 0._dp
       znihomte(1:kproma,jk) = 0._dp
       zninucte(1:kproma,jk) = 0._dp
       znidette(1:kproma,jk) = 0._dp

       ! FALLING ICE
       zqisflx(1:kproma) = zqisflx(1:kproma) + zdpg(1:kproma,jk)* &
                           (-zqismlt(1:kproma) - zqissub(1:kproma))
       zqisflx(1:kproma) = MAX(0._dp, zqisflx(1:kproma))

       znisflx(1:kproma) = znisflx(1:kproma) + zdpg(1:kproma,jk)* &
                           (-znismlt(1:kproma) - znissub(1:kproma))
       znisflx(1:kproma) = MAX(0._dp, znisflx(1:kproma))

       ! SNOW
       zqsflx(1:kproma) = zqsflx(1:kproma) + zdpg(1:kproma,jk)* &
                          (-zqssub(1:kproma) - zqsmlt(1:kproma) &
                           +zsaci(1:kproma,jk) + zsaut(1:kproma,jk) + zsacl(1:kproma,jk))
       zqsflx(1:kproma) = MAX(0._dp, zqsflx(1:kproma))
       zqsnow(1:kproma) = zstmst*zqsflx(1:kproma)/zdpg(1:kproma,jk)

       ! RAIN
       zqrflx(1:kproma)       = zqrflx(1:kproma) + zdpg(1:kproma,jk) &
                                *(-zqrevp(1:kproma) + zqimlt_rain(1:kproma) &
                                  +zqismlt(1:kproma) + zqsmlt(1:kproma) + zrpr(1:kproma,jk))
       zqrflx(1:kproma)   = MAX(0._dp, zqrflx(1:kproma))
    END IF !iprog

    ztmp1(1:kproma) = 1000._dp*(zt(1:kproma,jk)+zstmst*ztte(1:kproma,jk))
    ztmp1(1:kproma) = NINT(ztmp1(1:kproma))
    ll1(1:kproma) = (ztmp1(1:kproma)<jptlucu1 .OR. ztmp1(1:kproma)>jptlucu2)
    IF(ANY(ll1(1:kproma))) THEN
       DO jl=1,kproma
          IF(ll1(jl)) THEN
             WRITE(*,*) 'level', jk, nstep
             WRITE(*,*) 'zxlevp_cc', zlvdcp(jl,jk)*zxlevp_cc(jl)*zstmst
             WRITE(*,*) 'zxlevp_cv', zlvdcp(jl,jk)*zxlevp_cv(jl,jk)*zstmst
             WRITE(*,*) 'zxlevp_tp', zlvdcp(jl,jk)*zxlevp_tp(jl,jk)*zstmst
             WRITE(*,*) 'zcnd_co', zlvdcp(jl,jk)*ABS(zcnd_co(jl))*zstmst
             WRITE(*,*) 'zcnd_aj', zlvdcp(jl,jk)*zcnd_aj(jl)*zstmst
             WRITE(*,*) '-----------------------------------------------'
             WRITE(*,*) 'zxisub_cc', zlsdcp(jl,jk)*zxisub_cc(jl)*zstmst
             WRITE(*,*) 'zxisub_cv', zlsdcp(jl,jk)*zxisub_cv(jl,jk)*zstmst
             WRITE(*,*) 'zxisub_tp', zlsdcp(jl,jk)*zxisub_tp(jl,jk)*zstmst
             WRITE(*,*) 'zdep_co', zlsdcp(jl,jk)*ABS(zdep_co(jl))*zstmst
             WRITE(*,*) 'zdep_aj', zlsdcp(jl,jk)*zdep_aj(jl)*zstmst
             WRITE(*,*) 'zdep_ci1', zlsdcp(jl,jk)*zdep_ci1(jl,jk)*zstmst
             WRITE(*,*) 'zdep_ci2', zlsdcp(jl,jk)*zdep_ci2(jl)*zstmst
             WRITE(*,*) 'zqisub', zlsdcp(jl,jk)*zqisub(jl)*zstmst
             WRITE(*,*) '-----------------------------------------------'
             WRITE(*,*) 'zt', zt(jl,jk)
             WRITE(*,*) 'zq', zq(jl,jk)
             WRITE(*,*) 'zxi', zxi(jl,jk)
             WRITE(*,*) 'zicnc', zicnc(jl,jk)
             WRITE(*,*) 'zrim', zrim(jl,jk)
             WRITE(*,*) 'zrimfrac', zrimfrac(jl,jk)
             WRITE(*,*) 'zrhop', zrhop(jl,jk)
             WRITE(*,*) 'zsupw', zsupw(jl,jk)
             WRITE(*,*) 'zsupi', zsupi(jl,jk)
             WRITE(*,*) 'zqsi/zqsw-1', zqsw(jl,jk)/zqsi(jl,jk)-1
             WRITE(*,*) 'SICE', sice(jl,jk,krow)
             WRITE(*,*) 'allowed', (zq(jl,jk)-zqsi(jl,jk))*zlsdcp(jl,jk)
             WRITE(*,*) 'zqte_ls', zqte_ls(jl,jk)
             WRITE(*,*) 'ztte_ls', ztte_ls(jl,jk)
          END IF
       END DO !jl
    END IF
       
    CALL subtimestep(&
            !--IN
            kbdim, kproma, zstmst, zntmst_rcp, &
            ll_cc(:,jk), ll_cci(:,jk), paclc(:,jk), zaclci(:,jk), &
            zrho(:,jk), zlsdcp(:,jk), zlvdcp(:,jk), &
            !--INOUT
            ! local variables:
            zq(:,jk), zt(:,jk), zxl(:,jk), zcdnc(:,jk), zxi(:,jk), &
            zicnc(:,jk), zqirim(:,jk), zbirim(:,jk), zqihet(:,jk), zqiliq(:,jk), zqioliq(:,jk), &
            zqsrc(:,jk), zqprc(:,jk), &
            znihet(:,jk), znihom(:,jk), zninuc(:,jk), znidet(:,jk), &
            ! local tendencies:
            zqte(:,jk), ztte(:,jk), zxlte(:,jk), zcdncte(:,jk), zxite(:,jk), &
            zicncte(:,jk), zqirimte(:,jk), zbirimte(:,jk), zqihette(:,jk), zqiliqte(:,jk), zqioliqte(:,jk), &
            zqsrcte(:,jk), zqprcte(:,jk), &
            znihette(:,jk), znihomte(:,jk), zninucte(:,jk), znidette(:,jk), &
            ! global tendencies:
            pqte(:,jk), ptte(:,jk), pxlte(:,jk), pxtte(:,jk,idt_cdnc), &
            pxite(:,jk), pxtte(:,jk,idt_icnc), pxtte(:,jk,idt_qirim), pxtte(:,jk,idt_birim),&
            pxtte(:,jk,idt_qihet), pxtte(:,jk,idt_qiliq), pxtte(:,jk,idt_qioliq), &
            pxtte(:,jk,idt_qsrc), pxtte(:,jk,idt_qprc), &
            pxtte(:,jk,idt_nihet), pxtte(:,jk,idt_nihom), pxtte(:,jk,idt_ninuc), pxtte(:,jk,idt_nidet), &
            !--OUT
            zxlb(:,jk), zcdncb(:,jk), &
            zxib(:,jk), zicncb(:,jk), zqirimb(:,jk), zbirimb(:,jk), &
            zqihetb(:,jk), zqiliqb(:,jk), zqioliqb(:,jk), ztv(:,jk))

    CALL update_saturation_values(&
            !--IN
            'update_saturation_values: sedimentation', tmin, tmax, &
            kproma, kbdim, zt(:,jk), papm1(:,jk), zq(:,jk),    &
            !--OUT
            zesw(:,jk), zqsw(:,jk), zsupw(:,jk), zesi(:,jk), &
            zqsi(:,jk), zsupi(:,jk), zqswp1(:,jk), zqsip1(:,jk), zqsw0(:,jk))

    tmax = MAXVAL(zt(1:kproma,:))
    tmin = MINVAL(zt(1:kproma,:))

    CALL riming_variables(&
            !--IN
            kbdim, kproma, &
            zxi(:,jk), zqirim(:,jk), zbirim(:,jk), &
            !--OUT
            zrimfrac(:,jk), zrhop(:,jk))

    CALL secondary_prognostics(&
            !--IN
            kbdim, kproma, &
            zxi(:,jk), zqihet(:,jk), zqiliq(:,jk), zqioliq(:,jk), zqsrc(:,jk), zqprc(:,jk), &
            zicnc(:,jk), znihet(:,jk), znihom(:,jk), zninuc(:,jk), znidet(:,jk), &
            !--OUT
            zhetfrac(:,jk), zliqfrac(:,jk), zoliqfrac(:,jk), zsosifrac(:,jk), &
            znihetfrac(:,jk), znihomfrac(:,jk), zninucfrac(:,jk), znidetfrac(:,jk))

    CALL _read_fall_velocity(&
            !--IN
            kbdim, kproma, &
            zxib(:,jk), zicncb(:,jk), zrimfrac(:,jk), zrhop(:,jk), zbirimb(:,jk), &
            zrho_rcp(:,jk), rhofaci(:,jk), vt_rdc, &
            !--OUT
            zvtin, zvtim, zrim)

    CALL secondary_ice_properties_2m(&
            !--IN
            kbdim, kproma, &
            zrho(:,jk), zxi(:,jk), zicnc(:,jk), &
            zaaa(:,jk), &
            !--OUT
            zriv_2m(:,jk), zrieff_2m(:,jk), zvtim_2m(:,jk), zvtin_2m(:,jk))

    ! diagnose liquid droplet size
    ztmp1(1:kproma) = zxl(1:kproma,jk)/MAX(zcdnc(1:kproma,jk), cqtmin)
    ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0._dp)
    ! mass-weighted mean cloud droplet radius:
    zrcm(1:kproma,jk) = (ztmp1(1:kproma)/(1.333*pi*rhoh2o))**0.333

!............................................................
! cold precip

   IF(lfalling_ice .AND. (lpiggy .OR. iprog == 2)) THEN
      CALL _read_fall_velocity(&
              !--IN
              kbdim, kproma, &
              zxi(:,jk), zicnc(:,jk), zrimfrac(:,jk), zrhop(:,jk), zbirim(:,jk), &
              zrho_rcp(:,jk), rhofaci(:,jk), vt_rdc, &
              !--OUT
              zvtin, zvtim, zrim)

      CALL falling_ice(&
              !--IN
              kbdim, kproma, zstmst, zstmst_rcp, &
              zxi(:,jk), zicnc(:,jk), zvtim(:,jk), zvtim(:,jk), &
              zqisflx, znisflx, &
              zdpg(:,jk), zrho(:,jk), &
              !--INOUT
              zvtimfal, zvtimfal, &
              !--OUT
              zqifal(:,jk), znifal(:,jk))

      ! adjust falling fluxes (note that zqifal is a source for zxi)
      zqisflx(1:kproma) = zqisflx(1:kproma) - zdpg(1:kproma,jk)*zqifal(1:kproma,jk)
      zqisflx(1:kproma) = MAX(0._dp, zqisflx(1:kproma))

      znisflx(1:kproma) = znisflx(1:kproma) - zdpg(1:kproma,jk)*znifal(1:kproma,jk)
      znisflx(1:kproma) = MAX(0._dp, znisflx(1:kproma))
   ELSE
      zqifal(1:kproma,jk) = 0._dp
      znifal(1:kproma,jk) = 0._dp
      zqisflx(1:kproma) = 0._dp
      znisflx(1:kproma) = 0._dp
   ENDIF !lfalling_ice

!.................................................................
! in-cloud scavenging: warm microphysics

    ! liquid
    ztmp1(1:kproma) = MERGE(1._dp/zaclci(1:kproma,jk), 0._dp, ll_cci(1:kproma,jk))
    zfrain(1:kproma,jk) = zfrain(1:kproma,jk) &
                        + (zqrflx(1:kproma) + zqrevp(1:kproma)*zdpg(1:kproma,jk)) &
                          *zntmst_rcp*ztmp1(1:kproma)
    zfevapr(1:kproma,jk) = zfevapr(1:kproma,jk) &
                         + zdpg(1:kproma,jk)*zqrevp(1:kproma)*zntmst_rcp*ztmp1(1:kproma)
    zmratepr(1:kproma,jk) = zmratepr(1:kproma,jk) &
                          + zrpr(1:kproma,jk)*zstmst*ztmp1(1:kproma)


!............................................................
! diagnostics

    dfrl_het(1:kproma,jk,krow) = dfrl_het(1:kproma,jk,krow) + zdtime*zfrl_het(1:kproma)*zntmst_rcp
    dfrl_hom(1:kproma,jk,krow) = dfrl_hom(1:kproma,jk,krow) + zdtime*zfrl_hom(1:kproma)*zntmst_rcp
    dfrln_hom(1:kproma,jk,krow) = dfrln_hom(1:kproma,jk,krow) + zdtime*zfrln_hom(1:kproma)*zntmst_rcp
    dfrln_het(1:kproma,jk,krow) = dfrln_het(1:kproma,jk,krow) + zdtime*zfrln_het(1:kproma)*zntmst_rcp
    dqimlt(1:kproma,jk,krow) = dqimlt(1:kproma,jk,krow) + zdtime*zqimlt(1:kproma)*zntmst_rcp
    dqimlt_rain(1:kproma,jk,krow) = dqimlt_rain(1:kproma,jk,krow) + zdtime*zqimlt_rain(1:kproma)*zntmst_rcp
    dnimlt(1:kproma,jk,krow) = dnimlt(1:kproma,jk,krow) + zdtime*znimlt(1:kproma)*zntmst_rcp
    dnimlt_rain(1:kproma,jk,krow) = dnimlt_rain(1:kproma,jk,krow) + zdtime*znimlt_rain(1:kproma)*zntmst_rcp
    dnislf(1:kproma,jk,krow) = dnislf(1:kproma,jk,krow) + zdtime*znislf(1:kproma)*zntmst_rcp
    dqrflx(1:kproma,jk,krow) = dqrflx(1:kproma,jk,krow) + zdtime*zqrflx(1:kproma)*zntmst_rcp
    dqrevp(1:kproma,jk,krow) = dqrevp(1:kproma,jk,krow) + zdtime*zqrevp(1:kproma)*zntmst_rcp
    dnccol(1:kproma,jk,krow) = dnccol(1:kproma,jk,krow) + zdtime*znccol(1:kproma)*zntmst_rcp
    dqccol(1:kproma,jk,krow) = dqccol(1:kproma,jk,krow) + zdtime*zqccol(1:kproma)*zntmst_rcp
    dqisub(1:kproma,jk,krow) = dqisub(1:kproma,jk,krow) + zdtime*zqisub(1:kproma)*zntmst_rcp
    dnisub(1:kproma,jk,krow) = dnisub(1:kproma,jk,krow) + zdtime*znisub(1:kproma)*zntmst_rcp

    dcnd_co(1:kproma,jk,krow) = dcnd_co(1:kproma,jk,krow) + zdtime*zcnd_co(1:kproma)*zntmst_rcp
    ddep_co(1:kproma,jk,krow) = ddep_co(1:kproma,jk,krow) + zdtime*zdep_co(1:kproma)*zntmst_rcp
    ddep_a(1:kproma,jk,krow) = ddep_a(1:kproma,jk,krow) + zdtime*zdep_a(1:kproma)*zntmst_rcp
    ddep_wbf(1:kproma,jk,krow) = ddep_wbf(1:kproma,jk,krow) + zdtime*zdep_wbf(1:kproma)*zntmst_rcp
    dcdnc_wbf(1:kproma,jk,krow) = dcdnc_wbf(1:kproma,jk,krow) + zdtime*zcdnc_wbf(1:kproma)*zntmst_rcp
    dnevp_co(1:kproma,jk,krow) = dnevp_co(1:kproma,jk,krow) + zdtime*znevp_co(1:kproma)*zntmst_rcp
    dnsub_co(1:kproma,jk,krow) = dnsub_co(1:kproma,jk,krow) + zdtime*znsub_co(1:kproma)*zntmst_rcp
    dcnd_aj(1:kproma,jk,krow) = dcnd_aj(1:kproma,jk,krow) + zdtime*zcnd_aj(1:kproma)*zntmst_rcp
    ddep_aj(1:kproma,jk,krow) = ddep_aj(1:kproma,jk,krow) + zdtime*zdep_aj(1:kproma)*zntmst_rcp
    dqcdif(1:kproma,jk,krow) = dqcdif(1:kproma,jk,krow) + zdtime*zqcdif(1:kproma)*zntmst_rcp
    ddep_ci2(1:kproma,jk,krow) = ddep_ci2(1:kproma,jk,krow) + zdtime*zdep_ci2(1:kproma)*zntmst_rcp
    dnc_aj(1:kproma,jk,krow) = dnc_aj(1:kproma,jk,krow) + zdtime*znc_aj(1:kproma)*zntmst_rcp

    dxisub_cc(1:kproma,jk,krow) = dxisub_cc(1:kproma,jk,krow) + zdtime*zxisub_cc(1:kproma)*zntmst_rcp
    dnisub_cc(1:kproma,jk,krow) = dnisub_cc(1:kproma,jk,krow) + zdtime*znisub_cc(1:kproma)*zntmst_rcp
    dxlevp_cc(1:kproma,jk,krow) = dxlevp_cc(1:kproma,jk,krow) + zdtime*zxlevp_cc(1:kproma)*zntmst_rcp
    dncevp_cc(1:kproma,jk,krow) = dncevp_cc(1:kproma,jk,krow) + zdtime*zncevp_cc(1:kproma)*zntmst_rcp

    dqisflx(1:kproma,jk,krow) = dqisflx(1:kproma,jk,krow) + zdtime*zqisflx(1:kproma)*zntmst_rcp
    dnisflx(1:kproma,jk,krow) = dnisflx(1:kproma,jk,krow) + zdtime*znisflx(1:kproma)*zntmst_rcp
    dqissub(1:kproma,jk,krow) = dqissub(1:kproma,jk,krow) + zdtime*zqissub(1:kproma)*zntmst_rcp
    dnissub(1:kproma,jk,krow) = dnissub(1:kproma,jk,krow) + zdtime*znissub(1:kproma)*zntmst_rcp
    dqismlt(1:kproma,jk,krow) = dqismlt(1:kproma,jk,krow) + zdtime*zqismlt(1:kproma)*zntmst_rcp
    dnismlt(1:kproma,jk,krow) = dnismlt(1:kproma,jk,krow) + zdtime*znismlt(1:kproma)*zntmst_rcp

    dqsnow(1:kproma,jk,krow) = dqsnow(1:kproma,jk,krow) + zdtime*zqsnow(1:kproma)*zntmst_rcp
    dqsflx(1:kproma,jk,krow) = dqsflx(1:kproma,jk,krow) + zdtime*zqsflx(1:kproma)*zntmst_rcp
    dsacl(1:kproma,jk,krow) = dsacl(1:kproma,jk,krow) + zdtime*zsacl(1:kproma,jk)*zntmst_rcp
    dsacln(1:kproma,jk,krow) = dsacln(1:kproma,jk,krow) + zdtime*zsacln(1:kproma)*zntmst_rcp
    dqssub(1:kproma,jk,krow) = dqssub(1:kproma,jk,krow) + zdtime*zqssub(1:kproma)*zntmst_rcp
    dqsmlt(1:kproma,jk,krow) = dqsmlt(1:kproma,jk,krow) + zdtime*zqsmlt(1:kproma)*zntmst_rcp
    dvtimfal(1:kproma,jk,krow) = dvtimfal(1:kproma,jk,krow) + zdtime*zvtimfal(1:kproma)*zntmst_rcp
    dvtinfal(1:kproma,jk,krow) = dvtinfal(1:kproma,jk,krow) + zdtime*zvtinfal(1:kproma)*zntmst_rcp

    dvtim(1:kproma,jk,krow) = dvtim(1:kproma,jk,krow) + zdtime*zvtim(1:kproma,jk)*zntmst_rcp
    dvtin(1:kproma,jk,krow) = dvtin(1:kproma,jk,krow) + zdtime*zvtin(1:kproma,jk)*zntmst_rcp
    drim(1:kproma,jk,krow) = drim(1:kproma,jk,krow) + zdtime*zrim(1:kproma,jk)*zntmst_rcp
    driv(1:kproma,jk,krow) = driv(1:kproma,jk,krow) + zdtime*zriv(1:kproma,jk)*zntmst_rcp
    drieff(1:kproma,jk,krow) = drieff(1:kproma,jk,krow) + zdtime*zrieff(1:kproma,jk)*zntmst_rcp

    driv_2m(1:kproma,jk,krow) = driv_2m(1:kproma,jk,krow) + zdtime*zriv_2m(1:kproma,jk)*zntmst_rcp
    drieff_2m(1:kproma,jk,krow) = drieff_2m(1:kproma,jk,krow) + zdtime*zrieff_2m(1:kproma,jk)*zntmst_rcp
    dvtim_2m(1:kproma,jk,krow) = dvtim_2m(1:kproma,jk,krow) + zdtime*zvtim_2m(1:kproma,jk)*zntmst_rcp
    dvtin_2m(1:kproma,jk,krow) = dvtin_2m(1:kproma,jk,krow) + zdtime*zvtin_2m(1:kproma,jk)*zntmst_rcp

    dqr(1:kproma,jk,krow) = dqr(1:kproma,jk,krow) + zdtime*zqr(1:kproma)*zntmst_rcp

    !>>DN
    zacc_stm(1:kproma,jk)=zacc_stm(1:kproma,jk)+zacc(1:kproma,jk)*zntmst_rcp
    zaut_stm(1:kproma,jk)=zaut_stm(1:kproma,jk)+zaut(1:kproma,jk)*zntmst_rcp
    zqccol_stm(1:kproma,jk)=zqccol_stm(1:kproma,jk)+zqccol(1:kproma)*zntmst_rcp
    !<<DN
    
END DO column_processes

   ! compute the potential fall-speeds
   zvmpot(1:kproma,:) = zvtim(1:kproma,:)
   zvnpot(1:kproma,:) = zvtin(1:kproma,:)
   ll1_2d(1:kproma,:) = zvtim(1:kproma,:) > 0._dp
   ll2_2d(1:kproma,:) = zvtin(1:kproma,:) > 0._dp
   dumm1(1:kproma) = 0._dp
   dumn1(1:kproma) = 0._dp
   dumm2(1:kproma) = 0._dp
   dumn2(1:kproma) = 0._dp
   DO jk=1,klev
      ! reset zvmpot to 0 if there is a cloud
      dumm1(1:kproma) = MERGE(0._dp, dumm1(1:kproma), ll1_2d(1:kproma,jk))
      dumn1(1:kproma) = MERGE(0._dp, dumn1(1:kproma), ll2_2d(1:kproma,jk))

      ! ... and set zvmpot to given value
      dumm1(1:kproma) = MAX(zvtim(1:kproma,jk), dumm1(1:kproma))
      dumn1(1:kproma) = MAX(zvtin(1:kproma,jk), dumn1(1:kproma))

      ! use that zvmpot if there is a cloud, use zvmpot from above if there is no cloud  
      dumm2(1:kproma) = MERGE(dumm1(1:kproma), dumm2(1:kproma), ll1_2d(1:kproma,jk))
      dumn2(1:kproma) = MERGE(dumn1(1:kproma), dumn2(1:kproma), ll2_2d(1:kproma,jk))

      zvmpot(1:kproma,jk) = dumm2(1:kproma)
      zvnpot(1:kproma,jk) = dumn2(1:kproma)
   END DO !jk

  IF(ltimer) CALL timer_start(timer_sed)

  IF(lfalling_ice) THEN
     CALL ice_sedimentation(&
             !--IN
             kbdim, klev, kproma, zstmst, zstmst_rcp,  &
             zdpg, zvtim, zvtin, zdz, zt, &
             zvmpot, zvnpot, zrho, rhofaci, &
             zxi, zicnc, zqirim, zbirim, zqihet, zqiliq, zqioliq, &
             znihet, znihom, zninuc, znidet, &
             !--OUT
             zqisten, zqristen, zbgsten, znisten, zqhsten, zqlsten, zqlosten, &
             znihetsten, znihomsten, zninucsten, znidetsten, &
             zniflx, zqiflx, zqiflx_2d, znsedi)
  ELSE
     zqisten = 0._dp
     zqristen = 0._dp
     zbgsten = 0._dp
     zqhsten = 0._dp
     zqlsten = 0._dp
     zqlosten = 0._dp
     znisten = 0._dp
     znihetsten = 0._dp
     znihomsten = 0._dp
     zninucsten = 0._dp
     znidetsten = 0._dp
     zniflx = 0._dp
     zqiflx = 0._dp
     zqiflx_2d = 0._dp
     znsedi = 0
  ENDIF !lfalling_ice

  IF(ltimer) CALL timer_stop(timer_sed)


  ! average sedimentation of ICs.
  znissrc(1:kproma,:) = znissrc(1:kproma,:) &
                      + zntmst_rcp*znisten(1:kproma,:)

!............................................................
! sedimentation diagnostics

  drpr(1:kproma,:,krow) = drpr(1:kproma,:,krow) + zdtime*zrpr(1:kproma,:)*zntmst_rcp
  daut(1:kproma,:,krow) = daut(1:kproma,:,krow) + zdtime*zaut(1:kproma,:)*zntmst_rcp
  dacc(1:kproma,:,krow) = dacc(1:kproma,:,krow) + zdtime*zacc(1:kproma,:)*zntmst_rcp
  drprn(1:kproma,:,krow) = drprn(1:kproma,:,krow) + zdtime*zrprn(1:kproma,:)*zntmst_rcp

  dsaci(1:kproma,:,krow) = dsaci(1:kproma,:,krow) + zdtime*zsaci(1:kproma,:)*zntmst_rcp
  dsacin(1:kproma,:,krow) = dsacin(1:kproma,:,krow) + zdtime*zsacin(1:kproma,:)*zntmst_rcp
  dsaut(1:kproma,:,krow) = dsaut(1:kproma,:,krow) + zdtime*zsaut(1:kproma,:)*zntmst_rcp
  dsautn(1:kproma,:,krow) = dsautn(1:kproma,:,krow) + zdtime*zsautn(1:kproma,:)*zntmst_rcp

  dqifal(1:kproma,:,krow) = dqifal(1:kproma,:,krow) + zdtime*zqifal(1:kproma,:)*zntmst_rcp
  dnifal(1:kproma,:,krow) = dnifal(1:kproma,:,krow) + zdtime*znifal(1:kproma,:)*zntmst_rcp
  
  dqisten(1:kproma,:,krow) = dqisten(1:kproma,:,krow) + zdtime*zqisten(1:kproma,:)*zntmst_rcp
  dqristen(1:kproma,:,krow) = dqristen(1:kproma,:,krow) + zdtime*zqristen(1:kproma,:)*zntmst_rcp
  dbgsten(1:kproma,:,krow) = dbgsten(1:kproma,:,krow) + zdtime*zbgsten(1:kproma,:)*zntmst_rcp
  dnisten(1:kproma,:,krow) = dnisten(1:kproma,:,krow) + zdtime*znisten(1:kproma,:)*zntmst_rcp
  dqhsten(1:kproma,:,krow) = dqhsten(1:kproma,:,krow) + zdtime*zqhsten(1:kproma,:)*zntmst_rcp

  ztmp1(1:kproma) = 1._dp
  dnsedi(1:kproma,krow) = dnsedi(1:kproma,krow) +zdtime*znsedi*ztmp1(1:kproma)*zntmst_rcp

!.................................................................
! in-cloud scavenging: ice sedimentation

  ! calculate a proxy for snow formation
  ! first limit to 'cloud forming' area where ll_cc == .TRUE.
  ztmp2_2d(1:kproma,:) = MERGE(zqisten(1:kproma,:), 0._dp, ll_cc(1:kproma,:))
  ! then only look at sinks
  ztmp2_2d(1:kproma,:) = -MIN(ztmp2_2d(1:kproma,:), 0._dp)

  ztmp1_2d(1:kproma,:) = MERGE(1._dp/zaclci(1:kproma,:), 0._dp, zaclci(1:kproma,:) > clc_min)
  zfsnow(1:kproma,:) = zfsnow(1:kproma,:) &
                     + zqiflx_2d(1:kproma,:)*zntmst_rcp*ztmp1_2d(1:kproma,:)
  zmrateps(1:kproma,:) = zmrateps(1:kproma,:) &
                       + ztmp2_2d(1:kproma,:)*zstmst*ztmp1_2d(1:kproma,:)


!.................................................................
! surface fluxes
  
  prsfl(1:kproma) = prsfl(1:kproma) + zqrflx(1:kproma)*zntmst_rcp


!.................................................................
! update

  zxlte(1:kproma,:) = 0._dp
  zcdncte(1:kproma,:) = 0._dp

  IF(iprog == 1) THEN
     ! surface fluxes
     pssfl(1:kproma) = pssfl(1:kproma) + zqiflx(1:kproma)*zntmst_rcp
     znisurf(1:kproma) = znisurf(1:kproma) + zniflx(1:kproma)*zntmst_rcp

     zxite(1:kproma,:) = zqisten(1:kproma,:)
     zicncte(1:kproma,:) = znisten(1:kproma,:)
     zqirimte(1:kproma,:) = zqristen(1:kproma,:)
     zbirimte(1:kproma,:) = zbgsten(1:kproma,:)
     zqihette(1:kproma,:) = zqhsten(1:kproma,:)
     zqiliqte(1:kproma,:) = zqlsten(1:kproma,:)
     zqioliqte(1:kproma,:) = zqlosten(1:kproma,:)

     zqsrcte(1:kproma,:) = 0._dp
     zqprcte(1:kproma,:) = MAX(-zqisten(1:kproma,:), 0._dp)

     znihette(1:kproma,:) = znihetsten(1:kproma,:)
     znihomte(1:kproma,:) = znihomsten(1:kproma,:)
     zninucte(1:kproma,:) = zninucsten(1:kproma,:)
     znidette(1:kproma,:) = znidetsten(1:kproma,:)
  ELSE IF(iprog == 2) THEN
     pssfl(1:kproma) = pssfl(1:kproma) + zqsflx(1:kproma)*zntmst_rcp &
                                       + zqisflx(1:kproma)*zntmst_rcp ! contribution from falling ice

     zxite(1:kproma,:) = zqifal(1:kproma,:)
     zicncte(1:kproma,:) = znifal(1:kproma,:)
     zqirimte(1:kproma,:) = 0._dp
     zbirimte(1:kproma,:) = 0._dp
     zqihette(1:kproma,:) = 0._dp
     zqiliqte(1:kproma,:) = 0._dp
     zqioliqte(1:kproma,:) = 0._dp

     zqsrcte(1:kproma,:) = 0._dp
     zqprcte(1:kproma,:) = 0._dp

     znihette(1:kproma,:) = 0._dp
     znihomte(1:kproma,:) = 0._dp
     zninucte(1:kproma,:) = 0._dp
     znidette(1:kproma,:) = 0._dp
  ENDIF !iprog

  ztte = 0._dp
  zqte = 0._dp

  CALL subtimestep(&
          !--IN
          kbdim, klev, kproma, zstmst, zntmst_rcp, &
          ll_cc, ll_cci, paclc, zaclci, &
          zrho, zlsdcp, zlvdcp, &
          !--INOUT
          ! local variables:
          zq, zt, zxl, zcdnc, zxi, &
          zicnc, zqirim, zbirim, zqihet, zqiliq, zqioliq, zqsrc, zqprc, &
          znihet, znihom, zninuc, znidet, &
          ! local tendencies:
          zqte, ztte, zxlte, zcdncte, zxite, &
          zicncte, zqirimte, zbirimte, zqihette, zqiliqte, zqioliqte, zqsrcte, zqprcte, &
          znihette, znihomte, zninucte, znidette, &
          ! global tendencies:
          pqte, ptte, pxlte, pxtte(:,:,idt_cdnc), &
          pxite, pxtte(:,:,idt_icnc), pxtte(:,:,idt_qirim), pxtte(:,:,idt_birim),&
          pxtte(:,:,idt_qihet), pxtte(:,:,idt_qiliq), pxtte(:,:,idt_qioliq), &
          pxtte(:,:,idt_qsrc), pxtte(:,:,idt_qprc), &
          pxtte(:,:,idt_nihet), pxtte(:,:,idt_nihom), pxtte(:,:,idt_ninuc), pxtte(:,:,idt_nidet), &
          !--OUT
          zxlb, zcdncb, &
          zxib, zicncb, zqirimb, zbirimb, zqihetb, zqiliqb, zqioliqb, ztv)

  CALL update_saturation_values(&
          !--IN
          'update_saturation_values: sedimentation', tmin, tmax, &
          kproma, kbdim, ktdia, klev, zt, papm1, zq,    &
          !--OUT
          zesw, zqsw, zsupw, zesi, zqsi, zsupi, &
          zqswp1, zqsip1, zqsw0)

  tmax = MAXVAL(zt(1:kproma,:))
  tmin = MINVAL(zt(1:kproma,:))
  zmlwc(1:kproma,:) = zxlb(1:kproma,:)
  zmiwc(1:kproma,:) = zxib(1:kproma,:)

  CALL riming_variables(&
          !--IN
          kbdim, klev, kproma, &
          zxi, zqirim, zbirim, &
          !--OUT
          zrimfrac, zrhop)

  CALL secondary_prognostics(&
          !--IN
          kbdim, klev, kproma, &
          zxi, zqihet, zqiliq, zqioliq, zqsrc, zqprc, &
          zicnc, znihet, znihom, zninuc, znidet, &
          !--OUT
          zhetfrac, zliqfrac, zoliqfrac, zsosifrac, &
          znihetfrac, znihomfrac, zninucfrac, znidetfrac)

  ! fraction of source and cloud mass
  ztmp1_2d(1:kproma,:) = zxi(1:kproma,:) + zxl(1:kproma,:)
  zsocmfrac(1:kproma,:) = MERGE(pxtm1(1:kproma,:,idt_qsrc)/ztmp1_2d(1:kproma,:), 0._dp, &
                                ztmp1_2d(1:kproma,:) > eps)

  CALL _secondary_ice_properties_p3(&
          !--IN
          kbdim, klev, kproma, &
          zrhop, zrimfrac, zxi, zicnc, &
          rhofaci, &
          !--OUT
          zrim, zrieff, zriv, &
          zvtin, zvtim, zrhoice, &
          zlkp_col, zlkp_slf, &
          zlkp_depx1, zlkp_depx2, &
          zlkp_upper, zlkp_lower)

  CALL secondary_ice_properties_2m(&
          !--IN
          kbdim, klev, kproma, &
          zrho, zxi, zicnc, &
          zaaa, &
          !--OUT
          zriv_2m, zrieff_2m, zvtim_2m, zvtin_2m)

  ! make sure mean ice size is in bounds (i.e. apply lambda limiters)
  ! if the lookup table used another value than picncb, we adjust to that one
  ! and act as if it was some kind of aggregation/multiplication
  CALL handle_overflow(&
              !--IN
              kbdim, klev, kproma, zstmst_rcp, &
              zaclci, zlkp_upper, zlkp_lower, &
              !--INOUT
              zicnc, zicncb, zni_lkp)
  zni_lkp(1:kproma,:) = MERGE(zni_lkp(1:kproma,:), 0._dp, zxi(1:kproma,:) > 0._dp)
  ! update tendency
  pxtte(1:kproma,:,idt_icnc) = pxtte(1:kproma,:,idt_icnc) + zni_lkp(1:kproma,:)*zntmst_rcp
  ! write diagnostics
  dni_lkp(1:kproma,:,krow) = dni_lkp(1:kproma,:,krow) + zdtime*zni_lkp(1:kproma,:)*zntmst_rcp

  ! diagnose liquid droplet size
  ztmp1_2d(1:kproma,:) = zxl(1:kproma,:)/MAX(zcdnc(1:kproma,:), cqtmin)
  ztmp1_2d(1:kproma,:) = MAX(ztmp1_2d(1:kproma,:), 0._dp)
  ! mass-weighted mean cloud droplet radius:
  zrcm(1:kproma,:)     = (ztmp1_2d(1:kproma,:)/(1.333*pi*rhoh2o))**0.333

END DO substep_column

IF(ltimer) CALL timer_stop(timer_micro)

if(l2moment) THEN
   zqirimflx(1:kproma) = 0._dp
   zbirimflx(1:kproma) = 0._dp
   zqirim(1:kproma,:) = 0._dp
   zbirim(1:kproma,:) = 0._dp
   pxtte(1:kproma,:,idt_qirim) = 0._dp
   pxtte(1:kproma,:,idt_birim) = 0._dp
   pxtm1(1:kproma,:,idt_qirim) = 0._dp
   pxtm1(1:kproma,:,idt_birim) = 0._dp
ENDIF

!-------------------------------------------------------------------------------------
!                                IN-CLOUD SCAVENGING
!-------------------------------------------------------------------------------------

  CALL cloud_subm_2( &

          !-- IN   
          kproma, kbdim, klev, ktdia, krow, & 
          zmlwc(:,:), zmiwc(:,:), zmratepr(:,:), zmrateps(:,:), &
          zfrain(:,:), zfsnow(:,:), zfevapr(:,:), zfsubls(:,:), & 
          zmsnowacl(:,:), paclc(:,:), ptm1(:,:), ptte(:,:), pxtm1(:,:,:), & 
          !-- INOUT
          pxtte(:,:,:), &
          !-- IN
          paphp1(:,:), papp1(:,:), zrho(:,:), zaclci(:,:) )

!-------------------------------------------------------------------------------------
!                                UPDATE COVER
!-------------------------------------------------------------------------------------

  ll1_2d(1:kproma,:) = (zxi(1:kproma,:) + zxl(1:kproma,:)) > cqtmin

  paclc(1:kproma,:) = MERGE(paclc(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))
  zaclci(1:kproma,:) = MERGE(zaclci(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))

  ll1_2d(1:kproma,:) = paclc(1:kproma,:) < clc_min
  paclc(1:kproma,:) = MERGE(0._dp, paclc(1:kproma,:), ll1_2d(1:kproma,:))
  ll1_2d(1:kproma,:) = zaclci(1:kproma,:) < clc_min
  zaclci(1:kproma,:) = MERGE(0._dp, zaclci(1:kproma,:), ll1_2d(1:kproma,:))

!-------------------------------------------------------------------------------------
!                                WRITE TENDENCIES
!-------------------------------------------------------------------------------------

!............................................................
! diagnostics

! read ice PSD params
  CALL get_ice_psd_params(&
          !--IN
          kbdim, klev, kproma, &
          zrhop, zrimfrac, zxi, zicnc, &
          rhofaci, &
          !--OUT
          zice_mu, zice_lam)

  dice_mu(1:kproma,:,krow) = dice_mu(1:kproma,:,krow) + zdtime*zice_mu(1:kproma,:)
  dice_lam(1:kproma,:,krow) = dice_lam(1:kproma,:,krow) + zdtime*zice_lam(1:kproma,:)

! PROCESS RATES
  dxisub_cv(1:kproma,:,krow) = dxisub_cv(1:kproma,:,krow) + zdtime*zxisub_cv(1:kproma,:)
  dximlt_cv(1:kproma,:,krow) = dximlt_cv(1:kproma,:,krow) + zdtime*zximlt_cv(1:kproma,:)
  dxlevp_cv(1:kproma,:,krow) = dxlevp_cv(1:kproma,:,krow) + zdtime*zxlevp_cv(1:kproma,:)
  dxisub_tp(1:kproma,:,krow) = dxisub_tp(1:kproma,:,krow) + zdtime*zxisub_tp(1:kproma,:)
  dnisub_tp(1:kproma,:,krow) = dnisub_tp(1:kproma,:,krow) + zdtime*znisub_tp(1:kproma,:)
  dxlevp_tp(1:kproma,:,krow) = dxlevp_tp(1:kproma,:,krow) + zdtime*zxlevp_tp(1:kproma,:)
  dncevp_tp(1:kproma,:,krow) = dncevp_tp(1:kproma,:,krow) + zdtime*zncevp_tp(1:kproma,:)
  dnifrz_ci(1:kproma,:,krow) = dnifrz_ci(1:kproma,:,krow) + zdtime*znifrz_ci(1:kproma,:)
  ddep_ci1(1:kproma,:,krow) = ddep_ci1(1:kproma,:,krow) + zdtime*zdep_ci1(1:kproma,:)

! store the changes due to transport (to diagnose the entire budget)
  dxite_ls(1:kproma,:,krow) = dxite_ls(1:kproma,:,krow) + zdtime*zaclci(1:kproma,:)*zxite_ls(1:kproma,:)
  dicncte_ls(1:kproma,:,krow) = dicncte_ls(1:kproma,:,krow) + zdtime*zaclci(1:kproma,:)*zicncte_ls(1:kproma,:)

  dxlte_ls(1:kproma,:,krow) = dxlte_ls(1:kproma,:,krow) + zdtime*paclc(1:kproma,:)*zxlte_ls(1:kproma,:)
  dcdncte_ls(1:kproma,:,krow) = dcdncte_ls(1:kproma,:,krow) + zdtime*paclc(1:kproma,:)*zcdncte_ls(1:kproma,:)

  dxite_cv(1:kproma,:,krow) = dxite_cv(1:kproma,:,krow) + zdtime*zxite_cv(1:kproma,:)
  dxlte_cv(1:kproma,:,krow) = dxlte_cv(1:kproma,:,krow) + zdtime*zxlte_cv(1:kproma,:)

  dicnc_cv(1:kproma,:,krow) = dicnc_cv(1:kproma,:,krow) + zdtime*zicnc_cv(1:kproma,:)
  dcdnc_cv(1:kproma,:,krow) = dcdnc_cv(1:kproma,:,krow) + zdtime*zcdnc_cv(1:kproma,:)
  dncnuc(1:kproma,:,krow) = dncnuc(1:kproma,:,krow) + zdtime*zncnuc(1:kproma,:)
  dcdncact(1:kproma,:,krow) = dcdncact(1:kproma,:,krow) + zdtime*zcdncact(1:kproma,:)

! PROPERTIES
  drcm(1:kproma,:,krow) = drcm(1:kproma,:,krow) + zdtime*zrcm(1:kproma,:)
  dxlb(1:kproma,:,krow) = dxlb(1:kproma,:,krow) + zdtime*zxlb(1:kproma,:)
  dcdnc(1:kproma,:,krow) = dcdnc(1:kproma,:,krow) + zdtime*zcdnc(1:kproma,:)
  dcdncb(1:kproma,:,krow) = dcdncb(1:kproma,:,krow) + zdtime*zcdncb(1:kproma,:)
  dxib(1:kproma,:,krow) = dxib(1:kproma,:,krow) + zdtime*zxib(1:kproma,:)
  dicnc(1:kproma,:,krow) = dicnc(1:kproma,:,krow) + zdtime*zicnc(1:kproma,:)
  dicncb(1:kproma,:,krow) = dicncb(1:kproma,:,krow) + zdtime*zicncb(1:kproma,:)
  dqirim(1:kproma,:,krow) = dqirim(1:kproma,:,krow) + zdtime*zqirim(1:kproma,:)
  dbirim(1:kproma,:,krow) = dbirim(1:kproma,:,krow) + zdtime*zbirim(1:kproma,:)
  dqihet(1:kproma,:,krow) = dqihet(1:kproma,:,krow) + zdtime*zqihet(1:kproma,:)
  dqiliq(1:kproma,:,krow) = dqiliq(1:kproma,:,krow) + zdtime*zqiliq(1:kproma,:)
  dqioliq(1:kproma,:,krow) = dqioliq(1:kproma,:,krow) + zdtime*zqioliq(1:kproma,:)
  dfr(1:kproma,:,krow) = dfr(1:kproma,:,krow) + zdtime*zrimfrac(1:kproma,:)
  dfr_het(1:kproma,:,krow) = dfr_het(1:kproma,:,krow) + zdtime*zhetfrac(1:kproma,:)
  dfr_liq(1:kproma,:,krow) = dfr_liq(1:kproma,:,krow) + zdtime*zliqfrac(1:kproma,:)
  dfr_oliq(1:kproma,:,krow) = dfr_oliq(1:kproma,:,krow) + zdtime*zoliqfrac(1:kproma,:)
  dfr_sosi(1:kproma,:,krow) = dfr_sosi(1:kproma,:,krow) + zdtime*zsosifrac(1:kproma,:)
  dfr_socm(1:kproma,:,krow) = dfr_socm(1:kproma,:,krow) + zdtime*zsocmfrac(1:kproma,:)
  dfr_nihet(1:kproma,:,krow) = dfr_nihet(1:kproma,:,krow) + zdtime*znihetfrac(1:kproma,:)
  dfr_nihom(1:kproma,:,krow) = dfr_nihom(1:kproma,:,krow) + zdtime*znihomfrac(1:kproma,:)
  dfr_ninuc(1:kproma,:,krow) = dfr_ninuc(1:kproma,:,krow) + zdtime*zninucfrac(1:kproma,:)
  dfr_nidet(1:kproma,:,krow) = dfr_nidet(1:kproma,:,krow) + zdtime*znidetfrac(1:kproma,:)
  drhop(1:kproma,:,krow) = drhop(1:kproma,:,krow) + zdtime*zrhop(1:kproma,:)
  drhoice(1:kproma,:,krow) = drhoice(1:kproma,:,krow) + zdtime*zrhoice(1:kproma,:)

  dnihet(1:kproma,:,krow) = dnihet(1:kproma,:,krow) + zdtime*znihet(1:kproma,:)
  dnihom(1:kproma,:,krow) = dnihom(1:kproma,:,krow) + zdtime*znihom(1:kproma,:)
  dninuc(1:kproma,:,krow) = dninuc(1:kproma,:,krow) + zdtime*zninuc(1:kproma,:)
  dnidet(1:kproma,:,krow) = dnidet(1:kproma,:,krow) + zdtime*znidet(1:kproma,:)

  daclci(1:kproma,:,krow) = daclci(1:kproma,:,krow) + zdtime*zaclci(1:kproma,:)

  dtke(1:kproma,:,krow) = dtke(1:kproma,:,krow) + zdtime*ptkem1(1:kproma,:)
  dupdraft(1:kproma,:,krow) = dupdraft(1:kproma,:,krow) + zdtime*zvervx(1:kproma,:)
  dupdraftmax(1:kproma,:,krow) = dupdraftmax(1:kproma,:,krow) &
                               + 100._dp*zdtime*zvervmax(1:kproma,:)
  drhoair(1:kproma,:,krow) = drhoair(1:kproma,:,krow) + zdtime*zrho(1:kproma,:)
  ddpg(1:kproma,:,krow) = ddpg(1:kproma,:,krow) + zdtime*zdpg(1:kproma,:)
  ddz(1:kproma,:,krow) = ddz(1:kproma,:,krow) + zdtime*zdz(1:kproma,:)
  dsedtm(1:kproma,:,krow) = dsedtm(1:kproma,:,krow) + zdtime*zsedtm(1:kproma,:)
  dsedtn(1:kproma,:,krow) = dsedtn(1:kproma,:,krow) + zdtime*zsedtn(1:kproma,:)
  ztmp1(1:kproma) = 1._dp
  dnmicro(1:kproma,krow) = dnmicro(1:kproma,krow) + zdtime*zntmst*ztmp1(1:kproma)

! scavenging
  dfrain(1:kproma,:,krow) = dfrain(1:kproma,:,krow) + zdtime*zfrain(1:kproma,:)
  dfsnow(1:kproma,:,krow) = dfsnow(1:kproma,:,krow) + zdtime*zfsnow(1:kproma,:)
  dfevapr(1:kproma,:,krow) = dfevapr(1:kproma,:,krow) + zdtime*zfevapr(1:kproma,:)
  dfsubls(1:kproma,:,krow) = dfsubls(1:kproma,:,krow) + zdtime*zfsubls(1:kproma,:)
  dmsnowacl(1:kproma,:,krow) = dmsnowacl(1:kproma,:,krow) + zdtime*zmsnowacl(1:kproma,:)
  dmratepr(1:kproma,:,krow) = dmratepr(1:kproma,:,krow) + zdtime*zmratepr(1:kproma,:)
  dmrateps(1:kproma,:,krow) = dmrateps(1:kproma,:,krow) + zdtime*zmrateps(1:kproma,:)
  dmlwc(1:kproma,:,krow) = dmlwc(1:kproma,:,krow) + zdtime*zmlwc(1:kproma,:)
  dmiwc(1:kproma,:,krow) = dmiwc(1:kproma,:,krow) + zdtime*zmiwc(1:kproma,:)

  lo2_2d(:,:) = ice_switch(kbdim, kproma, klev, csecfrl, zt, zxi, zvervx, zvervmax)

  zqs(1:kproma,:) = MERGE(zqsi(1:kproma,:), zqsw(1:kproma,:),lo2_2d(1:kproma,:))
  prelhum(1:kproma,:) = zq(1:kproma,:)/zqs(1:kproma,:)


!-------------------------------------------------------------------------------------
!                                STANDARD DIAGNOSTICS
!-------------------------------------------------------------------------------------
  
! ice and liquid switches
  ll1_2d(1:kproma,:) = zxl(1:kproma,:) > cqtmin .AND. &
                       zcdnc(1:kproma,:) > cqtmin
  ll2_2d(1:kproma,:) = zxi(1:kproma,:) > cqtmin .AND. &
                       zicnc(1:kproma,:) > cqtmin

! LIQUID CLOUD TIME
  ztmp1_2d(1:kproma,:)        = cloud_time(1:kproma,:,krow) + zdtime
  cloud_time(1:kproma,:,krow) = MERGE(ztmp1_2d(1:kproma,:), cloud_time(1:kproma,:,krow), ll1_2d(1:kproma,:))

! ICE CLOUD TIME
  ztmp1_2d(1:kproma,:)        = cliwc_time(1:kproma,:,krow) + zdtime
  cliwc_time(1:kproma,:,krow) = MERGE(ztmp1_2d(1:kproma,:), cliwc_time(1:kproma,:,krow), ll2_2d(1:kproma,:))

! INTEGRATED NUMBERS
  DO jk = ktdia,klev
     ztmp1(1:kproma) = cdnc_burden(1:kproma,krow)+zdtime*zcdnc(1:kproma,jk)*zdpg(1:kproma,jk)
     cdnc_burden(1:kproma,krow) = MERGE(ztmp1(1:kproma), cdnc_burden(1:kproma,krow), ll1_2d(1:kproma,jk))

     ztmp1(1:kproma) = icnc_burden(1:kproma,krow)+zdtime*zicnc(1:kproma,jk)*zdpg(1:kproma,jk)
     icnc_burden(1:kproma,krow) = MERGE(ztmp1(1:kproma), icnc_burden(1:kproma,krow), ll2_2d(1:kproma,jk))
  END DO !jk

! CLOUD WEIGHTED NUMBERS
  ztmp1_2d(1:kproma,:) = cdnc(1:kproma,:,krow) &
                       + zdtime*zcdnc(1:kproma,:)*zrho(1:kproma,:)*paclc(1:kproma,:)
  cdnc(1:kproma,:,krow) = MERGE(ztmp1_2d(1:kproma,:), cdnc(1:kproma,:,krow), ll1_2d(1:kproma,:))
  ztmp1_2d(1:kproma,:) = icnc(1:kproma,:,krow) &
                       + zdtime*zicnc(1:kproma,:)*zrho(1:kproma,:)*zaclci(1:kproma,:)
  icnc(1:kproma,:,krow) = MERGE(ztmp1_2d(1:kproma,:), icnc(1:kproma,:,krow), ll2_2d(1:kproma,:))


! TOTAL ACCUMULATED NUMBERS
  ztmp1_2d(1:kproma,:) = cdnc_acc(1:kproma,:,krow) &
                       + zdtime*zcdnc(1:kproma,:)*zrho(1:kproma,:)
  cdnc_acc(1:kproma,:,krow) = MERGE(ztmp1_2d(1:kproma,:), cdnc_acc(1:kproma,:,krow), ll1_2d(1:kproma,:))
  ztmp1_2d(1:kproma,:) = icnc_acc(1:kproma,:,krow) &
                       + zdtime*zicnc(1:kproma,:)*zrho(1:kproma,:)
  icnc_acc(1:kproma,:,krow) = MERGE(ztmp1_2d(1:kproma,:), icnc_acc(1:kproma,:,krow), ll2_2d(1:kproma,:))

! CLOUD TOP NUMBER AND EFFECTIVE RADII
!>>DN bugfix
!  zrleff(1:kproma,:) = effective_liquid_radius(kbdim, kproma, klev, zcdncb, zxlb)
  zrleff(:,:) = effective_liquid_radius(kbdim, kproma, klev, zcdncb, zxlb)
!<<DN bugfix
  reffl(1:kproma,:,krow) = zrleff(1:kproma,:)
  reffl_acc(1:kproma,:,krow) = reffl_acc(1:kproma,:,krow) + zdtime*zrleff(1:kproma,:)
  reffi(1:kproma,:,krow) = zrieff(1:kproma,:)
  reffi_acc(1:kproma,:,krow) = reffi_acc(1:kproma,:,krow) + zdtime*zrieff(1:kproma,:)
  DO jk=ktdia,klev
     ll3(1:kproma) = ll1_2d(1:kproma,jk) &
               .AND. (itop(1:kproma,jk)  == jk ) &
               .AND. (zt(1:kproma,jk)   >  tmelt) &
               .AND. (zreffct(1:kproma)   <  4._dp) &
               .AND. (zrleff(1:kproma,jk) >= 4._dp)

     ztmp1(1:kproma) = reffl_ct(1:kproma,krow) + zdtime*zrleff(1:kproma,jk)
     reffl_ct(1:kproma,krow) = MERGE(ztmp1(1:kproma), reffl_ct(1:kproma,krow), ll3(1:kproma))

     ztmp1(1:kproma)        = cdnc_ct(1:kproma,krow) + zdtime*zcdnc(1:kproma,jk)*paclc(1:kproma,jk)
     cdnc_ct(1:kproma,krow) = MERGE(ztmp1(1:kproma), cdnc_ct(1:kproma,krow), ll3(1:kproma))

     ztmp1(1:kproma)           = reffl_time(1:kproma,krow) + zdtime
     reffl_time(1:kproma,krow) = MERGE(ztmp1(1:kproma), reffl_time(1:kproma,krow), ll3(1:kproma))

     zreffct(1:kproma) = MERGE(zrleff(1:kproma,jk), zreffct(1:kproma), ll3(1:kproma))
  END DO !jk

! NUMBERS OUT
  pacdnc(1:kproma,:) = MAX(1._dp, zcdnc(1:kproma,:)*zrho(1:kproma,:))
  picnc(1:kproma,:) = MAX(1._dp, zicnc(1:kproma,:)*zrho(1:kproma,:))


 !       10.2   Total cloud cover
 !
  zclcov(1:kproma) = 1.0_dp-paclc(1:kproma,1)

  DO 923 jk = 2,klev
     ztmp1(1:kproma) = MAX(paclc(1:kproma,jk), paclc(1:kproma,jk-1))
     ztmp2(1:kproma) = MIN(paclc(1:kproma,jk-1), xsec)

     zclcov(1:kproma) = zclcov(1:kproma) * (1._dp - ztmp1(1:kproma))  &
                      / (1._dp - ztmp2(1:kproma))
923 END DO

  zclcov(1:kproma)  = 1.0_dp-zclcov(1:kproma)
  paclcov(1:kproma) = paclcov(1:kproma) + zdtime*zclcov(1:kproma)

! vertical integrals
  zxi_cloud(1:kproma,:) = MERGE(zxi(1:kproma,:), 0._dp, ll_cc(1:kproma,:))
  zxi_snow(1:kproma,:) = MERGE(0._dp, zxi(1:kproma,:), ll_cc(1:kproma,:))
  dxi_cloud(1:kproma,:,krow) = dxi_cloud(1:kproma,:,krow) + zdtime*zxi_cloud(1:kproma,:)
  dxi_snow(1:kproma,:,krow) = dxi_snow(1:kproma,:,krow) + zdtime*zxi_snow(1:kproma,:)

  zqvi(1:kproma)  = 0.0_dp
  zxlvi(1:kproma) = 0.0_dp
  zxivi(1:kproma) = 0.0_dp
   !
  DO 933 jk = ktdia,klev
     zqvi(1:kproma)  = zqvi(1:kproma)  + zq(1:kproma,jk) *zdpg(1:kproma,jk)
     zxlvi(1:kproma) = zxlvi(1:kproma) + zxl(1:kproma,jk)*zdpg(1:kproma,jk) 
     zxivi(1:kproma) = zxivi(1:kproma) + zxi(1:kproma,jk)*zdpg(1:kproma,jk)
     dmete(1:kproma,krow) = dmete(1:kproma,krow) + &
                                (ptte(1:kproma,jk) - ztte_0(1:kproma,jk))&
                                *(cpd+zcons1*MAX(pqm1(1:kproma,jk),0.0_dp))*zdtime*zdpg(1:kproma,jk)
     dncdiff(1:kproma,krow) = dncdiff(1:kproma,krow) &
                            + (pxtte(1:kproma,jk,idt_cdnc)-zdnc_0(1:kproma,jk))*zdtime*zdpg(1:kproma,jk)
     dnidiff(1:kproma,krow) = dnidiff(1:kproma,krow) &
                            + (pxtte(1:kproma,jk,idt_icnc)-zdni_0(1:kproma,jk))*zdtime*zdpg(1:kproma,jk)
     dncdiff2(1:kproma,krow) = dncdiff2(1:kproma,krow) &
                            + (zcdnc(1:kproma,jk)-zdnc2_0(1:kproma,jk))*zdpg(1:kproma,jk)*zdt
     dnidiff2(1:kproma,krow) = dnidiff2(1:kproma,krow) &
                            + (zicnc(1:kproma,jk)-zdni2_0(1:kproma,jk))*zdpg(1:kproma,jk)*zdt
     ztmp1(1:kproma) = (pxtte(1:kproma,jk,idt_icnc)-zdni_0(1:kproma,jk))*zdtime*zdpg(1:kproma,jk)
     ztmp2(1:kproma) = (zicnc(1:kproma,jk)-zdni2_0(1:kproma,jk))*zdt*zdpg(1:kproma,jk)
     dnictrl(1:kproma,krow) = dnictrl(1:kproma,krow) + (ztmp2(1:kproma) - ztmp1(1:kproma))

     ztmp1(1:kproma) = (pxtte(1:kproma,jk,idt_cdnc)-zdnc_0(1:kproma,jk))*zdtime*zdpg(1:kproma,jk)
     ztmp2(1:kproma) = (zcdnc(1:kproma,jk)-zdnc2_0(1:kproma,jk))*zdpg(1:kproma,jk)*zdt
     dncctrl(1:kproma,krow) = dncctrl(1:kproma,krow) + (ztmp2(1:kproma) - ztmp1(1:kproma))
     
     dxivi_cloud(1:kproma,krow) = dxivi_cloud(1:kproma,krow) + zdtime*zxi_cloud(1:kproma,jk)*zdpg(1:kproma,jk)
     dxivi_snow(1:kproma,krow) = dxivi_snow(1:kproma,krow) + zdtime*zxi_snow(1:kproma,jk)*zdpg(1:kproma,jk)
933 END DO
  pqvi(1:kproma)  = pqvi(1:kproma)  + zdtime*zqvi(1:kproma)
  pxlvi(1:kproma) = pxlvi(1:kproma) + zdtime*zxlvi(1:kproma)
  pxivi(1:kproma) = pxivi(1:kproma) + zdtime*zxivi(1:kproma)

  paprl(1:kproma) = paprl(1:kproma) + zdtime * (prsfl(1:kproma)+pssfl(1:kproma))
  paprs(1:kproma) = paprs(1:kproma) + zdtime * pssfl(1:kproma)

  dnisurf(1:kproma,krow) = dnisurf(1:kproma,krow) + zdtime*znisurf(1:kproma)

  paclcac(1:kproma,:) = paclcac(1:kproma,:) + zdtime*paclc(1:kproma,:)

  icnt = icnt + 1
  
  !-- Some more data exchange to other model parts
  CALL cloud_subm_3( &
       !-- IN   
       kproma, kbdim, klev, krow, itop(:,:), ll1_2d(:,:), &
       zclcov(:), zcdnc(:,:), zrho(:,:), paclc(:,:), zrleff(:,:), prsfl(:), pssfl(:), &
       zaut_stm(:,:), zacc_stm(:,:), zdpg(:,:), zqccol_stm(:,:), zfrain(:,:), zfevapr(:,:),&
       zaclci(:,:))

END SUBROUTINE cloud_micro_interface_p3

SUBROUTINE cloud_formation(&
              !--IN
              kbdim, kproma, ptmst, ptmst_rcp, &
              pxl, pxi, pcdnc, picnc, pt, pq, &
              plvdcp, plsdcp, pqte, ptte, pdep, &
              pqsi, pqsip1, pqsw, pqswp1, psupi, paclc, &
              !--OUT
              pdep_co, pcnd_co, pnsub_co, pnevp_co, &
              pdep_wbf, pcdnc_wbf, pqcdif)

  INTEGER, INTENT(IN) :: kbdim, kproma
  REAL(dp), INTENT(IN) :: ptmst, ptmst_rcp

  REAL(dp), INTENT(IN), DIMENSION(kbdim) :: pxl, pxi, pcdnc, picnc, pt, pq, &
                                            plvdcp, plsdcp, pqte, ptte, pdep, &
                                            pqsi, pqsip1, pqsw, pqswp1, psupi, paclc

  REAL(dp), INTENT(OUT), DIMENSION(kbdim) :: pdep_co, pcnd_co, pnsub_co, pnevp_co, &
                                             pdep_wbf, pcdnc_wbf, pqcdif

  REAL(dp), DIMENSION(kbdim) :: zmix, zmixp1, zlc, zqs, zqsp1, &
                                zdqsdt, zdqsat, zdep, zdep_vap
  REAL(dp), DIMENSION(kbdim) :: ztmp1, ztmp2, ztmp3, ztmp4
  LOGICAL, DIMENSION(kbdim) :: ll1, ll2

! use water saturation throughout
  zlc(1:kproma)   = plvdcp(1:kproma)
  zqs(1:kproma)   = pqsw(1:kproma)
  zqsp1(1:kproma) = pqswp1(1:kproma)

  zdqsdt(1:kproma) = 1000._dp*(zqsp1(1:kproma)-zqs(1:kproma))

! calculate change in saturation specific humidity
  zdqsat(1:kproma) = ptte(1:kproma)*ptmst &
       + paclc(1:kproma)*ptmst*zlc(1:kproma)*pqte(1:kproma)

  zdqsat(1:kproma) = zdqsat(1:kproma) * zdqsdt(1:kproma) &
                     / (1._dp+paclc(1:kproma)*zlc(1:kproma)*zdqsdt(1:kproma))

! amount of condensation/deposition that is required for cloud scheme
  pqcdif(1:kproma) = (ptmst*pqte(1:kproma)-zdqsat(1:kproma))*paclc(1:kproma)

! amount missing to reach ice saturation (positive if supersaturated)
  zdep_vap(1:kproma) = (pq(1:kproma)-pqsi(1:kproma))/  &
                       (1+als**2*pqsi(1:kproma)/(cpd*rv*pt(1:kproma)**2))
  zdep_vap(1:kproma) = MERGE(MAX(0._dp, zdep_vap(1:kproma)), 0._dp, lallow_si)

  ! divergence
  ll1(1:kproma) = (pqcdif(1:kproma) < 0._dp)

  ! Deposition term only needs to be considered when ice is growing
  zdep(1:kproma) = MAX(0._dp, pdep(1:kproma))
  ! convert pdep to mass per timestep
  zdep(1:kproma) = ptmst*zdep(1:kproma)

  ! calculation of the rates according to Morrison-Gettelman (2008) scheme
  ztmp1(1:kproma) = MIN(paclc(1:kproma)*zdep(1:kproma),pqcdif(1:kproma)+zdep_vap(1:kproma)) !< (dq_i/dt)_dep
  ztmp2(1:kproma) = MAX(pqcdif(1:kproma)-ztmp1(1:kproma),0._dp)                             !< (dq_c/dt)_cond
  ztmp4(1:kproma) = MAX(pqcdif(1:kproma),-pxl(1:kproma))                                    !< (dq_c/dt)_evap
  ztmp3(1:kproma) = MAX(pqcdif(1:kproma)-ztmp4(1:kproma),-pxi(1:kproma))                    !< (dq_i/dt)_sub

  pdep_co(1:kproma) = MERGE(ztmp3(1:kproma),ztmp1(1:kproma),ll1(1:kproma)) 
  pcnd_co(1:kproma) = MERGE(ztmp4(1:kproma),ztmp2(1:kproma),ll1(1:kproma))

  ztmp1(1:kproma) = MIN(pxl(1:kproma), &
                         MAX(paclc(1:kproma)*zdep(1:kproma)-pqcdif(1:kproma),0._dp))
  pdep_wbf(1:kproma) = MERGE(0._dp, ztmp1(1:kproma), ll1(1:kproma))

  ! convergence deposition only in mixed phase as this will not allow RH_i > 1
  ! For cirrus we use the explicit deposition rate as f(RH_i) and only apply
  ! the sublimation case here.
  ll2(1:kproma) = (pt(1:kproma) > cthomi) .OR. &      ! mixed phase form/diss
                  (pt(1:kproma) <= cthomi .AND. &      ! cirrus dissipation
                   ll1(1:kproma) .AND. psupi(1:kproma) < 0._dp)

  pdep_co(1:kproma) = MERGE(pdep_co(1:kproma), 0._dp, ll2(1:kproma))
  pcnd_co(1:kproma) = MERGE(pcnd_co(1:kproma), 0._dp, pt(1:kproma) > cthomi)

  ! calculate the associated number changes
  ! liquid
  ll2(1:kproma) = pxl(1:kproma) > cqtmin
  ztmp1(1:kproma) = -MIN(0._dp, pcnd_co(1:kproma))
  ztmp1(1:kproma) = MERGE(ztmp1(1:kproma)/pxl(1:kproma), 0._dp, ll2(1:kproma))
  ztmp1(1:kproma) = MIN(1._dp, ztmp1(1:kproma))
  pnevp_co(1:kproma) = pcdnc(1:kproma)*ztmp1(1:kproma)

  ! due to WBF
  ztmp1(1:kproma) = MERGE(pdep_wbf(1:kproma)/pxl(1:kproma), 0._dp, ll2(1:kproma))
  ztmp1(1:kproma) = MIN(1._dp, ztmp1(1:kproma))
  pcdnc_wbf(1:kproma) = pcdnc(1:kproma)*ztmp1(1:kproma)

  ! ice
  ll2(1:kproma) = pxi(1:kproma) > cqtmin
  ztmp1(1:kproma) = -MIN(0._dp, pdep_co(1:kproma))
  ztmp1(1:kproma) = MERGE(ztmp1(1:kproma)/pxi(1:kproma), 0._dp, ll2(1:kproma))
  ztmp1(1:kproma) = MIN(1._dp, ztmp1(1:kproma))
  pnsub_co(1:kproma) = picnc(1:kproma)*ztmp1(1:kproma)

  ! convert to a rate
  pcnd_co(1:kproma) = ptmst_rcp*pcnd_co(1:kproma)
  pdep_co(1:kproma) = ptmst_rcp*pdep_co(1:kproma)
  pnevp_co(1:kproma) = ptmst_rcp*pnevp_co(1:kproma)
  pnsub_co(1:kproma) = ptmst_rcp*pnsub_co(1:kproma)
  pdep_wbf(1:kproma) = ptmst_rcp*pdep_wbf(1:kproma)
  pcdnc_wbf(1:kproma) = ptmst_rcp*pcdnc_wbf(1:kproma)

END SUBROUTINE cloud_formation

FUNCTION conservation_reduction_2d(kbdim, klev, kproma, psrcs, psnks) RESULT(pfrac)

  INTEGER :: kbdim, klev, kproma
  REAL(dp) :: ptmst

  REAL(dp), DIMENSION(kbdim,klev):: psrcs, psnks
  REAL(dp), DIMENSION(kbdim,klev):: pfrac

  LOGICAL :: ll1(kbdim,klev)

  ll1(1:kproma,:) = (psnks(1:kproma,:) .gt. psrcs(1:kproma,:)) .and.                     &
                    (psnks(1:kproma,:) .gt. epsec)
  pfrac(1:kproma,:) = MERGE(psrcs(1:kproma,:)/psnks(1:kproma,:), 1._dp, ll1(1:kproma,:))

END FUNCTION conservation_reduction_2d

FUNCTION conservation_reduction_1d(kbdim, kproma, psrcs, psnks) RESULT(pfrac)

  INTEGER :: kbdim, kproma
  REAL(dp) :: ptmst

  REAL(dp), DIMENSION(kbdim) :: psrcs, psnks
  REAL(dp), DIMENSION(kbdim) :: pfrac

  LOGICAL :: ll1(kbdim)

  ll1(1:kproma) = (psnks(1:kproma) .gt. psrcs(1:kproma)) .and.                     &
                  (psnks(1:kproma) .gt. epsec)
  pfrac(1:kproma) = MERGE(psrcs(1:kproma)/psnks(1:kproma), 1._dp, ll1(1:kproma))

END FUNCTION conservation_reduction_1d 

SUBROUTINE sublimation_snow(&
              !--IN
              kbdim, kproma, pdpg, &
              pq, pt, paclci, ll_cci, pqsflx, &
              psupi, pqsi, prho_rcp, plsdcp, &
              !--OUT
              pqssub)
  
  INTEGER, INTENT(IN) :: kbdim, kproma
  LOGICAL, INTENT(IN), DIMENSION(kbdim)  :: ll_cci
  REAL(dp), INTENT(IN), DIMENSION(kbdim) :: pt, pq, paclci, pqsflx, psupi, &
                                            pqsi, prho_rcp, plsdcp, pdpg
  REAL(dp), INTENT(OUT), DIMENSION(kbdim) :: pqssub

  REAL(dp), DIMENSION(kbdim) :: zqrho, zsubi, zclcpre
  REAL(dp), DIMENSION(kbdim) :: ztmp1, ztmp2, ztmp3, ztmp4
  LOGICAL, DIMENSION(kbdim) :: ll1

  zqrho(1:kproma) = 1.3_dp*prho_rcp(1:kproma)
  zsubi(1:kproma) = -MIN(0._dp, psupi)
  zclcpre(1:kproma) = MERGE(paclci(1:kproma), 1._dp, ll_cci(1:kproma))

  ll1(1:kproma) = (pqsflx(1:kproma) > cqtmin) .AND. &
                   ll_cci(1:kproma)

  ! sublimation parametrization
  ztmp1(1:kproma) = 1._dp/(2.43e-2_dp*rv)*plsdcp(1:kproma)**2/pt(1:kproma)**2
  ztmp1(1:kproma) = ztmp1(1:kproma) &
                  + 1._dp/0.211e-4_dp*prho_rcp(1:kproma)/pqsi(1:kproma)
  ztmp1(1:kproma) = 3.e6_dp*2._dp*pi*zsubi(1:kproma)*prho_rcp(1:kproma)/ztmp1(1:kproma)
  
  ztmp2(1:kproma) = zcons3*(MAX(cqtmin,pqsflx(1:kproma)) / zclcpre(1:kproma))**(0.25_dp/1.16_dp)
  ztmp2(1:kproma) = 0.78_dp  *ztmp2(1:kproma)**2 & 
                  + 232.19_dp*zqrho(1:kproma)**0.25_dp * ztmp2(1:kproma)**2.625_dp

  ! plug together
  pqssub(1:kproma) = ztmp2(1:kproma) * ztmp1(1:kproma)
                    
  ! check bounds.
  pqssub(1:kproma) = MIN(pqsflx(1:kproma)/pdpg(1:kproma), pqssub(1:kproma))
  pqssub(1:kproma)  = MERGE(pqssub(1:kproma), 0._dp, ll1(1:kproma))

END SUBROUTINE sublimation_snow

SUBROUTINE handle_overflow(&
              !--IN
              kbdim, klev, kproma, ptmst_rcp, &
              paclci, lkp_upper, lkp_lower, &
              !--INOUT
              picnc, picncb, pni_lkp)

  INTEGER, INTENT(IN)                            :: kbdim, klev, kproma
  REAL(dp), INTENT(IN)                           :: ptmst_rcp 
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev)    :: paclci, lkp_upper, lkp_lower
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim,klev) :: picnc, picncb
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev)   :: pni_lkp

  REAL(dp), DIMENSION(kbdim,klev) :: dum, delta

! compute under-/overflow
  dum(1:kproma,:) = MAX(MIN(lkp_upper(1:kproma,:), picncb(1:kproma,:)), lkp_lower(1:kproma,:))
  delta(1:kproma,:) = dum(1:kproma,:) - picncb(1:kproma,:)

! adjust numbers
  picnc(1:kproma,:) = picnc(1:kproma,:) + delta(1:kproma,:)*paclci(1:kproma,:)
  picncb(1:kproma,:) = picncb(1:kproma,:) + delta(1:kproma,:)

! save change for next update (to avoid double counting, this may be arranged in a better way...)
  pni_lkp(1:kproma,:) = delta(1:kproma,:)*paclci(1:kproma,:)*ptmst_rcp

END SUBROUTINE handle_overflow

SUBROUTINE get_ice_psd_params(&
              !--IN
              kbdim, klev, kproma, &
              prhop, primfrac, pxi, picnc, &
              prhofaci, &
              !--OUT
              pmu, plam)

  INTEGER, INTENT(IN) :: kbdim,klev, kproma
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: prhop, primfrac, pxi, picnc
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: prhofaci

  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: pmu, plam

  REAL(dp), DIMENSION(kbdim,klev) :: qnorm
  REAL(dp), DIMENSION(kbdim,klev) :: iqnorm, irhor, ifr
  REAL(dp), DIMENSION(kbdim,klev) :: zlkp1, zlkp2
  LOGICAL, DIMENSION(kbdim,klev) :: ll_ice

!>>DN bugfix
!  qnorm = MERGE(pxi(1:kproma,:)/picnc(1:kproma,:), 0._dp, picnc(1:kproma,:) > cqtmin)
  qnorm(1:kproma,:) = MERGE(pxi(1:kproma,:)/picnc(1:kproma,:), 0._dp, picnc(1:kproma,:) > cqtmin)
!<<DN bugfix

  CALL calculate_my_lookup_table_indices(&
          !--IN
          kproma, kbdim, klev, &
          qnorm, primfrac, prhop, &
          !--OUT
          iqnorm, irhor, ifr)

  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 1, &
          !--OUT
          zlkp1)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 2, &
          !--OUT
          zlkp2)

! ice present logical
  ll_ice(1:kproma,:) = pxi(1:kproma,:) .gt. cqtmin

  plam(1:kproma,:) = MERGE(zlkp1(1:kproma,:), -1._dp, ll_ice(1:kproma,:))
  pmu(1:kproma,:) = MERGE(zlkp2(1:kproma,:), -1._dp, ll_ice(1:kproma,:))

END SUBROUTINE get_ice_psd_params

SUBROUTINE my_secondary_ice_properties_p3(&
              !--IN
              kbdim, klev, kproma, &
              prhop, primfrac, pxi, picnc, &
              prhofaci, &
              !--OUT
              prim, prieff, priv, &
              pvtin, pvtim, prhoice, &
              plkp_col, plkp_slf, &
              plkp_depx1, plkp_depx2, &
              plkp_upper, plkp_lower)

  INTEGER, INTENT(IN) :: kbdim,klev, kproma
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: prhop, primfrac, pxi, picnc
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: prhofaci

  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: prim, prieff, priv, pvtin, &
                                                 pvtim, prhoice
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: plkp_upper, plkp_lower
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: plkp_col, plkp_slf
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: plkp_depx1, plkp_depx2

  REAL(dp), DIMENSION(kbdim,klev) :: qnorm
  REAL(dp), DIMENSION(kbdim,klev) :: iqnorm, irhor, ifr
  REAL(dp), DIMENSION(kbdim,klev) :: zlkp3, zlkp5, zlkp6, &
                                     zlkp7, zlkp8, zlkp10
  LOGICAL, DIMENSION(kbdim,klev) :: ll_ice

!>>DN bugfix
!  qnorm = MERGE(pxi(1:kproma,:)/picnc(1:kproma,:), 0._dp, picnc(1:kproma,:) > cqtmin)
  qnorm(1:kproma,:) = MERGE(pxi(1:kproma,:)/picnc(1:kproma,:), 0._dp, picnc(1:kproma,:) > cqtmin)
!<<DN bugfix

  CALL get_ni_limits(&
          !--IN
          kproma, kbdim, klev, &
          pxi, picnc, &
          !--OUT
          plkp_lower, plkp_upper)

  CALL calculate_my_lookup_table_indices(&
          !--IN
          kproma, kbdim, klev, &
          qnorm, primfrac, prhop, &
          !--OUT
          iqnorm, irhor, ifr)

  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 3, &
          !--OUT
          zlkp3)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 5, &
          !--OUT
          zlkp5)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 6, &
          !--OUT
          zlkp6)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 7, &
          !--OUT
          zlkp7)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 8, &
          !--OUT
          plkp_slf)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 9, &
          !--OUT
          plkp_col)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 10, &
          !--OUT
          zlkp10)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 11, &
          !--OUT
          plkp_depx1)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 12, &
          !--OUT
          plkp_depx2)

! ice present logical
  ll_ice(1:kproma,:) = pxi(1:kproma,:) .gt. cqtmin

  prim(1:kproma,:) = MERGE(0.5*zlkp3(1:kproma,:), 0._dp, ll_ice(1:kproma,:))
  prieff(1:kproma,:) = MERGE(zlkp5(1:kproma,:), 0._dp, ll_ice(1:kproma,:))
  priv(:,:) = effective_2_volmean_radius_param_Schuman_2011( &
                         kbdim, kproma, klev, prieff(:,:))
  priv(1:kproma,:) = MERGE(priv(1:kproma,:), 0._dp, ll_ice(1:kproma,:))
  pvtim(1:kproma,:) = MERGE(prhofaci(1:kproma,:)*zlkp6(1:kproma,:)*vt_rdc, 0._dp, ll_ice(1:kproma,:))
  pvtin(1:kproma,:) = MERGE(prhofaci(1:kproma,:)*zlkp7(1:kproma,:)*vt_rdc, 0._dp, ll_ice(1:kproma,:))
  prhoice(1:kproma,:) = MERGE(zlkp10(1:kproma,:), 0._dp, ll_ice(1:kproma,:))

  plkp_slf(1:kproma,:) = plkp_slf(1:kproma,:)*picnc(1:kproma,:)**2
  plkp_col(1:kproma,:) = plkp_col(1:kproma,:)*picnc(1:kproma,:)

END SUBROUTINE my_secondary_ice_properties_p3

SUBROUTINE secondary_ice_properties_p3(&
              !--IN
              kbdim, klev, kproma, &
              prhop, primfrac, pxi, picnc, &
              prhofaci, &
              !--OUT
              prim, prieff, priv, &
              pvtin, pvtim, prhoice, &
              plkp_col, plkp_slf, &
              plkp_depx1, plkp_depx2, &
              plkp_upper, plkp_lower)

  USE mo_p3_fields, ONLY: vt_rdc

  INTEGER, INTENT(IN) :: kbdim,klev, kproma
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: prhop, primfrac, pxi, picnc
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: prhofaci

  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: prim, prieff, priv, pvtin, &
                                                 pvtim, prhoice
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: plkp_upper, plkp_lower
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: plkp_col, plkp_slf
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: plkp_depx1, plkp_depx2

  REAL(dp), DIMENSION(kbdim,klev) :: dum1, dum2, dum4, dum5
  INTEGER, DIMENSION(kbdim,klev) :: dumi, dumk, dumii, dumjj
  REAL(dp), DIMENSION(kbdim,klev) :: zlkp1, zlkp2, zlkp6, &
                                     zlkp7, zlkp8, zlkp11, zlkp12, &
                                     zlkp5, zlkp10
  LOGICAL, DIMENSION(kbdim,klev)  :: ll_ice
  INTEGER :: jk, jl

!-------------------------------------------------------------------------------------
!                READ ICE MICROPHYSICAL VALUES FROM LOOKUP TABLE
!-------------------------------------------------------------------------------------

  CALL calculate_lookup_table_indices(&
       ! IN
       kbdim,klev,kproma,                                                              &
       prhop(:,:),primfrac(:,:),pxi(:,:),picnc(:,:),     &
       ! OUT
       dum1(:,:), dum2(:,:), dum4(:,:), dum5(:,:),      &
       dumi(:,:), dumk(:,:), dumii(:,:), dumjj(:,:))

  DO jl=1,kproma
    DO jk=1,klev
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),1,    &
           ! OUT
           zlkp1(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),2,    &
           ! OUT
           zlkp2(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),3,    &
           ! OUT
           plkp_slf(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),4,    &
           ! OUT
           plkp_col(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),6,    &
           ! OUT
           zlkp6(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),7,    &
           ! OUT
           plkp_upper(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),8,    &
           ! OUT
           plkp_lower(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),11,    &
           ! OUT
           zlkp11(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),12,    &
           ! OUT
           zlkp12(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),10,    &
           ! OUT
           zlkp10(jl,jk))
      CALL access_lookup_table(&
           ! IN
           dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
           dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),5,    &
           ! OUT
           zlkp5(jl,jk))
    ENDDO
  ENDDO

! ice present logical
  ll_ice(1:kproma,:) = pxi(1:kproma,:) .gt. cqtmin

  prim(1:kproma,:) = MERGE(0.5*zlkp11(1:kproma,:), 0._dp, ll_ice(1:kproma,:))
  prieff(1:kproma,:) = MERGE(zlkp6(1:kproma,:), 0._dp, ll_ice(1:kproma,:))
  priv(:,:) = effective_2_volmean_radius_param_Schuman_2011( &
                         kbdim, kproma, klev, prieff(:,:))
  priv(1:kproma,:) = MERGE(priv(1:kproma,:), 0._dp, ll_ice(1:kproma,:))
  pvtin(1:kproma,:) = MERGE(prhofaci(1:kproma,:)*zlkp1(1:kproma,:)*vt_rdc, 0._dp, ll_ice(1:kproma,:))
  pvtim(1:kproma,:) = MERGE(prhofaci(1:kproma,:)*zlkp2(1:kproma,:)*vt_rdc, 0._dp, ll_ice(1:kproma,:))
  prhoice(1:kproma,:) = MERGE(zlkp12(1:kproma,:), 0._dp, ll_ice(1:kproma,:))

  ! normalize Morrison lookup values to be the same as my lookup.
  plkp_depx1(1:kproma,:) = MERGE(zlkp5(1:kproma,:)/0.65/picnc(1:kproma,:), 0._dp, picnc(1:kproma,:)>cqtmin)
  plkp_depx2(1:kproma,:) = MERGE(zlkp10(1:kproma,:)/0.44/picnc(1:kproma,:), 0._dp, picnc(1:kproma,:)>cqtmin)

END SUBROUTINE secondary_ice_properties_p3

SUBROUTINE secondary_ice_properties_2m_2d(&
              !--IN
              kbdim, klev, kproma, &
              prho, pxi, picnc, &
              paaa, &
              !--OUT
              priv, prieff, pvtim, pvtin)

  USE mo_cloud_micro_2m, ONLY : eff_ice_crystal_radius

  INTEGER, INTENT(IN) :: kbdim, klev, kproma
  REAL(dp), DIMENSION(kbdim,klev), INTENT(IN) :: prho, pxi, picnc
  REAL(dp), DIMENSION(kbdim,klev), INTENT(IN) :: paaa
  
  REAL(dp), DIMENSION(kbdim,klev), INTENT(OUT) :: priv, prieff, pvtim, pvtin

  INTEGER :: jk
  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1, ztmp2, ztmp3, zicnc_m3
  REAL(dp), DIMENSION(kbdim,klev) :: zmmean, zalfased, zbetased
  LOGICAL, DIMENSION(kbdim,klev)  :: ll1, ll2

!------------------------------------------------------
! effective & volume mean ice crystal radius

  ! set up number per m3
  zicnc_m3(:,:) = 1._dp
  zicnc_m3(1:kproma,:) = prho(1:kproma,:)*picnc(1:kproma,:)
  zicnc_m3(:,:) = MAX(1._dp, zicnc_m3(:,:))

  ztmp1(:,:) = 0._dp
  ztmp1(1:kproma,:) = MAX(0._dp, 1000._dp*pxi(1:kproma,:)*prho(1:kproma,:))
  DO jk=1,klev
     ! rieff in (m)
     ztmp2(:,jk) = 1e-6_dp*eff_ice_crystal_radius(kbdim, kproma, ztmp1(:,jk), zicnc_m3(:,jk))
     ! riv in (m)
     ztmp3(:,jk) = effective_2_volmean_radius_param_Schuman_2011( &
                                              kbdim, kproma, ztmp2(:,jk))
  END DO !jk
  
  prieff(1:kproma,:) = MERGE(ztmp2(1:kproma,:), 0._dp, pxi(1:kproma,:) > cqtmin)
  priv(1:kproma,:) = MERGE(ztmp3(1:kproma,:), 0._dp, pxi(1:kproma,:) > cqtmin)

!------------------------------------------------------
! number and mass-weighted fallvelocities
  
  !mean ice crystal mass
  zmmean(1:kproma,:) = pxi(1:kproma,:)/MAX(picnc(1:kproma,:), 1._dp)
  zmmean(1:kproma,:) = MAX(zmmean(1:kproma,:), mi)
  zmmean(1:kproma,:) = MERGE(zmmean(1:kproma,:), 0._dp, picnc(1:kproma,:) > 1._dp)

  !get regime for fall velocity parametrization
  ll1(1:kproma,:) = (zmmean(1:kproma,:) < ri_vol_mean_1 )
  ll2(1:kproma,:) = (.NOT. ll1(1:kproma,:)) .AND. (zmmean(1:kproma,:) < ri_vol_mean_2 )

  zalfased(1:kproma,:) = MERGE(alfased_1, alfased_2, ll1(1:kproma,:))
  zalfased(1:kproma,:) = MERGE(alfased_3, zalfased(1:kproma,:), ll2(1:kproma,:))

  zbetased(1:kproma,:) = MERGE(betased_1, betased_2, ll1(1:kproma,:))
  zbetased(1:kproma,:) = MERGE(betased_3, zbetased(1:kproma,:), ll2(1:kproma,:))

  pvtim(1:kproma,:)    = fall*zalfased(1:kproma,:) & 
                         *(zmmean(1:kproma,:)**zbetased(1:kproma,:))*paaa(1:kproma,:)

  !limit fall velocity to 1 cm/s - 2 m/s:
  pvtim(1:kproma,:) = MAX(0.001_dp,pvtim(1:kproma,:))
  pvtim(1:kproma,:) = MIN(2._dp,pvtim(1:kproma,:))

  !assume number- and mass-weighted fallvelocities are equal
  pvtin(1:kproma,:) = pvtim(1:kproma,:)

END SUBROUTINE secondary_ice_properties_2m_2d

SUBROUTINE secondary_ice_properties_2m_1d(&
              !--IN
              kbdim, kproma, &
              prho, pxi, picnc, &
              paaa, &
              !--OUT
              priv, prieff, pvtim, pvtin)

  USE mo_cloud_micro_2m, ONLY : eff_ice_crystal_radius

  INTEGER, INTENT(IN) :: kbdim, kproma
  REAL(dp), DIMENSION(kbdim), INTENT(IN) :: prho, pxi, picnc
  REAL(dp), DIMENSION(kbdim), INTENT(IN) :: paaa
  
  REAL(dp), DIMENSION(kbdim), INTENT(OUT) :: priv, prieff, pvtim, pvtin

  INTEGER :: jk
  REAL(dp), DIMENSION(kbdim) :: ztmp1, ztmp2, ztmp3, zicnc_m3
  REAL(dp), DIMENSION(kbdim) :: zmmean, zalfased, zbetased
  LOGICAL, DIMENSION(kbdim)  :: ll1, ll2

!------------------------------------------------------
! effective & volume mean ice crystal radius

  ! set up number per m3
  zicnc_m3(:) = 1._dp
  zicnc_m3(1:kproma) = prho(1:kproma)*picnc(1:kproma)
  zicnc_m3(:) = MAX(1._dp, zicnc_m3(:))

  ztmp1(:) = 0._dp
  ztmp1(1:kproma) = MAX(0._dp, 1000._dp*pxi(1:kproma)*prho(1:kproma))
  ! rieff in (m)
  ztmp2(:) = 1e-6_dp*eff_ice_crystal_radius(kbdim, kproma, ztmp1(:), zicnc_m3(:))
  ! riv in (m)
  ztmp3(:) = effective_2_volmean_radius_param_Schuman_2011( &
                                              kbdim, kproma, ztmp2(:))
 
  ! prieff(1:kproma) = MIN(MAX(ztmp2(1:kproma), 1e-6*ceffmin), 1e-6*ceffmax)
  prieff(1:kproma) = MERGE(ztmp2(1:kproma), 0._dp, pxi(1:kproma) > cqtmin)
  priv(1:kproma) = MERGE(ztmp3(1:kproma), 0._dp, pxi(1:kproma) > cqtmin)

!------------------------------------------------------
! number and mass-weighted fallvelocities
  
  !mean ice crystal mass
  zmmean(1:kproma) = pxi(1:kproma)/MAX(picnc(1:kproma), 1._dp)
  zmmean(1:kproma) = MAX(zmmean(1:kproma), mi)
  zmmean(1:kproma) = MERGE(zmmean(1:kproma), 0._dp, picnc(1:kproma) > 1._dp)

  !get regime for fall velocity parametrization
  ll1(1:kproma) = (zmmean(1:kproma) < ri_vol_mean_1 )
  ll2(1:kproma) = (.NOT. ll1(1:kproma)) .AND. (zmmean(1:kproma) < ri_vol_mean_2 )

  zalfased(1:kproma) = MERGE(alfased_1, alfased_2, ll1(1:kproma))
  zalfased(1:kproma) = MERGE(alfased_3, zalfased(1:kproma), ll2(1:kproma))

  zbetased(1:kproma) = MERGE(betased_1, betased_2, ll1(1:kproma))
  zbetased(1:kproma) = MERGE(betased_3, zbetased(1:kproma), ll2(1:kproma))

  pvtim(1:kproma)    = fall*zalfased(1:kproma) & 
                         *(zmmean(1:kproma)**zbetased(1:kproma))*paaa(1:kproma)

  !limit fall velocity to 1 cm/s - 2 m/s:
  pvtim(1:kproma) = MAX(0.001_dp,pvtim(1:kproma))
  pvtim(1:kproma) = MIN(2._dp,pvtim(1:kproma))

  !assume number- and mass-weighted fallvelocities are equal
  pvtin(1:kproma) = pvtim(1:kproma)

END SUBROUTINE secondary_ice_properties_2m_1d

SUBROUTINE secondary_prognostics_1d(&
              !--IN
              kbdim, kproma, pxi, &
              pqihet, pqiliq, pqioliq, pqsrc, pqprc, &
              picnc, pnihet, pnihom, pninuc, pnidet, &
              !--OUT
              phetfrac, pliqfrac, poliqfrac, psosifrac, &
              pnihetfrac, pnihomfrac, pninucfrac, pnidetfrac)

 INTEGER, INTENT(IN) :: kbdim, kproma
 REAL(dp), INTENT(IN) :: pxi(kbdim)
 REAL(dp), INTENT(IN) :: pqihet(kbdim), pqiliq(kbdim), pqioliq(kbdim), &
                         pqsrc(kbdim), pqprc(kbdim), &
                         picnc(kbdim), pnihet(kbdim), pnihom(kbdim), pninuc(kbdim), pnidet(kbdim)

 REAL(dp), INTENT(OUT) :: phetfrac(kbdim), pliqfrac(kbdim), poliqfrac(kbdim), &
                          psosifrac(kbdim), &
                          pnihetfrac(kbdim), pnihomfrac(kbdim), pninucfrac(kbdim), pnidetfrac(kbdim)

 LOGICAL, DIMENSION(kbdim) :: ll1
 REAL(dp), DIMENSION(kbdim) :: ztmp1
 
! heterogeneously formed fraction
 ll1(1:kproma) = pxi(1:kproma) > eps
 phetfrac(1:kproma) = MERGE(pqihet(1:kproma)/pxi(1:kproma), 0._dp, ll1(1:kproma))
 phetfrac(1:kproma) = MAX(MIN(1._dp, phetfrac(1:kproma)), 0._dp)

! fraction of ice formed by frozen liquid
 ll1(1:kproma) = pxi(1:kproma) > eps
 pliqfrac(1:kproma) = MERGE(pqiliq(1:kproma)/pxi(1:kproma), 0._dp, ll1(1:kproma))
 pliqfrac(1:kproma) = MAX(MIN(1._dp, pliqfrac(1:kproma)), 0._dp)

! fraction of ice that formed through frozen liquid water
 ll1(1:kproma) = pxi(1:kproma) > eps
 poliqfrac(1:kproma) = MERGE(pqioliq(1:kproma)/pxi(1:kproma), 0._dp, ll1(1:kproma))
 poliqfrac(1:kproma) = MAX(MIN(1._dp, poliqfrac(1:kproma)), 0._dp)

! initialization by freezing can be very small, therefore we assume
! that the origin fractions can only be 0 or 1 for ice contents smaller than 1e-10 to avoid
! rounding errors of 2 very small numbers
 ll1(1:kproma) = pxi(1:kproma) < 1e-10_dp
 ztmp1(1:kproma) = MERGE(0._dp, 1._dp, phetfrac(1:kproma) < 0.5_dp)
 phetfrac(1:kproma) = MERGE(ztmp1(1:kproma), phetfrac(1:kproma), ll1(1:kproma))
 ztmp1(1:kproma) = MERGE(0._dp, 1._dp, poliqfrac(1:kproma) < 0.5_dp)
 poliqfrac(1:kproma) = MERGE(ztmp1(1:kproma), poliqfrac(1:kproma), ll1(1:kproma))

! fraction of source and sink terms
 ll1(1:kproma) = pqsrc(1:kproma) > eps .AND. pxi(1:kproma) > eps
 psosifrac(1:kproma) = MERGE(pqprc(1:kproma)/pqsrc(1:kproma), 0._dp, ll1(1:kproma))
 psosifrac(1:kproma) = MAX(0._dp, MIN(1._dp, psosifrac(1:kproma)))

! heterogeneously formed number fraction
 ll1(1:kproma) = picnc(1:kproma) > eps .AND. pxi(1:kproma) > eps
 pnihetfrac(1:kproma) = MERGE(pnihet(1:kproma)/picnc(1:kproma), 0._dp, ll1(1:kproma))
 pnihetfrac(1:kproma) = MAX(MIN(1._dp, pnihetfrac(1:kproma)), 0._dp)

! homogeneously formed number fraction
 ll1(1:kproma) = picnc(1:kproma) > eps .AND. pxi(1:kproma) > eps
 pnihomfrac(1:kproma) = MERGE(pnihom(1:kproma)/picnc(1:kproma), 0._dp, ll1(1:kproma))
 pnihomfrac(1:kproma) = MAX(MIN(1._dp, pnihomfrac(1:kproma)), 0._dp)

! nucleated number fraction
 ll1(1:kproma) = picnc(1:kproma) > eps .AND. pxi(1:kproma) > eps
 pninucfrac(1:kproma) = MERGE(pninuc(1:kproma)/picnc(1:kproma), 0._dp, ll1(1:kproma))
 pninucfrac(1:kproma) = MAX(MIN(1._dp, pninucfrac(1:kproma)), 0._dp)

! heterogeneously formed number fraction
 ll1(1:kproma) = picnc(1:kproma) > eps .AND. pxi(1:kproma) > eps
 pnidetfrac(1:kproma) = MERGE(pnidet(1:kproma)/picnc(1:kproma), 0._dp, ll1(1:kproma))
 pnidetfrac(1:kproma) = MAX(MIN(1._dp, pnidetfrac(1:kproma)), 0._dp)

END SUBROUTINE secondary_prognostics_1d

SUBROUTINE secondary_prognostics_2d(&
              !--IN
              kbdim, klev, kproma, pxi, &
              pqihet, pqiliq, pqioliq, pqsrc, pqprc, &
              picnc, pnihet, pnihom, pninuc, pnidet, &
              !--OUT
              phetfrac, pliqfrac, poliqfrac, psosifrac, &
              pnihetfrac, pnihomfrac, pninucfrac, pnidetfrac)

 INTEGER, INTENT(IN) :: kbdim, klev, kproma
 REAL(dp), INTENT(IN) :: pxi(kbdim,klev)
 REAL(dp), INTENT(IN) :: pqihet(kbdim,klev), pqiliq(kbdim,klev), pqioliq(kbdim,klev), &
                         pqsrc(kbdim,klev), pqprc(kbdim,klev), &
                         picnc(kbdim,klev), pnihet(kbdim,klev), pnihom(kbdim,klev), pninuc(kbdim,klev), &
                         pnidet(kbdim,klev)

 REAL(dp), INTENT(OUT) :: phetfrac(kbdim,klev), pliqfrac(kbdim,klev), poliqfrac(kbdim,klev), &
                          psosifrac(kbdim,klev), &
                          pnihetfrac(kbdim,klev), pnihomfrac(kbdim,klev), pninucfrac(kbdim,klev), &
                          pnidetfrac(kbdim,klev)

 LOGICAL, DIMENSION(kbdim,klev) :: ll1
 REAL(dp), DIMENSION(kbdim,klev) :: ztmp1
 
! heterogeneously formed fraction
 ll1(1:kproma,:) = pxi(1:kproma,:) > eps
 phetfrac(1:kproma,:) = MERGE(pqihet(1:kproma,:)/pxi(1:kproma,:), 0._dp, ll1(1:kproma,:))
 phetfrac(1:kproma,:) = MAX(MIN(1._dp, phetfrac(1:kproma,:)), 0._dp)

! fraction of ice formed by frozen liquid
 ll1(1:kproma,:) = pxi(1:kproma,:) > eps
 pliqfrac(1:kproma,:) = MERGE(pqiliq(1:kproma,:)/pxi(1:kproma,:), 0._dp, ll1(1:kproma,:))
 pliqfrac(1:kproma,:) = MAX(MIN(1._dp, pliqfrac(1:kproma,:)), 0._dp)

! fraction of ice that formed through frozen liquid water
 ll1(1:kproma,:) = pxi(1:kproma,:) > eps
 poliqfrac(1:kproma,:) = MERGE(pqioliq(1:kproma,:)/pxi(1:kproma,:), 0._dp, ll1(1:kproma,:))
 poliqfrac(1:kproma,:) = MAX(MIN(1._dp, poliqfrac(1:kproma,:)), 0._dp)

! initialization by freezing can be very small, therefore we assume
! that the origin fractions can only be 0 or 1 for ice contents smaller than 1e-10 to avoid
! rounding errors of 2 very small numbers
 ll1(1:kproma,:) = pxi(1:kproma,:) < 1e-10_dp
 ztmp1(1:kproma,:) = MERGE(0._dp, 1._dp, phetfrac(1:kproma,:) < 0.5_dp)
 phetfrac(1:kproma,:) = MERGE(ztmp1(1:kproma,:), phetfrac(1:kproma,:), ll1(1:kproma,:))
 ztmp1(1:kproma,:) = MERGE(0._dp, 1._dp, poliqfrac(1:kproma,:) < 0.5_dp)
 poliqfrac(1:kproma,:) = MERGE(ztmp1(1:kproma,:), poliqfrac(1:kproma,:), ll1(1:kproma,:))

! fraction of source and sink terms
 ll1(1:kproma,:) = pqsrc(1:kproma,:) > eps
 psosifrac(1:kproma,:) = MERGE(pqprc(1:kproma,:)/pqsrc(1:kproma,:), 0._dp, ll1(1:kproma,:))
 psosifrac(1:kproma,:) = MAX(0._dp, MIN(1._dp, psosifrac(1:kproma,:)))

! fraction of source and sink terms
 ll1(1:kproma,:) = pqsrc(1:kproma,:) > eps
 psosifrac(1:kproma,:) = MERGE(pqprc(1:kproma,:)/pqsrc(1:kproma,:), 0._dp, ll1(1:kproma,:))
 psosifrac(1:kproma,:) = MAX(0._dp, MIN(1._dp, psosifrac(1:kproma,:)))

! heterogeneously formed number fraction
 ll1(1:kproma,:) = picnc(1:kproma,:) > eps
 pnihetfrac(1:kproma,:) = MERGE(pnihet(1:kproma,:)/picnc(1:kproma,:), 0._dp, ll1(1:kproma,:))
 pnihetfrac(1:kproma,:) = MAX(MIN(1._dp, pnihetfrac(1:kproma,:)), 0._dp)

! homogeneously formed number fraction
 ll1(1:kproma,:) = picnc(1:kproma,:) > eps
 pnihomfrac(1:kproma,:) = MERGE(pnihom(1:kproma,:)/picnc(1:kproma,:), 0._dp, ll1(1:kproma,:))
 pnihomfrac(1:kproma,:) = MAX(MIN(1._dp, pnihomfrac(1:kproma,:)), 0._dp)

! nucleated number fraction
 ll1(1:kproma,:) = picnc(1:kproma,:) > eps
 pninucfrac(1:kproma,:) = MERGE(pninuc(1:kproma,:)/picnc(1:kproma,:), 0._dp, ll1(1:kproma,:))
 pninucfrac(1:kproma,:) = MAX(MIN(1._dp, pninucfrac(1:kproma,:)), 0._dp)

! heterogeneously formed number fraction
 ll1(1:kproma,:) = picnc(1:kproma,:) > eps
 pnidetfrac(1:kproma,:) = MERGE(pnidet(1:kproma,:)/picnc(1:kproma,:), 0._dp, ll1(1:kproma,:))
 pnidetfrac(1:kproma,:) = MAX(MIN(1._dp, pnidetfrac(1:kproma,:)), 0._dp)

END SUBROUTINE secondary_prognostics_2d

SUBROUTINE deposition(&
              !--IN
              kproma, kbdim,  &
              pxi, pt, psupi,      & 
              pastbsti, paclc, prho, paaa, pviscos, &
              pdv, pvtim, prim, picnc,                    &
              !--OUT
              pdep)

  INTEGER, INTENT(IN) :: kproma, kbdim
  REAL(dp), INTENT(IN), DIMENSION(kbdim) :: pxi, pt, psupi, pastbsti, &
                                            paclc, prho, paaa, pviscos, pdv, pvtim, prim, picnc
  REAL(dp), INTENT(OUT), DIMENSION(kbdim) :: pdep

  LOGICAL, DIMENSION(kbdim) :: ll_ice, ll1, ll2
  REAL(dp), DIMENSION(kbdim) :: ztmp1
  REAL(dp), DIMENSION(kbdim) :: zb2, zgtp, zvth, zfuchs, zfre, zre !< hydrodynamics

  ll_ice(1:kproma) = (pxi(1:kproma) > 0._dp) .AND. &
                     (paclc(1:kproma) > clc_min)

  zgtp(1:kproma)    = 1._dp/(prho(1:kproma)*pastbsti(1:kproma))
  zvth(1:kproma)    = SQRT( 8._dp*kb*pt(1:kproma) / (pi*xmw) )
  zb2(1:kproma)     = 0.25_dp * alpha * zvth(1:kproma)   / pdv(1:kproma)
  zfuchs(1:kproma)  = 1._dp/(1._dp+zb2(1:kproma)*prim(1:kproma))
  zre(1:kproma)     = 2._dp*prho(1:kproma)*prim(1:kproma)*pvtim(1:kproma)/pviscos(1:kproma)
  zfre(1:kproma)    = 1._dp + 0.229_dp*SQRT(zre(1:kproma))
  zfre(1:kproma)    = MERGE(zfre(1:kproma), 1._dp, ll_ice(1:kproma))

  pdep(1:kproma)  = 4._dp*pi*prim(1:kproma)*psupi(1:kproma)*picnc(1:kproma) &
                    *zfre(1:kproma)*zgtp(1:kproma)*zfuchs(1:kproma)*alpha

END SUBROUTINE deposition

SUBROUTINE deposition_l_1d(&
              !--IN
              kproma, kbdim,  &
              pxi, pt, psupi, pesi, papm1,   & 
              pastbsti, plsdcp, pvt, &
              prim, picnc, pfr, prho, pviscos, &
              plkp_depx1, plkp_depx2, prhofaci, psc, &
              !--OUT
              pdep)

  INTEGER, INTENT(IN) :: kproma, kbdim
  REAL(dp), INTENT(IN), DIMENSION(kbdim) :: pxi, pt, psupi, pastbsti, plsdcp, &
                                            prim, picnc, papm1, pesi
  REAL(dp), INTENT(IN), DIMENSION(kbdim) :: pfr, pvt, prho, pviscos
  REAL(dp), INTENT(IN), DIMENSION(kbdim) :: plkp_depx1, plkp_depx2, prhofaci, psc
  REAL(dp), INTENT(OUT), DIMENSION(kbdim) :: pdep

  REAL(dp), DIMENSION(kbdim) :: zfk, zfd, zcap, zfvs, zfvg, zfv, zreys

  !F_k
  zfk(1:kproma) = (plsdcp(1:kproma)/(rv*pt(1:kproma))-1)*plsdcp(1:kproma) &
                    /(con0_h*pt(1:kproma))
  !F_d
  zfd(1:kproma) = 2.11*1e-5*(pt(1:kproma)/tmelt)**1.94 & ! thermal conductivity of
                    *(101325/papm1(1:kproma))              ! water vapor in air
  zfd(1:kproma) = rv*pt(1:kproma)/(zfd(1:kproma)*pesi(1:kproma))

  !Sqrt(Reynolds number) of ice in air
  !plkp_depx2 = Integrate[ (D*v(D))^1/2*C(D)*m(D)*N(D) ] / Integrate[ m(D)*N(D) ]
  zreys(1:kproma) = (prhofaci(1:kproma)*prho(1:kproma)/pviscos(1:kproma))**0.5*plkp_depx2(1:kproma)
  zreys(1:kproma) = MAX(zreys(1:kproma), 0._dp)

  !Ventilation coefficient for snow aggregates (s) and graupel (g) (Straka 6.18)
  !plkp_depx1 = Integrate[ C(D)*m(D)*N(D) ] / Integrate[ m(D)*N(D) ]
  zfvs(1:kproma) = 0.65_dp*plkp_depx1(1:kproma) + 0.44_dp*psc(1:kproma)**(1._dp/3._dp)*zreys(1:kproma)
  zfvg(1:kproma) = 0.78_dp*plkp_depx1(1:kproma) + 0.308_dp*psc(1:kproma)**(1._dp/3._dp)*zreys(1:kproma)
  !Get total ventilation coefficient by weighting by rimed fraction
  ! zfv(1:kproma) = pfr(1:kproma)*zfvg(1:kproma) + (1._dp-pfr(1:kproma))*zfvs(1:kproma)
  zfv(1:kproma) = zfvs(1:kproma)

  !Sublimation rate: 0.5 is alpha_m == accomodation coefficient 
  !            (prob. of molecule attaching to lattice)
  !note that the capacitance is size dependant and thus included in lookup values
  pdep(1:kproma) = picnc(1:kproma)*0.5*4*pi*psupi(1:kproma)*zfv(1:kproma) &
                    /(zfk(1:kproma)+zfd(1:kproma))

END SUBROUTINE deposition_l_1d

SUBROUTINE deposition_l_2d(&
              !--IN
              kproma, kbdim, klev,  &
              pxi, pt, psupi, pesi, papm1,   & 
              pastbsti, plsdcp, pvt, &
              prim, picnc, pfr, prho, pviscos, &
              plkp_depx1, plkp_depx2, prhofaci, psc, &
              !--OUT
              pdep)

  INTEGER, INTENT(IN) :: kproma, kbdim, klev
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: pxi, pt, psupi, pastbsti, plsdcp, &
                                            prim, picnc, papm1, pesi
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: pfr, pvt, prho, pviscos
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: plkp_depx1, plkp_depx2, prhofaci, psc
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev) :: pdep

  REAL(dp), DIMENSION(kbdim,klev) :: zfk, zfd, zcap, zfvs, zfvg, zfv, zreys

  !F_k
  zfk(1:kproma,:) = (plsdcp(1:kproma,:)/(rv*pt(1:kproma,:))-1)*plsdcp(1:kproma,:) &
                    /(con0_h*pt(1:kproma,:))
  !F_d
  zfd(1:kproma,:) = 2.11*1e-5*(pt(1:kproma,:)/tmelt)**1.94 & ! thermal conductivity of
                    *(101325/papm1(1:kproma,:))              ! water vapor in air
  zfd(1:kproma,:) = rv*pt(1:kproma,:)/(zfd(1:kproma,:)*pesi(1:kproma,:))

  !Sqrt(Reynolds number) of ice in air
  !plkp_depx2 = Integrate[ (D*v(D))^1/2*C(D)*m(D)*N(D) ] / Integrate[ m(D)*N(D) ]
  zreys(1:kproma,:) = (prhofaci(1:kproma,:)*prho(1:kproma,:)/pviscos(1:kproma,:))**0.5*plkp_depx2(1:kproma,:)
  zreys(1:kproma,:) = MAX(zreys(1:kproma,:), 0._dp)

  !Ventilation coefficient for snow aggregates (s) and graupel (g) (Straka 6.18)
  !plkp_depx1 = Integrate[ C(D)*m(D)*N(D) ] / Integrate[ m(D)*N(D) ]
  zfvs(1:kproma,:) = 0.65_dp*plkp_depx1(1:kproma,:) + 0.44_dp*psc(1:kproma,:)**(1._dp/3._dp)*zreys(1:kproma,:)
  zfvg(1:kproma,:) = 0.78_dp*plkp_depx1(1:kproma,:) + 0.308_dp*psc(1:kproma,:)**(1._dp/3._dp)*zreys(1:kproma,:)
  !Get total ventilation coefficient by weighting by rimed fraction
  ! zfv(1:kproma,:) = pfr(1:kproma,:)*zfvg(1:kproma,:) + (1._dp-pfr(1:kproma,:))*zfvs(1:kproma,:)
  zfv(1:kproma,:) = zfvs(1:kproma,:)

  !Sublimation rate: 0.5 is alpha_m == accomodation coefficient 
  !            (prob. of molecule attaching to lattice)
  !note that the capacitance is size dependant and thus included in lookup values
  pdep(1:kproma,:) = picnc(1:kproma,:)*0.5*4*pi*psupi(1:kproma,:)*zfv(1:kproma,:) &
                    /(zfk(1:kproma,:)+zfd(1:kproma,:))

END SUBROUTINE deposition_l_2d

SUBROUTINE dep_from_cirrus(&
              !--IN
              kproma, kbdim, klev, ktdia, ptmst_rcp,     &
              papm1, pxi, picnc_new, pq, pt, &
              plsdcp, pastbsti, psupi, pesi, pqsi, & 
              pri, paclc, prho, paaa, pviscos, &
              pdv, pvtim, prim, picnc,                    &
              plkp_depx1, plkp_depx2, prhofaci, psc, &
              !--OUT
              pdep_ci)

  INTEGER, INTENT(IN) :: kproma, kbdim, klev, ktdia
  REAL(dp), INTENT(IN) :: ptmst_rcp
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: pxi, picnc_new, pq, pt, psupi, pqsi, pri, &
                                                 paclc, prho, paaa, pviscos, pdv, pvtim, prim, picnc,&
                                                 plkp_depx1, plkp_depx2, prhofaci, psc, &
                                                 papm1, pesi, plsdcp, pastbsti
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev) :: pdep_ci

  LOGICAL, DIMENSION(kbdim,klev) :: ll_ice, ll1, ll2
  REAL(dp), DIMENSION(kbdim,klev) :: zri, zrid, zicnc_tot
  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1, ztmp2
  REAL(dp) :: zb2, zgtp, zvth, zfuchs, zfre, zre !< for cirrus calculations
  REAL(dp) :: zmmean(kbdim,klev)      !< Volume mean ice crystal radius [m]
  REAL(dp) :: zalfased(kbdim,klev)    !< Parameter needed for the ice crystal fall velocity
  REAL(dp) :: zbetased(kbdim,klev)    !< Parameter needed for the ice crystal fall velocity
  REAL(dp) :: zxifallmc(kbdim,klev)   !< Fall velocity of ice crystal mass
  INTEGER :: jk, jl

  zicnc_tot(1:kproma,:) = picnc_new(1:kproma,:) + picnc(1:kproma,:)

  ll_ice(1:kproma,:) = (pxi(1:kproma,:) > 0._dp) .AND. &
                       (paclc(1:kproma,:) > clc_min)


  ! ice radius proxy

  ztmp1(1:kproma,:) = pt(1:kproma,:)-tmelt
  ztmp1(1:kproma,:) = MIN(0._dp,ztmp1(1:kproma,:))

  ztmp2(1:kproma,:) = 0.015_dp*ztmp1(1:kproma,:)
  ztmp2(1:kproma,:) = EXP(ztmp2(1:kproma,:))
  ztmp2(1:kproma,:) = 23.2_dp*ztmp2(1:kproma,:)
  ztmp2(1:kproma,:) = MAX(ztmp2(1:kproma,:),1.0_dp)
  zrid(:,:) = effective_2_volmean_radius_param_Schuman_2011( &
                         kbdim, kproma, klev, ztmp2(:,:))

  ! handle the case where zri is too small

  ll1(1:kproma,:) = (pri(1:kproma,:) >= epsec) 

  zri(1:kproma,:) = MERGE(pri(1:kproma,:), zrid(1:kproma,:), ll1(1:kproma,:))
  zri(1:kproma,:) = MAX(zri(1:kproma,:), 1.e-6_dp)
  ! take radius of pre-exising if they are bigger
  zri(1:kproma,:) = MAX(zri(1:kproma,:), prim(1:kproma,:))

  ! mean ice mass

  ztmp1(1:kproma,:) = pxi(1:kproma,:) & 
                      / (MAX(picnc(1:kproma,:), icemin) * MAX(paclc(1:kproma,:),clc_min))
  ztmp1(1:kproma,:) = MAX(ztmp1(1:kproma,:), mi)

  zmmean(1:kproma,:)   = MERGE(ztmp1(1:kproma,:), mi, ll_ice(1:kproma,:))

  ! fall speed

  ll1(1:kproma,:) = (zmmean(1:kproma,:) < ri_vol_mean_1 )
  ll2(1:kproma,:) = ( .NOT. ll1(1:kproma,:) ) .AND. &
                       (zmmean(1:kproma,:) < ri_vol_mean_2 )

  zalfased(1:kproma,:) = MERGE(alfased_1, alfased_2, ll1(1:kproma,:))
  zalfased(1:kproma,:) = MERGE(alfased_3, zalfased(1:kproma,:), ll2(1:kproma,:))

  zbetased(1:kproma,:) = MERGE(betased_1, betased_2, ll1(1:kproma,:))
  zbetased(1:kproma,:) = MERGE(betased_3, zbetased(1:kproma,:), ll2(1:kproma,:))

  zxifallmc(1:kproma,:) = fall*zalfased(1:kproma,:) &
                        * (zmmean(1:kproma,:)**zbetased(1:kproma,:))*paaa(1:kproma,:)

  zxifallmc(1:kproma,:) = MAX(zxifallmc(1:kproma,:), pvtim(1:kproma,:))
 
  ! rime fraction (all snow at cirrus temperatures)
  ztmp1(1:kproma,:) = 0._dp

  CALL deposition_l(&
          !--IN
          kproma, kbdim, klev,  &
          pxi, pt, psupi, pesi, papm1,   & 
          pastbsti, plsdcp, zxifallmc, &
          zri, zicnc_tot, ztmp1, prho, pviscos, &
          plkp_depx1, plkp_depx2, prhofaci, psc, &
          !--OUT
          pdep_ci)

  pdep_ci(1:kproma,:) = MIN(pdep_ci(1:kproma,:), (pq(1:kproma,:)-pqsi(1:kproma,:))*ptmst_rcp)

END SUBROUTINE dep_from_cirrus

SUBROUTINE vertical_velocity(&
              !--IN
              kproma, kbdim, klev, ktdia, krow, &
              ptkem1, paphm1, papm1, pgeo, pt, pq, &
              prho_rcp, pvervel, pqsi, &
              !--OUT
              pvervx)

  USE mo_memory_g3b,         ONLY: orostd,oromea,orogam,orothe, &
                                   oropic,oroval,orosig
  USE mo_memory_g2a,         ONLY: vm1, um1
  USE mo_cloud_utils,        ONLY: fact_tke
  USE mo_physical_constants, ONLY: rgrav
#ifdef HAMMOZ
  USE mo_orocirrus,          ONLY: orocirrus_w, orocirrus_cc
#endif
  
  INTEGER, INTENT(IN)                          :: kproma, kbdim, klev, ktdia, krow
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev)  :: ptkem1, paphm1, papm1, pgeo, pt, pq, &
                                                  prho_rcp, pvervel, pqsi
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev) :: pvervx

  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1

  ztmp1(1:kproma,:)    = 100._dp*fact_tke*SQRT(ptkem1(1:kproma,:))
  pvervx(1:kproma,:)   = -100._dp*rgrav*pvervel(1:kproma,:)*prho_rcp(1:kproma,:) + ztmp1(1:kproma,:)

#ifdef HAMMOZ
!  Update updraft velocity to take into account orography, if relevant:

   IF (lorocirrus) &
      CALL orocirrus_w( &
              !-- IN
              kproma, kbdim, klev, ktdia, &
              paphm1(:,:), &
              papm1(:,:), pgeo(:,:), pt(:,:), pq(:,:), um1(:,:,krow), &
              vm1(:,:,krow), orostd(:,krow), &
              orosig(:,krow), oromea(:,krow), &
              orogam(:,krow), orothe(:,krow), &
              oropic(:,krow), oroval(:,krow), &
              rgrav, prho_rcp(:,:), pvervel(:,:), pqsi(:,:), &
              !-- INOUT
              pvervx(:,:))
#endif

END SUBROUTINE vertical_velocity

SUBROUTINE aerosol_activation(&
              !--IN
              kbdim, kproma, klev, ktdia, &
              ibas, icl_minusbas, paclc, &
              pt, pxl, pcdnc, pcdncact, &
              !--OUT
              pncnuc)
  
  INTEGER, INTENT(IN)                         :: kbdim, kproma, klev, ktdia
  INTEGER, INTENT(IN), DIMENSION(kbdim,klev)  :: ibas, icl_minusbas
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: paclc
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: pt, pxl, pcdnc, pcdncact
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev):: pncnuc

  INTEGER :: jk, jl, jkk
  LOGICAL, DIMENSION(kbdim,klev)  :: ll1, ll2
  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1
  REAL(dp), DIMENSION(kbdim,klev) :: zncnuc_bas, zcdnc_bas, zcdnc

  zcdnc_bas(1:kproma,:) = 0._dp
  zcdnc(1:kproma,:) = 0._dp
  zncnuc_bas(1:kproma,:) = 0._dp
  pncnuc(1:kproma,:) = 0._dp

!--- Convert the aerosol activation into the number of newly formed cloud droplets
  ll1(1:kproma,:) = ibas(1:kproma,:) > 0 & 
                    .AND. pt(1:kproma,:) > cthomi

  !SF: first computes newly formed cloud droplets at cloud bases:
  ztmp1(1:kproma,:) = pcdncact(1:kproma,:) - pcdnc(1:kproma,:)
  pncnuc(1:kproma,:) = MAX(0._dp, ztmp1(1:kproma,:))
  pncnuc(1:kproma,:) = MERGE(pncnuc(1:kproma,:), 0._dp, ll1(1:kproma,:))

  !RD: array with updated values at cloud bases
  zcdnc(1:kproma,:) = pcdnc(1:kproma,:) + pncnuc(1:kproma,:)
 
  !SF: then computes newly formed cloud droplets above cloud base
  !    by assuming that the number of nucleated cloud droplets is constant above cloud base
  !    (adiabatic parcel theory)

  DO jk=ktdia,klev
     DO jl=1,kproma
        jkk = MAX(jk,icl_minusbas(jl,jk))  !sets the level index to either relevant cloud base or itself
        zncnuc_bas(jl,jk) = pncnuc(jl,jkk) !holds in each cloud level the number of
                                           !newly formed cloud droplets at the base

        zcdnc_bas(jl,jk) = zcdnc(jl,jkk)   !propagate cloud base values throughout cloud
     ENDDO !end jl
  ENDDO !end jk

  ll2(1:kproma,:) = (icl_minusbas(1:kproma,:) > 0)     .AND. &  !within a cloud
                    (zncnuc_bas(1:kproma,:)   > 0._dp)          !positive nucleation at base

  ! get number of newly activated droplets
  pncnuc(1:kproma,:) = MERGE(zcdnc_bas(1:kproma,:) - pcdnc(1:kproma,:), pncnuc(1:kproma,:), ll2(1:kproma,:))
  ! do not allow unphysical decrease in cloud droplets due to activation
  ! pncnuc(1:kproma,:) = MAX(0._dp, pncnuc(1:kproma,:))
  ! do not nucleate droplets in cirrus clouds
  pncnuc(1:kproma,:) = MERGE(pncnuc(1:kproma,:), 0._dp, pt(1:kproma,:) > cthomi)

END SUBROUTINE aerosol_activation

SUBROUTINE subtimestep_1d(&
              !--IN
              kbdim, kproma, ptmst, pntmst_rcp, &
              ll_cc, ll_cci, paclc, paclci, &
              prho, plsdcp, plvdcp, &
              !--INOUT
              ! local variables:
              pq, pt, pxl, pcdnc, pxi, picnc, pqirim, pbirim, &
              pqihet, pqiliq, pqioliq, pqsrc, pqprc, &
              pnihet, pnihom, pninuc, pnidet, &
              ! local tendencies:
              pqte, ptte, pxlte, pcdncte, pxite, &
              picncte, pqirimte, pbirimte, pqihette, pqiliqte, pqioliqte, &
              pqsrcte, pqprcte, &
              pnihette, pnihomte, pninucte, pnidette, &
              ! global tendencies:
              ppqte, pptte, ppxlte, ppcdncte, &
              ppxite, ppicncte, ppqirimte, ppbirimte,&
              ppqihette, ppqiliqte, ppqioliqte, &
              ppqsrcte, ppqprcte, &
              ppnihette, ppnihomte, ppninucte, ppnidette, &
              !--OUT
              pxlb, pcdncb, &
              pxib, picncb, pqirimb, pbirimb, pqihetb, pqiliqb, pqioliqb, ptv)

  INTEGER,  INTENT(in)                      :: kbdim, kproma
  REAL(dp), INTENT(in)                      :: ptmst, pntmst_rcp
  LOGICAL,  INTENT(in), DIMENSION(kbdim)    :: ll_cc, ll_cci
  REAL(dp), INTENT(in), DIMENSION(kbdim)    :: paclc, paclci, prho, &
                                               plsdcp, plvdcp

! local variables
  REAL(dp), INTENT(inout), DIMENSION(kbdim) :: pq, pt, pxl, pcdnc, pxi,  &
                                               picnc, pqirim, pbirim, pqihet, pqiliq, &
                                               pqioliq, pqsrc, pqprc, &
                                               pnihet, pnihom, pninuc, pnidet
! local tendencies
  REAL(dp), INTENT(inout), DIMENSION(kbdim) :: pqte, ptte, pxlte, pcdncte, &
                                               pxite, picncte, pqirimte, &
                                               pbirimte, pqihette, pqiliqte, pqioliqte, &
                                               pqsrcte, pqprcte, &
                                               pnihette, pnihomte, pninucte, pnidette
! global tendencies
  REAL(dp), INTENT(inout), DIMENSION(kbdim) :: ppqte, pptte, ppxlte, ppcdncte, &
                                               ppxite, ppicncte, ppqirimte, ppbirimte, &
                                               ppqihette, ppqiliqte, ppqioliqte, &
                                               ppqsrcte, ppqprcte, &
                                               ppnihette, ppnihomte, ppninucte, ppnidette

! output
  REAL(dp), INTENT(out), DIMENSION(kbdim)   :: pxlb, pcdncb, &
                                               pxib, picncb, pqirimb, pbirimb, &
                                               pqihetb, pqiliqb, pqioliqb, ptv

! LOCAL
  REAL(dp), DIMENSION(kbdim) :: ztmp1, ztmp2
  LOGICAL, DIMENSION(kbdim) :: ll_ice, ll_liq, ll

! assert that there are no negative secondary tracers
! due to conservation laws, this should not be applied to liquid/ice numbers and masses
  pqirimte(1:kproma) = MAX(pqirimte(1:kproma), -pqirim(1:kproma)/ptmst)
  pbirimte(1:kproma) = MAX(pbirimte(1:kproma), -pbirim(1:kproma)/ptmst)
  pqihette(1:kproma) = MAX(pqihette(1:kproma), -pqihet(1:kproma)/ptmst)
  pqiliqte(1:kproma) = MAX(pqiliqte(1:kproma), -pqiliq(1:kproma)/ptmst)
  pqioliqte(1:kproma) = MAX(pqioliqte(1:kproma), -pqioliq(1:kproma)/ptmst)
  pqsrcte(1:kproma) = MAX(pqsrcte(1:kproma), -pqsrc(1:kproma)/ptmst)
  pqprcte(1:kproma) = MAX(pqprcte(1:kproma), -pqprc(1:kproma)/ptmst)

  pnihette(1:kproma) = MAX(pnihette(1:kproma), -pnihet(1:kproma)/ptmst)
  pnihomte(1:kproma) = MAX(pnihomte(1:kproma), -pnihom(1:kproma)/ptmst)
  pninucte(1:kproma) = MAX(pninucte(1:kproma), -pninuc(1:kproma)/ptmst)
  pnidette(1:kproma) = MAX(pnidette(1:kproma), -pnidet(1:kproma)/ptmst)

  ! check if diabatic heating is allowed (turned off for numerical tests in SCM)
  pqte(1:kproma) = MERGE(pqte(1:kproma), 0._dp, ldiabheat)
  ptte(1:kproma) = MERGE(ptte(1:kproma), 0._dp, ldiabheat)

  IF(l2moment) THEN
     pqirimte(1:kproma) = 0._dp
     pbirimte(1:kproma) = 0._dp
  ENDIF !l2moment

! tendencies
  ppxlte(1:kproma) = ppxlte(1:kproma) + pxlte(1:kproma)*pntmst_rcp
  ppcdncte(1:kproma) = ppcdncte(1:kproma) + pcdncte(1:kproma)*pntmst_rcp

  ppxite(1:kproma) = ppxite(1:kproma) + pxite(1:kproma)*pntmst_rcp
  ppicncte(1:kproma) = ppicncte(1:kproma) + picncte(1:kproma)*pntmst_rcp
  ppqirimte(1:kproma) = ppqirimte(1:kproma) + pqirimte(1:kproma)*pntmst_rcp
  ppbirimte(1:kproma) = ppbirimte(1:kproma) + pbirimte(1:kproma)*pntmst_rcp
  ppqihette(1:kproma) = ppqihette(1:kproma) + pqihette(1:kproma)*pntmst_rcp
  ppqiliqte(1:kproma) = ppqiliqte(1:kproma) + pqiliqte(1:kproma)*pntmst_rcp
  ppqioliqte(1:kproma) = ppqioliqte(1:kproma) + pqioliqte(1:kproma)*pntmst_rcp

  ppqsrcte(1:kproma) = ppqsrcte(1:kproma) + pqsrcte(1:kproma)*pntmst_rcp
  ppqprcte(1:kproma) = ppqprcte(1:kproma) + pqprcte(1:kproma)*pntmst_rcp

  ppnihette(1:kproma) = ppnihette(1:kproma) + pnihette(1:kproma)*pntmst_rcp
  ppnihomte(1:kproma) = ppnihomte(1:kproma) + pnihomte(1:kproma)*pntmst_rcp
  ppninucte(1:kproma) = ppninucte(1:kproma) + pninucte(1:kproma)*pntmst_rcp
  ppnidette(1:kproma) = ppnidette(1:kproma) + pnidette(1:kproma)*pntmst_rcp

  ppqte(1:kproma) = ppqte(1:kproma) + pqte(1:kproma)*pntmst_rcp
  pptte(1:kproma) = pptte(1:Kproma) + ptte(1:kproma)*pntmst_rcp

! grid-mean
  pxi(1:kproma) = pxi(1:kproma) + ptmst*pxite(1:kproma)
  picnc(1:kproma) = picnc(1:kproma) + ptmst*picncte(1:kproma)
  pqirim(1:kproma) = pqirim(1:kproma) + ptmst*pqirimte(1:kproma)
  pbirim(1:kproma) = pbirim(1:kproma) + ptmst*pbirimte(1:kproma)
  pqihet(1:kproma) = pqihet(1:kproma) + ptmst*pqihette(1:kproma)
  pqiliq(1:kproma) = pqiliq(1:kproma) + ptmst*pqiliqte(1:kproma)
  pqioliq(1:kproma) = pqioliq(1:kproma) + ptmst*pqioliqte(1:kproma)

  pqsrc(1:kproma) = pqsrc(1:kproma) + ptmst*pqsrcte(1:kproma)
  pqprc(1:kproma) = pqprc(1:kproma) + ptmst*pqprcte(1:kproma)

  pnihet(1:kproma) = pnihet(1:kproma) + ptmst*pnihette(1:kproma)
  pnihom(1:kproma) = pnihom(1:kproma) + ptmst*pnihomte(1:kproma)
  pninuc(1:kproma) = pninuc(1:kproma) + ptmst*pninucte(1:kproma)
  pnidet(1:kproma) = pnidet(1:kproma) + ptmst*pnidette(1:kproma)

  pxl(1:kproma) = pxl(1:kproma) + ptmst*pxlte(1:kproma)
  pcdnc(1:kproma) = pcdnc(1:kproma) + ptmst*pcdncte(1:kproma)

  pq(1:kproma)  = pq(1:kproma) + ptmst*pqte(1:kproma)
  pt(1:kproma)  = pt(1:kproma) + ptmst*ptte(1:kproma)

! in-cloud
  ztmp1(1:kproma) = 1._dp/MAX(paclc(1:kproma),clc_min)
  pxlb(1:kproma) = MERGE(pxl(1:kproma)*ztmp1(1:kproma), 0._dp, ll_cc(1:kproma))
  pcdncb(1:kproma) = MERGE(pcdnc(1:kproma)*ztmp1(1:kproma), 1._dp, ll_cc(1:kproma))

  ztmp1(1:kproma) = 1._dp/MAX(paclci(1:kproma),clc_min)
  pxib(1:kproma) = MERGE(pxi(1:kproma)*ztmp1(1:kproma), 0._dp, ll_cci(1:kproma))
  picncb(1:kproma) = MERGE(picnc(1:kproma)*ztmp1(1:kproma), 0._dp, ll_cci(1:kproma))
  pqirimb(1:kproma) = MERGE(pqirim(1:kproma)*ztmp1(1:kproma), 0._dp, ll_cci(1:kproma))
  pbirimb(1:kproma) = MERGE(pbirim(1:kproma)*ztmp1(1:kproma), 0._dp, ll_cci(1:kproma))
  pqihetb(1:kproma) = MERGE(pqihet(1:kproma)*ztmp1(1:kproma), 0._dp, ll_cci(1:kproma))
  pqiliqb(1:kproma) = MERGE(pqiliq(1:kproma)*ztmp1(1:kproma), 0._dp, ll_cci(1:kproma))
  pqioliqb(1:kproma) = MERGE(pqioliq(1:kproma)*ztmp1(1:kproma), 0._dp, ll_cci(1:kproma))

! virtual temperature
  ptv(1:kproma) = pt(1:kproma)*(1._dp+vtmpc1*pq(1:kproma)-(pxl(1:kproma)+pxi(1:kproma)))


! reset summed tendencies
  pxite(1:kproma) = 0._dp
  picncte(1:kproma) = 0._dp
  pqirimte(1:kproma) = 0._dp
  pbirimte(1:kproma) = 0._dp
  pqihette(1:kproma) = 0._dp
  pqiliqte(1:kproma) = 0._dp
  pqioliqte(1:kproma) = 0._dp
  pqsrcte(1:kproma) = 0._dp
  pqprcte(1:kproma) = 0._dp
  pnihette(1:kproma) = 0._dp
  pnihomte(1:kproma) = 0._dp
  pninucte(1:kproma) = 0._dp
  pxlte(1:kproma) = 0._dp
  pcdncte(1:kproma) = 0._dp
  pqte(1:kproma) = 0._dp
  ptte(1:kproma) = 0._dp
  pnidette(1:kproma) = 0._dp

END SUBROUTINE subtimestep_1d

SUBROUTINE subtimestep_2d(&
              !--IN
              kbdim, klev, kproma, ptmst, pntmst_rcp, &
              ll_cc, ll_cci, paclc, paclci, &
              prho, plsdcp, plvdcp, &
              !--INOUT
              ! local variables:
              pq, pt, pxl, pcdnc, pxi, picnc, pqirim, pbirim, &
              pqihet, pqiliq, pqioliq, pqsrc, pqprc, &
              pnihet, pnihom, pninuc, pnidet, &
              ! local tendencies:
              pqte, ptte, pxlte, pcdncte, pxite, &
              picncte, pqirimte, pbirimte, pqihette, pqiliqte, pqioliqte, &
              pqsrcte, pqprcte, &
              pnihette, pnihomte, pninucte, pnidette, &
              ! global tendencies:
              ppqte, pptte, ppxlte, ppcdncte, &
              ppxite, ppicncte, ppqirimte, ppbirimte,&
              ppqihette, ppqiliqte, ppqioliqte, &
              ppqsrcte, ppqprcte, &
              ppnihette, ppnihomte, ppninucte, ppnidette, &
              !--OUT
              pxlb, pcdncb, &
              pxib, picncb, pqirimb, pbirimb, pqihetb, pqiliqb, pqioliqb, ptv)

  INTEGER,  INTENT(in)                           :: kbdim, kproma, klev
  REAL(dp), INTENT(in)                           :: ptmst, pntmst_rcp
  LOGICAL,  INTENT(in), DIMENSION(kbdim,klev)    :: ll_cc, ll_cci
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)    :: paclc, paclci, prho, &
                                                    plsdcp, plvdcp

! local variables
  REAL(dp), INTENT(inout), DIMENSION(kbdim,klev) :: pq, pt, pxl, pcdnc, pxi,  &
                                                    picnc, pqirim, pbirim, &
                                                    pqihet, pqiliq, pqioliq, &
                                                    pqsrc, pqprc, &
                                                    pnihet, pnihom, pninuc, pnidet
! local tendencies
  REAL(dp), INTENT(inout), DIMENSION(kbdim,klev) :: pqte, ptte, pxlte, pcdncte, &
                                                    pxite, picncte, pqirimte, &
                                                    pbirimte, pqihette, pqiliqte, pqioliqte, &
                                                    pqsrcte, pqprcte, &
                                                    pnihette, pnihomte, pninucte, pnidette
! global tendencies
  REAL(dp), INTENT(inout), DIMENSION(kbdim,klev) :: ppqte, pptte, ppxlte, ppcdncte, &
                                                    ppxite, ppicncte, ppqirimte, ppbirimte, &
                                                    ppqihette, ppqiliqte, ppqioliqte, &
                                                    ppqsrcte, ppqprcte, &
                                                    ppnihette, ppnihomte, ppninucte, ppnidette

! output
  REAL(dp), INTENT(out), DIMENSION(kbdim,klev)   :: pxlb, pcdncb, &
                                                    pxib, picncb, pqirimb, pbirimb, &
                                                    pqihetb, pqiliqb, pqioliqb, ptv

! LOCAL
  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1, ztmp2
  LOGICAL, DIMENSION(kbdim,klev) :: ll_ice, ll_liq, ll

! assert that there are no negative secondary tracers
! due to conservation laws, this should not be applied to liquid/ice numbers and masses
  pqirimte(1:kproma,:) = MAX(pqirimte(1:kproma,:), -pqirim(1:kproma,:)/ptmst)
  pbirimte(1:kproma,:) = MAX(pbirimte(1:kproma,:), -pbirim(1:kproma,:)/ptmst)
  pqihette(1:kproma,:) = MAX(pqihette(1:kproma,:), -pqihet(1:kproma,:)/ptmst)
  pqiliqte(1:kproma,:) = MAX(pqiliqte(1:kproma,:), -pqiliq(1:kproma,:)/ptmst)
  pqioliqte(1:kproma,:) = MAX(pqioliqte(1:kproma,:), -pqioliq(1:kproma,:)/ptmst)

  pqsrcte(1:kproma,:) = MAX(pqsrcte(1:kproma,:), -pqsrc(1:kproma,:)/ptmst)
  pqprcte(1:kproma,:) = MAX(pqprcte(1:kproma,:), -pqprc(1:kproma,:)/ptmst)

  pnihette(1:kproma,:) = MAX(pnihette(1:kproma,:), -pnihet(1:kproma,:)/ptmst)
  pnihomte(1:kproma,:) = MAX(pnihomte(1:kproma,:), -pnihom(1:kproma,:)/ptmst)
  pninucte(1:kproma,:) = MAX(pninucte(1:kproma,:), -pninuc(1:kproma,:)/ptmst)
  pnidette(1:kproma,:) = MAX(pnidette(1:kproma,:), -pnidet(1:kproma,:)/ptmst)

  ! check if diabatic heating is allowed (turned off for numerical tests)
  pqte(1:kproma,:) = MERGE(pqte(1:kproma,:), 0._dp, ldiabheat)
  ptte(1:kproma,:) = MERGE(ptte(1:kproma,:), 0._dp, ldiabheat)

  IF(l2moment) THEN
     pqirimte(1:kproma,:) = 0._dp
     pbirimte(1:kproma,:) = 0._dp
  ENDIF !l2moment

! tendencies
  ppxlte(1:kproma,:) = ppxlte(1:kproma,:) + pxlte(1:kproma,:)*pntmst_rcp
  ppcdncte(1:kproma,:) = ppcdncte(1:kproma,:) + pcdncte(1:kproma,:)*pntmst_rcp

  ppxite(1:kproma,:) = ppxite(1:kproma,:) + pxite(1:kproma,:)*pntmst_rcp
  ppicncte(1:kproma,:) = ppicncte(1:kproma,:) + picncte(1:kproma,:)*pntmst_rcp
  ppqirimte(1:kproma,:) = ppqirimte(1:kproma,:) + pqirimte(1:kproma,:)*pntmst_rcp
  ppbirimte(1:kproma,:) = ppbirimte(1:kproma,:) + pbirimte(1:kproma,:)*pntmst_rcp
  ppqihette(1:kproma,:) = ppqihette(1:kproma,:) + pqihette(1:kproma,:)*pntmst_rcp
  ppqiliqte(1:kproma,:) = ppqiliqte(1:kproma,:) + pqiliqte(1:kproma,:)*pntmst_rcp
  ppqioliqte(1:kproma,:) = ppqioliqte(1:kproma,:) + pqioliqte(1:kproma,:)*pntmst_rcp

  ppqsrcte(1:kproma,:) = ppqsrcte(1:kproma,:) + pqsrcte(1:kproma,:)*pntmst_rcp
  ppqprcte(1:kproma,:) = ppqprcte(1:kproma,:) + pqprcte(1:kproma,:)*pntmst_rcp

  ppnihette(1:kproma,:) = ppnihette(1:kproma,:) + pnihette(1:kproma,:)*pntmst_rcp
  ppnihomte(1:kproma,:) = ppnihomte(1:kproma,:) + pnihomte(1:kproma,:)*pntmst_rcp
  ppninucte(1:kproma,:) = ppninucte(1:kproma,:) + pninucte(1:kproma,:)*pntmst_rcp
  ppnidette(1:kproma,:) = ppnidette(1:kproma,:) + pnidette(1:kproma,:)*pntmst_rcp

  ppqte(1:kproma,:) = ppqte(1:kproma,:) + pqte(1:kproma,:)*pntmst_rcp
  pptte(1:kproma,:) = pptte(1:kproma,:) + ptte(1:kproma,:)*pntmst_rcp

! grid-mean
  pxi(1:kproma,:) = pxi(1:kproma,:) + ptmst*pxite(1:kproma,:)
  picnc(1:kproma,:) = picnc(1:kproma,:) + ptmst*picncte(1:kproma,:)
  pqirim(1:kproma,:) = pqirim(1:kproma,:) + ptmst*pqirimte(1:kproma,:)
  pbirim(1:kproma,:) = pbirim(1:kproma,:) + ptmst*pbirimte(1:kproma,:)
  pqihet(1:kproma,:) = pqihet(1:kproma,:) + ptmst*pqihette(1:kproma,:)
  pqiliq(1:kproma,:) = pqiliq(1:kproma,:) + ptmst*pqiliqte(1:kproma,:)
  pqioliq(1:kproma,:) = pqioliq(1:kproma,:) + ptmst*pqioliqte(1:kproma,:)

  pqsrc(1:kproma,:) = pqsrc(1:kproma,:) + ptmst*pqsrcte(1:kproma,:)
  pqprc(1:kproma,:) = pqprc(1:kproma,:) + ptmst*pqprcte(1:kproma,:)

  pnihet(1:kproma,:) = pnihet(1:kproma,:) + ptmst*pnihette(1:kproma,:)
  pnihom(1:kproma,:) = pnihom(1:kproma,:) + ptmst*pnihomte(1:kproma,:)
  pninuc(1:kproma,:) = pninuc(1:kproma,:) + ptmst*pninucte(1:kproma,:)
  pnidet(1:kproma,:) = pnidet(1:kproma,:) + ptmst*pnidette(1:kproma,:)

  pxl(1:kproma,:) = pxl(1:kproma,:) + ptmst*pxlte(1:kproma,:)
  pcdnc(1:kproma,:) = pcdnc(1:kproma,:) + ptmst*pcdncte(1:kproma,:)

  pq(1:kproma,:)  = pq(1:kproma,:) + ptmst*pqte(1:kproma,:)
  pt(1:kproma,:)  = pt(1:kproma,:) + ptmst*ptte(1:kproma,:)

! in-cloud
  ztmp1(1:kproma,:) = 1._dp/MAX(paclc(1:kproma,:),clc_min)
  pxlb(1:kproma,:) = MERGE(pxl(1:kproma,:)*ztmp1(1:kproma,:), 0._dp, ll_cc(1:kproma,:))
  pcdncb(1:kproma,:) = MERGE(pcdnc(1:kproma,:)*ztmp1(1:kproma,:), 1._dp, ll_cc(1:kproma,:))

  ztmp1(1:kproma,:) = 1._dp/MAX(paclci(1:kproma,:),clc_min)
  pxib(1:kproma,:) = MERGE(pxi(1:kproma,:)*ztmp1(1:kproma,:), 0._dp, ll_cci(1:kproma,:))
  picncb(1:kproma,:) = MERGE(picnc(1:kproma,:)*ztmp1(1:kproma,:), 0._dp, ll_cci(1:kproma,:))
  pqirimb(1:kproma,:) = MERGE(pqirim(1:kproma,:)*ztmp1(1:kproma,:), 0._dp, ll_cci(1:kproma,:))
  pbirimb(1:kproma,:) = MERGE(pbirim(1:kproma,:)*ztmp1(1:kproma,:), 0._dp, ll_cci(1:kproma,:))
  pqihetb(1:kproma,:) = MERGE(pqihet(1:kproma,:)*ztmp1(1:kproma,:), 0._dp, ll_cci(1:kproma,:))
  pqiliqb(1:kproma,:) = MERGE(pqiliq(1:kproma,:)*ztmp1(1:kproma,:), 0._dp, ll_cci(1:kproma,:))
  pqioliqb(1:kproma,:) = MERGE(pqioliq(1:kproma,:)*ztmp1(1:kproma,:), 0._dp, ll_cci(1:kproma,:))

! virtual temperature
  ptv(1:kproma,:) = pt(1:kproma,:)*(1._dp+vtmpc1*pq(1:kproma,:)-(pxl(1:kproma,:)+pxi(1:kproma,:)))

! reset summed tendencies
  pxite(1:kproma,:) = 0._dp
  picncte(1:kproma,:) = 0._dp
  pqirimte(1:kproma,:) = 0._dp
  pbirimte(1:kproma,:) = 0._dp
  pqihette(1:kproma,:) = 0._dp
  pqiliqte(1:kproma,:) = 0._dp
  pqioliqte(1:kproma,:) = 0._dp
  pqsrcte(1:kproma,:) = 0._dp
  pqprcte(1:kproma,:) = 0._dp
  pnihette(1:kproma,:) = 0._dp
  pnihomte(1:kproma,:) = 0._dp
  pninucte(1:kproma,:) = 0._dp
  pnidette(1:kproma,:) = 0._dp
  pxlte(1:kproma,:) = 0._dp
  pcdncte(1:kproma,:) = 0._dp
  pqte(1:kproma,:) = 0._dp
  ptte(1:kproma,:) = 0._dp

END SUBROUTINE subtimestep_2d

SUBROUTINE read_fall_velocity_2d(&
              !--IN
              kbdim, kproma, klev,  &
              pxi, picnc, pqifrac, prhop, pbirim, &
              prho_rcp, prhofaci, vt_rdc, &
              !--OUT
              pvtin, pvtim, prim)

  INTEGER, INTENT(in) :: kbdim, kproma, klev
  REAL(dp), INTENT(in) :: vt_rdc
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev) :: pxi, picnc, pqifrac, prhop, pbirim ! icnc in m-3
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev) :: prho_rcp, prhofaci
  REAL(dp), INTENT(out), DIMENSION(kbdim,klev) :: pvtim, pvtin, prim

  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1, zlkp1, zlkp2, zlkp11
  REAL(dp), DIMENSION(kbdim,klev) :: dum1,dum2,dum3,dum4,dum5
  INTEGER,  DIMENSION(kbdim,klev) :: dumi,dumk,dumii,dumjj
  LOGICAL :: ll1(kbdim,klev)
  INTEGER :: jl,jk

  CALL calculate_lookup_table_indices(&
       ! IN
       kbdim,klev,kproma,                                                           &
       prhop(:,:),pqifrac(:,:),pxi(:,:),picnc(:,:),     &
       ! OUT
       dum1(:,:), dum2(:,:), dum4(:,:), dum5(:,:),      &
       dumi(:,:), dumk(:,:), dumii(:,:), dumjj(:,:))

  DO jl=1,kproma
     DO jk=1,klev
        CALL access_lookup_table(&
             ! IN
             dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
             dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),1,    &
             ! OUT
             zlkp1(jl,jk))
        CALL access_lookup_table(&
             ! IN
             dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
             dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),2,    &
             ! OUT
             zlkp2(jl,jk))
        CALL access_lookup_table(&
             ! IN
             dum1(jl,jk),dumi(jl,jk),dum2(jl,jk),dumk(jl,jk),         &
             dum4(jl,jk),dumii(jl,jk),dum5(jl,jk),dumjj(jl,jk),11,    &
             ! OUT
             zlkp11(jl,jk))
    ENDDO !jk
  ENDDO !jl

! ice present logical
  ll1(1:kproma,:) = pxi(1:kproma,:) .gt. cqtmin
  pvtin(1:kproma,:) = MERGE(prhofaci(1:kproma,:)*zlkp1(1:kproma,:)*vt_rdc, 0._dp, ll1(1:kproma,:))
  pvtim(1:kproma,:) = MERGE(prhofaci(1:kproma,:)*zlkp2(1:kproma,:)*vt_rdc, 0._dp, ll1(1:kproma,:))
  prim(1:kproma,:) = MERGE(0.5*zlkp11(1:kproma,:), 0._dp, ll1(1:kproma,:))


END SUBROUTINE read_fall_velocity_2d

SUBROUTINE read_fall_velocity(&
              !--IN
              kbdim, kproma,  &
              pxi, picnc, pqifrac, prhop, pbirim, &
              prho_rcp, prhofaci, vt_rdc, &
              !--OUT
              pvtin, pvtim, prim)

  INTEGER, INTENT(in) :: kbdim, kproma
  REAL(dp), INTENT(in) :: vt_rdc
  REAL(dp), INTENT(in), DIMENSION(kbdim) :: pxi, picnc, pqifrac, prhop, pbirim ! icnc in m-3
  REAL(dp), INTENT(in), DIMENSION(kbdim) :: prho_rcp, prhofaci
  REAL(dp), INTENT(out), DIMENSION(kbdim) :: pvtim, pvtin, prim

  REAL(dp), DIMENSION(kbdim) :: ztmp1, zlkp1, zlkp2, zlkp11
  REAL(dp), DIMENSION(kbdim) :: dum1,dum2,dum3,dum4,dum5
  INTEGER,  DIMENSION(kbdim) :: dumi,dumk,dumii,dumjj
  LOGICAL :: ll1(kbdim)
  INTEGER :: jl

  CALL calculate_lookup_table_indices_1d(&
       ! IN
       kbdim,kproma,                                                        &
       prhop(:),pqifrac(:),pxi(:),picnc(:),     &
       ! OUT
       dum1(:), dum2(:), dum4(:), dum5(:),      &
       dumi(:), dumk(:), dumii(:), dumjj(:))

  DO jl=1,kproma
    CALL access_lookup_table(&
         ! IN
         dum1(jl),dumi(jl),dum2(jl),dumk(jl),         &
         dum4(jl),dumii(jl),dum5(jl),dumjj(jl),1,    &
         ! OUT
         zlkp1(jl))
    CALL access_lookup_table(&
         ! IN
         dum1(jl),dumi(jl),dum2(jl),dumk(jl),         &
         dum4(jl),dumii(jl),dum5(jl),dumjj(jl),2,    &
         ! OUT
         zlkp2(jl))
    CALL access_lookup_table(&
         ! IN
         dum1(jl),dumi(jl),dum2(jl),dumk(jl),         &
         dum4(jl),dumii(jl),dum5(jl),dumjj(jl),11,    &
         ! OUT
         zlkp11(jl))
  ENDDO

! ice present logical
  ll1(1:kproma) = pxi(1:kproma) .gt. cqtmin
  pvtin(1:kproma) = MERGE(prhofaci(1:kproma)*zlkp1(1:kproma)*vt_rdc, 0._dp, ll1(1:kproma))
  pvtim(1:kproma) = MERGE(prhofaci(1:kproma)*zlkp2(1:kproma)*vt_rdc, 0._dp, ll1(1:kproma))
  prim(1:kproma) = MERGE(0.5*zlkp11(1:kproma), 0._dp, ll1(1:kproma))

END SUBROUTINE read_fall_velocity

SUBROUTINE read_my_fall_velocity_2d(&
              !--IN
              kbdim, kproma, klev, &
              pxi, picnc, primfrac, prhop, pbirim, &
              prho_rcp, prhofaci, vt_rdc, &
              !--OUT
              pvtin, pvtim, prim)

  INTEGER, INTENT(in) :: kbdim, kproma, klev
  REAL(dp), INTENT(in) :: vt_rdc
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev) :: pxi, picnc, primfrac, prhop, pbirim ! icnc in m-3
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev) :: prho_rcp, prhofaci
  REAL(dp), INTENT(out), DIMENSION(kbdim,klev) :: pvtim, pvtin, prim

  REAL(dp), DIMENSION(kbdim,klev) :: qnorm, iqnorm, irhor, ifr
  REAL(dp), DIMENSION(kbdim,klev) :: zlkp3, zlkp6, zlkp7
  LOGICAL :: ll1(kbdim,klev)
  INTEGER :: jl

!>>DN bugfix
!  qnorm = MERGE(pxi(1:kproma,:)/picnc(1:kproma,:), 0._dp, picnc(1:kproma,:) > cqtmin)
  qnorm(1:kproma,:) = MERGE(pxi(1:kproma,:)/picnc(1:kproma,:), 0._dp, picnc(1:kproma,:) > cqtmin)
!<<DN bugfix

  CALL calculate_my_lookup_table_indices(&
          !--IN
          kproma, kbdim, klev, &
          qnorm, primfrac, prhop, &
          !--OUT
          iqnorm, irhor, ifr)

  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 3, &
          !--OUT
          zlkp3)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 6, &
          !--OUT
          zlkp6)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, klev, &
          iqnorm, ifr, irhor, 7, &
          !--OUT
          zlkp7)

! ice present logical
  ll1(1:kproma,:) = pxi(1:kproma,:) .gt. cqtmin
  pvtim(1:kproma,:) = MERGE(prhofaci(1:kproma,:)*zlkp6(1:kproma,:)*vt_rdc, 0._dp, ll1(1:kproma,:))
  pvtin(1:kproma,:) = MERGE(prhofaci(1:kproma,:)*zlkp7(1:kproma,:)*vt_rdc, 0._dp, ll1(1:kproma,:))
  prim(1:kproma,:) = MERGE(0.5*zlkp3(1:kproma,:), 0._dp, ll1(1:kproma,:))

END SUBROUTINE read_my_fall_velocity_2d

SUBROUTINE read_my_fall_velocity(&
              !--IN
              kbdim, kproma,  &
              pxi, picnc, primfrac, prhop, pbirim, &
              prho_rcp, prhofaci, vt_rdc, &
              !--OUT
              pvtin, pvtim, prim)

  INTEGER, INTENT(in) :: kbdim, kproma
  REAL(dp), INTENT(in) :: vt_rdc
  REAL(dp), INTENT(in), DIMENSION(kbdim) :: pxi, picnc, primfrac, prhop, pbirim ! icnc in m-3
  REAL(dp), INTENT(in), DIMENSION(kbdim) :: prho_rcp, prhofaci
  REAL(dp), INTENT(out), DIMENSION(kbdim) :: pvtim, pvtin, prim

  REAL(dp), DIMENSION(kbdim) :: qnorm, iqnorm, irhor, ifr
  REAL(dp), DIMENSION(kbdim) :: zlkp3, zlkp6, zlkp7
  LOGICAL :: ll1(kbdim)
  INTEGER :: jl

!>>DN bugfix
!  qnorm = MERGE(pxi(1:kproma)/picnc(1:kproma), 0._dp, picnc(1:kproma) > cqtmin)
  qnorm(1:kproma) = MERGE(pxi(1:kproma)/picnc(1:kproma), 0._dp, picnc(1:kproma) > cqtmin)
!<<DN bugfix

  CALL calculate_my_lookup_table_indices(&
          !--IN
          kproma, kbdim, &
          qnorm, primfrac, prhop, &
          !--OUT
          iqnorm, irhor, ifr)

  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, &
          iqnorm, ifr, irhor, 3, &
          !--OUT
          zlkp3)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, &
          iqnorm, ifr, irhor, 6, &
          !--OUT
          zlkp6)
  CALL access_my_lookup_table(&
          !--IN
          kproma, kbdim, &
          iqnorm, ifr, irhor, 7, &
          !--OUT
          zlkp7)

! ice present logical
  ll1(1:kproma) = pxi(1:kproma) .gt. cqtmin
  pvtim(1:kproma) = MERGE(prhofaci(1:kproma)*zlkp6(1:kproma)*vt_rdc, 0._dp, ll1(1:kproma))
  pvtin(1:kproma) = MERGE(prhofaci(1:kproma)*zlkp7(1:kproma)*vt_rdc, 0._dp, ll1(1:kproma))
  prim(1:kproma) = MERGE(0.5*zlkp3(1:kproma), 0._dp, ll1(1:kproma))

END SUBROUTINE read_my_fall_velocity

SUBROUTINE ice_sedimentation(&
              !--IN
              kbdim, klev, kproma, ptmst, ptmst_rcp, &
              pdpg, pvtim, pvtin, pdz, pt, &
              pvmpot, pvnpot, prho, prhofaci, &
              pxi, picnc, pqirim, pbirim, pqihet, pqiliq, pqioliq, &
              pnihet, pnihom, pninuc, pnidet, &
              !--OUT
              qisten, qristen, bgsten, nisten, qhsten, qlsten, qlosten, &
              nihetsten, nihomsten, ninucsten, nidetsten, &
              nisflx, qisflx, qiflx, nstep)

  INTEGER, INTENT(in)                          :: kbdim, klev, kproma
  REAL(dp), INTENT(in)                         :: ptmst, ptmst_rcp
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)  :: pdpg, pvmpot, pvnpot, prho, &
                                                  pxi, picnc, pqirim, pbirim, pqihet, pqiliq, pqioliq, &
                                                  pnihet, pnihom, pninuc, pnidet, &
                                                  pvtim, pvtin, pdz, pt
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)  :: prhofaci
  REAL(dp), INTENT(out), DIMENSION(kbdim,klev) :: qisten, qristen, bgsten, qhsten, qlsten, nisten, qiflx, qlosten, &
                                                  nihetsten, nihomsten, ninucsten, nidetsten
  REAL(dp), INTENT(out), DIMENSION(kbdim)      :: nisflx, qisflx
  INTEGER,  INTENT(out)                        :: nstep !< number of substeps

  REAL(dp) :: amax          !< highest sedimentation velocity in column (m and n)
  REAL(dp) :: nstep_rcp     !< reciprocal of nstep

  REAL(dp), DIMENSION(kbdim,klev)   :: zxi, zicnc, zqirim, zbirim, zqihet, & !< local copies of sed. vars.
                                       zqiliq, zqioliq, &
                                       znihet, znihom, zninuc, znidet
  REAL(dp), DIMENSION(kbdim,klev+1) :: falouti, faloutni, faloutri, &
                                       faloutbg, faloutqh, faloutqwc
  REAL(dp), DIMENSION(kbdim,klev)   :: faltndi, faltndni, faltndri, &
                                       faltndbg, faltndqh, faltndqwc              !< fraction of sedimentation
  REAL(dp), DIMENSION(kbdim,klev)   :: q_m, q_n
  REAL(dp), DIMENSION(kbdim,klev)   :: dumm, dumn

  INTEGER :: jl, jk, n  !< horizontal, vertical and step indices
  REAL(dp), DIMENSION(kbdim,klev) :: zdpg_rcp, zdz_rcp, zrho_rcp
  REAL(dp), DIMENSION(kbdim) :: dumsurf ! we only care about surface flux of qi and ni
  REAL(dp), DIMENSION(kbdim) :: qisurfk, nisurfk ! used to add up surface flux each step
!>>DN bugfix
  REAL(dp), DIMENSION(kbdim,klev+1) :: zout, dummout
!<<DN bugfix

  zdpg_rcp(1:kproma,:) = 1._dp/pdpg(1:kproma,:)
  zdz_rcp(1:kproma,:) = 1._dp/pdz(1:kproma,:)
  zrho_rcp(1:kproma,:) = 1._dp/prho(1:kproma,:)

  qiflx(1:kproma,:) = 0._dp
  qisflx(1:kproma) = 0._dp
  nisflx(1:kproma) = 0._dp

  qisten(1:kproma,:) = 0._dp
  qristen(1:kproma,:) = 0._dp
  bgsten(1:kproma,:) = 0._dp
  nisten(1:kproma,:) = 0._dp
  qhsten(1:kproma,:) = 0._dp
  qlsten(1:kproma,:) = 0._dp
  qlosten(1:kproma,:) = 0._dp
!>>DN bugfix
  qiflx(1:kproma,:) = 0._dp
!<<DN bugfix

  zxi(1:kproma,:) = pxi(1:kproma,:)
  zicnc(1:kproma,:) = picnc(1:kproma,:)
  zqirim(1:kproma,:) = pqirim(1:kproma,:)
  zbirim(1:kproma,:) = pbirim(1:kproma,:)
  zqihet(1:kproma,:) = pqihet(1:kproma,:)
  zqiliq(1:kproma,:) = pqiliq(1:kproma,:)
  zqioliq(1:kproma,:) = pqioliq(1:kproma,:)

  znihet(1:kproma,:) = pnihet(1:kproma,:)
  znihom(1:kproma,:) = pnihom(1:kproma,:)
  zninuc(1:kproma,:) = pninuc(1:kproma,:)
  znidet(1:kproma,:) = pnidet(1:kproma,:)

  amax = MAXVAL(ptmst*pvmpot(1:kproma,:)*zdz_rcp(1:kproma,:))*4
  nstep = MAX(CEILING(amax), 1)
  nstep = MERGE(nsedi, nstep, nsedi > 0)
  nstep_rcp = 1._dp/REAL(nstep)

  q_m(1:kproma,:) = pvmpot(1:kproma,:)*zdz_rcp(1:kproma,:)
  q_m(1:kproma,:) = MAX(MIN(ptmst_rcp*nstep, q_m(1:kproma,:)), 0._dp)
  q_n(1:kproma,:) = pvnpot(1:kproma,:)*zdz_rcp(1:kproma,:)
  q_n(1:kproma,:) = MAX(MIN(ptmst_rcp*nstep, q_n(1:kproma,:)), 0._dp)
  
  DO n=1,nstep

     IF(iintscheme==0) THEN
!>>DN bugfix: this bugfix should be extended to the implicit_euler scheme if it is used
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zxi, qisurfk)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, zicnc, nisurfk)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqirim, dumsurf)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zbirim, dumsurf)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqihet, dumsurf)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqiliq, dumsurf)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqioliq, dumsurf)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, znihet, dumsurf)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, znihom, dumsurf)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, zninuc, dumsurf)
!!$        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, znidet, dumsurf)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zxi, qisurfk,zout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, zicnc, nisurfk,dummout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqirim, dumsurf,dummout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zbirim, dumsurf,dummout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqihet, dumsurf,dummout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqiliq, dumsurf,dummout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqioliq, dumsurf,dummout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, znihet, dumsurf,dummout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, znihom, dumsurf,dummout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, zninuc, dumsurf,dummout)
        CALL euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, znidet, dumsurf,dummout)
!<<DN bugfix
     ELSE IF(iintscheme==1) THEN
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zxi, qisurfk)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, zicnc, nisurfk)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqirim, dumsurf)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zbirim, dumsurf)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqihet, dumsurf)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqiliq, dumsurf)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_m, pdpg, zqioliq, dumsurf)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, znihet, dumsurf)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, znihom, dumsurf)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, zninuc, dumsurf)
        CALL implicit_euler(kbdim, klev, kproma, ptmst*nstep_rcp, q_n, pdpg, znidet, dumsurf)
     ENDIF

     qisflx(1:kproma) = qisflx(1:kproma) + qisurfk(1:kproma)*nstep_rcp
     nisflx(1:kproma) = nisflx(1:kproma) + nisurfk(1:kproma)*nstep_rcp
!>>DN bugfix
     qiflx(1:kproma,1:klev) = qiflx(1:kproma,1:klev) + &
          zout(1:kproma,2:klev+1)*nstep_rcp!use sedimentation flux, not net outflux
!<<DN bugfix
  ENDDO !n

  ! since the outer loop uses an explicit euler method, convert back to
  ! the corresponding derivative
  qisten(1:kproma,:) = (zxi(1:kproma,:) - pxi(1:kproma,:))/ptmst
  nisten(1:kproma,:) = (zicnc(1:kproma,:) - picnc(1:kproma,:))/ptmst
  qristen(1:kproma,:) = (zqirim(1:kproma,:) - pqirim(1:kproma,:))/ptmst
  bgsten(1:kproma,:) = (zbirim(1:kproma,:) - pbirim(1:kproma,:))/ptmst
  qhsten(1:kproma,:) = (zqihet(1:kproma,:) - pqihet(1:kproma,:))/ptmst
  qlsten(1:kproma,:) = (zqiliq(1:kproma,:) - pqiliq(1:kproma,:))/ptmst
  qlosten(1:kproma,:) = (zqioliq(1:kproma,:) - pqioliq(1:kproma,:))/ptmst
  nihetsten(1:kproma,:) = (znihet(1:kproma,:) - pnihet(1:kproma,:))/ptmst
  nihomsten(1:kproma,:) = (znihom(1:kproma,:) - pnihom(1:kproma,:))/ptmst
  ninucsten(1:kproma,:) = (zninuc(1:kproma,:) - pninuc(1:kproma,:))/ptmst
  nidetsten(1:kproma,:) = (znidet(1:kproma,:) - pnidet(1:kproma,:))/ptmst

!>>DN bugfix
!  ! use net outflux as proxy for scavenging
!  qiflx(1:kproma,:) = MIN(0._dp, qisten(1:kproma,:))
!<<DN bugfix

END SUBROUTINE ice_sedimentation

SUBROUTINE implicit_euler(&
              !--IN
              kbdim, klev, kproma, ptmst, pquot, pdpg, &
              !--INOUT
              pval, &
              !--OUT
              psurf)

  INTEGER, INTENT(IN) :: kbdim, klev, kproma
  REAL(dp), INTENT(IN) :: ptmst
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: pquot, pdpg
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim,klev) :: pval
  REAL(dp), INTENT(OUT), DIMENSION(kbdim) :: psurf

  ! use 'z' values for absolute masses and convert back in the end
  ! note that psurf is given in absolute mass
  REAL(dp), DIMENSION(kbdim,klev) :: zval, zvalp1

  INTEGER :: jk

  ! convert mixing ratio to absolute mass
  zval(1:kproma,:) = pval(1:kproma,:)*pdpg(1:kproma,:)

  zvalp1(1:kproma,:) = 0._dp
  ! top level does not have influx from above
  zvalp1(1:kproma,1) = zval(1:kproma,1)/(1+ptmst*pquot(1:kproma,1))
  ! loop from top to bottom
  DO jk=2,klev
     ! implicit method depends on influx from above at p1
     zvalp1(1:kproma,jk) = (zval(1:kproma,jk)+ptmst*pquot(1:kproma,jk-1)*zvalp1(1:kproma,jk-1)) &
                           /(1+ptmst*pquot(1:kproma,jk))
  END DO ! jk
  ! surface is outflux of lowest level
  psurf(1:kproma) = pquot(1:kproma,klev)*zvalp1(1:kproma,klev)

  ! convert back to mixing ratio
  pval(1:kproma,:) = zvalp1(1:kproma,:)/pdpg(1:kproma,:)

END SUBROUTINE implicit_euler

SUBROUTINE euler(&
              !--IN
              kbdim, klev, kproma, ptmst, pquot, pdpg, &
              !--INOUT
              pval, &
              !--OUT
!>>DN bugfix
!              psurf)
              psurf, zout)
!<<DN bugfix

  INTEGER, INTENT(IN) :: kbdim, klev, kproma
  REAL(dp), INTENT(IN) :: ptmst
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: pquot, pdpg
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim,klev) :: pval
  REAL(dp), INTENT(OUT), DIMENSION(kbdim) :: psurf

!>>DN bugfix
!  REAL(dp), DIMENSION(kbdim,klev+1) :: zout
  REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev+1) :: zout
!<<DN bugfix

  ! use 'z' values for absolute masses and convert back in the end
  ! note that psurf is given in absolute mass
  REAL(dp), DIMENSION(kbdim,klev) :: zval, zvalp1

  INTEGER :: jk

  ! convert mixing ratio to absolute mass
  zval(1:kproma,:) = pval(1:kproma,:)*pdpg(1:kproma,:)

  ! array storing the outflux per level
  ! note that the array is shifted, such that klev+1 is the lowest level
  ! and 1 is above the uppermost level
  zout(1:kproma,2:klev+1) = zval(1:kproma,:)*pquot(1:kproma,:)
  zout(1:kproma,1) = 0._dp
  ! the derivative 
  zvalp1(1:kproma,:) = zval(1:kproma,:) + ptmst*(zout(1:kproma,1:klev)-zout(1:kproma,2:klev+1))

  ! surface is the outflux at lowest level. Note that this is in kg/m2/s
  psurf(1:kproma) = zout(1:kproma,klev+1)

  ! convert back to mixing ratios
  pval(1:kproma,:) = zvalp1(1:kproma,:)/pdpg(1:kproma,:)

END SUBROUTINE euler

SUBROUTINE update_saturation_values_1d(&
              !--IN
              name, tmin, tmax, &
              kproma, kbdim, ptm1, papm1, pqm1,      &
              !--OUT
              pesw, pqsw, psupw, pesi, pqsi, psupi,               &
              pqswp1, pqsip1, pqsw0)

  CHARACTER(len=*)                          :: name
  REAL(dp)                                  :: tmin, tmax
  INTEGER, INTENT(in)                       :: kproma, kbdim
  REAL(dp), DIMENSION(kbdim), INTENT(in)    :: ptm1, papm1, pqm1
  REAL(dp), DIMENSION(kbdim), INTENT(out)   :: pesw, pqsw, psupw, pesi,   &
                                               pqsi, psupi,pqswp1, pqsip1, &
                                               pqsw0

  INTEGER                    :: jl
  REAL(dp), DIMENSION(kbdim) :: ztmp1
  LOGICAL, DIMENSION(kbdim)  :: ll_look
  INTEGER, DIMENSION(kbdim)  :: itm1_look, itm1p1_look, it0_look !< lookup table indices
  
  ztmp1(1:kproma)      = 1000._dp*ptm1(1:kproma)
  itm1_look(1:kproma) = NINT(ztmp1(1:kproma))

  ll_look(1:kproma) = (itm1_look(1:kproma)<jptlucu1 .OR. itm1_look(1:kproma)>jptlucu2)
  IF (ANY(ll_look(1:kproma))) THEN
     WRITE(*,*) 'max temperature: ', MAXVAL(ptm1(1:kproma)), tmax
     WRITE(*,*) 'min temperature: ', MINVAL(ptm1(1:kproma)), tmin
     CALL lookuperror (name)
  ENDIF
  IF (ANY(ll_look(1:kproma))) CALL lookuperror (name)

  itm1_look(1:kproma) = MAX(MIN(itm1_look(1:kproma),jptlucu2),jptlucu1)

  itm1p1_look(1:kproma) = itm1_look(1:kproma) + 1
  itm1p1_look(1:kproma) = MAX(MIN(itm1p1_look(1:kproma),jptlucu2),jptlucu1)

  it0_look(1:kproma) = NINT(1000._dp*tmelt)

  DO jl=1,kproma
!SF water:
     ztmp1(jl) = tlucuaw(itm1_look(jl))/papm1(jl)     !<  e_s,w*Rd/Rv /p
     ztmp1(jl) = MIN(ztmp1(jl),0.5_dp)
     pqsw(jl)  = ztmp1(jl)/(1._dp-vtmpc1*ztmp1(jl))
     psupw(jl) = pqm1(jl)/pqsw(jl)-1.0_dp
     pesw(jl)  = ztmp1(jl)*papm1(jl)*rv/rd

     ztmp1(jl) = tlucuaw(it0_look(jl))/papm1(jl)     !<  e_s,w*Rd/Rv /p
     ztmp1(jl) = MIN(ztmp1(jl),0.5_dp)
     pqsw0(jl)  = ztmp1(jl)/(1._dp-vtmpc1*ztmp1(jl))

     !SF for later use in 5.:
     pqswp1(jl) = tlucuaw(itm1p1_look(jl))/papm1(jl)
     pqswp1(jl) = MIN(pqswp1(jl),0.5_dp)
     pqswp1(jl) = pqswp1(jl)/(1._dp-vtmpc1*pqswp1(jl))
     !SFend for later use in 5.

!SF ice:
     ztmp1(jl) = tlucua(itm1_look(jl))/papm1(jl)       !<  e_s,i*Rd/Rv /p
     ztmp1(jl) = MIN(ztmp1(jl),0.5_dp)
     pqsi(jl)  = ztmp1(jl)/(1._dp-vtmpc1*ztmp1(jl))
     psupi(jl) = pqm1(jl)/pqsi(jl)-1.0_dp
     pesi(jl)  = ztmp1(jl)*papm1(jl)*rv/rd

     !SF for later use in 5.:
     pqsip1(jl) = tlucua(itm1p1_look(jl))/papm1(jl)
     pqsip1(jl) = MIN(pqsip1(jl),0.5_dp)
     pqsip1(jl) = pqsip1(jl)/(1._dp-vtmpc1*pqsip1(jl))
     !SFend for later use in 5.

  END DO !jl
 
END SUBROUTINE update_saturation_values_1d

SUBROUTINE update_saturation_values_2d(&
              !--IN
              name, tmin, tmax, &
              kproma, kbdim, ktdia, klev, ptm1, papm1, pqm1,      &
              !--OUT
              pesw, pqsw, psupw, pesi, pqsi, psupi,               &
              pqswp1, pqsip1, pqsw0)

  CHARACTER(len=*)                               :: name
  REAL(dp)                                       :: tmin, tmax
  INTEGER, INTENT(in)                            :: kproma, kbdim, ktdia, klev
  REAL(dp), DIMENSION(kbdim,klev), INTENT(in)    :: ptm1, papm1, pqm1
  REAL(dp), DIMENSION(kbdim,klev), INTENT(out)   :: pesw, pqsw, psupw, pesi,   &
                                                    pqsi, psupi,pqswp1, pqsip1, &
                                                    pqsw0

  INTEGER                            :: jk, jl
  REAL(dp), DIMENSION(kbdim,klev)    :: ztmp1
  LOGICAL, DIMENSION(kbdim,klev)     :: ll_look
  INTEGER, DIMENSION(kbdim,klev)     :: itm1_look, itm1p1_look, it0_look !< lookup table indices
  
  ztmp1(1:kproma,:)      = 1000._dp*ptm1(1:kproma,:)
  itm1_look(1:kproma,:) = NINT(ztmp1(1:kproma,:))

  ll_look(1:kproma,:) = (itm1_look(1:kproma,:)<jptlucu1 .OR. itm1_look(1:kproma,:)>jptlucu2)
  IF (ANY(ll_look(1:kproma,ktdia:klev))) THEN
     WRITE(*,*) 'max temperature: ', MAXVAL(ptm1(1:kproma,:)), tmax
     WRITE(*,*) 'min temperature: ', MINVAL(ptm1(1:kproma,:)), tmin
     CALL lookuperror (name)
  ENDIF

  itm1_look(1:kproma,:) = MAX(MIN(itm1_look(1:kproma,:),jptlucu2),jptlucu1)

!SF for later use in 5.:
  itm1p1_look(1:kproma,:) = itm1_look(1:kproma,:) + 1
  itm1p1_look(1:kproma,:) = MAX(MIN(itm1p1_look(1:kproma,:),jptlucu2),jptlucu1)
!SFend for later use in 5.

  it0_look(1:kproma,:) = NINT(1000._dp*tmelt)

  DO jk=klev,ktdia,-1
     DO jl=1,kproma
!SF water:
        ztmp1(jl,jk) = tlucuaw(itm1_look(jl,jk))/papm1(jl,jk)     !<  e_s,w*Rd/Rv /p
        ztmp1(jl,jk) = MIN(ztmp1(jl,jk),0.5_dp)
        pqsw(jl,jk)  = ztmp1(jl,jk)/(1._dp-vtmpc1*ztmp1(jl,jk))
        psupw(jl,jk) = pqm1(jl,jk)/pqsw(jl,jk)-1.0_dp
        pesw(jl,jk)  = ztmp1(jl,jk)*papm1(jl,jk)*rv/rd

        ztmp1(jl,jk) = tlucuaw(it0_look(jl,jk))/papm1(jl,jk)     !<  e_s,w*Rd/Rv /p
        ztmp1(jl,jk) = MIN(ztmp1(jl,jk),0.5_dp)
        pqsw0(jl,jk)  = ztmp1(jl,jk)/(1._dp-vtmpc1*ztmp1(jl,jk))

        !SF for later use in 5.:
        pqswp1(jl,jk) = tlucuaw(itm1p1_look(jl,jk))/papm1(jl,jk)
        pqswp1(jl,jk) = MIN(pqswp1(jl,jk),0.5_dp)
        pqswp1(jl,jk) = pqswp1(jl,jk)/(1._dp-vtmpc1*pqswp1(jl,jk))
        !SFend for later use in 5.

!SF ice:

        ztmp1(jl,jk) = tlucua(itm1_look(jl,jk))/papm1(jl,jk)       !<  e_s,i*Rd/Rv /p
        ztmp1(jl,jk) = MIN(ztmp1(jl,jk),0.5_dp)
        pqsi(jl,jk)  = ztmp1(jl,jk)/(1._dp-vtmpc1*ztmp1(jl,jk))
        psupi(jl,jk) = pqm1(jl,jk)/pqsi(jl,jk)-1.0_dp
        pesi(jl,jk)  = ztmp1(jl,jk)*papm1(jl,jk)*rv/rd

        !SF for later use in 5.:
        pqsip1(jl,jk) = tlucua(itm1p1_look(jl,jk))/papm1(jl,jk)
        pqsip1(jl,jk) = MIN(pqsip1(jl,jk),0.5_dp)
        pqsip1(jl,jk) = pqsip1(jl,jk)/(1._dp-vtmpc1*pqsip1(jl,jk))
        !SFend for later use in 5.
 
     END DO !jl
  END DO !jk
 
END SUBROUTINE update_saturation_values_2d

SUBROUTINE sat_spec_hum_1d( &
              !-- IN
              kbdim, kproma, pap, ptlucu, &
              !-- OUT
              pes, pcor, pq)

  !-- Subroutine arguments
  INTEGER, INTENT(in) :: kbdim, kproma

  REAL(dp), INTENT(in) :: pap(kbdim)    !< Pressure at full levels [Pa]
  REAL(dp), INTENT(in) :: ptlucu(kbdim) !< Lookup table (e_s*Rd/Rv)

  REAL(dp), INTENT(out) :: pes(kbdim)  !< Saturation vapor pressure 
  REAL(dp), INTENT(out) :: pcor(kbdim)
  REAL(dp), INTENT(out) :: pq(kbdim)   !< Saturation specific humidity [kg/kg]

  pes(1:kproma) = ptlucu(1:kproma) / pap(1:kproma)
  pes(1:kproma) = MIN(pes(1:kproma), 0.5_dp)

  pcor(1:kproma) = 1._dp / (1._dp-vtmpc1*pes(1:kproma))

  pq(1:kproma) = pes(1:kproma) * pcor(1:kproma)

END SUBROUTINE sat_spec_hum_1d

SUBROUTINE sat_spec_hum_2d( &
              !-- IN
              kbdim, kproma, klev, pap, ptlucu, &
              !-- OUT
              pes, pcor, pq)

  !-- Subroutine arguments
  INTEGER, INTENT(in) :: kbdim, kproma, klev

  REAL(dp), INTENT(in) :: pap(kbdim,klev)    !< Pressure at full levels [Pa]
  REAL(dp), INTENT(in) :: ptlucu(kbdim,klev) !< Lookup table (e_s*Rd/Rv)

  REAL(dp), INTENT(out) :: pes(kbdim,klev)  !< Saturation vapor pressure
  REAL(dp), INTENT(out) :: pcor(kbdim,klev)
  REAL(dp), INTENT(out) :: pq(kbdim,klev)   !< Saturation specific humidity [kg/kg]

  pes(1:kproma,:) = ptlucu(1:kproma,:) / pap(1:kproma,:)
  pes(1:kproma,:) = MIN(pes(1:kproma,:), 0.5_dp)

  pcor(1:kproma,:) = 1._dp / (1._dp-vtmpc1*pes(1:kproma,:))

  pq(1:kproma,:) = pes(1:kproma,:) * pcor(1:kproma,:)

END SUBROUTINE sat_spec_hum_2d

SUBROUTINE freezing_below_238K( &
              !-- IN
              kbdim, kproma, ptmst, ptmst_rcp, &
              ld_frz_below_238K, paclc,        & !DN #475 (min cdnc)
              pxl, pcdnc,            &
              !-- OUT
              pfrl, pfrln)

  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL,INTENT(in) :: ld_frz_below_238K(kbdim) !< physical condition for freezing below 238K to occur

  REAL(dp), INTENT(in) :: ptmst, ptmst_rcp !< timestep
  REAL(dp), INTENT(in) :: paclc(kbdim)   !< Cloud cover
  REAL(dp), INTENT(in)  :: pxl(kbdim)  !< Cloud liquid water in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(in)  :: pcdnc(kbdim) !< Ice crystal number concentration (ICNC) [1/m3]

  REAL(dp), INTENT(out) :: pfrln(kbdim) !< cloud droplet freezing rate [m-3 s-1]
  REAL(dp), INTENT(out) :: pfrl(kbdim)  !< mass freezing rate [kg kg-1 s-1]

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim)

  pfrl(1:kproma)  = MERGE(pxl(1:kproma), 0._dp, ld_frz_below_238K(1:kproma))*ptmst_rcp

!--- Included for prognostic CDNC/IC scheme ----------------------------
  ztmp1(1:kproma) = MAX(pcdnc(1:kproma), 0._dp)
  ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, pxl(1:kproma) > cqtmin)
  pfrln(1:kproma) = MERGE(ztmp1(1:kproma)*ptmst_rcp, 0._dp, ld_frz_below_238K(1:kproma))

END SUBROUTINE freezing_below_238K

SUBROUTINE homogeneous_freezing( &
              !-- IN
              kbdim, kproma, ptmst, ptmst_rcp, &
              pt, paclc, pxl, pcdnc, prcm,     &
!>>DN bugfix
              ld_frz_below_238K,               &
!<<DN bugfix
              !-- OUT
              pfrl, pfrln)

  INTEGER, INTENT(in) :: kbdim, kproma

  REAL(dp), INTENT(in) :: ptmst, ptmst_rcp
  REAL(dp), DIMENSION(kbdim), INTENT(in) :: paclc, pxl, pcdnc, pt, prcm

!>>DN bugfix
  LOGICAL,INTENT(in) :: ld_frz_below_238K(kbdim) !< physical condition for freezing below 238K to occur
!<<DN bugfix

  REAL(dp), INTENT(out) :: pfrln(kbdim) !< cloud droplet freezing rate [m-3 s-1]
  REAL(dp), INTENT(out) :: pfrl(kbdim)  !< mass freezing rate [kg kg-1 s-1]

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim)
  REAL(dp) :: j_hom(kbdim), f_frz(kbdim)
  REAL(dp) :: t_c(kbdim)

  ! parameterization is given in celsius, so avoid confusion by new temperature variable
  t_c(1:kproma) = pt(1:kproma) - tmelt

  ! compute homogeneous freezing rate Jeffrey und Austin (1997)
  ztmp1(1:kproma) = 1.0e6_dp * 10._dp**(-7.63_dp-2.996_dp*(t_c(1:kproma)+30._dp))      ! T > 30 deg C
  ztmp2(1:kproma) = 1.0e6_dp * 10.**(-243.4_dp &                                       ! T < 30 deg C
                                    - 14.75_dp*t_c(1:kproma) &
                                    - 0.307_dp*t_c(1:kproma)**2 &
                                    - 0.00287_dp*t_c(1:kproma)**3 &
                                    - 0.0000102*t_c(1:kproma)**4)
  j_hom(1:kproma) = MERGE(ztmp1(1:kproma), ztmp2(1:kproma), t_c(1:kproma) > -30._dp)

  ! analytically integrated frozen fraction (~ probability of freezing * volume of droplet)
  f_frz(1:kproma) = 1._dp - EXP(-j_hom(1:kproma)*pi/6._dp*prcm(1:kproma)**3*ptmst)

  ! convert back to mass and number mixing ratio rates
  pfrl(1:kproma) = f_frz(1:kproma)*pxl(1:kproma)*ptmst_rcp
  pfrln(1:kproma) = f_frz(1:kproma)*pcdnc(1:kproma)*ptmst_rcp

!>>DN bugfix
  pfrl(1:kproma)  = MERGE(pfrl(1:kproma), 0._dp, ld_frz_below_238K(1:kproma))
  pfrln(1:kproma) = MERGE(pfrln(1:kproma), 0._dp, ld_frz_below_238K(1:kproma))
!<<DN bugfix

END SUBROUTINE homogeneous_freezing

SUBROUTINE het_mxphase_freezing( &
              !-- IN
              kbdim, kproma, ptmst, ptmst_rcp, &
              ld_mxphase_frz, papp1, ptkem1, pvervel, paclc, &
              pfracbcsol, pfracbcinsol, pfracdusol, pfracduai, pfracduci, &
              prho, prho_rcp, &
              prwetki, prwetai, prwetci, ptp1tmp, &
              picnc, pcdnc, pxib, pxlb, &
              !-- OUT
              pfrl, pfrln )

  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL, INTENT(in)  :: ld_mxphase_frz(kbdim) !< Physical conditions for het. mixed phase freezing to occur

  REAL(dp), INTENT(in) :: ptmst, ptmst_rcp    !< timestep
  REAL(dp), INTENT(in) :: papp1(kbdim)        !< pressure at full levels (t-1)
  REAL(dp), INTENT(in) :: ptkem1(kbdim)       !< turbulent kinetic energy (t-1)
  REAL(dp), INTENT(in) :: pvervel(kbdim)      !< large scale vertical velocity [m s-1]
  REAL(dp), INTENT(in) :: paclc(kbdim)        !< Cloud cover
  REAL(dp), INTENT(in) :: pfracbcsol(kbdim)   !< Fraction of BC in all soluble mixed modes
  REAL(dp), INTENT(in) :: pfracbcinsol(kbdim) !< Fraction of BC in all insoluble modes
  REAL(dp), INTENT(in) :: pfracdusol(kbdim)   !< Fraction of dust aerosols in all soluble mixed modes
  REAL(dp), INTENT(in) :: pfracduai(kbdim)    !< Fraction of dust in the insoluble accumulation mode
  REAL(dp), INTENT(in) :: pfracduci(kbdim)    !< Fraction of dust in the insoluble coarse mode
  REAL(dp), INTENT(in) :: prho(kbdim)         !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: prho_rcp(kbdim)     !< Inverse air density
  REAL(dp), INTENT(in) :: prwetki(kbdim)      !< wet radius, Aitken insoluble mode [m]
  REAL(dp), INTENT(in) :: prwetai(kbdim)      !< wet radius, accumulation insoluble mode [m]
  REAL(dp), INTENT(in) :: prwetci(kbdim)      !< wet radius, coarse insoluble mode [m]
  REAL(dp), INTENT(in) :: ptp1tmp(kbdim)      !< temporary value of the updated temperature (t) [K]
  REAL(dp), INTENT(in) :: picnc(kbdim)        !< Ice crystal number concentration (ICNC) [1/m3]
  REAL(dp), INTENT(in) :: pcdnc(kbdim)        !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp), INTENT(in) :: pxib(kbdim)         !< Cloud ice in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(in) :: pxlb(kbdim)         !< Cloud liquid water in the cloudy part of the grid box [kg/kg]

  REAL(dp), INTENT(out) :: pfrl(kbdim)        !< freezing rate [kg/kg]
  REAL(dp), INTENT(out) :: pfrln(kbdim)       !< Freezing rate for number conc. [1/m3]
  !local vars:
  INTEGER :: jl

  REAL(dp) :: zdfarbcki(kbdim) !< Aerosol diffusivity due to Brownian motion for insoluble BC
  REAL(dp) :: zdfarduai(kbdim) !< Aerosol diffusivity due to Brownian motion for insoluble accum. mode dust
  REAL(dp) :: zdfarduci(kbdim) !< Aerosol diffusivity due to Brownian motion for insoluble coarse mode dust
  REAL(dp) :: zetaair(kbdim)   !< Dynamic viscocity [kg/m/s]
  REAL(dp) :: zf1              !< temporary variable
  REAL(dp) :: zfrzcnt 
  REAL(dp) :: zfrzcntbc
  REAL(dp) :: zfrzcntdu
  REAL(dp) :: zfrzimm
  REAL(dp) :: znaimmbc
  REAL(dp) :: znaimmdu
  REAL(dp) :: zradl            !< mean volume radius of cloud droplets [m]
  REAL(dp) :: ztte             !< local temperature tendency
  REAL(dp) :: zomega           !< vertical velocity in the p-system [Pa/s]

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim)
  LOGICAL  :: ll1(kbdim), ll2(kbdim), ll3(kbdim)

 !--------- new freezing parameterisations (Lohmann & Diehl, 2006)
    ! corinna: included for contact/immersion freezing by dust and soot
 
  ztmp1(1:kproma) = 1._dp + 1.26_dp * 6.6E-8_dp / (prwetki(1:kproma) + eps) &
                                    * (p0sl_bg / papp1(1:kproma)) & !SF ccbcki
                                    * (ptp1tmp(1:kproma) / tmelt)
  ztmp2(1:kproma) = 1._dp + 1.26_dp * 6.6E-8_dp / (prwetai(1:kproma)+eps) &
                                    * (p0sl_bg / papp1(1:kproma)) & !SF ccduai
                                    * (ptp1tmp(1:kproma) / tmelt)
  ztmp3(1:kproma) = 1._dp + 1.26_dp * 6.6E-8_dp / (prwetci(1:kproma)+eps) &
                                    * (p0sl_bg / papp1(1:kproma)) & !SF ccduci
                                    * (ptp1tmp(1:kproma) / tmelt)
 
  zetaair(1:kproma) = 1.e-5_dp  &
                    * ( 1.718_dp + 0.0049_dp*(ptp1tmp(1:kproma)-tmelt) &
                      - 1.2e-5_dp*(ptp1tmp(1:kproma)-tmelt)*(ptp1tmp(1:kproma)-tmelt) )
 
  ll1(1:kproma) = (prwetki(1:kproma) < eps)

  zdfarbcki(1:kproma) = ak * ptp1tmp(1:kproma) * ztmp1(1:kproma) &
                      / ( 6._dp*pi*zetaair(1:kproma)*(prwetki(1:kproma)+eps) )
  zdfarbcki(1:kproma) = MERGE(0._dp, zdfarbcki(1:kproma), ll1(1:kproma))
 
  ll2(1:kproma) = (prwetai(1:kproma) < eps)

  zdfarduai(1:kproma) = ak * ptp1tmp(1:kproma) * ztmp2(1:kproma) &
                      / ( 6._dp*pi*zetaair(1:kproma)*(prwetai(1:kproma)+eps) )
  zdfarduai(1:kproma) = MERGE(0._dp, zdfarduai(1:kproma), ll2(1:kproma))
 
  ll3(1:kproma) = (prwetci(1:kproma) < eps)

! aerosol diffusivity
  zdfarduci(1:kproma) = ak * ptp1tmp(1:kproma) * ztmp3(1:kproma) &
                      / ( 6._dp*pi*zetaair(1:kproma)*(prwetci(1:kproma)+eps) )
  zdfarduci(1:kproma) = MERGE(0._dp, zdfarduci(1:kproma), ll3(1:kproma))
 
  ztmp1(1:kproma) = 0._dp
  ztmp2(1:kproma) = 0._dp
  DO jl=1,kproma
     IF(pxlb(jl) < cqtmin .OR. pcdnc(jl) < 1._dp) CYCLE
     ! volume mean water droplet radius
     zradl = (0.75_dp* pxlb(jl) * prho(jl) &
           / (pi*rhoh2o*pcdnc(jl)))**(1._dp/3._dp)

     zf1 = 4._dp*pi*zradl*pcdnc(jl)*prho_rcp(jl)

     zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1014_dp*(ptp1tmp(jl)-tmelt)+0.3277_dp)))  ! montmorillonite
     !zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1007_dp*(ptp1tmp(jl)-tmelt)+0.6935_dp)))  ! kaolinite

     !zfrzcntbc = MIN(1._dp,MAX(0._dp,-(0.0614_dp*(ptp1tmp(jl)-tmelt)+0.5730_dp)))
     zfrzcntbc = 0._dp ! disable BC contact freezing

     zfrzcnt = pxlb(jl) / pcdnc(jl) * prho(jl) * zf1 &
             * ( zfrzcntdu * ( zdfarduai(jl)*pfracduai(jl) &
                             + zdfarduci(jl)*pfracduci(jl) ) &
               + zfrzcntbc *   zdfarbcki(jl)*pfracbcinsol(jl) ) &
             * ( pcdnc(jl)*prho(jl)+picnc(jl)*prho(jl) )

     zfrzcnt = pxlb(jl)*(1._dp-EXP(-zfrzcnt/MAX(pxlb(jl), cqtmin)*ptmst))

     znaimmdu = 32.3_dp*pfracdusol(jl)     ! montmorillonite 
     !znaimmdu  = 6.15E-2_dp*pfracdusol((jl))  ! kaolinite 

     znaimmbc = 2.91E-3_dp*pfracbcsol(jl)

     !>>SF ToDo: create a function for computing zomega (see vervx calc in cloud_micro_interface)
     zomega = pvervel(jl) - fact_tke*SQRT(ptkem1(jl))*prho(jl)*grav !SF #345 changed TKE prefactor
     !<<SF ToDo
     ztte   = zomega / cpd *prho_rcp(jl)

     zfrzimm = -(znaimmdu+znaimmbc)*prho(jl)/rhoh2o*EXP(tmelt-ptp1tmp(jl))*MIN(ztte,0._dp) 
     zfrzimm = pxlb(jl)*(1._dp-EXP(-zfrzimm*pxlb(jl)/pcdnc(jl)*ptmst))

     ztmp1(jl) = zfrzcnt + zfrzimm

     ztmp1(jl) = MAX(0.0_dp,MIN(ztmp1(jl),pxlb(jl))) !SF pfrl surrogate

     ztmp2(jl) = pcdnc(jl)*ztmp1(jl)/(pxlb(jl)+eps)

     ztmp2(jl) = MAX(ztmp2(jl), 0._dp)  !SF pfrln surrogate
 
  ENDDO !SF end loop jl
 
  pfrl(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_mxphase_frz(1:kproma))
 
  ztmp1(1:kproma) = MIN(ztmp2(1:kproma), pcdnc(1:kproma))
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0._dp)
  pfrln(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_mxphase_frz(1:kproma))
 
  ztmp1(1:kproma) = MIN(pfrl(1:kproma), pxlb(1:kproma))
  pfrl(1:kproma)  = MERGE(ztmp1(1:kproma), 0._dp, ld_mxphase_frz(1:kproma))

  pfrl(1:kproma) = ptmst_rcp*pfrl(1:kproma)
  pfrln(1:kproma) = ptmst_rcp*pfrln(1:kproma)
 
END SUBROUTINE het_mxphase_freezing

SUBROUTINE falling_ice(&
              !--IN
              kbdim, kproma, ptmst, ptmst_rcp, &
              pxi, picnc, pvtim, pvtin, &
              pxiflx, pniflx, &
              pdpg, prho, &
              !--INOUT
              pvtimfal, pvtinfal, &
              !--OUT
              pxifal, pnifal)

  
  INTEGER, INTENT(IN) :: kbdim, kproma
  REAL(dp), INTENT(IN):: ptmst, ptmst_rcp

  REAL(dp), DIMENSION(kbdim), INTENT(INOUT) :: pvtimfal, pvtinfal
              
  REAL(dp), DIMENSION(kbdim), INTENT(IN)  :: pxi, picnc, pvtim, pvtin, &
                                             pxiflx, pniflx, &
                                             pdpg, prho

  REAL(dp), DIMENSION(kbdim), INTENT(OUT) :: pxifal, pnifal

  REAL(dp), DIMENSION(kbdim) :: zkm, zkn      ! gridbox fraction of vertical displacement in-cloud
  REAL(dp), DIMENSION(kbdim) :: zkmfal, zknfal! gridbox fraction of vertical displacement falling
  REAL(dp), DIMENSION(kbdim) :: zxiv, zniv    ! virtual number and mass mixing ratios for 
                                              ! falling crystals
  REAL(dp), DIMENSION(kbdim) :: zmflx_icr     ! increase in mass flux
  REAL(dp), DIMENSION(kbdim) :: znflx_icr     ! increase in number flux
  REAL(dp), DIMENSION(kbdim) :: zmflx_lft     ! left over mass flux from above
  REAL(dp), DIMENSION(kbdim) :: znflx_lft     ! left over number flux from above
  LOGICAL, DIMENSION(kbdim)  :: ll1, ll2

  zkm(1:kproma) = ptmst*pvtim(1:kproma)*prho(1:kproma)/pdpg(1:kproma)
  zkn(1:kproma) = ptmst*pvtin(1:kproma)*prho(1:kproma)/pdpg(1:kproma)
  zkmfal(1:kproma) = ptmst*pvtimfal(1:kproma)*prho(1:kproma)/pdpg(1:kproma)
  zknfal(1:kproma) = ptmst*pvtinfal(1:kproma)*prho(1:kproma)/pdpg(1:kproma)

  ! direct
  ! zkm(1:kproma) = MAX(MIN(1._dp, zkm(1:kproma)), 0._dp)
  ! zkn(1:kproma) = MAX(MIN(1._dp, zkn(1:kproma)), 0._dp)
  ! zkmfal(1:kproma) = MAX(MIN(1._dp, zkmfal(1:kproma)), 0._dp)
  ! zknfal(1:kproma) = MAX(MIN(1._dp, zknfal(1:kproma)), 0._dp)

  ! semi-analytic
  zkm(1:kproma) = 1._dp - EXP(-zkm(1:kproma))
  zkmfal(1:kproma) = 1._dp - EXP(-zkmfal(1:kproma))
  zkn(1:kproma) = 1._dp - EXP(-zkn(1:kproma))
  zknfal(1:kproma) = 1._dp - EXP(-zknfal(1:kproma))

  zxiv(1:kproma) = ptmst*pxiflx(1:kproma)/pdpg(1:kproma)
  zniv(1:kproma) = ptmst*pniflx(1:kproma)/pdpg(1:kproma)

  pxifal(1:kproma) = - pxi(1:kproma)*zkm(1:kproma) &            ! sink due to sedimentation from level
                     + zxiv(1:kproma)*(1._dp-zkmfal(1:kproma))  ! source due to sedimentation from above

  pnifal(1:kproma) = - picnc(1:kproma)*zkn(1:kproma) &          ! sink due to sedimentation from level
                     + zniv(1:kproma)*(1._dp-zknfal(1:kproma))  ! source due to sedimentation from above

  ! mass increase due to newly formed ice flux from level
  zmflx_icr(1:kproma) = pxi(1:kproma)*zkm(1:kproma)
  znflx_icr(1:kproma) = picnc(1:kproma)*zkn(1:kproma)

  ! left mass flux from above
  zmflx_lft(1:kproma) = zxiv(1:kproma)*zkmfal(1:kproma)
  znflx_lft(1:kproma) = zniv(1:kproma)*zknfal(1:kproma)

  ! weighted new flux velocity
  pvtimfal(1:kproma) = (zmflx_lft(1:kproma)*pvtimfal(1:kproma) + zmflx_icr(1:kproma)*pvtim(1:kproma)) &
                       / MAX(zmflx_lft(1:kproma) + zmflx_icr(1:kproma), cqtmin)
  pvtinfal(1:kproma) = (znflx_lft(1:kproma)*pvtinfal(1:kproma) + znflx_icr(1:kproma)*pvtin(1:kproma)) &
                       / MAX(znflx_lft(1:kproma) + znflx_icr(1:kproma), cqtmin)

  pxifal(1:kproma) = pxifal(1:kproma)*ptmst_rcp
  pnifal(1:kproma) = pnifal(1:kproma)*ptmst_rcp


END SUBROUTINE falling_ice

SUBROUTINE precip_formation_cold( &
              !--IN
              kbdim, kproma, ptmst, ptmst_rcp, &
              prho, prho_rcp, ll_cc, paclc, paclci, &
              pt, pqsnow, pxi, pxl, picnc, pcdnc, &
              prim, &
              !--OUT
              psaci, psacin, psaut, psautn)

  INTEGER, INTENT(IN) :: kbdim, kproma
  REAL(dp), INTENT(IN):: ptmst, ptmst_rcp
  LOGICAL, INTENT(IN) :: ll_cc(kbdim)
  REAL(dp), DIMENSION(kbdim), INTENT(IN)  :: prho, prho_rcp, paclc, paclci, &
                                             pt, pqsnow, pxi, pxl, picnc, pcdnc, &
                                             prim

  REAL(dp), DIMENSION(kbdim), INTENT(OUT) :: psaci, psaut, psacin, psautn

  REAL(dp) :: zcolleffi(kbdim) !< Collision efficiency for aggregation
  REAL(dp) :: zqrho(kbdim)     !< inverse air density (m3 kg-1)
  REAL(dp) :: zrim(kbdim)

  REAL(dp), DIMENSION(kbdim) :: ztmp1, ztmp2, ztmp3
  LOGICAL, DIMENSION(kbdim)  :: ll1, ll2

  ll1(1:kproma) = ll_cc .AND. pxi(1:kproma) > cqtmin
  ll2(1:kproma) = ll1(1:kproma) .AND. (pqsnow(1:kproma) > cqtmin) 

  zqrho(1:kproma) = 1.3_dp*prho_rcp(1:kproma)

!------------------------------------------------------
! collision efficiency

  ztmp1(1:kproma)     = fact_coll_eff*(pt(1:kproma)-tmelt)
  ztmp1(1:kproma)     = EXP(ztmp1(1:kproma))
  zcolleffi(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll1(1:kproma))

!------------------------------------------------------
! autoconversion

  ztmp1(1:kproma) = 17.5_dp / crhoi * prho(1:kproma) * zqrho(1:kproma)**0.33_dp
  ! ztmp1(1:kproma) = prho(1:kproma)*700*zcolleffi(1:kproma)*0.1 / crhoi * zqrho(1:kproma)**0.33_dp

  zrim(1:kproma) = MAX(1e-9_dp, prim(1:kproma))
  ztmp1(1:kproma) = MERGE(-6._dp/ztmp1(1:kproma)*LOG10(zrim(1:kproma)/ccrsnow), 0._dp, ll1(1:kproma))
  ztmp1(1:kproma) = MAX(1._dp, ztmp1(1:kproma)*ccftau) ! set conversion timescale to 1 sec if 
                                                       ! crystals are larger than 100 um
  ztmp1(1:kproma) = 1._dp / ztmp1(1:kproma)

  ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll1(1:kproma))
  psaut(1:kproma) = pxi(1:kproma)*(1._dp - 1._dp/(1._dp+pxi(1:kproma)*ptmst*ztmp1(1:kproma)))
  psaut(1:kproma) = MAX(0._dp, psaut(1:kproma))

  !------------------------------------------------------
! accretion

  ztmp1(1:kproma) = cons4*MAX(cqtmin,pqsnow(1:kproma))**0.8125_dp
  ztmp1(1:kproma) = pi*cn0s*3.078_dp*ztmp1(1:kproma)*zqrho(1:kproma)**0.5_dp
  ztmp1(1:kproma) = -ptmst*ztmp1(1:kproma)*zcolleffi(1:kproma)
  ztmp1(1:kproma) = EXP(ztmp1(1:kproma))
  ztmp1(1:kproma) = pxi(1:kproma) * (1._dp-ztmp1(1:kproma))

  psaci(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma))
  psaci(1:kproma) = MAX(0._dp, psaci(1:kproma))

!------------------------------------------------------
! assert ice mass conservation
  ztmp1(1:kproma) = psaci(1:kproma) + psaut(1:kproma)
  ztmp2(1:kproma) = pxi(1:kproma)*paclc(1:kproma)
  ztmp3(1:kproma) = conservation_reduction(kbdim, kproma, ztmp2(1:kproma), ztmp1(1:kproma))

  psaci(1:kproma) = psaci(1:kproma)*ztmp3(1:kproma)
  psaut(1:kproma) = psaut(1:kproma)*ztmp3(1:kproma)
  

!------------------------------------------------------
! associated number change

  ztmp1(1:kproma) = MERGE(picnc(1:kproma)/pxi(1:kproma), 0._dp, pxi(1:kproma) > 0._dp)

  psacin(1:kproma) = psaci(1:kproma)*ztmp1(1:kproma)
  psautn(1:kproma) = psaut(1:kproma)*ztmp1(1:kproma)

!------------------------------------------------------
! convert to rates

  psaut(1:kproma) = psaut(1:kproma)*ptmst_rcp
  psautn(1:kproma) = psautn(1:kproma)*ptmst_rcp

  psaci(1:kproma) = psaci(1:kproma)*ptmst_rcp
  psacin(1:kproma) = psacin(1:kproma)*ptmst_rcp

END SUBROUTINE precip_formation_cold

SUBROUTINE snow_cloud_coll( &
              !--IN
              kbdim, kproma, ptmst, ptmst_rcp, &
              prho, prho_rcp, ll_cc, paclc, paclci, &
              pt, pqsnow, pxi, pxl, picnc, pcdnc, &
              pviscos, &
              !--OUT
              psacl, psacln)

  INTEGER, INTENT(IN) :: kbdim, kproma
  REAL(dp), INTENT(IN):: ptmst, ptmst_rcp
  LOGICAL, INTENT(IN) :: ll_cc(kbdim)
  REAL(dp), DIMENSION(kbdim), INTENT(IN)  :: prho, prho_rcp, paclc, paclci, &
                                             pt, pqsnow, pxi, pxl, picnc, pcdnc, &
                                             pviscos

  REAL(dp), DIMENSION(kbdim), INTENT(OUT) :: psacl, psacln

  REAL(dp) :: zqrho(kbdim)     !< inverse air density (m3 kg-1)
  REAL(dp) :: zdplanar(kbdim)  !< constant snow volume radius
  REAL(dp) :: zstokes(kbdim)   !< stokes number for droplet and snow flake
  REAL(dp) :: zstcrit(kbdim)   !< Critical Stokes number
  REAL(dp) :: zrey(kbdim)      !< reynolds number for droplet and snow flake
  REAL(dp) :: zusnow(kbdim)    !< snow fall velocity
  REAL(dp) :: zudrop(kbdim)    !< droplet fall velocity
  REAL(dp) :: zcsacl(kbdim)    !< Temporary variable needed for accretion of snow with cloud droplets

  REAL(dp), DIMENSION(kbdim) :: ztmp1, ztmp2, ztmp3
  LOGICAL, DIMENSION(kbdim)  :: ll1, ll2, ll3, ll4, ll5, ll6, ll7, ll8

  ll2(1:kproma) = ll_cc(1:kproma)               .AND. &
                  (pqsnow(1:kproma)  >  cqtmin) .AND. &
                  (pxl(1:kproma)  >  cqtmin) .AND. &
                  (pcdnc(1:kproma) > cqtmin)

  zqrho(1:kproma) = 1.3_dp*prho_rcp(1:kproma)

!------------------------------------------------------
! drop fall velocity

  ztmp1(1:kproma) = ( 6._dp*pirho_rcp*prho(1:kproma)*MAX(0._dp, pxl(1:kproma)) &
                    /MAX(pcdnc(1:kproma), 1._dp) )**(1._dp/3._dp)
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 1.e-6_dp) !SF zdw

  zudrop(1:kproma) = 1.19e4_dp*2500._dp*ztmp1(1:kproma)**2*(1.3_dp*prho_rcp(1:kproma))**0.35_dp

!------------------------------------------------------
! drop fall velocity

  zdplanar(1:kproma) = 447.e-6_dp !DN #454: constant 100 mum volume radius -> constant max dimension of 447 mum
                                  !         after Table 2.2a in Pruppacher&Klett 1997
  zusnow(1:kproma) = 2.34_dp * (100._dp*zdplanar(1:kproma))**0.3_dp &
                   * (1.3_dp*prho_rcp(1:kproma))**0.35_dp

!------------------------------------------------------
! stokes number

  zstokes(1:kproma) = 2._dp*rgrav*(zusnow(1:kproma)-zudrop(1:kproma))*zudrop(1:kproma)/zdplanar(1:kproma)
  zstokes(1:kproma) = MAX(zstokes(1:kproma), cqtmin)

!------------------------------------------------------
! reynolds number

  zrey(1:kproma) = prho(1:kproma)*zdplanar(1:kproma)*zusnow(1:kproma)/pviscos(1:kproma)
  zrey(1:kproma) = MAX(zrey(1:kproma),cqtmin)

!------------------------------------------------------
! flow regime dependant accretion rate

  ll3(1:kproma) = (zrey(1:kproma) <=  5._dp)
  ll4(1:kproma) = (zrey(1:kproma) >   5._dp) .AND. (zrey(1:kproma) <  40._dp)
  ll5(1:kproma) = (zrey(1:kproma) >= 40._dp)

  ztmp1(1:kproma)   = 5.52_dp*zrey(1:kproma)**(-1.12_dp)
  ztmp2(1:kproma)   = 1.53_dp*zrey(1:kproma)**(-0.325_dp)
  zstcrit(1:kproma) = 1._dp
  zstcrit(1:kproma) = MERGE(ztmp1(1:kproma), zstcrit(1:kproma), ll3(1:kproma))
  zstcrit(1:kproma) = MERGE(ztmp2(1:kproma), zstcrit(1:kproma), ll4(1:kproma))

  zcsacl(1:kproma) = 0.2_dp * ( LOG10(zstokes(1:kproma)) - LOG10(zstcrit(1:kproma)) - 2.236_dp )**2
  zcsacl(1:kproma) = MIN(zcsacl(1:kproma), 1._dp-cqtmin)
  zcsacl(1:kproma) = MAX(zcsacl(1:kproma), 0._dp)
  zcsacl(1:kproma) = SQRT(1._dp - zcsacl(1:kproma))

  ll6(1:kproma) = ll5(1:kproma) .AND. (zstokes(1:kproma) <= 0.06_dp)
  ll7(1:kproma) = ll5(1:kproma) .AND. (zstokes(1:kproma) >  0.06_dp) .AND. (zstokes(1:kproma) <= 0.25_dp)
  ll8(1:kproma) = ll5(1:kproma) .AND. (zstokes(1:kproma) >  0.25_dp) .AND. (zstokes(1:kproma) <= 1.00_dp)

  ztmp1(1:kproma)  = (zstokes(1:kproma)+1.1_dp)**2/(zstokes(1:kproma)+1.6_dp)**2
  zcsacl(1:kproma) = MERGE(ztmp1(1:kproma), zcsacl(1:kproma), ll5(1:kproma))

  ztmp1(1:kproma)  = 1.034_dp*zstokes(1:kproma)**1.085_dp
  zcsacl(1:kproma) = MERGE(ztmp1(1:kproma), zcsacl(1:kproma), ll6(1:kproma))

  ztmp1(1:kproma)  = 0.787_dp*zstokes(1:kproma)**0.988_dp
  zcsacl(1:kproma) = MERGE(ztmp1(1:kproma), zcsacl(1:kproma), ll7(1:kproma))

  ztmp1(1:kproma)  = 0.7475_dp*LOG10(zstokes(1:kproma))+0.65_dp
  zcsacl(1:kproma) = MERGE(ztmp1(1:kproma), zcsacl(1:kproma), ll8(1:kproma))

  zcsacl(1:kproma) = MAX(MIN(zcsacl(1:kproma), 1._dp), 0.01_dp)
  zcsacl(1:kproma) = MERGE(zcsacl(1:kproma), 0._dp, ll2(1:kproma))

  ztmp1(1:kproma) = cons4*MAX(cqtmin,pqsnow(1:kproma))**0.8125_dp !SF zlamsm
  ztmp1(1:kproma) = pi*cn0s*3.078_dp*ztmp1(1:kproma)*zqrho(1:kproma)**0.5_dp
  ztmp2(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma)) 

  ztmp1(1:kproma) = -ptmst*ztmp2(1:kproma)*zcsacl(1:kproma)
  ztmp1(1:kproma) = EXP(ztmp1(1:kproma))
  ztmp1(1:kproma) = pxl(1:kproma)*(1._dp-ztmp1(1:kproma))

  ztmp2(1:kproma) = paclc(1:kproma)*ztmp1(1:kproma)
  psacl(1:kproma)   = MERGE(ztmp2(1:kproma), 0._dp, ll2(1:kproma))

!------------------------------------------------------
! associated number change

  ztmp1(1:kproma) = MERGE(pcdnc(1:kproma)/pxl(1:kproma), 0._dp, pxl(1:kproma) > 0._dp)
  psacln(1:kproma) = psacl(1:kproma)*ztmp1(1:kproma)

!------------------------------------------------------
! convert to rates

  psacl(1:kproma) = psacl(1:kproma)*ptmst_rcp
  psacln(1:kproma) = psacln(1:kproma)*ptmst_rcp 

END SUBROUTINE snow_cloud_coll

SUBROUTINE precip_formation_warm_2( &
              !--IN
              kbdim, kproma, ptmst, ptmst_rcp, &
              paclc_pre, prho, prho_rcp,       &
              pcdnc, pxl, pqr, prcm,           &
              !--OUT
              prpr, prprn, paut, pacc)

  !>>DN
  USE mo_param_switches,       ONLY: ac_scale_autoconversion, ac_scale_KK_LWP_exponent,&
       ac_scale_accretion, lMMPPE
  !<<DN
  
  INTEGER, INTENT(in) :: kbdim, kproma
  
  REAL(dp), INTENT(in) :: ptmst, ptmst_rcp
  REAL(dp), DIMENSION(kbdim), INTENT(IN) :: paclc_pre, prho, prho_rcp, pcdnc, pxl, pqr, prcm

  REAL(dp), DIMENSION(kbdim), INTENT(OUT):: prpr, prprn, paut, pacc

  REAL(dp), DIMENSION(kbdim) :: zcdnc_cm3      ! cdnc in 1/cm3
  REAL(dp), DIMENSION(kbdim) :: zcmass         ! mass of a mass-mean droplet
  LOGICAL, DIMENSION(kbdim)  :: ll1

! numbers in 1/cm3
  zcdnc_cm3(1:kproma) = MAX(1._dp, pcdnc(1:kproma))*prho(1:kproma)*1e-6_dp

! Autoconversion/Accretion rate from Khairoutdinov and Kogan, 2000
!>>DN  
!  paut(1:kproma) = ccraut*1350*MAX(0._dp,pxl(1:kproma))**(2.47)*zcdnc_cm3(1:kproma)**(-1.79)
  IF (.NOT.lMMPPE) THEN
     paut(1:kproma) = ccraut*1350*MAX(0._dp,pxl(1:kproma))**(2.47)*zcdnc_cm3(1:kproma)**(-1.79)
  ELSE
     paut(1:kproma) = ccraut*1350*MAX(0._dp,pxl(1:kproma))**(ac_scale_KK_LWP_exponent)*&
          zcdnc_cm3(1:kproma)**(ac_scale_autoconversion)     
  END IF
!<<DN

!>>DN  
!  pacc(1:kproma) = 3.7*pxl(1:kproma)*pqr(1:kproma)
  IF (.NOT.lMMPPE) THEN
     pacc(1:kproma) = 3.7*pxl(1:kproma)*pqr(1:kproma)
  ELSE
     pacc(1:kproma) = ac_scale_accretion*3.7*pxl(1:kproma)*pqr(1:kproma)     
  END IF
!<<DN

  prpr(1:kproma) = paut(1:kproma) + pacc(1:kproma)

  ll1(1:kproma) = pxl(1:kproma) > cqtmin .AND. pcdnc(1:kproma) > 1._dp

! ! mass of a mean droplet
!   zcmass(1:kproma) = 4._dp/3._dp*pi*rhoh2o*prcm(1:kproma)**3

! ! calculate number from mass weighted mean radius
!   prprn(1:kproma) = prpr(1:kproma)/MAX(cqtmin, zcmass(1:kproma))
!   prprn(1:kproma) = MERGE(prprn(1:kproma), pcdnc(1:kproma)*ptmst_rcp, ll1(1:kproma))

! rain out any water below cqtmin
  prprn(1:kproma) = MERGE(prpr(1:kproma)/pxl(1:kproma)*pcdnc(1:kproma), pcdnc(1:kproma)*ptmst_rcp, ll1(1:kproma))
  prpr(1:kproma) = MERGE(prpr(1:kproma), pxl(1:kproma)*ptmst_rcp, ll1(1:kproma))

END SUBROUTINE precip_formation_warm_2

SUBROUTINE precip_formation_warm( &
              !--IN
              kbdim, kproma, ptmst, ptmst_rcp, &
              ld_prcp_warm, pauloc, paclc, pclcstar, prho, prho_rcp, pxrp1, &
              pcdnc, pxlb, &
              !-- OUT
              prpr, prprn, paut, pacc)

!DN --> SF ToDo: divide into two subroutines: one for Khairoutdinov and Kogan, 2000; and one for Beheng (1994)
!DN --> SF ToDo: autoconversion and accretion cloud be separated as well

  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL, INTENT(in) :: ld_prcp_warm(kbdim) !< physical condition for warm precip to form

  REAL(dp), INTENT(in) :: pauloc(kbdim)   !< Part of the grid box allowed to particiation in accretion
                                          !< with newly formed condensate
  REAL(dp), INTENT(in) :: ptmst, ptmst_rcp!< timestep
  REAL(dp), INTENT(in) :: paclc(kbdim)    !< Cloud cover
  REAL(dp), INTENT(in) :: pclcstar(kbdim) !< Minimum of cloud cover and precipitation cover
  REAL(dp), INTENT(in) :: prho(kbdim)     !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: prho_rcp(kbdim) !< Inverse air density
  REAL(dp), INTENT(in) :: pxrp1(kbdim)    !< Rain mixing ratio (t) [kg/kg]
  REAL(dp), INTENT(in) :: pcdnc(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp), INTENT(in) :: pxlb(kbdim)  !< Cloud liquid water in the cloudy part of the grid box [kg/kg]

  REAL(dp), INTENT(out) :: prpr(kbdim)     !< rain formation rate [kg/kg]
  REAL(dp), INTENT(out) :: paut(kbdim)     !< for diagnostics [kg/kg]
  REAL(dp), INTENT(out) :: pacc(kbdim)     !< for diagnostics [kg/kg]
  REAL(dp), INTENT(out) :: prprn(kbdim)    !< Rain formation rate for number conc. [1/m3]

  !local vars:
  REAL(dp) :: zrac1(kbdim)     !< Accretion of cloud water with rain from above [kg/kg]
  REAL(dp) :: zrac2(kbdim)     !< Accretion of cloud water with rain formed inside the grid box [kg/kg]
  REAL(dp) :: zraut(kbdim)     !< Autoconversion of cloud droplets [kg/kg]
  REAL(dp) :: zrautself(kbdim) !< Sum of autoconversion and self collection [kg/kg]
  REAL(dp) :: zxlb(kbdim)      !< local copy of pxlb. RD Hack..
  REAL(dp) :: zcdnc(kbdim)     !< local copy of pcdnc. RD Hack..
  REAL(dp) :: zxlb0(kbdim)      !< local copy of pxlb. RD Hack..
  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim)
  LOGICAL  :: ll1(kbdim)

! Autoconversion rate from Khairoutdinov and Kogan, 2000

  zxlb(1:kproma) = pxlb(1:kproma)
  zxlb0(1:kproma) = pxlb(1:kproma)
  zcdnc(1:kproma) = MAX(1._dp, pcdnc(1:kproma))
  ztmp1(1:kproma) = ccraut*1350._dp*(1.e-6_dp*zcdnc(1:kproma)*prho(1:kproma))**(-1.79_dp)

  ztmp1(1:kproma) = zxlb(1:kproma) * (  1._dp &
                                     - (1._dp + ptmst*exm1_1*ztmp1(1:kproma)*&
                                                    MAX(cqtmin, zxlb(1:kproma))**exm1_1)**exp_1)

  ztmp1(1:kproma) = MIN(zxlb(1:kproma), ztmp1(1:kproma))
  zraut(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))

  ztmp1(1:kproma) = zxlb(1:kproma) - zraut(1:kproma)
  ztmp2(1:kproma) = zxlb(1:kproma) !SF keeps zxlb for later use
  zxlb(1:kproma)  = MERGE(ztmp1(1:kproma), zxlb(1:kproma), ld_prcp_warm(1:kproma))

!--- zrac1 is formed by accretion with rain from above
!--- zrac2 is formed by accretion with newly formed rain inside the grid box

  ztmp1(1:kproma) = -3.7_dp*ptmst*pxrp1(1:kproma)
  ztmp1(1:kproma) = zxlb(1:kproma)*(1._dp-EXP(ztmp1(1:kproma)))
  zrac1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))

  zxlb(1:kproma) = zxlb(1:kproma) - zrac1(1:kproma)

  ztmp1(1:kproma) = -3.7_dp*ptmst*pauloc(1:kproma)*prho(1:kproma)*zraut(1:kproma)
  ztmp1(1:kproma) = zxlb(1:kproma)*(1._dp-EXP(ztmp1(1:kproma)))
  zrac2(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))

  zxlb(1:kproma) = zxlb(1:kproma) - zrac2(1:kproma)

  prpr(1:kproma) = paclc(1:kproma)    * (zraut(1:kproma)+zrac2(1:kproma)) &
                 + pclcstar(1:kproma) *  zrac1(1:kproma)
  prpr(1:kproma) = MAX(0._dp, prpr(1:kproma))

!--- Autoconversion also changes the number of cloud droplets (prprn)

  ztmp1(1:kproma) = (zraut(1:kproma)+zrac1(1:kproma)+zrac2(1:kproma)) &
                  / (zxlb0(1:kproma)+eps)
  ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))

  ztmp2(1:kproma) = zcdnc(1:kproma)*ztmp1(1:kproma)*prho_rcp(1:kproma)
  prprn(1:kproma) = MERGE(ztmp2(1:kproma), 0._dp, ld_prcp_warm(1:kproma))

!--- save accretion and autoconversion rates for diagnostics
  pacc(1:kproma) = (zrac1(1:kproma) + zrac2(1:kproma))*ptmst_rcp
  paut(1:kproma) = zraut(1:kproma)*ptmst_rcp

END SUBROUTINE precip_formation_warm

SUBROUTINE set_lookup_index_1d( &
              !-- IN
              kbdim, kproma, pt, message, &
              !-- OUT
              kindex)

  !-- Subroutine arguments
  INTEGER, INTENT(in)          :: kbdim, kproma
  REAL(dp), INTENT(in)         :: pt(kbdim) !< Temperature [K]
  CHARACTER(len=*), INTENT(in) :: message   !< String to identify the location where a lookuptable overflow occured

  INTEGER, INTENT(out) :: kindex(kbdim) !< Corresponding index in lookup table

  !-- Local vars
  REAL(dp) :: zt_1000(kbdim)
  LOGICAL  :: ll_look(kbdim)

  zt_1000(1:kproma) = 1000._dp*pt(1:kproma)
  kindex(1:kproma)  = NINT(zt_1000(1:kproma))

  ll_look(1:kproma) = (kindex(1:kproma)<jptlucu1 .OR. kindex(1:kproma)>jptlucu2)

  IF (ANY(ll_look(1:kproma)))  CALL lookuperror (message)

END SUBROUTINE set_lookup_index_1d

SUBROUTINE set_lookup_index_2d( &
              !-- IN
              kbdim, kproma, klev, pt, message, &
              !-- OUT
              kindex)

  !-- Subroutine arguments
  INTEGER, INTENT(in)          :: kbdim, kproma, klev
  REAL(dp), INTENT(in)         :: pt(kbdim,klev) !< Temperature [K]
  CHARACTER(len=*), INTENT(in) :: message        !< String to identify the location where a lookuptable overflow occured

  INTEGER, INTENT(out) :: kindex(kbdim,klev) !< Corresponding index in lookup table

  !-- Local vars
  REAL(dp) :: zt_1000(kbdim,klev)
  LOGICAL  :: ll_look(kbdim,klev)

  zt_1000(1:kproma,:) = 1000._dp*pt(1:kproma,:)
  kindex(1:kproma,:)  = NINT(zt_1000(1:kproma,:))

  ll_look(1:kproma,:) = (kindex(1:kproma,:)<jptlucu1 .OR. kindex(1:kproma,:)>jptlucu2)

  IF (ANY(ll_look(1:kproma,:))) CALL lookuperror (message) 

END SUBROUTINE set_lookup_index_2d

FUNCTION ice_switch_1d(kbdim, kproma, pximin, pt, pxi, pvervx, pvervmax) RESULT(ll_ice)

  INTEGER :: kbdim, kproma
  REAL(dp), DIMENSION(kbdim) :: pt, pxi, pvervx, pvervmax
  LOGICAL :: ll_ice(kbdim)
  REAL(dp) :: pximin

  LOGICAL :: ll_up(kbdim)

  ! updraft condition
  IF(lvdiff) THEN
     ll_up(1:kproma) = 0.01_dp*pvervx(1:kproma) < pvervmax(1:kproma)
  ELSE
     ll_up(1:kproma) = 0.01_dp*csubw < pvervmax(1:kproma)
  ENDIF
  ! allow to turn off this condition
  ll_up(1:kproma) = ll_up(1:kproma) .OR. ldisableupdraftcond
  
  ll_ice(1:kproma)  = (pt(1:kproma) < cthomi)       .OR. &
                      (pt(1:kproma) < tmelt       .AND. &
                        ll_up(1:kproma))

END FUNCTION ice_switch_1d

FUNCTION ice_switch_2d(kbdim, kproma, klev, pximin, pt, pxi, pvervx, pvervmax) RESULT(ll_ice)

  INTEGER :: kbdim, kproma, klev
  REAL(dp), DIMENSION(kbdim,klev) :: pt, pxi, pvervx, pvervmax
  LOGICAL :: ll_ice(kbdim,klev)
  REAL(dp) :: pximin

  LOGICAL :: ll_up(kbdim,klev)

  ! updraft condition
  IF(lvdiff) THEN
     ll_up(1:kproma,:) = 0.01_dp*pvervx(1:kproma,:) < pvervmax(1:kproma,:)
  ELSE
     ll_up(1:kproma,:) = 0.01_dp*csubw < pvervmax(1:kproma,:)
  ENDIF
  ! allow to turn off this condition
  ll_up(1:kproma,:) = ll_up(1:kproma,:) .OR. ldisableupdraftcond

  ll_ice(1:kproma,:)  = (pt(1:kproma,:) < cthomi)      .OR. &
                        (pt(1:kproma,:) < tmelt      .AND. &
                          ll_up(1:kproma,:))

END FUNCTION ice_switch_2d

FUNCTION threshold_vert_vel_1d(kbdim, kproma, pesw, pesi, picnc, prho, price, peta) RESULT(pvervmax)

  !-- SF: this function computes the threshold vertical velocity needed for the WBF criterion

  !-- Function arguments
  INTEGER :: kbdim, kproma

  REAL(dp) :: pesw(kbdim)     !< Saturation vapor pressure w.r.t. water [Pa]
  REAL(dp) :: pesi(kbdim)     !< Saturation vapor pressure w.r.t. ice [Pa]
  REAL(dp) :: picnc(kbdim)    !< Ice crystal number concentration (ICNC) [1/kg]
  REAL(dp) :: prho(kbdim)     !< air density [kg/m3]
  REAL(dp) :: price(kbdim)    !< Volume mean ice crystal radius [m]
  REAL(dp) :: peta(kbdim)     !< ??
  REAL(dp) :: pvervmax(kbdim) !< Threshold vertical velocity

  pvervmax(1:kproma) = (pesw(1:kproma) - pesi(1:kproma)) / pesi(1:kproma) &
                     * picnc(1:kproma) * prho(1:kproma) * price(1:kproma) * peta(1:kproma)

END FUNCTION threshold_vert_vel_1d

FUNCTION threshold_vert_vel_2d(kbdim, kproma, klev, pesw, pesi, picnc, prho, price, peta) RESULT(pvervmax)

  !-- SF: this function computes the threshold vertical velocity needed for the WBF criterion

  !-- Function arguments
  INTEGER :: kbdim, kproma, klev

  REAL(dp) :: pesw(kbdim,klev)     !< Saturation vapor pressure w.r.t. water [Pa]
  REAL(dp) :: pesi(kbdim,klev)     !< Saturation vapor pressure w.r.t. ice [Pa]
  REAL(dp) :: picnc(kbdim,klev)    !< Ice crystal number concentration (ICNC) [1/m3]
  REAL(dp) :: prho(kbdim,klev)     !< air density [kg/m3]
  REAL(dp) :: price(kbdim,klev)    !< Volume mean ice crystal radius [m]
  REAL(dp) :: peta(kbdim,klev)     !< ??
  REAL(dp) :: pvervmax(kbdim,klev) !< Threshold vertical velocity

  pvervmax(1:kproma,:) = (pesw(1:kproma,:) - pesi(1:kproma,:)) / pesi(1:kproma,:) &
                       * picnc(1:kproma,:) * prho(1:kproma,:) * price(1:kproma,:) * peta(1:kproma,:)

END FUNCTION threshold_vert_vel_2d

FUNCTION effective_2_volmean_radius_param_Schuman_2011_1d(kbdim, kproma, prieff) RESULT(prvolmean)

  ! Simple param of r/re approximating the Schumann et al. 2011 data:

  !-- Function arguments
  INTEGER :: kbdim, kproma

  REAL(dp) :: prieff(kbdim)    !< effective ice crystal radius [m] <-- beware of units!
  REAL(dp) :: prvolmean(kbdim) !< volume mean ice crystal radius [m]     <-- beware of units!

  prvolmean(:) = 0._dp
  prvolmean(1:kproma) = MAX(1.e-6_dp,conv_effr2mvr*prieff(1:kproma))

END FUNCTION effective_2_volmean_radius_param_Schuman_2011_1d

FUNCTION effective_2_volmean_radius_param_Schuman_2011_2d(kbdim, kproma, klev, prieff) RESULT(prvolmean)

  ! Simple param of r/re approximating the Schumann et al. 2011 data:

  !-- Function arguments
  INTEGER :: kbdim, kproma, klev

  REAL(dp) :: prieff(kbdim,klev)    !< effective ice crystal radius [m] <-- beware of units!
  REAL(dp) :: prvolmean(kbdim,klev) !< volume mean ice crystal radius [m]     <-- beware of units!

  prvolmean(:,:) = 0._dp
  prvolmean(1:kproma,:) = MAX(1.e-6_dp,conv_effr2mvr*prieff(1:kproma,:))

END FUNCTION effective_2_volmean_radius_param_Schuman_2011_2d

END MODULE mo_cloud_micro_p3
