!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! @brief cloud_micro_interface computes all stratiform microphysics
!!
!! This routine computes the tendencies of the four prognostic
!! variables (temperature t, specific humidity q, cloud liquid
!! water xl, cloud ice xi) due to phase changes (condensation/
!! deposition, evaporation/sublimation of rain/snow falling
!! into the unsaturated part of the grid box, melting of snow,
!! melting/freezing of cloud ice/cloud water, sedimentation of
!! cloud ice, and precipitation formation in warm, cold and
!! mixed phase clouds).
!! The precipitation at the surface (rain and snow) is used in
!! later for computing the land surface hydrology in *surf*.
!! The cloud parameters (cloud cover, cloud liquid water and
!! cloud ice) are used for the calculation of radiation at the
!! next timestep.
!! Attention: 
!! In the current version the advective tendencies of skewness 
!! and variance are set to zero.
!!
!! @par References
!!     Lohmann and Roeckner, 1996: Clim. Dyn. 557-572
!!     Levkov et al., 1992: Beitr. Phys. Atm. 35-58.          (ice phase)
!!     Beheng, 1994: Atmos. Res. 193-206.                    (warm phase)
!!     Lenderink et al., 1998; KNMI-REPORT NO. 98-13       (condensation)
!!     Tompkins 2002, J. Atmos. Sci.                        (cloud cover)
!!
!! @author 
!! <ol>
!! <li> U.Lohmann    MPI-Hamburg  1995
!! <li> G.Lenderink  KNMI, de Bilt 1998
!! <li> M.Esch       MPI-Hamburg  1999
!! </ol>
!! @par Revision History
!! <ol>
!!<li>E.Roeckner    MPI-Hamburg  2000
!!<li>A.Tompkins    MPI-Hamburg  2000
!!<li>U.Schlese     MPI-Hamburg  2003
!!<li>U.Lohmann     Dalhousie University 2002-2006: Prognostic CDNC/IC scheme
!!<li>P.Stier       MPI-Hamburg          2002-2006: Prognostic CDNC/IC scheme
!!<li>                                              Scavenging parameters
!!<li>J.Zhang       Dalhousie University      2004: Prognostic CDNC/IC scheme
!!<li>S. Ferrachat  ETH Zuerich  2008: complete re-writing to allow vectorization 
!!<li>                                 cleanup of the code 
!!<li>S. Ferrachat ETH Zuerich 2010-04: removed HAM dependencies
!!<li>S. Ferrachat ETH Zuerich 2015-06: process-splitting refactoring (calls to subroutines with strict intent)
!! </ol>
!!
!! @par This module is used by
!! physc
!!
!! @par Responsible coder
!! sylvaine.ferrachat@env.ethz.ch
!!
!! @par Copyright
!! 2009 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ECHAM is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!! violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!! copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!! an according license agreement with MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_cloud_micro_2m

  USE mo_kind,               ONLY : dp
  USE mo_math_constants,     ONLY : pi
  USE mo_physical_constants, ONLY : cpd, grav, rgrav, rd, alv, als, rv,   &
                                    vtmpc1, vtmpc2, rhoh2o, ak, tmelt, p0sl_bg
  USE mo_time_control,   ONLY : delta_time, time_step_len
  USE mo_param_switches, ONLY : nauto, &   !++mgs
                                nic_cirrus, &
#ifdef HAMMOZ
                                lorocirrus, &
#endif
                                lsecprod, lconv, nclmi_progn, &
                                !>>UP #782.1
                                ldetr_liquid, &
                                !<<UP #782.1
                                !>>UP #797
                                lsedfix, &
                                !<<UP #797
                                !>>UP emulator: switches and values for phasing
                                !CMP processes, see e.g. #757
                                lemuphase_nic_cirrus, eta_emu_nic_cirrus, &
                                lemuphase_riming, eta_emu_riming, &
                                lemuphase_wbf, eta_emu_wbf, &
                                lemuphase_icaggr, eta_emu_icaggr, &
                                lemuphase_icaccr, eta_emu_icaccr, &
                                lemuphase_cdaccr, eta_emu_cdaccr, &
                                lemuphase_cdautc, eta_emu_cdautc, &
                                lemuphase_sci, eta_emu_sci, &
                                lemuphase_dep, eta_emu_dep, &
                                lemuphase_icnucl, eta_emu_icnucl, &
                                lemuphase_sub_evp_falling, eta_emu_subfis, &
                                eta_emu_evpr, &
                                lemuphase_cdnuc, eta_emu_cdnuc, &
                                lemuphase_sip, eta_emu_sip, &
                                lemuphase_mlt, eta_emu_mlt, &
                                !<<UP
                                !>>UP #821
                                lcmpsimpl_prescr, &
                                ncmpsimpl_prescr_rime, &
                                ncmpsimpl_prescr_icnucl, &
                                ncmpsimpl_prescr_subfis, &
                                ncmpsimpl_prescr_sci, &
                                ncmpsimpl_prescr_icaccr, &
                                ncmpsimpl_prescr_mlt, &
                                !<<UP
                                lpr_corr, & !UP, #783.3
                                lnewdiags, & !UP: switch for new, additional
                                !diagnostics, see e.g. #783
                                !>>UP #782.2
                                ldetr_convicfactor, convicfactor, &
                                !<<UP #782.2
                                !>>MA #790
                                lemuphase_hetfrz, eta_emu_hetfrz, &
                                !<<MA #790
                                !>>MA #795
                                lemuphase_fvic, eta_emu_fvic, &
                                !<<MA #795
                                !>>MA switch WBF, #765
                                lctrl_WBF, &
                                !<<MA, #765
                                lslf !UP #839 
  USE mo_echam_cloud_params, ONLY : cqtmin, cvtfall, crhosno, cn0s, ccwmin      &
                                  , cthomi,  clmax, clmin, jbmin, jbmax, lonacc &
                                  , ccraut, ceffmin, ceffmax, crhoi, ccsaut
  USE mo_cloud_utils,    ONLY: epsec, xsec, qsec, eps, mi, &
                               ri_vol_mean_1, ri_vol_mean_2, &
                               alfased_1, alfased_2, alfased_3, &
                               betased_1, betased_2, betased_3, &
                               icemin, &
                               cdi, mw0, mi0, mi0_rcp, ka, kb, &
                               alpha, xmw, fall, rhoice, conv_effr2mvr, clc_min, icemax, &
                               dw0, exm1_1, exp_1, exm1_2, &
                               exp_2, pirho_rcp, cap, cons4, cons5, &
                               fact_PK, pow_PK, & !SF
                               get_util_var, get_cloud_bounds, fact_coll_eff, fact_tke
  USE mo_boundary_condition, ONLY: bc_find, bc_apply
  USE mo_echam_convect_tables, ONLY: lookuperror, jptlucu1    &
                                   , jptlucu2, tlucua, tlucub, tlucuaw
  !>>UP: adding additional timers for CMP processes
  USE mo_control,          ONLY: ltimer
  USE mo_cmp_timer,        ONLY: timer_start,    timer_stop,                   &
                                 timer_cmp_icenucl
  USE mo_hammoz_timer,     ONLY: timer_ham_ifdef 
  !UP: below is debugging tries
  USE mo_submodel,         ONLY: print_value
  !>>UP #821
  USE mo_cmp_simpl,        ONLY: ibc_cmpsimpl_rime, ibc_cmpsimpl_icnucl, ibc_cmpsimpl_icaccr, &
                                 ibc_cmpsimpl_subfs, ibc_cmpsimpl_subfi, ibc_cmpsimpl_sci,    &
                                 ibc_cmpsimpl_mlt
  USE mo_mpi,              ONLY: p_io, p_pe
  ! UP: old emulator tries at phasing between two cirrus schemes
  !USE mo_cmp_emulator,     ONLY: lemuphase_nic_cirrus, eta_emu_nic_cirrus
  !<<UP
!davidn
  USE mo_activ,          ONLY :&
!>>SF Kasja diags
                              qnucl, qcnd, qlwc_detr, qevp_lwc, qautn, qracl, &
                              qracln, qsacl, qsacln, qfrz, qfrzn, qnuci, qdep, &
                              qiwc_detr, qsub_iwc, qagg, qaggn, qsaci, qsacin, &
                              qselfn, qsecprod, qsecprodn, qsedi, qsedin, qmlt, qmltn, &
                              qevp_rain, qsub_snow,qrprn,qsprn, &
!<<SF Kasja diags              
!>>UP #783     
                              qfrznhet, qfrznhom, &
!<<UP #783
!>>DN: new diags
                              qxmlt,qevabfn,qevabf,qdepbf,qcdnc_detr,qicnc_detr,qgentl,qgenti,&
                              qcorric1,qcorric2,qcorric3,qcorric4,qcorric5,qcorric6,qcorric7,&
                              qcorrcd1,qcorrcd2,qcorrcd3,qcorrcd4,qcorrcd5,qcorrcd6,qcorrcd7,qcorrcd8,&
                              qcorrcd9,qcorrcd10,qcorrxi,qcorrxi2,qcorrxl,qcorrxl2,qxlte,qxite,&
                              qxttecdnc,qxtteicnc,qspr,qrpr,qmltn2,qtestCD,qtestIC,qtestLWC,qtestIWC,&
                              qsub_ice,qmlt_snow,qmlts_atm,qmlt_conv,&
!<<DN: new diags
!>>UP: new diags for activation correction terms
! and split correction terms further, #783
                              qcorrcd2unphys_2d,dcorrcd2unphys_2d,&
                              qcorrcd1_1, qcorrcd1_2, qcorrcd1_3, qcorrcd1_4,&
                              qcorrcd2_1, qcorrcd2_2, qcorric1_1,&                
                              dcorrcd1_1, dcorrcd1_2, dcorrcd1_3, dcorrcd1_4,&
                              dcorrcd2_1, dcorrcd2_2, dcorric1_1,&                
!UP #783.3
                              qcorrcd2_3, qcorric1_2, qcorric1_3, &
                              dcorrcd2_3, dcorric1_2, dcorric1_3, &
!UP new diags for 1CCN
                              diag_cnt_cqtmins1, diag_cnt_cqtminl1, &
                              diag_cnt_cqtmins2, diag_cnt_cqtminl2, &
                              diag_cnt_cqtmins3, diag_cnt_cqtminl3, &
                              diag_cnt_ll12d_1, diag_cnt_ll12d_2, &
!<<UP: new diags
!>>DN: burden
                              daut,dfre,dnucl, dcnd, dlwc_detr, devp_lwc, dautn, dracl, &
                              dracln, dsacl, dsacln, dfrz, dfrzn, dnuci, ddep, &
                              diwc_detr, dsub_iwc, dagg, daggn, dsaci, dsacin, &
                              dselfn, dsecprod, dsecprodn, dsedi, dsedin, dmlt, dmltn, &
                              devp_rain, dsub_snow,drprn,dsprn, &
                              dxmlt,devabfn,devabf,ddepbf,dcdnc_detr,dicnc_detr,dgentl,dgenti,&
                              dcorric1,dcorric2,dcorric3,dcorric4,dcorric5,dcorric6,dcorric7,&
                              dcorrcd1,dcorrcd2,dcorrcd3,dcorrcd4,dcorrcd5,dcorrcd6,dcorrcd7,dcorrcd8,&
                              dcorrcd9,dcorrcd10,dcorrxi,dcorrxi2,dcorrxl,dcorrxl2,dxlte,dxite,&
                              dxttecdnc,dxtteicnc,dspr,drpr,dmltn2,dtestCD,dtestIC,dtestLWC,dtestIWC,&
                              dsub_ice,dmlt_snow,dsedi_sfc,dmlts_atm,dmlts_sfc,dsnow,dcldtte,dmlt_conv, &
!<<DN: burden
!>>UP #783     
                              dfrznhet, dfrznhom, &
!<<UP #783
!>>UP new lifetime diagnostics
                              q_lwc_sources, q_lwc_sinks, q_iwc_sources, q_iwc_sinks, &
                              q_cdnc_sources, q_cdnc_sinks, q_icnc_sources, q_icnc_sinks, &
                              d_lwc_lifetime_sources, d_lwc_lifetime_sinks, &
                              d_iwc_lifetime_sources, d_iwc_lifetime_sinks, &
                              d_cdnc_lifetime_sources, d_cdnc_lifetime_sinks, &
                              d_icnc_lifetime_sources, d_icnc_lifetime_sinks

!davidn

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cloud_micro_interface, eff_ice_crystal_radius, breadth_factor

  INTEGER :: ibc_cvcbot, ibc_wcape, ibc_tconv, ibc_detr_cond !< indices for vars passed with the 
                                                             !< boundary condition scheme

  LOGICAL :: ll_het !SF now set by cloud_subm_1

  REAL(dp) ::  zcons1, zcons2, zcons3      , & !< various constants
               zdtime, zdtime_rcp, ztmst, ztmst_rcp,ztmstdt, zdt !< various delta t's and reciprocal 

  INTERFACE set_lookup_index
     MODULE PROCEDURE set_lookup_index_1d
     MODULE PROCEDURE set_lookup_index_2d
  END INTERFACE set_lookup_index

  INTERFACE sat_spec_hum
     MODULE PROCEDURE sat_spec_hum_1d
     MODULE PROCEDURE sat_spec_hum_2d
  END INTERFACE sat_spec_hum

  INTERFACE threshold_vert_vel
     MODULE PROCEDURE threshold_vert_vel_1d
     MODULE PROCEDURE threshold_vert_vel_2d
  END INTERFACE threshold_vert_vel

  INTERFACE effective_2_volmean_radius_param_Schuman_2011
     MODULE PROCEDURE effective_2_volmean_radius_param_Schuman_2011_1d
     MODULE PROCEDURE effective_2_volmean_radius_param_Schuman_2011_2d
  END INTERFACE effective_2_volmean_radius_param_Schuman_2011

  CONTAINS

SUBROUTINE cloud_micro_interface( &
              !-- IN:
              kproma, kbdim, klev, klevp1, ktrac, ktdia, krow, knvb, &
              paphm1, papm1, papp1, pqm1, ptm1, ptvm1, pxlm1, &
              pxim1, &
              pvervel, pgeo, pxtm1, paphp1, ptkem1, &
              !-- INOUT:
              pqtec, paclc, paclcac, &
              paclcov, paprl, pqvi, pxlvi, pxivi, pqte, ptte, &
              pxlte, pxite, pxtte, paprs, pacdnc, &
              !-- OUT:
              picnc, prelhum, pssfl, prsfl)

USE mo_vphysc,             ONLY : set_vphysc_var      !++mgs
USE mo_activ,              ONLY : swat,            &
                                  qnuc, qaut, qacc, qfre, qmel,                      &
                                  cdnc_acc, lwc_acc, cloud_time, cloud_cover_duplic, &
                                  cdnc_burden_acc, reffl_acc, burden_time,           &
                                  icnc_burden_acc, reffi_acc, icnc_acc,              &
                                  cliwc_time, burdic_time, iwc_acc,                  &
                                  cdnc_burden, icnc_burden, cdnc, icnc,              &
                                  sice, reffl_ct, reffl_time, cdnc_ct,               &
                                  reffi_tovs, reffi_time, iwp_tovs,                  &
                                  idt_cdnc, idt_icnc, nfrzmod, reffl, reffi,         &              
                                  !>>UP #844
                                  slfdiag_iwp_mxpT, slfdiag_lwp_mxpT,                &
                                  slfdiag_iwp_l35, slfdiag_lwp_g0,                   & 
                                  slfdiag_icnc_mxpT, slfdiag_cdnc_mxpT,              &
                                  slfdiag_icnc_l35, slfdiag_cdnc_g0,                 & 
                                  !<<UP #844
                                  cdnc_bins, lwc_bins, cdnc_lwc_binned,              &
                                  prcp_bins, prcp_binned, cdnc_binned
USE mo_conv,               ONLY : cdncact_cv,     &
                                  conv_time, twc_conv
USE mo_cirrus,             ONLY : xfrzmstr
USE mo_submodel_interface, ONLY: cloud_subm_1, cloud_subm_2
USE mo_memory_g2a,         ONLY: vm1, um1
USE mo_memory_g3b,         ONLY: orostd,oromea,orogam,orothe, &
                                 oropic,oroval,orosig
#ifdef HAMMOZ
USE mo_orocirrus,          ONLY: orocirrus_w, orocirrus_cc
#endif

USE mo_time_control,       ONLY: lstart, lresume
!>>SF
USE mo_debugs
!<<SF
!>>UP: new CMP timers and new diagnostics for investigating the freezing/ice
!nucleation scheme
USE mo_cmp_timer,          ONLY: & ! already imported for the module: timer_start,timer_stop  
                                 timer_cmp_0inits, timer_cmp_1cbasetop,        &
                                 timer_cmp_1airdens, timer_cmp_2zerolocaltend, &
                                 timer_cmp_3melting, timer_cmp_3subandevap,    &
                                 timer_cmp_4cicesed, timer_cmp_4incloudwi,     &
                                 timer_cmp_5conddepevapsubl,                   &
                                 timer_cmp_5clearairevap,                      &
                                 timer_cmp_5incloudconddepevapsubl,            &
                                 timer_cmp_6freezingTs238,                     &
                                 timer_cmp_6freezingTl238,                     &
                                 timer_cmp_7cloudphyscandsurfprecip,           &
                                 timer_cmp_7warmclouds, timer_cmp_7coldclouds, &
                                 timer_cmp_7precipflux, timer_cmp_8updating,   &
                                 timer_cmp_9wetchem, timer_cmp_diagnostics,    &
                                 timer_cmp_7ccpfcoutside,             & !UP #777
                                 timer_cmp_inside, timer_cmp_cmp        !UP #777
 USE mo_submodel,          ONLY: lcmptend
 USE mo_cmp_diagn,         ONLY: tend_sacl, tend_inucl, icnc_bfrs238, icnc_pfrs238, &
                                 icnc_bfrl238, icnc_pfrl238, &
                                 icnc_dfrl238, icnc_dfrs238, &
                                 icnc_ctrfrl238, icnc_ctrfrs238, &
                                 tend_testdummy_1, & !, tend_testdummy_2
                                 !>>UP #821
                                 diag_delta_icnucl, diag_cntr_icnucl 
                                 !<<UP #821
USE mo_param_switches,     ONLY: cdnc_min_fixed ! added to get nclmi=0 to work
 !<<UP

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
  REAL(dp), INTENT(in)    :: pxtm1   (kbdim,klev,ktrac) !< tracer mmr                    (t-1)
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
  REAL(dp), INTENT(inout) :: pacdnc  (kbdim,klev)    !< cloud droplet number concentration (specified)

! Output:
  REAL(dp), INTENT(out), TARGET:: picnc(kbdim,klev)!< Ice crystal number concentration (ICNC) [1/m3]
  REAL(dp), INTENT(out) :: prelhum (kbdim,klev)    !< relative humidity
  REAL(dp), INTENT(out) :: pssfl   (kbdim)         !< surface snow flux (kg/m2/s)
  REAL(dp), INTENT(out) :: prsfl   (kbdim)         !< surface rain flux (kg/m2/s)

! Local variables
!
  INTEGER :: jl, jk, jm, jkk, jfrzmod, iqidx, ixidx, ierr

  INTEGER :: itop(kbdim,klev),         & !< flag for cloud tops
             ibas(kbdim,klev),         & !< flag for cloud bases
             icl_minusbas(kbdim,klev), & !< flag for all cloud levels excepted their base
             icl_minustop(kbdim,klev), & !<  flag for all cloud levels excepted their top (useless for now)
             iclbas(kbdim,klev),       & !< conv. cloud base
             itm1_look(kbdim,klev),    & !< index for temperature lookup table
             itp1_look(kbdim,klev)       !< itm1_look + 1

  REAL(dp) :: zb2, zgtp, zvth, zfuchs, zfre, zre !< for cirrus calculations

  REAL(dp):: zesw_tmp(kbdim,klev)    !< Temporary saturation vapor pressure w.r.t. water
  REAL(dp):: zes_tmp(kbdim,klev)     !< Temporary saturation vapor pressure w.r.t. water or ice
  REAL(dp):: zcorw(kbdim,klev)       !< Temporary variable needed for the calculation of the sat. vapor pressure
  REAL(dp):: zcor(kbdim,klev)        !< Temporary variable needed for the calculation of the sat. vapor pressure
  REAL(dp):: zlucua(kbdim,klev)      !< Temporary variable needed for the calculation of the sat. vapor pressure (t-1)
  REAL(dp):: zlucuaw(kbdim,klev)     !< Temporary variable needed for the calculation of the sat. vapor pressure (t-1)
  REAL(dp):: zlucuawp1(kbdim,klev)   !< Temporary variable needed for the calculation of the sat. vapor pressure (t)
  REAL(dp):: zlucuap1(kbdim,klev)    !< Temporary variable needed for the calculation of the sat. vapor pressure (t)
  REAL(dp):: zclcpre_2d(kbdim,klev)  !< fraction of grid box covered by precip (2D)
  REAL(dp):: zclcpre(kbdim)          !< fraction of grid box covered by precip (column-wise)
  REAL(dp):: zclcfi(kbdim)           !< fraction of grid box covered by sedimenting ice (column-wise)
  REAL(dp):: zcnd(kbdim)             !< condensation rate [kg/kg]
  REAL(dp):: zdep(kbdim)             !< deposition rate [kg/kg]
  REAL(dp):: zevp(kbdim)             !< evaporation of rain [kg/kg]
  REAL(dp):: zxievap(kbdim)          !< evaporation of cloud ice [kg/kg]
  REAL(dp):: zxlevap(kbdim)          !< evaporation of cloud water [kg/kg]
  REAL(dp):: zfrl(kbdim,klev)        !< freezing rate [kg/kg]
  REAL(dp):: zimlt(kbdim)            !< melting of ice if T>273 K [kg/kg]
  REAL(dp):: zsmlt(kbdim)            !< melting of snow [kg/kg]
  REAL(dp):: zrpr(kbdim)             !< rain formation rate [kg/kg]
  REAL(dp):: zspr(kbdim)             !< snow formation rate [kg/kg]
  REAL(dp):: zsub(kbdim)             !< sublimation of snow [kg/kg]
  REAL(dp):: zxtec(kbdim,klev)       !< tendency for detr. conv. cloud liq water or cloud ice (t)
  REAL(dp):: zxlte(kbdim)            !< tendency of cloud liquid water from detrainment [kg/kg/s]
  REAL(dp):: zxite(kbdim)            !< tendency of cloud ice from detrainment [kg/kg/s]
  REAL(dp):: zxiflux(kbdim)          !< flux of ice crystals falling into the grid box from above
  REAL(dp):: zsacl(kbdim)            !< accretion of snow flakes with cloud droplets [kg/kg], UP: i.e. riming
  REAL(dp):: zxlte2(kbdim)           !< temporary variable for liq water tendency from detrainment [kg/kg/s]
  REAL(dp):: zxite2(kbdim)           !< temporary variable for ice tendency from detrainment [kg/kg/s]
  REAL(dp):: zlsdcp(kbdim,klev)      !< latent heat of sublimation divided by
                                     !< the specific heat at constant pressure
  REAL(dp):: zlvdcp(kbdim,klev)      !< latent heat of vaporization divided by
                                     !< the specific heat at constant pressure
  REAL(dp):: zximlt(kbdim)           !< melting of the ice falling from above [kg/kg]
  REAL(dp):: ztp1tmp(kbdim,klev)     !< temporary value of the updated temperature (t) [K]
  REAL(dp):: zqp1tmp(kbdim)          !< temporary value of the updated specific humidity (t) [kg/kg]
  REAL(dp):: zxisub(kbdim)           !< sublimation of cloud ice [kg/kg]
  REAL(dp):: zxlb(kbdim,klev)        !< Cloud liquid water in the cloudy part of the grid box [kg/kg]
  REAL(dp):: zxib(kbdim,klev)        !< Cloud ice in the cloudy part of the grid box [kg/kg]
  REAL(dp):: zrho(kbdim,klev)        !< Air density [kg/m3]
  REAL(dp):: zclcov(kbdim)           !< Total cloud cover (overlap considered)
  REAL(dp):: zrho_rcp(kbdim,klev)    !< Inverse air density
  REAL(dp):: zqvi(kbdim)             !< Vertically integrated specific humidity [kg/m2]
  REAL(dp):: zxlvi(kbdim)            !< Vertically integrated cloud water [kg/m2]
  REAL(dp):: zxivi(kbdim)            !< Vertically integrated cloud ice [kg/m2]
  REAL(dp):: zbetaqt(kbdim)          !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zwide(kbdim)            !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zbetacl(kbdim)          !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zturbvar(kbdim)         !< Variable related to Tompkins cloud cover scheme 
  REAL(dp):: zturbskew(kbdim)        !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zconvvar(kbdim)         !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zconvskew(kbdim)        !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zvartg(kbdim)           !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zmicroskew(kbdim)       !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zgenti(kbdim)           !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zgentl(kbdim)           !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zxvarte(kbdim)          !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zxskewte(kbdim)         !< Variable related to Tompkins cloud cover scheme
  REAL(dp):: zgeoh(kbdim,klevp1)     !< Geopotential height at half levels
  REAL(dp):: zcdncact(kbdim,klev)    !< Number of newly activated cloud droplets
  REAL(dp):: zrprn(kbdim,klev)       !< Rain formation rate for number conc. [1/m3]
  REAL(dp):: zsprn(kbdim,klev)       !< Snow formation rate for number conc. [1/m3] !UP warning: self-collection wrongly included
  REAL(dp):: zsacln(kbdim,klev)      !< Accretion rate for number conc. [1/m3]
  REAL(dp):: zfrln(kbdim,klev)       !< Freezing rate for number conc. [1/m3]
  !>>SF #318
  REAL(dp),TARGET:: zcdnc(kbdim,klev) !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp),TARGET:: zcdnc_blackhole(kbdim,klev) !< Dummy 'black hole' variable
  REAL(dp),POINTER:: zcdnc_upd(:,:) !< Cloud droplet number concentration (CDNC) [1/m3]

  REAL(dp),TARGET:: zicnc_blackhole(kbdim,klev) !< Ice crystal number concentration (ICNC) [1/m3]
  REAL(dp),POINTER:: zicnc_upd(:,:) !< Dummy 'black hole' variable
  !SF note: I distinguish zcdnc_blackhole and zicnc_blackhole (rather than one single 'zblackhole' var)
  !         just in case someone wants later to do some intermediate diags...
  !         The primary idea is nevertheless that these vars are, as their name says, totally
  !         opaque and are simply necessary for handling smartly the 1m vs 2m cases...
  !<<SF #318
  REAL(dp):: zcdnc_burden(kbdim)     !< Vertically integrated CDNC [1/m2]
  REAL(dp):: zqlnuc(kbdim,klev)      !< Nucleation rate of CDNC [1/m3]
  REAL(dp):: zqlnuccvh(kbdim,klev)   !< Nucleated CDNC at base of convective clouds - temporary var.
  REAL(dp):: zqlnuccv(kbdim,klev)    !< Nucleated CDNC at base of convective clouds
  REAL(dp):: zqlnuc_bas(kbdim,klev)  !< Nucleated CDNC at base of stratiform clouds
  REAL(dp):: zcdnc_bas(kbdim,klev)   !< CDNC at base of stratiform clouds
  REAL(dp):: zicncq(kbdim,klev)      !< Temporary value for ICNC [1/m3]
  REAL(dp):: zicnc_burden(kbdim)     !< Vertically integrated ICNC [1/m2]
  REAL(dp):: zascs(kbdim,klev)       !< Number of soluble aerosols needed for homogeneous freezing [1/m3]
  REAL(dp):: zrid(kbdim,klev)        !< Mean volume radius of ice crystals [m]
  REAL(dp):: zqsi(kbdim,klev)        !< Saturation specific humidity w.r.t. ice [kg/kg]
  REAL(dp):: zninucl(kbdim, klev)    !< number conc. of newly nucleated IC [1/m3]
  REAL(dp):: zqinucl(kbdim, klev)    !< mixing ratio of newly nucleated IC [kg/kg]
  REAL(dp):: zri(kbdim, klev)        !< size of newly nucleated IC [m]
  REAL(dp):: zreffct(kbdim)          !< cloud top effective radius [m]
  REAL(dp):: ztau1i(kbdim)           !< ice cloud optical depth - visible wavelength
  REAL(dp):: znidetr(kbdim, klev)    !< ICNC from detrainment [1/m3]
  REAL(dp):: zfracdusol(kbdim,klev)  !< Fraction of dust aerosols in all soluble mixed modes
  REAL(dp):: zfracduai(kbdim,klev)   !< Fraction of dust in the insoluble accumulation mode 
  REAL(dp):: zfracduci(kbdim,klev)   !< Fraction of dust in the insoluble coarse mode 
  REAL(dp):: zfracbcsol(kbdim,klev)  !< Fraction of BC in all soluble mixed modes
  REAL(dp):: zfracbcinsol(kbdim,klev)!< Fraction of BC in all insoluble modes
  REAL(dp):: zxifluxn(kbdim)         !< flux of ice crystals number falling into the grid box from above
  REAL(dp):: zrwetki(kbdim,klev)     !< wet radius, Aitken insoluble mode [m]
  REAL(dp):: zrwetai(kbdim,klev)     !< wet radius, accumulation insoluble mode [m]
  REAL(dp):: zrwetci(kbdim,klev)     !< wet radius, coarse insoluble mode [m]  
  REAL(dp):: zdz(kbdim,klev)         !< layer thickness [m]
  REAL(dp):: zdp(kbdim,klev)         !< pressure difference of the layer [Pa]
  REAL(dp):: zdpg(kbdim,klev)        !< delta p over g [kg/m2]
  REAL(dp):: zaaa(kbdim,klev)        !< Air density correction needed for the ice crystal fall velocity
  REAL(dp):: zviscos(kbdim,klev)     !< Dynamic viscosity of water in air
  REAL(dp):: zqswp1(kbdim,klev)      !< Saturation specific humidity w.r.t. water (t) [kg/kg]
  REAL(dp):: zqsip1(kbdim,klev)      !< Saturation specific humidity w.r.t. ice (t) [kg/kg]
  REAL(dp):: zqsw(kbdim,klev)        !< Saturation specific humidity w.r.t. water (t-1) [kg/kg]
  REAL(dp):: zastbstw(kbdim,klev)    !< Thermodynamic term needed for water nucleation
  REAL(dp):: zastbsti(kbdim,klev)    !< Thermodynamic term needed for ice nucleation
  REAL(dp):: zmmean(kbdim,klev)      !< Volume mean ice crystal radius [m]
  REAL(dp):: zxifallmc(kbdim,klev)   !< Fall velocity of ice crystal mass
  REAL(dp):: zalfased(kbdim,klev)    !< Parameter needed for the ice crystal fall velocity
  REAL(dp):: zbetased(kbdim,klev)    !< Parameter needed for the ice crystal fall velocity
  REAL(dp):: zesw_2d(kbdim,klev)     !< Saturation vapor pressure w.r.t. water [Pa]
  REAL(dp):: zesi(kbdim,klev)        !< Saturation vapor pressure w.r.t. ice [Pa]
  REAL(dp):: zsusatw(kbdim,klev)     !< Supersat. with respect to water
  REAL(dp):: zsusatw_evap(kbdim,klev)!< Subsat. w.r.t. water
  REAL(dp):: zsusatix(kbdim,klev)    !< Supersat. with respect to ice
  REAL(dp):: zvervx(kbdim,klev)      !< Updraft velocity [cm/s]
  REAL(dp):: zicesub(kbdim,klev)     !< Subsat. w.r.t. ice
  REAL(dp):: zqrho(kbdim,klev)       !< Inverse air density [m3/kg]
  REAL(dp):: zxsp1(kbdim)            !< DOPPELT ??? - see line 355
  REAL(dp):: zxidtstar(kbdim)        !< ???
  REAL(dp):: zlc(kbdim)              !< Latent heat of vaporation/sublimation over specific heat depending on temp.
  REAL(dp):: zxldtstar(kbdim)        !< ???
  REAL(dp):: zxlm1evp(kbdim)         !< Evaporated cloud water ???
  REAL(dp):: zxim1evp(kbdim)         !< Evaporated cloud ice???
  REAL(dp):: zxldt(kbdim)            !< Local tendency of cloud water needed for condensational growth ???
  REAL(dp):: zxidt(kbdim)            !< Local tendency of cloud water needed for depositional growth ???
  REAL(dp):: zqsm1(kbdim,klev)       !< Saturation specific humidity (t-1) over ice/water depending on temp. [kg/kg]
  REAL(dp):: zqp1(kbdim)             !< Specific humidity [kg/kg] (t)
  REAL(dp):: zdtdt(kbdim)            !< Local temperature tendency [K]
  REAL(dp):: zdqsat(kbdim)           !< ???
  REAL(dp):: ztp1(kbdim)             !< Temperature [K] (t)
  REAL(dp):: zdqsdt(kbdim)           !< Difference in sat. specific humidity between t and t-1
  REAL(dp):: zqst1(kbdim)            !< Saturation specific humidity (t) over ice/water depending on temp. [kg/kg]
  REAL(dp):: zqvdt(kbdim)            !< Local specific humidity tendency [kg/kg]
  REAL(dp):: zxip1(kbdim,klev)       !< Cloud ice mass mixing ratio (t) [kg/kg]
  REAL(dp):: zqsp1tmp(kbdim)         !< Temporary new sat. specific humidity over ice/water depending on temp. [kg/kg]
  REAL(dp):: zclcstar(kbdim)         !< Minimum of cloud cover and precipitation cover
  REAL(dp):: zxrp1(kbdim)            !< Rain mixing ratio (t) [kg/kg]
  REAL(dp):: zauloc(kbdim)           !< Part of the grid box allowed to participation in accretion with
                                     !< newly formed condensate, UP: i.e. riming
  REAL(dp):: iqidx_1d(kbdim)         !< !SF why real ???
  REAL(dp):: zbap1(kbdim)            !< ???
  REAL(dp):: zgent(kbdim)            !< ???
  REAL(dp):: ixidx_1d(kbdim)         !< !SF why real ???
  REAL(dp):: zqtau(kbdim)            !< ???
  REAL(dp):: zxilb(kbdim)            !< Sum of in-cloud water and ice mass mixing ratio [kg/kg]
  REAL(dp):: zbbap1(kbdim)           !< ???
  REAL(dp):: zbqp1(kbdim)            !< ???
  REAL(dp):: zqcdif(kbdim)           !< ???
  REAL(dp):: zkair(kbdim,klev)       !< Thermal conductivity of air ???
  REAL(dp):: zvervmax(kbdim,klev)    !< Threshold vertical velocity
  REAL(dp):: zrice(kbdim,klev)       !< Volume mean ice crystal radius
  REAL(dp):: zeta(kbdim,klev)        !< Variable needed for the Bergeron-Findeisen process
  REAL(dp):: zdv(kbdim,klev)         !< ???
  REAL(dp):: zlwc_strat(kbdim,klev)  !< Liquid water content in the stratiform part
  REAL(dp):: zlwc_detr(kbdim,klev)   !< Liquid water content in the detrained part
  REAL(dp):: zlwc_tot(kbdim,klev)    !< Liquid water content in the strat+detrained parts
  REAL(dp):: zmratepr(kbdim,klev)    !< Rain formation rate in cloudy part of the grid box [kg/kg]
  REAL(dp):: zmrateps(kbdim,klev)    !< Snow formation rate in cloudy part of the grid box  [kg/kg]
  REAL(dp):: zfrain(kbdim,klev)      !< Rain flux before evaporation [kg/m2/s]
  REAL(dp):: zfsnow(kbdim,klev)      !< Snow flux before sublimation [kg/m2/s]
  REAL(dp):: zfevapr(kbdim,klev)     !< Evaporation of rain [kg/m2/s]
  REAL(dp):: zfsubls(kbdim,klev)     !< Sublimation of snow [kg/m2/s]
  REAL(dp):: zmlwc(kbdim,klev)       !< In-cloud liquid water mass mixing ratio before rain formation [kg/kg]
  REAL(dp):: zmiwc(kbdim,klev)       !< In-cloud ice mass mixing ratio before snow formation [kg/kg]
  REAL(dp):: zmsnowacl(kbdim,klev)   !< Accretion rate of snow with cl. drop. in cloudy part of the grid box [kg/kg] 
                                     !< UP: i.e. riming
  REAL(dp) :: zcd2ic(kbdim)          !< amount of cloud droplets erroneously present at temp < cthomi that
                                     !< should be credited to ice crystal number concentration
  REAL(dp) :: zcd2unphys(kbdim)      !< amount of cloud droplets erroneously present at temp < cthomi that
                                     !< are unphysical
  REAL(dp) :: zcd2ic_2d(kbdim,klev)  !< same as zcd2ic, 2D version
  REAL(dp) :: zcd2unphys_2d(kbdim,klev)  !< same as zcd2unphys, 2D version
  REAL(dp) :: zic2cd(kbdim)          !< amount of ice crystals erroneously present at temp > tmelt that
                                     !< should be credited to cloud droplet number concentration
  REAL(dp) :: znicex(kbdim, klev)         !< number of newly formed ice crystals 
  REAL(dp) :: zapnx(kbdim,klev,nfrzmod)   !< aerosol number available for homogeneous freezing [1/cm3]
  REAL(dp) :: zaprx(kbdim,klev,nfrzmod)   !< radius of aerosols available for homogeneous freezing  [cm]
  REAL(dp) :: zapsigx(kbdim,klev,nfrzmod) !< std. dev. of aerosols available for homogeneous freezing 
  REAL(dp) :: zap(kbdim,klev)             !< total number of aerosols available
  REAL(dp) :: zaclc_tm1(kbdim,klev)       !< cloud cover at t-1 [1]
  REAL(dp) :: zcvcbot (kbdim)             !< conv. cloud base index
  REAL(dp) :: zwcape  (kbdim)             !< CAPE contrib. to conv. vertical velocity [m s- 1]
  REAL(dp) :: ztconv  (kbdim,klev)        !< temperature as it was in the conv scheme [K]
!>>DN #475 (min CDNC)
  REAL(dp) :: zcdnc_min  (kbdim,klev)     !< Minimum CDNC concentration computed from maximum radius [1/m3] 
!<<DN #475 (min CDNC)


!>>UP SLF #839
  REAL(dp) :: zxlp1_slf(kbdim)     !< temp for SLF diagnostic: liquid water content [kg/kg]
  REAL(dp) :: zxip1_slf(kbdim)     !< temp for SLF diagnostic: ice water content [kg/kg]
  REAL(dp) :: zxtotp1_slf(kbdim)   !< temp for SLF diagnostic: total water content [kg/kg]
!<<UP SLF
  LOGICAL :: ll_mlt(kbdim)         !< switch that traces if temperature is above melting point
  LOGICAL :: ll_precip(kbdim)      !< switch that traces the presence of precipitation
  LOGICAL :: ll_falling_ice(kbdim) !< switch that traces the presence of falling ice
  LOGICAL :: ll_icecl(kbdim,klev)  !< switch that traces the presence of ice cloud
  LOGICAL :: ll_liqcl(kbdim,klev)  !< switch that traces the presence of liquid cloud
  LOGICAL :: ll_cc(kbdim)          !< switch that traces the presence of cloud cover
  LOGICAL :: ll_frz_below_238K(kbdim) !< physical condition for freezing below 238K to occur
  LOGICAL :: ll_mxphase_frz(kbdim) !< physical condition for mixed-phase freezing to occur
  LOGICAL :: ll_WBF(kbdim)         !< physical condition for WBF process to occur
  LOGICAL :: ll_prcp_warm(kbdim)   !< physical condition for warm precip to form

  LOGICAL, PARAMETER :: nosize = .true. !< .true. ---> aerosol size effects are ignored for homogeneous freezing

  !Temporary vars that do not have any single physical meaning          
  REAL(dp) :: ztmp1_2d(kbdim,klev), ztmp2_2d(kbdim,klev), ztmp3_2d(kbdim,klev), ztmp4_2d(kbdim,klev)
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim), ztmp4(kbdim), ztmp5(kbdim)
  !>>UP emulator nic_cirrus
  REAL(dp) :: zicncq_1(kbdim, klev), zicncq_2(kbdim, klev), zninucl_1(kbdim, klev), zninucl_2(kbdim, klev)
  REAL(dp) :: zdep_1(kbdim), zdep_2(kbdim)
  !<<UP
  LOGICAL  :: ll1(kbdim), ll2(kbdim), ll3(kbdim), lo2(kbdim)

  LOGICAL :: ll_cv(kbdim,klev), &
             ll_ice(kbdim,klev), &
             ll1_2d(kbdim,klev), ll2_2d(kbdim,klev), ll3_2d(kbdim,klev), ll4_2d(kbdim,klev), &
             lo2_2d(kbdim,klev)

!SF ToDo: lo2 and lo2_2d should be renamed with something meaningful, and a function could be created to avoid 
!         code duplication each time they are computed

  REAL(dp), PARAMETER :: zdummy = 1._dp !SF this is a dummy var to handle special non-physical cases 
                                        !   in a MERGE statement. This var is solely needed to prevent
                                        !   spurious calculations like div by zero for example

!>>SF Kasja diags
  REAL(dp):: zsacin(kbdim), zselfn(kbdim), zsecprodn(kbdim), zrautn(kbdim), zracl(kbdim), &
             zracln(kbdim), zaggn(kbdim), zmltn(kbdim), zsedi(kbdim), zsedin(kbdim), zfrzn(kbdim), zraut(kbdim),&
             zsaut(kbdim), zsaci(kbdim), &
!<<SF Kasja diags
!>>UP #783
             zfrznhet(kbdim), zfrznhom(kbdim)
!<<UP #783
!>>DN: new diags
  REAL(dp):: zevabf(kbdim),zlwcbf(kbdim),zncdetr(kbdim,klev),&
             zcorric1(kbdim,klev),zcorric2(kbdim,klev),zcorric3(kbdim,klev),zcorric4(kbdim),zcorric5(kbdim),zcorric6(kbdim),&
             zcorric7(kbdim),&
             zcorrcd1(kbdim,klev),zcorrcd2(kbdim,klev),zcorrcd3(kbdim,klev),zcorrcd4(kbdim,klev),zcorrcd5(kbdim,klev),&
             zcorrcd6(kbdim),zcorrcd7(kbdim),zcorrcd8(kbdim),zcorrcd9(kbdim),zcorrcd10(kbdim),&
             zcorrxi(kbdim),zcorrxi2(kbdim),zcorrxl(kbdim),zcorrxl2(kbdim),zactxlte(kbdim,klev),zactxite(kbdim,klev),&
             zxttecdnc(kbdim,klev),zxtteicnc(kbdim,klev),zsecprod(kbdim),zmltn2(kbdim),ztestcd(kbdim,klev),ztestic(kbdim,klev),&
             ztestlwc(kbdim,klev),ztestiwc(kbdim,klev)
!<<DN: new diags  
!>>UP: new diags for activation correction terms
! and divide correction terms further, #783
  REAL(dp):: zcorrcd2unphys_2d(kbdim,klev)
  REAL(dp):: zcorrcd1_1(kbdim,klev) 
  REAL(dp):: zcorrcd1_2(kbdim,klev) 
  REAL(dp):: zcorrcd1_3(kbdim,klev) 
  REAL(dp):: zcorrcd1_4(kbdim,klev) 
  REAL(dp):: zcorrcd2_1(kbdim,klev) 
  REAL(dp):: zcorrcd2_2(kbdim,klev) 
  REAL(dp):: zcorric1_1(kbdim,klev) 
!>>UP #783.3
  REAL(dp):: zcorric1_2(kbdim,klev)
  REAL(dp):: zcorric1_3(kbdim,klev)
  REAL(dp):: zcorrcd2_3(kbdim,klev)
  REAL(dp):: zcd2ic_2d_ic(kbdim,klev) ! zcd2ic, but ... (see below)
  REAL(dp):: zcd2ic_ic(kbdim) ! zcd2ic, but only those that go to IC (less than
!zcd2ic, because of the maximum IC condition 
  REAL(dp):: qnucl_upd(kbdim,klev)
!<<UP #783.3
!<<UP: new diags
!>>UP #797
  REAL(dp):: zicnc_diff(kbdim)
  REAL(dp):: zxip1_diff(kbdim)
!<<UP #797
!>>UP #821
  REAL(dp) :: zcmpsimpl_icnucl(kbdim,klev)   ! receives input for prescribed icnucl
  REAL(dp) :: zcmpsimpl_rime(kbdim,klev)     ! receives input for prescribed rime
  REAL(dp) :: zcmpsimpl_icaccr(kbdim,klev)   ! receives input for prescribed icaccr
  REAL(dp) :: zcmpsimpl_subfi(kbdim,klev)    ! receives input for prescribed subfi
  REAL(dp) :: zcmpsimpl_subfs(kbdim,klev)    ! receives input for prescribed subfs
  REAL(dp) :: zcmpsimpl_sci(kbdim,klev)      ! receives input for prescribed sci
  REAL(dp) :: zcmpsimpl_mlt(kbdim,klev)      ! receives input for prescribed mlt
!<<UP #821
!>>UP lifetime diagnostics
  REAL(dp) :: z_lwc_sources(kbdim, klev)
  REAL(dp) :: z_lwc_sinks(kbdim, klev)
  REAL(dp) :: z_iwc_sources(kbdim, klev)
  REAL(dp) :: z_iwc_sinks(kbdim, klev)
  REAL(dp) :: z_cdnc_sources(kbdim, klev)
  REAL(dp) :: z_cdnc_sinks(kbdim, klev)
  REAL(dp) :: z_icnc_sources(kbdim, klev)
  REAL(dp) :: z_icnc_sinks(kbdim, klev)
!<<UP

!--- This switches control most subprocesses.
! ToDo: their value is so far set here, but could be placed in a 
!       namelist if necessary
  LOGICAL :: lctrl_mlt_snow_ice    = .TRUE. !< Controls melting of snow and ice
  LOGICAL :: lctrl_subl_evap       = .TRUE. !< Controls sublim. of snow and ice, and evap. of rain
  LOGICAL :: lctrl_sed_ice         = .TRUE. !< Controls sedimentation of ice
  LOGICAL :: lctrl_upd_incl_wat    = .TRUE. !< Controls in-cloud water changes processes
  LOGICAL :: lctrl_frz_below_238K  = .TRUE. !< Controls freezing below 238K 
  LOGICAL :: lctrl_het_mxphase_frz = .TRUE. !< Controls heterogeneous freezing between 238K and 273K
  !LOGICAL :: lctrl_WBF             = .TRUE. !< Controls WBF process  <<MA orig. uncommented, #765
  LOGICAL :: lctrl_prcp_warm       = .TRUE. !< Controls warm cloud microphysics processes
  LOGICAL :: lctrl_prcp_cold       = .TRUE. !< Controls cold cloud microphysics processes
  LOGICAL :: lctrl_upd_flx         = .TRUE. !< Controls the update of precip fluxes
  LOGICAL :: lctrl_upd_tend        = .TRUE. !< Controls the update of tendencies and crucial variables
  LOGICAL :: lctrl_diags           = .TRUE. !< Controls the call to (most) diagnostics
!davidn SLF
  REAL(dp):: ztp1h(kbdim,klevp1) ! stores the temperature interpolated to layer boundaries
  LOGICAL  :: ll4(kbdim),ll5(kbdim),ll6(kbdim),ll7(kbdim),ll8(kbdim)
!davidn

!-----------------------------------------------------------------------------------
! Executable statements

!-- 0. Initialisations:

  !>>UP timers
  IF (ltimer) THEN
     CALL timer_start(timer_cmp_inside)
     CALL timer_start(timer_cmp_0inits)
     CALL timer_start(timer_cmp_cmp)
  END IF
  !<<UP

  zqinucl(1:kproma,:) = 0._dp
  zninucl(1:kproma,:) = 0._dp
  zreffct(1:kproma)   = 0._dp
  ztau1i(1:kproma)    = 0._dp

  zcdnc_burden(1:kproma) = 0._dp
  zcdnc_bas(1:kproma,:)  = 0._dp 
  zqlnuc(1:kproma,:)     = 0._dp
  zqlnuc_bas(1:kproma,:) = 0._dp
  zicnc_burden(1:kproma) = 0._dp
  zri(1:kproma,:)        = 1.e-6_dp

  zsusatw(1:kproma,:)      = 0._dp
  zsusatw_evap(1:kproma,:) = 0._dp
  zastbstw(1:kproma,:)     = 0._dp
  zastbsti(1:kproma,:)     = 0._dp

  zfrl(1:kproma,:)  = 0._dp
  zfrln(1:kproma,:) = 0._dp
  zsprn(1:kproma,:) = 0._dp

  zclcpre(1:kproma)  = 0._dp 
  zclcfi(1:kproma)   = 0._dp 
  zxiflux(1:kproma)  = 0._dp 
  zxifluxn(1:kproma) = 0._dp

!UP #783.3
! _upd is used here so that the new nucleation can be corrected before being
! added to the rate
  qnucl_upd(1:kproma,:) = 0._dp
!<<UP
!>>UP #797
  zicnc_diff(1:kproma) = 0._dp
  zxip1_diff(1:kproma) = 0._dp
!<<UP #797

  zaclc_tm1(1:kproma,:) = cloud_cover_duplic(1:kproma,:,krow)

  !>>SF initialization of most intent(out) variables compulsory over the whole array
  !     because most of them are passed back to stream elements in physc, which would cause
  !     a crash by printing out some undefined parts of the stream element array otherwise
    picnc(1:kbdim,klev)     = 0._dp
    !UP: so why is picnc above only initialized for klev, not 1:klev?
    prelhum(1:kbdim,1:klev) = 0._dp
    pssfl(1:kbdim)          = 0._dp
    prsfl(1:kbdim)          = 0._dp 
  !<<SF
  
!>>SF ToDo: put that in an initialization routine, that is called only once before the time loop
!---Computational constants:
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
!<<SF

!--- Get several utility variables:
  CALL get_util_var( &
          !-- IN
          kproma, kbdim, ktdia, klev, klevp1, &
          paphm1(:,:), pgeo(:,:), papm1(:,:), ptm1(:,:), &
          !-- OUT
          zgeoh(:,:), zdp(:,:), zdpg(:,:), &
          zdz(:,:), zaaa(:,:), zviscos(:,:) )

  !--- Retrieve some quantities from boundary condition scheme 
    IF (lstart .OR. lresume) THEN
       IF (lconv) THEN !csld #455
          CALL bc_find('Convective cloud base index', ibc_cvcbot, ierr=ierr)
          CALL bc_find('CAPE contrib. to conv. vertical velocity', ibc_wcape, ierr=ierr)
          CALL bc_find('Temperature in convective scheme', ibc_tconv, ierr=ierr) !SF #368
          CALL bc_find('Detrained condensate', ibc_detr_cond, ierr=ierr) !SF #518
       ENDIF !csld #455
    ENDIF

      !>>UP #821
    IF (lcmpsimpl_prescr) THEN
      IF (ncmpsimpl_prescr_icaccr /= 0) THEN
             CALL bc_apply(ibc_cmpsimpl_icaccr, kproma, krow, zcmpsimpl_icaccr(1:kproma,:))
      ENDIF
      IF (ncmpsimpl_prescr_rime /= 0) THEN
             CALL bc_apply(ibc_cmpsimpl_rime, kproma, krow, zcmpsimpl_rime(1:kproma,:))
      ENDIF
      IF (ncmpsimpl_prescr_icnucl /= 0) THEN
             CALL bc_apply(ibc_cmpsimpl_icnucl, kproma, krow, zcmpsimpl_icnucl(1:kproma,:))
      ENDIF
      IF (ncmpsimpl_prescr_subfis /= 0) THEN
             CALL bc_apply(ibc_cmpsimpl_subfi, kproma, krow, zcmpsimpl_subfi(1:kproma,:))
             CALL bc_apply(ibc_cmpsimpl_subfs, kproma, krow, zcmpsimpl_subfs(1:kproma,:))
      ENDIF
      IF (ncmpsimpl_prescr_sci /= 0) THEN
             CALL bc_apply(ibc_cmpsimpl_sci, kproma, krow, zcmpsimpl_sci(1:kproma,:))
      ENDIF
      IF (ncmpsimpl_prescr_mlt /= 0) THEN
             CALL bc_apply(ibc_cmpsimpl_mlt, kproma, krow, zcmpsimpl_mlt(1:kproma,:))
      ENDIF
    ENDIF
    !<<UP #821

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

    !>>SF #318
    NULLIFY(zcdnc_upd)
    IF (IAND(nclmi_progn,1) /= 0) THEN !prognostic CDNC (nclmi_progn==1 .or. 3)
        zcdnc_upd => zcdnc
    ELSE ! diagnostic CDNC
        ! Initialization of CDNC to the land vs ocean vertical profile as defined in physc
        ! zcdnc(1:kproma,:) = pacdnc(1:kproma,:) ! UP: original line
        ! Give a minimal CDNC value, see Issue #759
        zcdnc(1:kproma,:) = MAX(pacdnc(1:kproma,:),cdnc_min_fixed*1.e6_dp) ! UP: I think the units are cm-3 from mo_activ.f90 so we
        ! UP TODO: use minimum_CDNC here instead, cdnc_min_fixed was just a
        ! quick fix
        ! need to convert to m-3 here
        ! Added above line instead of the one before to fix the nclmi=0 run
        zcdnc_upd => zcdnc_blackhole
    ENDIF

    NULLIFY(zicnc_upd)
    IF (IAND(nclmi_progn,2) /= 0) THEN !prognostic ICNC (nclmi_progn==2 .or. 3)
        zicnc_upd => picnc !UPcomment: pointer assignment
    ELSE ! diagnostic ICNC
        ! Initialisation of ice nuclei number concentration (1/M**3) acc. to Pruppacher/Klett (1997)
        ! See Seinfeld/Pandis p.808
        !picnc(1:kproma,:) = MAX(MIN(1.e3_dp*exp(0.6_dp*(253.2_dp-ptm1(1:kproma,:))),1.e6_dp),cqtmin)
        ! UP fix issue #758>>
        ! Initialisation of ice nuclei number concentration acc. to Cooper (1986) ["Ice Initiation in Natural Clouds, 1986, William
        ! A. Cooper. Chapter 4 in: Roscoe, Braham, Jr. “Precipitation Enhancement - a Scientific Challenge.” Meteorological
        ! Monographs 21, no. 43 (1986): 178.]
        !picnc(1:kproma,:) = MAX(MIN(1e3*10**(-2.35-0.135*(ptm1(1:kproma,:) - 273.15)),1.e6_dp),cqtmin)
        !UP: Change maximal ICNC for IssueID #770
        picnc(1:kproma,:) = MAX(MIN(1e3*10**(-2.35-0.135*(ptm1(1:kproma,:) - 273.15)),1.e5_dp),cqtmin)
        ! << UP fix issue #758
        zicnc_upd => zicnc_blackhole
    ENDIF
    !<<SF #318
      
    !>>UP timers, #777
    IF (ltimer) THEN
       CALL timer_stop(timer_cmp_0inits)
       CALL timer_start(timer_cmp_1airdens)
    END IF
    !<<UP

!------------------------------------------------------------------
!       1.   Top boundary conditions, air density

!       1.2   Air density

  zrho(1:kproma,:)     = papm1(1:kproma,:)/(rd*ptvm1(1:kproma,:))
  zrho_rcp(1:kproma,:) = 1._dp / zrho(1:kproma,:)          
  zqrho(1:kproma,:)    = 1.3_dp * zrho_rcp(1:kproma,:)

!SF Detrained liq water/ice calc.: 

!SF ToDo: convert this loop to array syntax
  DO 122 jk = 1,klev  

     ll1(1:kproma) = (zxtec(1:kproma,jk) > 0.0_dp)

     ztmp1(1:kproma)            = twc_conv(1:kproma,jk,krow) + ztmstdt*zxtec(1:kproma,jk)
     twc_conv(1:kproma,jk,krow) = MERGE(ztmp1(1:kproma), twc_conv(1:kproma,jk,krow), ll1(1:kproma))

     ztmp1(1:kproma)             = conv_time(1:kproma,jk,krow) + zdtime
     conv_time(1:kproma,jk,krow) = MERGE(ztmp1(1:kproma), conv_time(1:kproma,jk,krow), ll1(1:kproma))
 
     pqtec(1:kproma,jk)  = MAX(pqtec(1:kproma,jk),0.0_dp)

!--- Included for prognostic CDNC/IC scheme ----------------------------

     ! calculate cloud droplet number concentration 
     zcdnc_upd(1:kproma,jk) = zrho(1:kproma,jk)*(pxtm1(1:kproma,jk,idt_cdnc) + ztmst*pxtte(1:kproma,jk,idt_cdnc))
     zcdnc_upd(1:kproma,jk) = MAX(zcdnc(1:kproma,jk),cqtmin)

     ! calculate ice crystal number
     zicnc_upd(1:kproma,jk)  = zrho(1:kproma,jk)*(pxtm1(1:kproma,jk,idt_icnc) + ztmst*pxtte(1:kproma,jk,idt_icnc))
     zicnc_upd(1:kproma,jk)  = MAX(picnc(1:kproma,jk),cqtmin)

!>>DN: new diags
     zxttecdnc(1:kproma,jk) = pxtte(1:kproma,jk,idt_cdnc)
     zxtteicnc(1:kproma,jk) = pxtte(1:kproma,jk,idt_icnc)
     zactxlte(1:kproma,jk) = pxlte(1:kproma,jk)
     zactxite(1:kproma,jk) = pxite(1:kproma,jk)
     ztmp2(1:kproma) = picnc(1:kproma,jk)
     ztmp3(1:kproma) = zcdnc(1:kproma,jk)
     ztestcd(1:kproma,jk) = ztmst*pxtte(1:kproma,jk,idt_cdnc)
     ztestic(1:kproma,jk) = ztmst*pxtte(1:kproma,jk,idt_icnc)
     ztestlwc(1:kproma,jk) = ztmst*pxlte(1:kproma,jk)
     ztestiwc(1:kproma,jk) = ztmst*pxite(1:kproma,jk)
!<<DN: new diags

     !>>SF correct for inconsistencies relative to phases and temperature
!>>DN #475
     ztmp1(1:kproma) = (pxlm1(1:kproma,jk)+ztmst*pxlte(1:kproma,jk))*zrho(1:kproma,jk)/MAX(paclc(1:kproma,jk),eps)
     zcdnc_min(:,jk)   = minimum_CDNC(kbdim, kproma, ztmp1(:))
!<<DN #475
     !DN / SF #477 refined correction for crediting to ice *only* when cl. cover and LWC are non-zero
     ll1(1:kproma) = (ptm1(1:kproma,jk) <= cthomi) .AND. (zcdnc(1:kproma,jk) > zcdnc_min(1:kproma,jk))
     ll2(1:kproma) = (paclc(1:kproma,jk) >= epsec) .AND. &
                     (pxlm1(1:kproma,jk) + ztmst*pxlte(1:kproma,jk) >= epsec)

     !UP: correction of CD activation for T and cloud cover, 1d
     zcd2ic(1:kproma)     = MERGE(zcdnc(1:kproma,jk), 0._dp, ll1(1:kproma) .AND.       ll2(1:kproma))
     zcd2unphys(1:kproma) = MERGE(zcdnc(1:kproma,jk), 0._dp, ll1(1:kproma) .AND. .NOT. ll2(1:kproma))

     ll1(1:kproma) = (ptm1(1:kproma,jk) > tmelt) .AND. (picnc(1:kproma,jk) > icemin)
     zic2cd(1:kproma) = MERGE(picnc(1:kproma,jk), 0._dp, ll1(1:kproma))

     zcdnc_upd(1:kproma,jk) = zcdnc(1:kproma,jk) - zcd2ic(1:kproma) - zcd2unphys(1:kproma) + zic2cd(1:kproma)
!>>UP #783.3: need to put in a maximum condition for icnc here
! don't use more zcd2ic than what brings you to maximum picnc
     IF (lpr_corr) THEN 
        zcd2ic_ic(1:kproma) = &
                MERGE(MAX(icemax - (picnc(1:kproma,jk) - zic2cd(1:kproma)),0._dp),&
!zcd2ic(1:kproma)-((picnc(1:kproma,jk) + zcd2ic(1:kproma) - zic2cd(1:kproma)) - icemax), &
                zcd2ic(1:kproma),&
                (picnc(1:kproma,jk) + zcd2ic(1:kproma) - zic2cd(1:kproma) > icemax))
     ELSE
        zcd2ic_ic(1:kproma) = zcd2ic(1:kproma)
     ENDIF
! zicnc_upd(1:kproma,jk) = picnc(1:kproma,jk) + zcd2ic(1:kproma) - zic2cd(1:kproma) ! orig. line,
! replaced by line below for #783.3
zicnc_upd(1:kproma,jk) = picnc(1:kproma,jk) + zcd2ic_ic(1:kproma) - zic2cd(1:kproma)
!<<UP #783.3

!>>UP: new diags, #783
! divide below diagnostics up into smaller pieces
IF (lnewdiags .OR. lpr_corr) THEN
     zcorrcd1_1(1:kproma,jk) = MAX(cqtmin - zcdnc(1:kproma,jk), 0._dp)
     zcorric1_1(1:kproma,jk) = MAX(cqtmin - picnc(1:kproma,jk), 0._dp)
ENDIF
!<<UP
     zcdnc_upd(1:kproma,jk) = MAX(zcdnc(1:kproma,jk),cqtmin)
     zicnc_upd(1:kproma,jk) = MAX(picnc(1:kproma,jk),cqtmin)
!UP question: Is it correct that picnc is put to minimal cqtmin and not some ice
!minimal value?
     !<<SF end correct for inconsistencies

!>>DN: new diags
     zcorric1(1:kproma,jk) = (picnc(1:kproma,jk)-ztmp2(1:kproma))
     zcorrcd1(1:kproma,jk) = (zcdnc(1:kproma,jk)-ztmp3(1:kproma))
!<<DN: new diags
!>>UP: new diags, #783
! divide above diagnostics up into smaller pieces
IF (lnewdiags .OR. lpr_corr) THEN
     zcorrcd1_2(1:kproma,jk) = zcd2ic(1:kproma)
     zcorrcd1_3(1:kproma,jk) = zcd2unphys(1:kproma)
     zcorrcd1_4(1:kproma,jk) = zic2cd(1:kproma)
     ! #783.3
     zcorric1_2(1:kproma,jk) = zcd2ic_ic(1:kproma)
     zcorric1_3(1:kproma,jk) = zic2cd(1:kproma) !UP #783.3, later gets distinguished from cd1_4
ENDIF
!<<UP

     zicncq(1:kproma,jk) = picnc(1:kproma,jk)

122 END DO

!       1.3   Calculate cloud base and cloud top

  !>>UP timers, #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_1airdens)
     CALL timer_start(timer_cmp_1cbasetop)
  END IF
  !<<UP

  CALL get_cloud_bounds( &
          !-- IN
          kproma, kbdim, ktdia, klev, paclc, &
          !- OUT
          itop, ibas, icl_minustop, icl_minusbas)

  CALL set_lookup_index( &
          !-- IN
          kbdim, kproma, klev, ptm1(:,:), '2-m cloud micro (1)', &
          !-- OUT
          itm1_look(:,:) )

  itp1_look(1:kproma,:) = itm1_look(1:kproma,:) + 1
  itp1_look(1:kproma,:) = MAX(MIN(itp1_look(1:kproma,:),jptlucu2),jptlucu1)

  DO jk=klev,ktdia,-1
     DO jl=1,kproma
        zlucuaw(jl,jk)   = tlucuaw(itm1_look(jl,jk)) 
        zlucuawp1(jl,jk) = tlucuaw(itp1_look(jl,jk)) 
        zlucua(jl,jk)    = tlucua(itm1_look(jl,jk)) 
        zlucuap1(jl,jk)  = tlucua(itp1_look(jl,jk)) 
     END DO !jl
  END DO !jk

!SF water:
  CALL sat_spec_hum( &
          !-- IN
          kbdim, kproma, klev, papm1(:,:), zlucuaw(:,:), &
          !-- OUT
          zesw_tmp(:,:), zcorw(:,:), zqsw(:,:) )

  ! Saturation water vapour pressure with respect to water:
  zesw_2d(1:kproma,:) = zesw_tmp(1:kproma,:)*papm1(1:kproma,:)*rv/rd

  zsusatw(1:kproma,:) = pqm1(1:kproma,:)/zqsw(1:kproma,:)-1.0_dp
  zsusatw(1:kproma,:) = MAX(zsusatw(1:kproma,:),0.0_dp)

  zsusatw_evap(1:kproma,:) = pqm1(1:kproma,:)/zqsw(1:kproma,:)-1.0_dp
  zsusatw_evap(1:kproma,:) = MIN(zsusatw_evap(1:kproma,:), 0._dp)

  CALL sat_spec_hum( &
          !-- IN
          kbdim, kproma, klev, papm1(:,:), zlucuawp1(:,:), &
          !-- OUT
          zesw_tmp(:,:), zcorw(:,:), zqswp1(:,:))

!SF ice:
  CALL sat_spec_hum( &
          !-- IN
          kbdim, kproma, klev, papm1(:,:), zlucua(:,:), &
          !-- OUT
          zesi(:,:), zcor(:,:), zqsi(:,:))

  ! Saturation water vapour pressure with respect to ice:
  zesi(1:kproma,:) = zesi(1:kproma,:)*papm1(1:kproma,:)*rv/rd

  sice(1:kproma,:,krow) = pqm1(1:kproma,:)/zqsi(1:kproma,:)-1.0_dp
  sice(1:kproma,:,krow) = MAX(sice(1:kproma,:,krow), 0._dp)

  zicesub(1:kproma,:) = pqm1(1:kproma,:)/zqsi(1:kproma,:)-1.0_dp
  zicesub(1:kproma,:) = MIN(zicesub(1:kproma,:), 0._dp) 

  CALL sat_spec_hum( &
          !-- IN
          kbdim, kproma, klev, papm1(:,:), zlucuap1(:,:), &
          !-- OUT
          zes_tmp(:,:), zcor(:,:), zqsip1(:,:))

  ! Store supersaturation with respect to water in stream:
  swat(1:kproma,:,krow)  = zsusatw(1:kproma,:)

!Set some more utility variables:
!SF ToDo: put that in get_util_var
  zastbstw(1:kproma,:) = &
                         !zast for water:
                         alv * (alv/(rv*ptm1(1:kproma,:)) - 1.0_dp) / (ka*ptm1(1:kproma,:)) &
                         !zbst for water:
                       + rv*ptm1(1:kproma,:)/(2.21_dp/papm1(1:kproma,:)*zesw_2d(1:kproma,:))

  zkair(1:kproma,:) = 4.1867e-3_dp * (5.69_dp + 0.017_dp*(ptm1(1:kproma,:)-tmelt)) ! eq. 13-18a Pruppacher & Klett 

  zastbsti(1:kproma,:) = &
                         !zast for ice:
                         als * (als/(rv*ptm1(1:kproma,:)) - 1.0_dp) / (zkair(1:kproma,:)*ptm1(1:kproma,:)) &
                         !zbst for ice:
                       + rv*ptm1(1:kproma,:)/(2.21_dp/papm1(1:kproma,:)*zesi(1:kproma,:))

  !--- Interface to aerosol calculations (avoids HAM dependencies in cloud_micro_interface):
  !    (includes activation)

!>>UP ccnclim
  ll1_2d(1:kproma,:) = ( ibas(1:kproma,:) > 0 )                               .AND. & !cloud base lev
                       ( ptm1(1:kproma,:) > cthomi )                          .AND. &
                       ( (zcdnc(1:kproma,:) <= zcdnc_min(1:kproma,:)) .OR.  & 
                         (pxlte(1:kproma,:) >  0._dp)                 .OR.  &
                         (paclc(1:kproma,:) >  zaclc_tm1(1:kproma,:)) .OR.  &
                         (zsusatw(1:kproma,:) >  eps)                     )
  ll2_2d(1:kproma,:) = ( ptm1(1:kproma,:) > cthomi )                          .AND. &
                       ( (zcdnc(1:kproma,:) <= zcdnc_min(1:kproma,:)) .OR.  & 
                         (pxlte(1:kproma,:) >  0._dp)                 .OR.  &
                         (paclc(1:kproma,:) >  zaclc_tm1(1:kproma,:)) .OR.  &
                         (zsusatw(1:kproma,:) >  eps)                     )
!<<UP ccnclim

  CALL cloud_subm_1( &
          !-- IN
          kproma, kbdim, klev, krow, ktdia, &
          ptkem1, zwcape, pvervel, zrho, & 
          zrho_rcp, pxtm1, pxtte, ptm1, papm1, pqm1, zesw_2d, &
          !>>UP
          zcdnc_min, &
          ll1_2d, ll2_2d, ibas, &
          zdpg, &
          !<<UP
          !-- OUT
          zcdncact, zrwetki, zrwetai, zrwetci, &
          zfracdusol, zfracduai, zfracduci, zfracbcsol, zfracbcinsol, &
          zascs, zapnx, zaprx, zapsigx, ll_het )

!DN --> SF ToDo: add a subroutine for cloud droplet activation after cloud_subm_1
!DN --> SF ToDo: add a subroutine for ice crystal nucleation
!DN --> SF ToDo: handling of detrained condensate should be in 2 subroutines:(CDNC/ICNC) here and (LWC/IWC) in loop 831

!--- Convert the aerosol activation into the number of newly formed cloud droplets

  ll1_2d(1:kproma,:) = ( ibas(1:kproma,:) > 0 )                               .AND. & !cloud base lev
                       ( ptm1(1:kproma,:) > cthomi )                          .AND. &
                       ( (zcdnc(1:kproma,:) <= zcdnc_min(1:kproma,:)) .OR.  & 
                         (pxlte(1:kproma,:) >  0._dp)                 .OR.  &
                         (paclc(1:kproma,:) >  zaclc_tm1(1:kproma,:)) .OR.  &
                         (zsusatw(1:kproma,:) >  eps)                     )

  !SF: first computes newly formed cloud droplets at cloud bases:
  !>>UP #805
  IF (lemuphase_cdnuc) THEN
          ztmp1_2d(1:kproma,:) = eta_emu_cdnuc*(zcdncact(1:kproma,:) - zcdnc(1:kproma,:))
  ELSE
  !UP original
          ztmp1_2d(1:kproma,:) = zcdncact(1:kproma,:) - zcdnc(1:kproma,:) 
  ENDIF
  !<<UP #805
  ztmp1_2d(1:kproma,:) = MAX(0._dp, ztmp1_2d(1:kproma,:))
  zqlnuc(1:kproma,:)   = MERGE(ztmp1_2d(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))

  ztmp1_2d(1:kproma,:)  = zcdnc(1:kproma,:) + zqlnuc(1:kproma,:)
  zcdnc_upd(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:), zcdnc(1:kproma,:), ll1_2d(1:kproma,:))
  !>>UP 1CCN add diags
  IF (lnewdiags) THEN
          diag_cnt_cqtmins1(1:kproma,:,krow) = diag_cnt_cqtmins1(1:kproma,:,krow) + zdtime*MERGE(1, 0, zcdnc(1:kproma,:) > cqtmin)
          diag_cnt_cqtminl1(1:kproma,:,krow) = diag_cnt_cqtminl1(1:kproma,:,krow) + zdtime*MERGE(1, 0, zcdnc(1:kproma,:) < cqtmin)
          diag_cnt_ll12d_1(1:kproma,:,krow) = diag_cnt_ll12d_1(1:kproma,:,krow) + zdtime*MERGE(1, 0, ll1_2d(1:kproma,:))
  ENDIF
  !<<UP

  qnuc(1:kproma,:,krow) = qnuc(1:kproma,:,krow) + zdt*zqlnuc(1:kproma,:)

!>>DN: new diags
  DO jk = 1,klev  
     !>>UP #783.3: separate updates to be able to correct them before adding
     !them
     IF (lpr_corr) THEN
         qnucl_upd(1:kproma,jk) = qnucl_upd(1:kproma,jk) + zdt*zqlnuc(1:kproma,jk)
         !UP TODO: is there a line for dnucl missing here? 
     ELSE
     ! original formulation
         qnucl(1:kproma,jk,krow) = qnucl(1:kproma,jk,krow) + zdt*zqlnuc(1:kproma,jk)
         dnucl(1:kproma,krow)    = dnucl(1:kproma,krow)    + zdt*zqlnuc(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
     ENDIF
     !<<UP #783.3
  END DO
!<<DN: new diags
 
  !SF: then computes newly formed cloud droplets above cloud base
  !    by assuming that the number of nucleated cloud droplets is constant above cloud base
  !    (adiabatic parcel theory)

  DO jk=ktdia,klev
     DO jl=1,kproma
        jkk = MAX(jk,icl_minusbas(jl,jk))  !sets the level index to either relevant cloud base or itself
        zqlnuc_bas(jl,jk) = zqlnuc(jl,jkk) !holds in each cloud level the number of
                                           !newly formed cloud droplets at the base

        zcdnc_bas(jl,jk) = zcdnc(jl,jkk)   !holds in each cloud level the number of
                                           !activated cloud droplets at the base
     ENDDO !end jl
  ENDDO !end jk

  ll1_2d(1:kproma,:) = (icl_minusbas(1:kproma,:) > 0)     .AND. & !all cloud levels minus base
                       (zqlnuc_bas(1:kproma,:)   > 0._dp)         !newly cloud droplets at base

!>>DN: new diags
  ztmp3_2d(1:kproma,:)     = MERGE(zcdnc_bas(1:kproma,:),0._dp, ll1_2d(1:kproma,:))
  ztmp1_2d(1:kproma,:)    = MERGE(zcdnc_bas(1:kproma,:)-zqlnuc_bas(1:kproma,:),0._dp, ll1_2d(1:kproma,:))
  !UP: replace zcdnc by zcdnc_upd below? No, left of = must have
  !_update, right of = without
  ztmp2_2d(1:kproma,:)    = MERGE(zcdnc(1:kproma,:)-zqlnuc(1:kproma,:),0._dp, ll1_2d(1:kproma,:))
!UP: corrcd5 is not a correction term, but additional info (correction is
!already applied to qnucl)
  zcorrcd5(1:kproma,:) = ztmp1_2d(1:kproma,:)-ztmp2_2d(1:kproma,:)
  DO jk = 1,klev  
     !>>UP #783.3: separate updates to be able to correct them before adding
     !them
     IF (lpr_corr) THEN
             qnucl_upd(1:kproma,jk) = qnucl_upd(1:kproma,jk) + &
                  zdt*(ztmp3_2d(1:kproma,jk)-ztmp2_2d(1:kproma,jk))
         !UP TODO: is there a line for dnucl missing here? (same as above) 
    ELSE
    ! original formulation
             qnucl(1:kproma,jk,krow) = qnucl(1:kproma,jk,krow) + &
                  zdt*(ztmp3_2d(1:kproma,jk)-ztmp2_2d(1:kproma,jk))
             dnucl(1:kproma,krow)    = dnucl(1:kproma,krow) + &
                  zdt*(ztmp3_2d(1:kproma,jk)-ztmp2_2d(1:kproma,jk))*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
    ENDIF
    !<<UP #783.3
  END DO
!<<DN: new diags

!UP #787: remove below line because it is useless
!  zqlnuc(1:kproma,:)   = MERGE(zqlnuc_bas(1:kproma,:), zqlnuc(1:kproma,:), ll1_2d(1:kproma,:))

  ztmp1_2d(1:kproma,:)  = qnuc(1:kproma,:,krow) + zdt*(zcdnc_bas(1:kproma,:) - zcdnc(1:kproma,:))
  qnuc(1:kproma,:,krow) = MERGE(ztmp1_2d(1:kproma,:), qnuc(1:kproma,:,krow), ll1_2d(1:kproma,:))

  zcdnc_upd(1:kproma,:) = MERGE(zcdnc_bas(1:kproma,:), zcdnc(1:kproma,:), ll1_2d(1:kproma,:))
  !>>UP 1CCN add diags
  IF (lnewdiags) THEN
          diag_cnt_cqtmins2(1:kproma,:,krow) = diag_cnt_cqtmins2(1:kproma,:,krow) + zdtime*MERGE(1, 0, zcdnc(1:kproma,:) > cqtmin)
          diag_cnt_cqtminl2(1:kproma,:,krow) = diag_cnt_cqtminl2(1:kproma,:,krow) + zdtime*MERGE(1, 0, zcdnc(1:kproma,:) < cqtmin)
          diag_cnt_ll12d_2(1:kproma,:,krow) = diag_cnt_ll12d_2(1:kproma,:,krow) + zdtime*MERGE(1, 0, ll1_2d(1:kproma,:))
  ENDIF
  !<<UP

  !>>SF correct for inconsistencies relative to phases and temperature
  !DN / SF #477 refined correction for crediting to ice *only* when cl. cover and LWC are non-zero
  !UP: I think this is where the activate CDs are converted to ice for T<-35C,
  !and removed where there is no cloud cover, 2d
  ll1_2d(1:kproma,:) = (ptm1(1:kproma,:) <= cthomi) .AND. (zcdnc(1:kproma,:) > zcdnc_min(1:kproma,:))
  ll2_2d(1:kproma,:) = (paclc(1:kproma,:) >= epsec) .AND. &
                       (pxlm1(1:kproma,:) + ztmst*pxlte(1:kproma,:) >= epsec)
  
  zcd2ic_2d(1:kproma,:)     = MERGE(zcdnc(1:kproma,:), 0._dp, ll1_2d(1:kproma,:) .AND.       ll2_2d(1:kproma,:))
  zcd2unphys_2d(1:kproma,:) = MERGE(zcdnc(1:kproma,:), 0._dp, ll1_2d(1:kproma,:) .AND. .NOT. ll2_2d(1:kproma,:))

!>>DN: new diags
  ztmp1_2d(1:kproma,:) = zcdnc(1:kproma,:)
!<<DN: new diags

  zcdnc_upd(1:kproma,:) = zcdnc(1:kproma,:) - zcd2ic_2d(1:kproma,:) - zcd2unphys_2d(1:kproma,:)
!>>UP #783.3: need to put in a maximum condition for icnc here
! don't use more zcd2ic than what brings you to maximum picnc
     IF (lpr_corr) THEN 
        zcd2ic_2d_ic(1:kproma,:) = &
                MERGE(MAX(icemax - (picnc(1:kproma,:)),0._dp),&
                zcd2ic_2d(1:kproma,:),&
                (picnc(1:kproma,:) + zcd2ic_2d(1:kproma,:)  > icemax))
     ELSE
        zcd2ic_2d_ic(1:kproma,:) = zcd2ic_2d(1:kproma,:)
     ENDIF
  ! zicnc_upd(1:kproma,:) = picnc(1:kproma,:) + zcd2ic_2d(1:kproma,:) 
  ! #783.3 original formulation above replace by line below to only use what
  ! really adds to icnc
  zicnc_upd(1:kproma,:) = picnc(1:kproma,:) + zcd2ic_2d_ic(1:kproma,:)
!<<UP #783.3

  !UP comment: zcdnc_upd is updated with zcdnc, because for prognostic they point at the same thing, and for diagnostic zcdnc_upd
  !goes to a black hole anyways
!>>UP: new diags, #783
! divide below diagnostics up into smaller pieces
IF (lnewdiags .OR. lpr_corr) THEN
  zcorrcd2_1(1:kproma,:) = MAX(cqtmin - zcdnc(1:kproma,:), 0._dp) 
ENDIF
!<<UP
  !>>UP 1CCN add diags
  IF (lnewdiags) THEN
          diag_cnt_cqtmins3(1:kproma,:,krow) = diag_cnt_cqtmins3(1:kproma,:,krow) + zdtime*MERGE(1, 0, zcdnc(1:kproma,:) > cqtmin)
          diag_cnt_cqtminl3(1:kproma,:,krow) = diag_cnt_cqtminl3(1:kproma,:,krow) + zdtime*MERGE(1, 0, zcdnc(1:kproma,:) < cqtmin)
  ENDIF
  !<<UP
  zcdnc_upd(1:kproma,:) = MAX(zcdnc(1:kproma,:),cqtmin)
  
!>>DN: new diags
  !UP: replace zcdnc by zcdnc_upd below? No, left of = must have
  !_update, right of = without
  zcorrcd2(1:kproma,:) = zcdnc(1:kproma,:)-ztmp1_2d(1:kproma,:)
!<<DN: new diags
!UP #783.3, change to _ic
  zcorric2(1:kproma,:) = zcd2ic_2d_ic(1:kproma,:)
!>>UP: new diags for CD activation correction
  IF (lnewdiags .OR. lpr_corr) THEN 
      !>>UP: new diags, #783
      ! divide above diagnostics up into smaller pieces
      zcorrcd2_2(1:kproma,:) = zcd2unphys_2d(1:kproma,:)
      !UP #783.3
      zcorrcd2_3(1:kproma,:) = zcd2ic_2d(1:kproma,:)
      !<<UP
      zcorrcd2unphys_2d(1:kproma,:) = zcd2unphys_2d(1:kproma,:)
  ENDIF
!<<UP: new diags

  !<<SF correct for inconsistencies

!--- obtain a number of cloud droplets for the detrained cloud water from convective anvils
!SF Todo: all calculation related to detrainment should be encapsulated into an IF (lconv) statement.
  
  DO jk=ktdia,klev
     iclbas(1:kproma,jk) = NINT(zcvcbot(1:kproma)) 
     ll1_2d(1:kproma,jk) = (jk == iclbas(1:kproma,jk))
  ENDDO

  ll2_2d(1:kproma,:) = (iclbas(1:kproma,:) > 0) 

!>>DN #286: consistent ice/liquid split for detrained condensates

!>>SF ToDo: move that into get_util_var and create a function for computing vervx (also used in het_mxphase_freezing)
!--- set up criterion for the Wegener-Bergeron-Findeisen process
  !-- Updraft velocity:
  ztmp1_2d(1:kproma,:)    = 100._dp*fact_tke*SQRT(ptkem1(1:kproma,:))
  ztmp1_2d(1:kproma,klev) = 0.0_dp
  zvervx(1:kproma,:)      = -100._dp*rgrav*pvervel(1:kproma,:)*zrho_rcp(1:kproma,:) + ztmp1_2d(1:kproma,:)
!<<SF ToDo

#ifdef HAMMOZ
!>>gf #65 (SFnote: lorocirrus is only possible when nic_cirrus == 2)
!SF ToDo: remove rgrav from argument list (this is a constant)

!  Update updraft velocity to take into account orography, if relevant:

!>>UP ham timers
IF (ltimer) CALL timer_start(timer_ham_ifdef)
!<<UP

! UP note: this seems not related to HAM, but to CMP (since also it says in mo_orocirrus.f90: "You can use this module only for
! Bernd's cirrus cloud scheme (nic_cirrus=2)"), so I won't add it to timer_ham_totsum
   IF (lorocirrus) &
      CALL orocirrus_w( &
              !-- IN
              kproma, kbdim, klev, ktdia, &
              paphm1(:,:), &
              papm1(:,:), pgeo(:,:), ptm1(:,:), pqm1(:,:), um1(:,:,krow), &
              vm1(:,:,krow), orostd(:,krow), &
              orosig(:,krow), oromea(:,krow), &
              orogam(:,krow), orothe(:,krow), &
              oropic(:,krow), oroval(:,krow), &
              rgrav, zrho_rcp(:,:), pvervel(:,:), zqsi(:,:), &
              !-- INOUT
              zvervx(:,:))
!<<gf
!>>UP ham timers
IF (ltimer) CALL timer_stop(timer_ham_ifdef)
!<<UP
#endif

!>>SF ToDo: move this into get_util_var
  !-- thermodynamic utilities:
  ztmp1_2d(1:kproma,:) = MAX(pqm1(1:kproma,:),0.0_dp)
  ztmp1_2d(1:kproma,:) = 1._dp/(cpd+zcons1*ztmp1_2d(1:kproma,:))
  zlvdcp(1:kproma,:)   = alv*ztmp1_2d(1:kproma,:)
  zlsdcp(1:kproma,:)   = als*ztmp1_2d(1:kproma,:)

  zdv(1:kproma,:) = 2.21_dp/papm1(1:kproma,:)

  ztmp2_2d(1:kproma,:) = 1._dp/MAX(pqm1(1:kproma,:),eps) &
                    + zlsdcp(1:kproma,:)*alv/(rv*ptm1(1:kproma,:)**2) !SF prot. against divide by zero 
  ztmp3_2d(1:kproma,:) = grav*(zlvdcp(1:kproma,:)*rd/rv/ptm1(1:kproma,:)-1._dp)/(rd*ptm1(1:kproma,:))
  ztmp4_2d(1:kproma,:) = 1._dp/(crhoi*als**2/(zkair(1:kproma,:)*rv*ptm1(1:kproma,:)**2) & 
                          + crhoi*rv*ptm1(1:kproma,:)/(zesi(1:kproma,:)*zdv(1:kproma,:)))

  zeta(1:kproma,:) = ztmp2_2d(1:kproma,:)/ztmp3_2d(1:kproma,:)*ztmp4_2d(1:kproma,:)*4._dp*pi*crhoi*cap/zrho(1:kproma,:)
!<<SF ToDo

  !-- cloud ice mass mixing ratio at t:
  zxip1(1:kproma,:) = pxim1(1:kproma,:) + ztmst*pxite(1:kproma,:)
  zxip1(1:kproma,:) = MAX(zxip1(1:kproma,:), 0._dp)

  !>>SF #441 (make the eff radii calculation consistent with the other similar ones)

  !SF conversion of ice mmr from grid-mean kg/kg to in-cloud g m-3
  ztmp2_2d(1:kproma,:) = 1000._dp*zxip1(1:kproma,:)*zrho(1:kproma,:)/MAX(paclc(1:kproma,:),clc_min)

  DO jk=ktdia,klev
     ztmp2_2d(:,jk) = eff_ice_crystal_radius(kbdim, kproma, ztmp2_2d(:,jk), picnc(:,jk))
  ENDDO
  !<<SF #441

  ztmp2_2d(1:kproma,:) = MIN(MAX(ztmp2_2d(1:kproma,:), ceffmin), ceffmax) !SF zrieff in micrometers

  ! Simple param of r/re approximating the Schumann et al. 2011 data:
  zrice(:,:) = effective_2_volmean_radius_param_Schuman_2011( &
                         kbdim, kproma, klev, ztmp2_2d(:,:))

  !-- Threshold vertical velocity for the Wegener-Bergeron-Findeisen process:
  zvervmax(:,:) = threshold_vert_vel(kbdim, kproma, klev, zesw_2d(:,:), &
                                            zesi(:,:), picnc(:,:), &
                                            zrice(:,:), zeta(:,:))

  !-- Final criterion used for the W-B-F process 
  lo2_2d(1:kproma,:) = (0.01_dp*zvervx(1:kproma,:) < zvervmax(1:kproma,:))
!<<DN #286

!>>DN 2013-07-05 #271
     zqlnuccvh(1:kproma,:) = MERGE(cdncact_cv(1:kproma,:,krow), 0._dp, ll1_2d(1:kproma,:))
!<<DN 2013-07-05

!>>UP #782: 1a
  IF (ldetr_liquid) THEN
          ll3_2d(1:kproma,:) = ll2_2d(1:kproma,:)              .AND. &
                               (zxtec(1:kproma,:) > 0._dp)     .AND. &
                               ( (ptm1(1:kproma,:) > tmelt)    .OR.  &
                                 (ptm1(1:kproma,:) > cthomi) &
                                ) !DN #286: modified condition for consistency with stratiform ice / liquid water split
  ELSE !original formulation
          ll3_2d(1:kproma,:) = ll2_2d(1:kproma,:)              .AND. &
                               (zxtec(1:kproma,:) > 0._dp)     .AND. &
                               ( (ptm1(1:kproma,:) > tmelt)    .OR.  &
                                 ( (ptm1(1:kproma,:) > cthomi) .AND. &
                                   (ptm1(1:kproma,:) < tmelt)  .AND. &
                                   (.NOT.lo2_2d(1:kproma,:))         &
                                 ) &
                                ) !DN #286: modified condition for consistency with stratiform ice / liquid water split
  ENDIF
!<<UP #782: 1a


!>>DN 2013-07-05 #271
  ztmp2_2d(1:kproma,:) = cdncact_cv(1:kproma,:,krow)
!<<DN 2013-07-05

  DO jk=ktdia,klev
     DO jl=1,kproma
        jkk = iclbas(jl,jk)
        jkk = MAX(1, jkk) !SF prevents cases where iclbas=0
        ztmp3_2d(jl,jk) = zqlnuccvh(jl,jkk) 
     ENDDO
  ENDDO

  ztmp1_2d(1:kproma,:) = MIN(ztmp2_2d(1:kproma,:),ztmp3_2d(1:kproma,:))
  ztmp1_2d(1:kproma,:) = MAX(0._dp,ztmp1_2d(1:kproma,:))

  zqlnuccv(1:kproma,:)  = MERGE(ztmp1_2d(1:kproma,:), 0._dp, ll3_2d(1:kproma,:))
  zlwc_detr(1:kproma,:) = MERGE(zxtec(1:kproma,:)*ztmst, 0._dp, ll3_2d(1:kproma,:))
  qnuc(1:kproma,:,krow) = qnuc(1:kproma,:,krow) + zdt*zqlnuccv(1:kproma,:)

  ll1_2d(1:kproma,:) = (paclc(1:kproma,:) >= epsec) .AND. &
                       (pxlm1(1:kproma,:)+ztmst*pxlte(1:kproma,:) >= epsec)

  ztmp2_2d(1:kproma,:) = MERGE(paclc(1:kproma,:), zdummy, ll1_2d(1:kproma,:))
  ztmp3_2d(1:kproma,:) = (pxlm1(1:kproma,:)+ztmst*pxlte(1:kproma,:))/ztmp2_2d(1:kproma,:)

  zlwc_strat(1:kproma,:) = MERGE(ztmp3_2d(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))

  ll2_2d(1:kproma,:) = ll1_2d(1:kproma,:) .OR. ll3_2d(1:kproma,:)

  zlwc_tot(1:kproma,:) = zlwc_strat(1:kproma,:) + zlwc_detr(1:kproma,:)
  zlwc_tot(1:kproma,:) = MERGE(zlwc_tot(1:kproma,:), zdummy, ll2_2d(1:kproma,:))

!>>UP, #783.3
  IF (lpr_corr) THEN

  !SF: weighs both stratiform and detrained contribs by their respective liquid water content:
  ztmp2_2d(1:kproma,:) = zcdnc(1:kproma,:)
  !UP: CDNC from detrainment are mixed with those of stratiform
  zcdnc_upd(1:kproma,:) = ( ztmp2_2d(1:kproma,:)*zlwc_strat(1:kproma,:)  &
                      + zqlnuccv(1:kproma,:)*zlwc_detr(1:kproma,:) ) &
                    / zlwc_tot(1:kproma,:)

  !UP: After mixing, CDNC should not be lower than before
  zcdnc_upd(1:kproma,:) = MAX(ztmp2_2d(1:kproma,:),zcdnc(1:kproma,:))
  !DN this ensures that the above weighting will never yield an unphysical decrease #271
  !>>DN: new diags
  zncdetr(1:kproma,:) = zcdnc(1:kproma,:) - ztmp2_2d(1:kproma,:)
  !<<DN
  !UP: true correction of having at least cqtmin
  !UP: when cqtmin is larger than zcdnc, the difference is added, checked by DN
  zcorrcd3(1:kproma,:) = - MIN(zcdnc(1:kproma,:)-cqtmin,0._dp)
  !UP: this is the correction 
  zcdnc_upd(1:kproma,:) = MAX(zcdnc(1:kproma,:),cqtmin)
!<<UP, #783.3
  ELSE

!>>DN: new diags
  !UP: replace zcdnc by zcdnc_upd below? No, left of = must have
  !_update, right of = without
  zncdetr(1:kproma,:)  = zqlnuccv(1:kproma,:) * zlwc_detr(1:kproma,:) &
                         / zlwc_tot(1:kproma,:)
  zcorrcd3(1:kproma,:) = zcdnc(1:kproma,:)*zlwc_strat(1:kproma,:)  &
                         / zlwc_tot(1:kproma,:) - zcdnc(1:kproma,:)
!<<DN: new diags

  !SF: weighs both stratiform and detrained contribs by their respective liquid
  !water content:
  ztmp2_2d(1:kproma,:) = zcdnc(1:kproma,:)
  zcdnc_upd(1:kproma,:) = ( ztmp2_2d(1:kproma,:)*zlwc_strat(1:kproma,:)  &
                      + zqlnuccv(1:kproma,:)*zlwc_detr(1:kproma,:) ) &
                    / zlwc_tot(1:kproma,:)

!>>DN: new diags
  zcorrcd3(1:kproma,:) = zcorrcd3(1:kproma,:) - MIN(zcdnc(1:kproma,:)-cqtmin,0._dp)
!<<DN: new diags

  zcdnc_upd(1:kproma,:) = MAX(ztmp2_2d(1:kproma,:),zcdnc(1:kproma,:),cqtmin)
  !DN this ensures that the above weighting will never yield an unphysical decrease
  !#271
  END IF

!  Ice nucleation

  !>>UP timers
  IF (ltimer) THEN
     CALL timer_start(timer_cmp_icenucl)
  END IF
  !<<UP

  ztmp1_2d(1:kproma,:) = ptm1(1:kproma,:)-tmelt
  ztmp1_2d(1:kproma,:) = MIN(0._dp,ztmp1_2d(1:kproma,:))

  ztmp2_2d(1:kproma,:) = 0.015_dp*ztmp1_2d(1:kproma,:)
  ztmp2_2d(1:kproma,:) = EXP(ztmp2_2d(1:kproma,:))
  ztmp2_2d(1:kproma,:) = 23.2_dp*ztmp2_2d(1:kproma,:)
  ztmp2_2d(1:kproma,:) = MAX(ztmp2_2d(1:kproma,:),1.0_dp)
     
!       Simple param of r/re approximating the Schumann et al. 2011 data

  zrid(:,:) = effective_2_volmean_radius_param_Schuman_2011( &
                         kbdim, kproma, klev, ztmp2_2d(:,:))

!>>UP #782: 1b
  IF (ldetr_liquid) THEN
          ll_cv(1:kproma,:) = (zxtec(1:kproma,:) > 0._dp)     .AND. &
                              (ptm1(1:kproma,:)  < cthomi)  
                              !DN #286: modified condition for consistency with stratiform ice / liquid water split
  ELSE !original formulation
          ll_cv(1:kproma,:) = (zxtec(1:kproma,:) > 0._dp)     .AND. &
                              ( (ptm1(1:kproma,:)  < cthomi)  .OR.  &
                                ( (ptm1(1:kproma,:)  < tmelt) .AND. &
                                  lo2_2d(1:kproma,:) &
                                ) &
                              ) !DN #286: modified condition for consistency with stratiform ice / liquid water split
  ENDIF
!<<UP #782: 1b

  ll1_2d(1:kproma,:) = (paclc(1:kproma,:) > clc_min)

!UL changed to consistenly use plates, 1.5.2012 [based on Lohmann et al, ERL 2008, expression (1)]
!
!>>SF #176 now this is done consistenly between cl. micro and radiation, with a general formula
  ztmp1_2d(1:kproma,:) = conv_effr2mvr*(0.5e-2_dp)**pow_PK*1000._dp/fact_PK*ztmst &
                       * zrho(1:kproma,:)*zxtec(1:kproma,:) &
                       / ( MAX(paclc(1:kproma,:),clc_min) * zrid(1:kproma,:)**pow_PK)
!<<SF
  ztmp1_2d(1:kproma,:) = MERGE(MAX(ztmp1_2d(1:kproma,:), 0._dp), 0._dp, ll1_2d(1:kproma,:))

!>>DN 2013-07-05 #271: removed the former IWC weigthing
  znidetr(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:),  0._dp, ll_cv(1:kproma,:))
!>>UP #782.2
  DO jl=1,kproma
             ! corresponds to SM's convUpdraft
             ztmp1_2d(jl,:) = MAX(0.01, zwcape(jl))
  ENDDO
  IF (ldetr_convicfactor) THEN
             znidetr(1:kproma,:) = MERGE(znidetr(1:kproma,:)/ &
                       MAX(1.0, 1.0 + convicfactor*ztmp1_2d(1:kproma,:)), &
                       znidetr(1:kproma,:), (znidetr(1:kproma,:) > 0._dp))
  ENDIF
!<<UP #782.2
  znidetr(1:kproma,:) = MAX(znidetr(1:kproma,:), cqtmin) !SF actually this could probably be rearranged more
                                                         !smartly with the above... (will be done when
                                                         !cleaning up all min values, see #241)

!>>UP #783.3: need to put in a maximum condition for icnc here
! don't use more znidetr than what brings you to maximum zicnc
     IF (lpr_corr) THEN 
        znidetr(1:kproma,:) = &
                MERGE(MAX(icemax - (zicncq(1:kproma,:)),0._dp),&
                znidetr(1:kproma,:),&
                (zicncq(1:kproma,:) + znidetr(1:kproma,:)  > icemax))
     ENDIF
!<<UP #783.3

  zicncq(1:kproma,:) = zicncq(1:kproma,:) + znidetr(1:kproma,:)
!<<DN 2013-07-05 #271: removed the former IWC weigthing

!>>UP emulator phase in/out cirrus scheme
  IF (lemuphase_nic_cirrus) THEN
      ! nic_cirrus=1 scheme
             ll_ice(1:kproma,:) = (sice(1:kproma,:,krow) > 0._dp) .AND. &
                                  (ptm1(1:kproma,:)      < cthomi)

             ztmp2_2d(1:kproma,:) = 0.75_dp*zrho(1:kproma,:)*sice(1:kproma,:,krow)*zqsi(1:kproma,:) &
                                  / (pi*rhoice*zrid(1:kproma,:)**3)-zicncq(1:kproma,:)
             ztmp3_2d(1:kproma,:) = zascs(1:kproma,:)-zicncq(1:kproma,:)
             ztmp1_2d(1:kproma,:) = MIN(ztmp2_2d(1:kproma,:),ztmp3_2d(1:kproma,:))
             ztmp1_2d(1:kproma,:) = MAX(ztmp1_2d(1:kproma,:),0._dp)

             zninucl_1(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:), 0._dp, ll_ice(1:kproma,:))

             !>>UP #783.3: need to put in a maximum condition for zninucl_1 here
             ! don't use more than what brings you to maximum zicnc
                  IF (lpr_corr) THEN 
                     zninucl_1(1:kproma,:) = &
                             MERGE(MAX(icemax - (zicncq(1:kproma,:)),0._dp),&
                             zninucl_1(1:kproma,:),&
                             (zicncq(1:kproma,:) + zninucl_1(1:kproma,:)  > icemax))
                  ENDIF
             !<<UP #783.3

             zicncq_1(1:kproma,:)  = zicncq(1:kproma,:) + zninucl_1(1:kproma,:)  

      ! nic_cirrus=2 scheme
             zsusatix(1:kproma,:) = sice(1:kproma,:,krow)
             znicex(1:kproma,:)   = 0._dp

        !--- the freezing rate is limited by the number of aerosols available in each mode. 
             DO jk = klev, ktdia, -1 !SFnote: is it really necessary to go upward?

                zapnx(1:kproma,jk,1) = 1.e-6_dp * (zapnx(1:kproma,jk,1) - zicncq(1:kproma,jk)) ![1/cm3]
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
                    CALL xfrzmstr( &
                            !-- IN
                            ll_het, nosize, ztmst, &
                            klev, kbdim, kproma, jk, nfrzmod, &
                            zsusatix(:,jk), zvervx(:,jk), &
                            zapnx(:,jk,:), zaprx(:,jk,:), zapsigx(:,jk,:), &
                            ptm1, tmelt, eps, papm1, cthomi, &
                            !-- OUT
                            zri, znicex )

             ENDDO !jk

        !--- Update ICNC

             ll1_2d(1:kproma,:) = (zri(1:kproma,:) >= epsec) 

             zri(1:kproma,:) = MERGE(zri(1:kproma,:), zrid(1:kproma,:), ll1_2d(1:kproma,:))
             zri(1:kproma,:) = MAX(zri(1:kproma,:), 1.e-6_dp)
             
             ll_ice(1:kproma,:) = (zsusatix(1:kproma,:) > 0._dp) .AND. &
                                  (ptm1(1:kproma,:)     < cthomi)

             ztmp2_2d(1:kproma,:) = 1.e6_dp*zap(1:kproma,:) 
             ztmp3_2d(1:kproma,:) = znicex(1:kproma,:)
             ztmp1_2d(1:kproma,:) = MIN(ztmp2_2d(1:kproma,:),ztmp3_2d(1:kproma,:)) 
             ztmp1_2d(1:kproma,:) = MAX(ztmp1_2d(1:kproma,:), 0._dp)

             zninucl_2(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:), 0._dp,ll_ice(1:kproma,:))
             !>>UP #783.3: need to put in a maximum condition for zninucl_2 here
             ! don't use more than what brings you to maximum zicnc
                  IF (lpr_corr) THEN 
                     zninucl_2(1:kproma,:) = &
                             MERGE(MAX(icemax - (zicncq(1:kproma,:)),0._dp),&
                             zninucl_2(1:kproma,:),&
                             (zicncq(1:kproma,:) + zninucl_2(1:kproma,:)  > icemax))
                  ENDIF
             !<<UP #783.3

             zicncq_2(1:kproma,:)  = zicncq(1:kproma,:) + zninucl_2(1:kproma,:)

      !combine
             zicncq(1:kproma,:)    = eta_emu_nic_cirrus*zicncq_1(1:kproma,:) + (1-eta_emu_nic_cirrus)*zicncq_2(1:kproma,:)
             ! make a new zninucl that combines the two schemes as well (so that it can be used safely after)
             zninucl(1:kproma,:)   = eta_emu_nic_cirrus*zninucl_1(1:kproma,:) + (1-eta_emu_nic_cirrus)*zninucl_2(1:kproma,:)
!UP #783.3: maximum condition for icnc to limit zninucl already implemented for
!the single zninucl_1,2 above

      !continue nic_cirrus=2 scheme
        !---calculate the deposition rate taking ventilation (zfre) into account

             ll_ice(1:kproma,:) = (pxim1(1:kproma,:) > 0._dp) .AND. &
                                  (paclc(1:kproma,:) > clc_min)

             ztmp1_2d(1:kproma,:) = MAX(zicncq(1:kproma,:), icemin)
             zicncq(1:kproma,:)   = MERGE(ztmp1_2d(1:kproma,:), zicncq(1:kproma,:), ll_ice(1:kproma,:)) 
      !UP note: I left this after the combination for now, because it looks like a safeguard.
      ! It could well be that the combination should be after this, though. TODO
!UP #783.3? No, there's nothing to be done here since zicncq is only used to
!compute rates

             ztmp1_2d(1:kproma,:) = zrho(1:kproma,:)*pxim1(1:kproma,:) & 
                                  / (zicncq(1:kproma,:) * MAX(paclc(1:kproma,:),clc_min))
             ztmp1_2d(1:kproma,:) = MAX(ztmp1_2d(1:kproma,:), mi)

             zmmean(1:kproma,:)   = MERGE(ztmp1_2d(1:kproma,:), mi, ll_ice(1:kproma,:))

             ll1_2d(1:kproma,:) = (zmmean(1:kproma,:) < ri_vol_mean_1 )
             ll2_2d(1:kproma,:) = ( .NOT. ll1_2d(1:kproma,:) ) .AND. &
                                  (zmmean(1:kproma,:) < ri_vol_mean_2 )

             zalfased(1:kproma,:) = MERGE(alfased_1, alfased_2, ll1_2d(1:kproma,:))
             zalfased(1:kproma,:) = MERGE(alfased_3, zalfased(1:kproma,:), ll2_2d(1:kproma,:))

             zbetased(1:kproma,:) = MERGE(betased_1, betased_2, ll1_2d(1:kproma,:))
             zbetased(1:kproma,:) = MERGE(betased_3, zbetased(1:kproma,:), ll2_2d(1:kproma,:))

             zxifallmc(1:kproma,:) = fall*zalfased(1:kproma,:) &
                                   * (zmmean(1:kproma,:)**zbetased(1:kproma,:))*zaaa(1:kproma,:)

             ztmp1_2d(1:kproma,:) = pqm1(1:kproma,:) - zqsi(1:kproma,:)

             DO jk = ktdia, klev
                DO jl = 1,kproma
                   zgtp    = 1._dp/(zrho(jl,jk)*zastbsti(jl,jk))
                   zvth    = SQRT( 8._dp*kb*ptm1(jl,jk) / (pi*xmw) )
                   zb2     = 0.25_dp * alpha * zvth   / zdv(jl,jk)
                   zfuchs  = 1._dp/(1._dp+zb2*zri(jl,jk))
                   zre     = 2._dp*zrho(jl,jk)*zri(jl,jk)*zxifallmc(jl,jk)/zviscos(jl,jk)
                   zfre    = 1._dp + 0.229_dp*SQRT(zre)
                   zfre    = MERGE(zfre, 1._dp, ll_ice(jl,jk))

                   zqinucl(jl,jk)  = 4._dp*pi*zri(jl,jk)*zsusatix(jl,jk)*zicncq(jl,jk) &
                                   *zfre*zgtp*zfuchs*alpha*ztmst
                   zqinucl(jl,jk) = MIN(zqinucl(jl,jk),ztmp1_2d(jl,jk))
                   zqinucl(jl,jk) = MAX(zqinucl(jl,jk),-pxim1(jl,jk))

                ENDDO !jl
             ENDDO !jk


#ifdef HAMMOZ
        !>>UP ham timers
        IF (ltimer) CALL timer_start(timer_ham_ifdef)
        !<<UP
             !--- update cloud cover with orographic cirrus contribs if relevant (#65)
             IF(lorocirrus) CALL orocirrus_cc( &
                                    !-- IN
                                    kproma, kbdim, klev, ktdia, krow, &
                                    zqinucl(1:kproma,:), zsusatix(1:kproma,:), &
                                    !-- INOUT
                                    paclc(1:kproma,:))
        !>>UP ham timers
        IF (ltimer) CALL timer_stop(timer_ham_ifdef)
        !<<UP
#endif


  ELSE
!<<UP
!SF modified the cond. statement
  IF ( nic_cirrus == 1 ) THEN ! Use ICNC scheme after Lohmann, JAS, 2002
        
     ll_ice(1:kproma,:) = (sice(1:kproma,:,krow) > 0._dp) .AND. &
                          (ptm1(1:kproma,:)      < cthomi)

     ztmp2_2d(1:kproma,:) = 0.75_dp*zrho(1:kproma,:)*sice(1:kproma,:,krow)*zqsi(1:kproma,:) &
                          / (pi*rhoice*zrid(1:kproma,:)**3)-zicncq(1:kproma,:)
     ztmp3_2d(1:kproma,:) = zascs(1:kproma,:)-zicncq(1:kproma,:)
     ztmp1_2d(1:kproma,:) = MIN(ztmp2_2d(1:kproma,:),ztmp3_2d(1:kproma,:))
     ztmp1_2d(1:kproma,:) = MAX(ztmp1_2d(1:kproma,:),0._dp)

     !>>UP #801
     IF (lemuphase_icnucl) THEN
             zninucl(1:kproma,:) = MERGE(eta_emu_icnucl*ztmp1_2d(1:kproma,:), 0._dp, ll_ice(1:kproma,:))
     ELSE
     !UP original line
             zninucl(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:), 0._dp, ll_ice(1:kproma,:))
     ENDIF
     !<<UP #801

     !>>UP #783.3: need to put in a maximum condition for zninucl here
     ! don't use more than what brings you to maximum zicnc
     IF (lpr_corr) THEN 
             zninucl(1:kproma,:) = &
              MERGE(MAX(icemax - (zicncq(1:kproma,:)),0._dp),&
              zninucl(1:kproma,:),&
              (zicncq(1:kproma,:) + zninucl(1:kproma,:)  > icemax))
     ENDIF
     !<<UP #783.3

     zicncq(1:kproma,:)  = zicncq(1:kproma,:) + zninucl(1:kproma,:)  

  ELSE IF ( nic_cirrus == 2 ) THEN !--- Use Bernd's cirrus scheme 

!
!--- Kaercher and Lohmann, JGR, 2002b; Lohmann et al., JGR, 2004)
!--- This implies that supersaturation with respect to ice is allowed, thus the depositional
!--- growth equation needs to be solved
!

     zsusatix(1:kproma,:) = sice(1:kproma,:,krow)
     znicex(1:kproma,:)   = 0._dp

!--- the freezing rate is limited by the number of aerosols available in each mode. 
     DO jk = klev, ktdia, -1 !SFnote: is it really necessary to go upward?

        zapnx(1:kproma,jk,1) = 1.e-6_dp * (zapnx(1:kproma,jk,1) - zicncq(1:kproma,jk)) ![1/cm3]
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

            CALL xfrzmstr( &
                    !-- IN
                    ll_het, nosize, ztmst, &
                    klev, kbdim, kproma, jk, nfrzmod, &
                    zsusatix(:,jk), zvervx(:,jk), &
                    zapnx(:,jk,:), zaprx(:,jk,:), zapsigx(:,jk,:), &
                    ptm1, tmelt, eps, papm1, cthomi, &
                    !-- OUT
                    zri, znicex )

         ENDDO !jk

!--- Update ICNC

     ll1_2d(1:kproma,:) = (zri(1:kproma,:) >= epsec) 

     zri(1:kproma,:) = MERGE(zri(1:kproma,:), zrid(1:kproma,:), ll1_2d(1:kproma,:))
     zri(1:kproma,:) = MAX(zri(1:kproma,:), 1.e-6_dp)
     
     ll_ice(1:kproma,:) = (zsusatix(1:kproma,:) > 0._dp) .AND. &
                          (ptm1(1:kproma,:)     < cthomi)

     ztmp2_2d(1:kproma,:) = 1.e6_dp*zap(1:kproma,:) 
     ztmp3_2d(1:kproma,:) = znicex(1:kproma,:)
     ztmp1_2d(1:kproma,:) = MIN(ztmp2_2d(1:kproma,:),ztmp3_2d(1:kproma,:)) 
     ztmp1_2d(1:kproma,:) = MAX(ztmp1_2d(1:kproma,:), 0._dp)

     !>>UP #821
     IF (lcmpsimpl_prescr) THEN
        SELECT CASE(ncmpsimpl_prescr_icnucl)
            CASE(0)
                !UP original line
                zninucl(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:), 0._dp,ll_ice(1:kproma,:))
            CASE(1,2)
                zninucl(1:kproma,:) = MERGE(zcmpsimpl_icnucl(1:kproma,:), 0._dp, ll_ice(1:kproma,:))
        END SELECT
     !<<UP #821
     ELSE !Only makes sense to phase if we don't prescribe
             !>>UP #801
             IF (lemuphase_icnucl) THEN
                     zninucl(1:kproma,:) = MERGE(eta_emu_icnucl*ztmp1_2d(1:kproma,:), 0._dp,ll_ice(1:kproma,:))
             ELSE
             !UP original line
                     zninucl(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:), 0._dp,ll_ice(1:kproma,:))
             ENDIF
             !<<UP #801
     ENDIF
     !>>UP #783.3: need to put in a maximum condition for zninucl here
     ! don't use more than what brings you to maximum zicnc
          IF (lpr_corr) THEN 
             zninucl(1:kproma,:) = &
                     MERGE(MAX(icemax - (zicncq(1:kproma,:)),0._dp),&
                     zninucl(1:kproma,:),&
                     (zicncq(1:kproma,:) + zninucl(1:kproma,:)  > icemax))
          ENDIF
     !<<UP #783.3
     !>>UP #821
     !This diagnostic needs to be here, after the maximum condition
     diag_delta_icnucl(1:kproma,:,krow) = diag_delta_icnucl(1:kproma,:,krow) + &
                                          zninucl(1:kproma,:)
     diag_cntr_icnucl(1:kproma,:,krow) = diag_cntr_icnucl(1:kproma,:,krow) + &
                                         MERGE(1._dp, 0._dp, ll_ice(1:kproma,:))
     !<<UP #821

     zicncq(1:kproma,:)  = zicncq(1:kproma,:) + zninucl(1:kproma,:)

!---calculate the deposition rate taking ventilation (zfre) into account

     ll_ice(1:kproma,:) = (pxim1(1:kproma,:) > 0._dp) .AND. &
                          (paclc(1:kproma,:) > clc_min)

     ztmp1_2d(1:kproma,:) = MAX(zicncq(1:kproma,:), icemin)
!UPtoDN: above is a correction that is not included in the diagnostics?
     zicncq(1:kproma,:)   = MERGE(ztmp1_2d(1:kproma,:), zicncq(1:kproma,:), ll_ice(1:kproma,:)) 
!UP #783.3: Again, nothing to be done here for zicncq

     ztmp1_2d(1:kproma,:) = zrho(1:kproma,:)*pxim1(1:kproma,:) & 
                          / (zicncq(1:kproma,:) * MAX(paclc(1:kproma,:),clc_min))
     ztmp1_2d(1:kproma,:) = MAX(ztmp1_2d(1:kproma,:), mi)

     zmmean(1:kproma,:)   = MERGE(ztmp1_2d(1:kproma,:), mi, ll_ice(1:kproma,:))

     ll1_2d(1:kproma,:) = (zmmean(1:kproma,:) < ri_vol_mean_1 )
     ll2_2d(1:kproma,:) = ( .NOT. ll1_2d(1:kproma,:) ) .AND. &
                          (zmmean(1:kproma,:) < ri_vol_mean_2 )

     zalfased(1:kproma,:) = MERGE(alfased_1, alfased_2, ll1_2d(1:kproma,:))
     zalfased(1:kproma,:) = MERGE(alfased_3, zalfased(1:kproma,:), ll2_2d(1:kproma,:))

     zbetased(1:kproma,:) = MERGE(betased_1, betased_2, ll1_2d(1:kproma,:))
     zbetased(1:kproma,:) = MERGE(betased_3, zbetased(1:kproma,:), ll2_2d(1:kproma,:))

     zxifallmc(1:kproma,:) = fall*zalfased(1:kproma,:) &
                           * (zmmean(1:kproma,:)**zbetased(1:kproma,:))*zaaa(1:kproma,:)

     ztmp1_2d(1:kproma,:) = pqm1(1:kproma,:) - zqsi(1:kproma,:)

     DO jk = ktdia, klev
        DO jl = 1,kproma
           zgtp    = 1._dp/(zrho(jl,jk)*zastbsti(jl,jk))
           zvth    = SQRT( 8._dp*kb*ptm1(jl,jk) / (pi*xmw) )
           zb2     = 0.25_dp * alpha * zvth   / zdv(jl,jk)
           zfuchs  = 1._dp/(1._dp+zb2*zri(jl,jk))
           zre     = 2._dp*zrho(jl,jk)*zri(jl,jk)*zxifallmc(jl,jk)/zviscos(jl,jk)
           zfre    = 1._dp + 0.229_dp*SQRT(zre)
           zfre    = MERGE(zfre, 1._dp, ll_ice(jl,jk))

           !>>UP #801 phasing deposition
           IF (lemuphase_dep) THEN
                   zqinucl(jl,jk)  = eta_emu_dep*4._dp*pi*zri(jl,jk)*zsusatix(jl,jk) &
                                   *zicncq(jl,jk)*zfre*zgtp*zfuchs*alpha*ztmst
           ELSE
           ! original formulation
                   zqinucl(jl,jk)  = 4._dp*pi*zri(jl,jk)*zsusatix(jl,jk) &
                                   *zicncq(jl,jk)*zfre*zgtp*zfuchs*alpha*ztmst
           ENDIF
           !<<UP
           !UP note: don't deposit more than the oversaturation allows
           zqinucl(jl,jk) = MIN(zqinucl(jl,jk),ztmp1_2d(jl,jk))
           !UP note: if the changes due to deposition are negative, at most all
           !the ice from the last timestep is taken away
           zqinucl(jl,jk) = MAX(zqinucl(jl,jk),-pxim1(jl,jk))

        ENDDO !jl
     ENDDO !jk


#ifdef HAMMOZ
!>>UP ham timers
IF (ltimer) CALL timer_start(timer_ham_ifdef)
!<<UP
     !--- update cloud cover with orographic cirrus contribs if relevant (#65)
     IF(lorocirrus) CALL orocirrus_cc( &
                            !-- IN
                            kproma, kbdim, klev, ktdia, krow, &
                            zqinucl(1:kproma,:), zsusatix(1:kproma,:), &
                            !-- INOUT
                            paclc(1:kproma,:))
!>>UP ham timers
IF (ltimer) CALL timer_stop(timer_ham_ifdef)
!<<UP
#endif
 
  ENDIF   !which cirrus scheme (nic_cirrus)
 ENDIF ! lemuphase_nic_cirrus
 
  !corinna: set zcdnc and picnc to minium now if nucleation is not strong enough

  ll4_2d(1:kproma,:) = ( paclc(1:kproma,:) >= epsec ) .AND. &
                       ( ptm1(1:kproma,:)  >  cthomi )

!>>DN: new diags
  ztmp1_2d(1:kproma,:)     = MAX(zcdnc_min(1:kproma,:)-zcdnc(1:kproma,:),0._dp)
!UP, #783.3: zcorrcd4 is a true correction term (zcdnc set to at least minimum 
! where you have a cloud, and to 0 outside of cloud)
! However, we do not have cloud water in here yet, so we need to correct corrcd4
! later, when we do have cloud water
  zcorrcd4(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:),0._dp,ll4_2d(1:kproma,:))
!<<DN: new diags

  ztmp1_2d(1:kproma,:) = MAX(zcdnc(1:kproma,:),zcdnc_min(1:kproma,:))
  zcdnc_upd(1:kproma,:) = MERGE(ztmp1_2d(1:kproma,:), zcdnc(1:kproma,:), ll4_2d(1:kproma,:))

  ll4_2d(1:kproma,:) = ( paclc(1:kproma,:) >= epsec ) .AND. &
                       ( ptm1(1:kproma,:)  <  tmelt  )

  ztmp1_2d(1:kproma,:) = MAX(zicncq(1:kproma,:), icemin)
  zicncq(1:kproma,:)   = MERGE(ztmp1_2d(1:kproma,:), zicncq(1:kproma,:), ll4_2d(1:kproma,:))
!UP question, #783: why did David not diagnose the correction of zicnc here?
!It's zicncq, which is only used for computing rates
!UP #783.3 max for icncq here not needed since zicncq is only used to compute
!rates

  !UPcomment: I guess ice nucleation goes on until here for now, there is more further down
  !>>UP timers
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_1cbasetop)
     CALL timer_stop(timer_cmp_icenucl)
  END IF
  !<<UP


!--- End included for CDNC/IC scheme -----------------------------------

column_processes:  DO jk=ktdia,klev

!     ------------------------------------------------------------------
!
!       2.    Set to zero local tendencies (increments)
!
     !>>UP timers
     IF (ltimer) THEN
        CALL timer_start(timer_cmp_2zerolocaltend)
     END IF
     !<<UP

     zcnd(1:kproma)       = 0.0_dp
     zdep(1:kproma)       = 0.0_dp
     zimlt(1:kproma)      = 0.0_dp
     zevp(1:kproma)       = 0.0_dp
     zsub(1:kproma)       = 0.0_dp
     zximlt(1:kproma)     = 0.0_dp
     zxisub(1:kproma)     = 0.0_dp
     zsmlt(1:kproma)      = 0.0_dp
     zsacl(1:kproma)      = 0.0_dp
     zgenti(1:kproma)     = 0.0_dp
     zgentl(1:kproma)     = 0.0_dp
     zxievap(1:kproma)    = 0.0_dp
     zxlevap(1:kproma)    = 0.0_dp
     zvartg(1:kproma)     = 0.0_dp
     zconvvar(1:kproma)   = 0.0_dp
     zconvskew(1:kproma)  = 0.0_dp
     zturbvar(1:kproma)   = 0.0_dp
     zturbskew(1:kproma)  = 0.0_dp
     zmicroskew(1:kproma) = 0.0_dp
     zxvarte(1:kproma)    = 0.0_dp
     zxskewte(1:kproma)   = 0.0_dp

     !>>UP debug issue #767
     zcorrcd10(1:kproma) = 0._dp
     !<<UP debug

!davidn
     zfrzn(1:kproma)  = 0._dp
     zevabf(1:kproma) = 0._dp
     zlwcbf(1:kproma) = 0._dp
     zmltn(1:kproma)  = 0._dp
     zmltn2(1:kproma) = 0._dp
!davidn
!>>UP #783
     zfrznhet(1:kproma)  = 0._dp
     zfrznhom(1:kproma)  = 0._dp
!<<UP #783

     !>>UP timers, #777
     IF (ltimer) THEN
        CALL timer_stop(timer_cmp_2zerolocaltend)
     END IF
     !<<UP

!     ------------------------------------------------------------------
!
!       3.   Modification of incoming precipitation fluxes by
!            melting, sublimation and evaporation
!
     IF (jk > 1) THEN

!       3.1   Melting of snow and ice

        !>>UP timers
        IF (ltimer) THEN
           CALL timer_start(timer_cmp_3melting)
        END IF
        !<<UP

        ll_mlt(1:kproma) = (ptm1(1:kproma,jk) > tmelt)

        IF (lctrl_mlt_snow_ice .AND. ANY(ll_mlt(1:kproma))) &

           CALL melting_snow_and_ice( &

                   !-- IN
                   kbdim, kproma, &
                   ll_mlt(:), &
                   ptm1(:,jk), pxim1(:,jk), &
                   zdp(:,jk), zicncq(:,jk), zlsdcp(:,jk), &
                   zlvdcp(:,jk), zcdnc(:,jk), picnc(:,jk), &     
                   jk, krow, zcmpsimpl_mlt(:,jk), & !UP #821
                   !-- INOUT
                   zicnc_upd(:,jk), qmel(:,jk,krow), zcdnc_upd(:,jk), &
                   prsfl(:), pssfl(:), zxiflux(:), zxifluxn(:), & 
                   pxite(:,jk), &
                   !-- OUT
                   zimlt(:), zsmlt(:), zximlt(:), zmltn(:), zmltn2(:))!davidn

!UP #783.3: need to put in a maximum condition for icnc here? No, as picnc_upd
!is only lowered inside the routine; DN agrees
        !>>UP timers
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_3melting)
           CALL timer_start(timer_cmp_3subandevap)
        END IF
        !<<UP

!>>DN: burden
        dmlts_atm(1:kproma,krow)   = dmlts_atm(1:kproma,krow)   + zdt*zsmlt(1:kproma)*zdp(1:kproma,jk)/grav
        qmlts_atm(1:kproma,jk,krow)= qmlts_atm(1:kproma,jk,krow)+ zdt*zsmlt(1:kproma)
!<<DN: burden

!       3.2   Sublimation of snow (zsub) and ice (zxisub) following (Lin et al., 1983)
!       and
!       3.3   Evaporation of rain (zevp following Rotstayn, 1997)

        ll_precip(1:kproma)      = (zclcpre(1:kproma) > 0._dp)
        ll_falling_ice(1:kproma) = (zclcfi(1:kproma) > 0._dp)

        IF (lctrl_subl_evap .AND. (ANY(ll_precip(1:kproma) .OR. ANY(ll_falling_ice(1:kproma))))) &

           CALL sublimation_snow_and_ice_evaporation_rain( &

                   !-- IN 
                   kbdim, kproma, &
                   ll_precip(:), ll_falling_ice(:), &
                   pqm1(:,jk), ptm1(:,jk), zclcpre(:), &
                   zdp(:,jk), zdpg(:,jk), zicesub(:,jk), &
                   zlsdcp(:,jk), zqrho(:,jk), zqsi(:,jk), &
                   zrho_rcp(:,jk), pssfl(:), zrho(:,jk), &    
                   zqsw(:,jk), prsfl(:), zsusatw_evap(:,jk), &
                   zastbstw(:,jk), zclcfi(:), &
                   jk, krow, & !UP #821
                   zcmpsimpl_subfi(:,jk), zcmpsimpl_subfs(:,jk), & !UP #821
                   !-- INOUT
                   zxiflux(:), zxifluxn(:), &
                   !-- OUT
                   zxisub(:), zsub(:), zevp(:) )

        !>>UP timers
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_3subandevap)
        END IF
        !<<UP

      ENDIF !SF end jk > 1
!
!     --- End included for CDNC/IC scheme -----------------------------------
!
!     ------------------------------------------------------------------
!       4.    Sedimentation of cloud ice from grid-mean values.
!             Updating the tendency 'pxite' to include sedimentation.
!             At jk=klev, the sedimentation sink is balanced by
!             precipitation at the surface (through 'zzdrs', see 7.3).
!             Finally: In-cloud cloud water/ice.
!
        !>>UP timers
        IF (ltimer) THEN
           CALL timer_start(timer_cmp_4cicesed)
        END IF
        !<<UP

        zxip1(1:kproma,jk) = pxim1(1:kproma,jk) + ztmst*pxite(1:kproma,jk)
        zxip1(1:kproma,jk) = MAX(zxip1(1:kproma,jk), EPSILON(1._dp))

!>>DN: new diags
        ztmp1(1:kproma) = picnc(1:kproma,jk)
!<<DN: new diags
!>>UP #797
        ztmp2(1:kproma) = zxip1(1:kproma,jk)
!<<UP #797
        IF (lctrl_sed_ice) THEN

           !>>UP #797
           IF (jk>1) THEN
                   CALL sedimentation_ice( &
                           !-- IN
                           kbdim, kproma, &
                           paclc(:,jk), zaaa(:,jk), zdp(:,jk), &
                           zrho(:,jk), zrho_rcp(:,jk), &
                           picnc(:,jk), &
                           !>>UP #797
                           zicnc_diff(:), &
                           zxip1_diff(:), &
                           pqm1(:,jk), zqsi(:,jk), &
                           !<<UP #797
                           !-- INOUT
                           zxip1(:,jk), &
                           zicnc_upd(:,jk), &
                           zxiflux(:), zxifluxn(:), zclcfi(:), &
                           zxisub(:), & !UP #797
                           !-- OUT
                           zmrateps(:,jk) )
           ELSE
                ! below is original, UP #797
                   CALL sedimentation_ice( &
                           !-- IN
                           kbdim, kproma, &
                           paclc(:,jk), zaaa(:,jk), zdp(:,jk), &
                           zrho(:,jk), zrho_rcp(:,jk), &
                           picnc(:,jk), &
                           !>>UP #797, in else case not used
                           picnc(:,jk), & 
                           zxip1(:,jk), &
                           pqm1(:,jk), zqsi(:,jk), &
                           !<<UP #797
                           !-- INOUT
                           zxip1(:,jk), &
                           zicnc_upd(:,jk), &
                           zxiflux(:), zxifluxn(:), zclcfi(:), &
                           zxisub(:), & !UP #797, in the else case not used
                           !-- OUT
                           zmrateps(:,jk) )
           ENDIF !jk>1
           !<<UP #797
        ELSE
            zmrateps(1:kproma,jk)  = 0._dp
        ENDIF
!UP #783.3: need to put in a maximum condition for icnc here? Looks like it
!might get raised in some cases, so better be safe
! But if I do it here, it does not register with any diagnostic, so rather do
! that check below
! DN: I wouldn't add a maximum icnc condition for sedimentation and rather keep
! it as it is with the check below

!>>UP #797
! safe the decreases in ice number and mass to input into the sedimentation
! routine in the next round as a restriction for how large any increase can be.
! The increase should be restricted to the size of the last available decrease
! in that grid box (going down the levels)
        zicnc_diff(1:kproma) = MERGE(picnc(1:kproma,jk) - ztmp1(1:kproma), &
                zicnc_diff(1:kproma), (picnc(1:kproma,jk) - ztmp1(1:kproma)) < 0)
! For mass we need to compute in-cloud values as this is where the restriction
! holds
! Need to use MAX with clc_min here to not divide by 0
! Could this be a problem by creating large numbers for 0 cloud cover? I think
! here that's not an issue because it only counts when the grid box is losing
! ice through sedimentation, which means that there should be cloud cover
        zxip1_diff(1:kproma) = MERGE( &
                (zxip1(1:kproma,jk) - ztmp2(1:kproma))/MAX(paclc(:,jk), clc_min), &
                zxip1_diff(1:kproma), (zxip1(1:kproma,jk) - ztmp2(1:kproma)) < 0)
!<<UP #797

!>>SF Kasja diags
        zsedi(1:kproma)  = pxim1(1:kproma,jk) + ztmst*pxite(1:kproma,jk) - zxip1(1:kproma,jk)
        zsedin(1:kproma) = (picnc(1:kproma,jk) - ztmp1(1:kproma))
!<<SF Kasja diags
        pxite(1:kproma,jk) = ztmst_rcp*(zxip1(1:kproma,jk)-pxim1(1:kproma,jk))

        !>>UP write out number conc. of newly nucleated IC for tendency analysis
        ! nucleation happens further above, but zninucl is not mentioned in between, 
        ! so I think it's okay to write it out here
        IF (lcmptend) THEN
                tend_inucl(1:kproma,jk,krow) = tend_inucl(1:kproma,jk,krow) + zninucl(1:kproma,jk) * zdtime
        ENDIF
        !<<UP

        !SF-- sums up icnc-after-sedimentation, detrainment and nucleation
        zicnc_upd(1:kproma,jk) = picnc(1:kproma,jk) + znidetr(1:kproma,jk) + zninucl(1:kproma,jk)
!>>DN: new diags
        ztmp1(1:kproma) = picnc(1:kproma,jk) ! UP: use picnc here, not zicnc_upd
!<<DN: new diags
!UP #783.3: separate below check into multiple right after icnc is changed, but
!leave this one here as last resort
        zicnc_upd(1:kproma,jk) = MIN(picnc(1:kproma,jk), icemax)
!UP #783.3: this is a true correction term and can remain here
        zicnc_upd(1:kproma,jk) = MAX(picnc(1:kproma,jk), icemin)
!>>DN: new diags
        zcorric7(1:kproma) = picnc(1:kproma,jk)-ztmp1(1:kproma)
!<<DN: new diags
     
        !>>SF correct for inconsistencies relative to phases and temperature
        ll1(1:kproma)    = (ptm1(1:kproma,jk)  > tmelt) .AND. &
                           (picnc(1:kproma,jk) > icemin)
        zic2cd(1:kproma) = MERGE(picnc(1:kproma,jk), 0._dp, ll1(1:kproma))
      
!>>DN: new diags
        ztmp1(1:kproma) = picnc(1:kproma,jk) 
!<<DN: new diags
        zicnc_upd(1:kproma,jk) = picnc(1:kproma,jk) - zic2cd(1:kproma)
        zcdnc_upd(1:kproma,jk) = zcdnc(1:kproma,jk) + zic2cd(1:kproma)
      
        zicnc_upd(1:kproma,jk) = MAX(picnc(1:kproma,jk),icemin)
!>>DN: new diags
        zcorrcd6(1:kproma) = zic2cd(1:kproma)
        zcorric4(1:kproma) = (picnc(1:kproma,jk)-ztmp1(1:kproma))
!<<DN: new diags
        !<<SF
    
        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_4cicesed)
           CALL timer_start(timer_cmp_4incloudwi)
        END IF
        !<<UP

!             In-cloud water/ice calculated from respective grid-means,
!             partial cloud cover, advective/diffusive tendencies,
!             detrained cloud water/ice and ice sedimentation.
!             In-cloud values are required for cloud microphysics.
!
!       THIS PART NEEDS TO BE COMMENTED BY ERICH/ADRIAN
!
!SF ToDo: this part needs refactoring. It actually does not belong to section 4. anymore
!        (sounds like prep for section 5.)

        ll_cc(1:kproma) = (paclc(1:kproma,jk) > clc_min) !DN / SF #473: (lower clc limit)
        
!>>SF #176 now this is done consistenly between cl. micro and radiation, with a general formula

        !SF conversion of ice mmr from grid-mean kg/kg to in-cloud g.m-3:
        ztmp1(1:kproma) = 1000._dp*zxip1(1:kproma,jk)*zrho(1:kproma,jk)/MAX(paclc(1:kproma,jk),clc_min)
   
        ztmp1(:) = eff_ice_crystal_radius(kbdim, kproma, ztmp1(:), picnc(:,jk))
!<<SF
        ztmp1(1:kproma) = MIN(MAX(ztmp1(1:kproma), ceffmin), ceffmax) !SF zrieff in micrometers

!       Simple param of r/re approximating the Schumann et al. 2011 data
        zrice(:,jk) = effective_2_volmean_radius_param_Schuman_2011( &
                                kbdim, kproma, ztmp1(:))

        zvervmax(:,jk) = threshold_vert_vel(kbdim, kproma, zesw_2d(:,jk), &
                                   zesi(:,jk), picnc(:,jk), &
                                   zrice(:,jk), zeta(:,jk))

        lo2(1:kproma)  = (ptm1(1:kproma,jk) < cthomi) .OR. &
                             (ptm1(1:kproma,jk) < tmelt .AND. &
                              0.01_dp*zvervx(1:kproma,jk) < zvervmax(1:kproma,jk) &
                             )         

!>>SF DN #368: energy consistency between conv and strat schemes
        IF (lconv) THEN  !csld#455
           IF (ldetr_liquid) THEN !UP, #782
              ll1(1:kproma) = (ztconv(1:kproma,jk) <= tmelt .AND. ptm1(1:kproma,jk) > cthomi)
           ELSE
              ll1(1:kproma) = (ztconv(1:kproma,jk) <= tmelt) &
                   .AND. .NOT. lo2(1:kproma)
           ENDIF

           ztmp1(1:kproma)   = ptte(1:kproma,jk)-(als-alv)*zxtec(1:kproma,jk)/cpd
           ptte(1:kproma,jk) = MERGE(ztmp1(1:kproma), ptte(1:kproma,jk), ll1(1:kproma))
        ENDIF !csld #455
!<<SF DN #368
!davidn
!>>UP, #782
        IF (ldetr_liquid) THEN !UP, #782
           WHERE((ztconv(1:kproma,jk)-tmelt)<=0._dp .AND. ptm1(1:kproma,jk) > cthomi)
              dmlt_conv(1:kproma,krow) = dmlt_conv(1:kproma,krow)+zdtime*&
                   (als-alv)*zxtec(1:kproma,jk)*zdpg(1:kproma,jk)
              qmlt_conv(1:kproma,jk,krow) = qmlt_conv(1:kproma,jk,krow)+zdtime*&
                   zxtec(1:kproma,jk)
           END WHERE
        ELSE
!<<UP, #782
           WHERE((ztconv(1:kproma,jk)-tmelt)<=0._dp.AND.(.NOT.lo2(1:kproma)))
              dmlt_conv(1:kproma,krow) = dmlt_conv(1:kproma,krow)+zdtime*&
                   (als-alv)*zxtec(1:kproma,jk)*zdpg(1:kproma,jk)
              qmlt_conv(1:kproma,jk,krow) = qmlt_conv(1:kproma,jk,krow)+zdtime*&
                   zxtec(1:kproma,jk)
           END WHERE
        ENDIF
!davidn

!>>UP #782: 1c
  IF (ldetr_liquid) THEN
        zxite(1:kproma) = MERGE(zxtec(1:kproma,jk), 0._dp, (ptm1(1:kproma,jk) < cthomi))
        zxlte(1:kproma) = MERGE(0._dp, zxtec(1:kproma,jk), (ptm1(1:kproma,jk) < cthomi))

        zxite2(1:kproma) = MERGE(zxite(1:kproma), 0._dp, (ptm1(1:kproma,jk) < cthomi))
        zxlte2(1:kproma) = MERGE(0._dp, zxlte(1:kproma), (ptm1(1:kproma,jk) < cthomi))
  ELSE !original formulation
        zxite(1:kproma) = MERGE(zxtec(1:kproma,jk), 0._dp, lo2(1:kproma))
        zxlte(1:kproma) = MERGE(0._dp, zxtec(1:kproma,jk), lo2(1:kproma))

        zxite2(1:kproma) = MERGE(zxite(1:kproma), 0._dp, lo2(1:kproma))
        zxlte2(1:kproma) = MERGE(0._dp, zxlte(1:kproma), lo2(1:kproma))
  ENDIF
!<<UP #782: 1c

        zxidt(1:kproma) = ztmst*(pxite(1:kproma,jk)+zxite2(1:kproma))
        zxldt(1:kproma) = ztmst*(pxlte(1:kproma,jk)+zxlte2(1:kproma)) + zximlt(1:kproma) + zimlt(1:kproma) 

        ll1(1:kproma) = (zxidt(1:kproma) > 0._dp) !fix Remo Dietlicher

!>>DN: new diags
        ztmp2(1:kproma) = pxite(1:kproma,jk)
!<<DN: new diags
!davidn code modif        
        ll1(1:kproma) = (zxidt(1:kproma) > 0._dp) !fix Remo
!davidn code modif

        ztmp1(1:kproma)    = -( ztmst_rcp*pxim1(1:kproma,jk) + zxite2(1:kproma) )
        ztmp1(1:kproma)    = MAX(pxite(1:kproma,jk), ztmp1(1:kproma))
        pxite(1:kproma,jk) = MERGE(pxite(1:kproma,jk), ztmp1(1:kproma), ll1(1:kproma))

!>>DN: new diags
        zcorrxi2(1:kproma)  = pxite(1:kproma,jk)-ztmp2(1:kproma)
!<<DN: new diags        

        ll1(1:kproma) = (zxldt(1:kproma) > 0._dp)

!>>DN: new diags
        ztmp2(1:kproma) = pxlte(1:kproma,jk)
!<<DN: new diags        

        ztmp1(1:kproma)    = -(pxlm1(1:kproma,jk)/ztmst+zxlte2(1:kproma))
        !limit negative water tendency to avoid negative cloud water:
        ztmp1(1:kproma)    = MAX(pxlte(1:kproma,jk), ztmp1(1:kproma))
        pxlte(1:kproma,jk) = MERGE(pxlte(1:kproma,jk), ztmp1(1:kproma), ll1(1:kproma))

!>>DN: new diags
        zcorrxl2(1:kproma)  = pxlte(1:kproma,jk)-ztmp2(1:kproma)
!<<DN: new diags        

!>>SF comment: evap / sublm section
        ztmp1(1:kproma) = pxim1(1:kproma,jk)/MAX(paclc(1:kproma,jk),clc_min) !DN / SF #473: (lower clc limit)
        ztmp2(1:kproma) = pxlm1(1:kproma,jk)/MAX(paclc(1:kproma,jk),clc_min) !DN / SF #473: (lower clc limit)
    
        zxib(1:kproma,jk) = MERGE(ztmp1(1:kproma), 0._dp, ll_cc(1:kproma))
        zxlb(1:kproma,jk) = MERGE(ztmp2(1:kproma), 0._dp, ll_cc(1:kproma))

        zxim1evp(1:kproma) = MERGE(0._dp, pxim1(1:kproma,jk), ll_cc(1:kproma))
        zxlm1evp(1:kproma) = MERGE(0._dp, pxlm1(1:kproma,jk), ll_cc(1:kproma)) 

        !-- ice cloud:

        ll1(1:kproma) = (zxidt(1:kproma) > 0._dp)

        zxidtstar(1:kproma) = MERGE(zxidt(1:kproma), 0._dp, ll1(1:kproma))
       
        ztmp1(1:kproma) = zxidt(1:kproma)/MAX(paclc(1:kproma,jk), clc_min) !DN / SF #473: (lower clc limit)
        ztmp1(1:kproma) = MAX(ztmp1(1:kproma), -zxib(1:kproma,jk))
        ztmp1(1:kproma) = MERGE(zxidt(1:kproma), ztmp1(1:kproma), ll1(1:kproma))
        ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll_cc(1:kproma))

        zxib(1:kproma,jk) = zxib(1:kproma,jk) + ztmp1(1:kproma)

        !--- evaporation of the negative ice tendency
        ll2(1:kproma) = (.NOT. ll_cc(1:kproma)) .AND. (.NOT. ll1(1:kproma))

        ztmp1(1:kproma) = ztmst * ( pxite(1:kproma,jk) + zxite2(1:kproma) )
        ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma))

        zxim1evp(1:kproma) = zxim1evp(1:kproma) + ztmp1(1:kproma)

        !-- water cloud:
        ll1(1:kproma) = (zxldt(1:kproma) > 0._dp)
        zxldtstar(1:kproma) = MERGE(zxldt(1:kproma), 0._dp, ll1(1:kproma))

        ztmp1(1:kproma) = zxldt(1:kproma)/MAX(paclc(1:kproma,jk), clc_min) !SF #460: add a missing conversion
                                                                           !         to in-cloud
                                                                           !+ DN / SF #473: (lower clc limit)
        !limit negative water tendency to avoid negative zxlb:
        !>>DN: follow-up of #460
        ztmp1(1:kproma) = MAX(ztmp1(1:kproma), -zxlb(1:kproma,jk))
        ztmp1(1:kproma) = MERGE(zxldt(1:kproma), ztmp1(1:kproma), ll1(1:kproma))
        !<<DN: follow-up of #460
        ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll_cc(1:kproma))

        zxlb(1:kproma,jk) = zxlb(1:kproma,jk) + ztmp1(1:kproma)

        !--- evaporation of the negative water tendency
        ll2(1:kproma) = (.NOT. ll_cc(1:kproma)) .AND. (.NOT. ll1(1:kproma))

        ztmp1(1:kproma)    = ztmst * ( pxlte(1:kproma,jk) + zxlte2(1:kproma) )
        ztmp1(1:kproma)    = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma))

        zxlm1evp(1:kproma) = zxlm1evp(1:kproma) + ztmp1(1:kproma)

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_4incloudwi)
           CALL timer_start(timer_cmp_5conddepevapsubl)
        END IF
        !<<UP

!<<SF comment evap / sublm section
!
!     ------------------------------------------------------------------
!       5.    Condensation/deposition and evaporation/sublimation
!
!             zlc       =  L_{v/s} / c_p
!
        zlc(1:kproma)      = MERGE(zlsdcp(1:kproma,jk), zlvdcp(1:kproma,jk), lo2(1:kproma))
        zqsm1(1:kproma,jk) = MERGE(zqsi(1:kproma,jk), zqsw(1:kproma,jk),     lo2(1:kproma))
        zqst1(1:kproma)    = MERGE(zqsip1(1:kproma,jk), zqswp1(1:kproma,jk), lo2(1:kproma))

        zdqsdt(1:kproma) = 1000._dp*(zqst1(1:kproma)-zqsm1(1:kproma,jk))

        zxievap(1:kproma) = (1.0_dp-paclc(1:kproma,jk)) * zxidtstar(1:kproma) + zxim1evp(1:kproma)
!DN 1st term on right hand side is for (partial) evaporation
!   of cloud condensate from convection/transport/precipitation (positive ice/water tendencies are handled here)
        zxlevap(1:kproma) = (1.0_dp-paclc(1:kproma,jk)) * zxldtstar(1:kproma) + zxlm1evp(1:kproma)

        zqvdt(1:kproma) = ztmst*pqte(1:kproma,jk) !RD/DN #556

        zdtdt(1:kproma) = ztmst*ptte(1:kproma,jk) &
                        -  zlvdcp(1:kproma,jk)*(zevp(1:kproma) + zxlevap(1:kproma)) &
                        - (zlsdcp(1:kproma,jk)-zlvdcp(1:kproma,jk)) &
                        * (zsmlt(1:kproma) + zximlt(1:kproma)+zimlt(1:kproma))

        zdtdt(1:kproma) = zdtdt(1:kproma) &
                        - zlsdcp(1:kproma,jk)*(zsub(1:kproma)+zxievap(1:kproma)+zxisub(1:kproma))

        zqp1(1:kproma) = pqm1(1:kproma,jk)+zqvdt(1:kproma)
        zqp1(1:kproma) = MAX(zqp1(1:kproma),0.0_dp)

        ztp1(1:kproma) = ptm1(1:kproma,jk)+zdtdt(1:kproma)

        zxib(1:kproma,jk) = MAX(zxib(1:kproma,jk),0.0_dp)
        zxlb(1:kproma,jk) = MAX(zxlb(1:kproma,jk),0.0_dp)
        zxilb(1:kproma)   = zxib(1:kproma,jk) + zxlb(1:kproma,jk)

        zdqsat(1:kproma) = zdtdt(1:kproma) &
                         + paclc(1:kproma,jk)*( ztmst*zlc(1:kproma)*pqte(1:kproma,jk) &
                         + zlvdcp(1:kproma,jk)*(zevp(1:kproma)+zxlevap(1:kproma)) &
                         + zlsdcp(1:kproma,jk)*(zsub(1:kproma)+zxievap(1:kproma)+zxisub(1:kproma)) )   !zdtdtstar 

        zdqsat(1:kproma) = zdqsat(1:kproma) * zdqsdt(1:kproma) &
                         / (1._dp+paclc(1:kproma,jk)*zlc(1:kproma)*zdqsdt(1:kproma))

        zqcdif(1:kproma) = (zqvdt(1:kproma)-zdqsat(1:kproma))*paclc(1:kproma,jk)

        ztmp1(1:kproma) = -zxilb(1:kproma)*paclc(1:kproma,jk)
        ztmp2(1:kproma) = qsec*zqp1(1:kproma)

        zqcdif(1:kproma) = MAX(zqcdif(1:kproma), ztmp1(1:kproma))
        zqcdif(1:kproma) = MIN(zqcdif(1:kproma), ztmp2(1:kproma))  ! limit to qv

        ll1(1:kproma) = (zqcdif(1:kproma) < 0._dp)     !SF cloud dissipation
        
        ztmp1(1:kproma) = MAX(epsec, zxilb(1:kproma))
        ztmp1(1:kproma) = zxib(1:kproma,jk) / ztmp1(1:kproma)
        ztmp1(1:kproma) = MAX(MIN(ztmp1(1:kproma), 1.0_dp), 0.0_dp) !zifrac

        ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 1._dp, ll1(1:kproma))

        ztmp2(1:kproma) = zqcdif(1:kproma)*(1.0_dp-ztmp1(1:kproma))
        zcnd(1:kproma)  = MERGE(ztmp2(1:kproma), 0._dp, ll1(1:kproma))

        !>>UP emulator phasing nic_cirrus 1/2
        IF (lemuphase_nic_cirrus) THEN
                zdep_1(1:kproma) = zqcdif(1:kproma)*ztmp1(1:kproma)  ! zqcdif gets computed for both before
                zdep_2(1:kproma) = zqinucl(1:kproma,jk)*ztmp1(1:kproma) ! zqinucl gets computed only for nic_cirrus=2 before,
                                                                        ! otherwise 0
                zdep(1:kproma)   = eta_emu_nic_cirrus * zdep_1 + (1 - eta_emu_nic_cirrus) * zdep_2
        ELSE
        !<<UP
        IF (nic_cirrus == 1) THEN
           zdep(1:kproma) = zqcdif(1:kproma)*ztmp1(1:kproma)
        ELSE IF (nic_cirrus == 2) THEN 
           zdep(1:kproma) = zqinucl(1:kproma,jk)*ztmp1(1:kproma)
        ENDIF
        ENDIF !UP lemuphase_nic_cirrus

        ll2(1:kproma) = (.NOT. ll1(1:kproma)) .AND. &  !SF no cloud dissipation
                        (.NOT. lo2(1:kproma))          !SF condensation

        zdep(1:kproma) = MERGE(0._dp, zdep(1:kproma), ll2(1:kproma))

!--- Included/changed for prognostic CDNC/IC scheme --------------------
!    Use saturation adjustment for condensation of water vapor to cloud water
!    (required when using the Sundqvist cloud cover scheme, which is currently the only available).
!    See detailed reasons in #485.

        zcnd(1:kproma) = MERGE(zqcdif(1:kproma), zcnd(1:kproma), ll2(1:kproma))

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_5conddepevapsubl)
           CALL timer_start(timer_cmp_5clearairevap)
        END IF
        !<<UP

!--- End included for CDNC/IC scheme -----------------------------------
!
!
!       5.4 Accounting for cloud evaporation in clear air and
!           checking for supersaturation
!

!SF Note/ToDo: against the prevalent philosophy in the rest of the code, I don't implement a logical switch
!              to offer the possibility to bypass these calculations. With a switch, I'd need to code surrogate values
!              for the 3 INTENT(out) vars, ie zqp1tmp(1:kproma), zqsp1tmp(1:kproma), ztp1tmp(1:kproma,jk), 
!              and this implies quite some code (with some code duplication with respect to the content of this
!              subroutine). Maybe reconsider later?

        CALL mixed_phase_deposition_and_corrections( &
         
                !-- IN
                kbdim, kproma, &
                papp1(:,jk), picnc(:,jk), pqm1(:,jk), &
                paclc(:,jk), zesi(:,jk), zesw_2d(:,jk), &
                zeta(:,jk), zgenti(:), zlsdcp(:,jk), &
                zlvdcp(:,jk), zqp1(:), zqsm1(:,jk), &
                zrho(:,jk), ztp1(:), zxievap(:), &
                zxip1(:,jk), zxite(:), zvervx(:,jk), &
                !-- INOUT 
                zcnd(:), zdep(:), & 
                !-- OUT
                zqp1tmp(:), zqsp1tmp(:), ztp1tmp(:,jk))

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_5clearairevap)
           CALL timer_start(timer_cmp_5incloudconddepevapsubl)
        END IF
        !<<UP

!
!       5.5 Change of in-cloud water due to deposition/sublimation and
!           condensation/evaporation (input for cloud microphysics)
!

!>>DN: new diags
        ztmp3(1:kproma) = picnc(1:kproma,jk) 
!<<DN: new diags
        IF (lctrl_upd_incl_wat) &

!DN activation/nucleation is also updated in this subroutine

           CALL update_in_cloud_water( &

                   !-- IN
                   kbdim, kproma, &
                   zap(:,jk), zcdncact(:,jk), zcnd(:), & 
                   zdep(:), zgenti(:), zgentl(:), &
                   znicex(:,jk), zqp1tmp(:), zqsp1tmp(:), &
                   zrho(:,jk), zrid(:,jk), ptm1(:,jk), zcdnc(:,jk), &
                   picnc(:,jk), &
                   !-- INOUT
                   ll_cc(:), zicnc_upd(:,jk), qnuc(:,jk,krow), &
                   zcdnc_upd(:,jk), paclc(:,jk), zxib(:,jk), zxlb(:,jk), &
                   zninucl(:,jk), & ! UP #783.3
                   !-- OUT
                   zcdnc_min(:,jk), ztmp1(:), ztmp2(:), & !DN #475 (min cdnc) !davidn
                   zcorric3) ! UP #783.3
    
        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_5incloudconddepevapsubl)
           CALL timer_stop(timer_cmp_cmp)
           CALL timer_start(timer_cmp_diagnostics)
        END IF
        !<<UP

!>>UP, #783.3
        IF (lpr_corr) THEN
                ! recreate ll1 from update_in_cloud_water
                ll1(1:kproma) = (paclc(1:kproma,jk) > 0._dp) &
                        .AND. (zxlb(1:kproma,jk) > cqtmin)
           
        !>>DN: new diags
                qnucl_upd(1:kproma,jk) = qnucl_upd(1:kproma,jk) + &
                     zdt*ztmp1(1:kproma)
                !UP: add the cdnc_min correction from inside the routine to zcorrcd4
                zcorrcd4(1:kproma, jk) = zcorrcd4(1:kproma,jk) + &
                     MAX(zcdnc_min(1:kproma,jk)-(ztmp1(1:kproma)+ztmp2(1:kproma)),0._dp)
                !UP: zcorrcd7 is the correction term of evaporating where there is
                ! no cloud
                ! We do not want to count corrections that add something outside
                ! a cloud (e.g. cd4) and then again count cd7 removing that
                ! Instead, we correct the other correction terms (below), and remove
                ! their contribution from cd7
                ! ztmp1 does not need to go into here, because outside of clouds
                ! we are not counting it (not for qnucl_upd above and not for
                ! pcdnc_upd inside the routine)
                zcorrcd7(1:kproma) = MERGE(0._dp, &
                MIN(zcdnc(1:kproma,jk) &
                -(ztmp2(1:kproma)-zcorrcd4(1:kproma,jk)-zcorrcd1_1(1:kproma,jk)&
                +zcorrcd1_2(1:kproma,jk)+zcorrcd1_3(1:kproma,jk)-zcorrcd1_4(1:kproma,jk)&
                -zcorrcd2_1(1:kproma,jk)+zcorrcd2_2(1:kproma,jk)+zcorrcd2_3(1:kproma,jk)&
                -zcorrcd6(1:kproma)-qnucl_upd(1:kproma,jk)/zdt), 0._dp),ll1)
                !UP: explanation of above: inside clouds we don't evaporate;
                ! outside clouds there has been evaporation where zdcnc_after -
                ! zcdnc_before < 0; correction terms are
                ! substracted or added depending on their sign and meaning (see above)
                ! e.g. zcorrcd1_1 raises to min cdnc -> needs to be substracted;
                ! zcorrcd1_2 cd2ic was substracted from cdnc before, now needs
                ! to be added
                ! As the nucleation is also corrected for wc=0 below, it also
                ! needs to be substracted here
                !UP: need to get the left-over change of
                ! setting to 0 where there is no cloud
                ! Inside clouds we don't evaporate
                ! Outside clouds we might have added something with corrections before, which we
                ! will remove from the corrections below -> now we also need to remove it
                ! from this correction term; else, it would only be listed as an
                ! evaporated sink (in CD7), but not as a source of CD (in other
                ! corrections)
                !
                !UP: correct nucleation rates for not happening where there is no cloud
                !water
                ! This does it in one go also for the nucleation that happens before
                ! update_in_cloud_water
                ! Need to correct the temporary update terms,
                ! and then they can be added to the real process rates
                qnucl_upd(1:kproma,jk) = MERGE(qnucl_upd(1:kproma,jk), 0._dp, ll1)
                ! The update terms contain only the change in this routine step
                ! (since they are initialised with 0 at the beginning of the
                ! routine),
                ! which is why they get added and do not replace dnucl and qnucl 
                qnucl(1:kproma,jk,krow) = qnucl(1:kproma,jk,krow) + qnucl_upd(1:kproma,jk)
                dnucl(1:kproma,krow)    = dnucl(1:kproma,krow) + qnucl_upd(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
                !UP: correct corrections for not happening where there is no cloud water
                !UP: setting to cqtmin or 0._dp for the below ll1 cases does not
                ! matter because both are so small
                !UP: correct correction rates for cw=0 here
                ! The thinking is that qcdnc and qnucl are getting corrected for
                ! cw=0, so also the correction terms that have been acting on
                ! them need to see correction for cw=0
                zcorrcd1_1(1:kproma,jk) = MERGE(zcorrcd1_1(1:kproma,jk),0._dp,ll1)
                zcorrcd1_2(1:kproma,jk) = MERGE(zcorrcd1_2(1:kproma,jk),0._dp,ll1)
                zcorrcd1_3(1:kproma,jk) = MERGE(zcorrcd1_3(1:kproma,jk),0._dp,ll1)
                zcorrcd1_4(1:kproma,jk) = MERGE(zcorrcd1_4(1:kproma,jk),0._dp,ll1)
                zcorrcd2_1(1:kproma,jk) = MERGE(zcorrcd2_1(1:kproma,jk),0._dp,ll1)
                zcorrcd2_2(1:kproma,jk) = MERGE(zcorrcd2_2(1:kproma,jk),0._dp,ll1)
                zcorrcd2_3(1:kproma,jk) = MERGE(zcorrcd2_3(1:kproma,jk),0._dp,ll1)
                zcorrcd4(1:kproma,jk) = MERGE(zcorrcd4(1:kproma,jk), 0._dp, ll1)
                zcorrcd6(1:kproma) = MERGE(zcorrcd6(1:kproma), 0._dp, ll1)

                ! ICNC corrections need a new logical (as in
                ! update_in_cloud_water):
                ll1(1:kproma) = (paclc(1:kproma,jk) > 0._dp) &
                        .AND. (zxib(1:kproma,jk) > cqtmin)

                zcorric5(1:kproma) = MERGE(0._dp, MIN(picnc(1:kproma,jk)&
                - (ztmp3(1:kproma) - zcorric1_1(1:kproma,jk) - zcorric1_2(1:kproma,jk) &
                - zcorric2(1:kproma,jk) -zcorric4(1:kproma) + zcorric1_3(1:kproma,jk) &
                - zcorric7(1:kproma)),&
                 0._dp), ll1)
                ! UP: correction terms from inside the routine don't need to be
                ! considered here because they are already only applied where
                ! ll1 is true
                !UP: correct IC correction terms
                zcorric1_1(1:kproma,jk) = MERGE(zcorric1_1(1:kproma,jk),0._dp,ll1)
                zcorric1_2(1:kproma,jk) = MERGE(zcorric1_2(1:kproma,jk),0._dp,ll1)
                zcorric1_3(1:kproma,jk) = MERGE(zcorric1_3(1:kproma,jk),0._dp,ll1)
                zcorric2(1:kproma,jk) = MERGE(zcorric2(1:kproma,jk),0._dp,ll1)
                zcorric4(1:kproma) = MERGE(zcorric4(1:kproma),0._dp,ll1)
                zcorric7(1:kproma) = MERGE(zcorric7(1:kproma),0._dp,ll1)
                !<<DN: new diags
!<<UP, #783.3
        ELSE
        !>>DN: new diags
                qnucl(1:kproma,jk,krow) = qnucl(1:kproma,jk,krow) + &
                     zdt*ztmp1(1:kproma)
                dnucl(1:kproma,krow)    = dnucl(1:kproma,krow) + &
                     zdt*ztmp1(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
                zcorrcd7(1:kproma) = zcdnc(1:kproma,jk)-ztmp2(1:kproma)
                zcorric5(1:kproma) = (picnc(1:kproma,jk)-ztmp3(1:kproma))
        !<<DN: new diags
        END IF


!
!     ------------------------------------------------------------------
!       6.    Freezing of cloud water
!
!       6.1   Freezing of all cloud water for T < 238 K
!
        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_diagnostics)
           CALL timer_start(timer_cmp_cmp)
           CALL timer_start(timer_cmp_6freezingTs238)
        END IF
        !<<UP

        !UP: original line; uncommented below to remain on linear branch, TODO: take
        !back in when sure of this change
        ! As detailed in r6806 this change requires the model to be retuned, 
        ! which is why we leave it for later
        ll_frz_below_238K(1:kproma) = (ztp1tmp(1:kproma,jk) <= cthomi)
        !>>UP #783.3
        ! freezing should only happen where there is liquid water mass present
!        IF (lpr_corr) THEN
!                ll_frz_below_238K(1:kproma) = (ztp1tmp(1:kproma,jk) <= cthomi) &
!                                      .AND. (zxlb(1:kproma,jk) > cqtmin)
!        ELSE
!                ! original line
!                ll_frz_below_238K(1:kproma) = (ztp1tmp(1:kproma,jk) <= cthomi)
!        ENDIF
        !<<UP

        !>>UP write out ICNC befor the freezing routine
        IF (lcmptend) THEN
                icnc_bfrs238(1:kproma,jk,krow) = icnc_bfrs238(1:kproma,jk,krow) + picnc(1:kproma,jk)*zdtime
                !tend_testdummy_1(1:kproma,jk,krow) = zicnc_upd(1:kproma,jk) ! watch out, this should also accumulate
        ENDIF
        !<<UP
!davidn
        ztmp1(1:kproma)=0._dp
!davidn

        !>>UP had to change the structure of this from a one-line if-statement to multiple-line
        IF (lctrl_frz_below_238K .AND. ANY(ll_frz_below_238K(1:kproma))) THEN

           !>>UP increase ctr for this freezing routine
           icnc_ctrfrs238(1:kproma,jk,krow) = icnc_ctrfrs238(1:kproma,jk,krow) + 1.0_dp*zdtime
           !<<UP

           CALL freezing_below_238K( &

                   !-- IN
                   kbdim, kproma, &
                   ll_frz_below_238K(:), paclc(:,jk), zcdnc_min(:,jk), & !DN #475 (min cdnc)
                   zcdnc(:,jk), picnc(:,jk), &
                   !-- INOUT
                   zicnc_upd(:,jk), qfre(:,jk,krow), zcdnc_upd(:,jk), &
                   zfrl(:,jk), zxib(:,jk), zxlb(:,jk), &
                   !-- OUT
                   zfrzn(:), ztmp1(:) )!davidn
        !>>UP #783
        zfrznhom(1:kproma)    = zfrznhom(1:kproma) + zfrzn(1:kproma)
        !<<UP #783
!>>DN: new diags
        zcorrcd10(1:kproma) = MERGE(-ztmp1(1:kproma),0._dp , ll_frz_below_238K(1:kproma))
!<<DN: new diags        
        !>>UPcomment: above is the thing with the prognostic/diagnostic switch: you input picnc, which is ICNC always, but output
        !zicnc_upd, which containts ICNC after this freezing. In the prognostic case, it points to picnc, so that gets updated. In
        !the diagnostic case, it points to a blackhole, so it does not get used.
        !<<UPcomment
        ENDIF
        !<<UP end of structrual change

        !>>UP write out ICNC after the freezing routine
        IF (lcmptend) THEN
                icnc_pfrs238(1:kproma,jk,krow) =  &
                        icnc_pfrs238(1:kproma,jk,krow) + zicnc_upd(1:kproma,jk)*zdtime
                icnc_dfrs238(1:kproma,jk,krow) = &
                        !icnc_dfrs238(1:kproma,jk,krow) + & don't add this because this would be double counting differences
                        (icnc_pfrs238(1:kproma,jk,krow) - icnc_bfrs238(1:kproma,jk,krow))
                ! here I don't need to multiply with zdtime, because
                ! that is already done for the values themselves
                !tend_testdummy_2(1:kproma,jk,krow) = picnc(1:kproma,jk) ! watch out, this should also accumulate
        ENDIF
        !<<UP

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_6freezingTs238)
           CALL timer_start(timer_cmp_6freezingTl238)
        END IF
        !<<UP

!--- End included for CDNC/IC scheme -----------------------------------
!
!       6.2   Freezing of cloud water between 238 and 273 K
!
        ll_mxphase_frz(1:kproma) = (zxlb(1:kproma,jk) > cqtmin)     &
                             .AND. (ztp1tmp(1:kproma,jk) < tmelt )  &
                             .AND. (ztp1tmp(1:kproma,jk) > cthomi ) & 
                             .AND. (zcdnc(1:kproma,jk) >= zcdnc_min(1:kproma,jk))   & !DN #475 (min cdnc)
                             .AND. ll_cc(1:kproma)

        !>>UP write out ICNC before the freezing routine
        IF (lcmptend) THEN
                icnc_bfrl238(1:kproma,jk,krow) = icnc_bfrl238(1:kproma,jk,krow) + picnc(1:kproma,jk) * zdtime
                !tend_testdummy_1(1:kproma,jk,krow) = picnc(1:kproma,jk) ! watch out, this should also accumulate
        ENDIF
        !<<UP

!>>DN: burden
        dcldtte(1:kproma,krow) = dcldtte(1:kproma,krow)-&
             ptte(1:kproma,jk)*(cpd+zcons1*MAX(pqm1(1:kproma,jk),0.0_dp))*zdtime*zdpg(1:kproma,jk)
!<<DN: burden

        IF (ANY(ll_mxphase_frz(1:kproma))) THEN

           !>>UP write out ICNC befor the freezing routine
           !IF (lcmptend) THEN
                   !icnc_bfrl238(1:kproma,jk,krow) = picnc(1:kproma,jk)
           !        tend_testdummy_1(1:kproma,jk,krow) = tend_testdummy_1(1:kproma,jk,krow) + zicnc_upd(1:kproma,jk) * zdtime
           !ENDIF
           !<<UP

           !-- Heterogeneous mixed-phase freezing
           !MA start here
           IF (lctrl_het_mxphase_frz) THEN

              !>>UP increase ctr for this freezing routine
              icnc_ctrfrl238(1:kproma,jk,krow) = icnc_ctrfrl238(1:kproma,jk,krow) + 1.0_dp*zdtime
              !<<UP

              CALL het_mxphase_freezing( &
   
                      !-- IN
                      kbdim, kproma, &
                      ll_mxphase_frz(:), papp1(:,jk), ptkem1(:,jk), &
                      pvervel(:,jk), paclc(:,jk), zfracbcsol(:,jk), &
                      zfracbcinsol(:,jk), zfracdusol(:,jk), zfracduai(:,jk), &
                      zfracduci(:,jk),  &
                      zrho(:,jk), zrho_rcp(:,jk), zrwetki(:,jk), &
                      zrwetai(:,jk), zrwetci(:,jk), ztp1tmp(:,jk), &
                      zcdnc_min(:,jk), & !DN #475 (min cdnc)
                      zcdnc(:,jk), picnc(:,jk), &
                      !-- INOUT
                      zicnc_upd(:,jk), zcdnc_upd(:,jk), zfrl(:,jk), &
                      zxib(:,jk), zxlb(:,jk), &
                      !-- OUT
                      zfrln(:,jk))
           ENDIF !lctrl_het_mxphase_frz
   
           !>>UP write out ICNC after the freezing routine
        ENDIF !ANY(ll_mxphase_frz(1:kproma)) 
        !UPcomment: need to end and restart the if statement here so that 
        ! the update of icnc_pfrl238 happens every time

           IF (lcmptend) THEN
                   icnc_pfrl238(1:kproma,jk,krow) = icnc_pfrl238(1:kproma,jk,krow) + zicnc_upd(1:kproma,jk) * zdtime
                   icnc_dfrl238(1:kproma,jk,krow) =  & ! icnc_dfrl238(1:kproma,jk,krow) + 
                           !& don't add this because this would be double counting differences
                           (icnc_pfrl238(1:kproma,jk,krow) - icnc_bfrl238(1:kproma,jk,krow))
           ENDIF

        IF (ANY(ll_mxphase_frz(1:kproma))) THEN !UPcomment: restart if statement
           !<<UP

           !>>UP #783
           zfrznhet(1:kproma)    = zfrznhet(1:kproma) + zfrln(1:kproma,jk)
           !<<UP #783
!>>SF Kasja diags
           zfrzn(1:kproma)    = zfrzn(1:kproma) + zfrln(1:kproma,jk)
!<<SF Kasja diags

           !-- Wegener-Bergeron-Findeisen process
           IF (lctrl_WBF) THEN
  
              !>>SF #427: need to recompute zvervmax:

              !SF conversion of in-cloud ice mmr from kg/kg to in-cloud g m-3:
              ztmp1(1:kproma) = 1000._dp*zxib(1:kproma,jk)*zrho(1:kproma,jk)

              ztmp2(:) = eff_ice_crystal_radius(kbdim, kproma, ztmp1(:), picnc(:,jk))
              ztmp2(1:kproma) = MIN(MAX(ztmp2(1:kproma), ceffmin), ceffmax) !SF zrieff in micrometers

              ! Simple param of r/re approximating the Schumann et al. 2011 data:
              zrice(:,jk) = effective_2_volmean_radius_param_Schuman_2011( &
                                     kbdim, kproma, ztmp2(:))

              zvervmax(:,jk) = threshold_vert_vel(kbdim, kproma, zesw_2d(:,jk), &
                                         zesi(:,jk), picnc(:,jk), &
                                         zrice(:,jk), zeta(:,jk))
              !<<SF #427

              ll_WBF(1:kproma) = ll_mxphase_frz(1:kproma) &
                           .AND. ll_cc(1:kproma)          &
                           .AND. (zdep(1:kproma) > 0._dp) &
                           .AND. (zxlb(1:kproma,jk) > 0._dp) &
                           .AND. (0.01_dp*zvervx(1:kproma,jk) < zvervmax(1:kproma,jk))      
    
              IF (ANY(ll_WBF(1:kproma))) THEN

                 CALL WBF_process( &
   
                         !-- IN
                         kbdim, kproma, &
                         ll_WBF(:), paclc(:,jk), zlsdcp(:,jk), zlvdcp(:,jk), &
                         zcdnc(:,jk), &
                         !-- INOUT
                         zcdnc_upd(:,jk), zxlb(:,jk), zxib(:,jk), &
                         pxlte(:,jk), pxite(:,jk), ptte(:,jk), zevabf(:), zlwcbf(:))!davidn

              ENDIF !ANY(ll_WBF(1:kproma))

           ENDIF !lctrl_WBF

        ENDIF !ANY(ll_mxphase_frz(1:kproma))
    
        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_6freezingTl238)
           CALL timer_start(timer_cmp_7cloudphyscandsurfprecip)
        END IF
        !<<UP

!     ------------------------------------------------------------------
!       7.  Cloud physics and precipitation fluxes at the surface
!

        zclcstar(1:kproma) = MIN(paclc(1:kproma,jk), zclcpre(1:kproma))
        zauloc(1:kproma)   = 3./5000._dp*zdz(1:kproma,jk)
        zauloc(1:kproma)   = MAX(MIN(zauloc(1:kproma), clmax), clmin)

        ll1(1:kproma) = (knvb(1:kproma) >= jbmin) .AND. &
                        (knvb(1:kproma) <= jbmax) .AND. &
                        (pvervel(1:kproma,jk) > 0._dp)

        ll2(1:kproma) = (jk == knvb(1:kproma)  ) .OR. &
                        (jk == knvb(1:kproma)+1)

        ll3(1:kproma) = ll1(1:kproma) .AND. ll2(1:kproma) .AND. lonacc

        zauloc(1:kproma) = MERGE(0._dp, zauloc(1:kproma), ll3(1:kproma))

        zxlb(1:kproma,jk) = MAX(zxlb(1:kproma,jk),1.e-20_dp)
        zxib(1:kproma,jk) = MAX(zxib(1:kproma,jk), 1.e-20_dp)

!--- for in-cloud scavenging
        zmlwc(1:kproma,jk) = zxlb(1:kproma,jk)
        zmiwc(1:kproma,jk) = zxib(1:kproma,jk)

!---  Calculate the rain and snow water content in kg/kg from the rain and snow flux
       ll1(1:kproma) = (zclcpre(1:kproma) > eps )
       ll2(1:kproma) = ll1(1:kproma) .AND. (prsfl(1:kproma) > cqtmin)
       ll3(1:kproma) = ll1(1:kproma) .AND. (pssfl(1:kproma) > cqtmin)

       ztmp1(1:kproma) = ( MAX(prsfl(1:kproma), cqtmin)/(12.45_dp*MAX(zclcpre(1:kproma),eps) &
                         * SQRT(zqrho(1:kproma,jk))) )**(8._dp/9._dp) 
                       !SF eq. 10.70 of Roeckner et al 2003 (ie MPI report 349),
                       ! 12.45 comes from:
                       ! a_10*(n_0r)**(-1/8) = 90.8*(8.e6)**(-1/8)

       ztmp2(1:kproma) = ( MAX(pssfl(1:kproma), cqtmin)/(cvtfall *MAX(zclcpre(1:kproma),eps)) )**(1._dp/1.16_dp) 
                       !SF eq. 10.74 of Roeckner et al 2003 (ie MPI report 349)

       zxrp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma))
       zxsp1(1:kproma) = MERGE(ztmp2(1:kproma), 0._dp, ll3(1:kproma))

!
!       7.1   Warm clouds: Coalescence processes after Beheng (1994):
!             Autoconversion of cloud droplets and collection of cloud
!             droplets by falling rain. Accretion of cloud droplets by
!             falling snow (zsacl) is calculated under 7.2 !UP: i.e. riming
!
        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_start(timer_cmp_7warmclouds)
        END IF
        !<<UP

        ll_prcp_warm(1:kproma) = ll_cc(1:kproma) &
                           .AND. (zxlb(1:kproma,jk) > cqtmin) &
                           .AND. (zcdnc(1:kproma,jk) >= zcdnc_min(1:kproma,jk)) !DN #475 (min cdnc)

        IF (lctrl_prcp_warm .AND. ANY(ll_prcp_warm(1:kproma))) THEN

           CALL precip_formation_warm( &

                   !-- IN
                   kbdim, kproma, &
                   ll_prcp_warm(:), zauloc(:), paclc(:,jk), &
                   zclcstar(:), zrho(:,jk), zrho_rcp(:,jk), &
                   zxrp1(:), zcdnc_min(:,jk), & !DN #475 (min cdnc) 
                   zcdnc(:,jk), &
                   !-- INOUT
                   zcdnc_upd(:,jk), zxlb(:,jk), &
                   !-- OUT
                   zmratepr(:,jk), zrpr(:), zrprn(:,jk) ,zracl(:), zracln(:), zraut(:), zrautn(:), zcorrcd9(:))!davidn
        ELSE
            zmratepr(1:kproma,jk) = 0._dp
            zrprn(1:kproma,jk)    = 0._dp
            zrpr(1:kproma)        = 0._dp
            !davidn
            zracl(1:kproma)=0._dp
            zracln(1:kproma)=0._dp
            zraut(1:kproma)=0._dp
            zrautn(1:kproma)=0._dp
            zcorrcd9(1:kproma)=0._dp
            !davidn
        ENDIF

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_7warmclouds)
           CALL timer_start(timer_cmp_7coldclouds)
           CALL timer_start(timer_cmp_7ccpfcoutside) !UP #777
        END IF
        !<<UP

!       7.2  Cold clouds:
!            Conversion of cloud ice to snow after Levkov et al. 1992:
!            Aggregation of ice crystals to snow assuming plates (zsaut) and accretion of ice
!            by falling snow. (zsaci)
!            Accrection of cloud droplets by falling snow. (zsacl) !UP: i.e.
!            riming
!            Effective radius of ice crystals assuming plates (Lohmann, ACPD, 2007)

        IF (lctrl_prcp_cold) &

           CALL precip_formation_cold( &

                   !-- IN
                   kbdim, kproma, &
                   ll_cc(:), zauloc(:), paclc(:,jk), &
                   zclcstar(:), zqrho(:,jk), zrho_rcp(:,jk), &
                   ztp1tmp(:,jk), zviscos(:,jk), zxsp1(:), &
                   zrho(:,jk), zcdnc_min(:,jk), & !DN #475 (min cdnc)
                   zcdnc(:,jk), picnc(:,jk), &
                   jk, krow, & !UP #821
                   zcmpsimpl_rime(:,jk), zcmpsimpl_icaccr(:,jk), zcmpsimpl_sci(:,jk), & !UP #821
                   !-- INOUT
                   zicnc_upd(:,jk), zcdnc_upd(:,jk), zmrateps(:,jk), &
                   zxib(:,jk), zxlb(:,jk), &
                   !-- OUT
                   zsprn(:,jk), zsacl(:), zsacln(:,jk), zmsnowacl(:,jk), &
                   zspr(:), zsaut(:), zsaci(:), zaggn(:), zsacin(:), zselfn(:), zsecprodn(:), zsecprod(:) )!davidn

        !>>UP timers #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_7ccpfcoutside)
        END IF
        !<<UP

        !--- Update CDNC for radiation:
        IF (IAND(nclmi_progn,1) /= 0) THEN !prognostic CDNC (nclmi_progn==1 .or. 3)
            pacdnc(1:kproma,jk) = zcdnc(1:kproma,jk)
        ENDIF

        !>>UP CMP tendencies
        IF (lcmptend) THEN
           tend_sacl(1:kproma,jk,krow) = tend_sacl(1:kproma,jk,krow) + zsacl(1:kproma) * zdtime
        ENDIF
        !<<UP

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_7coldclouds)
           CALL timer_start(timer_cmp_7precipflux)
        END IF
        !<<UP

!       7.3 Updating precipitation fluxes. In the lowest layer (klev),
!           the sedimentation sink of cloud ice is balanced
!           by precipitation at the surface (through 'zzdrs').
!           Fraction of precipitating clouds (zclcpre) used for the
!           calculation of evaporation/sublimation of rain/snow in
!           the next layer

        IF (lctrl_upd_flx) THEN

           CALL update_precip_fluxes( &

                   !-- IN
                   kbdim, kproma, &
                   jk,klev, &
                   paclc(:,jk), zdp(:,jk), zevp(:), &
                   zlsdcp(:,jk), zlvdcp(:,jk), zrpr(:), &
                   zsacl(:), zspr(:), zsub(:), &
                   ztp1tmp(:,jk), zxiflux(:), &
                   !-- INOUT
                   zclcpre(:), prsfl(:), pssfl(:), zsmlt(:), &
                   !-- OUT
                   zfevapr(:,jk), zfrain(:,jk), zfsnow(:,jk), zfsubls(:,jk), ztmp2(:) ) !davidn
        ELSE
           zfevapr(1:kproma,jk) = 0._dp 
           zfrain(1:kproma,jk)  = 0._dp
           zfsnow(1:kproma,jk)  = 0._dp
           zfsubls(1:kproma,jk) = 0._dp
        ENDIF 
!>>DN: burden
        IF (jk .EQ. klev) THEN
           dsedi_sfc(1:kproma,krow)   = dsedi_sfc(1:kproma,krow)   + zdtime*zxiflux(1:kproma)
           dmlts_sfc(1:kproma,krow)   = dmlts_sfc(1:kproma,krow)   + zdtime*ztmp2(1:kproma)
        END IF
!<<DN: burden
                                 
!SF for scavenging: create a 2d array of clcpre:
        zclcpre_2d(1:kproma,jk) = zclcpre(1:kproma)

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_7precipflux)
           CALL timer_stop(timer_cmp_7cloudphyscandsurfprecip)
           CALL timer_start(timer_cmp_8updating)
        END IF
        !<<UP

!     ------------------------------------------------------------------
!       8.    Updating tendencies of t, q, xl, xi and final cloud cover
!
!

!>>SF auxiliary variables that are necessary for both 'update_tendencies_and_important_vars' and
!     'diagnostics'    
!     These var track the presence of liquid-, respectively ice clouds

        ll_liqcl(1:kproma,jk) = (zxlb(1:kproma,jk) > eps)      .AND. &
                                (zcdnc(1:kproma,jk) >= zcdnc_min(1:kproma,jk)) !DN #475 (min cdnc)

        ll_icecl(1:kproma,jk) = (zxib(1:kproma,jk) > eps)  .AND. &
                                (picnc(1:kproma,jk) >= icemin) !SF #463: this condition is new - was added to make this
                                                               !   flag consistent with its counterpart for liq phase.

!<<SF

!
!       8.3 + 8.4   Tendencies of thermodynamic variables (+ corrections to avoid negative liq water/ice)
!             Attn: The terms zxisub and zximlt do not appear in
!                   pxite because these processes have already been
!                   included in pxite via changes in cloud ice
!                   sedimentation (see 3.1, 3.2 and 4)
! 

!SFNote: in contrast to diagnostics, all the variable updated here *are* important for the following physics and 
!        dynamics

        IF (lctrl_upd_tend) &

           CALL update_tendencies_and_important_vars( &

                   !-- IN
                   kbdim, kproma, &
                   picnc(:,jk), zcdnc(:,jk), pxim1(:,jk), &
                   pxlm1(:,jk), pxtm1(:,jk,idt_cdnc), &
                   pxtm1(:,jk,idt_icnc), zcnd(:), &
                   zdep(:), zevp(:), zfrl(:,jk), &
                   zgenti(:), zgentl(:), zimlt(:), &
                   zlsdcp(:,jk), zlvdcp(:,jk), zrho(:,jk), zrho_rcp(:,jk), &
                   zrpr(:), zsacl(:), zspr(:), &
                   zxievap(:), zximlt(:), zxite(:), &
                   zxlevap(:), zxlte(:), zxisub(:), &
                   zsub(:), zsmlt(:), &
                   zxib(:,jk), zxlb(:,jk), ztp1tmp(:,jk), &
                   ll_liqcl(:,jk), ll_icecl(:,jk), &
                   !-- INOUT
                   paclc(:,jk), pqte(:,jk), &
                   ptte(:,jk), pxite(:,jk), pxlte(:,jk), &
                   pxtte(:,jk,idt_cdnc), pxtte(:,jk,idt_icnc), &
                   zmlwc(:,jk), zmiwc(:,jk), &
                   !-- OUT
                   reffl(:,jk,krow), reffi(:,jk,krow), zcorrcd8(:), zcorrxl(:), zcorric6(:), zcorrxi(:) )!davidn

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_cmp)
           CALL timer_stop(timer_cmp_8updating)
           CALL timer_start(timer_cmp_diagnostics)
        END IF
        !<<UP
!>>DN: burden
        dcldtte(1:kproma,krow) = dcldtte(1:kproma,krow)+&
             ptte(1:kproma,jk)*(cpd+zcons1*MAX(pqm1(1:kproma,jk),0.0_dp))*zdtime*zdpg(1:kproma,jk)
!<<DN: burden

!>>SF additional diags for Kasja

        qcnd(1:kproma,jk,krow)       = qcnd(1:kproma,jk,krow)      + zdt*zcnd(1:kproma)
!>>DN
        qlwc_detr(1:kproma,jk,krow)  = qlwc_detr(1:kproma,jk,krow) + zdtime*zxlte(1:kproma)
!<<DN
        qevp_lwc(1:kproma,jk,krow)   = qevp_lwc(1:kproma,jk,krow)  - zdt*zxlevap(1:kproma)
        qaut(1:kproma,jk,krow)       = qaut(1:kproma,jk,krow)      - zdt*zraut(1:kproma)
        qautn(1:kproma,jk,krow)      = qautn(1:kproma,jk,krow)     - zdt*zrautn(1:kproma)
        qracl(1:kproma,jk,krow)      = qracl(1:kproma,jk,krow)     - zdt*zracl(1:kproma)
        qracln(1:kproma,jk,krow)     = qracln(1:kproma,jk,krow)    - zdt*zracln(1:kproma)
        qsacl(1:kproma,jk,krow)      = qsacl(1:kproma,jk,krow)     - zdt*zsacl(1:kproma)
        qsacln(1:kproma,jk,krow)     = qsacln(1:kproma,jk,krow)    - zdt*zsacln(1:kproma,jk)
        qfrz(1:kproma,jk,krow)       = qfrz(1:kproma,jk,krow)      - zdt*zfrl(1:kproma,jk)
        qfrzn(1:kproma,jk,krow)      = qfrzn(1:kproma,jk,krow)     - zdt*zfrzn(1:kproma)
!>>UP #783
        qfrznhet(1:kproma,jk,krow)      = qfrznhet(1:kproma,jk,krow)     - zdt*zfrznhet(1:kproma)
        qfrznhom(1:kproma,jk,krow)      = qfrznhom(1:kproma,jk,krow)     - zdt*zfrznhom(1:kproma)
!<<UP #783

        qnuci(1:kproma,jk,krow)     = qnuci(1:kproma,jk,krow)     + zdt*zninucl(1:kproma,jk)
        qdep(1:kproma,jk,krow)      = qdep(1:kproma,jk,krow)      + zdt*zdep(1:kproma)
!>>DN
        qiwc_detr(1:kproma,jk,krow) = qiwc_detr(1:kproma,jk,krow) + zdtime*zxite(1:kproma)
!<<DN
        qsub_iwc(1:kproma,jk,krow)  = qsub_iwc(1:kproma,jk,krow)  - zdt*zxievap(1:kproma)
        qagg(1:kproma,jk,krow)      = qagg(1:kproma,jk,krow)      - zdt*zsaut(1:kproma)
        qaggn(1:kproma,jk,krow)     = qaggn(1:kproma,jk,krow)     - zdt*zaggn(1:kproma)
        qsaci(1:kproma,jk,krow)     = qsaci(1:kproma,jk,krow)     - zdt*zsaci(1:kproma)
        qsacin(1:kproma,jk,krow)    = qsacin(1:kproma,jk,krow)    - zdt*zsacin(1:kproma)
        qselfn(1:kproma,jk,krow)    = qselfn(1:kproma,jk,krow)    - zdt*zselfn(1:kproma)
!>>DN
        IF (lsecprod) THEN
           qsecprod(1:kproma,jk,krow)  = qsecprod(1:kproma,jk,krow)  + zdt*zsecprod(1:kproma)
!<<DN
           qsecprodn(1:kproma,jk,krow) = qsecprodn(1:kproma,jk,krow) + zdt*zsecprodn(1:kproma)
        END IF
        qsedi(1:kproma,jk,krow)     = qsedi(1:kproma,jk,krow)     - zdt*zsedi(1:kproma)
        qsedin(1:kproma,jk,krow)    = qsedin(1:kproma,jk,krow)    - zdt*zsedin(1:kproma)
        qmlt(1:kproma,jk,krow)      = qmlt(1:kproma,jk,krow)      - zdt*zimlt(1:kproma)
        qxmlt(1:kproma,jk,krow)     = qxmlt(1:kproma,jk,krow)     - zdt*zximlt(1:kproma)
        qmltn(1:kproma,jk,krow)     = qmltn(1:kproma,jk,krow)     - zdt*zmltn(1:kproma)

        qevp_rain(1:kproma,jk,krow)  = qevp_rain(1:kproma,jk,krow) - zdt*zevp(1:kproma)
        qsub_snow(1:kproma,jk,krow)  = qsub_snow(1:kproma,jk,krow) - zdt*zsub(1:kproma)
        qrprn(1:kproma,jk,krow)      = qrprn(1:kproma,jk,krow)     - zdt*zrprn(1:kproma,jk)
        !UP warning: self-collection wrongly included below
        qsprn(1:kproma,jk,krow)      = qsprn(1:kproma,jk,krow)     - zdt*zsprn(1:kproma,jk)
!<<SF additional diags for Kasja

!>>DN: new diags
        qevabf(1:kproma,jk,krow)      = qevabf(1:kproma,jk,krow)      - zdt*zlwcbf(1:kproma)
        qdepbf(1:kproma,jk,krow)      = qdepbf(1:kproma,jk,krow)      + zdt*zlwcbf(1:kproma)
        qevabfn(1:kproma,jk,krow)     = qevabfn(1:kproma,jk,krow)     - zdt*zevabf(1:kproma)
        qicnc_detr(1:kproma,jk,krow)  = qicnc_detr(1:kproma,jk,krow)  + zdt*znidetr(1:kproma,jk) 
        qcdnc_detr(1:kproma,jk,krow)  = qcdnc_detr(1:kproma,jk,krow)  + zdt*zncdetr(1:kproma,jk)
        qcorric1(1:kproma,jk,krow)    = qcorric1(1:kproma,jk,krow)    + zdt*zcorric1(1:kproma,jk)
        qcorric2(1:kproma,jk,krow)    = qcorric2(1:kproma,jk,krow)    + zdt*zcorric2(1:kproma,jk)
!>>UP: new diags for activation correction terms
! and split correction terms further
        IF (lnewdiags) THEN
            qcorrcd2unphys_2d(1:kproma,jk,krow) = &
qcorrcd2unphys_2d(1:kproma,jk,krow) + zdt*zcorrcd2unphys_2d(1:kproma,jk)
            qcorrcd1_1(1:kproma,jk,krow) = &
qcorrcd1_1(1:kproma,jk,krow) + zdt*zcorrcd1_1(1:kproma,jk)
            qcorrcd1_2(1:kproma,jk,krow) = &
qcorrcd1_2(1:kproma,jk,krow) + zdt*zcorrcd1_2(1:kproma,jk)
            qcorrcd1_3(1:kproma,jk,krow) = &
qcorrcd1_3(1:kproma,jk,krow) + zdt*zcorrcd1_3(1:kproma,jk)
            qcorrcd1_4(1:kproma,jk,krow) = &
qcorrcd1_4(1:kproma,jk,krow) + zdt*zcorrcd1_4(1:kproma,jk)
            qcorrcd2_1(1:kproma,jk,krow) = &
qcorrcd2_1(1:kproma,jk,krow) + zdt*zcorrcd2_1(1:kproma,jk)
            qcorrcd2_2(1:kproma,jk,krow) = &
qcorrcd2_2(1:kproma,jk,krow) + zdt*zcorrcd2_2(1:kproma,jk)
!UP #783.3
            qcorrcd2_3(1:kproma,jk,krow) = &
qcorrcd2_3(1:kproma,jk,krow) + zdt*zcorrcd2_3(1:kproma,jk)
            qcorric1_1(1:kproma,jk,krow) = &
qcorric1_1(1:kproma,jk,krow) + zdt*zcorric1_1(1:kproma,jk)
!UP #783.3
            qcorric1_2(1:kproma,jk,krow) = &
qcorric1_2(1:kproma,jk,krow) + zdt*zcorric1_2(1:kproma,jk)
            qcorric1_3(1:kproma,jk,krow) = &
qcorric1_3(1:kproma,jk,krow) + zdt*zcorric1_3(1:kproma,jk)
!UP #783.3
        ENDIF
!<<UP: new diags
! UP comment: below needs to be commented out, see issue #767
!        qcorric3(1:kproma,jk,krow)    = qcorric3(1:kproma,jk,krow)    + zdt*zcorric3(1:kproma,jk)
        qcorric4(1:kproma,jk,krow)    = qcorric4(1:kproma,jk,krow)    + zdt*zcorric4(1:kproma)
        qcorric5(1:kproma,jk,krow)    = qcorric5(1:kproma,jk,krow)    + zdt*zcorric5(1:kproma)
        qcorric6(1:kproma,jk,krow)    = qcorric6(1:kproma,jk,krow)    + zdt*zcorric6(1:kproma)
        qcorric7(1:kproma,jk,krow)    = qcorric7(1:kproma,jk,krow)    + zdt*zcorric7(1:kproma)
        qcorrcd1(1:kproma,jk,krow)    = qcorrcd1(1:kproma,jk,krow)    + zdt*zcorrcd1(1:kproma,jk)
        qcorrcd2(1:kproma,jk,krow)    = qcorrcd2(1:kproma,jk,krow)    + zdt*zcorrcd2(1:kproma,jk)
        qcorrcd3(1:kproma,jk,krow)    = qcorrcd3(1:kproma,jk,krow)    + zdt*zcorrcd3(1:kproma,jk)
        qcorrcd4(1:kproma,jk,krow)    = qcorrcd4(1:kproma,jk,krow)    + zdt*zcorrcd4(1:kproma,jk)
        qcorrcd5(1:kproma,jk,krow)    = qcorrcd5(1:kproma,jk,krow)    + zdt*zcorrcd5(1:kproma,jk)
        qcorrcd6(1:kproma,jk,krow)    = qcorrcd6(1:kproma,jk,krow)    + zdt*zcorrcd6(1:kproma)
        qcorrcd7(1:kproma,jk,krow)    = qcorrcd7(1:kproma,jk,krow)    + zdt*zcorrcd7(1:kproma)
        qcorrcd8(1:kproma,jk,krow)    = qcorrcd8(1:kproma,jk,krow)    + zdt*zcorrcd8(1:kproma)
        qcorrcd9(1:kproma,jk,krow)    = qcorrcd9(1:kproma,jk,krow)    + zdt*zcorrcd9(1:kproma)
        qcorrcd10(1:kproma,jk,krow)   = qcorrcd10(1:kproma,jk,krow)   + zdt*zcorrcd10(1:kproma)
        qcorrxl(1:kproma,jk,krow)     = qcorrxl(1:kproma,jk,krow)     + zdt*zcorrxl(1:kproma)
        qcorrxl2(1:kproma,jk,krow)    = qcorrxl2(1:kproma,jk,krow)    + zdtime*zcorrxl2(1:kproma)
        qcorrxi(1:kproma,jk,krow)     = qcorrxi(1:kproma,jk,krow)     + zdt*zcorrxi(1:kproma)
        qcorrxi2(1:kproma,jk,krow)    = qcorrxi2(1:kproma,jk,krow)    + zdtime*zcorrxi2(1:kproma)
        qxlte(1:kproma,jk,krow)       = qxlte(1:kproma,jk,krow)       + zdtime*zactxlte(1:kproma,jk)
        qxite(1:kproma,jk,krow)       = qxite(1:kproma,jk,krow)       + zdtime*zactxite(1:kproma,jk)
        qxttecdnc(1:kproma,jk,krow)   = qxttecdnc(1:kproma,jk,krow)   + zdtime*zxttecdnc(1:kproma,jk)
        qxtteicnc(1:kproma,jk,krow)   = qxtteicnc(1:kproma,jk,krow)   + zdtime*zxtteicnc(1:kproma,jk)
        qspr(1:kproma,jk,krow)        = qspr(1:kproma,jk,krow)        + zdt*zspr(1:kproma)
        qrpr(1:kproma,jk,krow)        = qrpr(1:kproma,jk,krow)        + zdt*zrpr(1:kproma)
        qmltn2(1:kproma,jk,krow)      = qmltn2(1:kproma,jk,krow)      - zdt*zmltn2(1:kproma)

        ztestcd(1:kproma,jk) = ztestcd(1:kproma,jk)-ztmst*pxtte(1:kproma,jk,idt_cdnc)
        ztestic(1:kproma,jk) = ztestic(1:kproma,jk)-ztmst*pxtte(1:kproma,jk,idt_icnc)
        ztestlwc(1:kproma,jk) = ztestlwc(1:kproma,jk)-ztmst*pxlte(1:kproma,jk)
        ztestiwc(1:kproma,jk) = ztestiwc(1:kproma,jk)-ztmst*pxite(1:kproma,jk)
        qtestCD(1:kproma,jk,krow)     = qtestCD(1:kproma,jk,krow)     + zdt*ztestcd(1:kproma,jk)
        qtestIC(1:kproma,jk,krow)     = qtestIC(1:kproma,jk,krow)     + zdt*ztestic(1:kproma,jk)
        qtestLWC(1:kproma,jk,krow)    = qtestLWC(1:kproma,jk,krow)    + zdt*ztestlwc(1:kproma,jk)
        qtestIWC(1:kproma,jk,krow)    = qtestIWC(1:kproma,jk,krow)    + zdt*ztestiwc(1:kproma,jk)
        qgentl(1:kproma,jk,krow)      = qgentl(1:kproma,jk,krow)      + zdt*zgentl(1:kproma)
        qgenti(1:kproma,jk,krow)      = qgenti(1:kproma,jk,krow)      + zdt*zgenti(1:kproma)
        qsub_ice(1:kproma,jk,krow)    = qsub_ice(1:kproma,jk,krow)    + zdt*zxisub(1:kproma)
        qmlt_snow(1:kproma,jk,krow)   = qmlt_snow(1:kproma,jk,krow)   + zdt*zsmlt(1:kproma)
!<<DN: new diags

!>>UP: new diags for lifetimes
        IF (lnewdiags) THEN
        z_lwc_sources(1:kproma,jk) = qcnd(1:kproma,jk,krow) + qlwc_detr(1:kproma,jk,krow)       &
                                        - qmlt(1:kproma,jk,krow) - qxmlt(1:kproma,jk,krow)           &
                                        - qcorrxl(1:kproma,jk,krow) + qcorrxl2(1:kproma,jk,krow)
        z_lwc_sinks(1:kproma,jk) = qevp_lwc(1:kproma,jk,krow) - qracl(1:kproma,jk,krow)         &
                                      - qaut(1:kproma,jk,krow) + qsacl(1:kproma,jk,krow)             &
                                      + qfrz(1:kproma,jk,krow) + qevabf(1:kproma,jk,krow)
        z_iwc_sources(1:kproma,jk) = qdep(1:kproma,jk,krow) - qfrz(1:kproma,jk,krow)            &
                                        + qiwc_detr(1:kproma,jk,krow) + qdepbf(1:kproma,jk,krow)     &
                                        + qcorrxi(1:kproma,jk,krow)
        z_iwc_sinks(1:kproma,jk) = - qagg(1:kproma,jk,krow) - qsaci(1:kproma,jk,krow)           &
                                        + qsedi(1:kproma,jk,krow) + qsub_iwc(1:kproma,jk,krow)       &
                                        + qsecprod(1:kproma,jk,krow) 
        z_icnc_sources(1:kproma,jk) = qnuci(1:kproma,jk,krow) + qcorric2(1:kproma,jk,krow)      &
                                         - qfrzn(1:kproma,jk,krow) + qcorric1_2(1:kproma,jk,krow)    &
                                         + qcorric2(1:kproma,jk,krow) + qicnc_detr(1:kproma,jk,krow) &
                                         + qsecprodn(1:kproma,jk,krow) + qcorric1_1(1:kproma,jk,krow)
        z_icnc_sinks(1:kproma,jk) = - qaggn(1:kproma,jk,krow) - qsacin(1:kproma,jk,krow)        &
                                         - qselfn(1:kproma,jk,krow) + qsedin(1:kproma,jk,krow)       &
                                         + qmltn2(1:kproma,jk,krow) - qcorric4(1:kproma,jk,krow)     &
                                         + qcorric1_3(1:kproma,jk,krow) - qcorric6(1:kproma,jk,krow) &
                                         - qcorric7(1:kproma,jk,krow) - qcorric5(1:kproma,jk,krow)
        z_cdnc_sources(1:kproma,jk) = qnucl(1:kproma,jk,krow) - qcorrcd1_3(1:kproma,jk,krow)    &
                                         - qcorrcd2_2(1:kproma,jk,krow) - qmltn(1:kproma,jk,krow)    &
                                         + qcorrcd6(1:kproma,jk,krow) + qcorrcd1_4(1:kproma,jk,krow) &
                                         + qcdnc_detr(1:kproma,jk,krow) + qcorrcd3(1:kproma,jk,krow) &
                                         + qcorrcd4(1:kproma,jk,krow) + qcorrcd9(1:kproma,jk,krow)   &
                                         + qcorrcd1_1(1:kproma,jk,krow) + qcorrcd2_1(1:kproma,jk,krow)
        z_cdnc_sinks(1:kproma,jk) = - qracln(1:kproma,jk,krow) - qautn(1:kproma,jk,krow)        &
                                         - qsacln(1:kproma,jk,krow) - qfre(1:kproma,jk,krow)         &
                                         + qcorrcd1_2(1:kproma,jk,krow)                              &
                                         + qcorrcd2_3(1:kproma,jk,krow) - qevabfn(1:kproma,jk,krow)  &
                                         - qcorrcd8(1:kproma,jk,krow) - qcorrcd7(1:kproma,jk,krow)
        q_lwc_sources(1:kproma,jk,krow) = q_lwc_sources(1:kproma,jk,krow) + zdt*z_lwc_sources(1:kproma,jk)
        q_lwc_sinks(1:kproma,jk,krow) = q_lwc_sinks(1:kproma,jk,krow) + zdt*z_lwc_sinks(1:kproma,jk)
        q_iwc_sources(1:kproma,jk,krow) = q_iwc_sources(1:kproma,jk,krow) + zdt*z_iwc_sources(1:kproma,jk)
        q_iwc_sinks(1:kproma,jk,krow) = q_iwc_sinks(1:kproma,jk,krow) + zdt*z_iwc_sinks(1:kproma,jk)
        q_icnc_sources(1:kproma,jk,krow) = q_icnc_sources(1:kproma,jk,krow) + zdt*z_icnc_sources(1:kproma,jk)
        q_icnc_sinks(1:kproma,jk,krow) = q_icnc_sinks(1:kproma,jk,krow) + zdt*z_icnc_sinks(1:kproma,jk)
        q_cdnc_sources(1:kproma,jk,krow) = q_cdnc_sources(1:kproma,jk,krow) + zdt*z_cdnc_sources(1:kproma,jk)
        q_cdnc_sinks(1:kproma,jk,krow) = q_cdnc_sinks(1:kproma,jk,krow) + zdt*z_cdnc_sinks(1:kproma,jk)
        ENDIF
!<<UP

!>>DN: burden
        dfre(1:kproma,krow)        = dfre(1:kproma,krow)        - zdt*zfrln(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcnd(1:kproma,krow)        = dcnd(1:kproma,krow)        + zdt*zcnd(1:kproma)*zdpg(1:kproma,jk)
        dlwc_detr(1:kproma,krow)   = dlwc_detr(1:kproma,krow)   + zdtime*zxlte(1:kproma)*zdpg(1:kproma,jk)
        devp_lwc(1:kproma,krow)    = devp_lwc(1:kproma,krow)    - zdt*zxlevap(1:kproma)*zdpg(1:kproma,jk)
        daut(1:kproma,krow)        = daut(1:kproma,krow)        - zdt*zraut(1:kproma)*zdpg(1:kproma,jk)
        dautn(1:kproma,krow)       = dautn(1:kproma,krow)       - zdt*zrautn(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dracl(1:kproma,krow)       = dracl(1:kproma,krow)       - zdt*zracl(1:kproma)*zdpg(1:kproma,jk)
        dracln(1:kproma,krow)      = dracln(1:kproma,krow)      - zdt*zracln(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dsacl(1:kproma,krow)       = dsacl(1:kproma,krow)       - zdt*zsacl(1:kproma)*zdpg(1:kproma,jk)
        dsacln(1:kproma,krow)      = dsacln(1:kproma,krow)      - zdt*zsacln(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dfrz(1:kproma,krow)        = dfrz(1:kproma,krow)        - zdt*zfrl(1:kproma,jk)*zdpg(1:kproma,jk)
        dfrzn(1:kproma,krow)       = dfrzn(1:kproma,krow)       - zdt*zfrzn(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
!>>UP #783
        dfrznhet(1:kproma,krow)    = dfrznhet(1:kproma,krow)    - zdt*zfrznhet(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dfrznhom(1:kproma,krow)    = dfrznhom(1:kproma,krow)    - zdt*zfrznhom(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
!<<UP #783
        dnuci(1:kproma,krow)       = dnuci(1:kproma,krow)       +&
             zdt*zninucl(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        ddep(1:kproma,krow)        = ddep(1:kproma,krow)        + zdt*zdep(1:kproma)*zdpg(1:kproma,jk)
        diwc_detr(1:kproma,krow)   = diwc_detr(1:kproma,krow)   + zdtime*zxite(1:kproma)*zdpg(1:kproma,jk)
        dsub_iwc(1:kproma,krow)    = dsub_iwc(1:kproma,krow)    - zdt*zxievap(1:kproma)*zdpg(1:kproma,jk)
        dagg(1:kproma,krow)        = dagg(1:kproma,krow)        - zdt*zsaut(1:kproma)*zdpg(1:kproma,jk)
        daggn(1:kproma,krow)       = daggn(1:kproma,krow)       - zdt*zaggn(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dsaci(1:kproma,krow)       = dsaci(1:kproma,krow)       - zdt*zsaci(1:kproma)*zdpg(1:kproma,jk)
        dsacin(1:kproma,krow)      = dsacin(1:kproma,krow)      - zdt*zsacin(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dselfn(1:kproma,krow)      = dselfn(1:kproma,krow)      - zdt*zselfn(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        IF (lsecprod) THEN
           dsecprod(1:kproma,krow)    = dsecprod(1:kproma,krow)    + zdt*zsecprod(1:kproma)*zdpg(1:kproma,jk)
           dsecprodn(1:kproma,krow)   = dsecprodn(1:kproma,krow)   + &
                zdt*zsecprodn(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        END IF
        dsedi(1:kproma,krow)       = dsedi(1:kproma,krow)       - zdt*zsedi(1:kproma)*zdpg(1:kproma,jk)
        dsedin(1:kproma,krow)      = dsedin(1:kproma,krow)      - zdt*zsedin(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dmlt(1:kproma,krow)        = dmlt(1:kproma,krow)        - zdt*zimlt(1:kproma)*zdpg(1:kproma,jk)
        dxmlt(1:kproma,krow)       = dxmlt(1:kproma,krow)       - zdt*zximlt(1:kproma)*zdpg(1:kproma,jk)
        dmltn(1:kproma,krow)       = dmltn(1:kproma,krow)       - zdt*zmltn(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        devp_rain(1:kproma,krow)   = devp_rain(1:kproma,krow)   - zdt*zevp(1:kproma)*zdpg(1:kproma,jk)
        dsub_snow(1:kproma,krow)   = dsub_snow(1:kproma,krow)   - zdt*zsub(1:kproma)*zdpg(1:kproma,jk)
        drprn(1:kproma,krow)       = drprn(1:kproma,krow)       - zdt*zrprn(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        !UP warning: self-collection wrongly included below
        dsprn(1:kproma,krow)       = dsprn(1:kproma,krow)       - zdt*zsprn(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        devabf(1:kproma,krow)      = devabf(1:kproma,krow)      - zdt*zlwcbf(1:kproma)*zdpg(1:kproma,jk)
        ddepbf(1:kproma,krow)      = ddepbf(1:kproma,krow)      + zdt*zlwcbf(1:kproma)*zdpg(1:kproma,jk)
        devabfn(1:kproma,krow)     = devabfn(1:kproma,krow)     - zdt*zevabf(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dicnc_detr(1:kproma,krow)  = dicnc_detr(1:kproma,krow)  +&
             zdt*znidetr(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk) 
        dcdnc_detr(1:kproma,krow)  = dcdnc_detr(1:kproma,krow)  +&
             zdt*zncdetr(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorric1(1:kproma,krow)    = dcorric1(1:kproma,krow)    +&
             zdt*zcorric1(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorric2(1:kproma,krow)    = dcorric2(1:kproma,krow)    +&
             zdt*zcorric2(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
!>>UP: new diags for activation correction terms
        IF (lnewdiags) THEN
            dcorrcd2unphys_2d(1:kproma,krow) = dcorrcd2unphys_2d(1:kproma,krow) +&
zdt*zcorrcd2unphys_2d(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
            dcorrcd1_1(1:kproma,krow) = dcorrcd1_1(1:kproma,krow) +&
zdt*zcorrcd1_1(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
            dcorrcd1_2(1:kproma,krow) = dcorrcd1_2(1:kproma,krow) +&
zdt*zcorrcd1_2(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
            dcorrcd1_3(1:kproma,krow) = dcorrcd1_3(1:kproma,krow) +&
zdt*zcorrcd1_3(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
            dcorrcd1_4(1:kproma,krow) = dcorrcd1_4(1:kproma,krow) +&
zdt*zcorrcd1_4(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
            dcorrcd2_1(1:kproma,krow) = dcorrcd2_1(1:kproma,krow) +&
zdt*zcorrcd2_1(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
            dcorrcd2_2(1:kproma,krow) = dcorrcd2_2(1:kproma,krow) +&
zdt*zcorrcd2_2(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
!UP #783.3
            dcorrcd2_3(1:kproma,krow) = dcorrcd2_3(1:kproma,krow) +&
zdt*zcorrcd2_3(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
            dcorric1_1(1:kproma,krow) = dcorric1_1(1:kproma,krow) +&
zdt*zcorric1_1(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
!UP #783.3
            dcorric1_2(1:kproma,krow) = dcorric1_2(1:kproma,krow) +&
zdt*zcorric1_2(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
            dcorric1_3(1:kproma,krow) = dcorric1_3(1:kproma,krow) +&
zdt*zcorric1_3(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
!UP #783.3
        ENDIF
!<<UP: new diags
! UP comment: below needs to be commented out, see issue #767
!        dcorric3(1:kproma,krow)    = dcorric3(1:kproma,krow)    +&
!             zdt*zcorric3(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorric4(1:kproma,krow)    = dcorric4(1:kproma,krow)    + zdt*zcorric4(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorric5(1:kproma,krow)    = dcorric5(1:kproma,krow)    + zdt*zcorric5(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorric6(1:kproma,krow)    = dcorric6(1:kproma,krow)    + zdt*zcorric6(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorric7(1:kproma,krow)    = dcorric7(1:kproma,krow)    + zdt*zcorric7(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd1(1:kproma,krow)    = dcorrcd1(1:kproma,krow)    +&
             zdt*zcorrcd1(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd2(1:kproma,krow)    = dcorrcd2(1:kproma,krow)    +&
             zdt*zcorrcd2(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd3(1:kproma,krow)    = dcorrcd3(1:kproma,krow)    +&
             zdt*zcorrcd3(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd4(1:kproma,krow)    = dcorrcd4(1:kproma,krow)    +&
             zdt*zcorrcd4(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd5(1:kproma,krow)    = dcorrcd5(1:kproma,krow)    +&
             zdt*zcorrcd5(1:kproma,jk)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd6(1:kproma,krow)    = dcorrcd6(1:kproma,krow)    + zdt*zcorrcd6(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd7(1:kproma,krow)    = dcorrcd7(1:kproma,krow)    + zdt*zcorrcd7(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd8(1:kproma,krow)    = dcorrcd8(1:kproma,krow)    + zdt*zcorrcd8(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd9(1:kproma,krow)    = dcorrcd9(1:kproma,krow)    + zdt*zcorrcd9(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrcd10(1:kproma,krow)   = dcorrcd10(1:kproma,krow)   + zdt*zcorrcd10(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dcorrxl(1:kproma,krow)     = dcorrxl(1:kproma,krow)     + zdt*zcorrxl(1:kproma)*zdpg(1:kproma,jk)
        dcorrxl2(1:kproma,krow)    = dcorrxl2(1:kproma,krow)    + zdtime*zcorrxl2(1:kproma)*zdpg(1:kproma,jk)
        dcorrxi(1:kproma,krow)     = dcorrxi(1:kproma,krow)     + zdt*zcorrxi(1:kproma)*zdpg(1:kproma,jk)
        dcorrxi2(1:kproma,krow)    = dcorrxi2(1:kproma,krow)    + zdtime*zcorrxi2(1:kproma)*zdpg(1:kproma,jk)
        dxlte(1:kproma,krow)       = dxlte(1:kproma,krow)       + zdtime*zactxlte(1:kproma,jk)*zdpg(1:kproma,jk)
        dxite(1:kproma,krow)       = dxite(1:kproma,krow)       + zdtime*zactxite(1:kproma,jk)*zdpg(1:kproma,jk)
        dxttecdnc(1:kproma,krow)   = dxttecdnc(1:kproma,krow)   + zdtime*zxttecdnc(1:kproma,jk)*zdpg(1:kproma,jk)
        dxtteicnc(1:kproma,krow)   = dxtteicnc(1:kproma,krow)   + zdtime*zxtteicnc(1:kproma,jk)*zdpg(1:kproma,jk)
        dspr(1:kproma,krow)        = dspr(1:kproma,krow)        + zdt*zspr(1:kproma)*zdpg(1:kproma,jk)
        drpr(1:kproma,krow)        = drpr(1:kproma,krow)        + zdt*zrpr(1:kproma)*zdpg(1:kproma,jk)
        dmltn2(1:kproma,krow)      = dmltn2(1:kproma,krow)      - zdt*zmltn2(1:kproma)*zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk)
        dtestCD(1:kproma,krow)     = dtestCD(1:kproma,krow)     + zdt*ztestcd(1:kproma,jk)*zdpg(1:kproma,jk)
        dtestIC(1:kproma,krow)     = dtestIC(1:kproma,krow)     + zdt*ztestic(1:kproma,jk)*zdpg(1:kproma,jk)
        dtestLWC(1:kproma,krow)    = dtestLWC(1:kproma,krow)    + zdt*ztestlwc(1:kproma,jk)*zdpg(1:kproma,jk)
        dtestIWC(1:kproma,krow)    = dtestIWC(1:kproma,krow)    + zdt*ztestiwc(1:kproma,jk)*zdpg(1:kproma,jk)
        dgentl(1:kproma,krow)      = dgentl(1:kproma,krow)      + zdt*zgentl(1:kproma)*zdpg(1:kproma,jk)
        dgenti(1:kproma,krow)      = dgenti(1:kproma,krow)      + zdt*zgenti(1:kproma)*zdpg(1:kproma,jk)
        dsub_ice(1:kproma,krow)    = dsub_ice(1:kproma,krow)    + zdt*zxisub(1:kproma)*zdpg(1:kproma,jk)
        dmlt_snow(1:kproma,krow)   = dmlt_snow(1:kproma,krow)   + zdt*zsmlt(1:kproma)*zdpg(1:kproma,jk)
!<<DN: burden

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_diagnostics)
           CALL timer_start(timer_cmp_cmp)
           ! Timer needs to be started here so that it's on for the next
           ! iteration of the loop and also after finishing the loop
        END IF
        !<<UP

!--- wet-chemistry and in-cloud scavenging ----------------
END DO column_processes
        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_start(timer_cmp_8updating)
        END IF
        !<<UP
!>>DN: burden
       dsnow(1:kproma,krow)       = dsnow(1:kproma,krow)       + zdtime*pssfl(1:kproma)
!<<DN: burden

!>>UP FOR-ICE
!davidn SLF
    DO jl=1,kproma ! temperature interpolated to layer boundaries
       ! set temperature for lowest and highest level
       ztp1h(jl,1)=ptm1(jl,1)+ztmst*ptte(jl,1)
       ztp1h(jl,klevp1)=ptm1(jl,klev)+ztmst*ptte(jl,klev)
       DO jk=2,klev
          ! Loop through levels and get some pressure-weighted mean of the
          ! temperature for each level
          ztp1h(jl,jk)=((ptm1(jl,jk-1)+ztmst*ptte(jl,jk-1))*papm1(jl,jk-1)                    &
               *(papm1(jl,jk)-paphm1(jl,jk))                     &
               +(ptm1(jl,jk)+ztmst*ptte(jl,jk))*papm1(jl,jk)                         &
               *(paphm1(jl,jk)-papm1(jl,jk-1)))                  &
               /(paphm1(jl,jk)*(papm1(jl,jk)-papm1(jl,jk-1)))
          ! ztp1h is later used to judge the temperature bins
       END DO
    END DO
    
    DO jk=1,klev
       !UP question: Why are you overwriting picnc here?
       !UP note: because picnc isn't used for the SLF diagnostics, with which
       ! these two lines were implemented in r3578, but they are used for the
       ! diagnostics routine further down, and because I find no hint on what
       ! they are for, I interpret them as a mistake and take them out
       !picnc(1:kproma,jk)  = zrho(1:kproma,jk)*(pxtm1(1:kproma,jk,idt_icnc) + ztmst*pxtte(1:kproma,jk,idt_icnc))
       !picnc(1:kproma,jk)  = MAX(picnc(1:kproma,jk),cqtmin)

       ! Liquid water contents @t:
       zxlp1_slf(1:kproma) = pxlm1(1:kproma,jk)+ztmst*pxlte(1:kproma,jk)
       ! Ice water contents @t:
       zxip1_slf(1:kproma) = pxim1(1:kproma,jk)+ztmst*pxite(1:kproma,jk) 
       ! Total water contents @t:
       zxtotp1_slf(1:kproma) = zxlp1_slf(1:kproma) + zxip1_slf(1:kproma)
       ! Is there any ice or water present?
       ! ll1: there is ice or water present
       ! ll1: there is ice or water present
       ll1(1:kproma) = (zxtotp1_slf(1:kproma) > eps)
       ! ll1 is overwritten later, so we make a copy here
       ll2(:) = ll1(:)

       ztmp1(1:kproma) = 1000._dp*zxlp1_slf(1:kproma)*zdpg(1:kproma,jk)  ! lwp in g/m2
       ztmp2(1:kproma) = 1000._dp*zxip1_slf(1:kproma)*zdpg(1:kproma,jk)  ! iwp in g/m2

       IF (lslf) THEN
            !>>UP #844
            ! MASS AND NUMBER BURDENS DIVIDED BY TEMPERATURE REGIMES
            ll5(1:kproma) = ll2(1:kproma) .AND. &
               ! Temperature condition
               (ztp1h(1:kproma,jk) < tmelt) .AND. (ztp1h(1:kproma,jk) > cthomi)
            slfdiag_lwp_mxpT(1:kproma,krow) = slfdiag_lwp_mxpT(1:kproma,krow)+ &
                                  MERGE(ztmp1(1:kproma)*zdtime, 0._dp, ll5(1:kproma))
            slfdiag_cdnc_mxpT(1:kproma,krow) = slfdiag_cdnc_mxpT(1:kproma,krow)+ &
                                  MERGE(zcdnc(1:kproma,jk)*zdtime*zdz(1:kproma,jk), 0._dp, ll5(1:kproma))
            slfdiag_iwp_mxpT(1:kproma,krow) = slfdiag_iwp_mxpT(1:kproma,krow)+ &
                                  MERGE(ztmp2(1:kproma)*zdtime, 0._dp, ll5(1:kproma))
            slfdiag_icnc_mxpT(1:kproma,krow) = slfdiag_icnc_mxpT(1:kproma,krow)+ &
                                  MERGE(picnc(1:kproma,jk)*zdtime*zdz(1:kproma,jk), 0._dp, ll5(1:kproma))
            ll5(1:kproma) = ll2(1:kproma) .AND. &
               ! Temperature condition
               (ztp1h(1:kproma,jk) > tmelt)
            slfdiag_lwp_g0(1:kproma,krow) = slfdiag_lwp_g0(1:kproma,krow) + &
                                  MERGE(ztmp1(1:kproma)*zdtime, 0._dp, ll5(1:kproma))
            slfdiag_cdnc_g0(1:kproma,krow) = slfdiag_cdnc_g0(1:kproma,krow) + &
                                  MERGE(zcdnc(1:kproma,jk)*zdtime*zdz(1:kproma,jk), 0._dp, ll5(1:kproma))
            ll5(1:kproma) = ll2(1:kproma) .AND. &
               ! Temperature condition
               (ztp1h(1:kproma,jk) < cthomi)
            slfdiag_iwp_l35(1:kproma,krow) = slfdiag_iwp_l35(1:kproma,krow)+ &
                                  MERGE(ztmp2(1:kproma)*zdtime, 0._dp, ll5(1:kproma))
            slfdiag_icnc_l35(1:kproma,krow) = slfdiag_icnc_l35(1:kproma,krow)+ &
                                  MERGE(picnc(1:kproma,jk)*zdtime*zdz(1:kproma,jk), 0._dp, ll5(1:kproma))

            !<<UP #844
       ENDIF

       !>>UP
       IF (lnewdiags) THEN
            DO jm = 1, size(cdnc_bins)-1
                   DO jl = 1, size(lwc_bins)-1
                        ll1(1:kproma) = (zcdnc(1:kproma,jk) .GE. cdnc_bins(jm)) .AND. &
                                  (zcdnc(1:kproma,jk) < cdnc_bins(jm + 1)) .AND. &
                                  (zxlp1_slf(1:kproma) .GE. lwc_bins(jl)) .AND.  &
                                  (zxlp1_slf(1:kproma) < lwc_bins(jl + 1))
                        cdnc_lwc_binned(jm*jl)%ptr(1:kproma,jk,krow) = cdnc_lwc_binned(jm*jl)%ptr(1:kproma,jk,krow) + &
                                                          MERGE(1._dp * delta_time, 0._dp, ll1(1:kproma))
                   ENDDO
            ENDDO
            DO jm = 1, size(cdnc_bins)-1
                        ll1(1:kproma) = (zcdnc(1:kproma,jk) .GE. cdnc_bins(jm)) .AND. &
                                  (zcdnc(1:kproma,jk) < cdnc_bins(jm + 1))
                        cdnc_binned(jm)%ptr(1:kproma,jk,krow) = cdnc_binned(jm)%ptr(1:kproma,jk,krow) + &
                                                          MERGE(1._dp * delta_time, 0._dp, ll1(1:kproma))
            ENDDO
            d_lwc_lifetime_sources(1:kproma,krow) = d_lwc_lifetime_sources(1:kproma,krow) +   &
                        MERGE(zdt*(zxlp1_slf(1:kproma) / (z_lwc_sources(1:kproma,jk)*1000._dp*zdpg(1:kproma,jk))), &
                               0._dp, z_lwc_sources(1:kproma,jk) > eps)
            d_lwc_lifetime_sinks(1:kproma,krow) = d_lwc_lifetime_sinks(1:kproma,krow) +       &
                        MERGE(zdt*(zxlp1_slf(1:kproma) / (z_lwc_sinks(1:kproma,jk)*1000._dp*zdpg(1:kproma,jk))), &
                               0._dp, z_lwc_sinks(1:kproma,jk) > eps)
            d_iwc_lifetime_sources(1:kproma,krow) = d_iwc_lifetime_sources(1:kproma,krow) +   &
                        MERGE(zdt*(zxip1_slf(1:kproma) / (z_iwc_sources(1:kproma,jk)*1000._dp*zdpg(1:kproma,jk))), &
                               0._dp, z_iwc_sources(1:kproma,jk) > eps)
            d_iwc_lifetime_sinks(1:kproma,krow) = d_iwc_lifetime_sinks(1:kproma,krow) +       &
                        MERGE(zdt*(zxip1_slf(1:kproma) / (z_iwc_sinks(1:kproma,jk)*1000._dp*zdpg(1:kproma,jk))), &
                               0._dp, z_iwc_sinks(1:kproma,jk) > eps)
            d_cdnc_lifetime_sources(1:kproma,krow) = d_cdnc_lifetime_sources(1:kproma,krow) + &
                        MERGE(zdt*(zcdnc(1:kproma,jk) / (z_cdnc_sources(1:kproma,jk)*1000._dp     &
                                                   *zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk))), &
                               0._dp, z_cdnc_sources(1:kproma,jk) > eps)
            d_cdnc_lifetime_sinks(1:kproma,krow) = d_cdnc_lifetime_sinks(1:kproma,krow) +     &
                        MERGE(zdt*(zcdnc(1:kproma,jk) / (z_cdnc_sinks(1:kproma,jk)*1000._dp       &
                                                   *zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk))), &
                               0._dp, z_cdnc_sinks(1:kproma,jk) > eps)
            d_icnc_lifetime_sources(1:kproma,krow) = d_icnc_lifetime_sources(1:kproma,krow) + &
                        MERGE(zdt*(picnc(1:kproma,jk) / (z_icnc_sources(1:kproma,jk)*1000._dp     &
                                                   *zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk))), &
                               0._dp, z_icnc_sources(1:kproma,jk) > eps)
            d_icnc_lifetime_sinks(1:kproma,krow) = d_icnc_lifetime_sinks(1:kproma,krow) +     &
                        MERGE(zdt*(picnc(1:kproma,jk) / (z_icnc_sinks(1:kproma,jk)*1000._dp       &
                                                   *zdpg(1:kproma,jk)*zrho_rcp(1:kproma,jk))), &
                               0._dp, z_icnc_sinks(1:kproma,jk) > eps)
       ENDIF
       !<<UP

    END DO
!davidn
    IF (lnewdiags) THEN
            DO jm = 1, size(prcp_bins)-1
                        ll1(1:kproma) = (prsfl(1:kproma) .GE. prcp_bins(jm)) .AND. &
                                  (prsfl(1:kproma) < prcp_bins(jm + 1))
                        prcp_binned(jm)%ptr(1:kproma,krow) = prcp_binned(jm)%ptr(1:kproma,krow) + &
                                                          MERGE(1._dp * delta_time, 0._dp, ll1(1:kproma))
            ENDDO
    ENDIF
!<<UP FOR-ICE

 !--- save stratiform precipitation for dust emissions
  CALL set_vphysc_var(kproma, -1, krow, prflstrat=prsfl, psflstrat=pssfl)
 
        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_8updating)
           CALL timer_stop(timer_cmp_cmp)
           CALL timer_start(timer_cmp_9wetchem)
        END IF
        !<<UP

!--- wet-chemistry and in-cloud scavenging ----------------
!
!       9.    Wet chemistry and in-cloud scavenging

  CALL cloud_subm_2( &

          !-- IN   
          kproma, kbdim, klev, ktdia, krow, & 
          zmlwc(:,:), zmiwc(:,:), zmratepr(:,:), zmrateps(:,:), &
          zfrain(:,:), zfsnow(:,:), zfevapr(:,:), zfsubls(:,:), & 
          zmsnowacl(:,:), paclc(:,:), ptm1(:,:), ptte(:,:), pxtm1(:,:,:), & 
          !-- INOUT
          pxtte(:,:,:), &
          !-- IN
          paphp1(:,:), papp1(:,:), zrho(:,:), zclcpre_2d(:,:) ) 

!SF Note: I'd prefer to reorder the arguments, so that the intent(inout) is after all the intent(in),
!         but this implies some code change in echam code, so I leave that off for the moment.

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_9wetchem)
           CALL timer_start(timer_cmp_diagnostics)
        END IF
        !<<UP

!       10.    Diagnostics
!
!SF note: all the following is pure diagnostics, in the sense that the value of the calculated variables below
!         *does not* affect the rest of the physics and dynamics. In other words, disabling the following
!         yields a strictly identical climate as compared to not disabling it.

  IF (lctrl_diags) THEN 

     DO jk=ktdia,klev
   
           !SF ToDo: should be converted to a 2D subroutine
           CALL diagnostics( &
   
                   !-- IN
                   kbdim, kproma, &
                   itop(:,jk), jk, &
                   zcdnc(:,jk), picnc(:,jk), paclc(:,jk), zdpg(:,jk), &
                   zdz(:,jk), zfrln(:,jk), zrho(:,jk), &
                   zrprn(:,jk), zsacln(:,jk), zxib(:,jk), zxlb(:,jk), &
                   ztp1tmp(:,jk), reffl(:,jk,krow), reffi(:,jk,krow), &
                   ll_liqcl(:,jk), ll_icecl(:,jk), &
                   !-- INOUT
                   cdnc(:,jk,krow), cdnc_acc(:,jk,krow), cdnc_burden(:,krow), &
                   cdnc_ct(:,krow), cliwc_time(:,jk,krow), cloud_time(:,jk,krow), &
                   icnc(:,jk,krow), icnc_acc(:,jk,krow), icnc_burden(:,krow), &
                   iwc_acc(:,jk,krow), iwp_tovs(:,krow), lwc_acc(:,jk,krow), &
                   qacc(:,jk,krow), qaut(:,jk,krow), qfre(:,jk,krow), &
                   reffi_acc(:,jk,krow), reffi_time(:,krow), reffi_tovs(:,krow), &
                   reffl_acc(:,jk,krow), reffl_ct(:,krow), reffl_time(:,krow), &
                   zcdnc_burden(:), zicnc_burden(:), ztau1i(:), &
                   zreffct(:), paclcac(:,jk) )
   
     ENDDO

     cloud_cover_duplic(1:kproma,:,krow) = paclc(1:kproma,:) !SF #462: renamed var and moved its setting in the diags
   
     prelhum(1:kproma,:) = pqm1(1:kproma,:)/zqsm1(1:kproma,:)
     prelhum(1:kproma,:) = MAX(MIN(prelhum(1:kproma,:),1._dp),0._dp)

   !--- Included for prognostic CDNC/IC scheme ----------------------------
     ll1(1:kproma) = (zcdnc_burden(1:kproma) > epsec)
   
     ztmp1(1:kproma) = cdnc_burden_acc(1:kproma,krow) + zdtime*zcdnc_burden(1:kproma)
     ztmp2(1:kproma) = burden_time(1:kproma,krow)     + zdtime
   
     cdnc_burden_acc(1:kproma,krow) = MERGE(ztmp1(1:kproma), cdnc_burden_acc(1:kproma,krow), ll1(1:kproma))
     burden_time(1:kproma,krow)     = MERGE(ztmp2(1:kproma), burden_time(1:kproma,krow), ll1(1:kproma))
   
     ll1(1:kproma) = (zicnc_burden(1:kproma) > epsec)
   
     ztmp1(1:kproma) = icnc_burden_acc(1:kproma,krow) + zdtime*zicnc_burden(1:kproma)
     ztmp2(1:kproma) = burdic_time(1:kproma,krow)     + zdtime
   
     icnc_burden_acc(1:kproma,krow) = MERGE(ztmp1(1:kproma), icnc_burden_acc(1:kproma,krow), ll1(1:kproma))
     burdic_time(1:kproma,krow)     = MERGE(ztmp2(1:kproma), burdic_time(1:kproma,krow), ll1(1:kproma))
   
   !--- End included for CDNC/IC scheme -----------------------------------
   
   !       10.1   Accumulated precipitation at the surface
   !
      paprl(1:kproma) = paprl(1:kproma) + zdtime * (prsfl(1:kproma)+pssfl(1:kproma))
      paprs(1:kproma) = paprs(1:kproma) + zdtime * pssfl(1:kproma)
   
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
   !
   !       10.3   Vertical integrals of humidity, cloud water and cloud ice
   !
       zqvi(1:kproma)  = 0.0_dp
       zxlvi(1:kproma) = 0.0_dp
       zxivi(1:kproma) = 0.0_dp
   !
       DO 933 jk = ktdia,klev
          zqvi(1:kproma)  = zqvi(1:kproma)  + pqm1(1:kproma,jk) *zdpg(1:kproma,jk)
          zxlvi(1:kproma) = zxlvi(1:kproma) + pxlm1(1:kproma,jk)*zdpg(1:kproma,jk) 
          zxivi(1:kproma) = zxivi(1:kproma) + pxim1(1:kproma,jk)*zdpg(1:kproma,jk)
       933 END DO
   !
       pqvi(1:kproma)  = pqvi(1:kproma)  + zdtime*zqvi(1:kproma)
       pxlvi(1:kproma) = pxlvi(1:kproma) + zdtime*zxlvi(1:kproma)
       pxivi(1:kproma) = pxivi(1:kproma) + zdtime*zxivi(1:kproma)

    ENDIF !lctrl_diags

        !>>UP timers, #777
        IF (ltimer) THEN
           CALL timer_stop(timer_cmp_diagnostics)
           CALL timer_stop(timer_cmp_inside)
        END IF
        !<<UP

  RETURN

END SUBROUTINE cloud_micro_interface

SUBROUTINE melting_snow_and_ice( &
              !-- IN
              kbdim, kproma, &
              ld_mlt, &
              ptm1, pxim1, pdp, picncq, plsdcp, plvdcp, pcdnc, picnc, &
              pjk, pkrow, & !UP #821
              pcmpsimpl_mlt, & !UP #821
              !-- INOUT
              picnc_upd, pqmel, pcdnc_upd, prfl, psfl, pxiflux, pxifluxn, pxite, &
              ! OUT
              pimlt, psmlt, pximlt, pmltn, pmltn2 )!davidn

  !SF *Important note:*
  !   In this subroutine, most of the output variables are changed without enforcing
  !   explicitely ld_mlt to be true, even though the call to this subroutine is conditioned by 
  !   'IF (ANY(ll_mlt(1:kproma)))'
  !   This is ok and *DO NOT* generate any side-effect (with nproma dependency), because by construction 
  !   all calculations are only effective when ld_mlt is true.
  !   I therefore don't enforce any strict formulation of the form:
  !   output_var = MERGE(some_new_value, output_var, ld_mlt)
  !   BEWARE to keep this in mind in future, if any new development is made here!

  !>>UP #821
  USE mo_cmp_diagn,           ONLY: diag_delta_mlt, diag_cntr_mlt
  !<<UP #821

  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL, INTENT(in) :: ld_mlt(kbdim) !< switch that traces if temperature is above melting point

  REAL(dp), INTENT(in) :: ptm1(kbdim)   !< temperature (t-1)
  REAL(dp), INTENT(in) :: pxim1(kbdim)  !< cloud ice (t-1) 
  REAL(dp), INTENT(in) :: pdp(kbdim)    !< pressure difference of the layer [Pa]
  REAL(dp), INTENT(in) :: picncq(kbdim) !< Temporary value for ICNC [1/m3]
  REAL(dp), INTENT(in) :: plsdcp(kbdim) !< latent heat of sublimation div. by the specific heat at constant pres.
  REAL(dp), INTENT(in) :: plvdcp(kbdim) !< latent heat of vaporization div. by the specific heat at constant pres.
  REAL(dp), INTENT(in) :: pcdnc(kbdim)  !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp), INTENT(in) :: picnc(kbdim)  !< Ice crystal number concentration (ICNC) [1/m3]
  !>>UP #821
  INTEGER, INTENT(in)  :: pjk !jk for use in this subroutine
  INTEGER, INTENT(in)  :: pkrow ! krow for use in this subroutine
  REAL(dp), INTENT(in) :: pcmpsimpl_mlt(kbdim) !< input for prescribed mlt
  !<<UP #821

  REAL(dp), INTENT(inout) :: picnc_upd(kbdim) !< Ice crystal number concentration (ICNC) [1/m3]
  !SFNote:  intent(inout) for picdc_upd is important, because the corresponding argument of the caller
  !         is a pointer, which may point to the argument corresponding to dummy 'pcdnc'
  !         Since fortran is passing arguments by ref, an intent(out) in this
  !         case is not possible
  REAL(dp), INTENT(inout) :: pqmel(kbdim)    !< cloud droplet source rate from melting ice [m-3 s-1]
  REAL(dp), INTENT(inout) :: pcdnc_upd(kbdim) !< Updated cloud droplet number concentration (CDNC) [1/m3]
  !SFNote:  intent(inout) for pcdnc_upd is important, because the corresponding argument of the caller
  !         is a pointer, which may point to the argument corresponding to dummy 'pcdnc'
  !         Since fortran is passing arguments by ref, an intent(out) in this
  !         case is not possible
  REAL(dp), INTENT(inout) :: prfl(kbdim)     !< rain flux [kg/m2/s] 
  REAL(dp), INTENT(inout) :: psfl(kbdim)     !< snow flux [kg/m2/s]
  REAL(dp), INTENT(inout) :: pxiflux(kbdim)  !< flux of ice crystals falling into the grid box from above
  REAL(dp), INTENT(inout) :: pxifluxn(kbdim) !< flux of ice crystals number falling into the grid box from above
  REAL(dp), INTENT(inout) :: pxite(kbdim)    !< tendency of cloud ice

  REAL(dp), INTENT(out) :: pimlt(kbdim)  !< melting of ice if T>273 K [kg/kg]
  REAL(dp), INTENT(out) :: psmlt(kbdim)  !< melting of snow [kg/kg]
  REAL(dp), INTENT(out) :: pximlt(kbdim) !< melting of the ice falling from above [kg/kg]

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim)
  LOGICAL  :: ll1(kbdim)
  REAL(dp) :: pmltn(kbdim),pmltn2(kbdim)!davidn

  ztmp1(1:kproma) = ptm1(1:kproma) - tmelt
 
  ztmp1(1:kproma) = MAX(0._dp, ztmp1(1:kproma))   !SF ztdif
  ztmp1(1:kproma) = zcons2*ztmp1(1:kproma)*pdp(1:kproma) &
                     / ( plsdcp(1:kproma)-plvdcp(1:kproma) ) !SF zcons*ztdif
  
  ztmp2(1:kproma) = xsec*psfl(1:kproma)
  !>>UP #821
  IF (lcmpsimpl_prescr) THEN
        SELECT CASE(ncmpsimpl_prescr_mlt)
            !CASE(0)
                ! nothing happens, ztmp1 stays ztmp1
            CASE(1,2)
                ztmp1(1:kproma) = MERGE(MAX(0._dp,pcmpsimpl_mlt(1:kproma)), 0._dp, ld_mlt(1:kproma))
        END SELECT
  ELSE IF (lemuphase_mlt) THEN
  !<<UP #821
  !>>UP #802
          ! Phase melting of snow to rain here -> affects prfl, psfl, psmlt
          ! Phase melting of sedimenting ice here -> affects pxiflux, ztmp3 (->
          ! pxifluxn), pximlt below where ztmp1 is used a second time
                ztmp1(1:kproma) = eta_emu_mlt * ztmp1(1:kproma)
  ENDIF
  !<<UP #802
  !>>UP #821
  diag_delta_mlt(1:kproma,pjk,pkrow) = diag_delta_mlt(1:kproma,pjk,pkrow) + &
                                     MERGE(ztmp1(1:kproma), 0._dp, ld_mlt(1:kproma))
  ! Should I merge with ld_mlt here? No, because for the melting of sedimenting
  ! ice and snow that we care about, it is not used (and IC melting is not
  ! phased)
  diag_cntr_mlt(1:kproma,pjk,pkrow) = diag_cntr_mlt(1:kproma,pjk,pkrow) + &
                                     MERGE(1._dp, 0._dp, ld_mlt(1:kproma))
  !<<UP #821
  ! original line
  ztmp2(1:kproma) = MIN(ztmp2(1:kproma), ztmp1(1:kproma)) !SFzximelt
     
  prfl(1:kproma)  = prfl(1:kproma) + ztmp2(1:kproma)
  psfl(1:kproma)  = psfl(1:kproma) - ztmp2(1:kproma)

  psmlt(1:kproma) = ztmst*grav*ztmp2(1:kproma) / pdp(1:kproma)

  ztmp2(1:kproma) = xsec*pxiflux(1:kproma)
  ! original line
  ztmp2(1:kproma) = MIN(ztmp2(1:kproma), ztmp1(1:kproma))

  ll1(1:kproma) = (pxiflux(1:kproma) > epsec)

  ztmp3(1:kproma) = pxifluxn(1:kproma)*ztmp2(1:kproma)/MAX(pxiflux(1:kproma),epsec)
  ztmp3(1:kproma) = MERGE(ztmp3(1:kproma), 0._dp, ll1(1:kproma)) !SF zxinmelt

  pxiflux(1:kproma)  = pxiflux(1:kproma)  - ztmp2(1:kproma)
  pxifluxn(1:kproma) = pxifluxn(1:kproma) - ztmp3(1:kproma)

  pxifluxn(:) = consistency_number_to_mass(kbdim, kproma, epsec, pxiflux(:), pxifluxn(:))

  pximlt(1:kproma) = ztmst*grav*ztmp2(1:kproma) / pdp(1:kproma)
  !UP note: this pximlt is later added to the CD mass, but the number of ICs
  !that is melting is never added anywhere. This is a problem as it means that
  !the budget is not closed. However, DN and SF think that this will contribute
  !only a small amount of imbalance and we therefore don't fix it now. 

  ! If T > tmelt melt all ice crystals and transfer to cloud droplets 
  !UP note: this is an easy formulation that does not make sense to be phased

  ztmp1(1:kproma) = pxim1(1:kproma) + ztmst*pxite(1:kproma)
  ztmp1(1:kproma) = MAX(0._dp, ztmp1(1:kproma))
  pimlt(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_mlt(1:kproma))
  pxite(1:kproma) = pxite(1:kproma) - ztmst_rcp*pimlt(1:kproma)

  ztmp1(1:kproma)     = MERGE(picncq(1:kproma), 0._dp, ld_mlt(1:kproma))
  pcdnc_upd(1:kproma) = pcdnc(1:kproma) + ztmp1(1:kproma)
!>>DN: new diags
  pmltn2(1:kproma) = MERGE(picnc(1:kproma)-icemin,0._dp,ld_mlt(1:kproma)) 
!<<DN: new diags

  picnc_upd(1:kproma) = MERGE(icemin, picnc(1:kproma), ld_mlt(1:kproma)) 
!>>SF Kasja diags
  pmltn(1:kproma) = ztmp1(1:kproma)
!<<SF Kasja diags

  pqmel(1:kproma) = pqmel(1:kproma) + zdt*ztmp1(1:kproma)

END SUBROUTINE melting_snow_and_ice

SUBROUTINE sublimation_snow_and_ice_evaporation_rain( &
              !-- IN
              kbdim, kproma, &
              ld_precip, ld_falling_ice, &
              pqm1, ptm1, pclcpre, pdp, pdpg, picesub,   &
              plsdcp, pqrho, pqsi, prho_rcp, psfl, prho, &
              pqsw, prfl, psusatw_evap, pastbstw, pclcfi, &
              pjk, pkrow, & !UP #821
              pcmpsimpl_subfi, pcmpsimpl_subfs, & !UP #821
              !-- INOUT
              pxiflux, pxifluxn, &
              !-- OUT
              pxisub, psub, pevp)
!DN: Please note that this subroutine handles sublimation of falling ice only; 
!    sublimation of cloud ice is done in 4.

  !>>UP #821
  USE mo_cmp_diagn,           ONLY: diag_delta_subfs, diag_cntr_subfs, &
                                    diag_delta_subfi, diag_cntr_subfi
  !<<UP #821
  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL, INTENT(in) :: ld_precip(kbdim)      !< switch that traces the presence of precipitation
  LOGICAL, INTENT(in) :: ld_falling_ice(kbdim) !< switch that traces the presence of falling ice

  REAL(dp), INTENT(in) :: pqm1(kbdim)         !< specific humidity (t-1)
  REAL(dp), INTENT(in) :: ptm1(kbdim)         !< temperature (t-1)
  REAL(dp), INTENT(in) :: pclcpre(kbdim)      !< fraction of grid box covered by precip
  REAL(dp), INTENT(in) :: pdp(kbdim)          !< pressure difference of the layer [Pa]
  REAL(dp), INTENT(in) :: pdpg(kbdim)         !< delta p over g [kg/m2]
  REAL(dp), INTENT(in) :: picesub(kbdim)      !< Subsat. w.r.t. ice
  REAL(dp), INTENT(in) :: plsdcp(kbdim)       !< latent heat of subl. div. by the specific heat at constant pres.
  REAL(dp), INTENT(in) :: pqrho(kbdim)        !< Inverse air density [m3/kg]
  REAL(dp), INTENT(in) :: pqsi(kbdim)         !< Saturation specific humidity w.r.t. ice [kg/kg]
  REAL(dp), INTENT(in) :: prho_rcp(kbdim)     !< Inverse air density
  REAL(dp), INTENT(in) :: psfl(kbdim)         !< snow flux [kg/m2/s]
  REAL(dp), INTENT(in) :: prho(kbdim)         !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: pqsw(kbdim)         !< Saturation specific humidity w.r.t. water (t-1) [kg/kg]    
  REAL(dp), INTENT(in) :: prfl(kbdim)         !< rain flux [kg/m2/s]
  REAL(dp), INTENT(in) :: psusatw_evap(kbdim) !< Subsat. w.r.t. water
  REAL(dp), INTENT(in) :: pastbstw(kbdim)     !< Thermodynamic term needed for water nucleation
  REAL(dp), INTENT(in) :: pclcfi(kbdim)       !< Fraction of grid box covered by falling ice
  !>>UP #821
  INTEGER, INTENT(in) :: pjk !jk for use in this subroutine
  INTEGER, INTENT(in) :: pkrow ! krow for use in this subroutine
  REAL(dp), INTENT(in) :: pcmpsimpl_subfi(kbdim)    !< input for prescribed subfi
  REAL(dp), INTENT(in) :: pcmpsimpl_subfs(kbdim)    !< input for prescribed subfs
  !<<UP #821

  REAL(dp), INTENT(inout) :: pxiflux(kbdim)  !< flux of ice crystals falling into the grid box from above
  REAL(dp), INTENT(inout) :: pxifluxn(kbdim) !< flux of ice crystals number falling into the grid box from above

  REAL(dp), INTENT(out) :: pxisub(kbdim) !< sublimation of cloud ice [kg/kg]
  REAL(dp), INTENT(out) :: psub(kbdim)   !< sublimation of snow [kg/kg]
  REAL(dp), INTENT(out) :: pevp(kbdim)   !< evaporation of rain [kg/kg]

  !-- Local vars
  REAL(dp) :: zclcpre(kbdim) !< surrogate for pclcpre, with a dummy value (1._dp) replacing zeros to avoid div by zero
  REAL(dp) :: zclcfi(kbdim)  !< surrogate for pclfi, with a dummy value (1._dp) replacing zeros to avoid div by zero

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim), ztmp4(kbdim)
  LOGICAL  :: ll1(kbdim)

  !SF ToDo: this subroutine needs to be generalized by a function applicable to snow, falling ice and rain,
  !         to avoid current code duplication

  ztmp1(1:kproma) = 1._dp/(2.43e-2_dp*rv)*plsdcp(1:kproma)**2/ptm1(1:kproma)**2
  ztmp1(1:kproma) = ztmp1(1:kproma) &
                  + 1._dp/0.211e-4_dp*prho_rcp(1:kproma)/pqsi(1:kproma)
  ztmp1(1:kproma) = 3.e6_dp*2._dp*pi*picesub(1:kproma)*prho_rcp(1:kproma)/ztmp1(1:kproma) !zcoeff

  zclcpre(1:kproma) = MERGE(pclcpre(1:kproma), 1._dp, ld_precip(1:kproma))
  zclcfi(1:kproma)  = MERGE(pclcfi(1:kproma),  1._dp, ld_falling_ice(1:kproma))

  !Snow:  
  ll1(1:kproma) = (psfl(1:kproma) > cqtmin) .AND. &
                  ld_precip(1:kproma)

  ztmp2(1:kproma) = zcons3*(MAX(psfl(1:kproma),cqtmin) / zclcpre(1:kproma))**(0.25_dp/1.16_dp)  !zclambs
  ztmp2(1:kproma) = 0.78_dp  *ztmp2(1:kproma)**2 & 
                  + 232.19_dp*pqrho(1:kproma)**0.25_dp * ztmp2(1:kproma)**2.625_dp     !zcfac4c
  ztmp2(1:kproma) = ztmp2(1:kproma) * ztmp1(1:kproma) * pdpg(1:kproma)                     

  ztmp3(1:kproma) = -xsec * psfl(1:kproma) / zclcpre(1:kproma)
  ztmp3(1:kproma) = MAX(ztmp3(1:kproma), ztmp2(1:kproma))                              !zzeps
  ztmp3(1:kproma) = -ztmst*ztmp3(1:kproma) / pdpg(1:kproma) * zclcpre(1:kproma)

  ztmp4(1:kproma) = xsec*(pqsi(1:kproma)-pqm1(1:kproma))
  ztmp4(1:kproma) = MAX(ztmp4(1:kproma),0._dp)

  ztmp3(1:kproma)  = MIN(ztmp3(1:kproma),ztmp4(1:kproma))
  !>>UP #802
  !>>UP #821
  IF (lcmpsimpl_prescr) THEN
        SELECT CASE(ncmpsimpl_prescr_subfis)
            CASE(0)
                ! original line
                ztmp3(1:kproma) = MAX(ztmp3(1:kproma),0._dp)
            CASE(1,2)
                ! ztmp3 is used only with ll1 later -> use ll1 everywhere here
                ! as well
                !psfl (computed later on from psub) shall never be negative from too much sublimation
                ztmp3(1:kproma) = MERGE(MIN(pcmpsimpl_subfs(1:kproma), psfl(1:kproma)/(zcons2 * pdp(1:kproma))), &
                                        0._dp, ll1(1:kproma))
        END SELECT
  !<<UP #821
  ELSE IF (lemuphase_sub_evp_falling) THEN
          ! phase sublimation of snow
          ! phase here because eta can be larger than 1 and then the maximal
          ! condition below (MIN) is needed to keep a check on it
          ztmp3(1:kproma)  = MAX(eta_emu_subfis*ztmp3(1:kproma),0._dp)
          !psfl (computed later on from psub) shall never be negative from too much sublimation
          ztmp3(1:kproma) = MIN(ztmp3(1:kproma), psfl(1:kproma)/(zcons2 * pdp(1:kproma))) 
  ELSE
          ! original line
          ztmp3(1:kproma) = MAX(ztmp3(1:kproma),0._dp)
  END IF 
  !<<UP #802
  !>>UP #821
  !UP TODO: doesn't below need the ll1 merge as well?
  diag_delta_subfs(1:kproma,pjk,pkrow) = diag_delta_subfs(1:kproma,pjk,pkrow) + &
                                  MERGE(ztmp3(1:kproma), 0._dp, ll1(1:kproma))
  diag_cntr_subfs(1:kproma,pjk,pkrow) = diag_cntr_subfs(1:kproma,pjk,pkrow) + &
                                         MERGE(1._dp, 0._dp, ll1(1:kproma))
  !<<UP #821
  psub(1:kproma)  = MERGE(ztmp3(1:kproma), 0._dp, ll1(1:kproma))

  !Ice:
  !SF Note: #459: replaced ld_precip and zclcpre by ld_falling_ice and zclcfi in the ice section
  ll1(1:kproma) = (pxiflux(1:kproma) > cqtmin) &
            .AND. ld_falling_ice(1:kproma)

  ztmp2(1:kproma) = zcons3*(MAX(pxiflux(1:kproma),cqtmin) / zclcfi(1:kproma))**(0.25_dp/1.16_dp) !zclambs
  ztmp2(1:kproma) = 0.78_dp  *ztmp2(1:kproma)**2 &
                  + 232.19_dp*pqrho(1:kproma)**0.25_dp * ztmp2(1:kproma)**2.625_dp   !zcfac4c
  ztmp2(1:kproma) = ztmp2(1:kproma) * ztmp1(1:kproma) * pdpg(1:kproma)

  ztmp3(1:kproma) = -xsec * pxiflux(1:kproma) / zclcfi(1:kproma)
  ztmp3(1:kproma) = MAX(ztmp3(1:kproma), ztmp2(1:kproma))                       !zzeps
  ztmp3(1:kproma) = -ztmst*ztmp3(1:kproma) / pdpg(1:kproma) * zclcfi(1:kproma)

  ztmp4(1:kproma) = xsec*(pqsi(1:kproma)-pqm1(1:kproma))
  ztmp4(1:kproma) = MAX(ztmp4(1:kproma),0._dp)

  ztmp3(1:kproma)  = MIN(ztmp3(1:kproma),ztmp4(1:kproma))
  !>>UP #802
  !>>UP #821
  IF (lcmpsimpl_prescr) THEN
        SELECT CASE(ncmpsimpl_prescr_subfis)
            CASE(0)
                ztmp3(1:kproma)  = MAX(ztmp3(1:kproma),0._dp)
            CASE(1,2)
                ! ztmp3 is used only with ll1 after -> use ll1 here everywhere
                ! as well
                !pxiflux (computed later from pxisub, which in turn is computed from
                !ztmp3) shall never be negative from too much sublimation
                ztmp3(1:kproma) = MERGE(MIN(pcmpsimpl_subfi(1:kproma), pxiflux(1:kproma)/(zcons2 * pdp(1:kproma))), &
                                        0._dp, ll1(1:kproma))
        END SELECT
  !<<UP #821
  ELSE IF (lemuphase_sub_evp_falling) THEN
          ! phase sublimation of falling ice; this automatically propagates to
          ! pxiflux and pxifluxn
          ztmp3(1:kproma)  = MAX(eta_emu_subfis*ztmp3(1:kproma),0._dp)
          !pxiflux (computed later from pxisub, which in turn is computed from
          !ztmp3) shall never be negative from too much sublimation
          ztmp3(1:kproma) = MIN(ztmp3(1:kproma), pxiflux(1:kproma)/(zcons2 * pdp(1:kproma)))
  ELSE
          ! original line
          ztmp3(1:kproma)  = MAX(ztmp3(1:kproma),0._dp)
  END IF 
  !<<UP #802
  !>>UP #821
  !UP TODO: doesn't below need the ll1 merge as well?
  diag_delta_subfi(1:kproma,pjk,pkrow) = diag_delta_subfi(1:kproma,pjk,pkrow) + &
                                  MERGE(ztmp3(1:kproma), 0._dp, ll1(1:kproma))
  diag_cntr_subfi(1:kproma,pjk,pkrow) = diag_cntr_subfi(1:kproma,pjk,pkrow) + &
                                         MERGE(1._dp, 0._dp, ll1(1:kproma))
  !<<UP #821

  pxisub(1:kproma) = MERGE(ztmp3(1:kproma), 0._dp, ll1(1:kproma))

  ztmp4(1:kproma) = pxisub(1:kproma) * pxifluxn(1:kproma) / MAX(pxiflux(1:kproma),cqtmin)  !SF zsubin
  ztmp4(1:kproma) = zcons2 * ztmp4(1:kproma) * pdp(1:kproma)
  ztmp4(1:kproma) = MERGE(ztmp4(1:kproma), 0._dp, ll1(1:kproma))

  pxifluxn(1:kproma) = pxifluxn(1:kproma) - ztmp4(1:kproma)

  pxiflux(1:kproma) = pxiflux(1:kproma) - zcons2*pxisub(1:kproma)*pdp(1:kproma)

  pxifluxn(:) = consistency_number_to_mass(kbdim, kproma, epsec, pxiflux(:), pxifluxn(:)) !SF #480

  !Rain
  ll1(1:kproma) = (prfl(1:kproma) > cqtmin) .AND. &
                  ld_precip(1:kproma)

  ztmp2(1:kproma) = 870._dp * psusatw_evap(1:kproma) * pdpg(1:kproma) &
                  * (MAX(prfl(1:kproma),cqtmin)/zclcpre(1:kproma))**0.61_dp &
                  / (SQRT(prho(1:kproma))*pastbstw(1:kproma))                 

  ztmp3(1:kproma) = -xsec*prfl(1:kproma)/zclcpre(1:kproma)
  ztmp3(1:kproma) = MAX(ztmp3(1:kproma), ztmp2(1:kproma))
  ztmp3(1:kproma) = -ztmst*ztmp3(1:kproma)*zclcpre(1:kproma)/pdpg(1:kproma)

  ztmp4(1:kproma) = xsec*(pqsw(1:kproma)-pqm1(1:kproma))
  ztmp4(1:kproma) = MAX(ztmp4(1:kproma), 0._dp)

  !>>UP #802
  IF (lemuphase_sub_evp_falling) THEN
          ! phase evaporation of rain
          ! phase here because eta can be larger than 1 and then the maximal
          ! condition below (MIN) is needed to keep a check on it
          ztmp3(1:kproma) = MIN(eta_emu_evpr*ztmp3(1:kproma),ztmp4(1:kproma))
          !prfl (computed later from pevp, which in turn is computed from
          !ztmp3 below) shall never be negative from too much sublimation
          ztmp3(1:kproma) = MIN(ztmp3(1:kproma), prfl(1:kproma)/(zcons2 * pdp(1:kproma)))
  ELSE
          ! original line
          ztmp3(1:kproma)  = MIN(ztmp3(1:kproma),ztmp4(1:kproma))
  END IF 
  !<<UP #802
  ztmp3(1:kproma) = MAX(ztmp3(1:kproma),0._dp)

  pevp(1:kproma) = MERGE(ztmp3(1:kproma), 0._dp, ll1(1:kproma))

END SUBROUTINE sublimation_snow_and_ice_evaporation_rain

SUBROUTINE sedimentation_ice( &
              !-- IN
              kbdim, kproma, &
              paclc, paaa, pdp, &
              prho, prho_rcp, picnc, &
              !>>UP #797
              picnc_diff_up, &
              pxip1_diff_up, & 
              pqm1, pqsi, &
              !<<UP #797
              !-- INOUT
              pxip1, picnc_upd, pxiflux, pxifluxn, pclcfi, &
              pxisub, & !UP #797
              !-- OUT
              pmrateps)

  INTEGER, INTENT(in) :: kbdim, kproma

  REAL(dp), INTENT(in) :: paclc(kbdim) !< cloud cover
  REAL(dp), INTENT(in) :: paaa(kbdim) !< Air density correction needed for the ice crystal fall velocity
  REAL(dp), INTENT(in) :: pdp(kbdim)    !< pressure difference of the layer [Pa]
  REAL(dp), INTENT(in) :: prho(kbdim) !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: prho_rcp(kbdim) !< Inverse air density
  REAL(dp), INTENT(in) :: picnc(kbdim)    !< Ice crystal number concentration (ICNC) [1/m3]
  !>>UP #797
  ! difference in in-cloud ICNC from the upper level due to sedimentation
  REAL(dp), INTENT(in) :: picnc_diff_up(kbdim)   
  ! difference in in-cloud IC mass from the upper level due to sedimentation
  REAL(dp), INTENT(in) :: pxip1_diff_up(kbdim)   
  REAL(dp), INTENT(in) :: pqm1(kbdim)         !< specific humidity (t-1)
  REAL(dp), INTENT(in) :: pqsi(kbdim)         !< Saturation specific humidity w.r.t. ice [kg/kg]
  !<<UP #797
  
  REAL(dp), INTENT(inout) :: pxip1(kbdim) !< Ice mass (grid-mean values)
  REAL(dp), INTENT(inout) :: picnc_upd(kbdim)!< Ice crystal number concentration (ICNC) [1/m3]
  !SF Note: intent(inout) important for picnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pxiflux(kbdim)  !< flux of ice crystals
                                             !< IN: falling into the grid box from above
                                             !< OUT: updated by flux from current level 
  REAL(dp), INTENT(inout) :: pxifluxn(kbdim) !< flux of ice crystals number
                                             !< IN: falling into the grid box from above
                                             !< OUT: updated by flux from current level 
  REAL(dp), INTENT(inout) :: pclcfi(kbdim)   !< fraction of grid box covered by sedimenting ice

  !>>UP #797
  REAL(dp), INTENT(inout) :: pxisub(kbdim) !< sublimation of cloud ice [kg/kg]
  !<<UP #797

  REAL(dp), INTENT(out) :: pmrateps(kbdim) !< Snow formation rate in cloudy part of the grid box [kg/kg]

  !local vars:
  REAL(dp) :: zalfased(kbdim)              !< Parameter needed for the ice crystal fall velocity
  REAL(dp) :: zbetased(kbdim)              !< Parameter needed for the ice crystal fall velocity
  REAL(dp) :: zicnc_gridmean(kbdim)        !< ice crystal number conc. (grid-mean values) [1/m3]
  REAL(dp) :: zicnc_gridmean_bf_sed(kbdim) !< duplicate of ice crystal number conc. (grid-mean values) [1/m3]
                                           !< to keep a record of the value of zicnc_gridmean before sedimentation 
  REAL(dp) :: zmmean(kbdim)     !< volume mean ice crystal radius [m]
  REAL(dp) :: zxifallmc(kbdim)  !< fall velocity of ice crystal mass [m/s]
  REAL(dp) :: zxifallnc(kbdim)  !< fall velocity of ice crystal number [m/s]
  REAL(dp) :: zxi_bf_sed(kbdim) !< duplicate of cloud ice mass mixing ratio [kg/kg]
                                !< to keep a record of the value of pxip1 before sedimentation
  REAL(dp) :: zxiflx_from_level(kbdim) !< Flux of sedimenting ice crystals from current level
  REAL(dp) :: zxi_delta(kbdim)         !< Sedimented ice mmr (grid mean) [kg/kg]
 
  !>>UP #797
  REAL(dp) :: zicnc_bf_sed(kbdim) !store icnc before the routine for use later
  ! in the routine
  !<<UP #797

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim), ztmp4(kbdim)
  LOGICAL  :: ll1(kbdim), ll2(kbdim)

!SF Note: I don't get why the whole subroutine uses grid means
!         it looks to me that everything could be performed as in-cloud
!         which would avoid back and forth conversions
!

  zxi_bf_sed(1:kproma) = pxip1(1:kproma)
!>>UP #797
  zicnc_bf_sed(1:kproma) = picnc(1:kproma)
!<<UP #797

!UP note: conversion to grid-mean
  zicnc_gridmean(1:kproma) = picnc(1:kproma)*paclc(1:kproma)
!UP note: setting to minimal icemin here; this may introduce new crystals,
!especially where there is no cloud, but it should be a small amount
  zicnc_gridmean(1:kproma) = MAX(zicnc_gridmean(1:kproma),icemin)

  zicnc_gridmean_bf_sed(1:kproma) = zicnc_gridmean(1:kproma)

!UP note: zmmean is supposed to be a radius [m]: units don't work out here? 
! They are multiplying density and mmr to get mass concentration per volume;
! dividing it by the number you get mass per number, not a radius
! mi is documented to be a mass of ICs as well, so maybe zmmean is just
! erroneously documented as a radius
  zmmean(1:kproma) = prho(1:kproma)*pxip1(1:kproma)/zicnc_gridmean(1:kproma)
  zmmean(1:kproma) = MAX(zmmean(1:kproma), mi)

!UP note: ri_vol_mean_1 and _2 say they are radii, but from their values
!(2.166e-9 and 4.264e-8) it seems more likely that they are masses like mi
  ll1(1:kproma) = (zmmean(1:kproma) < ri_vol_mean_1 )
  ll2(1:kproma) = (.NOT. ll1(1:kproma)) .AND. (zmmean(1:kproma) < ri_vol_mean_2 )

!UP note: IC fall differently fast because of their radii
  zalfased(1:kproma) = MERGE(alfased_1, alfased_2, ll1(1:kproma))
  zalfased(1:kproma) = MERGE(alfased_3, zalfased(1:kproma), ll2(1:kproma))

  zbetased(1:kproma) = MERGE(betased_1, betased_2, ll1(1:kproma))
  zbetased(1:kproma) = MERGE(betased_3, zbetased(1:kproma), ll2(1:kproma))

!UP note: same equation as given in the "Documentation of the ECHAM5 
! double-moment cloud scheme (cloud cdnc icnc.f90)" from 2008
  zxifallmc(1:kproma) = fall*zalfased(1:kproma) & 
                      * (zmmean(1:kproma)**zbetased(1:kproma))*paaa(1:kproma)

  !limit fall velocity to 1 cm/s - 2 m/s:
  zxifallmc(1:kproma) = MAX(0.001_dp,zxifallmc(1:kproma))
  zxifallmc(1:kproma) = MIN(2._dp,zxifallmc(1:kproma))
  
  !>>MA Issue#795 Lemuphase_fvic
  IF (lemuphase_fvic) THEN
          zxifallmc(1:kproma) =  eta_emu_fvic * zxifallmc(1:kproma)
  ENDIF
  !<<MA

!UP note: fall velocity of number equal to that of mass
  zxifallnc(1:kproma) = zxifallmc(1:kproma) 

!UP note: here mass and number zxifall are equal except for zxifallmc,
!zxifallnc, and the prho_rcp
  ztmp1(1:kproma) = ztmst*grav*zxifallmc(1:kproma)*prho(1:kproma)/pdp(1:kproma) !SF zal1
  ll1(1:kproma)   = (zxifallmc(1:kproma) > eps)
  ztmp2(1:kproma) = pxiflux(1:kproma)*prho_rcp(1:kproma)/MAX(zxifallmc(1:kproma), eps)  
  ztmp2(1:kproma) = MERGE(ztmp2(1:kproma), 0._dp, ll1(1:kproma)) !SF zal2
  
!UP note: here mass and number zxifall are equal except for the prho_rcp
  ztmp3(1:kproma) = grav*ztmst*zxifallnc(1:kproma)*prho(1:kproma)/pdp(1:kproma) !SF zal3
  ll1(1:kproma)   = (zxifallnc(1:kproma) > eps)
  ztmp4(1:kproma) = pxifluxn(1:kproma)/MAX(zxifallnc(1:kproma), eps)
  ztmp4(1:kproma) = MERGE(ztmp4(1:kproma), 0._dp, ll1(1:kproma)) !SF zal4

  pxip1(1:kproma) = pxip1(1:kproma) &
                   * EXP(-ztmp1(1:kproma)) + ztmp2(1:kproma)*(1._dp-EXP(-ztmp1(1:kproma)))
  zicnc_gridmean(1:kproma) = zicnc_gridmean(1:kproma) &
                   * EXP(-ztmp3(1:kproma)) + ztmp4(1:kproma)*(1._dp-EXP(-ztmp3(1:kproma)))

  ll1(1:kproma) = (paclc(1:kproma) > clc_min)

!UP note: convert back to cloudy values
  ztmp1(1:kproma) = zicnc_gridmean(1:kproma) / MAX(paclc(1:kproma), clc_min)
!UP note: this is where icnc is getting updated; is ztmp1 always smaller than
!what IC conc. was before? Is zicnc_gridmean?
  !>>UP #797
  IF (lsedfix) THEN
        ! don't let the ICNC increase more than what it has decreased in the
        ! level above
        ! only do this where it is decreasing in the level above
        ! When it is increasing due to sedimentation also in the level above it
        ! is not clear to me what the restriction for the level below should be
        ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), MIN(ztmp1(1:kproma), &
                zicnc_bf_sed(1:kproma) - picnc_diff_up(1:kproma)), &
                picnc_diff_up(1:kproma) > 0)
        ! If the upper level is losing ICs in sedimentation, the lower level may
        ! not gain more than what the upper level is losing
        ! Same for number, but here we need an explicit conversion to
        ! in-cloud values first, because only those are constrained from the
        ! decreases above
        ! 0. save old pxip1 for later use
        ztmp3(1:kproma) = pxip1(1:kproma)
        ! 1. convert to in-cloud
        ztmp2(1:kproma) = pxip1(1:kproma) / MAX(paclc(1:kproma), clc_min)
        ! 2. restrict the increase to the decrease above
        ztmp2(1:kproma) = MERGE(ztmp2(1:kproma), MIN(ztmp2(1:kproma), &
                (zxi_bf_sed(1:kproma) / MAX(paclc(1:kproma), clc_min)) &
                - pxip1_diff_up(1:kproma)), &
                pxip1_diff_up(1:kproma) > 0)
        ! 3. convert back to gridmean
        pxip1(1:kproma) = ztmp2(1:kproma) * MAX(paclc(1:kproma), clc_min)
  ENDIF
  !<<UP #797
  !>>UP #798: change to epsilon where ll1 false
  IF (lsedfix) THEN
          picnc_upd(1:kproma) = MERGE(ztmp1(1:kproma), eps, ll1(1:kproma))
  ELSE
          ! original line
          picnc_upd(1:kproma) = MERGE(ztmp1(1:kproma), zicnc_gridmean(1:kproma), ll1(1:kproma))
  ENDIF
  !SF Note: ToDo: is it intended that picnc is set to a non-minimum value when paclc <= clc_min?
  !               I diagnosed the values of zicnc_gridmean(1:kproma) when paclc <= clc_min,
  !               and the result shows some high spots (up to ~ 1.e5)
  !UP note: Without cloud cover there should not be any ice, no? -> #798
  !<<UP

!UP note: below mass and number are used to calculate the fluxes
  zxi_delta(1:kproma) = zxi_bf_sed(1:kproma) - pxip1(1:kproma)

!UP note: this looks like some sort of level conversion
! The mmr is getting converted into a mass flux [kg/m^2/s]
  zxiflx_from_level(1:kproma) = zcons2*zxi_delta(1:kproma)*pdp(1:kproma)

!---for in-cloud scavenging
  !>>SF #453 convert to in-cloud value when cloud cover is non zero
  ztmp1(1:kproma)    = zxi_delta(1:kproma) / MAX(paclc(1:kproma), clc_min)
!UP note: sedimenting ice mass gets put to snow; but that would mean that there
!is no ice reaching a lower grid box?
  !>>UP #798: change to epsilon where ll1 false
  IF (lsedfix) THEN
          pmrateps(1:kproma) = MERGE(ztmp1(1:kproma), eps, ll1(1:kproma))
  ELSE
          ! original line
          pmrateps(1:kproma) = MERGE(ztmp1(1:kproma), zxi_delta(1:kproma), ll1(1:kproma))
  ENDIF
  !SF Note: it is not fully clear to me what should be the value in the absence of cloud cover (2nd argument of MERGE).
  !         Based on a discussion with David, and in consistency with what is done to icnc for the same pb (a few
  !         lines above), David and myself decided to implement it as above, ie leaving the grid-mean value where
  !         the cloud cover is below minimum.
  !         --> to be checked with Ulrike
  !<<UP #798
  !<<SF #453

  !SFNote #459: introducing a new variable for the grid box fraction covered by falling ice
  pclcfi(:) = gridbox_frac_falling_hydrometeor(kbdim, kproma, &
                        pxiflux(:), pclcfi(:), &
                        zxiflx_from_level(:), paclc(:))

!UP note: here the flux is getting updated
!UP: this should be using the corrected difference and the difference to
!the un-corrected one needs to be properly sublimated
  pxiflux(1:kproma)  = pxiflux(1:kproma)  + zxiflx_from_level(1:kproma)
!>>UP #797
! We attribute the difference due to the correction to instantaneous sublimation
! The update-routine will take care of the associated heat release and humidity
! changes
  IF (lsedfix) THEN
        ! copied and adapted from the sublimation routine
          ztmp4(1:kproma) = xsec*(pqsi(1:kproma)-pqm1(1:kproma))
          ztmp4(1:kproma) = MAX(ztmp4(1:kproma),0._dp)

          ztmp1(1:kproma) = &
                ! un-corrected zxiflx_from_level, as a mass difference
                (zxi_bf_sed(1:kproma) - ztmp3(1:kproma)) &
                ! substract the corrected zxiflx_from_level, convert to simple
                ! mass difference (same unit as pxisub)
                - zxiflx_from_level(1:kproma)/(zcons2*pdp(1:kproma))

          ! don't sublimate more (overall) than humidity allows
          ! MAX as safeguard
          pxisub(1:kproma) = pxisub(1:kproma) + MAX(ztmp1(1:kproma), 0._dp)
          !UP TODO: what if this condition strikes and what we've already
          !subtracted above cannot be sublimated?
          pxisub(1:kproma)  = MIN(pxisub(1:kproma),ztmp4(1:kproma))
          pxisub(1:kproma)  = MAX(pxisub(1:kproma),0._dp)
  ENDIF
!<<UP #797
! For the IC number flux the old zicnc_gridmean (not corrected for in-cloud
! restriction) is used and the old increase is subtracted from the flux -> the
! difference between the old and new increase is implicitly sublimated
  pxifluxn(1:kproma) = pxifluxn(1:kproma) + zcons2*(zicnc_gridmean_bf_sed(1:kproma) - &
                       zicnc_gridmean(1:kproma))*pdp(1:kproma)*prho_rcp(1:kproma)

!UP note: now this correction only captures the fluxes, not the icnc changes
!above; pxifluxn is getting corrected here, but it's actually not used after;
!pxiflux is only used to be converted to precip at the lowest level
  pxifluxn(:) = consistency_number_to_mass(kbdim, kproma, epsec, pxiflux(:), pxifluxn(:))

END SUBROUTINE sedimentation_ice

SUBROUTINE mixed_phase_deposition_and_corrections( &
              !-- IN
              kbdim, kproma, &
              papp1, picnc, pqm1, paclc, pesi, pesw, peta, pgenti, &
              plsdcp, plvdcp, pqp1, pqsm1, prho, ptp1, pxievap, & 
              pxip1, pxite, pvervx, & 
              !-- INOUT
              pcnd, pdep, &
              !-- OUT
              pqp1tmp, pqsp1tmp, ptp1tmp)

  INTEGER, INTENT(in) :: kbdim, kproma

  REAL(dp), INTENT(in) :: papp1(kbdim)   !< pressure at full levels (t-1)
  REAL(dp), INTENT(in) :: picnc(kbdim)   !< Ice crystal number concentration (ICNC) [1/m3]
  REAL(dp), INTENT(in) :: pqm1(kbdim)    !< specific humidity               (t-1)
  REAL(dp), INTENT(in) :: paclc(kbdim)   !< cloud cover
  REAL(dp), INTENT(in) :: pesi(kbdim)    !< Saturation vapor pressure w.r.t. ice [Pa]
  REAL(dp), INTENT(in) :: pesw(kbdim)    !< Saturation vapor pressure w.r.t. water [Pa]
  REAL(dp), INTENT(in) :: peta(kbdim)    !< Variable needed for the Bergeron-Findeisen process
  REAL(dp), INTENT(in) :: pgenti(kbdim)  !< Variable related to Tompkins cloud cover scheme
  REAL(dp), INTENT(in) :: plsdcp(kbdim)  !< latent heat of sublimation div. by the specific heat at constant pres.
  REAL(dp), INTENT(in) :: plvdcp(kbdim)  !< latent heat of vaporization div. by the specific heat at constant pres.
  REAL(dp), INTENT(in) :: pqp1(kbdim)    !< Specific humidity [kg/kg] (t)
  REAL(dp), INTENT(in) :: pqsm1(kbdim)   !< Sat. specific humidity (t-1) over ice/water dep. on temp. [kg/kg]
  REAL(dp), INTENT(in) :: prho(kbdim)    !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: ptp1(kbdim)    !< Temperature [K] (t)
  REAL(dp), INTENT(in) :: pxievap(kbdim) !< evaporation of cloud ice [kg/kg]
  REAL(dp), INTENT(in) :: pxip1(kbdim)   !< Ice mass sedimentation (grid-mean values)
  REAL(dp), INTENT(in) :: pxite(kbdim)  !< tendency of cloud ice from detrainment [kg/kg/s]
  REAL(dp), INTENT(in) :: pvervx(kbdim)  !< Updraft velocity [cm/s]

  REAL(dp), INTENT(inout) :: pcnd(kbdim) !< condensation rate [kg/kg]
  REAL(dp), INTENT(inout) :: pdep(kbdim) !< deposition rate [kg/kg]

  REAL(dp), INTENT(out) :: pqp1tmp(kbdim)  !< temporary value of the updated specific humidity (t) [kg/kg]
  REAL(dp), INTENT(out) :: pqsp1tmp(kbdim) !< Temporary new sat. specific humid. o. ice/water dep. on temp. [kg/kg]
  REAL(dp), INTENT(out) :: ptp1tmp(kbdim)  !< temporary value of the updated temperature (t) [K]

  !local vars:
  INTEGER :: jl
  INTEGER :: it(kbdim) !< index for temperature lookup table

  REAL(dp) :: zcor(kbdim)        !< Temporary variable needed for the calculation of the sat. specific humidity
  REAL(dp) :: zcorw(kbdim)       !< Temporary variable needed for the calculation of the sat. specific humidity
  REAL(dp) :: zlucua(kbdim)      !< Temporary variable needed for the calculation of the sat. vapor pressure (t-1)
  REAL(dp) :: zlucuap1(kbdim)    !< Temporary variable needed for the calculation of the sat. vapor pressure (t)
  REAL(dp) :: zlucuaw(kbdim)     !< Temporary variable needed for the calculation of the sat. vapor pressure (t-1)
  REAL(dp) :: zlucuawp1(kbdim)   !< Temporary variable needed for the calculation of the sat. vapor pressure (t)
  REAL(dp) :: zlucub(kbdim)      !< ???
  REAL(dp) :: zvervmax(kbdim)    !< Threshold vertical velocity needed for the Bergeron-Findeisen process
  REAL(dp) :: zdqsdt(kbdim)      !< Difference in sat. specific humidity between t and t-1
  REAL(dp) :: zes(kbdim)         !< Saturation vapor pressure w.r.t. water or ice
  REAL(dp) :: zesw(kbdim)        !< Saturation vapor pressure w.r.t. water
  REAL(dp) :: zlc(kbdim)         !< Latent heat of vaporation/sublimation over specific heat depending on temp.
  REAL(dp) :: zlcdqsdt(kbdim)    !< ???
  REAL(dp) :: zoversat(kbdim)    !< Security supersaturation over ice/water depending on temp. [kg/kg]
  REAL(dp) :: zoversatw(kbdim)   !< Security supersaturation w.r.t. water
  REAL(dp) :: zqcon(kbdim)       !< Temporary variable needed for the final calc. of condensation/deposition
  REAL(dp) :: zqsp1tmphet(kbdim) !< Onset critical specific humidity for heterogeneously formed cirrus clouds
  REAL(dp) :: zqsp1tmpw(kbdim)   !< Temporary new sat. specific humidity w.r.t. water [kg/kg]
  REAL(dp) :: zqst1(kbdim)       !< Saturation specific humidity (t) over ice/water depending on temp. [kg/kg]
  REAL(dp) :: zrhtest(kbdim)     !< Temporary variable needed for the final calc. of condensation/deposition
  REAL(dp) :: zrice(kbdim)       !< Volume mean ice crystal radius needed for the Bergeron-Findeisen process
  REAL(dp) :: zxip1(kbdim)       !< Cloud ice mass mixing ratio (t) [kg/kg]

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim), ztmp4(kbdim), ztmp5(kbdim)
  LOGICAL  :: ll1(kbdim), ll2(kbdim), ll3(kbdim), ll4(kbdim), ll5(kbdim), ll6(kbdim), lo2(kbdim)
  ! UP emulator phasing nic_cirrus
  REAL(dp) :: ztmp4_1(kbdim), ztmp4_2(kbdim)

  ptp1tmp(1:kproma) = ptp1(1:kproma)+plvdcp(1:kproma)*pcnd(1:kproma) &
                    + plsdcp(1:kproma)*pdep(1:kproma)
  pqp1tmp(1:kproma) = pqp1(1:kproma)-pcnd(1:kproma)-pdep(1:kproma)

  zxip1(1:kproma) = pxip1(1:kproma) + ztmst*pxite(1:kproma) &
                  - pxievap(1:kproma) + pgenti(1:kproma) + pdep(1:kproma)
  zxip1(1:kproma) = MAX(zxip1(1:kproma), 0._dp)

!>>SF #176 now this is done consistenly between cl. micro and radiation, with a general formula
  !SF conversion of ice mmr from grid-mean kg/kg to in-cloud g.m-3:
  ztmp1(1:kproma) = 1000._dp*zxip1(1:kproma)*prho(1:kproma)/MAX(paclc(1:kproma),clc_min)

  ztmp1(:) = eff_ice_crystal_radius(kbdim, kproma, ztmp1(:), picnc(:))
!<<SF
  ztmp1(1:kproma) = MIN(MAX(ztmp1(1:kproma), ceffmin), ceffmax) !SF zrieff in micrometers

  ! Simple param of r/re approximating the Schumann et al. 2011 data
  zrice(:) = effective_2_volmean_radius_param_Schuman_2011( &
                       kbdim, kproma, ztmp1(:))

  zvervmax(:) = threshold_vert_vel(kbdim, kproma, pesw(:), &
                          pesi(:), picnc(:), &
                          zrice(:), peta(:))

  lo2(1:kproma)   = (ptp1tmp(1:kproma) < cthomi) &
               .OR. (ptp1tmp(1:kproma) < tmelt .AND. 0.01_dp*pvervx(1:kproma) < zvervmax(1:kproma))

  CALL set_lookup_index( &
          !-- IN
          kbdim, kproma, ptp1tmp(:), &
          '2-m cloud micro (mixed_phase_deposition_and_corrections)', &
          !-- OUT
          it(:))

  DO jl=1,kproma
     zlucua(jl)  = tlucua(it(jl))
     zlucuaw(jl) = tlucuaw(it(jl))

     zlucuap1(jl)  = tlucua(it(jl)+1)
     zlucuawp1(jl) = tlucuaw(it(jl)+1)

     zlucub(jl) = tlucub(it(jl))
  ENDDO

  zlucua(1:kproma) = MERGE(zlucua(1:kproma), zlucuaw(1:kproma), lo2(1:kproma))

  CALL sat_spec_hum( &
          !-- IN
          kbdim, kproma, papp1(:), zlucua(:), &
          !-- OUT
          zes(:), zcor(:), pqsp1tmp(:))

  ll1(1:kproma) = (zes(1:kproma) < 0.4_dp)  !SF LO

  zoversat(1:kproma) = pqsp1tmp(1:kproma)*0.01_dp
 
  zrhtest(1:kproma) = pqm1(1:kproma)/pqsm1(1:kproma)
  zrhtest(1:kproma) = MIN(zrhtest(1:kproma), 1._dp)
  zrhtest(1:kproma) = zrhtest(1:kproma)*pqsp1tmp(1:kproma)

  CALL sat_spec_hum( &
          !-- IN
          kbdim, kproma, papp1(:), zlucuaw(:), &
          !-- OUT
          zesw(:), zcorw(:), zqsp1tmpw(:))

  zoversatw(1:kproma) = 0.01_dp*zqsp1tmpw(1:kproma)

  zqst1(1:kproma) = MERGE(zlucuap1(1:kproma),zlucuawp1(1:kproma),lo2(1:kproma))
  zqst1(1:kproma) = zqst1(1:kproma)/papp1(1:kproma)
  zqst1(1:kproma) = MIN(zqst1(1:kproma),0.5_dp)
  zqst1(1:kproma) = zqst1(1:kproma)/(1._dp-vtmpc1*zqst1(1:kproma))

  zdqsdt(1:kproma) = 1000._dp*(zqst1(1:kproma)-pqsp1tmp(1:kproma))

  zlc(1:kproma) = MERGE(plsdcp(1:kproma),plvdcp(1:kproma),lo2(1:kproma))

  ztmp1(1:kproma)    = zlc(1:kproma)*zdqsdt(1:kproma)
  ztmp2(1:kproma)    = pqsp1tmp(1:kproma)*zcor(1:kproma)*zlucub(1:kproma)
  zlcdqsdt(1:kproma) = MERGE(ztmp1(1:kproma), ztmp2(1:kproma), ll1(1:kproma))

  zqcon(1:kproma) = 1._dp/(1._dp+zlcdqsdt(1:kproma))

  ztmp1(1:kproma)       = zqsp1tmpw(1:kproma)+zoversatw(1:kproma)
  ztmp1(1:kproma)       = MIN(ztmp1(1:kproma),pqsp1tmp(1:kproma)*1.3_dp)
  zqsp1tmphet(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, lo2(1:kproma))

  !>>UP emulator phasing nic_cirrus
  IF (lemuphase_nic_cirrus) THEN
          ! case nic_cirrus = 1
                  ll1(1:kproma) = .true. 

                  ll2(1:kproma) = (pqp1tmp(1:kproma) > (pqsp1tmp(1:kproma)  + zoversat(1:kproma) ) )
                  ll3(1:kproma) = (pqp1tmp(1:kproma) > (zqsp1tmpw(1:kproma) + zoversatw(1:kproma)) )
                  ll4(1:kproma) = (pqp1tmp(1:kproma) > zqsp1tmphet(1:kproma))
                  ll5(1:kproma) = (ptp1tmp(1:kproma) >= cthomi)

                  ztmp1(1:kproma) = (pqp1tmp(1:kproma) - pqsp1tmp(1:kproma)  - zoversat(1:kproma) )*zqcon(1:kproma)
                  ztmp2(1:kproma) = (pqp1tmp(1:kproma) - zqsp1tmpw(1:kproma) - zoversatw(1:kproma))*zqcon(1:kproma)
                  ztmp3(1:kproma) = (pqp1tmp(1:kproma) - zqsp1tmphet(1:kproma)                )*zqcon(1:kproma)

                  !SF ice cloud cases:         
                  ll6(1:kproma) = (lo2(1:kproma) .AND. ll1(1:kproma) .AND. ll2(1:kproma)) &
                              .OR. &
                                  (lo2(1:kproma) .AND. (.NOT. ll1(1:kproma)) .AND. ll2(1:kproma) .AND. ll5(1:kproma))

                  ztmp4(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll6(1:kproma))

                  ll6(1:kproma) = lo2(1:kproma) &
                            .AND. (.NOT. ll1(1:kproma)) &
                            .AND. ll3(1:kproma) .AND. (.NOT. ll5(1:kproma)) .AND. (.NOT. ll_het)

                  ztmp4(1:kproma) = MERGE(ztmp2(1:kproma), ztmp4(1:kproma), ll6(1:kproma))

                  ll6(1:kproma) = lo2(1:kproma) &
                            .AND. (.NOT. ll1(1:kproma)) &
                            .AND. ll4(1:kproma) .AND. (.NOT. ll5(1:kproma)) .AND. ll_het

                  ztmp4_1(1:kproma) = MERGE(ztmp3(1:kproma), ztmp4(1:kproma), ll6(1:kproma))

          ! case nic_cirrus = 2
                  ll1(1:kproma) = .false. 

                  ll2(1:kproma) = (pqp1tmp(1:kproma) > (pqsp1tmp(1:kproma)  + zoversat(1:kproma) ) )
                  ll3(1:kproma) = (pqp1tmp(1:kproma) > (zqsp1tmpw(1:kproma) + zoversatw(1:kproma)) )
                  ll4(1:kproma) = (pqp1tmp(1:kproma) > zqsp1tmphet(1:kproma))
                  ll5(1:kproma) = (ptp1tmp(1:kproma) >= cthomi)

                  ztmp1(1:kproma) = (pqp1tmp(1:kproma) - pqsp1tmp(1:kproma)  - zoversat(1:kproma) )*zqcon(1:kproma)
                  ztmp2(1:kproma) = (pqp1tmp(1:kproma) - zqsp1tmpw(1:kproma) - zoversatw(1:kproma))*zqcon(1:kproma)
                  ztmp3(1:kproma) = (pqp1tmp(1:kproma) - zqsp1tmphet(1:kproma)                )*zqcon(1:kproma)

                  !SF ice cloud cases:         
                  ll6(1:kproma) = (lo2(1:kproma) .AND. ll1(1:kproma) .AND. ll2(1:kproma)) &
                              .OR. &
                                  (lo2(1:kproma) .AND. (.NOT. ll1(1:kproma)) .AND. ll2(1:kproma) .AND. ll5(1:kproma))

                  ztmp4(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll6(1:kproma))

                  ll6(1:kproma) = lo2(1:kproma) &
                            .AND. (.NOT. ll1(1:kproma)) &
                            .AND. ll3(1:kproma) .AND. (.NOT. ll5(1:kproma)) .AND. (.NOT. ll_het)

                  ztmp4(1:kproma) = MERGE(ztmp2(1:kproma), ztmp4(1:kproma), ll6(1:kproma))

                  ll6(1:kproma) = lo2(1:kproma) &
                            .AND. (.NOT. ll1(1:kproma)) &
                            .AND. ll4(1:kproma) .AND. (.NOT. ll5(1:kproma)) .AND. ll_het

                  ztmp4_2(1:kproma) = MERGE(ztmp3(1:kproma), ztmp4(1:kproma), ll6(1:kproma))

          ! combine
                  pdep(1:kproma) = pdep(1:kproma) + eta_emu_nic_cirrus*ztmp4_1(1:kproma) + (1-eta_emu_nic_cirrus) * &
                  ztmp4_2(1:kproma)

  ELSE
  !<<UP

  ll1(1:kproma) = (nic_cirrus == 1)

  ll2(1:kproma) = (pqp1tmp(1:kproma) > (pqsp1tmp(1:kproma)  + zoversat(1:kproma) ) )
  ll3(1:kproma) = (pqp1tmp(1:kproma) > (zqsp1tmpw(1:kproma) + zoversatw(1:kproma)) )
  ll4(1:kproma) = (pqp1tmp(1:kproma) > zqsp1tmphet(1:kproma))
  ll5(1:kproma) = (ptp1tmp(1:kproma) >= cthomi)

  ztmp1(1:kproma) = (pqp1tmp(1:kproma) - pqsp1tmp(1:kproma)  - zoversat(1:kproma) )*zqcon(1:kproma)
  ztmp2(1:kproma) = (pqp1tmp(1:kproma) - zqsp1tmpw(1:kproma) - zoversatw(1:kproma))*zqcon(1:kproma)
  ztmp3(1:kproma) = (pqp1tmp(1:kproma) - zqsp1tmphet(1:kproma)                )*zqcon(1:kproma)

  !SF ice cloud cases:         
  ll6(1:kproma) = (lo2(1:kproma) .AND. ll1(1:kproma) .AND. ll2(1:kproma)) &
              .OR. &
                  (lo2(1:kproma) .AND. (.NOT. ll1(1:kproma)) .AND. ll2(1:kproma) .AND. ll5(1:kproma))
  !UP comment: difference is that first has nic_cirrus=1 and second has nic_cirrus!=1 with T > homogeneous freezing

  !UP comment: above ll6 is apllied here to update ztmp4
  ztmp4(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll6(1:kproma))

  ll6(1:kproma) = lo2(1:kproma) &
            .AND. (.NOT. ll1(1:kproma)) &
            .AND. ll3(1:kproma) .AND. (.NOT. ll5(1:kproma)) .AND. (.NOT. ll_het)

  !UP comment: new update of ztmp4 only if nic_cirrus != 1
  ztmp4(1:kproma) = MERGE(ztmp2(1:kproma), ztmp4(1:kproma), ll6(1:kproma))

  ll6(1:kproma) = lo2(1:kproma) &
            .AND. (.NOT. ll1(1:kproma)) &
            .AND. ll4(1:kproma) .AND. (.NOT. ll5(1:kproma)) .AND. ll_het

  !UP comment: new update of ztmp4 only if nic_cirrus != 1
  ztmp4(1:kproma) = MERGE(ztmp3(1:kproma), ztmp4(1:kproma), ll6(1:kproma))

  pdep(1:kproma) = pdep(1:kproma) + ztmp4(1:kproma)

  ENDIF !UP lemuphase_nic_cirrus
  !------------------------------- UP until here -----------------------

  !SF water cloud cases:
  ll6(1:kproma) = (.NOT. lo2(1:kproma)) .AND. ll2(1:kproma)

  ztmp4(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll6(1:kproma))
 
  pcnd(1:kproma) = pcnd(1:kproma) + ztmp4(1:kproma)

  !SF final corrections for pdep and pcnd:
  ztmp5(1:kproma) = pqp1(1:kproma)-zrhtest(1:kproma)
  ztmp5(1:kproma) = MAX(ztmp5(1:kproma), 0._dp)

  ll1(1:kproma) = (pdep(1:kproma)     > 0._dp            )
  ll2(1:kproma) = (pcnd(1:kproma)     > 0._dp            )
  ll3(1:kproma) = (pqp1tmp(1:kproma)  < zrhtest(1:kproma))
  ll4(1:kproma) = (pqsp1tmp(1:kproma) <= pqsm1(1:kproma) )

  ll5(1:kproma) = lo2(1:kproma) .AND. ll1(1:kproma) .AND. ll3(1:kproma) .AND. ll4(1:kproma)

  !UP comment: new update of pdep only if nic_cirrus == 1 - no, this is not related to nic_cirrus directly anymore
  pdep(1:kproma) = MERGE(ztmp5(1:kproma), pdep(1:kproma), ll5(1:kproma)) 

  ll5(1:kproma) = (.NOT. lo2(1:kproma)) .AND. ll2(1:kproma) .AND. ll3(1:kproma) .AND. ll4(1:kproma)

  pcnd(1:kproma) = MERGE(ztmp5(1:kproma), pcnd(1:kproma), ll5(1:kproma))

  ptp1tmp(1:kproma) = ptp1(1:kproma)+plvdcp(1:kproma)*pcnd(1:kproma) &
                    + plsdcp(1:kproma)*pdep(1:kproma)
  pqp1tmp(1:kproma) = pqp1(1:kproma)-pcnd(1:kproma)-pdep(1:kproma) !SF #436 (update pqp1tmp)
  
END SUBROUTINE mixed_phase_deposition_and_corrections

SUBROUTINE update_in_cloud_water( &
              !-- IN
              kbdim, kproma, &
              pap, pcdncact, pcnd, pdep, pgenti, pgentl, &
              pnicex, pqp1tmp, pqsp1tmp, prho, prid, ptm1, &
              pcdnc, picnc, &
              !-- INOUT
              ld_cc, picnc_upd, pqnuc, pcdnc_upd, paclc, pxib, pxlb, &
              pninucl, & ! UP #783.3
              !-- OUT
              pcdnc_min, ptmp1, ptmp2, & !DN #475 (min cdnc) !davidn
              pcorric3) ! UP #783.3

!DN: update of in cloud water/ice, CDNC/ICNC(activation/nucleation), cloud cover

  INTEGER, INTENT(in) :: kbdim, kproma

  REAL(dp), INTENT(in) :: pap(kbdim)      !< total number of aerosols available
  REAL(dp), INTENT(in) :: pcdncact(kbdim) !< Number of newly activated cloud droplets
  REAL(dp), INTENT(in) :: pcnd(kbdim)     !< condensation rate [kg/kg]
  REAL(dp), INTENT(in) :: pdep(kbdim)     !< deposition rate [kg/kg]
  REAL(dp), INTENT(in) :: pgenti(kbdim)   !< Variable related to Tompkins cloud cover scheme
  REAL(dp), INTENT(in) :: pgentl(kbdim)   !< Variable related to Tompkins cloud cover scheme
  REAL(dp), INTENT(in) :: pnicex(kbdim)   !< number of newly formed ice crystals
  REAL(dp), INTENT(in) :: pqp1tmp(kbdim)  !< temporary value of the updated specific humidity (t) [kg/kg]
  REAL(dp), INTENT(in) :: pqsp1tmp(kbdim) !< Temporary new sat. specific humid. o. ice/water dep. on temp. [kg/kg]
  REAL(dp), INTENT(in) :: prho(kbdim)     !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: prid(kbdim)     !< Mean volume radius of ice crystals [m]
  REAL(dp), INTENT(in) :: ptm1(kbdim)     !< temperature (t-1)
  REAL(dp), INTENT(in) :: pcdnc(kbdim)    !< Cloud droplet number concentration (CDNC) [1/m3] 
  REAL(dp), INTENT(in) :: picnc(kbdim)    !< Ice crystal number concentration (ICNC) [1/m3]

  LOGICAL,  INTENT(inout) :: ld_cc(kbdim) !< True when cloud cover in a given layer is greater than eps
  REAL(dp), INTENT(inout) :: picnc_upd(kbdim) !< Ice crystal number concentration (ICNC) [1/m3]
  !SF Note: intent(inout) important for picnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pqnuc(kbdim) !< cloud droplet nucleation rate [m-3 s-1]
  REAL(dp), INTENT(inout) :: pcdnc_upd(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3] 
  !SF Note: intent(inout) important for pcdnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: paclc(kbdim) !< Cloud cover
  REAL(dp), INTENT(inout) :: pxib(kbdim)  !< Cloud ice in the cloudy part of the grid box [kg/kg] 
  REAL(dp), INTENT(inout) :: pxlb(kbdim)  !< Cloud liquid water in the cloudy part of the grid box [kg/kg] 
  REAL(dp), INTENT(inout) :: pninucl(kbdim) ! UP #783.3; number conc. of newly
                                            ! nucleated IC [1/m3], gets updated
                                            ! here, cirrus

!>>DN #475 (min cdnc)
  REAL(dp), INTENT(out) :: pcdnc_min(kbdim) !< Minimum CDNC concentration computed from maximum radius [1/m3]
!>>DN #475 (min cdnc)
  REAL(dp), INTENT(out) :: ptmp1(kbdim),ptmp2(kbdim)!davidn
  REAL(dp), INTENT(out) :: pcorric3(kbdim) ! UP #783.3; correction term for min.
                                           ! icemin

  !local vars:
  !UP #787: remove below line because it is useless
  !REAL(dp) :: zqlnuc(kbdim)  !< Nucleation rate of CDNC [1/m3]
  REAL(dp) :: zrelhum(kbdim) !< Relative humidity

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim), ztmp4(kbdim)
  LOGICAL  :: ll1(kbdim), ll2(kbdim)
  !>>UP emulator phasing nic_cirrus
  REAL(dp) :: ztmp1_1(kbdim), ztmp1_2(kbdim)
  !<<UP

  !UP update ice nucleation timer to include the following
  !>>UP timers
  IF (ltimer) THEN
     CALL timer_start(timer_cmp_icenucl)
  END IF
  !<<UP

  zrelhum(1:kproma) = pqp1tmp(1:kproma)/pqsp1tmp(1:kproma)
  !SFnote: security for division?

  ztmp1(1:kproma) = pdep(1:kproma) + pgenti(1:kproma)  
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0._dp) !SF deposition zdepos
 
  ztmp2(1:kproma) = pcnd(1:kproma) + pgentl(1:kproma)
  ztmp2(1:kproma) = MAX(ztmp2(1:kproma), 0._dp) !SF condensation zcond 

  ztmp3(1:kproma) = pxib(1:kproma)+pdep(1:kproma)/MAX(paclc(1:kproma), clc_min) !DN / SF #473: (lower clc limit) 
  ztmp3(1:kproma) = MAX(ztmp3(1:kproma), 0._dp)

  ztmp4(1:kproma) = pxlb(1:kproma)+pcnd(1:kproma)/MAX(paclc(1:kproma), clc_min) !DN / SF #473: (lower clc limit)
  ztmp4(1:kproma) = MAX(ztmp4(1:kproma), 0._dp)

  pxib(1:kproma) = MERGE(ztmp3(1:kproma), pxib(1:kproma), ld_cc(1:kproma))
  pxlb(1:kproma) = MERGE(ztmp4(1:kproma), pxlb(1:kproma), ld_cc(1:kproma))

  ll1(1:kproma) = (.NOT. ld_cc(1:kproma)) &
            .AND. ( (ztmp1(1:kproma) > 0._dp) .OR. (ztmp2(1:kproma) > 0._dp) )        

  ztmp3(1:kproma) = MAX(MIN(zrelhum(1:kproma), 1.0_dp), 0.01_dp)
  paclc(1:kproma) = MERGE(ztmp3(1:kproma), paclc(1:kproma), ll1(1:kproma))

  ztmp3(1:kproma) = ztmp1(1:kproma) / MAX(paclc(1:kproma), clc_min) !DN / SF #473: (lower clc limit)
  ztmp4(1:kproma) = ztmp2(1:kproma) / MAX(paclc(1:kproma), clc_min) !DN / SF #473: (lower clc limit)

  pxib(1:kproma) = MERGE(ztmp3(1:kproma), pxib(1:kproma), ll1(1:kproma))
  pxlb(1:kproma) = MERGE(ztmp4(1:kproma), pxlb(1:kproma), ll1(1:kproma))

  !UP end update ice nucleation timer
  !>>UP timers
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_icenucl)
  END IF
  !<<UP

!>>DN #475 (min cdnc)
  ztmp1(1:kproma) = pxlb(1:kproma)*prho(1:kproma)
  pcdnc_min(:)      = minimum_CDNC(kbdim, kproma, ztmp1(:))
!<<DN #475 (min cdnc)

  !SF re-define ld_cc which was no longer true:
  !SF ToDo WARNING : this definition is inconsistent with its previous def
  !                  shouldn't the lower limit set tp clc_min also?
  ld_cc(1:kproma) = (paclc(1:kproma) > 0._dp)

  !SF update cdnc:
  ll1(1:kproma) = ld_cc(1:kproma) .AND. (pxlb(1:kproma) > cqtmin)
  ll2(1:kproma) = ll1(1:kproma) &
            .AND. (pcdnc(1:kproma) <= pcdnc_min(1:kproma)) & !if there was no previous nucleation
            .AND. (ptm1(1:kproma) > cthomi)     !SF added a temperature condition for consistency

  !>>UP #805
  IF (lemuphase_cdnuc) THEN
          ztmp1(1:kproma)  = pcdncact(1:kproma) - pcdnc(1:kproma)
  ELSE
  !UP original
          ztmp1(1:kproma)  = pcdncact(1:kproma) - pcdnc(1:kproma)
  ENDIF
  !<<UP #805
  ztmp1(1:kproma)  = MAX(0._dp, ztmp1(1:kproma))
  !UP #787: remove below line because it is useless
  !zqlnuc(1:kproma) = MERGE(ztmp1(1:kproma), zqlnuc(1:kproma), ll2(1:kproma))

  ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma))
  ! UP: merging DN's r3593, he has the following line here that did not get merged 
  ! pcdnc(1:kproma) = pcdnc(1:kproma) + ztmp1(1:kproma)
  pqnuc(1:kproma) = pqnuc(1:kproma) + zdt*ztmp1(1:kproma)

  ptmp1(1:kproma) = ztmp1(1:kproma)!davidn
  ptmp2(1:kproma) = pcdnc(1:kproma)!davidn 

!UP: nucleation
  ztmp2(1:kproma) = pcdnc(1:kproma) + ztmp1(1:kproma)
!UP: cdnc min
  ztmp2(1:kproma) = MAX(ztmp2(1:kproma), pcdnc_min(1:kproma))
!UP: evaporation where there is no cloud water
!UP: with the check below, the correction for pcdnc_min is also only applied
!inside clouds
  pcdnc_upd(1:kproma) = MERGE(ztmp2(1:kproma), cqtmin, ll1(1:kproma)) 

  !UP update ice nucleation timer to include the following
  !>>UP timers
  IF (ltimer) THEN
     CALL timer_start(timer_cmp_icenucl)
  END IF
  !<<UP

  !SF update icnc:
  ll1(1:kproma) = ld_cc(1:kproma) .AND. (pxib(1:kproma)  >  cqtmin) 
  ll2(1:kproma) = ll1(1:kproma)   .AND. (picnc(1:kproma) <= icemin)

!SF modified cond. statement:
  !>>UP emulator phasing nic_cirrus
  IF (lemuphase_nic_cirrus) THEN
          ! vars used here are mostly passed to the subroutine or computed the same independently of nic_cirrus
          ! so I can simply compute both and then combine
          ztmp1_1(1:kproma) = 0.75_dp / (pi*rhoice) * prho(1:kproma) * pxib(1:kproma) / prid(1:kproma)**3
          ztmp1_2(1:kproma) = MIN(pnicex(1:kproma), pap(1:kproma)*1.e6_dp)
          ztmp1(1:kproma) = eta_emu_nic_cirrus*ztmp1_1(1:kproma) + (1-eta_emu_nic_cirrus)*ztmp1_2(1:kproma)
  ELSE
  !<<UP
  IF (nic_cirrus == 1) THEN
     ztmp1(1:kproma) = 0.75_dp / (pi*rhoice) * prho(1:kproma) * pxib(1:kproma) / prid(1:kproma)**3
  ELSE IF (nic_cirrus == 2) THEN
     ztmp1(1:kproma) = MIN(pnicex(1:kproma), pap(1:kproma)*1.e6_dp)
  ENDIF
  ENDIF !UP lemuphase_nic_cirrus

  ! UP comment: nucleation
  ! UP #783.3: include update to nucleation rate here
  pninucl(1:kproma)   = pninucl(1:kproma) + MERGE(MAX(ztmp1(1:kproma)-picnc(1:kproma),0._dp), 0._dp,ll2(1:kproma))
  ztmp1(1:kproma)     = MERGE(ztmp1(1:kproma), picnc(1:kproma), ll2(1:kproma))
  ! UP comment: ice min
  pcorric3(1:kproma)  = MAX((icemin - ztmp1(1:kproma)),0._dp) !UP #783.3
  ztmp1(1:kproma)     = MAX(ztmp1(1:kproma), icemin)
  picnc_upd(1:kproma) = MERGE(ztmp1(1:kproma), cqtmin, ll1(1:kproma)) 
  
  !UP end ice nucleation timer update
  !>>UP timers
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_icenucl)
  END IF
  !<<UP

END SUBROUTINE update_in_cloud_water

SUBROUTINE freezing_below_238K( &
              !-- IN
              kbdim, kproma, &
              ld_frz_below_238K, paclc, pcdnc_min, & !DN #475 (min cdnc)
              pcdnc, picnc, &
              !-- INOUT
              picnc_upd, pqfre, pcdnc_upd, pfrl, pxib, pxlb, pfrzn, ptmp1)!davidn

  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL,INTENT(in) :: ld_frz_below_238K(kbdim) !< physical condition for freezing below 238K to occur

  REAL(dp), INTENT(in) :: paclc(kbdim)   !< Cloud cover 
!>>DN #475 (min cdnc)
  REAL(dp), INTENT(in) :: pcdnc_min(kbdim) !< Minimum CDNC concentration computed from maximum radius [1/m3] 
!<<DN #475 (min cdnc)
  REAL(dp), INTENT(in) :: pcdnc(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp), INTENT(in) :: picnc(kbdim) !< Ice crystal number concentration (ICNC) [1/m3]

  REAL(dp), INTENT(inout) :: picnc_upd(kbdim) !< Ice crystal number concentration (ICNC) [1/m3]
  !SF Note: intent(inout) important for picnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pqfre(kbdim) !< cloud droplet freezing rate [m-3 s-1]
  REAL(dp), INTENT(inout) :: pcdnc_upd(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]
  !SF Note: intent(inout) important for pcdnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pfrl(kbdim)  !< freezing rate [kg/kg]
  REAL(dp), INTENT(inout) :: pxib(kbdim)  !< Cloud ice in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(inout) :: pxlb(kbdim)  !< Cloud liquid water in the cloudy part of the grid box [kg/kg]

  REAL(dp), INTENT(out) :: pfrzn(kbdim),ptmp1(kbdim)  !davidn

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim)

  !>>UPcomment: freezing rate is liquid water per cloudy * cloudy, wherever the condition for freezing below 238K is fulfilled
  ztmp1(1:kproma) = pfrl(1:kproma) + pxlb(1:kproma)*paclc(1:kproma)
  pfrl(1:kproma)  = MERGE(ztmp1(1:kproma), pfrl(1:kproma), ld_frz_below_238K(1:kproma))
  !<<UPcomment

  !>>UPcomment: below all liquid gets converted to ice where the condition for freezing below 238K is fulfilled
  ztmp1(1:kproma) = pxib(1:kproma)+pxlb(1:kproma)
  pxib(1:kproma)  = MERGE(ztmp1(1:kproma), pxib(1:kproma), ld_frz_below_238K(1:kproma))

  pxlb(1:kproma)  = MERGE(0._dp, pxlb(1:kproma), ld_frz_below_238K(1:kproma))
  !<<UPcomment

!--- Included for prognostic CDNC/IC scheme ----------------------------
  ztmp1(1:kproma) = pcdnc(1:kproma)-pcdnc_min(1:kproma)
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0._dp)
!>>DN: included in zcorrcd10
!  ztmp2(1:kproma) = pqfre(1:kproma) - zdt*ztmp1(1:kproma)
!  pqfre(1:kproma) = MERGE(ztmp2(1:kproma), pqfre(1:kproma), ld_frz_below_238K(1:kproma))
!<<DN
!>>SF Kasja diags
  pfrzn(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp,  ld_frz_below_238K(1:kproma))
!<<SF Kasja diags
 
  ztmp2(1:kproma)     = picnc(1:kproma) + ztmp1(1:kproma)

  !UP: original line; uncommented below to remain on linear branch, TODO: take
  !back in when sure of this change
  ! see r6806, this requires retuning and is thus excluded for now
  picnc_upd(1:kproma) = MERGE(ztmp2(1:kproma), picnc(1:kproma), ld_frz_below_238K(1:kproma))
  !>> UP need to put a maximum condition for picnc here (don't create more than
  ! max), #783.3
  ! This counts as sth. like immediate self-collection
!  IF (lpr_corr) THEN
!          picnc_upd(1:kproma) = MERGE(MIN(ztmp2(1:kproma), icemax), picnc(1:kproma), ld_frz_below_238K(1:kproma))
!  ELSE
  ! original line
!          picnc_upd(1:kproma) = MERGE(ztmp2(1:kproma), picnc(1:kproma), ld_frz_below_238K(1:kproma))
!  ENDIF
  !<<UP

!>>DN: new diags
  ptmp1(1:kproma)  = MAX(pcdnc(1:kproma)-cqtmin,0._dp)
!<<DN: new diags        
  
  pcdnc_upd(1:kproma) = MERGE(cqtmin, pcdnc(1:kproma), ld_frz_below_238K(1:kproma)) 

END SUBROUTINE freezing_below_238K

SUBROUTINE het_mxphase_freezing( &
              !-- IN
              kbdim, kproma, &
              ld_mxphase_frz, papp1, ptkem1, pvervel, paclc, &
              pfracbcsol, pfracbcinsol, pfracdusol, pfracduai, pfracduci, &
              prho, prho_rcp, &
              prwetki, prwetai, prwetci, ptp1tmp, pcdnc_min, pcdnc, &
              picnc, &
              !-- INOUT
              picnc_upd, pcdnc_upd, &
              pfrl, pxib, pxlb, &
              !-- OUT
              pfrln )

  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL, INTENT(in)  :: ld_mxphase_frz(kbdim) !< Physical conditions for het. mixed phase freezing to occur

  REAL(dp), INTENT(in) :: papp1(kbdim)        !< pressure at full levels (t-1)
  REAL(dp), INTENT(in) :: ptkem1(kbdim)       !< turbulent kinetic energy (t-1)
  REAL(dp), INTENT(in) :: pvervel(kbdim)      !< large scale vertical velocity [m s-1]
  REAL(dp), INTENT(in) :: paclc(kbdim)      !< Cloud cover
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
  REAL(dp), INTENT(in) :: pcdnc_min(kbdim)      !< Minimum CDNC concentration computed from maximum radius [1/m3] 
  REAL(dp), INTENT(in) :: pcdnc(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp), INTENT(in) :: picnc(kbdim) !< Ice crystal number concentration (ICNC) [1/m3]

  REAL(dp), INTENT(inout) :: picnc_upd(kbdim) !< Ice crystal number concentration (ICNC) [1/m3]
  !SF Note: intent(inout) important for picnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pcdnc_upd(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]
  !SF Note: intent(inout) important for pcdnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pfrl(kbdim)  !< freezing rate [kg/kg]
  REAL(dp), INTENT(inout) :: pxib(kbdim)  !< Cloud ice in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(inout) :: pxlb(kbdim)  !< Cloud liquid water in the cloudy part of the grid box [kg/kg]

  REAL(dp), INTENT(out) :: pfrln(kbdim) !< Freezing rate for number conc. [1/m3]
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

  zdfarduci(1:kproma) = ak * ptp1tmp(1:kproma) * ztmp3(1:kproma) &
                      / ( 6._dp*pi*zetaair(1:kproma)*(prwetci(1:kproma)+eps) )
  zdfarduci(1:kproma) = MERGE(0._dp, zdfarduci(1:kproma), ll3(1:kproma))
 
  DO jl=1,kproma
 
     zradl = (0.75_dp*pxlb(jl) * prho(jl) &
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
             * ( pcdnc(jl)+picnc(jl) )

     zfrzcnt = pxlb(jl)*(1._dp-EXP(-zfrzcnt/MAX(pxlb(jl), cqtmin)*ztmst))

     !UP: I think these are the temperature dependencies of the individual
     !species (from Fig. 1 in Lohmann and Diehl (2006))
     znaimmdu = 32.3_dp*pfracdusol(jl)     ! montmorillonite 
     !znaimmdu  = 6.15E-2_dp*pfracdusol((jl))  ! kaolinite 

     znaimmbc = 2.91E-3_dp*pfracbcsol(jl)

     !>>SF ToDo: create a function for computing zomega (see vervx calc in cloud_micro_interface)
     zomega = pvervel(jl) - fact_tke*SQRT(ptkem1(jl))*prho(jl)*grav !SF #345 changed TKE prefactor
     !<<SF ToDo
     ztte   = zomega / cpd *prho_rcp(jl)

     ! Eq. 3 (Lohmann and Diehl (2006))
     zfrzimm = -(znaimmdu+znaimmbc)*prho(jl)/rhoh2o*EXP(tmelt-ptp1tmp(jl))*MIN(ztte,0._dp) 
     zfrzimm = pxlb(jl)*(1._dp-EXP(-zfrzimm*pxlb(jl)/pcdnc(jl)*ztmst))

     !>>MA ToDo Issue#790 Lemuphase_hetfrz 
     IF (lemuphase_hetfrz) THEN 
         ztmp1(jl) =  eta_emu_hetfrz * (zfrzcnt + zfrzimm)      
     ELSE
         ztmp1(jl) = zfrzcnt + zfrzimm                      !MA orig.
     ENDIF
     !<< MA ToDo

     ztmp1(jl) = MAX(0.0_dp,MIN(ztmp1(jl),pxlb(jl))) !SF pfrl surrogate

     ztmp2(jl) = pcdnc(jl)*ztmp1(jl)/(pxlb(jl)+eps)

     ztmp2(jl) = MAX(ztmp2(jl), 0._dp)  !SF pfrln surrogate
 
  ENDDO !SF end loop jl
 
  pfrl(1:kproma) = MERGE(ztmp1(1:kproma), pfrl(1:kproma) , ld_mxphase_frz(1:kproma))
 
  ztmp1(1:kproma) = MIN(ztmp2(1:kproma), pcdnc(1:kproma)-pcdnc_min(1:kproma))
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0._dp)
  pfrln(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_mxphase_frz(1:kproma))
 
  ztmp1(1:kproma) = pcdnc(1:kproma)-pfrln(1:kproma)
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), cqtmin)
  pcdnc_upd(1:kproma) = MERGE(ztmp1(1:kproma), pcdnc(1:kproma), ld_mxphase_frz(1:kproma))
 
  ztmp1(1:kproma) = picnc(1:kproma)+pfrln(1:kproma)
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), cqtmin)
  picnc_upd(1:kproma) = MERGE(ztmp1(1:kproma), picnc(1:kproma), ld_mxphase_frz(1:kproma))
 
  ztmp1(1:kproma) = pxlb(1:kproma)-pfrl(1:kproma)
  pxlb(1:kproma)  = MERGE(ztmp1(1:kproma), pxlb(1:kproma), ld_mxphase_frz(1:kproma))
 
  ztmp1(1:kproma) = pxib(1:kproma)+pfrl(1:kproma)
  pxib(1:kproma)  = MERGE(ztmp1(1:kproma), pxib(1:kproma), ld_mxphase_frz(1:kproma))
 
  ztmp1(1:kproma) = pfrl(1:kproma)*paclc(1:kproma)
  pfrl(1:kproma)  = MERGE(ztmp1(1:kproma), pfrl(1:kproma), ld_mxphase_frz(1:kproma))
 
END SUBROUTINE het_mxphase_freezing

SUBROUTINE WBF_process( &
              !-- IN
              kbdim, kproma, &
              ld_WBF, paclc, plsdcp, plvdcp, pcdnc, &
              !-- INOUT
              pcdnc_upd, pxlb, pxib, pxlte, pxite, ptte, pevabf, plwcbf)

  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL, INTENT(in) :: ld_WBF(kbdim) !< Physical conditions for WBF to occur

  REAL(dp), INTENT(in) :: paclc(kbdim)  !< Cloud cover
  REAL(dp), INTENT(in) :: plsdcp(kbdim) !< lat heat of sublimation div. by the specific heat at const pres.
  REAL(dp), INTENT(in) :: plvdcp(kbdim) !< lat heat of vaporization div. by the specific heat at const pres.
  REAL(dp), INTENT(in) :: pcdnc(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]

  REAL(dp), INTENT(inout) :: pcdnc_upd(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]
  !SF Note: intent(inout) important for pcdnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pxlb(kbdim)  !< Cloud liquid water in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(inout) :: pxib(kbdim)  !< Cloud ice in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(inout) :: pxlte(kbdim) !< tendency of cloud liquid water
  REAL(dp), INTENT(inout) :: pxite(kbdim) !< tendency of cloud ice
  REAL(dp), INTENT(inout) :: ptte(kbdim)  !< tendency of temperature

  REAL(dp), INTENT(out) :: pevabf(kbdim),plwcbf(kbdim) !davidn

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim)

  !>>UP #765
  IF (lemuphase_wbf) THEN
          ztmp1(1:kproma) = ztmst_rcp*pxlb(1:kproma)*paclc(1:kproma) !SF zzevp
               
          ztmp2(1:kproma) = pxlte(1:kproma)-eta_emu_wbf*ztmp1(1:kproma)
          pxlte(1:kproma) = MERGE(ztmp2(1:kproma), pxlte(1:kproma), ld_WBF(1:kproma))
         
          ztmp2(1:kproma) = pxite(1:kproma)+eta_emu_wbf*ztmp1(1:kproma)
          pxite(1:kproma) = MERGE(ztmp2(1:kproma), pxite(1:kproma), ld_WBF(1:kproma))
         
          ztmp2(1:kproma) = ptte(1:kproma)+(plsdcp(1:kproma)-plvdcp(1:kproma))*ztmp1(1:kproma)*eta_emu_wbf
          ptte(1:kproma)  = MERGE(ztmp2(1:kproma), ptte(1:kproma), ld_WBF(1:kproma))
           
        !>>DN: new diags
          ztmp1(1:kproma) = MAX(0._dp,eta_emu_wbf*(pcdnc(1:kproma)-cqtmin)) 
          pevabf(1:kproma)   = MERGE(ztmp1(1:kproma),0._dp, ld_WBF(1:kproma))
          plwcbf(1:kproma)   = MERGE(eta_emu_wbf*pxlb(1:kproma)*paclc(1:kproma),0._dp, ld_WBF(1:kproma))
        !<<DN: new diags

          pcdnc_upd(1:kproma) = MERGE((1 - eta_emu_wbf)*pcdnc(1:kproma)+eta_emu_wbf*cqtmin, &
          pcdnc(1:kproma), ld_WBF(1:kproma))
          !UP: formulate it like above, not like this:
          !pcdnc(1:kproma)-eta_emu_wbf*(pcdnc(1:kproma) - cqtmin)
          ! to not run into the numerical problem from #766
         
          ztmp2(1:kproma) = pxib(1:kproma)+eta_emu_wbf*pxlb(1:kproma)
          pxib(1:kproma)  = MERGE(ztmp2(1:kproma), pxib(1:kproma), ld_WBF(1:kproma))

          !>>UPcomment: here WBF eats up all water within the cloud
          pxlb(1:kproma) = MERGE((1 - eta_emu_wbf)*pxlb(1:kproma), pxlb(1:kproma), ld_WBF(1:kproma))
          !<<UPcomment
  ELSE
  !UP original
          ztmp1(1:kproma) = ztmst_rcp*pxlb(1:kproma)*paclc(1:kproma) !SF zzevp
               
          ztmp2(1:kproma) = pxlte(1:kproma)-ztmp1(1:kproma)
          pxlte(1:kproma) = MERGE(ztmp2(1:kproma), pxlte(1:kproma), ld_WBF(1:kproma))
         
          ztmp2(1:kproma) = pxite(1:kproma)+ztmp1(1:kproma)
          pxite(1:kproma) = MERGE(ztmp2(1:kproma), pxite(1:kproma), ld_WBF(1:kproma))
         
          ztmp2(1:kproma) = ptte(1:kproma)+(plsdcp(1:kproma)-plvdcp(1:kproma))*ztmp1(1:kproma)
          ptte(1:kproma)  = MERGE(ztmp2(1:kproma), ptte(1:kproma), ld_WBF(1:kproma))
           
        !>>DN: new diags
          ztmp1(1:kproma) = MAX(0._dp,pcdnc(1:kproma)-cqtmin) 
          ztmp2(1:kproma) = ztmp1(1:kproma)
          pevabf(1:kproma)   = MERGE(ztmp2(1:kproma),0._dp, ld_WBF(1:kproma))
          plwcbf(1:kproma)   = MERGE(pxlb(1:kproma)*paclc(1:kproma),0._dp, ld_WBF(1:kproma))
        !<<DN: new diags

          pcdnc_upd(1:kproma) = MERGE(cqtmin, pcdnc(1:kproma), ld_WBF(1:kproma))
         
          ztmp2(1:kproma) = pxib(1:kproma)+pxlb(1:kproma)
          pxib(1:kproma)  = MERGE(ztmp2(1:kproma), pxib(1:kproma), ld_WBF(1:kproma))

          !>>UPcomment: here WBF eats up all water within the cloud
          pxlb(1:kproma) = MERGE(0._dp, pxlb(1:kproma), ld_WBF(1:kproma))
          !<<UPcomment
  ENDIF
  !<<UP #765

END SUBROUTINE WBF_process

SUBROUTINE precip_formation_warm( &
              !--IN
              kbdim, kproma, &
              ld_prcp_warm, pauloc, paclc, pclcstar, prho, prho_rcp, pxrp1, pcdnc_min, &
              pcdnc, &
              !-- INOUT
              pcdnc_upd, pxlb, &
              !-- OUT
              pmratepr, prpr, prprn, pracl, pracln, praut, prautn, pcorrcd9 )!davidn

!DN --> SF ToDo: divide into two subroutines: one for Khairoutdinov and Kogan, 2000; and one for Beheng (1994)
!DN --> SF ToDo: autoconversion and accretion could be separated as well

  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL, INTENT(in) :: ld_prcp_warm(kbdim) !< physical condition for warm precip to form

  REAL(dp), INTENT(in) :: pauloc(kbdim)   !< Part of the grid box allowed to participate in accretion
                                          !< with newly formed condensate
  REAL(dp), INTENT(in) :: paclc(kbdim)    !< Cloud cover
  REAL(dp), INTENT(in) :: pclcstar(kbdim) !< Minimum of cloud cover and precipitation cover
  REAL(dp), INTENT(in) :: prho(kbdim)     !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: prho_rcp(kbdim) !< Inverse air density
  REAL(dp), INTENT(in) :: pxrp1(kbdim)    !< Rain mixing ratio (t) [kg/kg]
  REAL(dp), INTENT(in) :: pcdnc_min(kbdim)  !< Minimum CDNC concentration computed from maximum radius [1/m3] 
  REAL(dp), INTENT(in) :: pcdnc(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]

  REAL(dp), INTENT(inout) :: pcdnc_upd(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]
  !SF Note: intent(inout) important for pcdnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pxlb(kbdim)  !< Cloud liquid water in the cloudy part of the grid box [kg/kg]

  REAL(dp), INTENT(out) :: pmratepr(kbdim) !< Rain formation rate in cloudy part
  REAL(dp), INTENT(out) :: prpr(kbdim)     !< rain formation rate [kg/kg]
  REAL(dp), INTENT(out) :: prprn(kbdim)    !< Rain formation rate for number conc. [1/m3]
  REAL(dp), INTENT(out) :: pracl(kbdim),pracln(kbdim),prautn(kbdim),pcorrcd9(kbdim),praut(kbdim)    !davidn

  !local vars:
  REAL(dp) :: zrac1(kbdim)     !< Accretion of cloud water with rain from above [kg/kg]
  REAL(dp) :: zrac2(kbdim)     !< Accretion of cloud water with rain formed inside the grid box [kg/kg]
  REAL(dp) :: zraut(kbdim)     !< Autoconversion of cloud droplets [kg/kg]
  REAL(dp) :: zrautself(kbdim) !< Sum of autoconversion and self collection [kg/kg]
  REAL(dp) :: zcdnc_tmp(kbdim) !< Temporary CDNC needed for Beheng scheme only [1/m3]

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim)
  LOGICAL  :: ll1(kbdim)

  IF (nauto == 2) THEN
!          Autoconversion rate from Khairoutdinov and Kogan, 2000

     ztmp1(1:kproma) = ccraut*1350._dp*(1.e-6_dp*pcdnc(1:kproma))**(-1.79_dp)
     
     ztmp1(1:kproma) = pxlb(1:kproma) * (  1._dp &
                                        - (1._dp + ztmst*exm1_1*ztmp1(1:kproma)*pxlb(1:kproma)**exm1_1)**exp_1)

     ztmp1(1:kproma) = MIN(pxlb(1:kproma), ztmp1(1:kproma))
     !>>UP #805
     IF (lemuphase_cdautc) THEN
             zraut(1:kproma) = MERGE(eta_emu_cdautc*ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))
             ! Don't convert more liquid than there is
             zraut(1:kproma) = MIN(pxlb(1:kproma), zraut(1:kproma))
     ELSE
     !UP original
             zraut(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))
     ENDIF
     !<<UP #805
  
     
     ztmp1(1:kproma) = pxlb(1:kproma) - zraut(1:kproma)
     ztmp2(1:kproma) = pxlb(1:kproma) !SF keeps pxlb for later use
     pxlb(1:kproma)  = MERGE(ztmp1(1:kproma), pxlb(1:kproma), ld_prcp_warm(1:kproma))
 
!--- zrac1 is formed by accretion with rain from above
!--- zrac2 is formed by accretion with newly formed rain inside the grid box

     ztmp1(1:kproma) = -3.7_dp*ztmst*pxrp1(1:kproma)
     ztmp1(1:kproma) = EXP(ztmp1(1:kproma))
     ztmp1(1:kproma) = pxlb(1:kproma)*(1._dp-ztmp1(1:kproma))
     !>>UP #805
     IF (lemuphase_cdaccr) THEN
             zrac1(1:kproma) = MERGE(eta_emu_cdaccr*ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))
             ! Don't convert more liquid than there is
             zrac1(1:kproma) = MIN(pxlb(1:kproma), zrac1(1:kproma))
     ELSE
     !UP original
             zrac1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))
     ENDIF
     !<<UP #805

     pxlb(1:kproma) = pxlb(1:kproma) - zrac1(1:kproma)

     ztmp1(1:kproma) = -3.7_dp*ztmst*pauloc(1:kproma)*prho(1:kproma)*zraut(1:kproma)
     ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))
     ztmp1(1:kproma) = pxlb(1:kproma)*(1._dp-EXP(ztmp1(1:kproma)))
     !>>UP #805
     IF (lemuphase_cdaccr) THEN
             zrac2(1:kproma) = MERGE(eta_emu_cdaccr*ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))
             ! Don't convert more liquid than there is
             zrac2(1:kproma) = MIN(pxlb(1:kproma), zrac2(1:kproma))
     ELSE
     !UP original
             zrac2(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))
     ENDIF
     !<<UP #805

     pxlb(1:kproma) = pxlb(1:kproma) - zrac2(1:kproma)

     prpr(1:kproma) = paclc(1:kproma)    * (zraut(1:kproma)+zrac2(1:kproma)) &
                    + pclcstar(1:kproma) *  zrac1(1:kproma)
!>>DN
     praut(1:kproma) = paclc(1:kproma)*zraut(1:kproma)
     pracl(1:kproma) = pclcstar(1:kproma)*zrac1(1:kproma)+paclc(1:kproma)*zrac2(1:kproma)
!<<DN

!--- for in-cloud scavenging
     ztmp1(1:kproma)    = zraut(1:kproma)+zrac1(1:kproma)+zrac2(1:kproma)
     pmratepr(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))

!--- Autoconversion also changes the number of cloud droplets (prprn)

     ztmp1(1:kproma) = (zraut(1:kproma)+zrac1(1:kproma)+zrac2(1:kproma)) &
                     / (ztmp2(1:kproma)+eps)
     prprn(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))

!>>DN
!UP comment: I don't think I need to change pcdnc to pcdnc_upd here, because that has not been filled yet
     pracln(1:kproma) = MERGE((zrac1(1:kproma)+zrac2(1:kproma)) / (ztmp2(1:kproma)+eps) &
          * pcdnc(1:kproma),0._dp, ld_prcp_warm(1:kproma))
     prautn(1:kproma) = MERGE(zraut(1:kproma) / (ztmp2(1:kproma)+eps) * pcdnc(1:kproma),&
          0._dp, ld_prcp_warm(1:kproma))
!<<DN

     ll1(1:kproma) = ld_prcp_warm(1:kproma) &
               .AND. (pxlb(1:kproma) > cqtmin)

     ztmp1(1:kproma) = MERGE(pcdnc_min(1:kproma), 0._dp, ll1(1:kproma))
     ztmp1(1:kproma) = pcdnc(1:kproma)-ztmp1(1:kproma)
     ztmp2(1:kproma) = pcdnc(1:kproma)*prprn(1:kproma)

     ztmp3(1:kproma) = MIN(ztmp1(1:kproma), ztmp2(1:kproma))
!>>DN: new diags
     pcorrcd9(1:kproma) = MAX(ztmp2(1:kproma)-ztmp1(1:kproma), 0._dp)
!UP comment: this does not seem right to me, UP TODO understand
!<<DN: new diags        
     prprn(1:kproma) = MERGE(ztmp3(1:kproma), 0._dp, ld_prcp_warm(1:kproma))

     ztmp1(1:kproma) = pcdnc(1:kproma)-prprn(1:kproma)
     ztmp1(1:kproma) = MAX(ztmp1(1:kproma), cqtmin)
     pcdnc_upd(1:kproma) = MERGE(ztmp1(1:kproma), pcdnc(1:kproma), ld_prcp_warm(1:kproma))

!--- End included alternative autoconversion parameterisation ----------
!--- Changed for alternative autoconversion parameterisation -----------

  ELSE   !SF( nauto ==1)
!          Beheng (1994) - ECHAM 5 standard

!---Changed for prognostic CDNC/IC scheme ------------------------------
!   (Replaced pacdnc by pcdnc which is set above)
     ztmp1(1:kproma) = ccraut*1.2e27_dp * prho_rcp(1:kproma) & 
                     * (pcdnc(1:kproma)*1.e-6_dp)**(-3.3_dp) &
                     * (prho(1:kproma)*1.e-3_dp)**4.7_dp
     zraut(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma))
!--- End changed for prognostic CDNC/IC scheme -------------------------

     zraut(1:kproma) = pxlb(1:kproma) * (  1._dp &
                                        - (1._dp + ztmst*exm1_2*zraut(1:kproma)*pxlb(1:kproma)**exm1_2)**exp_2)

     !>>UP #805
     IF (lemuphase_cdautc) THEN
             zraut(1:kproma) = MIN(pxlb(1:kproma), eta_emu_cdautc*zraut(1:kproma))
     ELSE
     !UP original
             zraut(1:kproma) = MIN(pxlb(1:kproma), zraut(1:kproma))
     ENDIF
     !<<UP #805
  
!--- Included for prognostic CDNC/IC scheme ----------------------------
     ztmp1(1:kproma) = 7.7e9_dp * zraut(1:kproma) * prho(1:kproma) !SF zrautn
     ztmp2(1:kproma) = 1.289e10_dp * 1.e-6_dp * ztmst * (prho(1:kproma)*pxlb(1:kproma))**2 !SF zself (1.e-6 comes
                                                                                           !   from a unit bug fix)
     ztmp3(1:kproma)     = ztmp1(1:kproma) + ztmp2(1:kproma)
     ztmp3(1:kproma)     = MIN(ztmp3(1:kproma),pcdnc(1:kproma))
     zrautself(1:kproma) = MERGE(ztmp3(1:kproma), 0._dp, ld_prcp_warm(1:kproma))

     ztmp1(1:kproma) = pcdnc(1:kproma)-zrautself(1:kproma)
     ztmp1(1:kproma) = MAX(ztmp1(1:kproma), cqtmin)
     zcdnc_tmp(1:kproma) = MERGE(ztmp1(1:kproma), pcdnc(1:kproma), ld_prcp_warm(1:kproma))

!--- End included for CDNC/IC scheme -----------------------------------

     ztmp1(1:kproma) = pxlb(1:kproma) - zraut(1:kproma)
     pxlb(1:kproma)  = MERGE(ztmp1(1:kproma), pxlb(1:kproma), ld_prcp_warm(1:kproma))
!
!--- zrac1 is formed by accretion with rain from above
!--- zrac2 is formed by accretion with newly formed rain inside the grid box
! !UP: warm cloud process
!
     zrac1(1:kproma) = -6._dp*ztmst*pxrp1(1:kproma)
     zrac1(1:kproma) = EXP(zrac1(1:kproma))
     !>>UP #805
     IF (lemuphase_cdaccr) THEN
             zrac1(1:kproma) = eta_emu_cdaccr*pxlb(1:kproma)*(1._dp-zrac1(1:kproma))
             ! Don't convert more liquid than there is
             zrac1(1:kproma) = MIN(pxlb(1:kproma), zrac1(1:kproma))
     ELSE
     !UP original
             zrac1(1:kproma) = pxlb(1:kproma)*(1._dp-zrac1(1:kproma))
     ENDIF
     !<<UP #805

     ztmp1(1:kproma) = pxlb(1:kproma) - zrac1(1:kproma)
     ztmp2(1:kproma) = pxlb(1:kproma)  !SF keeps pxlb for later use
     pxlb(1:kproma)  = MERGE(ztmp1(1:kproma), pxlb(1:kproma), ld_prcp_warm(1:kproma))
 
     zrac2(1:kproma) = -6._dp*ztmst*pauloc(1:kproma)*prho(1:kproma)*zraut(1:kproma)
     zrac2(1:kproma) = EXP(zrac2(1:kproma))
     !>>UP #805
     IF (lemuphase_cdaccr) THEN
             zrac2(1:kproma) = eta_emu_cdaccr*pxlb(1:kproma)*(1._dp-zrac2(1:kproma))
             ! Don't convert more liquid than there is
             zrac2(1:kproma) = MIN(pxlb(1:kproma), zrac2(1:kproma))
     ELSE
     !UP original
             zrac2(1:kproma) = pxlb(1:kproma)*(1._dp-zrac2(1:kproma))
     ENDIF
     !<<UP #805

     ztmp1(1:kproma) = pxlb(1:kproma)-zrac2(1:kproma)
     pxlb(1:kproma)  = MERGE(ztmp1(1:kproma), pxlb(1:kproma), ld_prcp_warm(1:kproma))

     ztmp1(1:kproma) = paclc(1:kproma)    * (zraut(1:kproma)+zrac2(1:kproma)) & 
                     + pclcstar(1:kproma) * zrac1(1:kproma)
     prpr(1:kproma)  = MERGE(ztmp1(1:kproma), prpr(1:kproma), ld_prcp_warm(1:kproma))

!--- for in-cloud scavenging
     ztmp1(1:kproma)    = zraut(1:kproma) + zrac1(1:kproma) + zrac2(1:kproma)
     pmratepr(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_prcp_warm(1:kproma)) 

!--- Included for prognostic CDNC/IC scheme ----------------------------
!--- Autoconversion also changes the number of cloud droplets (prprn)
     ztmp1(1:kproma) = (zrac1(1:kproma)+zrac2(1:kproma))/(ztmp2(1:kproma)+eps) !SF ztmp2=pxlb+zrac1+zrac2
     ztmp1(1:kproma) = zcdnc_tmp(1:kproma)*ztmp1(1:kproma)
     ztmp1(1:kproma) = MIN(ztmp1(1:kproma), zcdnc_tmp(1:kproma)) !SF zraccn

     ztmp2(1:kproma) = zcdnc_tmp(1:kproma)-ztmp1(1:kproma)
     ztmp2(1:kproma) = MAX(ztmp2(1:kproma), cqtmin)
     pcdnc_upd(1:kproma) = MERGE(ztmp2(1:kproma), zcdnc_tmp(1:kproma), ld_prcp_warm(1:kproma))

     ztmp2(1:kproma) = zrautself(1:kproma) + ztmp1(1:kproma)
     prprn(1:kproma) = MERGE(ztmp2(1:kproma), 0._dp, ld_prcp_warm(1:kproma))
!--- End included for CDNC/IC scheme -----------------------------------
  ENDIF 

END SUBROUTINE precip_formation_warm

SUBROUTINE precip_formation_cold( &
              !-- IN
              kbdim, kproma, &
              ld_cc, pauloc, paclc, pclcstar, pqrho, prho_rcp, ptp1tmp, &
              pviscos, pxsp1, prho, pcdnc_min, pcdnc, picnc, &
              pjk, pkrow, & !UP #821
              pcmpsimpl_rime, pcmpsimpl_icaccr, pcmpsimpl_sci, & !UP #821
              !-- INOUT
              picnc_upd, pcdnc_upd, pmrateps, &
              pxib, pxlb, & 
              !-- OUT
              psprn, psacl, psacln, pmsnowacl, pspr, psaut, psaci, paggn, psacin, pselfn, psecprodn, psecprod)!davidn

!>>UP #777
USE mo_cmp_timer,          ONLY: & ! already imported for the module: timer_start, timer_stop
                                 timer_cmp_7ccpfcaggr, timer_cmp_7ccpfcaccr,    &
                                 timer_cmp_7ccpfcrime, timer_cmp_7ccpfcsepr,    &
                                 timer_cmp_7ccpfcsci,                           &
                                 timer_cmp_7ccpfcfinish, timer_cmp_7ccpfcinside
!<<UP #777
!>>UP #821
USE mo_cmp_diagn,           ONLY: diag_delta_sci, diag_cntr_sci, &
                                  diag_delta_rime, diag_cntr_rime, &
                                  diag_delta_icaccr, diag_cntr_icaccr
!<<UP #821

  INTEGER, INTENT(in) :: kbdim, kproma

  LOGICAL, INTENT(in)  :: ld_cc(kbdim)

  REAL(dp), INTENT(in) :: pauloc(kbdim) !< Part of the grid box allowed to particiation in accretion with newly formed condensate
  !UP: i.e. riming
  REAL(dp), INTENT(in) :: paclc(kbdim)    !< Cloud cover
  REAL(dp), INTENT(in) :: pclcstar(kbdim) !< Minimum of cloud cover and precipitation cover
  REAL(dp), INTENT(in) :: pqrho(kbdim)    !< Inverse air density [m3/kg]
  REAL(dp), INTENT(in) :: prho_rcp(kbdim) !< Inverse air density
  REAL(dp), INTENT(in) :: ptp1tmp(kbdim)  !< temporary value of the updated temperature (t) [K]
  REAL(dp), INTENT(in) :: pviscos(kbdim)  !< Dynamic viscosity of water in air
  REAL(dp), INTENT(in) :: pxsp1(kbdim)    !< ??
  REAL(dp), INTENT(in) :: prho(kbdim)     !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: pcdnc_min(kbdim)  !< Minimum CDNC concentration computed from maximum radius [1/m3] 
  REAL(dp), INTENT(in) :: pcdnc(kbdim)    !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp), INTENT(in) :: picnc(kbdim)    !< Ice crystal number concentration (ICNC) [1/m3]
  !>>UP #821
  INTEGER, INTENT(in) :: pjk !jk for use in this subroutine
  INTEGER, INTENT(in) :: pkrow ! krow for use in this subroutine
  REAL(dp), INTENT(in) :: pcmpsimpl_rime(kbdim)     ! receives input for prescribed riming
  REAL(dp), INTENT(in) :: pcmpsimpl_icaccr(kbdim)   ! receives input for prescribed ICaccr
  REAL(dp), INTENT(in) :: pcmpsimpl_sci(kbdim)      ! receives input for prescribed sci
  !<<UP #821

  REAL(dp), INTENT(inout) :: picnc_upd(kbdim)    !< Ice crystal number concentration (ICNC) [1/m3]
  !SF Note: intent(inout) important for picnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pcdnc_upd(kbdim)    !< Cloud droplet number concentration (CDNC) [1/m3]
  !SF Note: intent(inout) important for pcdnc_upd, see my note in sub 'melting_snow_and_ice'
  REAL(dp), INTENT(inout) :: pmrateps(kbdim) !< Snow formation rate in cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(inout) :: pxib(kbdim)     !< Cloud ice in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(inout) :: pxlb(kbdim)     !< Cloud liquid water in the cloudy part of the grid box [kg/kg]

  REAL(dp), INTENT(out) :: psprn(kbdim)     !< Snow formation rate for number conc. [1/m3]
  REAL(dp), INTENT(out) :: psacl(kbdim)     !< accretion of snow flakes with cloud droplets [kg/kg] !UP: i.e. riming
  REAL(dp), INTENT(out) :: psacln(kbdim)    !< Accretion rate for number conc. [1/m3]
  REAL(dp), INTENT(out) :: pmsnowacl(kbdim) !< Accretion rate of snow with cloud droplets in cloudy 
  !UP: i.e. riming 
  REAL(dp), INTENT(out) :: pspr(kbdim)      !< snow formation rate [kg/kg]
  REAL(dp), INTENT(out) :: psaut(kbdim),psaci(kbdim),paggn(kbdim),psacin(kbdim), &!davidn
                           pselfn(kbdim),psecprodn(kbdim),psecprod(kbdim) !davidn

  !local vars:
  REAL(dp) :: zc1(kbdim)       !< Temporary variable needed for aggregation
  REAL(dp) :: zcolleffi(kbdim) !< Collision efficiency for aggregation
  !UP comment: to me it looks like the above zcolleffi is only used for accretion, not for aggregation
  REAL(dp) :: zcsacl(kbdim)    !< Temporary variable needed for accretion of snow with cloud droplets !UP: i.e. riming
  REAL(dp) :: zdplanar(kbdim)  !< Diameter of planar ice crystals
  REAL(dp) :: zrey(kbdim)      !< Reynolds number
  REAL(dp) :: zris(kbdim)      !< Size of ice crystals
  REAL(dp) :: zsaci(kbdim)     !< Accretion of snow with ice crystals [kg/kg]
!UP: i.e. real accretion
  REAL(dp) :: zsaclin(kbdim)   !< Accretion of snow crystals with cloud droplets [1/m3] !UP: i.e. riming
  REAL(dp) :: zsaut(kbdim)     !< Aggregation of ice crystals [kg/kg]
  REAL(dp) :: zsecprod(kbdim)  !< Secondary ice crystal production
  REAL(dp) :: zstcrit(kbdim)   !< Critical Stokes number
  REAL(dp) :: zstokes(kbdim)   !< Stokes number
  REAL(dp) :: zudrop(kbdim)    !< Cloud droplet vertical velocity [m/s]
  REAL(dp) :: zusnow(kbdim)    !< Fall velocity of snow [m/s]
  REAL(dp) :: zxibold(kbdim)   !< Temporary value of cloud ice [kg/kg]
  REAL(dp) :: zxsp(kbdim)      !< Snow mass mixing ratio [kg/kg]
  REAL(dp) :: zxsp2(kbdim)     !< Snow being formed inside the grid box [kg/kg]

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim), ztmp4(kbdim), ztmp5(kbdim)
  LOGICAL  :: ll1(kbdim), ll2(kbdim), ll3(kbdim), ll4(kbdim), ll5(kbdim), ll6(kbdim), &
              ll7(kbdim), ll8(kbdim) 
  !>>UP #821
  INTEGER :: ierr
  !<<UP #821

!SFToDo: this code needs thorough refactoring:
!          - clean separation between each process
!          - fixes for consistency
!        The fixes are however subject to discussion (e.g. #454)
! The refactoring will be done when the whole topic becomes clearer

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_start(timer_cmp_7ccpfcinside)
     CALL timer_start(timer_cmp_7ccpfcfinish)
  END IF
  !<<UP

  zxibold(1:kproma) = MAX(pxib(1:kproma),eps) !SF store pxib (with security) before any 
                                              !   cold-precip-related change for later use

  !>>SF needed for the temporary speedup right below
  pspr(1:kproma)      = 0._dp
  psprn(1:kproma)     = 0._dp
  pmsnowacl(1:kproma) = 0._dp
  psacl(1:kproma)     = 0._dp
  psacln(1:kproma)    = 0._dp
  !<<SF
!davidn
  psaut(1:kproma) = 0._dp
  psaci(1:kproma) = 0._dp
  psecprod(1:kproma) = 0._dp
  paggn(1:kproma)     = 0._dp
  psacin(1:kproma)    = 0._dp
  pselfn(1:kproma)    = 0._dp
  psecprodn(1:kproma) = 0._dp
!<<SFdiag
!davidn

  ll1(1:kproma) = ld_cc(1:kproma) .AND. (pxib(1:kproma) > cqtmin)
  !UP comment: needed for aggregation

!>>SF temporary speedup in waiting for the refactoring (the condition might then be moved before the
!     call to this subroutine)
  IF (.NOT. ANY(ll1(1:kproma))) return
!<<SF

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_7ccpfcfinish)
     CALL timer_start(timer_cmp_7ccpfcaggr)
  END IF
  !<<UP

  !UP comment: calculation of ice crystal size (zris) needed for aggregation starts here
!>>SF #176 now this is done consistenly between cl. micro and radiation, with a general formula
  !SF conversion of in-cloud ice mmr from kg/kg to in-cloud g m-3:
  ztmp1(1:kproma) = 1000._dp*pxib(1:kproma)*prho(1:kproma)

  ztmp1(:) = eff_ice_crystal_radius(kbdim, kproma, ztmp1(:), picnc(:))
!<<SF

  ztmp1(1:kproma) = MIN(MAX(ztmp1(1:kproma), ceffmin), ceffmax) !SF zrieff

  ztmp2(1:kproma) = 5113188._dp+2809._dp*ztmp1(1:kproma)**3
  ztmp2(1:kproma) = SQRT(ztmp2(1:kproma))
  ztmp2(1:kproma) = -2261._dp + ztmp2(1:kproma)  !SF zrih
  
  ztmp3(1:kproma) = 1.e-6_dp*ztmp2(1:kproma)**(1._dp/3._dp)
  zris(1:kproma)  = MERGE(ztmp3(1:kproma), 1._dp, ll1(1:kproma))   !SF 1. could be whatever, just not 0.
  !UP comment: end of size computation needed for aggregation
 
!--- temperature dependent collision efficiency (zcolleffi)

  ztmp1(1:kproma)     = fact_coll_eff*(ptp1tmp(1:kproma)-tmelt) !SF See #471
  ztmp1(1:kproma)     = EXP(ztmp1(1:kproma))
  zcolleffi(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll1(1:kproma)) 

  !UP comment: aggregation of ice crystals starts here
  zc1(1:kproma) = 17.5_dp / crhoi * prho(1:kproma) * pqrho(1:kproma)**0.33_dp

  ztmp1(1:kproma) = -6._dp / zc1(1:kproma) * LOG10(1.e4_dp*zris(1:kproma)) !SF zdt2
!SF Note: 1.e-4 = minimum size of snow flake
  ztmp1(1:kproma) = ccsaut / ztmp1(1:kproma)
 
  ztmp1(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll1(1:kproma))
  !>>UP lemuphase_icaggr: phase aggregation of ice cyrstals
  IF (lemuphase_icaggr) THEN
     zsaut(1:kproma) = MIN(eta_emu_icaggr * pxib(1:kproma) * &
        (1._dp - 1._dp/(1._dp+ztmp1(1:kproma)*ztmst*pxib(1:kproma))), pxib(1:kproma))
     !UP #776: need a maximal condition (aggregation only as much as there are ice crystals)
  ELSE
        zsaut(1:kproma) = pxib(1:kproma) * (1._dp - 1._dp/(1._dp+ztmp1(1:kproma)*ztmst*pxib(1:kproma))) ! original line
  ENDIF !lemuphase_icaggr
  !<<UP lemuphase_icaggr: phase aggregation of ice cyrstals

  ztmp1(1:kproma)   = pxib(1:kproma) - zsaut(1:kproma)
  zxibold(1:kproma) = pxib(1:kproma) !SF store pxib for later use at the end of 7.2
  pxib(1:kproma)    = MERGE(ztmp1(1:kproma), pxib(1:kproma), ll1(1:kproma))

  ztmp1(1:kproma) = pauloc(1:kproma)*prho(1:kproma)*zsaut(1:kproma)
  zxsp2(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll1(1:kproma))  !snow that is formed inside the grid box

  zxsp(1:kproma) = pxsp1(1:kproma) + zxsp2(1:kproma)
  
  !UP comment: aggregation of ice crystals at least until here

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_7ccpfcaggr)
     CALL timer_start(timer_cmp_7ccpfcrime)
  END IF
  !<<UP

  ll2(1:kproma) = ll1(1:kproma)               .AND. &
              (zxsp(1:kproma)  >  cqtmin) .AND. &
              (pxlb(1:kproma)  >  cqtmin) .AND. &
              (pcdnc(1:kproma) >= pcdnc_min(1:kproma))

!--- The riming of snow with droplets is calculated assuming planar snow flakes (Lohmann, JAS, 2004)
!--- It depends on the droplet (zudrop) and snow flake (zusnow) fall velocity, the Stokes
!--- number (zstokes) and the Reynolds number (zrey)

  ztmp1(1:kproma) = ( 6._dp*pirho_rcp*prho(1:kproma)*pxlb(1:kproma)/pcdnc(1:kproma) )**(1._dp/3._dp)
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 1.e-6_dp) !SF zdw

  zudrop(1:kproma) = 1.19e4_dp*2500._dp*ztmp1(1:kproma)**2*(1.3_dp*prho_rcp(1:kproma))**0.35_dp

  zdplanar(1:kproma) = 447.e-6_dp !DN #454: constant 100 mum volume radius -> constant max dimension of 447 mum
                                  !         after Table 2.2a in Pruppacher&Klett 1997
     
  zusnow(1:kproma) = 2.34_dp * (100._dp*zdplanar(1:kproma))**0.3_dp &
                   * (1.3_dp*prho_rcp(1:kproma))**0.35_dp

  zstokes(1:kproma) = 2._dp*rgrav*(zusnow(1:kproma)-zudrop(1:kproma))*zudrop(1:kproma)/zdplanar(1:kproma)
  zstokes(1:kproma) = MAX(zstokes(1:kproma), cqtmin)

  zrey(1:kproma) = prho(1:kproma)*zdplanar(1:kproma)*zusnow(1:kproma)/pviscos(1:kproma)
  zrey(1:kproma) = MAX(zrey(1:kproma),cqtmin)

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

  ztmp1(1:kproma) = cons4*zxsp(1:kproma)**0.8125_dp !SF zlamsm
  ztmp1(1:kproma) = pi*cn0s*3.078_dp*ztmp1(1:kproma)*pqrho(1:kproma)**0.5_dp
  ztmp2(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma)) 

  ztmp1(1:kproma) = -ztmst*ztmp2(1:kproma)*zcsacl(1:kproma)
  ztmp1(1:kproma) = EXP(ztmp1(1:kproma))
  !>>UP #763
  !>>UP #821
  IF (lcmpsimpl_prescr) THEN
      SELECT CASE(ncmpsimpl_prescr_rime)
          CASE(0)
              ! original line    
              ztmp1(1:kproma) = pxlb(1:kproma)*(1._dp-ztmp1(1:kproma))
          CASE(1,2)
              !UP #776: need to put in a maximum condition
              !         (not more riming that what is there as liquid)
              ! ztmp1 is used only with ll2 after -> use with ll2 here
              ztmp1(1:kproma) = MERGE(MIN(pcmpsimpl_rime(1:kproma), pxlb(1:kproma)), &
                                        0._dp, ll2(1:kproma))
      END SELECT
  !<<UP #821
  ! by phasing ztmp1 here I can get all the other cases in one go
  ELSE IF (lemuphase_riming) THEN
          ztmp1(1:kproma) = MIN(eta_emu_riming*pxlb(1:kproma)*(1._dp-ztmp1(1:kproma)), &
                                pxlb(1:kproma))
          !UP #776: need to put in a maximum condition (not more riming that what is there as liquid)
  ELSE
          ! original line    
          ztmp1(1:kproma) = pxlb(1:kproma)*(1._dp-ztmp1(1:kproma))
  ENDIF
  !<<UP #763
  !>>UP #821
  !UP TODO: don't I need the ll2 here as well?
  diag_delta_rime(1:kproma,pjk,pkrow) = diag_delta_rime(1:kproma,pjk,pkrow) + &
                                     MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma))
  ! only diagnose a count when it's also really used -> add logical merge
  diag_cntr_rime(1:kproma,pjk,pkrow) = diag_cntr_rime(1:kproma,pjk,pkrow) + MERGE(1._dp, 0._dp, ll2(1:kproma))
  !<<UP #821

  ! UP: automatically phased through ztmp1 being phased
  zsaclin(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma))

  ! UP: automatically phased through ztmp1 being phased
  ztmp2(1:kproma) = pxlb(1:kproma)-ztmp1(1:kproma)

  ztmp3(1:kproma) = pxlb(1:kproma)  !SF keeps pxlb for later use
  ! UP: riming gets phased automatically here, because ztmp1/ztmp2 is phased
  pxlb(1:kproma)     = MERGE(ztmp2(1:kproma), pxlb(1:kproma), ll2(1:kproma)) ! UP: original line

  ztmp2(1:kproma) = paclc(1:kproma)*ztmp1(1:kproma)
  ! UP: automatically phased through ztmp1 being phased
  psacl(1:kproma)   = MERGE(ztmp2(1:kproma), 0._dp, ll2(1:kproma)) ! UP: original line

  ll2(1:kproma) = (pxlb(1:kproma) > cqtmin) 

  ztmp1(1:kproma) = pcdnc(1:kproma)*zsaclin(1:kproma)/(ztmp3(1:kproma)+eps) !SF ztmp3 =pxlb+zsacl2in
  ztmp1(1:kproma) = MIN(ztmp1(1:kproma), pcdnc(1:kproma)-pcdnc_min(1:kproma))
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0._dp) 

  ! UP: riming phased here through zsaclin, which enters ztmp1; before I was
  ! phasing ztmp1, but I think that is less physical
  psacln(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma)) ! UP: original line 

!SF end ll2 condition (pxlb > cqtmin)

  ztmp1(1:kproma)     = pcdnc(1:kproma)-psacln(1:kproma)
  pcdnc_upd(1:kproma) = MERGE(ztmp1(1:kproma), pcdnc(1:kproma), ll1(1:kproma)) ! UP: phase is achieved here already by psacln above
  ! UP: riming is phased here through zsaclin
  pmsnowacl(1:kproma) = MERGE(zsaclin(1:kproma), 0._dp, ll1(1:kproma)) ! UP: original line

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_7ccpfcrime)
     CALL timer_start(timer_cmp_7ccpfcaccr)
  END IF
  !<<UP

!UP: I think this is where accretion of snow and ice starts, Issue #775: remove this from the riming phasing and make its own
!SF end ll2 condition:
!         (ll1(1:kproma) .AND. (zxsp(1:kproma) > cqtmin) .AND. &
!         (pxlb(1:kproma) > cqtmin) .AND. (pcdnc(1:kproma) > pcdnc_min(1:kproma)))     

  ll1(1:kproma) = ld_cc(1:kproma) .AND. (pxib(1:kproma) > cqtmin)

  ll2(1:kproma) = ll1(1:kproma) .AND. (zxsp(1:kproma) > cqtmin) 

  ztmp1(1:kproma) = cons4*zxsp(1:kproma)**0.8125_dp
  ztmp1(1:kproma) = pi*cn0s*3.078_dp*ztmp1(1:kproma)*pqrho(1:kproma)**0.5_dp
  ztmp1(1:kproma) = -ztmst*ztmp1(1:kproma)*zcolleffi(1:kproma)
  ztmp1(1:kproma) = EXP(ztmp1(1:kproma))
  ztmp1(1:kproma) = pxib(1:kproma) * (1._dp-ztmp1(1:kproma))

     !>>UP #821
     IF (lcmpsimpl_prescr) THEN
        SELECT CASE(ncmpsimpl_prescr_icaccr)
            CASE(0)
                !UP original line
                zsaci(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma)) 
            CASE(1,2)
                !UP #776: need to put in a maximum condition (not larger than pxib)
                zsaci(1:kproma) = MERGE(MIN(pcmpsimpl_icaccr(1:kproma), pxib(1:kproma)), 0._dp, ll2(1:kproma))
        END SELECT
     !<<UP #821
     ELSE !Only makes sense to phase if we don't prescribe
          ! UP: could phase zsaci here directly -> easier/more transparent than trying to
          ! find every instance where it is used
          !>>UP lemuphase_icaccr
          IF (lemuphase_icaccr) THEN
                  zsaci(1:kproma) = MERGE(MIN(eta_emu_icaccr*ztmp1(1:kproma),pxib(1:kproma)), 0._dp, ll2(1:kproma)) 
                  !UP #776: need to put in a maximum condition (not larger than pxib)
          ELSE
                  zsaci(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll2(1:kproma)) 
          ENDIF
          !<<UP lemuphase_icaccr
     ENDIF
  !>>UP #821
  diag_delta_icaccr(1:kproma,pjk,pkrow) = diag_delta_icaccr(1:kproma,pjk,pkrow) + zsaci(1:kproma)
! above contains the logical implicitly in zsaci
  diag_cntr_icaccr(1:kproma,pjk,pkrow) = diag_cntr_icaccr(1:kproma,pjk,pkrow) + &
                                         MERGE(1._dp, 0._dp, ll2(1:kproma))
  !<<UP #821

  !>>UP lemuphase_icaccr automatically through zsaci
  pxib(1:kproma) = pxib(1:kproma)-zsaci(1:kproma) ! UP: original line

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_7ccpfcaccr)
     CALL timer_start(timer_cmp_7ccpfcfinish)
  END IF
  !<<UP

!SF end ll2 condition:
!       ll2(1:kproma) = ll1(1:kproma).AND.(zxsp(1:kproma) > cqtmin) 

  ztmp1(1:kproma) = paclc(1:kproma)*( zsaut(1:kproma)+zsaci(1:kproma) )
  !UP #776: do I need a maximum condition here, so that the two processes together don't use more ice than is there?
  !UP: I don't think so, because zsaci acts only on the pxib that is still there after zsaut
  !>>UP lemuphase_icaccr and lemuphase_icaggr automatically through zsaut and zsaci
  pspr(1:kproma)  = MERGE(ztmp1(1:kproma), 0._dp, ll1(1:kproma)) ! UP: original line 

!--- for in-cloud scavenging
  ztmp1(1:kproma)    = zsaut(1:kproma)+zsaci(1:kproma)
  !>>UP lemuphase_icaggr and lemuphase_icaccr: phasing of ztmp1 automatic through
  !zsaut and zsaci, but phasing of pmrateps more tricky
  IF (lemuphase_icaggr .AND. lemuphase_icaccr) THEN
        ! New formulation after issue 766:
          pmrateps(1:kproma) = &
              MERGE(MAX((1._dp-(eta_emu_icaccr+eta_emu_icaggr)/2)*pmrateps(1:kproma) &
              + ztmp1(1:kproma), 0._dp), &
              pmrateps(1:kproma), ll1(1:kproma))
      !UP #776: need to put in a maximum condition so that (1-eta*) does not get smaller than 0
      !UP: for the double pmrateps there's no perfect solution for phasing as
      !ztmp1 is influenced by both but it's not clear how much at this point
      !(depends on the contribution from zsaut and zsaci)
      !1/2 is a pragmatic choice and we've also seen that this particular part
      !does not influence the result much
  ELSEIF (lemuphase_icaccr) THEN
          pmrateps(1:kproma) = &
              MERGE(MAX((1._dp-eta_emu_icaccr)*pmrateps(1:kproma) + ztmp1(1:kproma), 0._dp), &
              pmrateps(1:kproma), ll1(1:kproma))
  ELSEIF (lemuphase_icaggr) THEN
          pmrateps(1:kproma) = &
              MERGE(MAX((1._dp-eta_emu_icaggr)*pmrateps(1:kproma) + ztmp1(1:kproma), 0._dp), &
              pmrateps(1:kproma), ll1(1:kproma))
  ELSE
          pmrateps(1:kproma) = MERGE(ztmp1(1:kproma), pmrateps(1:kproma), ll1(1:kproma)) ! UP: original line
  ENDIF 
  !<<UP lemuphase_icaccr and lemuphase_icaggr

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_7ccpfcfinish)
     CALL timer_start(timer_cmp_7ccpfcsepr)
  END IF
  !<<UP

!       secondary ice crystal production (zsecprod) after Levkov et al. 1992
!       sink for snow, source for ice crystals
!uls    included size depedent accretion rate (Lohmann, JAS, 2004)

  zsecprod(1:kproma) = 0._dp

  IF (lsecprod) THEN !SF switch on/off secondary ice prod (see #251)

!SF ToDo: code block to refactor to avoid code duplication as compared to riming calcs...

     ll2(1:kproma) = (zxsp(1:kproma) > epsec) &
               .AND. (pxlb(1:kproma)    > epsec)  &
               .AND. (ptp1tmp(1:kproma) > 265.2_dp) &
               .AND. (ptp1tmp(1:kproma) < 270.2_dp)
   
     ztmp1(1:kproma) = ( 6._dp*pirho_rcp*prho(1:kproma) &
                       * pxlb(1:kproma)/pcdnc(1:kproma) )**(1._dp/3._dp)
     ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 1.e-6_dp) !SF zdw
   
     zudrop(1:kproma) = 1.19e4_dp * (50._dp*ztmp1(1:kproma))**2 &
                      * (1.3_dp*prho_rcp(1:kproma))**0.35_dp
    
     zstokes(1:kproma) = 2._dp * rgrav * (zusnow(1:kproma)-zudrop(1:kproma)) &
                       * zudrop(1:kproma) / zdplanar(1:kproma)
     zstokes(1:kproma) = MAX(zstokes(1:kproma), cqtmin)
   
     ztmp1(1:kproma) = 0.2_dp * (LOG10(zstokes(1:kproma))-LOG10(zstcrit(1:kproma))- 2.236_dp )**2
     ztmp1(1:kproma) = MIN(ztmp1(1:kproma), 1._dp-cqtmin)
     ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0._dp)
     ztmp1(1:kproma) = SQRT(1._dp - ztmp1(1:kproma))
   
     ll6(1:kproma) = ll5(1:kproma) .AND. (zstokes(1:kproma) <= 0.06_dp)
     ll7(1:kproma) = ll5(1:kproma) .AND. (zstokes(1:kproma) >  0.06_dp) &
                                   .AND. (zstokes(1:kproma) <= 0.25_dp)
     ll8(1:kproma) = ll5(1:kproma) .AND. (zstokes(1:kproma) >  0.25_dp) &
                                   .AND. (zstokes(1:kproma) <= 1.00_dp)
   
!SF Note: WHERE statements should be eliminated, because they are extremely costly
!         I nevertheless leave these for now, because this block (lsecprod==.true.) will be completely refactored
!         and beside it is almost never used.

     WHERE (ll6(1:kproma))
           ztmp1(1:kproma) = 1.034_dp*zstokes(1:kproma)**1.085_dp
     ELSEWHERE (ll7(1:kproma))
           ztmp1(1:kproma) = 0.787_dp*zstokes(1:kproma)**0.988_dp
     ELSEWHERE (ll8(1:kproma))
           ztmp1(1:kproma) = 0.7475_dp*LOG10(zstokes(1:kproma))+0.65_dp
     ELSEWHERE (ll5(1:kproma))
           ztmp1(1:kproma) = (zstokes(1:kproma)+1.1_dp)**2/(zstokes(1:kproma)+1.6_dp)**2
     ENDWHERE
   
     ztmp1(1:kproma)  = MAX(MIN(ztmp1(1:kproma), 1._dp), 0.01_dp)
     zcsacl(1:kproma) = MERGE(ztmp1(1:kproma), zcsacl(1:kproma), ll2(1:kproma))
   
     ztmp1(1:kproma) = cons5*zxsp(1:kproma)**0.875_dp !SF zlams2
   
     ztmp2(1:kproma) = cn0s * 0.831_dp * pi / mw0 &
                     * zcsacl(1:kproma) * prho(1:kproma) * pxlb(1:kproma) * ztmp1(1:kproma) &
                     * ( grav*crhosno / (0.75_dp*cdi*prho(1:kproma)) )**0.5_dp   !SF zj
 
     ztmp2(1:kproma) = MAX(0.00285_dp*ztmp2(1:kproma), 0._dp)  !SF zpn
   
     ztmp2(1:kproma) = ztmst*mi0*ztmp2(1:kproma)*prho_rcp(1:kproma)
     ztmp3(1:kproma) = zxsp(1:kproma)*prho_rcp(1:kproma)
     ztmp2(1:kproma) = MIN(ztmp3(1:kproma), ztmp2(1:kproma))
     ztmp2(1:kproma) = MAX(ztmp2(1:kproma), 0._dp)
   
     !>>UP #805
     IF (lemuphase_sip) THEN
             zsecprod(1:kproma) = MERGE( MAX(((1._dp - eta_emu_sip)*zsecprod(1:kproma) &
                                                + eta_emu_sip*ztmp2(1:kproma)), 0._dp), &
                                                zsecprod(1:kproma), ll2(1:kproma))
     ELSE
             zsecprod(1:kproma) = MERGE(ztmp2(1:kproma), zsecprod(1:kproma), ll2(1:kproma))
     END IF
     !<<UP #805
   
     ztmp1(1:kproma) = pxib(1:kproma)+zsecprod(1:kproma)
     pxib(1:kproma)  = MERGE(ztmp1(1:kproma), pxib(1:kproma), ll2(1:kproma))
   
     ztmp1(1:kproma) = pspr(1:kproma)-pclcstar(1:kproma)*zsecprod(1:kproma)
     ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0._dp)
     pspr(1:kproma)  = MERGE(ztmp1(1:kproma), pspr(1:kproma), ll2(1:kproma))
   
     ztmp1(1:kproma)    = pmrateps(1:kproma)-zsecprod(1:kproma)
     pmrateps(1:kproma) = MERGE(ztmp1(1:kproma), pmrateps(1:kproma), ll2(1:kproma))
   
!SF end ll2 condition (secondary ice production)
  ENDIF !SF end lsecprod

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_7ccpfcsepr)
     CALL timer_start(timer_cmp_7ccpfcfinish)
  END IF
  !<<UP

!davidn
  psaut(1:kproma) = paclc(1:kproma)*zsaut(1:kproma)
  psaci(1:kproma) = paclc(1:kproma)*zsaci(1:kproma)
  psecprod(1:kproma) = zsecprod(1:kproma)
!davidn

!SF end ll1 condition (ld_cc and pxib(jl) > cqtmin) 
!
!--- Also change the number of ice crystals due to the break-up of snow flakes
!
  ll1(1:kproma) = ld_cc(1:kproma) &
            .AND. (pxib(1:kproma)  >  epsec) &
            .AND. (picnc(1:kproma) >= icemin)
 
  ztmp1(1:kproma) = zxibold(1:kproma)
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0._dp) !SF zxibold

! UP: phasing accretion here indirectly via zsaci (same as aggregation via
! zsaut)
  ztmp2(1:kproma) = picnc(1:kproma) * (zsaci(1:kproma)+zsaut(1:kproma)) / (ztmp1(1:kproma)+eps) !SF zsprn1

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_7ccpfcfinish)
     CALL timer_start(timer_cmp_7ccpfcsci)
  END IF
  !<<UP

!UP: phase self-collection here, Issue #779
  !>>UP #821
  IF (lcmpsimpl_prescr) THEN
        SELECT CASE(ncmpsimpl_prescr_sci)
                CASE(0)
                  ztmp3(1:kproma) = 0.5_dp * ztmst * zc1(1:kproma) * picnc(1:kproma) * pxib(1:kproma) !SF zself !UP: original line
                CASE(1,2)
                  ! ztmp3 is used only with ll1 after -> use ll1 everywhere here
                  ztmp3(1:kproma) = MERGE(MIN(pcmpsimpl_sci(1:kproma), picnc(1:kproma)), 0._dp, ll1(1:kproma))
                  !UP #776: need to put in a maximum condition (not more self-collection than there is picnc)
        END SELECT
  !<<UP #821
  !>>UP lemuphase_sci
  ELSE IF (lemuphase_sci) THEN
     ztmp3(1:kproma) = MIN(eta_emu_sci*0.5_dp*ztmst*zc1(1:kproma)*picnc(1:kproma)*pxib(1:kproma), &
       picnc(1:kproma)) !SF zself
          !UP #776: need to put in a maximum condition (not more self-collection than there is picnc)
  ELSE
          ztmp3(1:kproma) = 0.5_dp * ztmst * zc1(1:kproma) * picnc(1:kproma) * pxib(1:kproma) !SF zself !UP: original line
  ENDIF
  !>>UP #821
  !UP TODO: doesn't below need the ll1 merge?
  diag_delta_sci(1:kproma,pjk,pkrow) = diag_delta_sci(1:kproma,pjk,pkrow) + MERGE(ztmp3(1:kproma), &
                                        0._dp, ll1(1:kproma))
  diag_cntr_sci(1:kproma,pjk,pkrow) = diag_cntr_sci(1:kproma,pjk,pkrow) + &
                                      MERGE(1._dp, 0._dp, ll1(1:kproma))
  !<<UP #821
  !<<UP lemuphase_sci

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_7ccpfcsci)
     CALL timer_start(timer_cmp_7ccpfcfinish)
  END IF
  !<<UP

  ztmp4(1:kproma) = mi0_rcp * prho(1:kproma) * zsecprod(1:kproma) !SF zsecprodn

  ztmp5(1:kproma) = ztmp2(1:kproma) + ztmp3(1:kproma) - ztmp4(1:kproma) 
  ztmp5(1:kproma) = MIN(ztmp5(1:kproma), picnc(1:kproma))    !SF zsprnself
!>>SFdiag
  paggn(1:kproma)     = picnc(1:kproma) * zsaut(1:kproma) / (ztmp1(1:kproma)+eps)
  paggn(1:kproma)     = MERGE(paggn(1:kproma),0._dp,ll1(1:kproma))
  psacin(1:kproma)    = picnc(1:kproma) * zsaci(1:kproma) / (ztmp1(1:kproma)+eps)
  psacin(1:kproma)    = MERGE(psacin(1:kproma),0._dp,ll1(1:kproma))
  pselfn(1:kproma)    = MERGE(ztmp3(1:kproma),0._dp,ll1(1:kproma))
  psecprodn(1:kproma) = MERGE(ztmp4(1:kproma),0._dp,ll1(1:kproma))
!<<SFdiag
   psprn(1:kproma) = MERGE(ztmp5(1:kproma), 0._dp, ll1(1:kproma)) !UP warning: self-collection wrongly included

  ztmp1(1:kproma)     = picnc(1:kproma) - psprn(1:kproma)
  ztmp1(1:kproma)     = MAX(ztmp1(1:kproma), cqtmin)
  picnc_upd(1:kproma) = MERGE(ztmp1(1:kproma), picnc(1:kproma), ll1(1:kproma))

  !>>UP timers #777
  IF (ltimer) THEN
     CALL timer_stop(timer_cmp_7ccpfcfinish)
     CALL timer_stop(timer_cmp_7ccpfcinside)
  END IF
  !<<UP

END SUBROUTINE precip_formation_cold

SUBROUTINE update_precip_fluxes( &
              !--IN
              kbdim, kproma, &
              kk, klev, paclc, pdp, pevp, plsdcp, plvdcp, prpr, psacl, pspr, psub, &
              ptp1tmp, pxiflux, & 
              !--INOUT
              pclcpre, prfl, psfl, psmlt, &
              !--OUT
              pfevapr, pfrain, pfsnow, pfsubls, ptmp2)!davidn

  INTEGER, INTENT(in) :: kbdim, kproma
  INTEGER, INTENT(in) :: kk
  INTEGER, INTENT(in) :: klev

  REAL(dp), INTENT(in) :: paclc(kbdim)   !< Cloud cover
  REAL(dp), INTENT(in) :: pdp(kbdim)     !< pressure difference of the layer [Pa]
  REAL(dp), INTENT(in) :: pevp(kbdim)    !< evaporation of rain [kg/kg]
  REAL(dp), INTENT(in) :: plsdcp(kbdim)  !< latent heat of sublimation div. by the specific heat at constant pres.
  REAL(dp), INTENT(in) :: plvdcp(kbdim)  !< latent heat of vaporization div. by the specific heat at constant pres.
  REAL(dp), INTENT(in) :: prpr(kbdim)    !< rain formation rate [kg/kg]
  REAL(dp), INTENT(in) :: psacl(kbdim)   !< accretion of snow flakes with cloud droplets [kg/kg] !UP: i.e. riming
  REAL(dp), INTENT(in) :: pspr(kbdim)    !< snow formation rate [kg/kg]
  REAL(dp), INTENT(in) :: psub(kbdim)    !< sublimation of snow [kg/kg]
  REAL(dp), INTENT(in) :: ptp1tmp(kbdim) !< temporary value of the updated temperature (t) [K]
  REAL(dp), INTENT(in) :: pxiflux(kbdim) !< flux of ice crystals falling into the grid box from above

  REAL(dp), INTENT(inout) :: pclcpre(kbdim) !< fraction of grid box covered by precip
  REAL(dp), INTENT(inout) :: prfl(kbdim)    !< rain flux [kg/m2/s]
  REAL(dp), INTENT(inout) :: psfl(kbdim)    !< snow flux [kg/m2/s]
  REAL(dp), INTENT(inout) :: psmlt(kbdim)   !< melting of snow [kg/kg]

  REAL(dp), INTENT(out) :: pfevapr(kbdim) !< Evaporation of rain [kg/m2/s]
  REAL(dp), INTENT(out) :: pfrain(kbdim)  !< Rain flux before evaporation [kg/m2/s]
  REAL(dp), INTENT(out) :: pfsnow(kbdim)  !< Snow flux before sublimation [kg/m2/s]
  REAL(dp), INTENT(out) :: pfsubls(kbdim) !< Sublimation of snow [kg/m2/s]
  REAL(dp), INTENT(out) :: ptmp2(kbdim) !davidn

  !local vars:
  REAL(dp) :: zpredel(kbdim) !< Precipitation (sum of rain and snow) formed in this level
  REAL(dp) :: zpretot(kbdim) !< Precipitation (sum of rain and snow) from above
  REAL(dp) :: zzdrr(kbdim)   !< Rain flux [kg/m2/s]
  REAL(dp) :: zzdrs(kbdim)   !< Snow flux [kg/m2/s]

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim), ztmp4(kbdim)
  LOGICAL  :: ll1(kbdim)

  zzdrr(1:kproma) = zcons2*pdp(1:kproma)*prpr(1:kproma)
  zzdrs(1:kproma) = zcons2*pdp(1:kproma)*(pspr(1:kproma)+psacl(1:kproma))

  IF (kk .EQ. klev) THEN
     zzdrs(1:kproma) = zzdrs(1:kproma)+pxiflux(1:kproma)
     ztmp1(1:kproma) = zcons2 * pdp(1:kproma) / (plsdcp(1:kproma)-plvdcp(1:kproma)) &
                     * MAX(0._dp, (ptp1tmp(1:kproma)-tmelt)) !SF zcons
     ztmp2(1:kproma) = MIN(xsec*zzdrs(1:kproma), ztmp1(1:kproma)) !SF zsnmlt
     ptmp2(1:kproma) = ztmp2(1:kproma)!davidn
     zzdrr(1:kproma) = zzdrr(1:kproma) + ztmp2(1:kproma)
     zzdrs(1:kproma) = zzdrs(1:kproma) - ztmp2(1:kproma)
     psmlt(1:kproma) = psmlt(1:kproma) + ztmp2(1:kproma) / (zcons2*pdp(1:kproma))
  END IF

  zpretot(1:kproma)  = prfl(1:kproma)  + psfl(1:kproma)
  zpredel(1:kproma)  = zzdrr(1:kproma) + zzdrs(1:kproma)

  pclcpre(:) = gridbox_frac_falling_hydrometeor(kbdim, kproma, &
                         zpretot(:), pclcpre(:), &
                         zpredel(:), paclc(:))

!--- for in-cloud scavenging
  ll1(1:kproma) = (pclcpre(1:kproma) > epsec)

  ztmp1(1:kproma) = ( prfl(1:kproma) + zzdrr(1:kproma)) / MAX(pclcpre(1:kproma), epsec)
  ztmp2(1:kproma) = ( psfl(1:kproma) + zzdrs(1:kproma)) / MAX(pclcpre(1:kproma), epsec)

  pfrain(1:kproma)  = MERGE(ztmp1(1:kproma), 0._dp, ll1(1:kproma))
  pfsnow(1:kproma)  = MERGE(ztmp2(1:kproma), 0._dp, ll1(1:kproma))

  ztmp3(1:kproma) = ( zcons2*pdp(1:kproma)*pevp(1:kproma) ) / MAX(pclcpre(1:kproma), epsec)
  ztmp4(1:kproma) = ( zcons2*pdp(1:kproma)*psub(1:kproma) ) / MAX(pclcpre(1:kproma), epsec)

  pfevapr(1:kproma) = MERGE(ztmp3(1:kproma), 0._dp, ll1(1:kproma))
  pfsubls(1:kproma) = MERGE(ztmp4(1:kproma), 0._dp, ll1(1:kproma))

  prfl(1:kproma) = prfl(1:kproma) + zzdrr(1:kproma) - zcons2 * pdp(1:kproma) * pevp(1:kproma)
  psfl(1:kproma) = psfl(1:kproma) + zzdrs(1:kproma) - zcons2 * pdp(1:kproma) * psub(1:kproma)

END SUBROUTINE update_precip_fluxes

SUBROUTINE update_tendencies_and_important_vars( &
              !-- IN
              kbdim, kproma, &
              picnc, pcdnc, pxim1, pxlm1, pxtm1_cdnc, pxtm1_icnc, &
              pcnd, pdep, pevp, pfrl, &
              pgenti, pgentl, pimlt, plsdcp, plvdcp, prho, prho_rcp, &
              prpr, psacl, pspr, pxievap, &
              pximlt, pxitec, pxlevap, pxltec, pxisub, psub, psmlt, &
              pxib, pxlb, ptp1tmp, ld_liqcl, ld_icecl, &
              !-- INOUT
              paclc, pqte, ptte, pxite, pxlte, &
              pxtte_cdnc, pxtte_icnc, pmlwc, pmiwc, &
              !-- OUT
              preffl, preffi, pcorrcd8, pcorrxl, pcorric6, pcorrxi)!davidn

  INTEGER, INTENT(in) :: kbdim, kproma

  REAL(dp), INTENT(in) :: picnc(kbdim)      !< Ice crystal number concentration (ICNC) [1/m3]
  REAL(dp), INTENT(in) :: pcdnc(kbdim)      !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp), INTENT(in) :: pxim1(kbdim)      !< cloud ice (t-1)
  REAL(dp), INTENT(in) :: pxlm1(kbdim)      !< cloud liquid water (t-1)
  REAL(dp), INTENT(in) :: pxtm1_cdnc(kbdim) !< tracer mmr for cdnc (t-1)
  REAL(dp), INTENT(in) :: pxtm1_icnc(kbdim) !< tracer mmr for icnc (t-1)
  REAL(dp), INTENT(in) :: pcnd(kbdim)       !< condensation rate [kg/kg]
  REAL(dp), INTENT(in) :: pdep(kbdim)       !< deposition rate [kg/kg]
  REAL(dp), INTENT(in) :: pevp(kbdim)       !< evaporation of rain [kg/kg]
  REAL(dp), INTENT(in) :: pfrl(kbdim)       !< freezing rate [kg/kg]
  REAL(dp), INTENT(in) :: pgenti(kbdim)     !< Variable related to Tompkins cloud cover scheme
  REAL(dp), INTENT(in) :: pgentl(kbdim)     !< Variable related to Tompkins cloud cover scheme
  REAL(dp), INTENT(in) :: pimlt(kbdim)      !< melting of ice if T>273 K [kg/kg]
  REAL(dp), INTENT(in) :: plsdcp(kbdim)     !< latent heat of sublimation div. by the specific heat at constant pres.
  REAL(dp), INTENT(in) :: plvdcp(kbdim)     !< latent heat of vaporization div. by the specific heat at constant pres.
  REAL(dp), INTENT(in) :: prho(kbdim)       !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: prho_rcp(kbdim)   !< Inverse air density
  REAL(dp), INTENT(in) :: prpr(kbdim)       !< rain formation rate [kg/kg]
  REAL(dp), INTENT(in) :: psacl(kbdim)      !< accretion of snow flakes with cloud droplets [kg/kg] !UP: i.e. riming
  REAL(dp), INTENT(in) :: pspr(kbdim)       !< snow formation rate [kg/kg]
  REAL(dp), INTENT(in) :: pxievap(kbdim)    !< evaporation of cloud ice [kg/kg]
  REAL(dp), INTENT(in) :: pximlt(kbdim)     !< melting of the ice falling from above [kg/kg]
  REAL(dp), INTENT(in) :: pxitec(kbdim)     !< tendency of cloud ice from detrainment [kg/kg/s]
  REAL(dp), INTENT(in) :: pxlevap(kbdim)    !< evaporation of cloud water [kg/kg]
  REAL(dp), INTENT(in) :: pxltec(kbdim)     !< tendency of cloud liquid water from detrainment [kg/kg/s]
  REAL(dp), INTENT(in) :: pxisub(kbdim)     !< sublimation of cloud ice [kg/kg]
  REAL(dp), INTENT(in) :: psub(kbdim)       !< sublimation of snow [kg/kg]
  REAL(dp), INTENT(in) :: psmlt(kbdim)      !< melting of snow [kg/kg]
  REAL(dp), INTENT(in) :: pxib(kbdim)       !< Cloud ice in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(in) :: pxlb(kbdim)       !< Cloud liquid water in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(in) :: ptp1tmp(kbdim)    !< updated temperature (t) [K]

  LOGICAL, INTENT(in) :: ld_liqcl(kbdim) !< switch that traces the presence of liquid cloud
  LOGICAL, INTENT(in) :: ld_icecl(kbdim) !< switch that traces the presence of ice cloud
  
  REAL(dp), INTENT(inout) :: paclc(kbdim)      !< cloud cover
  REAL(dp), INTENT(inout) :: pqte(kbdim)       !< tendency of specific humidity
  REAL(dp), INTENT(inout) :: ptte(kbdim)       !< tendency of temperature
  REAL(dp), INTENT(inout) :: pxite(kbdim)      !< tendency of cloud ice
  REAL(dp), INTENT(inout) :: pxlte(kbdim)      !< tendency of cloud liquid water
  REAL(dp), INTENT(inout) :: pxtte_cdnc(kbdim) !< tracer tendency for cdnc
  REAL(dp), INTENT(inout) :: pxtte_icnc(kbdim) !< tracer tendency for icnc
  REAL(dp), INTENT(inout) :: pmlwc(kbdim)      !< in-cloud liq water mmr before rain formation [kg/kg]
  REAL(dp), INTENT(inout) :: pmiwc(kbdim)      !< in-cloud ice water mmr before snow formation [kg/kg]

  REAL(dp), INTENT(out) :: preffl(kbdim) !< ice crystal effectiv radius [um]
  REAL(dp), INTENT(out) :: preffi(kbdim) !< cloud drop effectiv radius [um]
  REAL(dp), INTENT(out) :: pcorrcd8(kbdim),pcorrxl(kbdim),pcorric6(kbdim),pcorrxi(kbdim) !davidn

  !local vars:
  REAL(dp) :: zdxicor(kbdim) !< Avoids negative cloud ice
  REAL(dp) :: zdxlcor(kbdim) !< Avoids negative cloud liquid water
  REAL(dp) :: zxip1(kbdim)   !< Cloud ice mass mixing ratio (t) [kg/kg]
  REAL(dp) :: zxlp1(kbdim)   !< Cloud water mass mixing ratio (t+1) [kg/kg]


  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim)
  LOGICAL  :: ll1(kbdim), ll2(kbdim)

  !>>SF initialization of intent(out) variables compulsory over the whole array
  !     they are passed back to stream elements, which would cause
  !     a crash by printing out some undefined parts of the stream element array otherwise
  preffi(1:kbdim) = 0._dp
  preffl(1:kbdim) = 0._dp
  !<<SF

  pqte(1:kproma)  = pqte(1:kproma) &
                  + ztmst_rcp * ( -pcnd(1:kproma)-pgentl(1:kproma)+pevp(1:kproma)+pxlevap(1:kproma)    &
                                -  pdep(1:kproma)-pgenti(1:kproma)+psub(1:kproma)+pxievap(1:kproma)    &
                                +  pxisub(1:kproma) )

  ptte(1:kproma)  = ptte(1:kproma)     &
                  + ztmst_rcp * ( plvdcp(1:kproma) * (pcnd(1:kproma)+pgentl(1:kproma)   &
                                                     -pevp(1:kproma)-pxlevap(1:kproma)) &
                                + plsdcp(1:kproma) * (pdep(1:kproma)+pgenti(1:kproma)   &
                                                     -psub(1:kproma)-pxievap(1:kproma)-pxisub(1:kproma)) &  
                                + (plsdcp(1:kproma)-plvdcp(1:kproma)) &
                                                   * (-psmlt(1:kproma)-pimlt(1:kproma)  &
                                                     -pximlt(1:kproma)+pfrl(1:kproma)+psacl(1:kproma)) )

  ztmp1(1:kproma) = pxltec(1:kproma) + pxlte(1:kproma)

  ztmp2(1:kproma) = pimlt(1:kproma) + pximlt(1:kproma) - pfrl(1:kproma)   - prpr(1:kproma) &
                  - psacl(1:kproma) + pcnd(1:kproma)   + pgentl(1:kproma) - pxlevap(1:kproma)

  zxlp1(1:kproma) = pxlm1(1:kproma) + ztmst * ztmp1(1:kproma) + ztmp2(1:kproma)

  pxlte(1:kproma) = ztmp1(1:kproma) + ztmst_rcp * ztmp2(1:kproma)

  ztmp1(1:kproma) = pxitec(1:kproma) + pxite(1:kproma)
   
  ztmp2(1:kproma) = pfrl(1:kproma) - pspr(1:kproma) + pdep(1:kproma) + pgenti(1:kproma) - pxievap(1:kproma)

  zxip1(1:kproma) = pxim1(1:kproma) + ztmst * ztmp1(1:kproma) + ztmp2(1:kproma) 

  pxite(1:kproma) = ztmp1(1:kproma) + ztmst_rcp * ztmp2(1:kproma)

!
!--- Included for prognostic CDNC/IC scheme ----------------------------

  !--- Calculate new total tendency of CDNC:
  pxtte_cdnc(1:kproma) = ztmst_rcp * (pcdnc(1:kproma)*prho_rcp(1:kproma)-pxtm1_cdnc(1:kproma))

  !--- Calculate new total tendency of ICNC:
  pxtte_icnc(1:kproma) = ztmst_rcp * (picnc(1:kproma)*prho_rcp(1:kproma)-pxtm1_icnc(1:kproma))

! Corrections: Avoid negative cloud water/ice
!
  ztmp1(1:kproma) = zxlp1(1:kproma)      !SF zxlold

  ll1(1:kproma) = (zxlp1(1:kproma) < ccwmin)

!>>DN: new diags
!UP comment: I think this can remain as pcdnc and should not be changed to _upd
!UP comment: remove pcdnc where cc too small
  pcorrcd8(1:kproma) = MERGE(-pcdnc(1:kproma), 0._dp, ll1(1:kproma))
  pcorrxl(1:kproma)  = MERGE(-zxlp1(1:kproma), 0._dp, ll1(1:kproma))
!<<DN: new diags        

  zdxlcor(1:kproma) = -ztmst_rcp * zxlp1(1:kproma)
  zdxlcor(1:kproma) = MERGE(zdxlcor(1:kproma), 0._dp, ll1(1:kproma))

  pxlte(1:kproma) = pxlte(1:kproma) + zdxlcor(1:kproma)

  ztmp1(1:kproma)      = pxtte_cdnc(1:kproma) - ztmst_rcp*pcdnc(1:kproma)*prho_rcp(1:kproma)
  pxtte_cdnc(1:kproma) = MERGE(ztmp1(1:kproma), pxtte_cdnc(1:kproma), ll1(1:kproma))

  ll2(1:kproma) = (zxip1(1:kproma) < ccwmin)

!>>DN: new diags
!UP comment: I think this can remain as picnc and should not be changed to _upd
  pcorric6(1:kproma) = MERGE(-picnc(1:kproma), 0._dp, ll2(1:kproma))
  pcorrxi(1:kproma)  = MERGE(-zxip1(1:kproma), 0._dp, ll2(1:kproma))
!<<DN: new diags        

  zdxicor(1:kproma)   = -ztmst_rcp * zxip1(1:kproma) 
  zdxicor(1:kproma)   = MERGE(zdxicor(1:kproma), 0._dp, ll2(1:kproma)) 

  pxite(1:kproma) = pxite(1:kproma) + zdxicor(1:kproma)

  ztmp1(1:kproma)      = pxtte_icnc(1:kproma) - ztmst_rcp*picnc(1:kproma)*prho_rcp(1:kproma)
  pxtte_icnc(1:kproma) = MERGE(ztmp1(1:kproma), pxtte_icnc(1:kproma), ll2(1:kproma))

  paclc(1:kproma) = MERGE(0.0_dp,paclc(1:kproma),ll1(1:kproma) .AND. ll2(1:kproma))

  ll1(1:kproma) = (paclc(1:kproma) < clc_min)
  paclc(1:kproma) = MERGE(0._dp, paclc(1:kproma), ll1(1:kproma))

  ll2(1:kproma) = ll1(1:kproma) .OR.            &
                 (pmlwc(1:kproma) < 1.e-20_dp) !SFnote: why this value and not epsec for example??
  pmlwc(1:kproma) = MERGE(0._dp, pmlwc(1:kproma), ll2(1:kproma))

  ll2(1:kproma) = ll1(1:kproma) .OR.            &
                 (pmiwc(1:kproma) < 1.e-20_dp) !SFnote: why this value and not epsec for example??
  pmiwc(1:kproma) = MERGE(0._dp, pmiwc(1:kproma), ll2(1:kproma))

  pqte(1:kproma)    = pqte(1:kproma) - zdxlcor(1:kproma) - zdxicor(1:kproma)
  ptte(1:kproma)    = ptte(1:kproma) + plvdcp(1:kproma)*zdxlcor(1:kproma) &
                                     + plsdcp(1:kproma)*zdxicor(1:kproma)

  !--- In-cloud effective radius [um]:
  ztmp1(:) = breadth_factor(kbdim, kproma, pcdnc(:))

  ztmp1(1:kproma) = 1.E6_dp * ztmp1(1:kproma) &   
                  * ( (3._dp/(4._dp*pi*rhoh2o)) * pxlb(1:kproma) &
                      * prho(1:kproma) / pcdnc(1:kproma) )**(1._dp/3._dp)

  preffl(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_liqcl(1:kproma)) !SF #217

  !--- Ice crystal radii (mixed-phase and cirrus) [um]:
  !SF #176 now this is done consistenly between cl. micro and radiation, with a general formula
  !SF conversion of in-cloud ice mmr from kg/kg to in-cloud g m-3:
  ztmp1(1:kproma) = 1000._dp*pxib(1:kproma)*prho(1:kproma)
  ztmp1(:) = eff_ice_crystal_radius(kbdim, kproma, ztmp1(:), picnc(:))

  !>>SF Correction of cirrus-only ice crystal radii in case of alternate cirrus scheme (#448)
  IF (nic_cirrus == 1) THEN
     ll1(1:kproma)   = (ptp1tmp(1:kproma) < cthomi)
     ztmp2(1:kproma) = 83.8_dp*(1.e3_dp*MAX(pxib(1:kproma),eps)*prho(1:kproma))**0.216_dp  !SF (see #143)

     ztmp1(1:kproma) = MERGE(ztmp2(1:kproma), ztmp1(1:kproma), ll1(1:kproma))
  ENDIF
  !<<SF #448

  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), ceffmin)
  ztmp1(1:kproma) = MIN(ztmp1(1:kproma), ceffmax)   !zrieff

  preffi(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ld_icecl(1:kproma))

END SUBROUTINE update_tendencies_and_important_vars

SUBROUTINE diagnostics( &
              !--IN
              kbdim, kproma, &
              ktop, kk, &
              pcdnc, picnc, paclc, pdpg, pdz, pfrln, prho, prprn, psacln, pxib, pxlb, ptp1tmp, &
              preffl, preffi, ld_liqcl, ld_icecl, &
              !-- INOUT
              pcdnc_ave, pcdnc_ave_acc, pcdnc_ave_burd, pcdnc_ct, pcliwc_time, pcloud_time, &
              picnc_ave, picnc_ave_acc, picnc_ave_burd, piwc_acc, piwp_tovs, plwc_acc, pqacc, pqaut, &
              pqfre, preffi_acc, preffi_time, preffi_tovs, preffl_acc, preffl_ct, preffl_time, &
              pcdnc_burden, picnc_burden, ptau1i, preffct, paclcac )

  INTEGER, INTENT(in) :: kbdim, kproma
  INTEGER, INTENT(in) :: ktop(kbdim) !< !< flag for cloud tops
  INTEGER, INTENT(in) :: kk !< current level

  REAL(dp), INTENT(in) :: pcdnc(kbdim)   !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp), INTENT(in) :: picnc(kbdim)   !< Ice crystal number concentration (CDNC) [1/m3]
  REAL(dp), INTENT(in) :: paclc(kbdim)   !< Cloud cover
  REAL(dp), INTENT(in) :: pdpg(kbdim)    !< delta p over grav [kg/m2]
  REAL(dp), INTENT(in) :: pdz(kbdim)     !< layer thickness [m]
  REAL(dp), INTENT(in) :: pfrln(kbdim)   !< Freezing rate for number conc. [1/m3]
  REAL(dp), INTENT(in) :: prho(kbdim)    !< Air density [kg/m3]
  REAL(dp), INTENT(in) :: prprn(kbdim)   !< Rain formation rate for number conc. [1/m3]
  REAL(dp), INTENT(in) :: psacln(kbdim)  !< Accretion rate for number conc. [1/m3]
  REAL(dp), INTENT(in) :: pxib(kbdim)    !< Cloud ice in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(in) :: pxlb(kbdim)    !< Cloud liquid water in the cloudy part of the grid box [kg/kg]
  REAL(dp), INTENT(in) :: ptp1tmp(kbdim) !< updated temperature (t) [K]
  REAL(dp), INTENT(in) :: preffl(kbdim)  !< ice crystal effectiv radius [um]
  REAL(dp), INTENT(in) :: preffi(kbdim)  !< cloud drop effectiv radius [um]

  LOGICAL, INTENT(in) :: ld_liqcl(kbdim) !< switch that traces the presence of liquid cloud
  LOGICAL, INTENT(in) :: ld_icecl(kbdim) !< switch that traces the presence of ice cloud


  REAL(dp), INTENT(inout) :: pcdnc_ave(kbdim)      !< CDNC averaged over cloudy and cloud-free periods [1/m3]
  REAL(dp), INTENT(inout) :: pcdnc_ave_acc(kbdim)  !< as pcdnc_ave_acc, but accumulated [1/m3]
  REAL(dp), INTENT(inout) :: pcdnc_ave_burd(kbdim) !< CDNC burden [1/m2]
  REAL(dp), INTENT(inout) :: pcdnc_ct(kbdim)       !< cloud top cloud droplet number conc. [1/cm3]
  REAL(dp), INTENT(inout) :: pcliwc_time(kbdim)    !< ice cloud occurence time fraction (accumulated) [1]
  REAL(dp), INTENT(inout) :: pcloud_time(kbdim)    !< liquid cloud occurence time fraction (accumulated) [1]
  REAL(dp), INTENT(inout) :: picnc_ave(kbdim)      !< ICNC averaged over cloudy and cloud-free periods [1/m3]
  REAL(dp), INTENT(inout) :: picnc_ave_acc(kbdim)  !< as picnc_ave_acc, but accumulated [1/m3]
  REAL(dp), INTENT(inout) :: picnc_ave_burd(kbdim) !< ICNC burden [1/m2]
  REAL(dp), INTENT(inout) :: piwc_acc(kbdim)       !< ice wat cont acc.+ cloud weighted [kg m-3]
  REAL(dp), INTENT(inout) :: piwp_tovs(kbdim)      !< IWP sampled a la TOVS [kg m-2]
  REAL(dp), INTENT(inout) :: plwc_acc(kbdim)       !< liq wat cont acc.+ cloud weighted
  REAL(dp), INTENT(inout) :: pqacc(kbdim)          !< CD accretion rate [m-3 s-1]
  REAL(dp), INTENT(inout) :: pqaut(kbdim)          !< CD autoconversion rate [m-3 s-1]
  REAL(dp), INTENT(inout) :: pqfre(kbdim)          !< CD freezing rate [m-3 s-1]
  REAL(dp), INTENT(inout) :: preffi_acc(kbdim)     !< ice crystal effectiv radius weighted [um]
  REAL(dp), INTENT(inout) :: preffi_time(kbdim)    !< accumulted semi-transp. cirrus time [1]
  REAL(dp), INTENT(inout) :: preffi_tovs(kbdim)    !< semi-transparent cirrus effectiv radius [um]
  REAL(dp), INTENT(inout) :: preffl_acc(kbdim)     !< cloud drop effectiv radius weighted [um]
  REAL(dp), INTENT(inout) :: preffl_ct(kbdim)      !< cloud top effectiv radius weighted [um]
  REAL(dp), INTENT(inout) :: preffl_time(kbdim)    !< cloud top effectiv radius occ.time [1]
  REAL(dp), INTENT(inout) :: pcdnc_burden(kbdim)   !< Vertically integrated CDNC [1/m2]
  REAL(dp), INTENT(inout) :: picnc_burden(kbdim)   !< Vertically integrated ICNC [1/m2]
  REAL(dp), INTENT(inout) :: ptau1i(kbdim)         !< ice cloud optical depth - visible wavelength
  REAL(dp), INTENT(inout) :: preffct(kbdim)        !< cloud top effective radius [m]
  REAL(dp), INTENT(inout) :: paclcac(kbdim)        !< cloud cover, accumulated

  !-- Local temporary vars with multiple meanings
  REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim), ztmp3(kbdim), ztmp4(kbdim)
  LOGICAL  :: ll1(kbdim), ll2(kbdim), ll3(kbdim)

  pqaut(1:kproma) = pqaut(1:kproma) - zdt*prprn(1:kproma)
  pqfre(1:kproma) = pqfre(1:kproma) - zdt*pfrln(1:kproma)
  pqacc(1:kproma) = pqacc(1:kproma) - zdt*psacln(1:kproma)

  ztmp1(1:kproma)         = pcdnc_ave_acc(1:kproma) + zdtime*pcdnc(1:kproma)
  pcdnc_ave_acc(1:kproma) = MERGE(ztmp1(1:kproma), pcdnc_ave_acc(1:kproma), ld_liqcl(1:kproma))

  ztmp1(1:kproma)    = plwc_acc(1:kproma) + zdtime*pxlb(1:kproma)*prho(1:kproma)
  plwc_acc(1:kproma) = MERGE(ztmp1(1:kproma), plwc_acc(1:kproma), ld_liqcl(1:kproma))

  ztmp1(1:kproma)       = pcloud_time(1:kproma) + zdtime
  pcloud_time(1:kproma) = MERGE(ztmp1(1:kproma), pcloud_time(1:kproma), ld_liqcl(1:kproma))

  !--- In-cloud cdnc burden:
  ztmp1(1:kproma)        = pcdnc_burden(1:kproma)+pcdnc(1:kproma)*pdz(1:kproma)
  pcdnc_burden(1:kproma) = MERGE(ztmp1(1:kproma), pcdnc_burden(1:kproma), ld_liqcl(1:kproma))

  !---- cdnc and burden averaged over cloudy and cloud-free periods
  ztmp1(1:kproma)     = pcdnc_ave(1:kproma) + zdtime*pcdnc(1:kproma)*paclc(1:kproma)
  pcdnc_ave(1:kproma) = MERGE(ztmp1(1:kproma), pcdnc_ave(1:kproma), ld_liqcl(1:kproma))

  ztmp1(1:kproma)          = pcdnc_ave_burd(1:kproma) + zdtime*pcdnc(1:kproma)*pdz(1:kproma)*paclc(1:kproma)
  pcdnc_ave_burd(1:kproma) = MERGE(ztmp1(1:kproma), pcdnc_ave_burd(1:kproma), ld_liqcl(1:kproma))

  !--- Accumulated in-cloud effective radius [um]:
  preffl_acc(1:kproma) = preffl_acc(1:kproma) + zdtime*preffl(1:kproma)

  ll1(1:kproma) = ld_liqcl(1:kproma) &
            .AND. (ktop(1:kproma)  == kk ) &
            .AND. (ptp1tmp(1:kproma)   >  tmelt) &
            .AND. (preffct(1:kproma)   <  4._dp) &
            .AND. (preffl(1:kproma) >= 4._dp)

  ztmp1(1:kproma)     = preffl_ct(1:kproma) + zdtime*preffl(1:kproma)
  preffl_ct(1:kproma) = MERGE(ztmp1(1:kproma), preffl_ct(1:kproma), ll1(1:kproma))

  ztmp1(1:kproma)     = pcdnc_ct(1:kproma) + zdtime*pcdnc(1:kproma)*paclc(1:kproma)
  pcdnc_ct(1:kproma)  = MERGE(ztmp1(1:kproma), pcdnc_ct(1:kproma), ll1(1:kproma))

  ztmp1(1:kproma)       = preffl_time(1:kproma) + zdtime
  preffl_time(1:kproma) = MERGE(ztmp1(1:kproma), preffl_time(1:kproma), ll1(1:kproma))

  preffct(1:kproma) = MERGE(preffl(1:kproma), preffct(1:kproma), ll1(1:kproma))

  ztmp1(1:kproma)         = picnc_ave_acc(1:kproma) + zdtime*picnc(1:kproma)
  picnc_ave_acc(1:kproma) = MERGE(ztmp1(1:kproma), picnc_ave_acc(1:kproma), ld_icecl(1:kproma))

  ztmp1(1:kproma)    = piwc_acc(1:kproma) + zdtime*pxib(1:kproma)*prho(1:kproma)
  piwc_acc(1:kproma) = MERGE(ztmp1(1:kproma), piwc_acc(1:kproma), ld_icecl(1:kproma))

  preffi_acc(1:kproma) = preffi_acc(1:kproma) + zdtime*preffi(1:kproma)

  ztmp3(1:kproma)       = pcliwc_time(1:kproma) + zdtime
  pcliwc_time(1:kproma) = MERGE(ztmp3(1:kproma), pcliwc_time(1:kproma), ld_icecl(1:kproma)) 

  !SF tovs diagnostics:
  ll2(1:kproma) = ld_icecl(1:kproma) .AND. .NOT. ll1(1:kproma)

  ztmp3(1:kproma)  = 1000._dp*pxib(1:kproma)*paclc(1:kproma)*pdpg(1:kproma)  ! iwp in g/m2
  ztmp4(1:kproma)  = ptau1i(1:kproma) + 1.9787_dp*ztmp3(1:kproma)*MAX(preffi(1:kproma),ceffmin)**(-1.0365_dp)
  ptau1i(1:kproma) = MERGE(ztmp4(1:kproma), ptau1i(1:kproma), ll2(1:kproma))
  
  ll3(1:kproma) = ll2(1:kproma) &
            .AND. (ptau1i(1:kproma) > 0.7_dp) &
            .AND. (ptau1i(1:kproma) < 3.8_dp)

  ztmp2(1:kproma)       = preffi_tovs(1:kproma) + zdtime*preffi(1:kproma)
  preffi_tovs(1:kproma) = MERGE(ztmp2(1:kproma), preffi_tovs(1:kproma), ll3(1:kproma)) 
              
  ztmp2(1:kproma)       = preffi_time(1:kproma) + zdtime
  preffi_time(1:kproma) = MERGE(ztmp2(1:kproma), preffi_time(1:kproma), ll3(1:kproma)) 

  ztmp2(1:kproma)       = piwp_tovs(1:kproma) + zdtime*ztmp3(1:kproma)
  piwp_tovs(1:kproma)   = MERGE(ztmp2(1:kproma), piwp_tovs(1:kproma), ll3(1:kproma))
  !SF end tovs diagnostics 

  !--- In-cloud picnc_ave burden:
  ztmp1(1:kproma)        = picnc_burden(1:kproma)+picnc(1:kproma)*pdz(1:kproma)
  picnc_burden(1:kproma) = MERGE(ztmp1(1:kproma), picnc_burden(1:kproma), ld_icecl(1:kproma)) 

  !---- picnc_ave and burden averaged over cloudy and cloud-free periods
  ztmp1(1:kproma)     = picnc_ave(1:kproma) + zdtime*picnc(1:kproma)*paclc(1:kproma)
  picnc_ave(1:kproma) = MERGE(ztmp1(1:kproma), picnc_ave(1:kproma), ld_icecl(1:kproma))

  ztmp1(1:kproma)          = picnc_ave_burd(1:kproma) + zdtime*picnc(1:kproma)*pdz(1:kproma)*paclc(1:kproma)
  picnc_ave_burd(1:kproma) = MERGE(ztmp1(1:kproma), picnc_ave_burd(1:kproma), ld_icecl(1:kproma))  
  !--- End included for pcdnc_ave/IC scheme -----------------------------------

  paclcac(1:kproma) = paclcac(1:kproma) + zdtime*paclc(1:kproma)

END SUBROUTINE diagnostics

FUNCTION threshold_vert_vel_1d(kbdim, kproma, pesw, pesi, picnc, price, peta) RESULT(pvervmax)

  !-- SF: this function computes the threshold vertical velocity needed for the WBF criterion

  !-- Function arguments
  INTEGER :: kbdim, kproma

  REAL(dp) :: pesw(kbdim)     !< Saturation vapor pressure w.r.t. water [Pa]
  REAL(dp) :: pesi(kbdim)     !< Saturation vapor pressure w.r.t. ice [Pa]
  REAL(dp) :: picnc(kbdim)    !< Ice crystal number concentration (ICNC) [1/m3]
  REAL(dp) :: price(kbdim)    !< Volume mean ice crystal radius [m]
  REAL(dp) :: peta(kbdim)     !< ??
  REAL(dp) :: pvervmax(kbdim) !< Threshold vertical velocity

  pvervmax(1:kproma) = (pesw(1:kproma) - pesi(1:kproma)) / pesi(1:kproma) &
                     * picnc(1:kproma) * price(1:kproma) * peta(1:kproma)

END FUNCTION threshold_vert_vel_1d

FUNCTION threshold_vert_vel_2d(kbdim, kproma, klev, pesw, pesi, picnc, price, peta) RESULT(pvervmax)

  !-- SF: this function computes the threshold vertical velocity needed for the WBF criterion

  !-- Function arguments
  INTEGER :: kbdim, kproma, klev

  REAL(dp) :: pesw(kbdim,klev)     !< Saturation vapor pressure w.r.t. water [Pa]
  REAL(dp) :: pesi(kbdim,klev)     !< Saturation vapor pressure w.r.t. ice [Pa]
  REAL(dp) :: picnc(kbdim,klev)    !< Ice crystal number concentration (ICNC) [1/m3]
  REAL(dp) :: price(kbdim,klev)    !< Volume mean ice crystal radius [m]
  REAL(dp) :: peta(kbdim,klev)     !< ??
  REAL(dp) :: pvervmax(kbdim,klev) !< Threshold vertical velocity

  pvervmax(1:kproma,:) = (pesw(1:kproma,:) - pesi(1:kproma,:)) / pesi(1:kproma,:) &
                       * picnc(1:kproma,:) * price(1:kproma,:) * peta(1:kproma,:)

END FUNCTION threshold_vert_vel_2d

FUNCTION eff_ice_crystal_radius(kbdim, kproma, pxice, picnc) RESULT(prieff)

  !-- SF: this function computes the effective ice crystal radius based on Lohmann et al, ERL 2008, expression (1)
  !       The necessary params fact_PK and pow_PK are taken from Pruppacher & Klett 1997 and stored in 
  !       mo_cloud_utils.f90

  !-- Function arguments
  INTEGER :: kbdim, kproma

  REAL(dp) :: pxice(kbdim)  !< in-cloud ice mmr [g m-3] <--- WARNING!! beware of unit!!!
  REAL(dp) :: picnc(kbdim)  !< ice crystal number concentration (ICNC) [1/m3]
  REAL(dp) :: prieff(kbdim) !< effective ice crystal radius [um]

  prieff(1:kproma) = 0.5e4_dp * ( pxice(1:kproma) / fact_PK / picnc(1:kproma) )**(1._dp/pow_PK)

END FUNCTION eff_ice_crystal_radius

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

FUNCTION gridbox_frac_falling_hydrometeor( &
              kbdim, kproma, &
              pflx_from_above, pfrac_from_above, &
              pflx_from_level, pfrac_from_level) &

              RESULT(pfrac_tot)

  !-- This function computes the grid box fraction covered by falling hydrometeor (e.g. rain+snow, sedimenting ice...)

  !-- Function arguments
  INTEGER :: kbdim, kproma

  REAL(dp) :: pflx_from_above(kbdim)  !< Flux of falling hydrometeor from above
  REAL(dp) :: pfrac_from_above(kbdim) !< Fraction of gridbox covered by falling hydrometeor from above
  REAL(dp) :: pflx_from_level(kbdim)  !< Flux of falling hydrometeor from current level
  REAL(dp) :: pfrac_from_level(kbdim) !< Fraction of gridbox covered by falling hydrometeor from current level
  REAL(dp) :: pfrac_tot(kbdim)        !< Total fraction of gridbox covered by falling hydrometeor

  !-- Local vars
  LOGICAL :: ll1(kbdim)

  REAL(dp) :: zflx_tot(kbdim)
  REAL(dp) :: ztmp1(kbdim)

  ll1(1:kproma) = (pflx_from_above(1:kproma) > pflx_from_level(1:kproma))

  pfrac_from_above(1:kproma) = MERGE(pfrac_from_above(1:kproma), pfrac_from_level(1:kproma), ll1(1:kproma))

  zflx_tot(1:kproma) = pflx_from_above(1:kproma) + pflx_from_level(1:kproma)
  
  ll1(1:kproma) = (zflx_tot(1:kproma) > cqtmin)

  ztmp1(1:kproma) = ( pfrac_from_level(1:kproma) * pflx_from_level(1:kproma) &
                    + pfrac_from_above(1:kproma) * pflx_from_above(1:kproma) ) &
                  / MAX(zflx_tot(1:kproma), cqtmin)
  ztmp1(1:kproma) = MIN(ztmp1(1:kproma), 1.0_dp)
  ztmp1(1:kproma) = MAX(ztmp1(1:kproma), 0.0_dp)

  pfrac_tot(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ll1(1:kproma))

END FUNCTION gridbox_frac_falling_hydrometeor

FUNCTION effective_2_volmean_radius_param_Schuman_2011_1d(kbdim, kproma, prieff) RESULT(prvolmean)

  ! Simple param of r/re approximating the Schumann et al. 2011 data:

  !-- Function arguments
  INTEGER :: kbdim, kproma

  REAL(dp) :: prieff(kbdim)    !< effective ice crystal radius [1.e-6 m] <-- beware of units!
  REAL(dp) :: prvolmean(kbdim) !< volume mean ice crystal radius [m]     <-- beware of units!

  prvolmean(1:kproma) = MAX(1.e-6_dp,conv_effr2mvr*1.e-6_dp*prieff(1:kproma))

END FUNCTION effective_2_volmean_radius_param_Schuman_2011_1d

FUNCTION effective_2_volmean_radius_param_Schuman_2011_2d(kbdim, kproma, klev, prieff) RESULT(prvolmean)

  ! Simple param of r/re approximating the Schumann et al. 2011 data:

  !-- Function arguments
  INTEGER :: kbdim, kproma, klev

  REAL(dp) :: prieff(kbdim,klev)    !< effective ice crystal radius [1.e-6 m] <-- beware of units!
  REAL(dp) :: prvolmean(kbdim,klev) !< volume mean ice crystal radius [m]     <-- beware of units!

  prvolmean(1:kproma,:) = MAX(1.e-6_dp,conv_effr2mvr*1.e-6_dp*prieff(1:kproma,:))

END FUNCTION effective_2_volmean_radius_param_Schuman_2011_2d

FUNCTION breadth_factor(kbdim, kproma, pcdnc) RESULT(pkap)

  ! Breadth factor as a function of cloud droplet number concentration

  !-- Function arguments
  INTEGER :: kbdim, kproma

  REAL(dp) :: pcdnc(kbdim) !< Cloud droplet number concentration (CDNC) [1/m3]
  REAL(dp) :: pkap(kbdim)  !< Breadth factor [m?]

  pkap(1:kproma) = 0.00045e-6_dp*pcdnc(1:kproma) + 1.18_dp ! Peng & Lohmann, GRL 2003,
                                                           ! doi:10.1029/2003GL017192
                                                           ! equation 6

END FUNCTION breadth_factor

FUNCTION consistency_number_to_mass(kbdim, kproma, pthreshold, pmass, pnumber) RESULT(pnumber_physical)

  ! Returns a *physical* state of the number concentration of a given quantity (e.g. flux of ice crystal number)
  ! in the sense that the input number concentration is reset to 0 whenever its mass counterpart is below a 
  ! certain (small) threshold

  !-- Subroutine arguments
  INTEGER :: kbdim, kproma

  REAL(dp) :: pthreshold
  REAL(dp) :: pmass(kbdim)   !< Placeholder for a mass-related quantity
  REAL(dp) :: pnumber(kbdim) !< Placeholder for the corresponding number-related quantity (input)

  REAL(dp) :: pnumber_physical(kbdim) !< pnumber, converted to a physical quantity 

  !-- Local vars
  LOGICAL :: ll1(kbdim)

  ll1(1:kproma) = (pmass(1:kproma) < pthreshold)

  pnumber_physical(1:kproma) = MERGE(0._dp, pnumber(1:kproma), ll1(1:kproma))

END FUNCTION consistency_number_to_mass

FUNCTION minimum_CDNC(kbdim, kproma, pxwat) RESULT(pcdnc_min)

  ! Sets the minimum cloud droplet number concentration, either statically (cdnc_min_fixed), or
  ! dynamically based on a maximum cloud droplet volume radius (rcd_vol_max)
  ! Author: David Neubauer (see #475)
  !
  ! SF: refactoring to allow easy switch between both solution

  USE mo_cloud_utils,        ONLY: cdnc_min_lower, cdnc_min_upper
  USE mo_param_switches,     ONLY: ldyn_cdnc_min, cdnc_min_fixed, &
                                   rcd_vol_max !UP: move from mo_cloud_utils to
                                                ! make this tunable from
                                                ! settings file

  !-- Function arguments
  INTEGER :: kbdim, kproma 

  REAL(dp) :: pxwat(kbdim)     !< in-cloud wat mmr [kg m-3] <--- WARNING!! beware of unit!!!
  REAL(dp) :: pcdnc_min(kbdim) !< minimum cloud droplet number concentration [m-3]

  IF (ldyn_cdnc_min) THEN !SF dynamical value for minimum CDNC

     !DN #475: minimum CDNC by using a maximum droplet volume radius
     pcdnc_min(1:kproma) = rcd_vol_max**(-3._dp)*(3._dp/(4._dp*pi*rhoh2o)) * pxwat(1:kproma)
     pcdnc_min(1:kproma) = MIN( cdnc_min_upper, &
                                MAX(cdnc_min_lower, pcdnc_min(1:kproma)) )

  ELSE !SF static minimum CDNC
     pcdnc_min(1:kproma) = cdnc_min_fixed * 1.e6_dp !SF As a side-effect of #489, I prefered to change unit 
                                                    !   of cdnc_min_fixed to cm-3, in order to be able to
                                                    ! handle it as an integer...
                                                    !UP comment: so you need to change it back to m-3 here.
  ENDIF

END FUNCTION minimum_CDNC

END MODULE mo_cloud_micro_2m
