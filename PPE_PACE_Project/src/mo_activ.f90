MODULE mo_activ

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  USE mo_submodel_diag, ONLY: vmem3d, vmem2d

  IMPLICIT NONE

  PUBLIC activ_initialize
  PUBLIC activ_updraft
  PUBLIC activ_lin_leaitch
  PUBLIC construct_activ_stream
  !>>UP
  PUBLIC ham_activ_diag_lin_leaitch_2
  !<<UP

  PRIVATE

  INTEGER,         PUBLIC :: idt_cdnc, idt_icnc, nfrzmod

  TYPE (t_stream), PUBLIC, POINTER :: activ

  INTEGER,         PUBLIC          :: nw ! actual number of updraft velocity (w) bins 
                                         ! (can be 1 if characteristic updraft is used)

  REAL(dp),        PUBLIC, POINTER :: swat(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: w_large(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: w_turb(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: w_cape(:,:)
  REAL(dp),        PUBLIC, POINTER :: w_sigma(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: reffi(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: na(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qnuc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qaut(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qacc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qfre(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qeva(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qmel(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc_instantan(:,:,:) ! Ice crystal number concentration (ICNC), actual instantaneous value [1/m3]
  REAL(dp),        PUBLIC, POINTER :: lwc_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: iwc_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cloud_time(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cliwc_time(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_burden_acc(:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_burden(:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc_burden_acc(:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc_burden(:,:)
  REAL(dp),        PUBLIC, POINTER :: burden_time(:,:)
  REAL(dp),        PUBLIC, POINTER :: burdic_time(:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: reffi_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cloud_cover_duplic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: sice(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl_ct(:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl_time(:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_ct(:,:)
  REAL(dp),        PUBLIC, POINTER :: reffi_tovs(:,:)
  REAL(dp),        PUBLIC, POINTER :: reffi_time(:,:)
  REAL(dp),        PUBLIC, POINTER :: iwp_tovs(:,:)
!>>SF Kasja diags
  REAL(dp),        PUBLIC, POINTER :: qsprn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qrprn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qnucl(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcnd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qlwc_detr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qevp_lwc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qautn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qracl(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qracln(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsacl(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsacln(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qfrz(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qfrzn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qnuci(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qdep(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qiwc_detr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsub_iwc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qagg(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qaggn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsaci(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsacin(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qselfn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsecprod(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsecprodn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsedi(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsedin(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qmlt(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qmltn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qevp_rain(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsub_snow(:,:,:)
!<<SF Kasja diags
!>>UP #783
  REAL(dp),        PUBLIC, POINTER :: qfrznhet(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qfrznhom(:,:,:)
!<<UP #783
!>>DN: new diags
  REAL(dp),        PUBLIC, POINTER :: qmltn2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qxmlt(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qevabfn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qevabf(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qdepbf(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcdnc_detr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qicnc_detr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric1(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric3(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric4(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric5(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric6(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric7(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd1(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd3(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd4(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd5(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd6(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd7(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd8(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd9(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd10(:,:,:)
!>>UP: new diags for CD activation correction
! and further split of correction terms, #783
  REAL(dp),        PUBLIC, POINTER :: qcorrcd2unphys_2d(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd1_1(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd1_2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd1_3(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd1_4(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd2_1(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrcd2_2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric1_1(:,:,:)
!UP #783.3
  REAL(dp),        PUBLIC, POINTER :: qcorrcd2_3(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric1_2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorric1_3(:,:,:)
!<<UP: new diags
  REAL(dp),        PUBLIC, POINTER :: qcorrxi(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrxi2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrxl(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcorrxl2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qxlte(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qxite(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qxttecdnc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qxtteicnc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qspr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qrpr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qtestCD(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qtestIC(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qtestLWC(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qtestIWC(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qgentl(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qgenti(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qmlt_snow(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsub_ice(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qmlts_atm(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qmlt_conv(:,:,:)
!<<DN: new diags
!>>DN: burden
  REAL(dp),        PUBLIC, POINTER :: daut(:,:)
  REAL(dp),        PUBLIC, POINTER :: dfre(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsprn(:,:)
  REAL(dp),        PUBLIC, POINTER :: drprn(:,:)
  REAL(dp),        PUBLIC, POINTER :: dnucl(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcnd(:,:)
  REAL(dp),        PUBLIC, POINTER :: dlwc_detr(:,:)
  REAL(dp),        PUBLIC, POINTER :: devp_lwc(:,:)
  REAL(dp),        PUBLIC, POINTER :: dautn(:,:)
  REAL(dp),        PUBLIC, POINTER :: dracl(:,:)
  REAL(dp),        PUBLIC, POINTER :: dracln(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsacl(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsacln(:,:)
  REAL(dp),        PUBLIC, POINTER :: dfrz(:,:)
  REAL(dp),        PUBLIC, POINTER :: dfrzn(:,:)
  REAL(dp),        PUBLIC, POINTER :: dnuci(:,:)
  REAL(dp),        PUBLIC, POINTER :: ddep(:,:)
  REAL(dp),        PUBLIC, POINTER :: diwc_detr(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsub_iwc(:,:)
  REAL(dp),        PUBLIC, POINTER :: dagg(:,:)
  REAL(dp),        PUBLIC, POINTER :: daggn(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsaci(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsacin(:,:)
  REAL(dp),        PUBLIC, POINTER :: dselfn(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsecprod(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsecprodn(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsedi(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsedin(:,:)
  REAL(dp),        PUBLIC, POINTER :: dmlt(:,:)
  REAL(dp),        PUBLIC, POINTER :: dmltn(:,:)
  REAL(dp),        PUBLIC, POINTER :: devp_rain(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsub_snow(:,:)
  REAL(dp),        PUBLIC, POINTER :: dmltn2(:,:)
  REAL(dp),        PUBLIC, POINTER :: dxmlt(:,:)
  REAL(dp),        PUBLIC, POINTER :: devabfn(:,:)
  REAL(dp),        PUBLIC, POINTER :: devabf(:,:)
  REAL(dp),        PUBLIC, POINTER :: ddepbf(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcdnc_detr(:,:)
  REAL(dp),        PUBLIC, POINTER :: dicnc_detr(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric1(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric2(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric3(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric4(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric5(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric6(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric7(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd1(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd2(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd3(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd4(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd5(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd6(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd7(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd8(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd9(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd10(:,:)
!>>UP: new diags for CD activation correction, #783
  REAL(dp),        PUBLIC, POINTER :: dcorrcd2unphys_2d(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd1_1(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd1_2(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd1_3(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd1_4(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd2_1(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrcd2_2(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric1_1(:,:)
!UP #783.3
  REAL(dp),        PUBLIC, POINTER :: dcorrcd2_3(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric1_2(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorric1_3(:,:)
!<<UP: new diags
!>>UP #783
  REAL(dp),        PUBLIC, POINTER :: dfrznhet(:,:)
  REAL(dp),        PUBLIC, POINTER :: dfrznhom(:,:)
!<<UP #783
!>>UP: new cc diags #783
  REAL(dp),        PUBLIC, POINTER :: lcc_mod(:,:)
  REAL(dp),        PUBLIC, POINTER :: mcc_mod(:,:)
  REAL(dp),        PUBLIC, POINTER :: icc_mod(:,:)
  REAL(dp),        PUBLIC, POINTER :: diag_cnt_cqtmins1(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: diag_cnt_cqtmins2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: diag_cnt_cqtmins3(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: diag_cnt_cqtminl1(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: diag_cnt_cqtminl2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: diag_cnt_cqtminl3(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: diag_cnt_ll12d_1(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: diag_cnt_ll12d_2(:,:,:)
!<<UP: new cc diags
!>>UP: new diags to evaluate/improve ccnclim
  REAL(dp),        PUBLIC, POINTER :: nact_strat_tmp(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: nact_conv_tmp(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_tmp(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_conv_tmp(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_cdncmin(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_cdncmin_ctr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_actccn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_actccn_ctr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_noactccn(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_noactccn_ctr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_actccnwobase(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_actccnwobase_ctr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_nosmall40(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_nosmall40_ctr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_nosmall1(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_nosmall1_ctr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_so4ks(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_max(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_max_ctr(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: ccn_strat_2d(:,:)
  REAL(dp),        PUBLIC, POINTER :: rwet_KS_diag(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: rwet_KS_diag_cloudbase(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: rwet_KS_diag_cloudbase_ctr(:,:,:)
!<<UP
  REAL(dp),        PUBLIC, POINTER :: dcorrxi(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrxi2(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrxl(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcorrxl2(:,:)
  REAL(dp),        PUBLIC, POINTER :: dxlte(:,:)
  REAL(dp),        PUBLIC, POINTER :: dxite(:,:)
  REAL(dp),        PUBLIC, POINTER :: dxttecdnc(:,:)
  REAL(dp),        PUBLIC, POINTER :: dxtteicnc(:,:)
  REAL(dp),        PUBLIC, POINTER :: dspr(:,:)
  REAL(dp),        PUBLIC, POINTER :: drpr(:,:)
  REAL(dp),        PUBLIC, POINTER :: dtestCD(:,:)
  REAL(dp),        PUBLIC, POINTER :: dtestIC(:,:)
  REAL(dp),        PUBLIC, POINTER :: dtestLWC(:,:)
  REAL(dp),        PUBLIC, POINTER :: dtestIWC(:,:)
  REAL(dp),        PUBLIC, POINTER :: dgentl(:,:)
  REAL(dp),        PUBLIC, POINTER :: dgenti(:,:)
  REAL(dp),        PUBLIC, POINTER :: dmlt_snow(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsub_ice(:,:)
  REAL(dp),        PUBLIC, POINTER :: dmlts_atm(:,:)
  REAL(dp),        PUBLIC, POINTER :: dmlts_sfc(:,:)
  REAL(dp),        PUBLIC, POINTER :: dmlt_conv(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsedi_sfc(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsnow(:,:)
  REAL(dp),        PUBLIC, POINTER :: dcldtte(:,:)
  REAL(dp),        PUBLIC, POINTER :: dconvtte(:,:)
!<<DN: burden
!>>UP #844
  REAL(dp),        PUBLIC, POINTER  :: slfdiag_iwp_mxpT(:,:)
  REAL(dp),        PUBLIC, POINTER  :: slfdiag_lwp_mxpT(:,:)
  REAL(dp),        PUBLIC, POINTER  :: slfdiag_iwp_l35(:,:)
  REAL(dp),        PUBLIC, POINTER  :: slfdiag_lwp_g0(:,:)
  REAL(dp),        PUBLIC, POINTER  :: slfdiag_icnc_mxpT(:,:)
  REAL(dp),        PUBLIC, POINTER  :: slfdiag_cdnc_mxpT(:,:)
  REAL(dp),        PUBLIC, POINTER  :: slfdiag_icnc_l35(:,:)
  REAL(dp),        PUBLIC, POINTER  :: slfdiag_cdnc_g0(:,:)
!<<UP #844
!>>UP lifetime diagnostics
  REAL(dp),        PUBLIC, POINTER  :: q_lwc_sources(:,:,:)
  REAL(dp),        PUBLIC, POINTER  :: q_lwc_sinks(:,:,:)
  REAL(dp),        PUBLIC, POINTER  :: q_iwc_sources(:,:,:)
  REAL(dp),        PUBLIC, POINTER  :: q_iwc_sinks(:,:,:)
  REAL(dp),        PUBLIC, POINTER  :: q_cdnc_sources(:,:,:)
  REAL(dp),        PUBLIC, POINTER  :: q_cdnc_sinks(:,:,:)
  REAL(dp),        PUBLIC, POINTER  :: q_icnc_sources(:,:,:)
  REAL(dp),        PUBLIC, POINTER  :: q_icnc_sinks(:,:,:)
  REAL(dp),        PUBLIC, POINTER  :: d_lwc_lifetime_sources(:,:)
  REAL(dp),        PUBLIC, POINTER  :: d_lwc_lifetime_sinks(:,:)
  REAL(dp),        PUBLIC, POINTER  :: d_iwc_lifetime_sources(:,:)
  REAL(dp),        PUBLIC, POINTER  :: d_iwc_lifetime_sinks(:,:)
  REAL(dp),        PUBLIC, POINTER  :: d_cdnc_lifetime_sources(:,:)
  REAL(dp),        PUBLIC, POINTER  :: d_cdnc_lifetime_sinks(:,:)
  REAL(dp),        PUBLIC, POINTER  :: d_icnc_lifetime_sources(:,:)
  REAL(dp),        PUBLIC, POINTER  :: d_icnc_lifetime_sinks(:,:)
!<<UP

!>>UP: new ccnclim diagnostics
  REAL(dp),        PUBLIC, POINTER  :: cloudbase_lifetime(:,:,:) ! Variable to track
  ! how long a cloud base has been at this point
  REAL(dp), PUBLIC, DIMENSION(7)    :: ccn_bins = (/ 0._dp, 3.E7_dp, 6.E7_dp, 1.E8_dp, 3.E8_dp, 6.E8_dp, 1.E9_dp /)
  CHARACTER(4), PUBLIC, DIMENSION(7) :: ccn_bins_string = (/ '0E00', '3E07', '6E07', '1E08', '3E08', &
                                                             '6E08', '1E09'/)
  REAL(dp), PUBLIC, DIMENSION(8)    :: w_bins = (/ 0._dp, 2.E-1_dp, 4.E-1_dp, 6.E-1_dp, &
                                                    8.E-1_dp, 1.E0_dp, 2.E0_dp, 2.5E1_dp /)
  CHARACTER(3), PUBLIC, DIMENSION(8) :: w_bins_string = (/ '000', '002', '004', '006', &
                                                            '008', '010', '020', '250'/)
  REAL(dp), PUBLIC, DIMENSION(9)   :: cbl_bins = (/ -0.5_dp, 0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp, 4.5_dp, 10.5_dp, 20.5_dp, 50.5_dp /)
  CHARACTER(4), PUBLIC, DIMENSION(9) :: cbl_bins_string = (/ '-005', '0005', '0015', '0025', &
                                                             '0035', '0045', '0105', '0205', '0505' /)
  REAL(dp), PUBLIC, DIMENSION(7)    :: cdnc_bins = (/ 0._dp, 1.E00_dp, 1.E7_dp, 3.E7_dp, 5.E7_dp, 7.E7_dp, 1.E8_dp /)
  CHARACTER(4), PUBLIC, DIMENSION(7) :: cdnc_bins_string = (/ '0E00', '1E00', '1E07', '3E07', '5E07', '7E07', '1E08'/)
  REAL(dp), PUBLIC, DIMENSION(9)    :: lwc_bins = (/ 0._dp, 4.5E-2_dp, 6.E-2_dp, 7.5E-2_dp, 9.E-2_dp, &
                                                      1.05E-1_dp, 1.2E-1_dp, 1.35E-1_dp, 1.5E-1_dp /)
  CHARACTER(6), PUBLIC, DIMENSION(9) :: lwc_bins_string = (/ '000E-3', '045E-3', '060E-3', '075E-3', &
                                                             '090E-3', '105E-3', '120E-3', '135E-3', '150E-3' /)
  REAL(dp), PUBLIC, DIMENSION(8)    :: prcp_bins = (/ 0._dp, 1.E-8_dp, 1.E-5_dp, 2.E-5_dp, 4.E-5_dp, 8.E-5_dp, &
                                                           1.5E-4_dp, 1.E-3_dp /)
  CHARACTER(6), PUBLIC, DIMENSION(8) :: prcp_bins_string = (/ '0-0E00', '1-0E-8', '1-0E-5', '2-0E-5', '4-0E-5', '8-0E-5', &
                                                             '1-5E-4', '1-0E-3'/)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: ccn_strat_binned(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: ccn_conv_binned(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: ccn_strat_w_binned(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: ccn_strat_w_actccn_binned(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: ccn_strat_cbl_binned(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: cdnc_lwc_binned(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: cdnc_binned(:)
  TYPE(vmem2d), PUBLIC, ALLOCATABLE :: prcp_binned(:)
!<<UP

  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: w(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: w_pdf(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: swat_max_strat(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: swat_max_conv(:)

  REAL(dp)            :: w_min = 0.0_dp       ! minimum characteristic w for activation [m s-1]
  REAL(dp), PARAMETER :: w_sigma_min = 0.1_dp ! minimum value of w standard deviation [m s-1]

  !--- Subroutines:

CONTAINS

  SUBROUTINE activ_updraft(kproma,   kbdim,  klev,    krow, &
                           ptkem1,   pwcape, pvervel, prho, &
                           pw,       pwpdf                  )

    ! *activ_updraft* calculates the updraft vertical velocity
    !                 as sum of large scale and turbulent velocities
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford                 2008
    !
    ! References:
    ! -----------
    ! Lohmann et al., ACP, (2008)
    !

    USE mo_physical_constants, ONLY: grav
    !>>SF #345
    USE mo_cloud_utils,        ONLY: fact_tke
    USE mo_param_switches,     ONLY: ncd_activ
    !<<SF #345
    USE mo_param_switches, ONLY: nactivpdf !ZK

    IMPLICIT NONE

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(out) :: pw(kbdim,klev,nw)        ! stratiform updraft velocity bins, large-scale+TKE (>0.0) [m s-1]
    REAL(dp), INTENT(out) :: pwpdf(kbdim,klev,nw)     ! stratiform updraft velocity PDF

    REAL(dp), INTENT(in)  :: prho(kbdim,klev),      & ! air density
                             ptkem1(kbdim,klev),    & ! turbulent kinetic energy
                             pvervel(kbdim,klev),   & ! large scale vertical velocity [Pa s-1]
                             pwcape(kbdim)            ! CAPE contribution to convective vertical velocity [m s-1]

    REAL(dp) :: zwlarge(kbdim, klev), & ! large-scale vertical velocity [m s-1]
                zwturb(kbdim, klev)     ! TKE-derived vertical velocity or st. dev. thereof [m s-1]

    !--- Large scale vertical velocity in SI units:

    zwlarge(1:kproma,:)      = -1._dp* pvervel(1:kproma,:)/(grav*prho(1:kproma,:))
    w_large(1:kproma,:,krow) = zwlarge(1:kproma,:)

    !--- Turbulent vertical velocity:

    w_turb(1:kproma,:,krow)  = fact_tke*SQRT(ptkem1(1:kproma,:))

    !>>SF #345: correction for the TKE prefactor, in case of Lin & Leaitch scheme only
    IF (ncd_activ == 1) THEN ! Lin & Leaitch scheme
       w_turb(1:kproma,:,krow)  = 1.33_dp*SQRT(ptkem1(1:kproma,:))
    ENDIF
    !<<SF #345

    !--- Convective updraft velocity from CAPE:

    w_cape(1:kproma,krow)  = pwcape(1:kproma) !SF although this is no longer used as a contribution to the
                                              ! convective updraft velocity, this is just kept here
                                              ! for recording it into the activ stream

    !--- Total stratiform updraft velocity:

    IF (nactivpdf == 0) THEN
       !--- Turbulent vertical velocity:
       pw(1:kproma,:,1) = MAX(w_min,w_large(1:kproma,:,krow)+w_turb(1:kproma,:,krow))
       w(1)%ptr(1:kproma,:,krow) = pw(1:kproma,:,1)
       ! Only one "bin", with probability of 1. The actual value doesn't
       ! matter so long as it's finite, since it cancels out of the CDNC
       ! calculation.
       pwpdf(1:kproma,:,1) = 1.0_dp
    ELSE
       CALL aero_activ_updraft_sigma(kproma,   kbdim,   klev,    krow, &
                                     ptkem1,  zwturb                   )

       CALL aero_activ_updraft_pdf(kproma,  kbdim,  klev,  krow, &
                                   zwlarge, zwturb, pw,    pwpdf )
    END IF

  END SUBROUTINE activ_updraft

  SUBROUTINE aero_activ_updraft_sigma(kproma,   kbdim,   klev,    krow, &
                                      ptkem1,   pwsigma                 )

    ! *aero_activ_updraft_sigma* calculates the standard deviation of the pdf of uf 
    !                            updraft vertical velocity
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford                 2013
    !
    ! References:
    ! -----------
    ! West et al., ACP, 2013. 
    !

    IMPLICIT NONE

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(in)  :: ptkem1(kbdim,klev)  ! turbulent kinetic energy

    REAL(dp), INTENT(out) :: pwsigma(kbdim,klev) ! st. dev. of vertical velocity

    !--- Large scale vertical velocity in SI units:

    pwsigma(1:kproma,:)      = MAX(w_sigma_min, ((2.0_dp/3.0_dp)*ptkem1(1:kproma,:))**0.5_dp) ! m/s
    w_sigma(1:kproma,:,krow) = pwsigma(1:kproma,:)

  END SUBROUTINE aero_activ_updraft_sigma

  SUBROUTINE aero_activ_updraft_pdf(kproma,  kbdim,   klev, krow, &
                                    pwlarge, pwsigma, pw,   pwpdf )

    ! *aero_activ_updraft_* calculates Gaussian pdf of  
    !                       updraft vertical velocity
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford                 2013
    !
    ! References:
    ! -----------
    ! West et al., ACP, 2013. 
    !

    USE mo_math_constants, ONLY: pi
    USE mo_param_switches, ONLY: nactivpdf

    IMPLICIT NONE

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(in)  :: pwlarge(kbdim,klev), & ! large-scale vertical velocity [m s-1]
                             pwsigma(kbdim,klev)    ! st. dev. of vertical velocity [m s-1]

    REAL(dp), INTENT(out) :: pw(kbdim,klev,nw),   & ! vertical velocity bins [m s-1]
                             pwpdf(kbdim,klev,nw)   ! vettical velocity PDF [s m-1]

    INTEGER               :: jl, jk, jw

    REAL(dp)              :: zw_width(kbdim,klev), &
                             zw_min(kbdim,klev), &
                             zw_max(kbdim,klev)

    zw_min(1:kproma,:)   = 0.0_dp
    zw_max(1:kproma,:)   = 4.0_dp*pwsigma(1:kproma,:)
    zw_width(1:kproma,:) = (zw_max(1:kproma,:) - zw_min(1:kproma,:)) / DBLE(nw)

    DO jw=1, nw
      pw(1:kproma,:,jw) = zw_min(1:kproma,:) + (DBLE(jw) - 0.5_dp) * zw_width(1:kproma,:)

      pwpdf(1:kproma,:,jw) = (1.0_dp / ((2.0_dp*pi)**0.5_dp))                         &
                             * (1.0_dp / pwsigma(1:kproma,:))                         &
                             * EXP( -((pw(1:kproma,:,jw) - pwlarge(1:kproma,:))**2_dp &
                                    / (2.0_dp*pwsigma(1:kproma,:)**2.0_dp)) )
    END DO

    IF (nactivpdf < 0) THEN
       DO jw=1, nw
          w(jw)%ptr(1:kproma,:,krow) = pw(1:kproma,:,jw)
          w_pdf(jw)%ptr(1:kproma,:,krow) = pwpdf(1:kproma,:,jw)
       END DO
    END IF

  END SUBROUTINE aero_activ_updraft_pdf

!-------------------------------------------

  SUBROUTINE activ_lin_leaitch(kproma,  kbdim,    klev,     krow, &
                               pw, pcdncact                       )

    ! *activ_lin_leaitch* calculates the number of activated aerosol 
    !                     particles from the aerosol number concentration
    !SF now independent of HAM, since HAM-specific calculation are computed in mo_ham_activ
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-MET                       2004
    !
    ! Method:
    ! -------
    ! The parameterisation follows the simple empirical relations of 
    ! Lin and Leaitch (1997).
    ! Updraft velocity is parameterized following Lohmann et al. (1999).
    !

    USE mo_kind,       ONLY: dp
    USE mo_conv,       ONLY: na_cv, cdncact_cv

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(IN)  :: pw(kbdim,klev)  ! stratiform updraft velocity, large-scale+TKE (>0.0) [m s-1]
    REAL(dp), INTENT(out) :: pcdncact(kbdim,klev)  ! number of activated particles

    REAL(dp), PARAMETER :: c2=2.3E-10_dp, & ! [m4 s-1]
                           c3=1.27_dp       ! [1]

    INTEGER  :: jl, jk
    REAL(dp) :: zNmax, zeps

    zeps=EPSILON(1.0_dp)

    pcdncact(:,:) = 0._dp
    cdncact_cv(:,:,krow) = 0._dp

    !--- Aerosol activation:

    DO jk=1, klev
       DO jl=1, kproma

          !--- Stratiform clouds:

          ! Activation occurs only in occurrence of supersaturation

          !>>SF note: 
          !     The previous temperature restriction (temp > homogeneous freezing temp)
          !     has been removed because it was preventing to diagnose the number of
          !     dust and BC particules in soluble modes where temp < hom. freezing.
          !     The rationale behind is that diagnosing this allows further
          !     devel to implement concurrent homogeneous vs heterogenous freezing processes
          !     (which is not yet part of this version, though).
          !
          !IMPORTANT: 
          !     This temperature condition removal is completely transparent for the sanity 
          !     of the current code, since relevant temperature ranges are now safely checked
          !     directly in cloud_cdnc_icnc
          !<<SF

          IF(pw(jl,jk)>zeps .AND. na(jl,jk,krow)>zeps) THEN

             !--- Maximum number of activated particles [m-3]:

             zNmax=(na(jl,jk,krow)*pw(jl,jk))/(pw(jl,jk)+c2*na(jl,jk,krow))

             ! Average number of activated particles [m-3]:
             ! zNmax need to be converted to [cm-3] and the
             ! result to be converted back to [m-3].

             pcdncact(jl,jk)=0.1E6_dp*(1.0E-6_dp*zNmax)**c3

          END IF

          !--- Convective clouds:

          IF(pw(jl,jk)>zeps .AND. na_cv(jl,jk,krow)>zeps) THEN

             zNmax=(na_cv(jl,jk,krow)*pw(jl,jk))/(pw(jl,jk)+c2*na_cv(jl,jk,krow))
             cdncact_cv(jl,jk,krow)=0.1E6_dp*(1.0E-6_dp*zNmax)**c3

          ENDIF

       END DO
    END DO

  END SUBROUTINE activ_lin_leaitch

 SUBROUTINE activ_initialize

    USE mo_control,            ONLY: nlev, nn, &
                                     lcouple !SF
    USE mo_exception,          ONLY: message, em_param
    USE mo_submodel,           ONLY: print_value, lham, lhammoz, lccnclim
    USE mo_echam_cloud_params, ONLY: ccsaut, ccraut
    USE mo_param_switches,     ONLY : icover, nauto, &          !++mgs
                                      ncd_activ, nactivpdf, nic_cirrus, lclmi_progn, &
                                      cdnc_min_fixed
    USE mo_tracer,             ONLY: get_tracer
!davidn
    USE mo_param_switches,     ONLY: tun47ccraut,tun47ccsaut,tun47cdncmin
!davidn    

    CHARACTER(len=24)      :: csubmname

    !--- Set number of updraft bins: 

    SELECT CASE(ABS(nactivpdf))
      CASE(0)
        nw = 1
      CASE(1)
        nw = 20
      CASE DEFAULT
        nw = ABS(nactivpdf)
    END SELECT

    IF (nactivpdf <= 0) THEN
      ! These are used either if not using a PDF, or if per-bin
      ! diagnostics are requested.
      ALLOCATE(w(nw))
      ALLOCATE(w_pdf(nw))
      ALLOCATE(swat_max_strat(nw))
      ALLOCATE(swat_max_conv(nw))
    END IF

    !>>UP: new diagnostics for ccnclim
    ALLOCATE(ccn_strat_binned(size(ccn_bins-1)))
    ALLOCATE(ccn_conv_binned(size(ccn_bins-1)))
    ALLOCATE(ccn_strat_w_binned(size(ccn_bins-1)*size(w_bins-1)))
    ALLOCATE(ccn_strat_w_actccn_binned(size(ccn_bins-1)*size(w_bins-1)))
    ALLOCATE(ccn_strat_cbl_binned(size(ccn_bins-1)*size(cbl_bins-1)))
    ALLOCATE(cdnc_lwc_binned(size(cdnc_bins-1)*size(lwc_bins-1)))
    ALLOCATE(prcp_binned(size(prcp_bins-1)))
    ALLOCATE(cdnc_binned(size(cdnc_bins-1)))
    !<<UP

    !
    !-- overwrite values for coupled CDNC/ICNC cloud scheme
    !
    IF (lclmi_progn)  THEN
      IF (nlev == 31) THEN
         IF (nn == 63) THEN
            SELECT CASE (ncd_activ)
               CASE(1) ! LL activtion
                  !SF: updated on 2015.02.25 (David Neubauer / Katty Huang, pure atm run, HAM-M7, LL activation)
                  ccsaut = 1200._dp
                  ccraut = 3.5_dp
               CASE(2) !AR&G activation
                  SELECT CASE(cdnc_min_fixed)
                     CASE(10)
                        !SF: updated on 2017.02.14 (David Neubauer, pure atm run, HAM-M7)
                        ccsaut = 900._dp
                        ccraut = 2.8_dp
                     CASE(40)
                        !SF: updated on 2017.02.14 (David Neubauer, pure atm run, HAM-M7)
                        ccsaut = 900._dp
                        ccraut = 10.6_dp
                        ! ccraut = 30._dp ! UP changed to 15 here for tuning nclmi
                        ! = 0
                  END SELECT
            END SELECT
         ENDIF
      ENDIF

      IF (nlev == 47) THEN
         IF (nn == 63) THEN
            SELECT CASE (ncd_activ)
               CASE(1) ! LL activtion
                  !SF: updated on 2015.02.19 (David Neubauer, pure atm run, HAM-M7, LL activation)
                  ccsaut = 800._dp
                  ccraut = 5._dp
               CASE(2) ! AR&G activtion
                  SELECT CASE(cdnc_min_fixed)
                     CASE(10)
                        !SF: updated on 2017.02.14 (David Neubauer, pure atm run, HAM-M7)
                        ccsaut = 900._dp
                        ccraut = 2.8_dp
                     CASE(40)
                        !SF: updated on 2017.02.14 (David Neubauer, pure atm run, HAM-M7)
                        ccsaut = 900._dp
                        ccraut = 10.6_dp 
                        !ccraut = 30._dp ! UP changed here for tuning
                  END SELECT
            END SELECT
            !davidn tuning
            ccsaut = tun47ccsaut
            ccraut = tun47ccraut
            !davidn tuning
         ENDIF
      ENDIF
    ENDIF

    IF (ncd_activ == 2) THEN
       ! cdnc_min_fixed = 10.0E6_dp ! UP: original line
       cdnc_min_fixed = 10.0_dp ! UP: I think we need to use units of cm-3 here
       !davidn tuning
       cdnc_min_fixed = tun47cdncmin
       !davidn
       !SF: updated on 2016.04.04 (David Neubauer, pure atm run, HAM-M7, AR&G activation)
    ENDIF

!>>SF
    !-- Define the cdnc and icnc tracer index to point to the correct tracer:
    CALL get_tracer('CDNC',idx=idt_cdnc)
    CALL get_tracer('ICNC',idx=idt_icnc)
!<<SF

    !
    !-- Write out new parameters
    !
    IF (ncd_activ>0 .OR. nic_cirrus>0) THEN

      csubmname = 'UNKNOWN'
      IF (lham) csubmname = 'HAM'
      IF (lhammoz) csubmname = 'HAMMOZ'
      IF (lccnclim) csubmname = 'CCNCLIM'

      CALL message('','')
      CALL message('','----------------------------------------------------------')
      CALL message('activ_initialize','Parameter settings for the ECHAM-'//TRIM(csubmname)  &
                   //' cloud microphysics scheme')
      CALL message('','---')
      CALL print_value('              ncd_activ                       = ', ncd_activ)
      CALL print_value('              nic_cirrus                       = ', nic_cirrus)
      CALL message('', ' => Parameter adjustments in mo_activ:', level=em_param)
      CALL print_value('              ccsaut =', ccsaut)
      CALL print_value('              ccraut =', ccraut)
      CALL message('','---')
      CALL message('','----------------------------------------------------------')

    ENDIF
  END SUBROUTINE activ_initialize

  SUBROUTINE construct_activ_stream

    ! *construct_stream_activ* allocates output streams
    !                          for the activation schemes
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-MET                       2004
    !

  USE mo_memory_base,    ONLY: new_stream, add_stream_element, AUTO,  &
                               default_stream_setting, add_stream_reference
  USE mo_filename,       ONLY: trac_filetype
  USE mo_linked_list,    ONLY: HYBRID
  USE mo_param_switches, ONLY: ncd_activ, nactivpdf, nic_cirrus,& !SF
                               lnewdiags, lslf, lccnclimdiags !UP
  USE mo_ham,            ONLY: lccnclim_diag !UP

  IMPLICIT NONE

  INTEGER           :: jw, jccn
  CHARACTER(len=10) :: cbin


  !--- Create new stream:

  CALL new_stream (activ ,'activ',filetype=trac_filetype)


  !--- Add standard fields for post-processing:

  CALL add_stream_reference (activ, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
  CALL add_stream_reference (activ, 'lsp'     ,'sp'    ,lpost=.TRUE.)
  CALL add_stream_reference (activ, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
  CALL add_stream_reference (activ, 'gboxarea','geoloc',lpost=.TRUE.)

  CALL default_stream_setting (activ, lpost     = .TRUE. , &
                                      lrerun    = .TRUE. , &
                                      leveltype = HYBRID , &
                                      table     = 199,     &
                                      code      = AUTO     )
  !--- 1) Cloud Properties:

  CALL add_stream_element (activ,   'SWAT',       swat,                                   &
                           longname='ECHAM supersaturation over water',   units='% [0-1]' )

  IF (ncd_activ==2) THEN

     IF (nactivpdf == 0) THEN
        CALL add_stream_element (activ,   'SWAT_MAX_STRAT', swat_max_strat(1)%ptr, &
                                 longname='maximum supersaturation stratiform', units='% [0-1]' )

        CALL add_stream_element (activ,   'SWAT_MAX_CONV',  swat_max_conv(1)%ptr, &
                                 longname='maximum supersaturation convective', units='% [0-1]' )
     ELSE IF (nactivpdf < 0) THEN
        DO jw=1,nw
           WRITE (cbin, "(I2.2)") jw
           CALL add_stream_element (activ,   'SWAT_MAX_STRAT_'//TRIM(cbin), swat_max_strat(jw)%ptr, &
                                    longname='maximum supersaturation stratiform, vertical velocity bin '//TRIM(cbin), &
                                    units='% [0-1]' )

           CALL add_stream_element (activ,   'SWAT_MAX_CONV_'//TRIM(cbin), swat_max_conv(jw)%ptr, &
                                    longname='maximum supersaturation convective, vertical velocity bin '//TRIM(cbin), &
                                    units='% [0-1]' )
        END DO
     END IF
  ENDIF

  IF (nactivpdf == 0) THEN
     CALL add_stream_element (activ,   'W',          w(1)%ptr, &
                              longname='total vertical velocity for activation',units='m s-1')
  ELSE IF (nactivpdf < 0) THEN
     DO jw=1, nw
       WRITE (cbin, "(I2.2)") jw
       CALL add_stream_element (activ,   'W_'//TRIM(cbin), w(jw)%ptr, &
                                longname='Vertical velocity bin '//TRIM(cbin)//' for activation', &
                                units='m s-1')

       CALL add_stream_element (activ,   'W_PDF_'//TRIM(cbin), w_pdf(jw)%ptr, &
                                longname='Vertical velocity PDF in bin '//TRIM(cbin)//' for activation', &
                                units='s m-1')
     END DO
  END IF

  CALL add_stream_element (activ,   'W_LARGE',    w_large,                                &
                           longname='large scale vertical velocity',      units='m s-1'   )

  IF (nactivpdf == 0) THEN
     CALL add_stream_element (activ, 'W_TURB',     w_turb,                                 &
                              longname='turbulent vertical velocity',      units='m s-1'   )
  ELSE
     CALL add_stream_element (activ, 'W_SIGMA',    w_sigma,                                    &
                              longname='sub-grid st. dev. of vertical velocity', units='m s-1' )
  END IF

  CALL add_stream_element (activ,   'W_CAPE',     w_cape,                                 &
                           longname='convective updraft velocity from CAPE', units='m s-1')

  CALL add_stream_element (activ,   'REFFL',      reffl,                                  &
                           longname='cloud drop effectiv radius',         units='um'      )

  IF (nic_cirrus>0) THEN

  CALL add_stream_element (activ,   'REFFI',      reffi,                                  &
                           longname='ice crystal effectiv radius',        units='um'      )
  END IF

  CALL add_stream_element (activ,   'NA',         na,                                     &
                           longname='aerosol number for activation',      units='m-3'     )

  CALL default_stream_setting (activ, laccu=.TRUE.)

  CALL add_stream_element (activ,   'QNUC',       qnuc,                                   &
                           longname='CD nucleation rate',                 units='m-3 s-1' )

  CALL add_stream_element (activ,   'QAUT',       qaut,                                   &
                           longname='CD autoconversion rate',             units='m-3 s-1' )

  CALL add_stream_element (activ,   'QACC',       qacc,                                   &
                           longname='CD accretion rate',                  units='m-3 s-1' )

  CALL add_stream_element (activ,   'QFRE',       qfre,                                   &
                           longname='CD freezing rate',                   units='m-3 s-1' )
  !>>dod deleted QEVA, not used anywhere
  !  CALL add_stream_element (activ,   'QEVA',       qeva,                                   &
  !                           longname='CD evaporation rate',                units='m-3 s-1' )

  CALL add_stream_element (activ,   'QMEL',       qmel,                                   &
                           longname='CD source rate from melting ice',    units='m-3 s-1' )

  CALL add_stream_element (activ,   'CDNC_ACC',   cdnc_acc,                               &
                           longname='CDNC occurence acc.+ cloud weighted',units='m-3'     )

  CALL add_stream_element (activ,   'CDNC',       cdnc,                                   &
                           longname='CDNC',units='m-3'                                    )

  CALL add_stream_element (activ,   'CDNC_BURDEN_ACC',cdnc_burden_acc,                    &
                           longname='CDNC burden occurence accumulated',  units='m-2'     )

  CALL add_stream_element (activ,   'CDNC_BURDEN',cdnc_burden,                            &
                           longname='CDNC burden',                        units='m-2'     )

  CALL add_stream_element (activ,   'BURDEN_TIME',burden_time,                            &
                           longname='acc. cdnc burden occ.time fraction', units='1'       )

  CALL add_stream_element (activ,   'LWC_ACC',    lwc_acc,                                &
                           longname='liq wat cont acc.+ cloud weighted',  units='kg m-3'  )

  CALL add_stream_element (activ,   'CLOUD_TIME', cloud_time,                             &
                           longname='acc. cloud occurence time fraction', units='1'       )

  CALL add_stream_element (activ,   'REFFL_ACC',  reffl_acc,                              &
                           longname='cloud drop effectiv radius weighted',units='um'      )

  CALL add_stream_element (activ,   'REFFL_CT',  reffl_ct,                                &
                           longname='cloud top effectiv radius weighted',units='um'       )

  CALL add_stream_element (activ,   'REFFL_TIME',  reffl_time,                            &
                           longname='cloud top effectiv radius occ.time',units='1'        )

  CALL add_stream_element (activ,   'CDNC_CT',  cdnc_ct,                                  &
                           longname='cloud top cloud droplet number conc.',units='cm-3'   )

  CALL add_stream_element (activ,   'IWC_ACC',    iwc_acc,                                &
                           longname='ice wat cont acc.+ cloud weighted',  units='kg m-3'  )

  CALL add_stream_element (activ,   'CLIWC_TIME', cliwc_time,                             &
                           longname='acc. cloud occurence time fraction', units='1'       )

 !>>SF Kasja diags
  CALL default_stream_setting (activ, lrerun=.FALSE.)
  CALL add_stream_element (activ,   'QNUCL',       qnucl,                                 &
                           longname='ND768: CD nucleation rate',                 units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCND',       qcnd,                                   &
                           longname='ND768: condensation rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QLWC_DETR',     qlwc_detr,                           &
                           longname='ND768: xl detrainment rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QEVP_CD', qevp_lwc,                                  &
                           longname='ND768: evaporation of CD',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QAUTN',       qautn,                                 &
                           longname='ND768: autoconversion rate (N)',       units='m-3 s-1' )

  CALL add_stream_element (activ,   'QRACL',       qracl,                                 &
                           longname='ND768: rain accretion rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QRACLN',      qracln,                                &
                           longname='ND768: rain accretion rate (N)',           units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSACL',      qsacl,                                  &
                           longname='ND768: CD accretion rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSACLN',      qsacln,                                  &
                           longname='ND768: CD accretion rate (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QFRZ',       qfrz,                                   &
                           longname='ND768: CD freezing rate',               units='kg kg-1 s-1' )

  !>>UP #783
  CALL add_stream_element (activ,   'QFRZNHET',       qfrznhet,                           &
                           longname='UP783: CD freezing rate, heterogenous (N)', units='m-3 s-1' )
  CALL add_stream_element (activ,   'QFRZNHOM',       qfrznhom,                           &
                           longname='UP783: CD freezing rate, homogenous (N)', units='m-3 s-1' )
  !<<UP #783

  CALL add_stream_element (activ,   'QFRZN',       qfrzn,                                 &
                           longname='ND768: CD freezing rate (N)',               units='m-3 s-1' )

  CALL add_stream_element (activ,   'QNUCI',       qnuci,                                 &
                           longname='ND768: IC nucleation rate',                 units='m-3 s-1' )

  CALL add_stream_element (activ,   'QDEP',       qdep,                                   &
                           longname='ND768: deposition rate',                units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QIWC_DETR',     qiwc_detr,                               &
                           longname='ND768: xi detrainment rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSUB_IWC', qsub_iwc,                                 &
                           longname='ND768: sublimation of IC',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QAGG',       qagg,                                   &
                           longname='ND768: snow aggregation rate',          units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QAGGN',       qaggn,                                   &
                           longname='ND768: snow aggregation rate (N)',          units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSACI',       qsaci,                                 &
                           longname='ND768: IC accretion rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSACIN',       qsacin,                                 &
                           longname='ND768: IC accretion rate (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSELFN',       qselfn,                                 &
                           longname='ND768: IC self collection (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSECPROD',       qsecprod,                                 &
                           longname='ND768: Secondary IC production rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSECPRODN',       qsecprodn,                                 &
                           longname='ND768: Secondary IC production rate (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSEDI',       qsedi,                                 &
                           longname='ND768: IC sedimentation rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSEDIN',       qsedin,                                 &
                           longname='ND768: IC sedimentation rate (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QMLT',   qmlt,                                     &
                           longname='ND768: IC melting rate',                    units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QMLTN',   qmltn,                                     &
                           longname='ND768: IC melting rate (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QEVP_RAIN', qevp_rain,                               &
                           longname='ND768: evaporation of rain',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSUB_SNOW',  qsub_snow,                              &
                           longname='ND768: sublimation of snow',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QRPRN',       qrprn,                                 &
                           longname='ND768: Rain production rate (N)',            units='m-3 s-1' )

  !UP warning: self-collection wrongly included
  CALL add_stream_element (activ,   'QSPRN',   qsprn,                                     &
                           longname='ND768: Snow production rate, wrong with self-collection (N)', &
                           units='m-3 s-1' )

!<<SF Kasja diags

!>>DN: new diags
  CALL add_stream_element (activ,   'QMLTN2',   qmltn2,                                     &
                           longname='ND768: IC melting rate (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QXMLT',   qxmlt,                                     &
                           longname='ND768: IC (falling from above) melting rate',       units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QEVABFN',   qevabfn,                                     &
                           longname='ND768: CD evaporation rate BF (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QEVABF',  qevabf,                              &
                           longname='ND768: evaporation of CD BF',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QDEPBF',  qdepbf,                              &
                           longname='ND768: deposition rate BF',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QCDNC_DETR',   qcdnc_detr,                                     &
                           longname='ND768: cdnc detrainment rate (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QICNC_DETR',   qicnc_detr,                                     &
                           longname='ND768: icnc detrainment rate (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC1',   qcorric1,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC2',   qcorric2,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC3',   qcorric3,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC4',   qcorric4,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC5',   qcorric5,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC6',   qcorric6,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC7',   qcorric7,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD1',   qcorrcd1,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD2',   qcorrcd2,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD3',   qcorrcd3,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD4',   qcorrcd4,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD5',   qcorrcd5,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD6',   qcorrcd6,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD7',   qcorrcd7,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD8',   qcorrcd8,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD9',   qcorrcd9,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD10',   qcorrcd10,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-3 s-1' )

!>>UP: new diags for CD activation correction
! and further split correction terms, #783
  IF (lnewdiags) THEN 
      CALL add_stream_element (activ,   'QCORRCD2UNPHYS2D',   qcorrcd2unphys_2d,  &
                               longname='CDNC correction at cc=0 (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'QCORRCD11',   qcorrcd1_1,  &
                               longname='CDNC correction for min. cqtmin (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'QCORRCD12',   qcorrcd1_2,  &
                               longname='CDNC correction: CD to IC (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'QCORRCD13',   qcorrcd1_3,  &
                               longname='CDNC correction: CD to unphys (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'QCORRCD14',   qcorrcd1_4,  &
                               longname='CDNC correction: IC to CD (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'QCORRCD21',   qcorrcd2_1,  &
                               longname='CDNC correction for min. cqtmin (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'QCORRCD22',   qcorrcd2_2,  &
                               longname='CDNC correction: CD to unphys (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'QCORRIC11',   qcorric1_1,  &
                               longname='ICNC correction for min. cqtmin (N)',            units='m-3 s-1' )
!>>UP #783.3
      CALL add_stream_element (activ,   'QCORRCD23',   qcorrcd2_3,  &
                               longname='CDNC correction: CD to IC (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'QCORRIC12',   qcorric1_2,  &
                               longname='ICNC correction: CD to IC, for IC only (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'QCORRIC13',   qcorric1_3,  &
                               longname='ICNC correction: IC to CD, for IC only (N)',            units='m-3 s-1' )
!<<UP #783.3
  ENDIF
!<<UP: new diags

  CALL add_stream_element (activ,   'QXTTECDNC',   qxttecdnc,                                     &
                           longname='ND768: cdnc tendency (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QXTTEICNC',   qxtteicnc,                                     &
                           longname='ND768: icnc tendency (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRXI',  qcorrxi,                              &
                           longname='ND768: xi correction',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QCORRXI2',  qcorrxi2,                              &
                           longname='ND768: xi correction',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QCORRXL',  qcorrxl,                              &
                           longname='ND768: xl correction',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QCORRXL2',  qcorrxl2,                              &
                           longname='ND768: xl correction',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QXLTE',  qxlte,                              &
                           longname='ND768: xl tendency',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QXITE',  qxite,                              &
                           longname='ND768: xi tendency',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSPR',       qspr,                                 &
                           longname='ND768: Snow production rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QRPR',       qrpr,                                 &
                           longname='ND768: Rain production rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QtestCD',   qtestCD,                                     &
                           longname='ND768: CD test (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QtestLWC',  qtestLWC,                              &
                           longname='ND768: LWC test',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QtestIC',   qtestIC,                                     &
                           longname='ND768: IC test (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QtestIWC',  qtestIWC,                              &
                           longname='ND768: IWC test',            units='kg kg-1 s-1' )

  CALL default_stream_setting (activ, lrerun=.TRUE.)
!<<DN: new diags

  CALL default_stream_setting (activ, laccu=.FALSE., lpost=.FALSE.)

  CALL add_stream_element (activ,   'CLOUD_COVER_DUPLIC', cloud_cover_duplic,             &
                           longname='cloud cover duplicate for record at t+1', units='1'  )


  IF (nic_cirrus>0) THEN

  CALL add_stream_element (activ, 'ICNC_instantaneous', icnc_instantan, &
                           longname='ICNC instantaneous', units='m-3',  &
                           laccu=.FALSE., lpost=.TRUE., lrerun=.TRUE.)

  CALL default_stream_setting (activ, laccu=.TRUE., lpost=.TRUE.)

  CALL add_stream_element (activ,   'ICNC_ACC',   icnc_acc,                               &
                           longname='ICNC occurence acc.+ cloud weighted',units='m-3'     )

  CALL add_stream_element (activ,   'ICNC',       icnc,                                   &
                           longname='ICNC',units='m-3'                                    )

  CALL add_stream_element (activ,   'ICNC_BURDEN_ACC',icnc_burden_acc,                    &
                           longname='ICNC burden occurence accumulated',  units='m-2'     )

  CALL add_stream_element (activ,   'ICNC_BURDEN',icnc_burden,                            &
                           longname='ICNC burden',                        units='m-2'     )

  CALL add_stream_element (activ,   'BURDIC_TIME',burdic_time,                            &
                           longname='acc. icnc burden occ.time fraction', units='1'       )

  CALL add_stream_element (activ,   'REFFI_ACC',  reffi_acc,                              &
                           longname='ice crystal effectiv radius weighted',units='um'     )

  CALL add_stream_element (activ,   'REFFI_TOVS',  reffi_tovs,                            &
                           longname='semi-transparent cirrus effectiv radius',units='um'  )

  CALL add_stream_element (activ,   'REFFI_TIME',  reffi_time,                            &
                           longname='accumulted semi-transp. cirrus time',units='1'       )

  CALL add_stream_element (activ,   'IWP_TOVS',  iwp_tovs,                                &
                           longname='IWP sampled a la TOVS',units='kg m-2'                )

 !>>DN: new diags
  CALL default_stream_setting (activ, lrerun=.FALSE.)

  CALL add_stream_element (activ,   'QGENTL',   qgentl,                                     &
                           longname='ND768: Tompkins cloud cover scheme',               units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QGENTI',   qgenti,                                     &
                           longname='ND768: Tompkins cloud cover scheme',               units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSUB_ICE',   qsub_ice,                                     &
                           longname='ND768: sublimation of IC(falling from above)',     units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QMLT_SNOW',   qmlt_snow,                                     &
                           longname='ND768: snow melting rate (total)',                    units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QMLTS_ATM',   qmlts_atm,                                     &
                           longname='ND768: snow melting rate (atmosphere)',                    units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QMLT_CONV',   qmlt_conv,                                     &
                           longname='ND768: Melting of detrained ice',                    units='kg kg-1 s-1' )

  CALL default_stream_setting (activ, lrerun=.TRUE.)
!<<DN: new diags

  CALL default_stream_setting (activ, laccu=.FALSE.)

  CALL add_stream_element (activ,   'SICE',       sice,                                   &
                           longname='ECHAM supersaturation over ice',     units='% [0-1]' )

  END IF

!>>DN: burden
  CALL default_stream_setting (activ, lrerun=.FALSE.)

  CALL default_stream_setting (activ, laccu=.TRUE.)

  CALL add_stream_element (activ,   'DAUT',       daut,                                   &
                           longname='ND768: CD autoconversion rate',             units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DFRE',       dfre,                                   &
                           longname='ND768: CD freezing rate',                   units='m-2 s-1' )

  CALL add_stream_element (activ,   'DNUCL',       dnucl,                                 &
                           longname='ND768: CD nucleation rate',                 units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCND',       dcnd,                                   &
                           longname='ND768: condensation rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DLWC_DETR',     dlwc_detr,                           &
                           longname='ND768: xl detrainment rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DEVP_CD', devp_lwc,                                  &
                           longname='ND768: evaporation of CD',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DAUTN',       dautn,                                 &
                           longname='ND768: autoconversion rate (N)',       units='m-2 s-1' )

  CALL add_stream_element (activ,   'DRACL',       dracl,                                 &
                           longname='ND768: rain accretion rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DRACLN',      dracln,                                &
                           longname='ND768: rain accretion rate (N)',           units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSACL',      dsacl,                                  &
                           longname='ND768: CD accretion rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSACLN',      dsacln,                                  &
                           longname='ND768: CD accretion rate (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DFRZ',       dfrz,                                   &
                           longname='ND768: CD freezing rate',               units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DFRZN',       dfrzn,                                 &
                           longname='ND768: CD freezing rate (N)',               units='m-2 s-1' )

  !>>UP #783
  CALL add_stream_element (activ,   'DFRZNHET',       dfrznhet,                           &
                           longname='UP783: CD freezing rate, heterogenous (N)', units='m-2 s-1' )
  CALL add_stream_element (activ,   'DFRZNHOM',       dfrznhom,                           &
                           longname='UP783: CD freezing rate, homogenous (N)', units='m-2 s-1' )
  !<<UP #783

  CALL add_stream_element (activ,   'DNUCI',       dnuci,                                 &
                           longname='ND768: IC nucleation rate',                 units='m-2 s-1' )

  CALL add_stream_element (activ,   'DDEP',       ddep,                                   &
                           longname='ND768: deposition rate',                units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DIWC_DETR',     diwc_detr,                               &
                           longname='ND768: xi detrainment rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSUB_IWC', dsub_iwc,                                 &
                           longname='ND768: sublimation of IC',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DAGG',       dagg,                                   &
                           longname='ND768: snow aggregation rate',          units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DAGGN',       daggn,                                   &
                           longname='ND768: snow aggregation rate (N)',          units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSACI',       dsaci,                                 &
                           longname='ND768: IC accretion rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSACIN',       dsacin,                                 &
                           longname='ND768: IC accretion rate (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSELFN',       dselfn,                                 &
                           longname='ND768: IC self collection (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSECPROD',       dsecprod,                                 &
                           longname='ND768: Secondary IC production rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSECPRODN',       dsecprodn,                                 &
                           longname='ND768: Secondary IC production rate (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSEDI',       dsedi,                                 &
                           longname='ND768: IC sedimentation rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSEDIN',       dsedin,                                 &
                           longname='ND768: IC sedimentation rate (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DMLT',   dmlt,                                     &
                           longname='ND768: IC melting rate',                    units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DMLTN',   dmltn,                                     &
                           longname='ND768: IC melting rate (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DEVP_RAIN', devp_rain,                               &
                           longname='ND768: evaporation of rain',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSUB_SNOW',  dsub_snow,                              &
                           longname='ND768: sublimation of snow',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DRPRN',       drprn,                                 &
                           longname='ND768: Rain production rate (N)',            units='m-2 s-1' )

  !UP warning: self-collection wrongly included
  CALL add_stream_element (activ,   'DSPRN',   dsprn,                                     &
                           longname='ND768: Snow production rate, wrong with self-collection (N)', &
                           units='m-2 s-1' )

  CALL add_stream_element (activ,   'DMLTN2',   dmltn2,                                     &
                           longname='ND768: IC melting rate (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DXMLT',   dxmlt,                                     &
                           longname='ND768: IC (falling from above) melting rate',       units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DEVABFN',   devabfn,                                     &
                           longname='ND768: CD evaporation rate BF (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DEVABF',  devabf,                              &
                           longname='ND768: evaporation of CD BF',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DDEPBF',  ddepbf,                              &
                           longname='ND768: deposition rate BF',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCDNC_DETR',   dcdnc_detr,                                     &
                           longname='ND768: cdnc detrainment rate (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DICNC_DETR',   dicnc_detr,                                     &
                           longname='ND768: icnc detrainment rate (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC1',   dcorric1,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC2',   dcorric2,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC3',   dcorric3,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC4',   dcorric4,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC5',   dcorric5,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC6',   dcorric6,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC7',   dcorric7,                                     &
                           longname='ND768: icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD1',   dcorrcd1,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD2',   dcorrcd2,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD3',   dcorrcd3,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD4',   dcorrcd4,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD5',   dcorrcd5,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD6',   dcorrcd6,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD7',   dcorrcd7,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD8',   dcorrcd8,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD9',   dcorrcd9,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD10',   dcorrcd10,                                     &
                           longname='ND768: cdnc correction (N)',                    units='m-2 s-1' )

!>>UP: new diags for CD activation correction
! and further split correction terms, #783
  IF (lnewdiags) THEN 
      CALL add_stream_element (activ,   'DCORRCD2UNPHYS2D',   dcorrcd2unphys_2d,  &
                               longname='CDNC correction at cc=0 (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'DCORRCD11',   dcorrcd1_1,  &
                               longname='CDNC correction for min. cqtmin (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'DCORRCD12',   dcorrcd1_2,  &
                               longname='CDNC correction: CD to IC (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'DCORRCD13',   dcorrcd1_3,  &
                               longname='CDNC correction: CD to unphys (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'DCORRCD14',   dcorrcd1_4,  &
                               longname='CDNC correction: IC to CD (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'DCORRCD21',   dcorrcd2_1,  &
                               longname='CDNC correction for min. cqtmin (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'DCORRCD22',   dcorrcd2_2,  &
                               longname='CDNC correction: CD to unphys (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'DCORRIC11',   dcorric1_1,  &
                               longname='ICNC correction for min. cqtmin (N)',            units='m-3 s-1' )
!>>UP #783.3
      CALL add_stream_element (activ,   'DCORRCD23',   dcorrcd2_3,  &
                               longname='CDNC correction: CD to IC (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'DCORRIC12',   dcorric1_2,  &
                               longname='ICNC correction: CD to IC, for IC only (N)',            units='m-3 s-1' )
      CALL add_stream_element (activ,   'DCORRIC13',   dcorric1_3,  &
                               longname='ICNC correction: IC to CD, for IC only (N)',            units='m-3 s-1' )
!<<UP #783.3
  ENDIF
!<<UP: new diags

  CALL add_stream_element (activ,   'DXTTECDNC',   dxttecdnc,                                     &
                           longname='ND768: cdnc tendency (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DXTTEICNC',   dxtteicnc,                                     &
                           longname='ND768: icnc tendency (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRXI',  dcorrxi,                              &
                           longname='ND768: xi correction',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRXI2',  dcorrxi2,                              &
                           longname='ND768: xi correction',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRXL',  dcorrxl,                              &
                           longname='ND768: xl correction',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRXL2',  dcorrxl2,                              &
                           longname='ND768: xl correction',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DXLTE',  dxlte,                              &
                           longname='ND768: xl tendency',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DXITE',  dxite,                              &
                           longname='ND768: xi tendency',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSPR',       dspr,                                 &
                           longname='ND768: Snow production rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DRPR',       drpr,                                 &
                           longname='ND768: Rain production rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DtestCD',   dtestCD,                                     &
                           longname='ND768: CD test (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DtestLWC',  dtestLWC,                              &
                           longname='ND768: LWC test',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DtestIC',   dtestIC,                                     &
                           longname='ND768: IC test (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DtestIWC',  dtestIWC,                              &
                           longname='ND768: IWC test',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DGENTL',   dgentl,                                     &
                           longname='ND768: Tompkins cloud cover scheme',               units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DGENTI',   dgenti,                                     &
                           longname='ND768: Tompkins cloud cover scheme',               units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSUB_ICE',   dsub_ice,                                     &
                           longname='ND768: sublimation of IC(falling from above)',     units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DMLT_SNOW',   dmlt_snow,                                     &
                           longname='ND768: snow melting rate (total)',                    units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DMLTS_ATM',   dmlts_atm,                                     &
                           longname='ND768: snow melting rate (atmosphere)',                    units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DMLTS_SFC',   dmlts_sfc,                                     &
                           longname='ND768: snow melting rate (surface)',                    units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSEDI_SFC',       dsedi_sfc,                                 &
                           longname='ND768: IC sedimentation rate (surface)',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSNOW',       dsnow,                                 &
                           longname='ND768: snow fall (stratiform)',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCLDTTE',       dcldtte,                                 &
                           longname='ND768: 2-moment cloud microphysics energy tendency',              units='W m-2' )

  CALL add_stream_element (activ,   'DCONVTTE',       dconvtte,                                 &
                           longname='ND768: cumulus parameterization energy tendency',              units='W m-2' )

  CALL add_stream_element (activ,   'DMLT_CONV',   dmlt_conv,                                     &
                           longname='ND768: Melting of detrained ice',                    units='kg m-2 s-1' )

!>>UP: new cc diags, #783
  IF (lnewdiags) THEN 
      CALL add_stream_element (activ,   'ICC_MOD',   icc_mod,  &
                               longname='ice_cloud_cover as diagnosed by model',            units='1', &
                               laccu=.FALSE. )

      CALL add_stream_element (activ,   'LCC_MOD',   lcc_mod,  &
                               longname='liquid_cloud_cover as diagnosed by model',            units='1', &
                               laccu=.FALSE. )

      CALL add_stream_element (activ,   'MCC_MOD',   mcc_mod,  &
                               longname='mixed_phase_cloud_cover as diagnosed by model',            units='1', &
                               laccu=.FALSE. )

      CALL add_stream_element (activ,   'CNT_CQTMINS1',   diag_cnt_cqtmins1,  &
                               longname='Count how often the min. CDNC condition is hit, 1', units='1', &
                               laccu=.TRUE. )

      CALL add_stream_element (activ,   'CNT_CQTMINS2',   diag_cnt_cqtmins2,  &
                               longname='Count how often the min. CDNC condition is hit, 2', units='1', &
                               laccu=.TRUE. )

      CALL add_stream_element (activ,   'CNT_CQTMINS3',   diag_cnt_cqtmins3,  &
                               longname='Count how often the min. CDNC condition is hit, 3', units='1', &
                               laccu=.TRUE. )

      CALL add_stream_element (activ,   'CNT_CQTMINL1',   diag_cnt_cqtminl1,  &
                               longname='Count how often the min. CDNC condition is hit, 1', units='1', &
                               laccu=.TRUE. )

      CALL add_stream_element (activ,   'CNT_CQTMINL2',   diag_cnt_cqtminl2,  &
                               longname='Count how often the min. CDNC condition is hit, 2', units='1', &
                               laccu=.TRUE. )

      CALL add_stream_element (activ,   'CNT_CQTMINL3',   diag_cnt_cqtminl3,  &
                               longname='Count how often the min. CDNC condition is hit, 3', units='1', &
                               laccu=.TRUE. )

      CALL add_stream_element (activ,   'CNT_LL12D_1',   diag_cnt_ll12d_1,  &
                               longname='Count how often activation 1 happens', units='1', &
                               laccu=.TRUE. )

      CALL add_stream_element (activ,   'CNT_LL12D_2',   diag_cnt_ll12d_2,  &
                               longname='Count how often activation 2 happens', units='1', &
                               laccu=.TRUE. )

  DO jccn = 1, size(cdnc_bins)-1
        DO jw = 1, size(lwc_bins)-1
          CALL add_stream_element (activ, 'cdnc_lwc_binned_cdnc'//cdnc_bins_string(jccn)//'-'//cdnc_bins_string(jccn+1)//&
                                          '_lwc'//lwc_bins_string(jw)//'-'//lwc_bins_string(jw+1), &
                    cdnc_lwc_binned(jccn*jw)%ptr, &
                    longname='Binned frequency of CDNC [m-3] vs LWC [kg/kg]', units='', &
                    laccu=.TRUE. )
        ENDDO
  ENDDO
  DO jccn = 1, size(prcp_bins)-1
          CALL add_stream_element (activ, 'prcp_binned_min'//prcp_bins_string(jccn)//'_max'//prcp_bins_string(jccn+1), &
                    prcp_binned(jccn)%ptr, &
                    longname='Precipitation, binned', units='kg m-2 s-1', &
                    laccu=.TRUE. )
  ENDDO
  DO jccn = 1, size(cdnc_bins)-1
          CALL add_stream_element (activ, 'cdnc_binned_min'//cdnc_bins_string(jccn)//'_max'//cdnc_bins_string(jccn+1), &
                    cdnc_binned(jccn)%ptr, &
                    longname='CDNC binned', units='m-3', &
                    laccu=.TRUE. )
  ENDDO
  ENDIF
!<<UP: new diags

!>>UP #844
  IF (lslf) THEN
       !IWP and LWP divided by temperature regime
       CALL add_stream_element (activ, 'slfdiag_iwp_mxpT', slfdiag_iwp_mxpT,     &
                                longname='IWP in mixed phase temperatures', &
                                units='g/m2', &
                                ! keep this as DN implemented it for now
                                laccu = .TRUE. )
       CALL add_stream_element (activ, 'slfdiag_lwp_mxpT', slfdiag_lwp_mxpT,     &
                                longname='LWP in mixed phase temperatures', &
                                units='g/m2', &
                                ! keep this as DN implemented it for now
                                laccu = .TRUE. )
       CALL add_stream_element (activ, 'slfdiag_iwp_l35', slfdiag_iwp_l35,     &
                                longname='IWP at T<-35C', &
                                units='g/m2', &
                                ! keep this as DN implemented it for now
                                laccu = .TRUE. )
       CALL add_stream_element (activ, 'slfdiag_lwp_g0', slfdiag_lwp_g0,     &
                                longname='LWP at T>0C', &
                                units='g/m2', &
                                ! keep this as DN implemented it for now
                                laccu = .TRUE. )
       ! NUMBERS divided by temperature regime
       CALL add_stream_element (activ, 'slfdiag_icnc_mxpT', slfdiag_icnc_mxpT,     &
                                longname='ICNC in mixed phase temperatures', &
                                units='/m2', &
                                laccu = .TRUE. )
       CALL add_stream_element (activ, 'slfdiag_cdnc_mxpT', slfdiag_cdnc_mxpT,     &
                                longname='CDNC in mixed phase temperatures', &
                                units='/m2', &
                                laccu = .TRUE. )
       CALL add_stream_element (activ, 'slfdiag_icnc_l35', slfdiag_icnc_l35,     &
                                longname='ICNC at T<-35C', &
                                units='/m2', &
                                laccu = .TRUE. )
       CALL add_stream_element (activ, 'slfdiag_cdnc_g0', slfdiag_cdnc_g0,     &
                                longname='ICNC at T>0C', &
                                units='/m2', &
                                laccu = .TRUE. )
  ENDIF
!<<UP #844
  IF (lccnclimdiags) THEN 
  !>>UP Additional diagnostics to evaluate/improve ccnclim
  CALL add_stream_element (activ, 'CBASE_LIFETIME', cloudbase_lifetime,     &
            longname='How many time steps has a cloud base been here previously', &
            units='', &
            laccu = .FALSE. )

  CALL add_stream_element (activ, 'CCN_STRAT_TMP', ccn_strat_tmp, &
            longname='Number of activated aerosols in stratiform clouds (not acc.)', units='m-3', &
            laccu=.FALSE. )

  CALL add_stream_element (activ, 'CCN_CONV_TMP', ccn_conv_tmp, &
            longname='Number of activated aerosols in convective clouds (not acc.)', units='m-3', &
            laccu=.FALSE. )

  CALL add_stream_element (activ, 'NACT_STRAT_TMP', nact_strat_tmp, &
            longname='Number of activated CDs in stratiform clouds (not acc.)', units='m-3', &
            laccu=.FALSE. )

  CALL add_stream_element (activ, 'NACT_CONV_TMP', nact_conv_tmp, &
            longname='Number of activated CDs in convective clouds (not acc.)', units='m-3', &
            laccu=.FALSE. )

  CALL add_stream_element (activ, 'CCN_STRAT_2D', ccn_strat_2d, &
            longname='Number of activated aerosols in stratiform clouds, vertically integrated', units='m-2', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'RWET_KS_DIAG', rwet_KS_diag, &
            longname='Wet radius of KS, acc.', units='m', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'RWET_KS_DIAG_CLOUDBASE', rwet_KS_diag_cloudbase, &
            longname='Wet radius of KS, acc., counting only at cloud base', units='m', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'RWET_KS_DIAG_CLOUDBASE_CTR', rwet_KS_diag_cloudbase_ctr, &
            longname='Wet radius of KS, acc., counting only at cloud base, counter', units='', &
            laccu=.TRUE. )

  DO jccn = 1, size(ccn_bins)-1
          CALL add_stream_element (activ, 'ccn_strat_binned_min'//ccn_bins_string(jccn)//'_max'//ccn_bins_string(jccn+1), &
                    ccn_strat_binned(jccn)%ptr, &
                    longname='Number of activated aerosols in stratiform clouds, binned', units='m-3', &
                    laccu=.TRUE. )
          CALL add_stream_element (activ, 'ccn_conv_binned_min'//ccn_bins_string(jccn)//'_max'//ccn_bins_string(jccn+1), &
                    ccn_conv_binned(jccn)%ptr, &
                    longname='Number of activated aerosols in stratiform clouds, binned', units='m-3', &
                    laccu=.TRUE. )
  ENDDO

  DO jccn = 1, size(ccn_bins)-1
        DO jw = 1, size(w_bins)-1
          CALL add_stream_element (activ, 'ccn_strat_w_binned_ccn'//ccn_bins_string(jccn)//'-'//ccn_bins_string(jccn+1)//&
                                          '_w'//w_bins_string(jw)//'-'//w_bins_string(jw+1), &
                    ccn_strat_w_binned(jccn*jw)%ptr, &
                    longname='Number of activated aerosols in stratiform clouds, binned by CCN and updraft', units='m-3', &
                    laccu=.TRUE. )

          CALL add_stream_element (activ, 'ccn_strat_w_actccn_binned_ccn'//ccn_bins_string(jccn)//'-'//ccn_bins_string(jccn+1)//&
                                          '_w'//w_bins_string(jw)//'-'//w_bins_string(jw+1), &
                    ccn_strat_w_actccn_binned(jccn*jw)%ptr, &
                    longname='Number of activated aerosols in stratiform clouds, where ACTCCN logical, binned by CCN and updraft', &
                    units='m-3', &
                    laccu=.TRUE. )
        ENDDO
        DO jw = 1, size(cbl_bins)-1
          CALL add_stream_element (activ, 'ccn_strat_cbl_binned_ccn'//ccn_bins_string(jccn)//'-'//ccn_bins_string(jccn+1)//&
                                          '_cbl'//cbl_bins_string(jw)//'-'//cbl_bins_string(jw+1), &
                    ccn_strat_cbl_binned(jccn*jw)%ptr, &
                    longname='Number of activated aerosols in stratiform clouds, binned by CCN and cloud base lifetime', &
                    units='m-3', &
                    laccu=.TRUE. )
        ENDDO
  ENDDO

  ENDIF
!<<UP: new diags

  IF (lccnclim_diag) THEN

  CALL add_stream_element (activ, 'CCN_STRAT_CDNCMIN', ccn_strat_cdncmin, &
            longname='Number of activated aerosols in stratiform clouds, only above cdnc_min', units='m-3', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_CDNCMIN_CTR', ccn_strat_cdncmin_ctr, &
            longname='Number of activated aerosols in stratiform clouds, only above cdnc_min - counter', units='', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_ACTCCN', ccn_strat_actccn, &
            longname='Number of activated aerosols in stratiform clouds, only where logical for activation', units='m-3', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_ACTCCN_CTR', ccn_strat_actccn_ctr, &
            longname='Number of activated aerosols in stratiform clouds, only where logical for activation - counter', units='', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_NOACTCCN', ccn_strat_noactccn, &
            longname='Number of activated aerosols in stratiform clouds, only where no logical for activation', units='m-3', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_NOACTCCN_CTR', ccn_strat_noactccn_ctr, &
            longname='Number of activated aerosols in stratiform clouds, only where no logical for activation - counter', units='', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_ACTCCNWOBASE', ccn_strat_actccnwobase, &
            longname='CCN STRAT acc., only where logical for activation (excluding cloud base check)', units='m-3', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_ACTCCNWOBASE_CTR', ccn_strat_actccnwobase_ctr, &
            longname='CCN STRAT acc., only where logical for activation (excluding cloud base check) - counter', &
            units='', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_NOSMALL40', ccn_strat_nosmall40, &
            longname='Number of activated aerosols in stratiform clouds, only above 40e6m-3', units='m-3', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_NOSMALL40_CTR', ccn_strat_nosmall40_ctr, &
            longname='Number of activated aerosols in stratiform clouds, only above 40e6m-3 - counter', units='', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_NOSMALL1', ccn_strat_nosmall1, &
            longname='Number of activated aerosols in stratiform clouds, only above 1e6m-3', units='m-3', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_NOSMALL1_CTR', ccn_strat_nosmall1_ctr, &
            longname='Number of activated aerosols in stratiform clouds, only above 1e6m-3 - counter', units='', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_SO4KS', ccn_strat_so4ks, &
            longname='Number of activated aerosols in stratiform clouds, SO4KS', units='m-3', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_MAX', ccn_strat_max, &
            longname='Number of activated aerosols in stratiform clouds, record maximum concentration', units='m-3', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'CCN_STRAT_MAX_CTR', ccn_strat_max_ctr, &
            longname='Number of activated aerosols in stratiform clouds, record maximum concentration - counter', units='', &
            laccu=.TRUE. )

  ENDIF
!<<UP: new ccnclim diags

!>>UP: new diags
  IF (lnewdiags) THEN
  CALL add_stream_element (activ, 'RWET_KS_DIAG', rwet_KS_diag, &
            longname='Wet radius of KS, acc.', units='m', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'RWET_KS_DIAG_CLOUDBASE', rwet_KS_diag_cloudbase, &
            longname='Wet radius of KS, acc., counting only at cloud base', units='m', &
            laccu=.TRUE. )

  CALL add_stream_element (activ, 'RWET_KS_DIAG_CLOUDBASE_CTR', rwet_KS_diag_cloudbase_ctr, &
            longname='Wet radius of KS, acc., counting only at cloud base, counter', units='', &
            laccu=.TRUE. )

!>>UP lifetime diagnostics
  CALL add_stream_element (activ, 'Q_LWC_SOURCES', q_lwc_sources, &
            longname='Sources of LWC', units='kg kg-1 per timestep', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'Q_LWC_SINKS', q_lwc_sinks, &
            longname='Sinks of LWC', units='kg kg-1 per timestep', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'D_LWC_LIFETIME_SOURCES', d_lwc_lifetime_sources, &
            longname='LWC lifetime computed from sources', units='timesteps', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'D_LWC_LIFETIME_SINKS', d_lwc_lifetime_sinks, &
            longname='LWC lifetime computed from sinks', units='timesteps', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'Q_IWC_SOURCES', q_iwc_sources, &
            longname='Sources of IWC', units='kg kg-1 per timestep', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'Q_IWC_SINKS', q_iwc_sinks, &
            longname='Sinks of IWC', units='kg kg-1 per timestep', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'D_IWC_LIFETIME_SOURCES', d_iwc_lifetime_sources, &
            longname='IWC lifetime computed from sources', units='timesteps', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'D_IWC_LIFETIME_SINKS', d_iwc_lifetime_sinks, &
            longname='IWC lifetime computed from sinks', units='timesteps', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'Q_ICNC_SOURCES', q_icnc_sources, &
            longname='Sources of ICNC', units='m-3 per timestep', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'Q_ICNC_SINKS', q_icnc_sinks, &
            longname='Sinks of ICNC', units='m-3 per timestep', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'D_ICNC_LIFETIME_SOURCES', d_icnc_lifetime_sources, &
            longname='ICNC lifetime computed from sources', units='timesteps', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'D_ICNC_LIFETIME_SINKS', d_icnc_lifetime_sinks, &
            longname='ICNC lifetime computed from sinks', units='timesteps', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'Q_CDNC_SOURCES', q_cdnc_sources, &
            longname='Sources of CDNC', units='m-3 per timestep', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'Q_CDNC_SINKS', q_cdnc_sinks, &
            longname='Sinks of CDNC', units='m-3 per timestep', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'D_CDNC_LIFETIME_SOURCES', d_cdnc_lifetime_sources, &
            longname='CDNC lifetime computed from sources', units='timesteps', &
            laccu=.TRUE. )
  CALL add_stream_element (activ, 'D_CDNC_LIFETIME_SINKS', d_cdnc_lifetime_sinks, &
            longname='CDNC lifetime computed from sinks', units='timesteps', &
            laccu=.TRUE. )
!<<UP


  ENDIF
!<<UP: new diags

  CALL default_stream_setting (activ, lrerun=.TRUE.)
!<<DN: burden

  CALL default_stream_setting (activ, laccu=.FALSE.)

END SUBROUTINE construct_activ_stream

!>>UP Additional diagnostics to evaluate/improve ccnclim
SUBROUTINE ham_activ_diag_lin_leaitch_2(kproma, kbdim, klev, krow, pcdncact, pw, ll_actccn, ll_bas, pdpg)

    USE mo_conv,         ONLY: cdncact_cv, na_cv
    USE mo_time_control, ONLY: delta_time
    USE mo_ham,            ONLY: lccnclim_diag !UP

    INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(IN) :: pcdncact(kbdim,klev) ! number of activated particles
    REAL(dp), INTENT(IN) :: pw(kbdim,klev) ! number of activated particles
    LOGICAL, INTENT(IN) :: ll_actccn(kbdim,klev) ! whether these CCN are
    ! potentially used for CD activation
    INTEGER, INTENT(IN)  :: ll_bas(kbdim,klev) ! logical for cloud base
    REAL(dp):: pdpg(kbdim,klev)        !< delta p over g [kg/m2]

    LOGICAL  :: ll1(kbdim,klev)
    REAL(dp) :: zeps
    INTEGER  :: jccn, jw, jl, jk
    REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim)

    zeps = EPSILON(1._dp)

    ll1(1:kproma,:) = (na(1:kproma,:,krow) > zeps)

    nact_strat_tmp(1:kproma,:,krow)  = MERGE(pcdncact(1:kproma,:), 0._dp, ll1(1:kproma,:))
    nact_conv_tmp(1:kproma,:,krow)   = MERGE(cdncact_cv(1:kproma,:,krow), 0._dp, ll1(1:kproma,:))
    !nact_strat_tmp(1:kproma,:,krow)  = MERGE(pcdncact(1:kproma,:), nact_strat_tmp(1:kproma,:,krow), ll1(1:kproma,:))
    !nact_conv_tmp(1:kproma,:,krow)   = MERGE(cdncact_cv(1:kproma,:,krow), nact_conv_tmp(1:kproma,:,krow), ll1(1:kproma,:))

    ccn_strat_tmp(1:kproma,:,krow) = na(1:kproma,:,krow)
    ccn_conv_tmp(1:kproma,:,krow)  = na_cv(1:kproma,:,krow)

    ! If this level is a cloud base (ll_bas is filled with its own level if it
    ! contains a cloud base, else 0), raise the counter, else set it to 0 and
    ! start counting from new 
    cloudbase_lifetime(1:kproma,:,krow) = MERGE(cloudbase_lifetime(1:kproma,:,krow) + 1._dp, 0._dp, ll_bas > 0._dp)

    IF (lccnclim_diag) THEN
            ccn_strat_actccn(1:kproma,:,krow) = ccn_strat_actccn(1:kproma,:,krow) + & 
                                                MERGE(na(1:kproma,:,krow) * delta_time, 0._dp, ll_actccn(1:kproma,:))
            ccn_strat_actccn_ctr(1:kproma,:,krow) = ccn_strat_actccn_ctr(1:kproma,:,krow) + &
                                                    MERGE(1._dp*delta_time, 0._dp, ll_actccn(1:kproma,:))
    ENDIF

    !>>UP Diagnose CCN histogram
    DO jccn = 1, size(ccn_bins)-1
            ll1(1:kproma,:) = (na(1:kproma,:,krow) .GE. ccn_bins(jccn)) .AND. (na(1:kproma,:,krow) < ccn_bins(jccn + 1))
            ccn_strat_binned(jccn)%ptr(1:kproma,:,krow) = ccn_strat_binned(jccn)%ptr(1:kproma,:,krow) + &
                                                          MERGE(1._dp * delta_time, 0._dp, ll1(1:kproma,:))
            ll1(1:kproma,:) = (na_cv(1:kproma,:,krow) .GE. ccn_bins(jccn)) .AND. (na_cv(1:kproma,:,krow) < ccn_bins(jccn + 1))
            ccn_conv_binned(jccn)%ptr(1:kproma,:,krow) = ccn_conv_binned(jccn)%ptr(1:kproma,:,krow) + &
                                                         MERGE(1._dp * delta_time, 0._dp, ll1(1:kproma,:))
            ! Diagnose CCN vs W
            DO jw = 1, size(w_bins)-1
                ll1(1:kproma,:) = (na(1:kproma,:,krow) .GE. ccn_bins(jccn)) .AND. (na(1:kproma,:,krow) < ccn_bins(jccn + 1)) .AND. &
                                  (pw(1:kproma,:) .GE. w_bins(jw)) .AND. (pw(1:kproma,:) < w_bins(jw + 1))
                ccn_strat_w_binned(jccn*jw)%ptr(1:kproma,:,krow) = ccn_strat_w_binned(jccn*jw)%ptr(1:kproma,:,krow) + &
                                                              MERGE(1._dp * delta_time, 0._dp, ll1(1:kproma,:))
                ll1(1:kproma,:) = ll1(1:kproma,:) .AND. ll_actccn(1:kproma,:)
                ccn_strat_w_actccn_binned(jccn*jw)%ptr(1:kproma,:,krow) = ccn_strat_w_actccn_binned(jccn*jw)%ptr(1:kproma,:,krow) + &
                                                              MERGE(1._dp * delta_time, 0._dp, ll1(1:kproma,:))
            ENDDO

            ! Diagnose CCN vs cloudbase lifetime
           DO jl = 1, size(cbl_bins)-1
                ll1(1:kproma,:) = (na(1:kproma,:,krow) .GE. ccn_bins(jccn)) .AND. (na(1:kproma,:,krow) < ccn_bins(jccn + 1)) .AND. &
                                  (cloudbase_lifetime(1:kproma,:,krow) .GE. cbl_bins(jl)) .AND.  &
                                  (cloudbase_lifetime(1:kproma,:,krow) < cbl_bins(jl + 1))
                ccn_strat_cbl_binned(jccn*jl)%ptr(1:kproma,:,krow) = ccn_strat_cbl_binned(jccn*jl)%ptr(1:kproma,:,krow) + &
                                                              MERGE(1._dp * delta_time, 0._dp, ll1(1:kproma,:))

           ENDDO
    ENDDO
    DO jk = 1,klev
         ccn_strat_2d(1:kproma,krow)    = ccn_strat_2d(1:kproma,krow) + delta_time*na(1:kproma,jk,krow)*pdpg(1:kproma,jk)
    END DO
    !<<UP
 
END SUBROUTINE ham_activ_diag_lin_leaitch_2

END MODULE mo_activ
