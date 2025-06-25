MODULE mo_activ

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  USE mo_submodel_diag, ONLY: vmem3d

  IMPLICIT NONE

  PUBLIC activ_initialize
  PUBLIC activ_updraft
  PUBLIC activ_lin_leaitch
  PUBLIC construct_activ_stream

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
!>>DN convective mass flux diagnostics
  REAL(dp),        PUBLIC, POINTER :: w_conv(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: mc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: w_conv_inst(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: mc_inst(:,:,:)
!<<DN
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
!>>DN: new diagnostics
  REAL(dp),        PUBLIC, POINTER :: cloud_tm2(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cloud_cover_begin_micro_2m(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: pnucnic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: pnucncd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: pcdncactevap(:,:,:)    
  REAL(dp),        PUBLIC, POINTER :: pevapncd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: pevapncdbf(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: pcdncactfrz(:,:,:)     
  REAL(dp),        PUBLIC, POINTER :: pfrzncd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: picncactmelt(:,:,:)    
  REAL(dp),        PUBLIC, POINTER :: pmeltncd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: picncactsub(:,:,:)     
  REAL(dp),        PUBLIC, POINTER :: psubnic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: dpg(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: pcdncactap(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: picncactap(:,:,:)
!!<<DN: new diagnostics
  REAL(dp),        PUBLIC, POINTER :: sice(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl_ct(:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl_time(:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_ct(:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_incl_ct(:,:)!davidn
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
                                      ncd_activ, nactivpdf, nic_cirrus, lcdnc_progn, &
                                      cdnc_min_fixed
    USE mo_tracer,             ONLY: get_tracer
!davidn
    USE mo_param_switches,     ONLY: tun47ccraut,tun47ccsaut,tun47cdncmin
    USE mo_param_switches,     ONLY: tun31ccraut,tun31ccsaut,tun31cdncmin
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

    !
    !-- overwrite values for coupled CDNC/ICNC cloud scheme
    !
    IF (lcdnc_progn)  THEN
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
                  END SELECT
            END SELECT
            !davidn tuning
            ccsaut = tun31ccsaut
            ccraut = tun31ccraut
            !davidn tuning
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
       !davidn tuning
       cdnc_min_fixed = tun47cdncmin
       cdnc_min_fixed = tun31cdncmin
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
  USE mo_param_switches, ONLY: ncd_activ, nactivpdf, nic_cirrus !SF

  IMPLICIT NONE

  INTEGER           :: jw
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

!>>DN convective mass flux diagnostics
  CALL add_stream_element (activ,   'mc_inst',    mc_inst,                                &
                           longname='grid-mean convective mass flux (inst.)',             &
                           units='kg m-2 s-1', lpost=.FALSE.   )

  CALL add_stream_element (activ,   'W_CONV_INST',    w_conv_inst,                        &
                           longname='grid-mean convective velocity (inst.)',              &
                           units='m s-1', lpost=.FALSE.   )
!<<DN

  CALL add_stream_element (activ,   'REFFL',      reffl,                                  &
                           longname='cloud drop effectiv radius',         units='um'      )

  IF (nic_cirrus>0) THEN

  CALL add_stream_element (activ,   'REFFI',      reffi,                                  &
                           longname='ice crystal effectiv radius',        units='um'      )
  END IF

  CALL add_stream_element (activ,   'NA',         na,                                     &
                           longname='aerosol number for activation',      units='m-3'     )

  CALL default_stream_setting (activ, laccu=.TRUE.)

!>>DN convective mass flux diagnostics
  CALL add_stream_element (activ,   'mc',     mc,                                         &
                           longname='grid-mean convective mass flux (upd+downd)',         &
                           units='kg m-2 s-1'   )

  CALL add_stream_element (activ,   'W_CONV',    w_conv,                                  &
                           longname='grid-mean convective velocity (upd+downd)',          &
                           units='m s-1'   )
!<<DN

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
!>>DN: new diagnostics
  CALL add_stream_element (activ,   'QEVA',       qeva,                                   &
                           longname='CD evaporation rate',                units='m-3 s-1' )
!<<DN: new diagnostics

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

!davidn
  CALL add_stream_element (activ,   'CDNC_INCL_CT',  cdnc_incl_ct,                        &
                           longname='cloud top (inc cloud) cloud droplet number conc.',units='cm-3'   )
!davidn

  CALL add_stream_element (activ,   'IWC_ACC',    iwc_acc,                                &
                           longname='ice wat cont acc.+ cloud weighted',  units='kg m-3'  )

  CALL add_stream_element (activ,   'CLIWC_TIME', cliwc_time,                             &
                           longname='acc. cloud occurence time fraction', units='1'       )

 !>>SF Kasja diags
  CALL default_stream_setting (activ, lrerun=.FALSE.)
  CALL add_stream_element (activ,   'QNUCL',       qnucl,                                 &
                           longname='CD nucleation rate',                 units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCND',       qcnd,                                   &
                           longname='condensation rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QLWC_DETR',     qlwc_detr,                           &
                           longname='xl detrainment rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QEVP_CD', qevp_lwc,                                  &
                           longname='evaporation of CD',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QAUTN',       qautn,                                 &
                           longname='autoconversion rate (N)',       units='m-3 s-1' )

  CALL add_stream_element (activ,   'QRACL',       qracl,                                 &
                           longname='rain accretion rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QRACLN',      qracln,                                &
                           longname='rain accretion rate (N)',           units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSACL',      qsacl,                                  &
                           longname='CD accretion rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSACLN',      qsacln,                                  &
                           longname='CD accretion rate (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QFRZ',       qfrz,                                   &
                           longname='CD freezing rate',               units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QFRZN',       qfrzn,                                 &
                           longname='CD freezing rate (N)',               units='m-3 s-1' )

  CALL add_stream_element (activ,   'QNUCI',       qnuci,                                 &
                           longname='IC nucleation rate',                 units='m-3 s-1' )

  CALL add_stream_element (activ,   'QDEP',       qdep,                                   &
                           longname='deposition rate',                units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QIWC_DETR',     qiwc_detr,                               &
                           longname='xi detrainment rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSUB_IWC', qsub_iwc,                                 &
                           longname='sublimation of IC',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QAGG',       qagg,                                   &
                           longname='snow aggregation rate',          units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QAGGN',       qaggn,                                   &
                           longname='snow aggregation rate (N)',          units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSACI',       qsaci,                                 &
                           longname='IC accretion rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSACIN',       qsacin,                                 &
                           longname='IC accretion rate (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSELFN',       qselfn,                                 &
                           longname='IC self collection (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSECPROD',       qsecprod,                                 &
                           longname='Secondary IC production rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSECPRODN',       qsecprodn,                                 &
                           longname='Secondary IC production rate (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSEDI',       qsedi,                                 &
                           longname='IC sedimentation rate',              units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSEDIN',       qsedin,                                 &
                           longname='IC sedimentation rate (N)',              units='m-3 s-1' )

  CALL add_stream_element (activ,   'QMLT',   qmlt,                                     &
                           longname='IC melting rate',                    units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QMLTN',   qmltn,                                     &
                           longname='IC melting rate (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QEVP_RAIN', qevp_rain,                               &
                           longname='evaporation of rain',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSUB_SNOW',  qsub_snow,                              &
                           longname='sublimation of snow',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QRPRN',       qrprn,                                 &
                           longname='Rain production rate (N)',            units='m-3 s-1' )

  CALL add_stream_element (activ,   'QSPRN',   qsprn,                                     &
                           longname='Snow production rate (N)',               units='m-3 s-1' )

!<<SF Kasja diags

!>>DN: new diags
  CALL add_stream_element (activ,   'QMLTN2',   qmltn2,                                     &
                           longname='IC melting rate (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QXMLT',   qxmlt,                                     &
                           longname='IC (falling from above) melting rate',       units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QEVABFN',   qevabfn,                                     &
                           longname='CD evaporation rate BF (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QEVABF',  qevabf,                              &
                           longname='evaporation of CD BF',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QDEPBF',  qdepbf,                              &
                           longname='deposition rate BF',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QCDNC_DETR',   qcdnc_detr,                                     &
                           longname='cdnc detrainment rate (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QICNC_DETR',   qicnc_detr,                                     &
                           longname='icnc detrainment rate (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC1',   qcorric1,                                     &
                           longname='icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC2',   qcorric2,                                     &
                           longname='icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC3',   qcorric3,                                     &
                           longname='icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC4',   qcorric4,                                     &
                           longname='icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC5',   qcorric5,                                     &
                           longname='icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC6',   qcorric6,                                     &
                           longname='icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRIC7',   qcorric7,                                     &
                           longname='icnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD1',   qcorrcd1,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD2',   qcorrcd2,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD3',   qcorrcd3,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD4',   qcorrcd4,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD5',   qcorrcd5,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD6',   qcorrcd6,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD7',   qcorrcd7,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD8',   qcorrcd8,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD9',   qcorrcd9,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRCD10',   qcorrcd10,                                     &
                           longname='cdnc correction (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QXTTECDNC',   qxttecdnc,                                     &
                           longname='cdnc tendency (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QXTTEICNC',   qxtteicnc,                                     &
                           longname='icnc tendency (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QCORRXI',  qcorrxi,                              &
                           longname='xi correction',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QCORRXI2',  qcorrxi2,                              &
                           longname='xi correction',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QCORRXL',  qcorrxl,                              &
                           longname='xl correction',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QCORRXL2',  qcorrxl2,                              &
                           longname='xl correction',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QXLTE',  qxlte,                              &
                           longname='xl tendency',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QXITE',  qxite,                              &
                           longname='xi tendency',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSPR',       qspr,                                 &
                           longname='Snow production rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QRPR',       qrpr,                                 &
                           longname='Rain production rate',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QtestCD',   qtestCD,                                     &
                           longname='CD test (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QtestLWC',  qtestLWC,                              &
                           longname='LWC test',            units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QtestIC',   qtestIC,                                     &
                           longname='IC test (N)',                    units='m-3 s-1' )

  CALL add_stream_element (activ,   'QtestIWC',  qtestIWC,                              &
                           longname='IWC test',            units='kg kg-1 s-1' )

  CALL default_stream_setting (activ, lrerun=.TRUE.)
!<<DN: new diags

  CALL default_stream_setting (activ, laccu=.FALSE., lpost=.FALSE.)

  CALL add_stream_element (activ,   'CLOUD_COVER_DUPLIC', cloud_cover_duplic,             &
                           longname='cloud cover duplicate for record at t+1', units='1'  )

!>DN: new diagnostics
  CALL add_stream_element (activ,   'CLOUD_TM2', cloud_tm2,                               &
                           longname='cloud cover from t-2', units='1'                     )

  CALL add_stream_element (activ,   'CLOUD_BEFORE_MICRO', cloud_cover_begin_micro_2m,     &
                           longname='cloud cover before cloud microphysics', units='1'    )
!<<DN: new diagnostics

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
                           longname='Tompkins cloud cover scheme',               units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QGENTI',   qgenti,                                     &
                           longname='Tompkins cloud cover scheme',               units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QSUB_ICE',   qsub_ice,                                     &
                           longname='sublimation of IC(falling from above)',     units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QMLT_SNOW',   qmlt_snow,                                     &
                           longname='snow melting rate (total)',                    units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QMLTS_ATM',   qmlts_atm,                                     &
                           longname='snow melting rate (atmosphere)',                    units='kg kg-1 s-1' )

  CALL add_stream_element (activ,   'QMLT_CONV',   qmlt_conv,                                     &
                           longname='Melting of detrained ice',                    units='kg kg-1 s-1' )

  CALL default_stream_setting (activ, lrerun=.TRUE.)
!<<DN: new diags

  CALL default_stream_setting (activ, laccu=.FALSE.)

  CALL add_stream_element (activ,   'SICE',       sice,                                   &
                           longname='ECHAM supersaturation over ice',     units='% [0-1]' )

  END IF

!>>DN: new diagnostics
    CALL default_stream_setting (activ, laccu=.FALSE., leveltype = HYBRID)

    CALL add_stream_element (activ, 'NUC_NCD', pnucncd)
    CALL add_stream_element (activ, 'NUC_NIC', pnucnic)
    CALL add_stream_element (activ, 'EVAP_NCD', pevapncd)
    CALL add_stream_element (activ, 'EVAP_NCD_BF', pevapncdbf)
    CALL add_stream_element (activ, 'FRZ_NCD', pfrzncd)
    CALL add_stream_element (activ, 'MELT_NCD', pmeltncd)
    CALL add_stream_element (activ, 'SUB_NIC', psubnic)
    CALL add_stream_element (activ, 'CDNC_ACT_EVAP', pcdncactevap)
    CALL add_stream_element (activ, 'CDNC_ACT_FRZ', pcdncactfrz)
    CALL add_stream_element (activ, 'ICNC_ACT_MELT', picncactmelt)
    CALL add_stream_element (activ, 'ICNC_ACT_SUB', picncactsub)
    CALL add_stream_element (activ, 'DPG', dpg)
    CALL add_stream_element (activ, 'CDNC_ACT_AP', pcdncactap)
    CALL add_stream_element (activ, 'ICNC_ACT_AP', picncactap)
!<<DN: new diagnostics

!>>DN: burden
  CALL default_stream_setting (activ, lrerun=.FALSE.)

  CALL default_stream_setting (activ, laccu=.TRUE.)

  CALL add_stream_element (activ,   'DAUT',       daut,                                   &
                           longname='CD autoconversion rate',             units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DFRE',       dfre,                                   &
                           longname='CD freezing rate',                   units='m-2 s-1' )

  CALL add_stream_element (activ,   'DNUCL',       dnucl,                                 &
                           longname='CD nucleation rate',                 units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCND',       dcnd,                                   &
                           longname='condensation rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DLWC_DETR',     dlwc_detr,                           &
                           longname='xl detrainment rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DEVP_CD', devp_lwc,                                  &
                           longname='evaporation of CD',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DAUTN',       dautn,                                 &
                           longname='autoconversion rate (N)',       units='m-2 s-1' )

  CALL add_stream_element (activ,   'DRACL',       dracl,                                 &
                           longname='rain accretion rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DRACLN',      dracln,                                &
                           longname='rain accretion rate (N)',           units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSACL',      dsacl,                                  &
                           longname='CD accretion rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSACLN',      dsacln,                                  &
                           longname='CD accretion rate (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DFRZ',       dfrz,                                   &
                           longname='CD freezing rate',               units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DFRZN',       dfrzn,                                 &
                           longname='CD freezing rate (N)',               units='m-2 s-1' )

  CALL add_stream_element (activ,   'DNUCI',       dnuci,                                 &
                           longname='IC nucleation rate',                 units='m-2 s-1' )

  CALL add_stream_element (activ,   'DDEP',       ddep,                                   &
                           longname='deposition rate',                units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DIWC_DETR',     diwc_detr,                               &
                           longname='xi detrainment rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSUB_IWC', dsub_iwc,                                 &
                           longname='sublimation of IC',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DAGG',       dagg,                                   &
                           longname='snow aggregation rate',          units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DAGGN',       daggn,                                   &
                           longname='snow aggregation rate (N)',          units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSACI',       dsaci,                                 &
                           longname='IC accretion rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSACIN',       dsacin,                                 &
                           longname='IC accretion rate (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSELFN',       dselfn,                                 &
                           longname='IC self collection (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSECPROD',       dsecprod,                                 &
                           longname='Secondary IC production rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSECPRODN',       dsecprodn,                                 &
                           longname='Secondary IC production rate (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSEDI',       dsedi,                                 &
                           longname='IC sedimentation rate',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSEDIN',       dsedin,                                 &
                           longname='IC sedimentation rate (N)',              units='m-2 s-1' )

  CALL add_stream_element (activ,   'DMLT',   dmlt,                                     &
                           longname='IC melting rate',                    units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DMLTN',   dmltn,                                     &
                           longname='IC melting rate (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DEVP_RAIN', devp_rain,                               &
                           longname='evaporation of rain',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSUB_SNOW',  dsub_snow,                              &
                           longname='sublimation of snow',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DRPRN',       drprn,                                 &
                           longname='Rain production rate (N)',            units='m-2 s-1' )

  CALL add_stream_element (activ,   'DSPRN',   dsprn,                                     &
                           longname='Snow production rate (N)',               units='m-2 s-1' )

  CALL add_stream_element (activ,   'DMLTN2',   dmltn2,                                     &
                           longname='IC melting rate (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DXMLT',   dxmlt,                                     &
                           longname='IC (falling from above) melting rate',       units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DEVABFN',   devabfn,                                     &
                           longname='CD evaporation rate BF (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DEVABF',  devabf,                              &
                           longname='evaporation of CD BF',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DDEPBF',  ddepbf,                              &
                           longname='deposition rate BF',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCDNC_DETR',   dcdnc_detr,                                     &
                           longname='cdnc detrainment rate (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DICNC_DETR',   dicnc_detr,                                     &
                           longname='icnc detrainment rate (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC1',   dcorric1,                                     &
                           longname='icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC2',   dcorric2,                                     &
                           longname='icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC3',   dcorric3,                                     &
                           longname='icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC4',   dcorric4,                                     &
                           longname='icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC5',   dcorric5,                                     &
                           longname='icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC6',   dcorric6,                                     &
                           longname='icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRIC7',   dcorric7,                                     &
                           longname='icnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD1',   dcorrcd1,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD2',   dcorrcd2,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD3',   dcorrcd3,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD4',   dcorrcd4,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD5',   dcorrcd5,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD6',   dcorrcd6,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD7',   dcorrcd7,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD8',   dcorrcd8,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD9',   dcorrcd9,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRCD10',   dcorrcd10,                                     &
                           longname='cdnc correction (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DXTTECDNC',   dxttecdnc,                                     &
                           longname='cdnc tendency (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DXTTEICNC',   dxtteicnc,                                     &
                           longname='icnc tendency (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRXI',  dcorrxi,                              &
                           longname='xi correction',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRXI2',  dcorrxi2,                              &
                           longname='xi correction',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRXL',  dcorrxl,                              &
                           longname='xl correction',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCORRXL2',  dcorrxl2,                              &
                           longname='xl correction',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DXLTE',  dxlte,                              &
                           longname='xl tendency',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DXITE',  dxite,                              &
                           longname='xi tendency',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSPR',       dspr,                                 &
                           longname='Snow production rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DRPR',       drpr,                                 &
                           longname='Rain production rate',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DtestCD',   dtestCD,                                     &
                           longname='CD test (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DtestLWC',  dtestLWC,                              &
                           longname='LWC test',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DtestIC',   dtestIC,                                     &
                           longname='IC test (N)',                    units='m-2 s-1' )

  CALL add_stream_element (activ,   'DtestIWC',  dtestIWC,                              &
                           longname='IWC test',            units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DGENTL',   dgentl,                                     &
                           longname='Tompkins cloud cover scheme',               units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DGENTI',   dgenti,                                     &
                           longname='Tompkins cloud cover scheme',               units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSUB_ICE',   dsub_ice,                                     &
                           longname='sublimation of IC(falling from above)',     units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DMLT_SNOW',   dmlt_snow,                                     &
                           longname='snow melting rate (total)',                    units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DMLTS_ATM',   dmlts_atm,                                     &
                           longname='snow melting rate (atmosphere)',                    units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DMLTS_SFC',   dmlts_sfc,                                     &
                           longname='snow melting rate (surface)',                    units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSEDI_SFC',       dsedi_sfc,                                 &
                           longname='IC sedimentation rate (surface)',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DSNOW',       dsnow,                                 &
                           longname='snow fall (stratiform)',              units='kg m-2 s-1' )

  CALL add_stream_element (activ,   'DCLDTTE',       dcldtte,                                 &
                           longname='2-moment cloud microphysics energy tendency',              units='W m-2' )

  CALL add_stream_element (activ,   'DCONVTTE',       dconvtte,                                 &
                           longname='cumulus parameterization energy tendency',              units='W m-2' )

  CALL add_stream_element (activ,   'DMLT_CONV',   dmlt_conv,                                     &
                           longname='Melting of detrained ice',                    units='kg m-2 s-1' )

  CALL default_stream_setting (activ, lrerun=.TRUE.)
!<<DN: burden

  CALL default_stream_setting (activ, laccu=.FALSE.)

END SUBROUTINE construct_activ_stream

END MODULE mo_activ
