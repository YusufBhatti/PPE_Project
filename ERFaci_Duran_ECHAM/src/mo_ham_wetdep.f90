!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_wetdep.f90
!!
!! \brief
!! Module to compute HAM-specific aerosol wet deposition.
!! Originally contained in xt_wetdep.f90 and reorganized here by processes:
!! in-cloud (liq/ice phase), below-cloud (liq/ice phase)
!!
!! \author Philip Stier (MPI-Met)
!! \author Betty Croft (Dalhousie University)
!! \author Sylvaine Ferrachat (ETH Zurich)
!!
!! \responsible_coder
!! Sylvaine Ferrachat, sylvaine.ferrachat@env.ethz.ch
!!
!! \revision_history
!!   -# Johann Feicher (MPI-MET) - original code (xt_wetdep) (2001)
!!   -# Claudia Timmreck (MPI-MET) (2001)
!!   -# Philip Stier (MPI-Met) (2001-2004)
!!   -# Johann Sebastian Rast (MPI-Met) - introduction of nwetdep (2004)
!!   -# Betty Croft (Dalhousie Uni.) - size-dependent below-cloud scavenging (2005-2008)
!!   -# Sylvaine Ferrachat (ETH Zurich) - complete code cleanup and reorganisation (2009-11-04)
!!   -# Grazia Frontoso (C2SM) - porting of the in-cloud size-dep scavenging from Betty Croft (2013-06)
!!   -# Sylvaine Ferrachat (ETH Zurich) - code refactoring to optimize for mode-wise only calculations (2013-08) 
!!
!! \limitations
!! None
!!
!! \details
!! Stratiform and convective cases are handled transparently (former code duplication in 
!! xt_wetdep was removed).
!!
!! \bibliographic_references
!!    -# Croft et al, 2009, Aerosol size-dependent below-cloud scavenging by rain and snow in the ECHAM5-HAM,
!! Atmos. Chem. Phys. ; 9 ; 4653-4675
!!    -# Croft et all, 2010, Influences of in-cloud aerosol scavenging parameterizations on aerosol concentrations
!! and wet deposition in ECHAM5-HAM, Atmos. Chem. Phys. ; 10 ; 1511-1543 
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

MODULE mo_ham_wetdep

  USE mo_kind,          ONLY: dp
  USE mo_physical_constants, ONLY: tmelt
  USE mo_exception,     ONLY: finish
  USE mo_tracdef,       ONLY: ntrac, trlist, AEROSOLNUMBER, AEROSOLMASS
  USE mo_echam_cloud_params, ONLY: cthomi
  USE mo_time_control,  ONLY: time_step_len
  USE mo_ham,           ONLY: nclass
  
  IMPLICIT none

  PRIVATE

  PUBLIC :: ham_wetdep
  PUBLIC :: ham_setscav             ! to set the scavenging flags
  PUBLIC :: ham_conv_lfraq_so2      ! ++mgs: calculate liquid fraction of SO2 for conv cases
                                    ! ++mgs: former code from cuflx.f90
  PUBLIC :: prep_ham_mode_init

  !--- Constants:
  REAL(dp), PARAMETER :: zmin      = 1.e-10_dp
  REAL(dp), PARAMETER :: UNDEF     = -999._dp
  REAL(dp), PARAMETER :: zeps      = EPSILON(1._dp)
  REAL(dp), PARAMETER :: zeps_mass = 1.e-30_dp

  !--- Mode-wise scavenging coefficients and related
  INTEGER, POINTER :: indexy1(:,:,:,:) => NULL()
  INTEGER, POINTER :: indexy2(:,:,:,:) => NULL()   ! indices necessary for lookup table searches
                                                   ! shape: (kbdim,klev,numb/mass,mode)

  REAL(dp), POINTER :: mr(:,:,:,:) => NULL()       ! median radius (wet)
                                                   ! shape: (kbdim,klev, numb/mass, mode)
  REAL(dp), POINTER :: rcritrad(:,:,:,:) => NULL() ! critical radius [m] (ie min bound) for scav of mode
                                                   ! shape: (kbdim,klev, liq/ice, mode)
  REAL(dp), POINTER :: sfnuc(:,:,:,:,:) => NULL()  ! in-cloud, nucleation scavenging fraction
                                                   ! shape: (kbdim,klev, liq/ice, numb/mass, mode)
  REAL(dp), POINTER :: sfimp(:,:,:,:,:) => NULL()  ! in-cloud, impaction scavenging fraction
                                                   ! shape:  (kbdim,klev, liq/ice, numb/mass, mode)
  REAL(dp), POINTER :: sfrain(:,:,:,:)  => NULL()  ! below-cloud, scavenging fraction by rain 
                                                   ! shape: (kbdim,klev,numb/mass, mode)
  REAL(dp), POINTER :: sfsnow(:,:,:,:)  => NULL()  ! below-cloud, scavenging fraction by snow 
                                                   ! (kbdim,klev,numb/mass, mode)

  INTERFACE init_var
    MODULE PROCEDURE init_var_i_4d
    MODULE PROCEDURE init_var_r_4d
    MODULE PROCEDURE init_var_r_5d
  END INTERFACE init_var

  CONTAINS 

  !!----------------------------------------------------------------------------
  SUBROUTINE ham_wetdep(kproma, kbdim, klev, krow, ktop, kt,                          & 
                        kscavICtype, kscavBCtype,                                     & 
                        kscavICphase, kscavBCphase,                                   &
                        lstrat, ptm1, pxtm1, pxtte,                                   &  
                        pxtp10, pxtp1c, pxtp1c_sav,                                   &
                        pfrain, pfsnow, paclc, pmfu,                                  &
                        pmfuxt, prhop1, pdpg, pice,                                   &
                        peffice, peffwat, pclc, prevap,                               &
                        pdepint, pdepintbc, pdepintbcr, pdepintbcs,                   &
                        pdepintic, pdepintic_nucw, pdepintic_nucm,                    &
                        pdepintic_nucc, pdepintic_impw, pdepintic_impm, pdepintic_impc)

  ! master routine for ham scavenging calculations

  USE mo_ham_streams,   ONLY: rwet
  USE mo_ham_m7ctl,     ONLY: cmedr2mmedr
  USE mo_ham,           ONLY: nham_subm, HAM_M7
  ! YAB <
  USE mo_hammoz_perturbations, ONLY: lo_hammoz_perturbations, &
                                     scale_wetdep_bc
  !> YAB 
  !>>DN
  USE mo_hammoz_aerocom_diags, ONLY: lHEaci
  USE mo_hammoz_aerocom_HEaci, ONLY: wet3Dso4_inst
  USE mo_species, ONLY: speclist
  !<<DN

  !--- arguments
  INTEGER, INTENT(in) :: kproma         ! geographic block number of locations
  INTEGER, INTENT(in) :: kbdim          ! geographic block maximum number of locations
  INTEGER, INTENT(in) :: klev           ! numer of levels
  INTEGER, INTENT(in) :: krow           ! geographic block number
  INTEGER, INTENT(in) :: ktop           ! top layer index
  INTEGER, INTENT(in) :: kt             ! tracer index
  INTEGER, INTENT(in) :: kscavICtype    ! indicates in-cloud scavenging scheme
  INTEGER, INTENT(in) :: kscavBCtype    ! indicates below-cloud scavenging scheme
  INTEGER, INTENT(in) :: kscavICphase   ! indicates in-cloud scavenging by water and/or ice
  INTEGER, INTENT(in) :: kscavBCphase   ! indicates below-cloud scavenging by water and/or ice

  LOGICAL, INTENT(in) :: lstrat   ! stratiform or convective clouds case

  REAL(dp), INTENT(in)   :: ptm1(kbdim,klev),             & ! temperature
                            pxtm1(kbdim,klev,ntrac),      & ! tracer mixing ratio
                            pfrain(kbdim,klev),           & ! rain rate
                            pfsnow(kbdim,klev),           & ! snow rate
                            paclc(kbdim,klev),            & ! cloud cover
                            pmfu(kbdim,klev),             & ! convective flux
                            pdpg(kbdim,klev),             & ! grid box thickness
                            pice(kbdim,klev),             & ! ice fraction
                            peffice(kbdim,klev),          & ! autoconversion rate (ice)
                            peffwat(kbdim,klev),          & ! autoconversion rate (liq water)
                            pclc(kbdim,klev),             & ! fraction of grid covered by precip
                            prevap(kbdim,klev),           &
                            pxtp1c_sav(kbdim,klev,ntrac), & ! cloudy mixing ratio, untouched by wetdep
                            prhop1(kbdim,klev)              ! air density (t-dt)

  REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ntrac),  & ! tracer tendency
                             pxtp10(kbdim,klev,ntrac), & ! cloud-free mixing ratio
                             pxtp1c(kbdim,klev,ntrac), & ! cloudy mixing ratio
                             pmfuxt(kbdim,klev,ntrac), & ! updraft mmr
                             pdepint(kbdim),           & ! global scavenged mr
                             pdepintbc(kbdim),         & ! below-cloud scavenged mr
                             pdepintbcr(kbdim),        & ! below-cloud scavenged by rain mr
                             pdepintbcs(kbdim),        & ! below-cloud scavenged by snow mr
                             pdepintic(kbdim),         & ! in-cloud scavenged mr
                             pdepintic_nucw(kbdim),    & ! in-cloud by nucleation (warm cl) scav. mr
                             pdepintic_nucm(kbdim),    & ! in-cloud by nucleation (mixed-phase cl) scav. mr
                             pdepintic_nucc(kbdim),    & ! in-cloud by nucleation (cold cl) scav. mr
                             pdepintic_impw(kbdim),    & ! in-cloud by impaction (warm cl) scav. mr
                             pdepintic_impm(kbdim),    & ! in-cloud by impaction (mixed-phase cl) scav. mr
                             pdepintic_impc(kbdim)       ! in-cloud by impaction (cold cl) scav. mr

  !--- local variables 

  INTEGER  :: jk, imod, itrac_phase, itmp1(kbdim,klev), itmp2(kbdim,klev)

  LOGICAL :: ll1(kbdim,klev)
  LOGICAL :: ll_wat(kbdim,klev), ll_mxp(kbdim,klev), ll_ice(kbdim,klev)

  REAL(dp) :: ztmst, zrad_fac      ! conversion factor for median radius (mass vs number)

  REAL(dp):: zxtice(kbdim,klev),          zxtwat(kbdim,klev),          &
             zxtp1(kbdim,klev),           zdxtwat(kbdim,klev),         &
             zdxtice(kbdim,klev),         zdxtcol(kbdim,klev),         &
             zdep(kbdim,klev),                                         &
             zdxtwat_nuc(kbdim,klev),     zdxtice_nuc(kbdim,klev),     &
             zdxtwat_imp(kbdim,klev),     zdxtice_imp(kbdim,klev),     &
             zdxtcolr(kbdim,klev),        zdxtcols(kbdim,klev),        &
             zdxtevapic(kbdim,klev),      zdxtevapbc(kbdim,klev),      &
             zdxtevapic_nucw(kbdim,klev), zdxtevapic_nucm(kbdim,klev), &
             zdxtevapic_nucc(kbdim,klev),                              &
             zdxtevapic_impw(kbdim,klev), zdxtevapic_impm(kbdim,klev), &
             zdxtevapic_impc(kbdim,klev),                              &
             zdxtevapbcr(kbdim,klev),                                  &
             zdxtevapbcs(kbdim,klev),     zdep_nuc(kbdim,klev),        &
             zdep_imp(kbdim,klev),                                     &
             zxtte(kbdim,klev),           zxtp10(kbdim,klev),          &
             zmf(kbdim,klev),                                          &
             zxtfrac_col(kbdim,klev),     zxtfrac_colr(kbdim,klev),    &
             zxtfrac_cols(kbdim,klev),                                 &
             zcoeffr(kbdim,klev),         zcoeffs(kbdim,klev),         &
             ztmp1(kbdim,klev),           ztmp2(kbdim,klev)

  REAL(dp), POINTER :: rwet_p(:,:,:)

  !--- 0/ Initializations:

  ztmst = time_step_len

  imod    = trlist%ti(kt)%mode

  rwet_p => rwet(imod)%ptr

  zrad_fac = 1._dp !SF #380: proper initialization
  
  IF (trlist%ti(kt)%nphase == AEROSOLNUMBER) THEN
     itrac_phase = 1 
  ELSE IF (trlist%ti(kt)%nphase == AEROSOLMASS) THEN
     itrac_phase = 2
     !TB: count to mode median calculatin needed only for M7 (see #380)
     !    
     IF (nham_subm == HAM_M7) THEN
        zrad_fac = cmedr2mmedr(imod)
     ENDIF
  ENDIF

  zdxtevapic(1:kproma,:) = 0._dp
  zdxtevapbc(1:kproma,:) = 0._dp

  zxtfrac_col(1:kproma,:)  = 0._dp
  zxtfrac_colr(1:kproma,:) = 0._dp
  zxtfrac_cols(1:kproma,:) = 0._dp

  zdxtwat(1:kproma,:) = 0._dp
  zdxtice(1:kproma,:) = 0._dp

  zdxtwat_nuc(1:kproma,:) = 0._dp
  zdxtice_nuc(1:kproma,:) = 0._dp
  zdxtwat_imp(1:kproma,:) = 0._dp
  zdxtice_imp(1:kproma,:) = 0._dp

  IF (lstrat) THEN !stratiform case

     !--- Weight mixing ratios with cloud fraction:

     pxtp1c(1:kproma,:,kt) = pxtp1c(1:kproma,:,kt)*paclc(1:kproma,:)
     pxtp10(1:kproma,:,kt) = pxtp10(1:kproma,:,kt)*(1._dp-paclc(1:kproma,:))
     zxtp10(1:kproma,:)    = pxtp10(1:kproma,:,kt)
     zmf(1:kproma,:)       = pdpg(1:kproma,:) / ztmst
     !SF note: zxtp10 is needed in order to transparently
     !         handle the strat and the conv cases without writing too many separate,
     !         but very similar, equations
  ELSE
     zxtp10(1:kproma,:) = 0._dp
     zmf(1:kproma,:)    = pmfu(1:kproma,:)
  ENDIF

  !--- Associate tracer masses in the cloud fraction to water/ice phase
  !    to the respective mass fractions: 
  zxtwat(1:kproma,:) = pxtp1c(1:kproma,:,kt)*(1._dp-pice(1:kproma,:)) 
  zxtice(1:kproma,:) = pxtp1c(1:kproma,:,kt)*pice(1:kproma,:) 

  !--- 1/ Process:

  IF (kscavBCtype == 3 .OR. kscavICtype == 3 ) THEN !only necessary for size-dep scavenging:

     IF (indexy1(1,1,itrac_phase,imod) == UNDEF) THEN ! mode-phase calculation required!

        !--- Select aerosol wet radius and limit it to maximal 50 um:
        !
        !    If tracer is:
        !
        !       aerosol number mixing ratio: use number median radius
        !       aerosol mass   mixing ratio: use mass   median radius
        !--- Convert radius from metres to micrometres

        mr(1:kproma,:,itrac_phase,imod) = MIN(rwet_p(1:kproma,:,krow)*zrad_fac, 50.E-6_dp)
        mr(1:kproma,:,itrac_phase,imod) = mr(1:kproma,:,itrac_phase,imod)*1.E+6_dp

        ll1(1:kproma,:) = (mr(1:kproma,:,itrac_phase,imod) > zeps) !SF#294 replaced 0. by zeps

        ztmp1(1:kproma,:) = MERGE(mr(1:kproma,:,itrac_phase,imod), 1._dp, ll1(1:kproma,:)) !1. is just dummy
        ztmp2(1:kproma,:) = FLOOR((3._dp*(log(1.e4_dp*ztmp1(1:kproma,:))/log(2._dp)))+1._dp)

        itmp1(1:kproma,:) = MAX(0, MIN(60, INT(ztmp2(1:kproma,:))))
        itmp2(1:kproma,:) = MAX(0, MIN(60, INT(1._dp+ztmp2(1:kproma,:))))

        indexy1(1:kproma,:,itrac_phase,imod) = MERGE(itmp1(1:kproma,:), 0, ll1(1:kproma,:))
        indexy2(1:kproma,:,itrac_phase,imod) = MERGE(itmp2(1:kproma,:), 0, ll1(1:kproma,:))

     ENDIF !indexy1(itrac_phase,imod) == UNDEF

  ENDIF !kscavBCtype == 3 .OR. kscavICtype == 3

  !--- 1.1/ In-cloud scavenging
  IF(kscavICtype > 0) THEN

     !--- Set logical for temperature ranges (needed several times below)
     ll_wat(1:kproma,:) = (ptm1(1:kproma,:) > tmelt)
     ll_mxp(1:kproma,:) = .NOT. ll_wat(1:kproma,:) .AND. (ptm1(1:kproma,:) > cthomi)
     ll_ice(1:kproma,:) = (ptm1(1:kproma,:) <= cthomi)

     !--- 1.1.1/ Phase-specific calculations:

     IF (IAND(kscavICphase,1) /= 0) THEN !water scavenging on (kscavICphase==1 .or. 3)

        CALL ic_scav(kproma, kbdim, klev, krow, ktop, kt, &
                     imod, 1, itrac_phase, kscavICtype, &
                     lstrat, ll_wat, ll_mxp, ll_ice, &
                     prhop1, pxtp1c, pxtp1c_sav, paclc, peffwat, &
                     zdxtwat_nuc, zdxtwat_imp, zdxtwat, zxtwat)

     ENDIF !water scavenging on

     IF (IAND(kscavICphase,2) /= 0) THEN !ice scavenging on (kscavICphase==2 .or. 3)

        CALL ic_scav(kproma, kbdim, klev, krow, ktop, kt, &
                     imod, 2, itrac_phase, kscavICtype, &
                     lstrat, ll_wat, ll_mxp, ll_ice, &
                     prhop1, pxtp1c, pxtp1c_sav, paclc, peffice, &
                     zdxtice_nuc, zdxtice_imp, zdxtice, zxtice)

     ENDIF !ice scavenging on

     !--- 1.1.2/ Put everything together:

     pxtp1c(1:kproma,:,kt) = zxtwat(1:kproma,:) + zxtice(1:kproma,:)

     !--- Local deposition mass-flux [grid-box mean kg m-2 s-1]:
     zdep(1:kproma,:)     = (zdxtwat(1:kproma,:)     + zdxtice(1:kproma,:)    )*zmf(1:kproma,:)
     zdep_nuc(1:kproma,:) = (zdxtwat_nuc(1:kproma,:) + zdxtice_nuc(1:kproma,:))*zmf(1:kproma,:)
     zdep_imp(1:kproma,:) = (zdxtwat_imp(1:kproma,:) + zdxtice_imp(1:kproma,:))*zmf(1:kproma,:)

     DO jk=ktop,klev

        !--- Integrated deposition mass flux:
        pdepintic(1:kproma) = pdepintic(1:kproma) + zdep(1:kproma,jk)

        !>>SF #458 (replacing where statements)
        pdepintic_nucw(1:kproma) = pdepintic_nucw(1:kproma) &
                                 + MERGE(zdep_nuc(1:kproma,jk), 0._dp, ll_wat(1:kproma,jk)) 
        pdepintic_impw(1:kproma) = pdepintic_impw(1:kproma) &
                                 + MERGE(zdep_imp(1:kproma,jk), 0._dp, ll_wat(1:kproma,jk))

        pdepintic_nucm(1:kproma) = pdepintic_nucm(1:kproma) &
                                 + MERGE(zdep_nuc(1:kproma,jk), 0._dp, ll_mxp(1:kproma,jk)) 
        pdepintic_impm(1:kproma) = pdepintic_impm(1:kproma) &
                                 + MERGE(zdep_imp(1:kproma,jk), 0._dp, ll_mxp(1:kproma,jk))

        pdepintic_nucc(1:kproma) = pdepintic_nucc(1:kproma) &
                                 + MERGE(zdep_nuc(1:kproma,jk), 0._dp, ll_ice(1:kproma,jk)) 
        pdepintic_impc(1:kproma) = pdepintic_impc(1:kproma) &
                                 + MERGE(zdep_imp(1:kproma,jk), 0._dp, ll_ice(1:kproma,jk))
        !<<SF #458 (replacing where statements)

        !--- Re-evaporation:
        zdxtevapic(1:kproma,jk)      = pdepintic(1:kproma)     *prevap(1:kproma,jk)

        zdxtevapic_nucw(1:kproma,jk) = pdepintic_nucw(1:kproma)*prevap(1:kproma,jk)
        zdxtevapic_nucm(1:kproma,jk) = pdepintic_nucm(1:kproma)*prevap(1:kproma,jk)
        zdxtevapic_nucc(1:kproma,jk) = pdepintic_nucc(1:kproma)*prevap(1:kproma,jk)

        zdxtevapic_impw(1:kproma,jk) = pdepintic_impw(1:kproma)*prevap(1:kproma,jk)
        zdxtevapic_impm(1:kproma,jk) = pdepintic_impm(1:kproma)*prevap(1:kproma,jk)
        zdxtevapic_impc(1:kproma,jk) = pdepintic_impc(1:kproma)*prevap(1:kproma,jk)

        !--- Reduce integrated deposition mass flux by re-evap
        pdepintic(1:kproma)      = pdepintic(1:kproma)      - zdxtevapic(1:kproma,jk)

        pdepintic_nucw(1:kproma) = pdepintic_nucw(1:kproma) - zdxtevapic_nucw(1:kproma,jk)
        pdepintic_nucm(1:kproma) = pdepintic_nucm(1:kproma) - zdxtevapic_nucm(1:kproma,jk)
        pdepintic_nucc(1:kproma) = pdepintic_nucc(1:kproma) - zdxtevapic_nucc(1:kproma,jk)

        pdepintic_impw(1:kproma) = pdepintic_impw(1:kproma) - zdxtevapic_impw(1:kproma,jk)
        pdepintic_impm(1:kproma) = pdepintic_impm(1:kproma) - zdxtevapic_impm(1:kproma,jk)
        pdepintic_impc(1:kproma) = pdepintic_impc(1:kproma) - zdxtevapic_impc(1:kproma,jk)

        !>>DN
        IF(lHEaci .AND. trlist%ti(kt)%nphase == AEROSOLMASS.AND.trlist%ti(kt)%spid>0)THEN
           IF (speclist(trlist%ti(kt)%spid)%shortname == 'SO4') &
                wet3Dso4_inst(1:kproma,jk,krow)=&
                wet3Dso4_inst(1:kproma,jk,krow)+zdep(1:kproma,jk)-zdxtevapic(1:kproma,jk)
        END IF
        !<<DN        
     ENDDO

     IF (.NOT. lstrat) THEN ! conv case

        zxtp10(1:kproma,:) = -zdep(1:kproma,:)/pdpg(1:kproma,:)*ztmst
        !SF note: previously, the tendency - instead of the ambient value - was set here.
        !   tendency is indeed the only relevant quantity for later,
        !   but this has been done to establish more symetry with regard to the stratiform case scavenging,
        !   so that things can be generalized more easily (cf my note up where I set zxtp10 in the strat case)

        !--- Updraft mass flux:

        pmfuxt(1:kproma,:,kt) = pxtp1c(1:kproma,:,kt)*pmfu(1:kproma,:)

     ENDIF

  ENDIF !end in-cloud scavenging

  !--- 1.2/ Below-cloud scavenging
  IF(kscavBCtype > 0) THEN

     !--- 1.2.1/ Phase-specific calculations:

     IF (IAND(kscavBCphase,1) /= 0) THEN !rain scavenging on (kscavBCphase==1 .or. 3)

        CALL bc_rain(kproma, kbdim, klev, krow, ktop, imod, itrac_phase, kscavBCtype, pfrain)

     ENDIF !rain scavenging on

     IF (IAND(kscavBCphase,2) /= 0) THEN !snow scavenging on (kscavBCphase==2 .or. 3)

        CALL bc_snow(kproma, kbdim, klev, krow, ktop, imod, itrac_phase, kscavBCtype, pfsnow)

     ENDIF !snow scavenging on

     !--- 1.2.2/ Put everything together:

     !--- Calculate fraction of below cloud scavenged tracer:
     ll1(1:kproma,:) = (paclc(1:kproma,:) < zmin)

     ztmp1(1:kproma,:) = -ztmst*MAX(sfrain(1:kproma,:,itrac_phase,imod),0._dp)
     ztmp2(1:kproma,:) = -ztmst*MAX(sfsnow(1:kproma,:,itrac_phase,imod),0._dp)
!SFnote: in the above two expressions, the MAX function is here only to rule out the cases where sfrain and/or
!        sfsnow is/are equal to UNDEF, ie in cases where scavenging by rain and/or snow is not relevant.
!        Initializing sfrain and sfsnow to 0 at the beginning would defeat the concept of having an UNDEF value,
!        which is necessary to flag cases where mode-wise (and trace phase-wise) calculations are necessary.
!        This concept allows to perform these potentially expensive calculations only once per mode 
!        and tracer phase.

     !>>SF #458 (replacing where statements)
     zcoeffr(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:))
     zcoeffs(1:kproma,:) = MERGE(ztmp2(1:kproma,:), 0._dp, ll1(1:kproma,:))

     zxtfrac_colr(1:kproma,:) = 1._dp - EXP(zcoeffr(1:kproma,:)                    )
     zxtfrac_cols(1:kproma,:) = 1._dp - EXP(                    zcoeffs(1:kproma,:))
     zxtfrac_col(1:kproma,:)  = 1._dp - EXP(zcoeffr(1:kproma,:)+zcoeffs(1:kproma,:))

     zxtfrac_colr(1:kproma,:) = MAX(0._dp, MIN(1._dp, zxtfrac_colr(1:kproma,:) ) ) 
     zxtfrac_cols(1:kproma,:) = MAX(0._dp, MIN(1._dp, zxtfrac_cols(1:kproma,:) ) ) 
     zxtfrac_col(1:kproma,:)  = MAX(0._dp, MIN(1._dp, zxtfrac_col(1:kproma,:) ) ) 

     ll1(1:kproma,:) = (pclc(1:kproma,:) > zmin)

     ztmp1(1:kproma,:) = pxtp10(1:kproma,:,kt)*pclc(1:kproma,:)*zxtfrac_col(1:kproma,:)
     zdxtcol(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:))

     ztmp1(1:kproma,:)    = pxtp10(1:kproma,:,kt)*pclc(1:kproma,:)*zxtfrac_colr(1:kproma,:)*pdpg(1:kproma,:)/ztmst
     zdxtcolr(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:))

     ztmp1(1:kproma,:)    = pxtp10(1:kproma,:,kt)*pclc(1:kproma,:)*zxtfrac_cols(1:kproma,:)*pdpg(1:kproma,:)/ztmst
     zdxtcols(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:))
     !<<SF #458 (replacing where statements)

     
     !>>YAB Adding perturbed physics scaling for below-cloud scavenging:
     IF (lo_hammoz_perturbations) THEN 
        ! Scale the BC separately
!        IF (trlist%ti(kt)%basename == "BC") THEN
!           zdxtcol(1:kproma,:) = zdxtcol(1:kproma,:) * scale_wetdep_bc_bc_only
!           zdxtcolr(1:kproma,:) = zdxtcolr(1:kproma,:) * scale_wetdep_bc_bc_only
!           zdxtcols(1:kproma,:) = zdxtcols(1:kproma,:) * scale_wetdep_bc_bc_only
!        ENDIF

        ! Scale everything BC
        zdxtcol(1:kproma,:) = zdxtcol(1:kproma,:) * scale_wetdep_bc
        zdxtcolr(1:kproma,:) = zdxtcolr(1:kproma,:) * scale_wetdep_bc
        zdxtcols(1:kproma,:) = zdxtcols(1:kproma,:) * scale_wetdep_bc
     ENDIF
     !<<YAB

     zxtp10(1:kproma,:) = zxtp10(1:kproma,:) - zdxtcol(1:kproma,:)
     zdxtcol(1:kproma,:) = zdxtcol(1:kproma,:)*pdpg(1:kproma,:)/ztmst

     DO jk=ktop,klev

        pdepintbc(1:kproma)  = pdepintbc(1:kproma)  + zdxtcol(1:kproma,jk)
        pdepintbcr(1:kproma) = pdepintbcr(1:kproma) + zdxtcolr(1:kproma,jk)
        pdepintbcs(1:kproma) = pdepintbcs(1:kproma) + zdxtcols(1:kproma,jk)

        !--- Re-evaporation:
        zdxtevapbc(1:kproma,jk)  = pdepintbc(1:kproma) *prevap(1:kproma,jk)
        zdxtevapbcr(1:kproma,jk) = pdepintbcr(1:kproma)*prevap(1:kproma,jk)
        zdxtevapbcs(1:kproma,jk) = pdepintbcs(1:kproma)*prevap(1:kproma,jk)

        !--- Reduce integrated deposition mass flux by re-evap
        pdepintbc(1:kproma)  = pdepintbc(1:kproma)  - zdxtevapbc(1:kproma,jk)
        pdepintbcr(1:kproma) = pdepintbcr(1:kproma) - zdxtevapbcr(1:kproma,jk)
        pdepintbcs(1:kproma) = pdepintbcs(1:kproma) - zdxtevapbcs(1:kproma,jk)
        !>>DN
        IF(lHEaci .AND. trlist%ti(kt)%nphase == AEROSOLMASS)THEN
           IF (speclist(trlist%ti(kt)%spid)%shortname == 'SO4') &
                wet3Dso4_inst(1:kproma,jk,krow)=&
                wet3Dso4_inst(1:kproma,jk,krow)+zdxtcol(1:kproma,jk)-zdxtevapbc(1:kproma,jk)
        END IF
        !<<DN        
     ENDDO !jk

  ENDIF !end below-cloud scavenging

  !--- put together all contributions to integrated deposition mass flux
  pdepint(1:kproma) = pdepintic(1:kproma) + pdepintbc(1:kproma)

  !--- Total tendency + update updraft flux in conv case:   
  IF (lstrat) THEN !SF stratiform case

     pxtp10(1:kproma,:,kt) = zxtp10(1:kproma,:)                              &
                           + (zdxtevapic(1:kproma,:)+zdxtevapbc(1:kproma,:)) &
                             / pdpg(1:kproma,:) * ztmst

     zxtp1(1:kproma,:) = pxtm1(1:kproma,:,kt)+pxtte(1:kproma,:,kt)*ztmst
     zxtte(1:kproma,:) = (pxtp10(1:kproma,:,kt)+pxtp1c(1:kproma,:,kt)-zxtp1(1:kproma,:)) / ztmst

  ELSE !SF conv case

     zxtte(1:kproma,:) = zxtp10(1:kproma,:) / ztmst                                           &
                       + (zdxtevapic(1:kproma,:) + zdxtevapbc(1:kproma,:)) / pdpg(1:kproma,:)

  ENDIF

  pxtte(1:kproma,:,kt) = pxtte(1:kproma,:,kt) + zxtte(1:kproma,:)

  END SUBROUTINE ham_wetdep
  
  !!----------------------------------------------------------------------------
  SUBROUTINE ic_scav(kproma, kbdim, klev, krow, ktop, kt,           & !in
                     kmod, kwat_phase, ktrac_phase, kscavICtype,    & !in
                     lstrat, ld_wat, ld_mxp, ld_ice,                & !in
                     prhop1, pxtp1c, pxtp1c_sav, paclc, peff,       & !in
                     pdxt_nuc, pdxt_imp, pdxt, pxt)                   !out

  ! In-cloud scavenging master routine

    !>>YAB  Added perturbed physics setup:
    USE mo_hammoz_perturbations, ONLY: lo_hammoz_perturbations, &
                                     scale_wetdep_ic 
    !<< YAB

    INTEGER, INTENT(in)   :: kproma, kbdim, klev, krow, ktop, kt
    INTEGER, INTENT(in)   :: kmod                         ! current tracer mode
    INTEGER, INTENT(in)   :: kwat_phase                   ! kwat_phase=1 --> water; kwat_phase=2 --> ice 
    INTEGER, INTENT(in)   :: ktrac_phase                  ! ktrac_phase=1 --> number; kwat_phase=2 --> mass
    INTEGER, INTENT(in)   :: kscavICtype                  ! in-cloud scavenging scheme

    LOGICAL, INTENT(in)   :: lstrat                       ! flag to stratiform or convective clouds 
    LOGICAL, INTENT(in)   :: ld_wat(kbdim,klev)           ! liq water: ptm1 > tmelt
    LOGICAL, INTENT(in)   :: ld_mxp(kbdim,klev)           ! mixed-phase clouds: cthomi < ptm1 <= tmelt
    LOGICAL, INTENT(in)   :: ld_ice(kbdim,klev)           ! ice: ptm1 <= cthomi

    REAL(dp), INTENT(in)  :: pxtp1c(kbdim,klev,ntrac)     ! in-cloud tracer concentration
    REAL(dp), INTENT(in)  :: pxtp1c_sav(kbdim,klev,ntrac) ! in-cloud tracer concentration, untouched by wetdep
    REAL(dp), INTENT(in)  :: prhop1(kbdim,klev)           ! air density (t-dt)
    REAL(dp), INTENT(in)  :: paclc(kbdim,klev)            ! cloud cover
    REAL(dp), INTENT(in)  :: peff(kbdim,klev)             ! autoconversion rate (liq or ice)

    REAL(dp), INTENT(out) :: pdxt_nuc(kbdim,klev)      ! change in tracer mass assoc. with nucleation scav 
                                                       ! (for relevant phase)
    REAL(dp), INTENT(out) :: pdxt_imp(kbdim,klev)      ! change in tracer mass assoc. with impaction scav 
                                                       ! (for relevant phase)
    REAL(dp), INTENT(out) :: pdxt(kbdim,klev)          ! change in tracer mass (total)
                                                       ! (for relevant phase)
    REAL(dp), INTENT(inout) :: pxt(kbdim,klev)         ! tracer mass (total)
                                                       ! (for relevant phase)

    !local vars:
    LOGICAL :: ll1(1:kproma,klev)

    REAL(dp) :: zxtfrac(kbdim,klev), zxtfrac_nuc(kbdim,klev), zxtfrac_imp(kbdim,klev)
    REAL(dp) :: ztmp1(1:kbdim,klev)

    !--- Get the proper scavenging fractions
    CALL get_icscavfrac(kproma, kbdim, klev, krow, ktop, kt,                &
                        kmod, kwat_phase, ktrac_phase, kscavICtype, lstrat, &
                        ld_wat, ld_mxp, ld_ice,                             &
                        prhop1, pxtp1c, pxtp1c_sav,                         &
                        zxtfrac, zxtfrac_nuc, zxtfrac_imp)  
  
    !--- Change in in-cloud (strat) or updraft (conv) tracer concentration:
    !>>SF #458 (replacing where statements)
    ll1(1:kproma,:) = (paclc(1:kproma,:) > zmin)

    ztmp1(1:kproma,:)    = pxt(1:kproma,:)*zxtfrac_nuc(1:kproma,:)*peff(1:kproma,:)
    pdxt_nuc(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:))

    ztmp1(1:kproma,:)    = pxt(1:kproma,:)*zxtfrac_imp(1:kproma,:)*peff(1:kproma,:)
    pdxt_imp(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:))

    ztmp1(1:kproma,:) = pxt(1:kproma,:)*zxtfrac(1:kproma,:)*peff(1:kproma,:)
    pdxt(1:kproma,:)  = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:))
    !<<SF #458 (replacing where statements)

    !>>YAB Adding perturbed physics scaling for in-cloud scavenging:
    IF (lo_hammoz_perturbations) THEN 
       ! Scale the tendency rather than the fraction


       ! Scale everything (including BC)
       pdxt(1:kproma,:) = pdxt(1:kproma,:) * scale_wetdep_ic
       pdxt_imp(1:kproma,:) = pdxt_imp(1:kproma,:) * scale_wetdep_ic
       pdxt_nuc(1:kproma,:) = pdxt_nuc(1:kproma,:) * scale_wetdep_ic
    ENDIF
    !<<YAB

    pxt(1:kproma,:) = pxt(1:kproma,:) - pdxt(1:kproma,:)

  END SUBROUTINE ic_scav

  !!----------------------------------------------------------------------------
  SUBROUTINE get_icscavfrac(kproma, kbdim, klev, krow, ktop, kt, kmod, & ! in
                            kwat_phase, ktrac_phase, kscavICtype,      & ! in
                            lstrat, ld_wat, ld_mxp, ld_ice,            & ! in
                            prhop1, pxtp1c, pxtp1c_sav,                & ! in
                            pfrac, pfrac_nuc, pfrac_imp)                 ! out

! Utility routine to compute the in-cloud scavenging fractions

    USE mo_ham_wetdep_data, ONLY: csr_strat_wat, csr_strat_mix, csr_strat_ice, &
                                  csr_conv
    USE mo_ham,             ONLY: nham_subm, HAM_BULK, HAM_M7, HAM_SALSA

    INTEGER, INTENT(in)   :: kproma, kbdim, klev, krow, ktop, kt, kmod
    INTEGER, INTENT(in)   :: kwat_phase  ! kwat_phase=1 --> water; kwat_phase=2 --> ice 
    INTEGER, INTENT(in)   :: ktrac_phase ! ktrac_phase=1 --> number; kwat_phase=2 --> mass 
    INTEGER, INTENT(in)   :: kscavICtype

    LOGICAL, INTENT(in)   :: lstrat             ! flag to stratiform or convective clouds 
    LOGICAL, INTENT(in)   :: ld_wat(kbdim,klev) ! liq water: ptm1 > tmelt
    LOGICAL, INTENT(in)   :: ld_mxp(kbdim,klev) ! mixed-phase clouds: cthomi < ptm1 <= tmelt
    LOGICAL, INTENT(in)   :: ld_ice(kbdim,klev) ! ice: ptm1 <= cthomi

    REAL(dp), INTENT(in)  :: pxtp1c(kbdim,klev,ntrac)     ! in-cloud tracer concentration multiplied by cloud fraction
    REAL(dp), INTENT(in)  :: pxtp1c_sav(kbdim,klev,ntrac) ! in-cloud tracer concentration multiplied by cloud fraction
                                                          !SFnote: same as pxtp1c, but untouched by wetdep.
                                                          !        this is necessary for nucleation scavenging 
    REAL(dp), INTENT(in)  :: prhop1(kbdim,klev)           ! air density (t-dt)

    REAL(dp), INTENT(out) :: pfrac(kbdim,klev),     &  ! total      scavenging fraction
                             pfrac_nuc(kbdim,klev), &  ! nucleation scavenging fraction
                             pfrac_imp(kbdim,klev)     ! impaction  scavenging fraction

    REAL(dp) :: zcoeff_warm, zcoeff_mix, zcoeff_cold

    LOGICAL :: ll1(kbdim,klev)

    pfrac_imp(1:kproma,:) = 0._dp

    SELECT CASE(kscavICtype)

      !CASE(1) ! prescribed fractions !not implemented here

      CASE(2)  ! mode-wise fractions

         !--- Nucleation scavenging params:

         IF (lstrat) THEN !SF stratiform case

            IF (kwat_phase == 1) THEN
               !water
               zcoeff_warm  = csr_strat_wat(kmod)
               zcoeff_mix   = csr_strat_mix(kmod)
               zcoeff_cold  = csr_strat_wat(kmod)
            ELSE
               !ice
               zcoeff_warm  = csr_strat_wat(kmod)
               zcoeff_mix   = csr_strat_mix(kmod)
               zcoeff_cold  = csr_strat_ice(kmod)
            ENDIF

            !>>SF #458 (replacing where statements)
            pfrac_nuc(1:kproma,:) = MERGE(zcoeff_warm, 0._dp                , ld_wat(1:kproma,:))
            pfrac_nuc(1:kproma,:) = MERGE(zcoeff_mix,  pfrac_nuc(1:kproma,:), ld_mxp(1:kproma,:))
            pfrac_nuc(1:kproma,:) = MERGE(zcoeff_cold, pfrac_nuc(1:kproma,:), ld_ice(1:kproma,:))
            !<<SF #458 (replacing where statements)

         ELSE !SF convective clouds case

            pfrac_nuc(1:kproma,:) = csr_conv(kmod)

         ENDIF

      CASE(3) ! aerosol size-dep fraction

         SELECT CASE(nham_subm)
             CASE(HAM_BULK)
                ! not implemented
             CASE(HAM_M7)
         
                 IF (lstrat) THEN !GF stratiform case
         
                    ! Nucleation and impaction scavenging:
                    IF (sfnuc(1,1,kwat_phase,ktrac_phase,kmod) == UNDEF) THEN ! mode-phase calculation required!
                       CALL ic_scav_nuc(kproma, kbdim, klev, krow, kwat_phase, &
                                        ktrac_phase, kmod, prhop1, pxtp1c_sav)
                    ENDIF !sfnuc(1,1,kwat_phase,ktrac_phase,kmod) == UNDEF
         
                    pfrac_nuc(1:kproma,:) = sfnuc(1:kproma,:,kwat_phase,ktrac_phase,kmod)
                    
                    ! Impaction scavenging:
                    IF (sfimp(1,1,kwat_phase,ktrac_phase,kmod) == UNDEF) THEN ! mode-phase calculation required!
                       CALL ic_scav_imp(kproma, kbdim, klev, krow, ktop, kwat_phase, &
                                        ktrac_phase, kmod, prhop1, pxtp1c)
                    ENDIF !sfimp(1,1,kwat_phase,ktrac_phase,kmod) == UNDEF
         
                    ll1(1:kproma,:) = (pxtp1c(1:kproma,:,kt) > zeps) !SF#294 replaced 0. by zeps
                    pfrac_imp(1:kproma,:) = MERGE(sfimp(1:kproma,:,kwat_phase,ktrac_phase,kmod), 0._dp, ll1(1:kproma,:))
         
                 ELSE !GF convective case
         
                    pfrac_nuc(1:kproma,:) = csr_conv(kmod)  !no size-dep nuc scavenging for convective clouds!
         
                 ENDIF
             CASE(HAM_SALSA)
                ! Not implemented. This option will halt the model in mo_ham_init
         END SELECT
     
      CASE default

         CALL finish('get_icscavfrac','wrong kscavICtype value')

    END SELECT !kscavICtype

    !--- Calculate the scavenging parameter as sum over the processes:
    pfrac(1:kproma,:) = pfrac_nuc(1:kproma,:) + pfrac_imp(1:kproma,:)

    !--- Confine the fraction between 0% and 100% :
    pfrac(1:kproma,:)     = MAX(0._dp, MIN(1._dp, pfrac(1:kproma,:)))
    pfrac_nuc(1:kproma,:) = MAX(0._dp, MIN(1._dp, pfrac_nuc(1:kproma,:)))
    pfrac_imp(1:kproma,:) = MAX(0._dp, MIN(1._dp, pfrac_imp(1:kproma,:)))

  END SUBROUTINE get_icscavfrac

  !! ---------------------------------------------------------------------------------------
  SUBROUTINE ic_scav_nuc(kproma, kbdim, klev, krow, kwat_phase, &
                         ktrac_phase, kmod, prhop1, pxtp1c_sav)

! Grazia Frontoso, C2SM-ETHZ, 2013 - compute in-cloud size dependent 
!                                    nucleation scavenging coefficients
!                                    from Croft et al. 2010
!
!SF: output of this subroutine: sfnuc
!    sfnuc is defined in the whole module instead of being passed as an intent(out) because
!    it needs to be kept from one instance to the next

    USE mo_ham_m7_trac,          ONLY: idt_nks, idt_nas, idt_ncs
    USE mo_activ,                ONLY: na, idt_cdnc, idt_icnc
    USE mo_ham_streams,          ONLY: frac, rdry, rwet
    USE mo_ham_tools,            ONLY: ham_m7_logtail, ham_m7_invertlogtail
    USE mo_param_switches,       ONLY: ncd_activ
    USE mo_tracdef,              ONLY: ntrac

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow
    INTEGER, INTENT(in) :: kwat_phase  ! kwat_phase=1 --> water; kwat_phase=2 --> ice 
    INTEGER, INTENT(in) :: ktrac_phase ! ktrac_phase=1 --> number; ktrac_phase=2 --> mass
    INTEGER, INTENT(in) :: kmod        ! mode index
   
    REAL(dp), INTENT(in)  :: pxtp1c_sav(kbdim,klev,ntrac) ! in-cloud tracer concentration as untouched by wetdep 
    REAL(dp), INTENT(in)  :: prhop1(kbdim,klev)           ! air density (t-dt)

    ! Local variables
    REAL(dp) :: zxie(kbdim,klev),         & ! factor for inverse error function calculation
                zxtp1c(kbdim,klev,ntrac), & ! in-cloud tracer concentration
                ztmp1(kbdim,klev)

    REAL(dp), POINTER :: zrad_p(:,:) ! Pointer to dry or wet radius field,
                                     ! as appropriate to activation scheme

    LOGICAL :: ll_trac_phase, ll1(kbdim,klev)

    ! Initialization

    !SF note: the following is M7-dependent! It could be potentially generalized
    !         by replacing this condition by a true size condition + some generalized way
    !         of distinguishing mixed modes from insoluble ones
    !         For now, this is kept so, and size-dep nucleation scavenging is made unuseable
    !         if an alternate aerosol microphysics scheme is used (e.g. SALSA)
    IF (kmod < 2 .OR. kmod > 4) THEN !do nothing for non-relevant modes
       sfnuc(1:kproma,:,kwat_phase,ktrac_phase,kmod) = 0._dp
       RETURN
    ENDIF

    IF (ncd_activ == 2) THEN
      zrad_p => rdry(kmod)%ptr(:,:,krow) ! ARG activation is based on dry radius
    ELSE
      zrad_p => rwet(kmod)%ptr(:,:,krow) ! Lin & Leaitch activation is based on wet radius
    END IF

    ! Critical radius for scavenged mode.
    ! This calculation is only needed *once per mode* because it is not 
    ! dependent of whether this is a mass or number calculation
    ! (ie not dependent on ktrac_phase)

    IF (rcritrad(1,1,kwat_phase,kmod) == UNDEF) THEN ! rcritrad calculation required
       zxtp1c(1:kproma,:,:) = pxtp1c_sav(1:kproma,:,:)
       zxtp1c(1:kproma,:,:) = MAX(0._dp,zxtp1c(1:kproma,:,:)) !SFnote: this should be already ensured elsewhere!! 
                                                              !SF--> to cleanup

       SELECT CASE(kwat_phase)
          CASE(1) !liq water
      
             ll1(1:kproma,:) = (zxtp1c(1:kproma,:,idt_cdnc) > zeps_mass) .AND. &
                               (na(1:kproma,:,krow) > zeps)
          
             ztmp1(1:kproma,:) = zxtp1c(1:kproma,:,idt_cdnc) * prhop1(1:kproma,:)                & 
                               * frac(kmod)%ptr(1:kproma,:,krow) / MAX(na(1:kproma,:,krow),zeps)
         
          CASE(2) !ice
   
             ll1(1:kproma,:) = (zxtp1c(1:kproma,:,idt_icnc) > zeps_mass)
   
!>>SF to refactor! (M7-dependency)
             IF (kmod == 4) ztmp1(1:kproma,:) = MIN(1._dp, &
                                                    zxtp1c(1:kproma,:,idt_icnc) / (zxtp1c(1:kproma,:,idt_ncs)+zeps))
             IF (kmod == 3) ztmp1(1:kproma,:) = MIN(1._dp, &
                                                    MAX(0._dp, &
                                                        zxtp1c(1:kproma,:,idt_icnc)-zxtp1c(1:kproma,:,idt_ncs)) &
                                                    /(zxtp1c(1:kproma,:,idt_nas)+zeps))
             IF (kmod == 2) ztmp1(1:kproma,:) = MIN(1._dp, &
                                                    MAX(0._dp, &
                                                        zxtp1c(1:kproma,:,idt_icnc)-zxtp1c(1:kproma,:,idt_ncs) &
                                                                                   -zxtp1c(1:kproma,:,idt_nas)) &
                                                    /(zxtp1c(1:kproma,:,idt_nks)+zeps))
!<<SF to refactor!
   
       END SELECT !kwat_phase
   
       ztmp1(1:kproma,:) = 1._dp - 2._dp*MAX(0._dp,MIN(1._dp,ztmp1(1:kproma,:))) 
       zxie(1:kproma,:) = MERGE(ztmp1(1:kproma,:),1._dp,ll1(1:kproma,:)) 
   
       CALL ham_m7_invertlogtail(kproma, kbdim, klev, krow, kmod, &
                                 zrad_p(:,:), zxie(:,:),          &
                                 rcritrad(:,:,kwat_phase,kmod) )
 
    ENDIF !rcritrad calculation required

    ! Final scavenged fraction calculation:

    ll_trac_phase = (ktrac_phase == 1) ! true when tracer is a number, false when it is a mass

    CALL ham_m7_logtail(kproma, kbdim, klev, krow, kmod,      &
                        ll_trac_phase, zrad_p(:,:),           &
                        rcritrad(:,:,kwat_phase,kmod),        &
                        sfnuc(:,:,kwat_phase,ktrac_phase,kmod))

  END SUBROUTINE ic_scav_nuc

  !!----------------------------------------------------------------------------
  SUBROUTINE ic_scav_imp(kproma, kbdim, klev, krow, ktop, kwat_phase, &
                         ktrac_phase, kmod, prhop1, pxtp1c)

! Grazia Frontoso, C2SM-ETHZ, 2013 - compute in-cloud size dependent 
!                                    impaction scavenging coefficients
!                                    from Croft et al. 2010
!SF: output of this subroutine: sfimp
!    sfimp is defined in the whole module instead of being passed as an intent(out) because
!    it needs to be kept from one instance to the next

    USE mo_ham_wetdep_data,   ONLY: cdroprad, caerorad,   &
                                    cplaterad,            &
                                    scavdropm, scavdropn, &
                                    scaviceplate
 
    USE mo_activ,                ONLY: idt_icnc, reffl, reffi
    USE mo_ham_tools,            ONLY: scavcoef_bilinterp

    INTEGER, INTENT(in)     :: kproma, kbdim, klev, krow, ktop
    INTEGER, INTENT(in)     :: kwat_phase  ! kwat_phase=1 --> water; kwat_phase=2 --> ice
    INTEGER, INTENT(in)     :: ktrac_phase ! ktrac_phase=1 --> number; ktrac_phase=2 --> mass
    INTEGER, INTENT(in)     :: kmod

    REAL(dp), INTENT(in)    :: pxtp1c(kbdim,klev,ntrac) ! in-cloud tracer mass
    REAL(dp), INTENT(in)    :: prhop1(kbdim,klev)       ! air density (t-dt)
  
    ! Local variables

    INTEGER  :: jk, jl
    INTEGER  :: indexdropx1(kbdim,klev), indexdropx2(kbdim,klev),   &
                indexplatex1(kbdim,klev), indexplatex2(kbdim,klev), &  
                itmp1(kbdim,klev), itmp2(kbdim,klev)

    REAL(dp) :: X1(kbdim,klev), X2(kbdim,klev),                 &
                Y1(kbdim,klev), Y2(kbdim,klev) ,                &
                Q11(kbdim,klev), Q12(kbdim,klev) ,              &
                Q21(kbdim,klev), Q22(kbdim,klev)

    REAL(dp) :: zscavcoefplate(kbdim,klev), & ! scavenging coefficient for plate (ice)
                zscavcoefdrop(kbdim,klev),  & ! scavenging coefficient for cloud droplets
                zicnc(kbdim,klev)             ! ice crystal number concentration [# m-3]

    LOGICAL :: ll1(kbdim,klev), ll2(kbdim,klev)

    ! Initialization
    Q11(1:kproma,:) = 0._dp
    Q12(1:kproma,:) = 0._dp
    Q21(1:kproma,:) = 0._dp
    Q22(1:kproma,:) = 0._dp
    X1(1:kproma,:)  = 0._dp
    X2(1:kproma,:)  = 0._dp
    Y1(1:kproma,:)  = 0._dp
    Y2(1:kproma,:)  = 0._dp

    zicnc(1:kproma,:) = pxtp1c(1:kproma,:,idt_icnc)*prhop1(1:kproma,:)
               
    ! Calculate collision rates
    SELECT CASE(kwat_phase) 

       CASE(1) !liq

          ll1(1:kproma,:) = (reffl(1:kproma,:,krow) > zeps) !SF#294 replaced 0. by zeps

          itmp1(1:kproma,:) = MAX(0, MIN(9,INT(FLOOR(reffl(1:kproma,:,krow)/5._dp))))
          itmp2(1:kproma,:) = MAX(0, MIN(9,INT(1._dp + FLOOR(reffl(1:kproma,:,krow)/5._dp))))

          indexdropx1(1:kproma,:) = MERGE(itmp1(1:kproma,:), 0, ll1(1:kproma,:))
          indexdropx2(1:kproma,:) = MERGE(itmp2(1:kproma,:), 0, ll1(1:kproma,:))

          ! Assign in-cloud impaction scavenging coefficients (droplets)
          SELECT CASE(ktrac_phase) 
             CASE(1) ! number
                 DO jk=ktop,klev
                    DO jl=1,kproma
                       Q11(jl,jk) = scavdropn(indexdropx1(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                       Q12(jl,jk) = scavdropn(indexdropx2(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                       Q21(jl,jk) = scavdropn(indexdropx1(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
                       Q22(jl,jk) = scavdropn(indexdropx2(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
                    ENDDO
                 ENDDO
             CASE(2) ! mass
                 DO jk=ktop,klev
                    DO jl=1,kproma
                       Q11(jl,jk) = scavdropm(indexdropx1(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                       Q12(jl,jk) = scavdropm(indexdropx2(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                       Q21(jl,jk) = scavdropm(indexdropx1(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
                       Q22(jl,jk) = scavdropm(indexdropx2(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
                    ENDDO
                 ENDDO
          END SELECT
  
          DO jk=ktop,klev
             DO jl=1,kproma
                X1(jl,jk) = cdroprad(indexdropx1(jl,jk))
                X2(jl,jk) = cdroprad(indexdropx2(jl,jk))
                Y1(jl,jk) = caerorad(indexy1(jl,jk,ktrac_phase,kmod))
                Y2(jl,jk) = caerorad(indexy2(jl,jk,ktrac_phase,kmod))
             ENDDO
          ENDDO
              
          CALL scavcoef_bilinterp (kproma, kbdim, klev, krow,  ktop,          &
                                   reffl(:,:,krow), mr(:,:,ktrac_phase,kmod), &
                                   X1, X2, Y1, Y2,                            &
                                   Q11, Q12, Q21, Q22,                        &
                                   zscavcoefdrop)
  
          sfimp(1:kproma,:,kwat_phase,ktrac_phase,kmod) = zscavcoefdrop(1:kproma,:)

       CASE(2) !ice

          ll1(1:kproma,:) = (zicnc(1:kproma,:) >= zeps)

          ll2(1:kproma,:) = ll1(1:kproma,:)                   .AND. &
                            (reffi(1:kproma,:,krow) < 50._dp) .AND. &
                            (reffi(1:kproma,:,krow) >= 1._dp)

          itmp1(1:kproma,:) = MAX(0, MIN(10,INT(FLOOR(reffi(1:kproma,:,krow)/5._dp))))
          itmp2(1:kproma,:) = MAX(0, MIN(10,INT(1._dp + FLOOR(reffi(1:kproma,:,krow)/5._dp))))

          indexplatex1(1:kproma,:) = MERGE(itmp1(1:kproma,:), 0, ll2(1:kproma,:))
          indexplatex2(1:kproma,:) = MERGE(itmp2(1:kproma,:), 0, ll2(1:kproma,:))

          ll2(1:kproma,:) = ll1(1:kproma,:)                    .AND. &
                            (reffi(1:kproma,:,krow) >= 50._dp)

          itmp1(1:kproma,:) = MAX(0, MIN(34,INT(8._dp+FLOOR(reffi(1:kproma,:,krow)/50._dp))))
          itmp2(1:kproma,:) = MAX(0, MIN(34,INT(9._dp + FLOOR(reffi(1:kproma,:,krow)/50._dp))))

          indexplatex1(1:kproma,:) = MERGE(itmp1(1:kproma,:), indexplatex1(1:kproma,:), ll2(1:kproma,:))
          indexplatex2(1:kproma,:) = MERGE(itmp2(1:kproma,:), indexplatex2(1:kproma,:), ll2(1:kproma,:))

          ! Assign in-cloud impaction scavenging coefficients (ice plates)
          DO jk=ktop,klev
             DO jl=1,kproma
                Q11(jl,jk) = scaviceplate(indexplatex1(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                Q12(jl,jk) = scaviceplate(indexplatex2(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                Q21(jl,jk) = scaviceplate(indexplatex1(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
                Q22(jl,jk) = scaviceplate(indexplatex2(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
             ENDDO
          ENDDO

          DO jk=ktop,klev
             DO jl=1,kproma
                X1(jl,jk) = cplaterad(indexplatex1(jl,jk))
                X2(jl,jk) = cplaterad(indexplatex2(jl,jk))
                Y1(jl,jk) = caerorad(indexy1(jl,jk,ktrac_phase,kmod))
                Y2(jl,jk) = caerorad(indexy2(jl,jk,ktrac_phase,kmod))
             ENDDO
          ENDDO

          CALL scavcoef_bilinterp (kproma, kbdim,  klev,                      &
                                   krow, ktop,                                &
                                   reffi(:,:,krow), mr(:,:,ktrac_phase,kmod), &
                                   X1, X2, Y1, Y2,                            &
                                   Q11, Q12, Q21, Q22,                        &
                                   zscavcoefplate)          

          sfimp(1:kproma,:,kwat_phase,ktrac_phase,kmod) = 1._dp - EXP(-zscavcoefplate(1:kproma,:)*1.e-6_dp     &
                                                               *zicnc(1:kproma,:)*time_step_len)
          
    END SELECT !kwat_phase

  END SUBROUTINE ic_scav_imp

  !!----------------------------------------------------------------------------
  SUBROUTINE bc_rain(kproma, kbdim, klev, krow, ktop, kmod, & !in
                     ktrac_phase, kscavBCtype,              & !in
                     pfrain)                                  !in

  ! Below-cloud rain scavenging routine


!SF: output of this subroutine: sfrain
!    sfrain is defined in the whole module instead of being passed as an intent(out) because
!    it needs to be kept from one instance to the next

  USE mo_ham_tools,       ONLY: scavcoef_bilinterp
  USE mo_ham_wetdep_data, ONLY: cbcr, cscavbcrn, cscavbcrm, crainrate, caerorad

  INTEGER, INTENT(in)    :: kproma, kbdim, klev, krow, ktop, kmod, kscavBCtype, ktrac_phase
  REAL(dp), INTENT(in)   :: pfrain(kbdim,klev)

  !local vars:
  INTEGER  :: jk, jl, indexbcrx1(kbdim,klev), indexbcrx2(kbdim,klev)

  LOGICAL :: ll1(1:kproma, klev)

  REAL(dp) :: ztmp1(kbdim,klev), &
              X1(kbdim,klev),   X2(kbdim,klev),   &
              Y1(kbdim,klev),   Y2(kbdim,klev) ,  &
              Q11(kbdim,klev),  Q12(kbdim,klev) , &
              Q21(kbdim,klev),  Q22(kbdim,klev)

  SELECT CASE (kscavBCtype)
     !CASE (1) ! prescribed scavenging ratio !!not implemented yet!!
         !sfrain(1:kproma,:,ktrac_phase,kmod) = 

     CASE (2) ! standard mode-wise scavenging ratio

         sfrain(1:kproma,:,ktrac_phase,kmod) = cbcr(kmod)*pfrain(1:kproma,:)

     CASE (3) ! aerosol size-dep scavenging coeff

         IF (sfrain(1,1,ktrac_phase,kmod) == UNDEF) THEN !calculation required!
           
            !>>SF #458 (replacing where statements)
            ll1(1:kproma,:) = (pfrain(1:kproma,:) > 0._dp)

            ztmp1(1:kproma,:) = MERGE(pfrain(1:kproma,:), 1._dp, ll1(1:kproma,:)) !SF 1._dp is a dummy value
            ztmp1(1:kproma,:) = FLOOR( 2._dp*log10( 3600._dp*ztmp1(1:kproma,:) ) + 5._dp )

            indexbcrx1(1:kproma,:) = MERGE( &
                                          MAX( 0, MIN( 9, INT( ztmp1(1:kproma,:) ))), &
                                          0, &
                                          ll1(1:kproma,:) )
            
            indexbcrx2(1:kproma,:) = MERGE( &
                                          MAX( 0, MIN( 9, INT( 1_dp + ztmp1(1:kproma,:) ))), &
                                          0, &
                                          ll1(1:kproma,:) )
            !<<SF #458 (replacing where statements)
   
            ! --- Assign below-cloud scavenging coefficients for rain:
            SELECT CASE (ktrac_phase)
   
               CASE (1) ! number
                  DO jk=ktop,klev
                     DO jl=1,kproma
                        Q11(jl,jk) = cscavbcrn(indexbcrx1(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                        Q12(jl,jk) = cscavbcrn(indexbcrx2(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                        Q21(jl,jk) = cscavbcrn(indexbcrx1(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
                        Q22(jl,jk) = cscavbcrn(indexbcrx2(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
                        X1(jl,jk)  = crainrate(indexbcrx1(jl,jk))
                        X2(jl,jk)  = crainrate(indexbcrx2(jl,jk))
                        Y1(jl,jk)  = caerorad(indexy1(jl,jk,ktrac_phase,kmod))
                        Y2(jl,jk)  = caerorad(indexy2(jl,jk,ktrac_phase,kmod))
                     ENDDO
                  ENDDO
   
               CASE (2) ! mass
                  DO jk=ktop,klev
                     DO jl=1,kproma
                        Q11(jl,jk) = cscavbcrm(indexbcrx1(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                        Q12(jl,jk) = cscavbcrm(indexbcrx2(jl,jk),indexy1(jl,jk,ktrac_phase,kmod))
                        Q21(jl,jk) = cscavbcrm(indexbcrx1(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
                        Q22(jl,jk) = cscavbcrm(indexbcrx2(jl,jk),indexy2(jl,jk,ktrac_phase,kmod))
                        X1(jl,jk)  = crainrate(indexbcrx1(jl,jk))
                        X2(jl,jk)  = crainrate(indexbcrx2(jl,jk))
                        Y1(jl,jk)  = caerorad(indexy1(jl,jk,ktrac_phase,kmod))
                        Y2(jl,jk)  = caerorad(indexy2(jl,jk,ktrac_phase,kmod))
                     ENDDO
                  ENDDO
   
            END SELECT !ktrac_phase
   
            CALL scavcoef_bilinterp (kproma, kbdim,  klev,  krow, ktop, &
                                     pfrain, mr(:,:,ktrac_phase,kmod),  &
                                     X1, X2, Y1, Y2,                    &
                                     Q11, Q12, Q21, Q22,                &
                                     sfrain(:,:,ktrac_phase,kmod))
         ENDIF !calculation required for sfrain

     CASE default

       CALL finish('bc_rain','wrong kscavBCtype value')

  END SELECT !kscavBCtype

  END SUBROUTINE bc_rain

  !!----------------------------------------------------------------------------
  SUBROUTINE bc_snow(kproma, kbdim, klev, krow, ktop, kmod, & !in
                     ktrac_phase, kscavBCtype,              & !in
                     pfsnow)                                  !in

  ! Below-cloud snow scavenging routine

!SF: output of this subroutine: sfsnow
!    sfsnow is defined in the whole module instead of being passed as an intent(out) because
!    it needs to be kept from one instance to the next

  USE mo_ham_tools,       ONLY: scavcoef_bilinterp
  USE mo_ham_wetdep_data, ONLY: cbcs, csnowcolleff, caerorad

  INTEGER, INTENT(in)  :: kproma, kbdim, klev, krow, ktop, kmod, ktrac_phase, kscavBCtype

  REAL(dp), INTENT(in) :: pfsnow(kbdim,klev)

  !local vars:
  INTEGER  :: jk, jl

  LOGICAL :: ll1(kbdim,klev)

  REAL(dp) :: ztmp1(kbdim,klev),                 &
              X1(kbdim,klev),   X2(kbdim,klev),  &
              Y1(kbdim,klev),   Y2(kbdim,klev),  &
              Q11(kbdim,klev),  Q12(kbdim,klev), &
              Q21(kbdim,klev),  Q22(kbdim,klev)

  SELECT CASE (kscavBCtype)
     !CASE (1) !SF prescribed scavenging ratio !!not implemented yet!!

         !sfsnow(1:kproma,:,ktrac_phase,kmod) = 

     CASE (2) !SF standard mode-wise scavenging ratio

         sfsnow(1:kproma,:,ktrac_phase,kmod) = cbcs(kmod)*pfsnow(1:kproma,:)

     CASE (3) !SF aerosol size-dep scavenging coeff

         IF (sfsnow(1,1,ktrac_phase,kmod) == UNDEF) THEN !calculation required for sfsnow

            DO jk=ktop,klev
               DO jl=1,kproma
                  X1(jl,jk)  = 1._dp  ! overwrite with dummy value to cause no interpolation in x
                  X2(jl,jk)  = 1._dp  ! overwrite with dummy value to cause no interpolation in x
                  Y1(jl,jk)  = caerorad(indexy1(jl,jk,ktrac_phase,kmod))
                  Y2(jl,jk)  = caerorad(indexy2(jl,jk,ktrac_phase,kmod))
                  Q11(jl,jk) = csnowcolleff(1,indexy1(jl,jk,ktrac_phase,kmod))
                  Q12(jl,jk) = csnowcolleff(1,indexy2(jl,jk,ktrac_phase,kmod))
                  Q21(jl,jk) = csnowcolleff(1,indexy1(jl,jk,ktrac_phase,kmod))
                  Q22(jl,jk) = csnowcolleff(1,indexy2(jl,jk,ktrac_phase,kmod))
               ENDDO
            ENDDO
   
            ztmp1(1:kproma,:) = 0._dp
            CALL scavcoef_bilinterp (kproma, kbdim,  klev,  krow, ktop, &
                                     ztmp1, mr(:,:,ktrac_phase,kmod), &
                                     X1, X2, Y1, Y2, &
                                     Q11, Q12, Q21, Q22, &
                                     sfsnow(:,:,ktrac_phase,kmod))
  
            !>>SF #458 (replacing where statements)
            ll1(1:kproma,:) = (pfsnow(1:kproma,:) > zeps) !SF#294 replaced 0. by zeps

            ztmp1(1:kproma,:) = (0.6_dp/0.027_dp)*pfsnow(1:kproma,:)*sfsnow(1:kproma,:,ktrac_phase,kmod)

            sfsnow(1:kproma,:,ktrac_phase,kmod) = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:)) 
            !<<SF #458 (replacing where statements)

         ENDIF !calculation required for sfsnow

     CASE default

         CALL finish('bc_snow','wrong kscavBCtype value')

  END SELECT !kscavBCtype

  END SUBROUTINE bc_snow

  !-------------------------------------------------------------------
  SUBROUTINE ham_conv_lfraq_so2(kproma, kbdim, klev, &
                                ptu,                 &
                                pxtu,                &
                                prhou,               &
                                pmlwc,               &
                                plfrac_so2)

  USE mo_tracdef,      ONLY: ntrac
  USE mo_ham,          ONLY: mw_so2, mw_so4, nham_subm, HAM_M7, HAM_SALSA
  USE mo_ham_m7_trac,  ONLY: idt_so2, idt_so4, idt_ms4ns, idt_ms4ks, idt_ms4as, idt_ms4cs
  USE mo_ham_salsa_trac,ONLY: idt_so2s=>idt_so2, idt_so4s=>idt_so4, idt_ms4 !TB (!SF #397)
  
  USE mo_species,      ONLY: speclist
  USE mo_ham_species,  ONLY: id_so2
  USE mo_ham_salsactl, ONLY: fn2b, in1a

  INTEGER, INTENT(in)  :: kproma, kbdim, klev
  REAL(dp), INTENT(in) :: ptu(kbdim,klev)
  REAL(dp), INTENT(in) :: pxtu(kbdim,klev,ntrac)
  REAL(dp), INTENT(in) :: prhou(kbdim,klev)
  REAL(dp), INTENT(in) :: pmlwc(kbdim,klev)

  REAL(dp), INTENT(out) :: plfrac_so2(kbdim,klev) !liquid fraction of SO2

  !--- local variables:
  LOGICAL :: ll1(kbdim, klev)

  REAL(dp) :: zqtp1(kbdim,klev), ze2(kbdim,klev),   &
              ze3(kbdim,klev), zfac(kbdim,klev),    &
              zfac4(kbdim,klev),                    &     !!mgs(S)!!
              zso4l(kbdim,klev), zso2l(kbdim,klev), &
              zza(kbdim,klev), zzb(kbdim,klev),     &
              zzp(kbdim,klev), zzq(kbdim,klev),     &
              zqhp(kbdim,klev),                     &
              zheneff(kbdim,klev), zhp(kbdim,klev), &
              zhenry_so2(2)

  INTEGER  :: jt, jn_idt_so4, jdx, idt_so2_loc
  INTEGER, ALLOCATABLE :: idt_so4_array(:)

  !--- Calculate the solubility of SO2. Total sulfate
  !    is only used to calculate the ph of cloud water:

  zhenry_so2(:) = speclist(id_so2)%henry(:) ! Henry's law constant and activation energy

  zqtp1(1:kproma,:) = 1._dp/ptu(1:kproma,:) - 1._dp/298._dp
  ze2(1:kproma,:)   = zhenry_so2(1)*EXP(zhenry_so2(2)*zqtp1(1:kproma,:))
  ze3(1:kproma,:)   = 1.2e-02_dp*EXP(2010._dp*zqtp1(1:kproma,:))

  !>>SF #458 (replacing where statements)
  ll1(1:kproma,:) = (pmlwc(1:kproma,:) > 1.E-15_dp)

  zfac(1:kproma,:)  = MERGE( &
                            1000._dp/(pmlwc(1:kproma,:)*mw_so2), &
                            0._dp, &
                            ll1(1:kproma,:) )

  zfac4(1:kproma,:)  = MERGE( &
                             1000._dp/(pmlwc(1:kproma,:)*mw_so4), &
                             0._dp, &
                             ll1(1:kproma,:) )
  !<<SF #458 (replacing where statements)

  !--- Set the list of relevant tracer id's (!SF #397: compliance for SALSA)
  IF (nham_subm == HAM_M7) THEN

     idt_so2_loc = idt_so2

     jn_idt_so4 = 5
     ALLOCATE(idt_so4_array(jn_idt_so4))
     idt_so4_array(1) = idt_so4
     idt_so4_array(2) = idt_ms4ns
     idt_so4_array(3) = idt_ms4ks
     idt_so4_array(4) = idt_ms4as
     idt_so4_array(5) = idt_ms4cs

  ELSEIF (nham_subm == HAM_SALSA) THEN

     idt_so2_loc = idt_so2s

     jn_idt_so4 = 1+fn2b
     ALLOCATE(idt_so4_array(jn_idt_so4))
     idt_so4_array(1) = idt_so4s
     DO jt=in1a,fn2b
        idt_so4_array(1+jt) = idt_ms4(jt)
     ENDDO
     
  ENDIF

  IF ( ((idt_so2_loc > 0) .AND. (.NOT. ANY( (idt_so4_array(:) <= 0) ))) .AND. &
       ((nham_subm == HAM_M7) .OR. (nham_subm == HAM_SALSA)) ) THEN
     
      zso4l(1:kproma,:) = 0._dp

      DO jt=1,jn_idt_so4
         jdx = idt_so4_array(jt)
         zso4l(1:kproma,:) = zso4l(1:kproma,:) + pxtu(1:kproma,:,jdx)
      ENDDO

      zso4l(1:kproma,:) = zso4l(1:kproma,:) * zfac4(1:kproma,:)    !!mgs(S)!!

      zso4l(1:kproma,:) = MAX(zso4l(1:kproma,:),0._dp)

      zso2l(1:kproma,:) = pxtu(1:kproma,:,idt_so2_loc)*zfac(1:kproma,:)
      zso2l(1:kproma,:) = MAX(zso2l(1:kproma,:),0._dp)

      zza(1:kproma,:) = ze2(1:kproma,:)*8.2e-02_dp*ptu(1:kproma,:)    &
                      * pmlwc(1:kproma,:)*prhou(1:kproma,:)*1.e-03_dp

      zzb(1:kproma,:) = 2.5e-06_dp + zso4l(1:kproma,:)
      zzp(1:kproma,:) = ( zza(1:kproma,:)*ze3(1:kproma,:) - zzb(1:kproma,:) &
                        - zza(1:kproma,:)*zzb(1:kproma,:)                   &
                        ) / (1._dp+zza(1:kproma,:)) * 0.5_dp

      zzq(1:kproma,:) = - zza(1:kproma,:)*ze3(1:kproma,:)                             &
                      * (zzb(1:kproma,:)+zso2l(1:kproma,:)) / (1._dp+zza(1:kproma,:))

      zhp(1:kproma,:) = -zzp(1:kproma,:) + SQRT(zzp(1:kproma,:)**2 - zzq(1:kproma,:))

      zqhp(1:kproma,:) = 1._dp / zhp(1:kproma,:)

      zheneff(1:kproma,:) = 1._dp + ze3(1:kproma,:)*zqhp(1:kproma,:)

      plfrac_so2(1:kproma,:) = zza(1:kproma,:)*zheneff(1:kproma,:)

  ELSE

      plfrac_so2(1:kproma,:) = 0._dp

  END IF

  END SUBROUTINE ham_conv_lfraq_so2

  !-----------------------------------------------------------------------

  SUBROUTINE prep_ham_mode_init(kproma, kbdim, klev)

  ! Sylvaine Ferrachat Sept. 2013:
  !
  ! This routine initializes all the necessary variables that are mode-dependent
  ! and not tracer-dependent.
  ! This allows to reduce significantly the computing load by computing these variables
  ! only once per mode and phase (number or mass), when relevant.
  !
  ! 'Mode' is meant to be general here (not M7-dependent)

    INTEGER, INTENT(in) :: kproma, kbdim, klev

    CALL init_var(kproma, kbdim, klev, indexy1)
    CALL init_var(kproma, kbdim, klev, indexy2)
    CALL init_var(kproma, kbdim, klev, mr)
    CALL init_var(kproma, kbdim, klev, rcritrad)
    CALL init_var(kproma, kbdim, klev, sfnuc)
    CALL init_var(kproma, kbdim, klev, sfimp)
    CALL init_var(kproma, kbdim, klev, sfrain)
    CALL init_var(kproma, kbdim, klev, sfsnow)

  END SUBROUTINE prep_ham_mode_init

  !-----------------------------------------------------------------------

  SUBROUTINE init_var_r_4d(kproma, kbdim, klev, pvar)

  ! small utility to allocate, if relevant, and initialize

    INTEGER, INTENT(in)     :: kproma, kbdim, klev
    REAL(dp), POINTER :: pvar(:,:,:,:)

    IF (.NOT. ASSOCIATED(pvar)) ALLOCATE(pvar(kbdim,klev,2,nclass))
    pvar(1:kproma,:,:,:) = UNDEF 

  END SUBROUTINE init_var_r_4d

  !-----------------------------------------------------------------------

  SUBROUTINE init_var_r_5d(kproma, kbdim, klev, pvar)

  ! small utility to allocate, if relevant, and initialize

    INTEGER, INTENT(in)     :: kproma, kbdim, klev
    REAL(dp), POINTER :: pvar(:,:,:,:,:)

    IF (.NOT. ASSOCIATED(pvar)) ALLOCATE(pvar(kbdim,klev,2,2,nclass))
    pvar(1:kproma,:,:,:,:) = UNDEF 

  END SUBROUTINE init_var_r_5d

  !-----------------------------------------------------------------------

  SUBROUTINE init_var_i_4d(kproma, kbdim, klev, kvar)

  ! small utility to allocate, if relevant, and initialize

    INTEGER, INTENT(in)     :: kproma, kbdim, klev
    INTEGER, POINTER :: kvar(:,:,:,:)

    IF (.NOT. ASSOCIATED(kvar)) ALLOCATE(kvar(kbdim,klev,2,nclass))
    kvar(1:kproma,:,:,:) = UNDEF 

  END SUBROUTINE init_var_i_4d

  !-----------------------------------------------------------------------
  SUBROUTINE ham_setscav(kt,                       &
                         kscavICtype, kscavBCtype, &
                         kscavICphase, kscavBCphase)

  ! this routine sets the wetdep flags to handle:
  !   - in-cloud and/or below-cloud scav
  !   - water and/or ice scav (resp. rain and/or snow)

  ! kscavICtype  = 0 no in-cloud scavenging
  !                1 prescribed in-cloud scavenging params
  !                2 standard in-cloud scav params (mode-wise step function)
  !                3 aerosol-size dependent in-cloud params
  !
  ! kscavBCtype  = 0 no below-cloud scavenging
  !                1 prescribed below-cloud scavenging params
  !                2 standard below-cloud scav params (mode-wise step function)
  !                3 aerosol-size dependent in-cloud params
  !
  ! kscavICphase = 0 no in-cloud scavenging
  !                1 in-cloud water-only scavenging
  !                2 in-cloud ice-only scavenging
  !                3 in-cloud water+ice scavenging
  !
  ! kscavBCphase = 0 no below-cloud scavenging
  !                1 below-cloud rain-only scavenging
  !                2 below-cloud snow-only scavenging
  !                3 below-cloud rain+snow scavenging          

  USE mo_tracdef, ONLY: trlist
  !SF demo - deactivated USE mo_species, ONLY: speclist

  INTEGER, INTENT(in)  :: kt
  INTEGER, INTENT(out) :: kscavICtype, kscavBCtype, &
                          kscavICphase, kscavBCphase

  ! local variables
  !SF demo - deactivated INTEGER :: ispec

  !--- Normal setting: 

  SELECT CASE (trlist%ti(kt)%nwetdep)
     CASE(0) !wet dep off
        kscavICtype  = 0
        kscavBCtype  = 0
        kscavICphase = 0
        kscavBCphase = 0
     CASE(1) !standard scavenging (mode-wise step function params)
        kscavICtype  = 2
        kscavBCtype  = 2
        kscavICphase = 3
        kscavBCphase = 3
     CASE(2) !standard in-cloud scavenging + aerosol size-dependent below-cloud scavenging
        kscavICtype  = 2
        kscavBCtype  = 3
        kscavICphase = 3
        kscavBCphase = 3
     CASE(3) !aerosol size-dependent in-cloud + below-cloud scavenging
        kscavICtype  = 3
        kscavBCtype  = 3
        kscavICphase = 3
        kscavBCphase = 3
!       CASE(4) !Define here your own additional setup
!          kscavICtype  = 
!          kscavBCtype  = 
!          kscavICphase = 
!          kscavBCphase = 
     CASE default
        CALL finish('ham_setscav','wrong nwetdep setting')
  END SELECT

  !--- Tracer- or species-specific setting (demo -activate it if necessary-)

!  ispec = trlist%ti(kt)%spid
!  IF (speclist(ispec)%shortname == 'BC') THEN
!        kscavICtype  = 
!        kscavBCtype  = 
!        kscavICphase = 
!        kscavBCphase = 
!  ENDIF

  END SUBROUTINE ham_setscav

END MODULE mo_ham_wetdep
