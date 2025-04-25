!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
module mo_psrad_interface
  USE mo_kind,            ONLY: wp
  USE mo_physical_constants, ONLY: grav, rd, avo,                               &
       &                        amco2, amch4, amn2o, amo3, amo2, amd, amw
  USE mo_exception,       ONLY: finish
  USE mo_radiation_parameters, ONLY: rad_perm, psctm, ssi_factor
  USE mo_rrtm_params,     ONLY: maxxsec, maxinpx, nbndsw, nbndlw
  USE mo_cloud_optics,    ONLY: cloud_optics
  USE mo_aero_kinne,      ONLY: set_aop_kinne  
  USE mo_aero_volc,       ONLY: add_aop_volc
  USE mo_aero_volc_tab,   ONLY: add_aop_volc_ham, add_aop_volc_crow
  USE mo_lrtm_setup,      ONLY: setup_lrtm
  USE mo_srtm_setup,      ONLY: setup_srtm
  USE mo_lrtm_driver,     ONLY: lrtm
  USE mo_srtm_driver,     ONLY: srtm
  USE mo_submodel,        ONLY: lanysubmodel
  USE mo_submodel_interface, &
       ONLY: radiation_subm_1, radiation_subm_2  
  USE mo_cosp_simulator,  ONLY: cosp_reffl, cosp_reffi,      &
                                locosp, cosp_f3d, Lisccp_sim,&
                                cisccp_cldtau3d, cisccp_cldemi3d
  USE mo_spec_sampling,   ONLY: spec_sampling_strategy, get_num_gpoints
  USE mo_random_numbers,  ONLY: seed_size_random

  IMPLICIT NONE
  PUBLIC :: setup_psrad, psrad_interface, & 
            lw_strat, sw_strat
  
  TYPE(spec_sampling_strategy), SAVE &
    &                   :: lw_strat, sw_strat !< Spectral sampling strategies for longwave, shortwave
  INTEGER, PARAMETER    :: rng_seed_size = 4
CONTAINS
  !---------------------------------------------------------------------------
  !>
  !! @brief  Sets up (initializes) radation routines
  !! 
  !! @remarks
  !!   Modify preset variables of module MO_RADIATION which control the 
  !!   configuration of the radiation scheme.
  !
  SUBROUTINE setup_psrad
    CALL setup_lrtm
    CALL setup_srtm
    IF(seed_size_random() > rng_seed_size) THEN 
      CALL finish('setup_psrad','Random number seed size (rng_seed_size) is smaller than true size')
    END IF
  END SUBROUTINE setup_psrad

  !-----------------------------------------------------------------------------
  !>
  !! @brief arranges input and calls rrtm sw and lw routines
  !! 
  !! @par Revision History
  !! Original Source Rewritten and renamed by B. Stevens (2009-08)
  !!
  !! @remarks
  !!   Because the RRTM indexes vertical levels differently than ECHAM a chief
  !!   function of thise routine is to reorder the input in the vertical.  In 
  !!   addition some cloud physical properties are prescribed, which are 
  !!   required to derive cloud optical properties
  !!
  !! @par The gases are passed into RRTM via two multi-constituent arrays: 
  !!   zwkl and wx_r. zwkl has maxinpx species and  wx_r has maxxsec species  
  !!   The species are identifed as follows.  
  !!     ZWKL [#/cm2]          WX_R [#/cm2]
  !!    index = 1 => H20     index = 1 => n/a
  !!    index = 2 => CO2     index = 2 => CFC11
  !!    index = 3 =>  O3     index = 3 => CFC12
  !!    index = 4 => N2O     index = 4 => n/a
  !!    index = 5 => n/a
  !!    index = 6 => CH4
  !!    index = 7 => O2
  !

  SUBROUTINE psrad_interface( &
       & iaero           ,kproma          ,kbdim           ,klev            ,&
       & krow            ,ktrac           ,ktype           ,nb_sw           ,&
       & laland          ,laglac          ,cemiss          ,pmu0            ,&
       & pgeom1          ,alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,&
       & alb_nir_dif     ,pp_fl           ,pp_hl           ,pp_sfc          ,&
       & tk_fl           ,tk_hl           ,tk_sfc          ,xm_vap          ,&
!>>SF
       & xm_liq          ,xm_ice          ,cdnc            , picnc, cld_frc         ,&
!       & xm_liq          ,xm_ice          ,cdnc            ,cld_frc         ,&
!<<SF
       & cld_cvr         ,xm_o3           ,xm_co2          ,xm_ch4          ,&
       & xm_n2o          ,xm_cfc          ,xm_o2           ,pxtm1           ,&
       & flx_uplw        ,flx_uplw_clr    ,flx_dnlw        ,flx_dnlw_clr    ,&
       & flx_upsw        ,flx_upsw_clr    ,flx_dnsw        ,flx_dnsw_clr    ,&
       & vis_frc_sfc     ,par_dn_sfc      ,nir_dff_frc     ,vis_dff_frc     ,&
       & par_dff_frc                                                         )

    INTEGER,INTENT(IN)  ::             &
         iaero,                        & !< aerosol control
         kproma,                       & !< number of longitudes
         kbdim,                        & !< first dimension of 2-d arrays
         krow,                         & !< first dimension of 2-d arrays
         klev,                         & !< number of levels
         ktrac,                        & !< number of tracers
         ktype(kbdim),                 & !< type of convection
         nb_sw                           !< number of shortwave bands

    LOGICAL,INTENT(IN) ::              &
         laland(kbdim),                & !< land sea mask, land=.true.
         laglac(kbdim)                   !< glacier mask, glacier=.true.

    REAL(WP),INTENT(IN)  ::            &
         cemiss,                       & !< surface emissivity
         pmu0(kbdim),                  & !< mu0 for solar zenith angle
         pgeom1(kbdim,klev),           & !< geopotential above ground
         alb_vis_dir(kbdim),           & !< surface albedo for vis range and dir light
         alb_nir_dir(kbdim),           & !< surface albedo for NIR range and dir light
         alb_vis_dif(kbdim),           & !< surface albedo for vis range and dif light
         alb_nir_dif(kbdim),           & !< surface albedo for NIR range and dif light
         pp_fl(kbdim,klev),            & !< full level pressure in Pa
         pp_hl(kbdim,klev+1),          & !< half level pressure in Pa
         pp_sfc(kbdim),                & !< surface pressure in Pa
         tk_fl(kbdim,klev),            & !< full level temperature in K
         tk_hl(kbdim,klev+1),          & !< half level temperature in K
         tk_sfc(kbdim),                & !< surface temperature in K
         xm_vap(kbdim,klev),           & !< specific humidity in g/g
         xm_liq(kbdim,klev),           & !< specific liquid water content
         xm_ice(kbdim,klev),           & !< specific ice content in g/g
         cdnc(kbdim,klev),             & !< cloud nuclei concentration
         picnc(kbdim,klev),            & !< ice crystal number concentration
         cld_frc(kbdim,klev),          & !< fractional cloud cover
         cld_cvr(kbdim),               & !< total cloud cover in m2/m2
         xm_o3(kbdim,klev),            & !< o3  mass mixing ratio
         xm_co2(kbdim,klev),           & !< co2 mass mixing ratio
         xm_ch4(kbdim,klev),           & !< ch4 mass mixing ratio
         xm_n2o(kbdim,klev),           & !< n2o mass mixing ratio
         xm_cfc(kbdim,klev,2),         & !< cfc volume mixing ratio
         xm_o2(kbdim,klev),            & !< o2  mass mixing ratio
         pxtm1(kbdim,klev,ktrac)         !< tracer mass mixing ratios

    REAL (wp), INTENT (OUT) ::         &
         flx_uplw    (kbdim,klev+1),   & !<   upward LW flux profile, all sky
         flx_uplw_clr(kbdim,klev+1),   & !<   upward LW flux profile, clear sky
         flx_dnlw    (kbdim,klev+1),   & !< downward LW flux profile, all sky
         flx_dnlw_clr(kbdim,klev+1),   & !< downward LW flux profile, clear sky
         flx_upsw    (kbdim,klev+1),   & !<   upward SW flux profile, all sky
         flx_upsw_clr(kbdim,klev+1),   & !<   upward SW flux profile, clear sky
         flx_dnsw    (kbdim,klev+1),   & !< downward SW flux profile, all sky
         flx_dnsw_clr(kbdim,klev+1),   & !< downward SW flux profile, clear sky
         vis_frc_sfc(kbdim)        ,   & !< Visible (250-680) fraction of net surface radiation
         par_dn_sfc (kbdim)        ,   & !< Downward Photosynthetically Active Radiation (PAR) at surface
         nir_dff_frc(kbdim)        ,   & !< Diffuse fraction of downward surface near-infrared radiation
         vis_dff_frc(kbdim)        ,   & !< Diffuse fraction of downward surface visible radiation 
         par_dff_frc(kbdim)              !< Diffuse fraction of downward surface PAR

    ! -------------------------------------------------------------------------------------
    INTEGER  :: jk, jl, jp, jkb,          & !< loop indicies
         icldlyr(kbdim,klev)                !< index for clear or cloudy

    REAL(wp) ::                           &
         zsemiss     (kbdim,nbndlw),      & !< LW surface emissivity by band
         ppd_hl      (kbdim,klev),        & !< pressure thickness in Pa
         pm_sfc      (kbdim),             & !< surface pressure in mb
         delta,                           & !< pressure thickness 
         zscratch                           !< scratch array
    !
    ! --- vertically reversed _vr variables
    !
    REAL(wp) ::                             &
         col_dry_vr(kbdim,klev),           & !< number of molecules/cm2 of
         pm_fl_vr  (kbdim,klev),           & !< full level pressure [mb] 
         pm_hl_vr  (kbdim,klev+1),         & !< half level pressure [mb] 
         tk_fl_vr  (kbdim,klev),           & !< full level temperature [K]
         tk_hl_vr  (kbdim,klev+1),         & !< half level temperature [K]
         cdnc_vr   (kbdim,klev),           & !< cloud nuclei concentration
!>>SF
         icnc_vr   (kbdim,klev),           & !< ice crystal number concentration
!<<SF
         cld_frc_vr(kbdim,klev),           & !< secure cloud fraction
         ziwgkg_vr (kbdim,klev),           & !< specific ice water content
         ziwc_vr   (kbdim,klev),           & !< ice water content per volume
         ziwp_vr   (kbdim,klev),           & !< ice water path in g/m2  
         zlwgkg_vr (kbdim,klev),           & !< specific liquid water content
         zlwp_vr   (kbdim,klev),           & !< liquid water path in g/m2  
         zlwc_vr   (kbdim,klev),           & !< liquid water content per
         re_drop   (kbdim,klev),           & !< effective radius of liquid
         re_cryst  (kbdim,klev),           & !< effective radius of ice
         wkl_vr       (kbdim,maxinpx,klev),& !< number of molecules/cm2 of
         wx_vr        (kbdim,maxxsec,klev),& !< number of molecules/cm2 of
         cld_tau_lw_vr(kbdim,klev,nbndlw), & !< LW optical thickness of clouds
         cld_tau_sw_vr(kbdim,klev,nb_sw),  & !< extincion
         cld_cg_sw_vr (kbdim,klev,nb_sw),  & !< asymmetry factor
         cld_piz_sw_vr(kbdim,klev,nb_sw),  & !< single scattering albedo
         aer_tau_lw_vr(kbdim,klev,nbndlw), & !< LW optical thickness of aerosols
         aer_tau_sw_vr(kbdim,klev,nb_sw),  & !< aerosol optical thickness
         aer_cg_sw_vr (kbdim,klev,nb_sw),  & !< aerosol asymmetry factor
         aer_piz_sw_vr(kbdim,klev,nb_sw)     !< aerosol single scattering albedo

    REAL(wp) ::                            &
         flx_uplw_vr(kbdim,klev+1),        & !< upward flux, total sky
         flx_uplw_clr_vr(kbdim,klev+1),    & !< upward flux, clear sky
         flx_dnlw_vr(kbdim,klev+1),        & !< downward flux, total sky
         flx_dnlw_clr_vr(kbdim,klev+1)       !< downward flux, clear sky

    !
    ! Random seeds for sampling. Needs to get somewhere upstream 
    !
    INTEGER :: rnseeds(kbdim,rng_seed_size)

    !
    ! Number of g-points per time step. Determine here to allow automatic array allocation in 
    !   lrtm, srtm subroutines. 
    !
    INTEGER :: n_gpts_ts

    ! 1.0 Constituent properties 
    !--------------------------------
!IBM* ASSERT(NODEPS)
    DO jk = 1, klev
      jkb = klev+1-jk
      DO jl = 1, kproma
        cld_frc_vr(jl,jk) = MAX(EPSILON(1.0_wp),cld_frc(jl,jkb))
        ziwgkg_vr(jl,jk)  = xm_ice(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
        zlwgkg_vr(jl,jk)  = xm_liq(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
      END DO
    END DO
    !
    ! --- control for zero, infintesimal or negative cloud fractions
    !
    WHERE (cld_frc_vr(1:kproma,:) > 2.0_wp*EPSILON(1.0_wp))
      icldlyr(1:kproma,:) = 1
    ELSEWHERE
      icldlyr(1:kproma,:)  = 0
      ziwgkg_vr(1:kproma,:) = 0.0_wp
      zlwgkg_vr(1:kproma,:) = 0.0_wp
    END WHERE
    !
    ! --- main constituent reordering
    !
!IBM* ASSERT(NODEPS)
    DO jl = 1, kproma
      pm_hl_vr(jl,klev+1) = 0.01_wp*pp_hl(jl,1)
      tk_hl_vr(jl,klev+1) = tk_hl(jl,1)
      pm_sfc(jl)          = 0.01_wp*pp_sfc(jl)
    END DO

!IBM* ASSERT(NODEPS)
    DO jk = 1, klev
      jkb = klev+1-jk
      DO jl = 1, kproma
        !
        ! --- thermodynamic arrays
        !
        pm_hl_vr(jl,jk) = 0.01_wp*pp_hl(jl,jkb+1)
        pm_fl_vr(jl,jk) = 0.01_wp*pp_fl(jl,jkb)
        tk_hl_vr(jl,jk) = tk_hl(jl,jkb+1)
        tk_fl_vr(jl,jk) = tk_fl(jl,jkb)
        !
        ! --- cloud properties
        !
        zscratch      = pp_fl(jl,jkb)/(tk_fl(jl,jkb) * rd)
        delta         = pp_hl(jl,jkb+1)-pp_hl(jl,jkb)
        ziwc_vr(jl,jk) = ziwgkg_vr(jl,jk)*zscratch
        ziwp_vr(jl,jk) = ziwgkg_vr(jl,jk)*delta/grav
        zlwc_vr(jl,jk) = zlwgkg_vr(jl,jk)*zscratch
        zlwp_vr(jl,jk) = zlwgkg_vr(jl,jk)*delta/grav
        cdnc_vr(jl,jk) = cdnc(jl,jkb)*1.e-6_wp
!>>SF
        icnc_vr(jl,jk) = picnc(jl,jkb)
!<<SF
        ! In ECHAM, pressure is related to mass of N2, O2, and Ar only. 
        ! Now col_dry is the number of particles (atoms, molecules) of dry air per cm^2
        col_dry_vr(jl,jk) = 0.1_wp*delta*avo/(grav*amd) 
      END DO
    END DO

!IBM* ASSERT(NODEPS)
    DO jk = 1, klev
      jkb = klev+1-jk
      DO jl = 1, kproma
        !
        ! --- radiatively active gases
        !
        wkl_vr(jl,:,jk)   = 0.0_wp
        ! Water vapor arrives as specific humidity (mass of water/mass of moist air) 
        !   Note that this equation fails if the atmosphere is water only
        wkl_vr(jl,1,jk)   = xm_vap(jl,jkb)/(1._wp - xm_vap(jl,jkb)) * amd/amw 
        ! Other quantities arrive as mass mixing ratios (mass of substance/mass of dry air) 
        wkl_vr(jl,2,jk)   = xm_co2(jl,jkb)*amd/amco2
        wkl_vr(jl,3,jk)   = xm_o3 (jl,jkb)*amd/amo3
        wkl_vr(jl,4,jk)   = xm_n2o(jl,jkb)*amd/amn2o
        wkl_vr(jl,6,jk)   = xm_ch4(jl,jkb)*amd/amch4
        wkl_vr(jl,7,jk)   = xm_o2 (jl,jkb)*amd/amo2
      END DO
    END DO
    !
    ! --- CFCs are in volume mixing ratio 
    ! 
    wx_vr(:,:,:) = 0._wp
!IBM* ASSERT(NODEPS)
    DO jk = 1, klev
      jkb = klev+1-jk
      DO jl = 1, kproma
        wx_vr(jl,2,jk) = xm_cfc(jl,jkb,1)
        wx_vr(jl,3,jk) = xm_cfc(jl,jkb,2)
      END DO
    END DO

    !
    ! -- Convert to molecules/cm^2
    !
    wkl_vr(1:kproma,1:7,1:klev) = wkl_vr(1:kproma,1:7,1:klev) * &
         SPREAD(col_dry_vr(1:kproma,1:klev), NCOPIES=7, DIM=2) 
    wx_vr (1:kproma,2:3,1:klev) = wx_vr (1:kproma,2:3,1:klev) * &
         SPREAD(col_dry_vr(1:kproma,1:klev), NCOPIES=2, DIM=2)
    !
    ! 2.0 Surface Properties
    ! --------------------------------
    zsemiss(1:kproma,:) = cemiss 
    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------
    ppd_hl(1:kproma,:) = pp_hl(1:kproma,2:klev+1)-pp_hl(1:kproma,1:klev)

!>>SF dod bugfix
    aer_tau_lw_vr(:,:,:) = 0.0_wp
    aer_tau_sw_vr(:,:,:) = 0.0_wp
    aer_piz_sw_vr(:,:,:) = 0.0_wp
    aer_cg_sw_vr(:,:,:)  = 0.0_wp
!<<SF dod bugfix

    SELECT CASE (iaero)
    CASE (0)
!>>SF dod bugfix
      !aer_tau_lw_vr(:,:,:) = 0.0_wp
      !aer_tau_sw_vr(:,:,:) = 0.0_wp
      !aer_piz_sw_vr(:,:,:) = 0.0_wp
      !aer_cg_sw_vr(:,:,:)  = 0.0_wp
!<<SF dod bugfix
    CASE (3)
      CALL set_aop_kinne( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndlw                ,nb_sw            ,&
           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
           & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
           & tk_fl            ,pgeom1 )
    CASE (5)
      CALL set_aop_kinne( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndlw                ,nb_sw            ,&
           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
           & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
           & tk_fl            ,pgeom1 )
      CALL add_aop_volc( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndlw                ,nb_sw            ,&
           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
           & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
           & tk_fl )
    CASE (6)
      CALL set_aop_kinne( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndlw                ,nb_sw            ,&
           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
           & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
           & tk_fl            ,pgeom1 )
      CALL add_aop_volc( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndlw                ,nb_sw            ,&
           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
           & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
           & tk_fl )
      CALL add_aop_volc_ham( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndlw                ,nb_sw            ,&
           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
           & aer_cg_sw_vr                                               )
    CASE (7)
      CALL set_aop_kinne( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndlw                ,nb_sw            ,&
           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
           & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
           & tk_fl            ,pgeom1 )
      CALL add_aop_volc_crow( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,nbndlw                ,nb_sw            ,&
           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
           & aer_cg_sw_vr                                               )
    CASE DEFAULT
    END SELECT

    CALL cloud_optics(                                                  &
         & laglac        ,laland        ,kproma        ,kbdim          ,& 
         & klev          , ktype        ,nbndlw        ,nbndsw         ,&
         & icldlyr       ,zlwp_vr       ,ziwp_vr       ,zlwc_vr        ,&
!>>SF
         & ziwc_vr       ,cdnc_vr       , icnc_vr, cld_tau_lw_vr ,cld_tau_sw_vr  ,&
!         & ziwc_vr       ,cdnc_vr       ,cld_tau_lw_vr ,cld_tau_sw_vr  ,&
!<<SF
         & cld_piz_sw_vr ,cld_cg_sw_vr  ,re_drop       ,re_cryst    )  

    IF (locosp) THEN
      DO jk=1,klev
        jkb = klev+1-jk
        DO jl=1,kproma
          !
          ! To fix - particle sizes are local variables that aren't set! 
          !
          cosp_reffl(jl,jk,krow) = re_drop(jl,jkb)
          cosp_reffi(jl,jk,krow) = re_cryst(jl,jkb)
          cosp_f3d(jl,jk,krow)   = cld_frc(jl,jk) ! Or cld_frc_vr?
        END DO
      END DO
      IF ( Lisccp_sim ) THEN   
        !don't really need vert. re-arrange ...
        DO jk=1,klev
          jkb = klev+1-jk
          DO jl=1,kproma
            cisccp_cldtau3d(jl,jk,krow) = &
                 cld_tau_sw_vr(jl,jkb,9)  ! band 9 is 625 - 778 nm, needed is 670 nm
          END DO
        END DO
        DO jk=1,klev
          jkb = klev+1-jk
          DO jl=1,kproma
            !this is ONLY o.k. as long as wp equals dp, else conversion needed
            cisccp_cldemi3d(jl,jk,krow) = &
                 1._wp - exp(-1._wp*cld_tau_lw_vr(jl,jkb,6)) ! band 6 is 820 - 980 cm-1, 
            !  we need 10.5 µm channel
          END DO
        END DO
      END IF
    END IF

    !
    ! 3.5 Interface for submodels that provide aerosol and/or cloud radiative properties:
    ! -----------------------------------------------------------------------------------
    IF (lanysubmodel)  CALL radiation_subm_1(                     &
         kproma           ,kbdim            ,klev         ,krow  ,&
         ktrac            ,iaero            ,nbndlw       ,nb_sw ,&
         aer_tau_sw_vr    ,aer_piz_sw_vr    ,aer_cg_sw_vr        ,&
         aer_tau_lw_vr    , ppd_hl           ,pxtm1               )
    !
    ! 4.0 Radiative Transfer Routines
    ! --------------------------------
    !
    ! Seeds for random numbers come from least significant digits of pressure field 
    !
    rnseeds(1:kproma,1:rng_seed_size) = (pm_fl_vr(1:kproma,1:rng_seed_size) -  &
         int(pm_fl_vr(1:kproma,1:rng_seed_size)))* 1E9 + rad_perm

    n_gpts_ts = get_num_gpoints(lw_strat)
    CALL lrtm(kproma                                                          ,&
         & kbdim           ,klev            ,pm_fl_vr        ,pm_sfc          ,&
         & tk_fl_vr        ,tk_hl_vr        ,tk_sfc          ,wkl_vr          ,&
         & wx_vr           ,col_dry_vr      ,zsemiss         ,cld_frc_vr      ,&
         & cld_tau_lw_vr   ,aer_tau_lw_vr   ,rnseeds         ,lw_strat        ,&
         & n_gpts_ts       ,flx_uplw_vr     ,flx_dnlw_vr     ,flx_uplw_clr_vr ,&
         & flx_dnlw_clr_vr )
    !
    ! Reset random seeds so SW doesn't depend on what's happened in LW but is also independent
    !
    rnseeds(1:kproma,1:rng_seed_size) = (pm_fl_vr(1:kproma,rng_seed_size:1:-1) - &
         int(pm_fl_vr(1:kproma,rng_seed_size:1:-1)))* 1E9 + rad_perm
    n_gpts_ts = get_num_gpoints(sw_strat)
    ! 
    ! Potential pitfall - we're passing every argument but some may not be present
    !
    CALL srtm(kproma                                                           , & 
         &  kbdim           ,klev            ,pm_fl_vr        ,tk_fl_vr        , &
         &  wkl_vr          ,col_dry_vr      ,alb_vis_dir     ,alb_vis_dif     , &
         &  alb_nir_dir     ,alb_nir_dif     ,pmu0            ,ssi_factor      , &
         &  psctm           ,cld_frc_vr      ,cld_tau_sw_vr   ,cld_cg_sw_vr    , &
         &  cld_piz_sw_vr   ,aer_tau_sw_vr   ,aer_cg_sw_vr    ,aer_piz_sw_vr   , & 
         &  rnseeds         ,sw_strat        ,n_gpts_ts       ,flx_dnsw        , &
         &  flx_upsw        ,flx_dnsw_clr    ,flx_upsw_clr    ,vis_frc_sfc     , &
         &  par_dn_sfc      ,nir_dff_frc     ,vis_dff_frc     ,par_dff_frc       )
    !
    ! 5.0 Post Processing
    ! --------------------------------
    !
    ! Lw fluxes are vertically revered but SW fluxes are not
    !
    flx_uplw    (1:kproma,1:klev+1) = flx_uplw_vr    (1:kproma,klev+1:1:-1) 
    flx_uplw_clr(1:kproma,1:klev+1) = flx_uplw_clr_vr(1:kproma,klev+1:1:-1)
    flx_dnlw    (1:kproma,1:klev+1) = flx_dnlw_vr    (1:kproma,klev+1:1:-1)
    flx_dnlw_clr(1:kproma,1:klev+1) = flx_dnlw_clr_vr(1:kproma,klev+1:1:-1)

    !
    ! 6.0 Interface for submodel diagnosics after radiation calculation:
    ! ------------------------------------------------------------------
    IF (lanysubmodel)  CALL radiation_subm_2(kproma, kbdim, krow, klev, ktrac, iaero, pxtm1)

  END SUBROUTINE psrad_interface
end module mo_psrad_interface
