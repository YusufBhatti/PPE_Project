!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide interface to rrtmg shortwave radiation
!!
!! @remarks
!!   This module contains routines that provide the interface between ECHAM
!!   and the AER RRTMG radiation code.  Mostly it organizes and calculates the 
!!   information necessary to call the radiative transfer solvers.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2010-08)
!!         Robert Pincus, U. Colorado, visiting MPI-M, Hamburg (2011-07)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard) 
!!   code previously provided by AER and copyrighed by them.  The authors of the
!!   original AER code are: Eli J. Mlawer, Jennifer S. Delamere, Michael J. 
!!   Iacono and Shepard A. Clough with acknowledgments to Steven J. Taubman, 
!!   Karen Cady-Pereira, Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert 
!!   Bergstrom. The rewrites were designed to better interface with the structure
!!   of the ICON family of models and elements of the ICON programming standard.
!!
!! @par Copyright
!!   The AER copyright
!!   on the original code is as follows: Copyright 2002-2009, Atmospheric and
!!   Environmental Research, Inc. (AER). This software may be used, copied, or
!!   redistributed as long as it is not sold and this copyright notice is
!!   reproduced on each copy made.  This model is provided as is without any
!!   express or implied warranties. (http://www.rtweb.aer.com/)               
!! 
!
MODULE mo_srtm_driver

  USE mo_kind,               ONLY : wp
  USE mo_rrtm_params,        ONLY : nbndsw, ngptsw, jpband, jpb1, jpb2
  USE mo_radiation_parameters, &
                             ONLY : i_overlap, l_do_sep_clear_sky, rad_undef
  USE mo_srtm_setup,         ONLY : ngb, wavenum2, ssi_default, delwave
  USE mo_srtm_gas_optics,    ONLY : gpt_taumol, ih2o, ich4, ico2, io2, io3
  USE mo_rrtm_coeffs,        ONLY : srtm_coeffs
  USE mo_srtm_solver,        ONLY : srtm_solver, delta_scale, two_stream
  USE mo_cld_sampling,       ONLY : sample_cld_state
  USE mo_spec_sampling,      ONLY : spec_sampling_strategy, get_gpoint_set

  IMPLICIT NONE

  PRIVATE 
  PUBLIC :: srtm

  REAL (wp), PARAMETER     :: nir_vis_boundary   = 14500._wp, &
       & par_lower_boundary = 14285.7143_wp,  & !< equivalent to 700nm wavelength
       & par_upper_boundary = 25000._wp         !< equivalent to 400nm wavelength
  REAL (wp), PARAMETER      :: zepsec = 1.e-06_wp     !< epsilon
  INTEGER :: i
#if defined (__PGI)
  REAL (wp)                 :: frc_par_array(nbndsw)
  REAL (wp)                 ::       frc_vis(nbndsw)
#else
  REAL (wp), PARAMETER      :: frc_par_array(nbndsw) = &
             (/ (0.0_wp, i = 1, 8), 0.533725_wp, 1.0_wp, 0.550164_wp, (0.0_wp, i = 12, nbndsw) /) 
  REAL (wp), PARAMETER      ::       frc_vis(1:nbndsw) = MAX(0.0_wp, MIN(1.0_wp, & 
           (/ ((wavenum2(i+jpb1-1) - nir_vis_boundary) / delwave(i+jpb1-1), i = 1, nbndsw) /) ))
#endif

CONTAINS

  SUBROUTINE srtm( kproma                                                    , &
       &  kbdim           ,klev            ,play            ,tlay            , &
       &  wkl             ,coldry          ,asdir           ,asdif           , &
       &  aldir           ,aldif           ,prmu0           ,ssi_factor      , &
       &  psctm           ,cld_frc         ,cld_tau_sw      ,cld_cg_sw       , &
       &  cld_piz_sw      ,aer_tau_sw      ,aer_cg_sw       ,aer_piz_sw      , &
       &  rnseeds         ,strategy        ,n_gpts_ts       ,flxd_sw         , &
       &  flxu_sw         ,flxd_sw_clr     ,flxu_sw_clr     ,vis_frc_sfc     , &
       &  par_dn_sfc      ,nir_dff_frc     ,vis_dff_frc     ,par_dff_frc       )

    ! This program is the driver for RRTMG_SW, the AER SW radiation model for 
    !  application to GCMs, that has been adapted from RRTM_SW for improved
    !  efficiency and to provide fractional cloudiness and cloud overlap
    !  capability using McICA.
    !
    !
    INTEGER, INTENT(in) :: kproma        !< Number of horizontal columns 
    INTEGER, INTENT(in) :: kbdim         !< Maximum number of columns as declared in calling (sub)prog.
    INTEGER, INTENT(in) :: klev          !< Number of model layers

    REAL(wp), INTENT(in) :: &
         & psctm,            & !< local solar constant
         & play(kbdim,klev), & !< Layer pressures [hPa, mb]
         & tlay(kbdim,klev), & !< Layer temperatures [K]
         & prmu0(kbdim),     & !< Solar zenith angle
         & wkl(:,:,:),       & !< Gas volume mixing ratios [mol/frac]
         & coldry(kbdim,klev)  !< Column dry amount

    REAL(wp), INTENT(in) :: asdir(kbdim) !< UV/vis direct sfc albedo
    REAL(wp), INTENT(in) :: aldir(kbdim) !< Near-IR direct sfc albedo
    REAL(wp), INTENT(in) :: asdif(kbdim) !< UV/vis diffuse sfc albedo 
    REAL(wp), INTENT(in) :: aldif(kbdim) !< Near-IR diffuse sfc albedo

    REAL(wp), INTENT(in) :: &
         & ssi_factor(nbndsw)              , & !< solar constant factor, dimension(nbndsw)
         & cld_frc(kbdim,klev)     , & !< dimensions
         & cld_tau_sw(kbdim,klev,nbndsw), & !< dimensions (kproma, klev, nbndsw)
         & cld_cg_sw (kbdim,klev,nbndsw), &
         & cld_piz_sw(kbdim,klev,nbndsw), &
         & aer_tau_sw(kbdim,klev,nbndsw), &
         & aer_cg_sw (kbdim,klev,nbndsw), &
         & aer_piz_sw(kbdim,klev,nbndsw)

    INTEGER, INTENT(IN   ) :: n_gpts_ts      !< Number of gpts to sample
    INTEGER, INTENT(INOUT) :: rnseeds(:, :)  !< Seeds for random number generator (kbdim, :) 

    TYPE(spec_sampling_strategy), INTENT(IN) :: strategy

    REAL(WP),   INTENT(OUT)   ::   &
         flxd_sw    (kbdim,klev+1)        , & !< downward flux total sky
         flxd_sw_clr(kbdim,klev+1)        , & !< downward flux clear sky
         flxu_sw    (kbdim,klev+1)        , & !< upward flux total sky
         flxu_sw_clr(kbdim,klev+1)        , & !< upward flux clear sky
         vis_frc_sfc(kbdim)               , & !< Visible (250-680) fraction of net surface radiation
         par_dn_sfc(kbdim)                , & !< Downward Photosynthetically Active Radiation (PAR) at surface
         nir_dff_frc(kbdim)               , & !< Diffuse fraction of downward surface near-infrared radiation
         vis_dff_frc(kbdim)               , & !< Diffuse fraction of downward surface visible radiation 
         par_dff_frc(kbdim)                   !< Diffuse fraction of downward surface PAR

    ! ----------------
    LOGICAL     :: sunUp(kbdim)         !< Mask for sunlit points
    INTEGER(WP) :: idxSunUp(kbdim)      !< Indicies of sunlit points

    INTEGER :: jl   !< column loop index
    INTEGER :: jk   !< level index
    INTEGER :: ib   !< band loop index (starting at 1)
    INTEGER :: jb   !< band loop index (starting at first solar band)
    INTEGER :: ig   !< cumulative g-point index
    INTEGER :: igpt !< g-point to use in this calculation
    INTEGER :: ibs(kbdim, n_gpts_ts), igs(kbdim, n_gpts_ts) ! bands and g-points used for samples

    REAL(wp) :: cossza(kbdim)       !< Cosine of solar zenith angle
    REAL(wp) :: adjflux(jpband)      !< adjustment for current Earth/Sun distance
    REAL(wp) :: albdir(kbdim,n_gpts_ts)       !< surface albedo, direct          
    REAL(wp) :: albdif(kbdim,n_gpts_ts)       !< surface albedo, diffuse         

    ! -----------------
    ! Variables for gas optics calculations
    INTEGER :: laytrop(kbdim)          !< 100 hPa layer index
    INTEGER :: jp (kbdim,klev)         !< 
    INTEGER :: jt (kbdim,klev)         !<
    INTEGER :: jt1(kbdim,klev)         !<

    REAL(wp) :: colh2o(kbdim,klev)     !< column amount (h2o)
    REAL(wp) :: colco2(kbdim,klev)     !< column amount (co2)
    REAL(wp) :: colo3 (kbdim,klev)     !< column amount (o3)
    REAL(wp) :: coln2o(kbdim,klev)     !< column amount (n2o)
    REAL(wp) :: colch4(kbdim,klev)     !< column amount (ch4)
    REAL(wp) :: colo2 (kbdim,klev)     !< column amount (o2)
    REAL(wp) :: colmol(kbdim,klev)     !< column amount
    REAL(wp) :: gases (kbdim,klev,5)

    INTEGER  :: indself (kbdim,klev)
    INTEGER  :: indfor  (kbdim,klev)
    REAL(wp) :: selffac (kbdim,klev)
    REAL(wp) :: selffrac(kbdim,klev)
    REAL(wp) :: forfac  (kbdim,klev)
    REAL(wp) :: forfrac (kbdim,klev)

    REAL(wp) :: &                
         fac00(kbdim,klev), fac01(kbdim,klev), &
         fac10(kbdim,klev), fac11(kbdim,klev) 
    ! -----------------

    REAL(WP) :: bnd_wght(nbndsw)
    !
    ! Cloud, aerosol, and total optical properties per sample
    !
    REAL(wp) :: ztauc(kbdim,klev,n_gpts_ts)   !< cloud optical depth
    REAL(wp) :: zasyc(kbdim,klev,n_gpts_ts)   !< cloud asymmetry parameter 
    REAL(wp) :: zomgc(kbdim,klev,n_gpts_ts)   !< cloud single scattering albedo
    REAL(wp) :: ztaua(kbdim,klev,n_gpts_ts)   !< total aerosol optical depth
    REAL(wp) :: zasya(kbdim,klev,n_gpts_ts)   !< total aerosol asymmetry parameter 
    REAL(wp) :: zomga(kbdim,klev,n_gpts_ts)   !< total aerosol single scattering albedo

    REAL(wp) :: ztaut (kbdim,klev,n_gpts_ts)   !< total optical depth
    REAL(wp) :: zasyt (kbdim,klev,n_gpts_ts)   !< total asymmetry parameter 
    REAL(wp) :: zomgt (kbdim,klev,n_gpts_ts)   !< total single scattering albedo
    !
    ! When computing separate clear-sky fluxes, variables need to compute new optical 
    !   properties for initally cloudy cells
    !
    REAL(wp) :: Rdir(kbdim,klev), Rdif(kbdim,klev), & !< Direct and diffuse reflectance, transmittance
         Tdir(kbdim,klev), Tdif(kbdim,klev)    !< from two-stream calculation
    !
    ! Band-by-band fluxes, for computing spectrally-resolved surface fluxes fluxes
    !
    REAL(wp) :: zbbfu   (kbdim,nbndsw) !< All-sky   upward surface sw flux [w/m2]
    REAL(wp) :: zbbfd   (kbdim,nbndsw) !< All-sky downward surface sw flux [w/m2]
    REAL(wp) :: zbbfddir(kbdim,nbndsw) !< All-sky downward surface direct sw flux [w/m2]
    
    !
    ! Narrow band (g-point, or gp) fluxes
    !
    REAL(wp) :: zgpfu(kbdim,klev+1) !< Fullsky upward sw flux [w/m2]
    REAL(wp) :: zgpfd(kbdim,klev+1) !< Fullsky downward sw flux [w/m2]
    REAL(wp) :: ztdbt(kbdim)        !< Fullsky downward direct sw flux [w/m2]
    REAL(wp) :: zgpcu(kbdim,klev+1) !< Clearsky upward sw flux [w/m2]
    REAL(wp) :: zgpcd(kbdim,klev+1) !< Clearsky downward sw flux [w/m2]
    REAL(wp) :: ztdbtc(kbdim)       !< Clearsky downward direct sw flux [w/m2]
    !
    ! parameters from gas optics
    !
    REAL(wp) :: ztaug(kbdim,klev,n_gpts_ts), ztaur(kbdim,klev,n_gpts_ts)
    REAL(wp) :: zsflxzen(kbdim,n_gpts_ts)

    REAL(wp) :: zincflx(kbdim)
    REAL(wp) :: zfvis(kbdim,nbndsw), zfnir(kbdim,nbndsw), zfpar(kbdim,nbndsw)

    ! Variables for sampling strategy 
    REAL(WP) :: gpt_scaling, clrSky_scaling(1:kbdim)

    LOGICAL  ::    cldMask       (kbdim, klev, n_gpts_ts), & !< cloud mask in each cell
                colcldMask       (kbdim,       n_gpts_ts)



    !
    ! Initializations
    ! --------------------------------
#if defined (__PGI)
    frc_par_array(:) = (/ (0.0_wp, i = 1, 8), 0.533725_wp, 1.0_wp, 0.550164_wp, (0.0_wp, i = 12, nbndsw) /) 
    frc_vis(:) = MAX(0.0_wp, MIN(1.0_wp, (/ ((wavenum2(i+jpb1-1) - nir_vis_boundary) / delwave(i+jpb1-1), i = 1, nbndsw) /) ))
#endif
    !
    ! --- weight radiation within a band for the solar cycle ---
    ! psctm contains TSI (the "solar constant") scaled with the
    ! Sun-Earth distance. ssi_factor contains the relative contribution
    ! of each band to TSI. ssi_default is the (originally only
    ! implicitly defined) solar flux in the 14 bands.
    !
    bnd_wght(:) = psctm*ssi_factor(:) / ssi_default(:)
    DO jb = jpb1,jpb2
      adjflux(jb) = bnd_wght(jb-jpb1+1)
    ENDDO

    ! Which input points are sunlit? 
    sunUp(1:kproma)    = prmu0(1:kproma) > 0._wp 
    idxSunUp(1:kproma) = PACK( (/ (jl, jl = 1, kproma) /),  sunUp(1:kproma))

    !
    ! ---  1.0 Choose a set of g-points to do consistent with the spectral sampling strategy
    ! --------------------------------
    gpt_scaling = REAL(ngptsw)/REAL(n_gpts_ts)
    igs(1:kproma,1:n_gpts_ts) = get_gpoint_set(kproma, kbdim, strategy, rnseeds)

    ! Save the band nunber associated with each gpoint
    DO ig = 1, n_gpts_ts
      DO jl = 1, kproma  
        ibs(jl, ig) = ngb(igs(jl, ig)) - 15 ! arrays start from 1 but band numbers from 16 (shared with LW) 
      END DO
    END DO
    !
    ! ---  2.0 Optical properties 
    !
    ! ---  2.1 Cloud, aerosol, surface optical properties. 
    ! --------------------------------
    !
    ! Cloud optical depth is only saved for the band associated with this g-point
    !   We sample clouds first because we may want to adjust water vapor based 
    !   on presence/absence of clouds
    !
    CALL sample_cld_state(kproma, kbdim, klev, n_gpts_ts, rnseeds(:,:), i_overlap, &
         cld_frc, cldMask)
!IBM* ASSERT(NODEPS)
    DO ig = 1, n_gpts_ts
      DO jl = 1, kproma
        ztauc(jl,:,ig) = MERGE(cld_tau_sw(jl,1:klev,ibs(jl,ig)), 0._wp, cldMask(jl,:,ig)) 
        zasyc(jl,:,ig) = MERGE(cld_cg_sw (jl,1:klev,ibs(jl,ig)), 0._wp, cldMask(jl,:,ig)) 
        zomgc(jl,:,ig) = MERGE(cld_piz_sw(jl,1:klev,ibs(jl,ig)), 1._wp, cldMask(jl,:,ig))
      END DO
    END DO ! Loop over samples - done with cloud optical depth calculations 

!IBM* ASSERT(NODEPS)
    DO ig = 1, n_gpts_ts
      DO jl = 1, kproma  
        ib   = ibs(jl, ig) 
        ! Aerosol optical properties
        ztaua(jl,1:klev,ig) = aer_tau_sw(jl,1:klev,ib)
        zasya(jl,1:klev,ig) = aer_cg_sw (jl,1:klev,ib)
        zomga(jl,1:klev,ig) = aer_piz_sw(jl,1:klev,ib)
        albdif(jl,ig) = asdif(jl)*frc_vis(ib) + aldif(jl)*(1.0_wp - frc_vis(ib))
        albdir(jl,ig) = asdir(jl)*frc_vis(ib) + aldir(jl)*(1.0_wp - frc_vis(ib))
      END DO
    END DO

    !
    ! Cloud masks for sorting out clear skies - by cell and by column
    !   These calculations are only needed if l_do_sep_clear_sky is false

    IF(.not. l_do_sep_clear_sky) THEN 
      !
      ! Are any layers cloudy? 
      !
      colcldMask(1:kproma, :) = ANY(cldMask(1:kproma,:,:), DIM=2)
      !
      ! Clear-sky scaling is gpt_scaling/frac_clr or 0 if all samples are cloudy 
      !
      clrSky_scaling(1:kproma) = gpt_scaling *                  &
           MERGE( REAL(n_gpts_ts,KIND=wp) / &
           (REAL(n_gpts_ts - COUNT(colCldMask(1:kproma,:),DIM=2),KIND=wp)), &
           0._wp,                   &
           ANY(.not. colCldMask(1:kproma,:),DIM=2))
    END IF
    !
    ! ---  2.2 Gas optical depth calculations
    ! 
    ! --------------------------------
    !
    ! 2.2.1  Calculate information needed by the radiative transfer routine
    ! that is specific to this atmosphere, especially some of the 
    ! coefficients and indices needed to compute the optical depths
    ! by interpolating data from stored reference atmospheres. 
    ! The coefficients are functions of temperature and pressure and remain the same
    ! for all g-point samples.
    ! If gas concentrations, temperatures, or pressures vary with sample (ig) 
    !   the coefficients need to be calculated inside the loop over samples
    ! --------------------------------
    ! 
    CALL srtm_coeffs(kproma,kbdim, klev                               , &
         & play          ,tlay          ,coldry        ,wkl           , &
         & laytrop       ,jp            ,jt            ,jt1           , &
         & colch4        ,colco2        ,                               &
         & colh2o        ,colmol        ,coln2o        ,colo2         , &
         & colo3         ,fac00         ,fac01         ,fac10         , &
         & fac11         ,selffac       ,selffrac      ,indself       , &
         & forfac        ,forfrac       ,indfor)

    gases(1:kproma,:,ih2o) = colh2o(1:kproma,1:klev)
    gases(1:kproma,:,ico2) = colco2(1:kproma,1:klev)
    gases(1:kproma,:,ich4) = colch4(1:kproma,1:klev)
    gases(1:kproma,:,io2 ) = colo2 (1:kproma,1:klev)
    gases(1:kproma,:,io3 ) = colo3 (1:kproma,1:klev)

    !
    !  -- 2.2.2 Loop over g-points calculating gas optical properties. 
    !
    ! --------------------------------
!IBM* ASSERT(NODEPS)
    DO ig = 1, n_gpts_ts 
      DO jl = 1, kproma
        igpt = igs(jl, ig)
        CALL gpt_taumol(klev      ,igpt          ,                              &
             & jp      (jl,:) ,jt      (jl,:),jt1   (jl,:)  ,laytrop(jl)  , &
             & indself (jl,:) ,indfor  (jl,:),gases (jl,:,:),colmol (jl,:), &
             & fac00   (jl,:) ,fac01   (jl,:),fac10 (jl,:)  ,fac11  (jl,:), &
             & selffac (jl,:) ,selffrac(jl,:),forfac(jl,:)  ,forfrac(jl,:), &
             & zsflxzen(jl,ig),ztaug(jl, 1:klev, ig), ztaur(jl, 1:klev, ig))
      END DO
    END DO
    !
    ! --- 2.3 Combine optical properties of aerosols, clouds, gases (absorption + Rayleigh scattering)
    !
    ztaut(1:kproma,:,:) = ztaur(1:kproma,:,:) &
         & + ztaua(1:kproma,:,:) &
         & + ztaug(1:kproma,:,:) & 
         & + ztauc(1:kproma,:,:)
    zomgt(1:kproma,:,:) = ztaur(1:kproma,:,:) * 1.0_wp &
         & + ztaua(1:kproma,:,:) * zomga(1:kproma,:,:) & 
         & + ztauc(1:kproma,:,:) * zomgc(1:kproma,:,:) 
    zasyt(1:kproma,:,:) =(zasya(1:kproma,:,:) * zomga(1:kproma,:,:) *  ztaua(1:kproma,:,:) & 
         & + zasyc(1:kproma,:,:) * zomgc(1:kproma,:,:) *  ztauc(1:kproma,:,:)) / zomgt(1:kproma,:,:)
    zomgt(1:kproma,:,:) = zomgt(1:kproma,:,:) / ztaut(1:kproma,:,:)

    ! 
    ! ---  3.0 Compute radiative transfer.
    ! --------------------------------

    ! ---  3.1 Initialization
    !
    ! ---  3.1.1 Accumulated flux arrays
    !
    flxu_sw    (1:kproma,1:klev+1) = 0.0_wp
    flxd_sw    (1:kproma,1:klev+1) = 0.0_wp
    flxu_sw_clr(1:kproma,1:klev+1) = 0.0_wp
    flxd_sw_clr(1:kproma,1:klev+1) = 0.0_wp

    zbbfu   (:,:) = 0._wp
    zbbfd   (:,:) = 0._wp
    zbbfddir(:,:) = 0._wp

    !
    ! ---  3.1.2 Solar illumination
    !
    cossza(1:kproma) = prmu0(1:kproma) 

    !
    ! --- 3.2 Compute fluxes for each set of samples in turn
    ! 
    DO ig = 1, n_gpts_ts
      zincflx(1:kproma) =  adjflux(ibs(1:kproma,ig)+15) * zsflxzen(1:kproma,ig) * cossza(1:kproma)
      !
      ! All (or full) sky calculation 
      !
      IF(l_do_sep_clear_sky) THEN
        !
        ! Delta scale values of tau, g, w0 for solver. This scaling is done inside the 
        !   solver when l_do_sep_clear_sky, so these variables have different values and the end
        !   of this if block depending on l_do_sep_clear_sky (they aren't used again, though) 
        ! 
        CALL delta_scale(kproma, kbdim, klev, ztaut(:,:,ig), zasyt(:,:,ig), zomgt(:,:,ig)) 
        CALL two_stream(kproma, kbdim, klev                                    ,&
             & cossza(:)      ,ztaut(:,:,ig)  ,zomgt(:,:,ig)  ,zasyt(:,:,ig)   ,&
             & Rdir(:,:)      ,Rdif(:,:)      ,Tdir(:,:)      ,Tdif(:,:)) 

        CALL srtm_solver(kproma, kbdim, klev                                   ,&
             & albdif(:,  ig)  ,albdir(:, ig)  ,cossza(:)      ,ztaut(:,:,ig)  ,&
             & Rdir(:,:)       ,Rdif(:,:)      , Tdir(:,:)     ,Tdif(:,:)      ,&
             & zgpfd(:,:)      ,zgpfu(:,:)     , ztdbt(:))
      ELSE
        ! 
        ! Compute fluxes directly from optical properties
        !
        CALL srtm_solver(kproma, kbdim, klev, albdif(:,  ig), albdir(:, ig), cossza(:), &
             & ztaut(:,:,ig), zasyt(:,:,ig), zomgt(:,:,ig),            &
             & zgpfd(:,:)   , zgpfu(:,:)   , ztdbt(:))
      END IF
      !
      ! Accumulate broadband (flx) clear-sky fluxes, and band-by-band (zbb) clear sky fluxes
      !
      flxu_sw(1:kproma,1:klev+1) = flxu_sw(1:kproma,1:klev+1) &
           + gpt_scaling * SPREAD(zincflx(1:kproma),DIM=2, NCOPIES=klev+1) * zgpfu(1:kproma,:)
      flxd_sw(1:kproma,1:klev+1) = flxd_sw(1:kproma,1:klev+1) &
           + gpt_scaling * SPREAD(zincflx(1:kproma),DIM=2, NCOPIES=klev+1) * zgpfd(1:kproma,:)

!IBM* ASSERT(NODEPS)
      DO jl = 1, kproma
        ib = ibs(jl,ig)
        zbbfu   (jl,ib) = zbbfu   (jl,ib) + gpt_scaling * zincflx(jl) * zgpfu(jl,klev+1)
        zbbfd   (jl,ib) = zbbfd   (jl,ib) + gpt_scaling * zincflx(jl) * zgpfd(jl,klev+1)
        zbbfddir(jl,ib) = zbbfddir(jl,ib) + gpt_scaling * zincflx(jl) * ztdbt(jl) 
      END DO

      IF(l_do_sep_clear_sky) THEN

        IF(ANY(cldMask(1:kproma,:,ig))) THEN 
          WHERE(cldMask(1:kproma,:,ig))
            ztaut(1:kproma,:,ig) = ztaur(1:kproma,:,ig) + ztaua(1:kproma,:,ig)+ ztaug(1:kproma,:,ig)
            zomgt(1:kproma,:,ig) = ztaur(1:kproma,:,ig) * 1.0_wp + ztaua(1:kproma,:,ig) * zomga(1:kproma,:,ig)
            zasyt(1:kproma,:,ig) =(zasya(1:kproma,:,ig) * zomga(1:kproma,:,ig) * ztaua(1:kproma,:,ig)) / zomgt(1:kproma,:,ig)
            zomgt(1:kproma,:,ig) = zomgt(1:kproma,:,ig) / ztaut(1:kproma,:,ig)
          END WHERE

          CALL delta_scale(kproma, kbdim, klev, ztaut(:,:,ig), zasyt(:,:,ig), zomgt(:,:,ig), update=cldMask(:,:,ig)) 
          CALL two_stream(kproma, kbdim       ,klev                         , &
               cossza(:)      ,ztaut(:,:,ig)  ,zomgt(:,:,ig)  ,zasyt(:,:,ig), &
               Rdir(:,:)      ,Rdif(:,:)      ,Tdir(:,:)      ,Tdif(:,:)    , update = cldMask(:,:,ig)) 
          !
          ! Compute fluxes using clear-sky optical properties
          ! 
          CALL srtm_solver(kproma, kbdim, klev, &
               & albdif(:,  ig), albdir(:, ig), cossza(:), ztaut(:,:,ig), &
               & Rdir(:,:), Rdif(:,:), Tdir(:,:), Tdif(:,:),              &
               & zgpcd(:,:)   , zgpcu(:,:)   , ztdbtc(:))

        ELSE 
          !
          ! Clear-sky and all-sky fluxes are the same
          !
          zgpcu(1:kproma,:) = zgpfu(1:kproma,:)
          zgpcd(1:kproma,:) = zgpfd(1:kproma,:)
        END IF
        !
        ! Accumulate broadband (flx) clear-sky fluxes
        !
        flxu_sw_clr(1:kproma,1:klev+1) = flxu_sw_clr(1:kproma,1:klev+1) &
             + gpt_scaling * SPREAD(zincflx(1:kproma),DIM=2, NCOPIES=klev+1) * zgpcu(1:kproma,1:klev+1)
        flxd_sw_clr(1:kproma,1:klev+1) = flxd_sw_clr(1:kproma,1:klev+1) &
             + gpt_scaling * SPREAD(zincflx(1:kproma),DIM=2, NCOPIES=klev+1) * zgpcd(1:kproma,1:klev+1)
      ELSE 
        !
        ! Accumulate broadband (flx) clear-sky fluxes but here we
        ! exclude cloudy subcolumns and weight to account for smaller sample size
        !
!IBM* ASSERT(NODEPS)
        DO jk = 1, klev+1
          flxu_sw_clr(1:kproma,jk) = flxu_sw_clr(1:kproma,jk)                           &
               &    + MERGE(0._wp,clrSky_scaling * zincflx(1:kproma) * zgpfu(1:kproma,jk),colCldMask(1:kproma,ig))
          flxd_sw_clr(1:kproma,jk) = flxd_sw_clr(1:kproma,jk)                           &
               &    + MERGE(0._wp,clrSky_scaling * zincflx(1:kproma) * zgpfd(1:kproma,jk),colCldMask(1:kproma,ig))
        END DO
      END IF ! Separate clear-sky calculation 

    END DO   ! Loop over spectral samples 

    ! ---  4.0 Derived calculations - diagnostics and heating rates
    ! 
    ! --------------------------------
    !
    ! Spectrally resolved fluxes of various kinds
    !
    zfvis(1:kproma,1:nbndsw) = SPREAD(bnd_wght(1:nbndsw) *           frc_vis(1:nbndsw) , DIM=1, NCOPIES=kproma)
    zfnir(1:kproma,1:nbndsw) = SPREAD(bnd_wght(1:nbndsw) * (1.0_wp - frc_vis(1:nbndsw)), DIM=1, NCOPIES=kproma)
    zfpar(1:kproma,1:nbndsw) = SPREAD(bnd_wght(1:nbndsw) * frc_par_array(1:nbndsw)     , DIM=1, NCOPIES=kproma)

    vis_frc_sfc(1:kproma) = SUM(zfvis(1:kproma,1:nbndsw)                             * &
         (zbbfd(1:kproma,1:nbndsw) - zbbfu(1:kproma,1:nbndsw)), DIM=2) / &
         (flxd_sw(1:kproma,klev+1) - flxu_sw(1:kproma,klev+1) + zepsec)

    par_dn_sfc(1:kproma) = SUM(zfpar(1:kproma,1:nbndsw)*(zbbfd(1:kproma,1:nbndsw)), DIM=2) 

    nir_dff_frc(1:kproma) = SUM(zfnir(1:kproma,1:nbndsw) * &
         (zbbfd(1:kproma,1:nbndsw) - zbbfddir(1:kproma,1:nbndsw)), DIM = 2)   / & 
         (SUM(zfnir(1:kproma,1:nbndsw) * zbbfd(1:kproma,1:nbndsw), DIM=2) + zepsec)

    vis_dff_frc(1:kproma) = SUM(zfvis(1:kproma,1:nbndsw) * &
         (zbbfd(1:kproma,1:nbndsw) - zbbfddir(1:kproma,1:nbndsw)), DIM = 2)   / &
         (SUM(zfvis(1:kproma,1:nbndsw) * zbbfd(1:kproma,1:nbndsw), DIM=2) + zepsec)

    par_dff_frc(1:kproma) = SUM(zfpar(1:kproma,1:nbndsw) * &
         (zbbfd(1:kproma,1:nbndsw) - zbbfddir(1:kproma,1:nbndsw)), DIM = 2)   / & 
         (SUM(zfpar(1:kproma,1:nbndsw) * zbbfd(1:kproma,1:nbndsw), DIM=2) + zepsec) 
    !
    ! ---  4.1 If computing clear-sky fluxes from samples, flag any columns where all samples were cloudy
    ! 
    ! --------------------------------
    IF(.not. l_do_sep_clear_sky) THEN 
!IBM* ASSERT(NODEPS)
      DO jl = 1, kproma
        IF(ALL(colCldMask(jl,:))) THEN
          flxu_sw_clr(jl,1:klev+1) = rad_undef
          flxd_sw_clr(jl,1:klev+1) = rad_undef
        END IF
      END DO
    END IF
  END SUBROUTINE srtm

END MODULE mo_srtm_driver
