!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide interface to rrtmg longwave radiation
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
MODULE mo_lrtm_driver

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: grav, amd, amw
  USE mo_rrtm_params,        ONLY: nbndlw, ngptlw
  USE mo_radiation_parameters, &
                             ONLY : i_overlap, l_do_sep_clear_sky, rad_undef
  USE mo_lrtm_setup,         ONLY: ngb, delwave
  USE rrlw_planck,           ONLY: totplanck
  USE mo_rrtm_coeffs,        ONLY: lrtm_coeffs
  USE mo_lrtm_gas_optics,    ONLY: gas_optics_lw
  USE mo_lrtm_solver,        ONLY: lrtm_solver, find_secdiff
  USE mo_cld_sampling,       ONLY: sample_cld_state
  USE mo_spec_sampling,      ONLY: spec_sampling_strategy, get_gpoint_set
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lrtm

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief Prepares information for radiation call
  !! 
  !! @remarks: This program is the driver subroutine for the longwave radiative
  !! transfer routine.  This routine is adapted from the AER LW RRTMG_LW model
  !! that itself has been adapted from RRTM_LW for improved efficiency.  Our 
  !! routine does the spectral integration externally (the solver is explicitly
  !! called for each g-point, so as to facilitate sampling of g-points
  !! This routine:
  !!    1) calls INATM to read in the atmospheric profile from GCM;
  !!       all layering in RRTMG is ordered from surface to toa. 
  !!    2) calls COEFFS to calculate various quantities needed for 
  !!       the radiative transfer algorithm.  This routine is called only once for
  !!       any given thermodynamic state, i.e., it does not change if clouds chanege
  !!    3) calls TAUMOL to calculate gaseous optical depths for each 
  !!       of the 16 spectral bands, this is updated band by band.
  !!    4) calls SOLVER (for both clear and cloudy profiles) to perform the
  !!       radiative transfer calculation with a maximum-random cloud
  !!       overlap method, or calls RTRN to use random cloud overlap.
  !!    5) passes the necessary fluxes and cooling rates back to GCM
  !!
  !
  SUBROUTINE lrtm     (kproma                     ,&
       & kbdim      ,klev     ,play     ,psfc     ,&
       & tlay       ,tlev     ,tsfc     ,wkl      ,&
       & wx         ,coldry   ,emis     ,cldfr    ,&
       & taucld     ,tauaer   ,                    &
       & rnseeds    ,strategy ,n_gpts_ts,          &   
       & uflx       ,dflx     ,uflxc    ,dflxc    )

    INTEGER, INTENT(in) :: &
         & kbdim,          & !< Maximum block length
         & kproma,         & !< Number of horizontal columns
         & klev              !< Number of model layers

    REAL(wp), INTENT(in) :: &
         & play(kbdim,klev),     & !< Layer pressures [hPa, mb] (kbdim,klev)
         & psfc(kbdim),       & !< Surface pressure [hPa, mb] (kbdim)
         & tlay(kbdim,klev),     & !< Layer temperatures [K] (kbdim,klev)
         & tlev(kbdim,klev+1),   & !< Interface temperatures [K] (kbdim,klev+1)
         & tsfc(kbdim),       & !< Surface temperature [K] (kbdim)
         & wkl(:,:,:),    & !< Gas volume mixing ratios
         & wx(:,:,:),     & !< CFC type gas volume mixing ratios
         & coldry(kbdim,klev),   & !< Column dry amount
         & emis(kbdim,nbndlw),   & !< Surface emissivity  (kbdim,nbndlw)
         & cldfr(kbdim,klev),    & !< Cloud fraction  (kbdim,klev)
         & taucld(kbdim,klev,nbndlw), & !< Coud optical depth (kbdim,klev,nbndlw)
         & tauaer(kbdim,klev,nbndlw)    !< Aerosol optical depth (kbdim,klev,nbndlw)

    ! Variables for sampling cloud state and spectral points
    INTEGER, INTENT(INOUT) :: rnseeds(:, :)    !< Seeds for random number generator (kbdim,:) 
    TYPE(spec_sampling_strategy), INTENT(IN) :: strategy
    INTEGER, INTENT(IN   ) :: n_gpts_ts


    REAL(wp), INTENT(out) :: & 
         uflx (kbdim,0:klev), & !< Tot sky longwave upward   flux [W/m2], (kbdim,0:klev)
         dflx (kbdim,0:klev), & !< Tot sky longwave downward flux [W/m2], (kbdim,0:klev)
         uflxc(kbdim,0:klev), & !< Clr sky longwave upward   flux [W/m2], (kbdim,0:klev)
         dflxc(kbdim,0:klev)    !< Clr sky longwave downward flux [W/m2], (kbdim,0:klev)

    REAL(wp) ::                  & !< Properties for one column at a time
         & taug(klev)              !< gas optical depth 
         
    
    REAL(wp) ::                            &
         & fracs(kbdim,klev,n_gpts_ts),   & !< Planck fraction per g-point
         & pwvcm(kbdim),                  & !< precipitable water vapor [cm]
         & secdiff(kbdim),                & !< diffusivity angle for RT calculation 
         & taut  (kbdim,klev,n_gpts_ts),  & !< gaseous + aerosol optical depths for all columns
         & tautot(kbdim,klev,n_gpts_ts)     !< cloud + gaseous + aerosol optical depths for all columns


    REAL(wp) ::                            & ! Properties for all bands
         & planklay(kbdim,  klev,nbndlw), & ! Planck function at mid-layer
         & planklev(kbdim,0:klev,nbndlw), & ! Planck function at level interfaces
         & plankbnd(kbdim,       nbndlw)    ! Planck function at surface

    REAL(wp) ::                    & ! Properties for a single set of columns/g-points
         & layPlnk(kbdim,  klev), & ! Planck function at mid-layer
         & levPlnk(kbdim,0:klev), & ! Planck function at level interfaces
         & bndPlnk(kbdim),        & ! Planck function at surface
         & srfEmis(kbdim)           ! Surface emission

    REAL(wp) ::                  &
         & zgpcd(kbdim,0:klev), & ! < gpoint clearsky downward flux
         & zgpcu(kbdim,0:klev), & ! < gpoint clearsky downward flux
         & zgpfd(kbdim,0:klev), & ! < gpoint fullsky downward flux
         & zgpfu(kbdim,0:klev)    ! < gpoint fullsky downward flux

    ! -----------------
    ! Variables for gas optics calculations
    INTEGER ::                    &
         & laytrop (kbdim     ), & !< tropopause layer index
         & jp      (kbdim,klev), & !< lookup table index 
         & jt      (kbdim,klev), & !< lookup table index 
         & jt1     (kbdim,klev), & !< lookup table index 
         & indself (kbdim,klev), &
         & indfor  (kbdim,klev), &
         & indminor(kbdim,klev) 

    REAL(wp) ::                       &
         & colh2o      (kbdim,klev), & !< column amount (h2o)
         & colco2      (kbdim,klev), & !< column amount (co2)
         & colo3       (kbdim,klev), & !< column amount (o3)
         & coln2o      (kbdim,klev), & !< column amount (n2o)
         & colco       (kbdim,klev), & !< column amount (co)
         & colch4      (kbdim,klev), & !< column amount (ch4)
         & colo2       (kbdim,klev), & !< column amount (o2)
         & colbrd      (kbdim,klev), & !< column amount (broadening gases)
         & selffac     (kbdim,klev), &
         & selffrac    (kbdim,klev), &
         & forfac      (kbdim,klev), &
         & forfrac     (kbdim,klev), &           
         & minorfrac   (kbdim,klev), &
         & scaleminor  (kbdim,klev), &
         & scaleminorn2(kbdim,klev), & 
         & wbrodl      (kbdim,klev) 
         
    REAL(wp) ::                       &
         & wx_loc(SIZE(wx, 2), SIZE(wx, 3)) !< Normalized CFC amounts (molecules/cm^2) 

    REAL(wp) ::      & 
         fac00(kbdim,klev),        &
         fac01(kbdim,klev),        &
         fac10(kbdim,klev),        &
         fac11(kbdim,klev)
         
    REAL(wp) ::                     & 
         rat_h2oco2  (kbdim,klev), &
         rat_h2oco2_1(kbdim,klev), &
         rat_h2oo3   (kbdim,klev), & 
         rat_h2oo3_1 (kbdim,klev), &
         rat_h2on2o  (kbdim,klev), &
         rat_h2on2o_1(kbdim,klev), &
         rat_h2och4  (kbdim,klev), &
         rat_h2och4_1(kbdim,klev), &
         rat_n2oco2  (kbdim,klev), & 
         rat_n2oco2_1(kbdim,klev), &
         rat_o3co2   (kbdim,klev), &
         rat_o3co2_1 (kbdim,klev)
    ! -----------------

    INTEGER :: jl, jk, ig ! loop indicies
    INTEGER :: ib, igpt, ibs(kbdim, n_gpts_ts), igs(kbdim, n_gpts_ts)

    REAL(wp),    PARAMETER :: cldmin = 1.e-20_wp ! minimum val for clouds

    ! Variables for sampling strategy 
    REAL(WP) :: gpt_scaling, clrSky_scaling(1:kbdim)
    REAL(WP) ::    smp_tau(kbdim, klev, n_gpts_ts) 
    LOGICAL  ::    cldMask(kbdim, klev, n_gpts_ts), & !< cloud mask in each cell
                colcldMask(kbdim,       n_gpts_ts)    !< cloud mask for each column

    !
    ! --------------------------------

    !
    ! 1.0 Choose a set of g-points to do consistent with the spectral sampling strategy
    ! 
    ! --------------------------------
    gpt_scaling = REAL(ngptlw,KIND=wp)/REAL(n_gpts_ts,KIND=wp)
    igs(1:kproma,1:n_gpts_ts) = get_gpoint_set(kproma, kbdim, strategy, rnseeds)
    
    ! Save the band nunber associated with each gpoint
    DO jl = 1, kproma  
      DO ig = 1, n_gpts_ts
        ibs(jl, ig) = ngb(igs(jl, ig))    
      END DO
    END DO  
    
    !
    ! ---  2.0 Optical properties 
    !
    ! ---  2.1 Cloud optical properties. 
    ! --------------------------------
    ! Cloud optical depth is only saved for the band associated with this g-point
    !   We sample clouds first because we may want to adjust water vapor based 
    !   on presence/absence of clouds
    !
    CALL sample_cld_state(kproma, kbdim, klev, n_gpts_ts,   &
                           rnseeds(:,:), i_overlap, & 
                           cldfr(:,:), cldMask(:,:,:))
!IBM* ASSERT(NODEPS)
    DO ig = 1, n_gpts_ts
      DO jl = 1, kproma
        smp_tau(jl,:,ig) = MERGE(taucld(jl,1:klev,ibs(jl,ig)), 0._wp, cldMask(jl,:,ig)) 
      END DO
    END DO ! Loop over samples - done with cloud optical depth calculations 
    
    !
    ! Cloud masks for sorting out clear skies - by cell and by column
    !
    IF(.not. l_do_sep_clear_sky) THEN
      !
      ! Are any layers cloudy? 
      !
	  colcldMask(1:kproma,       1:n_gpts_ts) = ANY(cldMask(1:kproma,1:klev,1:n_gpts_ts), DIM=2)
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
    ! ---  2.2. Gas optical depth calculations
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
    ! Broadening gases -- the number of molecules per cm^2 of all gases not specified explicitly 
    !   (water is excluded) 
    wbrodl(1:kproma,1:klev) = coldry(1:kproma,1:klev) - SUM(wkl(1:kproma,2:,1:klev), DIM=2)
    CALL lrtm_coeffs(kproma           ,kbdim        ,klev         ,&
         & play         ,tlay         ,coldry       ,wkl          ,&
         & wbrodl       ,laytrop      ,jp           ,jt           ,&
         & jt1          ,                                          &
         & colh2o       ,colco2       ,colo3        ,coln2o       ,&
         & colco        ,colch4       ,colo2        ,colbrd       ,&
         & fac00        ,fac01        ,fac10        ,fac11        ,&
         & rat_h2oco2   ,rat_h2oco2_1 ,rat_h2oo3    ,rat_h2oo3_1  ,&
         & rat_h2on2o   ,rat_h2on2o_1 ,rat_h2och4   ,rat_h2och4_1 ,& 
         & rat_n2oco2   ,rat_n2oco2_1 ,rat_o3co2    ,rat_o3co2_1  ,&
         & selffac      ,selffrac     ,indself      ,forfac       ,&
         & forfrac      ,indfor       ,minorfrac    ,scaleminor   ,&
         & scaleminorn2 ,indminor     )
         
      !
      !  2.2.2 Loop over g-points calculating gas optical properties. 
      !
      ! --------------------------------
!IBM* ASSERT(NODEPS)
    DO ig = 1, n_gpts_ts
      DO jl = 1, kproma  
        ib   = ibs(jl, ig) 
        igpt = igs(jl, ig) 
        !
        ! Gas concentrations in colxx variables are normalized by 1.e-20_wp in lrtm_coeffs
        !   CFC gas concentrations (wx) need the same normalization
        !   Per Eli Mlawer the k values used in gas optics tables have been multiplied by 1e20
        wx_loc(:,:) = 1.e-20_wp * wx(jl,:,:)
        CALL gas_optics_lw  (                                                            &
             & klev              ,igpt              ,play        (jl,:),wx_loc    ( :,:),&
             & coldry      (jl,:),laytrop     (jl)  ,jp          (jl,:),jt        (jl,:),&
             & jt1         (jl,:),colh2o      (jl,:),colco2      (jl,:),colo3     (jl,:),&
             & coln2o      (jl,:),colco       (jl,:),colch4      (jl,:),colo2     (jl,:),&
             & colbrd      (jl,:),fac00       (jl,:),fac01       (jl,:),fac10     (jl,:),&
             & fac11       (jl,:),rat_h2oco2  (jl,:),rat_h2oco2_1(jl,:),rat_h2oo3 (jl,:),&
             & rat_h2oo3_1 (jl,:),rat_h2on2o  (jl,:),rat_h2on2o_1(jl,:),rat_h2och4(jl,:),&
             & rat_h2och4_1(jl,:),rat_n2oco2  (jl,:),rat_n2oco2_1(jl,:),rat_o3co2 (jl,:),&
             & rat_o3co2_1 (jl,:),selffac     (jl,:),selffrac    (jl,:),indself   (jl,:),&
             & forfac      (jl,:),forfrac     (jl,:),indfor      (jl,:),minorfrac (jl,:),&
             & scaleminor  (jl,:),scaleminorn2(jl,:),indminor    (jl,:),fracs     (jl,:,ig),&
             & taug         )
        DO jk = 1, klev
          taut (jl,jk,ig) = taug(jk) + tauaer(jl,jk,ib)
        ENDDO
      END DO  ! Loop over columns
    END DO ! Loop over g point samples - done with gas optical depth calculations 

    tautot(1:kproma,:,:) = taut(1:kproma,:,:) + smp_tau(1:kproma,:,:) ! All-sky optical depth. Mask for 0 cloud optical depth? 
    
    ! 
    ! ---  3.0 Compute radiative transfer.
    ! --------------------------------
    !
    ! Initialize fluxes to zero
    !
    uflx (1:kproma,0:klev) = 0.0_wp
    dflx (1:kproma,0:klev) = 0.0_wp
    uflxc(1:kproma,0:klev) = 0.0_wp
    dflxc(1:kproma,0:klev) = 0.0_wp
    
    !
    ! Planck function in each band at layers and boundaries
    !
!IBM* ASSERT(NODEPS)
    DO ig = 1, nbndlw
      planklay(1:kproma,1:klev,ig) = planckFunction(tlay(1:kproma,1:klev  ),ig) 
      planklev(1:kproma,0:klev,ig) = planckFunction(tlev(1:kproma,1:klev+1),ig) 
      plankbnd(1:kproma,       ig) = planckFunction(tsfc(1:kproma         ),ig) 
    END DO
        
    !
    ! Precipitable water vapor in each column - this can affect the integration angle secdiff
    !
    pwvcm(1:kproma) = ((amw * SUM(wkl(1:kproma,1,1:klev), DIM=2)) / & 
                       (amd * SUM(coldry(1:kproma,1:klev) + wkl(1:kproma,1,1:klev), DIM=2))) * & 
                      (1.e3_wp * psfc(1:kproma)) / (1.e2_wp * grav)

    !
    ! Compute radiative transfer for each set of samples
    !
    DO ig = 1, n_gpts_ts
      secdiff(1:kproma) = find_secdiff(ibs(1:kproma, ig), pwvcm(1:kproma))
!IBM* ASSERT(NODEPS)
      DO jl = 1, kproma
        ib = ibs(jl,ig)
        layPlnk(jl,1:klev) = planklay(jl,1:klev,ib) 
        levPlnk(jl,0:klev) = planklev(jl,0:klev,ib)
        bndPlnk(jl       ) = plankbnd(jl,       ib)
        srfEmis(jl       ) = emis    (jl,       ib) 
      END DO 
      
      !
      ! All sky fluxes
      ! 
      CALL lrtm_solver(kproma, kbdim, klev, &
          & tautot(:,:,ig),         &               
          & layPlnk      , levPlnk, &
          & fracs(:,:,ig), secdiff, &
          & bndPlnk      , srfEmis, & 
          & zgpfu    , zgpfd)
       
      uflx (1:kproma,0:klev) = uflx (1:kproma,0:klev) &
                             + zgpfu(1:kproma,0:klev) * gpt_scaling
      dflx (1:kproma,0:klev) = dflx (1:kproma,0:klev) &
                             + zgpfd(1:kproma,0:klev) * gpt_scaling
	                           
      !
      ! Clear-sky fluxes
      !
      IF(l_do_sep_clear_sky) THEN
        !
        ! Remove clouds and do second RT calculation
        !
        CALL lrtm_solver(kproma, kbdim, klev, &
            & taut (:,:,ig),                  &
            & layPlnk      , levPlnk,         &
            & fracs(:,:,ig), secdiff,         &
            & bndPlnk      , srfEmis,         & 
            & zgpcu    , zgpcd)
        uflxc(1:kproma,0:klev) = uflxc(1:kproma,0:klev) + zgpcu(1:kproma,0:klev) * gpt_scaling
        dflxc(1:kproma,0:klev) = dflxc(1:kproma,0:klev) + zgpcd(1:kproma,0:klev) * gpt_scaling
      ELSE
        !
        ! Accumulate fluxes by excluding cloudy subcolumns, weighting to account for smaller sample size
        !
!IBM* ASSERT(NODEPS)
        DO jk = 0, klev
		  uflxc(1:kproma,jk) = uflxc(1:kproma,jk)                             &
							 + MERGE(0._wp,                                   &
									 zgpfu(1:kproma,jk) * clrSky_scaling(1:kproma), &
									 colCldMask(1:kproma,ig))
		  dflxc(1:kproma,jk) = dflxc(1:kproma,jk)                             &
							 + MERGE(0._wp,                                   &
									 zgpfd(1:kproma,jk) * clrSky_scaling(1:kproma), &
									 colCldMask(1:kproma,ig))
        END DO 
      END IF 
    END DO ! Loop over samples

    !
    ! ---  3.1 If computing clear-sky fluxes from samples, flag any columns where all samples were cloudy
    ! 
    ! --------------------------------
	IF(.not. l_do_sep_clear_sky) THEN 
!IBM* ASSERT(NODEPS)
	  DO jl = 1, kproma
	    IF(ALL(colCldMask(jl,:))) THEN
	      uflxc(jl,0:klev) = rad_undef
	      dflxc(jl,0:klev) = rad_undef
	    END IF
	  END DO 
	END IF 
  END SUBROUTINE lrtm

  !----------------------------------------------------------------------------
  ELEMENTAL FUNCTION planckFunction(temp, band)
    !
    ! Compute the blackbody emission in a given band as a function of temperature
    !
    REAL(WP), INTENT(IN) :: temp
    INTEGER,  INTENT(IN) :: band 
    REAL(WP)             :: planckFunction
    
    INTEGER  :: index
    REAL(WP) :: fraction 
    
    index = MIN(MAX(1, INT(temp - 159._wp)),180)
    fraction = temp - 159._wp - float(index)
    
    planckFunction = totplanck(index, band) &
                   + fraction * (totplanck(index+1, band) - totplanck(index, band))
    planckFunction = planckFunction * delwave(band)
  END FUNCTION planckFunction

END MODULE mo_lrtm_driver

