!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_surface
  !---------------------------------------------------------------------
  ! USE statements 
  !---------------------------------------------------------------------
  USE mo_kind,                 ONLY: wp
  USE mo_io_units,             ONLY: nout
  USE mo_surface_types
  USE mo_surface_memory,       ONLY: land
  USE mo_surface_memory,       ONLY: ocean
  USE mo_surface_memory,       ONLY: ice 
  USE mo_surface_memory,       ONLY: box
  USE mo_surface_memory,       ONLY: surface
  USE mo_jsbach_interface,     ONLY: jsbach_inter_1d
! for station diagnostic
  USE mo_memory_g3b,           ONLY: evap_na, ahfl_na, ahfs_na, wind10_na, &
                                     ustr_na, vstr_na, tsurf_na
  
  !---------------------------------------------------------------------
  IMPLICIT NONE 
  ! 
  !---------------------------------------------------------------------
CONTAINS
  !---------------------------------------------------------------------
  ! SUBROUTINES: update_surface is one and only routine here 
  !---------------------------------------------------------------------
  SUBROUTINE update_surface (                         & ! logistic parameters
       current_nproma,                                & ! vector length
       jrow,                                          & ! krow
       levels,                                        & ! number of levels (lowest level)
       levels_plus_1,                                 & ! one up
       levels_minus_1,                                & ! one down
       lpnorth,                                       & ! true in northern hemisphere
  !!-----------------------------!! atmospheric conditions lowest level-
       cloud_water,                                   & ! total amount of water in liquid droplets
       cloud_ice,                                     & ! total amount of water in cloud ice
       geopotential,                                  & ! geopot.
       atm_temp,                                      & ! atmosphere temperature (lowest layer)
       atm_spec_humidity,                             & ! specific humidity atmos. (lowest layer)
       atm_full_lev_press,                            & ! full level pressure
       atm_half_lev_press,                            & ! half level pressure
       wind_u,                                        & ! wind speed zonal
       wind_v,                                        & ! wind speed meridional
       rain,                                          & ! total rain fall
       snow,                                          & ! total snow fall
       longwave_net,                                  & ! long wave net radiation
       sw_vis,                                        & ! net surface visible
       sw_vis_frac_diffuse,                           & ! fraction of diffuse visible
       sw_nir,                                        & ! net surface NIR
       sw_nir_frac_diffuse,                           & ! fraction of diffuse NIR
       sw_par,                                        & ! downward surface PAR
       sw_par_frac_diffuse,                           & ! fraction of diffuse PAR
       cos_zenith_angle,                              & ! zenith angle for rad. between l_trigrad
       cloud_cover,                                   & ! integrated cloud cover [0,...,1]
       !!-----------------------------------!! -------------------------
       ocean_temp,                                    & ! SST is read by echam and handed over here
       ocean_u,                                       & ! ocean_u_velocity
       ocean_v,                                       & ! ocean_v_velocity
       ice_depth,                                     & ! ICE depth comes from clsst 
       pameltdepth,                                   & ! melt pond depth
       pameltfrac,                                    & ! melt pond fraction
       pamlcorr,                                      & ! mixed layer ocean
       pamlcorac,                                     & ! mixed layer ocean
       pamlheatac,                                    & ! mixed layer ocean
       !!-----------------------------------!! -------------------------
       zcfh,                 zebsh,                   & !
       zqdif_pre,            ztdif_pre,               & !
       zudif,                zvdif,                   & !
       zghabl,               pi0,                     & !       
       ptrsol,                                        & !
       pco2_concentration,                            & ! CO2 concentration in lowest level 
                                                        ! (mass mixing ratio)
       !!-------OUTPUT--------------------------------------------------
       palbedo,              palbedo_vis,             &
       palbedo_nir,                                   &
       palbedo_vis_dir,      palbedo_nir_dir,         &
       palbedo_vis_dif,      palbedo_nir_dif,         &
       ptrflw,               ptrfli,                  &
       psofll,               psoflw,                  &
       psofli,               ptrfllac,                &
       ptrflwac,             ptrfliac,                &
       psofllac,             psoflwac,                &
       psofliac,             palsol,                  &
       palsoi,               palsow,                  &
       palsobs,              palsom,                  &
       ptaus,                pfage,                   &
       psnifrac,             pbarefrac,               &
       ustarm,               momentum_exchange_coeff, &
       heat_exchange_coeff,                           &
       tkevn_cond,           ztdif_new,               &
       zqdif_new,            ztvh,                    &
       zqsurf,               zth,                     &
       pwind10w ,            pu10 ,                   &
       pv10 ,                pwimax ,                 &
       pwind10,              pdew2,                   &
       ptemp2,               pt2max,                  &
       pt2min,               pevaplac,                &
       pevapwac,             pevapiac,                &
       pevap,                pahfllac,                &
       pahflwac,             pahfliac,                &
       pahfl,                pahfslac,                &
       pahfswac,             pahfsiac,                &
       pahfs,                pqhfla,                  &
       pevapw,                                        &
       pevapi,               pahfsl,                  &
       pahfsw,               pahfsi,                  &
       pevapot,                                       &
       pahflw,                                        &
       pahfli,               pahfll,                  &
       psni,                 pahfice,                 &
       pfluxres,             pqres,                   &
       pahfcon,              pahfres,                 &
       ptsi,                 ptslnew,                 &
       pzti,                 pzteffl4,                &
       pztsnew,              ptsurf,                  &
       paz0w,                paz0i,                   &
       paz0l,                paz0,                    &
       pustrl,               pvstrl,                  &
       pustrw,               pvstrw,                  &
       pustri,               pvstri,                  &
       pustr,                pvstr,                   &
       ptte_corr,            ptsw_new,                &
       pseaice,              pseaice_new,             &
       psiced_new,           pradtemp_old,            &
#ifdef HAMMOZ /* SF */
!SF gf #78 add some additional separate land/sea/ice diags in case of HAMMOZ
              zcfml,  zcfmw,  zcfmi,                  &
       zcfnc, zcfncl, zcfncw, zcfnci,                 &
       zri,   zril,   zriw,   zrii,                   &
              ztvl,   ztvw,   ztvi,                   &
              pfrl,   pfrw,   pfri,                   &
       palake,                                        &
       pcvs,  zsrfll,                                 &
       zcdn,  zcdnl,  zcdnw,  zcdni,                  &
       velo10m,                                       &
#else
       zcfnc, zri, ztvl, pfrl, pfrw, pfri, palake,    &
       pcvs, zsrfll, zcdn, velo10m,                   & 
#endif
       pco2_flux_ocean,                               &
       pco2_flux_land,                                &
       pco2_flux,                                     &
       pco2_flux_npp, pco2_flux_soilresp,             &
       pco2_flux_herbivory, pco2_flux_dynveg,         &
       pco2_emission_lcc, pco2_emission_harvest)

  !---------------------------------------------------------------------
  ! USE statements 
  !---------------------------------------------------------------------
    USE mo_memory_g3b, ONLY: &
         slf, slm,           &   !! Sea-land fraction (contains fraction of land [0,....,1]
         friac, tslm1, ws, sn, wl, &
         gld, wsmx

    USE mo_surface_land
    USE mo_surface_ocean
    USE mo_surface_ice
    USE mo_surface_boundary
    USE mo_time_control,         ONLY: lstart, delta_time
    USE mo_control,              ONLY: lmeltpond, lmlo, lcouple
    USE mo_submodel,             ONLY: lanysubmodel
    USE mo_vphysc,               ONLY: vphysc ! smelt for vphysc stream (dry deposition modules)
    !-------------------------------------------------------------------
    !                         :: NEW VARIABLE NAME          :: OLD ECHAM NAME AND DIMENSIONS
    !-------------------------------------------------------------------
    INTEGER, INTENT(in)       :: current_nproma             !! kbdim
    INTEGER, INTENT(IN)       :: jrow                       !! krow
    INTEGER, INTENT(IN)       :: levels                     !! klev
    INTEGER, INTENT(in)       :: levels_plus_1              !! klevp1
    INTEGER, INTENT(in)       :: levels_minus_1             !! klevm1
    LOGICAL, INTENT(in)       :: lpnorth(:)                 !! true in northern hemisphere
    REAL(wp),    INTENT(in)       :: geopotential(:)            !! pgeom1(kdim,klev)
    REAL(wp),    INTENT(in)       :: wind_u(:)                  !! pum1
    REAL(wp),    INTENT(in)       :: wind_v(:)                  !! pvm1
    REAL(wp),    INTENT(in)       :: atm_temp(:)                !! ptm1
    REAL(wp),    INTENT(in)       :: atm_spec_humidity(:)       !! pqm1
    REAL(wp),    INTENT(in)       :: rain(:)                    !! rsfl + rsfc from physc
    REAL(wp),    INTENT(in)       :: snow (:)                   !! ssfl + ssfc from physc
    REAL(wp),    INTENT(in)       :: longwave_net(:)            !! pemter
    REAL(wp),    INTENT(in)       :: sw_vis(:)                  !! net surface visible
    REAL(wp),    INTENT(in)       :: sw_vis_frac_diffuse(:)     !! fraction of diffuse visible
    REAL(wp),    INTENT(in)       :: sw_nir(:)                  !! net surface NIR
    REAL(wp),    INTENT(in)       :: sw_nir_frac_diffuse(:)     !! fraction of diffuse visible
    REAL(wp),    INTENT(in)       :: sw_par(:)                  !! downward surface PAR
    REAL(wp),    INTENT(in)       :: sw_par_frac_diffuse(:)     !! fraction of diffuse PAR
    REAL(wp),    INTENT(in)       :: atm_full_lev_press(:)      !! papm1
    REAL(wp),    INTENT(in)       :: cos_zenith_angle(:)        !! zm0
    REAL(wp),    INTENT(in)       :: atm_half_lev_press(:,:)    !! paphm1(kbdim, klevp1)
    REAL(wp),    INTENT(in)       :: cloud_water(:)             !! pxlm1
    REAL(wp),    INTENT(in)       :: cloud_ice(:)               !! pxim1
    REAL(wp),    INTENT(in)       :: cloud_cover(:)             !! paclc
    REAL(wp),    INTENT(in)       :: ocean_temp(:)              !! ptsw
    REAL(wp),    INTENT(in)       :: ocean_u(:)                 !! ocu
    REAL(wp),    INTENT(in)       :: ocean_v(:)                 !! ocv
    REAL(wp),    INTENT(in)       :: zcfh(:,:)                  !! zcfh: dimless exch coeff. from 
                                                                !!  diffusion scheme
    REAL(wp),    INTENT(in)       :: zebsh(:,:)                 !! zebsh from heat diffusion scheme
    REAL(wp),    INTENT(in)       :: zqdif_pre(:,:)             !! zqdif from heat diffusion scheme
    REAL(wp),    INTENT(in)       :: ztdif_pre(:,:)             !! ztdif from heat diffusion scheme
    REAL(wp),    INTENT(in)       :: zudif(:)                   !! zudif from momentum transfer 
                                                                !!  scheme
    REAL(wp),    INTENT(in)       :: zvdif(:)                   !! zvdif from momentum transfer 
                                                                !!  scheme
    REAL(wp),    INTENT(in)       :: zghabl(:)                  !! zghabl from pbl-height
    REAL(wp),    INTENT(in)       :: pi0(:)                     !! solar incidence factor for solar
                                                                !!  radiation
    REAL(wp),    INTENT(in)       :: ptrsol(:)                  !! ptrsol solar radiation
    REAL(wp),    INTENT(in)       :: ice_depth(:)               !! psiced comes from clsst
    REAL(wp),    INTENT(in)       :: pseaice(:)                 !! seaice comes from clsst
    REAL(wp),    INTENT(in)       :: pameltdepth(:)             !! melt pond depth
    REAL(wp),    INTENT(in)       :: pameltfrac(:)              !! melt pond fraction
    REAL(wp),    INTENT(inout)    :: pamlcorr(:)                !! mixed layer ocean
    REAL(wp),    INTENT(inout)    :: pamlcorac(:)               !! mixed layer ocean
    REAL(wp),    INTENT(inout)    :: pamlheatac(:)              !! mixed layer ocean
    REAL(wp),    INTENT(in)       :: pco2_concentration(:)      !! CO2 tracer
    !! - OUTPUT --------------------------------------------------------
    REAL(wp),    INTENT(out)      :: ustarm(:)                  !! ref. to zustarm in vdiff (341)
    REAL(wp),    INTENT(out)      :: momentum_exchange_coeff(:) !! zcdum for lowest layer-soil
    REAL(wp),    INTENT(out)      :: heat_exchange_coeff(:)     !! cduh for lowest layer-soil
    REAL(wp),    INTENT(out)      :: tkevn_cond(:)              !! ptkevn
    REAL(wp),    INTENT(out)      :: ztdif_new(:)               !! ztdif
    REAL(wp),    INTENT(out)      :: zqdif_new(:)               !! zqdif
    REAL(wp),    INTENT(out)      :: ztvh(:)                    !! ztvh
    REAL(wp),    INTENT(out)      :: zqsurf(:)                  !! zqsl*zcsat*pfrl+pfrw*zqsw
                                                                !!  +pfri*zqsi
    REAL(wp),    INTENT(out)      :: zth(:)
    REAL(wp),    INTENT(out)      :: pwind10w(:)
    REAL(wp),    INTENT(out)      :: pu10(:)
    REAL(wp),    INTENT(out)      :: pv10(:)
    REAL(wp),    INTENT(out)      :: pwimax(:)
    REAL(wp),    INTENT(out)      :: pwind10(:)
    REAL(wp),    INTENT(out)      :: pdew2(:)
    REAL(wp),    INTENT(out)      :: ptemp2(:)
    REAL(wp),    INTENT(out)      :: pt2max(:)
    REAL(wp),    INTENT(out)      :: pt2min(:)
    REAL(wp),    INTENT(out)      :: pevaplac(:)
    REAL(wp),    INTENT(out)      :: pevapwac(:)
    REAL(wp),    INTENT(out)      :: pevapiac(:)
    REAL(wp),    INTENT(out)      :: pevap(:)
    REAL(wp),    INTENT(out)      :: pahfllac(:)
    REAL(wp),    INTENT(out)      :: pahflwac(:)
    REAL(wp),    INTENT(out)      :: pahfliac(:)
    REAL(wp),    INTENT(out)      :: pahfl(:)
    REAL(wp),    INTENT(out)      :: pahfslac(:)
    REAL(wp),    INTENT(out)      :: pahfswac(:)
    REAL(wp),    INTENT(out)      :: pahfsiac(:)
    REAL(wp),    INTENT(out)      :: pahfs(:)
    REAL(wp),    INTENT(out)      :: pqhfla(:)
    REAL(wp),    INTENT(out)      :: pevapw(:)
    REAL(wp),    INTENT(out)      :: pevapi(:)
    REAL(wp),    INTENT(out)      :: pahfsl(:)
    REAL(wp),    INTENT(out)      :: pahfsw(:)
    REAL(wp),    INTENT(out)      :: pahfsi(:)
    REAL(wp),    INTENT(out)      :: pevapot(:)
    REAL(wp),    INTENT(out)      :: pahflw(:)
    REAL(wp),    INTENT(out)      :: pahfli(:)
    REAL(wp),    INTENT(out)      :: pahfll(:)
    REAL(wp),    INTENT(inout)    :: psni(:)
    REAL(wp),    INTENT(out)      :: pahfice(:)
    REAL(wp),    INTENT(out)      :: pfluxres(:)
    REAL(wp),    INTENT(out)      :: pqres(:)
    REAL(wp),    INTENT(out)      :: pahfcon(:)
    REAL(wp),    INTENT(out)      :: pahfres (:)
    REAL(wp),    INTENT(out)      :: ptsi (:)
    REAL(wp),    INTENT(out)      :: ptslnew(:)
    REAL(wp),    INTENT(out)      :: pzti(:)
    REAL(wp),    INTENT(out)      :: pzteffl4(:)
    REAL(wp),    INTENT(out)      :: pztsnew(:)
    REAL(wp),    INTENT(out)      :: ptsurf(:)
    REAL(wp),    INTENT(out)      :: paz0w(:)
    REAL(wp),    INTENT(out)      :: paz0i(:)
    REAL(wp),    INTENT(out)      :: paz0l(:)
    REAL(wp),    INTENT(out)      :: paz0(:)
    REAL(wp),    INTENT(out)      :: pustrl(:)
    REAL(wp),    INTENT(out)      :: pvstrl(:)
    REAL(wp),    INTENT(out)      :: pustrw(:)
    REAL(wp),    INTENT(out)      :: pvstrw(:)
    REAL(wp),    INTENT(out)      :: pustri(:)
    REAL(wp),    INTENT(out)      :: pvstri(:)
    REAL(wp),    INTENT(out)      :: pustr(:)
    REAL(wp),    INTENT(out)      :: pvstr(:)
    REAL(wp),    INTENT(out)      :: ptte_corr(:)
    REAL(wp),    INTENT(out)      :: ptsw_new(:)
    REAL(wp),    INTENT(out)      :: pseaice_new(:)
    REAL(wp),    INTENT(out)      :: psiced_new(:)
    REAL(wp),    INTENT(out)      :: pradtemp_old(:)
    REAL(wp),    INTENT(out)      :: palbedo(:)
    REAL(wp),    INTENT(out)      :: palbedo_vis(:)
    REAL(wp),    INTENT(out)      :: palbedo_nir(:)
    REAL(wp),    INTENT(inout)    :: palbedo_vis_dir(:)
    REAL(wp),    INTENT(inout)    :: palbedo_nir_dir(:)
    REAL(wp),    INTENT(inout)    :: palbedo_vis_dif(:)
    REAL(wp),    INTENT(inout)    :: palbedo_nir_dif(:)
    REAL(wp),    INTENT(out)      :: palsol(:) 
    REAL(wp),    INTENT(inout)      :: palsoi(:)  
    REAL(wp),    INTENT(inout)      :: palsow(:)
    REAL(wp),    INTENT(inout)      :: palsobs(:)
    REAL(wp),    INTENT(inout)      :: palsom(:)
    REAL(wp),    INTENT(inout)      :: ptaus(:)
    REAL(wp),    INTENT(inout)      :: pfage(:)
    REAL(wp),    INTENT(inout)      :: psnifrac(:)
    REAL(wp),    INTENT(inout)      :: pbarefrac(:)
    REAL(wp),    INTENT(out)      :: ptrflw(:)
    REAL(wp),    INTENT(out)      :: ptrfli(:) 
    REAL(wp),    INTENT(out)      :: psofll(:)
    REAL(wp),    INTENT(out)      :: psoflw(:)
    REAL(wp),    INTENT(out)      :: psofli(:)               
    REAL(wp),    INTENT(out)      :: ptrfllac(:)
    REAL(wp),    INTENT(out)      :: ptrflwac(:)
    REAL(wp),    INTENT(out)      :: ptrfliac(:)             
    REAL(wp),    INTENT(out)      :: psofllac(:)
    REAL(wp),    INTENT(out)      :: psoflwac(:)
    REAL(wp),    INTENT(out)      :: psofliac(:)  
    REAL(wp),    INTENT(in)       :: pco2_flux_ocean(:)         
    REAL(wp),    INTENT(out)      :: pco2_flux_land(:), &
                                     pco2_flux_npp(:), pco2_flux_soilresp(:), &
                                     pco2_flux_herbivory(:), pco2_flux_dynveg(:)
    REAL(wp),    INTENT(out)      :: pco2_emission_lcc(:), pco2_emission_harvest(:)             
    REAL(wp),    INTENT(out)      :: pco2_flux(:)             

!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------

!>>SF gf #78 add some additional separate land/sea/ice diags in case of HAMMOZ
#ifdef HAMMOZ
    REAL(wp),    INTENT(out)      :: zcfml(:)      !< mom. exch. coeff. over land 
    REAL(wp),    INTENT(out)      :: zcfmw(:)      !< mom. exch. coeff. over water       
    REAL(wp),    INTENT(out)      :: zcfmi(:)      !< mom. exch. coeff. over ice
    REAL(wp),    INTENT(out)      :: zcfncl(:)     !< func. of the heat transfer coeff. over land
    REAL(wp),    INTENT(out)      :: zcfncw(:)     !< func. of the heat transfer coeff. over water
    REAL(wp),    INTENT(out)      :: zcfnci(:)     !< func. of the heat transfer coeff. over ice
    REAL(wp),    INTENT(inout)    :: zril(:)       !< Richardson number over land
    REAL(wp),    INTENT(inout)    :: zriw(:)       !< Richardson number over water      
    REAL(wp),    INTENT(inout)    :: zrii(:)       !< Richardson number over ice
    REAL(wp),    INTENT(out)      :: ztvw(:)       !< virtual potential temp. over water
    REAL(wp),    INTENT(out)      :: ztvi(:)       !< virtual potential temp. over ice
    REAL(wp),    INTENT(out)      :: zcdnl(:)      
    REAL(wp),    INTENT(out)      :: zcdnw(:)             
    REAL(wp),    INTENT(out)      :: zcdni(:)             
#endif
!<<SF gf #78
    REAL(wp),    INTENT(out)      :: zcfnc(:)             
    REAL(wp),    INTENT(out)      :: zri(:)
    REAL(wp),    INTENT(out)      :: ztvl(:)
    REAL(wp),    INTENT(in)      :: pfrl(:)
    REAL(wp),    INTENT(in)      :: pfrw(:)
    REAL(wp),    INTENT(in)      :: pfri(:)
    REAL(wp),    INTENT(in)      :: palake(:)
    REAL(wp),    INTENT(out)      :: pcvs(:)
    REAL(wp),    INTENT(out)      :: zsrfll(:)
    REAL(wp),    INTENT(out)      :: zcdn(:)
    REAL(wp),    INTENT(out)      :: velo10m(:)
!>>SF gf #78 do not define these local vars in case of HAMMOZ because they are already arguments of the subroutine
#ifndef HAMMOZ
    REAL(wp)                      :: ztvw  (current_nproma)
    REAL(wp)                      :: ztvi  (current_nproma)
    REAL(wp)                      :: zcdnw (current_nproma)
    REAL(wp)                      :: zcdni (current_nproma)
    REAL(wp)                      :: zcfncw(current_nproma)
    REAL(wp)                      :: zcfnci(current_nproma)
#endif
!<<SF gf #78

!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! LOCAL VARIABLES
    !-------------------------------------------------------------------
    REAL(wp)       :: zqdif(current_nproma, levels)
    REAL(wp)       :: ztdif(current_nproma, levels)
    REAL(wp)       :: radtemp(current_nproma)
    REAL(wp)       :: longwave_down(current_nproma)
    REAL(wp)       :: radtemp_old(current_nproma)
    REAL(wp)       :: soil_wetness(current_nproma)
    REAL(wp)       :: snow_depth(current_nproma)
    REAL(wp)       :: skin_reservoir(current_nproma)
    REAL(wp)       :: tte_corr(current_nproma)
!>>SF gf #78 do not define these local vars in case of HAMMOZ because they are already arguments of the subroutine
#ifndef HAMMOZ
    REAL(wp)       :: zcdnl(current_nproma)
    REAL(wp)       :: zcfncl(current_nproma)
#endif
!<<SF gf #78
    REAL(wp)       :: zsnow_melt(current_nproma)
    REAL(dp)       :: sw_vis_land(current_nproma)
    REAL(dp)       :: sw_nir_land(current_nproma)
    !-------------------------------------------------------------------
    ! MASK FOR MASKOUT LAND,OCEAN,ICE
    !-------------------------------------------------------------------
    REAL(wp) ::  zero(current_nproma)
    INTEGER  ::  nproma
#if (! defined (__prism)) && (defined (__PGI))
    INTEGER  :: np
#endif
    !-------------------------------------------------------------------
    ! GET CURRENT VECTOR AND ROW DEFINITION OF DOMAIN 
    !-------------------------------------------------------------------
    nproma   = current_nproma
    zero(:)  = 0._wp
    !-------------------------------------------------------------------
    ! INITIALISATION NEW START OF MODEL 
    !-------------------------------------------------------------------
    IF (lstart) THEN
       ! New start of model
       land%wind_10_meter(1:nproma,jrow)=                              &
                                       SQRT(wind_u(:)**2 + wind_v(:)**2)
    END IF
    !-------------------------------------------------------------------
    ! PERMANENT INITIALISATION 
    !-------------------------------------------------------------------
    zqdif(:,:)                                  = zqdif_pre(:,:)
    ztdif(:,:)                                  = ztdif_pre(:,:)
    ocean%surface_temperature(1:nproma,jrow)    = ocean_temp(:) ! comes from clsst every time step
    box%seaice(1:nproma,jrow)                   = pseaice(:)    ! comes from clsst every time step
    ice%ice_depth(1:nproma,jrow)                = ice_depth(:)  ! comes from clsst every time step
    ice%snow_water_equivalent(1:nproma,jrow)    = psni(:)       ! comes from ocean every time step
                                                                !  in coupled mode (only over ice)
    ! FRACTIONAL COVER and MASK DEFINITION
    ! Note: seaice cover changes over time, therefore this isn't in init_surface
    !-------------------------------------------------------------------
    box%land_fract(1:nproma,jrow)   = pfrl(1:nproma)
    box%ocean_fract(1:nproma,jrow)  = pfrw(1:nproma)
    box%seaice_fract(1:nproma,jrow) = pfri(1:nproma)
    surface%is_land(1:nproma,jrow)  = pfrl(1:nproma) > 0.0_wp

    box%frac_ice_cover_acc(1:nproma,jrow)   =                          &
                        box%frac_ice_cover_acc(1:nproma,jrow)          &
                        + delta_time * box%seaice_fract(1:nproma,jrow)
    surface%is_seaice(1:nproma,jrow) = box%seaice_fract(1:nproma,jrow) > 0.0_wp
    surface%is_ocean (1:nproma,jrow) = box%ocean_fract (1:nproma,jrow) > 0.0_wp

    WHERE (.NOT. surface%is_land(1:nproma,jrow))  
       ice%zcvsi(1:nproma,jrow) =                                      &
                 TANH(ice%snow_water_equivalent(1:nproma,jrow) * 100._wp)
    ELSEWHERE
       ice%zcvsi(1:nproma,jrow) = 0._wp
    END WHERE
    !-------------------------------------------------------------------
    ! additional grid points for ocean coupling

    surface%is_ocean_for_coup(1:nproma,jrow) = .false.
    surface%is_ice_for_coup(1:nproma,jrow)   = .false.
    IF(lcouple) THEN
     surface%is_ocean_for_coup(1:nproma,jrow) =                         &
          slf(1:nproma,jrow) .lt. 1.0_wp                                &
          .and.(.not. surface%is_ocean(1:nproma,jrow))                  
     surface%is_ice_for_coup(1:nproma,jrow)   =                         &
          (surface%is_ocean_for_coup(1:nproma,jrow)  .or.               &
           surface%is_ocean(1:nproma,jrow)         ) .and.              &
          (.not. surface%is_seaice(1:nproma,jrow))
     surface%is_ocean (1:nproma,jrow) = surface%is_ocean (1:nproma,jrow)&
                     .OR. surface%is_ocean_for_coup(1:nproma,jrow)
     surface%is_seaice(1:nproma,jrow) = surface%is_seaice(1:nproma,jrow)&
                     .OR. surface%is_ice_for_coup(1:nproma,jrow)
    END IF

    !-------------------------------------------------------------------
    ! long wave radiation
    CALL longwave_down_rad( nproma                                     &
                          , longwave_net(1:nproma)                     &
                          , land%surface_temperature(1:nproma,jrow)    &
                          , ocean%surface_temperature(1:nproma,jrow)   &
                          , ice%surface_temperature(1:nproma,jrow)     &
                          , box%land_fract(1:nproma,jrow)              &
                          , box%ocean_fract(1:nproma,jrow)             &
                          , box%seaice_fract(1:nproma,jrow)            & 
                          , longwave_down(1:nproma) )
         
    box%longwave_down_acc(1:nproma,jrow)   =                           &
                        box%longwave_down_acc(1:nproma,jrow)           &
                        + delta_time * longwave_down(1:nproma)

    ! ATMOSPHERIC LOWEST LEVEL PARAMERTERS
    !-------------------------------------------------------------------
    CALL atm_conditions( nproma                                        &
                       , cloud_water(:)                                &
                       , cloud_ice(:)                                  &
                       , geopotential(:)                               &
                       , atm_temp(:)                                   &
                       , atm_spec_humidity(:)                          &
                       , atm_full_lev_press(:)                         &
                       , box%atm_tot_cloud_water(1:nproma,jrow)        &
                       , box%atm_dry_stat_energy(1:nproma,jrow)        &
                       , box%atm_pot_temp(1:nproma,jrow)               &
                       , box%atm_vir_pot_temp(1:nproma,jrow)           &
                       , box%atm_lat_heat_fact(1:nproma,jrow)          &
                       , box%atm_liq_wat_pot_temp(1:nproma,jrow)       &
                       , box%atm_sat_spec_hum(1:nproma,jrow)          )

    ! calculate surface and lowest level parameters for grid 
    ! and grid fractions of land, water, ice
    !-------------------------------------------------------------------
    CALL precalc_land( nproma                                          &
                     , surface%is_land(1:nproma,jrow)                  &  !! in
                     , land%surface_temperature(1:nproma,jrow)         &
                     , atm_half_lev_press(:,levels_plus_1)             &  !! in
                     , wind_u(:)                                       &
                     , wind_v(:)                                       &  !! in
                     , atm_spec_humidity(:)                            &
                     , box%atm_tot_cloud_water(1:nproma,jrow)          &  !! in
                     , box%atm_sat_spec_hum(1:nproma,jrow)             &
                     , box%atm_pot_temp(1:nproma,jrow)                 &  !! in
                     , box%atm_vir_pot_temp(1:nproma,jrow)             &
                     , box%atm_lat_heat_fact(1:nproma,jrow)            &  !! in
                     , cloud_cover(:)                                  &
                     , box%atm_liq_wat_pot_temp(1:nproma,jrow)         &  !! in
                     , geopotential(:)                                 &
                     , land%roughness_heat(1:nproma,jrow)              &  !! in
                     , land%roughness_momentum(1:nproma,jrow)          &  !! in
                     , atm_temp(:)                                     &
                     , zghabl(:)                                       &  !! in
                     , land%ZDQSL(1:nproma,jrow)                       &
                     , land%ZRIL(1:nproma,jrow)                        &  !! out, in  
                     , land%ZQSL(1:nproma,jrow)                        &
                     , zcfncl(:)                                       &  !! out, out
                     , land%ZCHL(1:nproma,jrow)                        &
                     , land%ZCFHL(1:nproma,jrow)                       &  !! out  
                     , land%zbnl(1:nproma,jrow)                        &
                     , land%zbhnl(1:nproma,jrow)                       &  !! out
                     , land%zbml(1:nproma,jrow)                        &
                     , land%zbhl(1:nproma,jrow)                        &  !! out
                     , land%zustarl(1:nproma,jrow)                     &
                     , land%ztkevl(1:nproma,jrow)                      &  !! out
                     , land%zcfml(1:nproma,jrow)                       &  !! in, out
                     , land%zustl(1:nproma,jrow)                       &
                     , land%zcptl(1:nproma,jrow)                       &
                     , land%zcpq(1:nproma,jrow)                        &  !! out
                     , land%zcair(1:nproma,jrow)                       &
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
                     , ztvl(1:nproma)                                  &  
                     , zcdnl(1:nproma)                                 &
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------
                     , land%zcsat(1:nproma,jrow)                      )

    CALL precalc_ocean( nproma                                         &
                      , surface%is_ocean(1:nproma,jrow)                &
                      , ocean%surface_temperature(1:nproma,jrow)       &
                      , atm_half_lev_press(:,levels_plus_1)            &
                      , atm_spec_humidity(:)                           &
                      , box%atm_tot_cloud_water(1:nproma,jrow)         &
                      , atm_temp(:)                                    &
                      , box%atm_sat_spec_hum(1:nproma,jrow)            &
                      , box%atm_pot_temp(1:nproma,jrow)                &
                      , box%atm_vir_pot_temp(1:nproma,jrow)            &
                      , box%atm_lat_heat_fact(1:nproma,jrow)           &
                      , cloud_cover(:)                                 &
                      , box%atm_liq_wat_pot_temp(1:nproma,jrow)        &
                      , wind_u(:)                                      &
                      , wind_v(:)                                      &
                      , ocean_u(:)                                     &
                      , ocean_v(:)                                     &
                      , ocean%roughness(1:nproma,jrow)                 &
                      , geopotential(:)                                &
                      , zghabl(:)                                      &
                      , ocean%zqsw(1:nproma,jrow)                      &
                      , ocean%zcptw(1:nproma,jrow)                     &
                      , ocean%zriw(1:nproma,jrow)                      &
                      , ocean%zcfhw(1:nproma,jrow)                     &
                      , ocean%zchw(1:nproma,jrow)                      &
                      , ocean%zbnw(1:nproma,jrow)                      &
                      , ocean%zbmw(1:nproma,jrow)                      &
                      , ocean%zbhw(1:nproma,jrow)                      &
                      , ocean%zustarw(1:nproma,jrow)                   &
                      , ocean%ztkevw(1:nproma,jrow)                    &
                      , ocean%zcfmw(1:nproma,jrow)                     &
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
                     , zcfncw(1:nproma)                                &
                     , ztvw(1:nproma)                                  &
                     , zcdnw(1:nproma)                                 &
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------
                      , ocean%zustw(1:nproma,jrow)                    )
    
    CALL precalc_ice( nproma                                           &
                    , surface%is_seaice(1:nproma,jrow)                 &
                    , ice%surface_temperature(1:nproma,jrow)           &
                    , atm_half_lev_press(:,levels_plus_1)              &
                    , atm_spec_humidity(:)                             &
                    , box%atm_tot_cloud_water(1:nproma,jrow)           &
                    , atm_temp(:)                                      &
                    , box%atm_sat_spec_hum(1:nproma,jrow)              &
                    , box%atm_pot_temp(1:nproma,jrow)                  &
                    , box%atm_vir_pot_temp(1:nproma,jrow)              &
                    , box%atm_lat_heat_fact(1:nproma,jrow)             &
                    , cloud_cover(:)                                   &
                    , box%atm_liq_wat_pot_temp(1:nproma,jrow)          &
                    , geopotential(:)                                  &
                    , wind_u(:)                                        &
                    , wind_v(:)                                        &
                    , ocean_u(:)                                       &
                    , ocean_v(:)                                       &
                    , ice%roughness(1:nproma,jrow)                     &
                    , zghabl(:)                                        &
                    , ice%zqsi(1:nproma,jrow)                          &
                    , ice%zcpti(1:nproma,jrow)                         &        
                    , ice%zrii(1:nproma,jrow)                          &
                    , ice%zcfmi(1:nproma,jrow)                         & 
                    , ice%zchi(1:nproma,jrow)                          &
                    , ice%zcfhi(1:nproma,jrow)                         & 
                    , ice%zbni(1:nproma,jrow)                          &
                    , ice%zbmi(1:nproma,jrow)                          & 
                    , ice%zbhi(1:nproma,jrow)                          &
                    , ice%zustari(1:nproma,jrow)                       & 
                    , ice%ztkevi(1:nproma,jrow)                        &
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
                     , zcfnci(1:nproma)                                &
                     , ztvi(1:nproma)                                  &
                     , zcdni(1:nproma)                                 &
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------
                    , ice%zusti(1:nproma,jrow)                        )
    !-------------------------------------------------------------------
    ! CALCULATE NEW ZO VALUES

    CALL update_z0_ocean( nproma                                       &
                        , levels_plus_1                                &
                        , surface%is_ocean(1:nproma,jrow)              &
                        , ocean%zcfmw(1:nproma,jrow)                   &
                        , zudif(1:nproma)                              &
                        , zvdif(1:nproma)                              &
                        , atm_temp(1:nproma)                           &
                        , atm_spec_humidity(1:nproma)                  &
                        , box%atm_tot_cloud_water(1:nproma,jrow)       &
                        , atm_half_lev_press(1:nproma,:)               &
                        , ocean%roughness(1:nproma,jrow)              )
    
    CALL update_z0_ice( nproma                                         &
                      , surface%is_seaice(1:nproma,jrow)               &
                      , ice%roughness(1:nproma,jrow)                  )
    !-------------------------------------------------------------------
    !   RICHTMYER-MORTON COEFFICIENTS FOR DIFFERENT SURFACE FRACTIONS
    !
    CALL richtmyer_land( nproma                                        &
                       , surface%is_land(1:nproma,jrow)                &
                       , levels, levels_plus_1                         &
                       , levels_minus_1                                &
                       , atm_half_lev_press(:,:)                       &
                       , zcfh(:,:)                                     &
                       , zebsh(:,:)                                    &
                       , zqdif(:,:)                                    & 
                       , ztdif(:,:)                                    &
                       , land%zcfhl(1:nproma,jrow)                     &
                       , land%zcair(1:nproma,jrow)                     &
                       , land%zcsat(1:nproma, jrow)                    &
                       , land%zetnl(1:nproma, jrow)                    &
                       , land%zftnl(1:nproma, jrow)                    &
                       , land%zeqnl(1:nproma, jrow)                    &
                       , land%zfqnl(1:nproma, jrow)                   )
    
    CALL richtmyer_ocean( nproma                                       &
                        , surface%is_ocean(1:nproma,jrow)              &
                        , levels, levels_plus_1                        &
                        , levels_minus_1                               &
                        , atm_half_lev_press(:,:)                      &
                        , zcfh(:,:)                                    &
                        , zebsh(:,:)                                   &
                        , zqdif(:,:)                                   &
                        , ztdif(:,:)                                   &
                        , ocean%zcfhw(1:nproma, jrow)                  &
                        , ocean%zetnw(1:nproma, jrow)                  &
                        , ocean%zftnw(1:nproma, jrow)                  &
                        , ocean%zeqnw(1:nproma, jrow)                  &
                        , ocean%zfqnw(1:nproma, jrow)                 )
    
    CALL richtmyer_ice( nproma                                         &
                      , surface%is_seaice(1:nproma,jrow)               &
                      , levels, levels_plus_1                          &
                      , levels_minus_1                                 &
                      , atm_half_lev_press(:,:)                        &
                      , zcfh(:,:)                                      &
                      , zebsh(:,:)                                     &
                      , zqdif(:,:)                                     &
                      , ztdif(:,:)                                     &
                      , ice%zcfhi(1:nproma, jrow)                      &
                      , ice%zetni(1:nproma, jrow)                      &
                      , ice%zftni(1:nproma, jrow)                      &
                      , ice%zeqni(1:nproma, jrow)                      &
                      , ice%zfqni(1:nproma, jrow)                     )
    !-------------------------------------------------------------------
    ! CALCULATE NEW ALBEDO VALUES

    sw_vis_land = 0._dp
    sw_nir_land = 0._dp
    WHERE (surface%is_land(1:nproma, jrow))
       sw_vis_land(:) = sw_vis(:) * (1._dp - land%albedo_vis(1:nproma,jrow)) / &
                        (1._dp - box%albedo_vis(1:nproma,jrow))
       sw_nir_land(:) = sw_nir(:) * (1._dp - land%albedo_nir(1:nproma,jrow)) / &
                        (1._dp - box%albedo_nir(1:nproma,jrow))
    END WHERE

    CALL update_albedo_land( nproma                                    & ! in
                    , surface%is_land(1:nproma, jrow)                  & ! in
                    , land%albedo_vis(1:nproma, jrow)                  & ! in
                    , land%albedo_nir(1:nproma, jrow)                  & ! in
                    , sw_vis_land(:)                                   & ! in
                    , sw_nir_land(:)                                   & ! in
                    , land%albedo(1:nproma, jrow)                     )  ! out

      CALL update_albedo_ocean(nproma                                    &
                        , surface%is_ocean(1:nproma,jrow)                &
                        , cos_zenith_angle(:)                            &
                        , sw_vis(:)                                      &
                        , sw_vis_frac_diffuse(:)                         &
                        , sw_nir(:)                                      &
                        , sw_nir_frac_diffuse(:)                         &
                        , palbedo_vis_dir(:)                             &
                        , palbedo_vis_dif(:)                             &
                        , palbedo_nir_dir(:)                             &
                        , palbedo_nir_dif(:)                             &
                        , ocean%albedo_vis_dir(1:nproma,jrow)            &
                        , ocean%albedo_nir_dir(1:nproma,jrow)            &
                        , ocean%albedo_vis_dif(1:nproma,jrow)            &
                        , ocean%albedo_nir_dif(1:nproma,jrow)            &
                        , ocean%albedo_vis(1:nproma,jrow)                &
                        , ocean%albedo_nir(1:nproma,jrow)                &
                        , ocean%albedo(1:nproma,jrow)                   )
    IF(.NOT.lmeltpond) THEN
      CALL update_albedo_ice(surface%is_seaice(1:nproma,jrow)            &
                            , ice%surface_temperature(1:nproma,jrow)     &
                            , ice%snow_water_equivalent(1:nproma,jrow)   &
                            , ice%albedo(1:nproma,jrow)                 )
      ice%albedo_vis(1:nproma,jrow)       = ice%albedo(1:nproma,jrow)
      ice%albedo_vis_dir(1:nproma,jrow)   = ice%albedo(1:nproma,jrow)
      ice%albedo_vis_dif(1:nproma,jrow)   = ice%albedo(1:nproma,jrow)
      ice%albedo_nir(1:nproma,jrow)       = ice%albedo(1:nproma,jrow)
      ice%albedo_nir_dir(1:nproma,jrow)   = ice%albedo(1:nproma,jrow)
      ice%albedo_nir_dif(1:nproma,jrow)   = ice%albedo(1:nproma,jrow)
   ! 
   ! bypass for meltponds
   !
      pfage(1:nproma)      = 0.0_wp
      psnifrac(1:nproma)   = 0.0_wp
      pbarefrac(1:nproma)  = 0.0_wp

    ELSE
      CALL update_albedo_ice_meltpond(nproma                             &
                        , surface%is_seaice(1:nproma,jrow)               &
                        , ice%surface_temperature(1:nproma,jrow)         &
                        , snow(:)                                        &
                        , ice%zcvsi(1:nproma,jrow)                       &
                        , ice%ice_depth(1:nproma,jrow)                   &
                        , cos_zenith_angle(:)                            &
                        , pameltdepth(1:nproma)                          &
                        , pameltfrac(1:nproma)                           &
                        , ptaus(1:nproma)                                &
                        , pfage(1:nproma)                                &
                        , psnifrac(1:nproma)                             &
                        , pbarefrac(1:nproma)                            &
                        , sw_vis(:)                                      &
                        , sw_vis_frac_diffuse(:)                         &
                        , sw_nir(:)                                      &
                        , sw_nir_frac_diffuse(:)                         &
                        , palbedo_vis_dir(:)                             &
                        , palbedo_vis_dif(:)                             &
                        , palbedo_nir_dir(:)                             &
                        , palbedo_nir_dif(:)                             &
                        , ice%albedo_vis_dir(1:nproma,jrow)              &
                        , ice%albedo_nir_dir(1:nproma,jrow)              &
                        , ice%albedo_vis_dif(1:nproma,jrow)              &
                        , ice%albedo_nir_dif(1:nproma,jrow)              &
                        , ice%albedo_vis(1:nproma,jrow)                  &
                        , ice%albedo_nir(1:nproma,jrow)                  &
                        , ice%alsom(1:nproma, jrow)                      &
                        , ice%alsobs(1:nproma, jrow)                     &
                        , ice%albedo(1:nproma,jrow)                   )
    END IF
    !------------------------------------------------------------------- ! Calculate windstress
    CALL update_stress_land( nproma                                    &
                           , surface%is_land(1:nproma,jrow)            &
                           , land%zcfml(1:nproma,jrow)                 &
                           , zudif(:)                                  &
                           , zvdif(:)                                  &
                           , land%u_stress(1:nproma,jrow)              &
                           , land%v_stress(1:nproma,jrow)             )
    CALL update_stress_ocean( nproma                                   &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , ocean%zcfmw(1:nproma,jrow)               &
                            , zudif(:)                                 &
                            , zvdif(:)                                 &
                            , ocean%u_stress(1:nproma,jrow)            &
                            , ocean%v_stress(1:nproma,jrow)           )
    CALL update_stress_ice( nproma                                     &
                          , surface%is_seaice(1:nproma,jrow)           &
                          , ice%zcfmi(1:nproma,jrow)                   &
                          , zudif(:)                                   &
                          , zvdif(:)                                   &
                          , ice%u_stress(1:nproma,jrow)                &
                          , ice%v_stress(1:nproma,jrow)               )
    !-------------------------------------------------------------------
    ! Derive snow cover fraction from snow water content
    !
    ice%snow_cover_fract(1:nproma, jrow) =                             &
                TANH(ice%snow_water_equivalent(1:nproma, jrow) * 100._wp)    
    !-------------------------------------------------------------------
    ! PRECALC FOR SURFACE LAYER ELIMINATION
    
    CALL update_ocean( nproma                                          &
                     , surface%is_ocean(1:nproma,jrow)                 &
                     , ocean%zetnw(1:nproma, jrow)                     &
                     , ocean%zftnw(1:nproma, jrow)                     &
                     , ocean%zeqnw(1:nproma, jrow)                     &
                     , ocean%zfqnw(1:nproma, jrow)                     &
                     , ocean%zcptw(1:nproma,jrow)                      &
                     , ocean%zqsw(1:nproma,jrow)                       &
                     , ocean%ztklevw(1:nproma,jrow)                    &
                     , ocean%zqklevw(1:nproma,jrow)                   )
    
    CALL update_ice( nproma                                            &
                   , surface%is_seaice(1:nproma,jrow)                  &
                   , ice%zetni(1:nproma, jrow)                         &
                   , ice%zftni(1:nproma, jrow)                         &
                   , ice%zeqni(1:nproma, jrow)                         &
                   , ice%zfqni(1:nproma, jrow)                         &
                   , ice%zcpti(1:nproma,jrow)                          &
                   , ice%zqsi(1:nproma,jrow)                           &
                   , ice%ztklevi(1:nproma,jrow)                        &
                   , ice%zqklevi(1:nproma,jrow)                       )

    land%zcair_old(1:nproma,jrow) = land%zcair(1:nproma,jrow)

    CALL jsbach_inter_1d( nproma,                                      &
                          COUNT(surface%is_land(1:nproma,jrow)),       &
         wind           = SQRT(wind_u(1:nproma)**2 + wind_v(1:nproma)**2),     &
         wind10         = land%wind_10_meter(1:nproma,jrow),           & 
         ! NOTE: This is the old 10m wind speed, we don't
         ! have the new one yet (as is the case in ECHAM)
         temp_air       = atm_temp(1:nproma),                          &
         qair           = atm_spec_humidity(1:nproma),                 &
         precip_rain    = rain,                                        &
         precip_snow    = snow,                                        &
         lwdown         = longwave_down,                               &
         sw_vis_net     = sw_vis_land(:),                              &
         sw_nir_net     = sw_nir_land(:),                              &
         sw_par_down    = sw_par(:),                                   &
         sw_par_frac_diffuse = sw_par_frac_diffuse(:),                 &
         pressure       = atm_half_lev_press(1:nproma,levels_plus_1),  &
         czenith        = cos_zenith_angle(1:nproma),                  &
         CO2_concentration = pco2_concentration(1:nproma),             &
         cdrag          = land%zcfhl(1:nproma,jrow),                   &
         etAcoef        = land%zetnl(1:nproma,jrow),                   &
         etBcoef        = land%zftnl(1:nproma,jrow),                   &
         eqAcoef        = land%zeqnl(1:nproma,jrow),                   &
         eqBcoef        = land%zfqnl(1:nproma,jrow),                   &
         cair           = land%zcair(1:nproma,jrow),                   &
         csat           = land%zcsat(1:nproma,jrow),                   &
         albedo_vis     = land%albedo_vis(1:nproma,jrow),              &
         albedo_nir     = land%albedo_nir(1:nproma,jrow),              &
         z0h            = land%roughness_heat(1:nproma,jrow),          &
         z0m            = land%roughness_momentum(1:nproma,jrow),      &
         evap_act       = land%evaporation_inst(1:nproma,jrow),        &
         evap_pot       = land%evaporation_pot(1:nproma,jrow),         &
         tsoil_rad      = land%surface_temperature_rad(1:nproma,jrow), &
         temp_soil_new  = land%surface_temperature_new(1:nproma,jrow), &
         qsurf          = land%surface_qsat_new(1:nproma,jrow),        &
         mask_land      = surface%is_land(1:nproma,jrow),              &
         land_fract     = box%land_fract(1:nproma,jrow),               &
         kblock         = jrow,                                        &
         latent         = land%latent_heat_flux_inst(1:nproma,jrow),   &
         sensible       = land%sensible_heat_flux_inst(1:nproma,jrow), &
         echam_zchl     = land%zchl(1:nproma,jrow),                    &
         surf_dry_static_energy = land%dry_static_energy_new(1:nproma,jrow), &
         soil_wetness   = soil_wetness(1:nproma),                      &
         snow_depth     = snow_depth(1:nproma),                        &
         skin_res       = skin_reservoir(1:nproma),                    & 
         tte_corr       = tte_corr(1:nproma),                          &
         glacier_depth = gld(1:nproma,jrow),                           &
         snow_melt     = zsnow_melt(1:nproma),                         &
         wsmax = wsmx(1:nproma,jrow),                                  &
         CO2_flux_npp         = pco2_flux_npp(1:nproma),               &
         CO2_flux_soilresp    = pco2_flux_soilresp(1:nproma),          &
         CO2_flux_herbivory   = pco2_flux_herbivory(1:nproma),         &
         CO2_flux_dynveg      = pco2_flux_dynveg(1:nproma),            &
         CO2_emission_lcc     = pco2_emission_lcc(1:nproma),           &
         CO2_emission_harvest = pco2_emission_harvest(1:nproma))

    IF (lanysubmodel) THEN
      IF (ASSOCIATED(vphysc%smelt)) THEN
        vphysc%smelt(1:nproma,jrow)=zsnow_melt(1:nproma)
      ENDIF
    ENDIF

    pco2_flux_land(1:nproma) = pco2_flux_npp(1:nproma) +               &
                               pco2_flux_soilresp(1:nproma) +          &
                               pco2_flux_herbivory(1:nproma) +         &
                               pco2_flux_dynveg(1:nproma)

    CALL update_land( nproma                                           &
                    , surface%is_land(1:nproma, jrow)                  &
                    , land%zetnl(1:nproma, jrow)                       &
                    , land%zftnl(1:nproma, jrow)                       &
                    , land%zeqnl(1:nproma, jrow)                       &
                    , land%zfqnl(1:nproma, jrow)                       &
                    , land%dry_static_energy_new(1:nproma,jrow)        &
                    , land%surface_qsat_new(1:nproma,jrow)             &
                    , land%ztklevl(1:nproma,jrow)                      &
                    , land%zqklevl(1:nproma,jrow)                     )

    !-------------------------------------------------------------------
    ! Copy land albedo to g3b memory stream
!    alsol(1:nproma,jrow) = land%albedo(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! update ztdif and zq dif due to new surface values for humidity and 
    ! temperature at the blending height (means lowest atmosphere level)

    CALL blend_zq_zt( nproma                                           &
                    , box%land_fract(1:nproma,jrow)                    &
                    , box%ocean_fract(1:nproma,jrow)                   &
                    , box%seaice_fract(1:nproma,jrow)                  &
                    , land%ztklevl(1:nproma,jrow)                      &
                    , ocean%ztklevw(1:nproma,jrow)                     &
                    , ice%ztklevi(1:nproma,jrow)                       &
                    , land%zqklevl(1:nproma,jrow)                      &
                    , ocean%zqklevw(1:nproma,jrow)                     &
                    , ice%zqklevi(1:nproma,jrow)                       &
                    , ztdif(:,levels)                                  &
                    , zqdif(:,levels)                                  &
                    , surface%is_land(1:nproma,jrow)                   &
                    , surface%is_ocean(1:nproma,jrow)                  &
                    , surface%is_seaice(1:nproma,jrow)                )
    
    ztdif_new(1:nproma) = ztdif(1:nproma,levels)
    zqdif_new(1:nproma) = zqdif(1:nproma,levels)
    !
    !-------------------------------------------------------------------
    CALL postproc_ocean( nproma                                        &
                       , levels_plus_1                                 &
                       , surface%is_ocean(1:nproma,jrow)               &
                       , ocean%zcfhw(1:nproma, jrow)                   &
                       , ocean%zqsw(1:nproma, jrow)                    &
                       , ocean%zqklevw(1:nproma,jrow)                  &
                       , atm_spec_humidity(:)                          &
                       , geopotential(:)                               &
                       , ocean%ztklevw(1:nproma,jrow)                  &
                       , box%atm_dry_stat_energy(1:nproma,jrow)        &
                       , ocean%zcptw(1:nproma, jrow)                   &
                       , ocean%surface_temperature(1:nproma,jrow)      &
                       , ocean%zbnw(1:nproma, jrow)                    &
                       , ocean%zbhw(1:nproma, jrow)                    &
                       , ocean%zriw(1:nproma, jrow)                    &
                       , atm_temp(:)                                   &
                       , atm_full_lev_press(:)                         &
                       , box%atm_tot_cloud_water(1:nproma, jrow)       &
                       , wind_u(:)                                     &
                       , wind_v(:)                                     &
                       , ocean_u(:)                                    &
                       , ocean_v(:)                                    &
                       , ocean%evaporation_inst(1:nproma, jrow)        &
                       , ocean%sensible_heat_flux_inst(1:nproma, jrow) &
                       , ocean%latent_heat_flux_inst(1:nproma, jrow)   &
                       , ocean%wind_10_meter(1:nproma, jrow)           &
                       , ocean%temp_2_meter(1:nproma, jrow)            &
                       , atm_half_lev_press(:,:)                       &
                       , ocean%zbmw(1:nproma,jrow)                     &
                       , ocean%dewpoint_2_meter(1:nproma,jrow)         &
                       , ocean%u_wind_10_meter(1:nproma,jrow)          &
                       , ocean%v_wind_10_meter(1:nproma,jrow)         )

    CALL postproc_ice( nproma                                          &
                     , levels_plus_1                                   &
                     , surface%is_seaice(1:nproma,jrow)                &
                     , ice%zcfhi(1:nproma, jrow)                       &
                     , ice%zqsi(1:nproma, jrow)                        &
                     , ice%zqklevi(1:nproma,jrow)                      &
                     , atm_spec_humidity(:)                            &
                     , geopotential(:)                                 &
                     , ice%ztklevi(1:nproma,jrow)                      &
                     , box%atm_dry_stat_energy(1:nproma,jrow)          &
                     , ice%zcpti (1:nproma, jrow)                      &
                     , ice%surface_temperature(1:nproma,jrow)          &
                     , ice%zbni(1:nproma, jrow)                        &
                     , ice%zbhi (1:nproma, jrow)                       &
                     , ice%zrii(1:nproma, jrow)                        &
                     , atm_temp(:)                                     &
                     , atm_full_lev_press(:)                           &
                     , box%atm_tot_cloud_water(1:nproma, jrow)         &
                     , wind_u(:)                                       &
                     , wind_v(:)                                       &
                     , ice%evaporation_inst(1:nproma, jrow)            &
                     , ice%sensible_heat_flux_inst(1:nproma, jrow)     &
                     , ice%latent_heat_flux_inst (1:nproma, jrow)      &
                     , ice%wind_10_meter(1:nproma, jrow)               &
                     , ice%temp_2_meter(1:nproma, jrow)                &
                     , atm_half_lev_press(:,:)                         &
                     , ice%zbmi(1:nproma, jrow)                        &
                     , ice%dewpoint_2_meter(1:nproma, jrow)            &
                     , ice%u_wind_10_meter(1:nproma,jrow)              &
                     , ice%v_wind_10_meter(1:nproma,jrow)             )

    CALL postproc_land( nproma                                         &
                      , surface%is_land(1:nproma,jrow)                 &
                      , geopotential(:)                                &
                      , land%zbnl(1:nproma,jrow)                       &
                      , land%zbml(1:nproma,jrow)                       &
                      , wind_u(:)                                      &
                      , wind_v(:)                                      &
                      , land%ZRIL(1:nproma,jrow)                       & 
                      , land%wind_10_meter(1:nproma,jrow)              &
                      , land%zbhnl(1:nproma,jrow)                      &
                      , land%zbhl(1:nproma,jrow)                       &
                      , land%zcptl(1:nproma,jrow)                      &
                      , box%atm_dry_stat_energy(1:nproma,jrow)         &
                      , atm_spec_humidity(:)                           &
                      , land%temp_2_meter(1:nproma,jrow)               &
                      , atm_temp(:)                                    &
                      , atm_full_lev_press(:)                          &
                      , atm_half_lev_press(:,:)                        &
                      , levels_plus_1                                  &
                      , box%atm_tot_cloud_water(1:nproma, jrow)        &
                      , land%dewpoint_2_meter(1:nproma, jrow)          &
                      , land%u_wind_10_meter(1:nproma,jrow)            &
                      , land%v_wind_10_meter(1:nproma,jrow)           )         
    !-------------------------------------------------------------------
    ! Average albedo
       CALL surface_box_average( nproma                                &
                               , land%albedo (1:nproma,jrow)           &
                               , ice%albedo  (1:nproma,jrow)           &
                               , ocean%albedo(1:nproma,jrow)           &
                               , box%albedo(1:nproma,jrow)             &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
       CALL surface_box_average( nproma                                &
                               , land%albedo_vis(1:nproma,jrow)        &
                               , ice%albedo_vis(1:nproma,jrow)         &
                               , ocean%albedo_vis(1:nproma,jrow)       &
                               , box%albedo_vis(1:nproma,jrow)         &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
       CALL surface_box_average( nproma                                &
                               , land%albedo_vis(1:nproma,jrow)        &
                               , ice%albedo_vis_dir(1:nproma,jrow)     &
                               , ocean%albedo_vis_dir(1:nproma,jrow)   &
                               , box%albedo_vis_dir(1:nproma,jrow)     &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
       CALL surface_box_average( nproma                                &
                               , land%albedo_vis(1:nproma,jrow)        &
                               , ice%albedo_vis_dif(1:nproma,jrow)     &
                               , ocean%albedo_vis_dif(1:nproma,jrow)   &
                               , box%albedo_vis_dif(1:nproma,jrow)     &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
       CALL surface_box_average( nproma                                &
                               , land%albedo_nir(1:nproma,jrow)        &
                               , ice%albedo_nir(1:nproma,jrow)         &
                               , ocean%albedo_nir(1:nproma,jrow)       &
                               , box%albedo_nir(1:nproma,jrow)         &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
       CALL surface_box_average( nproma                                &
                               , land%albedo_nir(1:nproma,jrow)        &
                               , ice%albedo_nir_dir(1:nproma,jrow)     &
                               , ocean%albedo_nir_dir(1:nproma,jrow)   &
                               , box%albedo_nir_dir(1:nproma,jrow)     &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
       CALL surface_box_average( nproma                                &
                               , land%albedo_nir(1:nproma,jrow)        &
                               , ice%albedo_nir_dif(1:nproma,jrow)     &
                               , ocean%albedo_nir_dif(1:nproma,jrow)   &
                               , box%albedo_nir_dif(1:nproma,jrow)     &
                               , surface%is_land(1:nproma,jrow)        &
                               , surface%is_ocean(1:nproma,jrow)       &
                               , surface%is_seaice(1:nproma,jrow)      &
                               , box%land_fract(1:nproma,jrow)         &
                               , box%ocean_fract(1:nproma,jrow)        &
                               , box%seaice_fract(1:nproma,jrow)      )
    !-------------------------------------------------------------------
    ! calculate surface lw and sw radiation for ice and lakeice

    CALL ice_rad( nproma                                               &
                , surface%is_seaice(1:nproma,jrow)                     &
                , longwave_down(:)                                     &
                , ice%surface_temperature(1:nproma,jrow)               &
                , pi0(:)                                               &
                , ptrsol(:)                                            &
                , ice%albedo  (1:nproma,jrow)                          &
                , box%albedo(1:nproma,jrow)                            &
                , ice%sofli(1:nproma,jrow)                             &
                , ice%trfli(1:nproma,jrow)                            )
    CALL ocean_rad( nproma                                             &
                  , surface%is_ocean(1:nproma,jrow)                    &
                  , longwave_down(:)                                   &
                  , ocean%surface_temperature(1:nproma,jrow)           &
                  , pi0(:)                                             &
                  , ptrsol(:)                                          &
                  , ocean%albedo  (1:nproma,jrow)                      &
                  , box%albedo(1:nproma,jrow)                          &
                  , ocean%soflw(1:nproma,jrow)                         &
                  , ocean%trflw(1:nproma,jrow)                        )
    CALL land_rad( nproma                                              &
                 , surface%is_land(1:nproma,jrow)                      &
                 , longwave_down(:)                                    &
                 , land%surface_temperature(1:nproma,jrow)             &
                 , pi0(:)                                              &
                 , ptrsol(:)                                           &
                 , land%albedo  (1:nproma,jrow)                        &
                 , box%albedo(1:nproma,jrow)                           &
                 , land%sofll(1:nproma,jrow)                           &
                 , land%trfll(1:nproma,jrow)                           &
                 , land%surface_temperature_rad(1:nproma,jrow)         &
                 , pzteffl4(1:nproma)                                 )
    land%sofllac(1:nproma,jrow)  =  land%sofllac(1:nproma,jrow)        &
             + box%land_fract(1:nproma,jrow) * delta_time              &
             * land%sofll(1:nproma,jrow)
    land%trfllac(1:nproma,jrow)  =  land%trfllac(1:nproma,jrow)        &
             + box%land_fract(1:nproma,jrow) * delta_time              &
             * land%trfll(1:nproma,jrow)
    ocean%soflwac(1:nproma,jrow) =  ocean%soflwac(1:nproma,jrow)       &
             + box%ocean_fract(1:nproma,jrow) * delta_time             &
             * ocean%soflw(1:nproma,jrow)
    ocean%trflwac(1:nproma,jrow) =  ocean%trflwac(1:nproma,jrow)       &
             + box%ocean_fract(1:nproma,jrow) * delta_time             &
             * ocean%trflw(1:nproma,jrow)
    ice%sofliac(1:nproma,jrow)   =  ice%sofliac(1:nproma,jrow)         &
             + box%seaice_fract(1:nproma,jrow) * delta_time            &
             * ice%sofli(1:nproma,jrow)
    ice%trfliac(1:nproma,jrow)   =  ice%trfliac(1:nproma,jrow)         &
             + box%seaice_fract(1:nproma,jrow) * delta_time            &
             * ice%trfli(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average roughness length
    CALL surface_box_average( nproma                                   &
                            , land%roughness_momentum(1:nproma,jrow)   &
                            , ice%roughness(1:nproma,jrow)             &
                            , ocean%roughness(1:nproma,jrow)           &
                            , box%roughness(1:nproma,jrow)             &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )

    CALL surface_box_average(nproma                                    &
                            , zcfncl(:)                                &
                            , zcfnci(:)                                &
                            , zcfncw(:)                                &
                            , box%ZCFNC(1:nproma,jrow)                 &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )

    CALL surface_box_average(nproma                                    &
                            , zcdnl(:)                                 &
                            , zcdni(:)                                 &
                            , zcdnw(:)                                 &
                            , zcdn(:)                                  &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )

    CALL surface_box_average(nproma                                    &
                            , land%ZRIL(1:nproma,jrow)                 &
                            , ice%zrii(1:nproma,jrow)                  &
                            , ocean%zriw(1:nproma,jrow)                &
                            , zri(:)                                   &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )

    !-------------------------------------------------------------------
    ! Average momentum exchange coefficient
    CALL surface_box_average( nproma                                   &
                            , land%zcfml(1:nproma,jrow)                &
                            , ice%zcfmi(1:nproma,jrow)                 &
                            , ocean%zcfmw(1:nproma,jrow)               &
                            , box%momentum_ex_coef(1:nproma,jrow)      &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         ) 
    momentum_exchange_coeff(1:nproma)    =                             &
                                   box%momentum_ex_coef(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average heat exchange coefficient
    CALL surface_box_average( nproma                                   &
                            , land%zcfhl(1:nproma,jrow)                &
                            , ice%zcfhi(1:nproma,jrow)                 &
                            , ocean%zcfhw(1:nproma,jrow)               &
                            , box%heat_ex_coef(1:nproma,jrow)          &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         ) 
        heat_exchange_coeff(1:nproma)    =                             &
                                   box%heat_ex_coef(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average ustar
    CALL average_ustar( nproma                                         &
                      , levels_plus_1                                  & !! in
                      , box%land_fract(1:nproma,jrow)                  &
                      , box%ocean_fract(1:nproma,jrow)                 & !! in
                      , box%seaice_fract(1:nproma,jrow)                &
                      , land%zustl(1:nproma,jrow)                      & !! in
                      , ocean%zustw(1:nproma,jrow)                     &
                      , ice%zusti(1:nproma,jrow)                       & !! in
                      , atm_temp(:)                                    &
                      , atm_spec_humidity(:)                           & !! in
                      , box%atm_tot_cloud_water(1:nproma,jrow)         &
                      , atm_half_lev_press(:,:)                        & !! in
                      , box%ustar(1:nproma,jrow)                       &
                      , surface%is_land(1:nproma,jrow)                 & !! out, in
                      , surface%is_ocean(1:nproma,jrow)                &
                      , surface%is_seaice(1:nproma,jrow)              )
    ustarm(1:nproma)                   = box%ustar(1:nproma,jrow) 
    !-------------------------------------------------------------------
    ! Average TKE condition
    CALL surface_box_average( nproma                                   &
                            , land%ztkevl(1:nproma,jrow)               &
                            , ice%ztkevi(1:nproma,jrow)                &
                            , ocean%ztkevw(1:nproma,jrow)              &
                            , box%ztkevn(1:nproma,jrow)                &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    tkevn_cond(:)              = box%ztkevn(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average wind stress
    CALL surface_box_average( nproma                                   &
                            , land%u_stress(1:nproma,jrow)             &
                            , ice%u_stress(1:nproma,jrow)              &
                            , ocean%u_stress(1:nproma,jrow)            &
                            , box%u_stress(1:nproma,jrow)              &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    CALL surface_box_average( nproma                                   &
                            , land%v_stress(1:nproma,jrow)             &
                            , ice%v_stress(1:nproma,jrow)              &
                            , ocean%v_stress(1:nproma,jrow)            &
                            , box%v_stress(1:nproma,jrow)              &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    box%u_stress_acc(1:nproma,jrow) =  box%u_stress_acc(1:nproma,jrow) &
                          + box%u_stress(1:nproma,jrow) * delta_time
    ustr_na(1:nproma,jrow) = box%u_stress(1:nproma,jrow)
    box%v_stress_acc(1:nproma,jrow) =  box%v_stress_acc(1:nproma,jrow) &
                          + box%v_stress(1:nproma,jrow) * delta_time
    vstr_na(1:nproma,jrow) = box%v_stress(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average 10 meter wind speed
    CALL surface_box_average( nproma                                   &
                            , land%wind_10_meter(1:nproma,jrow)        &
                            , ice%wind_10_meter(1:nproma,jrow)         &
                            , ocean%wind_10_meter(1:nproma,jrow)       &
                            , box%wind_10_meter(1:nproma,jrow)         &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    CALL surface_box_average( nproma                                   &
                            , land%u_wind_10_meter(1:nproma,jrow)      &
                            , ice%u_wind_10_meter(1:nproma,jrow)       &
                            , ocean%u_wind_10_meter(1:nproma,jrow)     &
                            , box%u_wind_10_meter(1:nproma,jrow)       &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )
    CALL surface_box_average( nproma                                   &
                            , land%v_wind_10_meter(1:nproma,jrow)      &
                            , ice%v_wind_10_meter(1:nproma,jrow)       &
                            , ocean%v_wind_10_meter(1:nproma,jrow)     &
                            , box%v_wind_10_meter(1:nproma,jrow)       &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )

    box%wind_10_meter_acc(1:nproma,jrow)    =                          &
                                  box%wind_10_meter_acc(1:nproma,jrow) &
                       + box%wind_10_meter(1:nproma,jrow) *delta_time
    wind10_na(1:nproma,jrow) = box%wind_10_meter(1:nproma,jrow)
    box%wind_10_max(1:nproma,jrow)          =                          &
                                  MAX(box%wind_10_max(1:nproma,jrow) , &
                        box%wind_10_meter(1:nproma,jrow) )
    !-------------------------------------------------------------------
    ! Average moisture flux (evaporation)
    CALL surface_box_average( nproma                                   &
                            , land%evaporation_inst(1:nproma,jrow)     &
                            , ice%evaporation_inst(1:nproma,jrow)      &
                            , ocean%evaporation_inst(1:nproma,jrow)    &
                            , box%evaporation_inst(1:nproma,jrow)      &
                            , surface%is_land(1:nproma,jrow)           &
                            , surface%is_ocean(1:nproma,jrow)          &
                            , surface%is_seaice(1:nproma,jrow)         &
                            , box%land_fract(1:nproma,jrow)            &
                            , box%ocean_fract(1:nproma,jrow)           &
                            , box%seaice_fract(1:nproma,jrow)         )

    ice%evaporation_acc(1:nproma,jrow) =                               &
                                 ice%evaporation_acc(1:nproma,jrow)    &
                       + MERGE((ice%evaporation_inst(1:nproma,jrow)    &
                       * box%seaice_fract(1:nproma,jrow)), zero,       &
                       surface%is_seaice(1:nproma,jrow)) * delta_time
    ocean%evaporation_acc(1:nproma,jrow) =                             &
                               ocean%evaporation_acc(1:nproma,jrow)    &
                     + MERGE((ocean%evaporation_inst(1:nproma,jrow)    &
                     * box%ocean_fract(1:nproma,jrow)), zero,          &
                     surface%is_ocean(1:nproma,jrow)) * delta_time
    land%evaporation_acc(1:nproma,jrow) =                              &
                               land%evaporation_acc(1:nproma,jrow)     &
                     + MERGE((land%evaporation_inst(1:nproma,jrow)     &
                     * box%land_fract(1:nproma,jrow)), zero,           &
                     surface%is_land(1:nproma,jrow)) *delta_time         
    box%evaporation_acc(1:nproma,jrow) =                               &
                               box%evaporation_acc(1:nproma,jrow)      &
                   + box%evaporation_inst(1:nproma,jrow) * delta_time
    evap_na(1:nproma,jrow) = box%evaporation_inst(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average sensible heat flux
    CALL surface_box_average( nproma                                   &
                      , land%sensible_heat_flux_inst(1:nproma,jrow)    &
                      , ice%sensible_heat_flux_inst(1:nproma,jrow)     &
                      , ocean%sensible_heat_flux_inst(1:nproma,jrow)   &
                      , box%sensible_heat_flux_inst(1:nproma,jrow)     &
                      , surface%is_land(1:nproma,jrow)                 &
                      , surface%is_ocean(1:nproma,jrow)                &
                      , surface%is_seaice(1:nproma,jrow)               &
                      , box%land_fract(1:nproma,jrow)                  &
                      , box%ocean_fract(1:nproma,jrow)                 &
                      , box%seaice_fract(1:nproma,jrow)               )

    land%sensible_flux_acc(1:nproma,jrow) =                            &
                                 land%sensible_flux_acc(1:nproma,jrow) & 
                + MERGE(land%sensible_heat_flux_inst(1:nproma,jrow),   &
                zero, surface%is_land(1:nproma,jrow))                  &
                * box%land_fract(1:nproma,jrow) * delta_time
    ocean%sensible_flux_acc(1:nproma,jrow) =                           &
                                ocean%sensible_flux_acc(1:nproma,jrow) & 
                + MERGE(ocean%sensible_heat_flux_inst(1:nproma,jrow),  &
                zero, surface%is_ocean(1:nproma,jrow))                 &
                * box%ocean_fract(1:nproma,jrow) * delta_time
    ice%sensible_flux_acc(1:nproma,jrow) =                             &
                                  ice%sensible_flux_acc(1:nproma,jrow) &
                + MERGE(ice%sensible_heat_flux_inst(1:nproma,jrow),    &
                zero,surface%is_seaice(1:nproma,jrow))                 &
                * box%seaice_fract(1:nproma,jrow) * delta_time
    box%sensible_heat_flux_acc(1:nproma,jrow) =                        &
               box%sensible_heat_flux_acc(1:nproma,jrow) +             &
               box%sensible_heat_flux_inst(1:nproma,jrow) * delta_time
    ahfs_na(1:nproma,jrow) = box%sensible_heat_flux_inst(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average latent heat flux
    CALL  surface_box_average( nproma                                  &
                      , land%latent_heat_flux_inst(1:nproma,jrow)      &
                      , ice%latent_heat_flux_inst(1:nproma,jrow)       &
                      , ocean%latent_heat_flux_inst(1:nproma,jrow)     &
                      , box%latent_heat_flux_inst(1:nproma,jrow)       &
                      , surface%is_land(1:nproma,jrow)                 &
                      , surface%is_ocean(1:nproma,jrow)                &
                      , surface%is_seaice(1:nproma,jrow)               &
                      , box%land_fract(1:nproma,jrow)                  &
                      , box%ocean_fract(1:nproma,jrow)                 &
                      , box%seaice_fract(1:nproma,jrow)               )

    ice%latent_heat_flux_acc(1:nproma,jrow)     =                      &
                               ice%latent_heat_flux_acc(1:nproma,jrow) &
                     + MERGE((ice%latent_heat_flux_inst(1:nproma,jrow) &
                     * box%seaice_fract(1:nproma,jrow)), zero,         &
                     surface%is_seaice(1:nproma,jrow)) * delta_time
    ocean%latent_heat_flux_acc(1:nproma,jrow)   =                      &
                             ocean%latent_heat_flux_acc(1:nproma,jrow) &
                   + MERGE((ocean%latent_heat_flux_inst(1:nproma,jrow) &
                   * box%ocean_fract(1:nproma,jrow)), zero,            &
                   surface%is_ocean(1:nproma,jrow)) * delta_time
    land%latent_heat_flux_acc(1:nproma,jrow)    =                      &
                              land%latent_heat_flux_acc(1:nproma,jrow) & 
                    + MERGE((land%latent_heat_flux_inst(1:nproma,jrow) &
                    * box%land_fract(1:nproma,jrow)), zero,            &
                    surface%is_land(1:nproma,jrow)) * delta_time
    box%latent_heat_flux_acc(1:nproma,jrow) =                          &
                             box%latent_heat_flux_acc(1:nproma,jrow) + &
                 box%latent_heat_flux_inst(1:nproma,jrow) * delta_time
    ahfl_na(1:nproma,jrow) = box%latent_heat_flux_inst(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Average 2 meter temperature
    CALL  surface_box_average( nproma                                  &
                             , land%temp_2_meter(1:nproma,jrow)        &
                             , ice%temp_2_meter(1:nproma,jrow)         &
                             , ocean%temp_2_meter(1:nproma,jrow)       &
                             , box%temp_2_meter(1:nproma,jrow)         &
                             , surface%is_land(1:nproma,jrow)          &
                             , surface%is_ocean(1:nproma,jrow)         &
                             , surface%is_seaice(1:nproma,jrow)        &
                             , box%land_fract(1:nproma,jrow)           &
                             , box%ocean_fract(1:nproma,jrow)          &
                             , box%seaice_fract(1:nproma,jrow)        )
    box%maxtemp_2_meter(1:nproma,jrow) =                               &
                           MAX( box%maxtemp_2_meter(1:nproma,jrow),    &
                                box%temp_2_meter(1:nproma,jrow)) 
    box%mintemp_2_meter(1:nproma,jrow) =                               &
                           MIN( box%mintemp_2_meter(1:nproma,jrow),    &
                                box%temp_2_meter(1:nproma,jrow)) 
    !-------------------------------------------------------------------
    ! Average 2 meter dew point
    CALL  surface_box_average( nproma                                  &
                             , land%dewpoint_2_meter(1:nproma,jrow)    &
                             , ice%dewpoint_2_meter(1:nproma,jrow)     &
                             , ocean%dewpoint_2_meter(1:nproma,jrow)   &
                             , box%dewpoint_2_meter(1:nproma,jrow)     &
                             , surface%is_land(1:nproma,jrow)          &
                             , surface%is_ocean(1:nproma,jrow)         &
                             , surface%is_seaice(1:nproma,jrow)        &
                             , box%land_fract(1:nproma,jrow)           &
                             , box%ocean_fract(1:nproma,jrow)          &
                             , box%seaice_fract(1:nproma,jrow)        )
    !-------------------------------------------------------------------
    ! Average virtual surface temperature and surface humidity
    CALL  average_tvh_qsurf( nproma                                    &
                           , land%surface_temperature(1:nproma,jrow)   &
                           , atm_spec_humidity(1:nproma)               &
                           , ocean%surface_temperature(1:nproma,jrow)  &
                           , ice%surface_temperature(1:nproma,jrow)    &
                           , box%land_fract(1:nproma,jrow)             &
                           , box%ocean_fract(1:nproma,jrow)            &
                           , box%seaice_fract(1:nproma,jrow)           &
                           , land%zcsat(1:nproma, jrow)                &
                           , land%zcair(1:nproma, jrow)                &
                           , ocean%zqsw(1:nproma,jrow)                 &
                           , ice%zqsi(1:nproma,jrow)                   &
                           , box%surf_vir_temp(1:nproma,jrow)          &
                           , box%surface_humidity(1:nproma,jrow)       &
                           , land%zqsl(1:nproma,jrow)                  &
                           , surface%is_land(1:nproma,jrow)            &
                           , surface%is_ocean(1:nproma,jrow)           &
                           , surface%is_seaice(1:nproma,jrow)         )
    !-------------------------------------------------------------------
    ! LAKE PHYSICS (old lake routine, 
    ! up to now located in mo_surface_ocean, mo_surface_ice)

    CALL s_lake( nproma                                                &
               , box%seaice(1:nproma,jrow)                             &
               , ice%ice_depth(1:nproma,jrow)                          &
               , palake(1:nproma)                         &
               , ice%surface_temperature(1:nproma,jrow)                &
               , ocean%surface_temperature(1:nproma,jrow)              &
               , ocean%latent_heat_flux_inst(1:nproma,jrow)            &
               , ocean%sensible_heat_flux_inst(1:nproma,jrow)          &
               , ice%qres(1:nproma,jrow)                               &
               , ocean%fluxres(1:nproma,jrow)                          &
               , ocean%trflw(1:nproma,jrow)                            &
               , ocean%soflw(1:nproma,jrow)                            &
               , ice%latent_heat_flux_inst(1:nproma,jrow)              &
               , ice%snow_water_equivalent(1:nproma,jrow)              &
               , ice%zcvsi(1:nproma,jrow)                      )

    CALL s_licetemp( nproma                                            &
                   , ice%ice_depth(1:nproma,jrow)                      &
                   , ice%snow_water_equivalent(1:nproma,jrow)          &
                   , palake(1:nproma)                     &
                   , ice%surface_temperature(1:nproma,jrow)            &
                   , ice%trfli(1:nproma,jrow)                          &
                   , ice%sofli(1:nproma,jrow)                          &
                   , ice%ahfice(1:nproma,jrow)                         &
                   , ice%qres(1:nproma,jrow)                           &
                   , ice%ahfcon(1:nproma,jrow)                         &
                   , ice%melting(1:nproma,jrow)                        &
                   , ice%evaporation_inst(1:nproma,jrow)               &
                   , snow(:)                                           &
                   , ice%sensible_heat_flux_inst(1:nproma,jrow)        &
                   , ice%latent_heat_flux_inst(1:nproma,jrow)          &
                   , ice%zcvsi(1:nproma,jrow)                          &
                   , box%seaice_fract(1:nproma,jrow)                   &
                   , ice%albedo(1:nproma,jrow)                         &
                   , ice%alsobs(1:nproma,jrow)                         &
                    )

  IF (lmlo) THEN    
    CALL ml_ocean ( nproma                                          &
                   , slf(1:nproma,jrow)                                &
                   , slm(1:nproma,jrow)                                &
!                   , box%land_fract(1:nproma,jrow)                     &
                   , lpnorth(1:nproma)                                 &
!                   , seaice(1,krow),    siced(1,krow),    alake(1,krow)      &
                   , box%seaice(1:nproma,jrow)                         &
                   , ice%ice_depth(1:nproma,jrow)                      &
                   , palake(1:nproma)                     &
!                   , tsi(1,krow),       tsw(1,krow)                          &
                   , ice%surface_temperature(1:nproma,jrow)            &
                   , ocean%surface_temperature(1:nproma,jrow)          &
!                   , zhflw,             zhfsw,            fluxres(1,krow)    &
                   , ocean%latent_heat_flux_inst(1:nproma,jrow)        &
                   , ocean%sensible_heat_flux_inst(1:nproma,jrow)      &
                   , ocean%fluxres(1:nproma,jrow)                      &
!                   , ztrflw,            zsoflw                               &
                   , ocean%trflw(1:nproma,jrow)                        &
                   , ocean%soflw(1:nproma,jrow)                        &
                   , ice%qres(1:nproma,jrow)                           &
                   , pamlcorr(1:nproma)                                &
                   , pamlcorac(1:nproma)                               &
                   , pamlheatac(1:nproma)                              &
!                   , zevapi,            sni(1,krow),      
                   , ice%evaporation_inst(1:nproma,jrow)               &
                   , ice%snow_water_equivalent(1:nproma,jrow)          &
!                   zcvsi              &
                   , ice%zcvsi(1:nproma,jrow)                          &
!                   , zfri           )
                   , box%seaice_fract(1:nproma,jrow) )
  END IF                 

    CALL s_sicetemp( nproma                                            &
                   , ice%ice_depth(1:nproma,jrow)                      &
                   , ice%snow_water_equivalent(1:nproma,jrow)          &
                   , palake(1:nproma)                     &
                   , slf(1:nproma,jrow)                                &
                   , ice%surface_temperature(1:nproma,jrow)            &
                   , ice%trfli(1:nproma,jrow)                          &
                   , ice%sofli(1:nproma,jrow)                          &
                   , ice%ahfice(1:nproma,jrow)                         &
                   , ice%qres(1:nproma,jrow)                           &
                   , ice%evaporation_inst(1:nproma,jrow)               &
                   , ice%ahfcon(1:nproma,jrow)                         &
                   , ice%melting(1:nproma,jrow)                        &
                   , ice%sensible_heat_flux_inst(1:nproma,jrow)        &
                   , ice%latent_heat_flux_inst(1:nproma,jrow)          &
                   , snow(:)                                           &
                   , ice%zcvsi(1:nproma,jrow)                          &
                   , box%seaice_fract(1:nproma,jrow)                   &
                   , ice%albedo(1:nproma,jrow)                        &
                   , ice%alsobs(1:nproma,jrow)                        &
                    )
    !-------------------------------------------------------------------
    ! Average surface temperature
    CALL surface_box_average( nproma                                   &
                       , land%surface_temperature_new(1:nproma,jrow)   &
                       , ice%surface_temperature(1:nproma,jrow)        &
                       , ocean%surface_temperature(1:nproma,jrow)      &
                       , box%surface_temperature(1:nproma,jrow)        &
                       , surface%is_land(1:nproma,jrow)                &
                       , surface%is_ocean(1:nproma,jrow)               &
                       , surface%is_seaice(1:nproma,jrow)              &
                       , box%land_fract(1:nproma,jrow)                 &
                       , box%ocean_fract(1:nproma,jrow)                &
                       , box%seaice_fract(1:nproma,jrow)              )
    box%surface_temperature_acc(1:nproma,jrow) =                       &
                  box%surface_temperature_acc(1:nproma,jrow)           &
                + box%surface_temperature(1:nproma,jrow) * delta_time
    tsurf_na(1:nproma,jrow) = box%surface_temperature(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! Radiative Temperature for subroutine RADIATION (from physc)
    CALL surface_box_average( nproma                                   &
                    , (land%surface_temperature_new(1:nproma,jrow)**4) &
                    , (ice%surface_temperature(1:nproma,jrow)**4)      &
                    , (ocean%surface_temperature(1:nproma,jrow)**4)    &
                    , radtemp(1:nproma)                                &
                    , surface%is_land(1:nproma,jrow)                   &
                    , surface%is_ocean(1:nproma,jrow)                  &
                    , surface%is_seaice(1:nproma,jrow)                 &
                    , box%land_fract(1:nproma,jrow)                    &
                    , box%ocean_fract(1:nproma,jrow)                   &
                    , box%seaice_fract(1:nproma,jrow)                 )
    !-------------------------------------------------------------------
    ! Radiative Temperature for subroutine RADHEAT (from physc)
    CALL surface_box_average( nproma                                   &
                       , (land%surface_temperature(1:nproma,jrow)**4)  &
                       , (ice%surface_temperature(1:nproma,jrow)**4)   &
                       , (ocean%surface_temperature(1:nproma,jrow)**4) &
                       , radtemp_old(1:nproma)                         &
                       , surface%is_land(1:nproma,jrow)                &
                       , surface%is_ocean(1:nproma,jrow)               &
                       , surface%is_seaice(1:nproma,jrow)              &
                       , box%land_fract(1:nproma,jrow)                 &
                       , box%ocean_fract(1:nproma,jrow)                &
                       , box%seaice_fract(1:nproma,jrow)              )
    land%surface_temperature(1:nproma,jrow)              =             &
                            land%surface_temperature_new(1:nproma,jrow)
    !-------------------------------------------------------------------
    ! MASK declaration for specific surface part
    ! (pingo masks for maskout)
    box%landpingo                                        = 0._wp
    box%oceanpingo                                       = 0._wp
    box%icepingo                                         = 0._wp
    WHERE(surface%is_land) box%landpingo                 = 1._wp
    WHERE(surface%is_ocean ) box%oceanpingo = 1._wp
    WHERE(surface%is_seaice ) box%icepingo   = 1._wp
    !-------------------------------------------------------------------
    ! Variables for ECHAM (updates Parameters and old stream vars ect.
    zth(1:nproma)    = (radtemp(1:nproma) )**0.25_wp
    pradtemp_old(1:nproma) = (radtemp_old(1:nproma))**0.25_wp
    ztvh(1:nproma)   = box%surf_vir_temp(1:nproma,jrow)
    zqsurf(1:nproma) = box%surface_humidity(1:nproma,jrow) 
    ptsw_new(:)      = ocean%surface_temperature(1:nproma,jrow)
    palsol(:)        = land%albedo (1:nproma,jrow)
    palsoi(:)        = ice%albedo  (1:nproma,jrow)
    palsow(:)        = ocean%albedo(1:nproma,jrow)
    palsobs(:)       = ice%alsobs(1:nproma,jrow)
    palsom(:)        = ice%alsom(1:nproma,jrow)
    palbedo(:)       = box%albedo(1:nproma,jrow)
    palbedo_vis(:)   = box%albedo_vis(1:nproma,jrow)
    palbedo_vis_dir(:) = box%albedo_vis_dir(1:nproma,jrow)
    palbedo_vis_dif(:) = box%albedo_vis_dif(1:nproma,jrow)
    palbedo_nir(:)   = box%albedo_nir(1:nproma,jrow)
    palbedo_nir_dir(:) = box%albedo_nir_dir(1:nproma,jrow)
    palbedo_nir_dif(:) = box%albedo_nir_dif(1:nproma,jrow)
    ptrflw(:)        = ocean%trflw(1:nproma,jrow)
    ptrfli(:)        = ice%trfli(1:nproma,jrow)
    psofll(:)        = land%sofll(1:nproma,jrow)
    psoflw(:)        = ocean%soflw(1:nproma,jrow)
    psofli(:)        = ice%sofli(1:nproma,jrow)   
    ptrfllac(:)      = land%trfllac(1:nproma,jrow)
    ptrflwac(:)      = ocean%trflwac(1:nproma,jrow)
    ptrfliac(:)      = ice%trfliac(1:nproma,jrow)
    psofllac(:)      = land%sofllac(1:nproma,jrow)
    psoflwac(:)      = ocean%soflwac(1:nproma,jrow)
    psofliac(:)      = ice%sofliac(1:nproma,jrow)
    pwind10w(:)      = ocean%wind_10_meter(1:nproma,jrow)
    pu10(:)          = box%u_wind_10_meter(1:nproma,jrow)   
    pv10(:)          = box%v_wind_10_meter(1:nproma,jrow)
    pwimax(:)        = box%wind_10_max(1:nproma,jrow)
    pwind10(:)       = box%wind_10_meter_acc(1:nproma,jrow)
    pdew2(:)         = box%dewpoint_2_meter(1:nproma,jrow)    
    ptemp2(:)        = box%temp_2_meter(1:nproma,jrow) 
    pt2max(:)        = box%maxtemp_2_meter(1:nproma,jrow)
    pt2min(:)        = box%mintemp_2_meter(1:nproma,jrow)
    pevaplac(:)      = land%evaporation_acc(1:nproma,jrow)
    pevapwac(:)      = ocean%evaporation_acc(1:nproma,jrow) 
    pevapiac(:)      = ice%evaporation_acc(1:nproma,jrow)
    pevap(:)         = box%evaporation_acc(1:nproma,jrow)
    pahfllac(:)      = land%latent_heat_flux_acc(1:nproma,jrow)  
    pahflwac(:)      = ocean%latent_heat_flux_acc(1:nproma,jrow)
    pahfliac(:)      = ice%latent_heat_flux_acc(1:nproma,jrow)   
    pahfl(:)         = box%latent_heat_flux_acc(1:nproma,jrow)
    pahfslac(:)      = land%sensible_flux_acc(1:nproma,jrow) 
    pahfswac(:)      = ocean%sensible_flux_acc(1:nproma,jrow) 
    pahfsiac(:)      = ice%sensible_flux_acc(1:nproma,jrow) 
    pahfs(:)         = box%sensible_heat_flux_acc(1:nproma,jrow) 
    pqhfla(:)        = box%evaporation_inst(1:nproma,jrow)  
    pevapw(:)        = ocean%evaporation_inst(1:nproma,jrow)  
    pevapi(:)        = ice%evaporation_inst(1:nproma,jrow) 
    pahfsl(:)        = land%sensible_heat_flux_inst(1:nproma,jrow) 
    pahfsw(:)        = ocean%sensible_heat_flux_inst(1:nproma,jrow)
    pahfsi(:)        = ice%sensible_heat_flux_inst(1:nproma,jrow)
    pevapot(:)       = land%evaporation_pot(1:nproma,jrow)
    pahflw(:)        = ocean%latent_heat_flux_inst(1:nproma,jrow)
    pahfli(:)        = ice%latent_heat_flux_inst(1:nproma,jrow)
    pahfll(:)        = land%latent_heat_flux_inst(1:nproma,jrow)
    psni(:)          = ice%snow_water_equivalent(1:nproma,jrow)
    pahfice(:)       = ice%ahfice(1:nproma,jrow)
    pfluxres(:)      = ocean%fluxres(1:nproma,jrow)
    pqres(:)         = ice%qres(1:nproma,jrow)
    pahfcon(:)       = ice%ahfcon(1:nproma,jrow)  
    pahfres(:)       = ice%melting(1:nproma,jrow)
    ptsi(:)          = ice%surface_temperature(1:nproma,jrow)
    ptslnew(:)       = land%surface_temperature_rad(1:nproma,jrow)
    pzti(:)          = radtemp(1:nproma)**0.25_wp
    pztsnew(:)       = ( box%land_fract  (1:nproma,jrow) * pzteffl4(1:nproma) + &
                  box%seaice_fract(1:nproma,jrow) * ice%surface_temperature(1:nproma,jrow)**4 + &
                  box%ocean_fract (1:nproma,jrow) * ocean%surface_temperature(1:nproma,jrow)**4 ) &
                  **0.25_wp
    ptsurf(:)        = box%surface_temperature_acc(1:nproma,jrow)
    paz0w(:)         = ocean%roughness(1:nproma,jrow)
    paz0i(:)         = ice%roughness(1:nproma,jrow)
    paz0l(:)         = land%roughness_momentum(1:nproma,jrow)
    paz0(:)          = box%roughness(1:nproma,jrow) 
    pustrl(:)        = land%u_stress(1:nproma,jrow)
    pvstrl(:)        = land%v_stress(1:nproma,jrow) 
    pustrw(:)        = ocean%u_stress(1:nproma,jrow)
    pvstrw(:)        = ocean%v_stress(1:nproma,jrow) 
    pustri(:)        = ice%u_stress(1:nproma,jrow)
    pvstri(:)        = ice%v_stress(1:nproma,jrow)
    pustr(:)         = box%u_stress_acc(1:nproma,jrow)
    pvstr(:)         = box%v_stress_acc(1:nproma,jrow)
    ptte_corr(:)     = tte_corr(:)/                                    &
                        (atm_half_lev_press(:,levels_plus_1)           &
                        -atm_half_lev_press(:,levels))
    pseaice_new(:)   = box%seaice(1:nproma,jrow)
    psiced_new(:)    = ice%ice_depth(1:nproma,jrow)
    pco2_flux(:)     = MERGE(box%land_fract (1:nproma,jrow) * pco2_flux_land(1:nproma), &
                             zero, surface%is_land (1:nproma,jrow)) + &
                       MERGE((box%ocean_fract(1:nproma,jrow)+box%seaice_fract(1:nproma,jrow)) * &
                        pco2_flux_ocean(1:nproma), zero, surface%is_ocean(1:nproma,jrow))

!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!   ztvl  !from subroutine precalc_land  
!   zcdnl !from subroutine precalc_land 

! TS    zcfncl(1:nproma) = land%ZCFNCL(1:nproma,jrow)
! TS    zril(1:nproma)   = land%ZRIL(1:nproma,jrow)
!>>SF gf #78
#ifdef HAMMOZ
    zril(1:nproma)   = land%ZRIL(1:nproma,jrow)
    zriw(1:nproma)   = ocean%ZRIW(1:nproma,jrow)
    zrii(1:nproma)   = ice%ZRII(1:nproma,jrow)
    zcfml(1:nproma)  = land%zcfml(1:nproma,jrow)
    zcfmw(1:nproma)  = ocean%zcfmw(1:nproma,jrow)
    zcfmi(1:nproma)  = ice%zcfmi(1:nproma,jrow)     
#endif
!<<SF gf #78
    zcfnc(1:nproma)  = box%ZCFNC(1:nproma,jrow)
!    pfrl(1:nproma)   = box%land_fract(1:nproma,jrow)  
!    pfrw(1:nproma)   = box%ocean_fract(1:nproma,jrow)  
!    pfri(1:nproma)   = box%seaice_fract(1:nproma,jrow)  
    pcvs(1:nproma)   = ice%snow_cover_fract(1:nproma,jrow)
    zsrfll(1:nproma) = land%sofll(1:nproma,jrow)   
    velo10m(1:nproma)= box%wind_10_meter(1:nproma,jrow)  

!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Variables for g3b stream only for diagnose !
    friac(1:nproma,jrow)    = box%frac_ice_cover_acc(1:nproma,jrow)
    tslm1(1:nproma,jrow)    = land%surface_temperature_new(1:nproma,jrow) 
    ws(1:nproma,jrow)       = soil_wetness(1:nproma)
    sn(1:nproma,jrow)       = snow_depth(1:nproma)
    wl(1:nproma,jrow)       = skin_reservoir(1:nproma)

  END SUBROUTINE update_surface

  SUBROUTINE init_surface

    USE mo_time_control, ONLY: lresume

    IF (.NOT. lresume) THEN
       land%roughness_momentum(:,:)   = 0.1_wp
       land%roughness_heat(:,:)       = 0.1_wp
       ocean%roughness(:,:)           = 0.001_wp
       ice%roughness(:,:)             = 0.001_wp

       ice%snow_water_equivalent(:,:) = 0.0_wp
       land%zcair(:,:)                = 0.0_wp
       land%zcsat(:,:)                = 0.0_wp
    END IF

  END SUBROUTINE init_surface

END MODULE mo_surface