!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!!++mgs: see lines 62-64 and 100
MODULE mo_surface_types 

  USE mo_kind,             ONLY: dp

PUBLIC
  
  TYPE grid_box_mean_type
     REAL(dp), POINTER, DIMENSION(:,:) ::  &        !! derived properties for the lowest atm level:
          atm_tot_cloud_water,             &        !! total cloud water[
          atm_sat_spec_hum,                &        !!
          atm_pot_temp,                    &        !!
          atm_vir_pot_temp,                &        !!
          atm_lat_heat_fact,               &        !!
          surf_vir_temp,                   &        !!
          evaporation_inst,                &        !! Moisture flux from surface to atmosphere
          momentum_ex_coef,                &        !! momentum exchange coefficient (zcfm -> zcdum in vdiff)
          heat_ex_coef,                    &        !! heat exchange coefficient (zcthv -> cduh in vdiff)
          ustar,                           &        !! u_star mean for tke calculation in vdiff
          ztkevn,                          &        !! tke boundary condition
          u_stress,                        &        !! pustr for accumulation in vdiff
          v_stress,                        &        !! pvstr for accumulation in vdiff
          wind_10_meter,                   &        !! wind speed 10 meter above displacement height
          sensible_heat_flux_inst,         &        !! sensible heat flux   
          latent_heat_flux_inst,           &        !! latent heat flux
          temp_2_meter,                    &        !! 2 meter temperature
          dewpoint_2_meter,                &        !! dew point 2 meter height
          surface_temp_rad,                &        !! Radiative surface temperature
          albedo,                          &        !! Surface albedo
          albedo_vis,                      &        !! Surface albedo visible range
          albedo_nir,                      &        !! Surface albedo NIR range
          albedo_vis_dir,                  &        !! Surface albedo visible range (dir)
          albedo_nir_dir,                  &        !! Surface albedo NIR range (dir)
          albedo_vis_dif,                  &        !! Surface albedo visible range (dif)
          albedo_nir_dif,                  &        !! Surface albedo NIR range (dif)
          surface_temperature,             &        !! Surface temperature
          lake_fract,                      &        !! Fraction of Lakes
          frac_ice_cover_acc,              &        !! Ice cover (fraction of grid box) (Accumulated, acc. to 'friac' in old echam5)
          longwave_down_acc,               &        !! Longwave radiation downward accumulated 
          surface_qsat,                    &        !! Saturated surface specific humidity
          surface_humidity,                &        !!
          atm_liq_wat_pot_temp,            &        !!
          land_fract,                      &        !! Fraction of land [fraction of grid box]
          seaice_fract,                    &        !! Fraction of sea ice [fraction of brid box] acc. to zfri
          seaice,                          &        !! seaice cover (1-SLM) acc. to seaice in echam
          ocean_fract,                     &        !! Fraction of ocean [fraction of grid box]
          atm_dry_stat_energy,             &
          sensible_heat_flux_acc,          &
          latent_heat_flux_acc,            &
          u_wind_10_meter,                 &
          v_wind_10_meter,                 &
          surface_temperature_acc,         &
          wind_10_meter_acc,               &
          roughness,                       &
          u_stress_acc,                    &
          v_stress_acc,                    &
          evaporation_acc,                 &
          maxtemp_2_meter,                 &
          mintemp_2_meter,                 &
          landpingo,                       &
          oceanpingo,                      &
          icepingo,                        &
          wind_10_max,                     &
          zcfnc,                           &        !! function of heat transfer coeff., TS new
          zri                                       !! moist Richardson number,          TS new
          
       END TYPE grid_box_mean_type

  TYPE surface_type
     LOGICAL, POINTER,  DIMENSION(:,:) ::  &        !! 
          is_land,                         &        !! Grid box contains land?
          is_seaice,                       &        !! Grid box contains sea ice?
          is_ocean,                        &        !! Grid box contains ocean?
          is_ocean_for_coup,               &        !! Gric box must additionally be calculated for ocean-coupling
          is_ice_for_coup                           !! Gric box must additionally be calculated for ocean-coupling
  END TYPE surface_type
  
  TYPE land2atmos_type
     REAL(dp), POINTER, DIMENSION(:,:) ::  &        !! Variables for land surface atmos exchange:
          surface_temperature_rad,         &        !! Radiative land surface temperature
          albedo,                          &        !! Land surface albedo
          albedo_vis,                      &        !! Land surface albedo visible range
          albedo_nir,                      &        !! Land surface albedo NIR range
          surface_temperature,             &        !! Surface temperature
          surface_temperature_new,         &        !! new surf_temp calculated from soil model within the time step
!! surf_temp according to ptslm1 in vdiff !!
          dry_static_energy_new,           &        !! 
          surface_qsat,                    &        !! Saturated surface specific humidity
          evaporation_inst,                &        !! Moisture flux from surface to atmosphere
          roughness_heat,                  &        !! roughness length of heat
          roughness_momentum,              &        !! roughness length of momentum (according to paz0l)
          evaporation_pot,                 &        !! potential evaporation
          sensible_heat_flux_inst,         &        !! sensible heat flux 
          latent_heat_flux_inst,           &        !! latent heat flux
          zcptlnew,                        &        !! new static energy land surface
          zqslnew,                         &        !! new satrated surface air moisture
          zdqsl,                           &        !! DERIV. OF SAT.SPEC.HUMIDITY AT THE OLD TEMP(dq_sat/dT)
          zril,                            &        !! moist richardson number
          zqsl,                            &        !! sat. mix. 'verhältnis' water vapour
!! TS     zcfncl,                          &        !! function of heat transfer coeff.
          zchl,                            &        !! heat transfer coeff. without function factor 
          zcfhl,                           &        !! stability dependent tranfer coeff. for heat
          zbnl,                            &        !! interpolation function for diagnose
          zbhnl,                           &        !! interpolation function for diagnose - maybe wrong ?
          zbml,                            &        !! interpolation function for diagose
          zbhl,                            &        !! interpolation function for diagnose
          zustarl,                         &        !! schubspannungsgeschwindigkeit (u*)
          zcfml,                           &        !! stability dependend transfer coeff. for momentum
          ztkevl,                          &        !! tke boundary condition
          zetnl,                           &        !! R-M-Coeff 
          zftnl,                           &        !! R-M-Coeff 
          zeqnl,                           &        !! R-M-Coeff 
          zfqnl,                           &        !! R-M-Coeff 
          zustl,                           &        !!
          zspeedl,                         &        !!
          zcair,                           &        !!
          zcsat,                           &        !!
          temp_2_meter,                    &        !!
          dewpoint_2_meter,                &        !!
          ztklevl,                         &        !!
          zqklevl,                         &        !!
          z0h,                             &        !! according to z0h in old vdiff
          zcptl,                           &
          zcpq,                            &
          surface_qsat_new,                &
          u_stress,                        &
          v_stress,                        &
          latent_heat_flux_acc,            &
          evaporation_acc,                 &
          sensible_flux_acc,               &
          u_wind_10_meter,                 &
          v_wind_10_meter,                 &
          wind_10_meter,                   &
          ZCAIR_OLD,                       &
          trfll,                           &
          sofll,                           &
          trfllac,                         &
          sofllac

  END TYPE land2atmos_type
  
  TYPE ocean2atmos_type
     REAL(dp), POINTER, DIMENSION(:,:) ::  &   !! Variables for ocean - atmos exchange:
          surface_temperature,             &        !! ocean surface temperature
          albedo,                          &        !! Ocean albedo
          albedo_vis,                      &        !! Ocean albedo visible range
          albedo_nir,                      &        !! Ocean albedo NIR range
          albedo_vis_dir,                  &        !! Surface albedo visible range (dir)
          albedo_nir_dir,                  &        !! Surface albedo NIR range (dir)
          albedo_vis_dif,                  &        !! Surface albedo visible range (dif)
          albedo_nir_dif,                  &        !! Surface albedo NIR range (dif)
          evaporation,                     &        !! Moisture flux from surface to atmosphere
          zqsw,                            &        !! sat. mix. 'verhältnis' water vapour
          zcptw,                           &        !! CpT
          zriw,                            &        !! moist richardson number
          zcfhw,                           &        !! stability dependent tranfer coeff. for heat
          zchw,                            &        !! heat transfer coeff. without function factor 
          zbnw,                            &        !! interpolation function for diagnose
          zbmw,                            &        !! interpolation function for diagnose
          zbhw,                            &        !! interpolation function for diagnose
          zustarw,                         &        !! schubspannungsgeschwindigkeit (u*)
          ztkevw,                          &        !! tke boundary condition
          ztklevw,                         &        !! 
          zqklevw,                         &        !!
          zcfmw,                           &        !! stability dependend transfer coeff. for momentum
          zetnw,                           &        !! R-M-Coeff 
          zftnw,                           &        !! R-M-Coeff 
          zeqnw,                           &        !! R-M-Coeff 
          zfqnw,                           &        !! R-M-Coeff 
          evaporation_inst,                &        !! For routine collect in physc:ocean coupling
          sensible_heat_flux_inst,         &        !!
          latent_heat_flux_inst,           &        !!
          temp_2_meter,                    &        !! Wind speed 2m over water
          roughness,                       &        !! 
          zustw,                           &        !!
          zspeedw,                         &        !!
          dewpoint_2_meter,                &        !! 
          fluxres,                         &        !! heat flux residuum from lake ice melting
          ahfres,                          &
          PHFLW,                           &
          trflw,                           &
          soflw,                           &
          trflwac,                         &
          soflwac,                         &
          PHFSW,                           &
          u_stress,                        &
          v_stress,                        &
          latent_heat_flux_acc,            &
          evaporation_acc,                 &
          sensible_flux_acc,               &
          u_wind_10_meter,                 &
          v_wind_10_meter,                 &
          wind_10_meter
     
  END TYPE ocean2atmos_type
  
  TYPE ice2atmos_type
     REAL(dp), POINTER, DIMENSION(:,:) ::  &        !! Variables for atmos-sea ice exchange:
          snow_cover_fract,                &        !! Fractional coverage of seaice by snow
          surface_temperature,             &        !! Surface temperature over ice [K]
          surface_temp_rad,                &        !! Radiative sea ice temperature
          snow_water_equivalent,           &        !! Water equivalent of snow on sea ice
          depth,                           &        !! Depth of sea ice
          albedo,                          &        !! Sea ice albedo
          albedo_vis,                      &        !! Sea ice albedo visible range
          albedo_nir,                      &        !! Sea ice albedo NIR range
          albedo_vis_dir,                  &        !! Surface albedo visible range (dir)
          albedo_nir_dir,                  &        !! Surface albedo NIR range (dir)
          albedo_vis_dif,                  &        !! Surface albedo visible range (dif)
          albedo_nir_dif,                  &        !! Surface albedo NIR range (dif)
          evaporation,                     &        !! Moisture flux from surface to atmosphere
          zqsi,                            &        !! sat. mix. 'verhältnis' water vapour
          zcpti,                           &        !! CpT
          zrii,                            &        !! moist richardson number
          zcfmi,                           &        !! stability dependend transfer coeff. for momentum
          zchi,                            &        !! heat transfer coeff. without function factor 
          zcfhi,                           &        !! stability dependent tranfer coeff. for heat
          zbni,                            &        !! interpolation function for diagnose
          zbmi,                            &        !! interpolation function for diagnose
          zbhi,                            &        !! interpolation function for diagnose
          zustari,                         &        !! schubspannungsgeschwindigkeit (u*)
          ztkevi,                          &        !! tke boundary condition
          ztklevi,                         &        !! 
          zqklevi,                         &        !!
          zetni,                           &        !! R-M-Coeff 
          zftni,                           &        !! R-M-Coeff 
          zeqni,                           &        !! R-M-Coeff 
          zfqni,                           &        !! R-M-Coeff
          evaporation_inst,                &        !!
          sensible_heat_flux_inst,         &        !!
          latent_heat_flux_inst,           &        !!
          temp_2_meter,                    &        !!
          roughness,                       &        !!
          zusti,                           &        !!
          dewpoint_2_meter,                &        !!
          ice_depth,                       &
          SNOW_ON_ICE,                     &
          PCVSI,                           &
          lw_flux,                         &
          sw_flux,                         &
          fraction,                        &
          v_stress,                        &
          u_stress,                        &
          latent_heat_flux_acc,            &
          evaporation_acc,                 &
          sensible_flux_acc,               &
          u_wind_10_meter,                 &
          v_wind_10_meter,                 &
          wind_10_meter,                   &
          qres,                            &
          melting ,                        &        !! ahfres in physc
          zcvsi,                           &
          ahfice,                          &
          ahfcon,                          &
          trfli,                           &
          sofli,                           &
          trfliac,                         &
          sofliac,                         &
          ameltdepth,                      &        !! melt pond depth
          ameltfrac,                       &        !! melt pond fraction
          alsom,                           &        !! albedo of melt ponds
          alsobs                                    !! albedo of bare ice and snow without ponds
  END TYPE ice2atmos_type

 END MODULE mo_surface_types
