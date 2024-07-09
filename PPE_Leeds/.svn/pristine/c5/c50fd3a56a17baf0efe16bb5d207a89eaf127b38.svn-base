!*    NAMELIST *FIRE_THONICKE_CTL* To control Thonicke fire algorithms.
!                                                                 
!      
!      NAME                TYPE       PURPOSE                                       
!      ----                ----       -------                                       
!                                                                       
!     *a_nd               REAL       propensity of humans to ignite a fire
!     *llight_ground      LOGICAL    file contains cloud to ground flashes only?
!     *lwind_speed_limit  LOGICAL    wind speed is limited
!     *human_para         REAL       parameter to tune the humen ignitions
!     *ign_para           REAL       parameter to tune the sum of ignitions
!     *wind_limit         REAL       parameter for the wind limitation
!     *SurfArea2Vol       REAL(3)    surface area to volume ration for 1hr 10hr and 100hr fuels
!     *moisture_scaling   REAL       factor to convert SAtVR to fuel drying parameter
!     *wind_max           REAL       maximum effective windspeed
!     *wind_slope         REAL       determines the reduction of high windspeed
!     *read_fuel_frac     LOGICAL    read fractions of fuel from file during initialization
!     *lduration_popd     LOGICAL    fire duration depends on population density
!     *lcalc_frp          LOGICAL    calculate total FRP of fires (requires extra input, default:
!                                                                  cosmos-landveg     : .false.
!                                                                  cosmos-landveg-fire: .true.
!     *Lethal2ResTime     REAL       Lethal heat time to fire residence time

NAMELIST /FIRE_THONICKE_CTL/                                            &
      llight_ground, lwind_speed_limit, a_nd, human_para, ign_para,     &
      wind_limit, SurfArea2Vol, moisture_scaling, wind_max, wind_slope, &
      read_fuel_frac, lduration_popd, lcalc_frp, Lethal2ResTime
!     ------------------------------------------------------------------
