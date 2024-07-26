!*    NAMELIST *WINDBREAK_JSBACH_CTL* To control JSBACH vegetation disturbance algorithms.
!                                                                 
!      
!      NAME                TYPE       PURPOSE                                       
!      ----                ----       -------                                       
!                                                                       
!     *WIND_THRESHOLD*     REAL       factor by which the previous day maximum wind speed must be larger
!                                     than the climatological daily maximum wind speed to allow any wind damage
!     *WIND_DAMAGE_SCALE*  REAL       scaling factor for wind damage
!

NAMELIST /WINDBREAK_JSBACH_CTL/    &
 wind_threshold, wind_damage_scale
!     ------------------------------------------------------------------
