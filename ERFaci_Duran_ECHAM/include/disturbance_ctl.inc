!*    NAMELIST *DISTURBANCE_CTL* To control JSBACH vegetation disturbance algorithms.
!                                                                 
!      
!      NAME                TYPE       PURPOSE                                       
!      ----                ----       -------                                       
!                                                                       
!     *FIRE_NAME*          CHARACTER  Name of used fire algorithm. Supported: none, jsbach, thonicke
!     *FUEL_NAME*          CHARACTER  Name of fuel consumption algorithm. Supported: none, jsbach, thonicke
!     *WINDBREAK_NAME*     CHARACTER  Name of windbreak algorithm. Suppored names: None, JSBach
!     *LBURN_PASTURE*      LOGICAL    Allow also pasture to be burned
!     *LDIAG*              LOGICAL    Produce additional diagnostic output if possible (thonicke fire only)
!     *FIRE_FRAC_WOOD_2_   REAL       fraction of wood carbon emitted to the atmosphere by fire
!      ATMOS
!
NAMELIST /DISTURBANCE_CTL/             &
 fire_name, fuel_name, windbreak_name, & ! Algorithm selection by name
 ldiag,                                & ! Diagnostic output
 fire_frac_wood_2_atmos, lburn_pasture
!     ------------------------------------------------------------------
