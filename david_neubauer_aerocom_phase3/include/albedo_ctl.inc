!*    NAMELIST *ALBEDO_CTL* To control land surface albedo.
!                                                                 
!      
!      NAME     TYPE      PURPOSE                                       
!      ----     ----      -------                                       
!                                                                       
!     *USE_ALBEDOCANOPY*    LOGICAL  true: use map of canopy albedo 
!                                    false: use PFT specific albedo values 
!     *USE_ALBEDOSOIL*      LOGICAL  true: calculate albedo of soil surface denpending on soil carbon and litter
!     *USE_ALBEDOSOILCONST* LOGICAL  true: the base albedo of the soil surface (without soil carbon and leaf litter)
!                                          is set to a global constant
!                                    false: the base albedo of the soil surface as read from the jsbach.nc file is used
!     *USE_SOC* CHARACTER(len=6)     linear: albedo of the soil is linearly reduced by soil carbon
!                                    log: logarithmic dependence of the albedo of the soil surface on soil carbon
!     *USE_LITER*         LOGICAL    true: albedo of the soil surface depends on leaf litter
!     *ALBEDO_AGE_WEIGHT  REAL(dp)  0: ECHAM5 scheme for albedo of snow is used
!                                   1: snow age scheme is used
!                                      0 < albedo_age_weight < 1: albedo of snow is calcualted by linearly weighting
!                                      the snow albedo of both schemes
   NAMELIST /ALBEDO_CTL/ USE_ALBEDOCANOPY, USE_ALBEDOSOIL, USE_ALBEDOSOILCONST, USE_SOC, USE_LITTER, ALBEDO_AGE_WEIGHT
!                                                                       
!     ------------------------------------------------------------------
