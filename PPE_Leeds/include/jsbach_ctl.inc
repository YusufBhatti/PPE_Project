!*    NAMELIST *JSBACH_CTL* To control the land surface scheme JSBACH.
!                                                                 
!      
!      NAME                 TYPE      PURPOSE                                       
!      ----                 ----      -------                                       
!                                                                       
!     *NTILES*              INTEGER   Number of tiles defined on each grid cell
!     *STANDALONE*          LOGICAL   Type of model run 
!                                       standalone: stand-alone jsbach run
!                                       coupled: jsbach driven by an atmosphere model
!     *COUPLING*            CHARACTER Type of coupling: explicit, semi, implicit
!     *LSS*                 CHARACTER Land surface sceme: ECHAM
!     *USE_BETHY*           LOGICAL   Use BETHY model (photosynthesis, respiration)
!     *USE_PHENOLOGY*       LOGICAL   Calculate LAI using the phenology module
!     *PHENO_SCHEME*        CHARACTER Phenonlgy scheme: LOGROP or KNORR
!     *USE_ALBEDO*          LOGICAL   Calculate albedo depending on vegetaion
!     *USE_ROUGHNESS_LAI*   LOGICAL   Calculate z0 depending on LAI
!     *USE_ROUGHNESS_ORO*   LOGICAL   Calculate z0 including subgrid-scale topographie
!     *WITH_NITROGEN*       LOGICAL   If .TRUE. nitrogen will be cycled in addition to carbon
!                                       Otherwise: only carbon cycle is simulated.
!     *WITH_YASSO*          LOGICAL      If .TRUE. the yasso model is used for litter and soil carbon decomposition 
!     *WITH_HD*             LOGICAL   Use hydrological discharge (HD) model
!     *USE_DYNVEG*          LOGICAL   Use the dynamic vegetation module
!     *USE_DISTURBANCE*     LOGICAL   Use of disturbance module independent of the dynamic vegetation
!     *LCC_FORCING_TYPE*    CHARACTER Scheme for landuse changes: none, maps, transitions
!     *LCC_SCHEME*          INTEGER   Pool scheme for landuse change: 1: litter (org) 2: anthropogenic
!     *FILE_TYPE*           INTEGER   Output format: 1=GRIB, 2=NETCDF, 6=NETCDF4
!     *FILE_ZTYPE*          INTEGER   Output compression type: GRIB SZIP=1, NETCDF4 ZIP=2
!     *OUT_STATE*           LOGICAL   Write the jsbach stream
!     *VEG_AT_1200*         LOGICAL   Write veg stream at 12:00 each day
!     *LPOST_ECHAM*         LOGICAL   Write jsbach output variables, even if they
!                                       are includeded in the echam output
!     *MISSING_VALUE*       REAL      Missing value for the output
!     *TEST_STREAM*         LOGICAL   Additional stream for model testing
!     *DEBUG*               LOGICAL   Additional output for debugging
!     *DEBUG_CCONSERVATION* LOGICAL   Additional output to debug carbon conservation
!     *TEST_CCONSERVATION*  LOGICAL   Activate carbon conservation test
!     *GRID_FILE*           CHARACTER File containing grid information
!     *LCTLIB_FILE*         CHARACTER Name of the land cover library file
!     *VEG_FILE*            CHARACTER File containig initial data for the vegetation 
!     *SURF_FILE*           CHARACTER File containig initial data for the surface
!     *SOIL_FILE*           CHARACTER File containig initial data for the soil
!     *READ_COVER_FRACT*    LOGICAL   read cover fractions from file
!     *INPUT_VERBOSE*       CHARACTER Debug flags for mo_input
!     *WRITE_INTERFACE_VARS*  LOGICAL   write all input variables goinig through the jsbach interface
!     *READ_INTERFACE_VARS*   LOGICAL   read stepwise interface variables from file
!     *INTERFACE_TEST*      LOGICAL   setup that allows for bit-identical results of echam/jsbach with 
!                                     WRITE_INTEFACE_VARS and jsbach offline with READ_INTERFACE_VARS

NAMELIST /JSBACH_CTL/ &
      NTILES, STANDALONE, LSS, USE_BETHY, USE_PHENOLOGY, USE_ALBEDO, &
      USE_ROUGHNESS_LAI, USE_ROUGHNESS_ORO, &
      WITH_NITROGEN, WITH_YASSO, WITH_HD, USE_DYNVEG, USE_DISTURBANCE, &
      LCC_FORCING_TYPE, LCC_SCHEME, PHENO_SCHEME, &
      FILE_TYPE, FILE_ZTYPE, OUT_STATE, VEG_AT_1200, LPOST_ECHAM, &
      READ_COVER_FRACT, MISSING_VALUE, TEST_STREAM, DEBUG, GRID_FILE, &
      LCTLIB_FILE, VEG_FILE, SURF_FILE, SOIL_FILE, DEBUG_CCONSERVATION, &
      TEST_CCONSERVATION, INPUT_VERBOSE, WRITE_INTERFACE_VARS, READ_INTERFACE_VARS, &
      INTERFACE_TEST
!                                                                       
!     ------------------------------------------------------------------
