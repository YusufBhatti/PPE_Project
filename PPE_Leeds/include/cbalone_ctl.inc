!*    NAMELIST *CBALONE_CTL* To control (stand alone) carbon balance
!                                                                 
!      
!    NAME     TYPE      PURPOSE                                       
!    ----     ----      -------                                       
!                                                                       
!   *DRIVER_DATA_PATH*      CHARACTER    path to the driver data
!   *EXPERIMENT*            CHARACTER    name of the experiment providing driving data
!   *SPECIFICNAME*          CHARACTER    driving data identifier: filenames are
!                                        assumed to be EXPERIMENT_yyyymm.SPECIFICNAME.nc
!   *CLIMATE_YEARSTART*     INTEGER      first year of the driving data 
!   *CLIMATE_YEAREND*       INTEGER      last year of the driving data
!   *RUN_YEAR_FIRST*        INTEGER      first year of the cbalance experiment
!   *RUN_YEAR_LAST*         INTEGER      last year of the cbalance experiment
!   *OUT_INTERVAL*          CHARACTER    output interval: YEAR / MONTH / DAY
!   *INPUT_SCALING*         LOGICAL      true for CO2-scaling
!   *INPUT_SCALING_NDEP*    LOGICAL      true for scaling of nitrogen deposition
!   *CO2_FILE*              CHARACTER    CO2 development file - only with input_scaling
!   *CLIMATOLOGY_DIFF_FILE* CHARACTER    Climatological differences of driver data
!                                         - only with input_scaling
!   *NDEPOSITION_DIFF_FILE* CHARACTER    Nitrogen deposition  differences of driver data
!                                         - only with input_scaling_ndep
!   *REF_YEAR_PAST*         INTEGER      reference years for input scaling
!   *REF_YEAR_RECENT*       INTEGER
! 
! !compare jsbach_ctl
!   *USE_DYNVEG*            LOGICAL      true to run the dynamic vegetation
!   *USE_DISTURBANCE*       LOGICAL      true to run disturbances independent of dynamic vegetation
!   *LCC_FORCING_TYPE*      LOGICAL      Scheme for landuse change: maps, transitions, none
!   *LCC_SCHEME*            INTEGER      Pool scheme for landuse change: 1 : litter (org) 2 : anthropognic
!   *DEBUG*                 LOGICAL      true for debugging
!   *DEBUG_CCONSERVATION*   LOGICAL      true to debug carbon conservation
!   *INPUT_VERBOSE*         CHARACTER    Flags for verbose level from mo_input
!
! !compare cbalance_ctl
!   *READ_CPOOLS*           LOGICAL      read inital values for carbon pools from file
!
!   *CPOOL_FILE_NAME*       CHARACTER    File name and path for initial data of carbon pools
!   *WITH_NITROGEN*         LOGICAL      If .TRUE. nitrogen will be cycled in addition to carbon
!                                        Otherwise: only carbon cycle is simulated.
!   *READ_NPOOLS*           LOGICAL      Initialize nitrogen pools with data from an 
!                                        external file 
!   *NPOOL_FILE_NAME*       CHARACTER    File name and path for initial data of nitrogen pools
!   *READ_NDEPO*            LOGICAL      Initialize nitrogen depositions from an external file
!   *NDEPO_FILE_NAME*       CHARACTER    File name and path for inital data of nitrogen depositions
!   *WITH_YASSO*            LOGICAL      If .TRUE. yasso mode will be used

NAMELIST /CBALONE_CTL/ &
      DRIVER_DATA_PATH, EXPERIMENT, SPECIFICNAME, CLIMATE_YEARSTART, &
      CLIMATE_YEAREND, RUN_YEAR_FIRST, RUN_YEAR_LAST, OUT_INTERVAL, &
      INPUT_SCALING, CO2_FILE, CLIMATOLOGY_DIFF_FILE, REF_YEAR_PAST, &
      REF_YEAR_RECENT, USE_DYNVEG, USE_DISTURBANCE, LCC_FORCING_TYPE, &
      LCC_SCHEME, &
      READ_CPOOLS, CPOOL_FILE_NAME, WITH_NITROGEN, READ_NPOOLS, NPOOL_FILE_NAME, &
      READ_NDEPO, NDEPO_FILE_NAME, INPUT_SCALING_NDEP, NDEPOSITION_DIFF_FILE, &
      DEBUG, DEBUG_CCONSERVATION, INPUT_VERBOSE, WITH_YASSO
!                                                                       
!     ------------------------------------------------------------------
