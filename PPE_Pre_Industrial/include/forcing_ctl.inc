!*    NAMELIST *FORCING_CTL* To control atmospheric forcing of jsbach in
!                            stand-alone configuration
!                                                                 
!      
!      NAME                      TYPE      PURPOSE                                       
!      ----                      ----      -------                                       
!                                                                       
!  Temperature
!  -----------
!  *FORCING_TEMP_FREQU*       CHARACTER Frequency of temperature forcing data
!                                       DAILY: Input data are daily
!                                       MONTHLY: Input data are monthly 
!                                       CONST: Constant values; need to be defined by keywords 
!                                              FORCING_TEMP_CONST_TMIN and FORCING_TEMP_CONST_TMAX
!  *FORCING_TEMP_FILE*        CHARACTER Name base of netcdf file with temperature forcing data [degC]
!                                       Field names 'tmin' and 'tmax'
!                                       Not needed for FORCING_TEMP_FREQU = CONST
!  *FORCING_TEMP_CONST_TMIN*  REAL      Constant value globally used for minimum daily temperature [degC]
!                                       Only needed for FORCING_TEMP_FREQU = CONST.
!  *FORCING_TEMP_CONST_TMAX*  REAL      Constant value globally used for maximum daily temperature [degC]
!                                       Only needed for FORCING_TEMP_FREQU = CONST.
!
!  Precipitation
!  -------------
!  *FORCING_PRECIP_FREQU*     CHARACTER Frequency of precipitation forcing data
!                                       DAILY: Input data are daily
!                                       MONTHLY: Input data are monthly 
!                                       CONST: Constant value; need to be defined by keyword 
!                                              FORCING_PRECIP_CONST_PRECIP
!  *FORCING_PRECIP_FILE*      CHARACTER Name base of netcdf file with precipitation forcing data [mm/day]
!                                       Field names: 'precip'(rain + snow in mm/day) and 'fwet' (number of wet 
!                                       days per month; not needed with daily input)
!                                       Not needed for FORCING_PRECIP_FREQU = CONST.
!  *FORCING_PRECIP_CONST_PRECIP*  REAL  Constant value globally used for precipitation [mm/day]
!                                       Only needed for FORCING_PRECIP_FREQU = CONST.
!  *FORCING_PRECIP_IN_MM_PER_DAY* LOGICAL  Unit of the input precipitation data:
!                                          true: in mm/day; false: in kg/m^2/s
!
!
!  Atmospheric Humidity
!  -------------
!  *FORCING_QAIR_TYPE*        CHARACTER Type of input data used for atmospheric humidity forcing:
!                                       NONE: specific humidity is calculated from air temperature
!                                       RH: relative humidity in percent
!                                       QAIR: mean daily or monthly specific humidity [kg/kg]
!  *FORCING_QAIR_FREQU*       CHARACTER Frequency of atmospheric humidity forcing data
!                                       DAILY: Input data are daily
!                                       MONTHLY: Input data are monthly
!                                       CONST: Constant value; need to be defined by keyword
!                                              FORCING_QAIR_CONST_RH
!  *FORCING_QAIR_CONST_RH*    REAL      Constant value globally used for relative humidity of the atmosphere [%]
!                                       Only needed for FORCING_QAIR_FREQU = CONST.
!  *FORCING_QAIR_FILE*        CHARACTER Name base of netcdf file with atmospheric humidity forcing data
!                                       Field names: 'rel_humidity' or 'qair', depending on FORCING_QAIR_TYPE
!                                       Not needed with FORCING_QAIR_FREQU = CONST.
!  Shortwave Radiation
!  -------------------
!  *FORCING_SW_TYPE*          CHARACTER Type of input data used for forcing by shortwave radiation: 
!                                       CLOUD: relative cloudcover in percent
!                                       MEAN_RAD: mean daily or monthly downward shortwave radiation [W/m^2]
!  *FORCING_SW_SCHEME*        CHARACTER Scheme to generate the  the components (PAR, NIR) of shortwave forcing
!                                       ORIG: Wolfgangs original scheme as described in his dissertation
!                                       CCDAS: The scheme currently used in CCDAS
!                                       NEW: Wolfgangs new algorithm 
!                                       NEW_ORIG_CF: Wolfgangs new algorithm with original conversion factor      
!  *FORCING_SW_FREQU*         CHARACTER Frequency of shortwave forcing data
!                                       DAILY: Input data are daily
!                                       MONTHLY: Input data are monthly 
!                                       CONST: constant value; need to be defined by keyword 
!                                              FORCING_SW_CONST_SHORTWAVE or FORCING_SW_CONST_CLOUD.
!  *FORCING_SW_FILE*          CHARACTER Name base of netcdf file with shortwave radiation data
!                                       Depending on the value for FORCING_SW_TYPE one has to
!                                       distinguish between:
!                                       CLOUD: field name 'cloud'; cloud cover [percent]
!                                       MEAN_RAD: field name 'shortwave'; downward shortwave rad. [W/m^2]
!                                       Not needed with FORCING_SW_FREQU = CONST
!  *FORCING_SW_CONST_CLOUD*   REAL  Constant value globally used for cloudcover [percent]
!                                         Only needed when FORCING_SW_FREQU = CONST and 
!                                         FORCING_SW_TYPE = CLOUD
!  *FORCING_SW_CONST_SHORTWAVE*  REAL  Constant value globally used for downward shortwave radiation [W/m^2]
!                                             Only needed when FORCING_SW_FREQU = CONST and 
!                                             FORCING_SW_TYPE = MEAN_RAD
!  *FORCING_TABLE_SW_POT_FREQU*  CHARACTER  Frequency of potential shortwave radiation
!                                           Only needed for FORCING_SW_TYPE = MEAN_RAD.
!                                          DAILY: Input data are daily
!                                          MONTHLY: Input data are monthly 
!                                          CONST: constant value; needs to be defined by keyword 
!                                                 FORCING_TABLE_SW_POT_CONST
!  *FORCING_TABLE_SW_POT_FILE*  CHARACTER  Name base of netcdf file containing table of potential sw-radiation [W/m^2]
!                                          Needed only for FORCING_SW_TYPE = MEAN_RAD.
!                                          field name 'mpot'; potential (i.e. clear-sky) shortwave radiation.
!  *FORCING_TABLE_SW_POT_CONST* CHARACTER  Globally constant value for potential shortwave radiation [W/m^2]
!                                          Only needed when FORCING_SW_TYPE = MEAN_RAD and 
!                                          FORCING_TABLE_SW_POT_FREQU = CONST.
!
!  Longwave Radiation
!  ------------------
!  *FORCING_LW_TYPE*          CHARACTER  Type of forcing data for longwave forcing
!                                        CLOUD: fractional cloud cover used as input data (in percent)
!                                        MEAN_RAD: Mean daily or monthly downward longwave radiation [W/m^2]
!  *FORCING_LW_FREQU*         CHARACTER  Frequency of longwave forcing data
!                                        DAILY: Input data are daily
!                                        MONTHLY: Input data are monthly 
!                                        CONST: constant value; needs to be defined by keyword 
!                                               FORCING_LW_CONST_LONGWAVE or FORCING_LW_CONST_CLOUD
!  *FORCING_LW_FILE*          CHARACTER  Name base of netcdf file with longwave radiation data [W/m^2]
!                                        Depending on the value for FORCING_LW_TYPE one has to
!                                        distinguish between:
!                                        CLOUD: field name 'cloud'; cloud cover [percent]
!                                        MEAN_RAD: field name 'longwave'; downward longwave radiation [W/m^2]
!                                        Not needed with FORCING_LW_FREQU = CONST
!  *FORCING_LW_CONST_CLOUD*   REAL  Constant value globally used for cloud cover [percent]
!                                         Only needed when FORCING_LW_FREQU = CONST and
!                                         FORCING_LW_TYPE = CLOUD
!  *FORCING_LW_CONST_LONGWAVE* REAL  Constant value globally used for downward longwave radiation [W/m^2]
!                                           Only needed when FORCING_LW_FREQU = CONST and 
!                                           FORCING_LW_TYPE = MEAN_RAD
!
!  CO2 Forcing
!  -----------
!  *FORCING_CO2_FREQU*        CHARACTER  Frequency of CO2 forcing data
!                                        DAILY: Input data are daily
!                                        MONTHLY: Input data are monthly 
!                                        CONST: constant value; needs to be defined by keyword 
!                                               FORCING_CO2_CONST_CO2 
!  *FORCING_CO2_FILE*         CHARACTER  Name base of netcdf file with CO2-concentration data 
!                                        (field name 'CO2')
!  *FORCING_CO2_UNIT*         CHARACTER  Unit of CO2 concentration 
!                                        PPMV: Input data are in ppmv CO2 
!                                        MOL_PER_MOL: Input data are in mol(CO2)/mol(DryAir) 
!                                        KG_PER_KG: Input data are in kg(CO2)/kg(DryAir) 
!  *FORCING_CO2_CONST_CO2*    REAL       global value for CO2-concentration 
!                                        Only needed when FORCING_CO2_FREQU = CONST
!
! Wind forcing
! ------------
!  *FORCING_WIND_FREQU*       CHARACTER  Frequency of windspeed forcing data
!                                        DAILY: Input data are daily
!                                        MONTHLY: Input data are monthly 
!                                        CONST: constant value; needs to be defined by keyword 
!                                               FORCING_WIND_CONST_WSPEED
!  *FORCING_WIND_FILE*        CHARACTER  Name base of netcdf file with windspeed data [m/s] (field name 'wspeed')
!  *FORCING_WIND_CONST_WSPEED*  REAL     Constant value globally used for windspeed [m/s]
!                                        Only needed when FORCING_WIND_FREQU = CONST
!
NAMELIST /FORCING_CTL/ &
      FORCING_TEMP_FREQU, FORCING_TEMP_FILE, FORCING_TEMP_CONST_TMIN, &
      FORCING_TEMP_CONST_TMAX, FORCING_PRECIP_FREQU, FORCING_PRECIP_FILE, &
      FORCING_PRECIP_CONST_PRECIP, FORCING_PRECIP_IN_MM_PER_DAY, &
      FORCING_SW_TYPE, FORCING_SW_SCHEME, FORCING_SW_FREQU, FORCING_SW_FILE, &
      FORCING_QAIR_CONST_RH,FORCING_QAIR_TYPE, FORCING_QAIR_FREQU, FORCING_QAIR_FILE, &
      FORCING_SW_CONST_CLOUD, FORCING_SW_CONST_SHORTWAVE, &
      FORCING_TABLE_SW_POT_FREQU, FORCING_TABLE_SW_POT_FILE, &
      FORCING_TABLE_SW_POT_CONST, FORCING_LW_TYPE, FORCING_LW_FREQU, &
      FORCING_LW_FILE, FORCING_LW_CONST_CLOUD, FORCING_LW_CONST_LONGWAVE, &
      FORCING_CO2_FREQU, FORCING_CO2_FILE, FORCING_CO2_CONST_CO2,FORCING_CO2_UNIT, &
      FORCING_WIND_FREQU, FORCING_WIND_FILE, FORCING_WIND_CONST_WSPEED
!                                                                       
!     ------------------------------------------------------------------
