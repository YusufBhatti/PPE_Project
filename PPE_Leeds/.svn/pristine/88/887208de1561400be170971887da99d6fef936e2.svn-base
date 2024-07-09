!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! This module contains substantial contributions by the Max Planck Institute for Biogeochemistry, Jena.
!!
!! ------------------------------------------------------------------------------------------------------
!!
!! This module belongs to the standalone version of JSBACH
!!
MODULE mo_jsbalone_forcing
  ! 
  !   Description: This module contains code to handle the forcing of the standalone JSBACH
  ! 
  ! Current Code Owner:  jsbach_admin
  !   Christian Reick, Thomas Raddatz (prior to 2010)
  !   Gregor Schuermann, Nuno Carvalhais, Soenke Zaehle (2010-)
  ! 
  ! History:
  !
  ! Version       Date                Comment
  ! -------       ------              -------
  ! 0.1           ????/??/??          Org. Code admin
  ! 0.2           2011/03/31          implementation of model forcing at model time step
  ! 0.3           2011/07/27          updates to code to correct errors in shortwave radiation calculation,  
  !                                   correction in the calculation of downward longwave radiation,
  !                                   implementation of atmospheric humidity forcing,
  !                                   implementation of a corresponding weather generator for atm. humidity
  !                                   removal of obsolete code w.r.t. the implicit hydrological coupling,
  !                                   the code fraction necessary to calculate actual / potential evapotransration 
  !                                   for unforced atm humidity is now part of the driving module mo_jsbalone.
  !                                   documentation and completion of forcing stream 
  ! 0.31          2012/04/03          correction of unit for CO2 concentration, which has to be in kg kg-1 at the
  !                                   end of this routine,SZ
  !
  ! Unresolved issues:
  !   CHR, 25.4.03: This module is incomplete in the following respects:
  !                 (1) The restriction to a rectangular section of the grid is not implemented. Therefore the program will
  !                     probably crash if it is not run on the full grid.
  !                 SZ, to my knowledge this is not the case, the code runs on single points and 
  !                      reduced domains for the cases tested
  !
  !
  !  GJS: requires input of air pressure? Potentially an asset if re-analysis / Echam data are used
  !  SZ:  different time step settings for forcing fields are allowed and should work, but may lead to inconsistent results
  !       The code would be cleaner if only timestep or daily forcing would be allowed for all variables, but this is left
  !       for now to help debugging and development of code
  !  SZ:  the potential radiation tables need to be removed or integrated into ForcingFields (as done now for diagnostics).
  !       To my understanding the radiation code based on fpar does not produce realistic radiation sums and hence should be deleted
  !       be deleted rather reworked

  USE mo_jsbach,        ONLY: options_type, debug
  USE mo_jsbach_constants, ONLY: molarMassCO2_kg, molarMassDryAir_kg, pi, p_sealevel, GasConstantDryAir, Gravity, Tmelt, &
                                 StefanBoltzmann, eps, solar_const, euler
  USE mo_jsbach_grid,   ONLY: grid_type, domain_type
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_util_string,   ONLY: int2string, real2string
  USE mo_mpi,           ONLY: p_parallel_io, p_parallel, p_io, p_bcast
  USE mo_echam_convect_tables, ONLY: tlucua
  USE mo_netcdf,        ONLY: nf_max_name, nf_open, NF_NOERR, nf_nowrite, nf_inq_dimid, nf_inq_dimlen, &
                              nf_get_var_double, nf_get_vara_double, nf_close, nf_inq_varid
  USE mo_kind,          ONLY: dp
  USE mo_time_control,  ONLY: lresume, lstart

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART =======================================================================================================

  PUBLIC :: init_forcing    ! Subroutine that reads the forcing namelist and allocates memory. 
  PUBLIC :: read_forcing    ! Subroutine that reads the forcing data. Has to be called for each (annual) forcing file.
  PUBLIC :: update_forcing  ! Subroutine that updates the forcing data. Has to be called at the beginning of each timestep.
  PUBLIC :: finish_forcing  ! Subroutine that deallocates memory that was allocated by init_forcing. To be used only when ending

  PUBLIC :: forcing_type    ! This is the central structure of the module on which all subroutines of the program operate.

  TYPE forcing_type
     REAL(dp), POINTER :: air_temp(:)        !! Air temperature at bottom (depending on local elevation) [Celsius]
     REAL(dp), POINTER :: precip_rain(:)     !! Liquid part of precipitation [kg/m^2/s]
     REAL(dp), POINTER :: precip_snow(:)     !! Snow part of precipitation [kg/m^2/s]
     REAL(dp), POINTER :: air_pressure(:)    !! air pressure at bottom (depending on local elevation) [N/m^2]
     REAL(dp), POINTER :: vapor_pressure(:)  !! vapor pressure [N/m^2]
     REAL(dp), POINTER :: spec_humidity(:)   !! specific humidity (from [0,1])
     REAL(dp), POINTER :: rad_sw_down(:)     !! total solar shortwave radiation flux: direct+diffuse, 280-3000 nm, (i.e. PAR+NIR),..
                                             !! ... that arrives at the bottom (includes cloud shading and scattering) [W/m^2]
     REAL(dp), POINTER :: rad_UV_down(:)     !! solar radiation flux (direct+diffuse) from the UV band 280-400 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
                                             !! Attention!!! This is not an accurate estimate of UV. It's simply a residual:
                                             !! rad_UV_down = rad_sw_down - rad_PAR_down - rad_NIR_down
     REAL(dp), POINTER :: rad_PAR_down(:)    !! solar radiation flux (direct+diffuse) from the visible band 400-700 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
     REAL(dp), POINTER :: rad_NIR_down(:)    !! solar radiation flux (direct+diffuse) from the near infrared band 700 - 3000 nm ...
                                             !! .. arriving at the bottom  (includes cloud shading and scattering) [W/m^2]
     REAL(dp), POINTER :: frac_PAR_diffuse(:)!! Fraction from  rad_PAR_down() that comes down as diffuse radiation (values in [0,1])
     REAL(dp), POINTER :: rad_sw_down_pot(:) !! potential shortwave downward radiation at the surface [W/m^2] 
     REAL(dp), POINTER :: rad_lw_down(:)     !! longwave downward flux [W/m^2] from thermal radiation of atmosphere and clouds
     REAL(dp), POINTER :: wind_speed(:)      !! wind speed of lowest atmosphere model layer [m/s]
     REAL(dp), POINTER :: wind_speed10(:)    !! 10m wind speed [m/s]
     REAL(dp), POINTER :: CO2_concentr(:)    !! CO2 concentration at bottom [kg(CO2)/kg(dry air)]
     REAL(dp), POINTER :: Tmin(:)            !! minimum temperature
     REAL(dp), POINTER :: Tmax(:)            !! maximum temperature
  END TYPE forcing_type

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ==============================================================================

  PRIVATE ! Make ALL following objects private

  ! === Private declarations =======================================================================================================

  ! --- time control
  LOGICAL ::  new_year_forcing     ! in contrast to new_year, new_month and new_day (stepon_jsbach), these switches are defined
  LOGICAL ::  new_month_forcing    !   comparing current_date with next_date. With stepwise forcing the dates in the forcing  
  LOGICAL ::  new_day_forcing      !   file shall correspond to the dates of stepwise output.

  ! --- parameters -----

  REAL(dp), PARAMETER :: gamma  = 6.0E-3_dp            ! fixed temperature gradient of atmosphere [K/m] (see Knorr, p. 31). Note that
                                                       ! for the standard ICAO-atmosphere this value is different, namely 6.5E-3 !!!
  REAL(dp), PARAMETER :: solar_const_vis = solar_const * 0.44_dp  ! solar constant corrected for visible rad. [W/m^2]
  REAL(dp), PARAMETER :: conv_prec_day   =  1._dp/86400._dp   ! Conversion factor for precipitation in [mm/day] to [kg/m^2/s]:
                                                              ! 1 mm/day is equivalent to 1 kg /86400 s / m^2
  REAL(dp), PARAMETER :: negligible = 1.0E-15_dp              ! a small number used to detect and prevent eventual division by zero

  ! ---- enumerations ---------------

  INTEGER, PARAMETER :: NONE_ = 0

  INTEGER,PARAMETER :: MONTHLY_ = 1
  INTEGER,PARAMETER :: DAILY_   = 2
  INTEGER,PARAMETER :: CONST_   = 3
  INTEGER,PARAMETER :: TIMESTEP_   = 4
  INTEGER,PARAMETER :: GHG_SCENARIO_  = 5
  ! enums for options controlling atmospheric humidity input data 
  INTEGER,PARAMETER :: RH_ = 1           ! relative humidity as forcing data 
  INTEGER,PARAMETER :: QAIR_ = 2         ! specific humidity as forcing data

  ! enums for options controlling radiation input data (used for longwave and shortwave)
  INTEGER,PARAMETER :: CLOUD_ = 1            ! cloudcover as forcing data 
  INTEGER,PARAMETER :: MEAN_RAD_ = 2         ! mean daily or monthly incoming radiation as input data

  ! enums for options controlling the shortwave-radiation algorithms used
  INTEGER,PARAMETER :: SW_SCHEME_ORIGINAL_    = 1 ! Use shortwave_from_fpar_orig()  (scheme from Wolfgangs Diss.)
  INTEGER,PARAMETER :: SW_SCHEME_CCDAS_       = 2 ! Use shortwave_from_fpar_ccdas() (scheme from CCDAS)
  INTEGER,PARAMETER :: SW_SCHEME_NEW_         = 3 ! Use shortwave_from_fsw()  (Wolfgangs new scheme)
  INTEGER,PARAMETER :: SW_SCHEME_NEW_ORIG_CF_ = 4 ! Use shortwave_from_fsw()  (Wolfgangs new scheme with original conversion factor)
  INTEGER,PARAMETER :: SW_SCHEME_DIRECT_        = 5 ! Use shortwave_from_direct_shortwave() (scheme from shortwave record)

  ! --- Namelist parameters
  CHARACTER(len=128)     :: forcing_temp_frequ           !! Frequency of temperature forcing data (DAILY/MONTHLY/CONST/TIMESTEP)
  CHARACTER(nf_max_name) :: forcing_temp_file            !! File with temperature forcing data [degC]
  REAL(dp)               :: forcing_temp_const_tmin      !! Constant value for minimum daily temperature [degC]
  REAL(dp)               :: forcing_temp_const_tmax      !! Constant value for maximum daily temperature [degC]
  CHARACTER(len=128)     :: forcing_precip_frequ         !! Frequency of precipitation forcing data (DAILY/MONTHLY/CONST/TIMESTEP)
  CHARACTER(nf_max_name) :: forcing_precip_file          !! File with precipitation forcing data (unit: mm/day or kg/m^2/s,
                                                         !! depending on 'forcing_precip_in_mm_per_day')
  REAL(dp)               :: forcing_precip_const_precip  !! Constant value for precipitation [mm/day]
  LOGICAL                :: forcing_precip_in_mm_per_day !! Precipitation input in mm/day (true) or kg/m^2/s (false)
  CHARACTER(len=128)     :: forcing_sw_type              !! Type of input data used for forcing by shortwave radiation
  CHARACTER(len=128)     :: forcing_sw_scheme            !! Scheme to generate the  the components (PAR, NIR) of sw forcing
  CHARACTER(len=128)     :: forcing_sw_frequ             !! Frequency of shortwave forcing data (DAILY/MONTHLY/CONST/TIMESTEP)
  CHARACTER(nf_max_name) :: forcing_sw_file              !! File with shortwave radiation data
  REAL(dp)               :: forcing_sw_const_cloud       !! Constant value for cloudcover [percent]
  REAL(dp)               :: forcing_sw_const_shortwave   !! Constant value for downward shortwave radiation [W/m^2]
  CHARACTER(len=128)     :: forcing_table_sw_pot_frequ   !! Frequency of potential sw radiation (DAILY/MONTHLY/CONST/TIMESTEP)
  CHARACTER(nf_max_name) :: forcing_table_sw_pot_file    !! File with potential short wave radiation data [W/m^2]
  REAL(dp)               :: forcing_table_sw_pot_const   !! Constant value for potential shortwave radiation [W/m^2]
  CHARACTER(len=128)     :: forcing_lw_type              !! Type of forcing data for longwave forcing
  CHARACTER(len=128)     :: forcing_lw_frequ             !! Frequency of longwave forcing data (DAILY/MONTHLY/CONST/TIMESTEP)
  CHARACTER(nf_max_name) :: forcing_lw_file              !! file with longwave radiation data [W/m^2]
  REAL(dp)               :: forcing_lw_const_cloud       !! Constant value for cloud cover [percent]
  REAL(dp)               :: forcing_lw_const_longwave    !! Constant value for downward longwave radiation [W/m^2]
  CHARACTER(len=128)     :: forcing_co2_frequ            !! Frequency of CO2 forcing data (DAILY/MONTHLY/CONST)
  CHARACTER(nf_max_name) :: forcing_co2_file             !! File with CO2-concentration data (unit: mol(CO2)/mol(air) or 
                                                         !! kg(CO2)/kg(air) or ppmv, depending on 'forcing_co2_unit'
  REAL(dp)               :: forcing_co2_const_co2        !! Constant value for CO2-concentration (unit: mol(CO2)/mol(air) or 
                                                         !! kg(CO2)/kg(air) or ppmv, depending on 'forcing_co2_unit')
  CHARACTER(len=128)     :: forcing_co2_unit             !! Unit of CO2 input (PPMV/MOL_PER_MOL/KG_PER_KG)
  CHARACTER(len=128)     :: forcing_wind_frequ           !! Frequency of windspeed forcing data (DAILY/MONTHLY/CONST/TIMESTEP)
  CHARACTER(nf_max_name) :: forcing_wind_file            !! File with windspeed data [m/s] 
  REAL(dp)               :: forcing_wind_const_wspeed    !! Constant value for windspeed [m/s]
  CHARACTER(len=128)     :: forcing_qair_type            !! Type of input data used for forcing by atm. humidity
  CHARACTER(len=128)     :: forcing_qair_frequ           !! Frequency of qair forcing data (DAILY/MONTHLY/CONST)
  CHARACTER(nf_max_name) :: forcing_qair_file            !! File with qair data 
  REAL(dp)               :: forcing_qair_const_rh        !! Constant value for relative humidity [%]

  ! --- Structures for handling input data for forcing -----------------------------------------------------------------------------

  !! The following structures collects all information on options relating to forcing (read in from the JSBACH run defiition file)

  TYPE forcing_options_type
     INTEGER :: type_of_qair_forcing           !! Possible values: see enumerations from above
     INTEGER :: type_of_shortwave_forcing      !! Possible values: see enumerations from above
     INTEGER :: type_of_shortwave_scheme       !! Possible values: see enumerations from above
     INTEGER :: type_of_longwave_forcing       !! Possible values: see enumerations from above
     LOGICAL :: precip_in_mm_per_day           !! Specifies units of precipitation input data:
                                               !! true: mm/day; false: kg/(m^2*s)
     REAL(dp) :: conv_CO2_2_MassRatio          !! Factor for the conversion of CO2 input data to kg(CO2)/kg(dry air)
  END TYPE forcing_options_type
  TYPE(forcing_options_type),SAVE :: forcing_options

  !! The following structure collects all information on a single input variable from a netcdf file
  TYPE input_variable_type 
     LOGICAL                  :: initialized = .FALSE.
     CHARACTER(nf_max_name)   :: variable_name="NONE"!! name of variable (name in netcdf-data file)
     CHARACTER(nf_max_name)   :: file_name = "NONE"  !! File in which variable is found
     INTEGER                  :: file_ID = 0         !! NetCDF file ID obtained from nf_open()
     INTEGER                  :: var_ID = 0          !! NetCDF variable ID obtained from nf_inq_varid ()
     INTEGER                  :: frequency =0        !! This number encodes the frequency of the input data:
                                                     !! 1: MONTHLY; 2: DAILY, 3: CONSTANT 4: TIMESTEP, 5: GHG_SCENARIO (CO2 only)
     REAL(dp),POINTER         :: data_daily(:)       !! pointer to daily data (if frequency=2 this points to the input data,
                                                     !! .. whereas for frequency=1 this points to data generated from monthly
                                                     !! .. data each time a day changes) (dim: domain%nland)
     REAL(dp),POINTER         :: data_monthly(:,:)   !! If input data are monthly, this points to those input data; otherwise
                                                     !! .. this pointer is not needed. (dim: domain%nland x 12)
     REAL(dp),POINTER         :: data_timestep(:,:)  !! pointer to timestep data 
                                                     !! .. (dim: domain%nland x #_of_timesteps)
     INTEGER                  :: length_time_series  !! The number of time steps found in the data, i.e. the number of months for ..
                                                     !! .. monthly data (only 12 allowed!) and the number of days for daily data.
     REAL(dp), POINTER        :: timevals(:)         !! The values of the time variable. Must be in YYYYMM format for monthly
                                                     !! .. data and YYYYMMDD for daily data
  END TYPE input_variable_type

  !! The following structure collects all POSSIBLE input data, i.e. some of the are usually not used. Which are used
  !! is controlled by options read in from the JSBACH runtime definition file.
  TYPE forcing_input_type
     TYPE(input_variable_type) :: air_temp_min  !! daily or monthly minimum of air temperature at surface [Celsius]
     TYPE(input_variable_type) :: air_temp_max  !! daily or monthly maximum of air temperature at surface [Celsius]
     TYPE(input_variable_type) :: air_temp      !! air temperature [Celsius]
     TYPE(input_variable_type) :: precipitation !! daily or monthly precipitation (including snow); units may be   
                                                !! [mm/day] or [kg/m^2/s] (see option  "precip_in_mm_per_day")
     TYPE(input_variable_type) :: frac_wet_days !! fraction of wet days per month (only needed for monthly precip. data!)
     TYPE(input_variable_type) :: cloud_cover   !! fraction of sky covered by clouds (from [0,1])
     TYPE(input_variable_type) :: shortwave     !! incident shortwave radiation [W/m^2]
     TYPE(input_variable_type) :: longwave      !! incident longwave radiation [W/m^2]
     TYPE(input_variable_type) :: CO2_concentr  !! CO2-concentration [mol(CO2)/mol(air)], [kg(CO2)/kg(dry Air)] or [ppmv] 
     TYPE(input_variable_type) :: wind_speed    !! wind speed [m/s]
     TYPE(input_variable_type) :: rel_humidity  !! relative humidity [%]
     TYPE(input_variable_type) :: qair          !! specific humidity [kg/kg]
  END TYPE forcing_input_type

  TYPE(forcing_input_type),SAVE :: forcing_input

  TYPE table_type
     TYPE(input_variable_type) :: shortwave_pot !! potential (i.e. clear sky) shortwave radiation incident to surface [W/m^2]
  END TYPE table_type

  TYPE(table_type),SAVE :: tables

  ! --- others ----------------------

  REAL(dp), ALLOCATABLE,TARGET     :: hlp_field3D_global(:,:,:)  ! temporary memory for reading global 3D-domain forcing data
  REAL(dp), ALLOCATABLE,TARGET     :: hlp_field2D_global(:,:) ! temporary memory for reading in global daily forcing data
  REAL(dp), ALLOCATABLE            :: hlp_field2D_local(:,:)  ! temporary memory for a 2D-domain on each processor

  ! --- output streams --------------

  ! All diagnostic output fields are collected in the following structure:

  TYPE forcing_diag_type
     TYPE(forcing_type) :: forcingFields       ! All forcing fields (mean values)
  END TYPE forcing_diag_type

  TYPE(forcing_diag_type),SAVE :: forcing_diag !! Collects all diagnosic fields, in particular all forcing fields

CONTAINS 

  ! --- init_forcing() -------------------------------------------------------------------------------------------------------------
  !
  ! This routine:
  !
  !   - reads in all forcing data (monthly means), or otherwise: 
  !   - checks (partly) consistence of forcing data with information from grid file
  !   - allocates memory for the forcing fields
  !
  SUBROUTINE init_forcing(grid, domain, global_options, forcingFields)

    USE mo_time_control, ONLY: delta_time
    USE mo_util_string, ONLY: tolower
    USE mo_namelist, ONLY: position_nml, POSITIONED
    USE mo_jsbach, ONLY: nml_unit
    USE mo_io_units, ONLY: nout
    USE mo_input_strings, ONLY: GetStringIndex
    USE mo_input,         ONLY: input_var_list, input_opt_list, &
                                INPUT_VAR_DIM_MIS_SPREAD, INPUT_VAR_VALID_MIDPOINT, &
                                InputAttrGet, InputVarAdd, InputOptNew

    TYPE(grid_type),    INTENT(in)    :: grid
    TYPE(domain_type),  INTENT(in)    :: domain
    TYPE(options_type), INTENT(inout) :: global_options  !! Part of the global options is read in here
    TYPE(forcing_type), INTENT(inout) :: forcingFields   !! This routine allocates memory for this structure 


    ! --- Other local parameters
    INTEGER  :: status
    INTEGER  :: read_status, f_unit
    INTEGER  :: time_steps_per_day ! number of forcing timesteps per day
    TYPE (input_var_list), POINTER :: CO2_var
    TYPE (input_opt_list), POINTER :: opt
    CHARACTER (len=64) :: CO2_unit

    INCLUDE 'forcing_ctl.inc'

    ! --- Nullify all pointer variables of forcing ---------------------------------------------------------------------------------

    NULLIFY(forcingFields%air_temp)
    NULLIFY(forcingFields%precip_rain)
    NULLIFY(forcingFields%precip_snow)
    NULLIFY(forcingFields%air_pressure)
    NULLIFY(forcingFields%vapor_pressure)
    NULLIFY(forcingFields%spec_humidity)
    NULLIFY(forcingFields%rad_sw_down)
    NULLIFY(forcingFields%rad_sw_down_pot)
    NULLIFY(forcingFields%rad_UV_down)
    NULLIFY(forcingFields%rad_PAR_down)
    NULLIFY(forcingFields%rad_NIR_down)
    NULLIFY(forcingFields%frac_PAR_diffuse)
    NULLIFY(forcingFields%rad_lw_down)
    NULLIFY(forcingFields%wind_speed)
    NULLIFY(forcingFields%CO2_concentr)
    NULLIFY(forcingFields%Tmin)
    NULLIFY(forcingFields%Tmax)

    ! --- Read namelist forcing_ctl ------------------------------------------------------------------------------------------------

    IF (p_parallel_io) THEN

       ! define default values

       forcing_temp_frequ = 'DAILY'
       forcing_temp_file = 'climate'
       forcing_temp_const_tmin = HUGE(0.0_dp)
       forcing_temp_const_tmax = HUGE(0.0_dp)
       forcing_precip_frequ = 'DAILY'
       forcing_precip_file = 'climate.nc'
       forcing_precip_const_precip = HUGE(0.0_dp)
       forcing_precip_in_mm_per_day = .TRUE.
       forcing_sw_type = 'MEAN_RAD'
       forcing_sw_scheme = 'NONE'
       forcing_sw_frequ = 'DAILY'
       forcing_sw_file = 'climate.nc'
       forcing_sw_const_cloud = HUGE(0.0_dp)
       forcing_sw_const_shortwave = HUGE(0.0_dp)
       forcing_table_sw_pot_frequ = 'DAILY'
       forcing_table_sw_pot_file = 'climate.nc'
       forcing_table_sw_pot_const = HUGE(0.0_dp)
       forcing_lw_type = 'MEAN_RAD'
       forcing_lw_frequ = 'DAILY'
       forcing_lw_file = 'climate.nc'
       forcing_lw_const_cloud = HUGE(0.0_dp)
       forcing_lw_const_longwave = HUGE(0.0_dp)
       forcing_co2_frequ = 'DAILY'
       forcing_co2_file = 'climate.nc'
       forcing_co2_const_co2 = 3.65e-4_dp 
       forcing_co2_unit = '' 
       forcing_wind_frequ = 'DAILY'
       forcing_wind_file = 'climate.nc'
       forcing_wind_const_wspeed = HUGE(0.0_dp)
       forcing_qair_type = 'NONE'
       forcing_qair_frequ = 'DAILY'
       forcing_qair_file = 'climate.nc'
       forcing_qair_const_rh = 100._dp 
       
       f_unit = position_nml ('FORCING_CTL', nml_unit, status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (f_unit, forcing_ctl)
          CALL message('init_forcing', 'Namelist FORCING_CTL: ')
          WRITE(nout, forcing_ctl)
       END SELECT

       forcing_options%precip_in_mm_per_day = forcing_precip_in_mm_per_day
       IF(forcing_options%precip_in_mm_per_day) THEN
          CALL message("init_forcing",'Precipitation forcing is in mm/day')
       ELSE
          CALL message("init_forcing",'Precipitation forcing is in kg/m^2/s')
       END IF

       SELECT CASE(TRIM(tolower(forcing_co2_unit)))
       CASE("ppmv")
           CALL message("init_forcing",'CO2 forcing is in ppmv')
           forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg * 0.000001_dp
       CASE("mol_per_mol")
           CALL message("init_forcing",'CO2 forcing is in mol(CO2)/mol(Dry Air)')
           forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg 
       CASE("kg_per_kg")
           CALL message("init_forcing",'CO2 forcing is in kg(CO2)/kg(Dry Air)')
           forcing_options%conv_CO2_2_MassRatio = 1.0_dp
       CASE default
           CALL message("init_forcing",'WARNING: "FORCING_CO2_IN_UNIT" missing in namelist forcing_ctl.')
           CALL message("init_forcing",'Recovering by assuming that the unit is mol(CO2)/mol(DryAir). Verify your settings!')
           forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg 
       END SELECT
       forcing_co2_const_co2 = forcing_co2_const_co2 * forcing_options%conv_CO2_2_MassRatio

       SELECT CASE(TRIM(tolower(forcing_qair_type)))
       CASE("rh")
          forcing_options%type_of_qair_forcing = RH_
       CASE("qair")
          forcing_options%type_of_qair_forcing = QAIR_
       CASE("none")
          forcing_options%type_of_qair_forcing = NONE_
          CALL message("init_forcing",'WARNING: QAIR will be diagnosed from Air temperature')
       CASE default
          forcing_options%type_of_qair_forcing = NONE_
          CALL message("init_forcing",'"FORCING_QAIR_TYPE" missing in namelist forcing_ctl.')
          CALL message("init_forcing",'WARNING: QAIR will be diagnosed from Air temperature')
       END SELECT
       CALL message("init_forcing","Type of air humidity forcing: "//TRIM(forcing_qair_type))

       SELECT CASE(TRIM(tolower(forcing_sw_type)))
       CASE("cloud")
          forcing_options%type_of_shortwave_forcing = CLOUD_
       CASE("mean_rad")
          forcing_options%type_of_shortwave_forcing = MEAN_RAD_
       CASE("NONE")
          CALL finish("init_forcing",'"FORCING_SW_TYPE" missing in namelist forcing_ctl')
       CASE default
          CALL finish("init_forcing",'Invalid value '//TRIM(forcing_sw_type)// &
               ' "FORCING_SW_TYPE" in namelist forcing_ctl.')
       END SELECT
       CALL message("init_forcing","Type of shortwave radiation forcing: "//TRIM(forcing_sw_type))

       SELECT CASE(TRIM(tolower(forcing_sw_scheme)))
       CASE("orig")
          forcing_options%type_of_shortwave_scheme= SW_SCHEME_ORIGINAL_
       CASE("ccdas")
          forcing_options%type_of_shortwave_scheme = SW_SCHEME_CCDAS_
       CASE("new")
          IF(forcing_options%type_of_shortwave_forcing == CLOUD_) THEN
             CALL finish("init_forcing",'Shortwave forcing scheme "NEW" not allowed for input type "CLOUD"')
          END IF
          forcing_options%type_of_shortwave_scheme = SW_SCHEME_NEW_
       CASE("new_orig_cf")
          IF(forcing_options%type_of_shortwave_forcing == CLOUD_) THEN
             CALL finish("init_forcing",'Shortwave forcing scheme "NEW_ORIG_CF" not allowed for input type "CLOUD"')
          END IF
          forcing_options%type_of_shortwave_scheme = SW_SCHEME_NEW_ORIG_CF_
       CASE("direct")
          forcing_options%type_of_shortwave_scheme = SW_SCHEME_DIRECT_
       CASE("NONE")
          CALL finish("init_forcing",'"FORCING_SW_SCHEME" missing in namelist.')
       CASE default
          CALL finish("init_forcing",'Invalid value '//TRIM(forcing_sw_scheme)//&
               ' of keyword "FORCING_SW_SCHEME" in namelist.')
       END SELECT
       CALL message("init_forcing","Radiation forcing sheme: "//TRIM(forcing_sw_scheme))

       IF(forcing_sw_frequ /= forcing_table_sw_pot_frequ ) THEN
          CALL finish("init_forcing","Invalid combination of shortwave ("//TRIM(forcing_sw_frequ)//&
               ") and mpot forcing frequency ("//TRIM(forcing_table_sw_pot_frequ)//"). These need to be the same.")
       ENDIF

       SELECT CASE(TRIM(tolower(forcing_lw_type)))
       CASE("cloud")
          forcing_options%type_of_longwave_forcing = CLOUD_
       CASE("mean_rad")
          forcing_options%type_of_longwave_forcing = MEAN_RAD_
       CASE("NONE")
          CALL finish("init_forcing",'"FORCING_LW_TYPE" missing in namelist.')
       CASE default
          CALL finish("init_forcing",'Invalid value '//TRIM(forcing_lw_type)//&
               ' of keyword "FORCING_LW_TYPE" in namlist.')
       END SELECT
       CALL message("init_forcing","Type of longwave radiation forcing: "//TRIM(forcing_lw_type))

       global_options%HeightTemperature = 2._dp
       global_options%HeightWind = -1._dp        ! -1: default to ECHAM lowest layer
       global_options%HeightHumidity = 10._dp

    ENDIF

    IF (p_parallel) THEN
       CALL p_bcast(global_options%HeightTemperature, p_io)
       CALL p_bcast(global_options%HeightWind,        p_io)
       CALL p_bcast(global_options%HeightHumidity,    p_io)
       CALL p_bcast(forcing_options%precip_in_mm_per_day, p_io)
       CALL p_bcast(forcing_options%type_of_qair_forcing, p_io)
       CALL p_bcast(forcing_options%type_of_shortwave_forcing, p_io)
       CALL p_bcast(forcing_options%type_of_shortwave_scheme, p_io)
       CALL p_bcast(forcing_options%type_of_longwave_forcing, p_io)
       CALL p_bcast(forcing_options%conv_CO2_2_MassRatio, p_io)
       CALL p_bcast(forcing_temp_const_tmin, p_io)
       CALL p_bcast(forcing_temp_const_tmax, p_io)
       CALL p_bcast(forcing_precip_const_precip, p_io)
       CALL p_bcast(forcing_qair_const_rh, p_io)
       CALL p_bcast(forcing_sw_const_cloud, p_io)
       CALL p_bcast(forcing_sw_const_shortwave, p_io)
       CALL p_bcast(forcing_table_sw_pot_const, p_io)
       CALL p_bcast(forcing_lw_const_cloud, p_io)
       CALL p_bcast(forcing_lw_const_longwave, p_io)
       CALL p_bcast(forcing_co2_const_co2, p_io)
       CALL p_bcast(forcing_wind_const_wspeed, p_io)
       CALL p_bcast(forcing_temp_frequ,p_io)
       CALL p_bcast(forcing_precip_frequ,p_io)
       CALL p_bcast(forcing_qair_frequ,p_io)
       CALL p_bcast(forcing_sw_frequ,p_io)
       CALL p_bcast(forcing_table_sw_pot_frequ,p_io)
       CALL p_bcast(forcing_lw_frequ,p_io)
       CALL p_bcast(forcing_co2_frequ,p_io)
       CALL p_bcast(forcing_wind_frequ,p_io)
    ENDIF

    ! --- initialization of constants ----------------------------------------------------------------------------------------------

    ! --- allocation of private fields that are used during the whole run, not only during initialization --------------------------
    IF(p_parallel_io) THEN                                                                                             !
      IF((forcing_temp_frequ == "TIMESTEP") .OR. & 
        (forcing_precip_frequ == "TIMESTEP") .OR. & 
        (forcing_qair_frequ == "TIMESTEP") .OR. &
        (forcing_sw_frequ == "TIMESTEP") .OR. & 
        (forcing_table_sw_pot_frequ == "TIMESTEP") .OR. & 
        (forcing_lw_frequ == "TIMESTEP") .OR. & 
        (forcing_co2_frequ == "TIMESTEP") .OR. & 
        (forcing_wind_frequ == "TIMESTEP")) THEN

          time_steps_per_day = 86400/int(delta_time)
          ALLOCATE(hlp_field3D_global(1:grid%nlon,1:grid%nlat,1:time_steps_per_day),STAT=status)  
          IF(status /= 0) CALL finish("init_forcing()", "ERROR: Field hlp_field3D_global() could not be allocated.")
          IF(debug) CALL message('init_forcing','hlp_field3D_global created')
      ENDIF  
  
      ALLOCATE(hlp_field2D_global(1:grid%nlon,1:grid%nlat),STAT=status)  
      IF(status /= 0) CALL finish("init_forcing()", "ERROR: Field hlp_field2D_global() could not be allocated.")
      IF(debug) CALL message ('init_forcing','hlp_field2D_global created')
 
    END IF

    ALLOCATE(hlp_field2D_local(domain%ndim,domain%nblocks),STAT=status) ! memory for a 2D-domain on each processor     
    IF(status /= 0) CALL finish("init_forcing()", "ERROR: Field hlp_field2D_local() could not be allocated.")

    !-- allocation of forcing fields

    ALLOCATE(forcingFields%air_temp(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%air_temp() could not be allocated.")

    IF(forcing_temp_frequ /= "TIMESTEP") THEN
      CALL init_variable_access(domain, grid, "FORCING_TEMP", "tmin", forcing_temp_frequ, forcing_temp_file, &
        forcing_temp_const_tmin, forcing_input%air_temp_min)
      CALL init_variable_access(domain, grid, "FORCING_TEMP", "tmax", forcing_temp_frequ, forcing_temp_file, &
        forcing_temp_const_tmax, forcing_input%air_temp_max)
        forcingFields%Tmin => forcing_input%air_temp_min%data_daily(:)
        forcingFields%Tmax => forcing_input%air_temp_max%data_daily(:)
    ELSE
      CALL init_variable_access(domain, grid, "FORCING_TEMP", "air_temp", forcing_temp_frequ, forcing_temp_file, &
         forcing_temp_const_tmin, forcing_input%air_temp)
      ALLOCATE(forcingFields%Tmin(domain%nland),STAT=status)
      IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%tmin() could not be allocated.")
      ALLOCATE(forcingFields%Tmax(domain%nland),STAT=status)
      IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%tmax() could not be allocated.")
    END IF

    ! --- initialize precipitation forcing
    ALLOCATE(forcingFields%precip_rain(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%precip_rain() could not be allocated.")
    ALLOCATE(forcingFields%precip_snow(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%precip_snow() could not be allocated.")

    ALLOCATE(forcingFields%vapor_pressure(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%vapor_pressure() could not be allocated.")
    ALLOCATE(forcingFields%spec_humidity(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%spec_humidity() could not be allocated.")

    ALLOCATE(forcingFields%rad_UV_down(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%rad_UV_down() could not be allocated.")
    ALLOCATE(forcingFields%rad_PAR_down(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%rad_PAR_down() could not be allocated.")
    ALLOCATE(forcingFields%rad_NIR_down(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%rad_NIR_down() could not be allocated.")
    ALLOCATE(forcingFields%frac_PAR_diffuse(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%frac_PAR_diffuse() could not be allocated.")
    ALLOCATE(forcingFields%rad_sw_down(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%rad_sw_down() could not be allocated.")
    ALLOCATE(forcingFields%rad_sw_down_pot(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%rad_sw_down_pot() could not be allocated.")

    IF (forcing_co2_frequ == "TIMESTEP") THEN
        ALLOCATE(forcingFields%CO2_concentr(domain%nland),STAT=status)
        IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%CO2_concentr() could not be allocated.")
    ELSEIF (TRIM(tolower(forcing_co2_frequ)) == "ghg_scenario" ) THEN
      CALL p_bcast(forcing_co2_file,p_io)
      CALL InputVarAdd("landpoints",forcingFields%CO2_concentr,"CO2",file_name=TRIM(forcing_co2_file), &
                       def_val=forcing_co2_const_co2,dt_unit="1d",                                     &
                       mul=forcing_options%conv_CO2_2_MassRatio, mis_dim=INPUT_VAR_DIM_MIS_SPREAD,     &
                       valid_time=INPUT_VAR_VALID_MIDPOINT,var_ref=CO2_var)
      CO2_unit = TRIM(InputAttrGet("units",var=CO2_var))
      IF (LEN_TRIM(CO2_unit) > 0) THEN
        status = GetStringIndex(TRIM(CO2_unit),"ppmv mol_per_mol kg_per_kg",ErrMsg="CO2 unit")
        SELECT CASE(status)
          CASE(1)
            forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg * 1.e-6_dp
          CASE(2)
            forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg 
          CASE(3)
            forcing_options%conv_CO2_2_MassRatio = 1.0_dp
        END SELECT
        IF ( status > 0 ) THEN ! If unit available and recognized: Override default unit conversion
          CALL InputOptNew(opt)
          opt%var_name = "CO2"
          opt%mul = forcing_options%conv_CO2_2_MassRatio
        ENDIF
      ENDIF
      forcing_input%CO2_concentr%frequency = GHG_SCENARIO_
    END IF                                                                                                                         
    IF (forcing_wind_frequ == "TIMESTEP") THEN
        ALLOCATE(forcingFields%wind_speed(domain%nland),STAT=status)
        IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%wind_speed() could not be allocated.")
    ENDIF
    ALLOCATE(forcingFields%wind_speed10(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%wind_speed10() could not be allocated.")

    ! --- initialization of diagnostic output streams ------------------------------------------------------------------------------

    CALL init_diagnostic_output(grid, domain, global_options%filetype)

    ! --- allocation of remaining forcing fields -----------------------------------------------------------------------------------

    ALLOCATE(forcingFields%air_pressure(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("init_forcing", "ERROR: Field forcingFields%air_pressure() could not be allocated.")

  END SUBROUTINE init_forcing

  ! --- read_forcing() -----------------------------------------------------------------------------
  !
  ! The code of this routine was formerly part of init_forcing. It is moved to this routine to work
  ! with multi-annual simulations and annual forcing files.
  !
  ! Routine to read forcing date from forcing data files.
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE read_forcing(grid, domain, forcingFields)

    USE mo_time_control,  ONLY: get_date_components, current_date, next_date, lstart, lresume
    USE mo_util_string, ONLY: tolower

    TYPE(grid_type),    INTENT(in)    :: grid
    TYPE(domain_type),  INTENT(in)    :: domain
    TYPE(forcing_type), INTENT(inout) :: forcingFields   !! This routine allocates memory for this structure

    INTEGER :: current_year, current_month, current_day, next_year, next_month, next_day
    INTEGER :: status
    REAL(dp) :: hlp_r

    ! We assume annual forcing files. Forcing data of next_date is read. We need to read from a new file, if
    ! next_date /= current_date, otherwise mothing needs to be done.
    CALL get_date_components (current_date, year=current_year, month=current_month, day=current_day)
    CALL get_date_components (next_date,  year=next_year, month=next_month, day=next_day)
    new_year_forcing =  (next_year /= current_year)
    new_month_forcing = (next_month /= current_month)
    new_day_forcing =   (next_day /= current_day)
    
    IF (.NOT. lstart .AND..NOT. lresume .AND..NOT. new_year_forcing) RETURN

    ! --- Temperature -----------------------------------------------------------------------------------------

    IF(forcing_temp_frequ /= "TIMESTEP") THEN
       CALL init_variable_access(domain, grid, "FORCING_TEMP", "tmin", forcing_temp_frequ, forcing_temp_file, &
            forcing_temp_const_tmin, forcing_input%air_temp_min)
       CALL init_variable_access(domain, grid, "FORCING_TEMP", "tmax", forcing_temp_frequ, forcing_temp_file, &
            forcing_temp_const_tmin, forcing_input%air_temp_max)
    ELSE
       CALL init_variable_access(domain, grid, "FORCING_TEMP", "air_temp", forcing_temp_frequ, forcing_temp_file, &
            forcing_temp_const_tmin, forcing_input%air_temp)
    END IF

    ! --- Precipitation ---------------------------------------------------------------------------------------

    CALL init_variable_access(domain, grid, "FORCING_PRECIP", "precip", forcing_precip_frequ, forcing_precip_file, &
         forcing_precip_const_precip, forcing_input%precipitation)
    IF(forcing_input%precipitation%frequency == MONTHLY_) &   !! Number of wet days only needed for monthly data
         CALL init_variable_access(domain, grid, "FORCING_PRECIP", "fwet", forcing_precip_frequ, forcing_precip_file, &
         HUGE(0.0_dp), forcing_input%frac_wet_days)

    ! eventually convert monthly precipitation data to kg/m^2/s:
    IF(forcing_options%precip_in_mm_per_day .AND. forcing_input%precipitation%frequency == MONTHLY_) &
         forcing_input%precipitation%data_monthly = forcing_input%precipitation%data_monthly * conv_prec_day

    ! --- Atmospheric humidity --------------------------------------------------------------------------------

    SELECT CASE(forcing_options%type_of_qair_forcing)
       CASE(NONE_)
          CALL message('read_forcing','Atmospheric humidity is diagnosed from Air temperature')
       CASE(RH_)
          ! need first call in allocate forcing_input%qair%data_daily, as this is needed to generate the diunal cycle of qair
          NULLIFY(forcing_input%qair%data_daily)
          ALLOCATE(forcing_input%qair%data_daily(domain%nland),STAT=status)
          IF(status /= 0) CALL finish('read_forcing', &
               'ERROR: For variable qair memory for daily data could not be allocated.')

          CALL init_variable_access(domain, grid, "FORCING_QAIR", "rel_humidity", forcing_qair_frequ, &
               forcing_qair_file, forcing_qair_const_rh, forcing_input%rel_humidity)

          SELECT CASE(forcing_input%rel_humidity%frequency)
             CASE(MONTHLY_)
                ! check data
                hlp_r = MAXVAL(forcing_input%rel_humidity%data_monthly(:,:))
                IF(hlp_r < 1._dp) THEN
                   WRITE (message_text,*) 'Maximum value of monthly relative humidity is ', hlp_r
                   CALL message('read_forcing', message_text)
                   CALL message('read_forcing','WARNING: rel. humidity data probably represent fractions instead of percentages !!')
                END IF
                !! Convert percentage data to fractional data:
                forcing_input%rel_humidity%data_monthly(:,:) = forcing_input%rel_humidity%data_monthly(:,:)/100._dp
             CASE(DAILY_)
                CALL message('read_forcing','Opening file for reading daily relative humidity')
             CASE(TIMESTEP_)
                CALL message('read_forcing','Opening file for reading timestep relative humidity')
             CASE(CONST_)
                hlp_r = MAXVAL(forcing_input%rel_humidity%data_daily(:))
                IF(hlp_r < 1._dp) THEN
                   WRITE (message_text,*) 'Constant value of relative humidity is ', hlp_r
                   CALL message('read_forcing', message_text)
                   CALL message('read_forcing',&
                        'WARNING: Constant value of relative humidity probably represents a fraction instead of a percentage!!')
                END IF
                !! Convert percentage data to fractional data:
                forcing_input%rel_humidity%data_daily(:) = forcing_input%rel_humidity%data_daily(:)/100._dp
             CASE default
                CALL finish('read_forcing','PROGRAMMING ERROR A54b1')
          END SELECT
       CASE(QAIR_)
          IF(forcing_input%qair%frequency == MONTHLY_) & 
               CALL finish('read_forcing','Qair forcing is currently not implemented on monthly forcing')
 
          CALL init_variable_access(domain, grid, "FORCING_QAIR", "qair", forcing_qair_frequ, &
               forcing_qair_file, forcing_qair_const_rh, forcing_input%qair)

       CASE default
          WRITE (message_text,*) 'PROGRAMMING ERROR: Value ', forcing_options%type_of_qair_forcing, ' should be impossible!!'
          CALL finish('read_forcing', message_text)
    END SELECT
    
    ! --- Shortwave radiation -----------------------------------------------------------------------------------

    SELECT CASE(forcing_options%type_of_shortwave_forcing)
       CASE(CLOUD_)
          CALL init_variable_access(domain, grid, "FORCING_SW_RADIAT", "cloud", forcing_sw_frequ, &
               forcing_sw_file, forcing_sw_const_cloud, forcing_input%cloud_cover)
          SELECT CASE(forcing_input%cloud_cover%frequency)
             CASE(MONTHLY_)
                ! check data
                hlp_r = MAXVAL(forcing_input%cloud_cover%data_monthly(:,:))
                IF(hlp_r < 1._dp) THEN
                   WRITE (message_text,*) 'Maximum value of monthly cloud cover is ', hlp_r
                   CALL message('read_forcing', message_text)
                   CALL message('read_forcing','WARNING: cloud cover data probably represent fractions instead of percentages !!')
                END IF
                !! Convert percentage data to fractional data:
                forcing_input%cloud_cover%data_monthly(:,:) = forcing_input%cloud_cover%data_monthly(:,:)/100._dp
             CASE(DAILY_)
                hlp_r = MAXVAL(forcing_input%cloud_cover%data_daily(:))
                IF(hlp_r < 1._dp) THEN
                   WRITE (message_text,*) 'Maximum value of daily cloud cover is ', hlp_r
                   CALL message('read_forcing', message_text)
                   CALL message('read_forcing','WARNING: cloud cover data probably represent fractionsinstead of percentages!!')
                END IF
                !! Convert percentage data to fractional data:
                forcing_input%cloud_cover%data_daily(:) = forcing_input%cloud_cover%data_daily(:)/100._dp
             CASE(CONST_)
                hlp_r = MAXVAL(forcing_input%cloud_cover%data_daily(:))
                IF(hlp_r < 1._dp) THEN
                   WRITE (message_text,*) 'Constant value of cloud cover is', hlp_r
                   CALL message('read_forcing', message_text)
                   CALL message('read_forcing',&
                        'WARNING: Constant value of cloud cover probably represents a fraction instead of a percentage!!')
                END IF
                !! Convert percentage data to fractional data:
                forcing_input%cloud_cover%data_daily(:) = forcing_input%cloud_cover%data_daily(:)/100._dp
             CASE(TIMESTEP_)
                CALL finish('read_forcing',"Cloud cover not implemented as time-step shortwave-forcing, force with radiation data")
             CASE default
                CALL finish('read_forcing',"PROGRAMMING ERROR A54b2")
          END SELECT
             
       CASE(MEAN_RAD_)
          ! read in shortwave data 
          CALL init_variable_access(domain, grid, "FORCING_SW_RADIAT", "shortwave", forcing_sw_frequ, &
               forcing_sw_file, forcing_sw_const_shortwave, forcing_input%shortwave)

          ! As a provisional solution we need here also a table of potential shortwave radiation 
          IF(forcing_options%type_of_shortwave_scheme /= SW_SCHEME_DIRECT_ ) THEN
                CALL init_variable_access(domain, grid, "FORCING_TABLE_SW_POT", "mpot", forcing_table_sw_pot_frequ, &
                   forcing_table_sw_pot_file, forcing_table_sw_pot_const, tables%shortwave_pot)
          ELSE
                ! not actually needed, but alloated here to avoid too many if statements in the main code
                ! evenually this tables options needs to be dropped anyhow
                CALL init_variable_access(domain, grid, "FORCING_TABLE_SW_POT", "mpot", forcing_table_sw_pot_frequ, &
                   forcing_table_sw_pot_file, forcing_table_sw_pot_const, tables%shortwave_pot,read_variable=.FALSE.)
          END IF

       CASE default
          WRITE (message_text,*) 'PROGRAMMING ERROR: Value ', forcing_options%type_of_shortwave_forcing, ' should be impossible!!'
          CALL finish('read_forcing', message_text)
    END SELECT
    ! read in longwave-radiation forcing options (type of radiation input)

    ! Allocate memory for longwave forcing

    SELECT CASE(forcing_options%type_of_longwave_forcing)
       CASE(CLOUD_) ! Use cloud data for longwave forcing
          ! Allocate memory for longwave forcing
          ALLOCATE(forcingFields%rad_lw_down(domain%nland),STAT=status)
          IF(status /= 0) CALL finish('read_forcing', "ERROR: Field forcingFields%rad_lw_down() could not be allocated.")

          !initialize cloud data
          IF(.NOT. forcing_input%cloud_cover%initialized) THEN! Initialize cloud data only if not already done for shortwave input 

             CALL init_variable_access(domain, grid, "FORCING_LW_RADIAT", "cloud", forcing_lw_frequ, &
                  forcing_lw_file, forcing_lw_const_cloud, forcing_input%cloud_cover)
             SELECT CASE(forcing_input%cloud_cover%frequency)
             CASE(MONTHLY_)
                ! check data
                hlp_r = MAXVAL(forcing_input%cloud_cover%data_monthly(:,:))
                IF(hlp_r < 1._dp) THEN
                   WRITE (message_text,*) 'Maximum value of monthly cloud cover is ', hlp_r
                   CALL message('read_forcing', message_text)
                   CALL message('read_forcing','WARNING: cloud cover data probably represent fractions instead of percentages !!')
                END IF
                !! Convert percentage data to fractional data:
                forcing_input%cloud_cover%data_monthly(:,:) = forcing_input%cloud_cover%data_monthly(:,:)/100._dp
             CASE(DAILY_)
                hlp_r = MAXVAL(forcing_input%cloud_cover%data_daily(:))
                IF(hlp_r < 1._dp) THEN
                   WRITE (message_text,*) 'Maximum value of daily cloud cover is ', hlp_r
                   CALL message('read_forcing', message_text)
                   CALL message('read_forcing','WARNING: cloud cover data probably represent fractions instead of percentages!!')
                END IF
                !! Convert percentage data to fractional data:
                forcing_input%cloud_cover%data_daily(:) = forcing_input%cloud_cover%data_daily(:)/100._dp
             CASE(CONST_)
                hlp_r = MAXVAL(forcing_input%cloud_cover%data_daily(:))
                IF(hlp_r < 1._dp) THEN
                   WRITE (message_text,*) 'Constant value of cloud cover is ', hlp_r
                   CALL message('read_forcing', message_text)
                   CALL message('read_forcing',&
                        'WARNING: Constant value of cloud cover probably represents a fraction instead of a percentage!!')
                END IF
                !! Convert percentage data to fractional data:
                forcing_input%cloud_cover%data_daily(:) = forcing_input%cloud_cover%data_daily(:)/100._dp
             CASE(TIMESTEP_)
                CALL finish('read_forcing',"Cloud cover not implemented as time-step longwave-forcing, force with radiation data")
             CASE default
                CALL finish('read_forcing',"PROGRAMMING ERROR F83-jj-23145F")
             END SELECT
          ELSE !! Cloudcover is used for both shortwave and longwave --> check frequency compatibility
             IF (p_parallel_io) THEN 
                IF (forcing_lw_frequ /= forcing_sw_frequ) THEN
                   CALL message('read_forcing', &
                        "Inconsistent frequencies for cloud data (keywords FORCING_LW_FREQU and FORCING_SW_FREQU in namelist)")
                END IF
             END IF
          END IF
       CASE(MEAN_RAD_) ! Use mean values for longwave forcing
          ! read in longwave data
          CALL init_variable_access(domain, grid, "FORCING_LW_RADIAT", "longwave", forcing_lw_frequ, &
               forcing_lw_file, forcing_lw_const_longwave, forcing_input%longwave)

          ALLOCATE(forcingFields%rad_lw_down(domain%nland),STAT=status)
          IF(status /= 0) CALL finish('read_forcing', "ERROR: Field forcingFields%rad_lw_down() could not be allocated.")
    END SELECT

    ! --- CO2 -------------------------------------------------------------------------------------------------

    IF (tolower(TRIM(forcing_co2_frequ)) /= "ghg_scenario") THEN
      CALL init_variable_access(domain, grid, "FORCING_CO2", "CO2", forcing_co2_frequ, forcing_co2_file, &
           forcing_co2_const_co2, forcing_input%CO2_concentr)
      IF (forcing_input%CO2_concentr%frequency /= TIMESTEP_) THEN ! TIMESTEP done in update_
          forcingFields%CO2_concentr => forcing_input%CO2_concentr%data_daily(:) !! For CO2-conc at all time steps the daily 
                                                                                 !! values are used
      END IF 
    ENDIF
                                                                                                                        
    ! --- Wind ----------------------------------------------------------------------------------------------------

    CALL init_variable_access(domain, grid, "FORCING_WIND", "wspeed", forcing_wind_frequ, forcing_wind_file, &
         forcing_wind_const_wspeed, forcing_input%wind_speed)
    IF (forcing_input%wind_speed%frequency /= TIMESTEP_) THEN
        forcingFields%wind_speed => forcing_input%wind_speed%data_daily(:) !! For wind speed at all time steps daily values are used
    ENDIF

  END SUBROUTINE read_forcing

  ! === update_forcing() ===========================================================================================================
  !
  ! Updates the central structure "forcing", i.e. it actualizes all forcing fields. This routine should be called every time step
  ! befor the land surface scheme is called.
  !
  ! More precisely
  !    (1) For monthly input data, each time, the time step falls in a new day, new daily mean values are computed,
  !        whereas for daily input data they are read in from file.
  !        ... at timestep forcing, at the beginning of each day, input files are read for the enitre day and stored in memory
  !    (2) At each time step the forcing fields are estimated from daily values.
  !
  SUBROUTINE update_forcing(grid, domain, evap_act2pot, forcingFields)

    USE mo_time_control,    ONLY: next_date, write_date, get_date_components, delta_time 
    USE mo_time_conversion, ONLY: time_days, day_in_year, IMerge_HMS2Sec
    USE mo_zenith,          ONLY: cos_zenith

    USE mo_atmosphere,   ONLY: sat_specific_humidity

    ! --- interface variables ---

    TYPE(domain_type),             INTENT(in)  :: domain        ! information on processor domain
    TYPE(grid_type),               INTENT(in)  :: grid          ! information on global grid
    REAL(dp),DIMENSION(1:domain%nland), INTENT(in) :: evap_act2pot ! ratio of actual to potential evapotranspiration from the
                                                                   ! previous day (needed for vapor pressure)
    TYPE(forcing_type), INTENT(inout)          :: forcingFields ! The forcing fields updated by this subroutine
    ! TYPE (time_days)   :: date                                ! the considered day in days-seconds format (actually, 
    !                                                            ! the seconds are not needed here)
    ! --- local variables ------------------------

    INTEGER          :: day, month, day_of_year          ! day of month, month, day of year at the actual time step

    REAL(dp),SAVE,ALLOCATABLE :: precip_generator_output(:,:) ! array to save the daily precipitation rate computed in the ...
                                                              ! ... precipitation_generator at the beginning of each month
    REAL(dp),SAVE,ALLOCATABLE :: relhum_generator_output(:,:) ! array to save the daily relative humidity computed in the ...
                                                              ! ... relhum_generator at the beginning of each month
    LOGICAL          :: preserve_monthly_precip = .TRUE. ! indicates, whether the daily precipitation rate is scaled to preserve ...
                                                         ! ... the monthly input values
    INTEGER          :: ith_timestep, t_temp    ! helpers for reading time step forcing
    INTEGER          :: hr, mn, se  ! integers indicating current time
    INTEGER          :: status

    REAL(dp),ALLOCATABLE,DIMENSION(:) :: hlp1,hlp2 ! helpers to swap between timestep and daily meteorological variables

    IF (debug) CALL message('update_forcing','Entering ...')
    IF (debug) CALL write_date(next_date,   'update_forcing: start - date: ')

    ! extract time information
    CALL get_date_components(next_date, MONTH=month, DAY=day, HOUR=hr, MINUTE=mn, SECOND=se)
    day_of_year=day_in_year(next_date)
    t_temp  = IMerge_HMS2Sec(hr, mn, se)
    ith_timestep = t_temp/int(delta_time)+1

    ! -- allocate helpers
    ALLOCATE(hlp1(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("update_forcing", "ERROR: Field hlp1() could not be allocated.")
    ALLOCATE(hlp2(domain%nland),STAT=status)
    IF(status /= 0) CALL finish("update_forcing", "ERROR: Field hlp2() could not be allocated.")

    ! --- temperature ------------------

    ! generate or read in daily data

    IF (new_day_forcing .OR. lstart .OR. lresume) THEN
       IF (forcing_input%air_temp_min%frequency == MONTHLY_) THEN !! monthly input data ==> linear interpolation to daily values
          CALL daily_from_monthly(forcing_input%air_temp_min%data_monthly,next_date,&
                                  -tmelt,0._dp,.TRUE.,.FALSE.,&
                                  forcing_input%air_temp_min%data_daily)
          CALL daily_from_monthly(forcing_input%air_temp_max%data_monthly,next_date,&
                                  -tmelt,0._dp,.TRUE.,.FALSE.,&
                                  forcing_input%air_temp_max%data_daily)

       ELSE IF( forcing_input%air_temp_min%frequency == DAILY_) THEN !! Daily input data ==> read in each new day
          CALL read_in_daily_data(grid,domain,forcing_input%air_temp_min)
          CALL read_in_daily_data(grid,domain,forcing_input%air_temp_max)

       ELSE IF( forcing_input%air_temp%frequency == TIMESTEP_) THEN
          CALL read_in_timestep_data(grid,domain,forcing_input%air_temp)
       END IF
    END IF

    ! estimate temperature at particular time step from daily temperature maximum and temperature minimum
    IF (forcing_input%air_temp%frequency /= TIMESTEP_) THEN
       CALL instantly_from_daily_temp_2(forcing_input%air_temp_min%data_daily(1:domain%nland),&
                                     forcing_input%air_temp_max%data_daily(1:domain%nland),&
                                     next_date, &
                                     forcingFields%air_temp(1:domain%nland))
    ELSE
       ! get the data for the ith time step from the the data_timestep field
       forcingFields%air_temp(1:domain%nland) = forcing_input%air_temp%data_timestep(1:domain%nland,ith_timestep)

    END IF

    ! --- compute bottom pressure --------------
    forcingFields%air_pressure(1:domain%nland)= bottom_pressure(domain%elev(1:domain%nland),forcingFields%air_temp(1:domain%nland))

    ! --- compute precipitation and snow ----------------

    ! For a whole day precipitation is kept constant (if not timestep). 
    ! Hence precipitation has to be updated only in case the day has changed from
    ! the previous to the current time step.
    ! in time step mode precipitation is changed every time

    IF (new_day_forcing .OR. lstart .OR. lresume) THEN ! each new day ...
       SELECT CASE(forcing_input%precipitation%frequency)
          CASE(MONTHLY_) ! Precipitation input data are monthly (and already converted to [kg/m^2/s])
             ! Precipitation is computed stochastically from the fraction of wet days and the average precipitation in the
             ! particular month (this is done in the subroutine precipitation_generator() at the first time step of each month)
             IF (new_month_forcing .OR. lstart .OR. lresume) THEN
                IF (.NOT. ALLOCATED(precip_generator_output)) ALLOCATE (precip_generator_output(domain%nland,31))
                CALL precipitation_generator(forcing_input%precipitation%data_monthly(:,month),&
                                             forcing_input%frac_wet_days%data_monthly(:,month),&
                                             next_date,preserve_monthly_precip,precip_generator_output(:,:))
             ENDIF
             forcing_input%precipitation%data_daily(:) = precip_generator_output(:,day)

          CASE(DAILY_) ! Precipitation input data are daily
             CALL read_in_daily_data(grid,domain,forcing_input%precipitation)
             !eventually convert precipitation data to kg/m^2/s:
             IF(forcing_options%precip_in_mm_per_day) &
                  forcing_input%precipitation%data_daily = forcing_input%precipitation%data_daily * conv_prec_day

          CASE(TIMESTEP_)
             CALL read_in_timestep_data(grid,domain,forcing_input%precipitation)

             ! In case precipitation is in mm/day, and convert it here to kg/m^2/s
             IF(forcing_options%precip_in_mm_per_day) &
                 forcing_input%precipitation%data_timestep = forcing_input%precipitation%data_timestep * conv_prec_day

          CASE(CONST_)
             ! Nothing to do

          CASE default
             CALL finish('update_forcing','PROGRAMMING ERROR P45Yb4-21')

       END SELECT

       ! the precipitation generator yields the TOTAL precipitation (i.e. rain+snow) --- so here, for the moment, 
       ! forcingFields%precip_rain is set to the total precipitation to save memory. 
       ! Divide total precipitation into rain and snow:
       IF (forcing_input%precipitation%frequency /= TIMESTEP_) THEN
          IF (forcing_input%air_temp%frequency == TIMESTEP_) THEN
             hlp1=MINVAL(forcing_input%air_temp%data_timestep(1:domain%nland,:), DIM=2)
             hlp2=MAXVAL(forcing_input%air_temp%data_timestep(1:domain%nland,:), DIM=2)
          ELSE
             hlp1=forcing_input%air_temp_min%data_daily(1:domain%nland)
             hlp2=forcing_input%air_temp_max%data_daily(1:domain%nland)
          ENDIF
          CALL snow_and_rain_from_precip(&
                  forcingFields%precip_rain(1:domain%nland), & ! output
                  forcingFields%precip_snow(1:domain%nland), & ! output
                  forcing_input%precipitation%data_daily(1:domain%nland), & ! input
                  0.5_dp*(hlp1(1:domain%nland)) + hlp2(1:domain%nland))  ! input 
       END IF
    END IF

    !! If precipitation inputs are time step we can change precip every time step...
    IF (forcing_input%precipitation%frequency == TIMESTEP_) THEN

       ! forcingFields%precip_rain is set to the total precipitation to save memory. 
       ! Divide total precipitation into rain and snow:
       IF (forcing_input%air_temp%frequency /= TIMESTEP_) THEN
          CALL snow_and_rain_from_precip(&
               forcingFields%precip_rain(1:domain%nland), & ! output
               forcingFields%precip_snow(1:domain%nland), & ! output
               forcing_input%precipitation%data_timestep(1:domain%nland,ith_timestep), & ! input
               0.5_dp*(forcing_input%air_temp_min%data_daily(1:domain%nland) + &
                       forcing_input%air_temp_max%data_daily(1:domain%nland))& 
               )
       ELSE
          CALL snow_and_rain_from_precip(&
               forcingFields%precip_rain(1:domain%nland), & ! output
               forcingFields%precip_snow(1:domain%nland), & ! output
               forcing_input%precipitation%data_timestep(1:domain%nland,ith_timestep), & ! input
               forcingFields%air_temp(1:domain%nland) ) ! input
       END IF

    END IF             

    ! --- computation of short wave radiation ---

    ! For a whole day cloud cover is kept constant.

    SELECT CASE(forcing_options%type_of_shortwave_forcing)
    CASE(CLOUD_)   ! short wave radiation is computed from cloudcover
          
       IF (new_day_forcing .OR. lstart .OR. lresume) THEN ! each new day
          SELECT CASE(forcing_input%cloud_cover%frequency)
             CASE(MONTHLY_)
                CALL daily_from_monthly(forcing_input%cloud_cover%data_monthly(:,:),next_date,&
                                        0._dp,100._dp,.TRUE.,.TRUE.,&
                                        forcing_input%cloud_cover%data_daily(:))! interpolation to a daily value  

             CASE(DAILY_)
                CALL read_in_daily_data(grid,domain,forcing_input%cloud_cover)

             CASE(TIMESTEP_)
                CALL finish('update_forcing','percentage of cloud cover at timestep subdaily not supported')

             CASE(CONST_)
             ! Nothing to do

             CASE default
                CALL finish('update_forcing','PROGRAMMING ERROR P24ghy/5')

          END SELECT
       END IF

       CALL shortwave_from_cloudcover(day_of_year, &                                        ! input
                                      cos_zenith                          (1:domain%nland),&  ! input
                                      forcingFields%air_pressure          (1:domain%nland),&  ! input
                                      forcing_input%cloud_cover%data_daily(1:domain%nland),&  ! input
                                      forcingFields%rad_UV_down           (1:domain%nland),&  ! output
                                      forcingFields%rad_PAR_down          (1:domain%nland),&  ! output
                                      forcingFields%rad_NIR_down          (1:domain%nland),&  ! output
                                      forcingFields%frac_PAR_diffuse      (1:domain%nland),&  ! output
                                      forcingFields%rad_sw_down           (1:domain%nland),&  ! output
                                      forcingFields%rad_sw_down_pot       (1:domain%nland) &  ! output
                                      )
    CASE(MEAN_RAD_) ! short wave radiation is computed from mean daily or monthly shortwave radiation incident to earth's surface

       IF (new_day_forcing .OR. lstart .OR. lresume) THEN
          ! update shortwave data
          SELECT CASE(forcing_input%shortwave%frequency)
             CASE(MONTHLY_) ! In case of monthly input data generate daily data
                CALL daily_from_monthly(forcing_input%shortwave%data_monthly(:,:),next_date, &
                                        0._dp,0._dp,.TRUE.,.FALSE.,&
                                        forcing_input%shortwave%data_daily(:))! interpolation to a daily value 

                IF(forcing_options%type_of_shortwave_scheme /= SW_SCHEME_DIRECT_) &
                   CALL daily_from_monthly(tables%shortwave_pot%data_monthly(:,:),next_date, &
                                        0._dp,0._dp,.TRUE.,.FALSE.,&
                                        tables%shortwave_pot%data_daily(:))! interpolation to a daily value 

             CASE(DAILY_) !  In case of daily input data read in the data for a new day
                CALL read_in_daily_data(grid,domain,forcing_input%shortwave)
                ! read potential radiation as well
                IF(forcing_options%type_of_shortwave_scheme /= SW_SCHEME_DIRECT_) &
                   CALL read_in_daily_data(grid,domain,tables%shortwave_pot)

             CASE(TIMESTEP_)
                CALL read_in_timestep_data(grid,domain,forcing_input%shortwave)
                ! read potential radiation as well
                IF(forcing_options%type_of_shortwave_scheme /= SW_SCHEME_DIRECT_) &
                   CALL read_in_timestep_data(grid,domain,tables%shortwave_pot)

             CASE(CONST_)
                ! Nothing to do

             CASE default
                CALL finish('update_forcing','PROGRAMMING ERROR Q4187bX-87')
          END SELECT

       END IF ! Preparations at a new day finished

       IF (forcing_input%shortwave%frequency /= TIMESTEP_) THEN
          hlp1(1:domain%nland) = forcing_input%shortwave%data_daily(1:domain%nland)
          hlp2(1:domain%nland) = tables%shortwave_pot%data_daily(1:domain%nland)
       ELSE
          hlp1(1:domain%nland) = forcing_input%shortwave%data_timestep(1:domain%nland,ith_timestep)
          hlp2(1:domain%nland) = tables%shortwave_pot%data_timestep(1:domain%nland,ith_timestep)
       ENDIF

          CALL shortwave_from_shortwaveMean(day_of_year, &                                       ! input
                                         cos_zenith                          (1:domain%nland),&  ! input
                                         domain%coslat                       (1:domain%nland),&  ! input
                                         domain%sinlat                       (1:domain%nland),&  ! input
                                         forcingFields%air_pressure          (1:domain%nland),&  ! input
                                         hlp1                                (1:domain%nland),&  ! input
                                         hlp2                                (1:domain%nland),&  ! input
                                         forcingFields%rad_UV_down           (1:domain%nland),&  ! output
                                         forcingFields%rad_PAR_down          (1:domain%nland),&  ! output
                                         forcingFields%rad_NIR_down          (1:domain%nland),&  ! output
                                         forcingFields%frac_PAR_diffuse      (1:domain%nland),&  ! output
                                         forcingFields%rad_sw_down           (1:domain%nland),&  ! output
                                         forcingFields%rad_sw_down_pot       (1:domain%nland) &  ! output
                                         )

    CASE default
       CALL finish('update_forcing',"PROGRAMMING ERROR Y4K7a9")
    END SELECT

    ! --- compute specific humidity ---
    SELECT CASE(forcing_options%type_of_qair_forcing)
    CASE(NONE_)   ! relative humidity is calculated from air temperature

       IF (forcing_input%air_temp%frequency /= TIMESTEP_) THEN
         hlp1=forcing_input%air_temp_min%data_daily(1:domain%nland)
       ELSE
         hlp1=MINVAL(forcing_input%air_temp%data_timestep(1:domain%nland,:), DIM=2)
       ENDIF

       ! --- compute vapour pressure --- 
       ! needed below to compute the specific humidity and (if not from input data) to compute the downward longwave radiation
       forcingFields%vapor_pressure(1:domain%nland) = &
               vapour_pressure_from_evapor(domain, &
                                     evap_act2pot(1:domain%nland), &
                                     hlp1(1:domain%nland), & ! input daily min temperature
                                     forcingFields%air_temp(1:domain%nland))

       ! --- compute specific humidity --- 
       forcingFields%spec_humidity(1:domain%nland) = &
            specific_humidity(forcingFields%vapor_pressure(1:domain%nland),& ! input
                              forcingFields%air_pressure  (1:domain%nland) & ! input
                                       )
        
    CASE(RH_)
       IF ( new_day_forcing .OR. lstart .OR. lresume ) THEN
          SELECT CASE(forcing_input%rel_humidity%frequency)
             CASE(MONTHLY_)
                IF (new_month_forcing .OR. lstart .OR. lresume) THEN
                IF (.NOT. ALLOCATED(relhum_generator_output)) & 
                     ALLOCATE (relhum_generator_output(domain%nland,31))
                CALL relhum_generator(forcing_input%rel_humidity%data_monthly(:,month),&
                                             forcing_input%frac_wet_days%data_monthly(:,month), &
                                             precip_generator_output(:,:),&
                                             next_date,relhum_generator_output(:,:))
                ENDIF
                forcing_input%rel_humidity%data_daily(:) = relhum_generator_output(:,day) 
             CASE(DAILY_)
                CALL read_in_daily_data(grid,domain,forcing_input%rel_humidity)
                forcing_input%rel_humidity%data_daily(1:domain%nland) = & 
                     forcing_input%rel_humidity%data_daily(1:domain%nland) / 100.
             CASE(TIMESTEP_)
                CALL read_in_timestep_data(grid,domain,forcing_input%rel_humidity)
                forcing_input%rel_humidity%data_timestep(1:domain%nland,:) = & 
                     forcing_input%rel_humidity%data_timestep(1:domain%nland,:) / 100.
             CASE(CONST_)
             ! Nothing to do
             CASE default
                CALL finish('update_forcing','PROGRAMMING ERROR XD4bHY-3')
          END SELECT

          IF (forcing_input%rel_humidity%frequency /= TIMESTEP_) THEN
             ! --- compute specific humidity --- 
             CALL spec_humidity_from_rel_humidity(forcing_input%rel_humidity%data_daily(1:domain%nland), &
                 forcing_input%air_temp_min%data_daily(1:domain%nland),&
                 forcing_input%air_temp_max%data_daily(1:domain%nland),&
                 domain%elev(1:domain%nland), &
                 delta_time, &
                 forcing_input%qair%data_daily(1:domain%nland))
          END IF
       END IF


       ! --- constrain daily specific humidity to be lower than or equal to the saturation spec. humidity
       hlp1=sat_specific_humidity(forcingFields%air_temp(1:domain%nland)+Tmelt, &
          forcingFields%air_pressure)

       IF (forcing_input%rel_humidity%frequency == TIMESTEP_) THEN
           forcingFields%spec_humidity(1:domain%nland) = & 
               forcing_input%rel_humidity%data_timestep(1:domain%nland,ith_timestep)*hlp1      
       ELSE
           forcingFields%spec_humidity(1:domain%nland) = & 
               MIN(forcing_input%qair%data_daily(1:domain%nland),hlp1)      
       ENDIF

       ! --- compute current vapour pressure --- 
       forcingFields%vapor_pressure(1:domain%nland) = &
           vapour_pressure(forcingFields%spec_humidity(1:domain%nland), & ! input
                           forcingFields%air_pressure  (1:domain%nland) & ! input
                                       )

        
    CASE(QAIR_)
       IF ( new_day_forcing .OR. lstart .OR. lresume ) THEN
          SELECT CASE(forcing_input%qair%frequency)
             CASE(MONTHLY_)
                CALL finish('update_forcing','QAIR at monthly time step not yet implemented')
             CASE(DAILY_)
                CALL read_in_daily_data(grid,domain,forcing_input%qair)
             CASE(TIMESTEP_)
                CALL read_in_timestep_data(grid,domain,forcing_input%qair)
             CASE(CONST_)
             ! Nothing to do
             CASE default
                CALL finish('update_forcing','PROGRAMMING ERROR XD4bHY-4')
          END SELECT
       END IF

       hlp1=sat_specific_humidity(forcingFields%air_temp(1:domain%nland)+Tmelt, &
          forcingFields%air_pressure)

       IF (forcing_input%shortwave%frequency /= TIMESTEP_) THEN
          forcingFields%spec_humidity(1:domain%nland) = & 
               MIN(forcing_input%qair%data_daily(1:domain%nland),hlp1)
       ELSE
          forcingFields%spec_humidity(1:domain%nland) = & 
               MIN(forcing_input%qair%data_timestep(1:domain%nland,ith_timestep),hlp1)
       ENDIF

       ! --- compute current vapour pressure --- 
       forcingFields%vapor_pressure(1:domain%nland) = &
           vapour_pressure(forcingFields%spec_humidity(1:domain%nland),& ! input
                           forcingFields%air_pressure  (1:domain%nland) & ! input
                                       )
 
    END SELECT


    ! --- compute downward longwave radiation ---

    SELECT CASE(forcing_options%type_of_longwave_forcing)
    CASE(CLOUD_)   ! long wave radiation is computed from cloudcover and vapour pressure
       IF ((new_day_forcing .OR. lstart .OR. lresume ).AND. (forcing_options%type_of_shortwave_forcing /= CLOUD_)) THEN
          SELECT CASE(forcing_input%cloud_cover%frequency)
            CASE(MONTHLY_)
                CALL daily_from_monthly(forcing_input%cloud_cover%data_monthly(:,:),next_date, &
                                        0._dp,100._dp,.TRUE.,.TRUE.,&
                                        forcing_input%cloud_cover%data_daily(:))! interpolation to a daily value  
             CASE(DAILY_)
                CALL read_in_daily_data(grid,domain,forcing_input%cloud_cover)
             CASE(CONST_)
             ! Nothing to do
             CASE(TIMESTEP_)
                CALL finish("update_forcing","CLOUD and TIMESTEP input not yet implemented for longwave rad")
             CASE default
                CALL finish('update_forcing','PROGRAMMING ERROR XD4bHY-3')
          END SELECT
       END IF

       forcingFields%rad_lw_down(1:domain%nland) = &
            longwave_rad_down(forcing_input%cloud_cover%data_daily(1:domain%nland), & ! input
                              forcingFields%vapor_pressure        (1:domain%nland), & ! input
                              forcingFields%air_temp              (1:domain%nland)  & ! input
                              )
      
    CASE(MEAN_RAD_) ! long wave radiation is taken from input data
       IF (new_day_forcing .OR. lstart .OR. lresume) THEN
          SELECT CASE(forcing_input%longwave%frequency)
          CASE(MONTHLY_) ! In case of monthly input data generate daily data
             CALL daily_from_monthly(forcing_input%longwave%data_monthly(:,:),next_date,&
                                     0._dp,0._dp,.TRUE.,.FALSE.,&
                                     forcing_input%longwave%data_daily(:))! interpolation to a daily value 
          CASE(DAILY_) !  In case of daily input data read in the data for a new day
             CALL read_in_daily_data(grid,domain,forcing_input%longwave)
          CASE(TIMESTEP_) 
             CALL read_in_timestep_data(grid,domain,forcing_input%longwave)
          CASE(CONST_)
             ! Nothing to do
          CASE default
             CALL finish('update_forcing','PROGRAMMING ERROR YPQwZ76b-89-4')
          END SELECT

          ! It is assumed that for the whole day the longwave radiation is constant, therefore in initialization ..
          ! .. forcingFields%rad_lw_down(:) is made pointing to forcing_input%longwave%data_daily(:) and no further ..
          ! .. computations are necessary
          ! .. this is differnt for TIMESTEP: longwave changes in each time step...

       END IF ! Preparations at a new day finished

       IF (forcing_input%longwave%frequency /= TIMESTEP_) THEN
           IF (forcing_input%shortwave%frequency /= TIMESTEP_ ) THEN
               hlp1 = forcing_input%air_temp_min%data_daily(1:domain%nland)
               hlp2 = forcing_input%air_temp_max%data_daily(1:domain%nland)
           ELSE 
               hlp1=MINVAL(forcing_input%air_temp%data_timestep(1:domain%nland,:), DIM=2)
               hlp2=MAXVAL(forcing_input%air_temp%data_timestep(1:domain%nland,:), DIM=2)
           ENDIF

           ! diurnal longwave downward radiation is generated such that the average radiation matches the 
           ! observations but the diurnal cycle follows temperature to help calculation of a correct
           ! net radiation flux. Note that reduces biases resulting from assuming a constant flux, but 
           ! this introduces uncertainty due to biases in the assumed dial course of temperature 
           forcingFields%rad_lw_down(1:domain%nland) = &
                longwave_from_daily_longwave(forcing_input%longwave%data_daily(1:domain%nland), & ! input
                                             hlp1(1:domain%nland),                              & ! input 
                                             hlp2(1:domain%nland),                              & ! input
                                             forcingFields%air_temp(1:domain%nland))              ! input               
       ELSE
          forcingFields%rad_lw_down(1:domain%nland) = forcing_input%longwave%data_timestep(1:domain%nland,ith_timestep)
       END IF

    CASE default
       CALL finish('update_forcing','PROGRAMMING ERROR Yppp4-OO-p3')
    END SELECT

    ! --- CO2-concentration -----------
    ! For CO2-concentration the daily values are used at each time step, i.e. the forcing field is declared as a
    ! pointer to the daily field: forcingFields%CO2_concentr => forcing_input%CO2_concentr%data_daily (see above)
    IF (new_day_forcing .OR. lstart .OR. lresume) THEN ! each new day ...
       SELECT CASE(forcing_input%CO2_concentr%frequency)
       CASE(MONTHLY_)
          CALL daily_from_monthly(forcing_input%CO2_concentr%data_monthly,next_date,&
                                  0._dp,0._dp,.TRUE.,.FALSE.,&
                                  forcing_input%CO2_concentr%data_daily)
          ! Convert CO2 to kg(CO2)/kg(dry air)
          forcing_input%CO2_concentr%data_daily = forcing_input%CO2_concentr%data_daily * & 
                                                  forcing_options%conv_CO2_2_MassRatio 
       CASE(DAILY_)
          CALL read_in_daily_data(grid,domain,forcing_input%CO2_concentr)
          ! Convert CO2 to kg(CO2)/kg(dry air)
          forcing_input%CO2_concentr%data_daily = forcing_input%CO2_concentr%data_daily * & 
                                                  forcing_options%conv_CO2_2_MassRatio 
       CASE(TIMESTEP_)
          CALL read_in_timestep_data(grid,domain,forcing_input%CO2_concentr)
          ! Convert CO2 to kg(CO2)/kg(dry air)
          forcing_input%CO2_concentr%data_timestep = forcing_input%CO2_concentr%data_timestep * & 
                                                     forcing_options%conv_CO2_2_MassRatio 
       CASE(CONST_)
          ! Nothing to do
       CASE(GHG_SCENARIO_)
          ! mo_input takes care
       CASE default
          CALL finish('update_forcing','PROGRAMMING ERROR: Program should never reach this point!!')       
       END SELECT
    END IF

    IF (forcing_input%CO2_concentr%frequency == TIMESTEP_) THEN
       forcingFields%CO2_concentr(1:domain%nland) = forcing_input%CO2_concentr%data_timestep(1:domain%nland,ith_timestep)
    END IF

    ! --- wind speed -----------
    ! For wind speed there is currently no interpolation from daily or monthly values to time step values, i.e. the forcing field 
    ! is declared as a pointer to the daily field: forcingFields%wind_speed => forcing_input%wind_speed%data_daily
    ! (see subroutine init_forcing())

    IF (new_day_forcing .OR. lstart .OR. lresume) THEN ! each new day
       select case(forcing_input%wind_speed%frequency)
          case(MONTHLY_)
             call daily_from_monthly(forcing_input%wind_speed%data_monthly,next_date,&
                                     0._dp,0._dp,.true.,.false.,&
                                     forcing_input%wind_speed%data_daily)
          case(DAILY_)
             call read_in_daily_data(grid,domain,forcing_input%wind_speed)
          case(CONST_)
             ! Nothing to do
          case(TIMESTEP_)
             CALL read_in_timestep_data(grid,domain,forcing_input%wind_speed)
          case default
             call finish('update_forcing','PROGRAMMING ERROR (wind_speed): Program should never reach this point!!')       
        end select
     end if

    IF (forcing_input%wind_speed%frequency == TIMESTEP_) THEN
       forcingFields%wind_speed(1:domain%nland) = forcing_input%wind_speed%data_timestep(1:domain%nland,ith_timestep)
    END IF

    ! --- update diagnostic output ---------
    CALL update_diagnostic_output(forcingFields)

    IF (debug) CALL message('update_forcing','... finished.')

    RETURN

    CONTAINS

    ! --- daily_from_monthly ----------------------------------------------------------------------------------------------
    ! 
    ! Computes for a given day of a year at each grid point an interpolated value.
    ! The interpolation is linear and preserving the monthly mean values.
    ! Thomas Raddatz 01.2004
    ! 
    ! ATTENTION: Interpolated values are not strictly in the range of the input values. So precipitation, wind speed might be
    !        negative, cloud cover larger than 1 or negative, etc. These out of range values should be excluded
    !        by setting lowerbound and/or upperbound, but then monthly mean values are not strictly preserved in every case.
    !        It can not be excluded, that tmin may be larger than tmax at some grid points. This should be checked after the
    !        interpolation (in the calling sequence). 
    !
    ! ATTENTION: This routine is faulty:
    !        The beginning of a year is computed here from the Julian/Gregorian calender, which is wrong by about 3 days in 10000 
    !        years. This may give wrong results for palaeo-computations. It would be better to detect the beginning of a year from
    !        the vernal equinox. To this end the EACHAM-Routine get_doy() from  the module mo_time_control could be used. --- 
    !        Particularly wrong is the year 1582, where the Julian switches to the Gregorian calender. Since in this year 
    !        October, 4, is immediately followed by October, 15, this routine gives discontinuous results for these days. 
    SUBROUTINE daily_from_monthly(monthly_values,date,lowerbound,upperbound,lowerbound_apply,upperbound_apply,daily_values)
      USE mo_time_conversion, ONLY: time_days
      USE mo_time_control, ONLY: get_date_components,get_month_len 
      IMPLICIT NONE

      REAL(dp), INTENT (IN)  :: monthly_values(domain%nland,12) ! one value for each month at each grid point
      TYPE (time_days)   :: date                                ! the considered day in days-seconds format (actually, 
                                                                ! the seconds are not needed here)
      REAL(dp), INTENT (IN)  :: lowerbound                      ! lower bound of the range of values
      REAL(dp), INTENT (IN)  :: upperbound                      ! upper bound of the range of values
      LOGICAL, INTENT (IN)   :: lowerbound_apply                ! determines whether the lower bound is regarded or not 
      LOGICAL, INTENT (IN)   :: upperbound_apply                ! determines whether the upper bound is regarded or not 
      REAL(dp), INTENT (OUT) :: daily_values(domain%nland)      ! the interpolated values for the considered day

      INTEGER  :: year,month,day,month_length,month_last_last,month_last,month_next,month_next_next
      INTEGER  :: days_first_half_month,days_second_half_month
      REAL(dp) :: diff_monthly(domain%nland)   ! difference in the value between two adjacent months
      REAL(dp) :: correct_fact(domain%nland)   ! factor to correct the start or end of the interpolation at the beginning or end of
                                               ! the current month
      REAL(dp) :: x1(domain%nland)             ! start of the interpolation at the beginning of the current month 
      REAL(dp) :: x2(domain%nland)             ! end of the interpolation at the end of the current month
      REAL(dp) :: xn(domain%nland)             ! value at the middle of the current month used for the interpolation

      CALL get_date_components(date,YEAR=year,MONTH=month,DAY=day)   ! extract month and day --> WRONG FOR PALAEO COMPUTATIONS!!!

      month_length = get_month_len (year, month)                     ! the length of the current month
      days_first_half_month = NINT(0.5_dp*REAL(month_length))        ! number of days from the beginning to the middle of the month
      days_second_half_month = month_length - days_first_half_month  ! number of days from the middle to the end of the month
      month_last = MOD(month+10,12) + 1                              ! ... determine the number of last month ...
      month_last_last = MOD(month_last+10,12) + 1                    ! ... determine the number of the last but one month ...
      month_next = MOD(month,12) + 1                                 ! ... determine the number of next month ...
      month_next_next = MOD(month_next,12) + 1                       ! ... determine the number of the next but one month ...
      diff_monthly(1:domain%nland) = monthly_values(1:domain%nland,month) - monthly_values(1:domain%nland,month_last)
      correct_fact(1:domain%nland) = (monthly_values(1:domain%nland,month_last_last) - &
                                      monthly_values(1:domain%nland,month_last) - &
                                      monthly_values(1:domain%nland,month) + &
                                      monthly_values(1:domain%nland,month_next)) / &
                                     (2._dp * MAX(ABS(diff_monthly(1:domain%nland)),negligible))
      x1(1:domain%nland) = ((monthly_values(1:domain%nland,month_last) + &
                           monthly_values(1:domain%nland,month)) / 2._dp) - &
                           ATAN(correct_fact(1:domain%nland)) * ABS(diff_monthly(1:domain%nland)) / pi
      diff_monthly(1:domain%nland) = monthly_values(1:domain%nland,month_next) - monthly_values(1:domain%nland,month)
      correct_fact(1:domain%nland) = (monthly_values(1:domain%nland,month_last) - &
                                      monthly_values(1:domain%nland,month) - &
                                      monthly_values(1:domain%nland,month_next) + &
                                      monthly_values(1:domain%nland,month_next_next)) / &
                                     (2._dp * MAX(ABS(diff_monthly(1:domain%nland)),negligible))
      x2(1:domain%nland) = ((monthly_values(1:domain%nland,month_next) + &
                           monthly_values(1:domain%nland,month)) / 2._dp) - &
                           ATAN(correct_fact(1:domain%nland)) * ABS(diff_monthly(1:domain%nland)) / pi
      xn(1:domain%nland) = (2._dp * REAL(month_length) * monthly_values(1:domain%nland,month) - &
                           x1(1:domain%nland) * REAL(days_first_half_month) - &
                           x2(1:domain%nland) * REAL(days_second_half_month)) / REAL(month_length) 
      IF (day <= days_first_half_month) THEN
         daily_values(1:domain%nland) = (REAL(day) - 0.5_dp) * &                               ! linear interpolation for the days .
                                        (xn(1:domain%nland) - x1(1:domain%nland)) / &          ! ... before the middle of the month
                                        REAL(days_first_half_month,dp) + x1(1:domain%nland)
      ELSE
         daily_values(1:domain%nland) = (REAL(day,dp) - 0.5_dp - REAL(days_first_half_month,dp)) * &
                                        (x2(1:domain%nland) - xn(1:domain%nland)) / &          ! ... after the middle of the month
                                        REAL(days_second_half_month,dp) + xn(1:domain%nland)
      END IF
      IF (lowerbound_apply) THEN
         daily_values(1:domain%nland) = MAX(daily_values(1:domain%nland),lowerbound)
      END IF
      IF (upperbound_apply) THEN
         daily_values(1:domain%nland) = MIN(daily_values(1:domain%nland),upperbound)
      END IF
      RETURN
    END SUBROUTINE DAILY_FROM_MONTHLY

    
    ! --- instantly_from_daily_temp_1() --------------------------------------------------------------------------------------
    !
    ! Given temperature mean and range for a particular day at every grid point,this routine 
    ! estimates the temperature for a particular timestep at that day for all grid points.
    ! In contrast to Marko's subroutine daytemp() (from which this code is derived) the routine
    ! correctly accounts for local time, i.e. temperatures are correct at all longitudes!
    !
    ! More precisely the temperature development is modelled as follows: 
    ! Let: T_min = T_mean - 0.5*T_range,  and: T_max = T_mean + 0.5*T_range. Then for a
    !     (i) normal day (more than 4, but less than 20 hours of daylight):
    !         --> during daylight: temperature behaves like a positive sine-wave, such that T_max is reached at 14:00 and 
    !             T_min at sunrise.
    !         --> during night: starting at sunset the temperature drops linearly until T_min at sunrise.
    !    (ii) short day (less than 4 hours of daylight):
    !         temperature is set to T_mean at all time steps, i.e. temperature is assumed to be constant.
    !   (iii) long day (more than 20 hours of daylight): 
    !         temperature is modelled as in the daylight-case of (i), but for the whole day, even during night hours.
    !
    ! REMARKS: 
    !      -- The time computations are based on the Julian/Gregorian calender. Hence there is a slight mismatch between the 
    !         correct begin of a year, as defined by the orbit around the sun, and the begin of the year according to the calender.
    !         The mismatch is about 3 days in 10000 years. This affects the computation of daylength in the current routine, and 
    !         thus the "phase" of the temperature curve, i.e. that e.g. the time when the temperature starts rising in the morning,
    !         may be somewhat shifted.
   
    SUBROUTINE instantly_from_daily_temp_1(dtmp,dtran,date,tmp)

      USE mo_kind,            ONLY: dp
      USE mo_time_conversion, ONLY: time_days
      USE mo_time_control,    ONLY: get_year_day
      USE mo_orbit,           ONLY: inquire_declination

      IMPLICIT NONE

      REAL(dp), PARAMETER :: pio180 = pi/180._dp   ! conversion factor from degrees to radians ("pi over 180")
      REAL(dp), PARAMETER :: rat1=14._dp/24._dp        ! day fraction at which maximum temperature shall occur
      REAL(dp), PARAMETER :: rat2=4._dp/24._dp         ! minimum length of normal day (as fraction of day)
      REAL(dp), PARAMETER :: rat3=20._dp/24._dp        ! maximum length of normal day (as fraction of day)
      REAL(dp), PARAMETER :: rat4=2._dp/24._dp         ! day fraction of dawn times (before sunrise + after sunset), used to set ...
                                                 ! ... temperature of a normal day to T_mean at sunrise-rat4/2 and sunset+rat4/2

      REAL(dp), INTENT (in) :: dtmp(1:domain%nland)   ! average temperature at the considered day for each grid point 
      REAL(dp), INTENT (in) :: dtran(1:domain%nland)  ! temperature range at that day for each grid point 
      TYPE (time_days), INTENT (in)  :: date          ! the considered simulation instant at that day in time-seconds format
      REAL(dp), INTENT (OUT) :: tmp(1:domain%nland)   ! estimated temperature at all grid points


      INTEGER    :: yearday           ! the number of the day of the current timestep in the current year
      REAL(dp)   :: year_instant      ! the time of the current time step in yearday-fraction-of-day format
      REAL(dp)   :: frac_timestep     ! the current time step expressed as the fraction of the current day
      REAL(dp)   :: frac_localtime    ! local time at a particular grid point expressed as the fraction of the current ...
                                      ! ... day since midnight (mod 1
      REAL(dp)   :: frac_UTC_offset   ! offset between UTC and local time (depends on longitude)
      REAL(dp)   :: declination       ! solar declination
      REAL(dp)   :: cpds, spds        ! cosine and sine contributions to zenith angle
      REAL(dp)   :: frac_daylight     ! fraction of day between sunrise and sunset
      REAL(dp)   :: frac_since_sunset ! time between current time step and sunset expressed as fraction of day
      REAL(dp)   :: frac_sunrise      ! time of sunrise expressed as fraction of day since midnight
      REAL(dp)   :: frac_sunset       ! time of sunset expressed as fraction of day since midnight
      REAL(dp)   :: temp_sunset       ! temperature at sunset
      REAL(dp)   :: temp_min          ! lowest day temperature
      REAL(dp)   :: sd                ! modulation factor for temperature
      REAL(dp)   :: arg
      INTEGER    :: i

      CALL inquire_declination(declination)

      year_instant=get_year_day(date)            ! conversion of current time step from day-seconds to yearday-dayfraction format
      yearday = INT(year_instant)                ! number of day of the current time step since the beginning of the current year
      frac_timestep=year_instant-REAL(yearday,dp)   ! fractional part of the day of the current time step since midnight in UTC
      DO i = 1,domain%nland
         frac_UTC_offset = domain%lon(i)/360._dp      ! difference between local time and UTC as fraction of a day 
         frac_localtime = MOD(frac_timestep+frac_UTC_offset+1,1._dp) ! local time as fraction of the local day (modulo 1)
         spds = domain%sinlat(i)*SIN(declination)  ! sine contribution to solar zenith angle
         cpds = domain%coslat(i)*COS(declination)  ! cosine contribution to solar zenith angle
         arg = -spds/cpds                        ! faulty?: for cpds=0 arg gets infinite
         IF (arg > 1._dp) THEN               ! polar night:
            frac_daylight = 0._dp            ! ... day length is zero
         ELSE IF (arg < -1._dp) THEN         ! polar day:
            frac_daylight = 1._dp            ! ... day length is the whole day
         ELSE                                ! normal day / night:
            frac_daylight = ACOS(arg)/pi    ! ... day length expressed as fraction of the whole day
         END IF         
         IF (frac_daylight>=rat2 .AND. frac_daylight<=rat3) THEN ! -------- normal day ---------------------------------------------
            frac_sunrise = 0.5_dp - frac_daylight/2._dp          ! time of sunrise expressed as fraction of day since midnight
            frac_sunset = frac_sunrise+frac_daylight             ! time of sunset expressed as fraction of day since midnight 
            IF (frac_localtime>frac_sunrise .AND. frac_localtime<frac_sunset) THEN  ! For time steps at daylight ...
               sd = COS(pi*(frac_localtime-rat1)/(frac_daylight/2._dp+rat4)) ! ... modulation factor for day-temperature such that
                                                                 ! .. temperature is at maximum at rat1 (2 pm) and at average at ...
                                                                 ! .. sunrise-rat4/2 and sunset+rat4/2
               tmp(i) = dtmp(i) + dtran(i)/2._dp*sd              ! temperature modulated by temperature range
            ELSE                                                 ! For time steps at night ...
               sd = COS(pi*(frac_sunset-rat1)/(frac_daylight/2._dp+rat4)) ! modulation factor for temperature at sunset
               temp_sunset = dtmp(i) + dtran(i)/2._dp*sd            ! temperature at sunset
               frac_since_sunset = MOD(frac_localtime-frac_sunset+1.,1._dp) ! day fraction since sunset
               temp_min = dtmp(i) - dtran(i)/2._dp                  ! lowest temperature
               tmp(i) = temp_min + (temp_sunset-temp_min)*(1.-frac_since_sunset/(1._dp-frac_daylight)) ! linear interpolation of ...
                                                                 ! ... temperature such that temperature is equal to the ...
                                                                 ! ... temperature at sunset and drops until the minimum ...
                                                                 ! ... temperature is met at sunrise.
            END IF
         ELSEIF (frac_daylight>rat3) THEN ! --- long day (close to polar day) ------------------------------------------------------
            sd = COS(pi*(frac_localtime-rat1)/(frac_daylight/2._dp+rat4)) ! modulation factor for day-temperature such that ...
                                                                 ! ... temperature is at maximum at rat1 (2 pm) and at average ...
                                                                 ! ... at sunrise-rat4/2 and sunset+rat4/2
            tmp(i) = dtmp(i) + dtran(i)/2._dp*sd                    ! temperature modulated by temperature range
         ELSE ! ------------------------------- short day (close to polar night) ---------------------------------------------------
            tmp(i) = dtmp(i)                                     ! temperature is set to mean temperture at all time steps
         END IF
      ENDDO
      
      RETURN
    END SUBROUTINE INSTANTLY_FROM_DAILY_TEMP_1
    
    ! --- instantly_from_daily_temp_2() --------------------------------------------------------------------------------------
    !
    ! The same as subroutine instantly_from_daily_temp_1(), but with temperature minimum and maximum as input instead of
    ! mean temperature and temperature range.
    ! In contrast to Marko's subroutine daytemp() (from which this code is derived) the routine
    ! correctly accounts for local time, i.e. temperatures are correct at all longitudes!
    !
    ! More precisely the temperature development is modelled as follows: 
    ! Let: T_min = T_mean - 0.5*T_range,  and: T_max = T_mean + 0.5*T_range. Then for a
    !     (i) normal day (more than 4, but less than 20 hours of daylight):
    !         --> during daylight: temperature behaves like a positive sine-wave, such that T_max is reached at 14:00 and 
    !             T_min at sunrise.
    !         --> during night: starting at sunset the temperature drops linearly until T_min at sunrise.
    !    (ii) short day (less than 4 hours of daylight):
    !         temperature is set to T_mean at all time steps, i.e. temperature is assumed to be constant.
    !   (iii) long day (more than 20 hours of daylight): 
    !         temperature is modelled as in the daylight-case of (i), but for the whole day, even during night hours.
    !
    ! REMARKS: 
    !      -- The time computations are based on the Julian/Gregorian calender. Hence there is a slight mismatch between the 
    !         correct begin of a year, as defined by the orbit around the sun, and the begin of the year according to the calender.
    !         The mismatch is about 3 days in 10000 years. This affects the computation of daylength in the current routine, and 
    !         thus the "phase" of the temperature curve, i.e. that e.g. the time when the temperature starts rising in the morning,
    !         may be somewhat shifted.  
    SUBROUTINE instantly_from_daily_temp_2(dtmp_min,dtmp_max,date,tmp)

      USE mo_kind,            ONLY: dp
      USE mo_time_conversion, ONLY: time_days
      USE mo_time_control,    ONLY: get_year_day
      USE mo_orbit,           ONLY: inquire_declination

      IMPLICIT NONE

      REAL(dp), PARAMETER :: pio180 = pi/180._dp    ! conversion factor from degrees to radians ("pi over 180")
      REAL(dp), PARAMETER :: rat1=14._dp/24._dp      ! day fraction at which maximum temperature shall occur
      REAL(dp), PARAMETER :: rat2=4._dp/24._dp       ! minimum length of normal day (as fraction of day)
      REAL(dp), PARAMETER :: rat3=24._dp/24._dp      ! maximum length of normal day (as fraction of day)
      REAL(dp), PARAMETER :: rat4=2._dp/24._dp       ! day fraction of dawn times (before sunrise + after sunset), used to set ...
                                                 ! ... temperature of a normal day to T_mean at sunrise-rat4/2 and sunset+rat4/2

      REAL(dp), INTENT (in) :: dtmp_min(1:domain%nland)  ! minimal temperature at the considered day for each grid point 
      REAL(dp), INTENT (in) :: dtmp_max(1:domain%nland)  ! maximal range at that day for each grid point 
      TYPE (time_days), INTENT (in)  :: date             ! the considered simulation instant at that day in time-seconds format
      REAL(dp), INTENT (OUT) :: tmp(1:domain%nland)      ! estimated temperature at all grid points


      INTEGER    :: yearday           ! the number of the day of the current timestep in the current year
      REAL(dp)   :: year_instant      ! the time of the current time step in yearday-fraction-of-day format
      REAL(dp)   :: frac_timestep     ! the current time step expressed as the fraction of the current day
      REAL(dp)   :: frac_localtime    ! local time at a particular grid point expressed as the fraction of the current ...
                                      ! ... day since midnight (mod 1
      REAL(dp)   :: frac_UTC_offset   ! offset between UTC and local time (depends on longitude)
      REAL(dp)   :: declination       ! solar declination
      REAL(dp)   :: cpds, spds        ! cosine and sine contributions to zenith angle
      REAL(dp)   :: frac_daylight     ! fraction of day between sunrise and sunset
      REAL(dp)   :: frac_since_sunset ! time between current time step and sunset expressed as fraction of day
      REAL(dp)   :: frac_sunrise      ! time of sunrise expressed as fraction of day since midnight
      REAL(dp)   :: frac_sunset       ! time of sunset expressed as fraction of day since midnight
      REAL(dp)   :: temp_sunset       ! temperature at sunset
      REAL(dp)   :: dtmp              ! mean day temperature
      REAL(dp)   :: dtran             ! day temperature range: maximum - minimum
      
      REAL(dp)   :: sd                ! modulation factor for temperature
      REAL(dp)   :: arg
      INTEGER    :: i

      CALL inquire_declination(declination)

      year_instant=get_year_day(date)            ! conversion of current time step from day-seconds to yearday-dayfraction format
      yearday = INT(year_instant)                ! number of day of the current time step since the beginning of the current year
      frac_timestep=year_instant-REAL(yearday,dp)   ! fractional part of the day of the current time step since midnight in UTC

      DO i = 1,domain%nland
         frac_UTC_offset = domain%lon(i)/360._dp   ! difference between local time and UTC as fraction of a day 
         frac_localtime = MOD(frac_timestep+frac_UTC_offset+2._dp,1._dp) ! local time as fraction of the local day (modulo 1)
         spds = domain%sinlat(i)*SIN(declination)  ! sine contribution to solar zenith angle
         cpds = domain%coslat(i)*COS(declination)  ! cosine contribution to solar zenith angle
         arg = -spds/cpds                        ! faulty?: for cpds=0 arg gets infinite
         IF (arg > 1._dp) THEN               ! polar night:
            frac_daylight = 0._dp            ! ... day length is zero
         ELSE IF (arg < -1._dp) THEN         ! polar day:
            frac_daylight = 1._dp            ! ... day length is the whole day
         ELSE                                ! normal day / night:
            frac_daylight = ACOS(arg)/pi    ! ... day length expressed as fraction of the whole day
         END IF

!! BLARPP empirical correction to account for the fact that the scheme does actually not 
!! consered dtmp=0.5*(max+min). The night-day turn temperature is adjusted such that the mean
!! temperature is within 0.05 K of the 0.5*(tmax+tmin) temperature
         IF (frac_daylight>=rat2 ) THEN ! -------- normal day ---------------------------------------------
           IF(frac_daylight < 0.36_dp)THEN
             arg=-8._dp * (frac_daylight-0.36_dp)**2 + 0.685_dp
           ELSE
             arg=-0.8_dp * (frac_daylight-0.36_dp)**2 + 0.685_dp
           ENDIF
         ELSE
           arg = 0.5_dp
         ENDIF

         dtmp = arg * dtmp_max(i) + (1._dp-arg) * dtmp_min(i) ! The dusk day temperature
         !dtmp = 0.5_dp*(dtmp_max(i)+dtmp_min(i)) ! The mean day temperature
         dtran = dtmp_max(i)- dtmp_min(i)        ! temperature range
         IF (frac_daylight>=rat2 .AND. frac_daylight<=rat3) THEN ! -------- normal day ---------------------------------------------
            frac_sunrise = 0.5_dp - frac_daylight/2._dp             ! time of sunrise expressed as fraction of day since midnight
            frac_sunset = frac_sunrise+frac_daylight                ! time of sunset expressed as fraction of day since midnight 
            IF (frac_localtime > frac_sunrise .AND. frac_localtime <= frac_sunset) THEN  ! For time steps at daylight ...
               sd = COS(pi*(frac_localtime-rat1)/(frac_daylight/2._dp+rat4)) ! ... modulation factor for day-temperature such that
                                                                    ! .. temperature is at maximum at rat1 (2 pm) and at average at
                                                                    ! .. sunrise-rat4/2 and sunset+rat4/2
               IF (sd >= 0) THEN
                   tmp(i) = dtmp + sd*(dtmp_max(i)-dtmp)                       ! temperature modulated by temperature range
               ELSE 
                   tmp(i) = dtmp + sd*(dtmp-dtmp_min(i))                       ! temperature modulated by temperature range
               ENDIF
               !tmp(i) = dtmp + sd*dtran/2._dp                       ! temperature modulated by temperature range
            ELSE                                                 ! For time steps at night ...
               sd = COS(pi*(frac_sunset-rat1)/(frac_daylight/2._dp+rat4)) ! modulation factor for temperature at sunset
               IF (sd >= 0) THEN
                   temp_sunset = dtmp + sd*(dtmp_max(i)-dtmp)                       ! temperature modulated by temperature range
               ELSE 
                   temp_sunset = dtmp + sd*(dtmp-dtmp_min(i))                       ! temperature modulated by temperature range
               ENDIF
               !temp_sunset = dtmp + dtran/2._dp*sd                     ! temperature at sunset
               frac_since_sunset = MOD(frac_localtime-frac_sunset+1._dp,1._dp) ! day fraction since sunset
               tmp(i) = dtmp_min(i)+(temp_sunset-dtmp_min(i))*(1._dp-frac_since_sunset/(1._dp-frac_daylight))! linear interpolation
                                                                       ! ... of temperature such that temperature is equal to the 
                                                                       ! ... temperature at sunset and drops until the minimum ...
                                                                       ! ... temperature is met at sunrise.
            END IF
         ELSEIF (frac_daylight>rat3) THEN ! --- long day (close to polar day) ------------------------------------------------------
           ! case is now obsolete
            sd = COS(pi*(frac_localtime-rat1)/(frac_daylight/2._dp+rat4)) ! modulation factor for day-temperature such that ...
                                                                       ! ... temperature is at maximum at rat1 (2 pm) and at average
                                                                       ! ... at sunrise-rat4/2 and sunset+rat4/2
            tmp(i) = dtmp + dtran/2._dp*sd                             ! temperature modulated by temperature range
         ELSE ! ------------------------------- short day (close to polar night) ---------------------------------------------------
            tmp(i) = dtmp                                        ! temperature is set to mean temperture at all time steps
         END IF
      ENDDO
      
      RETURN
    END SUBROUTINE INSTANTLY_FROM_DAILY_TEMP_2

    ! === precipitation_generator   ================================================================================================
    !
    ! This subroutine generates stochastically values for daily precipitation  The method consists in first determining by a first
    ! order Markov-chain whether the considered day is a wet or a dry day, and second, in case the day is a wet day, in determining
    ! the amount of precipitation from a Gamma-distribution. In principle the model has four parameters: the transition probabili-
    ! ties from a wet to a wet day p(W|W) and from a dry to wet day p(W|D), plus the two parameters alpha and beta of the Gamma-
    ! distribution. As was shown by Geng et al. all these parameters can be estimated from the fraction of wet days in the parti-
    ! cular month of the considered day (frac_wet_days) and the average precipitation on a wet day (precip_wetday_mean), which are
    ! the only inputs of the routine.
    !
    ! The implementation follows
    !     (1) S. Geng, F.W.T. Penning de Vries, I. Supit, "A simple method for generating daily rainfall data",
    !         Agricultural and Forest Meteorology, 36 (1986) 363-376
    !     (?) G.A. Larsen and R.B. Pense, "Stochastic simulation of daily climatic data for agronomic models",
    !         Agronomy Journal, 74 (1982) 510-514.
    !     (?) C.W. Richardson, "Stochastic simulation of daily precipitation, temperature, and solar-radiation"
    !        Water Resources Research, 17 (1981) 182-190.
    !
    ! ATTENTION: 
    !     This routine should be called only at the beginning of each month!
    !
    ! last modified by Thomas Raddatz 01.2004; preserving monthly mean value
    !
    SUBROUTINE precipitation_generator(precip_daily_mean,frac_wet_days,date,preserve_monthly_precip,precip_out)
      USE mo_time_control, ONLY: get_date_components,get_month_len
      IMPLICIT NONE
      INTEGER, PARAMETER :: ia = 16807, im = 2147483647, iq = 127773
      INTEGER, PARAMETER :: ir = 2836
      REAL(dp), PARAMETER:: am = 1._dp/im

      REAL(dp), INTENT (IN)  :: precip_daily_mean(domain%nland) ! mean daily precipitation in the considered month [kg/m^2/s]
      REAL(dp), INTENT (IN)  :: frac_wet_days(domain%nland)     ! average fraction of wet days in the considered month
      TYPE (time_days)   :: date                                ! the considered day in days-seconds format
      LOGICAL, INTENT (IN)   :: preserve_monthly_precip         ! determines whether the output values are scaled at the end of the
                                                                ! routine to preserve the monthly values of the input data
      REAL(dp), INTENT (OUT) :: precip_out(domain%nland,31)     ! the stochastically generated daily precipitation for a whole
                                                                ! month [kg/m^2/s]
      REAL(dp)  :: precip_sum(domain%nland)              ! array to sum up the precipitation within the current month

      LOGICAL,SAVE :: initialized=.FALSE.                ! to remember whether routine was intialized
      INTEGER,SAVE,ALLOCATABLE :: wet_state_yesterday(:) ! stores last days wet state (1:wet, 0:dry)
      INTEGER,SAVE :: idum                               ! a quasi random number 

      REAL(dp) :: ppw     ! average precipitation per wet day [mm/wet day]
      REAL(dp) :: pwd     ! transition probality for a wet day following a dry day p(W|D)
      REAL(dp) :: pww     ! transition probality for a wet day following a wet day p(W|W)
      REAL(dp) :: pwet    ! probability for the current day to be wet
      REAL(dp) :: x       ! random value from [0,1]
      REAL(dp) :: alpha   ! Gamma-distribution parameter, usually called gamma! [unitless]
      REAL(dp) :: beta    ! Gamma-distribution parameter [mm/day]  (<x>=alpha*beta, s^2=alpha*beta^2)
      REAL(dp) :: c, rr   ! helpers
      INTEGER :: i, j, k ! helpers
      INTEGER :: status ! for allocation status
      INTEGER :: year,month,month_length

      IF(.NOT. initialized) THEN ! quasirandom initialization of wet_state_yesterday and idum
         ALLOCATE(wet_state_yesterday(domain%nland),STAT=status)
         IF(status /= 0) CALL finish("precipitation_generator(", "ERROR: Field wet_state_yesterday() could not be allocated.")
         !idum = 10248098
         idum = 13458010
         DO i = 1,domain%nland
            wet_state_yesterday(i) = 0
            k = idum/iq
            idum = ia*(idum-k*iq) - ir*k
            IF (idum<0) idum = idum + im
            IF (am*idum<=frac_wet_days(i)) wet_state_yesterday(i) = 1
         END DO
         initialized = .TRUE.
      END IF
      CALL get_date_components(date,YEAR=year,MONTH=month)   ! extract year and month
                                                             ! --> WRONG FOR PALAEO COMPUTATIONS (s. daily_from_monthly)!!!
      month_length = get_month_len (year, month)             ! the length of the current month
      precip_sum(:) = 0._dp
      DO j = 1,month_length      
      DO i = 1,domain%nland
         ! Determine probability of a wet day from input values by the method of Geng et al.
         pwd = 0.75_dp*frac_wet_days(i) ! transition probability p(W|D) that a wet day is followed by a dry day is proportional ...
                                        ! ... to the fraction of wet days during one month (Eq. (3))
         pww = 0.25_dp + pwd            ! transition probability p(W|W) that a wet day is followed by a wet day (Eq. (4))
         pwet = pww*wet_state_yesterday(i) + pwd*REAL(1-wet_state_yesterday(i),dp) ! probability for a wet day today
         
         ! compute random number idum, so that am*idum is evenly distributed in [0,1]
         k = idum/iq
         idum = ia*(idum-k*iq) - ir*k
         IF (idum<0) idum = idum + im
         
         ! decide whether day is a wet day and eventually compute precipitation from Gamma-distribution
         IF (am*idum<=pwet .AND. precip_daily_mean(i)>0._dp) THEN     ! if wet day ...
            ppw = precip_daily_mean(i)/MAX(frac_wet_days(i),0.033_dp) ! ... compute average precipitation per wet day from mean ...
                                                                      ! ... daily precipitation and assume that precipitation ...
                                                                      ! ... at a wet day is at most 30 times the daily mean.
            ppw = ppw/conv_prec_day          ! convert precipitation from [kg/m^2/s] to [mm/day]
            IF(ppw >= 3._dp) THEN
               beta = -2.16_dp + 1.83_dp*ppw ! compute beta of the Gamma-distribution from mean precipitation at wet day (Eq. (5))
            ELSE
               beta = 1.11_dp*ppw            ! ... but take a corrected value for low precipitation (beta should stay positive!)
            END IF
            alpha = ppw/beta
            c = euler / (alpha+euler)

            ! determine gamma distribution by an appropriate random process
            gammaDistibution: DO    
               k = idum/iq
               idum = ia*(idum-k*iq) - ir*k
               IF (idum<0) idum = idum + im
               x = am*idum
               
               IF (x<=c) THEN
                  x = (x/c)**(1._dp/alpha)
                  rr = EXP(-x)
               ELSE
                  x = -LOG((x-c)/(1._dp-c)/euler)
                  rr = x**(alpha-1._dp)
               END IF
               
               k = idum/iq
               idum = ia*(idum-k*iq) - ir*k
               IF (idum<0) idum = idum + im
               
               IF (am*idum < rr) EXIT gammaDistibution
            END DO gammaDistibution
            
            precip_out(i,j) = beta*x                            ! final value for precipitation in [mm/day]
            precip_sum(i) = precip_sum(i) + precip_out(i,j)     ! sum up the precipitation within the month
            precip_out(i,j) = conv_prec_day * precip_out(i,j)   ! final value for precipitation in [kg/m^2/s]
            wet_state_yesterday(i) = 1                          ! remember wet state
            
         ELSE                                                   ! dry day
            precip_out(i,j) = 0._dp                             ! final value for precipitation
            wet_state_yesterday(i) = 0                          ! remember wet state
         END IF
      END DO
      END DO
      IF (preserve_monthly_precip) THEN
         ! set precip_sum to [kg/m^2/s]
         precip_sum(:) = conv_prec_day * precip_sum(:) / REAL(month_length,dp)
         ! set whole monthly precipitation to the 15th of the month, if precip_out is zero throughout the month
         DO i = 1,domain%nland
            IF (precip_sum(i) < negligible .AND. precip_daily_mean(i) > negligible) THEN
               precip_out(i,:) = 0._dp
               precip_out(i,15) = precip_daily_mean(i) * REAL(month_length,dp)
               precip_sum(i) = precip_daily_mean(i)
            END IF
         END DO
         ! preserve monthly mean precipitation rate
         DO j = 1,month_length
            precip_out(:,j) =  precip_out(:,j) * precip_daily_mean(:) / MAX(precip_sum(:),negligible)
         END DO
      END IF
      RETURN
    END SUBROUTINE PRECIPITATION_GENERATOR

    ! === relhum_generator   ================================================================================================
    !
    ! This subroutine generates stochastically values for daily relative humidity.  The method consists of taking into account
    ! dry and wet days (generated from precipitation generator, and second, in determining
    ! the amount of relative humidity from a Gamma-distribution. Implemented as in EPIC
    !
    ! ATTENTION: 
    !     This routine should be called only at the beginning of each month!
    !
    ! last modified by Soenke Zaehle 08.2010
    !
    SUBROUTINE relhum_generator(rel_hum_daily_mean,frac_wet_days,precip_daily,date,rel_hum_out)
      USE mo_time_control, ONLY: get_date_components,get_month_len
      IMPLICIT NONE

      REAL(dp), PARAMETER:: omqd = 0.5 ! from EPIC weather generator

      REAL(dp), INTENT (IN)  :: rel_hum_daily_mean(domain%nland) ! mean daily air humidity in the considered month [fraction]
      REAL(dp), INTENT (IN)  :: frac_wet_days(domain%nland)      ! average fraction of wet days in the considered month
      REAL(dp), INTENT (IN)  :: precip_daily(domain%nland,31)    ! daily precipitation in the considered month
      TYPE (time_days)   :: date                                 ! the considered day in days-seconds format
      REAL(dp), INTENT (OUT) :: rel_hum_out(domain%nland,31)     ! the stochastically generated relative humidity for a whole month
                                                                 ! [fraction]
      REAL(dp),DIMENSION(domain%nland) :: qdd,qdw           ! dry and wet day relative humidities [fraction]
      REAL(dp) :: qde,qd,qdup,qdlow,amn ! helpers
      REAL(dp) :: xx,y,x1 ! helpers
      REAL(dp) :: b1,b2,b3   ! parameters of the triangular distribution

      INTEGER :: i, j ! helpers
      INTEGER :: year,month,month_length

      REAL(dp) :: rn(domain%nland,31)

      INTEGER :: iseed1 = 13458010
      INTEGER :: iseed2 = 14352334
      INTEGER, ALLOCATABLE :: seed(:)
      INTEGER :: nseed

      !find out minimum size of the random_seeds
      CALL RANDOM_SEED(size=nseed)
      ALLOCATE(seed(nseed))
      seed(:) = 0
      seed(1) = iseed1
      seed(2) = iseed2

      !initialise and obtain pseudo random numbers [0..1]
      CALL RANDOM_SEED(put=seed)
      CALL RANDOM_NUMBER(rn)
      DEALLOCATE(seed)

      ! extract year and month
      ! --> WRONG FOR PALAEO COMPUTATIONS (s. daily_from_monthly)!!
      CALL get_date_components(date,YEAR=year,MONTH=month)   
      month_length = get_month_len (year, month)             

      !--- calculate mean dry and wet day relative humidity
      rel_hum_out = 0._dp
      DO i = 1,domain%nland
         if (frac_wet_days(i) > 0._dp) then
            qdd(i) = (rel_hum_daily_mean(i)-frac_wet_days(i)*omqd) / & 
                 (1.0-frac_wet_days(i)*omqd)
            if (qdd(i) < 0.2) then
               qdd(i) = 0.2
            endif   
            qdd(i) = MIN(1._dp,qdd(i))
            qdw(i) = (rel_hum_daily_mean(i)-(1.0-frac_wet_days(i))*qdd(i)) / &
                 frac_wet_days(i)
         else
            qdd(i) = rel_hum_daily_mean(i)
            qdw(i) = rel_hum_daily_mean(i)
         endif
      ENDDO

      DO j = 1,month_length      
         DO i = 1,domain%nland
     
            if (precip_daily(i,j) <= negligible ) then
               qde = qdd(i) ! dry day
            else
               qde = qdw(i) ! wet day
            endif
            !-----
            !---- estimate lower and upper bounds of humidity distribution function
            !---- following logic of the EPIC weather generator code
            !-----
            xx = exp(qde)
            qdup  = qde+(1._dp-qde)*xx/exp(1._dp)
            qdlow = qde*(1._dp-1._dp/xx)
            !-----
            !---- randomly select humidity from triangular distribution function
            !---- following logic of the EPIC weather generator code
            !-----
            y  = 2._dp/(qdup-qdlow)
            !-----
            b3 = qde-qdlow
            b2 = qdup-qde
            b1 = rn(i,j)/y
            !-----
            x1 = y*b3/2._dp
            !-----
            if (rn(i,j) > x1) then
               qd = qdup-sqrt (b2*b2-2._dp*b2*(b1-0._dp*b3))
            ELSE
               qd = qdlow+sqrt (2._dp*b1*b3)
            ENDIF
            !-----
            !---- adjust daily humidity to conserve monthly mean values
            !-----
            !---- note that this adjustment sometimes gives rise to humidity
            !---- values greater than 1.0 -- which is corrected below
            !-----
            amn = (qdup+qde+qdlow)/3.0
            qd  = qd*qde/amn
            !-----
            !---- constrain daily average relative humidity
            !-----
            rel_hum_out(i,j) = MAX(0.15_dp,MIN(qd,0.99_dp))

         ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RELHUM_GENERATOR

    SUBROUTINE spec_humidity_from_rel_humidity(rel_humidity,Tmin,Tmax,elevation,delta_time,qair)

      ! This subroutine converts daily relative humidity to specific humidity considering the daily course of temperature
      ! ops! Using this routine only makes sense together this the generated diurnal cycle of temperature as in 
      ! SUBROUTINE instantly_from_daily_temp_2 
      ! Otherwise the results will be faulty!

      USE mo_atmosphere,   ONLY: sat_specific_humidity
      USE mo_orbit,           ONLY: inquire_declination
      
      IMPLICIT NONE
      
      REAL(dp),         INTENT(in) :: rel_humidity(domain%nland)  !! mean daily relative humidity
      REAL(dp),         INTENT(in) :: Tmin(domain%nland)       !! minimum temperature (= temperature at sunrise)
      REAL(dp),         INTENT(in) :: Tmax(domain%nland)       !! maximum temperature (= at 2pm)
      REAL(dp),         INTENT(in) :: elevation(domain%nland)  !! height above sea leavel
      REAL(dp),         INTENT(in) :: delta_time !! time step of model
      REAL(dp),         INTENT(out) :: qair(domain%nland) !! the specific humidity of the day
      
      ! constants
      REAL(dp), PARAMETER :: pio180 = pi/180._dp    ! conversion factor from degrees to radians ("pi over 180")
      REAL(dp), PARAMETER :: rat1=14._dp/24._dp      ! day fraction at which maximum temperature shall occur
      REAL(dp), PARAMETER :: rat2=4._dp/24._dp       ! minimum length of normal day (as fraction of day)
      REAL(dp), PARAMETER :: rat3=20._dp/24._dp      ! maximum length of normal day (as fraction of day)
      REAL(dp), PARAMETER :: rat4=2._dp/24._dp       ! day fraction of dawn times (before sunrise + after sunset), used to set ...
                                                     ! ... temperature of a normal day to T_mean at sunrise-rat4/2 and sunset+rat4/2

      ! Temperature calculations
      REAL(dp)   :: declination       ! solar declination
      REAL(dp)   :: cpds, spds        ! cosine and sine contributions to zenith angle
      REAL(dp)   :: frac_daylight     ! fraction of day between sunrise and sunset
      REAL(dp)   :: frac_since_sunset ! time between current time step and sunset expressed as fraction of day
      REAL(dp)   :: frac_sunrise      ! time of sunrise expressed as fraction of day since midnight
      REAL(dp)   :: frac_sunset       ! time of sunset expressed as fraction of day since midnight
      REAL(dp)   :: temp_sunset       ! temperature at sunset
      REAL(dp)   :: dtmp              ! mean day temperature
      REAL(dp)   :: dtran             ! day temperature range: maximum - minimum
      REAL(dp)   :: sd                ! modulation factor for temperature
      REAL(dp)   :: arg,nstep
      REAL(dp),ALLOCATABLE :: tmp(:)  ! diurnal temperature 
      
      ! Qair calculations
      REAL(dp)   :: rh_test,rh_last,qsat_last ! temporary RH estimates
      REAL(dp),ALLOCATABLE :: qsat(:)      ! diurnal saturated specific humidity
      REAL(dp),ALLOCATABLE :: inv_qsat(:)  ! inverse diurnal saturated specific humidity
      REAL(dp),ALLOCATABLE :: qair_test(:) ! diurnal actual specific humidity, temporary
      REAL(dp),ALLOCATABLE :: pressure(:)  ! diurnal air pressure
      REAL(dp)   :: tstep        ! fractional time step of the day 
      
      ! helpers
      INTEGER    :: i,j,jj,time_steps_per_day
      INTEGER    :: offset(1)

      CALL inquire_declination(declination)

      nstep = (86400.0_dp/delta_time)
      time_steps_per_day = 86400/int(delta_time)

      IF(.NOT.(ALLOCATED(tmp))) ALLOCATE(tmp(time_steps_per_day))
      IF(.NOT.(ALLOCATED(pressure))) ALLOCATE(pressure(time_steps_per_day))
      IF(.NOT.(ALLOCATED(qsat))) ALLOCATE(qsat(time_steps_per_day))
      IF(.NOT.(ALLOCATED(inv_qsat))) ALLOCATE(inv_qsat(time_steps_per_day))
      IF(.NOT.(ALLOCATED(qair_test))) ALLOCATE(qair_test(time_steps_per_day))

    DO i = 1,domain%nland
       spds = domain%sinlat(i)*SIN(declination)  ! sine contribution to solar zenith angle
       cpds = domain%coslat(i)*COS(declination)  ! cosine contribution to solar zenith angle
       arg = -spds/cpds                        ! faulty?: for cpds=0 arg gets infinite
       IF (arg > 1._dp) THEN               ! polar night:
          frac_daylight = 0._dp            ! ... day length is zero
       ELSE IF (arg < -1._dp) THEN         ! polar day:
          frac_daylight = 1._dp            ! ... day length is the whole day
       ELSE                                ! normal day / night:
          frac_daylight = ACOS(arg)/pi    ! ... day length expressed as fraction of the whole day
       END IF
       dtmp = 0.5_dp*(Tmax(i)+Tmin(i)) ! The mean day temperature
       dtran = Tmax(i)- Tmin(i)        ! temperature range
       DO j = 1,time_steps_per_day
          tstep=j/nstep
          ! -------- normal day ---------------------------------------------
          IF (frac_daylight>=rat2 .AND. frac_daylight<=rat3) THEN 
             frac_sunrise = 0.5_dp - frac_daylight/2._dp             ! time of sunrise expressed as fraction of day since midnight
             frac_sunset = frac_sunrise+frac_daylight                ! time of sunset expressed as fraction of day since midnight 

             ! For time steps at daylight ...
             IF (tstep>frac_sunrise .AND. tstep<frac_sunset) THEN  
                sd = COS(pi*(tstep-rat1)/(frac_daylight/2._dp+rat4)) ! modulation factor for day-temperature such that
                                                                      !   temperature is at maximum at rat1 (2 pm) and at average at
                                                                      !   sunrise-rat4/2 and sunset+rat4/2
                tmp(j) = dtmp + sd*dtran/2._dp                        ! temperature modulated by temperature range
 
            ! For time steps at night ...
             ELSE                                                 
                sd = COS(pi*(frac_sunset-rat1)/(frac_daylight/2._dp+rat4)) ! modulation factor for temperature at sunset
                temp_sunset = dtmp + dtran/2._dp*sd                    ! temperature at sunset
                frac_since_sunset = MOD(tstep-frac_sunset+1._dp,1._dp) ! day fraction since sunset
                tmp(j) = Tmin(i)+(temp_sunset-Tmin(i))*(1._dp-frac_since_sunset/(1._dp-frac_daylight))! linear interpolation
                                                                       ! ... of temperature such that temperature is equal to the 
                                                                       ! ... temperature at sunset and drops until the minimum ...
                                                                       ! ... temperature is met at sunrise.
             END IF
          ! --- long day (close to polar day) ------------------------------------------------------
          ELSEIF (frac_daylight>rat3) THEN 
             sd = COS(pi*(tstep-rat1)/(frac_daylight/2._dp+rat4)) ! modulation factor for day-temperature such that ...
                                                                       ! ... temperature is at maximum at rat1 (2 pm) and at average
                                                                       ! ... at sunrise-rat4/2 and sunset+rat4/2
             tmp(j) = dtmp + dtran/2._dp*sd                             ! temperature modulated by temperature range

          ! ------------------------------- short day (close to polar night) ---------------------------------------------------
          ELSE 
             tmp(j) = dtmp                                        ! temperature is set to mean temperture at all time steps
          END IF
       ENDDO
 
       ! from function bottom pressure
       pressure(:) = &
            p_sealevel * ( 1._dp - elevation(i)*gamma /(Tmelt+tmp(:) + elevation(i)*gamma))**(Gravity/GasConstantDryAir/gamma)

       qsat(:)=sat_specific_humidity(tmp(:)+Tmelt, pressure(:)) ! present saturated qair    

       ! implied RH from minimum saturated specific humidity
       rh_test=SUM(MINVAL(qsat(:))/qsat(:))/nstep
       IF( rel_humidity(i) <= rh_test ) THEN
          ! Simply case, relative humitidy that assuming a constant qair there is never 
          ! saturation during the day. In this case, qair can be inferred by inverting
          ! RH = 1/n * sum ( QAIR / QSAT ) = QAIR / n * sum ( 1 / QSAT ) 
          ! -> QAIR = n * RH / sum ( 1 / QSAT )
          inv_qsat(:) = 1._dp / qsat(:) 
          qair(i) = nstep * rel_humidity(i) / SUM( inv_qsat(:) )
       ELSE
          ! The case where RH implies saturation at night. Under these conditions
          ! search for the qair value which is closed to giving the prescribed RH
          ! taking account of the implied period during the day where saturation occurs.
          ! Surely there's a better way to do this, but I cannot be bothered now.
          offset(:)=MAXLOC(qsat(:))
          IF( rel_humidity(i) >= 1._dp ) THEN
             qair(i)=qsat(offset(1))
          ELSE
             j=0
             jj=offset(1)
             rh_test=1._dp
             DO WHILE ( rh_test > rel_humidity(i) .AND. j < time_steps_per_day ) 
                ! save value from last time
                qsat_last=qsat(jj)
                rh_last=rh_test

                ! find tome step after hottest time of the day
                IF(j+offset(1)>time_steps_per_day) THEN
                   jj=j+offset(1)-time_steps_per_day
                ELSE
                   jj=j+offset(1)
                ENDIF

                ! create qair values constrained by qsat at low temperatures
                WHERE(qsat(:) <= qsat(jj))
                   qair_test(:)=qsat(:)
                ELSEWHERE
                   qair_test(:)=qsat(jj)
                ENDWHERE

                ! test whether target RH is reached
                rh_test=SUM(qair_test(:)/qsat(:))/nstep
                j=j+1
             ENDDO
             ! The above solution is inprecise. Attempt to reconcile this by linear approximation
             qair(i)=qsat(jj)+(qsat_last-qsat(jj))*(rel_humidity(i)-rh_test)/(rh_last-rh_test)
          ENDIF
       ENDIF
    ENDDO

    IF((ALLOCATED(tmp))) DEALLOCATE(tmp)
    IF((ALLOCATED(pressure))) DEALLOCATE(pressure)
    IF((ALLOCATED(qsat))) DEALLOCATE(qsat)
    IF((ALLOCATED(inv_qsat))) DEALLOCATE(inv_qsat)
    IF((ALLOCATED(qair_test))) DEALLOCATE(qair_test)

  END SUBROUTINE spec_humidity_from_rel_humidity

  END SUBROUTINE update_forcing

  ! === finish_forcing   ===========================================================================================================
  !
  !
  ! This subroutine deallocates memory that initially was allocated in init_forcing. Should be called when ending the program.
  !
  SUBROUTINE finish_forcing(options, forcingFields)

    TYPE(options_type), INTENT(in)      :: options 
    TYPE(forcing_type), INTENT(inout)   :: forcingFields   !! This routine deallocates memory for this structure 

    ! --- general ---
    DEALLOCATE(hlp_field2D_local)
    IF(p_parallel_io.AND.ALLOCATED(hlp_field2D_global)) DEALLOCATE(hlp_field2D_global) 
    IF(p_parallel_io.AND.ALLOCATED(hlp_field3D_global)) DEALLOCATE(hlp_field3D_global)

    IF (options%ReadInterfaceVars) RETURN

    ! --- temperature ---
    
    IF( (forcing_input%air_temp_min%initialized) ) THEN ! If one of them is present, also the other has been called
      CALL finish_variable_access(forcing_input%air_temp_min)
      CALL finish_variable_access(forcing_input%air_temp_max)
    ENDIF
    IF( (forcing_input%air_temp%initialized) ) THEN
      CALL finish_variable_access(forcing_input%air_temp) 
    ENDIF
    DEALLOCATE(forcingFields%air_temp)

    ! --- precipitation ---

    CALL finish_variable_access(forcing_input%precipitation)
    IF(forcing_input%frac_wet_days%initialized) CALL finish_variable_access(forcing_input%frac_wet_days)
    DEALLOCATE(forcingFields%precip_rain,forcingFields%precip_snow)

    ! --- radiation ---

    DEALLOCATE(forcingFields%rad_sw_down,forcingFields%rad_UV_down,forcingFields%rad_PAR_down,forcingFields%rad_NIR_down)
    DEALLOCATE(forcingFields%frac_PAR_diffuse)
    SELECT CASE(forcing_options%type_of_shortwave_forcing)
       CASE(CLOUD_)
          CALL finish_variable_access(forcing_input%cloud_cover)
          DEALLOCATE(forcingFields%rad_lw_down)
       CASE(MEAN_RAD_)
          CALL finish_variable_access(forcing_input%shortwave)
          CALL finish_variable_access(forcing_input%longwave)
       CASE default
          ! nothing to do
    END SELECT

    ! --- atmospheric humidity ---

    IF( (forcing_input%qair%initialized) ) THEN 
       CALL finish_variable_access(forcing_input%qair)
    ENDIF  
    IF( (forcing_input%rel_humidity%initialized) ) THEN 
       CALL finish_variable_access(forcing_input%rel_humidity)
    ENDIF  
  
    ! --- CO2 ----

    IF (.NOT. forcing_input%CO2_concentr%initialized) &
       CALL finish_variable_access(forcing_input%CO2_concentr)

    ! --- wind speed ---

    CALL finish_variable_access(forcing_input%wind_speed)

    ! --- remaining forcing fields

    DEALLOCATE(forcingFields%air_pressure,forcingFields%vapor_pressure,forcingFields%spec_humidity)

  END SUBROUTINE finish_forcing

  ! === bottom pressure() =======================================================================================================
  !
  ! To estimate bottom pressure a so called "polytrope atmosphere" is assumed for the troposhere. Here a fixed temperature
  ! gradient "gamma" is prescribed, i.e. air temperature at height z is set to 
  !     (1)     T(z)=T0-gamma*z, 
  ! where T0 is the temperature at the bottom. Assuming further hydrostatic conditions the vertical pressure gradient is given by
  !              dp
  !     (2)     ---- = -g*rho,
  !              dz
  ! where rho is the air density and g the acceleration due to gravity. Pressure and Temperature are related by the state 
  ! equation of an ideal gas
  !     (3)      p = rho*R*T,
  ! where R is the specific gas constant and T the temperature. By eliminating rho from (3) and (2) and entering (1) for T one 
  ! obtains
  !              dp        p*g
  !     (4)     ---- = - --------------
  !              dz      R(T0-gamma*z)
  ! By integration one obtains the final result
  !                              gamma                              g
  !     (5)      p = p0 * (1 - z*-----)^k , with the exponent k= ------- ,
  !                               T0                             R*gamma
  ! where p0 and T0 are the bottom pressure and temperature. It should be noted that the expression in brackets can get negative 
  ! only for extreme values of T0 and z. Assuming e.g. T0 > -80 Celsius (approx. 200 Kelvin), the bracket stays positive for 
  ! z<30km, which is true even at the Mount Everest. (Eq. (5) coincides with Eq. (31) from Knorr for z << 1!)
  !
  ! If temperature T(z) at height z is known instead of T0, using Eq. (1) Eq. (5) translates to:
  !                                  gamma      
  !     (6)      p = p0 * (1 - z*------------)^k
  !                              T(z)+z*gamma   
  !
  
  PURE elemental FUNCTION bottom_pressure(elevation,air_temperature)
    REAL(dp)              :: bottom_pressure
    REAL(dp), INTENT(in)  :: elevation
    REAL(dp), INTENT(in)  :: air_temperature

    bottom_pressure = &
        p_sealevel * ( 1._dp - elevation*gamma /(Tmelt+air_temperature + elevation*gamma))**(Gravity/GasConstantDryAir/gamma)
  END FUNCTION bottom_pressure


  ! === snow_and_rain_from_precip() ==============================================================================================
  !
  ! Precipitation is divided into rain and snow according to a method by M.S. Wigmosta, L. Vail and D.P. Lettenmaier, "A 
  ! distributed hydrology-vegetation model for complex terrain", Water Resources Research 30 (1994) 1665-1679: depending on 
  ! average day temperature "T", the snow part of precipitation "P_sn" is computed from the total precipitation "P" as:
  !               / P               for T < -1.1 Celsius
  ! (1)   P_sn = |  P*(3.3 - T)/4.4 for -1.1 Celsius <= T <= 3.3 Celsius
  !               \ 0               for T > 3.3 Celsius
  ! Hence the rain part ist then
  ! (2)   P_rain = P-P_sn.
  !  
  PURE ELEMENTAL SUBROUTINE snow_and_rain_from_precip(rain,snow,precipitation,air_temp_daily)
    REAL(dp),INTENT(out) :: rain          !! rain [kg/m^2/s]
    REAL(dp),INTENT(out) :: snow          !! snow [kg/m^2/s]
    REAL(dp),INTENT(in)  :: precipitation !! total precipitation (rain+snow) [kg/m^2/s]
    REAL(dp),INTENT(in)  :: air_temp_daily!! daily mean air_temperature at surface [Celsius]

    IF(air_temp_daily < -1.1_dp) THEN
       snow = precipitation
    ELSE IF(air_temp_daily > 3.3_dp) THEN
       snow = 0.0_dp 
    ELSE
       snow = precipitation*(3.3_dp-air_temp_daily)/4.4_dp
    END IF
    ! set the rain part according to Eq. (2):
    rain = precipitation - snow
  END SUBROUTINE snow_and_rain_from_precip

  ! === shortwave_from_cloudcover() ================================================================================================
  !
  ! Following the strategy described in [1], this routine computes PAR and NIR, together with the diffuse part of PAR, from 
  ! cloudcover, by converting cloudcover "n" to "fpar" (fpar is the fraction of actual to potential radiation flux in the PAR 
  ! band incident at the earth's surface): 
  ! (1)           fpar = 0.5 + 0.4 * n.
  ! Then the main computations are done in the subroutine "shortwave_from_fpar()" (see there for a description of the method).
  !
  ! Literature:
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49, 
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  !
  ELEMENTAL SUBROUTINE shortwave_from_cloudcover(day_of_year, cos_zenith, air_pressure, cloud_cover,                    & !! input
                                                 rad_UV_down, rad_PAR_down, rad_NIR_down, frac_PAR_diffuse, rad_sw_down,& !! output
                                                 rad_sw_down_pot)                                                         !! output
    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(dp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(dp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(dp),INTENT(in)  :: cloud_cover      !! fraction of sky covered by clouds (from [0,1]) 
    REAL(dp),INTENT(out) :: rad_UV_down      !! solar radiation flux (direct+diffuse) from the UV band 280-400 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: rad_PAR_down     !! solar radiation flux (direct+diffuse) from the visible band 400-700 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: rad_NIR_down     !! solar radiation flux (direct+diffuse) from the near infrared band 700-3000 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: frac_PAR_diffuse !! Fraction from  rad_PAR_down() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(dp),INTENT(out) :: rad_sw_down      !! total solar shortwave radiation : direct+diffuse, 280-3000 nm, (i.e. PAR+NIR), ..
                                             !! ... that actually arrives at the bottom [W/m^2]
    REAL(dp),INTENT(out) :: rad_sw_down_pot  !! total solar shortwave radiation : direct+diffuse, 280-3000 nm, (i.e. PAR+NIR), ..
                                             !! ... that actually arrives at the bottom [W/m^2]
 
    ! locals 

    REAL(dp) :: fpar

    ! Go ...
    
    fpar = 0.5_dp + 0.4_dp * (1.0_dp - cloud_cover)

    SELECT CASE(forcing_options%type_of_shortwave_scheme)
       CASE(SW_SCHEME_ORIGINAL_) 
          CALL shortwave_from_fpar_orig(day_of_year, cos_zenith, air_pressure, fpar,                            & !! inputs
                                        rad_UV_down, rad_PAR_down, rad_NIR_down, frac_PAR_diffuse, rad_sw_down  & !! outputs
                                        ) 
       CASE(SW_SCHEME_CCDAS_)
          CALL shortwave_from_fpar_ccdas(day_of_year, cos_zenith, air_pressure, fpar,                           & !! inputs
                                        rad_UV_down, rad_PAR_down, rad_NIR_down, frac_PAR_diffuse, rad_sw_down  & !! outputs
                                        ) 
    END SELECT

    rad_sw_down_pot = rad_sw_down / fpar

  END SUBROUTINE shortwave_from_cloudcover

  ! === shortwave_from_fpar_orig() ===============================================================================================
  !
  ! This routine follows exactly [1]. It computes PAR and NIR, together with the diffuse part of PAR, from "fpar", which is the 
  ! quotient of the actual PAR flux "Rpar_act" (i.e. at cloudy sky) to the potential PAR flux "Rpar_pot" (i.e. at clear sky), both 
  ! measured at the ground: 
  !                       Rpar_act
  ! (1)           fpot := -------- ,
  !                       Rpar_pot
  ! Following [1] the various components of the total shortwave radiation are computed as follows: The total flux in the PAR band 
  ! incident to the outer atmosphere "Rpar_top" is computed from
  ! (2)           Rpar_top = S0*cos(theta)*E,
  ! where S0=0.44*1360 W/m^2 =~= 600 W/m^2 is the PAR part of the solar constant (0.4-0.7 nm), "theta" is the zenith angle
  ! and
  ! (3)           E = (r0/r)^2 = 1.000110+0.034221 cos(psi)+0.001280*sin(psi)+0.000719*cos(2*psi)+0.000077*sin(2*psi)
  ! the inverse squared distance earth-sun "r" in units of the average distance "r0", where "psi" is the day-angle, i.e. the 
  ! day of the year expressed as an angle (see e.g. [4], Eq. (3.3), [1] Eq. (24)).
  ! The direct part of potential PAR reaching ground "Rpar_pot_dir" is according to [2] approximately given by
  ! (4)          Rpar_pot_dir = Rpar_top*exp(-0.185*(p/p0)/cos(theta)),
  ! where "p" is the current and "p0" the sealevel pressure (see [1] Eqs. (27)-(29) and [2] Eqs. (1)-(3)), and the exponential 
  ! term accounts for direct radiation reduced by Rayleigh scattering ( p/p0/cos(theta) is a measure for the mass of air along 
  ! the beam ). According to [2] the diffuse part of potential PAR ("Rpar_pot_diff") reaching ground is approximately 40% of the 
  ! difference between Rpar_top and the direct radiation (Equation (3) n [2] contains a printing error?):
  ! (5)           Rpar_pot_diff = 0.4*(Rpar_top - Rpar_pot_dir)
  ! so that potential PAR is
  ! (6)           Rpar_pot = Rpar_pot_dir + Rpar_pot_diff = Rpar_top*[0.4 + 0.6*exp(-0.185*(p/p0)/cos(theta))].
  ! By (1) the actual PAR then follows as
  ! (7)           Rpar_act = fpar * Rpar_pot.
  ! Next the near infrared (NIR) part of shortwave radiation has to be accounted for. Following [3] this is done by introducing
  ! a correction factor "F", such that the total downwards shortwave radiation "Rtot" reaching bottom is obtained from Rpar_act by
  ! (8)           Rtot = Rpar_act/F,
  ! where, using mfpar=1-fpar, the inverse factor is given by
  ! (9)           1/F = 1 + (1.185-0.437*mfpar-0.494*mfpar^2)* exp[(0.0305-0.208*mfpar+0.495*mfpar^2)/cos(theta)]
  ! Therefore the downward radiation reaching the surface in the near infrared band "Rnir" is
  ! (10)          Rnir = Rtot - Rpar_act
  ! Finally the fraction of direct radiation in PAR at the actual PAR 
  !                           Rpar_act_dir
  ! (11)          fdirPar =   ------------
  !                             Rpar_act
  ! is estimated from the direct part in potential PAR 
  ! (11)          Rpar_pot_dir = S0 * exp(-0.185*(p/p0)/cos(theta)) 
  ! total potential PAR "Rpar_pot" and "fpar" by (see Knorr Eq. (30))
  !                            /  0 for fpar < 0.2
  !                           /                   2/3
  !                          /      / 0.9 - fpar \      Rpar_pot_dir
  ! (12)          fdirPar = (  (1- (--------------)   ) ------------       for 0.2 <= fpar <= 0.9
  !                          \      \    0.7     /         Rpar_pot
  !                           \
  !                            \  1 for fpar > 0.9
  ! where the quotient showing up in (12) is
  !
  !               Rpar_pot_dir         exp(-0.185*(p/p0)/cos(theta))
  ! (13)          ------------ =   ---------------------------------------
  !                 Rpar_pot       0.4 + 0.6*exp(-0.185*(p/p0)/cos(theta))
  !
  ! It follows that the fraction of diffuse radiation "fdiffPar" in PAR is
  !
  ! (14)         fdiffPar = 1 - fdirPar.
  !
  ! REMARK: Eq. (12) is taken from [2] but contains an additional approximation: In [2] fdirPar does not depend on fpar, but on 
  ! the transmission ratio for TOTAL SHORTWAVE flux, whereas fpar is defined as the transmission ratio for PAR only (see Eq. (1)).
  ! Hence, (12) includes the additional approximation of the shortwave transmission ratio by the PAR transmission ratio (fpar).
  !
  ! Literature:
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49, 
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  ! [2] A. Weiss & J.M. Norman, "Partitioning solar radiation into direct and diffuse visible and near-infrared components",
  ! Agricultural and Forest Meteorology, 34 (1985) 205-213.
  ! [3] R.T. Pinker & I. Laszlo, "Global distribution of photosynthetically active radiation as observed from satellites",
  ! J. Climate 5 (1992) 56-65.
  ! [4] G.W. Paltridge and C.M.R. Platt, "Radiative Processes in Meteorology and Climatology" (Elsevier, Amsterdam, 1976).


  PURE ELEMENTAL SUBROUTINE shortwave_from_fpar_orig(day_of_year,cos_zenith,air_pressure,fpar, & !! inputs
                                                     rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw   & !! outputs
                                                     )

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(dp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(dp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(dp),INTENT(in)  :: fpar             !! daily mean fraction of actual to potential PAR at bottom (from [0,1])
    REAL(dp),INTENT(out) :: rad_UV           !! "RUV": solar radiation flux (direct+diffuse) from the UV band 280-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(dp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(dp),INTENT(out) :: rad_sw           !! "Rtot": total solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
 
    ! local variables

    REAL(dp) :: day_angle     !! day of year expressed as angle in radians
    REAL(dp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(dp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(dp) :: F_invers      !! factor 1/F to obtain total shortwave from visible radiation (see Eq. (9))
    REAL(dp) :: hlp_r,hlp2_r,hlp3_r,hlp4_r

    REAL(dp), PARAMETER  ::  frac_UV = 0._dp ! UV fraction of (UV + PAR)
                                             ! to be checked, if the ratio rad_PAR/rad_sw is within observed range (0.42-0.51) with this value of frac_UV

    day_angle = 2._dp*pi*(day_of_year-1._dp)/365._dp !! ATTENTION: this may be incorrect for paleo situations
    hlp_r=2._dp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
         1.000110_dp+3.4221E-2_dp*COS(day_angle)+1.280E-3_dp*SIN(day_angle)+7.19E-4_dp*COS(hlp_r)+7.7E-5_dp*SIN(hlp_r)

    IF(cos_zenith > 0.017452406_dp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
          Rpar_top   = solar_const_vis*cos_zenith*r_sun_earth_2       ! Eq. (2) from above 
          hlp3_r = EXP(-0.185_dp*(air_pressure/p_sealevel)/cos_zenith)
          hlp4_r = 0.4_dp+0.6_dp*hlp3_r
          rad_UV = fpar * Rpar_top * hlp4_r * frac_UV
          rad_PAR = fpar * Rpar_top * hlp4_r * (1._dp - frac_UV)   ! This is Rpar_act; Eq. (7) from above
          hlp_r=1._dp-fpar
          hlp2_r=hlp_r*hlp_r
          F_invers=1._dp + (1.185_dp-0.437_dp*hlp_r-0.494_dp*hlp2_r)* &
                         EXP(MIN(2.534_dp,(0.0305_dp-0.208_dp*hlp_r+0.495_dp*hlp2_r)/cos_zenith)) ! Eq. (9) from above
          rad_sw = (rad_PAR + rad_UV) *F_invers   ! This is "Rtot"; Eq. (8) from above
          rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
          IF(fpar <= 0.2_dp) THEN 
             frac_PAR_diffuse = 1.                                                     ! Eq. (12)-(15) from above
          ELSE IF(fpar >= 0.9_dp) THEN
             frac_PAR_diffuse = 0._dp                                                  ! Eq. (12)-(15) from above
          ELSE
             frac_PAR_diffuse = 1._dp - &
                  (1._dp- ((0.9_dp -fpar)/0.7_dp)**(2._dp/3._dp))*hlp3_r/hlp4_r ! Eq. (12)-(15) from above
          END IF
    ELSE ! zenith angle larger than 89 degrees
          rad_UV = 0._dp 
          rad_PAR = 0._dp
          rad_NIR = 0._dp
          frac_PAR_diffuse = 1._dp ! (could also be another value in this case)
          rad_sw  = 0._dp
    END IF
  END SUBROUTINE shortwave_from_fpar_orig

  PURE ELEMENTAL SUBROUTINE shortwave_from_fpar_sitelevel(day_of_year, cos_zenith, air_pressure, & !! input
                                                          rad_sw,                                & !! input and output
                                                          rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse& !! outputs
                                                          )

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(dp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(dp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(dp),INTENT(inout)  :: rad_sw           !! "Rtot": total solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
    REAL(dp),INTENT(out) :: rad_UV           !! "RUV": solar radiation flux (direct+diffuse) from the UV band 280-400 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(dp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
 
    ! local variables

    REAL(dp) :: day_angle     !! day of year expressed as angle in radians
    REAL(dp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(dp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(dp) :: Rpar_pot,fpar,rad_sw_pot
    REAL(dp) :: F_invers      !! factor 1/F to obtain total shortwave from visible radiation (see Eq. (9))
    REAL(dp) :: hlp_r,hlp2_r,hlp3_r,hlp4_r

    REAL(dp), PARAMETER  ::  frac_UV = 0._dp ! UV fraction of (UV + PAR)
                                             ! to be checked, if the ratio rad_PAR/rad_sw is within observed range (0.42-0.51) with this value of frac_UV

    day_angle = 2._dp*pi*(day_of_year-1._dp)/365._dp !! ATTENTION: this may be incorrect for paleo situations
    hlp_r=2._dp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
         1.000110_dp+3.4221E-2_dp*COS(day_angle)+1.280E-3_dp*SIN(day_angle)+7.19E-4_dp*COS(hlp_r)+7.7E-5_dp*SIN(hlp_r)

    IF(cos_zenith > 0.017452406_dp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
          Rpar_top   = solar_const_vis*cos_zenith*r_sun_earth_2       ! Eq. (2) from above 
          hlp3_r = EXP(-0.185_dp*(air_pressure/p_sealevel)/cos_zenith)
          hlp4_r = 0.4_dp+0.6_dp*hlp3_r
          Rpar_pot = Rpar_top*hlp4_r
          rad_sw_pot = Rpar_pot * solar_const / solar_const_vis
          !! rad_PAR = fpar*Rpar_top*hlp4_r    ! This is Rpar_act; Eq. (7) from above
          fpar = MAX(MIN(rad_sw/rad_sw_pot,1._dp),0._dp)
          hlp_r=1._dp-fpar
          hlp2_r=hlp_r*hlp_r
          F_invers=1._dp + (1.185_dp-0.437_dp*hlp_r-0.494_dp*hlp2_r)* &
                         EXP(MIN(2.534_dp,(0.0305_dp-0.208_dp*hlp_r+0.495_dp*hlp2_r)/cos_zenith)) ! Eq. (9) from above
          !! rad_sw = rad_PAR*F_invers   ! This is "Rtot"; Eq. (8) from above
          rad_UV = rad_sw * frac_UV / F_invers
          rad_PAR = rad_sw * (1._dp - frac_UV) / F_invers ! This is "Rtot"; Eq. (8) from above
          rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
          IF(fpar <= 0.2_dp) THEN 
             frac_PAR_diffuse = 1.                                                     ! Eq. (12)-(15) from above
          ELSE IF(fpar >= 0.9_dp) THEN
             frac_PAR_diffuse = 0._dp                                                  ! Eq. (12)-(15) from above
          ELSE
             frac_PAR_diffuse = 1._dp - &
                  (1._dp- ((0.9_dp -fpar)/0.7_dp)**(2._dp/3._dp))*hlp3_r/hlp4_r ! Eq. (12)-(15) from above
          END IF
    ELSE ! zenith angle larger than 89 degrees 
          rad_UV = 0._dp
          rad_PAR = 0._dp
          rad_NIR = 0._dp
          frac_PAR_diffuse = 1._dp ! (could also be another value in this case)
          rad_sw  = 0._dp
          rad_sw_pot  = 0._dp
    END IF
  END SUBROUTINE shortwave_from_fpar_sitelevel

  PURE ELEMENTAL SUBROUTINE shortwave_from_direct_shortwave(day_of_year,cos_zenith,cos_lat,sin_lat,air_pressure, & !! inputs
                                                     rad_sw_in, & !! inputs
                                                     rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw,rad_sw_pot   & !! outputs
                                                     )
    ! similar concept as above, however, the observed radiation is used to derive the remaining
    ! parameters. In daily version, the instantaneous flux at noon required to calculate irradiation at any 
    ! time of the day is approximated by solving the integral over the positive sun angle, following Prentice et al. Ecol. Mod. 1993


    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(dp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(dp),INTENT(in)  :: cos_lat          !! cosine of latitude 
    REAL(dp),INTENT(in)  :: sin_lat          !! sine of latitude 
    REAL(dp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(dp),INTENT(in)  :: rad_sw_in        !! "Rtot": daily average total solar shortwave radiation : direct+diffuse, 280-3000 nm,
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
    REAL(dp),INTENT(out) :: rad_UV           !! "RUV": solar radiation flux (direct+diffuse) from the UV band 280-400 nm
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(dp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(dp),INTENT(out) :: rad_sw,rad_sw_pot

    ! local variables

    REAL(dp) :: day_angle     !! day of year expressed as angle in radians
    REAL(dp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(dp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(dp) :: fsw           !! "fsw": fraction of potential radiation solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                              !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
    REAL(dp) :: F_invers      !! factor 1/F to obtain total shortwave from visible radiation (see Eq. (9))
    REAL(dp) :: hlp_r,hlp2_r,hlp3_r,hlp4_r
    REAL(dp) :: rad_sw_in_corr
    REAL(dp) :: delta,v,u,hh,sinehh,K

    REAL(dp), PARAMETER  ::  frac_UV = 0.1_dp ! UV fraction of (UV + PAR)

    day_angle = 2._dp*pi*(day_of_year+10._dp)/365._dp !! ATTENTION: this may be incorrect for paleo situations

    hlp_r=2._dp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
         1.000110_dp+3.4221E-2_dp*COS(day_angle)+1.280E-3_dp*SIN(day_angle)+7.19E-4_dp*COS(hlp_r)+7.7E-5_dp*SIN(hlp_r)
    K=13750.98708_dp

    IF (forcing_input%shortwave%frequency /= TIMESTEP_) THEN

      ! today's sum of top of the atmosphere radiation
      delta=-23.4_dp*pi/180._dp*COS(day_angle)
      u=sin_lat*SIN(delta)
      v=cos_lat*COS(delta)
      IF (u.GE.v) THEN
         hh=pi
      ELSE IF (u.LE.(0._dp-v)) THEN
         hh=0.0_dp
      ELSE
         hh=ACOS(-u/v)
      ENDIF
      sinehh=SIN(hh)

      !estimate of noon radiation giving rise to the recorded radiation sum
      rad_sw_in_corr = 0._dp
      IF(u*hh+v*sinehh.GT.1e-5_dp)THEN
         rad_sw_in_corr = rad_sw_in / (2._dp*(u*hh+v*sinehh) * K) * 86400._dp
      ENDIF

      rad_sw = rad_sw_in_corr * cos_zenith 

    ELSE
 
      ! time step forcing, use radiation directly
      rad_sw = rad_sw_in
    
    ENDIF
 
    IF(cos_zenith > 0.017452406_dp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
          Rpar_top   = solar_const*cos_zenith*r_sun_earth_2       ! Eq. (2) from above 
          hlp3_r = EXP(-0.185_dp*(air_pressure/p_sealevel)/cos_zenith)
          hlp4_r = 0.4_dp+0.6_dp*hlp3_r
          rad_sw_pot = Rpar_top * hlp4_r 
          fsw = MAX(MIN(rad_sw/rad_sw_pot,1._dp),0._dp)
          hlp_r=1._dp-fsw
          hlp2_r=hlp_r*hlp_r
          F_invers=1._dp + (1.185_dp-0.437_dp*hlp_r-0.494_dp*hlp2_r)* &
                         EXP(MIN(2.534_dp,(0.0305_dp-0.208_dp*hlp_r+0.495_dp*hlp2_r)/cos_zenith)) ! Eq. (9) from above
          rad_UV = rad_sw * frac_UV / F_invers
          rad_PAR = rad_sw * (1._dp - frac_UV) / F_invers   ! This is "Rtot"; Eq. (8) from above
          rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
          IF(fsw <= 0.2_dp) THEN 
             frac_PAR_diffuse = 1.                                                     ! Eq. (12)-(15) from above
          ELSE IF(fsw >= 0.9_dp) THEN
             frac_PAR_diffuse = 0.08_dp                                                ! Eq. (12)-(15) from above
          ELSE
             frac_PAR_diffuse = 1._dp - &
                  (1._dp- ((0.9_dp -fsw)/0.7_dp)**(2._dp/3._dp))*hlp3_r/hlp4_r ! Eq. (12)-(15) from above
          END IF
    ELSE ! zenith angle larger than 89 degrees
          rad_UV = 0._dp 
          rad_PAR = 0._dp
          rad_NIR = 0._dp
          frac_PAR_diffuse = 1._dp ! (could also be another value in this case)
          rad_sw  = 0._dp
          rad_sw_pot  = 0._dp
    END IF

  END SUBROUTINE shortwave_from_direct_shortwave

  ! === shortwave_from_fpar_ccdas() ================================================================================================
  !
  ! This routine implements the shortwave radiation scheme of ccdas. As far as I understand it is almost identical to 
  ! to shortwave_from_fpar_orig(), with the  only exceptions that Equation (9) is replaced by the more simple relationsship  !
  ! 
  ! (9')        F= 0.43 +0.25*(1.-fpar)*cos(theta)
  ! (whereever this comes from).
  !
  PURE ELEMENTAL SUBROUTINE shortwave_from_fpar_ccdas(day_of_year,cos_zenith,air_pressure,fpar, & !! inputs
                                                      rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw   & !! outputs
                                                       )

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(dp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(dp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(dp),INTENT(in)  :: fpar             !! daily mean fraction of actual to potential PAR at bottom (from [0,1])
    REAL(dp),INTENT(out) :: rad_UV           !! "RUV_act": solar radiation flux (direct+diffuse) from the UV band 280-400 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(dp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(dp),INTENT(out) :: rad_sw           !! "Rtot": total solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
 
    ! local variables

    REAL(dp) :: day_angle     !! day of year expressed as angle in radians
    REAL(dp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(dp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(dp) :: F             !! factor F to obtain total shortwave from visible radiation (see Eq. (9'))
    REAL(dp) :: hlp_r, hlp3_r,hlp4_r

    REAL(dp), PARAMETER  ::  frac_UV = 0._dp ! UV fraction of (UV + PAR)
                                             ! to be checked, if the ratio rad_PAR/rad_sw is within observed range (0.42-0.51) with this value of frac_UV

    day_angle = 2._dp*pi*(day_of_year-1._dp)/365._dp !! ATTENTION: this may be incorrect for paleo situations
    hlp_r=2._dp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
         1.000110_dp+3.4221E-2_dp*COS(day_angle)+1.280E-3_dp*SIN(day_angle)+7.19E-4_dp*COS(hlp_r)+7.7E-5_dp*SIN(hlp_r)
    IF(cos_zenith > 0.017452406_dp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
          Rpar_top   = solar_const_vis*cos_zenith*r_sun_earth_2       ! Eq. (2) from above 
          hlp3_r = EXP(-0.185_dp*(air_pressure/p_sealevel)/cos_zenith)
          hlp4_r = 0.4_dp+0.6_dp*hlp3_r
          rad_UV = fpar * Rpar_top * hlp4_r * frac_UV
          rad_PAR = fpar * Rpar_top * hlp4_r * (1._dp - frac_UV)    ! This is Rpar_act; Eq. (7) from above
          hlp_r=1.-fpar
          F = 0.43_dp +0.25_dp*(1.-fpar)*cos_zenith ! Eq. (9') from above
          rad_sw = (rad_PAR + rad_UV)/ F   ! This is "Rtot"; Eq. (8) from above
          rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
          IF(fpar <= 0.2_dp) THEN 
             frac_PAR_diffuse = 1._dp                                                     ! Eq. (12)-(14) from above
          ELSE IF(fpar >= 0.9_dp) THEN
             frac_PAR_diffuse = 0._dp                                                     ! Eq. (11)-(14) from above
          ELSE
             frac_PAR_diffuse = 1._dp - &
                  (1._dp- ((0.9_dp -fpar)/0.7_dp)**(2._dp/3._dp))*hlp3_r/hlp4_r ! Eq. (11)-(14) from above
          END IF
    ELSE ! zenith angle larger than 89 degrees
          rad_UV = 0._dp 
          rad_PAR = 0._dp
          rad_NIR = 0._dp
          frac_PAR_diffuse = 1._dp ! (could also be another value in this case)
          rad_sw  = 0._dp
    END IF
  END SUBROUTINE shortwave_from_fpar_ccdas

  ! === shortwave_from_fsw() =====================================================================================================
  !
  ! This routine is identical to the shortwave routine from CCDAS shortwave_from_fpar_ccdas(), but with a modification proposed by 
  ! Wolfgang to use not fpar as input, but the transmission ratio fsw for TOTAL SHORTWAVE radiation instead.
  !
  ! Here 
  !                     Rtot
  ! (A)        fsw := ----------
  !                     R_pot
  !
  ! the ratio between actual "Rtot" to potential "R_pot" shortwave radiation (daily means). "R_pot" is obtailed from the 
  ! potential PAR radiation by
  !                      Rpar_pot
  ! (B)        R_pot = ---------
  !                       F_pot
  ! where the conversion factor "F_pot" is taken as
  ! (C)        F_pot = 0.43.
  ! With (1) then total incident shortwave radiation follows from
  !                                       Rpar_pot
  ! (D)        Rtot = fsw * R_pot = fsw * ---------
  !                                         F_pot
  ! The PAR part of total actual radiation is then obtained by
  ! (E)        Rpar_act = Rtot * F,
  ! where the conversion factor "F" is
  ! (F)        F= 0.43 +0.25*(1.-fpar)*cos(theta)
  ! Now by equations (10) to (14) from shortwave_from_fpar_orig() all output fields can be computed.
  !
  PURE elemental SUBROUTINE shortwave_from_fsw(day_of_year,cos_zenith,air_pressure,fsw, & !! inputs
                                               rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw  & !! outputs
                                               )

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(dp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(dp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(dp),INTENT(in)  :: fsw              !! daily mean fraction of actual to potential TOTAL SHORTWAVE flux at bottom ..
                                             !! .. (from [0,1])
    REAL(dp),INTENT(out) :: rad_UV           !! "RUV_act": solar radiation flux (direct+diffuse) from the UV band 280-400 nm
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(dp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(dp),INTENT(out) :: rad_sw           !! "Rtot": total solar shortwave radiation : direct+diffuse, 400-3000 nm, ,..
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
 
    ! local variables

    REAL(dp) :: day_angle     !! day of year expressed as angle in radians
    REAL(dp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(dp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(dp) :: F             !! factor F to obtain total shortwave from visible radiation (see Eq. (9'))
    REAL(dp) :: hlp_r,hlp3_r,hlp4_r

    ! parameters

    REAL(dp),PARAMETER :: F_pot = 0.43_dp !! See Eq. (C) from above
    REAL(dp),PARAMETER :: frac_UV = 0._dp !! UV fraction of (UV + PAR)
                                          !! to be checked, if the ratio rad_PAR/rad_sw is within observed range (0.42-0.51) with this value of frac_UV


    day_angle = 2._dp*pi*(day_of_year-1._dp)/365._dp !! ATTENTION: this may be incorrect for paleo situations
    hlp_r=2._dp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
         1.000110_dp+3.4221E-2_dp*COS(day_angle)+1.280E-3_dp*SIN(day_angle)+7.19E-4_dp*COS(hlp_r)+7.7E-5_dp*SIN(hlp_r)
    IF(cos_zenith > 0.017452406_dp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
       Rpar_top   = solar_const_vis*cos_zenith*r_sun_earth_2       ! Eq. (2) from above 
       hlp3_r = EXP(-0.185_dp*(air_pressure/p_sealevel)/cos_zenith)
       hlp4_r = 0.4_dp+0.6_dp*hlp3_r
       rad_sw = fsw*Rpar_top*hlp4_r/F_pot   ! This is "Rtot"; Eq. (D) from above
       hlp_r=1.-fsw
       F = 0.43_dp + 0.25_dp*(1.-fsw)*cos_zenith ! Eq. (9') from shortwave_from_fpar_ccdas(), but with the correct ratio
                                                 ! fsw instead fpar
       rad_UV = rad_sw * F * frac_UV
       rad_PAR = rad_sw * F * (1._dp - frac_UV)  ! This is Rpar_act; Eq. (E) from above
       rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
       IF(fsw <= 0.2_dp) THEN 
          frac_PAR_diffuse = 1._dp                                                     ! Eq. (12)-(14) from above
       ELSE IF(fsw >= 0.9_dp) THEN
          frac_PAR_diffuse = 0._dp                                                  ! Eq. (11)-(14) from above
       ELSE
          frac_PAR_diffuse = 1._dp - &
               (1._dp - ((0.9_dp - fsw)/0.7_dp)**(2._dp/3._dp))*hlp3_r/hlp4_r ! Eq. (11)-(14) from above
       END IF
    ELSE ! zenith angle larger than 89 degrees 
       rad_UV = 0._dp
       rad_PAR = 0._dp
       rad_NIR = 0._dp
       frac_PAR_diffuse = 1._dp ! (could also be another value in this case)
       rad_sw  = 0._dp
    END IF
  END SUBROUTINE shortwave_from_fsw

  ! === shortwave_from_fsw_orig() ==================================================================================================
  !
  ! This routine is identical to the new shortwave routine from Wolfgang shortwave_from_fsw, except that the original conversion 
  ! 
  ! factor
  ! (9)           1/F = 1 + (1.185-0.437*mfpar-0.494*mfpar^2)* exp[(0.0305-0.208*mfpar+0.495*mfpar^2)/cos(theta)]
  ! with mfpar=1-fpar is used (compare Eq. (9) in subroutine shortwave_from_fpar_orig()).
  PURE elemental SUBROUTINE shortwave_from_fsw_orig(day_of_year,cos_zenith,air_pressure,fsw, & !! inputs
                                                    rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw  & !! outputs
                                                    )

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(dp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(dp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(dp),INTENT(in)  :: fsw              !! daily mean fraction of actual to potential TOTAL SHORTWAVE flux at bottom ..
                                             !! .. (from [0,1])
    REAL(dp),INTENT(out) :: rad_UV           !! "RUV_act": solar radiation flux (direct+diffuse) from the UV band 280-400 nm
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(dp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(dp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(dp),INTENT(out) :: rad_sw           !! "Rtot": total solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
 
    ! local variables

    REAL(dp) :: day_angle     !! day of year expressed as angle in radians
    REAL(dp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(dp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(dp) :: F             !! factor F to obtain total shortwave from visible radiation (see Eq. (9'))
    REAL(dp) :: hlp_r,hlp2_r,hlp3_r,hlp4_r

    ! parameters

    REAL(dp),PARAMETER :: F_pot = 0.43_dp !! See Eq. (C) from above
    REAL(dp),PARAMETER :: frac_UV = 0._dp !! UV fraction of (UV + PAR)
                                          !! to be checked, if the ratio rad_PAR/rad_sw is within observed range (0.42-0.51) with this value of frac_UV

    day_angle = 2._dp*pi*(day_of_year-1._dp)/365._dp !! ATTENTION: this may be incorrect for paleo situations
    hlp_r=2._dp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
         1.000110_dp+3.4221E-2_dp*COS(day_angle)+1.280E-3_dp*SIN(day_angle)+7.19E-4_dp*COS(hlp_r)+7.7E-5_dp*SIN(hlp_r)
    IF(cos_zenith > 0.017452406_dp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
       Rpar_top   = solar_const_vis*cos_zenith*r_sun_earth_2       ! Eq. (2) from above 
       hlp3_r = EXP(-0.185_dp*(air_pressure/p_sealevel)/cos_zenith)
       hlp4_r = 0.4_dp+0.6_dp*hlp3_r
       rad_sw = fsw*Rpar_top*hlp4_r/F_pot   ! This is "Rtot"; Eq. (D) from above
       hlp_r=1._dp-fsw
       hlp2_r=hlp_r*hlp_r
       F = 1._dp/(1._dp + (1.185_dp-0.437_dp*hlp_r-0.494_dp*hlp2_r) * &
            EXP(MIN(2.534_dp,(0.0305_dp-0.208_dp*hlp_r+0.495_dp*hlp2_r)/cos_zenith))) ! Eq. (9) from above
       rad_UV = rad_sw * F * frac_UV
       rad_PAR = rad_sw * F * (1._dp - frac_UV)    ! This is Rpar_act; Eq. (E) from above
       rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
       IF(fsw <= 0.2_dp) THEN 
          frac_PAR_diffuse = 1._dp                                                     ! Eq. (12)-(14) from above
       ELSE IF(fsw >= 0.9_dp) THEN
          frac_PAR_diffuse = 0._dp                                                     ! Eq. (11)-(14) from above
       ELSE
          frac_PAR_diffuse = 1._dp - &
               (1._dp - ((0.9_dp -fsw)/0.7_dp)**(2._dp/3._dp))*hlp3_r/hlp4_r ! Eq. (11)-(14) from above
       END IF
    ELSE ! zenith angle larger than 89 degrees
       rad_UV = 0._dp 
       rad_PAR = 0._dp
       rad_NIR = 0._dp
       frac_PAR_diffuse = 1._dp ! (could also be another value in this case)
       rad_sw  = 0._dp
    END IF
  END SUBROUTINE shortwave_from_fsw_orig

  ! === shortwave_from_shortwaveMean() =============================================================================================
  !
  ! This routine computes the various shortwave parts (PAR,NIR,fraction of diffuse PAR) from the mean daily shortwave flux that
  ! arrives at the bottom. The method is the same as that used in the subroutine "shortwave_from_cloudcover()", except that
  ! "fpar" (see description of  "shortwave_from_cloudcover()"), which is the fraction of actual to potential PAR at the earth's
  ! surface, is not computed from cloudcover, but estimated from the ratio of actual to potential mean daily shortwave flux
  !                       R0_act
  ! (1')           fpar = ------ ,
  !                       R0_pot
  ! where the actual mean daily shortwave flux "R0_act" is input to the routine and the potential mean daily shortwave flux
  ! is computed from the average shortwave flux reaching the bottom at a clear day:
  !
  !                                    1                                                             -0.185*p
  ! (A)           R0_pot = S * E *  ----- INTEGRAL[t_up:t_down] dt cos(theta(t)){ 0.4 + 0.6*exp(------------------) },
  !                                   24h                                                         p0*cos(theta(t))
  !
  ! where the term "0.4 + ..." accounts for Rayleigh scattering (which is insufficient!), "E" reflects the varying distance
  ! from the sun, "S=1360 W/m^2" is the solar constant and "t_up" and "t_down" are the times of sunrise and sunset 
  ! (for other symbol names see subroutine "shortwave_from_fpar()").
  ! To compute the integral one needs the explicit representation of the cosine of the zenith-angle theta(t):
  ! (B)          cos(theta(t)) = sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(omega(t))
  ! where "lat" is the latitude, "dec" the declination and "omega(t)" is the hour-angle
  ! (C)          omega(t) = t*2*pi/24h ;
  ! here the time "t" is measured since/to the culmination of the sun at the particular day, i.e. in the morning "t" is positive
  ! and negative in the afternoon. 
  ! The hour angles of sunrise and sunset are obtained from setting "theta(t) = pi/2" so that
  ! (D)          omega_sunrise = -omega_sunset = |arccos(tan(lat)*tan(dec))|.
  ! Since the cos(theta(t)) is symmetric the integral from (A) can be rewritten as
  !                                  1                                                       -0.185*p
  ! (A')          R0_pot = S * E *  -- INTEGRAL[0:omega_sunrise] dx cos(x){ 0.4 + 0.6*exp(--------------) },
  !                                 pi                                                       p0*cos(x)
  ! This value has to be provided when calling the subroutine (input variable "rad_sw_down_pot")
  !
  elemental PURE SUBROUTINE shortwave_from_shortwaveMean(day_of_year,cos_zenith,cos_lat,sin_lat,air_pressure,  & !! inputs
                                                         rad_sw_down_act,rad_sw_down_pot,                    & !! inputs
                                                         rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse, & !! outputs
                                                         rad_sw_down,rad_sw_down_pot_out)
    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(dp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(dp),INTENT(in)  :: cos_lat          !! cosine of latitude 
    REAL(dp),INTENT(in)  :: sin_lat          !! sine of latitude 
    REAL(dp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(dp),INTENT(in)  :: rad_sw_down_act  !! daily mean of total solar shortwave radiation : direct+diffuse, 280-3000 nm, ..
                                             !! .. (i.e. PAR+NIR), that actually arrives at the bottom  [W/m^2]
    REAL(dp),INTENT(in)  :: rad_sw_down_pot  !! daily mean of total solar shortwave radiation : direct+diffuse, 280-3000 nm, ..
                                             !! ..(i.e. PAR+NIR), that potentially arrives at the bottom (i.e. at clear sky) [W/m^2]
    REAL(dp),INTENT(out) :: rad_UV_down      !! solar radiation flux (direct+diffuse) from the UV band 280-400 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: rad_PAR_down     !! solar radiation flux (direct+diffuse) from the visible band 400-700 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: rad_NIR_down     !! solar radiation flux (direct+diffuse) from the near infrared band 700-3000 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(dp),INTENT(out) :: frac_PAR_diffuse !! Fraction from  rad_PAR_down() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(dp),INTENT(out) :: rad_sw_down      !! total solar shortwave radiation flux: direct+diffuse, 280-3000 nm, (i.e. PAR+NIR),..
                                             !! ... that actually arrives at the bottom for the particular zenith angle [W/m^2]
    REAL(dp),INTENT(out) :: rad_sw_down_pot_out !! total solar shortwave radiation flux: direct+diffuse, 280-3000 nm, (i.e. PAR+NIR)
                                             !! ... that potentially arrives at the bottom for the particular zenith angle [W/m^2]

    ! parameters

    REAL(dp),PARAMETER :: negligible_rad = 1e-8_dp !! Below this value (in [W/m^2]) radiation is assumed to be negligible
    ! locals

    REAL(dp) :: fsw         !! Fraction of actual out of potential radiation flux in the full shortwave band at surface

    IF(rad_sw_down_pot < negligible_rad ) THEN
       fsw = 0.5_dp
    ELSE 
       fsw = MIN(rad_sw_down_act/rad_sw_down_pot,1._dp)
    END IF
 
    SELECT CASE(forcing_options%type_of_shortwave_scheme)
       CASE(SW_SCHEME_ORIGINAL_)
          CALL shortwave_from_fpar_orig(day_of_year,cos_zenith,air_pressure,fsw, &   !! inputs
                                        rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse,rad_sw_down)   !! outputs
          rad_sw_down_pot_out = rad_sw_down_pot

       CASE(SW_SCHEME_CCDAS_)
          CALL shortwave_from_fpar_ccdas(day_of_year,cos_zenith,air_pressure,fsw, &   !! inputs
                                        rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse,rad_sw_down)   !! outputs
          rad_sw_down_pot_out = rad_sw_down_pot

       CASE(SW_SCHEME_NEW_)
          CALL shortwave_from_fsw(day_of_year,cos_zenith,air_pressure,fsw, &   !! inputs
                                  rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse,rad_sw_down)   !! outputs
          rad_sw_down_pot_out = rad_sw_down_pot

       CASE(SW_SCHEME_NEW_ORIG_CF_)
          CALL shortwave_from_fsw_orig(day_of_year,cos_zenith,air_pressure,fsw, &   !! inputs
                                       rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse,rad_sw_down)   !! outputs
          rad_sw_down_pot_out = rad_sw_down_pot
 
       CASE(SW_SCHEME_DIRECT_)
          CALL shortwave_from_direct_shortwave(day_of_year,cos_zenith,cos_lat,sin_lat,air_pressure, & !! inputs
                                        rad_sw_down_act, &   !! inputs
                                        rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse, & !! outputs
                                        rad_sw_down,rad_sw_down_pot_out)   !! outputs

   END SELECT

  END SUBROUTINE shortwave_from_shortwaveMean

  ! === vapour_pressure_from_evapor() =============================================================================================
  !
  ! Vapor pressure is computed according to a method described in [1] from actual and potential evapotranspiration. 
  ! Typically the mean vapor pressure at a day is approximately equal to the saturation vapor pressure at sunrise 
  ! "e_sat(T_min)", where temperature is at its minimum. Only in arid regions this approximation gets bad; in these regions 
  ! the saturation vapor pressure is not reached at sunrise. The following model accounts for this observation by using the 
  ! ratio of actual to potential evapotranspiration as a measure of aridity. Let "E_act(d)" and "E_pot(d)" denote the total 
  ! actual and total potential evapotranspiration at day "d" and "fe(d)" the ratio of these values (Eq. (22) in [1]):
  ! (1) fe(d) = E_act(d)/E_pot(d).
  ! For draught situations  E_pot(d)=infty so that fe(d)=0. In the other extreme, fe(d)=1, E_act(d)=E_pot(d), i.e. actual transpi-
  ! ration is only limited by potential evapotranspiration and this means the climate is extremely wet (sufficient soil moisture 
  ! and finite capacity of the atmosphere to absorb water vapor). Let the vapor pressure at sunrise be (Eq. (21) in [1])
  ! (2) e0_vap(d) = (h0 + (1-h0)*fe(d-1))*e_sat(T_min),
  ! where "e_sat(T_min)" is the saturation vapor pressure at sunrise temperature "T_min", "fe(d-1)" the yesterdays
  ! evapotranspiration ratio and "h0" a tuning parameter that has the meaning of the humidity at sunrise for total draught
  ! (fe=0). Here the yesterdays value of the evapotranspiration ratio "fe" is taken. The vapor pressure "e_vap(t)" in the actual
  ! timestep "t" is then estimated from the current temperature "T" by (Eq. (20) in [1])
  ! (3) e_vap(t) = e0_vap(d) + h1*fe(d-1)*(e_sat(T)-e0_vap(d)),
  ! where "h1" is another tuning constant that has the meaning of a daily amplitude of vapor pressure under extreme moist 
  ! conditions (fe=1). According to [1] one finds h0=0.96 and h1=0.49 (see pp. 65/66).
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49, 
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  !
  !vg>>
  ! Above values for h0 and h1 lead to extremely moist conditions, especially in desert areas.
  ! A comparison of qair from echam/jsbach output with qair calculated in this routine using different h0 values shows best
  ! agreement for h0=0.6. The same h0 is found when comparing CRUNCEP qair with qair calculated from CRUNCEP temperatures in
  ! this module. For desert studies it might be good to reduce h0 even further. As the specific humidity does not have a 
  ! significant diurnal cycle, h1 is set to 0. 
  !vg<<

  FUNCTION vapour_pressure_from_evapor(domain, evapo_ratio_mean, Tmin, Tair)

    TYPE(domain_type),INTENT(in) :: domain                     !! information on processor domain
    REAL(dp),         INTENT(in) :: evapo_ratio_mean(1:domain%nland) !! fe(d-1) of Eq. (1)
    REAL(dp),         INTENT(in) :: Tmin(1:domain%nland)       !! minimum temperature (= temperature at sunrise)
    REAL(dp),         INTENT(in) :: Tair(1:domain%nland)       !! temperature at current time step

    REAL(dp), PARAMETER :: h0         = 0.6_dp  !! (0.96_dp) tuning constant for computing vapor pressure
    REAL(dp), PARAMETER :: h1         = 0.0_dp  !! (0.49_dp) tuning constant for computing vapor pressure

    REAL(dp) :: vapour_pressure_from_evapor(1:domain%nland) !! the vapour pressure

    ! locals

    REAL(dp) :: evapo_sat        ! saturation pressure at actual temperature
    REAL(dp) :: evapo_sat_min    ! saturation pressure at minimum temperature of day (= at sunrise)
    REAL(dp) :: e0_vap           ! vapor pressure at sunrise 
    INTEGER :: i

    DO i=1,domain%nland
       evapo_sat_min = tlucua(INT(1000._dp * (Tmin(i)+tmelt)))/eps ! saturation pressure at minimum temperature
       evapo_sat     = tlucua(INT(1000._dp * (Tair(i)+tmelt)))/eps ! saturation pressure at actual temperature
       e0_vap = evapo_sat_min*(h0 +(1._dp-h0)*evapo_ratio_mean(i))                               ! Eq. (2) from above
       vapour_pressure_from_evapor(i) = e0_vap +h1*evapo_ratio_mean(i)*(evapo_sat-evapo_sat_min) ! Eq. (3) from above
    END DO
  END FUNCTION vapour_pressure_from_evapor

  ! === specific_humidity() =======================================================================================================
  !
  ! Specific humidity is computed here by assuming that each component of a mixture of dry air and vapor can be considered as an
  ! ideal gas. So, let m_d and m_v the masses of the two components in a fixed volume V. Then the specific humidity "q" is 
  ! defined as
  ! (1) q = m_v/(m_v+m_d)
  ! Closely related is the moisture mixing ratio "r"
  ! (2) r = m_v/m_d
  ! Hence, 
  ! (3) q = r/(1+r)
  ! Next it is assumed that both, dray air and vapor behave like ideal gases:
  ! (4) p_d*V = m_d*R_d*T and p_v*V = m_v*R_v*T,
  ! where p_d is the partial pressure of the dry air and p_v the vapor pressure, R_d and R_v the gas constants and T the 
  ! temperature. Noting that air pressure "p" is the sum of the two partial pressures (p=p_d+p_v) the division of the two gas 
  ! equation by each other gives by using (2):
  ! (5) r = eps*p_v/(p-p_v),
  ! where eps=R_d/R_v is the ratio of gas constants. Entering ths into (3) one finally obtains for the specific humidity:
  !                 p_v
  ! (6) q = eps -------------
  !             p-(1-eps)*p_v
  !
  PURE elemental FUNCTION specific_humidity(vapour_pressure, air_pressure)
    REAL(dp)            :: specific_humidity !! (from [0,1])
    REAL(dp),INTENT(in) :: vapour_pressure   !! [N/m^2]
    REAL(dp),INTENT(in) :: air_pressure      !! air pressure at bottom [N/m^2]
 
    specific_humidity = eps*vapour_pressure/(air_pressure - (1._dp-eps)*vapour_pressure)

  END FUNCTION specific_humidity

  ! === vapour_pressure() =========================================================================================================
  !
  ! Inverse of specific_humidity()
  !
  PURE elemental FUNCTION vapour_pressure(specific_humidity, air_pressure)
    REAL(dp) :: vapour_pressure   !! [N/m^2]
    REAL(dp),INTENT(in) :: specific_humidity !! (from [0,1])
    REAL(dp),INTENT(in) :: air_pressure      !! air pressure at bottom [N/m^2]

    vapour_pressure = air_pressure * specific_humidity / ( eps + specific_humidity * ( 1._dp - eps ) )

  END FUNCTION vapour_pressure

  ! === longwave_from_daily_longwave() ====================================================================================
  !
  ! Downward long wave radiation flux "R_d" [W/m^2] is according to [1],[2] computed by
  ! (1) R_d = r_cloud * epsA * sigma * T^4,
  ! where "T" is the air temperature [K] (in the bottom layer?), sigma = 5.6703e-8 W/(m^2 K^4) is the Stefan-Boltzmann constant,
  ! The mean daily longwave is known from observations, thus we are concerned only with the diurnal course associated with T
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49, 
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  ! [2] W. Brutsaert, "Evaporation into the Atmosphere", (Reidel, Dordrecht, 1982), pp. 299.
  !
  PURE elemental FUNCTION longwave_from_daily_longwave(longwave_down_daily,Tmin,Tmax,air_temp)
    REAL(dp)            :: longwave_from_daily_longwave !! hourly longwave downward radiation [W/m^2]
    REAL(dp),INTENT(in) :: longwave_down_daily          !! daily average downward longwave radiation [W/m^2] 
    REAL(dp),INTENT(in) :: Tmin                         !! daily mimimum temperature [Celsius] 
    REAL(dp),INTENT(in) :: Tmax                         !! daily maximum temperature [Celsius] 
    REAL(dp),INTENT(in) :: air_temp                     !! air temperature at bottom [Celsius]

    ! locals
    REAL(dp) :: dtmp    ! daily average temperature 
    REAL(dp) :: hlp_r   ! helper

    ! Go ...
    hlp_r=air_temp+Tmelt                                ! Celsius --> Kelvin
    dtmp=0.5_dp*(Tmax+Tmin)+Tmelt
    longwave_from_daily_longwave = longwave_down_daily*(hlp_r**4)/(dtmp**4)      ! lw-radiation flux Eq. (1)

  END FUNCTION longwave_from_daily_longwave 

  ! === longwave_rad_down() =======================================================================================================
  !
  ! Downward long wave radiation flux "R_d" [W/m^2] is according to [1],[2] computed by
  ! (1) R_d = r_cloud * epsA * sigma * T^4,
  ! where "T" is the air temperature [K] (in the bottom layer?), sigma = 5.6703e-8 W/(m^2 K^4) is the Stefan-Boltzmann constant,
  ! "epsA" is the emissivity of the cloudless atmosphere given by
  ! (2) epsA = epsA0 * (e_A/T)^(1/7).
  ! Here "e_vap" is the vapor pressure in [Pa] from above, and epsA0 = 0.64 (units: [(K/Pa)^(1/7)]) a constant (see [2]). 
  ! The factor "r_cloud" in (1) corrects for clouds and is given by
  ! (3) r_cloud = 1+0.22*n^2,
  ! where "n" is the cloudcover.
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49, 
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  ! [2] W. Brutsaert, "Evaporation into the Atmosphere", (Reidel, Dordrecht, 1982), pp. 299.
  !
  PURE elemental FUNCTION longwave_rad_down(cloud_cover,vapor_pressure,air_temp)
    REAL(dp)            :: longwave_rad_down !! longwave downward radiation [W/m^2]
    REAL(dp),INTENT(in) :: cloud_cover       !! fraction of sky covered by clouds (from [0,1])
    REAL(dp),INTENT(in) :: vapor_pressure    !! vapor pressure [N/m^2]
    REAL(dp),INTENT(in) :: air_temp          !! air temperature at bottom [Celsius]

    ! locals

    REAL(dp) :: r_cloud ! factor correcting longwave radiation for clouds
    REAL(dp) :: epsA    ! emissivity of the cloudless atmosphere 
    REAL(dp) :: hlp_r   ! helper

    ! Go ...

    r_cloud=1._dp+0.22_dp*cloud_cover*cloud_cover           ! cloud correction Eq. (3)
    hlp_r=air_temp+Tmelt                                ! Celsius --> Kelvin
    epsA=0.64_dp*(vapor_pressure/hlp_r)**0.14285714_dp   ! emissivity Eq. (2)
    longwave_rad_down = r_cloud*epsA*StefanBoltzmann*hlp_r**4      ! lw-radiation flux Eq. (1)
  END FUNCTION longwave_rad_down

  ! === init_variable_access() =====================================================================================================
  !
  ! This subroutine initializes the access to data of a particular variable: In the case of monthly data these are read in so that 
  ! in update_forcing() daily data can be generated. Or, in the case of daily data, the netcdf-file containing the respective
  ! file is opened, so that in subsequent netcdf-calls daily data can be read in. The whole information on the access to the data
  ! is after return found in a structure of type "input_variable_type". 
  !
  ! In particular the routine does the following:
  !
  !   *** The structure "variable_info" is filled in so that after return of this routine it contains all information
  !       to access the data, either by subsequent netcdf-input (in the case of daily data), or (in case of monthly data)
  !       by recomputing daily data from them.
  !   *** It does some consistency checks with respect to the underlying grid, and finishes the program in case of inconsistency
  !   *** In case of daily data it opens the netcdf-file for further access of the data.
  !   *** It allocates memory for the daily data
  !   *** In case of monthly data it allocates also memory for them and reads them in. The associated netcdf-file is closed
  !       after reading in the data.
  !
  ! The routine scans the JSBACH runtime definition file ("run.def") for keywords to obtain information on input data. Depending 
  ! on the values of the interface variables "keyword_stem" and "variable_name" the routine expects to find certain keywords in 
  ! run.def, i.e. there are certain rules what keywords are expected to be found. These rules are as follows:
  !
  ! Let (for example) be keyword stem = "FORCING_TEMP" and variable_name = "tmin". Then run.def is expected to contain the 
  ! keyword "FORCING_TEMP_FREQU". The value given to this keyword can be either "MONTHLY", "DAILY" or "CONST", that specify
  ! whether input data are monthly, daily or given by a constant value. In the case of  "MONTHLY" and "DAILY" in addition
  ! the file that contains the data has to be specified: here the keyword "FORCING_TEMP_FILE" is used. Instead, in the case
  ! of "CONST" the constant value is specified with "FORCING_TEMP_CONST_TMIN"; here for the construction of the keyword
  ! also the name of the variable ("tmin") is used for construction. In summary, there are three possibilities how the
  ! entries in "run.def" could look like:
  !
  ! Monthly input data:
  !                      FORCING_TEMP_FREQU = MONTHLY
  !                      FORCING_TEMP_FILE  = <fileName (absolute or relative to run.def)>
  !
  ! Daily input data:   
  !                      FORCING_TEMP_FREQU = DAILY
  !                      FORCING_TEMP_FILE  = <fileName (absolute or relative to run.def)>
  !
  ! Constant input data:
  !                      FORCING_TEMP_FREQU = CONST
  !                      FORCING_TEMP_CONST_TMIN = <real value>
  !  
    SUBROUTINE init_variable_access(domain, grid, keyword_stem, variable_name, frequ, forcing_file_base, constant, & 
        input_variable, read_variable)

    USE mo_util_string,   ONLY: tolower
    USE mo_tr_scatter,    ONLY: scatter_field
    USE mo_time_control,  ONLY: delta_time, next_date, get_date_components 
    
    TYPE(domain_type),INTENT(IN)             :: domain        ! Information on processor domain
    TYPE(grid_type) ,INTENT(IN)              :: grid          ! Information on global grid
    CHARACTER(len=*),INTENT(IN)              :: keyword_stem  ! Used to to construct namelist keywords
    CHARACTER(len=*),INTENT(IN)              :: frequ         ! data frequency (DAILY/MONTHLY/CONST)
    REAL(dp),        INTENT(IN)              :: constant      ! constant (only used if frequ=CONST)
    CHARACTER(nf_max_name),INTENT(IN)        :: forcing_file_base  ! name base of file containing the forcing data
    CHARACTER(len=*),INTENT(IN)              :: variable_name ! For monthly and daily input data: name of variable in netcdf file;
                                                              ! For constant value: used to construct keyword giving constant value
    TYPE(input_variable_type),INTENT(INOUT)  :: input_variable ! structure that on return contains all information
    LOGICAL,OPTIONAL,INTENT(IN)  :: read_variable

    INTEGER                  :: nc_file_id            ! NetCDF file ID
    INTEGER                  :: status                ! For return status of netcdf operations
    INTEGER                  :: no_of_lats,no_of_lons ! number of latitudes and longitudes in file opened here
    INTEGER                  :: nc_dim_id             ! ID of dimension variable in file opened here
    INTEGER                  :: nc_var_id             ! NetCDF variable ID
    REAL(dp), ALLOCATABLE,TARGET :: hlp_field3D_global(:,:,:) ! temporary memory for global monthly forcing data
    REAL(dp), ALLOCATABLE        :: hlp_field2D_locally(:,:)  ! temporary memory for a 2D-domain on each processor
    REAL(dp), POINTER,DIMENSION(:,:) :: hlp_real2D_ptr    ! pointer to a real 2D-field

    CHARACTER(len=128)       :: keyword
    CHARACTER(len=8)         :: year_label
    INTEGER                  :: n
    INTEGER                  :: zstart(1), zcount(1)
    INTEGER                  :: time_steps_per_day
    INTEGER                  :: next_year
    LOGICAL                  :: read_input

    time_steps_per_day = 86400/int(delta_time) 

    ! default is to read variable on initialisation
    read_input = .TRUE.
    IF(PRESENT(read_variable)) read_input=read_variable

    ! Go ...

    NULLIFY(input_variable%data_daily)
    NULLIFY(input_variable%data_monthly)
    NULLIFY(input_variable%data_timestep)

    ! Save variable name

    input_variable%variable_name = TRIM(variable_name)

    ! determine frequency of input data from (monthly, daily or timestep) from definition file

    keyword=TRIM(keyword_stem)//'_FREQU'
    SELECT CASE(TRIM(tolower(frequ)))
    CASE("monthly")
       input_variable%frequency = MONTHLY_
    CASE("daily") 
        input_variable%frequency = DAILY_
    CASE("const") 
        input_variable%frequency = CONST_
    CASE("timestep") 
        input_variable%frequency = TIMESTEP_
    CASE("-")
       CALL finish("init_variable_access()",'Keyword "'//TRIM(keyword)//'" missing in namelist.')
    CASE default
       CALL finish("init_variable_access()", &
            'Invalid value '//TRIM(frequ)//' of keyword "'//TRIM(keyword)//'" in namelist.')
    END SELECT
    CALL message("init_variable_access()", &
         'For forcing variable "'//variable_name//'" frequency of data is: '//TRIM(frequ))

    ! allocate memory for daily or timestep data on all processors (only memory for a single day is needed, because at the beginning
    ! of each day these data are either (i) newly generated from monthly data, (ii) input from file, 
    ! (iii) initialised for the time step of this day or (iv) set to a constant value)

    IF (input_variable%frequency == TIMESTEP_) THEN
         ALLOCATE(input_variable%data_timestep(1:domain%nland,1:time_steps_per_day),STAT=status)
         IF(status /= 0) CALL finish("init_variable_access()", &
               'ERROR: For variable  "'//TRIM(variable_name)//'" memory for timestep data could not be allocated.')
         input_variable%data_timestep(:,:)=0._dp
    ELSE
         ALLOCATE(input_variable%data_daily(1:domain%nland),STAT=status)
         IF(status /= 0) CALL finish("init_variable_access()", &
            'ERROR: For variable  "'//TRIM(variable_name)//'" memory for daily data could not be allocated.')
         input_variable%data_daily(:)=0._dp
    ENDIF

    IF (p_parallel_io) THEN

       ! Prepare for reading in data (only if data have to be input from a netcdf file)
       IF(input_variable%frequency .NE. CONST_ .AND. read_input ) THEN

          ! Name of file containing the forcing data (namelist forcing_ctl)

          CALL get_date_components(next_date, YEAR=next_year)
          WRITE(year_label,'(i8.4)') next_year
          WRITE(input_variable%file_name,'(a,a,a)') TRIM(forcing_file_base)//'_', TRIM(ADJUSTL(year_label)), '.nc'
          CALL message("init_variable_access()", &
               'For forcing variable "'//TRIM(variable_name)//'" using data from file: '//TRIM(input_variable%file_name))

          ! open forcing file and determine netcdf file-id

          status = nf_open(input_variable%file_name,nf_nowrite,nc_file_id)
          IF (status /= NF_NOERR) THEN
             CALL finish("init_variable_access()"," ERROR: Could not open forcing data file: "//TRIM(input_variable%file_name))
          END IF
          input_variable%file_ID = nc_file_id

          ! check whether dimensions are correct

          !determine dimension id for "lat" or "latitude"
          status = nf_inq_dimid(nc_file_id,"lat",nc_dim_id)
          IF (status /= NF_NOERR) THEN ! try dimension name "latitude
             status = nf_inq_dimid(nc_file_id,"latitude",nc_dim_id)
             IF (status /= NF_NOERR) CALL finish("init_variable_access()", &
               'Latitude dimension name ("lat" or "latitude") is missing in forcing file '//TRIM(input_variable%file_name))
          END IF
          ! determine and check the number of latitudes
          status = nf_inq_dimlen(nc_file_id,nc_dim_id,no_of_lats)
          IF (status /= NF_NOERR) CALL finish("init_variable_access()", &
               "Dimension of lat could not be read from forcing file "//TRIM(input_variable%file_name))
          IF(no_of_lats /= grid%nlat) THEN
             CALL message("init_variable_access()","Inconsistency between grid file and forcing file " &
                  //TRIM(input_variable%file_name))
             CALL finish("init_variable_access()", &
                  "Number of latitudes inconsistent: grid file: "//TRIM(int2string(grid%nlat))//&
                  " forcing file: "//int2string(no_of_lats))
          END IF
    
          !determine dimension id for "lon" or "longitude"
          status = nf_inq_dimid(nc_file_id,"lon",nc_dim_id)
          IF (status /= NF_NOERR) THEN
             status = nf_inq_dimid(nc_file_id,"longitude",nc_dim_id)
             IF (status /= NF_NOERR) CALL finish("init_variable_access()", &
                            'Longitude dimension name ("lon" or "longitude") is missing in forcing file ' &
                            //TRIM(input_variable%file_name))
          END IF
          ! determine and check the number of longitudes
          status = nf_inq_dimlen(nc_file_id,nc_dim_id,no_of_lons)
          IF (status /= NF_NOERR) &
               CALL finish("init_variable_access()","Dimension of lon could not be read from forcing file " &
               //TRIM(input_variable%file_name))
          IF(no_of_lons /= grid%nlon) THEN
             CALL message("init_variable_access()","Inconsistency between grid file and forcing file "// &
                  TRIM(input_variable%file_name)) 
             CALL finish("init_variable_access()", &
                  "Number of longtitudes inconsistent: grid file: "//TRIM(int2string(grid%nlon))//&
                  " forcing file: "//int2string(no_of_lons))
          END IF

          ! determine variable ID

          status = nf_inq_varid (nc_file_id,TRIM(variable_name), nc_var_id)
          IF (status /= NF_NOERR) CALL finish("init_variable_access()",  &
               'ERROR: No variable "'//TRIM(variable_name)//'" found in forcing data file '//TRIM(input_variable%file_name))
          input_variable%var_ID = nc_var_id

          ! determine the length of the time series
          
          status = nf_inq_dimid(nc_file_id,"time",nc_dim_id) ! determine dimension id for "time"
          IF (status /= NF_NOERR) &
               CALL finish("init_variable_access()","time dimension is missing in forcing file "//TRIM(input_variable%file_name))
          status = nf_inq_dimlen(nc_file_id,nc_dim_id, &
                                          input_variable%length_time_series)! determine the number of time steps
          IF (status /= NF_NOERR) CALL finish("init_variable_access()", &
                                         'Dimension of "time" could not be read from forcing file '//TRIM(input_variable%file_name))
          CALL message('init_variable_access()', 'Found '//TRIM(int2string(input_variable%length_time_series))//&
             ' records for variable '//TRIM(variable_name))

          ! get times of time series
          status = nf_inq_varid(nc_file_id, 'time', nc_var_id)
          ALLOCATE(input_variable%timevals(input_variable%length_time_series))
          zstart(1) = 1
          zcount(1) = input_variable%length_time_series
          status = nf_get_vara_double(nc_file_id, nc_var_id, zstart(1),zcount(1), input_variable%timevals)
          IF (status /= NF_NOERR) &
               CALL finish("init_variable_access()","time variable is missing in forcing file "//TRIM(input_variable%file_name))

          ! check length of timeseries

          SELECT CASE(input_variable%frequency)
          CASE(MONTHLY_)
             IF((input_variable%length_time_series .LE. 0) .OR. (MODULO(input_variable%length_time_series,12) .NE. 0)) THEN
                CALL finish("init_variable_access()",&
                     "Number of months ("//TRIM(int2String(input_variable%length_time_series))//") for variable "// &
                     TRIM(variable_name)//" not multiple of 12")  
             END IF
          CASE(TIMESTEP_)
             IF((input_variable%length_time_series .LE. 0) .OR. &
                (MODULO(input_variable%length_time_series,time_steps_per_day) .NE. 0)) THEN
                CALL finish("init_variable_access()",&
                     "Number of time steps ("//TRIM(int2String(input_variable%length_time_series))//") for variable "// &
                     TRIM(variable_name)//" not multiple of "//TRIM(int2String(time_steps_per_day))//")")  
             END IF
          END SELECT

       END IF

    END IF

    ! Broadcasting information (partially) to other processors

    IF (p_parallel) THEN
       CALL p_bcast(input_variable%file_name, p_io)
       CALL p_bcast(input_variable%frequency, p_io)
       CALL p_bcast(input_variable%length_time_series, p_io)
    END IF
    
    ! Additional preparations depending on frequency:

    SELECT CASE(input_variable%frequency)

    CASE(CONST_)!! A single value is used for the whole globe

       !! Check the namelist value of the constant

       IF (p_parallel_io) THEN
          keyword = TRIM(keyword_stem)//'_CONST_'//TRIM(variable_name)
          IF(constant == HUGE(constant)) &
               CALL finish("init_variable_access()", &
               'For variable '//TRIM(variable_name)//' keyword "'//TRIM(keyword)//'" is missing in namelist forcing_ctl')
       END IF
       !! set on all processors the daily field to constant value 
       input_variable%data_daily(:) = constant 

       CALL message("init_variable_access()",&
            'Variable '//TRIM(variable_name)//' is set globally to '//real2string(constant))

    CASE(MONTHLY_)!! additional preparations for monthly input data

       ! Check the number of months found in the data:

       IF(input_variable%length_time_series /= 12 .AND. read_input ) THEN
          CALL message("init_variable_access()", &
               'Variable '//TRIM(variable_name)//': data are taken from '//TRIM(input_variable%file_name)) 
          CALL finish("init_variable_access()", &
               'Variable '//TRIM(variable_name)//': '//TRIM(int2string(input_variable%length_time_series))//&
               ' months found in forcing file!')
       END IF

       ! allocate memory for monthly fields on all processors

       ALLOCATE(input_variable%data_monthly(1:domain%nland,1:12),STAT=status)
       IF(status /= 0) CALL finish("init_variable_access()", &
            'ERROR: For variable  "'//TRIM(variable_name)//'" memory for monthly data could not be allocated.')
       input_variable%data_monthly(:,:)=0._dp

       ! Read in monthly data
       IF( read_input ) THEN
       IF (p_parallel_io ) THEN 

          ! allocate temporary memory for reading in data

          ALLOCATE(hlp_field3D_global(1:grid%nlon,1:grid%nlat,1:12),STAT=status)
          IF(status /= 0) CALL finish("init_variable_access()", "ERROR: Field hlp_field3D_global() could not be allocated.")

          ! Read monthly data into temporary memory hlp_field3D_global() on io-processor
          status = nf_get_var_double(input_variable%file_ID,input_variable%var_ID, hlp_field3D_global)
          IF (status /= NF_NOERR) CALL finish('init_variable_access()', &
               'For variable "'//TRIM(variable_name)//'" reading data from forcing file ' &
               //TRIM(input_variable%file_name)//' failed.')

       END IF

       ! allocate temporary memory for 2D-fields

       ALLOCATE(hlp_field2D_locally(domain%ndim,domain%nblocks),STAT=status) ! memory for a 2D-domain on each processor
       IF(status /= 0) CALL finish("init_variable_access()()", "ERROR: Field hlp_field2D_locally() could not be allocated.")
       
       ! scatter data from io-processor to all processors and pack them, i.e. keep data only at landpoints

       DO n=1,12 ! go through all months
          NULLIFY(hlp_real2D_ptr)
          IF (p_parallel_io) hlp_real2D_ptr => hlp_field3D_global(:,:,n) ! Let hlp_real2D_ptr point to the global field ...
                                                                         ! ... of a particular month 
          CALL scatter_field(hlp_real2D_ptr,hlp_field2D_locally)         ! scatter global field to all processors
          NULLIFY(hlp_real2D_ptr)
          input_variable%data_monthly(:,n) = &
               PACK(hlp_field2D_locally,MASK=domain%mask) ! pack data to keep them only at landpoints
       END DO

       ! Deallocation of temporary memory

       IF (p_parallel_io) DEALLOCATE(hlp_field3D_global)
       DEALLOCATE(hlp_field2D_locally)

       ENDIF

       ! Close access to data file, because for monthly data no further data have to be read in

    END SELECT
  
    input_variable%initialized = .TRUE.

  END SUBROUTINE init_variable_access

  ! === finish_variable_access() ==================================================================================================

  ! This routine deallocates memory initially allocated for a variable and closes associated input files.
  SUBROUTINE finish_variable_access(input_variable)

    TYPE(input_variable_type),INTENT(INOUT)  :: input_variable ! the variable to be considered

    INTEGER :: status

    if (debug) CALL message('finish_variable_access','Finish variable '//TRIM(input_variable%variable_name))

    IF(.NOT. input_variable%initialized) THEN
       CALL message('finish_variable_access',&
            'WARNING: attemp to finish uninitialized variable '//TRIM(input_variable%variable_name))
       RETURN
    END IF

    IF (p_parallel_io) THEN
       status =  nf_close(input_variable%file_ID)
    END IF


    IF (ASSOCIATED(input_variable%data_timestep)) THEN 
       DEALLOCATE(input_variable%data_timestep)
       NULLIFY(input_variable%data_timestep)
       IF(input_variable%frequency .NE. TIMESTEP_) THEN
         CALL finish('finish_variable_access',&
           'PROGRAMMING ERROR: attempted to deallocate unallocated memory for timestep variable '// &
           TRIM(input_variable%variable_name))
       END IF
    END IF

    IF(ASSOCIATED(input_variable%data_daily)) THEN 
       IF(input_variable%frequency == TIMESTEP_) THEN
         CALL finish('finish_variable_access',&
           'PROGRAMMING ERROR: attempted to deallocate unallocated memory for daily variable '//TRIM(input_variable%variable_name))
       END IF
       DEALLOCATE(input_variable%data_daily)
       NULLIFY(input_variable%data_daily)
    END IF

    IF(ASSOCIATED(input_variable%data_monthly)) THEN 
       DEALLOCATE(input_variable%data_monthly)
       NULLIFY(input_variable%data_monthly)
       IF(input_variable%frequency .NE. MONTHLY_) THEN
          CALL finish('finish_variable_access',&
          'ERROR: For variable '//TRIM(input_variable%variable_name)//' there is erroneously memory  allocated for monthly data')
       END IF
    END IF
    
    input_variable%initialized = .FALSE.
    
  END SUBROUTINE finish_variable_access

  ! === read_in_daily_data() ======================================================================================================
  !
  ! This routine reads in daily data into the field "data_daily" of a structure of type "input_variable_type".
  !
  SUBROUTINE read_in_daily_data(grid,domain,input_variable)

    USE mo_tr_scatter,       ONLY: scatter_field
    USE mo_time_control,     ONLY: next_date, get_date_components

    TYPE(grid_type),          INTENT(in)    :: grid            !! information on global grid
    TYPE(domain_type),        INTENT(in)    :: domain          !! information on processor domain
    TYPE(input_variable_type),INTENT(inout) :: input_variable  !! variable to be read in

    INTEGER                     :: start_index(3) !! For the first index in the netcdf file where next days data begin
    INTEGER                     :: count_index(3) !! Length of data block along each dimension to be read in from netcdf-file 
    INTEGER                     :: tsID           !! Index of current day in netcdf file
    INTEGER                     :: status
    REAL(dp),POINTER,DIMENSION(:,:) :: hlp_real2D_ptr !! pointer to a real 2D-field
    INTEGER :: yr, mo, dy, hr, mn, se, ydate

    IF(input_variable%frequency /= DAILY_) RETURN !! Return for non-daily data input

    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)
    ydate = ISIGN(1,yr)*(IABS(yr)*10000+mo*100+dy)

    if (debug) CALL message('read_in_daily_data()',&
               'Read in variable '//TRIM(input_variable%variable_name)//' for day '// &
               TRIM(int2string(ydate))//' from '//TRIM(input_variable%file_name))
    
    !! read in global field

    IF (p_parallel_io) THEN 

       ! Check day number
       tsID = 1
       DO WHILE (tsID .LE. input_variable%length_time_series)
          IF ( ydate == INT(input_variable%timevals(tsID)) ) EXIT
          tsID = tsID + 1
       END DO
       IF (debug) CALL message('read_in_daily_data()', &
            'Index of day '//TRIM(int2string(ydate))//' in netCDF file: '//TRIM(int2string(tsID)))
       IF(tsID > input_variable%length_time_series) THEN 
          CALL message('read_in_daily_data()',&
               'Variable '//TRIM(input_variable%variable_name)//': To be read in from '//TRIM(input_variable%file_name))
          CALL finish('read_in_daily_data()', &
               'Day '//TRIM(int2string(ydate))//' not found')
       END IF

       ! Read next days data into temporary memory hlp_field2D_global() on io-processor

       start_index(1) = 1
       start_index(2) = 1
       start_index(3) = tsID
       count_index(1) = grid%nlon     
       count_index(2) = grid%nlat
       count_index(3) = 1

       status = nf_get_vara_double(input_variable%file_ID, input_variable%var_ID, &
                                              start_index(:),count_index(:), hlp_field2D_global(:,:))
       IF (status /= NF_NOERR) THEN
          CALL message('read_in_daily_data()',&
               'Variable '//TRIM(input_variable%variable_name)//': To be read in from '//TRIM(input_variable%file_name))
          CALL message('read_in_daily_data()',&
               'Reading in day '//TRIM(int2string(ydate))//' (record no. '//TRIM(int2string(tsID))//') failed')
       END IF

    END IF
    
    ! scatter data from io-processor to all processors and pack them, i.e. keep data only at landpoints

    NULLIFY(hlp_real2D_ptr)
    IF (p_parallel_io) hlp_real2D_ptr => hlp_field2D_global(:,:) ! Let hlp_real2D_ptr point to the global field
    CALL scatter_field(hlp_real2D_ptr,hlp_field2D_local)         ! scatter global field to all processors
    NULLIFY(hlp_real2D_ptr)
    input_variable%data_daily(:) = PACK(hlp_field2D_local,MASK=domain%mask) ! pack data to keep them only at landpoints

  END SUBROUTINE read_in_daily_data
  

  ! === read_in_timestep_data() ====================================================================================================
  !
  ! This routine reads in timestep data into the field "data_timestep" of a structure of type "input_variable_type".
  ! NB: This routine only works if model is run at subdaily timesteps...
  !
  SUBROUTINE read_in_timestep_data(grid,domain,input_variable)

    USE mo_tr_scatter,    ONLY: scatter_field
    USE mo_time_control,  ONLY: next_date, get_date_components, delta_time  

    TYPE(grid_type),          INTENT(in)    :: grid            !! information on global grid
    TYPE(domain_type),        INTENT(in)    :: domain          !! information on processor domain
    TYPE(input_variable_type),INTENT(inout) :: input_variable  !! variable to be read in

    INTEGER                         :: start_index(3) !! For the first index in the netcdf file where next days data begin
    INTEGER                         :: count_index(3) !! Length of data block along each dimension to be read in from netcdf-file 
    INTEGER                         :: tsID           !! Index of current day in netcdf file
    INTEGER                         :: tsIDstart, tsIDend   !! start and end indexes...
    INTEGER                         :: status
    REAL(dp),POINTER,DIMENSION(:,:) :: hlp_real2D_ptr !! pointer to a real 2D-field
    INTEGER                         :: yr, mo, dy, hr, mn, se, ydate
    INTEGER                         :: i, time_steps_per_day, ith_timestep_from_seconds 

    IF   (input_variable%frequency /= TIMESTEP_) RETURN 

    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)
    ydate = yr*10000+mo*100+dy

    time_steps_per_day = 86400/INT(delta_time) 

    IF(debug) THEN
       WRITE(message_text,*) 'Read in variable ', TRIM(input_variable%variable_name), ' from ', TRIM(input_variable%file_name)
       CALL message('read_in_timestep_data', message_text)
       WRITE(message_text,'(a,I8.4,5(a,I2.2))') ' date: ', yr, '-', mo, '-', dy, ' time: ', hr, ':', mn, ':', se
       CALL message('read_in_timestep_data', message_text)
    END IF 
    
    
    !! read in global field
    tsIDstart = 0 
    tsIDend = 0  
    IF (p_parallel_io) THEN 

       ! Check subdaily timestep number
       tsID = 1
       DO WHILE (tsID .LE. input_variable%length_time_series)
          IF ( ydate == INT(input_variable%timevals(tsID)) ) THEN
            IF (tsIDstart == 0) THEN
               IF(debug) THEN
                  CALL message('read_in_timestep_data()',&
                    TRIM(input_variable%variable_name)//' timevals @ tsIDstart: '// & 
                    TRIM(int2string(INT(input_variable%timevals(tsID))))// &
                    ' for '//TRIM(input_variable%variable_name))
               END IF
               tsIDstart = tsID
            END IF
            tsIDend = tsID
          ELSE IF ( ydate < INT(input_variable%timevals(tsID)) ) THEN
            EXIT 
          END IF
          tsID = tsID + 1
       END DO

       IF(debug) CALL message('read_in_timestep_data()',&
                    'tsID: '//TRIM(int2string(tsID))//' tsIDstart: '//TRIM(int2string(tsIDstart))// & 
                    'tsIDend: '//TRIM(int2string(tsIDend))// &
                    ' for '//TRIM(input_variable%variable_name))

       IF(tsIDstart > input_variable%length_time_series) THEN 
          CALL message('read_in_timestep_data()',&
               'Variable '//TRIM(input_variable%variable_name)//': To be read in from '//TRIM(input_variable%file_name))
          CALL finish('read_in_timestep_data()', &
               'Day '//TRIM(int2string(ydate))//' not found')
       END IF

       ! Read next days data into temporary memory hlp_field3D_global() on io-processor 
       
       ! specify how many records to extract from the file... 
       IF (time_steps_per_day .NE. (tsIDend - tsIDstart + 1)) THEN
          CALL finish('read_in_timestep_data()', &
            '# of records for the day '//TRIM(int2string(ydate))//' ('//TRIM(int2string(tsIDend-tsIDstart+1))// &
            ') is not correct ('//TRIM(int2string(time_steps_per_day))//') for the current temporal resolution (' &
            //TRIM(int2string(int(delta_time)))//' seconds).')
       END IF 

       start_index(1) = 1
       start_index(2) = 1
       start_index(3) = tsIDstart
       count_index(1) = grid%nlon     
       count_index(2) = grid%nlat
       count_index(3) = time_steps_per_day

       IF(debug) CALL message('read_in_timestep_data()','counter_index(3) is: '//TRIM(int2string(count_index(3)))) 

       status = nf_get_vara_double(input_variable%file_ID, input_variable%var_ID, &
                                              start_index(:),count_index(:), hlp_field3D_global(:,:,:)) 
       
       IF (status /= NF_NOERR) THEN
          CALL message('read_in_timestep_data()',&
               'Variable '//TRIM(input_variable%variable_name)//': To be read in from '//TRIM(input_variable%file_name))
          CALL message('read_in_timestep_data()',&
               'Reading in day '//TRIM(int2string(ydate))//' (record no. '//TRIM(int2string(tsIDstart))//') failed')
       END IF

    END IF
    
    ! scatter data from io-processor to all processors and pack them, i.e. keep data only at landpoints
    NULLIFY(hlp_real2D_ptr)

    DO i=1,time_steps_per_day 
      IF (p_parallel_io) THEN
        ith_timestep_from_seconds = NINT(MOD(input_variable%timevals(tsIDstart+i-1),1._dp)*86400._dp)/INT(delta_time)+1
        IF (ith_timestep_from_seconds /= i) THEN
          !!-nc: if it not matching time step... we could have a way to go around it... but instead just finish!
          !!-nc: only full days of data and ordered records are possible, which is cool, because then we do not 
          !!-nc: need to worry about the time step and time variables in the update forcing at timestep resolutions
          CALL finish('read_in_timestep_data()', &
            'the records for the day '//TRIM(int2string(ydate))// &
            ' are not ordered... ith_timestep_from_seconds is '//TRIM(int2string(ith_timestep_from_seconds))// &
            'while "i" is ('//TRIM(int2string(i))//').')
        END IF
      END IF
      IF (p_parallel_io) hlp_real2D_ptr => hlp_field3D_global(:,:,i) ! Let hlp_real2D_ptr point to the global field
      CALL scatter_field(hlp_real2D_ptr,hlp_field2D_local)         ! scatter global field to all processors
      NULLIFY(hlp_real2D_ptr)
      input_variable%data_timestep(:,i) = PACK(hlp_field2D_local,MASK=domain%mask) ! pack data to keep them only at landpoints
    END DO

  END SUBROUTINE read_in_timestep_data

  ! === init_diagnostic_output() ===================================================================================================
  !
  ! This routine initializes the forcing stream for diagnostic output. Implicitely memory is allocated for all fields of the
  ! structure "forcing_diag".
  !
  SUBROUTINE init_diagnostic_output(grid, domain, fileformat)

    USE mo_jsbach,        ONLY : missing_value
    USE mo_jsbach_grid,   ONLY: grid_type, domain_type
    USE mo_linked_list,   ONLY: LAND, t_stream
    USE mo_memory_base,   ONLY: new_stream, default_stream_setting, &
                                add =>add_stream_element

    TYPE(grid_type),  INTENT(in) :: grid
    TYPE(domain_type),INTENT(in) :: domain
    INTEGER,          INTENT(in) :: fileformat

    ! Local variables
    TYPE(t_stream),POINTER,SAVE  :: forcingStream
    INTEGER :: dim1_glob(1)
    INTEGER :: dim1_loc(1)

    IF (debug) CALL message('init_diagnostic_output','defining stream for forcing variables')

    ! Add new stream
    CALL new_stream(forcingStream, 'forcing', filetype=fileformat)
    ! Set default stream options
    CALL default_stream_setting(forcingStream, repr=LAND, lpost=.TRUE., lrerun=.TRUE.)

    dim1_glob(1) = grid%nland
    dim1_loc(1)  = domain%nland

    CALL add(forcingStream,'air_temp',        forcing_diag%forcingFields%air_temp,         dim1_loc, dim1_glob, &
         longname = 'Air Temperature (2m)',                   units = 'deg C',                                  &
         code = 1, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'precip_rain',     forcing_diag%forcingFields%precip_rain,      dim1_loc, dim1_glob, &
         longname = 'liquid water precipiation above canopy', units = 'kg m-2 s-1',                             & 
         code = 2, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'precip_snow',     forcing_diag%forcingFields%precip_snow,      dim1_loc, dim1_glob, &
         longname = 'frozen water precipitaion above canopy', units = 'kg m-2 s-1',                             &
         code = 3, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'air_pressure',    forcing_diag%forcingFields%air_pressure,     dim1_loc, dim1_glob, &
         longname = 'Air Pressure (2m)',                      units = 'Pa',                                     &
         code = 4, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'vapor_pressure',  forcing_diag%forcingFields%vapor_pressure,   dim1_loc, dim1_glob, &
         longname = 'Air Water Vapour Pressure (2m)',         units = 'Pa',                                     &
         code = 5, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'spec_humidity',   forcing_diag%forcingFields%spec_humidity,    dim1_loc, dim1_glob, &
         longname = 'Specific Air Humidity (2m)',             units = 'kg kg-1',                                &
         code = 6, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'rad_sw_down',     forcing_diag%forcingFields%rad_sw_down,      dim1_loc, dim1_glob, &
         longname = 'incident downward shortwave radiation',  units = 'W m-2',                                  & 
         code = 7, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'rad_PAR_down',    forcing_diag%forcingFields%rad_PAR_down,     dim1_loc, dim1_glob, &
         longname = 'photosynthetically active downward radiation', units = 'W m-2',                            & 
         code = 8, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'rad_NIR_down',    forcing_diag%forcingFields%rad_NIR_down,     dim1_loc, dim1_glob, &
         longname = 'near infrared downward radiation',       units = 'W m-2',                                  &
         code = 9, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'frac_PAR_diffuse',forcing_diag%forcingFields%frac_PAR_diffuse, dim1_loc, dim1_glob, &
         longname = 'diffuse fraction of photosynthetically active radiation', units = '--',                    & 
         code = 10, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'rad_lw_down',     forcing_diag%forcingFields%rad_lw_down,      dim1_loc, dim1_glob, &
         longname = 'incident longwave downward radiation',   units = 'W m-2',                                  &
         code = 11,  laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'wind_speed',      forcing_diag%forcingFields%wind_speed,       dim1_loc, dim1_glob, &
         longname = 'Wind Speed (10m)',                       units = 'm s-1',                                  & 
         code = 12, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'CO2_concentr',    forcing_diag%forcingFields%CO2_concentr,     dim1_loc, dim1_glob, &
         longname = 'Air CO2 concentration', units = 'kg kg-1',                                                 & 
         code = 13, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(forcingStream,'rad_sw_down_pot',  forcing_diag%forcingFields%rad_sw_down_pot, dim1_loc, dim1_glob, &
         longname = 'incident potential downward shortwave radiation',  units = 'W m-2',                        & 
         code = 14, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)

!!$    CALL add(forcingStream,'rad_UV_down',     forcing_diag%forcingFields%rad_UV_down,      dim1_loc, dim1_glob, &
!!$         longname = 'UV downward radiation', units = 'W m-2',                                                   & 
!!$         code = 15, laccu = .TRUE., lmiss=.TRUE., missval=missing_value)

    forcing_diag%forcingFields%air_temp(:)        = 0._dp 
    forcing_diag%forcingFields%precip_rain(:)     = 0._dp
    forcing_diag%forcingFields%precip_snow(:)     = 0._dp
    forcing_diag%forcingFields%air_pressure(:)    = 0._dp
    forcing_diag%forcingFields%vapor_pressure(:)  = 0._dp
    forcing_diag%forcingFields%spec_humidity(:)   = 0._dp
    forcing_diag%forcingFields%rad_sw_down(:)     = 0._dp
!!$    forcing_diag%forcingFields%rad_UV_down(:)    = 0._dp
    forcing_diag%forcingFields%rad_PAR_down(:)    = 0._dp
    forcing_diag%forcingFields%rad_NIR_down(:)    = 0._dp
    forcing_diag%forcingFields%frac_PAR_diffuse(:) = 0._dp
    forcing_diag%forcingFields%rad_lw_down(:)     = 0._dp
    forcing_diag%forcingFields%wind_speed(:)      = 0._dp
    forcing_diag%forcingFields%CO2_concentr(:)    = 0._dp
    forcing_diag%forcingFields%rad_sw_down_pot(:) = 0._dp

  END SUBROUTINE init_diagnostic_output

  ! === update_diagnostic_output() =================================================================================================
  !
  ! This routine updates the diagnostic output fields
  !
  SUBROUTINE update_diagnostic_output(forcingFields)

    USE mo_time_control, ONLY: dt => delta_time

    TYPE(forcing_type), INTENT(in)          :: forcingFields     ! The current forcing fields

    forcing_diag%forcingFields%air_temp(:)        = &
                                                forcing_diag%forcingFields%air_temp(:)        + dt*forcingFields%air_temp(:)
    forcing_diag%forcingFields%precip_rain(:)     = &
                                                forcing_diag%forcingFields%precip_rain(:)     + dt*forcingFields%precip_rain(:)
    forcing_diag%forcingFields%precip_snow(:)     = &
                                                forcing_diag%forcingFields%precip_snow(:)     + dt*forcingFields%precip_snow(:)
    forcing_diag%forcingFields%air_pressure(:)    = &
                                                forcing_diag%forcingFields%air_pressure(:)    + dt*forcingFields%air_pressure(:)
    forcing_diag%forcingFields%vapor_pressure(:)  = &
                                                forcing_diag%forcingFields%vapor_pressure(:)  + dt*forcingFields%vapor_pressure(:)
    forcing_diag%forcingFields%spec_humidity(:)   = &
                                                forcing_diag%forcingFields%spec_humidity(:)   + dt*forcingFields%spec_humidity(:)
    forcing_diag%forcingFields%rad_sw_down(:)     = &
                                                forcing_diag%forcingFields%rad_sw_down(:)     + dt*forcingFields%rad_sw_down(:)
!!$    forcing_diag%forcingFields%rad_UV_down(:)    = &
!!$                                                forcing_diag%forcingFields%rad_UV_down(:)     + dt*forcingFields%rad_UV_down(:)
    forcing_diag%forcingFields%rad_PAR_down(:)    = &
                                                forcing_diag%forcingFields%rad_PAR_down(:)    + dt*forcingFields%rad_PAR_down(:)
    forcing_diag%forcingFields%rad_NIR_down(:)    = &
                                                forcing_diag%forcingFields%rad_NIR_down(:)    + dt*forcingFields%rad_NIR_down(:)
    forcing_diag%forcingFields%frac_PAR_diffuse(:) = &
                                                forcing_diag%forcingFields%frac_PAR_diffuse(:)+ dt*forcingFields%frac_PAR_diffuse(:)
    forcing_diag%forcingFields%rad_lw_down(:)     = &
                                                forcing_diag%forcingFields%rad_lw_down(:)     + dt*forcingFields%rad_lw_down(:)
    forcing_diag%forcingFields%wind_speed(:)      = &
                                                forcing_diag%forcingFields%wind_speed(:)      + dt*forcingFields%wind_speed(:)
    forcing_diag%forcingFields%CO2_concentr(:)    = &
                                                forcing_diag%forcingFields%CO2_concentr(:)    + dt*forcingFields%CO2_concentr(:)
    forcing_diag%forcingFields%rad_sw_down_pot(:)     = &
                                                forcing_diag%forcingFields%rad_sw_down_pot(:) + dt*forcingFields%rad_sw_down_pot(:)
 
  END SUBROUTINE update_diagnostic_output
END MODULE mo_jsbalone_forcing
