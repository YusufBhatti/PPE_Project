!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_phenology

  ! 
  ! This module contains the new JSBACH-phenology LOGRO-P, based on logistic growth. It's main purpose is to compute the leaf 
  ! area index (LAI) at all land-gridpoints for every phenology type (this is done by the subroutine "update_phenology"). 
  !
  ! For a complete description see the LaTeX-file "phenology.tex" under jsbach/doc.
  !
  ! Author: C.H. Reick, MPI-BGC, 15.6.2003
  !         C.H. Reick, MPI-BGC, 17.11.2003, including NPP-criterion for leaf sheddin-rates of raingreens and grasses
  ! Remarks:
  !    -- In this module all times are measured in fractions of days!!
  !    -- Currently only natural vegetation is considered
  !    -- so far no restart mechanism
  !    -- A daylength criterion for assuring the begin of the rest phase for summer- and evergreen has not been implemented yet.
  !    -- How PFTs are related to landcover types is coded in the landcover library file, whose data are provided by
  !       the structure lctlib in module mo_jsbach_lctlib.
  !
  ! Abbreviations:
  !    EG ... evergreen
  !    SG ... summergreen
  !    RG ... raingreen
  !    GRS .. grasses

  USE mo_jsbach_grid, ONLY: kstart, kend, nidx
  USE mo_kind,        ONLY: dp 
  USE mo_jsbach,      ONLY: debug, new_day
  USE mo_exception,   ONLY: finish, message
  USE mo_util_string, ONLY: int2string, real2string
  USE mo_mpi,         ONLY: p_parallel_io

  IMPLICIT NONE  

  REAL(dp),PARAMETER :: NaN = -HUGE(1.0_dp)           ! Take minus 0.1 times the largest possible real as "not-a-number"
  REAL(dp),PARAMETER :: LARGE = 1.e+10_dp             ! A large vale

  ! === BEGIN OF PUBLIC PART =======================================================================================================

  ! --- public subroutines ---

  PUBLIC :: init_phenology             ! Allocates and initializes memory, sets standard parameters and reads in phenology-type data
  PUBLIC :: update_phenology           ! Recomputes the LAIs for all phenology types at all grid points.
  PUBLIC :: read_phenology_parameters  ! Reads phenology parameters from file
  PUBLIC :: print_phenology_parameters ! Prints out all parameters of the phenology model
  PUBLIC :: set_phenology_parameters   ! By this routine phenology parameters can be prescribed externally from a parameter set 
                                       ! given as structure
  PUBLIC :: set_pheno_params_serial    ! By this routine phenology parameters can be prescribed externally from a parameter set
                                       ! given as array
  PUBLIC :: get_phenology_parameters   ! This routine returns the currently used parameters of the phenology model
  PUBLIC :: get_pheno_params_serial    ! Returns the currently used phenology parameters as an array of rdata
  PUBLIC :: copy_phenology_parameters  ! Copies a set of phenology parameters
  PUBLIC :: init_pheno_params          ! Initializes phenology parameters; especially: sets names of parameters and default values
  PUBLIC :: init_status_pheno_params   ! Returns whether the phenology parameters used are initialized (.true.) or not (.false.)
  PUBLIC :: set_rdata                  ! Sets the values of a structure of type rdata.
  PUBLIC :: check_rdata_bounds         ! Checks whether a value of type rdata lies within bounds
  PUBLIC :: print_rdata_bounds         ! Writes the boundary values of number of type rdata to a string

  ! --- phenology type coding --- NEVER NEVER CHANGE THIS !!!! ----------

  INTEGER, PARAMETER, PUBLIC :: no_of_phenotypes = 5 !! Number of phenology types (not including "unvegetated")
  
  INTEGER, PARAMETER, PUBLIC :: unvegetated = 0   !! no vegetation present
  INTEGER, PARAMETER, PUBLIC :: evergreen   = 1
  INTEGER, PARAMETER, PUBLIC :: summergreen = 2
  INTEGER, PARAMETER, PUBLIC :: raingreen   = 3
  INTEGER, PARAMETER, PUBLIC :: grasses     = 4 
  INTEGER, PARAMETER, PUBLIC :: crops       = 5

  ! --- number of parameters ------------------

  INTEGER, PARAMETER, PUBLIC :: noOfPhenParams = 39 !! this is the total number of independent parameters of LoGro-P. This 
                                                    !! .. number derives from counting the parameters kept in 
                                                    !! .. PhenologyParameters_type. WARNING: A wrong values leads to crashes!

  ! --- public types used for parameters -----------------------------------------------------------------------------------------

  TYPE bounds_type  !! For saving parameter boundary values
     LOGICAL  ::  exists = .FALSE.   !! indicates whether a boundary value exists
     REAL(dp)  ::  value  =  NaN     !! the boundary value if it exists
  END TYPE bounds_type
  PUBLIC :: bounds_type

  TYPE rdata_type  !! To have a real value of a variable, its name and possible upper and lower range together. Should be set by
                   !! routine set_rdata().
     REAL(dp)          :: v =  NaN   !! the value of the variable
     CHARACTER(len=32) :: n =  ""    !! the name of the variable
     TYPE(bounds_type) :: uBound     !! the upper bound
     TYPE(bounds_type) :: lBound     !! the lower bound
     CHARACTER(len=32) :: units =""  !! the physical units of the stored values
  END TYPE rdata_type
  PUBLIC :: rdata_type

  TYPE GeneralParams_type                            !! Collects parameters that apply to all phenology types 
     TYPE(rdata_type),POINTER :: LAI_negligible      !! Below this value an LAI is considered to be zero.
     TYPE(rdata_type),POINTER :: laiSeed             !! Seed value used for the LAI to have a non-zero starting value 
                                                     !! .. whenever the vegetation starts growing.
     TYPE(rdata_type),POINTER :: wilt_point          !! Critical fraction of soil water bucket level: below this value growth of 
                                                     !! .. of grasses, raingreens and crops stops (wilting point).
  END TYPE GeneralParams_type
  PUBLIC :: GeneralParams_type

  TYPE EG_SG_common_params_type                      !! Parameters that are common to the evergreen and raingreen phenology type
     TYPE(rdata_type),POINTER :: tau_pseudo_soil     !! Characteristic time for the memory loss for computing the pseudo soil 
                                                     !! .. temperature from air temperature [days] (see subroutine 
                                                     !! .. "update_pseudo_soil_temp")
     TYPE(rdata_type),POINTER :: max_chill_days      !! This is an upper limit for the number of chill days (field chill_days(:)). 
                                                     !! .. It prevents the number of chill days to grow beyond any limit in regions
                                                     !! .. where the temperature is permanently below the alternation_temp (i.e.
                                                     !! .. especially in polar regions)
  END TYPE EG_SG_common_params_type
  PUBLIC :: EG_SG_common_params_type

  TYPE evergreen_params_type                         !! Parameters applying to the evergreen phenology type only
     TYPE(rdata_type), POINTER :: alternation_temp   !! Critical temperature of the alternating model, above which temperature
                                                     !! .. contributes to the heat sum, and below which days are considered as 
                                                     !! .. chilling days (see subroutine "update_growth_phase") [Celsius]
     TYPE(rdata_type), POINTER :: heat_sum_min       !! Minimum value of critical heat sum [degree days] (see "update_growth_phase")
     TYPE(rdata_type), POINTER :: heat_sum_range     !! Range of critical heat sum [degree days](see "update_growth_phase") 
     TYPE(rdata_type), POINTER :: chill_decay_const  !! Number of chill days at which chilling influence on critical heat sum
                                                     !! .. drops to 1/e (see "update_growth_phase")
     TYPE(rdata_type), POINTER :: growthPhaseLength  !! Length of growth phase [days]
     TYPE(rdata_type), POINTER :: shedRate_rest      !! Leaf shedding rate of evergreen during the vegetative phase [1/days]
     TYPE(rdata_type), POINTER :: growthRate         !! Fraction of NPP maximally allocated to leaves during growth     
  END TYPE evergreen_params_type
  PUBLIC :: evergreen_params_type

  TYPE summergreen_params_type                       !! Parameters applying to the summergreen phenology type only
     TYPE(rdata_type), POINTER :: alternation_temp   !! Critical temperature of the alternating model, above which temperature
                                                     !! .. contributes to the heat sum, and below which days are considered as 
                                                     !! .. chilling days (see subroutine "update_growth_phase") [Celsius]
     TYPE(rdata_type), POINTER :: heat_sum_min       !! Minimum value of critical heat sum [degree days] (see "update_growth_phase")
     TYPE(rdata_type), POINTER :: heat_sum_range     !! Range of critical heat sum [degree days](see "update_growth_phase") 
     TYPE(rdata_type), POINTER :: chill_decay_const  !! Number of chill days at which chilling influence on critical heat sum
                                                     !! .. drops to 1/e (see "update_growth_phase")
     TYPE(rdata_type), POINTER :: growthPhaseLength  !! Length of growth phase [years]
     TYPE(rdata_type), POINTER :: autumn_event_temp  !! Critical pseudo-soil-temperature that determines the autumn event, i.e. the
                                                     !! .. date at which rest phase begins [Celsius] (see "update_growth_phase")
                                                     !! .. [Celsius]
     TYPE(rdata_type), POINTER :: shedRate_veget     !! Leaf shedding rate of summergreen (SG) during the vegetative phase [1/days]
     TYPE(rdata_type), POINTER :: shedRate_rest      !! Leaf shedding rate of summergreen (SG) during the rest phase [1/days]
     TYPE(rdata_type), POINTER :: growthRate         !! Fraction of NPP maximally allocated to leaves during growth
     TYPE(rdata_type), POINTER :: maxLength_gvPhase  !! Number of days that growth plus vegetative phase maximally can have [days]
                                                     !! .. This parameter is introduced for technical reasons to assure that the
                                                     !! next growth pahse is not missed.
     TYPE(rdata_type), POINTER :: minLength_gvPhase  !! Number of days that growth plus vegetative phase minimally should last
                                                     !!   [days]. This parameter is introduced to prevent leaf shedding when early
                                                     !!   after the growth phase there is a cold snap 

  END TYPE summergreen_params_type
  PUBLIC :: summergreen_params_type

  TYPE raingreen_params_type !! Parameters applying to the raingreen phenology type only
     TYPE(rdata_type), POINTER :: shedRate_drySeason !! Leaf shedding rate (fast) in dry season [1/days]
     TYPE(rdata_type), POINTER :: shedRate_aging     !! Minimal leaf shedding rate by leaf aging (inverse leaf longevity) [1/days]
     TYPE(rdata_type), POINTER :: growthRate         !! Growth rate (only active during wet season; modified by leaf shedding)
                                                     !! [1/days]
     TYPE(rdata_type), POINTER :: bucketFill_critical !! If bucket filling drops below this value, the shedding rate is increased
                                                      !!   so that plant growth is reduced until it gets zero at the wilting point
                                                      !!   [fraction]
     TYPE(rdata_type), POINTER :: bucketFill_leafout !! The critical bucket filling at which plant start growing leaves 
                                                     !! .. (>= wiltingt point) [fraction]
  END TYPE raingreen_params_type
  PUBLIC :: raingreen_params_type

  TYPE grass_params_type     !! Parameters applying to the grass phenology type only
     TYPE(rdata_type), POINTER :: crit_temp          !! Critical temperature for growth of grasses [Celsius]: Below this 
                                                     !! .. temperature growth of grasses stops.
     TYPE(rdata_type), POINTER :: shedRate_growth    !! Leaf shedding rate in the growth phase when growth conditions are not so 
                                                     !!   good [1/days]
     TYPE(rdata_type), POINTER :: growthRate         !! Growth rate for good climate conditions [1/days]  
     TYPE(rdata_type), POINTER :: shedRate_drySeason !! Leaf shedding rate in dry season [1/days]
  END TYPE grass_params_type
  PUBLIC :: grass_params_type

  TYPE crop_param_type        !! Parameters applying to the crop phenology only
     TYPE(rdata_type), POINTER :: crit_temp          !! Critical temperature for growth of crops [Celsius]: Below this 
                                                     !! .. temperature growth of crops stops.
     TYPE(rdata_type), POINTER :: gdd_temp           !! Critical (base) temperature for counting growing degree days (i.e. heat sum)
                                                     !!   of crops
     TYPE(rdata_type), POINTER :: sproud             !! Critical fraction of soil water bucket level for setting LAI to the seed
                                                     !!   value
     TYPE(rdata_type), POINTER :: heat_sum_harvest   !! Heat sum (degree days) at which crops are harvested
     TYPE(rdata_type), POINTER :: shedRate_growth    !! Leaf shedding rate in the growth phase [1/days]
     TYPE(rdata_type), POINTER :: shedRate_rest      !! Shed rate for cold season of crops [1/days]
     TYPE(rdata_type), POINTER :: leafAlloc_frac     !! Fraction of NPP maximally allocated to leaves during growth     
  END TYPE crop_param_type
  PUBLIC :: crop_param_type

  TYPE PhenologyParameters_type                      !! Collects all parameters in a single structure
     TYPE(GeneralParams_type)       :: all           !! Parameters common to all phenology types
     TYPE(EG_SG_common_params_type) :: EG_SG         !! Parameters common to the evergreen (EG) and summergreen (SG) phenology types
     TYPE(evergreen_params_type)    :: EG            !! Parameters applying to the evergreen (EG) phenology type only
     TYPE(summergreen_params_type)  :: SG            !! Parameters applying to the summergreen (SG) phenology type only 
     TYPE(raingreen_params_type)    :: RG            !! Parameters applying to the raingreen (RG) phenology type only
     TYPE(grass_params_type)        :: GRS           !! Parameters applying to the grass (GRS) phenology type only 
     TYPE(crop_param_type)          :: CRP           !! Parameters applying to the crop (CRP) phenology type only
     TYPE(rdata_type), POINTER      :: memory(:)  => NULL()   !! pointer to serial array containing the parameter values physically
     logical                        :: initialized =.FALSE.   !! .true. indicates that phenology parameters have been initialized.
  END TYPE PhenologyParameters_type
  PUBLIC :: PhenologyParameters_type

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE ! Make ALL following objects private

  ! === Private declarations =======================================================================================================

  TYPE(PhenologyParameters_type),SAVE,TARGET :: PhenParams           !! Structure containing pointers to all parameters of the
                                                                     !!   phenology model

  INTEGER, SAVE :: ntiles  !! number of tiles
  INTEGER, SAVE :: nsoil   !! number of soil layers to be used


  ! --- time information (recomputed with each call of "update_phenology"), needed to detect midnight between two time steps

  REAL(dp), SAVE          :: timeStep_in_days          !! the time step expressed in days (set in "initPhenology")


  ! --- private fields (state variables and fields for intra-module communication)

  REAL(dp), SAVE, POINTER :: growth_phase_SG(:)   !! For summergreen only: values are:
                                                  !!   -1.0  .. during rest phase (i.e. from autumn event to spring event)
                                                  !!    0.0  .. during growth phase (i.e. some weeks after spring event)
                                                  !!    1.0  .. in vegetative phase (i.e. from end of growth phase to autum event)
                                                  !! Updated every day by the subroutine "update_growth_phase()".

  REAL(dp), SAVE, POINTER :: growth_phase_EG(:)   !! For evergreen only: values are:
                                                  !!    0.0  .. during growth phase (i.e. some weeks after spring event)
                                                  !!    1.0  .. in vegetative phase (i.e. from ent of growth phase to spring event)
                                                  !! Updated every day by the subroutine "update_growth_phase()".

  REAL(dp), SAVE, POINTER :: growth_phase_CRP(:)  !! For extra-tropical crops only: values are:
                                                  !!     -2.0  .. during rest phase of winter crops (i.e. from spring to autumn)
                                                  !!     -1.0  .. during rest phase of summer crops (i.e. from autumn to spring)
                                                  !!      0.0  .. during growth phase of summer crops (i.e. from spring to autumn)
                                                  !!      1.0  .. during growth phase of second crop (i.e. from summer to autumn)
                                                  !!      2.0  .. during growth phase of winter crops (i.e. from autumn to spring)
                                                  !! Updated every day by the subroutine "update_growth_phase()".

  REAL(dp), SAVE, POINTER :: days_since_growth_begin_EG(:) !! For summergreen and evergreen only: serves to remember the number of
                                                           !! days elapsed since the growth phase began 
  REAL(dp), SAVE, POINTER :: days_since_growth_begin_SG(:) !! (as real numbers to save it in a stream). 

  REAL(dp), SAVE, POINTER :: previous_day_temp(:) !! mean day air-temperature [Celsius] of previous day, because not known for ..
                                                  !! current day; updated at the beginning of each new day (stream)

  REAL(dp), SAVE, POINTER :: previous_day_temp_min(:) !! minimum air-temperature [Celsius] of previous day, because not known for ..
                                                      !! current day; updated at the beginning of each new day (stream)

  REAL(dp), SAVE, POINTER :: previous_day_temp_max(:) !! maximum air-temperature [Celsius] of previous day, because not known for ..
                                                      !! current day; updated at the beginning of each new day (stream)

  REAL(dp), SAVE, POINTER :: pseudo_soil_temp(:)  !! Pseudo-soil-temperature [Celsius]: a weighted running mean of the land ...
                                                  !! surface temperature, to simulate soil temperature. Updated every time step ..
                                                  !! by the subroutine "pseudo_soil_temp()". (stream)

  REAL(dp), SAVE, POINTER :: previous_day_NPP_rate(:,:)!! The mean primary production rate [mol(CO2)/(m^2 s)] at each grid point
                                                  !! for each PFT at previous day, because not known for current day; updated at ..
                                                  !! the beginning of each new day (stream) 


  REAL(dp), SAVE, POINTER :: laiMax_dyn(:,:)      !! Dynamical maximum LAI [m^2(leaf)/m^2(ground)]-- not larger than the maximum lai
                                                  !! coming from JSBACH, but smaller than that value when NPP is negative for more
                                                  !! than 1 year. (stream) 
  REAL(dp), SAVE, POINTER :: year_NPP_sum(:,:)    !! This field sums the NPP-rate during the whole year, starting on Jan,1, in the
                                                  !! northern hemisphere and at July,1, in the southern hemisphere. A negative value
                                                  !! indicates that summation has not yet started.

  ! --- internal fields of subroutine "update_previous_day_variables"

  REAL(dp), SAVE, POINTER :: day_temp_sum(:)      !! sum of air temperatures since midnight (needed to compute mean air 
                                                  !! temperature) [Celsius * 1] (stream)

  REAL(dp), SAVE, POINTER :: day_temp_min(:)      !! minimum of air temperatures since midnight (needed to compute growing 
                                                  !! degree days GDD) [Celsius * 1] (stream)

  REAL(dp), SAVE, POINTER :: day_temp_max(:)      !! maximum of air temperatures since midnight (needed to compute growing 
                                                  !! degree days GDD) [Celsius * 1] (stream)

  REAL(dp), SAVE, POINTER :: day_NPP_sum(:,:)     !! sum of NPP since midnight (needed to compute previous_day_NPP_rate(:,:))
                                                  !! [mol(CO2)/m^2]

  ! --- internal fields of subroutine "update_growth_phase"

  REAL(dp), SAVE, POINTER :: heat_sum_EG(:),    & !! Heat sums --- needed for the phenology of summergreens, evergreens and crops.
                             heat_sum_SG(:),    & !! Where heat summation has not yet started, "heat_sum" is set to a value < -1.
                             heat_sum_CRP(:),   & !! For summer crops heat summation starts with "0".
                             heat_sum_winter(:)   !! For winter crops heat summation starts with EPSILON(1.) and is positive in the
                                                  !! growing period (autumn to spring) and negative during summer. (stream)

  REAL(dp), SAVE, POINTER :: chill_days_EG(:),  & !! Number of chill days --- needed for the phenology of summer- and evergreens 
                             chill_days_SG(:)     !! (stream)

  ! --- constants needed in "update_pseudo_soil_temp" (initialized in "init_phenology")

  REAL(dp), SAVE          :: N_pseudo_soil_temp   !! Normalization for computing the pseudo soil temperature
  REAL(dp), SAVE          :: F_pseudo_soil_temp   !! exponential factor used for updating the pseudo soil temperature 
  INTEGER, SAVE           :: time_steps_per_day   !! Number of time steps per day

CONTAINS

  ! --- init_phenology() -----------------------------------------------------------------------------------------------------------
  !
  ! This routine initializes the phenology module. It has to be called before the first time step.
  !
  SUBROUTINE init_phenology(g_nland, l_nland, zntiles, n_soil_layers,isRestart, fileformat, fileztype, stream, parameters_serial)
    USE mo_time_control,ONLY: delta_time !! time step in seconds
    USE mo_memory_base, ONLY: new_stream,default_stream_setting, &
                              add =>add_stream_element
    USE mo_linked_list, ONLY: t_stream, LAND, TILES
    USE mo_netCDF,      ONLY: max_dim_name
    USE mo_output,      ONLY: veg_table
    USE mo_jsbach,      ONLY: missing_value

    INTEGER, INTENT(in)                   :: g_nland       ! number of global landpoints
    INTEGER, INTENT(in)                   :: l_nland       ! number of landpoints in domain
    INTEGER, INTENT(in)                   :: zntiles       ! The number of plant functional types to be used by the phenology module
    INTEGER, INTENT(in)                   :: n_soil_layers ! The number of soil layers to be used by the phenology module
    LOGICAL, INTENT(in)                   :: isRestart
    INTEGER, INTENT(in)                   :: fileformat    ! output file format
    INTEGER, INTENT(in)                   :: fileztype     ! output file compression
    TYPE(rdata_type),intent(in),optional  :: parameters_serial(:) ! The parameters (as serial array %memory()) with which the
                                                                  ! .. phenology shall be run. If this argument is missing standard
                                                                  ! .. parameters are used. 
    TYPE(t_stream), POINTER, OPTIONAL     :: stream

    ! local variables

    TYPE(t_stream),POINTER      :: phenStream

    INTEGER                     :: dim2p(2), dim2(2)
    CHARACTER(LEN=max_dim_name) :: dim2n(2)

    ! Go ...

    IF (debug) CALL message('init_phenology','Start initialization of PHENOLOGY')

    ! --- save number of soil layers and number of plant functional types to be used by the phenology module

    nsoil = n_soil_layers
    if(nsoil .le. 0) CALL finish('init_phenology','Initialization with non positive number of soil layers')
    ntiles = zntiles

    ! --- compute the number of time steps per day

    IF( mod(86400,INT(delta_time)) /= 0) THEN ! For computing the mean day temperature (see below) it is assumed that each day ..
                                              ! is computed with a fixed number of time steps. Therefore the program is  ..
                                              ! stopped wheen day length is not an integer multiple of the time step.
       CALL finish("init_phenology","ERROR: Day length is not an integer multiple of the time step!")
    ELSE 
       time_steps_per_day = 86400/int(delta_time)
    END IF

    ! --- express the time step in days

    timeStep_in_days = delta_time / 86400._dp

    ! --- set the phenology parameters

    IF(PRESENT(parameters_serial)) THEN
       CALL set_pheno_params_serial(parameters_serial) !! Use external parameters
    ELSE
       CALL set_phenology_parameters()                 !! Use standard parameters
    END IF

    ! --- define phenology stream

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'pheno', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       phenStream => stream
    ELSE
       ! Add new stream
       CALL new_stream(phenStream, 'pheno', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(phenStream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    ! --- Define state variables as stream elements

    dim2p = (/ l_nland, zntiles /)
    dim2  = (/ g_nland, zntiles /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'tiles'

    ! code numbers of this routine are 181-197, 199 and 205-208 using the GRIB veg_table
    CALL add(phenStream,'growth_phase_SG',            growth_phase_SG,   longname='Growth Phase of Summergreens', &
             units='', code=181, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'growth_phase_EG',            growth_phase_EG,   longname='Growth Phase of Evergreens', &
             units='', code=182, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'growth_phase_CRP',           growth_phase_CRP,  longname='Growth Phase of extra-tropical crops', &
             units='', code=196, lmiss=.TRUE., missval=missing_value)    
    CALL add(phenStream,'days_since_growth_begin_EG', days_since_growth_begin_EG, &
             longname='Number of Days since Beginning of Growth Phase - Evergreens',   units='', code=183, lpost=.FALSE.)
    CALL add(phenStream,'days_since_growth_begin_SG', days_since_growth_begin_SG, &
             longname='Number of Days since Beginning of Growth Phase - Summergreens', units='', code=184, lpost=.FALSE.)
    CALL add(phenStream,'previous_day_temp',          previous_day_temp, longname='Air Temperature of the Previous Day', &
             units='degC', code=185, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'previous_day_temp_min',  previous_day_temp_min, longname='Minimum air Temp. of the Previous Day', &
             units='degC', code=205, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'previous_day_temp_max',  previous_day_temp_max, longname='Maximum air Temp. of the Previous Day', &
             units='degC', code=206, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'pseudo_soil_temp',           pseudo_soil_temp,  longname='Pseudo Soil Temperature', &
             units='degC', code=186, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'day_temp_sum',               day_temp_sum,      longname='Sum of Air Temperature since Midnight', &
             units='degC', code=187, lpost=.FALSE.)
    CALL add(phenStream,'day_temp_min',               day_temp_min,      longname='Minimum of Air Temp. since Midnight', &
             units='degC', code=207, lpost=.FALSE.)
    CALL add(phenStream,'day_temp_max',               day_temp_max,      longname='Maximum of Air Temp. since Midnight', &
             units='degC', code=208, lpost=.FALSE.)
    CALL add(phenStream,'heat_sum_EG',                heat_sum_EG,       longname='Heat Sum - Evergreens', &
             units='degC d', code=188, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'heat_sum_SG',                heat_sum_SG,       longname='Heat Sum - Summergreens', &
             units='degC d', code=189, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'heat_sum_CRP',               heat_sum_CRP,      longname='Heat Sum - Crops', &
             units='degC d', code=197, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'heat_sum_winter',            heat_sum_winter,   longname='Heat Sum - Winter Crops', &
             units='degC d', code=199, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'chill_days_EG',              chill_days_EG,     longname='Number of Chill Days - Evergreens', &
             units='',       code=190, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'chill_days_SG',              chill_days_SG,     longname='Number of Chill Days - Summergreens', &
             units='',       code=191, lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'previous_day_NPP_rate',  previous_day_NPP_rate, longname='NPP Rate of the Previous Day', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=192, lpost=.FALSE., &
             lmiss=.TRUE., missval=missing_value)
    CALL add(phenStream,'day_NPP_sum',                day_NPP_sum,       longname='Sum of NPP since Midnight', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=193, lpost=.FALSE.)
    CALL add(phenStream,'laiMax_dyn',                 laiMax_dyn,        longname='Maximum Leaf Area Index', &
             units='',                         ldims=dim2p, gdims=dim2, dimnames=dim2n, code=194, lpost=.FALSE.)
    CALL add(phenStream,'year_NPP_sum',               year_NPP_sum,      longname='Sum of NPP since beginning of the Year', &
             units='mol(CO2) m-2(canopy) s-1', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=195, lpost=.FALSE.) 

    ! --- for computation of the pseudo soil temperature:

    F_pseudo_soil_temp = EXP(-timeStep_in_days/PhenParams%EG_SG%tau_pseudo_soil%v)
    N_pseudo_soil_temp = 1._dp / (1._dp - F_pseudo_soil_temp)

    laiMax_dyn(:,:) = -1._dp

    ! If this is a restart run we can exit now since the model state variables are read from restart file
    IF (isRestart) THEN
       RETURN
    ENDIF

    ! --- (ad hoc) initialization of private fields (lateron these fields should be initialized from external data!)

    growth_phase_SG(:) = -1._dp             !! We are during the SG-rest phase
    growth_phase_EG(:) = -1._dp             !! We are during the EG_vegetative phase (= not in growth phase)
    growth_phase_CRP(:) = -1._dp            !! We are during the CRP-rest phase
    days_since_growth_begin_EG(:) = 0._dp   !! We are during the EG-rest phase -> no growth days so far
    days_since_growth_begin_SG(:) = 0._dp   !! We are during the SG-rest phase -> no growth days so far
    pseudo_soil_temp(:) = 10._dp            !! An extremely warm winter temperature
    previous_day_temp(:) = 0._dp            !! We are during winter
    previous_day_temp_min(:) = 0._dp        !! We are during winter
    previous_day_temp_max(:) = 0._dp        !! We are during winter
    day_temp_sum(:) = 0._dp                 !! Computation of mean temperature has not started yet
    day_temp_min(:) = 0._dp                 !! Computation of minimum temperature has not started yet
    day_temp_max(:) = 0._dp                 !! Computation of maximum temperature has not started yet
    heat_sum_EG(:) = -99._dp                !! Heat summation has not yet started
    heat_sum_SG(:) = -99._dp                !! Heat summation has not yet started
    heat_sum_CRP(:) = 9999._dp              !! Heat sum is beyond harvest
    heat_sum_winter(:) = -EPSILON(1._dp)    !! Heat summation has not yet started
    chill_days_EG(:) = 0._dp                !! Counting of chill days has not yet started
    chill_days_SG(:) = 0._dp                !! Counting of chill days has not yet started
    previous_day_NPP_rate(:,:) = 0._dp      !! Computing of previous days NPP-rate has not yet started
    day_NPP_sum(:,:) = 0._dp                !! Computing of previous days NPP-rate has not yet started
    year_NPP_sum(:,:) = 2._dp * LARGE       !! Yearly summation has not yet started.

    ! --- debug message 

    IF (p_parallel_io)  CALL print_phenology_parameters
    CALL message('init_phenology()','Initialization of PHENOLOGY finished.')


  END SUBROUTINE init_phenology


  ! --- update_phenology() ---------------------------------------------------------------------------------------------------------
  ! 
  ! This is the main routine of the phenology module. 
  !
  ! Remark: This routine has to be called every time step. In case the whole grid is not processed in a single call, it may be 
  ! called several times during one time step --- in that case the interface variables do not cover all landpoints associated 
  ! with the particular processor, but only only those of a single block of the processor domain. Therefore, the routine has to 
  ! be called with the appropriate section of the processor fields. Example: Let  procArray(1:domain%nland) be a field covering all 
  ! landpoints a particular processor handles (e.g. the LAI). The block currently processed by the processor is 
  ! procArray(kstart:kend) --- this is what has to be handed over to update_phenology().
  !
  SUBROUTINE update_phenology(kidx,pheno_type_of_tile,lai,laiMax_stat,air_temp,relative_moisture,NPP_rate,specLeafArea,lat)

    INTEGER, INTENT(IN)           :: kidx                      !! the number of elements in the current call (= kend - kstart +1)
    INTEGER, INTENT(in)           :: pheno_type_of_tile(:,:)   !! The phenology types of the tiles of grid points (including
                                                               !! "unvegetated"); dimension: (1:nidx,1:ntiles)
    REAL(dp),INTENT(inout),TARGET :: lai(:,:)                  !! Leaf area index to be recomputed. Dimension: (1:nidx,1:ntiles)
    REAL(dp),INTENT(in)           :: laiMax_stat(:,:)          !! Prescribed (=static) maximum LAI for the current year.
                                                               !! Dimension: (1:nidx,1:ntiles)
    REAL(dp),INTENT(in)           :: air_temp(:)               !! Air temperature at current time step in lowest atmospheric layer
                                                               !! for each grid point. Dimension: (1:nidx)
    REAL(dp),INTENT(in)           :: relative_moisture(:,:)    !! Filling of the soil water bucket. Dimension: (1:nidx,1:ntiles)
    REAL(dp),INTENT(in)           :: NPP_rate(:,:)             !! Net primary production rate density [mol(CO2)/m^2/s] for each PFT.
                                                               !! Needed to determine the shedding rate of raingreens and grasses.
                                                               !! Dimension: (1:nidx,1:ntiles) 
    REAL(dp),INTENT(in)           :: specLeafArea(:,:)         !! Specific leaf area in molar units [m^2(leaf)/mol(Carbon)]
                                                               !! Dimension: (1:nidx,1:ntiles) 
    REAL(dp), INTENT(in)          :: lat(:)                    !! latitude of grid points(1:nidx)

    ! --- local variables ------------------------

    REAL(dp)          :: soil_moist_filling(kidx,nsoil,ntiles) !! Fractional filling of the soil water buckets in the different soil
                                                               !! layers. It is assumed that the soil layers are counted from the
                                                               !! surface down (index 1: topmost layer).
    INTEGER           :: i,j,k,kidx0,kidx1                     !! index counters
    REAL(dp)          :: grRate                                !! growth rate
    REAL(dp)          :: average_filling(kidx)                 !! fractional average filling of soil water buckets    
    LOGICAL           :: positive_NPP(kidx)                    !! flag to indicate the existence of a grass or crop type with a
                                                               !! positive NPP
    REAL(dp)          :: xtmp                                  !! dummy argument to calls of letItGrow and letItDie
                                                               !! (necessary to allow for inlining of these functions)
    ! --- set block range indices
    kidx0 = kstart    ! first value of kindex() (the index of the first element of a block in the processor domain)
    kidx1 = kend      ! last value of kindex()  (the index of the last element of a block in the processor domain)

    ! --- set initially the dynamical maximum LAI to the static maximum LAI
    IF (ANY(laiMax_dyn(kidx0:kidx1,:) < 0._dp)) laiMax_dyn(kidx0:kidx1,:) = laiMax_stat(1:kidx,:)

    ! --- initialize soil_moist_filling: as soon as a we have a full multi-layer soil hydrology soil_moist_filling should be
    !        provided from the calling routine. 
    DO i=1,nsoil
       soil_moist_filling(:,i,:) = relative_moisture(:,:)
    ENDDO
    
    
    IF (debug) THEN
       !! --- Check LAI <= LaiMax_dyn:
       DO i=1,kidx
          k=kidx0+i-1
          DO j=1,ntiles
             IF(lai(i,j) > laiMax_dyn(k,j) + EPSILON(1.0_dp)) THEN
                CALL message("update_phenology", "lai: "//real2string(lai(i,j))//" laiMax: "//real2string(laiMax_dyn(k,j)))
                CALL finish("update_phenology()","ERROR: On input one LAI is larger than maximum LAI.")
             END IF
          END DO
       END DO
    ELSE
       IF (ANY(lai(1:kidx,1:ntiles) > laiMax_dyn(kidx0:kidx1,1:ntiles))) &
            CALL finish("update_phenology","ERROR: LAI larger than maximum LAI on input")
    END IF

    ! --- sum up the day temperatures and NPP-rates

    CALL update_previous_day_variables(air_temp,NPP_rate)

    ! --- recompute pseudo-soil temperature

    CALL update_pseudo_soil_temp(air_temp)

    ! --- update growth phase (needed for evergreen, summergreen, and crop)
    positive_NPP(:) = .false.
    DO i = 1,ntiles
       DO j = 1,nidx
          IF (pheno_type_of_tile(j,i) == 4 .OR. pheno_type_of_tile(j,i) == 5) THEN
             IF (previous_day_NPP_rate(kidx0 + j - 1,i) > 0.0_dp) positive_NPP(j) = .true.
          END IF
       END DO
    END DO
    CALL update_growth_phase(lat,positive_NPP)

    ! --- update maximum LAI (once a year) and if needed reduce also actual LAI

!!$    CALL update_maximum_lai(lat(1:nidx),laiMax_stat(1:nidx,1:ntiles),lai(1:nidx,1:ntiles))

    ! --- Determine new LAI --------------------------

    DO i=1,ntiles

       DO j=1,nidx
          k=kidx0+j-1
          IF(laiMax_dyn(k,i) < PhenParams%all%LAI_negligible%v) then    !! ==> LAI is negligible ==> pheno state is essentialy zero 
             laiMax_dyn(k,i) = PhenParams%all%LAI_negligible%v
             lai(j,i) = PhenParams%all%LAI_negligible%v - EPSILON(1.0_dp)
          END IF
       END DO

       DO j=1,nidx
          k=kidx0+j-1
          SELECT CASE(pheno_type_of_tile(j,i))
             CASE(unvegetated) ! --- bare land ------------------------
                !lai(j,i) =  PhenParams%all%LAI_negligible%v 
                lai(j,i) = 0._dp
             CASE(evergreen)
                IF(ABS(growth_Phase_EG(k)) < EPSILON(1.0_dp)) THEN  !! we are during the growth phase 
                   IF(lai(j,i) < PhenParams%all%laiSeed%v) THEN
                      !! start growth at least with seed value smaller than laiMax
                      lai(j,i) = MIN(laiMax_dyn(k,i)-EPSILON(1.0_dp),PhenParams%all%laiSeed%v)
                   ENDIF
                   IF(lai(j,i) > PhenParams%all%LAI_negligible%v) THEN                      !! otherwise nothing to do
                      xtmp = PhenParams%EG%growthRate%v
                      lai(j,i) = letItGrow(lai(j,i), xtmp, 0.0_dp, laiMax_dyn(k,i))
                   END IF
                ELSE                                                                        !! we are during the vegetative phase
                   xtmp = PhenParams%EG%shedRate_rest%v
                   lai(j,i) = letItDie(lai(j,i),xtmp)
                END IF
             CASE(summergreen) ! --- summergreen -----------------------
                IF(growth_Phase_SG(k) > -0.5_dp) THEN                                       !! we are in the growth/vegetative phase
                   IF(growth_Phase_SG(k) > 0.5_dp) THEN                                     !! we are in the vegetative phase
                      xtmp = PhenParams%SG%shedRate_veget%v
                      lai(j,i) =  letItDie(lai(j,i),xtmp)
                   ELSE                                                                     !! we are in the growth phase
                      IF(lai(j,i) < PhenParams%all%laiSeed%v) THEN
                         !! start growth at least with seed value smaller than laiMax
                         lai(j,i) = MIN(laiMax_dyn(k,i)-EPSILON(1.0_dp),PhenParams%all%laiSeed%v)
                      ENDIF
                      IF(lai(j,i) > PhenParams%all%LAI_negligible%v) THEN                   !! (otherwise nothing to do)
                         xtmp = PhenParams%SG%growthRate%v
                         lai(j,i) = letItGrow(lai(j,i), xtmp, 0.0_dp, laiMax_dyn(k,i))
                      END IF
                   END IF
                ELSE                                                                        !! we are in the rest phase
                   xtmp = PhenParams%SG%shedRate_rest%v
                   lai(j,i) =  letItDie(lai(j,i),xtmp)
                END IF
             CASE(raingreen)  ! --- raingreen ------------------------
                average_filling(j) = SUM(soil_moist_filling(j,1:nsoil,i)) / REAL(nsoil)     !! average filling of the soil buckets
                                                                                           
                !! If LAI dropped below seed value and there is sufficient water, then start growth at least with seed value but
                !! smaller than laiMax. (Note: addition of 1E-3_dp for numerical reasons)
                IF(average_filling(j) > PhenParams%all%wilt_point%v + 1E-3_dp) THEN         
                   IF(lai(j,i) < PhenParams%all%laiSeed%v                                 &
                         .and.                                                            &
                      average_filling(j) > PhenParams%RG%bucketFill_leafout%v ) THEN
                      lai(j,i) = MIN(laiMax_dyn(k,i)-EPSILON(1.0),PhenParams%all%laiSeed%v)
                   ENDIF                                                                    
                   IF(previous_day_NPP_rate(k,i) > 0._dp) THEN                           !! If yesterdays net production is positive
                      grRate = PhenParams%RG%growthRate%v
                      xtmp   = shedRate_RG(average_filling(j))
                      lai(j,i) =  letItGrow(lai(j,i), grRate, xtmp, laiMax_dyn(k,i))     !! start growing.
                   END IF 
                ELSE                                                                     !! In case of bad growth conditions
                   xtmp = PhenParams%RG%shedRate_drySeason%v
                   lai(j,i) = letItDie(lai(j,i), xtmp)                                   !! shed leaves as fast as possible.
                END IF
             CASE(grasses) ! --- grasses ------------------------
                ! Surprisingly the current concept for grasses is more restrictive than that for raingreens!

                !! If temperature is fine and there is sufficient moisture in upper soil layer start growth at least with seed
                !! value, but smaller than laiMax. (Note: addition of 1E-3_dp for numerical reasons)
                IF(air_temp(j) > PhenParams%GRS%crit_temp%v &
                     .AND. &
                      soil_moist_filling(j,1,i) > PhenParams%all%wilt_point%v  + 1E-3_dp &
                   ) THEN
                   IF(lai(j,i) < PhenParams%all%laiSeed%v) THEN
                       lai(j,i) = MIN(laiMax_dyn(k,i)-EPSILON(1.0),PhenParams%all%laiSeed%v)
                   ENDIF 
                   IF( previous_day_NPP_rate(k,i) > 0._dp ) THEN                         !! If yesterdays net production is positive
                      xtmp = PhenParams%GRS%growthRate%v
                      lai(j,i) = letItGrow(lai(j,i), xtmp, 0.0_dp, laiMax_dyn(k,i))      !! let it grow
                   ELSE 
                      xtmp = PhenParams%GRS%shedRate_growth%v
                      lai(j,i) = letItGrow(lai(j,i), 0.0_dp, xtmp, laiMax_dyn(k,i))      !! no growth, only a bit of leaf shedding
                   END IF
                ELSE                                                                     !! otherwise ...
                   xtmp = PhenParams%GRS%shedRate_drySeason%v
                   lai(j,i) = letItDie(lai(j,i),xtmp)                                    !! let it die as fast as possible
                END IF                                                              
             CASE(crops) ! --- crops ------------------------

                IF (ABS(lat(j)) > 30._dp) THEN                                             !! We are outside the tropics

                !! If crops are in the growth phase and there is sufficient water in upper soil layer, then there are good growth
                !! conditions. (Note: addition of 0.001 for numerical reasons)
                IF(growth_phase_CRP(k) > -0.5_dp                                        &
                     .AND.                                                              &
                      soil_moist_filling(j,1,i) > PhenParams%all%wilt_point%v + 1E-3_dp &
                  ) THEN

                   !! If also temperature is fine growth should start at least with seed value ( < laiMax) if upper soil layer 
                   !! has at least a little bit water
                   IF(pseudo_soil_temp(k) > PhenParams%CRP%crit_temp%v) THEN
                      IF(lai(j,i) < PhenParams%all%laiSeed%v .AND.                      &
                         soil_moist_filling(j,1,i) > PhenParams%CRP%sproud%v) THEN
                         lai(j,i) = MIN(laiMax_dyn(k,i)-EPSILON(1.0_dp),PhenParams%all%laiSeed%v)
                      ENDIF

                      !! If also yesterdays net production is positive let it grow
                      IF (previous_day_NPP_rate(k,i) > 0._dp) THEN                        
                         grRate = PhenParams%CRP%leafAlloc_frac%v * specLeafArea(j,i) * &  !! growth rate is high (depending on NPP)
                                    MAX(0.0_dp,NPP_rate(j,i)) * 86400.0_dp / lai(j,i)      !! .. and is in units of 1/days
                         lai(j,i) = letItGrow(lai(j,i), grRate, 0.0_dp, laiMax_dyn(k,i))
                      ELSE                                                                 !! else ..
                         xtmp = PhenParams%CRP%shedRate_growth%v
                         lai(j,i) = letItDie(lai(j,i),xtmp)                                !! .. do a tiny bit of leaf shedding
                      END IF
                   ELSE
                      xtmp = PhenParams%CRP%shedRate_growth%v
                      lai(j,i) = letItDie(lai(j,i),xtmp)                                   !! .. do a tiny bit of leaf shedding
                   END IF
                ELSE                                                                       !! otherwise ...
                   xtmp = PhenParams%CRP%shedRate_rest%v
                   lai(j,i) = letItDie(lai(j,i),xtmp)                                      !! .. shed leaves fast
                END IF

                ELSE                              !! We are in the tropics
                !! If temperature is fine and there is sufficient water in upper soil layer then there are good growth conditions
                !!  (addition of 0.001 for numerical reasons)
                IF(pseudo_soil_temp(k) > PhenParams%CRP%crit_temp%v                     &
                     .AND.                                                              &
                      soil_moist_filling(j,1,i) > PhenParams%all%wilt_point%v + 1E-3_dp &
                  ) THEN
                   !! growth should start at least with with seed value ( < laiMax) if upper soil layer has at least a little bit
                   !! water
                   IF(lai(j,i) < PhenParams%all%laiSeed%v .AND.                         &
                      soil_moist_filling(j,1,i) > PhenParams%CRP%sproud%v) THEN 
                      lai(j,i) = MIN(laiMax_dyn(k,i)-EPSILON(1.0_dp),PhenParams%all%laiSeed%v)
                   ENDIF
                                                            
                    !! If yesterdays net production is positive let it grow
                   IF(previous_day_NPP_rate(k,i) > 0._dp   ) THEN
                      grRate = PhenParams%CRP%leafAlloc_frac%v * specLeafArea(j,i) *    &  !! growth rate 
                                    MAX(0.0_dp,NPP_rate(j,i)) * 86400.0_dp / lai(j,i)      !! .. is in units of 1/days
                      lai(j,i) = letItGrow(lai(j,i), grRate, 0.0_dp, laiMax_dyn(k,i))
                   ELSE                                                                    !! else ..
                      xtmp = PhenParams%CRP%shedRate_growth%v
                      lai(j,i) = letItDie(lai(j,i),xtmp)                                   !! .. do a tiny bit of leaf shedding
                   END IF
                ELSE                                                                       !! otherwise ...
                   xtmp = PhenParams%CRP%shedRate_rest%v
                   lai(j,i) = letItDie(lai(j,i),xtmp)                                      !! .. shed leaves fast
                END IF
                END IF
             CASE default
                ! if program comes here a missing or a non-allowed phenology type is detected
                CALL finish("update_phenology()","Non-allowed phenology type. PROGRAM SHOULD NEVER COME TO THIS POINT!!")
             END SELECT
       END DO
    END DO

  END SUBROUTINE update_phenology


  ! --- update_pseudo_soil_temp() --------------------------------------------------------------------------------------------------
  !
  ! This routine computes a weighted running mean of the air temperature, which is interpreted here as a pseudo soil temperature:
  ! (1)   T_ps(t) = N^(-1) * SUM(t'=-infty,t) T(t')*exp(-(t-t')*delta/tau_soil),
  ! where "T(t)" is the air temperature at time step with number "t", "delta" the length of the time step (in days) and 
  ! "tau_soil" is the characteristic time for loosing the memory of temperature in the soil (also in days; this is a tuning
  ! parameter! The normalization "N" is 
  ! (2)   N = SUM(t'=-infty,t) exp(-(t-t')*delta/tau_soil) = 1/(1 - exp(-delta/tau_soil)).
  ! This normalization constant (called "N_pseudo_soil_temp") is computed during initialization of this phenology module.
  ! Computation of T_ps(t) is performed iteratively: it follows from (1)
  !                    T(t+1)          delta
  ! (3)   T_ps(t+1) =  ------ + exp(- --------) * T_ps(t).
  !                      N            tau_soil
  ! The exponential factor (called F_pseudo_soil_temp) is computed during initialization of this phenology module.
  !
  ! Technically the only effect of this routine is an update of the field "pseudo_soil_temp(:)". The routine has to be called
  ! once every time step for every grid point.
  !
  ! Remark: Instead of air-temperature one could try to use the bottom temperature to compute the pseudo_soil_temp.
  !
  SUBROUTINE update_pseudo_soil_temp(air_temp)
    REAL(dp) :: air_temp(:) ! air temperature at current time step
    INTEGER :: kidx0, kidx1

    kidx0 = kstart
    kidx1 = kend

    pseudo_soil_temp(kidx0:kidx1) = &
         air_temp(1:nidx)/N_pseudo_soil_temp + F_pseudo_soil_temp * pseudo_soil_temp(kidx0:kidx1)

  END SUBROUTINE update_pseudo_soil_temp

  ! --- update_previous_day_variables() -------------------------------------------------------------------------------------------
  !
  ! Updates the mean day temperature (private field "previous_day_temp(:)") from air temperature and the mean day NPP-rate
  ! (private field previous_day_NPP(:,:)
  !
  SUBROUTINE update_previous_day_variables(air_temp,NPP_rate)

    USE mo_time_control, ONLY: delta_time   !! time step in seconds

    REAL(dp) :: air_temp(1:nidx) !! air temperature at current time step in lowest atmospheric layer for each grid point
    REAL(dp) :: NPP_rate(1:nidx,1:ntiles) !! NPP-rate during current time step for each tile (PFT)
    INTEGER :: kidx0, kidx1, i

    kidx0 = kstart
    kidx1 = kend

    ! --- update mean day values

    ! Note that updating is done globally at the same time step, i.e. for different longitudes the updating happens at different 
    ! local times.

    IF(new_day) THEN !! day has changed --> recompute mean, min, max day temperature of previous day
                             !! and reinitialize day_temp_sum(), day_temp_min(), day_temp_max()
       previous_day_temp(kidx0:kidx1) = day_temp_sum(kidx0:kidx1)/time_steps_per_day
       day_temp_sum(kidx0:kidx1) = air_temp(1:nidx)
       previous_day_temp_min(kidx0:kidx1) = day_temp_min(kidx0:kidx1)
       day_temp_min(kidx0:kidx1) = air_temp(1:nidx)
       previous_day_temp_max(kidx0:kidx1) = day_temp_max(kidx0:kidx1)
       day_temp_max(kidx0:kidx1) = air_temp(1:nidx)
       previous_day_NPP_rate(kidx0:kidx1,1:ntiles) = day_NPP_sum(kidx0:kidx1,1:ntiles)/86400.
       day_NPP_sum(kidx0:kidx1,1:ntiles) = NPP_rate(1:nidx,1:ntiles)*delta_time
    ELSE  !! day has not changed
       day_temp_sum(kidx0:kidx1) = day_temp_sum(kidx0:kidx1) + air_temp(1:nidx)
       day_NPP_sum(kidx0:kidx1,1:ntiles) = day_NPP_sum(kidx0:kidx1,1:ntiles) + NPP_rate(1:nidx,1:ntiles)*delta_time
       DO i = 1,nidx
          IF (day_temp_min(kidx0+i-1) > air_temp(i)) day_temp_min(kidx0+i-1) = air_temp(i)
          IF (day_temp_max(kidx0+i-1) < air_temp(i)) day_temp_max(kidx0+i-1) = air_temp(i)
       END DO
    END IF

  END SUBROUTINE update_previous_day_variables

  ! --- update_maximum_lai() ------------------------------------------------------------------------------------------------------
  !
  ! Sums NPP through the year and reduces maximum LAI (field "laiMax_dyn") at points where NPP is negative during a whole year, and 
  ! increases maximum LAI slightly, when when then a year with positive LAI comes. Updates are done on Jan 1 in the northern
  ! hemisphere and on July 1 in thesouthern hemisphere.
  !
  SUBROUTINE update_maximum_lai(lat,laiMax_stat,lai)
    USE mo_time_control,ONLY: current_date !! the current date
    USE mo_time_control,ONLY: get_year_day !! returns the day in a year, where the seconds make up the fractional part

    REAL(dp), INTENT(in)    :: lat(:)           !! latitude of grid points
    REAL(dp), INTENT(in)    :: laiMax_stat(:,:) !! maximum static LAI
    REAL(dp), INTENT(inout) :: lai(:,:)         !! actual leaf area index (eventually reduced if laiMax_stat gets smaller than lai)

    !! locals

    INTEGER :: day_in_year
    INTEGER :: k,j,n,kidx0,kidx1

    kidx0 = kstart
    kidx1 = kend

    IF(new_day) THEN !! day has changed --> update yearly_NPP_sum() is needed
       day_in_year = int(get_year_day(current_date) + EPSILON(1._dp))
       DO n=1,ntiles

          ! sum NPP
          WHERE(year_NPP_sum(kidx0:kidx1,n) < LARGE)                                    !! Where NPP summation has started,
             year_NPP_sum(kidx0:kidx1,n) = year_NPP_sum(kidx0:kidx1,n) + previous_day_NPP_rate(kidx0:kidx1,n) !! .. sum it!
          END WHERE

          ! check whether maximum LAI has to be changed 

          IF(day_in_year == 1) THEN                                              !! It is Jan. 1,
             DO j=1,nidx
                k=kidx0+j-1
                IF(year_NPP_sum(k,n) < LARGE .AND. lat(j) >= 0._dp) THEN         !! NPP summation has started in northern hemisphere
                   IF(year_NPP_sum(k,n) <= 0._dp) THEN                              !! If yearly NPP sum is negative
                      laiMax_dyn(k,n) = 0.8_dp * laiMax_dyn(k,n)                        !! .. then reduce laiMax
                      lai(j,n) = MIN(lai(j,n),laiMax_dyn(k,n) - EPSILON(1.0_dp))        !! .. and eventually also the actual LAI
                   ELSE                                                                 !! else increase LAImax slightly
                      laiMax_dyn(k,n) = min(1.1_dp * laiMax_dyn(k,n),laiMax_stat(j,n))  !! .. but not above maximum static LAI
                   END IF
                   year_NPP_sum(kidx0:kidx1,n) = 0._dp                                  !! .. and reinitialize NPP sum by zero.
                END IF
             END DO
          END IF

          IF(day_in_year == 181) THEN                                            !! It is  Jul. 1 (ignoring leap years)
             DO j=1,nidx
                k=kidx0+j-1
                IF(year_NPP_sum(k,n) < LARGE .AND. lat(j) < 0._dp) THEN          !! NPP summation has started in northern hemisphere
                   IF(year_NPP_sum(k,n) <= 0._dp) THEN                                  !! .. yearly NPP sum is negative
                      laiMax_dyn(k,n) = 0.8_dp * laiMax_dyn(k,n)                        !! .. then reduce laiMax
                      lai(j,n) = MIN(lai(j,n),laiMax_dyn(k,n) - EPSILON(1.0_dp))        !! .. and eventually also the actual LAI
                   ELSE                                                                 !! else increase LAImax slightly
                      laiMax_dyn(k,n) = MIN(1.1_dp * laiMax_dyn(k,n),laiMax_stat(j,n))  !! .. but not above maximum static LAI
                   END IF
                   year_NPP_sum(kidx0:kidx1,n) = 0._dp                                  !! .. and reinitialize NPP sum by zero.
                END IF
             end do
          end if

          ! Check start of NPP sum
          IF(day_in_year == 1) THEN                                                     !! It is Jan. 1,
             WHERE(year_NPP_sum(kidx0:kidx1,n) > LARGE .AND. lat(kidx0:kidx1) >= 0._dp) !! .. and NPP summation has not started ..
                year_NPP_sum(kidx0:kidx1,n) = 0._dp                                     !! .. then initialize NPP sum by zero.
             END WHERE
          END IF
          if(day_in_year == 181) THEN                                                   !! It is Jul. 1 (ignoring leap years)
             WHERE(year_NPP_sum(kidx0:kidx1,n) > LARGE .and. lat(kidx0:kidx1) < 0._dp)  !! .. and NPP summation has not started ..
                year_NPP_sum(kidx0:kidx1,n) = 0._dp                                     !! .. then initialize NPP sum by zero.
             END WHERE
          END IF
       END DO
    END IF

  END SUBROUTINE update_maximum_lai

  ! --- update_growth_phase() ------------------------------------------------------------------------------------------------------
  !
  ! This subroutine contains the two models for determining dates of budburst ("spring event") and leaf fall ("autumn event") for
  ! the summergreen and evergreen phenologies. For operating the dynamical phenology it needs only be known in what "phase" it
  ! currently is, i.e. whether the current time step falls into the time from spring to autumn event ("growth"), or into the time
  ! from autumn to spring event ("rest") --- the budburst and leaf fall dates itself are not needed. Hence, this routine does not
  ! provide these dates. Instead, it updates the field "growth_phase": A positive entry means that at the particular grid
  ! cell the evergreens and summergreens are growing, whereas a negative entry means that they are at rest.
  ! 
  ! For the spring event the "alternating model" of Murray et al. [1] is used: 
  ! Let "S(d)" denote the value of the heat sum at day "d":
  ! (1) S(d) = SUM(d'=d0,d) MAX(T(d')-T_alt,0),
  ! where "T(d)" is the mean day temperature at day "d", "T_alt" is the "alternating temperature" (which has the function of a
  ! cutoff temperature in the heat sum), and "d0" is the starting date for temperture summation. This starting date is determined
  ! by a temperature criterion (see below) --- once  more we need not know the date, but only whether summation has started: 
  ! this information is kept also in the field "heat_sum", which is set to a value smaller than -1 during times where no 
  ! heat summation takes place. 
  ! Another key quantity of the alternating model the number of chill days "C(d)": this is the number of days with a mean day 
  ! temperature below the alternating temperature "T_alt", where counting is started here at the day "d_a" of the last autumn event:
  ! (2) C(d) = SUM(d'=d_a,d) STEP(T(d)-T_alt);
  ! here STEP() is the Heaviside step function. From C(d) a critical heatsum "S_crit(d)" is computed:
  ! (3) S_crit(d) = S_crit_min + S_crit_range * exp(-C(d)/C_decay),
  ! where "S_crit_min", "S_crit_range" and "C_decay" are parameters: "S_crit_min" and "S_crit_range" define minimum value and
  ! maximum range of the critical heatsum, whereas "C_decay" determines how fast "S_crit(d)" decreases with increasing number of
  ! chill days. Finally, the spring event happens when first
  ! (4)  S(d) >= S_crit(d).
  ! Technically, at this date for the considered grid point the associated entry of the array growth_phase() is set to +1 and 
  ! "springEvent_flag" is set to "true" NO!!. Moreover the chill days count is reset to zero and heat_sum set to -99 to indicate
  ! that summation has stopped.
  !
  ! The "autumn event" is calculated from a pseudo soil temperature: It happens when during the growth phase first the pseudo soil 
  ! temperature Ts(d) falls below the critical soil temperature "Ts_crit" (variable "autumn_event_temp"); to prevent that this
  ! event is detected in spring the condition is added that the mean day air temperature T(d) is smaller than the soil temperature, 
  ! i.e. the autum event happens, when first
  ! (5) T(d) < Ts(d) < Ts_crit 
  ! At this event for the considered grid point the associated entries of the array growth_phase() are set to "-1".
  !  --- An additional daylength-criterion has not implemented yet!
  ! 
  ! It remains to determine the date "d0" for the start of heat summation: The idea is that heat summation starts when
  ! first the pseudo soil temperature during the rest phase gets larger than a critical temperature "heat_summation_start_temp"
  ! (a tuning parameter). Technically this means to initialize "heat_sum" at the particular grid point by zero; since this
  ! value is larger than -1 it indicates that heat summation has started. 
  !
  ! This routine requires that the soil temperature is updated before calling it.
  !
  !
  ! [1] M.B. murray, M.G.R. Cannell and R.I. Smith, "Date of budburst of fifteen tree species in Britain following climatic 
  !     warming", J. Appl. Ecology 26 (1989) 693-700). See also: A. Botta, N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray,
  !     "A global prognostic scheme of leaf onset using satellite data", Global Change Biol. 6 (2000) 709-725.
  !
  SUBROUTINE update_growth_phase(lat,positive_NPP)

    USE mo_time_control, ONLY: current_date, get_year_day

    REAL(dp), INTENT(in)          :: lat(:)                    !! latitude of grid points(1:nidx)
    LOGICAL, INTENT(in)           :: positive_NPP(:)           !! flag to indicate the existence of a grass or a crop type with
                                                               !! positive NPP (1:nidx)

    INTEGER  ::  day_of_year   ! contains the day number, counted from first of january of current year
    INTEGER  ::  k,n,kidx0,kidx1
    REAL(dp) ::  temp_min_max

    kidx0 = kstart
    kidx1 = kend

    day_of_year = INT(get_year_day(current_date))

    ! --- determine where heat summation for evergreen has to start

    DO k=kidx0,kidx1
       n=k - kidx0 +1
       IF( growth_phase_EG(k) < -0.5_dp                & ! we are during rest phase ..(Note: there is no vegetative phase for EG) 
           .AND.                                       & ! and .. 
           heat_sum_EG(k) < -1._dp                     & ! heat summation has not yet started
           .AND.                                       & ! and ..
           (                                           & ! ..
             (lat(n) >= 0._dp .AND. day_of_year == 1)  & ! we are at northern hemisphere at January 1 ..
             .OR.                                      & ! or ..
             (lat(n) < 0._dp .and. day_of_year == 183) & ! we are at southern hemisphere at July 2 (July 1 in leap years) .. 
           )                                           & ! ..
         ) THEN                                          ! --> heat summation has to start, ..
           heat_sum_EG(k) = 0._dp                        ! so initialize it by zero
       END IF
    END DO

    ! --- determine where heat summation for summergreen has to start

    DO k=kidx0,kidx1
       n= k- kidx0 + 1
       IF( growth_phase_SG(k) < -0.5_dp                & ! we are during rest .. 
           .AND.                                       & ! and ..                
           heat_sum_SG(k) < -1._dp                     & ! heat summation has not yet started
           .AND.                                       & ! and ..
           (                                           & ! ..
             (lat(n) >= 0._dp .AND. day_of_year == 1)  & ! we are at northern hemisphere at January 1 ..
             .OR.                                      & ! or ..
             (lat(n) < 0._dp .AND. day_of_year == 183) & ! we are at southern hemisphere at July 2 (July 1 in leap years) .. 
           )                                           & ! ..
         ) THEN                                          ! --> heat summation has to start, ..
           heat_sum_SG(k) = 0._dp                        ! so initialize it by zero
       END IF
    END DO

    IF(new_day) THEN                                       ! only at first time step of a day

    ! --- determine where heat summation for extra-tropical crops starts

       DO k=kidx0,kidx1
          n= k- kidx0 + 1
                                                                ! At the beginning of spring:
          IF(lat(n) >= 0._dp .AND. day_of_year == 70          & ! we are at northern hemisphere at March 11 (March 10 in leap years)
             .OR.                                             & ! or ..
             lat(n) < 0._dp .AND. day_of_year == 252          & ! we are at southern hemisphere at September 9 (September 8 in 
             ) THEN                                             ! leap years) ..
             heat_sum_CRP(k) = 0._dp                            ! --> heat summation has to start, so initialize it by zero
             IF (ABS(growth_phase_CRP(k)) < 1.5_dp) THEN        ! if there are no winter crops ..
                growth_phase_CRP(k) = -1._dp                    ! .. crop phase is set to rest phase of summer crops
             END IF
          END IF
                                                                ! At the beginning of summer:
          IF(lat(n) >= 0._dp .AND. day_of_year == 172         & ! we are at northern hemisphere at June 21 (June 20 in leap years)
             .OR.                                             & ! or ..
             lat(n) < 0._dp .AND. day_of_year == 354          & ! we are at southern hemisphere at December 21 (December 20 in leap
             ) THEN                                             ! years)

             ! If there are no winter crops and crops already have been harvested try double cropping.
             IF (ABS(growth_phase_CRP(k)) < 1.5_dp .AND.      &
                 heat_sum_CRP(k) > PhenParams%CRP%heat_sum_harvest%v) THEN
                heat_sum_CRP(k) = 0._dp                               ! then heat summation has to start, so initialize it by zero
                growth_phase_CRP(k) = 1._dp                           ! and crop phase is set to growth phase of double cropping
             ELSE IF (growth_phase_CRP(k) > 1.5_dp) THEN              ! if there are winter crops still in the growing phase ..
                growth_phase_CRP(k) = -2._dp                          ! .. then switch to the rest phase of winter crops 
             END IF
             IF (heat_sum_winter(k) > 0._dp)                        & ! Heat sum for winter is multiplied by -1 to indicate ..
                heat_sum_winter(k) = -1._dp * heat_sum_winter(k)      ! .. that there is no heat summation during summer
          END IF
                                                                ! In autumn:
          IF(lat(n) >= 0._dp .AND. day_of_year == 289         & ! we are at northern hemisphere at October 16 (October 15 in leap
             .OR.                                             & ! years) or ..
             lat(n) < 0._dp .AND. day_of_year == 105          & ! we are at southern hemisphere at April 15 (April 14 in leap years)
             ) THEN

             ! If there are winter crops and they were not harvested last winter or spring and summer crops would have been
             ! harvested already  --> crop phase is set to rest phase of summer crops
             IF (ABS(growth_phase_CRP(k)) > 1.5_dp .AND.      &
                ABS(heat_sum_winter(k)) < PhenParams%CRP%heat_sum_harvest%v .AND. &
                heat_sum_CRP(k) > PhenParams%CRP%heat_sum_harvest%v) THEN   
                growth_phase_CRP(k) = -1._dp
             END IF

             ! If there are winter crops and summer crops would have been harvested already twice
             ! --> crop phase is set to rest phase of summer crops
             END IF
             IF (ABS(growth_phase_CRP(k)) > 1.5_dp .AND. &
                heat_sum_CRP(k) > 2.0_dp * PhenParams%CRP%heat_sum_harvest%v) THEN
                growth_phase_CRP(k) = -1._dp

             ! If there is no double cropping and summer crops are still not harvested and winter crops would have been harvested
             ! last spring and it is still warm then try winter crops --> crop phase is set to rest phase of winter crops
             IF (ABS(growth_phase_CRP(k) - 1._dp) > 0.5_dp .AND. &
                 heat_sum_CRP(k) < PhenParams%CRP%heat_sum_harvest%v .AND. &
                 ABS(heat_sum_winter(k)) > PhenParams%CRP%heat_sum_harvest%v .AND. &
                 pseudo_soil_temp(k) > PhenParams%CRP%crit_temp%v) THEN
                 growth_phase_CRP(k) = -2._dp
             END IF
             heat_sum_winter(k) = EPSILON(1._dp) ! heat summation for winter crops has to start anyway
          END IF
       END DO

       ! --- update growth_phase

       ! --- update counter for days since last begin of growth phase

       DO k=kidx0,kidx1
          IF(days_since_growth_begin_EG(k) > 0.0_dp .AND. days_since_growth_begin_EG(k) < 365.0_dp) THEN
             days_since_growth_begin_EG(k) = days_since_growth_begin_EG(k) + 1.0_dp
          END IF
          IF(days_since_growth_begin_SG(k) > 0.0_dp .AND. days_since_growth_begin_SG(k) < 365.0_dp) THEN
             days_since_growth_begin_SG(k) = days_since_growth_begin_SG(k) + 1.0_dp
          END IF
       END DO

       ! --- count chill days for evergreens or update heat sum
       
       DO k=kidx0,kidx1
          IF(heat_sum_EG(k) > -1.0_dp ) THEN                                            ! heat summation has started
             IF(pseudo_soil_temp(k) < PhenParams%EG%alternation_temp%v ) THEN           ! the day is a chill day ..
                chill_days_EG(k) = chill_days_EG(k) + 1._dp                             ! --> increment chill days counter
                chill_days_EG(k) = &                                                    ! limit the number of chill days to prevent
                     MIN(REAL(chill_days_EG(k)),REAL(PhenParams%EG_SG%max_chill_days%v))! values beyond any limit in polar regions
             ELSE                                                                       ! increase heat sum:
                heat_sum_EG(k) = heat_sum_EG(k) +                                     & ! --> add mean day temperature excess
                     pseudo_soil_temp(k)-PhenParams%EG%alternation_temp%v               !     above alternating temperature
             END IF
          END IF
       END DO

       ! --- count chill days for summergreens or increase heat sum

       DO k=kidx0,kidx1
          IF(heat_sum_SG(k) > -1._dp) THEN                                              ! heat summation has started ..
             IF(pseudo_soil_temp(k) < PhenParams%SG%alternation_temp%v) THEN            ! the day is a chill day ..
                chill_days_SG(k) = chill_days_SG(k) + 1._dp                             ! --> increment chill days counter
                chill_days_SG(k) = &                                                    ! limit the number of chill days to prevent
                     MIN(REAL(chill_days_SG(k)),REAL(PhenParams%EG_SG%max_chill_days%v))! values beyond any limit in polar regions
             ELSE                                                                       ! increase heat sum:
                heat_sum_SG(k) = heat_sum_SG(k) +                                     & ! --> add mean day temperature excess
                     pseudo_soil_temp(k)-PhenParams%SG%alternation_temp%v               !     above alternating temperature
             END IF
          END IF
       END DO

       ! --- check begin and end of growth phase for evergreens (no rest phase for evergreens)

       DO k=kidx0,kidx1
          IF(growth_phase_EG(k) > -0.5_dp ) THEN                                        ! evergreens are during growth
             IF(days_since_growth_begin_EG(k) > PhenParams%EG%growthPhaseLength%v) THEN ! growth phase of EGs has ended
                growth_phase_EG(k) = -1.0_dp                                            ! --> rest phase starts
             END IF
          ELSE                                                                          ! everergreen are in rest phase
             IF(heat_sum_EG(k) > -1._dp                                               & ! if heat summation has started ..
                  .AND.                                                               & ! and ..
                heat_sum_EG(k) - PhenParams%EG%heat_sum_min%v >                       & ! criticality condition is fulfilled ..
                PhenParams%EG%heat_sum_range%v * EXP(-chill_days_EG(k) /              & ! ..
                PhenParams%EG%chill_decay_const%v)) THEN                                ! ..
                growth_phase_EG(k) = 0._dp                                              ! --> growth phase starts
                days_since_growth_begin_EG(k) = 1._dp                                   ! --> reinitialize counter for days since
                                                                                        !     growth phase begin
                heat_sum_EG(k) = -99._dp                                                ! --> set heat_sum to a value smaller than
                                                                                        !     -1 to indicate that heat summation is
                                                                                        !     no more needed
                chill_days_EG(k) = 0._dp                                                ! --> reset chill days count 
             END IF
          END IF
       END DO


       ! --- check begin and end of growth and vegetative phase for summergreens

       DO k=kidx0,kidx1

          ! check begin of growth phase
          IF(growth_phase_SG(k) < -0.5_dp ) THEN                               ! summergreens are during rest
             IF(heat_sum_SG(k) > -1._dp                                      & ! if heat summation has started ..
                  .AND.                                                      & ! and ..
                heat_sum_SG(k) - PhenParams%SG%heat_sum_min%v >              & ! criticality condition is fulfilled ..
                PhenParams%SG%heat_sum_range%v * EXP(-chill_days_SG(k) /     & ! ..
                PhenParams%SG%chill_decay_const%v)) THEN                       ! ..
                growth_phase_SG(k) = 0._dp                                     ! --> growth phase starts
                days_since_growth_begin_SG(k) = 1.0_dp                         !     start counting days since begin of growth phase
                heat_sum_SG(k) = -99._dp                                       ! --> set heat_sum to a value smaller than -1 to
                                                                               !     indicate that heat summation is no more needed
                chill_days_SG(k) = 0._dp                                       ! --> reset chill days count 
             END IF
          ELSE !! Summergreens are in growth or vegetative phase
             IF(growth_phase_SG(k) < 0.5_dp) THEN                                       ! Summergreens are in growth phase
             IF(days_since_growth_begin_SG(k) > PhenParams%SG%growthPhaseLength%v) THEN ! growth phase of EGs has ended
                   growth_phase_SG(k) = 1.0_dp                                          ! --> vegetative phase starts
                END IF
             ELSE                                                                  ! summergreen are in vegetative phase
                IF(                                                              & ! Check end of vegetative phase
                   previous_day_temp(k) < pseudo_soil_temp(k)                    & ! mean day temperature is smaller than soil temp.
                     .AND.                                                       & ! and ..
                   pseudo_soil_temp(k) < PhenParams%SG%autumn_event_temp%v       & ! soil temperature falls below critical ..
                     .AND.                                                       & ! and ..
                   days_since_growth_begin_SG(k)                                 & ! ..
                                  > PhenParams%SG%minLength_gvPhase%v) THEN        ! growth + vegetative phase have lasted the 
                                                                                   ! minimum time ..
                   growth_phase_SG(k) = -1._dp                                     ! --> rest phase starts
                ELSE                                  
                   IF(days_since_growth_begin_SG(k)                                &
                                     > PhenParams%SG%maxLength_gvPhase%v) THEN     ! if growth + vegetative phase has reached
                                                                                   ! maximum length then 
                      growth_phase_SG(k) = -1._dp                                  ! --> rest phase starts
                   END IF
                END IF
             END IF
          END IF
       END DO

       ! --- set growth phase for crops and increase heat sum

       DO k=kidx0,kidx1
          n= k- kidx0 + 1
          temp_min_max = (previous_day_temp_min(k) + previous_day_temp_max(k)) / 2.0_dp
          IF(temp_min_max > PhenParams%CRP%gdd_temp%v .AND.               & ! it is warm enough to count gdd ..
             positive_NPP(n)) THEN                                          ! and yesterdays productivity of a crop or grass type is
                                                                            ! positive
             heat_sum_CRP(k) = heat_sum_CRP(k) +                          & ! --> add mean day temperature excess above critical
                  temp_min_max - PhenParams%CRP%gdd_temp%v                  !     temperature to the heat sum for summer crops
             IF (heat_sum_winter(k) > 0._dp) THEN                           ! additionally, if we are in the growing period of 
                                                                            ! winter crops
                heat_sum_winter(k) = heat_sum_winter(k) +                 & ! --> add mean day temperature excess above critical
                     temp_min_max - PhenParams%CRP%gdd_temp%v               !     temperature to the heat sum for winter crops
             END IF
          END IF
          IF(pseudo_soil_temp(k) > PhenParams%CRP%crit_temp%v) THEN            ! it is warm enough for crop growth ..
             IF (ABS(growth_phase_CRP(k)) < 1.5_dp) THEN                       ! .. and if there are no winter crops ..
                IF (heat_sum_CRP(k) > PhenParams%CRP%heat_sum_harvest%v) THEN  ! .. and the heat sum for summer crops is high enough
                   growth_phase_CRP(k) = -1._dp                                ! --> crops are harvested,  i.e. rest phase starts,
                ELSE                                                           ! otherwise ..
                   IF (growth_phase_CRP(k) < 0.5_dp) growth_phase_CRP(k) = 0._dp  ! .. crops are still in the growth phase ..
                END IF
             ELSE                                                                 ! .. but if there are winter crops and
                IF (heat_sum_winter(k) > PhenParams%CRP%heat_sum_harvest%v .OR. & ! .. the heat sum for winter crops is high enough
                    heat_sum_winter(k) < 0._dp) THEN                              ! .. or their rest period (summer) started
                   growth_phase_CRP(k) = -2._dp                                   ! --> crops are harvested, i.e. rest phase starts,
                ELSE                                                              ! otherwise ..
                   growth_phase_CRP(k) = 2._dp                                    ! .. crops are still in the growth phase
                END IF
             END IF
          END IF
       END DO

    END IF
  END SUBROUTINE update_growth_phase

  ! --- letItGrow() ---------------------------------------------------------------------------------------------------------------
  !
  ! Performs one time step for the growth of the phenology state. 
  ! The growth model for the phenology state (x=LAI/LaiMax_dyn) involves two elements: logistic growth and exponential leaf
  ! shedding. Logistic growth guarantees that growth stops at the carrying capacity (here equal 1, because this is the maximum
  ! value of the pheno state). Together with a term accounting for leaf shedding the growth equation is
  !            dx
  ! (1)       ---- = k*(1-x)*x -p*x,
  !            dt
  ! where "k'" is the growth rate scaled by the maximum LAI
  !
  !                                             x(t)
  ! (2)      x(t+tau) = (k-p) ----------------------------------------.
  !                            k*x(t) + exp(-(k-p)*tau)*(k-p-k*x(t))
  !
  REAL(dp) FUNCTION letItGrow(x,k,p,z)  !! returns result of integration x(t+tau)

    REAL(dp),INTENT(IN) :: x            !! LAI at time t (0<= x <= 1)
    REAL(dp),INTENT(IN) :: k            !! growth rate          (units: 1/days)
    REAL(dp),INTENT(IN) :: p            !! leaf shedding rate   (units: 1/days)
    REAL(dp),INTENT(IN) :: z            !! maximum LAI
    REAL(dp)            :: hlp1,hlp2,numerator,denominator

    hlp1 = k-p
    hlp2 = k * x / z   ! x/z is the phenoState
    ! Note: make sure that x <= z on input. Call to finish in case of x>z  has been removed from here to allow for inlining.

    numerator = hlp1
    denominator = hlp2 + EXP(-hlp1*timeStep_in_days) * (hlp1 - hlp2)
    IF(ABS(ABS(numerator) - ABS(denominator)) <= 100._dp * EPSILON(1.0_dp)) THEN !! This prevents zero/zero divisions
       letItGrow = x
    ELSE
       letItGrow = x * numerator / denominator
    ENDIF
  END FUNCTION letItGrow
    
  ! --- letItDie() ----------------------------------------------------------------------------------------------------------------
  !
  ! This is the function letItGrow() for growth rate k=0, i.e. it describes only leaf shedding:
  !      x(t+tau) = x(t)*exp(-p*tau), 
  ! where "tau" is the time step and "p" the leaf shedding rate. This function is introduced only to save computation time.
  !
  REAL(dp) elemental FUNCTION letItDie(x,p)  !! returns x(t+tau)

    REAL(dp),INTENT(IN) :: x                 !! pheno state at time t (0<= x <= 1)
    REAL(dp),INTENT(IN) :: p                 !! leaf shedding rate (units: 1/days)

    letItDie = EXP(-p * timeStep_in_days) * x

  END FUNCTION letItDie

  ! --- shedRate_RG ---------------------------------------------------------------------------------------------------------------
  !
  ! Computes the shedding rate for the raingreens during wet season as a function of the soil water gauge
  ! The idea is, that for high water availability (bucket filling > bucketFill_critical) leafs are shedded only because of their
  ! aging. Below this critical value the shedding rate increases because the plants adapt to lower water availability. At the 
  ! wilting point the shedding is so large, that it is larger than the assumed fixed growth rate.
  !
  REAL(dp) elemental FUNCTION shedRate_RG(bucket_filling)
    REAL(dp),INTENT(IN) :: bucket_filling     !! average filling of the soil water buckets

    REAL(dp) shedRate,slope

    slope = PhenParams%RG%growthRate%v / (PhenParams%RG%bucketFill_critical%v - PhenParams%all%wilt_point%v) 
    shedRate = PhenParams%RG%growthRate%v + PhenParams%RG%shedRate_aging%v - slope * (bucket_filling - PhenParams%all%wilt_point%v)
    shedRate_RG = MIN(MAX(PhenParams%RG%shedRate_aging%v,shedRate),PhenParams%RG%growthRate%v+PhenParams%RG%shedRate_aging%v)
    
  END FUNCTION shedRate_RG

  ! --- set_rdata() ---------------------------------------------------------------------------------------------------------------
  !
  ! Sets the values of a structure of type rdata_type
  !
  FUNCTION set_rdata(value,name,ubound,lbound,units)
    TYPE(rdata_type)    :: set_rdata
    REAL(dp),intent(in)     :: value
    CHARACTER(len=*),INTENT(in) :: name
    REAL(dp)   ,INTENT(in),optional :: ubound     !! If present: The upper bound for value
    REAL(dp)   ,INTENT(in),optional :: lbound     !! If present: The lower bound for value
    CHARACTER(len=*),INTENT(in),optional :: units !! If present: the units of value

    set_rdata%n = name
    set_rdata%v = value

    IF(present(ubound)) THEN
       IF(value > ubound) THEN
!          CALL finish("set_rdata","ERROR: variable "//trim(name)//" = "//real2string(real(value))//&
!                                  " is larger than upper bound = "//real2string(real(ubound)))
       END IF
       set_rdata%ubound%exists = .TRUE.
       set_rdata%ubound%value = ubound
    ELSE
       set_rdata%ubound%exists = .FALSE.
       set_rdata%ubound%value = NaN
    END IF

    IF(present(lbound)) THEN
       IF(value < lbound) THEN
!          CALL finish("set_rdata","ERROR: variable "//trim(name)//" = "//real2string(real(value))//&
!                                  " is smaller than lower bound = "//real2string(real(lbound)))
       END IF
       set_rdata%lbound%exists = .TRUE.
       set_rdata%lbound%value = lbound
    ELSE
       set_rdata%lbound%exists = .FALSE.
       set_rdata%lbound%value = NaN
    END IF

    IF(present(units)) THEN
       set_rdata%units = trim(adjustl(units))
    ELSE
       set_rdata%units = ""
    END IF
  END FUNCTION set_rdata

  ! --- init_pheno_params() --------------------------------------------------------------------------------------------------------
  !
  ! Initializes a variable of type PhenologyParameters_type by setting all values to default parameters.  Here also unique verbose
  ! names of the parameters are defined. Whenever a variable of type PhenologyParameters_type is used, it should first be
  ! initialized by this routine. This routine also allocates memory for the parameters.
  !
  SUBROUTINE init_pheno_params(params)
    TYPE(PhenologyParameters_type),intent(inout) :: params

    INTEGER :: status
    INTEGER :: no

    !! Allocate memory if needed

    IF(.NOT. associated(params%memory)) THEN
       ALLOCATE(params%memory(1:noOfPhenParams),STAT=status)
       IF(status /= 0) THEN
          CALL finish("connectParams2memory()",&
               "ERROR: Allocation of params%memory(1:"//trim(int2string(noOfPhenParams))//") failed.")
       END IF

       IF (debug) CALL message("init_pheno_params()", &
                  "allocation of memory(1:"//trim(int2string(noOfPhenParams))//") for phenology parameters")

    END IF

    !! Start setting of default parameters

    no=0

    !! set general parameters

    no=no+1; params%all%LAI_negligible  => params%memory(no)
             params%all%LAI_negligible  =  set_rdata(1.0e-05_dp,"LAI_negligible",LBOUND=0._dp,UBOUND=0.01_dp,UNITS="--")
    no=no+1; params%all%laiSeed         => params%memory(no)
             params%all%laiSeed         =  set_rdata(4.0e-01_dp,"laiSeed",LBOUND=0.005_dp,UBOUND=1.0_dp,UNITS="--")
    no=no+1; params%all%wilt_point      => params%memory(no)
             params%all%wilt_point      =  set_rdata(0.35_dp,"wilt_point",LBOUND=0.34999_dp,UBOUND=0.35001_dp,UNITS="--")

    !! set parameters common to evergreen and summergreen

    no=no+1; params%EG_SG%tau_pseudo_soil  => params%memory(no)
             params%EG_SG%tau_pseudo_soil  =  set_rdata(10.0_dp,"tau_pseudo_soil",LBOUND=3.0_dp,UBOUND=100.0_dp,UNITS="days")
    no=no+1; params%EG_SG%max_chill_days   => params%memory(no)
             params%EG_SG%max_chill_days   = set_rdata(365._dp,"max_chill_days",LBOUND=0._dp,UBOUND=366.0_dp,UNITS="days")

    !! Parameters of evergreen only

    no=no+1; params%EG%alternation_temp => params%memory(no)
             params%EG%alternation_temp = set_rdata(4.0_dp,"EG%alternation_temp",LBOUND=0.0_dp,UBOUND=15.0_dp,UNITS="Celsius")
    no=no+1; params%EG%heat_sum_min     => params%memory(no)
             params%EG%heat_sum_min     =  set_rdata(10.0_dp,"EG%heat_sum_min",LBOUND=0.5_dp,UBOUND=200.0_dp,UNITS="degree days")
    no=no+1; params%EG%heat_sum_range   => params%memory(no)
             params%EG%heat_sum_range   &
                   = set_rdata(1.5e+02_dp,"EG%heat_sum_range",LBOUND=0._dp,UBOUND=2400._dp,UNITS="degree days")
    no=no+1; params%EG%chill_decay_const => params%memory(no)
             params%EG%chill_decay_const &
                   = set_rdata(1.5e+01_dp,"EG%chill_decay_const",LBOUND=0._dp,UBOUND=365._dp,UNITS="days")
    no=no+1; params%EG%growthPhaseLength => params%memory(no)
             params%EG%growthPhaseLength  &
                   = set_rdata(60.0_dp,"EG%growthPhaseLength",LBOUND=3._dp,UBOUND=180._dp,UNITS="days")
    no=no+1; params%EG%shedRate_rest     => params%memory(no)
             params%EG%shedRate_rest     &
                   = set_rdata(8.0e-04_dp,"EG%shedRate_rest",LBOUND=1._dp / 1460._dp,UBOUND=1._dp / 90._dp,UNITS="1/days")
    no=no+1; params%EG%growthRate        => params%memory(no)
             params%EG%growthRate        =  set_rdata(1.5e-02_dp, "EG%growthRate", LBOUND=0.01_dp,UBOUND=0.33_dp,UNITS="1/days")

    !! Parameters of summergreen only
    no=no+1; params%SG%alternation_temp => params%memory(no)
             params%SG%alternation_temp =  set_rdata( 4.0_dp,"SG%alternation_temp",LBOUND=0.0_dp,UBOUND=15.0_dp,UNITS="Celsius")
    no=no+1; params%SG%heat_sum_min     => params%memory(no)
             params%SG%heat_sum_min     =  set_rdata(30.0_dp,"SG%heat_sum_min",LBOUND=0.1_dp,UBOUND=200.0_dp,UNITS="degree days")
    no=no+1; params%SG%heat_sum_range   => params%memory(no)
             params%SG%heat_sum_range   = set_rdata(2.0e+02_dp,"SG%heat_sum_range",LBOUND=0._dp,UBOUND=2400._dp,UNITS="degree days")
    no=no+1; params%SG%chill_decay_const  => params%memory(no)
             params%SG%chill_decay_const  = set_rdata(2.5e+01_dp,"SG%chill_decay_const",LBOUND=0._dp,UBOUND=365._dp,UNITS="days")
    no=no+1; params%SG%growthPhaseLength  => params%memory(no)
             params%SG%growthPhaseLength  &
                   = set_rdata(60._dp,"SG%growthPhaseLength",LBOUND=3._dp,UBOUND=90._dp,UNITS="days")
    no=no+1; params%SG%autumn_event_temp  => params%memory(no)
             params%SG%autumn_event_temp  =  set_rdata(10._dp,"SG%autumn_event_temp",LBOUND=0._dp,UBOUND=15._dp,UNITS="Celsius")
    no=no+1; params%SG%shedRate_veget   => params%memory(no)
             params%SG%shedRate_veget   &
                   = set_rdata(4.0e-03_dp,"SG%shedRate_veget",LBOUND=1._dp / 730._dp,UBOUND=1._dp / 90._dp,UNITS="1/days")
    no=no+1; params%SG%shedRate_rest      => params%memory(no)
             params%SG%shedRate_rest      =  set_rdata(1.0e-01_dp,"SG%shedRate_rest",LBOUND=0._dp,UBOUND=10._dp,UNITS="1/days")
    no=no+1; params%SG%growthRate         => params%memory(no)
             params%SG%growthRate         =  set_rdata(8.7195e-02_dp,"SG%growthRate",LBOUND=0.035_dp,UBOUND=0.46_dp,UNITS="1/days")
    no=no+1; params%SG%maxLength_gvPhase  => params%memory(no)
             params%SG%maxLength_gvPhase  =  set_rdata(270._dp,"SG%maxLength_gvPhase", &
                                                          LBOUND=params%SG%growthPhaseLength%v,UBOUND=365._dp,UNITS="days")
    no=no+1; params%SG%minLength_gvPhase  => params%memory(no)
             params%SG%minLength_gvPhase  =  set_rdata(params%SG%growthPhaseLength%v + 30._dp,  &
                                                       "SG%minLength_gvPhase",                          &
                                                        LBOUND=params%SG%growthPhaseLength%v,   &
                                                        UBOUND=365._dp,UNITS="days")
    !! Parameters of raingreens

    no=no+1; params%RG%shedRate_drySeason => params%memory(no)
             params%RG%shedRate_drySeason  &
                   = set_rdata(1.20e-01_dp,"RG%shedRate_drySeason",LBOUND=1._dp / 180._dp,UBOUND=0.5_dp,UNITS="1/days")
    no=no+1; params%RG%shedRate_aging       => params%memory(no)
             params%RG%shedRate_aging       &
                   = set_rdata(5.0e-03_dp,"RG%shedRate_aging",LBOUND=1._dp / 730._dp,UBOUND=1._dp / 90._dp,UNITS="1/days")
    no=no+1; params%RG%growthRate         => params%memory(no)
             params%RG%growthRate         =  set_rdata(8.0e-02_dp,"RG%growthRate",LBOUND=0.001_dp,UBOUND=0.46_dp,UNITS="1/days")
    no=no+1; params%RG%bucketFill_critical => params%memory(no)
             params%RG%bucketFill_critical &
                   = set_rdata(6.50e-01_dp,"RG%bucketFill_critical",LBOUND=params%all%wilt_point%v,UBOUND=1.0_dp,UNITS="--")
    no=no+1; params%RG%bucketFill_leafout => params%memory(no)
             params%RG%bucketFill_leafout &
                   = set_rdata(4.00e-01_dp,"RG%bucketFill_leafout",LBOUND=params%all%wilt_point%v,UBOUND=1.0_dp,UNITS="--")
    !! Parameters of grasses

    no=no+1; params%GRS%crit_temp  => params%memory(no)
             params%GRS%crit_temp  =  set_rdata(4.0_dp, "GRS%crit_temp",LBOUND=0._dp,UBOUND=20._dp,UNITS="Celsius")
    no=no+1; params%GRS%shedRate_growth  => params%memory(no)
             params%GRS%shedRate_growth   &
                   = set_rdata(1.50e-02_dp,"GRS%shedRate_growth",LBOUND=0._dp,UBOUND=1._dp / 10._dp,UNITS="1/days")
    no=no+1; params%GRS%growthRate  => params%memory(no)
             params%GRS%growthRate  =  set_rdata(9.e-02_dp,"GRS%growthRate",LBOUND=0.05_dp,UBOUND=1.15_dp,UNITS="1/days")
    no=no+1; params%GRS%shedRate_drySeason => params%memory(no)
             params%GRS%shedRate_drySeason &
                   =set_rdata(1.5e-02_dp,"GRS%shedRate_drySeason",LBOUND=1._dp / 180._dp,UBOUND=0.5_dp,UNITS="1/days")

    ! parameters of crops only

    no=no+1; params%CRP%crit_temp  => params%memory(no)
             params%CRP%crit_temp  = set_rdata(10.0_dp,"CRP%crit_temp",LBOUND=0._dp,UBOUND=20._dp,UNITS="Celsius")
    no=no+1; params%CRP%gdd_temp  => params%memory(no)
             params%CRP%gdd_temp  = set_rdata(6.0_dp,"CRP%gdd_temp",LBOUND=0._dp,UBOUND=20._dp,UNITS="Celsius")
    no=no+1; params%CRP%sproud  => params%memory(no)
             params%CRP%sproud  = set_rdata(0.37_dp,"CRP%sproud",LBOUND=params%all%wilt_point%v,UBOUND=1._dp,UNITS="--")
    no=no+1; params%CRP%heat_sum_harvest  => params%memory(no)
             params%CRP%heat_sum_harvest = set_rdata(1300.0_dp,"CRP%heat_sum_harvest", &
             LBOUND=0._dp,UBOUND=3000._dp,UNITS="degree days")
    no=no+1; params%CRP%shedRate_growth  => params%memory(no)
             params%CRP%shedRate_growth  = set_rdata(0.03333_dp,"CRP%shedRate_growth",LBOUND=0._dp,UBOUND=0.04_dp,UNITS="1/days")
    no=no+1; params%CRP%shedRate_rest  => params%memory(no)
             params%CRP%shedRate_rest  =  set_rdata(0.1428_dp,"CRP%shedRate_rest",LBOUND=0._dp,UBOUND=1._dp,UNITS="1/days")
    no=no+1; params%CRP%leafAlloc_frac  => params%memory(no)
             params%CRP%leafAlloc_frac = set_rdata(0.8_dp,"CRP%leafAlloc_frac",LBOUND=0._dp,UBOUND=1._dp,UNITS="--")

    IF(no /= noOfPhenParams) THEN
       CALL finish("init_pheno_params()",&
            "PROGRAMMING ERROR: no (="//trim(int2string(no))//&
            ") different from noOfPhenParams (="//trim(int2string(noOfPhenParams))//")")
    END IF

    ! indicate that parameters have been initialized:
    
    params%initialized = .TRUE.
  END SUBROUTINE init_pheno_params

  ! --- copy_phenology_parameters() ------------------------------------------------------------------------------------------------
  !
  ! Copies phenology parameters paramsIn to paramsOut. If paramsCopy is not initialized, this routine will do so, i.e. it will also
  ! allocate memory.
  !
  SUBROUTINE copy_phenology_parameters(paramsIn,paramsOut)
    TYPE(PhenologyParameters_type),intent(in)    :: paramsIn   
    TYPE(PhenologyParameters_type),intent(inout) :: paramsOut  ! Contains on return the copy

    ! Check initioalization
    
    IF(.NOT. paramsIn%initialized) THEN
       CALL finish("copy_phenology_parameters()","PROGRAMMING ERROR:  paramsIn is not initialized.")
    END IF

    IF(.NOT. paramsOut%initialized) CALL init_pheno_params(paramsOut)

   ! Copy all parameters

    paramsOut%memory(:) = paramsIn%memory(:)
    
  END SUBROUTINE copy_phenology_parameters

  ! --- set_pheno_params_serial() --------------------------------------------------------------------------------------------------
  !
  ! Sets and checks the parameters "PhenParams" of the phenology model by using a serialized parameter set as input.
  !
  SUBROUTINE set_pheno_params_serial(params_serial)
    TYPE(rdata_type), intent(in) :: params_serial(:)

    INTEGER :: i

    !! Check initialization status 

    IF(.NOT. PhenParams%initialized) THEN   ! If not initialized, initialize 
       CALL init_pheno_params(PhenParams)   ! .. first by standard parameters
    END IF

    !! Overwrite with new values

    PhenParams%memory(:) = params_serial(:)

    !! Check bounds and terminate if bounds are not correct

    DO i=1,noOfPhenParams
       CALL check_rdata_bounds(PhenParams%memory(i),TERMINATE=.TRUE.)
    END DO

  END SUBROUTINE set_pheno_params_serial

  ! --- set_phenology_parameters() -------------------------------------------------------------------------------------------------
  !
  ! Sets and checks the parameters "PhenParams" of the phenology model. If the routine is called without arguments, it sets
  ! "PhenParams" standard values.
  !
  SUBROUTINE set_phenology_parameters(params)
    TYPE(PhenologyParameters_type),INTENT(in),optional :: params

    ! --- Set parameters ---

    IF(.NOT. present(params)) THEN
       CALL init_pheno_params(PhenParams)  ! Initialize by standard parameters
    ELSE
       IF(.NOT. associated(params%memory)) THEN
          CALL finish("set_phenology_parameters()","PROGRAMMING ERROR: pointer to physical memory of parameters not associated.")
       END IF
       CALL copy_phenology_parameters(params,PhenParams)
    END IF

    ! --- Check general parameters ---- Note that standard boundary values for the parameters are set in init_pheno_params() -------

    CALL check_rdata_bounds(PhenParams%all%LAI_negligible,TERMINATE=.TRUE.)    ! LAI_negligible
    CALL check_rdata_bounds(PhenParams%all%laiSeed,       TERMINATE=.TRUE.)    ! laiSeed
    CALL check_rdata_bounds(PhenParams%all%wilt_point,    TERMINATE=.TRUE.)    ! wilting point

    ! --- Check parameters common to ever- and summergreen ---

    CALL check_rdata_bounds(PhenParams%EG_SG%tau_pseudo_soil,TERMINATE=.TRUE.) ! tau_pseudo_soil
    CALL check_rdata_bounds(PhenParams%EG_SG%max_chill_days,TERMINATE=.TRUE.)  ! max_chill_days

    ! --- Check parameters of evergreens

    CALL check_rdata_bounds(PhenParams%EG%alternation_temp, TERMINATE=.TRUE.)  ! alternation_temp
    CALL check_rdata_bounds(PhenParams%EG%heat_sum_min,     TERMINATE=.TRUE.)  ! heat_sum_min
    CALL check_rdata_bounds(PhenParams%EG%heat_sum_range,   TERMINATE=.TRUE.)  ! heat_sum_range 
    CALL check_rdata_bounds(PhenParams%EG%chill_decay_const,TERMINATE=.TRUE.)  ! chill_decay_const
    CALL check_rdata_bounds(PhenParams%EG%growthPhaseLength,TERMINATE=.TRUE.)  ! growthPhaseLength
    CALL check_rdata_bounds(PhenParams%EG%shedRate_rest,    TERMINATE=.TRUE.)  ! shedRate_rest
    CALL check_rdata_bounds(PhenParams%EG%growthRate,       TERMINATE=.TRUE.)  ! growthRate

    ! --- Check parameters of summergreens

    CALL check_rdata_bounds(PhenParams%SG%alternation_temp, TERMINATE=.TRUE.)  ! alternation_temp
    CALL check_rdata_bounds(PhenParams%SG%heat_sum_min,     TERMINATE=.TRUE.)  ! heat_sum_min
    CALL check_rdata_bounds(PhenParams%SG%heat_sum_range,   TERMINATE=.TRUE.)  ! heat_sum_range 
    CALL check_rdata_bounds(PhenParams%SG%chill_decay_const,TERMINATE=.TRUE.)  ! chill_decay_const
    CALL check_rdata_bounds(PhenParams%SG%growthPhaseLength,TERMINATE=.TRUE.)  ! growthPhaseLength
    CALL check_rdata_bounds(PhenParams%SG%autumn_event_temp,TERMINATE=.TRUE.)  ! autumn_event_temp
    CALL check_rdata_bounds(PhenParams%SG%shedRate_veget,   TERMINATE=.TRUE.)  ! shedRate_veget
    CALL check_rdata_bounds(PhenParams%SG%shedRate_rest,    TERMINATE=.TRUE.)  ! shedRate_rest
    CALL check_rdata_bounds(PhenParams%SG%growthRate,       TERMINATE=.TRUE.)  ! growthRate

    ! shedRate_rest --- shedRate_veget
    IF(PhenParams%SG%shedRate_veget%v >= PhenParams%SG%shedRate_rest%v) &
!!$         CALL message("#set_phenology_parameters()","PhenParams%SG%shedRate_veget (="//&
!!$                              trim(real2string(real(PhenParams%SG%shedRate_veget%v)))//&
!!$                              " days^(-1)) should be smaller than PhenParams%SG%shedRate_rest (="//&
!!$                              trim(real2string(real(PhenParams%SG%shedRate_rest%v)))//&
!!$                              " days^(-1)).")

    ! --- Check parameters of raingreens

    CALL check_rdata_bounds(PhenParams%RG%shedRate_drySeason,   TERMINATE=.TRUE.)      ! shedRate_drySeason
    CALL check_rdata_bounds(PhenParams%RG%shedRate_aging,       TERMINATE=.TRUE.)      ! minimum shedding rate
    CALL check_rdata_bounds(PhenParams%RG%growthRate,           TERMINATE=.TRUE.)      ! growthRate
    CALL check_rdata_bounds(PhenParams%RG%bucketFill_critical,  TERMINATE=.TRUE.)      ! bucketFill_critical
    CALL check_rdata_bounds(PhenParams%RG%bucketFill_leafout,   TERMINATE=.TRUE.)      ! bucketFill_leafout
    ! --- Check parameters of grasses

    CALL check_rdata_bounds(PhenParams%GRS%crit_temp,           TERMINATE=.TRUE.)      ! crit_temp
    CALL check_rdata_bounds(PhenParams%GRS%shedRate_growth,     TERMINATE=.TRUE.)      ! shedRate_growth
    CALL check_rdata_bounds(PhenParams%GRS%growthRate,          TERMINATE=.TRUE.)      ! growthRate
    CALL check_rdata_bounds(PhenParams%GRS%shedRate_drySeason,  TERMINATE=.TRUE.)      ! shedRate_drySeason

    ! --- Check parameters of crops

    CALL check_rdata_bounds(PhenParams%CRP%crit_temp,           TERMINATE=.TRUE.)      ! crit_temp
    CALL check_rdata_bounds(PhenParams%CRP%gdd_temp,            TERMINATE=.TRUE.)      ! gdd_temp
    CALL check_rdata_bounds(PhenParams%CRP%sproud,              TERMINATE=.TRUE.)      ! sproud
    CALL check_rdata_bounds(PhenParams%CRP%heat_sum_harvest,    TERMINATE=.TRUE.)      ! heat_sum_harvest
    CALL check_rdata_bounds(PhenParams%CRP%shedRate_growth,     TERMINATE=.TRUE.)      ! shedRate_growth
    CALL check_rdata_bounds(PhenParams%CRP%leafAlloc_frac,      TERMINATE=.TRUE.)      ! leafAlloc_frac

  END SUBROUTINE set_phenology_parameters

  ! --- init_status_pheno_params() -------------------------------------------------------------------------------------------------
  !
  ! Returns the initialization status of the current phenology parameters
  !
  FUNCTION init_status_pheno_params()

    LOGICAL :: init_status_pheno_params
    init_status_pheno_params = PhenParams%initialized
  END FUNCTION init_status_pheno_params

  ! --- get_phenology_parameters() -------------------------------------------------------------------------------------------------
  !
  ! Returns the parameters currently used (if initialized; otherwise program is stopped)
  !
  SUBROUTINE get_phenology_parameters(params)
    TYPE(PhenologyParameters_type),INTENT(inout) :: params
    
    IF(PhenParams%initialized) THEN
       CALL copy_phenology_parameters(PhenParams,params)
    ELSE
       CALL finish("get_phenology_parameters()","ERROR: Phenology parameters are not initialized.")
    END IF
  END SUBROUTINE get_phenology_parameters

  ! --- get_pheno_params_serial() --------------------------------------------------------------------------------------------------
  !
  ! Returns the parameters currently used as an array of rdata (if initialized; otherwise program is stopped)
  !
  FUNCTION get_pheno_params_serial()
    TYPE(rdata_type) :: get_pheno_params_serial(1:noOfPhenParams)
    
    IF(PhenParams%initialized) THEN
       get_pheno_params_serial(:) = PhenParams%memory(:)
    ELSE
       CALL finish("get_pheno_params_serial()","ERROR: Phenology parameters are not initialized.")
    END IF
  END FUNCTION get_pheno_params_serial

  ! --- read_phenology_parameters() ------------------------------------------------------------------------------------------------
  !
  ! Reads the parameters for the phenology model. If the interface variable "PhenParamsOut" is present the parameters are returned
  ! in this structure. Otherwise by this routine the currently used parameter values are set to those read in from the file.
  !
  ! The file must conform to the following format:
  ! 
  ! There are two types of lines: 
  ! Comment line: 
  !        A line with "#" in column 1. Such lines are ignored.
  ! Data line:
  !        Data lines contain the data to be read. The general structure of such lines is
  !
  !                     "x keyword = value yyyyy"
  !
  !              x: This is an arbirary character (except '#') in column 1 that is ignored (this may be used for extra information, 
  !                 processed by other routines; usually left blank)
  !        keyword: This is one of the parameter names collected in the global fields "paramNames_???()"
  !              =: This equal sign separates the keyword from the data 
  !          yyyyy: Text of arbitrary length; is not processed
  !
  ! An example file can be generated by print_phenology_parameters().
  !
  SUBROUTINE read_phenology_parameters(paramFileName,PhenParamsOut)
    USE mo_filename, ONLY: find_next_free_unit

    CHARACTER(len=*)              ,INTENT(in)           :: paramFileName !! Name of file where to write the parameters
    TYPE(PhenologyParameters_type),INTENT(out),optional :: PhenParamsOut !! On return this structure contains all the parameters

    CHARACTER(len=1024) :: line,keyword
    LOGICAL             :: fileExists
    INTEGER             :: unitNo,ioerror
    INTEGER             :: posEqual,lineNo
    LOGICAL             :: paramPresent(1:noOfPhenParams)     !! For checking of parameter presence in category PhenParams%all
    CHARACTER(len=32)   :: paramNames(1:noOfPhenParams)       !! List of variable names
    INTEGER             :: i,n
    LOGICAL             :: keywordFound

    TYPE(PhenologyParameters_type) :: PhenParamsTmp

    !! inits

    CALL init_pheno_params(PhenParamsTmp)        !! initialize PhenParamsTmp with standard values. Will be overwritten below.

    paramPresent(:) = .FALSE.

    !! check existence of parameter file

    inquire(FILE=trim(adjustl(paramFileName)),EXIST=fileExists)
    IF(.NOT. fileExists) THEN
       CALL finish("read_phenology_parameters()","Could not find LoGro-P parameter input file "//trim(paramFileName))
    END IF

    !! open parameter file

    unitNo = find_next_free_unit(51,100)                               !! determine free unit number for output
    OPEN(UNIT=unitNo, FILE=paramFileName,ACTION='READ',IOSTAT=ioerror) !! open output file
    IF(ioerror /= 0) THEN 
       CALL finish('read_phenology_parameters()','ERROR when trying to open '//trim(paramFileName))
    END IF

    !! Make list of variable names

    DO n=1,noOfPhenParams
       paramNames(n)  = PhenParamsTmp%memory(n)%n
    END DO

    !! read file and set parameters"

    lineNo=0
    DO
       lineNo = lineNo+1
       READ(unitNo,'(A)',IOSTAT=ioerror) line
       IF(ioerror > 0) THEN
          CALL finish('read_phenology_parameters()','ERROR when reading line '//trim(int2string(lineNo))//&
                                                                               ' of file '//trim(paramFileName))
       ELSE IF(ioerror < 0) THEN
          EXIT !! End of file has been reached
       END IF

       !! Check whether this line is a comment line (starts with a '#')
       IF(line(1:1) == '#') cycle !! If it is a comment line ==> read next line

       !! determine position of "=" because keyword is found before this sign and numerical value after

       posEqual = index(line,'=') 
       IF(posEqual <= 0) THEN
          CALL finish('read_phenology_parameters()','ERROR: "=" sign missing in line '//trim(int2string(lineNo))//&          
                                                                               ' of file '//trim(paramFileName))
       END IF

       !! read keyword
       READ(line(2:posEqual-1),*) keyword !! read keyword with leading blanks removed; column 1 is ignored
       keyword = adjustl(keyword)

       !! read parameters

       keywordFound = .FALSE.
       DO n=1,noOfPhenParams
          IF(trim(keyword) == trim(paramNames(n))) THEN
             READ(line(posEqual+1:),*) PhenParamsTmp%memory(n)%v
             paramPresent(n) = .TRUE.
             keywordFound = .TRUE.
          END IF
       END DO
       
       IF(.NOT. keywordFound) THEN
          CALL finish('read_phenology_parameters()','ERROR: false keyword "'//trim(keyword)//'" in line '//&
                                                     trim(int2string(lineNo))//' of file "'//trim(paramFileName)//&
                                                     '". Note: column 1 is ignored!')
       END IF
    END DO

    !! close file

    CLOSE(unitNo)

    !! check completeness

    DO i=1,noOfPhenParams
       IF(.NOT. paramPresent(i)) THEN
          CALL finish("read_phenology_parameters()",'ERROR: keyword "'//trim(paramNames(i))//&
                                                   '" missing in LoGro-P parameter file "'//trim(paramFileName)//'"')
       END IF
    END DO

    IF(PRESENT(PhenParamsOut)) THEN
       CALL copy_phenology_parameters(PhenParamsTmp,PhenParamsOut)
    ELSE
       CALL set_phenology_parameters(PhenParamsTmp)
    END IF

  END SUBROUTINE read_phenology_parameters

  ! --- check_rdata_bounds()  ----------------------------------------------------------------------------------------------------
  !
  ! Checks whether a variable of type rdata_type lies within the range of it's lower and upper value. 
  !
  ! If optional input variable terminate is true, the program is stopped when "data" is out of its bounds.
  ! 
  ! If variable checkState is present the following return values are possible:
  !       0 ... value is OK, i.e. value is within bounds
  !      -1 ... value is below lower bound
  !      +1 ... value is above upper bound
  !
  ! At least one of the optional variables (TERMINATE,CHECKSTATE) has to be present, otherwise the routine is without effect.
  !
  SUBROUTINE check_rdata_bounds(data,terminate,checkState)
    TYPE(rdata_type),INTENT(in)  :: data
    LOGICAL,optional,INTENT(in)  :: terminate
    INTEGER,optional,INTENT(out) :: checkState

    IF(.NOT. PRESENT(terminate) .AND. .NOT. PRESENT(checkState)) THEN
       CALL finish("check_rdata_bounds()",&
            "PROGRAMMING ERROR: At least one of the optional variables 'TERMINATE'  or 'CHECKSTATE' has to be present")
    END IF

    !! Check lower bound
    IF(data%lBound%exists) THEN
       IF(data%v < data%lBound%value) THEN !! data is smaller than lower bound
          if(PRESENT(terminate)) THEN      !! Finish program if terminate=.true.
             IF(terminate) THEN 
!!$                CALL finish("check_rdata_bounds()",&
!!$                            "ERROR: Parameter "//trim(data%n)//" (="//trim(real2string(real(data%v)))//")"//&
!!$                            " has to be larger than lower bound "//trim(real2string(real(data%lBound%value))))
             END IF
          END IF
          checkState = -1  !! Note that if "terminate" is not present, then "checkState" has to be present
          RETURN
       END IF
    END IF
    !! Check upper bound
    IF(data%uBound%exists) THEN 
       IF(data%v > data%uBound%value) THEN !! data is larger than upper bound
          IF(present(terminate)) THEN      !! Finish program if terminate=.true.
             IF(terminate) THEN 
!!$                CALL finish("check_rdata_bounds()",&
!!$                            "ERROR: Parameter "//trim(data%n)//" (="//trim(real2string(real(data%v)))//")"//&
!!$                            " has to be smaller than upper bound "//trim(real2string(real(data%uBound%value))))
             END IF
          END IF
          checkState = +1 !! Note that if "terminate" is not present, then "checkState" has to be present
          RETURN
       END IF
    END IF
    !! if program comes here: data is within bounds
    IF(PRESENT(checkState)) checkState = 0
  END SUBROUTINE check_rdata_bounds

  ! --- print_phenology_parameters() ----------------------------------------------------------------------------------------------
  !
  ! If optional variable "params_in" is missing this routine prints out the values of all parameters currently used by the
  ! phenology model. But if "params_in" is present it prints out the set of phenology parameters found in this structure.
  ! Whether output is to file or standard out is decided by the presence of the optional variable "paramFileName". 
  ! The file written here can be read by read_phenology_parameters().
  !
  SUBROUTINE print_phenology_parameters(params_in,paramFileName)
    USE mo_filename, ONLY: find_next_free_unit

    TYPE(PhenologyParameters_type),TARGET,optional  :: params_in
    CHARACTER(len=*),           intent(in),optional :: paramFileName !! Name of file where to write the parameters

    INTEGER  ::  unitNo,ioerror
    TYPE(PhenologyParameters_type),POINTER :: params
    CHARACTER(len=8)  ::  date
    CHARACTER(len=10) ::  time
    CHARACTER(len=1)  ::  col1Char !! A single character plotted into column 1, being either '#' or blank ' '
    CHARACTER(len=32) ::  bstring

    !! Check presence of file name

    IF(PRESENT(paramFileName)) THEN !! output to file
       unitNo = find_next_free_unit(51,100)                                !! determine free unit number for output
       OPEN(UNIT=unitNo, FILE=paramFileName,ACTION='WRITE',IOSTAT=ioerror) !! open output file
       IF(ioerror /= 0) THEN 
          CALL finish('print_phenology_parameters()','ERROR when trying to open '//trim(paramFileName))
       END IF
       col1Char = ' ' !! --> when writing to file plot a blank into column 1
    ELSE   !! output to stdout
       unitNo = 6     !! Assume that unit 6 is stdout
       col1Char = '#' !! --> when writing to stdout plot a hash into column 1
    END IF

    !! Check presence of "params" argument in call

    IF(present(params_in)) THEN
       params => params_in   !! output of external phenology parameters
    ELSE
       params => PhenParams  !! output of internal phenology parameters
    END IF
    
    !! Get date and time

    CALL date_and_time(date,time)

    !! Write header

    WRITE(unitNo,'("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PARAMETERS OF LOGRO-PHENOLOGY <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")')
    WRITE(unitNo,'(A)') "# date: "//date(7:8)//"."//date(5:6)//"."//date(1:4)//&
                         " time: "//time(1:2)//":"//time(3:4)//":"//time(5:6)
    !! Write Parameters

    IF(params%initialized) THEN
       WRITE(unitNo,'("# GENERAL (params%all):")')
       bstring=print_rdata_bounds(params%all%LAI_negligible)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%all%LAI_negligible%n,params%all%LAI_negligible%v,&
                                           "              ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%all%laiSeed)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%all%laiSeed%n,params%all%laiSeed%v,&
                                        "              ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%all%wilt_point)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)')   col1Char,params%all%wilt_point%n,params%all%wilt_point%v,&
                                        "              ","   range: ",trim(bstring)

       WRITE(unitNo,'("# EVERGREEN and SUMMERGREEN (params%EG_SG):")')
       bstring=print_rdata_bounds(params%EG_SG%tau_pseudo_soil)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%EG_SG%tau_pseudo_soil%n,params%EG_SG%tau_pseudo_soil%v,&
                                        " days,        ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%EG_SG%max_chill_days)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%EG_SG%max_chill_days%n,params%EG_SG%max_chill_days%v,&
                                        " days,        ","   range: ",trim(bstring)

       WRITE(unitNo,'("# EVERGREEN ONLY (params%EG):")')
       bstring=print_rdata_bounds(params%EG%alternation_temp)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%EG%alternation_temp%n,params%EG%alternation_temp%v,&
                                        " Celsius,     ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%EG%heat_sum_min)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%EG%heat_sum_min%n,params%EG%heat_sum_min%v,&
                                        " degree days, ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%EG%heat_sum_range)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%EG%heat_sum_range%n,params%EG%heat_sum_range%v,&
                                        " degree days, ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%EG%chill_decay_const)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%EG%chill_decay_const%n,params%EG%chill_decay_const%v,&
                                        " days,        ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%EG%growthPhaseLength)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%EG%growthPhaseLength%n,params%EG%growthPhaseLength%v,&
                                        " days,        ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%EG%shedRate_rest)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%EG%shedRate_rest%n,params%EG%shedRate_rest%v,&
                                        " 1/days,      ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%EG%growthRate)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%EG%growthRate%n,params%EG%growthRate%v,&
                                        "              ","   range: ",trim(bstring)

       WRITE(unitNo,'( "# SUMMERGREEN ONLY (params%SG):")')
       bstring=print_rdata_bounds(params%SG%alternation_temp)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%alternation_temp%n,params%SG%alternation_temp%v,&
                                        " Celsius,     ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%heat_sum_min)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%heat_sum_min%n,params%SG%heat_sum_min%v,&
                                        " degree days, ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%heat_sum_range)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%heat_sum_range%n,params%SG%heat_sum_range%v,&
                                        " degree days, ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%chill_decay_const)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%chill_decay_const%n,params%SG%chill_decay_const%v,&
                                        " days,        ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%growthPhaseLength)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%growthPhaseLength%n,params%SG%growthPhaseLength%v,&
                                        " days,        ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%autumn_event_temp)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%autumn_event_temp%n,params%SG%autumn_event_temp%v,&
                                        " Celsius,     ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%shedRate_veget)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%shedRate_veget%n,params%SG%shedRate_veget%v,&
                                        " 1/days,      ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%shedRate_rest)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%shedRate_rest%n,params%SG%shedRate_rest%v,&
                                        " 1/days,      ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%growthRate)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%growthRate%n,params%SG%growthRate%v,&
                                        "              ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%maxLength_gvPhase)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%maxLength_gvPhase%n,params%SG%maxLength_gvPhase%v,&
                                        "              ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%SG%minLength_gvPhase)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%SG%minLength_gvPhase%n,params%SG%minLength_gvPhase%v,&
                                        "              ","   range: ",trim(bstring)

       WRITE(unitNo,'("# RAINGREEN ONLY (params%RG):")')
       bstring=print_rdata_bounds(params%RG%shedRate_drySeason)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%RG%shedRate_drySeason%n,params%RG%shedRate_drySeason%v,&
                                        " 1/days,      ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%RG%shedRate_aging)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%RG%shedRate_aging%n,params%RG%shedRate_aging%v,&
                                        " 1/days,      ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%RG%growthRate)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%RG%growthRate%n,params%RG%growthRate%v,&
                                        "              ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%RG%bucketFill_critical)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%RG%bucketFill_critical%n,params%RG%bucketFill_critical%v,&
                                        "              ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%RG%bucketFill_leafout)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%RG%bucketFill_leafout%n,params%RG%bucketFill_leafout%v,&
                                        "              ","   range: ",trim(bstring)

       WRITE(unitNo,'("# GRASSES ONLY (params%GRS):")')
       bstring=print_rdata_bounds(params%GRS%crit_temp)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%GRS%crit_temp%n,params%GRS%crit_temp%v,&
                                        " Celsius,     ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%GRS%shedRate_growth)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%GRS%shedRate_growth%n,params%GRS%shedRate_growth%v,&
                                        " 1/days,      ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%GRS%growthRate)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%GRS%growthRate%n,params%GRS%growthRate%v,&
                                        "              ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%GRS%shedRate_drySeason)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%GRS%shedRate_drySeason%n,params%GRS%shedRate_drySeason%v,&
                                        " m^2 s/g(C),  ","   range: ",trim(bstring)

       WRITE(unitNo,'("# CROPS ONLY (params%CRP):")')
       bstring=print_rdata_bounds(params%CRP%crit_temp)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%CRP%crit_temp%n,params%CRP%crit_temp%v,&
                                        " Celsius,     ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%CRP%gdd_temp)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%CRP%gdd_temp%n,params%CRP%gdd_temp%v,&
                                        " Celsius,     ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%CRP%sproud)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%CRP%sproud%n,params%CRP%sproud%v,&
                                        " Celsius,     ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%CRP%heat_sum_harvest)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%CRP%heat_sum_harvest%n,params%CRP%heat_sum_harvest%v,&
                                        " degree days, ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%CRP%shedRate_growth)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%CRP%shedRate_growth%n,params%CRP%shedRate_growth%v,&
                                        " 1/days,      ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%CRP%shedRate_rest)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%CRP%shedRate_rest%n,params%CRP%shedRate_rest%v,&
                                        "              ","   range: ",trim(bstring)
       bstring=print_rdata_bounds(params%CRP%leafAlloc_frac)
       WRITE(unitNo,'(A1,6X,A25,"= ",ES11.4,A,A,A)') col1Char,params%CRP%leafAlloc_frac%n,params%CRP%leafAlloc_frac%v,&
                                        "              ","   range: ",trim(bstring)
    ELSE
       CALL finish("print_phenology_parameters()","Parameters not initialized!! Use function init_pheno_params() first!")
    END IF
    WRITE(unitNo,'("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")')

    !! Eventually close file
    IF(PRESENT(paramFileName)) CLOSE(unitNo)
  END SUBROUTINE print_phenology_parameters

  ! --- print_rdata_bounds() -------------------------------------------------------------------------------------------------------
  !
  ! Returns the boundary values of a variable of type rdata as string "[lbound:ubound]", where lbound and ubound are the lower and
  ! upper boundary values (if they exist).
  !
  CHARACTER(len=32) PURE FUNCTION print_rdata_bounds(data)
    TYPE(rdata_type),INTENT(in) :: data

    IF(data%lBound%exists) THEN
       WRITE(print_rdata_bounds,'("[",ES11.4," :")') data%lBound%value
    ELSE
       WRITE(print_rdata_bounds,'("[   -inf     :")')
    END IF
    IF(data%uBound%exists) THEN
       WRITE(print_rdata_bounds,'(A,ES11.4," ]")') trim(print_rdata_bounds),data%uBound%value
    ELSE
       WRITE(print_rdata_bounds,'(A,"  +inf      ]")') trim(print_rdata_bounds)
    END IF
  END FUNCTION print_rdata_bounds

END MODULE mo_phenology
