!
MODULE mo_cosp_rttov_types
  ! Description:
  ! defines all derived types for RTTOV
  !
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0   01/12/2002  New F90 code with structures (P Brunel A Smith)
  !  1.1   29/01/2003  Add CO2 variable gaz to profile sturcture (P Brunel)
  !                    Add rain and solid precip. to profile cloud structure
  !  1.2   13/05/2003  Add structure for transmissions and optical depths (F Chevallier)
  !  1.3      08/2003  Add scattering facility (F Chevallier)
  !  1.4   18/09/2003  Add kice and kradip to profile_cloud_type (P Francis)
  !  1.5   09/12/2003  Change type for mclayer to INTEGER (R Saunders)
  !  1.6   06/01/2004  Add CO2 to ref profile (R Saunders)
  !  1.7   02/06/2004  Add fast model version compatibility level in coef type (P. Brunel)
  !  1.8   17/05/2005  Add q to profile_cloud_type ( U O'Keeffe)
  !  1.9   12/01/2006  Marco Matricardi (ECMWF):
  !           --       Added raytracing_type structure  for computation of variable
  !           --       local zenith angle.
  !           --       Added N2O,CO and CH4 to profile_type structure.
  !           --       Modified type rttov_coef.
  ! 1.10   30/01/2007  Removed frequency tagged variables R Saunders
  ! 1.11   15/03/2007  Added single stream transmittances for RTTOV
  !                    o/p R Saunders
  ! 1.12   11/10/2007  Move iaernum & iaertyp profile members to profile_aux
  !                    P.Marguinaud
  ! 1.13   05/07/2007  Added extra radiance outputs from cloud (R Saunders)
  ! 1.14   22/11/2007  RTTOV-SCATT version 9 - remove cld_radiance_type; slim
  !                    down profile_cloud_type ( A Geer )
  ! 1.15   12/12/2007  add IncTop (R Saunders)
  ! 1.16   12/08/2008  add IncZeeman, Be, cosbk (P. Rayer)
  ! 1.16   03/10/2008  RTTOV-SCATT revised cloud partitioning (A Geer)
  ! 1.17   01/09/2009  Added salinity for FASTEM-4 (Mark Liu)
  ! 1.19   03/11/2009  Transmittances / optical depths on levels (A Geer)
  ! 1.20   02/12/2009  Added principal component capability (Marco Matricard)
  ! 1.21   02/12/2009  Introduced new variables for the mixed cloud scheme (Marco Matricardi)
  ! 1.22   17/06/2010  Introduced spacetop flag to zero opdeps at user's
  !                    model-top in Zeeman channels (P Rayer)
  ! 1.23   05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
  ! 1.24   14/10/2010  Remove rt8_mode (J Hocking)
  ! 1.25   09/11/2012  Add Lambertian mode (R Saunders)
  ! 1.26   10/01/2013  Add PMC shifts (P Rayer)
  !
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
  ! 2011/11 Introduced PC multiple band and cloudy computations, Marco Matricardi, ECMWF
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:

  ! Imported Parameters:
  use mo_cosp_rttov_const, only: &
        & fastem_sp,     &
        & ncldtyp,       &
        & interp_rochon

!  Use parkind1, Only : jpim, jprb, jplm
  USE mo_kind,  only:wp
  Implicit None

  ! Surface skin
  Type sskin_type
  Integer*4 :: surftype        ! 0=land, 1=sea, 2=sea-ice
  Integer*4 :: watertype       ! 0=fresh water, 1=ocean water
  Real(wp)    :: t               ! radiative skin temperature (K)
  Real(wp)    :: salinity        ! practical salinity unit %o
  Real(wp)    :: fastem(fastem_sp)  ! land/sea-ice surface parameters for fastem-2
  End Type sskin_type

  ! Surface 2m
  Type s2m_type
  Real(wp) :: t                  ! temperature (K)
  Real(wp) :: q                  ! water vapour (ppmv)
  Real(wp) :: o                  ! ozone (ppmv)
  Real(wp) :: p                  ! surface pressure (hPa)
  Real(wp) :: u                  ! U 10m wind component (m/s)
  Real(wp) :: v                  ! V 10m wind component (m/s)
  Real(wp) :: wfetc              ! Wind fetch (metres)
  End Type s2m_type


  ! structure for atmospheric profiles on model pressure levels
  Type profile_type
  Character(Len=128) :: id
  Integer*4 :: date(3)  ! Year, Month, Day
  Integer*4 :: time(3)  ! Hour, Minute, Second
     ! number of atmospheric levels/layers
  Integer*4 :: nlevels
  Integer*4 :: nlayers
     ! ozone, CO2 and cloud liquid water profiles available
     ! atmosphere defined on nlevels
  Real(wp), Pointer :: p(:)      ! pressure (hPa)
  Real(wp), Pointer :: t(:)      ! temperature (K)
  Real(wp), Pointer :: q(:)      ! water vapour (ppmv)
  Real(wp), Pointer :: o3(:)     ! ozone (ppmv)
  Real(wp), Pointer :: co2(:)    ! carbon dioxide (ppmv)
  Real(wp), Pointer :: n2o(:)    ! n2o(ppmv)
  Real(wp), Pointer :: co(:)     ! co(ppmv)
  Real(wp), Pointer :: ch4(:)    ! ch4(ppmv)
  Real(wp), Pointer :: clw(:)    ! cloud liquid water (kg/kg)
  Real(wp), Pointer :: aerosols(:,:)
  Real(wp), Pointer :: cloud(:,:)
  Real(wp), Pointer :: cfrac(:)
  Real(wp), Pointer :: icede(:)  ! ice particle effective diameter (microns)
  Integer*4       :: idg
  Integer*4       :: ish
     ! surface
     Type(sskin_type) :: skin
     Type(s2m_type)   :: s2m
     !angles
  Real(wp) :: zenangle
  Real(wp) :: azangle
  Real(wp) :: sunzenangle
  Real(wp) :: sunazangle
  Real(wp) :: elevation
  Real(wp) :: latitude
  Real(wp) :: longitude
  Real(wp) :: snow_frac     ! snow coverage fraction for IR emissivity atlas (0 - 1)
  Real(wp) :: soil_moisture ! soil moisture (m^3/m^3)
  ! Be - Earth magnetic field strength (Gauss)
  ! cosbk - cosine of the angle between the Earth magnetic field and wave
  !         propagation direction
  Real(wp) :: Be
  Real(wp) :: cosbk
     ! Black body cloud
  Real(wp) :: ctp               ! cloud top pressure  (hPa)
  Real(wp) :: cfraction         ! cloud fraction (0 - 1) 1 for 100% cloud cover
     !CNam:
    logical :: ozone_data = .false.
    logical :: co2_data   = .false.
    logical :: n2o_data   = .false.
    logical :: co_data    = .false.
    logical :: ch4_data   = .false.

    logical :: addrefrac  = .false.
    logical :: clw_data   = .false.
    logical :: aer_data   = .false.
    logical :: cld_data   = .false.

    logical :: addsolar   = .false.

  End Type profile_type

  ! structure for atmospheric profile additional input for RTTOV-SCATT,
  ! with information on clouds and rain for each level. Full level pressures,
  ! t, q etc. should be placed in the profile_type.

  Type profile_cloud_type

  Integer*4 :: nlevels ! number of atmospheric levels (nlevels+1 for ph)
  logical :: use_totalice ! False => separate ice and snow  True => total ice
  logical :: mmr_snowrain ! snow and rain input units are: False => kg/m2/s  True => kg/kg
  Real(wp) :: cfrac      ! Average cloud fraction (only used if lusercfrac = .true.)

  Real(wp), Pointer :: ph(:)        ! nlevels+1 of half-level model pressures (hPa)
  Real(wp), Pointer :: cc(:)        ! nlevels of cloud cover
  Real(wp), Pointer :: clw(:)       ! nlevels of cloud liquid water (kg/kg)
  Real(wp), Pointer :: ciw(:)       ! nlevels of cloud ice water (kg/kg)
  Real(wp), Pointer :: totalice(:)  ! nlevels of total ice (kg/kg)
  Real(wp), Pointer :: rain(:)      ! nlevels of rain (units: see mmr_snowrain)
  Real(wp), Pointer :: sp(:)        ! nlevels of solid precipitation (units: see mmr_snowrain)

  End Type profile_cloud_type

  ! satellite geometry
  TYPE geometry_type
    Real(wp) :: sinzen
    Real(wp) :: sinzen_sq
    Real(wp) :: coszen
    Real(wp) :: coszen_sq
    Real(wp) :: seczen
    Real(wp) :: seczen_sq
    Real(wp) :: seczen_sqrt
    Real(wp) :: seczen_minus1
    Real(wp) :: seczen_minus1_sq
    Real(wp) :: sinview
    Real(wp) :: sinview_sq
    Real(wp) :: cosview_sq
    Real(wp) :: normzen
    Real(wp) :: viewang

    Real(wp) :: sinzen_sun
    Real(wp) :: sinlat
    Real(wp) :: coslat
  END TYPE geometry_type
  
  ! The actual predictor arrays
  Type rttov_path_pred
    Real(wp), Pointer     :: mixedgas(:,:,:)          ! (nmixed,  nlevels, nprofiles )
    Real(wp), Pointer     :: watervapour(:,:,:)       ! (nwater,  nlevels, nprofiles)
    Real(wp), Pointer     :: ozone(:,:,:)             ! (nozone,  nlevels)
    Real(wp), Pointer     :: wvcont(:,:,:)            ! (nwvcont, nlevels)
    Real(wp), Pointer     :: co2(:,:,:)               ! (nco2,    nlevels)
    Real(wp), Pointer     :: n2o(:,:,:)               ! (nn2o,    nlevels)
    Real(wp), Pointer     :: co(:,:,:)                ! (nco,     nlevels)
    Real(wp), Pointer     :: ch4(:,:,:)               ! (nch4,    nlevels)
    Real(wp), Pointer     :: clw(:,:)                 ! (         nlevels)  
    Real(wp), Pointer     :: pmc(:,:,:,:)    ! pressure modulated cell (npmc,nlevels,nprofiles.nchannels)  
  End Type rttov_path_pred
  
  ! Predictors
  Type predictors_type
     ! the nxxxx could be set to 0 to indicate the abscence
     ! of the predictor, in that case there is no need to
     ! allocate the corresponding predictor
  Integer*4 :: nlevels   ! number of levels for predictors (all same)
  Integer*4 :: nmixed    ! number of variables for Mixed Gases
  Integer*4 :: nwater    ! number of variables for Water Vapour
  Integer*4 :: nozone    ! number of variables for Ozone
  Integer*4 :: nwvcont   ! number of variables for WV Continuum
  Integer*4 :: nco2      ! number of variables for CO2
  Integer*4 :: nn2o      ! number of variables for N2O
  Integer*4 :: nco       ! number of variables for CO
  Integer*4 :: nch4      ! number of variables for CH4
  Integer*4 :: ncloud    ! number of variables for MW Cloud
  Integer*4 :: npmc      ! number of variables for pressure modulated cell correction
  
  Type(rttov_path_pred) :: path1  ! Predictors for surface-satellite path (always required)
  Type(rttov_path_pred) :: path2  ! Predictors for sun-surface-satellite path (only required for solar)
  End Type predictors_type

  ! Fast coefficients - separate structure to allow for non-PW and PW if necessary
  Type rttov_fast_coef
    ! separate arrays to allow different number of variables for each gaz
    Real(wp), Pointer      :: mixedgas(:,:,:)     ! Mixed gases coefs  (levels, channels, variables)
    Real(wp), Pointer      :: watervapour(:,:,:)  ! Water vapour coefs (levels, channels, variables)
    Real(wp), Pointer      :: ozone(:,:,:)        ! Ozone coefs        (levels, channels, variables)
    Real(wp), Pointer      :: wvcont(:,:,:)       ! WV Cont coefs      (levels, channels, variables)
    Real(wp), Pointer      :: co2(:,:,:)          ! CO2 coefs          (levels, channels, variables)
    Real(wp), Pointer      :: n2o(:,:,:)          ! N2O coefs          (levels, channels, variables)
    Real(wp), Pointer      :: co(:,:,:)           ! CO coefs           (levels, channels, variables)
    Real(wp), Pointer      :: ch4(:,:,:)          ! CH4 coefs          (levels, channels, variables)
  End Type rttov_fast_coef

  TYPE rttov_nlte_coef
    REAL(wp), POINTER :: coef(:,:,:,:) ! ncoef x nsat x nsol x nchan
    REAL(wp), POINTER :: sol_zen_angle(:), sat_zen_angle(:)
    REAL(wp), POINTER :: cos_sol(:), sec_sat(:)    
    INTEGER*4 :: ncoef, nsol, nsat, nchan
    INTEGER*4 :: start_chan, end_chan
    REAL(wp)    :: max_sat_angle
  END TYPE rttov_nlte_coef

  Type rttov_coef
     ! Structure for the storage of RTTOV coefficients
     ! this may differ from what is stored in the coefficient files especially
     ! for the units (ie kg/kg to ppmv)
     ! Gases are separated in MxG WV O3
     ! Number of levels is the same for all gases (taken from MxG).
     !
  Integer*4 :: id_platform  ! platform   (see documentation or MOD_CPARAM)
  Integer*4 :: id_sat    ! satellite  (.....)
  Integer*4 :: id_inst    ! instrument (.....)
  Integer*4 :: id_sensor  ! sensor
     !  1 = Infrared
     !  2 = Micro Wave
     !  3 = High resolution
  Integer*4 :: id_comp_lvl  ! RTTOV coefficient file version number
  Integer*4 :: id_comp_pc   ! Principal component coefficient file version number
  Integer*4 ,Dimension(3) :: id_creation_date  ! YYYY MM DD
  Character (len=80)    :: id_creation    ! Creation comment
  Character (len=32)    :: id_Common_name  ! usual name of the satellite
  Character (len=132)   :: line_by_line(100)
  Character (len=132)   :: readme_srf(100) ! readme Spectral Response Function


     !FAST_MODEL_VARIABLES section
  Character (len=32)    :: fmv_model_def  ! FMV definition (RTTOV6 OPTRAN RTTOV7)
  Integer*4               :: fmv_model_ver  ! fast model version compatibility level
  Integer*4               :: fmv_ori_nchn   ! number of channels in original file
  Integer*4               :: fmv_chn        ! number of channels read from file into coef structure
  Integer*4               :: fmv_gas        ! number of gases in file
  Integer*4, pointer      :: fmv_gas_id(:)    ! gas id. number i gas_id list (fmv_gas)
  Integer*4, Pointer      :: fmv_gas_pos(:)   ! respective position of each gas of gas_id list (ngases_max)
  Integer*4, Pointer      :: fmv_var(:)       ! number of variables/predictors by gaz (fmv_gas)
  Integer*4, Pointer      :: fmv_coe(:)       ! number of coefficients by gaz (fmv_gas)
  Integer*4, Pointer      :: fmv_int(:)       ! number of spectral intervals by gaz (fmv_gas)
  Integer*4, Pointer      :: fmv_lvl(:)       ! number of levels(pres/absorber) by gaz (fmv_gas)

  Integer*4               :: nmixed         ! number of variables/predictors for Mixed Gases
  Integer*4               :: nwater         ! number of variables/predictors for Water Vapour
  Integer*4               :: nozone         ! number of variables/predictors for Ozone
  Integer*4               :: nwvcont        ! number of variables/predictors for WV continuum
  Integer*4               :: nco2           ! number of variables/predictors for CO2
  Integer*4               :: nn2o           ! number of variables/predictors for N2O
  Integer*4               :: nco            ! number of variables/predictors for CO
  Integer*4               :: nch4           ! number of variables/predictors for CH4
  Integer*4               :: nlevels        ! number of levels(pres/absorber) same for all gases
  Integer*4               :: nlayers        ! number of layers(pres/absorber) nlevels-1
  logical               :: IncZeeman      ! Flag to include Zeeman effect for this sensor
  Integer*4               :: ncmixed        ! number of coefficients for Mixed Gases
  Integer*4               :: ncwater        ! number of coefficients for Water Vapour
  Integer*4               :: ncozone        ! number of coefficients for Ozone
  Integer*4               :: ncwvcont       ! number of coefficients for WV continuum
  Integer*4               :: ncco2          ! number of coefficients for CO2
  Integer*4               :: ncn2o          ! number of coefficients for N2O
  Integer*4               :: ncco           ! number of coefficients for CO
  Integer*4               :: ncch4          ! number of coefficients for CH4
  
  ! JAH - to remove?
  Integer*4               :: nintmixed
  Integer*4               :: nintwater
  Integer*4               :: nintozone
  Integer*4               :: nintwvcont
  Integer*4               :: nintco2
  Integer*4               :: nintn2o
  Integer*4               :: nintco
  Integer*4               :: nintch4

     !GAZ_UNITS section
     ! gases are in the order of gas id codes
  Integer*4, Pointer      :: gaz_units(:)  ! unit of gaz concentration for each gaz
                                             ! default value is specific conc. (kg/kg)
                                             ! value inside RTTOV calculations (ppmv)
     !FILTER_FUNCTIONS section  array size is fmv_chn
  logical :: ff_val_bc    ! are any band corrections to be applied?
  logical :: ff_val_gam   ! any gamma corrections?
  Integer*4 ,Pointer :: ff_ori_chn(:)   ! original chan number
  Integer*4 ,Pointer :: ff_val_chn(:)   ! validity of the channel (1=OK)
  Real(wp) ,Pointer :: ff_cwn (:)         ! cental wave number (cm-1)
  Real(wp) ,Pointer :: ff_bco (:)         ! band correction offset (K)
  Real(wp) ,Pointer :: ff_bcs (:)         ! band correction slope (K/K)
  Real(wp) ,Pointer :: ff_gam (:)         ! gamma factor transm. correction

     !TRANSMITTANCE_TRESHOLD section  array size is fmv_chn
  Integer*4 ,Pointer :: tt_chn(:)
  Integer*4 ,Pointer :: tt_val_chn(:)
  Real(wp)    ,Pointer :: tt_cwn (:)
  Real(wp)    ,Pointer :: tt_a0(:)
  Real(wp)    ,Pointer :: tt_a1(:)

     !PLANCK_WEIGHTED section array size if fmv_chn
  Integer*4 ,Pointer :: pw_chn(:)
  Integer*4 ,Pointer :: pw_val_chn(:)   ! 0 => non-PW thermal coefs, 1 => PW thermal coefs

     !SOLAR_SPECTRUM section array size is fmv_chn
  Integer*4 ,Pointer :: ss_chn(:)
  Integer*4 ,Pointer :: ss_val_chn(:)
  Real(wp)    ,Pointer :: ss_cwn (:)
  Real(wp)    ,Pointer :: ss_solar_spectrum(:)

     !WATER_OPTICAL_CONSTANT section array size is fmv_chn
  Integer*4 ,Pointer :: woc_chn(:)
  Real(wp)    ,Pointer :: woc_cwn (:)
  Complex(Kind=wp) ,Pointer :: woc_waopc_ow(:)
  Complex(Kind=wp) ,Pointer :: woc_waopc_fw(:)

     !WAVE_SPECTRUM section array size is ws_nomega
     !Data used to compute the frequency spectrum of the JONSWAP
     !wave model surface wave.
  Integer*4          :: ws_nomega
  Real(wp)  ,Pointer   :: ws_npoint(:)
  Real(wp)  ,Pointer   :: ws_k_omega(:)

     !FUNDAMENTAL_CONSTANTS section
  Real(wp) :: fc_speedl         ! speed of light (cm/s)
  Real(wp) :: fc_planck_c1      ! first radiation constant (mW/(m2*sr*cm-4))
  Real(wp) :: fc_planck_c2      ! second radiation constant (cm*K)
  Real(wp) :: fc_sat_height     ! satellite nominal altitude (km)

     !FASTEM section
  Integer*4 :: fastem_ver      ! fastem version number
  Integer*4, Pointer :: fastem_polar(:)  ! polarisation of each channel
     ! 0 = 0.5 V+H
     ! 1 = 90 - incident angle
     ! 2 = incident angle
     ! 3 = vertical
     ! 4 = horizontal
     ! 5 = V+H
     ! Full stokes vector

     !SSIREM section     array size is fmv_chn
     ! ems =   ssirem_a0
     !       - ssirem_a1*(zen**ssirem_xzn1)
     !       - ssirem_a2*(zen**ssirem_xzn2)
     ! where zen is satellite zenith angle in degrees, divided by 60.
  Integer*4 :: ssirem_ver                ! version number
  Integer*4,  Pointer  :: ssirem_chn(:)   ! original chan number
  Real(wp),  Pointer     :: ssirem_a0(:)    ! constant coef
  Real(wp),  Pointer     :: ssirem_a1(:)    ! first coef
  Real(wp),  Pointer     :: ssirem_a2(:)    ! second coef
  Real(wp),  Pointer     :: ssirem_xzn1(:)  ! 1st exponent on zenith angle
  Real(wp),  Pointer     :: ssirem_xzn2(:)  ! 2nd exponent on zenith angle

     !REFERENCE_PROFILE section  defined on Mixed gases pressure levels
     ! Not working for OPTRAN gas absorber levels
     ! gases are in the order of gas id codes
     ! unit for mr in coeff file is kg/kg or ppmv (see gaz_units section)
     ! unit for mr for optical depth calculations is ppmv
  Real(wp), Pointer      :: ref_prfl_p(:)     ! pressure  (hPa)       (levels)
  Real(wp), Pointer      :: ref_prfl_t(:,:)   ! temperature (K)       (levels, gases)
  Real(wp), Pointer      :: ref_prfl_mr(:,:)  ! mixing ratio (ppmv)   (levels, gases)
     !PROFILE_LIMITS section
     ! gases are in the order of gas id codes
     ! unit for mr in coeff file is kg/kg or ppmv (see gaz_units section)
     ! unit for mr for optical depth calculations is ppmv
  Real(wp), Pointer      :: lim_prfl_p(:)       ! pressure  (hPa)       (levels)
  Real(wp), Pointer      :: lim_prfl_tmax(:)    ! max temperature (K)   (levels)
  Real(wp), Pointer      :: lim_prfl_tmin(:)    ! min temperature (K)   (levels)
  Real(wp), Pointer      :: lim_prfl_gmax(:,:)  ! max mixing r (ppmv) (levels, gases)
  Real(wp), Pointer      :: lim_prfl_gmin(:,:)  ! min mixing r (ppmv) (levels, gases)


     !FAST_COEFFICIENTS/SOLAR_FAST_COEFFICIENTS section
     ! For non-PW instruments, "solar" will point to "thermal" coefs
     ! For instruments with solar-affected PW channels, both thermal and solar 
     ! structures will be populated from coef file
  Type(rttov_fast_coef), Pointer :: thermal             ! FAST_COEFFICIENTS
  Type(rttov_fast_coef), Pointer :: solar               ! SOLAR_FAST_COEFFICIENTS when present
  logical             :: solarcoef           ! .TRUE. if solar fast coefs present in file: need
                                                        ! to know whether solar points to thermal or is
                                                        ! allocated separately.
  logical             :: nltecoef = .FALSE.  ! .TRUE. if nlte coefs present
  TYPE(rttov_nlte_coef), POINTER :: nlte_coef           ! nlte_coef
  
       !PRESSURE_MODULATED_CELL section
  logical       :: pmc_shift = .FALSE.  ! .TRUE. if pmc shift coefs present
  Real(wp)          :: pmc_lengthcell       ! cell length (cm)
  Real(wp), Pointer :: pmc_pnominal(:)      ! nominal cell pressure (hPa) - nchannels
  Real(wp)          :: pmc_tempcell          ! cell temperature (K)
  Real(wp)          :: pmc_betaplus1        ! co2 band-average: self-HW/air-HW
  Integer*4       :: pmc_nlay             ! number of layers used
  Integer*4       :: pmc_nvar             ! number of variables used
  Real(wp), Pointer :: pmc_coef(:,:,:)     ! pressure moodulated cell corrections - nlevels, nchannels, nvariables
  Real(wp), Pointer :: pmc_ppmc(:)          ! actual cell pressure (hPa) - nchannels

   ! JAH - to remove?
  Real(wp), Pointer      :: mixedgasint(:,:)
  Real(wp), Pointer      :: watervapourint(:,:)
  Real(wp), Pointer      :: ozoneint(:,:)
  Real(wp), Pointer      :: wvcontint(:,:)
  Real(wp), Pointer      :: co2int(:,:)
  Real(wp), Pointer      :: n2oint(:,:)
  Real(wp), Pointer      :: coint(:,:)
  Real(wp), Pointer      :: ch4int(:,:)
      
     ! Auxillary variables
  Real(wp)               :: ratoe       ! ratio (H+R)/R  H=sat height, R=Earth radius
  Integer*4            :: mwcldtop    ! Upper layer for MW LWP calcs
  Real(wp), pointer      :: planck1(:)        ! C1 * Nu**3
  Real(wp), pointer      :: planck2(:)        ! C2 * Nu
  Real(wp), pointer      :: frequency_ghz(:)  ! frequency in GHz

     ! other predictor variables see Science and Validation report
  Real(wp), pointer      :: dp(:)        ! interval between standard p levels (hPa)
  Real(wp), pointer      :: dpp(:)       ! pressure based variable (hPa**2)
  Real(wp), pointer      :: tstar(:)     ! layer temp (K)
  Real(wp), pointer      :: to3star(:)   ! layer temp for O3 calculations (K)
  Real(wp), pointer      :: wstar(:)     ! layer WV  (ppmv)
  Real(wp), pointer      :: ostar(:)     ! layer O3  (ppmv)
  Real(wp), pointer      :: co2star(:)   ! layer co2 (ppmv)
  Real(wp), pointer      :: n2ostar(:)   ! layer n2o (ppmv)
  Real(wp), pointer      :: costar(:)    ! layer co  (ppmv)
  Real(wp), pointer      :: ch4star(:)   ! layer ch4 (ppmv)
  End Type rttov_coef

  Type rttov_scatt_coef
     ! Structure for the storage of RTTOV_SCATT coefficients
  Integer*4 :: nhydro ! Number of hydrometeors in computation
  Integer*4 :: mtype  ! Number of hydrometeors     in Mie tables
  Integer*4 :: mfreqm ! Number of frequencies      in Mie tables
  Integer*4 :: mtemp  ! Number of temperature bins in Mie tables
  Integer*4 :: mwc    ! Number of water bins       in Mie tables
  Real(wp)    :: offset_temp_rain       ! temperature offset in table for rain type
  Real(wp)    :: offset_temp_sp         ! temperature offset in table for solid prec. type
  Real(wp)    :: offset_temp_liq        ! temperature offset in table for cloud water type
  Real(wp)    :: offset_temp_ice        ! temperature offset in table for cloud ice type
  Real(wp)    :: offset_temp_totalice   ! temperature offset in table for total ice type
  Real(wp)    :: offset_water           ! liquid/ice water offset in table
  Real(wp)    :: scale_water            ! log10(liquid/ice water) scaling factor in table
  Real(wp)    :: from_scale_water       ! 10**(1._wp/scale_water)
  Real(wp)    :: conv_rain(2)           ! coefficients for rain unit conversion (mm.h-1 to g.m-3)
  Real(wp)    :: conv_sp  (2)           ! coefficients for solid prec. unit conversion (mm.h-1 to g.m-3)
  Real(wp)    :: conv_liq (2)           ! coefficients for cloud water conversion (not used)
  Real(wp)    :: conv_ice (2)           ! coefficients for cloud ice conversion   (not used)
  Real(wp)    :: conv_totalice (2)      ! coefficients for total ice conversion   (not used)
  Real(wp), pointer :: mie_freq(:)      ! list of frequencies in Mie table
  Real(wp), pointer :: ext(:,:,:,:)     ! extinction coefficent table
  Real(wp), pointer :: ssa(:,:,:,:)     ! single scattering albedo table
  Real(wp), pointer :: asp(:,:,:,:)     ! assymetry parameter table

  End Type rttov_scatt_coef


  TYPE profile_aux_s
    Integer*4 :: nearestlev_surf ! nearest model level above surface
    Real(wp)    :: pfraction_surf  ! pressure fraction of surface in model layer (hPa)
    Integer*4 :: nearestlev_ctp  ! nearest model level above cloud top
    Real(wp) :: pfraction_ctp   ! pressure fraction of cloud top pressure in layer (hPa)
    Real(wp) :: cfraction       ! cloud fraction (0 - 1) 1 for 100% cloud cover
  END TYPE profile_aux_s

  ! Auxillary profile variables
  ! variables calculated by the model from profile
  TYPE profile_aux
    logical :: on_coef_levels
    TYPE(profile_aux_s), POINTER :: s(:)
    Real(wp), POINTER :: debye_prof(:,:,:)  ! Debye terms
    Real(wp), POINTER :: relhum(:,:)        !Relative humidity
    Real(wp), POINTER :: relhumref(:,:)
    Real(wp), POINTER :: dg(:,:)            !Generalized effective diameter
    Real(wp), POINTER :: fac1_dg(:,:)       !Intermediate variables used to compute the
    Real(wp), POINTER :: fac2_dg(:,:)       !generalized diameter
    Real(wp), POINTER :: fac3_dg(:,:)
    Integer*4, POINTER :: iaertyp(:,:,:)
    Integer*4, POINTER :: iaernum(:,:)
    
    Real(wp), POINTER :: t_layer(:,:)       ! avg layer temperature
    Real(wp), POINTER :: w_layer(:,:)       ! avg layer humidity
    Real(wp), POINTER :: o3_layer(:,:)      ! avg layer humidity
    Real(wp), POINTER :: dt(:,:)            ! deviation from ref prof
    Real(wp), POINTER :: dto(:,:)           ! deviation from ref prof
    Real(wp), POINTER :: tr(:,:), tr_r(:,:)  ! ratio t / ref_t
    Real(wp), POINTER :: wr(:,:), wr_sqrt(:,:), wr_rsqrt(:,:)           ! 
    Real(wp), POINTER :: or(:,:), or_sqrt(:,:)
    Real(wp), POINTER :: tw(:,:), tw_sqrt(:,:), tw_4rt(:,:)            ! ratio t / ref_t
    Real(wp), POINTER :: ww(:,:), ww_r(:,:)            ! 
    Real(wp), POINTER :: ow(:,:), ow_r(:,:), ow_sqrt(:,:), ow_rsqrt(:,:)           ! 
    Real(wp), POINTER :: SUM(:,:)
    
  END TYPE profile_aux

  ! Auxillary profile variables for RTTOV_SCATT
  ! variables calculated by the model from profile
  Type profile_scatt_aux
  Real(wp), pointer :: cfrac(:)        ! horizontal cloud fraction (one value used for all layers)
  Real(wp), pointer :: ems_bnd(:)      ! surface emissivity for boundary conditions
  Real(wp), pointer :: ref_bnd(:)      ! surface emissivity for boundary conditions
  Real(wp), pointer :: ems_cld(:)      ! surface emissivity taking into account cloud/rain impact on od
  Real(wp), pointer :: ref_cld(:)      ! surface reflectivity taking into account cloud/rain impact on od
  Real(wp), pointer :: dz(:,:)         ! layer depth   [km]
  Real(wp), pointer :: tbd(:,:)        ! temperature at layer boundary [K]
  Real(wp), Pointer :: clw(:,:)        ! cloud liquid water (g/m3)
  Real(wp), Pointer :: ciw(:,:)        ! cloud ice water (g/m3)
  Real(wp), Pointer :: totalice(:,:)   ! total ice (g/m3)
  Real(wp), Pointer :: rain(:,:)       ! rain (g/m3)
  Real(wp), Pointer :: sp(:,:)         ! solid precipitation (g/m3)
!RWS  Real(wp), pointer :: mclayer(:)  ! upper level cloud layer
  Integer*4, pointer :: mclayer(:)   ! upper level cloud layer
  Real(wp), pointer :: delta(:,:)      ! (= ext*dz/coszen)
  Real(wp), pointer :: tau(:,:)        ! optical depths (= exp(-delta))
  Real(wp), pointer :: ext(:,:)        ! extinction coefficient integreated over hydrometeor types
  Real(wp), pointer :: ssa(:,:)        ! single scattering albedo integreated over hydrometeor types
  Real(wp), pointer :: asm(:,:)        ! asymetry parameter integreated over hydrometeor types [-1,1]
  Real(wp), pointer :: lambda(:,:)     ! eddington approx. variable
                                  ! (= sqrt( 3*ext*ext*(1-ssa)*(1-ssa*asm) )
  Real(wp), pointer :: h (:,:)         ! boundary condition variable (= 1.5_wp*ext(1-ssa*asm))
  Real(wp), pointer :: b0(:,:)         ! lower level temperature
  Real(wp), pointer :: b1(:,:)         ! temperature gradient
  Real(wp), pointer :: bn(:,:)         ! upper level temperature
  end type profile_scatt_aux

  type opdp_path_type
    ! path optical depths as predicted or interpolated (unitless)
    Real(wp), pointer :: atm_level(:,:) ! neg optical depth for thermal radiation (levels to space),
                                               ! size (levels, channels)
    Real(wp), pointer :: sun_level_path2(:,:) ! neg optical depth for solar radiation (levels to space) for
                                                     ! combined sun-surface-satellite path, size (levels, channels)
  end type opdp_path_type

  type transmission_type
    Real(wp), pointer  :: tau_total(:)    ! surface-satellite transmittance (channels)
    Real(wp), pointer  :: tau_levels(:,:) ! level-satellite transmittance (levels,channels)
    Real(wp), pointer  :: tausun_total_path2(:)    ! sun-surface-satellite solar transmittance
    Real(wp), pointer  :: tausun_levels_path2(:,:) ! sun-level-satellite solar transmittance for each level
    Real(wp), pointer  :: tausun_total_path1(:)    ! surface-satellite solar transmittance
    Real(wp), pointer  :: tausun_levels_path1(:,:) ! level-satellite solar transmittance for each level  
  end type

  type rttov_path_transmission
    ! Transmissions and optical depths (unitless)
    Real(wp), pointer  :: tau_surf_p(:,:)         ! transmittance from surface (streams,channels)
    Real(wp), pointer  :: tau_surf_p_r(:,:)       ! reciprocal transmittance from surface (streams,channels)
    Real(wp), pointer  :: tau_surf(:,:)           ! transmittance from surface (streams,channels)
    Real(wp), pointer  :: tau_surf_r(:,:)         ! reciprocal transmittance from surface (streams,channels)
    Real(wp), pointer  :: tau_level(:,:,:)        ! transmittance from each standard pressure level
                                                         ! (levels,streams,channels)
    Real(wp), pointer  :: tau_level_r(:,:,:)      ! reciprocal transmittance from each standard pressure level
                                                         ! (levels,streams,channels)
    Real(wp), pointer  :: tau_level_p(:,:,:)      ! transmittance from each standard pressure level
                                                         ! (levels,streams,channels)
    Real(wp), pointer  :: tau_level_p_r(:,:,:)    ! reciprocal transmittance from each standard pressure level
                                                         ! (levels,streams,channels)
    Real(wp), pointer  :: od_singlelayer(:,:,:)   ! single-layer optical depth
    Real(wp), pointer  :: od_singlelayer_r(:,:,:) ! reciprocal single-layer optical depth
    Real(wp), pointer  :: od_sfrac(:,:)
    Real(wp), pointer  :: od_sfrac_r(:,:)
    Real(wp), pointer  :: od_frac_ac(:,:)
    Real(wp), pointer  :: tau_surf_ac(:,:)

    Real(wp), pointer  :: fac2(:,:,:)    ! Mask for integration calculation: thermal and solar path1
  end type
  
  type transmission_type_aux
    Real(wp), pointer  :: fac1(:,:,:)    ! Mask for integration calculation
    Real(wp), pointer  :: surf_fac(:,:)
    Real(wp) :: anynegtau ! used to store information about presence of any negative transmittances

    Type(rttov_path_transmission), POINTER :: thermal_path1
    Type(rttov_path_transmission), POINTER :: solar_path2
    Type(rttov_path_transmission), POINTER :: solar_path1
  end type transmission_type_aux

  type radiance_type
     ! Primary radiance and corresponding brightness temperature/reflectance outputs
     ! Array size is (nchannels) or (nlayers,nchannels)
     ! unit for radiance is mw/cm-1/ster/sq.m
     ! unit for temperature is Kelvin
     ! reflectances are bi-directional reflectance factors (BRF, unitless)
     !
    Real(wp), pointer  :: clear(:)       ! clear sky radiance
    Real(wp), pointer  :: total(:)       ! cloudy radiance for given cloud
    Real(wp), pointer  :: bt_clear(:)    ! Brightness temp equivalent to clear radiance
    Real(wp), pointer  :: bt(:)          ! Brightness temp equivalent to total radiance
    Real(wp), pointer  :: refl_clear(:)  ! Reflectance calculated from clear radiance
    Real(wp), pointer  :: refl(:)        ! Reflectance calculated from total radiance
    Real(wp), pointer  :: overcast(:,:)  ! overcast radiance for opaque cloud at level bounding
                                                !   bottom of each layer
    Real(wp), pointer  :: cloudy(:)      ! 100% cloudy radiance for given cloud (simple cloud scheme)
                                                !   or same as total (addclouds/addaerosl true)
  end type radiance_type

  type radiance2_type
    ! Secondary radiances optionally calculated in direct model only for clear-sky with no solar contribution
    ! Array size is (nchannels) or (nlayers,nchannels)
    !
    Real(wp), pointer  :: upclear(:)     ! clear sky upwelling radiance without reflection term
    Real(wp), pointer  :: dnclear(:)     ! clear sky downwelling radiance
    Real(wp), pointer  :: refldnclear(:) ! reflected clear sky downwelling radiance
    Real(wp), pointer  :: up(:,:)        ! sum( B * dT ) above cloud upwelling radiance from each layer
    Real(wp), pointer  :: down(:,:)      ! sum ( B / T**2 dT ) above cloud downwelling radiance from
                                                !   each layer
    Real(wp), pointer  :: surf(:,:)      ! radiance at surface emitted from a black cloud
  end type radiance2_type

  type radiance_aux
     ! auxillary calculation arrays for RTE integration
     ! Direct model arrays need to be passed to TL AD and K codes
     ! array size is of (nchannels) or (nlevels, nchannels)
  Real(wp), pointer :: air(:,:)
  Real(wp), pointer :: surfair(:)
  Real(wp), pointer :: skin(:)
  Real(wp), pointer :: cosmic(:)
  Real(wp), pointer :: air_t_eff(:,:)
  Real(wp), pointer :: surf_t_eff(:)
  Real(wp), pointer :: skin_t_eff(:)
  Real(wp), pointer :: cosmic_t_eff(:)

  Real(wp), pointer :: up(:,:,:)               ! sum( B * dT )
  Real(wp), pointer :: down(:,:,:)             ! sum ( B / T**2 dT )
  Real(wp), pointer :: up_solar(:,:,:)         ! sum( B * dT )
  Real(wp), pointer :: down_solar(:,:,:)       ! sum ( B / T**2 dT )
  Real(wp), pointer :: meanrad_up(:,:)
  Real(wp), pointer :: meanrad_down(:,:)
  Real(wp), pointer :: meanrad_up_solar(:,:)
  Real(wp), pointer :: meanrad_down_solar(:,:)
  Real(wp), pointer :: down_ref(:,:,:)
  Real(wp), pointer :: down_ref_solar(:,:,:)
  Real(wp), pointer :: FAC1_2(:,:,:)
  Real(wp), pointer :: FAC2_2(:,:,:)
  Real(wp), pointer :: FAC3_2(:,:,:)
  Real(wp), pointer :: FAC4_2(:,:,:)
  Real(wp), pointer :: FAC5_2(:,:,:)
  Real(wp), pointer :: FAC6_2(:,:,:)
  Real(wp), pointer :: FAC7_2(:,:,:)
  Real(wp), pointer :: FAC1_3(:,:)
  Real(wp), pointer :: FAC2_3(:,:)
  Real(wp), pointer :: FAC3_3(:,:)
  Real(wp), pointer :: FAC4_3(:,:)
  Real(wp), pointer :: FAC5_3(:,:)
  Real(wp), pointer :: FAC6_3(:,:)
  Real(wp), pointer :: FAC7_3(:,:)
  Real(wp), pointer :: cloudy(:,:)
  end type radiance_aux

  type raytracing_type
  Real(wp), pointer :: LTICK   (:,:)   !(levels,profiles)
  Real(wp), pointer :: HGPL    (:,:)   !(levels)
  Real(wp), pointer :: DAIR    (:,:)   !(levels,profiles)
  Real(wp), pointer :: DMAIR   (:,:)   !(levels,profiles)
  Real(wp), pointer :: DMAIR_R (:,:)   !(levels,profiles)
  Real(wp), pointer :: REFRACTIVITY(:,:)   !(levels,profiles)
  Real(wp), pointer :: R       (:,:)   !(levels,profiles)
  Real(wp), pointer :: R_R     (:,:)   !(levels,profiles)
  Real(wp), pointer :: Z_R     (:,:)   !(levels,profiles)
  Real(wp), pointer :: RATOESUN(:,:)   !(levels,profiles)
  Real(wp), pointer :: RATOESAT(:,:)   !(levels,profiles)
  Real(wp), pointer :: ZASUN   (:,:)   !(levels,profiles)
  Real(wp), pointer :: ZASAT   (:,:)   !(levels,profiles)
  Real(wp), pointer :: INT     (:,:)   !(levels,profiles)
  Real(wp), POINTER :: ZTEMP   (:,:)   !(levels,profiles)
  Real(wp), pointer :: HL      (:,:)   !(levels,profiles)
  Real(wp), pointer :: PPW     (:,:)   !(levels,profiles)
  Real(wp), pointer :: DISPCO2 (:,:)   !(levels,profiles)
  Real(wp), pointer :: PATHSAT (:,:)   !(levels,profiles)
  Real(wp), pointer :: PATHSAT_rsqrt (:,:)   !(levels,profiles)
  Real(wp), pointer :: PATHSAT_sqrt (:,:)   !(levels,profiles)
  Real(wp), pointer :: PATHSUN (:,:)   !(levels,profiles)
  Real(wp), pointer :: PATHEFF (:,:)
  Real(wp), pointer :: CO2_CM  (:)     !(profiles)
  end type raytracing_type

  type sunglint_type_s
  Real(wp)          :: CSI
  Real(wp)          :: ALFA
  Real(wp)          :: C_SHAD
  Real(wp)          :: P_PRIME
  Real(wp)          :: PXY_GAMMAXY
  Real(wp)          :: GAMMA_O
  Real(wp)          :: GAMMA_P
  Real(wp)          :: G_SHAD
  Real(wp)          :: GAMMAX
  Real(wp)          :: Q_SHAD
  Real(wp)          :: ZENSAT
  Real(wp)          :: ZENSUN
  Real(wp)          :: DAZNG
  Real(wp)          :: FAC1
  Real(wp)          :: A_SHAD
  Real(wp)          :: B_SHAD
  Real(wp)          :: LAMBDA_A
  Real(wp)          :: LAMBDA_B
  Real(wp)          :: X_U
  Real(wp)          :: ALFA1
  Real(wp)          :: OMEGA_M
  Real(wp)          :: WINDSP
  Real(wp)          :: WANGL
  Real(wp)          :: GAMMA_SQ
  Real(wp)          :: GLINT
  Real(wp)          :: OMEGA
  end type sunglint_type_s

  type sunglint_type
  type(sunglint_type_s), pointer :: s(:)
  Real(wp),pointer  :: BETA (:,:)
  Real(wp),pointer  :: PSI  (:,:)
  end type sunglint_type

  type transmission_scatt_ir_type
  Real(wp)   , pointer  :: opdps         (:,:)
  Real(wp)   , pointer  :: opdpa         (:,:)
  Real(wp)   , pointer  :: gpar          (:,:)
  Real(wp)   , pointer  :: gpartot       (:,:)
  Real(wp)   , pointer  :: opdpscls      (:,:,:)
  Real(wp)   , pointer  :: opdpacls      (:,:,:)
  Real(wp)   , pointer  :: gparcls       (:,:,:)
  Real(wp)   , pointer  :: opdpaerla     (:,:)
  Real(wp)   , pointer  :: opdpcldla     (:,:)
  Real(wp)   , pointer  :: opdpsaer      (:,:)
  Real(wp)   , pointer  :: opdpaaer      (:,:)
  Real(wp)   , pointer  :: gparaera      (:,:)
  Real(wp)   , pointer  :: gparaer       (:,:)
  Real(wp)   , pointer  :: azphup        (:,:)
  Real(wp)   , pointer  :: azphdo        (:,:)
  Real(wp)   , pointer  :: azphupcls     (:,:,:)
  Real(wp)   , pointer  :: azphdocls     (:,:,:)
  Real(wp)   , pointer  :: azphuptot     (:,:)
  Real(wp)   , pointer  :: azphdotot     (:,:)
  Real(wp)   , pointer  :: azphaerup     (:,:)
  Real(wp)   , pointer  :: azphaerdo     (:,:)
  Real(wp)   , pointer  :: azphaerupa    (:,:)
  Real(wp)   , pointer  :: azphaerdoa    (:,:)
  Real(wp)   , pointer  :: phasintupref  (:,:,:)
  Real(wp)   , pointer  :: phasintdoref  (:,:,:)
  Real(wp)   , pointer  :: opdpabs       (:,:,:)
  Real(wp)   , pointer  :: opdpsca       (:,:,:)
  Real(wp)   , pointer  :: opdpac        (:,:,:)
  Real(wp)   , pointer  :: opdpacl       (:,:,:)
  Real(wp)   , pointer  :: opdpacsun     (:,:,:)
  Real(wp)   , pointer  :: opdpaclsun    (:,:,:)
  Real(wp)   , pointer  :: azphacup      (:,:,:)
  Real(wp)   , pointer  :: azphacdo      (:,:,:)
  Real(wp)   , pointer  :: opdpext       (:,:,:)
  Real(wp)   , pointer  :: ssa           (:,:,:)
  end type transmission_scatt_ir_type

  type rttov_opt_param
    Real(wp), pointer :: abs(:,:)     ! absorption coef (nchannels,nlayers)
    Real(wp), pointer :: sca(:,:)     ! scattering coef (nchannels,nlayers)
    Real(wp), pointer :: bpr(:,:)     ! b parameter (nchannels,nlayers)
    Real(wp), pointer :: pha(:,:,:)   ! phase function (nchannels,nlayers,nphangle)
    Real(wp), pointer :: phangle(:)   ! angles over which phase fns defined (nphangle)

    ! The following are for RTTOV internal purposes
    Real(wp)             :: minphadiff     ! minimum difference between phase angles
    Real(wp),    pointer :: cosphangle(:)  ! cosine of phase angles (nphangle)
    Integer*4, pointer :: iphangle(:)    ! array indexes (size depends on minphadiff)
  end type rttov_opt_param

  type rttov_pccomp
    Real(wp), pointer  :: pcscores(:)    ! Principal component scores
    Real(wp), pointer  :: bt_pccomp(:)   ! Brightness temp equivalent to radiances
                                                ! reconstructed using principal components
    Real(wp), pointer  :: total_pccomp(:)! Radiances reconstructed using principal
                                                ! components
  end type rttov_pccomp

  type rttov_coef_pccomp1
    Integer*4           :: fmv_pc_npred          ! Number of predictors in the regression set
    Integer*4, pointer  :: predictindex  (:)     ! Predictors channel indices
    Real(wp)   , pointer  :: coefficients  (:,:)   ! Regression coefficients
  end type rttov_coef_pccomp1

  type rttov_coef_pccomp2
    Real(wp)   , pointer  :: eigenvectors  (:,:)   ! Eigenvectors
  end type rttov_coef_pccomp2

  type rttov_coef_pccomp
    Integer*4           :: fmv_pc_comp_pc
    Integer*4           :: fmv_pc_cld
    Integer*4           :: fmv_pc_msets          ! Maximum number of regression sets
    Integer*4           :: fmv_pc_bands          ! Number of bands
    Integer*4           :: fmv_pc_mnum           ! Maximum number of eigenvectors
    Integer*4           :: fmv_pc_mchn           ! Maximum number of channels
    Integer*4           :: fmv_pc_nchn           ! Number of channels
    Integer*4           :: fmv_pc_nchn_noise     ! Number of channels for which instrument noise is available
    Integer*4           :: fmv_pc_nche           ! Number of channels for which emissisity coefs are available
    Integer*4           :: fmv_pc_gas            ! Number of gases for which a reference profile is given
    Integer*4, pointer  :: fmv_pc_sets   (:)     ! Number of regression sets in each band
    Integer*4, pointer  :: emiss_chn     (:)     ! Number of channels for which emissivity coefficients are
    Real(wp), pointer  :: emiss_c1      (:)     ! Emissivity coefficient
    Real(wp), pointer  :: emiss_c2      (:)     ! Emissivity coefficient
    Real(wp), pointer  :: emiss_c3      (:)     ! Emissivity coefficient
    Real(wp), pointer  :: emiss_c4      (:)     ! Emissivity coefficient
    Real(wp), pointer  :: emiss_c5      (:)     ! Emissivity coefficient
    Real(wp), pointer  :: emiss_c6      (:)     ! Emissivity coefficient
    Real(wp), pointer  :: emiss_c7      (:)     ! Emissivity coefficient
    Real(wp), pointer  :: emiss_c8      (:)     ! Emissivity coefficient
    Real(wp), pointer  :: emiss_c9      (:)     ! Emissivity coefficient
    Integer*4           :: fmv_pc_nlev           ! Number of reference profile levels
    Real(wp), Pointer     :: ref_pc_prfl_p (:)     ! pressure  (hPa)       (levels)
    Real(wp), Pointer     :: ref_pc_prfl_mr(:,:)   ! mixing ratio (ppmv)   (levels)
    Real(wp), Pointer     :: lim_pc_prfl_tmin(:)   ! Profile limit :temperature
    Real(wp), Pointer     :: lim_pc_prfl_tmax(:)   ! Profile limit :temperature
    Real(wp), Pointer     :: lim_pc_prfl_qmin(:)   ! Profile limit :water vapour
    Real(wp), Pointer     :: lim_pc_prfl_qmax(:)   ! Profile limit :water vapour
    Real(wp), Pointer     :: lim_pc_prfl_ozmin(:)  ! Profile limit :ozone
    Real(wp), Pointer     :: lim_pc_prfl_ozmax(:)  ! Profile limit :ozone
    Real(wp)              :: lim_pc_prfl_pmin      ! Surface pressure
    Real(wp)              :: lim_pc_prfl_pmax      ! Surface pressure
    Real(wp)              :: lim_pc_prfl_tsmin     ! Surface temperature
    Real(wp)              :: lim_pc_prfl_tsmax     ! Surface temperature
    Real(wp)              :: lim_pc_prfl_skmin     ! Skin temperature
    Real(wp)              :: lim_pc_prfl_skmax     ! Skin temperature
    Real(wp)              :: lim_pc_prfl_wsmin     ! 10m wind speed
    Real(wp)              :: lim_pc_prfl_wsmax     ! 10m wind speed
    Real(wp), Pointer     :: co2_pc_ref    (:)     ! Fixed co2 profile to be used in the computation of PC's
    Real(wp), Pointer     :: n2o_pc_ref    (:)     ! Fixed n2o profile to be used in the computation of PC's
    Real(wp), Pointer     :: co_pc_ref     (:)     ! Fixed co  profile to be used in the computation of PC's
    Real(wp), Pointer     :: ch4_pc_ref    (:)     ! Fixed ch4 profile to be used in the computation of PC's
    Real(wp), Pointer     :: noise_in      (:)     ! Noise values for the channels whose radiances are
                                                          ! reconstrucetd using principal components
    Real(wp), Pointer     :: noise         (:)     ! Noise values for the channels whose radiances are
                                                          ! used as predictors in the computation of principal components
    Integer*4, Pointer  :: ff_ori_chn_in(:)
    Real(wp),    Pointer  :: ff_cwn_in(:)          ! central wave number of reconstructed radiances
    Real(wp),    Pointer  :: ff_bco_in (:)         ! band correction offset (K)
    Real(wp),    Pointer  :: ff_bcs_in (:)         ! band correction slope (K/K)
    Real(wp),    pointer  :: planck1_in(:)         ! C1 * Nu**3
    Real(wp),    pointer  :: planck2_in(:)         ! C2 * Nu
    type(rttov_coef_pccomp1), pointer:: pcreg     (:,:)
    type(rttov_coef_pccomp2), pointer:: eigen     (:)
  end type rttov_coef_pccomp

  type rttov_coef_scatt_ir
  Integer*4           :: fmv_aer_chn    ! number of channels for which aerosol optical parameters are stored
  Integer*4           :: fmv_wcl_chn
  Integer*4           :: fmv_icl_chn
  Integer*4           :: fmv_aer_pha_chn
  Integer*4           :: fmv_wcl_pha_chn
  Integer*4           :: fmv_icl_pha_chn
  Integer*4           :: fmv_aer_comp
  Integer*4           :: fmv_wcl_comp
  Integer*4           :: fmv_icl_comp
  Integer*4           :: fmv_icl_ishp
  Integer*4           :: fmv_aer_pha_ioff
  Integer*4           :: fmv_wcl_pha_ioff
  Integer*4           :: fmv_icl_pha_ioff
  Integer*4           :: fmv_aer_ph
  Integer*4           :: fmv_wcl_ph
  Integer*4           :: fmv_icl_ph
  Integer*4           :: icl_nabs
  Integer*4           :: icl_nsca
  Integer*4           :: icl_nbpr
  Character(Len=4),   Pointer  :: fmv_aer_comp_name(:)
  Character(Len=4),   Pointer  :: fmv_wcl_comp_name(:)
  Character(Len=16),  Pointer  :: fmv_icl_comp_name(:,:)
  Integer*4, Pointer  :: fmv_aer_rh    (:)
  Integer*4, Pointer  :: fmv_wcl_rh    (:)
  Real(wp), Pointer  :: fmv_aer_rh_val(:)
  Real(wp), Pointer  :: fmv_wcl_rh_val(:)
  Integer*4, Pointer  :: ifmv_aer_ph_val(:)
  Integer*4, Pointer  :: ifmv_wcl_ph_val(:)
  Integer*4, Pointer  :: ifmv_icl_ph_val(:)
  Real(wp), Pointer  :: fmv_aer_ph_val(:)
  Real(wp), Pointer  :: fmv_wcl_ph_val(:)
  Real(wp), Pointer  :: fmv_icl_ph_val(:)
  Real(wp), Pointer  :: fmv_aer_ph_val_cos(:)
  Real(wp), Pointer  :: fmv_wcl_ph_val_cos(:)
  Real(wp), Pointer  :: fmv_icl_ph_val_cos(:)
  Real(wp)           :: fmv_aer_ph_val_min
  Real(wp)           :: fmv_wcl_ph_val_min
  Real(wp)           :: fmv_icl_ph_val_min
  Real(wp), Pointer  :: fmv_icl_dg    (:,:)
  Integer*4, Pointer  :: aer_pha_chanlist(:)
  Integer*4, Pointer  :: wcl_pha_chanlist(:)
  Integer*4, Pointer  :: icl_pha_chanlist(:)
  Integer*4, Pointer  :: aer_pha_index(:)
  Integer*4, Pointer  :: wcl_pha_index(:)
  Integer*4, Pointer  :: icl_pha_index(:)
  Real(wp)   , pointer  :: abs           (:,:)
  Real(wp)   , pointer  :: sca           (:,:)
  Real(wp)   , pointer  :: bpr           (:,:)
  Real(wp)   , pointer  :: pha           (:,:,:)
  Real(wp)   , pointer  :: confac        (:)
  end type rttov_coef_scatt_ir
  
  type rttov_coef_optpiclb
    ! interpolation factors for frequency
    Integer*4, pointer  :: iwn           (:) ! channel
    Integer*4, pointer  :: jwn           (:)
    Real(wp)   , pointer  :: dx_dwn        (:)
  end type rttov_coef_optpiclb
  
  type rttov_optpar_ir
    type(rttov_coef_scatt_ir), pointer :: optpaer(:)
    type(rttov_coef_scatt_ir), pointer :: optpwcl(:)
    type(rttov_coef_scatt_ir), pointer :: optpicl(:)
    type(rttov_coef_optpiclb), pointer :: optpiclb
  end type rttov_optpar_ir

  type ircld_type
  Integer*4, pointer  :: nstream(:)
  Integer*4, pointer  :: nstreamref(:)
  Integer*4, pointer  :: iloop(:)
  Integer*4, pointer  :: icount(:)
  Integer*4, pointer  :: icounstr(:)
  Integer*4, pointer  :: icount1(:)
  Real(wp)   , pointer  :: xstrclr(:)
  Integer*4, pointer  :: icldarr   (:,:,:)
  Real(wp)   , pointer  :: xstrref1  (:,:,:)
  Real(wp)   , pointer  :: xstrref2  (:,:,:)
  Integer*4, pointer  :: cldtyp    (:,:,:)
  Integer*4, pointer  :: indexstr  (:,:)
  Integer*4, pointer  :: icount1ref(:,:)
  Integer*4, pointer  :: iloopin   (:,:)
  Integer*4, pointer  :: iflag     (:,:)
  Real(wp)   , pointer  :: xstr      (:,:)
  Real(wp)   , pointer  :: xstrminref(:,:)
  Real(wp)   , pointer  :: xstrref   (:,:)
  Real(wp)   , pointer  :: cldcfr    (:,:)
  Real(wp)   , pointer  :: maxcov    (:,:)
  Real(wp)   , pointer  :: xstrmax   (:,:)
  Real(wp)   , pointer  :: xstrmin   (:,:)
  Real(wp)   , pointer  :: a         (:,:)
  Real(wp)   , pointer  :: ntotref   (:,:)

  Real(wp)   , pointer  :: tave      (:,:)
  Real(wp)   , pointer  :: wmixave   (:,:)
  Real(wp)   , pointer  :: xpresave  (:,:)
  Real(wp)   , pointer  :: ppv       (:,:)
  Real(wp)   , pointer  :: esw       (:,:)
  Real(wp)   , pointer  :: esi       (:,:)
  logical, pointer  :: flag      (:,:)
  end type ircld_type


  Type rttov_coefs
    logical         :: initialised = .false.
    Type(rttov_coef)           :: coef
    Type(rttov_coef_scatt_ir)  :: coef_scatt_ir
    Type(rttov_optpar_ir)      :: optp
    type(rttov_coef_pccomp)    :: coef_pccomp
  End Type


  ! Configuration options that make sense across all flavours of RTTOV
  Type rttov_opts_config
    logical :: apply_reg_limits = .false.
    logical :: verbose          = .true.
    logical :: do_checkinput    = .true.
  End Type

  ! PC options
  Type rttov_opts_pc
    logical :: addpc     = .false.
    Integer*4 :: ipcbnd    = -1
    Integer*4 :: ipcreg    = -1
    logical :: addradrec = .false.
  End Type

  ! General radiative transfer options
  Type rttov_opts_rt_all
    logical :: addrefrac  = .false.
    logical :: switchrad  = .false.
    logical :: use_q2m    = .true.
  End Type

  ! VIS/IR-only radiative transfer options
  Type rttov_opts_rt_ir
    Type(rttov_opts_pc) :: pc

    logical :: addsolar           = .false.
    logical :: do_nlte_correction = .false.
    logical :: addaerosl          = .false.
    logical :: addclouds          = .false.
    logical :: user_aer_opt_param = .false.
    logical :: user_cld_opt_param = .false.
    Real(wp)    :: cldstr_threshold   = -1.0_wp ! Recommend set this negative unless only running direct

    logical :: ozone_data = .false.
    logical :: co2_data   = .false.
    logical :: n2o_data   = .false.
    logical :: co_data    = .false.
    logical :: ch4_data   = .false.
  End Type

  ! MW-only radiative transfer options
  Type rttov_opts_rt_mw
    Integer*4 :: fastem_version = 5 ! Valid range: 1-5, otherwise version taken from coef file
    logical :: clw_data       = .false.
    logical :: do_lambertian  = .false.  ! Flag for setting MW Lambertian reflection
  End Type

  ! Options for internal vertical interpolation and vertical grid setup
  Type rttov_opts_interp
    logical :: addinterp   = .false.
    Integer*4 :: interp_mode = interp_rochon
    logical :: lgradp      = .false.
    logical :: spacetop    = .true.
  End Type

  ! RTTOV options
  Type rttov_options
    Type(rttov_opts_config)  :: config
    Type(rttov_opts_rt_all)  :: rt_all
    Type(rttov_opts_rt_ir)   :: rt_ir
    Type(rttov_opts_rt_mw)   :: rt_mw
    Type(rttov_opts_interp)  :: interpolation
  End Type

  ! RTTOV-SCATT options
  ! Note that RTTOV-SCATT deliberately does not give user control over core R/T options.
  Type rttov_options_scatt
    Type(rttov_opts_config) :: config
    logical      :: lusercfrac = .false.
    logical      :: use_q2m    = .true.
    Integer*4      :: fastem_version = 5
  End Type


  Type rttov_traj
!
! Hold RTTOV trajectory; these variables have counterparts in TL, AD, K,
! and their dimensions are known before running RTTOV (nlevels, nprofiles, nchannels)
! it is possible to allocate these variables from outside RTTOV
!
    Type(profile_Type), Pointer :: profiles_COEF(:)
    Type(predictors_Type)       :: predictors
    Type(raytracing_type)       :: raytracing
    Type(raytracing_type)       :: raytracing_COEF
    Type(ircld_type)            :: ircld
    Type(opdp_path_type)        :: opdp_path
    Type(opdp_path_type)        :: opdp_path_COEF

    Real(wp), Pointer :: thermrefl(:)    ! Surface refl for thermal calcs (nchanprof)
    Real(wp), Pointer :: fresnrefl(:)    ! Fresnel reflection coefficients (nchanprof)
    Type(profile_aux)  :: aux_prof
    Type(profile_aux)  :: aux_prof_COEF
    Type(transmission_scatt_ir_type)  :: transmission_scatt_ir
    Type(sunglint_type):: sunglint

    Type(rttov_coefs),         Pointer :: coefs
    Integer*4   :: nchannels
    Integer*4   :: nlevels
    Integer*4   :: nlayers
    Type(rttov_options)  :: opts
  End Type

  Type rttov_traj_dyn
!
! Hold RTTOV trajectory; these variables have counterparts in TL, AD, K,
! but their dimensions are known when RTTOV starts running (nstreams)
!
    Integer*4               :: nstreams
    TYPE(radiance_aux              ) :: auxrad_stream
    TYPE(transmission_scatt_ir_type) :: transmission_scatt_ir_stream
    TYPE(transmission_type_aux     ) :: transmission_aux
  End Type

  Type rttov_path_traj_trans

    ! Structure to hold optical depth and transmittance data
    ! within static trajectory.

    Real(wp),     POINTER :: tau_ref       (:,:)
    Real(wp),     POINTER :: tau_ref_surf  (:)
    Real(wp),     POINTER :: tau_level     (:,:)! sat to level transmittance
    Real(wp),     POINTER :: tau_surf      (:)
    Real(wp),     POINTER :: od_level      (:,:)! sat to level optical depth
    Real(wp),     POINTER :: opdp_ref_COEF (:,:)! layer optical depth before threshold

    Real(wp),     POINTER :: od_singlelayer(:,:)! single layer optical depth
    Real(wp),     POINTER :: od_frac       (:)
  End Type
 
  Type rttov_traj_sta
!
! Hold RTTOV trajectory; these variables do not have counterparts in TL, AD, K
!
    logical,  POINTER         :: thermal(:)        ! switch for thermal calculations (nchanprof)
    logical,  POINTER         :: solar(:)          ! switch for solar calculations (nchanprof)
    logical                   :: dothermal         ! flag to indicate thermal calculations required
    logical                   :: dosolar           ! flag to indicate solar calculations required

    logical, POINTER          :: do_lambertian(:)
    Real(wp), POINTER             :: solar_spec_esd(:) ! Solar spectrum adjusted for esd (nchanprof)
    Real(wp), POINTER             :: refl_norm(:)      ! Normalisation factor for solar surface reflectance
    TYPE(rttov_path_traj_trans), POINTER :: thermal_path1
    TYPE(rttov_path_traj_trans), POINTER :: solar_path2
    TYPE(rttov_path_traj_trans), POINTER :: solar_path1
    TYPE(geometry_Type), POINTER         :: angles(:)         ! geometry angles
    TYPE(geometry_Type), POINTER         :: angles_COEF(:)    ! geometry angles
    TYPE(profile_Type),  POINTER         :: profiles_COEF_ref(:)
    TYPE(radiance_aux)                   :: auxrad
    TYPE(rttov_chanprof), POINTER        :: chanprof_in(:)
    TYPE(rttov_chanprof), POINTER        :: chanprof_pc(:)
  End Type

!
  Type rttov_chanprof
    Integer*4 :: chan
    Integer*4 :: prof
  End Type

  Type rttov_emissivity
    Real(wp)    :: emis_in
    Real(wp)    :: emis_out
  End Type
  
  Type rttov_reflectance
    Real(wp)    :: refl_in
    Real(wp)    :: refl_out
    Real(wp)    :: refl_cloud_top
  End Type
  
  Type rttov_lbl_check
    Real(wp), Pointer :: atm_layer(:,:)
    logical :: plane_geometry
  End Type


END MODULE mo_cosp_rttov_types
