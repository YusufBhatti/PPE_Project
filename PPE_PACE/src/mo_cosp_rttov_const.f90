!
MODULE mo_cosp_rttov_const
  ! Description:
  ! Definition of all parameters (constants) for RTTOV
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
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0   01/12/2002  New F90 code with structures (P Brunel A Smith)
  !  1.1   29/01/2003  New platforms and instruments (P Brunel)
  !                    Hard limits for input profiles
  !  1.2   19/02/2003  Some changes to limits and comments (R Saunders)
  !  1.3   06/05/2003  Change version number to 7.3.1
  !                    and add references for physical constants (P Brunel)
  !  1.4      08/2003  Added variables for MW scattering (F Chevallier)
  !  1.5   18/09/2003  Added coefficients for cloud absorption properties (P Francis)
  !  1.6   15/10/2003  Added new sections in parameter files for scatt   (F Chevallier)
  !  1.7   23/11/2003  Added new definitions of polarisations 2.1 (S English)
  !  1.8   25/08/2005  Made inst_name a parameter (R Saunders)
  !  1.9   11/01/2006  Added logical flag for surface humidity use (R Saunders)
  !  1.10  12/01/2006  Marco Matricardi (ECMWF):
  !           --       Added variables for CO2,CO,N2O and CH4 molecules.
  !           --       Added parameters for the computation of the refractive index
  !           --       of air.
  !  1.11  06/02/2006  Added logical flag for linear in tau approx (R Saunders)
  !  1.12  06/04/2006  Added Meghatropiques (R. Saunders)
  !  1.13  14/03/2007  Added units conversion constants
  !  1.14  16/05/2007  Added polarimetric sensor type (R Saunders)
  !  1.15  25/09/2007  Added maximum number of warnings for checkinput (P Brunel)
  !  1.16  11/10/2007  Remove zhusta* and zice* constants ( P.Marguinaud )
  !  1.17  07/12/2007  Remove maximum number of warnings for checkinput (P Brunel)
  !  1.18  12/12/2007  Added hard limits for trace gases (R Saunders)
  !  1.19  13/12/2007  Renamed linear_tau (R Saunders)
  !  1.20  01/11/2007  Added parameters for section length and AD/K code (A. Geer)
  !  1.21  16/01/2008  Facility to apply regression limits  (N. Bormann)
  !  1.22  04/03/2008  Made min hard limit > zero (R Saunders)
  !  1.23  14/04/2008  Added SSM/T2 (R Saunders)
  !  1.24  02/06/2008  Changed mixing ratio for CO (R Saunders)
  !  1.25  12/08/2008  Added SSMISZ for SSMIS chan19-22 - Zeeman (P. Rayer)
  !  1.26  29/01/2009  Add Kalpana and FY-3 (R Saunders)
  !  1.27  26/05/2009  Add more platforms and sensors (R Saunders)
  !  1.28  02/12/2009  Add principal component capability (Marco matricardi)
  !  1.29  15/01/2010  Add rttov9 intervals constants (P Marguinaud)
  !  1.30  05/07/2010  Add maximum solar zenith angle constant (J Hocking)
  !  1.31  01/02/2011  Updates to platform and sensor lists (J Hocking)
  !  1.32  09/11/2012  Add theta_eff for Lambertian refl (R SAunders)
  !  1.33  07/05/2013  Add ice crystal diameter and ice water content hard limits (J Vidot)
  !
  ! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
  !

!  Use parkind1, Only : jpim     ,jprb
  USE mo_kind,  only:wp
  Implicit None

  !1.0 Precision and numerical constants
  ! Try to ensure this is large enough to avoid overflows in reciprocals but small enough to not affect the calculations.
  ! these parameters are defined at bottom of module (because they use 
  Real(wp), parameter :: max_exp_exponent = 50._wp ! approx 1e22 which should be sufficiently big for most purposes
  Real(wp), parameter :: min_exponent     = 1e-16_wp ! approx log_10(1+2^-52) - anything raised to this power or smaller 
                                                         ! should be approx equal to 1
! small_val is defined in rttov_transmit to avoid compiler incompatibility
!   ! small_val is used in rttov_transmit to ensure small values do not result in underflows. In subsequent calculations
!   ! these values are multiplied together hence the exponent of 1/3.
!   Real(wp)            :: small_val = (tiny(min_exponent)) ** (0.333333_wp) ! XLF doesn't like 1/3

  !1.1 general
  !-----------
  ! Version number of the current code

  Integer*4, Parameter :: version = 11
  Integer*4, Parameter :: release = 1
  Integer*4, Parameter :: minor_version = 0

  Integer*4, Parameter :: version_compatible_min = 10 ! minimum version number
  Integer*4, Parameter :: version_compatible_max = 11 ! maximum version number
          ! compatible for coefficients.
          ! coef files with "id_comp_lvl" outside range will be rejected

  Character (len=16), Parameter :: rttov_magic_string = '%RTTOV_COEFF    '
  Real(wp),    Parameter :: rttov_magic_number = 1.2345E+12_wp

  Integer*4, Parameter :: default_err_unit = 0  ! standard error unit number
                              ! standard error unit number is 7 for HPUX

  !1.2 physical constants
  !----------------------
  ! Molecular weights  (g/mole) are calculated by adding NIST Standard Atomic Weights
  ! Molecular weight of dry air refers to US standard atmosphere 1976
  ! NIST  Standard Atomic Weight are:
  ! H    1.00794   (7)
  ! C   12.0107    (8)
  ! N   14.0067    (2)
  ! O   15.9994    (3)
  Real(wp), Parameter :: mair = 28.9644_wp
  Real(wp), Parameter :: mh2o = 18.01528_wp
  Real(wp), Parameter :: mo3  = 47.9982_wp
  Real(wp), Parameter :: mco2 = 44.0095_wp
  Real(wp), Parameter :: mch4 = 16.04246_wp
  Real(wp), Parameter :: mn2o = 44.0128_wp
  Real(wp), Parameter :: mco  = 28.0101_wp

  ! Avogadro constant from NIST (mol-1)
  Real(wp), Parameter :: na = 6.02214129E23_wp
  
  ! Gravity from NIST 9.80665 ms-1 (exact)
  Real(wp), Parameter :: gravity = 9.80665_wp

  !
  ! Kaye & Laby latest library edition is 16e 1995, and gives
  ! * standard value  g = 9.80665 ms-1 exactly (p.191)
  ! * earth mean radius r= 6371.00 km (p191)
  !    [defined as [(r_equator)^2 (r_pole)]^1/3]
  Real(wp), Parameter :: pi      = 3.1415926535_wp
  Real(wp), Parameter :: deg2rad = pi/180.0_wp
  Real(wp), Parameter :: earthradius = 6371.00_wp
  Real(wp), Parameter :: flatt       = 3.3528107E-3_wp
  Real(wp), Parameter :: omega       = 7292115E-11_wp
  Real(wp), Parameter :: eqrad       = 6378.137_wp
  Real(wp), Parameter :: grave       = 9.7803267715_wp
  Real(wp), Parameter :: z4pi_r      = 0.0795774715_wp
  Real(wp), Parameter :: pi_r        = 0.3183098862_wp
  Real(wp), Parameter :: theta_eff   = 55.0_wp
  Real(wp), Parameter :: sec_theta_eff   = 1.743446796_wp

  ! The Cosmic Microwave Background Spectrum from the Full COBE FIRAS Data Set
  ! Fixsen D.J. et all
  ! Astrophysical Journal v.473, p.576 December 1996
  ! CMBR = 2.728 +- 0.004K
  Real(wp), Parameter :: tcosmic     = 2.728_wp
  !  Real(wp), Parameter :: tcosmic     = 0.1_wp !used for ECMWF tests

  ! Universal gas constant R = 8.314510 J/mol/K
  Real(wp), Parameter :: rgp = 8.314510_wp
  Real(wp), Parameter :: rgc = 8.314472_wp

  ! mean molar mass of dry air rm = 0.0289644 kg.mol^-1
  Real(wp), Parameter :: rm = 0.0289644_wp

  ! units conversion from  mixing ratio to ppmv
  Real(wp), Parameter :: q_mixratio_to_ppmv  = 1.60771704e+6_wp
  Real(wp), Parameter :: o3_mixratio_to_ppmv = 6.03504e+5_wp
  Real(wp), Parameter :: co2_mixratio_to_ppmv= 6.58114e+5_wp
  Real(wp), Parameter :: co_mixratio_to_ppmv = 1.0340699e+6_wp
  Real(wp), Parameter :: n2o_mixratio_to_ppmv= 6.58090e+5_wp
  Real(wp), Parameter :: ch4_mixratio_to_ppmv= 1.80548e+6_wp

  ! zero temperature(K)
  Real(wp), Parameter :: t0 = 273.15_wp
  ! standard pressure
  Real(wp), Parameter :: p0 = 1013.25_wp

  !1.3 satellite and instrument information
  !----------------------------------------

  !platform id codes
  Integer*4, Parameter :: nplatforms = 35
  Integer*4, Parameter :: &
       & platform_id_noaa      = 1, &
       & platform_id_dmsp      = 2, &
       & platform_id_meteosat  = 3, &
       & platform_id_goes      = 4, &
       & platform_id_gms       = 5, &
       & platform_id_fy2       = 6, &
       & platform_id_trmm      = 7, &
       & platform_id_ers       = 8, &
       & platform_id_eos       = 9, &
       & platform_id_metop     = 10, &
       & platform_id_envisat   = 11, &
       & platform_id_msg       = 12, &
       & platform_id_fy1       = 13, &
       & platform_id_adeos     = 14, &
       & platform_id_mtsat     = 15, &
       & platform_id_coriolis  = 16, &
       & platform_id_jpss      = 17, &
       & platform_id_gifts     = 18, &
       & platform_id_sentinel3 = 19, &
       & platform_id_meghatr   = 20, &
       & platform_id_kalpana   = 21, &
       & platform_id_insat_3d  = 22, &
       & platform_id_fy3       = 23, &
       & platform_id_coms      = 24, &
       & platform_id_meteorm   = 25, &
       & platform_id_gosat     = 26, &
       & platform_id_calipso   = 27, &
       & platform_id_dummy     = 28, &
       & platform_id_gcomw     = 29, &
       & platform_id_nimbus    = 30, &
       & platform_id_himawari  = 31, &
       & platform_id_mtg       = 32, &
       & platform_id_saral     = 33, &
       & platform_id_metopng   = 34, &
       & platform_id_landsat   = 35

  !platform names
  Character (len=9), Parameter :: platform_name(nplatforms) = &
       & (/ 'noaa     ', 'dmsp     ', 'meteosat ', 'goes     ', 'gms      ', &
          & 'fy2      ', 'trmm     ', 'ers      ', 'eos      ', 'metop    ', &
          & 'envisat  ', 'msg      ', 'fy1      ', 'adeos    ', 'mtsat    ', &
          & 'coriolis ', 'jpss     ', 'gifts    ', 'sentinel3', 'meghatr  ', &
          & 'kalpana  ', 'insat_3d ', 'fy3      ', 'coms     ', 'meteor-m ', &
          & 'gosat    ', 'calipso  ', 'dummy    ', 'gcom-w   ', 'nimbus   ', &
          & 'himawari ', 'mtg      ', 'saral    ', 'metopng  ', 'landsat  ' /)

  !instrument id codes
  Integer*4, Parameter :: &
       & inst_id_hirs   =  0, inst_id_msu    =  1, inst_id_ssu    =  2, inst_id_amsua   =  3, &
       & inst_id_amsub  =  4, inst_id_avhrr  =  5, inst_id_ssmi   =  6, inst_id_vtpr1   =  7, &
       & inst_id_vtpr2  =  8, inst_id_tmi    =  9, inst_id_ssmis  = 10, inst_id_airs    = 11, &
       & inst_id_hsb    = 12, inst_id_modis  = 13, inst_id_atsr   = 14, inst_id_mhs     = 15, &
       & inst_id_iasi   = 16, inst_id_amsre  = 17, inst_id_gmsim  = 18, inst_id_atms    = 19, &
       & inst_id_mviri  = 20, inst_id_seviri = 21, inst_id_goesim = 22, inst_id_goessd  = 23, &
       & inst_id_mtsatim= 24, inst_id_vissr  = 25, inst_id_mvisr  = 26, inst_id_cris    = 27, &
       & inst_id_cmis   = 28, inst_id_viirs  = 29, inst_id_windsat= 30, inst_id_gifts   = 31, &
       & inst_id_ssmt1  = 32, inst_id_ssmt2  = 33, inst_id_saphir = 34, inst_id_madras  = 35, &
       & inst_id_ssmisz = 36, inst_id_kavhrr = 37, inst_id_iimager= 38, inst_id_isoundr = 39, &
       & inst_id_mwts   = 40, inst_id_mwhs   = 41, inst_id_iras   = 42, inst_id_mwri    = 43, &
       & inst_id_abi    = 44, inst_id_mi     = 45, inst_id_msumr  = 46, inst_id_tansofts= 47, &
       & inst_id_iir    = 48, inst_id_mwr    = 49, inst_id_dummyir= 50, inst_id_dummymw = 51, &
       & inst_id_dummyhi= 52, inst_id_dummypo= 53, inst_id_scams  = 54, inst_id_smmr    = 55, &
       & inst_id_ahi    = 56, inst_id_irs    = 57, inst_id_altika = 58, inst_id_iasing  = 59, &
       & inst_id_tm     = 60, inst_id_fci    = 61, inst_id_amsr1  = 62, inst_id_amsr2   = 63, &
       & inst_id_vissr2 = 64, inst_id_slstr  = 65

  Integer*4, Parameter :: ninst = 66
  ! List of instruments  !!!! HIRS is number 0
  Character (len=8), Dimension(0:ninst-1),parameter :: inst_name =       &
        & (/ 'hirs    ', 'msu     ', 'ssu     ', 'amsua   ', 'amsub   ',  &
           & 'avhrr   ', 'ssmi    ', 'vtpr1   ', 'vtpr2   ', 'tmi     ',  &
           & 'ssmis   ', 'airs    ', 'hsb     ', 'modis   ', 'atsr    ',  &
           & 'mhs     ', 'iasi    ', 'amsre   ', 'imager  ', 'atms    ',  &
           & 'mviri   ', 'seviri  ', 'imager  ', 'sounder ', 'imager  ',  &
           & 'vissr   ', 'mvisr   ', 'cris    ', 'cmis    ', 'viirs   ',  &
           & 'windsat ', 'gifts   ', 'ssmt1   ', 'ssmt2   ', 'saphir  ',  &
           & 'madras  ', 'ssmisz  ', 'kavhrr  ', 'iimager ', 'isoundr ',  &
           & 'mwts    ', 'mwhs    ', 'iras    ', 'mwri    ', 'abi     ',  &
           & 'mi      ', 'msumr   ', 'tansofts', 'iir     ', 'mwr     ',  &
           & 'dummyir ', 'dummymw ', 'dummyhi ', 'dummypo ', 'scams   ',  &
           & 'smmr    ', 'ahi     ', 'irs     ', 'altika  ', 'iasing  ',  &
           & 'tm      ', 'fci     ', 'amsr    ', 'amsr2   ', 'vissr   ',  &
           & 'slstr   ' /)


  !1.4 Coefficient file Section names
  !----------------------------------
  Integer*4, Parameter :: nsections = 43
  Integer*4, Parameter :: lensection = 34
  Character(len=lensection), Parameter :: section_types(nsections) = &
    & (/ 'IDENTIFICATION                    ', 'LINE-BY-LINE                      ', &
       & 'FAST_MODEL_VARIABLES              ', 'FILTER_FUNCTIONS                  ', &
       & 'FUNDAMENTAL_CONSTANTS             ', 'SSIREM                            ', &
       & 'FASTEM                            ', 'REFERENCE_PROFILE                 ', &
       & 'PROFILE_LIMITS                    ', 'FAST_COEFFICIENTS                 ', &
       & 'COEF_SUB_FILES                    ', 'GAZ_UNITS                         ', &
       & 'DIMENSIONS                        ', 'FREQUENCIES                       ', &
       & 'HYDROMETEOR                       ', 'CONVERSIONS                       ', &
       & 'EXTINCTION                        ', 'ALBEDO                            ', &
       & 'ASYMMETRY                         ', 'GAS_SPECTRAL_INTERVAL             ', &
       & 'TRANSMITTANCE_TRESHOLD            ', 'SOLAR_SPECTRUM                    ', &
       & 'WATER_OPTICAL_CONSTANT            ', 'WAVE_SPECTRUM                     ', &
       & 'AEROSOLS_PARAMETERS               ', 'AEROSOLS_COMPONENTS               ', &
       & 'WATERCLOUD_TYPES                  ', 'WATERCLOUD_PARAMETERS             ', &
       & 'ICECLOUD_TYPES                    ', 'HEXAGONAL_PARAMETERS              ', &
       & 'AGGREGATE_PARAMETERS              ', 'PRINCOMP_PREDICTORS               ', &
       & 'PRINCOMP_EIGENVECTORS             ', 'PRINCOMP_COEFFICIENTS             ', &
       & 'EMISSIVITY_COEFFICIENTS           ', 'PC_REFERENCE_PROFILE              ', &
       & 'PC_PROFILE_LIMITS                 ', 'INSTRUMENT_NOISE                  ', &
       & 'PLANCK_WEIGHTED                   ', 'SOLAR_FAST_COEFFICIENTS           ', &
       & 'README_SPECTRAL_RESPONSE_FUNCTION ', 'NLTE_RADIANCE_COEFS               ', &
       & 'PRESSURE_MODULATED_CELL           '/)

  !sensors id codes
  Integer*4, Parameter :: nsensors = 4
  Integer*4, Parameter :: &
       & sensor_id_ir     = 1, &
       & sensor_id_mw     = 2, &
       & sensor_id_hi     = 3, &
       & sensor_id_po     = 4

  !sensors names
  Character (len=2), Parameter :: sensor_name(nsensors) = &
       & (/ 'ir', 'mw', 'hi', 'po' /)

  ! these codes are for the instrument from the inst_name array
  Integer*4, Parameter :: sensor_id(0:ninst-1) = (/ &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_mw,  &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_mw,  &
    sensor_id_mw, sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_ir,  &
    sensor_id_mw, sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_mw,  &
    sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_ir,  &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_mw, sensor_id_ir,  &
    sensor_id_po, sensor_id_hi, sensor_id_mw, sensor_id_mw, sensor_id_mw,  &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_ir,  &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_ir,  &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_ir, sensor_id_mw,  &
    sensor_id_ir, sensor_id_mw, sensor_id_hi, sensor_id_po, sensor_id_mw,  &
    sensor_id_mw, sensor_id_ir, sensor_id_hi, sensor_id_mw, sensor_id_hi,  &
    sensor_id_ir, sensor_id_ir, sensor_id_mw, sensor_id_mw, sensor_id_ir,  &
    sensor_id_ir /)

  !gas id codes
  Integer*4, Parameter :: ngases_max = 8
  Integer*4, Parameter :: &
        & gas_id_mixed       = 1, &
        & gas_id_watervapour = 2, &
        & gas_id_ozone       = 3, &
        & gas_id_wvcont      = 4, &
        & gas_id_co2         = 5, &
        & gas_id_n2o         = 6, &
        & gas_id_co          = 7, &
        & gas_id_ch4         = 8

  !gas names
  Character (len=12), Parameter :: gas_name(ngases_max) = &
        & (/ 'Mixed_gases ', &
           & 'Water_vapour', &
           & 'Ozone       ', &
           & 'WV_Continuum', &
           & 'CO2         ', &
           & 'N2O         ', &
           & 'CO          ', &
           & 'CH4         ' /)

  !gas units
  Integer*4, Parameter :: ngases_unit = 2
  Integer*4, Parameter :: &
        & gas_unit_specconc  = 1, &
        & gas_unit_ppmv      = 2
  Character (len=12), Parameter :: gas_unit_name(ngases_unit) = &
        & (/ 'spec. concen', &
           & 'ppmv        '  /)


  !1.5 error reporting
  !-------------------
  !error status values
  Integer*4, Parameter :: errorstatus_success = 0
  Integer*4, Parameter :: errorstatus_fatal   = 1


  !1.6 surface types
  !-----------------
  Integer*4, Parameter :: nsurftype = 2
  Integer*4, Parameter :: surftype_land = 0
  Integer*4, Parameter :: surftype_sea = 1
  Integer*4, Parameter :: surftype_seaice = 2

  !1.7 water types
  !---------------
  Integer*4, Parameter :: nwatertype = 1
  Integer*4, Parameter :: watertype_fresh_water = 0
  Integer*4, Parameter :: watertype_ocean_water = 1

  !1.8 cloud emissivity
  !---------------------
  Integer*4, Parameter :: overlap_scheme = 2    ! overlap scheme
  ! 1 => Geleyn and Hollingsworth (1979)
  ! 2 => Raisanen (1998)

  !
  !1.9 Hard limits for control of input profile
  !--------------------------------------------
  ! Temperature
  Real(wp), Parameter :: tmax   = 400.0_wp       ! degK
  Real(wp), Parameter :: tmin   = 90.0_wp        ! degK
  ! Water Vapour
  Real(wp), Parameter :: qmax   = 0.60E+06_wp    ! ppmv 0.373_wp kg/kg
  Real(wp), Parameter :: qmin   = 0.1E-10_wp     ! ppmv
  ! Ozone
  Real(wp), Parameter :: o3max  = 1000.0_wp      ! ppmv  1.657E-3_wp kg/kg
  Real(wp), Parameter :: o3min  = 0.1E-10_wp     ! ppmv
  ! CO2
  Real(wp), Parameter :: co2max = 1000.0_wp      ! ppmv
  Real(wp), Parameter :: co2min = 0.1E-10_wp     ! ppmv
  ! CO
  Real(wp), Parameter :: comax  = 10.0_wp        ! ppmv
  Real(wp), Parameter :: comin  = 0.1E-10_wp     ! ppmv
  ! N2O
  Real(wp), Parameter :: n2omax = 10.0_wp        ! ppmv
  Real(wp), Parameter :: n2omin = 0.1E-10_wp     ! ppmv
  ! CH4
  Real(wp), Parameter :: ch4max = 50.0_wp        ! ppmv
  Real(wp), Parameter :: ch4min = 0.1E-10_wp     ! ppmv
  ! Cloud Liquid Water
  Real(wp), Parameter :: clwmax = 1.0_wp         ! kg/kg
  Real(wp), Parameter :: clwmin = 0.0_wp         ! kg/kg
  ! Surface Pressure
  Real(wp), Parameter :: pmax   = 1100.0_wp      ! surface pressure hPa
  Real(wp), Parameter :: pmin   = 400.0_wp       ! hPa
  ! Surface Wind
  Real(wp), Parameter :: wmax   =  100.0_wp      ! surface wind speed (m/s)
  ! Zenith Angle
  Real(wp), Parameter :: zenmax = 75.0_wp        ! zenith angle (Deg) = secant 3.86_wp
  Real(wp), Parameter :: zenmaxv9 = 84.0_wp      ! larger zenmax for v9 predictors
  ! Cloud Top Pressure
  Real(wp), Parameter :: ctpmax = 1100.0_wp      ! (hPa)
  Real(wp), Parameter :: ctpmin =   50.0_wp      ! (hPa)
  ! Magnetic field strength
  Real(wp), Parameter :: bemax = 0.7_wp          ! (Gauss)
  Real(wp), Parameter :: bemin = 0.2_wp          ! (Guass)
  ! Ice Crystal Diameter
  Real(wp), Parameter :: dgmin_hex =  12.2_wp	  ! (micron)
  Real(wp), Parameter :: dgmax_hex =  118.29_wp  ! (micron)
  Real(wp), Parameter :: dgmin_agg =  5.61_wp	  ! (micron)
  Real(wp), Parameter :: dgmax_agg =  166.46_wp  ! (micron)  
  ! Ice Water Content
  Real(wp), Parameter :: iwcmin_hex =  0.000608_wp ! (g.m-3)
  Real(wp), Parameter :: iwcmax_hex =  0.254639_wp ! (g.m-3)
  Real(wp), Parameter :: iwcmin_agg =  0.000235_wp ! (g.m-3)
  Real(wp), Parameter :: iwcmax_agg =  0.489046_wp ! (g.m-3)  


  !1.10  Maximum Optical Depth
  !--------------------------
  ! maximum value of optical depth for transmittance calculation
  ! e(-30) -> 10**-14
  ! e(-50) -> 10**-22
  Real(wp), Parameter  :: max_optical_depth = 50._wp

  !1.11  Maximum solar zenith angle for which to apply solar calculation
  !---------------------------------------------------------------------
  Real(wp), Parameter  :: max_sol_zen = 84._wp
  
  
  !2 RTTOV7 aux parameters
  !-------------------------
  Integer*4, Parameter :: fastem_sp = 5  ! max. number of fastem surface parameters
  Real(wp), Parameter    :: mwcldtp = 322.0_wp  ! Upper pressure level (HPa) for lwp calcs
  Real(wp), Parameter    :: pressure_top = 0.004985_wp ! Pressure of top level for
                                                ! Line/Line calculations (hPa)
  Real(wp) , Dimension(8), Parameter :: dcoeff =        &! Debye coefs
        & (/ 17.1252_wp, 134.2450_wp, 310.2125_wp,  5.667_wp,   &
          & 188.7979_wp,  80.5419_wp,   0.1157_wp,  4.8417_wp/)

  !2.1 Polarisation definitions
  !----------------------------
  ! == pol_id +1
  !   1 average of vertical and horizontal
  !   2 nominal vertical at nadir, rotating
  !      with view angle
  !   3 nominal horizontal at nadir, rotating
  !      with view angle
  !   4 vertical
  !   5 horizontal
  !   6 + 45 minus -45 (3rd stokes vector)
  !   7 left circular - right circular (4th stokes vector)
  Integer*4, Dimension(7), Parameter :: npolar_compute = &
   & (/ 2, 2, 2, 1, 1, 2, 4/)
  Integer*4, Dimension(7), Parameter :: npolar_return = &
   & (/ 1, 1, 1, 1, 1, 2, 4/)

  ! pol_v and pol_h give proportion of v and h pol to use in emissivity calculation
  ! pol_s3 adds the 3rd/4th stokes vectors
  Real(wp), Parameter :: pol_v(3,7) = Reshape( &
    & (/ 0.5_wp, 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, 1.0_wp, &
       & 0.0_wp, 1.0_wp, 0.0_wp, &
       & 1.0_wp, 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, 0.0_wp  /), (/3,7/) )
  Real(wp), Parameter :: pol_h(3,7) = Reshape( &
    & (/ 0.5_wp, 0.0_wp, 0.0_wp, &
       & 0.0_wp, 1.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, 1.0_wp, &
       & 0.0_wp, 0.0_wp, 0.0_wp, &
       & 1.0_wp, 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, 0.0_wp  /), (/3,7/) )
  Real(wp), Parameter :: pol_s3(0:1,7) = Reshape( &
    & (/ 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, &
       & 0.0_wp, 0.0_wp, &
       & 1.0_wp, 0.0_wp, &
       & 0.0_wp, 1.0_wp  /), (/2,7/) )

  !3 RTTOVSCATT aux parameters
  !---------------------------
  ! Minimum cloud cover processed by rttov_scatt
  Real(wp), Parameter :: ccthres = 0.05_wp
  ! Minimum single scattering albedo processed by rttov_scatt
  Real(wp), Parameter :: min_ssa = 1.0E-03_wp
  ! Rain density (g.cm-3)
  Real(wp), Parameter :: rho_rain = 1.0_wp
  ! Snow density (g.cm-3)
  Real(wp), Parameter :: rho_snow = 0.1_wp

  ! Flags to identify function in shared K/Adjoint routines
  Integer*4, Parameter :: adk_adjoint = 0
  Integer*4, Parameter :: adk_k       = 1

  !4 Parameters to compute refractive index of air
  !--------------------------------------------------------
  Real(wp), Parameter :: D1   =8341.87_wp
  Real(wp), Parameter :: D2   =2405955.0_wp
  Real(wp), Parameter :: D3   =130.0_wp
  Real(wp), Parameter :: D4   =15996.0_wp
  Real(wp), Parameter :: D5   =38.9_wp
  Real(wp), Parameter :: DCO2 =0.540_wp
  Real(wp), Parameter :: ED1  =96095.43_wp
  Real(wp), Parameter :: ED2  =0.601_wp
  Real(wp), Parameter :: ED3  =0.00972_wp
  Real(wp), Parameter :: ED4  =0.003661_wp
  Real(wp), Parameter :: EW1  =3.7345_wp
  Real(wp), Parameter :: EW2  =0.0401_wp
  Real(wp), Parameter :: HTOP =100.0_wp
  Real(wp), Parameter :: CTOM =1.0E-4_wp
  Real(wp), Parameter :: WAVER=1700.0_wp

  !5 RTTOV8_M_SCATT
  !--------------------------------------------------------
  Integer*4, Parameter :: naer_max = 13

  Integer*4, Parameter :: &
        & aer_id_inso       = 1, &
        & aer_id_waso       = 2, &
        & aer_id_soot       = 3, &
        & aer_id_ssam       = 4, &
        & aer_id_sscm       = 5, &
        & aer_id_minm       = 6, &
        & aer_id_miam       = 7, &
        & aer_id_micm       = 8, &
        & aer_id_mitr       = 9, &
        & aer_id_suso       =10, &
        & aer_id_vola       =11, &
        & aer_id_vapo       =12, &
        & aer_id_asdu       =13

  Character (len=4), Parameter :: aer_name(naer_max) = &
        & (/ 'inso', &
           & 'waso', &
           & 'soot', &
           & 'ssam', &
           & 'sscm', &
           & 'minm', &
           & 'miam', &
           & 'micm', &
           & 'mitr', &
           & 'suso', &
           & 'vola', &
           & 'vapo', &
           & 'asdu'  /)

  Integer*4, Parameter :: nphangle = 208

  Real(wp), Parameter :: phangle(nphangle) = &
             (/   0.0_wp,   0.1_wp,   0.2_wp,   0.3_wp,   0.4_wp,   0.5_wp,   0.6_wp, &
              &   0.7_wp,   0.8_wp,   0.9_wp,   1.0_wp,   1.1_wp,   1.2_wp,   1.3_wp, &
              &   1.4_wp,   1.5_wp,   1.6_wp,   1.7_wp,   1.8_wp,   1.9_wp,   2.0_wp, &
              &   2.1_wp,   2.2_wp,   2.3_wp,   2.4_wp,   2.5_wp,   2.6_wp,   2.7_wp, &
              &   2.8_wp,   2.9_wp,   3.0_wp,   4.0_wp,   5.0_wp,   6.0_wp,   7.0_wp, &
              &   8.0_wp,   9.0_wp,  10.0_wp,  11.0_wp,  12.0_wp,  13.0_wp,  14.0_wp, &
              &  15.0_wp,  16.0_wp,  17.0_wp,  18.0_wp,  19.0_wp,  20.0_wp,  21.0_wp, &
              &  22.0_wp,  23.0_wp,  24.0_wp,  25.0_wp,  26.0_wp,  27.0_wp,  28.0_wp, &
              &  29.0_wp,  30.0_wp,  31.0_wp,  32.0_wp,  33.0_wp,  34.0_wp,  35.0_wp, &
              &  36.0_wp,  37.0_wp,  38.0_wp,  39.0_wp,  40.0_wp,  41.0_wp,  42.0_wp, &
              &  43.0_wp,  44.0_wp,  45.0_wp,  46.0_wp,  47.0_wp,  48.0_wp,  49.0_wp, &
              &  50.0_wp,  51.0_wp,  52.0_wp,  53.0_wp,  54.0_wp,  55.0_wp,  56.0_wp, &
              &  57.0_wp,  58.0_wp,  59.0_wp,  60.0_wp,  61.0_wp,  62.0_wp,  63.0_wp, &
              &  64.0_wp,  65.0_wp,  66.0_wp,  67.0_wp,  68.0_wp,  69.0_wp,  70.0_wp, &
              &  71.0_wp,  72.0_wp,  73.0_wp,  74.0_wp,  75.0_wp,  76.0_wp,  77.0_wp, &
              &  78.0_wp,  79.0_wp,  80.0_wp,  81.0_wp,  82.0_wp,  83.0_wp,  84.0_wp, &
              &  85.0_wp,  86.0_wp,  87.0_wp,  88.0_wp,  89.0_wp,  90.0_wp,  91.0_wp, &
              &  92.0_wp,  93.0_wp,  94.0_wp,  95.0_wp,  96.0_wp,  97.0_wp,  98.0_wp, &
              &  99.0_wp, 100.0_wp, 101.0_wp, 102.0_wp, 103.0_wp, 104.0_wp, 105.0_wp, &
              & 106.0_wp, 107.0_wp, 108.0_wp, 109.0_wp, 110.0_wp, 111.0_wp, 112.0_wp, &
              & 113.0_wp, 114.0_wp, 115.0_wp, 116.0_wp, 117.0_wp, 118.0_wp, 119.0_wp, &
              & 120.0_wp, 121.0_wp, 122.0_wp, 123.0_wp, 124.0_wp, 125.0_wp, 126.0_wp, &
              & 127.0_wp, 128.0_wp, 129.0_wp, 130.0_wp, 131.0_wp, 132.0_wp, 133.0_wp, &
              & 134.0_wp, 135.0_wp, 136.0_wp, 137.0_wp, 138.0_wp, 139.0_wp, 140.0_wp, &
              & 141.0_wp, 142.0_wp, 143.0_wp, 144.0_wp, 145.0_wp, 146.0_wp, 147.0_wp, &
              & 148.0_wp, 149.0_wp, 150.0_wp, 151.0_wp, 152.0_wp, 153.0_wp, 154.0_wp, &
              & 155.0_wp, 156.0_wp, 157.0_wp, 158.0_wp, 159.0_wp, 160.0_wp, 161.0_wp, &
              & 162.0_wp, 163.0_wp, 164.0_wp, 165.0_wp, 166.0_wp, 167.0_wp, 168.0_wp, &
              & 169.0_wp, 170.0_wp, 171.0_wp, 172.0_wp, 173.0_wp, 174.0_wp, 175.0_wp, &
              & 176.0_wp, 177.0_wp, 178.0_wp, 179.0_wp, 180.0_wp /)

  Integer*4, Parameter :: nwcl_max = 5

  Integer*4, Parameter :: &
        & wcl_id_stco       = 1, &
        & wcl_id_stma       = 2, &
        & wcl_id_cucc       = 3, &
        & wcl_id_cucp       = 4, &
        & wcl_id_cuma       = 5

  Character (len=4), Parameter :: wcl_name(nwcl_max) = &
        & (/ 'stco', &
           & 'stma', &
           & 'cucc', &
           & 'cucp', &
           & 'cuma' /)

  Integer*4, Parameter:: ncldtyp = 6

  Real(wp), Parameter :: E00       = 611.21_wp
  Real(wp), Parameter :: T00       = 273.16_wp
  Real(wp), Parameter :: TI        = T00 - 23.0_wp

  Real(wp), Parameter :: min_tau = 1.0e-8_wp
  Real(wp), Parameter :: min_od  = 1.0e-5_wp

!
! These are the RTTOV9 wavenumbers that make intervals
!
  Real(wp), Parameter :: rttov9_wv0690_50 =  690.50_wp, &
                                rttov9_wv1050_00 = 1050.00_wp, &
                                rttov9_wv1095_25 = 1095.25_wp, &
                                rttov9_wv1100_25 = 1100.25_wp, &
                                rttov9_wv1350_25 = 1350.25_wp, &
                                rttov9_wv1750_25 = 1750.25_wp, &
                                rttov9_wv1900_25 = 1900.25_wp, &
                                rttov9_wv1995_00 = 1995.00_wp, &
                                rttov9_wv2000_00 = 2000.00_wp, &
                                rttov9_wv2250_00 = 2250.00_wp, &
                                rttov9_wv2295_25 = 2295.25_wp, &
                                rttov9_wv2360_00 = 2360.00_wp, &
                                rttov9_wv2380_25 = 2380.25_wp, &
                                rttov9_wv2660_25 = 2660.25_wp, &
                                rttov9_wv2760_25 = 2760.25_wp

!
! Parameters for solar overcast radiance calculation
!

Real(wp), Parameter :: overcast_albedo_wvn = 10000._wp ! Wavenumber (cm-1) at which albedo changes
Real(wp), Parameter :: overcast_albedo1    = 0.7_wp    ! Overcast albedo for wvn > limit
Real(wp), Parameter :: overcast_albedo2    = 0.6_wp    ! Overcast albedo for wvn <= limit

!
! Parameters for Rayleigh cross-section parameterization taken from Bucholzt 1995.
!
Real(wp), Parameter :: ray_min_wvn = 5000.0_wp,     & ! Min wavenumber (cm-1) for which Rayleigh is calculated
                              ray_scs_wlm = 0.5_wp,        & ! Wavelength limit: below 0.5um the
                              ray_scs_a1 = 3.01577E-28_wp, & !   first set of parameters a1-d1 are used
                              ray_scs_b1 = -3.55212_wp,    & !   while above this a2-d2 are used.
                              ray_scs_c1 = -1.35579_wp,    &
                              ray_scs_d1 = -0.11563_wp,    &
                              ray_scs_a2 = 4.01061E-28_wp, &
                              ray_scs_b2 = -3.99668_wp,    &
                              ray_scs_c2 = -1.10298E-3_wp, &
                              ray_scs_d2 = -2.71393E-2_wp

!
! Interpolation modes
!
Integer*4, Parameter :: ninterp_modes           = 3 ! Number of valid interpolation options
Integer*4, Parameter :: interp_rochon           = 1 ! Rochon interpolation everywhere
Integer*4, Parameter :: interp_loglinear        = 2 ! Log-linear interpolation everywhere
Integer*4, Parameter :: interp_rochon_loglinear = 3 ! Rochon for user->coef, log-lin. for coef->user

END MODULE mo_cosp_rttov_const
