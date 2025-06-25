!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_jsbach_constants

  ! 
  !! Description: 
  !!   Declaration of constants used in JSBACH
  ! 
  !! Current Code Owner: jsbach_admin
  ! 
  !! History: 
  !  
  !! Version   Date        Comment 
  !! -------   ----        ------- 
  !! 0.1       2001/06/28  Original code. Reiner Schnur
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  ! Modules used: 
  ! 
  USE mo_kind,                 ONLY: dp
  USE mo_physical_constants,   ONLY: alf, als, alv, amco2, amd, amn2o, amw, argas, cpd, cpv, grav, &
                                     stbo, rd, rdv, rhoh2o, rv, tmelt, vtmpc1, vtmpc2, secperday, p0sl_bg
  USE mo_math_constants,       ONLY: pi, euler
  USE mo_radiation_parameters, ONLY: cemiss

  IMPLICIT NONE

  !! Global Parameters:

! REAL(dp),  PARAMETER :: Pi = pi                                !!  3.14159265358979323846264338327950288419717_dp
! REAL(dp),  PARAMETER :: euler = euler                          !!  2.71828182845904523536028747135266249775725_dp

! REAL(dp),  PARAMETER :: Tmelt  = tmelt                         !!   273.15_dp - Melting temperature of snow/ice [K]
! REAL(dp),  PARAMETER :: RhoH2O = rhoh2o                        !!    1000._dp - Density of liquid water (kg/m^3)

  REAL(dp),  PARAMETER :: Gravity = grav                         !!  9.80665_dp - av. gravitational acceleration at surface [m/s^2]

  REAL(dp),  PARAMETER :: UniversalGasConstant = argas           !! 8.314472_dp - Universal gas constant [J/(mole*K)]
  REAL(dp),  PARAMETER :: GasConstantDryAir = rd                 !!   287.04_dp - gas constant for dry air [J/(K*kg)]
  REAL(dp),  PARAMETER :: GasConstantWaterVapor = rv             !!   461.51_dp - gas constant for water vapor [J/(K*kg)]
  REAL(dp),  PARAMETER :: SpecificHeatDryAirConstPressure = cpd  !!  1004.64_dp - specific heat of dry air at constant pressure
  REAL(dp),  PARAMETER :: SpecificHeatVaporConstPressure =  cpv  !!  1869.46_dp - specific heat of water vapor at constant pressure
  REAL(dp),  PARAMETER :: eps = rdv                              !!             - rd/rv             [ ]
! REAL(dp),  PARAMETER :: vtmpc1 = vtmpc1                        !!             - rv/rd   - 1.0_dp  [ ]
! REAL(dp),  PARAMETER :: vtmpc2 = vtmpc2                        !!             - cpv/cpd - 1.0_dp  [ ]

  REAL(dp),  PARAMETER :: LatentHeatVaporization = alv           !! 2.5008e6_dp - latent heat for vaporization [J/kg]
  REAL(dp),  PARAMETER :: LatentHeatSublimation  = als           !! 2.8345e6_dp - latent heat for sublimation [J/kg]
  REAL(dp),  PARAMETER :: LatentHeatFusion = alf                 !!   als-alv   - latent heat of fusion [J/kg]

  REAL(dp),  PARAMETER :: solar_const =  1361.371_dp             !!              default solar constant (AMIP) [W/m2]

  REAL(dp),  PARAMETER :: StefanBoltzmann = stbo                 !! 5.6704e-8_dp - Stefan-Boltzmann constant [W/(m^2 K^4)] 
  REAL(dp),  PARAMETER :: vonKarman = 0.4_dp                     !!               von Karman constant for evapotranspiration

  REAL(dp),  PARAMETER :: Emissivity = cemiss                    !! 0.996_dp    - surface emissivity

  REAL(dp),  PARAMETER :: molarMassDryAir_kg = amd   * 1.e-3_dp  !!  28.970e-3_dp   -  Mass of 1 mol of dry air in kg
  REAL(dp),  PARAMETER :: molarMassC_kg      = 12.01_dp * 1.e-3_dp  !! 12.01e-3_dp  -  Mass of 1 mod C   in kg
  REAL(dp),  PARAMETER :: molarMassCO2_kg    = amco2 * 1.e-3_dp  !! 44.0095e-3_dp   -  Mass of 1 mol CO2 in kg
  REAL(dp),  PARAMETER :: molarMassN2O_kg    = amn2o * 1.e-3_dp  !!  44.013e-3_dp   -  Mass of 1 mol CO2 in kg 
  REAL(dp),  PARAMETER :: molarMassH2O_kg    = amw   * 1.e-3_dp  !! 18.0154e-3_dp   -  Mass of 1 mol H2O in kg

  REAL(dp),  PARAMETER :: p_sealevel  = p0sl_bg                  !! 101325._dp    - sea level pressure [Pa]

  REAL(dp),  PARAMETER :: snow_density = 330.0_dp                !!              snow density  [kg/m**3]
 
  ! non-physical model parameters
  REAL(dp),  PARAMETER :: cvinter = 0.25_dp                      !! Efficiency of interception of precipitation as rain

  ! Parameters used for snow cover fraction
  REAL(dp),  PARAMETER :: zepsec  = 1.E-12_dp
  REAL(dp),  PARAMETER :: zsigfac = 0.15_dp
  REAL(dp),  PARAMETER :: zqsncr  = 0.95_dp                      !! inverse of equivalent water height when snow is considered to
                                                                 !!    completely cover the ground

  ! Ratio of C mass to biomass (Spitfire and emission of chemical species by fire)
  REAL(dp),  PARAMETER :: C_2_biomass_ratio = 0.5_dp

END MODULE mo_jsbach_constants
