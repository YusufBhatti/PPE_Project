!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_bethy_constants

  ! 
  !! Description: 
  !!   Declaration of additional constants used in BETHY
  !!   The constants from mo_jsbach_constants are also available from this module
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
  USE mo_kind, ONLY : dp

  IMPLICIT NONE 

  !! Global (i.e. public) Declarations: 

  !! Global Parameters: 

  REAL(dp),    PARAMETER :: LaiMax = 8._dp                  !! Maximum LAI (used for nitrogen scaling) 
  REAL(dp),    PARAMETER :: LaiLimit = 3._dp                !! Minimum LAI: used for nitrogen scaling and estimation of fractional
                                                            !!    cover
  REAL(dp),    PARAMETER :: LaiMin = 1.E-9_dp               !! Minimum Lai in PAR computation
  REAL(dp),    PARAMETER :: EPar = 2.2E5_dp                 !! Energy content of PAR [J / mol(photons)]=(4.6 mol/MJ PAR)**-1
  REAL(dp),    PARAMETER :: FcMax = 1.0_dp                  !! Maximum fractional vegetation cover
  REAL(dp),    PARAMETER :: FcMin = 1.E-3_dp                !! Minimum fractional vegetation cover
  REAL(dp),    PARAMETER :: ZenithMin = 0.0174524_dp        !! Check for solar zenith angle > 89 degrees
  REAL(dp),    PARAMETER :: ZenithMinPar = 1.E-3_dp         !! Minimum cos of zenith angle for which PAR is calculated
  REAL(dp),    PARAMETER :: SoilReflectivityParMin = 0.0_dp !! Minimum soil reflectivity in PAR region

  REAL(dp),    PARAMETER :: minOfMaxCarboxrate = 1.0e-12_dp !! Minimum of maximum carboxylation rate [10^(-6) mol/(m^2 s)]
  REAL(dp),    PARAMETER :: minStomaConductance = 0.0_dp    !! Minimum stomatal conductance [mol H2O /(m^2 s) ??]

  !! Factors that relates leaf internal CO2-concentration to CO2-concentration of ambient air: 
  REAL(dp), PARAMETER    :: FCI1C3        = 0.87_dp         !! For C3 plants
  REAL(dp), PARAMETER    :: FCI1C4        = 0.67_dp         !! For C4 plants


END MODULE mo_bethy_constants
