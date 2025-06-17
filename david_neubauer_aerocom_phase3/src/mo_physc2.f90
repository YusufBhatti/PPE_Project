!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_physc2

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_physc2* constants to communicate between the main program
  !                    and the physical subroutines (except radiation ones).
  !
  ! ----------------------------------------------------------------

  REAL(dp) :: clam            !  *asymptotic mixing length for momentum.
  REAL(dp) :: ckap            !  *karman constant.
  REAL(dp) :: cb              !  *stability parameter near neutrality.
  REAL(dp) :: cc              !  *stability parameter for unstable cases.
  REAL(dp) :: cd              !  *stability parameter for stable cases.
  REAL(dp) :: cchar           !  *charnock constant.
  REAL(dp) :: cfreec          !  *free convection parameter
  REAL(dp) :: cgam            !  *free convection parameter
  REAL(dp) :: crgcg           !  *heat capacity x density of the soil.
  REAL(dp) :: cdif            !  *thermal diffusivity of soil.
  REAL(dp) :: cwsat           !  *soil water content at saturation
  REAL(dp) :: cwcap           !  *soil water content at field capacity (=wsmax)
  REAL(dp) :: cwcrit          !  *soil water content above which evaporation
                              !   occurs at potential rate
  REAL(dp) :: cwpwp           !  *soil water content at wilting point
  REAL(dp) :: cgh2o           !  *heat capacity of water
  REAL(dp) :: cqwevap         !  *inverse of (cwcrit-cwpwp)
  REAL(dp) :: cqwsbcr         !  *inverse of critical wetness for bare ground
  REAL(dp) :: cqsncr          !  *inverse of equivalent water height when snow
                              !  is considered to cover completely the ground
                              !  in the box.
  REAL(dp) :: csncri          !  CRITICAL SNOW DEPTH FOR SOIL COMPUTATIONS
  REAL(dp) :: cwlmax          !  *maximum moisture content of
                              !   the skin reservoir
  REAL(dp) :: cvdifts         !   factor for timestep weighting
                              !   in *rhs* of *vdiff* and *scv*.
  REAL(dp), ALLOCATABLE :: cevapcu(:) !   *evaporation coefficient for kuo0.
  REAL(dp) :: ctfreez         !   temperature at which sea
                              !   starts freezing/melting
  REAL(dp) :: cz0ice          !   roughness over sea-ice
  REAL(dp) :: corvari         !  *minimum sub-grid scale orography variance
                              !   to allow surface runoff
  REAL(dp) :: corvars         !  *sub-grid scale orography variance
                              !   calibration in surface runoff

END MODULE mo_physc2
