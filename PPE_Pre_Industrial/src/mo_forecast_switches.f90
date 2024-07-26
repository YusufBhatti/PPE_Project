!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_forecast_switches

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_forecast_switches* switches related to the dynamics
  !                         and the general control of the forecast.
  !
  ! ----------------------------------------------------------------

  LOGICAL :: lsimdt   !  *true for semi implicit time scheme for
                      !  divergence,temperature and surface pressure
                      !  equations.
  LOGICAL :: lsimzq   !  *true for semi implicit time scheme for
                      !  vorticity and humidity equations.
  LOGICAL :: lvtmpc1  !  *true for virtual temperature.
                      !  i.e:*rv*.ne.*rd*.
  LOGICAL :: lvtmpc2  !  *true for influence of humidity on *cp*.
                      !  i.e:*cpv*.ne.*cpd*.
  LOGICAL :: lumax    !  *true to compute and print information on
                      !  maximum wind.
  LOGICAL :: lzondia

END MODULE mo_forecast_switches
