!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_advection

  ! L. Kornblueh, MPI, August 2001, initial coding
  ! L. Kornblueh, MPI, Februrary 2012, remove spitfire
  
  IMPLICIT NONE

  PUBLIC

  ! defined advection schemes for passive tracer (in advection sense)

  INTEGER, PARAMETER :: no_advection    = 0
  INTEGER, PARAMETER :: semi_lagrangian = 1
  INTEGER, PARAMETER :: tpcore          = 3

  ! select variable, set in namelist runctl 

  INTEGER :: iadvec

  ! general dimension parameters to be set in the init section of
  ! the diferent advection modules 

  INTEGER :: nalatd, nalond, nalev, nacnst, nalat, nalon

  INTEGER, PARAMETER :: jps     =   3 ! basic advected variables without tracers:
                                      ! water vapour, cloud water, cloud ice
  
END MODULE mo_advection
