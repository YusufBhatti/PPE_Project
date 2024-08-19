!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_parameters

  IMPLICIT NONE

  ! parameters controlling array sizes.

!LK  INTEGER, PARAMETER :: jpm     = 106 ! max zonal wave number
!LK  INTEGER, PARAMETER :: jpn     = jpm ! max meridional wave number for m=0
!LK  INTEGER, PARAMETER :: jpk     = jpm ! max meridional wave number
                                      ! for ntrn only, f90 namelist restriction
  INTEGER, PARAMETER :: jpgrnd  = 5   ! number of soil layers
                                      
!LK  INTEGER, PARAMETER :: jptrac  = 100 ! maximum number of prognostic tracers
!LK  INTEGER, PARAMETER :: jps     = 3   ! basic Spitfire variables without tracers
                                      ! i.e. humidity, cloud water and cloud ice
!LK  INTEGER, PARAMETER :: jpmp1=jpm+1

END MODULE mo_parameters
