!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_helper
  
  USE mo_kind,        ONLY: dp 
  IMPLICIT NONE
  
  ! PARAMETERS FOR AUXILLIARY FUNCTIONS
  REAL(dp), PARAMETER :: zmin = 1e-18
!  REAL, PARAMETER :: eta = 0.99999            ! curvature parameter for mins/maxs
  
CONTAINS
    
  !*********************************************************
  !* FUNCTION errf
  !* the (cumulative) error function
  !* numerical recipes in Fortran 77, Chapter 6.2
  !*********************************************************
  
  REAL(DP) FUNCTION errf (x)
    REAL(dp) :: x, z, t
    z=ABS(x) 
    t=1./(1.+0.5*z) 
    errf=1.-0.5*t*EXP(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
         t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
         t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
    IF (x<0) errf=1.-errf
  END FUNCTION errf

  !*********************************************************
  !*  FUNCTION mins
  !*  smoothed minimum function
  !*********************************************************

!  REAL(DP) FUNCTION mins (x, y)
!    REAL :: x, y
!    REAL :: z
!    z = (x+y)**2 - 4.*eta*x*y
!    IF (z.GE.zmin) THEN
!       mins = (x + y - SQRT(z)) / (2.*eta)
!    ELSE
!       mins = 0.
!    ENDIF
!  END FUNCTION mins
  REAL(DP) FUNCTION mins (x, y, eta)
    REAL(dp) :: x, y, eta
    REAL(dp) :: z
    z = (x+y)**2 - 4._dp*eta*x*y
    z = max (z, zmin)
    mins = (x + y - SQRT(z)) / (2.*eta)
  END FUNCTION mins

  !*********************************************************
  !*  FUNCTION maxs
  !*  smoothed maximum function
  !*********************************************************

!  REAL(DP) FUNCTION maxs (x, y)
!    REAL :: x, y
!    REAL :: z
!    z = (x+y)**2 - 4.*eta*x*y
!    IF (z.GE.zmin) THEN
!       maxs = (x + y + SQRT(z)) / (2.*eta)
!    ELSE
!       maxs = 0.
!    ENDIF
!  END FUNCTION maxs
  REAL(DP) FUNCTION maxs (x, y, eta)
    REAL(dp) :: x, y, eta
    REAL(dp) :: z
    z = (x+y)**2 - 4.*eta*x*y
    z = max (z, zmin)
    maxs = (x + y + SQRT(z)) / (2.*eta)
  END FUNCTION maxs

  !*********************************************************
  !*  FUNCTION minx
  !*  minimum function with exponential transition
  !*********************************************************

  REAL(DP) FUNCTION minx (x, y, x0)
    REAL(dp) :: x, y, x0
    IF (x.LE.y+x0) THEN
       minx = x - x0*EXP((x-y)/x0-1.)
    ELSE
       minx = y
    ENDIF
  END FUNCTION minx

  !*********************************************************
  !*  FUNCTION maxx
  !*  maximum function with exponential transition
  !*********************************************************

  REAL(DP) FUNCTION maxx (x, y, x0)
    REAL(dp) :: x, y, x0
    IF (x.GE.y-x0) THEN
       maxx = x + x0*EXP(-(x-y)/x0-1.)
    ELSE
       maxx = y
    ENDIF
  END FUNCTION maxx

END MODULE mo_helper
