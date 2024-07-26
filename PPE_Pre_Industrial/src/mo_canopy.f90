!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_canopy

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! Parameters for computation of unstressed canopy resistance using Eq. 3.3.2.12 in ECHAM3 manual
  REAL(dp), PARAMETER :: &
       k = 0.9_dp,      &
       a = 5000._dp,    &
       b = 10._dp,      &
       c = 100._dp

CONTAINS
  
  SUBROUTINE unstressed_canopy_cond_par(lai, par, resistance)


    ! Computes the canopy conductance without water stress limitation according to ECHAM formalism, Eq. 3.3.2.12 in 
    ! ECHAM3 Manual

    REAL(dp), INTENT(in), DIMENSION(:) :: &
         lai,      &
         par
    REAL(dp), INTENT(out) :: resistance(SIZE(lai,DIM=1))
  

    ! Local variables
    REAL(dp) :: d(SIZE(par))    !! Variable d in Eq. 3.3.2.12 in ECHAM3 manual
    REAL(dp) :: par_echam(SIZE(par,DIM=1))
    par_echam = max(1.e-10_dp, par)

    WHERE(lai > EPSILON(1._dp))
       d = (a + b*c) / (c * par_echam)
       resistance = (LOG((d * EXP(k*lai) + 1._dp) / (d+1._dp)) * b / (d*par_echam) - LOG((d + EXP(-k*lai)) / (d+1._dp)))/(k*c)
    ELSEWHERE
       resistance = 1.e-20_dp
    END WHERE
    
  END SUBROUTINE unstressed_canopy_cond_par
!!$  SUBROUTINE unstressed_canopy_resistance_par(lai, par, resistance)
!!$
!!$
!!$    ! Computes the canopy resistance without water stress limitation according to ECHAM formalism, Eq. 3.3.2.12 in 
!!$    ! ECHAM3 Manual
!!$
!!$    REAL, INTENT(in), DIMENSION(:) :: &
!!$         lai,      &
!!$         par
!!$    REAL, INTENT(out) :: resistance(SIZE(lai,DIM=1))
!!$
!!$
!!$    ! Local variables
!!$    REAL :: d(SIZE(par))    !! Variable d in Eq. 3.3.2.12 in ECHAM3 manual
!!$
!!$    WHERE(par > 10.e-10 .AND. lai > EPSILON(1.))
!!$       d = (a + b*c) / (c * par)
!!$       resistance = ( k * c ) / (LOG((d * EXP(k*lai) + 1.) / (d+1.)) * b / (d*par) - LOG((d + EXP(-k*lai)) / (d+1)))
!!$    ELSEWHERE
!!$       resistance = HUGE(1.)
!!$    END WHERE
!!$    
!!$  END SUBROUTINE unstressed_canopy_resistance_par

!!$  SUBROUTINE unstressed_canopy_resistance_rmin(lai, rmin, resistance)
!!$
!!$    ! Computes the canopy resistance without water stress limitation according to VIC formalism, Eq. 2.13 in Xu Liang's thesis
!!$
!!$    REAL, INTENT(in), DIMENSION(:) :: &
!!$         lai,      &
!!$         rmin
!!$    REAL, INTENT(out) :: resistance(:)
!!$
!!$    WHERE(lai > EPSILON(1.))
!!$       resistance = rmin / lai
!!$    ELSEWHERE
!!$       resistance = HUGE(1.)
!!$    END WHERE
!!$
!!$  END SUBROUTINE unstressed_canopy_resistance_rmin

END MODULE mo_canopy
