!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_rad_fastmath
  USE mo_kind, ONLY : wp, dp
  
  IMPLICIT NONE
  PRIVATE 
  PUBLIC :: setup_rad_fastmath, expon, inv_expon, transmit, tautrans
  
  INTEGER,  PARAMETER :: ntbl = 10000
  REAL(wp), PARAMETER :: tblint = REAL(ntbl, wp)
  REAL(wp) :: tau_tbl(0:ntbl) !< Optical depth 
  REAL(wp) :: exp_tbl(0:ntbl) !< Exponential lookup table (EXP(-tau))
  REAL(wp) :: tfn_tbl(0:ntbl) !< Tau transition function
    ! i.e. the transition of the Planck function from that for the mean layer temperature 
    ! to that for the layer boundary temperature as a function of optical depth. 
    ! The "linear in tau" method is used to make the table.

  REAL(wp), PARAMETER :: od_lo  = 0.06_wp     !< Value of tau below which expansion is used
  REAL(wp), PARAMETER :: expeps = 1.e-20_wp   !< Smallest value for exponential table
  REAL(wp), PARAMETER :: pade   = 0.278_wp    !< Pade approximation constant
  REAL(wp), PARAMETER :: bpade  = 1._wp/pade
  REAL(wp), PARAMETER :: rec_6  = 1._wp/6._wp
  
  !
  ! The RRTMG tables are indexed with INT(tblint * x(i)/(bpade + x(i)) + 0.5_wp) 
  !   But these yield unstable values in the SW solver for some parameter sets, so 
  !   we'll remove this option (though the tables are initialized if people want them). 
  ! RRTMG table lookups are approximated second-order Taylor series expansion 
  !   (e.g. exp(-x) = 1._wp - x(1:n) + 0.5_wp * x(1:n)**2, tautrans = x/6._wp) for x < od_lo
  !

  LOGICAL, SAVE :: initialized = .false. 
  
CONTAINS
  ! ------------------------------------------------------------
  SUBROUTINE setup_rad_fastmath
    !
    ! Initialize lookup tables used by RRTMG
    !
    INTEGER  :: itr
    REAL(WP) :: tfn

    IF(.not. initialized) THEN 
	  ! Compute lookup tables for transmittance, tau transition function,
	  ! and clear sky tau (for the cloudy sky radiative transfer).  Tau is 
	  ! computed as a function of the tau transition function, transmittance 
	  ! is calculated as a function of tau, and the tau transition function 
	  ! is calculated using the linear in tau formulation at values of tau 
	  ! above 0.01.  TF is approximated as tau/6 for tau < 0.01.  All tables 
	  ! are computed at intervals of 0.001.  The inverse of the constant used
	  ! in the Pade approximation to the tau transition function is set to b.
  
	  tau_tbl(0) = 0.0_wp
	  exp_tbl(0) = 1.0_wp
	  tfn_tbl(0) = 0.0_wp
	  
	  DO itr = 1, ntbl-1
		tfn = float(itr) / float(ntbl)
		tau_tbl(itr) = bpade * tfn / (1._wp - tfn)
		exp_tbl(itr) = MAX(EXP(-tau_tbl(itr)), expeps)
		IF (tau_tbl(itr) .LT. od_lo) THEN
		  tfn_tbl(itr) = tau_tbl(itr) * rec_6
		ELSE
		  tfn_tbl(itr)=1._wp-2._wp*((1._wp/tau_tbl(itr))-(exp_tbl(itr)/(1.-exp_tbl(itr))))
		ENDIF
	  ENDDO
	  tau_tbl(ntbl) = 1.e10_wp
	  exp_tbl(ntbl) = expeps
	  tfn_tbl(ntbl) = 1.0_wp
      initialized = .true. 
    END IF 

  END SUBROUTINE setup_rad_fastmath
  ! ------------------------------------------------------------
  FUNCTION expon(x, n) 

    INTEGER,  INTENT(in) :: n
    REAL(wp), INTENT(in) :: x(n)
    REAL(wp)             :: expon(n)

#ifndef HAVE_FAST_MATH_LIB
    expon(1:n) = exp(x(1:n))
#else
#if (defined HAVE_ACML)
    CALL vrda_exp(n, x(1:n), expon(:))
#elif (defined HAVE_MASS)
    CALL vexp(expon(:), x(1:n), n)
#elif (defined HAVE_MKL)
    CALL vdexp(n, x(1:n), expon(:))
#else 
    ! no fast math library call available, use native Fortran ...
    expon(:) = exp(x(1:n))
#endif
#endif
  END FUNCTION expon
  ! ------------------------------------------------------------
  FUNCTION inv_expon(x,n) 
    !
    ! Compute EXP(-x) - but do it fast
    !
    INTEGER,  INTENT(in) :: n
    REAL(DP), INTENT(IN) :: x(n) 
    REAL(DP)             :: inv_expon(n) 
    
#ifndef HAVE_FAST_MATH_LIB
      inv_expon(1:n) = EXP(-x(1:n))
#else
#if (defined HAVE_ACML)
      CALL vrda_exp(n, -x(1:n), inv_expon(1:n))
#elif (defined HAVE_MASS)
      CALL vexp(inv_expon(1:n), -x(1:n), n)
#elif (defined HAVE_MKL)
      CALL vdexp(n, -x(1:n), inv_expon(1:n))
#else
      inv_expon(1:n) = EXP(-x(1:n))
     ! exp - no fast math library call available, use native Fortran ...
#endif
#endif

  END FUNCTION inv_expon
  ! ------------------------------------------------------------
  FUNCTION transmit(x, n) 
    !
    ! Compute transmittance 1 - EXP(-x) 
    ! 
    INTEGER,  INTENT(in) :: n
    REAL(DP), INTENT(IN) :: x(n) 
    REAL(DP)             :: transmit(n) 
    
    !
    ! MASS and MKL libraries have exp(x) - 1 functions; we could 
    !   use those here
    !
    transmit(1:n) = 1._wp - inv_expon(x,n)
    
  END FUNCTION transmit
  ! ------------------------------------------------------------
  FUNCTION tautrans(x, n) 
    !
    ! Compute "tau transition" using linear-in-tau approximation
    !
    INTEGER,  INTENT(in) :: n
    REAL(DP), INTENT(IN) :: x(n) 
    REAL(DP)             :: tautrans(n) 
    
    LOGICAL :: thin(n)
    REAL(DP) :: y(n) 
    ! 
    ! Default calculation is unstable (NaN) for the very lowest value of tau (3.6e-4) 
    !
    y(:) = inv_expon(x,n)
    tautrans(:)= MERGE(1._wp - 2._wp*(1._wp/x(1:n) - y(:)/(1._WP-y(:))), & 
                       x * rec_6,                                        & 
                       x > 1.E-3_wp)
  END FUNCTION tautrans
  ! ------------------------------------------------------------

END MODULE mo_rad_fastmath
