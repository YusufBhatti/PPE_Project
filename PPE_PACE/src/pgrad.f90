#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE pgrad(pgrd,pt,kdim,plnps,klen)

  ! Description:
  !
  ! Calculate linearized geopotential and pressure gradient terms
  ! for specified temperature and surface pressure fields.
  !
  ! Method:
  !
  ! *pgrad* is called from *inhysi* in the initial calculation
  ! of the gravity-wave matrix *bb,* and from subroutine for
  ! the calculation of semi-implicit corrections. 
  !
  ! Parameters are:
  !   *pgrd*      *computed sum.
  !   *pt*        *specified temperature.
  !   *kdim*      *first dimension of 2-d arrays *pgrd* and *pt*
  !   *plnps*     *specified logarithm of surface pressure.
  !   *klen*      *number of points for which calculation is performed.
  !
  ! Required constants are obtained from module *mo_hyb*
  !
  ! Results:
  ! Results are computed for *klen* consecutive points at
  ! each model level.
  !
  ! Reference:
  ! External documentation of the model equations and
  ! the organization of the vertical calculation.
  !
  ! Authors:
  !
  ! A. J. Simmons, ECMWF, November 1981, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,    ONLY: dp
  USE mo_control, ONLY: nlev
  USE mo_hyb,     ONLY: nlevm1, ralphr, rdtr, rlnmar

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kdim, klen

  !  Array arguments 
  REAL(dp) :: pgrd(kdim,*), plnps(*), pt(kdim,*)

  !  Local scalars: 
  REAL(dp) :: za, zb
  INTEGER :: ikp, jk, jl


  !  Executable statements 

!-- 1. Integrate hydrostatic equation

  za = ralphr(nlev)
  DO jl = 1, klen
    pgrd(jl,nlev) = za*pt(jl,nlev) + rdtr*plnps(jl)
  END DO

  DO jk = nlevm1, 1, -1
    ikp = jk + 1
    za = ralphr(jk)
    zb = rlnmar(ikp)

    DO jl = 1, klen
      pgrd(jl,jk) = za*pt(jl,jk) + zb*pt(jl,ikp) + pgrd(jl,ikp)

    END DO
  END DO

  RETURN
END SUBROUTINE pgrad
