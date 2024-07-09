!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE conteq(pdt,pdlnps,pdiv,kbdim,klen,lpperm)

  ! Description:
  !
  ! Computes temperature and surface pressure increments for semi-implict
  ! scheme.
  !
  ! Method:
  !
  ! To calculate linearized increments of the temperature
  ! and the logarithm of surface pressure for a specified divergence.
  !
  ! *conteq* is called from *inhysi* in the initial calculation of the 
  !          gravity-wave matrix *bb,* and from subroutines and in the 
  !          calculation of semi-implicit corrections. 
  !
  ! Parameters are:
  !
  !   *pdt*       *computed temperature increment.
  !   *kbdim*     *first dimension of 2-d array.
  !   *pdlnps*    *computed increment of ln(surface pressure).
  !   *pdiv*      *specified divergence.
  !   *klen*      *number of points for which calculation is
  !                performed.
  !   *lpperm*    *logical switch to permute the order of the
  !                indexes in the 2 dimensional arrays (if .t.).
  !
  ! Required constants are obtained from module *mo_hyb*.
  ! 
  ! Results are computed for *klen* consecutive points or
  ! spectral coefficients and (for temperature) at each model level.
  !
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

  USE mo_kind,       ONLY: dp
  USE mo_control,    ONLY: nlev
  USE mo_hyb,        ONLY: aktlrd, altrcp, delpr, nlevm1, rpr

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kbdim, klen
  LOGICAL :: lpperm

  !  Array arguments 
  REAL(dp) :: pdiv(kbdim,*), pdlnps(*), pdt(kbdim,*)

  !  Local scalars: 
  REAL(dp) :: za, zdp, zl
  INTEGER :: ikp, jk, jl

  !  Executable statements 

!-- 1. Integrate continuity equation

  IF (lpperm) THEN
    DO jl = 1, klen
      pdt(1,jl) = 0._dp
    END DO
  ELSE
    DO jl = 1, klen
      pdt(jl,1) = 0._dp
    END DO
  END IF

  DO jk = 1, nlevm1
    ikp = jk + 1
    zdp = delpr(jk)
    zl = aktlrd(jk)
    za = altrcp(jk)

    IF (lpperm) THEN
      DO jl = 1, klen
        pdt(ikp,jl) = zdp*pdiv(jk,jl) + pdt(jk,jl)
        pdt(jk,jl) = za*pdiv(jk,jl) + zl*pdt(jk,jl)
      END DO
    ELSE
      DO jl = 1, klen
        pdt(jl,ikp) = zdp*pdiv(jl,jk) + pdt(jl,jk)
        pdt(jl,jk) = za*pdiv(jl,jk) + zl*pdt(jl,jk)
      END DO

    END IF
  END DO

  zdp = delpr(nlev)
  zl = aktlrd(nlev)
  za = altrcp(nlev)
  IF (lpperm) THEN
    DO jl = 1, klen
      pdlnps(jl) = -rpr*(zdp*pdiv(nlev,jl)+pdt(nlev,jl))
      pdt(nlev,jl) = za*pdiv(nlev,jl) + zl*pdt(nlev,jl)
    END DO
  ELSE
    DO jl = 1, klen
      pdlnps(jl) = -rpr*(zdp*pdiv(jl,nlev)+pdt(jl,nlev))
      pdt(jl,nlev) = za*pdiv(jl,nlev) + zl*pdt(jl,nlev)
    END DO
  END IF

  RETURN
END SUBROUTINE conteq
