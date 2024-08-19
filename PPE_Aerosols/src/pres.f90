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
SUBROUTINE pres(ph,kdimp,ps,klen)

  ! Description:
  !
  ! Calculate half-level pressures at all model levels
  ! for a given surface pressure.
  !
  ! Method:
  !
  ! Calculations are performed separately for pressure,
  ! hybrid and sigma levels.
  !
  ! This subroutine is called from many points within the
  ! forecasting system. 
  !
  ! Parameters are:
  !    *ph*        computed half-level pressures.
  !    *kdimp*     first dimension of 2-d array *ph.*
  !    *ps*        surface pressure.
  !    *klen*      number of points for which calculation is
  !                performed.
  !
  ! Required constants are obtained from module *mo_hyb*. 
  ! The latter must have been initialiazed
  ! by a call of subroutine *inihyb*.
  !
  ! Results:
  ! Results are computed for *klen* consecutive points at
  ! each model half level.
  !
  ! Reference:
  ! External documentation of the model equations and the
  ! organization of the vertical calculation.
  !
  ! Authors:
  !
  ! A. J. Simmons, ECMWF, November 1981, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,    ONLY: wp
  USE mo_control, ONLY: nlevp1, nvclev, vct
  USE mo_hyb,     ONLY: nlmsgl, nlmslp, nplvp1, nplvp2

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER ,INTENT(in)  :: kdimp, klen

  !  Array arguments 
  REAL(wp),INTENT(in)  :: ps(*)
  REAL(wp),INTENT(inout) :: ph(kdimp,*)

  !  Local scalars: 
  REAL(wp):: zb, zp
  INTEGER :: jk, jl


  !  Executable statements 

!-- 1. Transfer pressure-level values

  DO jk = 1, nplvp1
    zp = vct(jk)

    DO jl = 1, klen
      ph(jl,jk) = zp

    END DO
  END DO

!-- 2. Compute hybrid-level values

  DO jk = nplvp2, nlmsgl
    zp = vct(jk)
    zb = vct(jk+nvclev)

    DO jl = 1, klen
      ph(jl,jk) = zp + zb*ps(jl)

    END DO
  END DO

!-- 3. Compute sigma-level values

  DO jk = nlmslp, nlevp1
    zb = vct(jk+nvclev)

    DO jl = 1, klen
      ph(jl,jk) = zb*ps(jl)

    END DO
  END DO

  RETURN
END SUBROUTINE pres
