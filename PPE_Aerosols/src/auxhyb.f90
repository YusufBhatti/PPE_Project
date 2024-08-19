#ifdef __xlC__
@PROCESS HOT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE auxhyb(pdelp,prdelp,plnpr,palpha,ph,kdim,klen)

  ! Description:
  !
  ! calculates auxiliary variables connected with the vertical 
  ! finite-difference scheme.
  !
  ! Method:
  !
  ! *To compute full-level values of auxiliary variables
  ! which depend on pre-computed half-level values of pressure.
  ! 
  ! *auxhyb* is called from several points within the forecasting system 
  ! in connection with the vertical discretization.
  ! Parameters are:
  !
  !   *pdelp*       *computed pressure difference across layers.
  !   *prdelp*      *reciprocal of *pdelp.*
  !   *plnpr*       *computed logarithm of ratio of pressures.
  !   *palpha*      *computed alphas for integration of the
  !                  hydrostatic equation and related terms.
  !   *ph*          *specified half-level pressures.
  !   *kdim*        *first dimension of 2-d arrays *pdelp,*
  !                  *plnpr,* *palpha,* and *ph.*
  !   *klen*        *number of points for which calculation is performed.
  !
  ! Required constants are obtained from modules
  ! *mo_physical_constants* and *mo_hyb*. The latter should have been
  ! initialized by a call of subroutine *inihyb*.
  ! 
  ! Results are computed for *klen* consecutive points for
  ! all required levels.
  ! 
  ! Calculations are performed separately for pressure,
  ! hybrid and sigma levels.
  !
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

  USE mo_kind,      ONLY: dp
  USE mo_physical_constants, ONLY: rd
  USE mo_control,   ONLY: nlev
  USE mo_hyb,       ONLY: delpr, nlmsgl, nlmsla, nplev, nplvpa,        &
                          ralpha, rdelpr, rlnpr

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER ,INTENT(in)  :: kdim, klen

  !  Array arguments 
  REAL(dp) ,INTENT(in)  :: ph(kdim, *)
  REAL(dp) ,INTENT(inout) :: palpha(kdim, *), pdelp(kdim, *),          &
                             plnpr(kdim, *),  prdelp(kdim, *)

  !  Local scalars: 
  REAL(dp):: za, zd, zl, zr, zrd
  INTEGER :: jk, jl

  !  Intrinsic functions 
  INTRINSIC LOG


  !  Executable statements 

!-- 1. Initialize variables

  zrd = rd

!-- 2. Set palpha and plnpr for top level

  palpha(1:klen,1) = ralpha(1)
  plnpr (1:klen,1) = rlnpr(1)

!-- 3. Set pressure-level values or other top-level values

  IF (nplev==0) THEN

!IBM* NOVECTOR
    DO jl = 1, klen
      pdelp(jl,1) = ph(jl,2) - ph(jl,1)
      prdelp(jl,1) = 1._dp/pdelp(jl,1)
    END DO

  ELSE

    DO jl = 1, klen
      pdelp(jl,1) = delpr(1)
      prdelp(jl,1) = rdelpr(1)
    END DO

    DO jk = 2, nplev
      zd = delpr(jk)
      zr = rdelpr(jk)
      zl = rlnpr(jk)
      za = ralpha(jk)

!IBM* NOVECTOR
      DO jl = 1, klen
        pdelp(jl,jk) = zd
        prdelp(jl,jk) = zr
        plnpr(jl,jk) = zl
        palpha(jl,jk) = za

      END DO
    END DO

  END IF

!-- 4. Calculate hybrid-level values

  DO jk = nplvpa, nlmsgl

    DO jl = 1, klen
      pdelp(jl,jk) = ph(jl,jk+1) - ph(jl,jk)
      prdelp(jl,jk) = 1._dp/pdelp(jl,jk)
      plnpr(jl,jk) = zrd*LOG(ph(jl,jk+1)/ph(jl,jk))
      palpha(jl,jk) = zrd - ph(jl,jk)*plnpr(jl,jk)*prdelp(jl,jk)

    END DO

  END DO

!-- 5. Set sigma-level values

  DO jk = nlmsla, nlev
    zl = rlnpr(jk)
    za = ralpha(jk)
    DO jl = 1, klen
      pdelp(jl,jk) = ph(jl,jk+1) - ph(jl,jk)
      prdelp(jl,jk) = 1._dp/pdelp(jl,jk)
      plnpr(jl,jk) = zl
      palpha(jl,jk) = za
    END DO
  END DO

  RETURN
END SUBROUTINE auxhyb
