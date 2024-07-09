#ifdef __xlC__
@PROCESS STRICT
#endif
!OCL NOALIAS

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE scctp

  ! Description:
  !
  ! Adds the implicit contribution of divergence to the temperature
  ! and surface pressure equations.
  !
  ! Method:
  !
  ! This subroutine adds the implicit contribution of divergence to the
  ! temperature and surface pressure equations.
  ! (see appendix *b-1,par 4-5-4 and 4-5-5)
  !
  ! *scctp* is called from *stepon*
  ! after the direct *legendre transforms of *t* and *p*
  !
  ! Results:
  ! The implicit contribution of divergence is provided by
  ! *conteq* in temporary arrays *zt* and *zp*
  ! *conteq*  now inlined !
  !
  ! Reference:
  ! 1-appendix b1:organisation of the spectral model.
  ! m.j       5/10/81.
  ! 2-sequence of calculations for the vertical hybrid scheme.
  ! a.s.     20/11/81.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, March 1982, original source
  ! U. Schlese, DKRZ, May 1993, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! T. Diehl, DKRZ, July 1999, parallel version 
  !
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_memory_sp,     ONLY: stp, sd
  USE mo_control,       ONLY: nlev, nlevp1
  USE mo_hyb,           ONLY: aktlrd, altrcp, delpr, nlevm1, rpr

  IMPLICIT NONE

  !  Local scalars: 
  REAL(dp):: za, zdp, zl
  INTEGER :: ikp, jj, jk, jlev, jt, snsp

  !  Local arrays: 
  REAL(dp):: zp(2,lc%snsp), zt(2,lc%snsp,nlev)


  !  Executable statements 

  snsp = lc%snsp

!-- 1. Compute temperature and surface pressure increments
!      and add them to the corresponding spectral fields.

  DO jj = 1, 2
    DO jt = 1, snsp
      zt(jj,jt,1) = 0._dp
    END DO
  END DO

  DO jk = 1, nlevm1
    ikp = jk + 1
    zdp = delpr(jk)
    zl = aktlrd(jk)
    za = altrcp(jk)
    DO jj = 1, 2
      DO jt = 1, snsp
        zt(jj,jt,ikp) = zdp*sd(jk,jj,jt) + zt(jj,jt,jk)
        zt(jj,jt,jk) = za*sd(jk,jj,jt) + zl*zt(jj,jt,jk)
      END DO
    END DO
  END DO

  zdp = delpr(nlev)
  zl = aktlrd(nlev)
  za = altrcp(nlev)
  DO jj = 1, 2
!DIR$ IVDEP
    DO jt = 1, snsp
      zp(jj,jt) = -rpr*(zdp*sd(nlev,jj,jt)+zt(jj,jt,nlev))
      stp(nlevp1,jj,jt) = stp(nlevp1,jj,jt) + zp(jj,jt)
      zt(jj,jt,nlev) = za*sd(nlev,jj,jt) + zl*zt(jj,jt,nlev)
    END DO
  END DO

  DO jlev = 1, nlev
    DO jj = 1, 2
      DO jt = 1, snsp
          stp(jlev,jj,jt) = stp(jlev,jj,jt) + zt(jj,jt,jlev)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE scctp
