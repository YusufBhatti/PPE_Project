#ifdef __xlC__
@PROCESS HOT
@PROCESS ALIAS(NOARYOVRLP,NOPTEOVRLP)
#endif

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE tf1

  ! Description:
  !
  ! It performs the first part of time filtering for next time step:
  ! xf=x+eps*(xm1-2*x).
  !
  ! Method:
  !
  ! *tf1* is called from *scan1sl* (subscan 2).
  !
  ! Reference:
  ! 1-appendix *b1:organisation of the spectral model.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, May 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! U. Schlese, DKRZ, Novenber 1999, cloud ice added
  ! I. Kirchner, MPI, December 2000, time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! D. Klocke, MPI, Nov, 2010, lnwp for NWP restarts
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_control,       ONLY: ltimer, lnwp
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_tf1
  USE mo_time_control,  ONLY: lstart
  USE mo_tracdef,       ONLY: trlist
  USE mo_memory_gl,     ONLY: q, xl, xi, xt
  USE mo_memory_g1a,    ONLY: alpsm1, dm1, qm1, tm1, vom1, xlm1, xim1, &
                              xtm1
  USE mo_memory_g1b,    ONLY: alpsf, df, qf, tf, vof, xlf, xif,        &
                              xtf
  USE mo_memory_g2a,    ONLY: um1, vm1, dudlm1, dvdlm1
  USE mo_memory_g2b,    ONLY: uf, vf, dudlf, dvdlf
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_semi_impl,     ONLY: eps
  USE mo_scan_buffer,   ONLY: vo, d, alps, t, u, v, dudl, dvdl

  IMPLICIT NONE

  !  Local scalars: 
  REAL(dp):: z1m2eps, zeps
  INTEGER :: jk, jl, jt, jrow
  INTEGER :: nlev, ngpblks, nproma


  !  Executable statements 

  nlev    = ldc% nlev    ! global: number of levels
  ngpblks = ldc% ngpblks ! number of rows

  !-- 1. First part of time filtering

  IF (lstart .OR. lnwp) THEN
    zeps = 0._dp
    z1m2eps = 1._dp
  ELSE
    zeps = eps
    z1m2eps = 1._dp - 2._dp*eps
  END IF

!CSD$ PARALLEL DO PRIVATE(nproma,jk,jl,jt)
!$OMP PARALLEL PRIVATE(nproma,jk,jl,jt,jrow)
  IF (ltimer) CALL timer_start(timer_tf1)
!$OMP DO
  DO jrow = 1, ngpblks

    IF ( jrow == ngpblks ) THEN
      nproma = ldc% npromz
    ELSE
      nproma = ldc% nproma
    END IF

    DO jk = 1, nlev
!DIR$ IVDEP
!DIR$ CONCURRENT
      DO jl = 1, nproma
        vof(jl,jk,jrow) = z1m2eps*vo(jl,jk,jrow) + zeps*vom1(jl,jk,jrow)
        df (jl,jk,jrow) = z1m2eps*d(jl,jk,jrow) + zeps*dm1(jl,jk,jrow)
        qf (jl,jk,jrow) = z1m2eps*q(jl,jk,jrow)  + zeps*qm1 (jl,jk,jrow)
        xlf(jl,jk,jrow) = z1m2eps*xl(jl,jk,jrow) + zeps*xlm1(jl,jk,jrow)
        xif(jl,jk,jrow) = z1m2eps*xi(jl,jk,jrow) + zeps*xim1(jl,jk,jrow)
        tf (jl,jk,jrow) = z1m2eps*t(jl,jk,jrow) + zeps*tm1(jl,jk,jrow)
        uf (jl,jk,jrow) = z1m2eps*u(jl,jk,jrow) + zeps*um1(jl,jk,jrow)
        vf (jl,jk,jrow) = z1m2eps*v(jl,jk,jrow) + zeps*vm1(jl,jk,jrow)
        dudlf(jl,jk,jrow) = z1m2eps*dudl(jl,jk,jrow) + zeps*dudlm1(jl,jk,jrow)
        dvdlf(jl,jk,jrow) = z1m2eps*dvdl(jl,jk,jrow) + zeps*dvdlm1(jl,jk,jrow)
      END DO
    END DO

    DO jt = 1, trlist% ntrac
      IF (trlist% ti(jt)% nint == 1) THEN
        DO jk = 1, nlev
!DIR$ IVDEP
!DIR$ CONCURRENT
          DO jl = 1, nproma
            xtf(jl,jk,jt,jrow) = z1m2eps*xt(jl,jk,jt,jrow) + zeps*xtm1(jl,jk,jt,jrow)
          END DO
        END DO
      ELSE
!DIR$ CONCURRENT
        xtf(:,:,jt,jrow) = xtm1(:,:,jt,jrow)
      ENDIF
    END DO

!DIR$ IVDEP
!DIR$ CONCURRENT
    DO jl = 1, nproma
      alpsf(jl,jrow) = z1m2eps*alps(jl,jrow) + zeps*alpsm1(jl,jrow)
    END DO

  END DO
!$OMP END DO NOWAIT
  IF (ltimer) CALL timer_stop(timer_tf1)
!$OMP END PARALLEL
!CSD$ END PARALLEL DO 

END SUBROUTINE tf1
