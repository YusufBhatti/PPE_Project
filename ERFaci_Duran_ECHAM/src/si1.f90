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
SUBROUTINE si1_extended_const

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: ldc => local_decomposition
  USE mo_semi_impl,     ONLY: betazq
  USE mo_control,       ONLY: nlev
  USE mo_time_control,  ONLY: time_step_len
  USE mo_scan_buffer,   ONLY: u0
  USE mo_geoloc,        ONLY: rz1u0
  USE mo_gaussgrid,     ONLY: gl_racst
  USE mo_transpose,     ONLY: reorder

  IMPLICIT NONE

  INTEGER :: jlev
  INTEGER :: jlat, lnlat, jglat
  REAL(dp):: zdt
  REAL(dp):: zrz1u0(ldc% nglon, nlev, ldc% nglat)

  zdt = 0.5_dp*time_step_len

  lnlat = ldc%nglat

  DO jlat = 1, lnlat
    jglat = ldc% glat(jlat) ! global latitude index N -> S
    DO jlev = 1, nlev
      zrz1u0(:,jlev,jlat) = (betazq*zdt*gl_racst(jglat))*u0(jlev,jlat)
    END DO
  END DO

  CALL reorder (rz1u0, zrz1u0)

END SUBROUTINE si1_extended_const

!OCL NOALIAS

SUBROUTINE si1(krow)

  ! Description:
  !
  ! 1st part of the semi-implicit scheme (done in grid point space). 
  !
  ! Method:
  !
  ! This subroutine computes the contribution in
  ! grid points to the semi-implicit scheme.
  !
  ! *si1* is called from *gpc*.
  !
  ! Reference:
  ! 1-appendix *b1:organisation of the spectral model.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! I. Kirchner, MPI, November 2000, date/time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! I. Kirchner, MPI, August 2002, bugfix tendency diagnostics
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
                              
  USE mo_scan_buffer,   ONLY: vol, vom,  qte, rh, xlte, xite, tte, xtte, &
                              alps, alpste, d, t, vo, dm
  USE mo_memory_g1a,    ONLY: vom1, qm1, xlm1, xim1, xtm1, dm1, tm1,   &
                              alpsm1
  USE mo_geoloc,        ONLY: rz1u0, racst_2d
  USE mo_control,       ONLY: nlev, ltimer
  USE mo_time_control,  ONLY: time_step_len
  USE mo_physical_constants,     ONLY: earth_radius
  USE mo_tracdef,       ONLY: trlist
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_si1

  IMPLICIT NONE

  INTEGER :: krow

  !  Local scalars: 
  REAL(dp) :: z1u0, z2, z3, zdt, ztwodt
  INTEGER :: jlev, jl, jrow, jt
  INTEGER :: nproma, nbdim
  LOGICAL :: loperm

  !  Local arrays: 
  REAL(dp), TARGET  :: zd(ldc% nproma, nlev)
  REAL(dp)          :: zp(ldc% nproma), zt(ldc% nproma, nlev)
  REAL(dp), POINTER :: zr(:,:)

  !  External subroutines 

  EXTERNAL :: conteq, pgrad

  !  Executable statements 

  IF (ltimer) CALL timer_start(timer_si1)

  !-- 1. Locate and allocate storage

  jrow  = krow        ! local  continuous index

  nbdim = ldc% nproma

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

  !   zp(:)   = 0.0_dp
  !   zt(:,:) = 0.0_dp
  !   zd(:,:) = 0.0_dp

  !-- 1.2 Equivalence arrays *zd* and *zr*

  zr => zd

  !-- 2. Skip over *si1* during initialisation iterations
  !      or prepare some temporary constants.

  !-- 2.1 Compute temporary constants

  zdt    = 0.5_dp*time_step_len
  ztwodt =        time_step_len

  z2 = ztwodt/earth_radius

  !-- 3. Semi implicit computations

  !-- 3.1 Vorticity and humidity equations
  !       and term dmu of divergence equation.
  DO jlev = 1, nlev

     ! get explicit part of time schema, (2*dt) included
     ! explicit part of vorticity

!DIR$ IVDEP
     DO jl = 1, nproma

        z1u0 = rz1u0(jl,jlev,jrow)
        z3   = ztwodt*racst_2d(jl,jrow)

        dm(jl,jlev,jrow)  =  z2 * vol (jl,jlev,jrow)

        vol(jl,jlev,jrow) =  z3 * vol (jl,jlev,jrow)                  &
                         -z1u0*(vom1(jl,jlev,jrow) - 2._dp*vo(jl,jlev,jrow))
        vom(jl,jlev,jrow) = -z2  * vom (jl,jlev,jrow)

        qm1(jl,jlev,jrow)  = qm1 (jl,jlev,jrow) + ztwodt*qte (jl,jlev,jrow)
        xlm1(jl,jlev,jrow) = xlm1(jl,jlev,jrow) + ztwodt*xlte(jl,jlev,jrow)
        xim1(jl,jlev,jrow) = xim1(jl,jlev,jrow) + ztwodt*xite(jl,jlev,jrow)
     END DO
  END DO

  DO jt = 1, trlist% ntrac
     if (trlist% ti(jt)% nint /= 1) CYCLE
     DO jlev = 1, nlev
        DO jl = 1, nproma
           xtm1(jl,jlev,jt,jrow) = xtm1(jl,jlev,jt,jrow) +             &
                                                ztwodt*xtte(jl,jlev,jt,jrow)
        END DO
     END DO
  END DO

  !-- 3.2 Compute implicit contribution of divergence
  !       to temperature and surface equations.

  DO jlev = 1, nlev
     DO jl = 1, nproma
        zd(jl,jlev) = .5_dp*dm1(jl,jlev,jrow) - d(jl,jlev,jrow)
     END DO
  END DO
  loperm = .FALSE.
  CALL conteq(zt,zp,zd,nbdim,nproma,loperm)

  !-- 3.3 Update *zt* and *zp* to compute the contribution
  !       of temperature and surface pressure to the
  !       divergence equation.

  DO jlev = 1, nlev
     DO jl = 1, nproma
        zt(jl,jlev) = zt(jl,jlev) + tm1(jl,jlev,jrow) - t(jl,jlev,jrow) +   &
                      zdt*tte(jl,jlev,jrow)
     END DO
  END DO
  DO jl = 1, nproma
     zp(jl) = zp(jl) + alpsm1(jl,jrow) - alps(jl,jrow) + zdt*alpste(jl,jrow)
  END DO
  CALL pgrad(zr,zt,nbdim,zp,nproma)

  ! explicit part of divergence
!DIR$ CONCURRENT

  !-- 3.4 Complete computation of the terms to be
  !       passed to *fftd*.

  DO jlev = 1, nlev
!DIR$ CONCURRENT
!DIR$ IVDEP
     DO jl = 1, nproma
        tm1(jl,jlev,jrow) = 2._dp*zt(jl,jlev) - tm1(jl,jlev,jrow) +    &
                                               2._dp*t(jl,jlev,jrow)
        rh(jl,jlev,jrow) = -ztwodt*(rh(jl,jlev,jrow)+zr(jl,jlev))
     END DO
  END DO

  DO jl = 1, nproma
     alpsm1(jl,jrow) = 2._dp*zp(jl) - alpsm1(jl,jrow) + 2._dp*alps(jl,jrow)
  END DO

  IF (ltimer) CALL timer_stop(timer_si1)

END SUBROUTINE si1
