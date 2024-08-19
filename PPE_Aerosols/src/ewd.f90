#ifdef __xlC__
@PROCESS STRICT
@PROCESS ALIAS(NOARYOVRLP,NOPTEOVRLP)
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE ewd

  ! Description:
  !
  ! Computes east west derivatives.
  !
  ! Method:
  !
  ! This subroutine computes east-west derivatives
  ! divided by(a*(1.-mu**2)).
  !
  ! *ewd* is called from *scan1sl*.
  !
  ! The east- west derivation is reduced in *fourier space to
  ! a multiplication by i*m.
  !
  ! 1-appendix *b1:organisation of the spectral model.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp 
  USE mo_buffer_fft,    ONLY: falps, ft, fdalpsl, fdtl, fu, fv, fdudl, fdvdl
  USE mo_gaussgrid,     ONLY: gl_racst
  USE mo_truncation,    ONLY: am
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_control,       ONLY: ltimer
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_ewd

  IMPLICIT NONE

  !  Local scalars: 
  REAL(dp):: zmracst, zracst
  INTEGER :: jglat, jlev, jm, jlat
  INTEGER :: n2mp1, nlev, nlevp1, nmp1, nlat
  INTEGER, POINTER :: glat (:)

  !  Executable statements 

!-- array bounds on local PE (Fourier space)

  nlev   =  lc% nflev
  nlevp1 =  lc% nflevp1
  nmp1   =  lc% nm + 1
  n2mp1  = (lc% nm + 1) * 2
  nlat   =  lc% nflat
  glat   => lc% glat

!-- 1. Compute east-west derivatives and blank end of the fields.

!$OMP PARALLEL 

  IF (ltimer) CALL timer_start(timer_ewd)

!$OMP DO PRIVATE(jlat, jglat, zracst, jlev, jm, zmracst)
  DO jlat = 1, nlat

    jglat = glat(jlat)

    zracst = gl_racst(jglat)

    DO jlev = 1, nlev

!DIR$ IVDEP
!OCL NOVREC
!IBM* ASSERT(NODEPS)
      DO jm = 1, nmp1
        zmracst = am(jm)*zracst
        fdtl (2*jm-1,jlev,jlat) = -zmracst*ft(2*jm  ,jlev,jlat)
        fdtl (2*jm  ,jlev,jlat) =  zmracst*ft(2*jm-1,jlev,jlat)
        fdudl(2*jm-1,jlev,jlat) = -am(jm) *fu(2*jm  ,jlev,jlat)
        fdudl(2*jm  ,jlev,jlat) =  am(jm) *fu(2*jm-1,jlev,jlat)
        fdvdl(2*jm-1,jlev,jlat) = -am(jm) *fv(2*jm  ,jlev,jlat)
        fdvdl(2*jm  ,jlev,jlat) =  am(jm) *fv(2*jm-1,jlev,jlat)
      END DO

      fdtl (n2mp1+1:,jlev,jlat) = 0.0_dp
      fdudl(n2mp1+1:,jlev,jlat) = 0.0_dp
      fdvdl(n2mp1+1:,jlev,jlat) = 0.0_dp

    END DO

    IF (nlevp1>nlev) THEN
!DIR$ IVDEP
!OCL NOVREC
!IBM* ASSERT(NODEPS)
      DO jm = 1, nmp1
        zmracst = am(jm)*zracst
        fdalpsl(2*jm-1,jlat) = -zmracst*falps(2*jm  ,jlat)
        fdalpsl(2*jm  ,jlat) =  zmracst*falps(2*jm-1,jlat)
      END DO
      fdalpsl(n2mp1+1:,jlat) = 0.0_dp
    END IF

  END DO
!$OMP END DO

  IF (ltimer) CALL timer_stop(timer_ewd)

!$OMP END PARALLEL

END SUBROUTINE ewd
