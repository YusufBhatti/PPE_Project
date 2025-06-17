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
SUBROUTINE legtri(psin,kcp,palp)

  ! Description:
  !
  ! Legendre functions for a triangular truncation.
  !
  ! Method:
  !
  ! This routine computes the values *palp* for the argument
  ! *psin* of the normalised *Legendre associated functions in the
  ! order ((jn1=jm1,kcp),jm1=1,kcp) for jn=jn1-1 and jm=jm1-1 .
  !
  ! *legtri* is called from *radmod*.
  !
  ! There are three dummy arguments:
  !          *psin* is the sine of latitude.
  !          *kcp* is one plus the limit wave number.
  !          *palp* is the array of the results.
  !
  ! Simple recurence formula.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, June 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !
  USE mo_kind,       ONLY: dp

  IMPLICIT NONE

  !  Scalar arguments 
  REAL(dp):: psin
  INTEGER :: kcp

  !  Array arguments 
  REAL(dp):: palp(*)

  !  Local scalars: 
  REAL(dp):: z2m, zcos, ze1, ze2, zf1m, zf2m, zm, zn, zn2, zre1, zsin
  INTEGER :: ic, icp, jj, jm, jm1, jm2, jn

  !  Intrinsic functions 
  INTRINSIC SQRT


  !  Executable statements 

!-- 1. Preliminary setting

  zsin = psin
  icp = kcp

!-- 2. Computations

  ic = icp - 1
  zcos = SQRT(1._dp-zsin**2)
  jj = 2
  palp(1) = 1._dp
  zf1m = SQRT(3._dp)
  palp(2) = zf1m*zsin
  DO jm1 = 1, icp
    jm = jm1 - 1
    zm = jm
    z2m = zm + zm
    zre1 = SQRT(z2m+3._dp)
    ze1 = 1._dp/zre1
    IF (jm/=0) THEN
      zf2m = zf1m*zcos/SQRT(z2m)
      zf1m = zf2m*zre1
      jj = jj + 1
      palp(jj) = zf2m
      IF (jm==ic) CYCLE
      jj = jj + 1
      palp(jj) = zf1m*zsin
      IF (jm1==ic) CYCLE
    END IF
    jm2 = jm + 2
    DO jn = jm2, ic
      zn = jn
      zn2 = zn**2
      ze2 = SQRT((4._dp*zn2-1._dp)/(zn2-zm**2))
      jj = jj + 1
      palp(jj) = ze2*(zsin*palp(jj-1)-ze1*palp(jj-2))
      ze1 = 1._dp/ze2
    END DO

  END DO

  RETURN
END SUBROUTINE legtri
