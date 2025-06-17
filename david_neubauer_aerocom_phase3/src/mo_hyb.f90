!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_hyb

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! ------------------------------------------------------------------
  !
  ! module *mo_hyb* - *loop indices and surface-pressure independent
  !        variables associated with the vertical finite-difference scheme.
  !
  ! a.j.simmons     e.c.m.w.f.     16/11/81.
  !
  ! ------------------------------------------------------------------

  INTEGER :: nlevm1         !  (number of levels)-1.
  INTEGER :: nplev          !  *number of pressure levels.
  INTEGER :: nplvp1         !  *nplev+1.*
  INTEGER :: nplvp2         !  *nplev+2.*
  INTEGER :: nplvpa         !  *nplvp1,* or 2 if *nplev=0.*
  INTEGER :: nlmsgl         !  *nlev* - (number of sigma levels).
  INTEGER :: nlmslp         !  *nlmsgl+1.*
  INTEGER :: nlmsla         !  *nlmslp,* or 2 if *nlmslp=1.*

  REAL(dp) :: apzero            !  *reference pressure for computation of the
                            !   hybrid vertical levels.
  REAL(dp) :: rpr               !  *reciprocal of reference surface pressure.
  REAL(dp) :: rdtr              !   rd*(reference temperature).
  REAL(dp) :: apsurf            !   fixed global mean of surface pressure
  REAL(dp) :: t0icao            !  *surface temperatur of reference atmosphere
  REAL(dp) :: tsticao           !  *stratospheric temperature of reference atmosphere
  REAL(dp) :: rdtstic           !  *rd*tsticao
  REAL(dp) :: rdlnp0i           !  *rd*ln(surface pressure) of reference atmosphere
  REAL(dp) :: alrrdic           !  *lapse-rate parameter of reference atmosphere
  REAL(dp) :: rdt0ral           !  *rd*t0icao/alphaic
  REAL(dp) :: ptricao           !  *tropopause pressure of reference atmosphere
  REAL(dp) :: rdlnpti           !  *rd*ln(ptricao)
  REAL(dp) :: gsticao           !  *constant used in geopotential calculation
  REAL(dp), ALLOCATABLE :: ralpha(:)    !   rd*alpha at pressure and sigma levels.
  REAL(dp), ALLOCATABLE :: rlnpr(:)     !   rd*ln(p(k+.5)/p(k-.5)) at pressure and sigma levels.
  REAL(dp), ALLOCATABLE :: dela(:)      !   a(k+.5)-a(k-.5).
  REAL(dp), ALLOCATABLE :: delb(:)      !   b(k+.5)-b(k-.5).
  REAL(dp), ALLOCATABLE :: rddelb(:)    !   rd*delb.
  REAL(dp), ALLOCATABLE :: cpg(:)       !   a(k+.5)*b(k-.5)-b(k+.5)*a(k-.5).
  REAL(dp), ALLOCATABLE :: delpr(:)     !   p(k+.5)-p(k-.5) for reference surface pressure.
  REAL(dp), ALLOCATABLE :: rdelpr(:)    !  *reciprocal of *delpr.*
  REAL(dp), ALLOCATABLE :: ralphr(:)    !  *constant array for use by pgrad.
  REAL(dp), ALLOCATABLE :: alpham(:)    !  *constant array for use by dyn.
  REAL(dp), ALLOCATABLE :: ardprc(:)    !  *constant array for use by dyn.
  REAL(dp), ALLOCATABLE :: rlnmar(:)    !  *constant array for use by pgrad.
  REAL(dp), ALLOCATABLE :: aktlrd(:)    !  *constant array for use by conteq.
  REAL(dp), ALLOCATABLE :: altrcp(:)    !  *constant array for use by conteq.
  REAL(dp), ALLOCATABLE :: ceta(:)      !  *full hybrid vertical levels.
  REAL(dp), ALLOCATABLE :: cetah(:)    !  *half hybrid vertical levels.
  REAL(dp), ALLOCATABLE :: bb(:,:) !  *gravity wave matrix

CONTAINS

  !+ initializes constants for vertical coordinate calculations.
  !+ $Id$

  SUBROUTINE inihyb

    ! Description:
    !
    ! Initializes constants for vertical coordinate calculations.
    !
    ! Method:
    !
    ! Compute loop indices and surface-pressure independent
    ! variables associated with the vertical finite-difference scheme.
    !
    ! Output is in module *mo_hyb*
    !
    ! Authors:
    !
    ! A. J. Simmons, ECMWF, November 1981, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, Jan 1999, subroutine inihyb -> module mo_hyb
    ! 
    ! for more details see file AUTHORS
    !
    
    USE mo_physical_constants, ONLY: grav, rd, rcpd
    USE mo_control,   ONLY: nlev, nlevp1, nvclev, vct

    !  Local scalars: 
    REAL(dp) :: za, zb, zetam, zetap, zp, zp0icao, zpp, zrd, zs, zsm
    INTEGER :: ilev, ilevp1, iplev, iplvp1, is, ism, ist, jk, jlev

    !  Intrinsic functions 
    INTRINSIC EXP, LOG


    !  Executable statements 

    !-- 1. Initialize variables
    
    apzero = 101325.0_dp
    zrd = rd
    ralpha(1) = zrd*LOG(2.0_dp)
    rlnpr(1) = 2.0_dp*ralpha(1)
    ilev = nlev
    ilevp1 = ilev + 1
    nlevp1 = ilevp1
    nlevm1 = ilev - 1
    iplev = 0
    iplvp1 = 1
    is = nvclev + ilevp1
    ism = is - 1
    zpp = vct(1)
    zsm = vct(is)
    apsurf = 98550.0_dp

    t0icao = 288.0_dp
    tsticao = 216.5_dp
    zp0icao = 101320.0_dp
    rdlnp0i = rd*LOG(zp0icao)
    rdtstic = rd*tsticao
    alrrdic = 0.0065_dp/grav
    rdt0ral = t0icao/alrrdic
    rdlnpti = rdlnp0i + (LOG(tsticao/t0icao))/alrrdic
    ptricao = EXP(rdlnpti/rd)
    gsticao = tsticao*(rdlnpti-1.0_dp/alrrdic)

    !-- 2. Calculate pressure-level values

10  CONTINUE

    zb = vct(nvclev+iplvp1+1)
    IF (zb > 0.0_dp) THEN
      nplev = iplev
      nplvp1 = iplvp1
      nplvp2 = iplvp1 + 1
      IF (iplev==0) THEN
        nplvpa = 2
      ELSE
        nplvpa = iplvp1
      END IF
      GO TO 20
    ELSE
      iplev = iplvp1
      iplvp1 = iplev + 1
      IF (iplvp1==ilevp1) GO TO 40
      zp = zpp
      zpp = vct(iplvp1)
      delpr(iplev) = zpp - zp
      rdelpr(iplev) = 1.0_dp/delpr(iplev)
      IF (iplev>1) THEN
        rlnpr(iplev) = zrd*LOG(zpp/zp)
        ralpha(iplev) = zrd - zp*rlnpr(iplev)/delpr(iplev)
      END IF
      alpham(iplev) = ralpha(iplev)*rcpd
      ardprc(iplev) = rlnpr(iplev)*rdelpr(iplev)*rcpd
      GO TO 10
    END IF

    !-- 3. Calculate sigma-level values

20  CONTINUE

    za = vct(ism-nvclev)
    IF (za > 0.0_dp) THEN
      nlmsgl = ism - nvclev
      nlmslp = nlmsgl + 1
      nlmsla = nlmslp
      GO TO 30
    ELSE
      is = ism
      ism = is - 1
      ist = is - nvclev
      zs = zsm
      zsm = vct(is)
      IF (ist==1) THEN
        nlmsgl = 0
        nlmslp = 1
        nlmsla = 2
        GO TO 30
      ELSE
        rlnpr(ist) = zrd*LOG(zs/zsm)
        ralpha(ist) = zrd - zsm*rlnpr(ist)/(zs-zsm)
      END IF
      GO TO 20
    END IF

    !-- 4. Calculate dela, delb, rddelb, cpg, and complete alphdb

30  CONTINUE

    DO jk = 1, nlev
      dela(jk) = vct(jk+1) - vct(jk)
      delb(jk) = vct(nvclev+jk+1) - vct(nvclev+jk)
      rddelb(jk) = rd*delb(jk)
      cpg(jk) = vct(nvclev+jk)*vct(jk+1) - vct(nvclev+jk+1)*vct(jk)
    END DO

    DO jk = nlmslp, nlev
      alpham(jk) = ralpha(jk)*delb(jk)
    END DO

    !-- 5. Compute full level values of the hybrid coordinate

40  CONTINUE

    zetam = vct(1)/apzero + vct(nvclev+1)
    cetah(1) = zetam

    DO jlev = 1, nlev
      zetap = vct(jlev+1)/apzero + vct(nvclev+1+jlev)
      ceta(jlev) = (zetam+zetap)*0.5_dp
      cetah(jlev+1) = zetap
      zetam = zetap
    END DO
    
  END SUBROUTINE inihyb

END MODULE mo_hyb
