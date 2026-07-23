!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_legendre
  !======================================================
  ! This module gathers all routines and coefficents
  ! required for performing Legendre transformations.
  !======================================================

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: dc => local_decomposition

  IMPLICIT NONE

  !================
  ! Public entities
  !================
  PRIVATE
  !------------------------------------------------------------------
  ! Quantities depending on the current truncation/spatial resolution
  !   set by subroutine inileg
  !------------------------------------------------------------------
  PUBLIC :: pnm     ! Legendre coefficients (scalar transform)
  PUBLIC :: anm     ! Legendre coefficients (?)
  !------------------------------------------------
  ! Quantities depending on the truncation/spatial,
  !   specific to the latitude currently processed  
  ! ..set by legmod
  !------------------------------------------------
  PUBLIC :: pnmd    ! Modified coefficients for direct legendre transform
  PUBLIC :: anmd    ! ...
  PUBLIC :: rnmd    !
  !----------------
  ! ..set by leginv
  !----------------
  PUBLIC :: pnmi    ! Modified coefficients for inverse legendre transform
  PUBLIC :: anmi    ! ...
  PUBLIC :: pnmiuv  !
  PUBLIC :: anmiuv  !
  !------------------
  ! module procedures
  !------------------
  PUBLIC :: inileg  ! Set up polynomials needed for the Legendre transforms.
  PUBLIC :: legmod  ! Calculate modified Legendre polynomials for direct tran.
  PUBLIC :: leginv  ! Calculate modified Legendre polynomials for inverse tra.
  PUBLIC :: cleanup_legendre ! deallocate module variables
  !============================================================================
  !==========================
  ! Declarations of variables
  !==========================
  !------------
  ! coefficents
  !------------
  REAL(dp)          ,ALLOCATABLE :: pnm(:,:) ! Legendre polinominals
  REAL(dp)          ,ALLOCATABLE :: anm(:,:)
  !--------------------------------------------------
  ! Modified legendre coefficients for one latitudude
  !   set by subroutine legmod (module mo_legendre)
  !--------------------------------------------------
  REAL(dp), TARGET,  ALLOCATABLE :: pnmd(:)   ! direct legendre transform
  REAL(dp), TARGET,  ALLOCATABLE :: anmd(:)
  REAL(dp), TARGET,  ALLOCATABLE :: rnmd(:)
  !------------------------------------------------
  !   set by subroutine leginv (module mo_legendre)
  !------------------------------------------------
  REAL(dp), TARGET,  ALLOCATABLE :: pnmi(:)   ! inverse legendre transform
  REAL(dp), TARGET,  ALLOCATABLE :: anmi(:)
  REAL(dp), TARGET,  ALLOCATABLE :: pnmiuv(:)
  REAL(dp), TARGET,  ALLOCATABLE :: anmiuv(:)

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE phcs(pnm, panm, kmp1, kkp1, pmu)

    ! Description:
    !
    ! *phcs* computes the values of the *Legendre polynomials and of their
    ! meridional derivatives at given latitudes for a rhomboidal truncation.
    !
    !
    ! 
    !        ^
    !   *N*  !
    !        !                            .
    !        !            +             .
    !        !          + +           .
    !        !        +   +         .
    !        !      +     +       .
    !        !    +       +     .
    !        !  +         +   .
    ! *KMAX* !+           + .
    !        !            .
    !        !          . .
    !        !        .   .
    !        !      .     .
    !        !    .       .
    !        !  .         .
    !        !________________________________________>
    !                  *MMAX*                      *M*
    !
    !
    ! Method:
    !
    ! *call* *phcs(pnm,panm,kmp1,kkp1,pmu)*
    !
    ! *pnm*    :*legendre polynomials values.
    ! *panm*   :(mu**2-1)*dpnm/dmu.
    ! *kmp1*   :mmax+1.
    ! *kkp1*   :kmax+1.
    ! *pmu*    :value at which *pnm* and *panm* are computed.
    !
    ! The *Legendre polynomials are defined as follows:
    !     * p(n,m)(mu)=sqrt((2n+1)*fact(n-m)/fact(n+m))    *
    !     *            /(fact(n)*2**(n+1))                 *
    !     *            *(1-mu**2)**(m/2)                   *
    !     *            *d**(n+m)(mu**2-1)**n/(dmu)**(n+m)  *
    !
    ! with *fact(k)=k*(k-1)**1
    !
    ! They are computed with the following numerically stable
    ! recurrence relation (Belousov,1962):
    !
    !     * p(n,m)(mu)=c(n,m)*p(n-2,m-2)                   *
    !     *           -d(n,m)*mu*p(n-1,m-2)                *
    !     *           +e(n,m)*p(n-1,m)                     *
    !
    ! with
    !     *c(n,m)=sqrt((2n+1)*(m+n-1)*(m+n-3)              *
    !     *           /(2n-3)/(m+n  )/(m+n-2))             *
    !
    !     *d(n,m)=sqrt((2n+1)*(m+n-1)*(n-m+1)              *
    !     *           /(2n-1)/(m+n  )/(m+n-2))             *
    !
    !     *e(n,m)=sqrt((2n+1)*(n-m)                        *
    !     *           /(2n-1)/(n+m))                       *
    !
    ! The derivatives *(panm)* are then computed as follows:
    !
    !     *pa(n,m)=n*f(n+1,m)*p(n+1,m)                     *
    !     *       -(n+1)*f(n,m)*p(n-1,m)                   *
    !
    ! with:
    !
    !     *f(n,m)=sqrt((n**2-m**2)/(4*n**2-1))             *
    !
    ! Results.:
    ! The *Legendre polynomials and their derivatives are stored
    ! column-wise. The following normalisation is used:
    ! Integral over [-1,+1] of *(pnm**2)* =.5
    !
    ! References:
    ! Belousov,S.L.,1962:Tables of normalised associated Legendre
    ! polynomials.(Mathematical tables series,
    ! Vol 18, PergamonPpress, New York, USA) 379pp
    !

    !  Scalar arguments 
    REAL(dp) :: pmu
    INTEGER :: kkp1, kmp1

    !  Array arguments 
    REAL(dp) :: panm(:), pnm(:)

    !  Local scalars: 
    REAL(dp) :: z2mm1, z2q2, zan, zateta, zcnm, zcos2, zcosfak, zcospar,  &
         zcostet, zdnm, zenm, zmm1, zn, zn2, zn2m1, znm1, zp, zp2, zq,    &
         zq2m1, zsinfak, zsinpar, zsintet, zsqp, zteta, zw, zwm2, zwm2q2, &
         zwnm1, zwq
    INTEGER :: ik, inmax, inmaxm, ito, iton, jk, jm, jn

    !  Local arrays: 
    REAL(dp) :: ztemp(3,kmp1+kkp1)

    !  Intrinsic functions 
    INTRINSIC ACOS, COS, SIN, SQRT


    !  Executable statements 

    ztemp(:,:) = 0.0_dp

    ! 1. Initiate recurrence by computing:
    !      *p(0,0)* *p(1,1)* *pa(0,0)* *pa(1,1)*.

    inmax = kkp1 + kmp1

    zcos2 = SQRT(1.0_dp-pmu**2)
    zteta = ACOS(pmu)
    zan = 1.0_dp

    ztemp(1,1) = 0.5_dp
    DO jn = 2, inmax
      zsinpar = 0.0_dp
      zcospar = 0.0_dp
      zp   = jn - 1.0_dp
      zp2  = zp*zp
      zsqp = 1.0_dp/SQRT(zp2+zp)
      zan  = zan*SQRT(1.0_dp-1.0_dp/(4*zp2))
      zcosfak = 1.0_dp
      zsinfak = zp*zsqp

      ik = jn
      DO jk = 1, ik, 2
        zq = jk - 1.0_dp
        zateta = (zp-zq)*zteta
        zcostet = COS(zateta)
        zsintet = SIN(zateta)
        IF (jn==jk) zcostet = zcostet*0.5_dp
        IF (jk/=1) THEN
           zcosfak = (zq-1.0_dp)/zq*(2*zp-zq+2.0_dp)/(2*zp-zq+1.0_dp)*zcosfak
           zsinfak = zcosfak*(zp-zq)*zsqp
        END IF
        zcospar = zcospar + zcostet*zcosfak
        zsinpar = zsinpar + zsintet*zsinfak
      END DO

      ztemp(1,jn) = zan*zcospar
      ztemp(2,jn-1) = zan*zsinpar
    END DO

    pnm(1) = 0.5_dp
    pnm(1+kkp1) = ztemp(2,1)
    panm(1) = 0.0_dp
    panm(1+kkp1) = pmu*ztemp(2,1)

    ! 2. Complete recurrence

    ! 2.1 First 2 columns.

    DO jn = 2, kkp1
      pnm(jn) = ztemp(1,jn)
      pnm(jn+kkp1) = ztemp(2,jn)
      zn2   = 2._dp*jn
      zn2m1 = zn2 - 1.0_dp
      zn    = REAL(jn,dp)
      znm1  = jn - 1.0_dp
      panm(jn) = znm1*(pmu*ztemp(1,jn)-SQRT(zn2m1/(zn2-3))*ztemp(1,jn-1))
      panm(jn+kkp1) = zn*pmu*ztemp(2,jn) &
                    - SQRT(((zn2+1)*(zn**2-1.0_dp))/zn2m1)*ztemp(2,jn-1)
    END DO

    ! 2.2 Other columns.

    DO jm = 3, kmp1
      zmm1  = jm - 1.0_dp
      z2mm1 = zmm1*2
      ztemp(3,1) = SQRT(1.0_dp+1.0_dp/z2mm1)*zcos2*ztemp(2,1)
      ito = (jm-1)*kkp1
      pnm(ito+1)  = ztemp(3,1)
      panm(ito+1) = zmm1*pmu*ztemp(3,1)
      inmaxm = inmax - jm

      DO jn = 2, inmaxm
        iton   = ito + jn
        znm1   = jn - 1.0_dp
        zq     = z2mm1 + znm1 - 1.0_dp
        zwm2   = zq + znm1
        zw     = zwm2 + 2.0_dp
        zwq    = zw*zq
        zq2m1  = zq*zq - 1.0_dp
        zwm2q2 = zwm2*zq2m1
        zwnm1  = zw*znm1
        z2q2   = zq2m1*2.0_dp
        zcnm   = SQRT((zwq*(zq-2.0_dp))/(zwm2q2-z2q2))
        zdnm   = SQRT((zwq*(znm1+1.0_dp))/zwm2q2)
        zenm   = SQRT(zwnm1/((zq+1.0_dp)*zwm2))
        ztemp(3,jn) = zcnm*ztemp(1,jn) - pmu*(zdnm*ztemp(1,jn+1)-zenm* &
                      ztemp(3,jn-1))
        pnm(iton)  = ztemp(3,jn)
        panm(iton) = (zmm1+znm1)*pmu*ztemp(3,jn) - SQRT(zwnm1*(zq+1.0_dp)/zwm2)* &
                      ztemp(3,jn-1)
      END DO

      DO jn = 1, inmax
        ztemp(1,jn) = ztemp(2,jn)
        ztemp(2,jn) = ztemp(3,jn)
      END DO

    END DO

    RETURN
  END SUBROUTINE phcs

!------------------------------------------------------------------------------
  SUBROUTINE inileg

    ! Description:
    !
    ! Set up polynomials needed for the Legendre transforms.
    !
    ! Method:
    !
    ! *pointer* dimensions made dynamic, and i/o unit numbers changed
    ! to comply with the i/o scheme.
    !
    ! This subroutine computes, normalises and writes suitably
    ! modified *Legendre polynomials needed for the *Legendre
    ! transforms.
    !
    ! *inileg* is called from *control*
    !
    ! Externals:
    !
    ! *phcs*      called to compute the *Legendre polynomials
    !             and their meridional derivatives.
    ! *reord*     reorders the *Legendre polynomials
    !             and their meridional derivatives.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, March 1982, original source
    ! J. K. Gibson, ECMWF, April 82, changed
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, September 1999, changes for parallel version
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,        ONLY: nhgl, nmp1, nnp1
    USE mo_gaussgrid,      ONLY: gl_cst, gl_gmu

    !  Local scalars: 
    REAL(dp) :: zmu, zrcst
    INTEGER :: jgl

    !  Local arrays: 
    REAL(dp) :: zanm(dc%lnsp), zanmt(nmp1*nnp1), zpnm(dc%lnsp), &
                zpnmt(nmp1*nnp1)


    !  Executable statements 

    IF (.NOT. ALLOCATED(pnm)) ALLOCATE (pnm(dc%lnsp,nhgl)); pnm=0.0_dp
    IF (.NOT. ALLOCATED(anm)) ALLOCATE (anm(dc%lnsp,nhgl)); anm=0.0_dp

    IF (.NOT. ALLOCATED(pnmd)) ALLOCATE (pnmd   (dc%lnsp))
    IF (.NOT. ALLOCATED(anmd)) ALLOCATE (anmd   (dc%lnsp))
    IF (.NOT. ALLOCATED(rnmd)) ALLOCATE (rnmd   (dc%lnsp))
    IF (.NOT. ALLOCATED(pnmi)) ALLOCATE (pnmi   (dc%lnsp))
    IF (.NOT. ALLOCATED(anmi)) ALLOCATE (anmi   (dc%lnsp))
    IF (.NOT. ALLOCATED(pnmiuv)) ALLOCATE (pnmiuv (dc%lnsp))
    IF (.NOT. ALLOCATED(anmiuv)) ALLOCATE (anmiuv (dc%lnsp))

    ! 1. Set up some constants and initiate scan

    DO jgl = 1, nhgl
      zmu   = gl_gmu(jgl)
      zrcst = 1.0_dp / gl_cst(jgl)

      ! 2. Compute, reorder and normalise Legendre
      !    polynomials and their meridional derivatives.

      ! 2.1 Compute *Legendre polynomials

      CALL phcs(zpnmt, zanmt, nmp1, nnp1, zmu)

      ! *phcs* computes polynomials for a rhomboidal truncation.

      ! 2.2 Reorder *Legendre polynomials

      CALL reord(zpnm, zpnmt, zanm, zanmt)

      ! 2.3 Normalise *Legendre polynomials

      zpnm(:) = zpnm(:) * 2.0_dp

      ! *anm* is divided by -(1.-mu**2) in order to get *d(pnm)/dmu.*

      zanm(:) = -zanm(:) * 2.0_dp * zrcst

      ! 3. Write Legendre polynomials to module fields

      pnm(:,jgl) = zpnm(:)
      anm(:,jgl) = zanm(:)

    END DO

  END SUBROUTINE inileg

!------------------------------------------------------------------------------
  SUBROUTINE reord(ppnm, ptp, phnm, pth)

    ! Description:
    !
    ! Reorders the *Legendre polynomials and their derivatives generated
    ! on a rhomboidal truncation by *phcs* to get them on the truncation
    ! used by the model.
    !
    ! Method:
    !
    ! *reord* is called from *inicom*. The polynomials are
    ! received *(ptp,pth)* and returned *(ppnm,phnm)* as subroutine
    ! arguments. The other parameters used are obtained from modules.
    !
    ! Results:
    ! *pnm* and *hnm* are returned in the natural order (i.e
    ! column wise ) on the truncation used by the model.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, March 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,       ONLY: nnp1

    !  Array arguments 
    REAL(dp) :: phnm(:), ppnm(:), pth(:), ptp(:)

    !  Local scalars: 
    INTEGER :: imp, inp, jm, jn
    INTEGER :: j


    !  Executable statements 

    ! 1. Reorders polynomials for a pentagonal truncation

    DO j = 1, dc%nlm    ! loop over zonal wavenumbers on this pe
      jm   = dc%lm(j)   ! zonal wave number
      inp  = dc%nlnp(j) ! lengt of m column
      imp  = dc%nlmp(j) ! offset of m column on this pe
      DO jn = 1, inp
        ppnm(imp+jn) = ptp((jm)*nnp1+jn)
        phnm(imp+jn) = pth((jm)*nnp1+jn)
      END DO
    END DO

  END SUBROUTINE reord
!------------------------------------------------------------------------------
  SUBROUTINE legmod(kpass, pnmt, anmt, rnmt)

    ! Description:
    !
    ! Calculate modified Legendre polynomials for direct tran.
    ! (grid-point to spherical harmonics)
    !
    ! Method:
    !
    ! *legmod* computes the modified Legendre polynomials 'pnmd', 
    ! 'anmd', 'rnmd' (one latitudinal index if kpass is present) for direc
    ! Legendre transforms using the normalised Legendre polynomials
    ! 'pnm', 'anm' stored in module 'mo_legendre'
    !
    ! *call* legmod(kpass)
    !
    ! *kpass*     The value of the main loop control variable
    !             (half latitudinal index)
    !             when *legmod* is called
    !
    ! Authors:
    !
    ! D. W. Dent, ECMWF, May 1984, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A.Rhodin        MPI, Sep 1999, changes for parallel version
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_physical_constants,     ONLY: earth_radius
    USE mo_gaussgrid,     ONLY: gl_gw

    !  Arguments 
    INTEGER , INTENT(in)  , OPTIONAL :: kpass      ! latitude index
    REAL(dp)    , INTENT(out) , OPTIONAL :: pnmt(:,:)  ! legendre coeffitients
    REAL(dp)    , INTENT(out) , OPTIONAL :: anmt(:,:)  !   present if kpass
    REAL(dp)    , INTENT(out) , OPTIONAL :: rnmt(:,:)  !   is not present

    !  Local scalars: 
    REAL(dp)    :: fnnp, z2w, zsa2
    INTEGER :: imp, innp, inp, jn, j

    !  Executable statements 

    !-- 1. Compute modified Legendre polynomials

    !      *pnmd*=2*w*pnmd
    !      *anmd*=2*w*anmd
    !      *rnmd*=2*w*(-n(n+1)/a**2)*pnmd

    IF (PRESENT (kpass)) THEN

      z2w = 2.0_dp*gl_gw(kpass)

      pnmd = pnm(:,kpass) * z2w
      anmd = anm(:,kpass) * z2w

      zsa2 = 1.0_dp/earth_radius**2

      DO j = 1, dc%nlm

        inp  = dc%nlnp(j)
        imp  = dc%nlmp(j)
        innp = dc%lm  (j)
        DO jn = 1, inp
          fnnp = REAL(innp+jn,dp)
          rnmd(imp+jn) = -zsa2*pnmd(imp+jn)*(fnnp-1.0_dp)*fnnp
        END DO
      END DO

    ELSE 

      DO j=1, dc%nlat/2    
        z2w = 2.0_dp * gl_gw(j)
        pnmt(j,:) = pnm(:,j) * z2w
        anmt(j,:) = anm(:,j) * z2w
      END DO

      zsa2 = 1.0_dp / earth_radius**2

      DO j = 1, dc%nlm
        inp  = dc%nlnp(j)
        imp  = dc%nlmp(j)
        innp = dc%lm  (j)
        DO jn = 1, inp
          fnnp = REAL(innp+jn,dp)
          rnmt(:,imp+jn) = -zsa2*pnmt(:,imp+jn)*(fnnp-1.0_dp)*(fnnp)
        END DO
      END DO

    ENDIF

  END SUBROUTINE legmod
!------------------------------------------------------------------------------
  SUBROUTINE leginv(kpass, pnmit, anmit, pnmiuvt, anmiuvt)

    ! Description:
    !
    ! Calculate modified Legendre polynomials for inverse tra.
    !
    ! Method:
    !
    ! *leginv* computes the modified Legendre polynomials for inver
    ! Legendre transforms using the normalised Legendre polynomials
    ! 'pnm', 'anm' stored in module 'mo_legendre'.
    !
    ! *call* *leginv(kpass)*
    !
    ! *kpass*      The value of the main loop control variable
    !              (half latitudinal index)
    !
    ! Authors:
    !
    ! D. W. Dent, ECMWF, May 1984, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A.Rhodin        MPI, Feb 2000, changes for parallel version
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_physical_constants,     ONLY: earth_radius
    USE mo_gaussgrid,     ONLY: gl_cst

    !  Arguments 
    INTEGER , INTENT(in)  , OPTIONAL :: kpass
    REAL(dp), INTENT(out) , OPTIONAL :: pnmit  (:,:)  !
    REAL(dp), INTENT(out) , OPTIONAL :: anmit  (:,:)  !
    REAL(dp), INTENT(out) , OPTIONAL :: pnmiuvt(:,:)  ! legendre coeffitients
    REAL(dp), INTENT(out) , OPTIONAL :: anmiuvt(:,:)  ! if kpass is not present

    !  Local scalars: 
    REAL(dp)    :: fnnp, zcst, zsa
    INTEGER :: imp, innp, inp, is, jm, jn, nhgl

    !  Executable statements 

    !-- 1. Compute modified Legendre polynomials
    !      suitable for inverse Legendre transforms of u and v

    !      *pnmiuv*=pnmi*a*m/(n(n+1)).
    !      *anmiuv*=anmi*a*(1.-mu**2)/(n(n+1)).

    IF(PRESENT(kpass)) THEN

      pnmi = pnm(:,kpass)
      anmi = anm(:,kpass)

      zcst = gl_cst(kpass)

      DO jm = 1, dc%nlm

        inp  = dc%nlnp(jm)
        imp  = dc%nlmp(jm)
        innp = dc%lm  (jm)

        IF (innp==0) THEN
          is = 2
        ELSE
          is = 1
        END IF

        DO jn = is, inp
          fnnp = REAL(innp+jn,dp)
          pnmiuv(imp+jn) = pnmi(imp+jn)*earth_radius*innp / ((fnnp-1._dp)*fnnp)
          anmiuv(imp+jn) = anmi(imp+jn)*earth_radius*zcst / ((fnnp-1._dp)*fnnp)
        END DO

      END DO

      IF (dc%nlnm0 > 0) THEN
        pnmiuv(dc%nlmp(1)+1) = 0.0_dp
        anmiuv(dc%nlmp(1)+1) = 0.0_dp
      ENDIF

      !-- 2. Compute inverse legendre polynomials *anmi*=*anmi/a*

      zsa = 1.0_dp/earth_radius

      anmi(:) = anmi(:)*zsa

    ELSE

      nhgl = dc%nlat/2

      DO jm = 1, dc%nlm

        inp  = dc%nlnp(jm)
        imp  = dc%nlmp(jm)
        innp = dc%lm  (jm)

        IF (innp==0) THEN
          is = 2
        ELSE
          is = 1
        END IF

        DO jn = is, inp
          fnnp = REAL(innp+jn,dp)
          pnmiuvt(:,imp+jn) = pnm(imp+jn,:)*earth_radius*innp         /((fnnp-1.0_dp)*fnnp)
          anmiuvt(:,imp+jn) = anm(imp+jn,:)*earth_radius*gl_cst(:nhgl)/((fnnp-1.0_dp)*fnnp)
        END DO
      END DO

      IF (dc%nlnm0 > 0) THEN
        pnmiuvt(:,dc%nlmp(1)+1) = 0.0_dp
        anmiuvt(:,dc%nlmp(1)+1) = 0.0_dp
      ENDIF

      !-- 2. Compute inverse legendre polynomials *anmi*=*anmi/a*

      zsa = 1.0_dp / earth_radius
      anmit(:,:) = TRANSPOSE(anm(:,:)) * zsa
      pnmit(:,:) = TRANSPOSE(pnm(:,:))

    ENDIF

  END SUBROUTINE leginv
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_legendre
    !
    ! deallocate module variables
    !
    IF (ALLOCATED(pnm))    DEALLOCATE (pnm)
    IF (ALLOCATED(anm))    DEALLOCATE (anm)
    IF (ALLOCATED(pnmd))   DEALLOCATE (pnmd)
    IF (ALLOCATED(anmd))   DEALLOCATE (anmd)
    IF (ALLOCATED(rnmd))   DEALLOCATE (rnmd)
    IF (ALLOCATED(pnmi))   DEALLOCATE (pnmi)
    IF (ALLOCATED(anmi))   DEALLOCATE (anmi)
    IF (ALLOCATED(pnmiuv)) DEALLOCATE (pnmiuv)
    IF (ALLOCATED(anmiuv)) DEALLOCATE (anmiuv)
  END SUBROUTINE cleanup_legendre
!------------------------------------------------------------------------------
END MODULE mo_legendre
