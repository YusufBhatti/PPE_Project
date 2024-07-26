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
SUBROUTINE ltd

  ! Description:
  !
  ! Direct Legendre transforms.
  !
  ! Method:
  !
  ! This subroutine performs direct *Legendre transforms for
  ! the divergence equation,
  ! the temperature and surface pressure equations,
  ! the vorticity equation,
  ! the mean wind.
  !
  ! *ltd* is called from *scan1sl*
  !
  ! Results:
  ! *ltd* adds in the spectral arrays:-
  !      *ld*  the contribution of the current latitude line
  !      *ltp* the contribution of the current latitude line
  !      *lvo* the contribution of the current latitude line
  !      *lu0* the contribution of the current latitude line
  !
  ! Authors:
  !
  ! D. W. Dent, ECMWF, December 1984, original source
  ! U. Schlese, DKRZ, in 1991, and 1994, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! L. Kornblueh, MPI, November 2002, optimization
  !
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_control,       ONLY: ltimer
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_ltd   
  USE mo_memory_ls,     ONLY: ld, ltp, lu0, lvo
  USE mo_memory_f,      ONLY: fadl, fadm, far, fatp1, faul, fazl, fazm, fsdl, &
                              fsdm, fsr, fstp1, fsul, fszl, fszm
  USE mo_legendre,      ONLY: legmod
  USE mo_decomposition, ONLY: lc => local_decomposition

  IMPLICIT NONE

  !  Local loop bounds

  INTEGER          :: nllev, nllevp1, nlmp1, nlnm0, lnsp
  INTEGER ,POINTER :: nlmp(:), nlnp(:)

  !  Global bounds
  INTEGER          :: nhgl

  !  Local scalars:

  INTEGER :: ims, ins, inu, iu, j, j0, jh, jhr, jj, kk
  INTEGER :: kf, km

  REAL(dp), POINTER :: ful (:,:)

  REAL(dp), SAVE, ALLOCATABLE :: pnml(:,:)
  REAL(dp), SAVE, ALLOCATABLE :: anml(:,:)
  REAL(dp), SAVE, ALLOCATABLE :: rnml(:,:)

  REAL(dp) :: fpstack((2*lc%nllev+lc%nllevp1)*2,lc%nlm,lc%nlat/2,2)
  REAL(dp) :: fastack(2*lc%nllev*2,lc%nlm,lc%nlat/2,2)
  REAL(dp) :: fr2stack(2*lc%nllev,lc%nlm,lc%nlat/2,2)
  REAL(dp) :: lpstack((2*lc%nllev+lc%nllevp1)*2+1,lc%lnsp)
  REAL(dp) :: lastack(2*lc%nllev*2+1,lc%lnsp)
  REAL(dp) :: lrstack(lc%nllev*2+1,lc%lnsp)

  !  External subroutines 

  !  Executable statements

  !-- Set local loop bounds

  nllev   =  lc% nllev    ! number of levels
  nllevp1 =  lc% nllevp1  ! number of levels + 1
  nlmp1   =  lc% nlm      ! number of m wave numbers
  nlnm0   =  lc% nlnm0    ! number of coefficients for m=0
  lnsp    =  lc% lnsp     ! number of complex spectral coefficients on this pe
  nlmp    => lc% nlmp     ! displacement of the first point of columns
  nlnp    => lc% nlnp     ! number of points on each column
  nhgl    =  lc% nlat/2   ! global halv number of gaussian latitudes

  IF ( .NOT. ALLOCATED(anml) ) THEN

    ALLOCATE( anml( lc% nlat/2, lc% lnsp))
    ALLOCATE( pnml( lc% nlat/2, lc% lnsp))
    ALLOCATE( rnml( lc% nlat/2, lc% lnsp))

    ! derive local wavenumber index
    ! calculate legendre coefficents for each latitude

    CALL legmod(pnmt=pnml,anmt=anml,rnmt=rnml)
  END IF

  !-- 1. Legendre transforms

    !-- 1.1 Transforms for d, vo, t and p

!$OMP PARALLEL PRIVATE(jj,jh,iu,j0,j,ims,ins,jhr,kk,kf,km,inu,ful)

  IF (ltimer) CALL timer_start(timer_ltd)

    jh = 1   ! northern hemisphere
    ! Loop over latitudes
!$OMP DO
    DO jhr = 1, nhgl

      DO j=1,nlmp1
        DO kf = 1,2
!DIR$ CONCURRENT
          DO kk = 1, nllev
            km = (kf-1)*nllev+kk 
            fpstack(km,j,jhr,jh)             = fsdl(kk,kf,j,jhr)
            fpstack(km+nllev*2,j,jhr,jh)     = fszl(kk,kf,j,jhr)
            fastack(km,j,jhr,jh)             = fadm(kk,kf,j,jhr)
            fastack(km+nllev*2,j,jhr,jh)     = fazm(kk,kf,j,jhr)
            fr2stack(km,j,jhr,jh)            = fsr(kk,kf,j,jhr)
          ENDDO
        ENDDO
        DO kf = 1,2
          DO kk=1,nllevp1
            km = (kf-1)*nllevp1+kk 
            fpstack(km+nllev*4,j,jhr,jh)     = fstp1(kk,kf,j,jhr)
          ENDDO
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

    jh = 2   ! southern hemisphere
    ! Loop over latitudes
!$OMP DO
    DO jhr = 1, nhgl

      DO j=1,nlmp1
        DO kf = 1,2
!DIR$ CONCURRENT
          DO kk = 1, nllev
            km = (kf-1)*nllev+kk 
            fpstack(km,j,jhr,jh)             = fadl(kk,kf,j,jhr)
            fpstack(km+nllev*2,j,jhr,jh)     = fazl(kk,kf,j,jhr)
            fastack(km,j,jhr,jh)             = fsdm(kk,kf,j,jhr)
            fastack(km+nllev*2,j,jhr,jh)     = fszm(kk,kf,j,jhr)
            fr2stack(km,j,jhr,jh)            = far(kk,kf,j,jhr)
          ENDDO
        ENDDO
        DO kf = 1,2
          DO kk=1,nllevp1
            km = (kf-1)*nllevp1+kk 
            fpstack(km+nllev*4,j,jhr,jh)     = fatp1(kk,kf,j,jhr)
          ENDDO
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

!$OMP DO
!CDIR NOVECTOR
    DO jj = 1, 2*nlmp1
      jh = 1 + (jj-1)/nlmp1
      iu = 1 - (jj-1)/nlmp1
      j0 = mod(jj-1,nlmp1) + 1
      j = (nlmp1 - j0/2) * MOD(j0,2) + j0/2 * MOD(j0+1,2)
      ims  = nlmp(j) - iu   ! offset to local m columns  (spectral coef.)
      ins  = nlnp(j) + iu   ! column length

      IF (ins>=2) then
        CALL dgemm("N","N",                                     &
                   4*nllev+2*nllevp1,ins/2,nhgl,                &
                   1.0_dp,                                      &
                   fpstack(1,j,1,jh),(4*nllev+2*nllevp1)*nlmp1, &
                   pnml(1,ims+2),     2*nhgl,                   &
                   0.0_dp,                                      &
                   lpstack(1,ims+2),  2*((2*nllev+nllevp1)*2+1))
        CALL dgemm("N","N",                          &
                   4*nllev,          ins/2,nhgl,     &
                   -1.0_dp,                          &
                   fastack(1,j,1,jh), 4*nllev*nlmp1, &
                   anml(1,ims+2),     2*nhgl,        &
                   0.0_dp,                           &
                   lastack(1,ims+2),  2*(2*nllev*2+1))
        CALL dgemm("N","N",                          &
                   2*nllev,          ins/2,nhgl,     &
                   1.0_dp,                           &
                   fr2stack(1,j,1,jh),2*nllev*nlmp1, &
                   rnml(1,ims+2),     2*nhgl,        &
                   0.0_dp,                           &
                   lrstack(1,ims+2),  2*(nllev*2+1))
      ENDIF

    END DO
!$OMP END DO

  !-- 1.2 Transforms for the mean wind

  IF ( nlnm0 > 0 ) THEN
!$OMP DO
!CDIR NOVECTOR
    DO jh = 1,2
      iu = 2 - jh
      
      inu = (nlnm0+iu)/2
      
      IF (jh==1) THEN
        ful  => fsul 
      ELSE
        ful  => faul
      END IF

      CALL dgemm("N","N",             &
                 nllev,inu,nhgl,      &
                 1.0_dp,              &
                 ful(1,1),    nllev,  &
                 pnml(1,2-iu),2*nhgl, &
                 1.0_dp,              &
                 lu0(1,2-iu), 2*nllev)
    END DO
!$OMP END DO
  END IF

!$OMP DO
  DO jhr = 1, lnsp
    DO kf = 1,2
!DIR$ CONCURRENT
      DO kk = 1, nllev
        km = (kf-1)*nllev+kk 
        ld(kk,kf,jhr)=lpstack(km,jhr)+lastack(km,jhr)+lrstack(km,jhr)
        lvo(kk,kf,jhr)=lpstack(2*nllev+km,jhr)+lastack(2*nllev+km,jhr)
      ENDDO
    ENDDO
    DO kf = 1,2
      DO kk=1,nllevp1
        km = (kf-1)*nllevp1+kk 
        ltp(kk,kf,jhr)=lpstack(4*nllev+km,jhr)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

  IF (ltimer) CALL timer_stop(timer_ltd)

!$OMP END PARALLEL

END SUBROUTINE ltd
