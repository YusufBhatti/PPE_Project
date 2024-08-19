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
SUBROUTINE lti

  ! Description:
  !
  ! Inverse Legendre transforms
  !
  ! Method:
  !
  ! This subroutine performs inverse *legendre transforms
  !
  ! *lti* is called from *scan2*
  !
  ! Results:
  ! *lti* computes the *fourier components
  ! for the current latitude line.
  !
  ! Authors:
  !
  ! D. W. Dent, ECMWF, December 1984, original source
  ! U. Schlese, DKRZ, December 1994, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! L. Kornblueh, MPI, November 2002, optimization
  ! L. Kornblueh, MPI, february 2004, optimization
  ! R. Smith, MPI, February 2004, optimization
  ! S. Shingu, NEC, March 2004, optimization (ES) 
  !
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp 
  USE mo_control,       ONLY: ltimer
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_lti
  USE mo_memory_ls,     ONLY: ld, ltp, lu0, lvo
  USE mo_memory_f,      ONLY: fad, fadu0, fatp, fatpm, fau, fau0, fav, favo, &
                              fsd, fsdu0, fstp, fstpm, fsu, fsu0, fsv, fsvo
  USE mo_legendre,      ONLY: leginv
  USE mo_decomposition, ONLY: lc => local_decomposition

  IMPLICIT NONE

  !  Local loop bounds

  INTEGER          :: nllev, nllevp1, nlmp1, nlnm0, lnsp
  INTEGER ,POINTER :: nlmp(:), nlnp(:)

  !  Global bounds
  INTEGER          :: nhgl

  !  Local scalars:

  INTEGER :: ims, ins, irow, iu, j, j0, jh, jhr, jj, jl
  INTEGER :: jm

  !  Local arrays:

  REAL(dp) :: pnmit   (lc% nlat/2 ,lc% lnsp) ,&
              anmit   (lc% nlat/2 ,lc% lnsp) ,&
              pnmiuvt (lc% nlat/2 ,lc% lnsp) ,&
              anmiuvt (lc% nlat/2 ,lc% lnsp)

  REAL(dp) :: ftmp1 ( lc% nlat /2 , lc% nllev  , 2, lc% nlm, 2 ), &
              ftmp2 ( lc% nlat /2 , lc% nllev  , 2, lc% nlm, 2 ), &
              ftmp3 ( lc% nlat    , lc% nllev  , 2, lc% nlm, 2 ), &
              ftmp4 ( lc% nlat    , lc% nllev  , 2, lc% nlm, 2 ), &
              ftmpp1( lc% nlat    , lc% nllevp1, 2, lc% nlm, 2 )

  REAL(dp), SAVE, ALLOCATABLE :: pnmil(:,:), pnmiuvl(:,:)

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

  !-- Set initial values for transforms (zero coefficients (1,m=0 n=0)))

  IF (nlnm0 > 0) lvo(:,1,1) = 0.0_dp
  IF (nlnm0 > 0) ld (:,1,1) = 0.0_dp

  IF ( .NOT. ALLOCATED(pnmil) ) THEN

    ALLOCATE( pnmil   (2*nhgl,lc% lnsp))
    ALLOCATE( pnmiuvl (2*nhgl,lc% lnsp))

    CALL leginv(pnmit=pnmit ,anmit=anmit ,pnmiuvt=pnmiuvt ,anmiuvt=anmiuvt)
!$OMP PARALLEL PRIVATE(irow,jhr)
!$OMP DO
    DO jhr = 1, lnsp
      DO irow = 1, nhgl
        pnmil  (     irow,jhr) = pnmit (irow,jhr)
        pnmil  (nhgl+irow,jhr) = anmit (irow,jhr)
        pnmiuvl(     irow,jhr) = pnmiuvt(irow,jhr)
        pnmiuvl(nhgl+irow,jhr) = anmiuvt(irow,jhr)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ENDIF

  !-- Inverse *Legendre transforms

!$OMP PARALLEL PRIVATE(jj,j0,j,jh,jm,ims,ins,iu,jl,irow)

  IF (ltimer) CALL timer_start(timer_lti)

!$OMP DO
!CDIR$ NOVECTOR
!CDIR NOVECTOR
  DO jj = 1, 2*nlmp1
    jh = 1 + (jj-1)/nlmp1
    iu = 1 - (jj-1)/nlmp1
    j0 = mod(jj-1,nlmp1) + 1
    j = (nlmp1 - j0/2) * MOD(j0,2) + j0/2 * MOD(j0+1,2)
    ims  = nlmp(j) - iu   ! offset to local m columns  (spectral coef.)
    ins  = nlnp(j) + iu   ! column length
    IF (ins>=2) THEN
      CALL dgemm('N','T',                      &
                 nhgl,  2*nllev,  ins/2,       &
                 1.0_dp,                       &
                 pnmil(1,ims+2),    nhgl*2*2,  &
                 ld(1,1,ims+2),     nllev*2*2, &
                 0.0_dp,                       &
                 ftmp1(1,1,1,j,jh), nhgl)
      CALL dgemm('N','T',                      &
                 nhgl,  2*nllev,  ins/2,       &
                 1.0_dp,                       &
                 pnmil(1,ims+2),    nhgl*2*2,  &
                 lvo(1,1,ims+2),    nllev*2*2, &
                 0.0_dp,                       &
                 ftmp2(1,1,1,j,jh), nhgl)
      CALL dgemm('N','T',                        &
                 2*nhgl,2*nllevp1,ins/2,         &
                 1.0_dp,                         &
                 pnmil(1,ims+2),    nhgl*2*2,    &
                 ltp(1,1,ims+2),    nllevp1*2*2, &
                 0.0_dp,                         &
                 ftmpp1(1,1,1,j,jh),2*nhgl)
      CALL dgemm('N','T',                      &
                 2*nhgl,2*nllev,ins/2,         &
                 1.0_dp,                       &
                 pnmiuvl(1,ims+2),  nhgl*2*2,  &
                 ld(1,1,ims+2),     nllev*2*2, &
                 0.0_dp,                       &
                 ftmp3(1,1,1,j,jh), 2*nhgl)
      CALL dgemm('N','T',                      &
                 2*nhgl,2*nllev,ins/2,         &
                 1.0_dp,                       &
                 pnmiuvl(1,ims+2),  nhgl*2*2,  &
                 lvo(1,1,ims+2),    nllev*2*2, &
                 0.0_dp,                       &
                 ftmp4(1,1,1,j,jh), 2*nhgl)
    END IF

  END DO
!$OMP END DO

  jh=1
!$OMP DO
!CDIR$ NOVECTOR
!CDIR NOVECTOR
  DO j0 = 1, nlmp1
    j = (nlmp1 - j0/2) * MOD(j0,2) + j0/2 * MOD(j0+1,2)
    ins=(nlnp(j)+1)/2

    IF (ins>0) THEN
      DO jl = 1, nllev
        DO jm = 1,2
          DO irow = 1, nhgl
            fsd (jl,jm,j,irow) = ftmp1(irow,jl,jm,j,jh)
            fsvo(jl,jm,j,irow) = ftmp2(irow,jl,jm,j,jh)
          ENDDO
        ENDDO
      ENDDO
      DO jl = 1, nllevp1
        DO jm = 1,2
          DO irow = 1, nhgl
            fstp (jl,jm,j,irow) = ftmpp1(     irow,jl,jm,j,jh)
            fatpm(jl,jm,j,irow) = ftmpp1(nhgl+irow,jl,jm,j,jh)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      fsd  (:,:,j,:) = 0.0_dp
      fsvo (:,:,j,:) = 0.0_dp
      fstp (:,:,j,:) = 0.0_dp
      fatpm(:,:,j,:) = 0.0_dp
      ftmp3(1:nhgl       ,:,:,j,jh) = 0.0_dp
      ftmp3(nhgl+1:2*nhgl,:,:,j,jh) = 0.0_dp
      ftmp4(1:nhgl       ,:,:,j,jh) = 0.0_dp
      ftmp4(nhgl+1:2*nhgl,:,:,j,jh) = 0.0_dp
    END IF

  END DO
!$OMP END DO

  jh=2
!$OMP DO
!CDIR$ NOVECTOR
!CDIR NOVECTOR
  DO j0 = 1, nlmp1
    j = (nlmp1 - j0/2) * MOD(j0,2) + j0/2 * MOD(j0+1,2)
    ins=(nlnp(j))/2
    IF (ins>0) THEN
      DO jl = 1, nllev
        DO jm = 1,2
          DO irow = 1, nhgl
            fad (jl,jm,j,irow) = ftmp1(irow,jl,jm,j,jh)
            favo(jl,jm,j,irow) = ftmp2(irow,jl,jm,j,jh)
          ENDDO
        ENDDO
      ENDDO
      DO jl = 1, nllevp1
        DO jm = 1,2
          DO irow = 1, nhgl
            fatp (jl,jm,j,irow) = ftmpp1(     irow,jl,jm,j,jh)
            fstpm(jl,jm,j,irow) = ftmpp1(nhgl+irow,jl,jm,j,jh)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      fad  (:,:,j,:) = 0.0_dp
      favo (:,:,j,:) = 0.0_dp
      fatp (:,:,j,:) = 0.0_dp
      fstpm(:,:,j,:) = 0.0_dp
      ftmp3(1:nhgl       ,:,:,j,jh) = 0.0_dp
      ftmp3(nhgl+1:2*nhgl,:,:,j,jh) = 0.0_dp
      ftmp4(1:nhgl       ,:,:,j,jh) = 0.0_dp
      ftmp4(nhgl+1:2*nhgl,:,:,j,jh) = 0.0_dp
    END IF

  END DO
!$OMP END DO

  !-- Combine rotational and divergent parts of u and v

!DIR$ CONCURRENT
!$OMP DO
  DO j = 1, nlmp1
    DO jl = 1, nllev
!DIR$ IVDEP
      DO irow = 1, nhgl
        fsu(jl,1,j,irow) =  ftmp4(nhgl+irow,jl,1,j,2) + ftmp3(irow,jl,2,j,1)
        fsu(jl,2,j,irow) =  ftmp4(nhgl+irow,jl,2,j,2) - ftmp3(irow,jl,1,j,1)
        fau(jl,1,j,irow) =  ftmp4(nhgl+irow,jl,1,j,1) + ftmp3(irow,jl,2,j,2)
        fau(jl,2,j,irow) =  ftmp4(nhgl+irow,jl,2,j,1) - ftmp3(irow,jl,1,j,2)
        fsv(jl,1,j,irow) =  ftmp4(irow,jl,2,j,1) - ftmp3(nhgl+irow,jl,1,j,2)
        fsv(jl,2,j,irow) = -ftmp4(irow,jl,1,j,1) - ftmp3(nhgl+irow,jl,2,j,2)
        fav(jl,1,j,irow) =  ftmp4(irow,jl,2,j,2) - ftmp3(nhgl+irow,jl,1,j,1)
        fav(jl,2,j,irow) = -ftmp4(irow,jl,1,j,2) - ftmp3(nhgl+irow,jl,2,j,1)
      END DO
    END DO
  END DO
!$OMP END DO
  
  !-- Transform mean wind

  IF ( nlnm0 > 0 ) THEN
!$OMP SECTIONS
!$OMP SECTION
    jh = 1
    ins = nlnm0+1
    CALL dgemm('N','T',              &
               nllev,nhgl,ins/2,     &
               1.0_dp,               &
               lu0(1,jh),  nllev*2,  &
               pnmil(1,jh),nhgl*2*2, &
               0.0_dp,               &
               fsu0 (1,1), nllev)
!$OMP SECTION
    jh = 1
    ins = nlnm0+1
    CALL dgemm('N','T',                   &
               nllev,nhgl,ins/2,          &
               1.0_dp,                    &
               lu0(1,jh),       nllev*2,  &
               pnmil(nhgl+1,jh),nhgl*2*2, &
               0.0_dp,                    &
               fadu0(1,1),      nllev)
!$OMP SECTION
    jh = 2
    ins = nlnm0
    CALL dgemm('N','T',              &
               nllev,nhgl,ins/2,     &
               1.0_dp,               &
               lu0(1,jh),  nllev*2,  &
               pnmil(1,jh),nhgl*2*2, &
               0.0_dp,               &
               fau0 (1,1), nllev)
!$OMP SECTION
    jh = 2
    ins = nlnm0
    CALL dgemm('N','T',                   &
               nllev,nhgl,ins/2,          &
               1.0_dp,                    &
               lu0(1,jh),       nllev*2,  &
               pnmil(nhgl+1,jh),nhgl*2*2, &
               0.0_dp,                    &
               fsdu0(1,1),      nllev)
!$OMP END SECTIONS
  END IF

  IF (ltimer) CALL timer_stop(timer_lti)
!$OMP END PARALLEL

END SUBROUTINE lti
