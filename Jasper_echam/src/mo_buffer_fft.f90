!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_buffer_fft

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: pe_decomposed

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fftz          ! Fourier space array
  PUBLIC :: fftl          ! Legendre
  PUBLIC :: fbm0          ! buffer with m=0 only (Fourier)
  PUBLIC :: lbm0          ! buffer with m=0 only (Legendre)
  PUBLIC :: construct_fft ! fourier space array deallocation routine
  PUBLIC :: destruct_fft  ! fourier space array allocation routine
  PUBLIC :: fd, ft, fu, fv, fvo, fdtl, fdtm, falps, fdalpsl, fdalpsm
  PUBLIC :: fdudl, fdvdl, fulz, fvmax
  PUBLIC :: ld, lt, lu, lv, lvo, ldtm, lalps, ldalpsm
  PUBLIC :: fu0, fdu0, ful, lu0, ldu0, lul
  PUBLIC :: ldl,  ldm, ltm1, lrh, lvol, lvom,        lalpsm1 !, lul
  PUBLIC :: fdm1, fdm, ftm1, frh, fvol, fvom, fvom1, falpsm1 !, ful
  PUBLIC :: nvar, nvar_fsls, nvar_lsfs, nvar0_fsls, nvar0_lsfs

  !
  ! common array for all variables
  !
  REAL(dp), ALLOCATABLE, TARGET :: fftz  (:,:,:,:) ! the fft is performed in fftz
  REAL(dp), ALLOCATABLE, TARGET :: fftl  (:,:,:,:) ! sym1,2  works in        fftl
  REAL(dp), ALLOCATABLE, TARGET :: fbm0  (:,:,:)   ! buffer with spectral c. m=0 only (Fourier)
  REAL(dp), ALLOCATABLE, TARGET :: lbm0  (:,:,:)   ! buffer with spectral c. m=0 only (Legendr)
  REAL(dp), ALLOCATABLE, TARGET :: fulz  (:,:)     ! buffer with zonal max of fu (Fourier)
  REAL(dp), ALLOCATABLE, TARGET :: fvmax (:)       ! buffer with horizonal max of !u!+!v! (Fourier)

  !
  ! pointers into fftz, fftl used by inverse fft
  !                inverse fft,    sym2
  !
  REAL(dp), POINTER :: fd     (:,:,:) ,ld   (:,:,:)
  REAL(dp), POINTER :: ft     (:,:,:) ,lt   (:,:,:)
  REAL(dp), POINTER :: fu     (:,:,:) ,lu   (:,:,:)
  REAL(dp), POINTER :: fdudl  (:,:,:) 
  REAL(dp), POINTER :: fv     (:,:,:) ,lv   (:,:,:)
  REAL(dp), POINTER :: fdvdl  (:,:,:) 
  REAL(dp), POINTER :: fvo    (:,:,:) ,lvo  (:,:,:)
  REAL(dp), POINTER :: fdtm   (:,:,:) ,ldtm (:,:,:)
  REAL(dp), POINTER :: fdtl   (:,:,:) 
  REAL(dp), POINTER :: falps  (:,:)   ,lalps  (:,:)
  REAL(dp), POINTER :: fdalpsl(:,:)
  REAL(dp), POINTER :: fdalpsm(:,:)   ,ldalpsm(:,:)

  REAL(dp), POINTER :: fu0    (:,:)   ,lu0    (:,:)
  REAL(dp), POINTER :: fdu0   (:,:)   ,ldu0   (:,:)
  REAL(dp), POINTER :: ful    (:,:)   ,lul    (:,:)

  !
  ! pointers into fftz, fftl used by direct fft
  !                direct fft,    sym1
  !
  REAL(dp), POINTER :: fdm1   (:,:,:), ldl    (:,:,:)
  REAL(dp), POINTER :: fdm    (:,:,:), ldm    (:,:,:)
  REAL(dp), POINTER :: ftm1   (:,:,:), ltm1   (:,:,:)
  REAL(dp), POINTER :: frh    (:,:,:), lrh    (:,:,:)
  REAL(dp), POINTER :: fvol   (:,:,:), lvol   (:,:,:)
  REAL(dp), POINTER :: fvom   (:,:,:), lvom   (:,:,:)
  REAL(dp), POINTER :: falpsm1(:,:),   lalpsm1(:,:)
  REAL(dp), POINTER :: fvom1  (:,:,:)

  INTEGER, PARAMETER :: nvar = 9       ! number of variables (4th index of zfft)
  INTEGER, PARAMETER :: nvar_fsls  = 6 ! number of variables spectral c. m=0 only  
  INTEGER, PARAMETER :: nvar_lsfs  = 9 ! number of variables spectral c. m=0 only
  INTEGER, PARAMETER :: nvar0_fsls = 1 ! number of variables spectral c. m=0 only  
  INTEGER, PARAMETER :: nvar0_lsfs = 3 ! number of variables spectral c. m=0 only
  

CONTAINS

  SUBROUTINE construct_fft (dc)

    TYPE (pe_decomposed), INTENT(in) :: dc   ! decomposition table

#ifdef _OPENMP
    INTEGER :: jlan, jlas
    INTEGER :: jlat
    INTEGER :: ivar
#endif

    !
    ! multi level arrays
    !
    ALLOCATE (fftz (dc% nlon+2, dc% nflevp1, dc% nflat, nvar))
    ALLOCATE (fftl (dc% nlm *2, dc% nflevp1, dc% nlat , nvar))

#ifdef _OPENMP
!$OMP PARALLEL
    ! sym2 affinity:
!$OMP DO PRIVATE(jlan, jlas)
    DO jlan = 1, dc%nlat/2
      jlas = dc%nlat-jlan+1
      fftl(:,:,jlan,:) = 0.0_dp
      fftl(:,:,jlas,:) = 0.0_dp
    ENDDO
!$OMP END DO
    ! ewd affinity:
!$OMP DO PRIVATE(jlat)
    DO jlat = 1, dc%nflat
      fftz(:,:,jlat,:) = 0.0_dp 
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
#else
    fftl(:,:,:,:) = 0.0_dp
    fftz(:,:,:,:) = 0.0_dp 
#endif

    !
    ! arrays with spectral coefficients m=0 only
    !
    ALLOCATE (lbm0 (            dc% nflev  , dc% nlat,  nvar0_lsfs))
    ALLOCATE (fbm0 (            dc% nflev  , dc% nflat, nvar0_lsfs))
    lbm0(:,:,:) = 0.0_dp
    fbm0(:,:,:) = 0.0_dp
    !
    ALLOCATE (fulz (            dc% nflev  , dc% nflat))
    ALLOCATE (fvmax(            dc% nflev))
    fulz(:,:)   = 0.0_dp
    fvmax(:)    = 0.0_dp
    !
    ! pointers for inverse transforms
    !
    fd    => fftz (:,1:dc% nflev,:,1) ;ld   => fftl (:,1:dc% nflev,:,1)
    ft    => fftz (:,1:dc% nflev,:,2) ;lt   => fftl (:,1:dc% nflev,:,2)
    fu    => fftz (:,1:dc% nflev,:,3) ;lu   => fftl (:,1:dc% nflev,:,3)
    fv    => fftz (:,1:dc% nflev,:,4) ;lv   => fftl (:,1:dc% nflev,:,4)
    fvo   => fftz (:,1:dc% nflev,:,5) ;lvo  => fftl (:,1:dc% nflev,:,5)
    fdtm  => fftz (:,1:dc% nflev,:,6) ;ldtm => fftl (:,1:dc% nflev,:,6)
    fdtl  => fftz (:,1:dc% nflev,:,7)
    fdudl => fftz (:,1:dc% nflev,:,8)
    fdvdl => fftz (:,1:dc% nflev,:,9)
    !
    ! single level arrays
    !
    IF (dc%nflev == dc%nflevp1) THEN
      NULLIFY (falps, fdalpsl, fdalpsm)
    ELSE
      fdalpsl => fftz (:,dc% nflevp1,:,1) 
      fdalpsm => fftz (:,dc% nflevp1,:,2) ;ldalpsm => fftl (:,dc% nflevp1,:,2)
      falps   => fftz (:,dc% nflevp1,:,3) ;lalps   => fftl (:,dc% nflevp1,:,3)
    ENDIF
    !
    ! zonal means (m=0 only)
    !
    ful => fbm0 (:,:,1)             ;lul  => lbm0 (:,:,1)
    fu0 => fbm0 (:,:,2)             ;lu0  => lbm0 (:,:,2)
    fdu0=> fbm0 (:,:,3)             ;ldu0 => lbm0 (:,:,3)
    !
    ! pointers for direct transforms
    !
    fdm1 => fftz (:,1:dc% nflev,:,1); ldl  => fftl (:,1:dc% nflev,:,1)
    fdm  => fftz (:,1:dc% nflev,:,2); ldm  => fftl (:,1:dc% nflev,:,2)
    ftm1 => fftz (:,1:dc% nflev,:,3); ltm1 => fftl (:,1:dc% nflev,:,3)
    frh  => fftz (:,1:dc% nflev,:,4); lrh  => fftl (:,1:dc% nflev,:,4)
    fvol => fftz (:,1:dc% nflev,:,5); lvol => fftl (:,1:dc% nflev,:,5)
    fvom => fftz (:,1:dc% nflev,:,6); lvom => fftl (:,1:dc% nflev,:,6)
    fvom1=> fftz (:,1:dc% nflev,:,7)
    ! ful, lul used for both direct and inverse transform
    !
    ! single level arrays
    !
    IF (dc%nflev == dc%nflevp1) THEN
      NULLIFY (falpsm1); NULLIFY (lalpsm1)
    ELSE
      falpsm1 => fftz(:,dc%nflevp1,:,3); lalpsm1 => fftl(:,dc%nflevp1,:,3)
    ENDIF

  END SUBROUTINE construct_fft

  SUBROUTINE destruct_fft
    DEALLOCATE (fftz)
    DEALLOCATE (fftl)
    DEALLOCATE (lbm0)
    DEALLOCATE (fbm0)
    DEALLOCATE (fulz)
    DEALLOCATE (fvmax)
  END SUBROUTINE destruct_fft

END MODULE mo_buffer_fft


