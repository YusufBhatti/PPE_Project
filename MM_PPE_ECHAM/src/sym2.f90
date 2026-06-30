#ifdef __xlC__
@PROCESS HOT
@PROCESS ALIAS(NOARYOVRLP,NOPTEOVRLP)
#endif

#if defined(__uxp__) || defined(__SX__) || defined(ES)
#define FAST_AND_DIRTY 1
#endif

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE sym2

  ! Description:
  !
  ! This subroutine computes Ffourier components from
  ! their symmetric and antisymmetric parts and
  ! retrieves u and v from their rotational and divergent parts.
  !
  ! Method:
  !
  ! sym2 is called from scan1.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, March 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! E. Tschirschnitz, NEC, February 2003, Optimization
  ! L. Kornblueh, MPI, February 2003, clean code for subscript checking
  ! 
  ! for more details see file AUTHORS
  !

#ifdef FAST_AND_DIRTY
  USE mo_kind,          ONLY: dp
#endif
  USE mo_memory_f,      ONLY: fad, fadu0, fatp, fatpm, fau, fau0, fav, favo, &
                              fsd, fsdu0, fstp, fstpm, fsu, fsu0, fsv, fsvo
  USE mo_buffer_fft,    ONLY: lalps, ld, lt, lu, lv, lvo, ldalpsm, ldtm, &
                              lu0, ldu0, lul
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_control,       ONLY: ltimer   
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_sym2

  IMPLICIT NONE

  !  Local scalars: 

  INTEGER :: jlev, jm, jlan, jlas, nlat
  INTEGER :: nlev, nlevp1, nmp1

#ifdef FAST_AND_DIRTY
  INTEGER , PARAMETER :: jlev1=1, jm1=1

  INTEGER :: ix

  INTEGER :: ixl(2*dc%nllev*dc% nlm)
  INTEGER :: ixr(2*dc%nllev*dc% nlm)
  INTEGER :: ixrp(2*dc%nllev*dc% nlm)
#endif

  !  Executable statements 

  ! Array bounds on local PE

  nmp1   = dc% nlm
  nlev   = dc% nllev
  nlevp1 = dc% nllevp1
  nlat   = dc% nlat

#ifdef FAST_AND_DIRTY
  DO jlev = 1, nlev
    DO jm = 1, 2*nmp1
      ixl(jm+(jlev-1)*2*nmp1)      = jm+(jlev-1)*2*nmp1
    ENDDO
!OCL NOVREC, NOALIAS
    DO jm = 1, nmp1
      ixr(2*jm-1+(jlev-1)*2*nmp1)  = jlev+(jm-1)*2*(nlev+0.5_dp)
      ixr(2*jm+(jlev-1)*2*nmp1)    = jlev+nlev+(jm-1)*2*(nlev+0.5_dp)
      ixrp(2*jm-1+(jlev-1)*2*nmp1) = jlev+(jm-1)*2*(nlevp1+0.5_dp)
      ixrp(2*jm+(jlev-1)*2*nmp1)   = jlev+nlevp1+(jm-1)*2*(nlevp1+0.5_dp)
    ENDDO
  ENDDO
#endif

!$OMP PARALLEL

    IF (ltimer) CALL timer_start(timer_sym2)

#ifdef FAST_AND_DIRTY
!$OMP DO PRIVATE(jlan,jlas,jlev,ix,jm)
#else
!$OMP DO PRIVATE(jlan,jlas,jlev,jm)
#endif
  DO jlan = 1, nlat/2         ! northern latitude or half row index
    jlas = nlat-jlan+1        ! southern latitude row index
    
    !-- 1. Compute Fourier components

#ifdef FAST_AND_DIRTY
!CDIR NODEP
!OCL NOVREC, NOALIAS
!IBM* ASSERT(NODEPS)
    DO ix = 1, nlev*2*nmp1
      lvo (ixl(ix),jlev1,jlan) = fsvo(ixr(ix),1,jm1,jlan) + favo(ixr(ix),1,jm1,jlan)
      ld  (ixl(ix),jlev1,jlan) = fsd(ixr(ix),1,jm1,jlan)  + fad(ixr(ix),1,jm1,jlan)
      lu  (ixl(ix),jlev1,jlan) = fsu(ixr(ix),1,jm1,jlan)  + fau(ixr(ix),1,jm1,jlan)
      lv  (ixl(ix),jlev1,jlan) = fsv(ixr(ix),1,jm1,jlan)  + fav(ixr(ix),1,jm1,jlan)
      lvo (ixl(ix),jlev1,jlas) = fsvo(ixr(ix),1,jm1,jlan) - favo(ixr(ix),1,jm1,jlan)
      ld  (ixl(ix),jlev1,jlas) = fsd(ixr(ix),1,jm1,jlan)  - fad(ixr(ix),1,jm1,jlan)
      lu  (ixl(ix),jlev1,jlas) = fsu(ixr(ix),1,jm1,jlan)  - fau(ixr(ix),1,jm1,jlan)
      lv  (ixl(ix),jlev1,jlas) = fsv(ixr(ix),1,jm1,jlan)  - fav(ixr(ix),1,jm1,jlan)
      lt  (ixl(ix),jlev1,jlan) = fstp(ixrp(ix),1,jm1,jlan)  + fatp(ixrp(ix),1,jm1,jlan)
      ldtm(ixl(ix),jlev1,jlan) = fstpm(ixrp(ix),1,jm1,jlan) + fatpm(ixrp(ix),1,jm1,jlan)
      lt  (ixl(ix),jlev1,jlas) = fstp(ixrp(ix),1,jm1,jlan)  - fatp(ixrp(ix),1,jm1,jlan)
      ldtm(ixl(ix),jlev1,jlas) = fstpm(ixrp(ix),1,jm1,jlan) - fatpm(ixrp(ix),1,jm1,jlan)
    ENDDO
#else
    !-- 1.1 Northern hemisphere

    DO jlev = 1, nlev
      
!DIR$ IVDEP
!OCL NOVREC, NOALIAS
!DIR$ VECTOR ALWAYS
      DO jm = 1, nmp1
        lvo (2*jm-1,jlev,jlan) = fsvo(jlev,1,jm,jlan) + favo(jlev,1,jm,jlan)
        lvo (2*jm  ,jlev,jlan) = fsvo(jlev,2,jm,jlan) + favo(jlev,2,jm,jlan)
        ld  (2*jm-1,jlev,jlan) = fsd(jlev,1,jm,jlan)  + fad(jlev,1,jm,jlan)
        ld  (2*jm  ,jlev,jlan) = fsd(jlev,2,jm,jlan)  + fad(jlev,2,jm,jlan)
        lt  (2*jm-1,jlev,jlan) = fstp(jlev,1,jm,jlan) + fatp(jlev,1,jm,jlan)
        lt  (2*jm  ,jlev,jlan) = fstp(jlev,2,jm,jlan) + fatp(jlev,2,jm,jlan)
        lu  (2*jm-1,jlev,jlan) = fsu(jlev,1,jm,jlan)  + fau(jlev,1,jm,jlan)
        lu  (2*jm  ,jlev,jlan) = fsu(jlev,2,jm,jlan)  + fau(jlev,2,jm,jlan)
        lv  (2*jm-1,jlev,jlan) = fsv(jlev,1,jm,jlan)  + fav(jlev,1,jm,jlan)
        lv  (2*jm,  jlev,jlan) = fsv(jlev,2,jm,jlan)  + fav(jlev,2,jm,jlan)
        ldtm(2*jm-1,jlev,jlan) = fstpm(jlev,1,jm,jlan)+ fatpm(jlev,1,jm,jlan)
        ldtm(2*jm  ,jlev,jlan) = fstpm(jlev,2,jm,jlan)+ fatpm(jlev,2,jm,jlan)
      END DO

      !-- 1.2 Southern hemisphere

!DIR$ IVDEP
!OCL NOVREC, NOALIAS
!DIR$ VECTOR ALWAYS
      DO jm = 1, nmp1
        lvo (2*jm-1,jlev,jlas) = fsvo(jlev,1,jm,jlan) - favo(jlev,1,jm,jlan)
        lvo (2*jm  ,jlev,jlas) = fsvo(jlev,2,jm,jlan) - favo(jlev,2,jm,jlan)
        ld  (2*jm-1,jlev,jlas) = fsd(jlev,1,jm,jlan)  - fad(jlev,1,jm,jlan)
        ld  (2*jm  ,jlev,jlas) = fsd(jlev,2,jm,jlan)  - fad(jlev,2,jm,jlan)
        lt  (2*jm-1,jlev,jlas) = fstp(jlev,1,jm,jlan) - fatp(jlev,1,jm,jlan)
        lt  (2*jm  ,jlev,jlas) = fstp(jlev,2,jm,jlan) - fatp(jlev,2,jm,jlan)
        lu  (2*jm-1,jlev,jlas) = fsu(jlev,1,jm,jlan)  - fau(jlev,1,jm,jlan)
        lu  (2*jm  ,jlev,jlas) = fsu(jlev,2,jm,jlan)  - fau(jlev,2,jm,jlan)
        lv  (2*jm-1,jlev,jlas) = fsv(jlev,1,jm,jlan)  - fav(jlev,1,jm,jlan)
        lv  (2*jm,  jlev,jlas) = fsv(jlev,2,jm,jlan)  - fav(jlev,2,jm,jlan)
        ldtm(2*jm-1,jlev,jlas) = fstpm(jlev,1,jm,jlan)- fatpm(jlev,1,jm,jlan)
        ldtm(2*jm  ,jlev,jlas) = fstpm(jlev,2,jm,jlan)- fatpm(jlev,2,jm,jlan)
      END DO
      
    END DO
#endif

    IF (nlevp1>nlev) THEN
!DIR$ IVDEP
!OCL NOVREC
!CDIR NODEP
!IBM* ASSERT(NODEPS)
!DIR$ VECTOR ALWAYS
      DO jm = 1, nmp1
#ifdef FAST_AND_DIRTY
         lalps  (2*jm-1,jlan) = fstp (nlevp1,1,jm,jlan) + fatp (nlevp1,1,jm,jlan)
         ldalpsm(2*jm-1,jlan) = fstpm(nlevp1,1,jm,jlan) + fatpm(nlevp1,1,jm,jlan)
         lalps  (2*jm-1,jlas) = fstp (nlevp1,1,jm,jlan) - fatp (nlevp1,1,jm,jlan)
         ldalpsm(2*jm-1,jlas) = fstpm(nlevp1,1,jm,jlan) - fatpm(nlevp1,1,jm,jlan)
         ldalpsm(2*jm  ,jlas) = fstpm(2*nlevp1,1,jm,jlan) - fatpm(2*nlevp1,1,jm,jlan)
         ldalpsm(2*jm  ,jlan) = fstpm(2*nlevp1,1,jm,jlan) + fatpm(2*nlevp1,1,jm,jlan)
         lalps  (2*jm  ,jlas) = fstp (2*nlevp1,1,jm,jlan) - fatp (2*nlevp1,1,jm,jlan)
         lalps  (2*jm  ,jlan) = fstp (2*nlevp1,1,jm,jlan) + fatp (2*nlevp1,1,jm,jlan)
#else
        lalps  (2*jm-1,jlan) = fstp (nlevp1,1,jm,jlan) + fatp (nlevp1,1,jm,jlan)
        ldalpsm(2*jm-1,jlan) = fstpm(nlevp1,1,jm,jlan) + fatpm(nlevp1,1,jm,jlan)
        lalps  (2*jm-1,jlas) = fstp (nlevp1,1,jm,jlan) - fatp (nlevp1,1,jm,jlan)
        ldalpsm(2*jm-1,jlas) = fstpm(nlevp1,1,jm,jlan) - fatpm(nlevp1,1,jm,jlan)
        ldalpsm(2*jm  ,jlas) = fstpm(nlevp1,2,jm,jlan) - fatpm(nlevp1,2,jm,jlan)
        ldalpsm(2*jm  ,jlan) = fstpm(nlevp1,2,jm,jlan) + fatpm(nlevp1,2,jm,jlan)
        lalps  (2*jm  ,jlas) = fstp (nlevp1,2,jm,jlan) - fatp (nlevp1,2,jm,jlan)
        lalps  (2*jm  ,jlan) = fstp (nlevp1,2,jm,jlan) + fatp (nlevp1,2,jm,jlan)
#endif

      END DO
    ENDIF
    
    ! arrays with spectral coefficients m=0 only

    IF (dc% nlnm0 > 0) THEN
      lu0 (:,jlas) = fsu0 (:,jlan) - fau0 (:,jlan)
      ldu0(:,jlas) = fsdu0(:,jlan) - fadu0(:,jlan)
      lu0 (:,jlan) = fsu0 (:,jlan) + fau0 (:,jlan)
      ldu0(:,jlan) = fsdu0(:,jlan) + fadu0(:,jlan)
      lul (:,jlan) = lu (1,:,jlan)
      lul (:,jlas) = lu (1,:,jlas)
    ENDIF
  END DO

!$OMP END DO

  IF (ltimer) CALL timer_stop(timer_sym2)

!$OMP END PARALLEL

END SUBROUTINE sym2
