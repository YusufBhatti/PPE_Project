#ifdef __xlC__
@PROCESS HOT
@PROCESS ALIAS(NOPTEOVRLP,NOARYOVRLP)
#endif

!OCL NOALIAS

#if defined(__uxp__) || defined(__SX__) || defined(ES) || defined(_UNICOSMP) || defined(__xlC__)
#define FAST_AND_DIRTY 1
#endif

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE sym1

  ! Description:
  !
  ! This subroutine computes symmetric and antisymmetric
  ! parts of Fourier components.
  !
  ! Method:
  !
  ! sym1 computes symetric and antisymmetric parts of
  ! Fourier components in 2 steps: The contribution of northern
  ! is added when nrow is odd and of southern hemisphere when
  ! nrow is even.
  !
  ! Results:
  ! The contributions are returned as symmetric and antisymmetric
  ! Fourier coefficients.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, March 1982, original source
  ! J. K. Gibson, ECMWF, August 1983, Adapted for multi-tasking
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! E. Tschirschnitz, NEC, February 2003, Optimization
  ! L. Kornblueh, MPI, February 2003, clean code for subscript checking
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_control,       ONLY: ltimer
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_sym1
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_memory_f,      ONLY: fadl, fadm, far, fatp1, faul, fazl, fazm, &
                              fsdl, fsdm, fsr, fstp1, fsul, fszl, fszm
  USE mo_buffer_fft,    ONLY: ldl, ldm, ltm1, lrh, lvol, lvom ,lul, lalpsm1


  IMPLICIT NONE

  !  Local scalars: 

  INTEGER :: nlev, nlevp1, nmp1, nlnm0, ngl, nhgl
  INTEGER :: jm
#ifdef FAST_AND_DIRTY
  INTEGER :: jl, jlev, js, ix

  INTEGER :: ixl(dc% nlm*dc% nflev*2)
  INTEGER :: ixr(dc% nlm*dc% nflev*2)
  INTEGER :: ixlp(dc% nlm*dc% nflev*2)
#endif

  !  Executable statements 

  ! local array bounds on this PE

  nlev   = dc% nflev   ! number of levels
  nlevp1 = dc% nflevp1 ! number of levels including highest level
  nmp1   = dc% nlm     ! number of spectral coefficients m
  nlnm0  = dc% nlnm0   ! number of spectral coefficients n for m=0
  ngl    = dc% nlat    ! number of latitudes
  nhgl   = ngl / 2

#ifdef FAST_AND_DIRTY
!DIR$ CONCURRENT
!DIR$ PREFERSTREAM
!$OMP PARALLEL PRIVATE(jm,js,jlev)
  IF (ltimer) CALL timer_start(timer_sym1)
!$OMP DO
  DO jm = 1, nmp1
    DO js = 1,2
      DO jlev = 1,nlev
        ixl(jlev+(js-1)*nlev+(jm-1)*nlev*2)  = jlev+(js-1)*nlev+(jm-1)*nlev*2
        ixlp(jlev+(js-1)*nlev+(jm-1)*nlev*2) = jlev+(js-1)*nlevp1+(jm-1)*nlevp1*2
        ixr(jlev+(js-1)*nlev+(jm-1)*nlev*2)  = 2*jm+(js-2)+(jlev-1)*2*nmp1
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO

!DIR$ CONCURRENT
!DIR$ PREFERSTREAM
!$OMP DO  PRIVATE(jl,ix)
  DO jl = 1, nhgl
!CDIR NODEP
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
    do ix = 1, nlev*2*nmp1
      fszl (ixl(ix),1,1,jl)  = (lvol(ixr(ix),1,jl)+lvol(ixr(ix),1,ngl+1-jl))*0.5_dp
      fazm (ixl(ix),1,1,jl)  = (lvom(ixr(ix),1,jl)-lvom(ixr(ix),1,ngl+1-jl))*0.5_dp
      fsdl (ixl(ix),1,1,jl)  = (ldl (ixr(ix),1,jl)+ldl (ixr(ix),1,ngl+1-jl))*0.5_dp
      fadm (ixl(ix),1,1,jl)  = (ldm (ixr(ix),1,jl)-ldm (ixr(ix),1,ngl+1-jl))*0.5_dp
      fsr  (ixl(ix),1,1,jl)  = (lrh (ixr(ix),1,jl)+lrh (ixr(ix),1,ngl+1-jl))*0.5_dp
      fazl (ixl(ix),1,1,jl)  = (lvol(ixr(ix),1,jl)-lvol(ixr(ix),1,ngl+1-jl))*0.5_dp
      fszm (ixl(ix),1,1,jl)  = (lvom(ixr(ix),1,jl)+lvom(ixr(ix),1,ngl+1-jl))*0.5_dp
      fadl (ixl(ix),1,1,jl)  = (ldl (ixr(ix),1,jl)-ldl (ixr(ix),1,ngl+1-jl))*0.5_dp
      fsdm (ixl(ix),1,1,jl)  = (ldm (ixr(ix),1,jl)+ldm (ixr(ix),1,ngl+1-jl))*0.5_dp
      far  (ixl(ix),1,1,jl)  = (lrh (ixr(ix),1,jl)-lrh (ixr(ix),1,ngl+1-jl))*0.5_dp
      fstp1(ixlp(ix),1,1,jl) = (ltm1(ixr(ix),1,jl)+ltm1(ixr(ix),1,ngl+1-jl))*0.5_dp
      fatp1(ixlp(ix),1,1,jl) = (ltm1(ixr(ix),1,jl)-ltm1(ixr(ix),1,ngl+1-jl))*0.5_dp
    END DO
  END DO
!$OMP END DO
#else
!$OMP PARALLEL PRIVATE(jm)
  IF (ltimer) CALL timer_start(timer_sym1)
!$OMP DO
  DO jm = 1, nmp1
!DIR$ CONCURRENT
    fszl (:    ,1,jm,:)=(lvol(2*jm-1,:,:nhgl)+lvol(2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fszl (:    ,2,jm,:)=(lvol(2*jm  ,:,:nhgl)+lvol(2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fazm (:    ,1,jm,:)=(lvom(2*jm-1,:,:nhgl)-lvom(2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fazm (:    ,2,jm,:)=(lvom(2*jm  ,:,:nhgl)-lvom(2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fsdl (:    ,1,jm,:)=(ldl (2*jm-1,:,:nhgl)+ldl (2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fsdl (:    ,2,jm,:)=(ldl (2*jm  ,:,:nhgl)+ldl (2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fadm (:    ,1,jm,:)=(ldm (2*jm-1,:,:nhgl)-ldm (2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fadm (:    ,2,jm,:)=(ldm (2*jm  ,:,:nhgl)-ldm (2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fsr  (:    ,1,jm,:)=(lrh (2*jm-1,:,:nhgl)+lrh (2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fsr  (:    ,2,jm,:)=(lrh (2*jm  ,:,:nhgl)+lrh (2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fstp1(:nlev,1,jm,:)=(ltm1(2*jm-1,:,:nhgl)+ltm1(2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fstp1(:nlev,2,jm,:)=(ltm1(2*jm  ,:,:nhgl)+ltm1(2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fazl (:    ,1,jm,:)=(lvol(2*jm-1,:,:nhgl)-lvol(2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fazl (:    ,2,jm,:)=(lvol(2*jm  ,:,:nhgl)-lvol(2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fszm (:    ,1,jm,:)=(lvom(2*jm-1,:,:nhgl)+lvom(2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fszm (:    ,2,jm,:)=(lvom(2*jm  ,:,:nhgl)+lvom(2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fadl (:    ,1,jm,:)=(ldl (2*jm-1,:,:nhgl)-ldl (2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fadl (:    ,2,jm,:)=(ldl (2*jm  ,:,:nhgl)-ldl (2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fsdm (:    ,1,jm,:)=(ldm (2*jm-1,:,:nhgl)+ldm (2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fsdm (:    ,2,jm,:)=(ldm (2*jm  ,:,:nhgl)+ldm (2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    far  (:    ,1,jm,:)=(lrh (2*jm-1,:,:nhgl)-lrh (2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    far  (:    ,2,jm,:)=(lrh (2*jm  ,:,:nhgl)-lrh (2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fatp1(:nlev,1,jm,:)=(ltm1(2*jm-1,:,:nhgl)-ltm1(2*jm-1,:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    fatp1(:nlev,2,jm,:)=(ltm1(2*jm  ,:,:nhgl)-ltm1(2*jm  ,:,ngl:nhgl+1:-1))*0.5_dp
  END DO
!$OMP END DO
#endif  

  IF (nlevp1 > nlev) THEN
!DIR$ CONCURRENT
!DIR$ PREFERSTREAM
!$OMP DO PRIVATE(jm)
    DO jm = 1, nmp1
!DIR$ CONCURRENT
      fstp1(nlevp1,1,jm,:) = (lalpsm1(2*jm-1,:nhgl)+lalpsm1(2*jm-1,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
      fstp1(nlevp1,2,jm,:) = (lalpsm1(2*jm  ,:nhgl)+lalpsm1(2*jm  ,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
      fatp1(nlevp1,1,jm,:) = (lalpsm1(2*jm-1,:nhgl)-lalpsm1(2*jm-1,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
      fatp1(nlevp1,2,jm,:) = (lalpsm1(2*jm  ,:nhgl)-lalpsm1(2*jm  ,ngl:nhgl+1:-1))*0.5_dp
    END DO
!$OMP END DO
  END IF

  IF (nlnm0 > 0) THEN
!$OMP WORKSHARE
!DIR$ CONCURRENT
    fsul(:,:) = (lul(:,:nhgl)+lul(:,ngl:nhgl+1:-1))*0.5_dp
!DIR$ CONCURRENT
    faul(:,:) = (lul(:,:nhgl)-lul(:,ngl:nhgl+1:-1))*0.5_dp
!$OMP END WORKSHARE
  ENDIF

  IF (ltimer) CALL timer_stop(timer_sym1)
!$OMP END PARALLEL

END SUBROUTINE sym1
