#ifdef __xlC__
@PROCESS STRICT
#endif
#if defined (__SX__)
#define VECTOR 1
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

SUBROUTINE ffti

#ifdef MKL_DFT
  USE mo_mkl_dft,       ONLY: mkl_dft
#elif ESSL_DFT
  USE mo_essl_dft,      ONLY: essl_dft
#else
  USE mo_fft992,        ONLY: fft992
#endif
  USE mo_buffer_fft,    ONLY: fftz, nvar
  USE mo_decomposition, ONLY: dc => local_decomposition
#if defined (VECTOR) && defined (_OPENMP)
  USE omp_lib,          ONLY: omp_get_thread_num, &
                              omp_get_max_threads
#endif
#ifdef _PROFILE
  USE mo_profile,       ONLY: trace_start, trace_stop
#endif

  IMPLICIT NONE

  INTEGER :: inc, isign
  INTEGER :: nlon, nlp2, nlev, nlat
#if defined (VECTOR) && defined (_OPENMP)
  INTEGER :: ivar, tid, nthreads, chunk, rest
  INTEGER, ALLOCATABLE, SAVE :: istart(:), icount(:)
#elif defined (MKL_DFT) || defined (ESSL_DFT)
#else
  INTEGER :: ilat, ivar
#endif
  LOGICAL :: col_1d

#ifdef _PROFILE
  CALL trace_start ('ffti', 31)
#endif

!-- 2. Inverse Fourier transforms

!-- 2.1 Set constants

  inc    = 1
  isign  = 1
  nlon   = dc% nlon
  nlp2   = nlon + 2
  nlev   = dc% nflevp1
  nlat   = dc% nflat
  col_1d = dc% col_1d

!-- 2.2 fft(vo, d, t, alps, u, v, dtl, dtm, dalpsl, dalpsm, dudl, dvdl)

  IF (.NOT.col_1d) THEN

#ifndef _OPENMP

#ifdef VECTOR

    CALL fft992(fftz,inc,nlp2,nlon,nvar*nlev*nlat,isign)

#elif MKL_DFT

    CALL mkl_dft(fftz,inc,nlp2,nlon,nvar*nlev*nlat,isign)

#elif ESSL_DFT

    CALL essl_dft(fftz,nlon,nlp2,nvar*nlev*nlat,isign)

#else

    DO ivar = 1, nvar
      DO ilat = 1, nlat
        CALL fft992(fftz(1,1,ilat,ivar),inc,nlp2,nlon,nlev,isign)
      ENDDO
    ENDDO

#endif

#else

#ifdef VECTOR

    IF (.NOT. ALLOCATED(istart)) THEN
      nthreads = omp_get_max_threads()
      ALLOCATE(istart(0:nthreads), icount(0:nthreads))
      istart(0) = 1
      DO tid = 0, nthreads-1
        chunk = nlat/nthreads 
        rest  = MOD(nlat, nthreads)
        if (tid < rest) chunk = chunk+1
        icount(tid) = chunk
        istart(tid+1) = istart(tid)+chunk
      ENDDO
    ENDIF

!$OMP PARALLEL PRIVATE(tid)
    tid = omp_get_thread_num()
    DO ivar = 1, nvar
        CALL fft992(fftz(1,1,istart(tid),ivar),inc,nlp2,nlon,nlev*icount(tid),isign)
    END DO
!$OMP END PARALLEL

#elif MKL_DFT

!$OMP PARALLEL PRIVATE(ilat)
    DO ivar = 1, nvar
!$OMP DO
      DO ilat = 1, nlat
        CALL mkl_dft(fftz(1,1,ilat,ivar),inc,nlp2,nlon,nlev,isign)
      ENDDO
!$OMP END DO
    ENDDO
!$OMP END PARALLEL

#elif ESSL_DFT

    CALL finish('ffti','OpenMP parallel ESSL DFT not supported yet.')

#else

!$OMP PARALLEL PRIVATE(ilat)
    DO ivar = 1, nvar
!$OMP DO
      DO ilat = 1, nlat
        CALL fft992(fftz(1,1,ilat,ivar),inc,nlp2,nlon,nlev,isign)
      ENDDO
!$OMP END DO
    ENDDO
!$OMP END PARALLEL

#endif

#endif

  ENDIF

#ifdef _PROFILE
  CALL trace_stop ('ffti', 31)
#endif

END SUBROUTINE ffti
