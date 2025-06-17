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
#ifdef __STANDALONE
SUBROUTINE init_decomposition(nproca, nprocb,  &
                              ngl, nlon, nlev, &
                              nm, nn, nk,      &
                              nproma,          &
                              lcolumn)
#else
SUBROUTINE init_decomposition
#endif

#ifndef __STANDALONE
  USE mo_control,       ONLY: nproca, nprocb, ngl, nlon, nlev, nm      &
                            , nn, nk, lcolumn, nproma, nprocio
  USE mo_column,        ONLY: lat_1d, lon_1d !, inicolumn
  USE mo_advection,     ONLY: iadvec
  USE mo_decomposition, ONLY: local_decomposition,                     & 
                              global_decomposition,                    &
                              decompose, debug_seriell,                &
                              print_decomposition
#else
  USE mo_decomposition, ONLY: local_decomposition,                     & 
                              global_decomposition,                    &
!                              print_decomposition,                     &
                              decompose, debug_seriell
#endif
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_mpi,           ONLY: p_nprocs, p_pe, p_io, p_parallel,        &
                              p_parallel_io

  USE mo_transpose,     ONLY: indx
  USE mo_util_string,   ONLY: separator
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_max_threads, omp_set_num_threads
#endif

  IMPLICIT NONE

#ifdef __STANDALONE
  ! mo_control
  INTEGER, INTENT(in) :: nproca, nprocb, ngl, nlon, nlev, nm, nn, nk, nproma
  LOGICAL, INTENT(in) :: lcolumn

  ! mo_column
  INTEGER :: lat_1d(1), lon_1d(1)
 
  ! mo_advection, = 3 selects tpcore
  INTEGER :: iadvec = 3
#endif

  INTEGER :: p,i
  
  LOGICAL :: lrot, lfull_m
  INTEGER :: debug

  lrot          = .TRUE.  ! true: no rotation of longitudes
  lfull_m       = .FALSE. ! true: full m-columns per PE
  debug_seriell = .FALSE. ! true: same results as ser.version if nprocb == 1

  ! debug = 0,1 : PE 0 takes full domain (always no rotation)
  !         -1  : no special treatment of PE 0
  !          0  : gather from PE 0
  !          1  : gather from PEs > 0

#ifdef NOMPI
  IF (nproca*nprocb > 1) THEN
     CALL message ('', &
          'Parallel execution selected with wrong compiler options')
     CALL message ('', 'Please, recompile without -DNOMPI') 
     CALL finish ('init_decomposition', 'Program aborted')
  END IF
#endif

  IF (p_nprocs == nproca*nprocb) THEN
     debug = -1
  ELSE IF ( p_nprocs == nproca*nprocb+1) THEN
    IF (nprocio > 0) THEN
      CALL finish ('init_decomposition',                                &
           'Parallel I/O and parallel debug not suported')
    ENDIF
    debug = 0
#ifdef _OPENMP
     ! reset number of threads to one so that OpenMP can be checked as well
     IF (p_pe == p_io) THEN
       CALL OMP_SET_NUM_THREADS(1)
     ENDIF
#endif
  ELSE
     CALL finish ('init_decomposition',                                &
          'Number of runtime PEs doesn''t fit nproca*nprocb(+1)')
  END IF

  IF (p_parallel .AND. p_parallel_io) THEN 
     CALL message ('', '')     
     WRITE (message_text,'(a,i4,a,i3,a,i3)')                           &
          ' Total number of PEs: ', p_nprocs,                          &
          ' set A: ', nproca, ' set B: ', nprocb
     CALL message ('', message_text)     
  ENDIF

  ALLOCATE (global_decomposition(1:p_nprocs))

  ! derive decomposition

  IF (lcolumn) THEN
    CALL decompose (global_decomposition, 0, nproca, nprocb, ngl,      &
         nlon, nlev, nm, nn, nk, iadvec, norot=lrot, debug=debug,      &
         lfull_m=lfull_m, lats_1d=lat_1d(1), lons_1d=lon_1d(1))
  ELSE
    CALL decompose (global_decomposition, nproma, nproca, nprocb,      &
         ngl, nlon, nlev, nm, nn, nk, iadvec, norot=lrot, debug=debug, &
         lfull_m=lfull_m)
  ENDIF

  ! keep index values, not id values

  DO p = 1, p_nprocs
     global_decomposition(p)%mapmesh(:,:) =                            &
        indx(global_decomposition(p)%mapmesh(:,:),global_decomposition)
  END DO

  ! copy global decomposition table entry to local decomposition

  DO i = 1, p_nprocs
     IF (global_decomposition(i)% pe == p_pe) THEN
        local_decomposition = global_decomposition(i)
     END IF
  END DO

  IF (p_parallel_io) THEN
    CALL message('', separator)
    CALL message('',' Blocking information:')
    WRITE (message_text,'(a,l5)') '   reordering (lreg)           = ', local_decomposition%lreg
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') '   number of blocks (ngplks)   = ', local_decomposition%ngpblks
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') '   size of blocks (nproma)     = ', local_decomposition%nproma
    CALL message('', message_text)
    WRITE (message_text,'(a,i5)') '   size of last block (npromz) = ', local_decomposition%npromz
    CALL message('', message_text)
    CALL message('', separator)
  END IF

!!$#ifdef __STANDALONE
!!$  ! decomposition printout
!!$
!!$  IF (p_parallel_io .AND. p_parallel) THEN
!!$     DO i = 1, p_nprocs
!!$        CALL message('', separator)
!!$        CALL print_decomposition(global_decomposition(i))
!!$     END DO
!!$     CALL message('', separator)
!!$     CALL message('', '')
!!$  END IF
!!$#endif

END SUBROUTINE init_decomposition
