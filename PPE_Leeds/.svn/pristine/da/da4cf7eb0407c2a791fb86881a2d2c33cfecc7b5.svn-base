!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_tr_omp_decomposition

#ifdef _OPENMP
  USE omp_lib
#endif

  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_decomposition, ONLY: ldc => local_decomposition

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: max_nthreads=64

  LOGICAL, PARAMETER :: ldebug = .FALSE.

  ! type for 1d-thread-space
  TYPE thspd1_type
    INTEGER :: begin
    INTEGER :: end
  END TYPE thspd1_type

  PUBLIC :: thspd1_type

  ! public data (for read only):

#ifdef _CRAYFTN
  ! number of available threads
  INTEGER, PUBLIC :: nthreads = 0

  ! thread space for full model latitudes-segments
  TYPE(thspd1_type), PUBLIC :: thsp_lat(2,0:max_nthreads-1) 
  LOGICAL,           PUBLIC :: uniform_thsp_lat
  
  ! thread space for mpi-task-local latitudes
  TYPE(thspd1_type), PUBLIC :: thsp_glat(0:max_nthreads-1) 
  LOGICAL,           PUBLIC :: uniform_thsp_glat

  ! thread space for mpi-task-local fourier latitudes
  TYPE(thspd1_type), PUBLIC :: thsp_flat(0:max_nthreads-1) 
  LOGICAL,           PUBLIC :: uniform_thsp_flat
  LOGICAL,           PUBLIC :: nonzero_thsp_flat

  ! thread space for nproma-blocks
  TYPE(thspd1_type), PUBLIC :: thsp_blk(0:max_nthreads-1) 
  LOGICAL,           PUBLIC :: uniform_thsp_blk
#else
  ! number of available threads
  INTEGER, PUBLIC, PROTECTED :: nthreads = 0

  ! thread space for full model latitudes-segments
  TYPE(thspd1_type), PUBLIC, PROTECTED :: thsp_lat(2,0:max_nthreads-1) 
  LOGICAL,           PUBLIC, PROTECTED :: uniform_thsp_lat
  
  ! thread space for mpi-task-local latitudes
  TYPE(thspd1_type), PUBLIC, PROTECTED :: thsp_glat(0:max_nthreads-1) 
  LOGICAL,           PUBLIC, PROTECTED :: uniform_thsp_glat

  ! thread space for mpi-task-local fourier latitudes
  TYPE(thspd1_type), PUBLIC, PROTECTED :: thsp_flat(0:max_nthreads-1) 
  LOGICAL,           PUBLIC, PROTECTED :: uniform_thsp_flat
  LOGICAL,           PUBLIC, PROTECTED :: nonzero_thsp_flat

  ! thread space for nproma-blocks
  TYPE(thspd1_type), PUBLIC, PROTECTED :: thsp_blk(0:max_nthreads-1) 
  LOGICAL,           PUBLIC, PROTECTED :: uniform_thsp_blk
#endif

  ! public units:

  PUBLIC :: init_omp_decomposition
  ! exported from either omp_lib or an internal replacement without omp usage
  PUBLIC :: omp_get_thread_num
  PUBLIC :: omp_get_max_threads
  PUBLIC :: omp_in_parallel

CONTAINS

#ifndef _OPENMP
  ! functions replacing openmp functions used 
  INTEGER FUNCTION omp_get_thread_num()
    omp_get_thread_num = 0
  END FUNCTION omp_get_thread_num

  INTEGER FUNCTION omp_get_max_threads()
    omp_get_max_threads = 1
  END FUNCTION omp_get_max_threads
  LOGICAL FUNCTION omp_in_parallel()
    omp_in_parallel = .FALSE.
  END FUNCTION omp_in_parallel
#endif

  SUBROUTINE init_omp_decomposition(nt)
    INTEGER, INTENT(in), OPTIONAL :: nt

    IF (PRESENT(nt)) THEN
      nthreads = nt
    ELSE
#ifdef _OPENMP
      nthreads = omp_get_max_threads()
#else
      nthreads = 1
#endif
    ENDIF

    IF (nthreads > max_nthreads) THEN
      CALL finish('init_omp_decomposition','max_nthreads too small',1)
    ENDIF

    CALL init_thsp

  END SUBROUTINE init_omp_decomposition

  SUBROUTINE init_thsp
    TYPE(thspd1_type) :: tmp_lat(0:max_nthreads-1) 
    INTEGER :: it, nlat

    nlat = ldc%nlat

    ! 1-d full model latitudes
    CALL decompose_thsp(nlat/2, tmp_lat, uniform_thsp_lat, 'lat')

    DO it = 0, nthreads-1
      ! first segment (northern latitudes)
      thsp_lat(1,it)%begin = tmp_lat(it)%begin
      thsp_lat(1,it)%end   = tmp_lat(it)%end
      ! second segment (southern latitudes)
      thsp_lat(2,it)%begin = nlat - thsp_lat(1,it)%end + 1
      thsp_lat(2,it)%end   = nlat - thsp_lat(1,it)%begin + 1
    ENDDO
    ! 1-d gridspace latitudes :
    CALL decompose_thsp(ldc%nglat, thsp_glat, uniform_thsp_glat, 'glat')

    ! 1-d fourierspace latitudes :
    IF (ldc%nflat == 0) THEN
      nonzero_thsp_flat = .FALSE.
      thsp_flat(:)      = thspd1_type(-1,-2)
      uniform_thsp_flat = .TRUE.
    ELSE
      IF ( ldc%nglat /= ldc%nflat ) THEN
        CALL finish('init_thsp','ldc%nglat /= ldc%nflat',1)
      ENDIF
      nonzero_thsp_flat =.TRUE.
      thsp_flat(:)      = thsp_glat(:)
      uniform_thsp_flat = uniform_thsp_glat
    ENDIF

    ! 1-d grid-block space:
    CALL decompose_thsp(ldc%ngpblks, thsp_blk, uniform_thsp_blk, 'blk')

  END SUBROUTINE init_thsp

  SUBROUTINE decompose_thsp(nx, thsp_x, uniform_thsp_x, label)
    INTEGER,                    INTENT(in)  :: nx
    TYPE(thspd1_type),          INTENT(out) :: thsp_x(0:)
    LOGICAL,                    INTENT(out) :: uniform_thsp_x
    CHARACTER(len=*), OPTIONAL, INTENT(in)  :: label

    INTEGER :: min_chunk, rest, it
    INTEGER :: chunk(0:nthreads-1)

    thsp_x(:) = thspd1_type(-1,-2)
    
    min_chunk = nx/nthreads
    rest      = MOD(nx,nthreads)

    IF (rest == 0) THEN
      chunk(:)       = min_chunk
      uniform_thsp_x = .TRUE.
    ELSE
      chunk(0:rest-1)        = min_chunk+1
      chunk(rest:nthreads-1) = min_chunk
      uniform_thsp_x         = .FALSE.
    ENDIF
    IF (ldebug .AND. PRESENT(label)) THEN
      WRITE(message_text,*) label, ' uniform_thsp:', uniform_thsp_x
      CALL message('decompose_thsp',message_text)
    ENDIF

    thsp_x(0)%begin = 1
    thsp_x(0)%end   = chunk(0)
    DO it = 1, nthreads-1
      thsp_x(it)%begin = thsp_x(it-1)%end+1
      thsp_x(it)%end   = thsp_x(it-1)%end+chunk(it)

    ENDDO
    IF (ldebug .AND. PRESENT(label)) THEN
      DO it = 0, nthreads-1
        WRITE(message_text,*) label,' (it, begin, end) ',it, thsp_x(it)%begin, thsp_x(it)%end
      CALL message('decompose_thsp',message_text)
      ENDDO
    ENDIF

  END SUBROUTINE decompose_thsp

END MODULE mo_tr_omp_decomposition
