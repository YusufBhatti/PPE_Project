!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_util_buffer_pool

  USE, INTRINSIC :: iso_c_binding

  IMPLICIT NONE

  INTERFACE
    INTEGER(c_int) FUNCTION init_buffer_pool(entries) &
         BIND(C, NAME='init_buffer_pool')
      IMPORT :: c_int
      INTEGER(c_int), VALUE :: entries
    END FUNCTION init_buffer_pool
  END INTERFACE

  INTERFACE
    INTEGER(c_int) FUNCTION setup_buffer_pool(pd) &
         BIND(C, NAME='setup_buffer_pool')
      IMPORT :: c_int
      INTEGER(c_int), VALUE :: pd
    END FUNCTION setup_buffer_pool
  END INTERFACE

  INTERFACE
    INTEGER(c_int) FUNCTION allocate_buffer_pool(pd, allocation_type) &
         BIND(C, NAME='allocate_buffer_pool')
      IMPORT :: c_int
      INTEGER(c_int), VALUE :: pd
      INTEGER(c_int), VALUE :: allocation_type
    END FUNCTION allocate_buffer_pool
  END INTERFACE

  INTERFACE
    INTEGER(c_int) FUNCTION private_copy_chunk_to_buffer_pool(pd, chunk_no, chunk) &
         BIND(C, NAME='copy_chunk_to_buffer_pool')
      IMPORT :: c_int, c_double
      INTEGER(c_int), VALUE :: pd
      INTEGER(c_int), VALUE :: chunk_no
      REAL(c_double), DIMENSION(*) :: chunk
    END FUNCTION private_copy_chunk_to_buffer_pool
  END INTERFACE

  INTERFACE
    INTEGER(c_int) FUNCTION buffer_pool_entry(pd, size, dims, chunk) &
         BIND(C, NAME='buffer_pool_entry') 
      IMPORT :: c_int, c_double
      INTEGER(c_int), VALUE :: pd
      INTEGER(c_int), VALUE :: size
      INTEGER(c_int), VALUE :: dims      
      REAL(c_double), DIMENSION(*) :: chunk
    END FUNCTION buffer_pool_entry
  END INTERFACE

  INTERFACE
    SUBROUTINE print_test() BIND(C, NAME='print_test')
    END SUBROUTINE print_test
  END INTERFACE

CONTAINS

  FUNCTION copy_chunk_to_buffer_pool(pd, chunk_no, chunk) RESULT(iret)
    INTEGER, INTENT(in) :: pd
    INTEGER, INTENT(in) :: chunk_no
    DOUBLE PRECISION, DIMENSION(*) :: chunk
    INTEGER :: iret
    iret = private_copy_chunk_to_buffer_pool(pd, chunk_no, chunk)
  END FUNCTION copy_chunk_to_buffer_pool

END MODULE mo_util_buffer_pool
