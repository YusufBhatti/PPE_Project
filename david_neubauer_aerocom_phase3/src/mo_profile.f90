!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_profile

  IMPLICIT NONE

  PRIVATE

#if defined(__xlC__) && defined(__hpm__)
#include "f_hpc.h"
#endif

  PUBLIC :: trace_init
  PUBLIC :: trace_finalize
  PUBLIC :: trace_start
  PUBLIC :: trace_stop

CONTAINS

  SUBROUTINE trace_init(name, pe) 
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER,          INTENT(in) :: pe

#if defined(__xlC__) && defined(__hpm__)
    call f_hpminit(pe, TRIM(name))
#endif

  END SUBROUTINE trace_init

  SUBROUTINE trace_finalize (pe)
    INTEGER, INTENT(in) :: pe

#if defined(__xlC__) && defined(__hpm__)
    call f_hpmterminate(pe)
#endif

  END SUBROUTINE trace_finalize

  SUBROUTINE trace_start(name, number)
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER,          INTENT(in) :: number

#if defined(__SX__) && defined(_FTRACE)
    call ftrace_region_begin(TRIM(name))
#endif
#if defined(__xlC__) && defined(__hpm__)
    CALL f_hpmstart(number, TRIM(name))
#endif

  END SUBROUTINE trace_start

  SUBROUTINE trace_stop(name, number)
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER,          INTENT(in) :: number

#if defined(__SX__) && defined(_FTRACE)
    call ftrace_region_end(TRIM(name))
#endif
#if defined(__xlC__) && defined(__hpm__)
    CALL f_hpmstop(number)
#endif

  END SUBROUTINE trace_stop

END MODULE mo_profile
