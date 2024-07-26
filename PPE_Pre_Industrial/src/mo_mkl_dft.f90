!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
#ifdef MKL_DFT

! Load Intel MKL definitions for DFT
#include "mkl_dfti.f90"

MODULE mo_mkl_dft

  USE mo_kind, ONLY: dp

  USE mkl_dft_type 
  USE mkl_dfti

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_mkl_dft
  PUBLIC :: cleanup_mkl_dft
  PUBLIC :: mkl_dft

  TYPE(DFTI_DESCRIPTOR), POINTER :: my_desc_handle => NULL()

  INTEGER :: status

CONTAINS

  SUBROUTINE init_mkl_dft(inc, jump, n, lot)

    INTEGER, INTENT(in) :: inc, jump, n, lot

    INTEGER :: inc_arr(2)

    inc_arr(1) = 0
    inc_arr(2) = inc

    status = DftiCreateDescriptor( my_desc_handle, DFTI_DOUBLE, &
         DFTI_REAL, 1, n )

    status = DftiSetValue(my_desc_handle, DFTI_NUMBER_OF_TRANSFORMS, lot)
    status = DftiSetValue(my_desc_handle, DFTI_INPUT_DISTANCE, jump)
    status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_DISTANCE, jump)
    status = DftiSetValue(my_desc_handle, DFTI_INPUT_STRIDES, inc_arr)
    status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, inc_arr)
    status = DftiSetValue(my_desc_handle, DFTI_FORWARD_SCALE, 1.0_dp/REAL(n,dp))

    status = DftiCommitDescriptor(my_desc_handle)

  END SUBROUTINE init_mkl_dft

  SUBROUTINE cleanup_mkl_dft

    IF (ASSOCIATED(my_desc_handle)) THEN
        status = DftiFreeDescriptor(my_desc_handle)
        NULLIFY(my_desc_handle)
    ENDIF

  END SUBROUTINE cleanup_mkl_dft

  SUBROUTINE mkl_dft(a, inc, jump, n, lot, isign)

    INTEGER, INTENT(in) :: inc, isign, jump, lot, n

    REAL(dp), INTENT(inout) :: a(*)

    INTEGER :: inc1(2), inc2(2)
    INTEGER :: jump1,jump2, lot1


    IF (.NOT. ASSOCIATED(my_desc_handle)) THEN
        CALL init_mkl_dft(inc, jump, n, lot)
    ELSE
        status = DftiGetValue(my_desc_handle, DFTI_NUMBER_OF_TRANSFORMS, lot1)
        status = DftiGetValue(my_desc_handle, DFTI_INPUT_DISTANCE, jump1)
        status = DftiGetValue(my_desc_handle, DFTI_OUTPUT_DISTANCE, jump2)
        status = DftiGetValue(my_desc_handle, DFTI_INPUT_STRIDES, inc1)
        status = DftiGetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, inc2)

        IF ( lot1 /= lot    .OR.                    &
             jump1 /= jump  .OR. jump2 /= jump .OR. &
             inc1(2) /= inc .OR. inc2(2) /= inc ) THEN
          CALL cleanup_mkl_dft
          CALL init_mkl_dft(inc, jump, n, lot)
        ENDIF
    ENDIF

    IF (isign == -1) THEN
        status = DftiComputeForward(my_desc_handle, a)
    ELSE
        status = DftiComputeBackward(my_desc_handle, a)
    ENDIF

  END SUBROUTINE mkl_dft

END MODULE mo_mkl_dft

#else

MODULE mo_mkl_dft
END MODULE mo_mkl_dft

#endif
