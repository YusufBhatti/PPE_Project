!>
!> @file xt_core_f.f90
!> @brief Fortran interface to yaxt core declarations
!>
!> @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
!>                                 Moritz Hanke <hanke@dkrz.de>
!>                                 Thomas Jahns <jahns@dkrz.de>
!>
!> @author Jörg Behrens <behrens@dkrz.de>
!>         Moritz Hanke <hanke@dkrz.de>
!>         Thomas Jahns <jahns@dkrz.de>
!>

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://doc.redmine.dkrz.de/yaxt/html/
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
MODULE xt_core
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char, c_int, &
       c_long, c_short, c_long_long
  USE xt_mpi, ONLY: XT_INT_FC_MPIDT
  IMPLICIT NONE
  PRIVATE
  INTEGER, PUBLIC, PARAMETER :: xt_int_kind   = XT_INT_FC_KIND
  INTEGER, PUBLIC, PARAMETER :: pi2 = 4
  INTEGER, PUBLIC, PARAMETER :: pi4 = 9
  INTEGER, PUBLIC, PARAMETER :: pi8 = 14
  INTEGER, PUBLIC, PARAMETER :: i2 = SELECTED_INT_KIND(pi2)
  INTEGER, PUBLIC, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)
  INTEGER, PUBLIC, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)
  PUBLIC :: xt_initialize, xt_finalize, xt_abort, xt_get_default_comm, char
  PUBLIC :: xt_initialized, xt_finalized
  PUBLIC :: xt_slice_c_loc
  PUBLIC :: OPERATOR(==), OPERATOR(/=)

  INTEGER, PUBLIC, PARAMETER :: xt_mpi_fint_kind = XT_MPI_FINT_FC_KIND
  INTEGER, PUBLIC, PARAMETER :: xt_int_mpidt = XT_INT_FC_MPIDT
  INTEGER(xt_int_kind), PARAMETER :: dummy = 0_xt_int_kind
  !> number of decimal places needed to print any variable of type
  !! INTEGER(xt_int_kind)
  INTEGER, PUBLIC, PARAMETER :: xt_int_dec_len &
       = CEILING(1.0 + REAL(DIGITS(dummy)) * LOG10(REAL(RADIX(dummy))))
  CHARACTER(9), PARAMETER :: xt_stripe_tag = 'xt_stripe'
  !> maximal length of string xt_stripe(a, b, c)
  INTEGER, PUBLIC, PARAMETER :: xt_stripe2s_len &
       = LEN(xt_stripe_tag) + 2 + 4 + 3 * xt_int_dec_len

  TYPE, BIND(C), PUBLIC :: xt_stripe
    INTEGER(xt_int_kind) :: start
    INTEGER(xt_int_kind) :: stride
    INTEGER(c_int) :: nstrides
  END TYPE xt_stripe

  TYPE, BIND(C), PUBLIC :: xt_bounds
    INTEGER(xt_int_kind) :: start, size
  END TYPE xt_bounds

  !> describes range of positions starting with start up to start + size - 1
  !! i.e. [start,start+size) if size is positive and down to start + size + 1
  !! i.e. (start+size,start] if size is negative
  TYPE, BIND(c), PUBLIC :: xt_pos_ext
    INTEGER(c_int) :: start, size
  END TYPE xt_pos_ext

  INTERFACE

    FUNCTION xt_get_default_comm() RESULT(comm) &
         BIND(c, name='xt_get_default_comm_f')
      IMPORT :: xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(xt_mpi_fint_kind) :: comm
    END FUNCTION xt_get_default_comm

    SUBROUTINE xt_initialize(default_comm) BIND(C, name='xt_initialize_f')
      IMPORT:: c_int, xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(xt_mpi_fint_kind), INTENT(in) :: default_comm
    END SUBROUTINE xt_initialize

    SUBROUTINE xt_finalize() BIND(C, name='xt_finalize')
    END SUBROUTINE xt_finalize

  END INTERFACE

  INTERFACE xt_abort
    MODULE PROCEDURE xt_abort4
    MODULE PROCEDURE xt_abort3
  END INTERFACE xt_abort

  INTERFACE
    SUBROUTINE xt_abort_c(comm, msg, source, line) BIND(c, name='xt_abort_f')
      IMPORT :: c_char, c_int, xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in):: comm
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: msg
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: source
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: line
    END SUBROUTINE xt_abort_c
  END INTERFACE

  INTERFACE char
    MODULE PROCEDURE xt_stripe2char
  END INTERFACE char

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE xt_pos_ext_eq
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE xt_pos_ext_ne
  END INTERFACE OPERATOR(/=)

  EXTERNAL :: xt_slice_c_loc

CONTAINS

  SUBROUTINE xt_abort4(comm, msg, source, line)
    INTEGER, INTENT(in) :: comm
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    CALL xt_abort_c(comm, TRIM(msg)//c_null_char, &
         TRIM(source)//c_null_char, INT(line, c_int))
  END SUBROUTINE xt_abort4

  SUBROUTINE xt_abort3(msg, source, line)
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    CALL xt_abort_c(xt_get_default_comm(), TRIM(msg)//c_null_char, &
         TRIM(source)//c_null_char, INT(line, c_int))
  END SUBROUTINE xt_abort3

  ELEMENTAL FUNCTION xt_stripe2char(stripe) RESULT(str)
    CHARACTER(len=xt_stripe2s_len) :: str
    TYPE(xt_stripe), INTENT(in) :: stripe
    WRITE (str, '(2a,3(i0,a))') xt_stripe_tag, '(', stripe%start, ', ', &
         stripe%stride, ', ', stripe%nstrides, ')'
  END FUNCTION xt_stripe2char

  FUNCTION xt_initialized() RESULT(is_initialized)
    LOGICAL :: is_initialized
    INTERFACE
      FUNCTION xt_initialized_c() BIND(c, name='xt_initialized') &
           RESULT(is_initialized)
        IMPORT :: c_int
        INTEGER(c_int) :: is_initialized
      END FUNCTION xt_initialized_c
    END INTERFACE
    is_initialized = xt_initialized_c() /= 0
  END FUNCTION xt_initialized

  FUNCTION xt_finalized() RESULT(is_finalized)
    LOGICAL :: is_finalized
    INTERFACE
      FUNCTION xt_finalized_c() BIND(c, name='xt_finalized') &
           RESULT(is_finalized)
        IMPORT :: c_int
        INTEGER(c_int) :: is_finalized
      END FUNCTION xt_finalized_c
    END INTERFACE
    is_finalized = xt_finalized_c() /= 0
  END FUNCTION xt_finalized

  ELEMENTAL FUNCTION xt_pos_ext_eq(a, b) RESULT(p)
    TYPE(xt_pos_ext), INTENT(in) :: a, b
    LOGICAL :: p
    p = a%start == b%start .AND. (a%size == b%size &
         .OR. (ABS(a%size) == 1 .AND. ABS(a%size) == ABS(b%size)))
  END FUNCTION xt_pos_ext_eq

  ELEMENTAL FUNCTION xt_pos_ext_ne(a, b) RESULT(p)
    TYPE(xt_pos_ext), INTENT(in) :: a, b
    LOGICAL :: p
    p = a%start /= b%start .OR. (a%size /= b%size &
         .AND. .NOT. (ABS(a%size) == 1 .AND. ABS(a%size) == ABS(b%size)))
  END FUNCTION xt_pos_ext_ne

END MODULE xt_core
