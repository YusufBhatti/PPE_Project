!>
!! @file xt_ut.f90
!!
!! @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
!!
!!                                 Thomas Jahns <jahns@dkrz.de>
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
!  Keywords:
!  Maintainer: Jörg Behrens <behrens@dkrz.de>
!              Thomas Jahns <jahns@dkrz.de>
!  URL: https://doc.redmine.dkrz.de/yaxt/html/
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are  permitted provided that the following conditions are
!  met:
!
!  Redistributions of source code must retain the above copyright notice,
!  this list of conditions and the following disclaimer.
!
!  Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!
!  Neither the name of the DKRZ GmbH nor the names of its contributors
!  may be used to endorse or promote products derived from this software
!  without specific prior written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!

MODULE xt_ut
  !
  ! unitrans interface for accessing yaxt
  !

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char, c_int, &
       c_long, c_short, c_long_long, c_ptr, c_loc

  ! \todo to be resolved later(dependency problem): USE yaxt, ONLY : xt_int_kind

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ut_init, ut_init_decomposition, &
       ut_init_oneway_transposition_template, &
       ut_init_transposition, ut_transpose, &
       ut_destroy_transposition_template, &
       ut_abort, ut_destroy_decomposition, ut_destroy_transposition, &
       ut_finalize

  ! \todo to be replaced by use yaxt....:
  INTEGER, PARAMETER, PUBLIC :: xt_int_kind   = XT_INT_FC_KIND

  !PUBLIC :: xt_int_kind

  INTEGER, PARAMETER :: inflate_inner = 1
  INTEGER, PARAMETER :: inflate_outer = 2

  INTEGER, PUBLIC, PARAMETER :: comm_forward  = 1
  INTEGER, PUBLIC, PARAMETER :: comm_backward = 2

  INTEGER, PUBLIC, PARAMETER :: ut_mode_dt_p2p        = 1
  INTEGER, PUBLIC, PARAMETER :: ut_mode_dt_alltoall   = 2
  INTEGER, PUBLIC, PARAMETER :: ut_mode_pack_p2p      = 3
  INTEGER, PUBLIC, PARAMETER :: ut_mode_pack_alltoall = 4

  INTERFACE ut_init_decomposition
    MODULE PROCEDURE ut_init_decomposition_1d
  END INTERFACE ut_init_decomposition


  INTERFACE ut_init_transposition
    MODULE PROCEDURE ut_init_transposition_simple
    MODULE PROCEDURE ut_init_transposition_with_offsets
  END INTERFACE ut_init_transposition

  INTERFACE ut_transpose
    MODULE PROCEDURE ut_transpose_int
  END INTERFACE ut_transpose

  INTERFACE
    SUBROUTINE xt_ut_abort(msg, source, line) BIND(C, name='xt_ut_abort')
      IMPORT:: C_CHAR, C_INT
      IMPLICIT NONE
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: msg
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: source
      INTEGER(C_INT), INTENT(in) :: line
    END SUBROUTINE xt_ut_abort
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_ut_init(decomp_size, comm_tmpl_size, comm_size, debug_lvl, &
         mode, idebug_unit) BIND(C, name='xt_ut_init')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: decomp_size, comm_tmpl_size, &
           comm_size, debug_lvl, mode, idebug_unit
    END SUBROUTINE xt_ut_init
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_ut_finalize() BIND(C, name='xt_ut_finalize')
    END SUBROUTINE xt_ut_finalize
  END INTERFACE

  INTERFACE
    INTEGER(C_INT) FUNCTION xt_ut_init_decomposition_1d(idx_vec, idx_vec_n) &
         BIND(C, name='xt_ut_init_decomposition_1d')
      IMPORT:: C_INT
      IMPORT:: xt_int_kind
      IMPLICIT NONE
      INTEGER(xt_int_kind), DIMENSION(*), INTENT(in) :: idx_vec
      INTEGER(c_int), VALUE, INTENT(in) :: idx_vec_n
    END FUNCTION xt_ut_init_decomposition_1d
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_ut_destroy_decomposition(handle) &
         BIND(C, name='xt_ut_destroy_decomposition')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: handle
    END SUBROUTINE xt_ut_destroy_decomposition
  END INTERFACE

  INTERFACE
    INTEGER(C_INT) FUNCTION xt_ut_init_oneway_transposition_template(&
         decomp_handle_in, decomp_handle_out, mpi_world, check_unique) &
         & BIND(C, name='xt_ut_init_oneway_transposition_template')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: decomp_handle_in, &
           decomp_handle_out, mpi_world, check_unique
    END FUNCTION xt_ut_init_oneway_transposition_template
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_ut_destroy_transposition_template(handle) &
         BIND(C, name='xt_ut_destroy_transposition_template')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: handle
    END SUBROUTINE xt_ut_destroy_transposition_template
  END INTERFACE

  INTERFACE
    INTEGER FUNCTION xt_ut_init_transposition_simple(itemplate, datatype) &
         BIND(C, name='xt_ut_init_transposition_simple')
      IMPLICIT NONE
      INTEGER, VALUE, INTENT(in) :: itemplate
      INTEGER, VALUE, INTENT(in) :: datatype
    END FUNCTION xt_ut_init_transposition_simple
  END INTERFACE

  INTERFACE
    INTEGER FUNCTION xt_ut_init_transposition(itemplate, offset_in, &
         offset_in_size, offset_out, offset_out_size, &
         & datatype)  BIND(C, name='xt_ut_init_transposition')
      IMPLICIT NONE
      INTEGER, VALUE, INTENT(in) :: itemplate, offset_in_size, offset_out_size
      INTEGER, DIMENSION(*), INTENT(in) :: offset_in, offset_out
      INTEGER, VALUE, INTENT(in) :: datatype
    END FUNCTION xt_ut_init_transposition
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_ut_destroy_transposition(handle) &
         BIND(C, name='xt_ut_destroy_transposition')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER(C_INT), VALUE, INTENT(in) :: handle
    END SUBROUTINE xt_ut_destroy_transposition
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_ut_transpose(pt_in, transposition_handle, &
         direction, pt_out) BIND(C, name='xt_ut_transpose')
      IMPORT:: C_INT, c_ptr
      IMPLICIT NONE
      TYPE(c_ptr) :: pt_in
      TYPE(c_ptr) :: pt_out
      INTEGER(C_INT), VALUE, INTENT(in) :: transposition_handle
      INTEGER(C_INT), VALUE, INTENT(in) :: direction
    END SUBROUTINE xt_ut_transpose
  END INTERFACE


  LOGICAL, PARAMETER :: debug = .TRUE.

CONTAINS

  SUBROUTINE ut_abort(msg, source, line)
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line

    CALL xt_ut_abort(TRIM(msg)//c_null_char, source, line)

  END SUBROUTINE ut_abort

  SUBROUTINE ut_finalize()
    CALL xt_ut_finalize()
  END SUBROUTINE ut_finalize

  SUBROUTINE ut_init(decomp_size, comm_tmpl_size, comm_size, debug_lvl, &
       mode, debug_unit)
    INTEGER, INTENT(in) :: decomp_size
    INTEGER, INTENT(in) :: comm_tmpl_size
    INTEGER, INTENT(in) :: comm_size
    INTEGER, INTENT(in) :: debug_lvl
    INTEGER, INTENT(in) :: mode
    INTEGER, INTENT(in) :: debug_unit

    CALL xt_ut_init(decomp_size, comm_tmpl_size, comm_size, debug_lvl, &
         mode, debug_unit)

  END SUBROUTINE ut_init

  SUBROUTINE ut_init_decomposition_1d(myindex, global_size, handle)
    INTEGER, INTENT(in)  :: myindex(:)
    INTEGER, INTENT(in)  :: global_size
    INTEGER, INTENT(out) :: handle

    IF (HUGE(1_xt_int_kind) < HUGE(myindex)) THEN
      IF (ANY(myindex > HUGE(1_xt_int_kind)) &
           .OR. ANY(myindex < -HUGE(1_xt_int_kind))) &
           CALL ut_abort('ut_init_decomposition_1d: &
           &index value not supported', &
           __FILE__, &
           __LINE__)
    END IF
    IF (HUGE(SIZE(myindex)) > HUGE(1_c_int)) THEN
      IF (SIZE(myindex) > HUGE(1_c_int)) &
           CALL ut_abort('ut_init_decomposition_1: &
           &array size unsupported', &
           __FILE__, &
           __LINE__)
    END IF
    handle = xt_ut_init_decomposition_1d(idx_vec=INT(myindex, xt_int_kind), &
         idx_vec_n=INT(SIZE(myindex), c_int))

  END SUBROUTINE ut_init_decomposition_1d

  SUBROUTINE ut_destroy_decomposition(handle)
    INTEGER, INTENT(in) :: handle

    CALL xt_ut_destroy_decomposition(handle)

  END SUBROUTINE ut_destroy_decomposition

  SUBROUTINE ut_init_oneway_transposition_template(decomp_handle_in, &
       decomp_handle_out, mpi_world, comm_tmpl_handle, check_unique)
    INTEGER, INTENT(in)  :: decomp_handle_in
    INTEGER, INTENT(in)  :: decomp_handle_out
    INTEGER, INTENT(in)  :: mpi_world
    INTEGER, INTENT(out) :: comm_tmpl_handle
    LOGICAL, OPTIONAL, INTENT(in) :: check_unique

    INTEGER :: icheck_unique

    icheck_unique = 0
    IF (PRESENT(check_unique)) THEN
      IF (check_unique) icheck_unique = 1
    ENDIF
    comm_tmpl_handle = xt_ut_init_oneway_transposition_template(&
         decomp_handle_in, decomp_handle_out, mpi_world, icheck_unique)

  END SUBROUTINE ut_init_oneway_transposition_template

  SUBROUTINE ut_destroy_transposition_template(handle)
    INTEGER, INTENT(in) :: handle

    CALL xt_ut_destroy_transposition_template(handle)

  END SUBROUTINE ut_destroy_transposition_template

  SUBROUTINE ut_init_transposition_simple(comm_template_handle, &
       datatype, comm_handle)
    INTEGER, INTENT(in)  :: comm_template_handle
    INTEGER, INTENT(in)  :: datatype
    INTEGER, INTENT(out) :: comm_handle

    comm_handle = xt_ut_init_transposition_simple(comm_template_handle, &
         datatype)

  END SUBROUTINE ut_init_transposition_simple

  SUBROUTINE ut_init_transposition_with_offsets(comm_template_handle, &
       offset_in, offset_out, datatype_in, datatype_out, comm_handle)
    INTEGER, INTENT(in)  :: comm_template_handle
    INTEGER, INTENT(in)  :: offset_in(:)
    INTEGER, INTENT(in)  :: offset_out(:)
    INTEGER, INTENT(in)  :: datatype_in
    INTEGER, INTENT(in)  :: datatype_out
    INTEGER, INTENT(out) :: comm_handle

    INTEGER :: datatype

    IF (datatype_in /= datatype_out) THEN
      CALL ut_abort('ut_init_transposition: &
           &(datatype_in /= datatype_out) not supported', &
           __FILE__, &
           __LINE__)
    ENDIF

    datatype = datatype_in

    comm_handle = xt_ut_init_transposition(comm_template_handle, &
         & offset_in, SIZE(offset_in), &
         & offset_out, SIZE(offset_out), &
         & datatype)

  END SUBROUTINE ut_init_transposition_with_offsets

  SUBROUTINE ut_destroy_transposition(handle)
    INTEGER, INTENT(in) :: handle

    CALL xt_ut_destroy_transposition(handle)

  END SUBROUTINE ut_destroy_transposition

  SUBROUTINE ut_transpose_int(field_in, transposition_handle, direction, &
       field_out)
    INTEGER, POINTER     :: field_in, field_out
    INTEGER, INTENT(in) :: transposition_handle
    INTEGER, INTENT(in) :: direction

    TYPE(c_ptr) :: pt_in, pt_out

    pt_in = C_LOC(field_in)
    pt_out = C_LOC(field_out)
    CALL xt_ut_transpose(pt_in, transposition_handle, direction, pt_out)

  END SUBROUTINE ut_transpose_int

End MODULE xt_ut
