!>
!! @file xt_redist_f.f90
!! @brief xt_redist-related procedures of Fortran interface
!!
!! @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

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
MODULE xt_redist_base
  USE xt_core, ONLY: xt_abort, xt_mpi_fint_kind, i2, i4, i8
  USE xt_xmap_abstract, ONLY: xt_xmap
  USE iso_c_binding, ONLY: c_int, c_null_ptr, c_ptr, c_associated
  USE xt_mpi, ONLY: mpi_address_kind
  IMPLICIT NONE
  PRIVATE
  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_redist
    PRIVATE
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_redist

  TYPE, BIND(c), PUBLIC :: xt_offset_ext
    INTEGER(c_int) :: start, size, stride
  END TYPE xt_offset_ext

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_redist_f2c(redist) BIND(c, name='xt_redist_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_redist
      IMPLICIT NONE
      TYPE(xt_redist), INTENT(in) :: redist
      TYPE(c_ptr) :: p
    END FUNCTION xt_redist_f2c
  END INTERFACE

  INTERFACE xt_redist_delete
    MODULE PROCEDURE xt_redist_delete_1
    MODULE PROCEDURE xt_redist_delete_a1d
  END INTERFACE xt_redist_delete

  INTERFACE
    SUBROUTINE xt_redist_delete_c(redist) &
         BIND(C, name='xt_redist_delete')
      IMPORT :: c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE :: redist
    END SUBROUTINE xt_redist_delete_c
  END INTERFACE

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_redist_is_null
  END INTERFACE xt_is_null

  INTERFACE xt_redist_s_exchange
    MODULE PROCEDURE xt_redist_s_exchange1
    MODULE PROCEDURE xt_redist_s_exchange_a1d
    MODULE PROCEDURE xt_redist_s_exchange_i2_a1d
    MODULE PROCEDURE xt_redist_s_exchange_i4_a1d
    MODULE PROCEDURE xt_redist_s_exchange_i8_a1d
  END INTERFACE xt_redist_s_exchange

  INTERFACE
    SUBROUTINE xt_redist_s_exchange_c(redist, num_ptr, src_data_cptr, &
         dst_data_cptr) BIND(C, name='xt_redist_s_exchange')
      IMPORT:: c_ptr, c_int
      TYPE(c_ptr), VALUE, INTENT(in) :: redist
      INTEGER(c_int), VALUE, INTENT(in) :: num_ptr
      TYPE(c_ptr) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    END SUBROUTINE xt_redist_s_exchange_c

    FUNCTION xt_redist_get_mpi_comm(redist) &
         BIND(c, name='xt_redist_get_mpi_comm_c2f') RESULT(comm)
      IMPORT :: xt_redist, xt_mpi_fint_kind
      TYPE(xt_redist), INTENT(in) :: redist
      INTEGER(xt_mpi_fint_kind) :: comm
    END FUNCTION xt_redist_get_mpi_comm

    FUNCTION xt_redist_p2p_ext_new_c2f(xmap, num_src_ext, src_extents, &
         num_dst_ext, dst_extents, datatype) &
         BIND(c, name='xt_redist_p2p_ext_new_c2f') RESULT(redist)
      IMPORT :: c_int, c_ptr, xt_offset_ext, xt_mpi_fint_kind, xt_xmap
      TYPE(xt_xmap), INTENT(in) :: xmap
      INTEGER(c_int), VALUE, INTENT(in) :: num_src_ext, num_dst_ext
      TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
           dst_extents(num_dst_ext)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
      TYPE(c_ptr) :: redist
    END FUNCTION xt_redist_p2p_ext_new_c2f
  END INTERFACE

  INTERFACE xt_redist_p2p_ext_new
    MODULE PROCEDURE xt_redist_p2p_ext_new_i2_a1d_i2_a1d
    MODULE PROCEDURE xt_redist_p2p_ext_new_i4_a1d_i4_a1d
    MODULE PROCEDURE xt_redist_p2p_ext_new_i8_a1d_i8_a1d
    MODULE PROCEDURE xt_redist_p2p_ext_new_a1d_a1d
  END INTERFACE xt_redist_p2p_ext_new

  PUBLIC :: xt_redist_c2f, xt_redist_f2c, xt_is_null, &
       xt_redist_delete, xt_redist_s_exchange1, xt_redist_s_exchange, &
       xt_redist_p2p_new, xt_redist_p2p_off_new, &
       xt_redist_p2p_blocks_new, xt_redist_p2p_blocks_off_new, &
       xt_redist_collection_static_new, xt_redist_collection_new, &
       xt_redist_repeat_new, xt_redist_get_mpi_comm, xt_redist_p2p_ext_new
CONTAINS

  FUNCTION xt_redist_is_null(redist) RESULT(p)
    TYPE(xt_redist), INTENT(in) :: redist
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(redist%cptr)
  END FUNCTION xt_redist_is_null

  FUNCTION xt_redist_c2f(redist) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: redist
    TYPE(xt_redist) :: p
    p%cptr = redist
  END FUNCTION xt_redist_c2f

  SUBROUTINE xt_redist_delete_1(redist)
    TYPE(xt_redist), INTENT(inout) :: redist
    CALL xt_redist_delete_c(redist%cptr)
    redist%cptr = c_null_ptr
  END SUBROUTINE xt_redist_delete_1

  SUBROUTINE xt_redist_delete_a1d(redists)
    TYPE(xt_redist), INTENT(inout) :: redists(:)
    INTEGER :: i, n
    n = SIZE(redists)
    DO i = 1, n
      CALL xt_redist_delete_c(redists(i)%cptr)
      redists(i)%cptr = c_null_ptr
    END DO
  END SUBROUTINE xt_redist_delete_a1d

  SUBROUTINE xt_redist_s_exchange1(redist, src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    TYPE(c_ptr) :: src_data_cptr, dst_data_cptr
    INTERFACE
      SUBROUTINE xt_redist_s_exchange1_c(redist, src_data_cptr, dst_data_cptr) &
           BIND(C, name='xt_redist_s_exchange1')
        IMPORT:: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: redist
        TYPE(c_ptr), VALUE :: src_data_cptr, dst_data_cptr
      END SUBROUTINE xt_redist_s_exchange1_c
    END INTERFACE
    CALL xt_redist_s_exchange1_c(xt_redist_f2c(redist), src_data_cptr, &
         dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange1

  SUBROUTINE xt_redist_s_exchange_a1d(redist, src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    TYPE(c_ptr) :: src_data_cptr(:), dst_data_cptr(:)
    INTEGER :: n
    n = SIZE(src_data_cptr)
    IF (n /= SIZE(dst_data_cptr) .OR. n > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of pointers", &
         __FILE__, &
         __LINE__)
    CALL xt_redist_s_exchange_c(xt_redist_f2c(redist), INT(n, c_int), &
         src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_a1d

  SUBROUTINE xt_redist_s_exchange_i2_a1d(redist, num_ptr, &
       src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(i2), INTENT(in) :: num_ptr
    TYPE(c_ptr) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    IF (num_ptr < 0_i2) &
         CALL xt_abort("invalid number of pointers", &
         __FILE__, &
         __LINE__)
    CALL xt_redist_s_exchange_c(xt_redist_f2c(redist), INT(num_ptr, c_int), &
         src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_i2_a1d

  SUBROUTINE xt_redist_s_exchange_i4_a1d(redist, num_ptr, &
       src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(i4), INTENT(in) :: num_ptr
    TYPE(c_ptr) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    IF (num_ptr < 0_i4 .OR. num_ptr > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of pointers", &
         __FILE__, &
         __LINE__)
    CALL xt_redist_s_exchange_c(xt_redist_f2c(redist), INT(num_ptr, c_int), &
         src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_i4_a1d

  SUBROUTINE xt_redist_s_exchange_i8_a1d(redist, num_ptr, &
       src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(i8), INTENT(in) :: num_ptr
    TYPE(c_ptr) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    IF (num_ptr < 0_i8 .OR. num_ptr > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of pointers", &
         __FILE__, &
         __LINE__)
    CALL xt_redist_s_exchange_c(xt_redist_f2c(redist), INT(num_ptr, c_int), &
         src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_i8_a1d

  FUNCTION xt_redist_p2p_new(xmap, datatype) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: res

    INTERFACE
      FUNCTION xt_redist_p2p_new_f(xmap, datatype) &
           BIND(C, name='xt_redist_p2p_new_f') RESULT(res_ptr)
        IMPORT:: xt_xmap, xt_redist, c_int, c_ptr, xt_mpi_fint_kind
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_new_f
    END INTERFACE

    res = xt_redist_c2f(xt_redist_p2p_new_f(xmap, datatype))

  END FUNCTION xt_redist_p2p_new

  FUNCTION xt_redist_p2p_off_new(xmap, src_offsets, dst_offsets, datatype) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: src_offsets(*)
    INTEGER, INTENT(in) :: dst_offsets(*)
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: res

    INTERFACE
      FUNCTION xt_redist_p2p_off_new_f(xmap, src_offsets, dst_offsets, &
           datatype) BIND(C, name='xt_redist_p2p_off_new_f') RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_redist, c_int, c_ptr, xt_mpi_fint_kind
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(xt_mpi_fint_kind), INTENT(in) :: src_offsets(*)
        INTEGER(xt_mpi_fint_kind), INTENT(in) :: dst_offsets(*)
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_off_new_f
    END INTERFACE

    res = xt_redist_c2f(&
         xt_redist_p2p_off_new_f(xmap, src_offsets, dst_offsets, datatype))

  END FUNCTION xt_redist_p2p_off_new

  FUNCTION xt_redist_p2p_blocks_new(xmap, src_block_sizes, src_block_num, &
       &                                  dst_block_sizes, dst_block_num, &
       &                                  datatype) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
    INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_blocks_new_f(xmap, &
           &                              src_block_sizes, src_block_num, &
           &                              dst_block_sizes, dst_block_num, &
           &                              datatype) &
           BIND(C, name='xt_redist_p2p_blocks_new_f') RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_mpi_fint_kind, xt_redist, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
        INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_new_f
    END INTERFACE

    res = xt_redist_c2f(&
         xt_redist_p2p_blocks_new_f(xmap, &
         &                          src_block_sizes, src_block_num, &
         &                          dst_block_sizes, dst_block_num, &
         &                          datatype))
  END FUNCTION xt_redist_p2p_blocks_new

  FUNCTION xt_redist_p2p_blocks_off_new(xmap, src_block_offsets, &
       src_block_sizes, src_block_num, &
       dst_block_offsets, dst_block_sizes, dst_block_num, &
       datatype) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(in) :: src_block_offsets(*)
    INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
    INTEGER(c_int), INTENT(in) :: dst_block_offsets(*)
    INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_blocks_off_new_f(xmap, src_block_offsets, &
           src_block_sizes, src_block_num, &
           dst_block_offsets, dst_block_sizes, dst_block_num, &
           datatype) BIND(C, name='xt_redist_p2p_blocks_off_new_f') &
           RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_redist, xt_mpi_fint_kind, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), INTENT(in) :: src_block_offsets(*)
        INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
        INTEGER(c_int), INTENT(in) :: dst_block_offsets(*)
        INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_off_new_f
    END INTERFACE

    res = xt_redist_c2f(&
         xt_redist_p2p_blocks_off_new_f(xmap, src_block_offsets, &
           src_block_sizes, src_block_num, &
           dst_block_offsets, dst_block_sizes, dst_block_num, &
           datatype))

  END FUNCTION xt_redist_p2p_blocks_off_new

  FUNCTION xt_redist_collection_static_new(redists, num_redists, &
       src_displacements, dst_displacements, comm) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(*)
    INTEGER, INTENT(in) :: num_redists, comm
    INTEGER(mpi_address_kind), INTENT(in) :: src_displacements(*)
    INTEGER(mpi_address_kind), INTENT(in) :: dst_displacements(*)
    TYPE(xt_redist) :: res

    INTERFACE
      FUNCTION xt_redist_collection_static_new_f(redists_f, num_redists, &
           src_displacements, dst_displacements, comm_f) &
           BIND(C, name='xt_redist_collection_static_new_f') RESULT(res)
        IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind

        IMPLICIT NONE
        TYPE(xt_redist), INTENT(in) :: redists_f(*)
        INTEGER(xt_mpi_fint_kind), VALUE :: num_redists
        INTEGER(mpi_address_kind), INTENT(in) :: src_displacements(*)
        INTEGER(mpi_address_kind), INTENT(in) :: dst_displacements(*)
        INTEGER(xt_mpi_fint_kind), VALUE :: comm_f
        TYPE(c_ptr) :: res
      END FUNCTION xt_redist_collection_static_new_f
    END INTERFACE

    res = xt_redist_c2f(xt_redist_collection_static_new_f(redists, &
         num_redists, src_displacements, dst_displacements, comm))

  END FUNCTION xt_redist_collection_static_new

  FUNCTION xt_redist_collection_new(redists, num_redists, cache_size, comm) &
       RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(*)
    INTEGER, INTENT(in) :: num_redists, cache_size, comm
    TYPE(xt_redist) :: res

    INTERFACE
      FUNCTION xt_redist_collection_new_f(redists_f, num_redists, cache_size, &
           comm_f) BIND(C, name='xt_redist_collection_new_f') RESULT(res)
        IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind
        IMPLICIT NONE
        TYPE(xt_redist), INTENT(in) :: redists_f(*)
        INTEGER(xt_mpi_fint_kind), VALUE :: num_redists, cache_size, comm_f
        TYPE(c_ptr) :: res
      END FUNCTION xt_redist_collection_new_f
    END INTERFACE

    res = xt_redist_c2f(xt_redist_collection_new_f(redists, num_redists, &
         cache_size, comm))

  END FUNCTION xt_redist_collection_new

  FUNCTION xt_redist_repeat_new(redist, src_extent, dst_extent, &
                                num_repetitions, displacements) &
       RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(mpi_address_kind), INTENT(in) :: src_extent, dst_extent
    INTEGER, INTENT(in) :: num_repetitions
    INTEGER(c_int), INTENT(in) :: displacements(*)
    TYPE(xt_redist) :: res

    INTERFACE
      FUNCTION xt_redist_repeat_new_c(redist_f, src_extent, dst_extent, &
           num_repetitions, displacements) &
           BIND(C, name='xt_redist_repeat_new') RESULT(res)
        IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in)  :: redist_f
        INTEGER(mpi_address_kind), VALUE :: src_extent, dst_extent
        INTEGER(c_int), VALUE :: num_repetitions
        INTEGER(c_int), INTENT(in) :: displacements(*)
        TYPE(c_ptr) :: res
      END FUNCTION xt_redist_repeat_new_c
    END INTERFACE

    res = xt_redist_c2f(xt_redist_repeat_new_c(xt_redist_f2c(redist), &
         src_extent, dst_extent, num_repetitions, displacements))

  END FUNCTION xt_redist_repeat_new

  FUNCTION xt_redist_p2p_ext_new_i2_a1d_i2_a1d(xmap, num_src_ext, src_extents, &
       num_dst_ext, dst_extents, datatype) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i2), INTENT(in) :: num_src_ext, num_dst_ext
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
         dst_extents(num_dst_ext)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: redist
    IF (num_src_ext < 0_i2 .OR. num_src_ext > HUGE(1_c_int) &
         .OR. num_dst_ext < 0_i2 .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", &
         __FILE__, &
         __LINE__)
    redist = xt_redist_c2f(xt_redist_p2p_ext_new_c2f(xmap, &
         INT(num_src_ext, c_int), src_extents, &
         INT(num_dst_ext, c_int), dst_extents, datatype))
  END FUNCTION xt_redist_p2p_ext_new_i2_a1d_i2_a1d

  FUNCTION xt_redist_p2p_ext_new_i4_a1d_i4_a1d(xmap, num_src_ext, src_extents, &
       num_dst_ext, dst_extents, datatype) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i4), INTENT(in) :: num_src_ext, num_dst_ext
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
         dst_extents(num_dst_ext)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: redist
    IF (num_src_ext < 0_i4 .OR. num_src_ext > HUGE(1_c_int) &
         .OR. num_dst_ext < 0_i4 .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", &
         __FILE__, &
         __LINE__)
    redist = xt_redist_c2f(xt_redist_p2p_ext_new_c2f(xmap, &
         INT(num_src_ext, c_int), src_extents, &
         INT(num_dst_ext, c_int), dst_extents, datatype))
  END FUNCTION xt_redist_p2p_ext_new_i4_a1d_i4_a1d

  FUNCTION xt_redist_p2p_ext_new_i8_a1d_i8_a1d(xmap, num_src_ext, src_extents, &
       num_dst_ext, dst_extents, datatype) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i8), INTENT(in) :: num_src_ext, num_dst_ext
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
         dst_extents(num_dst_ext)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: redist
    IF (num_src_ext < 0_i8 .OR. num_src_ext > HUGE(1_c_int) &
         .OR. num_dst_ext < 0_i8 .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", &
         __FILE__, &
         __LINE__)
    redist = xt_redist_c2f(xt_redist_p2p_ext_new_c2f(xmap, &
         INT(num_src_ext, c_int), src_extents, &
         INT(num_dst_ext, c_int), dst_extents, datatype))
  END FUNCTION xt_redist_p2p_ext_new_i8_a1d_i8_a1d

  FUNCTION xt_redist_p2p_ext_new_a1d_a1d(xmap, src_extents, dst_extents, &
       datatype) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(:), &
         dst_extents(:)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: redist
    INTEGER :: num_src_ext, num_dst_ext

    num_src_ext = SIZE(src_extents)
    num_dst_ext = SIZE(dst_extents)

    IF (num_src_ext > HUGE(1_c_int) .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", &
         __FILE__, &
         __LINE__)
    redist = xt_redist_c2f(xt_redist_p2p_ext_new_c2f(xmap, &
         INT(num_src_ext, c_int), src_extents, &
         INT(num_dst_ext, c_int), dst_extents, datatype))
  END FUNCTION xt_redist_p2p_ext_new_a1d_a1d

END MODULE xt_redist_base
