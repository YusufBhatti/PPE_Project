!>
!! @file xt_xmap_f.f90
!! @brief Fortran interface to yaxt xmap declarations
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
MODULE xt_xmap_abstract
  USE iso_c_binding, ONLY: c_int, c_ptr, c_null_ptr, &
       c_associated, c_f_pointer
  USE xt_core, ONLY: xt_mpi_fint_kind, xt_pos_ext
  USE xt_idxlist_abstract, ONLY: xt_idxlist
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xt_xmap_c2f, xt_xmap_f2c, xt_is_null
  PUBLIC :: xt_xmap_delete, xt_xmap_get_num_destinations, &
       xt_xmap_get_num_sources, xt_xmap_get_destination_ranks, &
       xt_xmap_get_source_ranks
  PUBLIC :: xt_xmap_all2all_new, xt_xmap_dist_dir_new
  PUBLIC :: xt_xmap_get_in_iterator, xt_xmap_get_out_iterator, &
       xt_xmap_iterator_next, xt_xmap_iterator_get_rank, &
       xt_xmap_iterator_get_transfer_pos, &
       xt_xmap_iterator_get_transfer_pos_ext, &
       xt_xmap_iterator_get_num_transfer_pos, &
       xt_xmap_iterator_get_num_transfer_pos_ext, &
       xt_xmap_iterator_delete


  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_xmap
    PRIVATE
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_xmap

  TYPE, BIND(c), PUBLIC :: xt_xmap_iter
    PRIVATE
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_xmap_iter

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_xmap_f2c(xmap) BIND(c, name='xt_xmap_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_xmap
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(in) :: xmap
      TYPE(c_ptr) :: p
    END FUNCTION xt_xmap_f2c

    SUBROUTINE xt_xmap_delete_c(xmap) BIND(C, name='xt_xmap_delete')
      IMPORT :: c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE :: xmap
    END SUBROUTINE xt_xmap_delete_c

  END INTERFACE

  INTERFACE xt_xmap_delete
    MODULE PROCEDURE xt_xmap_delete_1
    MODULE PROCEDURE xt_xmap_delete_a1d
  END INTERFACE xt_xmap_delete

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_xmap_is_null
    MODULE PROCEDURE xt_xmap_iterator_is_null
  END INTERFACE xt_is_null

  INTERFACE
    FUNCTION xt_xmap_iterator_get_num_transfer_pos_c(iter) RESULT(num) &
         BIND(c, name='xt_xmap_iterator_get_num_transfer_pos')
      IMPORT :: c_int, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iter
      INTEGER(c_int) :: num
    END FUNCTION xt_xmap_iterator_get_num_transfer_pos_c

    FUNCTION xt_xmap_iterator_get_num_transfer_pos_ext_c(iter) RESULT(num) &
         BIND(c, name='xt_xmap_iterator_get_num_transfer_pos_ext')
      IMPORT :: c_int, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iter
      INTEGER(c_int) :: num
    END FUNCTION xt_xmap_iterator_get_num_transfer_pos_ext_c
  END INTERFACE
CONTAINS

  FUNCTION xt_xmap_is_null(xmap) RESULT(p)
    TYPE(xt_xmap), INTENT(in) :: xmap
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(xmap%cptr)
  END FUNCTION xt_xmap_is_null


  FUNCTION xt_xmap_c2f(xmap) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: xmap
    TYPE(xt_xmap) :: p
    p%cptr = xmap
  END FUNCTION xt_xmap_c2f

  SUBROUTINE xt_xmap_delete_1(xmap)
    TYPE(xt_xmap), INTENT(inout) :: xmap
    CALL xt_xmap_delete_c(xt_xmap_f2c(xmap))
    xmap%cptr = c_null_ptr
  END SUBROUTINE xt_xmap_delete_1

  SUBROUTINE xt_xmap_delete_a1d(xmaps)
    TYPE(xt_xmap), INTENT(inout) :: xmaps(:)
    INTEGER :: i, n
    n = SIZE(xmaps)
    DO i = 1, n
      CALL xt_xmap_delete_c(xt_xmap_f2c(xmaps(i)))
      xmaps(i)%cptr = c_null_ptr
    END DO
  END SUBROUTINE xt_xmap_delete_a1d

  FUNCTION xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, VALUE, INTENT(in) :: comm
    TYPE(xt_xmap) :: res

    INTERFACE
      FUNCTION xt_xmap_all2all_new_f(src_idxlist, dst_idxlist, comm) &
           BIND(C, name='xt_xmap_all2all_new_f') RESULT(res_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, c_ptr
        IMPLICIT NONE
        TYPE(Xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(Xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_xmap_all2all_new_f
    END INTERFACE

    res = xt_xmap_c2f(xt_xmap_all2all_new_f(src_idxlist, dst_idxlist, comm))
  END FUNCTION xt_xmap_all2all_new

  FUNCTION xt_xmap_dist_dir_new(src_idxlist, dst_idxlist, comm) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, VALUE, INTENT(in) :: comm
    TYPE(xt_xmap) :: res

    INTERFACE
      FUNCTION xt_xmap_dist_dir_new_f(src_idxlist, dst_idxlist, comm) &
           BIND(C, name='xt_xmap_dist_dir_new_f') RESULT(res_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, c_ptr
        IMPLICIT NONE
        TYPE(Xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(Xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_xmap_dist_dir_new_f
    END INTERFACE

    res = xt_xmap_c2f(xt_xmap_dist_dir_new_f(src_idxlist, dst_idxlist, comm))
  END FUNCTION xt_xmap_dist_dir_new

  FUNCTION xt_xmap_get_num_destinations(xmap) RESULT(num)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER :: num
    INTERFACE
      FUNCTION xt_xmap_get_num_destinations_c(xmap) RESULT(num) &
           BIND(c, name='xt_xmap_get_num_destinations')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int) :: num
      END FUNCTION xt_xmap_get_num_destinations_c
    END INTERFACE
    num = INT(xt_xmap_get_num_destinations_c(xmap%cptr))
  END FUNCTION xt_xmap_get_num_destinations

  FUNCTION xt_xmap_get_num_sources(xmap) RESULT(num)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER :: num
    INTERFACE
      FUNCTION xt_xmap_get_num_sources_c(xmap) RESULT(num) &
           BIND(c, name='xt_xmap_get_num_sources')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int) :: num
      END FUNCTION xt_xmap_get_num_sources_c
    END INTERFACE
    num = INT(xt_xmap_get_num_sources_c(xmap%cptr))
  END FUNCTION xt_xmap_get_num_sources

  SUBROUTINE xt_xmap_get_destination_ranks(xmap, ranks)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(out) :: ranks(*)
    INTERFACE
      SUBROUTINE xt_xmap_get_destination_ranks_c(xmap, ranks) &
           BIND(c, name='xt_xmap_get_destination_ranks')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int), INTENT(out) :: ranks(*)
      END SUBROUTINE xt_xmap_get_destination_ranks_c
    END INTERFACE
    CALL xt_xmap_get_destination_ranks_c(xmap%cptr, ranks)
  END SUBROUTINE xt_xmap_get_destination_ranks

  SUBROUTINE xt_xmap_get_source_ranks(xmap, ranks)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(out) :: ranks(*)
    INTERFACE
      SUBROUTINE xt_xmap_get_source_ranks_c(xmap, ranks) &
           BIND(c, name='xt_xmap_get_source_ranks')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int), INTENT(out) :: ranks(*)
      END SUBROUTINE xt_xmap_get_source_ranks_c
    END INTERFACE
    CALL xt_xmap_get_source_ranks_c(xmap%cptr, ranks)
  END SUBROUTINE xt_xmap_get_source_ranks

  FUNCTION xt_xmap_get_out_iterator(xmap) RESULT(iter)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(xt_xmap_iter) :: iter
    INTERFACE
      FUNCTION xt_xmap_get_out_iterator_c(xmap) RESULT(cptr) &
           BIND(c, name='xt_xmap_get_out_iterator')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        TYPE(c_ptr) :: cptr
      END FUNCTION xt_xmap_get_out_iterator_c
    END INTERFACE
    iter%cptr = xt_xmap_get_out_iterator_c(xmap%cptr)
  END FUNCTION xt_xmap_get_out_iterator

  FUNCTION xt_xmap_get_in_iterator(xmap) RESULT(iter)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(xt_xmap_iter) :: iter
    INTERFACE
      FUNCTION xt_xmap_get_in_iterator_c(xmap) RESULT(cptr) &
           BIND(c, name='xt_xmap_get_in_iterator')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        TYPE(c_ptr) :: cptr
      END FUNCTION xt_xmap_get_in_iterator_c
    END INTERFACE
    iter%cptr = xt_xmap_get_in_iterator_c(xmap%cptr)
  END FUNCTION xt_xmap_get_in_iterator

  FUNCTION xt_xmap_iterator_is_null(iter) RESULT(p)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(iter%cptr)
  END FUNCTION xt_xmap_iterator_is_null

  FUNCTION xt_xmap_iterator_next(iter) RESULT(avail)
    TYPE(xt_xmap_iter), INTENT(inout) :: iter
    LOGICAL :: avail
    INTERFACE
      FUNCTION xt_xmap_iterator_next_c(iter) RESULT(avail) &
           BIND(c, name='xt_xmap_iterator_next')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        INTEGER(c_int) :: avail
      END FUNCTION xt_xmap_iterator_next_c
    END INTERFACE
    avail = xt_xmap_iterator_next_c(iter%cptr) /= 0
  END FUNCTION xt_xmap_iterator_next

  FUNCTION xt_xmap_iterator_get_rank(iter) RESULT(rank)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    INTEGER :: rank
    INTERFACE
      FUNCTION xt_xmap_iterator_get_rank_c(iter) RESULT(rank) &
           BIND(c, name='xt_xmap_iterator_get_rank')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        INTEGER(c_int) :: rank
      END FUNCTION xt_xmap_iterator_get_rank_c
    END INTERFACE
    rank = INT(xt_xmap_iterator_get_rank_c(iter%cptr))
  END FUNCTION xt_xmap_iterator_get_rank

  !> note: result is read-only
  FUNCTION xt_xmap_iterator_get_transfer_pos(iter) RESULT(transfer_pos)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    INTEGER(c_int), POINTER :: transfer_pos(:)

    INTERFACE
      FUNCTION xt_xmap_iterator_get_transfer_pos_c(iter) RESULT(transfer_pos) &
           BIND(c, name='xt_xmap_iterator_get_transfer_pos')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        TYPE(c_ptr) :: transfer_pos
      END FUNCTION xt_xmap_iterator_get_transfer_pos_c
    END INTERFACE
    INTEGER :: n(1)
    TYPE(c_ptr) :: transfer_pos_cptr
    NULLIFY(transfer_pos)
    n(1) = INT(xt_xmap_iterator_get_num_transfer_pos_c(iter%cptr))
    transfer_pos_cptr = xt_xmap_iterator_get_transfer_pos_c(iter%cptr)
    CALL C_F_POINTER(transfer_pos_cptr, transfer_pos, n)
  END FUNCTION xt_xmap_iterator_get_transfer_pos

  FUNCTION xt_xmap_iterator_get_num_transfer_pos(iter) RESULT(num)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    INTEGER :: num
    num = INT(xt_xmap_iterator_get_num_transfer_pos_c(iter%cptr))
  END FUNCTION xt_xmap_iterator_get_num_transfer_pos

  !> note: result is read-only
  FUNCTION xt_xmap_iterator_get_transfer_pos_ext(iter) RESULT(transfer_pos_ext)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    TYPE(xt_pos_ext), POINTER :: transfer_pos_ext(:)

    INTERFACE
      FUNCTION xt_xmap_iterator_get_transfer_pos_ext_c(iter) &
           RESULT(transfer_pos_ext) &
           BIND(c, name='xt_xmap_iterator_get_transfer_pos_ext')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        TYPE(c_ptr) :: transfer_pos_ext
      END FUNCTION xt_xmap_iterator_get_transfer_pos_ext_c
    END INTERFACE
    INTEGER :: n(1)
    TYPE(c_ptr) :: transfer_pos_ext_cptr
    NULLIFY(transfer_pos_ext)
    n(1) = INT(xt_xmap_iterator_get_num_transfer_pos_ext_c(iter%cptr))
    transfer_pos_ext_cptr = xt_xmap_iterator_get_transfer_pos_ext_c(iter%cptr)
    CALL C_F_POINTER(transfer_pos_ext_cptr, transfer_pos_ext, n)
  END FUNCTION xt_xmap_iterator_get_transfer_pos_ext

  FUNCTION xt_xmap_iterator_get_num_transfer_pos_ext(iter) RESULT(num)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    INTEGER :: num
    num = INT(xt_xmap_iterator_get_num_transfer_pos_ext_c(iter%cptr))
  END FUNCTION xt_xmap_iterator_get_num_transfer_pos_ext

  SUBROUTINE xt_xmap_iterator_delete(iter)
    TYPE(xt_xmap_iter), INTENT(inout) :: iter
    INTERFACE
      SUBROUTINE xt_xmap_iterator_delete_c(iter) &
           BIND(c, name='xt_xmap_iterator_delete')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
      END SUBROUTINE xt_xmap_iterator_delete_c
    END INTERFACE
    CALL xt_xmap_iterator_delete_c(iter%cptr)
    iter%cptr = c_null_ptr
  END SUBROUTINE xt_xmap_iterator_delete
END MODULE xt_xmap_abstract
