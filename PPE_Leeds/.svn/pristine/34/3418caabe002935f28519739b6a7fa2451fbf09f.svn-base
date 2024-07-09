!>
!! @file test_xmap_all2all_fail_f.f90
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
PROGRAM test_xmap_all2all_fail
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, &
       xt_idxlist, xt_idxlist_delete, xt_idxvec_new, &
       xt_xmap, xt_xmap_delete, xt_xmap_all2all_new, &
       xt_xmap_get_num_destinations, xt_xmap_get_num_sources, &
       xt_xmap_get_destination_ranks, xt_xmap_get_source_ranks
  IMPLICIT NONE
  INTEGER, PARAMETER :: xi = xt_int_kind
  INTEGER :: my_rank, ierror
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  CALL mpi_comm_rank(mpi_comm_world, my_rank, ierror)

  CALL test_xmap1

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE shift_idx(idx, offset)
    INTEGER(xt_int_kind), INTENT(inout) :: idx(:)
    INTEGER(xt_int_kind), INTENT(in) :: offset
    INTEGER :: i
    DO i = 1, SIZE(idx)
      idx(i) = idx(i) + INT(my_rank, xi) * offset
    END DO
  END SUBROUTINE shift_idx

  SUBROUTINE test_xmap(src_index_list, dst_index_list)
    INTEGER(xt_int_kind), INTENT(in) :: src_index_list(:), dst_index_list(:)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    INTEGER :: rank(1)
    src_idxlist = xt_idxvec_new(src_index_list)
    dst_idxlist = xt_idxvec_new(dst_index_list)

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    IF (xt_xmap_get_num_destinations(xmap) /= 1) &
         CALL test_abort("error in xmap construction", &
         __FILE__, &
         __LINE__)

    IF (xt_xmap_get_num_sources(xmap) /= 1) &
         CALL test_abort("error in xt_xmap_get_num_sources", &
         __FILE__, &
         __LINE__)
    CALL xt_xmap_get_destination_ranks(xmap, rank)
    IF (rank(1) /= my_rank) &
         CALL test_abort("error in xt_xmap_get_destination_ranks", &
         __FILE__, &
         __LINE__)

    CALL xt_xmap_get_source_ranks(xmap, rank)
    IF (rank(1) /= my_rank) &
         CALL test_abort("error in xt_xmap_get_source_ranks", &
         __FILE__, &
         __LINE__)
    CALL xt_xmap_delete(xmap)
  END SUBROUTINE test_xmap

  SUBROUTINE test_xmap1
    INTEGER(xt_int_kind) :: i
    INTEGER(xt_int_kind), PARAMETER :: num_src_idx = 7, num_dst_idx = 7
    INTEGER(xt_int_kind) :: src_index_list(num_src_idx), &
         dst_index_list(num_dst_idx)
    DO i = 1_xi, num_src_idx
      src_index_list(i) = i
    END DO
    CALL shift_idx(src_index_list, num_src_idx)
    DO i = 1_xi, num_dst_idx
      dst_index_list(i) = num_dst_idx - i + 2_xi
    END DO
    CALL shift_idx(dst_index_list, num_src_idx)
    ! note: this should fail because dst/src indices don't match
    CALL test_xmap(src_index_list, dst_index_list)
  END SUBROUTINE test_xmap1

END PROGRAM test_xmap_all2all_fail
