!>
!! @file test_xmap_all2all_parallel_f.f90
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
MODULE test_xmap_common_parallel
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, xt_stripe, &
       xt_idxlist, xt_idxlist_delete, xt_idxvec_new, &
       xt_idxstripes_new, xt_idxempty_new, &
       xt_xmap, xt_xmap_delete, &
       xt_xmap_get_num_destinations, xt_xmap_get_num_sources, &
       xt_xmap_get_destination_ranks, xt_xmap_get_source_ranks
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: xi = xt_int_kind
  INTEGER :: comm_rank, comm_size
  PUBLIC :: xmap_parallel_test_main
CONTAINS
  SUBROUTINE xmap_parallel_test_main(xmap_new)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER :: ierror
    CALL init_mpi
    CALL xt_initialize(mpi_comm_world)
    CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
    CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
    IF (comm_size > HUGE(xi)) &
         CALL test_abort("number of ranks exceeds test limit", &
         __FILE__, &
         __LINE__)

    CALL test_allgather_analog(xmap_new)
    IF (comm_size > 2) CALL test_ring_1d(xmap_new)
    IF (comm_size == 2) CALL test_pair(xmap_new)
    IF (comm_size > 1) CALL test_ping_pong(xmap_new)

    IF (test_err_count() /= 0) &
         CALL test_abort("non-zero error count!", &
         __FILE__, &
         __LINE__)
    CALL xt_finalize
    CALL finish_mpi
  END SUBROUTINE xmap_parallel_test_main

  SUBROUTINE test_allgather_analog(xmap_new)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER(xt_int_kind) :: i, src_index_list(1)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_stripe) :: comm_rank_stripe(1)
    INTEGER, ALLOCATABLE :: ranks(:)
    src_index_list(1) = INT(comm_rank, xi)
    src_idxlist = xt_idxvec_new(src_index_list)

    comm_rank_stripe(1) = xt_stripe(0, 1, comm_size)
    dst_idxlist = xt_idxstripes_new(comm_rank_stripe)

    xmap = xmap_new(src_idxlist, dst_idxlist, mpi_comm_world)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    IF (xt_xmap_get_num_destinations(xmap) /= INT(comm_size, xi)) &
         CALL test_abort("error in xmap construction", &
         __FILE__, &
         __LINE__)

    IF (xt_xmap_get_num_sources(xmap) /= INT(comm_size, xi)) &
         CALL test_abort("error in xt_xmap_get_num_sources", &
         __FILE__, &
         __LINE__)

    ALLOCATE(ranks(comm_size))

    CALL xt_xmap_get_destination_ranks(xmap, ranks)
    IF (ANY(ranks /= (/ (i, i=0,INT(comm_size-1, xi)) /))) &
         CALL test_abort("error in xt_xmap_get_destination_ranks", &
         __FILE__, &
         __LINE__)

    CALL xt_xmap_get_source_ranks(xmap, ranks)
    IF (ANY(ranks /= (/ (i, i=0,INT(comm_size-1, xi)) /))) &
         CALL test_abort("error in xt_xmap_get_source_ranks", &
         __FILE__, &
         __LINE__)
    DEALLOCATE(ranks)
    CALL xt_xmap_delete(xmap)
  END SUBROUTINE test_allgather_analog

  SUBROUTINE test_ring_1d(xmap_new)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    ! test in which each process talks WITH two other processes
    INTEGER(xt_int_kind) :: src_index_list(1), dst_index_list(2), temp
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(Xt_xmap) :: xmap
    INTEGER :: ranks(2)

    src_index_list(1) = INT(comm_rank, xi)
    src_idxlist = xt_idxvec_new(src_index_list)

    ! destination index list
    dst_index_list(1) = INT(MOD(comm_rank + comm_size - 1, comm_size), xi)
    dst_index_list(2) = INT(MOD(comm_rank             + 1, comm_size), xi)
    IF (dst_index_list(1) > dst_index_list(2)) THEN
      temp = dst_index_list(1)
      dst_index_list(1) = dst_index_list(2)
      dst_index_list(2) = temp
    END IF
    dst_idxlist = xt_idxvec_new(dst_index_list, 2)

    ! test of exchange map
    xmap = xmap_new(src_idxlist, dst_idxlist, mpi_comm_world)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist);

    ! test results
    IF (xt_xmap_get_num_destinations(xmap) /= 2) &
         CALL test_abort("error in xt_xmap_get_num_destinations", &
         __FILE__, &
         __LINE__)

    IF (xt_xmap_get_num_sources(xmap) /= 2) &
         CALL test_abort("error in xt_xmap_get_num_sources", &
         __FILE__, &
         __LINE__)

    CALL xt_xmap_get_destination_ranks(xmap, ranks)
    IF (ANY(ranks /= dst_index_list)) &
         CALL test_abort("error in xt_xmap_get_destination_ranks", &
         __FILE__, &
         __LINE__)

    CALL xt_xmap_get_source_ranks(xmap, ranks)
    IF (ANY(ranks /= dst_index_list)) &
         CALL test_abort("error in xt_xmap_get_source_ranks", &
         __FILE__, &
         __LINE__)

    ! clean up
    CALL xt_xmap_delete(xmap)

  END SUBROUTINE test_ring_1d

  SUBROUTINE test_pair(xmap_new)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    !src_index_list(index, rank)
    INTEGER(xt_int_kind) :: i, j, k
#ifdef __xlC__
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(20, 0:1) = RESHAPE((/ &
         &  1_xi, 2_xi, 3_xi, 4_xi, 5_xi, &
         &  9_xi, 10_xi, 11_xi, 12_xi, 13_xi, &
         & 17_xi, 18_xi, 19_xi, 20_xi, 21_xi, &
         & 25_xi, 26_xi, 27_xi, 28_xi, 29_xi, &
         &  4_xi,  5_xi,  6_xi,  7_xi,  8_xi, &
         & 12_xi, 13_xi, 14_xi, 15_xi, 16_xi, &
         & 20_xi, 21_xi, 22_xi, 23_xi, 24_xi, &
         & 28_xi, 29_xi, 30_xi, 31_xi, 32_xi /),  &
         (/ 20, 2 /))
#else
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(20, 0:1) = RESHAPE((/ &
         (((i + j * 8_xi + k * 3_xi, i = 1_xi, 5_xi), j = 0_xi,3_xi), &
         k = 0_xi,1_xi) /), (/ 20, 2 /))
#endif
    ! dst_index_list(index,rank)
    INTEGER(xt_int_kind), PARAMETER :: dst_index_list(20, 0:1) = RESHAPE((/ &
         10_xi, 15_xi, 14_xi, 13_xi, 12_xi, &
         15_xi, 10_xi, 11_xi, 12_xi, 13_xi, &
         23_xi, 18_xi, 19_xi, 20_xi, 21_xi, &
         31_xi, 26_xi, 27_xi, 28_xi, 29_xi, &
         13_xi, 12_xi, 11_xi, 10_xi, 15_xi, &
         12_xi, 13_xi, 14_xi, 15_xi, 10_xi, &
         20_xi, 21_xi, 22_xi, 23_xi, 18_xi, &
         28_xi, 29_xi, 30_xi, 31_xi, 26_xi /),  &
         (/ 20, 2 /))
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    INTEGER :: ranks(2)

    src_idxlist = xt_idxvec_new(src_index_list(:, comm_rank))

    ! destination index list
    dst_idxlist = xt_idxvec_new(dst_index_list(:, comm_rank))

    ! test of exchange map
    xmap = xmap_new(src_idxlist, dst_idxlist, mpi_comm_world)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    ! test results
    IF (xt_xmap_get_num_destinations(xmap) /= 2) &
         CALL test_abort("error in xt_xmap_get_num_destinations", &
         __FILE__, &
         __LINE__)

    IF (xt_xmap_get_num_sources(xmap) /= 2) &
         CALL test_abort("error in xt_xmap_get_num_sources", &
         __FILE__, &
         __LINE__)

    CALL xt_xmap_get_destination_ranks(xmap, ranks)
    IF (ranks(1) /= 0 .OR. ranks(2) /= 1) &
         CALL test_abort("error in xt_xmap_get_destination_ranks", &
         __FILE__, &
         __LINE__)

    CALL xt_xmap_get_source_ranks(xmap, ranks)
    IF (ranks(1) /= 0 .OR. ranks(2) /= 1) &
         CALL test_abort("error in xt_xmap_get_source_ranks", &
         __FILE__, &
         __LINE__)

    ! clean up
    CALL xt_xmap_delete(xmap)
  END SUBROUTINE test_pair

  SUBROUTINE test_ping_pong(xmap_new)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER(xt_int_kind) :: index_list(5) = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi /)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    CHARACTER(len=80) :: msg
    INTEGER :: expect, dst_rank(1), src_rank(1)

    IF (comm_rank == 0) THEN
      src_idxlist = xt_idxvec_new(index_list, 5)
    ELSE
      src_idxlist = xt_idxempty_new()
    END IF


    IF (comm_rank == comm_size-1) THEN
      dst_idxlist = xt_idxvec_new(index_list, 5)
    ELSE
      dst_idxlist = xt_idxempty_new()
    END IF

    ! test of exchange map

    xmap = xmap_new(src_idxlist, dst_idxlist, mpi_comm_world)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)


    WRITE (msg, '(a,i0,a)') "error in xt_xmap_get_num_destinations (rank == ", &
         comm_rank, ")"
    ! test results
    expect = MERGE(0, 1, comm_rank /= 0)
    IF (xt_xmap_get_num_destinations(xmap) /= expect) &
         CALL test_abort(TRIM(msg), &
         __FILE__, &
         __LINE__)

    expect = MERGE(0, 1, comm_rank /= comm_size - 1)
    IF (xt_xmap_get_num_sources(xmap) /= expect) &
         CALL test_abort(msg, &
         __FILE__, &
         __LINE__)

    IF (comm_rank == 0) THEN
      CALL xt_xmap_get_destination_ranks(xmap, dst_rank)
      IF (dst_rank(1) /= comm_size - 1) &
           CALL test_abort("error in xt_xmap_get_destination_ranks", &
           __FILE__, &
           __LINE__)
    ELSE IF (comm_rank == comm_size) THEN
      CALL xt_xmap_get_source_ranks(xmap, src_rank)
      IF (src_rank(1) /= 0) &
           CALL test_abort("error in xt_xmap_get_source_ranks", &
           __FILE__, &
           __LINE__)
    END IF

    ! clean up
    CALL xt_xmap_delete(xmap)
  END SUBROUTINE test_ping_pong
END MODULE test_xmap_common_parallel
