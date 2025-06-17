!>
!! @file test_redist_p2p.c
!!
!! @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!
!
!  Keywords:
!  Maintainer: Jörg Behrens <behrens@dkrz.de>
!              Moritz Hanke <hanke@dkrz.de>
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
PROGRAM test_redist_p2p_f
  USE iso_c_binding, ONLY: c_ptr, c_loc
  USE mpi
  USE yaxt, ONLY: xt_int_kind, xt_xmap, xt_idxlist, xt_redist, xt_offset_ext, &
       xi => xt_int_kind, xt_int_mpidt, xt_initialize, xt_finalize, &
       xt_idxvec_new, xt_idxlist_delete, &
       xt_redist_p2p_new, xt_redist_p2p_off_new, xt_redist_p2p_ext_new, &
       xt_redist_s_exchange, xt_redist_s_exchange1, &
       xt_redist_delete,  xt_redist_get_mpi_comm, &
       xt_xmap_all2all_new, xt_xmap_delete
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  IMPLICIT NONE

  ! init mpi
  CALL init_mpi

  CALL xt_initialize(mpi_comm_world)

  ! offset-free test:
  ! source index list
  CALL test_without_offsets
  CALL test_with_offsets


  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)

  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE test_without_offsets
    INTEGER, PARAMETER :: src_num_indices = 14, dst_num_indices = 13, &
         num_arr = 1
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(src_num_indices) &
         = (/  5_xi, 67_xi,  4_xi,  5_xi, 13_xi, &
         &     9_xi,  2_xi,  1_xi,  0_xi, 96_xi, &
         &    13_xi, 12_xi,  1_xi,  3_xi /), &
         dst_index_list(dst_num_indices) = &
         & (/  5_xi,  4_xi,  3_xi, 96_xi,  1_xi, &
         &     5_xi,  4_xi,  5_xi,  4_xi,  3_xi, &
         &    13_xi,  2_xi,  1_xi /)
    INTEGER :: i
    DOUBLE PRECISION, TARGET :: src_data(src_num_indices) = &
         (/ (DBLE(i), i=0,13) /)
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_num_indices) &
         = (/ 0.0d0,  2.0d0, 13.0d0,  9.0d0,  7.0d0, &
         &    0.0d0,  2.0d0,  0.0d0,  2.0d0, 13.0d0, &
         &    4.0d0,  6.0d0,  7.0d0 /)

    DOUBLE PRECISION, TARGET :: dst_data(dst_num_indices)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist
    TYPE(c_ptr) :: src_data_p(num_arr), dst_data_p(num_arr)

    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices)

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices)

    ! xmap
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)

    ! redist_p2p
    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    ! test communicator of redist
    IF (.NOT. test_communicator(xt_redist_get_mpi_comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_MPI_Comm", &
         __FILE__, &
         __LINE__)

    ! test exchange
    dst_data(:) = -1.0d0

    src_data_p(1) = C_LOC(src_data)
    dst_data_p(1) = C_LOC(dst_data)

    CALL xt_redist_s_exchange(redist, num_arr, src_data_p, dst_data_p)

    IF (ANY(ref_dst_data(:) /= dst_data(:))) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    ! clean up
    CALL xt_redist_delete(redist)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE test_without_offsets

  SUBROUTINE test_with_offsets
    ! source index list
    INTEGER, PARAMETER :: src_num = 14, dst_num = 13
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(src_num) = &
         (/  5_xi, 67_xi,  4_xi,  5_xi, 13_xi, &
         &   9_xi,  2_xi,  1_xi,  0_xi, 96_xi, &
         &  13_xi, 12_xi,  1_xi,  3_xi /), &
         dst_index_list(dst_num) = &
         (/  5_xi,  4_xi,  3_xi, 96_xi,  1_xi, &
         &   5_xi,  4_xi,  5_xi,  4_xi,  3_xi, &
         &  13_xi,  2_xi,  1_xi /)
    INTEGER :: i
    INTEGER, PARAMETER :: src_pos(src_num) = (/ (i, i = 0, src_num - 1) /), &
         dst_pos(dst_num) = (/ ( dst_num - i, i = 1, dst_num ) /)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist
    DOUBLE PRECISION, TARGET :: src_data(src_num) = &
         (/ ( DBLE(i), i = 0, 13 ) /)
    DOUBLE PRECISION, TARGET :: dst_data(dst_num)
    TYPE(c_ptr) :: src_data_p, dst_data_p
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_num) = &
         (/  0.0d0,  2.0d0, 13.0d0,  9.0d0,  7.0d0, &
         &   0.0d0,  2.0d0,  0.0d0,  2.0d0, 13.0d0, &
         &   4.0d0,  6.0d0,  7.0d0 /)


    src_idxlist = xt_idxvec_new(src_index_list)

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num)

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)

    ! redist_p2p with offsets
    redist = xt_redist_p2p_off_new(xmap, src_pos, dst_pos, mpi_double_precision)

    ! test communicator of redist
    IF (.NOT. test_communicator(xt_redist_get_mpi_comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_MPI_Comm", &
         __FILE__, &
         __LINE__)

    ! test exchange
    dst_data(:) = -1.0d0

    src_data_p = C_LOC(src_data)
    dst_data_p = C_LOC(dst_data)

    CALL xt_redist_s_exchange1(redist, src_data_p, dst_data_p)


    IF (ANY(ref_dst_data(dst_pos(:) + 1) /= dst_data(:))) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    ! clean up
    CALL xt_redist_delete(redist)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE test_with_offsets

  SUBROUTINE test_offset_extents
    ! source/destination index lists
    INTEGER, PARAMETER :: src_num = 14, dst_num = 13
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(src_num) = &
         (/  5_xi, 67_xi,  4_xi,  5_xi, 13_xi, &
         &   9_xi,  2_xi,  1_xi,  0_xi, 96_xi, &
         &  13_xi, 12_xi,  1_xi,  3_xi /), &
         dst_index_list(dst_num) = &
         (/  5_xi,  4_xi,  3_xi, 96_xi,  1_xi, &
         &   5_xi,  4_xi,  5_xi,  4_xi,  3_xi, &
         &  13_xi, 2_xi, 1_xi /)
    INTEGER :: i
    INTEGER(xt_int_kind), TARGET :: src_data(src_num) &
         = (/ (INT(i, xi), i = 0, 13) /), dst_data(dst_num)
    INTEGER(xt_int_kind), PARAMETER :: ref_dst_data(dst_num) = &
         (/  7_xi,  6_xi,  4_xi, 13_xi,  2_xi, &
         &   0_xi,  2_xi,  0_xi,  7_xi,  9_xi, &
         &  13_xi,  2_xi,  0_xi /)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist
    TYPE(xt_offset_ext), PARAMETER :: &
         src_pos(1) = (/ xt_offset_ext(0, src_num, 1) /), &
         dst_pos(1) = (/ xt_offset_ext(dst_num - 1, dst_num, -1) /)
    TYPE(c_ptr) :: src_data_p, dst_data_p

    src_idxlist = xt_idxvec_new(src_index_list)
    dst_idxlist = xt_idxvec_new(dst_index_list)

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)

    ! redist_p2p with extents of offsets
    redist = xt_redist_p2p_ext_new(xmap, &
         src_pos, dst_pos, xt_int_mpidt)
    ! test communicator of redist
    IF (.NOT. test_communicator(xt_redist_get_MPI_Comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_MPI_Comm", &
         __FILE__, &
         __LINE__)

    dst_data(:) = -1_xi

    ! test exchange
    src_data_p = C_LOC(src_data)
    dst_data_p = C_LOC(dst_data)

    CALL xt_redist_s_exchange1(redist, src_data_p, dst_data_p)


    IF (ANY(ref_dst_data(:) /= dst_data(:))) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    ! clean up
    CALL xt_redist_delete(redist)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE test_offset_extents

  FUNCTION test_communicator(comm1, comm2) RESULT(congruent)
    INTEGER, INTENT(in) :: comm1, comm2
    LOGICAL :: congruent

    INTEGER :: ierror, rcode

    CALL mpi_comm_compare(comm1, comm2, rcode, ierror)
    congruent = ((rcode == mpi_ident) .OR. (rcode == mpi_congruent))
  END FUNCTION test_communicator
END PROGRAM test_redist_p2p_f
