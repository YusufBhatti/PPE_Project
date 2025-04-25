!>
!! @file test_idxstripes_f.f90
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
PROGRAM test_idxstripes_f
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE mpi
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count, &
       idxlist_pack_unpack_copy
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_bounds, xt_pos_ext, &
       xt_stripe, xt_idxlist, xt_idxlist_delete, xt_idxstripes_new, &
       xt_idxvec_from_stripes_new, xt_int_kind, xt_idxlist_copy, &
       xt_idxlist_get_index_stripes, xt_idxlist_get_intersection, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_bounding_box, OPERATOR(/=), &
       xt_idxlist_get_pos_exts_of_index_stripes, &
       xt_idxlist_get_num_indices, xt_idxvec_new, &
       xt_idxstripes_from_idxlist_new
  USE iso_c_binding, ONLY: c_int
  IMPLICIT NONE
  INTEGER, PARAMETER :: xi = xt_int_kind

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  CALL stripe_test_general1
  CALL stripe_test_general2
  CALL stripe_test_general3
  CALL stripe_test_general4
  CALL stripe_test_general5
  CALL test_intersection1
  CALL test_intersection2
  CALL test_intersection3
  CALL test_intersection4
  CALL test_intersection5
  CALL test_intersection6
  CALL test_intersection7
  CALL test_intersection8
  CALL test_intersection9
  CALL test_intersection10
  CALL test_intersection11
  CALL test_get_pos1
  CALL test_get_pos2
  CALL test_get_pos3
  CALL test_get_pos4
  CALL test_stripe_overlap
  CALL test_stripe_bb1
  CALL test_stripe_bb2
  CALL check_pos_ext1
  CALL check_pos_ext2
  CALL check_pos_ext3
  CALL check_pos_ext4
  CALL check_pos_ext5
  CALL check_pos_ext6
  CALL check_pos_ext7
  CALL check_pos_ext8
  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi

CONTAINS
  SUBROUTINE stripe_test_general(stripes, ref_indices)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    TYPE(xt_idxlist) :: idxstripes, idxvec
    INTEGER :: num_ext, num_unmatched, num_pos, i
    INTEGER(c_int) :: ext_size
    TYPE(xt_pos_ext), ALLOCATABLE :: pos_ext(:)

    idxstripes = xt_idxstripes_new(stripes, SIZE(stripes))
    CALL do_tests(idxstripes, ref_indices)

    num_unmatched = xt_idxlist_get_pos_exts_of_index_stripes(idxstripes, &
         SIZE(stripes), stripes, num_ext, pos_ext, .TRUE.)
    IF (num_unmatched /= 0) &
         CALL test_abort("stripes not found", &
         __FILE__, &
         __LINE__)

    num_pos = 0
    DO i = 1, num_ext
      ext_size = pos_ext(i)%size
      IF (num_pos /= pos_ext(i)%start) &
           CALL test_abort("position/start mismatch", &
           __FILE__, &
           __LINE__)
      num_pos = num_pos + ext_size
    END DO
    IF (num_pos /= xt_idxlist_get_num_indices(idxstripes)) &
         CALL test_abort("index list length/positions overlap mismatch", &
         __FILE__, &
         __LINE__)

    DEALLOCATE(pos_ext)
    CALL xt_idxlist_delete(idxstripes)

    ! test recreation of stripes from reference vector
    idxvec = xt_idxvec_new(ref_indices)
    idxstripes = xt_idxstripes_from_idxlist_new(idxvec)
    CALL check_idxlist(idxstripes, ref_indices)
    CALL xt_idxlist_delete(idxvec)
    CALL xt_idxlist_delete(idxstripes)
  END SUBROUTINE stripe_test_general

  SUBROUTINE stripe_test_general1
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 1, 5), xt_stripe(20, 1, 5) /);
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0_xi,  1_xi,  2_xi,  3_xi,  4_xi, &
         &   10_xi, 11_xi, 12_xi, 13_xi, 14_xi, &
         &   20_xi, 21_xi, 22_xi, 23_xi, 24_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general1

  SUBROUTINE stripe_test_general2
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 2, 5), xt_stripe(20, 3, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0_xi,  1_xi,  2_xi,  3_xi,  4_xi, &
         &   10_xi, 12_xi, 14_xi, 16_xi, 18_xi, &
         &   20_xi, 23_xi, 26_xi, 29_xi, 32_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general2

  SUBROUTINE stripe_test_general3
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, 6, 5), &
         xt_stripe(1, 3, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, 6_xi, 12_xi, 18_xi, 24_xi, &
         &    1_xi, 4_xi,  7_xi, 10_xi, 13_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general3

  SUBROUTINE stripe_test_general4
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, -1, 5), &
         xt_stripe(1, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, -1_xi, -2_xi, -3_xi, -4_xi, &
         &    1_xi,  2_xi,  3_xi,  4_xi,  5_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general4

  SUBROUTINE stripe_test_general5
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(9, -2, 5), &
         xt_stripe(0, 2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 9_xi,  7_xi,  5_xi,  3_xi,  1_xi, &
         &    0_xi,  2_xi,  4_xi,  6_xi,  8_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general5

  SUBROUTINE test_intersection(stripes_a, stripes_b, ref_indices_a, ref_indices_b)
    TYPE(xt_stripe), INTENT(in) :: stripes_a(:), stripes_b(:)
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices_a(:)
    INTEGER(xt_int_kind), OPTIONAL, INTENT(in) :: ref_indices_b(:)
    TYPE(xt_idxlist) :: idxstripes_a, idxstripes_b, intersection(2)

    idxstripes_a = xt_idxstripes_new(stripes_a)
    idxstripes_b = xt_idxstripes_new(stripes_b)
    intersection(1) = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b)
    intersection(2) = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a)
    CALL do_tests(intersection(1), ref_indices_a)
    IF (PRESENT(ref_indices_b)) THEN
      CALL do_tests(intersection(2), ref_indices_b)
    ELSE
      CALL do_tests(intersection(2), ref_indices_a)
    END IF
    CALL xt_idxlist_delete(intersection(2))
    CALL xt_idxlist_delete(intersection(1))
    CALL xt_idxlist_delete(idxstripes_a)
    CALL xt_idxlist_delete(idxstripes_b)
  END SUBROUTINE test_intersection

  SUBROUTINE test_intersection1
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 1, 4), &
         xt_stripe(6, 1, 4) /), &
         stripes_b(1) = (/ xt_stripe(1, 1, 8) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(6) &
         = (/ 1_xi, 2_xi, 3_xi, 6_xi, 7_xi, 8_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection1

  SUBROUTINE test_intersection2
    TYPE(xt_stripe), PARAMETER :: stripes_a(3) = (/ xt_stripe(0, 1, 4), &
         xt_stripe(6, 1, 4), xt_stripe(11, 1, 4) /), &
         stripes_b(2) = (/ xt_stripe(1, 1, 7), xt_stripe(9, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(9) &
         = (/ 1_xi, 2_xi, 3_xi, 6_xi, 7_xi, 9_xi, 11_xi, 12_xi, 13_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection2

  SUBROUTINE test_intersection3
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 1, 3), &
         xt_stripe(8, 1, 3) /), &
         stripes_b(2) = (/ xt_stripe(3, 1, 5), xt_stripe(11, 1, 3) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ -1_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices(1:0))
  END SUBROUTINE test_intersection3

  SUBROUTINE test_intersection4
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 1, 10) /), &
         stripes_b(2) = (/ xt_stripe(0, 2, 5), xt_stripe(9, -2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 8_xi, 9_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection4

  SUBROUTINE test_intersection5
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 3, 5), &
         xt_stripe(1, 7, 5) /), &
         stripes_b(2) = (/ xt_stripe(0, 2, 7), xt_stripe(24, -1, 10) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(6) &
         = (/ 0_xi, 6_xi, 8_xi, 12_xi, 15_xi, 22_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection5

  SUBROUTINE test_intersection6
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 1, 10) /), &
         stripes_b(2) = (/ xt_stripe(5, 1, 5), xt_stripe(4, -1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 8_xi, 9_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection6

  SUBROUTINE test_intersection7
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 1, 10) , &
            xt_stripe(20, 1, 5) /), &
         stripes_b(2) = (/ xt_stripe(3, 1, 5), xt_stripe(17, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(7) &
         = (/ 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 20_xi, 21_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection7

  SUBROUTINE test_intersection8
    TYPE(xt_stripe), PARAMETER :: stripes_a(10) = (/ xt_stripe(0, 1, 2), &
         xt_stripe(3, 1, 2), xt_stripe(5, 1, 2), xt_stripe(8, 1, 2), &
         xt_stripe(10, 1, 2), xt_stripe(14, 1, 2), xt_stripe(17, 1, 2), &
         xt_stripe(20, 1, 2), xt_stripe(23, 1, 2), xt_stripe(25, 1, 2) /), &
         stripes_b(5) = (/ xt_stripe(5, 1, 3), xt_stripe(8, 1, 2), &
         xt_stripe(19, 1, 1), xt_stripe(20, 1, 2), xt_stripe(30, 1, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(6) &
         = (/ 5_xi, 6_xi, 8_xi, 9_xi, 20_xi, 21_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection8

  SUBROUTINE test_intersection9
    TYPE(xt_stripe), PARAMETER :: stripes_a(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(1, 1, 5), xt_stripe(2, 1, 5) /), &
         stripes_b(1) = (/ xt_stripe(-2, 1, 10) /)
    INTEGER(xt_int_kind) :: i
    INTEGER(xt_int_kind), PARAMETER :: ref_indices_a(7) &
         = (/ (i, i=0_xi,6_xi) /), &
         ref_indices_b(15) = (/ 0_xi, 1_xi, 1_xi, 2_xi, 2_xi, 2_xi, 3_xi, &
         &                      3_xi, 3_xi, 4_xi, 4_xi, 4_xi, 5_xi, 5_xi, 6_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices_a, ref_indices_b)
  END SUBROUTINE test_intersection9

  SUBROUTINE test_intersection10
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 2, 5) /), &
         stripes_b(1) = (/ xt_stripe(1, 2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: dummy(1) = (/ -1_xi /)
    CALL test_intersection(stripes_a, stripes_b, dummy(1:0))
  END SUBROUTINE test_intersection10

  SUBROUTINE test_intersection11
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 5, 20) /), &
         stripes_b(1) = (/ xt_stripe(1, 7, 15) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(3) = (/ 15_xi, 50_xi, 85_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection11

  ! both ranges overlap in range but have no
  ! indices in common because of stride
  SUBROUTINE test_intersection12
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(34, 29, 12) /), &
         stripes_b(1) = (/ xt_stripe(36, 7, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: dummy(1) = (/ -1_xi /)

    CALL test_intersection(stripes_a, stripes_b, dummy(1:0))
  END SUBROUTINE test_intersection12

  ! same as test_intersection12 but with negative stride
  SUBROUTINE test_intersection13
    TYPE(xt_stripe), PARAMETER :: &
         stripes_a(1) = (/ xt_stripe(353, -29, 12) /), &
         stripes_b(1) = (/ xt_stripe(36, 7, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: dummy(1) = (/ -1_xi /)

    CALL test_intersection(stripes_a, stripes_b, dummy(1:0))
  END SUBROUTINE test_intersection13

  SUBROUTINE test_intersection14
    TYPE(xt_stripe), PARAMETER :: &
         stripes_a(1) = (/ xt_stripe(95, -29, 2) /), &
         stripes_b(1) = (/ xt_stripe(81, 14, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ 95_xi /)

    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection14

  SUBROUTINE test_get_pos(stripes, pos)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER, INTENT(in) :: pos(:)
    INTEGER(xt_int_kind), PARAMETER :: dummy = 1_xi
    INTEGER(xt_int_kind) :: ref_sel_idx(SIZE(pos)), sel_idx(SIZE(pos))
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -HUGE(dummy)
    INTEGER :: num_pos, ip, p, ref_undef_count, undef_count
    TYPE(xt_idxlist) :: idxlist
    idxlist = xt_idxstripes_new(stripes)
    num_pos = SIZE(pos)
    ref_undef_count = 0
    DO ip = 1, num_pos
      p = pos(ip)
      IF (xt_idxlist_get_index_at_position(idxlist, p, ref_sel_idx(ip))) THEN
        ref_sel_idx(ip) = undef_idx
        ref_undef_count = ref_undef_count + 1
      END IF
    END DO
    undef_count = xt_idxlist_get_indices_at_positions(idxlist, pos, sel_idx, &
         undef_idx)
    IF (undef_count /= ref_undef_count) &
         CALL test_abort("inequal undef count!", &
         __FILE__, &
         __LINE__)
    IF (ANY(sel_idx /= ref_sel_idx)) &
         CALL test_abort("incorrect index returned for position!", &
         __FILE__, &
         __LINE__)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_get_pos

  SUBROUTINE test_get_pos1
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 1, 5), xt_stripe(20, -1, 5) /)
    INTEGER, PARAMETER :: pos(13) = &
         (/   0,   2,   7,   9,  11, &
         &  100,  11, 200,   9, 300, &
         &   18, 400,   5 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos1

  SUBROUTINE test_get_pos2
    TYPE(xt_stripe), PARAMETER :: stripes(4) = (/ xt_stripe(0, 1, 3), &
         xt_stripe(10, 1, 2), xt_stripe(20, -1, 6), xt_stripe(30, -1, 7) /)
    INTEGER, PARAMETER :: pos(19) = &
         (/   -1,    0,    1,    2,    3,    4,   23,    5,    6,    7, &
         &     8,    9,   10,   11,   12,    0,    2,  100, 2000 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos2

  SUBROUTINE test_get_pos3
    TYPE(xt_stripe), PARAMETER :: stripes(4) = (/ xt_stripe(0, 1, 3), &
         xt_stripe(10, 1, 2), xt_stripe(20, -1, 6), xt_stripe(30, -1, 7) /)
    INTEGER, PARAMETER :: pos(13) = &
         (/    4,    7,    2,    5,    9,    0,   10,    6,   11,    8, &
         &    12,    1,    3 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos3

  SUBROUTINE test_get_pos4
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 1, 5), xt_stripe(20, -1, 5) /)
    INTEGER, PARAMETER :: pos(7) = &
         (/  -10,  200,  700,   90,   90,   18,  141 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos4

  SUBROUTINE test_stripe_overlap
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(1, 1, 5) /)
    INTEGER(xt_int_kind) :: i, j
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ ((i + j, i=0,4), j = 0, 1) /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE test_stripe_overlap

  SUBROUTINE test_stripe_bb(stripes, global_size, global_start_index, bounds_ref)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER(xt_int_kind), INTENT(in) :: global_size(:), global_start_index
    TYPE(xt_bounds), INTENT(in) :: bounds_ref(:)

    TYPE(xt_bounds) :: bounds(SIZE(global_size))
    TYPE(xt_idxlist) :: idxstripes

    IF (SIZE(global_size) /= SIZE(bounds_ref)) &
         CALL test_abort("size mismatch for bounding-box", &
         __FILE__, &
         __LINE__)
    idxstripes = xt_idxstripes_new(stripes, SIZE(stripes))

    bounds = xt_idxlist_get_bounding_box(idxstripes, global_size, &
         global_start_index)
    IF (ANY(bounds /= bounds_ref)) &
         CALL test_abort("boundary box doesn't match reference", &
         __FILE__, &
         __LINE__)
    CALL xt_idxlist_delete(idxstripes)
  END SUBROUTINE test_stripe_bb

  SUBROUTINE test_stripe_bb1
    TYPE(xt_stripe), PARAMETER :: stripes(1) = (/ xt_stripe(-1, -1, -1) /)
    INTEGER(xt_int_kind), PARAMETER :: global_size(3) = 4_xi, &
         global_start_index = 0
    TYPE(xt_bounds), PARAMETER :: bounds_ref(3) = xt_bounds(0, 0)
    CALL test_stripe_bb(stripes(1:0), global_size, global_start_index, bounds_ref)
  END SUBROUTINE test_stripe_bb1

  SUBROUTINE test_stripe_bb2
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(47, -12, 2), &
         xt_stripe(32, 12, 2), xt_stripe(36, 12, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: global_size(3) = (/ 5_xi, 4_xi, 3_xi /), &
         global_start_index = 1
    TYPE(xt_bounds), PARAMETER :: bounds_ref(3) = (/ xt_bounds(2, 2), &
         xt_bounds(2, 2), xt_bounds(1, 2) /)
    CALL test_stripe_bb(stripes, global_size, global_start_index, bounds_ref)
  END SUBROUTINE test_stripe_bb2

  SUBROUTINE do_tests(idxlist, ref_indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    TYPE(xt_stripe), ALLOCATABLE :: stripes(:)
    TYPE(xt_stripe), PARAMETER :: dummy(1) = (/ xt_stripe(0,0,0) /)
    INTEGER :: num_stripes
    TYPE(xt_idxlist) :: temp_idxlist, idxlist_copy

    CALL check_idxlist(idxlist, ref_indices)
    CALL xt_idxlist_get_index_stripes(idxlist, stripes)
    IF (ALLOCATED(stripes)) THEN
      num_stripes = SIZE(stripes)
      temp_idxlist = xt_idxvec_from_stripes_new(stripes, num_stripes)
    ELSE
      num_stripes = 0
      temp_idxlist = xt_idxvec_from_stripes_new(dummy, num_stripes)
    END IF
    CALL check_idxlist(temp_idxlist, ref_indices)

    CALL xt_idxlist_delete(temp_idxlist)

    IF (ALLOCATED(stripes)) DEALLOCATE(stripes)

    ! test packing and unpacking
    idxlist_copy = idxlist_pack_unpack_copy(idxlist)

    ! check copy
    CALL check_idxlist(idxlist_copy, ref_indices)

    CALL xt_idxlist_delete(idxlist_copy)

    ! test copying
    idxlist_copy = xt_idxlist_copy(idxlist)

    ! check copy
    CALL check_idxlist(idxlist_copy, ref_indices)

    ! clean up
    CALL xt_idxlist_delete(idxlist_copy)
  END SUBROUTINE do_tests

  SUBROUTINE check_pos_ext(stripes, search_stripes, ref_pos_ext, &
       single_match_only, ref_unmatched, test_desc)
    TYPE(xt_stripe), INTENT(in) :: stripes(:), search_stripes(:)
    TYPE(xt_pos_ext), intent(in) :: ref_pos_ext(:)
    LOGICAL, INTENT(in) :: single_match_only
    INTEGER, INTENT(in) :: ref_unmatched
    CHARACTER(len=*) :: test_desc

    INTEGER :: num_stripes, num_search_stripes, num_ref_pos_ext, num_ext, &
         unmatched
    TYPE(xt_idxlist) :: idxstripes
    TYPE(xt_pos_ext), ALLOCATABLE :: pos_ext(:)

    num_stripes = SIZE(stripes)
    num_search_stripes = SIZE(search_stripes)
    num_ref_pos_ext = SIZE(ref_pos_ext)

    idxstripes = xt_idxstripes_new(stripes)
    unmatched = xt_idxlist_get_pos_exts_of_index_stripes(idxstripes, &
         num_search_stripes, search_stripes, &
         num_ext, pos_ext, single_match_only)
    IF (unmatched /= ref_unmatched) &
         CALL test_abort("error in number of unmatched indices for " &
         // test_desc, &
         __FILE__, &
         __LINE__)
    IF (num_ext < 0 .OR. num_ext /= num_ref_pos_ext) &
         CALL test_abort("error finding " // test_desc, &
         __FILE__, &
         __LINE__)
    IF (ANY(pos_ext /= ref_pos_ext)) &
         CALL test_abort("incorrect position extent length found in "&
         // test_desc, &
         __FILE__, &
         __LINE__)
    DEALLOCATE(pos_ext)
    CALL xt_idxlist_delete(idxstripes)
  END SUBROUTINE check_pos_ext

  SUBROUTINE check_pos_ext1
    INTEGER, PARAMETER :: num_stripes = 1, num_ref_pos_ext = 1, &
         num_ref_unmatched = 0

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10) /), &
         search_stripes(1) = (/ xt_stripe(10_xi, -1_xi, 5) /)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(9, -5) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "simple inverted stripe")
  END SUBROUTINE check_pos_ext1

  SUBROUTINE check_pos_ext2
    INTEGER, PARAMETER :: num_stripes = 1, num_ref_pos_ext = 1, &
         num_ref_unmatched = 5

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10) /), &
         search_stripes(2) = xt_stripe(10_xi, -1_xi, 5)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(9, -5) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "simple inverted stripe")
  END SUBROUTINE check_pos_ext2

  SUBROUTINE check_pos_ext3
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10) /), &
         search_stripes(1) = xt_stripe(10_xi, 1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(9, 2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search inc stripe over inc gap")
  END SUBROUTINE check_pos_ext3

  SUBROUTINE check_pos_ext4
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(25_xi, -1_xi, 11), xt_stripe(10_xi, -1_xi, 10) /), &
         search_stripes(1) = xt_stripe(10_xi, 1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(11, -2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search inc stripe over dec gap")
  END SUBROUTINE check_pos_ext4

  SUBROUTINE check_pos_ext5
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(25_xi, -1_xi, 11), xt_stripe(10_xi, -1_xi, 10) /), &
         search_stripes(1) = xt_stripe(15_xi, -1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(10, 2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over dec gap")
  END SUBROUTINE check_pos_ext5

  SUBROUTINE check_pos_ext6
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10) /), &
         search_stripes(1) = xt_stripe(15_xi, -1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(10, -2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over inc gap")
  END SUBROUTINE check_pos_ext6

  SUBROUTINE check_pos_ext7
    INTEGER, PARAMETER :: num_stripes = 3, num_ref_pos_ext = 1, &
         num_ref_unmatched = 8

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10), &
         &    xt_stripe(29_xi, 1_xi, 10) /), &
         search_stripes(1) = xt_stripe(32_xi, -1_xi, 30)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(23, -22) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over 2 inc gap")
  END SUBROUTINE check_pos_ext7

  SUBROUTINE check_pos_ext8
    INTEGER, PARAMETER :: num_stripes = 5, num_ref_pos_ext = 5, &
         num_ref_unmatched = 0

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10), &
         &    xt_stripe(29_xi, 1_xi, 10), xt_stripe(14_xi, -1_xi, 4), &
         &    xt_stripe(28_xi, -1_xi, 4) /), &
         search_stripes(1) = xt_stripe(32_xi, -1_xi, 30)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(23, -4), xt_pos_ext(34, 4), xt_pos_ext(19, -10), &
         &    xt_pos_ext(30, 4), xt_pos_ext(9, -8) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over jumbled stripes")
  END SUBROUTINE check_pos_ext8

END PROGRAM test_idxstripes_f
