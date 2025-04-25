!>
!! @file test_redist_collection_f.f90
!! @brief Fortran test of redist_collection class
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
PROGRAM test_redist_collection
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, xi => xt_int_kind, &
       xt_xmap, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_collection_new, &
       xt_redist_delete, xt_redist_s_exchange1, xt_redist_s_exchange, &
       xt_slice_c_loc
  USE test_redist_common, ONLY: build_odd_selection_xmap
  USE iso_c_binding, ONLY: c_loc, c_ptr
  IMPLICIT NONE
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL simple_test
  CALL test_repeated_redist(-1)
  CALL test_repeated_redist(0)
  CALL test_displacement_variations

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE simple_test
    ! general test with one redist
    ! set up data
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_coll
    INTEGER, PARAMETER :: src_slice_len = 5, dst_slice_len = 3
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_slice_len) &
         = (/ 1.0d0, 3.0d0, 5.0d0 /)
    DOUBLE PRECISION, TARGET :: &
         src_data(src_slice_len) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /), &
         dst_data(dst_slice_len)

    dst_data = -1.0d0

    xmap = build_odd_selection_xmap(src_slice_len)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_collection
    redist_coll = xt_redist_collection_new((/ redist /), 1, -1, mpi_comm_world)

    CALL xt_redist_delete(redist)

    ! test exchange
    CALL xt_redist_s_exchange1(redist_coll, C_LOC(src_data), C_LOC(dst_data))


    IF (ANY(ref_dst_data /= dst_data)) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    ! clean up
    CALL xt_redist_delete(redist_coll)
  END SUBROUTINE simple_test

  SUBROUTINE test_repeated_redist_ds1(redist_coll)
    TYPE(xt_redist), INTENT(in) :: redist_coll
    INTEGER :: i, j
    DOUBLE PRECISION, SAVE, TARGET :: src_data(5, 3) = RESHAPE((/&
         (DBLE(i), i = 1, 15)/), (/ 5, 3 /))
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(3, 3) &
         = RESHAPE((/ ((DBLE(i + j), i = 1,5,2), j = 0,10,5) /), (/ 3, 3 /))
    DOUBLE PRECISION, TARGET :: dst_data(3, 3)
    TYPE(c_ptr) :: src_data_p(3), dst_data_p(3)
    dst_data = -1.0d0
    DO i = 1, 3
      CALL xt_slice_c_loc(src_data(:, i), src_data_p(i))
      CALL xt_slice_c_loc(dst_data(:, i), dst_data_p(i))
    END DO
    CALL xt_redist_s_exchange(redist_coll, 3, src_data_p, dst_data_p)

    IF (ANY(ref_dst_data /= dst_data)) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
  END SUBROUTINE test_repeated_redist_ds1

  SUBROUTINE test_repeated_redist_ds2(redist_coll)
    TYPE(xt_redist), INTENT(in) :: redist_coll
    INTEGER :: i, j
    DOUBLE PRECISION, SAVE, TARGET :: src_data(5, 3) = RESHAPE((/&
         (DBLE(i), i = 1, 15)/), (/ 5, 3 /))
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(3, 3) &
         = RESHAPE((/ ((DBLE(i + j), i = 1,5,2), j = 0,10,5) /), (/ 3, 3 /))
    DOUBLE PRECISION, TARGET :: dst_data(3, 3)
    TYPE(c_ptr) :: src_data_p(3), dst_data_p(3)
    dst_data = -1.0d0
    CALL xt_slice_c_loc(src_data(:, 2), src_data_p(1))
    CALL xt_slice_c_loc(src_data(:, 1), src_data_p(2))
    CALL xt_slice_c_loc(src_data(:, 3), src_data_p(3))
    CALL xt_slice_c_loc(dst_data(:, 2), dst_data_p(1))
    CALL xt_slice_c_loc(dst_data(:, 1), dst_data_p(2))
    CALL xt_slice_c_loc(dst_data(:, 3), dst_data_p(3))
    CALL xt_redist_s_exchange(redist_coll, 3, src_data_p, dst_data_p)

    IF (ANY(ref_dst_data /= dst_data)) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
  END SUBROUTINE test_repeated_redist_ds2

  SUBROUTINE test_repeated_redist(cache_size)
    INTEGER, INTENT(in) :: cache_size
    ! test with one redist used three times (with two different input data
    ! displacements -> test of cache) (with default cache size)
    ! set up data
    INTEGER, PARAMETER :: num_slice = 3
    INTEGER, PARAMETER :: src_slice_len = 5
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redists(num_slice), redist_coll

    xmap = build_odd_selection_xmap(src_slice_len)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_collection

    redists = redist
    redist_coll = xt_redist_collection_new(redists, 3, cache_size, &
         mpi_comm_world)

    CALL xt_redist_delete(redist)

    ! test exchange
    CALL test_repeated_redist_ds1(redist_coll)
    ! test exchange with changed displacements
    CALL test_repeated_redist_ds2(redist_coll)
    ! test exchange with original displacements
    CALL test_repeated_redist_ds1(redist_coll)
    ! clean up
    CALL xt_redist_delete(redist_coll)
  END SUBROUTINE test_repeated_redist

  ! test with one redist used three times (with different input
  ! data displacements until the cache is full)
  ! set up data
  SUBROUTINE test_displacement_variations
    INTEGER(xt_int_kind) :: i, j
    INTEGER, PARAMETER :: cache_size = 16
    INTEGER(xt_int_kind), PARAMETER :: num_slice = 3_xi, dst_step = 2_xi
    INTEGER, PARAMETER :: src_slice_len = 5
    INTEGER, PARAMETER :: dst_slice_len &
         = (src_slice_len + dst_step - 1)/dst_step
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redists(num_slice), redist_coll
    DOUBLE PRECISION, TARGET, SAVE :: src_data(src_slice_len, num_slice) &
         = RESHAPE((/ (DBLE(i), i = 1_xi, src_slice_len*num_slice) /), &
         (/ INT(src_slice_len), INT(num_slice) /))
    DOUBLE PRECISION, TARGET :: dst_data(dst_slice_len, num_slice)
    DOUBLE PRECISION, TARGET, ALLOCATABLE :: src_data_(:), dst_data_(:)
    TYPE(c_ptr) :: src_data_p(num_slice), dst_data_p(num_slice)
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_slice_len, num_slice) = &
         RESHAPE((/ ((DBLE(i + j * src_slice_len), &
         &            i = 1_xi, src_slice_len, dst_step), &
         &           j = 0_xi, num_slice - 1_xi) /), &
         &       (/ INT(dst_slice_len), INT(num_slice) /))
    INTEGER :: k

    xmap = build_odd_selection_xmap(src_slice_len)
    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_collection
    redists = redist

    redist_coll = xt_redist_collection_new(redists, INT(num_slice), &
         cache_size, mpi_comm_world)

    CALL xt_redist_delete(redist)

    dst_data = -1.0d0


    ALLOCATE(src_data_(src_slice_len + cache_size + 2), &
         dst_data_(dst_slice_len + cache_size + 2))

    DO i = 1, num_slice - 1
      CALL xt_slice_c_loc(src_data(:, i), src_data_p(i))
      CALL xt_slice_c_loc(dst_data(:, i), dst_data_p(i))
    END DO

    ! test exchange
    DO k = 1, cache_size + 2
      src_data_(k:k+src_slice_len-1) = src_data(:,num_slice)
      dst_data_(k:k+dst_slice_len-1) = dst_data(:,num_slice)

      CALL xt_slice_c_loc(src_data_(k:k+src_slice_len-1), src_data_p(3))
      CALL xt_slice_c_loc(dst_data_(k:k+dst_slice_len-1), dst_data_p(3))

      CALL xt_redist_s_exchange(redist_coll, INT(num_slice), src_data_p, &
           dst_data_p)

      IF (ANY(ref_dst_data(:, 1:num_slice-1) /= dst_data(:, 1:2)) &
           .OR. ANY(ref_dst_data(:,3) /= dst_data_(k:k+dst_slice_len-1))) &
           CALL test_abort("error in xt_redist_s_exchange", &
           __FILE__, &
           __LINE__)
    END DO

    ! clean up
    DEALLOCATE(src_data_, dst_data_)
    CALL xt_redist_delete(redist_coll)
  END SUBROUTINE test_displacement_variations

END PROGRAM test_redist_collection
