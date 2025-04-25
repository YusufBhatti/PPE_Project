!>
!! @file test_redist_collection_static_parallel_f.f90
!! @brief Fortran test of redist_collection_static class
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
PROGRAM test_redist_collection_static_parallel
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, xi => xt_int_kind, &
       xt_idxlist, xt_idxlist_delete, xt_stripe, xt_idxvec_new, &
       xt_idxsection_new, xt_idxlist_collection_new, xt_idxstripes_new, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_collection_static_new, &
       xt_redist_delete, xt_redist_s_exchange, xt_redist_s_exchange1, &
       xt_idxlist_get_indices, xt_int_mpidt
  USE iso_c_binding, ONLY: c_loc, c_ptr
  IMPLICIT NONE
  INTEGER(xt_int_kind) :: dummy
  INTEGER :: rank, size, ierror, dt_xt_int
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL mpi_comm_rank(mpi_comm_world, rank, ierror)
  IF (ierror /= MPI_SUCCESS) &
       CALL test_abort('mpi_comm_rank failed', &
       __FILE__, &
       __LINE__)
  CALL mpi_comm_size(mpi_comm_world, size, ierror)
  IF (ierror /= MPI_SUCCESS) &
       CALL test_abort('mpi_comm_size failed', &
       __FILE__, &
       __LINE__)
  IF (xt_int_mpidt == mpi_datatype_null) THEN
    CALL mpi_type_create_f90_integer(RANGE(dummy), dt_xt_int, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort('mpi_type_create_f90_integer failed', &
         __FILE__, &
         __LINE__)
  ELSE
    dt_xt_int = xt_int_mpidt
  END IF

  IF (size > 1) THEN
    CALL test_4redist
    CALL test_rr_exchange
  END IF

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE build_idxlists(indices_a, indices_b, indices_all)
    ! redist test with four different redists
    TYPE(xt_idxlist), INTENT(out) :: indices_a, indices_b, indices_all

    TYPE(xt_idxlist) :: indices_a_(2)
    INTEGER :: i
    INTEGER(xt_int_kind), PARAMETER :: start = 0
    INTEGER(xt_int_kind) :: global_size(2), local_start(2, 2)
    INTEGER :: local_size(2)

    TYPE(xt_stripe) :: stripe

    global_size(1) = INT(2 * size, xi)
    global_size(2) = INT(size**2, xi)
    local_size = size
    local_start = RESHAPE((/ 0_xi, INT(rank*size, xi), INT(size, xi), &
         INT(size*size-(rank+1)*size, xi) /), (/ 2, 2 /))

    DO i = 1, 2
      indices_a_(i) = xt_idxsection_new(start, global_size, local_size, &
           local_start(:, i))
    END DO
    indices_a = xt_idxlist_collection_new(indices_a_)

    CALL xt_idxlist_delete(indices_a_(1))
    CALL xt_idxlist_delete(indices_a_(2))

    stripe = xt_stripe(INT(rank * 2 * size**2, xi), 1_xi, 2*size**2)
    indices_b = xt_idxstripes_new(stripe)

    stripe = xt_stripe(0_xi, 1_xi, 2*size**3)
    indices_all = xt_idxstripes_new(stripe)
  END SUBROUTINE build_idxlists

  SUBROUTINE test_4redist
    INTEGER, PARAMETER :: num_tx = 4
    TYPE(xt_idxlist) :: indices_a, indices_b, indices_all
    INTEGER(xt_int_kind), TARGET :: index_vector_a(2*size**2), &
         index_vector_b(2*size**2), index_vector_all(2*size**3)
    TYPE(xt_xmap) :: xmaps(num_tx)
    TYPE(xt_redist) :: redists(num_tx), redist
    INTEGER(xt_int_kind), TARGET :: results_1(2*size**2), &
         results_2(2*size**2), results_3(2*size**3), results_4(2*size**3)
    TYPE(c_ptr) :: results(1), input(1)
    INTEGER(mpi_address_kind) :: src_displacements(num_tx), &
         dst_displacements(num_tx)
    INTEGER :: i, ierror

    CALL build_idxlists(indices_a, indices_b, indices_all)

    CALL xt_idxlist_get_indices(indices_a, index_vector_a)
    CALL xt_idxlist_get_indices(indices_b, index_vector_b)
    CALL xt_idxlist_get_indices(indices_all, index_vector_all)

    xmaps(1) = xt_xmap_all2all_new(indices_a, indices_b, mpi_comm_world)
    xmaps(2) = xt_xmap_all2all_new(indices_b, indices_a, mpi_comm_world)
    xmaps(3) = xt_xmap_all2all_new(indices_a, indices_all, mpi_comm_world)
    xmaps(4) = xt_xmap_all2all_new(indices_b, indices_all, mpi_comm_world)

    CALL xt_idxlist_delete(indices_a)
    CALL xt_idxlist_delete(indices_b)
    CALL xt_idxlist_delete(indices_all)

    DO i = 1, num_tx
      redists(i) = xt_redist_p2p_new(xmaps(i), dt_xt_int)
      CALL xt_xmap_delete(xmaps(i))
    END DO

    CALL mpi_get_address(index_vector_a, src_displacements(1), ierror)
    CALL mpi_get_address(index_vector_b, src_displacements(2), ierror)
    CALL mpi_get_address(index_vector_a, src_displacements(3), ierror)
    CALL mpi_get_address(index_vector_b, src_displacements(4), ierror)

    src_displacements = src_displacements - src_displacements(1)

    CALL mpi_get_address(results_1, dst_displacements(1), ierror)
    CALL mpi_get_address(results_2, dst_displacements(2), ierror)
    CALL mpi_get_address(results_3, dst_displacements(3), ierror)
    CALL mpi_get_address(results_4, dst_displacements(4), ierror)

    dst_displacements = dst_displacements - dst_displacements(1)

    redist = xt_redist_collection_static_new(redists, num_tx, &
         src_displacements, dst_displacements, mpi_comm_world)

    ! test communicator of redist
    ! if (!test_communicator(xt_redist_get_MPI_Comm(redist), MPI_COMM_WORLD))
    !   PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    CALL xt_redist_delete(redists)

    input(1) = C_LOC(index_vector_a)
    results(1) = C_LOC(results_1)

    CALL xt_redist_s_exchange(redist, 1, input, results)

    ! check results
    IF (ANY(results_1(:) /= index_vector_b)) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    IF (ANY(results_2(:) /= index_vector_a)) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    IF (ANY(results_3(:) /= index_vector_all)) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    IF (ANY(results_3(:) /= index_vector_all)) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    ! clean up

    CALL xt_redist_delete(redist)
  END SUBROUTINE test_4redist

  ! redist test with two redists that do a round robin exchange in
  ! different directions
  SUBROUTINE test_rr_exchange
    TYPE(xt_idxlist) :: src_indices, dst_indices(2)
    INTEGER, PARAMETER :: num_local_indices = 5
    INTEGER(xt_int_kind), TARGET :: src_indices_(num_local_indices)
    INTEGER(xt_int_kind) :: i, temp, dst_indices_(num_local_indices, 2)
    TYPE(xt_xmap) :: xmaps(2)
    TYPE(xt_redist) :: redists(2), redist
    INTEGER(xt_int_kind), TARGET :: results_1(num_local_indices)
    INTEGER(xt_int_kind), VOLATILE, TARGET ::results_2(num_local_indices)
    INTEGER(mpi_address_kind) :: src_displacements(2), dst_displacements(2), &
         addr_temp
    INTEGER :: ierror

    DO i = 1_xi, INT(num_local_indices, xi)
      src_indices_(i) = INT(rank, xi) * INT(num_local_indices, xi) + (i - 1_xi)
      dst_indices_(i, 1) = MOD(src_indices_(i) + 1_xi, &
           &                   INT(size, xi) * INT(num_local_indices, xi))
      temp = src_indices_(i) - 1_xi
      dst_indices_(i, 2) = MERGE(INT(size, xi) * INT(num_local_indices, xi) &
           &                     - 1_xi, &
           &                     temp, temp < 0_xi)
    END DO

    src_indices = xt_idxvec_new(src_indices_, num_local_indices)
    dst_indices(1) = xt_idxvec_new(dst_indices_(:, 1))
    dst_indices(2) = xt_idxvec_new(dst_indices_(:, 2))

    xmaps(1) = xt_xmap_all2all_new(src_indices, dst_indices(1), mpi_comm_world)
    xmaps(2) = xt_xmap_all2all_new(src_indices, dst_indices(2), mpi_comm_world)

    CALL xt_idxlist_delete(src_indices)
    CALL xt_idxlist_delete(dst_indices)

    redists(1) = xt_redist_p2p_new(xmaps(1), dt_xt_int)
    redists(2) = xt_redist_p2p_new(xmaps(2), dt_xt_int)

    CALL xt_xmap_delete(xmaps)

    src_displacements = 0_mpi_address_kind
    dst_displacements(1) = 0_mpi_address_kind
    CALL mpi_get_address(results_2, dst_displacements(2), ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error in mpi_get_address", &
         __FILE__, &
         __LINE__)
    CALL mpi_get_address(results_1, addr_temp, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error in mpi_get_address", &
         __FILE__, &
         __LINE__)
    dst_displacements(2) = dst_displacements(2) - addr_temp

    redist = xt_redist_collection_static_new(redists, 2, src_displacements, &
         dst_displacements, mpi_comm_world)

    ! test communicator of redist
    ! IF (!test_communicator(xt_redist_get_MPI_Comm(redist), MPI_COMM_WORLD))
    !     PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    CALL xt_redist_delete(redists)

    results_1 = -1
    results_2 = -1

    CALL xt_redist_s_exchange1(redist, C_LOC(src_indices_), C_LOC(results_1))

    ! check results
    IF (ANY(results_1(:) /= dst_indices_(:, 1))) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
    IF (ANY(results_2(:) /= dst_indices_(:, 2))) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    ! clean up
    CALL xt_redist_delete(redist)
  END SUBROUTINE test_rr_exchange

END PROGRAM test_redist_collection_static_parallel
