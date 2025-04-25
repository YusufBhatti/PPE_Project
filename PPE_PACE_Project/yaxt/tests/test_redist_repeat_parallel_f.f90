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
       xt_idxlist, xt_idxlist_delete, xt_stripe, &
       xt_idxsection_new, xt_idxlist_collection_new, xt_idxstripes_new, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_repeat_new, &
       xt_redist_delete, xt_redist_s_exchange, &
       xt_idxlist_get_indices, xt_int_mpidt
  USE iso_c_binding, ONLY: c_loc, c_ptr, c_int
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
  END IF

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE build_idxlists(indices_a, indices_b)
    ! redist test with four different redists
    TYPE(xt_idxlist), INTENT(out) :: indices_a, indices_b

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

    stripe = xt_stripe(start = INT(rank * 2 * size**2, xi), stride = 1_xi, &
                       nstrides = INT(2*size**2, c_int))
    indices_b = xt_idxstripes_new(stripe)
  END SUBROUTINE build_idxlists

  SUBROUTINE test_4redist
    TYPE(xt_idxlist) :: indices_a, indices_b
    INTEGER(xt_int_kind) :: index_vector_a(2*size**2), &
                            index_vector_b(2*size**2)
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist_repeat, redist_repeat_2, redist_p2p
    INTEGER(xt_int_kind), TARGET :: results_1(2*size**2,4), &
                                    results_2(2*size**2,9)
    INTEGER(xt_int_kind), TARGET :: input_data(2*size**2,9)
    INTEGER(xt_int_kind) :: ref_results_1(2*size**2,4), &
                            ref_results_2(2*size**2,9)
    TYPE(c_ptr) :: results(1), input(1)
    INTEGER(mpi_address_kind) :: extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int) :: displacements(4), displacements_2(4)
    INTEGER :: i, ierror
    INTEGER(xt_int_kind) :: j

    CALL build_idxlists(indices_a, indices_b)

    CALL xt_idxlist_get_indices(indices_a, index_vector_a)
    CALL xt_idxlist_get_indices(indices_b, index_vector_b)

    xmap = xt_xmap_all2all_new(indices_a, indices_b, mpi_comm_world)

    CALL xt_idxlist_delete(indices_a)
    CALL xt_idxlist_delete(indices_b)

    redist_p2p = xt_redist_p2p_new(xmap, dt_xt_int)
    CALL xt_xmap_delete(xmap)

    CALL mpi_get_address(input_data(1,1), base_address, ierror)
    CALL mpi_get_address(input_data(1,2), temp_address, ierror)
    extent = temp_address - base_address

    displacements = (/0,1,2,3/)
    displacements_2 = (/1,2,4,8/)

    redist_repeat = xt_redist_repeat_new(redist_p2p, extent, extent, &
         4, displacements)
    redist_repeat_2 = xt_redist_repeat_new(redist_p2p, extent, extent, &
         4, displacements_2)

    CALL xt_redist_delete(redist_p2p)

    DO j = 1, 9
      DO i = 1, 2*size**2
        input_data(i, j) = index_vector_a(i) + (j - 1) * INT(4*size**2, xi)
      END DO
    END DO
    results_1 = -1
    results_2 = -1

    input(1) = C_LOC(input_data)
    results(1) = C_LOC(results_1)

    CALL xt_redist_s_exchange(redist_repeat, 1, input, results)
    results(1) = C_LOC(results_2)
    CALL xt_redist_s_exchange(redist_repeat_2, 1, input, results)

    DO j = 1, 4
      DO i = 1, 2*size**2
        ref_results_1(i, j) = index_vector_b(i) + (j - 1) * INT(4*size**2, xi)
      END DO
    END DO
    DO j = 1, 9
      DO i = 1, 2*size**2
        ref_results_2(i, j) = index_vector_b(i) + (j - 1) * INT(4*size**2, xi)
      END DO
    END DO

    ref_results_2(:,1:4:3) = -1
    ref_results_2(:,6:8:1) = -1

    ! check results
    IF (ANY(results_1 /= ref_results_1)) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
    IF (ANY(results_2 /= ref_results_2)) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    ! clean up

    CALL xt_redist_delete(redist_repeat)
    CALL xt_redist_delete(redist_repeat_2)
  END SUBROUTINE test_4redist

END PROGRAM test_redist_collection_static_parallel
