!
! @file test_redist_p2p_parallel_f.f90
!
! @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
!                                 Moritz Hanke <hanke@dkrz.de>
!                                 Thomas Jahns <jahns@dkrz.de>
!
! @author Jörg Behrens <behrens@dkrz.de>
!         Moritz Hanke <hanke@dkrz.de>
!         Thomas Jahns <jahns@dkrz.de>
!
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
PROGRAM test_redist_p2p_parallel
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE iso_c_binding, ONLY: c_ptr, c_loc
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, &
       xt_int_kind, xi => xt_int_kind, &
       xt_idxlist, xt_idxvec_new, xt_idxlist_delete, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_get_mpi_comm, &
       xt_redist_s_exchange1, xt_redist_delete, &
       xt_redist_p2p_blocks_off_new, xt_redist_p2p_blocks_new
  USE test_idxlist_utils, ONLY: test_err_count
  IMPLICIT NONE

  INTEGER :: rank, size, ierror

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL mpi_comm_rank(mpi_comm_world, rank, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("MPI error!", &
       __FILE__, &
       __LINE__)

  CALL mpi_comm_size(mpi_comm_world, size, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("MPI error!", &
       __FILE__, &
       __LINE__)

  CALL simple_test
  CALL nonuniform_test
  CALL block_redist_test

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi


CONTAINS
  SUBROUTINE simple_test
    INTEGER, PARAMETER :: data_size = 10
    INTEGER, PARAMETER :: src_num_indices = data_size, &
         dst_num_indices = data_size
    INTEGER(xt_int_kind) :: src_index_list(data_size), &
         dst_index_list(data_size)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist
    DOUBLE PRECISION, TARGET :: src_data(data_size), dst_data(data_size)
    INTEGER :: i

    ! source index list
    DO i = 1, src_num_indices
      src_index_list(i) = INT(rank * data_size + (i - 1), xi)
    END DO

    src_idxlist = xt_idxvec_new(src_index_list)
    ! destination index list
    DO i = 1, dst_num_indices
      dst_index_list(i) &
           = INT(MOD(rank * data_size + i + 1, size * data_size), xi)
    END DO
    dst_idxlist = xt_idxvec_new(dst_index_list)
    ! xmap
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)
    ! redist_p2p
    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    ! test communicator of redist
    IF (.NOT. test_communicator(xt_redist_get_mpi_comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", &
         __FILE__, &
         __LINE__)

    ! test exchange
    dst_data(:) = -1.0d0

    DO i = 1, src_num_indices
      src_data(i) = DBLE(rank * data_size + i - 1)
    END DO

    CALL xt_redist_s_exchange1(redist, C_LOC(src_data), C_LOC(dst_data))

    IF (ANY(DBLE(dst_index_list(:)) /= dst_data(:))) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
    ! clean up
    CALL xt_redist_delete(redist)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE simple_test

  ! test nonuniform numbers of send and receive partners
  SUBROUTINE nonuniform_test
    ! source index list
    INTEGER(xt_int_kind), ALLOCATABLE :: src_index_list(:), dst_index_list(:)
    DOUBLE PRECISION, ALLOCATABLE, TARGET :: src_data(:), dst_data(:)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(Xt_redist) :: redist
    INTEGER :: i, src_num_indices, dst_num_indices

    ALLOCATE(src_index_list(size), dst_index_list(size), &
         src_data(size), dst_data(size))
    src_num_indices = MERGE(size, 0, rank == 0)
    DO i = 1, src_num_indices
      src_index_list(i) = INT(i - 1, xi)
    END DO

    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices)

    ! destination index list
    dst_num_indices = size
    DO i = 1, dst_num_indices
      dst_index_list(i) = INT(i - 1, xi)
    END DO

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices)

    ! xmap
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)

    ! redist_p2p
    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    ! test communicator of redist
    IF (.NOT. test_communicator(xt_redist_get_mpi_comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", &
         __FILE__, &
         __LINE__)

    ! test exchange
    IF (rank == 0) THEN
      DO i = 1, size
        src_data(i) = DBLE(i - 1)
      END DO
    ELSE
      src_data(:) = -2.0d0
    END IF
    dst_data(:) = -1.0d0

    CALL xt_redist_s_exchange1(redist, C_LOC(src_data), C_LOC(dst_data))

    DO i = 1, size
      IF (dst_data(i) /= DBLE(i - 1)) &
           CALL test_abort("error in xt_redist_s_exchange", &
           __FILE__, &
           __LINE__)
    END DO

    ! clean up
    CALL xt_redist_delete(redist)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE nonuniform_test

  ! test redist with blocks
  SUBROUTINE block_redist_test
    ! gvol_size: volume of deep ocean
    INTEGER :: ngdom, gvol_size, i, nwin, ig0, ig, j, p, qa, qb, &
         a_vol_size, b_vol_size
    ! gdepth: ocean depth of an one dim. ocean
    INTEGER, ALLOCATABLE :: gdoma(:), gdomb(:), gsurfdata(:), &
         gdepth(:), ig2col_off(:), b_surfdata_ref(:), gvoldata(:), &
         src_block_offsets(:), src_block_sizes(:), dst_block_offsets(:), &
         dst_block_sizes(:), b_voldata_ref(:)
    INTEGER, TARGET, ALLOCATABLE :: a_surfdata(:), b_surfdata(:), &
         a_voldata(:), b_voldata(:)
    INTEGER(xi), ALLOCATABLE :: iveca(:), ivecb(:)
    INTEGER(xi) :: ia, ib
    TYPE(Xt_idxlist) :: idxlist_a, idxlist_b
    TYPE(xt_xmap) :: xmap
    TYPE(Xt_redist) :: redist, block_redist, block_redist2

    IF (2 * size > HUGE(1_xt_int_kind)) &
         CALL test_abort('too large number of tasks', &
         __FILE__, &
         __LINE__)
    ! the global index domain (1dim problem):
    ngdom = 2 * size
    ! start state (index distribution) of global domain
    ALLOCATE(gdoma(ngdom), gdomb(ngdom))
    ! end state ""
    ALLOCATE(gsurfdata(ngdom), gdepth(ngdom))
    ALLOCATE(ig2col_off(ngdom)) ! offset of surface DATA within vol
    gvol_size = 0
    DO i = 1, ngdom
      gdoma(i) = i - 1
      gdomb(i) = ngdom - i
      gsurfdata(i) = 99 + i
      gdepth(i) = i
      ig2col_off(i) = gvol_size
      gvol_size = gvol_size + gdepth(i)
    END DO

    nwin = ngdom / size ! my local window size of the global surface domain
    ! start of my window within global index domain (== global offset)
    ig0 = rank * nwin
    IF (nwin * size /= ngdom) &
         CALL test_abort("internal error", &
         __FILE__, &
         __LINE__)

    ! local index
    ALLOCATE(iveca(nwin), ivecb(nwin))
    DO i = 1, nwin
      ig = ig0 + i
      iveca(i) = INT(gdoma(ig), xi)
      ivecb(i) = INT(gdomb(ig), xi)
    END DO

    idxlist_a = xt_idxvec_new(iveca, nwin)
    idxlist_b = xt_idxvec_new(ivecb, nwin)

    xmap = xt_xmap_all2all_new(idxlist_a, idxlist_b, mpi_comm_world)

    ! simple redist
    redist = xt_redist_p2p_new(xmap, mpi_integer)

    ! test communicator of redist
    IF (.NOT. test_communicator(xt_redist_get_mpi_comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", &
         __FILE__, &
         __LINE__)

    ALLOCATE(a_surfdata(nwin), b_surfdata(nwin), b_surfdata_ref(nwin))
    DO i = 1, nwin
      a_surfdata(i) = gsurfdata(iveca(i) + 1)
      b_surfdata(i) = -1
      b_surfdata_ref(i) = gsurfdata(ivecb(i) + 1)
    END DO

    CALL xt_redist_s_exchange1(redist, C_LOC(a_surfdata), C_LOC(b_surfdata))
    IF (ANY(b_surfdata(:) /= b_surfdata_ref(:))) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
    CALL xt_redist_delete(redist)

    ! generate global volume data
    ALLOCATE(gvoldata(gvol_size))
    DO i = 1, ngdom
      DO j = 1, gdepth(i)
        p = ig2col_off(i) + j
        gvoldata(p) = (i - 1) * 100 + j - 1
      END DO
    END DO

    ! generate blocks
    ALLOCATE(src_block_offsets(nwin), src_block_sizes(nwin), &
         dst_block_offsets(nwin), dst_block_sizes(nwin))
    a_vol_size = 0 ! state a volume of my proc
    b_vol_size = 0 ! state b volume of my proc
    ! we only need local size but simply oversize here
    ALLOCATE(a_voldata(gvol_size), b_voldata(gvol_size), &
         b_voldata_ref(gvol_size))
    a_voldata(:) = -1
    b_voldata(:) = -1
    b_voldata_ref(:) = -1

    qa = 1
    DO i = 1, nwin
      ia = iveca(i)
      IF (i > 1) THEN
        src_block_offsets(i) = src_block_offsets(i - 1) + src_block_sizes(i - 1)
      ELSE
        src_block_offsets(i) = 0
      END IF
      src_block_sizes(i) = gdepth(INT(ia) + 1)
      DO j = 1, gdepth(INT(ia) + 1)
        p = ig2col_off(INT(ia) + 1) + j
        a_voldata(qa) = gvoldata(p)
        qa = qa + 1
      END DO
      a_vol_size = a_vol_size + src_block_sizes(i)
    END DO

    qb = 1
    DO i = 1, nwin
      ib = ivecb(i)
      IF (i > 1) THEN
        dst_block_offsets(i) = dst_block_offsets(i - 1) + dst_block_sizes(i - 1)
      ELSE
        dst_block_offsets(i) = 0
      END IF
      dst_block_sizes(i) = gdepth(INT(ib) + 1)
      DO j = 1, gdepth(INT(ib) + 1)
        p = ig2col_off(INT(ib) + 1) + j
        b_voldata_ref(qb) = gvoldata(p)
        qb = qb + 1
      END DO
      b_vol_size = b_vol_size + dst_block_sizes(i)
    END DO

    ! redist with blocks
    block_redist &
         = xt_redist_p2p_blocks_off_new(xmap, &
         src_block_offsets, src_block_sizes, nwin, &
         dst_block_offsets, dst_block_sizes, nwin, mpi_integer)
    ! test communicator of redist
    IF (.NOT. test_communicator(xt_redist_get_mpi_comm(block_redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", &
         __FILE__, &
         __LINE__)

    CALL xt_redist_s_exchange1(block_redist, &
         C_LOC(a_voldata), C_LOC(b_voldata))

    IF (ANY(b_voldata(:) /= b_voldata_ref(:))) &
         CALL test_abort("error in xt_redist_s_exchange (1) for volume data", &
         __FILE__, &
         __LINE__)

    ! redist with blocks but without explicit offsets:
    block_redist2 = xt_redist_p2p_blocks_new(xmap, &
         src_block_sizes, nwin, dst_block_sizes, nwin, mpi_integer)
    ! test communicator of redist

    IF (.NOT. test_communicator(xt_redist_get_mpi_comm(block_redist2), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", &
         __FILE__, &
         __LINE__)

    b_voldata(:) = -1

    CALL xt_redist_s_exchange1(block_redist2, &
         C_LOC(a_voldata), C_LOC(b_voldata))

    IF (ANY(b_voldata(:) /= b_voldata_ref(:))) &
         CALL test_abort("error in xt_redist_s_exchange (2) for volume data", &
         __FILE__, &
         __LINE__)

    ! cleanup
    CALL xt_redist_delete(block_redist2)
    CALL xt_redist_delete(block_redist)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(idxlist_a)
    CALL xt_idxlist_delete(idxlist_b)
  END SUBROUTINE block_redist_test

  FUNCTION test_communicator(comm1, comm2) RESULT(p)
    INTEGER, INTENT(in) :: comm1, comm2
    LOGICAL :: p
    INTEGER :: cmpres, ierror
    CALL mpi_comm_compare(comm1, comm2, cmpres, ierror)
    IF (ierror /= mpi_success) &
       CALL test_abort("MPI error!", &
       __FILE__, &
       __LINE__)
    p = ((cmpres == MPI_IDENT) .OR. (cmpres == MPI_CONGRUENT))
  END FUNCTION test_communicator

END PROGRAM test_redist_p2p_parallel
