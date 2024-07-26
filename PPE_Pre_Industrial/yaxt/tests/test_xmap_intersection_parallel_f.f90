!>
!! @file test_xmap_intersection_parallel_f.f90
!!
!! @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!
!!
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
PROGRAM test_xmap_intersection_parallel
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE mpi
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, &
       xt_idxlist, xt_idxvec_new, xt_idxlist_delete, xt_xmap, &
       xt_idxempty_new, xi => xt_int_kind, &
       xt_xmap_intersection_new, xt_xmap_intersection_ext_new, xt_com_list, &
       xt_xmap_delete, xt_xmap_iter, &
       xt_xmap_get_in_iterator, xt_xmap_get_out_iterator, &
       xt_xmap_get_num_destinations, xt_xmap_get_num_sources, &
       xt_xmap_iterator_get_rank, xt_xmap_iterator_get_num_transfer_pos, &
       xt_xmap_iterator_get_transfer_pos, xt_xmap_iterator_next, &
       xt_xmap_iterator_delete
#if defined __PGI && ( __PGIC__ < 12 || (__PGIC__ ==  11 && __PGIC_MINOR__ <= 7))
  ! PGI Fortran 12.7 and older has a bug that prevents proper passing of
  ! generic interfaces through multiple modules, direct USE instead
  USE xt_xmap_abstract, ONLY: xt_is_null
#else
  USE yaxt, ONLY: xt_is_null
#endif
  IMPLICIT NONE

  TYPE test_message
    INTEGER :: rank       ! rank of communication partner
    INTEGER, POINTER :: pos(:)   ! positions to be sent/received
  END TYPE test_message
  INTERFACE
    SUBROUTINE posix_exit(code) BIND(c, name='exit')
      USE iso_c_binding, ONLY: c_int
      INTEGER(c_int), VALUE, INTENT(in) :: code
    END SUBROUTINE posix_exit
  END INTERFACE

  INTEGER, PARAMETER :: xmi_type_base = 0, xmi_type_ext = 1
  INTEGER :: xmi_type

  INTEGER :: ierror
  INTEGER :: my_rank, comm_size

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  xmi_type = xmi_type_base
  CALL parse_options

  CALL mpi_comm_rank(mpi_comm_world, my_rank, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("MPI error!", &
       __FILE__, &
       __LINE__)

  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("MPI error!", &
       __FILE__, &
       __LINE__)

  IF (comm_size /= 3) THEN
    CALL xt_finalize
    CALL finish_mpi
    CALL posix_exit(77)
  END IF

  ! parse_options(&argc, &argv);
  CALL simple_rr_test
  CALL elimination_test
  CALL one_to_one_comm_test
  CALL full_comm_matrix_test
  CALL dedup_test

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi

CONTAINS
  ! simple test (round robin)
  SUBROUTINE simple_rr_test
    ! setup
    INTEGER(xi) :: src_index(1), dst_index(1)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    INTEGER, PARAMETER :: num_src_intersections = 1, &
         num_dst_intersections = 1, num_sends = 1, num_recvs = 1
    INTEGER, SAVE, TARGET :: send_pos(num_sends) = (/ 0 /), &
         recv_pos(num_recvs) = (/ 0 /)
    TYPE(xt_com_list) :: src_com(num_src_intersections), &
         dst_com(num_dst_intersections)
    TYPE(xt_xmap) :: xmap
    TYPE(test_message) :: send_messages(1), recv_messages(1)

    src_index(1) = INT(my_rank, xi)
    dst_index(1) = INT(MOD(my_rank + 1, comm_size), xi)
    src_idxlist = xt_idxvec_new(src_index)
    dst_idxlist = xt_idxvec_new(dst_index)
    src_com(1) = xt_com_list(src_idxlist, MOD(my_rank+1, comm_size))
    dst_com(1) = xt_com_list(dst_idxlist, MOD(my_rank+comm_size-1, comm_size))

    xmap = xmi_new(src_com(1:num_src_intersections), &
         dst_com(1:num_dst_intersections), &
         src_idxlist, dst_idxlist, mpi_comm_world)

    ! test
    send_messages(1)%rank = MOD(my_rank+1, comm_size)
    send_messages(1)%pos => send_pos
    recv_messages(1)%rank = MOD(my_rank+comm_size-1, comm_size)
    recv_messages(1)%pos => recv_pos

    CALL test_xmap(xmap, send_messages, recv_messages)

    ! cleanup
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(dst_idxlist)
    CALL xt_idxlist_delete(src_idxlist)
  END SUBROUTINE simple_rr_test

  ! rank 0 receives the same point from rank 1 and 2
  SUBROUTINE elimination_test
    INTEGER(xi), PARAMETER :: src_index(1) = (/ 0_xi /), &
         dst_index(1) = (/ 0_xi /)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    INTEGER :: num_src_intersections,  num_dst_intersections, num_sends, &
         num_recvs
    INTEGER, TARGET :: send_pos(1), recv_pos(1)
    TYPE(xt_com_list) :: src_com(1), dst_com(2)
    TYPE(xt_xmap) :: xmap
    TYPE(test_message) :: send_messages(1), recv_messages(1)
    ! setup
    IF (my_rank == 0) THEN
      src_idxlist = xt_idxempty_new()
      dst_idxlist = xt_idxvec_new(dst_index)
    ELSE
      src_idxlist = xt_idxvec_new(src_index)
      dst_idxlist = xt_idxempty_new()
    END IF
    num_src_intersections = MERGE(1, 0, my_rank /= 0)
    src_com = xt_com_list(src_idxlist, 0)
    num_dst_intersections = MERGE(0, 2, my_rank /= 0)
    dst_com(1) = xt_com_list(dst_idxlist, 1)
    dst_com(2) = xt_com_list(dst_idxlist, 2)

    xmap = xmi_new(src_com(1:num_src_intersections), &
         dst_com(1:num_dst_intersections), &
         src_idxlist, dst_idxlist, mpi_comm_world)

    ! test
    send_pos(1) = 0
    num_sends = MERGE(1, 0, my_rank == 1)
    send_messages(1)%rank = 0
    send_messages(1)%pos => send_pos
    recv_pos(1) = 0;
    num_recvs = MERGE(1, 0, my_rank == 0)
    recv_messages(1)%rank = 1
    recv_messages(1)%pos => recv_pos

    CALL test_xmap(xmap, send_messages(1:num_sends), recv_messages(1:num_recvs))

    ! cleanup

    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(dst_idxlist)
    CALL xt_idxlist_delete(src_idxlist)
  END SUBROUTINE elimination_test

  ! all ranks can receive data from one of the others
  SUBROUTINE one_to_one_comm_test
    ! rank               |  0  |  1  |  2  |
    ! source indices     | 1,2 | 2,0 | 0,1 |
    ! destination indice |  0  |  1  |  2  |

    INTEGER(xi) :: src_indices(2), dst_index(1)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist, src_intersection_idxlist(2)
    INTEGER, PARAMETER :: num_src_intersections(3) = (/ 2, 1, 0 /)
    INTEGER :: num_sends, num_recvs, s_s, s_e, i
    INTEGER, TARGET :: send_pos(2), recv_pos(1)
    TYPE(xt_com_list) :: src_com(2), dst_com(1)
    TYPE(xt_xmap) :: xmap
    TYPE(test_message) :: send_messages(2), recv_messages(1)
    ! setup
    dst_index(1) = INT(my_rank, xi)
    DO i = 1, 2
      src_indices(i) = INT(MOD(my_rank+i, comm_size), xi)
      src_intersection_idxlist(i) = xt_idxvec_new(src_indices(i:i), 1)
    END DO
    src_idxlist = xt_idxvec_new(src_indices, 2)
    dst_idxlist = xt_idxvec_new(dst_index, 1)
    src_com(1) = xt_com_list(src_intersection_idxlist(1), 1)
    src_com(2) = xt_com_list(src_intersection_idxlist(2), &
         MERGE(2, 0, my_rank == 0))
    dst_com = xt_com_list(dst_idxlist, MERGE(1, 0, my_rank == 0))
    s_s = MERGE(my_rank + 1, 1, my_rank /= 2)
    s_e = s_s + num_src_intersections(my_rank + 1) - 1
    xmap = xmi_new(src_com(s_s:s_e), dst_com(:), src_idxlist, dst_idxlist, &
         mpi_comm_world)

    ! test
    recv_pos(1) = 0
    num_recvs = 1
    SELECT CASE (my_rank)
    CASE (0)
      send_pos(1) = 0
      send_pos(2) = 1
      num_sends = 2
      send_messages(1)%rank = 1
      send_messages(1)%pos => send_pos(1:1)
      send_messages(2)%rank = 2
      send_messages(2)%pos => send_pos(2:2)
      recv_messages(1)%rank = 1
    CASE (1)
      send_pos = 1
      num_sends = 1
      send_messages(1)%rank = 0
      send_messages(1)%pos => send_pos(1:1)
      recv_messages(1)%rank = 0
    CASE default
      num_sends = 0
      recv_messages(1)%rank = 0
    END SELECT
    recv_messages(1)%pos => recv_pos(1:1)
    CALL test_xmap(xmap, send_messages(1:num_sends), recv_messages(1:num_recvs))

    ! cleanup
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_intersection_idxlist(2))
    CALL xt_idxlist_delete(src_intersection_idxlist(1))
    CALL xt_idxlist_delete(dst_idxlist)
    CALL xt_idxlist_delete(src_idxlist)
  END SUBROUTINE one_to_one_comm_test

  ! all ranks receive data from each of the others
  SUBROUTINE full_comm_matrix_test
    !rank               |        0        |        1        |        2
    !source indices     |    0,1,2,3,4    |    3,4,5,6,7    |    6,7,8,0,1
    !destination indices|0,1,2,3,4,5,6,7,8|0,1,2,3,4,5,6,7,8|0,1,2,3,4,5,6,7,8


    INTEGER(xi), PARAMETER :: src_indices(5,0:2) &
         = RESHAPE((/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi, &
         &            3_xi, 4_xi ,5_xi, 6_xi, 7_xi, &
         &            6_xi, 7_xi, 8_xi, 0_xi, 1_xi /), (/ 5, 3 /)), &
         dst_indices(9) &
         = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 8_xi /)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_com_list) :: src_com(0:2), dst_com(0:2)
    TYPE(xt_xmap) :: xmap
    INTEGER, SAVE, TARGET :: send_pos(5, 0:2) &
         = RESHAPE((/ 0,1,2,3,4, 2,3,4,-1,-1, 2,-1,-1,-1,-1 /), (/ 5, 3 /)), &
         num_send_pos(0:2) = (/ 5, 3, 1 /), &
         recv_pos(5, 0:2) &
         = RESHAPE((/ 0,1,2,3,4, 5,6,7,-1,-1, 8,-1,-1,-1,-1 /), (/ 5, 3 /)), &
         num_recv_pos(0:2) = (/ 5, 3, 1 /)
    TYPE(test_message) :: send_messages(0:2), recv_messages(0:2)
    INTEGER :: i

    ! setup
    src_idxlist = xt_idxvec_new(src_indices(:, my_rank))
    dst_idxlist = xt_idxvec_new(dst_indices, 9)
    DO i = 0, 2
      src_com(i) = xt_com_list(src_idxlist, i)
      dst_com(i) = xt_com_list(xt_idxvec_new(src_indices(:, i)), i)
    END DO
    xmap = xmi_new(src_com, dst_com, src_idxlist, dst_idxlist, &
         mpi_comm_world)

    ! test
    DO i = 0, 2
      send_messages(i)%rank = i
      send_messages(i)%pos => send_pos(1:num_send_pos(my_rank), my_rank)
      recv_messages(i)%rank = i
      recv_messages(i)%pos => recv_pos(1:num_recv_pos(i), i)
    END DO
    CALL test_xmap(xmap, send_messages, recv_messages)

    ! cleanup
    CALL xt_xmap_delete(xmap)
    DO i = 2, 0, -1
      CALL xt_idxlist_delete(dst_com(i)%list)
    END DO
    CALL xt_idxlist_delete(dst_idxlist)
    CALL xt_idxlist_delete(src_idxlist)
  END SUBROUTINE full_comm_matrix_test

  ! one rank receives data from the other two, that have duplicated indices
  ! (this provokes a bug found by Joerg Behrens)
  SUBROUTINE dedup_test
    ! rank                |  0  |  1  |   2   |
    ! source indices      | 0,2 | 1,2 |       |
    ! destination indices |     |     | 0,1,2 |

    INTEGER(xt_int_kind), PARAMETER :: src_indices(2, 0:1) &
         = RESHAPE((/ 0_xi,2_xi, 1_xi,2_xi /), (/ 2, 2 /))
    TYPE(xt_com_list) ::  src_com(1), dst_com(2)
    INTEGER :: num_src_intersections, num_dst_intersections
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    INTEGER(xt_int_kind), PARAMETER :: dst_indices(3) = (/ 0_xi, 1_xi, 2_xi /)
    TYPE(xt_xmap) :: xmap
    INTEGER :: i, num_recv_messages, num_send_messages
    INTEGER, PARAMETER :: num_recv_pos(2) = (/ 2, 1 /), &
         num_send_pos(2) = (/ 2, 1 /)
    INTEGER, SAVE, TARGET :: &
         recv_pos(2, 2) = RESHAPE((/ 0, 2, 1, -1 /), (/ 2, 2 /)), &
         send_pos(2, 2) = RESHAPE((/ 0, 1, 0, -1 /), (/ 2, 2 /))
    TYPE(test_message) :: recv_messages(2), send_messages(1)
    ! setup
    IF (my_rank == 2) THEN
      num_src_intersections = 0
      num_dst_intersections = 2
      DO i = 0, 1
        dst_com(i+1)%list = xt_idxvec_new(src_indices(:, i))
        dst_com(i+1)%rank = i
      END DO
      src_idxlist = xt_idxempty_new()
      dst_idxlist = xt_idxvec_new(dst_indices(:))
    ELSE
      num_src_intersections = 1
      src_com(1)%list = xt_idxvec_new(src_indices(:, my_rank))
      src_com(1)%rank = 2
      num_dst_intersections = 0;
      src_idxlist = xt_idxvec_new(src_indices(:, my_rank))
      dst_idxlist = xt_idxempty_new()
    END IF
    xmap = xmi_new(src_com(1:num_src_intersections), &
         dst_com(1:num_dst_intersections), &
         src_idxlist, dst_idxlist, mpi_comm_world)

    ! test
    IF (my_rank == 2) THEN
      num_recv_messages = 2
      num_send_messages = 0
      DO i = 1, 2
        recv_messages(i)%rank = i - 1
        recv_messages(i)%pos => recv_pos(1:num_recv_pos(i), i)
      END DO
    ELSE
      num_recv_messages = 0
      num_send_messages = 1
      send_messages(1)%rank = 2
      send_messages(1)%pos => send_pos(1:num_send_pos(my_rank + 1), my_rank + 1)
    END IF
    CALL test_xmap(xmap, send_messages(1:num_send_messages), &
         recv_messages(1:num_recv_messages))

    ! cleanup
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(dst_idxlist)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_com(1:num_dst_intersections)%list)
    CALL xt_idxlist_delete(src_com(1:num_src_intersections)%list)
  END SUBROUTINE dedup_test

  SUBROUTINE test_xmap_iter(iter, msgs)
    TYPE(xt_xmap_iter), INTENT(inout) :: iter
    TYPE(test_message), INTENT(in) :: msgs(:)

    INTEGER :: num_msgs, num_pos, i, j
    INTEGER, POINTER :: pos(:)
    LOGICAL :: iter_is_null

    num_msgs = SIZE(msgs)
    iter_is_null = xt_is_null(iter)
    IF (num_msgs == 0) THEN
      IF (.NOT. iter_is_null) &
           CALL test_abort('ERROR: xt_xmap_get_*_iterator (non-null when ' &
           // 'iter should be null)', &
           __FILE__, &
           __LINE__)
    ELSE IF (iter_is_null) THEN
      CALL test_abort('ERROR: xt_xmap_get_*_iterator ' &
           // '(iter should not be NULL)', &
           __FILE__, &
           __LINE__)
    ELSE
      i = 1
      DO WHILE(.TRUE.)
        IF (xt_xmap_iterator_get_rank(iter) /= msgs(i)%rank) &
             CALL test_abort('ERROR: xt_xmap_iterator_get_rank', &
             __FILE__, &
             __LINE__)
        num_pos = SIZE(msgs(i)%pos)
        IF (xt_xmap_iterator_get_num_transfer_pos(iter) /= num_pos) THEN
          CALL test_abort("ERROR: xt_xmap_iterator_get_num_transfer_pos", &
               __FILE__, &
               __LINE__)
        END IF
        pos => xt_xmap_iterator_get_transfer_pos(iter)
        DO j = 1, num_pos
          IF (pos(j) /= msgs(i)%pos(j)) &
               CALL test_abort('ERROR: xt_xmap_iterator_get_transfer_pos', &
               __FILE__, &
               __LINE__)
        END DO
        IF (.NOT. xt_xmap_iterator_next(iter)) EXIT
        i = i + 1
      END DO
      IF (i /= num_msgs) &
           CALL test_abort('ERROR: xt_xmap_iterator_next &
           &(wrong number of messages)', &
           __FILE__, &
           __LINE__)
    END IF
  END SUBROUTINE test_xmap_iter

  SUBROUTINE test_xmap(xmap, send_messages, recv_messages)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(test_message), INTENT(in) :: send_messages(:), recv_messages(:)

    INTEGER :: num_sends, num_recvs
    TYPE(xt_xmap_iter) :: send_iter, recv_iter

    num_sends = SIZE(send_messages)
    num_recvs = SIZE(recv_messages)
    IF (xt_xmap_get_num_destinations(xmap) /= num_sends) &
         CALL test_abort('ERROR: xt_xmap_get_num_destinations', &
         __FILE__, &
         __LINE__)
    IF (xt_xmap_get_num_sources(xmap) /= num_recvs) &
         CALL test_abort('ERROR: xt_xmap_get_num_sources', &
         __FILE__, &
         __LINE__)
    send_iter = xt_xmap_get_out_iterator(xmap)
    recv_iter = xt_xmap_get_in_iterator(xmap)

    CALL test_xmap_iter(send_iter, send_messages)
    CALL test_xmap_iter(recv_iter, recv_messages)

    IF (.NOT. xt_is_null(recv_iter)) CALL xt_xmap_iterator_delete(recv_iter)
    IF (.NOT. xt_is_null(send_iter)) CALL xt_xmap_iterator_delete(send_iter)
  END SUBROUTINE test_xmap

  SUBROUTINE parse_options
    INTEGER :: i, num_cmd_args, arg_len
    INTEGER, PARAMETER :: max_opt_arg_len = 80
    CHARACTER(max_opt_arg_len) :: optarg
    num_cmd_args = COMMAND_ARGUMENT_COUNT()
    i = 1
    DO WHILE (i < num_cmd_args)
      CALL GET_COMMAND_ARGUMENT(i, optarg, arg_len)
      IF (optarg(1:2) == '-m' .AND. i < num_cmd_args .AND. arg_len == 2) THEN
        CALL GET_COMMAND_ARGUMENT(i + 1, optarg, arg_len)
        IF (arg_len > max_opt_arg_len) &
             CALL test_abort('incorrect argument to command-line option -m', &
             __FILE__, &
             __LINE__)
        IF (optarg(1:arg_len) == "xt_xmap_intersection_new") THEN
          xmi_type = xmi_type_base
        ELSE IF (optarg(1:arg_len) == "xt_xmap_intersection_ext_new") THEN
          xmi_type = xmi_type_ext
        ELSE
          WRITE (0, *) 'arg to -m: ', optarg(1:arg_len)
          CALL test_abort('incorrect argument to command-line option -m', &
               __FILE__, &
               __LINE__)
        END IF
        i = i + 2
      ELSE
        WRITE (0, *) 'unexpected command-line argument parsing error: ', &
             TRIM(optarg)
        FLUSH(0)
        CALL test_abort('unexpected command-line argument -m', &
             __FILE__, &
             __LINE__)
      END IF
    END DO
  END SUBROUTINE parse_options

  FUNCTION xmi_new(src_com, dst_com, src_idxlist, dst_idxlist, comm) &
       RESULT(xmap)
    TYPE(xt_com_list), INTENT(in) :: src_com(:), dst_com(:)
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap
    SELECT CASE(xmi_type)
    CASE(xmi_type_base)
      xmap = xt_xmap_intersection_new(src_com, dst_com, &
           src_idxlist, dst_idxlist, comm)
    CASE(xmi_type_ext)
      xmap = xt_xmap_intersection_ext_new(src_com, dst_com, &
           src_idxlist, dst_idxlist, comm)
    END SELECT
  END FUNCTION xmi_new

END PROGRAM test_xmap_intersection_parallel
