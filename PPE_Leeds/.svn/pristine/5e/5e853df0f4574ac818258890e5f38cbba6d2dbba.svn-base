! recreates whole array from variously displaced chunks
PROGRAM perf_sparse_array_gather
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_abort, &
       xt_int_kind, xt_stripe, xt_abort, xt_idxlist_collection_new, &
       xt_idxlist, xt_idxfsection_new, xt_idxstripes_new, xt_idxlist_delete, &
       xt_redist, xt_redist_p2p_ext_new, xt_offset_ext, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist_delete, xt_redist_s_exchange1, &
       xt_idxstripes_from_idxlist_new
  USE iso_c_binding, ONLY: c_int, c_loc, c_ptr
  USE mpi
  IMPLICIT NONE

  INTEGER :: gshape(3), gstart(3), nparts(2), num_parts, sshape(3), gsize, &
       section_shape(3), section_start(3), section_size, ofs, ierror
  DOUBLE PRECISION, TARGET, ALLOCATABLE :: gathered_data(:, :, :)
  DOUBLE PRECISION, TARGET, ALLOCATABLE :: scattered_data(:, :, :)
  TYPE(c_ptr) :: scattered_data_p, gathered_data_p
  TYPE(xt_offset_ext) :: gather_ext(1)
  TYPE(xt_offset_ext), ALLOCATABLE :: scatter_exts(:,:)
  TYPE(xt_idxlist) :: gather_list
  TYPE(xt_idxlist), ALLOCATABLE :: scatter_lists(:,:)
  TYPE(xt_idxlist) :: scatter_list, scatter_list_stripes
  TYPE(xt_xmap) :: gather_xmap
  TYPE(xt_redist) :: gather_redist
  REAL :: r
  INTEGER :: i, j, k, t
  CHARACTER(len=132) :: msg
  LOGICAL :: check_correctness_of_exchange
  INTERFACE
    SUBROUTINE POSIX_EXIT(code) BIND(c, name="exit")
      IMPORT :: c_int
      INTEGER(c_int), VALUE, INTENT(in) :: code
    END SUBROUTINE POSIX_EXIT
  END INTERFACE

  CALL mpi_init(ierror)
  IF (ierror /= mpi_success) THEN
    WRITE (0, *) "MPI initialization failed"
    CALL POSIX_EXIT(1)
  END IF
  CALL xt_initialize(mpi_comm_world)

  check_correctness_of_exchange = .FALSE.
  CALL parse_options

  gstart(:) = 1
  gshape(1) = 192
  gshape(2) = 96
  gshape(3) = 49
  gsize = PRODUCT(gshape)

  nparts(1) = 8
  nparts(2) = 8
  num_parts = PRODUCT(nparts)
  sshape(1) = (gsize + num_parts - 1) / num_parts * 20
  sshape(2) = nparts(1)
  sshape(3) = nparts(2)

  ALLOCATE(gathered_data(gshape(1), gshape(2), gshape(3)), &
       scattered_data(sshape(1), sshape(2), sshape(3)), &
       scatter_exts(nparts(1), nparts(2)), scatter_lists(nparts(1), nparts(2)))

  gather_ext(1) = xt_offset_ext(0, gsize, 1)

  gathered_data = -HUGE(1.0d0)

  ! perform 2d deco of cube
  ! divide gshape(1:2) by nparts(1:2)
  section_start(3) = 1
  section_shape(3) = gshape(3)
  scattered_data_p = C_LOC(scattered_data)
  gathered_data_p = C_LOC(gathered_data)
  DO t = 1, 100
    gather_list = xt_idxstripes_new(xt_stripe(0, 1, gsize))
    CALL RANDOM_NUMBER(r)
    DO j = 1, nparts(2)
      section_start(2) = (gshape(2) * (j - 1))/nparts(2) + 1
      section_shape(2) = (gshape(2) * j)/nparts(2) - section_start(2) + 1
      DO i = 1, nparts(1)
        section_start(1) = (gshape(1) * (i - 1))/nparts(1) + 1
        section_shape(1) = (gshape(1) * i)/nparts(1) - section_start(1) + 1
        section_size = PRODUCT(section_shape(:))
        ofs = INT(r * REAL(sshape(1) - 1 - section_size))
        scatter_exts(i, j) &
             = xt_offset_ext(ofs + ((i-1) + (j-1) * sshape(2)) * sshape(1), &
             section_size, 1)
        scatter_lists(i, j) = xt_idxfsection_new(0_xt_int_kind, &
             INT(gshape, xt_int_kind), &
             section_shape, INT(section_start, xt_int_kind))
        scattered_data(1:ofs, i, j) = -9.0d99
        scattered_data(ofs+1:ofs + section_size, i, j) &
             = DBLE(i + (j - 1) * nparts(1))
        scattered_data(ofs + section_size + 1:, i, j) = -9.0d99
      END DO
    END DO
    scatter_list = xt_idxlist_collection_new(scatter_lists)
    CALL xt_idxlist_delete(scatter_lists)
    scatter_list_stripes = xt_idxstripes_from_idxlist_new(scatter_list)
    CALL xt_idxlist_delete(scatter_list)
    gather_xmap = xt_xmap_all2all_new(scatter_list_stripes, gather_list, &
         mpi_comm_self)
    CALL xt_idxlist_delete(gather_list)
    gather_redist = xt_redist_p2p_ext_new(gather_xmap, &
         RESHAPE(scatter_exts, (/ num_parts /)), gather_ext, &
         mpi_double_precision)
    CALL xt_xmap_delete(gather_xmap)
    CALL xt_redist_s_exchange1(gather_redist, scattered_data_p, gathered_data_p)
    IF (check_correctness_of_exchange) THEN
      DO k = 1, gshape(3)
        DO j = 1, gshape(2)
          DO i = 1, gshape(1)
            IF (gathered_data(i, j, k) < 1.0d0 &
                 .OR. gathered_data(i, j, k) > DBLE(num_parts)) THEN
              WRITE (0, "(a,3(i0,a),g28.16)") "gathered_data(", &
                   i, ", ", j, ", ", k, ") = ", gathered_data(i, j, k)
              WRITE (msg, "(a,3(', ',i0))") "error in data", i, j, k
              CALL xt_abort(TRIM(msg), &
                   __FILE__, &
                   __LINE__)
            END IF
          END DO
        END DO
      END DO
    END IF
    CALL xt_redist_delete(gather_redist)
  END DO
  CALL xt_finalize
  CALL mpi_finalize(ierror)
  IF (ierror /= mpi_success) THEN
    WRITE (0, *) "MPI finalization failed"
    CALL posix_exit(1)
  END IF

CONTAINS
  SUBROUTINE parse_options
    INTEGER :: i, num_cmd_args, arg_len
    INTEGER, PARAMETER :: max_opt_arg_len = 80
    CHARACTER(max_opt_arg_len) :: optarg
    CHARACTER(len=18), PARAMETER :: &
         check_correctness_argstr = '-check-correctness'
    num_cmd_args = COMMAND_ARGUMENT_COUNT()
    i = 1
    DO WHILE (i < num_cmd_args)
      CALL GET_COMMAND_ARGUMENT(i, optarg, arg_len)
      IF (optarg(1:arg_len) == check_correctness_argstr &
           .AND. arg_len == LEN(check_correctness_argstr)) THEN
        check_correctness_of_exchange = .TRUE.
      ELSE
        WRITE (0, *) 'unexpected command-line argument parsing error: ', &
             TRIM(optarg)
        FLUSH(0)
        CALL xt_abort(mpi_comm_world, 'unexpected command-line argument "' &
             // optarg(1:arg_len) // '"', &
             __FILE__, &
             __LINE__)
      END IF
      i = i + 1
    END DO

  END SUBROUTINE parse_options

END PROGRAM perf_sparse_array_gather
