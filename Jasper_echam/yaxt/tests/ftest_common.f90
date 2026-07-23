!>
!> @file ftest_common.f90
!>
!> @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
!>                                 Moritz Hanke <hanke@dkrz.de>
!>                                 Thomas Jahns <jahns@dkrz.de>
!>
!> @author Jörg Behrens <behrens@dkrz.de>
!>         Moritz Hanke <hanke@dkrz.de>
!>         Thomas Jahns <jahns@dkrz.de>
!>

!
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
MODULE ftest_common
  USE mpi
  USE xt_core, ONLY: i2, i4, i8
  IMPLICIT NONE
  PRIVATE
  INTEGER, PUBLIC, PARAMETER :: dp = SELECTED_REAL_KIND(12, 307)

  TYPE timer
    CHARACTER(len=20) :: label = 'undef'
    INTEGER  :: istate  = -1
    REAL(dp) :: t0      = 0.0_dp
    REAL(dp) :: dt_work = 0.0_dp
  END TYPE timer

  INTERFACE test_abort
    MODULE PROCEDURE test_abort_cmsl_f
    MODULE PROCEDURE test_abort_msl_f
  END INTERFACE test_abort

  INTERFACE icmp
    MODULE PROCEDURE icmp_2d
    MODULE PROCEDURE icmp_3d
  END INTERFACE icmp

  INTERFACE id_map
    MODULE PROCEDURE id_map_i2, id_map_i4, id_map_i8
  END INTERFACE id_map

  REAL(dp) :: sync_dt_sum = 0.0_dp
  LOGICAL, PARAMETER :: debug = .FALSE.
  LOGICAL :: verbose = .FALSE.

  PUBLIC :: init_mpi, finish_mpi
  PUBLIC :: timer, treset, tstart, tstop, treport, mysync
  PUBLIC :: id_map, icmp, factorize, regular_deco
  PUBLIC :: test_abort, set_verbose, get_verbose
CONTAINS

  SUBROUTINE init_mpi
    CHARACTER(len=*), PARAMETER :: context = 'init_mpi: '
    INTEGER :: ierror
    CALL mpi_init(ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_INIT failed', &
         __FILE__, __LINE__)
  END SUBROUTINE init_mpi

  SUBROUTINE finish_mpi
    CHARACTER(len=*), PARAMETER :: context = 'finish_mpi: '
    INTEGER :: ierror
    CALL MPI_FINALIZE(ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_FINALIZE failed', &
         __FILE__, &
         __LINE__)
  END SUBROUTINE finish_mpi

  SUBROUTINE set_verbose(verb)
    LOGICAL, INTENT(in) :: verb
    verbose = verb
  END SUBROUTINE set_verbose

  SUBROUTINE get_verbose(verb)
    LOGICAL, INTENT(out) :: verb
    verb = verbose
  END SUBROUTINE get_verbose

  PURE SUBROUTINE treset(t, label)
    TYPE(timer), INTENT(inout) :: t
    CHARACTER(len=*), INTENT(in) :: label
    t%label   = label
    t%istate  = 0
    t%t0      = 0.0_dp
    t%dt_work = 0.0_dp
  END SUBROUTINE treset

  SUBROUTINE tstart(t)
    TYPE(timer), INTENT(inout) :: t
    IF (debug) WRITE(0,*) 'tstart: ',t%label
    CALL mysync
    t%istate = 1
    t%t0 = work_time()
  END SUBROUTINE tstart

  SUBROUTINE tstop(t)
    TYPE(timer), INTENT(inout) :: t
    REAL(dp) :: t1
    IF (debug) WRITE(0,*) 'tstop: ',t%label
    t1 = work_time()
    t%dt_work = t%dt_work + (t1 - t%t0)
    t%istate = 0
    CALL mysync

  END SUBROUTINE tstop

  SUBROUTINE treport(t,extra_label,comm)
    TYPE(timer), INTENT(in) :: t
    CHARACTER(len=*), INTENT(in) :: extra_label
    INTEGER, INTENT(in) :: comm

    CHARACTER(len=*), PARAMETER :: context = 'treport: '
    REAL(dp) :: work_sum, work_max, work_avg, e
    REAL(dp) :: sbuf
    REAL(dp), ALLOCATABLE :: rbuf(:)
    INTEGER :: nprocs, rank, ierror

    sbuf = t%dt_work
    CALL mpi_comm_rank(comm, rank, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort(context//'MPI_COMM_RANK failed', &
         __FILE__, &
         __LINE__)
    CALL mpi_comm_size(comm, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort(context//'MPI_COMM_RANK failed', &
         __FILE__, &
         __LINE__)
    ALLOCATE(rbuf(0:nprocs-1))
    rbuf = -1.0_dp
    CALL mpi_gather(sbuf, 1, MPI_DOUBLE_PRECISION, &
         &  rbuf, 1, MPI_DOUBLE_PRECISION, &
         &  0, comm, ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_GATHER failed', &
         __FILE__, &
         __LINE__)

    IF (rank == 0) THEN
      IF (rbuf(0) /= sbuf) CALL test_abort(context//'internal error (1)', &
           __FILE__, &
           __LINE__)
      IF (ANY(rbuf < 0.0_dp)) CALL test_abort(context//'internal error (2)', &
           __FILE__, &
           __LINE__)
      work_sum = SUM(rbuf)
      work_max = MAXVAL(rbuf)
      work_avg = work_sum / REAL(nprocs, dp)
      e = work_avg / (work_max + 1.e-20_dp)

      IF (verbose) WRITE(0,'(A,I4,2X,A16,3F18.8)') &
           'nprocs, label, wmax, wavg, e =', &
           nprocs, extra_label//':'//t%label, &
           work_max, work_avg, e
    ENDIF

  END SUBROUTINE treport

  SUBROUTINE mysync
    CHARACTER(len=*), PARAMETER :: context = 'mysync: '
    INTEGER :: ierror
    REAL(dp) :: t0, dt

    t0 = mpi_wtime()

    CALL mpi_barrier(MPI_COMM_WORLD, IERROR)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_BARRIER failed', &
         __FILE__, &
         __LINE__)

    dt = (MPI_WTIME() - t0)
    sync_dt_sum = sync_dt_sum + dt

  END SUBROUTINE mysync

  REAL(dp) FUNCTION work_time()
    work_time = MPI_WTIME() - sync_dt_sum
    RETURN
  END FUNCTION work_time

  PURE SUBROUTINE id_map_i2(map)
    INTEGER(i2), INTENT(out) :: map(:,:)

    INTEGER :: i, j, m, n

    m = SIZE(map, 1)
    n = SIZE(map, 2)
    DO j = 1, n
      DO i = 1, m
        map(i,j) = INT((j - 1) * m + i, i2)
      ENDDO
    ENDDO

  END SUBROUTINE id_map_i2

  PURE SUBROUTINE id_map_i4(map)
    INTEGER(i4), INTENT(out) :: map(:,:)

    INTEGER :: i, j, m, n

    m = SIZE(map, 1)
    n = SIZE(map, 2)
    DO j = 1, n
      DO i = 1, m
        map(i,j) = INT((j - 1) * m + i, i4)
      ENDDO
    ENDDO

  END SUBROUTINE id_map_i4

  PURE SUBROUTINE id_map_i8(map)
    INTEGER(i8), INTENT(out) :: map(:,:)

    INTEGER :: i, j, m, n

    m = SIZE(map, 1)
    n = SIZE(map, 2)
    DO j = 1, n
      DO i = 1, m
        map(i,j) = INT((j - 1) * m + i, i8)
      ENDDO
    ENDDO

  END SUBROUTINE id_map_i8

  SUBROUTINE test_abort_msl_f(msg, source, line)
    CHARACTER(*), INTENT(in) :: msg
    CHARACTER(*), INTENT(in) :: source
    INTEGER, VALUE, INTENT(in) :: line
    CALL test_abort_cmsl_f(mpi_comm_world, msg, source, line)
  END SUBROUTINE test_abort_msl_f

  SUBROUTINE test_abort_cmsl_f(comm, msg, source, line)
    INTEGER, INTENT(in):: comm
    CHARACTER(*), INTENT(in) :: msg
    CHARACTER(*), INTENT(in) :: source
    INTEGER, VALUE, INTENT(in) :: line

    INTERFACE
      SUBROUTINE c_abort() BIND(c, name='abort')
      END SUBROUTINE c_abort
    END INTERFACE

    INTEGER :: ierror
    LOGICAL :: flag

    WRITE (0, '(3a,i0,2a)') 'Fatal error in ', source, ', line ', line, &
         ': ', msg
    FLUSH(0)
    CALL mpi_initialized(flag, ierror)
    IF (ierror == mpi_success .AND. flag) &
         CALL mpi_abort(comm, 1, ierror)
    CALL c_abort
  END SUBROUTINE test_abort_cmsl_f

  SUBROUTINE icmp_2d(label, f,g, rank)
    CHARACTER(len=*), PARAMETER :: context = 'ftest_common::icmp_2d: '
    CHARACTER(len=*), INTENT(in) :: label
    INTEGER, INTENT(in) :: f(:,:)
    INTEGER, INTENT(in) :: g(:,:)
    INTEGER, INTENT(in) :: rank

    INTEGER :: i, j, n1, n2

    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    IF (SIZE(g,1) /= n1 .OR. SIZE(g,2) /= n2) &
         CALL test_abort(context//'shape mismatch error', &
         __FILE__, &
         __LINE__)

    DO j = 1, n2
      DO i = 1, n1
        IF (f(i,j) /= g(i,j)) THEN
          WRITE(0,'(2a,4(a,i0))') context, label, ' test failed: i=', &
               i, ', j=', j, ', f(i,j)=', f(i,j), ', g(i,j)=', g(i,j)
          CALL test_abort(context//label//' test failed', &
               __FILE__, &
               __LINE__)
        ENDIF
      ENDDO
    ENDDO
    IF (verbose) WRITE(0,*) rank,':',context//label//' passed'
  END SUBROUTINE icmp_2d

  SUBROUTINE icmp_3d(label, f,g, rank)
    CHARACTER(len=*), PARAMETER :: context = 'ftest_common::icmp_3d: '
    CHARACTER(len=*), INTENT(in)  :: label
    INTEGER, INTENT(in)  :: f(:,:,:)
    INTEGER, INTENT(in)  :: g(:,:,:)
    INTEGER, INTENT(in) :: rank

    INTEGER :: i1, i2, i3, n1, n2, n3

    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    n3 = SIZE(f,3)
    IF (SIZE(g,1) /= n1 .OR. SIZE(g,2) /= n2 .OR. SIZE(g,3) /= n3) &
      CALL test_abort(context//label//'shape mismatch', __FILE__, __LINE__)

    DO i3 = 1, n3
      DO i2 = 1, n2
        DO i1 = 1, n1
          IF (f(i1,i2,i3) /= g(i1,i2,i3)) THEN
            WRITE(0,*) context,label,&
                 ' test failed: i1, i2, i3, f(i1,i2,i3), g(i1,i2,i3) =', &
                 i1, i2, i3, f(i1,i2,i3), g(i1,i2,i3)
            CALL test_abort(context//label//' test failed', &
                 __FILE__, &
                 __LINE__)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    IF (verbose) WRITE(0,*) rank,':',context//label//' passed'
  END SUBROUTINE icmp_3d

  SUBROUTINE factorize(c, a, b)
    INTEGER, INTENT(in) :: c
    INTEGER, INTENT(out) :: a, b ! c = a*b

    INTEGER :: x0, i

    IF (c<1) CALL test_abort('factorize: invalid process space', &
         __FILE__, &
         __LINE__)
    IF (c <= 3 .OR. c == 5 .OR. c == 7) THEN
      a = c
      b = 1
      RETURN
    ENDIF

    ! simple approach, we try to be near c = (2*x) * x
    x0 = INT(SQRT(0.5 * REAL(c)) + 0.5)
    a = 2*x0
    f_loop: DO i = a, 1, -1
      IF (MOD(c,i) == 0) THEN
        a = i
        b = c/i
        EXIT f_loop
      ENDIF
    ENDDO f_loop

  END SUBROUTINE factorize

  SUBROUTINE regular_deco(g_cn, c0, cn)
    INTEGER, INTENT(in) :: g_cn
    INTEGER, INTENT(out) :: c0(0:), cn(0:)

    ! convention: process space coords start at 0, grid point coords start at 1

    integer :: tn
    INTEGER :: d, m
    INTEGER :: it

    tn = SIZE(c0)
    IF (tn<0) CALL test_abort('(tn<0)', __FILE__, __LINE__)
    IF (tn>g_cn) CALL test_abort('regular_deco: too many task for such a core&
         & region', &
         __FILE__, &
         __LINE__)

    d = g_cn/tn
    m = MOD(g_cn, tn)

    DO it = 0, m-1
      cn(it) = d + 1
    ENDDO
    DO it = m, tn-1
      cn(it) = d
    ENDDO

    c0(0)=0
    DO it = 1, tn-1
      c0(it) = c0(it-1) + cn(it-1)
    ENDDO
    IF (c0(tn-1)+cn(tn-1) /= g_cn) &
         CALL test_abort('regular_deco: internal error 1', &
         __FILE__, &
         __LINE__)
  END SUBROUTINE regular_deco

END MODULE ftest_common
