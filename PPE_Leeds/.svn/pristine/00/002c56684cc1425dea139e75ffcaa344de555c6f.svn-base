!>
!> @file test_perf.f90
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

PROGRAM test_perf
  USE mpi
  USE xt_ut, ONLY: ut_init, ut_mode_dt_p2p, comm_forward, ut_abort, &
       & ut_transpose, ut_init_decomposition, &
       & ut_init_oneway_transposition_template, ut_destroy_decomposition, &
       & ut_destroy_transposition_template, ut_destroy_transposition, &
       & ut_init_transposition, xt_int_kind
  USE ftest_common, ONLY: finish_mpi, init_mpi, treset, tstart, tstop, &
       treport, timer, id_map, factorize, regular_deco, set_verbose
  USE yaxt, ONLY: xt_finalize

  IMPLICIT NONE
  ! global extents including halos:

  INTEGER, PARAMETER :: nlev = 30
  INTEGER, PARAMETER :: undef_int = HUGE(undef_int)/2 - 1
  INTEGER(xt_int_kind), PARAMETER :: undef_index = -1
  INTEGER, PARAMETER :: nhalo = 1 ! 1dim. halo border size

  INTEGER, PARAMETER :: grid_kind_test = 1
  INTEGER, PARAMETER :: grid_kind_toy  = 2
  INTEGER, PARAMETER :: grid_kind_tp10 = 3
  INTEGER, PARAMETER :: grid_kind_tp04 = 4
  INTEGER, PARAMETER :: grid_kind_tp6M = 5
  INTEGER :: grid_kind = grid_kind_test

  CHARACTER(len=10) :: grid_label
  INTEGER :: g_ie, g_je ! global domain extents
  INTEGER :: ie, je, ke ! local extents, including halos
  INTEGER :: p_ioff, p_joff ! offsets within global domain
  INTEGER :: nprocx, nprocy ! process space extents
  INTEGER :: nprocs ! == nprocx*nprocy
  ! process rank, process coords within (0:, 0:) process space
  INTEGER :: mype, mypx, mypy
  LOGICAL :: lroot ! true only for proc 0

  INTEGER, ALLOCATABLE :: g_id(:,:) ! global id
  ! global "tripolar-like" toy bounds exchange
  INTEGER, ALLOCATABLE :: g_tpex(:, :)
  INTEGER :: template_tpex_2d, trans_tpex_2d
  INTEGER :: template_tpex_3d, trans_tpex_3d

  INTEGER, ALLOCATABLE :: loc_id_2d(:,:), loc_tpex_2d(:,:)
  INTEGER, ALLOCATABLE :: loc_id_3d(:,:,:), loc_tpex_3d(:,:,:)
  INTEGER, ALLOCATABLE :: fval_2d(:,:), gval_2d(:,:)
  INTEGER, ALLOCATABLE :: fval_3d(:,:,:), gval_3d(:,:,:)
  INTEGER, ALLOCATABLE :: id_pos(:,:), pos3d_surf(:,:)
  INTEGER :: surf_trans_tpex_2d
  LOGICAL :: verbose

  TYPE(timer) :: t_all, t_surf_trans, t_exch_surf
  TYPE(timer) :: t_template_2d, t_trans_2d, t_exch_2d
  TYPE(timer) :: t_template_3d, t_trans_3d, t_exch_3d

  CALL treset(t_all, 'all')
  CALL treset(t_surf_trans, 'surf_trans')
  CALL treset(t_exch_surf, 'exch_surf')
  CALL treset(t_template_2d, 'template_2d')
  CALL treset(t_trans_2d, 'trans_2d')
  CALL treset(t_exch_2d, 'exch_2d')
  CALL treset(t_template_3d, 'template_3d')
  CALL treset(t_trans_3d, 'trans_3d')
  CALL treset(t_exch_3d, 'exch_3d')

  CALL init_mpi

  CALL tstart(t_all)

  ! mpi & decomposition & allocate mem:
  CALL init_all

  ! full global index space:
  CALL id_map(g_id)

  ! local window of global index space:
  CALL get_window(g_id, loc_id_2d)


  ! define bounds exchange for full global index space
  CALL def_exchange()

  ! local window of global bounds exchange:
  CALL get_window(g_tpex, loc_tpex_2d)

  ! template: loc_id_2d -> loc_tpex_2d
  CALL tstart(t_template_2d)
  CALL gen_template_2d(loc_id_2d, loc_tpex_2d, template_tpex_2d)
  CALL tstop(t_template_2d)

  ! transposition: loc_id_2d:data -> loc_tpex_2d:data
  CALL tstart(t_trans_2d)
  CALL gen_trans(template_tpex_2d, MPI_INTEGER, MPI_INTEGER, trans_tpex_2d)
  CALL tstop(t_trans_2d)


  ! test 2d-to-2d transposition:
  fval_2d = loc_id_2d
  CALL tstart(t_exch_2d)
  CALL exchange2d_2d(trans_tpex_2d, fval_2d, gval_2d)
  CALL tstop(t_exch_2d)
  CALL icmp_2d('2d to 2d check', gval_2d, loc_tpex_2d)

  ! define positions of surface elements within (i,k,j) array
  CALL id_map(id_pos)
  CALL id_map(pos3d_surf)
  CALL gen_pos3d_surf(pos3d_surf)

  ! generate surface transposition:
  CALL tstart(t_surf_trans)
  CALL gen_off_trans(template_tpex_2d, MPI_INTEGER, id_pos(:,:)-1, &
       MPI_INTEGER, pos3d_surf(:,:)-1, surf_trans_tpex_2d)
  CALL tstop(t_surf_trans)


  ! 2d to surface boundsexchange:
  ALLOCATE(gval_3d(ie,nlev,je))
  gval_3d = -1
  CALL tstart(t_exch_surf)
  CALL exchange2d_3d(surf_trans_tpex_2d, fval_2d, gval_3d)
  CALL tstop(t_exch_surf)

  CALL icmp_2d('surface check', gval_3d(:,1,:), loc_tpex_2d)

  ! check sub surface:
  CALL icmp_2d('sub surface check', gval_3d(:,2,:), loc_tpex_2d*0-1)

  ! cleanup 2d:
  CALL ut_destroy_transposition_template(template_tpex_2d)
  CALL ut_destroy_transposition(trans_tpex_2d)
  CALL ut_destroy_transposition(surf_trans_tpex_2d)

  ! inflate 2d -> 3d
  CALL inflate_idx(3, loc_id_2d, loc_id_3d)
  CALL inflate_idx(3, loc_tpex_2d, loc_tpex_3d)

  ! template: loc_id_3d -> loc_tpex_3d
  CALL tstart(t_template_3d)
  CALL gen_template_3d(loc_id_3d, loc_tpex_3d, template_tpex_3d)
  CALL tstop(t_template_3d)

  ! transposition: loc_id_3d:data -> loc_tpex_3d:data
  CALL tstart(t_trans_3d)
  CALL gen_trans(template_tpex_3d, MPI_INTEGER, MPI_INTEGER, trans_tpex_3d)
  CALL tstop(t_trans_3d)

  ! test 3d-to-3d transposition:
  DEALLOCATE(gval_3d)
  ALLOCATE(fval_3d(ie,je,nlev), gval_3d(ie,je,nlev))
  fval_3d = loc_id_3d
  CALL tstart(t_exch_3d)
  CALL exchange3d_3d(trans_tpex_3d, fval_3d, gval_3d)
  CALL tstop(t_exch_3d)
  CALL icmp_3d('3d to 3d check', gval_3d, loc_tpex_3d)

  ! cleanup 3d:
  CALL ut_destroy_transposition_template(template_tpex_3d)
  CALL ut_destroy_transposition(trans_tpex_3d)

  CALL tstop(t_all)

  IF (verbose) WRITE(0,*) 'timer report for nprocs=',nprocs

  CALL treport(t_all, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_surf_trans, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_exch_surf, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_template_2d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_trans_2d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_exch_2d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_template_3d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_trans_3d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_exch_3d, TRIM(grid_label), mpi_comm_world)

  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE msg(s)
    CHARACTER(len=*), INTENT(in) :: s
    IF (verbose) WRITE(0,*) s
  END SUBROUTINE msg

  SUBROUTINE inflate_idx(inflate_pos, idx_2d, idx_3d)
    CHARACTER(len=*), PARAMETER :: context = 'test_perf::inflate_idx: '
    INTEGER, INTENT(in) :: inflate_pos
    INTEGER, INTENT(in) :: idx_2d(:,:)
    INTEGER, ALLOCATABLE, INTENT(out) :: idx_3d(:,:,:)

    INTEGER :: i, j, k

    IF (ALLOCATED(idx_3d)) DEALLOCATE(idx_3d)

    IF (inflate_pos == 3) THEN
      ALLOCATE(idx_3d(ie, je, ke))
      DO k=1,ke
        DO j=1,je
          DO i=1,ie
            idx_3d(i,j,k) = idx_2d(i,j) + (k-1) * g_ie * g_je
          ENDDO
        ENDDO
      ENDDO
    ELSE
      CALL ut_abort(context//' unsupported inflate position', &
           __FILE__, &
           __LINE__)
    ENDIF

  END SUBROUTINE inflate_idx

  SUBROUTINE gen_pos3d_surf(pos)
    INTEGER, INTENT(inout) :: pos(:,:)
    ! positions for zero based arrays (ECHAM grid point dim order)
    ! old pos = i + j*ie
    ! new pos = i + k*ie + j*ie*nlev
    INTEGER :: ii,jj, i,j,k, p,q

    k = 0 ! surface
    DO jj=1,je
      DO ii=1,ie
        p = pos(ii,jj) - 1 ! shift to 0-based index
        j = p/ie
        i = MOD(p,ie)
        q = i + k*ie + j*ie*nlev
        pos(ii,jj) = q + 1 ! shift to 1-based index
      ENDDO
    ENDDO

  END SUBROUTINE gen_pos3d_surf

  SUBROUTINE icmp_2d(label, f,g)
    CHARACTER(len=*), PARAMETER :: context = 'test_perf::icmp_2d: '
    CHARACTER(len=*), INTENT(in)  :: label
    INTEGER, INTENT(in)  :: f(:,:)
    INTEGER, INTENT(in)  :: g(:,:)

    INTEGER :: i, j, n1, n2

    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    IF (SIZE(g,1) /= n1 .OR. SIZE(g,2) /= n2) &
         CALL ut_abort(context//'internal error', &
         __FILE__, &
         __LINE__)

    DO j = 1, n2
      DO i = 1, n1
        IF (f(i,j) /= g(i,j)) THEN
          WRITE(0,*) context,label,' test failed: i, j, f(i,j), g(i,j) =', &
               i, j, f(i,j), g(i,j)
          CALL ut_abort(context//label//' test failed', &
               __FILE__, &
               __LINE__)
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE icmp_2d

  SUBROUTINE icmp_3d(label, f,g)
    CHARACTER(len=*), PARAMETER :: context = 'test_perf::icmp_3d: '
    CHARACTER(len=*), INTENT(in)  :: label
    INTEGER, INTENT(in)  :: f(:,:,:)
    INTEGER, INTENT(in)  :: g(:,:,:)

    INTEGER :: i, j, k, n1, n2, n3

    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    n3 = SIZE(f,3)
    IF (SIZE(g,1) /= n1 .OR. SIZE(g,2) /= n2 .OR. SIZE(g,3) /= n3) &
         CALL ut_abort(context//'internal error', &
         __FILE__, &
         __LINE__)

    DO k = 1, n3
      DO j = 1, n2
        DO i = 1, n1
          IF (f(i,j,k) /= g(i,j,k)) THEN
             WRITE(0,*) context,label, &
                  ' test failed: i, j, f(i,j,k), g(i,j,k) =', &
                  i, j, k, f(i,j,k), g(i,j,k)
            CALL ut_abort(context//label//' test failed', &
                 __FILE__, &
                 __LINE__)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE icmp_3d

  SUBROUTINE init_all
    CHARACTER(len=*), PARAMETER :: context = 'init_all: '
    INTEGER :: ierror
    CHARACTER(len=20) :: grid_str

    CALL get_environment_variable('YAXT_TEST_PERF_GRID', grid_str)

    verbose = .TRUE.

    SELECT CASE (TRIM(ADJUSTL(grid_str)))
    CASE('TOY')
      grid_kind = grid_kind_toy
      grid_label = 'TOY'
      g_ie = 66
      g_je = 36
    CASE('TP10')
      grid_kind = grid_kind_tp10
      grid_label = 'TP10'
      g_ie = 362
      g_je = 192
    CASE('TP04')
      grid_kind = grid_kind_tp04
      grid_label = 'TP04'
      g_ie = 802
      g_je = 404
    CASE('TP6M')
      grid_kind = grid_kind_tp6m
      grid_label = 'TP6M'
      g_ie = 3602
      g_je = 2394
    CASE default
      grid_kind = grid_kind_test
      grid_label = 'TEST'
      g_ie = 32
      g_je = 12
      verbose = .FALSE.
    END SELECT

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) CALL ut_abort(context//'MPI_COMM_SIZE failed', &
         __FILE__, &
         __LINE__)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierror)
    IF (ierror /= MPI_SUCCESS) CALL ut_abort(context//'MPI_COMM_RANK failed', &
         __FILE__, &
         __LINE__)
    IF (mype==0) THEN
      lroot = .true.
    ELSE
      lroot = .FALSE.
      verbose = .FALSE.
    ENDIF
    CALL set_verbose(verbose)
    CALL factorize(nprocs, nprocx, nprocy)
    IF (lroot .AND. verbose) WRITE(0,*) 'nprocx, nprocy=',nprocx, nprocy
    IF (lroot .AND. verbose) WRITE(0,*) 'g_ie, g_je=',g_ie, g_je
    mypy = mype / nprocx
    mypx = MOD(mype, nprocx)

    CALL deco
    ke = nlev

    ALLOCATE(g_id(g_ie, g_je), g_tpex(g_ie, g_je))

    ALLOCATE(fval_2d(ie,je), gval_2d(ie,je))
    ALLOCATE(loc_id_2d(ie,je), loc_tpex_2d(ie,je))
    ALLOCATE(id_pos(ie,je), pos3d_surf(ie,je))

    fval_2d = undef_int
    gval_2d = undef_int
    loc_id_2d = undef_int
    loc_tpex_2d = undef_int
    id_pos = undef_int
    pos3d_surf = undef_int

    CALL ut_init(decomp_size=30, comm_tmpl_size=30, comm_size=30, &
         &       debug_lvl=0, mode=ut_mode_dt_p2p, debug_unit=0)

  END SUBROUTINE init_all

  SUBROUTINE exchange2d_2d(itrans, f, g)
    INTEGER, INTENT(in) :: itrans
    INTEGER, TARGET, INTENT(in) :: f(:,:)
    INTEGER, TARGET, INTENT(out) :: g(:,:)

    INTEGER, POINTER :: p_in, p_out

    p_in => f(1,1)
    p_out=> g(1,1)

    CALL ut_transpose(p_in, itrans, comm_forward, p_out)

  END SUBROUTINE exchange2d_2d

  SUBROUTINE exchange3d_3d(itrans, f, g)
    INTEGER, INTENT(in) :: itrans
    INTEGER, TARGET, INTENT(in) :: f(:,:,:)
    INTEGER, TARGET, INTENT(out) :: g(:,:,:)

    INTEGER, POINTER :: p_in, p_out

    p_in => f(1,1,1)
    p_out=> g(1,1,1)

    CALL ut_transpose(p_in, itrans, comm_forward, p_out)

  END SUBROUTINE exchange3d_3d


  SUBROUTINE exchange2d_3d(itrans, f, g)
    INTEGER, INTENT(in) :: itrans
    INTEGER, TARGET, INTENT(in) :: f(:,:)
    INTEGER, TARGET, INTENT(out) :: g(:,:,:)

    INTEGER, POINTER :: p_in, p_out

    p_in => f(1,1)
    p_out=> g(1,1,1)

    CALL ut_transpose(p_in, itrans, comm_forward, p_out)

  END SUBROUTINE exchange2d_3d

  SUBROUTINE gen_trans(itemp, send_dt, recv_dt, itrans)
    INTEGER,INTENT(in) :: itemp, send_dt, recv_dt
    INTEGER,INTENT(out) :: itrans

    INTEGER :: dt

    IF (send_dt /= recv_dt) &
         CALL ut_abort('gen_trans: (send_dt /= recv_dt) unsupported', &
         __FILE__, &
         __LINE__)
    dt = send_dt
    CALL ut_init_transposition(itemp, dt, itrans)

  END SUBROUTINE gen_trans

  SUBROUTINE gen_off_trans(itemp, send_dt, send_off, recv_dt, recv_off, itrans)
    INTEGER,INTENT(in) :: itemp, send_dt, recv_dt
    INTEGER,INTENT(in) :: send_off(:,:), recv_off(:,:)
    INTEGER,INTENT(out) :: itrans

    INTEGER :: send_offsets(SIZE(send_off)), recv_offsets(SIZE(recv_off))

    send_offsets = RESHAPE(send_off, (/SIZE(send_off)/) )
    recv_offsets = RESHAPE(recv_off, (/SIZE(recv_off)/) )

    CALL ut_init_transposition(itemp, send_offsets, recv_offsets, &
         send_dt, recv_dt, itrans)

  END SUBROUTINE gen_off_trans

  SUBROUTINE get_window(gval, win)
    INTEGER, INTENT(in) :: gval(:,:)
    INTEGER, INTENT(out) :: win(:,:)

    INTEGER :: i, j, ig, jg

    DO j = 1, je
      jg = p_joff + j
      DO i = 1, ie
        ig = p_ioff + i
        win(i,j) =  gval(ig,jg)
      ENDDO
    ENDDO

  END SUBROUTINE get_window

  SUBROUTINE gen_template_2d(local_src_idx, local_dst_idx, ihandle)
    INTEGER, INTENT(in) :: local_src_idx(:,:)
    INTEGER, INTENT(in) :: local_dst_idx(:,:)
    INTEGER, INTENT(out) :: ihandle

    INTEGER :: src(SIZE(local_src_idx)), dst(size(local_dst_idx))

    INTEGER :: src_handle, dst_handle

    src = RESHAPE( local_src_idx, (/SIZE(local_src_idx)/) )
    dst = RESHAPE( local_dst_idx, (/SIZE(local_dst_idx)/) )

    CALL ut_init_decomposition(src, g_ie * g_je, src_handle)
    CALL ut_init_decomposition(dst, g_ie * g_je, dst_handle)

    CALL ut_init_oneway_transposition_template(src_handle, dst_handle, &
         MPI_COMM_WORLD, ihandle)

    CALL ut_destroy_decomposition(dst_handle)
    CALL ut_destroy_decomposition(src_handle)

  END SUBROUTINE gen_template_2d

  SUBROUTINE gen_template_3d(local_src_idx, local_dst_idx, ihandle)
    INTEGER, INTENT(in) :: local_src_idx(:,:,:)
    INTEGER, INTENT(in) :: local_dst_idx(:,:,:)
    INTEGER, INTENT(out) :: ihandle

    INTEGER, ALLOCATABLE :: src(:), dst(:)
    INTEGER :: src_handle, dst_handle

    ALLOCATE(src(SIZE(local_src_idx)), dst(SIZE(local_dst_idx)))

    src = RESHAPE( local_src_idx, (/SIZE(local_src_idx)/) )
    dst = RESHAPE( local_dst_idx, (/SIZE(local_dst_idx)/) )

    CALL ut_init_decomposition(src, g_ie * g_je, src_handle)
    CALL ut_init_decomposition(dst, g_ie * g_je, dst_handle)

    CALL ut_init_oneway_transposition_template(src_handle, dst_handle, &
         MPI_COMM_WORLD, ihandle)

    CALL ut_destroy_decomposition(dst_handle)
    CALL ut_destroy_decomposition(src_handle)

  END SUBROUTINE gen_template_3d


  SUBROUTINE def_exchange()

    LOGICAL, PARAMETER :: increased_north_halo = .FALSE.
    LOGICAL, PARAMETER :: with_north_halo = .true.
    INTEGER :: i, j
    INTEGER :: g_core_is, g_core_ie, g_core_js, g_core_je
    INTEGER :: north_halo

    ! global core domain:
    g_core_is = nhalo + 1
    g_core_ie = g_ie-nhalo
    g_core_js = nhalo + 1
    g_core_je = g_je-nhalo

    ! global tripolar boundsexchange:
    g_tpex = undef_index
    g_tpex(g_core_is:g_core_ie, g_core_js:g_core_je) &
         = g_id(g_core_is:g_core_ie, g_core_js:g_core_je)

    IF (with_north_halo) THEN

      ! north inversion, (maybe with increased north halo)
      IF (increased_north_halo) THEN
        north_halo = nhalo+1
      ELSE
        north_halo = nhalo
      ENDIF

      IF (2*north_halo > g_core_je) &
           CALL ut_abort('def_exchange: grid too small (or halo too large) &
           &for tripolar north exchange', &
           __FILE__, &
           __LINE__)
      DO j = 1, north_halo
        DO i = g_core_is, g_core_ie
          g_tpex(i,j) = g_tpex(g_core_ie + (g_core_is-i), 2*north_halo + (1-j))
        ENDDO
      ENDDO

    ELSE

      DO j = 1, nhalo
        DO i = nhalo+1, g_ie-nhalo
          g_tpex(i,j) = g_id(i,j)
        ENDDO
      ENDDO

    ENDIF

    ! south: no change
    DO j = g_core_je+1, g_je
      DO i = nhalo+1, g_ie-nhalo
        g_tpex(i,j) = g_id(i,j)
      ENDDO
    ENDDO

    ! PBC
    DO j = 1, g_je
      DO i = 1, nhalo
        g_tpex(g_core_is-i,j) = g_tpex(g_core_ie+(1-i),j)
      ENDDO
      DO i = 1, nhalo
        g_tpex(g_core_ie+i,j) = g_tpex(nhalo+i,j)
      ENDDO
    ENDDO

    CALL check_g_idx (g_tpex)

  END SUBROUTINE def_exchange

  SUBROUTINE check_g_idx(gidx)
    INTEGER,INTENT(in) :: gidx(:,:)

    IF (ANY(gidx == undef_index)) THEN
      CALL ut_abort('check_g_idx: check failed', __FILE__, __LINE__)
    ENDIF
  END SUBROUTINE check_g_idx

  SUBROUTINE deco
    INTEGER :: cx0(0:nprocx-1), cxn(0:nprocx-1)
    INTEGER :: cy0(0:nprocy-1), cyn(0:nprocy-1)

    cx0 = 0
    cxn = 0
    CALL regular_deco(g_ie-2*nhalo, cx0, cxn)

    cy0 = 0
    cyn = 0
    CALL regular_deco(g_je-2*nhalo, cy0, cyn)

    ! process local deco variables:
    ie = cxn(mypx) + 2*nhalo
    je = cyn(mypy) + 2*nhalo
    p_ioff = cx0(mypx)
    p_joff = cy0(mypy)

  END SUBROUTINE deco

END PROGRAM test_perf
