!>
!> @file test_yaxt.f90
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

PROGRAM test_yaxt
  USE iso_c_binding, ONLY: c_loc, c_int
  USE mpi
  USE yaxt, ONLY: xt_idxlist, xt_idxlist_delete, xt_idxvec_new, &
       &          xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       &          xt_redist, xt_redist_p2p_new, xt_redist_delete, &
       &          xt_initialize, xt_finalize, xt_idxlist_get_num_indices, &
       &          xt_idxlist_get_indices, xt_modifier, xt_idxmod_new, &
       xt_redist_p2p_off_new, xt_redist_s_exchange1, xt_int_kind, &
       xt_idxfsection_new, xt_redist_collection_static_new
  USE ftest_common, ONLY: test_abort, icmp, id_map, factorize, regular_deco, &
       finish_mpi
  USE iso_c_binding, ONLY: c_loc, c_int

  IMPLICIT NONE

  INTEGER, PARAMETER :: g_ie = 8, g_je = 4! global extents including halos
  LOGICAL, PARAMETER :: verbose = .FALSE.
  INTEGER, PARAMETER :: nlev = 3
  INTEGER, PARAMETER :: undef_int = -1
  INTEGER(xt_int_kind), PARAMETER :: undef_index = -1
  INTEGER, PARAMETER :: nhalo = 1 ! 1dim. halo border size
  LOGICAL, PARAMETER :: increased_north_halo = .FALSE.
  LOGICAL, PARAMETER :: with_north_halo = .TRUE.

  INTEGER :: ie, je ! local extents, including halos
  INTEGER :: p_ioff, p_joff ! offsets within global domain
  INTEGER :: nprocx, nprocy ! process space extents
  INTEGER :: nprocs ! == nprocx*nprocy
  INTEGER :: mype, mypx, mypy ! process rank, process coords within (0:, 0:) process space
  LOGICAL :: lroot ! true only for proc 0

  INTEGER(xt_int_kind) :: g_id(g_ie, g_je) ! global id

  ! global "tripolar-like" toy bounds exchange
  INTEGER(xt_int_kind) :: g_tpex(g_ie, g_je)
  TYPE(xt_xmap) :: xmap_tpex
  TYPE(xt_redist) :: redist_tpex
  TYPE(xt_redist) :: redist_surf_tpex

  INTEGER(xt_int_kind), ALLOCATABLE :: loc_id(:,:), loc_tpex(:,:)
  INTEGER(c_int), ALLOCATABLE :: mstate(:,:)
  INTEGER, ALLOCATABLE :: fval(:,:), gval(:,:)
  INTEGER, ALLOCATABLE :: gval3d(:,:,:)
  INTEGER, ALLOCATABLE :: id_pos(:,:), pos3d_surf(:,:)

  ! mpi & decomposition & allocate mem:
  CALL init_all

  ! full global index space:
  CALL id_map(g_id)

  ! local window of global index space:
  CALL get_window(g_id, loc_id)

  ! define bounds exchange for full global index space
  CALL def_exchange(g_id, g_tpex)

  ! local window of global bounds exchange:
  CALL get_window(g_tpex, loc_tpex)

  ! check interface to idxsection:
  CALL general_fsection_test

  ! compare current index construction with modifier results
  CALL check_modifiers

  ! template: loc_id -> loc_tpex
  CALL gen_template(loc_id, loc_tpex, xmap_tpex) ! todo rename template to xmap

  ! transposition: loc_id:data -> loc_tpex:data
  CALL gen_trans(xmap_tpex, MPI_INTEGER, MPI_INTEGER, redist_tpex)

  ! test 2d-to-2d transposition:
  fval = INT(loc_id)
  CALL exchange_int(redist_tpex, fval, gval)

  CALL icmp('2d to 2d check', gval, INT(loc_tpex), mype)

  CALL check_redist_collection_static

  ! define positions of surface elements within (i,k,j) array
  CALL gen_id_pos(id_pos)
  CALL gen_id_pos(pos3d_surf)
  CALL gen_pos3d_surf(pos3d_surf)

  ! generate surface transposition:
  CALL gen_off_trans(xmap_tpex, MPI_INTEGER, INT(id_pos(:,:)) - 1, &
       MPI_INTEGER, INT(pos3d_surf(:,:)) - 1, redist_surf_tpex)

  ! 2d to surface boundsexchange:
  gval3d = -1
  CALL exchange_int(redist_surf_tpex, fval, gval3d)

  CALL icmp('surface check', gval3d(:,1,:), INT(loc_tpex), mype)
  ! check sub surface:
  CALL icmp('sub surface check', gval3d(:,2,:), INT(loc_tpex)*0-1, mype)

  ! cleanup:
  CALL xt_xmap_delete(xmap_tpex)

  CALL xt_redist_delete(redist_tpex)

  CALL xt_redist_delete(redist_surf_tpex)

  CALL xt_finalize()
  CALL finish_mpi

CONTAINS

  SUBROUTINE check_redist_collection_static
    INTEGER, PARAMETER :: nr = 2
    TYPE(xt_redist) :: rvec(nr), rcol
    INTEGER, TARGET :: f(ie,je,nr), g(ie,je,nr), ref_g(ie,je,nr)
    INTEGER(mpi_address_kind) :: f_addr(nr), g_addr(nr)
    INTEGER(mpi_address_kind) :: f_disp(nr), g_disp(nr)

    INTEGER :: ir, ierror
    rvec(:) = redist_tpex
    DO ir = 1, nr
      CALL MPI_GET_ADDRESS(f(1,1,ir), f_addr(ir), ierror)
      IF (ierror /= MPI_SUCCESS) CALL my_abort('MPI_GET_ADDRESS failed', __LINE__)
      CALL MPI_GET_ADDRESS(g(1,1,ir), g_addr(ir), ierror)
      IF (ierror /= MPI_SUCCESS) CALL my_abort('MPI_GET_ADDRESS failed', __LINE__)
      f_disp(ir) = f_addr(ir) - f_addr(1)
      g_disp(ir) = g_addr(ir) - g_addr(1)
    ENDDO

    rcol = xt_redist_collection_static_new(rvec, nr, f_disp, g_disp, MPI_COMM_WORLD)
    DO ir = 1, nr
      f(:,:,ir) = INT(loc_id) + (ir-1) * ie*je
    ENDDO

    ref_g = 0
    DO ir = 1, nr
      CALL xt_redist_s_exchange1(rvec(ir), C_LOC(f(1,1,ir)), C_LOC(ref_g(1,1,ir)));
    ENDDO

    g = 0
    CALL xt_redist_s_exchange1(rcol, C_LOC(f), C_LOC(g))
    IF (ANY(g /= ref_g)) CALL my_abort('(g /= ref_g)', __LINE__)
    CALL xt_redist_delete(rcol)

  END SUBROUTINE check_redist_collection_static

  SUBROUTINE check_modifiers()
    TYPE(xt_modifier) :: m_tpex(5)
    INTEGER :: m_tpex_num
    TYPE(xt_idxlist) :: loc_id_idxlist
    INTEGER(xt_int_kind) :: loc_tpex2(ie,je)
    TYPE(xt_idxlist) :: loc_tpex2_idxlist

    loc_id_idxlist = xt_idxvec_new(loc_id, SIZE(loc_id))

    ! use one simple modifier to define index transfer:
    CALL def_tpex_mod_via_idxvec(m_tpex, m_tpex_num)

    loc_tpex2_idxlist = xt_idxmod_new(loc_id_idxlist, m_tpex, m_tpex_num, mstate)
    loc_tpex2 = -1
    CALL xt_idxlist_get_indices(loc_tpex2_idxlist, loc_tpex2)
    IF (ANY(loc_tpex2 /= loc_tpex)) &
         CALL my_abort('idx copy does not match', __LINE__)
    CALL xt_idxlist_delete(loc_tpex2_idxlist)

    ! test call without mstate
    loc_tpex2_idxlist = xt_idxmod_new(loc_id_idxlist, m_tpex, m_tpex_num)
    loc_tpex2 = -1
    CALL xt_idxlist_get_indices(loc_tpex2_idxlist, loc_tpex2)
    IF (ANY(loc_tpex2 /= loc_tpex)) &
         CALL my_abort('idx copy does not match', __LINE__)
    CALL xt_idxlist_delete(loc_tpex2_idxlist)
    CALL delete_modifiers(m_tpex(1:m_tpex_num))

    ! use compact modifiers to define index transfer:
    CALL def_tpex_mod_via_sections(m_tpex, m_tpex_num)
    loc_tpex2_idxlist = xt_idxmod_new(loc_id_idxlist, m_tpex, m_tpex_num, mstate)
    loc_tpex2 = -1
    CALL xt_idxlist_get_indices(loc_tpex2_idxlist, loc_tpex2)

    IF (ANY(loc_tpex2 /= loc_tpex)) &
         CALL my_abort('idx copy does not match', __LINE__)
    CALL xt_idxlist_delete(loc_tpex2_idxlist)
    CALL delete_modifiers(m_tpex(1:m_tpex_num))

    ! cleanup:
    CALL xt_idxlist_delete(loc_id_idxlist)
  END SUBROUTINE check_modifiers

  SUBROUTINE delete_modifiers(m)
    TYPE(xt_modifier), INTENT(inout) :: m(:)

    INTEGER :: i

    DO i = 1, SIZE(m)
      CALL xt_idxlist_delete(m(i)%extract)
      CALL xt_idxlist_delete(m(i)%subst)
    ENDDO

  END SUBROUTINE delete_modifiers

  SUBROUTINE my_abort(msg, line)
    CHARACTER(*), INTENT(in) :: msg
    INTEGER, VALUE, INTENT(in) :: line
    CALL test_abort(msg, &
         __FILE__, &
         line)
  END SUBROUTINE my_abort

  SUBROUTINE general_fsection_test
    INTEGER(xt_int_kind), PARAMETER :: gdx = 10_xt_int_kind, gdy=5_xt_int_kind
    INTEGER, PARAMETER :: ldx = 4,  ldy=2
    INTEGER(xt_int_kind), PARAMETER :: gstart = 1
    INTEGER(xt_int_kind), PARAMETER :: gsize(2) = (/ gdx, gdy /)
    TYPE(xt_idxlist) :: global_section, local_section
    INTEGER(xt_int_kind) :: indices(gdx*gdy), lstart(2)
    INTEGER :: egis(gdx, gdy)
    INTEGER :: i, j, idx, p

    ! prepare explicit global index space
    idx = gstart - 1
    DO j = 1, gdy
      DO i = 1, gdx
        idx = idx + 1
        egis(i,j) = idx
      ENDDO
    ENDDO

    lstart = (/ 1_xt_int_kind, 1_xt_int_kind /)

    ! check case: local section == global section
    global_section = xt_idxfsection_new(gstart, gsize, INT(gsize), lstart)
    indices = -1
    CALL xt_idxlist_get_indices(global_section, indices)
    p = 0
    DO j = 1, gdy
      DO i = 1, gdx
        p = p + 1
        IF (egis(i,j) /= indices(p)) CALL my_abort('(1) bad indices', __LINE__)
      ENDDO
    ENDDO
    CALL xt_idxlist_delete(global_section)

    ! check case: simple subsection
    local_section = xt_idxfsection_new(gstart, gsize, (/ ldx, ldy /), lstart)
    indices = -1
    CALL xt_idxlist_get_indices(local_section, indices)
    p = 0
    DO j = 1, ldy
      DO i = 1, ldx
        p = p + 1
        IF (egis(i,j) /= indices(p)) CALL my_abort('(2) bad indices', __LINE__)
      ENDDO
    ENDDO
    CALL xt_idxlist_delete(local_section)

    ! check case: i-reverse subsection
    local_section = xt_idxfsection_new(gstart, gsize, &
         (/ -ldx, ldy /), lstart)
    indices = -1
    CALL xt_idxlist_get_indices(local_section, indices)
    p = 0
    DO j = 1, ldy
      DO i = ldx, 1, -1
        p = p + 1
        IF (egis(i,j) /= indices(p)) CALL my_abort('(3) bad indices', __LINE__)
      ENDDO
    ENDDO
    CALL xt_idxlist_delete(local_section)

    ! check case: j-reverse subsection
    local_section = xt_idxfsection_new(gstart, gsize, (/ ldx, -ldy /), lstart)
    indices = -1
    CALL xt_idxlist_get_indices(local_section, indices)
    p = 0
    DO j = ldy, 1, -1
      DO i = 1, ldx
        p = p + 1
        IF (egis(i,j) /= indices(p)) CALL my_abort('(4) bad indices', __LINE__)
      ENDDO
    ENDDO
    CALL xt_idxlist_delete(local_section)

    ! check case: ij-reverse subsection
    local_section = xt_idxfsection_new(gstart, gsize, &
         (/ -ldx, -ldy /), lstart)
    indices = -1
    CALL xt_idxlist_get_indices(local_section, indices)
    p = 0
    DO j = ldy, 1, -1
      DO i = ldx, 1, -1
        p = p + 1
        IF (egis(i,j) /= indices(p)) CALL my_abort('(5) bad indices', __LINE__)
      ENDDO
    ENDDO
    CALL xt_idxlist_delete(local_section)
  END SUBROUTINE general_fsection_test

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

  SUBROUTINE init_all
    CHARACTER(len=*), PARAMETER :: context = 'init_all: '
    INTEGER :: ierror

    CALL MPI_INIT(ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL my_abort(context//'MPI_INIT failed', __LINE__)

    CALL xt_initialize(MPI_COMM_WORLD)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL my_abort(context//'MPI_COMM_SIZE failed', __LINE__)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL my_abort(context//'MPI_COMM_RANK failed', __LINE__)
    IF (mype==0) THEN
      lroot = .true.
    ELSE
      lroot = .FALSE.
    ENDIF

    CALL factorize(nprocs, nprocx, nprocy)
    IF (verbose .AND. lroot) WRITE(0,*) 'nprocx, nprocy=',nprocx, nprocy
    mypy = mype / nprocx
    mypx = MOD(mype, nprocx)

    !CALL ut_init(decomp_size=30, comm_tmpl_size=30, comm_size=30, &
    !     &       debug_lvl=0, mode=ut_mode_dt_p2p, debug_unit=0)

    CALL deco

    ALLOCATE(fval(ie,je), gval(ie,je))
    ALLOCATE(loc_id(ie,je), loc_tpex(ie,je), mstate(ie,je))
    ALLOCATE(id_pos(ie,je), gval3d(ie,nlev,je), pos3d_surf(ie,je))

    fval = undef_int
    gval = undef_int
    loc_id = INT(undef_int, xt_int_kind)
    loc_tpex = INT(undef_int, xt_int_kind)
    id_pos = undef_int
    gval3d = undef_int
    pos3d_surf = undef_int

  END SUBROUTINE init_all

  SUBROUTINE gen_id_pos(pos)
    INTEGER, INTENT(out) :: pos(:,:)

    INTEGER :: i,j,p

    p = 0
    DO j = 1, SIZE(pos,2)
      DO i = 1, SIZE(pos,1)
        p = p + 1
        pos(i,j) = p
      ENDDO
    ENDDO

  END SUBROUTINE gen_id_pos

  SUBROUTINE exchange_int(redist, f, g)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, TARGET, INTENT(in) :: f(*)
    INTEGER, TARGET, VOLATILE, INTENT(out) :: g(*)

    CALL xt_redist_s_exchange1(redist, c_loc(f), c_loc(g));

  END SUBROUTINE exchange_int

  SUBROUTINE gen_trans(xmap, send_dt, recv_dt, redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER,INTENT(in) :: send_dt, recv_dt
    TYPE(xt_redist),INTENT(out) :: redist

    INTEGER :: dt

    IF (send_dt /= recv_dt) &
         CALL my_abort('gen_trans: (send_dt /= recv_dt) unsupported', __LINE__)
    dt = send_dt
    redist = xt_redist_p2p_new(xmap, dt)
    !CALL ut_init_transposition(itemp, dt, itrans)

  END SUBROUTINE gen_trans

  SUBROUTINE gen_off_trans(xmap, send_dt, send_off, recv_dt, recv_off, redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER,INTENT(in) :: send_dt, recv_dt
    INTEGER(c_int),INTENT(in) :: send_off(:,:), recv_off(:,:)
    TYPE(xt_redist),INTENT(out) :: redist

    !INTEGER :: send_offsets(SIZE(send_off)), recv_offsets(SIZE(recv_off))

    !send_offsets = RESHAPE(send_off, (/SIZE(send_off)/) )
    !recv_offsets = RESHAPE(recv_off, (/SIZE(recv_off)/) )
    IF (recv_dt /= send_dt) &
         CALL my_abort('(datatype_in /= datatype_out) not supported', &
         __LINE__)

    redist = xt_redist_p2p_off_new(xmap, send_off, recv_off, send_dt);
    !CALL ut_init_transposition(itemp, send_offsets, recv_offsets, send_dt, recv_dt, itrans)

  END SUBROUTINE gen_off_trans

  SUBROUTINE get_window(gval, win)
    INTEGER(xt_int_kind), INTENT(in) :: gval(:,:)
    INTEGER(xt_int_kind), INTENT(out) :: win(:,:)

    INTEGER :: i, j, ig, jg

    DO j = 1, je
      jg = p_joff + j
      DO i = 1, ie
        ig = p_ioff + i
        win(i,j) =  gval(ig,jg)
      ENDDO
    ENDDO

  END SUBROUTINE get_window

  SUBROUTINE gen_template(local_src_idx, local_dst_idx, xmap)
    INTEGER(xt_int_kind), INTENT(in) :: local_src_idx(:,:)
    INTEGER(xt_int_kind), INTENT(in) :: local_dst_idx(:,:)
    TYPE(xt_xmap), INTENT(out) :: xmap

    TYPE(Xt_idxlist) :: src_idxlist, dst_idxlist
    INTEGER :: src_num, dst_num
    INTEGER(xt_int_kind) :: cp_src_idx(g_ie, g_je)
    INTEGER(xt_int_kind) :: cp_dst_idx(g_ie, g_je)

    src_idxlist = xt_idxvec_new(local_src_idx, g_ie * g_je)
    src_num = INT(xt_idxlist_get_num_indices(src_idxlist))
    IF (src_num /= g_ie*g_je) CALL my_abort('unexpected src_num', __LINE__)
    CALL xt_idxlist_get_indices(src_idxlist, cp_src_idx)
    IF (ANY(cp_src_idx /= local_src_idx)) CALL my_abort('idx copy does not match', &
         __LINE__)

    dst_idxlist = xt_idxvec_new(local_dst_idx, g_ie * g_je)
    dst_num = INT(xt_idxlist_get_num_indices(dst_idxlist))
    IF (dst_num /= g_ie*g_je) CALL my_abort('unexpected dst_num', __LINE__)
    CALL xt_idxlist_get_indices(dst_idxlist, cp_dst_idx)
    IF (ANY(cp_dst_idx /= local_dst_idx)) CALL my_abort('idx copy does not match', &
         __LINE__)

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist,  MPI_COMM_WORLD)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

  END SUBROUTINE gen_template

  SUBROUTINE def_tpex_mod_via_idxvec(mvec, mvec_num)
    TYPE(xt_modifier), INTENT(out) :: mvec(:)
    INTEGER, INTENT(out) :: mvec_num

    INTEGER(xt_int_kind) :: g_start_indices(g_ie, g_je)
    INTEGER(xt_int_kind) :: g_end_indices(g_ie, g_je)
    TYPE(xt_idxlist) :: g_start_idxlist
    TYPE(xt_idxlist) :: g_end_idxlist

    IF (SIZE(mvec)<1) CALL my_abort('def_tpex_mod_via_idxvec mvec too small', &
         __LINE__)

    CALL id_map(g_start_indices)
    g_start_idxlist = xt_idxvec_new(g_start_indices, SIZE(g_start_indices))

    CALL def_exchange(g_start_indices, g_end_indices)
    g_end_idxlist = xt_idxvec_new(g_end_indices, SIZE(g_end_indices))

    mvec(1)%extract = g_start_idxlist
    mvec(1)%subst = g_end_idxlist
    mvec(1)%mask = 1
    mvec_num = 1

  END SUBROUTINE def_tpex_mod_via_idxvec

  SUBROUTINE def_tpex_mod_via_sections(mvec, mvec_num)
    TYPE(xt_modifier), INTENT(out) :: mvec(:)
    INTEGER, INTENT(out) :: mvec_num

    INTEGER(xt_int_kind), PARAMETER :: gstart_idx = 1_xt_int_kind
    INTEGER(xt_int_kind), PARAMETER :: gsize(2) &
         = (/ INT(g_ie, xt_int_kind), INT(g_je, xt_int_kind) /)
    INTEGER :: ldx, ldy

    INTEGER(xt_int_kind) :: g_core_is, g_core_ie, g_core_je
    INTEGER(xt_int_kind) :: north_halo, im

    ! global core domain:
    g_core_is = nhalo + 1
    g_core_ie = g_ie-nhalo
    g_core_je = g_je-nhalo

    im = 0

    ! global tripolar boundsexchange:
    IF (with_north_halo) THEN

      ! north inversion, (maybe with increased north halo)
      IF (increased_north_halo) THEN
        north_halo = nhalo+1
      ELSE
        north_halo = nhalo
      ENDIF

      IF (2*north_halo > g_core_je) &
           CALL my_abort('def_tpex_mod_via_sections: grid too small (or halo too large)' // &
           'for tripolar north exchange', __LINE__)

      im = im + 1_xt_int_kind
      IF (SIZE(mvec)<im) CALL my_abort('(SIZE(mvec)<im)', __LINE__)
      ! north border exchange without ew-halos
      ldx = INT(g_core_ie - g_core_is + 1)
      ldy = INT(north_halo)
      mvec(im)%extract = xt_idxfsection_new(gstart_idx, gsize, &
           (/ ldx, ldy /),   (/g_core_is, 1_xt_int_kind/))
      mvec(im)%subst   = xt_idxfsection_new(gstart_idx, gsize, &
           (/ -ldx, -ldy /), (/g_core_is, north_halo+1_xt_int_kind/))
      mvec(im)%mask    = 1

      ! 1. north edge:
      im = im + 1_xt_int_kind
      IF (SIZE(mvec)<im) CALL my_abort('(SIZE(mvec)<im)', __LINE__)
      ldx = 1
      ldy = INT(north_halo)
      mvec(im)%extract = xt_idxfsection_new(gstart_idx, gsize, &
           (/ ldx, ldy /),   (/1_xt_int_kind, 1_xt_int_kind/))
      mvec(im)%subst   = xt_idxfsection_new(gstart_idx, gsize, &
           (/ -ldx, -ldy /), (/2_xt_int_kind, north_halo+1_xt_int_kind/))
      mvec(im)%mask    = 1
      ! 2. north edge:
      im = im + 1_xt_int_kind
      IF (SIZE(mvec)<im) CALL my_abort('(SIZE(mvec)<im)', __LINE__)
      ldx = 1
      ldy = INT(north_halo)
      mvec(im)%extract = xt_idxfsection_new(gstart_idx, gsize, &
           (/ ldx, ldy /),   (/INT(g_ie, xt_int_kind), 1_xt_int_kind/))
      mvec(im)%subst   = xt_idxfsection_new(gstart_idx, gsize, &
           (/ -ldx, -ldy /), &
           (/INT(g_ie - 1, xt_int_kind), north_halo+1_xt_int_kind/))
      mvec(im)%mask    = 1

    ELSE

      ! nothing to do at the north border

    ENDIF

    ! PBC below north border
    ldx = nhalo
    ldy = INT(INT(g_je, xt_int_kind) - north_halo)
    im = im + 1_xt_int_kind
    IF (SIZE(mvec)<im) CALL my_abort('(SIZE(mvec)<im)', __LINE__)
    mvec(im)%extract = xt_idxfsection_new(gstart_idx, gsize, &
         (/ ldx, ldy /), (/1_xt_int_kind, north_halo+1_xt_int_kind/))
    mvec(im)%subst   = xt_idxfsection_new(gstart_idx, gsize, &
         (/ ldx, ldy /), &
         (/ g_core_ie - INT(ldx, xt_int_kind) + 1_xt_int_kind, &
           north_halo + 1_xt_int_kind/))
    mvec(im)%mask    = 1

    im = im + 1_xt_int_kind
    IF (SIZE(mvec)<im) CALL my_abort('(SIZE(mvec)<im)', __LINE__)
    mvec(im)%extract = xt_idxfsection_new(gstart_idx, gsize, &
         (/ ldx, ldy /), (/g_core_ie+1_xt_int_kind, north_halo+1_xt_int_kind/))
    mvec(im)%subst   = xt_idxfsection_new(gstart_idx, gsize, &
         (/ ldx, ldy /), (/ INT(ldx, xt_int_kind) + 1_xt_int_kind, &
         &                  north_halo+1_xt_int_kind/))
    mvec(im)%mask    = 1

    mvec_num = INT(im)

  END SUBROUTINE def_tpex_mod_via_sections

  SUBROUTINE def_exchange(id_in, id_out)
    INTEGER(xt_int_kind), INTENT(in) :: id_in(:,:)
    INTEGER(xt_int_kind), INTENT(out) :: id_out(:,:)

    INTEGER :: i, j
    INTEGER :: g_core_is, g_core_ie, g_core_js, g_core_je
    INTEGER :: north_halo

    ! global core domain:
    g_core_is = nhalo + 1
    g_core_ie = g_ie-nhalo
    g_core_js = nhalo + 1
    g_core_je = g_je-nhalo

    ! global tripolar boundsexchange:
    id_out = undef_index
    id_out(g_core_is:g_core_ie, g_core_js:g_core_je) &
         = id_in(g_core_is:g_core_ie, g_core_js:g_core_je)

    IF (with_north_halo) THEN

      ! north inversion, (maybe with increased north halo)
      IF (increased_north_halo) THEN
        north_halo = nhalo+1
      ELSE
        north_halo = nhalo
      ENDIF

      IF (2*north_halo > g_core_je) &
           CALL my_abort('def_exchange: grid too small (or halo too large)' // &
           'for tripolar north exchange', __LINE__)
      DO j = 1, north_halo
        DO i = g_core_is, g_core_ie
          id_out(i,j) = id_out(g_core_ie + (g_core_is-i), 2*north_halo + (1-j))
        ENDDO
      ENDDO

    ELSE

      DO j = 1, nhalo
        DO i = nhalo+1, g_ie-nhalo
          id_out(i,j) = id_in(i,j)
        ENDDO
      ENDDO

    ENDIF

    ! south: no change
    DO j = g_core_je+1, g_je
      DO i = nhalo+1, g_ie-nhalo
        id_out(i,j) = id_in(i,j)
      ENDDO
    ENDDO

    ! PBC
    DO j = 1, g_je
      DO i = 1, nhalo
        id_out(g_core_is-i,j) = id_out(g_core_ie+(1-i),j)
      ENDDO
      DO i = 1, nhalo
        id_out(g_core_ie+i,j) = id_out(nhalo+i,j)
      ENDDO
    ENDDO

    CALL check_g_idx(id_out)

  END SUBROUTINE def_exchange

  SUBROUTINE check_g_idx(gidx)
    INTEGER(xt_int_kind), INTENT(in) :: gidx(:,:)

    IF (ANY(gidx == undef_index)) THEN
      CALL my_abort('check_g_idx: check failed', __LINE__)
    ENDIF
  END SUBROUTINE check_g_idx

  SUBROUTINE deco
    INTEGER :: cx0(0:nprocx-1), cxn(0:nprocx-1)
    INTEGER :: cy0(0:nprocy-1), cyn(0:nprocy-1)

    CALL regular_deco(g_ie-2*nhalo, cx0, cxn)
    CALL regular_deco(g_je-2*nhalo, cy0, cyn)

    ! process local deco variables:
    ie = cxn(mypx) + 2*nhalo
    je = cyn(mypy) + 2*nhalo
    p_ioff = cx0(mypx)
    p_joff = cy0(mypy)

  END SUBROUTINE deco

END PROGRAM test_yaxt
