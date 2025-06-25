!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide YAXT support.
!!
!! @remarks
!!   This module contains routines for description of domain
!!   decomposition for parallel output with CDI-PIO and global data
!!   transpositions using YAXT.
!!
!! @author Jörg Behrens, DKRZ, Hamburg (2013-06-20): Initial version
!!         Irina Fast,   DKRZ, Hamburg (2013-07-18): Revised version
!!
MODULE mo_echam_yaxt
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_loc, c_null_ptr, c_ptr
  USE mo_kind,          ONLY: dp
  USE mo_mpi,           ONLY: p_real_dp, p_pe, p_io
  USE mo_decomposition, ONLY: ldc => local_decomposition, debug_parallel

  USE mo_exception,     ONLY: message, message_text, finish
  USE mo_control,       ONLY: lyaxt_transposition
#ifdef HAVE_YAXT
  USE mo_control,       ONLY: control_nsp => nsp
  USE mo_advection,     ONLY: advection_jps => jps
  USE mo_tracdef,       ONLY: tracdef_trlist => trlist, tracdef_ntrac => ntrac
  USE yaxt,             ONLY: xt_initialize, xt_finalize, xt_abort,       &
                              xt_idxlist, xt_idxvec_new, xt_idxempty_new, &
                              xt_xmap, xt_xmap_all2all_new, xt_redist, &
                              xt_redist_p2p_off_new, xt_redist_s_exchange1, &
                              xt_redist_collection_new, xt_redist_s_exchange, &
                              xt_stripe, xt_idxstripes_new, xt_is_null, &
                              xt_idxlist_c2f, xt_xmap_delete, &
                              xt_mpi_comm_mark_exclusive
#endif
  USE mo_transpose,     ONLY: tr_gp_ffsl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: yaxt_initialize, yaxt_finalize, setup_yaxt_decomposition, &
            generate_yaxt_redist, add_yaxt_gp_nlevs
  PUBLIC :: yaxt_tr_gp_to_ffsl, yaxt_tr_ffsl_to_gp, yaxt_tr_tracer_gp_to_ffsl, &
       &    yaxt_tr_gp_to_fs, yaxt_tr_fs_to_gp, &
       &    yaxt_tr_fs_to_ls, yaxt_tr_ls_to_fs, &
       &    yaxt_tr_ls_to_sp, yaxt_tr_sp_to_ls

#ifdef HAVE_YAXT
  PUBLIC :: gp_gdeco, sp_sdeco
#endif


! switch for runtime-selftest
!#define YAXT_SELFTEST

! check memory layout expectations of ECHAM data
!#define CHECK_MEM_STRIDES


  ! methods for generation of index lists
  INTEGER, PARAMETER :: idx_descr_vector  = 1
  INTEGER, PARAMETER :: idx_descr_stripes = 2
  INTEGER, PARAMETER :: idx_descr_section = 3

#ifdef HAVE_YAXT

  ! agregation switches (aggregation leads to bigger and fewer messages):
  LOGICAL, PARAMETER :: aggregate_default    = .TRUE.
  LOGICAL, PARAMETER :: aggregate_gp_to_fs   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_fs_to_gp   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_fs_to_ls   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_ls_to_fs   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_ls_to_sp   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_sp_to_ls   = aggregate_default
  LOGICAL, PARAMETER :: aggregate_gp_to_ffsl = aggregate_default
  LOGICAL, PARAMETER :: aggregate_ffsl_to_gp = aggregate_default

  ! tie each decomposition (idxlist) in echam to a fixed data layout (offsets):
  TYPE yaxt_deco
    TYPE(xt_idxlist)     :: idxlist
    INTEGER, ALLOCATABLE :: offset(:)
  END TYPE yaxt_deco

  ! Note:
  ! We distinguish 4 kinds of index spaces: gdeco, sdeco, zmdeco, and m0deco.
  ! The general decomposition naming is: {short decomposition name}_{index space name}
  ! Only decompositions belonging to the same index space are combined within a redist object.
  ! The decompositions are tied to offsets as used in ECHAM.
  TYPE(yaxt_deco), ALLOCATABLE, SAVE :: gp_gdeco(:) ! standard gridpoint decomposition of grid point data
  TYPE(yaxt_deco), ALLOCATABLE, SAVE :: sp_sdeco(:) ! fourier-decomposition of grid point data
  TYPE(yaxt_deco), ALLOCATABLE, SAVE :: fs_gdeco(:) ! fourier-decomposition of grid point data
  TYPE(yaxt_deco), ALLOCATABLE, SAVE :: ls_sdeco(:) ! legendre decomposition of spectral point data

  ! standard gridpoint decomposition for tracers:
  TYPE(yaxt_deco), SAVE:: gp_4d_fb_gdeco ! offsets as in mo_tpcore::fb_gp
  TYPE(yaxt_deco), SAVE:: gp_4d_xtm1_gdeco ! offsets as in mo_tpcore::xtm1

  ! FFSL decompositions:
  TYPE(yaxt_deco), SAVE :: ffsl_2d_gdeco ! decomposition description for FFSL 2d
  TYPE(yaxt_deco), SAVE :: ffsl_3d_gdeco ! decomposition description for FFSL 3d
  TYPE(yaxt_deco), SAVE :: ffsl_4d_fb_gdeco ! deco with offsets as in mo_tpcore::fb
  TYPE(yaxt_deco), SAVE :: ffsl_4d_tcfb_gdeco ! like ffsl_4d_fb_gdeco but only trlist_count tracers

  ! subset of fs_gdeco(nlev+1), used for fs <--> ls:
  TYPE(yaxt_deco), SAVE :: fs_subset_gdeco

  ! 2d FS decomposition:
  TYPE(yaxt_deco), SAVE :: fs_2d_gdeco

  ! decomposition of GP zonal means for all levels:
  TYPE(yaxt_deco), SAVE :: gp_zmdeco

  ! decomposition of GP zonal means for selected levels (fs levels) only
  TYPE(yaxt_deco), SAVE :: gp_sel_zmdeco

  ! decomposition of FS zonal means:
  TYPE(yaxt_deco), SAVE :: fs_zmdeco

  ! 3d LS decomposition:
  TYPE(yaxt_deco), SAVE :: ls_3d_gdeco

  ! decomposition of LS zonal means:
  TYPE(yaxt_deco), SAVE :: ls_zmdeco

  ! decomposition of LS coeff. with m = 0:
  TYPE(yaxt_deco), SAVE :: ls_m0deco

  ! decomposition of SP coeff. with m = 0:
  TYPE(yaxt_deco), SAVE :: sp_m0deco

#endif

  ! global dimensions
  INTEGER :: nlon, nlat, nlev

  ! gp-deco info:
  INTEGER :: nproca, nprocb, nproma, ngpblks
  INTEGER :: glons(2), glone(2), glats(2), glate(2)

  ! ffsl-deco info:
  INTEGER :: ffsl_nlat, ffsl_nlev
  INTEGER :: ffsl_gp_lat1, ffsl_gp_lat2 ! ffsl lat-borders in gp-oriented coords

  ! fs-deco info:
  INTEGER :: flats(2), flate(2), nflat, flevs, fleve, nflev, nflevp1

  ! tracer info:
  INTEGER :: pcnst, trdim_size, active_trnum
  INTEGER, ALLOCATABLE :: active_trpos(:) ! position of active tracer within trlist

#ifdef HAVE_YAXT
  ! single var redists:

  ! gp -> ffsl:
  TYPE(xt_redist), SAVE :: gp2ffsl_4d_xtm1_redist
  TYPE(xt_redist), SAVE :: gp2ffsl_3d_redist
  TYPE(xt_redist), SAVE :: gp2ffsl_2d_redist
  TYPE(xt_redist), SAVE :: gp2ffsl_all_redist

  ! ffsl -> gp:
  TYPE(xt_redist), SAVE :: ffsl2gp_4d_fb_redist
  TYPE(xt_redist), SAVE :: ffsl2gp_3d_redist
  TYPE(xt_redist), SAVE :: ffsl2gp_2d_redist
  TYPE(xt_redist), SAVE :: ffsl2gp_all_redist

  ! gp -> fs:
  TYPE(xt_redist), SAVE :: gp2fs_3d_redist
  TYPE(xt_redist), SAVE :: gp2fs_2d_redist
  TYPE(xt_redist), SAVE :: gp2fs_zm_redist
  TYPE(xt_redist), SAVE :: gp2fs_all_redist

  ! fs -> gp:
  TYPE(xt_redist), SAVE :: fs2gp_3d_redist
  TYPE(xt_redist), SAVE :: fs2gp_2d_redist
  TYPE(xt_redist), SAVE :: fs2gp_zm_redist
  TYPE(xt_redist), SAVE :: fs2gp_all_redist

  ! fs -> ls:
  TYPE(xt_redist), SAVE :: fs2ls_3d_redist
  TYPE(xt_redist), SAVE :: fs2ls_zm_redist
  TYPE(xt_redist), SAVE :: fs2ls_all_redist

  ! ls -> fs:
  TYPE(xt_redist), SAVE :: ls2fs_3d_redist
  TYPE(xt_redist), SAVE :: ls2fs_zm_redist
  TYPE(xt_redist), SAVE :: ls2fs_all_redist

  ! ls -> sp:
  TYPE(xt_redist), ALLOCATABLE, SAVE :: ls2sp_3d_redist(:)
  TYPE(xt_redist), SAVE :: ls2sp_m0_redist
  TYPE(xt_redist), SAVE :: ls2sp_all_redist

  ! sp -> ls:
  TYPE(xt_redist), ALLOCATABLE, SAVE :: sp2ls_3d_redist(:)
  TYPE(xt_redist), SAVE :: sp2ls_m0_redist
  TYPE(xt_redist), SAVE :: sp2ls_all_redist

#ifdef CHECK_MEM_STRIDES
  INTERFACE check_mem_strides
    MODULE PROCEDURE check_mem_strides_4d
    MODULE PROCEDURE check_mem_strides_3d
    MODULE PROCEDURE check_mem_strides_2d
  END INTERFACE check_mem_strides
#endif

#endif /* HAVE_YAXT */

  ! nproca x nprocb domain:
  INTEGER :: ab_comm, ab_rank, ab_size

  LOGICAL, SAVE :: decos_are_finalized = .FALSE.
  LOGICAL, SAVE, ALLOCATABLE :: gp_nlevs(:)

CONTAINS

  SUBROUTINE yaxt_initialize(comm)
    INTEGER, INTENT(IN) :: comm     ! mpi communicator
#ifndef HAVE_YAXT
  IF (lyaxt_transposition) THEN
    WRITE (message_text, *) 'For array transpositions with YAXT '// &
       &'specification of -DHAVE_YAXT is required. Please, recompile '// &
       &'and link with YAXT libary OR set lyaxt_transposition=.FALSE. '// &
       &'in namelist PARCTL.'
    CALL finish('yaxt_initialize', message_text)
  ENDIF
#else
  CALL xt_initialize(comm)
#endif
  END SUBROUTINE yaxt_initialize

!-------------------------------------------------------------------------------

  SUBROUTINE yaxt_finalize()
#ifdef HAVE_YAXT
    CALL xt_finalize()
#endif
  END SUBROUTINE yaxt_finalize

!-------------------------------------------------------------------------------


  SUBROUTINE add_yaxt_gp_nlevs(levels)
    INTEGER, INTENT(in) :: levels(:)
#ifdef HAVE_YAXT
    INTEGER :: min_lev, nlevs_size_new, nlevs_size_old, i
    LOGICAL, ALLOCATABLE :: tmp(:)

    ! For every 3d-array that is to be written via mo_output, this subroutine must
    ! have been called with the size of third dimension as element of the argument vector.
    ! Multiple calls are possible but must precede the call of setup_yaxt_decomposition.
    ! This routine has to be called at least once if setup_yaxt_decomposition is called


    IF (decos_are_finalized) THEN
      WRITE (message_text, *) 'Late call of add_yaxt_gp_nlevs. Yaxt decos are already finalized.'
      CALL finish('mo_echam_yaxt: add_yaxt_gp_nlevs', 'Decos are already finalized.')
    ENDIF

    min_lev = MINVAL(levels)
    IF (min_lev<1) THEN
      WRITE (message_text, *) 'invalid level',MINVAL(levels)
      CALL finish('mo_echam_yaxt: add_yaxt_gp_nlevs', message_text)
    ENDIF

    nlevs_size_new = MAXVAL(levels)
    IF ( .NOT. ALLOCATED(gp_nlevs)) THEN
      ! alloc:
      ALLOCATE(gp_nlevs(nlevs_size_new))
      gp_nlevs = .FALSE.
    END IF

    nlevs_size_old = SIZE(gp_nlevs)
    IF (nlevs_size_new > nlevs_size_old) THEN
      ! realloc:
      ALLOCATE(tmp(nlevs_size_old))
      tmp = gp_nlevs
      DEALLOCATE(gp_nlevs)
      ALLOCATE(gp_nlevs(nlevs_size_new))
      gp_nlevs(1:nlevs_size_old) = tmp
      gp_nlevs(nlevs_size_old+1:) = .FALSE.
    END IF

    DO i = 1, size(levels)
      gp_nlevs(levels(i)) = .TRUE.
    ENDDO
#endif
  END SUBROUTINE add_yaxt_gp_nlevs

  SUBROUTINE setup_yaxt_decomposition()
#ifdef HAVE_YAXT
    INTEGER :: itracer, nt

    ! Gaussian grid decomposition data
    nproca  = ldc%nproca
    nprocb  = ldc%nprocb
    nproma  = ldc%nproma
    ngpblks = ldc%ngpblks

    nlon  = ldc%nlon
    nlat  = ldc%nlat
    nlev  = ldc%nlev

    glons = ldc%glons
    glone = ldc%glone
    glats = ldc%glats
    glate = ldc%glate


    ! FFSL decomposition data
    ! translate ffsl region boundaries to gp latitudes
    ffsl_gp_lat1 = nlat + 1 - ldc%ffsl%latn
    ffsl_gp_lat2 = nlat + 1 - ldc%ffsl%lats
    ffsl_nlat = ffsl_gp_lat2 - ffsl_gp_lat1 + 1

    ! FS decomposition data
    flats = ldc%flats
    flate = ldc%flate
    nflat = ldc%nflat
    flevs = ldc%flevs
    fleve = ldc%fleve
    nflev = ldc%nflev
    nflevp1 = ldc%nflevp1

    ! tracer data:
    trdim_size = tracdef_ntrac
    active_trnum = COUNT(tracdef_trlist%ti(1:tracdef_trlist%ntrac)%ntran /= 0)
    pcnst = advection_jps + active_trnum
    ALLOCATE(active_trpos(active_trnum))
    nt = 0
    DO itracer = 1, tracdef_trlist%ntrac
      IF (tracdef_trlist%ti(itracer)%ntran /= 0) THEN
        nt = nt + 1
        active_trpos(nt) = itracer
      ENDIF
    ENDDO

    ! define index lists and offsets for Gaussian grid
    CALL setup_gp_gdeco(idx_descr_stripes)

    ! spectral coeff. space decompositions:
    CALL setup_sp_sdeco(idx_descr_stripes)
    CALL setup_ls_sdeco

    CALL setup_fs_2d_gdeco
    CALL setup_fs_3d_gdeco
    CALL setup_ls_3d_gdeco

    ! define index lists and offsets for FFSL 2D data
    CALL setup_ffsl_2d_gdeco(idx_descr_stripes)
    ! define index lists and offsets for FFSL 3D data
    CALL setup_ffsl_3d_gdeco(idx_descr_stripes)
    CALL setup_ffsl_4d_gdeco

    ! zonal-means space decompositions:
    CALL setup_gp_zmdeco
    CALL setup_fs_zmdeco
    CALL setup_ls_zmdeco

    ! spectral m==0 coeff. space decompositions:
    CALL setup_sp_m0deco
    CALL setup_ls_m0deco

    decos_are_finalized = .TRUE.
#endif /* HAVE_YAXT */
  END SUBROUTINE setup_yaxt_decomposition

  SUBROUTINE generate_yaxt_redist(model_comm)
    INTEGER, INTENT(IN) :: model_comm  ! communicator of application without I/O servers
#ifdef HAVE_YAXT
    INTEGER, PARAMETER :: max_num_redists = 15
    TYPE(xt_redist) :: redists(max_num_redists)
    INTEGER, PARAMETER :: redist_cs = 1 ! best setting for static aggregation
    INTEGER :: num_redists, ir, n

    CALL set_comms(model_comm)

    ! generate gp <--> ffsl redistributions:

    gp2ffsl_2d_redist = new_redist(gp_gdeco(1), ffsl_2d_gdeco)
    ffsl2gp_2d_redist = new_redist(ffsl_2d_gdeco, gp_gdeco(1))
    gp2ffsl_3d_redist = new_redist(gp_gdeco(nlev), ffsl_3d_gdeco)
    ffsl2gp_3d_redist = new_redist(ffsl_3d_gdeco, gp_gdeco(nlev))

    gp2ffsl_3d_redist = new_redist(gp_gdeco(nlev), ffsl_3d_gdeco)
    ffsl2gp_3d_redist = new_redist(ffsl_3d_gdeco, gp_gdeco(nlev))
    ffsl2gp_4d_fb_redist = new_redist(ffsl_4d_fb_gdeco, gp_4d_fb_gdeco)
    IF (active_trnum>0) gp2ffsl_4d_xtm1_redist = new_redist(gp_4d_xtm1_gdeco, ffsl_4d_tcfb_gdeco)

    IF (aggregate_gp_to_ffsl) THEN
      redists(1:5) = gp2ffsl_3d_redist
      redists(6:7) = gp2ffsl_2d_redist
      num_redists = 7
      gp2ffsl_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    IF (aggregate_ffsl_to_gp) THEN
      redists(1) = ffsl2gp_4d_fb_redist
      redists(2) = ffsl2gp_2d_redist
      redists(3) = ffsl2gp_3d_redist
      num_redists = 3
      ffsl2gp_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    ! generate gp <--> fs redistributions:
    gp2fs_3d_redist = new_redist(gp_gdeco(nlev), fs_gdeco(nlev))
    fs2gp_3d_redist = new_redist(fs_gdeco(nlev), gp_gdeco(nlev))

    gp2fs_2d_redist = new_redist(gp_gdeco(1), fs_2d_gdeco)
    fs2gp_2d_redist = new_redist(fs_2d_gdeco, gp_gdeco(1))

    gp2fs_zm_redist = new_redist(gp_sel_zmdeco, fs_zmdeco)
    fs2gp_zm_redist = new_redist(fs_zmdeco, gp_zmdeco)

    IF (aggregate_gp_to_fs) THEN
      DO ir = 1, 7
        redists(ir) = gp2fs_3d_redist
      ENDDO
      redists(8) = gp2fs_2d_redist
      DO ir = 9, 11
        redists(ir) = gp2fs_zm_redist
      ENDDO
      num_redists = 11
      gp2fs_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    IF (aggregate_fs_to_gp) THEN
      DO ir = 1, 9
        redists(ir) = fs2gp_3d_redist
      ENDDO
      DO ir = 10, 12
        redists(ir) = fs2gp_2d_redist
      ENDDO
      DO ir = 13, 15
        redists(ir) = fs2gp_zm_redist
      ENDDO
      num_redists = 15
      fs2gp_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    ! generate fs <--> ls redistributions:
    fs2ls_3d_redist = new_redist(fs_subset_gdeco, ls_3d_gdeco)
    ls2fs_3d_redist = new_redist(ls_3d_gdeco, fs_subset_gdeco)

    fs2ls_zm_redist = new_redist(fs_zmdeco, ls_zmdeco)
    ls2fs_zm_redist = new_redist(ls_zmdeco, fs_zmdeco)

    IF (aggregate_fs_to_ls) THEN
      DO ir = 1, 6
        redists(ir) = fs2ls_3d_redist
      ENDDO
      redists(7) = fs2ls_zm_redist
      num_redists = 7
      fs2ls_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    IF (aggregate_ls_to_fs) THEN
      DO ir = 1, 9
        redists(ir) = ls2fs_3d_redist
      ENDDO
      DO ir = 10, 12
        redists(ir) = ls2fs_zm_redist
      ENDDO
      num_redists = 12
      ls2fs_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    ! generate ls <--> sp redistributions:
    ALLOCATE(ls2sp_3d_redist(nlev:nlev+1), sp2ls_3d_redist(nlev:nlev+1))
    DO n = nlev, nlev+1
      ls2sp_3d_redist(n) = new_redist(ls_sdeco(n), sp_sdeco(n))
      sp2ls_3d_redist(n) = new_redist(sp_sdeco(n), ls_sdeco(n))
    ENDDO
    ls2sp_m0_redist = new_redist(ls_m0deco, sp_m0deco)
    sp2ls_m0_redist = new_redist(sp_m0deco, ls_m0deco)

    IF (aggregate_ls_to_sp) THEN
      redists(1) = ls2sp_3d_redist(nlev)
      redists(2) = ls2sp_3d_redist(nlev)
      redists(3) = ls2sp_3d_redist(nlev+1)
      redists(4) = ls2sp_m0_redist
      num_redists = 4
      ls2sp_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

    IF (aggregate_sp_to_ls) THEN
      redists(1) = sp2ls_3d_redist(nlev)
      redists(2) = sp2ls_3d_redist(nlev)
      redists(3) = sp2ls_3d_redist(nlev+1)
      redists(4) = sp2ls_m0_redist
      num_redists = 4
      sp2ls_all_redist = xt_redist_collection_new(redists, num_redists, redist_cs, ab_comm)
    ENDIF

#ifdef YAXT_SELFTEST
    CALL selftest
#endif

  CONTAINS

    TYPE(xt_redist) FUNCTION new_redist(src, dst)
      TYPE(yaxt_deco), INTENT(in) :: src, dst
      TYPE(xt_xmap) :: x
      x = xt_xmap_all2all_new(src%idxlist, dst%idxlist, ab_comm)
      new_redist = xt_redist_p2p_off_new(x, src%offset, dst%offset, p_real_dp)
      CALL xt_xmap_delete(x)
    END FUNCTION new_redist

#endif /* HAVE_YAXT */
  END SUBROUTINE generate_yaxt_redist

  SUBROUTINE yaxt_tr_gp_to_ffsl(gp3d1, gp3d2, gp3d3, gp3d4, gp3d5, gp2d1, gp2d2, &
       &                        ffsl3d1, ffsl3d2, ffsl3d3, ffsl3d4, ffsl3d5, ffsl2d1, ffsl2d2)
    REAL(dp), TARGET, INTENT(in)  :: gp3d1(:,:,:)   ! decomposed with 3d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(in)  :: gp3d2(:,:,:)   ! ""
    REAL(dp), TARGET, INTENT(in)  :: gp3d3(:,:,:)   ! ""
    REAL(dp), TARGET, INTENT(in)  :: gp3d4(:,:,:)   ! ""
    REAL(dp), TARGET, INTENT(in)  :: gp3d5(:,:,:)   ! ""
    REAL(dp), TARGET, INTENT(in)  :: gp2d1(:,:)     ! decomposed with 2d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(in)  :: gp2d2(:,:)     ! ""
    REAL(dp), TARGET, INTENT(out) :: ffsl3d1(:,:,:) ! decomposed with 3d ffsl deco.
    REAL(dp), TARGET, INTENT(out) :: ffsl3d2(:,:,:) ! ""
    REAL(dp), TARGET, INTENT(out) :: ffsl3d3(:,:,:) ! ""
    REAL(dp), TARGET, INTENT(out) :: ffsl3d4(:,:,:) ! ""
    REAL(dp), TARGET, INTENT(out) :: ffsl3d5(:,:,:) ! ""
    REAL(dp), TARGET, INTENT(out) :: ffsl2d1(:,:)   ! decomposed with 2d ffsl gridpoint deco.
    REAL(dp), TARGET, INTENT(out) :: ffsl2d2(:,:)   ! ""
#ifdef HAVE_YAXT
    TYPE(c_ptr) :: gp_c_ptr(7), ffsl_c_ptr(7)

#ifdef CHECK_MEM_STRIDES
    INTEGER :: mem_strides_3d(3), mem_strides_2d(2)
    mem_strides_3d = (/ 1, nproma, nproma*nlev /)
    CALL check_mem_strides(gp3d1, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp3d2, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp3d3, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp3d4, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp3d5, mem_strides_3d, 8, __LINE__)
    mem_strides_2d = (/ 1, nproma /)
    CALL check_mem_strides(gp2d1, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(gp2d2, mem_strides_2d, 8, __LINE__)
    mem_strides_3d = (/ 1, nlon, nlon*ffsl_nlat /)
    CALL check_mem_strides(ffsl3d1, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(ffsl3d2, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(ffsl3d3, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(ffsl3d4, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(ffsl3d5, mem_strides_3d, 8, __LINE__)
    mem_strides_2d = (/ 1, nlon /)
    CALL check_mem_strides(ffsl2d1, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(ffsl2d2, mem_strides_2d, 8, __LINE__)
#endif

    IF (aggregate_gp_to_ffsl) THEN
      gp_c_ptr = (/ C_LOC(gp3d1(1,1,1)), C_LOC(gp3d2(1,1,1)), C_LOC(gp3d3(1,1,1)), &
           &        C_LOC(gp3d4(1,1,1)), C_LOC(gp3d5(1,1,1)), &
           &        C_LOC(gp2d1(1,1)), C_LOC(gp2d2(1,1)) /)

      ffsl_c_ptr = (/ C_LOC(ffsl3d1(1,1,1)), C_LOC(ffsl3d2(1,1,1)), C_LOC(ffsl3d3(1,1,1)), &
           &          C_LOC(ffsl3d4(1,1,1)), C_LOC(ffsl3d5(1,1,1)), &
           &          C_LOC(ffsl2d1(1,1)), C_LOC(ffsl2d2(1,1)) /)
      CALL xt_redist_s_exchange(gp2ffsl_all_redist, gp_c_ptr, ffsl_c_ptr)
    ELSE
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d1(1,1,1)), C_LOC(ffsl3d1(1,1,1)) )
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d2(1,1,1)), C_LOC(ffsl3d2(1,1,1)) )
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d3(1,1,1)), C_LOC(ffsl3d3(1,1,1)) )
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d4(1,1,1)), C_LOC(ffsl3d4(1,1,1)) )
      CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(gp3d5(1,1,1)), C_LOC(ffsl3d5(1,1,1)) )
      CALL xt_redist_s_exchange1(gp2ffsl_2d_redist, C_LOC(gp2d1(1,1)),   C_LOC(ffsl2d1(1,1))   )
      CALL xt_redist_s_exchange1(gp2ffsl_2d_redist, C_LOC(gp2d2(1,1)),   C_LOC(ffsl2d2(1,1))   )
    ENDIF
#endif
  END SUBROUTINE yaxt_tr_gp_to_ffsl

  SUBROUTINE yaxt_tr_tracer_gp_to_ffsl(xtm1, fb)
    REAL(dp), TARGET, INTENT(in)  :: xtm1(:,:,:,:)   ! gp tracer data
    REAL(dp), TARGET, INTENT(inout) :: fb(:,:,:,:)       ! ffsl tracer data
#ifdef HAVE_YAXT
    INTEGER :: iit, it
    LOGICAL :: match

    IF (active_trnum == 0) RETURN

#ifdef CHECK_MEM_STRIDES
    CALL check_mem_strides(xtm1, (/ 1, nproma, nproma*nlev, nproma*nlev*tracdef_ntrac /), 8, __LINE__)
    CALL check_mem_strides(fb, (/ 1, nlon, nlon*ffsl_nlat, nlon*ffsl_nlat*ffsl_nlev /), 8, __LINE__)
#endif

    ! check if our redist object still matches the current active tracer list of the model
    match = .TRUE.
    check: DO

      IF (active_trnum /= COUNT(tracdef_trlist%ti(1:tracdef_trlist%ntrac)%ntran /= 0) ) THEN
        match = .FALSE.
        EXIT check
      ENDIF

      DO iit = 1, active_trnum
        it = active_trpos(iit)
        IF (tracdef_trlist%ti(it)%ntran == 0) THEN
          match = .FALSE.
          EXIT check
        ENDIF
      ENDDO

      EXIT check
    ENDDO check

    IF (match) THEN
      CALL xt_redist_s_exchange1(gp2ffsl_4d_xtm1_redist, C_LOC(xtm1(1,1,1,1)),C_LOC(fb(1,1,1,1)))
    ELSE
      CALL fallback
    ENDIF

  CONTAINS

    SUBROUTINE fallback
      LOGICAL, SAVE :: first_call = .TRUE.
      INTEGER :: ib
      REAL(dp), ALLOCATABLE, TARGET :: tmp(:,:,:)

      ! fallback solution: work on copies of 3d parts

      IF (first_call) THEN
        WRITE(0,*) 'WARNING: (yaxt_tr_tracer_gp_to_ffsl) invalid initial active tracerlist; use fallback solution'
        first_call = .FALSE.
      ENDIF

      ALLOCATE(tmp(nproma, nlev, ngpblks))
      ib = advection_jps
      DO it = 1, tracdef_trlist%ntrac
        IF (tracdef_trlist%ti(it)%ntran /= 0) THEN
          tmp = xtm1(:,:,it,:)
          ib = ib + 1
          CALL xt_redist_s_exchange1(gp2ffsl_3d_redist, C_LOC(tmp(1,1,1)), C_LOC(fb(1,1,1,ib)) )
        ENDIF
      ENDDO

    END SUBROUTINE fallback
#endif
  END SUBROUTINE yaxt_tr_tracer_gp_to_ffsl

  SUBROUTINE yaxt_tr_ffsl_to_gp(gp4d, gp2d, gp3d, ffsl4d, ffsl2d, ffsl3d)
    REAL(dp), TARGET, INTENT(out) :: gp4d(:,:,:,:)   ! decomposed with 4d standard (gridpoint&tracer) deco.
    REAL(dp), TARGET, INTENT(out) :: gp2d(:,:)       ! decomposed with 2d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(out) :: gp3d(:,:,:)     ! decomposed with 3d standard gridpoint deco.
    REAL(dp), TARGET, INTENT(in)  :: ffsl4d(:,:,:,:) ! decomposed with 4d ffsl deco.
    REAL(dp), TARGET, INTENT(in)  :: ffsl2d(:,:)     ! decomposed with 2d ffsl deco.
    REAL(dp), TARGET, INTENT(in)  :: ffsl3d(:,:,:)   ! decomposed with 3d ffsl deco.
#ifdef HAVE_YAXT
    TYPE(c_ptr) :: gp_c_ptr(3), ffsl_c_ptr(3)

#ifdef CHECK_MEM_STRIDES
    CALL check_mem_strides(gp4d, (/ 1, nproma, nproma*nlev, nproma*nlev*pcnst /), 8, __LINE__)
    CALL check_mem_strides(gp3d, (/ 1, nproma, nproma*nlev /), 8, __LINE__)
    CALL check_mem_strides(gp2d, (/ 1, nproma /), 8, __LINE__)
    CALL check_mem_strides(ffsl4d, (/ 1, nlon, nlon*ffsl_nlat, nlon*ffsl_nlat*ffsl_nlev /), 8, __LINE__)
    CALL check_mem_strides(ffsl3d, (/ 1, nlon, nlon*ffsl_nlat /), 8, __LINE__)
    CALL check_mem_strides(ffsl2d, (/ 1, nlon /), 8, __LINE__)
#endif

    IF (SIZE(ffsl4d,4) /= SIZE(gp4d,3)) THEN
      WRITE(0,*) 'SIZE(gp4d,3), SIZE(ffsl4d,4)=',SIZE(gp4d,3), SIZE(ffsl4d,4)
      CALL die('yaxt_tr_ffsl_to_gp: tracer dimensions do not match', __LINE__)
    ENDIF

    IF (aggregate_ffsl_to_gp) THEN
      ffsl_c_ptr = (/ C_LOC(ffsl4d(1,1,1,1)), C_LOC(ffsl2d(1,1)), C_LOC(ffsl3d(1,1,1)) /)
      gp_c_ptr = (/  C_LOC(gp4d(1,1,1,1)), C_LOC(gp2d(1,1)), C_LOC(gp3d(1,1,1)) /)
      CALL xt_redist_s_exchange(ffsl2gp_all_redist, ffsl_c_ptr, gp_c_ptr)
    ELSE
      CALL xt_redist_s_exchange1(ffsl2gp_4d_fb_redist, C_LOC(ffsl4d(1,1,1,1)), C_LOC(gp4d(1,1,1,1)))
      CALL xt_redist_s_exchange1(ffsl2gp_2d_redist, C_LOC(ffsl2d(1,1)), C_LOC(gp2d(1,1)))
      CALL xt_redist_s_exchange1(ffsl2gp_3d_redist, C_LOC(ffsl3d(1,1,1)), C_LOC(gp3d(1,1,1)))
    ENDIF
    IF (ldc%npromz /= nproma) THEN
      CALL zero_fractional_block_4d(gp4d)
      CALL zero_fractional_block_3d(gp3d)
      CALL zero_fractional_block_2d(gp2d)
    ENDIF
#endif
  END SUBROUTINE yaxt_tr_ffsl_to_gp


  SUBROUTINE yaxt_tr_gp_to_fs(gp1, gp2, gp3, gp4, gp5, gp6, gp7, &
       &                      sf3, zm1, zm2, zm3, fs, fs0)
    !
    !   grid point space  -> Fourier space
    !
    REAL(dp), TARGET, INTENT(in)    :: gp1(:,:,:)   ! gridpoint space 3d
    REAL(dp), TARGET, INTENT(in)    :: gp2(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp3(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp4(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp5(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp6(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: gp7(:,:,:)   !
    REAL(dp), TARGET, INTENT(in)    :: sf3(:,:)     ! gridpoint space 2d
    REAL(dp), TARGET, INTENT(inout) :: zm1(:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(inout) :: zm2(:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(inout) :: zm3(:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(inout) :: fs (:,:,:,:) ! Fourier space
    REAL(dp), TARGET, INTENT(inout) :: fs0(:,:,:)   ! zonal mean, Four.
#ifdef HAVE_YAXT
    TYPE(c_ptr) :: gp_c_ptr(7)
    TYPE(c_ptr) :: zm_c_ptr(3)
    TYPE(c_ptr) :: all_src_c_ptr(11), all_dst_c_ptr(11)
    INTEGER :: m

#ifdef CHECK_MEM_STRIDES
    INTEGER :: mem_strides_3d(3), mem_strides_2d(2)
    mem_strides_3d = (/ 1, nproma, nproma*nlev /)
    CALL check_mem_strides(gp1, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp2, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp3, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp4, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp5, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp6, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp7, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(sf3, (/ 1, nproma /), 8, __LINE__)
    mem_strides_2d = (/ 1, nlev /)
    CALL check_mem_strides(zm1, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(zm2, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(zm3, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(fs0, (/ 1, nflev, nflev*nflat /), 8, __LINE__)
    CALL check_mem_strides(fs, (/ 1, nlon+2, (nlon+2)*nflevp1, (nlon+2)*nflevp1*nflat /), 8, __LINE__)
#endif

    gp_c_ptr = (/ C_LOC(gp1(1,1,1)), C_LOC(gp2(1,1,1)), C_LOC(gp3(1,1,1)), &
         &        C_LOC(gp4(1,1,1)), C_LOC(gp5(1,1,1)), C_LOC(gp6(1,1,1)), &
         &        C_LOC(gp7(1,1,1))  /)

    zm_c_ptr = (/ C_LOC(zm1(1,1)), C_LOC(zm2(1,1)), C_LOC(zm3(1,1)) /)

    ! set some values to zero as in original echam implementation:
    CALL zero_fs_4d(fs, 1, 7)

    IF (aggregate_gp_to_fs) THEN
      DO m = 1,7
        all_src_c_ptr(m) =  gp_c_ptr(m)
        all_dst_c_ptr(m) =  C_LOC(fs(1,1,1,m))
      ENDDO
      all_src_c_ptr(8) =  C_LOC(sf3(1,1))
      all_dst_c_ptr(8) =  C_LOC(fs(1,nflevp1,1,3))
      DO m = 1,3
        all_src_c_ptr(8+m) =  zm_c_ptr(m)
        all_dst_c_ptr(8+m) =  C_LOC(fs0(1,1,m))
      ENDDO
      CALL xt_redist_s_exchange(gp2fs_all_redist, all_src_c_ptr, all_dst_c_ptr)
    ELSE
      DO m = 1, 7
        CALL xt_redist_s_exchange1(gp2fs_3d_redist, gp_c_ptr(m), C_LOC(fs(1,1,1,m)));
      ENDDO
      CALL xt_redist_s_exchange1(gp2fs_2d_redist, C_LOC(sf3(1,1)), C_LOC(fs(1,nflevp1,1,3)))
      DO m = 1, 3
        CALL xt_redist_s_exchange1(gp2fs_zm_redist, zm_c_ptr(m), C_LOC(fs0(1,1,m)))
      ENDDO
    ENDIF
#endif
  END SUBROUTINE yaxt_tr_gp_to_fs

  SUBROUTINE yaxt_tr_fs_to_gp(gp1, gp2, gp3, gp4, gp5, gp6, gp7, gp8, gp9, &
       &                      sf1, sf2, sf3, zm1, zm2, zm3, fs, fs0)
    !
    !   grid point space  -> Fourier space
    !
    REAL(dp), TARGET, INTENT(inout)  :: gp1    (:,:,:)   ! gridpoint space 3d
    REAL(dp), TARGET, INTENT(inout)  :: gp2    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp3    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp4    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp5    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp6    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp7    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp8    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: gp9    (:,:,:)   !
    REAL(dp), TARGET, INTENT(inout)  :: sf1    (:,:)     ! gridpoint space 2d
    REAL(dp), TARGET, INTENT(inout)  :: sf2    (:,:)     !
    REAL(dp), TARGET, INTENT(inout)  :: sf3    (:,:)     !
    REAL(dp), TARGET, INTENT(inout)  :: zm1    (:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(inout)  :: zm2    (:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(inout)  :: zm3    (:,:)     ! zonal mean
    REAL(dp), TARGET, INTENT(in)     :: fs     (:,:,:,:) ! Fourier space
    REAL(dp), TARGET, INTENT(in)     :: fs0    (:,:,:)   ! zonal mean, Four.
#ifdef HAVE_YAXT
    !
    TYPE(c_ptr) :: gp_c_ptr(9)
    TYPE(c_ptr) :: sf_c_ptr(3)
    TYPE(c_ptr) :: zm_c_ptr(3)
    TYPE(c_ptr) :: all_src_c_ptr(15), all_dst_c_ptr(15)
    INTEGER :: m, npromz

#ifdef CHECK_MEM_STRIDES
    INTEGER :: mem_strides_3d(3), mem_strides_2d(2)
    mem_strides_3d = (/ 1, nproma, nproma*nlev /)
    CALL check_mem_strides(gp1, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp2, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp3, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp4, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp5, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp6, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp7, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp8, mem_strides_3d, 8, __LINE__)
    CALL check_mem_strides(gp9, mem_strides_3d, 8, __LINE__)
    mem_strides_2d = (/ 1, nproma /)
    CALL check_mem_strides(sf1, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(sf2, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(sf3, mem_strides_2d, 8, __LINE__)
    mem_strides_2d = (/ 1, nlev /)
    CALL check_mem_strides(zm1, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(zm2, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(zm3, mem_strides_2d, 8, __LINE__)
    CALL check_mem_strides(fs, (/ 1, nlon+2, (nlon+2)*nflevp1, (nlon+2)*nflevp1*nflat /), 8, __LINE__)
    CALL check_mem_strides(fs0, (/ 1, nflev, nflev*nflat /), 8, __LINE__)
#endif

    gp_c_ptr = (/ C_LOC(gp1(1,1,1)), C_LOC(gp2(1,1,1)), C_LOC(gp3(1,1,1)), &
         &        C_LOC(gp4(1,1,1)), C_LOC(gp5(1,1,1)), C_LOC(gp6(1,1,1)), &
         &        C_LOC(gp7(1,1,1)), C_LOC(gp8(1,1,1)), C_LOC(gp9(1,1,1))  /)

    sf_c_ptr = (/ C_LOC(sf1(1,1)), C_LOC(sf2(1,1)), C_LOC(sf3(1,1)) /)
    zm_c_ptr = (/ C_LOC(zm1(1,1)), C_LOC(zm2(1,1)), C_LOC(zm3(1,1)) /)

    npromz = ldc%npromz

    IF (aggregate_fs_to_gp) THEN
      DO m = 1,9
        all_src_c_ptr(m) =  C_LOC(fs(1,1,1,m))
        all_dst_c_ptr(m) =  gp_c_ptr(m)
      ENDDO
      DO m = 1,3
        all_src_c_ptr(9+m) =  C_LOC(fs(1,nflevp1,1,m))
        all_dst_c_ptr(9+m) =  sf_c_ptr(m)
      ENDDO
      DO m = 1,3
        all_src_c_ptr(12+m) =  C_LOC(fs0(1,1,m))
        all_dst_c_ptr(12+m) =  zm_c_ptr(m)
      ENDDO
      CALL xt_redist_s_exchange(fs2gp_all_redist, all_src_c_ptr, all_dst_c_ptr)
    ELSE
      DO m = 1, 9
        CALL xt_redist_s_exchange1(fs2gp_3d_redist, C_LOC(fs(1,1,1,m)), gp_c_ptr(m))
      ENDDO
      DO m = 1, 3
        CALL xt_redist_s_exchange1(fs2gp_2d_redist, C_LOC(fs(1,nflevp1,1,m)), sf_c_ptr(m) )
      ENDDO
      DO m = 1, 3
        CALL xt_redist_s_exchange1(fs2gp_zm_redist, C_LOC(fs0(1,1,m)), zm_c_ptr(m))
      ENDDO
    ENDIF

    IF (npromz /= nproma) THEN
      CALL zero_fractional_block_3d(gp1)
      CALL zero_fractional_block_3d(gp2)
      CALL zero_fractional_block_3d(gp3)
      CALL zero_fractional_block_3d(gp4)
      CALL zero_fractional_block_3d(gp5)
      CALL zero_fractional_block_3d(gp6)
      CALL zero_fractional_block_3d(gp7)
      CALL zero_fractional_block_3d(gp8)
      CALL zero_fractional_block_3d(gp9)
      CALL zero_fractional_block_2d(sf1)
      CALL zero_fractional_block_2d(sf2)
      CALL zero_fractional_block_2d(sf3)
    ENDIF

#endif
  END SUBROUTINE yaxt_tr_fs_to_gp

  SUBROUTINE yaxt_tr_fs_to_ls(fs, ls, fs0, ls0)
    REAL(dp), TARGET, INTENT(in)  :: fs(:,:,:,:)
    REAL(dp), TARGET, INTENT(out) :: ls(:,:,:,:)
    REAL(dp), TARGET, INTENT(in)  :: fs0(:,:,:)
    REAL(dp), TARGET, INTENT(out) :: ls0(:,:,:)
#ifdef HAVE_YAXT
    INTEGER, PARAMETER :: n4_expected = 6
    INTEGER :: m, n4
    TYPE(c_ptr) :: fs_c_ptr(7)
    TYPE(c_ptr) :: ls_c_ptr(7)

#ifdef CHECK_MEM_STRIDES
    CALL check_mem_strides(fs, (/ 1, nlon+2, (nlon+2)*nflevp1, (nlon+2)*nflevp1*nflat /), 8, __LINE__)
    CALL check_mem_strides(fs0, (/ 1, nflev, nflev*nflat /), 8, __LINE__)
    CALL check_mem_strides(ls, (/ 1, 2*ldc%nlm, 2*ldc%nlm*nflevp1, 2*ldc%nlm*nflevp1*nlat /), 8, __LINE__)
    CALL check_mem_strides(ls0, (/ 1, nflev, nflev*nlat /), 8, __LINE__)
#endif

    n4 = MIN(SIZE(fs,4), SIZE(ls,4))
    IF (aggregate_fs_to_ls .AND. n4 == n4_expected) THEN
      DO m = 1, 6
        fs_c_ptr(m) =  C_LOC(fs(1,1,1,m))
        ls_c_ptr(m) =  C_LOC(ls(1,1,1,m))
      ENDDO
      m = 1
      fs_c_ptr(7) =  C_LOC(fs0(1,1,m))
      ls_c_ptr(7) =  C_LOC(ls0(1,1,m))
      CALL xt_redist_s_exchange(fs2ls_all_redist, fs_c_ptr, ls_c_ptr)
    ELSE
      DO m = 1, n4
        CALL xt_redist_s_exchange1(fs2ls_3d_redist, C_LOC(fs(1,1,1,m)), C_LOC(ls(1,1,1,m)) )
      ENDDO
      m=1
      CALL xt_redist_s_exchange1(fs2ls_zm_redist, C_LOC(fs0(1,1,m)), C_LOC(ls0(1,1,m)) )
    ENDIF
#endif
  END SUBROUTINE yaxt_tr_fs_to_ls

  SUBROUTINE yaxt_tr_ls_to_fs(fs, ls, fs0, ls0)
    REAL(dp), TARGET, INTENT(out)  :: fs(:,:,:,:)
    REAL(dp), TARGET, INTENT(in)   :: ls(:,:,:,:)
    REAL(dp), TARGET, INTENT(out)  :: fs0(:,:,:)
    REAL(dp), TARGET, INTENT(in)   :: ls0(:,:,:)
#ifdef HAVE_YAXT
    TYPE(c_ptr) :: fs_c_ptr(12)
    TYPE(c_ptr) :: ls_c_ptr(12)
    INTEGER, PARAMETER :: n1_expected = 9
    INTEGER, PARAMETER :: n2_expected = 3
    INTEGER :: m, n1, n2

#ifdef CHECK_MEM_STRIDES
    CALL check_mem_strides(fs, (/ 1, nlon+2, (nlon+2)*nflevp1, (nlon+2)*nflevp1*nflat /), 8, __LINE__)
    CALL check_mem_strides(fs0, (/ 1, nflev, nflev*nflat /), 8, __LINE__)
    CALL check_mem_strides(ls, (/ 1, 2*ldc%nlm, 2*ldc%nlm*nflevp1, 2*ldc%nlm*nflevp1*nlat /), 8, __LINE__)
    CALL check_mem_strides(ls0, (/ 1, nflev, nflev*nlat /), 8, __LINE__)
#endif

    n1 = MIN(SIZE(fs,4), SIZE(ls,4))
    n2 = MIN(SIZE(fs0,3), SIZE(ls0,3))

    IF (aggregate_ls_to_fs .AND. n1==n1_expected .AND. n2==n2_expected) THEN
      DO m = 1, n1
        fs_c_ptr(m) =  C_LOC(fs(1,1,1,m))
        ls_c_ptr(m) =  C_LOC(ls(1,1,1,m))
      ENDDO
      DO m = 1, n2
        fs_c_ptr(n1+m) =  C_LOC(fs0(1,1,m))
        ls_c_ptr(n1+m) =  C_LOC(ls0(1,1,m))
      ENDDO
      CALL xt_redist_s_exchange(ls2fs_all_redist, ls_c_ptr, fs_c_ptr)
      DO m = 1, n1
        CALL zero_fs_3d(fs(:,:,:,m))
      ENDDO
    ELSE
      DO m = 1, n1
        CALL xt_redist_s_exchange1(ls2fs_3d_redist, C_LOC(ls(1,1,1,m)), C_LOC(fs(1,1,1,m)) )
        CALL zero_fs_3d(fs(:,:,:,m))
      ENDDO
      DO m = 1, n2
        CALL xt_redist_s_exchange1(ls2fs_zm_redist, C_LOC(ls0(1,1,m)), C_LOC(fs0(1,1,m)) )
      ENDDO
    ENDIF

  CONTAINS

    ! 3d fs post processing
    SUBROUTINE zero_fs_3d(f)
      REAL(dp), INTENT(inout) :: f(:,:,:)
      INTEGER :: i, jj, kk

      DO jj = 1, nflat
        DO kk = 1, nflevp1
          DO i = 2*(ldc%nm+1)+1, nlon+2
            f(i,kk,jj) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE zero_fs_3d
#endif
  END SUBROUTINE yaxt_tr_ls_to_fs

  SUBROUTINE yaxt_tr_ls_to_sp(ls1, sp1, ls2, sp2, ls3, sp3, ls0, sp0)
    REAL(dp), TARGET, INTENT(in)    :: ls1(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(inout) :: sp1(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(in)    :: ls2(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(inout) :: sp2(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(in)    :: ls3(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(inout) :: sp3(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(in)    :: ls0(:,:)   ! Legendre (m=0 only)
    REAL(dp), TARGET, INTENT(inout) :: sp0(:,:)   ! spectral (m=0 only)

#ifdef HAVE_YAXT
    TYPE(c_ptr) :: ls_c_ptr(4), sp_c_ptr(4), ls0_c_ptr, sp0_c_ptr
    ! use to prevent operations on null pointer
    REAL(dp), TARGET, SAVE :: dummy(1)

#ifdef CHECK_MEM_STRIDES
    INTEGER :: ls_nk(nlev:nlev+1)
    ls_nk(nlev)   = MIN(ldc%lleve,nlev)   - ldc%llevs + 1
    ls_nk(nlev+1) = MIN(ldc%lleve,nlev+1) - ldc%llevs + 1
    CALL check_mem_strides(ls1, (/ 1, ls_nk(nlev), 2*ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(ls3, (/ 1, ls_nk(nlev), 2*ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(ls2, (/ 1, ls_nk(nlev+1), 2*ls_nk(nlev+1) /), 8, __LINE__)
    CALL check_mem_strides(sp1, (/ 1, nlev, 2*nlev /), 8, __LINE__)
    CALL check_mem_strides(sp3, (/ 1, nlev, 2*nlev /), 8, __LINE__)
    CALL check_mem_strides(sp2, (/ 1, nlev+1, 2*(nlev+1)/), 8, __LINE__)
    CALL check_mem_strides(ls0, (/ 1, ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(sp0, (/ 1, nlev /), 8, __LINE__)
#endif

    IF (SIZE(ls0)>0) THEN
      ls0_c_ptr = C_LOC(ls0(1,1))
    ELSE
      ls0_c_ptr = C_LOC(dummy)
    ENDIF
    IF (SIZE(sp0)>0) THEN
      sp0_c_ptr = C_LOC(sp0(1,1))
    ELSE
      sp0_c_ptr = C_LOC(dummy)
    ENDIF

    IF (aggregate_ls_to_sp) THEN
      ls_c_ptr(1) = C_LOC(ls1(1,1,1))
      sp_c_ptr(1) = C_LOC(sp1(1,1,1))

      ls_c_ptr(2) = C_LOC(ls3(1,1,1))
      sp_c_ptr(2) = C_LOC(sp3(1,1,1))

      ls_c_ptr(3) = C_LOC(ls2(1,1,1))
      sp_c_ptr(3) = C_LOC(sp2(1,1,1))

      ls_c_ptr(4) = ls0_c_ptr
      sp_c_ptr(4) = sp0_c_ptr
      CALL xt_redist_s_exchange(ls2sp_all_redist, ls_c_ptr, sp_c_ptr)
    ELSE
      CALL xt_redist_s_exchange1(ls2sp_3d_redist(nlev), C_LOC(ls1(1,1,1)), C_LOC(sp1(1,1,1)) )
      CALL xt_redist_s_exchange1(ls2sp_3d_redist(nlev), C_LOC(ls3(1,1,1)), C_LOC(sp3(1,1,1)) )
      CALL xt_redist_s_exchange1(ls2sp_3d_redist(nlev+1), C_LOC(ls2(1,1,1)), C_LOC(sp2(1,1,1)) )
      CALL xt_redist_s_exchange1(ls2sp_m0_redist, ls0_c_ptr, sp0_c_ptr)
    ENDIF
#endif
  END SUBROUTINE yaxt_tr_ls_to_sp

  SUBROUTINE yaxt_tr_sp_to_ls(ls1, sp1, ls2, sp2, ls3, sp3, ls0, sp0)
    REAL(dp), TARGET, INTENT(inout) :: ls1(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(in)    :: sp1(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(inout) :: ls2(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(in)    :: sp2(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(inout) :: ls3(:,:,:) ! Legendre space
    REAL(dp), TARGET, INTENT(in)    :: sp3(:,:,:) ! spectral space
    REAL(dp), TARGET, INTENT(inout) :: ls0(:,:)   ! Legendre (m=0 only)
    REAL(dp), TARGET, INTENT(in)    :: sp0(:,:)   ! spectral (m=0 only)
#ifdef HAVE_YAXT
    TYPE(c_ptr) :: ls_c_ptr(4), sp_c_ptr(4), ls0_c_ptr, sp0_c_ptr
    ! use to prevent operations on null pointer
    REAL(dp), TARGET, SAVE :: dummy(1)

#ifdef CHECK_MEM_STRIDES
    INTEGER :: ls_nk(nlev:nlev+1)
    ls_nk(nlev)   = MIN(ldc%lleve,nlev)   - ldc%llevs + 1
    ls_nk(nlev+1) = MIN(ldc%lleve,nlev+1) - ldc%llevs + 1
    CALL check_mem_strides(ls1, (/ 1, ls_nk(nlev), 2*ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(ls3, (/ 1, ls_nk(nlev), 2*ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(ls2, (/ 1, ls_nk(nlev+1), 2*ls_nk(nlev+1) /), 8, __LINE__)
    CALL check_mem_strides(sp1, (/ 1, nlev, 2*nlev /), 8, __LINE__)
    CALL check_mem_strides(sp3, (/ 1, nlev, 2*nlev /), 8, __LINE__)
    CALL check_mem_strides(sp2, (/ 1, nlev+1, 2*(nlev+1)/), 8, __LINE__)
    CALL check_mem_strides(ls0, (/ 1, ls_nk(nlev) /), 8, __LINE__)
    CALL check_mem_strides(sp0, (/ 1, nlev /), 8, __LINE__)
#endif

    IF (SIZE(ls0)>0) THEN
      ls0_c_ptr = C_LOC(ls0(1,1))
    ELSE
      ls0_c_ptr = C_LOC(dummy)
    ENDIF
    IF (SIZE(sp0)>0) THEN
      sp0_c_ptr = C_LOC(sp0(1,1))
    ELSE
      sp0_c_ptr = C_LOC(dummy)
    ENDIF

    IF (aggregate_sp_to_ls) THEN
      ls_c_ptr(1) = C_LOC(ls1(1,1,1))
      sp_c_ptr(1) = C_LOC(sp1(1,1,1))

      ls_c_ptr(2) = C_LOC(ls3(1,1,1))
      sp_c_ptr(2) = C_LOC(sp3(1,1,1))

      ls_c_ptr(3) = C_LOC(ls2(1,1,1))
      sp_c_ptr(3) = C_LOC(sp2(1,1,1))

      ls_c_ptr(4) = ls0_c_ptr
      sp_c_ptr(4) = sp0_c_ptr
      CALL xt_redist_s_exchange(sp2ls_all_redist, sp_c_ptr, ls_c_ptr)
    ELSE
      CALL xt_redist_s_exchange1(sp2ls_3d_redist(nlev), C_LOC(sp1(1,1,1)), C_LOC(ls1(1,1,1)) )
      CALL xt_redist_s_exchange1(sp2ls_3d_redist(nlev), C_LOC(sp3(1,1,1)), C_LOC(ls3(1,1,1)) )
      CALL xt_redist_s_exchange1(sp2ls_3d_redist(nlev+1), C_LOC(sp2(1,1,1)), C_LOC(ls2(1,1,1)) )
      CALL xt_redist_s_exchange1(sp2ls_m0_redist, sp0_c_ptr, ls0_c_ptr)
    ENDIF
#endif
  END SUBROUTINE yaxt_tr_sp_to_ls

#ifdef HAVE_YAXT
  !
  ! start of yaxt-only section until end of module
  !

  ! zero selected 4d fs data
  SUBROUTINE zero_fs_4d(f,m1,m2)
    REAL(dp), INTENT(inout) :: f(:,:,:,:)
    INTEGER, INTENT(in) :: m1, m2
    INTEGER :: i, jj, kk, m, f_shape(4)

    f_shape = SHAPE(f)
    DO m = m1, m2
      DO jj = 1, nflat
        DO kk = 1, nflevp1
          DO i = nlon+1, nlon+2
            f(i,kk,jj,m) = 0.0_dp
          ENDDO
        ENDDO
        DO kk = nflev+1,nflevp1
          DO i = 1, nlon+2
            f(i,kk,jj,m) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE zero_fs_4d

  ! zero tail of fractional nproma block for 4d arrays
  SUBROUTINE zero_fractional_block_4d(g)
    REAL(dp), INTENT(inout) :: g(:,:,:,:)
    INTEGER :: i, k, it, nt, npromz
    nt = SIZE(g,3)
    npromz = ldc%npromz
    DO k=1, nlev
      DO it = 1, nt
        DO i = npromz+1, nproma
          g(i,k,it,ngpblks) = 0.0_dp
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE zero_fractional_block_4d

  SUBROUTINE zero_fractional_block_3d(g)
    REAL(dp), INTENT(inout) :: g(:,:,:)
    INTEGER :: i, k, npromz
    npromz = ldc%npromz
    DO k=1, nlev
      DO i = npromz+1, nproma
        g(i,k,ngpblks) = 0.0_dp
      ENDDO
    ENDDO
  END SUBROUTINE zero_fractional_block_3d

  SUBROUTINE zero_fractional_block_2d(g)
    REAL(dp), INTENT(inout) :: g(:,:)
    INTEGER :: i, npromz
    npromz = ldc%npromz
    DO i = npromz+1, nproma
      g(i,ngpblks) = 0.0_dp
    ENDDO
  END SUBROUTINE zero_fractional_block_2d

  SUBROUTINE setup_gp_gdeco(idx_description_type)
    INTEGER, INTENT(IN) :: idx_description_type
    INTEGER :: nglon, nglat

    nglon = ldc%nglon
    nglat = ldc%nglat

    CALL setup_gp

  CONTAINS

    !
    ! Define global index list and offsets for Gaussian grid
    !
    SUBROUTINE setup_gp
      INTEGER :: nidx, nstr, nk, i
      INTEGER :: trpos(pcnst)

      ALLOCATE(gp_gdeco(SIZE(gp_nlevs)))
      ! the next line is only needed for some buggy compilers
      DO nk = 1, SIZE(gp_gdeco)
        gp_gdeco(nk)%idxlist = xt_idxlist_c2f(c_null_ptr)
      ENDDO

      DO nk = 1, SIZE(gp_nlevs)
        ! generate index list
        IF (gp_nlevs(nk)) THEN
          IF (.NOT. xt_is_null(gp_gdeco(nk)%idxlist)) &
               & CALL finish('mo_echam_yaxt: setup_gp', 'Internal error.')
          nidx = nglon * nk * nglat
          SELECT CASE (idx_description_type)
            ! index list as vector
            CASE (idx_descr_vector)
              gp_gdeco(nk)%idxlist = new_gp_idxvec(nidx, nk)
            ! index list by stripes
            CASE (idx_descr_stripes)
              nstr = nk * nglat
              gp_gdeco(nk)%idxlist = new_gp_idxstripes(nstr, nk)
            CASE DEFAULT
              WRITE (message_text, '(a,i0,a)') &
              'Unsupported global index desciption type: ', &
              idx_description_type, ' (Gaussain grid)'
              CALL finish('setup_yaxt_decomposition', message_text)
          END SELECT

          ! generate offsets
          ALLOCATE(gp_gdeco(nk)%offset(nidx))
          CALL set_gp_offset(nidx, nk, gp_gdeco(nk)%offset)
        ENDIF
      ENDDO

      ! gp deco for tracer:
      nidx = nglon * nlev * nglat * pcnst
      gp_4d_fb_gdeco%idxlist = new_gp_4d_idxvec(nidx, nlev, pcnst)
      ALLOCATE(gp_4d_fb_gdeco%offset(nidx))
      DO i = 1, pcnst
        trpos(i) = i
      ENDDO

      CALL set_gp_4d_offset(nidx, nlev, trpos, pcnst, pcnst, gp_4d_fb_gdeco%offset)

      IF (active_trnum /= 0) THEN
        nidx = nglon * nlev * nglat * active_trnum
        gp_4d_xtm1_gdeco%idxlist = new_gp_4d_idxvec(nidx, nlev, active_trnum)
        ALLOCATE(gp_4d_xtm1_gdeco%offset(nidx))
        CALL set_gp_4d_offset(nidx, nlev, active_trpos, active_trnum, trdim_size, gp_4d_xtm1_gdeco%offset)
      ENDIF

    END SUBROUTINE setup_gp

    !
    ! Gaussian grid
    ! use index vector to describe decomposition
    TYPE(xt_idxlist) FUNCTION new_gp_idxvec(nidx, nlv)
      INTEGER, INTENT(IN) :: nidx, nlv
      INTEGER :: idx(nidx)
      INTEGER :: i, j, k, p, r

      p = 0
      DO k = 1, nlv
        DO r = 1, 2
          DO j = glats(r), glate(r)
            DO i = glons(r), glone(r)
              p = p + 1
              idx(p) = i + nlon * ( (j-1) + nlat * (k-1) ) - 1
            ENDDO !i
          ENDDO !j
        ENDDO !r
      ENDDO !k
      new_gp_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_gp_idxvec

    !
    ! Gaussian grid
    ! use index stripes to describe decomposition
    TYPE(xt_idxlist) FUNCTION new_gp_idxstripes(nstr, nlv)
      INTEGER, INTENT(IN) :: nstr, nlv
      TYPE(xt_stripe) :: s(nstr)
      INTEGER :: p, k, r, j

      s(:)%stride = 1
      p = 0 ! strides index
      DO k = 1, nlv
        DO r = 1, 2
          DO j = glats(r), glate(r)
            p = p + 1
            s(p)%start = glons(r) + nlon * ( (j-1) + nlat * (k-1) ) - 1
            s(p)%nstrides = nglon
          ENDDO !j
        ENDDO !r
      ENDDO !k
      new_gp_idxstripes = xt_idxstripes_new(s, nstr)
    END FUNCTION new_gp_idxstripes

    !
    ! Gaussian grid
    ! set offsets
    SUBROUTINE set_gp_offset(nidx, nlv, offset)
      INTEGER, INTENT(IN)  :: nidx, nlv
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: i, j, k, p, r, ib, ia

      p = 0
      DO k = 1, nlv
        ib = 1
        ia = 0
        DO r = 1, 2
          DO j = glats(r), glate(r)
            DO i = glons(r), glone(r)
              p = p + 1
              ia = ia +1
              IF (ia > nproma) THEN
                ib = ib + 1
                ia = 1
              ENDIF
              offset(p) = ia + nproma * (k - 1) + nproma * nlv * (ib - 1) - 1
            ENDDO !i
          ENDDO !j
        ENDDO !r
      ENDDO !k
    END SUBROUTINE set_gp_offset

    !
    ! Gaussian grid with tracer dimension
    ! use index vector to describe decomposition
    TYPE(xt_idxlist) FUNCTION new_gp_4d_idxvec(nidx, nlv, nt)
      INTEGER, INTENT(IN) :: nidx, nlv, nt
      INTEGER :: idx(nidx)
      INTEGER :: i, j, k, it, p, r

      p = 0
      DO it = 1, nt
        DO k = 1, nlv
          DO r = 1, 2
            DO j = glats(r), glate(r)
              DO i = glons(r), glone(r)
                p = p + 1
                idx(p) = (i-1) + nlon * ( (j-1) + nlat * ( (k-1) + nlv * (it-1) ) )
              ENDDO !i
            ENDDO !j
          ENDDO !r
        ENDDO !k
      ENDDO
      IF (p /= nidx) CALL die('bad case', __LINE__)
      new_gp_4d_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_gp_4d_idxvec

    !
    ! Gaussian grid with tracer dimension
    ! set offsets
    SUBROUTINE set_gp_4d_offset(nidx, nlv, trpos, trnum, trsize, offset)
      INTEGER, INTENT(IN)  :: nidx, nlv, trnum, trsize
      INTEGER, INTENT(IN)  :: trpos(:)
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: i, j, k, iit,it, p, r, ib, ia
      INTEGER :: w_ia, w_ib, w_k, w_it

      ! offsets weights for different loops:
      w_ia = 1
      w_k  = w_ia * nproma
      w_it = w_k * nlv
      w_ib = w_it * trsize

      p = 0
      DO iit = 1, trnum
        it = trpos(iit)
        DO k = 1, nlv
          ib = 1
          ia = 0
          DO r = 1, 2
            DO j = glats(r), glate(r)
              DO i = glons(r), glone(r)
                p = p + 1
                ia = ia +1
                IF (ia > nproma) THEN
                  ib = ib + 1
                  ia = 1
                ENDIF
                offset(p) = (ia-1)*w_ia + (ib-1)*w_ib + (k-1)*w_k + (it-1)*w_it
              ENDDO !i
            ENDDO !j
          ENDDO !r
        ENDDO !k
      ENDDO !it
    END SUBROUTINE set_gp_4d_offset

  END SUBROUTINE setup_gp_gdeco

  ! decomposition & offsets of zonal means
  SUBROUTINE setup_gp_zmdeco
    INTEGER :: nglat

    nglat = ldc%nglat
    CALL setup_gp_zm

  CONTAINS

    SUBROUTINE setup_gp_zm
      INTEGER :: nidx
      ! zonal means for all levels:
      nidx = nlev*nglat
      gp_zmdeco%idxlist = new_gp_zm_idxvec(nidx)
      ALLOCATE(gp_zmdeco%offset(nidx))
      CALL set_gp_zm_offset(nidx, gp_zmdeco%offset)

      ! zonal means for fs levels:
      nidx = nflev*nglat
      gp_sel_zmdeco%idxlist = new_gp_sel_zm_idxvec(nidx)
      ALLOCATE(gp_sel_zmdeco%offset(nidx))
      CALL set_gp_sel_zm_offset(nidx, gp_sel_zmdeco%offset)

    END SUBROUTINE setup_gp_zm

    ! indices for grid point zonal means for all levels
    TYPE(xt_idxlist) FUNCTION new_gp_zm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)
      INTEGER :: k, j, r, p

      p = 0
      DO r = 1, 2
        DO j = glats(r), glate(r)
          DO k = 1, nlev
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = k + nlev * (j-1)  - 1
          ENDDO
        ENDDO
      ENDDO
      new_gp_zm_idxvec = xt_idxvec_new(idx, nidx)

    END FUNCTION new_gp_zm_idxvec

    ! offsets for grid point zonal means for all levels
    SUBROUTINE set_gp_zm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: jj, k, p

      ! data shape: (nlev, nglat)
      ! see: variables ul, u0, du0 in module mo_scan_buffer

      p = 0
      DO jj = 1, nglat
        DO k = 1, nlev
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          offset(p) = k - 1 + nlev*(jj-1)
        ENDDO
      ENDDO

    END SUBROUTINE set_gp_zm_offset

    ! indices for grid point zonal means within fs levels
    TYPE(xt_idxlist) FUNCTION new_gp_sel_zm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)
      INTEGER :: k, j, r, p, confined_fleve

      ! indices belong to zonal index space
      confined_fleve = MIN(fleve, nlev)

      p = 0
      DO r = 1, 2
        DO j = glats(r), glate(r)
          DO k = flevs, confined_fleve
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = k + nlev * (j-1)  - 1
          ENDDO
        ENDDO
      ENDDO
      new_gp_sel_zm_idxvec = xt_idxvec_new(idx, nidx)

    END FUNCTION new_gp_sel_zm_idxvec

    ! offsets for grid point zonal means within fs levels
    SUBROUTINE set_gp_sel_zm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: jj, k, p, confined_fleve

      ! data shape: (nlev, nglat)
      ! see: variables ul, u0, du0 in module mo_scan_buffer
      confined_fleve = MIN(fleve, nlev)

      p = 0
      DO jj = 1, nglat
        DO k = flevs, confined_fleve
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          offset(p) = k - 1 + nlev*(jj-1)
        ENDDO
      ENDDO

    END SUBROUTINE set_gp_sel_zm_offset

  END SUBROUTINE setup_gp_zmdeco


  SUBROUTINE setup_fs_2d_gdeco
    INTEGER :: num_indices

    ! the surface case is only used within the extra fs layer:
    IF (nflevp1 > nflev) THEN
      num_indices = nflat*nlon
      fs_2d_gdeco%idxlist = new_fs_2d_idxvec(num_indices)
      ALLOCATE(fs_2d_gdeco%offset(num_indices))
      CALL set_fs_2d_offset(num_indices, fs_2d_gdeco%offset)
    ELSE
      ! this pe has no extra layer
      fs_2d_gdeco%idxlist = xt_idxempty_new()
      ! we still need a valid address for offsets
      ALLOCATE(fs_2d_gdeco%offset(1))
      fs_2d_gdeco%offset(1) = 0
    ENDIF

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_fs_2d_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)

      INTEGER :: i, j, r, p

      p = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          DO i = 1, nlon
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = i + nlon * (j-1)  - 1
          ENDDO
        ENDDO
      ENDDO

      new_fs_2d_idxvec = xt_idxvec_new( idx, nidx)
    END FUNCTION new_fs_2d_idxvec

    SUBROUTINE set_fs_2d_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: r, i, j, jj, p

      ! fs 2d-data: (1:nlon+2, nflevp1, 1:nflat)

      p = 0
      jj = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          jj = jj + 1
          DO i = 1, nlon
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            offset(p) = i + (nlon+2) * nflevp1 * (jj-1) - 1
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE set_fs_2d_offset

  END SUBROUTINE setup_fs_2d_gdeco

  SUBROUTINE setup_fs_3d_gdeco

    CALL setup_fs_3d

  CONTAINS

    SUBROUTINE setup_fs_3d
      INTEGER :: nidx, confined_fleve, global_kmax, nk(nlev:nlev+1)

      ! std fs deco:
      nk = (/ nflev, nflevp1 /)
      ALLOCATE(fs_gdeco(nlev:nlev+1))
      DO global_kmax = nlev, nlev+1
        confined_fleve = MIN(fleve, global_kmax)
        nidx = nk(global_kmax)*nflat*nlon
        fs_gdeco(global_kmax)%idxlist = new_fs_3d_idxvec(nidx, nlon, confined_fleve)
        ALLOCATE(fs_gdeco(global_kmax)%offset(nidx))
        CALL set_fs_3d_offset(nidx, fs_gdeco(global_kmax)%offset, nlon, confined_fleve)
      ENDDO

      ! special subset deco for transpositions between fs and ls:
      nidx = nflevp1*nflat*2*(ldc%nm+1)
      fs_subset_gdeco%idxlist = new_fs_3d_idxvec(nidx, 2*(ldc%nm+1), fleve)
      ALLOCATE(fs_subset_gdeco%offset(nidx))
      CALL set_fs_3d_offset(nidx, fs_subset_gdeco%offset, 2*(ldc%nm+1), fleve)

    END SUBROUTINE setup_fs_3d

    ! fs 3d indices
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_fs_3d_idxvec(nidx,active_nlon,active_fleve)
      INTEGER, INTENT(IN) :: nidx, active_nlon, active_fleve
      INTEGER :: idx(nidx)

      INTEGER :: i, j, k, r, p

      ! assumptions:
      !
      ! in fs-space echam distributes nlev+1 levels in block form,
      ! these blocks are described by (flevs:fleve), we have nflevp1 = fleve-flevs+1
      ! and nflev= min(fleve,nlev)-flevs, empty blocks have nlevp1 <= 0,
      ! the transpositions from/to gp-deco only touch levels <= nlev,
      ! the level nlev+1 carries extra data (outside the normal grid point space)

      IF (fleve-flevs+1 /= nflevp1) CALL die('bad case',__LINE__)

      p = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          DO k = flevs, active_fleve
            DO i = 1, active_nlon
              p = p + 1
              IF (p > nidx) CALL die('bad case', __LINE__)
              idx(p) = i + nlon * ( (j-1) + nlat * (k-1) ) - 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      new_fs_3d_idxvec = xt_idxvec_new( idx, nidx)
    END FUNCTION new_fs_3d_idxvec

    SUBROUTINE set_fs_3d_offset(nidx, offset, active_nlon, active_fleve)
      INTEGER, INTENT(IN) :: nidx, active_nlon, active_fleve
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: r, i, j, k, jc, kc, p

      ! assumptions:
      !
      ! fs data layout:
      ! shape(fftz) =  (dc% nlon+2, dc% nflevp1, dc% nflat, nvar)

      p = 0
      jc = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          jc = jc + 1
          DO k = flevs, active_fleve
            kc = k - flevs + 1
            DO i = 1, active_nlon
              p = p + 1
              IF (p > nidx) CALL die('bad case', __LINE__)
              offset(p) = i-1 + (nlon+2) * ( (kc-1)  + nflevp1 * (jc-1))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE set_fs_3d_offset

  END SUBROUTINE setup_fs_3d_gdeco

  ! fs zonal means:
  SUBROUTINE setup_fs_zmdeco
    INTEGER :: num_indices

    num_indices = nflev*nflat
    fs_zmdeco%idxlist = new_fs_zm_idxvec(num_indices)
    ALLOCATE(fs_zmdeco%offset(num_indices))
    CALL set_fs_zm_offset(num_indices, fs_zmdeco%offset)

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_fs_zm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)

      INTEGER :: k, j, r, p, confined_fleve

      confined_fleve = MIN(fleve, nlev)

      p = 0
      DO r = 1, 2
        DO j = flats(r), flate(r)
          DO k = flevs, confined_fleve
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = k + nlev * (j-1)  - 1
          ENDDO
        ENDDO
      ENDDO
      new_fs_zm_idxvec = xt_idxvec_new( idx, nidx)
    END FUNCTION new_fs_zm_idxvec

    SUBROUTINE set_fs_zm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: k, jj, p

      ! offsets for zonal means

      ! data shape: (nflev, nflat)
      ! see: first two dimensions of fbm0 in mo_buffer_fft

      p = 0
      DO jj = 1, nflat
        DO k = 1, nflev
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          offset(p) = p -1
        ENDDO
      ENDDO
    END SUBROUTINE set_fs_zm_offset

  END SUBROUTINE setup_fs_zmdeco

  !
  ! Define global index list and offsets for spectral coefficients
  !
  SUBROUTINE setup_sp_sdeco(idx_description_type)
    INTEGER, INTENT(IN) :: idx_description_type
    INTEGER :: nm, nsm, snsp
    INTEGER, ALLOCATABLE :: sm(:), nmp(:), snn0(:), snnp(:)

    ! Spectral decomposition data
    nm  = ldc%nm   ! spectral truncation
    nsm = ldc%nsm  ! number of m wave numbers on PE
    snsp= ldc%snsp ! number of spectral (complex) coefficients on PE
    ALLOCATE(sm(nsm), snn0(nsm), snnp(nsm), nmp(nm+2))
    sm  = ldc%sm
    nmp = ldc%nmp  ! displacement of the first point of
                   ! m-columns with respect to the first point
                   ! of the first m-column
    snn0 = ldc%snn0
    snnp = ldc%snnp

    CALL setup_sp

  CONTAINS

    SUBROUTINE setup_sp
      INTEGER :: nlev_list(3), nidx, nstr, i

      ! array with numbers of vertical levels / tiles
      nlev_list = (/1, nlev, nlev+1/)

      ALLOCATE(sp_sdeco(MAXVAL(nlev_list)))
      ! the next line is only needed for some buggy compilers
      DO i = 1, SIZE(sp_sdeco)
        sp_sdeco(i)%idxlist = xt_idxlist_c2f(c_null_ptr)
      ENDDO

      DO i = 1, SIZE(nlev_list)
        ! generate index list
        IF (xt_is_null(sp_sdeco(nlev_list(i))%idxlist)) THEN
          nidx = 2 * snsp * nlev_list(i)
          SELECT CASE (idx_description_type)
            ! index list as vector
          CASE (idx_descr_vector)
            sp_sdeco(nlev_list(i))%idxlist = new_sp_idxvec(nidx, nlev_list(i))
            ! index list by stripes
          CASE (idx_descr_stripes)
            nstr = nsm * nlev_list(i)
            sp_sdeco(nlev_list(i))%idxlist = new_sp_idxstripes(nstr, nlev_list(i))
          CASE DEFAULT
            WRITE (message_text, '(a,i0,a)') &
                 'Unsupported global index desciption type: ', &
                 idx_description_type, ' (spectral data)'
            CALL finish('setup_yaxt_decomposition', message_text)
          END SELECT

          ! generate offsets
          ALLOCATE(sp_sdeco(nlev_list(i))%offset(nidx))
          CALL set_sp_offset(nidx, nlev_list(i), sp_sdeco(nlev_list(i))%offset)
        ENDIF
      ENDDO

    END SUBROUTINE setup_sp

    !
    ! spectral data
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_sp_idxvec(nidx, nlv)
      INTEGER, INTENT(IN) :: nidx, nlv
      INTEGER :: idx(nidx)
      INTEGER :: nsp, mp1, spgls, spgle
      INTEGER :: p, k, im, n, c

      ! global 3d spectral index space: sp3_index(1:2, 0:nsp-1, 1:nlev)
      ! sp3_index(c,n,k) := c-1 + 2*n + 2*nsp*(k-1)

      nsp = (nm + 1) * ( nm + 2) / 2

      p = 0
      DO k = 1, nlv
        DO im = 1, nsm
          mp1  = sm(im) + 1
          spgls = nmp(mp1) + snn0(im)
          spgle = spgls + snnp(im) - 1
          DO n = spgls, spgle
            DO c = 1, 2
              p = p + 1
              idx(p) = c + 2 * n + 2 * nsp * (k - 1) - 1
            ENDDO
          ENDDO
        END DO
      END DO
      new_sp_idxvec = xt_idxvec_new(idx, nidx)

    END FUNCTION new_sp_idxvec

    !
    ! spectral data
    ! describe decomposition by index stripes
    TYPE(xt_idxlist) FUNCTION new_sp_idxstripes(nstr, nlv)
      INTEGER, INTENT(IN) :: nstr, nlv
      TYPE(xt_stripe) :: s(nstr)
      INTEGER :: nsp, p, k, im, mp1, spgls

      nsp = (nm + 1) * ( nm + 2) / 2
      s(:)%stride = 1
      p = 0 ! strides index
      DO k = 1, nlv
        DO im = 1, nsm
          mp1  = sm(im) + 1
          spgls = nmp(mp1) + snn0(im)
          p = p + 1
          s(p)%start = 2 * spgls + 2 * nsp * (k - 1)
          s(p)%nstrides = 2 * snnp(im)
        END DO
      END DO
      new_sp_idxstripes = xt_idxstripes_new(s, nstr)

    END FUNCTION new_sp_idxstripes
    !
    ! spectral data
    ! set offsets
    SUBROUTINE set_sp_offset(nidx, nlv, offset)
      INTEGER, INTENT(IN)  :: nidx, nlv
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: p, k, n, c
      INTEGER :: im, mp1, spgls, spgle

      !
      ! local echam 3d spectral offsets: sp3_offset(1:nlv,1:2,1:snsp)
      ! sp3_offset(k,c,n) = k-1 + nlv*(c-1) + nlv*2*(n-1)
      !

      p = 0
      DO k = 1, nlv
        DO n = 1, snsp
          DO c = 1, 2
            p = p + 1
            offset(p) = k + nlv * (c - 1) + nlv * 2 * (n - 1) - 1
          ENDDO !k
        ENDDO !n
      ENDDO !c

    END SUBROUTINE set_sp_offset

  END SUBROUTINE setup_sp_sdeco

  !
  ! decomposition & offsets for spectral coeff. with m==0
  SUBROUTINE setup_sp_m0deco()
    INTEGER :: num_indices

    num_indices = nlev*ldc%nsnm0
    sp_m0deco%idxlist = new_sp_m0_idxvec(num_indices, nlev)
    ALLOCATE(sp_m0deco%offset(num_indices))
    CALL set_sp_m0_offset(num_indices, nlev, sp_m0deco%offset)

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_sp_m0_idxvec(nidx, global_kmax)
      INTEGER, INTENT(IN) :: nidx, global_kmax
      INTEGER :: idx(nidx)
      INTEGER :: n, p, k, nlev, nsnm0, snn0

      ! see scatter_sp0 for details

      nlev   = global_kmax
      nsnm0  = ldc%nsnm0
      snn0   = ldc%snn0(1)

      p = 0
      DO n = snn0+1,snn0+nsnm0
        DO k = 1, nlev
          p = p + 1
          idx(p) = k + nlev * (n-1)
        ENDDO
      END DO

      new_sp_m0_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_sp_m0_idxvec

    SUBROUTINE set_sp_m0_offset(nidx, global_kmax, offset)
      INTEGER, INTENT(IN)  :: nidx, global_kmax
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: p, k, n

      p = 0
      DO n = 1, ldc%nsnm0
        DO k = 1, global_kmax
          p = p + 1
          offset(p) = p-1
        ENDDO
      ENDDO

    END SUBROUTINE set_sp_m0_offset

  END SUBROUTINE setup_sp_m0deco

  !
  ! Define global index list and offsets for FFSL 2D-grid
  !
  SUBROUTINE setup_ffsl_2d_gdeco(idx_description_type)
    INTEGER, INTENT(IN) :: idx_description_type

    CALL setup_ffsl_2d

  CONTAINS

    SUBROUTINE setup_ffsl_2d
      INTEGER :: nidx, nstr

      nidx = nlon * ffsl_nlat
      ! generate index list
      SELECT CASE(idx_description_type)
        CASE (idx_descr_vector)
          ffsl_2d_gdeco%idxlist = new_ffsl_2d_idxvec(nidx)
        CASE (idx_descr_stripes)
          nstr = ffsl_nlat
          ffsl_2d_gdeco%idxlist = new_ffsl_2d_idxstripes(nstr)
        CASE DEFAULT
          CALL finish('setup_yaxt_decomposition', &
              'Index list description type is not valid')
      END SELECT

      ! generate offsets
      ALLOCATE(ffsl_2d_gdeco%offset(nidx))
      CALL set_ffsl_2d_offset(nidx, ffsl_2d_gdeco%offset)
    END SUBROUTINE setup_ffsl_2d

    !
    ! FFSL 2D data
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_ffsl_2d_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)
      INTEGER :: i, j, p

      p = 0
      DO j = ffsl_gp_lat1, ffsl_gp_lat2
        DO i = 1, nlon
          p = p + 1
          idx(p) =i + nlon * (j-1) - 1
        ENDDO
      ENDDO
      new_ffsl_2d_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_ffsl_2d_idxvec

    !
    ! FFSL 2D data
    ! describe decomposition by index stripes
    TYPE(xt_idxlist) FUNCTION new_ffsl_2d_idxstripes(nstr)
      INTEGER, INTENT(IN) :: nstr
      TYPE(xt_stripe)     :: s(nstr)
      INTEGER :: j, p

      s(:)%stride = 1
      s(:)%nstrides = nlon
      p = 0
      DO j = ffsl_gp_lat1, ffsl_gp_lat2
        p = p + 1
        s(p)%start = 1 + nlon * (j-1) - 1
      ENDDO
      new_ffsl_2d_idxstripes = xt_idxstripes_new(s, nstr)
    END FUNCTION new_ffsl_2d_idxstripes

    !
    ! FFSL 2D data
    ! set offsets
    SUBROUTINE set_ffsl_2d_offset(nidx, offset)
      INTEGER, INTENT(IN)  :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: coords_off(nlon, ffsl_nlat)
      INTEGER :: i, j, ic, jc, p

      ! data offsets - here we take care of the reversed latitute iteration in ffsl space
      p = 0
      DO jc = ffsl_nlat, 1, -1 ! change in j-orientation
        DO ic = 1, nlon
          coords_off(ic, jc) = p
          p = p + 1
        ENDDO
      ENDDO

      p = 0
      DO j = ffsl_gp_lat1, ffsl_gp_lat2
        jc = j - ffsl_gp_lat1 + 1
        DO i = 1, nlon
          ic = i
          p = p + 1
          offset(p) = coords_off(ic, jc)
        ENDDO
      ENDDO
    END SUBROUTINE set_ffsl_2d_offset

  END SUBROUTINE setup_ffsl_2d_gdeco

  !
  ! Define global index list and offsets for FFSL 3D-grid
  !
  SUBROUTINE setup_ffsl_3d_gdeco(idx_description_type)
    INTEGER, INTENT(IN) :: idx_description_type
    INTEGER, ALLOCATABLE :: ffsl_kstack(:)

    CALL setup_ffsl_3d

  CONTAINS

    SUBROUTINE setup_ffsl_3d
      INTEGER :: nidx, nstr
      INTEGER :: kstack(ldc%nlev)

      CALL get_ffsl_kstack(kstack, ffsl_nlev)
      ALLOCATE(ffsl_kstack(ffsl_nlev))
      ffsl_kstack = kstack(1:ffsl_nlev)
      nidx = nlon * ffsl_nlat * ffsl_nlev

      ! generate index list
      SELECT CASE(idx_description_type)
        CASE (idx_descr_vector)
          ffsl_3d_gdeco%idxlist = new_ffsl_3d_idxvec(nidx)
        CASE (idx_descr_stripes)
          nstr = ffsl_nlev * ffsl_nlat
          ffsl_3d_gdeco%idxlist = new_ffsl_3d_idxstripes(nstr)
        CASE DEFAULT
          CALL finish('setup_yaxt_decomposition', &
              'Index list description type is not valid')
      END SELECT

      ! generate offsets
      ALLOCATE(ffsl_3d_gdeco%offset(nidx))
      CALL set_ffsl_3d_offset(nidx, ffsl_3d_gdeco%offset)
    END SUBROUTINE setup_ffsl_3d

    !
    ! FFSL 3D data
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_ffsl_3d_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)

      INTEGER :: i, j, k, kk, p

      p = 0
      DO kk = 1, ffsl_nlev
        k = ffsl_kstack(kk)
        DO j = ffsl_gp_lat1, ffsl_gp_lat2
          DO i = 1, nlon
            p = p + 1
            idx(p) = i + nlon * ( (j-1) + nlat * (k-1) ) - 1
          ENDDO
        ENDDO
      ENDDO
      new_ffsl_3d_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_ffsl_3d_idxvec
    !
    ! FFSL 3D data
    ! describe decomposition by index stripes
    TYPE(xt_idxlist) FUNCTION new_ffsl_3d_idxstripes(nstr)
      INTEGER, INTENT(IN) :: nstr
      TYPE(xt_stripe) :: s(nstr)

      INTEGER :: j, k, kk, p

      s(:)%stride = 1
      s(:)%nstrides = nlon
      p = 0
      DO kk = 1, ffsl_nlev
        k = ffsl_kstack(kk)
        DO j = ffsl_gp_lat1, ffsl_gp_lat2
          p = p + 1
          s(p)%start = 1 + nlon * ( (j-1) + nlat * (k-1) ) - 1
        ENDDO
      ENDDO
      new_ffsl_3d_idxstripes = xt_idxstripes_new(s, nstr)
    END FUNCTION new_ffsl_3d_idxstripes
    !
    ! FFSL 3D data
    ! set offsets
    SUBROUTINE set_ffsl_3d_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: coords_off(nlon, ffsl_nlat, ffsl_nlev)
      INTEGER :: i, j, kk, ic, jc, kc, p

      ! acces to data offsets via coords
      p = 0
      DO kc = 1, ffsl_nlev
        DO jc = ffsl_nlat, 1, -1 ! change in j-orientation
          DO ic = 1, nlon
            coords_off(ic, jc, kc) = p
            p = p + 1
          ENDDO
        ENDDO
      ENDDO

      ! offsets in index order (must match order in idxlist definition)
      p = 0
      DO kk = 1, ffsl_nlev
        kc = kk
        DO j = ffsl_gp_lat1, ffsl_gp_lat2
          jc = j-ffsl_gp_lat1+1
          DO i = 1, nlon
            ic = i
            p = p + 1
            offset(p) = coords_off(ic, jc, kc)
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE set_ffsl_3d_offset

  END SUBROUTINE setup_ffsl_3d_gdeco

  !
  ! Define global index list and offsets for FFSL 4D-grid
  !
  SUBROUTINE setup_ffsl_4d_gdeco
    INTEGER, ALLOCATABLE :: ffsl_kstack(:)


    CALL setup_ffsl_4d

  CONTAINS

    SUBROUTINE setup_ffsl_4d
      INTEGER :: nidx
      INTEGER :: kstack(nlev)

      CALL get_ffsl_kstack(kstack, ffsl_nlev)
      ALLOCATE(ffsl_kstack(ffsl_nlev))
      ffsl_kstack = kstack(1:ffsl_nlev)

      nidx = nlon * ffsl_nlat * ffsl_nlev * pcnst
      ! generate index list:
      ffsl_4d_fb_gdeco%idxlist = new_ffsl_4d_idxvec(nidx,pcnst)
      ! generate offsets:
      ALLOCATE(ffsl_4d_fb_gdeco%offset(nidx))
      CALL set_ffsl_4d_offset(nidx, 0, pcnst, ffsl_4d_fb_gdeco%offset)

      IF (active_trnum /= 0) THEN
        nidx = nlon * ffsl_nlat * ffsl_nlev * active_trnum
        ! generate index list:
        ffsl_4d_tcfb_gdeco%idxlist = new_ffsl_4d_idxvec(nidx,active_trnum)
        ! generate offsets:
        ALLOCATE(ffsl_4d_tcfb_gdeco%offset(nidx))
        CALL set_ffsl_4d_offset(nidx, advection_jps, active_trnum, ffsl_4d_tcfb_gdeco%offset)
      ENDIF
    END SUBROUTINE setup_ffsl_4d
    !
    ! FFSL 4 data
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_ffsl_4d_idxvec(nidx,nt)
      INTEGER, INTENT(IN) :: nidx, nt
      INTEGER :: idx(nidx)

      INTEGER :: i, j, k, it, kk, p

      p = 0
      DO it = 1, nt
        DO kk = 1, ffsl_nlev
          k = ffsl_kstack(kk)
          DO j = ffsl_gp_lat1, ffsl_gp_lat2
            DO i = 1, nlon
              p = p + 1
              idx(p) = i-1 + nlon * ( (j-1) + nlat * ( (k-1) + nlev * (it-1) ) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF (p /= nidx) CALL die('bad case', __LINE__)
      new_ffsl_4d_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_ffsl_4d_idxvec
    !
    ! FFSL 4D data
    ! set offsets
    SUBROUTINE set_ffsl_4d_offset(nidx, trshift, nt, offset)
      INTEGER, INTENT(IN) :: nidx, trshift, nt
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: coords_off(nlon, ffsl_nlat, ffsl_nlev, nt)
      INTEGER :: i, j, it, kk, ic, jc, kc, tc, p, offshift

      offshift = trshift * ffsl_nlev * ffsl_nlat * nlon

      ! access to data offsets via coords
      p = 0
      DO tc = 1, nt
        DO kc = 1, ffsl_nlev
          DO jc = ffsl_nlat, 1, -1 ! change in j-orientation
            DO ic = 1, nlon
              coords_off(ic, jc, kc, tc) = p + offshift
              p = p + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      ! offsets in index order (must match order in idxlist definition)
      p = 0
      DO it = 1, nt
        tc = it
        DO kk = 1, ffsl_nlev
          kc = kk
          DO j = ffsl_gp_lat1, ffsl_gp_lat2
            jc = j-ffsl_gp_lat1+1
            DO i = 1, nlon
              ic = i
              p = p + 1
              offset(p) = coords_off(ic, jc, kc, tc)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF (p /= nidx) CALL die('bad case', __LINE__)
    END SUBROUTINE set_ffsl_4d_offset

  END SUBROUTINE setup_ffsl_4d_gdeco

  !
  ! Define global index list and offsets for LS 3D
  !
  SUBROUTINE setup_ls_3d_gdeco()
    INTEGER :: num_indices

    num_indices = nlat * nflevp1 * 2 * ldc%nlm

    ! generate index list
    ls_3d_gdeco%idxlist = new_ls_3d_idxvec(num_indices)

    ! generate offsets
    ALLOCATE(ls_3d_gdeco%offset(num_indices))
    CALL set_ls_3d_offset(num_indices, ls_3d_gdeco%offset)

  CONTAINS

    ! Legendre space indices for 3d data,
    ! describe decomposition by index vector
    TYPE(xt_idxlist) FUNCTION new_ls_3d_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)

      INTEGER :: i, j, k, ii, p
      INTEGER, POINTER :: intr(:) ! index array

      intr => ldc%intr
      p = 0
      DO j = 1, nlat
        DO k = flevs, fleve
          DO ii = 1, 2*ldc%nlm
            i = intr(ii)
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            idx(p) = i + nlon * ( (j-1) + nlat * (k-1) ) - 1
            IF (idx(p) <0 .OR. idx(p) > nlon*nlat*(nlev+1)) CALL die("idx out of range",__LINE__)
          ENDDO
        ENDDO
      ENDDO
      new_ls_3d_idxvec = xt_idxvec_new( idx, nidx)

    END FUNCTION new_ls_3d_idxvec

    !
    ! Legendre space 3D data
    ! set offsets
    SUBROUTINE set_ls_3d_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: ii, j, k, p

      ! assumption:
      ! fftl shape = (nlm*2, nflevp1, nlat , nvar)

      p = 0
      DO j = 1, nlat
        DO k = flevs, fleve
          DO ii = 1, 2*ldc%nlm
            p = p + 1
            IF (p > nidx) CALL die('bad case', __LINE__)
            offset(p) = p - 1
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE set_ls_3d_offset

  END SUBROUTINE setup_ls_3d_gdeco

  !
  ! Define global index list and offsets for LS zm
  !
  SUBROUTINE setup_ls_zmdeco()
    INTEGER :: num_indices

    IF (ldc%nlnm0 > 0) THEN
      num_indices = nflev*nlat
      ls_zmdeco%idxlist = new_ls_zm_idxvec(num_indices)
      ALLOCATE(ls_zmdeco%offset(num_indices))
      CALL set_ls_zm_offset(num_indices, ls_zmdeco%offset)
    ELSE
      ls_zmdeco%idxlist = xt_idxempty_new()
      ALLOCATE(ls_zmdeco%offset(1))
      ls_zmdeco%offset(1)=0
    ENDIF

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_ls_zm_idxvec(nidx)
      INTEGER, INTENT(IN) :: nidx
      INTEGER :: idx(nidx)
      INTEGER :: j, k, p, confined_fleve

      confined_fleve = MIN(fleve, nlev)

      p = 0
      DO j = 1, nlat
        DO k = flevs, confined_fleve
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          idx(p) =  k + nlev * (j-1)  - 1
        ENDDO
      ENDDO

      new_ls_zm_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_ls_zm_idxvec

    SUBROUTINE set_ls_zm_offset(nidx, offset)
      INTEGER, INTENT(IN) :: nidx
      INTEGER, INTENT(OUT) :: offset(nidx)
      INTEGER :: j, k, p

      ! offsets for zonal means
      ! data shape: (nflev, nlat)
      ! see: first two dimensions of lbm0 in mo_buffer_fft

      p = 0
      DO j = 1, nlat
        DO k = 1, nflev
          p = p + 1
          IF (p > nidx) CALL die('bad case', __LINE__)
          offset(p) = p -1
        ENDDO
      ENDDO
    END SUBROUTINE set_ls_zm_offset

  END SUBROUTINE setup_ls_zmdeco

  ! legendre decomposition of spectral coeff.:
  SUBROUTINE setup_ls_sdeco()

    CALL setup_ls

  CONTAINS

    SUBROUTINE setup_ls
      INTEGER :: nidx, global_kmax, nk, ke

      ALLOCATE(ls_sdeco(nlev:nlev+1))

      DO global_kmax = nlev, nlev+1

        ke     = MIN(ldc%lleve, global_kmax)
        nk     = ke - ldc%llevs + 1
        nidx = nk * 2 * ldc%lnsp

        ! generate index list
        ls_sdeco(global_kmax)%idxlist = new_ls3_idxvec(nidx, global_kmax)

        ! generate offsets
        ALLOCATE(ls_sdeco(global_kmax)%offset(nidx))
        CALL set_ls3_offset(nidx, ls_sdeco(global_kmax)%offset, global_kmax)

      ENDDO
    END SUBROUTINE setup_ls

    TYPE(xt_idxlist) FUNCTION new_ls3_idxvec(nidx, global_kmax)
      INTEGER, INTENT(IN) :: nidx, global_kmax
      INTEGER :: idx(nidx)

      ! global 3d legendre index space: ls3_index(1:2, 0:nsp-1, 1:nlev)
      ! ls3_index(c,n,k) := c-1 + 2*n + 2*nsp*(k-1)
      ! Note: sp3 and ls3 use the same global index space

      INTEGER :: mp, nsp
      INTEGER :: c, p, k, ke, im, mp1, mpgl, np, n

      nsp = control_nsp
      ke     = MIN(ldc%lleve, global_kmax)

      p = 0
      DO k = ldc%llevs, ke
        mp=0
        DO im=1,ldc%nlm
          mp1  = ldc%lm(im)+1
          np   = ldc%nnp(mp1)
          mpgl = ldc%nmp(mp1)
          DO n = mpgl+1, mpgl+np
            DO c = 1, 2
              p = p + 1
              idx(p)= c-1 + 2 * (n-1) + 2 * nsp * (k-1)
            ENDDO
          ENDDO
          mp = mp + np
        ENDDO
      ENDDO

      new_ls3_idxvec = xt_idxvec_new( idx, nidx)
    END FUNCTION new_ls3_idxvec

    SUBROUTINE set_ls3_offset(nidx, offset, global_kmax)
      INTEGER, INTENT(IN) :: nidx, global_kmax
      INTEGER, INTENT(out) :: offset(nidx)
      INTEGER :: mp1, llevs, lleve, np, ke, nk
      INTEGER :: p, k, im, n, c

      llevs  = ldc%llevs
      lleve  = ldc%lleve
      ke = MIN (lleve,global_kmax)
      nk = ke - llevs + 1

      p = 0
      DO k = 1, nk
        DO n = 0, ldc%lnsp-1
          DO c = 1, 2
            p = p + 1
            offset(p) = k + nk*(c-1) + (n-0)*2*nk - 1
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE set_ls3_offset

  END SUBROUTINE setup_ls_sdeco

  SUBROUTINE setup_ls_m0deco()
    INTEGER :: num_indices, global_kmax, nk, ke

    global_kmax = nlev
    ke     = MIN(ldc%lleve, global_kmax)
    nk     = ke - ldc%llevs + 1
    num_indices = nk * ldc%nlnm0
    ls_m0deco%idxlist = new_ls_m0_idxvec(num_indices)
    ALLOCATE(ls_m0deco%offset(num_indices))
    CALL set_ls_m0_offset(num_indices,ls_m0deco%offset)

  CONTAINS

    TYPE(xt_idxlist) FUNCTION new_ls_m0_idxvec(nidx)
      INTEGER, INTENT(in) :: nidx
      INTEGER :: idx(nidx), p, k, n
      p = 0
      DO n = 1, ldc%nlnm0
        DO k = ldc%llevs,ke
          p = p + 1
          idx(p) = k + global_kmax*(n-1)
        ENDDO
      ENDDO

      new_ls_m0_idxvec = xt_idxvec_new(idx, nidx)
    END FUNCTION new_ls_m0_idxvec

    SUBROUTINE set_ls_m0_offset(nidx, offset)
      INTEGER, INTENT(in) :: nidx
      INTEGER, INTENT(out) :: offset(nidx)
      INTEGER :: p, k, n
      p = 0
      DO n = 1, ldc%nlnm0
        DO k = ldc%llevs,ke
          p = p + 1
          offset(p) = p-1
        ENDDO
      ENDDO
    END SUBROUTINE set_ls_m0_offset

  END SUBROUTINE setup_ls_m0deco


    SUBROUTINE get_ffsl_kstack(kstack, kstack_n)
      INTEGER, INTENT(OUT) :: kstack(:)
      INTEGER, INTENT(OUT) :: kstack_n
      INTEGER :: dest_idx(2), dest_set_b
      INTEGER :: k, ia, my_idx

      my_idx = ldc%pe + 1

      ! echam::mo_transpose send-logic
      kstack_n = 0
      kstack = 0
      DO k = 1, ldc%nlev
        !
        ! for PEs with odd set_a:
        !   Northern Hemisphere sent to same set_a
        !   Southern Hemisphere sent to set_a+1 (unless set_a  ==  nproca)
        ! for PEs with even set_a:
        !   Northern Hemisphere sent to set_a-1
        !   Southern Hemisphere sent to same set_a
        !
        dest_set_b  =  MOD(k-1,nprocb)+1 ! set_b of receiving PE for this k
        IF( MOD(ldc%set_a,2)  ==  1 ) THEN
          dest_idx(1) = ldc%mapmesh(dest_set_b,ldc%set_a) ! N-region goes to my nprocb proc-family
          ia = MIN(ldc%set_a+1,ldc%nproca)
          dest_idx(2) = ldc%mapmesh(dest_set_b,ia) ! S-region goes preferably to northern nprocb proc-family
        ELSE
          dest_idx(1) = ldc%mapmesh(dest_set_b,ldc%set_a-1) ! N-region goes to southern nprocb proc-family
          dest_idx(2) = ldc%mapmesh(dest_set_b,ldc%set_a) ! S-region goes to my nprocb proc-family
        ENDIF

        IF (dest_idx(1) == my_idx .OR. dest_idx(2) == my_idx) THEN
          kstack_n = kstack_n + 1
          kstack(kstack_n) = k
        ENDIF

      ENDDO
    END SUBROUTINE get_ffsl_kstack

  ! set some required sub communicators of the given model communicator
  SUBROUTINE set_comms(model_comm)
    USE mpi
    INTEGER, INTENT(IN) :: model_comm  ! communicator of application without I/O servers
    INTEGER :: ierror, color, key

    IF (debug_parallel < 0) THEN
      ! use own echam communicator to avoid tag collision:
      CALL mpi_comm_dup(model_comm, ab_comm, ierror)
      IF (ierror /= MPI_SUCCESS) CALL die('mpi_comm_dup failed', __LINE__)
    ELSE
      ! If echam runs in debug_parallel mode then the process space size is nproca x nprocb + 1.
      ! In that case we want to split the process space into two spaces of sizes nproca x nprocb and 1 x 1.
      color = ldc%spe
      key = ldc%pe
      CALL mpi_comm_split(model_comm, color, key, ab_comm, ierror)
      IF (ierror /= MPI_SUCCESS) CALL die('mpi_comm_split failed', __LINE__)
    ENDIF

    CALL xt_mpi_comm_mark_exclusive(ab_comm)

    ! nproca x nprocb space:
    CALL mpi_comm_size(ab_comm, ab_size, ierror)
    CALL mpi_comm_rank(ab_comm, ab_rank, ierror)

    IF (ab_size /= nproca * nprocb) CALL die('set_comms: internal error (1)',__LINE__)

  END SUBROUTINE set_comms

#ifdef YAXT_SELFTEST
  SUBROUTINE selftest()
    REAL(dp), ALLOCATABLE :: gp_2d_testdata(:,:), gp_3d_testdata(:,:,:)
    REAL(dp), ALLOCATABLE :: ffsl_2d_testdata(:,:), ffsl_3d_testdata(:,:,:)

    ! test transpositions:
    ALLOCATE( gp_2d_testdata(nproma, ngpblks), ffsl_2d_testdata(nlon, ffsl_nlat) )
    CALL set_gp_2d_testdata(gp_2d_testdata)
    CALL exchange(gp2ffsl_2d_redist, gp_2d_testdata, ffsl_2d_testdata);
    CALL check_gp2ffsl_2d(gp_2d_testdata, ffsl_2d_testdata)

    ALLOCATE( gp_3d_testdata(nproma, nlev, ngpblks), ffsl_3d_testdata(nlon, ffsl_nlat, ffsl_nlev) )
    CALL set_gp_3d_testdata(gp_3d_testdata)
    CALL exchange(gp2ffsl_3d_redist, gp_3d_testdata, ffsl_3d_testdata)
    CALL check_gp2ffsl_3d(gp_3d_testdata, ffsl_3d_testdata)

    IF (ab_rank == 0) WRITE(0,*) '(mo_echam_yaxt) Selftest passed.'

  END SUBROUTINE selftest

  SUBROUTINE check_gp2ffsl_2d(gp, ffsl)
    REAL(dp), INTENT(IN) :: gp(:,:)
    REAL(dp), INTENT(IN) :: ffsl(:,:)
    REAL(dp) :: ref_gp(SIZE(gp,1), SIZE(gp,2))
    REAL(dp) :: ref_ffsl(SIZE(ffsl,1), SIZE(ffsl,2))

    ref_gp = gp
    CALL tr_gp_ffsl(ldc , 1 ,ref_gp, ref_ffsl)
    IF (ANY(ffsl /= ref_ffsl)) THEN
      CALL FINISH('mo_echam_yaxt::check_gp2ffsl_2d', &
                  'selftest for gp2ffsl 2d failed')
    ELSE
      CALL message('mo_echam_yaxt::check_gp2ffsl_2d', &
                   'selftest for gp2ffsl 2d passed')
    ENDIF
  END SUBROUTINE check_gp2ffsl_2d

  SUBROUTINE check_gp2ffsl_3d(gp, ffsl)
    REAL(dp), INTENT(IN) :: gp(:,:,:)
    REAL(dp), INTENT(IN) :: ffsl(:,:,:)
    REAL(dp) :: ref_gp(SIZE(gp,1), SIZE(gp,2), SIZE(gp,3))
    REAL(dp) :: ref_ffsl(SIZE(ffsl,1), SIZE(ffsl,2), SIZE(ffsl,3))

    ref_gp = gp
    CALL tr_gp_ffsl(ldc , 1 ,ref_gp, ref_ffsl)
    IF (ANY(ffsl /= ref_ffsl)) THEN
      CALL FINISH('mo_echam_yaxt::check_gp2ffsl_3d', &
                  'selftest for gp2ffsl 3d failed')
    ELSE
      CALL message('mo_echam_yaxt::check_gp2ffsl_3d', &
                   'selftest for gp2ffsl 3d passed')
    ENDIF
  END SUBROUTINE check_gp2ffsl_3d
#endif /* YAXT_SELFTEST */

  SUBROUTINE exchange(redist, f, g)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(dp), TARGET, INTENT(in) :: f(*)
    REAL(dp), TARGET, INTENT(out) :: g(*)
    CALL xt_redist_s_exchange1(redist, C_LOC(f), C_LOC(g));
  END SUBROUTINE exchange

  SUBROUTINE set_gp_2d_testdata(f)
    REAL(dp), INTENT(OUT) :: f(:,:) !f(ia,ib)
    INTEGER :: i, j, ia, ib, r
    INTEGER :: ival

    ib = 1
    ia = 0
    DO r = 1, 2
      DO j = glats(r), glate(r)
        DO i = glons(r), glone(r)
          ia = ia +1
          IF (ia > nproma) THEN
            ib = ib + 1
            ia = 1
          ENDIF
          ival = i + (j-1)*1000
          f(ia,ib) = REAL(ival, dp)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE set_gp_2d_testdata

  SUBROUTINE set_gp_3d_testdata(f)
    REAL(dp), INTENT(OUT) :: f(:,:,:) !f(ia,k,ib)
    INTEGER :: i, j, k, ia, ib, r
    INTEGER :: ival

    DO k = 1, nlev
      ib = 1
      ia = 0
      DO r = 1, 2
        DO j = glats(r), glate(r)
          DO i = glons(r), glone(r)
            ia = ia +1
            IF (ia > nproma) THEN
              ib = ib + 1
              ia = 1
            ENDIF
            ival = i + ( (j-1) + (k-1)*1000)*1000
            f(ia,k,ib) = REAL(ival, dp)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE set_gp_3d_testdata

  ! module specific termination subroutine,
  ! simplifies code
  SUBROUTINE die(reason, line)
    USE mo_exception,  ONLY: finish
    USE yaxt,          ONLY: xt_finalize, xt_abort
    CHARACTER(len=*), INTENT(in) :: reason
    INTEGER, INTENT(in) :: line
    CHARACTER(len=*), PARAMETER :: file = __FILE__

    WRITE(message_text,*) reason,', file ',file,', line',line
    CALL xt_finalize()
    CALL xt_abort(TRIM(reason), __FILE__, line)
    CALL finish('mo_echam_yaxt', message_text)
    WRITE(0,*) 'model_abort: fallback stop'
    STOP
  END SUBROUTINE die

#ifdef CHECK_MEM_STRIDES
  SUBROUTINE check_mem_strides_4d(f, memstride, element_size, line)
    USE mpi
    CHARACTER(len=*), PARAMETER :: prefix = 'check_mem_strides_4d: '
    REAL(dp), INTENT(in) :: f(:,:,:,:)
    INTEGER, INTENT(in) :: memstride(:) ! vector of memory strides in units of element_size
    INTEGER, INTENT(in) :: element_size, line
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a, b
    INTEGER :: ierror

    IF (SIZE(memstride) /= 4) CALL die(prefix//'bad usage',__LINE__)

    IF (SIZE(f)<2) RETURN
    CALL mpi_get_address(f(1,1,1,1), a, ierror)

    IF (SIZE(f,1) >= 2) THEN
      CALL mpi_get_address(f(2,1,1,1), b, ierror)
      IF ( b-a /= INT(memstride(1)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 1', line)
    ENDIF

    IF (SIZE(f,2) >= 2) THEN
      CALL mpi_get_address(f(1,2,1,1), b, ierror)
      IF ( b-a /= INT(memstride(2)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 2', line)
    ENDIF

    IF (SIZE(f,3) >= 2) THEN
      CALL mpi_get_address(f(1,1,2,1), b, ierror)
      IF ( b-a /= INT(memstride(3)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 3', line)
    ENDIF

    IF (SIZE(f,4) >= 2) THEN
      CALL mpi_get_address(f(1,1,1,2), b, ierror)
      IF ( b-a /= INT(memstride(4)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 4', line)
    ENDIF

  END SUBROUTINE check_mem_strides_4d

  SUBROUTINE check_mem_strides_3d(f, memstride, element_size, line)
    USE mpi
    CHARACTER(len=*), PARAMETER :: prefix = 'check_mem_strides_3d: '
    REAL(dp), INTENT(in) :: f(:,:,:)
    INTEGER, INTENT(in) :: memstride(:) ! vector of memory strides in units of element_size
    INTEGER, INTENT(in) :: element_size, line
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a, b
    INTEGER :: ierror

    IF (SIZE(memstride) /= 3) CALL die(prefix//'bad usage',__LINE__)

    IF (SIZE(f)<2) RETURN
    CALL mpi_get_address(f(1,1,1), a, ierror)

    IF (SIZE(f,1) >= 2) THEN
      CALL mpi_get_address(f(2,1,1), b, ierror)
      IF ( b-a /= INT(memstride(1)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 1', line)
    ENDIF

    IF (SIZE(f,2) >= 2) THEN
      CALL mpi_get_address(f(1,2,1), b, ierror)
      IF ( b-a /= INT(memstride(2)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 2', line)
    ENDIF

    IF (SIZE(f,3) >= 2) THEN
      CALL mpi_get_address(f(1,1,2), b, ierror)
      IF ( b-a /= INT(memstride(3)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 3', line)
    ENDIF
  END SUBROUTINE check_mem_strides_3d

  SUBROUTINE check_mem_strides_2d(f, memstride, element_size, line)
    USE mpi
    CHARACTER(len=*), PARAMETER :: prefix = 'check_mem_strides_2d: '
    REAL(dp), INTENT(in) :: f(:,:)
    INTEGER, INTENT(in) :: memstride(:) ! vector of memory strides in units of element_size
    INTEGER, INTENT(in) :: element_size, line
    INTEGER(KIND=MPI_ADDRESS_KIND) :: a, b
    INTEGER :: ierror

    IF (SIZE(memstride) /= 2) CALL die(prefix//'bad usage',__LINE__)

    IF (SIZE(f)<2) RETURN
    CALL mpi_get_address(f(1,1), a, ierror)

    IF (SIZE(f,1) >= 2) THEN
      CALL mpi_get_address(f(2,1), b, ierror)
      IF ( b-a /= INT(memstride(1)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 1', line)
    ENDIF

    IF (SIZE(f,2) >= 2) THEN
      CALL mpi_get_address(f(1,2), b, ierror)
      IF ( b-a /= INT(memstride(2)*element_size, MPI_ADDRESS_KIND)) &
           & CALL die(prefix//'bad memory stride in dim 2', line)
    ENDIF
  END SUBROUTINE check_mem_strides_2d

#endif /* CHECK_MEM_STRIDES */

#endif /* HAVE_YAXT */

END MODULE mo_echam_yaxt
