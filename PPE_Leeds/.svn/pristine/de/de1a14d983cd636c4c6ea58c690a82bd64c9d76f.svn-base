!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!-----------------------------------------------------------------------------
! Control of gather: unpacking and receive maby overlapped. 
! 
! The Cray XT3 and Portland Group f95 compiler has significant 
! problems with the existance of buffers with MPI_Isend and MPI_Irecv,
! so we use the old blocking MPI_Send/MPI_Recv scheme. Buffers are copied
! to temporary arrays without any need.
! 
#if defined (__XT3__) || defined (__PGI)
#define __ASYNC_GATHER 1
#define __ASYNC_GATHER_ANY 1
#endif
!
MODULE mo_jsbach_transpose
  !
  ! this module holds the global distribution and transposition routines:
  !   global distribution routines (global <-> pe)
  !   transposition routines (pe <-> pe)
  !
  ! Authors:
  !
  ! A. Rhodin, MPI,     August    1999, original source
  ! A. Rhodin, DWD/MPI, September 2001, FFSL transposition
  ! A. Rhodin, DWD/MPI, April     2002, blocking (nproma)
  ! L. Kornblueh, MPI, July 2004, replaced buffer allocation methode to reduce 
  !      allocates and allow for later usage with NEC SX GMEM and MPI-2 single
  !      sided communication to improve performance 
  ! R. Smith, and L. Kornblueh, MPI, October 2004, improvement in gather_gp by 
  !      implementing MPI_Wait_any functionality 
  ! L. Kornblueh, MPI, November 2004, removed some minor bugs related to the
  !      optimizations affecting code parts not used in standard ECHAM5
  ! L. Kornblueh, MPI, February 2005, work on buffering problems appearing 
  !      with different compiler.
  ! L. Kornblueh, MPI, February 2008, incorporated optimizations from
  !      NEC and IBM. Additional cleanup of the code.
  !
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish
  USE mo_mpi,           ONLY: p_nprocs,         &! number of processors
                              p_pe,             &! this processor
                              p_io,             &! processor which performs I/O
#ifdef __XT3__
                              p_barrier,        &! MPI_barrrier routine
#endif
#ifdef __ASYNC_GATHER
                              p_irecv,          &! MPI_irecv    routine
                              p_wait,           &! MPI_wait     routine
#endif
#ifdef __ASYNC_GATHER_ANY
                              p_wait_any,       &! MPI_wait_any routine
#endif
                              p_send,           &! MPI_send     routine
                              p_recv             ! MPI_recv     routine

  USE mo_decomposition, ONLY: pe_decomposed,    &! decomposition table data type
                              debug_parallel     ! debug flag
  USE mo_transpose,     ONLY: indx, reorder

  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ! public routines:
  !
  !
  !   scatter : global arrays -> local arrays
  !
  PUBLIC :: jsb_scatter_gp3 ! global grid point field -> local pe's (nlon,nlev,nlat)
  PUBLIC :: jsb_scatter_gp2 ! global grid point field -> local pe's (nlon,nlat)
  !
  !
  !   gather  : global arrays <- local arrays
  !
  PUBLIC :: jsb_gather_gp3 ! global grid point field <- local pe's (nlon,nlev,nlat)
  PUBLIC :: jsb_gather_gp2 ! global grid point field <- local pe's (nlon,nlat)
  !
  !
  ! define tags
  !
  INTEGER, PARAMETER :: tag_scatter_gp   = 100
  INTEGER, PARAMETER :: tag_gather_gp    = 101
  INTEGER, PARAMETER :: tag_scatter_gpc  = 103
  INTEGER, PARAMETER :: tag_gather_gpc   = 104
  !
  !====================================================================== 
  !
  ! gather buffer
  !
#ifdef __ASYNC_GATHER
  !
  TYPE gather_buffer
    REAL(dp), ALLOCATABLE :: receive2d(:,:)
    REAL(dp), ALLOCATABLE :: receive3d(:,:,:)
  END TYPE gather_buffer
  !
  TYPE(gather_buffer), ALLOCATABLE :: gather_array(:) 
  !
#endif
  !
  !====================================================================== 
  !
  ! transpose buffer
  !
  TYPE transpose_buffer
    REAL(dp), ALLOCATABLE :: send_buffer(:,:,:,:)
    REAL(dp), ALLOCATABLE :: recv_buffer(:,:,:,:)
    ! 3 dims for gp-fs and fs-gp, and fs-ls and ls-fs
    ! 2 dims for ls-sp and sp-ls (set last index to 1 for allocate)
    REAL(dp), ALLOCATABLE :: send_buffer0(:,:,:)
    REAL(dp), ALLOCATABLE :: recv_buffer0(:,:,:)
  END TYPE transpose_buffer
  !
!==============================================================================
CONTAINS

  SUBROUTINE jsb_scatter_gp2 (gl, lc)
    !
    ! send global 2D grid point field to local pe's (nlon,nlat)
    !
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    !
    REAL(dp), POINTER                :: gl   (:,:) ! global field
    REAL(dp), TARGET ,INTENT(out)    :: lc   (:,:) ! local  field
    !
    ! Data structure:
    !
    !   global grid point field GL: (1:nlon (+2), 1:nlat)
    !
    !   local  grid point field LC: (1:nglon(+2), 1:nglat)
    !
    !   The actual size of the first index may or may be not larger than 
    !   NLON or NGLON
    !
    !   The local grid point field LC covers two distinct areas at opposite
    !   sides of the globe. The array elements correspond with the global
    !   field as follows:
    !
    !   LC (1:nglon,1:nglat/2 ) = GL (glons(1):glone(1),glats(1):glate(1))
    !   LC (1:nglon,nglat/2+1:) = GL (glons(2):glone(2),glats(2):glate(2))
    !
    ! local variables
    !
    REAL(dp), ALLOCATABLE :: buf (:,:) ! buffer
    INTEGER               :: imype     ! index of this PE
    INTEGER               :: i         ! loop index
    INTEGER               :: pe        ! processor to communicate with
    INTEGER               :: nlon      ! global number of longitudes
    LOGICAL               :: lreg      ! flag for regular grid
    REAL(dp), POINTER     :: lcb (:,:) ! pointer/temporary buffer
    REAL(dp)        ,POINTER :: gl_comp(:)  ! global compressed field

    imype = indx (p_pe, gl_dc)

    !
    ! Compressed decomposition
    IF (ASSOCIATED(gl_dc(imype)%mask)) THEN
       IF (p_pe == p_io) THEN
          ALLOCATE(gl_comp(gl_dc(imype)%npts))
          gl_comp = PACK(gl, MASK=gl_dc(imype)%mask)
       ENDIF
       CALL jsb_scatter_gpc1to2(gl_comp, lc, gl_dc)
       IF (p_pe == p_io) DEALLOCATE(gl_comp)
       RETURN
    ENDIF

    !
    ! set, allocate local variables
    !
    nlon  = gl_dc(imype)% nlon
    lreg  = gl_dc(imype)% lreg

    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% nglon, &
                    gl_dc(imype)% nglat))
    ENDIF
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        ALLOCATE (buf (gl_dc(i)% nglon, gl_dc(i)% nglat))
        !
        ! pack first segment
        !
        buf(:,:gl_dc(i)% nglh(1)) =                    &
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1), &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1))
        !
        ! pack second segment
        !
        IF (gl_dc(i)% nglh(2)>0) THEN
          IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
            buf(:,gl_dc(i)% nglat/2+1:) =                  &
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2), &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ELSE
             !
            ! pack second segment, split in longitudes
            !
            buf(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:) = &
              gl (gl_dc(i)% glons(2) : nlon,                       &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
            buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,             &
                  gl_dc(i)% nglat/2+1:) =                          &
              gl (1: gl_dc(i)% glone(2),                           &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ENDIF
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_gp)
        ELSE
          lcb(:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      CALL p_recv (lcb(:,:), p_io, tag_scatter_gp)
    END IF
    IF (.NOT.lreg) THEN
      CALL reorder (lc,lcb)
      DEALLOCATE (lcb)
    ENDIF
  END SUBROUTINE jsb_scatter_gp2
!------------------------------------------------------------------------------
  SUBROUTINE jsb_scatter_gp3 (gl, lc)
    !
    ! send global 3D grid point field to local pe's (nlon,nlev,nlat)
    !
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    !
    REAL(dp), POINTER                :: gl   (:,:,:) ! global field
    REAL(dp), TARGET  ,INTENT(out)   :: lc   (:,:,:) ! local  field
    !
    ! Data structure:
    !
    !   global grid point field GL: (1:nlon (+2), 1:nlev, 1:nlat)
    !
    !   local  grid point field LC: (1:nglon(+2), 1:nlev, 1:nglat)
    !
    !   The actual size of the first index may or may be not larger than 
    !   NLON or NGLON
    !
    !   The local grid point field LC covers two distinct areas at opposite
    !   sides of the globe. The array elements correspond with the global
    !   field as follows:
    !
    !   LC (1:nglon,:,1:nglat/2 ) = GL (glons(1):glone(1),:,glats(1):glate(1))
    !   LC (1:nglon,:,nglat/2+1:) = GL (glons(2):glone(2),:,glats(2):glate(2))
    !
    ! local variables
    !
    REAL(dp), ALLOCATABLE :: buf (:,:,:) ! buffer
    INTEGER               :: i           ! loop index
    INTEGER               :: pe          ! processor to communicate with
    INTEGER               :: nlon        ! global number of longitudes
!    INTEGER               :: nglon       ! local
    INTEGER               :: imype       ! index of this pe
    LOGICAL               :: lreg        ! flag for regular grid
    REAL(dp), POINTER     :: lcb (:,:,:) ! pointer/temporary buffer
    REAL(dp)        ,POINTER :: gl_comp(:,:)  ! global compressed field

    imype = indx (p_pe, gl_dc)
    !
    ! Compressed decomposition
    IF (ASSOCIATED(gl_dc(imype)%mask)) THEN
       IF (p_pe == p_io) THEN
          ALLOCATE(gl_comp(gl_dc(imype)%npts,SIZE(gl,2)))
          DO i=1,SIZE(gl,2)
             gl_comp(:,i) = PACK(gl(:,i,:), MASK=gl_dc(imype)%mask)
          ENDDO
       ENDIF
       CALL jsb_scatter_gpc2to3(gl_comp, lc, gl_dc)
       IF (p_pe == p_io) DEALLOCATE(gl_comp)
       RETURN
    ENDIF

    !
    ! set, allocate local variables
    !
    nlon = gl_dc(1)% nlon
    lreg  = gl_dc(imype)% lreg
    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% nglon,&
                    SIZE(lc,2),         &
                    gl_dc(imype)% nglat))
    ENDIF
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        ALLOCATE (buf (gl_dc(i)% nglon, SIZE(gl,2), gl_dc(i)% nglat))
        !
        ! pack first segment
        !
        buf(:,:,:gl_dc(i)% nglh(1)) =                   &
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1), : , &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1))
        !
        ! pack second segment
        !
        IF (gl_dc(i)% nglh(2)>0) THEN
          IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
            buf(:,:,gl_dc(i)% nglh(1)+1:) =                 &
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2), : , &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ELSE
            !
            ! pack second segment, split in longitudes
            !
            buf(:nlon-gl_dc(i)% glons(2)+1,:,gl_dc(i)% nglh(1)+1:) = &
              gl (gl_dc(i)% glons(2) : nlon, : ,                     &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
            buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,:,&
                  gl_dc(i)% nglh(1)+1:) = &
              gl (1: gl_dc(i)% glone(2), : ,          &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ENDIF
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_gp)
        ELSE
          lcb(:,:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      CALL p_recv (lcb(:,:,:), p_io, tag_scatter_gp)
    END IF
    IF (.NOT.lreg) THEN
      CALL reorder (lc,lcb)
      DEALLOCATE (lcb)
    ENDIF
  END SUBROUTINE jsb_scatter_gp3
!------------------------------------------------------------------------------
  SUBROUTINE jsb_scatter_gpc1to2 (gl, lc, gl_dc)
  !
  ! send global compressed grid point field (npts) to local pe's (nproma,ngpblks)
  !
  REAL(dp)             ,POINTER     :: gl   (:)   ! global field
  REAL(dp)     ,TARGET ,INTENT(out) :: lc   (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition
    !
    ! Data structure:
    !
    !   global grid point field GL: (1:nlon (+2), 1:nlat)
    !
    !   local  grid point field LC: (1:nglon(+2), 1:nglat)
    !
    !   The actual size of the first index may or may be not larger than
    !   NLON or NGLON
    !
    !   The local grid point field LC covers two distinct areas at opposite
    !   sides of the globe. The array elements correspond with the global
    !   field as follows:
    !
    !   LC (1:nglon,1:nglat/2 ) = GL (glons(1):glone(1),glats(1):glate(1))
    !   LC (1:nglon,nglat/2+1:) = GL (glons(2):glone(2),glats(2):glate(2))
    !
    ! local variables
    !
    REAL(dp),ALLOCATABLE :: buf (:) ! buffer
    INTEGER              :: imype     ! index of this PE
    INTEGER              :: i         ! loop index
    INTEGER              :: pe        ! processor to communicate with
    INTEGER              :: npts      ! global number of longitudes
    REAL(dp)    ,POINTER :: lcb (:) ! pointer/temporary buffer
    !
    ! set, allocate local variables
    !
    imype = indx (p_pe, gl_dc)
    npts  = gl_dc(imype)% npts
    IF (gl_dc(imype)%lreg) &
         CALL finish('jsb_scatter_gpc12', 'Regular grid not allowed')
    IF (npts < 0) &
         CALL finish('jsb_scatter_gpc12', 'This is not a compressed decomposition')

    ALLOCATE (lcb(gl_dc(imype)% ngpts))
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        ALLOCATE (buf (gl_dc(i)% ngpts))

        buf(:) = gl (gl_dc(i)% gptss : gl_dc(i)% gptse)
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_gpc)
        ELSE
          lcb(:) = buf
!         lcb(nglon+1:     ,:) = 0._dp
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      CALL p_recv (lcb(:), p_io, tag_scatter_gpc)
!     lcb(nglon+1:,:) = 0._dp
    END IF
    CALL reorder (lc,lcb)
    DEALLOCATE (lcb)
  END SUBROUTINE jsb_scatter_gpc1to2
!------------------------------------------------------------------------------
  SUBROUTINE jsb_scatter_gpc2to3 (gl, lc, gl_dc)
  !
  ! send global compressed grid point field (npts) to local pe's (nproma,ngpblks)
  !
  REAL(dp)             ,POINTER     :: gl   (:,:)   ! global field
  REAL(dp)     ,TARGET ,INTENT(out) :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
    !
    ! Data structure:
    !
    !   global grid point field GL: (1:nlon (+2), 1:nlat)
    !
    !   local  grid point field LC: (1:nglon(+2), 1:nglat)
    !
    !   The actual size of the first index may or may be not larger than
    !   NLON or NGLON
    !
    !   The local grid point field LC covers two distinct areas at opposite
    !   sides of the globe. The array elements correspond with the global
    !   field as follows:
    !
    !   LC (1:nglon,1:nglat/2 ) = GL (glons(1):glone(1),glats(1):glate(1))
    !   LC (1:nglon,nglat/2+1:) = GL (glons(2):glone(2),glats(2):glate(2))
    !
    ! local variables
    !
    REAL(dp),ALLOCATABLE :: buf (:,:) ! buffer
    INTEGER              :: imype     ! index of this PE
    INTEGER              :: i         ! loop index
    INTEGER              :: pe        ! processor to communicate with
    INTEGER              :: npts      ! global number of longitudes
    REAL(dp)    ,POINTER :: lcb (:,:) ! pointer/temporary buffer
    !
    ! set, allocate local variables
    !
    imype = indx (p_pe, gl_dc)
    npts  = gl_dc(imype)% npts
    IF (gl_dc(imype)%lreg) &
         CALL finish('jsb_scatter_gpc12', 'Regular grid not allowed')
    IF (npts < 0) &
         CALL finish('jsb_scatter_gpc12', 'This is not a compressed decomposition')

    ALLOCATE (lcb(gl_dc(imype)% ngpts, SIZE(gl,2)))
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        ALLOCATE (buf (gl_dc(i)% ngpts, SIZE(gl,2)))

        buf(:,:) = gl (gl_dc(i)% gptss : gl_dc(i)% gptse,:)
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_gpc)
        ELSE
          lcb(:,:) = buf(:,:)
!         lcb(nglon+1:     ,:) = 0._dp
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      CALL p_recv (lcb(:,:), p_io, tag_scatter_gpc)
!     lcb(nglon+1:,:) = 0._dp
    END IF
    CALL reorder (lc,lcb)
    DEALLOCATE (lcb)
  END SUBROUTINE jsb_scatter_gpc2to3
!------------------------------------------------------------------------------
  SUBROUTINE jsb_gather_gp2 (gl, lc, source)
    !
    ! receive global grid point field from local pe's (nlon,nlat)
    !
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    USE mo_jsbach,        ONLY: missing_value
    !
    REAL(dp), POINTER                :: gl   (:,:)   ! global field
    REAL(dp), TARGET ,INTENT(in)     :: lc   (:,:)   ! local  field
    INTEGER, OPTIONAL    ,INTENT(in) :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
#ifndef __ASYNC_GATHER
    REAL(dp), ALLOCATABLE :: buf (:,:)   ! buffer
    INTEGER               :: nglon       ! local number of longitudes
#endif
    INTEGER               :: nlon        ! global number of longitudes
    INTEGER               :: i           ! loop index
#ifdef __ASYNC_GATHER_ANY
    INTEGER               :: ii          ! loop index
#endif
    INTEGER               :: pe          ! processor to communicate with
    INTEGER               :: imype       ! index of this pe
    INTEGER               :: src         ! source 
    LOGICAL               :: lreg        ! flag for regular grid
    REAL(dp), POINTER     :: lcb (:,:)   ! pointer/temporary buffer
#ifdef __XT3__
    INTEGER               :: iplow, iphigh
    INTEGER, PARAMETER    :: ipstrip=128
#endif 
    REAL(dp)        ,POINTER :: gl_comp(:)  ! global compressed field

    imype = indx (p_pe, gl_dc)

    !
    ! Compressed decomposition
    IF (ASSOCIATED(gl_dc(imype)%mask)) THEN
       IF (p_pe == p_io) ALLOCATE(gl_comp(gl_dc(imype)%npts))
       CALL gather_gpc2to1(gl_comp, lc, gl_dc, source)
       IF (p_pe == p_io) THEN
          gl = UNPACK(gl_comp, gl_dc(imype)%mask, missing_value)
          DEALLOCATE(gl_comp)
       ENDIF
       RETURN
    ENDIF
    !
    ! set and allocate local variables
    !
    ! for parallel debugging mode
    !
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    nlon  = gl_dc(imype)% nlon
    !
    lreg  = gl_dc(imype)% lreg
    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% nglon,&
                    gl_dc(imype)% nglat))
      CALL reorder (lcb,lc)
    ENDIF
    !
#ifdef __XT3__
    do iplow = 0, p_nprocs-1, ipstrip
    iphigh = min(iplow+ipstrip-1,p_nprocs-1)
    call p_barrier

#endif 
    IF (p_pe /= p_io) THEN
      !
      ! send if pe /= p_io
      !
#ifdef __XT3__
      if ( gl_dc(p_pe)%pe >= iplow .and. gl_dc(p_pe)%pe <= iphigh ) &
#endif 
      CALL p_send (lcb(:,:), p_io, tag_gather_gp)
      !
    ELSE
      !
      ! receive
      !
#ifdef __ASYNC_GATHER
      IF (.NOT. ALLOCATED(gather_array)) THEN
        ALLOCATE(gather_array(0:p_nprocs-1))
      ENDIF
      !
      IF (.NOT. ALLOCATED(gather_array(0)%receive2d)) THEN
      !
        DO i = 1, p_nprocs
#ifdef __XT3__
           if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) &
#endif 
          ALLOCATE(gather_array(gl_dc(i)%pe)%receive2d(gl_dc(i)%nglon, gl_dc(i)%nglat))
        ENDDO
      ENDIF
      !
      DO i = 1, p_nprocs
#ifdef __XT3__
        if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) then
#endif 
        IF (gl_dc(i)%pe /= gl_dc(imype)%pe) THEN
          CALL p_irecv(gather_array(gl_dc(i)%pe)%receive2d, gl_dc(i)%pe, tag_gather_gp)
        ELSE
          gather_array(gl_dc(i)%pe)%receive2d(:,:) = lcb(:,:)
        ENDIF
#ifdef __XT3__
        endif
#endif 
      ENDDO
       
#ifdef __ASYNC_GATHER_ANY
      DO ii = 1, p_nprocs
        IF (ii == 1) THEN
          pe = p_pe
        ELSE
          CALL p_wait_any(pe)
        END IF
        i = indx(pe, gl_dc)
#else
      CALL p_wait
       
      DO i = 1, p_nprocs
        pe = gl_dc(i)%pe
#endif
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            gather_array(pe)%receive2d(:,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                gather_array(pe)%receive2d(:,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% glons(2) : nlon,                       &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) =       &
                  gather_array(pe)%receive2d(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:)
              gl (1: gl_dc(i)% glone(2),                        &
                     gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                  gather_array(pe)%receive2d(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
      END DO
#else
      DO i = 1, p_nprocs
#ifdef __XT3__
        if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) then
#endif 
        pe    = gl_dc(i)% pe
        nglon = gl_dc(i)% nglon
        ALLOCATE (buf (nglon, gl_dc(i)% nglat))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_gp)
        ELSE
          buf = lcb(:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            buf(:,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(:,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% glons(2) : nlon,                       &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) =       &
                buf(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:)
              gl (1: gl_dc(i)% glone(2),                        &
                     gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,      &
                    gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
        DEALLOCATE (buf)
#ifdef __XT3__
      ENDIF
#endif 
      END DO
#endif
    ENDIF
#ifdef __XT3__
    ENDDO
#endif 
    !
    IF (lreg) THEN
      NULLIFY (lcb)
    ELSE
      DEALLOCATE (lcb)
    ENDIF
    !
  END SUBROUTINE jsb_gather_gp2
!------------------------------------------------------------------------------
  SUBROUTINE jsb_gather_gp3 (gl, lc, source)
    !
    ! receive global grid point field from local pe's (nlon,nlev,nlat)
    !
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    USE mo_jsbach,        ONLY: missing_value
    !
    REAL(dp), POINTER                :: gl   (:,:,:) ! global field
    REAL(dp), TARGET ,INTENT(in)     :: lc   (:,:,:) ! local  field
    INTEGER, OPTIONAL    ,INTENT(in) :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
#ifndef __ASYNC_GATHER
    REAL(dp), ALLOCATABLE :: buf (:,:,:) ! buffer
    INTEGER               :: nglon       ! local number of longitudes
#endif 
    INTEGER               :: nlon        ! global number of longitudes
    INTEGER               :: i           ! loop index
#ifdef __ASYNC_GATHER_ANY
    INTEGER               :: ii          ! loop index
#endif
    INTEGER               :: nk          ! local number of levels
    INTEGER, SAVE         :: nks = -1    ! save for control reasons
    INTEGER               :: pe          ! processor to communicate with
    INTEGER               :: imype       ! index of this pe
    INTEGER               :: src         ! source 
    LOGICAL               :: lreg        ! flag for regular grid
    REAL(dp), POINTER     :: lcb (:,:,:) ! pointer/temporary buffer
#ifdef __XT3__
    INTEGER               :: iplow, iphigh
    INTEGER, PARAMETER    :: ipstrip=128
#endif 
    REAL(dp)        ,POINTER :: gl_comp(:,:) ! global compressed field

    imype = indx (p_pe, gl_dc)
    !
    ! Compressed decomposition
    IF (ASSOCIATED(gl_dc(imype)%mask)) THEN
       IF (p_pe == p_io) ALLOCATE(gl_comp(gl_dc(imype)%npts, SIZE(gl,2)))
       CALL gather_gpc3to2(gl_comp, lc, gl_dc, source)
       IF (p_pe == p_io) THEN
          DO i=1,SIZE(gl,2)
             gl(:,i,:) = UNPACK(gl_comp(:,i), gl_dc(imype)%mask, missing_value)
          ENDDO
          DEALLOCATE(gl_comp)
       ENDIF
       RETURN
    ENDIF
    !
    ! set, allocate local variables
    !
    ! for parallel debugging mode
    !
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! second dimension adjustment
    !
    nk = SIZE(lc,2)
    if (nks == -1) nks = nk
    !
    lreg  = gl_dc(imype)% lreg
    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% nglon,&
                    nk,                 & 
                    gl_dc(imype)% nglat))
      CALL reorder (lcb,lc)
    ENDIF
    !
    ! send if pe /= p_io
    !
#ifdef __XT3__
    do iplow = 0, p_nprocs-1, ipstrip
    iphigh = min(iplow+ipstrip-1,p_nprocs-1)
    call p_barrier

#endif 
    IF (p_pe /= p_io) THEN
#ifdef __XT3__
      if ( gl_dc(p_pe)%pe >= iplow .and. gl_dc(p_pe)%pe <= iphigh ) &
#endif 
      CALL p_send (lcb(:,:,:), p_io, tag_gather_gp)
    ELSE
      !
      ! receive
      !
#ifdef __ASYNC_GATHER
      IF (.NOT. ALLOCATED(gather_array)) THEN
        ALLOCATE(gather_array(0:p_nprocs-1))
      ENDIF
      !
      IF (.NOT. ALLOCATED(gather_array(0)%receive3d)) THEN
        DO i = 1, p_nprocs
          ALLOCATE(gather_array(gl_dc(i)%pe)%receive3d(gl_dc(i)%nglon, &
                                nk,                                    &
                                gl_dc(i)%nglat))
        ENDDO
      ENDIF  
      !
      ! eventually adjust size
      !
      IF (nks /= nk) THEN
        DO i = 1, p_nprocs
          DEALLOCATE (gather_array(gl_dc(i)%pe)%receive3d)
          ALLOCATE(gather_array(gl_dc(i)%pe)%receive3d(gl_dc(i)%nglon, &
                                nk,                                    &
                                gl_dc(i)%nglat))
        ENDDO
        nks = nk
      ENDIF
      !
      ! receive
      !
      IF (gl_dc(imype)%pe == p_io) THEN
        gather_array(gl_dc(imype)%pe)%receive3d(:,1:nk,:) = lcb(:,1:nk,:)
      ENDIF
      DO i = 1, p_nprocs
        IF (gl_dc(i)%pe /= gl_dc(imype)%pe) THEN
          CALL p_irecv(gather_array(gl_dc(i)%pe)%receive3d, gl_dc(i)%pe, &
                       tag_gather_gp)
        ENDIF
      ENDDO

#ifdef __ASYNC_GATHER_ANY
      DO ii = 1, p_nprocs
        IF (gl_dc(ii)%pe == p_io) THEN
          pe = p_pe
        ELSE
          CALL p_wait_any(pe)
        END IF
        i = indx(pe, gl_dc)
#else
      CALL p_wait
      
      DO i = 1, p_nprocs
        pe = gl_dc(i)%pe
#endif
        nlon = gl_dc(i)% nlon
        IF(src == -1 .OR. (src == 0 .AND. p_io == pe)      &
                     .OR. (src == 1 .AND. p_io /= pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),1:nk,   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            gather_array(pe)%receive3d(:,1:nk,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),1:nk,   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                gather_array(pe)%receive3d(:,1:nk,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl(gl_dc(i)%glons(2):nlon,                                      &
                 1:nk,                                                        &
                 gl_dc(i)%glats(2):gl_dc(i)%glate(2))                         &
              =                                                               &
              gather_array(pe)%receive3d(:nlon-gl_dc(i)%glons(2)+1,           &
                                         1:nk,                                &
                                         gl_dc(i)% nglat/2+1:)
              !
              gl (1:gl_dc(i)%glone(2),                                        &
                  1:nk,                                                       &
                  gl_dc(i)%glats(2):gl_dc(i)%glate(2))                        &
              =                                                               &
              gather_array(pe)%receive3d(gl_dc(i)%nglon-gl_dc(i)%glone(2)+1:, &
                                         1:nk,                                &
                                         gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
      END DO
#else
      nlon  = gl_dc(imype)% nlon
      DO i = 1, p_nprocs
#ifdef __XT3__
        if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) then
#endif 
        pe    = gl_dc(i)% pe
        nglon = gl_dc(i)% nglon
        ALLOCATE (buf (nglon, SIZE(gl,2), gl_dc(i)% nglat))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_gp)
        ELSE
          buf = lcb(:,:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),:,   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            buf(:,:,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),:,   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(:,:,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% glons(2) : nlon, : ,                     &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) =       &
                buf(:nlon-gl_dc(i)% glons(2)+1,:,gl_dc(i)% nglat/2+1:)
              gl (1: gl_dc(i)% glone(2), : ,                      &
                     gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,:,      &
                    gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
        DEALLOCATE (buf)
#ifdef __XT3__
        ENDIF
#endif 
      END DO
#endif
    ENDIF
#ifdef __XT3__
    ENDDO
#endif 
    IF (lreg) THEN
      NULLIFY (lcb)
    ELSE
      DEALLOCATE (lcb)
    ENDIF
  END SUBROUTINE jsb_gather_gp3
!------------------------------------------------------------------------------
  SUBROUTINE gather_gpc2to1 (gl, lc, gl_dc, source)
  !
  ! receive global compressed grid point field (npts) from local pe's (nproma,ngpblks)
  !
  REAL(dp)             ,POINTER     :: gl   (:)     ! global field
  REAL(dp)     ,TARGET ,INTENT(in)  :: lc   (:,:)   ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
    REAL(dp),ALLOCATABLE :: buf (:)     ! buffer
    INTEGER              :: i           ! loop index
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: npts        ! global number of compressed points
    INTEGER              :: ngpts       ! local number of compressed points
    INTEGER              :: imype       ! index of this pe
    INTEGER              :: src         ! source
    REAL(dp)    ,POINTER :: lcb (:)     ! pointer/temporary buffer
    !
    ! set, allocate local variables
    !
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    imype = indx (p_pe, gl_dc)
    npts  = gl_dc(1)% npts
    IF (gl_dc(imype)%lreg) &
         CALL finish('gather_gpc2to1', 'Regular grid not allowed')
    IF (npts < 0) &
         CALL finish('gather_gpc2to1', 'This is not a compressed decomposition')

    ALLOCATE (lcb(gl_dc(imype)% ngpts))
    CALL reorder(lcb,lc)
    !
    ! send if pe /= p_io
    !
    IF (p_pe /= p_io) THEN
      ngpts = gl_dc(imype)% ngpts
      CALL p_send (lcb(:), p_io, tag_gather_gpc)
    ELSE
      DO i = 1, p_nprocs
        pe    = gl_dc(i)% pe
        ngpts = gl_dc(i)% ngpts
        ALLOCATE (buf (ngpts))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_gpc)
        ELSE
          buf = lcb(:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
           gl (gl_dc(i)% gptss : gl_dc(i)% gptse) = buf(:)
        ENDIF
        DEALLOCATE (buf)
      END DO
      !
      ! set elements with i>nlon to zero
      !
!     gl (nlon+1:,:) = 0._dp
    ENDIF
    DEALLOCATE (lcb)
  END SUBROUTINE gather_gpc2to1
!------------------------------------------------------------------------------
  SUBROUTINE gather_gpc3to2 (gl, lc, gl_dc, source)
  !
  ! receive global compressed grid point field (npts) from local pe's (nproma,ngpblks)
  !
  REAL(dp)             ,POINTER     :: gl   (:,:)     ! global field
  REAL(dp)     ,TARGET ,INTENT(in)  :: lc   (:,:,:)   ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)       ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source         ! source to gather from
    !                                                 ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
    REAL(dp),ALLOCATABLE :: buf (:,:)   ! buffer
    INTEGER              :: i           ! loop index
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: npts        ! global number of compressed points
    INTEGER              :: ngpts       ! local number of compressed points
    INTEGER              :: imype       ! index of this pe
    INTEGER              :: src         ! source
    REAL(dp)    ,POINTER :: lcb (:,:)   ! pointer/temporary buffer
    !
    ! set, allocate local variables
    !
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    imype = indx (p_pe, gl_dc)
    npts  = gl_dc(1)% npts
    IF (gl_dc(imype)%lreg) &
         CALL finish('gather_gpc3to2', 'Regular grid not allowed')
    IF (npts < 0) &
         CALL finish('gather_gpc3to2', 'This is not a compressed decomposition')

    ALLOCATE (lcb(gl_dc(imype)% ngpts,SIZE(lc,2)))
    CALL reorder(lcb,lc)
    !
    ! send if pe /= p_io
    !
    IF (p_pe /= p_io) THEN
      ngpts = gl_dc(imype)% ngpts
      CALL p_send (lcb(:,:), p_io, tag_gather_gpc)
    ELSE
      DO i = 1, p_nprocs
        pe    = gl_dc(i)% pe
        ngpts = gl_dc(i)% ngpts
        ALLOCATE (buf (ngpts,SIZE(gl,2)))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf(:,:), pe, tag_gather_gpc)
        ELSE
          buf(:,:) = lcb(:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
           gl (gl_dc(i)% gptss : gl_dc(i)% gptse, :) = buf(:,:)
        ENDIF
        DEALLOCATE (buf)
      END DO
      !
      ! set elements with i>nlon to zero
      !
!     gl (nlon+1:,:) = 0._dp
    ENDIF
    DEALLOCATE (lcb)
  END SUBROUTINE gather_gpc3to2
!------------------------------------------------------------------------------
END MODULE mo_jsbach_transpose
