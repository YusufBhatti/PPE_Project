#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_global_op
!
! This module holds the global operations (max, min, sum, etc.)
! on gridpoint fields. The operators return the same value in
! parallel and seriell mode.
!
! Authors:
!
! A. Rhodin, MPI, August 1999, original source
!
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish
  USE mo_mpi,           ONLY: p_send,          &! MPI_send routine
                              p_recv,          &! MPI_recv routine
                              p_max, p_min,    &! reduce operators
                              p_communicator_a,&! communication groups
                              p_communicator_b,&!
                              p_communicator_d
  USE mo_decomposition, ONLY: dc               => local_decomposition, &
                              gl_dc            => global_decomposition
  USE mo_transpose,     ONLY: reorder

  IMPLICIT NONE

  PRIVATE

  !
  ! public global operations:
  !
  PUBLIC :: maxval_zonal  ! maximum value, reduce 1st  index in zonal direction
  PUBLIC :: maxval_latit  ! maximum value, reduce last index in latit.direction
  PUBLIC :: maxval_global ! maximum value, reduce 1st(zonal),last(latit.) index
  PUBLIC :: minval_zonal  ! minimum value, reduce 1st  index in zonal direction
  PUBLIC :: minval_latit  ! minimum value, reduce last index in latit.direction
  PUBLIC :: minval_global ! minimum value, reduce 1st(zonal),last(latit.) index
  PUBLIC :: sum_zonal     ! sum,           reduce 1st  index in zonal direction
  PUBLIC :: sum_latit     ! sum,           reduce last index in latit.direction
  PUBLIC :: sum_global    ! sum,           reduce 1st(zonal),last(latit.) index
  PUBLIC :: sum_zonal_sl  ! zonal sum used in mass fixer: reduce both indices 
  PUBLIC :: sum_latit_sl  ! latitudinal sum used in mass fixer: reduce 2nd ind.
  !
  ! interfaces
  !                                 ! index bounds of argument x
  INTERFACE maxval_zonal
    MODULE PROCEDURE maxval_zonal1  ! (nlon)
    MODULE PROCEDURE maxval_zonal2  ! (nlon, any)
    MODULE PROCEDURE maxval_zonal3  ! (nlon, any, any)
  END INTERFACE

  INTERFACE maxval_latit
    MODULE PROCEDURE maxval_latit1  ! (     nlat)
    MODULE PROCEDURE maxval_latit2  ! (any, nlat)
  END INTERFACE

  INTERFACE maxval_global
    MODULE PROCEDURE maxval_global2 ! (nlon,      nlat)
    MODULE PROCEDURE maxval_global3 ! (nlon, any, nlat)
  END INTERFACE

  INTERFACE minval_zonal
    MODULE PROCEDURE minval_zonal1  ! (nlon)
    MODULE PROCEDURE minval_zonal2  ! (nlon, any)
    MODULE PROCEDURE minval_zonal3  ! (nlon, any, any)
  END INTERFACE

  INTERFACE minval_latit
    MODULE PROCEDURE minval_latit1  ! (     nlat)
    MODULE PROCEDURE minval_latit2  ! (any, nlat)
  END INTERFACE

  INTERFACE minval_global
    MODULE PROCEDURE minval_global2 ! (nlon,      nlat)
    MODULE PROCEDURE minval_global3 ! (nlon, any, nlat)
  END INTERFACE

  INTERFACE sum_zonal
    MODULE PROCEDURE sum_zonal1     ! (nlon)
    MODULE PROCEDURE sum_zonal2     ! (nlon,      nlat or any)
    MODULE PROCEDURE sum_zonal3     ! (nlon, any, nlat or any)
  END INTERFACE

  INTERFACE sum_latit
    MODULE PROCEDURE sum_latit1     ! (     nlat)
    MODULE PROCEDURE sum_latit2     ! (any, nlat)
  END INTERFACE

  INTERFACE sum_global
    MODULE PROCEDURE sum_global2    ! (nlon,      nlat)
    MODULE PROCEDURE sum_global3    ! (nlon, any, nlat)
  END INTERFACE
  !
  ! define tags
  !
  INTEGER, PARAMETER :: tag_maxval_zonal  = 300
  INTEGER, PARAMETER :: tag_maxval_latit  = 301
  INTEGER, PARAMETER :: tag_maxval_global = 302
  INTEGER, PARAMETER :: tag_minval_zonal  = 310
  INTEGER, PARAMETER :: tag_minval_latit  = 311
  INTEGER, PARAMETER :: tag_minval_global = 312
  INTEGER, PARAMETER :: tag_sum_zonal     = 320
  INTEGER, PARAMETER :: tag_sum_latit     = 321
!==============================================================================
CONTAINS
!==============================================================================
  FUNCTION maxval_zonal1 (x)
  !
  ! reduce to zonal maximum (within set_a communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:)
  REAL(dp)              :: maxval_zonal1
    maxval_zonal1 = p_max(MAXVAL(x(:dc%nglon)),p_communicator_a)
  END FUNCTION maxval_zonal1
!------------------------------------------------------------------------------
  FUNCTION maxval_zonal2 (x)
  !
  ! reduce to zonal maximum (within set_a communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:)
  REAL(dp)              :: maxval_zonal2 (SIZE(x,2))
    maxval_zonal2 = p_max(MAXVAL(x(:dc%nglon,:),dim=1),p_communicator_a)
  END FUNCTION maxval_zonal2
!------------------------------------------------------------------------------
  FUNCTION maxval_zonal3 (x)
  !
  ! reduce to zonal maximum (within set_a communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:,:)
  REAL(dp)              :: maxval_zonal3 (SIZE(x,2), SIZE(x,3))
    maxval_zonal3 = p_max(MAXVAL(x(:dc%nglon,:,:),dim=1),p_communicator_a)
  END FUNCTION maxval_zonal3
!------------------------------------------------------------------------------
  FUNCTION maxval_latit1 (x)
  !
  ! reduce to latit maximum (within set_b communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:)
  REAL(dp)              :: maxval_latit1
    maxval_latit1 = p_max(MAXVAL(x(:)),p_communicator_b)
  END FUNCTION maxval_latit1
!------------------------------------------------------------------------------
  FUNCTION maxval_latit2 (x)
  !
  ! reduce to latit maximum (within set_b communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:)
  REAL(dp)              :: maxval_latit2 (SIZE(x,1))
    maxval_latit2 = p_max(MAXVAL(x,dim=2),p_communicator_b)
  END FUNCTION maxval_latit2
!------------------------------------------------------------------------------
  FUNCTION maxval_global2 (x)
  !
  ! reduce to global maximum (within model communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:)
  REAL(dp)              :: maxval_global2
    maxval_global2 = p_max(MAXVAL(x(:dc%nglon,:)),p_communicator_d)
  END FUNCTION maxval_global2
!------------------------------------------------------------------------------
  FUNCTION maxval_global3 (x)
  !
  ! reduce to global maximum (within model communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:,:)
  REAL(dp)              :: maxval_global3 (SIZE(x,2))
    INTEGER :: i
    DO i=1,SIZE(x,2)
      maxval_global3(i) = MAXVAL(x(:dc%nglon,i,:))
    END DO
    maxval_global3 = p_max(maxval_global3,p_communicator_d)
  END FUNCTION maxval_global3
!==============================================================================
  FUNCTION minval_zonal1 (x)
  !
  ! reduce to zonal minimum (within set_a communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:)
  REAL(dp)              :: minval_zonal1
    minval_zonal1 = p_min(MINVAL(x(:dc%nglon)),p_communicator_a)
  END FUNCTION minval_zonal1
!------------------------------------------------------------------------------
  FUNCTION minval_zonal2 (x)
  !
  ! reduce to zonal minimum (within set_a communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:)
  REAL(dp)              :: minval_zonal2 (SIZE(x,2))
    minval_zonal2 = p_min(MINVAL(x(:dc%nglon,:),dim=1),p_communicator_a)
  END FUNCTION minval_zonal2
!------------------------------------------------------------------------------
  FUNCTION minval_zonal3 (x)
  !
  ! reduce to zonal minimum (within set_a communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:,:)
  REAL(dp)              :: minval_zonal3 (SIZE(x,2), SIZE(x,3))
    minval_zonal3 = p_min(MINVAL(x(:dc%nglon,:,:),dim=1),p_communicator_a)
  END FUNCTION minval_zonal3
!------------------------------------------------------------------------------
  FUNCTION minval_latit1 (x)
  !
  ! reduce to latit minimum (within set_b communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:)
  REAL(dp)              :: minval_latit1
    minval_latit1 = p_min(MINVAL(x(:)),p_communicator_b)
  END FUNCTION minval_latit1
!------------------------------------------------------------------------------
  FUNCTION minval_latit2 (x)
  !
  ! reduce to latit minimum (within set_b communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:)
  REAL(dp)              :: minval_latit2 (SIZE(x,1))
    minval_latit2 = p_min(MINVAL(x,dim=2),p_communicator_b)
  END FUNCTION minval_latit2
!------------------------------------------------------------------------------
  FUNCTION minval_global2 (x)
  !
  ! reduce to global minimum (within model communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:)
  REAL(dp)              :: minval_global2
    minval_global2 = p_min(MINVAL(x(:dc%nglon,:)),p_communicator_d)
  END FUNCTION minval_global2
!------------------------------------------------------------------------------
  FUNCTION minval_global3 (x)
  !
  ! reduce to global minimum (within model communicator group)
  !
  REAL(dp) ,INTENT (in) :: x(:,:,:)
  REAL(dp)              :: minval_global3 (SIZE(x,2))
    INTEGER :: i
    DO i=1,SIZE(x,2)
      minval_global3(i) = MINVAL(x(:dc%nglon,i,:))
    END DO
    minval_global3 = p_min(minval_global3,p_communicator_d)
  END FUNCTION minval_global3
!==============================================================================
  FUNCTION sum_zonal3 (x, jlat) RESULT (y)
  !
  ! zonal sum: reduce first index, bounds: (nlon[+x], any, any)
  !            jlat is reqired if 3rd index is not latitude index
  !
  REAL(dp)    ,INTENT (in)           :: x(:,:,:) ! array to sum over 1st index
  INTEGER ,INTENT (in), OPTIONAL :: jlat     ! local continuous row index
  REAL(dp) :: y (SIZE(x,2),SIZE(x,3))
    INTEGER :: i, k, set_a, pe, pe_root, glons(2), glone(2), shift
    LOGICAL :: first, root
    REAL(dp)    :: buf (dc%nlon,SIZE(x,2),SIZE(x,3))
    set_a = dc% set_a  ! local set a
    pe    = dc% pe     ! local pe
    shift = dc% nlon/2
    first =.TRUE.
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)% set_a == set_a) THEN
        !
        ! if root receive x
        !
        glons = gl_dc(i)% glons
        glone = gl_dc(i)% glone
        IF (first) THEN
          pe_root      = gl_dc(i)% pe  ! the first pe in the group is 'root'
          root         = (pe_root == pe)
          first        = .FALSE.
          IF(.NOT.root) EXIT
          buf(glons(1):glone(1),:,:) = x
        ELSE
          CALL p_recv (buf(glons(1):glone(1),:,:), &
                       gl_dc(i)% pe, tag_sum_zonal)  ! receve
        ENDIF
      ENDIF
    END DO
    IF (root) THEN
      !
      ! rotate longitudes
      !
      IF (PRESENT(jlat)) THEN
        IF (dc%glon(jlat) /= dc%glon(1)) buf = CSHIFT (buf, shift=shift, dim=1)
      ELSE IF (SIZE(x,3) == dc%nglat) THEN
        IF (glons(1)/=glons(2))    &
          buf(:,:,dc%nglat/2+1:) = &
          CSHIFT (buf(:,:,dc%nglat/2+1:), shift=shift, dim=1)
      ELSE
        CALL finish ('sum_zonal3','size(x,3) /= nglat')
      ENDIF
      !
      ! sum
      !
      y = 0.0_dp
      DO k = 1, dc% nlon
        y = y + buf(k,:,:)
      END DO
      !
      ! send result
      !
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_a == set_a) THEN
          IF (gl_dc(i)% pe /= pe) THEN
            CALL p_send (y, gl_dc(i)% pe, tag_sum_zonal)
          ENDIF
        ENDIF
      END DO
    ELSE
      !
      ! if not root send x and receive sum
      !
      CALL p_send (x(:dc%nglon,:,:), pe_root, tag_sum_zonal)
      CALL p_recv (y               , pe_root, tag_sum_zonal)
    ENDIF
  END FUNCTION sum_zonal3
!------------------------------------------------------------------------------
  FUNCTION sum_zonal2 (x, jlat) RESULT (y)
  !
  ! zonal sum: reduce first index, bounds: (nlon[+x], any)
  !            jlat is reqired if 2nd index is not latitude index
  !
  REAL(dp)    ,INTENT (in)           :: x(:,:) ! array to sum over 1st index
  INTEGER ,INTENT (in), OPTIONAL :: jlat   ! continuous local row index
  REAL(dp) :: y (SIZE(x,2))
    INTEGER :: i, k, set_a, pe, pe_root, glons(2), glone(2), shift
    LOGICAL :: first, root
    REAL(dp)    :: buf (dc%nlon,SIZE(x,2))
    set_a = dc% set_a  ! local set a
    pe    = dc% pe     ! local pe
    shift = dc% nlon/2
    first =.TRUE.
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)% set_a == set_a) THEN
        !
        ! if root receive x
        !
        glons = gl_dc(i)% glons
        glone = gl_dc(i)% glone
        IF (first) THEN
          pe_root      = gl_dc(i)% pe  ! the first pe in the group is 'root'
          root         = (pe_root == pe)
          first        = .FALSE.
          IF(.NOT.root) EXIT
          buf(glons(1):glone(1),:) = x
        ELSE
          CALL p_recv (buf(glons(1):glone(1),:), &
                       gl_dc(i)% pe, tag_sum_zonal)  ! receve
        ENDIF
      ENDIF
    END DO
    IF (root) THEN
      !
      ! rotate longitudes
      !
      IF (PRESENT(jlat)) THEN
        IF (dc%glon(jlat) /= dc%glon(1)) buf = CSHIFT (buf, shift=shift, dim=1)
      ELSE IF (SIZE(x,2) == dc%nglat) THEN
        IF (glons(1)/=glons(2))    &
          buf(:,dc%nglat/2+1:) = &
          CSHIFT (buf(:,dc%nglat/2+1:), shift=shift, dim=1)
      ELSE
        CALL finish ('sum_zonal2','size(x,2) /= nglat')
      ENDIF
      !
      ! sum
      !
      y = 0.0_dp
      DO k = 1, dc% nlon
        y = y + buf(k,:)
      END DO
      !
      ! send result
      !
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_a == set_a) THEN
          IF (gl_dc(i)% pe /= pe) THEN
            CALL p_send (y, gl_dc(i)% pe, tag_sum_zonal)
          ENDIF
        ENDIF
      END DO
    ELSE
      !
      ! if not root send x and receive sum
      !
      CALL p_send (x(:dc%nglon,:), pe_root, tag_sum_zonal)
      CALL p_recv (y             , pe_root, tag_sum_zonal)
    ENDIF
  END FUNCTION sum_zonal2
!------------------------------------------------------------------------------
  FUNCTION sum_zonal1 (x, jlat) RESULT (y)
  !
  ! zonal sum: sum x, index bound: (nlon[+x])
  !            jlat is the continuous local row index
  !
  REAL(dp)    ,INTENT (in) :: x(:)
  INTEGER ,INTENT (in) :: jlat
  REAL(dp) :: y
    INTEGER :: i, k, set_a, pe, pe_root, glons(2), glone(2), shift
    LOGICAL :: first, root
    REAL(dp)    :: buf (dc%nlon)
    set_a = dc% set_a  ! local set a
    pe    = dc% pe     ! local pe
    shift = dc% nlon/2
    first =.TRUE.
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)% set_a == set_a) THEN
        !
        ! if root receive x
        !
        glons = gl_dc(i)% glons
        glone = gl_dc(i)% glone
        IF (first) THEN
          pe_root      = gl_dc(i)% pe  ! the first pe in the group is 'root'
          root         = (pe_root == pe)
          first        = .FALSE.
          IF(.NOT.root) EXIT
          buf(glons(1):glone(1)) = x
        ELSE
          CALL p_recv (buf(glons(1):glone(1)), &
                       gl_dc(i)% pe, tag_sum_zonal)  ! receve
        ENDIF
      ENDIF
    END DO
    IF (root) THEN
      !
      ! rotate longitudes
      !
      IF (dc%glon(jlat) /= dc%glon(1)) buf = CSHIFT (buf, shift=shift)
      !
      ! sum
      !
      y = 0.0_dp
      DO k = 1, dc% nlon
        y = y + buf(k)
      END DO
      !
      ! send result
      !
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_a == set_a) THEN
          IF (gl_dc(i)% pe /= pe) THEN
            CALL p_send (y, gl_dc(i)% pe, tag_sum_zonal)
          ENDIF
        ENDIF
      END DO
    ELSE
      !
      ! if not root send x and receive sum
      !
      CALL p_send (x(:dc%nglon), pe_root, tag_sum_zonal)
      CALL p_recv (y           , pe_root, tag_sum_zonal)
    ENDIF
  END FUNCTION sum_zonal1
!==============================================================================
  FUNCTION sum_zonal_sl (x, jlat) RESULT (y)
  !
  ! zonal sum: reduce first index, bounds: (nlon[+x], any)
  !            jlat is reqired if 2nd index is not latitude index
  !
  REAL(dp)    ,INTENT (in)           :: x(:,:) ! array to sum over 1st index
  INTEGER ,INTENT (in)           :: jlat   ! continuous local row index
  REAL(dp) :: y
    INTEGER :: i, k, l, set_a, pe, pe_root, glons(2), glone(2), shift
    LOGICAL :: first, root
    REAL(dp)    :: buf (dc%nlon,SIZE(x,2))
    set_a = dc% set_a  ! local set a
    pe    = dc% pe     ! local pe
    shift = dc% nlon/2
    first =.TRUE.
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)% set_a == set_a) THEN
        !
        ! if root receive x
        !
        glons = gl_dc(i)% glons
        glone = gl_dc(i)% glone
        IF (first) THEN
          pe_root      = gl_dc(i)% pe  ! the first pe in the group is 'root'
          root         = (pe_root == pe)
          first        = .FALSE.
          IF(.NOT.root) EXIT
          buf(glons(1):glone(1),:) = x
        ELSE
          CALL p_recv (buf(glons(1):glone(1),:), &
                       gl_dc(i)% pe, tag_sum_zonal)  ! receve
        ENDIF
      ENDIF
    END DO
    IF (root) THEN
      !
      ! rotate longitudes
      !
      IF (dc%glon(jlat) /= dc%glon(1)) buf = CSHIFT (buf, shift=shift, dim=1)
      !
      ! sum
      !
      y = 0.0_dp
      DO l=1,SIZE(buf,2)
        DO k = 1, dc% nlon
          y = y + buf(k,l)
        END DO
      END DO
      !
      ! send result
      !
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_a == set_a) THEN
          IF (gl_dc(i)% pe /= pe) THEN
            CALL p_send (y, gl_dc(i)% pe, tag_sum_zonal)
          ENDIF
        ENDIF
      END DO
    ELSE
      !
      ! if not root send x and receive sum
      !
      CALL p_send (x(:dc%nglon,:), pe_root, tag_sum_zonal)
      CALL p_recv (y             , pe_root, tag_sum_zonal)
    ENDIF
  END FUNCTION sum_zonal_sl
!------------------------------------------------------------------------------
  FUNCTION sum_latit_sl (x) RESULT (y)
  !
  ! latitudinal sum: reduce last index, bounds: (any, nlat)
  !
  REAL(dp) ,INTENT (in) :: x(:,:)
  REAL(dp)              :: y (SIZE(x,1))
    INTEGER              :: i, k, set_b, pe, pe_root, glats(2), glate(2), nglat
    LOGICAL              :: first, root
    REAL(dp)                 :: buf (SIZE(x,1),dc%nlat)
    REAL(dp)    ,ALLOCATABLE :: tmp (:,:)
    set_b = dc% set_b  ! local set b
    pe    = dc% pe     ! local pe
    first =.TRUE.
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)% set_b == set_b) THEN
        !
        ! if root receive x
        !
        glats = gl_dc(i)% glats
        glate = gl_dc(i)% glate
        nglat = gl_dc(i)% nglat
        IF (first) THEN
          pe_root      = gl_dc(i)% pe  ! the first pe in the group is 'root'
          root         = (pe_root == pe)
          first        = .FALSE.
          IF(.NOT.root) EXIT
          buf(:,glats(1):glate(1)) = x(:,         :nglat/2)
          buf(:,glats(2):glate(2)) = x(:,nglat/2+1:       )
        ELSE
          ALLOCATE (tmp(SIZE(x,1),nglat))
          CALL p_recv (tmp, gl_dc(i)% pe, tag_sum_zonal)  ! receive
          buf(:,glats(1):glate(1)) = tmp(:,         :nglat/2)
          buf(:,glats(2):glate(2)) = tmp(:,nglat/2+1:       )
          DEALLOCATE (tmp)
        ENDIF
      ENDIF
    END DO
    IF (root) THEN
      !
      ! sum
      !
      y = 0.0_dp
      DO k = dc% nlat, 1, -1
        y = y + buf(:,k)
      END DO
      !
      ! send result
      !
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_b == set_b) THEN
          IF (gl_dc(i)% pe /= pe) THEN
            CALL p_send (y, gl_dc(i)% pe, tag_sum_zonal)
          ENDIF
        ENDIF
      END DO
    ELSE
      !
      ! if not root send x and receive sum
      !
      CALL p_send (x, pe_root, tag_sum_zonal)
      CALL p_recv (y, pe_root, tag_sum_zonal)
    ENDIF
  END FUNCTION sum_latit_sl
!==============================================================================
  FUNCTION sum_latit2 (x) RESULT (y)
  !
  ! latitudinal sum: reduce last index, bounds: (any, nlat)
  !
  REAL(dp) ,INTENT (in) :: x(:,:)
  REAL(dp)              :: y (SIZE(x,1))
    INTEGER              :: i, k, set_b, pe, pe_root, glats(2), glate(2), nglat
    LOGICAL              :: first, root
    REAL(dp)                 :: buf (SIZE(x,1),dc%nlat)
    REAL(dp)    ,ALLOCATABLE :: tmp (:,:)
    set_b = dc% set_b  ! local set b
    pe    = dc% pe     ! local pe
    first =.TRUE.
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)% set_b == set_b) THEN
        !
        ! if root receive x
        !
        glats = gl_dc(i)% glats
        glate = gl_dc(i)% glate
        nglat = gl_dc(i)% nglat
        IF (first) THEN
          pe_root      = gl_dc(i)% pe  ! the first pe in the group is 'root'
          root         = (pe_root == pe)
          first        = .FALSE.
          IF(.NOT.root) EXIT
          buf(:,glats(1):glate(1)) = x(:,         :nglat/2)
          buf(:,glats(2):glate(2)) = x(:,nglat/2+1:       )
        ELSE
          ALLOCATE (tmp(SIZE(x,1),nglat))
          CALL p_recv (tmp, gl_dc(i)% pe, tag_sum_zonal)  ! receive
          buf(:,glats(1):glate(1)) = tmp(:,         :nglat/2)
          buf(:,glats(2):glate(2)) = tmp(:,nglat/2+1:       )
          DEALLOCATE (tmp)
        ENDIF
      ENDIF
    END DO
    IF (root) THEN
      !
      ! sum
      !
      y = 0.0_dp
      DO k = 1, dc% nlat
        y = y + buf(:,k)
      END DO
      !
      ! send result
      !
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_b == set_b) THEN
          IF (gl_dc(i)% pe /= pe) THEN
            CALL p_send (y, gl_dc(i)% pe, tag_sum_zonal)
          ENDIF
        ENDIF
      END DO
    ELSE
      !
      ! if not root send x and receive sum
      !
      CALL p_send (x, pe_root, tag_sum_zonal)
      CALL p_recv (y, pe_root, tag_sum_zonal)
    ENDIF
  END FUNCTION sum_latit2
!------------------------------------------------------------------------------
  FUNCTION sum_latit1 (x) RESULT (y)
  !
  ! latitudinal sum: index, bound: (nlat)
  !
  REAL(dp) ,INTENT (in) :: x(:)
  REAL(dp)              :: y
    INTEGER              :: i, k, set_b, pe, pe_root, glats(2), glate(2), nglat
    LOGICAL              :: first, root
    REAL(dp)                 :: buf (dc%nlat)
    REAL(dp)    ,ALLOCATABLE :: tmp (:)
    set_b = dc% set_b  ! local set b
    pe    = dc% pe     ! local pe
    first =.TRUE.
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)% set_b == set_b) THEN
        !
        ! if root receive x
        !
        glats = gl_dc(i)% glats
        glate = gl_dc(i)% glate
        nglat = gl_dc(i)% nglat
        IF (first) THEN
          pe_root      = gl_dc(i)% pe  ! the first pe in the group is 'root'
          root         = (pe_root == pe)
          first        = .FALSE.
          IF(.NOT.root) EXIT
          buf(glats(1):glate(1)) = x(         :nglat/2)
          buf(glats(2):glate(2)) = x(nglat/2+1:       )
        ELSE
          ALLOCATE (tmp(nglat))
          CALL p_recv (tmp, gl_dc(i)% pe, tag_sum_zonal)  ! receve
          buf(glats(1):glate(1)) = tmp(         :nglat/2)
          buf(glats(2):glate(2)) = tmp(nglat/2+1:       )
          DEALLOCATE (tmp)
        ENDIF
      ENDIF
    END DO
    IF (root) THEN
      !
      ! sum
      !
      y = 0.0_dp
      DO k = 1, dc% nlat
        y = y + buf(k)
      END DO
      !
      ! send result
      !
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_b == set_b) THEN
          IF (gl_dc(i)% pe /= pe) THEN
            CALL p_send (y, gl_dc(i)% pe, tag_sum_zonal)
          ENDIF
        ENDIF
      END DO
    ELSE
      !
      ! if not root send x and receive sum
      !
      CALL p_send (x, pe_root, tag_sum_zonal)
      CALL p_recv (y, pe_root, tag_sum_zonal)
    ENDIF
  END FUNCTION sum_latit1
!==============================================================================
  FUNCTION sum_global3 (x) RESULT (y)
  !
  ! global sum (reduce index 1: lon, index 3 : lat)
  !
  REAL(dp) ,INTENT (in) :: x(:,:,:)
  REAL(dp)              :: y (SIZE(x,2))
   
    REAL(dp), ALLOCATABLE :: z (:,:,:)

    if (dc% lreg) then
      y = sum_latit (sum_zonal (x, 1))
    else
      allocate (z (dc% nglon, size(x,2), dc% nglat))
      call reorder (z,x)
      y = sum_latit (sum_zonal (z, 1))
      deallocate (z)
    endif

  END FUNCTION sum_global3
!------------------------------------------------------------------------------
  FUNCTION sum_global2 (x) RESULT (y)
  !
  ! global sum (reduce index 1: lon, index 2 : lat)
  !
  REAL(dp) ,INTENT (in) :: x(:,:)
  REAL(dp)              :: y

    REAL(dp), ALLOCATABLE :: z (:,:)

    if (dc% lreg) then    
      y = sum_latit (sum_zonal (x, 1))
    else
      allocate (z (dc% nglon, dc% nglat))
      call reorder (z,x)
      y = sum_latit (sum_zonal (z, 1))
      deallocate (z)
    endif
                           
  END FUNCTION sum_global2
!==============================================================================
END MODULE mo_global_op
