!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_test_trans
  !
  ! This module holds the routines to compare the decomposed fields
  ! on PE 1.. with the global fields calculated on PE 0 in debug mode
  !
  ! Authors:
  !
  ! A. Rhodin, MPI, August 1999, original source
  !

  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish
  USE mo_mpi,           ONLY: p_pe, p_io, p_send, p_recv, p_barrier
  USE mo_decomposition, ONLY: dcg => global_decomposition, &
                              dcl => local_decomposition,  &
                              debug_parallel, any_col_1d
  USE mo_transpose,     ONLY: gather_sp, gather_ls, gather_gp, &
                              gather_sa, tag_gather_gp
  USE mo_io_units,        ONLY: nerr
  USE mo_linked_list,   ONLY: t_stream
  USE mo_memory_base,   ONLY: new_stream,                          &
                              add_stream_element ,get_stream_element, &
                              remove_stream_element
  IMPLICIT NONE

  PRIVATE

  ! Generic routines to compare fields 

  PUBLIC :: test_spectral   ! in spectral  space (generic)
  PUBLIC :: test_legendre   ! in Legendre  space (generic)
  PUBLIC :: test_gridpoint  ! in gridpoint space (generic)
  PUBLIC :: test_symasym    ! in Legendre  space (generic, [a]sym Fourier comp)
  PUBLIC :: test_zonmean    ! zonal mean in gridpoint space
  PUBLIC :: test_scalar     ! compare scalar value
  PUBLIC :: test_row        ! compare latitude row
  PUBLIC :: test_slrow      ! compare latitude row (SLT conform index)

  ! Specific routines

  INTERFACE test_spectral
     MODULE PROCEDURE test_spectral2     ! (nlev, nsp for m=0)
     MODULE PROCEDURE test_spectral3     ! (nlev, nsp)
  END INTERFACE

  INTERFACE test_legendre
     MODULE PROCEDURE test_legendre0     ! (nlev,    nsp for m=0)
     MODULE PROCEDURE test_legendre3     ! (nlev, 2, nsp)
     MODULE PROCEDURE test_legendre4     ! (nmp1, nlevp1, nlat, nvar)
  END INTERFACE

  INTERFACE test_gridpoint
     MODULE PROCEDURE test_gridpoint2    ! (nlon,      nlat)
     MODULE PROCEDURE test_gridpoint3    ! (nlon, nlev,nlat)
     MODULE PROCEDURE test_gridpoint4
  END INTERFACE

  INTERFACE test_symasym
     MODULE PROCEDURE test_symasym2      ! (nlev[+1],          nhgl)
     MODULE PROCEDURE test_symasym4      ! (nlev[+1], 2, nmp1, nhgl)   
  END INTERFACE

  INTERFACE test_row
     ! real
     MODULE PROCEDURE test_row0          ! scalar
     MODULE PROCEDURE test_row1          ! (nlon[+x])
     MODULE PROCEDURE test_row2          ! (nlon[+x], nlev[+1])
     MODULE PROCEDURE test_row3          ! (nlon[+x], any, any)
     ! integer
     MODULE PROCEDURE itest_row0         ! scalar
     MODULE PROCEDURE itest_row1         ! (nlon[+x])
     MODULE PROCEDURE itest_row2         ! (nlon[+x], nlev[+1])
     MODULE PROCEDURE itest_row3         ! (nlon[+x], any, any)
     ! logical
     MODULE PROCEDURE ltest_row0         ! scalar
     MODULE PROCEDURE ltest_row1         ! (nlon[+x])
     MODULE PROCEDURE ltest_row2         ! (nlon[+x], nlev[+1])
     MODULE PROCEDURE ltest_row3         ! (nlon[+x], any, any)
  END INTERFACE

  INTERFACE test_slrow
     MODULE PROCEDURE test_slrow0        ! scalar
     MODULE PROCEDURE test_slrow1        ! (nlon[+x])
     MODULE PROCEDURE test_slrow2        ! (nlon[+x], nlev[+1])
     MODULE PROCEDURE test_slrow3        ! (nlon[+x], any, any)
  END INTERFACE

  TYPE (t_stream) ,POINTER :: p_test =>NULL()

CONTAINS

  !===========================================================================
  SUBROUTINE test_spectral3 (sp, name)
    REAL(dp)          ,INTENT(in) :: sp (:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL(dp) ,POINTER :: tmp (:,:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) ALLOCATE (tmp(SIZE(sp,1),SIZE(sp,2),SIZE(sp,3)))
       CALL gather_sp(tmp,sp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=sp)) THEN
             CALL finish('test_spectral3', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_spectral3
  !---------------------------------------------------------------------------
  SUBROUTINE test_spectral2 (sp, name)
    REAL(dp)          ,INTENT(in) :: sp (:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL(dp) ,POINTER :: tmp (:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) ALLOCATE (tmp(SIZE(sp,1),SIZE(sp,2)))
       CALL gather_sp(tmp,sp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=sp)) THEN
             CALL finish('test_spectral2', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_spectral2
  !===========================================================================
  SUBROUTINE test_legendre4 (ls, name)
    REAL(dp) ,INTENT(in)              :: ls (:,:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL(dp) ,POINTER :: tmp (:,:,:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) &
            ALLOCATE (tmp(SIZE(ls,1),SIZE(ls,2),SIZE(ls,3),SIZE(ls,4)))
       CALL gather_ls (tmp,ls,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=ls)) THEN
             CALL finish('test_legendre4', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_legendre4
  !---------------------------------------------------------------------------
  SUBROUTINE test_legendre3 (ls, name)
    REAL(dp) ,INTENT(in)              :: ls (:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL(dp) ,POINTER :: tmp (:,:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) ALLOCATE (tmp(SIZE(ls,1),SIZE(ls,2),SIZE(ls,3)))
       CALL gather_ls(tmp,ls,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=ls)) THEN
             CALL finish('test_legendre3','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_legendre3
  !---------------------------------------------------------------------------
  SUBROUTINE test_legendre0 (ls, name)
    REAL(dp)              ,INTENT(in) :: ls (:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL(dp) ,POINTER :: tmp (:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) ALLOCATE (tmp(SIZE(ls,1),SIZE(ls,2)))
       CALL gather_ls(tmp,ls,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=ls)) THEN
             CALL finish('test_legendre0', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_legendre0
  !===========================================================================
  SUBROUTINE test_symasym4 (sa, name)
    REAL(dp) ,INTENT(in)              :: sa (:,:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL(dp) ,POINTER :: tmp (:,:,:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) &
            ALLOCATE (tmp(SIZE(sa,1),SIZE(sa,2),SIZE(sa,3),SIZE(sa,4)))
       CALL gather_sa (tmp,sa,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=sa)) THEN
!print *,'tmp:shape,minval,maxval=',shape(tmp),minval(tmp),maxval(tmp)
!print *,'sa :shape,minval,maxval=',shape(sa ),minval(sa ),maxval(sa )
             CALL finish('test_symasym4', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_symasym4
  !---------------------------------------------------------------------------
  SUBROUTINE test_symasym2 (sa, name)
    REAL(dp) ,INTENT(in)              :: sa (:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL(dp) ,POINTER :: tmp (:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) &
            ALLOCATE (tmp(SIZE(sa,1),SIZE(sa,2)))
       CALL gather_sa (tmp,sa,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=sa)) THEN
             CALL finish('test_symasym4', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_symasym2
  !===========================================================================
  SUBROUTINE test_gridpoint4 (gp, name, abort)
  REAL(dp) ,INTENT(in)              :: gp (:,:,:,:)
  CHARACTER (len=*) ,INTENT(in) :: name
  LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL(dp) ,POINTER :: tmp (:,:,:,:)
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(gp,1),SIZE(gp,2),SIZE(gp,3),SIZE(gp,4)))
         tmp = gp
       ENDIF
       CALL gather_gp (tmp,gp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp(:dcl%nlon,:,:,:)/=gp(:dcl%nlon,:,:,:))) THEN
             WRITE(nerr,*)'test_gridpoint4: decomposition test failed for '&
                    //name
             IF(ab) CALL finish('test_gridpoint4', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_gridpoint4
  !----------------------------------------------------------------------------
  SUBROUTINE test_gridpoint3 (gp, name, abort)
  REAL(dp) ,INTENT(in)              :: gp (:,:,:)
  CHARACTER (len=*) ,INTENT(in) :: name
  LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL(dp) ,POINTER :: tmp (:,:,:)
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(gp,1),SIZE(gp,2),SIZE(gp,3)))
         tmp = gp
       ENDIF
       CALL gather_gp (tmp,gp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp(:dcl%nlon,:,:)/=gp(:dcl%nlon,:,:))) THEN
             WRITE(nerr,*)'test_gridpoint3: decomposition test failed for '&
                    //name
             IF(ab) CALL finish('test_gridpoint3', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_gridpoint3
  !----------------------------------------------------------------------------
  SUBROUTINE test_gridpoint2 (gp, name, abort)
  REAL(dp) ,INTENT(in)              :: gp (:,:)
  CHARACTER (len=*) ,INTENT(in) :: name
  LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL(dp) ,POINTER :: tmp (:,:)
    LOGICAL       :: ab
    IF (debug_parallel>=0) THEN
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(gp,1),SIZE(gp,2)))
         tmp = gp
       ENDIF
       CALL gather_gp(tmp,gp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp(:dcl%nlon,:)/=gp(:dcl%nlon,:))) THEN
             WRITE(nerr,*)'test_gridpoint2: decomposition test failed for '&
                    //name
             IF(ab)              CALL finish('test_gridpoint2', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
    CALL p_barrier
  END SUBROUTINE test_gridpoint2
  !===========================================================================
  SUBROUTINE test_zonmean (zm, name, abort)
    REAL(dp) ,INTENT(in)              :: zm (:,:) ! zonal mean (nlev,nlat)
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL(dp) ,POINTER :: tmp (:,:)
    INTEGER       :: i
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             ALLOCATE (tmp( SIZE(zm,1), dcg(i)%nglat))
             CALL p_recv (tmp, dcg(i)%pe, tag_gather_gp)
             IF (ANY(tmp(:,                :dcg(i)%nglh(1))   &
                  /=zm(:,dcg(i)%glats(1) :dcg(i)%glate(1)))) THEN
                  WRITE(nerr,*)'test_zonmean: decomposition test failed for '&
                    //name
                  IF (ab) CALL finish('test_zonmean',&
                    'decomposition test failed for '//name)
             ENDIF
             IF (ANY(tmp(:,dcg(i)%nglh(1)+1:)                 &
                  /=zm(:,dcg(i)%glats(2) :dcg(i)%glate(2)))) THEN
                  WRITE(nerr,*)'test_zonmean: decomposition test failed for '&
                    //name
                  IF (ab) CALL finish('test_zonmean',&
                    'decomposition test failed for '//name)
             ENDIF
             DEALLOCATE (tmp)
          END DO
       ELSE
          CALL p_send (zm, p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE test_zonmean
  !===========================================================================
  SUBROUTINE test_scalar (sc, name)
    REAL(dp) ,INTENT(in)              :: sc   ! scalar
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL(dp)    :: tmp
    INTEGER :: i
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             CALL p_recv (tmp, dcg(i)%pe, tag_gather_gp)
             IF (tmp/=sc) THEN
                WRITE (nerr,*)'test_scalar: decomposition test failed'
                WRITE (nerr,*)'test_scalar: PE,value=',p_io,sc
                WRITE (nerr,*)'test_scalar: PE,value=',dcg(i)%pe,tmp
                CALL finish('test_scalar','decomposition test failed for '//name)
             ENDIF
          END DO
       ELSE
          CALL p_send (sc, p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE test_scalar
  !===========================================================================
  !===========================================================================
  ! real tests
  SUBROUTINE test_row1 (rw, j, name, abort)
    REAL(dp)          ,INTENT(in) :: rw   (:) ! row (nlon[+x])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL(dp) ,POINTER :: rtmp (:), buf(:), z(:,:)
    INTEGER       :: i, k, jg
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
      IF(.NOT. ALL(dcg% lreg)) THEN
        !---------------------------------------------------------------
        ! nproma-blocking is enabled, fields are stored in stream p_test
        !---------------------------------------------------------------
        IF (.NOT.ASSOCIATED(p_test)) CALL construct_p_test
        IF (j==1) THEN
          CALL add_stream_element (p_test, name, z)
        ELSE
          CALL get_stream_element (p_test, name, z)
        ENDIF
        z (1:SIZE(rw,1),j) = rw(:)
        IF (j==dcl% ngpblks) THEN
          CALL test_gridpoint (z, name, abort)
          CALL remove_stream_element (p_test, name)
        ENDIF
      ELSE
        !-------------------------------------------------
        ! nproma-blocking is disabled, compare single rows
        !-------------------------------------------------
        IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE(rtmp(dcl%nlon))
          rtmp = rw
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon)) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   rtmp(dcg(i)%glons(1):dcg(i)%glone(1)) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      rtmp(dcg(i)%glons(2):dcg(i)%glone(2)) = buf
                   ELSE
                      rtmp(:dcg(i)%glone(2) )=buf( dcg(i)%nglon-dcg(i)%glone(2)+1:)
                      rtmp( dcg(i)%glons(2):)=buf(:dcg(i)%nglon-dcg(i)%glone(2))
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(rtmp/=rw(:dcl%nglon))) THEN
             DO k=1,dcl%nglon
                IF(rtmp(k)/=rw(k)) WRITE(nerr,*) k, rw(k), rtmp(k), rw(k)-rtmp(k)
             END DO
             IF (ab) &
                  CALL finish('test_row1','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(rtmp)
        ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (rw(:dcl%nglon), p_io, tag_gather_gp)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE test_row1
  !---------------------------------------------------------------------------
  SUBROUTINE test_row2 (rw, j, name, abort)
  REAL(dp)          ,INTENT(in) :: rw (:,:) ! row (nlon[+x], nlev[+1])
  INTEGER           ,INTENT(in) :: j        ! row index
  CHARACTER (len=*) ,INTENT(in) :: name
  LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL(dp) ,POINTER :: rtmp (:,:), buf(:,:), z(:,:,:)
    INTEGER       :: i, jg, k,l
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
      IF(.NOT. ALL(dcg% lreg)) THEN
        !---------------------------------------------------------------
        ! nproma-blocking is enabled, fields are stored in stream p_test
        !---------------------------------------------------------------
        IF (.NOT.ASSOCIATED(p_test)) CALL construct_p_test
        IF (j==1) THEN
          CALL add_stream_element (p_test, name, z, klev = SIZE(rw,2))
        ELSE
          CALL get_stream_element (p_test, name, z)
        ENDIF
        z (1:SIZE(rw,1),:,j) = rw(:,:)
        IF (j==dcl% ngpblks) THEN
          CALL test_gridpoint (z, name, abort)
          CALL remove_stream_element (p_test, name)
        ENDIF
      ELSE
        !-------------------------------------------------
        ! nproma-blocking is disabled, compare single rows
        !-------------------------------------------------
        IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE (rtmp (dcl%nlon, SIZE(rw,2)))
          rtmp = rw
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon, SIZE(rw,2))) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   rtmp(dcg(i)%glons(1):dcg(i)%glone(1),:) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      rtmp(dcg(i)%glons(2):dcg(i)%glone(2),:) = buf
                   ELSE
                      rtmp(:dcg(i)%glone(2) ,:)=buf( dcg(i)%nglon-dcg(i)%glone(2)+1:,:)
                      rtmp( dcg(i)%glons(2):,:)=buf(:dcg(i)%nglon-dcg(i)%glone(2),:)
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(rtmp/=rw(:dcl%nglon,:))) THEN
             WRITE(nerr,*)'test_row2: FAIL: jg=',jg,name
             DO l=1,SIZE(rtmp,2)
                DO k=1,SIZE(rtmp,1)
                   IF(rtmp(k,l)/=rw(k,l)) WRITE(nerr,*) &
                        'k,l, value(k,l) (1PE vs nPEs)',k,l,rw(k,l),rtmp(k,l)
                END DO
             END DO
             IF (ab) &
                  CALL finish('test_row2','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(rtmp)
        ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (rw(:dcl%nglon,:), p_io, tag_gather_gp)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE test_row2
  !---------------------------------------------------------------------------
  SUBROUTINE test_row3 (rw, j, name, abort)
    REAL(dp)          ,INTENT(in) :: rw (:,:,:) ! row (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL(dp) ,POINTER :: rtmp (:,:,:), buf(:,:,:), z(:,:,:,:)
    INTEGER       :: i, jg, k, l, m
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
      IF(.NOT. ALL(dcg% lreg)) THEN
        !---------------------------------------------------------------
        ! nproma-blocking is enabled, fields are stored in stream p_test
        !---------------------------------------------------------------
        IF (.NOT.ASSOCIATED(p_test)) CALL construct_p_test
        IF (j==1) THEN
          CALL add_stream_element (p_test, name, z, klev  = SIZE(rw,2),&
                                                    ktrac = SIZE(rw,3))
        ELSE
          CALL get_stream_element (p_test, name, z)
        ENDIF
        z (1:SIZE(rw,1),:,:,j) = rw(:,:,:)
        IF (j==dcl% ngpblks) THEN
          CALL test_gridpoint (z, name, abort)
          CALL remove_stream_element (p_test, name)
        ENDIF
      ELSE
        !-------------------------------------------------
        ! nproma-blocking is disabled, compare single rows
        !-------------------------------------------------
        IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE (rtmp (dcl%nlon, SIZE(rw,2), SIZE(rw,3)))
          rtmp = 0.0_dp
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon, SIZE(rw,2), SIZE(rw,3))) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   rtmp(dcg(i)%glons(1):dcg(i)%glone(1),:,:) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      rtmp(dcg(i)%glons(2):dcg(i)%glone(2),:,:) = buf
                   ELSE
                      rtmp(:dcg(i)%glone(2) ,:,:) = &
                           buf( dcg(i)%nglon-dcg(i)%glone(2)+1:,:,:)
                      rtmp( dcg(i)%glons(2):,:,:) = &
                           buf(:dcg(i)%nglon-dcg(i)%glone(2),:,:)
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(rtmp/=rw(:dcl%nglon,:,:))) THEN
             WRITE(nerr,*)'test_row3: FAIL: jg=',jg,name
             DO m=1,SIZE(rtmp,3)
                DO l=1,SIZE(rtmp,2)
                   DO k=1,SIZE(rtmp,1)
                      IF(rtmp(k,l,m)/=rw(k,l,m)) WRITE(nerr,*) &
                           'k,l,m, value(k,l,m) (1PE vs nPEs)',k,l,m,rw(k,l,m),rtmp(k,l,m)
                   END DO
                END DO
             END DO
             IF (ab) &
                  CALL finish('test_row3','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(rtmp)
        ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (rw(:dcl%nglon,:,:), p_io, tag_gather_gp)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE test_row3
  !---------------------------------------------------------------------------
  SUBROUTINE test_row0 (rw, jl, name, abort)
    REAL(dp)          ,INTENT(in) :: rw    ! scalar to compare
    INTEGER           ,INTENT(in) :: jl     ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL(dp)          :: rtmp
    INTEGER       :: i, jg
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF(.NOT. dcl% lreg) &
         CALL finish ('test_row','does not work with nproma enabled')
       IF (p_pe==p_io) THEN
          jg = dcl%glat(jl)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          rtmp = 0.0_dp
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive row0 only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                CALL p_recv (rtmp, dcg(i)%pe, tag_gather_gp)
                !--------
                ! compare
                !--------
                IF (rtmp/=rw) THEN
                   WRITE(nerr,*)'test_row0: FAIL: pe, jg=',dcg(i)%pe,jg,name
                   WRITE(nerr,*)'test_row0: pe, value',p_pe,rw
                   WRITE(nerr,*)'test_row0: pe, value',dcg(i)%pe,rtmp
                   IF (ab) &
                        CALL finish('test_row0','decomposition test failed for '//name)
                ENDIF
             ENDIF
          END DO
       ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (rw, p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE test_row0
  !===========================================================================
  !===========================================================================
  ! integer tests
  SUBROUTINE itest_row1 (iw, j, name, abort)
    INTEGER           ,INTENT(in) :: iw   (:) ! row (nlon[+x])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    INTEGER ,POINTER :: itmp (:), buf(:)
    INTEGER       :: i, k, jg
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
      IF(.NOT. ALL(dcg% lreg)) THEN
        CALL test_row1 (REAL(iw,dp), j, name, abort)
      ELSE
        IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE(itmp(dcl%nlon))
          itmp = iw
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon)) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   itmp(dcg(i)%glons(1):dcg(i)%glone(1)) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      itmp(dcg(i)%glons(2):dcg(i)%glone(2)) = buf
                   ELSE
                      itmp(:dcg(i)%glone(2) )=buf( dcg(i)%nglon-dcg(i)%glone(2)+1:)
                      itmp( dcg(i)%glons(2):)=buf(:dcg(i)%nglon-dcg(i)%glone(2))
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(itmp/=iw(:dcl%nglon))) THEN
             DO k=1,dcl%nglon
                IF(itmp(k)/=iw(k)) WRITE(nerr,*) k, iw(k), itmp(k)
             END DO
             IF (ab) &
                  CALL finish('test_row1','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(itmp)
        ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (iw(:dcl%nglon), p_io, tag_gather_gp)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE itest_row1
  !---------------------------------------------------------------------------
  SUBROUTINE itest_row2 (iw, j, name, abort)
    INTEGER           ,INTENT(in) :: iw (:,:) ! row (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    INTEGER ,POINTER :: itmp (:,:), buf(:,:)
    INTEGER       :: i, jg, k,l
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
      IF(.NOT. ALL(dcg% lreg)) THEN
        CALL test_row2 (REAL(iw,dp), j, name, abort)
      ELSE
        IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE (itmp (dcl%nlon, SIZE(iw,2)))
          itmp = iw
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon, SIZE(iw,2))) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   itmp(dcg(i)%glons(1):dcg(i)%glone(1),:) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      itmp(dcg(i)%glons(2):dcg(i)%glone(2),:) = buf
                   ELSE
                      itmp(:dcg(i)%glone(2) ,:)=buf( dcg(i)%nglon-dcg(i)%glone(2)+1:,:)
                      itmp( dcg(i)%glons(2):,:)=buf(:dcg(i)%nglon-dcg(i)%glone(2),:)
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(itmp/=iw(:dcl%nglon,:))) THEN
             WRITE(nerr,*)'test_row2: FAIL: jg=',jg,name
             DO l=1,SIZE(itmp,2)
                DO k=1,SIZE(itmp,1)
                   IF(itmp(k,l)/=iw(k,l)) WRITE(nerr,*) &
                        'k,l, value(k,l) (1PE vs nPEs)',k,l,iw(k,l),itmp(k,l)
                END DO
             END DO
             IF (ab) &
                  CALL finish('test_row2','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(itmp)
        ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (iw(:dcl%nglon,:), p_io, tag_gather_gp)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE itest_row2
  !---------------------------------------------------------------------------
  SUBROUTINE itest_row3 (iw, j, name, abort)
    INTEGER           ,INTENT(in) :: iw (:,:,:) ! row (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    INTEGER ,POINTER :: itmp (:,:,:), buf(:,:,:)
    INTEGER       :: i, jg, k, l, m
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
      IF(.NOT. ALL(dcg% lreg)) THEN
        CALL test_row3 (REAL(iw,dp), j, name, abort)
      ELSE
        IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE (itmp (dcl%nlon, SIZE(iw,2), SIZE(iw,3)))
          itmp = 0
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon, SIZE(iw,2), SIZE(iw,3))) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   itmp(dcg(i)%glons(1):dcg(i)%glone(1),:,:) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      itmp(dcg(i)%glons(2):dcg(i)%glone(2),:,:) = buf
                   ELSE
                      itmp(:dcg(i)%glone(2) ,:,:) = &
                           buf( dcg(i)%nglon-dcg(i)%glone(2)+1:,:,:)
                      itmp( dcg(i)%glons(2):,:,:) = &
                           buf(:dcg(i)%nglon-dcg(i)%glone(2),:,:)
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(itmp/=iw(:dcl%nglon,:,:))) THEN
             WRITE(nerr,*)'test_row3: FAIL: jg=',jg,name
             DO m=1,SIZE(itmp,3)
                DO l=1,SIZE(itmp,2)
                   DO k=1,SIZE(itmp,1)
                      IF(itmp(k,l,m)/=iw(k,l,m)) WRITE(nerr,*) &
                           'k,l,m, value(k,l,m) (1PE vs nPEs)',k,l,m,iw(k,l,m),itmp(k,l,m)
                   END DO
                END DO
             END DO
             IF (ab) &
                  CALL finish('test_row3','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(itmp)
        ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (iw(:dcl%nglon,:,:), p_io, tag_gather_gp)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE itest_row3
  !---------------------------------------------------------------------------
  SUBROUTINE itest_row0 (iw, jl, name, abort)
    INTEGER           ,INTENT(in) :: iw    ! scalar to compare
    INTEGER           ,INTENT(in) :: jl     ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    INTEGER       :: itmp
    INTEGER       :: i, jg
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF(.NOT. dcl% lreg) &
         CALL finish ('test_row','does not work with nproma enabled')
       IF (p_pe==p_io) THEN
          jg = dcl%glat(jl)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          itmp = 0
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive row0 only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                CALL p_recv (itmp, dcg(i)%pe, tag_gather_gp)
                !--------
                ! compare
                !--------
                IF (itmp/=iw) THEN
                   WRITE(nerr,*)'test_row0: FAIL: pe, jg=',dcg(i)%pe,jg,name
                   WRITE(nerr,*)'test_row0: pe, value',p_pe,iw
                   WRITE(nerr,*)'test_row0: pe, value',dcg(i)%pe,itmp
                   IF (ab) &
                        CALL finish('test_row0','decomposition test failed for '//name)
                ENDIF
             ENDIF
          END DO
       ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (iw, p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE itest_row0
  !===========================================================================
  !===========================================================================
  ! logical tests
  SUBROUTINE ltest_row1 (lw, j, name, abort)
    LOGICAL           ,INTENT(in) :: lw   (:) ! row (nlon[+x])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    LOGICAL ,POINTER :: ltmp (:), buf(:)
    REAL(dp)    ,ALLOCATABLE  :: z (:)
    INTEGER       :: i, k, jg
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
      IF(.NOT. ALL(dcg% lreg)) THEN
        ALLOCATE (z(SIZE(lw)))
        z = 0.0_dp; WHERE (lw) z = 1.0_dp
        CALL test_row1 (z, j, name, abort)
        DEALLOCATE (z)
      ELSE
        IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE(ltmp(dcl%nlon))
          ltmp = lw
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon)) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   ltmp(dcg(i)%glons(1):dcg(i)%glone(1)) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      ltmp(dcg(i)%glons(2):dcg(i)%glone(2)) = buf
                   ELSE
                      ltmp(:dcg(i)%glone(2) )=buf( dcg(i)%nglon-dcg(i)%glone(2)+1:)
                      ltmp( dcg(i)%glons(2):)=buf(:dcg(i)%nglon-dcg(i)%glone(2))
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(ltmp .NEQV. lw(:dcl%nglon))) THEN
             DO k=1,dcl%nglon
                IF(ltmp(k) .NEQV. lw(k)) WRITE(nerr,*) k, lw(k), ltmp(k)
             END DO
             IF (ab) &
                  CALL finish('test_row1','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(ltmp)
        ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (lw(:dcl%nglon), p_io, tag_gather_gp)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE ltest_row1
  !---------------------------------------------------------------------------
  SUBROUTINE ltest_row2 (lw, j, name, abort)
    LOGICAL           ,INTENT(in) :: lw (:,:) ! row (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    LOGICAL ,POINTER :: ltmp (:,:), buf(:,:)
    REAL(dp)    ,ALLOCATABLE  :: z (:,:)
    INTEGER       :: i, jg, k,l
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
      IF(.NOT. ALL(dcg% lreg)) THEN
        ALLOCATE (z(SIZE(lw,1),SIZE(lw,2)))
        z = 0.0_dp; WHERE (lw) z = 1.0_dp
        CALL test_row2 (z, j, name, abort)
        DEALLOCATE (z)
      ELSE
        IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE (ltmp (dcl%nlon, SIZE(lw,2)))
          ltmp = lw
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon, SIZE(lw,2))) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   ltmp(dcg(i)%glons(1):dcg(i)%glone(1),:) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      ltmp(dcg(i)%glons(2):dcg(i)%glone(2),:) = buf
                   ELSE
                      ltmp(:dcg(i)%glone(2) ,:)=buf( dcg(i)%nglon-dcg(i)%glone(2)+1:,:)
                      ltmp( dcg(i)%glons(2):,:)=buf(:dcg(i)%nglon-dcg(i)%glone(2),:)
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(ltmp .NEQV. lw(:dcl%nglon,:))) THEN
             WRITE(nerr,*)'test_row2: FAIL: jg=',jg,name
             DO l=1,SIZE(ltmp,2)
                DO k=1,SIZE(ltmp,1)
                   IF(ltmp(k,l) .NEQV. lw(k,l)) WRITE(nerr,*) &
                        'k,l, value(k,l) (1PE vs nPEs)',k,l,lw(k,l),ltmp(k,l)
                END DO
             END DO
             IF (ab) &
                  CALL finish('test_row2','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(ltmp)
        ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (lw(:dcl%nglon,:), p_io, tag_gather_gp)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE ltest_row2
  !---------------------------------------------------------------------------
  SUBROUTINE ltest_row3 (lw, j, name, abort)
    LOGICAL           ,INTENT(in) :: lw (:,:,:) ! row (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    LOGICAL ,POINTER :: ltmp (:,:,:), buf(:,:,:)
    REAL(dp)    ,ALLOCATABLE  :: z (:,:,:)
    INTEGER       :: i, jg, k, l, m
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
      IF(.NOT. ALL(dcg% lreg)) THEN
        ALLOCATE (z(SIZE(lw,1),SIZE(lw,2),SIZE(lw,3)))
        z = 0.0_dp; WHERE (lw) z = 1.0_dp
        CALL test_row3 (z, j, name, abort)
        DEALLOCATE (z)
      ELSE
        IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE (ltmp (dcl%nlon, SIZE(lw,2), SIZE(lw,3)))
          ltmp = .FALSE.
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon, SIZE(lw,2), SIZE(lw,3))) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   ltmp(dcg(i)%glons(1):dcg(i)%glone(1),:,:) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      ltmp(dcg(i)%glons(2):dcg(i)%glone(2),:,:) = buf
                   ELSE
                      ltmp(:dcg(i)%glone(2) ,:,:) = &
                           buf( dcg(i)%nglon-dcg(i)%glone(2)+1:,:,:)
                      ltmp( dcg(i)%glons(2):,:,:) = &
                           buf(:dcg(i)%nglon-dcg(i)%glone(2),:,:)
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(ltmp .NEQV. lw(:dcl%nglon,:,:))) THEN
             WRITE(nerr,*)'test_row3: FAIL: jg=',jg,name
             DO m=1,SIZE(ltmp,3)
                DO l=1,SIZE(ltmp,2)
                   DO k=1,SIZE(ltmp,1)
                      IF(ltmp(k,l,m) .NEQV. lw(k,l,m)) WRITE(nerr,*) &
                           'k,l,m, value(k,l,m) (1PE vs nPEs)',k,l,m,lw(k,l,m),ltmp(k,l,m)
                   END DO
                END DO
             END DO
             IF (ab) &
                  CALL finish('test_row3','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(ltmp)
        ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (lw(:dcl%nglon,:,:), p_io, tag_gather_gp)
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE ltest_row3
  !---------------------------------------------------------------------------
  SUBROUTINE ltest_row0 (lw, jl, name, abort)
    LOGICAL           ,INTENT(in) :: lw    ! scalar to compare
    INTEGER           ,INTENT(in) :: jl     ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    LOGICAL       :: ltmp
    INTEGER       :: i, jg
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF(.NOT. dcl% lreg) &
         CALL finish ('test_row','does not work with nproma enabled')
       IF (p_pe==p_io) THEN
          jg = dcl%glat(jl)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ltmp = .FALSE.
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive row0 only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                CALL p_recv (ltmp, dcg(i)%pe, tag_gather_gp)
                !--------
                ! compare
                !--------
                IF (ltmp .NEQV. lw) THEN
                   WRITE(nerr,*)'test_row0: FAIL: pe, jg=',dcg(i)%pe,jg,name
                   WRITE(nerr,*)'test_row0: pe, value',p_pe,lw
                   WRITE(nerr,*)'test_row0: pe, value',dcg(i)%pe,ltmp
                   IF (ab) &
                        CALL finish('test_row0','decomposition test failed for '//name)
                ENDIF
             ENDIF
          END DO
       ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (lw, p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE ltest_row0
  !============================================================================
  SUBROUTINE test_slrow3 (rw, jl, name, abort)
    REAL(dp)          ,INTENT(in) :: rw (:,:,:) ! slrow (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: jl         ! slt row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    CALL test_row3 (rw, dcl%nglat-jl+1, name, abort)
  END SUBROUTINE test_slrow3
  !---------------------------------------------------------------------------
  SUBROUTINE test_slrow2 (rw, jl, name, abort)
    REAL(dp)          ,INTENT(in) :: rw (:,:) ! row (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: jl       ! slt row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    CALL test_row2 (rw, dcl%nglat-jl+1, name, abort)
  END SUBROUTINE test_slrow2
  !---------------------------------------------------------------------------
  SUBROUTINE test_slrow1 (rw, jl, name, abort)
    REAL(dp)          ,INTENT(in) :: rw   (:) ! row (nlon[+x])
    INTEGER           ,INTENT(in) :: jl       ! slt row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    CALL test_row1 (rw, dcl%nglat-jl+1, name, abort)
  END SUBROUTINE test_slrow1
  !---------------------------------------------------------------------------
  SUBROUTINE test_slrow0 (rw, jl, name, abort)
    REAL(dp)          ,INTENT(in) :: rw     ! scalar to compare
    INTEGER           ,INTENT(in) :: jl     ! slt row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    CALL test_row0 (rw, dcl%nglat-jl+1, name, abort)
  END SUBROUTINE test_slrow0
  !============================================================================
  SUBROUTINE construct_p_test
    CALL new_stream (p_test ,'p_test',lrerun=.FALSE.,lpost=.FALSE.)
  END SUBROUTINE construct_p_test
  !============================================================================
END MODULE mo_test_trans
