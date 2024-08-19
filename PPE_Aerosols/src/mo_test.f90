!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_test

  USE mo_kind,        ONLY: dp
  USE mo_jsbach_grid, ONLY: grid_type, domain_type
  USE mo_jsbach,      ONLY: debug, missing_value
  USE mo_linked_list, ONLY: t_stream, GAUSSIAN, NETCDF, GRIB
  USE mo_exception,   ONLY: message, message_text, finish
  USE mo_util_string, ONLY: int2string, real2string
  USE mo_memory_g3b,  ONLY: g3b

#if defined (__SX__) && defined (_OPENMP)
  USE omp_lib,          ONLY: omp_get_thread_num, &
                              omp_get_num_threads
#endif

  IMPLICIT NONE

  REAL(dp), POINTER, DIMENSION(:) :: test_hack_var1, test_hack_var2, test_hack_var3

  REAL(dp), POINTER, DIMENSION(:,:) :: &
       test0, test1, test2, test3, test4, test5, test6, test7, test8, test9, &
       test10, test11, test12, test13, test14, test15, test16, test17, test18, test19

  REAL(dp), POINTER, DIMENSION(:,:) :: &
       etest0, etest1, etest2, etest3, etest4, etest5, etest6, etest7, etest8, etest9, &
       etest10, etest11, etest12, etest13, etest14, etest15, etest16, etest17, etest18, etest19

  REAL(dp), POINTER, DIMENSION(:,:) :: &
       diff0, diff1, diff2, diff3, diff4, diff5, diff6, diff7, diff8, diff9, &
       diff10, diff11, diff12, diff13, diff14, diff15, diff16, diff17, diff18, diff19

  TYPE(t_stream), SAVE, POINTER :: IO_test
  LOGICAL, SAVE, POINTER :: domain_mask(:,:)
  INTEGER, SAVE :: domain_ndim, domain_nblocks, domain_nland
  REAL(dp), POINTER :: cover_fract(:,:)
  LOGICAL, POINTER :: cover_mask(:,:)
  INTEGER :: ntiles
  INTEGER, SAVE :: ilon0, ilon1, ilat0, ilat1
  INTEGER :: kidx0, kidx1, kblock
!$OMP THREADPRIVATE(kidx0, kidx1, kblock)

  INTERFACE test
     MODULE PROCEDURE test_0d_name
     MODULE PROCEDURE test_1d
     MODULE PROCEDURE test_1d_name
     MODULE PROCEDURE test_2d
     MODULE PROCEDURE test_2d_name
  END INTERFACE

  INTERFACE put_test_jvar
     MODULE PROCEDURE put_test_jvar_0d
  END INTERFACE

  INTERFACE put_test_var
     MODULE PROCEDURE put_test_var_1d
     MODULE PROCEDURE put_test_var_0d
  END INTERFACE

  INTERFACE put_test_var_nomask
     MODULE PROCEDURE put_test_var_nomask_1d
     MODULE PROCEDURE put_test_var_nomask_0d
  END INTERFACE

CONTAINS

  SUBROUTINE test_init_memory(grid, domain, stream)

    USE mo_netCDF,      ONLY : max_dim_name
    USE mo_memory_base, ONLY: new_stream, default_stream_setting, add =>add_stream_element
    USE mo_time_event,   ONLY: io_time_event, TIME_INC_MINUTES, TRIG_FIRST
    USE mo_time_control, ONLY: delta_time

    TYPE(grid_type), INTENT(in) :: grid
    TYPE(domain_type), INTENT(in) :: domain
    TYPE(t_stream), POINTER, OPTIONAL :: stream

    INTEGER  :: dim2p(2), dim2(2)
    CHARACTER(LEN=max_dim_name) :: dim2n(2)

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'test', filetype=NETCDF, lpost=.TRUE., lrerun=.FALSE., &
               interval = io_time_event(INT(delta_time/60),TIME_INC_MINUTES, TRIG_FIRST, 0))
          ! Set default stream options
          CALL default_stream_setting(stream, table  = 198, repr=GAUSSIAN, lpost=.TRUE., lrerun=.FALSE.)
       ENDIF
       IO_test => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_test, 'test', filetype=NETCDF, lpost=.TRUE., lrerun=.FALSE., &
            interval = io_time_event(INT(delta_time/60),TIME_INC_MINUTES, TRIG_FIRST, 0))
       ! Set default stream options
       CALL default_stream_setting(IO_test, table  = 198, repr=GAUSSIAN, lpost=.TRUE., lrerun=.FALSE.)
    ENDIF

    dim2p = (/ domain%ndim, domain%nblocks /)
    dim2  = (/ grid%nlon, grid%nlat /)
    dim2n(1) = 'lon'
    dim2n(2) = 'lat'

    CALL add(IO_test, 'test0', test0, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=1)
    CALL add(IO_test, 'test1', test1, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=2)
    CALL add(IO_test, 'test2', test2, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=3)
    CALL add(IO_test, 'test3', test3, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=4)
    CALL add(IO_test, 'test4', test4, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=5)
    CALL add(IO_test, 'test5', test5, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=6)
    CALL add(IO_test, 'test6', test6, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=7)
    CALL add(IO_test, 'test7', test7, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=8)
    CALL add(IO_test, 'test8', test8, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=9)
    CALL add(IO_test, 'test9', test9, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=10)
    CALL add(IO_test, 'test10', test10, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=11)
    CALL add(IO_test, 'test11', test11, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=12)
    CALL add(IO_test, 'test12', test12, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=13)
    CALL add(IO_test, 'test13', test13, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=14)
    CALL add(IO_test, 'test14', test14, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=15)
    CALL add(IO_test, 'test15', test15, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=16)
    CALL add(IO_test, 'test16', test16, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=17)
    CALL add(IO_test, 'test17', test17, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=18)
    CALL add(IO_test, 'test18', test18, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=19)
    CALL add(IO_test, 'test19', test19, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=20)

    CALL add(IO_test, 'etest0', etest0, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=51)
    CALL add(IO_test, 'etest1', etest1, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=52)
    CALL add(IO_test, 'etest2', etest2, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=53)
    CALL add(IO_test, 'etest3', etest3, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=54)
    CALL add(IO_test, 'etest4', etest4, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=55)
    CALL add(IO_test, 'etest5', etest5, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=56)
    CALL add(IO_test, 'etest6', etest6, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=57)
    CALL add(IO_test, 'etest7', etest7, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=58)
    CALL add(IO_test, 'etest8', etest8, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=59)
    CALL add(IO_test, 'etest9', etest9, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=60)
    CALL add(IO_test, 'etest10', etest10, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=61)
    CALL add(IO_test, 'etest11', etest11, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=62)
    CALL add(IO_test, 'etest12', etest12, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=63)
    CALL add(IO_test, 'etest13', etest13, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=64)
    CALL add(IO_test, 'etest14', etest14, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=65)
    CALL add(IO_test, 'etest15', etest15, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=66)
    CALL add(IO_test, 'etest16', etest16, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=67)
    CALL add(IO_test, 'etest17', etest17, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=68)
    CALL add(IO_test, 'etest18', etest18, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=69)
    CALL add(IO_test, 'etest19', etest19, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=70)

    CALL add(IO_test, 'diff0', diff0, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=101)
    CALL add(IO_test, 'diff1', diff1, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=102)
    CALL add(IO_test, 'diff2', diff2, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=103)
    CALL add(IO_test, 'diff3', diff3, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=104)
    CALL add(IO_test, 'diff4', diff4, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=105)
    CALL add(IO_test, 'diff5', diff5, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=106)
    CALL add(IO_test, 'diff6', diff6, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=107)
    CALL add(IO_test, 'diff7', diff7, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=108)
    CALL add(IO_test, 'diff8', diff8, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=109)
    CALL add(IO_test, 'diff9', diff9, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=110)
    CALL add(IO_test, 'diff10', diff10, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=111)
    CALL add(IO_test, 'diff11', diff11, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=112)
    CALL add(IO_test, 'diff12', diff12, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=113)
    CALL add(IO_test, 'diff13', diff13, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=114)
    CALL add(IO_test, 'diff14', diff14, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=115)
    CALL add(IO_test, 'diff15', diff15, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=116)
    CALL add(IO_test, 'diff16', diff16, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=117)
    CALL add(IO_test, 'diff17', diff17, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=118)
    CALL add(IO_test, 'diff18', diff18, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=119)
    CALL add(IO_test, 'diff19', diff19, dim2p, dim2, dimnames=dim2n,lmiss=.TRUE., missval=missing_value, &
         bits=64, code=120)

    test0   = missing_value
    test1   = missing_value
    test2   = missing_value
    test3   = missing_value
    test4   = missing_value
    test5   = missing_value
    test6   = missing_value
    test7   = missing_value
    test8   = missing_value
    test9   = missing_value
    test10  = missing_value
    test11  = missing_value
    test12  = missing_value
    test13  = missing_value
    test14  = missing_value
    test15  = missing_value
    test16  = missing_value
    test17  = missing_value
    test18  = missing_value
    test19  = missing_value
    etest0  = missing_value
    etest1  = missing_value
    etest2  = missing_value
    etest3  = missing_value
    etest4  = missing_value
    etest5  = missing_value
    etest6  = missing_value
    etest7  = missing_value
    etest8  = missing_value
    etest9  = missing_value
    etest10 = missing_value
    etest11 = missing_value
    etest12 = missing_value
    etest13 = missing_value
    etest14 = missing_value
    etest15 = missing_value
    etest16 = missing_value
    etest17 = missing_value
    etest18 = missing_value
    etest19 = missing_value
    diff0   = missing_value
    diff1   = missing_value
    diff2   = missing_value
    diff3   = missing_value
    diff4   = missing_value
    diff5   = missing_value
    diff6   = missing_value
    diff7   = missing_value
    diff8   = missing_value
    diff9   = missing_value
    diff10  = missing_value
    diff11  = missing_value
    diff12  = missing_value
    diff13  = missing_value
    diff14  = missing_value
    diff15  = missing_value
    diff16  = missing_value
    diff17  = missing_value
    diff18  = missing_value
    diff19  = missing_value

    ALLOCATE(test_hack_var1(domain%ndim))
    ALLOCATE(test_hack_var2(domain%ndim))
    ALLOCATE(test_hack_var3(domain%ndim))

  END SUBROUTINE test_init_memory

  SUBROUTINE init_test(grid, domain, itiles, IO_stream)

    USE mo_jsbach_grid, ONLY: grid_type, domain_type

    TYPE(grid_type), INTENT(in)          :: grid
    TYPE(domain_type), INTENT(in)        :: domain
    INTEGER, INTENT(in)                  :: itiles
    TYPE(t_stream), POINTER, OPTIONAL    :: IO_stream

    INTEGER :: i

    CALL test_init_memory(grid, domain, IO_stream)

    ALLOCATE(domain_mask(domain%ndim,domain%nblocks))
    domain_mask = domain%mask
    domain_ndim = domain%ndim
    domain_nblocks = domain%nblocks
    domain_nland = domain%nland
    ntiles = itiles
    ALLOCATE(cover_fract(grid%nland,ntiles))
    ALLOCATE(cover_mask(grid%nland,ntiles))

    DO i=1,grid%nlon
       IF (ABS(domain%lon(1)-grid%lon(i))<EPSILON(1._dp)) ilon0 = i
       IF (ABS(domain%lon(domain%nland)-grid%lon(i))<EPSILON(1._dp)) ilon1 = i
    END DO
    DO i=1,grid%nlat
       IF (ABS(domain%lat(1)-grid%lat(i))<EPSILON(1._dp)) ilat0 = i
       IF (ABS(domain%lat(domain%nland)-grid%lat(i))<EPSILON(1._dp)) ilat1 = i
    END DO

  END SUBROUTINE init_test

  SUBROUTINE init_test_2(i1,i2, iblock, fract, mask, domain)

    INTEGER, INTENT(in) :: i1, i2, iblock
    REAL(dp), INTENT(in)                     :: fract(:,:)
    LOGICAL, INTENT(in)                  :: mask(:,:)
    TYPE(domain_type), INTENT(in) :: domain

    kidx0 = i1
    kidx1 = i2
    kblock = iblock

    cover_fract(kidx0:kidx1,:) = fract
    cover_mask(kidx0:kidx1,:) = mask

    domain_mask(:,iblock) = domain%mask(:,iblock)
    domain_ndim = domain%ndim
    domain_nblocks = domain%nblocks
    domain_nland = domain%nland

  END SUBROUTINE init_test_2

  SUBROUTINE test_2d_name(jvar, index, ecode, text, mask)

    USE mo_utils, ONLY: average_tiles
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info

    REAL(dp), INTENT(in) :: jvar(:,:)
    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ecode
    LOGICAL, INTENT(in), OPTIONAL :: mask(:,:)
    
    REAL(dp) :: jvar_avg(SIZE(jvar,1))
    REAL(dp), POINTER :: evar(:,:), evar2(:,:)
    REAL(dp), POINTER :: testvar(:,:), diffvar(:,:)
    CHARACTER(LEN=8) :: cname

    IF (PRESENT(mask)) THEN
       CALL average_tiles(jvar, cover_mask(kidx0:kidx1,:).AND.mask, cover_fract(kidx0:kidx1,:), jvar_avg)
    ELSE
       CALL average_tiles(jvar, cover_mask(kidx0:kidx1,:), cover_fract(kidx0:kidx1,:), jvar_avg)
    END IF
    cname = 'test'//int2string(index)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF
    CALL get_stream_element(IO_test, cname, testvar)
    testvar(:,kblock) = UNPACK(jvar_avg, domain_mask(:,kblock), missing_value)

    IF (PRESENT(ecode)) THEN
       IF (LEN_TRIM(ecode)>5 .AND. ecode(1:MIN(5,LEN(ecode))) == 'etest') THEN
          CALL get_stream_element(IO_test, ecode, evar)
       ELSE
          CALL get_stream_element(g3b, ecode, evar)
          cname = 'etest'//int2string(index)
          CALL get_stream_element(IO_test, cname, evar2)
          evar2(:,kblock) = MERGE(evar(:,kblock), missing_value, domain_mask(:,kblock))
       END IF
       cname = 'diff'//int2string(index)
       CALL get_stream_element(IO_test, cname, diffvar)
       diffvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:,kblock), missing_value, domain_mask(:,kblock))
    END IF

    IF (debug .AND. kblock==domain_nblocks) THEN
       IF (PRESENT(ecode)) THEN
          PRINT*, text, ' (difference) ', MINVAL(diffvar), &
               MAXVAL(MERGE(diffvar,-missing_value,diffvar<0.95_dp*missing_value))
       ELSE
          PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       END IF
    END IF

  END SUBROUTINE test_2d_name

  SUBROUTINE test_2d(jvar, evar, index, text, mask)

    USE mo_utils, ONLY: average_tiles
    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info

    REAL(dp), INTENT(in) :: jvar(:,:)
    REAL(dp), INTENT(in) :: evar(:)
    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    LOGICAL, INTENT(in), OPTIONAL :: mask(:,:)
    
    REAL(dp) :: jvar_avg(SIZE(jvar,1))
    REAL(dp), POINTER :: testvar(:,:), evar2(:,:)
    CHARACTER(LEN=8) :: cname

    IF (PRESENT(mask)) THEN
       CALL average_tiles(jvar, cover_mask(kidx0:kidx1,:).AND.mask, cover_fract(kidx0:kidx1,:), jvar_avg)
    ELSE
       CALL average_tiles(jvar, cover_mask(kidx0:kidx1,:), cover_fract(kidx0:kidx1,:), jvar_avg)
    END IF
    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF
    testvar(:,kblock) = UNPACK(jvar_avg, domain_mask(:,kblock), missing_value)
    cname = 'etest'//int2string(index)
    CALL get_stream_element(IO_test, cname, evar2)
    evar2(:,kblock) = MERGE(evar, missing_value, domain_mask(:,kblock))

    testvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:), missing_value, domain_mask(:,kblock))

    IF (debug .AND. kblock==domain_nblocks) &
         PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))

  END SUBROUTINE test_2d

  SUBROUTINE test_1d_name(jvar, index, ecode, text)

    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info

    REAL(dp), INTENT(in) :: jvar(:)
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ecode
    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    
    REAL(dp), POINTER :: evar(:,:), evar2(:,:)
    REAL(dp), POINTER :: testvar(:,:), diffvar(:,:)
    CHARACTER(LEN=8) :: cname

#if defined (__SX__) && defined (_OPENMP)
    INTEGER :: tid
#endif

#if defined (__SX__) && defined (_OPENMP)
!!$    IF (debug) THEN
!!$       tid = omp_get_thread_num()
!!$       CALL message('test_1d_name', 'OpenMP thread #'//int2string(tid)//' '//int2string(kblock))
!!$    END IF
#endif

    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF
    testvar(:,kblock) = UNPACK(jvar, domain_mask(:,kblock), missing_value)


    IF (PRESENT(ecode)) THEN
       IF (LEN_TRIM(ecode)>5 .AND. ecode(1:MIN(5,LEN(ecode))) == 'etest') THEN
          CALL get_stream_element(IO_test, ecode, evar)
       ELSE
          CALL get_stream_element(g3b, ecode, evar)
          cname = 'etest'//int2string(index)
          CALL get_stream_element(IO_test, cname, evar2)
          evar2(:,kblock) = MERGE(evar(:,kblock), missing_value, domain_mask(:,kblock))
       END IF
       cname = 'diff'//int2string(index)
       CALL get_stream_element(IO_test, cname, diffvar)
       diffvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:,kblock), missing_value, domain_mask(:,kblock))
    END IF

    IF (debug .AND. kblock==domain_nblocks) THEN
       IF (PRESENT(ecode)) THEN
          PRINT*, text, ' (difference) ', MINVAL(diffvar), &
               MAXVAL(MERGE(diffvar,-missing_value,testvar<0.95_dp*missing_value))
       ELSE
          PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       END IF
    END IF

  END SUBROUTINE test_1d_name

  SUBROUTINE test_0d_name(index, ecode, text)

    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info

    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ecode
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    
    REAL(dp), POINTER :: evar(:,:), evar2(:,:)
    REAL(dp), POINTER :: testvar(:,:), diffvar(:,:)
    CHARACTER(LEN=8) :: cname

!    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF

    IF (PRESENT(ecode)) THEN
       IF (LEN_TRIM(ecode)>5 .AND. ecode(1:MIN(5,LEN(ecode))) == 'etest') THEN
          CALL get_stream_element(IO_test, ecode, evar)
       ELSE
          CALL get_stream_element(g3b, ecode, evar)
          cname = 'etest'//int2string(index)
          CALL get_stream_element(IO_test, cname, evar2)
          evar2(:,kblock) = MERGE(evar(:,kblock), missing_value, domain_mask(:,kblock))
       END IF
       cname = 'diff'//int2string(index)
       CALL get_stream_element(IO_test, cname, diffvar)
       diffvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:,kblock), missing_value, domain_mask(:,kblock))
    END IF

    IF (debug .AND. kblock==domain_nblocks) THEN
       IF (PRESENT(ecode)) THEN
          PRINT*, text, ' (difference) ', MINVAL(diffvar), &
               MAXVAL(MERGE(diffvar,-missing_value,testvar<0.95_dp*missing_value))
       ELSE
          PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       END IF
    END IF

  END SUBROUTINE test_0d_name

  SUBROUTINE put_test_jvar_0d(jvar, index, iland)

    USE mo_memory_base, ONLY: get_stream_element
 
    REAL(dp), INTENT(in) :: jvar
    INTEGER, INTENT(in) :: index, iland

    REAL(dp) :: var(domain_nland)
    REAL(dp), POINTER :: testvar(:,:)
    CHARACTER(LEN=8) :: cname

    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    var = PACK(testvar(:,kblock), domain_mask(:,kblock))
    var(iland) = jvar
    testvar(:,kblock) = UNPACK(var, domain_mask(:,kblock), missing_value)

  END SUBROUTINE put_test_jvar_0d

  SUBROUTINE test_1d(jvar, evar, index, text)

    USE mo_memory_base, ONLY: get_stream_element, set_stream_element_info

    REAL(dp), INTENT(in) :: jvar(:)
    REAL(dp), INTENT(in) :: evar(:)
    INTEGER, INTENT(in) :: index
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: text
    
    REAL(dp), POINTER :: testvar(:,:), evar2(:,:)
    CHARACTER(LEN=8) :: cname

    cname = 'test'//int2string(index)
    CALL get_stream_element(IO_test, cname, testvar)
    IF (PRESENT(text)) THEN
       CALL set_stream_element_info(IO_test, cname, longname=text)
    END IF
    testvar(:,kblock) = UNPACK(jvar, domain_mask(:,kblock), missing_value)
    cname = 'etest'//int2string(index)
    CALL get_stream_element(IO_test, cname, evar2)
    evar2(:,kblock) = MERGE(evar, missing_value, domain_mask(:,kblock))

    testvar(:,kblock) = MERGE(testvar(:,kblock) - evar(:), missing_value, domain_mask(:,kblock))

    IF (debug .AND. kblock==domain_nblocks) THEN
       IF (PRESENT(text)) THEN
          PRINT*, text, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       ELSE
          PRINT*, MINVAL(testvar), MAXVAL(MERGE(testvar,-missing_value,testvar<0.95_dp*missing_value))
       END IF
    END IF

   
  END SUBROUTINE test_1d


  SUBROUTINE put_test_var_1d(jrow, evar, index)

    USE mo_memory_base, ONLY: get_stream_element

    REAL(dp), INTENT(in) :: evar(:)
    INTEGER, INTENT(in) :: index, jrow
!    INTEGER, INTENT(in) :: iblock

    REAL(dp), POINTER :: testvar(:,:)

    CHARACTER(LEN=9) :: cname

    cname = 'etest'//int2string(index)

    CALL get_stream_element(IO_test, cname, testvar)

    testvar(:,jrow) = MERGE(evar, missing_value, domain_mask(:,jrow))
    

  END SUBROUTINE put_test_var_1d

  SUBROUTINE put_test_var_0d(jrow, evar, index, iproma)

    USE mo_memory_base, ONLY: get_stream_element

    REAL(dp), INTENT(in) :: evar
    INTEGER, INTENT(in) :: index
    INTEGER, INTENT(in) :: iproma, jrow

    REAL(dp), POINTER :: testvar(:,:)

    CHARACTER(LEN=9) :: cname

    cname = 'etest'//int2string(index)

    CALL get_stream_element(IO_test, cname, testvar)

    IF (domain_mask(iproma,jrow)) THEN
       testvar(iproma,jrow) = evar
    ELSE
       testvar(iproma,jrow) = missing_value
    END IF
    
  END SUBROUTINE put_test_var_0d

  SUBROUTINE put_test_var_nomask_1d(jrow, evar, index)

    USE mo_memory_base, ONLY: get_stream_element

    REAL(dp), INTENT(in) :: evar(:)
    INTEGER, INTENT(in) :: index, jrow
!    INTEGER, INTENT(in) :: iblock

    REAL(dp), POINTER :: testvar(:,:)
    INTEGER :: nproma

    CHARACTER(LEN=9) :: cname

    nproma = SIZE(evar)

    cname = 'etest'//int2string(index)

!!$    CALL message('put_test_var_nomask_1d',TRIM(cname)//' - '//TRIM(int2string(nproma))//', '//&
!!$         TRIM(int2string(jrow))//' - '//TRIM(real2string(MINVAL(evar(1:nproma))))//', '//&
!!$         TRIM(real2string(MAXVAL(evar(1:nproma)))))

    CALL get_stream_element(IO_test, cname, testvar)

    testvar(1:nproma,jrow) = evar(:)

  END SUBROUTINE put_test_var_nomask_1d

  SUBROUTINE put_test_var_nomask_0d(jrow, evar, index, iproma)

    USE mo_memory_base, ONLY: get_stream_element

    INTEGER, INTENT(in) :: jrow
    REAL(dp), INTENT(in) :: evar
    INTEGER, INTENT(in) :: index
    INTEGER, INTENT(in) :: iproma

    REAL(dp), POINTER :: testvar(:,:)

    CHARACTER(LEN=9) :: cname

    cname = 'etest'//int2string(index)

    CALL get_stream_element(IO_test, cname, testvar)

    testvar(iproma,jrow) = evar
    
  END SUBROUTINE put_test_var_nomask_0d
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE write_interface_variables(grid, domain, wind, wind10, temp_air, qair, &
       etAcoef, etBcoef, eqAcoef, eqBcoef, cdrag, precip_rain, precip_snow, CO2_concentration, lwdown, &
       sw_vis_net, sw_nir_net, sw_par_frac_diffuse, pressure, &
       echam_zchl, sw_par_down, czenith)
  !----------------------------------------------------------------------------

    USE mo_time_control, ONLY: lstart, lresume, next_date, get_date_components
    USE mo_filename,     ONLY: compose_filenames, standard_output_file
    USE mo_jsbach_grid,  ONLY: grid_type, domain_type
    USE mo_netcdf,       ONLY: nf_check

    TYPE(grid_type),    INTENT(in) :: grid
    TYPE(domain_type),  INTENT(in) :: domain

    REAL(dp), INTENT(in) :: wind(:)
    REAL(dp), INTENT(in) :: wind10(:)
    REAL(dp), INTENT(in) :: temp_air(:)
    REAL(dp), INTENT(in) :: qair(:)
    REAL(dp), INTENT(in) :: etAcoef(:)
    REAL(dp), INTENT(in) :: etBcoef(:)
    REAL(dp), INTENT(in) :: eqAcoef(:)
    REAL(dp), INTENT(in) :: eqBcoef(:)
    REAL(dp), INTENT(in) :: cdrag(:)
    REAL(dp), INTENT(in) :: precip_rain(:)
    REAL(dp), INTENT(in) :: precip_snow(:)
    REAL(dp), INTENT(in) :: CO2_concentration(:)
    REAL(dp), INTENT(in) :: lwdown(:)
    REAL(dp), INTENT(in) :: sw_vis_net(:)
    REAL(dp), INTENT(in) :: sw_nir_net(:)
    REAL(dp), INTENT(in) :: sw_par_frac_diffuse(:)
    REAL(dp), INTENT(in) :: pressure(:)
    REAL(dp), INTENT(in) :: echam_zchl(:)
    REAL(dp), INTENT(in), OPTIONAL :: sw_par_down(:)
    REAL(dp), INTENT(in), OPTIONAL :: czenith(:)

    INCLUDE 'netcdf.inc'

    !-- local variables

    INTEGER, SAVE :: ncid                   ! netcdf file id
    INTEGER :: lon_dim, lat_dim, time_dim   ! dimension IDs
    INTEGER :: lon_id, lat_id               ! dimension variable IDs
    INTEGER :: dim3(3)
    INTEGER :: tstart(1), tcount(1), fstart(3), fcount(3)

    CHARACTER(len=nf_max_name), SAVE :: filename
  
    INTEGER, SAVE :: icount = 0             ! time step counter

    INTEGER  :: year, month, day, hour, minute, second
    REAL(dp) :: date

    INTEGER, SAVE :: time_id
    INTEGER, SAVE :: wind_id
    INTEGER, SAVE :: wind10_id
    INTEGER, SAVE :: temp_air_id
    INTEGER, SAVE :: qair_id
    INTEGER, SAVE :: etAcoef_id
    INTEGER, SAVE :: etBcoef_id
    INTEGER, SAVE :: eqAcoef_id
    INTEGER, SAVE :: eqBcoef_id
    INTEGER, SAVE :: cdrag_id
    INTEGER, SAVE :: precip_rain_id
    INTEGER, SAVE :: precip_snow_id
    INTEGER, SAVE :: CO2_concentration_id
    INTEGER, SAVE :: lwdown_id
    INTEGER, SAVE :: sw_vis_net_id
    INTEGER, SAVE :: sw_nir_net_id
    INTEGER, SAVE :: sw_par_down_id
    INTEGER, SAVE :: sw_par_frac_diffuse_id
    INTEGER, SAVE :: czenith_id
    INTEGER, SAVE :: pressure_id
    INTEGER, SAVE :: echam_zchl_id


    !-- routine does not work in parallel runs nor with nproma smaller than nland

    IF (grid%nland > domain%nland) THEN
       CALL finish ('write_interface_variables','it is not possible to write the jsbach forcing variables in parallel runs')
    END IF

    !-- Create netcdf file

    IF (lstart .OR. lresume) THEN

       CALL compose_filenames
       filename = TRIM(standard_output_file)//'_interface.nc'

       ! enter define mode

       CALL nf_check(nf_create(TRIM(filename), NF_CLOBBER, ncid))

       WRITE (message_text,*) 'JSBACH forcing variables are written to: ', TRIM(filename)
       CALL message('write_interface_variables', message_text)


       ! define dimensions

       CALL nf_check(nf_def_dim(ncid, 'lon', grid%nlon, lon_dim))
       CALL nf_check(nf_def_dim(ncid, 'lat', grid%nlat, lat_dim))
       CALL nf_check(nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim))

       ! define variables

       CALL nf_check(nf_def_var(ncid, 'lon',  NF_DOUBLE, 1, (/lon_dim/),  lon_id))
       CALL nf_check(nf_def_var(ncid, 'lat',  NF_DOUBLE, 1, (/lat_dim/),  lat_id))
       CALL nf_check(nf_def_var(ncid, 'time', NF_DOUBLE, 1, (/time_dim/), time_id))

       dim3(1:3) = (/ lon_dim, lat_dim, time_dim /)
       CALL nf_check(nf_def_var(ncid, 'wind',        NF_DOUBLE, 3, dim3, wind_id))
       CALL nf_check(nf_def_var(ncid, 'wind10',      NF_DOUBLE, 3, dim3, wind10_id))
       CALL nf_check(nf_def_var(ncid, 'temp_air',    NF_DOUBLE, 3, dim3, temp_air_id))
       CALL nf_check(nf_def_var(ncid, 'qair',        NF_DOUBLE, 3, dim3, qair_id))
       CALL nf_check(nf_def_var(ncid, 'etAcoef',     NF_DOUBLE, 3, dim3, etAcoef_id))
       CALL nf_check(nf_def_var(ncid, 'etBcoef',     NF_DOUBLE, 3, dim3, etBcoef_id))
       CALL nf_check(nf_def_var(ncid, 'eqAcoef',     NF_DOUBLE, 3, dim3, eqAcoef_id))
       CALL nf_check(nf_def_var(ncid, 'eqBcoef',     NF_DOUBLE, 3, dim3, eqBcoef_id))
       CALL nf_check(nf_def_var(ncid, 'cdrag',       NF_DOUBLE, 3, dim3, cdrag_id))
       CALL nf_check(nf_def_var(ncid, 'precip_rain', NF_DOUBLE, 3, dim3, precip_rain_id))
       CALL nf_check(nf_def_var(ncid, 'precip_snow', NF_DOUBLE, 3, dim3, precip_snow_id))
       CALL nf_check(nf_def_var(ncid, 'CO2_concentration', NF_DOUBLE, 3, dim3, CO2_concentration_id))
       CALL nf_check(nf_def_var(ncid, 'lwdown',      NF_DOUBLE, 3, dim3, lwdown_id))
       CALL nf_check(nf_def_var(ncid, 'sw_vis_net',  NF_DOUBLE, 3, dim3, sw_vis_net_id))
       CALL nf_check(nf_def_var(ncid, 'sw_nir_net',  NF_DOUBLE, 3, dim3, sw_nir_net_id))
       CALL nf_check(nf_def_var(ncid, 'sw_par_frac_diffuse', NF_DOUBLE, 3, dim3, sw_par_frac_diffuse_id))
       CALL nf_check(nf_def_var(ncid, 'pressure',    NF_DOUBLE, 3, dim3, pressure_id))
       CALL nf_check(nf_def_var(ncid, 'echam_zchl',  NF_DOUBLE, 3, dim3, echam_zchl_id))
       CALL nf_check(nf_def_var(ncid, 'sw_par_down', NF_DOUBLE, 3, dim3, sw_par_down_id))
       CALL nf_check(nf_def_var(ncid, 'czenith',     NF_DOUBLE, 3, dim3, czenith_id))

       ! assign global attributes
       CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'conventions', 6, 'CF-1.0'))
       CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'source', 15, 'ECHAM6'))
       CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'institution', 36, 'Max-Planck-Institute for Meteorology'))

       ! leave define mode
       CALL nf_check(nf_enddef(ncid))

       CALL nf_check(nf_put_var_double(ncid, lat_id, grid%lat(:)))
       CALL nf_check(nf_put_var_double(ncid, lon_id, grid%lon(:)))

    ELSE

       CALL nf_check(nf_open(TRIM(filename), NF_WRITE, ncid))

    END IF

    !-- write data of the current time step

    icount = icount+1

    CALL get_date_components(next_date, year, month, day, hour, minute, second)
    date = ISIGN(1,year)*(IABS(year))*10000 + month*100 + day  &
             + (hour*3600 + minute*60 + second)/86400._dp

    tstart(1) = icount
    tcount(1) = 1
    CALL nf_check(nf_put_vara_double(ncid, time_id, tstart, tcount, date))

    fstart(1:3) = (/ 1, 1, icount /)
    fcount(1:3) = (/ grid%nlon, grid%nlat, 1 /)
    CALL write_var(ncid, wind_id,       fstart, fcount, wind,       grid)
    CALL write_var(ncid, wind10_id,     fstart, fcount, wind10,     grid)
    CALL write_var(ncid, temp_air_id,   fstart, fcount, temp_air,   grid)
    CALL write_var(ncid, qair_id,       fstart, fcount, qair,       grid)
    CALL write_var(ncid, etAcoef_id,    fstart, fcount, etAcoef,    grid)
    CALL write_var(ncid, etBcoef_id,    fstart, fcount, etBcoef,    grid)
    CALL write_var(ncid, eqAcoef_id,    fstart, fcount, eqAcoef,    grid)
    CALL write_var(ncid, eqBcoef_id,    fstart, fcount, eqBcoef,    grid)
    CALL write_var(ncid, cdrag_id,      fstart, fcount, cdrag,      grid)
    CALL write_var(ncid, precip_rain_id, fstart, fcount, precip_rain, grid)
    CALL write_var(ncid, precip_snow_id, fstart, fcount, precip_snow, grid)
    CALL write_var(ncid, CO2_concentration_id, fstart, fcount, CO2_concentration, grid)
    CALL write_var(ncid, lwdown_id,     fstart, fcount, lwdown,     grid)
    CALL write_var(ncid, sw_vis_net_id, fstart, fcount, sw_vis_net, grid)
    CALL write_var(ncid, sw_nir_net_id, fstart, fcount, sw_nir_net, grid)
    CALL write_var(ncid, sw_par_frac_diffuse_id, fstart, fcount, sw_par_frac_diffuse, grid)
    CALL write_var(ncid, pressure_id,   fstart, fcount, pressure,   grid)
    CALL write_var(ncid, echam_zchl_id, fstart, fcount, echam_zchl, grid)
    IF (PRESENT(sw_par_down)) CALL write_var(ncid, sw_par_down_id, fstart, fcount, sw_par_down, grid)
    IF (PRESENT(czenith))     CALL write_var(ncid, czenith_id,    fstart, fcount, czenith, grid)


    !-- close the file

    CALL nf_check(nf_close(ncid))

  END SUBROUTINE write_interface_variables
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE read_interface_variables(grid, domain, forcing, driving, cos_zenith)
  !----------------------------------------------------------------------------

    USE mo_time_control,      ONLY: lstart, lresume, l_putrerun, next_date, get_date_components
    USE mo_mpi,               ONLY: p_parallel_io
    USE mo_jsbach_grid,       ONLY: grid_type, domain_type
    USE mo_jsbalone_forcing,  ONLY: forcing_type
    USE mo_jsbalone,          ONLY: drive_type
    USE mo_netcdf,            ONLY: nf_check
    USE mo_jsbach_constants,  ONLY: Tmelt

    TYPE(grid_type),     INTENT(in)    :: grid
    TYPE(domain_type),   INTENT(in)    :: domain
    TYPE(forcing_type),  INTENT(inout) :: forcing
    TYPE(drive_type),    INTENT(inout) :: driving
    REAL(dp),            INTENT(out)   :: cos_zenith(:)

    INCLUDE 'netcdf.inc'

    !-- local variables

    INTEGER, SAVE :: ncid                         ! netcdf file id
    INTEGER, SAVE :: lon_dim, lat_dim, time_dim   ! dimensions
    INTEGER, SAVE :: lon_id, lat_id               ! dimension IDs

    CHARACTER(len=nf_max_name), SAVE :: filename

    INTEGER  :: ts 
    INTEGER  :: year, month, day, hour, minute, second
    REAL(dp) :: date
    REAL(dp), SAVE, ALLOCATABLE :: timevals(:)

    INTEGER, SAVE :: time_id
    INTEGER, SAVE :: wind_id
    INTEGER, SAVE :: wind10_id
    INTEGER, SAVE :: temp_air_id
    INTEGER, SAVE :: qair_id
    INTEGER, SAVE :: etAcoef_id
    INTEGER, SAVE :: etBcoef_id
    INTEGER, SAVE :: eqAcoef_id
    INTEGER, SAVE :: eqBcoef_id
    INTEGER, SAVE :: cdrag_id
    INTEGER, SAVE :: precip_rain_id
    INTEGER, SAVE :: precip_snow_id
    INTEGER, SAVE :: CO2_concentration_id
    INTEGER, SAVE :: lwdown_id
    INTEGER, SAVE :: sw_vis_net_id
    INTEGER, SAVE :: sw_nir_net_id
    INTEGER, SAVE :: sw_par_down_id
    INTEGER, SAVE :: sw_par_frac_diffuse_id
    INTEGER, SAVE :: czenith_id
    INTEGER, SAVE :: pressure_id
    INTEGER, SAVE :: echam_zchl_id


    IF (debug) WRITE (message_text,*) 'lstart: ', lstart, ', lresume: ', lresume 
    IF (debug) CALL message ('read_interface_variables', message_text)

    IF (p_parallel_io) THEN
       IF (lstart .OR. lresume) THEN

          ! open netcdf file

          filename = 'interface_variables.nc'
          CALL nf_check(nf_open(TRIM(filename), NF_NOWRITE, ncid))

          WRITE (message_text,*) 'JSBACH interface variables are read from: ', TRIM(filename)
          CALL message('read_interface_variables', message_text)

          ! check dimensions

          CALL nf_check(nf_inq_dimid(ncid, 'lon', lon_id))
          CALL nf_check(nf_inq_dimlen(ncid, lon_id, lon_dim))
          IF (lon_dim /= grid%nlon) THEN
             WRITE (message_text,*) 'longitude dimensions do not match: grid%nlon=', grid%nlon, ', nlon in file: ', lon_dim
             CALL finish ('read_interface_variables', message_text)
          END IF
          CALL nf_check(nf_inq_dimid(ncid, 'lat', lat_id))
          CALL nf_check(nf_inq_dimlen(ncid, lat_id, lat_dim))
          IF (lat_dim /= grid%nlat) THEN
             WRITE (message_text,*) 'latitude dimensions do not match: grid%nlat=', grid%nlat, ', nlat in file: ', lat_dim
             CALL finish ('read_interface_variables', message_text)
          END IF
          CALL nf_check(nf_inq_dimid(ncid, 'time', time_id))
          CALL nf_check(nf_inq_dimlen(ncid, time_id, time_dim))
          CALL nf_check(nf_inq_varid(ncid, 'time', time_id))
          ALLOCATE(timevals(time_dim))
          CALL nf_check(nf_get_var_double(ncid, time_id, timevals))

          ! get variable ids

          CALL nf_check(nf_inq_varid(ncid,  'wind',         wind_id))
          CALL nf_check(nf_inq_varid(ncid,  'wind10',       wind10_id))
          CALL nf_check(nf_inq_varid(ncid,  'temp_air',     temp_air_id))
          CALL nf_check(nf_inq_varid(ncid,  'qair',         qair_id))
          CALL nf_check(nf_inq_varid(ncid,  'etAcoef',      etAcoef_id))
          CALL nf_check(nf_inq_varid(ncid,  'etBcoef',      etBcoef_id))
          CALL nf_check(nf_inq_varid(ncid,  'eqAcoef',      eqAcoef_id))
          CALL nf_check(nf_inq_varid(ncid,  'eqBcoef',      eqBcoef_id))
          CALL nf_check(nf_inq_varid(ncid,  'cdrag',        cdrag_id))
          CALL nf_check(nf_inq_varid(ncid,  'precip_rain',  precip_rain_id))
          CALL nf_check(nf_inq_varid(ncid,  'precip_snow',  precip_snow_id))
          CALL nf_check(nf_inq_varid(ncid,  'CO2_concentration',  CO2_concentration_id))
          CALL nf_check(nf_inq_varid(ncid,  'lwdown',       lwdown_id))
          CALL nf_check(nf_inq_varid(ncid,  'sw_vis_net',   sw_vis_net_id))
          CALL nf_check(nf_inq_varid(ncid,  'sw_nir_net',   sw_nir_net_id))
          CALL nf_check(nf_inq_varid(ncid,  'sw_par_frac_diffuse', sw_par_frac_diffuse_id))
          CALL nf_check(nf_inq_varid(ncid,  'pressure',     pressure_id))
          CALL nf_check(nf_inq_varid(ncid,  'echam_zchl',   echam_zchl_id))
          CALL nf_check(nf_inq_varid(ncid,  'sw_par_down',  sw_par_down_id))
          CALL nf_check(nf_inq_varid(ncid,  'czenith',      czenith_id))
       END IF
    END IF


    !-- allocate forcing arrays
    NULLIFY(forcing%wind_speed);       ALLOCATE(forcing%wind_speed(domain%nland))
    NULLIFY(forcing%wind_speed10);     ALLOCATE(forcing%wind_speed10(domain%nland))
    NULLIFY(forcing%air_temp);         ALLOCATE(forcing%air_temp(domain%nland))
    NULLIFY(forcing%spec_humidity);    ALLOCATE(forcing%spec_humidity(domain%nland))
    NULLIFY(forcing%precip_rain);      ALLOCATE(forcing%precip_rain(domain%nland))
    NULLIFY(forcing%precip_snow);      ALLOCATE(forcing%precip_snow(domain%nland))
    NULLIFY(forcing%CO2_concentr);     ALLOCATE(forcing%CO2_concentr(domain%nland))
    NULLIFY(forcing%rad_lw_down);      ALLOCATE(forcing%rad_lw_down(domain%nland))
    NULLIFY(forcing%frac_PAR_diffuse); ALLOCATE(forcing%frac_PAR_diffuse(domain%nland))
    NULLIFY(forcing%air_pressure);     ALLOCATE(forcing%air_pressure(domain%nland))
    NULLIFY(forcing%rad_PAR_down);     ALLOCATE(forcing%rad_PAR_down(domain%nland))


    !-- get time step in file corresponding to current date
    IF (p_parallel_io) THEN
       CALL get_date_components(next_date, year, month, day, hour, minute, second)
       date = REAL(ISIGN(1,year)*(IABS(year)*10000+month*100+day), dp)  &
            + REAL(hour*3600 + minute*60 + second, dp)/86400._dp

       IF (debug) WRITE (message_text,*) 'date: ', year, '-', month, '-', day, ' time: ', hour, ':', minute, ':', second
       IF (debug) CALL message ('read_interface_variables', message_text)
       IF (debug) WRITE (message_text,*) 'date: ', date
       IF (debug) CALL message ('read_interface_variables', message_text)

       ts = 1
       DO WHILE (ts .LE. time_dim)
          IF ( date == timevals(ts) ) EXIT
          ts = ts + 1
       END DO
       IF (ts > time_dim) THEN 
          WRITE (message_text,*) 'date ', date, ' not available in ', filename 
          CALL finish('read_interface_variables', message_text)
       END IF
    ELSE
       ts = 0  ! not used in read_far
    END IF
 
    !-- read data of the current time step

    CALL read_var(grid, domain, ncid, wind_id,                ts, forcing%wind_speed)
    CALL read_var(grid, domain, ncid, wind10_id,              ts, forcing%wind_speed10)
    CALL read_var(grid, domain, ncid, temp_air_id,            ts, forcing%air_temp)
    forcing%air_temp = forcing%air_temp - Tmelt
    CALL read_var(grid, domain, ncid, qair_id,                ts, forcing%spec_humidity)
    CALL read_var(grid, domain, ncid, etAcoef_id,             ts, driving%etAcoef)
    CALL read_var(grid, domain, ncid, etBcoef_id,             ts, driving%etBcoef)
    CALL read_var(grid, domain, ncid, eqAcoef_id,             ts, driving%eqAcoef)
    CALL read_var(grid, domain, ncid, eqBcoef_id,             ts, driving%eqBcoef)
    CALL read_var(grid, domain, ncid, cdrag_id,               ts, driving%cdrag)
    CALL read_var(grid, domain, ncid, precip_rain_id,         ts, forcing%precip_rain)
    CALL read_var(grid, domain, ncid, precip_snow_id,         ts, forcing%precip_snow)
    CALL read_var(grid, domain, ncid, CO2_concentration_id,   ts, forcing%CO2_concentr)
    CALL read_var(grid, domain, ncid, lwdown_id,              ts, forcing%rad_lw_down)
    CALL read_var(grid, domain, ncid, sw_vis_net_id,          ts, driving%vis_net)
    CALL read_var(grid, domain, ncid, sw_nir_net_id,          ts, driving%nir_net)
    CALL read_var(grid, domain, ncid, sw_par_frac_diffuse_id, ts, forcing%frac_PAR_diffuse)
    CALL read_var(grid, domain, ncid, pressure_id,            ts, forcing%air_pressure)
    CALL read_var(grid, domain, ncid, echam_zchl_id,          ts, driving%zchl)
    CALL read_var(grid, domain, ncid, sw_par_down_id,         ts, forcing%rad_PAR_down)
    CALL read_var(grid, domain, ncid, czenith_id,             ts, cos_zenith)

    !-- close the file

    IF (p_parallel_io .AND. l_putrerun) DEALLOCATE (timevals)
    IF (p_parallel_io .AND. l_putrerun) CALL nf_check(nf_close(ncid))

  END SUBROUTINE read_interface_variables
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE write_var(ncid, varid, fstart, fcount, var, grid)
  !----------------------------------------------------------------------------

    USE mo_jsbach,   ONLY: missing_value
    USE mo_netcdf,   ONLY: nf_check

    INCLUDE 'netcdf.inc'

    INTEGER,  INTENT(in) :: ncid, varid
    INTEGER,  INTENT(in) :: fstart(3), fcount(3)
    REAL(dp), INTENT(in) :: var(:)
    TYPE(grid_type), INTENT(in) :: grid

    REAL(dp), DIMENSION(grid%nlon,grid%nlat) :: var_masked


    IF (SIZE(var) == grid%nland) THEN

       !-- unpack array (necessary in jsbach offline runs)

       var_masked = UNPACK(var, grid%mask, missing_value)

    ELSE

       !-- fill ocean grid cells with missing value

       WHERE (grid%mask(:,:))
          var_masked(:,:) = RESHAPE(var(:),(/grid%nlon,grid%nlat/))
       ELSEWHERE
          var_masked(:,:) = missing_value
       END WHERE
    END IF

    !-- write the variable

    CALL nf_check(nf_put_vara_double(ncid, varid, fstart, fcount, var_masked))

  END SUBROUTINE write_var

  !----------------------------------------------------------------------------
  SUBROUTINE read_var(grid, domain, ncid, varid, ts, var_local)
  !----------------------------------------------------------------------------

    USE mo_mpi,        ONLY: p_parallel_io
    USE mo_netcdf,     ONLY: nf_check
    USE mo_tr_scatter, ONLY: scatter_field
    USE mo_temp,       ONLY: zreal2d, zzreal2d, zreal2d_ptr

    INCLUDE 'netcdf.inc'

    TYPE(grid_type),   INTENT(in)    :: grid
    TYPE(domain_type), INTENT(in)    :: domain
    INTEGER,           INTENT(in)    :: ncid, varid
    INTEGER,           INTENT(in)    :: ts
    REAL(dp),          INTENT(out)   :: var_local(:)

    INTEGER            :: fstart(3), fcount(3)

    IF (debug) WRITE (message_text,*) 'reading variable id ', varid
    IF (debug) CALL message ('read_var', message_text)

    ALLOCATE(zreal2d(domain%ndim,domain%nblocks))

    IF (p_parallel_io) THEN
       ALLOCATE(zzreal2d(grid%nlon,grid%nlat))
       !-- read 2d variable

       fstart(1:3) = (/ 1, 1, ts /)
       fcount(1:3) = (/ grid%nlon, grid%nlat, 1 /)

       CALL nf_check(nf_get_vara_double(ncid, varid, fstart, fcount, zzreal2d))
    END IF
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zzreal2d(:,:)
    CALL scatter_field(zreal2d_ptr, zreal2d)
    var_local(:) = PACK(zreal2d, MASK=domain%mask)

    DEALLOCATE(zreal2d)
    IF (p_parallel_io) DEALLOCATE(zzreal2d)

  END SUBROUTINE read_var

END MODULE mo_test
