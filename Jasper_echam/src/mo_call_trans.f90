#if defined(__SX__) || defined(ES)
#define FAST_AND_DIRTY 1
#endif

!#define MEASURE_LOAD_IMBALANCE

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_call_trans
  !
  ! This module holds the routines to invoke the transpositions
  ! with the fields to be transformed as actual parameters
  !
  ! Aditionally the following tests are made
  !
  !  In debug mode (PE 0 handles the whole domain, the other PE's
  !  handle the decomposed domain) the fields on PE0 are compared to those
  !  gathered from the other PE's .
  !
  ! Authors:
  !
  ! A. Rhodin, MPI, August 1999, original source
  !
  USE mo_decomposition, ONLY: gdc  => global_decomposition, &
                              debug_parallel, any_col_1d
  USE mo_exception,     ONLY: message
  USE mo_test_trans,    ONLY: test_spectral, test_legendre, test_gridpoint, &
                              test_symasym,  test_zonmean
  USE mo_control,       ONLY: ltimer, lyaxt_transposition
  USE mo_timer,         ONLY: timer_start, timer_stop, &
                              timer_s2l,   timer_l2s,  &
                              timer_l2f,   timer_f2l,  &
                              timer_f2g,   timer_g2f
#ifdef MEASURE_LOAD_IMBALANCE
  USE mo_real_timer,    ONLY: new_timer
  USE mo_mpi,           ONLY: p_barrier
#endif

  IMPLICIT NONE

  PRIVATE
  !
  ! inverse transpositions
  !
  PUBLIC :: spectral_to_legendre
  PUBLIC :: legendre_to_fourier
  PUBLIC :: fourier_to_gridpoint
  !
  ! direct transpositions
  !
  PUBLIC :: gridpoint_to_fourier
  PUBLIC :: fourier_to_legendre
  PUBLIC :: legendre_to_spectral
  !
  ! test routines
  !
  PUBLIC :: test_memory_f
  PUBLIC :: test_memory_gp
  PUBLIC :: test_scan_buffer

CONTAINS
  !============================================================================
  SUBROUTINE legendre_to_spectral
    USE mo_memory_ls, ONLY: ld, ltp, lvo, lu0 ! Legendre space
    USE mo_memory_sp, ONLY: sd, stp, svo, su0 ! spectral space
    USE mo_transpose, ONLY: tr_ls_sp          ! transposition routine
    USE mo_echam_yaxt,ONLY: yaxt_tr_ls_to_sp
#ifdef MEASURE_LOAD_IMBALANCE
    INTEGER, SAVE :: t_sync = 0
#endif
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
      CALL test_legendre (ld, 'ld')
      CALL test_legendre (ltp,'ltp')
      CALL test_legendre (lvo,'lvo')
      CALL test_legendre (lu0,'lu0')
      CALL message('','Test before legendre_to_spectral suceeded.')
    ENDIF
    !
    ! transposition: Legendre -> spectral
    !
#ifdef MEASURE_LOAD_IMBALANCE
    IF (ltimer) THEN
      IF (t_sync == 0) t_sync = new_timer('sync_ls2sp')
      CALL timer_start(t_sync)
      CALL p_barrier
      CALL timer_stop(t_sync)
    ENDIF
#endif
    IF (ltimer) CALL timer_start(timer_l2s)
    IF (lyaxt_transposition) THEN
      CALL yaxt_tr_ls_to_sp(ld, sd, ltp, stp, lvo, svo, lu0, su0)
    ELSE
      CALL tr_ls_sp (gdc, 1, ld, sd, ltp, stp, lvo, svo, lu0, su0)
    ENDIF
    IF (ltimer) CALL timer_stop(timer_l2s)
    !
    ! compare spectral fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
      CALL test_spectral (sd, 'sd')
      CALL test_spectral (stp,'stp')
      CALL test_spectral (svo,'svo')
      CALL test_spectral (su0,'su0')
      CALL message('','Test after legendre_to_spectral suceeded.')
    ENDIF
  END SUBROUTINE legendre_to_spectral
  !----------------------------------------------------------------------------
  SUBROUTINE spectral_to_legendre
    USE mo_memory_ls, ONLY: ld, ltp, lvo, lu0 ! Legendre space
    USE mo_memory_sp, ONLY: sd, stp, svo, su0 ! spectral space
    USE mo_transpose, ONLY: tr_ls_sp          ! transposition routine
    USE mo_echam_yaxt,ONLY: yaxt_tr_sp_to_ls
#ifdef MEASURE_LOAD_IMBALANCE
    INTEGER, SAVE :: t_sync = 0
#endif
    !
    ! compare spectral fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
      CALL test_spectral (sd, 'sd')
      CALL test_spectral (stp,'stp')
      CALL test_spectral (svo,'svo')
      CALL test_spectral (su0,'su0')
      CALL message('','Test before spectral_to_legendre suceeded.')
    ENDIF
    !
    ! transposition: spectral -> Legendre
    !
#ifdef MEASURE_LOAD_IMBALANCE
    IF (ltimer) THEN
      IF (t_sync == 0) t_sync = new_timer('sync_sp2ls')
      CALL timer_start(t_sync)
      CALL p_barrier
      CALL timer_stop(t_sync)
    ENDIF
#endif
    IF (ltimer) CALL timer_start(timer_s2l)
    IF (lyaxt_transposition) THEN
      CALL yaxt_tr_sp_to_ls(ld, sd, ltp, stp, lvo, svo, lu0, su0)
    ELSE
      CALL tr_ls_sp (gdc, -1, ld, sd, ltp, stp, lvo, svo, lu0, su0)
    ENDIF
    IF (ltimer) CALL timer_stop(timer_s2l)
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
      CALL test_legendre (ld, 'ld')
      CALL test_legendre (ltp,'ltp')
      CALL test_legendre (lvo,'lvo')
      CALL test_legendre (lu0,'lu0')
      CALL message('','Test after spectral_to_legendre suceeded.')
    ENDIF
  END SUBROUTINE spectral_to_legendre
  !============================================================================
  SUBROUTINE legendre_to_fourier
    USE mo_buffer_fft,  ONLY: fftz, fftl,& ! Fourier and Legendre space buffer
                              fbm0, lbm0   ! buffer for zonal means (m=0)
    USE mo_transpose,   ONLY: tr_fs_ls     ! transposition routine
!    USE mo_tr_alltoall, ONLY: transpose_fs_ls 
    USE mo_echam_yaxt,  ONLY: yaxt_tr_ls_to_fs
#ifdef MEASURE_LOAD_IMBALANCE
    INTEGER, SAVE :: t_sync = 0
#endif
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
      CALL test_legendre (fftl, 'fftl')
      CALL message('','Test before legendre_to_fourier suceeded.')
    ENDIF
    !
    ! transposition Fourier <- Legendre
    !
#ifdef MEASURE_LOAD_IMBALANCE
    IF (ltimer) THEN
      IF (t_sync == 0) t_sync = new_timer('sync_ls2fs')
      CALL timer_start(t_sync)
      CALL p_barrier
      CALL timer_stop(t_sync)
    ENDIF
#endif
    IF (ltimer) CALL timer_start(timer_l2f)
    IF (lyaxt_transposition) THEN
      CALL yaxt_tr_ls_to_fs(fftz, fftl, fbm0, lbm0)
    ELSE

!    IF (debug_parallel >= 0) THEN
      CALL tr_fs_ls (gdc, -1, fftz, fftl, fbm0, lbm0)
!    ELSE
!      CALL transpose_fs_ls (gdc, -1, fftz, fftl, fbm0, lbm0)
!    ENDIF

    ENDIF
    IF (ltimer) CALL timer_stop(timer_l2f)
  END SUBROUTINE legendre_to_fourier
  !----------------------------------------------------------------------------
  SUBROUTINE fourier_to_legendre
    USE mo_buffer_fft,  ONLY: fftz, fftl,& ! Fourier and Legendre space buffer
                              fbm0, lbm0   ! buffer for zonal means (m=0)
    USE mo_transpose,   ONLY: tr_fs_ls     ! transposition routine
!    USE mo_tr_alltoall, ONLY: transpose_fs_ls 
    USE mo_echam_yaxt,  ONLY: yaxt_tr_fs_to_ls
#ifdef MEASURE_LOAD_IMBALANCE
    INTEGER, SAVE :: t_sync = 0
#endif
    !
    ! transposition Fourier <- Legendre
    !
#ifdef MEASURE_LOAD_IMBALANCE
    IF (ltimer) THEN
      IF (t_sync == 0) t_sync = new_timer('sync_fs2ls')
      CALL timer_start(t_sync)
      CALL p_barrier
      CALL timer_stop(t_sync)
    ENDIF
#endif
    IF (ltimer) CALL timer_start(timer_f2l)
    IF (lyaxt_transposition) THEN
      CALL yaxt_tr_fs_to_ls(fftz(:,:,:,:6), fftl(:,:,:,:6), &
           &                fbm0(:,:,1:1),  lbm0(:,:,1:1))
    ELSE

!    IF (debug_parallel >= 0) THEN
      CALL tr_fs_ls (gdc, 1, fftz(:,:,:,:6), fftl(:,:,:,:6), &
                     fbm0(:,:,1:1),  lbm0(:,:,1:1))
!    ELSE
!      CALL transpose_fs_ls (gdc, 1, fftz(:,:,:,:6), fftl(:,:,:,:6), &
!                            fbm0(:,:,1:1),  lbm0(:,:,1:1))
!    ENDIF

    ENDIF
    IF (ltimer) CALL timer_stop(timer_f2l)
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
      CALL test_legendre (fftl(:,:,:,:), 'fftl')
      CALL message('','Test after fourier_to_legendre suceeded.')
    ENDIF
  END SUBROUTINE fourier_to_legendre
  !============================================================================
  SUBROUTINE fourier_to_gridpoint
    USE mo_buffer_fft,             ONLY: fftz, fbm0
    USE mo_scan_buffer,            ONLY: d, t, u, v, vo, dtm,  &
                                         dudl, dvdl, dtl,      &
                                         alps, dalpsl, dalpsm, &
                                         u0, du0, ul
    USE mo_transpose,              ONLY: tr_gp_fs
    USE mo_tr_alltoall,            ONLY: transpose_gp_fs
    USE mo_echam_yaxt,             ONLY: yaxt_tr_fs_to_gp
#ifdef MEASURE_LOAD_IMBALANCE
    INTEGER, SAVE :: t_sync = 0
#endif
    !
    ! transposition grid point <- Fourier
    !
#ifdef MEASURE_LOAD_IMBALANCE
    IF (ltimer) THEN
      IF (t_sync == 0) t_sync = new_timer('sync_fs2gp')
      CALL timer_start(t_sync)
      CALL p_barrier
      CALL timer_stop(t_sync)
    ENDIF
#endif
    IF (ltimer) CALL timer_start(timer_f2g)
      IF (debug_parallel >= 0) THEN
        CALL tr_gp_fs (gdc,-1,d,t,u,v,vo,dtm,dtl,        &
                       gp8=dudl, gp9=dvdl,               & 
                       sf1=dalpsl, sf2=dalpsm, sf3=alps, &
                       zm1=ul, zm2=u0, zm3=du0,          &
                       fs=fftz, fs0=fbm0)
      ELSE if (lyaxt_transposition) then
        CALL yaxt_tr_fs_to_gp(d,t,u,v,vo,dtm,dtl,        &
             &                gp8=dudl, gp9=dvdl,               &
             &                sf1=dalpsl, sf2=dalpsm, sf3=alps, &
             &                zm1=ul, zm2=u0, zm3=du0,          &
             &                fs=fftz, fs0=fbm0)
      ELSE
        CALL transpose_gp_fs (gdc,-1,d,t,u,v,vo,dtm,dtl,        &
                              gp8=dudl, gp9=dvdl,               & 
                              sf1=dalpsl, sf2=dalpsm, sf3=alps, &
                              zm1=ul, zm2=u0, zm3=du0,          &
                              fs=fftz, fs0=fbm0)
      ENDIF
    !
    IF (ltimer) CALL timer_stop(timer_f2g)
    !
    ! compare gridpoint fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
      IF(any_col_1d) THEN
        CALL message('','Test after fourier_to_gridpoint skipped')
        CALL message('','not tested: d t u v vo dtm')
        CALL message('','            dtl dalpsl dalpsm alps')
      ELSE
        CALL test_gridpoint (d,     'd')
        CALL test_gridpoint (t,     't')
        CALL test_gridpoint (u,     'u')
        CALL test_gridpoint (v,     'v')
        CALL test_gridpoint (vo,    'vo')
        CALL test_gridpoint (dtm,   'dtm')
        CALL test_gridpoint (dtl,   'dtl')
        CALL test_gridpoint (dalpsl,'dalpsl')
        CALL test_gridpoint (dalpsm,'dalpsm')
        CALL test_gridpoint (alps,  'alps')
        CALL message('','Test after fourier_to_gridpoint suceeded.')
      ENDIF
    ENDIF
  END SUBROUTINE fourier_to_gridpoint
  !------------------------------------------------------------------------------
  SUBROUTINE gridpoint_to_fourier
    USE mo_buffer_fft,             ONLY: fftz, fbm0
    USE mo_scan_buffer,            ONLY: rh, dm, vom, vol, u0, du0, ul
    USE mo_memory_g1a,             ONLY: alpsm1, dm1, tm1, vom1
    USE mo_transpose,              ONLY: tr_gp_fs
    USE mo_tr_alltoall,            ONLY: transpose_gp_fs
    USE mo_echam_yaxt,             ONLY: yaxt_tr_gp_to_fs
#ifdef MEASURE_LOAD_IMBALANCE
    INTEGER, SAVE :: t_sync = 0
#endif
    !
    ! transposition grid point -> Fourier
    !
#ifdef MEASURE_LOAD_IMBALANCE
    IF (ltimer) THEN
      IF (t_sync == 0) t_sync = new_timer('sync_gp2fs')
      CALL timer_start(t_sync)
      CALL p_barrier
      CALL timer_stop(t_sync)
    ENDIF
#endif
    IF (ltimer) CALL timer_start(timer_g2f)
    IF (debug_parallel >= 0) THEN
      CALL tr_gp_fs (gdc, 1,dm1,dm,tm1,rh,vol,vom,vom1, &
           sf3=alpsm1,                        &
           zm1=ul, zm2=u0, zm3=du0,           &
           fs=fftz, fs0=fbm0)
    ELSE IF(lyaxt_transposition) THEN
      CALL yaxt_tr_gp_to_fs(dm1,dm,tm1,rh,vol,vom,vom1,   &
           &                sf3=alpsm1,                   &
           &                zm1=ul, zm2=u0, zm3=du0,      &
           &                fs=fftz, fs0=fbm0)
    ELSE
      CALL transpose_gp_fs (gdc, 1,dm1,dm,tm1,rh,vol,vom,vom1, &
           sf3=alpsm1,                        &
           zm1=ul, zm2=u0, zm3=du0,           &
           fs=fftz, fs0=fbm0)
    ENDIF
    IF (ltimer) CALL timer_stop(timer_g2f)
  END SUBROUTINE gridpoint_to_fourier
  !==============================================================================
  SUBROUTINE test_scan_buffer (text)
    USE mo_scan_buffer,   ONLY: dtm, dtl, dalpsl, dalpsm, vo, d, t,   &
         alps, u, v, vol, vom, rh, qte, xlte,  &
         xite, tte, alpste, u0, du0, ul,       &
         alnpr, alpha, vervel
    CHARACTER (len=*) ,INTENT(in) :: text
    IF (debug_parallel>=0) THEN
      CALL test_gridpoint (dtm,   'dtm')  
      CALL test_gridpoint (dtl,   'dtl')
      CALL test_gridpoint (dalpsl,'dalpsl')
      CALL test_gridpoint (dalpsm,'dalpsm')
      CALL test_gridpoint (vo,    'vo')
      CALL test_gridpoint (d,     'd')
      CALL test_gridpoint (t,     't'    ,abort=.FALSE.)! column model
      CALL test_gridpoint (alps,  'alps' ,abort=.FALSE.)!   modifies
      CALL test_gridpoint (u,     'u'    ,abort=.FALSE.)!   these
      CALL test_gridpoint (v,     'v'    ,abort=.FALSE.)!   value
      CALL test_gridpoint (vol,   'vol')
      CALL test_gridpoint (vom,   'vom')
      CALL test_gridpoint (rh,    'rh')
      CALL test_gridpoint (qte,   'qte')
      CALL test_gridpoint (xlte,  'xlte')
      CALL test_gridpoint (xite,  'xite')
      CALL test_gridpoint (tte,   'tte')
      CALL test_gridpoint (alpste,'alpste')
      CALL test_zonmean   (u0,    'u0', abort=.FALSE.)! not
      CALL test_zonmean   (du0,   'du0',abort=.FALSE.)!(lon,[lev],lat)
      CALL test_zonmean   (ul,    'ul', abort=.FALSE.)! but (lev,lat)
      CALL test_gridpoint (alnpr, 'alnpr')
      CALL test_gridpoint (alpha, 'alpha')
      CALL test_gridpoint (vervel,'vervel')
      CALL message('','Test on scan_buffer suceeded '//text)
    ENDIF
  END SUBROUTINE test_scan_buffer
  !------------------------------------------------------------------------------
  SUBROUTINE test_memory_f (text)
    USE mo_memory_f,      ONLY: f
    USE mo_linked_list,   ONLY: list_element
    CHARACTER (len=*) ,INTENT(in) :: text
    TYPE (list_element) ,POINTER :: e
#ifdef FAST_AND_DIRTY
    IF (debug_parallel>=0) THEN
      CALL message('','Test on memory_f disabled!')
    ENDIF
#else
    IF (debug_parallel>=0) THEN
      e => f% first_list_element
      DO
        IF(.NOT.ASSOCIATED(e)) EXIT
        CALL test_symasym (e% field% ptr, e% field% info% name)
        e => e% next_list_element
      END DO
      CALL message('','Test on memory_f suceeded '//text)
    ENDIF
#endif    
  END SUBROUTINE test_memory_f
  !------------------------------------------------------------------------------
  SUBROUTINE test_memory_gp (gp, text)
    USE mo_linked_list,   ONLY: t_stream, list_element
    TYPE(t_stream)      ,INTENT(in) :: gp
    CHARACTER (len=*) ,INTENT(in) :: text
    TYPE (list_element) ,POINTER :: e
    IF (debug_parallel>=0) THEN
      e => gp% first_list_element
      DO
        IF(.NOT.ASSOCIATED(e)) EXIT
        CALL test_gridpoint (e% field% ptr(:,:,:,:), e% field% info% name)
        CALL message('','Test on memory_gp suceeded: '//e% field% info% name)
        e => e% next_list_element
      END DO
      CALL message('','Test on memory_gp suceeded '//text)
    ENDIF
  END SUBROUTINE test_memory_gp
  !==============================================================================
END MODULE mo_call_trans
