! Reshape is not that powerful implemented in some compiler
!
#if defined (__crayx1) || defined (sun) || defined (__SX__) || defined (ES) || defined (__PGI) || defined (__GFORTRAN__)
#define __REPLACE_RESHAPE 1
#endif
!
! Switch on explicit buffer packing and unpacking, allowing vectorization
!
#if defined (__crayx1) || defined (__SX__) || defined (ES) || defined(__PGI) || (defined __xlC__) || defined (_CRAYFTN)
#define __EXPLICIT 1
#endif
!
! Select communication type, if not defined NON NBLOCKING communication is selected
!
! do not use at all
#undef __FS_LS_ALLTOALLV

#ifdef NOMPI
#undef __FS_LS_ALLTOALLV
#endif

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_tr_alltoall
  !
  ! this module holds the transposition routines:
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
  USE mo_exception,     ONLY: message_text, message, finish
  USE mo_mpi,           ONLY: p_pe,             &! this processor
                              p_send,           &! MPI_send     routine
                              p_isend,          &! MPI_isend    routine
                              p_recv,           &! MPI_recv     routine
                              p_irecv,          &! MPI_irecv    routine
                              p_sendrecv,       &! MPI_sendrecv routine
                              p_wait,           &! MPI_wait     routine
                              p_real_dp, p_all_comm
  USE mo_decomposition, ONLY: pe_decomposed,    &! decomposition table data type
                              dc=>local_decomposition, &
                              gdc=>global_decomposition
  USE mo_buffer_fft,    ONLY: nvar_fsls, nvar_lsfs, nvar0_fsls, nvar0_lsfs

  USE mo_tr_omp_decomposition, ONLY: init_omp_decomposition,         &
                                     omp_get_thread_num,             &
                                     max_nthreads, nthreads,         &
                                     thsp_glat, thsp_lat, thsp_blk,  &
                                     thsp_flat
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ! public routines:
  !
  ! transpositions:
  !
  PUBLIC :: transpose_gp_fs   ! transpose grid point space <-> Fourier  space
  PUBLIC :: transpose_fs_ls   ! transpose Fourier    space <-> Legendre space
  PUBLIC :: transpose_ls_sp   ! transpose Legendre   space <-> spectral space
  PUBLIC :: transpose_gp_ffsl ! transpose grid point space <-> ffsl decomposition
  !
  ! reorder arrays in gridpoint space 
  !
  PUBLIC :: reorder
  !
  ! interfaces (specific routines)
  !
  ! flux-form semi-lagrangian transport scheme decomposition
  !
  INTERFACE transpose_gp_ffsl
    MODULE PROCEDURE transpose_gp_ffsl_2
    MODULE PROCEDURE transpose_gp_ffsl_3
    MODULE PROCEDURE transpose_gp_ffsl_4
  END INTERFACE
  !
  ! reorder array elements for blocksize = nproma
  !
  INTERFACE reorder
    MODULE PROCEDURE reorder12
    MODULE PROCEDURE reorder21
    MODULE PROCEDURE reorder23
    MODULE PROCEDURE reorder32
    MODULE PROCEDURE reorder2
    MODULE PROCEDURE reorder3
    MODULE PROCEDURE reorder4

    MODULE PROCEDURE subreorder2
    MODULE PROCEDURE subreorder3
  END INTERFACE
  !
  INTERFACE indx
    MODULE PROCEDURE indx0
    MODULE PROCEDURE indx2
  END INTERFACE
  !
  ! define tags
  !
  INTEGER, PARAMETER :: tag_tr_gp_fs     = 200
  INTEGER, PARAMETER :: tag_tr_fs_ls     = 210
  INTEGER, PARAMETER :: tag_tr_ls_sp     = 220
  INTEGER, PARAMETER :: tag_tr_gp_ffsl   = 230
  !
  ! define communication types
  ! 
  INTEGER, PARAMETER :: sendrecv     = 0
  INTEGER, PARAMETER :: non_blocking = 1
  INTEGER, PARAMETER :: onesided     = 2
  INTEGER, PARAMETER :: alltoall     = 3
  !
  !====================================================================== 
  !
  ! transpose buffer
  !
  ! 3 dims for gp-fs and fs-gp, and fs-ls and ls-fs
  ! 2 dims for ls-sp and sp-ls (set last index to 1 for allocate)
  !
  TYPE transpose_buffer
    SEQUENCE
    REAL(dp), ALLOCATABLE :: send_buffer(:,:,:,:)
    REAL(dp), ALLOCATABLE :: recv_buffer(:,:,:,:)
    REAL(dp), ALLOCATABLE :: send_buffer0(:,:,:)
    REAL(dp), ALLOCATABLE :: recv_buffer0(:,:,:)
  END TYPE transpose_buffer
  !
  TYPE(transpose_buffer), ALLOCATABLE :: fs_ls(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: ls_fs(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: gp_fs(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: fs_gp(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: ls_sp(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: sp_ls(:) 
  !
  INTEGER, ALLOCATABLE :: plan_a(:)
  INTEGER, ALLOCATABLE :: plan_b(:)
  !
  LOGICAL, SAVE :: require_init = .TRUE.
  ! 
  INTEGER, PARAMETER, PUBLIC :: foreward =  1
  INTEGER, PARAMETER, PUBLIC :: backward = -1
  !
  !==============================================================================
CONTAINS
  !==============================================================================
  !
  SUBROUTINE init_transpose

    IF (require_init) THEN
      CALL init_omp_decomposition
      require_init = .FALSE.
    ENDIF

  END SUBROUTINE init_transpose
  !
  !==============================================================================
  !
  SUBROUTINE transpose_gp_fs (gl_dc, sign, gp1, gp2, gp3, gp4, gp5, gp6, gp7, gp8, gp9, &
                                           sf1, sf2, sf3, zm1, zm2, zm3, fs, fs0)
    !
    ! transpose:
    !   sign=foreward : grid point space to Fourier space
    !   sign=backward : grid point space from  Fourier space
    !
    TYPE (pe_decomposed) ,INTENT(in)     :: gl_dc  (:)       ! decomposition
    INTEGER              ,INTENT(in)     :: sign             ! foreward:gp>fs; backward:gp<fs
    REAL(dp)             ,INTENT(inout)  :: gp1    (:,:,:)   ! gridpoint space 3d
    REAL(dp)             ,INTENT(inout)  :: gp2    (:,:,:)   !
    REAL(dp)             ,INTENT(inout)  :: gp3    (:,:,:)   !
    REAL(dp)             ,INTENT(inout)  :: gp4    (:,:,:)   !
    REAL(dp)             ,INTENT(inout)  :: gp5    (:,:,:)   !
    REAL(dp)             ,INTENT(inout)  :: gp6    (:,:,:)   !
    REAL(dp)             ,INTENT(inout)  :: gp7    (:,:,:)   !
    REAL(dp) ,OPTIONAL   ,INTENT(inout)  :: gp8    (:,:,:)   ! for u wind deriv.
    REAL(dp) ,OPTIONAL   ,INTENT(inout)  :: gp9    (:,:,:)   ! for v wind deriv.
    REAL(dp) ,OPTIONAL   ,INTENT(inout)  :: sf1    (:,:)     ! gridpoint space 2d
    REAL(dp) ,OPTIONAL   ,INTENT(inout)  :: sf2    (:,:)     ! gridpoint space 2d
    REAL(dp) ,OPTIONAL   ,INTENT(inout)  :: sf3    (:,:)     ! gridpoint space 2d
    REAL(dp) ,OPTIONAL   ,INTENT(inout)  :: zm1    (:,:)     ! zonal mean
    REAL(dp) ,OPTIONAL   ,INTENT(inout)  :: zm2    (:,:)     ! zonal mean
    REAL(dp) ,OPTIONAL   ,INTENT(inout)  :: zm3    (:,:)     ! zonal mean
    REAL(dp)             ,INTENT(inout)  :: fs     (:,:,:,:) ! Fourier space
    REAL(dp) ,OPTIONAL   ,INTENT(inout)  :: fs0    (:,:,:)   ! zonal mean, Fourier space
    !
    ! Data structures:
    !
    ! local variables
    !
    INTEGER :: i, k, n
    INTEGER :: imype         ! index of this pe
    INTEGER :: nprocb        ! number of PEs in set A
    INTEGER :: ks, ke, nk    ! vertical range of buffer
    INTEGER :: nk0           ! vertical range of buffer bu0
    INTEGER :: nglat, nglon  ! gridpoint space no. lons,lats
    INTEGER :: glons(2)      ! first longitudes in gridspace
    INTEGER :: glone(2)      ! last  longitudes in gridspace
    INTEGER :: nglh(2)       ! number of lats in each domains
    INTEGER :: nlon          ! global number of longitudes
    INTEGER :: nvar, nvarmax ! number of variables in fft buffer
    INTEGER :: nva0          ! number of variables (zonal mean)
    LOGICAL :: lreg          ! regular lon/lat ordering

    INTEGER, ALLOCATABLE, SAVE :: idest(:) ! destination of data
    INTEGER, ALLOCATABLE, SAVE :: isrc(:)  ! source of data

    LOGICAL, ALLOCATABLE, SAVE :: lm0r(:)  ! receive m0
    LOGICAL, ALLOCATABLE, SAVE :: lm0s(:)  ! send m0

    LOGICAL, SAVE :: first_call = .TRUE.

    INTEGER, SAVE :: ilats(0:max_nthreads-1),ilate(0:max_nthreads-1)
    INTEGER, SAVE :: iblks(0:max_nthreads-1),iblke(0:max_nthreads-1)

    INTEGER :: it, il, iv, ils, ile, ibs, ibe

    IF (require_init) THEN
      CALL init_transpose
    ENDIF

    IF (first_call) THEN
       first_call=.FALSE.

       imype  = indx (p_pe, gl_dc)

       DO i = 0, nthreads-1
         ilats(i) = thsp_glat(i)%begin
         ilate(i) = thsp_glat(i)%end
         iblks(i) = thsp_blk(i)%begin
         iblke(i) = thsp_blk(i)%end
       ENDDO

    ENDIF

    IF (sign == foreward) nvar = SIZE (fs,4)-2
    IF (sign == backward) nvar = SIZE (fs,4)
    nvarmax =  SIZE (fs,4)
    imype  = indx (p_pe, gl_dc)
    nprocb = gl_dc(imype)%nprocb
    lreg   = gl_dc(imype)%lreg
    IF (gl_dc(imype)%col_1d) RETURN
    !
    IF (.NOT. ALLOCATED(plan_b)) THEN
      ALLOCATE(plan_b(0:nprocb-1))
    ENDIF

    k = 0
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)%set_a /=  gl_dc(imype)%set_a) CYCLE
      plan_b(k) = i ! = gl_dc(i)%pe
      IF (i == imype) n = k
      k = k + 1
    END DO
    plan_b = CSHIFT (plan_b,n)

    nva0 = 0 
    IF (PRESENT(fs0)) nva0 = SIZE (fs0,3)
    
    IF (.NOT. ALLOCATED(idest)) THEN
      ALLOCATE(idest(0:nprocb-1))
      ALLOCATE(isrc(0:nprocb-1))
    ENDIF
      
    DO k = 0, nprocb-1
      idest(k) = plan_b(           k        ) ! PE index (send)
      isrc(k)  = plan_b(MOD(nprocb-k,nprocb)) ! PE index (recv)
    ENDDO
    
    IF (.NOT. ALLOCATED(lm0r)) THEN
       ALLOCATE(lm0r(0:nprocb-1))
       ALLOCATE(lm0s(0:nprocb-1))
     ENDIF
    !
    !------------------------------------------------------------------------
    !
    IF (.NOT. ALLOCATED(gp_fs)) THEN
      ALLOCATE (gp_fs(0:nprocb-1))
      
      DO k = 0, nprocb-1
        
        ! GridPoint -> Fourier Space
        
        nk    = gl_dc(idest(k))%nflevp1
        nk0   = gl_dc(idest(k))%nflev
        nglat = gl_dc(imype)%nglat
        nglon = gl_dc(imype)%nglon
        
        ALLOCATE (gp_fs(k)%send_buffer(nglon, nk, nglat, nvar) )
        ALLOCATE (gp_fs(k)%send_buffer0(nk0, nglat, nva0))

!$OMP PARALLEL PRIVATE(it,ils,ile)
        it = omp_get_thread_num()
        ils = ilats(it)
        ile = ilate(it)
        gp_fs(k)%send_buffer(:,:,ils:ile,:) = 0.0_dp
        gp_fs(k)%send_buffer0(:,ils:ile,:)  = 0.0_dp
!$OMP END PARALLEL

        nk    = gl_dc(imype)%nflevp1
        nk0   = gl_dc(imype)%nflev
        nglat = gl_dc(isrc(k))%nglat
        nglon = gl_dc(isrc(k))%nglon
        
        ALLOCATE (gp_fs(k)%recv_buffer(nglon, nk, nglat, nvar) )
        ALLOCATE (gp_fs(k)%recv_buffer0(nk0, nglat, nva0))
        
      ENDDO
      
    ENDIF
    
    IF (.NOT. ALLOCATED(fs_gp)) THEN
      ALLOCATE (fs_gp(0:nprocb-1))
      
      DO k = 0, nprocb-1
        
        ! Fourier Space -> GridPoint
        
        nk    = gl_dc(imype)%nflevp1
        nk0   = gl_dc(imype)%nflev
        nglat = gl_dc(idest(k))%nglat
        nglon = gl_dc(idest(k))%nglon

        ALLOCATE (fs_gp(k)%send_buffer(nglon, nk, nglat, nvarmax) )
        ALLOCATE (fs_gp(k)%send_buffer0(nk0, nglat, nva0))

!$OMP PARALLEL PRIVATE(it,ils,ile)
        it = omp_get_thread_num()
        ils = ilats(it)
        ile = ilate(it)
        fs_gp(k)%send_buffer(:,:,ils:ile,:) = 0.0_dp
        fs_gp(k)%send_buffer0(:,ils:ile,:)  = 0.0_dp
!$OMP END PARALLEL
       
        nk    = gl_dc(isrc(k))%nflevp1
        nk0   = gl_dc(isrc(k))%nflev
        nglat = gl_dc(imype)%nglat
        nglon = gl_dc(imype)%nglon

        ALLOCATE (fs_gp(k)%recv_buffer(nglon, nk, nglat, nvarmax) )
        ALLOCATE (fs_gp(k)%recv_buffer0(nk0, nglat, nva0))

      ENDDO
    
    ENDIF
    !
    !------------------------------------------------------------------------
    !
    SELECT CASE (sign)

    CASE (foreward)
      ! 
      ! GridPoint -> Fourier Space
      !
      DO k = 0, nprocb-1
         
        lm0s(k)  = (imype == idest(k) .OR. sign == backward ) .AND. PRESENT(fs0)
        IF (imype == idest(k)) THEN
          lm0r(k) = lm0s(k)
        ELSE
          lm0r(k) = (imype == isrc(k) .OR. sign == backward) .AND. PRESENT(fs0)
        ENDIF
        
        ks    = gl_dc(idest(k))%flevs
        ke    = gl_dc(idest(k))%fleve
        nk    = gl_dc(idest(k))%nflevp1
        nk0   = gl_dc(idest(k))%nflev
        nglat = gl_dc(imype)%nglat
        nglon = gl_dc(imype)%nglon
        glons = gl_dc(imype)%glons
        glone = gl_dc(imype)%glone
        nlon  = gl_dc(imype)%nlon
        nglh  = gl_dc(imype)%nglh
        
        CALL pack_gp_buf(k)
        
      ENDDO
      
      CALL sendrecv_gpfs(gp_fs)

      DO k = 0, nprocb-1
         
        ks    = gl_dc(imype)%flevs
        ke    = gl_dc(imype)%fleve
        nk    = gl_dc(imype)%nflevp1
        nk0   = gl_dc(imype)%nflev
        nglat = gl_dc(isrc(k))%nglat
        nglon = gl_dc(isrc(k))%nglon
        glons = gl_dc(isrc(k))%glons
        glone = gl_dc(isrc(k))%glone
        nlon  = gl_dc(isrc(k))%nlon
        nglh  = gl_dc(isrc(k))%nglh
        
        CALL unpack_buf_fs(k)
        
      ENDDO

!$OMP PARALLEL PRIVATE(it,ils,ile,iv,il,k)
      it = omp_get_thread_num()
      ils = ilats(it)
      ile = ilate(it)
      IF (SIZE(fs,1) == gl_dc(imype)%nlon+2) THEN
        DO iv = 1, SIZE(fs,4)
          DO il = ils, ile
            DO k = 1, SIZE(fs,2)
              fs (gl_dc(imype)%nlon+1,k,il,iv) = 0.0_dp
              fs (gl_dc(imype)%nlon+2,k,il,iv) = 0.0_dp
            ENDDO
          ENDDO
        ENDDO
      ELSE
        fs (gl_dc(imype)%nlon+1:,:,ils:ile,:) = 0.0_dp
      ENDIF
!$OMP END PARALLEL

    CASE (backward)
      ! 
      ! Fourier Space -> GridPoint
      !
      DO k = 0, nprocb-1
         
        lm0s(k)  = (imype == idest(k) .OR. sign == backward) .AND. PRESENT(fs0)
        IF (imype == idest(k)) THEN
          lm0r(k) = lm0s(k)
        ELSE
          lm0r(k) = (imype == isrc(k) .OR. sign == backward) .AND. PRESENT(fs0)
        ENDIF
        
        ks    = gl_dc(imype)%flevs
        ke    = gl_dc(imype)%fleve
        nk    = gl_dc(imype)%nflevp1
        nk0   = gl_dc(imype)%nflev
        nglat = gl_dc(idest(k))%nglat
        nglon = gl_dc(idest(k))%nglon
        glons = gl_dc(idest(k))%glons
        glone = gl_dc(idest(k))%glone
        nlon  = gl_dc(idest(k))%nlon
        nglh  = gl_dc(idest(k))%nglh
        
        CALL pack_fs_buf(k)
        
      ENDDO

      CALL sendrecv_gpfs(fs_gp)

      DO k = 0, nprocb-1
        
        ks    = gl_dc(isrc(k))%flevs
        ke    = gl_dc(isrc(k))%fleve
        nk    = gl_dc(isrc(k))%nflevp1
        nk0   = gl_dc(isrc(k))%nflev
        nglat = gl_dc(imype)%nglat
        nglon = gl_dc(imype)%nglon
        glons = gl_dc(imype)%glons
        glone = gl_dc(imype)%glone
        nlon  = gl_dc(imype)%nlon
        nglh  = gl_dc(imype)%nglh

        CALL unpack_buf_gp(k)

      ENDDO

    CASE default

      CALL finish ('transpose_fs_gp', 'invalid transposition parameter (not foreward,backward)')

    END SELECT
    !------------------------------------------------------------------------------
    !
  CONTAINS
    !
    !------------------------------------------------------------------------------
    SUBROUTINE pack_gp_buf(kl)
      INTEGER, INTENT(in) :: kl
      !
      ! pack message to send/recv buffer buf
      !
!$OMP PARALLEL PRIVATE(it,ils,ile,ibs,ibe), FIRSTPRIVATE(ke,nk)
      it = omp_get_thread_num()
      ils = ilats(it)
      ile = ilate(it)
      ibs = iblks(it)
      ibe = iblke(it)
      !
      IF (ke == gl_dc(imype)%nlev+1) THEN
        !
        gp_fs(kl)%send_buffer (:,nk,ils:ile,:) = 0.0_dp
        !
        IF (lreg) THEN
          IF(PRESENT(sf1)) gp_fs(kl)%send_buffer (:,nk,ils:ile,1) = sf1(:nglon,ils:ile)
          IF(PRESENT(sf2)) gp_fs(kl)%send_buffer (:,nk,ils:ile,2) = sf2(:nglon,ils:ile)
          IF(PRESENT(sf3)) gp_fs(kl)%send_buffer (:,nk,ils:ile,3) = sf3(:nglon,ils:ile)
        ELSE
          IF(PRESENT(sf1)) CALL reorder (gp_fs(kl)%send_buffer(:,nk,ils:ile,1), sf1(:,:), ils, ile)
          IF(PRESENT(sf2)) CALL reorder (gp_fs(kl)%send_buffer(:,nk,ils:ile,2), sf2(:,:), ils, ile)
          IF(PRESENT(sf3)) CALL reorder (gp_fs(kl)%send_buffer(:,nk,ils:ile,3), sf3(:,:), ils, ile)
        ENDIF
        ke = ke - 1
        nk = nk - 1
      ENDIF
      !
      IF (nk > 0) THEN
        IF (lreg) THEN
          gp_fs(kl)%send_buffer (:,:nk,ils:ile,1) = gp1 (:nglon,ks:ke,ils:ile)
          gp_fs(kl)%send_buffer (:,:nk,ils:ile,2) = gp2 (:nglon,ks:ke,ils:ile)
          gp_fs(kl)%send_buffer (:,:nk,ils:ile,3) = gp3 (:nglon,ks:ke,ils:ile)
          gp_fs(kl)%send_buffer (:,:nk,ils:ile,4) = gp4 (:nglon,ks:ke,ils:ile)
          gp_fs(kl)%send_buffer (:,:nk,ils:ile,5) = gp5 (:nglon,ks:ke,ils:ile)
          gp_fs(kl)%send_buffer (:,:nk,ils:ile,6) = gp6 (:nglon,ks:ke,ils:ile)
          gp_fs(kl)%send_buffer (:,:nk,ils:ile,7) = gp7 (:nglon,ks:ke,ils:ile)
        ELSE
          CALL reorder(gp_fs(kl)%send_buffer (:,:nk,ils:ile,1), gp1 (:,ks:ke,:), ils, ile)
          CALL reorder(gp_fs(kl)%send_buffer (:,:nk,ils:ile,2), gp2 (:,ks:ke,:), ils, ile)
          CALL reorder(gp_fs(kl)%send_buffer (:,:nk,ils:ile,3), gp3 (:,ks:ke,:), ils, ile)
          CALL reorder(gp_fs(kl)%send_buffer (:,:nk,ils:ile,4), gp4 (:,ks:ke,:), ils, ile)
          CALL reorder(gp_fs(kl)%send_buffer (:,:nk,ils:ile,5), gp5 (:,ks:ke,:), ils, ile)
          CALL reorder(gp_fs(kl)%send_buffer (:,:nk,ils:ile,6), gp6 (:,ks:ke,:), ils, ile)
          CALL reorder(gp_fs(kl)%send_buffer (:,:nk,ils:ile,7), gp7 (:,ks:ke,:), ils, ile)
        ENDIF
      ENDIF
      !
      ! pack zonal mean
      !
      IF (lm0s(kl)) THEN
        IF(PRESENT(zm1)) gp_fs(kl)%send_buffer0 (:,ils:ile,1) = zm1 (ks:ke,ils:ile)
        IF(PRESENT(zm2)) gp_fs(kl)%send_buffer0 (:,ils:ile,2) = zm2 (ks:ke,ils:ile)
        IF(PRESENT(zm3)) gp_fs(kl)%send_buffer0 (:,ils:ile,3) = zm3 (ks:ke,ils:ile)
      ENDIF
!$OMP END PARALLEL
      !
    END SUBROUTINE pack_gp_buf
    !------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_gp(kl)
      INTEGER, INTENT(in) :: kl
      !
      ! unpack grid point space from send/recv buffer buf
      !   
!$OMP PARALLEL PRIVATE(it,ils,ile,ibs,ibe), FIRSTPRIVATE(ke,nk)
      it = omp_get_thread_num()
      ils = ilats(it)
      ile = ilate(it)
      ibs = iblks(it)
      ibe = iblke(it)
      !
      IF (ke == gl_dc(imype)%nlev+1) THEN
        !
        IF (lreg) THEN
          IF(PRESENT(sf1)) sf1(:,ils:ile) = fs_gp(kl)%recv_buffer (:,nk,ils:ile,1)
          IF(PRESENT(sf2)) sf2(:,ils:ile) = fs_gp(kl)%recv_buffer (:,nk,ils:ile,2)
          IF(PRESENT(sf3)) sf3(:,ils:ile) = fs_gp(kl)%recv_buffer (:,nk,ils:ile,3)
        ELSE
          IF(PRESENT(sf1)) CALL reorder (sf1(:,ibs:ibe), fs_gp(kl)%recv_buffer (:,nk,:,1), ibs, ibe)
          IF(PRESENT(sf2)) CALL reorder (sf2(:,ibs:ibe), fs_gp(kl)%recv_buffer (:,nk,:,2), ibs, ibe)
          IF(PRESENT(sf3)) CALL reorder (sf3(:,ibs:ibe), fs_gp(kl)%recv_buffer (:,nk,:,3), ibs, ibe)
        ENDIF
        ke = ke - 1
        nk = nk - 1
        !
      ENDIF
      !
      IF (nk > 0) THEN
        IF (lreg) THEN
          gp1 (:,ks:ke,ils:ile) = fs_gp(kl)%recv_buffer (:,:nk,ils:ile,1) 
          gp2 (:,ks:ke,ils:ile) = fs_gp(kl)%recv_buffer (:,:nk,ils:ile,2)
          gp3 (:,ks:ke,ils:ile) = fs_gp(kl)%recv_buffer (:,:nk,ils:ile,3)
          gp4 (:,ks:ke,ils:ile) = fs_gp(kl)%recv_buffer (:,:nk,ils:ile,4)
          gp5 (:,ks:ke,ils:ile) = fs_gp(kl)%recv_buffer (:,:nk,ils:ile,5)
          gp6 (:,ks:ke,ils:ile) = fs_gp(kl)%recv_buffer (:,:nk,ils:ile,6)
          gp7 (:,ks:ke,ils:ile) = fs_gp(kl)%recv_buffer (:,:nk,ils:ile,7)
          IF (PRESENT(gp8)) gp8 (:,ks:ke,ils:ile) = fs_gp(kl)%recv_buffer (:,:nk,ils:ile,8)
          IF (PRESENT(gp9)) gp9 (:,ks:ke,ils:ile) = fs_gp(kl)%recv_buffer (:,:nk,ils:ile,9)
        ELSE
          CALL reorder (gp1 (:,ks:ke,ibs:ibe), fs_gp(kl)%recv_buffer (:,:nk,:,1), ibs, ibe)
          CALL reorder (gp2 (:,ks:ke,ibs:ibe), fs_gp(kl)%recv_buffer (:,:nk,:,2), ibs, ibe)
          CALL reorder (gp3 (:,ks:ke,ibs:ibe), fs_gp(kl)%recv_buffer (:,:nk,:,3), ibs, ibe)
          CALL reorder (gp4 (:,ks:ke,ibs:ibe), fs_gp(kl)%recv_buffer (:,:nk,:,4), ibs, ibe)
          CALL reorder (gp5 (:,ks:ke,ibs:ibe), fs_gp(kl)%recv_buffer (:,:nk,:,5), ibs, ibe)
          CALL reorder (gp6 (:,ks:ke,ibs:ibe), fs_gp(kl)%recv_buffer (:,:nk,:,6), ibs, ibe)
          CALL reorder (gp7 (:,ks:ke,ibs:ibe), fs_gp(kl)%recv_buffer (:,:nk,:,7), ibs, ibe)
          IF (PRESENT(gp8)) CALL reorder (gp8 (:,ks:ke,ibs:ibe), fs_gp(kl)%recv_buffer (:,:nk,:,8), ibs, ibe)
          IF (PRESENT(gp9)) CALL reorder (gp9 (:,ks:ke,ibs:ibe), fs_gp(kl)%recv_buffer (:,:nk,:,9), ibs, ibe)
        ENDIF
      ENDIF
      ! 
      ! unpack zonal mean
      !   
      IF (lm0r(kl)) THEN
        IF(PRESENT(zm1)) zm1 (ks:ke,ils:ile) = fs_gp(kl)%recv_buffer0 (:,ils:ile,1)
        IF(PRESENT(zm2)) zm2 (ks:ke,ils:ile) = fs_gp(kl)%recv_buffer0 (:,ils:ile,2)
        IF(PRESENT(zm3)) zm3 (ks:ke,ils:ile) = fs_gp(kl)%recv_buffer0 (:,ils:ile,3)
      ENDIF
!$OMP END PARALLEL
      !
    END SUBROUTINE unpack_buf_gp
    !------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_fs(kl)
      INTEGER, INTENT(in) :: kl
      !
      ! unpack message to fourier buffer fs
      !
      INTEGER :: ils1, ile1, ils2, ile2, it, ils, ile
      !
!$OMP PARALLEL PRIVATE(it, ils, ile, ils1, ile1, ils2, ile2)
      it = omp_get_thread_num()
      ils = ilats(it)
      ile = ilate(it)
      ils1 = ils ! start-latitude in first segment for this thread
      ile1 = ile ! end- " "
      ils2 = ils ! start-latitude in second segment for this thread
      ile2 = ile ! end- " "
      IF (ile <= nglat/2) THEN
        ile2 = -1
      ELSE IF (ils > nglat/2) THEN
        ile1 = -1
      ELSE
        ile1 = nglat/2
        ils2 = ile1 + 1
      ENDIF
      !
      ! unpack first segment
      !
      IF (ile1 > 0) THEN
        fs(glons(1):glone(1),:,ils1:ile1,:nvar) = gp_fs(kl)%recv_buffer(:,:,ils1:ile1,:nvar)
      ENDIF
      ! 
      ! unpack second segment
      !
      IF (ile2 > 0) THEN
        IF (glone(2) > glons(2)) THEN
          fs(glons(2):glone(2), :, ils2:ile2, :nvar) = gp_fs(kl)%recv_buffer(:, :, ils2:ile2, :nvar)
        ELSE
          ! 
          ! unpack second segment, split into longitudes
          !
          fs(glons(2):nlon, :, ils2:ile2,:nvar) = &
               gp_fs(kl)%recv_buffer(:nlon-glons(2)+1, :, ils2:ile2, :nvar)
          fs(1:glone(2), :, ils2:ile2, :nvar) = &
               gp_fs(kl)%recv_buffer(nglon-glone(2)+1:, :, ils2:ile2, :nvar)          
        ENDIF
      ENDIF
      ! 
      ! unpack zonal mean
      !
      IF (lm0r(kl)) THEN
        fs0 (:,ils:ile,:) = gp_fs(kl)%recv_buffer0 (:,ils:ile,:)
      ENDIF
!$OMP END PARALLEL
      !
    END SUBROUTINE unpack_buf_fs
    !------------------------------------------------------------------------------
    SUBROUTINE pack_fs_buf(kl)
      INTEGER, INTENT(in) :: kl
      ! 
      ! pack fourier buffer fs to buffer
      ! 
      INTEGER :: ils1, ile1, ils2, ile2, it, ils, ile
      !
      IF (nglh(1) < 1) THEN
        CALL finish ('pack_fs_buf','unexpected case: nglh(1) < 1',1)
      ENDIF
      !
!$OMP PARALLEL private(it, ils,ile, ils1,ile1, ils2,ile2)
      it = omp_get_thread_num()
      ils = ilats(it)
      ile = ilate(it)
      ils1 = ils ! start-latitude in first segment for this thread
      ile1 = ile ! end- " "
      ils2 = ils ! start-latitude in second segment for this thread
      ile2 = ile ! end- " "
      !
      IF (ile <= nglh(1)) THEN
        ile2 = -1
      ELSE IF (ils > nglh(1)) THEN
        ile1 = -1
      ELSE
        ile1 = nglh(1)
        ils2 = ile1 + 1
      ENDIF
      ! 
      ! pack first segment
      !
      IF (ile1 > 0) THEN
        fs_gp(kl)%send_buffer(:,:,ils1:ile1,:) = fs(glons(1):glone(1),:,ils1:ile1,:)
      ENDIF
      ! 
      ! pack second segment
      !
      IF (ile2 > 0) THEN
        IF (glone(2) > glons(2)) THEN
          fs_gp(kl)%send_buffer(:,:,ils2:ile2,:) = fs(glons(2):glone(2),:,ils2:ile2,:)
        ELSE
          ! 
          ! pack second segment, split into longitudes
          !
          fs_gp(kl)%send_buffer  (                :nlon-glons(2)+1,:,ils2:ile2,:) = &
               fs (glons(2)        :nlon           ,:,ils2:ile2,:) 
          fs_gp(kl)%send_buffer  (nglon-glone(2)+1:               ,:,ils2:ile2,:) = &
               fs (1               :glone(2)       ,:,ils2:ile2,:)
        ENDIF
      ENDIF
      ! 
      ! pack zonal mean
      !
      IF (lm0s(kl)) THEN
        fs_gp(kl)%send_buffer0  (:,ils:ile,:) = fs0 (:,ils:ile,:)
      ENDIF
!$OMP END PARALLEL
      !
    END SUBROUTINE pack_fs_buf
    !------------------------------------------------------------------------------
    SUBROUTINE sendrecv_gpfs(trbuf)
      !
      ! send and receive buffer
      ! deallocate send buffer
      !
      TYPE(transpose_buffer), INTENT(inout) :: trbuf(0:)
      !
      INTEGER :: it, ils, ile, kl
      !
#ifdef __SENDRECV
      INTEGER :: communication_type = sendrecv
#else
      INTEGER :: communication_type = non_blocking
#endif
      !
      SELECT CASE(communication_type)
        
      CASE(sendrecv)
        !
        DO kl = 0, nprocb-1
          !
          IF (imype /= idest(kl)) THEN
            CALL p_sendrecv (trbuf(kl)%send_buffer, gl_dc(idest(kl))%pe, &
                             trbuf(kl)%recv_buffer, gl_dc(isrc(kl))%pe,  &
                             tag_tr_gp_fs)
            IF (lm0s(kl) .AND. lm0r(kl)) THEN
              CALL p_sendrecv (trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                               trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe,  &
                               tag_tr_gp_fs)
            ELSE IF (lm0s(kl)) THEN
              CALL p_send (trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                           tag_tr_gp_fs)
            ELSE IF (lm0r(k)) THEN
              CALL p_recv (trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe, &
                           tag_tr_gp_fs)
            ENDIF
          ELSE
!$OMP PARALLEL PRIVATE(it,ils,ile)
            it = omp_get_thread_num()
            ils = ilats(it)
            ile = ilate(it)
            trbuf(kl)%recv_buffer(:,:,ils:ile,:) = trbuf(kl)%send_buffer(:,:,ils:ile,:)
            IF (lm0r(kl)) THEN
              trbuf(kl)%recv_buffer0(:,ils:ile,:) = trbuf(kl)%send_buffer0(:,ils:ile,:)
            ENDIF
!$OMP END PARALLEL
          ENDIF
        ENDDO
        !
      CASE (non_blocking)
        !
        DO kl = 0, nprocb-1
          IF (imype /= idest(kl)) THEN
            CALL p_isend (trbuf(kl)%send_buffer, gl_dc(idest(kl))%pe, &
                          tag_tr_gp_fs)
            CALL p_irecv (trbuf(kl)%recv_buffer, gl_dc(isrc(kl))%pe,  &
                          tag_tr_gp_fs)
            IF (lm0s(kl) .AND. lm0r(kl)) THEN
              CALL p_isend(trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                           tag_tr_gp_fs)
              CALL p_irecv (trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe,  &
                            tag_tr_gp_fs)
            ELSE IF (lm0s(kl)) THEN
              CALL p_isend (trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                            tag_tr_gp_fs)
            ELSE IF (lm0r(kl)) THEN
              CALL p_irecv (trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe, &
                            tag_tr_gp_fs)
            ENDIF
          ELSE
!$OMP PARALLEL PRIVATE(it,ils,ile)
            it = omp_get_thread_num()
            ils = ilats(it)
            ile = ilate(it)
            trbuf(kl)%recv_buffer(:,:,ils:ile,:) = trbuf(kl)%send_buffer(:,:,ils:ile,:)
            IF (lm0r(kl)) THEN
              trbuf(kl)%recv_buffer0(:,ils:ile,:) = trbuf(kl)%send_buffer0(:,ils:ile,:)
            ENDIF
!$OMP END PARALLEL
          ENDIF
        ENDDO
        !
        CALL p_wait
        !
      END SELECT
      !
    END SUBROUTINE sendrecv_gpfs
    !
  END SUBROUTINE transpose_gp_fs
  !==============================================================================
  SUBROUTINE transpose_fs_ls (gl_dc, sign, fs, ls, fs0, ls0)
    !
    ! transpose
    !   sign=foreward : Fourier space  -> Legendre space
    !   sign=backward : Fourier space <-  Legendre space
    !
    TYPE (pe_decomposed) ,INTENT(in)     :: gl_dc  (:)       ! decomposition
    INTEGER              ,INTENT(in)     :: sign             ! foreward:fs>ls; backward:gs<ls
    !
    ! Assumed shape array association:
    !
    REAL(dp),           INTENT(inout)  :: fs   (:,:,:,:)   ! fs
    REAL(dp),           INTENT(inout)  :: ls   (:,:,:,:)   ! ls
    REAL(dp), OPTIONAL, INTENT(inout)  :: fs0  (:,:,:)     ! fs, zonal means
    REAL(dp), OPTIONAL, INTENT(inout)  :: ls0  (:,:,:)     ! ls, zonal means
    !
    ! Array element sequence association to speed up indexing:
    !
    ! local variables
    !
    INTEGER :: i, k, n        ! loop indices
    INTEGER :: nlev           ! number of levels
    INTEGER :: nlev0          ! number of levels (m=0 only)
    INTEGER :: nflat          ! number of latitudes (2*nhgl)
    INTEGER :: flats(2)       ! first latitude  in Fourier space
    INTEGER :: flate(2)       ! last  latitude  in Fourier space
    INTEGER :: nlm            ! no. of m coeff. in Legendre space
    INTEGER :: nproca         ! number of PEs in set A
    INTEGER :: n2mp1          ! total number of coeff. from lgti 
    INTEGER :: nb,nf,n2       ! explicit shapes
    !
    INTEGER, POINTER :: intr(:)        ! index array
    !
    LOGICAL, ALLOCATABLE, SAVE :: lm0r(:)        ! receive m0
    LOGICAL, ALLOCATABLE, SAVE :: lm0s(:)        ! send m0
    !
    INTEGER                    :: imype    ! index of this pe
    INTEGER, ALLOCATABLE, SAVE :: idest(:) ! destination of data
    INTEGER, ALLOCATABLE, SAVE :: isrc(:)  ! source of data
    LOGICAL, SAVE :: first_call=.TRUE.
    LOGICAL, SAVE :: simple_intr
    INTEGER, SAVE :: intrs,intre
    INTEGER :: ialloc
    !
    if (require_init) then
      call init_transpose
    endif
    !
    imype  = indx (p_pe, gl_dc)                 ! pe index (my)
    nproca =  gl_dc(imype)%nproca
    n2mp1  = (gl_dc(imype)%nm + 1) * 2
    !
    if (first_call) then
      intrs = gl_dc(imype)%intr(1)
      intre = gl_dc(imype)%intr(2*gl_dc(imype)%nlm)
      simple_intr=.true.
      do i = 1, 2*gl_dc(imype)%nlm
        if (gl_dc(imype)%intr(i) /= intrs-1+i) then
          simple_intr = .FALSE.
          EXIT
        ENDIF
      ENDDO
      first_call = .FALSE.
    ENDIF

!!$#ifdef _OPEMP
!!$    IF (nproca == 1) THEN
!!$      intr  => gl_dc(imype)%intr
!!$      CALL pure_omp()
!!$      RETURN
!!$    ENDIF
!!$#endif
    !
    IF (.NOT. ALLOCATED(plan_a)) THEN

      ALLOCATE(plan_a(0:nproca-1))

      k = 0
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)%set_b /= gl_dc(imype)%set_b) CYCLE
        plan_a(k) = i          ! gl_dc(i)%pe
        IF (i == imype) n = k
        k = k + 1
      ENDDO

      plan_a = CSHIFT (plan_a,n)

      ALLOCATE(idest(0:nproca-1))
      ALLOCATE(isrc(0:nproca-1))

      DO k = 0, nproca-1
        idest(k) = plan_a(           k        ) ! PE index (send)
        isrc(k)  = plan_a(MOD(nproca-k,nproca)) ! PE index (recv)
      ENDDO

      ALLOCATE(lm0r(0:nproca-1))
      ALLOCATE(lm0s(0:nproca-1))
    ENDIF
    !
    !------------------------------------------------------------------------
    !
    IF (.NOT. ALLOCATED(fs_ls)) THEN

      ALLOCATE (fs_ls(0:nproca-1))
      ALLOCATE (ls_fs(0:nproca-1))

      DO ialloc = 1, 8
 
        DO k = 0, nproca-1

          ! Fourier -> Legendre space

          nlm    = gl_dc(idest(k))%nlm
          nlev   = gl_dc(imype)%nflevp1
          nlev0  = gl_dc(imype)%nflev
          nflat  = gl_dc(imype)%nflat
          
          IF (ialloc == 1) ALLOCATE (fs_ls(k)%send_buffer(2*nlm,nlev,nflat,nvar_fsls))
          IF (ialloc == 3) ALLOCATE (fs_ls(k)%send_buffer0(nlev0,nflat,nvar0_fsls))
          
          nlm    = gl_dc(imype)%nlm
          nlev   = gl_dc(isrc(k))%nflevp1
          nlev0  = gl_dc(isrc(k))%nflev
          nflat  = gl_dc(isrc(k))%nflat
          
          IF (ialloc == 2) ALLOCATE (fs_ls(k)%recv_buffer(2*nlm,nlev,nflat,nvar_fsls))
          IF (ialloc == 4) ALLOCATE (fs_ls(k)%recv_buffer0(nlev0,nflat,nvar0_fsls))
          
          ! Legendre -> Fourier space
          
          nlm    = gl_dc(imype)%nlm
          nlev   = gl_dc(idest(k))%nflevp1
          nlev0  = gl_dc(idest(k))%nflev
          nflat  = gl_dc(idest(k))%nflat
          
          IF (ialloc == 5) ALLOCATE (ls_fs(k)%send_buffer(2*nlm,nlev,nflat,nvar_lsfs))
          IF (ialloc == 7) ALLOCATE (ls_fs(k)%send_buffer0(nlev0,nflat,nvar0_lsfs))
          
          nlm    = gl_dc(isrc(k))%nlm
          nlev   = gl_dc(imype)%nflevp1
          nlev0  = gl_dc(imype)%nflev
          nflat  = gl_dc(imype)%nflat
          
          IF (ialloc == 6) ALLOCATE (ls_fs(k)%recv_buffer(2*nlm,nlev,nflat,nvar_lsfs))
          IF (ialloc == 8) ALLOCATE (ls_fs(k)%recv_buffer0(nlev0,nflat,nvar0_lsfs))
          
        ENDDO

      ENDDO

    ENDIF
    !------------------------------------------------------------------------
    SELECT CASE (sign)

    CASE (foreward)
      !
      ! Fourier -> Legendre space
      !
      nlev  =  gl_dc(imype)%nflevp1
      nlev0 =  gl_dc(imype)%nflev
      nflat =  gl_dc(imype)%nflat
      flats =  gl_dc(imype)%flats
      flate =  gl_dc(imype)%flate
      nf = SIZE (fs,1)
      n2 = nlev * nflat * nvar_fsls

      DO k = 0, nproca-1
        lm0s(k)  = gl_dc(idest(k))%nlnm0 > 0 .AND. PRESENT(fs0)
        IF (imype == idest(k)) THEN 
          lm0r(k) = lm0s(k)
        ELSE
          lm0r(k) = gl_dc(imype)%nlnm0 > 0 .AND. PRESENT(fs0)
        ENDIF
      ENDDO

!$OMP PARALLEL PRIVATE(k)
      DO k = 0, nproca-1
        CALL pack_fs_buf(k)
      ENDDO
!$OMP END PARALLEL

#ifdef __FS_LS_ALLTOALLV
      CALL alltoallv_fsls(fs_ls)
#else
      CALL sendrecv_fsls(fs_ls)
#endif

      nlm    = gl_dc(imype)%nlm
      intr  => gl_dc(imype)%intr
      nb = 2*nlm
      nf = SIZE (fs,1)
      n2 = nlev * nflat * nvar_fsls

!$OMP PARALLEL PRIVATE(k)
      DO k = 0, nproca-1
        CALL unpack_buf_ls(k)
      END DO
!$OMP END PARALLEL

    CASE (backward)
      !
      ! Legendre -> Fourier
      !
      nlm    = gl_dc(imype)%nlm
      nb = 2*nlm
      nf = SIZE (fs,1)
      intr  => gl_dc(imype)%intr

      DO k = 0, nproca-1
        lm0s(k)  = gl_dc(imype)%nlnm0 > 0 .AND. PRESENT(fs0)
        IF (imype == idest(k)) THEN 
          lm0r(k) = lm0s(k)
        ELSE
          lm0r(k) = gl_dc(isrc(k))%nlnm0 > 0 .AND. PRESENT(fs0)
        ENDIF
      ENDDO

!$OMP PARALLEL PRIVATE(k)
      DO k = 0, nproca-1
        CALL pack_ls_buf(k)
      ENDDO
!$OMP END PARALLEL

#ifdef __FS_LS_ALLTOALLV
      CALL alltoallv_fsls(ls_fs)
#else
      CALL sendrecv_fsls(ls_fs)
#endif

      DO k = 0, nproca-1

        nlm    = gl_dc(isrc(k))%nlm
        nlev   = gl_dc(imype)%nflevp1
        nlev0  = gl_dc(imype)%nflev
        nflat =  gl_dc(imype)%nflat
        flats =  gl_dc(imype)%flats
        flate =  gl_dc(imype)%flate

        intr  => gl_dc(isrc(k))%intr

        nb = 2*nlm
        nf = SIZE (fs,1)
        n2 = nlev * nflat * nvar_lsfs

        CALL unpack_buf_fs(k)

      END DO

      ! set coefficients not provided by inverse Legendre transform zero
      fs(n2mp1+1:,:,:,:) = 0.0_dp

    CASE default

      CALL finish ('transpose_fs_ls', 'invalid transposition parameter (not foreward,backward)')

    END SELECT
    !------------------------------------------------------------------------------
  CONTAINS
    !------------------------------------------------------------------------------
!!$    SUBROUTINE pure_omp
!!$
!!$      INTEGER :: it, ivar, ilat, ilats, ilate, is
!!$      LOGICAL :: lm0r_k0
!!$
!!$      nflat =  gl_dc(imype)%nflat
!!$      flats =  gl_dc(imype)%flats
!!$      flate =  gl_dc(imype)%flate
!!$      lm0r_k0  = gl_dc(imype)%nlnm0 > 0 .AND. PRESENT(fs0)
!!$
!!$!$OMP PARALLEL PRIVATE(it,ivar,ilat,ilats,ilate)
!!$      it = omp_get_thread_num()
!!$      DO is = 1, 2
!!$        ilats = thsp_lat(is,it)%begin
!!$        ilate = thsp_lat(is,it)%end
!!$
!!$        SELECT CASE (sign)
!!$
!!$        CASE (foreward)
!!$          ! f2l:
!!$          DO ivar = 1, nvar_fsls
!!$            IF (simple_intr) THEN
!!$              DO ilat = ilats, ilate
!!$                ls(:,:,ilat,ivar) = fs(intrs:intre,:,ilat,ivar)
!!$              ENDDO
!!$            ELSE
!!$              DO ilat = ilats, ilate
!!$                ls(:,:,ilat,ivar) = fs(intr,:,ilat,ivar)
!!$              ENDDO
!!$            ENDIF
!!$          ENDDO
!!$          IF (lm0r_k0) THEN
!!$            DO ivar = 1, nvar0_lsfs
!!$              DO ilat = ilats, ilate
!!$                ls0 (:,ilat,ivar) = fs0(:,ilat,ivar)
!!$              ENDDO
!!$            ENDDO
!!$          ENDIF
!!$        CASE (backward)
!!$          ! l2f:
!!$          DO ivar = 1, nvar_fsls
!!$            IF (simple_intr) THEN
!!$              DO ilat = ilats, ilate
!!$                fs (intrs:intre,:,ilat,ivar)     = ls(:,:,ilat,ivar)
!!$                fs (n2mp1+1:,:,ilat,ivar) = 0.0_dp
!!$              ENDDO
!!$            ELSE
!!$              DO ilat = ilats, ilate
!!$                fs (intr,:,ilat,ivar)     = ls(:,:,ilat,ivar)
!!$                fs (n2mp1+1:,:,ilat,ivar) = 0.0_dp
!!$              ENDDO
!!$            ENDIF
!!$          ENDDO
!!$          IF (lm0r_k0) THEN
!!$            DO ivar = 1, nvar0_lsfs
!!$              DO ilat = ilats, ilate
!!$                fs0(:,ilat,ivar)   = ls0 (:,ilat,ivar)
!!$              ENDDO
!!$            ENDDO
!!$          ENDIF
!!$        CASE default
!!$          CALL finish ('transpose_fs_ls', 'invalid transposition parameter (not foreward,backward)')
!!$        END SELECT
!!$      ENDDO
!!$!$OMP END PARALLEL
!!$
!!$    END SUBROUTINE pure_omp

    !------------------------------------------------------------------------------

    SUBROUTINE sendrecv_fsls(trbuf)
      !
      ! send and receive buffer
      ! deallocate send buffer
      !
      TYPE(transpose_buffer), INTENT(inout) :: trbuf(0:)
      INTEGER :: kl
      !
#ifdef __SENDRECV
      INTEGER :: communication_type = sendrecv
#else
      INTEGER :: communication_type = non_blocking
#endif
      !
      SELECT CASE (communication_type)

      CASE(sendrecv)

        DO kl = 0, nproca-1
          !
          IF(imype /= idest(kl)) THEN
            CALL p_sendrecv (trbuf(kl)%send_buffer, gl_dc(idest(kl))%pe, &
                             trbuf(kl)%recv_buffer, gl_dc(isrc(kl))%pe,  &
                             tag_tr_fs_ls)
            IF (lm0s(kl) .AND. lm0r(kl)) THEN
              CALL p_sendrecv (trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                               trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe,  &
                               tag_tr_fs_ls)
            ELSE IF (lm0s(kl)) THEN
              CALL p_send (trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                           tag_tr_fs_ls)
            ELSE IF (lm0r(kl)) THEN
              CALL p_recv (trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe, &
                           tag_tr_fs_ls)
            ENDIF
          ELSE
            trbuf(kl)%recv_buffer = trbuf(kl)%send_buffer
            IF (lm0r(kl)) trbuf(kl)%recv_buffer0 = trbuf(kl)%send_buffer0
          ENDIF

        ENDDO

      CASE(non_blocking)

        DO kl = 0, nproca-1
          !
          IF(imype /= idest(kl)) THEN
            CALL p_isend(trbuf(kl)%send_buffer, gl_dc(idest(kl))%pe, &
                         tag_tr_fs_ls)
            CALL p_irecv(trbuf(kl)%recv_buffer, gl_dc(isrc(kl))%pe,  &
                         tag_tr_fs_ls)
            IF (lm0s(kl) .AND. lm0r(kl)) THEN
              CALL p_isend(trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                           tag_tr_fs_ls)
              CALL p_irecv(trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe,  &
                           tag_tr_fs_ls)
            ELSE IF (lm0s(kl)) THEN
              CALL p_isend (trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                            tag_tr_fs_ls)
            ELSE IF (lm0r(kl)) THEN
              CALL p_irecv (trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe, &
                            tag_tr_fs_ls)
            ENDIF
          ELSE
            trbuf(kl)%recv_buffer = trbuf(kl)%send_buffer
            IF (lm0r(kl)) trbuf(kl)%recv_buffer0 = trbuf(kl)%send_buffer0
          ENDIF
        ENDDO

        CALL p_wait

      END SELECT

    END SUBROUTINE sendrecv_fsls
    !------------------------------------------------------------------------------
#ifdef __FS_LS_ALLTOALLV
    SUBROUTINE alltoallv_fsls(trbuf)

      USE mpi, ONLY: MPI_ADDRESS_KIND
      USE mo_mpi, ONLY: p_communicator_a, p_communicator_b
      USE mo_util_cpu

      TYPE(transpose_buffer), INTENT(inout) :: trbuf(0:)

      INTEGER(kind=MPI_ADDRESS_KIND) :: send_base, recv_base
      INTEGER(kind=MPI_ADDRESS_KIND) :: send_base0, recv_base0

      INTEGER :: sd(0:dc%d_nprocs-1), rd(0:dc%d_nprocs-1)   ! SDISPLS, RDISPLS
      INTEGER :: sd0(0:dc%d_nprocs-1), rd0(0:dc%d_nprocs-1) !

      INTEGER :: sc(0:dc%d_nprocs-1), rc(0:dc%d_nprocs-1)   ! SENDCOUNTS, RECVCOUNTS
      INTEGER :: sc0(0:dc%d_nprocs-1), rc0(0:dc%d_nprocs-1) !

      INTEGER:: comall_mype, coma_mype, comb_mype
      INTEGER :: ia, jdest, jsrc, ierror, p_error
      INTEGER(kind=MPI_ADDRESS_KIND) :: addr
      INTEGER :: idiff
      LOGICAL, SAVE :: first_call=.TRUE.

      sc(:)  = 0
      sd(:)  = 0
      sc0(:) = 0
      sd0(:) = 0

      CALL MPI_GET_ADDRESS( trbuf(0)%send_buffer(1,1,1,1), send_base, ierror) 
      IF (ierror > 0) CALL finish('alltoallv_fsls','bad case (1a)')
      CALL MPI_GET_ADDRESS( trbuf(0)%send_buffer0(1,1,1), send_base0, ierror) 
      IF (ierror > 0) CALL finish('alltoallv_fsls','bad case (1b)')

      DO ia = 0, nproca-1
        jdest = gl_dc(idest(ia))%pe

        CALL MPI_GET_ADDRESS( trbuf(ia)%send_buffer(1,1,1,1), addr, ierror) 
        IF (ierror > 0) CALL finish('alltoallv_fsls','bad case (3)')
        idiff = INT(addr-send_base)
        sd(jdest) = idiff/8
        IF (sd(jdest)*8 /= idiff) CALL finish('alltoallv_fsls','bad alignment')
        sc(jdest) = SIZE(trbuf(ia)%send_buffer)

        CALL MPI_GET_ADDRESS( trbuf(ia)%send_buffer0(1,1,1), addr, ierror) 
        IF (ierror > 0) CALL finish('alltoallv_fsls','bad case (4)')
        idiff = INT(addr-send_base0)
        sd0(jdest) = idiff/8
        IF (sd0(jdest)*8 /= idiff) CALL finish('alltoallv_fsls','bad alignment')
        sc0(jdest) = SIZE(trbuf(ia)%send_buffer0)

        IF (.NOT. lm0s(ia)) THEN
          sc0(jdest) = 0
          sd0(jdest) = 0
        ENDIF
      ENDDO

      rc(:)  = 0
      rd(:)  = 0
      rc0(:) = 0
      rd0(:) = 0
      CALL MPI_GET_ADDRESS( trbuf(0)%recv_buffer(1,1,1,1), recv_base, ierror) 
      IF (ierror > 0) CALL finish('alltoallv_fsls','bad case (5)')
      CALL MPI_GET_ADDRESS( trbuf(0)%recv_buffer0(1,1,1), recv_base0, ierror) 
      IF (ierror > 0) CALL finish('alltoallv_fsls','bad case (6)')

      DO ia=0, nproca-1
        jsrc=gl_dc(isrc(ia))%pe

        CALL MPI_GET_ADDRESS( trbuf(ia)%recv_buffer(1,1,1,1), addr, ierror) 
        IF (ierror > 0) CALL finish('alltoallv_fsls','bad case (7)')
        idiff=INT(addr-recv_base)
        rd(jsrc)=idiff/8
        IF ( rd(jsrc)*8 /= idiff ) CALL finish('alltoallv_fsls','bad alignment')
        rc(jsrc)=SIZE(trbuf(ia)%recv_buffer)

        CALL MPI_GET_ADDRESS( trbuf(ia)%recv_buffer0(1,1,1), addr, ierror) 
        IF (ierror > 0) CALL finish('alltoallv_fsls','bad case (8)')
        idiff=INT(addr-recv_base0)
        rd0(jsrc)=idiff/8
        IF ( rd0(jsrc)*8 /= idiff ) CALL finish('alltoallv_fsls','bad alignment')
        rc0(jsrc)=SIZE(trbuf(ia)%recv_buffer0)

        IF (.NOT. lm0r(ia)) THEN
          rc0(jsrc)=0
          rd0(jsrc)=0
        ENDIF

      ENDDO

      CALL MPI_ALLTOALLV(trbuf(0)%send_buffer(1,1,1,1), sc, sd, p_real_dp, & !SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPE,
           &             trbuf(0)%recv_buffer(1,1,1,1), rc, rd, p_real_dp, & !RECVBUF, RECVCOUNTS, RDISPLS, RECVTYPE
           &             p_all_comm, ierror)
      
      IF (ierror > 0) CALL finish('alltoallv_fsls','MPI_ALLTOALLV (1) failed')

      CALL MPI_ALLTOALLV(trbuf(0)%send_buffer0, sc0, sd0, p_real_dp, & !SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPE,
           &             trbuf(0)%recv_buffer0, rc0, rd0, p_real_dp, & !RECVBUF, RECVCOUNTS, RDISPLS, RECVTYPE
           &             p_all_comm, ierror)

      IF (ierror > 0) CALL finish('alltoallv_fsls','MPI_ALLTOALLV (2) failed')

    END SUBROUTINE alltoallv_fsls
#endif
    !------------------------------------------------------------------------------
    SUBROUTINE pack_fs_buf(kl)
      INTEGER, INTENT(in) :: kl
      INTEGER :: nlm, nb
      INTEGER, POINTER :: intr(:)

      INTEGER :: it, ifls, ifle
#ifdef __EXPLICIT
      INTEGER :: ilat, ivar, ilev, iq
#endif

      it = omp_get_thread_num()
      ifls = thsp_flat(it)%begin
      ifle = thsp_flat(it)%end

      nlm = gl_dc(idest(kl))%nlm
      nb  = 2*nlm
      
      intr => gl_dc(idest(kl))%intr

#ifdef __EXPLICIT
      DO ivar = 1, nvar_fsls
        DO ilat = ifls, ifle
          DO ilev = 1, nlev
            DO iq = 1, nb              
              fs_ls(kl)%send_buffer(iq,ilev,ilat,ivar) = fs(intr(iq),ilev,ilat,ivar)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
#else
      fs_ls(kl)%send_buffer(:,:,ifls:ifle,:) = fs(intr,:,ifls:ifle,:)
#endif
      IF (lm0s(kl)) THEN
        fs_ls(kl)%send_buffer0(:,ifls:ifle,:) = fs0(:,ifls:ifle,:)
      ENDIF

    END SUBROUTINE pack_fs_buf
    !------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_fs(kl)
      INTEGER, INTENT(in) :: kl
      !     
      INTEGER :: it, ifls, ifle
#ifdef __EXPLICIT
      INTEGER :: ilat, ivar, ilev, iq
#endif

#ifdef __EXPLICIT
!$OMP PARALLEL PRIVATE(it,ifls,ifle,ivar,ilat,ilev,iq)
#else
!$OMP PARALLEL PRIVATE(it,ifls,ifle)
#endif
      it = omp_get_thread_num()
      ifls = thsp_flat(it)%begin
      ifle = thsp_flat(it)%end

#ifdef __EXPLICIT
      DO ivar = 1, nvar_fsls
        DO ilat = ifls, ifle
          DO ilev = 1, nlev
            DO iq = 1, nb              
              fs (intr(iq),ilev,ilat,ivar) = ls_fs(kl)%recv_buffer(iq,ilev,ilat,ivar) 
            ENDDO
          ENDDO
        ENDDO
      ENDDO
#else
      fs(intr,:,ifls:ifle,:) = ls_fs(kl)%recv_buffer(:,:,ifls:ifle,:)
#endif
      IF (lm0r(kl)) THEN
        fs0(:,ifls:ifle,:) = ls_fs(kl)%recv_buffer0(:,ifls:ifle,:)
      ENDIF
!$OMP END PARALLEL

    END SUBROUTINE unpack_buf_fs
    !------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_ls(kl)
      INTEGER, INTENT(in) :: kl
      INTEGER :: nlev, nlev0, nflat, flats(2), flate(2)
      !
#ifdef _OPENMP
      INTEGER :: it, is, ilats,ilate, ilat, ivar, ilev, iq, offset(2)
#endif
      
      nlev   = gl_dc(isrc(kl))%nflevp1
      nlev0  = gl_dc(isrc(kl))%nflev
      nflat =  gl_dc(isrc(kl))%nflat
      flats =  gl_dc(isrc(kl))%flats
      flate =  gl_dc(isrc(kl))%flate
      
#ifdef _OPENMP
      it = omp_get_thread_num()
      offset(1) = 1-flats(1)
      offset(2) = nflat/2+1-flats(2)
      DO is = 1, 2
        ilats=thsp_lat(is,it)%begin
        ilate=thsp_lat(is,it)%end
        DO ivar = 1, nvar_lsfs
          DO ilat = ilats, ilate
            DO ilev = 1, nlev
              DO iq= 1, 2*nlm
                ls(iq,ilev,ilat,ivar) = fs_ls(kl)%recv_buffer(iq,ilev,offset(is)+ilat,ivar)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF (lm0r(kl)) THEN
          DO ivar = 1, nvar0_lsfs
            DO ilat = ilats, ilate
              DO ilev = 1, nlev
                ls0(ilev,ilat,ivar) = fs_ls(kl)%recv_buffer0 (ilev,offset(is)+ilat,ivar)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
#else
#ifdef __EXPLICIT
      CALL explicit_unpack_buf_ls(ls(1,1,1,1),fs_ls(kl)%recv_buffer(1,1,1,1), &
           2*nlm*nlev,SIZE(LS,3),nflat,nvar_fsls,flats(1),flate(1),1,nflat/2)
      CALL explicit_unpack_buf_ls(ls(1,1,1,1),fs_ls(kl)%recv_buffer(1,1,1,1), &
           2*nlm*nlev,SIZE(LS,3),nflat,nvar_fsls,flats(2),flate(2),nflat/2+1,nflat)
      IF (lm0r(kl)) THEN
        ls0 (:,flats(1):flate(1),:) = fs_ls(kl)%recv_buffer0 (:,         :nflat/2,:)
        ls0 (:,flats(2):flate(2),:) = fs_ls(kl)%recv_buffer0 (:,nflat/2+1:       ,:)
      ENDIF
#else
      ls(:,:,flats(1):flate(1),:) = fs_ls(kl)%recv_buffer (:,:,         :nflat/2,:)
      ls(:,:,flats(2):flate(2),:) = fs_ls(kl)%recv_buffer (:,:,nflat/2+1:       ,:)
      IF (lm0r(kl)) THEN
        ls0 (:,flats(1):flate(1),:) = fs_ls(kl)%recv_buffer0 (:,         :nflat/2,:)
        ls0 (:,flats(2):flate(2),:) = fs_ls(kl)%recv_buffer0 (:,nflat/2+1:       ,:)
      ENDIF
#endif
#endif

    END SUBROUTINE unpack_buf_ls
    !------------------------------------------------------------------------------
    SUBROUTINE pack_ls_buf(kl)
      INTEGER, INTENT(in) :: kl

      INTEGER :: nlev, nlev0, nflat, flats(2), flate(2), n2
#ifdef __EXPLICIT
      INTEGER :: i1, i2, i3
#ifndef _OPENMP
      INTEGER :: it
#endif
#endif
#ifdef _OPENMP
      INTEGER :: it, ifls, ifle, ilat1s, ilat1e, ilat2s, ilat2e
      INTEGER :: ilat, ivar, ilat10, ilat20
      INTEGER :: icount1, icount2
#endif

      nlev   = gl_dc(idest(kl))%nflevp1
      nlev0  = gl_dc(idest(kl))%nflev
      nflat =  gl_dc(idest(kl))%nflat
      flats =  gl_dc(idest(kl))%flats
      flate =  gl_dc(idest(kl))%flate      
      n2 = nlev * nflat * nvar_lsfs

#ifdef _OPENMP
      it = omp_get_thread_num()
      ifls = thsp_flat(it)%begin
      ifle = thsp_flat(it)%end

      ilat1s = MAX(ifls,1)
      ilat1e = MIN(ifle,nflat/2)
      ilat2s = MAX(ifls,nflat/2+1)
      ilat2e = MIN(ifle,nflat)
      icount1 = ilat1e-ilat1s+1
      icount2 = ilat2e-ilat2s+1

      ilat10=flats(1)-1
      ilat20=flats(2)-1-nflat/2

      DO ivar = 1, nvar_lsfs
        DO ilat = ilat1s,ilat1e
          ls_fs(kl)%send_buffer(:,:,ilat,ivar) = ls(:,:,ilat10+ilat,ivar)
        ENDDO
        DO ilat = ilat2s,ilat2e
          ls_fs(kl)%send_buffer(:,:,ilat,ivar) = ls(:,:,ilat20+ilat,ivar)
        ENDDO
      ENDDO
      IF (lm0s(kl)) THEN
        DO ivar=1,nvar0_lsfs
          DO ilat = ilat1s,ilat1e
            ls_fs(kl)%send_buffer0(:,ilat,ivar) = ls0 (:,ilat10+ilat,ivar)
          ENDDO
          DO ilat = ilat2s,ilat2e
            ls_fs(kl)%send_buffer0(:,ilat,ivar) = ls0 (:,ilat20+ilat,ivar)
          ENDDO
        ENDDO
      ENDIF
#else
#ifdef __EXPLICIT
      CALL explicit_pack_ls_buf(ls_fs(kl)%send_buffer(1,1,1,1),ls(1,1,1,1), &
           2*nlm*nlev,SIZE(LS,3),nflat,nvar_lsfs,flats(1),flate(1),1,nflat/2)
      CALL explicit_pack_ls_buf(ls_fs(kl)%send_buffer(1,1,1,1),ls(1,1,1,1), &
           2*nlm*nlev,SIZE(LS,3),nflat,nvar_lsfs,flats(2),flate(2),nflat/2+1,nflat)
      IF (lm0s(kl)) THEN
        DO i3 = 1, nvar0_lsfs
          DO i2 = 1, nflat/2
            it = flats(1)-1+i2
            DO i1 = 1, nlev0
              ls_fs(kl)%send_buffer0 (i1,i2,i3) = ls0(i1,it,i3)
            ENDDO
          ENDDO
          DO i2 = nflat/2+1, nflat
            it = flats(2)+i2-(nflat/2+1)
            DO i1 = 1, nlev0
              ls_fs(kl)%send_buffer0 (i1,i2,i3) = ls0(i1,it,i3)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
#else
      ls_fs(kl)%send_buffer (:,:,         :nflat/2,:) = ls(:,:,flats(1):flate(1),:)
      ls_fs(kl)%send_buffer (:,:,nflat/2+1:       ,:) = ls(:,:,flats(2):flate(2),:)
      IF (lm0s(kl)) THEN
        ls_fs(kl)%send_buffer0 (:,         :nflat/2,:) = ls0 (:,flats(1):flate(1),:)
        ls_fs(kl)%send_buffer0 (:,nflat/2+1:       ,:) = ls0 (:,flats(2):flate(2),:)
      ENDIF
#endif
#endif

    END SUBROUTINE pack_ls_buf

  END SUBROUTINE transpose_fs_ls
  !==============================================================================
  SUBROUTINE transpose_ls_sp (gl_dc, sign, ls1, sp1, ls2, sp2, ls3, sp3, ls0, sp0)
    !
    ! transpose
    !   sign=foreward : Legendre space  -> spectral space
    !   sign=backward : Legendre space <-  spectral space
    !
    TYPE (pe_decomposed) ,INTENT(in)     :: gl_dc (:)     ! decomposition
    INTEGER              ,INTENT(in)     :: sign          ! foreward:ls>sp; backward:ls<sp
    REAL(dp)             ,INTENT(inout)  :: ls1   (:,:,:) ! Legendre space 
    REAL(dp)             ,INTENT(inout)  :: sp1   (:,:,:) ! spectral space
    REAL(dp)             ,INTENT(inout)  :: ls2   (:,:,:) ! Legendre space
    REAL(dp)             ,INTENT(inout)  :: sp2   (:,:,:) ! spectral space
    REAL(dp)             ,INTENT(inout)  :: ls3   (:,:,:) ! Legendre space
    REAL(dp)             ,INTENT(inout)  :: sp3   (:,:,:) ! spectral space
    REAL(dp), OPTIONAL   ,INTENT(inout)  :: ls0   (:,:)   ! Legendre (m=0 only)
    REAL(dp), OPTIONAL   ,INTENT(inout)  :: sp0   (:,:)   ! spectral (m=0 only)
    !
    ! local variables
    !
    INTEGER :: k, i, n       ! loop indices
    INTEGER :: imype         ! decomposition table index of this pe
    INTEGER :: nllevp1       ! number of levels in Legendre space
    INTEGER :: llevs         ! first level in Legendre space
    INTEGER :: lleve         ! last level in Legendre space
    INTEGER :: nlnm0         ! number of coeff. with m=0 (Legendre)
    INTEGER :: snsp          ! number of coefficients in sp. space
    INTEGER :: ssps          ! first coefficients in spectral space
    INTEGER :: sspe          ! last coefficients in spectral space
    INTEGER :: nsnm0         ! number of coeff. with m=0 (spectral)
    INTEGER :: snn0          ! first n for m=0 in spectral space
    INTEGER :: ke, nk        ! actual last level, number of levels
    INTEGER :: nprocb        ! number of PEs in set A

    INTEGER, ALLOCATABLE, SAVE :: idest(:) ! destination of data
    INTEGER, ALLOCATABLE, SAVE :: isrc(:)  ! source of data

    LOGICAL, ALLOCATABLE, SAVE :: lm0r(:)  ! receive m0
    LOGICAL, ALLOCATABLE, SAVE :: lm0s(:)  ! send m0

    !LOGICAL, SAVE :: first_call = .TRUE.
    !-------------------------------------------------------------------
    imype  = indx (p_pe, gl_dc)
    nprocb = gl_dc(imype)%nprocb
    IF (gl_dc(imype)%col_1d) RETURN
    !
    IF (.NOT. ALLOCATED(plan_b)) THEN
      ALLOCATE(plan_b(0:nprocb-1))
    ENDIF

    k = 0
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)%set_a /=  gl_dc(imype)%set_a) CYCLE
      plan_b(k) = i ! gl_dc(i)%pe
      IF (i == imype) n = k
      k = k + 1
    END DO
    plan_b = CSHIFT (plan_b,n)

    IF (.NOT. ALLOCATED(idest)) THEN
      ALLOCATE(idest(0:nprocb-1))
      ALLOCATE(isrc(0:nprocb-1))
    ENDIF

    DO k = 0, nprocb-1
      idest(k) = plan_b(k)                    ! PE index (send)
      isrc(k)  = plan_b(MOD(nprocb-k,nprocb)) ! PE index (recv)
    ENDDO

    IF (.NOT. ALLOCATED(lm0r)) THEN
       ALLOCATE(lm0r(0:nprocb-1))
       ALLOCATE(lm0s(0:nprocb-1))
    ENDIF

    IF (.NOT. ALLOCATED(ls_sp)) THEN
      ALLOCATE (ls_sp(0:nprocb-1))

      DO k = 0, nprocb-1

        !  Legendre Space -> Spectral space

        nllevp1 = gl_dc(imype)%nllevp1
        snsp    = gl_dc(idest(k))%snsp
        nsnm0   = gl_dc(idest(k))%nsnm0
        
        ALLOCATE (ls_sp(k)%send_buffer(nllevp1, 2, snsp, 3))
        ALLOCATE (ls_sp(k)%send_buffer0(nllevp1, nsnm0,1))

        ls_sp(k)%send_buffer(:,:,:,:) = 0.0_dp
        ls_sp(k)%send_buffer0(:,:,:) = 0.0_dp
        
        nllevp1 = gl_dc(isrc(k))%nllevp1
        snsp    = gl_dc(imype)%snsp
        nsnm0   = gl_dc(imype)%nsnm0
        
        ALLOCATE (ls_sp(k)%recv_buffer(nllevp1, 2, snsp, 3))
        ALLOCATE (ls_sp(k)%recv_buffer0(nllevp1, nsnm0,1))
 
      ENDDO
    ENDIF

    IF (.NOT. ALLOCATED(sp_ls)) THEN
      ALLOCATE (sp_ls(0:nprocb-1))

      DO k = 0, nprocb-1

        ! Spectral -> Legendre Space

        nllevp1 = gl_dc(idest(k))%nllevp1
        snsp    = gl_dc(imype)%snsp
        nsnm0   = gl_dc(imype)%nsnm0

        ALLOCATE (sp_ls(k)%send_buffer(nllevp1, 2, snsp, 3))
        ALLOCATE (sp_ls(k)%send_buffer0(nllevp1, nsnm0,1))

        sp_ls(k)%send_buffer(:,:,:,:) = 0.0_dp
        sp_ls(k)%send_buffer0(:,:,:) = 0.0_dp

        nllevp1 = gl_dc(imype)%nllevp1
        snsp    = gl_dc(isrc(k))%snsp
        nsnm0   = gl_dc(isrc(k))%nsnm0

        ALLOCATE (sp_ls(k)%recv_buffer(nllevp1, 2, snsp, 3))
        ALLOCATE (sp_ls(k)%recv_buffer0(nllevp1, nsnm0,1))

      ENDDO

    ENDIF
    !
    !------------------------------------------------------------------------
    !
    SELECT CASE (sign)
    CASE (foreward)
      !
      ! Legendre space -> spectral space
      !
      DO k = 0, nprocb-1

         lm0s(k)  = gl_dc(imype)%nlnm0>0 .AND. gl_dc(idest(k))%nsnm0>0.AND. PRESENT(ls0).AND. PRESENT(sp0)
         IF (imype == idest(k)) THEN
           lm0r(k) = lm0s(k)
         ELSE
           lm0r(k) =  gl_dc(isrc(k))%nlnm0>0 .AND. gl_dc(imype)%nsnm0>0.AND. PRESENT(ls0).AND. PRESENT(sp0)
         ENDIF

         nllevp1 = gl_dc(imype)%nllevp1
         llevs   = gl_dc(imype)%llevs  
         lleve   = gl_dc(imype)%lleve  
         nlnm0   = gl_dc(imype)%nlnm0  
         snsp    = gl_dc(idest(k))%snsp   
         ssps    = gl_dc(idest(k))%ssps   
         sspe    = gl_dc(idest(k))%sspe   
         nsnm0   = gl_dc(idest(k))%nsnm0  
          
         IF (lm0s(k)) snn0 = gl_dc(idest(k))%snn0(1) 

        CALL pack_ls_buf(k)

      END DO

      DO k = 0, nprocb-1

        nllevp1 = gl_dc(isrc(k))%nllevp1
        llevs   = gl_dc(isrc(k))%llevs  
        lleve   = gl_dc(isrc(k))%lleve  
        nlnm0   = gl_dc(isrc(k))%nlnm0  
        snsp    = gl_dc(imype)%snsp   
        ssps    = gl_dc(imype)%ssps   
        sspe    = gl_dc(imype)%sspe   
        nsnm0   = gl_dc(imype)%nsnm0  
    
        IF (lm0r(k)) snn0 = gl_dc(imype)%snn0(1) 

        CALL unpack_buf_sp(k)

      END DO

    CASE (backward)
      !
      ! Legendre space <- spectral space
      !
      DO k = 0, nprocb-1

         lm0s(k)  = gl_dc(idest(k))%nlnm0>0 .AND. gl_dc(imype)%nsnm0>0.AND. PRESENT(ls0).AND. PRESENT(sp0)
         IF (imype == idest(k)) THEN
           lm0r(k) = lm0s(k)
         ELSE
           lm0r(k) =  gl_dc(imype)%nlnm0>0 .AND. gl_dc(isrc(k))%nsnm0>0.AND. PRESENT(ls0).AND. PRESENT(sp0)
         ENDIF

        nllevp1 = gl_dc(idest(k))%nllevp1
        llevs   = gl_dc(idest(k))%llevs  
        lleve   = gl_dc(idest(k))%lleve  
        nlnm0   = gl_dc(idest(k))%nlnm0  
        snsp    = gl_dc(imype)%snsp   
        ssps    = gl_dc(imype)%ssps   
        sspe    = gl_dc(imype)%sspe   
        nsnm0   = gl_dc(imype)%nsnm0  
    
        IF(lm0s(k)) snn0 = gl_dc(imype)%snn0(1) 

        CALL pack_sp_buf(k)

      END DO

      DO k = 0, nprocb-1

        nllevp1 = gl_dc(imype)%nllevp1
        llevs   = gl_dc(imype)%llevs  
        lleve   = gl_dc(imype)%lleve  
        nlnm0   = gl_dc(imype)%nlnm0  
        snsp    = gl_dc(isrc(k))%snsp   
        ssps    = gl_dc(isrc(k))%ssps   
        sspe    = gl_dc(isrc(k))%sspe   
        nsnm0   = gl_dc(isrc(k))%nsnm0  
    
        IF(lm0r(k)) snn0 = gl_dc(isrc(k))%snn0(1) 

        CALL unpack_buf_ls(k)

      END DO

    CASE default

      CALL finish ('transpose_ls_sp', 'invalid transposition parameter (not foreward,backward)')

    END SELECT

    !------------------------------------------------------------------------------
  CONTAINS
    !------------------------------------------------------------------------------
    SUBROUTINE pack_ls_buf(kl)
      INTEGER, INTENT(in) :: kl

#ifdef __EXPLICIT
      INTEGER :: i1, i2, i3, it
!$OMP PARALLEL PRIVATE(i1,i2,i3,it)
!$OMP DO
      DO i3 = 1, snsp
        it = ssps-1+i3
        DO i2 = 1, 2
          DO i1 = 1, SIZE(ls1,1)
            ls_sp(kl)%send_buffer(i1,i2,i3,1) = ls1(i1,i2,it)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO
      DO i3 = 1, snsp
        it = ssps-1+i3
        DO i2 = 1, 2
          DO i1 = 1, SIZE(ls2,1)
            ls_sp(kl)%send_buffer(i1,i2,i3,2) = ls2(i1,i2,it)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO
      DO i3 = 1, snsp
        it = ssps-1+i3
        DO i2 = 1, 2
          DO i1 = 1, SIZE(ls3,1)
            ls_sp(kl)%send_buffer(i1,i2,i3,3) = ls3(i1,i2,it)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
      IF (lm0s(kl)) THEN
!$OMP DO
        DO i2 = 1, nsnm0
          it = snn0+i2
          DO i1 = 1, SIZE(ls0,1)
            ls_sp(kl)%send_buffer0 (i1,i2,1) = ls0(i1,it)
          ENDDO
        ENDDO 
!$OMP END DO NOWAIT
      ENDIF
!$OMP END PARALLEL
#else
!$OMP PARALLEL
!$OMP WORKSHARE
      ls_sp(kl)%send_buffer (:SIZE(ls1,1),:,:,1)  = ls1 (:,:,ssps:sspe)
      ls_sp(kl)%send_buffer (:SIZE(ls2,1),:,:,2)  = ls2 (:,:,ssps:sspe)
      ls_sp(kl)%send_buffer (:SIZE(ls3,1),:,:,3)  = ls3 (:,:,ssps:sspe)
!$OMP END WORKSHARE
      IF (lm0s(kl)) THEN
!$OMP WORKSHARE
        ls_sp(kl)%send_buffer0 (:SIZE(ls0,1),:,1) = ls0 (:,snn0+1:snn0+nsnm0)
!$OMP END WORKSHARE
      ENDIF
!$OMP END PARALLEL
#endif
    END SUBROUTINE pack_ls_buf
    !------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_sp(kl)
      INTEGER, INTENT(in) :: kl
!$OMP PARALLEL
      ke = MIN(SIZE(sp1,1), lleve)
      nk = ke-llevs+1
      IF (nk > 0) THEN
!$OMP WORKSHARE
        sp1 (llevs:ke, :, :) = ls_sp(kl)%recv_buffer (:nk,:,:,1)
!$OMP END WORKSHARE
      ENDIF
      ke = MIN(SIZE(sp2,1), lleve)
      nk = ke-llevs+1
      IF (nk > 0) THEN
!$OMP WORKSHARE
        sp2 (llevs:ke, :, :) = ls_sp(kl)%recv_buffer (:nk,:,:,2)
!$OMP END WORKSHARE
      ENDIF
      ke = MIN(SIZE(sp3,1), lleve)
      nk = ke-llevs+1
      IF (nk > 0) THEN
!$OMP WORKSHARE
        sp3 (llevs:ke, :, :) = ls_sp(kl)%recv_buffer (:nk,:,:,3)
!$OMP END WORKSHARE
      ENDIF
      IF (lm0r(kl)) THEN
        ke = MIN(SIZE(sp0,1), lleve)
        nk = ke-llevs+1
        IF(nk > 0) THEN
!$OMP WORKSHARE
          sp0 (llevs:ke,:) = ls_sp(kl)%recv_buffer0 (:nk,:,1)
!$OMP END WORKSHARE
        ENDIF
      ENDIF
!$OMP END PARALLEL
    END SUBROUTINE unpack_buf_sp
    !------------------------------------------------------------------------------
    SUBROUTINE pack_sp_buf(kl)
      INTEGER, INTENT(in) :: kl
#ifdef __EXPLICIT
      INTEGER :: i1, i2, i3, it
!$OMP PARALLEL PRIVATE(i1,i2,i3,it)
      ke = MIN(SIZE(sp1,1), lleve)
      nk = MAX(0, ke-llevs+1)
!$OMP DO
      DO i3 = 1, snsp
        DO i2 = 1, 2
          DO i1 = 1, nk
            it = llevs-1+i1
            sp_ls(kl)%send_buffer (i1,i2,i3,1) = sp1(it,i2,i3) 
          ENDDO
          DO i1 = nk+1, nllevp1
            sp_ls(kl)%send_buffer (i1,i2,i3,1) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
      ke = MIN(SIZE(sp2,1), lleve) 
      nk = MAX(0, ke-llevs+1)
!$OMP DO
      DO i3 = 1, snsp
        DO i2 = 1, 2
          DO i1 = 1, nk
            it = llevs-1+i1
            sp_ls(kl)%send_buffer (i1,i2,i3,2) = sp2(it,i2,i3) 
          ENDDO
          DO i1 = nk+1, nllevp1
            sp_ls(kl)%send_buffer (i1,i2,i3,2) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
      ke = MIN(SIZE(sp3,1), lleve)
      nk = MAX(0, ke-llevs+1)
!$OMP DO
      DO i3 = 1, snsp
        DO i2 = 1, 2
          DO i1 = 1, nk
            it = llevs-1+i1
            sp_ls(kl)%send_buffer (i1,i2,i3,3) = sp3(it,i2,i3) 
          ENDDO
          DO i1 = nk+1, nllevp1
            sp_ls(kl)%send_buffer (i1,i2,i3,3) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
      IF (lm0s(kl)) THEN
        ke = MIN(SIZE(sp0,1), lleve)
        nk = MAX(0, ke-llevs+1)
!$OMP DO
        DO i2 = 1, nsnm0
          DO i1 = 1, nk
            it = llevs-1+i1
            sp_ls(kl)%send_buffer0 (i1,i2,1) = sp0 (it,i2)
          ENDDO
          DO i1 = nk+1, nllevp1
            sp_ls(kl)%send_buffer0 (i1,i2,:) = 0.0_dp
          ENDDO
        ENDDO
!$OMP ENDDO
      ENDIF
!$OMP END PARALLEL
#else
!$OMP PARALLEL
      ke = MIN(SIZE(sp1,1), lleve)
      nk = MAX(0, ke-llevs+1)
!$OMP WORKSHARE
      sp_ls(kl)%send_buffer (:nk,:,:,1) = sp1 (llevs:ke, :, :) 
      sp_ls(kl)%send_buffer (nk+1:,:,:,1) = 0.0_dp
!$OMP END WORKSHARE
      ke = MIN(SIZE(sp2,1), lleve) 
      nk = MAX(0, ke-llevs+1)
!$OMP WORKSHARE
      sp_ls(kl)%send_buffer (:nk,:,:,2) = sp2 (llevs:ke, :, :) 
      sp_ls(kl)%send_buffer (nk+1:,:,:,2) = 0.0_dp
!$OMP END WORKSHARE
      ke = MIN(SIZE(sp3,1), lleve)
      nk = MAX(0, ke-llevs+1)
!$OMP WORKSHARE
      sp_ls(kl)%send_buffer (:nk,:,:,3) = sp3 (llevs:ke, :, :) 
      sp_ls(kl)%send_buffer (nk+1:,:,:,3) = 0.0_dp
!$OMP END WORKSHARE
      IF (lm0s(kl)) THEN
        ke = MIN(SIZE(sp0,1), lleve)
        nk = MAX(0, ke-llevs+1)
!$OMP WORKSHARE
        sp_ls(kl)%send_buffer0 (:nk,:,1) = sp0 (llevs:ke,:)
        sp_ls(kl)%send_buffer0 (nk+1:,:,:) = 0.0_dp
!$OMP END WORKSHARE
      ENDIF
!$OMP END PARALLEL
#endif
    END SUBROUTINE pack_sp_buf
    !------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_ls(kl)
      INTEGER, INTENT(in) :: kl
      INTEGER :: sls0, sls1, sls2, sls3
      sls1 = SIZE(ls1,1) 
      sls2 = SIZE(ls2,1) 
      sls3 = SIZE(ls3,1) 
!$OMP PARALLEL
!$OMP WORKSHARE
      ls1 (:,:,ssps:sspe) = sp_ls(kl)%recv_buffer (:sls1,:,:,1)
      ls2 (:,:,ssps:sspe) = sp_ls(kl)%recv_buffer (:sls2,:,:,2)
      ls3 (:,:,ssps:sspe) = sp_ls(kl)%recv_buffer (:sls3,:,:,3)
!$OMP END WORKSHARE
      IF (lm0r(kl)) THEN
        sls0 = SIZE(ls0,1) 
!$OMP WORKSHARE
        ls0 (:,snn0+1:snn0+nsnm0) = sp_ls(kl)%recv_buffer0 (:sls0,:,1)
!$OMP END WORKSHARE
      ENDIF
!$OMP END PARALLEL
    END SUBROUTINE unpack_buf_ls
    !------------------------------------------------------------------------------
   SUBROUTINE sendrecv_lssp(trbuf)
      !
      ! send and receive buffer
      ! deallocate send buffer
      !
      TYPE(transpose_buffer), INTENT(inout) :: trbuf(0:)
      INTEGER :: kl
      !
#ifdef __SENDRECV
      INTEGER :: communication_type = SENDRECV
#else
      INTEGER :: communication_type = NON_BLOCKING
#endif
      !
      SELECT CASE(communication_type)
      CASE (sendrecv)
        DO kl = 0, nprocb-1
          IF(imype /= idest(kl)) THEN
            CALL p_sendrecv (trbuf(kl)%send_buffer, gl_dc(idest(kl))%pe, &
                             trbuf(kl)%recv_buffer, gl_dc(isrc(kl))%pe,  &
                             tag_tr_ls_sp)
            IF (lm0s(kl) .AND. lm0r(kl)) THEN
              CALL p_sendrecv (trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                               trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe,  &
                               tag_tr_ls_sp)
            ELSE IF (lm0s(kl)) THEN
              CALL p_send (trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                           tag_tr_ls_sp)
            ELSE IF (lm0r(kl)) THEN
              CALL p_recv (trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe, &
                           tag_tr_ls_sp)
            ENDIF
          ELSE
            trbuf(kl)%recv_buffer = trbuf(kl)%send_buffer
            IF (lm0r(kl)) trbuf(kl)%recv_buffer0 = trbuf(kl)%send_buffer0
          ENDIF
        ENDDO
      CASE (non_blocking)
        DO kl = 0, nprocb-1
          IF(imype /= idest(kl)) THEN
            CALL p_isend (trbuf(kl)%send_buffer, gl_dc(idest(kl))%pe, &
                          tag_tr_ls_sp)
            CALL p_irecv (trbuf(kl)%recv_buffer, gl_dc(isrc(kl))%pe,  &
                          tag_tr_ls_sp)
            IF (lm0s(kl) .AND. lm0r(kl)) THEN
              CALL p_isend(trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                           tag_tr_ls_sp)
              CALL p_irecv (trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe,  &
                            tag_tr_ls_sp)
            ELSE IF (lm0s(kl)) THEN
              CALL p_isend (trbuf(kl)%send_buffer0, gl_dc(idest(kl))%pe, &
                            tag_tr_ls_sp)
            ELSE IF (lm0r(kl)) THEN
              CALL p_irecv (trbuf(kl)%recv_buffer0, gl_dc(isrc(kl))%pe, &
                            tag_tr_ls_sp)
            ENDIF
          ELSE
            trbuf(kl)%recv_buffer = trbuf(kl)%send_buffer
            IF (lm0r(kl)) trbuf(kl)%recv_buffer0 = trbuf(kl)%send_buffer0
          ENDIF
        ENDDO
        CALL p_wait
      END SELECT
    END SUBROUTINE sendrecv_lssp
    !------------------------------------------------------------------------------
  END SUBROUTINE transpose_ls_sp
  !==============================================================================
  SUBROUTINE transpose_gp_ffsl_2 (dc, sign, x_gp, x_ffsl)
    TYPE(pe_decomposed),INTENT(in)    :: dc           ! decomposition
    INTEGER            ,INTENT(in)    :: sign         ! direction
    REAL(dp), TARGET   ,INTENT(inout) :: x_gp   (:,:) ! grid point field
    REAL(dp)           ,INTENT(inout) :: x_ffsl (:,:) ! ffsl decomposition
    !
    ! transpose to/from decomposition required by the ffsl
    !   (Flux-Form Semi-Lagrangian) transport scheme
    !
    !   ordering required by ffsl:
    !     Latitudes are sorted South to North
    !     No splitting of hemispheres
    !
    !   Task of this routine:
    !     sign=foreward : Gridpoint space  -> ffsl      space
    !     sign=backward : ffsl      space  -> Gridpoint space
    !
    !     Andreas Rhodin, DWD/MPI, Aug 2001
    !
    REAL(dp), ALLOCATABLE :: buf   (:,:) ! receive buffer
    REAL(dp)     ,POINTER :: lx_gp (:,:) ! pointer/temporary buffer

    LOGICAL :: lreg        ! flag for regular grid
    INTEGER :: nlat2, nlatx, nlatf, nlatg, nlong
    INTEGER :: n, idx

    IF (require_init) CALL init_transpose
    !
    !-------------------------------------------------------------------
    nlong = dc%      nglon   !      number of gp   longitudes on this  pe
    nlatg = dc%      nglat   !      number of gp   latitudes  on this  pe
    nlat2 = dc%      nglat/2 ! half number of gp   latitudes  on this  pe
    nlatx = dc%ffsl%nlatx   ! half number of gp   latitudes  on other pe
    nlatf = dc%ffsl%nlat    !      number of ffsl latitudes  on this  pe
    lreg  = dc%lreg
    !
    ! gridpoint -> ffsl
    !
    IF (sign == foreward) THEN
      !
      ! blocking (nproma)
      !
      IF (lreg) THEN
        lx_gp => x_gp
      ELSE
        ALLOCATE (lx_gp(dc%nglon,dc%nglat))
        CALL reorder (lx_gp,x_gp)
      ENDIF
      !
      ! nprocb == 1
      !
      IF(dc%nprocb==1) THEN
        IF (dc%pe == dc%ffsl%pe_x) THEN          ! just reorder
          x_ffsl (:,1:nlatf) = lx_gp (:,nlatg:1:-1)
        ELSE                                        ! send / receive
          ALLOCATE (buf (dc%nlon, dc%ffsl%nlatx))
          IF (dc%ffsl%latn == dc%glate(2)) THEN  ! exchange southern hemisp.
            CALL p_sendrecv                                    &
              (lx_gp(:,nlat2+1:), dc%ffsl%pe_x,               &
                buf,              dc%ffsl%pe_x, tag_tr_gp_ffsl)
            x_ffsl (:,:nlatx   ) = buf   (:,nlatx:1:-1)
            x_ffsl (:, nlatx+1:) = lx_gp (:,nlat2:1:-1)
          ELSE                                      ! exchange northern hemisp.
            CALL p_sendrecv                                  & 
              (lx_gp(:,:nlat2), dc%ffsl%pe_x,               &
               buf,             dc%ffsl%pe_x, tag_tr_gp_ffsl)
            x_ffsl (:,:nlatx   ) = buf   (:,nlatx:      1:-1)
            x_ffsl (:, nlatx+1:) = lx_gp (:,nlatg:nlat2+1:-1)
          ENDIF
          DEALLOCATE (buf)
        ENDIF
        !
        ! nprocb /= 1
        !
      ELSE
        ALLOCATE (buf (dc%nglon, dc%nglat))
        buf(:,:) = lx_gp(:,dc%nglat:1:-1)
        !
        ! Distribute our part to ALL PEs with identical set_a
        ! This is kind of a broadcast, data is replicated on all PEs with 
        ! identical set_a. send:
        !
        DO n=1,dc%nprocb
          !
          ! For PEs with odd set_a:
          !   Northern Hemisphere sent to same set_a
          !   Southern Hemisphere sent to set_a+1 (unless set_a == nproca)
          ! For PEs with even set_a:
          !   Northern Hemisphere sent to set_a-1
          !   Southern Hemisphere sent to same set_a
          !
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_isend(buf(1,nlat2+1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
            idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
            CALL p_isend(buf(1,1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
          ELSE
            idx = dc%mapmesh(n,dc%set_a-1)
            CALL p_isend(buf(1,nlat2+1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_isend(buf(1,1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
          ENDIF
        ENDDO
        !
        ! recv:
        !
        DO n=1,dc%nprocb
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),nlatx+1:), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
            idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
            CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),1:nlatx), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
          ELSE
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),nlatx+1:), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
            idx = dc%mapmesh(n,dc%set_a-1)
            CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),1:nlatx), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
          ENDIF
        ENDDO
        !
        ! wait, deallocate
        !
        CALL p_wait
        DEALLOCATE (buf)
        !
      ENDIF
      IF (.NOT.lreg) DEALLOCATE (lx_gp)
      !
      ! ffsl -> gridpoint
      !
    ELSE
      !
      ! blocking (nproma)
      !
      IF (lreg) THEN
        lx_gp => x_gp
      ELSE
        ALLOCATE (lx_gp(dc%nglon,dc%nglat))
      ENDIF
      !
      ! nprocb == 1
      !
      IF(dc%nprocb==1) THEN
        IF (dc%pe == dc%ffsl%pe_x) THEN          ! just reorder
          lx_gp (:,nlatg:1:-1) = x_ffsl (:,1:nlatf)
        ELSE                                        ! send / receive
          ALLOCATE (buf (dc%nlon, dc%ffsl%nlatx))
          IF (dc%ffsl%latn == dc%glate(2)) THEN  ! exchange southern hemisp.
            buf   (:,nlatx:1:-1) = x_ffsl (:,:nlatx   )
            lx_gp (:,nlat2:1:-1) = x_ffsl (:, nlatx+1:)
            CALL p_sendrecv                                    &
              (buf,               dc%ffsl%pe_x,               &
               lx_gp(:,nlat2+1:), dc%ffsl%pe_x, tag_tr_gp_ffsl)
          ELSE                                      ! exchange northern hemisp.
            lx_gp (:,nlatg:nlat2+1:-1) = x_ffsl (:, nlatx+1:)
            buf   (:,nlatx:      1:-1) = x_ffsl (:,:nlatx   )
            CALL p_sendrecv                                  &
              (buf,             dc%ffsl%pe_x,               &
               lx_gp(:,:nlat2), dc%ffsl%pe_x, tag_tr_gp_ffsl)
          ENDIF
          DEALLOCATE (buf)
        ENDIF
        !
        ! nproca /= 1
        !
      ELSE
        ! 
        ! This is exactly like gp -> ffsl with sends/recvs exchanged 
        ! and the reordering step put to the end
        !
        ALLOCATE (buf (dc%nglon, dc%nglat))
        !
        DO n = 1, dc%nprocb
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_irecv(buf(1,nlat2+1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
            idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))            
            CALL p_irecv(buf(1,1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
          ELSE
            idx = dc%mapmesh(n,dc%set_a-1)
            CALL p_irecv(buf(1,nlat2+1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_irecv(buf(1,1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
          ENDIF
        ENDDO
        !
        DO n=1,dc%nprocb
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_send(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),nlatx+1:), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
            idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
            CALL p_send(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),1:nlatx), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
          ELSE
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_send(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),nlatx+1:), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
            idx = dc%mapmesh(n,dc%set_a-1)
            CALL p_send(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),1:nlatx), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
          ENDIF
        ENDDO
        !
        CALL p_wait
        lx_gp(:,dc%nglat:1:-1) = buf(:,:)
        DEALLOCATE (buf)
      ENDIF
      !
      ! blocking (nproma)
      !
      IF (.NOT.lreg) THEN
        CALL reorder (x_gp,lx_gp)
        DEALLOCATE (lx_gp)
      ENDIF
    ENDIF
  END SUBROUTINE transpose_gp_ffsl_2
  !------------------------------------------------------------------------------
  SUBROUTINE transpose_gp_ffsl_3 (dc, sign, x_gp, x_ffsl)
    TYPE(pe_decomposed),INTENT(in)    :: dc             ! decomposition
    INTEGER            ,INTENT(in)    :: sign           ! direction
    REAL(dp), TARGET   ,INTENT(inout) :: x_gp   (:,:,:) ! grid point field
    REAL(dp)           ,INTENT(inout) :: x_ffsl (:,:,:) ! ffsl decomposition
    !
    ! transpose to/from decomposition required by the ffsl
    !   (Flux-Form Semi-Lagrangian) transport scheme
    !
    !   ordering required by ffsl:
    !     Latitudes are sorted South to North
    !     No splitting of hemispheres
    !
    !   Task of this routine:
    !     sign=foreward : Gridpoint space  -> ffsl      space
    !     sign=backward : ffsl      space  -> Gridpoint space
    !
    !     Andreas Rhodin, DWD/MPI, Aug 2001
    !
    REAL(dp), ALLOCATABLE :: buf (:,:,:)   ! receive buffer
    REAL(dp)     ,POINTER :: lx_gp (:,:,:) ! pointer/temporary buffer
    
    LOGICAL :: lreg        ! flag for regular grid
    INTEGER :: nlat2, nlatx, nlatf, nlatg, nlong
    INTEGER :: j, k, k2, n, idx, dest_set_b, src_set_b
    
    !-------------------------------------------------------------------
    ! additional barrier to catch load imbalance
    !
    IF (require_init) CALL init_transpose
    !
    !-------------------------------------------------------------------
    nlong = dc%      nglon   !      number of gp   longitudes on this  pe
    nlatg = dc%      nglat   !      number of gp   latitudes  on this  pe
    nlat2 = dc%      nglat/2 ! half number of gp   latitudes  on this  pe
    nlatx = dc%ffsl%nlatx   ! half number of gp   latitudes  on other pe
    nlatf = dc%ffsl%nlat    !      number of ffsl latitudes  on this  pe
    lreg  = dc%lreg
    !
    ! gridpoint -> ffsl
    !
    IF (sign==foreward) THEN
      !
      ! blocking (nproma)
      !
      IF (lreg) THEN
        lx_gp => x_gp
      ELSE
        ALLOCATE (lx_gp(dc%nglon,dc%nlev,dc%nglat))
        !PRINT*,'tr_gp_ffsl_3:lreg=',lreg
        !PRINT*,'size(x_gp)=',SIZE(x_gp,1),SIZE(x_gp,2),SIZE(x_gp,3)
        !PRINT*,'dc%nglon, dc%nlev, dc%nglat=',dc%nglon, dc%nlev, dc%nglat
        CALL reorder (lx_gp,x_gp)
      ENDIF
      !
      ! nprocb == 1
      !
      IF(dc%nprocb==1) THEN
        IF (dc%pe == dc%ffsl%pe_x) THEN          ! just reorder
          DO j=1,nlatf
            x_ffsl (:,j,:) = lx_gp (:,:,nlatg+1-j)
          END DO
        ELSE                                        ! send / receive
          ALLOCATE (buf (dc%nlon, dc%nlev, dc%ffsl%nlatx))
          IF (dc%ffsl%latn == dc%glate(2)) THEN  ! exchange southern hemisp.
            CALL p_sendrecv                                      &
              (lx_gp(:,:,nlat2+1:), dc%ffsl%pe_x,               &
               buf,                 dc%ffsl%pe_x, tag_tr_gp_ffsl)
            DO j=1,nlatx
              x_ffsl (:,      j,:) = buf   (:,:,nlatx+1-j)
            END DO
            DO j=1,nlat2
              x_ffsl (:,nlatx+j,:) = lx_gp (:,:,nlat2+1-j)
            END DO
          ELSE                                      ! exchange northern hemisp.
            CALL p_sendrecv                                    & 
              (lx_gp(:,:,:nlat2), dc%ffsl%pe_x,               &
               buf,               dc%ffsl%pe_x, tag_tr_gp_ffsl)
            DO j=1,nlatx
              x_ffsl (:,       j,:) = buf   (:,:,nlatx+1-j)
            END DO
            DO j=1,nlat2
              x_ffsl (:, nlatx+j,:) = lx_gp (:,:,nlatg+1-j)
            END DO
          ENDIF
          DEALLOCATE (buf)
        ENDIF
      !
      ! nprocb /= 1
      !
      ELSE
        ALLOCATE (buf (dc%nglon, dc%nglat, dc%nlev))
        !
        ! send:
        !
        DO k=1,dc%nlev
          dest_set_b = MOD(k-1,dc%nprocb)+1 ! set_b of receiving PE for this k
          DO j=1,dc%nglat
            buf(:,j,k) = lx_gp(:,k,dc%nglat+1-j)
          ENDDO
          !
          ! For PEs with odd set_a:
          !   Northern Hemisphere sent to same set_a
          !   Southern Hemisphere sent to set_a+1 (unless set_a == nproca)
          ! For PEs with even set_a:
          !   Northern Hemisphere sent to set_a-1
          !   Southern Hemisphere sent to same set_a
          !
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(dest_set_b,dc%set_a)
            !PRINT*,'(a):k,dest_set_b,dc%set_a,idx=',k,dest_set_b,dc%set_a,idx

            CALL p_isend(buf(1,nlat2+1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
            idx = dc%mapmesh(dest_set_b,MIN(dc%set_a+1,dc%nproca))
            CALL p_isend(buf(1,1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
          ELSE
            idx = dc%mapmesh(dest_set_b,dc%set_a-1)
            !PRINT*,'(b):k,dest_set_b,dc%set_a-1,idx=',k,dest_set_b,dc%set_a-1,idx
            CALL p_isend(buf(1,nlat2+1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
            idx = dc%mapmesh(dest_set_b,dc%set_a)
            CALL p_isend(buf(1,1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
          ENDIF
        ENDDO
        !
        ! recv:
        !
        DO k=1,(dc%nlev-dc%set_b)/dc%nprocb+1 !dc%set_b,dc%nlev,dc%nprocb
          k2 = (k-1)*dc%nprocb+dc%set_b
          DO n=1,dc%nprocb
            IF( MOD(dc%set_a,2) == 1 ) THEN
              idx = dc%mapmesh(n,dc%set_a)
              CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),&
                          nlatx+1:,k), gdc(idx)%pe,1000+k2)
              idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
              CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),&
                          1:nlatx,k), gdc(idx)%pe,1000+k2)
            ELSE
              idx = dc%mapmesh(n,dc%set_a)
              CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),&
                          nlatx+1:,k), gdc(idx)%pe,1000+k2)
              idx = dc%mapmesh(n,dc%set_a-1)
              CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),&
                          1:nlatx,k), gdc(idx)%pe,1000+k2)
            ENDIF
          ENDDO
        ENDDO
        !
        ! wait, deallocate
        !
        CALL p_wait
        DEALLOCATE (buf)
      ENDIF
      IF (.NOT.lreg) DEALLOCATE (lx_gp)
    !
    ! ffsl -> gridpoint
    !
    ELSE
      !
      ! blocking (nproma)
      !
      IF (lreg) THEN
        lx_gp => x_gp
      ELSE
        ALLOCATE (lx_gp(dc%nglon,dc%nlev,dc%nglat))
      ENDIF
      !
      ! nprocb == 1
      !
      IF(dc%nprocb==1) THEN

        IF (dc%pe == dc%ffsl%pe_x) THEN          ! just reorder
          DO j=1,nlatf
            lx_gp (:,:,nlatg+1-j) = x_ffsl (:,j,:)
          END DO
        ELSE                                        ! send / receive
          ALLOCATE (buf (dc%nlon, dc%nlev, dc%ffsl%nlatx))
          IF (dc%ffsl%latn == dc%glate(2)) THEN  ! exchange southern hemisp.
            DO j=1,nlatx
              buf   (:,:,nlatx+1-j) = x_ffsl (:,      j,:)
            END DO
            DO j=1,nlat2
              lx_gp (:,:,nlat2+1-j) = x_ffsl (:,nlatx+j,:)
            END DO
            CALL p_sendrecv                                      &
              (buf,                 dc%ffsl%pe_x,               &
               lx_gp(:,:,nlat2+1:), dc%ffsl%pe_x, tag_tr_gp_ffsl)
          ELSE                                      ! exchange northern hemisp.
            DO j=1,nlat2
              lx_gp (:,:,nlatg+1-j) = x_ffsl (:,nlatx+j,:)
            END DO
            DO j=1,nlatx
              buf   (:,:,nlatx+1-j) = x_ffsl (:,      j,:)
            END DO
            CALL p_sendrecv                                    &
              (buf,               dc%ffsl%pe_x,               &
               lx_gp(:,:,:nlat2), dc%ffsl%pe_x, tag_tr_gp_ffsl)
          ENDIF
          DEALLOCATE (buf)
        ENDIF
      !
      ! nprocb /= 1
      !
      ELSE
        ! 
        ! This is exactly like gp -> ffsl with sends/recvs exchanged 
        ! and the reordering step put to the end
        !
        ALLOCATE (buf (dc%nglon, dc%nglat, dc%nlev))
        !
        ! irecv:
        !
        DO k=1,dc%nlev
          src_set_b = MOD(k-1,dc%nprocb)+1 ! set_b of receiving PE for this k
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(src_set_b,dc%set_a)
            CALL p_irecv(buf(1,nlat2+1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
            idx = dc%mapmesh(src_set_b,MIN(dc%set_a+1,dc%nproca))
            CALL p_irecv(buf(1,1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
          ELSE
            idx = dc%mapmesh(src_set_b,dc%set_a-1)
            CALL p_irecv(buf(1,nlat2+1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
            idx = dc%mapmesh(src_set_b,dc%set_a)
            CALL p_irecv(buf(1,1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
          ENDIF
        ENDDO
        !
        ! send:
        !
        DO k=1,(dc%nlev-dc%set_b)/dc%nprocb+1 !dc%set_b,dc%nlev,dc%nprocb
          k2 = (k-1)*dc%nprocb+dc%set_b
          DO n=1,dc%nprocb
            IF( MOD(dc%set_a,2) == 1 ) THEN
              idx = dc%mapmesh(n,dc%set_a)
              CALL p_send(x_ffsl(gdc(idx)%glons(1):             &
                                 gdc(idx)%glone(1),nlatx+1:,k), &
                          gdc(idx)%pe,1000+k2)
              idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
              CALL p_send(x_ffsl(gdc(idx)%glons(1):            &
                                 gdc(idx)%glone(1),1:nlatx,k), &
                          gdc(idx)%pe,1000+k2)
            ELSE
              idx = dc%mapmesh(n,dc%set_a)
              CALL p_send(x_ffsl(gdc(idx)%glons(1):             &
                                 gdc(idx)%glone(1),nlatx+1:,k), &
                          gdc(idx)%pe,1000+k2)
              idx = dc%mapmesh(n,dc%set_a-1)
              CALL p_send(x_ffsl(gdc(idx)%glons(1):            &
                                 gdc(idx)%glone(1),1:nlatx,k), &
                          gdc(idx)%pe,1000+k2)
            ENDIF
          ENDDO
        ENDDO
        !
        ! wait, reorder, deallocate
        !
        CALL p_wait
        DO k=1,dc%nlev
          DO j=1,dc%nglat
            lx_gp(:,k,dc%nglat+1-j) = buf(:,j,k)
          ENDDO
        ENDDO
        DEALLOCATE (buf)
      ENDIF
      !
      ! blocking (nproma)
      !
      IF (.NOT.lreg) THEN
        CALL reorder (x_gp,lx_gp)
        DEALLOCATE (lx_gp)
      ENDIF
    ENDIF
    !
  END SUBROUTINE transpose_gp_ffsl_3
  !------------------------------------------------------------------------------
  SUBROUTINE transpose_gp_ffsl_4 (dc, sign, x_gp, x_ffsl)
    TYPE(pe_decomposed), INTENT(in)    :: dc               ! decomposition
    INTEGER,             INTENT(in)    :: sign             ! direction
    REAL(dp),            INTENT(inout) :: x_gp   (:,:,:,:) ! grid point field
    REAL(dp),            INTENT(inout) :: x_ffsl (:,:,:,:) ! ffsl decomposition
    !
    ! transpose to/from decomposition required by the ffsl
    !   (Flux-Form Semi-Lagrangian) transport scheme
    !
    !   ordering required by ffsl:
    !     Latitudes are sorted South to North
    !     No splitting of hemispheres
    !
    !   Task of this routine:
    !     sign=foreward : Gridpoint space  -> ffsl      space
    !     sign=backward : ffsl      space  -> Gridpoint space
    !
    !     Andreas Rhodin, DWD/MPI, Aug 2001
    !
    INTEGER :: jt
    IF (require_init) CALL init_transpose
    DO jt = 1, SIZE(x_ffsl,4)
      CALL transpose_gp_ffsl_3(dc, sign, x_gp(:,:,jt,:), x_ffsl(:,:,:,jt))
    END DO
    !
  END SUBROUTINE transpose_gp_ffsl_4
  !==============================================================================
  FUNCTION indx0 (pe, gl_dc)
    INTEGER              ,INTENT(in) :: pe       ! processor id
    TYPE (pe_decomposed) ,INTENT(in) :: gl_dc(:) ! global decomposition
    INTEGER                          :: indx0    ! index
    !
    ! returns the index of a given PE in the global decomposition table
    !
    INTEGER :: i
    indx0 = -1
    DO i = 1, SIZE(gl_dc)
      IF(gl_dc(i)%pe == pe) THEN
        indx0 = i
        RETURN
      ENDIF
    END DO
    WRITE (message_text,*) 'index not found in decomposition table'
    CALL message ('', TRIM(message_text))
    WRITE (message_text,*) '  required:',pe
    CALL message ('', TRIM(message_text))
    WRITE (message_text,*) '  found  :',gl_dc%pe
    CALL message ('', TRIM(message_text))
    CALL finish('mo_transpose:indx', 'index not found in decomposition table')
  END FUNCTION indx0
  !------------------------------------------------------------------------------
  FUNCTION indx2 (pe, gl_dc)
    INTEGER              ,INTENT(in) :: pe(:,:)       ! processor id
    TYPE (pe_decomposed) ,INTENT(in) :: gl_dc(:)      ! global decomposition
    INTEGER                          :: indx2(SIZE(pe,1),SIZE(pe,2))     ! index
    !
    ! returns the index of a given PE in the global decomposition table
    !
    INTEGER :: i
    indx2 = -1
    DO i = 1, SIZE(gl_dc)
      WHERE(gl_dc(i)%pe == pe)
        indx2 = i
      endwhere
    END DO
    IF(ANY(indx2==-1))THEN
      WRITE (message_text,*) 'index not found in decomposition table'
      CALL message ('', TRIM(message_text))
      WRITE (message_text,*) '  required:',PACK(pe,mask=(indx2==-1))
      CALL message ('', TRIM(message_text))
      WRITE (message_text,*) '  found  :',gl_dc%pe
      CALL message ('', TRIM(message_text))
      CALL finish('mo_transpose:indx', &
           'index not found in decomposition table')
    END IF
  END FUNCTION indx2
  !==============================================================================
  FUNCTION gen_inv_map (gl_dc)
    TYPE (pe_decomposed), INTENT(in) :: gl_dc(:)
    INTEGER :: gen_inv_map(0:SIZE(gl_dc)-1) 
    INTEGER :: i
    gen_inv_map = -1
    DO i = 1, SIZE(gl_dc)
      gen_inv_map(gl_dc(i)%pe) = i
    END DO
    IF(ANY(gen_inv_map==-1))THEN
      WRITE (message_text,*) 'index not found in decomposition table'
      CALL message ('', TRIM(message_text))
      WRITE (message_text,*) '  required:',i
      CALL message ('', TRIM(message_text))
      WRITE (message_text,*) '  found  :',gl_dc%pe
      CALL message ('', TRIM(message_text))
      CALL finish('mo_transpose:gen_inv_map', &
           'index not found in decomposition table')
    END IF
  END FUNCTION gen_inv_map
  !==============================================================================
  SUBROUTINE reorder21 (y,x)
    REAL(dp) ,INTENT(out) :: y (:)
    REAL(dp) ,INTENT(in)  :: x (:,:)

#ifdef __REPLACE_RESHAPE
    CALL util_reshape2(y, x, SIZE(y), SIZE(x,1)*SIZE(x,2))
#else
    y = RESHAPE (x,(/SIZE(y)/),(/0.0_dp/))
#endif

  END SUBROUTINE reorder21
  !------------------------------------------------------------------------------
  SUBROUTINE reorder12 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:)
    REAL(dp) ,INTENT(in)  :: x (:)

#ifdef __REPLACE_RESHAPE
    CALL util_reshape2(y, x, SIZE(y,1)*SIZE(y,2), SIZE(x))
#else
    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2)/),(/0.0_dp/))
#endif

  END SUBROUTINE reorder12
  !------------------------------------------------------------------------------
  SUBROUTINE reorder32 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:,:)

#ifdef __REPLACE_RESHAPE
    CALL util_reshape2(y, x, SIZE(y,1)*SIZE(y,2), SIZE(x,1)*SIZE(x,2)*SIZE(x,3))
#else
    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2)/),(/0.0_dp/))
#endif

  END SUBROUTINE reorder32
  !------------------------------------------------------------------------------
  SUBROUTINE reorder23 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:)

#ifdef __REPLACE_RESHAPE
    CALL util_reshape2(y, x, SIZE(y,1)*SIZE(y,2)*SIZE(y,3), SIZE(x,1)*SIZE(x,2))
#else
    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2),SIZE(y,3)/),(/0.0_dp/))
#endif
    
  END SUBROUTINE reorder23
  !------------------------------------------------------------------------------
  SUBROUTINE reorder2 (y,x)
    REAL(dp), INTENT(out) :: y (:,:)
    REAL(dp), INTENT(in)  :: x (:,:)

#ifdef __REPLACE_RESHAPE
    CALL util_reshape2(y, x, SIZE(y,1)*SIZE(y,2), SIZE(x,1)*SIZE(x,2))
#else
    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2)/),(/0.0_dp/))
#endif

  END SUBROUTINE reorder2
  !------------------------------------------------------------------------------
  SUBROUTINE reorder3 (y,x)
    REAL(dp), INTENT(out) :: y (:,:,:)
    REAL(dp), INTENT(in)  :: x (:,:,:)
    INTEGER :: k
#ifdef __SX__
    REAL(dp) :: tmp(MAX(SIZE(x,1)*SIZE(x,3),SIZE(y,1)*SIZE(y,3)))
    INTEGER :: j0, j, l, nx1, nx3, ny1, ny3

    nx1 = SIZE(x,1)
    nx3 = SIZE(x,3)
    ny1 = SIZE(y,1)
    ny3 = SIZE(y,3)
    j0 = nx1*nx3+1
    IF (nx1*nx3 /= ny1*ny3) stop 'reorder3: unexpected field sizes'
      
!$-disabled-OMP PARALLEL PRIVATE(j,k,l,tmp)
#else
!$-disabled-OMP PARALLEL PRIVATE(k)
#endif

!$-disabled-OMP DO SCHEDULE(STATIC,1)
    DO k = 1, SIZE(x,2)
#ifdef __SX__
      DO l = 1, nx3
        DO j = 1, nx1
          tmp(j+(l-1)*nx1) = x(j,k,l)
        ENDDO
      ENDDO
      DO j = j0, ny1*ny3
        tmp(j) = 0.0_dp
      ENDDO

      DO l = 1, ny3
        DO j = 1, ny1
          y(j,k,l) = tmp(j+(l-1)*ny1)
        ENDDO
      ENDDO
#elif defined (__REPLACE_RESHAPE)
      CALL util_reshape2(y(:,k,:), x(:,k,:), SIZE(y,1)*SIZE(y,3), SIZE(x,1)*SIZE(x,3))
#else
      y(:,k,:) = RESHAPE (x(:,k,:),(/SIZE(y,1),SIZE(y,3)/),(/0.0_dp/))
#endif
    END DO
!$-disabled-OMP END DO
!$-disabled-OMP END PARALLEL

  END SUBROUTINE reorder3
  !------------------------------------------------------------------------------
  SUBROUTINE reorder4 (y,x)
    REAL(dp), INTENT(out) :: y (:,:,:,:)
    REAL(dp), INTENT(in)  :: x (:,:,:,:)
    INTEGER :: k, l

    DO l = 1, SIZE(x,3)
      DO k = 1, SIZE(x,2)
#ifdef __REPLACE_RESHAPE
        CALL util_reshape2(y(:,k,l,:), x(:,k,l,:), SIZE(y,1)*SIZE(y,4), SIZE(x,1)*SIZE(x,4))
#else
        y(:,k,l,:) = RESHAPE (x(:,k,l,:),(/SIZE(y,1),SIZE(y,4)/),(/0.0_dp/))
#endif
      END DO
    END DO

  END SUBROUTINE reorder4
  !==============================================================================
  ! reorder for a sub range of 2d-array y
  SUBROUTINE subreorder2(y, x, js, je)
    REAL(dp), INTENT(out) :: y(:,:)
    REAL(dp), INTENT(in)  :: x(:,:)
    INTEGER,  INTENT(in)  :: js, je ! needed to calculate how many source-elements we skip
                                    ! we don't skip dest-elements
                                    ! js, je could be replaced by just one argument joff=(je-js+1)

    INTEGER :: ny1, ny2, nx1, nx2, ncopy, nskip, nsource

    nx1 = SIZE(x,1)
    nx2 = SIZE(x,2)
    ny1 = SIZE(y,1)
    ny2 = SIZE(y,2)

    nskip   = (js-1)*ny1      ! number of source-elements we skip
    ncopy   = (je-js+1)*ny1   ! number of elements we want to copy
    nsource = nx1*nx2         ! total size of source

    CALL util_subreshape2(y, x, ncopy, nsource, nskip)

  END SUBROUTINE subreorder2
  !------------------------------------------------------------------------------
  ! reorder for a subrange of 3d-array y
  ! reference implementation, not optimized yet
  SUBROUTINE subreorder3(y,x, iy3s, iy3e)
    REAL(dp), INTENT(out) :: y(:,:,:) 
    REAL(dp), INTENT(in)  :: x(:,:,:)
    INTEGER,  INTENT(in)  :: iy3s, iy3e 

    INTEGER :: ny2, k
    
!-$-omp critical
    ny2 = SIZE(y,2)
    DO k = 1, ny2
      CALL subreorder2(y(:,k,:), x(:,k,:), iy3s, iy3e)
    ENDDO
!-$-omp end critical

  END SUBROUTINE subreorder3

END MODULE mo_tr_alltoall
!------------------------------------------------------------------------------
#ifdef __EXPLICIT
SUBROUTINE explicit_unpack_buf_ls(ls,buf,n1,nfls,nf,nv,i1,i2,j1,j2)
  !
  ! explicit shape size array argument version of unpack_buf_ls
  ! (enforces vectorization)
  !
  USE mo_kind,          ONLY: dp
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n1, nfls, nf, nv
  REAL(dp) :: ls(n1,nfls,nv)
  REAL(dp) :: buf (n1,nf,nv)
    !
  INTEGER :: i1,i2,j1,j2
  INTEGER :: i,j,k
  
!$OMP PARALLEL PRIVATE(i,j,k)
!$OMP DO
  DO k=1,nv
    DO j=i1,i2
!IBM* ASSERT(NODEPS)
      DO i=1,n1
        ls(i,j,k)=buf(i,j-i1+j1,k)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL
END SUBROUTINE explicit_unpack_buf_ls
!------------------------------------------------------------------------------
SUBROUTINE explicit_pack_ls_buf(buf,ls,n1,nfls,nf,nv,i1,i2,j1,j2)
  !
  ! explicit shape size array argument version of pack_ls_buf
  ! (enforces vectorization)
  !
  USE mo_kind,          ONLY: dp
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n1, nfls, nf, nv
  REAL(dp) :: ls(n1,nfls,nv)
  REAL(dp) :: buf (n1,nf,nv)
  
  INTEGER :: i1,i2,j1,j2
  INTEGER :: i,j,k
  
!$OMP PARALLEL PRIVATE(i,j,k)
!$OMP DO
  DO k=1,nv
    DO j=i1,i2
!IBM* ASSERT(NODEPS)
      DO i=1,n1
        buf(i,j-i1+j1,k)=ls(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL
END SUBROUTINE explicit_pack_ls_buf
#endif
!==============================================================================
