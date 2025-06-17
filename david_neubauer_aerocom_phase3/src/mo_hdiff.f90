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
MODULE mo_hdiff

  USE mo_kind,          ONLY: dp
  USE mo_mpi,           ONLY: p_pe
  USE mo_exception,     ONLY: finish, message
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_memory_sp,     ONLY: sd, stp, svo
  USE mo_control,       ONLY: lmidatm, ltdiag, nn, nk, nkp1,          &
                              nlev, nlevp1
  USE mo_truncation,    ONLY: ntrn
  USE mo_semi_impl,     ONLY: hdamp, vcrit, vmax
  USE mo_physical_constants,     ONLY: earth_radius
  USE mo_diag_tendency, ONLY: pdvor, pdtem, pddiv
  USE mo_time_control,  ONLY: time_step_len, delta_time, lfirst_day,  &
                              get_time_step
  USE mo_filename,      ONLY: find_next_free_unit

  IMPLICIT NONE

  INTEGER,  ALLOCATABLE ::  ncdif(:)
  INTEGER,  ALLOCATABLE ::  iq(:)

  REAL(dp), ALLOCATABLE :: diftcor(:) !  correction profile for temperature

  REAL(dp) :: damhih      !  for extra diffusion in middle atmosphere
  REAL(dp) :: dampth
  REAL(dp) :: difvo       !  coefficient for vorticity.
  REAL(dp) :: difd        !  coefficient for divergence.
  REAL(dp) :: dift        !  coefficient for temperature.

  REAL(dp) :: enstdif     !  factor by which stratospheric
                          !  horizontal diffusion is increased from one
                          !  level to next level above.
  INTEGER :: nlvstd1      !  last (uppermost) layer at which
                          !  stratospheric horizontal diffusion is
                          !  enhanced.
  INTEGER :: nlvstd2      !  first (lowest) layer at which
                          !  stratospheric horizontal diffusion is
                          !  enhanced.

  LOGICAL :: ldiahdf      !  true for statistics of horizontal diffusion
  INTEGER :: ndiahdf = -1 ! I/O unit for hdiff diagnostics

CONTAINS

  SUBROUTINE init_hdiff_diag

    CHARACTER(len=*), PARAMETER :: cbasename = 'hdiffdiag'
    CHARACTER(len=64) :: cfilename

    IF (ndiahdf < 0) THEN
      WRITE(cfilename,'(a,i4.4)') cbasename//'.', p_pe
      ndiahdf = find_next_free_unit(80,89)
      OPEN (unit=ndiahdf,file=TRIM(cfilename))
    ELSE
      CALL finish('init_hdiff_diag','Diagnostics of hdiff already in use.')
    ENDIF

  END SUBROUTINE init_hdiff_diag

  SUBROUTINE finalize_hdiff_diag

    IF (ndiahdf > 0) THEN
      CLOSE(ndiahdf)
    ELSE
      CALL finish('finalize_hdiff_diag','Diagnostics of hdiff not in use.')
    ENDIF

  END SUBROUTINE finalize_hdiff_diag

  SUBROUTINE sudif

    ! Description:
    !
    ! Initializes module mo_hdiff
    !
    ! Authors:
    !
    ! M. Esch, MPI, September 1993, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! M. Esch, MPI, August 1999, modifications for ECHAM5
    ! 
    ! for more details see file AUTHORS
    !
    !  Local arrays: 
    INTEGER :: ihq255(nlev)
    INTEGER :: ihq127(nlev)
    INTEGER :: ihq63(nlev)
    INTEGER :: ihq31(nlev)

    !  Executable statements 

    !-- 1. Set parameter

    IF (lmidatm) THEN
      IF (nn == 255) THEN
        IF (nlev == 95) THEN
          ihq255(1:10)   = 1  !l=10 => ~0.15 hPa  del2
          ihq255(11:25)  = 2  !l=25 => ~1.5 hPa   del4
          ihq255(26:)    = 3  !                   del6
        ELSE
          CALL finish ('sudif', 'Truncation not supported.')
        ENDIF
      ELSE IF (nn == 127) THEN
        IF (nlev == 95) THEN
          ihq127(1:10)   = 1  ! l=10 => ~0.15 hPa  del2
          ihq127(11:25)  = 2  ! l=25 => ~1.50 hPa  del4 
          ihq127(26:)    = 3  !                    del6
        ELSE
          CALL finish ('sudif', 'Truncation not supported.')
        ENDIF
      ELSE IF (nn == 63) THEN
        IF (nlev == 95) THEN
          ihq63(1:10)   = 1  ! l=10 => ~0.15 hPa  del2
          ihq63(11:20)  = 2  ! l=20 => ~0.77 hPa  del4 
          ihq63(21:25)  = 3  ! l=25 => ~1.50 hPa  del6
          ihq63(26:)    = 4  !                    del8
        ELSE IF (nlev == 47) THEN
          ihq63(1:4)   = 1     ! l=4 => ~0.23 hPa  del2 
          ihq63(5:7)   = 2     ! l=7 => ~1.22 hPa  del4 
          ihq63(8:9)   = 3     ! l=9 => ~2.96 hPa  del6 
          ihq63(10:)   = 4 
        ELSE
          CALL finish('sudif','Truncation not supported.')
        ENDIF
      ELSE IF (nn == 31) THEN
        IF (nlev == 47) THEN
          ihq31(1:4)   = 1     ! l=4  => ~0.23 hPa  del2
          ihq31(5:7)   = 2     ! l=7  => ~1.22 hPa  del4 
          ihq31(8:9)   = 3     ! l=9  => ~2.96 hPa  del6 
          ihq31(10:11) = 4     ! l=11 => ~6.37 hPa  del8 
          ihq31(12:)   = 5
        ELSE
          CALL finish('sudif','Truncation not supported.')
        ENDIF
      ELSE
        CALL finish('sudif','Truncation not supported.')
      ENDIF

    ELSE  ! not lmidatm

      IF (nn == 31) THEN
        IF (nlev == 31) THEN
          ihq31(1:3)  = 1
          ihq31(4)    = 2
          ihq31(5)    = 3
          ihq31(6)    = 4
          ihq31(7:)   = 5
        ELSE
          CALL finish('sudif','Truncation not supported.')
        ENDIF
!>>SF
#ifdef HAMMOZ
      ELSE IF (nn == 63) THEN
        IF (nlev == 31) THEN
           ihq63(1:3)  = 1
           ihq63(4)    = 2
           ihq63(5)    = 3
           ihq63(6:)   = 4
        ENDIF
#endif
!<<SF
      ELSE
        CALL finish('sudif','Truncation not supported.')
      ENDIF

    ENDIF

    ncdif(1:nlev) = 0

    !-- 2. Copy to iq

    SELECT CASE (nn)
    CASE (255)
      iq(:) = ihq255(:)      
    CASE (127)
      iq(:) = ihq127(:)      
    CASE (63) 
      iq(:) = ihq63(:)
    CASE (31) 
      iq(:) = ihq31(:)
    CASE DEFAULT
      CALL finish('sudif','This model resolution is not supported.')
    END SELECT

  END SUBROUTINE sudif

  SUBROUTINE hdiff

    ! Description:
    !
    ! Horizontal diffusion.
    !
    ! Method:
    !
    ! This subroutine performs horizontal diffusion.
    !
    ! *hdiff* is called from *stepon*
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, June 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! I. Kirchner, MPI, August 1998, tendency diagnostics
    ! T. Diehl, DKRZ, July 1999, parallel version 
    ! I. Kirchner, MPI, December 2000, time control
    ! 
    ! for more details see file AUTHORS
    !
    !  Local scalars: 
    REAL(dp):: zaa, zcons, zcor, zdamp, zdifd, zdift, zdifvo, znn       &
             , ztwodt, zzcor, zztemp, zzzn
    INTEGER :: in, is, jk, jlev, jn, jndif, jr, ic, i, snsp, nns

    !  Local arrays: 
    REAL(dp):: zcor1(nkp1,nlev), zcor2(nkp1,nlev), znd(nkp1,2*nlev)     &
             , zndifd(nkp1,nlev), zndift(nkp1,nlev), zndifvo(nkp1,nlev) &
             , znt(nkp1,nlev), znvo(nkp1,2*nlev), zzus(nkp1,nlev)       &
             , zhigh(nkp1,nlev)

    INTEGER :: np1(lc%snsp), nindex(lc%nns)
    INTEGER :: mymsp(lc%snsp)
    INTEGER :: istep
    !  Executable statements 

    istep = get_time_step()

    snsp   = lc%snsp
    nns    = lc%nns
    np1    = lc%np1
    nindex = lc%nindex
    mymsp  = lc%mymsp

    ! Each PE loops only over n's owned by itself 

    DO i = 1, nns
      jn = nindex(i)
      zzus(jn,:) = 0._dp
    ENDDO

    !-- 0. Calculate diffusion coefficients 

    difvo = earth_radius*earth_radius/(nk*nkp1*3600._dp*dampth)
    IF (lfirst_day .AND. .NOT. lmidatm) &
         difvo = earth_radius*earth_radius/(nk*nkp1*3600._dp*3._dp)

    difd = 5._dp *difvo
    dift = 0.4_dp*difvo

    !-- 1. Compute diffusion correction factor

    ztwodt = time_step_len
    zaa    = 1._dp/(earth_radius*earth_radius)
    zcons  = 2.5_dp*delta_time/earth_radius

    DO i = 1, nns
      jn = nindex(i)
      IF (jn >= 2) THEN
        in = jn - 1
        zdamp = 1._dp
        IF (lmidatm) THEN
          DO jk = nlev, 1, -1
            IF (in > ntrn(jk)) zdamp = zdamp*hdamp
            IF (jk <= nlvstd2 .AND. jk >= nlvstd1) zdamp=zdamp*enstdif
            zzus(jn,jk) = zdamp
          END DO
        ELSE
          IF (jn <= nk/3) THEN
            DO jk = nlev, 1, -1
              IF (in > ntrn(jk)) zdamp = zdamp*hdamp
              IF (jk <= nlvstd2 .AND. jk >= nlvstd1) zdamp=zdamp*enstdif
              zzus(jn,jk) = zdamp
            END DO
          ELSE
            DO jk = 1, nlev
              zzus(jn,jk) = 1._dp
            END DO
          END IF
        END IF
      ENDIF
    END DO

    IF (lmidatm) THEN
      DO i = 1, nns
        jn = nindex(i)
        DO jk = 1,nlev
          zhigh(jn,jk) = 1._dp
        END DO
      END DO

      IF (.NOT.ldiahdf) THEN
        DO i = 1, nns
          jn = nindex(i)
          IF (jn >= 3) THEN
            znn = REAL(jn-1,dp)
            DO jk = 1,nlev
              zcor = (1._dp+(znn*vmax(jk)-vcrit)) 
              IF (zcor > 1._dp) zhigh(jn,jk) = damhih
            END DO
          ENDIF
        END DO
      ELSE
        ! same loop but WITH statistics
        DO i = 1, nns
          jn = nindex(i)
          IF (jn >= 3) THEN
            znn = REAL(jn-1,dp)
            DO jk = 1,nlev
              zcor = (1._dp+(znn*vmax(jk)-vcrit)) 
              IF (zcor > 1._dp) THEN
                zhigh(jn,jk) = damhih
                WRITE(ndiahdf,'(i12,i5,i5,f18.12)') &
                     istep, jk, jn-1, vmax(jk)   
              ENDIF
            END DO
          ENDIF
        END DO
      END IF
    END IF

    ! Vertical loop

    DO jk = 1, nlev

      jndif = ncdif(jk)
      zztemp = (nkp1*nk*zaa)**(1-iq(jk))*ztwodt
      zdifd = difd*zztemp
      zdifvo = difvo*zztemp
      zdift = dift*zztemp

      DO i = 1, nns
        jn = nindex(i)
        zndifd(jn,jk) = 1._dp
        zndifvo(jn,jk) = 1._dp
        zndift(jn,jk) = 1._dp
      END DO

      DO i = 1, nns
        jn = nindex(i)
        IF (jn >= jndif+1) THEN
          zzzn = REAL(jn - 1 - jndif,dp)
          zztemp = (zaa*zzzn*(zzzn+1._dp))**iq(jk)*zzus(jn,jk)
          IF (lmidatm) THEN
            zndifd(jn,jk)  = 1._dp + zhigh(jn,jk)*zdifd*zztemp
            zndifvo(jn,jk) = 1._dp + zhigh(jn,jk)*zdifvo*zztemp
            zndift(jn,jk)  = 1._dp + zhigh(jn,jk)*zdift*zztemp
          ELSE
            zndifd(jn,jk)  = 1._dp + zdifd*zztemp
            zndifvo(jn,jk) = 1._dp + zdifvo*zztemp
            zndift(jn,jk)  = 1._dp + zdift*zztemp
          ENDIF
        ENDIF
      END DO

      IF ( .NOT. lmidatm) THEN
        IF ( .NOT. ldiahdf) THEN
          DO i = 1, nns
            jn = nindex(i)
            IF (jn >= 3) THEN
              znn = REAL(jn - 1,dp)
              zcor = (1._dp+zcons*(znn*vmax(jk)-vcrit))
              IF (zcor > 1._dp) THEN
                zndifd(jn,jk) = zcor*zndifd(jn,jk)
                zndifvo(jn,jk) = zcor*zndifvo(jn,jk)
              END IF
            ENDIF
          END DO
        ELSE

          ! Same loop but with statistics

          DO i = 1, nns
            jn = nindex(i)
            IF (jn >= 3) THEN
              znn = REAL(jn - 1,dp)
              zcor = (1._dp+zcons*(znn*vmax(jk)-vcrit))
              IF (zcor > 1._dp) THEN
                WRITE (ndiahdf,'(i12,i5,i5,2f18.12)') &
                     istep, jk, jn - 1, zcor, zndifd(jn,jk)
                zndifd(jn,jk) = zcor*zndifd(jn,jk)
                zndifvo(jn,jk) = zcor*zndifvo(jn,jk)
              END IF
            ENDIF
          END DO
        END IF
      END IF

      DO i = 1, nns
        jn = nindex(i)
        IF (jn >= 3) THEN
          znd(jn,jk) = 1._dp/zndifd(jn,jk)
          znd(jn,jk+nlev) = znd(jn,jk)
          znvo(jn,jk) = 1._dp/zndifvo(jn,jk)
          znvo(jn,jk+nlev) = znvo(jn,jk)
        ELSE
          znd(jn,jk) = 1._dp
          znd(jn,jk+nlev) = 1._dp
          znvo(jn,jk) = 1._dp
          znvo(jn,jk+nlev) = 1._dp
        ENDIF
      END DO

      ! Temperature

      DO i = 1, nns
        jn = nindex(i)
        IF (jn >= 2) THEN
          zzcor = 1._dp
          IF (jn-1 > ntrn(jk)) zzcor = 0._dp
          znn = REAL(jn - 1,dp)
          IF (lmidatm) THEN
            zcor = 1._dp
          ELSE
            zcor = (1._dp+zcons*(znn*vmax(jk)-vcrit))
            IF (zcor > 1._dp) zndift(jn,jk) = zcor*zndift(jn,jk)
            IF (zcor <= 1._dp) zcor = 1._dp
          END IF
          zcor2(jn,jk) = zzcor*diftcor(jk)
          zcor1(jn,jk) = zcor2(jn,jk)/zcor
        ELSE
          zcor1(jn,jk) = 0._dp
          zcor2(jn,jk) = 0._dp
        ENDIF
      END DO

      DO i = 1, nns
        jn = nindex(i)
        znt(jn,jk) = 1._dp/zndift(jn,jk)
      END DO
    END DO

    !-- 2. Modify fields

    DO is = 1, snsp

      IF (lmidatm .AND. mymsp(is)==0) CYCLE

      ic = np1(is)

      sd (:,1,is) = sd (:,1,is)*znd (ic,:nlev)
      sd (:,2,is) = sd (:,2,is)*znd (ic, nlev+1:)
      svo(:,1,is) = svo(:,1,is)*znvo(ic,:nlev)
      svo(:,2,is) = svo(:,2,is)*znvo(ic, nlev+1:)

      DO jr = 1, 2

!DIR$ IVDEP
!CDIR NODEP
!OCL NOVREC,NOALIAS

        DO jlev = 1, nlev
          stp(jlev,jr,is) = stp(nlevp1,jr,is)*zcor1(ic,jlev) +         &
               (stp(jlev,jr,is)-zcor2(ic,jlev)*stp(nlevp1,jr,is))*     &
               znt(ic,jlev)
        END DO

      END DO

    END DO

  END SUBROUTINE hdiff
  
END MODULE mo_hdiff
