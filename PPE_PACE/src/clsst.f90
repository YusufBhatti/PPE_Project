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
SUBROUTINE clsst(krow)

  ! Description:
  !
  ! Passes climate sea-surface-temperatures and sea ice to atmosphere
  !
  ! Method:
  !
  ! This subroutine interpolates the sea-surface temperatures and
  ! sea-ice concentration at each time step and updates *tsw* and *siced*.
  !
  ! *clsst* is called from *gpc*.
  !
  ! Authors: 
  !
  ! U. Schlese, DKRZ, January 1993, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ and M. Esch, MPI, Oct. 1999, modifications for ECHAM5
  ! R. Voss, U. Schlese, MPI, Jan 2000, mods for fract. surface coverage
  ! A. Rhodin, MPI, prescribe sst in SCM run
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! for more details see file AUTHORS
  ! 
  !

  USE mo_kind,          ONLY: dp
  USE mo_memory_g3b,    ONLY: slf, tsw, tsi, seaice, siced, alake
  USE mo_sst,           ONLY: sst, aice, dailysst, dailyice
  USE mo_control,       ONLY: lcouple, ldailysst
  USE mo_physc2,        ONLY: ctfreez
  USE mo_physical_constants,     ONLY: tmelt
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2, ndw1, ndw2, wgtd1, wgtd2
  USE mo_column,        ONLY: sst_1d ! if > 0 : prescribed SST
  USE mo_geoloc,        ONLY: philat_2d
  USE mo_radiation_parameters,     ONLY: nmonth

  IMPLICIT NONE

  INTEGER :: krow

  ! Local variables

  INTEGER :: jl      ! longitude loop index
  INTEGER :: jrow    ! local latitude index
  INTEGER :: nproma  ! number of longitudes on PE
  REAL(dp):: zic(ldc%nproma)
  REAL(dp):: zts(ldc%nproma)
  LOGICAL :: lonorth ! .true. for northern latitude
  LOGICAL :: lsst_1d

  !  Executable statements

  jrow    = krow        ! local latitude index

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

!-- 1. Update temperatures and sea ice

!-- 1.1.1 Uncoupled mode

  IF(.NOT.lcouple) THEN

  lsst_1d = ALLOCATED(sst_1d)

!
!-- 1.1 Annual cycle
!
  IF (nmonth == 0) THEN
    IF ( ldailysst ) THEN
!DIR$ CONCURRENT
      DO jl=1,nproma
        zts(jl)=wgtd1*dailysst(jl,jrow,ndw1)+wgtd2*dailysst(jl,jrow,ndw2)
        zic(jl)=wgtd1*dailyice(jl,jrow,ndw1)+wgtd2*dailyice(jl,jrow,ndw2)
      END DO
    ELSE
!DIR$ CONCURRENT
      DO jl=1,nproma
        zts(jl)=wgt1*sst(jl,jrow,nmw1)+wgt2*sst(jl,jrow,nmw2)
        zic(jl)=wgt1*aice(jl,jrow,nmw1)+wgt2*aice(jl,jrow,nmw2)
      END DO
    ENDIF
!
!-- 1.2 Perpetual month
!
  ELSE
!DIR$ CONCURRENT
    DO jl=1,nproma
      zts(jl)=sst(jl,jrow,nmonth)
      zic(jl)=aice(jl,jrow,nmonth)
    END DO
  END IF
!

  DO jl=1,nproma

    IF(alake(jl,jrow).EQ.0._dp) THEN
      IF(slf(jl,jrow).LT.1._dp) THEN

        zic(jl)=zic(jl)*0.01_dp  !!!  assuming input data is in percent
        seaice(jl,jrow)=MAX(0._dp,MIN(0.99_dp,zic(jl)))
        IF (seaice(jl,jrow).LE.0.01_dp) seaice(jl,jrow)=0.0_dp

        IF (seaice(jl,jrow).GT.0._dp) THEN

          tsw(jl,jrow)=ctfreez

          lonorth = philat_2d(jl,jrow).GT.0 ! true in northern hemisphere
          IF (lonorth) THEN
            siced(jl,jrow)=2._dp
          ELSE
            siced(jl,jrow)=1._dp
          END IF

        ELSE                               ! water

          tsw(jl,jrow)=MAX(zts(jl),ctfreez)
          siced(jl,jrow)=0._dp

        END IF
        !
        ! prescribe sst in SCM
        !
        IF (lsst_1d) THEN
          IF (sst_1d(jl,jrow) > 0._dp) tsw(jl,jrow) = sst_1d(jl,jrow)
        ENDIF

      ELSE                                 ! land
        seaice(jl,jrow)=0._dp
        siced(jl,jrow)=0._dp
        tsw(jl,jrow)=tmelt  !! dummy setting to some reasonable value
        tsi(jl,jrow)=tmelt
      END IF
    END IF
  END DO

  END IF

  RETURN
END SUBROUTINE clsst
