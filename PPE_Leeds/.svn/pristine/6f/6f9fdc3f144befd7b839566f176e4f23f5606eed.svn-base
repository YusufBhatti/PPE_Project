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
SUBROUTINE initemp(krow)
  !
  ! Description:
  !
  ! Initialize surface temperatures 
  !
  ! Method:
  !
  !  Initializes surface temperatures and ice over land and lakes
  !  Sea surface temperatures and sea ice 
  !  are set in *clsst* in case of an uncoupled run 
  !
  ! *initemp* is called from *scan1* at the first timestep only
  !
  ! 
  ! Input data is taken from array *tslclim*.
  ! 
  !
  ! Computation of  soil temperatures:
  ! 
  ! Starting from the tslclim field over land points temperatures are set
  ! in relation to depth of the soil layer and position of the initial
  ! day in the annual cycle. Over sea all levels are set to *tmelt*.
  ! tsl is at 0.07 m
  ! thickness of layers 0.065, 0.254, 0.913, 2.902, 5.700 m
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, April 2000, original source
  ! L. Dumenil, MPI, June 1989, original source of *initemp*
  !             which is now part of this routine
  ! I. Kirchner, MPI, November 2000, date/time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! E. Roeckner, MPI, Sep 2002, initialization of mixed layer ocean
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! V. Gayler, MPI, August 2006, take out soil temperature initialization.
  !             It happens in jsbach (mo_soil: ini_soil_temp).
  ! 
  !
  USE mo_kind,          ONLY: dp
  USE mo_memory_g3b,    ONLY: slf, tsw, tsi, seaice, siced             &
                            , alake, tsl, tslm, tslm1, radtemp
  USE mo_sst,           ONLY: sst, aice, dailysst, dailyice
  USE mo_clim,          ONLY: tslclim
  USE mo_math_constants,ONLY: pi
  USE mo_physical_constants,     ONLY: tmelt
  USE mo_physc2,        ONLY: ctfreez
  USE mo_radiation_parameters,     ONLY: nmonth
  USE mo_control,       ONLY: lcouple, ldailysst, lcolumn, lfractional_mask
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2, ndw1, ndw2, wgtd1, wgtd2
  USE mo_time_control,  ONLY: get_date_components, NDAYLEN,            &
                              get_month_len, start_date, get_year_len
  USE mo_surface_memory, ONLY: land, ice
  USE mo_column,        ONLY: sst_1d_input, ml_input

  IMPLICIT NONE

  INTEGER :: krow

  !  Local scalars: 
  REAL(dp):: zdmax, zkap, zmax, zmin, zsqrt, zday, zdelti, yearl       &
           , dayl, zd
  INTEGER :: jrow, iyday, jl, im, iim, jmax, jmin, jmonth              &
           , nmomid(12), nmonthl(12)
  INTEGER :: yr, mo, dy, hr, mn, se
  INTEGER :: nproma

  !  Local arrays: 
  REAL(dp):: zcount(ldc%nproma), ztsl(ldc%nproma)                      &
           , zic(ldc%nproma), zmth(12), znmea(ldc%nproma)              &
           , zrange(ldc%nproma), zts(ldc%nproma)
          

  !  Intrinsic functions 
  INTRINSIC COS, EXP, SQRT


  !  Executable Statements 

  jrow = krow   ! local continuous latitude index

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

  CALL get_date_components(start_date, yr, mo, dy, hr, mn, se)

  yearl = get_year_len(yr)

  dayl  = REAL(NDAYLEN,dp)

  ! Computational constants

  zkap = 7.5E-7_dp

  ! set length of month
  DO im = 1,12

  ! use variable length of month
     nmonthl(im) = get_month_len(yr,im)
     nmomid(im)  = INT(nmonthl(im)*0.5_dp)

  END DO

  ! Year day for which interpolation is done
  ! zdmax= day of local annual maximum
  ! layer depths

!
  IF(nmonth == 0) THEN
    DO jl=1,nproma
      ztsl(jl)=wgt1*tslclim(jl,jrow,nmw1)+wgt2*tslclim(jl,jrow,nmw2)
      IF (lcolumn .AND. ml_input) ztsl=sst_1d_input(1)
    END DO
  ELSE
    DO jl=1,nproma
      ztsl(jl)=tslclim(jl,jrow,nmonth)
    END DO
  END IF

  IF(lcouple) THEN
    DO jl=1,nproma
      zic(jl)=seaice(jl,jrow)
      zts(jl)=tsw(jl,jrow)
    END DO
  ELSE
    IF ( ldailysst ) THEN
      DO jl=1,nproma
        zts(jl)= wgtd1*dailysst(jl,jrow,ndw1)+wgtd2*dailysst(jl,jrow,ndw2)
        zic(jl)=(wgtd1*dailyice(jl,jrow,ndw1)+wgtd2*dailyice(jl,jrow,ndw2))*0.01_dp
        IF(zic(jl).LE.0.01_dp) zic(jl)=0._dp
      END DO
    ELSE
      DO jl=1,nproma
        IF(nmonth == 0) THEN
          zts(jl)=wgt1*sst(jl,jrow,nmw1)+wgt2*sst(jl,jrow,nmw2)
          zic(jl)=(wgt1*aice(jl,jrow,nmw1)+wgt2*aice(jl,jrow,nmw2))*0.01_dp 
        ELSE
          zts(jl)=sst(jl,jrow,nmonth)
          zic(jl)=aice(jl,jrow,nmonth)*0.01_dp 
        END IF
        IF(zic(jl).LE.0.01_dp) zic(jl)=0._dp
      END DO
    END IF
  END IF
!
!
!   Initialize temperatures and ice
!
!
  DO jl = 1,nproma
    IF ((alake(jl,jrow).GT.0.0_dp .AND. lfractional_mask) .OR.       &! lakes in slf
        (alake(jl,jrow).GE.0.5_dp .AND. .NOT.lfractional_mask)) THEN  ! lakes in slm
      tsi(jl,jrow)=ztsl(jl)
      zdelti=tsi(jl,jrow)-tmelt
      IF (zdelti.LT.0._dp) THEN
        siced(jl,jrow)=1._dp-EXP(-0.005_dp*zdelti**2)
      ELSE
        siced(jl,jrow)=0._dp
      END IF
      IF (siced(jl,jrow).GE.0.1_dp) THEN
        seaice(jl,jrow)=1._dp
        tsw(jl,jrow)=tmelt
      ELSE
        siced(jl,jrow)=0._dp
        seaice(jl,jrow)=0._dp
        tsw(jl,jrow)=MAX(ztsl(jl),tmelt)
        tsi(jl,jrow)=MIN(ztsl(jl),tmelt)
      END IF
    ELSE
      siced(jl,jrow)=0._dp
      seaice(jl,jrow)=0._dp
      tsw(jl,jrow)=MAX(ztsl(jl),tmelt)
      tsi(jl,jrow)=MIN(ztsl(jl),tmelt)
    END IF
    IF (alake(jl,jrow).EQ.0.0_dp .AND. slf(jl,jrow).LT.1.0_dp) THEN         !  ocean
      seaice(jl,jrow)=MAX(0._dp,MIN(0.99_dp,zic(jl)))
      IF (seaice(jl,jrow).GT.0._dp) THEN
        tsi(jl,jrow)=MIN(ztsl(jl),tmelt)
        tsw(jl,jrow)=ctfreez
        zdelti=tsi(jl,jrow)-tmelt
        IF (zdelti.LT.0._dp) THEN
          siced(jl,jrow)=2._dp*(1._dp-EXP(-0.005_dp*zdelti**2))+0.1_dp
        ELSE
          siced(jl,jrow)=0.1_dp
        END IF
      ELSE
        tsi(jl,jrow)=MIN(ztsl(jl),tmelt)
        tsw(jl,jrow)=MAX(ctfreez,zts(jl))
        siced(jl,jrow)=0.1_dp
      END IF
    END IF
  END DO
!
!    Initialize soil temperatures
!

!-- Calculate annual mean temperature

  DO jl = 1, nproma
    znmea(jl) = 0._dp
  END DO

  DO jmonth = 1, 12
    DO jl = 1, nproma
      IF (slf(jl,jrow).GT.0._dp) THEN
        znmea(jl) = znmea(jl) + tslclim(jl,jrow,jmonth) 
      END IF
    END DO
  END DO

  DO jl = 1, nproma
    IF (slf(jl,jrow).GT.0._dp) THEN
      znmea(jl) = znmea(jl)/12._dp
    END IF
  END DO

!--  Month of annual maximum/minimum

  DO jl = 1, nproma
    IF (slf(jl,jrow).GT.0._dp) THEN
      DO jmonth = 1, 12
        zmth(jmonth) = tslclim(jl,jrow,jmonth)
      END DO
      jmax = MAXLOC(zmth,DIM=1)
      jmin = MINLOC(zmth,DIM=1)
      zmax = zmth(jmax)
      zmin = zmth(jmin)
      zrange(jl) = zmax - zmin
      zcount(jl) = REAL(jmax,dp)
    END IF
  END DO

!--  Soil temperature initialization happens in jsbach (mo_soil).
!--  The algorithm below is only used to initialize the upper most soil layer,
!--  i.e. the land surface temperature.

  IF (nmonth == 0) THEN
     ! annual cycle run
     iyday = dy
     DO im = 1, mo-1
        iyday = iyday + nmonthl(im)
     END DO
  ELSE
     ! perpetual months
     iyday = nmomid(nmonth)
     DO im = 1, nmonth-1
        iyday = iyday + nmonthl(im)
     END DO
  END IF
  zday = REAL(iyday,dp)
  zsqrt = SQRT(zkap*yearl*dayl/pi)

  zd = (-0.07_dp) + 0.5*0.065_dp
  DO jl = 1, nproma
     IF (slf(jl,jrow).GT.0._dp) THEN 
        im = INT(zcount(jl)+0.0001_dp)
        zdmax = 0.0_dp
        DO iim = 1,im
           zdmax = zdmax + REAL(nmonthl(iim),dp)
        END DO
        zdmax = zdmax-REAL(nmomid(im),dp)
        tsl(jl,jrow) = znmea(jl)+0.5_dp*zrange(jl)*EXP(-zd/zsqrt) &
                       *COS(2._dp*pi*(zday-zdmax)/yearl-zd/zsqrt)
     ELSE
        tsl(jl,jrow) = tmelt
     END IF
  END DO
!
!  Set all time levels of surface temperature to uppermost soil temp.
!
  DO jl = 1,nproma
    tslm(jl,jrow)  = tsl(jl,jrow)
    tslm1(jl,jrow) = tsl(jl,jrow)
! INIT MO_SURFACE FOR JSBACH-------------------------------------
    land%surface_temperature(jl,jrow) = tslm1(jl,jrow)
    ice%surface_temperature(jl,jrow)  = tsi(jl,jrow)
!    radtemp(jl,jrow)                  = ztsl(jl)
  END DO
! Set surface temperature for radiation calculation according to slf and seaice
  DO jl = 1,nproma
    radtemp(jl,jrow)=slf(jl,jrow)*ztsl(jl)+(1._dp-slf(jl,jrow))* &
                     (seaice(jl,jrow)*tsi(jl,jrow)+(1._dp-seaice(jl,jrow))*tsw(jl,jrow))
  END DO

  RETURN
END SUBROUTINE initemp
