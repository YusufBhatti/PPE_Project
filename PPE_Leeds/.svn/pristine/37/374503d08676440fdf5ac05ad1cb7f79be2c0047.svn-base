!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_geopack

  ! This module provides the Tsyganenko GEOPACK-2005 routines, modified and
  ! tested for Fortran 90, and a few auxiliary routines (Jan Kazil 2007-09-03 11:52:45).

  ! This collection of subroutines is a result of several upgrades of the original package
  ! written by N. A. Tsyganenko in 1978-1979. This version is dated May 04, 2005. On that
  ! date, the IGRF coefficients were updated according to the recently published table of
  ! IGRF-10 coefficients, so that the main field model now extends through 2010 (a linear
  ! extrapolation is used for 2005 - 2010, based on the table of secular velocities). For
  ! more details, see  

  ! http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html (revision of 03/22/2005)

  ! Prefatory notes to the version of April 22, 2003:

  ! This package represents an in-depth revision of the previous version, with significant
  ! changes in the format of calling statements. Users should familiarize themselves with

  ! the new formats and rules, and accordingly adjust their source codes, as specified
  !  below. Please consult the documentation file

  ! http://modelweb.gsfc.nasa.gov/magnetos/data-based/Geopack-2005.doc for detailed
  ! descriptions of individual subroutines.

  ! The following changes were made to the previous release of GEOPACK (of Jan 5, 2001).

  ! (1) Subroutine IGRF, calculating the Earth's main field:

  !  (a) Two versions of this subroutine are provided here. In the first one (IGRF_GSM)
  !    both input (position) and output (field components) are in the Geocentric Solar-
  !    Magnetospheric Cartesian coordinates, while the second one (IGRF_GEO) uses sphe-
  !    rical geographical (geocentric) coordinates, as in the older releases.

  !  (b) updating of all expansion coefficients is now made separately in the s/r RECALC,
  !    which also takes into account the secular change of the coefficients within
  !    a given year (at the Earth's surface, the rate of the change can reach 7 nT/month).

  !  (c) the optimal length of spherical harmonic expansions is now automatically set
  !    inside the code, based on the radial distance, so that the deviation from the
  !    full-length approximation does not exceed 0.01 nT. (In the previous versions,
  !    the upper limit NM of the order of harmonics had to be specified by users),

  ! (2) Subroutine DIP, calculating the Earth's field in the dipole approximation:

  !  (a) no longer accepts the tilt angle via the list of formal parameters. Instead,
  !    the sine SPS and cosine CPS of that angle are now forwarded into DIP via the
  !    first common block /GEOPACK1/.  Accordingly, there are two options: (i) to
  !    calculate SPS and CPS by calling RECALC before calling DIP, or (ii) to specify
  !    them explicitly. In the last case, SPS and CPS should be specified AFTER the
  !    invocation of RECALC (otherwise they will be overridden by those returned by
  !    RECALC).

  !  (b) the Earth's dipole moment is now calculated by RECALC, based on the table of
  !    the IGRF coefficients and their secular variation rates, for a given year and
  !    the day of the year, and the obtained value of the moment is forwarded into DIP
  !    via the second common block /GEOPACK2/. (In the previous versions, only a single
  !    fixed value was provided for the geodipole moment, corresponding to the most
  !    recent epoch).

  ! (3) Subroutine RECALC now consolidates in one module all calculations needed to
  !    initialize and update the values of coefficients and quantities that vary in
  !    time, either due to secular changes of the main geomagnetic field or as a result
  !    of Earth's diurnal rotation and orbital motion around Sun. That allowed us to
  !    simplify the codes and make them more compiler-independent.

  ! (4) Subroutine GEOMAG is now identical in its structure to other coordinate trans-
  !    formation subroutines. It no longer invokes RECALC from within GEOMAG, but uses
  !    precalculated values of the rotation matrix elements, obtained by a separate
  !    external invocation of RECALC. This eliminates possible interference of the
  !    two subroutines in the old version of the package.

  ! (5) Subroutine TRACE (and the subsidiary modules STEP and RHAND):

  !  (a) no longer needs to specify the highest order of spherical harmonics in the
  !    main geomagnetic field expansion - it is now calculated automatically inside the
  !    IGRF_GSM (or IGRF_GEO) subroutine.

  !  (b) the internal field model can now be explicitly chosen by specifying the para-
  !     meter INNAME (either IGRF_GSM or DIP).

  ! (6) A new subroutine BCARSP was added, providing a conversion of Cartesian field
  !    components into spherical ones (operation, inverse to that performed by the sub-
  !    routine  BSPCAR).

  ! (7) Two new subroutines were added, SHUETAL_MGNP and T96_MGNP, providing the position
  !    of the magnetopause, according to the model of Shue et al. [1998] and the one
  !    used in the T96 magnetospheric magnetic field model.

  USE mo_kind, ONLY: wp
  USE mo_math_constants, ONLY : pi

  IMPLICIT NONE
!  INTEGER, PARAMETER, PRIVATE :: wp = kind(0.0D0)

  INTEGER, PARAMETER :: dp = 2

CONTAINS

  SUBROUTINE geo2mag(lat,lon,mag_lat,mag_lon)

    ! *geo2mag* calculates the magnetic coordinates from given geographic
    ! coordinates. The routine operates for years 1965-2005, for which the IGRF
    ! (International Geomagnetic Reference Field) coefficients are known.

    IMPLICIT NONE

    ! Input arguments:
    REAL (wp) :: lat, lon ! degrees

    ! Output arguments:
    REAL (wp) :: mag_lat, mag_lon ! degrees

    ! Local variables:

    INTEGER :: j

    REAL (wp) :: r, theta, phi
    REAL (wp) :: xgeo, ygeo, zgeo
    REAL (wp) :: xmag, ymag, zmag

    ! Transform the latitude and longitude (degrees) into colatitude theta
    ! and longitude phi (radians):
    theta = pi*(0.5_wp-lat/180.0_wp)
    phi = pi*lon/180.0_wp

    ! We are computing everything on the Earth's surface:
    r = 1.0_wp ! Earth radii

    ! Transform spherical geographic coordinates into cartesian coordinates:

    j = 1

    CALL sphcar(r,theta,phi,xgeo,ygeo,zgeo,j)

    ! Transform the geographic cartesian coordinates into magnetic cartesian
    ! coordinates:

    j = 1

    CALL geomag(xgeo,ygeo,zgeo,xmag,ymag,zmag,j)

    ! Transform cartesian magnetic coordinates into spherical coordinates:

    j = -1

    CALL sphcar(r,theta,phi,xmag,ymag,zmag,j)

    ! Transform colatitude theta and longitude phi (radians) into latitude
    ! and longitude (degrees) :
    mag_lat = 180.0_wp*(0.5_wp-real(theta,kind=wp)/pi)
    mag_lon = 180.0_wp*phi/pi

  END SUBROUTINE geo2mag


  !=============================================================================


  FUNCTION vertical_cutoff_rigidity(mag_lat)

    ! *vertical_cutoff_rigidity* calculates the vertical cutoff rigidity (GV)
    ! from a given magnetic latitude (degrees).

    IMPLICIT NONE
    REAL (wp) :: vertical_cutoff_rigidity

    ! Input arguments:
    REAL (wp) :: mag_lat ! degrees

    vertical_cutoff_rigidity = 14.9_wp*(cos(pi*mag_lat/180.0_wp))**4.0_wp ! GV

  END FUNCTION vertical_cutoff_rigidity


  !=============================================================================


  SUBROUTINE igrf_gsm(xgsm,ygsm,zgsm,hxgsm,hygsm,hzgsm)

    ! CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE GEOCENTRIC SOLAR
    ! MAGNETOSPHERIC COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL
    ! COEFFICIENTS (e.g., http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html Revised: 22 March, 2005)


    ! BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR IF THE DATE/TIME (IYEAR,IDAY,IHOUR,MIN,ISEC)
    ! WAS CHANGED, THE MODEL COEFFICIENTS AND GEO-GSM ROTATION MATRIX ELEMENTS SHOULD BE UPDATED
    ! BY CALLING THE SUBROUTINE RECALC

    !-----INPUT PARAMETERS:

    !    XGSM,YGSM,ZGSM - CARTESIAN GSM COORDINATES (IN UNITS RE=6371.2 KM)

    !-----OUTPUT PARAMETERS:

    !    HXGSM,HYGSM,HZGSM - CARTESIAN GSM COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA

    !    LAST MODIFICATION:  MAY 4, 2005.
    !    THIS VERSION OF THE CODE ACCEPTS DATES FROM 1965 THROUGH 2010.

    !    AUTHOR: N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: xgsm, ygsm, zgsm, hxgsm, hygsm, hzgsm

    REAL (wp) :: xgeo, ygeo, zgeo

    REAL (wp) :: a(14), b(14), g(105), h(105), rec(105)

    COMMON /geopack2/g, h, rec

    REAL (wp) :: rho2, r, c, rho, s, cf, sf, pp, p, d, bbr, bbt, bbf, w, &
         x, y, q, z, bi, p2, d2, an, e, hh, qq, xk, dk, pm, br, bt, bf, he, &
         hxgeo, hygeo, hzgeo

    INTEGER :: irp3, nm, k, n, m, mm, mn

    CALL geogsm(xgeo,ygeo,zgeo,xgsm,ygsm,zgsm,-1)
    rho2 = xgeo**2.0_wp + ygeo**2.0_wp
    r = sqrt(rho2+zgeo**2.0_wp)
    c = zgeo/r
    rho = sqrt(rho2)
    s = rho/r
    IF (s<1.E-5_wp) THEN
      cf = 1._wp
      sf = 0._wp
    ELSE
      cf = xgeo/rho
      sf = ygeo/rho
    END IF

    pp = 1._wp/r
    p = pp

    ! IN THIS NEW VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
    !   HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
    !     ON THE VALUE OF THE RADIAL DISTANCE R:

    irp3 = r + 2
    nm = 3 + 30/irp3
    IF (nm>13) nm = 13

    k = nm + 1
    DO n = 1, k
      p = p*pp
      a(n) = p
      b(n) = p*n
    END DO

    p = 1._wp
    d = 0._wp
    bbr = 0._wp
    bbt = 0._wp
    bbf = 0._wp

    DO m = 1, k
      IF (m==1) THEN
        x = 0._wp
        y = 1._wp
      ELSE
        mm = m - 1
        w = x
        x = w*cf + y*sf
        y = y*cf - w*sf
      END IF
      q = p
      z = d
      bi = 0._wp
      p2 = 0._wp
      d2 = 0._wp
      DO n = m, k
        an = a(n)
        mn = n*(n-1)/2 + m
        e = g(mn)
        hh = h(mn)
        w = e*y + hh*x
        bbr = bbr + b(n)*w*q
        bbt = bbt - an*w*z
        IF (m/=1) THEN
          qq = q
          IF (s<1.E-5_wp) qq = z
          bi = bi + an*(e*x-hh*y)*qq
        END IF
        xk = rec(mn)
        dk = c*z - s*q - xk*d2
        pm = c*q - xk*p2
        d2 = z
        p2 = q
        z = dk
        q = pm
      END DO
      d = s*d + c*p
      p = s*p
      IF (m==1) CYCLE
      bi = bi*mm
      bbf = bbf + bi
    END DO

    br = bbr
    bt = bbt
    IF (s<1.E-5_wp) THEN
      IF (c<0._wp) bbf = -bbf
      bf = bbf
    ELSE
      bf = bbf/s
    END IF

    he = br*s + bt*c
    hxgeo = he*cf - bf*sf
    hygeo = he*sf + bf*cf
    hzgeo = br*c - bt*s

    CALL geogsm(hxgeo,hygeo,hzgeo,hxgsm,hygsm,hzgsm,1)


  END SUBROUTINE igrf_gsm


  !=============================================================================


  SUBROUTINE igrf_geo(r,theta,phi,br,btheta,bphi)

    ! CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE SPHERICAL GEOGRAPHIC
    ! (GEOCENTRIC) COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL
    ! COEFFICIENTS  (e.g., http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html, revised: 22 March, 2005)

    ! BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR IF THE DATE (IYEAR AND IDAY) WAS CHANGED,
    ! THE MODEL COEFFICIENTS SHOULD BE UPDATED BY CALLING THE SUBROUTINE RECALC

    !-----INPUT PARAMETERS:

    !  R, THETA, PHI - SPHERICAL GEOGRAPHIC (GEOCENTRIC) COORDINATES:
    !  RADIAL DISTANCE R IN UNITS RE=6371.2 KM, COLATITUDE THETA AND LONGITUDE PHI IN RADIANS

    !-----OUTPUT PARAMETERS:

    !    BR, BTHETA, BPHI - SPHERICAL COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
    !     (POSITIVE BR OUTWARD, BTHETA SOUTHWARD, BPHI EASTWARD)

    !    LAST MODIFICATION:  MAY 4, 2005.
    !    THIS VERSION OF THE  CODE ACCEPTS DATES FROM 1965 THROUGH 2010.

    !    AUTHOR: N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: r, theta, phi, br, btheta, bphi

    REAL (wp) :: a(14), b(14), g(105), h(105), rec(105)

    COMMON /geopack2/g, h, rec

    REAL (wp) :: c, s, cf, sf, pp, p, d, bbr, bbt, bbf, w, x, y, q, z, bi, &
         p2, d2, an, e, hh, qq, xk, dk, pm

    INTEGER :: irp3, nm, k, n, m, mm, mn

    c = cos(theta)
    s = sin(theta)
    cf = cos(phi)
    sf = sin(phi)

    pp = 1._wp/r
    p = pp

    ! IN THIS NEW VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
    !   HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
    !     ON THE VALUE OF THE RADIAL DISTANCE R:

    irp3 = r + 2
    nm = 3 + 30/irp3
    IF (nm>13) nm = 13

    k = nm + 1
    DO n = 1, k
      p = p*pp
      a(n) = p
      b(n) = p*n
    END DO

    p = 1._wp
    d = 0._wp
    bbr = 0._wp
    bbt = 0._wp
    bbf = 0._wp

    DO m = 1, k
      IF (m==1) THEN
        x = 0._wp
        y = 1._wp
      ELSE
        mm = m - 1
        w = x
        x = w*cf + y*sf
        y = y*cf - w*sf
      END IF
      q = p
      z = d
      bi = 0._wp
      p2 = 0._wp
      d2 = 0._wp
      DO n = m, k
        an = a(n)
        mn = n*(n-1)/2 + m
        e = g(mn)
        hh = h(mn)
        w = e*y + hh*x
        bbr = bbr + b(n)*w*q
        bbt = bbt - an*w*z
        IF (m/=1) THEN
          qq = q
          IF (s<1.E-5_wp) qq = z
          bi = bi + an*(e*x-hh*y)*qq
        END IF
        xk = rec(mn)
        dk = c*z - s*q - xk*d2
        pm = c*q - xk*p2
        d2 = z
        p2 = q
        z = dk
        q = pm
      END DO
      d = s*d + c*p
      p = s*p
      IF (m==1) CYCLE
      bi = bi*mm
      bbf = bbf + bi
    END DO

    br = bbr
    btheta = bbt
    IF (s<1.E-5_wp) THEN
      IF (c<0._wp) bbf = -bbf
      bphi = bbf

    ELSE
      bphi = bbf/s
    END IF

  END SUBROUTINE igrf_geo


  !=============================================================================


  SUBROUTINE dip(xgsm,ygsm,zgsm,bxgsm,bygsm,bzgsm)

    ! CALCULATES GSM COMPONENTS OF A GEODIPOLE FIELD WITH THE DIPOLE MOMENT
    ! CORRESPONDING TO THE EPOCH, SPECIFIED BY CALLING SUBROUTINE RECALC (SHOULD BE
    ! INVOKED BEFORE THE FIRST USE OF THIS ONE AND IN CASE THE DATE/TIME WAS CHANGED).

    !--INPUT PARAMETERS: XGSM,YGSM,ZGSM - GSM COORDINATES IN RE (1 RE = 6371.2 km)

    !--OUTPUT PARAMETERS: BXGSM,BYGSM,BZGSM - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.

    ! LAST MODIFICATION: MAY 4, 2005

    ! AUTHOR: N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: xgsm, ygsm, zgsm, bxgsm, bygsm, bzgsm
    REAL (wp) :: aaa(10), sps, cps, bbb(23)
    REAL (wp) :: g(105), h(105), rec(105)

    COMMON /geopack1/aaa, sps, cps, bbb
    COMMON /geopack2/g, h, rec

    REAL (wp) :: dipmom, p, u, v, t, q

    dipmom = sqrt(g(2)**2.0_wp+g(3)**2.0_wp+h(3)**2.0_wp)

    p = xgsm**2.0_wp
    u = zgsm**2.0_wp
    v = 3._wp*zgsm*xgsm
    t = ygsm**2.0_wp
    q = dipmom/sqrt(p+t+u)**5.0_wp
    bxgsm = q*((t+u-2._wp*p)*sps-v*cps)
    bygsm = -3._wp*ygsm*q*(xgsm*sps+zgsm*cps)
    bzgsm = q*((p+t-2._wp*u)*cps-v*sps)


  END SUBROUTINE dip


  !=============================================================================


  SUBROUTINE sun(iyear,iday,ihour,min,isec,gst,slong,srasn,sdec)

    ! CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
    ! WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)

    !-------  INPUT PARAMETERS:
    ! IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
    !   AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).

    !-------  OUTPUT PARAMETERS:
    ! GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
    ! SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
    ! ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
    ! RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.

    ! LAST MODIFICATION:  MARCH 31, 2003 (ONLY SOME NOTATION CHANGES)

    !    ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead

    IMPLICIT NONE

    INTEGER :: iyear, iday, ihour, min, isec
    REAL (wp) :: gst, slong, srasn, sdec

    REAL (wp) :: dj, fday, rad

    REAL (wp) :: t, vl, g, obliq, sob, slp, sind, cosd, sc

    DATA rad/57.295779513_wp/

    IF (iyear>=1901 .AND. iyear<=2099) THEN
      fday = real(ihour*3600+min*60+isec,kind=wp)/86400._wp
      dj = 365*(iyear-1900) + (iyear-1901)/4 + iday - 0.5_wp + fday
      t = dj/36525._wp
      vl = mod(279.696678_wp+0.9856473354_wp*dj,360._wp)
      gst = mod(279.690983_wp+.9856473354_wp*dj+360._wp*fday+180._wp, &
           360._wp)/rad
      g = mod(358.475845_wp+0.985600267_wp*dj,360._wp)/rad
      slong = (vl+(1.91946_wp-0.004789_wp*t)*sin(g)+ &
           0.020094_wp*sin(2._wp*g))/rad
      IF (slong>6.2831853_wp) slong = slong - 6.2831853_wp
      IF (slong<0._wp) slong = slong + 6.2831853_wp
      obliq = (23.45229_wp-0.0130125_wp*t)/rad
      sob = sin(obliq)
      slp = slong - 9.924E-5_wp

      !  THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
      !  THE ORBITAL MOTION OF THE EARTH

      sind = sob*sin(slp)
      cosd = sqrt(1._wp-sind**2.0_wp)
      sc = sind/cosd
      sdec = atan(sc)
      srasn = pi - atan2(cos(obliq)/sob*sc,-cos(slp)/cosd)

    END IF

  END SUBROUTINE sun


  !=============================================================================


  SUBROUTINE sphcar(r,theta,phi,x,y,z,j)

    !  CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
    !   (THETA AND PHI IN RADIANS).

    !                 J>0            J<0
    !-----INPUT:   J,R,THETA,PHI     J,X,Y,Z
    !----OUTPUT:      X,Y,Z        R,THETA,PHI

    ! NOTE: AT THE POLES (X=0 AND Y=0) WE ASSUME PHI=0 (WHEN CONVERTING
    !       FROM CARTESIAN TO SPHERICAL COORDS, I.E., FOR J<0)

    !  LAST MOFIFICATION:  APRIL 1, 2003 (ONLY SOME NOTATION CHANGES AND MORE
    !                        COMMENTS ADDED)

    !  AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: r, theta, phi, x, y, z

    REAL (wp) :: sq

    INTEGER :: j

    IF (j>0) THEN
      sq = r*sin(theta)
      x = sq*cos(phi)
      y = sq*sin(phi)
      z = r*cos(theta)

    ELSE
      sq = x**2.0_wp + y**2.0_wp
      r = sqrt(sq+z**2.0_wp)
      IF (sq/=0._wp) THEN
        sq = sqrt(sq)
        phi = atan2(y,x)
        theta = atan2(sq,z)
        IF (phi<0._wp) phi = phi + 2 * pi
      ELSE
        phi = 0._wp
        IF (z<0._wp) THEN
          theta = pi
        ELSE
          theta = 0._wp
        END IF
      END IF
    END IF

  END SUBROUTINE sphcar


  !=============================================================================


  SUBROUTINE bspcar(theta,phi,br,btheta,bphi,bx,by,bz)

    !  CALCULATES CARTESIAN FIELD COMPONENTS FROM SPHERICAL ONES
    !-----INPUT:   THETA,PHI - SPHERICAL ANGLES OF THE POINT IN RADIANS
    !             BR,BTHETA,BPHI -  SPHERICAL COMPONENTS OF THE FIELD
    !-----OUTPUT:  BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD

    !  LAST MOFIFICATION:  APRIL 1, 2003 (ONLY SOME NOTATION CHANGES)

    !  WRITTEN BY:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: theta, phi, br, btheta, bphi, bx, by, bz

    REAL (wp) :: s, c, sf, cf, be

    s = sin(theta)
    c = cos(theta)
    sf = sin(phi)
    cf = cos(phi)
    be = br*s + btheta*c
    bx = be*cf - bphi*sf
    by = be*sf + bphi*cf
    bz = br*c - btheta*s


  END SUBROUTINE bspcar


  !=============================================================================


  SUBROUTINE bcarsp(x,y,z,bx,by,bz,br,btheta,bphi)

    ! CALCULATES SPHERICAL FIELD COMPONENTS FROM THOSE IN CARTESIAN SYSTEM

    !-----INPUT:   X,Y,Z  - CARTESIAN COMPONENTS OF THE POSITION VECTOR
    !             BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD VECTOR
    !-----OUTPUT:  BR,BTHETA,BPHI - SPHERICAL COMPONENTS OF THE FIELD VECTOR

    ! NOTE: AT THE POLES (THETA=0 OR THETA=PI) WE ASSUME PHI=0,
    !       AND HENCE BTHETA=BX, BPHI=BY

    !  WRITTEN AND ADDED TO THIS PACKAGE:  APRIL 1, 2003,
    !  AUTHOR:   N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: x, y, z, bx, by, bz, br, btheta, bphi
    REAL (wp) :: rho2, r, rho, cphi, sphi, ct, st

    rho2 = x**2.0_wp + y**2.0_wp
    r = sqrt(rho2+z**2.0_wp)
    rho = sqrt(rho2)

    IF (rho/=0._wp) THEN
      cphi = x/rho
      sphi = y/rho
    ELSE
      cphi = 1._wp
      sphi = 0._wp
    END IF

    ct = z/r
    st = rho/r

    br = (x*bx+y*by+z*bz)/r
    btheta = (bx*cphi+by*sphi)*ct - bz*st
    bphi = by*cphi - bx*sphi


  END SUBROUTINE bcarsp


  !=============================================================================


  SUBROUTINE recalc(iyear,iday,ihour,min,isec)

    ! 1. PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN
    !    SEVERAL COORDINATE SYSTEMS, MOST FREQUENTLY USED IN SPACE PHYSICS.

    ! 2. PREPARES COEFFICIENTS USED IN THE CALCULATION OF THE MAIN GEOMAGNETIC FIELD
    !     (IGRF MODEL)

    ! THIS SUBROUTINE SHOULD BE INVOKED BEFORE USING THE FOLLOWING SUBROUTINES:
    !   IGRF_GEO, IGRF_GSM, DIP, GEOMAG, GEOGSM, MAGSM, SMGSM, GSMGSE, GEIGEO.

    ! THERE IS NO NEED TO REPEATEDLY INVOKE RECALC, IF MULTIPLE CALCULATIONS ARE MADE
    !   FOR THE SAME DATE AND TIME.

    !-----INPUT PARAMETERS:

    !    IYEAR   -  YEAR NUMBER (FOUR DIGITS)
    !    IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
    !    IHOUR -  HOUR OF DAY (00 TO 23)
    !    MIN   -  MINUTE OF HOUR (00 TO 59)
    !    ISEC  -  SECONDS OF MINUTE (00 TO 59)

    !-----OUTPUT PARAMETERS:   NONE (ALL OUTPUT QUANTITIES ARE PLACED
    !                        INTO THE COMMON BLOCKS /GEOPACK1/ AND /GEOPACK2/)

    !   OTHER SUBROUTINES CALLED BY THIS ONE: SUN

    !   AUTHOR:  N.A. TSYGANENKO
    !   DATE:    DEC.1, 1991

    !   CORRECTION OF MAY 9, 2006:  INTERPOLATION OF THE COEFFICIENTS (BETWEEN
    !    LABELS 50 AND 105) IS NOW MADE THROUGH THE LAST ELEMENT OF THE ARRAYS
    !    G(105)  AND H(105) (PREVIOUSLY MADE ONLY THROUGH N=66, WHICH IN SOME
    !    CASES CAUSED RUNTIME ERRORS)

    !   REVISION OF MAY 3, 2005:
    !    The table of IGRF coefficients was extended to include those for the epoch 2005
    !      the maximal order of spherical harmonics was also increased up to 13
    !        (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)

    !   REVISION OF APRIL 3, 2003:
    !    The code now includes preparation of the model coefficients for the subroutines
    !      IGRF and GEOMAG. This eliminates the need for the SAVE statements, used in the
    !       old versions, making the codes easier and more compiler-independent.

    IMPLICIT NONE

    INTEGER :: iyear, iday, ihour, min, isec

    REAL (wp) :: g65(105), h65(105), g70(105), h70(105), g75(105), &
         h75(105), g80(105), h80(105), g85(105), h85(105), g90(105), &
         h90(105), g95(105), h95(105), g00(105), h00(105), g05(105), &
         h05(105), dg05(45), dh05(45)

    REAL (wp) :: st0, ct0, sl0, cl0, ctcl, stcl, ctsl, stsl, sfi, cfi, &
         sps, cps, shi, chi, hi, psi, xmut, a11, a21, a31, a12, a22, a32, &
         a13, a23, a33, ds3, cgst, sgst, ba(6)

    REAL (wp) :: g(105), h(105), rec(105)

    COMMON /geopack1/st0, ct0, sl0, cl0, ctcl, stcl, ctsl, stsl, sfi, cfi, &
         sps, cps, shi, chi, hi, psi, xmut, a11, a21, a31, a12, a22, a32, &
         a13, a23, a33, ds3, cgst, sgst, ba

    ! THE COMMON BLOCK /GEOPACK1/ CONTAINS ELEMENTS OF THE ROTATION MATRICES AND OTHER
    !  PARAMETERS RELATED TO THE COORDINATE TRANSFORMATIONS PERFORMED BY THIS PACKAGE

    COMMON /geopack2/g, h, rec

    ! THE COMMON BLOCK /GEOPACK2/ CONTAINS COEFFICIENTS OF THE IGRF FIELD MODEL, CALCULATED
    !   FOR A GIVEN YEAR AND DAY FROM THEIR STANDARD EPOCH VALUES. THE ARRAY REC CONTAINS
    !   COEFFICIENTS USED IN THE RECURSION RELATIONS FOR LEGENDRE ASSOCIATE POLYNOMIALS.

    INTEGER :: iy, n, n2, m, mn, mnn

    REAL (wp) :: dt, f2, f1, s, p, aa, g10, g11, h11, sq, sqq, sqr, gst, &
         slong, srasn, sdec, s1, s2, s3, dip1, dip2, dip3, y1, y2, y3, y, z1, &
         z2, z3, dj, t, obliq, dz1, dz2, dz3, dy1, dy2, dy3, exmagx, exmagy, &
         exmagz, eymagx, eymagy

    DATA g65/0._wp, -30334._wp, -2119._wp, -1662._wp, 2997._wp, 1594._wp, &
         1297._wp, -2038._wp, 1292._wp, 856._wp, 957._wp, 804._wp, 479._wp, &
         -390._wp, 252._wp, -219._wp, 358._wp, 254._wp, -31._wp, -157._wp, &
         -62._wp, 45._wp, 61._wp, 8._wp, -228._wp, 4._wp, 1._wp, -111._wp, &
         75._wp, -57._wp, 4._wp, 13._wp, -26._wp, -6._wp, 13._wp, 1._wp, &
         13._wp, 5._wp, -4._wp, -14._wp, 0._wp, 8._wp, -1._wp, 11._wp, 4._wp, &
         8._wp, 10._wp, 2._wp, -13._wp, 10._wp, -1._wp, -1._wp, 5._wp, 1._wp, &
         -2._wp, -2._wp, -3._wp, 2._wp, -5._wp, -2._wp, 4._wp, 4._wp, 0._wp, &
         2._wp, 2._wp, 0._wp, 39*0.0_wp/
    DATA h65/0._wp, 0._wp, 5776._wp, 0._wp, -2016._wp, 114._wp, 0._wp, &
         -404._wp, 240._wp, -165._wp, 0._wp, 148._wp, -269._wp, 13._wp, &
         -269._wp, 0._wp, 19._wp, 128._wp, -126._wp, -97._wp, 81._wp, 0._wp, &
         -11._wp, 100._wp, 68._wp, -32._wp, -8._wp, -7._wp, 0._wp, -61._wp, &
         -27._wp, -2._wp, 6._wp, 26._wp, -23._wp, -12._wp, 0._wp, 7._wp, &
         -12._wp, 9._wp, -16._wp, 4._wp, 24._wp, -3._wp, -17._wp, 0._wp, &
         -22._wp, 15._wp, 7._wp, -4._wp, -5._wp, 10._wp, 10._wp, -4._wp, &
         1._wp, 0._wp, 2._wp, 1._wp, 2._wp, 6._wp, -4._wp, 0._wp, -2._wp, &
         3._wp, 0._wp, -6._wp, 39*0.0_wp/

    DATA g70/0._wp, -30220._wp, -2068._wp, -1781._wp, 3000._wp, 1611._wp, &
         1287._wp, -2091._wp, 1278._wp, 838._wp, 952._wp, 800._wp, 461._wp, &
         -395._wp, 234._wp, -216._wp, 359._wp, 262._wp, -42._wp, -160._wp, &
         -56._wp, 43._wp, 64._wp, 15._wp, -212._wp, 2._wp, 3._wp, -112._wp, &
         72._wp, -57._wp, 1._wp, 14._wp, -22._wp, -2._wp, 13._wp, -2._wp, &
         14._wp, 6._wp, -2._wp, -13._wp, -3._wp, 5._wp, 0._wp, 11._wp, 3._wp, &
         8._wp, 10._wp, 2._wp, -12._wp, 10._wp, -1._wp, 0._wp, 3._wp, 1._wp, &
         -1._wp, -3._wp, -3._wp, 2._wp, -5._wp, -1._wp, 6._wp, 4._wp, 1._wp, &
         0._wp, 3._wp, -1._wp, 39*0.0_wp/
    DATA h70/0._wp, 0._wp, 5737._wp, 0._wp, -2047._wp, 25._wp, 0._wp, &
         -366._wp, 251._wp, -196._wp, 0._wp, 167._wp, -266._wp, 26._wp, &
         -279._wp, 0._wp, 26._wp, 139._wp, -139._wp, -91._wp, 83._wp, 0._wp, &
         -12._wp, 100._wp, 72._wp, -37._wp, -6._wp, 1._wp, 0._wp, -70._wp, &
         -27._wp, -4._wp, 8._wp, 23._wp, -23._wp, -11._wp, 0._wp, 7._wp, &
         -15._wp, 6._wp, -17._wp, 6._wp, 21._wp, -6._wp, -16._wp, 0._wp, &
         -21._wp, 16._wp, 6._wp, -4._wp, -5._wp, 10._wp, 11._wp, -2._wp, &
         1._wp, 0._wp, 1._wp, 1._wp, 3._wp, 4._wp, -4._wp, 0._wp, -1._wp, &
         3._wp, 1._wp, -4._wp, 39*0.0_wp/

    DATA g75/0._wp, -30100._wp, -2013._wp, -1902._wp, 3010._wp, 1632._wp, &
         1276._wp, -2144._wp, 1260._wp, 830._wp, 946._wp, 791._wp, 438._wp, &
         -405._wp, 216._wp, -218._wp, 356._wp, 264._wp, -59._wp, -159._wp, &
         -49._wp, 45._wp, 66._wp, 28._wp, -198._wp, 1._wp, 6._wp, -111._wp, &
         71._wp, -56._wp, 1._wp, 16._wp, -14._wp, 0._wp, 12._wp, -5._wp, &
         14._wp, 6._wp, -1._wp, -12._wp, -8._wp, 4._wp, 0._wp, 10._wp, 1._wp, &
         7._wp, 10._wp, 2._wp, -12._wp, 10._wp, -1._wp, -1._wp, 4._wp, 1._wp, &
         -2._wp, -3._wp, -3._wp, 2._wp, -5._wp, -2._wp, 5._wp, 4._wp, 1._wp, &
         0._wp, 3._wp, -1._wp, 39*0.0_wp/
    DATA h75/0._wp, 0._wp, 5675._wp, 0._wp, -2067._wp, -68._wp, 0._wp, &
         -333._wp, 262._wp, -223._wp, 0._wp, 191._wp, -265._wp, 39._wp, &
         -288._wp, 0._wp, 31._wp, 148._wp, -152._wp, -83._wp, 88._wp, 0._wp, &
         -13._wp, 99._wp, 75._wp, -41._wp, -4._wp, 11._wp, 0._wp, -77._wp, &
         -26._wp, -5._wp, 10._wp, 22._wp, -23._wp, -12._wp, 0._wp, 6._wp, &
         -16._wp, 4._wp, -19._wp, 6._wp, 18._wp, -10._wp, -17._wp, 0._wp, &
         -21._wp, 16._wp, 7._wp, -4._wp, -5._wp, 10._wp, 11._wp, -3._wp, &
         1._wp, 0._wp, 1._wp, 1._wp, 3._wp, 4._wp, -4._wp, -1._wp, -1._wp, &
         3._wp, 1._wp, -5._wp, 39*0.0_wp/

    DATA g80/0._wp, -29992._wp, -1956._wp, -1997._wp, 3027._wp, 1663._wp, &
         1281._wp, -2180._wp, 1251._wp, 833._wp, 938._wp, 782._wp, 398._wp, &
         -419._wp, 199._wp, -218._wp, 357._wp, 261._wp, -74._wp, -162._wp, &
         -48._wp, 48._wp, 66._wp, 42._wp, -192._wp, 4._wp, 14._wp, -108._wp, &
         72._wp, -59._wp, 2._wp, 21._wp, -12._wp, 1._wp, 11._wp, -2._wp, &
         18._wp, 6._wp, 0._wp, -11._wp, -7._wp, 4._wp, 3._wp, 6._wp, -1._wp, &
         5._wp, 10._wp, 1._wp, -12._wp, 9._wp, -3._wp, -1._wp, 7._wp, 2._wp, &
         -5._wp, -4._wp, -4._wp, 2._wp, -5._wp, -2._wp, 5._wp, 3._wp, 1._wp, &
         2._wp, 3._wp, 0._wp, 39*0.0_wp/
    DATA h80/0._wp, 0._wp, 5604._wp, 0._wp, -2129._wp, -200._wp, 0._wp, &
         -336._wp, 271._wp, -252._wp, 0._wp, 212._wp, -257._wp, 53._wp, &
         -297._wp, 0._wp, 46._wp, 150._wp, -151._wp, -78._wp, 92._wp, 0._wp, &
         -15._wp, 93._wp, 71._wp, -43._wp, -2._wp, 17._wp, 0._wp, -82._wp, &
         -27._wp, -5._wp, 16._wp, 18._wp, -23._wp, -10._wp, 0._wp, 7._wp, &
         -18._wp, 4._wp, -22._wp, 9._wp, 16._wp, -13._wp, -15._wp, 0._wp, &
         -21._wp, 16._wp, 9._wp, -5._wp, -6._wp, 9._wp, 10._wp, -6._wp, &
         2._wp, 0._wp, 1._wp, 0._wp, 3._wp, 6._wp, -4._wp, 0._wp, -1._wp, &
         4._wp, 0._wp, -6._wp, 39*0.0_wp/

    DATA g85/0._wp, -29873._wp, -1905._wp, -2072._wp, 3044._wp, 1687._wp, &
         1296._wp, -2208._wp, 1247._wp, 829._wp, 936._wp, 780._wp, 361._wp, &
         -424._wp, 170._wp, -214._wp, 355._wp, 253._wp, -93._wp, -164._wp, &
         -46._wp, 53._wp, 65._wp, 51._wp, -185._wp, 4._wp, 16._wp, -102._wp, &
         74._wp, -62._wp, 3._wp, 24._wp, -6._wp, 4._wp, 10._wp, 0._wp, &
         21._wp, 6._wp, 0._wp, -11._wp, -9._wp, 4._wp, 4._wp, 4._wp, -4._wp, &
         5._wp, 10._wp, 1._wp, -12._wp, 9._wp, -3._wp, -1._wp, 7._wp, 1._wp, &
         -5._wp, -4._wp, -4._wp, 3._wp, -5._wp, -2._wp, 5._wp, 3._wp, 1._wp, &
         2._wp, 3._wp, 0._wp, 39*0.0_wp/
    DATA h85/0._wp, 0._wp, 5500._wp, 0._wp, -2197._wp, -306._wp, 0._wp, &
         -310._wp, 284._wp, -297._wp, 0._wp, 232._wp, -249._wp, 69._wp, &
         -297._wp, 0._wp, 47._wp, 150._wp, -154._wp, -75._wp, 95._wp, 0._wp, &
         -16._wp, 88._wp, 69._wp, -48._wp, -1._wp, 21._wp, 0._wp, -83._wp, &
         -27._wp, -2._wp, 20._wp, 17._wp, -23._wp, -7._wp, 0._wp, 8._wp, &
         -19._wp, 5._wp, -23._wp, 11._wp, 14._wp, -15._wp, -11._wp, 0._wp, &
         -21._wp, 15._wp, 9._wp, -6._wp, -6._wp, 9._wp, 9._wp, -7._wp, 2._wp, &
         0._wp, 1._wp, 0._wp, 3._wp, 6._wp, -4._wp, 0._wp, -1._wp, 4._wp, &
         0._wp, -6._wp, 39*0.0_wp/

    DATA g90/0._wp, -29775._wp, -1848._wp, -2131._wp, 3059._wp, 1686._wp, &
         1314._wp, -2239._wp, 1248._wp, 802._wp, 939._wp, 780._wp, 325._wp, &
         -423._wp, 141._wp, -214._wp, 353._wp, 245._wp, -109._wp, -165._wp, &
         -36._wp, 61._wp, 65._wp, 59._wp, -178._wp, 3._wp, 18._wp, -96._wp, &
         77._wp, -64._wp, 2._wp, 26._wp, -1._wp, 5._wp, 9._wp, 0._wp, 23._wp, &
         5._wp, -1._wp, -10._wp, -12._wp, 3._wp, 4._wp, 2._wp, -6._wp, 4._wp, &
         9._wp, 1._wp, -12._wp, 9._wp, -4._wp, -2._wp, 7._wp, 1._wp, -6._wp, &
         -3._wp, -4._wp, 2._wp, -5._wp, -2._wp, 4._wp, 3._wp, 1._wp, 3._wp, &
         3._wp, 0._wp, 39*0._wp/

    DATA h90/0._wp, 0._wp, 5406._wp, 0._wp, -2279._wp, -373._wp, 0._wp, &
         -284._wp, 293._wp, -352._wp, 0._wp, 247._wp, -240._wp, 84._wp, &
         -299._wp, 0._wp, 46._wp, 154._wp, -153._wp, -69._wp, 97._wp, 0._wp, &
         -16._wp, 82._wp, 69._wp, -52._wp, 1._wp, 24._wp, 0._wp, -80._wp, &
         -26._wp, 0._wp, 21._wp, 17._wp, -23._wp, -4._wp, 0._wp, 10._wp, &
         -19._wp, 6._wp, -22._wp, 12._wp, 12._wp, -16._wp, -10._wp, 0._wp, &
         -20._wp, 15._wp, 11._wp, -7._wp, -7._wp, 9._wp, 8._wp, -7._wp, &
         2._wp, 0._wp, 2._wp, 1._wp, 3._wp, 6._wp, -4._wp, 0._wp, -2._wp, &
         3._wp, -1._wp, -6._wp, 39*0._wp/

    DATA g95/0._wp, -29692._wp, -1784._wp, -2200._wp, 3070._wp, 1681._wp, &
         1335._wp, -2267._wp, 1249._wp, 759._wp, 940._wp, 780._wp, 290._wp, &
         -418._wp, 122._wp, -214._wp, 352._wp, 235._wp, -118._wp, -166._wp, &
         -17._wp, 68._wp, 67._wp, 68._wp, -170._wp, -1._wp, 19._wp, -93._wp, &
         77._wp, -72._wp, 1._wp, 28._wp, 5._wp, 4._wp, 8._wp, -2._wp, 25._wp, &
         6._wp, -6._wp, -9._wp, -14._wp, 9._wp, 6._wp, -5._wp, -7._wp, 4._wp, &
         9._wp, 3._wp, -10._wp, 8._wp, -8._wp, -1._wp, 10._wp, -2._wp, &
         -8._wp, -3._wp, -6._wp, 2._wp, -4._wp, -1._wp, 4._wp, 2._wp, 2._wp, &
         5._wp, 1._wp, 0._wp, 39*0._wp/

    DATA h95/0._wp, 0._wp, 5306._wp, 0._wp, -2366._wp, -413._wp, 0._wp, &
         -262._wp, 302._wp, -427._wp, 0._wp, 262._wp, -236._wp, 97._wp, &
         -306._wp, 0._wp, 46._wp, 165._wp, -143._wp, -55._wp, 107._wp, 0._wp, &
         -17._wp, 72._wp, 67._wp, -58._wp, 1._wp, 36._wp, 0._wp, -69._wp, &
         -25._wp, 4._wp, 24._wp, 17._wp, -24._wp, -6._wp, 0._wp, 11._wp, &
         -21._wp, 8._wp, -23._wp, 15._wp, 11._wp, -16._wp, -4._wp, 0._wp, &
         -20._wp, 15._wp, 12._wp, -6._wp, -8._wp, 8._wp, 5._wp, -8._wp, &
         3._wp, 0._wp, 1._wp, 0._wp, 4._wp, 5._wp, -5._wp, -1._wp, -2._wp, &
         1._wp, -2._wp, -7._wp, 39*0._wp/

    DATA g00/0._wp, -29619.4_wp, -1728.2_wp, -2267.7_wp, 3068.4_wp, &
         1670.9_wp, 1339.6_wp, -2288._wp, 1252.1_wp, 714.5_wp, 932.3_wp, &
         786.8_wp, 250._wp, -403._wp, 111.3_wp, -218.8_wp, 351.4_wp, &
         222.3_wp, -130.4_wp, -168.6_wp, -12.9_wp, 72.3_wp, 68.2_wp, 74.2_wp, &
         -160.9_wp, -5.9_wp, 16.9_wp, -90.4_wp, 79.0_wp, -74.0_wp, 0._wp, &
         33.3_wp, 9.1_wp, 6.9_wp, 7.3_wp, -1.2_wp, 24.4_wp, 6.6_wp, -9.2_wp, &
         -7.9_wp, -16.6_wp, 9.1_wp, 7.0_wp, -7.9_wp, -7._wp, 5._wp, 9.4_wp, &
         3._wp, -8.4_wp, 6.3_wp, -8.9_wp, -1.5_wp, 9.3_wp, -4.3_wp, -8.2_wp, &
         -2.6_wp, -6._wp, 1.7_wp, -3.1_wp, -0.5_wp, 3.7_wp, 1._wp, 2._wp, &
         4.2_wp, 0.3_wp, -1.1_wp, 2.7_wp, -1.7_wp, -1.9_wp, 1.5_wp, -0.1_wp, &
         0.1_wp, -0.7_wp, 0.7_wp, 1.7_wp, 0.1_wp, 1.2_wp, 4.0_wp, -2.2_wp, &
         -0.3_wp, 0.2_wp, 0.9_wp, -0.2_wp, 0.9_wp, -0.5_wp, 0.3_wp, -0.3_wp, &
         -0.4_wp, -0.1_wp, -0.2_wp, -0.4_wp, -0.2_wp, -0.9_wp, 0.3_wp, &
         0.1_wp, -0.4_wp, 1.3_wp, -0.4_wp, 0.7_wp, -0.4_wp, 0.3_wp, -0.1_wp, &
         0.4_wp, 0._wp, 0.1_wp/


    DATA h00/0._wp, 0._wp, 5186.1_wp, 0._wp, -2481.6_wp, -458.0_wp, 0._wp, &
         -227.6_wp, 293.4_wp, -491.1_wp, 0._wp, 272.6_wp, -231.9_wp, &
         119.8_wp, -303.8_wp, 0._wp, 43.8_wp, 171.9_wp, -133.1_wp, -39.3_wp, &
         106.3_wp, 0._wp, -17.4_wp, 63.7_wp, 65.1_wp, -61.2_wp, 0.7_wp, &
         43.8_wp, 0._wp, -64.6_wp, -24.2_wp, 6.2_wp, 24._wp, 14.8_wp, &
         -25.4_wp, -5.8_wp, 0.0_wp, 11.9_wp, -21.5_wp, 8.5_wp, -21.5_wp, &
         15.5_wp, 8.9_wp, -14.9_wp, -2.1_wp, 0.0_wp, -19.7_wp, 13.4_wp, &
         12.5_wp, -6.2_wp, -8.4_wp, 8.4_wp, 3.8_wp, -8.2_wp, 4.8_wp, 0.0_wp, &
         1.7_wp, 0.0_wp, 4.0_wp, 4.9_wp, -5.9_wp, -1.2_wp, -2.9_wp, 0.2_wp, &
         -2.2_wp, -7.4_wp, 0.0_wp, 0.1_wp, 1.3_wp, -0.9_wp, -2.6_wp, 0.9_wp, &
         -0.7_wp, -2.8_wp, -0.9_wp, -1.2_wp, -1.9_wp, -0.9_wp, 0.0_wp, &
         -0.4_wp, 0.3_wp, 2.5_wp, -2.6_wp, 0.7_wp, 0.3_wp, 0.0_wp, 0.0_wp, &
         0.3_wp, -0.9_wp, -0.4_wp, 0.8_wp, 0.0_wp, -0.9_wp, 0.2_wp, 1.8_wp, &
         -0.4_wp, -1.0_wp, -0.1_wp, 0.7_wp, 0.3_wp, 0.6_wp, 0.3_wp, -0.2_wp, &
         -0.5_wp, -0.9_wp/


    DATA g05/0._wp, -29556.8_wp, -1671.8_wp, -2340.5_wp, 3047._wp, &
         1656.9_wp, 1335.7_wp, -2305.3_wp, 1246.8_wp, 674.4_wp, 919.8_wp, &
         798.2_wp, 211.5_wp, -379.5_wp, 100.2_wp, -227.6_wp, 354.4_wp, &
         208.8_wp, -136.6_wp, -168.3_wp, -14.1_wp, 72.9_wp, 69.6_wp, 76.6_wp, &
         -151.1_wp, -15.0_wp, 14.7_wp, -86.4_wp, 79.8_wp, -74.4_wp, -1.4_wp, &
         38.6_wp, 12.3_wp, 9.4_wp, 5.5_wp, 2.0_wp, 24.8_wp, 7.7_wp, -11.4_wp, &
         -6.8_wp, -18.0_wp, 10.0_wp, 9.4_wp, -11.4_wp, -5.0_wp, 5.6_wp, &
         9.8_wp, 3.6_wp, -7.0_wp, 5.0_wp, -10.8_wp, -1.3_wp, 8.7_wp, -6.7_wp, &
         -9.2_wp, -2.2_wp, -6.3_wp, 1.6_wp, -2.5_wp, -0.1_wp, 3.0_wp, 0.3_wp, &
         2.1_wp, 3.9_wp, -0.1_wp, -2.2_wp, 2.9_wp, -1.6_wp, -1.7_wp, 1.5_wp, &
         -0.2_wp, 0.2_wp, -0.7_wp, 0.5_wp, 1.8_wp, 0.1_wp, 1.0_wp, 4.1_wp, &
         -2.2_wp, -0.3_wp, 0.3_wp, 0.9_wp, -0.4_wp, 1.0_wp, -0.4_wp, 0.5_wp, &
         -0.3_wp, -0.4_wp, 0.0_wp, -0.4_wp, 0.0_wp, -0.2_wp, -0.9_wp, 0.3_wp, &
         0.3_wp, -0.4_wp, 1.2_wp, -0.4_wp, 0.7_wp, -0.3_wp, 0.4_wp, -0.1_wp, &
         0.4_wp, -0.1_wp, -0.3_wp/

    DATA h05/0._wp, 0.0_wp, 5080.0_wp, 0.0_wp, -2594.9_wp, -516.7_wp, &
         0.0_wp, -200.4_wp, 269.3_wp, -524.5_wp, 0.0_wp, 281.4_wp, -225.8_wp, &
         145.7_wp, -304.7_wp, 0.0_wp, 42.7_wp, 179.8_wp, -123.0_wp, -19.5_wp, &
         103.6_wp, 0.0_wp, -20.2_wp, 54.7_wp, 63.7_wp, -63.4_wp, 0.0_wp, &
         50.3_wp, 0.0_wp, -61.4_wp, -22.5_wp, 6.9_wp, 25.4_wp, 10.9_wp, &
         -26.4_wp, -4.8_wp, 0.0_wp, 11.2_wp, -21.0_wp, 9.7_wp, -19.8_wp, &
         16.1_wp, 7.7_wp, -12.8_wp, -0.1_wp, 0.0_wp, -20.1_wp, 12.9_wp, &
         12.7_wp, -6.7_wp, -8.1_wp, 8.1_wp, 2.9_wp, -7.9_wp, 5.9_wp, 0.0_wp, &
         2.4_wp, 0.2_wp, 4.4_wp, 4.7_wp, -6.5_wp, -1.0_wp, -3.4_wp, -0.9_wp, &
         -2.3_wp, -8.0_wp, 0.0_wp, 0.3_wp, 1.4_wp, -0.7_wp, -2.4_wp, 0.9_wp, &
         -0.6_wp, -2.7_wp, -1.0_wp, -1.5_wp, -2.0_wp, -1.4_wp, 0.0_wp, &
         -0.5_wp, 0.3_wp, 2.3_wp, -2.7_wp, 0.6_wp, 0.4_wp, 0.0_wp, 0.0_wp, &
         0.3_wp, -0.8_wp, -0.4_wp, 1.0_wp, 0.0_wp, -0.7_wp, 0.3_wp, 1.7_wp, &
         -0.5_wp, -1.0_wp, 0.0_wp, 0.7_wp, 0.2_wp, 0.6_wp, 0.4_wp, -0.2_wp, &
         -0.5_wp, -1.0_wp/

    DATA dg05/0.0_wp, 8.8_wp, 10.8_wp, -15.0_wp, -6.9_wp, -1.0_wp, &
         -0.3_wp, -3.1_wp, -0.9_wp, -6.8_wp, -2.5_wp, 2.8_wp, -7.1_wp, &
         5.9_wp, -3.2_wp, -2.6_wp, 0.4_wp, -3.0_wp, -1.2_wp, 0.2_wp, -0.6_wp, &
         -0.8_wp, 0.2_wp, -0.2_wp, 2.1_wp, -2.1_wp, -0.4_wp, 1.3_wp, -0.4_wp, &
         0.0_wp, -0.2_wp, 1.1_wp, 0.6_wp, 0.4_wp, -0.5_wp, 0.9_wp, -0.2_wp, &
         0.2_wp, -0.2_wp, 0.2_wp, -0.2_wp, 0.2_wp, 0.5_wp, -0.7_wp, 0.5_wp/

    DATA dh05/0.0_wp, 0.0_wp, -21.3_wp, 0.0_wp, -23.3_wp, -14.0_wp, &
         0.0_wp, 5.4_wp, -6.5_wp, -2.0_wp, 0.0_wp, 2.0_wp, 1.8_wp, 5.6_wp, &
         0.0_wp, 0.0_wp, 0.1_wp, 1.8_wp, 2.0_wp, 4.5_wp, -1.0_wp, 0.0_wp, &
         -0.4_wp, -1.9_wp, -0.4_wp, -0.4_wp, -0.2_wp, 0.9_wp, 0.0_wp, 0.8_wp, &
         0.4_wp, 0.1_wp, 0.2_wp, -0.9_wp, -0.3_wp, 0.3_wp, 0.0_wp, -0.2_wp, &
         0.2_wp, 0.2_wp, 0.4_wp, 0.2_wp, -0.3_wp, 0.5_wp, 0.4_wp/


    iy = iyear

    ! WE ARE RESTRICTED BY THE INTERVAL 1965-2010, FOR WHICH THE IGRF COEFFICIENTS
    !   ARE KNOWN; IF IYEAR IS OUTSIDE THIS INTERVAL, THEN THE SUBROUTINE USES THE
    !     NEAREST LIMITING VALUE AND PRINTS A WARNING:

    IF (iy<1965) THEN
      iy = 1965
      WRITE (*,'(/,/,1x,a,i4,/,6x,a,i4,/)') &
           '**** RECALC WARNS: YEAR IS OUT OF INTERVAL 1965-2010: IYEAR=', &
           iyear, 'CALCULATIONS WILL BE DONE FOR IYEAR=', iy
    END IF

    IF (iy>2010) THEN
      iy = 2010
      WRITE (*,'(/,/,1x,a,i4,/,6x,a,i4,/)') &
           '**** RECALC WARNS: YEAR IS OUT OF INTERVAL 1965-2010: IYEAR=', &
           iyear, 'CALCULATIONS WILL BE DONE FOR IYEAR=', iy
    END IF


    ! CALCULATE THE ARRAY REC, CONTAINING COEFFICIENTS FOR THE RECURSION RELATIONS,
    ! USED IN THE IGRF SUBROUTINE FOR CALCULATING THE ASSOCIATE LEGENDRE POLYNOMIALS
    ! AND THEIR DERIVATIVES:

    CONSTRUCT_1: DO n = 1, 14
      n2 = 2*n - 1
      n2 = n2*(n2-2)
      DO m = 1, n
        mn = n*(n-1)/2 + m
        rec(mn) = real((n-m)*(n+m-2),kind=wp)/real(n2,kind=wp)
      END DO
    END DO CONSTRUCT_1

    IF (iy<1970) THEN !INTERPOLATE BETWEEN 1965 - 1970

      !      INTERPOLATE BETWEEEN 1965 - 1970:

      f2 = (real(iy,kind=wp)+real(iday-1,kind=wp)/365.25_wp-1965._wp)/ &
           5._wp
      f1 = 1._wp - f2
      DO n = 1, 105
        g(n) = g65(n)*f1 + g70(n)*f2
        h(n) = h65(n)*f1 + h70(n)*f2
      END DO
    ELSE IF (iy<1975) THEN !INTERPOLATE BETWEEN 1970 - 1975

      !      INTERPOLATE BETWEEN 1970 - 1975:

      f2 = (real(iy,kind=wp)+real(iday-1,kind=wp)/365.25_wp-1970._wp)/ &
           5._wp
      f1 = 1._wp - f2
      DO n = 1, 105
        g(n) = g70(n)*f1 + g75(n)*f2
        h(n) = h70(n)*f1 + h75(n)*f2
      END DO
    ELSE IF (iy<1980) THEN !INTERPOLATE BETWEEN 1975 - 1980

      !      INTERPOLATE BETWEEN 1975 - 1980:

      f2 = (real(iy,kind=wp)+real(iday-1,kind=wp)/365.25_wp-1975._wp)/ &
           5._wp
      f1 = 1._wp - f2
      DO n = 1, 105
        g(n) = g75(n)*f1 + g80(n)*f2
        h(n) = h75(n)*f1 + h80(n)*f2
      END DO
    ELSE IF (iy<1985) THEN !INTERPOLATE BETWEEN 1980 - 1985

      !      INTERPOLATE BETWEEN 1980 - 1985:

      f2 = (real(iy,kind=wp)+real(iday-1,kind=wp)/365.25_wp-1980._wp)/ &
           5._wp
      f1 = 1._wp - f2
      DO n = 1, 105
        g(n) = g80(n)*f1 + g85(n)*f2
        h(n) = h80(n)*f1 + h85(n)*f2
      END DO
    ELSE IF (iy<1990) THEN !INTERPOLATE BETWEEN 1985 - 1990

      !      INTERPOLATE BETWEEN 1985 - 1990:

      f2 = (real(iy,kind=wp)+real(iday-1,kind=wp)/365.25_wp-1985._wp)/ &
           5._wp
      f1 = 1._wp - f2
      DO n = 1, 105
        g(n) = g85(n)*f1 + g90(n)*f2
        h(n) = h85(n)*f1 + h90(n)*f2
      END DO
    ELSE IF (iy<1995) THEN !INTERPOLATE BETWEEN 1990 - 1995

      !      INTERPOLATE BETWEEN 1990 - 1995:

      f2 = (real(iy,kind=wp)+real(iday-1,kind=wp)/365.25_wp-1990._wp)/ &
           5._wp
      f1 = 1._wp - f2
      DO n = 1, 105
        g(n) = g90(n)*f1 + g95(n)*f2
        h(n) = h90(n)*f1 + h95(n)*f2
      END DO
    ELSE IF (iy<2000) THEN !INTERPOLATE BETWEEN 1995 - 2000

      !      INTERPOLATE BETWEEN 1995 - 2000:

      f2 = (real(iy,kind=wp)+real(iday-1,kind=wp)/365.25_wp-1995._wp)/ &
           5._wp
      f1 = 1._wp - f2
      DO n = 1, 105 !  THE 2000 COEFFICIENTS (G00) GO THROUGH THE ORDER 13, NOT 10
        g(n) = g95(n)*f1 + g00(n)*f2
        h(n) = h95(n)*f1 + h00(n)*f2
      END DO
    ELSE IF (iy<2005) THEN !INTERPOLATE BETWEEN 2000 - 2005

      !      INTERPOLATE BETWEEN 2000 - 2005:

      f2 = (real(iy,kind=wp)+real(iday-1,kind=wp)/365.25_wp-2000._wp)/ &
           5._wp
      f1 = 1._wp - f2
      DO n = 1, 105
        g(n) = g00(n)*f1 + g05(n)*f2
        h(n) = h00(n)*f1 + h05(n)*f2
      END DO
    ELSE

      !      EXTRAPOLATE BEYOND 2005:


      dt = real(iy,kind=wp) + real(iday-1,kind=wp)/365.25_wp - 2005._wp
      DO n = 1, 105
        g(n) = g05(n)
        h(n) = h05(n)
        IF (n>45) CYCLE
        g(n) = g(n) + dg05(n)*dt
        h(n) = h(n) + dh05(n)*dt
      END DO
    END IF

    !  COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
    !  THEM BY SCHMIDT NORMALIZATION FACTORS:

    s = 1._wp
    CONSTRUCT_2: DO n = 2, 14
      mn = n*(n-1)/2 + 1
      s = s*real(2*n-3,kind=wp)/real(n-1,kind=wp)
      g(mn) = g(mn)*s
      h(mn) = h(mn)*s
      p = s
      DO m = 2, n
        aa = 1._wp
        IF (m==2) aa = 2._wp
        p = p*sqrt(aa*real(n-m+1,kind=wp)/real(n+m-2,kind=wp))
        mnn = mn + m - 1
        g(mnn) = g(mnn)*p
        h(mnn) = h(mnn)*p
      END DO
    END DO CONSTRUCT_2

    g10 = -g(2)
    g11 = g(3)
    h11 = h(3)

    ! NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD.SYSTEM:
    !  SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
    !        ST0 * CL0                ST0 * SL0                CT0

    sq = g11**2.0_wp + h11**2.0_wp
    sqq = sqrt(sq)
    sqr = sqrt(g10**2.0_wp+sq)
    sl0 = -h11/sqq
    cl0 = -g11/sqq
    st0 = sqq/sqr
    ct0 = g10/sqr
    stcl = st0*cl0
    stsl = st0*sl0
    ctsl = ct0*sl0
    ctcl = ct0*cl0

    CALL sun(iy,iday,ihour,min,isec,gst,slong,srasn,sdec)

    ! S1,S2, AND S3 ARE THE COMPONENTS OF THE UNIT VECTOR EXGSM=EXGSE IN THE
    !  SYSTEM GEI POINTING FROM THE EARTH'S CENTER TO THE SUN:

    s1 = cos(srasn)*cos(sdec)
    s2 = sin(srasn)*cos(sdec)
    s3 = sin(sdec)
    cgst = cos(gst)
    sgst = sin(gst)

    ! DIP1, DIP2, AND DIP3 ARE THE COMPONENTS OF THE UNIT VECTOR EZSM=EZMAG
    !  IN THE SYSTEM GEI:

    dip1 = stcl*cgst - stsl*sgst
    dip2 = stcl*sgst + stsl*cgst
    dip3 = ct0

    ! NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EYGSM IN THE SYSTEM GEI
    !  BY TAKING THE VECTOR PRODUCT D x S AND NORMALIZING IT TO UNIT LENGTH:

    y1 = dip2*s3 - dip3*s2
    y2 = dip3*s1 - dip1*s3
    y3 = dip1*s2 - dip2*s1
    y = sqrt(y1*y1+y2*y2+y3*y3)
    y1 = y1/y
    y2 = y2/y
    y3 = y3/y

    !  THEN IN THE GEI SYSTEM THE UNIT VECTOR Z = EZGSM = EXGSM x EYGSM = S x Y
    !   HAS THE COMPONENTS:

    z1 = s2*y3 - s3*y2
    z2 = s3*y1 - s1*y3
    z3 = s1*y2 - s2*y1

    !   THE VECTOR EZGSE (HERE DZ) IN GEI HAS THE COMPONENTS (0,-SIN(DELTA),
    !    COS(DELTA)) = (0.,-0.397823,0.917462); HERE DELTA = 23.44214 DEG FOR
    !  THE EPOCH 1978 (SEE THE BOOK BY GUREVICH OR OTHER ASTRONOMICAL HANDBOOKS).
    !   HERE THE MOST ACCURATE TIME-DEPENDENT FORMULA IS USED:

    dj = real(365*(iy-1900)+(iy-1901)/4+iday,kind=wp) - 0.5_wp + &
         real(ihour*3600+min*60+isec,kind=wp)/86400._wp
    t = dj/36525._wp
    obliq = (23.45229_wp-0.0130125_wp*t)/57.2957795_wp
    dz1 = 0._wp
    dz2 = -sin(obliq)
    dz3 = cos(obliq)

    ! THEN THE UNIT VECTOR EYGSE IN GEI SYSTEM IS THE VECTOR PRODUCT DZ x S :

    dy1 = dz2*s3 - dz3*s2
    dy2 = dz3*s1 - dz1*s3
    dy3 = dz1*s2 - dz2*s1

    !  THE ELEMENTS OF THE MATRIX GSE TO GSM ARE THE SCALAR PRODUCTS:
    !  CHI=EM22=(EYGSM,EYGSE), SHI=EM23=(EYGSM,EZGSE), EM32=(EZGSM,EYGSE)=-EM23,
    !    AND EM33=(EZGSM,EZGSE)=EM22

    chi = y1*dy1 + y2*dy2 + y3*dy3
    shi = y1*dz1 + y2*dz2 + y3*dz3
    hi = asin(shi)

    !   TILT ANGLE: PSI=ARCSIN(DIP,EXGSM)

    sps = dip1*s1 + dip2*s2 + dip3*s3
    cps = sqrt(1._wp-sps**2.0_wp)
    psi = asin(sps)

    !   THE ELEMENTS OF THE MATRIX MAG TO SM ARE THE SCALAR PRODUCTS:
    ! CFI=GM22=(EYSM,EYMAG), SFI=GM23=(EYSM,EXMAG); THEY CAN BE DERIVED AS FOLLOWS:

    ! IN GEO THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
    ! AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THE COMPONENTS ARE:
    ! EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
    !           CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
    !           -ST0
    ! EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
    !           -SL0*SIN(GST)+CL0*COS(GST)
    !            0
    ! THE COMPONENTS OF EYSM IN GEI WERE FOUND ABOVE AS Y1, Y2, AND Y3;
    ! NOW WE ONLY HAVE TO COMBINE THE QUANTITIES INTO SCALAR PRODUCTS:

    exmagx = ct0*(cl0*cgst-sl0*sgst)
    exmagy = ct0*(cl0*sgst+sl0*cgst)
    exmagz = -st0
    eymagx = -(sl0*cgst+cl0*sgst)
    eymagy = -(sl0*sgst-cl0*cgst)
    cfi = y1*eymagx + y2*eymagy
    sfi = y1*exmagx + y2*exmagy + y3*exmagz

    xmut = (atan2(sfi,cfi)+pi)*3.8197186342_wp

    ! THE ELEMENTS OF THE MATRIX GEO TO GSM ARE THE SCALAR PRODUCTS:

    !  A11=(EXGEO,EXGSM), A12=(EYGEO,EXGSM), A13=(EZGEO,EXGSM),
    !  A21=(EXGEO,EYGSM), A22=(EYGEO,EYGSM), A23=(EZGEO,EYGSM),
    !  A31=(EXGEO,EZGSM), A32=(EYGEO,EZGSM), A33=(EZGEO,EZGSM),

    !  ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:

    ! EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
    ! EXGSM=(S1,S2,S3),  EYGSM=(Y1,Y2,Y3),   EZGSM=(Z1,Z2,Z3)
    !                                                          AND  THEREFORE:

    a11 = s1*cgst + s2*sgst
    a12 = -s1*sgst + s2*cgst
    a13 = s3
    a21 = y1*cgst + y2*sgst
    a22 = -y1*sgst + y2*cgst
    a23 = y3
    a31 = z1*cgst + z2*sgst
    a32 = -z1*sgst + z2*cgst
    a33 = z3



  END SUBROUTINE recalc


  !=============================================================================


  SUBROUTINE geomag(xgeo,ygeo,zgeo,xmag,ymag,zmag,j)

    !   CONVERTS GEOGRAPHIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICA VERSA.

    !                   J>0                       J<0
    !-----INPUT:  J,XGEO,YGEO,ZGEO           J,XMAG,YMAG,ZMAG
    !-----OUTPUT:    XMAG,YMAG,ZMAG           XGEO,YGEO,ZGEO


    ! ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEOMAG IN TWO CASES:
    !    /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
    !    /B/  IF THE VALUES OF IYEAR AND/OR IDAY HAVE BEEN CHANGED


    !  LAST MOFIFICATION:  MARCH 30, 2003 (INVOCATION OF RECALC INSIDE THIS S/R WAS REMOVED)

    !  AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: xgeo, ygeo, zgeo, xmag, ymag, zmag

    INTEGER :: j

    REAL (wp) :: st0, ct0, sl0, cl0, ctcl, stcl, ctsl, stsl, ab(19), bb(8)

    COMMON /geopack1/st0, ct0, sl0, cl0, ctcl, stcl, ctsl, stsl, ab, bb

    IF (j>0) THEN
      xmag = xgeo*ctcl + ygeo*ctsl - zgeo*st0
      ymag = ygeo*cl0 - xgeo*sl0
      zmag = xgeo*stcl + ygeo*stsl + zgeo*ct0
    ELSE
      xgeo = xmag*ctcl - ymag*sl0 + zmag*stcl
      ygeo = xmag*ctsl + ymag*cl0 + zmag*stsl
      zgeo = zmag*ct0 - xmag*st0
    END IF


  END SUBROUTINE geomag


  !=============================================================================


  SUBROUTINE geigeo(xgei,ygei,zgei,xgeo,ygeo,zgeo,j)

    !  CONVERTS EQUATORIAL INERTIAL (GEI) TO GEOGRAPHICAL (GEO) COORDS
    !  OR VICA VERSA.
    !                   J>0                J<0
    !----INPUT:  J,XGEI,YGEI,ZGEI    J,XGEO,YGEO,ZGEO
    !----OUTPUT:   XGEO,YGEO,ZGEO      XGEI,YGEI,ZGEI

    ! ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEIGEO IN TWO CASES:
    !    /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
    !    /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED

    !    LAST MODIFICATION:  MARCH 31, 2003

    !    AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: xgei, ygei, zgei, xgeo, ygeo, zgeo

    INTEGER :: j

    REAL (wp) :: a(27), cgst, sgst, b(6)

    COMMON /geopack1/a, cgst, sgst, b

    IF (j>0) THEN
      xgeo = xgei*cgst + ygei*sgst
      ygeo = ygei*cgst - xgei*sgst
      zgeo = zgei
    ELSE
      xgei = xgeo*cgst - ygeo*sgst
      ygei = ygeo*cgst + xgeo*sgst
      zgei = zgeo
    END IF


  END SUBROUTINE geigeo


  !=============================================================================


  SUBROUTINE magsm(xmag,ymag,zmag,xsm,ysm,zsm,j)

    ! CONVERTS DIPOLE (MAG) TO SOLAR MAGNETIC (SM) COORDINATES OR VICA VERSA

    !                   J>0              J<0
    !----INPUT: J,XMAG,YMAG,ZMAG     J,XSM,YSM,ZSM
    !----OUTPUT:    XSM,YSM,ZSM       XMAG,YMAG,ZMAG

    ! ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE MAGSM IN TWO CASES:
    !    /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
    !    /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED

    !    LAST MODIFICATION:  MARCH 31, 2003

    !    AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: xmag, ymag, zmag, xsm, ysm, zsm

    INTEGER :: j

    REAL (wp) :: a(8), sfi, cfi, b(7), ab(10), ba(8)

    COMMON /geopack1/a, sfi, cfi, b, ab, ba

    IF (j>0) THEN
      xsm = xmag*cfi - ymag*sfi
      ysm = xmag*sfi + ymag*cfi
      zsm = zmag
    ELSE
      xmag = xsm*cfi + ysm*sfi
      ymag = ysm*cfi - xsm*sfi
      zmag = zsm
    END IF


  END SUBROUTINE magsm


  !=============================================================================


  SUBROUTINE gsmgse(xgsm,ygsm,zgsm,xgse,ygse,zgse,j)

    ! CONVERTS GEOCENTRIC SOLAR MAGNETOSPHERIC (GSM) COORDS TO SOLAR ECLIPTIC (GSE) ONES
    !  OR VICA VERSA.
    !                   J>0                J<0
    !-----INPUT: J,XGSM,YGSM,ZGSM    J,XGSE,YGSE,ZGSE
    !----OUTPUT:   XGSE,YGSE,ZGSE      XGSM,YGSM,ZGSM

    ! ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GSMGSE IN TWO CASES:
    !    /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
    !    /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED

    !    LAST MODIFICATION:  MARCH 31, 2003

    !    AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: xgsm, ygsm, zgsm, xgse, ygse, zgse

    INTEGER :: j

    REAL (wp) :: a(12), shi, chi, ab(13), ba(8)

    COMMON /geopack1/a, shi, chi, ab, ba

    IF (j>0) THEN
      xgse = xgsm
      ygse = ygsm*chi - zgsm*shi
      zgse = ygsm*shi + zgsm*chi
    ELSE
      xgsm = xgse
      ygsm = ygse*chi + zgse*shi
      zgsm = zgse*chi - ygse*shi
    END IF


  END SUBROUTINE gsmgse


  !=============================================================================


  SUBROUTINE smgsm(xsm,ysm,zsm,xgsm,ygsm,zgsm,j)

    ! CONVERTS SOLAR MAGNETIC (SM) TO GEOCENTRIC SOLAR MAGNETOSPHERIC
    !  (GSM) COORDINATES OR VICA VERSA.
    !                 J>0                 J<0
    !-----INPUT: J,XSM,YSM,ZSM        J,XGSM,YGSM,ZGSM
    !----OUTPUT:  XGSM,YGSM,ZGSM       XSM,YSM,ZSM

    ! ATTENTION:  SUBROUTINE RECALC  MUST BE INVOKED BEFORE SMGSM IN TWO CASES:
    !    /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
    !    /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED

    !    LAST MODIFICATION:  MARCH 31, 2003

    !    AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: xsm, ysm, zsm, xgsm, ygsm, zgsm

    INTEGER :: j

    REAL (wp) :: a(10), sps, cps, b(15), ab(8)

    COMMON /geopack1/a, sps, cps, b, ab

    IF (j>0) THEN
      xgsm = xsm*cps + zsm*sps
      ygsm = ysm
      zgsm = zsm*cps - xsm*sps
    ELSE
      xsm = xgsm*cps - zgsm*sps
      ysm = ygsm
      zsm = xgsm*sps + zgsm*cps
    END IF


  END SUBROUTINE smgsm


  !=============================================================================


  SUBROUTINE geogsm(xgeo,ygeo,zgeo,xgsm,ygsm,zgsm,j)

    ! CONVERTS GEOGRAPHIC (GEO) TO GEOCENTRIC SOLAR MAGNETOSPHERIC (GSM) COORDINATES
    !  OR VICA VERSA.

    !                  J>0                   J<0
    !----- INPUT:  J,XGEO,YGEO,ZGEO    J,XGSM,YGSM,ZGSM
    !---- OUTPUT:    XGSM,YGSM,ZGSM      XGEO,YGEO,ZGEO

    ! ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEOGSM IN TWO CASES:
    !    /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
    !    /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC  HAVE BEEN CHANGED

    !    LAST MODIFICATION: MARCH 31, 2003

    !    AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: xgeo, ygeo, zgeo, xgsm, ygsm, zgsm

    INTEGER :: j

    REAL (wp) :: aa(17), a11, a21, a31, a12, a22, a32, a13, a23, a33, d, &
         b(8)

    COMMON /geopack1/aa, a11, a21, a31, a12, a22, a32, a13, a23, a33, d, b

    IF (j>0) THEN
      xgsm = a11*xgeo + a12*ygeo + a13*zgeo
      ygsm = a21*xgeo + a22*ygeo + a23*zgeo
      zgsm = a31*xgeo + a32*ygeo + a33*zgeo
    ELSE
      xgeo = a11*xgsm + a21*ygsm + a31*zgsm
      ygeo = a12*xgsm + a22*ygsm + a32*zgsm
      zgeo = a13*xgsm + a23*ygsm + a33*zgsm
    END IF


  END SUBROUTINE geogsm


  !=============================================================================


  SUBROUTINE rhand(x,y,z,r1,r2,r3,iopt,parmod,exname,inname)

    ! CALCULATES THE COMPONENTS OF THE RIGHT HAND SIDE VECTOR IN THE GEOMAGNETIC FIELD
    !   LINE EQUATION  (a subsidiary subroutine for the subroutine STEP)

    !    LAST MODIFICATION:  MARCH 31, 2003

    !    AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: x, y, z, r1, r2, r3, parmod(10)

    INTEGER :: iopt

    REAL (wp) :: a(15), psi, aa(10), ds3, bb(8)

    EXTERNAL exname, inname

    !    EXNAME AND INNAME ARE NAMES OF SUBROUTINES FOR THE EXTERNAL AND INTERNAL
    !    PARTS OF THE TOTAL FIELD

    COMMON /geopack1/a, psi, aa, ds3, bb

    REAL (wp) :: bx, bxgsm, hxgsm, by, bygsm, hygsm, bz, bzgsm, hzgsm, b

    CALL exname(iopt,parmod,psi,x,y,z,bxgsm,bygsm,bzgsm)
    CALL inname(x,y,z,hxgsm,hygsm,hzgsm)

    bx = bxgsm + hxgsm
    by = bygsm + hygsm
    bz = bzgsm + hzgsm
    b = ds3/sqrt(bx**2.0_wp+by**2.0_wp+bz**2.0_wp)
    r1 = bx*b
    r2 = by*b
    r3 = bz*b


  END SUBROUTINE rhand


  !=============================================================================


  SUBROUTINE step(x,y,z,ds,errin,iopt,parmod,exname,inname)

    !  RE-CALCULATES {X,Y,Z}, MAKING A STEP ALONG A FIELD LINE.
    !  DS IS THE STEP SIZE, ERRIN IS PERMISSIBLE ERROR VALUE, IOPT SPECIFIES THE EXTERNAL
    !  MODEL VERSION, THE ARRAY PARMOD CONTAINS INPUT PARAMETERS FOR THAT MODEL
    !  EXNAME IS THE NAME OF THE EXTERNAL FIELD SUBROUTINE
    !  INNAME IS THE NAME OF THE INTERNAL FIELD SUBROUTINE (EITHER DIP OR IGRF)

    !  ALL THE PARAMETERS ARE INPUT ONES; OUTPUT IS THE RENEWED TRIPLET X,Y,Z

    !    LAST MODIFICATION:  MARCH 31, 2003

    !    AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: x, y, z, ds, errin, parmod(10)

    INTEGER :: iopt

    REAL (wp) :: a(26), ds3, b(8)

    COMMON /geopack1/a, ds3, b

    REAL (wp) :: r11, r12, r13, r21, r22, r23, r31, r32, r33, r41, r42, &
         r43, r51, r52, r53, errcur

    EXTERNAL exname, inname

    DO

      ds3 = -ds/3._wp
      CALL rhand(x,y,z,r11,r12,r13,iopt,parmod,exname,inname)
      CALL rhand(x+r11,y+r12,z+r13,r21,r22,r23,iopt,parmod,exname,inname)
      CALL rhand(x+.5_wp*(r11+r21),y+.5_wp*(r12+r22),z+.5_wp*(r13+r23), &
           r31,r32,r33,iopt,parmod,exname,inname)
      CALL rhand(x+.375_wp*(r11+3._wp*r31),y+.375_wp*(r12+3._wp*r32), &
           z+.375_wp*(r13+3._wp*r33),r41,r42,r43,iopt,parmod,exname,inname)
      CALL rhand(x+1.5_wp*(r11-3._wp*r31+4._wp*r41), &
           y+1.5_wp*(r12-3._wp*r32+4._wp*r42),z+1.5_wp*(r13-3._wp*r33+4._wp* &
           r43),r51,r52,r53,iopt,parmod,exname,inname)
      errcur = abs(r11-4.5_wp*r31+4._wp*r41-.5_wp*r51) + &
           abs(r12-4.5_wp*r32+4._wp*r42-.5_wp*r52) + &
           abs(r13-4.5_wp*r33+4._wp*r43-.5_wp*r53)
      IF (errcur<errin) THEN
        EXIT
      ELSE
        ds = ds*.5_wp
      END IF
    END DO
    x = x + .5_wp*(r11+4._wp*r41+r51)
    y = y + .5_wp*(r12+4._wp*r42+r52)
    z = z + .5_wp*(r13+4._wp*r43+r53)
    IF (errcur<errin*.04_wp .AND. abs(ds)<1.33_wp) ds = ds*1.5_wp


  END SUBROUTINE step


  !=============================================================================


  SUBROUTINE trace(xi,yi,zi,dir,rlim,r0,iopt,parmod,exname,inname,xf,yf, &
       zf,xx,yy,zz,l)

    ! TRACES A FIELD LINE FROM AN ARBITRARY POINT OF SPACE TO THE EARTH'S
    ! SURFACE OR TO A MODEL LIMITING BOUNDARY.

    ! THE HIGHEST ORDER OF SPHERICAL HARMONICS IN THE MAIN FIELD EXPANSION USED
    ! IN THE MAPPING IS CALCULATED AUTOMATICALLY. IF INNAME=IGRF_GSM, THEN AN IGRF MODEL
    ! FIELD WILL BE USED, AND IF INNAME=DIP, A PURE DIPOLE FIELD WILL BE USED.

    ! IN ANY CASE, BEFORE CALLING TRACE, ONE SHOULD INVOKE RECALC, TO CALCULATE CORRECT
    ! VALUES OF THE IGRF COEFFICIENTS AND ALL QUANTITIES NEEDED FOR TRANSFORMATIONS
    ! BETWEEN COORDINATE SYSTEMS INVOLVED IN THIS CALCULATIONS.

    ! ALTERNATIVELY, THE SUBROUTINE RECALC CAN BE INVOKED WITH THE DESIRED VALUES OF
    ! IYEAR AND IDAY (TO SPECIFY THE DIPOLE MOMENT), WHILE THE VALUES OF THE DIPOLE
    ! TILT ANGLE PSI (IN RADIANS) AND ITS SINE (SPS) AND COSINE (CPS) CAN BE EXPLICITLY
    ! SPECIFIED AND FORWARDED TO THE COMMON BLOCK GEOPACK1 (11th, 12th, AND 16th ELEMENTS, RESP.)

    !------------- INPUT PARAMETERS:

    !  XI,YI,ZI - GSM COORDS OF INITIAL POINT (IN EARTH RADII, 1 RE = 6371.2 km),

    !  DIR - SIGN OF THE TRACING DIRECTION: IF DIR=1.0 THEN WE MOVE ANTIPARALLEL TO THE
    !    FIELD VECTOR (E.G. FROM NORTHERN TO SOUTHERN CONJUGATE POINT),
    !    AND IF DIR=-1.0 THEN THE TRACING GOES IN THE OPPOSITE DIRECTION.

    !  R0 -  RADIUS OF A SPHERE (IN RE) FOR WHICH THE FIELD LINE ENDPOINT COORDINATES
    !    XF,YF,ZF  SHOULD BE CALCULATED.

    !  RLIM - UPPER LIMIT OF THE GEOCENTRIC DISTANCE, WHERE THE TRACING IS TERMINATED.

    !  IOPT - A MODEL INDEX; CAN BE USED FOR SPECIFYING AN OPTION OF THE EXTERNAL FIELD
    !      MODEL (E.G., INTERVAL OF THE KP-INDEX). ALTERNATIVELY, ONE CAN USE THE ARRAY
    !      PARMOD FOR THAT PURPOSE (SEE BELOW); IN THAT CASE IOPT IS JUST A DUMMY PARAMETER.

    !  PARMOD -  A 10-ELEMENT ARRAY CONTAINING MODEL PARAMETERS, NEEDED FOR A UNIQUE
    !     SPECIFICATION OF THE EXTERNAL FIELD. THE CONCRETE MEANING OF THE COMPONENTS
    !     OF PARMOD DEPENDS ON A SPECIFIC VERSION OF THE EXTERNAL FIELD MODEL.

    !  EXNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE EXTERNAL MAGNETIC FIELD
    !   (E.G., T96_01).
    !  INNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE INTERNAL MAGNETIC FIELD
    !   (EITHER DIP OR IGRF_GSM).

    !-------------- OUTPUT PARAMETERS:

    !  XF,YF,ZF - GSM COORDS OF THE LAST CALCULATED POINT OF A FIELD LINE
    !  XX,YY,ZZ - ARRAYS, CONTAINING COORDS OF FIELD LINE POINTS. HERE THEIR MAXIMAL LENGTH WAS
    !     ASSUMED EQUAL TO 999.
    !  L - ACTUAL NUMBER OF THE CALCULATED FIELD LINE POINTS. IF L EXCEEDS 999, TRACING
    !    TERMINATES, AND A WARNING IS DISPLAYED.


    !    LAST MODIFICATION:  MARCH 31, 2003.

    !    AUTHOR:  N. A. TSYGANENKO

    IMPLICIT NONE

    REAL (wp) :: xi, yi, zi, dir, rlim, r0, parmod(10), xf, yf, zf

    REAL (wp) :: xx(1000), yy(1000), zz(1000)

    INTEGER :: iopt, l

    REAL (wp) :: aa(26), dd, bb(8)

    COMMON /geopack1/aa, dd, bb

    REAL (wp) :: err, ds, x, y, z, al, r1, r2, r3, ad, rr, ryz, r, fc, xr, &
         yr, zr

    EXTERNAL exname, inname

    err = 0.0001_wp
    l = 0
    ds = 0.5_wp*dir
    x = xi
    y = yi
    z = zi
    dd = dir
    al = 0._wp

    ! here we call RHAND just to find out the sign of the radial component of the field
    !  vector, and to determine the initial direction of the tracing (i.e., either away
    !  or towards Earth):

    CALL rhand(x,y,z,r1,r2,r3,iopt,parmod,exname,inname)
    ad = 0.01_wp
    IF (x*r1+y*r2+z*r3<0._wp) ad = -0.01_wp

    !    |AD|=0.01 and its sign follows the rule:
    ! (1) if DIR=1 (tracing antiparallel to B vector) then the sign of AD is the same as of Br
    ! (2) if DIR=-1 (tracing parallel to B vector) then the sign of AD is opposite to that of Br
    !    AD is defined in order to initialize the value of RR (radial distance at previous step):

    rr = sqrt(x**2.0_wp+y**2.0_wp+z**2.0_wp) + ad
    DO
      l = l + 1
      IF (l>999_dp) THEN
        GO TO 6000
      ELSE
        xx(l) = x
        yy(l) = y
        zz(l) = z
        ryz = y**2.0_wp + z**2.0_wp
        r2 = x**2.0_wp + ryz
        r = sqrt(r2)

        ! check if the line hit the outer tracing boundary; if yes, then terminate
        !  the tracing (label 8):

        IF (r>rlim .OR. ryz>1600._wp .OR. x>20._wp) THEN
          GO TO 6010
        ELSE

          ! check whether or not the inner tracing boundary was crossed from outside,
          ! if yes, then calculate the footpoint position by interpolation (go to label 6):

          IF (r<r0 .AND. rr>r) THEN
            EXIT
          ELSE

            ! check if (i) we are moving outward, or (ii) we are still sufficiently
            !   far from Earth (beyond R=5Re); if yes, proceed further:

            IF (r<rr .AND. r<=5._wp) THEN

              ! now we moved closer inward (between R=3 and R=5); go to 3 and begin logging
              ! previous values of X,Y,Z, to be used in the interpolation (after having
              ! crossed the inner tracing boundary):

              IF (r>=3._wp) THEN
                ds = dir
              ELSE

                ! we entered inside the sphere R=3: to avoid too large steps (and hence inaccurate
                ! interpolated position of the footpoint), enforce the progressively smaller
                ! stepsize values as we approach the inner boundary R=R0:

                fc = 0.2_wp
                IF (r-r0<0.05_wp) fc = 0.05_wp
                al = fc*(r-r0+0.2_wp)
                ds = dir*al
              END IF
              xr = x
              yr = y
              zr = z
            END IF
            rr = r
            CALL step(x,y,z,ds,err,iopt,parmod,exname,inname)
          END IF
        END IF
      END IF
    END DO

    ! find the footpoint position by interpolating between the current and previous
    !  field line points:

    r1 = (r0-r)/(rr-r)
    x = x - (x-xr)*r1
    y = y - (y-yr)*r1
    z = z - (z-zr)*r1
    GO TO 6010
6000 CONTINUE
    WRITE (*,'(/,/,a,a,/,/)') &
         '**** COMPUTATIONS IN THE SUBROUTINE TRACE ARE', &
         ' TERMINATED: THE CURRENT NUMBER OF POINTS EXCEEDED 1000 ****'
    l = 999
6010 CONTINUE
    xf = x
    yf = y
    zf = z
  END SUBROUTINE trace


  !=============================================================================


  SUBROUTINE shuetal_mgnp(xn_pd,vel,bzimf,xgsm,ygsm,zgsm,xmgnp,ymgnp, &
       zmgnp,dist,id)

    ! FOR ANY POINT OF SPACE WITH COORDINATES (XGSM,YGSM,ZGSM) AND SPECIFIED CONDITIONS
    ! IN THE INCOMING SOLAR WIND, THIS SUBROUTINE:

    ! (1) DETERMINES IF THE POINT (XGSM,YGSM,ZGSM) LIES INSIDE OR OUTSIDE THE
    !     MODEL MAGNETOPAUSE OF SHUE ET AL. (JGR-A, V.103, P. 17691, 1998).

    ! (2) CALCULATES THE GSM POSITION OF A POINT {XMGNP,YMGNP,ZMGNP}, LYING AT THE MODEL
    !     MAGNETOPAUSE AND ASYMPTOTICALLY TENDING TO THE NEAREST BOUNDARY POINT WITH
    !     RESPECT TO THE OBSERVATION POINT {XGSM,YGSM,ZGSM}, AS IT APPROACHES THE MAGNETO-
    !     PAUSE.

    ! INPUT: XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
    !                   OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
    !        BZIMF - IMF BZ IN NANOTESLAS

    !        VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
    !                 OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
    !                    FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY

    !        XGSM,YGSM,ZGSM - GSM POSITION OF THE OBSERVATION POINT IN EARTH RADII

    ! OUTPUT: XMGNP,YMGNP,ZMGNP - GSM POSITION OF THE BOUNDARY POINT
    !         DIST - DISTANCE (IN RE) BETWEEN THE OBSERVATION POINT (XGSM,YGSM,ZGSM)
    !                AND THE MODEL NAGNETOPAUSE
    !         ID -  POSITION FLAG:  ID=+1 (-1) MEANS THAT THE OBSERVATION POINT
    !         LIES INSIDE (OUTSIDE) OF THE MODEL MAGNETOPAUSE, RESPECTIVELY.

    ! OTHER SUBROUTINES USED: T96_MGNP

    !         AUTHOR:  N.A. TSYGANENKO,
    !         DATE:    APRIL 4, 2003.

    IMPLICIT NONE

    REAL (wp) :: xn_pd, vel, bzimf, xgsm, ygsm, zgsm, xmgnp, ymgnp, zmgnp, &
         dist

    REAL (wp) :: pd, phi, r0, alpha, r, rm, xmt96, ymt96, zmt96, rho2, st, &
         ct, t, f, gradf_r, gradf_t, gradf, dr, dt, ds, rho

    INTEGER :: id96, nit, id

    IF (vel<0._wp) THEN
      pd = xn_pd
    ELSE
      pd = 1.94E-6_wp*xn_pd*vel**2.0_wp ! PD IS THE SOLAR WIND DYNAMIC PRESSURE (IN nPa)
    END IF


    ! DEFINE THE ANGLE PHI, MEASURED DUSKWARD FROM THE NOON-MIDNIGHT MERIDIAN PLANE;
    ! IF THE OBSERVATION POINT LIES ON THE X AXIS, THE ANGLE PHI CANNOT BE UNIQUELY
    ! DEFINED, AND WE SET IT AT ZERO:

    IF (ygsm/=0._wp .OR. zgsm/=0._wp) THEN
      phi = atan2(ygsm,zgsm)
    ELSE
      phi = 0._wp
    END IF

    ! FIRST, FIND OUT IF THE OBSERVATION POINT LIES INSIDE THE SHUE ET AL BDRY
    ! AND SET THE VALUE OF THE ID FLAG:

    id = -1
    r0 = (10.22_wp+1.29_wp*tanh(0.184_wp*(bzimf+ &
         8.14_wp)))*pd**(-.15151515_wp)
    alpha = (0.58_wp-0.007_wp*bzimf)*(1._wp+0.024_wp*log(pd))
    r = sqrt(xgsm**2.0_wp+ygsm**2.0_wp+zgsm**2.0_wp)
    rm = r0*(2._wp/(1._wp+xgsm/r))**alpha
    IF (r<=rm) id = + 1

    ! NOW, FIND THE CORRESPONDING T96 MAGNETOPAUSE POSITION, TO BE USED AS
    ! A STARTING APPROXIMATION IN THE SEARCH OF A CORRESPONDING SHUE ET AL.
    ! BOUNDARY POINT:

    CALL t96_mgnp(pd,-1._wp,xgsm,ygsm,zgsm,xmt96,ymt96,zmt96,dist,id96)

    rho2 = ymt96**2.0_wp + zmt96**2.0_wp
    r = sqrt(rho2+xmt96**2.0_wp)
    st = sqrt(rho2)/r
    ct = xmt96/r

    ! NOW, USE NEWTON'S ITERATIVE METHOD TO FIND THE NEAREST POINT AT THE
    !  SHUE ET AL.'S BOUNDARY:

    nit = 0
    DO

      t = atan2(st,ct)
      rm = r0*(2._wp/(1._wp+ct))**alpha

      f = r - rm
      gradf_r = 1._wp
      gradf_t = -alpha/r*rm*st/(1._wp+ct)
      gradf = sqrt(gradf_r**2.0_wp+gradf_t**2.0_wp)

      dr = -f/gradf**2.0_wp
      dt = dr/r*gradf_t

      r = r + dr
      t = t + dt
      st = sin(t)
      ct = cos(t)

      ds = sqrt(dr**2.0_wp+(r*dt)**2.0_wp)

      nit = nit + 1

      IF (nit>1000) PRINT *, &
           ' BOUNDARY POINT COULD NOT BE FOUND; ITERATIONS DO NOT CONVERGE'

      IF (ds<=1.E-4_wp) EXIT
    END DO

    xmgnp = r*cos(t)
    rho = r*sin(t)

    ymgnp = rho*sin(phi)
    zmgnp = rho*cos(phi)

    dist = sqrt((xgsm-xmgnp)**2.0_wp+(ygsm-ymgnp)**2.0_wp+ &
         (zgsm-zmgnp)**2.0_wp)


  END SUBROUTINE shuetal_mgnp


  !=============================================================================


  SUBROUTINE t96_mgnp(xn_pd,vel,xgsm,ygsm,zgsm,xmgnp,ymgnp,zmgnp,dist,id)

    ! FOR ANY POINT OF SPACE WITH GIVEN COORDINATES (XGSM,YGSM,ZGSM), THIS SUBROUTINE DEFINES
    ! THE POSITION OF A POINT (XMGNP,YMGNP,ZMGNP) AT THE T96 MODEL MAGNETOPAUSE, HAVING THE
    ! SAME VALUE OF THE ELLIPSOIDAL TAU-COORDINATE, AND THE DISTANCE BETWEEN THEM.  THIS IS
    ! NOT THE SHORTEST DISTANCE D_MIN TO THE BOUNDARY, BUT DIST ASYMPTOTICALLY TENDS TO D_MIN,
    ! AS THE OBSERVATION POINT GETS CLOSER TO THE MAGNETOPAUSE.

    ! INPUT: XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
    !                   OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
    !        VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
    !                 OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
    !                    FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY

    !        XGSM,YGSM,ZGSM - COORDINATES OF THE OBSERVATION POINT IN EARTH RADII

    ! OUTPUT: XMGNP,YMGNP,ZMGNP - GSM POSITION OF THE BOUNDARY POINT, HAVING THE SAME
    !         VALUE OF TAU-COORDINATE AS THE OBSERVATION POINT (XGSM,YGSM,ZGSM)
    !         DIST -  THE DISTANCE BETWEEN THE TWO POINTS, IN RE,
    !         ID -    POSITION FLAG; ID=+1 (-1) MEANS THAT THE POINT (XGSM,YGSM,ZGSM)
    !         LIES INSIDE (OUTSIDE) THE MODEL MAGNETOPAUSE, RESPECTIVELY.

    ! THE PRESSURE-DEPENDENT MAGNETOPAUSE IS THAT USED IN THE T96_01 MODEL
    ! (TSYGANENKO, JGR, V.100, P.5599, 1995; ESA SP-389, P.181, OCT. 1996)

    !  AUTHOR:  N.A. TSYGANENKO
    !  DATE:    AUG.1, 1995, REVISED APRIL 3, 2003.


    ! DEFINE SOLAR WIND DYNAMIC PRESSURE (NANOPASCALS, ASSUMING 4% OF ALPHA-PARTICLES),
    !  IF NOT EXPLICITLY SPECIFIED IN THE INPUT:

    IMPLICIT NONE

    REAL (wp) :: xn_pd, vel, xgsm, ygsm, zgsm, xmgnp, ymgnp, zmgnp, dist

    INTEGER :: id

    REAL (wp) :: pd, rat, rat16, a0, s00, x00, a, s0, x0, xm, phi, rho, &
         rhomgnp, xksi, xdzt, sq1, sq2, sigma, tau, arg

    IF (vel<0._wp) THEN
      pd = xn_pd
    ELSE
      pd = 1.94E-6_wp*xn_pd*vel**2.0_wp

    END IF

    ! RATIO OF PD TO THE AVERAGE PRESSURE, ASSUMED EQUAL TO 2 nPa:

    rat = pd/2.0_wp
    rat16 = rat**0.14_wp

    ! (THE POWER INDEX 0.14 IN THE SCALING FACTOR IS THE BEST-FIT VALUE OBTAINED FROM DATA
    !   AND USED IN THE T96_01 VERSION)

    ! VALUES OF THE MAGNETOPAUSE PARAMETERS FOR  PD = 2 nPa:

    a0 = 70._wp
    s00 = 1.08_wp
    x00 = 5.48_wp

    !  VALUES OF THE MAGNETOPAUSE PARAMETERS, SCALED BY THE ACTUAL PRESSURE:

    a = a0/rat16
    s0 = s00
    x0 = x00/rat16
    xm = x0 - a

    ! (XM IS THE X-COORDINATE OF THE "SEAM" BETWEEN THE ELLIPSOID AND THE CYLINDER)

    !    (FOR DETAILS ON THE ELLIPSOIDAL COORDINATES, SEE THE PAPER:
    !     N.A.TSYGANENKO, SOLUTION OF CHAPMAN-FERRARO PROBLEM FOR AN
    !     ELLIPSOIDAL MAGNETOPAUSE, PLANET.SPACE SCI., V.37, P.1037, 1989).

    IF (ygsm/=0._wp .OR. zgsm/=0._wp) THEN
      phi = atan2(ygsm,zgsm)
    ELSE
      phi = 0._wp
    END IF

    rho = sqrt(ygsm**2.0_wp+zgsm**2.0_wp)

    IF (xgsm<xm) THEN
      xmgnp = xgsm
      rhomgnp = a*sqrt(s0**2.0_wp-1)
      ymgnp = rhomgnp*sin(phi)
      zmgnp = rhomgnp*cos(phi)
      dist = sqrt((xgsm-xmgnp)**2.0_wp+(ygsm-ymgnp)**2.0_wp+ &
           (zgsm-zmgnp)**2.0_wp)
      IF (rhomgnp>rho) id = + 1
      IF (rhomgnp<=rho) id = -1
    ELSE

      xksi = (xgsm-x0)/a + 1._wp
      xdzt = rho/a
      sq1 = sqrt((1._wp+xksi)**2.0_wp+xdzt**2.0_wp)
      sq2 = sqrt((1._wp-xksi)**2.0_wp+xdzt**2.0_wp)
      sigma = 0.5_wp*(sq1+sq2)
      tau = 0.5_wp*(sq1-sq2)

      ! NOW CALCULATE (X,Y,Z) FOR THE CLOSEST POINT AT THE MAGNETOPAUSE

      xmgnp = x0 - a*(1._wp-s0*tau)
      arg = (s0**2.0_wp-1._wp)*(1._wp-tau**2.0_wp)
      IF (arg<0._wp) arg = 0._wp
      rhomgnp = a*sqrt(arg)
      ymgnp = rhomgnp*sin(phi)
      zmgnp = rhomgnp*cos(phi)

      ! NOW CALCULATE THE DISTANCE BETWEEN THE POINTS {XGSM,YGSM,ZGSM} AND {XMGNP,YMGNP,ZMGNP}:
      !  (IN GENERAL, THIS IS NOT THE SHORTEST DISTANCE D_MIN, BUT DIST ASYMPTOTICALLY TENDS
      !   TO D_MIN, AS WE ARE GETTING CLOSER TO THE MAGNETOPAUSE):

      dist = sqrt((xgsm-xmgnp)**2.0_wp+(ygsm-ymgnp)**2.0_wp+ &
           (zgsm-zmgnp)**2.0_wp)

      IF (sigma>s0) THEN !  ID=-1 MEANS THAT THE POINT LIES OUTSIDE

        id = -1 !  ID=-1 MEANS THAT THE POINT LIES OUTSIDE
      END IF
      IF (sigma<=s0) THEN !  ID=+1 MEANS THAT THE POINT LIES INSIDE
        id = + 1 !  ID=+1 MEANS THAT THE POINT LIES INSIDE
      END IF
      !                                          THE MAGNETOSPHERE

    END IF

  END SUBROUTINE t96_mgnp

END MODULE mo_geopack
