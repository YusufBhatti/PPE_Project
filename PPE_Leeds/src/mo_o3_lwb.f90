!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_o3_lwb

!- Description:
!
!  This module contains a global ozone climatology following London et al. (1976)
!  for ozone amount and Wilcox and Belmont (1977) for altitude of the maximum 
!  concentration. Both fields are given by T5 spectral representations.
!  The spectral fields are given for 4 seasons. These spectral data are inter-
!  polated to the actual model time, done by
!  subroutine pre_o3_lwb called from prerad. The conversion to grid space
!  is done by the subroutine o3_lwb called from radint.
!
!  This module contains:
!
!  A) spectral component fields
!  B) the subroutine pre_o3_lwb to provide time dependent spectral components
!  C) the subroutine o3_lwb to provide time gridded fields
!
!
!  This module replaces the following parts of the ECHAM4.f90 model:
!   - part of mo_rad1 for ozone
!   - subroutine ozone
!   - parts of radint for ozone field computation
!
!- Author:
!
!  Marco Giorgetta, MPI, April 1999
!  U. Schulzweida,  MPI, May 2002, blocking (nproma)

!=======================================================================

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pre_o3_lwb, o3_lwb

!=======================================================================
!
! A)
  REAL(dp), dimension(21) :: cozqc ! ozone amount cosine coefficients
  REAL(dp), dimension(15) :: cozqs ! ozone amount sine coefficients
  REAL(dp), dimension(21) :: cozhc ! altitude of ozone max. cosine coefficients
  REAL(dp), dimension(15) :: cozhs ! altitude of ozone max. sine coefficients

!=======================================================================

  CONTAINS

!=======================================================================
!
! B)

SUBROUTINE pre_o3_lwb(pytime) ! former ozone,
                              ! provides time dependent spectral components


  ! Description:
  !
  ! Computes instantaneous values for the yearly cycle of the 
  ! ozone distibution.
  !
  ! Method:
  !
  ! This routine computes instantaneous values of a t5 spetral index
  ! distribution for two ozone parameters (total quantity and pressure
  ! at the maximum of concentration,both in *pascal) from the time
  ! of the year (see *orbit*).
  !
  ! Staightforward, a second order *Fourier development for the
  ! time of the year.
  !
  ! *ozone*   is called from *physc* at the first latitude row at
  !           the time of a full radiation computation.
  ! There are five dummy arguments: 
  !         *pytime*   is the time of the year (in radians).
  !         *cozqc*, *pozqs*, *pozhc* and *pozhs* are arrays for the 
  !                    t5 distributions (*q for quantity, 
  !                    *h for height, *c for cosine and *s for sine).
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, June 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! M.A. Giorgetta, MPI, Aril 1999, modification for new module mo_o3_london
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_gaussgrid,      ONLY: coslon, sinlon, gl_twomu
  USE mo_geoloc,         ONLY: zozq_x, zozh_x
  USE mo_decomposition,  ONLY: ldc => local_decomposition
  USE mo_transpose,      ONLY: reorder

  !  Scalar arguments 
  REAL(dp) :: pytime

  !  Module array arguments 
  !REAL(dp) :: cozhc(21), cozhs(15), cozqc(21), cozqs(15)

  !  Local scalars: 
  REAL(dp)    :: zc1yt, zc2yt, zs1yt, zs2yt, zytime
  REAL(dp)    :: zzozh
  REAL(dp)    :: zcos1, zcos2, zcos3, zcos4, zcos5,  &
             zsin1, zsin2, zsin3, zsin4, zsin5
  REAL(dp)    :: zsin
  INTEGER :: jmn
  INTEGER :: jmm,imm,imnc,imns,jnn,jl

  INTEGER :: nglon, nglat, nlof
  INTEGER :: jlat, jglat

  !  Local arrays: 
  REAL(dp) :: zozhc0(21), zozhc1(21), zozhc2(21), zozhc3(21), zozhc4(21), &
          zozhs0(15), zozhs1(15), zozhs2(15), zozhs3(15), zozhs4(15), zozqc0(21), &
          zozqc1(21), zozqc2(21), zozqc3(21), zozqc4(21), zozqs0(15), zozqs1(15), &
          zozqs2(15), zozqs3(15), zozqs4(15)
  REAL(dp),DIMENSION(11)           :: zfozh,zfozq
  REAL(dp),DIMENSION(66)           :: zalp
  REAL(dp)    :: zozq(ldc% nglon, ldc% nglat)
  REAL(dp)    :: zozh(ldc% nglon, ldc% nglat)

  !  Intrinsic functions 
  INTRINSIC COS, SIN

  !  External subroutines
  EXTERNAL  legtri  ! for the Legendre transform at a given latitude

  !  Data statements 
  ! *zoz q/h c/s n* (n=0,4) corresponds to the *poz q/h c/s*
  ! (see above) and to the five terms of the *fourier development.
  DATA zozqc0/  0.6012E-01_dp,  0.1887E-02_dp,  0.7410E-02_dp,  0.9950E-03_dp, &
               -0.1426E-02_dp, -0.2072E-03_dp, -0.4954E-03_dp,  0.7955E-05_dp, -0.3701E-03_dp, &
                0.4116E-04_dp, -0.4163E-04_dp, -0.2933E-03_dp,  0.2154E-04_dp, -0.2849E-03_dp, &
               -0.1604E-03_dp, -0.1054E-03_dp,  0.4974E-03_dp,  0.1047E-03_dp,  0.8323E-04_dp, &
                0.2874E-03_dp,  0.1333E-03_dp/
  DATA zozqs0/  0.4210E-03_dp, -0.9591E-03_dp,  0.2811E-03_dp, -0.2257E-03_dp, -0.1713E-03_dp, &
               -0.3538E-03_dp,  0.1095E-03_dp, -0.4390E-03_dp, -0.5605E-05_dp,  0.1478E-03_dp, &
                0.2849E-03_dp,  0.3430E-03_dp,  0.8248E-04_dp,  0.1442E-03_dp, -0.1375E-04_dp/
  DATA zozhc0/  0.3166E+04_dp,  0.8663E+02_dp,  0.9401E+03_dp,  0.1999E+02_dp, &
               -0.3530E+03_dp, -0.3311E+02_dp, -0.4903E+02_dp, -0.4015E+00_dp, -0.1333E+02_dp, &
                0.5675E+01_dp,  0.7221E+01_dp, -0.3001E+02_dp,  0.7570E+01_dp, -0.1142E+02_dp, &
               -0.1365E+02_dp, -0.1502E+02_dp,  0.4911E+02_dp,  0.1425E+02_dp,  0.8983E+01_dp, &
                0.3064E+02_dp,  0.1693E+02_dp/
  DATA zozhs0/  0.4231E+02_dp, -0.7391E+02_dp,  0.1273E+02_dp,  0.2086E+02_dp, -0.1597E+02_dp, &
               -0.3591E+02_dp,  0.1059E+02_dp, -0.2779E+02_dp, -0.6923E+01_dp,  0.1397E+02_dp, &
                0.2387E+02_dp,  0.2883E+02_dp,  0.8626E+01_dp,  0.1607E+02_dp, -0.2676E+01_dp/
  DATA zozqc1/  0.7090E-04_dp,  0.4930E-05_dp,  0.6829E-03_dp,  0.1897E-03_dp, &
                0.7226E-04_dp, -0.2807E-03_dp,  0.4970E-04_dp, -0.1753E-03_dp, -0.7843E-04_dp, &
               -0.1649E-03_dp, -0.1037E-03_dp, -0.4830E-04_dp, -0.6304E-04_dp, -0.1100E-03_dp, -0.7952E-04_dp, &
                0.1326E-04_dp,  0.2599E-04_dp,  0.9926E-05_dp, -0.9247E-05_dp, -0.3521E-05_dp, &
               -0.1780E-04_dp/
  DATA zozqs1/  0.6333E-04_dp,  0.1145E-03_dp,  0.1192E-03_dp,  0.4934E-04_dp, &
                0.2699E-04_dp,  0.3684E-04_dp, -0.2395E-05_dp,  0.2045E-04_dp, -0.8684E-04_dp, &
                0.5301E-04_dp, -0.4176E-05_dp,  0.4103E-04_dp,  0.2783E-04_dp,  0.1754E-04_dp, &
                0.1116E-04_dp/
  DATA zozhc1/ -0.3450E+02_dp,  0.2148E+03_dp,  0.3376E+02_dp,  0.6535E+02_dp, -0.1564E+02_dp, &
               -0.4273E+02_dp,  0.9553E+01_dp, -0.4647E+01_dp, -0.6129E+01_dp, -0.6727E+01_dp, &
               -0.6761E+01_dp, -0.2467E+01_dp, -0.2181E+01_dp, -0.5361E+01_dp, -0.2395E+01_dp, &
                0.5952E+00_dp,  0.2106E+01_dp, -0.1367E+01_dp, -0.2349E+01_dp,  0.3532E+00_dp, &
               -0.3169E+01_dp/
  DATA zozhs1/  0.3977E+01_dp,  0.5032E+01_dp,  0.6226E+01_dp, -0.3625E+00_dp, -0.1373E+01_dp, &
                0.4600E+01_dp,  0.4312E+01_dp,  0.2882E+01_dp, -0.6351E+01_dp,  0.5731E+01_dp, &
               -0.2574E+01_dp,  0.3235E+00_dp,  0.2806E+01_dp,  0.8133E+00_dp,  0.2032E+01_dp/
  DATA zozqc2/  0.8571E-03_dp,  0.3086E-02_dp,  0.9287E-03_dp,  0.2787E-03_dp, &
                0.1826E-03_dp, -0.1006E-03_dp,  0.1092E-03_dp, -0.1266E-03_dp,  0.5372E-04_dp, &
               -0.1188E-03_dp, -0.3285E-04_dp, -0.1783E-04_dp, -0.3018E-05_dp, -0.8709E-04_dp, -0.8707E-04_dp, &
                0.8633E-04_dp,  0.3530E-04_dp,  0.4863E-04_dp,  0.3917E-05_dp, -0.3252E-04_dp, &
               -0.1936E-06_dp/
  DATA zozqs2/ -0.8822E-04_dp,  0.1341E-03_dp,  0.3095E-04_dp,  0.8230E-04_dp, &
                0.2735E-04_dp,  0.1714E-04_dp, -0.9406E-04_dp,  0.1912E-04_dp, -0.5402E-04_dp, &
                0.3571E-04_dp,  0.3897E-04_dp,  0.4487E-04_dp,  0.3079E-04_dp,  0.3196E-04_dp, &
               -0.2391E-05_dp/
  DATA zozhc2/  0.5216E+02_dp,  0.1613E+03_dp,  0.3284E+02_dp, -0.7670E+02_dp, -0.9548E+01_dp, &
                0.1608E+02_dp,  0.1023E+02_dp, -0.1090E+02_dp,  0.2748E+01_dp, -0.3846E+01_dp, &
               -0.4135E+01_dp,  0.1255E+01_dp, -0.3301E-01_dp, -0.5273E+01_dp, -0.7247E+01_dp, &
                0.1387E+02_dp,  0.4184E+01_dp,  0.6495E+01_dp,  0.2944E+01_dp, -0.1947E+01_dp, &
                0.1132E+01_dp/
  DATA zozhs2/ -0.1968E+02_dp,  0.1192E+02_dp, -0.1194E+01_dp,  0.1084E+01_dp,  0.2946E+01_dp, &
                0.2630E+01_dp, -0.1256E+02_dp,  0.1395E+01_dp, -0.2222E+01_dp,  0.4864E+01_dp, &
                0.6450E+01_dp,  0.5568E+01_dp,  0.5292E+01_dp,  0.4876E+01_dp, -0.7579E+00_dp/
  DATA zozqc3/ -0.2759E-03_dp, -0.2781E-03_dp, -0.1087E-03_dp, -0.1633E-03_dp, -0.3627E-04_dp, &
               -0.4242E-04_dp,  0.6045E-05_dp, -0.1703E-04_dp,  0.4562E-04_dp, -0.1009E-04_dp, &
                0.2663E-04_dp, -0.1786E-04_dp,  0.1550E-04_dp, -0.9135E-06_dp,  0.2372E-04_dp, &
                0.1100E-05_dp,  0.2299E-04_dp,  0.4659E-05_dp,  0.2423E-05_dp,  0.7321E-05_dp, &
                0.8852E-05_dp/
  DATA zozqs3/ -0.3678E-04_dp, -0.2219E-04_dp, -0.3911E-04_dp, -0.4398E-04_dp, -0.1142E-04_dp, &
               -0.9121E-05_dp, -0.2011E-04_dp,  0.4711E-06_dp, -0.3775E-05_dp,  0.3866E-05_dp, &
                0.2400E-04_dp,  0.2043E-04_dp, -0.1824E-05_dp, -0.5550E-05_dp,  0.2506E-05_dp/
  DATA zozhc3/ -0.1534E+03_dp, -0.2095E+02_dp, -0.1006E+03_dp, -0.7385E+01_dp,  0.5203E+01_dp, &
                0.9434E+00_dp, -0.3814E+00_dp, -0.3175E+01_dp,  0.3366E+01_dp,  0.3378E+00_dp, &
                0.2740E+00_dp, -0.2669E+01_dp,  0.8452E+00_dp,  0.3498E+00_dp,  0.2192E+01_dp, &
               -0.4024E+00_dp,  0.1544E+01_dp, -0.4588E+00_dp,  0.6998E+00_dp,  0.6263E+00_dp, &
                0.1228E+01_dp/
  DATA zozhs3/ -0.3588E+01_dp,  0.2076E+00_dp, -0.2088E+01_dp, -0.4159E+01_dp,  0.2244E+00_dp, &
               -0.7751E+00_dp, -0.2749E+01_dp,  0.7234E+00_dp,  0.4390E+00_dp, -0.1646E+00_dp, &
                0.1700E+01_dp,  0.1046E+01_dp, -0.7856E+00_dp, -0.1644E+01_dp,  0.2648E+00_dp/
  DATA zozqc4/ -0.1460E-03_dp,  0.3422E-03_dp, -0.3529E-04_dp,  0.1791E-03_dp, -0.1917E-03_dp, &
               -0.2558E-04_dp,  0.6547E-04_dp,  0.6401E-04_dp,  0.4823E-04_dp,  0.7084E-05_dp, &
                0.2895E-04_dp, -0.1561E-04_dp,  0.8179E-06_dp,  0.1028E-04_dp, -0.7667E-05_dp, &
               -0.4347E-05_dp,  0.7293E-05_dp, -0.5735E-05_dp,  0.7838E-05_dp, -0.2933E-05_dp, &
                0.3686E-05_dp/
  DATA zozqs4/ -0.4560E-05_dp, -0.5292E-04_dp, -0.1252E-04_dp,  0.1850E-04_dp, -0.2273E-04_dp, &
                0.6552E-05_dp,  0.1422E-04_dp, -0.6545E-05_dp,  0.7998E-06_dp,  0.2845E-04_dp, &
                0.2497E-04_dp,  0.2844E-04_dp,  0.3855E-06_dp, -0.1487E-04_dp,  0.1954E-05_dp/
  DATA zozhc4/  0.9260E+01_dp, -0.9055E+01_dp,  0.5460E+01_dp, -0.7603E+01_dp, -0.3329E+02_dp, &
               -0.1048E+02_dp,  0.9328E+01_dp,  0.4597E+01_dp,  0.3827E+01_dp, -0.3201E+01_dp, &
                0.1708E+01_dp, -0.1548E+01_dp, -0.5323E+00_dp,  0.3039E+01_dp,  0.5740E+00_dp, &
                0.1353E+00_dp, -0.2354E+01_dp,  0.2818E+00_dp,  0.1113E+01_dp, -0.1891E+01_dp, &
               -0.3074E+00_dp/
  DATA zozhs4/ -0.2446E+01_dp,  0.4199E+01_dp, -0.2571E+01_dp,  0.8194E+01_dp,  0.4206E+00_dp, &
                0.3856E+01_dp,  0.1159E+01_dp,  0.2547E+01_dp, -0.1314E+01_dp,  0.2331E+01_dp, &
                0.1144E+01_dp, -0.4408E+00_dp, -0.6797E+00_dp, -0.2598E+01_dp,  0.8953E+00_dp/


  !  Executable statements 

!-- 1. Preliminary setting

  zytime = pytime

!-- 2. Computations

  zc1yt = COS(zytime)
  zs1yt = SIN(zytime)
  zc2yt = zc1yt**2 - zs1yt**2
  zs2yt = 2._dp*zs1yt*zc1yt
  DO jmn = 1, 21
    cozqc(jmn) = zozqc0(jmn) + 2._dp*(zozqc1(jmn)*zc1yt+zozqc2(jmn)*zs1yt+zozqc3 &
&        (jmn)*zc2yt+zozqc4(jmn)*zs2yt)
    cozhc(jmn) = zozhc0(jmn) + 2._dp*(zozhc1(jmn)*zc1yt+zozhc2(jmn)*zs1yt+zozhc3 &
&        (jmn)*zc2yt+zozhc4(jmn)*zs2yt)
  END DO
  DO jmn = 1, 15
    cozqs(jmn) = zozqs0(jmn) + 2._dp*(zozqs1(jmn)*zc1yt+zozqs2(jmn)*zs1yt+zozqs3 &
&        (jmn)*zc2yt+zozqs4(jmn)*zs2yt)
    cozhs(jmn) = zozhs0(jmn) + 2._dp*(zozhs1(jmn)*zc1yt+zozhs2(jmn)*zs1yt+zozhs3 &
&        (jmn)*zc2yt+zozhs4(jmn)*zs2yt)

  END DO


  !  Local array bounds

  nglon = ldc% nglon      ! local number of longitudes
  nglat = ldc% nglat      ! local number of latitudes

  DO jlat = 1, nglat

    jglat = ldc% glat(jlat)
    nlof  = ldc% glon(jlat) ! longitude offset to global field
     
    zsin = .5_dp*gl_twomu(jglat)

    ! Prepare Legendre transform at given latitute

    CALL legtri(zsin,6,zalp)

    ! Legendre transform

    DO jmm = 1, 11
      zfozq(jmm) = 0._dp
      zfozh(jmm) = 0._dp
    END DO
    imm = 0
    imnc = 0
    imns = 0
    DO jmm = 1, 6
      imm = imm + 1
      DO jnn = jmm, 6
        imnc = imnc + 1
        zfozq(imm) = zfozq(imm) + zalp(imnc)*cozqc(imnc)
        zfozh(imm) = zfozh(imm) + zalp(imnc)*cozhc(imnc)
      END DO
      IF (jmm/=1) THEN
        imm = imm + 1
        DO jnn = jmm, 6
          imns = imns + 1
          zfozq(imm) = zfozq(imm) + zalp(imns+6)*cozqs(imns)
          zfozh(imm) = zfozh(imm) + zalp(imns+6)*cozhs(imns)
        END DO
      END IF
    END DO

    ! Fourier transform

    DO jl = 1, nglon

      zcos1 = coslon(jl+nlof)
      zsin1 = sinlon(jl+nlof)
      zcos2 = zcos1*zcos1 - zsin1*zsin1
      zsin2 = zsin1*zcos1 + zcos1*zsin1
      zcos3 = zcos2*zcos1 - zsin2*zsin1
      zsin3 = zsin2*zcos1 + zcos2*zsin1
      zcos4 = zcos3*zcos1 - zsin3*zsin1
      zsin4 = zsin3*zcos1 + zcos3*zsin1
      zcos5 = zcos4*zcos1 - zsin4*zsin1
      zsin5 = zsin4*zcos1 + zcos4*zsin1

      zozq(jl,jlat) = zfozq(1) + 2._dp*(zfozq(2)*zcos1+zfozq(3)*zsin1+zfozq(4)*zcos2+ &
                      zfozq(5)*zsin2+zfozq(6)*zcos3+zfozq(7)*zsin3+zfozq(8)*zcos4+ &
                      zfozq(9)*zsin4+zfozq(10)*zcos5+zfozq(11)*zsin5)
      zzozh = zfozh(1) + 2._dp*(zfozh(2)*zcos1+zfozh(3)*zsin1+zfozh(4)*zcos2+ &
              zfozh(5)*zsin2+zfozh(6)*zcos3+zfozh(7)*zsin3+zfozh(8)*zcos4+ &
              zfozh(9)*zsin4+zfozh(10)*zcos5+zfozh(11)*zsin5)
      zozh(jl,jlat) = SQRT(zzozh)**3

    END DO

  END DO

  CALL reorder (zozq_x, zozq)
  CALL reorder (zozh_x, zozh)

END SUBROUTINE pre_o3_lwb

!=======================================================================
!
! C)

FUNCTION o3_lwb(jrow,pdp,pph)

  ! follows ECHAM4 radint

  USE mo_geoloc,         ONLY: zozq_x, zozh_x

  INTEGER ,INTENT(in) :: jrow ! local latitude index
  REAL(dp), INTENT(in) , DIMENSION(:,:)        :: pdp,pph
  REAL(dp), DIMENSION(SIZE(pdp,1),SIZE(pdp,2)) :: o3_lwb

  INTEGER             :: klon,klev
  INTEGER             :: jl,jk

  REAL(dp)                :: zcphn3,zcpho3,zsdpn3,zsdpo3,zqofn,zqofo


  klon = SIZE(pdp,1)
  klev = SIZE(pdp,2)

  ! output: ozone longitude height section in o3_lwb

  DO jk = 1, klev
    DO jl = 1, klon
      zcpho3 = pph(jl,jk)**3
      zcphn3 = pph(jl,jk+1)**3
      zsdpo3 = SQRT(zcpho3)
      zsdpn3 = SQRT(zcphn3)
      zqofo = zozq_x(jl,jrow)*zsdpo3/(zsdpo3+zozh_x(jl,jrow))
      zqofn = zozq_x(jl,jrow)*zsdpn3/(zsdpn3+zozh_x(jl,jrow))
      o3_lwb(jl,jk) = zqofn - zqofo
    END DO
  END DO

    o3_lwb=o3_lwb/pdp

END FUNCTION o3_lwb

!=======================================================================

END MODULE mo_o3_lwb
