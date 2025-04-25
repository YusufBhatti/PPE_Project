!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_srtm_setup

  USE mo_kind, ONLY : wp
  USE mo_rrtm_params, ONLY : nbndsw, mg, ngptsw, jpb1, jpb2
  USE mo_srtm_netcdf, ONLY : srtm_read
  USE mo_rad_fastmath,ONLY : setup_rad_fastmath

  IMPLICIT NONE

  ! ------- Definitions -------
  !     Arrays for the g-point reduction from 224 to 112 for the 16 LW bands:
  !     This mapping from 224 to 112 points has been carefully selected to 
  !     minimize the effect on the resulting fluxes and cooling rates, and
  !     caution should be used if the mapping is modified.  The full 224
  !     g-point set can be restored with ngpt=224, ngc=16*16, ngn=224*1., etc.
  !     ngpt    The total number of new g-points
  !     ngc     The number of new g-points in each band
  !     ngs     The cumulative sum of new g-points for each band
  !     ngm     The index of each new g-point relative to the original
  !             16 g-points for each band.  
  !     ngn     The number of original g-points that are combined to make
  !             each new g-point in each band.
  !     ngb     The band index for each new g-point.
  !     wt      

  INTEGER, PARAMETER :: ng(jpb1:jpb2)   = (/&
       & 16,16,16,16,16,16,16,16,16,16,16,16, 16, 16/)
  INTEGER, PARAMETER :: nspa(jpb1:jpb2) = (/&
       &  9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1,  9,  1/)
  INTEGER, PARAMETER :: nspb(jpb1:jpb2) = (/&
       &  1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1,  5,  1/)
  INTEGER, PARAMETER :: ngc(nbndsw)=      (/&
       &  6,12, 8, 8,10,10, 2,10, 8, 6, 6, 8,  6, 12/)
  INTEGER, PARAMETER :: ngs(nbndsw) =     (/&
       &  6,18,26,34,44,54,56,66,74,80,86,94,100,112/)

  INTEGER, PARAMETER :: ngm(nbndsw*mg) = (/ &
       & 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, & ! band 16
       & 1, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9,10,10,11,12,12, & ! band 17
       & 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, & ! band 18
       & 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, & ! band 19
       & 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,10,10,10,10,10, & ! band 20
       & 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,10,10,10,10,10, & ! band 21
       & 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, & ! band 22
       & 1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,10,10,10, & ! band 23
       & 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, & ! band 24
       & 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! band 25
       & 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! band 26
       & 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, & ! band 27
       & 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! band 28
       & 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9,10,11,12 /) ! band 29

  INTEGER, PARAMETER :: ngn(ngptsw) = (/ &
       2,2,2,2,4,4,             & ! band 16
       1,1,1,1,1,2,1,2,1,2,1,2, & ! band 17
       1,1,1,1,2,2,4,4,         & ! band 18
       1,1,1,1,2,2,4,4,         & ! band 19
       1,1,1,1,1,1,1,1,2,6,     & ! band 20
       1,1,1,1,1,1,1,1,2,6,     & ! band 21
       8,8,                     & ! band 22
       2,2,1,1,1,1,1,1,2,4,     & ! band 23
       2,2,2,2,2,2,2,2,         & ! band 24
       1,1,2,2,4,6,             & ! band 25
       1,1,2,2,4,6,             & ! band 26
       1,1,1,1,1,1,4,6,         & ! band 27
       1,1,2,2,4,6,             & ! band 28
       1,1,1,1,2,2,2,2,1,1,1,1 /) ! band 29

  INTEGER, PARAMETER :: ngb(ngptsw) = (/    &
       16,16,16,16,16,16,                   & ! band 16
       17,17,17,17,17,17,17,17,17,17,17,17, & ! band 17
       18,18,18,18,18,18,18,18,             & ! band 18
       19,19,19,19,19,19,19,19,             & ! band 19
       20,20,20,20,20,20,20,20,20,20,       & ! band 20
       21,21,21,21,21,21,21,21,21,21,       & ! band 21
       22,22,                               & ! band 22
       23,23,23,23,23,23,23,23,23,23,       & ! band 23
       24,24,24,24,24,24,24,24,             & ! band 24
       25,25,25,25,25,25,                   & ! band 25
       26,26,26,26,26,26,                   & ! band 26
       27,27,27,27,27,27,27,27,             & ! band 27
       28,28,28,28,28,28,                   & ! band 28
       29,29,29,29,29,29,29,29,29,29,29,29 /) ! band 29

  ! Shortwave spectral band limits (wavenumbers)
  REAL(wp), PARAMETER :: wavenum1(jpb1:jpb2) = (/ &
       2600._wp, 3250._wp, 4000._wp, 4650._wp, 5150._wp, 6150._wp, 7700._wp, &
       8050._wp,12850._wp,16000._wp,22650._wp,29000._wp,38000._wp,  820._wp/)
  REAL(wp), PARAMETER :: wavenum2(jpb1:jpb2) = (/ &
       3250._wp, 4000._wp, 4650._wp, 5150._wp, 6150._wp, 7700._wp, 8050._wp, &
       12850._wp,16000._wp,22650._wp,29000._wp,38000._wp,50000._wp, 2600._wp/)
  REAL(wp), PARAMETER :: delwave(jpb1:jpb2)  = (/ &
       650._wp,  750._wp,  650._wp,  500._wp, 1000._wp, 1550._wp,  350._wp, &
       4800._wp, 3150._wp, 6650._wp, 6350._wp, 9000._wp,12000._wp, 1780._wp/)

  REAL(wp) :: wt(mg) = (/ & !< RRTM weights for 16 g-points.
       0.1527534276_wp, 0.1491729617_wp, 0.1420961469_wp, &
       0.1316886544_wp, 0.1181945205_wp, 0.1019300893_wp, &
       0.0832767040_wp, 0.0626720116_wp, 0.0424925000_wp, &
       0.0046269894_wp, 0.0038279891_wp, 0.0030260086_wp, &
       0.0022199750_wp, 0.0014140010_wp, 0.0005330000_wp, &
       0.0000750000_wp /)

  REAL(wp) :: rwgt(nbndsw*mg)

  REAL(wp), PARAMETER :: ssi_default(14) =  (/ & !< SRTM default solar flux (W/m2) in 14 SW bands
    & 12.1095682699999987_wp   , 20.3650825429849398_wp   , 23.7297328613475429_wp,    &
    & 22.4276934179066849_wp   , 55.6266126199999960_wp   , 1.0293153082385436E+02_wp, &
    & 24.2936128100154392_wp   , 3.4574251380000004E+02_wp, 2.1818712729860400E+02_wp, &
    & 3.4719231470000005E+02_wp, 1.2949501812000000E+02_wp, 50.1522503011060650_wp,    &
    & 3.0799387010047838_wp    , 12.8893773299999985_wp  /)
    ! sum of 14 bands is: 1.3682223735968237E+03

  REAL(wp), PARAMETER :: ssi_amip(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for
                           ! AMIP-type CMIP5 simulation (average from 1979-1988)
    & 11.95053_wp, 20.14766_wp, 23.40394_wp, 22.09458_wp, 55.41401_wp,  &
    & 102.5134_wp, 24.69814_wp, 347.5362_wp, 217.2925_wp, 343.4221_wp,  &
    & 129.403_wp, 47.14264_wp, 3.172126_wp, 13.18075_wp /)
    ! sum of 14 bands is: 1361.371

  REAL(wp), PARAMETER :: ssi_preind(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for
                           ! preindutrial CMIP5 simulation (average from 1944-1856)
    & 11.95005_wp, 20.14612_wp, 23.40302_wp, 22.09443_wp, 55.41679_wp,  &
    & 102.512_wp , 24.69536_wp, 347.4719_wp, 217.2217_wp, 343.2816_wp,  &
    & 129.3001_wp, 47.07624_wp, 3.130199_wp, 13.17521_wp /)
    ! sum of 14 bands is: 1360.875

  ! ssi_RCEdiurnOn added (diurnal cycle on)
  REAL(wp), PARAMETER :: ssi_RCEdiurnOn(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for
                           ! RCE simulations with diurnal cycle. global mean insolation = 340.3 W/m2
    & 9.386766_wp,  15.82535_wp,  18.38306_wp,   17.3546_wp,  43.52597_wp,  &
    & 80.52106_wp,  19.39961_wp,  272.9788_wp,  170.6764_wp,  269.7473_wp,  &
    & 101.642_wp,   37.02906_wp,  2.491606_wp,  10.35307_wp/)
    ! sum of 14 bands is: 1069.315 for diurnal cycle on.

  ! ssi_RCEdiurnOFF added (diurnal cycle off)
  REAL(wp), PARAMETER :: ssi_RCEdiurnOFF(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for
                             ! RCE simulations with diurnal cycle switched off. global mean insolation = 340.3 W/m2
                             ! rescaled from ssi_amip above, with constant factor of app. 0.3183092
    & 3.803964_wp,  6.413186_wp,  7.44969_wp,   7.032908_wp,  17.63879_wp, &
    & 32.63096_wp,  7.861645_wp,  110.624_wp,   69.16621_wp,  109.3144_wp, &
    & 41.19017_wp,  15.00594_wp, 1.009717_wp,   4.195554_wp  /)
    ! sum of 14 bands is: 433.3371

CONTAINS

  ! **************************************************************************
  SUBROUTINE setup_srtm
    ! **************************************************************************
    !
    !  Original version:   Michael J. Iacono; February, 2004
    !  Revision for F90 formatting:  M. J. Iacono, July, 2006
    !
    !  This subroutine performs calculations necessary for the initialization
    !  of the shortwave model.  Lookup tables are computed for use in the SW
    !  radiative transfer, and input absorption coefficient data for each
    !  spectral band are reduced from 224 g-point intervals to 112.
    ! **************************************************************************

    ! ------- Local -------

    INTEGER :: ibnd, igc, ig, ind, ipr
    INTEGER :: igcsm, iprsm

    REAL(wp) :: wtsum, wtsm(mg)

    ! Initialize model data

    CALL srtm_read
    CALL setup_rad_fastmath

    ! Perform g-point reduction from 16 per band (224 total points) to
    ! a band dependent number (112 total points) for all absorption
    ! coefficient input data and Planck fraction input data.
    ! Compute relative weighting for new g-point combinations.

    igcsm = 0
    DO ibnd = 1,nbndsw
      iprsm = 0
      IF (ngc(ibnd).LT.mg) THEN
        DO igc = 1,ngc(ibnd)
          igcsm = igcsm + 1
          wtsum = 0.
          DO ipr = 1, ngn(igcsm)
            iprsm = iprsm + 1
            wtsum = wtsum + wt(iprsm)
          ENDDO
          wtsm(igc) = wtsum
        ENDDO
        DO ig = 1, ng(ibnd+15)
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = wt(ig)/wtsm(ngm(ind))
        ENDDO
      ELSE
        DO ig = 1, ng(ibnd+15)
          igcsm = igcsm + 1
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = 1.0_wp
        ENDDO
      ENDIF
    ENDDO

    ! Reduce g-points for absorption coefficient data in each LW spectral band.

    CALL cmbgb16s
    CALL cmbgb17
    CALL cmbgb18
    CALL cmbgb19
    CALL cmbgb20
    CALL cmbgb21
    CALL cmbgb22
    CALL cmbgb23
    CALL cmbgb24
    CALL cmbgb25
    CALL cmbgb26
    CALL cmbgb27
    CALL cmbgb28
    CALL cmbgb29

  END SUBROUTINE setup_srtm

  !***************************************************************************
  SUBROUTINE cmbgb16s
    !***************************************************************************
    !
    !  Original version:       MJIacono; July 1998
    !  Revision for RRTM_SW:   MJIacono; November 2002
    !  Revision for RRTMG_SW:  MJIacono; December 2003
    !  Revision for F90 reformatting:  MJIacono; July 2006
    !
    !  The subroutines CMBGB16->CMBGB29 input the absorption coefficient
    !  data for each band, which are defined for 16 g-points and 14 spectral
    !  bands. The data are combined with appropriate weighting following the
    !  g-point mapping arrays specified in RRTMG_SW_INIT.  Solar source 
    !  function data in array SFLUXREF are combined without weighting.  All
    !  g-point reduced data are put into new arrays for use in RRTMG_SW.
    !
    !  band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
    !
    !-----------------------------------------------------------------------

    USE rrsw_kg16, ONLY : kao, kbo, selfrefo, forrefo, sfluxrefo, &
         ka, kb, selfref, forref, sfluxref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(1)
            sumk = 0.
            DO ipr = 1, ngn(igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(1)
          sumk = 0.
          DO ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(1)
        sumk = 0.
        DO ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(1)
        sumk = 0.
        DO ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(1)
      sumf = 0.
      DO ipr = 1, ngn(igc)
        iprsm = iprsm + 1
        sumf = sumf + sfluxrefo(iprsm)
      ENDDO
      sfluxref(igc) = sumf
    ENDDO

  END SUBROUTINE cmbgb16s

  !***************************************************************************
  SUBROUTINE cmbgb17
    !***************************************************************************
    !
    !     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
    !-----------------------------------------------------------------------

    USE rrsw_kg17, ONLY : kao, kbo, selfrefo, forrefo, sfluxrefo, &
         ka, kb, selfref, forref, sfluxref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(2)
            sumk = 0.
            DO ipr = 1, ngn(ngs(1)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+16)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,5
      DO jt = 1,5
        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(2)
            sumk = 0.
            DO ipr = 1, ngn(ngs(1)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+16)
            ENDDO
            kb(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(2)
        sumk = 0.
        DO ipr = 1, ngn(ngs(1)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+16)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(2)
        sumk = 0.
        DO ipr = 1, ngn(ngs(1)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+16)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,5
      iprsm = 0
      DO igc = 1,ngc(2)
        sumf = 0.
        DO ipr = 1, ngn(ngs(1)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb17

  !***************************************************************************
  SUBROUTINE cmbgb18
    !***************************************************************************
    !
    !     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
    !-----------------------------------------------------------------------

    USE rrsw_kg18, ONLY : kao, kbo, selfrefo, forrefo, sfluxrefo, &
         ka, kb, selfref, forref, sfluxref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(3)
            sumk = 0.
            DO ipr = 1, ngn(ngs(2)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+32)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(3)
          sumk = 0.
          DO ipr = 1, ngn(ngs(2)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+32)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(3)
        sumk = 0.
        DO ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+32)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(3)
        sumk = 0.
        DO ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+32)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(3)
        sumf = 0.
        DO ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb18

  !***************************************************************************
  SUBROUTINE cmbgb19
    !***************************************************************************
    !
    !     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
    !-----------------------------------------------------------------------

    USE rrsw_kg19, ONLY : kao, kbo, selfrefo, forrefo, sfluxrefo, &
         ka, kb, selfref, forref, sfluxref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(4)
            sumk = 0.
            DO ipr = 1, ngn(ngs(3)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+48)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(4)
          sumk = 0.
          DO ipr = 1, ngn(ngs(3)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+48)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(4)
        sumk = 0.
        DO ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+48)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(4)
        sumk = 0.
        DO ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+48)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(4)
        sumf = 0.
        DO ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb19

  !***************************************************************************
  SUBROUTINE cmbgb20
    !***************************************************************************
    !
    !     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
    !-----------------------------------------------------------------------

    USE rrsw_kg20, ONLY : kao, kbo, selfrefo, forrefo, sfluxrefo, absch4o, &
         ka, kb, selfref, forref, sfluxref, absch4

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(5)
          sumk = 0.
          DO ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+64)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(5)
          sumk = 0.
          DO ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+64)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(5)
        sumk = 0.
        DO ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+64)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(5)
        sumk = 0.
        DO ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+64)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(5)
      sumf1 = 0.
      sumf2 = 0.
      DO ipr = 1, ngn(ngs(4)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + absch4o(iprsm)*rwgt(iprsm+64)
      ENDDO
      sfluxref(igc) = sumf1
      absch4(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb20

  !***************************************************************************
  SUBROUTINE cmbgb21
    !***************************************************************************
    !
    !     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
    !-----------------------------------------------------------------------

    USE rrsw_kg21, ONLY : kao, kbo, selfrefo, forrefo, sfluxrefo, &
         ka, kb, selfref, forref, sfluxref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(6)
            sumk = 0.
            DO ipr = 1, ngn(ngs(5)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+80)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,5
      DO jt = 1,5
        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(6)
            sumk = 0.
            DO ipr = 1, ngn(ngs(5)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+80)
            ENDDO
            kb(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(6)
        sumk = 0.
        DO ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+80)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(6)
        sumk = 0.
        DO ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+80)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(6)
        sumf = 0.
        DO ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb21

  !***************************************************************************
  SUBROUTINE cmbgb22
    !***************************************************************************
    !
    !     band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)
    !-----------------------------------------------------------------------

    USE rrsw_kg22, ONLY : kao, kbo, selfrefo, forrefo, sfluxrefo, &
         ka, kb, selfref, forref, sfluxref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(7)
            sumk = 0.
            DO ipr = 1, ngn(ngs(6)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+96)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(7)
          sumk = 0.
          DO ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+96)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(7)
        sumk = 0.
        DO ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+96)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(7)
        sumk = 0.
        DO ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+96)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(7)
        sumf = 0.
        DO ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb22

  !***************************************************************************
  SUBROUTINE cmbgb23
    !***************************************************************************
    !
    !     band 23:  8050-12850 cm-1 (low - h2o; high - nothing)
    !-----------------------------------------------------------------------

    USE rrsw_kg23, ONLY : kao, selfrefo, forrefo, sfluxrefo, raylo, &
         ka, selfref, forref, sfluxref, rayl

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(8)
          sumk = 0.
          DO ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+112)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(8)
        sumk = 0.
        DO ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+112)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(8)
        sumk = 0.
        DO ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+112)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(8)
      sumf1 = 0.
      sumf2 = 0.
      DO ipr = 1, ngn(ngs(7)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + raylo(iprsm)*rwgt(iprsm+112)
      ENDDO
      sfluxref(igc) = sumf1
      rayl(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb23

  !***************************************************************************
  SUBROUTINE cmbgb24
    !***************************************************************************
    !
    !     band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)
    !-----------------------------------------------------------------------

    USE rrsw_kg24, ONLY : kao, kbo, selfrefo, forrefo, sfluxrefo, &
         abso3ao, abso3bo, raylao, raylbo, &
         ka, kb, selfref, forref, sfluxref, &
         abso3a, abso3b, rayla, raylb

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2, sumf3


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(9)
            sumk = 0.
            DO ipr = 1, ngn(ngs(8)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+128)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(9)
          sumk = 0.
          DO ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+128)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(9)
        sumk = 0.
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+128)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,3
      iprsm = 0
      DO igc = 1,ngc(9)
        sumk = 0.
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+128)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(9)
      sumf1 = 0.
      sumf2 = 0.
      sumf3 = 0.
      DO ipr = 1, ngn(ngs(8)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + raylbo(iprsm)*rwgt(iprsm+128)
        sumf2 = sumf2 + abso3ao(iprsm)*rwgt(iprsm+128)
        sumf3 = sumf3 + abso3bo(iprsm)*rwgt(iprsm+128)
      ENDDO
      raylb(igc) = sumf1
      abso3a(igc) = sumf2
      abso3b(igc) = sumf3
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(9)
        sumf1 = 0.
        sumf2 = 0.
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumf1 = sumf1 + sfluxrefo(iprsm,jp)
          sumf2 = sumf2 + raylao(iprsm,jp)*rwgt(iprsm+128)
        ENDDO
        sfluxref(igc,jp) = sumf1
        rayla(igc,jp) = sumf2
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb24

  !***************************************************************************
  SUBROUTINE cmbgb25
    !***************************************************************************
    !
    !     band 25:  16000-22650 cm-1 (low - h2o; high - nothing)
    !-----------------------------------------------------------------------

    USE rrsw_kg25, ONLY : kao, sfluxrefo, &
         abso3ao, abso3bo, raylo, &
         ka, sfluxref, &
         abso3a, abso3b, rayl

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2, sumf3, sumf4


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(10)
          sumk = 0.
          DO ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+144)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(10)
      sumf1 = 0.
      sumf2 = 0.
      sumf3 = 0.
      sumf4 = 0.
      DO ipr = 1, ngn(ngs(9)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + abso3ao(iprsm)*rwgt(iprsm+144)
        sumf3 = sumf3 + abso3bo(iprsm)*rwgt(iprsm+144)
        sumf4 = sumf4 + raylo(iprsm)*rwgt(iprsm+144)
      ENDDO
      sfluxref(igc) = sumf1
      abso3a(igc) = sumf2
      abso3b(igc) = sumf3
      rayl(igc) = sumf4
    ENDDO

  END SUBROUTINE cmbgb25

  !***************************************************************************
  SUBROUTINE cmbgb26
    !***************************************************************************
    !
    !     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)
    !-----------------------------------------------------------------------

    USE rrsw_kg26, ONLY : sfluxrefo, raylo, &
         sfluxref, rayl

    ! ------- Local -------
    INTEGER :: igc, ipr, iprsm
    REAL(wp) :: sumf1, sumf2


    iprsm = 0
    DO igc = 1,ngc(11)
      sumf1 = 0.
      sumf2 = 0.
      DO ipr = 1, ngn(ngs(10)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + raylo(iprsm)*rwgt(iprsm+160)
        sumf2 = sumf2 + sfluxrefo(iprsm)
      ENDDO
      rayl(igc) = sumf1
      sfluxref(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb26

  !***************************************************************************
  SUBROUTINE cmbgb27
    !***************************************************************************
    !
    !     band 27:  29000-38000 cm-1 (low - o3; high - o3)
    !-----------------------------------------------------------------------

    USE rrsw_kg27, ONLY : kao, kbo, sfluxrefo, raylo, &
         ka, kb, sfluxref, rayl

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(12)
          sumk = 0.
          DO ipr = 1, ngn(ngs(11)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+176)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(12)
          sumk = 0.
          DO ipr = 1, ngn(ngs(11)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+176)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(12)
      sumf1 = 0.
      sumf2 = 0.
      DO ipr = 1, ngn(ngs(11)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + raylo(iprsm)*rwgt(iprsm+176)
      ENDDO
      sfluxref(igc) = sumf1
      rayl(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb27

  !***************************************************************************
  SUBROUTINE cmbgb28
    !***************************************************************************
    !
    !     band 28:  38000-50000 cm-1 (low - o3,o2; high - o3,o2)
    !-----------------------------------------------------------------------

    USE rrsw_kg28, ONLY : kao, kbo, sfluxrefo, &
         ka, kb, sfluxref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(13)
            sumk = 0.
            DO ipr = 1, ngn(ngs(12)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+192)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,5
      DO jt = 1,5
        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(13)
            sumk = 0.
            DO ipr = 1, ngn(ngs(12)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+192)
            ENDDO
            kb(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jp = 1,5
      iprsm = 0
      DO igc = 1,ngc(13)
        sumf = 0.
        DO ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumf = sumf + sfluxrefo(iprsm,jp)
        ENDDO
        sfluxref(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb28

  !***************************************************************************
  SUBROUTINE cmbgb29
    !***************************************************************************
    !
    !     band 29:  820-2600 cm-1 (low - h2o; high - co2)
    !-----------------------------------------------------------------------

    USE rrsw_kg29, ONLY : kao, kbo, selfrefo, forrefo, sfluxrefo, &
         absh2oo, absco2o, &
         ka, kb, selfref, forref, sfluxref, &
         absh2o, absco2

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm
    REAL(wp) :: sumk, sumf1, sumf2, sumf3


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(14)
          sumk = 0.
          DO ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+208)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(14)
          sumk = 0.
          DO ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+208)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(14)
        sumk = 0.
        DO ipr = 1, ngn(ngs(13)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+208)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(14)
        sumk = 0.
        DO ipr = 1, ngn(ngs(13)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+208)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(14)
      sumf1 = 0.
      sumf2 = 0.
      sumf3 = 0.
      DO ipr = 1, ngn(ngs(13)+igc)
        iprsm = iprsm + 1
        sumf1 = sumf1 + sfluxrefo(iprsm)
        sumf2 = sumf2 + absco2o(iprsm)*rwgt(iprsm+208)
        sumf3 = sumf3 + absh2oo(iprsm)*rwgt(iprsm+208)
      ENDDO
      sfluxref(igc) = sumf1
      absco2(igc) = sumf2
      absh2o(igc) = sumf3
    ENDDO

  END SUBROUTINE cmbgb29

END MODULE mo_srtm_setup


