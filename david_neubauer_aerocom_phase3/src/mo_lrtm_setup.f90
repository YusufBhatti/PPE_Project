!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_lrtm_setup

  USE mo_kind, ONLY : wp
  USE mo_rrtm_params, ONLY : mg, nbndlw, ngptlw, maxinpx
  USE mo_lrtm_netcdf, ONLY : lrtm_read
  USE mo_rad_fastmath,ONLY : setup_rad_fastmath

  IMPLICIT NONE
  !
  ! spectra information that is entered at run time
  !
  REAL(wp) :: rwgt(nbndlw*mg) !< Weights for combining original gpts into reduced gpts

  INTEGER :: nxmol            !< Number of cross-section molecules
  INTEGER :: ixindx(maxinpx)  !< Flag for active cross-sections in calculation

  INTEGER, PARAMETER :: ngc(nbndlw) = (/    & !< The number of new g-intervals in each band
       10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2/)
  INTEGER, PARAMETER :: ngs(nbndlw) = (/    & !< The cumulative sum of new g-intervals for each band
       10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,140/)
  INTEGER, PARAMETER :: ngm(nbndlw*mg) = (/ & !< The index of each new gpt relative to the orignal
       1,2,3,3,4,4,5,5,6,6,7,7,8,8,9,10, &          ! band 1
       1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 2
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 3
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,14, &    ! band 4
       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 5
       1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 6
       1,1,2,2,3,4,5,6,7,8,9,10,11,11,12,12, &      ! band 7
       1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 8
       1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 9
       1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 10
       1,2,3,3,4,4,5,5,6,6,7,7,7,8,8,8, &           ! band 11
       1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 12
       1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4, &           ! band 13
       1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 14
       1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 15
       1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2/)            ! band 16
  INTEGER, PARAMETER :: ngn(ngptlw) = (/    & !< The number of original gs combined to make new pts
       1,1,2,2,2,2,2,2,1,1, &                       ! band 1
       1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 2
       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 3
       1,1,1,1,1,1,1,1,1,1,1,1,1,3, &               ! band 4
       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 5
       2,2,2,2,2,2,2,2, &                           ! band 6
       2,2,1,1,1,1,1,1,1,1,2,2, &                   ! band 7
       2,2,2,2,2,2,2,2, &                           ! band 8
       1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 9
       2,2,2,2,4,4, &                               ! band 10
       1,1,2,2,2,2,3,3, &                           ! band 11
       1,1,1,1,2,2,4,4, &                           ! band 12
       3,3,4,6, &                                   ! band 13
       8,8, &                                       ! band 14
       8,8, &                                       ! band 15
       4,12/)                                       ! band 16
  INTEGER, PARAMETER :: ngb(ngptlw) = (/ & !< The band index for each new g-interval
       1,1,1,1,1,1,1,1,1,1, &                       ! band 1
       2,2,2,2,2,2,2,2,2,2,2,2, &                   ! band 2
       3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, &           ! band 3
       4,4,4,4,4,4,4,4,4,4,4,4,4,4, &               ! band 4
       5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, &           ! band 5
       6,6,6,6,6,6,6,6, &                           ! band 6
       7,7,7,7,7,7,7,7,7,7,7,7, &                   ! band 7
       8,8,8,8,8,8,8,8, &                           ! band 8
       9,9,9,9,9,9,9,9,9,9,9,9, &                   ! band 9
       10,10,10,10,10,10, &                         ! band 10
       11,11,11,11,11,11,11,11, &                   ! band 11
       12,12,12,12,12,12,12,12, &                   ! band 12
       13,13,13,13, &                               ! band 13
       14,14, &                                     ! band 14
       15,15, &                                     ! band 15
       16,16/)                                      ! band 16

  REAL(wp), PARAMETER :: wt(mg) = (/ & !< RRTM weights for the original 16 g-intervals
       0.1527534276_wp, 0.1491729617_wp, 0.1420961469_wp, &
       0.1316886544_wp, 0.1181945205_wp, 0.1019300893_wp, &
       0.0832767040_wp, 0.0626720116_wp, 0.0424925000_wp, &
       0.0046269894_wp, 0.0038279891_wp, 0.0030260086_wp, &
       0.0022199750_wp, 0.0014140010_wp, 0.0005330000_wp, &
       0.0000750000_wp/)

  INTEGER, PARAMETER :: nspa(nbndlw) = (/ & !< Number of reference atmospheres for lower atmosphere
       1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9/)  
  INTEGER, PARAMETER :: nspb(nbndlw) = (/ & !< Number of reference atmospheres for upper atmosphere
       1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0/)
  INTEGER, PARAMETER :: ng(nbndlw)   = (/ &! < Number of g intervals in each band
       16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/)

  REAL(wp), PARAMETER :: wavenum1(nbndlw) = (/ & !< Spectral band lower boundary in wavenumbers
       &   10._wp, 350._wp, 500._wp, 630._wp, 700._wp, 820._wp, &
       &  980._wp,1080._wp,1180._wp,1390._wp,1480._wp,1800._wp, &
       & 2080._wp,2250._wp,2380._wp,2600._wp/)
  REAL(wp), PARAMETER :: wavenum2(nbndlw) = (/ & !< Spectral band upper boundary in wavenumbers
       &  350._wp, 500._wp, 630._wp, 700._wp, 820._wp, 980._wp, &
       & 1080._wp,1180._wp,1390._wp,1480._wp,1800._wp,2080._wp, &
       & 2250._wp,2380._wp,2600._wp,3250._wp/)
  REAL(wp), PARAMETER :: delwave(nbndlw)  = (/ & !< Spectral band width in wavenumbers
       340._wp, 150._wp, 130._wp,  70._wp, 120._wp, 160._wp, &
       100._wp, 100._wp, 210._wp,  90._wp, 320._wp, 280._wp, &
       170._wp, 130._wp, 220._wp, 650._wp/)

CONTAINS

  ! **************************************************************************
  SUBROUTINE setup_lrtm
    ! **************************************************************************
    !
    !  Original version:       Michael J. Iacono; July, 1998
    !  First revision for GCMs:   September, 1998
    !  Second revision for RRTM_V3.0:  September, 2002
    !  M. ESCH, MPI-M, 30-May-2012: Write keywords in capital letters
    !
    !  This subroutine performs calculations necessary for the initialization
    !  of the longwave model.  Lookup tables are computed for use in the LW
    !  radiative transfer, and input absorption coefficient data for each
    !  spectral band are reduced from 256 g-point intervals to 140.
    ! **************************************************************************

    ! ------- Local -------

    INTEGER :: ibnd, igc, ig, ind, ipr 
    INTEGER :: igcsm, iprsm

    REAL(wp) :: wtsum, wtsm(mg)        !

    ! Initialize model data
    !
    !     nxmol     - number of cross-sections input by user
    !     ixindx(i) - index of cross-section molecule corresponding to Ith
    !                 cross-section specified by user
    !
    nxmol = 4
    ixindx(1) = 1 ! ccl4
    ixindx(2) = 2 ! cfc11
    ixindx(3) = 3 ! cfc12
    ixindx(4) = 4 ! cfc22
    ixindx(5:maxinpx) = 0

    CALL lrtm_read              ! molecular absorption coefficients
    CALL setup_rad_fastmath

    ! Perform g-point reduction from 16 per band (256 total points) to
    ! a band dependant number (140 total points) for all absorption
    ! coefficient input data and Planck fraction input data.
    ! Compute relative weighting for new g-point combinations.

    igcsm = 0
    DO ibnd = 1,nbndlw
      iprsm = 0
      IF (ngc(ibnd).LT.mg) THEN
        DO igc = 1,ngc(ibnd) 
          igcsm = igcsm + 1
          wtsum = 0._wp
          DO ipr = 1, ngn(igcsm)
            iprsm = iprsm + 1
            wtsum = wtsum + wt(iprsm)
          ENDDO
          wtsm(igc) = wtsum
        ENDDO
        DO ig = 1, ng(ibnd)
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = wt(ig)/wtsm(ngm(ind))
        ENDDO
      ELSE
        DO ig = 1, ng(ibnd)
          igcsm = igcsm + 1
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = 1.0_wp
        ENDDO
      ENDIF
    ENDDO

    ! Reduce g-points for absorption coefficient data in each LW spectral band.

    CALL cmbgb1
    CALL cmbgb2
    CALL cmbgb3
    CALL cmbgb4
    CALL cmbgb5
    CALL cmbgb6
    CALL cmbgb7
    CALL cmbgb8
    CALL cmbgb9
    CALL cmbgb10
    CALL cmbgb11
    CALL cmbgb12
    CALL cmbgb13
    CALL cmbgb14
    CALL cmbgb15
    CALL cmbgb16

  END SUBROUTINE setup_lrtm

  !***************************************************************************
  SUBROUTINE cmbgb1
    !  Original version:    MJIacono; July 1998
    !  Revision for GCMs:   MJIacono; September 1998
    !  Revision for RRTMG:  MJIacono, September 2002
    !  Revision for F90 reformatting:  MJIacono, June 2006
    !  M. ESCH, MPI-M, 30-May-2012: Write keywords in capital letters
    !
    !  The subroutines CMBGB1->CMBGB16 input the absorption coefficient
    !  data for each band, which are defined for 16 g-points and 16 spectral
    !  bands. The data are combined with appropriate weighting following the
    !  g-point mapping arrays specified in RRTMINIT.  Plank fraction data
    !  in arrays FRACREFA and FRACREFB are combined without weighting.  All
    !  g-point reduced data are put into new arrays for use in RRTM.
    !
    !  band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
    !                       (high key - h2o; high minor - n2)
    !  note: previous versions of rrtm band 1: 
    !        10-250 cm-1 (low - h2o; high - h2o)
    !***************************************************************************

    USE rrlw_kg01, ONLY: fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, ka_mn2, kb_mn2, &
         selfref, forref

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumk1, sumk2, sumf1, sumf2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(1)
          sumk = 0.0_wp
          DO ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(1)
          sumk = 0.0_wp
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
        sumk = 0.0_wp
        DO ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(1)
        sumk = 0.0_wp
        DO ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,19
      iprsm = 0
      DO igc = 1,ngc(1)
        sumk1 = 0.0_wp
        sumk2 = 0.0_wp
        DO ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk1 = sumk1 + kao_mn2(jt,iprsm)*rwgt(iprsm)
          sumk2 = sumk2 + kbo_mn2(jt,iprsm)*rwgt(iprsm)
        ENDDO
        ka_mn2(jt,igc) = sumk1
        kb_mn2(jt,igc) = sumk2
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(1)
      sumf1 = 0.0_wp
      sumf2 = 0.0_wp
      DO ipr = 1, ngn(igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      ENDDO
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb1

  !***************************************************************************
  SUBROUTINE cmbgb2
    !***************************************************************************
    !
    !     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
    !
    !     note: previous version of rrtm band 2: 
    !           250 - 500 cm-1 (low - h2o; high - h2o)
    !***************************************************************************

    USE rrlw_kg02, ONLY: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, selfref, forref

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf1, sumf2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(2)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+16)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(2)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+16)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(2)
        sumk = 0.0_wp
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
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(1)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+16)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(2)
      sumf1 = 0.0_wp
      sumf2 = 0.0_wp
      DO ipr = 1, ngn(ngs(1)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      ENDDO
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb2

  !***************************************************************************
  SUBROUTINE cmbgb3
    !***************************************************************************
    !
    !     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
    !                           (high key - h2o,co2; high minor - n2o)
    !
    ! old band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
    !***************************************************************************

    USE rrlw_kg03, ONLY: fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, ka_mn2o, kb_mn2o, &
         selfref, forref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(3)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(2)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+32)
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
          DO igc = 1,ngc(3)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(2)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+32)
            ENDDO
            kb(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,9
      DO jt = 1,19
        iprsm = 0
        DO igc = 1,ngc(3)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(2)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
          ENDDO
          ka_mn2o(jn,jt,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,5
      DO jt = 1,19
        iprsm = 0
        DO igc = 1,ngc(3)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(2)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
          ENDDO
          kb_mn2o(jn,jt,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(3)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+32)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(3)
        sumk = 0.0_wp
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
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        ENDDO
        fracrefa(igc,jp) = sumf
      ENDDO
    ENDDO

    DO jp = 1,5
      iprsm = 0
      DO igc = 1,ngc(3)
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefbo(iprsm,jp)
        ENDDO
        fracrefb(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb3

  !***************************************************************************
  SUBROUTINE cmbgb4
    !***************************************************************************
    !
    !     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
    !
    ! old band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
    !***************************************************************************

    USE rrlw_kg04, ONLY: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, selfref, forref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(4)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(3)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+48)
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
          DO igc = 1,ngc(4)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(3)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+48)
            ENDDO
            kb(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(4)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+48)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(4)
        sumk = 0.0_wp
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
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        ENDDO
        fracrefa(igc,jp) = sumf
      ENDDO
    ENDDO

    DO jp = 1,5
      iprsm = 0
      DO igc = 1,ngc(4)
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefbo(iprsm,jp)
        ENDDO
        fracrefb(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb4

  !***************************************************************************
  SUBROUTINE cmbgb5
    !***************************************************************************
    !
    !     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
    !                           (high key - o3,co2)
    !
    ! old band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
    !***************************************************************************

    USE rrlw_kg05, ONLY: fracrefao, fracrefbo, kao, kbo, kao_mo3, ccl4o, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, ka_mo3, ccl4, &
         selfref, forref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(5)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(4)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+64)
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
          DO igc = 1,ngc(5)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(4)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+64)
            ENDDO
            kb(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,9
      DO jt = 1,19
        iprsm = 0
        DO igc = 1,ngc(5)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mo3(jn,jt,iprsm)*rwgt(iprsm+64)
          ENDDO
          ka_mo3(jn,jt,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(5)
        sumk = 0.0_wp
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
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+64)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(5)
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        ENDDO
        fracrefa(igc,jp) = sumf
      ENDDO
    ENDDO

    DO jp = 1,5
      iprsm = 0
      DO igc = 1,ngc(5)
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefbo(iprsm,jp)
        ENDDO
        fracrefb(igc,jp) = sumf
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(5)
      sumk = 0.0_wp
      DO ipr = 1, ngn(ngs(4)+igc)
        iprsm = iprsm + 1
        sumk = sumk + ccl4o(iprsm)*rwgt(iprsm+64)
      ENDDO
      ccl4(igc) = sumk
    ENDDO

  END SUBROUTINE cmbgb5

  !***************************************************************************
  SUBROUTINE cmbgb6
    !***************************************************************************
    !
    !     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
    !                           (high key - nothing; high minor - cfc11, cfc12)
    !
    ! old band 6:  820-980 cm-1 (low - h2o; high - nothing)
    !***************************************************************************

    USE rrlw_kg06, ONLY: fracrefao, kao, kao_mco2, cfc11adjo, cfc12o, &
         selfrefo, forrefo, &
         fracrefa, ka, ka_mco2, cfc11adj, cfc12, &
         selfref, forref

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf, sumk1, sumk2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(6)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(5)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+80)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,19
      iprsm = 0
      DO igc = 1,ngc(6)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumk = sumk + kao_mco2(jt,iprsm)*rwgt(iprsm+80)
        ENDDO
        ka_mco2(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(6)
        sumk = 0.0_wp
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
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+80)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(6)
      sumf = 0.0_wp
      sumk1= 0.0_wp
      sumk2= 0.0_wp
      DO ipr = 1, ngn(ngs(5)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefao(iprsm)
        sumk1= sumk1+ cfc11adjo(iprsm)*rwgt(iprsm+80)
        sumk2= sumk2+ cfc12o(iprsm)*rwgt(iprsm+80)
      ENDDO
      fracrefa(igc) = sumf
      cfc11adj(igc) = sumk1
      cfc12(igc) = sumk2
    ENDDO

  END SUBROUTINE cmbgb6

  !***************************************************************************
  SUBROUTINE cmbgb7
    !***************************************************************************
    !
    !     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
    !                            (high key - o3; high minor - co2)
    !
    ! old band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
    !***************************************************************************

    USE rrlw_kg07, ONLY: fracrefao, fracrefbo, kao, kbo, kao_mco2, kbo_mco2, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, ka_mco2, kb_mco2, &
         selfref, forref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(7)
            sumk = 0.0_wp
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
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+96)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,9
      DO jt = 1,19
        iprsm = 0
        DO igc = 1,ngc(7)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+96)
          ENDDO
          ka_mco2(jn,jt,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,19
      iprsm = 0
      DO igc = 1,ngc(7)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumk = sumk + kbo_mco2(jt,iprsm)*rwgt(iprsm+96)
        ENDDO
        kb_mco2(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(7)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+96)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(7)
        sumk = 0.0_wp
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
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        ENDDO
        fracrefa(igc,jp) = sumf
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(7)
      sumf = 0.0_wp
      DO ipr = 1, ngn(ngs(6)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefbo(iprsm)
      ENDDO
      fracrefb(igc) = sumf
    ENDDO

  END SUBROUTINE cmbgb7

  !***************************************************************************
  SUBROUTINE cmbgb8
    !***************************************************************************
    !
    !     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
    !                             (high key - o3; high minor - co2, n2o)
    !
    ! old band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
    !***************************************************************************

    USE rrlw_kg08, ONLY: fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
         kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
         cfc12o, cfc22adjo, &
         fracrefa, fracrefb, ka, ka_mco2, ka_mn2o, &
         ka_mo3, kb, kb_mco2, kb_mn2o, selfref, forref, &
         cfc12, cfc22adj

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumk1, sumk2, sumk3, sumk4, sumk5, sumf1, sumf2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(8)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+112)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO
    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(8)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+112)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(8)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+112)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(8)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+112)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,19
      iprsm = 0
      DO igc = 1,ngc(8)
        sumk1 = 0.0_wp
        sumk2 = 0.0_wp
        sumk3 = 0.0_wp
        sumk4 = 0.0_wp
        sumk5 = 0.0_wp
        DO ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk1 = sumk1 + kao_mco2(jt,iprsm)*rwgt(iprsm+112)
          sumk2 = sumk2 + kbo_mco2(jt,iprsm)*rwgt(iprsm+112)
          sumk3 = sumk3 + kao_mo3(jt,iprsm)*rwgt(iprsm+112)
          sumk4 = sumk4 + kao_mn2o(jt,iprsm)*rwgt(iprsm+112)
          sumk5 = sumk5 + kbo_mn2o(jt,iprsm)*rwgt(iprsm+112)
        ENDDO
        ka_mco2(jt,igc) = sumk1
        kb_mco2(jt,igc) = sumk2
        ka_mo3(jt,igc) = sumk3
        ka_mn2o(jt,igc) = sumk4
        kb_mn2o(jt,igc) = sumk5
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(8)
      sumf1= 0.0_wp
      sumf2= 0.0_wp
      sumk1= 0.0_wp
      sumk2= 0.0_wp
      DO ipr = 1, ngn(ngs(7)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
        sumk1= sumk1+ cfc12o(iprsm)*rwgt(iprsm+112)
        sumk2= sumk2+ cfc22adjo(iprsm)*rwgt(iprsm+112)
      ENDDO
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
      cfc12(igc) = sumk1
      cfc22adj(igc) = sumk2
    ENDDO

  END SUBROUTINE cmbgb8

  !***************************************************************************
  SUBROUTINE cmbgb9
    !***************************************************************************
    !
    !     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
    !                             (high key - ch4; high minor - n2o)!

    ! old band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
    !***************************************************************************

    USE rrlw_kg09, ONLY: fracrefao, fracrefbo, kao, kao_mn2o, &
         kbo, kbo_mn2o, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, ka_mn2o, &
         kb, kb_mn2o, selfref, forref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(9)
            sumk = 0.0_wp
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
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+128)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,9
      DO jt = 1,19
        iprsm = 0
        DO igc = 1,ngc(9)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+128)
          ENDDO
          ka_mn2o(jn,jt,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,19
      iprsm = 0
      DO igc = 1,ngc(9)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + kbo_mn2o(jt,iprsm)*rwgt(iprsm+128)
        ENDDO
        kb_mn2o(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(9)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+128)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(9)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+128)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(9)
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        ENDDO
        fracrefa(igc,jp) = sumf
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(9)
      sumf = 0.0_wp
      DO ipr = 1, ngn(ngs(8)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefbo(iprsm)
      ENDDO
      fracrefb(igc) = sumf
    ENDDO

  END SUBROUTINE cmbgb9

  !***************************************************************************
  SUBROUTINE cmbgb10
    !***************************************************************************
    !
    !     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
    !
    ! old band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
    !***************************************************************************

    USE rrlw_kg10, ONLY: fracrefao, fracrefbo, kao, kbo, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, &
         selfref, forref

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf1, sumf2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(10)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+144)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(10)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+144)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(10)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(9)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+144)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(10)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(9)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+144)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(10)
      sumf1= 0.0_wp
      sumf2= 0.0_wp
      DO ipr = 1, ngn(ngs(9)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      ENDDO
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb10

  !***************************************************************************
  SUBROUTINE cmbgb11
    !***************************************************************************
    !
    !     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
    !                              (high key - h2o; high minor - o2)
    !
    ! old band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
    !                              (high key - h2o; high minor - o2)
    !***************************************************************************

    USE rrlw_kg11, ONLY: fracrefao, fracrefbo, kao, kao_mo2, &
         kbo, kbo_mo2, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, ka_mo2, &
         kb, kb_mo2, selfref, forref

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumk1, sumk2, sumf1, sumf2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(11)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(10)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+160)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO
    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(11)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(10)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+160)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,19
      iprsm = 0
      DO igc = 1,ngc(11)
        sumk1 = 0.0_wp
        sumk2 = 0.0_wp
        DO ipr = 1, ngn(ngs(10)+igc)
          iprsm = iprsm + 1
          sumk1 = sumk1 + kao_mo2(jt,iprsm)*rwgt(iprsm+160)
          sumk2 = sumk2 + kbo_mo2(jt,iprsm)*rwgt(iprsm+160)
        ENDDO
        ka_mo2(jt,igc) = sumk1
        kb_mo2(jt,igc) = sumk2
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(11)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(10)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+160)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(11)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(10)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+160)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(11)
      sumf1= 0.0_wp
      sumf2= 0.0_wp
      DO ipr = 1, ngn(ngs(10)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      ENDDO
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb11

  !***************************************************************************
  SUBROUTINE cmbgb12
    !***************************************************************************
    !
    !     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
    !
    ! old band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
    !***************************************************************************

    USE rrlw_kg12, ONLY: fracrefao, kao, selfrefo, forrefo, &
         fracrefa, ka, selfref, forref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(12)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(11)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+176)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(12)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(11)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+176)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(12)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(11)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+176)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(12)
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(11)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        ENDDO
        fracrefa(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb12

  !***************************************************************************
  SUBROUTINE cmbgb13
    !***************************************************************************
    !
    !     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
    !
    ! old band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
    !***************************************************************************

    USE rrlw_kg13, ONLY: fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
         kbo_mo3, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, ka_mco2, ka_mco, &
         kb_mo3, selfref, forref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumk1, sumk2, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(13)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(12)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+192)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,9
      DO jt = 1,19
        iprsm = 0
        DO igc = 1,ngc(13)
          sumk1 = 0.0_wp
          sumk2 = 0.0_wp
          DO ipr = 1, ngn(ngs(12)+igc)
            iprsm = iprsm + 1
            sumk1 = sumk1 + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+192)
            sumk2 = sumk2 + kao_mco(jn,jt,iprsm)*rwgt(iprsm+192)
          ENDDO
          ka_mco2(jn,jt,igc) = sumk1
          ka_mco(jn,jt,igc) = sumk2
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,19
      iprsm = 0
      DO igc = 1,ngc(13)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumk = sumk + kbo_mo3(jt,iprsm)*rwgt(iprsm+192)
        ENDDO
        kb_mo3(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(13)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+192)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(13)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+192)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(13)
      sumf = 0.0_wp
      DO ipr = 1, ngn(ngs(12)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefbo(iprsm)
      ENDDO
      fracrefb(igc) = sumf
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(13)
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        ENDDO
        fracrefa(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb13

  !***************************************************************************
  SUBROUTINE cmbgb14
    !***************************************************************************
    !
    !     band 14:  2250-2380 cm-1 (low - co2; high - co2)
    !
    ! old band 14:  2250-2380 cm-1 (low - co2; high - co2)
    !***************************************************************************

    USE rrlw_kg14, ONLY: fracrefao, fracrefbo, kao, kbo, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, &
         selfref, forref

    ! ------- Local -------
    INTEGER :: jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf1, sumf2


    DO jt = 1,5
      DO jp = 1,13
        iprsm = 0
        DO igc = 1,ngc(14)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+208)
          ENDDO
          ka(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(14)
          sumk = 0.0_wp
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
        sumk = 0.0_wp
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
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(13)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+208)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(14)
      sumf1= 0.0_wp
      sumf2= 0.0_wp
      DO ipr = 1, ngn(ngs(13)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      ENDDO
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    ENDDO

  END SUBROUTINE cmbgb14

  !***************************************************************************
  SUBROUTINE cmbgb15
    !***************************************************************************
    !
    !     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
    !                              (high - nothing)
    !
    ! old band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
    !***************************************************************************

    USE rrlw_kg15, ONLY: fracrefao, kao, kao_mn2, selfrefo, forrefo, &
         fracrefa, ka, ka_mn2, selfref, forref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(15)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(14)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+224)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jn = 1,9
      DO jt = 1,19
        iprsm = 0
        DO igc = 1,ngc(15)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(14)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mn2(jn,jt,iprsm)*rwgt(iprsm+224)
          ENDDO
          ka_mn2(jn,jt,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(15)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(14)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+224)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(15)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(14)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+224)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(15)
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(14)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        ENDDO
        fracrefa(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb15

  !***************************************************************************
  SUBROUTINE cmbgb16
    !***************************************************************************
    !
    !     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
    !
    ! old band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
    !***************************************************************************

    USE rrlw_kg16, ONLY: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, selfref, forref

    ! ------- Local -------
    INTEGER :: jn, jt, jp, igc, ipr, iprsm 
    REAL(wp) :: sumk, sumf


    DO jn = 1,9
      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(16)
            sumk = 0.0_wp
            DO ipr = 1, ngn(ngs(15)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+240)
            ENDDO
            ka(jn,jt,jp,igc) = sumk
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,5
      DO jp = 13,59
        iprsm = 0
        DO igc = 1,ngc(16)
          sumk = 0.0_wp
          DO ipr = 1, ngn(ngs(15)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+240)
          ENDDO
          kb(jt,jp,igc) = sumk
        ENDDO
      ENDDO
    ENDDO

    DO jt = 1,10
      iprsm = 0
      DO igc = 1,ngc(16)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(15)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+240)
        ENDDO
        selfref(jt,igc) = sumk
      ENDDO
    ENDDO

    DO jt = 1,4
      iprsm = 0
      DO igc = 1,ngc(16)
        sumk = 0.0_wp
        DO ipr = 1, ngn(ngs(15)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+240)
        ENDDO
        forref(jt,igc) = sumk
      ENDDO
    ENDDO

    iprsm = 0
    DO igc = 1,ngc(16)
      sumf = 0.0_wp
      DO ipr = 1, ngn(ngs(15)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefbo(iprsm)
      ENDDO
      fracrefb(igc) = sumf
    ENDDO

    DO jp = 1,9
      iprsm = 0
      DO igc = 1,ngc(16)
        sumf = 0.0_wp
        DO ipr = 1, ngn(ngs(15)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        ENDDO
        fracrefa(igc,jp) = sumf
      ENDDO
    ENDDO

  END SUBROUTINE cmbgb16
  !***************************************************************************

END MODULE mo_lrtm_setup
