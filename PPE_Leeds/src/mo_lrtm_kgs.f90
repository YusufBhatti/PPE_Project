!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE rrlw_planck
  USE mo_kind,        ONLY : wp 
  USE mo_rrtm_params, ONLY : nbndlw
  
  REAL(wp) :: chi_mls(7,59)
  REAL(wp) :: totplanck(181,nbndlw) !< planck function for each band
  REAL(wp) :: totplanck16(181)      !< for band 16

END MODULE rrlw_planck

!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 1 Parameters: 10-250 cm-1 (low - h2o; high - h2o)
!
MODULE rrlw_kg01

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no1  = 16 !< original abs coefficients

  REAL(wp) :: fracrefao(no1)  , fracrefbo(no1)
  REAL(wp) :: kao(5,13,no1)
  REAL(wp) :: kbo(5,13:59,no1)
  REAL(wp) :: kao_mn2(19,no1) , kbo_mn2(19,no1)
  REAL(wp) :: selfrefo(10,no1), forrefo(4,no1)

  INTEGER, PARAMETER :: ng1  = 10 !< combined abs. coefficients

  REAL(wp) :: fracrefa(ng1)  , fracrefb(ng1)
  REAL(wp) :: ka(5,13,ng1)   , absa(65,ng1)
  REAL(wp) :: kb(5,13:59,ng1), absb(235,ng1)
  REAL(wp) :: ka_mn2(19,ng1) , kb_mn2(19,ng1)
  REAL(wp) :: selfref(10,ng1), forref(4,ng1)

  EQUIVALENCE (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrlw_kg01
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 2 Parameters 250-500 cm-1 (low - h2o; high - h2o)
!
MODULE rrlw_kg02

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no2  = 16

  REAL(wp) :: fracrefao(no2)   , fracrefbo(no2)
  REAL(wp) :: kao(5,13,no2)
  REAL(wp) :: kbo(5,13:59,no2)
  REAL(wp) :: selfrefo(10,no2) , forrefo(4,no2)

  INTEGER, PARAMETER :: ng2  = 12

  REAL(wp) :: fracrefa(ng2)  , fracrefb(ng2)
  REAL(wp) :: ka(5,13,ng2)   , absa(65,ng2)
  REAL(wp) :: kb(5,13:59,ng2), absb(235,ng2)
  REAL(wp) :: selfref(10,ng2), forref(4,ng2)
  REAL(wp) :: refparam(13)

  EQUIVALENCE (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE rrlw_kg02
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 3 Parameters 500-630 cm-1 (low - h2o,co2; high - h2o,co2)
!
MODULE rrlw_kg03

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no3  = 16

  REAL(wp) :: fracrefao(no3,9) ,fracrefbo(no3,5)
  REAL(wp) :: kao(9,5,13,no3)
  REAL(wp) :: kbo(5,5,13:59,no3)
  REAL(wp) :: kao_mn2o(9,19,no3), kbo_mn2o(5,19,no3)
  REAL(wp) :: selfrefo(10,no3)
  REAL(wp) :: forrefo(4,no3)

  INTEGER, PARAMETER :: ng3  = 16

  REAL(wp) :: fracrefa(ng3,9) ,fracrefb(ng3,5)
  REAL(wp) :: ka(9,5,13,ng3)  ,absa(585,ng3)
  REAL(wp) :: kb(5,5,13:59,ng3),absb(1175,ng3)
  REAL(wp) :: ka_mn2o(9,19,ng3), kb_mn2o(5,19,ng3)
  REAL(wp) :: selfref(10,ng3)
  REAL(wp) :: forref(4,ng3)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

END MODULE rrlw_kg03
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 4 Parameters 630-700 cm-1 (low - h2o,co2; high - o3,co2)
!
MODULE rrlw_kg04

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no4  = 16

  REAL(wp) :: fracrefao(no4,9)  ,fracrefbo(no4,5)
  REAL(wp) :: kao(9,5,13,no4)
  REAL(wp) :: kbo(5,5,13:59,no4)
  REAL(wp) :: selfrefo(10,no4)  ,forrefo(4,no4)

  INTEGER, PARAMETER :: ng4  = 14

  REAL(wp) :: fracrefa(ng4,9)  ,fracrefb(ng4,5)
  REAL(wp) :: ka(9,5,13,ng4)   ,absa(585,ng4)
  REAL(wp) :: kb(5,5,13:59,ng4),absb(1175,ng4)
  REAL(wp) :: selfref(10,ng4)  ,forref(4,ng4)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

END MODULE rrlw_kg04
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 5 Parameters 700-820 cm-1 (low - h2o,co2; high - o3,co2)
!
MODULE rrlw_kg05

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no5  = 16

  REAL(wp) :: fracrefao(no5,9) ,fracrefbo(no5,5)
  REAL(wp) :: kao(9,5,13,no5)
  REAL(wp) :: kbo(5,5,13:59,no5)
  REAL(wp) :: kao_mo3(9,19,no5)
  REAL(wp) :: selfrefo(10,no5)
  REAL(wp) :: forrefo(4,no5)
  REAL(wp) :: ccl4o(no5)

  INTEGER, PARAMETER :: ng5  = 16

  REAL(wp) :: fracrefa(ng5,9) ,fracrefb(ng5,5)
  REAL(wp) :: ka(9,5,13,ng5)   ,absa(585,ng5)
  REAL(wp) :: kb(5,5,13:59,ng5),absb(1175,ng5)
  REAL(wp) :: ka_mo3(9,19,ng5)
  REAL(wp) :: selfref(10,ng5)
  REAL(wp) :: forref(4,ng5)
  REAL(wp) :: ccl4(ng5)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

END MODULE rrlw_kg05
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 6 Parameters 820-980 cm-1 (low - h2o; high - nothing)
!
MODULE rrlw_kg06

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no6  = 16

  REAL(wp) , DIMENSION(no6) :: fracrefao
  REAL(wp) :: kao(5,13,no6)
  REAL(wp) :: kao_mco2(19,no6)
  REAL(wp) :: selfrefo(10,no6)
  REAL(wp) :: forrefo(4,no6)

  REAL(wp) , DIMENSION(no6) :: cfc11adjo
  REAL(wp) , DIMENSION(no6) :: cfc12o

  INTEGER, PARAMETER :: ng6  = 8

  REAL(wp) , DIMENSION(ng6) :: fracrefa
  REAL(wp) :: ka(5,13,ng6),absa(65,ng6)
  REAL(wp) :: ka_mco2(19,ng6)
  REAL(wp) :: selfref(10,ng6)
  REAL(wp) :: forref(4,ng6)

  REAL(wp) , DIMENSION(ng6) :: cfc11adj
  REAL(wp) , DIMENSION(ng6) :: cfc12

  EQUIVALENCE (ka(1,1,1),absa(1,1))

END MODULE rrlw_kg06
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 7 Parameters 980-1080 cm-1 (low - h2o,o3; high - o3)
!
MODULE rrlw_kg07

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no7  = 16

  REAL(wp) , DIMENSION(no7) :: fracrefbo
  REAL(wp) :: fracrefao(no7,9)
  REAL(wp) :: kao(9,5,13,no7)
  REAL(wp) :: kbo(5,13:59,no7)
  REAL(wp) :: kao_mco2(9,19,no7)
  REAL(wp) :: kbo_mco2(19,no7)
  REAL(wp) :: selfrefo(10,no7)
  REAL(wp) :: forrefo(4,no7)

  INTEGER, PARAMETER :: ng7  = 12

  REAL(wp) , DIMENSION(ng7) :: fracrefb
  REAL(wp) :: fracrefa(ng7,9)
  REAL(wp) :: ka(9,5,13,ng7) ,absa(585,ng7)
  REAL(wp) :: kb(5,13:59,ng7),absb(235,ng7)
  REAL(wp) :: ka_mco2(9,19,ng7)
  REAL(wp) :: kb_mco2(19,ng7)
  REAL(wp) :: selfref(10,ng7)
  REAL(wp) :: forref(4,ng7)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE rrlw_kg07
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 8 Parameters 1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!
MODULE rrlw_kg08

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no8  = 16

  REAL(wp) , DIMENSION(no8) :: fracrefao
  REAL(wp) , DIMENSION(no8) :: fracrefbo
  REAL(wp) , DIMENSION(no8) :: cfc12o
  REAL(wp) , DIMENSION(no8) :: cfc22adjo

  REAL(wp) :: kao(5,13,no8)
  REAL(wp) :: kao_mco2(19,no8)
  REAL(wp) :: kao_mn2o(19,no8)
  REAL(wp) :: kao_mo3(19,no8)
  REAL(wp) :: kbo(5,13:59,no8)
  REAL(wp) :: kbo_mco2(19,no8)
  REAL(wp) :: kbo_mn2o(19,no8)
  REAL(wp) :: selfrefo(10,no8)
  REAL(wp) :: forrefo(4,no8)

  INTEGER, PARAMETER :: ng8  = 8

  REAL(wp) , DIMENSION(ng8) :: fracrefa
  REAL(wp) , DIMENSION(ng8) :: fracrefb
  REAL(wp) , DIMENSION(ng8) :: cfc12
  REAL(wp) , DIMENSION(ng8) :: cfc22adj

  REAL(wp) :: ka(5,13,ng8)    ,absa(65,ng8)
  REAL(wp) :: kb(5,13:59,ng8) ,absb(235,ng8)
  REAL(wp) :: ka_mco2(19,ng8)
  REAL(wp) :: ka_mn2o(19,ng8)
  REAL(wp) :: ka_mo3(19,ng8)
  REAL(wp) :: kb_mco2(19,ng8)
  REAL(wp) :: kb_mn2o(19,ng8)
  REAL(wp) :: selfref(10,ng8)
  REAL(wp) :: forref(4,ng8)

  EQUIVALENCE (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE rrlw_kg08
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 9 parameters 1180-1390 cm-1 (low - h2o,ch4; high - ch4)
!
MODULE rrlw_kg09

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no9  = 16

  REAL(wp) , DIMENSION(no9) :: fracrefbo

  REAL(wp) :: fracrefao(no9,9)
  REAL(wp) :: kao(9,5,13,no9)
  REAL(wp) :: kbo(5,13:59,no9)
  REAL(wp) :: kao_mn2o(9,19,no9)
  REAL(wp) :: kbo_mn2o(19,no9)
  REAL(wp) :: selfrefo(10,no9)
  REAL(wp) :: forrefo(4,no9)

  INTEGER, PARAMETER :: ng9  = 12

  REAL(wp) , DIMENSION(ng9) :: fracrefb
  REAL(wp) :: fracrefa(ng9,9)
  REAL(wp) :: ka(9,5,13,ng9) ,absa(585,ng9)
  REAL(wp) :: kb(5,13:59,ng9) ,absb(235,ng9)
  REAL(wp) :: ka_mn2o(9,19,ng9)
  REAL(wp) :: kb_mn2o(19,ng9)
  REAL(wp) :: selfref(10,ng9)
  REAL(wp) :: forref(4,ng9)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE rrlw_kg09
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 10 Parameters 1390-1480 cm-1 (low - h2o; high - h2o)
!
MODULE rrlw_kg10

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no10 = 16

  REAL(wp) , DIMENSION(no10) :: fracrefao
  REAL(wp) , DIMENSION(no10) :: fracrefbo

  REAL(wp) :: kao(5,13,no10)
  REAL(wp) :: kbo(5,13:59,no10)
  REAL(wp) :: selfrefo(10,no10)
  REAL(wp) :: forrefo(4,no10)

  INTEGER, PARAMETER :: ng10 = 6

  REAL(wp) , DIMENSION(ng10) :: fracrefa
  REAL(wp) , DIMENSION(ng10) :: fracrefb

  REAL(wp) :: ka(5,13,ng10)   , absa(65,ng10)
  REAL(wp) :: kb(5,13:59,ng10), absb(235,ng10)
  REAL(wp) :: selfref(10,ng10)
  REAL(wp) :: forref(4,ng10)

  EQUIVALENCE (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE rrlw_kg10
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 11 Parameters 1480-1800 cm-1 (low - h2o; high - h2o)
!
MODULE rrlw_kg11

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no11 = 16

  REAL(wp) , DIMENSION(no11) :: fracrefao
  REAL(wp) , DIMENSION(no11) :: fracrefbo

  REAL(wp) :: kao(5,13,no11)
  REAL(wp) :: kbo(5,13:59,no11)
  REAL(wp) :: kao_mo2(19,no11)
  REAL(wp) :: kbo_mo2(19,no11)
  REAL(wp) :: selfrefo(10,no11)
  REAL(wp) :: forrefo(4,no11)

  INTEGER, PARAMETER :: ng11 = 8

  REAL(wp) , DIMENSION(ng11) :: fracrefa
  REAL(wp) , DIMENSION(ng11) :: fracrefb

  REAL(wp) :: ka(5,13,ng11)   , absa(65,ng11)
  REAL(wp) :: kb(5,13:59,ng11), absb(235,ng11)
  REAL(wp) :: ka_mo2(19,ng11)
  REAL(wp) :: kb_mo2(19,ng11)
  REAL(wp) :: selfref(10,ng11)
  REAL(wp) :: forref(4,ng11)

  EQUIVALENCE (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE rrlw_kg11
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 12 Parameters 1800-2080 cm-1 (low - h2o,co2; high - nothing)
!
MODULE rrlw_kg12

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no12 = 16

  REAL(wp) :: fracrefao(no12,9)
  REAL(wp) :: kao(9,5,13,no12)
  REAL(wp) :: selfrefo(10,no12)
  REAL(wp) :: forrefo(4,no12)

  INTEGER, PARAMETER :: ng12 = 8

  REAL(wp) :: fracrefa(ng12,9)
  REAL(wp) :: ka(9,5,13,ng12) ,absa(585,ng12)
  REAL(wp) :: selfref(10,ng12)
  REAL(wp) :: forref(4,ng12)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1))

END MODULE rrlw_kg12
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 13 Parameters 2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!
MODULE rrlw_kg13

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no13 = 16

  REAL(wp) , DIMENSION(no13) :: fracrefbo

  REAL(wp) :: fracrefao(no13,9)
  REAL(wp) :: kao(9,5,13,no13)
  REAL(wp) :: kao_mco2(9,19,no13)
  REAL(wp) :: kao_mco(9,19,no13)
  REAL(wp) :: kbo_mo3(19,no13)
  REAL(wp) :: selfrefo(10,no13)
  REAL(wp) :: forrefo(4,no13)

  INTEGER, PARAMETER :: ng13 = 4

  REAL(wp) , DIMENSION(ng13) :: fracrefb

  REAL(wp) :: fracrefa(ng13,9)
  REAL(wp) :: ka(9,5,13,ng13) ,absa(585,ng13)
  REAL(wp) :: ka_mco2(9,19,ng13)
  REAL(wp) :: ka_mco(9,19,ng13)
  REAL(wp) :: kb_mo3(19,ng13)
  REAL(wp) :: selfref(10,ng13)
  REAL(wp) :: forref(4,ng13)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1))

END MODULE rrlw_kg13
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 14 Parameters 2250-2380 cm-1 (low - co2; high - co2)
!
MODULE rrlw_kg14

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no14 = 16

  REAL(wp) , DIMENSION(no14) :: fracrefao
  REAL(wp) , DIMENSION(no14) :: fracrefbo

  REAL(wp) :: kao(5,13,no14)
  REAL(wp) :: kbo(5,13:59,no14)
  REAL(wp) :: selfrefo(10,no14)
  REAL(wp) :: forrefo(4,no14)

  INTEGER, PARAMETER :: ng14 = 2

  REAL(wp) , DIMENSION(ng14) :: fracrefa
  REAL(wp) , DIMENSION(ng14) :: fracrefb

  REAL(wp) :: ka(5,13,ng14)   ,absa(65,ng14)
  REAL(wp) :: kb(5,13:59,ng14),absb(235,ng14)
  REAL(wp) :: selfref(10,ng14)
  REAL(wp) :: forref(4,ng14)

  EQUIVALENCE (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrlw_kg14
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 15 Parameters 2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
MODULE rrlw_kg15

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no15 = 16

  REAL(wp) :: fracrefao(no15,9)
  REAL(wp) :: kao(9,5,13,no15)
  REAL(wp) :: kao_mn2(9,19,no15)
  REAL(wp) :: selfrefo(10,no15)
  REAL(wp) :: forrefo(4,no15)

  INTEGER, PARAMETER :: ng15 = 2

  REAL(wp) :: fracrefa(ng15,9)
  REAL(wp) :: ka(9,5,13,ng15) ,absa(585,ng15)
  REAL(wp) :: ka_mn2(9,19,ng15)
  REAL(wp) :: selfref(10,ng15)
  REAL(wp) :: forref(4,ng15)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1))

END MODULE rrlw_kg15
!-----------------------------------------------------------------------------
!>
!! @brief LRTM Band 16 Parameters 2600-3000 cm-1 (low - h2o,ch4; high - nothing)
!
MODULE rrlw_kg16

  USE mo_kind ,ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no16 = 16

  REAL(wp) , DIMENSION(no16) :: fracrefbo

  REAL(wp) :: fracrefao(no16,9)
  REAL(wp) :: kao(9,5,13,no16)
  REAL(wp) :: kbo(5,13:59,no16)
  REAL(wp) :: selfrefo(10,no16)
  REAL(wp) :: forrefo(4,no16)

  INTEGER, PARAMETER :: ng16 = 2

  REAL(wp) , DIMENSION(ng16) :: fracrefb

  REAL(wp) :: fracrefa(ng16,9)
  REAL(wp) :: ka(9,5,13,ng16) ,absa(585,ng16)
  REAL(wp) :: kb(5,13:59,ng16), absb(235,ng16)
  REAL(wp) :: selfref(10,ng16)
  REAL(wp) :: forref(4,ng16)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrlw_kg16
