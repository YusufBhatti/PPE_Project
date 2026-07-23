!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE rrsw_kg16

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no16 = 16
  REAL(wp) :: kao(9,5,13,no16)
  REAL(wp) :: kbo(5,13:59,no16)
  REAL(wp) :: selfrefo(10,no16), forrefo(3,no16)
  REAL(wp) :: sfluxrefo(no16)
  REAL(wp) :: rayl

  INTEGER, PARAMETER :: ng16 = 6
  REAL(wp) :: ka(9,5,13,ng16) , absa(585,ng16)
  REAL(wp) :: kb(5,13:59,ng16), absb(235,ng16)
  REAL(wp) :: selfref(10,ng16), forref(3,ng16)
  REAL(wp) :: sfluxref(ng16)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrsw_kg16

MODULE rrsw_kg17

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no17 = 16
  REAL(wp) :: kao(9,5,13,no17)
  REAL(wp) :: kbo(5,5,13:59,no17)
  REAL(wp) :: selfrefo(10,no17), forrefo(4,no17)
  REAL(wp) :: sfluxrefo(no17,5)
  REAL(wp) :: rayl

  INTEGER, PARAMETER :: ng17 = 12
  REAL(wp) :: ka(9,5,13,ng17) , absa(585,ng17)
  REAL(wp) :: kb(5,5,13:59,ng17), absb(1175,ng17)
  REAL(wp) :: selfref(10,ng17), forref(4,ng17)
  REAL(wp) :: sfluxref(ng17,5)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

END MODULE rrsw_kg17

MODULE rrsw_kg18

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no18 = 16
  REAL(wp) :: kao(9,5,13,no18)
  REAL(wp) :: kbo(5,13:59,no18)
  REAL(wp) :: selfrefo(10,no18), forrefo(3,no18)
  REAL(wp) :: sfluxrefo(no18,9)
  REAL(wp) :: rayl

  INTEGER, PARAMETER :: ng18 = 8
  REAL(wp) :: ka(9,5,13,ng18), absa(585,ng18)
  REAL(wp) :: kb(5,13:59,ng18), absb(235,ng18)
  REAL(wp) :: selfref(10,ng18), forref(3,ng18)
  REAL(wp) :: sfluxref(ng18,9)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrsw_kg18

MODULE rrsw_kg19

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no19 = 16
  REAL(wp) :: kao(9,5,13,no19)
  REAL(wp) :: kbo(5,13:59,no19)
  REAL(wp) :: selfrefo(10,no19), forrefo(3,no19)
  REAL(wp) :: sfluxrefo(no19,9)
  REAL(wp) :: rayl

  INTEGER, PARAMETER :: ng19 = 8
  REAL(wp) :: ka(9,5,13,ng19), absa(585,ng19)
  REAL(wp) :: kb(5,13:59,ng19), absb(235,ng19)
  REAL(wp) :: selfref(10,ng19), forref(3,ng19)
  REAL(wp) :: sfluxref(ng19,9)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrsw_kg19

MODULE rrsw_kg20

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no20 = 16

  REAL(wp) :: kao(5,13,no20)
  REAL(wp) :: kbo(5,13:59,no20)
  REAL(wp) :: selfrefo(10,no20), forrefo(4,no20)
  REAL(wp) :: sfluxrefo(no20)
  REAL(wp) :: absch4o(no20)

  REAL(wp) :: rayl 

  INTEGER, PARAMETER :: ng20 = 10
  REAL(wp) :: ka(5,13,ng20), absa(65,ng20)
  REAL(wp) :: kb(5,13:59,ng20), absb(235,ng20)
  REAL(wp) :: selfref(10,ng20), forref(4,ng20)
  REAL(wp) :: sfluxref(ng20)
  REAL(wp) :: absch4(ng20)

  EQUIVALENCE (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrsw_kg20

MODULE rrsw_kg21

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no21 = 16
  REAL(wp) :: kao(9,5,13,no21)
  REAL(wp) :: kbo(5,5,13:59,no21)
  REAL(wp) :: selfrefo(10,no21), forrefo(4,no21)
  REAL(wp) :: sfluxrefo(no21,9)
  REAL(wp) :: rayl

  INTEGER, PARAMETER :: ng21 = 10
  REAL(wp) :: ka(9,5,13,ng21), absa(585,ng21)
  REAL(wp) :: kb(5,5,13:59,ng21), absb(1175,ng21)
  REAL(wp) :: selfref(10,ng21), forref(4,ng21)
  REAL(wp) :: sfluxref(ng21,9)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

END MODULE rrsw_kg21

MODULE rrsw_kg22

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no22 = 16

  REAL(wp) :: kao(9,5,13,no22)
  REAL(wp) :: kbo(5,13:59,no22)
  REAL(wp) :: selfrefo(10,no22), forrefo(3,no22)
  REAL(wp) :: sfluxrefo(no22,9)

  REAL(wp) :: rayl

  INTEGER, PARAMETER :: ng22 = 2
  REAL(wp) :: ka(9,5,13,ng22), absa(585,ng22)
  REAL(wp) :: kb(5,13:59,ng22), absb(235,ng22)
  REAL(wp) :: selfref(10,ng22), forref(3,ng22)
  REAL(wp) :: sfluxref(ng22,9)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrsw_kg22

MODULE rrsw_kg23

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no23 = 16
  REAL(wp) :: kao(5,13,no23)
  REAL(wp) :: selfrefo(10,no23), forrefo(3,no23)
  REAL(wp) :: sfluxrefo(no23)
  REAL(wp) :: raylo(no23)

  INTEGER, PARAMETER :: ng23 = 10
  REAL(wp) :: ka(5,13,ng23), absa(65,ng23)
  REAL(wp) :: selfref(10,ng23), forref(3,ng23)
  REAL(wp) :: sfluxref(ng23), rayl(ng23)

  EQUIVALENCE (ka(1,1,1),absa(1,1))

END MODULE rrsw_kg23

MODULE rrsw_kg24

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no24 = 16
  REAL(wp) :: kao(9,5,13,no24)
  REAL(wp) :: kbo(5,13:59,no24)
  REAL(wp) :: selfrefo(10,no24), forrefo(3,no24)
  REAL(wp) :: sfluxrefo(no24,9)
  REAL(wp) :: abso3ao(no24), abso3bo(no24)
  REAL(wp) :: raylao(no24,9), raylbo(no24)

  INTEGER, PARAMETER :: ng24 = 8
  REAL(wp) :: ka(9,5,13,ng24), absa(585,ng24)
  REAL(wp) :: kb(5,13:59,ng24), absb(235,ng24)
  REAL(wp) :: selfref(10,ng24), forref(3,ng24)
  REAL(wp) :: sfluxref(ng24,9)
  REAL(wp) :: abso3a(ng24), abso3b(ng24)
  REAL(wp) :: rayla(ng24,9), raylb(ng24)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrsw_kg24

MODULE rrsw_kg25

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no25 = 16
  REAL(wp) :: kao(5,13,no25)
  REAL(wp) :: sfluxrefo(no25)
  REAL(wp) :: abso3ao(no25), abso3bo(no25)
  REAL(wp) :: raylo(no25)

  INTEGER, PARAMETER :: ng25 = 6
  REAL(wp) :: ka(5,13,ng25), absa(65,ng25)
  REAL(wp) :: sfluxref(ng25)
  REAL(wp) :: abso3a(ng25), abso3b(ng25)
  REAL(wp) :: rayl(ng25)

  EQUIVALENCE (ka(1,1,1),absa(1,1))

END MODULE rrsw_kg25

MODULE rrsw_kg26

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no26 = 16
  REAL(wp) :: sfluxrefo(no26)
  REAL(wp) :: raylo(no26)

  INTEGER, PARAMETER :: ng26 = 6
  REAL(wp) :: sfluxref(ng26)
  REAL(wp) :: rayl(ng26)

END MODULE rrsw_kg26

MODULE rrsw_kg27

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no27 = 16
  REAL(wp) :: kao(5,13,no27)
  REAL(wp) :: kbo(5,13:59,no27)
  REAL(wp) :: sfluxrefo(no27)
  REAL(wp) :: raylo(no27)

  INTEGER, PARAMETER :: ng27 = 8
  REAL(wp) :: ka(5,13,ng27), absa(65,ng27)
  REAL(wp) :: kb(5,13:59,ng27), absb(235,ng27)
  REAL(wp) :: sfluxref(ng27)
  REAL(wp) :: rayl(ng27)

  EQUIVALENCE (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrsw_kg27

MODULE rrsw_kg28

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no28 = 16
  REAL(wp) :: kao(9,5,13,no28)
  REAL(wp) :: kbo(5,5,13:59,no28)
  REAL(wp) :: sfluxrefo(no28,5)
  REAL(wp) :: rayl

  INTEGER, PARAMETER :: ng28 = 6
  REAL(wp) :: ka(9,5,13,ng28), absa(585,ng28)
  REAL(wp) :: kb(5,5,13:59,ng28), absb(1175,ng28)
  REAL(wp) :: sfluxref(ng28,5)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

END MODULE rrsw_kg28

MODULE rrsw_kg29

  USE mo_kind, ONLY : wp
  IMPLICIT NONE

  INTEGER, PARAMETER :: no29 = 16
  REAL(wp) :: kao(5,13,no29)
  REAL(wp) :: kbo(5,13:59,no29)
  REAL(wp) :: selfrefo(10,no29), forrefo(4,no29)
  REAL(wp) :: sfluxrefo(no29)
  REAL(wp) :: absh2oo(no29), absco2o(no29)
  REAL(wp) :: rayl

  INTEGER, PARAMETER :: ng29 = 12
  REAL(wp) :: ka(5,13,ng29), absa(65,ng29)
  REAL(wp) :: kb(5,13:59,ng29), absb(235,ng29)
  REAL(wp) :: selfref(10,ng29), forref(4,ng29)
  REAL(wp) :: sfluxref(ng29)
  REAL(wp) :: absh2o(ng29), absco2(ng29)

  EQUIVALENCE (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE rrsw_kg29

