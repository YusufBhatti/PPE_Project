!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_methox

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE

  REAL(KIND=wp), PARAMETER :: RQLIM=4.25E-6_wp
  REAL(KIND=wp), PARAMETER :: RPBOTOX=10000._wp
  REAL(KIND=wp), PARAMETER :: RPBOTPH=20._wp
  REAL(KIND=wp), PARAMETER :: RPTOPOX=50._wp
  REAL(KIND=wp), PARAMETER :: RPTOPPH=0.1_wp

  REAL(KIND=wp) :: RALPHA1
  REAL(KIND=wp) :: RALPHA2
  REAL(KIND=wp) :: RALPHA3
  REAL(KIND=wp) :: RLOGPPH

  CONTAINS

  SUBROUTINE init_methox

!**** *init_methox*   - INITIALIZE MODULE mo_methox CONTROLLING *METHOX*

!     PURPOSE.
!     --------
!           INITIALIZE methox variables

!**   INTERFACE.
!     ----------
!        CALL *init_methox* FROM *init_subm*
!              --------        ------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "INTEGRATED FORECASTING SYSTEM"

!     AUTHOR.
!     -------
!        C.JAKOB   *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 98-04-07
!        Modified : 02-01-29  A.Simmons: increase RQLIM from 3.75e-6
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.G. Fogli & C. Cagnazzo (CMCC): July 2010, ECHAM5 implementation
!        H. Schmidt: Aug 2010, ECHAM6 implementation
!     ------------------------------------------------------------------


   IMPLICIT NONE

    RALPHA1 = (19._wp*LOG(10._wp))/(LOG(20._wp)**4)
    RALPHA2 = LOG(1.0_wp/3._wp+0.01_wp)
    RALPHA3 = 0.5_wp*(LOG(100._wp)+RALPHA2)
    RLOGPPH = LOG(RPTOPPH/RPBOTPH)

  END SUBROUTINE init_methox






  SUBROUTINE methox(kproma, kbdim, klev, pq, ptenq, pap)

!**** *METHOX*   - Calculate humidity tendencies from methane
!                  oxidation and photolysis

!**   INTERFACE.
!     ----------
!        CALL *METHOX* FROM *physc*
!              ------        -------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KPROMA*       geographic block number of locations
!    *KBDIM*        geographic block maximum number of locations
!    *KLEV*         NUMBER OF LEVELS

!     INPUT PARAMETERS (REAL):

!    *PAP*          PRESSURE                                      PA
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG

!     UPDATED PARAMETERS (REAL):

!    *PTENQ*        TENDENCY OF SPECIFIC HUMIDITY                 KG/(KG*S)

!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        MODULE mo_kind
!        MODULE mo_math_constants

!     PURPOSE.
!     --------
!       In GCMs without chemistry coupling, the water vapor concentration in the startosphere and mesosphere
!       is unrealistic because the major production (oxidation of methane) and loss (photolysis) mechanisms are missing.
!       This submodel parameterizes both effects and should lead to a vertical profile of water vapor mixing ratio with 
!       a secondary maximum near the stratopause. 

!     METHOD.
!     -------
!        SEE RD-MEMO R60.1/AJS/31

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        SEE RD-MEMO R60.1/AJS/31

!     AUTHOR.
!     -------
!        C.JAKOB   *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 98-04-07
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.G. Fogli & C. Cagnazzo (CMCC): July 2010, ECHAM5 implementation
!        H. Schmidt: Aug 2010, ECHAM6 implementation
!     ------------------------------------------------------------------

    USE mo_kind,   ONLY : wp

    USE mo_math_constants   , ONLY :  pi

    IMPLICIT NONE

    INTEGER,INTENT(IN)    :: KBDIM
    INTEGER,INTENT(IN)    :: KPROMA
    INTEGER,INTENT(IN)    :: KLEV
    REAL(KIND=wp)   ,INTENT(IN)    :: PQ(KBDIM,KLEV)
    REAL(KIND=wp)   ,INTENT(INOUT) :: PTENQ(KBDIM,KLEV)
    REAL(KIND=wp)   ,INTENT(IN)    :: PAP(KBDIM,KLEV)

    LOGICAL :: LLOXID,         LLPHOTO

    INTEGER :: JK, JL

    REAL(KIND=wp) :: ZARG, ZPRATIO, ZTAU1, ZTAU2, ZTDAYS

    DO JK=1,KLEV
      DO JL=1,KPROMA        

        LLOXID=PAP(JL,JK) < RPBOTOX.AND.PQ(JL,JK) < RQLIM
        LLPHOTO=PAP(JL,JK) < RPBOTPH

!       METHANE OXIDATION

        IF(LLOXID) THEN
          IF(PAP(JL,JK) <= RPTOPOX) THEN
            ZTDAYS=100._wp
          ELSE
            ZPRATIO=(LOG(PAP(JL,JK)/RPTOPOX))**4._wp/LOG(RPBOTOX/PAP(JL,JK))
            ZTDAYS=100._wp*(1._wp+RALPHA1*ZPRATIO)
          ENDIF
          ZTAU1=86400._wp*ZTDAYS
          PTENQ(JL,JK)=PTENQ(JL,JK)+(RQLIM-PQ(JL,JK))/ZTAU1
        ENDIF


!       PHOTOLYSIS

        IF(LLPHOTO) THEN
          IF(PAP(JL,JK) <= RPTOPPH) THEN
            ZTDAYS=3._wp
          ELSE
            ZARG=RALPHA2-RALPHA3*(1._wp+COS((PI*LOG(PAP(JL,JK)/RPBOTPH))/RLOGPPH))
            ZTDAYS=1.0_wp/(EXP(ZARG)-0.01_wp)
          ENDIF
          ZTAU2=86400._wp*ZTDAYS
          PTENQ(JL,JK)=PTENQ(JL,JK)-PQ(JL,JK)/ZTAU2
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE methox


END MODULE mo_methox


