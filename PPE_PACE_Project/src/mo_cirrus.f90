MODULE mo_cirrus

  ! Author:  Axel Lauer, DLR-IPA Oberpfaffenhofen
  ! last modified: 01/2006
  ! Modified by Ulrike Lohmann, ETH Zurich for use with ECHAM5-HAM

  USE mo_kind, ONLY: dp

  INTEGER :: nfrzmod

  CONTAINS

! ------------------------------------------------------------------------------

  REAL(dp) FUNCTION SCRHOM(T)

     ! ***** CRITICAL ICE SATURATION RATIO FOR HOMOGENEOUS FREEZING
     ! ***** EVALUATED FOR AN INTERMEDIATE PARTICLE SIZE OF 0.25 MUM
     ! ***** TH.KOOP, B.P.LUO, A.TSIAS, TH.PETER, NATURE 406, 611-614, 2000

     IMPLICIT NONE

     ! function parameters

     REAL(dp), INTENT(in) :: T

     ! local variables

     REAL(dp) :: TEMP

     TEMP = MIN(240.0_dp,MAX(170.0_dp,T))
     SCRHOM = 2.418_dp - (TEMP/245.68_dp)

     RETURN

  END FUNCTION SCRHOM

! ------------------------------------------------------------------------------

  REAL(dp) FUNCTION PISAT(T)

     ! ***** VAPOR PRESSURE OVER ICE IN MBAR (T IN K)
     ! ***** J.MARTI AND K.MAUERSBERGER, GRL 20(5), 363-366, 1993

     IMPLICIT NONE

     ! function parameters

     REAL(dp), INTENT(in) :: T

     ! parameters

     REAL(dp), PARAMETER :: A = 0.01_dp
     REAL(dp), PARAMETER :: B = 12.537_dp
     REAL(dp), PARAMETER :: C = -2663.5_dp

     PISAT = A * 10.0_dp**( B + C/T )

     RETURN

  END FUNCTION PISAT

! ------------------------------------------------------------------------------

  REAL(dp) FUNCTION XEERFC(Y)

     ! ***** PRODUCT EXP(1/Y) * ERFC(1/SQRT(Y)) AND ASYMPTOTES

     IMPLICIT NONE

     ! function parameters

     REAL(dp), INTENT(in) :: Y

     ! parameters

     REAL(dp), PARAMETER :: SQPI = 1.7724539_dp

     ! local variables

     REAL(dp) :: Y1,SQY1,PPROD,POLY
     INTEGER :: K

!     ! functions
!
!     REAL(dp) :: XERF ! error function with rational approximation

     Y1   = 1.0_dp / Y
     SQY1 = SQRT(Y1)
     IF (Y.LE.0.2_dp) THEN
        PPROD = 1.0_dp
        POLY  = 1.0_dp
        DO K = 1, 4
           PPROD = PPROD * REAL(2*K-1,dp)
!           POLY  = POLY  + PPROD * (-0.5_dp*Y)**REAL(K,dp)
           POLY  = POLY  + PPROD * (-0.5_dp*Y)**K
        END DO
        XEERFC = POLY / SQY1 / SQPI
     ELSEIF (Y.GE.2.0_dp) THEN
        PPROD  = 1.0_dp
        POLY   = 1.0_dp
        DO K = 1, 4
           PPROD = PPROD * REAL(2*K+1,dp)
           POLY  = POLY  + REAL(2**K,dp) * Y1**K / PPROD
        END DO
        XEERFC = EXP(Y1) - 2.0_dp*SQY1/SQPI * POLY
     ELSE
        XEERFC = ( 1.0_dp - XERF(SQY1) ) * EXP(Y1)
     END IF

     RETURN

  END FUNCTION XEERFC

! ------------------------------------------------------------------------------

  REAL(dp) FUNCTION XERF(X)

     ! ***** ERROR FUNCTION WITH RATIONAL APPROXIMATION
     ! ***** M.ABRAMOWITZ, I.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS
     ! ***** DOVER, 10th PRINTING 1972, p.299, # 7.1.26

     IMPLICIT NONE

     ! function parameters

     REAL(dp), INTENT(in) :: X

     ! parameters

     REAL(dp), PARAMETER :: P  = 0.3275911_dp
     REAL(dp), PARAMETER :: A1 = 0.254829592_dp
     REAL(dp), PARAMETER :: A2 =-0.284496736_dp
     REAL(dp), PARAMETER :: A3 = 1.421413741_dp
     REAL(dp), PARAMETER :: A4 =-1.453152027_dp
     REAL(dp), PARAMETER :: A5 = 1.061405429_dp

     ! local variables

     REAL(dp) :: T,ARG

     T    = 1.0_dp / ( 1.0_dp + P * ABS(X) )
     ARG  = MIN( X**2, 75.0_dp )
     XERF = 1.0_dp - ( A1*T + A2*T**2 + A3*T**3 + A4*T**4 + A5*T**5 ) &
            * EXP( -ARG )
     XERF = SIGN( XERF, X )

     RETURN

     END FUNCTION XERF

! ------------------------------------------------------------------------------

 REAL(dp) FUNCTION SCRHET(T, NFRZ)

     IMPLICIT NONE

     ! function parameters

     REAL(dp), INTENT(in) :: T
     INTEGER,  INTENT(in) :: NFRZ

     ! ***** CRITICAL ICE SATURATION RATIO FOR HETEROGENEOUS FREEZING

     IF (NFRZ.EQ.0) THEN

        ! ***** CONSTANT THRESHOLD (FREEZING OF 1 PTCL PER SEC, RADIUS 0.25 MUM)
        SCRHET = 1.3_dp

     ELSE
        STOP 'VALUE IFRZ NOT IMPLEMENTED.'
     END IF

     RETURN

  END FUNCTION SCRHET

! ------------------------------------------------------------------------------

  REAL(dp) FUNCTION TAUG(B,Y,X0)

     ! ***** DIMENSIONLESS GROWTH TIME SCALE
     ! ***** B.KÄRCHER AND S.SOLOMON, JGR 104(D22), 27441-27459, 1999

     IMPLICIT NONE

     REAL(dp), INTENT(in) :: B, Y, X0

     REAL(dp), PARAMETER :: SQ31 = 0.577350269_dp
     REAL(dp), PARAMETER :: SIX1 = 0.166666666_dp

     REAL(dp) :: X, F1, F10, F2, F20

     TAUG = 0.0_dp

     IF (Y.LE.X0) RETURN

     X    = MIN( Y, 0.999_dp )
     F1   = SIX1 * LOG( (1.0_dp+X +X**2)  / (1.0_dp-X)**2 )
     F10  = SIX1 * LOG( (1.0_dp+X0+X0**2) / (1.0_dp-X0)**2 )
     F2   = SQ31 * ATAN( SQ31*(1.0_dp+2.0_dp*X) )
     F20  = SQ31 * ATAN( SQ31*(1.0_dp+2.0_dp*X0) )
     TAUG = (B+1.0_dp)*(F1-F10) + (B-1.0_dp)*(F2-F20)

     RETURN

  END FUNCTION TAUG

! ------------------------------------------------------------------------------

  SUBROUTINE XFRZMSTR(LHET,NOSIZE,ZTMST,                   &
                      KLEV,kbdim,kproma,JK, nfrzmod,       &
                      ZSUSATIX,ZVERVX,ZAPNX,               &
                      ZAPRX,ZAPSIGX,PTM1,TMELT,ZMIN,PAPM1, &
                      ZTHOMI,ZRI,ZNICEX)

    ! This is the calling loop for XFRZHOM/XFRZHET extracted from CLOUD

    IMPLICIT NONE

    ! SUBROUTINE parameters

    LOGICAL,  INTENT(in)  :: NOSIZE  ! .true. ---> aerosol size effects are
                                     !             ignored for homogeneous
                                     !             freezing
    LOGICAL,  INTENT(in)  :: LHET    ! .true. ---> heterogeneous freezing of
                                     !             aerosol particles below 235K
                                     !             is considered
                                     !             aerosol particles below 235K
                                     !             is considered

    INTEGER,  INTENT(in)  :: KLEV    ! number of levels
    INTEGER,  INTENT(in)  :: KBDIM   ! size of arrays
    INTEGER,  INTENT(in)  :: KPROMA  ! number of cells in arrays
    INTEGER,  INTENT(in)  :: JK      ! current level
    INTEGER,  INTENT(in)  :: NFRZMOD  ! number of aerosol modes
    REAL(dp), INTENT(in)  :: ZSUSATIX(kbdim)  ! ice supersaturation pressure
    REAL(dp), INTENT(in)  :: ZVERVX(kbdim)    ! updraft [cm/s]
    REAL(dp), INTENT(in)  :: ZAPNX(kbdim,NFRZMOD)  ! aerosol number conc. [1/cm3]
    REAL(dp), INTENT(in)  :: ZAPRX(kbdim,NFRZMOD)  ! aerosol radius [cm]
    REAL(dp), INTENT(in)  :: ZAPSIGX(kbdim,NFRZMOD)! aerosol standard deviation
                                                  ! (log-normal distribution)
    REAL(dp), INTENT(in)  :: PTM1(kbdim,KLEV) ! temperature [K]
    REAL(dp), INTENT(in)  :: PAPM1(kbdim,KLEV)! pressure [Pa]
    REAL(dp), INTENT(in)  :: ZTMST            ! time step [s]
    REAL(dp), INTENT(in)  :: TMELT            ! melting temperature
    REAL(dp), INTENT(in)  :: ZMIN             ! "epsilon"
    REAL(dp), INTENT(in)  :: ZTHOMI           ! temperature of hom. freezing

    REAL(dp) ZNICEX(kbdim,klev)    ! number of ice nuclei
    REAL(dp) ZRI(kbdim,KLEV)  ! 

    ! local variables

    REAL(dp) :: ZAPN(NFRZMOD)    ! aerosol number concentration [1/cm3]
    REAL(dp) :: ZAPR(NFRZMOD)    ! aerosol radius [cm]
    REAL(dp) :: ZAPSIG(NFRZMOD)  ! aerosol standard deviation (log-normal distr.)

    INTEGER  :: IXLIST(kbdim)
    INTEGER  :: IX2LIST(kbdim)
    INTEGER  :: KLIST(kbdim)
    REAL(dp) :: PWX(kbdim), COOLRX(kbdim)

    INTEGER  :: IXC,JL,IX,I,IX0,IXC2,K
    REAL(dp) :: ZSUSATI,ZVERV,ZSATI,COOLR,T,DT,SI,TMIN,SCR
    REAL(dp) :: PW,PWCR,SL,PISATTMIN,ZNICE,ZRI1,PICE,ZPHPA,TEMP,CTOT

    ! parameters

    INTEGER,  PARAMETER :: KMAX   = 120

    REAL(dp), PARAMETER :: RHOICE = 0.925_dp
    REAL(dp), PARAMETER :: PI     = 3.1415927_dp
    REAL(dp), PARAMETER :: ALPHA  = 0.5_dp
    REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
    REAL(dp), PARAMETER :: RGAS   = 8.3145E7_dp
    REAL(dp), PARAMETER :: BK     = 1.3807E-16_dp
    REAL(dp), PARAMETER :: CPAIR  = 1.00467E7_dp
    REAL(dp), PARAMETER :: HEAT   = 2830.3E7_dp
    REAL(dp), PARAMETER :: AVOG   = 6.02213E23_dp
    REAL(dp), PARAMETER :: GRAV   = 981.0_dp
    REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
    REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
    REAL(dp), PARAMETER :: XWW    = 18.016_dp
    REAL(dp), PARAMETER :: XWA    = 28.966_dp

!    ! functions
!
!    REAL(dp) :: PISAT  ! vapor pressure over ice [hPa]
!    REAL(dp) :: SCRHOM ! critical ice saturation ratio for homogeneous freezinh

    IXC=0

    DO 361 JL = 1, KPROMA
       ZSUSATI = ZSUSATIX(JL)
       ZVERV   = ZVERVX(JL)
       ZAPN(:) = ZAPNX(JL,:) !gf !SF related to #160: the previous MIN statement here was 1. buggy
                                 !                    and 2. unnecessary since it's taken care of right before
                                 ! the call to xfrzmstr in cloud_cdnc_icnc
       zri(jl,jk)=0._dp
       znicex(jl,jk)=0._dp

       IF (ZSUSATI.GT.0.0_dp.AND.PTM1(JL,JK).LT.ZTHOMI.AND.ZVERV.GT.ZMIN) THEN
          ZSATI=ZSUSATI+1.0_dp
          ZNICE=0.0_dp
          ZRI1=0.0_dp
          IF (PTM1(JL,JK).LT.ZTHOMI) THEN
             CTOT=0.0_dp
             DO I=1,NFRZMOD
                CTOT=CTOT+ZAPN(I)
             END DO
             IF ( CTOT > 0.0_dp .AND. ZVERV > 0.0_dp ) THEN

                !*****CHECK IF FREEZING CONDITIONS ARE MET WITHIN DT

                COOLR  = GRAV * ZVERV / CPAIR
                T = PTM1(JL,JK)
                DT = ZTMST
                SI = ZSATI
                TMIN   = T    - COOLR * DT
                TMIN   = MAX( TMIN, 170.0_dp )
                SCR    = SCRHOM(TMIN)
                PW     = SI * PISAT(T)
                PWCR   = PW * ( TMIN / T )**3.5_dp
                PISATTMIN=PISAT(TMIN)
                IF ((PWCR/PISATTMIN).GE. SCR) THEN
                   COOLRX(JL)=COOLR
                   PWX(JL)=PW
                   IXC=IXC+1
                   IXLIST(IXC)=JL
                END IF
             END IF
          END IF
       END IF
 361 END DO ! JL-loop

    ! *****FREEZING TEMPERATURE logic

    IX0=IXC
    IXC=0
    IXC2=0
    DO IX=1,IX0
       JL=IXLIST(IX)
       ZSUSATI=ZSUSATIX(JL)
       ZSATI=ZSUSATI+1.0_dp
       T = PTM1(JL,JK)
       SI=ZSATI
       SCR    = SCRHOM(T)
       IF (SI.GE.SCR) THEN
          ! Need not to figure out a K
          IXC=IXC+1
          IXLIST(IXC)=JL
       ELSE
          ! Need to figure out a K
          IXC2=IXC2+1
          IX2LIST(IXC2)=JL
       END IF
    END DO ! IX-loop

    DO K=1,KMAX+1
       IX0=IXC2
       IXC2=0
       IF ( IX0 > 0 ) THEN
          DO IX=1,IX0
             JL=IX2LIST(IX)
             T = PTM1(JL,JK)
             PW = PWX(JL)
             SL    = ( T - 170.0_dp ) / REAL(KMAX,dp)
             TEMP = T  - SL * REAL(K-1,dp)
             SCR  = SCRHOM(TEMP)
             PICE = PISAT(TEMP)
             PWCR = PW * (TEMP/T)**3.5_dp
             IF ((PWCR/PICE).GE.SCR) THEN
                ! The K condition has been reached.
                ! Keep the K in mind and add this JL point to the list
                ! of points to be passed to XFRZHOM
                KLIST(JL)=K
                IXC=IXC+1
                IXLIST(IXC)=JL
             ELSE
                ! Keep for further K iteration
                IXC2=IXC2+1
                IX2LIST(IXC2)=JL
             END IF
          END DO ! IX-loop
       END IF
    END DO ! K-loop

    DO IX=1,IXC
       JL=IXLIST(IX)
       ZSUSATI=ZSUSATIX(JL)
       ZVERV=ZVERVX(JL)
       ZAPN(:)=ZAPNX(JL,:)
       ZAPR(:)=ZAPRX(JL,:)
       ZAPSIG(:)=ZAPSIGX(JL,:)
       COOLR=COOLRX(JL)
       PW=PWX(JL)

       ZPHPA = PAPM1(JL,JK)*0.01_dp
       ZSATI = ZSUSATI+1.0_dp
       ZNICE = 0.0_dp
       ZRI1  = 0.0_dp

       T = PTM1(JL,JK)
       SI = ZSATI

       ! ***** FREEZING TEMPERATURE

       SCR    = SCRHOM(T)
       IF (SI.GE.SCR) THEN
          TEMP  = T
          PICE  = PISAT(TEMP)
          PW    = SCR * PICE
          PWCR  = PW
       ELSE
          SL   = ( T - 170.0_dp ) / REAL(KMAX,dp)
          K    = KLIST(JL)
          TEMP = T  - SL * REAL(K-1,dp)
          SCR  = SCRHOM(TEMP)
          PICE = PISAT(TEMP)
          PWCR = PW * (TEMP/T)**3.5_dp
       END IF

       IF (LHET) THEN
          CALL XFRZHET (NOSIZE,nfrzmod,ZTMST,ZAPN,ZAPR,ZAPSIG,ZPHPA,    &
                        PTM1(JL,JK),ZVERV,ZSATI,ZNICE,ZRI1, COOLR, PW, &
                        SCR, TEMP, PICE, PWCR)
       ELSE
          CALL XFRZHOM (NOSIZE,nfrzmod,ZTMST,ZAPN,ZAPR,ZAPSIG,ZPHPA,    &
                        PTM1(JL,JK),ZVERV,ZSATI,ZNICE,ZRI1, COOLR, PW, &
                        SCR, TEMP, PICE, PWCR)
       END IF

       ZRI(JL,JK) = ZRI1 *1.0E-2_dp
       ZNICEX(JL,JK) = ZNICE*1.0E6_dp

    END DO ! IX-loop

    RETURN

  END SUBROUTINE XFRZMSTR

! ------------------------------------------------------------------------------

  SUBROUTINE XFRZHET (NOSIZE,nfrzmod,DT,C,R,SIG,P,T,V,SI,CI,RI,COOLR, PW, &
                      SCR, TEMP, PICE, PWCR)

     ! ***** HETEROGENEOUS FREEZING (ADIABATIC ASCENT)
     !
     ! ***** BERND KÄRCHER  APR 14  2003
     ! ***** bernd.kaercher@dlr.de  http://www.op.dlr.de/~pa1c/
     !
     ! ***** Ulrike Lohmann      Dalhousie Univ   Apr 14 03
     !       Ulrike.Lohmann@Dal.Ca
     ! ***** Johannes Hendricks  DLR-IPA          Apr 14 03
     !       Johannes.Hendricks@dlr.de
     !
     ! ***** References
     !       Kärcher, B. and U. Lohmann
     !       A parameterization of cirrus cloud formation: Heterogeneous
     !       freezing.
     !       J. Geophys. Res. 108, doi:10.1029/2002JD003220, in press, 2003.

     IMPLICIT NONE

     ! SUBROUTINE parameters

     LOGICAL,  INTENT(in)  :: NOSIZE     ! .true. ---> aerosol size effects are
                                         !             ignored for homogeneous
                                         !             freezing
     INTEGER,  INTENT(in)  :: NFRZMOD     ! number of aerosol modes
     REAL(dp), INTENT(in)  :: DT         ! time step [s]
     REAL(dp), INTENT(in)  :: C(NFRZMOD)  ! aerosol number concentration [1/cm3]
     REAL(dp), INTENT(in)  :: R(NFRZMOD)  ! aerosol radius [cm]
     REAL(dp), INTENT(in)  :: SIG(NFRZMOD)! aerosol standard deviation
                                         ! (log-normal distribution)
     REAL(dp), INTENT(in)  :: P          ! pressure [hPa]
     REAL(dp), INTENT(in)  :: T          ! temperature [K]
     REAL(dp), INTENT(in)  :: V          ! updraft [cm/s]

     REAL(dp), INTENT(inout) :: SI       ! saturation pressure (ice)

     REAL(dp), INTENT(out) :: COOLR      ! cooling rate
     REAL(dp), INTENT(out) :: RI         ! 
     REAL(dp), INTENT(out) :: PICE       ! 
     REAL(dp), INTENT(out) :: TEMP       ! 
     REAL(dp), INTENT(out) :: CI         ! ice nuclei
     REAL(dp), INTENT(out) :: SCR        ! 
     REAL(dp), INTENT(out) :: PW         ! 
     REAL(dp), INTENT(out) :: PWCR       ! 

     ! local variables

     INTEGER, SAVE :: NFRZ = 0

     INTEGER :: IFRZ,N,K

     REAL(dp) :: A(3)
     REAL(dp) :: B(2)
     REAL(dp) :: CCR(nfrzmod)

     REAL(dp) :: TMIN,SISVE,SL,VOLF,PCR,SUPSCR,DLNJDT,TAU,DIFFC
     REAL(dp) :: BKT,VTH,CISAT,THETA,XMI0,XMIMAX,RIMAX,XMISAT,CTAU
     REAL(dp) :: TGROW,ZF,XMFP,BETA,RIHAT,X0,Z,XMI,YK,CVF,CTOT,DSCRMIN
     REAL(dp) :: RS,X

     ! parameters

     INTEGER,  PARAMETER :: KMAX   = 120

     REAL(dp), PARAMETER :: RHOICE = 0.925_dp
     REAL(dp), PARAMETER :: PI     = 3.1415927_dp
     REAL(dp), PARAMETER :: SQPI   = 1.7724539_dp
     REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
     REAL(dp), PARAMETER :: RGAS   = 8.3145E7_dp
     REAL(dp), PARAMETER :: ALPHA  = 0.5_dp
     REAL(dp), PARAMETER :: BK     = 1.3807E-16_dp
     REAL(dp), PARAMETER :: CPAIR  = 1.00467E7_dp
     REAL(dp), PARAMETER :: HEAT   = 2830.3E7_dp
     REAL(dp), PARAMETER :: AVOG   = 6.02213E23_dp
     REAL(dp), PARAMETER :: GRAV   = 981.0_dp
     REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
     REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
     REAL(dp), PARAMETER :: XWW    = 18.016_dp
     REAL(dp), PARAMETER :: XWA    = 28.966_dp

!     ! functions
!
!     REAL(dp) :: SCRHET ! critical ice saturation ratio for hom. freezing
!     REAL(dp) :: TAUG   ! dimensionless growth time scale
!     REAL(dp) :: PISAT  ! vapor pressure over ice [hPa]

     IFRZ    = 0
     CVF     = 0.99_dp
     CTAU    = 50.0_dp
     DSCRMIN = 0.05_dp

     ! ***** DO NOTHING CRITERIA

     CI    = 0.0_dp
     RI    = 0.0_dp
     XMI   = 0.0_dp
     CTOT  = 0.0_dp
     DO N  = 1, NFRZMOD
        CTOT = CTOT + C(N)
     END DO
     IF ((CTOT.LE.0.0_dp).OR.(V.LE.0.0_dp)) RETURN
     NFRZ  = IFRZ

     ! ***** CHECK IF FREEZING CONDITIONS ARE MET WITHIN DT

     COOLR  = GRAV * V     / CPAIR
     TMIN   = T    - COOLR * DT
     TMIN   = MAX( TMIN, 170.0_dp )
     SCR    = SCRHET(TMIN,NFRZ)
     PW     = SI * PISAT(T)
     PWCR   = PW * ( TMIN / T )**3.5_dp
     IF ((PWCR/PISAT(TMIN)).LT.SCR) RETURN

     ! ***** FREEZING TEMPERATURE

     SCR    = SCRHET(T,NFRZ)
     IF (SI.GE.SCR) THEN
        SISVE = SCR
        TEMP  = T
        PICE  = PISAT(TEMP)
        PW    = SCR * PICE
        PWCR  = PW
        GOTO 10
     ELSE
        SISVE = SI
        SL    = ( T - 170.0_dp ) / REAL(KMAX,dp)
        DO K = 1, KMAX+1
           TEMP = T  - SL * REAL(K-1,dp)
           SCR  = SCRHET(TEMP,NFRZ)
           PICE = PISAT(TEMP)
           PWCR = PW * (TEMP/T)**3.5_dp
           IF ((PWCR/PICE).GE.SCR) GOTO 10
        END DO ! K-loop
     END IF

     RETURN

 10  CONTINUE

     VOLF    = PI * RHOICE / 0.75_dp
     PCR     = P  * PWCR   / PW
     SCR     = MAX( SCR, 1.001_dp )
     SUPSCR  = SCR - 1.0_dp + DSCRMIN
     DO N    = 1, NFRZMOD
        CCR(N) = C(N) * PWCR / PW
     END DO ! N-loop

     ! ***** TIMESCALE OF THE FREEZING EVENT (FROM HOMOGENEOUS RATE)

     DLNJDT = ABS( 4.37_dp - 0.03_dp*TEMP )
     TAU    = 1.0_dp / ( CTAU*DLNJDT*COOLR )

     ! ***** ICE CRYSTAL PROPERTIES AFTER THE FREEZING EVENT

     DIFFC  = 4.0122E-3_dp * TEMP**1.94_dp / PCR
     BKT    = BK   * TEMP
     VTH    = SQRT( 8.0_dp*BKT / (PI*XMW) )
     CISAT  = 1.0E3_dp * PICE / BKT
     THETA  = HEAT * XWW / (RGAS*TEMP)
     A(1)   = ( THETA/CPAIR - XWA/RGAS ) * ( GRAV/TEMP )
     A(2)   = 1.0_dp   / CISAT
     A(3)   = 0.001_dp*XWW**2*HEAT**2 / ( AVOG*XWA*PCR*TEMP*CPAIR )
     B(1)   = SVOL * 0.25_dp  * ALPHA * VTH * CISAT * SUPSCR
     B(2)   = 0.25_dp * ALPHA * VTH   / DIFFC

     CALL XICEHET (NFRZMOD,V,TEMP,TAU,SCR,A,B,CCR,R,SIG,CI,RIHAT,RS,YK)

     XMI0   = VOLF * CI  * RIHAT**3  * CVF
     XMIMAX = XMI0 + XMW * CISAT * ( SCR - 1.0_dp )
     RIMAX  = ( XMIMAX / (VOLF*CI) )**THIRD

     ! ***** VAPOR RELAXATION: ICE CRYSTAL PROPERTIES AFTER DT

     XMISAT = XMW   * CISAT
     TGROW  = 0.75_dp  / ( PI*DIFFC*CI*RIMAX )
     ZF     = TGROW / DT
     XMFP   = 3.0_dp    * DIFFC / VTH
     BETA   = XMFP  / ( 0.75_dp*ALPHA*RIMAX )
     X0     = RIHAT / RIMAX
!     DO X   = 1.0_dp, X0, -0.01_dp
     X = 1.0_dp
     DO WHILE (X >= X0)
        Z     = ZF * TAUG(BETA,X,X0)
        IF (Z.LE.1.0_dp) EXIT
        X = X - 0.01_dp
     END DO

     RI     = X    * RIMAX
     XMI    = VOLF * CI * RI**3
     SI     = SCR  - ( XMI - XMI0 )/XMISAT

     RETURN

  END SUBROUTINE XFRZHET

! ------------------------------------------------------------------------------

  SUBROUTINE XICEHET (NFRZMOD,V,T,TAU,SCR,A,B,C,R,SIG,CI,RIHAT,RS,YK)

    ! ***** ICE CRYSTAL CONCENTRATION AND SIZE AFTER FREEZING EVENT

    IMPLICIT NONE

    ! subroutine parameters

    INTEGER,  INTENT(in)  :: NFRZMOD      ! number of aerosol modes
    REAL(dp), INTENT(in)  :: V           ! updraft [cm/s]
    REAL(dp), INTENT(in)  :: T           ! temperature [K]
    REAL(dp), INTENT(in)  :: TAU         ! time scale of the freezing event
    REAL(dp), INTENT(in)  :: SCR         ! 
    REAL(dp), INTENT(in)  :: A(3)        ! 
    REAL(dp), INTENT(in)  :: B(2)        ! 
    REAL(dp), INTENT(in)  :: C(nfrzmod)   ! 
    REAL(dp), INTENT(in)  :: R(nfrzmod)   ! aerosol radius [cm]
    REAL(dp), INTENT(in)  :: SIG(nfrzmod) ! aerosol standard deviation
                                         ! (log-normal distribution)
    REAL(dp), INTENT(out) :: CI          ! ice nuclei
    REAL(dp), INTENT(out) :: RIHAT       ! 
    REAL(dp), INTENT(out) :: RS          ! 
    REAL(dp), INTENT(out) :: YK          ! 

    ! parameters

    INTEGER,  PARAMETER :: IBIN   = 80

    REAL(dp), PARAMETER :: RMIN   = 1.0E-7_dp
    REAL(dp), PARAMETER :: RMAX   = 1.0E-3_dp
    REAL(dp), PARAMETER :: VRAT   = 1.5_dp
    REAL(dp), PARAMETER :: PI     = 3.1415927_dp
    REAL(dp), PARAMETER :: SQPI   = 1.7724539_dp
    REAL(dp), PARAMETER :: TWOPI  = 6.2831853_dp
    REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
    REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
    REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
    REAL(dp), PARAMETER :: RHOICE = 0.925_dp

    ! local variables

    REAL(dp) :: SHAPE(IBIN),   R0(IBIN+1)
    REAL(dp) :: CRINT(IBIN+1), CIINT(IBIN+1)
    REAL(dp) :: VOLF,PHI,RIMFC,TBBT,CTOT,DELTA,DELP1,SYK,EERFC,RIM,ARG,RMEAN
    REAL(dp) :: SLOPE,SLOPER,RLOGRAT,CSH,SHAPFC1,SHAPFC2,SIGL,RIRAT,SUMSH
    REAL(dp) :: XMIHAT,DLOGR0

    INTEGER :: NBIN,N,I,II

!    ! functions
!
!    REAL(dp) :: XEERFC

    ! ***** CONSTANTS

    NBIN  = 1 + INT( LOG( (RMAX/RMIN)**3 ) / LOG(VRAT) )
    VOLF  = PI * RHOICE    / 0.75_dp
    PHI   = V  * A(1)*SCR / ( A(2) + A(3)*SCR )
    RIMFC = 4.0_dp * PI   * B(1)/B(2)**2 / SVOL
    TBBT  = 2.0_dp * B(1) * B(2) * TAU
    CTOT  = 0.0_dp
    DO N  = 1, NFRZMOD
       CTOT = CTOT + C(N)
    END DO

    ! ***** MONODISPERSE APPROXIMATION (SINGLE MODE ONLY)

    IF ((NFRZMOD.EQ.1).AND.(SIG(NFRZMOD).LT.1.1)) THEN
       RS     = R(NFRZMOD)
       DELTA  = B(2) * RS
       DELP1  = 1.0_dp   + DELTA
       YK     = TBBT / DELP1**2
       SYK    = SQRT(YK)
       EERFC  = XEERFC(YK)
       RIM    = RIMFC / DELP1 * ( DELTA**2 - 1.0_dp &
                                  + (1.0_dp+0.5_dp*YK*DELP1**2)*SQPI*EERFC/SYK )
       CI     = PHI / RIM
       RIRAT  = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
       RIHAT  = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
       XMIHAT = VOLF * CI * RIHAT**3
       CI     = MIN( CI, CTOT )
       RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
       RETURN
    END IF

    ! ***** SIZE DISTRIBUTION PROPERTIES

    R0(NBIN+1)    = RMAX * VRAT**THIRD
    CIINT(NBIN+1) = 1.0E-35_dp
    CRINT(NBIN+1) = 1.0E-25_dp
    DLOGR0     = 2.0_dp**THIRD * (VRAT**THIRD-1.0_dp) / (VRAT+1.0_dp)**THIRD
    SUMSH      = 0.0_dp
    DO I       = 1, NBIN
       CIINT(I)  = 0.0_dp
       CRINT(I)  = 0.0_dp
       SHAPE(I)  = 0.0_dp
       R0(I)     = RMIN * VRAT**( THIRD*REAL(I-1,dp) )
       DO N      = 1, NFRZMOD
          SIGL     = LOG( MAX( SIG(N), 1.1_dp ) )
          SHAPFC1  = 1.0_dp   / ( SQRT(TWOPI) * SIGL )
          SHAPFC2  = 0.5_dp / SIGL**2
          ARG      = SHAPFC2  * ( LOG(R0(I)/R(N)) )**2
          ARG      = MIN( ARG, 75.0_dp )
          SHAPE(I) = SHAPE(I) + DLOGR0 * C(N) * SHAPFC1 * EXP( -ARG )
       END DO ! N-loop
       SUMSH     = SUMSH    + SHAPE(I)
    END DO ! I-loop

    ! ***** ICE CRYSTAL PROPERTIES

    RMEAN     = 0.0_dp
    DO I  = NBIN, 1, -1
       DELTA    = B(2) * R0(I)
       DELP1    = 1.0_dp   + DELTA
       YK       = TBBT / DELP1**2
       SYK      = SQRT(YK)
       EERFC    = XEERFC(YK)
       RIM      = RIMFC / DELP1 * ( DELTA**2 - 1.0_dp &
                                  + (1.0_dp+0.5_dp*YK*DELP1**2)*SQPI*EERFC/SYK )
       CSH      = SHAPE(I) / SUMSH   * CTOT
       CRINT(I) = CRINT(I+1) + RIM   * CSH
       CIINT(I) = CIINT(I+1) +         CSH
       RMEAN    = RMEAN      + R0(I) * CSH
       IF (CRINT(I).GE.PHI) GOTO 10
    END DO ! I-loop
    RS        = R0(1)
    RMEAN     = RMEAN / CTOT
    CI        = CTOT  * PHI   / CRINT(1)
    DELP1     = 1.0_dp    + B(2)  * RMEAN
    YK        = TBBT  / DELP1**2
    SYK       = SQRT(YK)
    EERFC     = XEERFC(YK)
    RIRAT     = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
    RIHAT     = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
    XMIHAT    = VOLF    * CI * RIHAT**3
    CI        = CTOT
    RIHAT     = ( XMIHAT / (VOLF*CI) )**THIRD
    RETURN
 10 CONTINUE
    II        = MAX( I, 1 )
    RLOGRAT   = LOG( R0(II) / R0(II+1) )
    SLOPER    = LOG( CRINT(II) / CRINT(II+1) ) / RLOGRAT
    RS        = R0(II+1)    * ( PHI / CRINT(II+1) )**(1.0_dp/SLOPER)
    SLOPE     = LOG( CIINT(II) / CIINT(II+1) ) / RLOGRAT
    CI        = CIINT(II+1) * ( RS / R0(II+1) )**SLOPE
    RMEAN     = RMEAN / CTOT
    DELP1     = 1.0_dp   + B(2) * MAX( RS, RMEAN )
    YK        = TBBT / DELP1**2
    SYK       = SQRT(YK)
    EERFC     = XEERFC(YK)
    RIRAT     = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
    RIHAT     = ( RIRAT * DELP1 - 1.0_dp ) / B(2)

    ! ***** EXIT
    RETURN

  END SUBROUTINE XICEHET

! ------------------------------------------------------------------------------

  SUBROUTINE XFRZHOM (NOSIZE,nfrzmod,DT,C,R,SIG,P,T,V,SI,CI,RI,  &
                      COOLR, PW, SCR, TEMP, PICE, PWCR)

     ! ***** HOMOGENEOUS FREEZING OF SUPERCOOLED AEROSOL (ADIABATIC ASCENT)
     !
     ! ***** BERND KÄRCHER  APR 30  2002
     ! ***** bernd.kaercher@dlr.de  http://www.op.dlr.de/~pa1c/
     ! ***** Ulrike Lohmann  Dalhousie Univ  Apr 02  Ulrike.Lohmann@Dal.Ca
     !
     ! ***** References
     ! Kärcher, B. and U. Lohmann
     ! A parameterization of cirrus cloud formation:
     ! Homogeneous freezing including effects of aerosol size.
     ! J. Geophys. Res. 107 (D ), in press, 2002.

     IMPLICIT NONE

     ! subroutine parameters
 
     LOGICAL,  INTENT(in)  :: NOSIZE      ! take into account aerosol size?
     INTEGER,  INTENT(in)  :: NFRZMOD      ! number of aerosol modes
     REAL(dp), INTENT(in)  :: DT          ! time step [s]
     REAL(dp), INTENT(in)  :: C(NFRZMOD)   !
     REAL(dp), INTENT(in)  :: R(NFRZMOD)   ! aerosol radius [cm]
     REAL(dp), INTENT(in)  :: SIG(NFRZMOD) ! aerosol standard deviation
                                          ! (log-normal size distribution)
     REAL(dp), INTENT(in)  :: P           ! pressure [hPa]
     REAL(dp), INTENT(in)  :: T           ! temperature [K]
     REAL(dp), INTENT(in)  :: V           ! updraft [cm/s]
     REAL(dp), INTENT(in)  :: COOLR       ! 
     REAL(dp), INTENT(in)  :: PW          ! 
     REAL(dp), INTENT(in)  :: SCR         ! 
     REAL(dp), INTENT(in)  :: TEMP        ! 
     REAL(dp), INTENT(in)  :: PICE        ! 
     REAL(dp), INTENT(in)  :: PWCR        ! 

     REAL(dp), INTENT(out) :: SI          ! saturation pressure (ice)
     REAL(dp), INTENT(out) :: CI          ! ice nuclei
     REAL(dp), INTENT(out) :: RI          ! 

     ! parameters

     INTEGER,  PARAMETER :: KMAX   = 120

     REAL(dp), PARAMETER :: RHOICE = 0.925_dp
     REAL(dp), PARAMETER :: PI     = 3.1415927_dp
     REAL(dp), PARAMETER :: ALPHA  = 0.5_dp
     REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
     REAL(dp), PARAMETER :: RGAS   = 8.3145E7_dp
     REAL(dp), PARAMETER :: BK     = 1.3807E-16_dp
     REAL(dp), PARAMETER :: CPAIR  = 1.00467E7_dp
     REAL(dp), PARAMETER :: HEAT   = 2830.3E7_dp
     REAL(dp), PARAMETER :: AVOG   = 6.02213E23_dp
     REAL(dp), PARAMETER :: GRAV   = 981.0_dp
     REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
     REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
     REAL(dp), PARAMETER :: XWW    = 18.016_dp
     REAL(dp), PARAMETER :: XWA    = 28.966_dp

     ! local variables

     REAL(dp) :: A(3)        ! 
     REAL(dp) :: B(2)        ! 
     REAL(dp) :: CCR(NFRZMOD) ! 

     REAL(dp) :: XMI,VOLF,PCR,CTAU,DLNJDT,TAU,DIFFC,BKT,VTH,CISAT,THETA
     REAL(dp) :: YK,RS,RIHAT,XMI0,XMIMAX,RIMAX,XMISAT,TGROW,ZF,XMFP
     REAL(dp) :: BETA,X0,X,Z

     INTEGER  :: N

!     ! functions
!
!     REAL(dp) :: PISAT  ! vapor pressure over ice [hPa]
!     REAL(dp) :: TAUG   ! dimensionless growth time scale

     ! ***** DO NOTHING CRITERIA

     CI    = 0.0_dp
     RI    = 0.0_dp
     XMI   = 0.0_dp

     ! ***** FREEZING TEMPERATURE

!!$CC      SCR    = SCRHOM(T)
!!$CC      IF (SI.GE.SCR) THEN
!!$CC       SISVE = SCR
!!$CC       TEMP  = T
!!$CC       PICE  = PISAT(TEMP)
!!$CC       PW    = SCR * PICE
!!$CC       PWCR  = PW
!!$CC       GOTO 10
!!$CC      ELSE
!!$CC       SISVE = SI
!!$CC       SL    = ( T - 170. ) / REAL(KMAX,dp)
!!$CC       DO  K = 1, KMAX+1
!!$CC        TEMP = T  - SL * REAL(K-1,dp)
!!$CC        SCR  = SCRHOM(TEMP)
!!$CC        PICE = PISAT(TEMP)
!!$CC        PWCR = PW * (TEMP/T)**3.5
!!$CC        IF ((PWCR/PICE).GE.SCR) GOTO 10
!!$CC       ENDDO
!!$CC      ENDIF
!!$CC      RETURN
!!$CC 10   CONTINUE

     VOLF    = PI * RHOICE / 0.75_dp
     PCR     = P  * PWCR   / PW
     DO N    = 1, NFRZMOD
        CCR(N) = C(N) * PWCR / PW
     END DO

     ! ***** TIMESCALE OF THE FREEZING EVENT

     IF (NOSIZE) THEN
        CTAU  = MAX( (2260.0_dp-10.0_dp*TEMP), 100.0_dp )
     ELSE
        CTAU  = 50.0_dp
     END IF
     DLNJDT = ABS( 4.37_dp - 0.03_dp*TEMP )
     TAU    = 1.0_dp / ( CTAU*DLNJDT*COOLR )

     ! ***** ICE CRYSTAL PROPERTIES AFTER THE FREEZING EVENT

     DIFFC  = 4.0122E-3_dp * TEMP**1.94_dp / PCR
     BKT    = BK   * TEMP
     VTH    = SQRT( 8.0_dp*BKT / (PI*XMW) )
     CISAT  = 1.E3_dp * PICE / BKT
     THETA  = HEAT * XWW / (RGAS*TEMP)
     A(1)   = ( THETA/CPAIR - XWA/RGAS ) * ( GRAV/TEMP )
     A(2)   = 1.0_dp   / CISAT
     A(3)   = 0.001_dp*XWW**2*HEAT**2 / ( AVOG*XWA*PCR*TEMP*CPAIR )
     B(1)   = SVOL * 0.25_dp  * ALPHA * VTH * CISAT * ( SCR - 1.0_dp )
     B(2)   = 0.25_dp * ALPHA * VTH   / DIFFC

     CALL XICEHOM (NOSIZE,NFRZMOD,V,TEMP,TAU,SCR,A,B,CCR,R,SIG,CI,RIHAT,RS,YK)

     XMI0   = VOLF * CI  * RIHAT**3
     XMIMAX = XMI0 + XMW * CISAT * ( SCR - 1.0_dp )
     RIMAX  = ( XMIMAX / (VOLF*CI) )**THIRD

     ! ***** VAPOR RELAXATION: ICE CRYSTAL PROPERTIES AFTER DT

     XMISAT = XMW   * CISAT
     TGROW  = 0.75_dp  / ( PI*DIFFC*CI*RIMAX )
     ZF     = TGROW / DT
     XMFP   = 3.0_dp    * DIFFC / VTH
     BETA   = XMFP  / ( 0.75_dp*ALPHA*RIMAX )
     X0     = RIHAT / RIMAX
!     DO X   = 1.0_dp, X0, -0.01_dp
     X = 1.0_dp
     DO WHILE (X >= X0)
        Z = ZF * TAUG(BETA,X,X0)
        IF (Z.LE.1.0_dp) EXIT
        X = X - 0.01_dp
     END DO

     RI     = X    * RIMAX
     XMI    = VOLF * CI * RI**3
     SI     = SCR  - ( XMI - XMI0 )/XMISAT

     RETURN

  END SUBROUTINE XFRZHOM

! ------------------------------------------------------------------------------

  SUBROUTINE XICEHOM (NOSIZE,NFRZMOD,V,T,TAU,SCR,A,B,C,R,SIG,CI,RIHAT,RS,YK)

     ! ***** ICE CRYSTAL CONCENTRATION AND SIZE AFTER FREEZING EVENT

     IMPLICIT NONE

     ! subroutine parameters

     LOGICAL,  INTENT(in)  :: NOSIZE      ! take into account aerosol size?
     INTEGER,  INTENT(in)  :: NFRZMOD      ! number of aerosol modes
     REAL(dp), INTENT(in)  :: V           ! updraft [cm/s]
     REAL(dp), INTENT(in)  :: T           ! temperature [K]
     REAL(dp), INTENT(in)  :: TAU         ! 
     REAL(dp), INTENT(in)  :: SCR         ! 
     REAL(dp), INTENT(in)  :: A(3)        ! 
     REAL(dp), INTENT(in)  :: B(2)        ! 
     REAL(dp), INTENT(in)  :: C(NFRZMOD)
     REAL(dp), INTENT(in)  :: R(NFRZMOD)   ! aerosol radius [cm]
     REAL(dp), INTENT(in)  :: SIG(NFRZMOD) ! aerosol standard deviation
                                          ! (log-normal size distribution)
     REAL(dp), INTENT(out) :: CI          ! ice nuclei
     REAL(dp), INTENT(out) :: RIHAT       ! 
     REAL(dp), INTENT(out) :: RS          ! 
     REAL(dp), INTENT(out) :: YK          ! 

     ! parameters

     INTEGER,  PARAMETER :: IBIN   = 80
     REAL(dp), PARAMETER :: RMIN   = 1.0E-7_dp
     REAL(dp), PARAMETER :: RMAX   = 1.0E-3_dp
     REAL(dp), PARAMETER :: VRAT   = 1.5_dp
     REAL(dp), PARAMETER :: PI     = 3.1415927_dp
     REAL(dp), PARAMETER :: SQPI   = 1.7724539_dp
     REAL(dp), PARAMETER :: TWOPI  = 6.2831853_dp
     REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
     REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
     REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
     REAL(dp), PARAMETER :: RHOICE = 0.925_dp

     ! local variables

     REAL(dp) :: SHAPE(IBIN), R0(IBIN+1)
     REAL(dp) :: CRINT(IBIN+1), CIINT(IBIN+1)
     REAL(dp) :: RIMX(IBIN+1), CSHX(IBIN+1)

     REAL(dp) :: VOLF,PHI,RIMFC,TBBT,CTOT,XMIHAT,DELTA,DELP1,SYK,RIM,RIRAT
     REAL(dp) :: DLOGR0,SUMSH,SIGL,SHAPFC1,SHAPFC2,ARG,RMEAN,CSH,RLOGRAT
     REAL(dp) :: SLOPER,SLOPE,EERFC

     INTEGER  :: NBIN,N,I,II

!     ! functions
!
!     REAL(dp) :: XEERFC

     ! ***** CONSTANTS

     NBIN  = 1 + INT( LOG( (RMAX/RMIN)**3 ) / LOG(VRAT) )
     VOLF  = PI * RHOICE    / 0.75_dp
     PHI   = V  * A(1)*SCR / ( A(2) + A(3)*SCR )
     RIMFC = 4.0_dp * PI   * B(1)/B(2)**2 / SVOL
     TBBT  = 2.0_dp * B(1) * B(2) * TAU
     CTOT  = 0.0_dp
     DO N  = 1, NFRZMOD
        CTOT = CTOT + C(N)
     END DO

     ! ***** NO SIZE EFFECTS

     IF (NOSIZE) THEN
        RS     = 0.25E-4_dp
        YK     = TAU  * ( B(1)/RS ) / ( 1.0_dp + B(2)*RS )
        CI     = SVOL * ( B(2) / (TWOPI*B(1)) )**1.5_dp * PHI / SQRT(TAU)
        CI     = MIN( CI, CTOT )
        XMIHAT = XMW * PI * PHI * TAU / 6.0_dp
        RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
        RETURN
     END IF

     ! ***** MONODISPERSE AEROSOL (SINGLE MODE ONLY)

     IF ((NFRZMOD.EQ.1).AND.(SIG(NFRZMOD).LT.1.1)) THEN
        RS     = R(NFRZMOD)
        DELTA  = B(2) * RS
        DELP1  = 1.0_dp   + DELTA
        YK     = TBBT / DELP1**2
        SYK    = SQRT(YK)
        EERFC  = XEERFC(YK)
        RIM    = RIMFC / DELP1 * ( DELTA**2 - 1.0_dp   &
                                  + (1.0_dp+0.5_dp*YK*DELP1**2)*SQPI*EERFC/SYK )
        CI     = PHI / RIM
        RIRAT  = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
        RIHAT  = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
        XMIHAT = VOLF * CI * RIHAT**3
        CI     = MIN( CI, CTOT )
        RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
        RETURN
     END IF

     ! ***** SIZE DISTRIBUTION PROPERTIES

     R0(NBIN+1)    = RMAX * VRAT**THIRD
     CIINT(NBIN+1) = 1.E-35_dp
     CRINT(NBIN+1) = 1.E-25_dp
     DLOGR0     = 2.0_dp**THIRD * (VRAT**THIRD-1.0_dp) / (VRAT+1.0_dp)**THIRD
     SUMSH      = 0.0_dp
     DO I       = 1, NBIN
        CIINT(I)  = 0.0_dp
        CRINT(I)  = 0.0_dp
        SHAPE(I)  = 0.0_dp
        R0(I)     = RMIN * VRAT**( THIRD*REAL(I-1,dp) )
        DO N      = 1, NFRZMOD
           SIGL     = LOG( MAX( SIG(N), 1.1_dp ) )
           SHAPFC1  = 1.0_dp   / ( SQRT(TWOPI) * SIGL )
           SHAPFC2  = 0.5_dp / SIGL**2
           ARG      = SHAPFC2  * ( LOG(R0(I)/R(N)) )**2
           ARG      = MIN( ARG, 75.0_dp )
           SHAPE(I) = SHAPE(I) + DLOGR0 * C(N) * SHAPFC1 * EXP( -ARG )
        END DO
        SUMSH     = SUMSH    + SHAPE(I)
     END DO

     ! ***** ICE CRYSTAL PROPERTIES

     RMEAN     = 0.0_dp
     DO     I  = NBIN, 1, -1
        DELTA    = B(2) * R0(I)
        DELP1    = 1.0_dp   + DELTA
        YK       = TBBT / DELP1**2
        SYK      = SQRT(YK)
        EERFC    = XEERFC(YK)
        RIMX(I)  = RIMFC / DELP1 * ( DELTA**2 - 1.0_dp  &
                                  + (1.0_dp+0.5_dp*YK*DELP1**2)*SQPI*EERFC/SYK )
        CSHX(I)  = SHAPE(I) / SUMSH   * CTOT
     END DO
     DO     I  = NBIN, 1, -1
        RIM = RIMX(I)
        CSH = CSHX(I)
        CRINT(I) = CRINT(I+1) + RIM   * CSH
        CIINT(I) = CIINT(I+1) +         CSH
        RMEAN    = RMEAN      + R0(I) * CSH
        IF (CRINT(I).GE.PHI) GOTO 10
     END DO
     RS        = R0(1)
     RMEAN     = RMEAN / CTOT
     CI        = CTOT  * PHI   / CRINT(1)
     DELP1     = 1.0_dp    + B(2)  * RMEAN
     YK        = TBBT  / DELP1**2
     SYK       = SQRT(YK)
     EERFC     = XEERFC(YK)
     RIRAT     = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
     RIHAT     = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
     XMIHAT    = VOLF    * CI * RIHAT**3
     CI        = CTOT
     RIHAT     = ( XMIHAT / (VOLF*CI) )**THIRD

     RETURN

 10  CONTINUE

     II        = MAX( I, 1 )
     RLOGRAT   = LOG( R0(II) / R0(II+1) )
     SLOPER    = LOG( CRINT(II) / CRINT(II+1) ) / RLOGRAT
     RS        = R0(II+1)    * ( PHI / CRINT(II+1) )**(1.0_dp/SLOPER)
     SLOPE     = LOG( CIINT(II) / CIINT(II+1) ) / RLOGRAT
     CI        = CIINT(II+1) * ( RS / R0(II+1) )**SLOPE
     RMEAN     = RMEAN / CTOT
     DELP1     = 1.0_dp   + B(2) * MAX( RS, RMEAN )
     YK        = TBBT / DELP1**2
     SYK       = SQRT(YK)
     EERFC     = XEERFC(YK)
     RIRAT     = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
     RIHAT     = ( RIRAT * DELP1 - 1.0_dp ) / B(2)

     ! **** EXIT
     RETURN

  END SUBROUTINE XICEHOM

! ------------------------------------------------------------------------------

END MODULE mo_cirrus
