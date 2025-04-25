!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE soilhyd(kidia, kfdia, klon, kdeep, lperma, zdiagt, ilog, jllog,  &
    wsi, dzsi, zwsat, zwsfc, fksat, vpor, bclapp, fmpot,  &
    zdrain, spsi, zwpwp, zdel, ice)
!     **************************************************************************

!     ******* SOIL HYDROLOGY - Routine to be build into REMO
!     *******   to calculate soil hydrology for 5 soil layers
!     ******* SOIL HYDROLOGY = percolation and vertical diffusion of water
!     ******* Type of diffusiont:  Richtmeyr/Morton

!     *** PROGRAMMING AND DEVELOPMENT: STEFAN HAGEMANN
!     ****
!     *** VS. 1.0 - Test version for soil processes without run-off and others
!              percolation dependent of IPERC type
!              REMO calculations in [m], but output is in [mm]/output interval.
!     ***
!     *** VS. 1.1 - Version for REMO
!     ***           Units must be made consistent.

!     *** VS. 1.2 - Limitation to IPERC=2, IDIFF=2
!     ***           --> alternatives removed

!     *** VS. 1.3 - Last loop splitted and loop order changed to get longest loop inside

!     *** VS. 1.4 - September 2007 - HAG
!     ***           Loop from Numerical Recipes changed to get longest loop inside.
!     ***               Therefore were ZTRI(KDEEP) and ZDC changed into fields
!     ***               ZTRI(KLON,KDEEP) and ZDC(KLON,KDEEP).
!     ***           Inconsistency for percolation overflow of layers eliminated
!     ***               now overflow exists if WSI > ZWSFC (instead of ZWSAT as used before)
!     ***           Test for 1E-10. was changed due to the calculation accuracy of the NEC
!     ***               at WSMX-1.E-10 < ZEPS (1E-15)
!     ***
!     ***
!     *** Fluxes are always calculated to be valid two time steps ahead.
!     *** Therefore the fluxes must be divided by 2 for accumulation.
!     *** The calculations themselves uses volumes, not fluxes.
!     *** The flux is calculated so that later the unit of the flux is the accumulated volume
!     *** per time step.
!     ***
!     *** VS. 1.5 - December 2007 - HAG
!     ***           Use of ZWSFC instead of ZWSAT for percolation

!     ***           Juli 2008 - S. HAGEMANN (R. PODZUN)
!     ***           Removal of all occurances of ZVFC in soilhyd.f and surf.f
!     ***           UNSAFE! Use may cause problems. on nec no obvious errors were caused.

!     *** VS. 1.6 - March 2009 - HAG
!     ***           Introduction of maximum percolation rate and solving the percolation
!     ***           equation with midpoint method instead of euler method.

!     *** VS. 1.7 - March-May 2009 - HAG
!     ***           Implementation of van genuchetn equation for percolation.
!     ***           preliminary using iperc to allow switching to previuos
!     ***               Clapp & Hornberger method.
!     ***           Allow ECHAM drainage for lowest layer above bedrock (IDRECH).
!     ***           Allow reduction of hydraulic conductivity with depth for percolation.

!     *** VS. 1.8 - October 2009 - HAG
!     ***           Percolation overflow is not allowed anymore from layers below the surface
!     ***               layer. This has led to a too strong drying
!     ***               of soil layers in wet areas, especially from the second lowest layer
!     ***               above the bedrock. Now it assumed that water piles upwards from
!     ***               SATURATEDsoil5layer.f90 layers, very similar to a rising water table. 
!     ***               Whether this may actually be used to describe a water table, has to be investigated.
!     ***               But this can be also done offline.

!     *** VS. 1.9 - June 2010 - HAG
!     ***            Initialization of zdrain with zero must happen in calling routine before
!     ***                soilchange is called.
!
!     *** vs. 1.10 - Jay 2011 - HAG
!     ***            Code converted using TO_F90 by Alan Miller
!     ***            plus REAL --> REAL(dp) and fine tuning afterwards
!
!     SOILHYD
!     ***  FKSAT = Saturated hydraulic conductivity: unit [m/s] -> Must be checked in ECMGBDT.F
!     ***  FMPOT = Matrix potential [m]
!     *** BCLAPP = Exponent B in Clapp and Hornberger
!     ***   VPOR = Volumetric soil porosity [m/m]
!     ***   SPSI = Soil pore size distribution index used in Van Genuchten

!     *** ZDIAGT = Time step: unit [s]   ##### org. was [day]

!     *** ZDRAIN = Drainage: Volumen during time step: unit [m]
!     ***  ZWSAT = Maximum storage in soil layer I (porosity) [m]
!     ***          Important for the take up of the infiltration in the soil layers
!     ***          Defined down to the depth of the soil - that is where the bedrock begins.
!     ***  ZWSFC = Field capacity for soil layer I (analog to ZWSAT) [m]
!     ***  ZWPWP = Texture depending wilting point of the soil layer I (analog to ZWSAT) [m]
!     ***
!     ***  ZPERC = Downward flux from soil layer I: First in m/day, then in [m].
!     ***          First because of gravitational drainage.
!     ***  ZDIFF = Diffusion of water between soil layer i and layer i+1 [m^2/day]
!     ***   ZMVG = Exponent M in Van Genuchten scheme = SPSI / (SPSI+1)
!     ***   ZDEL = Thickness of the 5 soil layers [m]
!     *** ZKREDU = Reduction factor of hydraulic conductivity with depth for percolation
!     ***  IREDU = Reduction of conductivity with depth yes/no (1/0)
!     ***
!     ***  NSOIL - Number of soil layers
!     *** WSI(I) = Soil moisture content in I. soil layer [m]
!     ***          Updated through the fluxes calculated in SOILHYD
!     ***  ILOG =  1 --> Log output to StdOut with keword SOILFLUX for grid cell JLLOG
!     ***                5 * DIFF. FLUXES, 5 PERC. FLUXES, ZDRAIN
!     ***          2 --> Log output to StdOut with keyword SOILHYD1 for grid cell JLLOG

!     ***  IPERC = Type of percolation
!     ***          1  ECHAM4 drainage (removed, not available)
!     ***          2  Clapp & Hornberger
!     ***          3  Van Genuchten methode 2, after these
!     ***
!     ***  IDRECH = ECHAM drainage
!     ***          0   No, drainage is percolation from lowest layer above bedrock.
!     ***          1   Yes, for lowest layer above bedrock (replaces percolation)
!     ***          2   For all layers above bedrock (in addition to percolation, except
!     ***              for lowest layer above bedrock where it replaces percolation.

USE mo_kind, ONLY: dp
USE mo_land_surface, ONLY: fract_small
IMPLICIT NONE

INTEGER,  INTENT(IN)                     :: kidia
INTEGER,  INTENT(IN)                     :: kfdia
INTEGER,  INTENT(IN)                     :: klon
INTEGER,  INTENT(IN)                     :: kdeep
LOGICAL,  INTENT(IN)                     :: lperma
REAL(dp), INTENT(IN)                     :: zdiagt
INTEGER,  INTENT(IN)                     :: ilog
INTEGER,  INTENT(INOUT)                  :: jllog
REAL(dp), INTENT(INOUT)                  :: wsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: dzsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: zwsat(klon,kdeep)
REAL(dp), INTENT(IN)                     :: zwsfc(klon,kdeep)
REAL(dp), INTENT(IN)                     :: fksat(klon)
REAL(dp), INTENT(IN)                     :: vpor(klon)
REAL(dp), INTENT(IN)                     :: bclapp(klon)
REAL(dp), INTENT(IN)                     :: fmpot(klon)
REAL(dp), INTENT(INOUT)                  :: zdrain(klon)
REAL(dp), INTENT(IN)                     :: spsi(klon)
REAL(dp), INTENT(IN)                     :: zwpwp(klon,kdeep)
REAL(dp), INTENT(IN)                     :: zdel(kdeep)
REAL(dp), INTENT(IN), OPTIONAL           :: ice(klon, kdeep)

REAL(dp), PARAMETER :: zeps   = 1.e-15_dp
INTEGER,  PARAMETER :: iperc  = 3
INTEGER,  PARAMETER :: idrech = 2
INTEGER,  PARAMETER :: iredu  = 0
REAL(dp), PARAMETER :: credu  = 1.0_dp

REAL(dp) :: zdum, zdum2, zdum3
REAL(dp) :: zperc(klon,kdeep), zdum4, zctsat
REAL(dp) :: zmvg(klon), zwrel
REAL(dp) :: zdiflog(10), wslog(kdeep)
!     *** ECHAM-Drainage-Parameter
REAL(dp) :: zdrmin, zdrmax, zdrexp, zwdmin, zwdcri

!     *** Diffusion
REAL(dp) :: zda(klon,kdeep), zdb(klon,kdeep), zdiff(klon,kdeep-1)
REAL(dp) :: ztri(klon,kdeep), zdc(klon,kdeep)

REAL(dp) :: zpmax, zp2         ! Maximum percolation rate
REAL(dp) :: zkredu(kdeep)      ! KSAT reduction factor with depth

!     *** Dummy variables, running indices, etc.
INTEGER  :: i, jl

  zwdmin = 0.05_dp
  zwdcri = 0.9_dp
  zdrmin = 0.001_dp / (3600._dp*1000._dp)
  zdrmax = 0.1_dp / (3600._dp*1000._dp)
  zdrexp = 1.5_dp
  IF (iredu == 1) THEN
    zkredu(1) = (1._dp + EXP(-credu*zdel(1)) )/2._dp
    DO i=2, kdeep
      zkredu(i) = (EXP(-credu*SUM(zdel(1:i-1))) +  &
          EXP(-credu*SUM(zdel(1:i))) )/2._dp
    END DO
  ELSE
    zkredu(1:kdeep) = 1._dp
  END IF

!     ***** Initialization of local fields
  zperc(kidia:kfdia,1:kdeep) = 0._dp
  zmvg(kidia:kfdia) = spsi(kidia:kfdia) / (spsi(kidia:kfdia)+1._dp)

!     *** Percolation = Gravitational drainage - calculated from bottom upward.
!     *** Maximum percoplation due to numeric instability at large time steps
!     *** and small spatial distances.
  IF (iperc == 2) THEN          ! Clapp & Hornberger-formulation
    DO i=1, kdeep
      DO jl=kidia, kfdia
!cc        IF (DZSI(JL,I).GT.0) THEN
        IF (zwsfc(jl,i)-1.e-10_dp > zeps) THEN
          
!         *** MAXIMUM PERCOLATION RATE ZPMAX
          zdum = ( zwsfc(jl,i) -  &
              fksat(jl)*zkredu(i) * zdiagt - zwpwp(jl,i) ) /  &
              ( zwsfc(jl,i) - zwpwp(jl,i))
          IF (zdum > zeps) THEN
            zp2 = fksat(jl)*zkredu(i) * zdum ** (2._dp*bclapp(jl)+3._dp)
          ELSE
            zp2=0._dp
          END IF
          zpmax  = MAX( (fksat(jl)*zkredu(i) + zp2)/2._dp , 0._dp)
          zctsat = MIN( MAX((wsi(jl,i)-zwsfc(jl,i))/zdiagt ,0._dp)  &
              / (fksat(jl)*zkredu(i)) /zdiagt , 1._dp)
          zpmax  = zctsat * fksat(jl)*zkredu(i) + (1._dp-zctsat) * zpmax
          
          zwrel  = MAX(wsi(jl,i)-zwpwp(jl,i), 0._dp) / (zwsfc(jl,i)-zwpwp(jl,i))
          zwrel  = MIN(zwrel, 1._dp)
          
          IF (zwrel-1.e-10_dp > zeps) THEN
            zp2   = fksat(jl)*zkredu(i) * zdiagt / 2._dp * zwrel ** (2._dp*bclapp(jl)+3._dp)
            zdum2 = MAX( MAX(wsi(jl,i)-zp2,0._dp)  -zwpwp(jl,i), 0._dp)/  &
                (zwsfc(jl,i)-zwpwp(jl,i))
            zdum2 = MIN(zdum2, 1._dp)
            zperc(jl,i) = fksat(jl)*zkredu(i) * zdum2 ** (2._dp*bclapp(jl)+3._dp)
            zperc(jl,i) = MIN(zperc(jl,i) , zpmax )
          ELSE
            zperc(jl,i) = 0._dp
          END IF
        END IF
      END DO
    END DO
  ELSE IF (iperc == 3) THEN          ! Van Genuchten-formulation
    DO i = 1, kdeep
      DO jl=kidia, kfdia
        IF (zwsfc(jl,i)-1.e-10_dp > zeps) THEN
          
!         *** MAXIMUM PERCOLATION RATE ZPMAX
          zdum = MAX(zwsfc(jl,i)-zwpwp(jl,i) -  &
              fksat(jl)*zkredu(i)*zdiagt, 0._dp) / (zwsfc(jl,i)-zwpwp(jl,i))
          IF (zdum > zeps) THEN
            zp2 = fksat(jl)*zkredu(i) * SQRT(zdum) * ( 1._dp - ( 1._dp -  &
                zdum**(1._dp/zmvg(jl)) )**zmvg(jl) )**2._dp
          ELSE
            zp2=0._dp
          END IF
          zpmax = MAX( (fksat(jl)*zkredu(i) + zp2)/2._dp , 0._dp)
          zctsat = MIN( MAX((wsi(jl,i)-zwsfc(jl,i))/zdiagt ,0._dp)  &
              / (fksat(jl)*zkredu(i)) /zdiagt , 1._dp)
          zpmax = zctsat * fksat(jl)*zkredu(i) + (1._dp-zctsat) * zpmax
          
!         *** Midpoint-Method
          
          zwrel = MAX(wsi(jl,i)-zwpwp(jl,i), 0._dp) / (zwsfc(jl,i)-zwpwp(jl,i))
          zwrel = MIN(zwrel, 1._dp)
          
          zdum2 = fksat(jl)*zkredu(i) * SQRT(zwrel) * ( 1._dp - ( 1._dp -  &
              zwrel**(1._dp/zmvg(jl)) )**zmvg(jl) )**2._dp
          zdum3 = wsi(jl,i) - MIN(zdum2,zpmax) * zdiagt / 2._dp
          
          zdum4 = MAX(zdum3-zwpwp(jl,i), 0._dp) / (zwsfc(jl,i)-zwpwp(jl,i))
          zdum4 = MIN(zdum4, 1._dp)
          
          zperc(jl,i) = fksat(jl)*zkredu(i) * SQRT(zdum4) * ( 1._dp -  &
              ( 1._dp - zdum4**(1._dp/zmvg(jl)) )**zmvg(jl) )**2._dp
          
          zperc(jl,i) = MIN(zperc(jl,i) , zpmax)
        ELSE
          zperc(jl,i) = 0._dp
        END IF
      END DO
    END DO
  END IF

! *** Using ECHAM4 drainage calculation? for the lowest layer above bedrock
! *** this replaces the percolation value.

  IF (idrech == 1 .OR. idrech == 2) THEN
    DO i = 1, kdeep-1
      DO jl=kidia, kfdia
        IF (dzsi(jl,1) > 0) THEN
          IF (dzsi(jl,i) > 0) THEN
            IF (dzsi(jl,i+1) <= 0) THEN           ! *** Lowest layer above bedrock
              
!cc              ZDUM = ZWSFC(JL,I) * ZWDMIN
              zdum = zwpwp(jl,i)
              zdum2 = zwsfc(jl,i) * zwdcri
              IF (wsi(jl,i) <= zdum) THEN
                zperc(jl,i) = 0._dp
              ELSE
                zdum3 = zdrmin * (wsi(jl,i)-zdum) / (zwsfc(jl,i)-zdum)
                IF (wsi(jl,i) > zdum2) THEN
                  zdum4 = (zdrmax-zdrmin)*( (wsi(jl,i)-zdum2)  &
                      / (zwsat(jl,i) - zdum2) )**zdrexp
                  zdum3 = zdum3 + MIN(zdum4, (wsi(jl,i)-zdum2)/zdiagt)
                END IF
                zperc(jl,i) = MIN(zdum3, (wsi(jl,i)-zdum)/zdiagt)
              END IF
            ELSE IF (idrech == 2) THEN            ! *** Other layers for IDRECH=2
              
!             *** Drainage is calculated with percolation removed from WSI
              zp2 = wsi(jl,i) - zperc(jl,i) * zdiagt
              zdum = zwpwp(jl,i)
              zdum2 = zwsfc(jl,i) * zwdcri
              IF (zp2 > zdum) THEN         ! Drainage occurs only above wilting point
                zdum3 = zdrmin * (zp2-zdum) / (zwsfc(jl,i)-zdum)
                IF (zp2 > zdum2) THEN
                  zdum4 = (zdrmax-zdrmin)*( (zp2-zdum2)  &
                      / (zwsat(jl,i) - zdum2) )**zdrexp
                  zdum3 = zdum3 + MIN(zdum4, (zp2-zdum2)/zdiagt)
                END IF
                
!               *** Except for the lowest layer above bedrock, drainage is removed from
!               *** soil moisture layer wsi before diffusion occurs
!               *** here, drainage is used as accumulated value over ZDIAGT
                zdum4 = MIN(zdum3 * zdiagt, (zp2-zdum))
                wsi(jl,i) = wsi(jl,i) - zdum4
                zdrain(jl) = zdrain(jl) + zdum4
              END IF
            END IF                                 ! Lowest layer / IDRECH=2
          END IF
        END IF
      END DO
    END DO
    
!   *** Lowest layer (=5)
    DO jl=kidia, kfdia
      IF (dzsi(jl,kdeep) > 0) THEN
!cc          ZDUM = ZWSFC(JL,KDEEP) * ZWDMIN
        zdum = zwpwp(jl,kdeep)
        zdum2 = zwsfc(jl,kdeep) * zwdcri
        IF (wsi(jl,kdeep) <= zdum) THEN
          zperc(jl,kdeep) = 0._dp
        ELSE
          zdum3 = zdrmin * (wsi(jl,kdeep)-zdum) / (zwsfc(jl,kdeep)-zdum)
          IF (wsi(jl,kdeep) > zdum2) THEN
            zdum4 = (zdrmax-zdrmin)*( (wsi(jl,kdeep)-zdum2)  &
                / (zwsat(jl,kdeep) - zdum2) )**zdrexp
            zdum3= zdum3 + MIN(zdum4, (wsi(jl,kdeep)-zdum2)/zdiagt)
          END IF
          zperc(jl,kdeep) = MIN(zdum3, (wsi(jl,kdeep)-zdum)/zdiagt)
        END IF
      END IF
    END DO
  END IF

! *** Original unit of percolation fluxes (=unit of FKSAT) was [m/day]
! *** Now it's accumulatied [m] over ZDIAGE (unit of ZDIAGT=[s]
! *** --> Take time step into account : * ZDIAGT
  zperc(kidia:kfdia,1:kdeep) = zperc(kidia:kfdia,1:kdeep) * zdiagt

! *** Diffusion of water *******************************************

! *** Numerical recipes/Richtmyer & -Morton diffusion  *************

  DO i = 1, kdeep-1
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0) THEN
        
!       ***  Soil moisture diffusivity [m^2/day]
!       ***  method c: diffusivity of both levels is averaged
!       ***            to calculate diffusivity between layers
        IF (wsi(jl,i+1) > 0) THEN
          zdum = wsi(jl,i+1)/dzsi(jl,i+1)
          
!         *** Calculating the diffusivity of layer i+1
          IF (lperma) THEN
            zdum4 = zdum / (MAX(0.001_dp, (vpor(jl) - (ice(jl,i+1) / MAX(fract_small,dzsi(jl,i+1))))))
          ELSE
            zdum4 = zdum / vpor(jl)
          ENDIF
          IF (zdum4-1.e-10_dp > zeps) THEN
            zdum2 = bclapp(jl) * fksat(jl) * fmpot(jl) /  &
                zdum * zdum4 ** (bclapp(jl)+3.)
          ELSE
            zdum2 = 0._dp
          END IF
        ELSE
          zdum2=0._dp
        END IF
        IF (wsi(jl,i) > 0._dp) THEN
          zdum = wsi(jl,i)/dzsi(jl,i)
          
!         *** Calculating the diffusivity of layer i
          IF (lperma) THEN
            zdum4 = zdum / (MAX(0.001_dp, (vpor(jl) - (ice(jl,i+1) / MAX(fract_small,dzsi(jl,i+1))))))
          ELSE
            zdum4 = zdum / vpor(jl)
          ENDIF
          IF (zdum4-1.e-10_dp > zeps) THEN
            zdum3 = bclapp(jl) * fksat(jl) * fmpot(jl) /  &
                zdum *  zdum4 ** (bclapp(jl)+3._dp)
          ELSE
            zdum3 = 0._dp
          END IF
        ELSE
          zdum3=0._dp
        END IF
        
!       *** Calculating the diffusivity at the bottom of layer i
        zdiff(jl,i) = (zdum2+zdum3)*0.5_dp
      END IF
    END DO
  END DO

! *** Calculation of diffusion coefficients
! *** Deepest layer KDEEP = 5
  DO jl=kidia, kfdia
    zda(jl, kdeep) = 0._dp
    IF (dzsi(jl,kdeep) > 0._dp) THEN
      
!     *** [ZDIAGT] = s
      zdb(jl, kdeep) = zdiff(jl,kdeep-1) *zdiagt /dzsi(jl,kdeep)  &
          / (dzsi(jl,kdeep)+dzsi(jl,kdeep-1)) * 2._dp
    ELSE
      zdb(jl, kdeep) = 0._dp
    END IF
  END DO

! *** Layers i = 4,3,2
  DO i=kdeep-1, 2,-1
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0) THEN
        IF (dzsi(jl,i+1) > 0._dp) THEN
          zda(jl,i) = zdiff(jl,i) * zdiagt / dzsi(jl,i) /  &
              (dzsi(jl,i)+dzsi(jl,i+1)) * 2._dp
        ELSE
          zda(jl,i) = 0._dp
        END IF
        zdb(jl,i) = zdiff(jl,i-1) * zdiagt / dzsi(jl,i) /  &
            (dzsi(jl,i)+dzsi(jl,i-1)) * 2._dp
      ELSE
        zda(jl, i) = 0._dp
        zdb(jl, i) = 0._dp
      END IF
    END DO
  END DO

! *** Layer 1
  DO jl=kidia, kfdia
    zdb(jl, 1) = 0._dp
    IF (dzsi(jl,1) > 0) THEN
      IF (dzsi(jl,2) > 0) THEN
        zda(jl,1) = zdiff(jl,1) * zdiagt / dzsi(jl,1) /  &
            (dzsi(jl,1)+dzsi(jl,2)) * 2._dp
      ELSE
        zda(jl,1) = 0._dp
      END IF
    ELSE
      zda(jl, 1) = 0._dp
    END IF
  END DO

! *** Temporary storage for log output
  IF (ilog == 1 .OR. ilog == 3) wslog(1:kdeep) = wsi(jllog, 1:kdeep)

! *** Routine TRIDIAG from Numerical Recipes, p. 43
! *** -ZDA = CI, -ZDB=AI, ZDA+ZDB+1=BI, WSI(T)=RI, WSI(T+1)=UI
  DO jl=kidia, kfdia
    IF (dzsi(jl,1) > 0) THEN
      zdc(jl,1) = zda(jl,1)+zdb(jl,1)+1._dp
      wsi(jl,1) = wsi(jl,1) / zdc(jl,1)
    END IF
  END DO
! *** Decomposition and forward substitution
  DO i=2, kdeep
    DO jl=kidia, kfdia
      IF (dzsi(jl,1) > 0) THEN
        
        IF (dzsi(jl,i) > 0) THEN
          ztri(jl,i) = -zda(jl,i-1) / zdc(jl,i-1)
          zdc(jl,i) = zda(jl,i)+zdb(jl,i)+1._dp + zdb(jl,i)*ztri(jl,i)
          wsi(jl,i) = (wsi(jl,i)/dzsi(jl,i) + zdb(jl,i)  &
              * wsi(jl,i-1)/dzsi(jl,i-1) ) / zdc(jl,i) * dzsi(jl,i)
        END IF
      END IF
    END DO
  END DO

  DO i=kdeep-1,1,-1       ! Backsubstitution
    DO jl=kidia, kfdia
      IF (dzsi(jl,1) > 0) THEN
        
        IF (dzsi(jl,i+1) > 0) THEN
          wsi(jl,i) = wsi(jl,i) - ztri(jl,i+1)*wsi(jl,i+1)*  &
              (dzsi(jl,i)/dzsi(jl,i+1))
        END IF
      END IF
    END DO
  END DO

! *** Log output of fluxes ?
  IF (ilog == 1 .OR. ilog == 3) THEN
    DO i=1, kdeep
      zdiflog(i) = wslog(i) - wsi(jllog, i)
    END DO
  END IF

! *** If layer overflow occurs after diffusion, it is added to drainage.
  DO i = 1, kdeep
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0) THEN
        IF (wsi(jl,i) > zwsfc(jl,i)) THEN
          zdrain(jl) = zdrain(jl) + wsi(jl,i) - zwsfc(jl,i)
          wsi(jl,i) = zwsfc(jl,i)
        END IF
      END IF
    END DO
  END DO

  IF (ilog == 2 .OR. ilog == 3) THEN
    WRITE(*,'(A10,1X, 11(E16.9E2, 1X) )')  &
        "SOILHYD1: ", fksat(jllog)*1000._dp, wsi(jllog,1:kdeep)*1000._dp,  &
        zwsfc(jllog,1:kdeep)*1000._dp
  END IF

! ******* NEW SOIL LAYER SOIL MOISTURE CONTENTS *********************

! *** Percolation per time step must be smaller than the usable soil 
! *** layer contents.
  DO i = 1, kdeep
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0 .AND. zperc(jl,i) > 0) THEN
        zdum = zwdmin*zwsfc(jl,i)
        IF (wsi(jl,i)-zperc(jl,i) < zdum) THEN
          IF (wsi(jl,i) > zdum) THEN
            zperc(jl,i) = wsi(jl,i) - zdum
            wsi(jl,i) = zdum
          ELSE
            zperc(jl,i) = 0._dp
          END IF
        ELSE
          wsi(jl,i) = wsi(jl,i) - zperc(jl,i)
        END IF
      END IF
    END DO
  END DO

! *** Percolation added to downward layer or drainage
! *** if lowest layer or bedrock is reached.
  DO i = 1, kdeep-1
    DO jl=kidia, kfdia
      IF (dzsi(jl,1) > 0) THEN
        IF (dzsi(jl,i+1) > 0) THEN
          wsi(jl,i+1) = wsi(jl,i+1) + zperc(jl,i)
        ELSE
          zdrain(jl) = zdrain(jl) + zperc(jl,i)
          zperc(jl,i) = 0._dp
        END IF
      END IF
    END DO
  END DO

  DO jl=kidia, kfdia
    IF (dzsi(jl,1) > 0) THEN
      zdrain(jl) = zdrain(jl) + zperc(jl,kdeep)
    END IF
  END DO

! *** Percolation overflow is not allowed from layers below the surface layer.
! *** Instead it is assumed that water piles upwards from saturated layers.
  DO i = kdeep,2, -1
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0) THEN
        IF (wsi(jl,i) > zwsfc(jl,i)) THEN
          zperc(jl,i-1) = MAX( zperc(jl,i-1) - wsi(jl,i) + zwsfc(jl,i), 0._dp)
          wsi(jl,i-1) = wsi(jl,i-1) + wsi(jl,i) - zwsfc(jl,i)
          wsi(jl,i) = zwsfc(jl,i)
        END IF
      END IF
    END DO
  END DO

! *** Overflow from the surface layer enters drainage. this should not happen,
! *** so this is just for security.
  DO jl=kidia, kfdia
    IF (dzsi(jl,1) > 0) THEN
      IF (wsi(jl,1) > zwsfc(jl,1)) THEN
        zdrain(jl) = zdrain(jl) + wsi(jl,1) - zwsfc(jl,1)
        wsi(jl,1) = zwsfc(jl,1)
      END IF
    END IF
  END DO

! *** Log output of the fluxes ?
  IF (ilog == 1 .OR. ilog == 3) THEN
    DO i=1, kdeep
      zdiflog(kdeep+i) = zperc(jllog,i)
    END DO
    WRITE(*, '(A10,1X,11(E15.6E3,1X) )') "SOILFLUX",zdiflog(1:10)*1000._dp, zdrain(jllog)*1000._dp
  END IF

! *** Log output of WSI!
  IF (ilog == 2 .OR. ilog == 3) THEN
    WRITE(*,'(A10,1X,10(E16.9E2,1X) )' ) "SOILHYD2: ", wsi(jllog,1:kdeep)*1000._dp, dzsi(jllog,1:kdeep)*1000._dp
  END IF

END SUBROUTINE soilhyd

SUBROUTINE soildef (kidia, kfdia, klon, kdeep, lperma, &
    dzr, dzs, vpor, vfc, vpwp, loland, zdel, dzrsi, dzsi, zwsat, zwsfc, zwpwp, ice)
!     ***********************************************************************

!     *** Routine defining various fields for the 5-layer soil scheme.

!     ******** Version 1.0 - Juli 2006
!              Programming and development: Stefan Hagemann
!     ******** Version 1.1 - Juni 2007
!              Optimization : Ralf Podzun
!     ******** Version 1.2 - September 2007
!              A loop optimized : Stefan Hagemann
!     ******** Version 1.3 - Oktober 2007
!              Minor technical correction without effect: Stefan Hagemann
!     ******** Version 1.4 - January 2008
!              Minor technical correction without effect: Stefan Hagemann
!     ******** Version 1.5 - February 2008 - Stefan Hagemann
!              Minor technical correction that may effect only 1-2 gridboxes
!              where DZR=0. (near lakes)
!     ******** Version 1.6 - March 2009 - Stefan Hagemann
!              Volumetric soil field capacity vfc and wilting point vpwp are given
!              via input list and do not need to be computed in soildef.
!              new output: ZWPWP, WSMX eliminated.

!     ***   ZDEL = Thickness of the 5 soil layers [m]
!     ***  DZRSI = Rooted depth per soil layer (until rooting depth DZR) [m]
!     ***   DZSI = Depth that can be saturated with water per soil layer
!     ***          (until bedrock DZS) [m]
!     ***    DZR = Rooting depth in [m]
!     ***    DZS = Soil depth derived from textures in [m]
!     ***   VPOR = Volumetric porosity [m water / m soil]
!     ***    VFC = Volumetric field capacity [m water / m soil]
!     ***   VPWP = Volumetric permanent wilting point [m water / m soil]
!     ***  KDEEP = Number of soil layers
!     ***  ZWSAT = Maximum storage capacity of soil layer I (porosity) [m]
!     ***          Important for the uptake of infiltration in the soil layers.
!     ***          Defined to the depth of the soil, that is, where the bedrock begins.
!     ***  ZWSFC = Field capacity of soil I (analog to ZWSAT) [m]
!     ***  ZWPWP = Texture depending wilting point of soil layer I (analog to ZWSAT) [m]
!
USE mo_kind, ONLY: dp
IMPLICIT NONE
!
INTEGER,  INTENT(IN)                     :: kidia
INTEGER,  INTENT(IN)                     :: kfdia
INTEGER,  INTENT(IN)                     :: klon
INTEGER,  INTENT(IN)                     :: kdeep
LOGICAL,  INTENT(IN)                     :: lperma
REAL(dp), INTENT(IN)                     :: dzr(klon)
REAL(dp), INTENT(IN)                     :: dzs(klon)
REAL(dp), INTENT(IN)                     :: vpor(klon)
REAL(dp), INTENT(IN)                     :: vfc(klon)
REAL(dp), INTENT(IN)                     :: vpwp(klon)
LOGICAL,  INTENT(IN)                     :: loland(klon)
REAL(dp), INTENT(IN)                     :: zdel(kdeep)

REAL(dp), INTENT(IN), OPTIONAL           :: ice(klon,kdeep)

REAL(dp), INTENT(OUT)                    :: dzrsi(klon, kdeep)
REAL(dp), INTENT(OUT)                    :: dzsi(klon, kdeep)
REAL(dp), INTENT(OUT)                    :: zwsat(klon,kdeep)
REAL(dp), INTENT(OUT)                    :: zwsfc(klon,kdeep)
REAL(dp), INTENT(OUT)                    :: zwpwp(klon,kdeep)

REAL(dp), PARAMETER :: zeps = 1.e-15_dp

LOGICAL  :: loflag1(klon), loflag2(klon)

INTEGER  :: i, jl
REAL(dp) :: zdels(kdeep), zdelsm1(kdeep)

! ************ Field initialization -
  DO  jl=1,klon
    loflag1(jl) = .true.
    loflag2(jl) = .true.
  END DO

  dzsi (:,:) = 0._dp
  dzrsi(:,:) = 0._dp
!  CALL setra(dzsi,klon*kdeep,0._dp)
!  CALL setra(dzrsi,klon*kdeep,0._dp)

!CDIR NOVECTOR
  DO  i=2,kdeep
    zdels(i)   = sum(zdel(1:i))
    zdelsm1(i) = sum(zdel(1:i-1))
  END DO

  zdels(1)   = zdel(1)
  zdelsm1(1) = 0._dp

  DO  i=1,kdeep
    DO  jl=kidia, kfdia
      IF (loland(jl).AND.loflag1(jl)) THEN
        IF (dzr(jl) > zdels(i)) THEN
          dzrsi(jl,i) = zdel(i)
        ELSE
          dzrsi(jl,i) = dzr(jl) - zdelsm1(i)
          loflag1(jl) = .false.
        END IF
      END IF
    END DO
  END DO

! *** Saturatable depth
  DO  i=1,kdeep
    DO  jl=kidia, kfdia
      IF (loland(jl) .AND. loflag2(jl)) THEN
        IF (dzs(jl) > zdels(i)) THEN
          dzsi(jl,i) = zdel(i)
        ELSE
          dzsi(jl,i)  = dzs(jl) - zdelsm1(i)
          loflag2(jl) = .false.
        END IF
      END IF
    END DO
  END DO

  IF (lperma) THEN
    DO  i=1,kdeep
      DO  jl=1,klon
        zwsat(jl,i) = MAX(vpor(jl) * dzsi(jl,i) - ice(jl,i),0._dp)
        IF (vpor(jl).GT. zeps) THEN
           zwsfc(jl,i) = zwsat(jl,i) * vfc(jl) / vpor(jl)
           zwpwp(jl,i) = zwsat(jl,i) * vpwp(jl) / vpor(jl)
        ELSE
           zwsfc(jl,i) = 0._dp
           zwpwp(jl,i) = 0._dp
        ENDIF
      END DO
    END DO
  ELSE
    DO  i=1,kdeep
      DO  jl=1,klon
        zwsat(jl,i) = vpor(jl) * dzsi(jl,i)
        zwsfc(jl,i) = vfc (jl) * dzsi(jl,i)
        zwpwp(jl,i) = vpwp(jl) * dzsi(jl,i)
      END DO
    END DO
  ENDIF

END SUBROUTINE soildef


SUBROUTINE soilchange (kidia, kfdia, klon, kdeep,  isch ,  &
    dt, zdiags, loland, &
!         *** - Soil moisture layers,  
    wsi   , wsim1m, dzr, dzrsi  , dzsi, zwsat, zwsfc, &
!         ***   Vertical water fluxes  -- evap. fluxes over land
    finfil, drain, aetrans, aebsoil,                  &
!         ***   Reduction of evaporation if necessary 
    redevap)

!     Intermediate routine for changing 5 soil moisture layers due to
!     evaporation fluxes, infiltration and drainage that are calculated
!     with the old bucket scheme. The soil moisture changes are calculated
!     for one time step, not two as previously done for ws using the
!     3 time step scheme

!     *** VS. 1.0 -- Juli 2006 -- Stefan Hagemann

!     *** VS 1.1 August 2006 -- Stefan Hagemann
!     *** Security check moved to ISCH =1 loop

!     *** VS 1.2 -- January 2007 -- Stefan Hagemann
!     *** IWORK renamed to ISCH .

!     *** VS 1.3 -- September 2007 -- Stefan Hagemann
!     *** Small error in transpiration loop eliminated that occured if transpiration
!     ***   was taken from a layer where WSI < wilting point. The error lead to the
!     ***   the result that the WSI of the corresponding layer was set to the wilting point.
!     ***   The error only causes problems in dry region gridboxes with small WSMX that are
!     ***   close to the bare soil evaporation water depth of 10 cm, and probably
!     ***   only if small time steps are used, e.g. at 0.25 degree.
!     ***   resolution.
!     *** Introduction of zeps as small increment to be used instead of zero in some
!     ***   comparison-if statements. small flux: 0.01 mm/day = 1.16e-10 m/s
!     ***   ZEPS = 1e-15 << against a small flux.
!     *** If WSI < WWILT for one layer, transpiration amount will be reduced
!     ***   instead of re-distribution to other layers.
!     ***   Thus, AETRANS is changed in SOILCHANGE, but also EVAP and EVAPL
!     ***   must be changed.
!     *** Bare soil evaporation will only be taken from the most upper layer.
!     ***   overshooting amounts will be disregarded and AEBSOIL will be reduced.
!     ***   thus, AEBSOIL is changed in SOILCHANGE, but also EVAP and EVAPL
!     ***   must be changed.
!     ***  --> New array submitted to routine: REDEVAP
!     ***
!     *** VS 1.4 -- March 2008 -- Stefan Hagemann
!     *** Reduction of evaporation is taken out again as there was no feedback
!     *** (technical) to the atmosphere. An implementation of this feedback would
!     *** require more complex changes to VDIFF and comprehensive testing. In
!     *** my opinion a reduction of transpiration in VDIFF would be desirable.
!     *** Field REDEVAP is kept for diagnostic purposes.
!     ***
!     *** VS 1.5 - October 2009 Stefan Hagemann
!     ***   Small error in soil moisture change for REDEVAP not zero was corrected.
!     ***   This has led to errors in dry areas.
!     ***   Infiltration fills the layers from above until field capacity. Initially
!     ***   it was filled up until saturation but this caused (likely) unreal(dp)istic
!     ***   high layer overflow values after diffusion(and thus drainage) in dry areas
!     ***   as congo.
!     ***
!     *** VS 1.6 - February 2010 Stefan Hagemann
!     *** Adding infiltration security check if soil for some reason
!     *** cannot take it all

!     ******** Input fields
!     ***AEBSOIL = Bare soil evaporation in [m], accumulated over TWODT/2
!     ***AETRANS = Transpiration in [m], accumulated over TWODT/2
!     ***  DRAIN = Drainage in [m], accumulated over TWODT
!     ***          --> Conversion to TWODT/2 needed
!     *** FINFIL = Infiltration in [m/s] --> multiplication with DT required.

!     ***    DZR = Rooting depth in [m]
!     ***  DZRSI = Rooted depth per soil layer (until rooting depth DZR)
!     ***   DZSI = Depth that can be saturated with water per soil layer
!     ***          (until bedrock DZS)
!     ***  ZWSAT = Maximum storage capacity of soil layer I (porosity) [m]
!     ***          Important for the uptake of infiltration in the soil layers.
!     ***          Defined to the depth of the soil, that is, where the bedrock begins.
!     ***  ZWSFC = Field capacity of soil I (analog to ZWSAT) [m]

!     ***    DT = time step, usually [s]
!     *** ISCH  = Switch which values will change the soil moisture
!     ***         0 = No change: set WSI = WSIM1M
!     ***         1 = change by evaporation fluxes
!     ***         2 = change by infiltration
!     ***         3 = change by infiltration and drainage
USE mo_kind, ONLY: dp
IMPLICIT NONE
!     ---------------------

INTEGER,  INTENT(IN)                     :: kidia
INTEGER,  INTENT(IN)                     :: kfdia
INTEGER,  INTENT(IN)                     :: klon
INTEGER,  INTENT(IN)                     :: kdeep
INTEGER,  INTENT(IN)                     :: isch
REAL(dp), INTENT(IN)                     :: dt
REAL(dp), INTENT(IN)                     :: zdiags
LOGICAL,  INTENT(IN)                     :: loland(klon) ! -- Land and ice mask      --
REAL(dp), INTENT(INOUT)                  :: wsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: wsim1m(klon, kdeep)
REAL(dp), INTENT(IN)                     :: dzr(klon)
REAL(dp), INTENT(IN)                     :: dzrsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: dzsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: zwsat(klon,kdeep)
REAL(dp), INTENT(IN)                     :: zwsfc(klon,kdeep)
REAL(dp), INTENT(IN)                     :: finfil(klon)
REAL(dp), INTENT(INOUT)                  :: drain(klon)
REAL(dp), INTENT(IN)                     :: aetrans(klon)
REAL(dp), INTENT(IN)                     :: aebsoil(klon)
REAL(dp), INTENT(INOUT)                  :: redevap(klon)

REAL(dp), PARAMETER :: zeps = 1.e-15_dp

INTEGER :: jl, jk

REAL(dp) :: zwdmin, zwilt
REAL(dp) :: zrebsoil(klon), zretrans(klon)

REAL(dp) :: zdum, zdum2, zrel
REAL(dp) :: zover(klon)

! *** Constants 
  zwdmin = 0.05_dp
  zwilt  = 0.35_dp ! Read WAVA later on

! ******* No change = set WSI to WSIM1M
  IF (isch == 0) THEN
    wsi(kidia:kfdia,1:kdeep) = wsim1m(kidia:kfdia,1:kdeep)
    RETURN
  END IF

! ******* Change by evaporation fluxes
  IF (isch == 1) THEN
    zrebsoil(kidia:kfdia) = 0._dp
    zretrans(kidia:kfdia) = 0._dp
    zover(kidia:kfdia) = aebsoil(kidia:kfdia)

!   *** Bare soil evaporation loop
    DO  jl=kidia, kfdia
      IF (loland(jl)) THEN
!        *** TAU is packed in SURF in precipitation
!        ***     --> Is contained in the infiltration!
        IF (zover(jl) < -zeps) THEN
          wsi(jl,1) = wsi(jl,1) + zover(jl)
          IF (wsi(jl,1) < 0._dp) THEN
            zover(jl) = wsi(jl,1)
            wsi(jl,1) = 0._dp
          ELSE
            zover(jl)=0._dp
          END IF
        END IF
      END IF
    END DO
    
!       *** Reduction of AEBSOIL necessary? --> diagnostic
    DO jl=kidia, kfdia
      IF (loland(jl)) THEN
        IF (zover(jl) < -zeps) THEN
          zrebsoil(jl) = zover(jl)
!!!            AEBSOIL(JL) = MIN(AEBSOIL(JL) - ZOVER(JL), 0.)
        END IF
      END IF
    END DO

!   *** transpiration loop
    zover(kidia:kfdia) = 0._dp
    DO  jk=1, kdeep
      DO  jl=kidia, kfdia
        IF (loland(jl)) THEN
          IF (aetrans(jl) < -zeps) THEN
            
!           *** Distribute transpiration over rooting depth
!           *** for this, calculate rooted water depth 
            IF (dzr(jl) <= zeps) THEN
              
!             *** Special case: sealed grid cell
              IF (jk == 1) zretrans(jl) = zretrans(jl) + aetrans(jl)
              
            ELSE IF (dzrsi(jl,jk) > 0) THEN
!             *** Note that transpiration should only access water until the
!             *** rooting depth --> pay regard where the relative term of the
!             *** rooted depth within the layer has to be applied
              zrel = MIN(dzrsi(jl,jk)/dzsi(jl,jk), 1._dp)
              
              zdum2 = zwsfc(jl,jk)*zwilt * zrel
              IF (wsi(jl,jk)*zrel >= zdum2) THEN
                zdum = wsi(jl,jk) * zrel + aetrans(jl) * dzrsi(jl,jk)/dzr(jl)
                IF (zdum < zdum2) THEN
                  zover(jl) = zdum - zdum2
                  wsi(jl,jk) = zdum2 + wsi(jl,jk) * MAX(1._dp-zrel, 0._dp)
                ELSE
                  wsi(jl,jk) = zdum + wsi(jl,jk) * MAX(1._dp-zrel, 0._dp)
                  zover(jl)=0._dp
                END IF
              ELSE
                zover(jl) = aetrans(jl) * dzrsi(jl,jk)/dzr(jl)
              END IF
              zretrans(jl) = zretrans(jl) + zover(jl)
            END IF
          END IF
        END IF
      END DO
    END DO

!   *** Adding ZRETRANS and ZREBSOIL to REDEVAP
    DO jl=kidia, kfdia
      IF (loland(jl)) THEN
!CC          IF (zretrans(JL).LT.-zeps) THEN
!CC            aetrans(JL) = MIN(aetrans(JL) - redevap(JL), 0.)
!CC          ENDIF
        redevap(jl) = redevap(jl) + zrebsoil(jl) + zretrans(jl)
      END IF
    END DO
    
!   *** If redevap not zero, soil moisture must be changed
!   *** to close the land surface water balance.
!   *** Reduction will be relatively distributed to the wetness
!   *** of each layer.
    DO  jl=kidia, kfdia
      IF (loland(jl)) THEN
        IF (redevap(jl) < -zeps) THEN
          zdum2 = SUM(wsi(jl,:))
          IF (zdum2 > 0._dp) THEN
            DO jk=1, kdeep
              zdum = wsi(jl,jk) + redevap(jl) * wsi(jl,jk)/zdum2
              wsi(jl,jk) = MAX(zdum, 0._dp)
            END DO
          END IF
        END IF
      END IF
    END DO
  END IF

! ******* Change by infiltration
  IF (isch == 2 .OR. isch == 3) THEN
    zover(kidia:kfdia) = finfil(kidia:kfdia)*dt
    
    DO  jk = 1, kdeep
      DO  jl=kidia, kfdia
        IF (loland(jl)) THEN
          
!         *** infiltration-loop for 5 layer routine

          wsi(jl,jk) = wsi(jl,jk) + zover(jl)
          IF (wsi(jl,jk) > zwsfc(jl,jk)) THEN
            zover(jl)  = wsi(jl,jk) - zwsfc(jl,jk)
            wsi(jl,jk) = zwsfc(jl,jk)
          ELSE
            zover(jl) = 0._dp
          END IF
        END IF
      END DO
    END DO
  
!   *** If infiltration cannot fully enter the soil
    DO jl=kidia, kfdia
      IF (loland(jl)) THEN
        IF (zover(jl) > 0) drain(jl) = drain(jl)+zover(jl)
      END IF
    END DO
  END IF

! ******* Change by drainage (unit m in SURF), used only for testing
  IF (isch == 3) THEN
    DO  jl=kidia, kfdia
      IF (loland(jl).AND.drain(jl) > 0._dp) THEN
        zover(jl) = -drain(jl) * zdiags
        DO jk=kdeep,1,-1
!         *** Start saturated(?) depth
          zdum = zwsfc(jl,jk)*zwdmin
          IF (zwsat(jl,jk) > 0 .AND. wsi(jl,jk)-zdum > 0) THEN
            wsi(jl,jk) =  wsi(jl,jk) + zover(jl)
            IF (wsi(jl,jk)-zdum < 0) THEN
              zover(jl) = wsi(jl,jk)-zdum
              wsi(jl,jk) = zdum
            ELSE
              zover(jl) =0._dp
              EXIT      ! Exit soil layer loop
            END IF
          END IF
        END DO
      END IF
    END DO
  END IF

! ****** Avoid negative values - shouldn't happen
  wsi(kidia:kfdia,1:kdeep) = MAX( wsi(kidia:kfdia,1:kdeep), 0._dp)

END SUBROUTINE soilchange

SUBROUTINE setra ( rar, idim, rwert )

!**** SETRA    -   Up: Set value for a REAL(dp)-field
!**
!**   Call     :   CALL SETRA ( RAR, IDIM, RWERT )
!**   Purpose  :   Fill REAL(dp)-field with given value
!**   Input parameters: 
!**                RAR  -  REAL(dp)-field to fill
!**                IDIM -  Number of elements to fill
!**                RWERT-  REAL(dp)-value to fill RAR with
!**
!**   Author   : Volker Wohlgemuth,F3
!
USE mo_kind, ONLY: dp
!
INTEGER, INTENT(IN)                      :: idim
REAL(dp), INTENT(IN)                     :: rwert
REAL(dp), INTENT(OUT)                    :: rar(idim)

INTEGER ii

  DO  ii = 1,idim
    rar(ii) = rwert
  END DO

END SUBROUTINE setra
