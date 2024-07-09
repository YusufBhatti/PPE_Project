!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_bethy_photosyn
  !
  ! This module contains the photosynthesis/stomata model of BETHY
  ! See W. Knorr, "Satellite Remote Sensing and Modelling of the Global CO2 Exchange of Land Vegetation", 
  !     Examensarbeit Nr. 49 (MPI for Meteorology, Hamburg,1998)
  !

  USE mo_kind,        ONLY: dp
  USE mo_exception,   ONLY: finish

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART =======================================================================================================

  PUBLIC :: photosyn      !! subroutine that models photosynthesis and stomata control

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE ! Make ALL following objects private

  ! --- parameters ----------------------------------------------------------------------------------------------------------------

  ! C3 PLANTS: FARQHUAR, G.D., S. VON CAEMMERER AND J.A. BERRY, 1980.           
  !            A BIOCHEMICAL MODEL OF PHOTOYNTHESIS IN LEAVES OF C3 SPECIES.    
  !            PLANTA 149, 78-90.                                               

  REAL(dp), parameter ::  ALPHA = 0.28_dp    ! EFFICIENCY OF OF PHOTON CAPTURE
  REAL(dp), parameter ::  OX    = 0.21_dp    ! OXYGEN CONCENTRATION [MOL(O2) / MOL(AIR)] 
  REAL(dp), parameter ::  KC0   = 460.E-6_dp ! MICHAELIS-MENTEN CONSTANT FOR CO2 AT 25C [MOL(CO2) / MOL(AIR)] 
  REAL(dp), parameter ::  KO0   = 330.E-3_dp ! MICHAELIS-MENTEN CONSTANT FOR O2 AT 25C [MOL(O2) / MOL(AIR)] 
  REAL(dp), parameter ::  EC    = 59356._dp  ! ACTIVATION ENERGY FOR KC [J / MOL]
  REAL(dp), parameter ::  EO    = 35948._dp  ! ACTIVATION ENERGY FOR KO [J / MOL]
  REAL(dp), parameter ::  EV    = 58520._dp  ! ACTIVATION ENERGY FOR VCMAX [J / MOL]
  REAL(dp), parameter ::  ER    = 45000._dp  ! ACTIVATION ENERGY FOR DARK RESPIRATION [J / MOL]
  REAL(dp), parameter ::  EK    = 50967._dp  !   = Q10=2 (Collatz et al. 1992)
  REAL(dp), parameter ::  FRDC3 = 0.011_dp   ! RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C3

  ! C4 PLANTS: COLLATZ, G.J., M. RIBAS-CARBO AND J.A. BERRY, 1992.              
  !            COUPLED PHOTOSYNTHESIS-STOMATAL CONDUCTANCE MODEL FOR LEAVES     
  !            OF C4 PLANTS. AUST. J. PLANT PHYSIOL. 19, 519-538.               

  REAL(dp), parameter ::  FRDC4 = 0.031_dp   ! RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C4
  REAL(dp), parameter ::  ALC4  = 0.04_dp    ! EFFECTIVE QUANTUM EFFICIENCY
  REAL(dp), parameter ::  THETA = 0.83_dp    ! CURVATURE PARAMETER

  ! --- private fields (for diagnostic output only) ---------


CONTAINS

! --- photosyn() -------------------------------------------------------------------------------------------------------------------
!
! This subroutine operates in two modes: 
!    (i) KFLG = 0: Computes carbon production and stomatal conductance without water stress at given leaf internal CO2-concentration
!   (ii) KFLG = 1: Computes carbon production and leaf internal CO2-concentration at given stomatal conductance
!
! All computations are performed for a single plant functional type and a single canopy layer !!!!
!
! Diagnosis of photorespiration is commented, preliminary and not tested (Thomas Raddatz, Sep. 2006)
! 
  SUBROUTINE photosyn(kidx, VegetationFlag, waterLimitationFlag, CarboxRate, ETransport, C4flag, &
       PAR, TC, P, NSCL, PAR_DOWN, atm_co2_conc, GS, CI, &
       gross_assimilation, dark_respiration, max_carbox_rate, max_e_transport_rate, &
       carbox_rate, e_transport_rate)
    
    USE mo_bethy_constants,   ONLY: minOfMaxCarboxrate, minStomaConductance
    USE mo_jsbach_constants,  ONLY: Tmelt, R_gas => UniversalGasConstant
    
    ! --- interface definition -----------------------------------------------------------------------------------------------------
    INTEGER, INTENT(in)      :: kidx
    logical, intent(in)      :: VegetationFlag (kidx)! Flag is true if Vegetation is present
    logical,  intent(in)     :: waterLimitationFlag   ! FLAG FOR DIFFERENT COMPUTATION MODES:
    !    false: water unlimited: 
    !           Compute photosynthetic production and stomatal conductance 'GS' at
    !           given leaf-internal carbon dioxide concentration 'CI'.
    !     true: water limited:
    !           Compute photosynthetic production and carbon dioxide concentration 
    !           'CI' at given stomatal conductance 'GS'.
    REAL(dp), INTENT(in)     :: CarboxRate(kidx), & ! Maximum carboxilation rate at 25 degrees Celsius [1.E-6 * Mol(CO2)/m^2/s] ...
         ETransport(1:kidx)    ! Maximum rate of electron transport at 25 degC for the PFT handled in this call [MOL(CO2) /m^2 s]
    LOGICAL,  INTENT(in)     :: C4flag(kidx) ! Photosynthetic pathway: C4=.true or C3=.false.
    REAL(dp), intent(in)     :: PAR(kidx)    ! ABSORBED PAR per layer [MOL(PHOTONS) / M^2 S]
    REAL(dp), intent(in)     :: TC(kidx)     ! Canopy (leaf) TEMPERATURE [Celsius]
    REAL(dp), intent(in)     :: P(kidx)        ! AIR PRESSURE [PA]
    REAL(dp), intent(in)     :: NSCL(kidx)     ! Nitrogen scaling factor at maximum carboxylation rate ('Vm') and maximum ..
    ! ..  electron transport rate ('Jm').
    REAL(dp), intent(in)     :: PAR_DOWN(kidx) ! Downward PAR [MOL(PHOTONS) / M^2 S]
    REAL(dp), intent(in)     :: atm_co2_conc(kidx) ! Total irridiance at the surface [mol/(m^2 s)]
    REAL(dp), intent(inout)  :: GS(kidx)     ! STOMATAL CONDUCTANCE FOR WATER VAPOUR [m/s]:
    !    input: for waterLimitationFlag=.true.
    !   output: for waterLimitationFlag=.false.
    REAL(dp), intent(inout)  :: CI(kidx)     ! CO2 concentration inside leaf  [MOL(CO2)/MOL(AIR)]: 
    !    input: for waterLimitationFlag=.false.
    !   output: for waterLimitationFlag=.true.
    REAL(dp), INTENT(out), DIMENSION(kidx), OPTIONAL :: &  ! Output variables
         gross_assimilation,       &                  ! GROSS PHOTOSYNTHESIS [MOL(CO2) / M^2 S]
         dark_respiration,         &                  ! DARK RESPIRATION OF LEAF [MOL(CO2) / M^2 S]
!!$         photo_respiration,        &                  ! PHOTORESPIRATION [MOL(CO2) / M^2 S] 
         max_carbox_rate,          &                  ! Maximum carboxylation rate (=VCmax) [MOL(CO2)/M^2 S]
         max_e_transport_rate,     &                  ! Maximum electron transport rate (=Jmax) (only C3 plants)[MOL(CO2)/M^2 S]
         carbox_rate,              &                  ! Actual carboxylation rate (=JC) [MOL(CO2)/M^2 S]
         e_transport_rate                             ! Actual electron transport rate (=JE) [MOL(CO2)/M^2 S]
    
    ! --- parameters -----------
    
    REAL(dp) ::  T1 =  25._dp + Tmelt                 ! 25 degress Celsius expressed in Kelvin 
    
    ! --- local variables -----------------------------------------------------------------------------------------------------------
    
    REAL(dp) :: VCMAX ! MAXIMUM CARBOXYLATION RATE [MOL(CO2)/M^2 S] 
    REAL(dp) :: GAM   ! COMPENSATION POINT WITHOUT DARK RESPIRATION [MOL(CO2) / MOL(AIR)]
    REAL(dp) :: JMAX  ! MAXIMUM RATE OF ELECTRON TRANSPORT [MOL(CO2)/M^2 S]
    REAL(dp) ::  K    ! CO2 SPECIFICITY OF PECASE  
    REAL(dp) ::  T    ! Canopy Temperature in Kelvin
    
    
    REAL(dp) :: GASS      ! gross assimilation (formerly called A)
    REAL(dp) :: DARK_RESP ! dark respiration (formerly called RD)
!!$    REAL(dp) :: PHOTO_RESP! photo respiration (formerly called RPH)  
    
    REAL(dp) :: KC        ! Michaelis-Menten constant for CO2
    REAL(dp) :: KO        ! Michaelis-Menten constant for O2
    REAL(dp) :: JE        ! Electron transport rate (Knorr, Eq. 102c)
    REAL(dp) :: JC        ! Carboxylation rate (Knorr, Eq. 102b)
    REAL(dp) :: J0, J1
    REAL(dp) :: T0        ! Temperature relative to 25 degree Celcius, i.e. T - 25
    REAL(dp) :: G0 
    REAL(dp) :: K1, W1, K2, W2, B, C ! helpers
    
    INTEGER  :: JL

    !---------------------------------------------------------------------------------
    !                  FIRST ENTRY, TC=TA, CI0 => GC0, AC0                 
    !---------------------------------------------------------------------------------  

    IF (.not. waterLimitationFlag) THEN   ! No water limitation
       
       
       IF (PRESENT(gross_assimilation)) &
            CALL finish('photosyn', 'water-unlimited case should not be called with diagnostic accumulation fields')
       DO JL = 1,kidx
          IF (VegetationFlag(jl)) THEN
             T = TC(JL) + Tmelt          !  Canopy or Vegetation Temperature in Kelvin
             T0 = T - T1                 !  relative T to 25 degree Celcius, means T - 25
             IF (.not. C4flag(JL)) THEN      ! C3
                !---------------------------------------------------------------------------------
                !   C3     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT          
                !
                ! Rate (with Aktivationenergy) vegetation temperature dependance is
                !  k = k(25C) * EXP((Veg Temp - 25) * aktivation energy
                !                    / 298.16 * R * (Veg Temp + 273.16))
                ! => k = k0 * EXP( T0 * E / T1 / R / T ),  WHERE R is the gas constant (8.314)
                ! This holds for OX-Oxygen partial pressure, KC-Michaelis-Menten constant for CO2,
                !    KO-Michaelis-Menten constant for O2, PVM-carboxylation capacity,
                !    RD-Dark respiration, K-PEPcase CO2 specivity
                !    Knorr (106)
                !---------------------------------------------------------------------------------
                KC = KC0 * EXP(EC * T0 / T1 / R_gas / T)
                KO = KO0 * EXP(EO * T0 / T1 / R_gas / T)
                !---------------------------------------------------------------------------------!
                ! CO2 compensation point without leaf respiration, Gamma* is assumed to be linearly 
                ! dependant on vegetation temperature, Gamma* = 1.7 * TC (IF Gamma* in microMol/Mol)
                ! Here, Gam in Mol/Mol,       Knorr (105)
                !---------------------------------------------------------------------------------
                GAM = MAX (1.7E-6_dp * TC(JL), 0._dp)
                !---------------------------------------------------------------------------------
                ! PVM and PJM are not only temperature dependant but also differ inside the canopy. 
                ! This is due to the fact that the plant distributes its Nitrogen content and 
                ! therefore Rubisco content so that, the place with the most incoming light got the 
                ! most Rubisco. Therefore, it is assumed that the Rubisco content falls 
                ! exponentially inside the canopy. This is reflected directly in the values of PVM 
                ! and PJM at 25 Celsius (PVM * nscl),  Knorr (107/108)
                !---------------------------------------------------------------------------------
                VCMAX = CarboxRate(JL) * NSCL(JL) &
                     * EXP(EV * T0 / T1 / R_gas / T)
                !---------------------------------------------------------------------------------
                ! The temperature dependance of the electron transport capacity follows
                ! Farqhuar(1988) with a linear temperature dependance according to the vegetation 
                ! temperature
                !  J = J(25C) * TC / 25 WHERE J(25) = J0 * NSCL
                ! ?????????????????????/ minOfMaxCarboxrate=1E-12 vielleicht etwas gering
                !---------------------------------------------------------------------------------
                JMAX = ETransport(JL) * NSCL(JL)*TC(JL)/25._dp
                JMAX = Max(JMAX,minOfMaxCarboxrate)
            
                !---------------------------------------------------------------------------------
                !                                                                                
                !  C3       GROSS PHOTOSYNTHESIS AT GIVEN CI                                     
                !                                                                                
                !---------------------------------------------------------------------------------
                !c  The assimilation follows the Farqhuar (1980) formulation for C3 plants
                !  A = min{JC, JE} - RD
                !  JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
                !  JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)      with
                !   J = alpha * I * PJM / sqrt(PJM^2 + alpha^2 * I^2) with I=PAR in Mol(Photons)
                !        Knorr (102a-c, 103)
                !  Here J = J1 and A is the gross photosynthesis (GASS), i.e. still including the
                !          respiratory part RD
                !---------------------------------------------------------------------------------
                IF ( JMAX .GT. minOfMaxCarboxrate) THEN
                   J1 = ALPHA * PAR(JL) * JMAX &
                        / SQRT(JMAX**2 + (ALPHA * PAR(JL))**2)
                ELSE
                   J1 = 0._dp
                ENDIF
                JE = J1 * (CI(JL) - GAM) / 4._dp / (CI(JL) + 2._dp * GAM)
                JC = VCMAX * (CI(JL) - GAM) / ( CI(JL) &
                     + KC * (1._dp + OX / KO) )
                GASS = MIN (JE, JC) * HITINHIB(TC(JL))

                !---------------------------------------------------------------------------------
                !                                                                                
                !  C3      COMPUTE 'DARK' RESPIRATION, PHOTORESPIRATION, STOMATAL CONDUCTANCE    
                !                                                                                
                !---------------------------------------------------------------------------------
                !---------------------------------------------------------------------------------
                ! Following Farqhuar et al. (1980), the dark respiration at 25C is proportional to
                ! PVM at 25C, therefore RD = const * PVM, but the temperature dependance goes with
                ! ER (for respiration) and not with EV (for PVM)
                !---------------------------------------------------------------------------------
                DARK_RESP = FRDC3 * CarboxRate(JL) * NSCL(JL) &
                     * EXP(ER * T0 / T1 / R_gas / T) &
                     * HITINHIB(TC(JL)) &
                     * DARKINHIB(PAR_DOWN(JL))
                !---------------------------------------------------------------------------------
                ! Diffusion equation Flux = (CA - CI) / resistence, rs
                !   conductance gs = 1 / rs  =>  Flux = (CA-CI) * gs
                !   (CA ... CO2mixingRatio)
                !   Flux of CO2 is A * amount, Assimilation rate * amount
                !   A is here, Gross Assimilation, though A-RD = (net) Assimilation rate
                !   the amount comes from the ideal gas equation pV=nRT => n/V = p / RT
                !   the stomatal conductance for CO2 is less THEN the conductance of H2O by
                !   the factor of 1.6: gs(CO2) = gs(H2O) / 1.6, due to its lower mobiblity due
                !   to its higher mass
                !   => A (net) = gs/1.6 * (CA-CI) * p/RT
                !   => gs = A(net)*1.6*RT/p/(CA-CI)
                !---------------------------------------------------------------------------------
                
                GS(JL) = MAX (&
                     1.6_dp * (GASS - DARK_RESP) / (Atm_co2_conc(Jl) - CI(JL)) * R_gas * T / P(JL), &
                     minStomaConductance)
 
                !---------------------------------------------------------------------------------
                !
                !                               C 4                                              
                !
                !---------------------------------------------------------------------------------
             ELSE IF (C4flag(JL)) THEN  ! C4
                !---------------------------------------------------------------------------------
                !                                                                                
                !   C4     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT          
                !           AND 'DARK' RESPIRATION                                               
                !                                                                                
                !---------------------------------------------------------------------------------
                ! For C4 plants the Farquhar equations are replaced by the set of equations of
                !  Collatz et al. 1992:
                !  A = min{JC, JE} - RD
                !  JC = k * CI
                !  JE = 1/2/Theta *[PVM + Ji - sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji)]      with
                !  Ji = alphai * Ipar / Epar with Ji=PAR in Mol(Photons)
                !        Knorr (114a-d)
                !  alphai is the integrated quantum efficiency for C4 plants (ALC4 = 0.04, 
                !    compared to the efficiency of C3 plants, ALPHA = 0.28)
                !  Theta is the curve PARAMETER (0.83) which makes the change between
                !   PVM and K limitation smooth
                !  K is the PECase CO2 specifity instead of the electron transport capacity
                !   within C3 plants
                !  Ci is the stomatal CO2 concentration = Cimin + (Ci0 - Cimin)* GC/GC0 with
                !    Cimin = 0.3 and 0.15 CA respectivly (CA is the CO2 mixing ratio)
                !
                ! The factor 1E3 comes that PJM for C3 is in microMol and K is in milliMol,
                !   which is not considered in INITVEGDATA
                ! K scales of course with EK
                !---------------------------------------------------------------------------------
                K = ETransport(JL) * 1.E3_dp * NSCL(JL) &
                     * EXP(EK * T0 / T1 / R_gas / T)
                !---------------------------------------------------------------------------------
                ! same as C3
                !---------------------------------------------------------------------------------
                VCMAX = CarboxRate(JL) &
                     * NSCL(JL) * EXP(EV * T0 / T1 / R_gas / T)
                !---------------------------------------------------------------------------------
                !  same as C3, just the 25 degree Celsius proportional factor is different
                !    0.011 for C3,  0.0042 for C4 
                !---------------------------------------------------------------------------------
                DARK_RESP = FRDC4 * CarboxRate(JL) * NSCL(JL) * EXP(ER * T0 / T1 / R_gas / T) &
                     * HITINHIB(TC(JL)) &
                     * DARKINHIB(PAR_DOWN(JL))
                !---------------------------------------------------------------------------------
                !
                !  C4       GROSS PHOTOSYNTHESIS AT GIVEN CI                                     C
                !
                !---------------------------------------------------------------------------------
                !  JE = 1/2/Theta *[PVM + Ji - sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji)]
                !    Ji = ALC4 * PAR
                !  J0 is the sum of the first two terms in JE
                !---------------------------------------------------------------------------------
                J0 = (ALC4 * PAR(JL) + VCMAX) /  2._dp / THETA
                !---------------------------------------------------------------------------------
                !  last 2 terms:  with J0^2 = 1/4/Theta^2*(PVM+Ji)^2
                !       sqrt(1/4/Theta^2)*sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji))
                !   = sqrt (J0^2 - PVM*Ji/Theta)
                !---------------------------------------------------------------------------------
                JE = J0 - SQRT (J0**2 - VCMAX &
                     * ALC4 * PAR(JL) / THETA)
                !---------------------------------------------------------------------------------
                !         see above
                !---------------------------------------------------------------------------------
                JC = K * CI(JL)
                !---------------------------------------------------------------------------------
                ! same as C3, Farquhar Assimilation
                !---------------------------------------------------------------------------------
                GASS = MIN(JE, JC) * HITINHIB(TC(JL))
                !---------------------------------------------------------------------------------
                !                                                                                
                !   C4     COMPUTE PHOTORESPIRATION, STOMATAL CONDUCTANCE                        
                !
                !---------------------------------------------------------------------------------
                ! same as C3 (diffusion equation)
                !---------------------------------------------------------------------------------
                GS(JL) = 1.6_dp * (GASS - DARK_RESP) / (Atm_co2_conc(Jl) - CI(JL)) * R_gas * T / P(JL)
                GS(JL) = Max(GS(JL), minStomaConductance)
                !---------------------------------------------------------------------------------
                ! Photorespiration is 0 for C4 plants
                !---------------------------------------------------------------------------------
!!$                PHOTO_RESP = 0._dp
           
             ENDIF

          ELSE ! not natural vegetation
             GS(JL) = 0._dp !not defined 
          END IF
       END DO
       !---------------------------------------------------------------------------------      
    ELSE   !  with water limitation (now gs is input and Ci is calculated)
       !---------------------------------------------------------------------------------
          
          IF (.NOT. PRESENT(gross_assimilation)) &
               CALL finish('photosyn', 'water-limited case must be called with diagnostic accumulation fields')
          DO JL = 1,kidx
             IF (VegetationFlag(jl)) THEN
                T = TC(JL) + Tmelt         !  Canopy or Vegetation Temperature in Kelvin
                T0 = T - T1                !  relative T to 25 degree Celcius, means T - 25
                IF (.not. C4flag(JL) ) THEN    ! C3
                   !---------------------------------------------------------------------------------
                   !   C3     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT          
                   !           And 'DARK' RESPIRATION                                               
                   !---------------------------------------------------------------------------------
                   KC = KC0 * EXP(EC * T0 / T1 / R_gas / T)
                   KO = KO0 * EXP(EO * T0 / T1 / R_gas / T)
                   GAM = MAX (1.7E-6_dp * TC(JL), 0._dp)           ! Gam=1.7*TC
                   VCMAX = CarboxRate(JL) * NSCL(JL) &
                        * EXP(EV * T0 / T1 / R_gas / T)           ! PVM=PVM0*NSCL*EXP(Energy scaled relative 25C)
                   DARK_RESP = FRDC3 * CarboxRate(JL) * NSCL(JL)     &  ! RD=RD0*NSCL*EXP(Energy scaled relative 25C)
                        * EXP(ER * T0 / T1 / R_gas / T) &  ! with RD0 proportional PVM0 at 25C
                        * HITINHIB(TC(JL)) &
                        * DARKINHIB(PAR_DOWN(JL))
                   JMAX = ETransport(JL) * NSCL(JL)     & ! PJM=PJM0*NSCL*EXP(Energy scaled relative 25C)
                        * TC(JL) / 25._dp
                   JMAX = Max(JMAX, minOfMaxCarboxrate)
                   !---------------------------------------------------------------------------------
                   !
                   !  C3       GROSS PHOTOSYNTHESIS AT GIVEN TC                                     
                   !
                   !---------------------------------------------------------------------------------
                   !  Remember:
                   !  A = min{JC, JE} - RD
                   !  JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
                   !  JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)      with
                   !   J = alpha * I * PJM / sqrt(PJM^2 + alpha^2 * I^2) with I=PAR in Mol(Photons)
                   !        Knorr (102a-c, 103)
                   ! J = J1
                   !---------------------------------------------------------------------------------
                   IF ( JMAX .GT. minOfMaxCarboxrate) THEN
                      J1 = ALPHA * PAR(JL) * JMAX &
                           / SQRT(JMAX**2 + (ALPHA * PAR(JL))**2)
                   ELSE
                      J1 = 0._dp
                   ENDIF
                   !---------------------------------------------------------------------------------
                   !         Helping friends K1, W1, W2, K2
                   !---------------------------------------------------------------------------------
                   K1 = 2._dp * GAM
                   W1 = J1 / 4._dp
                   W2 = VCMAX
                   K2 = KC * (1._dp + OX / KO)
                   !---------------------------------------------------------------------------------
                   ! A = gs / 1.6 * (CA - Ci) * p / R / T
                   ! <=> Ci = CA - 1.6 * R * T / p / gs * A = CA - A / G0
                   !---------------------------------------------------------------------------------
                   !if (R_gas < 0.1 .OR. T < 0.1 .OR. P(jl) < 0.1) &
                   !print*,jl,R_gas,T,P(jl)
                   G0 = GS(JL) / 1.6_dp / R_gas / T * P(JL)
                   !---------------------------------------------------------------------------------
                   ! A = min{JC, JE} - RD
                   ! => A = JC - RD ^ A = JE - RD
                   ! Set this (A =) in Ci formula above
                   ! Set Ci in
                   !  JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)
                   ! => quadratic formula in JE
                   ! 0 = JE^2 -(RD+J/4+G0*(CA+2*GAMMA))*JE +J/4*G0*(CA-GAMMA)+J/4*RD
                   !---------------------------------------------------------------------------------
                   B = DARK_RESP + W1 + G0 * (Atm_co2_conc(JL) + K1)
                   C = W1 * G0 * (Atm_co2_conc(Jl) - GAM) + W1 * DARK_RESP
                   !---------------------------------------------------------------------------------
                   ! with 0 = x^2 - bx + c
                   !      x1/2 = b/2 +/- sqrt(b^2/4 - c)
                   !  take '-' as minimum value of formula
                   !---------------------------------------------------------------------------------
                   IF ( JMAX .GT. minOfMaxCarboxrate) THEN
                      JE = B / 2._dp - SQRT ( MAX (B**2 / 4._dp - C, 0._dp))
                   ELSE   ! In this case, C = 0
                      JE = 0._dp
                   END IF
                   !---------------------------------------------------------------------------------
                   ! Set Ci in
                   !  JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
                   ! WRITE JC = PVM * (Ci - Gam) / (Ci + K2)
                   ! => quadratic formula in JC
                   ! 0 = JC^2 -(RD+PVM+G0*(CA+K2))*JC +PVM*G0*(CA-GAMMA)+RD*PVM
                   !---------------------------------------------------------------------------------
                   B = DARK_RESP + W2 + G0 * (Atm_co2_conc(Jl) + K2)
                   C = W2 * G0 * (Atm_co2_conc(Jl) - GAM) + W2 * DARK_RESP
                   JC = B / 2._dp - SQRT ( MAX(B**2 / 4._dp - C, 0._dp))
                   !---------------------------------------------------------------------------------
                   ! A = min{JC, JE} - RD
                   !  but here A = Gross Photosynthesis = GPP
                   !  => A = min{JC, JE}
                   !---------------------------------------------------------------------------------
                   GASS = MIN(JE, JC) * HITINHIB(TC(JL))
                   !---------------------------------------------------------------------------------
                   !
                   !   C3     COMPUTE PHOTORESPIRATION AND INTER STOMATE CO2 CONCENTRATION          
                   !
                   !---------------------------------------------------------------------------------
                   ! A = gs / 1.6 * (CA - Ci) * p / R / T
                   ! <=> Ci = CA - 1.6 * R * T / p / gs * A = CA - A / G0
                   !   with A = Net assimilation = NPP
                   ! (CA is the CO2 mixing ratio)
                   !---------------------------------------------------------------------------------
                   CI(JL) = Atm_co2_conc(Jl) - &
                        MAX((GASS - DARK_RESP) / MAX(G0, 1.E-6_dp), 0._dp)
                   !---------------------------------------------------------------------------------
                   ! Photorespiration 
                   ! Carboxylilation controlled assimilation:
                   ! JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
                   ! Photorespiration = PVM * Gam / (Ci + KC * (1 + OX/KO))
                   ! Light limited Assimilation:
                   ! JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)
                   ! Photorespiration = J * Gam / 4 / (Ci + 2 * Gam)
!!$                   PHOTO_RESP = MIN(VCMAX * GAM / ( CI(JL) + KC * ( 1._dp + OX / KO )), &
!!$                        J1 * GAM / 4._dp / (CI(JL) + 2._dp * GAM))                      &
!!$                        * HITINHIB(TC(JL))
                   !---------------------------------------------------------------------------------
                ELSE                          ! C3 -> C4 Plants
                   !---------------------------------------------------------------------------------
                   !
                   !                               C 4                                              
                   !
                   !                                                                                
                   !   C4     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT          
                   !           AND 'DARK' RESPIRATION                                               
                   !
                   !---------------------------------------------------------------------------------
                   ! Remind:
                   !  Collatz et al. 1992:
                   !  A = min{JC, JE} - RD
                   !  JC = k * Ci
                   !  JE = 1/2/Theta *[PVM + Ji - sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji)]      with
                   !   Ji = alphai * Ipar / Epar with Ji=PAR in Mol(Photons)           and
                   !   PAR = Ipar / Epar;  ALC4=alphai; J0=1/2/Theta *(PVM + Ji);
                   !   Ci = CA - 1.6 * R * T / p / gs * A = CA - A / G0                and
                   !   => A = JC - RD ^ A = JE - RD
                   !---------------------------------------------------------------------------------
                   
                   K = ETransport(JL) * 1.E3_dp * NSCL(JL)   &             ! K=K0*NSCL*EXP(Energie scaled relative 25C)
                        * EXP(EK * T0 / T1 / R_gas / T)              !   in mmol CO2 Pecase for C4
                   VCMAX = CarboxRate(JL) * NSCL(JL) &                  ! PVM=PVM0*NSCL*EXP(Energie scaled relative 25C)
                        * EXP(EV * T0 / T1 / R_gas / T)
                   DARK_RESP = FRDC4 * CarboxRate(JL) * NSCL(JL)  &     ! DARK_RESP=RD0*NSCL*EXP(Energie scaled relative 25C)
                        * EXP(ER * T0 / T1 / R_gas / T)   & ! with RD0 proportional PVM0 at 25C
                        * HITINHIB(TC(JL)) &
                        * DARKINHIB(PAR_DOWN(JL))
                   !---------------------------------------------------------------------------------
                   !
                   !  C4       GROSS PHOTOSYNTHESIS AT GIVEN CI                                     
                   !
                   !---------------------------------------------------------------------------------
                   !---------------------------------------------------------------------------------
                   !  Ci = CA - 1.6 * R * T / p / gs * A = CA - A / G0    
                   !--------------------------------------------------------------------------------- 
                   G0 = GS(JL) / 1.6_dp / R_gas / T * P(JL)
                   !---------------------------------------------------------------------------------
                   !  J0=1/2/Theta *(PVM + Ji) = (alphai * PAR + PVM) / 2 / Theta
                   !---------------------------------------------------------------------------------
                   J0 = (ALC4 * PAR(JL) + VCMAX) /  2._dp / THETA
                   !---------------------------------------------------------------------------------
                   !  JE = J0 - sqrt( J0^2 - PVM*alphai*PAR/Theta)
                   !---------------------------------------------------------------------------------
                   JE = J0 - SQRT (J0**2 - VCMAX * ALC4 * PAR(JL) / THETA)
                   !---------------------------------------------------------------------------------
                   !  JC = (G0*CA + RD) / (1 + G0/K)
                   !---------------------------------------------------------------------------------
                   
                   JC = (G0 * Atm_co2_conc(Jl) + DARK_RESP) / (1._dp + G0 / K) 
                   !---------------------------------------------------------------------------------
                   !  A = min{JC, JE} - RD, but here: A = GPP
                   !---------------------------------------------------------------------------------
                   GASS = MIN(JE, JC) * HITINHIB(TC(JL))
                   !---------------------------------------------------------------------------------
                   !   C4     COMPUTE PHOTORESPIRATION AND INTER STOMATE CO2 CONCENTRATION          
                   !---------------------------------------------------------------------------------
                   !  Ci = CA - 1.6 * R * T / p / gs * A = CA - A / G0, with A = NPP
                   !---------------------------------------------------------------------------------
                   CI(JL) = Atm_co2_conc(Jl) - &
                        MAX((GASS-DARK_RESP) / MAX(G0, 1.E-6_dp), 0._dp)
                   !---------------------------------------------------------------------------------
                   ! No Photorespiration with C4 plants
                   !---------------------------------------------------------------------------------
!!$                   PHOTO_RESP = 0._dp
                   JMAX = 0._dp
                ENDIF                !  C3 Plants = 1
                
                !------------------------------------------------------------------------------------------------------------------
                ! Diagnostic output for water limited case (waterLimitationFlag = .true.)
                !------------------------------------------------------------------------------------------------------------------
                gross_assimilation(jl)    = GASS             ! gross assimilation
                dark_respiration(jl)      = DARK_RESP        ! dark respiration
!!$                photo_respiration(jl)     = PHOTO_RESP       ! photo respiration
                max_carbox_rate(jl)       = VCMAX            ! max carboxylation rate
                max_e_transport_rate(jl)  = JMAX             ! max e-transport rate
                carbox_rate(jl)           = JC               ! carboxylation rate
                e_transport_rate(jl)      = JE               ! e-transport rate
                
             ELSE !not natural vegetation
                
                gross_assimilation(jl)    = 0.0_dp           ! gross assimilation
                dark_respiration(jl)      = 0.0_dp           ! dark respiration
!!$                photo_respiration(jl)     = 0.0_dp           ! photo respiration
                max_carbox_rate(jl)       = 0.0_dp           ! max carboxylation rate
                max_e_transport_rate(jl)  = 0.0_dp           ! max e-transport rate
                carbox_rate(jl)           = 0.0_dp           ! carboxylation rate
                e_transport_rate(jl)      = 0.0_dp           ! e-transport rate
             END IF
          END DO                  !  JL
       ENDIF                      !  IPHLG=0

       RETURN
       
     END SUBROUTINE photosyn

!=================================================================================
! Hight temperature inhibition
!=================================================================================
!-----------------------------------------
! FUNCTION which inhibits Assimilation and Respiration at temperatures above 
! 55 Celsius from Collatz et al., Physiological and environmental regulation
! of stomatal conductance, photosynthesis and transpiration: a model that 
! includes a laminar boundary layer,
! Agricultural and Forest Meteorology, 54, pp. 107-136, 1991
!------------------------------------------

!!$REAL(dp) elemental FUNCTION HITINHIB(T)
REAL(dp) FUNCTION HITINHIB(T)
  REAL(dp),intent(in) :: T  !! Canopy temperature in Celsius
  HITINHIB = 1._dp / ( 1._dp + EXP( 1.3_dp * ( T - 55._dp ) ) )
END FUNCTION HITINHIB

!=================================================================================
! Dark inhibition
!=================================================================================
!-----------------------------------------
! FUNCTION which inhibits Dark-Respiration in light
! after Brooks and Farquhar, Effect of temperature on the CO2/O2 specifity on RBisCO
! and the rate of respiration in the light, Planta 165, 397-406, 1985
!
! It is fitted to inhibit the dark-respiration to 50% of it's uninhibited value 
! up from 50 umol/m^2s total solar downward radiation, i.e. 25 umol/m^2s downward PAR
! (idea for improvement: darkinhib should be closer to 1 at very high latitudes)
! The FUNCTION is an own creation.
!---------------------------------------------

!!$REAL(dp) elemental FUNCTION DARKINHIB(PAR_DOWN)
REAL(dp) FUNCTION DARKINHIB(PAR_DOWN)
  REAL(dp),intent(in) :: PAR_DOWN !! downward PAR at the surface [mol phot./(m^2 s)]
  DARKINHIB = 0.5_dp + 0.5_dp * EXP(-2.e5_dp * MAX(PAR_DOWN,0.0_dp))
END FUNCTION DARKINHIB

END module mo_bethy_photosyn
