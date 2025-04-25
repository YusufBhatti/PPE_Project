module mo_cirrus_4

  ! This file is the cirrus box model. It calculates the freezing of aerosol "freezing modes".
  ! This is done by integrating the saturation in timesteps for the length of the simulation time.
  ! It is assuming an adiabatic ascend of the air parcel, this results in an cooling which increases the saturation.
  ! The water vapor can then be deposited on existing ice crystals (frozen freezing modes), or if the critical saturation
  ! of an unfrozen freezing mode is reached, this can freeze.
  !
  ! It can be used as independed box model (see BoxModelExample) or called from a global model (see ModelInterface).
  ! 
  ! It is based on Kärcher et al. (2006), therefore you will find a lot of references to this paper in the code.
  !
  ! In this file the units used units are cm-3 for concentrations and cm for radii.
  ! However, it can be called by the subroutine ModelInterfaceUnitsLiterMicron which takes l-1 and µm as input and output units.
  !
  ! Written in 2016 by Steffen Münch at ETH Zürich
  ! Before it was rewritten, it was used by Miriam Kuebbeler (MK) and Blaz Gasparini (BG). I kept their comments when useful.

  !use mo_kind, only: dp
  
  implicit none
 
  logical :: lCheck = .false.
  logical :: lIgnoreAerosolSizeEffects = .false.
  
  logical :: lDontUseFreezRad = .true.
  logical :: lDustDepRecalcActiveFraction = .true.
  logical :: lActivatedFractionUllrich = .false.
  logical :: lNewDeposition = .true.
  logical :: lVentilation = .true.
  
  integer :: nFreezingModes
  integer :: debugOutput = 0

  !--- Description of freezing modes --------------------------------
  
  type :: FreezingMode
    
    integer :: freezingType = -1
    !-1: freezing mode will be ignored when lFrozen = false (standard unchanged)
    ! 0: preexisting ice; no nucleation but growth by deposition of vapor
    ! 1: homogeneously freezing aerosol
    ! 2: heterogeneous freezing aerosol: Immersion freezing of dust
    ! 3: heterogeneous freezing aerosol: Deposition nucleation
    ! 4: seeding particle
    ! 5: pore condensation freezing, accumulation mode
    ! 6: pore condensation freezing, coarse mode
    ! 7: …(add own here)…
    ! this has to be consistent with CriticalSaturationForFreezing
    
    logical :: lFrozen = .false.
    ! if lFrozen = .false. then nucleation events are calculated that give the iceConc and iceRad after event
    ! before iceConc should be 0.0
    ! if lFrozen = .true. then nucleation events are not calculated and any values in aeroConc are ignored

    real :: aeroConc = 0.0
    real :: aeroRad = 0.0
    real :: aeroStdDev = 0.0
    real :: iceConc = 0.0
    real :: iceRad = 0.0
    
    ! Output: If the mode freezes during cirrus calculation, the time after which if froze is saved here
    real :: freezingTime = -1.0
    real :: depositedWaterVapor = 0.0
    
    ! For recalculate activated fraction
    real :: subtractIce = 0.0
    
    ! Internal: This one is set by the code, it is used to calculate concentration changes due to adiabatic compressions
    real :: tempOfCurrentConc = 0.0

  end type FreezingMode
  
  integer, dimension(3), parameter :: dustDepFreezingTypes = (/3, 5, 6/)
  
  real, parameter :: PI   = 3.1415927
  real, parameter :: iceDens = 925.0 !kg/m3
  ! for deposition
  real, parameter :: depositionCoefficient = 0.5

  contains
    
  !--- Interface functions for global and box model ---------------------------

  subroutine BoxModelExample()

    implicit none

    integer, parameter :: nModes = 1

    real :: simulationTime ! s
    real :: pres ! hPa
    real :: temp ! K
    real :: updraftVelocity ! cm/s
    real :: sat
    type(FreezingMode) :: freezingModes(nModes)

    !lCheck = .false.
    !lIgnoreAerosolSizeEffects = .false.
  
    simulationTime = 30 * 60 ! 30 min
  
    pres = 100 ! hPa
    temp = 190 ! K
    sat = 1.2
    updraftVelocity = 0.01 ! cm/s
  
    freezingModes(1)%freezingType = 1 ! homogeneous
    freezingModes(1)%aeroRad = 0.1 * 1e-6 * 1e2 ! 0.1µm in cm
    freezingModes(1)%aeroConc = 500 ! cm-3
    freezingModes(1)%aeroStdDev = 1
  
    !freezingModes(2)%freezingType = 3 ! heterogeneous
    !freezingModes(2)%aeroRad = 0.5 * 1e-6 * 1e2 ! 0.1µm in cm
    !freezingModes(2)%aeroConc = 0 ! cm-3
    !freezingModes(2)%aeroStdDev = 2
    
    !freezingModes(3)%lFrozen = .true. ! preexisting ice
    !freezingModes(3)%iceRad = 50 * 1e-6 * 1e2 ! cm
    !freezingModes(3)%iceConc = 10e-3 ! cm-3
  
    call ModelInterface(simulationTime, pres, temp, updraftVelocity, sat, freezingModes)
  
    print*, "After calculation..."
    print*, "Pres", pres
    print*, "Temp", temp
    print*, "Sat", sat
    print*, "hom.iceRad", freezingModes(1)%iceRad
    print*, "hom.iceConc", freezingModes(1)%iceConc
    !print*, "het.iceRad", freezingModes(2)%iceRad
    !print*, "het.iceConc", freezingModes(2)%iceConc
    !print*, "pre.iceRad", freezingModes(3)%iceRad
    !print*, "pre.iceConc", freezingModes(3)%iceConc
  
  end subroutine BoxModelExample
  
  subroutine ModelInterface(simulationTime, presIn, tempIn, updraftVelocityIn, sat, freezingModes)

    implicit none

    real, intent(in) :: simulationTime ! s
    real, intent(in) :: presIn ! hPa
    real, intent(in) :: tempIn ! K
    real, intent(in) :: updraftVelocityIn ! cm/s
    real, intent(inout) :: sat
    type(FreezingMode), intent(inout) :: freezingModes(:) ! all sizes and concentrations in cm and cm-3
      
    integer :: i
    
    nFreezingModes = size(freezingModes)
    
    ! set temp of concentrations to start temp; this is used to calculate adiabatic concentration changes
    do i = 1, nFreezingModes
      freezingModes(i)%tempOfCurrentConc = tempIn
    end do
  
    call CalculateCirrus(simulationTime, presIn, tempIn, updraftVelocityIn, sat, freezingModes)
  
    call AdjustConcentrationsToNewConditions(freezingModes, tempIn)
  
  end subroutine ModelInterface
  
  subroutine ModelInterfaceUnitsLiterMicron(simulationTime, pres, temp, updraftVelocity, sat, freezingModes)

    implicit none

    real, intent(in) :: simulationTime ! s
    real, intent(in) :: pres ! hPa
    real, intent(in) :: temp ! K
    real, intent(in) :: updraftVelocity ! cm/s
    real, intent(inout) :: sat
    type(FreezingMode), intent(inout) :: freezingModes(:) ! all sizes in µm and concentrations in l-1
      
    integer :: i
    
    nFreezingModes = size(freezingModes)
        
    ! convert the units to cirrus box model units
    do i = 1, nFreezingModes
      freezingModes(i)%aeroConc = freezingModes(i)%aeroConc * 1e-3 ! l-1 -> cm-3
      freezingModes(i)%aeroRad = freezingModes(i)%aeroRad * 1e-4 ! µm -> cm
      freezingModes(i)%iceConc = freezingModes(i)%iceConc * 1e-3 ! l-1 -> cm-3
      freezingModes(i)%iceRad = freezingModes(i)%iceRad * 1e-4 ! µm -> cm
      
      freezingModes(i)%subtractIce = freezingModes(i)%subtractIce * 1e-3 ! l-1 -> cm-3
    end do
    
    ! now call the actual interface
    call ModelInterface(simulationTime, pres, temp, updraftVelocity, sat, freezingModes)
  
    ! convert the units back
    do i = 1, nFreezingModes
      freezingModes(i)%aeroConc = freezingModes(i)%aeroConc * 1e3 ! cm-3 -> l-1
      freezingModes(i)%aeroRad = freezingModes(i)%aeroRad * 1e4 ! cm -> µm
      freezingModes(i)%iceConc = freezingModes(i)%iceConc * 1e3 ! cm-3 -> l-1
      freezingModes(i)%iceRad = freezingModes(i)%iceRad * 1e4 ! cm -> µm
      
      freezingModes(i)%subtractIce = freezingModes(i)%subtractIce * 1e3 ! cm-3 -> l-1
    end do
  
  end subroutine ModelInterfaceUnitsLiterMicron

  !--- Description of freezing modes: critical saturations --------------------------------

  function CriticalSaturationForFreezing(mode, temp)
    
    ! IMPORTANT LIMITATION: The Critical saturations do not be the same as the model can not freeze two types at ones
    ! If they are the same, one of them freezes first (the one that is the first in the freezing modes array)
    
    implicit none
    
    ! in/out
    type(FreezingMode), intent(in) :: mode
    real, intent(in) :: temp
    real :: CriticalSaturationForFreezing

    if(mode%freezingType == 1) then
      
      ! Critical ice saturation ratio for homogeneous freezing
      ! Evaluated for an intermediate particle size of 0.25 µm
      ! According to: Th.Koop, B.P.Luo, A.Tsias, Th.Peter, Nature 406, 611-614, 2000
      !if( temp < 170.0 .or. temp > 240.0 )then
      !  write(debugOutput,*) &
      !  "Warning: Input temperature not in range 170<temp<240 for homogeneous freezing...Adjusting temp...", temp
      !end if
      CriticalSaturationForFreezing = 2.418 - ( min(240.0, max(170.0, temp)) /245.68)
      
    else if(mode%freezingType == 2) then
      
      ! Immersion freezing of dust (Moehler et al., 2008)
      CriticalSaturationForFreezing = 1.3
      
    else if(mode%freezingType == 3) then
    
      ! Deposition nucleation
      ! In the old code the temp was used, with which the simulation was started, not the current temp
      ! Makes no sense in my oppinion and would need more implementations...
      if(temp >= 220.0) then
        CriticalSaturationForFreezing = 1.2
      else
        CriticalSaturationForFreezing = 1.1
      end if
      
    else if(mode%freezingType == 4) then
    
      ! Seeding particle, use 1.05 for now...
      CriticalSaturationForFreezing = 1.05
      !CriticalSaturationForFreezing = 2.418 - ( min(240.0, max(170.0, temp)) /245.68) - 0.1
      
    else if(mode%freezingType == 5 .or. mode%freezingType == 6) then
  
      ! Pore condensation freezing
      ! Starts already at 100% RHi (see activated fraction)
      CriticalSaturationForFreezing = 1.0
      
    else
      
      ! Unknown freezingType
      CriticalSaturationForFreezing = -1.0
      write(debugOutput,*) "ERROR: Unknown freezingType in CriticalSaturationForFreezing in mo_cirrus_4: ", mode%freezingType
      
    end if
    
  end function CriticalSaturationForFreezing
  
  function ActiveFractionDustDep(tempIn, sat, freezingType, dustRad)
    
    implicit none
    
    ! in/out
    real, intent(in) :: tempIn, sat, dustRad
    integer, intent(in) ::  freezingType
    real :: ActiveFractionDustDep
    real :: temp
    
    real :: a, S0 ! vars for Kuebbeler
    
    ! vars for Ullrich
    real :: ns,surfDust
    real, parameter :: aUll = 285.692 ! paper: alpha
    real, parameter :: bUll = 0.017   ! paper: beta
    real, parameter :: cUll = 256.692 ! paper: gamma
    real, parameter :: dUll = 0.080   ! paper: kappa
    real, parameter :: eUll = 200.745 ! paper: lambda
    
    ! vars for PCF
    integer, parameter :: version = 3
    real, parameter :: A1=-0.0225
    real, parameter :: A2=1.0092
    real, parameter :: x01=122.6554
    real, parameter :: x02=150.7881
    real, parameter :: h1=0.09673
    real, parameter :: h2=0.04779
    real, parameter :: p=0.35949
    
    temp = tempIn

    if(freezingType == 3) then
      if(.not.lActivatedFractionUllrich) then
        
        !MK: calculate active fraction of dust particles following Moehler et al., 2006 (ACP) equation 3
        if(sat >= 1.0 .and. temp <= 238.0) then
          if(temp <= 220.0) then
            a = 2.0
            S0 = 1.1
          else
            a = 0.5
            S0 = 1.2
          end if    
          ActiveFractionDustDep = exp(a * (sat - S0)) - 1.0
          ActiveFractionDustDep = min(ActiveFractionDustDep, 1.0)
          ActiveFractionDustDep = max(ActiveFractionDustDep, 0.0)
        else
          ActiveFractionDustDep = 0.0
        end if
        
      else
        
        ! Activated fraction following Ullrich et al. (2017)
    
        ! only use it for temp > 200 K
        if(temp < 200.0) temp = 200.0
    
        ActiveFractionDustDep = 0.0
        if(sat > 1.0) then

          ! INAS densitiy
          ns = exp( aUll*(sat-1.0)**0.25 * cos(bUll*(temp-cUll))**2 * (PI/2.0 - atan(dUll*(temp-eUll)))/PI ) ! m-2

          ! dust surface in m2
          surfDust = 4. * PI * (dustRad*1e-2)**2

          ActiveFractionDustDep = 1.0 - exp(-ns*surfDust)
          ActiveFractionDustDep = min(ActiveFractionDustDep, 1.0)
          ActiveFractionDustDep = max(ActiveFractionDustDep, 0.0)
      
        end if
        
      end if
    else if(freezingType == 5) then
      ! Pore condensation freezing for accumulation mode after Rob's parameterisation
      if(sat >= 1.0 .and. temp <= 235.0) then
        if(version==1) then
          ActiveFractionDustDep = -0.7/EXP(sat*100/10-10)+0.7
        else if(version==2) then
          ActiveFractionDustDep = -0.7/EXP(sat*100/10.75-10)+0.7
        else if(version==3) then
          ActiveFractionDustDep = A1 + (A2 - A1) * ( p/(1+10**( (x01-sat*100)*h1 )) + (1-p)/(1+10**( (x02-sat*100)*h2 )) )
        end if
      else
        ActiveFractionDustDep = 0.0
      end if
    else if(freezingType == 6) then
      ! Pore condensation freezing for coarse mode after Rob's parameterisation
      if(sat >= 1.0 .and. temp <= 235.0) then
        if(version==1) then
          ActiveFractionDustDep = -1/EXP(sat*100/10-10)+1
        else if(version==2) then
          ActiveFractionDustDep = -1/EXP(sat*100/10.75-10)+1
        else if(version==3) then
          ActiveFractionDustDep = -1/EXP(sat*100/10.75-10)+1
        end if
      else
        ActiveFractionDustDep = 0.0
      end if
    end if
      
    
  end function ActiveFractionDustDep

  !--- Main functions -------------------------------------------------------------------------

  subroutine CalculateCirrus(simulationTimeIn, presStart, tempStart, updraftVelocity, sat, freezingModes)
    
    implicit none
    
    ! input
    real, intent(in) :: simulationTimeIn ! s
    real, intent(in) :: presStart ! hPa
    real, intent(in) :: tempStart ! K
    real, intent(in) :: updraftVelocity ! cm/s
    real, intent(inout) :: sat
    type(FreezingMode), intent(inout) :: freezingModes(nFreezingModes) ! all sizes and concentrations in cm and cm-3

    ! local
    real :: simulationTime
    real :: pres, temp
    real :: satCrit = -1.0 ! random number, always use satCrit in combination with a test for .not.AllModesFrozen(freezingModes) !!!
    integer :: iFrzMode ! for calculating nucleation event
    real :: coolRate
    integer :: i
    logical :: lAllFrozen
    
    ! for stepping back
    real :: stepBackTime
    
    ! for time integration of saturation
    real :: time, dt
    real :: satPres
    real :: a1_Paper, a2_Paper, a3_Paper, b1_Paper, b2_Paper
    ! for deposition on existing ice crystals
    real :: airDens, viscosity, deposOnIce
    real :: fallVel, reynoldsNumber, ventilationCoeff
    real :: depositionOnIceAsDowndraft, effectiveVel
    real :: alpha, vaporAtSatIce, beta
    
    ! for adjusting dt
    logical :: freezeNextTimestep, afterFreezingSmallDt
    real :: dSat, oldDt
    ! for calculating new ice crystal radius
    real :: deltaPlus1, eta
    
    integer :: nTimeSteps
    logical :: lFreezingOccured
    
    ! for deposition
    logical :: lNoIce
    real :: K_thermConduc, D_diffusionCoeff
    real :: dMass_sameForAllIC, dMass_dt, depositedWaterPerIC(nFreezingModes)
    real :: massOld, massNew
    
    ! for recalculating dust activated fraction during ascend
    integer :: nDustDepModes, dustDepIndex(nFreezingModes)
    real :: actFrac, diffIceConc, dustDepAeroVol(nFreezingModes), concNew
    logical :: recalcActFrac
    
    ! new deposition
    !Ls_latHeat = (2834.1 - 0.29*(temp-273.15) - 0.004*(temp-273.15)^2) * 1e3
    real, parameter :: Ls_latHeat = 2839.3*1e3 ! J/kg
    real, parameter :: Rv_specGasConst = 461.9144 ! J/kgK
    
    real, parameter :: relativeSatChangePerTimestep = 1e-2
    
    ! --- Initialization --------------------------------------------------------------------------
    
    ! check input
    if(lCheck) then
      ! just check for unphysical input
      if(simulationTimeIn <= 0.0) then
        write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid simulationTime", simulationTimeIn
        return
      end if
      if(presStart <= 0.0) then
        write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid pressure", presStart
        return
      end if
      if(tempStart <= 0.0) then
        write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid temperature", tempStart
        return
      end if
      if(updraftVelocity <= 0.0) then
        write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid updraftVelocity", updraftVelocity
        return
      end if
      if(sat < 0.0) then
        write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid sat", sat
        return
      end if
      do i = 1, nFreezingModes
        if(.not.freezingModes(i)%lFrozen) then
          if(freezingModes(i)%aeroRad < 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%aeroRad", freezingModes(i)%aeroRad
            return
          end if
          if(freezingModes(i)%aeroConc < 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%aeroConc", freezingModes(i)%aeroConc
            return
          end if
          if(freezingModes(i)%aeroStdDev < 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%aeroStdDev", freezingModes(i)%aeroStdDev
            return
          end if
          if(freezingModes(i)%iceRad .ne. 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%iceRad", freezingModes(i)%iceRad
            return
          end if
          if(freezingModes(i)%iceConc .ne. 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%iceConc", freezingModes(i)%iceConc
            return
          end if
        else
          if(freezingModes(i)%aeroRad .ne. 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%aeroRad", freezingModes(i)%aeroRad
            return
          end if
          if(freezingModes(i)%aeroConc .ne. 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%aeroConc", freezingModes(i)%aeroConc
            return
          end if
          if(freezingModes(i)%aeroStdDev .ne. 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%aeroStdDev", freezingModes(i)%aeroStdDev
            return
          end if
          if(freezingModes(i)%iceRad < 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%iceRad", freezingModes(i)%iceRad
            return
          end if
          if(freezingModes(i)%iceConc < 0.0) then
            write(debugOutput,*) "ERROR in mo_cirrus_4: Invalid freezingModes(",i,")%iceConc", freezingModes(i)%iceConc
            return
          end if
        end if
        
      end do
    end if
    
    simulationTime = simulationTimeIn
    pres = presStart
    temp = tempStart
    
    ! set temp of concentrations to start temp; this is used to calculate adiabatic concentration changes
    ! already done above but for safety reasons here again...
    do i = 1, nFreezingModes
      freezingModes(i)%tempOfCurrentConc = temp
      freezingModes(i)%depositedWaterVapor = 0.0
    end do
    
    coolRate = 9.80999976e-05 * updraftVelocity ! number is dryAdiabLapseRate
    
    lAllFrozen = AllFrozen(freezingModes)
    lNoIce = .true.
    do i = 1, nFreezingModes
      if(freezingModes(i)%iceConc > 0.0) lNoIce = .false.
    end do
    
    if(lDustDepRecalcActiveFraction) then
      call GetDustDepIndex(freezingModes,nDustDepModes,dustDepIndex)    
      if(nDustDepModes > 0) then
        do i = 1, nDustDepModes
          dustDepAeroVol(i) = 4./3.*PI*(freezingModes(dustDepIndex(i))%aeroRad)**3
        end do
      end if
    end if
    
    ! --- Go back in time? --------------------------------------------------------------------------
    
    ! Check if sat is already too high so that freezing would have occured before
    ! If so, step back in time (desend air parcel) so that the nucleation is included in the simulation
    ! Should this be done in global-model-case???
    if(.false.) then
      if(.not.lAllFrozen) then
        call NextCriticalSaturation(freezingModes, temp, satCrit, iFrzMode)
      
        if(sat > satCrit) then
        
          ! go to the point where sat was 95% of satCrit
          stepBackTime = 1.0e5 * log(sat/(0.95*satCrit)) / updraftVelocity
          sat = 0.95 * satCrit
          simulationTime = simulationTime + stepBackTime ! stepping back is the same as increasing the simulation time
      
          temp = temp + 1.0e-4 * updraftVelocity * stepBackTime
          pres = presStart * (temp/tempStart)**3.5
      
          call AdjustConcentrationsToNewConditions(freezingModes, temp)
                  
        end if
      
      end if
    end if
    
    ! --- Are we already above satCrit? ------------------------------------------------------
    
    ! Maybe implement some check here that freezes INP if we already are above the critical saturation
    ! Right now this just happens in the first step of the time loop
    
    ! --- Time loop --------------------------------------------------------------------------
    
    ! init time loop
    nTimeSteps = 0
    time = 0.0
    dt = 0.0005 * (1.0e5/updraftVelocity) ! not needed for new algorithm
    dt = min(dt, simulationTime)
    freezeNextTimestep = .false.
    afterFreezingSmallDt = .false.
    recalcActFrac = .true.
    lFreezingOccured = .false.
    dMass_sameForAllIC = 0.0
    
    ! do time integration
    do while(time < simulationTime)
      
      ! get next freezing saturation
      if(.not.lAllFrozen) call NextCriticalSaturation(freezingModes, temp, satCrit, iFrzMode)
      
      ! --- Calculate deposition --------------------------------------------------------------------------
      
      ! a1 is used for deposition and the calculation of dt
      ! see Kärcher et al. 2006 for the formulas of a1,a2,a3,b1,b2
      a1_Paper = ( 0.601272523 / temp - 0.000342181855 ) / temp
      
      depositionOnIceAsDowndraft = 0.0
      if(.not.lNoIce) then
        ! the next block (deposition) is necessary to calculate the time derivation of sat dSat
        ! this is then used to calculate a timestep with an euler forward algorithm
        ! for the euler forward, the derivation has to be calculated at the current time (a1,a2,..., and deposition at current time)
        ! therefore just the current temp, pres, sat, iceConc,... can be used

        ! the next block (deposition) is also necessary for the calculation of a freezing event
        ! so for both cases these calculations have to be done first
      
        ! calculate deposition on existing ice crystals
        airDens = 0.35 * pres / temp
        viscosity = 1e-5 * ( 1.512 + 0.0052*(temp-233.15) )
        if(.not.lNewDeposition) then
        
          call CalculateA1A2A3B1B2(temp, pres, sat, a1_Paper, a2_Paper, a3_Paper, b1_Paper, b2_Paper)

          deposOnIce = 0.0
        
          do i = 1, nFreezingModes
            if(freezingModes(i)%lFrozen .and. freezingModes(i)%iceConc > 0.0) then
          
              ventilationCoeff = 1.0
              if(lVentilation) then
                fallVel = 14.0 * freezingModes(i)%iceRad * (1.3/airDens)**0.35
                ! BG: What kind of fall velocity parametrization is that??????
                !     not mentioned in article neither in Pruppacher and Klett
                reynoldsNumber = airDens * 2.0 * freezingModes(i)%iceRad*1e-2 * fallVel / viscosity
                ventilationCoeff = 1.0 + 0.229 * sqrt(reynoldsNumber)
                ! BG: ventilation coeff for preexisting ice particles: because they may have grown > 25 microns, VENT>1
              end if
            
              deposOnIce = deposOnIce + 3.89051704e23 * freezingModes(i)%iceConc * ventilationCoeff & ! number is fourPI_SpecVol
                * b1_Paper * freezingModes(i)%iceRad**2.0 / ( 1.0 + b2_Paper * freezingModes(i)%iceRad )
                ! according to Kärcher et al. (2006)

            end if
          end do
        
          ! deposition can be formulated as a downdraft which decreases the updraftVelocity
          depositionOnIceAsDowndraft = ( a2_Paper + a3_Paper * sat ) * deposOnIce / ( a1_Paper * sat )
        
        else ! new deposition
        
          if(.not.lFreezingOccured .or. dMass_sameForAllIC == 0.0) then
            ! if this condition applies, temp and pres changed and the following values have to be recalculated
      
            K_thermConduc = 4.1868e-3 * (5.69 + 0.017*(temp-273.15)) ! W/mK
            D_diffusionCoeff = 2.11e-5 * (temp/273.15)**1.94 * (101325 / (pres*1e2)) ! m2/s
            satPres = 3.6e10 * exp(-(6145.0/temp))

            dMass_sameForAllIC = 4.0 * PI * (sat - 1.0) * depositionCoefficient / &
              ( Ls_latHeat/(K_thermConduc*temp) * ( Ls_latHeat/(Rv_specGasConst*temp) - 1.0 ) &
                + Rv_specGasConst*temp/( D_diffusionCoeff *  satPres*1e2) ) ! kg per m, s, and per IC
                ! note, that to get kg/s this has to be multiplied by the IC capacitance
                ! this capacitance depends on the IC shape and for a sphere is just the radius
                ! however for plates it would be 2r/PI
          end if
          
          dMass_dt = 0.0
          depositedWaterPerIC = 0.0
          do i = 1, nFreezingModes
            if(freezingModes(i)%lFrozen .and. freezingModes(i)%iceConc > 0.0 .and. freezingModes(i)%iceRad > 0.0) then
          
              ventilationCoeff = 1.0
              if(lVentilation) then
                fallVel = 14.0 * freezingModes(i)%iceRad * (1.3/airDens)**0.35
                ! BG: What kind of fall velocity parametrization is that??????
                !     not mentioned in article neither in Pruppacher and Klett
                reynoldsNumber = airDens * 2.0 * freezingModes(i)%iceRad*1e-2 * fallVel / viscosity
                ventilationCoeff = 1.0 + 0.229 * sqrt(reynoldsNumber)
                ! BG: ventilation coeff for preexisting ice particles: because they may have grown > 25 microns, VENT>1
              end if
          
              depositedWaterPerIC(i) = dMass_sameForAllIC * freezingModes(i)%iceRad*1e-2 * ventilationCoeff ! kg/s
          
              dMass_dt = dMass_dt + depositedWaterPerIC(i) * freezingModes(i)%iceConc*1e6 ! kg/(m3 s)
          
            end if
          end do
        
          ! deposition can be formulated as a downdraft which decreases the updraftVelocity
          
          ! if(satPres > 0.0 .and. sat > 0.0) then
          !   depositionOnIceAsDowndraft = dMass_dt * 287.0 * temp / (satPres*1e2) / (a1_Paper * sat)
          ! else
          !   depositedWaterPerIC = 0.0
          !   depositionOnIceAsDowndraft = 0.0
          ! end if

          alpha = 59.7839 / (temp**2) - 0.03417642 / temp
          vaporAtSatIce = 0.62 * 611.2 * exp(22.46*(temp-273.15) / (272.62 + temp - 273.15)) / pres
          beta = 1.0/vaporAtSatIce + 17273953.0 / (temp**2)

          depositionOnIceAsDowndraft = beta/alpha * dMass_dt/airDens * 1e2 ! m/s -> cm/s
          
        end if
      end if
      
      effectiveVel = updraftVelocity - depositionOnIceAsDowndraft
      
      ! --- Freezing event --------------------------------------------------------------------------
      
      lFreezingOccured = .true.
      
      ! see if freezing event has to be calculated
      if(.not.lAllFrozen .and. (sat >= satCrit .or. freezeNextTimestep)) then
        
        freezeNextTimestep = .false.
        
        ! The freezing parameterization only works for positive velocities
        if(effectiveVel>0.0) then
        
          ! check input
          if(freezingModes(iFrzMode)%aeroConc > 0.0 .and. freezingModes(iFrzMode)%aeroRad > 0.0 &
            .and. freezingModes(iFrzMode)%aeroStdDev > 0.0) then
          
            ! calculate freezing of that mode: done in seperate subroutine
            call CalculateFreezingEvent(pres, temp, effectiveVel, sat, freezingModes(iFrzMode))
              
            ! Save time of freezing
            freezingModes(iFrzMode)%freezingTime = time
            
            lNoIce = .false.
            
            ! After a freezing event: do shorter timesteps starting with 1 second
            afterFreezingSmallDt = .true.
            dt = 1
          
          else
          
            ! Just set the mode to frozen so that it does not disturb the ongoing calculations
            freezingModes(iFrzMode)%lFrozen = .true.
            freezingModes(iFrzMode)%iceConc = 0.0
            freezingModes(iFrzMode)%iceRad = 0.0
            
          end if
          
          ! Recalculate everything after freezing event as ice crystals have changed
          lAllFrozen = AllFrozen(freezingModes)
          cycle
          
        else !if(effectiveVel>0.0) then
          
          ! Unfortunately as the freezing parameterization only works for positive velocities there is nothing
          ! that can be done right now.
          ! Just wait and see if we get positive velocity later...
          
        end if
        
      end if
      
      ! --- Freeze activated fraction of dustDep --------------------------------------------------------------------------
      
      ! if lDustDepRecalcActiveFraction then dust deposition is not calculated above
      if(lDustDepRecalcActiveFraction .and. recalcActFrac .and. nDustDepModes > 0) then
        
        if((sat > 1.1 .and. temp <= 220.0) .or. (sat > 1.2 .and. temp > 220.0)) then
          
          recalcActFrac = .true.
          
          do i=1,nDustDepModes
          
            actFrac = ActiveFractionDustDep(temp, sat, freezingModes(dustDepIndex(i))%freezingType, &
              freezingModes(dustDepIndex(i))%aeroRad)
            diffIceConc = freezingModes(dustDepIndex(i))%aeroConc*actFrac - freezingModes(dustDepIndex(i))%iceConc &
              - freezingModes(dustDepIndex(i))%subtractIce
          
            if(diffIceConc > 0.0) then
                        
              massOld = 4./3.*PI*(freezingModes(dustDepIndex(i))%iceRad)**3 * & ! actally this is volume not mass
                freezingModes(dustDepIndex(i))%iceConc
              massNew = dustDepAeroVol(i) * diffIceConc
              concNew = freezingModes(dustDepIndex(i))%iceConc + diffIceConc
              freezingModes(dustDepIndex(i))%iceRad = ((massOld+massNew) / &
                concNew / (4./3.*PI))**(1./3.)
              
              freezingModes(dustDepIndex(i))%iceConc = concNew
            
              freezingModes(dustDepIndex(i))%lFrozen = .true.
            
              lNoIce = .false.
            
              recalcActFrac = .false.
            
            end if
            
          end do
          
          if(.not.recalcActFrac) cycle
          
        end if
        
      end if
      recalcActFrac = .true.
      
      lFreezingOccured = .false.
      
      ! --- Time integration: find timestep length --------------------------------------------------------------------------
      
      !if(lAllFrozen) effectiveVel = effectiveVel - updraftVelocity ! no further cooling ?
            
      ! new timestep algorithm, the timestep length can be calculated to get a desired satRelativeChange
      ! formula for satRelativeChange: satRelativeChange = abs(a1_Paper * ( effectiveVel ) * dt)
      ! this is solved for dt to get the timestep length
      oldDt = dt
      if(effectiveVel /= 0.0) then
        dt = relativeSatChangePerTimestep / (a1_Paper * abs( effectiveVel ))
        
        ! if there just was a freezing event: do smaller timesteps starting with 1 s
        ! but increase dt every step by a factor of 2
        if(afterFreezingSmallDt) then
          if(oldDt*2<dt) then
            dt = oldDt * 2
          else
            afterFreezingSmallDt = .false.
          end if
        end if
        
      else
        dt = simulationTime * 1e-4
      end if
      
      ! check if timestep is too short -> speed it up a little
      if(dt < (simulationTime*1e-4)) dt = simulationTime * 1e-4
      
      ! include a little mechanism that the timestep does not become too big?
      ! without that the steadystate criteria will never work
      ! however, maybe it is not necessary...
      ! if you want to speed up the model code - remove this if clause
      !if(dt > simulationTime * 1e-3) dt = simulationTime * 1e-3
      
      ! check if timestep is too big
      if(time+dt > simulationTime) dt = simulationTime - time
      
      ! now calculate saturation change
      dSat = a1_Paper * sat * ( effectiveVel ) * dt
            
      ! check if timestep is too big so that it exceeds satCrit
      if(.not.lAllFrozen .and. sat+dSat > satCrit) then
        ! This can also occur in some strange situations when satCrit changed between timesteps (e.g. dustDep < 220 K)
        ! and freezing can not be calculated as effectiveVel <= 0.0.
        ! But in this special situation dt should not be adjusted, as the calculated dt would be wrong.
        ! Therefore only:
        if(effectiveVel > 0.0 .and. satCrit > sat) then
          dSat = satCrit - sat
          dt = dSat / (a1_Paper * sat * effectiveVel )
          ! Unfortunately this can be not 100% accurate as satCrit can change in next timestep (e.g. homo is temp dependend).
          ! Therefore use additional switch:
          freezeNextTimestep = .true.
        end if
      end if
      
      ! --- Time integration: calculate values for next timestep -----------------------------------------------
      
      time = time + dt
      sat = sat + dSat
      temp  = temp - dt * coolRate
      !if(lAllFrozen) temp  = temp + dt * coolRate ! no further cooling?
      pres = presStart * ( temp / tempStart )**3.5
      ! Concentration change due to adiabatic compression has to be calculated:
      call AdjustConcentrationsToNewConditions(freezingModes, temp) 
      
      ! calc ice crystal growth by water vapor deposition -> new ice radius:
      if(.not.lNoIce) then
        if(.not.lNewDeposition) then
        
          do i = 1, nFreezingModes
            if(freezingModes(i)%lFrozen .and. freezingModes(i)%iceConc > 0.0 .and. freezingModes(i)%iceRad > 0.0) then
          
              ventilationCoeff = 1.0
              if(lVentilation) then
                fallVel = 14.0* freezingModes(i)%iceRad * (1.3/airDens)**0.35
                ! BG: What kind of fall velocity parametrization is that??????
                !     not mentioned in article neither in Pruppacher and Klett
                reynoldsNumber = airDens * 2.0 * freezingModes(i)%iceRad*1e-2 * fallVel / viscosity
                ventilationCoeff = 1.0 + 0.229 * sqrt(reynoldsNumber)
                ! BG: ventilation coeff for preexisting ice particles: because they may have grown > 25 microns, VENT>1
              end if
            
              deltaPlus1 = b2_Paper * freezingModes(i)%iceRad + 1.0
              eta = 1.0 + 2.0 * ventilationCoeff * b1_Paper * b2_Paper * dt / deltaPlus1**2.0
              if(eta >= 0.0) then
                freezingModes(i)%iceRad = ( deltaPlus1 * sqrt(eta) - 1.0 ) / b2_Paper ! using eq 18, Kärcher et al. 2006
              else
                ! I guess in this case the particles evaporate BUT I HAVE NOT CHECKED THAT YET !!!
                ! I think it only happens at saturations close to 0.0
                freezingModes(i)%iceRad = 0.0
                freezingModes(i)%iceConc = 0.0
              end if
          
            end if
          end do
        
        else ! new deposition
        
          ! Calculate deposition like in the books Lohmann et al. (2016) and Pruppracher and Klett (2010)
        
          do i = 1, nFreezingModes
            if(depositedWaterPerIC(i) > 0.0) then
            
              dMass_dt = depositedWaterPerIC(i) ! kg/s
              massOld = 4./3.*PI*(freezingModes(i)%iceRad*1e-2)**3 * iceDens
              massNew = massOld + dMass_dt * dt
              freezingModes(i)%iceRad = (massNew / (4./3.*PI * iceDens))**(1./3.) * 1e2
            
              freezingModes(i)%depositedWaterVapor = freezingModes(i)%depositedWaterVapor &
                + dMass_dt * dt * freezingModes(i)%iceConc*1e6 ! kg/m3
            
            end if
          end do
        
        end if ! deposition
      end if
      
      ! This still occurs (sat very rarly; temp gets small at very high updrafts). Fix this later!
      if(sat <= 0.0 .or. temp <= 100.0) then
        !write(debugOutput,*) "sat <= 0.0", time, simulationTime, dt, sat, dSat, temp, effectiveVel
        exit
      end if
      ! acctually with the minimum time step of simulationTime * 1e-4 this should never happen
      ! I leave it here anyway as a security check to avoid an endless loop 
      nTimeSteps = nTimeSteps + 1
      if(nTimeSteps>10500) then
        !write(debugOutput,*) "Cirrus scheme: 10500 time steps reached...this should not happen",&
        !"(time, simulationTime)", time, simulationTime
        exit
      end if
      
    end do
    
  end subroutine CalculateCirrus

  subroutine CalculateFreezingEvent(pres, temp, effectiveVelocity, sat, frzMode)

    ! Right now this calculates the new iceConc, iceRad, aeroConc, and aeroRad
    ! sat is not changed, also the calculation of aeroRad is not that good...
    !
    ! Most of this was taken from the old code (with some restructuring).
    !
    ! At the moment there are three freezing parameterizations:
    ! 1. size effects are completely ignored, all aerosols are assumed to be 0.25 µm, according to Kärcher and Lohmann 2002a
    ! 2. size effects are considered, but all aerosols are assumed to have the same size, according to Kärcher and Lohmann 2002b
    ! 3. full size effects are considered, aerosols are sorted into size bins and then aerosols freeze, first the biggest ones
    !    then going down the bins, until all available water is frozen, thereby not all aerosol have to freeze

    real, intent(inout) :: pres ! hPa
    real, intent(inout) :: temp ! K
    real, intent(in) :: effectiveVelocity ! cm/s
    real, intent(inout) :: sat
    type(FreezingMode), intent(inout) :: frzMode
      
    real :: a1_Paper, a2_Paper, a3_Paper, b1_Paper, b2_Paper
    real :: coolRate
    real :: freezingParameterization_C, freezingTimeScale
    real :: availableWaterForNucleation
    
    real :: kappa, delta, deltaPlus1, TBBT
    real :: errorFuncApprox, growthTermFactors, R_growthTerm_monoDisperse
    
    real :: minIceRad
    real :: newIceRad, newIceMass, newIceConc
    
    integer :: i

    real, parameter :: SVOL = 3.23e-23 ! specVolumeWater
    !real, parameter :: PI   = 3.1415927
    real, parameter :: TWOPI = 6.2831853 ! 2 * PI
    real, parameter :: XMW = 2.992e-23 ! some Mass of water
    real, parameter :: RHOICE = 0.925
    
    real, parameter :: VRAT = 2.0
    real, parameter :: THIRD = 1.0/3.0
    
    integer, parameter :: nBinsMax = 41
    integer :: nBins, lastFrozenBin
    real, parameter :: radMax = 1.E-3
    real, parameter :: radMin = 1.E-7
    real :: sumFrozenAero, sumDeposWaterMole, nAllAeroInBins
    real :: sumFrozenAeroLast, sumDeposWaterMoleLast
    real :: DLOGR0, SIGL, SHAPFC1, SHAPFC2, ARG
    real :: binsRad(nBinsMax), nAeroInBin(nBinsMax)
    real :: radMean, RIRAT, RLOGRAT, SLOPE, SLOPER
    
    ! --- Initialization --------------------------------------------------------------------------
    
    call CalculateA1A2A3B1B2(temp, pres, sat, a1_Paper, a2_Paper, a3_Paper, b1_Paper, b2_Paper)
    
    if(frzMode%freezingType>1) then
      ! in case of heterogeneous nucleation b1_Paper is calculated with a slightly higher saturation
      ! don't know why at the moment...
      !b1_Paper = waterToCrystalFlux * 3.23e-23 * vaporConcAtSatPres * (sat-1.0) ! number is specVolumeWater
      b1_Paper = b1_Paper / (sat-1.0) * (sat-1.0+0.05)
    end if
    
    coolRate = -9.80999976e-05 * effectiveVelocity ! number is dryAdiabLapseRate
    freezingParameterization_C = -0.004*temp**2 + 2.0*temp - 304.4
    freezingTimeScale = 1.0 / ( freezingParameterization_C * coolRate )
    
    availableWaterForNucleation = effectiveVelocity * a1_Paper*sat / ( a2_Paper + a3_Paper*sat )
    
    ! --- Freezing parameterization type 1: all size effects are ignored ---------------------------------------
    if(lIgnoreAerosolSizeEffects .and. frzMode%freezingType==1) then
      ! according to Kärcher and Lohmann 2002a
      ! all crystals are assumed to be 0.25 µm for the calculations, only at the end the radius is calculated from the mass
      minIceRad = 0.25e-4 ! 0.25 µm
      kappa = freezingTimeScale * ( b1_Paper/minIceRad ) / ( 1.0 + b2_Paper*minIceRad )
      newIceConc = SVOL * ( b2_Paper / (TWOPI*b1_Paper) )**1.5 * availableWaterForNucleation / sqrt(freezingTimeScale)
      newIceConc = min(newIceConc, frzMode%aeroConc)
      newIceMass = XMW * PI * availableWaterForNucleation * freezingTimeScale / 6.0
      newIceRad = ( newIceMass / (PI * RHOICE / 0.75 * newIceConc) )**(1.0/3.0)
      
    ! --- Freezing parameterization type 2: all aerosol are assumed to have the same size --------------------------------------
    else if(frzMode%aeroStdDev < 1.1) then ! use this if aerosols can be assumed to have the same size: aeroStdDev is <= 1.1
      ! according to Kärcher and Lohmann 2002b
      ! this includes the size effect but assumes that all aerosol have the same size, so no sizebins are used
      
      minIceRad = frzMode%aeroRad
      
      growthTermFactors = 4.0 * PI * b1_Paper/b2_Paper**2.0 / SVOL
      TBBT = 2.0 * b1_Paper * b2_Paper * freezingTimeScale
      delta = b2_Paper * minIceRad            !BG: delta from KL 2002b, eq.7
      deltaPlus1 = 1.0 + delta
      kappa = TBBT / deltaPlus1**2.0       !BG: = kappa from KL 2002b, eq.8
      
      errorFuncApprox = 3.0 * sqrt(kappa) / ( 2.0 + sqrt(1.0+9.0*kappa/PI) ) ! according to Ren and Mackenzie 2005
      R_growthTerm_monoDisperse = growthTermFactors / deltaPlus1 * ( delta**2.0 - 1.0 &
        + (1.0 + 0.5*kappa*deltaPlus1**2.0) * errorFuncApprox/sqrt(kappa) )
        
      newIceConc = availableWaterForNucleation / R_growthTerm_monoDisperse
      newIceRad  = ( ( 1.0 + 0.5 * sqrt(kappa) * errorFuncApprox) * deltaPlus1 - 1.0 ) / b2_Paper
      newIceMass = PI * RHOICE / 0.75 * newIceConc * newIceRad**3.0
      
      ! if there are not enough aerosols, distribute the complete mass on the existing ones...is this a good assumption?
      newIceConc = min( newIceConc, frzMode%aeroConc )
      newIceRad  = ( newIceMass / ( PI * RHOICE / 0.75*newIceConc) )**(1.0/3.0)
      
    ! --- Freezing parameterization type 3: use size bins to include size effect -------------------------------
    else
      ! according to Kärcher and Lohmann 2002b
      ! now calculate freezing with considering an aerosol size distribution
      ! therefore size bins are created. the aerosols are distritributed into the sizebins
      ! then going from the biggest bin to the lowest, the aerosols are frozen
      ! until all available water vapor has been used, the aerosols below this point are not frozen
      
      nBins = 1 + INT( log( (radMax/radMin)**3.0 ) / log(VRAT) )
      if(lCheck) then
        if(nBins >= nBinsMax) write(debugOutput,*) "ERROR in mo_cirrus_4: nBinsMax too low: ", nBins
        return
      else
        if(nBins >= nBinsMax) nBins = nBinsMax - 1
      end if

      ! --- Calculate distribution of aerosols into the size bins ----------------------------------------------
      
      DLOGR0 = 2.0**THIRD * (VRAT**THIRD - 1.0) / (VRAT+1.0)**THIRD
      nAllAeroInBins = 0.0
      binsRad(nBins+1) = radMax * VRAT**THIRD

      do i = 1, nBins
        
        binsRad(i) = radMin * VRAT**( THIRD*FLOAT(i-1) )
        
        SIGL     = log( max( frzMode%aeroStdDev, 1.1 ) )
        SHAPFC1  = 1.0   / ( SQRT(TWOPI) * SIGL )
        SHAPFC2  = 0.5 / SIGL**2.0
        ARG      = SHAPFC2  * ( log(binsRad(i)/frzMode%aeroRad) )**2.0
        ARG      = min( ARG, 75.0 )
        
        nAeroInBin(i) = DLOGR0 * frzMode%aeroConc * SHAPFC1 * EXP(-ARG)
        nAllAeroInBins = nAllAeroInBins + nAeroInBin(i)
        
      enddo
      
      nAeroInBin = nAeroInBin / nAllAeroInBins * frzMode%aeroConc
      
      ! --- Go through size bins, biggest to smallest --------------------------------------------------
      
      sumFrozenAero = 1e-35     ! to avoid division by 0
      sumDeposWaterMole = 1e-25 ! not good style: understand algorithm and improve...
      radMean = 0.0
      
      growthTermFactors = 4.0 * PI * b1_Paper/b2_Paper**2.0 / SVOL
      TBBT = 2.0 * b1_Paper * b2_Paper * freezingTimeScale
      
      lastFrozenBin = -1
      do i = nBins, 1, -1
        
        ! calculate the water molecules that condense on the aerosol in this size bin
        delta = b2_Paper * binsRad(i)
        deltaPlus1 = 1.0 + delta
        kappa = TBBT / deltaPlus1**2.0
        
        errorFuncApprox = 3.0 * sqrt(kappa) / ( 2.0 + sqrt(1.0+9.0*kappa/PI) ) ! according to Ren and Mackenzie 2005
        R_growthTerm_monoDisperse = growthTermFactors / deltaPlus1 * ( delta**2.0 - 1.0 &
          + (1.0 + 0.5*kappa*deltaPlus1**2.0) * errorFuncApprox/sqrt(kappa) )
          
        sumDeposWaterMoleLast = sumDeposWaterMole
        sumFrozenAeroLast = sumFrozenAero
        
        sumDeposWaterMole = sumDeposWaterMole + R_growthTerm_monoDisperse * nAeroInBin(i)
        sumFrozenAero = sumFrozenAero + nAeroInBin(i)
        radMean = radMean + binsRad(i) * nAeroInBin(i)
        
        ! if all available water is consumed
        if(sumDeposWaterMole >= availableWaterForNucleation) then
          
          lastFrozenBin = i
          exit ! all available water frozen therefore exit loop over size bins
          
        end if
        
      enddo
      
      ! --- Calculate new iceConc and iceRad --------------------------------------------------
      
      ! loop is finished, either because all water is consumed (lastFrozenBin>0)
      ! or because all aerosol are frozen (lastFrozenBin==-1)
      if(lastFrozenBin>0) then
        
        ! lastFrozenBin is the last bin that freezes, however only partly
        ! therefore calculate the iceRad and iceConc that is somewhere between this bin and the last one
        RLOGRAT = log( binsRad(lastFrozenBin) / binsRad(lastFrozenBin+1) )
        SLOPER = log( sumDeposWaterMole / sumDeposWaterMoleLast ) / RLOGRAT
        minIceRad = binsRad(lastFrozenBin+1) * ( availableWaterForNucleation / sumDeposWaterMoleLast )**(1.0/SLOPER)
        
        SLOPE = log( sumFrozenAero / sumFrozenAeroLast ) / RLOGRAT
        newIceConc = sumFrozenAeroLast * ( minIceRad / binsRad(lastFrozenBin+1) )**SLOPE
        
        radMean = radMean / frzMode%aeroConc
        
        deltaPlus1 = 1.0 + b2_Paper * max( minIceRad, radMean )
        kappa = TBBT / deltaPlus1**2.0
        errorFuncApprox = 3.0 * sqrt(kappa) / ( 2.0 + sqrt(1.0+9.0*kappa/PI) )
        RIRAT = 1.0 + 0.5 * sqrt(Kappa) * errorFuncApprox
        newIceRad = ( RIRAT * deltaPlus1 - 1.0 ) / b2_Paper
        
        ! newIceMass not calculated --- but also not needed...
        
      else if(lastFrozenBin==-1) then ! smallest bin
        
        ! all aerosol freeze; calculate their iceRad and iceConc
        
        minIceRad = binsRad(1)
        radMean = radMean / frzMode%aeroConc
        
        deltaPlus1 = 1.0 + b2_Paper * radMean
        kappa = TBBT / deltaPlus1**2.0
        errorFuncApprox = 3.0 * sqrt(kappa) / ( 2.0 + sqrt(1.0+9.0*kappa/PI) )
        RIRAT = 1.0 + 0.5 * sqrt(Kappa) * errorFuncApprox
        newIceRad = ( RIRAT * deltaPlus1 - 1.0 ) / b2_Paper
        
        newIceConc = frzMode%aeroConc * availableWaterForNucleation / sumDeposWaterMole ! sure that this is correct???
        newIceMass = PI * RHOICE / 0.75 * newIceConc * newIceRad**3.0
        newIceConc = frzMode%aeroConc
        newIceRad = ( newIceMass / (PI * RHOICE / 0.75 * newIceConc) )**(1.0/3.0)
        
      end if
      
    end if ! --- End Freezing parameterization types --------------------------------------------------
        
    ! --- Set new values of the freezingMode --------------------------------------------------
    
    ! for heterogeneous freezing the newIceMass is multiplied by CVF = 0.99
    ! however, as we do not use newIceMass further, it seems not necessary
    !if(frzMode%freezingType>1) newIceMass = newIceMass * 0.99
    
    if(lDontUseFreezRad) frzMode%iceRad = frzMode%aeroRad
    
    ! adjust frzMode properties
    frzMode%aeroConc = frzMode%aeroConc - newIceConc
    frzMode%aeroRad = minIceRad ! this should be done better!!!
    !frzMode%aeroStdDev ! not changed at the moment
    if(.not.lDontUseFreezRad) frzMode%iceRad = newIceRad
    frzMode%iceConc = newIceConc
    frzMode%lFrozen = .true.

    frzMode%tempOfCurrentConc = temp
    
    ! change in saturation should be calculated here, is not done in old code

  end subroutine CalculateFreezingEvent
  
  !--- Helper function --------------------------------------------------------
    
  function AllFrozen(modes)
    
    implicit none
    
    ! in/out
    type(FreezingMode), intent(in) :: modes(nFreezingModes)
    logical :: AllFrozen
    
    integer :: i
    
    AllFrozen = .true.
    
    do i = 1, nFreezingModes
      if(.not.modes(i)%lFrozen .and. modes(i)%freezingType /= -1) then
        if(lDustDepRecalcActiveFraction .and. any(dustDepFreezingTypes == modes(i)%freezingType)) cycle
        AllFrozen = .false.
        return
      end if
    end do
    
  end function AllFrozen
  
  subroutine NextCriticalSaturation(modes, temp, nextSatCrit, nextFreezingIndex)
    
    implicit none
    
    ! in/out
    type(FreezingMode), intent(in) :: modes(nFreezingModes)
    real, intent(in) :: temp
    real, intent(out) :: nextSatCrit
    integer, intent(out) :: nextFreezingIndex
    
    integer :: i
    real :: satCrit
    
    nextSatCrit = -1.0
    nextFreezingIndex = -1
    
    do i = 1, nFreezingModes
      if(.not.modes(i)%lFrozen .and. modes(i)%freezingType /= -1) then
        
        if(lDustDepRecalcActiveFraction .and. any(dustDepFreezingTypes == modes(i)%freezingType)) cycle
        
        satCrit = CriticalSaturationForFreezing(modes(i), temp)
        
        if(nextSatCrit==-1.0 .or. satCrit < nextSatCrit) then
          nextSatCrit = satCrit
          nextFreezingIndex = i
        end if

      end if
    end do
    
  end subroutine NextCriticalSaturation
  
  subroutine GetDustDepIndex(modes, number, indices)
    
    implicit none
    
    ! in/out
    type(FreezingMode), intent(in) :: modes(nFreezingModes)
    integer, intent(out) :: number, indices(nFreezingModes)
    
    integer :: i
    
    number = 0
    
    do i = 1, nFreezingModes
      if(.not.modes(i)%lFrozen .and. any(dustDepFreezingTypes == modes(i)%freezingType)) then
        
        if(modes(i)%aeroConc > 0.0 .and. modes(i)%aeroRad > 0.0 &
          .and. modes(i)%aeroStdDev > 0.0) then
        
          number = number + 1
          indices(number) = i
          
        end if

      end if
    end do
    
  end subroutine GetDustDepIndex
  
  subroutine AdjustConcentrationsToNewConditions(modes, newTemp)
    ! Adjusts aerosol and ice concentrations of all modes to the
    ! adiabatic compression calculated from the temperature change.
    
    implicit none
    
    ! in/out
    type(FreezingMode), intent(inout) :: modes(nFreezingModes)
    real, intent(in) :: newTemp
    
    integer :: i
    real :: conversionFactor
    
    do i = 1, nFreezingModes
      
      conversionFactor = (newTemp/modes(i)%tempOfCurrentConc)**2.5
      
      modes(i)%aeroConc = modes(i)%aeroConc * conversionFactor
      modes(i)%iceConc = modes(i)%iceConc * conversionFactor
      
      modes(i)%tempOfCurrentConc = newTemp
      
    end do
    
  end subroutine AdjustConcentrationsToNewConditions
  
  function SatPresOverIce(temp)

    ! ***** VAPOR PRESSURE OVER ICE IN MBAR (temp IN K)

    implicit none

    real, intent(in) :: temp
    real :: SatPresOverIce

    real, parameter :: A = 0.01
    real, parameter :: B = 12.537
    real, parameter :: C = -2663.5

    ! Version 1: J.MARTI and K.MAUERSBERGER, GRL 20(5), 363-366, 1993
    SatPresOverIce = A * 10.0**( B + C/temp )
    
    ! Version 2: Rogers and Yau (1989) ! maybe in different unit! check that!
    !SatPresOverIce = 3.6e10 * exp(-(6145.0/temp))

  end function SatPresOverIce
  
  subroutine CalculateA1A2A3B1B2(temp, pres, sat, a1_Paper, a2_Paper, a3_Paper, b1_Paper, b2_Paper)

    ! ***** VAPOR PRESSURE OVER ICE IN MBAR (temp IN K)

    implicit none

    real, intent(in) :: temp, pres, sat
    real, intent(out) :: a1_Paper, a2_Paper, a3_Paper, b1_Paper, b2_Paper
    
    real :: satPres, waterToCrystalFlux, vaporConcAtSatPres

    ! calculate a1,a2,a3,b1,b2 according to Kärcher et al. (2006)
    satPres = 3.6e10 * exp(-(6145.0/temp))
    waterToCrystalFlux = depositionCoefficient / 4 * sqrt(11713803.0*temp) ! sqrt calculates thermal velocity of molecules
    vaporConcAtSatPres = 7.24637701e18 * satPres / temp ! ideal gas law; number is inverse boltzmann constant with cm in unit

    ! see Kärcher et al. 2006 for the formulas of a1,a2,a3,b1,b2
    a1_Paper = ( 0.601272523 / temp - 0.000342181855 ) / temp
    a2_Paper = 1.0 / vaporConcAtSatPres
    a3_Paper = 1.49236645e-12 / (temp * pres)
    b1_Paper = waterToCrystalFlux * 3.23e-23 * vaporConcAtSatPres * (sat-1.0) ! number is specVolumeWater
    b2_Paper = waterToCrystalFlux * 249.239822 * pres / temp**1.94 ! 249.239822 * pres / temp**1.94 is diffusion coefficient ! C = 1

  end subroutine CalculateA1A2A3B1B2

end module 