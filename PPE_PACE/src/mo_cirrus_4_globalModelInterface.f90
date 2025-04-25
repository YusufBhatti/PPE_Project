module mo_cirrus_4_globalModelInterface

  ! This file is the interface between the 2 moment cloud routine of ECHAM-HAM and the cirrus box model.
  ! It sets the concentrations for the different freezing modes, calls the cirrus model for the freezing calculation,
  ! and calculates the deposition of water vapor on existing and new ice crystals.
  !
  ! In this file the units used units are l-1 for concentrations and µm for radii (also especially in the freezing modes).
  ! Values that are from or given to the global model usually are in SI unit (m-3 and m, with exceptions).
  !
  ! Written in 2016 by Steffen Münch at ETH Zürich
  ! Before it was rewritten, it was used by Miriam Kuebbeler (MK) and Blaz Gasparini (BG). I kept their comments when useful.

!#define CIRRUS_MODULE_DEBUG

  use mo_kind, only: dp, wp
  use mo_linked_list, only: t_stream
  use mo_echam_cloud_params, only : cloudParams_homogeneousIceTemperature=>cthomi
  use mo_cirrus_4, only: FreezingMode, ActiveFractionDustDep
  use mo_physical_constants, only: als, cpd, rv

  use mo_ham_m7ctl, only: iacci, icoai, iaccs, icoas
  use mo_ham_species, only: id_du
  use mo_ham, only : sizeclass

#ifdef CIRRUS_MODULE_DEBUG
  use CustomOutputStream
#endif
  
  implicit none
  
  logical :: lDustDepRecalcActiveFraction = .true.
  logical :: lUseAscendForDeposition = .false. ! this needs the new deposition calculation activated in the cirrus scheme
  logical :: lRemoveFrozenDustFromTracer = .false.
  logical :: lPoreCondensationFreezing = .false.
  
  logical :: lSeeding = .false.
  real, parameter :: seedingConc = 10.0 ! l-1
  real, parameter :: seedingRad = 50.0 ! µm
  
  !>> --- used for lDebug = true ------------------------------------------------------------------------------
#ifdef CIRRUS_MODULE_DEBUG
  
  logical :: lDebug = .true.
  integer :: outputFileUnit = 103
  type(t_stream), pointer :: cirrusStream
  
  type(NormalSE) :: geopotHeightStream
  type(NormalSE) :: presStream
  type(NormalSE) :: tempStream
  type(NormalSE) :: updraftStream
  type(NormalSE) :: satStream
  type(NormalSE) :: satAfterStream
  
  type(TimedValueSE) :: satClearSkyStream, satCloudySkyStream
  
  type(ConcRadSE) :: homoAeroInStream, homoAeroOutStream, homoIceStream
  type(ConcRadSE) :: dustImmAeroInStream, dustImmAeroOutStream, dustImmIceStream
  type(ConcRadSE) :: dustDepAAeroInStream, dustDepAAeroOutStream, dustDepAIceStream
  type(ConcRadSE) :: dustDepCAeroInStream, dustDepCAeroOutStream, dustDepCIceStream
    
  type(ConcRadSE) :: homoAeroAllStream, dustImmAeroAllStream, dustDepAeroAllStream, dustDepIceAllStream
  type(NormalSE) :: dustImmVolFractionSE
  type(ConcRadSE) :: dustImmAccAeroInStream, dustImmCoaAeroInStream
    
  type(NormalSE) :: homoIceFreezingTime, dustImmIceFreezingTime, dustDepAIceFreezingTime, dustDepCIceFreezingTime
    
  type(ConcRadSE) :: stratiformIceStream, detrainedIceStream
  type(NormalSE) :: stratiformIceRadOutStream, detrainedIceRadOutStream
  type(ConcRadSE) :: stratiformIceLarger150
    
  type(ConcRadSE) :: newIceCrystalsStream

  type(NormalSE) :: depositedWaterVaporStream
  
  type(NormalSE) :: calculateCirrusFreqStream
  
  type(ConcRadSE) :: seedAeroInStream, seedAeroOutStream, seedIceStream
  type(NormalSE) :: seedIceFreezingTime
    
  real(wp), pointer :: radForcingLW(:,:), radForcingSW(:,:), radForcingNet(:,:)
  real(wp), pointer :: radForcingColdBaseLW(:,:), radForcingColdBaseSW(:,:), radForcingColdBaseNet(:,:)
  real(wp), pointer :: radForcingColdBase2LW(:,:), radForcingColdBase2SW(:,:), radForcingColdBase2Net(:,:)
  real(wp), pointer :: radForcingWarmBaseLW(:,:), radForcingWarmBaseSW(:,:), radForcingWarmBaseNet(:,:)
  real(wp), pointer :: radForcingWarmBase2LW(:,:), radForcingWarmBase2SW(:,:), radForcingWarmBase2Net(:,:)
  
  type(TimedValueSE) :: satAboveHomoSatCrit
  type(ConcRadSE) :: iceByOvershoot
  type(TimedValueSE) :: tempBelow238
    
  real(wp), pointer :: iceCloudTime(:,:,:), iceCloudIWC(:,:,:), iceCloudEffRad(:,:,:), iceCloudICNC(:,:,:), &
    iceCloudCover(:,:,:), iceCloudOpticalDepth(:,:,:)
  real(wp), pointer :: cirrusCod3Time(:,:,:), cirrusCod3IWC(:,:,:), cirrusCod3EffRad(:,:,:), cirrusCod3ICNC(:,:,:), &
    cirrusCod3Cover(:,:,:), cirrusCod3OpticalDepth(:,:,:)
  
#endif
  !<< --- END used for lDebug = true --------------------------------------------------------------------------
  
  public :: CirrusModelInterface, InitCirrusModel, ConstructCirrusModelStreams
#ifdef CIRRUS_MODULE_DEBUG
  public :: iceCloudTime, iceCloudIWC, iceCloudEffRad, iceCloudICNC, iceCloudCover, iceCloudOpticalDepth
  public :: cirrusCod3Time, cirrusCod3IWC, cirrusCod3EffRad, cirrusCod3ICNC, cirrusCod3Cover, cirrusCod3OpticalDepth
#endif
  
  contains
  
  subroutine CirrusModelInterface( bp,lev,bn, & ! multi processing block position
    simulationTimeIn, geopotHeight, pres, temp, updraftVelocityIn, & ! cirrus scheme in
    specificHumidity, satSpecificHumidityIce, aeroConcHomo, aeroRadHomo, & ! cirrus scheme in
    cloudIceIn, airDensity, stratiformIC, detrainedIC, detrainedRad, &
    pxtm1, pxtte, &
    cloudCover, airDensCorrectionForICfallVelocity, & ! for deposition calculation
    themodynTermForIceNucleation, diffusionCoeff, viscosity, & ! for deposition calculation
    newIceCrystals, depositedWaterVapor & ! output
    )
    
    use mo_ham_streams, only: extMixedAccuDust=>nduinsolai, &       !ext mixed accu. dust AP for frz in stream_ham
                                  extMixedCoarseDust=>nduinsolci, & !ext mixed coarse dust AP for frz in stream_ham 
                                  intMixedDust=>ndusol_strat ,&     !int mixed dust aerosols from strat. activ. part.  
                                  dryRadius=>rdry, wetRadius=>rwet
                                  ! BG: put these in cloud_subm_1=>submodelCloud1???
    use mo_ham_freezing, only: aero_massvolratio, aero_nc_surfw
    use mo_cloud_utils, only: conv_effr2mvr, kb, xmw, alpha, clc_min, fall, &
          sedimen_alpha1=>alfased_1, sedimen_alpha2=>alfased_2, sedimen_alpha3=>alfased_3, &
          sedimen_beta1=>betased_1, sedimen_beta2=>betased_2, sedimen_beta3=>betased_3, &
          icMass_border1=>ri_vol_mean_1, icMass_border2=>ri_vol_mean_2, &
          icMassAssumed=>mi, icemin, cloudParams_totalWaterMin=>cqtmin
    use mo_physical_constants, only : meltingTemp=>tmelt
    use mo_math_constants, only : PI=>pi
    
    USE mo_ham, only: naerocomp, aerocomp
    
    use mo_cirrus_4, only: ModelInterfaceUnitsLiterMicron

    implicit none
    
    integer, intent(in) :: bp,lev,bn ! multi processing position in output streams: bp = position in block, lev = grid point level, bn = block number

    real(dp), intent(in) :: simulationTimeIn ! s
    real(dp), intent(in) :: geopotHeight
    real(dp), intent(in) :: pres ! hPa
    real(dp), intent(in) :: temp ! K
    real(dp), intent(in) :: updraftVelocityIn ! cm/s
    real(dp), intent(in) :: specificHumidity ! kg/kg
    real(dp), intent(in) :: satSpecificHumidityIce ! kg/kg
    !real(dp), intent(in) :: satIn ! as sat is calculated from superSat this should be > 1.0 if cirrus scheme should be called
    
    real(dp), intent(in) :: aeroConcHomo ! m-3
    real(dp), intent(in) :: aeroRadHomo ! m
    real(dp), intent(in) :: cloudIceIn ! kg/kg
    real(dp) :: cloudIce ! kg/kg
    real(dp), intent(in) :: airDensity ! kg/m3
    real(dp), intent(in) :: stratiformIC, detrainedIC, detrainedRad ! m-3 and m
    
    REAL(dp), INTENT(in) :: pxtm1(:,:,:) !< tracer mmr (t-1)
    REAL(dp), INTENT(inout) :: pxtte(:,:,:) !< tracer tendency
    
    ! for deposition calculation
    real(dp), intent(in) :: cloudCover, airDensCorrectionForICfallVelocity, &
    themodynTermForIceNucleation, diffusionCoeff, viscosity
    
    real(dp), intent(out) :: newIceCrystals ! m-3
    real(dp), intent(out) :: depositedWaterVapor ! kg/kg

    integer, parameter :: nFreezingModes = 7
    
    type(FreezingMode) :: homoMode, dustImmMode, dustDepAMode, dustDepCMode, stratiformIceMode, detrainedIceMode
    type(FreezingMode) :: seedMode
    type(FreezingMode) :: freezingModes(nFreezingModes)
      
    real :: simulationTime ! s
    real :: updraftVelocity ! cm/s
    real :: sat
    
    real :: dustSolAcc, dustSolCoa
    real(dp) :: volRatio(1,1), numConc(1,1), tracerConc(1,1)
    real(dp) :: massAcc, massCoa
    
    real :: activeFraction, a, S0, param
    
    real :: newIceCrystalsRad

    ! for deposition calculations
    real, parameter :: RHOICE = 925.0 ! kg/m3
    real :: iceMassBefore, icMass ! kg
    real :: icFallVel_Param1, icFallVel_Param2, icMassFallVelocity
    real :: gtp, thermalVelocity, b2, fuchs, reynold, ventilation
    real :: depositedWaterVaporInKgPerL
    !real :: massOut
    
    !real :: detrainedFraction, icRad
    
    real :: allDepositedWater
    integer :: i
    
    real :: overshootIceConc, overshootIceRad
    
    logical :: lPotentialCirrusConditions

    integer :: jn, jclass, jspec, jt
    integer :: iTracerDustDepAccMass, iTracerDustDepCoaMass, iTracerDustDepAccNumber, iTracerDustDepCoaNumber
    
    
    simulationTime = simulationTimeIn
    updraftVelocity = updraftVelocityIn
    if(updraftVelocity > 1000) updraftVelocity = 1000
    
    ! limit to 1 km upward movement
    ! CHECK HOW MUCH THIS EFFECTS GLOBAL MODEL RESULTS
    !if(updraftVelocity * simulationTime > 100000) then
    !  simulationTime = 100000 / updraftVelocity
      !if(lDebug) write(outputFileUnit,*) "Simulation time limited due to high updraft:", updraftVelocity, simulationTime / simulationTimeIn
      !print*, "Reduce simulation time to not exceed 1 km to", simulationTime
    !end if

    ! convert gridmean ice mass to incloud ice mass
    if(cloudCover >= clc_min) then
      cloudIce = cloudIceIn / cloudCover
    else
      cloudIce = 0.0
    end if
    cloudIce = max(cloudIce, 0.0)
    
    sat = specificHumidity / satSpecificHumidityIce
    
    
    lPotentialCirrusConditions = (updraftVelocity > 0.01 & !releaseLimit .and. updraftVelocity > iceDowndraft & !not checking that at the moment
      .and. sat > 1.0 .and. temp < cloudParams_homogeneousIceTemperature)
      !.and. (stratiformIC+detrainedIC) < cloudParams_totalWaterMin
    
    !>> --- used for lDebug = true ------------------------------------------------------------------------------
#ifdef CIRRUS_MODULE_DEBUG
    if(lDebug) then
      ! write output of input values
      
      call SaveNormalValue(bp,lev,bn, geopotHeightStream, geopotHeight)
      call SaveNormalValue(bp,lev,bn, presStream, pres)
      call SaveNormalValue(bp,lev,bn, tempStream, temp)
      call SaveNormalValue(bp,lev,bn, updraftStream, updraftVelocity)
      call SaveNormalValue(bp,lev,bn, satStream, sat)

      if(temp < cloudParams_homogeneousIceTemperature) then
        call SaveTimedValue(bp,lev,bn, tempBelow238, real(temp))
        
        if(cloudIce < 1e-7) then
          ! stream clear sky sat
          call SaveTimedValue(bp,lev,bn, satClearSkyStream, sat)
        else
          ! stream cloudy sat
          call SaveTimedValue(bp,lev,bn, satCloudySkyStream, sat)
        end if
      end if
      
      if(sat > 2.418 - ( min(240.0, max(170.0, temp)) /245.68)) then
        call SaveTimedValue(bp,lev,bn, satAboveHomoSatCrit, sat)
      end if
      
      if(lPotentialCirrusConditions) call SaveNormalValue(bp,lev,bn, calculateCirrusFreqStream, 1.0)
      
    end if
#endif
    !<< --- END used for lDebug = true --------------------------------------------------------------------------
    
#ifndef CIRRUS_MODULE_DEBUG
    ! if not in debug mode, the homo and dust modes only have to be set if we have potential cirrus conditions (saves computational time)
    if(lPotentialCirrusConditions) then
#endif
    
      ! homogeneous -----------------------------------------------------------------------------------------------
      homoMode%freezingType = 1
      homoMode%aeroConc = (aeroConcHomo) * 1e-3 ! m-3 -> l-1
      homoMode%aeroRad = aeroRadHomo * 1e6 ! m -> µm
      homoMode%aeroStdDev = 1.0

      !releaseLimit homoMode%aeroConc = max(homoMode%aeroConc, 1.e-2) !BG: due to the numerical stability of cirrus scheme
      !releaseLimit homoMode%aeroRad = max(homoMode%aeroRad, 1.01e-3) !due to num stability
  
  
      ! heterogeneous - dust immersion mode ------------------------------------------------------------------------
      dustImmMode%freezingType = 2
      
      !dustImmMode%aeroConc = (intMixedDust(bp,lev,bn)) * 1e-3 ! m-3 -> l-1
      ! For dust immersion freezing: either take all soluble concentration acc+coa
      !dustSolAcc = pxtm1(bp,lev,sizeclass(iaccs)%idt_no) * airDensity ! soluble dust accumulation mode
      !dustSolCoa = pxtm1(bp,lev,sizeclass(icoas)%idt_no) * airDensity ! soluble dust coarse mode
      
      ! tracer index for dustDep accumulation number: sizeclass(iaccs)%idt_no
      ! tracer index for dustDep coarse number: sizeclass(icoas)%idt_no
      ! tracer index for dustDep accumulation mass: aerocomp(jn)%idt
      ! tracer index for dustDep coarse mass: 
      
      ! Or weight the total concentration by surface
      tracerConc(1,1) = pxtm1(bp,lev,sizeclass(iaccs)%idt_no)*airDensity
      CALL aero_massvolratio(1, 1, 1, iaccs, id_du, 'volu', pxtm1(bp,lev,:), volRatio)
      CALL aero_nc_surfw(1, 1, 1, volRatio, tracerConc, numConc)
      dustSolAcc = numConc(1,1)
      tracerConc(1,1) = pxtm1(bp,lev,sizeclass(icoas)%idt_no)*airDensity
      CALL aero_massvolratio(1, 1, 1, icoas, id_du, 'volu', pxtm1(bp,lev,:), volRatio)
      CALL aero_nc_surfw(1, 1, 1, volRatio, tracerConc, numConc)
      dustSolCoa = numConc(1,1)
      
      ! Sum of soluble dust accu+coar
      dustImmMode%aeroConc = (dustSolAcc + dustSolCoa) * 1e-3 ! m-3 -> l-1
      
      dustImmMode%aeroConc = dustImmMode%aeroConc * 0.05 !MK: only 5% of DU acts as IN
      
      !dustImmMode%aeroRad = aeroRadHomo * 1e6 ! m -> µm
      if(dustImmMode%aeroConc > 0.0) &
        dustImmMode%aeroRad = (wetRadius(iaccs)%ptr(bp,lev,bn)*dustSolAcc + &
                wetRadius(icoas)%ptr(bp,lev,bn)*dustSolCoa) / (dustSolAcc+dustSolCoa) * 1e6 ! m -> µm
                ! careful! Should here radius or mass be weighted?
      
      dustImmMode%aeroStdDev = 1.0
    
      !releaseLimit dustImmMode%aeroConc =  max(dustImmMode%aeroConc, 1.e-2) !MK: to avoid concentrations less than 0.01/L
      !releaseLimit dustImmMode%aeroRad = max(dustImmMode%aeroRad, 1.01e-3) !MK: to avoid particles smaller than 1nm
  
  
      ! heterogeneous - dust deposition nucleation mode --------------------------------------------------------------------------------------------------------
      
      !dustDepMode%freezingType = 3
      !dustDepMode%aeroConc = (extMixedAccuDust(bp,lev,bn) + extMixedCoarseDust(bp,lev,bn) & !dust accum+coarse insoluble
      !  ) * 1e-3 ! m-3 -> l-1
      !if(dustDepMode%aeroConc > 0.0) &
      !  dustDepMode%aeroRad = (dryRadius(iacci)%ptr(bp,lev,bn)*extMixedAccuDust(bp,lev,bn) + &
      !    dryRadius(icoai)%ptr(bp,lev,bn)*extMixedCoarseDust(bp,lev,bn)) &
      !    / (extMixedAccuDust(bp,lev,bn) + extMixedCoarseDust(bp,lev,bn)) * 1e6
          ! careful! Should here radius or mass be weighted?
          !BG: ACCUM DRY INSOLUBLE radius, added 18.12.2014
      !dustDepMode%aeroStdDev = 1.0
      
      ! dustDep accumulation mode
      if(lPoreCondensationFreezing) then
        dustDepAMode%freezingType = 5
      else
        dustDepAMode%freezingType = 3
      end if
      dustDepAMode%aeroConc = extMixedAccuDust(bp,lev,bn) * 1e-3 ! m-3 -> l-1
      dustDepAMode%aeroRad = dryRadius(iacci)%ptr(bp,lev,bn) * 1e6
      dustDepAMode%aeroStdDev = 1.0
      
      ! dustDep coarse mode
      if(lPoreCondensationFreezing) then
        dustDepCMode%freezingType = 6
      else
        dustDepCMode%freezingType = 3
      end if
      dustDepCMode%aeroConc = extMixedCoarseDust(bp,lev,bn) * 1e-3 ! m-3 -> l-1
      dustDepCMode%aeroRad = dryRadius(icoai)%ptr(bp,lev,bn) * 1e6
      dustDepCMode%aeroStdDev = 1.0
      
      if(lDustDepRecalcActiveFraction) then
        
        if(.not.lRemoveFrozenDustFromTracer) then
          dustDepAMode%subtractIce = stratiformIC*1e-3
          dustDepCMode%subtractIce = stratiformIC*1e-3
        end if
        
      else
        
        activeFraction = ActiveFractionDustDep(real(temp), sat, dustDepAMode%freezingType, dustDepAMode%aeroRad)
        if(activeFraction < 0.001) activeFraction = 0.001
        dustDepAMode%aeroConc = dustDepAMode%aeroConc * activeFraction !MK: only active fraction of min. dust acts as IN
        
        activeFraction = ActiveFractionDustDep(real(temp), sat, dustDepCMode%freezingType, dustDepCMode%aeroRad)
        if(activeFraction < 0.001) activeFraction = 0.001
        dustDepCMode%aeroConc = dustDepCMode%aeroConc * activeFraction !MK: only active fraction of min. dust acts as IN
        
      end if
  
      !releaseLimit dustDepMode%aeroConc = max(dustDepMode%aeroConc, 1.e-2) !MK: to avoid concentrations less than 0.01/L
      !releaseLimit dustDepMode%aeroRad = max(dustDepMode%aeroRad, 1.01e-3) !MK: to avoid particles smaller than 1nm
      
      ! seeding particle --------------------------------------------------------------------------------------------
      if(lSeeding) then
        seedMode%freezingType = 4
        seedMode%aeroConc = seedingConc
        seedMode%aeroRad = seedingRad
        seedMode%aeroStdDev = 1.0
      end if
  
#ifndef CIRRUS_MODULE_DEBUG
    end if
#endif
  
    ! existing ice --------------------------------------------------------------------------------------------------------
    if(cloudIce > 0.0) then
      stratiformIceMode%lFrozen = .true. ! preexisting ice
      stratiformIceMode%iceConc = (stratiformIC) * 1e-3 ! m-3 -> l-1
      !BG: simple parametrisation based on temperature for IC size
      !iceRad = max(23.2 * exp(0.015 * min(0.0_dp, temp - meltingTemp)), 1.0_dp)
      !detrainedIceRad = max(1.0e-6_dp,conv_effr2mvr*1.0e-6*iceRad)
      ! This already has been calculated in the cloud routine: detrainedRad
      ! Old version:
      !existingIceMode%iceRad = detrainedRad * 1e6 ! m -> µm
      ! New version (but wrong weighting):
      !existingIceMode%iceRad = (detrainedRad*detrainedConc + existingIceRad*(existingIceCrystals-detrainedConc)) &
      !   / existingIceCrystals * 1e6 ! m -> µm
      ! New version ( calculate mean radius from total volume ):
!>>DN bugfix
!      if(stratiformIC > 0.0) then
      if(stratiformIC > EPSILON(1._dp)) then
!>>DN bugfix
        stratiformIceMode%iceRad = ((cloudIce*airDensity/RHOICE) / (stratiformIC * 4./3.*PI))**(1./3.) * 1e6 ! m -> µm
        
        if(stratiformIceMode%iceRad > 150) then
#ifdef CIRRUS_MODULE_DEBUG
          call SaveConcRad(bp,lev,bn, stratiformIceLarger150, stratiformIceMode%iceConc, stratiformIceMode%iceRad)
#endif
          stratiformIceMode%iceRad = 150
        end if
      end if
    
      if(stratiformIceMode%iceConc*1e3<cloudParams_totalWaterMin) stratiformIceMode%iceConc = 0.0 ! cloudParams_totalWaterMin = 1e-12
      !releaseLimit existingIceMode%iceConc = min(existingIceMode%iceConc, 1.0e4) !BG: max preex ice conc = 1e+7/m3 -> 10000 /l
      !releaseLimit existingIceMode%iceRad = min(existingIceMode%iceRad, 1.0e3)  !BG: max radius of preex ice = 1 mm
    end if
    if(detrainedIC > cloudParams_totalWaterMin) then ! cloudParams_totalWaterMin = 1e-12
      detrainedIceMode%lFrozen = .true. ! preexisting ice
      detrainedIceMode%iceConc = (detrainedIC) * 1e-3 ! m-3 -> l-1
      detrainedIceMode%iceRad = detrainedRad * 1e6 ! m -> µm
    end if
    
    
    !>> --- used for lDebug = true ------------------------------------------------------------------------------
#ifdef CIRRUS_MODULE_DEBUG
    if(lDebug) then
      
      call SaveConcRad(bp,lev,bn, homoAeroInStream, homoMode%aeroConc, homoMode%aeroRad)
      call SaveConcRad(bp,lev,bn, dustImmAeroInStream, dustImmMode%aeroConc, dustImmMode%aeroRad)
      call SaveConcRad(bp,lev,bn, dustDepAAeroInStream, dustDepAMode%aeroConc, dustDepAMode%aeroRad)
      call SaveConcRad(bp,lev,bn, dustDepCAeroInStream, dustDepCMode%aeroConc, dustDepCMode%aeroRad)
      
      call SaveNormalValue(bp,lev,bn, dustImmVolFractionSE, volRatio(1,1))

      call SaveConcRad(bp,lev,bn, stratiformIceStream, stratiformIceMode%iceConc, stratiformIceMode%iceRad)
      call SaveConcRad(bp,lev,bn, detrainedIceStream, detrainedIceMode%iceConc, detrainedIceMode%iceRad)
      
      call SaveConcRad(bp,lev,bn, homoAeroAllStream, real(aeroConcHomo*1e-3), real(aeroRadHomo*1e6))
      call SaveConcRad(bp,lev,bn, dustImmAeroAllStream, real((dustSolAcc + dustSolCoa) * 1e-3), dustImmMode%aeroRad)
      call SaveConcRad(bp,lev,bn, dustDepAeroAllStream, real(extMixedAccuDust(bp,lev,bn) + extMixedCoarseDust(bp,lev,bn))*1e-3, &
        real(dryRadius(iacci)%ptr(bp,lev,bn)*1e6))
        
      call SaveConcRad(bp,lev,bn, dustImmAccAeroInStream, real(dustSolAcc*1e-3), real(wetRadius(iaccs)%ptr(bp,lev,bn)*1e6))
      call SaveConcRad(bp,lev,bn, dustImmCoaAeroInStream, real(dustSolCoa*1e-3), real(wetRadius(icoas)%ptr(bp,lev,bn)*1e6))
      
      if(lSeeding) call SaveConcRad(bp,lev,bn, seedAeroInStream, seedMode%aeroConc, seedMode%aeroRad)
        
    end if
#endif
    !<< --- END used for lDebug = true --------------------------------------------------------------------------

    ! Assume that in all existing IC there is also one INP in there
    homoMode%aeroConc = homoMode%aeroConc - stratiformIC*1e-3
    dustImmMode%aeroConc = dustImmMode%aeroConc - stratiformIC*1e-3
    if(.not.lDustDepRecalcActiveFraction .and. .not.lRemoveFrozenDustFromTracer) then
      dustDepAMode%aeroConc = dustDepAMode%aeroConc - stratiformIC*1e-3
      dustDepCMode%aeroConc = dustDepCMode%aeroConc - stratiformIC*1e-3
    end if
    !if(lSeeding) seedMode%aeroConc = seedMode%aeroConc - stratiformIC*1e-3
    
    ! Security check
    if(homoMode%aeroConc < 0.0) homoMode%aeroConc = 0.0
    if(dustImmMode%aeroConc < 0.0) dustImmMode%aeroConc = 0.0
    if(dustDepAMode%aeroConc < 0.0) dustDepAMode%aeroConc = 0.0
    if(dustDepCMode%aeroConc < 0.0) dustDepCMode%aeroConc = 0.0
    if(seedMode%aeroConc < 0.0) seedMode%aeroConc = 0.0
    if(stratiformIceMode%iceConc < 0.0) stratiformIceMode%iceConc = 0.0
    if(detrainedIceMode%iceConc < 0.0) detrainedIceMode%iceConc = 0.0
    
    ! Combine all modes in input variable - order does not matter for cirrus scheme --------------------------------------
    freezingModes(1) = homoMode
    freezingModes(2) = dustImmMode
    freezingModes(3) = dustDepAMode
    freezingModes(4) = dustDepCMode
    freezingModes(5) = stratiformIceMode
    freezingModes(6) = detrainedIceMode
    if(lSeeding) freezingModes(7) = seedMode

    !iceMassBefore = 4.0/3.0 * PI * (existingIceMode%iceRad*1e-6)**3 * existingIceMode%iceConc * RHOICE ! in kg/l
    
    newIceCrystals = 0.0
    newIceCrystalsRad = 0.0
    depositedWaterVapor = 0.0
    
    
    ! --- 1. calculate newIceCrystals ----------------------------------------------------------------------------
    
    if (lPotentialCirrusConditions) then ! test for sat below 1.0
    
      ! Call cirrus box model
      call ModelInterfaceUnitsLiterMicron(simulationTime, real(pres), real(temp), real(updraftVelocity), sat, freezingModes)
      
      homoMode = freezingModes(1)
      dustImmMode = freezingModes(2)
      dustDepAMode = freezingModes(3)
      dustDepCMode = freezingModes(4)
      stratiformIceMode = freezingModes(5)
      detrainedIceMode = freezingModes(6)
      if(lSeeding) seedMode = freezingModes(7)
      
      newIceCrystals = (homoMode%iceConc + dustImmMode%iceConc + dustDepAMode%iceConc + dustDepCMode%iceConc) * 1e3 ! l-1 -> m-3
      if(lSeeding) newIceCrystals = newIceCrystals + seedMode%iceConc * 1e3 ! l-1 -> m-3

#ifdef CIRRUS_MODULE_DEBUG
      ! right now this is not used by anything except debug output
      if(newIceCrystals > 0.0) then
        if(.not.lSeeding) then
          newIceCrystalsRad = (homoMode%iceRad*homoMode%iceConc + dustImmMode%iceRad*dustImmMode%iceConc &
          + dustDepAMode%iceRad*dustDepAMode%iceConc + dustDepCMode%iceRad*dustDepCMode%iceConc) / (newIceCrystals*1e-3)
        else
          newIceCrystalsRad = (homoMode%iceRad*homoMode%iceConc + dustImmMode%iceRad*dustImmMode%iceConc &
          + dustDepAMode%iceRad*dustDepAMode%iceConc + dustDepCMode%iceRad*dustDepCMode%iceConc + seedMode%iceRad*seedMode%iceConc) &
          / (newIceCrystals*1e-3)
        end if
      end if
#endif
      
      if(lUseAscendForDeposition) then
        ! For old deposition calculation (doesn't work right now):
        !depositedWaterVaporInKgPerL = &
        !  4.0/3.0 * PI * ((homoMode%iceRad*1e-6)**3 - (homoMode%aeroRad*1e-6)**3) * homoMode%iceConc * RHOICE &
        !+ 4.0/3.0 * PI * ((dustImmMode%iceRad*1e-6)**3 - (dustImmMode%aeroRad*1e-6)**3) * dustImmMode%iceConc * RHOICE &
        !+ 4.0/3.0 * PI * ((dustDepMode%iceRad*1e-6)**3 - (dustDepMode%aeroRad*1e-6)**3) * dustDepMode%iceConc * RHOICE &
        !+ 4.0/3.0 * PI * (existingIceMode%iceRad*1e-6)**3 * existingIceMode%iceConc * RHOICE - iceMassBefore
        !if(lSeeding) depositedWaterVaporInKgPerL = depositedWaterVaporInKgPerL &
        !+ 4.0/3.0 * PI * ((seedMode%iceRad*1e-6)**3 - (seedMode%aeroRad*1e-6)**3) * seedMode%iceConc * RHOICE
        !depositedWaterVapor = depositedWaterVaporInKgPerL / (airDensity*1e-3) ! kg/l -> kg/kg
        !if(depositedWaterVapor > 0.0) depositedWaterVapor = min(depositedWaterVapor, specificHumidity - satSpecificHumidityIce )
        !if(depositedWaterVapor < 0.0) depositedWaterVapor = max(depositedWaterVapor,-cloudIce)

        ! For new deposition calculation:
        depositedWaterVapor = 0.0
        do i = 1, nFreezingModes
          depositedWaterVapor = depositedWaterVapor + freezingModes(i)%depositedWaterVapor !kg/m3
        end do
        depositedWaterVapor = depositedWaterVapor / airDensity !kg/kg
        ! think about that min statement below...
        if(depositedWaterVapor > 0.0) depositedWaterVapor = min(depositedWaterVapor, specificHumidity - satSpecificHumidityIce )
        if(depositedWaterVapor < 0.0) depositedWaterVapor = max(depositedWaterVapor,-cloudIce)
      
      else

        depositedWaterVapor = homoMode%depositedWaterVapor + dustImmMode%depositedWaterVapor &
          + dustDepAMode%depositedWaterVapor + dustDepCMode%depositedWaterVapor !kg/m3

        if(lSeeding) depositedWaterVapor = depositedWaterVapor + seedMode%depositedWaterVapor

        depositedWaterVapor = depositedWaterVapor / airDensity !kg/kg
        ! think about that min statement below...
        if(depositedWaterVapor > 0.0) depositedWaterVapor = min(depositedWaterVapor, specificHumidity - satSpecificHumidityIce )
        if(depositedWaterVapor < 0.0) depositedWaterVapor = max(depositedWaterVapor,-cloudIce)

      end if
    
    end if
    
    ! --- 2. calculate depositedWaterVapor ----------------------------------------------------------------------------
    
    if(.false.) then

    if (.not.lUseAscendForDeposition .or. .not.lPotentialCirrusConditions) then
    
      ! The whole deposition calculation is in SI units (m-3 and m). FreezingMode values have to be converted!
    
      sat = specificHumidity / satSpecificHumidityIce ! go back to sat before cirrus scheme
      allDepositedWater = 0.0
    
      do i = 1, nFreezingModes
    
        !if(freezingModes(i)%iceRad < 1.0) freezingModes(i)%iceRad = 1.0 ! CHANGE THIS LATER
    
        if(freezingModes(i)%iceConc > 0.0 .and. freezingModes(i)%iceRad > 0.0) then
      
          !if(cloudIce == 0.0 .or. cloudCover <= clc_min) existingIceCrystals = icemin
 
          ! air Density kg/m3
          ! cloudIce kg/kg
          ! existingIceCrystals 1/m3
          !icMass = max(airDensity*cloudIce / (freezingModes(i)%iceConc*1e3 * max(cloudCover,clc_min)), icMassAssumed)
          icMass = max( 4.0/3.0 * PI * (freezingModes(i)%iceRad*1e-6)**3 * RHOICE, icMassAssumed)
          if(freezingModes(i)%iceConc*1e3 < 1e-12) icMass = icMassAssumed

          if(icMass < icMass_border1) then
            icFallVel_Param1 = sedimen_alpha1
            icFallVel_Param2 = sedimen_beta1
          else if(icMass < icMass_border2) then
            icFallVel_Param1 = sedimen_alpha3
            icFallVel_Param2 = sedimen_beta3
          else ! > icMass_border2
            icFallVel_Param1 = sedimen_alpha2
            icFallVel_Param2 = sedimen_beta2
          end if
          icMassFallVelocity = fall*icFallVel_Param1 * icMass**icFallVel_Param2 * airDensCorrectionForICfallVelocity
 
          gtp = 1.0/(airDensity*themodynTermForIceNucleation)
          thermalVelocity = sqrt( 8.0*kb*temp / (PI*xmw) )
          b2 = 0.25 * alpha * thermalVelocity / diffusionCoeff
          fuchs = 1.0/(1.0+b2*freezingModes(i)%iceRad*1e-6)
          reynold = 2.0 * airDensity * freezingModes(i)%iceRad*1e-6 * icMassFallVelocity / viscosity
          ventilation = 1.0 + 0.229*sqrt(reynold)
          if(freezingModes(i)%iceConc*1e3 < 1e-12) ventilation = 1.0

          depositedWaterVapor = 4.0*PI*freezingModes(i)%iceRad*1e-6*(sat-1.0)*freezingModes(i)%iceConc*1e3 &
            *ventilation*gtp*fuchs*alpha*simulationTimeIn
          !if(freezingModes(i)%freezingTime > 0.0) then
            ! Calculate deposition only after crystal formation
          !  depositedWaterVapor = depositedWaterVapor * (simulationTime - freezingModes(i)%freezingTime) / simulationTime
          !end if
          if(depositedWaterVapor > 0.0) depositedWaterVapor = min(depositedWaterVapor, specificHumidity - satSpecificHumidityIce )
          if(depositedWaterVapor < 0.0) depositedWaterVapor = max(depositedWaterVapor,-cloudIce)
      
          allDepositedWater = allDepositedWater + depositedWaterVapor
        
          !sat = (specificHumidity - depositedWaterVapor) / satSpecificHumidityIce ! for sat after stream
      
          ! repair this later, also units...deactivated for now...but it is anyway only for debug output...
          !if(i==4 .and. .false.) then
            ! for existingIce%iceRad out stream
            !iceMassBefore = 4.0/3.0 * PI * existingIceMode%iceRad**3 * existingIceMode%iceConc * RHOICE*1e-3 ! in g as RHOICE = 0.925 ! g/cm3
          !  massOut = iceMassBefore + depositedWaterVapor*1e3*airDensity*1e-6 ! kg/kg -> g/cm3
          !  if(massOut>0.0) then
          !    existingIceMode%iceRad = (massOut / existingIceMode%iceConc / (RHOICE*1e-3) * 3.0/4.0 / PI)**(1.0/3.0)
          !  else
              !if(lDebug) write(outputFileUnit,*) "more mass evaporating than available", massOut, depositedWaterVapor, cloudIce
          !    existingIceMode%iceRad = 0.0
          !    existingIceMode%iceConc = 0.0
          !  end if
          !endif
 
        end if
      end do
      depositedWaterVapor = allDepositedWater

      ! Check again
      if(depositedWaterVapor > 0.0) depositedWaterVapor = min(depositedWaterVapor, specificHumidity - satSpecificHumidityIce )
      if(depositedWaterVapor < 0.0) depositedWaterVapor = max(depositedWaterVapor,-cloudIce)
      
    end if ! (.not.lUseAscendForDeposition .or. .not.lPotentialCirrusConditions)
    
    if(depositedWaterVapor > 0.0) depositedWaterVapor = min(depositedWaterVapor, &
    (specificHumidity - satSpecificHumidityIce) / ( 1+als**2*satSpecificHumidityIce / (cpd*rv*temp**2) ) )

    if(depositedWaterVapor /= depositedWaterVapor) &
    print*,"cirrus nan",lPotentialCirrusConditions,temp,updraftVelocity,cloudCover,cloudIce,stratiformIC
    
    ! The following statement was in Miriam's and Blaz's code. I am not sure why.
    ! The saturation inhibites negative deposition (=evaporation) which is allowed in nic_cirrus=2
    ! The temperature - 1 makes no sense for me, so I comment this statement out for now...
    !if(sat <= 1.0 .or. temp >= (cloudParams_homogeneousIceTemperature-1.0)) then
    !  depositedWaterVapor = 0.0
    !  newIceCrystals = 0.0
    !end if
          
    end if

    if(lRemoveFrozenDustFromTracer) then
      DO jn = 1,naerocomp    !loop over all mode-species

         jclass = aerocomp(jn)%iclass
         jspec  = aerocomp(jn)%spid
         jt     = aerocomp(jn)%idt

         if( jclass == iacci .and. jspec == id_du ) iTracerDustDepAccMass = jt
         if( jclass == icoai .and. jspec == id_du ) iTracerDustDepCoaMass = jt

      ENDDO !end loop over all mode-species
    
      iTracerDustDepAccNumber = sizeclass(iacci)%idt_no
      iTracerDustDepCoaNumber = sizeclass(icoai)%idt_no
    end if
    if(lRemoveFrozenDustFromTracer) then
      if(dustDepAMode%iceConc > 0.0) then
        massAcc = pxtm1(bp,lev,iTracerDustDepAccMass) / pxtm1(bp,lev,iTracerDustDepAccNumber)
      
        pxtte(bp,lev,iTracerDustDepAccNumber) = pxtte(bp,lev,iTracerDustDepAccNumber) &
        - dustDepAMode%iceConc*1e3/airDensity / simulationTime
        pxtte(bp,lev,iTracerDustDepAccMass) = pxtte(bp,lev,iTracerDustDepAccMass) &
        - dustDepAMode%iceConc*1e3/airDensity * massAcc / simulationTime
      end if
      if(dustDepCMode%iceConc > 0.0) then
        massCoa = pxtm1(bp,lev,iTracerDustDepCoaMass) / pxtm1(bp,lev,iTracerDustDepCoaNumber)
      
        pxtte(bp,lev,iTracerDustDepCoaNumber) = pxtte(bp,lev,iTracerDustDepCoaNumber) &
        - dustDepCMode%iceConc*1e3/airDensity / simulationTime
        pxtte(bp,lev,iTracerDustDepCoaMass) = pxtte(bp,lev,iTracerDustDepCoaMass) &
        - dustDepCMode%iceConc*1e3/airDensity * massCoa / simulationTime
      end if
    end if
        
    !>> --- used for lDebug = true ------------------------------------------------------------------------------
#ifdef CIRRUS_MODULE_DEBUG
    if(lDebug) then
      
	    if(lPotentialCirrusConditions) call SaveNormalValue(bp,lev,bn, satAfterStream, sat)
      
      call SaveConcRad(bp,lev,bn, homoAeroOutStream, homoMode%aeroConc, homoMode%aeroRad)
      call SaveConcRad(bp,lev,bn, dustImmAeroOutStream, dustImmMode%aeroConc, dustImmMode%aeroRad)
      call SaveConcRad(bp,lev,bn, dustDepAAeroOutStream, dustDepAMode%aeroConc, dustDepAMode%aeroRad)
      call SaveConcRad(bp,lev,bn, dustDepCAeroOutStream, dustDepCMode%aeroConc, dustDepCMode%aeroRad)
      
      call SaveNormalValue(bp,lev,bn, stratiformIceRadOutStream, stratiformIceMode%iceRad)
      call SaveNormalValue(bp,lev,bn, detrainedIceRadOutStream, detrainedIceMode%iceRad)
      
      if(homoMode%iceConc > 0.0) then
        call SaveConcRad(bp,lev,bn, homoIceStream, homoMode%iceConc, homoMode%iceRad)
        call SaveNormalValue(bp,lev,bn, homoIceFreezingTime, homoMode%freezingTime)
      end if
      if(dustImmMode%iceConc > 0.0) then
        call SaveConcRad(bp,lev,bn, dustImmIceStream, dustImmMode%iceConc, dustImmMode%iceRad)
        call SaveNormalValue(bp,lev,bn, dustImmIceFreezingTime, dustImmMode%freezingTime)
      end if
      if(dustDepAMode%iceConc > 0.0) then
        call SaveConcRad(bp,lev,bn, dustDepAIceStream, dustDepAMode%iceConc, dustDepAMode%iceRad)
        call SaveNormalValue(bp,lev,bn, dustDepAIceFreezingTime, dustDepAMode%freezingTime)
      end if
      if(dustDepCMode%iceConc > 0.0) then
        call SaveConcRad(bp,lev,bn, dustDepCIceStream, dustDepCMode%iceConc, dustDepCMode%iceRad)
        call SaveNormalValue(bp,lev,bn, dustDepCIceFreezingTime, dustDepCMode%freezingTime)
      end if
      if(dustDepAMode%iceConc > 0.0 .or. dustDepCMode%iceConc > 0.0) then
        call SaveConcRad(bp,lev,bn, dustDepIceAllStream, dustDepAMode%iceConc+dustDepCMode%iceConc, &
        ((dustDepAMode%iceRad)**3+(dustDepCMode%iceRad)**3)**(1./3.) )
      end if
      if(newIceCrystals > 0.0) then
        call SaveConcRad(bp,lev,bn, newIceCrystalsStream, real(newIceCrystals*1e-3), newIceCrystalsRad)
      end if

      call SaveNormalValue(bp,lev,bn, depositedWaterVaporStream, depositedWaterVapor*1e6) ! kg/kg -> g/kg
      
      if(lSeeding) then
        call SaveConcRad(bp,lev,bn, seedAeroOutStream, seedMode%aeroConc, seedMode%aeroRad)
        if(seedMode%iceConc > 0.0) then
          call SaveConcRad(bp,lev,bn, seedIceStream, seedMode%iceConc, seedMode%iceRad)
          call SaveNormalValue(bp,lev,bn, seedIceFreezingTime, seedMode%freezingTime)
        end if
      end if
      
      overshootIceConc = 0.0
      overshootIceRad = 0.0
      do i = 1, nFreezingModes
        if(freezingModes(i)%lFrozen .and. freezingModes(i)%freezingTime >= 0.0 .and. freezingModes(i)%iceConc > 0.0) then
          
          ! If freezing occured after a rise of 1 km
          if(freezingModes(i)%freezingTime > 100000 / updraftVelocity) then
            
            overshootIceConc = overshootIceConc + freezingModes(i)%iceConc
            overshootIceRad = overshootIceRad + freezingModes(i)%iceConc*freezingModes(i)%iceRad
            
          end if
                    
        end if
      end do
      if(overshootIceConc > 0.0) overshootIceRad = overshootIceRad / overshootIceConc
      call SaveConcRad(bp,lev,bn, iceByOvershoot, overshootIceConc, overshootIceRad)
      
    end if
#endif
    !<< --- END used for lDebug = true --------------------------------------------------------------------------
  
  end subroutine CirrusModelInterface

  !>> --- used for lDebug = true ------------------------------------------------------------------------------
  
  subroutine InitCirrusModel()
    
    ! call this from 
    ! init_subm in mo_submodel_interface or 
    ! activ_initialize in mo_activ or
    ! init_cloud_micro_2m in mo_cloud_utils
    
    use mo_exception, only: finish
    use mo_cirrus_4, only: debugOutput, lNewDeposition, depRecalcInCirrus=>lDustDepRecalcActiveFraction
    
    implicit none
    
    logical :: lOpened
    
    if(lUseAscendForDeposition .and. .not.lNewDeposition) &
      call finish("","ERROR: Cirrus module. Use ascend for deposition is only working if lNewDeposition is true!")
      
    if((lDustDepRecalcActiveFraction .and. .not.depRecalcInCirrus) .or. &
    (.not.lDustDepRecalcActiveFraction .and. depRecalcInCirrus)) &
      call finish("","ERROR: Cirrus module. lDustDepRecalcActiveFraction in global model interface and cirrus model don't match!")
    
#ifdef CIRRUS_MODULE_DEBUG
    if(lDebug) then
    
      ! open cirrus output file for messages
    
      INQUIRE (UNIT=outputFileUnit,OPENED=lOpened)
    
      if(.not.lOpened) then
        open(unit=outputFileUnit,file="cirrusOutput.txt")
      else
        call finish("","Could not open cirrusOutput.txt, unit already used...")
      end if
      
      debugOutput = outputFileUnit
    
    end if
    
#endif
  end subroutine

  subroutine ConstructCirrusModelStreams()
#ifdef CIRRUS_MODULE_DEBUG
    
    ! call this from 
    ! init_subm_memory in mo_submodel_interface or
    ! construct_activ_stream in mo_activ
    
    ! construct cirrus output stream to save all relevant variables
  
    USE mo_memory_base, ONLY: add_stream_element
    use mo_radiation_parameters, only: lradforcing
    
    implicit none
    
    if(lDebug) then
      
      write(outputFileUnit,*) "Building output streams"
    
      call CreateNewOutputStream(cirrusStream, 'cirrus')
    
      call AddNormalSE (cirrusStream, geopotHeightStream, 'geopotHeight', longname='Geopotential Height',units='m2 s-2' )
      call AddNormalSE (cirrusStream, presStream, 'pres', longname='Pressure',units='hPa' )
      call AddNormalSE (cirrusStream, tempStream, 'temp', longname='Temperature',units='K' )
      call AddNormalSE (cirrusStream, updraftStream, 'updraft', longname='Updraft velocity',units='cm/s' )
      call AddNormalSE (cirrusStream, satStream, 'sat', longname='Saturation',units='1' )
      call AddNormalSE (cirrusStream, satAfterStream, 'satAfter', longname='Saturation After Cirrus Scheme',units='1' )

      call AddTimedValueSE (cirrusStream, satAboveHomoSatCrit, 'satAboveHomoSatCrit')
      call AddTimedValueSE (cirrusStream, satClearSkyStream, 'satClearSky')
      call AddTimedValueSE (cirrusStream, satCloudySkyStream, 'satCloudySky')
      
      call AddConcRadSE (cirrusStream, homoAeroInStream, 'homoAeroIn', &
      longname='Homogeneous Aerosol In', timed = .false.)
      call AddConcRadSE (cirrusStream, homoAeroOutStream, 'homoAeroOut', &
      longname='Homogeneous Aerosol Out', timed = .false.)
      call AddConcRadSE (cirrusStream, homoIceStream, 'homoIce', longname='Homogeneous Ice')
      call AddNormalSE (cirrusStream, homoIceFreezingTime, 'homoIceTimeFreezing', &
      longname='Homogeneous Ice Freezing Time',units='s')
      
      call AddConcRadSE (cirrusStream, dustImmAeroInStream, 'dustImmAeroIn', &
      longname='Dust immersion Aerosol In', timed = .false.)
      call AddConcRadSE (cirrusStream, dustImmAeroOutStream, 'dustImmAeroOut', &
      longname='Dust immersion Aerosol Out', timed = .false.)
      call AddConcRadSE (cirrusStream, dustImmIceStream, 'dustImmIce', longname='Dust immersion Ice')
      call AddNormalSE (cirrusStream, dustImmIceFreezingTime, 'dustImmIceTimeFreezing', &
      longname='DustImm Ice Freezing Time',units='s')
      call AddNormalSE (cirrusStream, dustImmVolFractionSE, 'dustImmVolFraction')
      call AddConcRadSE (cirrusStream, dustImmAccAeroInStream, 'dustImmAccAeroIn', timed = .false.)
      call AddConcRadSE (cirrusStream, dustImmCoaAeroInStream, 'dustImmCoaAeroIn', timed = .false.)
      
      call AddConcRadSE (cirrusStream, dustDepAAeroInStream, 'dustDepAAeroIn', &
      longname='Dust deposition accumulation mode Aerosol In', timed = .false.)
      call AddConcRadSE (cirrusStream, dustDepAAeroOutStream, 'dustDepAAeroOut', &
      longname='Dust deposition accumulation mode Aerosol Out', timed = .false.)
      call AddConcRadSE (cirrusStream, dustDepAIceStream, 'dustDepAIce', longname='Dust deposition accumulation mode Ice')
      call AddNormalSE (cirrusStream, dustDepAIceFreezingTime, 'dustDepAIceTimeFreezing', &
      longname='DustDep Ice accumulation mode Freezing Time',units='s')
      
      call AddConcRadSE (cirrusStream, dustDepCAeroInStream, 'dustDepCAeroIn', &
      longname='Dust deposition coarse mode Aerosol In', timed = .false.)
      call AddConcRadSE (cirrusStream, dustDepCAeroOutStream, 'dustDepCAeroOut', &
      longname='Dust deposition coarse mode Aerosol Out', timed = .false.)
      call AddConcRadSE (cirrusStream, dustDepCIceStream, 'dustDepCIce', longname='Dust deposition coarse mode Ice')
      call AddNormalSE (cirrusStream, dustDepCIceFreezingTime, 'dustDepCIceTimeFreezing', &
      longname='DustDep Ice coarse mode Freezing Time',units='s')
      
      call AddConcRadSE (cirrusStream, homoAeroAllStream, 'homoAeroAll', timed = .false.)
      call AddConcRadSE (cirrusStream, dustImmAeroAllStream, 'dustImmAeroAll', timed = .false.)
      call AddConcRadSE (cirrusStream, dustDepAeroAllStream, 'dustDepAeroAll', timed = .false.)
      call AddConcRadSE (cirrusStream, dustDepIceAllStream, 'dustDepIceAll')
      
      call AddConcRadSE (cirrusStream, stratiformIceStream, 'iceStratiform', longname='Stratiform Ice', timed = .false.)
      call AddNormalSE (cirrusStream, stratiformIceRadOutStream, 'iceStratiformRadOut', &
      longname='Stratiform Ice Out Radius',units='µm' )
      call AddConcRadSE (cirrusStream, detrainedIceStream, 'iceDetrained', longname='Detrained Ice')
      call AddNormalSE (cirrusStream, detrainedIceRadOutStream, 'iceDetrainedRadOut', &
      longname='Detrained Ice Out Radius',units='µm' )
      call AddConcRadSE (cirrusStream, stratiformIceLarger150, 'iceZStratiformLarger150')
      
      call AddConcRadSE (cirrusStream, newIceCrystalsStream, 'newIceCrystals', longname='New Formed Ice Crystals')
      call AddNormalSE (cirrusStream, depositedWaterVaporStream, 'depositedWaterVapor', units='g/kg')
      
      call AddNormalSE (cirrusStream, calculateCirrusFreqStream, 'calculateCirrusFreq')

      call AddConcRadSE (cirrusStream, iceByOvershoot, 'iceZZOvershoot')
      call AddTimedValueSE (cirrusStream, tempBelow238, 'tempBelow238')
      
      call add_stream_element (cirrusStream, 'ziceCloudTime', iceCloudTime)
      call add_stream_element (cirrusStream, 'ziceCloud_IWC', iceCloudIWC, units='mg/m3')
      call add_stream_element (cirrusStream, 'ziceCloud_EffRad', iceCloudEffRad, units='µm')
      call add_stream_element (cirrusStream, 'ziceCloud_ICNC', iceCloudICNC, units='l-1')
      call add_stream_element (cirrusStream, 'ziceCloud_Cover', iceCloudCover)
      call add_stream_element (cirrusStream, 'ziceCloud_OpticalDepth', iceCloudOpticalDepth)
      
      call add_stream_element (cirrusStream, 'zcirrusCod3Time', cirrusCod3Time)
      call add_stream_element (cirrusStream, 'zcirrusCod3_IWC', cirrusCod3IWC, units='mg/m3')
      call add_stream_element (cirrusStream, 'zcirrusCod3_EffRad', cirrusCod3EffRad, units='µm')
      call add_stream_element (cirrusStream, 'zcirrusCod3_ICNC', cirrusCod3ICNC, units='l-1')
      call add_stream_element (cirrusStream, 'zcirrusCod3_Cover', cirrusCod3Cover)
      call add_stream_element (cirrusStream, 'zcirrusCod3_OpticalDepth', cirrusCod3OpticalDepth)
      
      if(lSeeding) then
        call AddConcRadSE (cirrusStream, seedAeroInStream, 'seedAeroIn', longname='Seeding Aerosol In', timed = .false.)
        call AddConcRadSE (cirrusStream, seedAeroOutStream, 'seedAeroOut', &
        longname='Seeding deposition Aerosol Out', timed = .false.)
        call AddConcRadSE (cirrusStream, seedIceStream, 'seedIce', longname='Seeding Ice')
        call AddNormalSE (cirrusStream, seedIceFreezingTime, 'seedIceTimeFreezing', longname='Seed Ice Freezing Time',units='s')
      end if
           
      if(lradforcing(1)) then
        call add_stream_element (cirrusStream, 'radForcingSW', radForcingSW, &
        longname='Ice cloud shortwave radiative forcing TOA', units='W m-2')
      end if
      if(lradforcing(2)) then
        call add_stream_element (cirrusStream, 'radForcingLW', radForcingLW, &
        longname='Ice cloud longwave radiative forcing TOA', units='W m-2')
      end if
      if(lradforcing(1) .and. lradforcing(2)) then
        call add_stream_element (cirrusStream, 'radForcingNet', radForcingNet, &
        longname='Ice cloud net radiative forcing TOA', units='W m-2')
      end if
      
      if(lradforcing(1)) then
        call add_stream_element (cirrusStream, 'radForcingColdBaseSW', radForcingColdBaseSW, &
        longname='Ice cloud shortwave radiative forcing TOA', units='W m-2')
        call add_stream_element (cirrusStream, 'radForcingColdBase2SW', radForcingColdBase2SW, &
        longname='Ice cloud shortwave radiative forcing TOA', units='W m-2')
        call add_stream_element (cirrusStream, 'radForcingWarmBaseSW', radForcingWarmBaseSW, &
        longname='Ice cloud shortwave radiative forcing TOA', units='W m-2')
        call add_stream_element (cirrusStream, 'radForcingWarmBase2SW', radForcingWarmBase2SW, &
        longname='Ice cloud shortwave radiative forcing TOA', units='W m-2')
      end if
      if(lradforcing(2)) then
        call add_stream_element (cirrusStream, 'radForcingColdBaseLW', radForcingColdBaseLW, &
        longname='Ice cloud longwave radiative forcing TOA', units='W m-2')
        call add_stream_element (cirrusStream, 'radForcingColdBase2LW', radForcingColdBase2LW, &
        longname='Ice cloud longwave radiative forcing TOA', units='W m-2')
        call add_stream_element (cirrusStream, 'radForcingWarmBaseLW', radForcingWarmBaseLW, &
        longname='Ice cloud longwave radiative forcing TOA', units='W m-2')
        call add_stream_element (cirrusStream, 'radForcingWarmBase2LW', radForcingWarmBase2LW, &
        longname='Ice cloud longwave radiative forcing TOA', units='W m-2')
      end if
      if(lradforcing(1) .and. lradforcing(2)) then
        call add_stream_element (cirrusStream, 'radForcingColdBaseNet', radForcingColdBaseNet, &
        longname='Ice cloud net radiative forcing TOA', units='W m-2')
        call add_stream_element (cirrusStream, 'radForcingColdBase2Net', radForcingColdBase2Net, &
        longname='Ice cloud net radiative forcing TOA', units='W m-2')
        call add_stream_element (cirrusStream, 'radForcingWarmBaseNet', radForcingWarmBaseNet, &
        longname='Ice cloud net radiative forcing TOA', units='W m-2')
        call add_stream_element (cirrusStream, 'radForcingWarmBase2Net', radForcingWarmBase2Net, &
        longname='Ice cloud net radiative forcing TOA', units='W m-2')
      end if
          
    end if
    
#endif
  end subroutine
  
  subroutine SaveRadiativeForcing(bs, bn, shortwave, longwave, iDiag)
    
    ! call this from calculate_forcing in mo_radiation_forcing
    
    use mo_radiation_parameters, only: lradforcing
    USE mo_time_control, ONLY: dTime=>delta_time
    
    implicit none
    
    integer, intent(in) :: bs,bn ! multi processing position in output streams: bs = block size, bn = block number
    real(wp), intent(in) :: shortwave(:), longwave(:)
    integer, intent(in) :: iDiag
    
#ifdef CIRRUS_MODULE_DEBUG
    
    if(lDebug) then
      
      if(iDiag == 1) then
        
        if(lradforcing(1)) then
          radForcingSW(1:bs,bn) = radForcingSW(1:bs,bn) + shortwave(1:bs)*dTime
        end if
        if(lradforcing(2)) then
          radForcingLW(1:bs,bn) = radForcingLW(1:bs,bn) + longwave(1:bs)*dTime
        end if
        if(lradforcing(1) .and. lradforcing(2)) then
          radForcingNet(1:bs,bn) = radForcingNet(1:bs,bn) + (shortwave(1:bs) + longwave(1:bs))*dTime
        end if
        
      else if(iDiag == 2) then
        
        if(lradforcing(1)) then
          radForcingColdBaseSW(1:bs,bn) = radForcingColdBaseSW(1:bs,bn) + shortwave(1:bs)*dTime
        end if
        if(lradforcing(2)) then
          radForcingColdBaseLW(1:bs,bn) = radForcingColdBaseLW(1:bs,bn) + longwave(1:bs)*dTime
        end if
        if(lradforcing(1) .and. lradforcing(2)) then
          radForcingColdBaseNet(1:bs,bn) = radForcingColdBaseNet(1:bs,bn) + (shortwave(1:bs) + longwave(1:bs))*dTime
        end if
        
      else if(iDiag == 3) then
      
        if(lradforcing(1)) then
          radForcingColdBase2SW(1:bs,bn) = radForcingColdBase2SW(1:bs,bn) + shortwave(1:bs)*dTime
        end if
        if(lradforcing(2)) then
          radForcingColdBase2LW(1:bs,bn) = radForcingColdBase2LW(1:bs,bn) + longwave(1:bs)*dTime
        end if
        if(lradforcing(1) .and. lradforcing(2)) then
          radForcingColdBase2Net(1:bs,bn) = radForcingColdBase2Net(1:bs,bn) + (shortwave(1:bs) + longwave(1:bs))*dTime
        end if
        
      else if(iDiag == 4) then
    
        if(lradforcing(1)) then
          radForcingWarmBaseSW(1:bs,bn) = radForcingWarmBaseSW(1:bs,bn) + shortwave(1:bs)*dTime
        end if
        if(lradforcing(2)) then
          radForcingWarmBaseLW(1:bs,bn) = radForcingWarmBaseLW(1:bs,bn) + longwave(1:bs)*dTime
        end if
        if(lradforcing(1) .and. lradforcing(2)) then
          radForcingWarmBaseNet(1:bs,bn) = radForcingWarmBaseNet(1:bs,bn) + (shortwave(1:bs) + longwave(1:bs))*dTime
        end if
        
      else if(iDiag == 5) then
  
        if(lradforcing(1)) then
          radForcingWarmBase2SW(1:bs,bn) = radForcingWarmBase2SW(1:bs,bn) + shortwave(1:bs)*dTime
        end if
        if(lradforcing(2)) then
          radForcingWarmBase2LW(1:bs,bn) = radForcingWarmBase2LW(1:bs,bn) + longwave(1:bs)*dTime
        end if
        if(lradforcing(1) .and. lradforcing(2)) then
          radForcingWarmBase2Net(1:bs,bn) = radForcingWarmBase2Net(1:bs,bn) + (shortwave(1:bs) + longwave(1:bs))*dTime
        end if
        
      end if
      
    end if
    
#endif
  end subroutine
    
  !<< --- END used for lDebug = true --------------------------------------------------------------------------
  
end module
