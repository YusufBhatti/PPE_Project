! submodelctl.inc
!
! This namelist defines general submodel switches which are needed in the interface
! layer (*mo_submodel_interface* or other parts of the standard ECHAM code). Also
! included are switches which define the coupling between various submodels (for
! example in HAMMOZ, which coupled aerosol and gas-phase chemical processes). Other
! submodel-specific switches should be defined and maintained in extra namelists
! which carry the name of the submodel itself (e.g. mozctl or hamctl).
!
!     ------------------------------------------------------------------
NAMELIST /submodelctl/ &
  lxt,              &  ! switch generic test tracer submodule on/off
  lmethox,          &  ! switch for upper atmospheric H2O production from methane
  ltransdiag,       &  ! switch to turn on atmospheric energy transport diagnostics
  lco2,             &  ! switch for CO2 submodel (JSBACH related)
  lham,             &  ! switch HAM aerosol module on/off
  lmoz,             &  ! switch MOZART on/off
  lhammoz,          &  ! switch HAM and MOZ on/off together with the coupling between the two
                       ! note: lhammoz overrides lham and lmoz
  lhammonia,        &  ! switch HAMMONIA on/off
  llght,            &  ! switch lightning emissions on/off
  lbioemi_stdalone, &  ! switch biogenic emissions model as a standalone submodel
                       ! (ie not embedded in HAM or MOZ) on/off
  losat,            &  ! satellite simulator on/off
  loisccp,          &  ! isccp diagnostics on/off
  lhmzphoto,        &  ! hammoz photolysis frequency coupling on/off
  lhmzoxi,          &  ! hammoz coupling of oxidant fields on/off
  lhmzhet,          &  ! hammoz heterogeneous chemistry coupling on/off
  lchemfeedback,    &  ! combi-switch for interactive chemistry:
                       ! for moz: lchemrad and linterh2o = true
                       ! for hammonia: as above plus lchemheat, atmospheric mass and cp
  lchemrad,         &  ! chemistry interacts with radiation on/off
  linterh2o,        &  ! feedback water content from MOZ to ECHAM and vice versa
  lchemheat,        &  ! chemical heating on/off
  lccnclim,         &  ! activate CCN climatology as submodel
  lflighttrack,     &  ! flight-track simulator
  linteram,         &  ! hammonia air mass from chemistry on/off
  lintercp,         &  ! hammonia specific heat from chemistry on/off
  lemissions,       &  ! switch emissions on/off
  lchemistry,       &  ! switch chemistry calculations on/off
  ldrydep,          &  ! switch dry deposition on/off
  lwetdep,          &  ! switch wet deposition on/off
  lsedimentation,   &  ! switch sedimentation on/off
  laero_micro,      &  ! switch aerosol microphysical processes on/off
  lburden,          &  ! activate burden (column integral) diagnostics for mass mixing ratio tracers
  emi_basepath,     &  ! path to emission files specified in emi_spec.txt
  emi_scenario,     &  ! RCP (Representative Concentration Pathway) to be modelled
  laoa                 ! age-of-air tracer submodel
!     ------------------------------------------------------------------

