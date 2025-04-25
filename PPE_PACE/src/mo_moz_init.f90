!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_init
!!
!! \brief
!!  Initialisation routines for MOZART chemistry
!!
!! \author Martin G. Schultz (FZ-Juelich)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - MGS: original version (2009-09-22)
!!  - MGS (2013-07-25): revision for MOZ 0.9
!!
!! \limitations
!!  none
!!
!! \details
!!  This module contains the general MOZ initialisation routine which is called 
!!  from init_subm.
!!  It first calls setmoz to initialize moz module variables and read the mozctl namelist.
!!  Then it initializes various MOZ components and defines the MOZ species and tracers.
!!
!! \bibliographic_references
!!  none
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE mo_moz_init

  USE mo_kind,            ONLY : dp

  IMPLICIT NONE

  PRIVATE

  SAVE

  PUBLIC  :: start_moz, moz_define_tracer, moz_initialize, init_moz_ic

! Variable declarations

  TYPE mozsp           ! MOZART species properties
    CHARACTER(len=24)   :: name
    REAL(dp)            :: henry(2)      ! 1=H(298K), 2=B; H = H(298K)*exp(B/T-B/298.), units mole atm-1
    REAL(dp)            :: f0
  END TYPE


  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! start_moz: high-level MOZ initialisation routine; interface for MOZART chemistry 
!! initialisation including species definition
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Externals:
!! <ol>
!! <li>none
!! </ol>
!!
!! @par Notes
!! 
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!++mgs #249: added lchemistry flag
SUBROUTINE start_moz (lchemistry, lcouple_oxidants)
!--mgs

  USE mo_kind,            ONLY : dp
  USE mo_exception,       ONLY : message, message_text, finish, em_error, em_debug
!++mgs #249: added lphotolysis flag
  USE mo_moz,             ONLY : setmoz, lphotolysis, uvalbedo_file, photovars, nphotovars
  USE mo_moz_eslookup,    ONLY : esinti
  USE mo_moz_mods,        ONLY : pcnstm1, chem_mods_inti  ! preprocessor generated!
  USE mo_moz_subs,        ONLY : set_sim_dat              ! preprocessor generated!
  USE mo_moz_tropo_rates, ONLY : tropo_rates_inti
  USE mo_moz_uvalbedo,    ONLY : moz_uvalbedo_init
  USE mo_moz_photo,       ONLY : photo_initialize
  USE mo_hammoz_het_tropo_rates, ONLY: bc_aerosol

!++mgs #249
  LOGICAL, INTENT(in)    :: lchemistry         ! run with chemistry
!--mgs
  LOGICAL, INTENT(in)    :: lcouple_oxidants   ! couple (HAM) oxidant fields (species)

! -- parameters for esinti 
  REAL(dp), PARAMETER   ::   &
       epsilo = 0.622_dp, &
       latvap = 2.5104e06_dp, &
       latice = 3.336e5_dp, &
       rh2o   = 461._dp, &
       cpair  = 1004.64_dp   
  LOGICAL, PARAMETER ::  &
       esice = .TRUE.

  CALL message('','Beginning execution of start_moz', level=em_debug)
! -- 0. sanity check: at least one "transported species" must be defined in MOZART mechanism
  if (pcnstm1 == 0) then
     write(message_text,'(a40)') 'no transported species defined (pcnstm1==0). Defunct MOZ mechanism!!'
     CALL message('start_moz', message_text, level=em_error)
     CALL finish('start_moz: Critical error - please ask Martin!')
  endif


! -- 1. set default values (incl. mo_moz_constants) and read mozctl namelist
!++mgs #249: added lchemistry argument
  CALL setmoz(lchemistry,bc_aerosol)
!--mgs

! -- 2. initialize MOZART components: simulation parameters, saturation vapour pressure,
!                                     module initialisation, constants
  CALL chem_mods_inti ()
  CALL message('start_moz','finished chem_mods_inti', level=em_debug)
  CALL set_sim_dat ()
  CALL message('start_moz','finished set_sim_dat', level=em_debug)
  CALL esinti( epsilo, latvap, latice, rh2o, cpair, esice )
  CALL message('start_moz','finished esinti', level=em_debug)

!!+++  IF (.NOT. ALLOCATED(latwts)) ALLOCATE(latwts(plat))
!!+++  IF (.NOT. ALLOCATED(phi)) ALLOCATE(phi(plat))
!!+++  IF (.NOT. ALLOCATED(lam)) ALLOCATE(lam(plon))

!!+++  latwts(1:plat) = gl_gw(1:plat)
!!+++  phi(1:plat)    = twopi*philat(1:plat)/360
!!+++  lam(1:plon)    = twopi*philon(1:plon)/360

!-- initialisation of UV albedo and photolysis scheme
!++mgs #249: if block
  IF (lphotolysis) THEN
     CALL moz_uvalbedo_init (uvalbedo_file)
     CALL message('start_moz','finished moz_uvalbedo_init', level=em_debug)
     CALL photo_initialize(photovars, nphotovars)
  END IF
!--mgs

!-- MOZART species definition
!   define all species in tracnam as ECHAM species
  CALL moz_species(lcouple_oxidants)

  CALL tropo_rates_inti ()
  CALL message('start_moz','finished tropo_rates_inti', level=em_debug)

  CALL message('start_moz','done.', level=em_debug)

END SUBROUTINE start_moz

!
!  moz_species_properties: Set gas-phase species properties which are not defined by
!  MOZART preprocessor.
!
!

SUBROUTINE moz_species_properties (kntrx, pmozprop)

  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish, message_text

  IMPLICIT NONE

! arguments
  INTEGER,          INTENT(out) :: kntrx
!++mgs
  TYPE(mozsp),      INTENT(out) :: pmozprop(:)
!--mgs
! local variables
  INTEGER                       :: ntrx          ! actual number of values defined below
  INTEGER                       :: kdim          ! dimension of pmozprop

  !ntrx = 128
  ntrx = 134
  kdim = SIZE(pmozprop(:)%f0)
! test if dimension is large enough
  IF (kdim < ntrx) THEN
    WRITE(message_text,'(a,i0,a,i0)') 'Internal error: field dimension too small. Should be ', ntrx, &
                          ', is ', kdim
    CALL finish('moz_species_properties', message_text)
  END IF

! define list of tracers and henry constants and f0 values
! species names according to MOZART preprocessor input
! JPL-report 02-25 (February, 2003)

!++mgs: NEW LIST
  ! ## ToDo ## : f0 values need review (species property "dryreac" - see lg_drydep scheme)

  !                     name        henry(1)      henry(2)        f0     
  ! section 1: inorganic and organic species
  pmozprop(  1) = mozsp( 'O3'        , (/1.03e-2_dp  , 2830._dp/)     , 1.0_dp  )  ! JPL15(2006)
  pmozprop(  2) = mozsp( 'H2'        , (/7.8e-4_dp   , 500._dp/)      , 0._dp   )  ! R. Sander (1999)
  pmozprop(  3) = mozsp( 'H2O2'      , (/8.44e4_dp   , 7600._dp/)     , 1.0_dp  )  ! JPL17(2011)
  pmozprop(  4) = mozsp( 'HO2'       , (/690._dp     , 0._dp/)        , 1.0_dp   )  ! JPL15(2006)
  pmozprop(  5) = mozsp( 'HO2NO2'    , (/1.2e4_dp    , 6900._dp/)     , 1.0_dp  )  ! R. Sander (1999)
  pmozprop(  6) = mozsp( 'HNO3'      , (/3.2e11_dp   , 8700._dp/)     , 1.0_dp  )  ! H=2.1E5 M/atm and Ka = 15.4 M (Schwartz and White 1981). At an average cloud droplet pH=5 => H* = H * (1+Ka/[H+]) = 2.1E5*(1+15.4/1E-5) = 3.2E11 
  pmozprop(  7) = mozsp( 'HONO'      , (/5.05e3_dp   , 4800._dp/)     , 0.1_dp  )  ! R. Sander (1999): H=49 and Ka=5.1E-4 => H* = 49(1+5.1E-4*1E7) = 5.05E3
  pmozprop(  8) = mozsp( 'NO'        , (/1.92e-3_dp  , 1790._dp/)     , 1.0_dp  )  ! JPL15(2006), comment
  pmozprop(  9) = mozsp( 'NO2'       , (/12.e-2_dp   , 2360._dp/)     , 1.0_dp  )  ! JPL17(2011)
  pmozprop( 10) = mozsp( 'NO3'       , (/3.8e-2_dp   , 0._dp/)        , 1.0_dp  )  ! JPL15(2006)
!  pmozprop( 11) = mozsp( 'N2O5'      , (/2.1_dp      , 3400._dp/)     , 1.0_dp  )  ! R. Sander (1999): 2 ref. give it as infinite and one as 2.1 M/atm
  pmozprop( 11) = mozsp( 'N2O5'      , (/3.2e11_dp   , 8700._dp/)     , 1.0_dp  )  ! HNO3 as proxy as done in GEOS-Chem
  ! Section: organic compounds (no deposition of organic RO2) 
  ! ROOH are assumed to have a surface reactivity of 1, in all other cases when oxygen present f0 = 0.1
  ! -OOH groups are considered similar to -OH in affecting the Henry's law constant
  !                1C
  pmozprop( 12) = mozsp( 'CH4'       , (/1.41e-3_dp  , 2040._dp/)     , 0._dp   )  ! JPL17(2011) SMILES: C
  pmozprop( 13) = mozsp( 'CO'        , (/9.81e-4_dp  , 1720._dp/)     , 0._dp   )  ! JPL17(2011) SMILES: [C-]#[O+] 
  pmozprop( 14) = mozsp( 'CH2O'      , (/3.23e3_dp   , 7100._dp/)     , 0.1_dp  )  ! JPL17(2011) SMILES: C=O
  pmozprop( 15) = mozsp( 'CH3OH'     , (/203._dp     , 5640._dp/)     , 0._dp   )  ! JPL17(2011) SMILES: CO
  pmozprop( 16) = mozsp( 'CH3OOH'    , (/300._dp     , 5280._dp/)     , 0.1_dp  )  ! JPL17(2011) SMILES: COO
  pmozprop( 17) = mozsp( 'HCN'       , (/12._dp      , 5000._dp/)     , 1.0_dp  )  ! R. Sander (1999) SMILES: C#N
  pmozprop( 18) = mozsp( 'HCOOH'     , (/8.9e3_dp    , 6100._dp/)     , 0._dp   )  ! R. Sander (1999): measured value by Johnson et al. 1996, surface reactivity as in Nguyen et al. PNAS 2015 SMILES: C(=O)O
  pmozprop( 19) = mozsp( 'CH3O2NO2'  , (/2.0_dp      , 4700._dp/)     , 0.1_dp  )  ! R. Sander (1999): methyl nitrate as proxy species SMILES: COO[N+](=O)[O-]
  !                2C
  pmozprop( 20) = mozsp( 'C2H5OH'    , (/190._dp     , 6660._dp/)     , 0._dp   )  ! JPL17(2011) SMILES: CCO
  pmozprop( 21) = mozsp( 'C2H5OOH'   , (/336._dp     , 5995._dp/)     , 1.0_dp  )  ! JPL17(2011) SMILES: CCOO
  pmozprop( 22) = mozsp( 'CH3CHO'    , (/12.9_dp     , 5890._dp/)     , 0.1_dp  )  ! JPL17(2011) SMILES: CC=O
  pmozprop( 23) = mozsp( 'CH3CN'     , (/52.8_dp     , 3970._dp/)     , 0._dp   )  ! JPL17(2011) SMILES: CC#N
  pmozprop( 24) = mozsp( 'PAN'       , (/2.8_dp      , 5730._dp/)     , 0.1_dp  )  ! JPL17(2011) SMILES: CC(=O)OON(=O)=O
  pmozprop( 25) = mozsp( 'CH3COOH'   , (/4.1e3_dp    , 6200._dp/)     , 0._dp   )  ! JPL17(2011) SMILES: CC(=O)O
  pmozprop( 26) = mozsp( 'CH3COOOH'  , (/837._dp     , 5310._dp/)     , 0.1_dp  )  ! JPL17(2011) SMILES: C(C(=O)O)O
  pmozprop( 27) = mozsp( 'GLYALD'    , (/4.1e4_dp    , 4600._dp/)     , 0.1_dp  )  ! R. Sander (1999) : 2-hydroxyethanal SMILES: C(C=O)O
  pmozprop( 28) = mozsp( 'GLYOXAL'   , (/4.19e5_dp   , 7480._dp/)     , 0.1_dp  )  ! JPL17(2011) : Heff SMILES: C(=O)C=O 
  pmozprop( 29) = mozsp( 'HCOCO2H'   , (/1.1e4_dp    , 4800._dp/)     , 0.1_dp  )  ! Sander (2015) : glyoxilic acid SMILES: O=CC(=O)O
  pmozprop( 30) = mozsp( 'HCOCO3H'   , (/1.1e4_dp    , 4800._dp/)     , 0.1_dp  )  ! Sander (2015) : glyoxilic acid SMILES: OOC(=O)C=O
  pmozprop( 31) = mozsp( 'HOCH2CO2H' , (/2.8e4_dp    , 4000._dp/)     , 0.1_dp  )  ! Sander (2015) : glycolic acid SMILES: OCC(=O)O
  pmozprop( 32) = mozsp( 'HOCH2CO3H' , (/2.8e4_dp    , 4000._dp/)     , 0.1_dp  )  ! Sander (2015) : glycolic acid SMILES: OCC(=O)OO
  !                3C
  pmozprop( 33) = mozsp( 'C3H7OOH'   , (/336._dp     , 5995._dp/)     , 1.0_dp  )  ! ** as C2H5OOH SMILES: OOC(C)C 
  pmozprop( 34) = mozsp( 'CH3COCH3'  , (/27.8_dp     , 5530._dp/)     , 0._dp   )  ! JPL17(2011) SMILES: CC(=O)C
  pmozprop( 35) = mozsp( 'CH3COCHO'  , (/3.7e3_dp    , 7500._dp/)     , 0.1_dp  )  ! R. Sander (1999) : propanonal SMILES: O=CC(=O)C
  pmozprop( 36) = mozsp( 'HYAC'      , (/7.7e3_dp    ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) SMILES: CCC(=O)O
  pmozprop( 37) = mozsp( 'NOA'       , (/1.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) SMILES: CC(=O)CON(=O)=O
  !                4C
  pmozprop( 38) = mozsp( 'MACR'      , (/6.5_dp      , 5300._dp/)     , 0.1_dp  )  ! R. Sander (1999) : 2-methylpropenal SMILES: CC(=C)C=O
  pmozprop( 39) = mozsp( 'MACROOH'   , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: OCC(C)(OO)C=O
  pmozprop( 40) = mozsp( 'MACROH'    , (/2.1e5_dp    , 9900._dp/)     , 0._dp   )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: CC(O)(CO)C=O
  pmozprop( 41) = mozsp( 'MEK'       , (/37._dp      , 8200._dp/)     , 0.1_dp  )  ! Sander (2015) SMILES: CCC(=O)C
  pmozprop( 42) = mozsp( 'MEKOOH'    , (/1.1e3_dp    , 7300._dp/)     , 1.0_dp  )  ! Sander (2015) : 2-butanol as proxy SMILES: C/C=C/C(=O)O
  pmozprop( 43) = mozsp( 'MPAN'      , (/1.7_dp      , 0._dp/)        , 0.1_dp  )  ! R. Sander (1999) SMILES: O=N(=O)OOC(=O)C(=C)C
  pmozprop( 44) = mozsp( 'MVK'       , (/4.1e1_dp    , 7800._dp/)     , 0.1_dp  )  ! R. Sander (1999) : 3-buten-2-one SMILES: CC(=O)C=C
  pmozprop( 45) = mozsp( 'MACO2H'    , (/1.5e3_dp    , 6800._dp/)     , 0.1_dp  )  ! Sander (2015) : methyl butanoic acid as proxy SMILES: CC(=C)C(=O)O
  pmozprop( 46) = mozsp( 'MACO3H'    , (/1.5e3_dp    , 6800._dp/)     , 1.0_dp  )  ! Sander (2015) : methyl butanoic acid as proxy SMILES: CC(=C)C(=O)OO
  pmozprop( 47) = mozsp( 'BIGALD1'   , (/2.5e5_dp    ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : TOL_EPOX as proxy MCM name MALDIAL SMILES: O=CC=CC=O  
  pmozprop( 48) = mozsp( 'CO2H3CHO'  , (/4.1e4_dp    , 4600._dp/)     , 0.1_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy SMILES: O=CC(O)C(=O)C
  pmozprop( 49) = mozsp( 'CO2H3CO3H' , (/4.1e4_dp    , 4600._dp/)     , 1.0_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy SMILES: CC(=O)C(O)C(=O)OO
  pmozprop( 50) = mozsp( 'BIACETOH'  , (/4.1e4_dp    , 4600._dp/)     , 0.1_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy SMILES: CC(=O)C(=O)CO
  pmozprop( 51) = mozsp( 'IBUTALOH'  , (/4.1e4_dp    , 4600._dp/)     , 0.1_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy SMILES: O=CC(C)(C)O
  pmozprop( 52) = mozsp( 'IBUTALOHOOH' , (/4.1e4_dp    , 4600._dp/)     , 1.0_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy MCM species IPRHOCO3H SMILES: OOC(=O)C(C)(C)O   
  pmozprop( 53) = mozsp( 'LHMVKABOOH', (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: CC(=O)C(O)COO 
  pmozprop( 54) = mozsp( 'MVKN'      , (/5.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! like for ISOPNO3 from Nguyen et al. 2015 SMILES: OCC(ON(=O)=O)C(=O)C 
  pmozprop( 55) = mozsp( 'MACRN'     , (/5.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! like for ISOPNO3 from Nguyen et al. 2015 SMILES: OCC(C)(C=O)ON(=O)=O
  !                5C
  pmozprop( 56) = mozsp( 'MBO'        , (/64._dp      ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : 2-methyl-3-buten-2-ol as proxy SMILES: C=CC(C)(C)O 
  pmozprop( 57) = mozsp( 'HPALD'      , (/4.1e4_dp    , 4600._dp/)     , 1.0_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy SMILES: C/C(=C/C=O)/COO
  pmozprop( 58) = mozsp( 'PACALD'     , (/4.1e4_dp    , 4600._dp/)     , 1.0_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy SMILES: OOC(=O)/C=C(\C=O)/C
  pmozprop( 59) = mozsp( 'ISOPBNO3'   , (/5.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! Nguyen et al. 2015 SMILES: OCC(C)(C=C)ON(=O)=O
  pmozprop( 60) = mozsp( 'ISOPDNO3'   , (/5.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! Nguyen et al. 2015 SMILES: OCC(ON(=O)=O)C(=C)C
  pmozprop( 61) = mozsp( 'LISOPACNO3' , (/5.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! Nguyen et al. 2015 SMILES: OCC=C(C)CON(=O)=O
  pmozprop( 62) = mozsp( 'LC5PAN1719' , (/5.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! Nguyen et al. 2015 SMILES: OCC(=CC(=O)OON(=O)=O)C
  pmozprop( 63) = mozsp( 'LHC4ACCO2H' , (/2.1e5_dp    , 9900._dp/)     , 0.1_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: OCC(=CC(=O)O)C
  pmozprop( 64) = mozsp( 'LHC4ACCO3H' , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: OCC(=CC(=O)OO)C
  pmozprop( 65) = mozsp( 'LIECHO'     , (/4.1e4_dp    , 4600._dp/)     , 0.1_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy SMILES: CC1(OC1CO)C=O
  pmozprop( 66) = mozsp( 'LIECO3H'    , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: CC(O)(C1CO1)C(=O)OO
  pmozprop( 67) = mozsp( 'LIEPOX'     , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: OCC1OC1(C)CO
  pmozprop( 68) = mozsp( 'LNISOOH'    , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: O=CC(O)C(C)(OO)CON(=O)=O
  pmozprop( 69) = mozsp( 'MBOOOH'     , (/3.e11_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2,3-butanetriol as proxy SMILES: OCC(OO)C(C)(C)O
  pmozprop( 70) = mozsp( 'NC4CHO'     , (/2.3_dp      ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : 2-methylbutanal as proxy SMILES: O=CC=C(C)CON(=O)=O
  pmozprop( 71) = mozsp( 'NISOPOOH'   , (/3.6e4_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : 5-nitrooxy-2-pentanol as proxy SMILES: OOCC=C(C)CON(=O)=O
  pmozprop( 72) = mozsp( 'LHC4ACCHO'  , (/4.1e4_dp    , 4600._dp/)     , 0.1_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy SMILES: CC(=CC=O)CO
  pmozprop( 73) = mozsp( 'HCOC5'      , (/4.1e4_dp    , 4600._dp/)     , 1.0_dp  )  ! R. Sander (1999) : 2-hydroxyethanal as proxy SMILES: CC(=C)C(=O)CO
  pmozprop( 74) = mozsp( 'ISOPBOOH'   , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: OCC(C)(OO)C=C
  pmozprop( 75) = mozsp( 'ISOPDOOH'   , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: OCC(OO)C(=C)C
  pmozprop( 76) = mozsp( 'LISOPACOOH' , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: OOCC=C(C)CO
  pmozprop( 77) = mozsp( 'ISOPBOH'    , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: CC(O)(CO)C=C
  pmozprop( 78) = mozsp( 'ISOPDOH'    , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: CC(=C)C(O)CO
  pmozprop( 79) = mozsp( 'ISOPAOH'    , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: OCC=C(C)CO
  pmozprop( 80) = mozsp( 'LISOPOOHOOH', (/2.e16_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2,3,4-butanetetrol as proxy SMILES: OC(C)(COO)C(CO)OO
  pmozprop( 81) = mozsp( 'LISOPNO3OOH', (/3.e11_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2,3-butanetriol as proxy SMILES: O=[N+]([O-])OC(C)(CO)C(CO)OO 
  pmozprop( 82) = mozsp( 'LISOPNO3NO3', (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: O=[N+]([O-])OC(C)(CO)C(O[N+]([O-])=O)CO
  pmozprop( 83) = mozsp( 'LC578OOH'   , (/3.e11_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2,3-butanetriol as proxy SMILES: OCC(O)C(C)(OO)C=O
  pmozprop( 84) = mozsp( 'C59OOH'     , (/3.e11_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2,3-butanetriol as proxy SMILES: OCC(=O)C(C)(CO)OO
  pmozprop( 85) = mozsp( 'BIGALD2'    , (/2.5e5_dp    ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : TOL_EPOX as proxy MCM name C5DICARB SMILES: O=CC=CC(=O)C 
  pmozprop( 86) = mozsp( 'BIGALD3'    , (/2.5e5_dp    ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : TOL_EPOX as proxy MCM name C4MDIAL SMILES:  O=CC=C(C)C=O  
  !                6C
  pmozprop( 87) = mozsp( 'BZALD'      , (/38._dp      , 5500._dp/)     , 0.1_dp  )  ! Sander (2015) SMILES:  O=Cc1ccccc1 
  pmozprop( 88) = mozsp( 'BZOOH'      , (/2.9e3_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : benzyl alcohol as proxy SMILES: OOCc1ccccc1
  pmozprop( 89) = mozsp( 'PHENOOH'    , (/3.e11_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2,3-butanetriol as proxy SMILES: OOC1(O)C=CC2OOC1C2O
  pmozprop( 90) = mozsp( 'C6H5OOH'    , (/3.e3_dp     , 5900._dp/)     , 1.0_dp  )  ! Sander (2015) : phenol as proxy see Harrison et al. (2001) SMILES: OOc1ccccc1
  pmozprop( 91) = mozsp( 'PHENOL'     , (/3.e3_dp     , 5900._dp/)     , 0.1_dp  )  ! Sander (2015) : Harrison et al. (2001) SMILES: Oc1ccccc1
  pmozprop( 92) = mozsp( 'CATECHOL'   , (/4.6e3_dp    ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : 1,2-dihydroxybenzene SMILES: Oc1ccccc1O
  pmozprop( 93) = mozsp( 'CATEC1OOH'  , (/4.6e3_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : CATECHOL as proxy SMILES: OOc1ccccc1O
  pmozprop( 94) = mozsp( 'BEPOMUC'    , (/2.5e5_dp    ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : TOL_EPOX MCM name BZEPOXMUC SMILES: O=CC=CC1OC1C=O
  pmozprop( 95) = mozsp( 'BENZOOH'    , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy MCM name BZBIPEROOH SMILES: OOC1C=CC2OOC1C2O 
  pmozprop( 96) = mozsp( 'PBZNIT'     , (/2.8_dp      , 5730._dp/)     , 0.1_dp  )  ! PAN as proxy MCM name PBZN SMILES: O=N(=O)OOC(=O)c1ccccc1 
  pmozprop( 97) = mozsp( 'BIGALD4'    , (/2.5e5_dp    ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : TOL_EPOX as proxy, according to G. Tyndall MCM name is C5MDICARB SMILES: O=CC(=CC(=O)C)C
  !                7C
  pmozprop( 98) = mozsp( 'TOLOOH'     , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy MCM name TLBIPEROOH SMILES: OOC1C=CC2(C)OOC1C2O  
  pmozprop( 99) = mozsp( 'CRESOL'     , (/1.1e3_dp    , 6700._dp/)     , 0.1_dp  )  ! Sander (2015) : Harrison et al. (2001) SMILES: Cc1ccccc1O
  pmozprop(100) = mozsp( 'TEPOMUC'    , (/2.5e5_dp    ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : TOL_EPOX MCM name TLEPOXMUC SMILES: O=CC1OC1C=CC(=O)C 
  pmozprop(108) = mozsp( 'TERPROD2'   , (/2.5e5_dp    ,    0._dp/)     , 0.1_dp  )  ! Sander (2015) : TOL_EPOX as proxy, oxygenated VOC with 7C (G. Tyndall pers. comm.) assumed MCM name SMILES: ??? 
  !               8C
  pmozprop(101) = mozsp( 'XYLENOOH'   , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy MCM name OXYBPEROOH SMILES: OOC1(C)C=CC2OOC1(C)C2O 
  pmozprop(102) = mozsp( 'XYLOL'      , (/1.1e3_dp    , 6700._dp/)     , 0.1_dp  )  ! cresol as proxy MCM name OXYLOL SMILES: Cc1c(C)cccc1O 
  pmozprop(103) = mozsp( 'XYLOLOOH'   , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy MCM name OXYOLOOH SMILES: OOC1(O)C=C(C)C2(C)OOC1C2O   
  !                9C
  !                10C+more
  pmozprop(104) = mozsp( 'TERPOOH'    , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy SMILES: OOC1(C)C(O)CC2CC1C2(C)C 
  pmozprop(105) = mozsp( 'TERP2OOH'   , (/2.1e5_dp    , 9900._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2-butanediol as proxy MCM name  SMILES: ???
  pmozprop(106) = mozsp( 'TERPNO3'    , (/5.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! like an isoprene nitrate MCM name APINANO3 SMILES:  O=N(=O)OC1(C)C(O)CC2CC1C2(C)C 
  pmozprop(107) = mozsp( 'TERPROD1'   , (/2.5e5_dp    ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : TOL_EPOX as proxy MCM name C109OOH second major APINENE ozonolysis comp. under low NOx SMILES: OOCC(=O)C1CC(CC=O)C1(C)C 
  pmozprop(109) = mozsp( 'NTERPNO3'   , (/5.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! like an isoprene nitrate MCM name NAPINAOOH SMILES: CC1(C)C2CC1C(C)(OO)C(O[N](=O)=O)C2
  pmozprop(110) = mozsp( 'MTHOM'     , (/2.e9_dp     ,    0._dp/)     , 1.0_dp  )  ! Sander (2015) : 1,2,6-hexanetriol as proxy SMILES:  CC(=O)C(OO)C1CC(OO)(C(=O)OO)C1(C)C 
  !                xC
  pmozprop(111) = mozsp( 'ALKOH'      , (/336._dp     , 5995._dp/)     , 0.1_dp  )  ! ** as C2H5OOH MCM name IPECOH SMILES: CCC(C)(C)O  
  pmozprop(112) = mozsp( 'ALKOOH'     , (/336._dp     , 5995._dp/)     , 1.0_dp  )  ! ** as C2H5OOH MCM name IPECOOH SMILES: CCC(C)(C)OO 
  pmozprop(113) = mozsp( 'POOH'       , (/300._dp     , 5280._dp/)     , 0.1_dp  )  ! ** as CH3OOH MCM name HYPROPO2H SMILES: OCC(C)OO 
  pmozprop(114) = mozsp( 'ROOH'       , (/300._dp     , 5280._dp/)     , 1.0_dp  )  ! ** as CH3OOH MCM name HYPERACT SMILES: CC(=O)COO  
  pmozprop(115) = mozsp( 'EOOH'       , (/300._dp     , 5280._dp/)     , 0.1_dp  )  ! ** as CH3OOH MCM name HYETHO2H SMILES: OCCOO 
  pmozprop(116) = mozsp( 'ALKNO3'     , (/5.e3_dp     ,    0._dp/)     , 0.1_dp  )  ! like for ISOPNO3 from Nguyen et al. 2015 MCM name IPECNO3 SMILES: CCC(C)(C)ON(=O)=O  
  ! section : inorganic Halogens
  pmozprop(117) = mozsp( 'HF'         , (/1.3e4_dp    ,    0._dp/)     , 0._dp   )  ! Sander (2015) 
  pmozprop(118) = mozsp( 'HCL'        , (/2.0e11_dp   ,  600._dp/)     , 0._dp   )  ! Fernandez et al. 2014
  pmozprop(119) = mozsp( 'HOCL'       , (/6.6e2_dp    , 5900._dp/)     , 0._dp   )  ! R. Sander (1999) 
  pmozprop(120) = mozsp( 'CL2'        , (/9.2e-2_dp   , 2000._dp/)     , 0._dp   )  ! R. Sander (1999) 
  pmozprop(121) = mozsp( 'CLONO2'     , (/2.e13_dp    ,    0._dp/)     , 1._dp   )  ! R. Sander (1999) : listed as "infinite"
  pmozprop(122) = mozsp( 'BR2'        , (/7.2e-1_dp   ,    0._dp/)     , 0._dp   )  ! Sander (2015) 
  pmozprop(123) = mozsp( 'HBR'        , (/7.2e13_dp   , 6100._dp/)     , 0._dp   )  ! Fernandez et al. 2014
  pmozprop(124) = mozsp( 'HOBR'       , (/6.1e3_dp    ,    0._dp/)     , 0._dp   )  ! R. Sander (1999) 
  pmozprop(125) = mozsp( 'BRONO2'     , (/2.e13_dp    ,    0._dp/)     , 1._dp   )  ! R. Sander (1999) : listed as "infinite"
  pmozprop(126) = mozsp( 'BRONO'      , (/3.e-1_dp    ,    0._dp/)     , 0._dp   )  ! Sander (2015) 
  pmozprop(127) = mozsp( 'BRNO2'      , (/3.e-1_dp    ,    0._dp/)     , 0._dp   )  ! Sander (2015) 
  pmozprop(128) = mozsp( 'BRCL'       , (/9.7e-1_dp   , 5600._dp/)     , 0._dp   )  ! Sander (2015) 
!  ! section : Sulfur
  pmozprop(129) = mozsp( 'DMS'        , (/0.54_dp     , 3460._dp/)     , 0._dp   )  ! JPL17(2011)
  pmozprop(130) = mozsp( 'SO2'        , (/2.45E5_dp   , 3100._dp/)     , 0._dp   )  ! Sander (2015): H(298)=1.2 M/atm and H*(298) calculated at pH=7 with K1a=1.23E-2 M and K2a=6.61E-8 M 
  pmozprop(131) = mozsp( 'H2SO4'      , (/1.3e15_dp   , 20000._dp/)    , 0._dp   )  ! Sander (2015)
  pmozprop(132) = mozsp( 'DMSO'       , (/9.8e4_dp    ,     0._dp/)    , 0.1_dp  )  ! JPL17(2011)
  pmozprop(133) = mozsp( 'CH3SO3H'    , (/1.e30_dp    ,     0._dp/)    , 1._dp   )  ! As in EMAC (ECHAM/MESSy)
  ! section NH3
  pmozprop(134) = mozsp( 'NH3'        , (/1.02e4_dp   , 4200._dp/)     , 0._dp   )  ! Sander (2015): H(298)=58 M/atm and H*(298) calculated at pH=7 with Kb=1.75E5 M

!--mgs (NEW LIST)
! Don't forget to update to last number in the first statement of this
! subroutine if you change the table!

!!Attention! See comment above!!!
  kntrx = ntrx

END SUBROUTINE moz_species_properties


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! moz_species: define species based on MOZART preprocessor output (mo_moz_subs)
!!
!! @author see module info
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! moz_init
!!

!! @par Notes
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE moz_species (lcouple_oxidants)

  USE mo_kind,          ONLY: dp
  USE mo_util_string,   ONLY: toupper
  USE mo_string_utls,   ONLY: st1_in_st2_idx
  USE mo_moz_mods,      ONLY: pcnstm1, adv_mass
  USE mo_moz,           ONLY: tracnam
  USE mo_tracdef,       ONLY: GAS, itrprog, ln, jptrac
  USE mo_species,       ONLY: new_species, query_species, speclist
  USE mo_exception,     ONLY: message, em_debug, em_warn, em_info
 
  LOGICAL, INTENT(in)     :: lcouple_oxidants   ! re-define existing species as prognostic
                                                ! or rename to '_clim'
  INTEGER            :: i, idd, idx, ispec, ntrx
  CHARACTER(len=24)  :: ztracnam(jptrac)
  CHARACTER(len=64)  :: clongname, cshortname
  REAL(dp)           :: zeps
  TYPE(mozsp)        :: pmozprop(jptrac)

  zeps=EPSILON(1.0_dp)

  CALL message('moz_species', 'starting...', level=em_debug)

  ! define henry and f0 values
  CALL moz_species_properties (ntrx, pmozprop)
  DO i = 1, ntrx
    ztracnam(i) = pmozprop(i)%name
  END DO

  ! loop over MOZ tracnam and define new species or re-use existing ones
  DO i = 1, pcnstm1
    ! -- special handling of H2SO4: Rename to SO4 as in HAM!
!! ### check again! MOZ has SO4 and H2SO4 defined!!
!!    IF (TRIM(tracnam(i)) == 'H2SO4') tracnam(i) = 'SO4'

    clongname  = tracnam(i)
    cshortname = tracnam(i)

    ! check if species undergoes dry and wet deposition
    idd = st1_in_st2_idx( tracnam(i), ztracnam )

    ! test if species has already been defined
    CALL query_species(shortname=tracnam(i), nphase=GAS, index=ispec)
    IF (ispec > 0) THEN
      ! re-define existing species as MOZ species if oxidants shall be coupled
      IF (lcouple_oxidants) THEN

        IF (TRIM(toupper(speclist(ispec)%tsubmname)) /= 'HAM') THEN
          CALL message('moz_species', 'Re-defining species '//TRIM(tracnam(i))//   &
                       ' although it is not from HAM!', level=em_warn)
        ELSE
          CALL message('moz_species', 'Re-defining species '//TRIM(tracnam(i))//   &
                       ' initially defined by HAM!', level=em_info)
        END IF
        speclist(ispec)%tsubmname = 'MOZ'
        speclist(ispec)%itrtype = itrprog
        speclist(ispec)%moleweight = adv_mass(i)
        IF (idd > 0) THEN
          speclist(ispec)%henry      = pmozprop(idd)%henry
          speclist(ispec)%dryreac    = pmozprop(idd)%f0
          speclist(ispec)%ldrydep    = .true.
!++mgs 2013/12/12: set wetdep true only if henry(1) >= 1.
          speclist(ispec)%lwetdep    = (pmozprop(idd)%henry(1) >= 1.0_dp )
!--mgs
        ELSE
          IF (speclist(ispec)%ldrydep) THEN
            CALL message('moz_species', 'Setting ldrydep to FALSE for species '//TRIM(tracnam(i)), &
                         level=em_info)
            speclist(ispec)%ldrydep   = .false.
          END IF
          IF (speclist(ispec)%lwetdep) THEN
            CALL message('moz_species', 'Setting lwetdep to FALSE for species '//TRIM(tracnam(i)), &
                         level=em_info)
            speclist(ispec)%lwetdep   = .false.
          END IF
        END IF

      ELSE      ! lcouple_oxidants

!++sschr: #283 Code inconsistent between lhammoz=true and lhammoz=false (species and tracer definitions) 
!for the lcouple_oxidants=false case: if a species is already defined, add a "_moz" to the species name [unchanged]

        clongname = TRIM(tracnam(i))//'_moz'
        cshortname = TRIM(tracnam(i))//'_moz'
        ispec = 0                               ! force definition of new species

!--sschr: #283

      END IF    ! lcouple_oxidants

    END IF

    ! define new species if no valid species is already in the list
    IF (ispec <= 0) THEN
      IF (idd > 0) THEN
        CALL new_species(GAS,                               &
                         clongname,                         &
                         cshortname,                        &
                         'mol mol-1',                       &
                         adv_mass(i),                       &
                         tsubmname = 'MOZ',                 &
                         itrtype   = itrprog,               &
                         henry     = pmozprop(idd)%henry,   &
                         dryreac   = pmozprop(idd)%f0,      &
                         ldrydep   = .true.,                &
!++mgs 2013/12/12: set wetdep true only if henry(1) >= 1.
                         lwetdep   = (pmozprop(idd)%henry(1) >= 1.0_dp ),  &
!--mgs
                         idx=idx)
      ELSE
        CALL new_species(GAS,                               &
                         clongname,                         &
                         cshortname,                        &
                         'mol mol-1',                       &
                         adv_mass(i),                       &
                         tsubmname = 'MOZ',                 &
                         itrtype   = itrprog,               &
                         henry     = (/ 0._dp, 0._dp /),    &
                         dryreac   = 0._dp,                 &
                         ldrydep   = .false.,               &
                         lwetdep   = .false.,               &
                         idx=idx)
      END IF
    END IF
  END DO

!### may be necessary to "clean up" the definition of drydep and wetdep species !!!

  CALL message('moz_species', 'done.', level=em_debug)

END SUBROUTINE moz_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! helper routine: expand_default 
! replace 'default' token with predefined list of names

SUBROUTINE expand_default(speclist, deflist, ndef)

  USE mo_exception,                ONLY: finish, message, message_text, em_debug
  USE mo_mpi,                      ONLY: p_parallel_io


  CHARACTER (LEN=32), INTENT(inout)     :: speclist(:)
  CHARACTER (LEN=32), INTENT(in)        :: deflist(:)
  INTEGER, INTENT(in)                  :: ndef

  INTEGER                  :: i, ii, nspec, nspeclist

  WRITE(message_text,'(a,2i0)') "size,shape=", size(speclist), shape(speclist)
  CALL message('expand_default', message_text, level=em_debug)
 
  nspeclist = SIZE(speclist)
  ! find number of species in speclist. We already know that 'default' is on position 1
  nspec = 0
  DO i=2,nspeclist
     IF (TRIM(speclist(i)) /= '') nspec = nspec + 1
  END DO
  ! shift species after 'default' to make room for default species
  IF (nspec+ndef > nspeclist) THEN         ! shouldn't happen...
     CALL finish('expand_default', 'Too many species in out_species or burden_species list!')
  END IF
  ! fill values from the end of list (don't overwrite entries on positions after ndef!)
  DO i=nspec,1,-1
     speclist(ndef+i) = speclist(i+1)
  END DO
  nspec = nspec+ndef
  ! insert default species
  DO i=1,ndef
     speclist(i) = deflist(i)
  END DO
  ! remove duplicates
  DO i=nspec,1,-1
     DO ii=1,i-1
        IF (TRIM(speclist(i)) == TRIM(speclist(ii))) THEN
           speclist(i:nspec-1) = speclist(i+1:nspec) 
           speclist(nspec) = ''
        END IF
     END DO
  END DO
  DO i=1,nspec
     WRITE(message_text,'(i0,a,a)') i,': ',TRIM(speclist(i))
     CALL message('expand_default', message_text, level=em_debug)
  END DO

END SUBROUTINE expand_default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! moz_define_tracer: create ECHAM tracers based on MOZ species
!!
!! @author see module info
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Externals:
!! <ol>
!! <li>none
!! </ol>
!!
!! @par Notes
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE moz_define_tracer

  USE mo_tracdef,           ONLY: OK, RESTART, &
                                  CONSTANT, INITIAL, &
                                  PERIODIC, ON, OFF
  USE mo_tracer,            ONLY: new_tracer, get_tracer, new_diag_burden
  USE mo_species,           ONLY: nspec, speclist
! USE mo_submodels,         ONLY: lhmzhet     ! heterogeneous chemistry coupling in HAMMOZ
  USE mo_moz_mods,          ONLY: pcnstm1
  USE mo_moz,               ONLY: tracnam, budget_species, out_species, burden_species, &
                                  ndrydep, nwetdep, e2m
  USE mo_moz_util,          ONLY: get_spc_ndx
  USE mo_util_string,       ONLY: tolower
  USE mo_string_utls,       ONLY: st1_in_st2_proof
  USE mo_exception,         ONLY: finish, message, message_text, em_error, em_warn
 
  IMPLICIT NONE

  INTEGER, PARAMETER            :: ndef_out=16
  CHARACTER(LEN=32)             :: cdef_out(1:ndef_out)= &
                                (/'C2H6    ', &
                                  'C3H8    ', &
                                  'C5H8    ', &
                                  'CH2O    ', &
                                  'CH3OOH  ', &
                                  'CO      ', &
                                  'H2O2    ', &
                                  'HNO3    ', &
                                  'HO2     ', &
                                  'HO2NO2  ', &
                                  'N2O5    ', &
                                  'NO      ', &
                                  'NO2     ', &
                                  'NO3     ', &
                                  'O3      ', &
                                  'PAN     '    /)
  INTEGER, PARAMETER            :: ndef_burden=4
  CHARACTER(LEN=32)             :: cdef_burden(1:ndef_burden)= &
                                (/'O3      ', &
                                  'NO      ', &
                                  'NO2     ', &
                                  'CO      '    /)
 
  INTEGER :: ierr, i, idt, init, iburdenid, ihelp, jt
  INTEGER :: nbud, nddep, nwdep, nwritet, ispc_ndx
  CHARACTER(LEN=32) :: helpname, mapname
 !                                                  ! nwdep: 2 for in-cloud scavenging only
 !                                                  !        3 for in-cloud and precipitation scavenging

  ! -- expand 'default' and 'all' strings in out_species list
  IF (tolower(TRIM(out_species(1))) == 'default') THEN
     CALL expand_default(out_species, cdef_out, ndef_out)
  END IF 
  IF (tolower(TRIM(out_species(1))) == 'all') THEN
    DO jt= 1, pcnstm1 ! only MOZART tracers
      out_species(jt) = TRIM(tracnam(jt))
    END DO
  END IF 
  ! -- expand 'default' and 'all' strings in burden_species list
  IF (tolower(TRIM(burden_species(1))) == 'default') THEN
     CALL expand_default(burden_species, cdef_burden, ndef_burden)
  END IF 
  IF (tolower(TRIM(burden_species(1))) == 'all') THEN
    DO jt= 1, pcnstm1 ! only MOZART tracers
      burden_species(jt) = TRIM(tracnam(jt))
    END DO
  END IF 
  ! -- check validity of tracer lists from mozctl namelist 
  IF (.NOT. st1_in_st2_proof(out_species, tracnam, ierr=ierr)) THEN
     IF (ierr > 0) &
     CALL message ('moz_define_tracer',' species '//TRIM(out_species(ierr))// &
                   ' in out_species does not exist in species list', level=em_warn)
  END IF
  IF (.NOT. st1_in_st2_proof(budget_species, tracnam, ierr=ierr)) THEN
     IF (ierr > 0) &
     CALL message ('moz_define_tracer',' species '//TRIM(budget_species(ierr))// &
                  ' in budget_species does not exist in species list', level=em_warn)
  END IF
  IF (.NOT. st1_in_st2_proof(burden_species, tracnam, ierr=ierr)) THEN
     IF (ierr > 0) &
     CALL message ('moz_define_tracer',' species '//TRIM(burden_species(ierr))// &
                   ' in burden_species does not exist in species list', level=em_warn)
  END IF

  ! -- preset flag values
  ! Note: File(s) containing initial values for tracers must be in MMR!!!
  init=RESTART+CONSTANT+INITIAL

  ! -- define tracers based on species list.
  ! +++ note: perhaps we can generalize this at some point. Now define only MOZ tracers in this way.

  DO i = 1,nspec
    ! look for MOZ species
    IF (TRIM(speclist(i)%tsubmname) /= 'MOZ') CYCLE

    nbud    = OFF
    nddep   = OFF
    nwdep   = OFF
    nwritet = OFF
    iburdenid = 0
    IF(st1_in_st2_proof(speclist(i)%shortname,out_species)) nwritet = ON
    IF(st1_in_st2_proof(speclist(i)%shortname,budget_species)) nbud = 1
    IF(st1_in_st2_proof(speclist(i)%shortname,burden_species)) &
      ! define burden diagnostics with total column mass and tropospheric mass
      iburdenid = new_diag_burden(speclist(i)%shortname, itype=3, lclobber=.false.)

!++sschr: #283 Code inconsistent between lhammoz=true and lhammoz=false (species and tracer definitions) 
! If get_tracer(stripped_named) returns idt/=0, then define the tracer with the "_moz",
!                                               otherwise use the stripped_name as tracer name.

    helpname = speclist(i)%shortname
    ihelp=index(speclist(i)%shortname,'_moz')
    IF (ihelp > 0) helpname = speclist(i)%shortname(:ihelp-1)
    mapname = helpname
  ! Don't use argument ierr=... ==> finish is called if tracers are not found
    CALL get_tracer(helpname, ierr=ierr)
    IF (ierr .eq. 0) helpname = speclist(i)%shortname
      
!--sschr: #283 

    IF (ihelp > 0) THEN
      IF(st1_in_st2_proof(speclist(i)%shortname(1:ihelp-1),out_species)) nwritet = ON
      IF(st1_in_st2_proof(speclist(i)%shortname(1:ihelp-1),budget_species)) nbud = 1
      IF(st1_in_st2_proof(speclist(i)%shortname(1:ihelp-1),burden_species)) &
        ! define burden diagnostics with total column mass and tropospheric mass
        iburdenid = new_diag_burden(helpname, itype=3, lclobber=.false.)
    ENDIF

    IF (speclist(i)%ldrydep) nddep = ndrydep
    IF (speclist(i)%lwetdep) nwdep = nwetdep
        ! note: only two gas-phase species should be subject to below-cloud scavenging: HNO3 and H2SO4


    CALL new_tracer(helpname, 'MOZ', units='kg kg-1',           &
                    moleweight = speclist(i)%moleweight,                     &
                    nphase     = speclist(i)%nphase,                         &
                    spid       = i,                                          &
                    nbudget    = nbud,                                       &
!                   nsemis     = nemis,                                      &
                    ndrydep    = nddep,                                      &
                    nwetdep    = nwdep,                                      &
                    ninit      = init,                                       &
                    nwrite     = nwritet,                                    &
                    code       = i,                                          &
                    burdenid   = iburdenid,                                  &
                    idx        = idt,                                        & 
                    ierr       = ierr       )

    IF (ierr == OK) THEN
      speclist(i)%idt = idt
      ! get_spc_ndx: get MOZ species number (from mozpp)
      ! this may create duplicate entries in e2m if lhammoz=.false.
      ! this doesn't matter, because X_moz is always defined after X
      ispc_ndx = get_spc_ndx(mapname)
      e2m(idt) = ispc_ndx
!     WRITE(message_text,'(3a,2i4)')"shortname,mapname,idt,ispc_ndx=",speclist(i)%shortname,mapname,idt,ispc_ndx
!     CALL message('moz_define_tracer', message_text)
    ELSE
      WRITE(message_text,'(a,i0)') 'new_tracer returned error code ',ierr
      CALL message('moz_define_tracer', message_text, level=em_error)
    END IF

  END DO
 
END SUBROUTINE moz_define_tracer


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! moz_initialize: perform remaining initialisation after tracers are defined
!!
!! @author see module info
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! moz_init
!!
!! @par Notes
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE moz_initialize

  USE mo_exception,       ONLY : finish, message, message_text, em_param, em_debug
  USE mo_control,         ONLY : nproma, nlev
  USE mo_tracer,          ONLY : get_tracer
!++mgs #249: added lchemsolv
  USE mo_moz,             ONLY : lchemsolv, lphotolysis, lfastj, lhammoniagases, &
                                 lbc_species, ubc_species, nlbc, nubc, bc_lbc, bc_ubc
!--mgs
  USE mo_moz_mods,        ONLY : plonl, plev, plevp, plnplv
  USE mo_moz_util,        ONLY : get_spc_ndx
  USE mo_moz_setinv,      ONLY : invariants_inti
  USE mo_moz_waccm_photo, ONLY : waccm_prate_inti
  USE mo_moz_lbc_ubc,     ONLY : moz_lbc_ubc_initialize
  USE mo_moz_chemdr,      ONLY : chemdr_inti
!++sschr 2014/09/03: added lchemistry
  USE mo_submodel,        ONLY : lchemistry, linterh2o
!--sschr 2014/09/03: added lchemistry


  INTEGER :: idt_h2o, ierr

  !-----------------------------------------------------------------------
  !     ... Initialize species and reaction indices for MOZART solver and diagnostics
  !-----------------------------------------------------------------------
  IF ( get_spc_ndx('O2') > 0 .AND.    &
       get_spc_ndx('CO2') > 0 .AND.   &
       get_spc_ndx('N2O') > 0 .AND.   &
       get_spc_ndx('CH4') > 0 ) THEN
    lhammoniagases = .TRUE.
  ENDIF
!++mgs #249: added if condition and moved call to invariants_inti
  IF (lchemistry) THEN
    CALL invariants_inti
    CALL message('moz_initialize','finished invariants_inti', level=em_debug)
  ELSE
    CALL message('moz_initialize','lchemistry=.false.', level=em_debug)
  END IF
!--mgs

  !-- evaluate namelist string lists for lower and upper boundary conditions and define
  !   boundary conditions
  CALL moz_lbc_ubc_initialize(lbc_species, ubc_species, nlbc, nubc, bc_lbc, bc_ubc)

  !-----------------------------------------------------------------------
  !       ... Initialize photorate module
  !-----------------------------------------------------------------------
  IF ( lphotolysis ) THEN
    CALL waccm_prate_inti
    CALL message('moz_initialize','finished waccm_prate_inti', level=em_debug)
  END IF

  !-----------------------------------------------------------------------
  !       ... Initialize chemistry solvers and stratospheric aerosol
  !           and identify moz species for diagnostics
  !-----------------------------------------------------------------------
!++mgs #249: added if condition
  IF (lchemistry) THEN
    CALL chemdr_inti
  END IF
!--mgs

  !-----------------------------------------------------------------------
  !       .. Check if water vapour tracer is defined when linterh2o = true
  !-----------------------------------------------------------------------
  CALL get_tracer( 'H2O', idx=idt_h2o, ierr=ierr )
  IF ( linterh2o ) THEN
     IF ( ierr == 0 ) THEN
        write(message_text,'(a,i0)') 'linterh2o = true. Tracer id for chemical water vapour = ',idt_h2o
        CALL message('moz_initialize', message_text, level=em_param)
     ELSE
        CALL finish('moz_initialize','linterh2o true, but no chemical water vapour tracer defined.')
     END IF
  END IF

  CALL message('moz_initialize','Done.', level=em_debug)

END SUBROUTINE moz_initialize

!---------------------------------------------------------------------------
SUBROUTINE init_moz_ic

! Pre-set initial conditions for MOZART. These will be overwritten by an 
! initial conditions file if the respective species are contained in the file.
! ### Ideally, the ics defined here should suffice to initialize a "technically correct"
! simulation. For now, we only focus on O2.

  USE mo_kind,                 ONLY: dp
  USE mo_exception,            ONLY: message, em_error
  USE mo_radiation_parameters, ONLY: o2vmr              ! O2 volume mixing ratio
! USE mo_moz_chemdr,           ONLY: lhas_o2
  USE mo_tracer,               ONLY: get_tracer
  USE mo_tracdef,              ONLY: trlist,         &  ! tracer info variable
                                     CONSTANT


  INTEGER          :: idt          ! tracer index

  ! --- initialize O2 mixing ratio
    CALL get_tracer('O2', idx=idt)
    IF (idt > 0) THEN
      trlist% ti(idt)% ninit = IOR(trlist% ti(idt)% ninit, CONSTANT)
      trlist% ti(idt)% vini  = o2vmr * 31.9988003_dp / 29.87_dp
    ELSE
      CALL message('init_moz_ic', 'No O2 tracer defined!', &
                   level=em_error)
    END IF
END SUBROUTINE init_moz_ic

END MODULE mo_moz_init
