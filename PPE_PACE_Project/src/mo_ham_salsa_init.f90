!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_salsa_init.f90
!!
!! \brief
!! Contains subroutines and functions that are used to initialize the particle size grid and aerosol 
!! processes for SALSA.         
!!
!! \author Harri Kokkola (FMI)
!!
!! \responsible_coder
!! Harri Kokkola, harri.kokkola@fmi.fi
!!
!! \revision_history
!!   -# H. Korhonen (FMI) - original code (2005)
!!   -# H. Kokkola (FMI) - original code (2006-2014)
!!   -# M. Niskanen(FMI) - change from 3 subregions to 2 subregions (2012)
!!   -# A. Laakso (FMI) - ECHAM6-HAMMOZ implementation with some changes (2013)
!!
!! \limitations
!! None
!!
!! \details
!! This module defines salsa sizebins and wetdep parameters. Salsa_initialize is called from
!! mo_ham_init when SALSA is used. Salsa_initialize calls set_sizebins.
!! This module contains the following subroutines
!!      set_sizebins
!!      set_coagc
!!      actcurve
!!      salsa_initialize
!! and function
!!      coagc
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_ham_salsa_init

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: coagc
  PUBLIC :: salsa_initialize

CONTAINS

  ! fxm: when dpmid is used for calculating coagulation coefficients
  ! (only?), would it make more sense to use approximate wet radii
  ! e.g. for sea salt particles?
  !********************************************************************
  !
  ! subroutine SET_SIZEBINS()
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Initializes particle size distribution grid by 
  ! calculating size bin limits and mid-size for
  ! *dry* particles in each bin
  !
  !
  ! Method:
  ! ------- 
  ! Size distribution described using
  !   1) moving center method (regimes 1 and 2)
  !   (Jacobson, Atmos. Env., 31, 131-144, 1997)
  !   2) fixed sectional method (regime 3)
  ! 
  ! Size bins in each regime are spaced logarithmically
  ! based on given regime size limits and bin number.
  !
  !
  ! Interface:
  ! ----------
  ! Called from model driver
  ! (only at the beginning of simulation)
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005 
  ! Harri Kokkola (FMI) 2006
  ! Matti Niskanen(FMI) 2012
  ! Anton Laakso  (FMI) 2013
  !
  !---------------------------------------------------------------------

  SUBROUTINE set_sizebins()

    USE mo_kind, ONLY : dp

    USE mo_math_constants, ONLY : &
         pi_6

    USE mo_ham_salsactl, ONLY : &
         reglim,      & ! diameter limits for size regimes [m]
         nbin,        & ! number of size bins in each (sub)regime
         nbin2,       & ! number of bins in former 2-region
         nbin3,       & !     - " -       3-region
         in1a, fn1a,  & ! size regime bin indices: 1a
         in2a, fn2a,  & !     - " -       2a
         in2b, fn2b,  & !     - " -       2b
         vlolim,      & ! low volume limit for bins in regimes 1 and 2 [fxm]
         vhilim,      & ! high    - " -
         vratiohi,    & ! volume ratio from dpmid to bounds
         vratiolo,    & ! volume ratio from dpmid to bounds
         dpmid          ! mid-diameter (dry) for bins in all regimes [m]  

    IMPLICIT NONE

    !-- local variables ----------
    INTEGER :: ii, jj ! loop indices

    REAL(dp) :: ratio ! ratio of regime upper and lower diameter

 
    !-- 1) size regime 1: --------------------------------------
    !  - minimum & maximum *dry* volumes [fxm] 
    !  - bin mid *dry* diameter [m]

    ratio = reglim(2)/reglim(1)   ! section spacing

    DO ii = in1a,fn1a
       vlolim(ii) = pi_6*(reglim(1)*ratio**(real(ii-1)/nbin(1)))**3
       vhilim(ii) = pi_6*(reglim(1)*ratio**(real(ii)/nbin(1)))**3
       dpmid(ii) = ((vhilim(ii) + vlolim(ii))/(2.*pi_6))**(1._dp/3._dp)
       vratiohi(ii)= (vhilim(ii))/(pi_6*dpmid(ii)**3)
       vratiolo(ii)= (vlolim(ii))/(pi_6*dpmid(ii)**3)

    END DO


    !-- 2) size regime 2: --------------------------------------
    !  - minimum & maximum *dry* volumes [fxm] 
    !  - bin mid *dry* diameter [m]

    !-- 2.1.1) subregime 2a beginning
    ratio = reglim(3)/reglim(2)   ! section spacing

    DO jj = in2a,in2a+(nbin2-1)
       ii = jj - in2a

       vlolim(jj) = pi_6*(reglim(2)*ratio**(real(ii)/real(nbin(2)-nbin3)))**3
       vhilim(jj) = pi_6*(reglim(2)*ratio**(real(ii+1)/real(nbin(2)-nbin3)))**3
       dpmid(jj) = ((vhilim(jj) + vlolim(jj))/(2._dp*pi_6))**(1._dp/3._dp)
       vratiohi(jj)= (vhilim(jj))/(pi_6*dpmid(jj)**3)
       vratiolo(jj)= (vlolim(jj))/(pi_6*dpmid(jj)**3)

    END DO

    !-- 2.1.2) subregime 2a end
    ratio = reglim(4)/reglim(3)   ! section spacing

    DO jj = fn2a-(nbin3-1),fn2a
       ii = jj - (fn2a-(nbin3-1))

       vlolim(jj) = pi_6*(reglim(3)*ratio**(real(ii)/real(nbin(2)-nbin2)))**3
       vhilim(jj) = pi_6*(reglim(3)*ratio**(real(ii+1)/real(nbin(2)-nbin2)))**3
       dpmid(jj) = ((vhilim(jj) + vlolim(jj))/(2._dp*pi_6))**(1._dp/3._dp)
       vratiohi(jj)= (vhilim(jj))/(pi_6*dpmid(jj)**3)
       vratiolo(jj)= (vlolim(jj))/(pi_6*dpmid(jj)**3)

    END DO

    !-- 2.2) same values for subregime 2b
    vlolim(in2b:fn2b)   = vlolim(in2a:fn2a)
    vhilim(in2b:fn2b)   = vhilim(in2a:fn2a)
    dpmid(in2b:fn2b)    = dpmid(in2a:fn2a)
    vratiohi(in2b:fn2b) = vratiohi(in2a:fn2a)
    vratiolo(in2b:fn2b) = vratiolo(in2a:fn2a)

  END SUBROUTINE set_sizebins


  ! fxm: Here we use a same approximate density for all bins; should be changed?
  ! fxm: Should we use a larger bin than in3a to represent regime 3
  !********************************************************************
  !
  ! subroutine SET_COAGC(klev,pres0,temp0)
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates particle coagulation coefficients for
  ! each bin pair (regime 3 particles lumped together)
  ! and each pressure level
  !
  !
  ! Method:
  ! ------- 
  ! The coefficients are calculated for particle diameters 
  ! *dpmid*; the coefficients are scaled at every timestep to 
  ! account for the actual wet size of the particles.
  !
  ! The coefficients are calculated separately for 
  ! each pressure level because atmospheric pressure 
  ! has a significant nonlinear effect on the
  ! coefficients. The effect of temperature
  ! is much weaker and 'typical' temperatures for each
  ! pressure level can be used.
  !
  !
  ! Interface:
  ! ----------
  ! Called from model driver
  ! (only at the beginning of simulation)
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005 
  !
  !---------------------------------------------------------------------

  SUBROUTINE set_coagc(klev,pres0,temp0)

    USE mo_math_constants, ONLY : pi_6

    USE mo_kind, ONLY : dp

    USE mo_ham_salsactl, ONLY : &
         dpmid,    & ! mid-diameter (dry) of each size bin [m]
         fn2b

    IMPLICIT NONE

    !-- input variables --------
    INTEGER, INTENT(in) :: &
         klev    ! number of vertical levels

    REAL(dp), INTENT(in) ::  &
         pres0(klev), & ! pressure at each vertical level [Pa]
         temp0(klev)   ! temperature at each vertical level [K]

    !-- local variables --------
    INTEGER :: ii, jj, kk  ! loop indices

    REAL(dp) :: &
         mpart(fn2b)   ! approximate mass of particles [kg]

    !------------------------------------------------------------------------------- 

    ! fxm: should we use something else than fixed density?
    !-- particle mass; density of 1500 kg/m3 assumed [kg] 
    mpart = pi_6*dpmid(1:fn2b)**3*1500.

  END SUBROUTINE set_coagc


  !*********************************************************************
  !
  !  function COAGC(diam1,diam2,mass1,mass2,temp,pres)
  !
  !*********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates the coagulation coefficient for two
  ! colliding particles. 
  !
  !
  ! Method: 
  ! -------  
  ! Only Brownian coagulation taken into account.
  ! Transition regime correction is done with Fuchs  
  ! flux matching.
  !
  !
  ! Interface:
  ! ----------
  ! Called from subroutine SET_COAG
  ! (only at the beginning of simulation)
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005 
  !
  !----------------------------------------------------------------------

  FUNCTION coagc(diam1,diam2,mass1,mass2,temp,pres)

    USE mo_kind, ONLY : dp

    USE mo_math_constants, ONLY : pi

    USE mo_physical_constants, ONLY: ak, p0sl_bg

    IMPLICIT NONE

    !-- Input variables ----------
    REAL(dp), INTENT(IN) :: &
         diam1,  &   ! diameters of colliding particles [m]
         diam2,  &   !
         mass1,  &   ! masses -"- [kg]
         mass2,  &
         temp,   &   ! ambient temperature [K]
         pres        ! ambient pressure [fxm]

    !-- Output variables ---------
    REAL(dp) ::  &
         coagc       ! coagulation coefficient of particles [m3/s]

    !-- Local variables ----------  
    REAL(dp) ::  &
         visc,   &   ! viscosity of air [kg/(m s)]
         mfp,    &   ! mean free path of air molecules [m]
         mdiam,  &   ! mean diameter of colliding particles [m]
         fmdist      ! distance of flux matching [m]

    REAL(dp), DIMENSION (2) :: &
         diam,   &   ! diameters of particles [m]
         mpart,  &   ! masses of particles [kg]
         knud,   &   ! particle knudsen number [1]
         beta,   &   ! Cunningham correction factor [1]
         dfpart, &   ! particle diffusion coefficient [m2/s]
         mtvel,  &   ! particle mean thermal velocity [m/s]
         omega,  &   !
         tva,    &   ! temporary variable [m]
         flux        ! flux in continuum and free molec. regime [m/s]

    !------------------------------------------------------------------------------- 

    !-- 0) Initializing particle and ambient air variables --------------------
    diam = (/ diam1, diam2 /)       ! particle diameters [m]
    mpart = (/ mass1, mass2 /)       ! particle masses [kg]

    visc = (7.44523e-3_dp*temp**1.5_dp)/ &
         (5093._dp*(temp+110.4_dp))                   ! viscosity of air [kg/(m s)]
    mfp = (1.656e-10_dp*temp+1.828e-8_dp)*p0sl_bg/pres ! mean free path of air [m]


    !-- 2) Slip correction factor for small particles -------------------------

    knud = 2._dp*mfp/diam                                    ! Knudsen number
    beta = 1._dp+knud*(1.142_dp+0.558_dp*exp(-0.999_dp/knud))! Cunningham correction factor
    ! (Allen and Raabe, Aerosol Sci. Tech. 4, 269)

    !-- 3) Particle properties ------------------------------------------------

    dfpart = beta*ak*temp/(3._dp*pi*visc*diam)  ! diffusion coefficient [m2/s]
    mtvel = sqrt((8._dp*ak*temp)/(pi*mpart))    ! mean thermal velocity [m/s]
    omega = 8._dp*dfpart/(pi*mtvel)

    mdiam = 0.5_dp*(diam(1)+diam(2))               ! mean diameter [m]


    !-- 4) Calculation of fluxes and flux matching ----------------------------

    flux(1) = 4._dp*pi*mdiam*(dfpart(1)+dfpart(2)  )    ! flux in continuum regime [m3/s]
    flux(2) = pi*sqrt(mtvel(1)**2+mtvel(2)**2)*mdiam**2 !  -"- in free molec. regime [m3/s]

    tva(1) = ((mdiam+omega(1))**3 - &              ! temporary variable [m]
         (mdiam**2+omega(1)**2)* &
         sqrt((mdiam**2+omega(1)**2)))/ &
         (3._dp*mdiam*omega(1)) - mdiam

    tva(2) = ((mdiam+omega(2))**3 - &              ! temporary variable [m]
         (mdiam**2+omega(2)**2)* &
         sqrt((mdiam**2+omega(2)**2)))/ &
         (3._dp*mdiam*omega(2)) - mdiam

    fmdist = sqrt(tva(1)**2+tva(2)**2)             ! flux matching distance [m]


    !-- 5) Coagulation coefficient [m3/s] -------------------------------------

    coagc = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2)) 


  END FUNCTION coagc


  SUBROUTINE salsa_initialize

    ! Purpose:
    ! ---------
    ! Initializes constants and parameters 
    ! used in the SALSA aerosol model.
    !
    !
    ! Interface:
    ! ---------
    ! *salsa_initialize*  is called from *start_ham* in mo_ham_init
    !

    USE mo_ham,              ONLY: sizeclass

    USE mo_ham_salsactl
    USE mo_ham_wetdep_data, ONLY: csr_strat_wat, csr_strat_mix, csr_strat_ice, &
                                  csr_conv, cbcr, cbcs
    USE mo_exception, ONLY: finish
    USE mo_ham, ONLY:naeroclass,nham_subm
    USE mo_kind, ONLY: dp

    IMPLICIT NONE

    !--- Local variables:

    INTEGER :: jclass, jj
    INTEGER :: ii, jbin
   
    REAL    :: k1

    INTEGER :: wetdepswitch
  
    CHARACTER(len=3) :: cbin(fn2b)

    LOGICAL          :: lofine(fn2b)

    ! Subrange 1a

    DO jclass = in1a, fn1a

       jj = jclass - in1a + 1

       sizeclass(jclass)%classname   = "Size bin 1a"//CHAR(jj+48)
       sizeclass(jclass)%shortname  = "1a"//CHAR(jj+48)
       sizeclass(jclass)%self       = 1
       sizeclass(jclass)%lsoluble   = .TRUE.
       sizeclass(jclass)%lsed       = .FALSE.
       sizeclass(jclass)%lsoainclass = .TRUE. ! thk: VBS
       sizeclass(jclass)%lactivation = .FALSE. !>>dod<<

    END DO

    ! Subrange 2a

    DO jclass = in2a, fn2a

       jj = jclass - in2a + 1

       sizeclass(jclass)%classname   = "Size bin 2a"//CHAR(jj+48)
       sizeclass(jclass)%shortname  = "2a"//CHAR(jj+48)
       sizeclass(jclass)%self       = 1
       sizeclass(jclass)%lsoluble   = .TRUE.
       sizeclass(jclass)%lsed       = .TRUE.
       sizeclass(jclass)%lsoainclass = .TRUE. ! thk: VBS
       sizeclass(jclass)%lactivation = .TRUE. !>>dod<<

    END DO

    ! Subrange 2b

    DO jclass = in2b, fn2b

       jj = jclass - in2b + 1

       sizeclass(jclass)%classname   = "Size bin 2b"//CHAR(jj+48)
       sizeclass(jclass)%shortname  = "2b"//CHAR(jj+48)
       sizeclass(jclass)%self       = 1
       sizeclass(jclass)%lsoluble   = .FALSE.
       sizeclass(jclass)%lsed       = .TRUE.
       sizeclass(jclass)%lsoainclass = .FALSE.
       sizeclass(jclass)%lactivation = .TRUE. !>>dod<<

    END DO

 !--------------------------------------------------------
    !-- Prescribed scavenging ratios  (mod_aero_salsa.f90)
    !-- Prescribed mean mass scavenging coefficients


    IF (.NOT. ALLOCATED(csr_strat_wat)) ALLOCATE(csr_strat_wat(naeroclass(nham_subm)))
    IF (.NOT. ALLOCATED(csr_strat_mix)) ALLOCATE(csr_strat_mix(naeroclass(nham_subm)))
    IF (.NOT. ALLOCATED(csr_strat_ice)) ALLOCATE(csr_strat_ice(naeroclass(nham_subm)))
    IF (.NOT. ALLOCATED(csr_conv)) ALLOCATE(csr_conv(naeroclass(nham_subm)))
    IF (.NOT. ALLOCATED(cbcr)) ALLOCATE(cbcr(naeroclass(nham_subm)))
    IF (.NOT. ALLOCATED(cbcs)) ALLOCATE(cbcs(naeroclass(nham_subm)))

    wetdepswitch=1

    cbcs(:)=5e-3_dp

    if (wetdepswitch==1) then !stier et al wetdep coefficients

       ii = in1a
       cbin(ii) = '1a'//CHAR(48+ii)
       lofine(ii) = .TRUE.
       csr_strat_wat(ii) = 0.1_dp
       csr_strat_mix(ii) = 0.1_dp
       csr_strat_ice(ii) = 0.1_dp
       csr_conv(ii)      = 0.2_dp
       cbcr(ii)          = 5.0E-4_dp

       

       DO ii = in1a+1, fn1a

          cbin(ii) = '1a'//CHAR(48+ii)
          lofine(ii) = .TRUE.
          csr_strat_wat(ii) = 0.25_dp
          csr_strat_mix(ii) = 0.4_dp
          csr_strat_ice(ii) = 0.1_dp
          csr_conv(ii)      = 0.6_dp
          cbcr(ii)          = 1.e-4_dp

       END DO

       DO ii = in2a, in2a+(nbin2-1)

          cbin(ii) = '2a'//CHAR(49+ii-in2a)
          lofine(ii) = .TRUE.
          csr_strat_wat(ii) = 0.85_dp
          csr_strat_mix(ii) = 0.75_dp
          csr_strat_ice(ii) = 0.1_dp
          csr_conv(ii)      = 0.99_dp
          cbcr(ii)          = 1.0E-3_dp

       END DO
       DO ii = fn2a-(nbin3-1), fn2a

          cbin(ii) = '2a'//CHAR(49+ii-in2a)
          lofine(ii) = .TRUE.
          csr_strat_wat(ii) = 0.99_dp
          csr_strat_mix(ii) = 0.75_dp
          csr_strat_ice(ii) = 0.1_dp
          csr_conv(ii)      = 0.99_dp
          cbcr(ii)          = 1.0E-1_dp

       END DO

       DO ii = in2b, in2b+(nbin2-1)

          cbin(ii) = '2b'//CHAR(49+ii-in2b)
          lofine(ii) = .FALSE.
          csr_strat_wat(ii) = 0.2_dp
          csr_strat_mix(ii) = 0.1_dp
          csr_strat_ice(ii) = 0.1_dp
          csr_conv(ii)      = 0.2_dp
          cbcr(ii)          = 1.0E-4_dp

       END DO
       DO ii = fn2b-(nbin3-1), fn2b

          cbin(ii) = '2b'//CHAR(49+ii-in2b)
          lofine(ii) = .FALSE.
          csr_strat_wat(ii) = 0.40_dp
          csr_strat_mix(ii) = 0.40_dp
          csr_strat_ice(ii) = 0.1_dp
          csr_conv(ii)      = 0.40_dp
          cbcr(ii)          = 1.0E-1_dp

       END DO

    elseif(wetdepswitch==0) then ! linear wetdeposition coefficients

       DO ii = in1a, fn1a

          cbin(ii) = '1a'//CHAR(48+ii)
          lofine(ii) = .TRUE.
          k1 = (0.25_dp - 0.06_dp)/(real(fn1a) - real(in1a))
          csr_strat_wat(ii) = k1 * (real(ii) - real(in1a)) + 0.06_dp
          !k1 = (0.40_dp - 0.10_dp)/(real(fn1a) - real(in1a))
          !csr_strat_mix(ii) = k1 * (real(ii) - real(in1a)) + 0.1_dp
          ! new coefficients by Bourgeois & Bey 2011 jgr
          csr_strat_mix(ii) = 0.06_dp
          csr_strat_ice(ii) = 0.06_dp
          k1 = (0.60_dp - 0.20_dp)/(real(fn1a) - real(in1a))
          csr_conv(ii)      = k1 * (real(ii) - real(in1a)) + 0.2_dp
          k1 = (1.0E-4_dp - 5.0E-4_dp)/(real(fn1a) - real(in1a))
          !--- Mean mass scavenging coefficients normalized by rain-rate [kg m-2]:
          !    Rain: Seinfeld & Pandis, Fig 20.15:
          cbcr(ii)          = k1 * (real(ii) - real(in1a)) + 5.0E-4_dp

       END DO

       DO ii = in2a, in2a+(nbin2-1)

          cbin(ii) = '2a'//CHAR(49+ii-in2a)
          lofine(ii) = .TRUE.
          k1 = (0.85_dp - 0.25_dp)/(real(fn2a) - real(in2a))
          csr_strat_wat(ii) = k1 * (real(ii) - real(in2a)) + 0.25_dp
          !k1 = (0.75_dp - 0.40_dp)/(real(fn2a) - real(in2a))
          !csr_strat_mix(ii) = k1 * (real(ii) - real(in2a)) + 0.4_dp
          csr_strat_mix(ii) = 0.06_dp
          csr_strat_ice(ii) = 0.06_dp
          k1 = (0.99_dp - 0.60_dp)/(real(fn2a) - real(in2a))
          csr_conv(ii)      = k1 * (real(ii) - real(in2a)) + 0.6_dp
          k1 = (1.0E-3_dp - 1.0E-4_dp)/(real(fn2a) - real(in2a))
          cbcr(ii)          = k1 * (real(ii) - real(in2a)) + 1.0E-4_dp

       END DO
       DO ii = fn2a-(nbin3-1), fn2a

          cbin(ii) = '2a'//CHAR(49+ii-in2a)
          lofine(ii) = .TRUE.
          k1 = (0.99_dp - 0.85_dp)/(real(fn2a) - real(fn2a-2))
          csr_strat_wat(ii) = k1 * (real(ii) - real(fn2a-2)) + 0.85_dp
          k1 = (0.75_dp - 0.75_dp)/(real(fn2a) - real(fn2a-2))
          csr_strat_mix(ii) = k1 * (real(ii) - real(fn2a-2)) + 0.75_dp
!          csr_strat_mix(ii) = 0.
          csr_strat_ice(ii) = 0.06_dp
          k1 = (0.40_dp - 0.99_dp)/(real(fn2a) - real(fn2a-2))
          csr_conv(ii)      = k1 * (real(ii) - real(fn2a-2)) + 0.99_dp
          k1 = (1.0E-1_dp - 1.0E-3_dp)/(real(fn2a) - real(fn2a-2))
          cbcr(ii)          = k1 * (real(ii) - real(fn2a-2)) + 1.0E-3_dp

       END DO

       DO ii = in2b, in2b+(nbin2-1)

          cbin(ii) = '2b'//CHAR(49+ii-in2b)
          lofine(ii) = .FALSE.
          k1 = (0.40_dp - 0.20_dp)/(real(fn2b) - real(in2b))
          csr_strat_wat(ii) = k1 * (real(ii) - real(in2b)) + 0.20_dp
          csr_strat_mix(ii) = 0.06_dp
          csr_strat_ice(ii) = 0.06_dp
          k1 = (0.40_dp - 0.20_dp)/(real(fn2b) - real(in2b))
          csr_conv(ii)      = k1 * (real(ii) - real(in2b)) + 0.20_dp
          k1 = (1.0E-3_dp - 1.0E-4_dp)/(real(fn2b) - real(in2b))
          cbcr(ii)          = k1 * (real(ii) - real(in2b)) + 1.0E-4_dp

       END DO
       DO ii = fn2b-(nbin3-1), fn2b

          cbin(ii) = '2b'//CHAR(49+ii-in2b)
          lofine(ii) = .TRUE.
          k1 = (0.40_dp - 0.40_dp)/(real(fn2b) - real(fn2b-2))
          csr_strat_wat(ii) = k1 * (real(ii) - real(fn2b-2)) + 0.40_dp
          k1 = (0.40_dp - 0.40_dp)/(real(fn2b) - real(fn2b-2))
          csr_strat_mix(ii) = k1 * (real(ii) - real(fn2b-2)) + 0.40_dp
          csr_strat_ice(ii) = 0.06_dp
          k1 = (0.40_dp - 0.40_dp)/(real(fn2b) - real(fn2b-2))
          csr_conv(ii)      = k1 * (real(ii) - real(fn2b-2)) + 0.40_dp
          k1 = (1.0E-1_dp - 1.0E-3_dp)/(real(fn2b) - real(fn2b-2))
          cbcr(ii)          = k1 * (real(ii) - real(fn2b-2)) + 1.0E-3_dp

       END DO

    elseif(wetdepswitch==2) then

       ii = in1a
       cbin(ii) = '1a'//CHAR(48+ii)
       lofine(ii) = .TRUE.
       csr_strat_wat(ii) = 0.06_dp
       csr_strat_mix(ii) = 0.06_dp
       csr_strat_ice(ii) = 0.06_dp
       csr_conv(ii)      = 0.2_dp
       cbcr(ii)          = 5.0E-4_dp

       DO ii = in1a+1, fn1a

          cbin(ii) = '1a'//CHAR(48+ii)
          lofine(ii) = .TRUE.
          csr_strat_wat(ii) = 0.25_dp
          csr_strat_mix(ii) = 0.06_dp
          csr_strat_ice(ii) = 0.06_dp
          csr_conv(ii)      = 0.6_dp
          cbcr(ii)          = 1.e-4_dp

       END DO

       DO ii = in2a, in2a+(nbin2-1)
          cbin(ii) = '2a'//CHAR(49+ii-in2a)
          lofine(ii) = .TRUE.
          csr_strat_wat(ii) = 0.85_dp
          csr_strat_mix(ii) = 0.06_dp
          csr_strat_ice(ii) = 0.06_dp
          csr_conv(ii)      = 0.99_dp
          cbcr(ii)          = 1.0E-3_dp

       END DO
       DO ii = fn2a-(nbin3-1), fn2a
          cbin(ii) = '2a'//CHAR(49+ii-in2a)
          lofine(ii) = .TRUE.
          csr_strat_wat(ii) = 0.99_dp
          csr_strat_mix(ii) = 0.75_dp
          csr_strat_ice(ii) = 0.06_dp
          csr_conv(ii)      = 0.99_dp
          cbcr(ii)          = 1.0E-1_dp

       END DO

       DO ii = in2b, in2b+(nbin2-1)

          cbin(ii) = '2b'//CHAR(49+ii-in2b)
          lofine(ii) = .FALSE.
          csr_strat_wat(ii) = 0.2_dp
          csr_strat_mix(ii) = 0.06_dp
          csr_strat_ice(ii) = 0.06_dp
          csr_conv(ii)      = 0.2_dp
          cbcr(ii)          = 1.0E-4_dp

       END DO
       DO ii = fn2b-(nbin3-1), fn2b

          cbin(ii) = '2b'//CHAR(49+ii-in2b)
          lofine(ii) = .TRUE.
          csr_strat_wat(ii) = 0.40_dp
          csr_strat_mix(ii) = 0.06_dp
          csr_strat_ice(ii) = 0.06_dp
          csr_conv(ii)      = 0.40_dp
          cbcr(ii)          = 1.0E-1_dp

       END DO

    else
       call finish('mo_ham_salsa_init','wetdep not defined correctly')
    endif

    CALL set_sizebins
        
  END SUBROUTINE salsa_initialize

END MODULE mo_ham_salsa_init
