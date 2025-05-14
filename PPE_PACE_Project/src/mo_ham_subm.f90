!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_subm.f90
!!
!! \brief
!! Module to provide interface to the aerosol microphysics schemes
!!
!! \author Martin G. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Sylvaine Ferrachat, sylvaine.ferrachat@env.ethz.ch
!!
!! \revision_history
!!   -# The original code is from J. Feichter, J. Wilson and E. Vignatti, JRC Ispra 
!! and was adapted for ECHAM by P. Stier, Oxford. Other contributions include 
!! D. O'Donnell, K. Zhang and others
!!   -# M.G. Schultz (FZ Juelich) - new code structure for integration into echam6-hammoz (2009-09-24)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!   -# H. Kokkola (FMI) - generalization to include also SALSA aerosol microphysics (2014)
!!
!! \limitations
!! None
!!
!! \details
!! This module contains the ham_subm_interface routine and common routines
!! for aerosol microphysics interfaces.
!! This module contains the following subroutines which used to be individual files.
!!       ham_subm_interface
!!       subm_mass_sum
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

MODULE mo_ham_subm

IMPLICIT NONE

PRIVATE

PUBLIC :: ham_subm_interface


CONTAINS

! ---------------------------------------------------------------------------
!  ham_subm_interface: called from mo_submodel_interface 
! ---------------------------------------------------------------------------

SUBROUTINE ham_subm_interface(kproma, kbdim, klev, krow, ktrac, &
                        pap,    paph,                     &
                        pt,     pq,                       &
                        pxtm1,  pxtte,                    &
                        paclc,  pgrvolm1, ppbl)
  !
  ! Authors:
  ! --------
  !
  ! 1997 J. Feichter and J. Wilson (original code)
  ! 2000 J. Feichter and E. Vignati (adopted to m7)
  ! 2001 P. Stier
  ! 2008 J. Kazil
  !   - Ionization rate calculation
  !   - Changed m7 interface and call
  !   - Extended comments
  ! 2008 D.O'Donnell, MPI-Met, Hamburg
  !   - generalisation of the species (no more hardcoded so4, bc, oc, etc.)
  !   - added sizeclass and new modules, removed dependency on mo_ham_m7_trac,
  ! 2009 Kai Zhang, MPI-Met, Hamburg
  !   - submodel interface
  ! 2013 Harri Kokkola, FMI, Kuopio
  !   - implementation of SALSA microphysics
  !   - M7 microphysics routines separated to mo_ham_m7
  !
  !
  ! Method:
  ! -------
  !
  ! The interface takes the tracer mixing ratios of the at the time t-dt, adds
  ! the tendencies (pxtte) of the preceeding processes, converts the units,
  ! calculates ambient and auxiliary quantities, and hands the appropriate
  ! quantities to microphysics modules. A mass conservation check is performed 
  ! (before and after calling the microphysics) if requested.
  !
  ! After the call of microfysical module, the returning concentrations are 
  ! reconverted to mixing ratios, the tendencies of the mixing ratios are 
  ! calculated and returned to ECHAM.
  
  USE mo_kind,               ONLY: dp
  USE mo_physical_constants, ONLY: rd, vtmpc1, grav, avo, argas
  USE mo_time_control,       ONLY: time_step_len
  USE mo_echam_convect_tables, ONLY: tlucuaw, jptlucu1, jptlucu2
  USE mo_exception,          ONLY: finish, message, message_text, em_warn
  USE mo_species,            ONLY: speclist
  USE mo_ham,                ONLY: lgcr,  &
                                   nsolact,                                        &
                                   subm_ngasspec, subm_gasspec, subm_gasunitconv,        &
                                   subm_naerospec, subm_aerospec, subm_aerounitconv,     &
                                   subm_aero_idx,                                    &
                                   immr2ug, immr2molec, ivmr2molec
  USE mo_ham_subm_species,   ONLY: isubm_so4g, isubm_ocnv
  USE mo_ham,                ONLY: naerocomp, aerocomp, aerowater, mw_so4, mw_oc, &
                                   sizeclass, nclass, &
                                   nham_subm,       &
                                   HAM_BULK,        &
                                   HAM_M7,          &
                                   HAM_SALSA, nsol, &
                                   lmass_diag
  USE mo_ham_streams,        ONLY: rdry, rwet, densaer, relhum, ipr 
  USE mo_tracer_processes,   ONLY: xt_borrow
  USE mo_ham_gcrion,         ONLY: solar_activity,gcr_ionization
  !>>dod timers
  USE mo_control,            ONLY: ltimer
  USE mo_hammoz_timer,       ONLY: timer_start, timer_stop, timer_ham_m7_main
  USE mo_ham_salsa,          ONLY: salsa
  USE mo_ham_m7,             ONLY: m7
  !<<dod

  IMPLICIT NONE 
   
! compulsory arguments
  
  INTEGER :: kproma                     ! geographic block number of locations
  INTEGER :: kbdim                      ! geographic block maximum number of locations
  INTEGER :: klev                       ! numer of levels
  INTEGER :: ktrac                      ! number of tracers
  
! optional arguments

  INTEGER,  OPTIONAL :: krow                       ! geographic block number
  REAL(dp), OPTIONAL :: pap     (kbdim,klev)       ! pressure [Pa], at full levels
  REAL(dp), OPTIONAL :: paph    (kbdim,klev+1)     ! pressure [Pa], at half levels
  REAL(dp), OPTIONAL :: pt      (kbdim,klev)       ! temperature [K]  
  REAL(dp), OPTIONAL :: pq      (kbdim,klev)       ! specific humidity [kg/kg] 
  REAL(dp), OPTIONAL :: pxtm1   (kbdim,klev,ktrac) ! tracer mass/number mixing ratio at 
  REAL(dp), OPTIONAL :: pxtte   (kbdim,klev,ktrac) ! tracer mass/number mixing ratio tendencies [kg/kg s-1 or #/kg s-1]
  REAL(dp), OPTIONAL :: paclc   (kbdim,klev)       ! cloud cover [0,1]
  REAL(dp), OPTIONAL :: pgrvolm1(kbdim,klev)       ! grid box volume [m3]
  REAL(dp), OPTIONAL :: ppbl    (kbdim)            ! Planetary boundary layer top level
 
  !
  ! Streams and stream elements:
  !
  
  REAL(dp), POINTER :: rdry_p(:,:,:)
  REAL(dp), POINTER :: rwet_p(:,:,:)
  REAL(dp), POINTER :: densaer_p(:,:,:)
  
  !
  ! Local variables:
  !
  
  INTEGER :: it,jl,jk,jc,jn,jt,jspec,jclass,ilevp1
  
  REAL(dp):: ztmst, zqtmst, zqs, zq_amb
  
  REAL(dp):: zfac,     zqfac,            &
             zfacm,    zqfacm,           &
             zfacc,    zqfacc,           &
             zfacn,    zqfacn,           &
             zeps,     zaclc

  REAL(dp) :: zfac_vmr                     ! Conversion of vmr to molecules cm-3   

  !>>dod changed index to naerospec
  REAL(dp):: zmass_pre(subm_naerospec),       & ! mass of aerosol compounds before and
             zmass_post(subm_naerospec)         ! after microphysics call, for mass conservation check
  !<<dod

  REAL(dp):: zgso4(kbdim,klev),          & ! [H2SO4(g)] [cm-3]
             zgso4m1(kbdim,klev),        & ! [H2SO4(g)] at t-1 [cm-3]
             zgso4p1(kbdim,klev)           ! [H2SO4(g)] at t+1 [cm-3]

  !>>dod soa
  REAL(dp) :: zgas(kbdim,klev,subm_ngasspec)
  REAL(dp) :: zunitfac, zqunitfac, zfac1
  !<<dod

  REAL(dp):: zdgso4(kbdim,klev)            ! d[H2SO4(g)]/dt [cm-3 s-1]

  REAL(dp):: zrh(kbdim,klev),            & ! relative humidity at t+1 [0,1]
             zrhoa(kbdim,klev),          & ! air mass density [kg m-3]
             zipr(kbdim,klev),           & ! ionization rate [cm-3 s-1]
             zdpg(kbdim,klev),           & ! mass of air column in layer [kg m-2]
             zdz(kbdim,klev)               ! layer thickness [m]
             
  REAL(dp):: zaerml(kbdim,klev,naerocomp), & ! aerosol mass for individual compounds [molec. cm-3 for sulfate and OCNV, 
                                             ! and ug m-3 for bc, oc, ss, and dust]
             zaernl(kbdim,klev,nclass),    & ! aerosol number for each mode [cm-3]
             zm6rp(kbdim,klev,nclass),     & ! mean mode actual radius (wet for soluble and dry for insoluble modes) [cm]
             zm6dry(kbdim,klev,nsol),    & ! dry radius for soluble modes [cm]
             zrhop(kbdim,klev,nclass),     & ! mean mode particle density [g cm-3]
             zww(kbdim,klev,nclass),       &   ! aerosol water content for each mode [kg(water) m-3(air)]
             zaervl(kbdim,klev,naerocomp) ! aerosol mass for individual compounds [molec. cm-3 for sulfate and ug m-3 for bc, oc, ss, and dust] !alaak
  
  REAL(dp):: zsolact                       ! Solar activity parameter [-1,1]
  
  !>>dod from array
  LOGICAL :: labort
  !<<dod

  !--- 0) Initialisations: -----------------------------------------------------
  
  zmass_pre  = 0._dp
  zmass_post = 0._dp
  
  ztmst      = time_step_len
  zqtmst     = 1._dp/time_step_len
  
  zeps       = 1.E-10_dp
  
  ilevp1     = klev+1
  
  !>>dod
  labort = .FALSE.
  !<<dod

  !--- Mass conserving correction of negative tracer values:
  
  CALL xt_borrow(kproma, kbdim,  klev, ilevp1, ktrac, &
                 pap,    paph,                        &
                 pxtm1,  pxtte                        )
  
  !--- 1) Calculate necessary parameters: --------------------------------------
  
!!mgs(S)!!  !--- Factor to transform mass sulfur in kg into molecules per kg:
!!mgs(S)!!  zfacm  = 6.022e+20_dp/32._dp
  !--- Factor to transform mass SO4 in kg into molecules per kg:
  zfacm  = 6.022e+20_dp/mw_so4

  !--- Factor to transform mass OC in kg into molecules per kg:
  zfacc  = 6.022e+20_dp/mw_oc

  !--- Factor to transform kg into micro gram:

  zfac   = 1.e09_dp

  !--- Factor to transform N/m**3 into N/cm**3:

  zfacn  = 1.0e-06_dp

  !--- Prefactor used when converting VMR into molec cm-3
  zfac_vmr = 1.E-6_dp*avo/argas

  zqfac  = 1.0_dp/zfac
  zqfacm = 1.0_dp/zfacm
  zqfacc = 1.0_dp/zfacc
  zqfacn = 1.0_dp/zfacn

  !--- 2) Calculate ambient properties: ----------------------------------------
  
  DO jk = 1,klev
    DO  jl = 1,kproma
      
      !--- 2.1) Calculate air density:
      !         (currently neglects volume occupied by liquid and ice water  = > physc)
      
      zrhoa(jl,jk) = pap(jl,jk)/(pt(jl,jk)*rd*(1._dp+vtmpc1*pq(jl,jk)))
      
      !--- 2.2) New calculation of the relative humidity (over water):
      
      it    = NINT(pt(jl,jk)*1000._dp)
      it    = MAX(MIN(it,jptlucu2),jptlucu1)
      
      zqs = tlucuaw(it)/pap(jl,jk)
      zqs = MIN(zqs,0.5_dp)
      zqs = zqs/(1._dp-vtmpc1*zqs)
      
      zaclc = MIN(paclc(jl,jk),1.0_dp-zeps)
      
      zq_amb = MAX( 0.0_dp , (pq(jl,jk)-zqs*zaclc)/(1._dp-zaclc) )
      
      zrh(jl,jk) = zq_amb/zqs
      
      zrh(jl,jk) = MAX(0.0_dp,MIN(zrh(jl,jk),1.0_dp))
      
      !--- Air mass auxiliary variable:
      
      zdpg(jl,jk) = (paph(jl,jk+1)-paph(jl,jk))/grav
      
      !--- Layer thickness zdz = dp/(rho*grav) [m]:
      
      zdz(jl,jk) = zdpg(jl,jk)/zrhoa(jl,jk)
      
    ENDDO
    
    !ham_ps:changed position to suppress loop-fusion on NEC-SX6
    
    relhum(1:kproma,jk,krow) = zrh(1:kproma,jk)
     
  ENDDO
  
  !--- 2.3) Ionization rate
  !
  
  IF (lgcr) then
    
    !--- Solar activity:
    
    IF (abs(nsolact) > 1.0_dp) then
      ! Parameterize solar activity as function of the current date:
      zsolact = solar_activity()
    ELSE
      ! Use solar activity parameter set by user:
      zsolact = nsolact
    ENDIF
    
    !--- Galactic cosmic ray ionization rate:
    
    CALL gcr_ionization(krow,kproma,kbdim,klev,zsolact,pt,pap,zipr)
    
    !--- Save the ionization rate:
    
    ipr(1:kproma,1:klev,krow) = 1.0e6_dp*zipr(1:kproma,1:klev)
    
  ENDIF
  
  !--- 3) Convert units and add the tendencies of preceeding proceses: ---------
  !       Convert:
!!mgs(S)!!  !         - mass of sulfur species from [mass(S)/mass(air)]
!!mgs(S)!!  !           to [molecules/cm+3]
  !         - mass of sulfur species from [mass(SO4)/mass(air)]
  !         - mass of all other species to micro-gram/cubic-meter
  !         - particle numbers are converted from [N/kg(air)] to [N/cm+3]
  !
  !--- 3.1) Gases:
  
!++mgs: Note that logic changed here - looping only over subm_gasspec !! 
!! this is equivalent to the former "IF (gasspec(jn)%lm7gas) THEN"
!! the former "gasspec(jn)%m7unitconv" is now replaced with subm_gasunitconv(jn)
!--mgs
!++mgs: initialize zgso4p1
  zgso4p1(:,:) = 0._dp

  DO jn=1,subm_ngasspec
    jt = speclist(subm_gasspec(jn))%idt

    SELECT CASE(subm_gasunitconv(jn))
    CASE(immr2ug)
      DO jk=1,klev
        DO jl=1,kproma
          zunitfac = zfac*zrhoa(jl,jk)
          zgas(jl,jk,jn) = zunitfac*(pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst)
        END DO
      END DO

    CASE(immr2molec)
      DO jk=1,klev
        DO jl=1,kproma
           ! >> thk #513
           zunitfac = 1e-3*zrhoa(jl,jk)*avo/speclist(subm_gasspec(jn))%moleweight
           ! << thk
           zgas(jl,jk,jn) = zunitfac*(pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst)
        END DO
      END DO
           
      !---need gas phase SO4 tendency for Jan Kazil's SO4 condensation scheme
      IF (jn == isubm_so4g) THEN
        DO jk=1,klev
          DO jl=1,kproma
            zunitfac = zfacm*zrhoa(jl,jk)
            zgso4m1(jl,jk) = zunitfac*pxtm1(jl,jk,jt)
            zdgso4(jl,jk) = zunitfac*pxtte(jl,jk,jt)
            zgso4p1(jl,jk) = zgso4m1(jl,jk)+zdgso4(jl,jk)*time_step_len
          END DO
        END DO
      END IF
              
    CASE(ivmr2molec)
      DO jk=1,klev
        DO jl=1,kproma
          zunitfac = zfac_vmr*pap(jl,jk)/pt(jl,jk)
          zgas(jl,jk,jn) = zunitfac*(pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst)
        END DO
      END DO
    END SELECT

  END DO
  
  !--- 3.2 Particle mass:
  
  DO jn=1, naerocomp
     jt    = aerocomp(jn)%idt       ! get tracer id
     jspec = aerocomp(jn)%spid      ! get species id
     jl    = subm_aero_idx(jspec)     ! get index to subm_aerospec list
     !!mgs=old code!!     IF (aerocomp(jn)%species%m7unitconv == immr2molec) THEN
     IF (jl <= 0) THEN
        WRITE(message_text,*) 'SUBM_AERO_IDX Mapping error !! No index for jspec=',jspec
        CALL finish('ham_subm_interface', message_text)
     END IF
     IF (subm_aerounitconv(jl) == immr2molec) THEN
        zfac1 = zfacm
     ELSE
        zfac1 = zfac
     END IF
     
     zaerml(1:kproma,:,jn) = zfac1*zrhoa(1:kproma,:)*(pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt)*ztmst)
     
  END DO
        
  !--- 3.3) Particle numbers:

  DO jn=1, nclass
     jt = sizeclass(jn)%idt_no     
     zaernl(1:kproma,:,jn) = zfacn*zrhoa(1:kproma,:)*(pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
  END DO
  !<<dod soa           

  !--- Discard potential remaining negative values:
  
  zgso4p1(1:kproma,:) = MAX(zgso4p1(1:kproma,:),0._dp)
  
  zaerml(1:kproma,:,:) = MAX(zaerml(1:kproma,:,:),0._dp)
  zaernl(1:kproma,:,:) = MAX(zaernl(1:kproma,:,:),0._dp)
  
  !--- Sum total mass of all compounds for mass diagnostics:
  !>>dod moved to separate subroutine 

  IF (lmass_diag) CALL subm_mass_sum(kbdim, kproma, klev, krow, zgas, zaerml, zdpg, zmass_pre)
  !<<dod

  !--- 4) Call of microphysics model (subm): ----------------------------------------------------------
  
  
  SELECT CASE(nham_subm)
     
  CASE(HAM_BULK)
     
     !  CALL bulk_model

  CASE(HAM_M7)

     ! The gas phase concentrations at t-1 and their derivatives for processes
     ! preceeding subm are passed to subm, as subm performs their time integration
     ! together with the processes it accounts for. The aerosol tracers (mass and
     ! number) at t+1 after processes preceeding subm are passed to microphysics module, 
     ! which applies its processes on these.
  
     zgso4(1:kproma,:) = zgso4m1(1:kproma,:)

     !>>dod timers
     IF (ltimer) CALL timer_start(timer_ham_m7_main)
     !<<dod

     CALL m7(kproma,  kbdim,   klev,    krow, &  ! ECHAM indices
             pap,     pt,      zrh,           &  ! Pressure, temperature, RH 
             zgas,                            &  ! Gases at t+dt
             zgso4,   zdgso4,                 &  ! [H2SO4(g)], derivative (production rate)
             zaerml,  zaernl,                 &  ! Aerosol mass and number
             zm6rp,   zm6dry,  zrhop,   zww,  &  ! Aerosol properties
             zipr,    paclc,                  &  ! Ionization rate, cloud cover
             zdz,     pgrvolm1,               &  ! Layer thickness, grid box volume 
             ppbl                              ) ! Planetary boundary layer top level

     !>>dod timers
     IF (ltimer) CALL timer_stop(timer_ham_m7_main)
     !<<dod

     !>>dod again special handling for so4 gas for Jan Kazil's SO4 condensation scheme
     
     zgas(1:kproma,:, isubm_so4g) = zgso4(1:kproma,:)

  CASE(HAM_SALSA)
     
     ! Number concentration converted from cm-3 to m-3 for SALSA
     zaernl(1:kproma,:,:) = zaernl(1:kproma,:,:) * 1.e6_dp

     ! >> thk: adapting for VBS
     ! Gas phase concentrations converted from cm-3 to m-3 for SALSA
     zgas(1:kproma,:,:) = zgas(1:kproma,:,:) * 1.e6_dp

     CALL salsa(kproma,  kbdim,   klev,    krow, &  ! ECHAM indices
                pap,     zrh,     pt,      ztmst,&  ! Pressure, RH, temperature, time step length
                !zgso4,   zgocnv,  zgocsv,        &  ! [H2SO4(g)], [OCNV(g)], [OCSV(g)]
                zgas,                            &  ! gas phase concentrations
                zaerml,  zaernl,                 &  ! Aerosol volume and number
                zm6rp,   zm6dry,  zrhop,   zww,  &  ! Aerosol properties
                ppbl                             &
                ) ! Planetary boundary layer top level

     ! Number concentration converted from m-3 to cm-3 for compatibility with M7
     zaernl(1:kproma,:,:) = zaernl(1:kproma,:,:) * 1.e-6_dp
     
     ! Gas phase concentrations converted from m-3 to cm-3 for compatibility with M7 
     zgas(1:kproma,:,:) = zgas(1:kproma,:,:) * 1.e-6_dp
     ! << thk

  END SELECT
  
  !--- 5) Reconvert masses and numbers into mixing ratios, other ---------------
  !       quantities to SI units and calculate the tendencies (xtte):
  !
  !>>dod soa
  !---5.1) Gases
!++mgs: changes to replace gasspec and m7unitconv (see above)
  DO jn = 1,subm_ngasspec
    jt=speclist(subm_gasspec(jn))%idt

    SELECT CASE(subm_gasunitconv(jn))    
    CASE(immr2ug)                            
      DO jk=1,klev
        DO jl=1,kproma
          zqunitfac = zqfac/zrhoa(jl,jk)
          pxtte(jl,jk,jt) = (zgas(jl,jk,jn)*zqunitfac-pxtm1(jl,jk,jt))*zqtmst
        END DO
      END DO
           
    CASE(immr2molec)
      DO jk=1,klev
        DO jl=1,kproma
           ! >> thk #513
           zqunitfac = 1e3*speclist(subm_gasspec(jn))%moleweight/(avo*zrhoa(jl,jk))
           ! << thk
           pxtte(jl,jk,jt) = (zgas(jl,jk,jn)*zqunitfac-pxtm1(jl,jk,jt))*zqtmst
        END DO
      END DO
    END SELECT

!>>csld #538
! csld : A more elegant way to resolve this bug could be to define an extra case
! for h2so4 (subm_gasunitconv(jn)). Another alternative would be to do the back unit conversion 
! inside the HAM_M7 case. 
! I let this cosmetic issue for the moment.
    IF ((nham_subm == HAM_M7) .AND. (jn == isubm_so4g)) THEN
        DO jk=1,klev
          DO jl=1,kproma
             zqunitfac = zqfacm/zrhoa(jl,jk)
             pxtte(jl,jk,jt) = (zgas(jl,jk,jn)*zqunitfac-pxtm1(jl,jk,jt))*zqtmst
          END DO
      END DO
    END IF
!<<csld 538

  END DO

  !--- 5.2) Particle mass:
  DO jn = 1,naerocomp
     jt    = aerocomp(jn)%idt       ! get tracer id
     jspec = aerocomp(jn)%spid      ! get secies id
     jl    = subm_aero_idx(jspec)     ! get index to subm_aerospec list

!!mgs-old code!!     SELECT CASE(aerocomp(jn)%species%m7unitconv)    

     SELECT CASE(subm_aerounitconv(jl))
     CASE(immr2ug)                            
        zqunitfac = zqfac
     CASE(immr2molec)
        zqunitfac = zqfacm
     END SELECT

     DO jk=1,klev
        DO jl=1,kproma
           pxtte(jl,jk,jt) = (zaerml(jl,jk,jn)*zqunitfac/zrhoa(jl,jk)-pxtm1(jl,jk,jt))*zqtmst
        END DO
     END DO
  END DO

  !--- 5.3 Particle numbers:

  DO jn=1,nclass
     jt = sizeclass(jn)%idt_no     
     DO jk=1,klev
        DO jl=1,kproma
           zqunitfac = zqfacn/zrhoa(jl,jk)
           pxtte(jl,jk,jt) = (zaernl(jl,jk,jn)*zqunitfac-pxtm1(jl,jk,jt))*zqtmst
        END DO
     END DO
  END DO
  !<<dod soa

  !--- 6) Convert microphysical model output quantities to SI units and store in streams: --------------
  
  DO jclass = 1, nclass
    
    rdry_p     => rdry(jclass)%ptr
    rwet_p     => rwet(jclass)%ptr
    densaer_p  => densaer(jclass)%ptr
    
    !--- Mean mode density from [g/cm3] to [kg/m3]:
    
    densaer_p(1:kproma,:,krow) = zrhop(1:kproma,:,jclass)*1.E3_dp
    
    !--- Ambient Count Median Radius from [cm] to [m]:
    
    rwet_p(1:kproma,:,krow) = zm6rp(1:kproma,:,jclass)/100._dp
    
    !--- Dry Count Median Radius from [cm] to [m]:
    
    IF (jclass <=  nsol) THEN 
       rdry_p(1:kproma,:,krow) = zm6dry(1:kproma,:,jclass)/100._dp
    ELSE
       rdry_p(1:kproma,:,krow) = zm6rp(1:kproma,:,jclass)/100._dp
    END IF
    
  END DO
  
  !--- Store diagnostic aerosol properties in pseudo-tracers:
  !>>dod soa
  DO jn=1,nsol
     pxtm1(1:kproma,:,aerowater(jn)%idt)=zww(1:kproma,:,jn)/zrhoa(1:kproma,:)
  END DO
  !<<dod

!gf
  !--- Mass conserving correction of negative tracer values:

  CALL xt_borrow(kproma, kbdim,  klev, ilevp1, ktrac, &
                 pap,    paph,                        &
                 pxtm1,  pxtte                        )
!gf

  !--- 7) Perform mass conservation check: -------------------------------------
  
  IF (lmass_diag) THEN
    
     !--- Sum total mass of all compounds for mass diagnostics:

     CALL subm_mass_sum(kbdim, kproma, klev, krow, zgas, zaerml, zdpg, zmass_post)

     !--- Perform mass conservation check:
     !>>dod soa rewritten to avoid hardcoding of individual species
     DO jc=1,subm_naerospec

        IF( ABS(zmass_pre(jc)-zmass_post(jc)) > 0.100_dp*ABS(MAX(zmass_pre(jc),zmass_post(jc)))) THEN

           CALL message('ham_subm_interface', 'microphysics module violates mass conservation by >10% for '// &
                        speclist(subm_aerospec(jc))%longname, level=em_warn)

           labort=.TRUE.        
        
        ELSE IF( ABS(zmass_pre(jc)-zmass_post(jc)) > 0.010_dp*ABS(MAX(zmass_pre(jc),zmass_post(jc)))) THEN

           CALL message('ham_subm_interface', 'microphysics module violates mass conservation by >1% for '// &
                        speclist(subm_aerospec(jc))%longname)

        ELSE IF( ABS(zmass_pre(jc)-zmass_post(jc)) > 0.001_dp*ABS(MAX(zmass_pre(jc),zmass_post(jc)))) THEN

           CALL message('ham_subm_interface', 'microphysics module violates mass conservation by >0.1% for '// &
                        speclist(subm_aerospec(jc))%longname)
        END IF

     END DO

  END IF
    
END SUBROUTINE ham_subm_interface

! ---------------------------------------------------------------------------

SUBROUTINE subm_mass_sum(kbdim, kproma, klev, krow, pgas, paerml, pdpg, pmasssum)
    
  ! Purpose: sum the total mass per species over the gridpoints on one processor
  ! This subroutine is called if the namelist variables nsoa == 1 and lomassdiag == .TRUE.
  ! It is called twice, from ham_subm_interface, once before and once after the 
  ! microphysical processes have been calculated.
  ! It tests the total mass and throws an exception if mass conservation is violated

  !---inherited functions, types and data
  USE mo_kind,             ONLY: dp
  USE mo_ham,              ONLY: naerocomp, aerocomp, nsoa, nsoaspec, nham_subm, HAM_SALSA
  USE mo_ham,              ONLY: subm_ngasspec, subm_naerospec
  USE mo_ham_subm_species, ONLY: isubm_so2, isubm_so4, isubm_so4g, isubm_oc, isubm_ocnv
  USE mo_geoloc,           ONLY: gboxarea_2d
  USE mo_ham_soa,          ONLY: soaprop

  IMPLICIT NONE

  !---subroutine interface
  INTEGER,  INTENT(IN) :: kbdim, kproma, klev, krow        ! grid parameters
  REAL(dp), INTENT(IN) :: pgas(kbdim,klev,subm_ngasspec)        ! gas concentrations
  REAL(dp), INTENT(IN) :: paerml(kbdim,klev,naerocomp)         ! aerosol concentrations
  REAL(dp), INTENT(IN) :: pdpg(kbdim,klev)                 ! air mass auxiliary variable
  REAL(dp), INTENT(OUT) :: pmasssum(subm_naerospec)             ! total mass

  !---Local data
  !   Parameters:
  !   -
  !   Local variables
  REAL(dp) :: zarea(kbdim)                     ! gridbox area
  INTEGER :: jspec                             ! species index
  INTEGER :: jl, jk, jn, jm                    ! loop counter


    !---executable procedure
    pmasssum(:) = 0._dp
    zarea(1:kproma) = gboxarea_2d(1:kproma,krow)

!!mgs(S)!!: ### update needed to obtain correct mass??

    !---sum aerosol masses per species
    DO jn=1,naerocomp
       !---get species index
       jspec = aerocomp(jn)%spid
       
       DO jk=1,klev
          DO jl=1,kproma
             pmasssum(jspec) = pmasssum(jspec) + paerml(jl,jk,jn)*pdpg(jl,jk)*zarea(jl)
          END DO
       END DO
       
    END DO
       
    !---sulphate: add SO2 and gas phase SO4
    DO jk=1,klev
       DO jl=1,kproma
          pmasssum(isubm_so4) = pmasssum(isubm_so4) + (pgas(jl,jk,isubm_so2)+pgas(jl,jk,isubm_so4g)) * &
                                                pdpg(jl,jk)*zarea(jl)
       END DO
    END DO

    !---SOA....
    IF (nsoa == 1) THEN
       DO jm=1,nsoaspec
          IF (soaprop(jm)%lvolatile) THEN
             !---gas species index
             jspec = soaprop(jm)%spid_soa

             DO jk=1,klev
                DO jl=1,kproma
                   pmasssum(jspec) = pmasssum(jspec)+pgas(jl,jk,jspec)*pdpg(jl,jk)*zarea(jl)
                END DO
             END DO
             
          END IF
       END DO
       
    END IF
    
  END SUBROUTINE subm_mass_sum
  !
  !*******************************************************************************
  !   


END MODULE mo_ham_subm
