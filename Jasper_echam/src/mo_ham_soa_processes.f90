!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_soa_processes.f90
!!
!! \brief
!! This file contains the heart of the SOA module: calculation of gas-phase and 
!! aerosol-phase equilibrium, 2-product model, etc.
!!
!! \author Declan O'Donnell (MPI-Met)
!!
!! \responsible_coder
!! Declan O'Donnell, declan.Odonnell@fmi.fi
!!
!! \revision_history
!!   -# Declan O'Donnell (MPI-Met) - original code (2009)
!!   -# Martin Schultz (FZ-Juelich) -  combined routines in module (2009-08-25)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!
!! \limitations
!! None
!!
!! \details
!! Note: These routines used to be separate subroutine files:
!!   - soa2prod.f90
!!   - soa_equi0.f90
!!   - soa_equi1.f90
!!   - soa_part.f90
!!   - soa_progn.f90
!!
!! \bibliographic_references
!! D. O'Donnell et al, Estimating the direct and indirect effects of secondary organic aerosols
!! using ECHAM5-HAM, Atmos. Chem. Phys., 11, 8635-8659, doi:10.5194/acp-11-8635-2011, (2011)
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

MODULE mo_ham_soa_processes

  IMPLICIT NONE

  !--- public member functions                                    
  PUBLIC :: soa_equi0               ! equilibrium SOA gas and aerosol masses, semi-volatile species
  PUBLIC :: soa2prod                ! implements a 2-product model of SOA formation from gas-phase precursors
  PUBLIC :: soa_progn               ! calculates the tendency of total SOA mass per species

  !---private member functions 
  PUBLIC :: soa_equi1               ! distribution of newly formed SOA among modes, non-volatile species
  PRIVATE :: soa_part               ! iterative calculation of equilibrium gas-aerosol masses 

CONTAINS

  SUBROUTINE soa2prod(kproma, kbdim, klev,  klevp1, krow, ktrac, paph, pap, &
                      pt,     pxtm1, pxtte)
    
    !   soa2prod implements a 2-product model of SOA formation from gas-phase precursors
    !
    !   Each SOA precursor species oxidises into 2 condensable SOA species of differing volatilities
    !   Precursor + oxidant -> a1[P1] + a2[P2] 
    !   where P1 and P2 are SOA species and a1, a2 are mass-based stoichiometric coefficients, which are
    !   obtained empirically. Here 'oxidant' may be O3, OH or NO3.
    !
    !   Author:
    !   ------
    !   Declan O'Donnell, MPI-M

    USE mo_kind,         ONLY: dp
    USE mo_physical_constants, ONLY: avo, argas, grav, amd 
    USE mo_time_control, ONLY: time_step_len, delta_time
    USE mo_geoloc,       ONLY: rdayl_x
    USE mo_boundary_condition, ONLY: bc_apply, bc_query
    USE mo_external_field_processor, ONLY: EF_FILE
    USE mo_species,      ONLY: speclist
    !>>dod simplified soa id handling
    USE mo_ham_soa,      ONLY: soaprop,                                              &
    !>>csld #404
                               id_apin, id_tbeta, id_bpin, id_lim, id_sab, id_myrc,  &
                               id_car, id_isop, id_tol, id_xyl, id_benz,             &
    !<<csld #404
                               d_prod_soa_mterp, d_prod_soa_tol, d_prod_soa_xyl,     &
                               d_prod_soa_isop, d_prod_soa_benz,                     &
                               d_chem_sink_mterp, d_chem_sink_isop, d_chem_sink_tol, &
                               d_chem_sink_xyl, d_chem_sink_benz
    !<<dod
    USE mo_ham,          ONLY: nsoalumping, nsoaspec, ibc_oh, ibc_o3, ibc_no3
    USE mo_ham_streams,  ONLY: daylength

    USE mo_ham_species,  ONLY: id_oh, id_o3, id_oc, &
                               id_no3 !gf see #146
    
    IMPLICIT NONE

    INTEGER,  INTENT(in)    :: kproma                    ! geographic block number of locations
    INTEGER,  INTENT(in)    :: kbdim                     ! geographic block maximum number of locations 
    INTEGER,  INTENT(in)    :: klev                      ! numer of levels 
    INTEGER,  INTENT(in)    :: klevp1                    ! numer of levels + 1  
    INTEGER,  INTENT(in)    :: ktrac                     ! number of tracers
    INTEGER,  INTENT(in)    :: krow                      ! geographic block number
    REAL(dp), INTENT(in)    :: paph(kbdim,klevp1)        ! air pressure at layer interface 
    REAL(dp), INTENT(in)    :: pap(kbdim,klev)           ! air pressure at layer center    
    REAL(dp), INTENT(in)    :: pt (kbdim,klev)           ! air temperature 
    REAL(dp), INTENT(inout) :: pxtm1(kbdim,klev,ktrac)   ! tracer mass/number mixing ratio
    REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ktrac)   ! tracer mass/number mixing ratio tendency
    
    !--- subroutine variables
    REAL(dp) :: zdpg(kbdim,klev)            ! air mass auxiliary factor (for diagnostics)
    
    REAL(dp) :: zdayfac(kbdim)              ! Increase oxidant concentrations (monthly averages) by day
    REAL(dp) :: znitefac(kbdim)             ! and night

    REAL(dp) :: znvte(kbdim,klev,nsoaspec)      ! tendency for non-volatile species (no gas tracer defined)
    
    REAL(dp) :: zvmrfac(kbdim,klev)         ! conversion factor VMR-> molecules cm-3
    
    REAL(dp) :: zoh(kbdim,klev)             ! OH concentration ### at t+dt
    REAL(dp) :: zo3(kbdim,klev)             ! O3 concentration ### at t+dt
    REAL(dp) :: zno3(kbdim,klev)            ! NO3 concentration ### at t+dt
    REAL(dp) :: zprecursor                  ! precursor concentration at t+dt
    ! Reaction rates [cm**3 molec-1 s-1]:
    REAL(dp) :: zrr_oh                      ! reaction rate precursor - OH
    REAL(dp) :: zrr_o3                      ! reaction rate precursor - O3
    REAL(dp) :: zrr_no3                     ! reaction rate precursor - NO3
    ! Loss/Formation rates [s-1]:
    REAL(dp) :: zrate_oh                    ! k(oh)*[OH]
    REAL(dp) :: zrate_o3                    ! k(o3)*[O3]
    REAL(dp) :: zrate_no3                   ! k(no3)*[NO3]
    REAL(dp) :: zratesum                    ! sum of 2 or more of above
    
    REAL(dp) :: zlossfac                    ! fractional precursor loss
    REAL(dp) :: zlossrate                   ! precursor loss rate [mmr s-1]
    REAL(dp) :: zprodfac                    ! SOA production factor
    REAL(dp) :: zprodfac1                   ! SOA production factor product 1
    REAL(dp) :: zprodfac2                   ! SOA production factor product 2
    REAL(dp) :: zprodrate1                  ! SOA production rate product 1
    REAL(dp) :: zprodrate2                  ! SOA production rate product 2
    REAL(dp) :: zeps                        ! EPSILON for type dp

    !>>csld #404
    REAL(dp) :: zlossfacsum  ! sum of fractional precursor losses (in part. for all different monoterpene species)
    REAL(dp) :: zlossratesum ! sum of precursor loss rate [mmr s-1] (in part. for all different monoterpene species)
    !<<csld #404
   
    INTEGER  :: idt_mt, idt_mt1, idt_mt2, & ! tracer indices
                idt_is, idt_is1, idt_is2, & 
                idt_tol, idt_xyl, idt_benz

    INTEGER  :: ix_tolpr, ix_xylpr, ix_benzpr      ! indices for the oxidation products of
                                                   ! toluene, xylene, benzene in the array znvte

    !>>csld #404
    INTEGER :: zarr_idt_mt(7)  ! array containing the indices of different monoterpenes 
    INTEGER :: jidt            ! loop counter inside the above array (zarr_idt_mt)
    !<<csld #404  

    INTEGER  :: jk, jl                      ! loop counters 
    INTEGER  :: ef_type
  
    !--- executable procedure

    !      Comment:
    !      Prescribed oxidant concentrations are used, we do not calculate changes
    !      in oxidant concentrations

    zeps = EPSILON(1._dp)
    znvte(:,:,:) = 0._dp

    !--- calculate day/night length 
    DO jl=1,kproma
        
       IF (daylength(jl,krow) > zeps) THEN
          zdayfac(jl) = 1._dp/daylength(jl,krow)
       ELSE
          zdayfac(jl) = 0._dp
       END IF

       IF ((1._dp-daylength(jl,krow)) > zeps) THEN
          znitefac(jl) = 1._dp/(1._dp-daylength(jl,krow))
       ELSE
          znitefac(jl) = 0._dp
       END IF

    END DO

    !---conversion factors
    DO jk=1,klev
       DO jl=1,kproma
          !---layer thickness*density (for diagnostics)
          zdpg(jl,jk) = (paph(jl,jk+1) - paph(jl,jk)) / grav

          !---for conversion of oxidant MMR to molecules cm-3
          zvmrfac(jl,jk) = 1.E-6_dp*avo*pap(jl,jk)/(argas*pt(jl,jk)) * amd

       END DO
    END DO

    !---SOA chemistry parameterisation

    !---precursor loss (with and without SOA production) and SOA production are calculated
    !   Reactions are (case with SOA production)
    !   PRE + OX -> (alpha1) P1 + (alpha2) P2 where P1, P2 are gas-phase SOA 
    !
    !   The reaction rate for a single reaction is:
    !   d[PRE]/dt = -k[PRE][OX]     where k=reaction rate
    !   since [OX] is prescribed (fixed) this has the simple solution
    !   [PRE] (final) = [PRE](initial)*EXP(-k[OX]t)    where t=time  
    ! 
    !   If the same precursor reacts with several oxidants (e.g. monoterpenes with O3, OH and NO3) we have
    !   d[PRE]/dt = SUM(i) [ -k(i)[PRE][OX(i)]
    !   leading to
    !   [PRE] (final) = [PRE](initial)*EXP( (sum(i) (-k(i)[OX(i)]t) ) ) = [PRE](initial)*(prod(i)( EXP(-k(i)[OX(i)]t) ) )
    !
    !   Whereas for a given SOA species we must only take into account only the reaction that produces that species,
    !   d[P1]/dt = alpha1 * 
    !   The method below calculates the loss factor EXP(-k[OX]t) for each oxidant (for the total
    !   precursor loss) and separately for precursor loss that leads to SOA production (for the 
    !   calculation of SOA formation)
    !
    !   OH and O3 are considered as daytime reactions, while NO3 is considered for nighttime.
    !   SOA yield from O3 oxidation differes between day and night due to different secondary 
    !   reaction chains.
    
    ! monoterpenes  
    !>>csld #404
    ! array containing the indices of the monoterpenes species
    zarr_idt_mt = (/ speclist(id_apin)%idt,  speclist(id_tbeta)%idt, &
                     speclist(id_bpin)%idt,  speclist(id_lim)%idt,   &
                     speclist(id_sab)%idt,   speclist(id_myrc)%idt,  &
                     speclist(id_car)%idt /)            
    !<<csld #404
    idt_mt1 = speclist(soaprop(1)%spid_tot)%idt
    idt_mt2 = speclist(soaprop(2)%spid_tot)%idt
    
    ! isoprene
    idt_is = speclist(id_isop)%idt
    idt_is1 = speclist(soaprop(3)%spid_tot)%idt
    idt_is2 = speclist(soaprop(4)%spid_tot)%idt
    
    ! toluene
    idt_tol = speclist(id_tol)%idt
    
    ! xylene
    idt_xyl = speclist(id_xyl)%idt
    
    ! benzene
    idt_benz = speclist(id_benz)%idt
    
    !!mgs!!  ! OH
    !!mgs!!  idt_oh  = speclist(id_oh)%idt
    !!mgs!!  ! O3
    !!mgs!!  idt_o3  = speclist(id_o3)%idt
    !!mgs!!  ! NO3
    !!mgs!!  idt_no3  = speclist(id_no3)%idt

    !---anthropogenics: subject to lumping
    !   indices are set to distinct values in the case of no lumping
    !   and to the same value in case of lumping.
    !   this avoids branching in the gridpoint loop

    SELECT CASE (nsoalumping) 
    CASE(0)                   ! no lumping all anthropogenics distinct
       ix_tolpr = 5
       ix_xylpr = 6
       ix_benzpr = 7
       
    CASE(1)                 ! lumping of anthropogenics into one SOA species
       ix_tolpr = 5
       ix_xylpr = ix_tolpr
       ix_benzpr = ix_tolpr

    CASE(2)                 ! lumping of anthropogenics into OC
       ix_tolpr = id_oc
       ix_xylpr = ix_tolpr
       ix_benzpr = ix_tolpr
    END SELECT

    ! obtain oxidant concentrations and convert from mass mixing ratio to concentration
    !Note: will use either climatologies or interactive fields from chemistry module
    CALL bc_apply(ibc_oh,  kproma, krow, zoh)
    CALL bc_apply(ibc_o3,  kproma, krow, zo3)
    CALL bc_apply(ibc_no3, kproma, krow, zno3)

    DO jk = 1,klev
       DO jl = 1,kproma
          zoh(jl,jk)  = zoh(jl,jk) * zvmrfac(jl,jk) / speclist(id_oh)%moleweight
          zo3(jl,jk)  = zo3(jl,jk) * zvmrfac(jl,jk) / speclist(id_o3)%moleweight
          zno3(jl,jk) = zno3(jl,jk) * zvmrfac(jl,jk) / speclist(id_no3)%moleweight
       END DO
       
       ! Why per level ?
       !-- apply daylength scaling if values are from climatology
       CALL bc_query(ibc_oh, ef_type=ef_type)
       IF (ef_type == EF_FILE) zoh(1:kproma,jk)  = zoh(1:kproma,jk) * zdayfac(1:kproma)

       CALL bc_query(ibc_o3, ef_type=ef_type)
       IF (ef_type == EF_FILE) zo3(1:kproma,jk)  = zo3(1:kproma,jk) * zdayfac(1:kproma)

       !-- apply lenght of night factor to no3
       CALL bc_query(ibc_no3, ef_type=ef_type)
       IF (ef_type == EF_FILE) zno3(1:kproma,jk)  = zno3(1:kproma,jk) * znitefac(1:kproma)
    END DO
    
    ! day time reactions
  
    DO jk = 1,klev
       DO jl = 1,kproma

          IF ((NINT(rdayl_x(jl,krow)) == 1)) THEN
!!mgs!!           !---convert oxidant species from VMR to molecules cm-3
!!mgs!!           zoh = zdayfac(jl) * zvmrfac(jl,jk) * &
!!mgs!!                 (pxtm1(jl,jk,idt_oh)+pxtte(jl,jk,idt_oh)*time_step_len)
!!mgs!!           zo3 = zdayfac(jl) * zvmrfac(jl,jk) * &
!!mgs!!                 (pxtm1(jl,jk,idt_o3)+pxtte(jl,jk,idt_o3)*time_step_len)

             !---monoterpenes

             !---reaction rates k=Aexp(B/T)...rate data from IUPAC
             zrr_oh = 1.2E-11_dp*EXP(440._dp/pt(jl,jk))
             zrr_o3 = 6.3E-16_dp*EXP(-580._dp/pt(jl,jk))

             zrate_oh = zrr_oh * zoh(jl,jk)
             zrate_o3 = zrr_o3 * zo3(jl,jk)

             zratesum = zrate_oh + zrate_o3

             !>> csld #404
             zlossratesum = 0._dp
             zlossfacsum = 0._dp  

             DO jidt = 1, size(zarr_idt_mt)
                ! indice of the monoterpene specie
                idt_mt = zarr_idt_mt(jidt)
                
                zprecursor = pxtm1(jl,jk,idt_mt) + pxtte(jl,jk,idt_mt)*time_step_len 
                zlossfac   = zprecursor * (1._dp - EXP(-zratesum*time_step_len))           
                zlossrate  = zlossfac / time_step_len

                pxtte(jl,jk,idt_mt) = pxtte(jl,jk,idt_mt) - zlossrate
               
                !sum over all monoterpenes
                zlossfacsum  = zlossfacsum  + zlossfac
                zlossratesum = zlossratesum + zlossrate 
             END DO
             !<<csld $404

             IF (zratesum > 1.E-20_dp) THEN
                zprodfac = (zrate_o3/zratesum) * zlossfacsum !csld #404 zlossfac -> zlosfacsum
             ELSE
                zprodfac = 0._dp
             END IF

             !---temperature dependent stoichiometric coefficients (Saathoff et al ACPD 2008)
             zprodfac1 = zprodfac * (0.715_dp - 0.002_dp*pt(jl,jk))
             zprodrate1 = zprodfac1 / time_step_len
             pxtte(jl,jk,idt_mt1) = pxtte(jl,jk,idt_mt1) + zprodrate1
             
             zprodfac2 = zprodfac * 1200._dp*EXP(-pt(jl,jk)/35._dp)
             zprodrate2 = zprodfac2 / time_step_len
             pxtte(jl,jk,idt_mt2) = pxtte(jl,jk,idt_mt2) + zprodrate2

             d_prod_soa_mterp(jl,krow) = d_prod_soa_mterp(jl,krow) + &
                                       (zprodrate1+zprodrate2) * delta_time * zdpg(jl,jk)
             d_chem_sink_mterp(jl,krow) = d_chem_sink_mterp(jl,krow) + &
                                          zlossratesum * delta_time * zdpg(jl,jk) !csld #404 zlossfac -> zlosfacsum

             !---isoprene
           
             !---reaction rates k=Aexp(B/T)...rate data from IUPAC
             zrr_oh = 2.7E-11_dp*EXP(390._dp/pt(jl,jk))
             zrr_o3 = 1.03E-14_dp*EXP(-1995._dp/pt(jl,jk))

             zrate_oh = zrr_oh * zoh(jl,jk)
             zrate_o3 = zrr_o3 * zo3(jl,jk)

             zratesum = zrate_oh + zrate_o3

             zprecursor = pxtm1(jl,jk,idt_is)+pxtte(jl,jk,idt_is)*time_step_len
             zlossfac = zprecursor * (1._dp - EXP(-zratesum*time_step_len))           
             zlossrate = zlossfac / time_step_len

             pxtte(jl,jk,idt_is) = pxtte(jl,jk,idt_is) - zlossrate

             IF (zratesum > 1.E-20_dp) THEN
                zprodfac = (zrate_oh/zratesum) * zlossfac
             ELSE
                zprodfac = 0._dp
             END IF

             !---isoprene 2-product data Henze and Seinfeld GRL 2007
             zprodfac1 = 0.232_dp * zprodfac 
             zprodrate1 = zprodfac1 / time_step_len
             pxtte(jl,jk,idt_is1) = pxtte(jl,jk,idt_is1) + zprodrate1

             zprodfac2 = 0.0288_dp * zprodfac
             zprodrate2 = zprodfac2 / time_step_len
             pxtte(jl,jk,idt_is2) = pxtte(jl,jk,idt_is2) + zprodrate2
             
             d_prod_soa_isop(jl,krow) = d_prod_soa_isop(jl,krow) + &
                                        (zprodrate1+zprodrate2) * delta_time * zdpg(jl,jk)

             d_chem_sink_isop(jl,krow) = d_chem_sink_isop(jl,krow) + &
                                         zlossrate * delta_time * zdpg(jl,jk)

             !---toluene           
             zrr_oh = 1.81E-12_dp*EXP(338._dp/pt(jl,jk))
             zrate_oh = zrr_oh * zoh(jl,jk)

             zprecursor = pxtm1(jl,jk,idt_tol)+pxtte(jl,jk,idt_tol)*time_step_len
             zlossfac = zprecursor * (1._dp - EXP(-zrate_oh*time_step_len))           
             zlossrate = zlossfac / time_step_len

             pxtte(jl,jk,idt_tol) = pxtte(jl,jk,idt_tol) - zlossrate

             zprodfac = zlossfac

             zprodfac1 = 0.36_dp * zprodfac 
             zprodrate1 = zprodfac1 / time_step_len
             znvte(jl,jk,ix_tolpr) = zprodrate1

             d_prod_soa_tol(jl,krow) = d_prod_soa_tol(jl,krow) + &
                                       zprodrate1 * delta_time * zdpg(jl,jk)
 
             d_chem_sink_tol(jl,krow) = d_chem_sink_tol(jl,krow) + &
                                        zlossrate * delta_time * zdpg(jl,jk)

             !---xylene
             zrr_oh = 2.31E-11_dp
             zrate_oh = zrr_oh * zoh(jl,jk)

             zprecursor = pxtm1(jl,jk,idt_xyl)+pxtte(jl,jk,idt_xyl)*time_step_len
             zlossfac = zprecursor * (1._dp - EXP(-zrate_oh*time_step_len))           
             zlossrate = zlossfac / time_step_len

             pxtte(jl,jk,idt_xyl) = pxtte(jl,jk,idt_xyl) - zlossrate

             zprodfac = zlossfac

             zprodfac1 = 0.30_dp * zprodfac 
             zprodrate1 = zprodfac1 / time_step_len
             znvte(jl,jk,ix_xylpr) =  znvte(jl,jk,ix_xylpr)+zprodrate1

             d_prod_soa_xyl(jl,krow) = d_prod_soa_xyl(jl,krow) + &
                                       zprodrate1 * delta_time * zdpg(jl,jk)

             d_chem_sink_xyl(jl,krow) = d_chem_sink_xyl(jl,krow) + &
                                        zlossrate * delta_time * zdpg(jl,jk)

             !---benzene           
             zrr_oh = 2.33E-12_dp*EXP(-193._dp/pt(jl,jk))
             zrate_oh = zrr_oh *zoh(jl,jk)

             zprecursor = pxtm1(jl,jk,idt_benz)+pxtte(jl,jk,idt_benz)*time_step_len
             zlossfac = zprecursor * (1._dp - EXP(-zrate_oh*time_step_len))           
             zlossrate = zlossfac / time_step_len

             pxtte(jl,jk,idt_benz) = pxtte(jl,jk,idt_benz) - zlossrate

             zprodfac = zlossfac

             zprodfac1 = 0.37_dp * zprodfac 
             zprodrate1 = zprodfac1 / time_step_len
             znvte(jl,jk,ix_benzpr) = znvte(jl,jk,ix_benzpr)+zprodrate1

             d_prod_soa_benz(jl,krow) = d_prod_soa_benz(jl,krow) + &
                                        zprodrate1 * delta_time * zdpg(jl,jk)

             d_chem_sink_benz(jl,krow) = d_chem_sink_benz(jl,krow) + &
                                         zlossrate * delta_time * zdpg(jl,jk)

          ELSE               ! night time reactions
             !---convert oxidant species from VMR to molecules cm-3
             !!mgs!!           zno3 = znitefac(jl) * zvmrfac(jl,jk) * &
             !!mgs!!                  (pxtm1(jl,jk,idt_no3)+pxtte(jl,jk,idt_no3)*time_step_len)

             !---monoterpenes
             zrr_no3 = 1.2E-12_dp*EXP(490._dp/pt(jl,jk))
             zrate_no3 = zrr_no3 * zno3(jl,jk)
             
             !>> csld #404
             zlossratesum = 0._dp

             DO jidt = 1, size(zarr_idt_mt)     
                ! indice of the monoterpene specie
                idt_mt = zarr_idt_mt(jidt)
            
                zprecursor = pxtm1(jl,jk,idt_mt) + pxtte(jl,jk,idt_mt)*time_step_len 
                zlossfac   =  zprecursor * (1._dp - EXP(-zrate_no3*time_step_len))
                zlossrate  = zlossfac / time_step_len
                
                pxtte(jl,jk,idt_mt) = pxtte(jl,jk,idt_mt) - zlossrate  

                !sum over all monoterpenes
                zlossratesum = zlossratesum + zlossrate 
             END DO
             !<<csld $404

            d_chem_sink_mterp(jl,krow) = d_chem_sink_mterp(jl,krow) + &
                                         zlossratesum * delta_time * zdpg(jl,jk) !csld #404 zlossfac -> zlosfacsum

             !---isoprene
             zrr_no3 = 3.15E-12_dp*EXP(-450._dp/pt(jl,jk))
             zrate_no3 = zrr_no3 * zno3(jl,jk)

             zprecursor = pxtm1(jl,jk,idt_is)+pxtte(jl,jk,idt_is)*time_step_len
             zlossfac = zprecursor * (1._dp - EXP(-zrate_no3*time_step_len))           
             zlossrate = zlossfac / time_step_len

             pxtte(jl,jk,idt_is) = pxtte(jl,jk,idt_is) - zlossrate

             d_chem_sink_isop(jl,krow) = d_chem_sink_isop(jl,krow) + &
                                         zlossrate * delta_time * zdpg(jl,jk)

             !---xylene
             zrr_no3 = 2.6E-16_dp
             zrate_no3 = zrr_no3 * zno3(jl,jk)

             zprecursor = pxtm1(jl,jk,idt_xyl)+pxtte(jl,jk,idt_xyl)*time_step_len
             zlossfac = zprecursor * (1._dp - EXP(-zrate_no3*time_step_len))           
             zlossrate = zlossfac / time_step_len

             pxtte(jl,jk,idt_xyl) = pxtte(jl,jk,idt_xyl) - zlossrate

             d_chem_sink_xyl(jl,krow) = d_chem_sink_xyl(jl,krow) + &
                                        zlossrate * delta_time * zdpg(jl,jk)

          END IF
       END DO
    END DO

    !---split the newly-formed non-volatile SOA among the modes
    ! CALL soa_progn(kproma, kbdim, klev, krow, ktrac, pxtm1, pxtte)
    CALL soa_equi1(kproma, kbdim, klev, krow, ktrac, pxtm1, pxtte, znvte)
    

  END SUBROUTINE soa2prod


  SUBROUTINE soa_equi0(kproma, kbdim, klev, krow,  ktrac, &
                       pap,    pt ,   pq,   pxtm1, pxtte  )

    ! This subroutine calculates the equilibrium SOA gas and aerosol masses 
    ! 
    ! In order to minimise the number of prognostic tracers in the model, for each volatile SOA species
    ! only a single prognostic tracer is used. That tracer contains the total mass 
    ! of the species (gas + all applicable aerosol modes). This subroutine calculates the equilibrium
    ! gas-aerosol phase partitioning according to the 2-product model of SOA formation and 
    ! the distribution of aerosol phase mass between the different modes. Diagnostic (not transported) 
    ! tracers hold the resulting gas and per-mode aerosol masses. These can subsequently be used
    ! in radiation and aerosol microphysical calculations.
    !
    ! Author
    ! ------
    ! Declan O'Donnell MPI-M 2007

    USE mo_kind,         ONLY: dp
    USE mo_time_control, ONLY: time_step_len
    USE mo_physical_constants, ONLY: rd, rv
    USE mo_species,      ONLY: speclist
    USE mo_ham_species,  ONLY: id_so4, id_oc
    USE mo_ham,          ONLY: aerocomp, nsoaspec, sizeclass, nclass
    USE mo_ham_soa,      ONLY: soaprop, lso4_in_m0 

    IMPLICIT NONE

    INTEGER,  INTENT(in)    :: kproma                    ! geographic block number of locations
    INTEGER,  INTENT(in)    :: kbdim                     ! geographic block maximum number of locations 
    INTEGER,  INTENT(in)    :: klev                      ! numer of levels 
    INTEGER,  INTENT(in)    :: ktrac                     ! number of tracers
    INTEGER,  INTENT(in)    :: krow                      ! geographic block number
    REAL(dp), INTENT(in)    :: pap(kbdim,klev)           ! air pressure 
    REAL(dp), INTENT(in)    :: pt(kbdim,klev)            ! air temperature 
    REAL(dp), INTENT(in)    :: pq(kbdim,klev)            ! specific humidity
    REAL(dp), INTENT(inout) :: pxtm1(kbdim,klev,ktrac)   ! tracer mass/number mixing ratio
    REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ktrac)   ! tracer mass/number mixing ratio tendency

    !--- subroutine variables
    REAL(dp) :: zmnvmode(kbdim,klev,nclass)          ! SOA absorbing non-volatile mass per mode
    REAL(dp) :: zmnv(kbdim,klev)                   ! total non-volatiles
    REAL(dp) :: zsoamass(kbdim,klev,nsoaspec)          ! total (gas+aerosol) mass
    REAL(dp) :: zm0(kbdim,klev)                    ! equilibrium SOA absorbing mass
    REAL(dp) :: zkp(kbdim,klev,nsoaspec)               ! partitioning coefficient per SOA compound
    REAL(dp) :: zmmr2ug(kbdim,klev)                ! conversion factor mmr -> ug m-3
    
    REAL(dp) :: zgasmass(kbdim,klev,nsoaspec)
    REAL(dp) :: zaeromass

    REAL(dp) :: zrmoist                            ! gas constant of moist air
    REAL(dp) :: ztmst                              ! model time step

    INTEGER :: ixso4, ixoc
    INTEGER :: jk, jl, jn, jm, jt, jspec, jaero 

    !--- executable procedure

    !---initialisations
    ztmst = time_step_len
    zmnv(:,:) = 0._dp
    zmnvmode(:,:,:) = 0._dp
    zsoamass(:,:,:) = 0._dp

    !---factor to convert mmr to ug m-3 = 10e9*air density
    
    DO jk=1,klev
       DO jl=1,kproma
          zrmoist = rd*(1._dp - pq(jl,jk)) + rv*pq(jl,jk)
          zmmr2ug(jl,jk) = 1.E9_dp*pap(jl,jk) / (pt(jl,jk)*zrmoist)
       END DO
    END DO

    !---total non-volatile SOA absorbing mass

    !---total non-volatile SOA absorbing mass
    DO jn=1,nclass
       IF (sizeclass(jn)%lsoainclass) THEN

          !---start with OC
          ixoc = speclist(id_oc)%iaerocomp(jn)       !!mgs!! im7table(jn,id_oc)
          IF (ixoc > 0) THEN
             jt = aerocomp(ixoc)%idt
             zmnvmode(1:kproma,:,jn) = zmnvmode(1:kproma,:,jn) + &
                                       (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*time_step_len)
          END IF

          !---add the non-volatile SOA
          DO jm=1,nsoaspec
             IF (.NOT. soaprop(jm)%lvolatile) THEN
                jspec = soaprop(jm)%spid_soa
                jaero = speclist(jspec)%iaerocomp(jn)
                jt = aerocomp(jaero)%idt
                zmnvmode(1:kproma,:,jn) = zmnvmode(1:kproma,:,jn) + &
                                          (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*time_step_len)
             END IF
          END DO
        
          !---sulphate if included according to model switch
          IF (lso4_in_m0) THEN
             ixso4 = speclist(id_so4)%iaerocomp(jn)       !!mgs!! im7table(jn,id_so4)
             IF (ixso4 > 0) THEN
                jt = aerocomp(ixso4)%idt
                zmnvmode(1:kproma,:,jn) = zmnvmode(1:kproma,:,jn) + &
                                          (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*time_step_len)
             END IF
          END IF

          zmnvmode(1:kproma,:,jn) = MAX(zmnvmode(1:kproma,:,jn),0._dp)
          zmnv(1:kproma,:) = zmnv(1:kproma,:) + zmnvmode(1:kproma,:,jn)
        
       END IF
     
    END DO

    !---convert to ug m-3
    zmnv(1:kproma,:) = zmnv(1:kproma,:) * zmmr2ug(1:kproma,:)
    DO jn=1,nclass
       zmnvmode(1:kproma,:,jn) = zmnvmode(1:kproma,:,jn) * zmmr2ug(1:kproma,:)
    END DO

    !---now compute the total SOA mass per volatile species (aerosol+gas)
    DO jm=1,nsoaspec
       IF (soaprop(jm)%lvolatile) THEN

          ! prognostic tracer
          jt = speclist(soaprop(jm)%spid_tot)%idt 

          zsoamass(1:kproma,:,jm) = zmmr2ug(1:kproma,:)*&
                                    (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       END IF
    END DO
  
    !---calculation of total SOA absorbing mass M0
    !
    CALL soa_part(kproma, kbdim, klev, krow, pt, zmnv, zsoamass, zm0, zkp) 

    DO jm=1,nsoaspec
       IF (soaprop(jm)%lvolatile) THEN
          !---convert total soa mass back to mmr
          zsoamass(1:kproma,:,jm) = zsoamass(1:kproma,:,jm)/zmmr2ug(1:kproma,:)

          ! gas phase mass G(i) = S(i) / (1+Kp*M0) where S(i) = total (gas+aerosol) SOA mass
          jt = speclist(soaprop(jm)%spid_soa)%idt

          zgasmass(1:kproma,:,jm) = zsoamass(1:kproma,:,jm) / &
                                    (1._dp+zkp(1:kproma,:,jm)*zm0(1:kproma,:))
          pxtm1(1:kproma,:,jt) = MAX(zgasmass(1:kproma,:,jm),0._dp)
          pxtte(1:kproma,:,jt) = 0._dp
       END IF
    END DO

    !---aerosol mass is partitioned in proportion to the non-volatile content of each mode
    DO jn=1,nclass
       IF (sizeclass(jn)%lsoainclass) THEN
          DO jm=1,nsoaspec

             IF (soaprop(jm)%lvolatile) THEN
                jspec = soaprop(jm)%spid_soa
                jaero = speclist(jspec)%iaerocomp(jn)
                jt = aerocomp(jaero)%idt
                
                DO jk=1,klev
                   DO jl=1,kproma
                      IF (zmnv(jl,jk) > 0._dp) THEN              
                         zaeromass = zsoamass(jl,jk,jm) - zgasmass(jl,jk,jm)
                         ! zaeromass = MAX(zaeromass, 0._dp)

                         pxtm1(jl,jk,jt) = zaeromass*zmnvmode(jl,jk,jn)/zmnv(jl,jk)
                         pxtte(jl,jk,jt) = 0._dp
                      END IF
                   END DO
                END DO            !---end do gridpoints
             END IF            !---end if volatile
             
          END DO            !---end do soa
       END IF            !---end if soa in mode
    END DO            !---end do modes


  END SUBROUTINE soa_equi0

  
  SUBROUTINE soa_equi1(kproma, kbdim, klev, krow, ktrac, pxtm1, pxtte, pnvte)

    ! This subroutine calculates the equilibrium SOA gas and aerosol masses 
    ! 
    ! In order to minimise the number of prognostic tracers in the model, for each volatile SOA species
    ! only a single prognostic tracer is used. That tracer contains the total mass 
    ! of the species (gas + all applicable aerosol modes). This subroutine calculates the equilibrium
    ! gas-aerosol phase partitioning according to the 2-product model of SOA formation and 
    ! the distribution of aerosol phase mass between the different modes. Diagnostic (not transported) 
    ! tracers hold the resulting gas and per-mode aerosol masses. These can subsequently be used
    ! in radiation and aerosol microphysical calculations.
    !
    ! Author
    ! ------
    ! Declan O'Donnell MPI-M 2007
    
    USE mo_kind,         ONLY: dp
    USE mo_time_control, ONLY: time_step_len
    USE mo_ham,          ONLY: aerocomp, nsoaspec, sizeclass, nclass
    USE mo_ham_soa,      ONLY: soaprop, lso4_in_m0
    USE mo_species,      ONLY: speclist
    USE mo_ham_species,  ONLY: id_oc, id_so4
    
    IMPLICIT NONE
    
    INTEGER,  INTENT(in)     :: kproma                    ! geographic block number of locations
    INTEGER,  INTENT(in)     :: kbdim                     ! geographic block maximum number of locations 
    INTEGER,  INTENT(in)     :: klev                      ! numer of levels 
    INTEGER,  INTENT(in)     :: ktrac                     ! number of tracers
    INTEGER,  INTENT(in)     :: krow                      ! geographic block number
    REAL(dp), INTENT(inout)  :: pxtm1(kbdim,klev,ktrac)   ! tracer mass/number mixing ratio
    REAL(dp), INTENT(inout)  :: pxtte(kbdim,klev,ktrac)   ! tracer tendency 
    REAL(dp), INTENT(in)     :: pnvte(kbdim,klev,nsoaspec)
    
    !--- subroutine variables
    REAL(dp) :: zmnvmode(kbdim,klev,nclass)                 ! SOA absorbing non-volatile mass per mode
    REAL(dp) :: zmnv(kbdim,klev)                          ! total non-volatiles
    
    INTEGER :: jn, jm, jspec, jaero, jt                   ! mode, SOA, species, species-in-mode and tracer indices

    INTEGER :: ixoc, ixso4
    
    !--- executable procedure

    !---initialisations
    zmnv(:,:) = 0._dp
    zmnvmode(:,:,:) = 0._dp

   !---total non-volatile SOA absorbing mass
   DO jn=1,nclass
      IF (sizeclass(jn)%lsoainclass) THEN

         !---start with OC
         ixoc = speclist(id_oc)%iaerocomp(jn)       !!mgs!! im7table(jn,id_oc)
         IF (ixoc > 0) THEN
            jt = aerocomp(ixoc)%idt
            zmnvmode(1:kproma,:,jn) = zmnvmode(1:kproma,:,jn) + &
                                      (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*time_step_len)
         END IF

         !---add the non-volatile SOA
         DO jm=1,nsoaspec
            IF (.NOT. soaprop(jm)%lvolatile) THEN
               jspec = soaprop(jm)%spid_soa
               jaero = speclist(jspec)%iaerocomp(jn)
               jt = aerocomp(jaero)%idt
               zmnvmode(1:kproma,:,jn) = zmnvmode(1:kproma,:,jn) + &
                                         (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*time_step_len)
            END IF
         END DO
        
         !---sulphate if included according to model switch
         IF (lso4_in_m0) THEN
            ixso4 = speclist(id_so4)%iaerocomp(jn)       !!mgs!! im7table(jn,id_so4)
            IF (ixso4 > 0) THEN
               jt = aerocomp(ixso4)%idt
               zmnvmode(1:kproma,:,jn) = zmnvmode(1:kproma,:,jn) + &
                                         (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*time_step_len)
            END IF
         END IF

         zmnvmode(1:kproma,:,jn) = MAX(zmnvmode(1:kproma,:,jn),0._dp)
         zmnv(1:kproma,:) = zmnv(1:kproma,:) + zmnvmode(1:kproma,:,jn)
        
      END IF
     
   END DO

    !---newly formed non-volatile SOA is partitioned according to the mass fraction of 
    !   existing organic mass in each mode
    DO jn = 1,nclass
       IF (sizeclass(jn)%lsoainclass) THEN
          DO jm = 1,nsoaspec
             IF (.NOT. soaprop(jm)%lvolatile) THEN
                jspec = soaprop(jm)%spid_soa
                jaero = speclist(jspec)%iaerocomp(jn)
                jt = aerocomp(jaero)%idt

                WHERE (zmnv(1:kproma,:) > 0._dp)
                   pxtte(1:kproma,:,jt) = pxtte(1:kproma,:,jt) + pnvte(1:kproma,:,jm) * &
                                          zmnvmode(1:kproma,:,jn)/zmnv(1:kproma,:)
                END WHERE
             END IF
          END DO
       END IF
    END DO

  END SUBROUTINE soa_equi1


  SUBROUTINE soa_part(kproma, kbdim, klev, krow, pt, pmnv, psoa, pm0, pkp) 

    !----------------------------------- SOA gas-aerosol partitioning -------------------------------------
    ! Author:
    ! ------
    ! Declan O'Donnell, MPI-M, 2007
    !
    
    ! soa_part handles the partitioning of secondary organics between the gas and aerosol phases
    
    USE mo_kind,         ONLY: dp
    USE mo_physical_constants, ONLY: argas
    USE mo_ham,          ONLY: nsoaspec
    USE mo_ham_soa,      ONLY: soaprop
    USE mo_ham_m7ctl,    ONLY: cmin_aerml
    
    IMPLICIT NONE

    !--- Subroutine interface
    INTEGER, INTENT(IN)     :: kproma, kbdim, klev, krow             ! Grid parameters for this processor
    REAL(dp), INTENT(IN)    :: pt(kbdim,klev)                        ! Air temperature
    REAL(dp), INTENT(IN)    :: pmnv(kbdim,klev)                      ! Non-volatile SOA absorbing mass concentration 
    REAL(dp), INTENT(IN)    :: psoa(kbdim,klev,nsoaspec)                 ! Total (gas+aerosol) SOA species concentration
    REAL(dp), INTENT(OUT)   :: pm0(kbdim,klev)                       ! SOA absorbing mass
    REAL(dp), INTENT(OUT)   :: pkp(kbdim,klev,nsoaspec)                  ! SOA partitioning coefficient
    !--- Local variables
    
    !--- parameters
    INTEGER, PARAMETER :: nmax = 30
    REAL(dp), PARAMETER :: zminmass = 1.E-25_dp
    
    !--- subroutine variables
    REAL(dp) :: ztol                              ! Tolerance in iterative solution
    REAL(dp) :: zt0  
    REAL(dp) :: zk0
    REAL(dp) :: zk0t0
    REAL(dp) :: zdh_fac                           ! for calculation of Kp_om
    REAL(dp) :: zexpdht0
    REAL(dp) :: zkpsoa(kbdim,klev,nsoaspec)           ! Kom,i*[SOA]i                
    REAL(dp) :: zx1(kbdim,klev)                   ! work areas for iterative solution
    REAL(dp) :: zx2(kbdim,klev)                   ! 
    REAL(dp) :: zx3(kbdim,klev)                   ! 
    REAL(dp) :: zy1(kbdim,klev)                   ! 
    REAL(dp) :: zy2(kbdim,klev)                   ! 
    
    LOGICAL :: lconverged(kbdim,klev)             ! Iterative solution found or not
    INTEGER :: nunsolved                          ! Number of gridboxes where we have not yet solved for M0
    
    INTEGER :: i, jk, jl, jm                      ! loop counters
  
    !--- executable procedure 

    DO jm = 1,nsoaspec

       IF (soaprop(jm)%lvolatile) THEN

          !---Calculation of temperature-dependent partitioning coefficient Kp,om
          !
          !                 T        [ dH  ( 1    1  ) ]
          !  Kp,om = Kref *---- * EXP[ ---*( - - --- ) ]   
          !                Tref      [  R  ( T   Tref) ]
          
          zt0 = soaprop(jm)%tref
          zk0 = soaprop(jm)%Kp
          zdh_fac = soaprop(jm)%dH/argas
          zk0t0 = zk0/zt0
          zexpdht0 = EXP(-zdh_fac/zt0)
          
          DO jk=1,klev
             DO jl=1,kproma
              
                pkp(jl,jk,jm) = zk0t0 * pt(jl,jk) * zexpdht0 * EXP(zdh_fac/pt(jl,jk)) 
                zkpsoa(jl,jk,jm) = pkp(jl,jk,jm) * psoa(jl,jk,jm)
             END DO
          END DO
           
       END IF

    END DO

  
    !--- 3.3 Iterative calculation of partition
    !
    !       Find M0 such that 
    !                     -- Kom,i * [SOA]i  
    !       M0 = Mnv + M0*\  --------------
    !                     /  1 + Kom,i*M0
    !                     --
    !       
    !       solved iteratively with secant method - 
    !       method chosen so that it can run on a vector machine
    !       without stopping the show
    
    !---initialisations
    zx1(1:kproma,:) = pmnv(1:kproma,:)
    zx2(1:kproma,:) = pmnv(1:kproma,:)
    zy1(:,:) = 0._dp
    zy2(:,:) = 0._dp
    pm0(:,:) = 0._dp
    nunsolved = 0
    lconverged(:,:) = .FALSE.
    ztol = cmin_aerml*100._dp
  
    !---initial guess: maximum mass
    DO jm=1,nsoaspec
       IF (soaprop(jm)%lvolatile) THEN
          DO jk=1,klev
             DO jl=1,kproma
                !---start with zx2=all mass in aerosol phase
                zx2(jl,jk) = zx2(jl,jk)+psoa(jl,jk,jm)
                
                !---evaluate y2 = f(x2) and y1 = f(x1)
                !   start with sum term
                zy1(jl,jk) = zy1(jl,jk) +  &
                             psoa(jl,jk,jm) * pkp(jl,jk,jm) / &
                             (1._dp+pkp(jl,jk,jm)*zx1(jl,jk))
                   
                zy2(jl,jk) = zy2(jl,jk) +  &
                             psoa(jl,jk,jm) * pkp(jl,jk,jm) / &
                             (1._dp+pkp(jl,jk,jm)*zx2(jl,jk))
             END DO
          END DO
       END IF
    END DO

    DO jk=1,klev
       DO jl=1,kproma
        
          IF (zx2(jl,jk) < zminmass) THEN          ! Very small total mass, nothing to do
             lconverged(jl,jk) = .TRUE.            ! Mark the gridbox solved
          ELSE                                     ! else we need to solve for M0
             zy1(jl,jk) = pmnv(jl,jk)-zx1(jl,jk)+zx1(jl,jk)*zy1(jl,jk)
             zy2(jl,jk) = pmnv(jl,jk)-zx2(jl,jk)+zx2(jl,jk)*zy2(jl,jk)
        
             !---evaluate if either y1 or y2 is a solution
             IF (ABS(zy1(jl,jk)) < ztol) THEN
                lconverged(jl,jk) = .TRUE.
                pm0(jl,jk) = zx1(jl,jk)
             ELSEIF ((ABS(zy2(jl,jk)) < ztol)) THEN
                lconverged(jl,jk) = .TRUE.
                pm0(jl,jk) = zx2(jl,jk)
             ELSEIF (ABS(zy2(jl,jk)-zy1(jl,jk)) < ztol) THEN
                lconverged(jl,jk) = .TRUE.           ! Where no difference between y2 and y1
                pm0(jl,jk) = zx2(jl,jk)              ! arbitrary choice of x2 over x1
             ELSE
                nunsolved = nunsolved + 1                             
             END IF

          END IF
       END DO
    END DO

    iterative_solution: &
       
    DO i=1,nmax
    WHERE ( .NOT. lconverged(1:kproma,:) )                 ! Find next x according to secant method
       zx3(1:kproma,:) = zy2(1:kproma,:)*(zx2(1:kproma,:)-zx1(1:kproma,:)) / &
                         (zy2(1:kproma,:)-zy1(1:kproma,:))
       zx1(1:kproma,:) = zx2(1:kproma,:)
       zy1(1:kproma,:) = zy2(1:kproma,:)
       zx2(1:kproma,:) = zx2(1:kproma,:) - zx3(1:kproma,:)
       
       ! Now evaluate y2 = f(x2). Before starting loop over SOA, initialise y2 = Mnv-x2
       zy2(1:kproma,:) = pmnv(1:kproma,:)-zx2(1:kproma,:)
    END WHERE

                        ! Now add x2*SUM(pkp*SOA[i]/(1+x2*pkp))       
    DO jm=1,nsoaspec
       IF (soaprop(jm)%lvolatile) THEN
          WHERE (.NOT.lconverged(1:kproma,:))
             zy2(1:kproma,:) = zy2(1:kproma,:)+zx2(1:kproma,:) * &
                               psoa(1:kproma,:,jm)*pkp(1:kproma,:,jm) / &
                              (1._dp+pkp(1:kproma,:,jm)*zx2(1:kproma,:))
          END WHERE
       END IF
    END DO

    DO jk = 1,klev               ! Test for convergence
       DO jl = 1,kproma
          IF (.NOT. lconverged(jl,jk)) THEN    

             IF ((ABS(zy2(jl,jk)) .LT. ztol).OR.(ABS(zy2(jl,jk)-zy1(jl,jk)).LT. ztol)) THEN
                lconverged(jl,jk) = .TRUE.         
                                           
                pm0(jl,jk) = zx2(jl,jk)
                nunsolved = nunsolved - 1
             END IF
          END IF
       END DO
    END DO
    ! Now check the number of gridboxes where convergence has been achieved

    IF (nunsolved == 0) EXIT iterative_solution
     
                                 
  END DO iterative_solution

  pm0(1:kproma,:) = MAX(0._dp,pm0(1:kproma,:))

 END SUBROUTINE soa_part


 SUBROUTINE soa_progn(kproma, kbdim, klev, krow, ktrac, pxtm1, pxtte)

   ! Purpose
   ! -------
   ! This subroutine calculates the tendency of total SOA mass per species. 
   ! In order to minimise CPU consumption, only the total SOA mass is transported
   ! (for semivolatile species). The gas and aerosol phase masses are diagnosed early
   ! in the gridpoint calculations (before radiation is calculated) in the subroutine
   ! soa_equi0. 
   ! Subsequent processes (dry deposition, wet deposition, aerosol microphysics...) 
   ! calculate the impact on the tracer tendencies of the gas and aerosol phase 
   ! of each species but leave the tendency of the total mass unchanged.
   ! This subroutine calculates the tendency (per species) of the total semivolatile mass
   ! at the end of each timestep.
   !
   ! Author
   ! ------
   ! Declan O'Donnell, MPI-MET Hamburg, 2008
   !
   
   !---inherited data, types and functions
   USE mo_kind,            ONLY: dp
   USE mo_species,         ONLY: speclist
   USE mo_ham,             ONLY: aerocomp, nsoaspec, sizeclass, nclass
   !!mgs!!  USE mo_ham_soa,         ONLY: soaprop, soaspecies, isoa_ix
   USE mo_ham_soa,         ONLY: soaprop 
   
   IMPLICIT NONE 

   !---subroutine interface
   INTEGER,  INTENT(in)    :: kproma                    ! geographic block number of locations
   INTEGER,  INTENT(in)    :: kbdim                     ! geographic block maximum number of locations 
   INTEGER,  INTENT(in)    :: klev                      ! numer of levels 
   INTEGER,  INTENT(in)    :: ktrac                     ! number of tracers
   INTEGER,  INTENT(in)    :: krow                      ! geographic block number
   REAL(dp), INTENT(inout) :: pxtm1(kbdim,klev,ktrac)   ! tracer mass/number mixing ratio
   REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ktrac)   ! tracer mass/number mixing ratio tendency
   
   !---local data
   INTEGER :: jn, jm, jt1, jt2, jt3, jspec, jaero 
   
   !---executable procedure
   
   DO jm=1,nsoaspec
      IF (soaprop(jm)%lvolatile) THEN

         !---index to diagnostic (not transported) tracer for SOA gases
         jt1 = speclist(soaprop(jm)%spid_soa)%idt

         !---index to prognostic (transported) tracer for total SOA mass
         jt2 = speclist(soaprop(jm)%spid_tot)%idt

         !---update total mass tendency from gas phase tendency
         pxtte(1:kproma,:,jt2) = pxtte(1:kproma,:,jt2) + pxtte(1:kproma,:,jt1)

         DO jn=1,nclass
            IF (sizeclass(jn)%lsoainclass) THEN
               jspec = soaprop(jm)%spid_soa
               jaero = speclist(jspec)%iaerocomp(jn)
               
               !---index to diagnostic (not transported) aerosol phase tracer per applicable mode
               jt3 = aerocomp(jaero)%idt

               !---update total mass tendency from aerosol phase tendency
               pxtte(1:kproma,:,jt2) = pxtte(1:kproma,:,jt2) + pxtte(1:kproma,:,jt3)
            END IF
         END DO
              
      END IF
   END DO

 END SUBROUTINE soa_progn

END MODULE mo_ham_soa_processes
