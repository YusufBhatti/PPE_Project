!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_vbs_partition.f90
!!
!! \brief
!! mo_ham_vbs contains parameters, switches and initialization routines 
!! for the volatility basis set 'add-on'
!!
!! \responsible_coder
!! Thomas Kühn (thk), thomas.h.kuhn@uef.fi
!!
!! \revision_history
!!   -# T. Kühn (UEF) - initial setup
!!   -# J. Merikanto (FMI) - aqueous phase SOA code
!! 
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
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

MODULE mo_ham_vbs_partition

  USE mo_kind, ONLY: &
       dp

  USE mo_ham_vbsctl, ONLY: &
       nvbs_setup,    &
       t_voc_prec,         &
       t_vbs_group,        &
       vbs_nvocs,     &
       vbs_voc_prec,  &
       vbs_ngroup,    &
       vbs_set,       &
       t_aq_soa,           &   
       laqsoa,             &
       aqsoa_ngroup,       &
       aqsoa_set


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
       vbs_gas_phase_chem,&         ! computes gas phase oxidation of VOCs
       vbs_condensation_salsa, &       ! computes condensation/evaporation of
       aqsoa_condensation_salsa             ! sVOCs onto/off particles for SALSA   


CONTAINS

  SUBROUTINE vbs_gas_phase_chem(&
       kproma, kbdim, klev,  klevp1, krow, ktrac, paph, pap, pt, pxtm1, pxtte&
       )

    ! -----------------------------------------------------------------------
    !
    ! SUBROUTINE vbs_gas_phase_chem
    !
    ! Computes the production of VBS species due to oxidation of VOC 
    ! precursors in the gas phase. Considered oxidants are OH, O3 and NO3
    !  
    ! Authors:
    ! --------
    ! Declan O'Donnell, MPI-Met
    ! Kai Zhang, MPI-Met
    ! Thomas Kühn, UEF    6/2015 --
    ! Joonas Merikanto, FMI 7/2015
    !
    ! vbs_gas_phase_chem is called 
    ! from *physc_subm_2* in *mo_submodel_interface*
    !
    ! -----------------------------------------------------------------------



    ! use statements
    USE mo_physical_constants, ONLY: avo, argas, grav, amd 
    USE mo_time_control, ONLY: time_step_len
    USE mo_boundary_condition, ONLY: bc_apply, bc_query
    USE mo_ham_streams,  ONLY: daylength
    USE mo_geoloc,       ONLY: rdayl_x
    USE mo_external_field_processor, ONLY: EF_FILE
    USE mo_species,      ONLY: speclist
    USE mo_ham,          ONLY: ibc_oh, ibc_o3, ibc_no3
    USE mo_ham_species,  ONLY: id_oh, id_o3, id_oc, id_no3


    ! input/output parameters
    INTEGER,  INTENT(in)    :: kproma                  ! geographic block number of locations
    INTEGER,  INTENT(in)    :: kbdim                   ! geographic block max number of locations 
    INTEGER,  INTENT(in)    :: klev                    ! numer of levels 
    INTEGER,  INTENT(in)    :: klevp1                  ! numer of levels + 1  
    INTEGER,  INTENT(in)    :: ktrac                   ! number of tracers
    INTEGER,  INTENT(in)    :: krow                    ! geographic block number
    REAL(dp), INTENT(in)    :: paph(kbdim,klevp1)      ! air pressure at layer interface 
    REAL(dp), INTENT(in)    :: pap(kbdim,klev)         ! air pressure at layer center    
    REAL(dp), INTENT(in)    :: pt (kbdim,klev)         ! air temperature 
    REAL(dp), INTENT(in)    :: pxtm1(kbdim,klev,ktrac) ! tracer mass/number mixing ratio
    REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ktrac) ! tracer mass/number mixing ratio tendency

    
    ! thk: change to ECHAM naming convention!
    ! local variables
    TYPE(t_voc_prec), POINTER :: zprec ! local copy for easier reading
    TYPE(t_vbs_group), POINTER :: zbase ! local copy for easier reading
    TYPE(t_aq_soa), POINTER :: zaqsoa ! local copy for easier reading

    ! conversion factors
    REAL(dp) ::     &
         zmmr_to_c(kbdim, klev) ! conversion factor from mass mixing ratio to #/cm3
         
    ! Gas phase reactant concentrations
    REAL(dp) ::     &
         zc_oh(kbdim, klev),  & ! OH concentration at t+dt
         zc_o3(kbdim, klev),  & ! O3 concentration at t+dt
         zc_no3(kbdim, klev), & ! NO3 concentration at t+dt
         zc_prec,             & ! VOC precursor concentation at t+dt
         zc_aqsoa_prec          ! aqsoa concentation at t+dt
    
    ! Gas phase reaction rates
    REAL(dp) ::     &
         zr_rate_oh_prec,   & ! voc precursor reaction rate with OH [m3/(mol*s)]
         zc_rate_oh_prec,   & ! voc precursor conversion rate due to OH [1/s]
         zr_rate_o3_prec,   & ! voc precursor reaction rate with O3 [m3/(mol*s)]
         zc_rate_o3_prec,   & ! voc precursor conversion rate due to O3 [1/s]
         zr_rate_no3_prec,  & ! voc precursor reaction rate with NO3 [m3/(mol*s)]
         zc_rate_no3_prec,  & ! voc precursor conversion rate due to NO3 [1/s]
         zc_rate_tot_prec,  & ! voc precursor total conversion rate
         zr_rate_oh_aqsoa,   & ! aqsoa reaction rate with OH [m3/(mol*s)]
         zc_rate_oh_aqsoa,   & ! aqsoa conversion rate due to OH [1/s]
         zr_rate_o3_aqsoa,   & ! aqsoa reaction rate with O3 [m3/(mol*s)]
         zc_rate_o3_aqsoa,   & ! aqsoa conversion rate due to O3 [1/s]
         zr_rate_no3_aqsoa,  & ! aqsoa reaction rate with NO3 [m3/(mol*s)]
         zc_rate_no3_aqsoa,  & ! aqsoa conversion rate due to NO3 [1/s]
         zc_rate_tot_aqsoa     ! aqsoa total conversion rate
    
    ! zloss/zproduction rates
    REAL(dp) :: &
         zloss_prec,           & ! total zloss of precursor VOC
         zloss_aqsoa,          & ! total zloss of aqsoa
         zlossrate_prec,       & ! zloss rate of precursor VOC
         zlossrate_aqsoa,      & ! zloss rate of aqsoa
         zprodrate_vbs,        & ! zproduction rate for vbs
         zprodrate_aqsoa,      & ! zproduction rate for aqsoa
         zlossoh_aqsoa,        & ! OH zloss of aqsoa
         zlosso3_aqsoa,        & ! O3 zloss of aqsoa
         zlossno3_aqsoa,       & ! NO3 zloss of aqsoa
         zlossphoto_aqsoa        ! photodissociation zloss of aqsoa

    !scaling factors
    REAL(dp) ::     &
         zdayfac(kbdim),      & ! day length scaling factor for climatoloty
         znightfac(kbdim),    & ! night length scaling factor for climatoloty
         zdayfac_oh(kbdim),   & ! day length scaling factor for climatoloty -- OH
         zdayfac_o3(kbdim),   & ! day length scaling factor for climatoloty -- O3
         znightfac_no3(kbdim)   ! night length scaling factor for climatoloty -- NO3
    

    !small number
    REAL(dp) :: &
         zeps

    INTEGER :: jk, jl, jb, jv
    INTEGER :: ef_type
    INTEGER :: idt_prec, idt_base, idt_aqsoa
    INTEGER :: zday_switch
    
    ! -----------------------------------------------------------------------
    ! executable procedure
    ! -----------------------------------------------------------------------
    

    ! >> from soa_processes:
    ! thk: change the day/night treatment to how it's done 
    !      in mo_ham_chemistry ?
    !--- calculate day/night length 
    zeps = EPSILON(1._dp)
    DO jl=1,kproma
        
       IF (daylength(jl,krow) > zeps) THEN
          zdayfac(jl) = 1._dp/daylength(jl,krow)
       ELSE
          zdayfac(jl) = 0._dp
       END IF

       IF ((1._dp-daylength(jl,krow)) > zeps) THEN
          znightfac(jl) = 1._dp/(1._dp-daylength(jl,krow))
       ELSE
          znightfac(jl) = 0._dp
       END IF

    END DO


!    ! >> from mo_ham_chemistry (not sure about this yet)
!    WHERE(daylength(1:kproma,krow) > zeps)
!       ! First we account for daylength, so that the integral of the OH
!       ! concentration over a month gives the monthly mean (this introduces an
!       ! artificial variation of the OH concentration over the month):
!       zdayfac(1:kproma) = 1._dp/daylength(1:kproma,krow)
!       ! Solar local time (as fraction of the day):
!       slt(1:kproma) = mod(1.15741E-5_dp*(isecond + philon_2d(1:kproma,krow)*240._dp),1._dp)
!       ! Then we account for the diurnal cycle so that the OH concentration
!       ! follows a cosine peak between sunrise and sunset:
!       zdayfac(1:kproma) = zdayfac(1:kproma)*0.5_dp*pi*cos(pi*(slt(1:kproma)-0.5_dp)*zdayfac(1:kproma))
!       !>>dod correction of negative zdayfac values:
!       !      check if solar local time is within half a day of local noon 
!       ldaylight(1:kproma) = (slt(1:kproma) > (0.5_dp-0.5_dp*daylength(1:kproma,krow)) ) .AND. &
!            (slt(1:kproma) < (0.5_dp+0.5_dp*daylength(1:kproma,krow)) )
!    ELSEWHERE
!       zdayfac(1:kproma)= 0._dp
!       ldaylight(1:kproma) = .FALSE.
!    END WHERE



    ! scaling is only needed if the oxidant concentrations come from a
    ! climatology. If they are computed online, we set the scaling factor
    ! to one, i.e. we don't scale

    ! OH
    CALL bc_query(ibc_oh, ef_type=ef_type)
    IF (ef_type == EF_FILE) THEN
       zdayfac_oh(1:kproma) = zdayfac(1:kproma)
    ELSE 
       zdayfac_oh(1:kproma) = 1.0_dp
    END IF
       
    ! O3
    CALL bc_query(ibc_o3, ef_type=ef_type)
    IF (ef_type == EF_FILE) THEN
       zdayfac_o3(1:kproma) = zdayfac(1:kproma)
    ELSE 
       zdayfac_o3(1:kproma) = 1.0_dp
    END IF

    ! NO3
    !-- apply length of night factor to no3
    CALL bc_query(ibc_no3, ef_type=ef_type)
    IF (ef_type == EF_FILE) THEN
       znightfac_no3(1:kproma) = znightfac(1:kproma)
    ELSE 
       znightfac_no3(1:kproma) = 1.0_dp
    END IF



    !---conversion factors
    DO jk=1,klev
       DO jl=1,kproma
          !---for conversion of oxidant MMR to molecules cm-3
          ! note that both amd and speclist%moleweight below are in [g/mol]
          zmmr_to_c(jl,jk) = 1.0e-6_dp*avo*pap(jl,jk)/(argas*pt(jl,jk))*amd

       END DO
    END DO



    ! obtain oxidant concentrations and convert from mass mixing ratio to concentration
    ! Note: will use either climatologies or interactive fields from chemistry module
    ! obtaining oxidant mass mixing ratios at t+dt
    CALL bc_apply(ibc_oh,  kproma, krow, zc_oh)
    CALL bc_apply(ibc_o3,  kproma, krow, zc_o3)
    CALL bc_apply(ibc_no3, kproma, krow, zc_no3)

    ! converting mass mixing ratios to molecule/cm3
    ! and scale for day/night time length, if values come from climatology
    DO jk = 1,klev
       DO jl = 1,kproma
          zc_oh(jl,jk) =zc_oh(jl,jk) *zmmr_to_c(jl,jk)/speclist(id_oh)%moleweight *zdayfac_oh(jl)
          zc_o3(jl,jk) =zc_o3(jl,jk) *zmmr_to_c(jl,jk)/speclist(id_o3)%moleweight *zdayfac_o3(jl)
          zc_no3(jl,jk)=zc_no3(jl,jk)*zmmr_to_c(jl,jk)/speclist(id_no3)%moleweight*znightfac_no3(jl)
       END DO
    END DO
    ! << from soa_processes
    
    
    DO jv = 1,vbs_nvocs

       ! marking the precursor we currently process
       zprec => vbs_voc_prec(jv)
       idt_prec = speclist(zprec%spid)%idt

       DO jk = 1,klev
          DO jl = 1,kproma

             ! the precursor concentration in this grid cell
             zc_prec = pxtm1(jl,jk,idt_prec)+pxtte(jl,jk,idt_prec)*time_step_len 
             
          
             IF ((NINT(rdayl_x(jl,krow)) == 1)) THEN
 
              !---reaction rates k=k_0*exp(E_p/T); E_p=E/R
              zr_rate_oh_prec  = zprec%k_0_OH  * EXP(zprec%Eact_p_OH/pt(jl,jk))!*zday_switch
              zr_rate_o3_prec  = zprec%k_0_O3  * EXP(zprec%Eact_p_O3/pt(jl,jk))!*zday_switch
              zr_rate_no3_prec = zprec%k_0_NO3 * EXP(zprec%Eact_p_NO3/pt(jl,jk))!*(1-zday_switch)
 
              !---conversion rates k_x*[x]
              zc_rate_oh_prec  = zr_rate_oh_prec  * zc_oh(jl,jk)
              zc_rate_o3_prec  = zr_rate_o3_prec  * zc_o3(jl,jk)
              zc_rate_no3_prec = zr_rate_no3_prec * zc_no3(jl,jk)
 
              !---total conversion rate
              zc_rate_tot_prec = zc_rate_oh_prec + zc_rate_o3_prec + zc_rate_no3_prec
              
              !---amount of reacted voc
              zloss_prec = zc_prec * (1._dp - EXP(-zc_rate_tot_prec*time_step_len))           
              
              !---and the rate of zloss
              zlossrate_prec = zloss_prec / time_step_len
              
              !---correcting the tendency
              pxtte(jl,jk,idt_prec) = pxtte(jl,jk,idt_prec) - zlossrate_prec
 
                 
              !---calculating the production of products
              DO jb = 1,vbs_ngroup

                  ! marking the product we currently process
                 zbase => vbs_set(jb)
                 idt_base = speclist(zbase%spid)%idt
                 !idt_base = vbs_set(jb)%idt
                
                 ! the production rate for jb
                 zprodrate_vbs = zlossrate_prec * zprec%stoich_coeff(jb)
 
                 !---correcting the tendency
                 pxtte(jl,jk,idt_base) = pxtte(jl,jk,idt_base) + zprodrate_vbs
  
              END DO !vbs_ngroup

              if(laqsoa) then
         
                 DO jb = 1 , aqsoa_ngroup  
                    zaqsoa => aqsoa_set(jb)                           
                    idt_aqsoa = speclist(zaqsoa%spid)%idt

                    ! the zproduction rate for jb
                    zprodrate_aqsoa = zc_prec*(1._dp-EXP(-zc_rate_oh_prec*time_step_len))/time_step_len  &
                                * zprec%stoich_coeff(jb+vbs_ngroup)
                    if(jb.eq.2.and.jv.eq.2) zprodrate_aqsoa =  zprodrate_aqsoa +  & !Hard coded for ISOP+O3 -> 0.01 GLYX
                               zc_prec*(1._dp-EXP(-zc_rate_o3_prec*time_step_len))/time_step_len & 
                                * zprec%stoich_coeff(jb+vbs_ngroup)*0.5_dp
                    if(jb.eq.2) zprodrate_aqsoa = zprodrate_aqsoa + 0.24*zlossoh_aqsoa/time_step_len !hard coded for EIPOX+OH->0.24*Glyoxal   

                    !---correcting the tendency
                    pxtte(jl,jk,idt_aqsoa) = pxtte(jl,jk,idt_aqsoa)+ zprodrate_aqsoa

                    ! the gas phase loss rates for aqSOA                         
                    !concentration                  
                    zc_aqsoa_prec=pxtm1(jl,jk,idt_aqsoa)+pxtte(jl,jk,idt_aqsoa)*time_step_len                   
                    !---reaction rates k=k_0*exp(E_p/T); E_p=E/R
                    zr_rate_oh_aqsoa   = zaqsoa%k_0_OH  * EXP(zaqsoa%Eact_p_OH /pt(jl,jk))!*zday_switch
                    zr_rate_o3_aqsoa   = zaqsoa%k_0_O3  * EXP(zaqsoa%Eact_p_O3 /pt(jl,jk))!*zday_switch
                    zr_rate_no3_aqsoa  = zaqsoa%k_0_NO3 * EXP(zaqsoa%Eact_p_NO3/pt(jl,jk))!*zday_switch
                    
                    zc_rate_oh_aqsoa  = zr_rate_oh_aqsoa  * zc_oh(jl,jk)
                    zc_rate_o3_aqsoa  = zr_rate_o3_aqsoa  * zc_o3(jl,jk)
                    zc_rate_no3_aqsoa = zr_rate_no3_aqsoa * zc_no3(jl,jk)
                
                    zlossoh_aqsoa    = zc_aqsoa_prec * (1._dp - EXP(-zc_rate_oh_aqsoa*time_step_len))   
                    zlosso3_aqsoa    = zc_aqsoa_prec * (1._dp - EXP(-zc_rate_o3_aqsoa*time_step_len)) 
                    zlossno3_aqsoa   = zc_aqsoa_prec * (1._dp - EXP(-zc_rate_no3_aqsoa*time_step_len)) 
                    zlossphoto_aqsoa = zc_aqsoa_prec * (1._dp - EXP(-zaqsoa%photodis *time_step_len/2.0)) !2.0*zday_switch (now average over whole day)
                    zloss_aqsoa = zlossoh_aqsoa + zlosso3_aqsoa + zlossno3_aqsoa + zlossphoto_aqsoa
                    
                    zlossrate_aqsoa = zloss_aqsoa / time_step_len

                    pxtte(jl,jk,idt_aqsoa) = pxtte(jl,jk,idt_aqsoa) - zlossrate_aqsoa                 
                
                 END DO
              END IF
           END IF
        END DO !kproma
     END DO !klev
  END DO !vbs_nvocs               



  END SUBROUTINE vbs_gas_phase_chem


  SUBROUTINE vbs_condensation_salsa(&
       ptemp,                               &
       pt_step,                             &
       pc_gas,                              &
       pnaero, pvols,                       &
       pd_wet, pcolrate                     &
       )

    ! -----------------------------------------------------------------------
    !
    ! SUBROUTINE vbs_condensation_salsa
    !
    ! Computes the exchange of volatility basis set (VBS) groups between 
    ! gas- and particle phases and their distribution among the aerosol 
    ! size bins. Changes are only calculated for one grid cell, looping 
    ! over grid and height happens in the calling routine ( *condensation*
    ! in *mo_ham_salsa_dynamics* ).
    ! 
    ! Authors:
    ! --------
    ! Thomas Kühn, UEF    6/2015 --
    !
    ! vbs_condensation_salsa is called 
    ! from *condensation* in *mo_ham_salsa_dynamics*
    !
    ! Related literature:
    ! -Jacobson (2005), Fundamentals of Atmospheric Modelling, 2nd Edition
    ! -Farina et. al (2010), Modeling global secondary organic aerosol
    !  formation and processing with the volatility basis set: Implications
    !  for anthropogenic secondary organic aerosol
    !
    !
    ! -----------------------------------------------------------------------

    USE mo_math_constants, ONLY: &
         pi_6

    USE mo_physical_constants, ONLY: &
         avo

    USE mo_species, ONLY: &
         speclist              ! list of chemical species in ham

    USE mo_ham_species, ONLY: &
         id_wat

    USE mo_ham_vbsctl, ONLY: &
         t_vbs_group,        & ! the VBS data structure
         vbs_ngroup,         & ! number of VBS bins
         vbs_set,            & ! VBS
         nclass_vbs            ! number of modes/bins that include VBS soa

    USE mo_ham_salsactl, ONLY: &
         surfw0                ! surface tension of water [J/m2]
         
    USE mo_ham, ONLY : &
         sizeclass,          & ! aerosol bin properties
         subm_aerospec,      & ! aerosol species list
         subm_ngasspec,      & ! number of gas species in HAM
         subm_naerospec,     & ! aerosol species list
         subm_naerospec_nowat,     &
         subm_aerospec_nowat,     &
         nclass                ! number of bins

    ! -----------------------------------------------------------------------

    ! input / output variables

    REAL(dp), INTENT(IN) :: &  
         pd_wet(nclass),                       & ! wet diameter of particles in each bin [m]
         pcolrate(nclass),                     & ! collision rate per bin
         ptemp,                                & ! ambient temperature [K]
         pt_step                                 ! timestep [s]

    REAL(dp), INTENT(INOUT) :: &
         pc_gas(subm_ngasspec),                & ! Gas phase tracer concentrations [molec/m3]
         pnaero(nclass),                       & ! aerosol number concentrations [#/m3]
         pvols(nclass,subm_naerospec_nowat)                ! aerosol volume concentrations per species [m3/m3]


    ! -----------------------------------------------------------------------

    ! local variables

    ! indexing VBS groups
    TYPE(t_vbs_group), POINTER :: &
         zgroup                                  ! local copy for easier reading

    ! intermediates for calculation
    REAL(dp) :: &
         zeps = 1.e-30_dp,                       & ! small number
         zlarge = 1e30_dp,                       & ! big number
         ztemp_sum,                              & ! temporary storage
         zv_tot,                                 & ! total dry aerosol volume concentraion
         zc_part_all(nclass_vbs, subm_naerospec),  & ! particle phase concentrations per species
                                                   ! for the vbs relevant bins + water [m3/m3]
                                                   ! note the different dimensions!
         zk_mass(nclass_vbs),                    & ! mass transfer coefficient
         zk_mass_tot,                            & ! mass transfer coefficient
         zs_prime(vbs_ngroup,nclass_vbs),        & ! equilibrium saturation ratio
         zc_sat(vbs_ngroup),                     & ! uncorrected saturation vapor concentration
         zc_sat_part(nclass_vbs),                & ! mole-fraction corrected sat. vap. conc. per class
         zc_gas(vbs_ngroup),                     & ! current gas phase concentration [mol/m3]
         zc_gas_new(vbs_ngroup),                 & ! new gas phase concentration [mol/m3]
         zc_part(vbs_ngroup,nclass_vbs),         & ! currect particle phase concentration [mol/m3]
         zc_part_new(vbs_ngroup,nclass_vbs),     & ! new particle phase concentration [mol/m3]
         zc_tot(vbs_ngroup)                        ! total concentration [mol/m3]



    ! we divide the time step into increasing 
    ! sub-intervals, with 
    ! sum(dt_n) = pt_step and
    ! dt_(n+1)=dt_n * ztime_fac
    INTEGER, PARAMETER :: zn_steps = 10          ! amount of sub-time steps
    REAL(dp),PARAMETER :: ztime_fac = 1.5_dp     ! sub-time step growth factor
    REAL(dp) ::           zh                     ! current sub-time step length

    ! loop indices
    INTEGER :: & 
         jt,                                   & ! time sub-step
         jc,                                   & ! class
         js,                                   & ! species
         jvc,                                  & ! vbs class
         jg                                      ! vbs group

    ! -----------------------------------------------------------------------
    ! executable procedure
    ! -----------------------------------------------------------------------


    ! Calculating the particle phase concentrations for all aerosol species plus 
    ! water for the VBS-relevant bins. 
    jvc = 0
    DO jc = 1,nclass
       
       ! treat only bins that are relevant for VBS:
       IF (sizeclass(jc)%lsoainclass) THEN

          ! increment vbs class counter
          jvc = jvc + 1

          ! total volume concentration
          zv_tot = 0.0_dp
          
          ! All aerosol species are taken from pvols:
          DO js = 1,subm_naerospec_nowat

             ! speclist%moleweight is in g/mol
             zc_part_all(jvc,js) = pvols(jc,js)*speclist(subm_aerospec_nowat(js))%density&
                  /speclist(subm_aerospec_nowat(js))%moleweight*1e3 
    
             ! adding to zv_tot
             zv_tot = zv_tot + pvols(jc,js)
          END DO ! js


          ! aersol water (from wet diameter, number concentration and dry aersol volume concentration)
          ! mwa is in kg/mol
          zc_part_all(jvc, subm_naerospec) =&
               ((pi_6*pd_wet(jc)**3)*pnaero(jc)-zv_tot)*speclist(id_wat)%density/&
               (speclist(id_wat)%moleweight/1000.)


          ! mass transfer coefficient per aerosol bin Jacobson (16.64) 
          zk_mass(jvc) = pcolrate(jc)

          ! equilibrium saturation ratio per aerosol bin Jacobson (16.47)
          DO jg = 1,vbs_ngroup
             ! easier access
             zgroup => vbs_set(jg)

             ! Kelvin effect only
             zs_prime(jg,jvc) = exp(4.0_dp*surfw0*zgroup%mv/(avo*ptemp*pd_wet(jc)))
          END DO

       END IF ! lsoainclass

    END DO ! jc


    ! total mass transfer coefficent
    zk_mass_tot = sum(zk_mass(:))


    ! mapping to VBS 
    DO jg = 1,vbs_ngroup

       ! easier access
       zgroup => vbs_set(jg)

       ! computing the uncorrected saturation vapor concentration
       ! after Farina et al (2010), Equation (4)
       zc_sat(jg) = zgroup%C0*(zgroup%T0/ptemp)*&
            exp(zgroup%Hvap_eff*(1._dp/zgroup%T0-1._dp/ptemp))

       ! copyping the particle volume concentrations
       DO jvc = 1, nclass_vbs
          zc_part(jg, jvc)=zc_part_all(jvc,zgroup%id_vols)
       END DO ! jc

       ! converting gas phase concentration
       zc_gas(jg) = pc_gas(zgroup%id_gasspec)/avo

       ! total concentration
       zc_tot(jg) = zc_gas(jg)+sum(zc_part(jg,:))


    END DO ! jg


    ! initial copy
    zc_gas_new(:)    = zc_gas(:)
    zc_part_new(:,:) = zc_part(:,:)

    ! Starting the sub-time step iteration
    zh=pt_step*(ztime_fac-1._dp)/(ztime_fac**real(zn_steps,dp)-1._dp)
    
    DO jt = 1, zn_steps ! loop over sub-time steps
    
       ! the particle phase concentrations [mol/m3] and
       ! the gas phase concentrations [mol/m3] before condensation:
       DO jg = 1,vbs_ngroup
          
          ! if concentrations are to small, we skip
          IF (zc_tot(jg) < zeps) CYCLE            
       
          ! easier access
          zgroup => vbs_set(jg)
       
          ! computing the new mole-fraction weighted, per-bin equilibrium vapor concentration
          DO jvc = 1,nclass_vbs
             IF (sum(zc_part_all(jvc,:)) > 0.0_dp) THEN
                zc_sat_part(jvc) = zc_sat(jg)*zc_part(jg,jvc)/sum(zc_part_all(jvc,:))
             ELSE
                zc_sat_part(jvc) = zc_sat(jg)
             END IF
          END DO ! jvc



          ! retrieve previous gas and particle phase concentration
          zc_gas(jg)    = zc_gas_new(jg)
          zc_part(jg,:) = zc_part_new(jg,:)
             
          ! computing the new gas phase concentration (16.71)
          zc_gas_new(jg) = min(&
               (zc_gas(jg)+zh*sum(zk_mass(:)*zs_prime(jg,:)*zc_sat_part(:)))/(1.0_dp+zh*zk_mass_tot),&
               zc_tot(jg)&
               )

          ! computing the new particle phase concentrations and correcting for negative
          ! concentrations (16.69) + (16.72a)
          DO jvc = 1, nclass_vbs
             zc_part_new(jg,jvc)=max(&
                  zc_part(jg,jvc)+zh*zk_mass(jvc)*(zc_gas_new(jg)-zs_prime(jg,jvc)*zc_sat_part(jvc)),&
                  0._dp&
                  )
          END DO ! jvc

          ! correcting for too high particle concentrations
          ztemp_sum = sum(zc_part_new(jg,:))
          IF (ztemp_sum < zeps) THEN
             zc_part_new(jg,:) = 0.0_dp
             zc_gas_new(jg)    = zc_tot(jg)
          ELSE
             zc_part_new(jg,:) = zc_part_new(jg,:)*(zc_tot(jg)-zc_gas_new(jg))/ztemp_sum
          END IF

          ! mapping back from VBS
          zc_part_all(:,zgroup%id_vols) = zc_part_new(jg, :)

       END DO ! jg


       ! increasing the time sub-step size
       zh = zh*ztime_fac

    END DO ! jt


    ! mapping back from VBS
    DO jg = 1,vbs_ngroup

       ! easier access
       zgroup => vbs_set(jg)

       ! converting gas phase concentration
       pc_gas(zgroup%id_gasspec) = zc_gas_new(jg)*avo

    END DO ! jg



    ! Calculating the new species volume fractions from concentrations
    jvc = 0
    DO jc = 1,nclass

       ! only VBS-relevant bins
       IF (sizeclass(jc)%lsoainclass) THEN

          !increment vbs class counter
          jvc = jvc+1
          
          DO jg = 1,vbs_ngroup
             zgroup => vbs_set(jg)
             ! speclist%moleweight is in g/mol
             pvols(jc,zgroup%id_vols) = zc_part_all(jvc,zgroup%id_vols)&
                  /speclist(zgroup%spid)%density&
                  *speclist(zgroup%spid)%moleweight*1e-3
          END DO ! js
          
       END IF ! lsoainclass
    END DO ! jc
          


  END SUBROUTINE vbs_condensation_salsa



  SUBROUTINE aqsoa_condensation_salsa(   &
    ptemp,                               & ! ambient conditions
    pc_gas,                              & ! gas phase concentrations
    pnaero,                              & ! aerosol number concentration
    pvols,                               & ! aerosol properties
    pd_wet                               &
    )

    ! -----------------------------------------------------------------------
    !
    ! SUBROUTINE aqsoa_condensation_salsa
    !
    ! Computes the exchange of aqsoa precursors between 
    ! gas- and particle phases and their distribution among the aerosol 
    ! size bins.
    ! 
    ! In order to keep changes to the original SALSA routines to a minimum, 
    ! we retrieve the gas phase concentrations for the aqsoa groups directly
    ! from pxtm1 and pxtte.
    !  
    ! Authors:
    ! --------
    ! Thomas Kühn, UEF    6/2015 --
    ! Joonas Merikanto, FMI    6/2015 --
    !
    ! aqsoa_condensation_salsa is called 
    ! from *condensation* in *mo_ham_salsa_dynamics*
    !
    ! NOTE: Currently a very simple scheme based on the use of effective Henry's coefficient
    !       Can be updated later on
    !
    ! Related literature:
    ! -Nguyen, T. B., et al. "Organic aerosol formation from the reactive uptake of isoprene 
    ! epoxydiols (IEPOX) onto non-acidified inorganic seeds." Atmospheric Chemistry and Physics 
    ! 14.7 (2014): 3497-3510.
    !
    !
    ! -----------------------------------------------------------------------

    USE mo_physical_constants, ONLY: &
         avo

    USE mo_math_constants, ONLY: &
         pi_6

    USE mo_ham_vbsctl, ONLY: &
         t_aq_soa,        & ! the VBS data structure
         aqsoa_ngroup,    & ! number of VBS bins
         aqsoa_set,       &  ! a
         nclass_aqsoa        ! number of modes/bins that include wet soa
         
         
    USE mo_ham, ONLY : &
         sizeclass,          & ! aerosol bin properties
         subm_aerospec,      & ! aerosol species list
         subm_ngasspec,      & ! number of gas species in HAM
         subm_aerospec_nowat, &! all but water aero species list
         subm_naerospec_nowat,&! number of non water aero species
         subm_naerospec,      &! number of aero species
         nclass                ! number of bins

    USE mo_species, ONLY: &
         speclist              ! list of chemical species in ham

    USE mo_ham_species, ONLY: &
         id_wat

    ! -----------------------------------------------------------------------

    ! input / output variables

         
    REAL(dp), INTENT(IN) :: &  
         pd_wet(nclass),                       & ! wet diameter of particles in each bin [m]
         ptemp                                ! ambient temperature [K]

    REAL(dp), INTENT(INOUT) ::&
         pc_gas(subm_ngasspec),         & ! gas phase tracer concentrations [molec/m3]
         pnaero(nclass),                       & ! aerosol number concentrations [#/m3]
         pvols(nclass,subm_naerospec_nowat)   ! aerosol volume concentrations per species [m3/m3]

    ! -----------------------------------------------------------------------

    ! local variables

    ! indexing aqsoa groups
    TYPE(t_aq_soa), POINTER :: zgroup      ! local copy for easier reading

    ! intermediates for calculation
    REAL(dp) :: &
         zv_tot,                                 & ! total dry aerosol volume concentration
         zc_part_all(nclass_aqsoa, subm_naerospec),& ! particle phase concentrations per species
                                                   ! for the vbs relevant bins + water [mol/m3]
         zc_gas(aqsoa_ngroup),                    & ! current gas phase concentration [mol/m3]
         zc_gas_new(aqsoa_ngroup),                & ! new gas phase concentration [mol/m3]
         zc_part(aqsoa_ngroup,nclass_aqsoa),            & ! currect particle phase concentration [mol/m3]
         zc_part_new(aqsoa_ngroup,nclass_aqsoa),        & ! new particle phase concentration [mol/m3]
         zc_tot(aqsoa_ngroup)                       ! total concentration [m3/m3]

    REAL(dp) :: &
         zeps = 1.e-30_dp                         ! small number

    INTEGER :: & 
         jc,                                & ! bin (nclass) counter
         ja,                                & ! aqsoa_ngroup counter
         jac,                               & ! aqsoa class counter
         js                                   ! non-water aerosol counter

    REAL(dp) ::    p_lwc                      ! Total aerosol water [g/m3]
    REAL(dp) ::    p_lwc_bin                  ! Aerosol water in specific bin [g/m3]

    ! -----------------------------------------------------------------------
    ! executable procedure
    ! -----------------------------------------------------------------------

   ! experimenetally derived partitioning based on aerosol water content

   ! the particle phase concentrations [mol/m3] and
    ! the gas phase concentrations [mol/m3] before condensation:

    jac = 0
    DO jc = 1,nclass

       ! treat only bins that are relevant for VBS:
       IF (sizeclass(jc)%lsoainclass) THEN

          ! increment aqsoa class counter
          jac = jac + 1
          ! total volume concentration
          zv_tot = 0.0_dp

          ! All aerosol species are taken from pvols:
          DO js = 1,subm_naerospec_nowat
             zc_part_all(jac,js) = pvols(jc,js)*speclist(subm_aerospec_nowat(js))%density&
                  /speclist(subm_aerospec_nowat(js))%moleweight*1e3 
             ! adding to zv_tot [m3/m3]
             zv_tot = zv_tot + pvols(jc,js)
          END DO ! js           

          ! aersol water (from wet diameter, number concentration and dry aersol volume concentration)
          ! mwa is in kg/mol
          zc_part_all(jac, subm_naerospec) =&
               ((pi_6*pd_wet(jc)**3)*pnaero(jc)-zv_tot)*speclist(id_wat)%density/&
               (speclist(id_wat)%moleweight/1000.)

       END IF ! lsoainclass

    END DO ! jc

    ! Total relevant aerosol water in [g/m3] 
    p_lwc=max(0.0_dp,sum(zc_part_all(:, subm_naerospec))*speclist(id_wat)%moleweight)

    IF(p_lwc.gt.zeps) THEN
       DO ja = 1,aqsoa_ngroup
          zgroup => aqsoa_set(ja)
          
          DO jac = 1, nclass_aqsoa
             zc_part(ja,jac)=zc_part_all(jac,zgroup%id_aqsoa)
          END DO

          zc_gas(ja) = pc_gas(zgroup%id_gasspec)/avo
          
          ! Total aqsoa concentration (gas+liquid phase) [mol/m3]
          zc_tot(ja) = zc_gas(ja)+sum(zc_part(ja,:))
       
          !First calculate new gas phase concentration according to total aerosol water content and
          !"effective Henry's constant", ref. Nguyen et al. (2014)  
          zc_gas_new(ja)=zc_tot(ja)*(1.0_dp/(1.0_dp+zgroup%Eff_henry_aerosol*1.E-6*0.08205736*ptemp*p_lwc))
          
          !Update gas phase concentration
          pc_gas(zgroup%id_gasspec) = zc_gas_new(ja)*avo
     
          !Then distribute remaining aqsoa according to aerosol water content 
          DO jac = 1, nclass_aqsoa        
             p_lwc_bin=max(0.0_dp,zc_part_all(jac, subm_naerospec)*speclist(id_wat)%moleweight)  
             zc_part_new(ja,jac)=(zc_tot(ja)-zc_gas_new(ja))*p_lwc_bin/p_lwc  
          END DO

          ! Update liquid phase concentration
                 zc_part_all(:,zgroup%id_aqsoa) = zc_part_new(ja, :)

       END DO

       ! Calculating the new species volume fractions from concentrations
       jac = 0

       DO jc = 1, nclass
          ! only aqsoa-relevant bins
          IF (sizeclass(jc)%lsoainclass) THEN
             !increment aqsoa class counter
             jac = jac+1
             DO ja = 1,aqsoa_ngroup
                zgroup => aqsoa_set(ja)
                pvols(jc,zgroup%id_aqsoa)= zc_part_all(jac,zgroup%id_aqsoa)&
                     /speclist(zgroup%spid)%density&
                     *speclist(zgroup%spid)%moleweight*1e-3
             END DO
          END IF ! lsoainclass
       END DO
  END IF 

  END SUBROUTINE aqsoa_condensation_salsa


END MODULE mo_ham_vbs_partition
