!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_vbs.f90
!!
!! \brief
!! mo_ham_vbs_species makes sure that all volatile and semi-volatile 
!! compounds are set up properly
!!
!! \responsible_coder
!! Thomas Kuehn (thk), thomas.h.kuhn@uef.fi
!!
!! \revision_history
!!   -# Declan O'Donnell (MPI-Met) - soa code; used as guideline here
!!   -# Kai Zhang (MPI-Met) - modifications to soa code
!!   -# T. Kuehn (UEF) - original code (2015)
!!   -# J. Merikanto (FMI) - added aqueous phase SOA
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
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of 
!! any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_vbs
  USE mo_kind, ONLY: &
       dp

  USE mo_ham_vbsctl, ONLY: &
       nvbs_setup,         &
       vbs_setup_3class_OC,&
       vbs_setup_3class,   &
       laqsoa,             &
       t_voc_prec,         &
       t_vbs_group,        &
       t_aq_soa,           &
       vbs_nvocs,          &
       vbs_voc_prec,       &
       vbs_ngroup,         &
       vbs_set,            &
       aqsoa_ngroup,       &
       aqsoa_set

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: vbs_species

CONTAINS

  SUBROUTINE vbs_species
    
    ! -----------------------------------------------------------------------
    !
    ! SUBROUTINE vbs_species
    !
    ! depending on nvbs_setup:
    ! - creates the species (new_species; GAS) needed for the VOC precursers
    ! - creates the species (new_species; GAS_OR_AEROSOL) for the basis set
    ! - sets up the control structures vbs_voc_prec and vbs_set
    !  
    ! Authors:
    ! --------
    ! Declan O'Donnell, MPI-Met
    ! Kai Zhang, MPI-Met
    ! Thomas Kuehn, UEF    6/2015 --
    !
    ! vbs_species is called from ham_species in mo_ham_species
    !
    ! -----------------------------------------------------------------------

    USE mo_exception,          ONLY: &
         message,                    &
         message_text,               &
         em_info,                    &
         finish
    
    USE mo_physical_constants, ONLY: &
         argas,                      & ! ideal gas constant [J/K/mol]
         avo

    USE mo_tracdef,            ONLY: &
         GAS,                        &
         GAS_OR_AEROSOL,             &
         ON,                         &
         OFF,                        &
         itrprog

    USE mo_species,            ONLY: &
         speclist,                   & ! list of species defined in ham
         new_species                   ! for registering new ham species

    USE mo_ham_vbsctl,         ONLY: &
         laqsoa

    USE mo_ham_species,        ONLY: &
         id_oc, id_ocnv                ! for cross-referencing

    USE mo_ham_subm_species,   ONLY: &
         isubm_oc, isubm_ocnv          ! for cross-referencing

    USE mo_ham_rad_data,       ONLY: &
         iradoc
    
    INTEGER :: jb, jv

    INTEGER :: spid_temp

    ! -----------------------------------------------------------------------
    ! executable procedure
    ! -----------------------------------------------------------------------


    ! the additional aqueous phase soa species, if laqsoa is on
    IF (laqsoa) THEN
       aqsoa_ngroup = 2
    ELSE
       aqsoa_ngroup = 0
    ENDIF

    SELECT CASE(nvbs_setup)
       
       
    CASE(vbs_setup_3class_OC) ! 3 class scheme with OC as non-volatile

       ! Defining the VOC precursor species
       ! we use five VOC precursors, Monoterpenes and Isoprene from Megan,
       ! and Toluene, Xylene, and Benzene as anthropogenics 
       ! hence nvbs_vocs = 5
       vbs_nvocs = 5

       ! the basic VBS has two base species + OC
       vbs_ngroup = 3

       ! allocating memory for the precursor properties:
       IF (.NOT. ALLOCATED(vbs_voc_prec)) THEN
          ALLOCATE(vbs_voc_prec(vbs_nvocs))

          DO jv=1,vbs_nvocs
             ! allocating memory for the stoichiometric coefficients
             IF (.NOT. ALLOCATED(vbs_voc_prec(jv)%stoich_coeff)) THEN
                ALLOCATE(vbs_voc_prec(jv)%stoich_coeff(vbs_ngroup+aqsoa_ngroup))
             END IF
          END DO
       END IF

       ! allocating memory for the basis set:
       IF (.NOT. ALLOCATED(vbs_set)) THEN
          ALLOCATE(vbs_set(vbs_ngroup))
       END IF



       ! Defining the VOC species

       ! Monoterpenes
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Monoterpenes',           &
            shortname   = 'APIN',                   &
            units       = 'kg kg-1',                &
            mw          = 136._dp,                  &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(1)%spid &
            ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       ! OH (data from IUPAC)
       vbs_voc_prec(1)%k_0_OH     = 1.2E-11_dp ! pre-factor [m3/(mol*s)]
       vbs_voc_prec(1)%Eact_p_OH  = 440._dp    ! reduced activation energy [K]: Eact_p=Eact/R
       ! O3 (data from IUPAC)
       vbs_voc_prec(1)%k_0_O3     = 6.3E-16_dp ! pre-factor [m3/(mol*s)]
       vbs_voc_prec(1)%Eact_p_O3  = -580._dp   ! reduced activation energy [K]: Eact_p=Eact/R
       ! NO3
       vbs_voc_prec(1)%k_0_NO3    = 1.2E-12_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(1)%Eact_p_NO3 = 490._dp !reduced activation energy [K]: Eact_p=Eact/R

       ! Stoichiometric Coefficients (different amount depending on laqsoa)
       IF (laqsoa) THEN
          ! from Harri (hi NOx)
          !                            (  OC   ,  VBS1   ,  VBS10  , IEPOX , Glyx)
          vbs_voc_prec(1)%stoich_coeff=(/0.1_dp, 0.037_dp, 0.088_dp, 0.0_dp, 0.0_dp/) 
       ELSE
          !                            (  OC      ,  VBS1   , VBS10    )  
          vbs_voc_prec(1)%stoich_coeff=(/0.1_dp, 0.037_dp, 0.088_dp/) ! from Harri (hi NOx)
          !vbs_voc_prec(1)%stoich_coeff=(/0.002_dp, 0.003_dp, 0.065_dp/) ! from Harri (low NOx)
       ENDIF

       ! Isoprene                           
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Isoprene',               &
            shortname   = 'C5H8',                   &
            units       = 'kg kg-1',                &
            mw          = 68._dp,                   &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(2)%spid &
          ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_prec(2)%k_0_OH     = 2.7E-11_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(2)%Eact_p_OH  = 390._dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_prec(2)%k_0_O3     = 1.03E-14_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(2)%Eact_p_O3  = -1995._dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_prec(2)%k_0_NO3    = 3.15E-12_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(2)%Eact_p_NO3 = -450._dp !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC   ,  VBS1    ,   VBS10  ,   IEPOX , Glyx)
          vbs_voc_prec(2)%stoich_coeff=(/0.0_dp, 0.0295_dp, 0.0453_dp, 0.525_dp, 0.04_dp/) 
       ELSE
          !                            (  OC   ,  VBS1    , VBS10     )  
          vbs_voc_prec(2)%stoich_coeff=(/0.0_dp, 0.0295_dp, 0.0453_dp/) 
       ENDIF

       ! Toluene
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Toluene',                &
            shortname   = 'TOL',                    &
            units       = 'kg kg-1',                &
            mw          = 92._dp,                   &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(3)%spid &
            ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_prec(3)%k_0_OH     = 1.81E-12_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(3)%Eact_p_OH  = 338._dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_prec(3)%k_0_O3     = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(3)%Eact_p_O3  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_prec(3)%k_0_NO3    = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(3)%Eact_p_NO3 = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC    ,  VBS1 , VBS10 , IEPOX , Glyx    )
          vbs_voc_prec(3)%stoich_coeff=(/0.36_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.24_dp/) ! to be reviewed
       ELSE
          !                            (  OC    ,  VBS1 ,   VBS10)
          vbs_voc_prec(3)%stoich_coeff=(/0.36_dp, 0.0_dp, 0.0_dp/) ! to be reviewed
       ENDIF

       ! Xylene
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Xylene',                 &
            shortname   = 'XYL',                    &
            units       = 'kg kg-1',                &
            mw          = 106._dp,                  &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(4)%spid &
            ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_prec(4)%k_0_OH     = 2.31E-11_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(4)%Eact_p_OH  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_prec(4)%k_0_O3     = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(4)%Eact_p_O3  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_prec(4)%k_0_NO3    = 2.6E-16_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(4)%Eact_p_NO3 = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC    ,  VBS1 , VBS10 ,  IEPOX,  Glyx   )
          vbs_voc_prec(4)%stoich_coeff=(/0.30_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.25_dp/) ! to be reviewed
       ELSE
          !                            (  OC    ,  VBS1 , VBS10  )
          vbs_voc_prec(4)%stoich_coeff=(/0.30_dp, 0.0_dp, 0.0_dp/) ! to be reviewed
       ENDIF
		      
       ! Benzene
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Benzene',                &
            shortname   = 'BENZ',                   &
            units       = 'kg kg-1',                &
            mw          = 66._dp,                   &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(5)%spid &
            ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_prec(5)%k_0_OH     = 2.33E-12_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(5)%Eact_p_OH  = -193._dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_prec(5)%k_0_O3     = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(5)%Eact_p_O3  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_prec(5)%k_0_NO3    = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(5)%Eact_p_NO3 = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC    ,  VBS1 , VBS10 , IEPOX ,  Glyx   )
          vbs_voc_prec(5)%stoich_coeff=(/0.37_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.35_dp/) ! to be reviewed
       ELSE
          !                            (  OC    ,  VBS1 , VBS10  )
          vbs_voc_prec(5)%stoich_coeff=(/0.37_dp, 0.0_dp, 0.0_dp/) ! to be reviewed
       ENDIF


       ! OC (C*= 0 mug/m3)
       vbs_set(1)%lcreateaero = .FALSE.
       vbs_set(1)%spid        = id_ocnv  ! thk: !!!only works, if SALSA is on!!!
       vbs_set(1)%spid_aero   = id_oc    
       vbs_set(1)%C0          = 0.0_dp               ! equ. vapor conc. at T0 [mol/m3]
       vbs_set(1)%T0          = 298_dp               ! T0 [K]
       vbs_set(1)%Hvap_eff    = 30e3_dp/argas        ! eff. evap. enthalpy [K]
       vbs_set(1)%mv          =                    & ! molecular volume [m3]
            speclist(id_oc)%moleweight/avo/speclist(id_oc)%density*1e-3 ! [speclist%moleweight] = g/mol

       ! VBS group 1 (C*= 1 mug/m3)
       CALL new_species(&
            nphase      = GAS_OR_AEROSOL,               &
            longname    = 'Volatility Basis set C*=1',  &
            shortname   = 'VBS1'   ,                    &
            units       = 'kg kg-1',                    &
            mw          = 186._dp,                      & ! to be reviewed
            tsubmname   = 'HAM',                        &
            itrtype     = itrprog,                      &
            ldrydep     = .TRUE.,                       & ! to be reviewed
            lwetdep     = .TRUE.,                       & ! to be reviewed
            dryreac     = 0._dp,                        & ! to be reviewed
            henry       = (/ 1.E5_dp, 0._dp /),         & ! to be reviewed
            density     = 1320._dp,                     & ! to be reviewed
            iaerorad    = iradoc,                       & ! to be reviewed
            lwatsol     = .TRUE.,                       & ! to be reviewed
            kappa       = 0.037_dp,                     & ! Petters and Kreidenweis (2007)
            lburden     = .TRUE.,                       &
            idx         = spid_temp                     &
            ) 

       vbs_set(2)%lcreateaero = .TRUE.
       vbs_set(2)%spid        = spid_temp
       vbs_set(2)%spid_aero   = 0
       vbs_set(2)%C0          = &                    ! equ. vapor con. at T0 [mol/m3]
            1.0e-6_dp/speclist(spid_temp)%moleweight      ! [speclist%moleweight] = g/mol
       vbs_set(2)%T0          = 298_dp               ! T0 [K]
       vbs_set(2)%Hvap_eff    = 30e3_dp/argas        ! eff. evap. enthalpy [K]
       vbs_set(2)%mv          =                    & ! molecular volume [m3]
            speclist(spid_temp)%moleweight/avo/speclist(spid_temp)%density*1e-3 ! [speclist%moleweight] = g/mol


       ! VBS group 10 (C* = 10 mug/m3)
       CALL new_species(&
            nphase      = GAS_OR_AEROSOL,               &
            longname    = 'Volatility Basis Set C*=10', &
            shortname   = 'VBS10',                      &
            units       = 'kg kg-1',                    &
            mw          = 186._dp,                      & ! to be reviewed
            tsubmname   = 'HAM',                        &
            itrtype     = itrprog,                      &
            ldrydep     = .TRUE.,                       & ! to be reviewed
            lwetdep     = .TRUE.,                       & ! to be reviewed
            dryreac     = 0._dp,                        & ! to be reviewed
            henry       = (/ 1.E5_dp, 0._dp /),         & ! to be reviewed
            density     = 1320._dp,                     & ! to be reviewed
            iaerorad    = iradoc,                       & ! to be reviewed
            lwatsol     = .TRUE.,                       & ! to be reviewed
            kappa       = 0.037_dp,                     & !Petters and Kreidenweis (2007)
            lburden     = .TRUE.,                       & ! to be reviewed
            idx         = spid_temp                     &
            ) 


       vbs_set(3)%lcreateaero = .TRUE.
       vbs_set(3)%spid        = spid_temp
       vbs_set(3)%spid_aero   = 0
       vbs_set(3)%C0          = &                    ! equ. vapor con. at T0 [mol/m3]
            10.0e-6_dp/speclist(spid_temp)%moleweight      ! [speclist%moleweight] = g/mol
       vbs_set(3)%T0          = 298_dp               ! T0 [K]
       vbs_set(3)%Hvap_eff    = 30e3_dp/argas        ! eff. evap. enthalpy [K]
       vbs_set(3)%mv          =                    & ! molecular volume [m3]
            speclist(spid_temp)%moleweight/avo/speclist(spid_temp)%density*1e-3 ! [speclist%moleweight] = g/mol





    CASE(vbs_setup_3class) ! 3 class scheme with volatilities 0, 1, and 10 ug/m3

       ! Defining the VOC precursor species
       ! we use five VOC precursors, Monoterpenes and Isoprene from Megan,
       ! and Toluene, Xylene, and Benzene as anthropogenics 
       ! hence nvbs_vocs = 5
       vbs_nvocs = 5

       ! the VBS has three volatility groups
       vbs_ngroup = 3

       ! allocating memory for the precursor properties:
       IF (.NOT. ALLOCATED(vbs_voc_prec)) THEN
          ALLOCATE(vbs_voc_prec(vbs_nvocs))

          DO jv=1,vbs_nvocs
             ! allocating memory for the stoichiometric coefficients
             IF (.NOT. ALLOCATED(vbs_voc_prec(jv)%stoich_coeff)) THEN
                ALLOCATE(vbs_voc_prec(jv)%stoich_coeff(vbs_ngroup+aqsoa_ngroup))
             END IF
          END DO
       END IF

       ! allocating memory for the basis set:
       IF (.NOT. ALLOCATED(vbs_set)) THEN
          ALLOCATE(vbs_set(vbs_ngroup))
       END IF



       ! Defining the VOC species

       ! Monoterpenes
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Monoterpenes',           &
            shortname   = 'APIN',                   &
            units       = 'kg kg-1',                &
            mw          = 136._dp,                  &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(1)%spid &
            ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       ! OH (data from IUPAC)
       vbs_voc_prec(1)%k_0_OH     = 1.2E-11_dp ! pre-factor [m3/(mol*s)]
       vbs_voc_prec(1)%Eact_p_OH  = 440._dp    ! reduced activation energy [K]: Eact_p=Eact/R
       ! O3 (data from IUPAC)
       vbs_voc_prec(1)%k_0_O3     = 6.3E-16_dp ! pre-factor [m3/(mol*s)]
       vbs_voc_prec(1)%Eact_p_O3  = -580._dp   ! reduced activation energy [K]: Eact_p=Eact/R
       ! NO3
       vbs_voc_prec(1)%k_0_NO3    = 1.2E-12_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(1)%Eact_p_NO3 = 490._dp !reduced activation energy [K]: Eact_p=Eact/R

       ! Stoichiometric Coefficients (different amount depending on laqsoa)
       IF (laqsoa) THEN
          ! from Harri (hi NOx)
          !                            (  OC   ,  VBS1   ,  VBS10  , IEPOX , Glyx)
          vbs_voc_prec(1)%stoich_coeff=(/0.1_dp, 0.037_dp, 0.088_dp, 0.0_dp, 0.0_dp/) 
       ELSE
          !                            (  OC      ,  VBS1   , VBS10    )  
          vbs_voc_prec(1)%stoich_coeff=(/0.1_dp, 0.037_dp, 0.088_dp/) ! from Harri (hi NOx)
          !vbs_voc_prec(1)%stoich_coeff=(/0.002_dp, 0.003_dp, 0.065_dp/) ! from Harri (low NOx)
       ENDIF

       ! Isoprene                           
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Isoprene',               &
            shortname   = 'C5H8',                   &
            units       = 'kg kg-1',                &
            mw          = 68._dp,                   &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(2)%spid &
          ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_prec(2)%k_0_OH     = 2.7E-11_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(2)%Eact_p_OH  = 390._dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_prec(2)%k_0_O3     = 1.03E-14_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(2)%Eact_p_O3  = -1995._dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_prec(2)%k_0_NO3    = 3.15E-12_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(2)%Eact_p_NO3 = -450._dp !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC   ,  VBS1    ,   VBS10  ,   IEPOX , Glyx)
          vbs_voc_prec(2)%stoich_coeff=(/0.0_dp, 0.0295_dp, 0.0453_dp, 0.525_dp, 0.04_dp/) 
       ELSE
          !                            (  OC   ,  VBS1    , VBS10     )  
          vbs_voc_prec(2)%stoich_coeff=(/0.0_dp, 0.0295_dp, 0.0453_dp/) 
       ENDIF

       ! Toluene
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Toluene',                &
            shortname   = 'TOL',                    &
            units       = 'kg kg-1',                &
            mw          = 92._dp,                   &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(3)%spid &
            ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_prec(3)%k_0_OH     = 1.81E-12_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(3)%Eact_p_OH  = 338._dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_prec(3)%k_0_O3     = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(3)%Eact_p_O3  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_prec(3)%k_0_NO3    = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(3)%Eact_p_NO3 = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC    ,  VBS1 , VBS10 , IEPOX , Glyx    )
          vbs_voc_prec(3)%stoich_coeff=(/0.36_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.24_dp/) ! to be reviewed
       ELSE
          !                            (  OC    ,  VBS1 ,   VBS10)
          vbs_voc_prec(3)%stoich_coeff=(/0.36_dp, 0.0_dp, 0.0_dp/) ! to be reviewed
       ENDIF

       ! Xylene
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Xylene',                 &
            shortname   = 'XYL',                    &
            units       = 'kg kg-1',                &
            mw          = 106._dp,                  &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(4)%spid &
            ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_prec(4)%k_0_OH     = 2.31E-11_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(4)%Eact_p_OH  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_prec(4)%k_0_O3     = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(4)%Eact_p_O3  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_prec(4)%k_0_NO3    = 2.6E-16_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(4)%Eact_p_NO3 = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC    ,  VBS1 , VBS10 ,  IEPOX,  Glyx   )
          vbs_voc_prec(4)%stoich_coeff=(/0.30_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.25_dp/) ! to be reviewed
       ELSE
          !                            (  OC    ,  VBS1 , VBS10  )
          vbs_voc_prec(4)%stoich_coeff=(/0.30_dp, 0.0_dp, 0.0_dp/) ! to be reviewed
       ENDIF
		      
       ! Benzene
       CALL new_species(&
            nphase      = GAS,                      &
            longname    = 'Benzene',                &
            shortname   = 'BENZ',                   &
            units       = 'kg kg-1',                &
            mw          = 66._dp,                   &
            tsubmname   = 'HAM',                    &
            itrtype     = itrprog,                  &
            ldrydep     = .FALSE.,                  &
            lwetdep     = .FALSE.,                  &
            lburden     = .TRUE.,                   &
            idx         = vbs_voc_prec(5)%spid &
            ) 
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_prec(5)%k_0_OH     = 2.33E-12_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(5)%Eact_p_OH  = -193._dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_prec(5)%k_0_O3     = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(5)%Eact_p_O3  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_prec(5)%k_0_NO3    = 0.0_dp !pre-factor [m3/(mol*s)]
       vbs_voc_prec(5)%Eact_p_NO3 = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC    ,  VBS1 , VBS10 , IEPOX ,  Glyx   )
          vbs_voc_prec(5)%stoich_coeff=(/0.37_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.35_dp/) ! to be reviewed
       ELSE
          !                            (  OC    ,  VBS1 , VBS10  )
          vbs_voc_prec(5)%stoich_coeff=(/0.37_dp, 0.0_dp, 0.0_dp/) ! to be reviewed
       ENDIF


       ! VBS group 0 (C*= 1 mug/m3)
       CALL new_species(&
            nphase      = GAS_OR_AEROSOL,               &
            longname    = 'Volatility Basis set C*=0',  &
            shortname   = 'VBS0'   ,                    &
            units       = 'kg kg-1',                    &
            mw          = 186._dp,                      & ! to be reviewed
            tsubmname   = 'HAM',                        &
            itrtype     = itrprog,                      &
            ldrydep     = .TRUE.,                       & ! to be reviewed
            lwetdep     = .TRUE.,                       & ! to be reviewed
            dryreac     = 0._dp,                        & ! to be reviewed
            henry       = (/ 1.E5_dp, 0._dp /),         & ! to be reviewed
            density     = 1320._dp,                     & ! to be reviewed
            iaerorad    = iradoc,                       & ! to be reviewed
            lwatsol     = .TRUE.,                       & ! to be reviewed
            kappa       = 0.037_dp,                     & ! Petters and Kreidenweis (2007)
            lburden     = .TRUE.,                       &
            idx         = spid_temp                     &
            ) 

       vbs_set(1)%lcreateaero = .TRUE.
       vbs_set(1)%spid        = spid_temp
       vbs_set(1)%spid_aero   = 0
       vbs_set(1)%C0          = &                    ! equ. vapor con. at T0 [mol/m3]
            0.0e-6_dp/speclist(spid_temp)%moleweight      ! [speclist%moleweight] = g/mol
       vbs_set(1)%T0          = 298_dp               ! T0 [K]
       vbs_set(1)%Hvap_eff    = 30e3_dp/argas        ! eff. evap. enthalpy [K]
       vbs_set(1)%mv          =                    & ! molecular volume [m3]
            speclist(spid_temp)%moleweight/avo/speclist(spid_temp)%density*1e-3 ! [speclist%moleweight] = g/mol


       ! VBS group 1 (C*= 1 mug/m3)
       CALL new_species(&
            nphase      = GAS_OR_AEROSOL,               &
            longname    = 'Volatility Basis set C*=1',  &
            shortname   = 'VBS1'   ,                    &
            units       = 'kg kg-1',                    &
            mw          = 186._dp,                      & ! to be reviewed
            tsubmname   = 'HAM',                        &
            itrtype     = itrprog,                      &
            ldrydep     = .TRUE.,                       & ! to be reviewed
            lwetdep     = .TRUE.,                       & ! to be reviewed
            dryreac     = 0._dp,                        & ! to be reviewed
            henry       = (/ 1.E5_dp, 0._dp /),         & ! to be reviewed
            density     = 1320._dp,                     & ! to be reviewed
            iaerorad    = iradoc,                       & ! to be reviewed
            lwatsol     = .TRUE.,                       & ! to be reviewed
            kappa       = 0.037_dp,                     & ! Petters and Kreidenweis (2007)
            lburden     = .TRUE.,                       &
            idx         = spid_temp                     &
            ) 

       vbs_set(2)%lcreateaero = .TRUE.
       vbs_set(2)%spid        = spid_temp
       vbs_set(2)%spid_aero   = 0
       vbs_set(2)%C0          = &                    ! equ. vapor con. at T0 [mol/m3]
            1.0e-6_dp/speclist(spid_temp)%moleweight      ! [speclist%moleweight] = g/mol
       vbs_set(2)%T0          = 298_dp               ! T0 [K]
       vbs_set(2)%Hvap_eff    = 30e3_dp/argas        ! eff. evap. enthalpy [K]
       vbs_set(2)%mv          =                    & ! molecular volume [m3]
            speclist(spid_temp)%moleweight/avo/speclist(spid_temp)%density*1e-3 ! [speclist%moleweight] = g/mol


       ! VBS group 10 (C* = 10 mug/m3)
       CALL new_species(&
            nphase      = GAS_OR_AEROSOL,               &
            longname    = 'Volatility Basis Set C*=10', &
            shortname   = 'VBS10',                      &
            units       = 'kg kg-1',                    &
            mw          = 186._dp,                      & ! to be reviewed
            tsubmname   = 'HAM',                        &
            itrtype     = itrprog,                      &
            ldrydep     = .TRUE.,                       & ! to be reviewed
            lwetdep     = .TRUE.,                       & ! to be reviewed
            dryreac     = 0._dp,                        & ! to be reviewed
            henry       = (/ 1.E5_dp, 0._dp /),         & ! to be reviewed
            density     = 1320._dp,                     & ! to be reviewed
            iaerorad    = iradoc,                       & ! to be reviewed
            lwatsol     = .TRUE.,                       & ! to be reviewed
            kappa       = 0.037_dp,                     & !Petters and Kreidenweis (2007)
            lburden     = .TRUE.,                       & ! to be reviewed
            idx         = spid_temp                     &
            ) 


       vbs_set(3)%lcreateaero = .TRUE.
       vbs_set(3)%spid        = spid_temp
       vbs_set(3)%spid_aero   = 0
       vbs_set(3)%C0          = &                    ! equ. vapor con. at T0 [mol/m3]
            10.0e-6_dp/speclist(spid_temp)%moleweight      ! [speclist%moleweight] = g/mol
       vbs_set(3)%T0          = 298_dp               ! T0 [K]
       vbs_set(3)%Hvap_eff    = 30e3_dp/argas        ! eff. evap. enthalpy [K]
       vbs_set(3)%mv          =                    & ! molecular volume [m3]
            speclist(spid_temp)%moleweight/avo/speclist(spid_temp)%density*1e-3 ! [speclist%moleweight] = g/mol





    CASE DEFAULT
       WRITE (message_text,'(a,i0)') &
            'ERROR: No scheme is implemented for nvbs_setup =',&
            nvbs_setup

       CALL message(&
            '',&
            message_text&
            )

       CALL finish(&
            'vbs_species',&
            'run terminated'&
            )
    END SELECT



    !>> thk: added aqueous phase soa directly
    IF (laqsoa) THEN

       ! allocating memory for the wet SOA set:
       IF (.NOT. ALLOCATED(aqsoa_set)) THEN
          ALLOCATE(aqsoa_set(aqsoa_ngroup))
       END IF

       ! Isoprene epoxide
       CALL new_species(&
            nphase      = GAS_OR_AEROSOL,               &
            longname    = 'Isoprene epoxide',         &
            shortname   = 'IEPOX',                       &
            units       = 'kg kg-1',                    &
            mw          = 118._dp,                      & ! to be reviewed
            tsubmname   = 'HAM',                        &
            itrtype     = itrprog,                      &
            ldrydep     = .TRUE.,                       & ! to be reviewed
            lwetdep     = .TRUE.,                       & ! to be reviewed
            dryreac     = 0._dp,                        & ! to be reviewed
            henry       = (/ 1.E5_dp, 0._dp /),         & ! to be reviewed
            density     = 1320._dp,                     & ! to be reviewed
            iaerorad    = iradoc,                       & ! to be reviewed
            lwatsol     = .TRUE.,                       & ! to be reviewed
            kappa       = 0.037_dp,                     & !Petters and Kreidenweis (2007)
            lburden     = .TRUE.,                       & ! to be reviewed
            idx         = spid_temp                     &
            ) 

       ! Physical paarameters
       aqsoa_set(1)%lcreateaero = .TRUE.
       aqsoa_set(1)%spid        = spid_temp
       aqsoa_set(1)%spid_aero   = 0
       aqsoa_set(1)%mv           = 1.13E-28_dp     ! molecular volume [m3] (to be reviewed!!) (assuming 0.3nm radius)              
       aqsoa_set(1)%Eff_henry_water   = 1.0E5_dp        ! Effective Henry's constant for cloud droplets          
       aqsoa_set(1)%Eff_henry_aerosol = 1.0E8_dp        ! Effective Henry's constant for aerosols           
       !OH
       !aqsoa_set(1)%k_0_OH     = 1.25E-11 !pre-factor [m3/(mol*s)] ((Bates et al.. Average between cis and trans isomers)
       aqsoa_set(1)%k_0_OH     = 3.56E-11_dp !pre-factor [m3/(mol*s)] ((Jacobs et al.. average between IEPOX1 and IEPOX4)
       aqsoa_set(1)%Eact_p_OH  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       aqsoa_set(1)%k_0_O3     = 0.0_dp !pre-factor [m3/(mol*s)]
       aqsoa_set(1)%Eact_p_O3  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       aqsoa_set(1)%k_0_NO3    = 0.0_dp !pre-factor [m3/(mol*s)]
       aqsoa_set(1)%Eact_p_NO3 = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !Photodissociation 
       aqsoa_set(1)%photodis   = 0.0_dp !desctruction rate [% s-1] with sunlight       


       ! Glyoxal
       CALL new_species(&
            nphase      = GAS_OR_AEROSOL,               &
            longname    = 'Glyoxal',         &
            shortname   = 'Glyx',                       &
            units       = 'kg kg-1',                    &
            mw          = 58._dp,                      & ! to be reviewed
            tsubmname   = 'HAM',                        &
            itrtype     = itrprog,                      &
            ldrydep     = .TRUE.,                       & ! to be reviewed
            lwetdep     = .TRUE.,                       & ! to be reviewed
            dryreac     = 0._dp,                        & ! to be reviewed
            henry       = (/ 4.19E5_dp, 0._dp /),        & ! from Nguyen et al. 2014. Value for pure water.
            density     = 1320._dp,                     & ! to be reviewed
            iaerorad    = iradoc,                       & ! to be reviewed
            lwatsol     = .TRUE.,                       & ! to be reviewed
            kappa       = 0.037_dp,                     & !Petters and Kreidenweis (2007)
            lburden     = .TRUE.,                       & ! to be reviewed
            idx         = spid_temp                     &
            ) 


       ! Physical paarameters
       aqsoa_set(2)%lcreateaero = .TRUE.
       aqsoa_set(2)%spid        = spid_temp
       aqsoa_set(2)%spid_aero   = 0
       aqsoa_set(2)%mv          = 1.13E-28_dp     ! molecular volume [m3] (to be reviewed!!) (assuming 0.3nm radius)
       aqsoa_set(2)%Eff_henry_water   = 4.19E5_dp  ! Effective Henry's constant for cloud droplets
       aqsoa_set(2)%Eff_henry_aerosol = 3.0E8_dp  ! Effective Henry's constant for aerosols (Kampf et al. 2013)
       !OH
       aqsoa_set(2)%k_0_OH     = 0.0_dp !pre-factor [m3/(mol*s)]
       aqsoa_set(2)%Eact_p_OH  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       aqsoa_set(2)%k_0_O3     = 0.0_dp !pre-factor [m3/(mol*s)]
       aqsoa_set(2)%Eact_p_O3  = 0.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       aqsoa_set(2)%k_0_NO3    = 6.E-13_dp !pre-factor [m3/(mol*s)]
       aqsoa_set(2)%Eact_p_NO3 = -1900.0_dp !reduced activation energy [K*mol]: Eact_p = Eact/R
       !Photodissociation 
       aqsoa_set(2)%photodis   = 3.141593_dp/2._dp*8.21E-5 !desctruction rate [% s-1] with sunlight        

    END IF !laqsoa

    !<< thk
    
    



  END SUBROUTINE vbs_species




END MODULE mo_ham_vbs
