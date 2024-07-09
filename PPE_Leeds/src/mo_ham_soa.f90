!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_soa.f90
!!
!! \brief
!! This module contains HAM routines related to the handling of secondary organic aerosols.
!!
!! \author  Declan O'Donnell (MPI-M)
!!
!! \responsible_coder
!! Declan O'Donnell, declan.Odonnell@fmi.fi 
!!
!! \revision_history
!!   -# Declan O'Donnell (MPI-Met) - original code (YYYY)
!!   -# Kai Zhang (MPI-Met) - rewrite species registeration part, change longname
!!                            and shortname to avoid duplicate longnames in species list (2009-08)
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

MODULE mo_ham_soa

  !--- inherited types, data and functions
  USE mo_kind,         ONLY: dp       
  USE mo_ham,         ONLY: nsoaspec, new_aerocomp
  USE mo_tracdef,      ONLY: ON, OFF, itrnone, itrdiag, itrprog, itrpresc
  USE mo_ham_species, ONLY: id_oc
  USE mo_decomposition,ONLY: lc => local_decomposition
  USE mo_species,      ONLY: t_species, new_species
  USE mo_tracdef,      ONLY: GAS, AEROSOL, GAS_OR_AEROSOL, ON, OFF

  IMPLICIT NONE  

  !---public member functions
  
  PUBLIC :: soa_species
  !>>dod deleted start_map_soaspec <<dod
  PUBLIC :: start_soa_aerosols
  PUBLIC :: construct_soa_streams
  !>>dod 
  PUBLIC :: set_soa_tracer_attr
  !<<dod

  !---public data types 
  !>>dod deleted dead code (t_soa_rxn)
  !<< dod
  TYPE, PUBLIC :: t_soa_prop                           ! SOA properties
     LOGICAL  :: lvolatile                             ! Volatile species T/F
     REAL(dp) :: Kp                                    ! 2-product model parameter: 
                                                       ! partitioning coefficient 
     REAL(dp) :: tref                                  ! 2-product model parameter: 
                                                       ! temperature at which Kp was derived
     REAL(dp) :: dH                                    ! 2-product model parameter:
                                                       ! Enthalpy of vaporisation
     INTEGER  :: spid_tot                              ! Species id, prognostic species (gas+aerosol total mass)
     INTEGER  :: spid_soa                              ! Species id, diagnostic species (gas or aerosol)
  END TYPE t_soa_prop                                  

!!mgs!! removed obsolete t_soa_species and replaced with soagas_idx

  !---public module data

  !---maximum number of SOA species
  INTEGER, PARAMETER, PUBLIC :: nmaxsoa = 7

  !---species identities for SOA precursors, gases and aerosols  ### now in mo_soa_species
  !   soa precursors
  !>>csld #404 
  !define new monoterpene species for use of MEGAN2.1
  INTEGER, PUBLIC :: id_apin, id_tbeta, id_bpin, id_lim, id_sab, id_myrc, id_car
  !<<csld #404 
  INTEGER, PUBLIC :: id_isop, id_tol, id_xyl, id_benz

  !---2-product model properties
  TYPE(t_soa_prop), PUBLIC, TARGET :: soaprop(nmaxsoa)         

  !---individual SOA species
!!mgs!!  TYPE(t_soa_species), PUBLIC :: soaspecies(nmaxsoa)

  INTEGER, PUBLIC :: nnvol                 ! number of non-volatile species
  !>>dod removed dead code (soa_rxn)
  !<< dod

  !---model functional switches
  LOGICAL, PUBLIC  :: lso4_in_m0 = .FALSE.                ! include sulphate in SOA absorbing mass?

  !---diagnostic fields...
  REAL(dp), PUBLIC, POINTER :: d_prod_soa_mterp(:,:)       
  REAL(dp), PUBLIC, POINTER :: d_prod_soa_isop(:,:)       
  REAL(dp), PUBLIC, POINTER :: d_prod_soa_tol(:,:)       
  REAL(dp), PUBLIC, POINTER :: d_prod_soa_xyl(:,:)       
  REAL(dp), PUBLIC, POINTER :: d_prod_soa_benz(:,:)       
  REAL(dp), PUBLIC, POINTER :: d_chem_sink_mterp(:,:)       
  REAL(dp), PUBLIC, POINTER :: d_chem_sink_isop(:,:)       
  REAL(dp), PUBLIC, POINTER :: d_chem_sink_tol(:,:)       
  REAL(dp), PUBLIC, POINTER :: d_chem_sink_xyl(:,:)       
  REAL(dp), PUBLIC, POINTER :: d_chem_sink_benz(:,:)

  !-----------------------------------------------------------------------------------------------

  CONTAINS
  !-----------------------------------------------------------------------------------------------
  
  
  
    SUBROUTINE soa_species
!### NOTE: species type shall be changed to cover both gas-phase and aerosol definitions in one structure
!### Will require nphase as bit flags (1 = gas, 2=aerosol [4=liquid...]) - only one species per oxidation product then.

      USE mo_ham,           ONLY: nsoalumping
!++mgs
      USE mo_ham_rad_data,  ONLY: iradsoa
!--mgs
      !>> dod
      ! USE mo_submodel,      ONLY: lham, lhammoz
      USE mo_tracdef,       ONLY: GAS, AEROSOL, GAS_OR_AEROSOL
      !<<
      IMPLICIT NONE

      !---local data:
      INTEGER :: i
      !---executable procedure
      

      !---create precursors
      !   Monoterpenes
      !>>dod changed lwetdep to FALSE on the basis that it is negligible in runs where it was used

!>>csld #404 new species for use of MEGAN2.1

      CALL new_species(nphase      = GAS,                 &
                       longname    = 'alpha pinene',      &
                       shortname   = 'APIN',              &
                       units       = 'kg kg-1',           &
                       mw          = 136._dp,             &
                       tsubmname   = 'HAM',               &
                       itrtype     = itrprog,             &
                       ldrydep     = .FALSE.,             &
                       lwetdep     = .FALSE.,             &
                       lburden     = .TRUE.,              &
                       idx         = id_apin              ) 

      CALL new_species(nphase      = GAS,                 &
                       longname    = 't beta ocimene',    &
                       shortname   = 'TBETAOCI',          &
                       units       = 'kg kg-1',           &
                       mw          = 136._dp,             &
                       tsubmname   = 'HAM',               &
                       itrtype     = itrprog,             &
                       ldrydep     = .FALSE.,             &
                       lwetdep     = .FALSE.,             &
                       lburden     = .TRUE.,              &
                       idx         = id_tbeta             ) 

      CALL new_species(nphase      = GAS,                 &
                       longname    = 'b pinene',          &
                       shortname   = 'BPIN',              &
                       units       = 'kg kg-1',           &
                       mw          = 136._dp,             &
                       tsubmname   = 'HAM',               &
                       itrtype     = itrprog,             &
                       ldrydep     = .FALSE.,             &
                       lwetdep     = .FALSE.,             &
                       lburden     = .TRUE.,              &
                       idx         = id_bpin              ) 

      CALL new_species(nphase      = GAS,                 &
                       longname    = 'limonene',          &
                       shortname   = 'LIMON',             &
                       units       = 'kg kg-1',           &
                       mw          = 136._dp,             &
                       tsubmname   = 'HAM',               &
                       itrtype     = itrprog,             &
                       ldrydep     = .FALSE.,             &
                       lwetdep     = .FALSE.,             &
                       lburden     = .TRUE.,              &
                       idx         = id_lim               ) 

      CALL new_species(nphase      = GAS,                 &
                       longname    = 'sabinene',          &
                       shortname   = 'SABIN',             &
                       units       = 'kg kg-1',           &
                       mw          = 136._dp,             &
                       tsubmname   = 'HAM',               &
                       itrtype     = itrprog,             &
                       ldrydep     = .FALSE.,             &
                       lwetdep     = .FALSE.,             &
                       lburden     = .TRUE.,              &
                       idx         = id_sab               ) 

      CALL new_species(nphase      = GAS,                 &
                       longname    = 'myrcene',           &
                       shortname   = 'MYRC',              &
                       units       = 'kg kg-1',           &
                       mw          = 136._dp,             &
                       tsubmname   = 'HAM',               &
                       itrtype     = itrprog,             &
                       ldrydep     = .FALSE.,             &
                       lwetdep     = .FALSE.,             &
                       lburden     = .TRUE.,              &
                       idx         = id_myrc              ) 

      CALL new_species(nphase      = GAS,                 &
                       longname    = '3-carene',          &
                       shortname   = 'CARENE3',           &
                       units       = 'kg kg-1',           &
                       mw          = 136._dp,             &
                       tsubmname   = 'HAM',               &
                       itrtype     = itrprog,             &
                       ldrydep     = .FALSE.,             &
                       lwetdep     = .FALSE.,             &
                       lburden     = .TRUE.,              &
                       idx         = id_car               ) 
!<<csld #404

      !   Isoprene                           
     
!++mgs: are you sure that isoprene is wet and dry deposited ?? 
     !>>no, but I am not sure that it is not, either. It is reactive stuff, so dry dep.
     !  is kept, wet dep is removed, as for monoterpenes...might delete drydep too, align with 
     !  the chemistry model
     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Isoprene',          &
                      shortname   = 'C5H8',              &
                      units       = 'kg kg-1',           &
                      mw          = 68._dp,              &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrprog,             &
                      ldrydep     = .FALSE.,             &
                      lwetdep     = .FALSE.,             &
                      lburden     = .TRUE.,              &
                      idx         = id_isop                 ) 
     !<<dod
		      
      !   Toluene
		
     !>> dod anthropogenics: removed henry since wetdep=F
     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Toluene',           &
                      shortname   = 'TOL',               &
                      units       = 'kg kg-1',           &
                      mw          = 92._dp,              &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrprog,             &
                      ldrydep     = .FALSE.,             &
                      lwetdep     = .FALSE.,             &
                      lburden     = .TRUE.,              &
                      idx         = id_tol                  ) 

      !   Xylene
			       
     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Xylene',            &
                      shortname   = 'XYL',               &
                      units       = 'kg kg-1',           &
                      mw          = 106._dp,             &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrprog,             &
                      ldrydep     = .FALSE.,             &
                      lwetdep     = .FALSE.,             &
                      lburden     = .TRUE.,              &
                      idx         = id_xyl                  ) 
		      
      !   Benzene
                 
     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Benzene',           &
                      shortname   = 'BENZ',              &
                      units       = 'kg kg-1',           &
                      mw          = 66._dp,              &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrprog,             &
                      ldrydep     = .FALSE.,             &
                      lwetdep     = .FALSE.,             &
                      lburden     = .TRUE.,              &
                      idx         = id_benz                 ) 
		      
     !<< dod

      !---create biogenic gas phase SOA species
!++mgs: note: ltrreq now moved to aerosol species to avoid duplicate subname 'gas'

     !>> dod
     ! semi-volatile species are handled in a special way. Such species may exist both in the
     ! gas and aerosol phases, and in many different aerosol size modes. Partitioning between
     ! gas and aerosol is done with an equilibrium model. This has the consequence that from the
     ! point of view of tracer transport, we do not actually need separate aerosol and gas phase
     ! tracers: since both gas and aerosol are transported in the same way, we can simply transport
     ! the total mass (gas+aerosol) and then diagnose the gas and aerosol masses using the 
     ! equilibrium model. This cuts the tracer transport requirements by a factor of 5 or so. 
     ! For purposes of calculating radiative effects, wet deposition,
     ! etc., it is convenient to also implement the separated gas and aerosols as tracers.
  
     ! First the prognostic species (total). Define it as a gas, otherwise HAM will
     ! allocate many tracers.
     CALL new_species(nphase      = GAS,                          &
                      longname    = 'Monoterpene SOA Total 1',    &
                      shortname   = 'SOA_MT_T1',                  &
                      units       = 'kg kg-1',                    &
                      mw          = 186._dp,                      &
                      tsubmname   = 'HAM',                        &
                      itrtype     = itrprog,                      &
                      lburden     = .FALSE.,                      &
                      idx         = soaprop(1)%spid_tot           ) 


     CALL new_species(nphase      = GAS,                          &
                      longname    = 'Monoterpene SOA Total 2',    &
                      shortname   = 'SOA_MT_T2',                  &
                      units       = 'kg kg-1',                    &
!gf                      mw          = 186._dp,                      &
                      mw          = 168._dp,                      &
                      tsubmname   = 'HAM',                        &
                      itrtype     = itrprog,                      &
                      lburden     = .FALSE.,                      &
                      idx         = soaprop(2)%spid_tot           ) 

     
     ! now the gases and aerosols...
     CALL new_species(nphase      = GAS_OR_AEROSOL,               &
                      longname    = 'Monoterpene SOA 1',          &
                      shortname   = 'SOA_MT1',                    &
                      units       = 'kg kg-1',                    &
                      mw          = 186._dp,                      &
                      tsubmname   = 'HAM',                        &
                      itrtype     = itrdiag,                      &
                      ldrydep     = .TRUE.,                       &
                      lwetdep     = .TRUE.,                       &
                      dryreac     = 0._dp,                        & 
!gf(#94)                      henry       = (/ 1.E4_dp, 0._dp /),         &
                      henry       = (/ 1.E5_dp, 0._dp /),         &
                      density     = 1320._dp,                     &
                      iaerorad    = iradsoa,                      &
                      lwatsol     = .TRUE.,                       &  
                      kappa       = 0.06_dp,                       &     !>>dod <<
                      lburden     = .TRUE.,                       &
                      idx         = soaprop(1)%spid_soa           ) 


     CALL new_species(nphase      = GAS_OR_AEROSOL,               &
                      longname    = 'Monoterpene SOA 2',          &
                      shortname   = 'SOA_MT2',                    &
                      units       = 'kg kg-1',                    &
!gf                      mw          = 186._dp,                      &     
                      mw          = 168._dp,                      &     
                      tsubmname   = 'HAM',                        &
                      itrtype     = itrdiag,                      &
                      ldrydep     = .TRUE.,                       &
                      lwetdep     = .TRUE.,                       &
                      dryreac     = 0._dp,                        & 
                      henry       = (/ 1.E5_dp, 0._dp /),         &
                      density     = 1320._dp,                     &
                      iaerorad    = iradsoa,                      &
                      lwatsol     = .TRUE.,                       &  
                      kappa       = 0.06_dp,                       &  !>>dod <<
                      lburden     = .TRUE.,                       &
                      idx         = soaprop(2)%spid_soa           ) 


     ! End of monoterpene definition

     ! Isoprene

     ! Prognostic species (total mass)
     CALL new_species(nphase      = GAS,                            &
                      longname    = 'Isoprene SOA Total 1',         &
                      shortname   = 'SOA_IS_T1',                    &
                      units       = 'kg kg-1',                      &
                      mw          = 125._dp,                        &
                      tsubmname   = 'HAM',                          &
                      itrtype     = itrprog,                        &
                      lburden     = .FALSE.,                        &
                      idx         = soaprop(3)%spid_tot             ) 

     CALL new_species(nphase      = GAS,                            &
                      longname    = 'Isoprene SOA Total 2',         &
                      shortname   = 'SOA_IS_T2',                    &
                      units       = 'kg kg-1',                      &
                      mw          = 125._dp,                        &
                      tsubmname   = 'HAM',                          &
                      itrtype     = itrprog,                        &
                      lburden     = .FALSE.,                        &
                      idx         = soaprop(4)%spid_tot             ) 
				      
     ! Diagnostic species. 
     CALL new_species(nphase      = GAS_OR_AEROSOL,                 &
                      longname    = 'Isoprene SOA 1',               &
                      shortname   = 'SOA_IS1',                      &
                      units       = 'kg kg-1',                      &
                      mw          = 125._dp,                        &
                      tsubmname   = 'HAM',                          &
                      itrtype     = itrdiag,                        &
                      ldrydep     = .TRUE.,                         &
                      lwetdep     = .TRUE.,                         &
                      dryreac     = 0._dp,                          & 
!gf(#94)                      henry       = (/ 1.E4_dp, 0._dp /),           &
                      henry       = (/ 1.E5_dp, 0._dp /),           &
                      density     = 1320._dp,                       &
                      iaerorad    = iradsoa,                        &
                      lwatsol     = .TRUE.,                         &  
                      kappa       = 0.06_dp,                         &  !>>dod <<
                      lburden     = .TRUE.,                         &
                      idx         = soaprop(3)%spid_soa             ) 

     CALL new_species(nphase      = GAS_OR_AEROSOL,                 &
                      longname    = 'Isoprene SOA 2',               &
                      shortname   = 'SOA_IS2',                      &
                      units       = 'kg kg-1',                      &
                      mw          = 125._dp,                        &
                      tsubmname   = 'HAM',                          &
                      itrtype     = itrdiag,                        &
                      ldrydep     = .TRUE.,                         &
                      lwetdep     = .TRUE.,                         &
                      dryreac     = 0._dp,                          & 
                      henry       = (/ 1.E5_dp, 0._dp /),           &
                      lburden     = .TRUE.,                         &
                      density     = 1320._dp,                       &
                      iaerorad    = iradsoa,                        &
                      lwatsol     = .TRUE.,                         &  
                      kappa       = 0.06_dp,                         & !>>dod <<
                      idx         = soaprop(4)%spid_soa             ) 

     ! End isoprene definitons

     !---create anthropogenic gas phase SOA species according to lumping parameter
     ! 
     SELECT CASE(nsoalumping)
     CASE(0)                              ! no lumping, all anthropogenics distinct 

        CALL new_species(nphase      = AEROSOL,                      &
                         longname    = 'Toluene SOA',                &
                         shortname   = 'SOA_TOL',                    &
                         units       = 'kg kg-1',                    &
                         mw          = 124._dp,                      & 
                         tsubmname   = 'HAM',                        &
                         itrtype     = itrprog,                      &
                         ldrydep     = .TRUE.,                       &
                         lwetdep     = .TRUE.,                       &
                         density     = 1450._dp,                     &
                         iaerorad    = iradsoa,                      &
                         lwatsol     = .TRUE.,                       &  
                         kappa       = 0.06_dp,                       &  !>>dod <<
                         lburden     = .TRUE.,                       &
                         idx         = soaprop(5)%spid_soa               ) 


        CALL new_species(nphase      = AEROSOL,                      &
                         longname    = 'Xylene SOA',                 &
                         shortname   = 'SOA_XYL',                    &
                         units       = 'kg kg-1',                    &
                         mw          = 138._dp,                      & 
                         tsubmname   = 'HAM',                        &
                         itrtype     = itrprog,                      &
                         ldrydep     = .TRUE.,                       &
                         lwetdep     = .TRUE.,                       &
                         density     = 1330._dp,                     &
                         iaerorad    = iradsoa,                      &
                         lwatsol     = .TRUE.,                       &  
                         kappa       = 0.06_dp,                       & !>>dod <<
                         lburden     = .TRUE.,                       &
                         idx         = soaprop(6)%spid_soa              ) 

        CALL new_species(nphase      = AEROSOL,                      &
                         longname    = 'Benzene SOA',                &
                         shortname   = 'SOA_BENZ',                   &
                         units       = 'kg kg-1',                    &
                         mw          = 98._dp,                       & 
                         tsubmname   = 'HAM',                        &
                         itrtype     = itrprog,                      &
                         ldrydep     = .TRUE.,                       &
                         lwetdep     = .TRUE.,                       &
                         density     = 1450._dp,                     &
                         iaerorad    = iradsoa,                      &
                         lwatsol     = .TRUE.,                       &  
                         kappa       = 0.06_dp,                       & !>>dod <<
                         lburden     = .TRUE.,                       &
                         idx         = soaprop(7)%spid_soa              ) 

     CASE(1)                               ! lump all anthropogenics into one distinct SOA species

        CALL new_species(nphase      = AEROSOL,                      &
                         longname    = 'Anthropogenic SOA',          &
                         shortname   = 'ASOA',                       &
                         units       = 'kg kg-1',                    &
                         mw          = 124._dp,                      & 
                         tsubmname   = 'HAM',                        &
                         itrtype     = itrprog,                      &
                         ldrydep     = .TRUE.,                       &
                         lwetdep     = .TRUE.,                       &
                         density     = 1450._dp,                     &
                         iaerorad    = iradsoa,                      &
                         lwatsol     = .TRUE.,                       &  
                         kappa       = 0.06_dp,                       &  !>>dod <<
                         lburden     = .TRUE.,                       &
                         idx         = soaprop(5)%spid_soa              ) 
        
     CASE(2)                              ! lump all anthropogenic SOA together with primary OC

        !! soaprop(5)%spid_soa = id_oc

     END SELECT

     !---SOA 2-product properties

     !---biogenics
     !   monoterpene SOA
     soaprop(1)%lvolatile = .TRUE.
         !soaprop(1)%Kp = 0.0637_dp                       ! Presto et al EST 2005
         ! soaprop(1)%tref = 295._dp                       ! Presto et al EST 2005
         !soaprop(1)%Kp = 0.038                                     ! Griffin et al JGR 1999
         !soaprop(1)%tref = 307._dp                                 ! Griffin et al JGR 1999
     soaprop(1)%Kp = 2.3_dp                                     ! Saathoff et al ACPD 2008
     soaprop(1)%tref = 293.3_dp                                 ! Saathoff et al ACPD 2008
     soaprop(1)%dH = 5.9E4_dp
     
     soaprop(2)%lvolatile = .TRUE.
         !soaprop(2)%Kp = 0.0026_dp                       ! Presto et al EST 2005
         !soaprop(2)%tref = 295._dp                       ! Presto et al EST 2005
         !soaprop(2)%Kp = 0.326                                     ! Griffin et al JGR 1999
         !soaprop(2)%tref = 307._dp                                  ! Griffin et al JGR 1999
     soaprop(2)%Kp = 0.028                                    ! Saathoff et al ACPD 2008
     soaprop(2)%tref = 293.3_dp                                  ! Saathoff et al ACPD 2008
     soaprop(2)%dH = 2.4E4_dp
     
     !   isoprene SOA
     soaprop(3)%lvolatile = .TRUE.
     soaprop(3)%Kp = 0.00862                                     ! Henze and Seinfeld GRL 2006
     soaprop(3)%tref = 295._dp
     soaprop(3)%dH = 4.2E4_dp

     soaprop(4)%lvolatile = .TRUE.
     soaprop(4)%Kp = 1.62_dp                                     ! Henze and Seinfeld GRL 2006
     soaprop(4)%tref = 295._dp
     soaprop(4)%dH = 4.2E4_dp


     !---anthropogenics according to lumping parameter
     SELECT CASE(nsoalumping) 
     CASE(0)                              ! no lumping, all anthropogenics distinct 
        !   toluene SOA
        soaprop(5)%lvolatile = .FALSE.
        soaprop(5)%Kp = 0._dp
        soaprop(5)%tref = 298._dp
        soaprop(5)%dH = 0._dp
        soaprop(5)%spid_tot = -1
        
        !   xylene SOA
        soaprop(6)%lvolatile = .FALSE.
        soaprop(6)%Kp = 0._dp
        soaprop(6)%tref = 298._dp
        soaprop(6)%dH = 0._dp
        soaprop(6)%spid_tot = -1

        !  benzene SOA
        soaprop(7)%lvolatile = .FALSE.
        soaprop(7)%Kp = 0._dp
        soaprop(7)%tref = 298._dp
        soaprop(7)%dH = 0._dp
        soaprop(7)%spid_tot = -1

        nsoaspec = 7

     CASE(1)                            ! lump all anthropogenics into one distinct SOA species

        soaprop(5)%lvolatile = .FALSE.
        soaprop(5)%Kp = 0._dp
        soaprop(5)%tref = 298._dp
        soaprop(5)%dH = 0._dp
        soaprop(5)%spid_tot = -1
        
        nsoaspec = 5

     CASE(2)
        nsoaspec = 4
     END SELECT

     !---set species id for unused array elements to 'undefined'
     IF (nsoaspec < nmaxsoa) THEN
        soaprop(nsoaspec+1:nmaxsoa)%spid_tot = -1
        soaprop(nsoaspec+1:nmaxsoa)%spid_soa = -1
     END IF

     !---count the number of non-volatile species (needed by soa2prod and soa_part
     !   see mo_ham_soa_processes).
     nnvol = 0
     DO i=1,nsoaspec
        IF (.NOT. soaprop(i)%lvolatile) nnvol = nnvol + 1
     END DO
     
     IF (nsoalumping == 2) nnvol = 1

   END SUBROUTINE soa_species
   
    !-----------------------------------------------------------------------------------------------

    SUBROUTINE start_soa_aerosols(nmod, lsoainclass)

      USE mo_ham,          ONLY: new_aerocomp
      USE mo_species,      ONLY: speclist

      IMPLICIT NONE

      INTEGER, INTENT(in)    :: nmod
      LOGICAL, INTENT(in)    :: lsoainclass(nmod)

      !---local data:
      INTEGER :: jn, jm, ispec 

      !---executable procedure

      DO jn = 1,nmod
         IF (lsoainclass(jn)) THEN

            DO jm=1,nsoaspec

               ispec = soaprop(jm)%spid_soa

               !---allocate memory for the aerosol phase tracers
               IF (.NOT. ALLOCATED(speclist(ispec)%iaerocomp) ) THEN
                 ALLOCATE(speclist(ispec)%iaerocomp(nmod))
                 speclist(ispec)%iaerocomp(:) = 0
               END IF

               speclist(ispec)%iaerocomp(jn) = new_aerocomp(jn, ispec, speclist(ispec)%itrtype)
            END DO
         END IF
      END DO

    END SUBROUTINE start_soa_aerosols
    !-----------------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------------
    SUBROUTINE construct_soa_streams

      ! adds input and diagnostic output streams for SOA
      ! called from: call_submodels

      ! Author:
      ! Declan O'Donnell, MPI-M, 2007

      !---inherited types, functions and data ---
      USE mo_memory_base,   ONLY: new_stream, default_stream_setting, add_stream_element, add_stream_reference
      USE mo_linked_list,   ONLY: t_stream, SURFACE
      USE mo_filename,      ONLY: trac_filetype

      IMPLICIT NONE

      !--- subroutine interface ---
      !    -
 
      !--- Local data ---
      TYPE(t_stream), POINTER   :: stream_soa

      !--- executable procedure ---
!++mgs ### Open a new soa stream (disentangling)
      CALL new_stream(stream_soa, 'soadiag', filetype=trac_filetype)

      !---diagnostics 
      CALL default_stream_setting (stream_soa, units='kg m-2 s-1',  lrerun=.TRUE.,  laccu=.TRUE.,  &
                                   lpost=.TRUE.,  leveltype=SURFACE )

      !>>dod redmine #113
      !--- Add standard fields for post-processing:
      CALL add_stream_reference (stream_soa, 'geosp'   ,'g3b',    lpost=.TRUE.)
      CALL add_stream_reference (stream_soa, 'lsp'     ,'sp',     lpost=.TRUE.)
      CALL add_stream_reference (stream_soa, 'aps'     ,'g3b',    lpost=.TRUE.)
      CALL add_stream_reference (stream_soa, 'gboxarea','geoloc', lpost=.TRUE.)
      !<<dod

      !   SOA excepted (for now) from emissions output switch

      !   SOA (pseudo-) chemistry
      CALL add_stream_element (stream_soa, 'D_PROD_SOA_MTERP',  d_prod_soa_mterp,  units='kg m-2 s-1',   &
                               longname='production of biogenic SOA from monoterpenes ')
      CALL add_stream_element (stream_soa, 'D_PROD_SOA_ISOP',  d_prod_soa_isop,  units='kg m-2 s-1',   &
                               longname='production of biogenic SOA from isoprene')
      CALL add_stream_element (stream_soa, 'D_PROD_SOA_TOL',  d_prod_soa_tol,  units='kg m-2 s-1',   &
                               longname='Production of anthropogenic SOA from toluene ')
      CALL add_stream_element (stream_soa, 'D_PROD_SOA_XYL',  d_prod_soa_xyl,  units='kg m-2 s-1',   &
                               longname='Production of anthropogenic SOA from xylene ')
      CALL add_stream_element (stream_soa, 'D_PROD_SOA_BENZ',  d_prod_soa_benz,  units='kg m-2 s-1', &
                               longname='Production of anthropogenic SOA from benzene ')

      CALL add_stream_element (stream_soa, 'D_CHEM_SINK_MTERP', d_chem_sink_mterp, units='kg m-2 s-1',   &
                               longname='SOA precursor destruction - monoterpenes ')
      CALL add_stream_element (stream_soa, 'D_CHEM_SINK_ISOP',  d_chem_sink_isop,  units='kg m-2 s-1',   &
                               longname='SOA precursor destruction - isoprene')
      CALL add_stream_element (stream_soa, 'D_CHEM_SINK_TOL',  d_chem_sink_tol,  units='kg m-2 s-1',   &
                               longname='SOA precursor destruction -  toluene ')
      CALL add_stream_element (stream_soa, 'D_CHEM_SINK_XYL',  d_chem_sink_xyl,  units='kg m-2 s-1',   &
                               longname='SOA precursor destruction -  xylene ')
      CALL add_stream_element (stream_soa, 'D_CHEM_SINK_BENZ',  d_chem_sink_benz,  units='kg(C) m-2 s-1',   &
                               longname='SOA precursor destruction -  benzene ')
      
    END SUBROUTINE construct_soa_streams

    !--------------------------------------------------------------------------------------------------

    SUBROUTINE set_soa_tracer_attr(kspid, init, itran, iwrite, iconv, ivdiff, iint, iwdep, iddep)

      USE mo_tracdef,           ONLY: ON, OFF, &
                                      RESTART, CONSTANT, INITIAL
      USE mo_advection,         ONLY: iadvec
      USE mo_ham,               ONLY: nsoaspec, nwetdep
      USE mo_exception,         ONLY: finish 
      USE mo_species,           ONLY: speclist

      IMPLICIT NONE
      
      !---subroutine interface
      INTEGER, INTENT(IN)  :: kspid                  ! species index
      INTEGER, INTENT(OUT) :: init                   ! tracer initialisation method
      INTEGER, INTENT(OUT) :: itran                  ! advection of the tracer
      INTEGER, INTENT(OUT) :: iwrite                 ! output of the tracer 
      INTEGER, INTENT(OUT) :: iconv                  ! convective transport of the tracer
      INTEGER, INTENT(OUT) :: ivdiff                 ! vertical diffusion of the tracer
      INTEGER, INTENT(OUT) :: iint                   ! integration flag
      INTEGER, INTENT(OUT) :: iwdep                  ! wet deposition flag
      INTEGER, INTENT(OUT) :: iddep                  ! dry deposition flag

      !---local data
      INTEGER :: jm
      LOGICAL :: lfound

      !--executable procedure
      lfound = .FALSE.
      DO jm=1,nsoaspec
         IF (kspid == soaprop(jm)%spid_tot) THEN
            lfound = .TRUE.
            init     = RESTART+CONSTANT+INITIAL
            itran    = iadvec  
            iwrite   = OFF
            iconv    = OFF
            ivdiff   = OFF
            iint     = ON
            iwdep    = OFF
            iddep    = OFF
         ELSE IF (kspid == soaprop(jm)%spid_soa) THEN
            lfound = .TRUE.
            IF (soaprop(jm)%lvolatile) THEN
               init     = -1
               itran    = OFF
               iwrite   = ON
               iconv    = ON
               ivdiff   = ON
               iint     = OFF
               iwdep    = nwetdep
               iddep    = 2
            ELSE
               init     = RESTART+CONSTANT+INITIAL
               itran    = iadvec
               iwrite   = ON
               iconv    = ON
               ivdiff   = ON
               iint     = ON
               iwdep    = nwetdep
               iddep    = 2
            END IF
         END IF
      END DO

      IF (.NOT. lfound) CALL finish('set_soa_tracer_attr', 'subroutine called for non-SOA species '//speclist(kspid)%shortname)
    END SUBROUTINE set_soa_tracer_attr

END MODULE mo_ham_soa
