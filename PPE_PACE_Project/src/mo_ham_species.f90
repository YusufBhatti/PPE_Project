!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_species.f90
!!
!! \brief
!! mo_ham_species assigns identities and properties to species in the HAM model
!!
!! \author Declan O'Donnell (MPI-Met)
!!
!! \responsible_coder
!! Declan O'Donnell, declan.Odonnell@fmi.fi
!!
!! \revision_history
!!   -# Declan O'Donnell (MPI-Met) - original code (2008)
!!   -# K. Zhang (MPI-Met) - seperate species list for ham; new_species (2009-07) 
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

MODULE mo_ham_species
  
  !--- inherited types, data and functions
  USE mo_kind,       ONLY: dp
  USE mo_tracdef,    ONLY: jptrac 
  USE mo_species,    ONLY: t_species 
  

  IMPLICIT NONE

  !--- public member functions
  
  PUBLIC :: ham_species 

  !--- public module data

  
  !   Maximum number of species in the model
  INTEGER, PARAMETER, PUBLIC :: nmaxspec = 100 

  !   tracer request types
  INTEGER, PARAMETER, PUBLIC :: inotr = 0       ! request no tracer   
  INTEGER, PARAMETER, PUBLIC :: idiagtr = 1     ! request a diagnostic (not transported) tracer
  INTEGER, PARAMETER, PUBLIC :: iprogtr = 2     ! request a prognostic (transported) tracer

  ! ### keys to refractive index tables now in mo_ham_rad_data

  !   basic HAM model, gas phase compounds
  INTEGER, PUBLIC :: id_dms,                  & ! Dimethyl sulphide
                     id_so2,                  & ! Sulphur dioxide
                     id_so4g,                 & ! Sulphate (gas)
                     id_oh,                   & ! Hydroxyl radical
                     id_h2o2,                 & ! Hydrogen Peroxide
                     id_o3,                   & ! Ozone
                     id_no2,                  & ! Nitrogen dioxide
!gf see #146
                     id_no3,                  & ! Nitrate radical
!gf
                     id_ocnv                   ! Organic carbon (gas)

  !   basic HAM model, aerosol phase compounds
  INTEGER, PUBLIC :: id_so4,                  & ! Sulphate (aerosol)
                     id_bc,                   & ! Black carbon
                     id_oc,                   & ! Organic carbon 
                     id_ss,                   & ! Sea Salt
                     id_du,                   & ! Dust
                     id_wat                     ! Aerosol water

  ! Number of model species
  INTEGER, PUBLIC :: ham_nspec      ! number of species in HAM
  INTEGER, PUBLIC :: ham_naerospec  ! number of aerosol species in HAM

  ! Model species. These arrays are overdimensioned: the sum of defined species 
  ! (gas plus aerosol) shall not exceed nmaxspec
  TYPE(t_species), PUBLIC, ALLOCATABLE, TARGET :: ham_aerospec(:)
  
  
  
  CONTAINS




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Register ham gas and aerosol species (except for SOA) 
!!
!! @author 
!! D. O'Dannel (MPI-Met) 
!! K. Zhang    (MPI-Met) 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! D. O'Dannel (MPI-Met) - original version - (2008-??) 
!! K. Zhang    (MPI-Met) - seperate species list for ham - (2009-07) 
!!
!! @par This subroutine is called by
!! to-be-filled
!!
!! @par Externals:
!! <ol>
!! <li>new_species 
!! </ol>
!!
!! @par Notes 
!!  
!! @par Responsible coder
!! kai.zhang@zmaw.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ham_species

  !---inherited types, data and functions
  USE mo_ham,             ONLY: nham_subm, HAM_SALSA, HAM_M7
  USE mo_submodel,        ONLY: lham
  USE mo_tracdef,         ONLY: GAS, AEROSOL, GAS_OR_AEROSOL, itrprog, itrpresc, itrnone
  USE mo_species,         ONLY: new_species
  USE mo_ham_rad_data,    ONLY: iradso4, iradbc, iradoc, iradss, iraddu, iradwat
!gf
  USE mo_ham,             ONLY: mw_so2, mw_so4, mw_dms
!gf

  IMPLICIT NONE

  !--- executable procedure ----

  !--------- 0. Initialisations    

  !--- instantiate the model species. start with the gas phase
  IF (lham) THEN

     !--------- 1. Dimethyl Sulphide (DMS)
     ! 

     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Dimethyl sulphide', &
                      shortname   = 'DMS',               &
                      units       = 'kg kg-1',           &
!gf                      mw          = 62.019_dp,           &
                      mw          = mw_dms,              &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrprog,             &
                      ldrydep     = .FALSE.,             &
                      lwetdep     = .FALSE.,             &
                      henry       = (/ 0.54_dp, 3460._dp /),   & !csld(#275)
                      idx         = id_dms                  )

     !--------- 2. Sulphur Dioxide 
     !

     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Sulphur Dioxide',   &
                      shortname   = 'SO2',               &
                      units       = 'kg kg-1',           &
!gf                      mw          = 63.692_dp,           &
                      mw          = mw_so2,              &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrprog,             &
                      ldrydep     = .TRUE.,              &
                      lwetdep     = .TRUE.,              &
                      dryreac     = 0._dp,               &
                      henry       = (/ 1.36_dp, 4250._dp /), & !csld(#275)
                      idx         = id_so2                  )


     !--------- 4. Hydroxyl Radical 
     !             

     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Hydroxyl radical',  &
                      shortname   = 'OH',                &
                      units       = 'VMR',               &
                      mw          = 17.003_dp,           &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrnone,             & !gf #57
                      ldrydep     = .FALSE.,             &
                      lwetdep     = .FALSE.,             &
                      dryreac     = 0._dp,               &
                      henry       = (/ 39._dp, 0._dp /), & !csld(#275)
                      idx         = id_oh                   )
  
     !--------- 5. Hydrogen Peroxide
     !

     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Hydrogen peroxide', &
                      shortname   = 'H2O2',              &
                      units       = 'VMR',               &
                      mw          = 34.005_dp,           &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrnone,             & !gf #57
                      ldrydep     = .FALSE.,             &
                      lwetdep     = .FALSE.,             &
                      dryreac     = 1._dp,               & 
                      henry       = (/ 8.44E4_dp, 7600._dp /), & !csld(#275)
                      idx         = id_h2o2                 )

     !--------- 6. Ozone
     !

     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Ozone',             &
                      shortname   = 'O3',                &
                      units       = 'VMR',               &
                      mw          = 48.0_dp,             &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrnone,             & !gf #57
                      ldrydep     = .FALSE.,             &
                      lwetdep     = .FALSE.,             &
                      dryreac     = 1._dp,               & 
                      henry       = (/ 1.03E-2_dp, 2830._dp /), & !csld(#275)
                      idx         = id_o3                   )

     !--------- 7. Nitrogen Dioxide
     !

     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Nitrogen dioxide',  &
                      shortname   = 'NO2',               &
                      units       = 'VMR',               &
                      mw          = 46.006_dp,           &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrnone,             & !gf #57
                      ldrydep     = .FALSE.,             &
                      lwetdep     = .FALSE.,             &
                      dryreac     = 1._dp,               & 
                      henry       = (/ 1.2E-2_dp, 2360._dp /), & !csld(#275)
                      idx         = id_no2                  )
!gf see #146

     !--------- 8. Nitrate radical
     !

     CALL new_species(nphase      = GAS,                 &
                      longname    = 'Nitrate radical',   &
                      shortname   = 'NO3',               &
                      units       = 'VMR',               &
                      mw          = 62.049_dp,           &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrnone,             & !gf #57
                      ldrydep     = .FALSE.,             &
                      henry       = (/ 3.8E-2_dp, 0._dp /), & !csld(#275)
                      lwetdep     = .FALSE.,             &
                      idx         = id_no3               )

!gf

     !--------- 3.a Gas phase sulphate
     !  (sulphate can be present in both the gas-phase and the aerosol phase)
     !  density here is that of H2SO4 
!>>DT
     CALL new_species(nphase      = GAS,      &
                      longname    = 'Sulphuric acid',    &
                      shortname   = 'H2SO4',             &
                      units       = 'kg kg-1',           &
                      mw          = mw_so4,              &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrprog,             &
                      ldrydep     = .TRUE.,              &
                      lwetdep     = .TRUE.,              &
                      henry       = (/ 1.3e15_dp, 20000._dp /),  & !SF #571
                      dryreac     = 0._dp,               &
                      idx          = id_so4g             )
!<<DT

     !--------- end of basic model gas phase species

     !--------- basic model aerosol species

     !--------- 3.b Aerosol phase sulphate
     !  (sulphate can be present in both the gas-phase and the aerosol phase)
     !  density here is that of H2SO4 

!>>DT
     CALL new_species(nphase      = AEROSOL,             &
                      longname    = 'Sulphate',          &
                      shortname   = 'SO4',               &
                      units       = 'kg kg-1',           &
                      mw          = mw_so4,              &
                      tsubmname   = 'HAM',               &
                      itrtype     = itrprog,             &
                      ldrydep     = .TRUE.,              &
                      lwetdep     = .TRUE.,              &
                      density      = 1841._dp,           &
                      iaerorad     = iradso4,            &
                      lwatsol      = .TRUE.,             & 
                      lelectrolyte = .TRUE.,             &
                      nion         = 2,                  &
                      osm          = 1._dp,              &
                      kappa        = 0.60_dp,            &
                      lburden      = .TRUE.,             &
                      idx          = id_so4              )
!<<DT

     IF(nham_subm == HAM_SALSA) THEN

        !--------- 10a. Organic Carbon (Primary organic aerosol)
        !
        
        CALL new_species(nphase       = GAS_OR_AEROSOL,      &
                         longname     = 'Organic carbon',    &
                         shortname    = 'OC',                &
                         units        = 'kg kg-1',           &
                         mw           = 180._dp,             &
                         tsubmname   = 'HAM',                &
                         itrtype      = itrprog,             &
                         density      = 2000._dp,            &
                         iaerorad     = iradoc,              &
                         lwatsol      = .TRUE.,              &  
                         kappa        = 0.06_dp,             &
                         ldrydep      = .TRUE.,              &
                         lwetdep      = .TRUE.,              &   
                         henry       = (/1.E5_dp, 0._dp/),   &
                         idx          = id_ocnv                   )

        ! equate "aerosol" species id with "gas" species id for OCNV
        id_oc = id_ocnv
        
     END IF

     !--------- 9. Black Carbon
     !

     CALL new_species(nphase       = AEROSOL,             &
                      longname     = 'Black carbon',      &
                      shortname    = 'BC',                &
                      units        = 'kg kg-1',           &
                      mw           = 12.010_dp,           &
                      tsubmname    = 'HAM',               &
                      itrtype      = itrprog,             &
                      density      = 2000._dp,            &
                      iaerorad     = iradbc,              &
                      lwatsol      = .FALSE.,             &  
                      ldrydep      = .TRUE.,              &
                      lwetdep      = .TRUE.,              &   
                      idx          = id_bc                   )

     IF(nham_subm == HAM_M7) THEN

        !--------- 10b. Organic Carbon (Primary organic aerosol)
        !
        
        CALL new_species(nphase       = AEROSOL,             &
                         longname     = 'Organic carbon',    &
                         shortname    = 'OC',                &
                         units        = 'kg kg-1',           &
                         mw           = 180._dp,             &
                         tsubmname   = 'HAM',                &
                         itrtype      = itrprog,             &
                         density      = 2000._dp,            &
                         iaerorad     = iradoc,              &
                         lwatsol      = .TRUE.,              &  
                         kappa        = 0.06_dp,             &
                         ldrydep      = .TRUE.,              &
                         lwetdep      = .TRUE.,              &   
                         idx          = id_oc                   )
        
     END IF

     !--------- 11. Sea salt
     !

     CALL new_species(nphase       = AEROSOL,             &
                      longname     = 'Sea salt',          &
                      shortname    = 'SS',                &
                      units        = 'kg kg-1',           &
                      mw           = 58.443_dp,           &
                      tsubmname    = 'HAM',               &
                      itrtype      = itrprog,             &
                      density      = 2165._dp,            &
                      iaerorad     = iradss,              &
                      lwatsol      = .TRUE.,              & 
                      lelectrolyte = .TRUE.,              &
                      nion         = 2,                   &
                      osm          = 1._dp,               & 
                      kappa        = 1._dp,               &    !>>dod<<
                      ldrydep      = .TRUE.,              &
                      lwetdep      = .TRUE.,              &   
                      idx          = id_ss                   )


     !--------- 12. Mineral Dust
     !

     CALL new_species(nphase       = AEROSOL,             &
                      longname     = 'Dust',              &
                      shortname    = 'DU',                &
                      units        = 'kg kg-1',           &
                      mw           = 250._dp,             &
                      tsubmname    = 'HAM',               &
                      itrtype      = itrprog,             &
                      density      = 2650._dp,            &
                      iaerorad     = iraddu,              &
                      lwatsol      = .FALSE.,             &  
                      ldrydep      = .TRUE.,              &
                      lwetdep      = .TRUE.,              &   
                      idx          = id_du                   )


     !--------- 13. Water (includes only water on aerosols, not water vapour, cloud water or cloud ice)
     !                    

     CALL new_species(nphase       = AEROSOL,             &
                      longname     = 'Aerosol water',     &
                      shortname    = 'WAT',               &
                      units        = 'kg kg-1',           &
                      mw           = 18.0_dp,             &
                      tsubmname    = 'HAM',               &
                      itrtype      = itrprog,             &
                      density      = 1000._dp,            &
                      iaerorad     = iradwat,             &
                      lwatsol      = .TRUE.,              &  
                      lemis        = .FALSE.,             &
                      ldrydep      = .FALSE.,             &
                      lwetdep      = .FALSE.,             &   
                      idx          = id_wat                  )


     !--------- end of basic model aerosol phase species
  END IF
  END SUBROUTINE ham_species

END MODULE mo_ham_species

