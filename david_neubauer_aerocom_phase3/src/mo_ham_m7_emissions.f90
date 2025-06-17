!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_m7_emissions.f90
!!
!! \brief
!! Module for HAM-M7 specific emissions that cannot be dealt with the standard
!! processing scheme.
!!
!! \author M.G. Schultz and S. Schroeder, (FZ Juelich)
!! \author D. O'Donnell, ETH-Z (ETH Zurich)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M.G. Schultz and S. Schroeder (FZ Juelich) - original code (2010-02-11)
!!   -# D. O'Donnell (ETH-Zurich) - added AEROCOM emissions (2010-04-21)
!!
!! \limitations
!! None
!!
!! \details
!!  Routine ham_m7_emissions is called from emi_interface and must return
!!  lprocessed=.TRUE. when HAM decided to do something for this sector and 
!!  species. In this case, ham_m7_emissions must set the boundary conditions
!!  for emission mass fluxes.
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

MODULE mo_ham_m7_emissions

  USE mo_kind,             ONLY: dp
  USE mo_ham,              ONLY: mw_so2, mw_so4, mw_s, sigma_fine, sigma_coarse

  IMPLICIT NONE

  PRIVATE

  ! public variables  (see declaration below)

  ! subprograms
  PUBLIC                       :: ham_m7_init_emissions
  PUBLIC                       :: ham_m7_init_emi_stream
  PUBLIC                       :: ham_m7_emissions

  !### temporary: probably obsolete when bgc_dust and megan use boundary conditions for input fields


  ! sector indices 
  ! IPCC emissions have separate sectors for various fossil and fire emissions
  ! for AEROCOM emissions use fossil and fire
  INTEGER           :: idsec_seasalt
  INTEGER           :: idsec_fire, idsec_ffire, idsec_gfire, idsec_awb
  INTEGER           :: idsec_fossil, idsec_dom, idsec_ene, idsec_ind, idsec_tra, idsec_wst
  INTEGER           :: idsec_ships
  INTEGER           :: idsec_biofuel             ! aerocom sector

  ! indices for weighting factors to distribute mass flux among modes
  INTEGER           :: idx_nbcki, idx_mbcki      ! black carbon number and mass
  INTEGER           :: idx_nocki, idx_nocks, &   ! organic carbon number and mass
                       idx_mocki, idx_mocks, idx_mocas
  INTEGER           :: idx_ns4ks, idx_ns4as, idx_ns4cs, &   ! sulfate
                       idx_ms4ks, idx_ms4as, idx_ms4cs
  INTEGER           :: idx_nduai, idx_nduci, &   ! dust
                       idx_nduas, idx_nducs, &
                       idx_mduai, idx_mduci, &
                       idx_mduas, idx_mducs
  INTEGER           :: idx_nssas, idx_nsscs, &   ! seasalt
                       idx_mssas, idx_msscs

  ! diagnostic pointers
  REAL(dp), POINTER :: d_emi_mduai(:,:), d_emi_nduai(:,:), d_emi_mduci(:,:), d_emi_nduci(:,:)
  REAL(dp), POINTER :: d_emi_mssas(:,:), d_emi_nssas(:,:), d_emi_msscs(:,:), d_emi_nsscs(:,:)

  ! Module parameters to distribute emission mass flux across modes and to convert
  ! mass flux in number flux
  REAL(dp), PARAMETER:: zbb_wsoc_perc  = 0.65_dp,       &! Biom. Burn. Percentage of Water Soluble OC (WS
                                                         ! (M.O. Andreae; Talk: Smoke and Climate)
                        zbg_wsoc_perc  = 0.65_dp,       &! Assume same Percentage of WSOC for biogenic OC
                        zom2oc         = 1.4_dp,        &! Mass ratio organic species to organic carbon
                                                         ! (Seinfeld and Pandis, 1998, p709;
                                                         !  Ferek et al., JGR, 1998)
                        zfacso2        = 0.975_dp,      &! factor to scale primary SO4 emissions
                                                         ! AEROCOM assumption 2.5 % of the SO2 emissions
                                                         ! in the from of SO4
                        zso2ts         = mw_s/mw_so2,   &! conversion factor SO2 to S
                        zso2tso4       = mw_so4/mw_so2, &! conversion mass SO2 to mass SO4
                        ref_sigma_dust_ai = sigma_fine,    &! source standard deviations
                        ref_sigma_dust_ci = sigma_coarse

  !---variables for mapping emission masses into aerosol numbers
  REAL(dp)           :: zm2n_bcki_ff, zm2n_bcki_bb, &
                        zm2n_bcks_bb, zm2n_ocki_ff, &
                        zm2n_ocki_bb, zm2n_ocki_bg, &
                        zm2n_ocks_bb, zm2n_ocks_bg, &
                        zm2n_s4ks_ff, zm2n_s4as_ff, &
                        zm2n_s4cs_ff,               &
                        zm2n_s4ks_bb, zm2n_s4as_bb, &
                        zm2n_s4cs_bb,               &
                        zm2n_duai, zm2n_duci

  CONTAINS

  !! ---------------------------------------------------------------------------
  !! subroutine to initialize HAM specific emissions: must set appropriate local sector
  !! indices idsec_xyz

  SUBROUTINE ham_m7_init_emissions(nsectors)

  USE mo_emi_matrix,             ONLY: em_get_sector_info
  USE mo_util_string,            ONLY: toupper
  USE mo_math_constants,         ONLY: pi
  USE mo_species,                ONLY: speclist, spec_idt, spec_ntrac
  USE mo_ham,                    ONLY: sizeclass, idsec_dust, idsec_biogenic, & 
                                       ibc_dust
  USE mo_ham_species,            ONLY: id_bc, id_oc, id_so4, id_du, id_ss
  USE mo_ham_m7_trac,            ONLY: idt_mbcki, idt_mocki, idt_mocks, idt_mocas,  &
                                       idt_ms4ks, idt_ms4as, idt_ms4cs, &
                                       idt_mduai, idt_mduci, &
                                       idt_mduas, idt_mducs, &
                                       idt_mssas, idt_msscs
  USE mo_ham_m7ctl,              ONLY: cmr2ram, iaiti, iaits, iaccs, icoas, iacci, &
                                       icoai
  USE mo_ham_dust,               ONLY: bgc_dust_initialize
  USE mo_boundary_condition,       ONLY: bc_query
  USE mo_external_field_processor, ONLY: EF_MODULE

  INTEGER, INTENT(in)      :: nsectors    ! number of sectors defined in emi matrix

  INTEGER           :: i, nvars, jt, ieftype
  CHARACTER(LEN=64) :: secname

  !-- parameters
  ! ### Note: these parameters should be further cleaned up -- there are some ad-hoc 
  ! assumptions in the code below (50% mode split) which are not entirely style-conform
  ! with the cmr_ parameters
  REAL(dp), PARAMETER:: cmr_ff         = 0.03E-6_dp,    &! Fossil fuel emissions:
                                                         ! assumed number median radius of the emitted
                                                         ! particles with the standard deviation given
                                                         ! in mo_ham_m7ctl [m]. Has to lie within the
                                                         ! Aitken mode for the current setup!
                        cmr_bb         = 0.075E-6_dp,   &! Biomass burning emissions:
                                                         ! Assumed number median radius of the emitted
                                                         ! particles with the standard deviation given
                                                         ! in mo_ham_m7ctl [m]. Has to lie within the
                                                         ! Aitken mode for the current setup!
                        cmr_bg         = 0.03E-6_dp,    &! Biogenic secondary particle formation:
                                                         ! Assumed number median radius of the emitted
                                                         ! particles with the standard deviation given
                                                         ! in mo_ham_m7ctl [m]. Has to lie within the
                                                         ! Aitken mode for the current setup!
                        cmr_sk         = 0.03E-6_dp,    &! SO4 primary emission  ---> aitken mode
                                                         ! Assumed number median radius of the emitted
                                                         ! particles with the standard deviation given
                                                         ! in mo_ham_m7ctl [m]. Has to lie within the
                                                         ! Aitken mode for the current setup!
                        cmr_sa         = 0.075E-6_dp,   &! SO4 primary emission  ---> accumulation mode
                                                         ! Assumed number median radius of the emitted
                                                         ! particles with the standard deviation given
                                                         ! in mo_ham_m7ctl [m]. Has to lie within the
                                                         ! Accumulation mode for the current setup!
                        cmr_sc         = 0.75E-6_dp,    &! SO4 primary emission  ---> coarse mode
                                                         ! Assumed number median radius of the emitted
                                                         ! particles with the standard deviation given
                                                         ! in mo_ham_m7ctl [m]. Has to lie within the
                                                         ! Coarse mode for the current setup!
                        mmr_dust_ai = 0.35E-6_dp,       &! source mass median radii [m]
                        mmr_dust_ci = 1.75E-6_dp
  REAL(dp) :: densdust

  !-- initialize
  idsec_dust      = -1
  idsec_seasalt   = -1
  idsec_fire      = -1    ! fire emissions all in one (AEROCOM)
  idsec_ffire     = -1    ! forest fire emissions (IPCC)
  idsec_gfire     = -1    ! grass fire emissions (IPCC)
  idsec_awb       = -1    ! agricultural waste burning (AEROCOM and IPCC)
  idsec_fossil    = -1    ! fossil fuel emissions all in one (AEROCOM)
  idsec_dom       = -1    ! domestic (fossil and biofuel) emissions (IPCC)
  idsec_ene       = -1    ! energy sector emissions (IPCC)
  idsec_ind       = -1    ! industry sector emissions (IPCC)
  idsec_tra       = -1    ! traffic emissions (IPCC)
  idsec_wst       = -1    ! waste handling emissions (IPCC)
  idsec_ships     = -1    ! international ship traffic emissions (IPCC)
  idsec_biogenic  = -1
  idsec_biofuel   = -1    ! biofuel emissions (AEROCOM)

  DO i = 1, nsectors
    CALL em_get_sector_info(i, secname, nvars)
    IF ((TRIM(toupper(secname)) == 'DUST') .AND. (nvars > 0))     idsec_dust = i
    IF ((TRIM(toupper(secname)) == 'SEASALT') .AND. (nvars > 0))  idsec_seasalt = i
    IF ((TRIM(toupper(secname)) == 'FIRE') .AND. (nvars > 0))     idsec_fire = i
    IF ((TRIM(toupper(secname)) == 'FFIRE') .AND. (nvars > 0))    idsec_ffire = i
    IF ((TRIM(toupper(secname)) == 'GFIRE') .AND. (nvars > 0))    idsec_gfire = i
    IF ((TRIM(toupper(secname)) == 'AWB') .AND. (nvars > 0))      idsec_awb = i
    IF ((TRIM(toupper(secname)) == 'FOSSIL') .AND. (nvars > 0))   idsec_fossil = i
    IF ((TRIM(toupper(secname)) == 'DOM') .AND. (nvars > 0))      idsec_dom = i
    IF ((TRIM(toupper(secname)) == 'ENE') .AND. (nvars > 0))      idsec_ene = i
    IF ((TRIM(toupper(secname)) == 'IND') .AND. (nvars > 0))      idsec_ind = i
    IF ((TRIM(toupper(secname)) == 'TRA') .AND. (nvars > 0))      idsec_tra = i
    IF ((TRIM(toupper(secname)) == 'WST') .AND. (nvars > 0))      idsec_wst = i
    IF ((TRIM(toupper(secname)) == 'SHIPS') .AND. (nvars > 0))    idsec_ships = i
    IF ((TRIM(toupper(secname)) == 'BIOGENIC') .AND. (nvars > 0)) idsec_biogenic = i
    IF ((TRIM(toupper(secname)) == 'BIOFUEL')  .AND. (nvars > 0)) idsec_biofuel = i
  END DO

  !-- set mass to number conversion factors for HAM
  !   Calculate factors to convert mass flux in number flux for
  !   given number median radii (cmr) and standard deviation
  !   (implicitly by the conversion factor cmr2ram) of the modes
  !
  !    N = M/m = M/(4/3 * pi * dens * R(averageMass)**3)
  !            = M * (3/(4*pi*dens*R(averageMass)))
  !            !
  !            = M * zm2n_xx

  ! sschr:
  ! inucs -> ns
  ! iaits -> ks
  ! iaccs -> as
  ! icoas -> cs
  ! iaiti -> ki
  ! iacci -> ai
  ! icoai -> ci

  zm2n_bcki_ff=3._dp/(4._dp*pi*speclist(id_bc)%density*(cmr_ff*cmr2ram(iaiti))**3)
  zm2n_bcki_bb=3._dp/(4._dp*pi*speclist(id_bc)%density*(cmr_bb*cmr2ram(iaiti))**3)
  zm2n_bcks_bb=3._dp/(4._dp*pi*speclist(id_bc)%density*(cmr_bb*cmr2ram(iaits))**3)

  zm2n_ocki_ff=3._dp/(4._dp*pi*speclist(id_oc)%density*(cmr_ff*cmr2ram(iaiti))**3)
  zm2n_ocki_bb=3._dp/(4._dp*pi*speclist(id_oc)%density*(cmr_bb*cmr2ram(iaiti))**3)
  zm2n_ocki_bg=3._dp/(4._dp*pi*speclist(id_oc)%density*(cmr_bg*cmr2ram(iaiti))**3)
  zm2n_ocks_bb=3._dp/(4._dp*pi*speclist(id_oc)%density*(cmr_bb*cmr2ram(iaits))**3)
  zm2n_ocks_bg=3._dp/(4._dp*pi*speclist(id_oc)%density*(cmr_bg*cmr2ram(iaits))**3)

!>>gf - Redmine #267
  zm2n_s4ks_ff=3._dp/((4._dp*pi*speclist(id_so4)%density*(cmr_sk*cmr2ram(iaits))**3))
  zm2n_s4as_ff=3._dp/((4._dp*pi*speclist(id_so4)%density*(cmr_sa*cmr2ram(iaccs))**3))
  zm2n_s4cs_ff=3._dp/((4._dp*pi*speclist(id_so4)%density*(cmr_sc*cmr2ram(icoas))**3))
  zm2n_s4ks_bb=3._dp/((4._dp*pi*speclist(id_so4)%density*(cmr_sk*cmr2ram(iaits))**3))
  zm2n_s4as_bb=3._dp/((4._dp*pi*speclist(id_so4)%density*(cmr_sa*cmr2ram(iaccs))**3))
  zm2n_s4cs_bb=3._dp/((4._dp*pi*speclist(id_so4)%density*(cmr_sc*cmr2ram(icoas))**3))
!<<gf

  IF (idsec_dust > 0) THEN
    densdust=speclist(id_du)%density
    !--- Calculate conversion factor from mass flux to number flux for log-normal distributions:
    zm2n_duai=3._dp/(4._dp*pi*densdust*mmr_dust_ai**3) * EXP(4.5_dp*LOG(ref_sigma_dust_ai)**2)
    zm2n_duci=3._dp/(4._dp*pi*densdust*mmr_dust_ci**3) * EXP(4.5_dp*LOG(ref_sigma_dust_ci)**2)
  !-- call submodel specific initialisations
    CALL bc_query(ibc_dust, ef_type=ieftype)
    !-- call submodel specific initialisations
    ! (nothing to be done if those are replaced by compound-specific file input)
    IF (ieftype == EF_MODULE) CALL bgc_dust_initialize 
!spec_ntrac(id_du) is set to four (DU_AS, DU_CS, DU_AI, DU_CI)!
!output of number tracers is done by BYNUMMODE
    nvars = spec_ntrac(id_du)+4
    spec_idt(id_du, nvars-3) = sizeclass(iaccs)%idt_no ! iaccs -> as
    spec_idt(id_du, nvars-2) = sizeclass(icoas)%idt_no ! icoas -> cs
    spec_idt(id_du, nvars-1) = sizeclass(iacci)%idt_no ! iacci -> ai
    spec_idt(id_du, nvars)   = sizeclass(icoai)%idt_no   ! icoai -> ci
    spec_ntrac(id_du) = nvars
    idx_nduas = nvars-3
    idx_nducs = nvars-2       
    idx_nduai = nvars-1
    idx_nduci = nvars         
    idx_mduas = -1
    idx_mducs = -1
    idx_mduai = -1
    idx_mduci = -1
    DO jt=1,nvars-4
      IF (spec_idt(id_du, jt) == idt_mduas) idx_mduas = jt
      IF (spec_idt(id_du, jt) == idt_mducs) idx_mducs = jt
      IF (spec_idt(id_du, jt) == idt_mduai) idx_mduai = jt
      IF (spec_idt(id_du, jt) == idt_mduci) idx_mduci = jt
    END DO
  END IF
    
  IF (idsec_seasalt > 0) THEN
    nvars = spec_ntrac(id_ss)+2
    spec_idt(id_ss, nvars-1) = sizeclass(iaccs)%idt_no ! iaccs -> as
    spec_idt(id_ss, nvars)   = sizeclass(icoas)%idt_no   ! icoas -> cs
    spec_ntrac(id_ss) = nvars
    idx_nssas = nvars-1
    idx_nsscs = nvars         
    idx_mssas = -1
    idx_msscs = -1
    DO jt=1,nvars-2
      IF (spec_idt(id_ss, jt) == idt_mssas) idx_mssas = jt
      IF (spec_idt(id_ss, jt) == idt_msscs) idx_msscs = jt
    END DO
  END IF

  !-- expand spec_idt list to include aerosol number tracers and find appropriate aerosol mass tracers
  !   for the distribution of emission among aerosol modes
  !-- black carbon
  nvars = spec_ntrac(id_bc)+1
  spec_idt(id_bc, nvars) = sizeclass(iaiti)%idt_no      ! aitken mode insoluble (former idt_nki)
  spec_ntrac(id_bc) = nvars
  idx_nbcki = nvars                                  ! index to weight array for aerosol number BC
  idx_mbcki = -1
  DO jt=1,nvars-1
    IF (spec_idt(id_bc, jt) == idt_mbcki) idx_mbcki = jt
  END DO

  !-- organic carbon
  nvars = spec_ntrac(id_oc)+2
  spec_idt(id_oc, nvars-1) = sizeclass(iaiti)%idt_no    ! aitken mode insoluble (former idt_nki)
  spec_idt(id_oc, nvars)   = sizeclass(iaits)%idt_no    ! aitken mode soluble (former idt_nks)
  spec_ntrac(id_oc) = nvars
  idx_nocki = nvars-1                                ! index to weight array for aerosol number OC
  idx_nocks = nvars         
  idx_mocki = -1
  idx_mocks = -1
  idx_mocas = -1
  DO jt=1,nvars-2
    IF (spec_idt(id_oc, jt) == idt_mocki) idx_mocki = jt
    IF (spec_idt(id_oc, jt) == idt_mocks) idx_mocks = jt
    IF (spec_idt(id_oc, jt) == idt_mocas) idx_mocas = jt
  END DO

  !-- sulfur
  nvars = spec_ntrac(id_so4)+3
  spec_idt(id_so4, nvars-2) = sizeclass(iaits)%idt_no  ! aitken mode soluble (former idt_nks)
  spec_idt(id_so4, nvars-1) = sizeclass(iaccs)%idt_no  ! accumulation mode soluble (former idt_nas)
  spec_idt(id_so4, nvars)   = sizeclass(icoas)%idt_no  ! coarse mode soluble (former idt_ncs)
  spec_ntrac(id_so4) = nvars
  idx_ns4ks = nvars-2                               ! index to weight array for aerosol number so4
  idx_ns4as = nvars-1   
  idx_ns4cs = nvars         
  idx_ms4ks = -1
  idx_ms4as = -1
  idx_ms4cs = -1
  DO jt=1,nvars-3
    IF (spec_idt(id_so4, jt) == idt_ms4ks) idx_ms4ks = jt
    IF (spec_idt(id_so4, jt) == idt_ms4as) idx_ms4as = jt
    IF (spec_idt(id_so4, jt) == idt_ms4cs) idx_ms4cs = jt
  END DO


  END SUBROUTINE ham_m7_init_emissions

  !! ---------------------------------------------------------------------------
  !! subroutine to define additional diagnostics for the emi diagnostic stream
  SUBROUTINE ham_m7_init_emi_stream(nsectors, emi_stream, ldiagdetail)

  ! add field to emis stream (init_emi_stream in mo_emi_interface must have been called!)

  USE mo_linked_list,   ONLY: SURFACE
  USE mo_memory_base,   ONLY: t_stream, add_stream_element, default_stream_setting, &
                              AUTO, SURFACE
  USE mo_ham_dust,      ONLY: bgc_dust_init_diag
  USE mo_ham,           ONLY: idsec_dust

  INTEGER, INTENT(in)      :: nsectors
  TYPE(t_stream), POINTER  :: emi_stream
  LOGICAL, INTENT(inout)   :: ldiagdetail(nsectors)

  ! -- Diagnostics for dust emissions
  IF (idsec_dust > 0) THEN

    !-- Dust emission details (accumulated)
    CALL default_stream_setting (emi_stream, lrerun    = .FALSE.,     &
                                       laccu     = .TRUE. ,     &
                                       lpost     = .TRUE. ,     &
                                       leveltype = SURFACE,     &
                                       table     = 199,         &
                                       code      = AUTO         )

    CALL add_stream_element (emi_stream, 'emi_DU_mai', d_emi_mduai,     &
                             longname='emission mass flux of dust in insoluble accumulation mode',   &
                             units='kg m-2 s-1'  )
    CALL add_stream_element (emi_stream, 'emi_DU_nai', d_emi_nduai,     &
                             longname='emission number flux of dust in insoluble accumulation mode', &
                             units='1 m-2 s-1'  )
    CALL add_stream_element (emi_stream, 'emi_DU_mci', d_emi_mduci,     &
                             longname='emission mass flux of dust in insoluble coarse mode',         &
                             units='kg m-2 s-1'  )
    CALL add_stream_element (emi_stream, 'emi_DU_nci', d_emi_nduci,     &
                             longname='emission number flux of dust in insoluble coarse mode',       &
                             units='1 m-2 s-1'  )
    !! add ustar_acrit field (and other dust diagnostics?)
    CALL bgc_dust_init_diag(emi_stream)
    !! return info that diag for dust is taken care off
!sschr????
!   ldiagdetail(idsec_dust) = .FALSE.
!sschr????
  END IF

  ! -- Diagnostics for seasalt emissions
  IF (idsec_seasalt > 0) THEN
    !-- Seasalt emission details (accumulated)
    CALL default_stream_setting (emi_stream, lrerun    = .FALSE.,     &
                                 laccu     = .TRUE. ,     &
                                 lpost     = .TRUE. ,     &
                                 leveltype = SURFACE,     &
                                 table     = 199,         &
                                 code      = AUTO         )

    CALL add_stream_element (emi_stream, 'emi_SS_mas', d_emi_mssas,     &
                             longname='emission mass flux of seasalt in soluble accumulation mode',   &
                             units='kg m-2 s-1'  )
    CALL add_stream_element (emi_stream, 'emi_SS_nas', d_emi_nssas,     &
                             longname='emission number flux of seasalt in soluble accumulation mode', &
                             units='1 m-2 s-1'  )
    CALL add_stream_element (emi_stream, 'emi_SS_mcs', d_emi_msscs,     &
                             longname='emission mass flux of seasalt in soluble coarse mode',         &
                             units='kg m-2 s-1'  )
    CALL add_stream_element (emi_stream, 'emi_SS_ncs', d_emi_nsscs,     &
                             longname='emission number flux of seasalt in soluble coarse mode',       &
                             units='1 m-2 s-1'  )
    !! return info that diag for dust is taken care off
    ldiagdetail(idsec_seasalt) = .FALSE.
  END IF

  END SUBROUTINE ham_m7_init_emi_stream

  !! ---------------------------------------------------------------------------
  !! subroutine to handle special emissions for HAM
  SUBROUTINE ham_m7_emissions(kproma, kbdim, klev, krow, ktrac, ksec, kspec,     &
                              ibc_extra, pfactor)


  USE mo_time_control,       ONLY: delta_time
  USE mo_submodel_diag,      ONLY: t_diag_list, get_diag_pointer
  USE mo_submodel_streams,   ONLY: emi_lpost
  USE mo_species,            ONLY: nmaxtrspec
  USE mo_ham,                ONLY: aerocomp, nsoa, sizeclass, nseasalt, idsec_dust, idsec_biogenic
  ! for dust and seasalt emissions
  USE mo_ham_species,        ONLY: id_du, id_ss
  USE mo_ham_m7ctl,          ONLY: iduai, iduci, iacci, icoai, &
                                   issas, isscs, iaccs, icoas
  ! for carbon emissions
  USE mo_ham_species,        ONLY: id_bc, id_oc
  ! for sulfur emissions
  USE mo_ham_species,        ONLY: id_so2, id_so4
  !>>dod
  USE mo_ham_soa,            ONLY: soaprop 
  !<<dod
  !>>dod (redmine #44) import seasalt emission schemes from HAM2
  USE mo_ham_m7_emi_seasalt, ONLY:seasalt_emissions_monahan, seasalt_emissions_lsce,    &
                                  seasalt_emissions_mh,      seasalt_emissions_guelle,  &
                                  seasalt_emissions_gong, seasalt_emissions_long, seasalt_emissions_gong_SST
  !<<dod
  USE mo_boundary_condition,       ONLY: bc_set, bc_query
  USE mo_external_field_processor, ONLY: EF_MODULE

  USE mo_ham,   ONLY: ibc_dust, ibc_seasalt !alaak +

  ! arguments

  INTEGER,  INTENT(in)    :: kproma                   ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                    ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: klev                     ! number of levels
  INTEGER,  INTENT(in)    :: krow                     ! geographic block number
  INTEGER,  INTENT(in)    :: ktrac                    ! number of tracers
  INTEGER,  INTENT(in)    :: ksec                     ! emission sector number (match with idsec_xyz)
  INTEGER,  INTENT(in)    :: kspec                    ! species number
  INTEGER,  INTENT(inout) :: ibc_extra(nmaxtrspec)    ! some sectors (dust, seasalt) use more than one bc
  REAL(dp), INTENT(inout) :: pfactor(nmaxtrspec)      ! weighting factor for size distribution

  ! local variables
  REAL(dp)          :: zmassf1(kbdim), zmassf2(kbdim), znumf1(kbdim), znumf2(kbdim)
  REAL(dp), POINTER :: fld2d(:,:)
  INTEGER             :: ierr, ieftype

  !-- Initialisation

  ibc_extra(:) = 0

  !-- Test for sector-species combinations where HAM should take action

  !--- Dust emissions
  IF (ksec == idsec_dust .AND. kspec == id_du) THEN
! note: id_du = 0 should not occur here. Let it crash if it does ...
    CALL bc_query(ibc_dust, ef_type=ieftype)
    ! calculate dust emissions and map onto HAM modes
    IF (ieftype == EF_MODULE) THEN
      CALL ham_m7_dust_emissions(kproma, kbdim, krow)
      ibc_extra(idx_mduai) = ibc_dust
      ibc_extra(idx_mduci) = ibc_dust + 1
      ibc_extra(idx_mduas) = 0
      ibc_extra(idx_mducs) = 0
      ibc_extra(idx_nduai) = ibc_dust
      ibc_extra(idx_nduci) = ibc_dust + 1
      ibc_extra(idx_nduas) = 0
      ibc_extra(idx_nducs) = 0
      ! next statement for efficiency reasons
      ibc_extra(1) = 0
      pfactor(:) = 0._dp                  ! reset weighting factors
      pfactor(idx_mduai) = 1.0_dp
      pfactor(idx_mduci) = 1.0_dp
    ELSE
      pfactor(:) = 0._dp                  ! reset weighting factors
      pfactor(idx_mduai) = 0.1_dp
      pfactor(idx_mduci) = 0.9_dp
    ENDIF
    pfactor(idx_mduas) = 0.0_dp
    pfactor(idx_mducs) = 0.0_dp
    pfactor(idx_nduai) = zm2n_duai
    pfactor(idx_nduci) = zm2n_duci
    pfactor(idx_nduas) = 0.0_dp
    pfactor(idx_nducs) = 0.0_dp
  END IF

  !--- Seasalt emissions
  IF (ksec == idsec_seasalt .AND. kspec == id_ss) THEN
! note: id_ss = 0 should not occur here. Let it crash if it does ...
    ! identify seasalt tracers
    ! calculate seasalt emissions for HAM modes
    SELECT CASE(nseasalt)
    CASE (1)
       CALL seasalt_emissions_monahan(kproma, kbdim, krow, zmassf1, zmassf2, znumf1, znumf2)
    CASE (2)
       CALL seasalt_emissions_lsce(kproma, kbdim, krow, zmassf1, zmassf2, znumf1, znumf2)
    !>>dod (redmine #44) import of seasalt schemes from HAM2   
    !   CASE(3)      
       !...not yet implemented
    CASE(4)
       CALL seasalt_emissions_mh(kproma, kbdim, krow, zmassf1, zmassf2, znumf1, znumf2)
    CASE(5)
       CALL seasalt_emissions_guelle(kproma, kbdim, krow, zmassf1, zmassf2, znumf1, znumf2)
    CASE(6)
       CALL seasalt_emissions_gong(kproma, kbdim, krow, zmassf1, zmassf2, znumf1, znumf2)
    CASE(7)
       CALL seasalt_emissions_long(kproma, kbdim, krow, zmassf1, zmassf2,znumf1, znumf2)
    CASE(8)
       CALL seasalt_emissions_gong_SST(kproma, kbdim, krow, zmassf1, zmassf2,znumf1,znumf2)


    END SELECT
    !<<dod

    ibc_extra(idx_mssas) = ibc_seasalt
    ibc_extra(idx_msscs) = ibc_seasalt + 1
    ibc_extra(idx_nssas) = ibc_seasalt + 2
    ibc_extra(idx_nsscs) = ibc_seasalt + 3
    ! next statement for efficiency reasons
    ibc_extra(1) = 0
    CALL bc_set(ibc_seasalt, kproma, krow, zmassf1)
    CALL bc_set(ibc_seasalt+1, kproma, krow, zmassf2)
    CALL bc_set(ibc_seasalt+2, kproma, krow, znumf1)
    CALL bc_set(ibc_seasalt+3, kproma, krow, znumf2)
    pfactor(:) = 0._dp                  ! reset weighting factors
    pfactor(idx_mssas) = 1.0_dp
    pfactor(idx_msscs) = 1.0_dp
    pfactor(idx_nssas) = 1.0_dp
    pfactor(idx_nsscs) = 1.0_dp
  END IF

  !--- Black carbon emissions
  IF (kspec == id_bc) THEN
    pfactor(:) = 0._dp                  ! reset weighting factors
    pfactor(idx_mbcki) = 1._dp          ! all black carbon mass goes into insoluble aitken mode
    ! apply conversion factor mass to number
    IF (ksec == idsec_fire .OR. ksec == idsec_ffire .OR. ksec == idsec_gfire    &
!gf(#140)        .OR. ksec == idsec_awb .OR. ksec == idsec_biofuel) THEN
        .OR. ksec == idsec_awb .OR.ksec == idsec_dom .OR. ksec == idsec_biofuel) THEN
      pfactor(idx_nbcki) = zm2n_bcki_bb   ! use biomass burning factor
    ELSE
      pfactor(idx_nbcki) = zm2n_bcki_ff   ! use fossil fuel factor
    END IF
  END IF

  !--- Organic carbon emissions
  IF (kspec == id_oc) THEN
    ! default mode distribution is for fossil fuel emissions
    pfactor(:) = 0._dp                  ! reset weighting factors
    pfactor(idx_mocki) = zom2oc         ! all organic carbon mass goes into insoluble aitken mode
    pfactor(idx_nocki) = zom2oc*zm2n_ocki_ff
    !--- biomass burning, agricultural waste burning and domestic
    IF (ksec == idsec_fire .OR. ksec == idsec_ffire .OR. ksec == idsec_gfire     &
!gf(#140)        .OR. ksec == idsec_awb .OR. ksec ==  idsec_biofuel) THEN
        .OR. ksec == idsec_awb .OR. ksec == idsec_dom .OR. ksec ==  idsec_biofuel) THEN
      pfactor(idx_mocki) = zom2oc*(1._dp-zbb_wsoc_perc)     ! insoluble fraction
      pfactor(idx_mocks) = zom2oc*zbb_wsoc_perc             ! soluble fraction
      pfactor(idx_nocki) = zom2oc*(1._dp-zbb_wsoc_perc)*zm2n_ocki_bb  ! number for insoluble aitken
      pfactor(idx_nocks) = zom2oc*zbb_wsoc_perc*zm2n_ocks_bb          ! number for soluble aitken
    END IF
    !--- biogenic
    !--- if SOA scheme from D. O'Donnell is active, the biogenic sector does
    !    not emit aerosols (only precursors)
    IF (ksec == idsec_biogenic .AND. nsoa /= 1) THEN
      !! no OM to OC conversion here ###
      ! 50% goes into aitken mode and 50% into accumulation mode
      ! note that soluble number concentration is not altered because it is assumed that
      ! soluble biogenic particles immediately condense on other particles (P. Stier, 2010)
      pfactor(idx_mocki) = (1._dp-zbg_wsoc_perc)               ! insoluble fraction
      pfactor(idx_mocks) = 0.5*zbg_wsoc_perc                   ! soluble fraction (aitken)
      pfactor(idx_mocas) = 0.5*zbg_wsoc_perc                   ! soluble fraction (accumulation)
      pfactor(idx_nocki) = (1._dp-zbg_wsoc_perc)*zm2n_ocki_bg  ! number for insoluble aitken
    END IF
  END IF

  !--- Sulfur emissions
! ### NOTE: lmoz=.true. will require the sulfur tracer in units of kg(SO2) rather than kg(S) !! ###
  ! Note: emi_matrix only contains SO2 as species, but a fraction (1-zso2frac) of SO2 
  ! emissions is assumed to go directly into the (SO4) aerosol phase
  IF (kspec == id_so2) THEN
    pfactor(:) = 0._dp
    pfactor(1) = zfacso2
  END IF
  IF (kspec == id_so4) THEN
    ! default mode distribution is for surface level fossil fuel emissions
    ! and also valid for volcanic emissions
    ! 50% aitken mode and 50% accumulation mode
    pfactor(:) = 0._dp                  ! reset weighting factors
    pfactor(idx_ms4ks) = 0.5_dp * zso2tso4 * (1._dp-zfacso2)
    pfactor(idx_ms4as) = 0.5_dp * zso2tso4 * (1._dp-zfacso2)
    pfactor(idx_ns4ks) = 0.5_dp * zso2tso4 * (1._dp-zfacso2) * zm2n_s4ks_ff
    pfactor(idx_ns4as) = 0.5_dp * zso2tso4 * (1._dp-zfacso2) * zm2n_s4as_ff
    !--- energy sector and industrial emissions (smoke stacks)
    ! 50% accumulation mode and 50% coarse mode
!gf(#140) - Add ships and remove industrial sector
!    IF (ksec == idsec_ene .OR. ksec == idsec_ind) THEN
!gf(#140)    IF (ksec == idsec_ene .OR. ksec == idsec_ind .OR. ksec == idsec_ships) THEN
    IF (ksec == idsec_ene .OR. ksec == idsec_ships) THEN
!gf
      pfactor(idx_ms4ks) = 0._dp
      pfactor(idx_ms4cs) = 0.5_dp * zso2tso4 * (1._dp-zfacso2)
      pfactor(idx_ns4ks) = 0._dp
      pfactor(idx_ns4cs) = 0.5_dp * zso2tso4 * (1._dp-zfacso2) * zm2n_s4cs_ff
    END IF
    !--- biomass burning, agricultural waste burning and domestic
    ! 50% aitken mode and 50% accumulation mode
!gf(#140)    IF (ksec == idsec_fire .OR. ksec == idsec_ffire .OR. ksec == idsec_gfire .OR. ksec == idsec_awb) THEN
    IF (ksec == idsec_fire .OR. ksec == idsec_ffire .OR. ksec == idsec_gfire .OR. ksec == idsec_awb .OR. ksec == idsec_dom) THEN
      pfactor(idx_ns4ks) = 0.5_dp * zso2tso4 * (1._dp-zfacso2) * zm2n_s4ks_bb
      pfactor(idx_ns4as) = 0.5_dp * zso2tso4 * (1._dp-zfacso2) * zm2n_s4as_bb
    END IF
  END IF
 
  END SUBROUTINE ham_m7_emissions

  !! ===========================================================================
  !! private subroutines to interface with external modules (dust, seasalt)

  !! ---------------------------------------------------------------------------
  SUBROUTINE ham_m7_dust_emissions(kproma, kbdim, krow)

  ! Description:
  ! ------------
  ! *ham_m7_dust_emissions* is the interface between the MPI BGC dust emission scheme
  ! and the HAM/M7 aerosol module. Dust number flux is calculated from the mass flux and
  ! both fluxes are stored as ECHAM5 tracers.
  !
  ! Author:
  ! -------
  ! M. Werner, MPI BGC, November 2002
  !
  ! Interface:
  ! ----------
  ! *ham_m7_dust_emissions* is called from *ham_m7_emissions* which is called from
  ! *emi_interface*.
  !
  ! Adaptation to ECHAM6-HAMMOZ code structure by Martin Schultz, FZ Juelich (2010-02-11)

  USE mo_memory_g3b,         ONLY: slm, glac
  USE mo_ham_species,        ONLY: id_du
  USE mo_exception,          ONLY: finish
  USE mo_ham_dust,           ONLY: bgc_dust_calc_emis, flux_6h
  USE mo_ham_m7ctl,          ONLY: sigma, iacci, icoai
  USE mo_boundary_condition, ONLY: bc_set
  USE mo_ham,                ONLY: ibc_dust

  IMPLICIT NONE

  INTEGER,  INTENT(in)    :: kproma                   ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                    ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: krow                     ! geographic block number

  !--- Local Variables:
  REAL(dp)            :: pmassf_ai(kbdim)             ! mass flux of dust in accumulation mode
  REAL(dp)            :: pmassf_ci(kbdim)             ! ..                in coarse mode
  REAL(dp)            :: ztmp1(kbdim)                 !SF #458 temporary variable
  LOGICAL             :: ll1(kbdim)                   !SF #458 temporary variable
  INTEGER,  PARAMETER :: min_ai = 1                   ! index boundaries for BGC_Dust bin numbers
  INTEGER,  PARAMETER :: max_ai = 1
  INTEGER,  PARAMETER :: min_ci = 2
  INTEGER,  PARAMETER :: max_ci = 4

  INTEGER             :: nn

  !--- Consistency check:
  IF( ABS(ref_sigma_dust_ai-sigma(iacci))>1.E-3_dp .OR. ABS(ref_sigma_dust_ci-sigma(icoai))>1.E-3_dp) THEN
    CALL finish('ham_m7_dust_emis','inconsistent source standard deviation')
  END IF

  !--- 1) Calculate BGC dust emissions: result is in g m-2 s-1
  CALL bgc_dust_calc_emis(kproma, krow)

  ! Sum over BGC_Dust internal tracer classes for HAM modal scheme:
  pmassf_ai(:) = 0._dp
  pmassf_ci(:) = 0._dp

  ! Accumulation mode: trc 1
  DO nn=min_ai, max_ai
    pmassf_ai(1:kproma) = pmassf_ai(1:kproma) + flux_6h(1:kproma,nn,krow)
  END DO

  ! Coarse mode: trc 2,3,4
  DO nn=min_ci, max_ci
    pmassf_ci(1:kproma) = pmassf_ci(1:kproma) + flux_6h(1:kproma,nn,krow)
  ENDDO 

  ! Make sure mass flux is positive !SF #458 (replacing where statements)
  ll1(1:kproma) = (pmassf_ai(1:kproma) <= 0._dp)

  pmassf_ai(1:kproma) = MERGE( &
                              0._dp, &
                              pmassf_ai(1:kproma), &
                              ll1(1:kproma))

  pmassf_ci(1:kproma) = MERGE( &
                              0._dp, &
                              pmassf_ci(1:kproma), &
                              ll1(1:kproma))

  ! Mask out glacier and ocean values on ECHAM grid !SF #458 (replacing where statements)
  ll1(1:kproma) = (glac(1:kproma,krow) > 0.5_dp) .OR. (slm(1:kproma,krow) < 0.5_dp)

  pmassf_ai(1:kproma) = MERGE( &
                              0._dp, &
                              pmassf_ai(1:kproma), &
                              ll1(1:kproma))

  pmassf_ci(1:kproma) = MERGE( &
                              0._dp, &
                              pmassf_ci(1:kproma), &
                              ll1(1:kproma))

  ! Reduce flux over glacier grids according to glacier fraction !SF #458 (replacing where statements)
  ll1(1:kproma)   = (glac(1:kproma,krow) > 0.1_dp)
  ztmp1(1:kproma) = (1._dp-glac(1:kproma,krow))**3

  pmassf_ai(1:kproma) = MERGE( &
                              pmassf_ai(1:kproma) * ztmp1(1:kproma), &
                              pmassf_ai(1:kproma), &
                              ll1(1:kproma))

  pmassf_ci(1:kproma) = MERGE( &
                              pmassf_ci(1:kproma) * ztmp1(1:kproma), &
                              pmassf_ci(1:kproma), &
                              ll1(1:kproma))
  
  ! Reduce flux over not-all-land grids according to land fraction !SF #458 (replacing where statements)
  ll1(1:kproma)   = (slm(1:kproma,krow) < 0.99_dp)
  ztmp1(1:kproma) = MIN(slm(1:kproma,krow),0.7_dp)**4

  pmassf_ai(1:kproma) = MERGE( &
                              pmassf_ai(1:kproma) * ztmp1(1:kproma), &
                              pmassf_ai(1:kproma), &
                              ll1(1:kproma))

  pmassf_ci(1:kproma) = MERGE( &
                              pmassf_ci(1:kproma) * ztmp1(1:kproma), &
                              pmassf_ci(1:kproma), &
                              ll1(1:kproma))
  
  ! Convert to kg m-2 s-1
  pmassf_ai(1:kproma) = MAX(pmassf_ai(1:kproma), 0._dp) / 1000._dp 
  pmassf_ci(1:kproma) = MAX(pmassf_ci(1:kproma), 0._dp) / 1000._dp 

  CALL bc_set(ibc_dust, kproma, krow, pmassf_ai)
  CALL bc_set(ibc_dust+1, kproma, krow, pmassf_ci)

  END SUBROUTINE ham_m7_dust_emissions

  !>>dod (redmine #44) import of seasalt emission schemes from HAM2:
  !      moved seasalt emission routines to own module (mo_ham_emi_seasalt)
END MODULE mo_ham_m7_emissions
