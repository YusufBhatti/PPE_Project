!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! hammoz_perturbations.inc
!!
!! \brief
!! Namelist for aerosol parameter perturbations
!!
!! \author Philip Stier (Oxford)
!!
!! \responsible_coder
!! Philip Stier, philip.stier@physics.ox.ac.uk
!!
!! \revision_history
!!   -# Philip Stier (Oxford) - original code (2013)
!!   -# Duncan Watson-Parris (Oxford) - Added wet depostion, fire height, BC ni and
!!                                      convection perturbations
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

NAMELIST /hammoz_perturbations/ &
!--- Sources:
  scale_nuc_bl,             & ! Scale factor for BL nucleation scheme (activation or kinetic)
  scale_nuc_ft,             & ! Scale factor for free tropospheric nucleation scheme (Vehkamaeki or Kazil & Lovejoy)
  scale_emi_cmr_ff,         & ! Scale factor for emission count median raduius for fossil fules
  scale_emi_cmr_bb,         & ! Scale factor for emission count median raduius for wildfires
  scale_emi_cmr_bf,         & ! Scale factor for emission count median raduius for biogenic emissions
  scale_emi_ff,             & ! Scale factor for fossil fule emissions
  scale_emi_bb,             & ! Scale factor for wildfire emissions
  scale_emi_bf,             & ! Scale factor for biofuel emissions
  scale_emi_bc,             & ! Scale factor for BC emissions (all sectors)
  scale_emi_dms,            & ! Scale factor for dms emissions (terrestrial)
  scale_emi_ssa,            & ! Scale factor for SSA emissions (all sectors) 
  scale_emi_du,             & ! Scale factor for DUST emissions (all sectors)
  scale_emi_so2,             & ! Scale factor for ANTH SO2 emissions (all sectors)
!  emi_prim_so4_frac,        & ! Absolute fraction of SO2 emitted as primary SO4
!  emi_prim_so4_cmr,         & ! Absolute count median radius of primary SO4 particles [m]
!  prop_fire_in_pbl_p1,      & ! The proportion of (G/F)FIRE tracer emitted into the layer imediately above the PBL (PBL+1)
!  prop_fire_in_pbl_p2,      & ! The proportion of (G/F)FIRE tracer emitted into PBL+2
!  prop_fire_in_pbl,         & ! The proportion of (G/F)FIRE tracer emitted into the PBL
!--- Sinks:
  scale_drydep_ait,         & !< Scale factors for dry deposition of Aitken modes
  scale_drydep_acc,         & !< Scale factors for dry deposition of accumulation modes
  scale_drydep_coa,         & !< Scale factors for dry deposition of coarse modes i
  scale_wetdep_ic,          & ! Scale factor for in-cloud wet deposition
  scale_wetdep_bc,          & ! Scale factor for below-cloud wet deposition
!  scale_wetdep_ic_bc_only,  & ! Scale factor for in-cloud wet deposition of BC (independantly of the above scalings)
!  scale_wetdep_bc_bc_only,  & ! Scale factor for below-cloud wet deposition of BC (independantly of the above scalings)
!  scale_sed
!--- Ageing:
  scale_so4_coating,        & ! Scale the coating thickness of SO4 required to 'age' particles
!  scale_intra_mode_coagulation,   & ! Scale the diagonal elements of the coagulation kernel matrix
!  scale_inter_mode_coagulation,   & ! Scale the off-diagonal elements of the coagulation kernel matrix
!--- Convection:
!  scale_tr_entrainment,     & ! Scale aerosol entrainment rate
!  scale_vertical_velocity,  & ! Scale (total) vertical velocity - ONLY for nactivpdf == 0
!--- Aerosol processes:
!---
!--- Radiative properties:
  bc_rad_ni,          & ! Scale the imaginary part of the refractive index (for BC)
!  oc_rad_ni,          & ! Scale the imaginary part of the refractive index (for OC)
  du_rad_ni,          & ! Scale the imaginary part of the refractive index (for DUST)
!--- Cloud scheme (ONLY for naut == 2):
!  KK_exponent,              & ! The exponent to use in the KK autoconversion scheme
!  KK_LWP_exponent,            & ! The LWP exponent in the KK autoconversion scheme.
!  scale_activation,         & ! Scale the number of ARG activated particles
!  scale_accretion,          & ! Scale the accretion
!  scale_solar_const,        & ! scale the solar constant (equally across all bands)
  scale_so2_reactions, &
!  scale_water, &
!  scale_RH,    &
  pH_pert,     & ! Scale the cloud-water pH
  scale_kappa_ss, &
  scale_kappa_so4
