!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! ham_dustctl.inc
!!
!! \brief
!! namelist for the ECHAM/HAM BGC dust emissions
!!
!! \author Martin G. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# Martin G. Schultz (FZ Juelich) - original code (2009-10-01)
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

NAMELIST /ham_dustctl/ &

  ndurough,               & !< Global surface roughness length parameter.
  nduscale_reg,           & !< Regional threshold wind friction velocity parameter #csld #433
                            !! This is an array of shape/dim (8)
                            !! The indices correspond to the following regions:
                            !!     1       All the other locations than the following regions
                            !!     2       North america
                            !!     3       South America 
                            !!     4       North Africa 
                            !!     5       South Africa 
                            !!     6       Middle East 
                            !!     7       Asia
                            !!     8       Australia 
  r_dust_lai,             & !< Parameter for the threshold lai value
  r_dust_umin,            & !< Minimum U* for dust mobilization 
  r_dust_z0s,             & !< Surface roughness length for smooth surfaces 
  r_dust_scz0,            & !< Scale factor of surface roughness length
  r_dust_z0min,           & !< Minimum surface roughness length 
  k_dust_smst,            & !< Effect of soil moisture on threshold wind velocity of dust emission 
  k_dust_easo,            & !< New soil type for East Asia region 
  r_dust_sf13,            & !< Parameter of the duscale value over Takelimakan desert 
  r_dust_sf14,            & !< Parameter of the duscale value over Loess
  r_dust_sf15,            & !< Parameter of the duscale value over Gobi desert 
  r_dust_sf16,            & !< Parameter of the duscale value over other mixture soils
  r_dust_sf17,            & !< Parameter of the duscale value over desert and sand land 
  r_dust_af13,            & !< Parameter of the G/F value over Takelimakan desert 
  r_dust_af14,            & !< Parameter of the G/F value over Loess
  r_dust_af15,            & !< Parameter of the G/F value over Gobi desert 
  r_dust_af16,            & !< Parameter of the G/F value over other mixture soils
  r_dust_af17               !< Parameter of the G/F value over desert and sand land 
