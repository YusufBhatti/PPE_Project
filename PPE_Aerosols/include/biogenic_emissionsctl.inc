!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! biogenic_emissionsctl.inc
!!
!! \brief
!! General namelist for biogenic emissions
!!
!! \author Grazia Frontoso (C2SM)
!!
!! \responsible_coder
!! Grazia Frontoso, grazia.frontoso@env.ethz.ch
!!
!! \revision_history
!!   -# none
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

NAMELIST /biogenic_emissionsctl/ &
         nlai_biogenic_ef_type, & !< Choice of lai external field type in the
                                  !< biogenic emission module (MEGAN)
         nef_pft, &               !< Choice of PFT fractions from MEGAN-CLM4 or JSBACH
         emifact_files_species, & !< Choice of emission factors
                                  !< read from specific file or calculated with PFT fractions
         ldebug_bioemi            !< Switch on detailed output for bioemi diagnostics
