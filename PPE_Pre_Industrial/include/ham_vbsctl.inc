!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! ham_salsactl
!!
!! \brief
!! namelist for the SALSA aerosol model
!!
!! \author Anton Laakso (FMI)
!!
!! \responsible_coder
!! Anton Laakso, anton.laakso@fmi.fi
!!
!! \revision_history
!!   -# Anton Laakso (FMI) +SALSA (2013-6)
!!   -# Thomas Kühn (UEF) + volatility basisi set (2015-06)
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

NAMELIST /ham_vbsctl/ nvbs_setup, &  !! Choice of setup
                 laqsoa              !! wet SOA scheme
                                         
