!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! set_flighttrack_element.inc
!!
!! \brief
!! Namelist for the flighttrack submodel (related to the flighttrack stream)
!!
!! \author Zak Kipling (Uni Oxford)
!!
!! \responsible_coder
!! Zak Kipling, zak.kipling@physics.ox.ac.uk
!!
!! \revision_history
!!   -# Zak Kipling (Uni Oxford) - ported from ECHAM5.5-HAM2 (2013)
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

NAMELIST /set_flighttrack_element/ &

  stream,   &      ! stream name
  element          ! element name
