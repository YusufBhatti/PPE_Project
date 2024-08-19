!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_data.f90
!!
!! \brief
!! Module for storing AeroCom utility variables (mostly constants?)
!!
!! \author D. Neubauer (ETH Zurich)
!!  adapted from AerChemMIP diagnostics from S. Ferrachat
!!
!! D. Neubauer, david.neubauer@env.ethz.ch
!!
!! \revision_history
!!   -# D. Neubauer (ETH Zurich) - original code (2018-08-14)
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
MODULE mo_hammoz_aerocom_data

  USE mo_filename, ONLY: NETCDF4, NETCDF, GRIB

  IMPLICIT NONE

  !debug INTEGER, PARAMETER :: AEROCOM_FILETYPE = NETCDF4
  INTEGER, PARAMETER :: AEROCOM_FILETYPE = NETCDF
  INTEGER, PARAMETER :: AEROCOM_FILETYPE_GRIB = GRIB

END MODULE mo_hammoz_aerocom_data
