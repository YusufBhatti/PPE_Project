!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_zenith
  ! 
  ! Description: 
  !   <Say what this module is for> 
  ! 
  ! Current Code Owner: <Name of person responsible for this code> 
  ! 
  ! History: 
  !  
  ! Version   Date     Comment 
  ! -------   ----     ------- 
  ! <version> <date>   Original code. <Your name> 
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 
  ! Modules used: 
  !
  USE mo_kind,         ONLY: dp
  USE mo_jsbach_grid , ONLY: domain_type

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: compute_orbit_and_solar, init_zenith, cos_zenith

  REAL(dp), ALLOCATABLE, SAVE :: cos_zenith(:)    !< Cosine of solar zenith angle for domain

  LOGICAL :: module_initialized = .FALSE.

CONTAINS 

  !=================================================================================================
  SUBROUTINE compute_orbit_and_solar(domain)

    ! Compute orbital parameters for current time step

    ! Called from *jsbalone_driver* at the beginning of each time step
    ! If coupled with ECHAM: 
    !    this is done in *pre_radiation* from *mo_radiation*, called in *scan1*
    !    the cosine of zenith angle is passed from ECHAM to the interface
    !    the interface copies the packed zenith angle to *cos_zenith*

    USE mo_time_control, ONLY: l_orbvsop87, get_orbit_times
    USE mo_orbit,        ONLY: orbit_kepler, orbit_vsop87
    USE mo_radiation_parameters,  ONLY: nmonth, solar_parameters, yr_perp, lyr_perp

    TYPE(domain_type), INTENT(in) :: domain

    LOGICAL :: l_rad_call = .FALSE.

    REAL(dp) :: rasc_sun, decl_sun, dist_sun, orbit_date, time_of_day
    REAL(dp) :: flx_ratio, cos_mu0(domain%nland,1), daylght_frc(domain%nland,1)
    REAL(dp), DIMENSION(domain%nland,1) :: sinlon, sinlat, coslon, coslat
    INTEGER :: nland

    nland = domain%nland

    IF (.NOT. module_initialized) CALL init_zenith(domain)

    CALL get_orbit_times(l_rad_call, lyr_perp, nmonth, yr_perp, time_of_day, orbit_date)

    IF (l_orbvsop87) THEN
       CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
    ELSE
       CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
    END IF

    sinlon(:,1) = domain%sinlon(1:nland)
    sinlat(:,1) = domain%sinlat(1:nland)
    coslon(:,1) = domain%coslon(1:nland)
    coslat(:,1) = domain%coslat(1:nland)

    CALL solar_parameters(decl_sun, dist_sun, time_of_day       &
                         ,sinlon, sinlat, coslon, coslat        &
                         ,flx_ratio, cos_mu0, daylght_frc)

    cos_zenith(1:nland) = cos_mu0(1:nland,1)

  END SUBROUTINE compute_orbit_and_solar
  !
  !=================================================================================================
  SUBROUTINE init_zenith(domain)

    ! Pre-compute some work quantities for domain
    
    TYPE(domain_type), INTENT(in) :: domain

    IF (module_initialized) RETURN

    ALLOCATE(cos_zenith(domain%nland))

    module_initialized = .TRUE.

  END SUBROUTINE init_zenith

END MODULE mo_zenith
