!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_atmosphere

  USE mo_jsbach_constants, ONLY: Gravity
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

CONTAINS

  ELEMENTAL FUNCTION height_above_surface(geopotential, elevation)
    !
    ! Converts geopotential to height above surface
    !
    REAL(dp), INTENT(in)  :: geopotential   ! Geopotential [m^2/s^2]
    REAL(dp), INTENT(in)  :: elevation      ! Elevation above sea level of surface
    REAL(dp)              :: height_above_surface

    height_above_surface = geopotential/Gravity - elevation

  END FUNCTION height_above_surface
  !
  !----------------------------------------------------------------------------------------------------------------------
  FUNCTION sat_specific_humidity(temp, pressure) RESULT(qsat)
    !
    ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
    ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
    ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)
    !
    USE mo_echam_convect_tables, ONLY: tlucua, &   ! Table for water vapor pressure multiplied by R_d/R_v
                                 jptlucu1, jptlucu2

    REAL(dp), INTENT(in) :: temp(:)      ! Air temperature at level [K]
    REAL(dp), INTENT(in) :: pressure(:)  ! Pressure at level [Pa]
    REAL(dp)             :: qsat(SIZE(temp,1))

    REAL(dp) :: vtmpc1
    INTEGER :: it(SIZE(temp,1))
    REAL(dp) :: tluc(SIZE(temp,1))

    vtmpc1 = 0.6077686814143877_dp
    it = NINT(temp*1000._dp)
    tluc = tlucua(it)

    WHERE (it >= jptlucu1 .AND. it <= jptlucu2)
       qsat = tluc / (pressure - vtmpc1*tluc)
    ELSEWHERE
       qsat = HUGE(1._dp)
    END WHERE

  END FUNCTION sat_specific_humidity

END MODULE mo_atmosphere
