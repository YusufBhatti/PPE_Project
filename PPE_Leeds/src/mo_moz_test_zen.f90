!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_moz_test_zen.f90
!!
!! \brief
!! this module tests the zenith angle calculation from MOZART vs. the ECHAM solar_parameters code.
!!
!! \author Martin Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (XXXX-XX-XX)
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

#if defined (NAG)
#define ARGCHECK 1
#endif

module mo_moz_test_zen

      use mo_kind,   only: dp, wp
      
      implicit none

      private
      public :: test_zen

      save

      real(dp) :: esfact                   ! relative earth sun distance factor (1/R^2)

      contains

  subroutine test_zen( kproma, kbdim, krow )


  USE mo_radiation_parameters, ONLY: lyr_perp, yr_perp, nmonth, isolrad, nb_sw,   &
                                     twopi, deg2rad, cemiss, diff, solc,          &
                                     solar_parameters

  USE mo_time_conversion, ONLY: TC_convert, TC_get, time_intern
  USE mo_time_control,    ONLY: current_date, ndaylen, l_orbvsop87, get_orbit_times
  USE mo_geoloc,          ONLY: amu0_x, rdayl_x, amu0m_x, rdaylm_x, coslon_2d, &
       &                        sinlon_2d, sinlat_2d, coslat_2d
  USE mo_orbit,           ONLY: orbit_vsop87, orbit_kepler
  USE mo_moz_waccm_photo, ONLY: diurnal_geom
  USE mo_exception,       ONLY : message, message_text, em_debug

!-----------------------------------------------------------------------
!     	... Local variables
!-----------------------------------------------------------------------
  integer, intent(in) :: kproma, kbdim, krow
  integer, parameter :: ip = 1    ! longitude tile index
  ! Timer
  TYPE(time_intern)  :: mydate
  INTEGER            :: yymmdd, hhmmss, sec, day
  ! Variables for diurnal_geometry. Not used except for zen_angle in waccm_photo.
  real(dp)           :: caldayn
  real(dp)           :: zen_angle(kbdim), loc_angle(kbdim)
  real(dp)           :: sunon(kbdim), sunoff(kbdim)
  logical            :: polar_night(kbdim), polar_day(kbdim)

  REAL(wp)           :: rasc_sun, decl_sun, dist_sun, orbit_date, time_of_day
  REAL(wp)           :: flx_ratio_cur
  logical            :: l_rad_call

  !-- Initialize
  CALL TC_convert(current_date, mydate)
  CALL TC_get(mydate, yymmdd, hhmmss)
  CALL TC_get(current_date, day, sec)


!--- Date and time
  caldayn = caldayr( yymmdd, sec )
  WRITE(message_text,'(a,e25.15,2i0)') 'caldayn, yymmdd, sec = ',caldayn, yymmdd, sec
  CALL message('test_zen', message_text, level=em_debug)

!--- Calculate parameters for diurnal geometry   (MOZART)
!        The only parameter which is really used lateron is zen_angle
!        in waccm_photo.
!        We keep the other parameters from diurnal_geometry to preserve
!        the calling sequence to waccm_photo. Note that they are no longer
!        calculated in the diurnal_geometry subroutine.

  call diurnal_geom( ip, krow, caldayn, polar_night, polar_day, &
                     sunon, sunoff, loc_angle, zen_angle )

!--- Calculate zenith angle (and other things) (ECHAM)
    !
    ! 1.0 Compute orbital parameters for current time step
    ! --------------------------------
    l_rad_call = .FALSE.
    CALL get_orbit_times(l_rad_call, lyr_perp, nmonth, yr_perp, time_of_day, &
         &               orbit_date)

    IF (l_orbvsop87) THEN
      CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
    ELSE
      CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
    END IF
    CALL solar_parameters(decl_sun, dist_sun, time_of_day, &
         &                sinlon_2d, sinlat_2d, coslon_2d, coslat_2d, &
         &                flx_ratio_cur, amu0_x, rdayl_x)

  end subroutine test_zen


      integer function DOY( mon, cday )
!-----------------------------------------------------------------------
!   ... Compute day of year ignoring leap years.
!           Returns values in the range [1,365].
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: mon, cday

      integer, save :: jdcon(12) = &
            (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)

      doy = jdcon(mon) + cday

      end function DOY


      real(dp) function CALDAYR( idate, isec )
!-----------------------------------------------------------------------
!   ... Calendar day with fractional part.  Returns values in the
!           range [1., 366.)
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: idate, isec

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
      integer :: mon, day

      mon = MOD( idate,10000 ) / 100
      day = MOD( idate,100 )

      CALDAYR = DOY( mon, day ) + REAL( isec, dp )/86400._dp

      end function CALDAYR

end module mo_moz_test_zen
