!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!++mgs: new module 14.03.2010
!++mgs: added decl_sun_cur (for MOZ photolysis) 02.06.2010
!! @brief Module to provide parameters to radiation routines and avoid circular dependencies.
!!
!! @remarks
!!   This module contains the public parameters provided by the radiation module
!!   mo_radiation.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19):
!!
!!         Martin Schultz, FZJ, Juelich (2010-04-13):
!!              extracted parameters from mo_radiation
!!         Dagmar Popke, MPI-M, Hamburg (2013-11-15):
!!              Implementation of RCE
!!
!! $ID: n/a$
!!
!! @par Origin
!!   see mo_radiation.f90
!!
!
MODULE mo_radiation_parameters

  USE mo_kind,            ONLY: wp
  USE mo_math_constants,  ONLY: pi
  USE mo_control,         ONLY: lrce

IMPLICIT NONE

  PRIVATE

  PUBLIC :: ldiur, lradforcing
  PUBLIC :: lyr_perp, yr_perp, nmonth, isolrad, nb_sw
  PUBLIC :: lw_spec_samp, sw_spec_samp, lw_gpts_ts, sw_gpts_ts, rad_perm
  PUBLIC :: i_overlap, l_do_sep_clear_sky
  PUBLIC :: ih2o, ico2, ich4, io3, io2, in2o, icfc, ighg, iaero, fco2
  PUBLIC :: co2vmr, ch4vmr, o2vmr, n2ovmr, cfcvmr
  PUBLIC :: co2mmr, ch4mmr, o2mmr, n2ommr            
  PUBLIC :: twopi, deg2rad, ch4_v, n2o_v
  PUBLIC :: cemiss, diff, rad_undef, solc
  PUBLIC :: psct, psctm, ssi_factor, flx_ratio_cur, flx_ratio_rad, decl_sun_cur 
  PUBLIC :: l_interp_rad_in_time, zepzen
  
  ! 1.0 NAMELIST global variables and parameters
  ! --------------------------------
  LOGICAL :: ldiur       = .TRUE.  !< diurnal cycle
  LOGICAL :: lradforcing(2) = (/.FALSE.,.FALSE./) !< &! switch on/off diagnostic 
             !of instantaneous aerosol solar (lradforcing(1)) and 
             !thermal (lradforcing(2)) radiation forcing
  LOGICAL :: lyr_perp    = .FALSE. !< switch to specify perpetual vsop87 year
  INTEGER :: yr_perp     = -99999  !< year if (lyr_perp == .TRUE.)
  INTEGER :: nmonth      =  0      !< 0=annual cycle; 1-12 for perpetual month
  ! nmonth currently works for zonal mean ozone and the orbit (year 1987) only
  INTEGER :: isolrad     =  3      !< mode of solar constant calculation
                                   !< default is rrtm solar constant
  INTEGER :: nb_sw      !< number of shortwave bands, set in setup
  
  ! Spectral sampling
  INTEGER :: lw_spec_samp = 1, sw_spec_samp = 1 ! 1 is broadband, 2 is MCSI, 3 and up are teams
  INTEGER :: lw_gpts_ts = 1, sw_gpts_ts = 1     ! Number of g-points per time step using MCSI
  INTEGER :: rad_perm = 0                       ! Integer for perturbing random number seeds
  ! Radiation driver
  LOGICAL :: l_do_sep_clear_sky = .TRUE. ! Compute clear-sky fluxes by removing clouds
  INTEGER :: i_overlap = 1               ! 1 = max-ran, 2 = max, 3 = ran
  LOGICAL :: l_interp_rad_in_time = .TRUE. 
  !
  ! --- Switches for radiative agents
  !
  INTEGER :: ih2o  = 1  !< water vapor, clouds and ice for radiation
  INTEGER :: ico2  = 2  !< carbon dioxide
  INTEGER :: ich4  = 3  !< methane
  INTEGER :: io3   = 3  !< ozone
  INTEGER :: io2   = 2  !< molecular oxygen
  INTEGER :: in2o  = 3  !< nitrous oxide
  INTEGER :: icfc  = 2  !< cfc11 and cfc12
  INTEGER :: ighg  = 0  !< greenhouse gase scenario
  INTEGER :: iaero = 2  !< aerosol model
  REAL(wp) :: fco2  = 1._wp !< factor for external co2 scenario (ico2=4)
  !
  ! --- Default gas volume mixing ratios - 1990 values (CMIP5)
  !
  REAL(wp) :: co2vmr    =  353.9e-06_wp !< CO2
  REAL(wp) :: ch4vmr    = 1693.6e-09_wp !< CH4
  REAL(wp) :: o2vmr     =    0.20946_wp !< O2
  REAL(wp) :: n2ovmr    =  309.5e-09_wp !< N20
  REAL(wp) :: cfcvmr(2) = (/252.8e-12_wp,466.2e-12_wp/)  !< CFC 11 and CFC 12
  !
  ! 2.0 Non NAMELIST global variables and parameters
  ! --------------------------------
  REAL(wp), PARAMETER :: twopi    = 2.0_wp*pi
  REAL(wp), PARAMETER :: deg2rad  = pi/180.0_wp
  REAL(wp), PARAMETER :: ch4_v(3) = (/1.25e-01_wp,  683.0_wp, -1.43_wp  /)
  REAL(wp), PARAMETER :: n2o_v(3) = (/1.20e-02_wp, 1395.0_wp, -1.43_wp/)
  !
  ! --- radiative transfer parameters
  !
  REAL(wp), PARAMETER :: cemiss = 0.996_wp  !< LW Emissivity Factor
  REAL(wp), PARAMETER :: diff   = 1.66_wp   !< LW Diffusivity Factor
  REAL(wp), PARAMETER :: rad_undef = -999._wp
  REAL(wp), PARAMETER :: zepzen = 1.e-10_wp 
  !
  !
  REAL(wp) :: solc = 1361.371_wp            !< default solar constant [W/m2] for
                                            !  AMIP-type CMIP5 simulation
  !++hs
  REAL(wp) :: psct                          !< local (orbit relative and possibly
  !                                            time dependent) solar constant
  REAL(wp) :: psctm                         !< orbit and time dependent solar constant for radiation time step
  REAL(wp) :: ssi_factor(14)                !< fraction of TSI in the 14 RRTM SW bands
  !--hs
  REAL(wp) :: flx_ratio_cur, flx_ratio_rad
  REAL(wp) :: decl_sun_cur                  !< solar declination at current time step
  !
  ! 3.0 Variables computed by routines in mo_radiation (export to submodels)
  ! --------------------------------
  !
  REAL(wp) :: co2mmr, ch4mmr, o2mmr, n2ommr                ! setup_radiation

  public solar_parameters

contains

  !---------------------------------------------------------------------------
  !>
  !! @brief Scans a block and fills with solar parameters
  !! 
  !! @remarks: This routine calculates the solar zenith angle for each
  !! point in a block of data.  For simulations with no dirunal cycle 
  !! the cosine of the zenith angle is set to its average value (assuming 
  !! negatives to be zero and for a day divided into nds intervals).  
  !! Additionally a field is set indicating the fraction of the day over 
  !! which the solar zenith angle is greater than zero.  Otherwise the field 
  !! is set to 1 or 0 depending on whether the zenith angle is greater or 
  !! less than 1. 
  !
  SUBROUTINE solar_parameters(decl_sun, dist_sun, time_of_day                  &
       &                     ,sinlon, sinlat, coslon, coslat                   &
       &                     ,flx_ratio, cos_mu0, daylght_frc)

    REAL(wp), INTENT(in)  :: &
         decl_sun,           & !< delination of the sun
         dist_sun,           & !< distance from the sun in astronomical units
         time_of_day,        & !< time_of_day (in radians)
         sinlon(:,:),        & !< sines of longitudes
         sinlat(:,:),        & !< and latitudes
         coslon(:,:),        & !< cosines of longitudes
         coslat(:,:)           !< and latitudes
    REAL(wp), INTENT(out) :: &
         flx_ratio,          & !< ratio of actual to average solar constant
         cos_mu0(:,:),       & !< cos_mu_0, cosine of the solar zenith angle
         daylght_frc(:,:)      !< daylight fraction (0 or 1) with diurnal cycle

    INTEGER     :: i, j
    REAL(wp)    :: zen1, zen2, zen3, z1, z2, z3, xx

    LOGICAL, SAVE      :: initialized = .FALSE.
    INTEGER, PARAMETER :: nds = 128 !< number of diurnal samples
    REAL (wp), SAVE :: cosrad(nds), sinrad(nds)
    REAL (wp)       :: xsmpl(nds), xnmbr(nds)

    IF (.NOT.initialized) THEN
      DO i = 1, nds
        xx = twopi*(i-1.0_wp)/nds
        sinrad(i) = SIN(xx)
        cosrad(i) = COS(xx)
      END DO
      initialized = .TRUE.
    END IF

    flx_ratio = 1.0_wp/dist_sun**2
    zen1 = SIN(decl_sun)
    zen2 = COS(decl_sun)*COS(time_of_day)
    zen3 = COS(decl_sun)*SIN(time_of_day)

    IF (lrce) THEN        ! radiation calculation in RCE configuration (same equatorial zenith
                          ! angle for all gridpoints)
       IF (ldiur) THEN                  ! all grid points have diurnal cycle of equator
          cos_mu0(:,:)     = -zen2
          daylght_frc(:,:) = 1.0_wp
          WHERE (cos_mu0(:,:) < 0.0_wp)
             cos_mu0(:,:)     = 0.0_wp
             daylght_frc(:,:) = 0.0_wp
          END WHERE
       ELSE
          DO j = 1, SIZE(cos_mu0,2)
             DO i = 1, SIZE(cos_mu0,1)

                xsmpl(:) = 0.7854_wp    ! solar weighted zenith angle of equator
                xnmbr(:) = 1.0_wp
                WHERE (xsmpl(:) < EPSILON(1.0_wp))
                   xsmpl(:) = 0.0_wp
                   xnmbr(:) = 0.0_wp
                END WHERE

                cos_mu0(i,j)     = SUM(xsmpl(:))
                daylght_frc(i,j) = SUM(xnmbr(:))
             END DO
          END DO

          WHERE (daylght_frc(:,:) > EPSILON(1.0_wp))
             cos_mu0(:,:)     = cos_mu0(:,:)/daylght_frc(:,:)
             daylght_frc(:,:) = daylght_frc(:,:)/nds
          END WHERE
       END IF
    ELSE           ! normal radiation calculation
     IF (ldiur) THEN
       cos_mu0(:,:)     =  zen1*sinlat(:,:)                &
            &             -zen2*coslat(:,:)*coslon(:,:)    &
            &             +zen3*coslat(:,:)*sinlon(:,:) 
       daylght_frc(:,:) = 1.0_wp
       WHERE (cos_mu0(:,:) < 0.0_wp)
         cos_mu0(:,:)     = 0.0_wp
         daylght_frc(:,:) = 0.0_wp
       END WHERE
     ELSE
       DO j = 1, SIZE(cos_mu0,2)
         DO i = 1, SIZE(cos_mu0,1)
 
           z1 =  zen1*sinlat(i,j)
           z2 = -zen2*coslat(i,j)
           z3 =  zen3*coslat(i,j)

           xsmpl(:) = z1 + z2*cosrad(:) + z3*sinrad(:)
           xnmbr(:) = 1.0_wp
           WHERE (xsmpl(:) < EPSILON(1.0_wp))
             xsmpl(:) = 0.0_wp
             xnmbr(:) = 0.0_wp
           END WHERE

           cos_mu0(i,j)     = SUM(xsmpl(:))
           daylght_frc(i,j) = SUM(xnmbr(:))
         END DO
       END DO
 
       WHERE (daylght_frc(:,:) > EPSILON(1.0_wp))
         cos_mu0(:,:)     = cos_mu0(:,:)/daylght_frc(:,:)
         daylght_frc(:,:) = daylght_frc(:,:)/nds
       END WHERE
     END IF
    END IF

  END SUBROUTINE solar_parameters

END MODULE mo_radiation_parameters
