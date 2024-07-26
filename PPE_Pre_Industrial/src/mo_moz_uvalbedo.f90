!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_uvalbedo
!!
!! \brief
!!  Compute UV albedo for MOZART photolysis frequency calculations and provide 
!!  it as boundary condition
!!
!! \author Martin G. Schultz (FZ-Juelich)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - MGS: original version based on T. Laepple's code (2009-08-29)
!!
!! \limitations
!!  none
!!
!! \details
!!  The UV albedo is parameterized according to Laepple et al., 2005 using two
!!  MODIS derived albedo maps for snow-covered and snow-free land.
!!  Technically, it is realized as a user-defined boundary condition. The two
!!  maps are read from file (monthly time resolution) and the combined value
!!  (depending on current snow cover) is stored as a "module" type boundary
!!  condition.
!!  User control: via uvalbedo_file in mozctl.
!!
!! \bibliographic_references
!!  - Laepple T; Schultz MG; Lamarque JF; et al., J. Geophys. Res., 110(D11), doi: 10.1029/2004JD005463, 2005.
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_moz_uvalbedo

  USE mo_kind,               ONLY: dp

  IMPLICIT NONE

  PRIVATE

  SAVE

! module routines
  PUBLIC :: moz_uvalbedo_init        ! define boundary conditions
  PUBLIC :: moz_uvalbedo             ! compute combined UV albedo and store as boundary condition

  PUBLIC :: ibc_uvalbedo             ! index of UV albedo boundary condition.
                                     ! ibc_uvalbedo+1 is used for green map
                                     ! ibc_uvalbedo+2 is used for white map

! Variable declarations
!++sschr #379: added missing initialization
  INTEGER            :: ibc_uvalbedo = -1
!++sschr #379


  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! moz_uvalbedo_init: initializes the boundary condition structures
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Externals:
!! <ol>
!! <li>none
!! </ol>
!!
!! @par Notes
!! 
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE moz_uvalbedo_init (uvalbedo_file)

  USE mo_boundary_condition, ONLY: bc_nml, bc_define, BC_EVERYWHERE, BC_BOTTOM, BC_REPLACE, BC_SPECIAL, BC_RELAX
  USE mo_external_field_processor, ONLY: EF_FILE, EF_MODULE, EF_LONLAT, EF_IGNOREYEAR

  CHARACTER (LEN=*), INTENT(in)       :: uvalbedo_file     ! filename from mozctl namelist

  TYPE (bc_nml)                       :: bc_struc
  INTEGER                             :: idum


  ! set parameters for UV albedo boundary condition
  bc_struc%ef_type      = EF_MODULE
  bc_struc%ef_geometry  = EF_LONLAT       !! only for the other two bcs ??
  bc_struc%bc_mode      = BC_SPECIAL
  bc_struc%bc_domain    = BC_EVERYWHERE
  ! define UV albedo boundary condition
  ibc_uvalbedo = bc_define('MOZART UV albedo', bc_struc,2,.TRUE.)

  ! set parameters for the green map
  bc_struc%ef_type      = EF_FILE
  bc_struc%ef_template  = TRIM(uvalbedo_file) 
  bc_struc%ef_varname   = 'agreen'
  bc_struc%ef_geometry  = EF_LONLAT
  bc_struc%ef_timedef   = EF_IGNOREYEAR
  bc_struc%bc_mode      = BC_REPLACE
  idum = bc_define('MOZART snow-free UV albedo', bc_struc,2,.TRUE.)

  ! set parameters for the white map
  bc_struc%ef_varname   = 'awhite'
  idum = bc_define('MOZART snow-covered UV albedo', bc_struc,2,.TRUE.)

END SUBROUTINE moz_uvalbedo_init


SUBROUTINE moz_uvalbedo (kbdim, kproma, krow, pfrl, pfrw, pfri, pfrg)

  USE mo_boundary_condition, ONLY: bc_set, bc_apply
  USE mo_memory_g3b,         ONLY: sn       !! deprecated: better to use vphysc stream...###
  USE mo_geoloc,             ONLY: philat_2d
  USE mo_time_control,       ONLY: time_step_len
  USE mo_moz_diag,           ONLY: dpalb    ! diagnostics pointer

  INTEGER, INTENT(in)       :: kbdim, kproma, krow
  REAL(dp), INTENT(in)      :: pfrl(kbdim), pfrw(kbdim),      &
                               pfri(kbdim), pfrg(kbdim)         ! land, water and ice fractions

  REAL(dp)                  :: zfrw(kbdim)                      ! total water fraction incl. seaice
  REAL(dp)                  :: zgreen(kproma), zwhite(kproma)   ! MODIS albedo values
  REAL(dp)                  :: zl(kbdim), zw(kbdim), zi(kbdim)  ! albedo for land, water, ice
  REAL(dp)                  :: zalb(kproma)                     ! combined albedo value

  REAL(dp), PARAMETER       :: alb_ocean = 0.07_dp              ! open ocean albedo
  REAL(dp), PARAMETER       :: alb_arctic = 0.78_dp             ! Arctic ice albedo
  REAL(dp), PARAMETER       :: alb_antarc = 0.89_dp             ! Antarctic ice albedo

  ! Note: we use the sea ice fraction on water as pfri here, because otherwise we produce
  ! "holes" of low albedo especially in Northern Canada, where the glacier and snow masks
  ! don't work well for grid boxes that are partially land and partially ocean.
  ! For simplicity we apply the seaice albedo value to the entire grid box wherever the 
  ! seaice fraction is > 0.5

  zalb(:) = 0._dp
  zfrw(1:kproma) = 1._dp - pfrl(1:kproma)

  ! get values for green and white albedo
  CALL bc_apply(ibc_uvalbedo+1, kproma, krow, zgreen)
  CALL bc_apply(ibc_uvalbedo+2, kproma, krow, zwhite)
  ! compute combined UV albedo: 
  ! first derive land albedo value from "green" or "white" map
  ! - threshold is snow depth of 1 cm -- old ECHAM5-MOZ had 0.008
  ! - use white map for glacier fraction of land
  zl(1:kproma) = zgreen(1:kproma)
  WHERE (sn(1:kproma,krow) > 0.00999_dp)  zl(1:kproma) = zwhite(1:kproma)
  zl(1:kproma) = pfrg(1:kproma)*zwhite(1:kproma) + (1._dp-pfrg(1:kproma))*zl(1:kproma)
  ! next, compute water albedo
  ! - ocean value: constant at 0.07
  zw(1:kproma) = alb_ocean
  ! - seaice value: 0.78 for Arctic, 0.89 for Antarctic
  zi(:) = alb_arctic
  WHERE (philat_2d(1:kproma,krow) < 0._dp)  zi(1:kproma) = alb_antarc
  ! - set water albedo to seaice value if seaice fraction exceeds 30%
  WHERE (pfri(1:kproma) > 0.3_dp)  zw(1:kproma) = zi(1:kproma)
  ! Combine albedo values
  ! Note: pfrl+zfrw = 1 per definition; differs from notation in physc!!
  zalb(1:kproma) = pfrl(1:kproma)*zl(1:kproma) + zfrw(1:kproma)*zw(1:kproma)
  ! Increase albedo for grid boxes that are predominantly land but have 
  ! a large seaice fraction ("Canada fix")
  WHERE (pfrl(1:kproma) > 0.5_dp .AND. pfri(1:kproma) > 0.75_dp)     &
     zalb(1:kproma) = 0.2_dp*zl(1:kproma) + 0.8_dp*zw(1:kproma)

  CALL bc_set(ibc_uvalbedo, kproma, krow, zalb)

  ! save as diagnostics
!+sschr #578
  IF (associated(dpalb)) dpalb(1:kproma, krow) = zalb(1:kproma)
!-sschr

END SUBROUTINE moz_uvalbedo

END MODULE mo_moz_uvalbedo
