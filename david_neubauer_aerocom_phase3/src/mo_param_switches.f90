!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_param_switches

  !>>DN AeroCom
  USE mo_kind,               ONLY : dp
  !<<DN AeroCom

  IMPLICIT NONE

  ! M.A. Giorgetta, March 2000, lrad added
  !
  ! ----------------------------------------------------------------
  !
  ! module *mo_param_switches* switches related to the parameterisations of
  !                            diabatic processes except for radiation.
  !
  ! ----------------------------------------------------------------

  LOGICAL :: lphys        !   *true for parameterisation of diabatic processes.
  LOGICAL :: lrad         !   *true for radiation.
  LOGICAL :: lvdiff       !   *true for vertical diffusion.
  LOGICAL :: lcond        !   *true for large scale condensation scheme.
  LOGICAL :: lsurf        !   *true for surface exchanges.
  LOGICAL :: lconv        !   *true to allow convection
  LOGICAL :: lgwdrag      !   *true for gravity wave drag scheme
  LOGICAL :: lice         !   *true for sea-ice temperature calculation
  LOGICAL :: lconvmassfix !   *false for switching off aerosol mass fixer in conv

  INTEGER :: iconv = 1    !   *1,2,3 for different convection schemes
  INTEGER :: icover = 1   !   *1 for default cloud cover scheme
                          !   *2 for prognostic cloud cover scheme

!++mgs : new switches for interactive cloud scheme
  LOGICAL :: lcdnc_progn  !   true for prognostic cloud droplet activation
  INTEGER :: ncd_activ    !   type of cloud droplet activation scheme (see physctl.inc)
  INTEGER :: nactivpdf    !   sub-grid scale PDF of updraft velocities (see physctl.inc) !ZK
!>>SF changed the former logical lice_supersat into an integer characterizing the cirrus scheme
  INTEGER :: nic_cirrus   !   type of cirrus scheme (see physctl.inc)
!<<SF
  INTEGER :: nauto        !   autoconversion scheme    (1,2)
  LOGICAL :: lsecprod     ! switch to take into account secondary ice production (see cloud_cdnc_icnc.f90) #251
                          ! lsecprod = .FALSE. (default, no secondary ice production)
                          !          = .TRUE.  (secondary ice prod)
  LOGICAL :: lorocirrus   ! switch to take into account gravity waves updraft velocity in 
                          ! ice crystal number concentration calculation (--> orographic cirrus clouds)
                          ! lorocirrus = .FALSE. (default, no orographic cirrus clouds)
                          !            = .TRUE.  (orographic cirrus clouds on)
!>>SF #475
  LOGICAL :: ldyn_cdnc_min ! switch to turn on the dynamical setting of the minimum cloud droplet number concentration
                           ! ldyn_cdnc_min = .FALSE. (default, fixed minimum CDNC)
                           !                 .TRUE.  (dynamical min CDNC)
!<<SF #475
!>>SF #589
 INTEGER :: cdnc_min_fixed ! fixed value for min CDNC in cm-3 (used when ldyn_cdnc_min is FALSE)
                           ! Warning! So far only values of 40 or 10 are accepted.
!<<SF #589

 !>>DN
 REAL(dp) :: ac_scale_activation
 REAL(dp) :: ac_scale_autoconversion
 REAL(dp) :: ac_scale_accretion
 REAL(dp) :: ac_scale_KK_LWP_exponent
 !<<DN
!>>DN ambient AOD per species
 LOGICAL  :: lnoh2oaod
!<<DN

 
END MODULE mo_param_switches
