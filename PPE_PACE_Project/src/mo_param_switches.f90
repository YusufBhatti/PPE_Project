!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_param_switches

 USE mo_kind,          ONLY: dp

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
  LOGICAL :: lclmi_progn  !   true for prognostic cloud droplet activation
!>>SF
  INTEGER :: nclmi_progn  ! allows to fine-tune the 2-moment cloud microphysics calculations in
!                           mo_cloud_micro_2m
!                           WARNING: this is meant for special sensitivity studies!!
!                                    Do not modify if you don't know about it
!                           0: off (both CDNC and ICNC are diagnostic variables)
!                              WARNING!! this is *not* equivalent to lclmi_progn = .false. because
!                              this still uses cloud_cdnc_icnc.f90 (with a different autoconversion
!                              scheme as compared to using cloud.f90).
!                           1: CDNC is a prognostic var, ICNC is a diagnostic var
!                           2: CDNC is a diagnotic var,  ICNC is a prognostic var
!                           3: CDNC and ICNC are prognostic vars [DEFAULT, normal case]         
!<<SF
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
!>>UP emulator
  LOGICAL :: lemuphase_nic_cirrus      ! Switch to phase the cirrus scheme in/out (between nic_cirrus = 1 and = 2)
  REAL(dp) :: eta_emu_nic_cirrus       ! Phase cirrus scheme 1/2 in/out by this much, between 0 and 1
  LOGICAL :: lemuphase_riming          ! Switch to phase riming in/out, issue #763
  REAL(dp) :: eta_emu_riming           ! Phase riming in/out by this much, between 0 and 2, issue #763
  LOGICAL :: lemuphase_wbf             ! Switch to phase WBF in/out, issue #765
  REAL(dp) :: eta_emu_wbf              ! Phase WBF in/out by this much, between 0 and 1 (!), issue #765
  LOGICAL :: lemuphase_icaggr          ! Switch to phase aggregation of ice crystals in/out, issue #757
  REAL(dp) :: eta_emu_icaggr           ! Phase aggregation in/out by this much, between 0 and 2, issue #757
  LOGICAL :: lemuphase_icaccr          ! Switch to phase accretion of ice crystals to snow in/out, issue #775
  REAL(dp) :: eta_emu_icaccr           ! Phase accretion in/out by this much, between 0 and 2, issue #775
  LOGICAL :: lemuphase_cdaccr          ! Switch to phase accretion of cloud droplets in/out, issue #805
  REAL(dp) :: eta_emu_cdaccr           ! Phase accretion in/out by this much, between 0 and 2, issue #805
  LOGICAL :: lemuphase_cdautc          ! Switch to phase autoconversion of cloud droplets in/out, issue #805
  REAL(dp) :: eta_emu_cdautc           ! Phase autoconversion in/out by this much, between 0 and 2, issue #805
  LOGICAL :: lemuphase_sci             ! Switch to phase self-collection of ice crystals to ice in/out, issue #779
  REAL(dp) :: eta_emu_sci              ! Phase self-collection of ice crystals in/out by this much, between 0 and 2, issue #779
  LOGICAL :: lemuphase_dep             ! Switch to phase deposition onto ice in/out, issue #801
  REAL(dp) :: eta_emu_dep              ! Phase deposition in/out by this much, between 0 and 2, issue #801
  LOGICAL :: lemuphase_icnucl          ! Switch to phase ice nucleation in the cirrus scheme in/out, issue #801
  REAL(dp) :: eta_emu_icnucl           ! Phase IC nucleation in/out by this much, between 0 and 2, issue #801
  LOGICAL :: lemuphase_sub_evp_falling ! Switch to phase sublimation/evaporation of falling ice, snow and rain in/out, issue #802
  REAL(dp) :: eta_emu_subfis           ! Phase sublimation of falling ice and snow in/out by this much, between 0 and 2, issue #802
  REAL(dp) :: eta_emu_evpr             ! Phase evaporation of rain in/out by this much, between 0 and 2, issue #802
  LOGICAL :: lemuphase_cdnuc           ! Switch to phase deposition onto ice in/out, issue #801
  REAL(dp) :: eta_emu_cdnuc            ! Phase deposition in/out by this much, between 0 and 2, issue #801
  LOGICAL :: lemuphase_sip             ! Switch to phase secondary ice production in/out, issue #805
  REAL(dp) :: eta_emu_sip              ! Phase SIP in/out by this much, between 0 and 2, issue #805
  LOGICAL :: lemuphase_mlt             ! Switch to phase melting of sedimenting ice and snow in/out, issue #802
  REAL(dp) :: eta_emu_mlt              ! Phase melting in/out by this much, between 0 and 2, issue #802
!<<UP
!>>MA 
  LOGICAL :: lemuphase_hetfrz         ! Switch to phase heterogeneous freezing in/out, issue #790
  REAL(dp) :: eta_emu_hetfrz          ! Phase heterogeneous freezing in/out by this much, between 0 and 1, issue #790
  LOGICAL :: lemuphase_fvic           ! Switch to phase fall velocity of ice crystal mass in/out, issue #795
  REAL(dp) :: eta_emu_fvic            ! Phase fall velocity of ice crystal mass in/out by this much, between 0 and 2, issue #795
  LOGICAL :: lctrl_WBF !switch WBF on/off
!<<MA
!>>UP #821
  LOGICAL :: lcmpsimpl_prescr ! Switch to prescribe various CMP deltas
                              ! instead of using the calucated values
  ! Following are integer switches for single processes:
  ! 0: computed values, default
  ! 1: one prescribed constant
  ! 2: prescried latlev profile
  INTEGER :: ncmpsimpl_prescr_rime
  INTEGER :: ncmpsimpl_prescr_icnucl
  INTEGER :: ncmpsimpl_prescr_subfis
  INTEGER :: ncmpsimpl_prescr_sci
  INTEGER :: ncmpsimpl_prescr_icaccr
  INTEGER :: ncmpsimpl_prescr_mlt
!<<UP


!>>UP #782.1
  LOGICAL :: ldetr_liquid ! switch to turn to only liquid detrainment at mixed-phase temperatures
  ! (otherwise (default), ice is detrained depending on WBF criterion)
!<<UP #782.1
!<<UP #783.3
  LOGICAL :: lpr_corr ! switch to turn on corrections of the process rates that
! serve to get reasonable correction terms
!<<UP #783.3
!>>UP #782.2
  REAL(dp) :: convicfactor
  LOGICAL :: ldetr_convicfactor
!<<UP #782.2
!>>UP #797
  LOGICAL :: lsedfix ! switch to turn on sedimentation fix for too large
! increases due to sedimentation, see #797
!<<UP #797
!>>UP: new diags for CD activation correction
  LOGICAL :: lnewdiags
  ! new diags
  LOGICAL :: lccnclimdiags
!<<UP: new diags
!>>UP FOR-ICE SLF #839
  LOGICAL :: lslf
!<<UP FOR-ICE SLF
  REAL(dp) :: rcd_vol_max         !maximum droplet volume radius [m]
                                  !UP: moved this here to make it
                                  !tunable via settings file
                                  ! Used for ldyn_cdnc_min

!davidn
  REAL(dp) :: tun47zinhomi
  REAL(dp) :: tun47zinhoml1
  REAL(dp) :: tun47zinhoml2
  REAL(dp) :: tun47zinhoml3
  REAL(dp) :: tun47zinpar
  REAL(dp) :: tun47cdncmin
  REAL(dp) :: tun47ccraut
  REAL(dp) :: tun47ccsaut
  REAL(dp) :: tun47entrpen
  REAL(dp) :: tun47entrscv
  REAL(dp) :: tun47cmfctop
  REAL(dp) :: tun47cprcon
  LOGICAL :: lepsfix !#680
!davidn

END MODULE mo_param_switches
