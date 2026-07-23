! gwsctl.inc
! 
!     ------------------------------------------------------------------
NAMELIST /gwsctl/ &
  lextro,         &! true for  hines' doppler spreading extrowave
                   !  parameterization (Hines, 1997a,b)
  lfront,         &! true for gw emerging from fronts and background
                   !  (Charron and Manzini, 2002)
  lozpr,          &! true for background enhancement associated with
                   !  precipitation (Manzini et al., 1997)
  iheatcal,       &! switch for activating upper atmosphere processes
                   !  iheatcal = 1 to calculate heating rates and diffusion
                   !               coefficient.
                   !  iheatcal = 0  only momentum flux deposition
  rmscon,         &! root mean square gravity wave wind at lowest level (m/s)
  lrmscon_lat,    &! true: latitude dependent rmscon as defined in 
                   !       setgws       
                   ! false: rmscon uniform over all latitudes as defined 
                   !        by default or in this namelist
                   ! attention: may be overwritten if lfront or lozpr are true 
  emiss_lev,      &! number of levels above the ground at which gw
                   !  are emitted (attention: this is vertical resolution
                   !  dependent. Must be generalized)
  kstar,          &! typical gravity wave horizontal wavenumber (1/m)
  m_min,          &! minimum bound in  vertical wavenumber (1/m)
  rms_front,      &! rms frontal gw wind at source level  (m/s)
  front_thres,    &! minimum value of the frontogenesis function
                   !  for which gw are emitted from fronts [(K/m)^2/hour]
  pcrit,          &! critical precipitation value (mm/d) above which
                   !  rms gw wind enhancement is applied
  pcons,          &! adimensional facto rfor  background enhancement
                   !  associated with precipitation
  lat_rmscon_lo,  &! latitude below which tropical GW source is used
                   ! (only used for lrmscon_lat=.true.)
  lat_rmscon_hi,  &! latitude above which tropical GW source is used
                   ! (only used for lrmscon_lat=.true.)
  rmscon_lo,      &! tropical GW source parameter
                   ! (only used for lrmscon_lat=.true.)
  rmscon_hi        ! extratropical GW source parameter
                   ! (only used for lrmscon_lat=.true.)
!     ------------------------------------------------------------------



