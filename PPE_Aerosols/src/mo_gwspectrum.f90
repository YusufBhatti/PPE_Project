!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_gwspectrum

  !========================================================================
  !
  ! module *mo_gwspectrum* contains:
  ! A) controls variables for the set up of the Hines parameterization 
  !    for a gravity wave spectrum and its source spectrum. 
  ! B) Some internal switches and constants specific for the 
  !    Hines gw parameterization.
  !
  !
  ! ----------------------------------------------------------------------
  !
  ! A) Parameters and switches in namelist gwsctl
  ! =============================================
  ! see subroutine setgws for setup of namelist and modifications
  ! read from job control file 
  !
  ! Switches:
  ! ---------
  !
  ! lextro     : true for  hines' doppler spreading extrowave 
  !              parameterization (Hines, 1997a,b)
  !
  ! lfront     : true for gw emerging from fronts and background
  !              (Charron and Manzini, 2002)  
  !
  ! lozpr      : true for background enhancement associated with 
  !              precipitation (Manzini et al., 1997)
  !
  ! iheatcal   : switch for activating upper atmosphere processes
  !              iheatcal = 1 to calculate heating rates and diffusion 
  !                           coefficient.
  !              iheatcal = 0  only momentum flux deposition
  !
  ! The following parameter is set in setgws depending on the model
  ! resolution, unless defined in the namelist.
  !
  ! lrmscon_lat: true for latitude dependent rmscon as defined
  !                  in setgws (default) or gwsctl through rmscon_lo, rmscon_hi,
  !                  lat_rmscon_lo, and lat_rmscon_hi
  !              false: rmscon uniform over all latitudes as defined
  !                  by default or in gwsctl  
  !              attention: may be overwritten if lfront or lozpr true
  ! 
  !
  ! Parameters:
  ! -----------
  !
  !
  ! rmscon        : root mean square gravity wave wind at lowest level (m/s)
  !
  ! emiss_lev     : number of levels above the ground at which gw
  !                 are emitted (attention: this is vertical resolution
  !                 dependent. Must be generalized)
  !
  ! kstar         : typical gravity wave horizontal wavenumber (1/m)
  !
  ! m_min         : minimum bound in  vertical wavenumber (1/m)
  !               
  ! rms_front     : rms frontal gw wind at source level  (m/s)
  !
  ! front_thres   : minimum value of the frontogenesis function
  !                 for which gw are emitted from fronts [(K/m)^2/hour]
  !
  !
  !  pcrit        : critical precipitation value (mm/d) above which 
  !                 rms gw wind enhancement is applied
  !
  !  pcons        : adimensional facto rfor  background enhancement 
  !                 associated with precipitation 
  !
  ! The following parameters are set in setgws depending on the model
  ! resolution, unless defined in the namelist.
  !
  !  lat_rmscon_lo : latitude below which tropical GW source is used
  !  lat_rmscon_hi : latitude above which tropical GW source is used
  !  rmscon_lo     : tropical GW source parameter
  !  rmscon_hi     : extratropical GW source parameter
  ! 
  ! B) Internal swithches and constants 
  ! ===================================
  ! 
  !  naz        : actual number of horizontal azimuths used
  !               (code set up presently for only naz = 4 or 8)
  !
  !  slope      : slope of incident vertical wavenumber spectrum
  !               (code set up presently for slope 1., 1.5 or 2.)
  !               (if m_min not zero, code set up only for ??)
  !
  !  f1         : "fudge factor" used in calculation of trial value of
  !                azimuthal cutoff wavenumber m_alpha (1.2 <= f1 <= 1.9)
  !
  !  f2         : "fudge factor" used in calculation of maximum
  !                permissible instabiliy-induced cutoff wavenumber 
  !                (0.1 <= f2 <= 1.4)
  !
  !  f3         : "fudge factor" used in calculation of maximum 
  !               permissible molecular viscosity-induced cutoff wavenumber 
  !                m_sub_m_mol 
  !
  !  f5         : "fudge factor" used in calculation of heating rate
  !                (1 <= f5 <= 3)
  !
  !  f6         : "fudge factor" used in calculation of turbulent 
  !                diffusivity coefficient 
  !
  !  ksmin      : additional parameter to define a latitudinally varying kstar
  !                (ksmin = 1.e-5 [1/m] used maecham4)        
  !  ksmax      : additional parameter to define a latitudinally varying kstar
  !                (ksmax = 1.e-4  [1/m] used maecham4) 
  !                kstar = ksmin/( coslat+(ksmin/ksmax) ) as for maecham4
  !
  !  icutoff    : switch  for activating exponenetial damp of gwd, 
  !               heating and diffusion arrays above alt_cutoff
  !               icutoff = 1 : The exponentially damp is on
  !               icutoff = 0 : The arrays are not modified (default)
  !
  !  alt_cutoff : altitude (meters) above which exponential decay applied
  !
  !  smco       : smoother used to smooth cutoff vertical wavenumbers
  !               and total rms winds before calculating drag or heating.
  !                 (==> a 1:smco:1 stencil used; smco >= 1.)
  !
  !  nsmax      : number of times smoother applied ( >= 1)
  !
  !===========================================================================

  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  PUBLIC

  ! Switches in namelist gwsctl
  ! ---------------------------

  LOGICAL :: lextro   = .TRUE.
  LOGICAL :: lfront   = .FALSE.
  LOGICAL :: lozpr    = .FALSE.
  INTEGER :: iheatcal =   1
  LOGICAL :: lrmscon_lat

  ! Parameters in namelist gwsctl
  ! ----------------------------

  REAL(wp) :: rmscon     = 1.0_wp

  INTEGER  :: emiss_lev

  REAL(wp) :: kstar      = 5.0e-5_wp        ! = 2*pi/(126 km)
  REAL(wp) :: m_min      = 0.0_wp

  REAL(wp) :: rms_front   = 2.0_wp
  REAL(wp) :: front_thres = 0.12_wp              ! default for T42

  REAL(wp) :: pcrit       = 5.0_wp
  REAL(wp) :: pcons       = 4.75_wp

  REAL(wp) :: lat_rmscon_lo
  REAL(wp) :: lat_rmscon_hi
  REAL(wp) :: rmscon_lo
  REAL(wp) :: rmscon_hi

  !---------------------------------------------------------------
  ! Internal switches and constants 
  ! --------------------------------

  INTEGER :: naz   = 8

  REAL(wp) :: slope = 1.0_wp

  REAL(wp) :: f1    = 1.5_wp 
  REAL(wp) :: f2    = 0.3_wp 
  REAL(wp) :: f3    = 1.0_wp 
  REAL(wp) :: f5    = 1.0_wp 
  REAL(wp) :: f6    = 0.5_wp   

  REAL(wp) :: ksmin = 1.e-5_wp       
  REAL(wp) :: ksmax = 1.e-4_wp       

  INTEGER  :: icutoff    = 0   
  REAL(wp) :: alt_cutoff = 105.e3_wp

  REAL(wp) :: smco       = 2.0_wp      !  (test value: smco = 1.0)
  INTEGER  :: nsmax      = 5           !  (test value: nsmax = 2)

  REAL(wp),DIMENSION(:),ALLOCATABLE :: rmscon_lat  ! latitude dependent rmscon

  !-----------------------------------------------------------------

END MODULE mo_gwspectrum

