! ndgctl.inc
! 
!     ------------------------------------------------------------------
NAMELIST /ndgctl/ &
         lnudgini,          &! data synchronisation
         lnudgcli,          &! use nudging data cyclic
         dt_nudg_start,     &! nudging start date
         dt_nudg_stop,      &! nudging stop date
         ndg_file_stp,      &! nudging data, temperature and log surface pressure
         ndg_file_div,      &! nudging data, divergence
         ndg_file_vor,      &! nudging data, vorticity
         nsstinc,           &! sst handling, read new sst after NSSTINC hours
         nsstoff,           &! read sst at hour NSSTOFF first
         ndg_freez,         &! sea ice mask detection value
         ndg_file_sst,      &! sst data file
         lnudgimp,          &! nudging method, implicit/explicit
         lnudgpat,          &! pattern nudging
         lnudgfrd,          &! define place for NMI filter
         nudglmin,          &! vertical level separation, index of uppermost layer
         nudglmax,          &! index of deepest layer
         nudgtrun,          &! spectral domain selection
         nudgsmin,          &! lowest nudged wavenumber
         nudgsmax,          &! highest nudged wavenumber
         nudgp,             &! relaxation time definition for surface parameter
         nudgt,             &! relaxation time definition for temperature
         nudgd,             &! relaxation time definition for divergence
         nudgv,             &! relaxation time definition for vorticity
         ltintlin,          &! time interpolation linear/cubic spline
         ldamplin,          &! change relaxation between synoptic times
         nudgdamp,          &! 1= no damping, 0= full damping
         nudgdsize,         &! environment radius for damping
         lnudgdbx,          &! diagnostics
         lnudgwobs,         &! store reference fields
         lsite,             &! tendency diagnostics,
         ndg_file_nc,       &! netcdf-file of the nudging data
         inudgformat         ! switch to specify format of nudging data
!     ------------------------------------------------------------------
