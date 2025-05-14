!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! hamctl.inc
!!
!! \brief
!! General namelist for HAM
!!
!! \author Philip Stier (MPI-Met)
!! \author Martin G. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# Philip Stier (MPI-Met) - original code (YYYY)
!!   -# Martin G. Schultz (FZ Juelich) (2009)
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

NAMELIST /hamctl/ &
!--- microphysics scheme
  nham_subm,              & !< Choice of aerosol microphysics scheme
!--- aerosol processes
  nseasalt,               & !< Choice of the Sea Salt emission scheme
  npist,                  & !< Choice of the air-sea exchange scheme
  ndrydep,                & !< choice of dry deposition scheme (if ldrydep == .true.)
  nwetdep,                & !< choice of wet deposition scheme (if lwetdep == .true.)
  ndust,                  & !< choice of parameter sets for the BGC dust scheme
  naerorad,               & !< Choice for radiatively active aerosols
  laerocom_diag,          & !< Extended aerosol diagnostics
  nrad,                   & !< Radiation calculation (for each mode)
  nradmix,                & !< Mixing scheme for refractive indices (for each mode)
  nraddiag,               & !< Extended radiation diagnostics
  nraddry,                & !< dry aerosol radiation diagnostics
  nradbeta,               & !< LIDAR diagnostics
  lhetfreeze,             & !< switch to set heterogeneous freezing below 235K (cirrus scheme)
  nsoa,                   & !< Choice for Secondary Organic Aerosols
                            !! 0: no SOA scheme
                            !! 1: SOA scheme from O'Donnell et all, ACP 2011
                            !! 2: SOA scheme with VBS approach from Farina et al, JGR 2010
                            !!    (curr. SALSA only)
  nsoalumping,            & !< SOA lumping scheme to apply
  nlai_drydep_ef_type,    & !< Choice of lai external field type in the drydep scheme  !gf #244
  lscond,                 & !< Condensation of H2SO4 
  lscoag,                 & !< Coagulation 
  lgcr,                   & !< Galactic cosmic ray ionization 
  nsolact,                & !< Solar activity parameter
  lmass_diag,             & !< Mass diagnostics switch
  nccndiag,               & !< (C)CN diagnostics at fixed supersaturations
                            !! 0: OFF
                            !! 1: 2D CCN diagnostics
                            !! 2: 3D CCN diagnostics
                            !! 3: 2D CCN + CN diagnostics
                            !! 4: 3D CCN + CN diagnostics
                            !! 5: 2D CCN + CN diagnostics + burdens
                            !! 6: 3D CCN + CN diagnostics + burdens
!--- output control
  burden_keytype,         & !< 0 : no output
                            !! 1 : by tracer
                            !! 2 : per species
                            !! 3 : per mode (default)
!--- boundary conditions (oxidant fields)
  bc_oh,                  & !< tracer boundary condition for OH
  bc_o3,                  & !< tracer boundary condition for O3
  bc_h2o2,                & !< tracer boundary condition for H2O2
  bc_no2,                 & !< tracer boundary condition for NO2
  bc_no3                    !< tracer boundary condition for NO3
