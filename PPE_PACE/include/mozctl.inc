!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mozctl
!!
!! \brief
!!  Controls settings for MOZART chemistry
!!
!! \author Sebastian Rast (MPI-M)
!! \author Martin G. Schultz (FZ-Juelich)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - SR and MGS: original version (2005)
!!
!! \limitations
!!  none
!!
!! \details
!!
!! \bibliographic_references
!!  none
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

NAMELIST /mozctl/ &
  lchemsolv,                & ! activate chemistry solver
  lphotolysis,              & ! activate photolysis calculation
  lfastj,                   & ! use fastJ photolysis code instead of WACCM
  lfastjaero,               & ! use fastJ photolysis code without aerosols
  lfastjcloud,              & ! use fastJ photolysis code without clouds
  lstrathet,                & ! heterogeneous chemistry in stratosphere on/off
  ndrydep,                  & ! choice of dry deposition scheme (if ldrydep == .true.)
  nwetdep,                  & ! choice of wet deposition scheme (if lwetdep == .true.)
  lbc_species,              & ! list of species for tropospheric lower boundary condition
  ubc_species,              & ! list of species for stratospheric/mesospheric upper boundary condition
  out_species,              & ! list of species which shall be output in _tracer file
  budget_species,           & ! list of species for which budget diag. shall be run
  burden_species,           & ! list of species for which burden diagnostics is done
  photovars,                & ! list of variables in photo stream (or "all")
  uvalbedo_file,            & ! file name for UV albedo data (green and white)
  bc_ch4,                   & ! (lower) boundary condition for methane
  bc_lbc,                   & ! (lower) boundary conditions for tracers (strato. run/HAMMONIA)
  bc_ubc,                   & ! (upper) boundary condition for tracers (tropospheric run)
  bc_aerosol,               & ! boundary condition for aerosol (trophet run)
  ltrophet                    ! Tropospheric heterogeneous chemistry (if ltrophet == .true.)
!     ------------------------------------------------------------------

