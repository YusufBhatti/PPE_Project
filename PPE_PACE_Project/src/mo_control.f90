!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_control

  ! Control variables for model housekeeping.
  !
  ! U. Schlese, DKRZ, December 1994
  ! A. Rhodin, MPI, January 1999,
  !      Subroutine m_control renamed to alloc_mods and moved from
  !      module mo_control to module m_alloc_mods.
  ! L. Kornblueh, MPI, June 1999,
  !      added nproca and nprocb for driving the parallel decomposition
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! I. Kirchner, MPI, December 2000, time control
  ! I. Kirchner, MPI, March 2001, revision
  ! M. Esch, MPI, September 2002, add switch for mixed layer ocean
  ! M. Esch, MPI, November  2003, add switch for scenario runs
  ! L. Kornblueh, MPI, April 2010, more switches for parallel I/O
  ! D. Popke, MPI, Nov 2013, Implementation of RCE (lrce)
  ! ------------------------------------------------------------------

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  SAVE    ! to assure that the scope of all variables is for the whole program run

  REAL(dp), POINTER :: vct(:)=>NULL() ! vertical coefficients table.

  INTEGER :: nproma   = -1  !   working dimension for grid-point computations
  INTEGER :: nproca   = 1   !   number of processors in set A
  INTEGER :: nprocb   = 1   !   number of processors in set B
  INTEGER :: nprocio  = 0   !   number of processors for I/O server

  INTEGER :: nm             !   max zonal wave number.
  INTEGER :: nn             !   max meridional wave number for m=0.
  INTEGER :: nk             !   max meridional wave number.
  INTEGER :: ngl            !   number of gaussian latitudes.
  INTEGER :: nlon           !   max number of points on each latitude line.
  INTEGER :: nlev           !   number of vertical levels.
  INTEGER :: nmp1           !   max zonal wave number + 1.
  INTEGER :: nnp1           !   max meridional wave number + 1.
  INTEGER :: nkp1
  INTEGER :: n2mp1          !   2 * (max zonal wave number + 1).
  INTEGER :: n4mp1          !   4 * (max zonal wave number + 1).
  INTEGER :: nlp2           !   max number of points per latitude line + 2.
  INTEGER :: nlevp1         !   *nlev+1.
  INTEGER :: nsp            !   number of spectral coefficients.
  INTEGER :: n2sp           !   2*number of spectral coefficients.
  INTEGER :: nhgl           !   (number of gaussian latitudes)/2.
  INTEGER :: nscan          !   current scan number.

  INTEGER :: nspace1        !   memory manager space for use of root task
  INTEGER :: nspace2        !   memory manager space for use of subtasks

  INTEGER :: nspadd         !   memory manager space increase

  INTEGER :: nvclev         !   number of levels with vertical coefficients.

  LOGICAL :: ldebug    = .FALSE. !   .true. for mass fixer diagnostics
  LOGICAL :: ldailysst = .FALSE. !   .true. for using daily SST and SIC
  LOGICAL :: lamip     = .FALSE. !   .true. for using variable sst
  LOGICAL :: ldiagamip = .FALSE. !   .true. for AMIP diagnostics
  LOGICAL :: lcouple   = .FALSE. !   .true. for a coupled run
  LOGICAL :: lcouple_co2 = .FALSE. !   .true. if coupled run includes co2
  LOGICAL :: lcouple_parallel = .FALSE. !   .true. if all PEs should communicate with oasis
  LOGICAL :: lnwp      = .FALSE. !   .false. for climate mode .true. for NWP mode
  LOGICAL :: lnudge    = .FALSE. !   .true. for Nudging mode
  LOGICAL :: lmidatm   = .TRUE.  !   .true. for middle atmosphere model version
  LOGICAL :: lmlo      = .FALSE. !   .true. for mixed layer ocean
  LOGICAL :: lmeltpond = .TRUE.  !   .true. for calculation of meltponds
  LOGICAL :: lmlo_ice  = .TRUE.  !   .true. for sea ice calculation in mixed layer ocean
  REAL(dp):: rmlo_depth = 50._dp !   default mixed layer depth
  REAL(dp):: earth_angular_velocity = 7.29212e-5_dp !! [rad/s] default value for model angular
                                                    !velocity, is set to namelist value once
                                                    !runctl is read in inictl.f90
  LOGICAL  :: lrce      = .FALSE. !   .true. if ECHAM is run in radiative-convective equilibrium
                                  !   configuration; switches to RCE radiation calculation,
                                  !   regulates pbl height calculation (vdiff.f90),
                                  !   sets oceans surface albedo of direct radiation to constant
  LOGICAL :: lfractional_mask = .FALSE.  !   .true. for fractional land sea mask

  LOGICAL :: lprint_m0 = .FALSE. !   .false. for printing t(m0) and timestep time in stepon 
  LOGICAL :: lnmi      = .FALSE. !   .true. normal mode initialisation
  LOGICAL :: ltdiag    = .FALSE. !   .true. run with additional tendency diagnostics
  LOGICAL :: lcolumn   = .FALSE. !   .true. column model mode
  LOGICAL :: lvctch    = .FALSE. !   .true. if column model has changed vct

  LOGICAL :: l_volc    = .FALSE. !   switch for volcanic forcing

  REAL(dp):: satoverpasstime = 10.5_dp ! satellite overpass time at 10.30 a.m. [s]

  LOGICAL :: lhd  ! activates hydrologic discharge (HD) model (defined in jsbach namelist jsbach_ctl)


  ! Spectral and grid point initial files.
  INTEGER, PARAMETER :: nisp  = 23   ! unit for initial spectral fields.
  INTEGER, PARAMETER :: nigp  = 24   ! unit for initial grid point fields.

  ! Climate sea surface temperature and sea ice annual cycle file
  INTEGER, PARAMETER :: nist  = 20   ! unit for surf.temp. file
  INTEGER, PARAMETER :: nice  = 96   ! unit for amip ice file

  ! Climate flux correction file
  INTEGER, PARAMETER :: nflu  = 42   ! unit for flux correction file

  ! Climate leaf area index and vegetation ratio annual cycle file
  INTEGER, PARAMETER :: nvltcl   = 90 ! unit for climate leaf area index
  INTEGER, PARAMETER :: nvgratcl = 91 ! unit for climate vegetation ratio

  ! Climate land surface temperature annual cycle file
  INTEGER, PARAMETER :: ntslcl   = 92 ! unit for climate land surface temperature

  INTEGER, PARAMETER :: nsurface = 79 ! unit for boundary atm-surface jsbach

  ! optional files
  INTEGER :: ndiahdf = 10 ! I/O unit for hdiff diagnostics.

  ! Switches
  LOGICAL :: ltimer            = .FALSE. ! to use timer
  LOGICAL :: ldebugio          = .FALSE. ! to debug IO
  LOGICAL :: ldebugmem         = .FALSE. ! to debug memory
  LOGICAL :: ltctest           = .FALSE. ! to test time control
  LOGICAL :: ldebugs           = .FALSE. ! to enable debug stream
  LOGICAL :: lindependent_read = .FALSE. ! to read initial/restart input on each MPI rank separatly 
  LOGICAL :: lcollective_write = .FALSE. ! to write restart in parallel on each MPI rank
  LOGICAL :: ldebugcpl         = .FALSE. ! to debug coupling
!>>SF: to switchable internal rerun
  LOGICAL :: lforcererun       = .TRUE. ! to switch on/off internal rerun:
                                        !  - true means internal rerun are forced (default in ECHAM6)
                                        !  - false means internal rerun are bypassed (compulsory in ECHAM-HAMMOZ)
!<<SF 
  LOGICAL :: lyaxt_transposition = .TRUE. ! to use YAXT for global transpositions

  ! output redirection
  INTEGER :: stdout_redir = 0
  INTEGER :: stderr_redir = 0

END MODULE mo_control
