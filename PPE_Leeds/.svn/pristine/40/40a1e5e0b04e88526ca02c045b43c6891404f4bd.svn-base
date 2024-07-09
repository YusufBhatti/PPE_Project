!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_nudging_constants
!BOP
  !
  ! !MODULE: mo_nudging_constants (layer 4)

  ! !DESCRIPTION: 
  !
  ! contains internal constants and control variables of nudging module,
  ! the following features can be modified using the namelist {\it \&NDGCTL}
  ! \begin{enumerate}
  ! \item data synchronisation
  ! \item nudging data
  ! \item sst handling
  ! \item nudging method
  ! \item vertical level separation
  ! \item spectral domain selection
  ! \item relaxation time definition
  ! \item time interpolation
  ! \item diagnostics
  ! \end{enumerate}

  ! !REVISION HISTORY: 
  ! Ingo Kirchner, MPI Hamburg, April-2001
  ! Ingo Kirchner, MPI Hamburg, Aug-2002, revision
  ! Aiko Voigt, MPI Hamburg, March 2011, add inudgformat to chose between netcdf and cray-format
  !                                      0 --> CRAY binary format input
  !                                      1 --> GRIB (not implemented)
  !                                      2 --> netcdf format input

  ! !USES:
  USE mo_kind,       ONLY: dp

!EOP
!BOC
!BOX
  IMPLICIT NONE
!EOX
  CHARACTER(len=*), PARAMETER :: &
       ndg_version = 'Nudging Version E5R1.04 August-2002 (kirchner@dkrz.de)'

  ! FEATURE: data synchronisation ==============================================

  ! EXTERNAL
  LOGICAL :: lnudgini = .FALSE.  ! adjust date/time from nudging data
  !                              true  --> use first date from nudging data
  !                              false --> use information from history data
  !
  ! the adjustment (lnudgini=.true.) of the date is dependent on LRESUME
  ! a) LRESUME=.TRUE. .... the NEXT_DATE must be coded in the file name
  ! b) LRESUME=.FALSE. ... the DT_NUDG_START date must be coded in the file name

  LOGICAL :: lnudgcli = .FALSE.  ! read input data cyclic
  !                              true  --> year ignored, only month/day/time 
  !                                        important for synchronization
  !                              false --> read data continuously

  !  DT_NUDG_START  nudging start date (see mo_time_control)
  !  DT_NUDG_STOP   nudging stop date (see mo_time_control)
!BOX

  ! INTERNAL
  LOGICAL :: lnudgstop = .FALSE. ! test nudging stop time

!EOX
  ! FEATURE: nudging data ======================================================

  ! the file name template should contain the date information in the
  ! form e.g. %y4%m2 for year with four and months with two digits

  ! EXTERNAL
  !  NDG_FILE_VOR    file name templates for vorticity (see mo_nudging_io)
  !  NDG_FILE_DIV    file name templates for divergence (see mo_nudging_io)
  !  NDG_FILE_STP    file name templates for temp and lnps (see mo_nudging_io)


  ! FEATURE: sst handling ======================================================

  ! EXTERNAL
  !  NSSTINC       read new sst after NSSTINC hours, zero means no usage of 
  !                nudging sst (see mo_nudging_sst)
  !  NSSTOFF       offset in hours relativ to 00 UTC (see mo_nudging_sst)
  !  NDG_FREEZ     temperatur point for ice mask coding (see mo_nudging_sst)
  !  NDG_FILE_SST  sst file name template (see mo_nudging_sst)


  ! FEATURE: nudging method ====================================================

  ! EXTERNAL
  LOGICAL :: lnudgimp = .TRUE.   ! choose between implicit or explicit method
  !                              true  --> implicit nudging (default)
  !                              false --> explicit nudging

  LOGICAL :: lnudgpat = .FALSE.  ! pattern nudging
  !                              true  --> assimilation of correlation pattern
  !                              false --> no effect
  !
  ! Remark 1: for pattern nudging only linear time interpolation is possible
  ! Remark 2: pattern nudging is not running in parallel mode

  LOGICAL :: lnudgfrd = .FALSE.  ! place for location of NMI filter
  !                              true  --> NMI filter performed in ReadOneBlock
  !                              false --> NMI filter after time interpolation
!BOX

  ! INTERNAL
  LOGICAL :: lnudg_run = .FALSE.  ! perform nudging

!EOX
  ! FEATURE: vertical level separation =========================================
  !
  ! the level domain will be defined by the external data amount
  ! negative values of nudging coefficients means no data in nudging file

  ! CONSTANT
  INTEGER, PARAMETER :: nmaxc = 80     ! maximal number of nudging weights for levels

  ! EXTERNAL 
  INTEGER :: nudglmin = 1        ! lowest level index (top->down)
  INTEGER :: nudglmax = nmaxc    ! highest level index
!BOX

  ! INTERNAL
  ! specific level selection
  INTEGER :: &
       ilev_v_min = 0 , ilev_v_max = 0, ino_v_lev = 0, &! vorticty
       ilev_d_min = 0 , ilev_d_max = 0, ino_d_lev = 0, &! divergence
       ilev_t_min = 0 , ilev_t_max = 0, ino_t_lev = 0   ! temperature
                                                        ! last level for lnps
  ! reduce number of input data
  LOGICAL :: &
       linp_vor = .FALSE., &! read vorticity
       linp_div = .FALSE., &! read divergence
       linp_tem = .FALSE., &! read temperatur
       linp_lnp = .FALSE.   ! read pressure

!EOX
  ! FEATURE: spectral domain selection =========================================

  ! definitions for the spectral domain which will be nudged

  ! CONSTANT
  !                    reduce the number of spectral coefficients nudged
  INTEGER, PARAMETER :: &
       NDG_WINDOW_ALL  = 0, &! use all coefficients of active wave no.
       NDG_WINDOW_CUT  = 1, &! cut off meri. index triangular  
       NDG_WINDOW_CUT0 = 2   ! cut off, except wave 0
  !
  INTEGER, PARAMETER :: max_wavenum = 106 !  highest modified wave number

  ! EXTERNAL
  INTEGER :: nudgtrun = NDG_WINDOW_ALL    ! type of wave index selection

  !                     boundaries of nudging window
  INTEGER :: nudgsmin = 0                 ! lowest modified wave number
  !                      0 .......... skip global average
  !                     -1 .......... incl. global mean
  INTEGER :: nudgsmax = max_wavenum       ! highest modified wave number
  !                     set level region
!BOX

  ! INTERNAL
  ! total spectral coefficient number used
  INTEGER :: nudgmin  = 1
  INTEGER :: nudgmax  = max_wavenum*max_wavenum  ! resolution dependent

!EOX
  ! FEATURE: relaxation time definition ========================================

  ! CONSTANT
  REAL(kind=dp), PARAMETER :: nfact = 1.e-5_dp  ! rescaling factor for coefficients
  !
  !                   unit of external nudging weigths [1.e-5/sec]
  !                   set nudging coefficients following Jeuken et al. (1996)
  !                   not optimal for all resolutions
  !

  ! EXTERNAL
  ! **** external defined nudging weights -> WEIGHT*NFACT
  !      for each level an individual weight can be defined

  ! coefficients from Jeuken et al. 1996
  !REAL(kind=dp)    :: &
  !     nudgd(nmaxc) =  5.0, &! nudging weight for divergence (5.6 hours)
  !     nudgv(nmaxc) = 10.0, &! nudging weight for vorticity (2.8 hours)
  !     nudgt(nmaxc) =  1.0, &! nudging weight for temperature (27.8 hours)
  !     nudgp        = 10.0   ! nudging weight for log sfc pressure (2.8 hours)
  !
  ! use DMI coefficients as default, corresponding relaxation time in braces
  REAL(kind=dp) :: &
       nudgd(nmaxc) = 0.5787_dp, &! nudging weight for divergence (48 hours)
       nudgv(nmaxc) = 4.6296_dp, &! nudging weight for vorticity (6 hours)
       nudgt(nmaxc) = 1.1574_dp, &! nudging weight for temperature (24 hours)
       nudgp        = 1.1574_dp   ! nudging weight for log sfc pressure (24 hours)


  ! FEATURE: time interpolation ================================================

  ! EXTERNAL
  LOGICAL :: ltintlin  = .TRUE.  ! time interpolation method
  !                              true  --> linear time interpolation
  !                              false --> non-linear time interpolation

  LOGICAL :: ldamplin  = .TRUE.  ! damping type of time interpolation
  !                              true  --> linear damping between two nudging times
  !                              false --> non-linear damping

  !          set time dependent damping/realaxation
  !          reduce nudging weight between the reference date/time (0...1)
  REAL(kind=dp) :: nudgdamp  = 1.0_dp  ! minimal fraction of the nudging weigth
                                    ! reached at times corresponding to the 
                                    ! influence radius

  !          influence radius near nudging times
  REAL(kind=dp) :: nudgdsize = 0.5_dp  ! nudging all the times


  ! FEATURE: diagnostics =======================================================

  ! EXTERNAL
  LOGICAL :: lnudgdbx = .FALSE.  ! debugging message control
  !                              false --> skip additional messages
  !                              true  --> additional debugging messages

  LOGICAL :: lnudgwobs = .FALSE. ! store reference fields
  !                              false --> no additional output of reference data
  !                              true  --> store reference data in model output

  LOGICAL :: lsite = .FALSE.     ! calculation of the Systematic Initial Tendency Error
  !                              true  --> additional output, memory
  !                              false --> no calculation
!BOX

  INTEGER, PUBLIC      :: ndunit  ! unit for special diagnostics
!EOX
!EOC
  INTEGER :: inudgformat = 0     ! switch to choose between CRAY format or netcdf input
  
END MODULE mo_nudging_constants
