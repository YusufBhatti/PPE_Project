!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_semi_lagrangian

  ! J. Olson, NCAR, unknown, Original version 
  ! J. Rosinski, NCAR, June 1992, Standardized
  ! D. Williamson, P. Rasch, NCAR, August 1992, Reviewed         
  ! L. Kornblueh, U. Schulzweida, MPI, May 1998, f90 version
  ! T. Diehl, DKRZ, July 1999, parallelization
  ! L. Kornblueh, MPI, July 2001, packing in a module

  USE mo_kind,          ONLY: dp
  USE mo_transpose,     ONLY: indx
  USE mo_decomposition, ONLY: gc=>global_decomposition,  &
                              lc => local_decomposition, &
                              debug_seriell
  USE mo_mpi,           ONLY: p_probe, p_recv, p_send, p_isend, p_wait, &
                              p_barrier, p_communicator_d
  USE mo_global_op,     ONLY: sum_zonal_sl
  USE mo_memory_base,   ONLY: set_stream_element_info
  USE mo_memory_gl,     ONLY: gl, lammp, phimp, sigmp
  USE mo_io_units,        ONLY: nout, nerr
  USE mo_exception,     ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_semi_lagrangian
  PUBLIC :: setup_semi_lagrangian
  PUBLIC :: semi_lagrangian_transport
  PUBLIC :: extend_semi_lagrangian
  PUBLIC :: cleanup_semi_lagrangian

  PUBLIC :: zonal_mass_integral
  PUBLIC :: global_mass_integral
  PUBLIC :: mass_fixer

  ! basic grid point resolution parameters
  !
  !  INTEGER, PARAMETER :: nxpt=1   ! no.of points outside active domain for
  !                                 ! interpolant
  !  INTEGER, PARAMETER :: jintmx=1 ! number of extra latitudes in polar region
  !
  ! choose the max. value for nxpt depending on resolution
  !
  ! replace by communication on demand!
  ! recommended: t21      : 18
  !              t106     : 60
  !              nprocb=1 : 2       ! 1 PE in E-W direction
  !
  INTEGER, PARAMETER :: nxpt=60
  ! INTEGER, PARAMETER :: nxpt=18
  ! INTEGER, PARAMETER :: nxpt=2

  ! Minimum value for joverlap is 2
  INTEGER, PARAMETER :: joverlap=3
  INTEGER, PARAMETER :: jintmx=joverlap-nxpt


  INTEGER, PARAMETER :: pmap = 20000 ! Dimension of artificial evenly spaced 
                                     ! vertical grid arrays

  !  integer pointers to physical data starting locations in 3-d arrays

  INTEGER, PARAMETER :: i1  = 1 + nxpt           ! model starting longitude
  INTEGER, PARAMETER :: j1  = 1 + nxpt + jintmx  ! model starting latitude

  INTEGER :: plon    ! number of local longitudes
  INTEGER :: plev    ! number of vertical levels
  INTEGER :: plat    ! number of local latitudes
  INTEGER :: plato2  ! number of local latitudes/2
  INTEGER :: pcnst   ! number of constituents (including water vapor)
  INTEGER :: plevp1  ! plev + 1
  INTEGER :: plond   ! local slt extended domain longitude
  INTEGER :: platd   ! local slt extended domain latitude per hemisphere
  INTEGER :: pgls    ! length of local latitude slice
  INTEGER :: jfirst  ! first index to be computed
  INTEGER :: jlast   ! last  index to be computed
  INTEGER :: i2pi    ! start of eastern long. extension
  INTEGER :: plono2  ! plon/2
  INTEGER :: istart  ! index to start computation
  INTEGER :: istop   ! index to stop  computation
  INTEGER :: js      ! index of southernmost model lat
  INTEGER :: jn      ! index of northernmost model lat
  INTEGER :: jstart  ! index for first model lat.
  INTEGER :: pbpts   ! (length of local latitude slice)*fields
  INTEGER :: plevm1

  REAL(dp) :: dphibr     ! reciprocal of maximum del phi

  REAL(dp), ALLOCATABLE :: sinlam(:)       ! Sine of longitudes in global grid 
                                           ! (no extension points).
  REAL(dp), ALLOCATABLE :: coslam(:)       ! Cosine of longitudes in global 
                                           ! grid (no extension points).
  REAL(dp), ALLOCATABLE :: detai(:)        ! Increment between model interfaces
                                           ! ("half" levels).
  REAL(dp), ALLOCATABLE :: detam(:)        ! Increment between model mid-levels
                                           ! ("full" levels)
  REAL(dp), ALLOCATABLE :: dphi(:,:)       ! Interval between latitudes in the
                                           ! extended grid
  REAL(dp)              :: dlam            ! Interval between longitudes
  REAL(dp), ALLOCATABLE :: etaint(:)       ! Half-index hybrid-levels from 
                                           ! sig(1/2) = etaint(1) = 0. to
                                           ! sig(plev+1/2) = 
                                           ! etaint(plevp1) = 1.
  REAL(dp), ALLOCATABLE :: etamid(:)       ! Full-index hybrid-levels in 
                                           ! vertical grid.
  REAL(dp), ALLOCATABLE :: gw(:)           ! Gauss weights for latitudes in 
                                           ! the global grid. 
                                           ! (These sum to 2.0.)
  REAL(dp), ALLOCATABLE :: phi(:,:)        ! Latitude values in the extended 
                                           ! grid.
  REAL(dp), ALLOCATABLE :: lam(:)          ! Longitude values in the extended 
                                           ! grid.
  REAL(dp), ALLOCATABLE :: lbasdy(:,:,:,:) ! Weights for Lagrange cubic 
                                           ! derivative estimates on the
                                           ! unequally spaced latitude grid
  REAL(dp), ALLOCATABLE :: lbasiy(:,:,:,:) ! Weights for Lagrange cubic 
                                           ! interpolation on the
                                           ! unequally spaced latitude grid
  REAL(dp), ALLOCATABLE :: lbasdz(:,:,:)   ! Weights for Lagrange cubic 
                                           ! derivative estimates on the
                                           ! unequally spaced vertical grid 
                                           ! (corresponding to model full 
                                           ! levels).
  REAL(dp), ALLOCATABLE :: lbassd(:,:,:)   ! Weights for Lagrange cubic 
                                           ! derivative estimates on the
                                           ! unequally spaced vertical grid 
                                           ! (corresponding to model half 
                                           ! levels).
  REAL(dp)              :: cwava           ! 1.0/(g*plon)

  ! Memory space for the semi Lagrangian scheme

  REAL(dp), ALLOCATABLE :: ub(:,:,:,:)
  REAL(dp), ALLOCATABLE :: vb(:,:,:,:)
  REAL(dp), ALLOCATABLE :: fb(:,:,:,:,:)

  INTEGER, ALLOCATABLE  :: kftype(:)

  INTEGER :: kdpmpf(pmap), kdpmph(pmap) 

  !  parameters common to many slt routines

  INTEGER, PARAMETER :: ppdy   = 4      ! length of interpolation grid stencil
  LOGICAL, PARAMETER :: plimdr = .TRUE. ! flag to limit derivatives
  
  ! Syntax for variables: [n][pe]{snd,rcv}_{eq,po,ea,we}
  ! n indicates "number of ..."
  ! pe indicates processor
  ! snd = send, rcv = receive, eq = equator, po = pole, ea = east, we = west

  INTEGER, ALLOCATABLE :: nxpt_a(:,:,:)  ! nxpt_a(plev,platd,2); array 

  REAL(dp), ALLOCATABLE :: pdel(:,:,:)

  REAL(dp), ALLOCATABLE :: qxl(:,:,:,:,:)
  REAL(dp), ALLOCATABLE :: qxr(:,:,:,:,:)
  REAL(dp), ALLOCATABLE :: uxl(:,:,:,:,:)
  REAL(dp), ALLOCATABLE :: uxr(:,:,:,:,:)

  ! send/receive "towards" equator
  INTEGER, ALLOCATABLE :: pe_snd_eq(:)   ! send lats to PEs with these id's
  INTEGER, ALLOCATABLE :: nsnd_eq(:)     ! # of lats to send to PEs
  INTEGER, ALLOCATABLE :: snd_eq(:,:,:)  ! my local lats to send
  INTEGER              :: npe_snd_eq     ! # of PEs to send to

  INTEGER, ALLOCATABLE :: nrcv_eq(:)     ! # of lats to receive from PEs
  INTEGER, ALLOCATABLE :: rcv_eq(:,:,:)  ! my local lats to receive
  INTEGER              :: npe_rcv_eq     ! # of PEs to receive from

  ! send/receive "towards" poles
  INTEGER, ALLOCATABLE :: pe_snd_po(:)
  INTEGER, ALLOCATABLE :: nsnd_po(:)
  INTEGER, ALLOCATABLE :: snd_po(:,:,:)
  INTEGER              :: npe_snd_po

  INTEGER, ALLOCATABLE :: nh(:,:,:)

  INTEGER, ALLOCATABLE :: nrcv_po(:)
  INTEGER, ALLOCATABLE :: rcv_po(:,:,:)
  INTEGER              :: npe_rcv_po
  INTEGER              :: nrcv_po_own
  INTEGER, ALLOCATABLE :: rcv_po_own(:,:)

  ! send to the east/receive from the west
  INTEGER, ALLOCATABLE :: blks_ea(:)

  ! send to the west/receive from the east
  INTEGER, ALLOCATABLE :: blks_we(:)

CONTAINS

!!============================================================================
!!
!! Section 1: init_semi_lagrangian is public interface 
!!
!! Contains: SUBROUTINE init_semi_lagrangian
!!           SUBROUTINE grdini
!!           SUBROUTINE grdxy
!!           SUBROUTINE basdy
!!           SUBROUTINE basdz
!!           SUBROUTINE basiy
!!           SUBROUTINE vrtmap
!!           SUBROUTINE lcdbas
!!           SUBROUTINE lcbas
!!
!!           SUBROUTINE init_mass_fixer is referenced in section mass_fixer
!!           SUBROUTINE setup_overlay is referenced in section overlay
!!

  SUBROUTINE init_semi_lagrangian

    ! Description:
    !
    ! Initializes slt-scheme
    !
    ! Authors:
    !
    ! M. Esch, MPI, March 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! T. Diehl, DKRZ, July 1999, parallel version
    ! 
    ! for more details see file AUTHORS

    USE mo_advection,    ONLY: nalatd, nalond, nalev, nacnst, nalat, nalon
    USE mo_physical_constants,    ONLY: grav
    USE mo_advection,    ONLY: jps
    USE mo_tracdef,      ONLY: trlist
    USE mo_hyb,          ONLY: ceta, cetah

    integer :: ih, ilat, k, i

    !  Executable statements 

    ! Set values

    plon   = lc%nglon
    plev   = lc%nlev
    plevp1 = lc%nlev + 1
    plat   = lc%nglat
    plato2 = lc%nglat/2
    
    pcnst  = jps + trlist% nadvec
    
    plond  = plon + 1 + 2*nxpt
    ! platd calculated per hemisphere
    platd  = plato2 + 2*nxpt + 2*jintmx

    pgls   = plon*plev
    !  jfirst = nxpt + 1
    !  jlast  = platd - nxpt - 1
    ! assume that the interpolant requires only one extra point
    jfirst = 2
    jlast  = platd - 2
    
    i2pi   = nxpt + plon + 1
    plono2 = plon/2
    istart = nxpt + 1
    istop  = nxpt + plon
    js     = 1 + nxpt + jintmx
    jn     = plato2 + nxpt + jintmx
    jstart = nxpt + jintmx + 1
    pbpts  = plond*plev*pcnst
    plevm1 = plev - 1

    nalev  = plev
    nacnst = pcnst
    nalat  = plat
    nalon  = plon
    nalatd = platd
    nalond = plond

    ! Arrays depending on platd get an extra dimension for the hemisphere 
    ! (since platd is the local number of latitudes on the extended grid 
    ! per hemisphere)

    ALLOCATE (coslam(plon), detai(plevp1), detam(plev), dphi(platd,2), &
              etaint(plevp1), etamid(plev), gw(plat), lam(plond),      &
              lbasdy(4,2,platd,2), lbasdz(4,2,plev),                   &
              lbasiy(4,2,platd,2), lbassd(4,2,plevp1), phi(platd,2),   &
              sinlam(plon))

    lbasdy(:,:,:,:) = 0.0_dp
    lbasiy(:,:,:,:) = 0.0_dp
    lbasdz(:,:,:)   = 0.0_dp
    lbassd(:,:,:)   = 0.0_dp

    ALLOCATE (kftype(pcnst))

    ! Memory space for the semi Lagrangian scheme

    ALLOCATE (ub(plond,plev,platd,2),       &
              vb(plond,plev,platd,2),       &
              fb(plond,plev,pcnst,platd,2))

    ALLOCATE (pdel(plond,plev,plat))

    ALLOCATE (qxl(plond,plev,pcnst,platd,2), &
              qxr(plond,plev,pcnst,platd,2), &
              uxl(plond,plev,2,platd,2),     &
              uxr(plond,plev,2,platd,2))    

    ub(:,:,:,:)   = 0.0_dp
    vb(:,:,:,:)   = 0.0_dp
    fb(:,:,:,:,:) = 0.0_dp

    ! Set mass fixer parameters

    CALL init_mass_fixer

    ! Define eta coordinates: Used for calculation etadot vertical velocity
    ! for slt.

    etamid(:) = ceta(:)
    etaint(1) = 0.0001_dp
    IF (plev > 19) etaint(1) = MIN(1.E-5_dp, etamid(1)*0.1_dp)
    etaint(2:) = cetah(2:)

    CALL grdini(grav)

    ! Initial guess for trajectory midpoints in spherical coords.
    ! lstart = T:  Use arrival points as initial guess for trajectory midpoin
    ! lstart = F:  Use calculated trajectory midpoints from previous time
    ! step as first guess.
    ! NOTE:  Reduce number of iters necessary for convergence after nstep=1

    CALL set_stream_element_info (gl,'lammp',lrerun=.TRUE.)
    CALL set_stream_element_info (gl,'phimp',lrerun=.TRUE.)
    CALL set_stream_element_info (gl,'sigmp',lrerun=.TRUE.)

    ! in case this fields are contained in the restart file the following
    ! set values are overwritten. In case not, they make life more expensive.
    ! Changing the advection scheme during a run disables reproducibility.    

    DO ih = 1,2
      DO ilat = 1, plato2
        DO k = 1, plev
          DO i = 1, plon
            lammp(i,k,(ih-1)*plato2 + ilat) = lam(i1 + i   - 1)
            phimp(i,k,(ih-1)*plato2 + ilat) = phi(j1 + ilat - 1,ih)
            sigmp(i,k,(ih-1)*plato2 + ilat) = etamid(k)
          END DO
        END DO
      END DO
    END DO
    
    ALLOCATE (nxpt_a(plev,platd,2))

    CALL setup_overlap

  END SUBROUTINE init_semi_lagrangian
!=============================================================================
  SUBROUTINE grdini(gravit)

    !
    ! Description:
    !
    ! Initialize model and extended grid parameters.
    ! Initialize weights for Lagrange cubic derivative estimates.
    ! Initialize weights for Lagrange cubic interpolant.
    !
    ! Method:
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    ! Parallel version:  T. Diehl, DKRZ, July 1999
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: gravit

    !  Local scalars:
    INTEGER :: j, k, ih

    !  Local arrays:
    REAL(dp) :: detailn(plevp1)    ! dlog(etaint)
    REAL(dp) :: detamln(plev)     ! dlog(etamid)
    REAL(dp) :: etailn(plevp1)     ! log(etaint)
    REAL(dp) :: etamln(plev)      ! log(etamid)

    !  Intrinsic functions
    INTRINSIC LOG

    !  Executable statements

    ! Initialize extended horizontal grid coordinates.

    CALL grdxy

    ! Basis functions for computing Lagrangian cubic derivatives
    ! on unequally spaced latitude and vertical grids.

    CALL basdy(phi(1,1),lbasdy(1,1,1,1))
    CALL basdy(phi(1,2),lbasdy(1,1,1,2))
    CALL basdz(plev,etamid,lbasdz)
    CALL basdz(plevp1,etaint,lbassd)

    ! Basis functions for computing weights for Lagrangian cubic
    ! interpolation on unequally spaced latitude grids.

    CALL basiy(phi(1,1),lbasiy(1,1,1,1))
    CALL basiy(phi(1,2),lbasiy(1,1,1,2))

    ! Compute interval lengths in latitudinal grid

    dphi(:,:)   = 0.0_dp
    dphibr      = 0.0_dp
    DO ih = 1,2
      DO j = 1, platd - 1
        dphi(j,ih) = phi(j+1,ih) - phi(j,ih)
        dphibr = MAX(ABS(dphi(j,ih)),dphibr)
      END DO
    END DO
    dphibr = 1.0_dp/dphibr

    ! Compute interval lengths in vertical grids.

detam(:) = 0.0_dp
detai(:) = 0.0_dp

    DO k = 1, plev
      etamln(k) = LOG(etamid(k))
    END DO
    DO k = 1, plevp1
      etailn(k) = LOG(etaint(k))
    END DO
    DO k = 1, plev - 1
      detam(k) = etamid(k+1) - etamid(k)
      detamln(k) = etamln(k+1) - etamln(k)
    END DO
    DO k = 1, plev
      detai(k) = etaint(k+1) - etaint(k)
      detailn(k) = etailn(k+1) - etailn(k)
    END DO

    ! Build artificial evenly spaced vertical grid for use in determining
    ! vertical position of departure point.
    ! Build one grid for full model levels and one for half levels.

    CALL vrtmap(plev,etamln,detamln,kdpmpf)
    CALL vrtmap(plevp1,etailn,detailn,kdpmph)

    ! Compute moisture integration constant

    cwava = 1.0_dp/(lc%nlon*gravit)

  END SUBROUTINE grdini
!=============================================================================
  SUBROUTINE grdxy

    ! Description:
    !
    ! Define the "extended" grid used in the semi-Lagrangian transport
    ! scheme.  The longitudes are equally spaced and the latitudes are
    ! Gaussian.  The global grid is extended to include "wraparound" points
    ! on all sides.
    !
    ! Method:
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    ! parallel version:  T. Diehl, DKRZ, July 1999
    !
    ! for more details see file AUTHORS
    !

    USE mo_gaussgrid,     ONLY: gauaw ! module subroutine
    USE mo_control,       ONLY: ngl

    ! dlam    Length of increment in longitude grid.
    ! lam     Longitude values in the extended grid.
    ! phi     Latitude values in the extended grid.
    ! gw      Gauss weights for latitudes in the global grid.
    !         (These sum to 2.0 like the ones in CCM1.)
    ! sinlam  Sine of longitudes in global grid (no extension points).
    ! coslam  Cosine of longitudes in global grid (no extension points).

    !  Local scalars:
    REAL(dp)    :: lam0, pi
    INTEGER :: i, ig, j
    INTEGER :: jg, jgn
    INTEGER :: jstartg, jstopg

    ! Local global arrays
    !  INTEGER, PARAMETER :: pglatd=lc%nlat + 2*jintmx + 2*nxpt
    !  REAL(dp) :: phig(pglatd), wrkg(pglatd), wg(lc%nlat)
    INTEGER :: pglatd
    REAL(dp) :: phig(ngl+2*jintmx+2*nxpt)
    REAL(dp) :: wrkg(ngl+2*jintmx+2*nxpt)
    REAL(dp) :: wg(ngl+2*jintmx+2*nxpt)

    !  Intrinsic functions
    INTRINSIC ASIN, ATAN, COS, SIN

    !  Executable statements

    pglatd = ngl + 2*jintmx + 2*nxpt
    lam0 = 0.0_dp
    pi = 4.0_dp*ATAN(1.0_dp)
    jstartg = jstart
    jstopg = jstart - 1 + lc%nlat

    ! Interval length in equally spaced longitude grid.

    dlam = 2.0_dp*pi/REAL(lc%nlon,dp)

    ! Longitude values on extended grid.

    !  ! compute offset for general case
    !  offset = 0
    !  ppid_current   = lc%pid
    !  DO ioff=1,lc%pid_ln-1
    !     plft_current         = gc(ppid_current)%pid_lft
    !     offset               = offset + gc(plft_current)%nglon
    !     ppid_current         = plft_current
    !  END DO

    DO i = 1, plond
      lam(i) = REAL(i - istart + lc%glons(1)-1,dp)*dlam + lam0
      !     lam(i) = REAL(dp)(i-istart + offset)*dlam + lam0
    END DO

    ! Compute Gauss latitudes and weights.  On return; phi contains the
    ! sine of the latitudes starting closest to the north pole and going
    ! toward the south; w contains the corresponding Gauss weights.

    CALL gauaw(phig,wg,lc%nlat)

    ! Reorder and compute latitude values.

    DO j = jstartg, jstopg
      wrkg(j) = ASIN(phig(jstopg-j+1))
    END DO
    phig(jstartg:jstopg)=wrkg(jstartg:jstopg)

    ! North and south poles.

    phig(jstartg-1) = -pi*0.5_dp
    phig(jstopg+1) = pi*0.5_dp

    ! Extend Gauss latitudes below south pole so that the spacing above
    ! the pole is symmetric, and phi is decreasing, i.e., phi < -pi/2

    IF (jstartg > 2) THEN
      DO j = 1, jstartg - 2
        phig(j) = -pi - phig(2*jstartg-2-j)
      END DO
    END IF

    ! Analogously for Northern Hemisphere

    IF (pglatd > jstopg+1) THEN
      DO j = jstopg + 2, pglatd
        phig(j) = pi - phig(2*jstopg+2-j)
      END DO
    END IF

    ! Pack into local phi and w arrays

    DO j=1,platd
      ! Southern Hemisphere
      jg = lc%glats(1) + j -1
      phi(j,1) = phig(jg)
      ! Northern hemisphere
      jgn = pglatd - lc%glats(1) - plato2 - 2*nxpt - 2*jintmx + j + 1
      phi(j,2) = phig(jgn)
    END DO
    DO j=1,plato2
      ! Southern Hemisphere
      jg = lc%glats(1) + j - 1
      gw(j) = wg(jg)
      ! Northern Hemisphere
      jgn = lc%nlat - lc%glats(1) - plato2 + j + 1
      gw(j+plato2) = wg(jgn)
    END DO

    ! Sine and cosine of longitude.

    ig = 0
    DO i = istart, istop
      ig = ig + 1
      sinlam(ig) = SIN(lam(i))
      coslam(ig) = COS(lam(i))
    END DO

  END SUBROUTINE grdxy
!=============================================================================
  SUBROUTINE basdy(pphi,pbasdy)

    ! Description:
    !
    ! Compute weights for the calculation of derivative estimates
    !
    ! Method:
    !
    ! Compute weights for the calculation of derivative estimates at the two
    ! center points of the four point stencil for each interval in the
    ! unequally spaced latitude grid. Estimates are from differentiating
    ! a Lagrange cubic polynomial through the four point stencil.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: pphi(platd)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: pbasdy(4,2,platd)

    ! pphi     Latitude values in the extended grid.
    ! plbasdy  Weights for derivative estimates based on Lagrange cubic
    !          polynomial on the unequally spaced latitude grid.
    !          If grid interval j (in extended grid) is surrounded by
    !          a 4 point stencil, then the derivative at the "bottom"
    !          of the interval uses the weights lbasdy(1,1,j),
    !          lbasdy(2,1,j), lbasdy(3,1,j), and lbasdy(4,1,j).
    !          The derivative at the "top" of the interval
    !          uses lbasdy(1,2,j), lbasdy(2,2,j), lbasdy(3,2,j),
    !          and lbasdy(4,2,j).

    !  Local scalars:
    INTEGER :: jj

    !  Executable statements

    DO jj = jfirst, jlast
      CALL lcdbas(pphi(jj-1),pbasdy(1,1,jj),pbasdy(1,2,jj))
    END DO

  END SUBROUTINE basdy
!=============================================================================
  SUBROUTINE basdz(pkdim,psig,pbasdz)

    ! Description:
    !
    ! Compute weights for the calculation of derivative estimates
    !
    ! Method:
    !
    ! Compute weights for the calculation of derivative estimates at two
    ! center points of the four point stencil for each interval in the
    ! unequally spaced vertical grid (as defined by the array sig).
    ! Estimates are from differentiating a Lagrange cubic polynomial
    ! through the four point stencil.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pkdim

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: psig(pkdim)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: pbasdz(4,2,pkdim)

    ! pkdim   Number of grid points in vertical grid.
    ! psig    Sigma values in the vertical grid.
    ! pbasdz  Weights for derivative estimates based on Lagrange cubic
    !         polynomial on the unequally spaced vertical grid.
    !         If grid interval j is surrounded by a 4 point stencil,
    !         then the derivative at the "top" of the interval (smaller
    !         sigma value) uses the weights lbasdz(1,1,j),lbasdz(2,1,j),
    !         lbasdz(3,1,j), and lbasdz(4,1,j).  The derivative at the
    !         "bottom" of the interval uses lbasdz(1,2,j), lbasdz(2,2,j),
    !         lbasdz(3,2,j), and lbasdz(4,2,j).  (Recall the vertical
    !         level indices increase from the top of the atmosphere
    !         towards the bottom.)

    !  Local scalars:
    INTEGER :: kk

    !  Executable statements

    DO kk = 2, pkdim - 2
      CALL lcdbas(psig(kk-1),pbasdz(:,1,kk),pbasdz(:,2,kk))
    END DO

  END SUBROUTINE basdz
!=============================================================================
  SUBROUTINE basiy(pphi,pbasiy)

    ! Description:
    !
    ! Compute weights used in Lagrange cubic polynomial interpolation
    !
    ! Method:
    !
    ! Compute weights used in Lagrange cubic polynomial interpolation in
    ! the central interval of a four point stencil.  Done for each interval
    ! in the unequally spaced latitude grid.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: pphi(platd)  ! grid values in extended gri

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (out) :: pbasiy(4,2,platd)  ! Weights for Lagrange cubic intp.

    !  Local scalars:
    INTEGER :: jj

    !  Executable statements

    DO jj = jfirst, jlast
      CALL lcbas(pphi(jj-1),pbasiy(1,1,jj),pbasiy(1,2,jj))
    END DO

  END SUBROUTINE basiy
!=============================================================================
  SUBROUTINE vrtmap(pkdim,sigln,dsigln,kdpmap)
    
    ! Description:
    !
    ! Map indices of an artificial evenly spaced (in log) vertical grid to
    ! the indices of the log of the model vertical grid
    !
    ! Method:
    !
    ! The resultant array of mapped indices will be used by "kdpfnd"
    ! to find the vertical location of any departure point relative
    ! to the model grid.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !
    
    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pkdim
    
    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: dsigln(pkdim)    ! intervals between model levels (log)
    REAL(dp), INTENT (in) :: sigln(pkdim)     ! model levels (log(eta))
    
    !  Array arguments with intent(Out):
    INTEGER, INTENT (out) :: kdpmap(pmap) ! array of mapped indices
    
    !  Local scalars:
    REAL(dp) :: del                ! artificial grid interval
    REAL(dp) :: dpa                ! artificial departure point
    REAL(dp) :: eps                ! epsilon factor
    INTEGER :: k, kk
    INTEGER :: newmap         ! estimated value of "pmap"
    
    !  Local arrays:
    INTEGER :: imin(1)
    
    !!LK  !  External functions
    !!LK  INTEGER, EXTERNAL :: ismin
    
    !  Executable statements
    
    eps = 1.0E-05_dp
    del = (sigln(pkdim)-sigln(1))/REAL(pmap,dp)
    !!LK  imin = ismin(pkdim-1,dsigln,1)
    imin(:) = MINLOC(dsigln(1:pkdim-1))
    IF (del+eps>=dsigln(imin(1))) THEN
      newmap = INT((sigln(pkdim)-sigln(1))/dsigln(imin(1)) + 1)
      WRITE(nout,'(A)') ' VRTMAP:  Not enough artificial grid intervals.'
      WRITE(nout,'(A,I20)') ' Currently, "pmap" is set to ', pmap
      WRITE(nout,'(A,I20)') ' Reset parameter "pmap" to at least ', newmap
      CALL finish('vrtmap','Run terminated.')
    END IF
    
    kdpmap(1) = 1
    !CDIR LOOPCNT=20000
    DO kk = 2, pmap
      dpa = sigln(1) + REAL(kk-1,dp)*del
      DO k = 1, pkdim - 1
        IF (dpa>sigln(k)+eps) THEN
          kdpmap(kk) = k
        END IF
      END DO
    END DO
    
  END SUBROUTINE vrtmap
!=============================================================================

  SUBROUTINE lcbas(grd,bas1,bas2)

    ! Description:
    !
    ! Evaluate the partial Lagrangian cubic basis functions (denominator
    ! only ) for the grid points and gather grid values
    !
    ! Authors:
    !
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) ::  grd(4)   ! Grid stencil

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: bas1(4)  ! Grid values on stencil
    REAL(dp), INTENT (out) :: bas2(4)  ! Lagrangian basis functions

    !  Local scalars:
    !  Grid value differences used in weights
    REAL(dp) :: x0mx1, x0mx2, x0mx3, x1mx2, x1mx3, x2mx3

    !  Executable statements

    x0mx1 = grd(1) - grd(2)
    x0mx2 = grd(1) - grd(3)
    x0mx3 = grd(1) - grd(4)
    x1mx2 = grd(2) - grd(3)
    x1mx3 = grd(2) - grd(4)
    x2mx3 = grd(3) - grd(4)

    bas1(1) = grd(1)
    bas1(2) = grd(2)
    bas1(3) = grd(3)
    bas1(4) = grd(4)

    bas2(1) =  1.0_dp/(x0mx1*x0mx2*x0mx3)
    bas2(2) = -1.0_dp/(x0mx1*x1mx2*x1mx3)
    bas2(3) =  1.0_dp/(x0mx2*x1mx2*x2mx3)
    bas2(4) = -1.0_dp/(x0mx3*x1mx3*x2mx3)

  END SUBROUTINE lcbas

  SUBROUTINE lcdbas(grd,dbas2,dbas3)

    ! Description:
    !
    ! Calculate weights.
    !
    ! Method:
    !
    ! Calculate weights used to evaluate derivative estimates at the
    ! inner grid points of a four point stencil based on Lagrange
    ! cubic polynomial through four unequally spaced points.
    !
    ! Authors:
    !
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) ::  grd(4)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: dbas2(4), dbas3(4)

    ! grd    Coordinate values of four points in stencil.
    ! dbas2  Derivatives of the four basis functions at grid point 2.
    ! dbas3  Derivatives of the four basis functions at grid point 3.

    !  Local scalars:
    REAL(dp) :: x1, x2, x3, x4                        ! - grid values
    REAL(dp) :: x1mx2,x1mx3,x1mx4,x2mx3,x2mx4,x3mx4   ! - differences of grid values

    !  Executable atatements

    x1 = grd(1)
    x2 = grd(2)
    x3 = grd(3)
    x4 = grd(4)

    x1mx2 = x1 - x2
    x1mx3 = x1 - x3
    x1mx4 = x1 - x4
    x2mx3 = x2 - x3
    x2mx4 = x2 - x4
    x3mx4 = x3 - x4

    dbas2(1) = x2mx3*x2mx4/(x1mx2*x1mx3*x1mx4)
    dbas2(2) = -1.0_dp/x1mx2 + 1.0_dp/x2mx3 + 1.0_dp/x2mx4
    dbas2(3) = -x1mx2*x2mx4/(x1mx3*x2mx3*x3mx4)
    dbas2(4) = x1mx2*x2mx3/(x1mx4*x2mx4*x3mx4)

    dbas3(1) = -x2mx3*x3mx4/(x1mx2*x1mx3*x1mx4)
    dbas3(2) = x1mx3*x3mx4/(x1mx2*x2mx3*x2mx4)
    dbas3(3) = -1.0_dp/x1mx3 - 1.0_dp/x2mx3 + 1.0_dp/x3mx4
    dbas3(4) = -x1mx3*x2mx3/(x1mx4*x2mx4*x3mx4)

  END SUBROUTINE lcdbas
!!
!!============================================================================

!!============================================================================
!!
!! Section 2: setup_semi_lagrangian is public interface 
!!
!! Contains: SUBROUTINE setup_semi_lagrangian

!OCL NOALIAS

  SUBROUTINE setup_semi_lagrangian

    ! Description:
    !
    ! Fills 3dim arrays for the semi lagrangian scheme
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, February 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! T. Diehl, DKRZ, July 1999, parallel version
    ! 
    ! for more details see file AUTHORS
    !
    
    USE mo_gaussgrid,     ONLY: gl_sqcst
    USE mo_tracdef,       ONLY: trlist
    USE mo_scan_buffer,   ONLY: u, v
    USE mo_memory_g1a,    ONLY: qm1, xlm1, xim1, xtm1, alpsm1

    !  Local scalars: 
    REAL(dp) :: zlimit, zrcst
    INTEGER :: i, ita, jk, jl, jt, k, l, jglat, ih
    
    !  Local arrays: 
    REAL(dp) :: zphm1(plon,plevp1,plat), zpsm1(plon,plat)
    
    !  External subroutines 
    EXTERNAL pres
    
    !  Intrinsic functions 
    INTRINSIC ABS
    
    !  Executable statements 
    
    ! Set scalar fields to zero if the absolute value
    ! is smaller than a hardware dependent value
    ! to prevent slt scheme from underflow.
    
    zlimit = 1.0E-200_dp
    
    WHERE (ABS(qm1(:,:,:)) < zlimit) 
      qm1(:,:,:) = 0.0_dp
    END WHERE
    
    WHERE (ABS(xlm1(:,:,:)) < zlimit)
      xlm1(:,:,:) = 0.0_dp
    END WHERE
    
    WHERE (ABS(xim1(:,:,:)) < zlimit)
      xim1(:,:,:) = 0.0_dp
    END WHERE

    zpsm1(:,:) = EXP(alpsm1(1:plon,:))
    
    ! Compute half level pressure at t-dt for pdel.
    
    DO i = 1, plat
      CALL pres(zphm1(:,:,i), plon, zpsm1(:,i), plon)
    END DO
    
    ! Fill up 3d arrays for slt
    ! and divide winds by cos(latitude)
    
    DO i = 1, plat
      IF (i <= plato2) THEN
        ih = 1
        k = i + js - 1
      ELSE
        ih = 2
        k = i + js - 1 - plato2
      END IF
      
      !     l = MOD((i+1),2)*(plat-(i-1)/2)+MOD(i,2)*(i+1)/2
      l = plat-i+1
      
      jglat  = lc%glat(l)   ! global index (continuous from north to south)
      zrcst  = 1.0_dp/gl_sqcst(jglat)                 ! 1./cos(latitude)
      
      ! arrays for SLT must have the format south -> north
      DO jk = 1, plev
        DO jl = 1, plon
          pdel(jl+nxpt,jk,i)    = zphm1(jl,jk+1,l)-zphm1(jl,jk,l)
          ub(jl+nxpt,jk,k,ih)   = u(jl,jk,l)*zrcst
          vb(jl+nxpt,jk,k,ih)   = v(jl,jk,l)*zrcst
          fb(jl+nxpt,jk,1,k,ih) = qm1(jl,jk,l)
          fb(jl+nxpt,jk,2,k,ih) = xlm1(jl,jk,l)
          fb(jl+nxpt,jk,3,k,ih) = xim1(jl,jk,l)
        END DO
      END DO
      
      ita = 0
      DO jt = 1, trlist% ntrac
        IF (trlist% ti(jt)% ntran /= 0 ) THEN
          ita = ita + 1
          DO jk = 1, plev
            DO jl = 1, plon
              IF (ABS(xtm1(jl,jk,jt,l))<zlimit) xtm1(jl,jk,jt,l) = 0.0_dp
              fb(jl+nxpt,jk,3+ita,k,ih) = xtm1(jl,jk,jt,l)
            END DO
          END DO
        END IF
      END DO
    END DO
    
  END SUBROUTINE setup_semi_lagrangian

!!
!!============================================================================

!!============================================================================
!!
!! Section 2: semi_lagrangian_transport is public interface 
!!
!! Contains: SUBROUTINE semi_lagrangian_transport


  SUBROUTINE cleanup_semi_lagrangian

    DEALLOCATE (blks_we)
    DEALLOCATE (blks_ea)
    DEALLOCATE (rcv_po_own)
    DEALLOCATE (rcv_po)
    DEALLOCATE (nrcv_po)
    DEALLOCATE (nh)
    DEALLOCATE (snd_po)
    DEALLOCATE (nsnd_po)
    DEALLOCATE (pe_snd_po)
    DEALLOCATE (rcv_eq)
    DEALLOCATE (nrcv_eq)
    DEALLOCATE (snd_eq)
    DEALLOCATE (nsnd_eq)
    DEALLOCATE (pe_snd_eq)

    DEALLOCATE (nxpt_a)
    DEALLOCATE (uxr, uxl, qxr, qxl)
    DEALLOCATE (pdel)
    DEALLOCATE (fb, vb, ub)
    DEALLOCATE (kftype)
    DEALLOCATE (sinlam, phi, lbassd, lbasiy, lbasdz, lbasdy, lam, gw, &
                etamid, etaint, dphi, detam, detai, coslam)

  END SUBROUTINE cleanup_semi_lagrangian

  SUBROUTINE bandij(dlam, phib, lamp, phip, iband, jband, nxpt_a)

    ! Description:
    !
    ! Calculate longitude and latitude indices
    !
    ! Method:
    !
    ! Calculate longitude and latitude indices that identify the
    ! intervals in the extended grid that contain the departure points.
    ! Upon entry, all dep. points should be within jintmx intervals of the
    ! Northern- and Southern-most model latitudes. Note: the algorithm
    ! relies on certain relationships of the intervals in the Gaussian grid.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    ! parallel version:  T. Diehl, July 1999
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    REAL(dp),    INTENT (in) :: dlam

    !  Array arguments with intent(In):
    REAL(dp),    INTENT (in) :: lamp(pgls), phib(platd), phip(pgls)
    INTEGER, INTENT (in) :: nxpt_a(plev)

    !  Array arguments with intent(Out):
    INTEGER, INTENT (out) :: iband(pgls), jband(pgls)

    ! dlam    Length of increment in equally spaced longitude grid (radians)
    ! phib    Latitude values for the extended grid.
    ! lamp    Longitude coordinates of the points.  It is assumed that
    !         0.0 .le. lamp(i) .lt. 2*pi .
    ! phip    Latitude coordinates of the points.
    ! iband   Longitude index of the points.  This index points into
    !         the extended arrays, e.g.,
    !         lam(iband(i)) .le. lamp(i) .lt. lam(iband(i)+1) .
    ! jband   Latitude index of the points.  This index points into
    !         the extended arrays, e.g.,
    !         phib(jband(i)) .le. phip(i) .lt. phib(jband(i)+1) .

    !  Local scalars:
    INTEGER :: i           ! index
    INTEGER :: k, ieloc(1)
    LOGICAL :: lerror(pgls)

    !  Intrinsic functions
    INTRINSIC INT, MERGE, ANY, MAXLOC

    !  Executable statements

    ! Longitude indices.

    DO i = 1, pgls

      !   iband(i) = istart + CEILING(lamp(i)/dlam) - lc%glons(1)
      iband(i) = istart + INT(lamp(i)/dlam) - lc%glons(1) &
           + MERGE(0,1,lamp(i) < 0.0_dp)

      ! For dynamic sltpads need to check the level maximum

      k = (i-1)/plon + 1  ! determine current level

      !   lerror(i) = (iband(i) <  istart-nxpt_a(k)) .OR. &
      !               (iband(i) >= plon+1+nxpt+nxpt_a(k))
      !
      ! the original check was not sufficient
      !   cf. subroutine cubxdr:
      !   arrays are indexed with i-1 to i+2
      !
      !   if only one PE is present in e-W direction the test is skipped.
      !   nxpt must be at least 2 in this case
      !
      lerror(i) = .false.
      IF (lc% nprocb > 1) THEN
        lerror(i) = (iband(i) <  istart-nxpt_a(k)+1) .OR. &
                    (iband(i) >  plon+nxpt+nxpt_a(k)-2)
      ENDIF

    END DO

    IF (ANY(lerror)) THEN
      ieloc = MAXLOC((/(i,i=1,pgls)/),lerror)
      WRITE (nerr,*) 'bandij: set B = ', lc%set_b, &
           ' set A = ', lc%set_a, &
           ' bandi = ', iband(ieloc(1))
      WRITE (nerr,*) ' Increase nxpt and/or isave (in sltb1)'
      CALL finish('bandij','Point out of bounds in E-W.')
    END IF

    ! Latitude indices.

    jband(:) = INT ((phip(:) - phib(1))*dphibr + 1._dp)

    DO i = 1, pgls
      IF (phip(i) >= phib(jband(i)+1)) jband(i) = jband(i)+1
    END DO

    !  lerror(:) = jband(:) < 1 .OR. jband(:) > platd
    !
    ! the original check was not sufficient
    !   cf. subroutine herxin:
    !   arrays are indexed with jband-1 to jband+2

    lerror(:) = jband(:) < 2 .OR. jband(:) > platd-2

    IF (ANY(lerror)) THEN
      ieloc = MAXLOC((/(i,i=1,pgls)/),lerror)
      WRITE (nerr,*) 'bandij: set B = ', lc%set_b, &
           ' set A = ', lc%set_a, &
           ' bandij = ',jband(ieloc(1))
      CALL finish('bandij','Point out of bounds in N-S.')
    END IF

  END SUBROUTINE bandij


  SUBROUTINE cubxdr(pidim,ibeg,len,dx,f,qxl,qxr)

    ! Description:
    !
    ! Compute Lagrangian cubic derivative estimates for data on an equally
    ! spaced grid.
    !
    ! Method:
    !
    ! Compute Lagrangian cubic derivative estimates for data on an equally
    ! spaced grid.  Suppose grid interval i is centered in a 4 point
    ! stencil consisting of grid points i-1, i, i+1, and i+2.  Then the
    ! derivative at the left edge of the interval (i.e., grid point i)
    ! is stored in qxl(i), and the derivative at the right edge of the
    ! interval (i.e., grid point i+1) is stored in qxr(i).  Note that
    ! qxl(i) is not necessarily equal to qxr(i-1) even though both of
    ! these values are estimates of the derivative at grid point i.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: dx
    INTEGER, INTENT (in) :: ibeg, len, pidim

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: f(pidim)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: qxl(pidim), qxr(pidim)

    ! pidim   Length of f, qxl, and qxr.
    ! ibeg    First interval of grid for which derivatives are computed.
    ! len     Number of grid intervals for which derivatives are computed.
    !         (There are pidim - 1 intervals between the pidim gridpoints
    !         represented in f, qxl, and qxr.)
    ! dx      Value of grid spacing.
    ! f       Values on equally spaced grid for which derivatives are
    !         computed.
    ! qxl     qxl(i) is the derivative at the left  edge of interval i.
    ! qxr     qxr(i) is the derivative at the right edge of interval i.

    !  Local scalars:
    REAL(dp) :: rdx6        ! normalization weight
    INTEGER :: i, iend

    !  Executable statements

    iend = ibeg + len - 1
    rdx6 = 1.0_dp/(6.0_dp*dx)

    DO i = ibeg, iend
      qxl(i) = (-2.0_dp*f(i-1)-3.0_dp*f(i)+6.0_dp*f(i+1)-f(i+2))*rdx6
      qxr(i) = (f(i-1)-6.0_dp*f(i)+3.0_dp*f(i+1)+2.0_dp*f(i+2))*rdx6
    END DO

  END SUBROUTINE cubxdr

  SUBROUTINE cubydr(pf,fint,wdy,jdp,fyb,fyt)

    ! Description:
    !
    ! Compute Lagrangian cubic derivative estimates.
    !
    ! Method:
    !
    ! Compute Lagrangian cubic derivative estimates at both ends of the
    ! intervals in the y coordinate (unequally spaced) containing the
    ! departure points for the latitude slice being forecasted.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! Modified:          U. Schlese,  DKRZ - Hamburg,  May 1994
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pf

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: fint(pgls,ppdy,pf), wdy(4,2,platd)
    INTEGER, INTENT (in) :: jdp(pgls)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: fyb(pgls,pf), fyt(pgls,pf)

    ! pf      Number of fields being interpolated.
    ! fint    (fint(i,k,j,m),j=1,ppdy) contains the x interpolants at each
    !         latitude needed for the y derivative estimates at the
    !         endpoints of the interval that contains the departure point
    !         for grid point (i,k).  The last index of fint allows for
    !         interpolation of multiple fields. fint is generated by a
    !         call to herxin.
    ! wdy     Weights for Lagrange cubic derivative estimates on the
    !         unequally spaced latitude grid. If grid interval j (in
    !         extended array) is surrounded by a 4 point stencil, then
    !         the derivative at the "bottom" of the interval uses the
    !         weights wdy(1,1,j),wdy(2,1,j), wdy(3,1,j), and wdy(4,1,j).
    !         The derivative at the "top" of the interval uses wdy(1,2,j),
    !         wdy(2,2,j), wdy(3,2,j), and wdy(4,2,j).
    ! jdp     jdp(i,k) is the index of the y-interval that contains the
    !         departure point corresponding to global grid point (i,k) in
    !         the latitude slice being forecasted.
    !         Suppose yb contains the y-coordinates of the extended array
    !         and ydp(i,k) is the y-coordinate of the departure point
    !         corresponding to grid point (i,k).  Then,
    !         yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
    ! fyb     fyb(i,k,.) is the derivative at the bottom of the y interval
    !         that contains the departure point of global grid point (i,k).
    ! fyt     fyt(i,k,.) is the derivative at the top of the y interval
    !         that contains the departure point of global grid point (i,k).

    !  Local scalars:
    INTEGER :: i, m

    !  Local arrays:
    REAL(dp) :: wdytem(pgls,4,2)  ! Work array to gather weights

    !  Executable statements

    ! Load temp arrays with weights for Lagrangian cubic derivative
    ! estimates on unequally spaced grid.

    DO i = 1, pgls
      wdytem(i,1,1) = wdy(1,1,jdp(i))
      wdytem(i,2,1) = wdy(2,1,jdp(i))
      wdytem(i,3,1) = wdy(3,1,jdp(i))
      wdytem(i,4,1) = wdy(4,1,jdp(i))
      wdytem(i,1,2) = wdy(1,2,jdp(i))
      wdytem(i,2,2) = wdy(2,2,jdp(i))
      wdytem(i,3,2) = wdy(3,2,jdp(i))
      wdytem(i,4,2) = wdy(4,2,jdp(i))
    END DO

    ! Loop over fields.

    DO m = 1, pf
      DO i = 1, pgls
        fyb(i,m) = wdytem(i,1,1)*fint(i,1,m) + wdytem(i,2,1)*fint(i,2,m) + &
             wdytem(i,3,1)*fint(i,3,m) + wdytem(i,4,1)*fint(i,4,m)

        fyt(i,m) = wdytem(i,1,2)*fint(i,1,m) + wdytem(i,2,2)*fint(i,2,m) + &
             wdytem(i,3,2)*fint(i,3,m) + wdytem(i,4,2)*fint(i,4,m)
      END DO
    END DO

  END SUBROUTINE cubydr

  SUBROUTINE cubzdr(pidim,pkdim,f,pbasdz,dfz1,dfz2)

    ! Description:
    !
    ! Vertical derivative estimates for a vertical slice using Lagrangian
    ! cubic formulas.
    !
    ! Method:
    !
    ! Vertical derivative estimates for a vertical slice using Lagrangian
    ! cubic formulas.  Derivatives are set to zero at the top and bottom.
    ! At the "inner nodes" of the top and bottom intervals, a "one sided"
    ! estimate is used.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pidim, pkdim

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: f(pidim,pkdim), pbasdz(4,2,pkdim)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: dfz1(pidim,pkdim), dfz2(pidim,pkdim)

    ! pidim   Horizontal dimension of arrays. (pidim = plon)
    ! pkdim   Vertical dimension of arrays.
    ! f       Vertical slice of data for which derivative estimates are made
    ! pbasdz  Lagrangian cubic basis functions for evaluating the
    !         derivatives on the unequally spaced vertical grid.
    ! dfz1    dfz1 contains derivative estimates at the "top" edges of the
    !         intervals in the f array.
    ! dfz2    dfz2 contains derivative estimates at the "bottom" edges of
    !         the intervals in the f array.

    !  Local scalars:
    INTEGER :: i, k   ! Indices

    !  Executable statements

dfz1 = 0.0_dp
dfz2 = 0.0_dp

    DO k = 2, pkdim - 2
      DO i = 1, pidim

        ! Lagrangian derivative estimates (cubic) for the two center nodes in a
        ! four node stencil.

        dfz1(i,k) = pbasdz(1,1,k)*f(i,k-1) + pbasdz(2,1,k)*f(i,k) + &
                    pbasdz(3,1,k)*f(i,k+1) + pbasdz(4,1,k)*f(i,k+2)

        dfz2(i,k) = pbasdz(1,2,k)*f(i,k-1) + pbasdz(2,2,k)*f(i,k) + &
                    pbasdz(3,2,k)*f(i,k+1) + pbasdz(4,2,k)*f(i,k+2)
      END DO
    END DO

    ! Constrain derivatives to zero at top and bottom of vertical grid.
    ! At the interior nodes of the intervals at the top and bottom of the
    ! vertical grid, use the derivative estimate at that same node for the
    ! adjacent interval.  (This is a "one-sided" estimate for that node.)

    DO i = 1, pidim
      dfz1(i,1) = 0.0_dp
      dfz2(i,1) = dfz1(i,2)
      dfz1(i,pkdim-1) = dfz2(i,pkdim-2)
      dfz2(i,pkdim-1) = 0.0_dp
    END DO

  END SUBROUTINE cubzdr

!OCL NOVREC

  SUBROUTINE extx(pkcnst,pkdim,fb)

    ! Description:
    !
    ! Copy data to the longitude extensions of the extended array
    !
    ! Method:
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !



    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pkcnst ! Dimension construct for 3-D arrays
    INTEGER, INTENT (in) :: pkdim  ! Vertical dimension

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: fb(plond,pkdim*pkcnst,platd) ! Constituents array

    !  Local scalars:
    INTEGER :: i, j, k, m


    !  Executable statements

    IF (lc%nprocb > 1) RETURN

    DO m=1,pkcnst

      ! Fill west edge points.

      IF (nxpt >= 1) THEN
        DO j = 1, platd
          DO k = 1, pkdim
            DO i = 1, nxpt
              fb(i,k+(m-1)*pkdim,j) = fb(i+plon,k+(m-1)*pkdim,j)
            END DO
          END DO
        END DO
      END IF

      ! Fill east edge points

      DO j = 1, platd
        DO k = 1, pkdim
          DO i = i2pi, plond
            fb(i,k+(m-1)*pkdim,j) = fb(i-plon,k+(m-1)*pkdim,j)
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE extx

  SUBROUTINE extys(pkcnst,pkdim,fb)

    ! Description:
    !
    ! Fill latitude extensions of a scalar extended array.
    !
    ! Method:
    !
    ! This is done in 2 steps:
    !
    ! 1) Interpolate to the pole points; use the mean field value on the
    !    Gaussian latitude closest to the pole.
    ! 2) Add latitude lines beyond the poles.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    ! Parallel version:  T. Diehl, July 1999
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pkcnst, pkdim

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: fb(plond,pkdim,pkcnst,platd,2)

    ! pkcnst  Number of tracers
    ! pkdim   Vertical dimension
    ! fb      Output is same as on entry except with the pole latitude
    !         and extensions beyond it filled.

    ! Local arrays:
    REAL(dp) :: fb_buf(plon,2,pkdim), zbuf(2,pkdim)
    REAL(dp) :: fb_buf_g((plon+1)*2*pkdim,lc%nprocb)
    REAL(dp) :: fb_g(plon+1,2,pkdim,lc%nprocb)
    ! "plon+1" in fb_buf_g: need max. local #of lons
    REAL(dp) :: fb_g_one(lc%nlon,2,pkdim)
    REAL(dp) :: fbuf1(plon*pkdim*(nxpt+jintmx-1)*2,lc%nprocb)
    REAL(dp) :: fbuf2(plon*pkdim*(nxpt+jintmx-1)*2,lc%nprocb)
    INTEGER :: tagtable(lc%nprocb), pe_id2(lc%nprocb), pe_id(lc%nprocb)
    INTEGER :: nlons(lc%nprocb), ilon(plon,lc%nprocb)

    ! Local scalars:
    REAL(dp) :: zave         ! accumulator for zonal averaging
    REAL(dp) :: fnlon        ! REAL(dp)(lc%nlon)
    INTEGER :: i, ig, j, k
    INTEGER :: tag_base, tagcount, jproc, isource, src_ln
    INTEGER :: src, ig_mirror, node, inode, il, npes, idest
    INTEGER :: imeslen, jj, ii, ip, m
    INTEGER :: tagr1, mcountr1, sourcer1, tags1, mcounts1, dests1
    INTEGER :: tagr2, mcountr2, sourcer2, tags2, mcounts2, dests2
    INTEGER :: tagr3, mcountr3, sourcer3, tags3, mcounts3, dests3
    INTEGER :: len

    !  Intrinsic functions
    INTRINSIC COUNT, PACK, SIZE


    !  Executable statements

    ! return, if I am not a pole PE
    ! This version currently works only if "nxpt+jintmx-1 <= plat/2";
    ! (nxpt+jintmx-1 = the lines beyond the pole line)
    ! otherwise, communication with non pole PEs will be necessary

    IF (lc%set_a /= 1) RETURN

    fnlon = REAL(lc%nlon,dp)
    tag_base = 400
    fb_buf_g(:,:) = 0.0_dp
    fb_buf(:,:,:) = 0.0_dp
    fb_g(:,:,:,:) = 0.0_dp

    ! loop over fields

    DO m=1,pkcnst

      ! Fill north and south pole line.

      IF (lc%set_a == 1) THEN  ! if I am a pole PE

        DO k = 1, pkdim
          ig = 0
          DO i = istart, istop
            ig = ig + 1

            IF (lc%set_b /= 1) THEN
              ! north
              fb_buf(ig,2,k) = fb(i,k,m,jn,2)
              ! south
              fb_buf(ig,1,k) = fb(i,k,m,js,1)
            ELSE
              ! north
              fb_g(ig,2,k,1) = fb(i,k,m,jn,2)
              ! south
              fb_g(ig,1,k,1) = fb(i,k,m,js,1)
            END IF

          END DO
        END DO

        IF (lc%nprocb > 1) THEN
          IF (lc%set_b /= 1) THEN
            tags1 = tag_base + m
            mcounts1 = SIZE(fb_buf)
            dests1 = gc(lc%mapmesh(1,lc%set_a))%pe
            CALL p_isend(fb_buf,dests1,tags1,p_count=mcounts1)
          ELSE
            tagcount = 1
            tagtable(1) = tag_base + m
            DO jproc = 1,lc%nprocb - 1
              CALL p_probe(tagcount,tagtable,sourcer1,tagr1,mcountr1)
              isource = indx(sourcer1,gc)
              src_ln = gc(isource)%set_b
              CALL p_recv(fb_buf_g(1,src_ln),sourcer1,tagr1,p_count=mcountr1)

              len=gc(isource)%nglon
              ig = 0
              DO k=1,pkdim
                fb_g(1:len,1,k,src_ln) = fb_buf_g(ig+1:ig+len,src_ln)
                ig =ig + len
                fb_g(1:len,2,k,src_ln) = fb_buf_g(ig+1:ig+len,src_ln)
                ig =ig + len
              END DO

            END DO
          END IF
        END IF

        IF (lc%set_b == 1) THEN
          DO k = 1,pkdim
            ii = 0
            DO src_ln = 1,lc%nprocb
              src = lc%mapmesh(src_ln,lc%set_a)
              DO i = 1,gc(src)%nglon
                ii = ii + 1
                fb_g_one(ii,2,k) = fb_g(i,2,k,src_ln)
                fb_g_one(ii,1,k) = fb_g(i,1,k,src_ln)
              END DO
            END DO
          END DO

          zbuf(:,:) = 0.0_dp
          DO k = 1,pkdim
            DO i = 1,lc%nlon
              ! north
              zbuf(2,k) = zbuf(2,k) + fb_g_one(i,2,k)
              ! south
              zbuf(1,k) = zbuf(1,k) + fb_g_one(i,1,k)
            END DO
          END DO

          IF (lc%nprocb > 1) THEN
            DO jproc = 1,lc%nprocb
              IF (gc(lc%mapmesh(jproc,lc%set_a))%pe /= lc%pe) THEN
                dests2 = gc(lc%mapmesh(jproc,lc%set_a))%pe
                tags2 = tag_base + pkcnst + m
                mcounts2 = SIZE(zbuf)
                CALL p_isend(zbuf,dests2,tags2,p_count=mcounts2)
              END IF
            END DO
          END IF
        ELSE
          tagcount = 1
          tagtable(1) = tag_base + pkcnst + m
          CALL p_probe(tagcount,tagtable,sourcer2,tagr2,mcountr2)
          CALL p_recv(zbuf,sourcer2,tagr2,p_count=mcountr2)
        END IF

        DO k = 1,pkdim
          ! north
          zave = zbuf(2,k)/fnlon
          DO i = istart, istop
            fb(i,k,m,jn+1,2) = zave
          END DO
          ! south
          zave = zbuf(1,k)/fnlon
          DO i = istart, istop
            fb(i,k,m,js-1,1) = zave
          END DO
        END DO

        call p_wait

      END IF

    END DO

    ! Fill lines beyond pole line

    ! Case only one pole PE
    IF (lc%nprocb == 1) THEN
      DO m=1,pkcnst  ! loop over fields
        ! Fill northern lines beyond pole line.
        DO j = jn + 2, platd
          DO k = 1, pkdim
!DIR$ IVDEP
!OCL NOVREC
            DO i = istart, istart + plono2 - 1
              fb(i,k,m,j,2) = fb(plono2+i,k,m,2*jn+2-j,2)
              fb(plono2+i,k,m,j,2) = fb(i,k,m,2*jn+2-j,2)
            END DO
          END DO
        END DO

        ! Fill southern lines beyond pole line.
        DO j = 1, js - 2
          DO k = 1, pkdim
!DIR$ IVDEP
!OCL NOVREC
            DO i = istart, istart + plono2 - 1
              fb(i,k,m,j,1) = fb(plono2+i,k,m,2*js-2-j,1)
              fb(plono2+i,k,m,j,1) = fb(i,k,m,2*js-2-j,1)
            END DO
          END DO
        END DO
      END DO

    ELSE  ! Case more than one pole PE

      pe_id2(:) = -999
      pe_id(:) = -999
      nlons(:) = 0
      ilon(:,:) = -999

      il = 0
      DO ig = lc%glons(1),lc%glone(1)
        il = il + 1
        IF ((ig >= 1) .AND. (ig <= lc%nlon/2)) THEN
          ig_mirror = lc%nlon/2 + ig
        ELSE
          ig_mirror = ig - lc%nlon/2
        END IF

        DO jproc = 1,lc%nprocb
          inode = lc%mapmesh(jproc,lc%set_a)
          node = gc(inode)%pe
          IF (node /= lc%pe) THEN
            IF ((ig_mirror >= gc(inode)%glons(1)) .AND. &
                 (ig_mirror <= gc(inode)%glone(1))) THEN
              pe_id2(jproc) = node
              nlons(jproc) = nlons(jproc) + 1
              ilon(nlons(jproc),jproc) = istart + il - 1
              EXIT
            END IF
          END IF
        END DO
      END DO

      npes = COUNT(pe_id2 /= -999)
      pe_id(:npes) = PACK(pe_id2,pe_id2 /= -999)

      fbuf1(:,:) = 0.0_dp
      fbuf2(:,:) = 0.0_dp

      DO m=1,pkcnst ! loop over fields

        ! send part
        DO ip = 1,npes
          dests3 = pe_id(ip)
          idest = indx(dests3,gc)
          jproc = gc(idest)%set_b
          imeslen = 0

          ! north pole
          jj = 0
          DO j = jn+2,platd
            jj = jj+1
            DO k=1,pkdim
!DIR$ IVDEP
!OCL NOVREC
              DO i=1,nlons(jproc)
                ii = ilon(i,jproc)
                imeslen = imeslen + 1
                fbuf1(imeslen,ip) = fb(ii,k,m,2*jn+2-j,2)
              END DO
            END DO
          END DO

          ! south pole
          jj = 0
          DO j = 1,js - 2
            jj = jj + 1
            DO k = 1,pkdim
!DIR$ IVDEP
!OCL NOVREC
              DO i=1,nlons(jproc)
                ii = ilon(i,jproc)
                imeslen = imeslen + 1
                fbuf1(imeslen,ip) = fb(ii,k,m,2*js - 2 -j,1)
              END DO
            END DO
          END DO

          tags3 = tag_base + 2*pkcnst + 10 + (dests3+1)*m
          mcounts3 = imeslen
          ! Note that mcount=0 if jn+2>platd or js-2<1
          CALL p_isend(fbuf1(1,ip),dests3,tags3,p_count=mcounts3)
        END DO

        ! receive part
        tagcount = npes
        DO ip = 1,npes
          tagtable(ip) = tag_base + 2*pkcnst + 10 + (lc%pe + 1)*m
        END DO

        DO ip = 1,npes
          CALL p_probe(tagcount,tagtable,sourcer3,tagr3,mcountr3)
          CALL p_recv(fbuf2(1,ip),sourcer3,tagr3,p_count=mcountr3)

          isource = indx(sourcer3,gc)
          jproc = gc(isource)%set_b
          imeslen = 0

          ! north pole
          jj = 0
          DO j = jn+2,platd
            jj = jj + 1
            DO k = 1,pkdim
!DIR$ IVDEP
!OCL NOVREC
              DO i = 1,nlons(jproc)
                ii = ilon(i,jproc)
                imeslen = imeslen + 1
                fb(ii,k,m,j,2) = fbuf2(imeslen,ip)
              END DO
            END DO
          END DO

          ! south pole
          jj = 0
          DO j = 1,js-2
            jj = jj + 1
            DO k = 1,pkdim
!DIR$ IVDEP
!OCL NOVREC
              DO i = 1,nlons(jproc)
                ii = ilon(i,jproc)
                imeslen = imeslen + 1
                fb(ii,k,m,j,1) = fbuf2(imeslen,ip)
              END DO
            END DO
          END DO

        END DO

        call p_wait

      END DO

    END IF

  END SUBROUTINE extys

  SUBROUTINE extyv(pkdim,coslam,sinlam,vb)

    ! Description:
    !
    ! Fill latitude extensions of a vector component extended array.
    !
    ! Method:
    !
    ! This is done in 2 steps:
    !
    ! 1) Interpolate to the pole points; use coefficients for zonal wave
    !    number 1 on the Gaussian latitude closest to the pole.
    ! 2) Add latitude lines beyond the poles.
    !
    ! Dimensioning construct for 3-D arrays
    ! Vertical dimension
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    ! Parallel version:  T. Diehl, July 1999
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pkdim

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: coslam(plon), sinlam(plon)

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: vb(plond,pkdim,platd,2)

    ! pkdim   Vertical dimension
    ! coslam  Cos of longitude at x-grid points (global grid).
    ! sinlam  Sin of longitude at x-grid points (global grid).
    ! vb      Output is same as on entry except with the pole latitude
    !         and extensions beyond it filled.

    ! Local arrays:
    REAL(dp) :: vb_buf(plon,2,2,pkdim), zbuf(2,2,pkdim)
    REAL(dp) :: vb_buf_g((plon+1)*2*2*pkdim,lc%nprocb)
    REAL(dp) :: vb_g(plon+1,2,2,pkdim,lc%nprocb)
    ! "plon+1" in vb_buf_g: need max. local #of lons
    REAL(dp) :: vb_g_one(lc%nlon,2,2,pkdim)
    REAL(dp) :: vbuf1(plon*pkdim*(nxpt+jintmx-1)*2,lc%nprocb)
    REAL(dp) :: vbuf2(plon*pkdim*(nxpt+jintmx-1)*2,lc%nprocb)
    INTEGER :: tagtable(lc%nprocb), pe_id2(lc%nprocb), pe_id(lc%nprocb)
    INTEGER :: nlons(lc%nprocb), ilon(plon,lc%nprocb)

    ! Local scalars:
    REAL(dp) :: fnlono2       ! REAL(dp)(lc%nlon/2)
    REAL(dp) :: zavec         ! accumulator for wavenumber 1
    REAL(dp) :: zaves         ! accumulator for wavenumber 1
    INTEGER :: i, ig, j, k
    INTEGER :: tag_base, tagcount, jproc, isource, src_ln
    INTEGER :: src, ig_mirror, node, inode, il, npes, idest
    INTEGER :: imeslen, jj, ii, ip
    INTEGER :: tagr1, mcountr1, sourcer1, tags1, mcounts1, dests1
    INTEGER :: tagr2, mcountr2, sourcer2, tags2, mcounts2, dests2
    INTEGER :: tagr3, mcountr3, sourcer3, tags3, mcounts3, dests3
    INTEGER :: len

    !  Intrinsic functions
    INTRINSIC COUNT, SIZE, PACK


    !  Executable statements

    ! return, if I am not a pole PE
    ! This version currently works only if "nxpt+jintmx-1 <= plat/2";
    ! (nxpt+jintmx-1 = the lines beyond the pole line)
    ! otherwise, communication with non pole PEs will be necessary

    IF (lc%set_a /= 1) RETURN

    fnlono2 = REAL(lc%nlon/2,dp)
    tag_base = 400
    vb_buf_g(:,:)   = 0.0_dp
    vb_buf(:,:,:,:) = 0.0_dp
    vb_g(:,:,:,:,:) = 0.0_dp

    ! Fill north and south pole line.

    IF (lc%set_a == 1) THEN   ! if I am a pole PE ...

      DO k = 1, pkdim
        ig = 0
        DO i = istart, istop
          ig = ig + 1

          IF (lc%set_b /= 1) THEN
            ! north
            vb_buf(ig,1,2,k) = vb(i,k,jn,2)*coslam(ig)
            vb_buf(ig,2,2,k) = vb(i,k,jn,2)*sinlam(ig)

            ! south
            vb_buf(ig,1,1,k) = vb(i,k,js,1)*coslam(ig)
            vb_buf(ig,2,1,k) = vb(i,k,js,1)*sinlam(ig)
          ELSE
            ! north
            vb_g(ig,1,2,k,1) = vb(i,k,jn,2)*coslam(ig)
            vb_g(ig,2,2,k,1) = vb(i,k,jn,2)*sinlam(ig)

            ! south
            vb_g(ig,1,1,k,1) = vb(i,k,js,1)*coslam(ig)
            vb_g(ig,2,1,k,1) = vb(i,k,js,1)*sinlam(ig)
          END IF

        END DO
      END DO

      IF (lc%nprocb > 1) THEN
        IF (lc%set_b /= 1) THEN
          tags1 = tag_base + 1
          mcounts1 = SIZE(vb_buf)
          dests1 = gc(lc%mapmesh(1,lc%set_a))%pe
          CALL p_isend(vb_buf,dests1,tags1,p_count=mcounts1)
        ELSE
          tagcount = 1
          tagtable(1) = tag_base + 1
          DO jproc = 1,lc%nprocb - 1
            CALL p_probe(tagcount,tagtable,sourcer1,tagr1,mcountr1)
            isource = indx(sourcer1,gc)
            src_ln = gc(isource)%set_b
            CALL p_recv(vb_buf_g(1,src_ln),sourcer1,tagr1,p_count=mcountr1)

            len=gc(isource)%nglon
            ig = 0
            DO k=1,pkdim
              vb_g(1:len,1,1,k,src_ln) = vb_buf_g(ig+1:ig+len,src_ln)
              ig =ig + len
              vb_g(1:len,2,1,k,src_ln) = vb_buf_g(ig+1:ig+len,src_ln)
              ig =ig + len
              vb_g(1:len,1,2,k,src_ln) = vb_buf_g(ig+1:ig+len,src_ln)
              ig =ig + len
              vb_g(1:len,2,2,k,src_ln) = vb_buf_g(ig+1:ig+len,src_ln)
              ig =ig + len
            END DO

          END DO
        END IF
      END IF

      IF (lc%set_b == 1) THEN
        DO k = 1,pkdim
          ii = 0
          DO src_ln = 1,lc%nprocb
            src = lc%mapmesh(src_ln,lc%set_a)
            DO i = 1,gc(src)%nglon
              ii = ii + 1
              vb_g_one(ii,1,2,k) = vb_g(i,1,2,k,src_ln)
              vb_g_one(ii,2,2,k) = vb_g(i,2,2,k,src_ln)
              vb_g_one(ii,1,1,k) = vb_g(i,1,1,k,src_ln)
              vb_g_one(ii,2,1,k) = vb_g(i,2,1,k,src_ln)
            END DO
          END DO
        END DO

        zbuf(:,:,:) = 0.0_dp
        DO k = 1,pkdim
          DO i = 1,lc%nlon
            ! north
            zbuf(1,2,k) = zbuf(1,2,k) + vb_g_one(i,1,2,k)
            zbuf(2,2,k) = zbuf(2,2,k) + vb_g_one(i,2,2,k)
            ! south
            zbuf(1,1,k) = zbuf(1,1,k) + vb_g_one(i,1,1,k)
            zbuf(2,1,k) = zbuf(2,1,k) + vb_g_one(i,2,1,k)
          END DO
        END DO

        IF (lc%nprocb > 1) THEN
          DO jproc = 1,lc%nprocb
            IF (gc(lc%mapmesh(jproc,lc%set_a))%pe /= lc%pe) THEN
              dests2 = gc(lc%mapmesh(jproc,lc%set_a))%pe
              tags2 = tag_base + 2
              mcounts2 = SIZE(zbuf)
              CALL p_isend(zbuf,dests2,tags2,p_count=mcounts2)
            END IF
          END DO
        END IF
      ELSE
        tagcount = 1
        tagtable(1) = tag_base + 2
        CALL p_probe(tagcount,tagtable,sourcer2,tagr2,mcountr2)
        !    source = gc(lc%mapmesh(1,lc%set_a))%pe
        !    tagr2 = tag_base + 2
        !    mcountr2 = SIZE(zbuf)
        CALL p_recv(zbuf,sourcer2,tagr2,p_count=mcountr2)
      END IF

      DO k = 1,pkdim
        ! north
        zavec = zbuf(1,2,k)/fnlono2
        zaves = zbuf(2,2,k)/fnlono2
        ig = 0
        DO i = istart, istop
          ig = ig + 1
          vb(i,k,jn+1,2) = zavec*coslam(ig) + zaves*sinlam(ig)
        END DO
        ! south
        zavec = zbuf(1,1,k)/fnlono2
        zaves = zbuf(2,1,k)/fnlono2
        ig = 0
        DO i = istart, istop
          ig = ig + 1
          vb(i,k,js-1,1) = zavec*coslam(ig) + zaves*sinlam(ig)
        END DO
      END DO

    END IF

    !call p_barrier

    ! Fill lines beyond pole line
    ! Case only one pole PE
    IF (lc%nprocb == 1) THEN
      !     IF (jn+2 <= platd) THEN
      ! Fill northern lines beyond pole line.
      DO j = jn + 2, platd
        DO k = 1, pkdim
!DIR$ IVDEP
!OCL NOVREC
          DO i = istart, istart + plono2 - 1
            vb(i,k,j,2) = -vb(plono2+i,k,2*jn+2-j,2)
            vb(plono2+i,k,j,2) = -vb(i,k,2*jn+2-j,2)
          END DO
        END DO
      END DO
      !     END IF

      !     IF (js-2 >= 1) THEN
      ! Fill southern lines beyond pole line.
      DO j = 1, js - 2
        DO k = 1, pkdim
!DIR$ IVDEP
!OCL NOVREC
          DO i = istart, istart + plono2 - 1
            vb(i,k,j,1) = -vb(plono2+i,k,2*js-2-j,1)
            vb(plono2+i,k,j,1) = -vb(i,k,2*js-2-j,1)
          END DO
        END DO
      END DO
      !     END IF

    ELSE  ! Case more than one pole PE

      pe_id2(:) = -999
      pe_id(:) = -999
      nlons(:) = 0
      ilon(:,:) = -999

      il = 0
      DO ig = lc%glons(1),lc%glone(1)
        il = il + 1
        IF ((ig >= 1) .AND. (ig <= lc%nlon/2)) THEN
          ig_mirror = lc%nlon/2 + ig
        ELSE
          ig_mirror = ig - lc%nlon/2
        END IF

        DO jproc = 1,lc%nprocb
          inode = lc%mapmesh(jproc,lc%set_a)
          node = gc(inode)%pe
          IF (node /= lc%pe) THEN
            IF ((ig_mirror >= gc(inode)%glons(1)) .AND. &
                 (ig_mirror <= gc(inode)%glone(1))) THEN
              pe_id2(jproc) = node
              nlons(jproc) = nlons(jproc) + 1
              ilon(nlons(jproc),jproc) = istart + il - 1
              EXIT
            END IF
          END IF
        END DO
      END DO

      npes = COUNT(pe_id2 /= -999)
      pe_id(:npes) = PACK(pe_id2,pe_id2 /= -999)

      vbuf1(:,:) = 0.0_dp
      vbuf2(:,:) = 0.0_dp

      ! send part
      DO ip = 1,npes
        dests3 = pe_id(ip)
        idest = indx(dests3,gc)
        jproc = gc(idest)%set_b
        imeslen = 0

        ! north pole
        jj = 0
        DO j = jn+2,platd
          jj = jj+1
          DO k=1,pkdim
!DIR$ IVDEP
!OCL NOVREC
            DO i=1,nlons(jproc)
              ii = ilon(i,jproc)
              imeslen = imeslen + 1
              vbuf1(imeslen,ip) = -vb(ii,k,2*jn+2-j,2)
            END DO
          END DO
        END DO

        ! south pole
        jj = 0
        DO j = 1,js - 2
          jj = jj + 1
          DO k = 1,pkdim
!DIR$ IVDEP
!OCL NOVREC
            DO i=1,nlons(jproc)
              ii = ilon(i,jproc)
              imeslen = imeslen + 1
              vbuf1(imeslen,ip) = -vb(ii,k,2*js - 2 -j,1)
            END DO
          END DO
        END DO

        !        tags3 = tag_base + 10 + ip
        tags3 = tag_base + 10 + dests3
        mcounts3 = imeslen
        ! Note that mcount=0 if jn+2>platd or js-2<1
        CALL p_isend(vbuf1(1,ip),dests3,tags3,p_count=mcounts3)

      END DO

      ! receive part
      tagcount = npes
      DO ip = 1,npes
        tagtable(ip) = tag_base + 10 + lc%pe
      END DO

      DO ip = 1,npes
        CALL p_probe(tagcount,tagtable,sourcer3,tagr3,mcountr3)
        CALL p_recv(vbuf2(1,ip),sourcer3,tagr3,p_count=mcountr3)

        isource = indx(sourcer3,gc)
        jproc = gc(isource)%set_b
        imeslen = 0

        ! north pole
        jj = 0
        DO j = jn+2,platd
          jj = jj + 1
          DO k = 1,pkdim
!DIR$ IVDEP
!OCL NOVREC
            DO i = 1,nlons(jproc)
              ii = ilon(i,jproc)
              imeslen = imeslen + 1
              vb(ii,k,j,2) = vbuf2(imeslen,ip)
            END DO
          END DO
        END DO

        ! south pole
        jj = 0
        DO j = 1,js-2
          jj = jj + 1
          DO k = 1,pkdim
!DIR$ IVDEP
!OCL NOVREC
            DO i = 1,nlons(jproc)
              ii = ilon(i,jproc)
              imeslen = imeslen + 1
              vb(ii,k,j,1) = vbuf2(imeslen,ip)
            END DO
          END DO
        END DO

      END DO

    END IF

    CALL p_wait

  END SUBROUTINE extyv

  SUBROUTINE g2spos(ih, dttmp, lam, phib, phi, cosphi, sinphi, upr, vpr,&
       lamgc, phigc, lamsc, phisc)

    ! Description:
    !
    ! Transform position coordinates for a set of points.
    !
    ! Method:
    !
    ! Transform position coordinates for a set of points, each of which is
    ! associated with a grid point in a global latitude slice, from local
    ! geodesic to spherical coordinates.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Local parameters:
    ! 1. - eps, eps from machine precision
    REAL(dp), PARAMETER :: fac = 1.0_dp-1.0E-12

    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: cosphi, dttmp, phi, sinphi

    INTEGER, intent (in) :: ih

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: lam(plon), lamgc(pgls), phib(platd), phigc(pgls), &
         upr(pgls), vpr(pgls)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: lamsc(plon,plev), phisc(pgls)

    ! dttmp  Time step over which midpoint/endpoint trajectory is
    !        calculated (seconds).
    ! lam    Longitude coordinates of the global grid points in spherical
    !        system.  The grid points in the global array are the reference
    ! points for the local geodesic systems.
    ! phib   Latitude values for the extended grid.
    ! phi    Latitude coordinate (in the global grid) of the current
    !        latitude slice.
    ! cosphi cos( phi )
    ! sinphi sin( phi )
    ! upr    zonal      velocity at departure point in local geodesic coord
    ! vpr    Meridional velocity at departure point in local geodesic coord
    ! lamgc  Longitude coordinate of points in geodesic coordinates.
    ! phigc  Latitude coordinate of points in geodesic coordinates.
    ! lamsc  Longitude coordinate of points in spherical coordinates.
    ! phisc  Latitude coordinate of points in spherical coordinates.

    !  Local scalars:
    REAL(dp) :: clamgc        ! cos(lamgc)
    REAL(dp) :: coeff         ! tmp variable
    REAL(dp) :: cphigc        ! cos(phigc)
    REAL(dp) :: distmx        ! max distance
    REAL(dp) :: phipi2        ! tmp variable
    REAL(dp) :: pi            ! 4.*atan(1.)
    REAL(dp) :: pi2           ! pi/2
    REAL(dp) :: sgnphi        ! holds sign of phi
    REAL(dp) :: slam2         ! sin(lamgc)**2
    REAL(dp) :: sphigc        ! sin(phigc)
    REAL(dp) :: twopi         ! 2.*pi
    INTEGER :: i, ii, k
    INTEGER :: nval       ! number of values returned from whenfgt

    !  Local arrays:
    REAL(dp) :: dist(pgls)    ! approx. distance traveled along traj.
    REAL(dp) :: dlam(pgls)    ! zonal extent of trajectory
    REAL(dp) :: slamgc(pgls)  ! sin(lamgc)
    INTEGER :: indx(pgls) ! index holder

    !!LK  !  External subroutines
    !!LK  EXTERNAL whenfgt

    !  Intrinsic functions
!    INTRINSIC ABS, ASIN, ATAN, COS, SIGN, SIN, SQRT
    INTRINSIC ASIN, ATAN, COS, SIGN, SIN, SQRT

    !  Executable statements

    pi     = 4.0_dp*ATAN(1.0_dp)
    twopi  = pi*2.0_dp
    pi2    = pi*0.5_dp
    coeff  = (1.1_dp*dttmp)**2
    distmx = (SIGN(pi2,phi)-phi)**2/coeff
    sgnphi = SIGN(1.0_dp,phi)

    DO i = 1, pgls
      sphigc = SIN(phigc(i))
      cphigc = COS(phigc(i))
      slamgc(i) = SIN(lamgc(i))
      clamgc = COS(lamgc(i))
      phisc(i) = ASIN((sphigc*cosphi+cphigc*sinphi*clamgc)*fac)

      ! IF (ABS(phisc(i))>=phib(j1+plat)*fac) phisc(i) = SIGN(phib(j1+plat), &
      !      phisc(i))*fac

      ! IF (ABS(phisc(i))>=phib(j1+plato2)*fac) phisc(i) = SIGN(phib(j1+plato2), &
      !      phisc(i))*fac

      ! IF(phisc(i)>=phib(j1+plato2))   phisc(i) = phib(j1+plato2)*fac

      ! parallel version
      IF (lc%set_a ==1) THEN  ! If I am a pole PE
        IF (ih == 2) THEN
          IF (phisc(i) >= phib(j1+plato2)*fac) phisc(i) = phib(j1+plato2)*fac
        ELSE
          IF (phisc(i) <= phib(j1-1)*fac) phisc(i) = phib(j1-1)*fac
        END IF
      END IF

      dlam(i) = ASIN((slamgc(i)*cphigc/COS(phisc(i)))*fac)

      ! Compute estimated trajectory distance based upon winds alone

      dist(i) = upr(i)**2 + vpr(i)**2
    END DO

    ! Determine which trajectories may have crossed over pole

    CALL whenfgt(pgls,dist,1,distmx,indx,nval)

    ! Check that proper branch of arcsine is used for calculation of
    ! dlam for those trajectories which may have crossed over pole.

!DIR$ IVDEP
!OCL NOVREC
    DO ii = 1, nval
      i = indx(ii)
      slam2 = slamgc(i)**2
      phipi2 = ASIN((SQRT((slam2-1.0_dp)/(slam2-1.0_dp/cosphi**2)))*fac)
      IF (sgnphi*phigc(i) > phipi2) THEN
        dlam(i) = SIGN(pi,lamgc(i)) - dlam(i)
      END IF
    END DO

    DO k = 1, plev
      DO i = 1, plon
        lamsc(i,k) = lam(i) + dlam((k-1)*plon+i)

        ! Restrict longitude to be in the range [0, twopi).
        ! But: Parallel code allows extension beyond twopi

        IF (debug_seriell .and. lc%nprocb == 1) THEN
          IF (lamsc(i,k) >= twopi) lamsc(i,k) = lamsc(i,k) - twopi
          IF (lamsc(i,k) < 0.0_dp) lamsc(i,k) = lamsc(i,k) + twopi
        END IF

      END DO
    END DO

  END SUBROUTINE g2spos

  SUBROUTINE herxin(pf,fb,fxl,fxr,xdp,idp,jdp,fint)

    ! Description:
    !
    ! Interpolate to its x value at each latitude
    !
    ! Method:
    !  
    ! For each departure point in the latitude slice being forecasted,
    ! interpolate (using equally spaced Hermite cubic formulas) to its
    ! x value at each latitude required for later interpolation in the y
    ! direction.
    !
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! Modified:          U. Schlese,  DKRZ - Hamburg, May 1994
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    ! 
    ! for more details see file AUTHORS
    !
    
    !  Scalar arguments with intent(In):
    INTEGER, INTENT (IN) :: pf

    !  Array arguments with intent(In):
    REAL(dp), INTENT (IN) :: fb(plond,plev,pf,platd), &
         fxl(plond,plev,pf,platd), fxr(plond,plev,pf,platd), &
         xdp(plon,plev)
    
    INTEGER, INTENT (IN) :: idp(plon,plev), jdp(plon,plev)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (OUT) :: fint(plon,plev,ppdy,pf)

    ! pf      Number of fields being interpolated.
    ! fb      extended array of data to be interpolated.
    ! fxl     x derivatives at the left edge of each interval containing
    !         the departure point
    ! fxr     x derivatives at the right edge of each interval containing
    !         the departure point
    ! xdp     xdp(i,k) is the x-coordinate (extended grid) of the
    !         departure point that corresponds to global grid point (i,k)
    !         in the latitude slice being forecasted.
    ! idp     idp(i,k) is the index of the x-interval (extended grid) that
    !         contains the departure point corresponding to global grid
    !         point (i,k) in the latitude slice being forecasted.
    !         Note that x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
    ! jdp     jdp(i,k) is the index of the y-interval (extended grid) that
    !         contains the departure point corresponding to global grid
    !         point (i,k) in the latitude slice being forecasted.
    !         Suppose yb contains the y-coordinates of the extended array
    !         and ydp(i,k) is the y-coordinate of the departure point
    !         corresponding to grid point (i,k).  Then,
    !         yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
    ! fint    (fint(i,k,j,n),j=1,ppdy) contains the x interpolants at each
    !         latitude needed for the y derivative estimates at the
    !         endpoints of the interval that contains the departure point
    !         for grid point (i,k).  The last index of fint allows for
    !         interpolation of multiple fields.
    
    !  Local scalars: 
    REAL(dp) :: dx        ! x-increment
    REAL(dp) :: rdx       ! 1./dx
    ! interpolation coefficient
    REAL(dp) :: xl, xr
    REAL(dp) :: pi, xx
    INTEGER  :: i, k, m
    
    !  Local arrays: 
    ! interpolation coefficient
    REAL(dp) :: dhl(plon,plev), dhr(plon,plev), hl(plon,plev), hr(plon,plev)
    INTEGER  :: id2(plon,plev), ioff

    !  Executable statements 
    ! dx = x(nxpt+2) - x(nxpt+1)
    ! parallel version

    pi = 4.0_dp*ATAN(1.0_dp)
    dx = 2.0_dp*pi/REAL(lc%nlon,dp)
    ioff = lc%glons(1) - istart

    rdx = 1.0_dp/dx
    DO k = 1, plev
      DO i = 1, plon
        xx = REAL(idp(i,k)+ioff,dp)*dx         ! same as in grdxy
        xl = ( xx           -xdp(i,k))*rdx
        !     xl = ( x(idp(i,k)+1)-xdp(i,k))*rdx  ! old version, relies on extended grid
        xr = 1.0_dp - xl
        hl(i,k) = (3.0_dp-2.0_dp*xl)*xl**2
        hr(i,k) = (3.0_dp-2.0_dp*xr)*xr**2
        dhl(i,k) = -dx*(xl-1.0_dp)*xl**2
        dhr(i,k) = dx*(xr-1.0_dp)*xr**2
      END DO
    END DO

    ! access inner domain only, if 1 PE in e-W direction

    id2 = idp
    IF (lc% nprocb == 1) THEN
      WHERE (id2 >  nxpt + plon) id2 = id2 - plon
      WHERE (id2 <= nxpt       ) id2 = id2 + plon
    ENDIF

    ! x interpolation at each latitude needed for y interpolation.
    ! Once for each field.

!vdir noloopchg
    DO m = 1, pf
      DO k = 1, plev
        DO i = 1, plon
          fint(i,k,1,m) = fb(id2(i,k),k,m,jdp(i,k)-1)*hl(i,k) +   &
                          fb(id2(i,k)+1,k,m,jdp(i,k)-1)*hr(i,k) + &
                          fxl(id2(i,k),k,m,jdp(i,k)-1)*dhl(i,k) + &
                          fxr(id2(i,k),k,m,jdp(i,k)-1)*dhr(i,k)
          fint(i,k,2,m) = fb(id2(i,k),k,m,jdp(i,k))*hl(i,k) +     &
                          fb(id2(i,k)+1,k,m,jdp(i,k))*hr(i,k) +   &
                          fxl(id2(i,k),k,m,jdp(i,k))* &
                          dhl(i,k) + fxr(id2(i,k),k,m,jdp(i,k))*dhr(i,k)
          fint(i,k,3,m) = fb(id2(i,k),k,m,jdp(i,k)+1)*hl(i,k) +   &
                          fb(id2(i,k)+1,k,m,jdp(i,k)+1)*hr(i,k) + &
                          fxl(id2(i,k),k,m,jdp(i,k)+1)*dhl(i,k) + &
                          fxr(id2(i,k),k,m,jdp(i,k)+1)*dhr(i,k)
          fint(i,k,4,m) = fb(id2(i,k),k,m,jdp(i,k)+2)*hl(i,k) +   &
                          fb(id2(i,k)+1,k,m,jdp(i,k)+2)*hr(i,k) + &
                          fxl(id2(i,k),k,m,jdp(i,k)+2)*dhl(i,k) + &
                          fxr(id2(i,k),k,m,jdp(i,k)+2)*dhr(i,k)
        END DO
      END DO
    END DO

  END SUBROUTINE herxin


  SUBROUTINE heryin(pf,fint,fyb,fyt,y,dy,ydp,jdp,fdp)

    ! Description:
    !
    ! Interpolate the x interpolants to the y value of the departure point
    !
    ! Method:
    !
    ! For each departure point in the latitude slice to be forecasted,
    ! interpolate (using unequally spaced Hermite cubic formulas) the
    ! x interpolants to the y value of the departure point.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! Modified:          U. Schlese,  DKRZ - Hamburg, May 1994
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pf

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: dy(platd), fint(pgls,ppdy,pf), fyb(pgls,pf), &
         &                      fyt(pgls,pf), y(platd), ydp(pgls)
    INTEGER, INTENT (in) :: jdp(pgls)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) ::  fdp(pgls,pf)

    ! pf      Number of fields being interpolated.
    ! fint    (fint(i,k,j,m),j=ppdy/2,ppdy/2 + 1) contains the x
    !         interpolants at the endpoints of the y-interval that
    !         contains the departure point for grid point (i,k).  The last
    !         index of fint allows for interpolation of multiple fields.
    !         fint is generated by a call to herxin.
    ! fyb     fyb(i,k,.) is the derivative at the "bottom" of the
    !         y-interval that contains the departure point of grid
    !         point (i,k).  fyb is generated by a call to cubydr.
    ! fyt     fyt(i,k,.) is the derivative at the "top" of the y-interval
    !         that contains the departure point of grid point (i,k).
    !         fyt is generated by a call to cubydr.
    ! y       y-coordinate (latitude) values in the extended array.
    ! dy      Increment in the y-coordinate value for each interval in the
    !         extended array.
    ! ydp     ydp(i,k) is the y-coordinate of the departure point that
    !         corresponds to global grid point (i,k) in the latitude slice
    !         being forecasted.
    ! jdp     jdp(i,k) is the index of the y-interval that contains the
    !         departure point corresponding to global grid point (i,k) in
    !         the latitude slice being forecasted.
    !         Note that
    !         y(jdp(i,k)) .le. ydp(i,k) .lt. y(jdp(i,k)+1) .
    ! fdp     Horizontally interpolated field values at the departure point
    !         for the latitude slice being forecasted.

    !  Local scalars:
    ! latitude interval containing dep. pt.
    REAL(dp) :: dyj, yb, yt
    INTEGER :: i, jb, jt, m

    !  Local arrays:
    ! interpolation coefficients
    REAL(dp) :: dhb(pgls), dht(pgls), hb(pgls), ht(pgls)


    !  Executable statements
    jb = ppdy/2
    jt = jb + 1

    DO i = 1, pgls
      dyj = dy(jdp(i))
      yb = (y(jdp(i)+1)-ydp(i))/dyj
      yt = 1._dp - yb
      hb(i) = (3.0_dp-2.0_dp*yb)*yb**2
      ht(i) = (3.0_dp-2.0_dp*yt)*yt**2
      dhb(i) = -dyj*(yb-1.0_dp)*yb**2
      dht(i) = dyj*(yt-1.0_dp)*yt**2
    END DO

    ! Loop over fields.

    DO m = 1, pf
      DO i = 1, pgls
        fdp(i,m) = fint(i,jb,m)*hb(i) + fyb(i,m)*dhb(i) &
                 + fint(i,jt,m)*ht(i) + fyt(i,m)*dht(i)
      END DO
    END DO

  END SUBROUTINE heryin

  SUBROUTINE herzin(pkdim,pf,f,fst,fsb,psig,pdsig,psigdp,kdp,fdp)

    ! Description:
    !
    ! Interpolate field on vertical slice to vertical departure point
    !
    ! Method:
    !
    ! Interpolate field on vertical slice to vertical departure point using
    ! Hermite cubic interpolation.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pf, pkdim

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: pdsig(pkdim), f(plon,pkdim,pf), fsb(plon,pkdim,pf),  &
         fst(plon,pkdim,pf), psig(pkdim), psigdp(plon,plev)
    INTEGER, INTENT (in) :: kdp(plon,plev)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: fdp(plon,plev,pf)

    ! pkdim   Vertical dimension of vertical slice arrays.
    ! pf      Number of fields being interpolated.
    ! f       Vertical slice of data to be interpolated.
    ! fst     z-derivatives at the top edge of each interval contained in f
    ! fsb     z-derivatives at the bot edge of each interval contained in f
    ! sig     Sigma values corresponding to the vertical grid
    ! dsig    Increment in sigma value for each interval in vertical grid.
    ! sigdp   Sigma value at the trajectory midpoint or endpoint for each
    !         gridpoint in a vertical slice from the global grid.
    ! kdp     Vertical index for each gridpoint.  This index points into a
    !         vertical slice array whose vertical grid is given by sig.
    ! E.g.,   sig(kdp(i,j)) .le. sigdp(i,j) .lt. sig(kdp(i,j)+1) .
    ! fdp     Value of field at the trajectory midpoints or endpoints.

    !  Local scalars:
    ! vert interval containing the dep. pt.
    REAL(dp) :: dzk, zb, zt
    INTEGER :: i, k, m

    !  Local arrays:
    ! interpolation coefficients
    REAL(dp) :: dhb(plon,plev), dht(plon,plev), hb(plon,plev), ht(plon,plev)

!DIR$ NOBOUNDS

    !  Executable statements

    DO k = 1, plev
      DO i = 1, plon
        dzk      = pdsig(kdp(i,k))
        zt       = (psig(kdp(i,k)+1)-psigdp(i,k))/dzk
        zb       = 1.0_dp - zt
        ht(i,k)  = (3.0_dp-2.0_dp*zt)*zt**2
        hb(i,k)  = (3.0_dp-2.0_dp*zb)*zb**2
        dht(i,k) = -dzk*(zt-1.0_dp)*zt**2
        dhb(i,k) =  dzk*(zb-1.0_dp)*zb**2
      END DO
    END DO

    ! Loop over fields.

    DO m = 1, pf
      DO k = 1, plev
        DO i = 1, plon
          fdp(i,k,m) = f(i,kdp(i,k)  ,m)*ht(i,k)+fst(i,kdp(i,k),m)*dht(i,k) &
                     + f(i,kdp(i,k)+1,m)*hb(i,k)+fsb(i,kdp(i,k),m)*dhb(i,k)
        END DO
      END DO
    END DO

  END SUBROUTINE herzin

  SUBROUTINE hrintp(pf,fb,qxl,qxr,y,dy,wdy,xdp,ydp,idp,jdp,limitd, &
                    fint,fyb,fyt,fdp)

    ! Description:
    !
    ! Interpolate 2-d field to departure point
    !
    ! Method:
    !
    ! Interpolate 2-d field to departure point using tensor product
    ! Hermite cubic interpolation.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pf
    LOGICAL, INTENT (in) :: limitd    ! flag for shape-preservation

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: dy(platd), fb(plond,plev,pf,platd),  &
         qxl(plond,plev,pf,platd), qxr(plond,plev,pf,platd),  &
         wdy(4,2,platd), xdp(plon,plev), y(platd), ydp(plon,plev)
    INTEGER, INTENT (in) :: idp(plon,plev), jdp(plon,plev)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: fint(plon,plev,ppdy,pf), fyb(plon,plev,pf),  &
         fyt(plon,plev,pf), fdp(plon,plev,pf)

    ! pf      Number of fields being interpolated.
    ! fb      Extended array of data to be interpolated.
    ! qxl     x-derivatives at the left  edge of each interval containing
    !         the departure point.
    ! qxr     x-derivatives at the right edge of each interval containing
    !         the departure point.
    ! y       y-coordinate (latitude) values in the extended array.
    ! dy      Increment in the y-coordinate value for each interval in the
    !         extended array.
    ! wdy     Weights for Lagrange cubic derivative estimates on the
    !         unequally spaced y-grid.  If grid interval j (in extended
    !         array is surrounded by a 4 point stencil, then the
    !         derivative at the "bottom" of the interval uses the weights
    !         wdy(1,1,j),wdy(2,1,j), wdy(3,1,j), and wdy(4,1,j).  The
    !         derivative at the "top" of the interval uses wdy(1,2,j),
    !         wdy(2,2,j), wdy(3,2,j) and wdy(4,2,j).
    ! xdp     xdp(i,k) is the x-coordinate of the departure point that
    !         corresponds to global grid point (i,k) in the latitude slice
    !         being forecasted.
    ! ydp     ydp(i,k) is the y-coordinate of the departure point that
    !         corresponds to global grid point (i,k) in the latitude slice
    !         being forecasted.
    ! idp     idp(i,k) is the index of the x-interval that contains the
    !         departure point corresponding to global grid point (i,k) in
    !         the latitude slice being forecasted.
    !         Note that
    !         x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
    ! jdp     jdp(i,k) is the index of the y-interval that contains the
    !         departure point corresponding to global grid point (i,k) in
    !         the latitude slice being forecasted.
    !         Suppose yb contains the y-coordinates of the extended array
    !         and ydp(i,k) is the y-coordinate of the departure point
    !         corresponding to grid point (i,k).  Then,
    !         yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
    !         limitd  Logical flag to specify whether or not the y-derivatives
    !         will be limited.
    ! fint    WORK ARRAY, results not used on return
    ! fyb     WORK ARRAY, results not used on return
    ! fyt     WORK ARRAY, results not used on return
    ! fdp     Value of field at the horizontal departure points.

    !!LK  !  External aubroutines
    !!LK  EXTERNAL cubydr, herxin, heryin, limdy


    !  Executable atatements

    ! Hermite cubic interpolation to the x-coordinate of each
    ! departure point at each y-coordinate required to compute the
    ! y-derivatives.

    CALL herxin(pf,fb,qxl,qxr,xdp,idp,jdp,fint)

    ! Compute y-derivatives.

    CALL cubydr(pf,fint,wdy,jdp,fyb,fyt)
    IF (limitd) THEN
      CALL limdy(pf,fint,dy,jdp,fyb,fyt)
    END IF

    ! Hermite cubic interpolation in the y-coordinate.

    CALL heryin(pf,fint,fyb,fyt,y,dy,ydp,jdp,fdp)

  END SUBROUTINE hrintp

  SUBROUTINE hrvint(pf,fb,qxl,qxr,wdy,xdp,ydp,idp,jdp,fint,fdp)

    ! Description:
    !
    ! Interpolate 2-d field to departure point
    !
    ! Method:
    !
    ! Interpolate 2-d field to departure point using equivalent of tensor
    ! product Lagrange cubic interpolation. For economy, code actually
    ! uses Hermite cubic with Lagrange cubic derivative estimate in x
    ! and Lagrange cubic in y.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pf

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: fb(plond,plev,pf,platd), qxl(plond,plev,pf,platd), &
         qxr(plond,plev,pf,platd), wdy(4,2,platd), xdp(plon,plev),ydp(plon,plev)
    INTEGER :: idp(plon,plev), jdp(plon,plev)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: fint(plon,plev,ppdy,pf), fdp(plon,plev,pf)

    ! pf      Number of fields being interpolated.
    ! fb      Extended array of data to be interpolated.
    ! qxl     x-derivatives at the left  edge of each interval containing
    !         the departure point.
    ! qxr     x-derivatives at the right edge of each interval containing
    !         the departure point.
    ! wdy     Weights for Lagrange cubic interpolation on the
    !         unequally spaced y-grid.
    ! xdp     xdp(i,k) is the x-coordinate of the departure point that
    !         corresponds to global grid point (i,k) in the latitude slice
    !         being forecasted.
    ! ydp     ydp(i,k) is the y-coordinate of the departure point that
    !         corresponds to global grid point (i,k) in the latitude slice
    !         being forecasted.
    ! idp     idp(i,k) is the index of the x-interval that contains the
    !         departure point corresponding to global grid point (i,k) in
    !         the latitude slice being forecasted.
    !         Note that
    !         x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
    ! jdp     jdp(i,k) is the index of the y-interval that contains the
    !         departure point corresponding to global grid point (i,k) in
    !         the latitude slice being forecasted.
    !         Suppose yb contains the y-coordinates of the extended array
    !         and ydp(i,k) is the y-coordinate of the departure point
    !         corresponding to grid point (i,k).  Then,
    !         yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
    ! fint    WORK ARRAY, results not used on return
    ! fdp     Value of field at the departure points.

    !!LK  !  External subroutines
    !!LK  EXTERNAL herxin, lagyin


    !  Executable statements

    ! Hermite cubic interpolation to the x-coordinate of each
    ! departure point at each y-coordinate required to compute the
    ! y-interpolants.

    CALL herxin(pf,fb,qxl,qxr,xdp,idp,jdp,fint)

    ! Lagrange cubic interpolation in y

    CALL lagyin(pf,fint,wdy,ydp,jdp,fdp)

  END SUBROUTINE hrvint

  SUBROUTINE kdpfnd(pkdim,pmap,sig,sigdp,kdpmap,kdp)

    ! Description:
    !
    ! Determine vertical departure point indices.
    !
    ! Method:
    !
    ! Determine vertical departure point indices that point into a grid
    ! containing the full or half sigma levels.  Use an artificial evenly
    ! spaced vertical grid to map into the true model levels.
    !
    ! Indices are computed assuming the the sigdp values have
    ! been constrained so that sig(1) .le. sigdp(i,j) .lt. sig(pkdim).
    !
    ! Authors:
    !
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pkdim         ! dimension of "sig"
    INTEGER, INTENT (in) :: pmap          ! dimension of "kdpmap"

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: sig(pkdim)       ! vertical grid coordinates
    REAL(dp), INTENT (in) ::  sigdp(pgls)     ! vertical coords. of departure points
    INTEGER, INTENT (in) ::  kdpmap(pmap) ! array of model grid indices which
    ! are mapped into the artificial grid.

    !  Array arguments with intent(Out):
    INTEGER, INTENT (out) :: kdp(pgls)    ! vertical index for each dep. pt.

    !  Local scalars:
    REAL(dp) :: rdel          ! Reciprocal of interval in artificial g
    REAL(dp) :: sig1ln        ! ln (sig(1))
    INTEGER :: i, ii      ! indices

    !  Intrinsic functions
    INTRINSIC INT, LOG, MAX, MIN


    !  Executable statements
    rdel = REAL(pmap,dp)/(LOG(sig(pkdim))-LOG(sig(1)))
    sig1ln = LOG(sig(1))

    DO i = 1, pgls

      ! First guess of the departure point's location in the model grid

      ii = MAX(1,MIN(pmap,INT((LOG(sigdp(i))-sig1ln)*rdel+1.0_dp)))
      kdp(i) = kdpmap(ii)

      ! Determine if location is in next interval

      IF (sigdp(i)>=sig(kdp(i)+1)) THEN
        kdp(i) = kdp(i) + 1
      END IF
    END DO

  END SUBROUTINE kdpfnd

  SUBROUTINE lagyin(pf,fint,wdy,ydp,jdp,fdp)

    ! Description:
    !
    ! Interpolate the x interpolants to the y value of the departure point
    !
    ! Mehotd:
    !
    ! For each departure point in the latitude slice to be forecasted,
    ! interpolate (using unequally spaced Lagrange cubic formulas) the
    ! x interpolants to the y value of the departure point.
    !
    ! Authors:
    !
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER :: pf

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: fint(pgls,ppdy,pf), wdy(4,2,platd), ydp(pgls)
    INTEGER, INTENT (in) :: jdp(pgls)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: fdp(pgls,pf)

    ! pf      Number of fields being interpolated.
    ! fint    (fint(i,k,j,m),j=ppdy/2,ppdy/2 + 1) contains the x
    !         interpolants at the endpoints of the y-interval that contains
    !         the departure point for grid point (i,k).  The last index of
    !         fint allows for interpolation of multiple fields.  fint is
    !         generated by a call to herxin.
    ! wdy     Grid values and weights for Lagrange cubic interpolation on
    !         the unequally spaced y-grid.
    ! ydp     ydp(i,k) is the y-coordinate of the departure point that
    !         corresponds to global grid point (i,k) in the latitude slice
    !         being forecasted.
    ! jdp     jdp(i,k) is the index of the y-interval that contains the
    !         departure point corresponding to global grid point (i,k) in
    !         the latitude slice being forecasted.
    !         Note that
    !         y(jdp(i,k)) .le. ydp(i,k) .lt. y(jdp(i,k)+1) .
    ! fdp     Horizontally interpolated field values at the departure point
    !         for the latitude slice being forecasted.

    !  Local scalars:
    ! -- interpolation weights/coeffs.
    REAL(dp) :: coef12, coef34, ymy1, ymy2, ymy3, ymy4
    INTEGER :: i, m

    !  Local arrays:
    REAL(dp) :: term1(pgls), term2(pgls), term3(pgls), term4(pgls)


    !  Executable statements

    IF (ppdy /= 4) THEN
      WRITE (nout,'(A)') 'LAGYIN:Error:  ppdy .ne. 4'
      CALL finish('lagyin','Run terminated.')
    END IF

    DO i = 1, pgls
      ymy3 = ydp(i) - wdy(3,1,jdp(i))
      ymy4 = ydp(i) - wdy(4,1,jdp(i))
      coef12 = ymy3*ymy4
      ymy2 = ydp(i) - wdy(2,1,jdp(i))
      term1(i) = coef12*ymy2*wdy(1,2,jdp(i))
      ymy1 = ydp(i) - wdy(1,1,jdp(i))
      term2(i) = coef12*ymy1*wdy(2,2,jdp(i))
      coef34 = ymy1*ymy2
      term3(i) = coef34*ymy4*wdy(3,2,jdp(i))
      term4(i) = coef34*ymy3*wdy(4,2,jdp(i))
    END DO

    ! Loop over fields.

    DO m = 1, pf
      DO i = 1, pgls
        fdp(i,m) = fint(i,1,m)*term1(i) + fint(i,2,m)*term2(i) + &
                   fint(i,3,m)*term3(i) + fint(i,4,m)*term4(i)
      END DO
    END DO

  END SUBROUTINE lagyin

  SUBROUTINE limdx(pidim,ibeg,len,dx,f,qxl,qxr)

    ! Description:
    !
    ! Limit the derivative estimates for data.
    !
    ! Method:
    !
    ! Limit the derivative estimates for data on an equally spaced grid
    ! so they satisfy the SCM0 condition, that is, the spline will be
    ! monotonic, but only C0 continuous on the domain
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: dx
    INTEGER, INTENT (in) :: ibeg, len, pidim

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: f(pidim)

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: qxl(pidim), qxr(pidim)

    ! pidim   Length of f, qxl, and qxr.
    ! ibeg    First interval of grid for which derivatives are computed.
    ! len     Number of grid intervals for which derivatives are computed.
    !         (There are pidim - 1 intervals between the pidim gridpoints
    !         represented in f, qxl, and qxr.)
    ! dx      Value of grid spacing.
    ! f       Values on equally spaced grid from which derivatives qxl and
    !         qxr were computed.
    ! qxl     qxl(i) is the limited derivative at the left  edge of interval
    ! qxr     qxr(i) is the limited derivative at the right edge of interval

    !  Local scalars:
    REAL(dp) :: rdx                ! 1./dx
    INTEGER :: i               ! index
    INTEGER :: iend            ! index to end work on vector

    !  Local arrays:
    REAL(dp) :: deli(pbpts)        ! simple linear derivative

    !!LK  !  External Subroutines
    !!LK  EXTERNAL scm0


    !  Executable Statements

    IF (pidim>pbpts) THEN
      WRITE (nout,'(A)') 'LIMDX: Local work array DELI not dimensioned large enough'
      WRITE (nout,'(A,I5)') '       Increase local parameter pbpts to ',pidim
      CALL finish('limdx','Run terminated.')
    END IF

    iend = ibeg + len - 1
    rdx = 1.0_dp/dx

    DO i = ibeg, iend
      deli(i) = (f(i+1)-f(i))*rdx
    END DO

    ! Limiter

    CALL scm0(len,deli(ibeg),qxl(ibeg),qxr(ibeg))

  END SUBROUTINE limdx

  SUBROUTINE limdy(pf,fint,dy,jdp,fyb,fyt)

    ! Description:
    !
    ! Limit the y-derivative estimates.
    !
    ! Method:
    !
    ! Limit the y-derivative estimates so they satisy the SCM0 for the
    ! x-interpolated data corresponding to the departure points of a single
    ! latitude slice in the global grid, that is, they are monotonic, but
    ! spline has only C0 continuity
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pf

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: fint(plon,plev,ppdy,pf), dy(platd)
    INTEGER, INTENT (in) :: jdp(plon,plev)

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: fyb(plon,plev,pf), fyt(plon,plev,pf)

    ! pf      Number of fields being interpolated.
    ! fint    (fint(i,k,j,m),j=1,ppdy) contains the x interpolants at each
    !         latitude needed for the y derivative estimates at the
    !         endpoints of the interval that contains the departure point
    !         for grid point (i,k).  The last index of fint allows for
    !         interpolation of multiple fields.  fint is generated by a
    !         call to herxin.
    ! dy      Increment in the y-coordinate value for each interval in the
    !         extended array.
    ! jdp     jdp(i,k) is the index of the y-interval that contains the
    !         departure point corresponding to global grid point (i,k) in
    !         the latitude slice being forecasted.
    !         Suppose yb contains the y-coordinates of the extended array
    !         and ydp(i,k) is the y-coordinate of the departure point
    !         corresponding to grid point (i,k).  Then,
    !         yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
    ! fyb     fyb(i,k,.) is the limited derivative at the bot of the y
    !         interval that contains the departure point of global grid
    !         point (i,k).
    ! fyt     fyt(i,k,.) is the limited derivative at the top of the y
    !         interval that contains the departure point of global grid
    !         point (i,k).

    !  Local scalars:
    INTEGER :: i, k, m          ! indices
    INTEGER :: jb               ! index for bottom of interval
    INTEGER :: jt               ! index for top    of interval

    !  Local arrays:
    REAL(dp) :: rdy(plon,plev)      ! 1./dy
    REAL(dp) :: deli(plon,plev)     ! simple linear derivative

    !!LK  !  External subroutines
    !!LK  EXTERNAL scm0


    !  Executable statements

    jb = ppdy/2
    jt = jb + 1

    DO k = 1, plev
      DO i = 1, plon
        rdy(i,k) = 1.0_dp/dy(jdp(i,k))
      END DO
    END DO

    ! Loop over fields.

    DO m = 1, pf
      DO k = 1, plev
        DO i = 1, plon
          deli(i,k) = (fint(i,k,jt,m)-fint(i,k,jb,m))*rdy(i,k)
        END DO
      END DO

      ! Limiter

      CALL scm0(plon*plev,deli,fyb(1,1,m),fyt(1,1,m))
    END DO

  END SUBROUTINE limdy

  SUBROUTINE limdz(f,dsig,fst,fsb)

    ! Description:
    !
    ! Apply SCMO limiter to vertical derivative estimates on a vertical
    ! slice.
    !
    ! Method:
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: f(plon,plev,pcnst), dsig(plev)

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: fst(plon,plev,pcnst), fsb(plon,plev,pcnst)

    ! f       Field values used to compute the discrete differences for
    !         each interval in the vertical grid.
    ! dsig    Increment in the sigma-coordinate value for each interval.
    ! fst     Limited derivative at the top of each interval.
    ! fsb     Limited derivative at the bottom of each interval.

    !  Local scalars:
    REAL(dp) :: rdsig                ! 1./dsig
    INTEGER :: i                 ! longitude   index
    INTEGER :: k                 ! vertical    index
    INTEGER :: m                 ! constituent index

    !  Local arrays:
    REAL(dp) :: deli(plon,plevm1)    ! simple linear derivative

    !!LK  !  External subroutines
    !!LK  EXTERNAL scm0


    !  Executable statements

    ! Loop over fields.

    DO m = 1, pcnst
      DO k = 1, plevm1
        rdsig = 1.0_dp/dsig(k)
        DO i = 1, plon
          deli(i,k) = (f(i,k+1,m)-f(i,k,m))*rdsig
        END DO
      END DO
      CALL scm0(plon*plevm1,deli,fst(1,1,m),fsb(1,1,m))
    END DO

  END SUBROUTINE limdz

  SUBROUTINE s2gphi(lam,cosphi,sinphi,lamsc,phisc,phigc)

    ! Description:
    !
    ! Calculate transformed local geodesic latitude coordinates.
    !
    ! Method:
    !
    ! Calculate transformed local geodesic latitude coordinates for a set
    ! of points, each of which is associated with a grid point in a global
    ! latitude slice. Transformation is spherical to local geodesic.
    ! (Williamson and Rasch, 1991)
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: cosphi, sinphi

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: lam(plon), lamsc(plon,plev), phisc(plon,plev)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: phigc(plon,plev)

    ! lam    longitude coordinates of the global grid points in spherical
    !        system. The grid points in the global array are the reference
    !        points for the local geodesic systems.
    ! cosphi cosine of the latitude of the global latitude slice.
    ! sinphi sine of the latitude of the global latitude slice.
    ! lamsc  longitude coordinate of dep. points in spherical coordinates.
    ! phisc  latitude  coordinate of dep. points in spherical coordinates.
    ! phigc  latitude  coordinate of dep. points in local geodesic coords.

    !  Local scalars:
    REAL(dp) :: clamsc, cphisc, sphisc
    INTEGER :: i, k

    !  Intrinsic functions
    INTRINSIC ASIN, COS, SIN


    !  Executable statements

    DO k = 1, plev
      DO i = 1, plon
        sphisc = SIN(phisc(i,k))
        cphisc = COS(phisc(i,k))
        clamsc = COS(lam(i)-lamsc(i,k))
        phigc(i,k) = ASIN(sphisc*cosphi-cphisc*sinphi*clamsc)
      END DO
    END DO

  END SUBROUTINE s2gphi

  SUBROUTINE s2gvel(udp,vdp,lam,cosphi,sinphi,lamdp,phidp,upr,vpr)

    ! Description:
    !
    ! Transform velocity components at departure points.
    !
    ! Method:
    !
    ! Transform velocity components at departure points associated with a
    ! single latitude slice from spherical coordinates to local geodesic
    ! coordinates. (Williamson and Rasch, 1991)
    !
    ! Authors
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: cosphi, sinphi

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: lam(plon), lamdp(plon,plev), phidp(plon,plev),  &
         &              udp(plon,plev), vdp(plon,plev)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: upr(plon,plev), vpr(plon,plev)

    ! udp    u-component of departure point velocity in spherical coords.
    ! vdp    v-component of departure point velocity in spherical coords.
    ! lam    Longitude of arrival point position (model grid point) in
    !        spherical coordinates.
    ! cosphi Cos of latitude of arrival point positions (model grid pt).
    ! sinphi Sin of latitude of arrival point positions (model grid pt).
    ! lamdp  Longitude of departure point position in spherical
    !        coordinates.
    ! phidp  Latitude  of departure point position in spherical
    !        coordinates.
    ! upr    u-component of departure point velocity in geodesic coords.
    ! vpr    v-component of departure point velocity in geodesic coords.

    !  Local scalars:
    REAL(dp) :: cdlam, clamp, cphid, cphip, dlam, sdlam, slamp, sphid, sphip
    INTEGER :: i, k

    !  Intrinsic functions
    INTRINSIC ASIN, COS, SIN


    !  Executable statements

    DO k = 1, plev
      DO i = 1, plon
        dlam = lam(i) - lamdp(i,k)
        sdlam = SIN(dlam)
        cdlam = COS(dlam)
        sphid = SIN(phidp(i,k))
        cphid = COS(phidp(i,k))
        sphip = sphid*cosphi - cphid*sinphi*cdlam
        cphip = COS(ASIN(sphip))
        slamp = -sdlam*cphid/cphip
        clamp = COS(ASIN(slamp))

        vpr(i,k) = (vdp(i,k)*(cphid*cosphi+sphid*sinphi*cdlam) &
             -udp(i,k)*sinphi*sdlam)/cphip
        upr(i,k) = (udp(i,k)*cdlam+vdp(i,k)*sphid*sdlam &
             +vpr(i,k)*slamp*sphip)/clamp
      END DO
    END DO

  END SUBROUTINE s2gvel

  SUBROUTINE scm0(n,deli,df1,df2)

    ! Description:
    !
    ! Apply SCM0 limiter to derivative estimates.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! Modified:          U. Schlese,  DKRZ - Hamburg,  May 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !

    !  Local parameters:
    ! machine precision parameter to prevent rounding
    ! negatives in interpolant
    REAL(dp), PARAMETER :: eps = 1.E-12_dp
    ! factor applied in limiter
    REAL(dp), PARAMETER :: fac = 3.0_dp*(1.0_dp-eps)

    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: n

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: deli(n)

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: df1(n), df2(n)

    ! n      Dimension of input arrays.
    ! deli   deli(i) is the discrete derivative on interval i, i.e.,
    ! deli(i) = ( f(i+1) - f(i) )/( x(i+1) - x(i) ).
    ! df1    df1(i) is the limited derivative at the left  edge of interval
    ! df2    df2(i) is the limited derivative at the right edge of interval

    !  Local scalars:
    REAL(dp) :: tmp1  ! derivative factor
    REAL(dp) :: tmp2  ! abs(tmp1)
    INTEGER :: i

    !  Intrinsic functions
    INTRINSIC ABS, MERGE

    !  Executable statements

    DO i = 1, n
      tmp1 = fac*deli(i)
      tmp2 = ABS(tmp1)

      df1(i) = MERGE(0.0_dp,df1(i),-(deli(i)*df1(i)) >= 0.0_dp)
      df2(i) = MERGE(0.0_dp,df2(i),-(deli(i)*df2(i)) >= 0.0_dp)
      df1(i) = MERGE(tmp1,df1(i),tmp2-ABS(df1(i)) < 0.0_dp)
      df2(i) = MERGE(tmp1,df2(i),tmp2-ABS(df2(i)) < 0.0_dp)

    END DO

  END SUBROUTINE scm0

#ifndef SLDIAG
  SUBROUTINE semi_lagrangian_transport(dt,ra,iterdp,wb, &
       fbout,hw1lat,hw2lat,hw3lat)
#else
  SUBROUTINE semi_lagrangian_transport (dt,ra,iterdp,wb,&
       fbout,hw1lat,hw2lat,hw3lat,diag)
#endif

    ! Description:
    !
    ! Drive the slt algorithm on a given latitude slice.
    !
    ! Method:
    !
    ! Drive the slt algorithm on a given latitude slice in the extended
    ! data arrays using information from the entire latitudinal extent
    ! of the arrays.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    ! parallel version:  T. Diehl, DKRZ, July 1999
    !
    ! for more details see file AUTHORS
    !
    
    !  Local parameters:
    REAL(dp), PARAMETER :: phigs = 1.221730_dp  ! cut-off latitude: about 70 degree
    
    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: dt, ra
    INTEGER, INTENT (in) :: iterdp
    
    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: wb(plon,plevp1,plat)

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: fbout(plond,plev,pcnst,platd,2)
#ifdef SLDIAG
    REAL(dp), INTENT (out) :: diag(plond,plev,pcnst)
#endif
    REAL(dp), INTENT (out) :: hw1lat(pcnst,plat), hw2lat(pcnst,plat), &
         hw3lat(pcnst,plat)

    ! pmap    Dimension of kdpmpX arrays
    ! jcen    Latitude index in extended grid corresponding to lat slice
    !         being forecasted.
    ! jgc     Latitude index in model    grid corresponding to lat slice
    !         being forecasted.
    ! dt      Time interval that parameterizes the parcel trajectory.
    ! ra      Reciprocal of radius of earth.
    ! iterdp  Number of iterations used for departure point calculation.
    ! ub      u-velocity component.
    ! uxl     uxl ( , , 1, ) = x-derivatives of u at the left  (west) edge
    !         of given interval
    ! uxl ( , , 2, ) = x-derivatives of v at the left  (west) edge
    !         of given interval
    ! uxr     uxr ( , , 1, ) = x-derivatives of u at the right (east) edge
    !         of given interval
    ! uxr ( , , 2, ) = x-derivatives of v at the right (east) edge
    !         of given interval
    ! wb      z-velocity component (eta-dot).
    ! fb      Scalar components to be transported.
    ! qxl     x-derivatives at the left  edge of each interval containing
    !         the departure point.
    ! qxr     x-derivatives at the right edge of each interval containing
    !         the departure point.
    ! lam     Longitude values for the extended grid.
    ! phib    Latitude  values for the extended grid.
    ! dphib   Interval between latitudes in the extended grid.
    ! sig     Hybrid eta values at the "full-index" levels.
    ! sigh    Half-index eta-levels including sigh(i,1) = eta(1/2) = 0.0
    !         and sigh(i,plev+1) = eta(plev+1/2) = 1.  Note that in general
    !         sigh(i,k) .lt. sig(i,k)  where sig(i,k) is the hybrid value
    !         at the k_th full-index level.
    ! dsig    Interval lengths in full-index hybrid level grid.
    ! dsigh   Interval lengths in half-index hybrid level grid.
    ! lbasdy  Weights for Lagrange cubic derivative estimates on the
    !         unequally spaced latitude grid.
    ! lbasdz  Weights for Lagrange cubic derivative estimates on the
    !         unequally spaced vertical grid (full levels).
    ! lbassd  Weights for Lagrange cubic derivative estimates on the
    !         unequally spaced vertical grid (half levels).
    ! lbasiy  Weights for Lagrange cubic interpolation on the unequally
    !         spaced latitude grid.
    ! kdpmpf  indices of artificial grid mapped into the full level grid
    ! kdpmph  indices of artificial grid mapped into the half level grid
    ! lammp   Longitude coordinates of the trajectory mid-points of the
    !         parcels that correspond to the global grid points contained
    !         in the latitude slice being forecasted.  On entry lammp
    !         is an initial guess.
    ! phimp   Latitude coordinates of the trajectory mid-points of the
    !         parcels that correspond to the global grid points contained
    !         in the latitude slice being forecasted.  On entry phimp
    !         is an initial guess.
    ! sigmp   Hybrid value at the trajectory midpoint for each gridpoint
    !         in a vertical slice from the global grid.  On entry sigmp is
    !         an initial guess.
    ! fbout   Extended array only one latitude of which, however, is filled
    !         with forecasted (transported) values.  This routine must be
    !         called multiple times to fill the entire array.  This is
    !         done to facilitate multi-tasking.
    ! diag    Time tendency due to horizontal advection for each scalar
    !         component at the global grid points that correspond to the
    !         latitude slice being forecasted.  As in "fbout", only one
    !         latitude slice at a time is filled.
    
    !  Local scalars:
    INTEGER :: m
    INTEGER :: i, j, k, ih, ilat, jlat, jcen
    INTEGER :: isafe, nxpt_max, iext ! for parallel version
    LOGICAL :: locgeo  ! flag indicating coordinate sys
    
    REAL(dp) :: fhr(plon,plev,pcnst)        ! horizontal interpolants
    REAL(dp) :: fhsb(plon,plev,pcnst)       ! derivative at bot of interval
    REAL(dp) :: fhst(plon,plev,pcnst)       ! derivative at top of interval
    REAL(dp) :: fint(plon,plev,ppdy,pcnst)  ! work space
    REAL(dp) :: fyb(plon,plev,pcnst)        ! work space
    REAL(dp) :: fyt(plon,plev,pcnst)        ! work space
    REAL(dp) :: lamdp(plon,plev)            ! zonal      departure pt. coord.
    REAL(dp) :: phidp(plon,plev)            ! meridional departure pt. coord.
    REAL(dp) :: sigdp(plon,plev)            ! vertical   departure pt. coord.
    REAL(dp) :: wsb(plon,plevp1)             ! w derivative at bot of interval
    REAL(dp) :: wst(plon,plevp1)             ! w derivative at top of interval
    
    REAL(dp) :: fdp(plon,plev,pcnst)
    
#ifdef SLDIAG
    REAL(dp) :: rdt                         ! 1./dt
#endif
    INTEGER :: idp(plon,plev)           ! zonal      dep point index
    INTEGER :: jdp(plon,plev)           ! meridional dep point index
    INTEGER :: kdp(plon,plev)           ! vertical   dep point index

    REAL(dp) :: uvb(plond,plev,2,platd,2)

    !  Intrinsic functions
    INTRINSIC ABS

    !  Executable statements
    
    hw1lat(:,:) = 0.0_dp
    hw2lat(:,:) = 0.0_dp
    hw3lat(:,:) = 0.0_dp

    DO ih = 1,2
      DO j = 1,platd
        DO k = 1,plev
          DO i = 1,plond
            uvb(i,k,1,j,ih) = ub(i,k,j,ih)
            uvb(i,k,2,j,ih) = vb(i,k,j,ih)
          END DO
        END DO
      END DO
    END DO

    DO ilat = 1, plat  ! continuous local latitude index for SLT (S->N)
      jlat = plat+1-ilat  ! latitude index for wfld (wfld structure: N->S)
      IF (ilat <= plato2) THEN
        ih = 1
        jcen = j1 - 1 + ilat  ! latitude index for S
      ELSE
        ih = 2
        jcen = j1 - 1 + ilat -plato2  ! latitude index for N
      END IF

      ! Calculate zonal mass of constituents before advected by slt.

      CALL zonal_mass_integral(cwava,gw(ilat),fb(1,1,1,jcen,ih), &
           pdel(1,1,ilat), hw1lat(:,jlat),jlat)


      ! Horizontal interpolation
    
      ! Compute departure points and corresponding indices.
      ! Poleward of latitude phigs (radians), perform the computation in
      ! local geodesic coordinates.
      ! Equatorward of latitude phigs, perform the computation in global
      ! spherical coordinates
      
      locgeo = ABS(phi(jcen,ih)) >= phigs
      
      CALL sphdep(ih,jcen,ilat,dt,ra,iterdp,locgeo,uvb(1,1,1,1,ih), &
           uxl(1,1,1,1,ih),uxr(1,1,1,1,ih),&
           lam,phi(1,ih),lbasiy(1,1,1,ih), &
           lammp(:,:,ilat),phimp(:,:,ilat),lamdp,phidp,idp,jdp, &
           nxpt_a(1,jcen,ih))
      
      ! Interpolate scalar fields to the departure points.
      
      CALL hrintp(pcnst,fb(1,1,1,1,ih), &
           qxl(1,1,1,1,ih),qxr(1,1,1,1,ih),&
           lam,dphi(1,ih),lbasdy(1,1,1,ih),lamdp,phidp,idp, &
           jdp,plimdr,fint,fyb,fyt,fhr)

#ifdef SLDIAG

      ! Compute time tendency due to horizontal advection.

      rdt = 1.0_dp/dt
      DO m = 1, pcnst
        DO k = 1, plev
          DO i = 1, plon
            diag(i,k,m) = (fhr(i,k,m)-fb(nxpt+i,k,m,jcen,ih))*rdt
          END DO
        END DO
      END DO
#endif

      ! Vertical interpolation.
      ! Compute vertical derivatives of vertical wind

      CALL cubzdr(plon,plevp1,wb(1,1,jlat),lbassd,wst,wsb)

      ! Compute departure points and corresponding indices.

      CALL vrtdep(dt,iterdp,wb(1,1,jlat),wst,wsb,etamid,etaint,detai,&
           kdpmpf,kdpmph,sigmp(:,:,ilat),sigdp,kdp)
      
      ! Vertical derivatives of scalar fields.
      ! Loop over constituents.
      
      DO m = 1, pcnst
        CALL cubzdr(plon,plev,fhr(1,1,m),lbasdz,fhst(1,1,m),fhsb(1,1,m))
      END DO
      IF (plimdr) THEN
        CALL limdz(fhr,detam,fhst,fhsb)
      END IF
      
      ! Vertical interpolation of scalar fields.

      CALL herzin(plev,pcnst,fhr,fhst,fhsb,etamid,detam,sigdp,kdp,fdp)

      ! Transfer transported values to extended array

      DO m = 1,pcnst
        DO k = 1,plev
          DO i = 1,plon
            fbout(nxpt+i,k,m, jcen,ih) = fdp(i,k,m)
          END DO
        END DO
      END DO
      
      ! Locally update the nxpt_p array using idp
      ! The array values over the row of processors in longitude
      ! are synchronized in sltini
      
      isafe = 4
      DO k = 1,plev
        nxpt_max = 0
        DO i = 1,plon
          iext= i + nxpt
          nxpt_max = MAX(nxpt_max,ABS(iext - idp(i,k)) + isafe)
        END DO
        nxpt_a(k,jcen,ih) = MIN(nxpt_max,nxpt)
      END DO

      ! Calculate zonal mass of constituents after advected by slt.

      ! Compute contribution to global integral hw3.
      
      CALL zonal_mass_integral(cwava,gw(ilat),fbout(1,1,1,jcen,ih), &
           pdel(1,1,ilat), hw2lat(:,jlat),jlat)
      
      ! Compute contribution to global integral hw3.
      
      CALL global_mass_integral(cwava,gw(ilat),fb(1,1,1,jcen,ih), &
           fbout(1,1,1,jcen,ih), pdel(1,1,ilat),etamid,kftype, &
           hw3lat(:,jlat),jlat)
    END DO

  END SUBROUTINE semi_lagrangian_transport

  SUBROUTINE extend_semi_lagrangian

    ! Description:
    !
    !...Prepare the extended arrays for use in the SLT routines
    !
    !   1)  Fill latitude extensions.
    !   2)  Fill longitude extensions.
    !   3)  Compute x-derivatives
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! Parallel version:  T. Diehl, DKRZ, July 1999
    !
    ! for more details see file AUTHORS
    !
       
    INTEGER       :: itmp(plev,platd,2)
    INTEGER       :: nxpt_a_g(plev,platd,2,lc%nprocb)
    INTEGER       :: puvpts, pqpts, ih, j, k, jjfirst, jjlast, nei_nxpt
    INTEGER       :: tag_base, tag, count, dest, source, tag2, count2
    INTEGER       :: tagcount, tagtable(lc%nprocb), jproc, isource, src_ln
    LOGICAL, SAVE :: first = .TRUE.
    
    !!LK  ! external subroutines
    !!LK  EXTERNAL extyv, extys, extx, cubxdr, limdx
    
    ! executable statements

    qxl(:,:,:,:,:) = 0.0_dp
    qxr(:,:,:,:,:) = 0.0_dp
    uxl(:,:,:,:,:) = 0.0_dp
    uxr(:,:,:,:,:) = 0.0_dp

    puvpts = plond*plev
    pqpts = plond*plev*pcnst
    
    ! Fill latitude extensions beyond the southern- and northern-most
    ! latitudes in the global grid
    
    CALL extyv(plev,coslam,sinlam,ub)
    CALL extyv(plev,coslam,sinlam,vb)
    
    ! loop over fields now in extys
    
    CALL extys(pcnst,plev,fb)

    !=============================
    !either call p_barrier or use different message tags in following section
    
    CALL p_barrier (p_communicator_d)
    
    ! Begin of parallel extensions
    !======================================================================
    ! Synchronize the nxpt_a array for this timestep
    ! (The local update takes place in sltb1)
    
    IF (first .OR. lc%nprocb == 1) THEN
      !  Initialize the nxpt_a array
      nxpt_a(:,:,:) = nxpt
      first = .FALSE.
    ELSE
      !     IF (.not. debug_seriell .or. lc%nprocb > 1) THEN
      
      !   Be sure the overlap extension is enough to support interpolation
      
      jjfirst = nxpt + jintmx + 1
      DO ih = 1,2
        DO j = 1,nxpt+jintmx
          DO k = 1,plev
            itmp(k,j,ih) = nxpt_a(k,jjfirst,ih)
          END DO
        END DO
      END DO
      jjlast = plato2 + nxpt + jintmx
      DO ih = 1,2
        DO j = jjlast+1,platd
          DO k = 1,plev
            itmp(k,j,ih) = nxpt_a(k,jjlast,ih)
          END DO
        END DO
      END DO
      
      DO ih = 1,2
        DO j = jjfirst,jjlast
          DO k = 1,plev
            nei_nxpt = 0
            IF (k > 1) nei_nxpt = MAX(nxpt_a(k-1,j,ih),nei_nxpt)
            nei_nxpt = MAX(nxpt_a(k,j,ih),nei_nxpt)
            IF (k <= plev-1) nei_nxpt = MAX(nxpt_a(k+1,j,ih),nei_nxpt)
            IF (k <= plev-2) nei_nxpt = MAX(nxpt_a(k+2,j,ih),nei_nxpt)
            nei_nxpt = MAX(nxpt_a(k,j-1,ih),nei_nxpt)
            nei_nxpt = MAX(nxpt_a(k,j,ih),nei_nxpt)
            nei_nxpt = MAX(nxpt_a(k,j+1,ih),nei_nxpt)
            nei_nxpt = MAX(nxpt_a(k,j+2,ih),nei_nxpt)
            itmp(k,j,ih) = nei_nxpt
          END DO
        END DO
      END DO
      DO ih = 1,2
        DO j = 1,platd
          DO k = 1,plev
            nxpt_a(k,j,ih) = itmp(k,j,ih)
          END DO
        END DO
      END DO
      
      ! Synchronize the nxpt_a across the processors in the row
      
      ! On entry nxpt_a contains local data; on exit it contains the
      ! vector max of the processor row
      
      tag_base = 400 + 10 + lc%nprocb  ! last change occurred in extys
      IF (lc%set_b /= 1) THEN
        tag = tag_base + 1
        count = SIZE(itmp)
        dest = gc(lc%mapmesh(1,lc%set_a))%pe
        ! send to column 1
        CALL p_isend(itmp,dest,tag,p_count=count)
        ! receive updated values from column 1
        source = gc(lc%mapmesh(1,lc%set_a))%pe
        tag2 = tag_base + 2
        count2 = SIZE(nxpt_a)
        CALL p_recv(nxpt_a,source,tag2,p_count=count2)
      ELSE
        ! collect
        nxpt_a_g(:,:,:,1) = itmp(:,:,:)
        tagcount = 1
        tagtable(1) = tag_base + 1
        DO jproc = 1,lc%nprocb - 1
          CALL p_probe(tagcount,tagtable,source,tag,count)
          isource = indx(source,gc)
          src_ln = gc(isource)%set_b
          CALL p_recv(nxpt_a_g(1,1,1,src_ln),source,tag,p_count=count)
        END DO
        
        nxpt_a = MAXVAL(nxpt_a_g, DIM=4)
        
        ! distribute
        
        DO jproc = 1,lc%nprocb
          IF (gc(lc%mapmesh(jproc,lc%set_a))%pe /= lc%pe) THEN
            dest = gc(lc%mapmesh(jproc,lc%set_a))%pe
            tag = tag_base + 2
            count = SIZE(nxpt_a)
            CALL p_send(nxpt_a,dest,tag,p_count=count)
          END IF
        END DO
      END IF
      
      !     END IF
    END IF


    CALL p_wait
    
    !  Update the overlap between processors
    !  The arrays have been extended as if poles were present
    
    CALL overlap(ub,vb,fb)
    
    ! end of parallel extensions
    !====================================================================
    
    ! Fill longitude extensions
    
    CALL extx(1,plev,ub(1,1,1,1))
    CALL extx(1,plev,ub(1,1,1,2))
    CALL extx(1,plev,vb(1,1,1,1))
    CALL extx(1,plev,vb(1,1,1,2))
    
    ! loop over fields now in extx
    
    CALL extx(pcnst,plev,fb(1,1,1,1,1))
    CALL extx(pcnst,plev,fb(1,1,1,1,2))
    
    ! Compute x-derivatives
    
    DO ih = 1,2
      DO j = 1, platd
        CALL cubxdr(puvpts,2,puvpts-3,dlam,ub(1,1,j,ih), &
             uxl(1,1,1,j,ih),uxr(1,1,1,j,ih))
        CALL cubxdr(puvpts,2,puvpts-3,dlam,vb(1,1,j,ih), &
             uxl(1,1,2,j,ih),uxr(1,1,2,j,ih))
        CALL cubxdr(pqpts,2,pqpts-3,dlam,fb(1,1,1,j,ih), &
             qxl(1,1,1,j,ih),qxr(1,1,1,j,ih))
        IF (plimdr) THEN
          CALL limdx(pqpts,2,pqpts-3,dlam,fb(1,1,1,j,ih), &
               qxl(1,1,1,j,ih),qxr(1,1,1,j,ih))
        END IF
      END DO
    END DO

    CALL p_barrier (p_communicator_d)

  END SUBROUTINE extend_semi_lagrangian

  SUBROUTINE sphdep(ih,jcen,jgc,dt,ra,iterdp,locgeo,uvb,uxl,uxr,lam,phib, &
       lbasiy,lammp,phimp,lamdp,phidp,idp,jdp,nxpt_a)

    ! Description:
    !
    ! Compute departure points for semi-Lagrangian transport on surface of sphere
    !
    ! Method:
    !
    ! Compute departure points for semi-Lagrangian transport on surface of
    ! sphere using midpoint quadrature.  Computations are done in:
    !
    !   1) "local geodesic"   coordinates for "locgeo" = .true.
    !   2) "global spherical" coordinates for "locgeo" = .false.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !
    
    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: dt, ra
    INTEGER, INTENT (in) :: iterdp, ih, jcen, jgc
    LOGICAL, INTENT (in) :: locgeo
    
    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: lam(plond), lbasiy(4,2,platd), phib(platd), &
         uvb(plond,plev,2,platd), uxl(plond,plev,2,platd), &
         uxr(plond,plev,2,platd)
    INTEGER, INTENT (in) :: nxpt_a(plev)

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: lammp(plon,plev), phimp(plon,plev)
    
    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: lamdp(plon,plev), phidp(plon,plev)
    INTEGER, INTENT (out) :: idp(plon,plev), jdp(plon,plev)
    
    
    ! jcen    Index in extended grid corresponding to latitude being
    !         forecast.
    ! jgc     Index in model    grid corresponding to latitude being
    !         forecast.
    ! dt      Time interval that parameterizes the parcel trajectory.
    ! ra      Reciprocal of radius of earth.
    ! iterdp  Number of iterations used for departure point calculation.
    ! locgeo  Logical flag to indicate computation in "local geodesic" or
    !         "global spherical" space.
    ! ub      Longitudinal velocity components in spherical coordinates.
    ! uxl     uxl ( , , 1, ) = x-derivatives of u at the left  (west) edge
    !         of given interval
    !         uxl ( , , 2, ) = x-derivatives of v at the left  (west) edge
    !         of given interval
    ! uxr     uxr ( , , 1, ) = x-derivatives of u at the right (east) edge
    !         of given interval
    !         uxr ( , , 2, ) = x-derivatives of v at the right (east) edge
    !         of given interval
    ! lam     Longitude values for the extended grid.
    ! phib    Latitude  values for the extended grid.
    ! lbasiy  Weights for Lagrange cubic interpolation on the unequally
    !         spaced latitude grid.
    ! lammp   Longitude coordinates of the trajectory mid-points of the
    !         parcels that correspond to the global grid points contained
    !         in the latitude slice being forecast.  On entry lammp
    !         is an initial guess.
    ! phimp   Latitude coordinates of the trajectory mid-points of the
    !         parcels that correspond to the global grid points contained
    !         in the latitude slice being forecast.  On entry phimp
    !         is an initial guess.
    ! lamdp   Longitude coordinates of the departure points that correspond
    !         to the global grid points contained in the latitude slice
    !         being forecast.  lamdp is constrained so that
    !         0.0 .le. lamdp(i) .lt. 2*pi .
    ! phidp   Latitude coordinates of the departure points that correspond
    !         to the global grid points contained in the latitude slice
    !         being forecast.  If phidp is computed outside the latitudinal
    !         domain of the extended grid, then an abort will be called by
    !         subroutine "trjgl".
    ! idp     Longitude index of departure points.  This index points into
    !         the extended arrays, e.g.,
    !         lam (idp(i,k)) .le. lamdp(i,k) .lt. lam (idp(i,k)+1).
    ! jdp     Latitude  index of departure points.  This index points into
    !         the extended arrays, e.g.,
    !         phib(jdp(i,k)) .le. phidp(i,k) .lt. phib(jdp(i,k)+1).
    
    !  Local scalars:
    REAL(dp) :: cphic       ! cos(phicen)
    REAL(dp) :: dlam        ! increment of grid in x-direction
    REAL(dp) :: dttmp       ! time step (seconds)
    REAL(dp) :: finc        ! time step factor
    REAL(dp) :: phicen      ! latitude coord of current lat slice
    REAL(dp) :: sphic       ! sin(phicen)
    REAL(dp) :: pi
    INTEGER :: i, iter, k
    
    !  Local arrays:
    REAL(dp) :: fint(plon,plev,ppdy,2) ! u/v x-interpolants
    REAL(dp) :: lampr(plon,plev)       ! relative long coord of dep pt
    REAL(dp) :: phipr(plon,plev)       ! relative lat  coord of dep pt
    REAL(dp) :: upr(plon,plev)         ! u in local geodesic coords
    REAL(dp) :: uvmp(plon,plev,2)      ! u/v (spherical) interpltd to dep pt
    REAL(dp) :: vpr(plon,plev)         ! v in local geodesic coords
    
    !!LK  !  External subroutines
    !!LK  EXTERNAL bandij, g2spos, hrvint, s2gphi, s2gvel, trajmp, trjgl, trjmps
    
    !  Intrinsic functions
    INTRINSIC COS, SIN, ATAN
    
    
    !  Executable statements
    
    !  dlam = lam(nxpt+2) - lam(nxpt+1)
    ! for parallel version
    pi = 4.0_dp*ATAN(1.0_dp)
    dlam = 2.0_dp*pi/REAL(lc%nlon,dp)
    
    phicen = phib(jcen)
    cphic = COS(phicen)
    sphic = SIN(phicen)
    
    ! Convert latitude coordinates of trajectory midpoints from spherical
    ! to local geodesic basis.
    
    IF (locgeo) CALL s2gphi(lam(i1),cphic,sphic,lammp,phimp,phipr)
    
    ! Loop over departure point iterates.
    
    DO iter = 1, iterdp
      
      ! Compute midpoint indices.
      
      CALL bandij(dlam,phib,lammp,phimp,idp,jdp,nxpt_a)
      
      ! Interpolate velocity fields to midpoint locations (tensor product
      ! Lagrange cubic interpolation).

      CALL hrvint(2,uvb,uxl,uxr,lbasiy,lammp,phimp,idp,jdp,fint,uvmp)
      
      ! Put u/v on unit sphere
      
      DO k = 1, plev
        DO i = 1, plon
          uvmp(i,k,1) = uvmp(i,k,1)*ra
          uvmp(i,k,2) = uvmp(i,k,2)*ra
        END DO
      END DO
      
      ! For local geodesic:
      
      !   a) Convert velocity coordinates at trajectory midpoints from
      !      spherical coordinates to local geodesic coordinates,
      !   b) Estimate midpoint parcel trajectory,
      !   c) Convert back to spherical coordinates
      
      ! Else, for global spherical
      
      !      Estimate midpoint trajectory with no conversions
      
      IF (locgeo) THEN
        CALL s2gvel(uvmp(1,1,1),uvmp(1,1,2),lam(i1),cphic,sphic,lammp,phimp, &
             upr,vpr)
        CALL trajmp(dt,upr,vpr,phipr,lampr)
        dttmp = 0.5_dp*dt
        
        CALL g2spos(ih,dttmp,lam(i1),phib,phicen,cphic,sphic,upr,vpr, &
             lampr,phipr, lammp,phimp)
        
      ELSE
        
        CALL trjmps(dt,uvmp(1,1,1),uvmp(1,1,2),phimp,lampr,phipr)
        finc = 1.0_dp
        
        CALL trjgl(jgc,finc,phicen,lam(i1),lampr,phipr,lammp,phimp)
      END IF
      
    END DO ! End of iter=1,iterdp loop
    
    ! Compute departure points in geodesic coordinates, and convert back
    ! to spherical coordinates.
    
    ! Else, compute departure points directly in spherical coordinates
    
    IF (locgeo) THEN
      DO k = 1, plev
        DO i = 1, plon
          lampr(i,k) = 2.0_dp*lampr(i,k)
          phipr(i,k) = 2.0_dp*phipr(i,k)
        END DO
      END DO
      dttmp = dt
      
      CALL g2spos(ih,dttmp,lam(i1),phib,phicen,cphic,sphic,upr,vpr, &
                          lampr,phipr,lamdp,phidp)
    ELSE
      finc = 2.0_dp
      CALL trjgl(jgc,finc,phicen,lam(i1),lampr,phipr,lamdp,phidp)
    END IF
    
    ! Compute departure point indicies.
    
    CALL bandij(dlam,phib,lamdp,phidp,idp,jdp,nxpt_a)
    
  END SUBROUTINE sphdep

  SUBROUTINE trajmp(dt,upr,vpr,phipr,lampr)
    
    ! Description:
    !
    ! Estimate mid-point of parcel trajectory (geodesic coordinates) based
    ! upon horizontal wind field.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !
    
    !  Scalar arguments with intent(In):
    REAL(dp) :: dt
    
    !  Array arguments with intent(In):
    REAL(dp) ::upr(pgls), vpr(pgls)
    
    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) ::  phipr(pgls)
    
    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: lampr(pgls)
    
    ! dt      Time interval that corresponds to the parcel trajectory.
    ! upr     u-coordinate of velocity corresponding to the most recent
    !         estimate of the trajectory mid-point (in geodesic system).
    ! vpr     v-coordinate of velocity corresponding to the most recent
    !         estimate of the trajectory mid-point (in geodesic system).
    ! phipr   Phi value at trajectory mid-point (geodesic coordinates).
    !         On entry this is the most recent estimate.
    ! lampr   Lambda value at trajectory mid-point (geodesic coordinates).
    
    !  Local scalars:
    INTEGER :: i
    
    !  Intrinsic functions
    INTRINSIC COS
    
    
    !  Executable statements
    
    DO i = 1, pgls
      lampr(i) = -0.5_dp*dt*upr(i)/COS(phipr(i))
      phipr(i) = -0.5_dp*dt*vpr(i)
    END DO
    
  END SUBROUTINE trajmp

  SUBROUTINE trjgl(jgc,finc,phicen,lam,lampr,phipr,lamp,phip)
    
    ! Description:
    !
    ! Map relative trajectory mid/departure point coordinates to global
    ! latitude/longitude coordinates and test limits
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !
    
    !  Scalar arguments with intent(in):
    REAL(dp), INTENT (in) :: finc, phicen
    INTEGER, INTENT (in) :: jgc
    
    !  Array arguments with intent(in):
    REAL(dp), INTENT (in) :: lam(plon), lampr(plon,plev), phipr(plon,plev)
    
    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: lamp(plon,plev), phip(plon,plev)
    
    
    ! jgc     Index in model grid corresponding to latitude being forecast.
    ! finc    Time step factor (1. for midpoint, 2. for dep. point)
    ! phicen  Latitude value for current latitude being forecast.
    ! lam     Longitude values for the extended grid.
    ! lampr   Longitude coordinates (relative to the arrival point) of the
    !         trajectory mid-points of the parcels that correspond to the
    ! global grid points contained in the latitude slice being
    !         forecast.
    ! phipr   Latitude coordinates (relative to the arrival point) of the
    !         trajectory mid-points of the parcels that correspond to the
    !         global grid points contained in the latitude slice being
    !         forecast.
    ! lamp    Longitude coordinates of the trajectory mid-points of the
    !         parcels that correspond to the global grid points contained
    !         in the latitude slice being forecast.
    ! phip    Latitude  coordinates of the trajectory mid-points of the
    !         parcels that correspond to the global grid points contained
    !         in the latitude slice being forecast.

    !  Local scalars:
    REAL(dp) :: pi, twopi
    INTEGER :: i, k
    
    !  Local arrays:
    INTEGER :: imax(2), imin(2)
    
    !  Intrinsic functions
    INTRINSIC ATAN
    
    
    !  Executable statements
    
    pi = 4.0_dp*ATAN(1.0_dp)
    twopi = pi*2.0_dp
    
    DO k = 1, plev
      DO i = 1, plon
        lamp(i,k) = lam(i) + finc*lampr(i,k)
        phip(i,k) = phicen + finc*phipr(i,k)
        
        IF (debug_seriell .AND. lc%nprocb == 1) THEN
          IF (lamp(i,k) >= twopi) lamp(i,k) = lamp(i,k) - twopi
          IF (lamp(i,k) < 0.0_dp) lamp(i,k) = lamp(i,k) + twopi
        END IF
        
      END DO
    END DO
      
    ! Test that the latitudinal extent of trajectory is NOT over the poles

    imax(:) = MAXLOC(phip)
    imin(:) = MINLOC(phip)
    ! Since trjgl is called only if ABS(lat) < 70(degrees):
    ! Test that the latitudinal extent of trajectory is NOT over the poles
    !
    IF (phip(imax(1),imax(2)) >= pi/2.0_dp)  THEN
      CALL errmsg
    ELSE IF (phip(imin(1),imin(2)) <= -pi/2.0_dp) THEN
      CALL errmsg
    END IF
    
    RETURN
    
  CONTAINS
    
    SUBROUTINE errmsg
      
      WRITE (nerr,'(/,A)') ' Model is blowing up ...'
      WRITE (nerr,'(/,A,I5,A,I5,A,I5,A)') &
           ' Parcel associated with longitude ',imax(1),', level ',imax(2), &
           ' and latitude ',jgc,' is outside the model domain.'
      CALL finish('trjgl','Run terminated.')

    END SUBROUTINE errmsg

  END SUBROUTINE trjgl

  SUBROUTINE trjmps(dt,upr,vpr,phimp,lampr,phipr)

    ! Description:
    !
    ! Estimate mid-point interval of parcel trajectory (global spherical
    ! coordinates).
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !
    
    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: dt
    
    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: phimp(pgls), upr(pgls), vpr(pgls)
    
    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: lampr(pgls), phipr(pgls)
    
    ! dt      Time interval that corresponds to the parcel trajectory.
    ! upr     u-coordinate of velocity corresponding to the most recent
    !         estimate of the trajectory mid-point.
    ! vpr     v-coordinate of velocity corresponding to the most recent
    !         estimate of the trajectory mid-point.
    ! phimp   Phi value of trajectory midpoint (most recent estimate).
    ! lampr   Longitude coordinate of trajectory mid-point relative to the
    !         arrival point.
    ! phipr   Latitude  coordinate of trajectory mid-point relative to the
    !         arrival point.
    
    !  Local scalars:
    INTEGER :: i
    
    !  Intrinsic functions
    INTRINSIC COS
    
    
    !  Executable statements
    
    DO i = 1, pgls
      lampr(i) = -0.5_dp*dt*upr(i)/COS(phimp(i))
      phipr(i) = -0.5_dp*dt*vpr(i)
    END DO
      
  END SUBROUTINE trjmps

  SUBROUTINE vdplim(pkdim,sig,sigdp)
    
    ! Description:
    !
    ! Restrict vertical departure points to be between the top and bottom
    ! sigma levels of the "full-" or "half-" level grid.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !
    
    !  Local parameters:
    ! machine precision to restrict dep
    ! point to be inside last grid point
    REAL(dp), PARAMETER :: eps = 1.E-12_dp
    
    !  Scalar arguments with intent(In):
    INTEGER, INTENT (in) :: pkdim
    
    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: sig(pkdim)
    
    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: sigdp(pgls)
    
    ! pkdim   Vertical dimension of "sig"
    ! sig     Sigma values at the "full" or "half" model levels
    ! sigdp   Sigma value at the trajectory endpoint or midpoint for each
    !         gridpoint in a vertical slice from the global grid.  This
    !         routine restricts those departure points to within the
    !         model's vertical grid.
    
    !  Local scalars:
    INTEGER :: i
    
    
    !  Executable statements
    
    DO i = 1, pgls
      IF (sigdp(i) < sig(1)) THEN
        sigdp(i) = sig(1)
      END IF
      IF (sigdp(i) >= sig(pkdim)) THEN
        sigdp(i) = sig(pkdim) - eps
      END IF
    END DO
    
  END SUBROUTINE vdplim

  SUBROUTINE vrtdep(dt,iterdp,wb,wst,wsb,sig,sigh,dsigh,kdpmpf,kdpmph, &
       sigmp,sigdp,kdp)

    ! Description:
    !
    ! Compute vertical departure point and departure point index.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          D. Williamson, P. Rasch, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    !
    ! for more details see file AUTHORS
    !
    
    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: dt
    INTEGER, INTENT (in) :: iterdp
    
    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: dsigh(plevp1), sig(plev), sigh(plevp1), &
         wb(plon,plevp1), wsb(plon,plevp1), wst(plon,plevp1)
    INTEGER, INTENT (in) :: kdpmpf(pmap), kdpmph(pmap)
    
    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: sigmp(plon,plev)
    
    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: sigdp(plon,plev)
    INTEGER, INTENT (out) :: kdp(plon,plev)
    
    ! pmap    Dimension of kdpmap arrays
    ! dt      Time interval that parameterizes the parcel trajectory.
    ! iterdp  Number of iterations used for departure point calculation.
    ! wb      Vertical velocity component (sigma dot).
    ! wst     z-derivs at the top edge of each interval contained in wb
    ! wsb     z-derivs at the bot edge of each interval contained in wb
    ! sig     Sigma values at the full-index levels.
    ! sigh    Half-index sigma levels including sigh(1) = sigma(1/2) = 0.0
    !         sigh(plev+1) = sigma(plev+1/2) = 1.0 .  Note that in general
    !         sigh(k) .lt. sig(k)  where sig(k) is the sigma value at the
    !         k_th full-index level.
    ! dsigh   Increment in half-index sigma levels.
    ! kdpmpf  Array of indices of the model full levels which are mapped
    !         into an artificial evenly spaced vertical grid.  Used to aid
    !         in search for vertical position of departure point
    ! kdpmph  Array of indices of the model half levels which are mapped
    !         into an artificial evenly spaced vertical grid.  Used to aid
    !         in search for vertical position of departure point
    ! sigmp   Sigma value at the trajectory midpoint for each gridpoint
    !         in a vertical slice from the global grid.  On entry sigmp is
    !         an initial guess.
    ! sigdp   Sigma value at the trajectory endpoint for each gridpoint
    !         in a vertical slice from the global grid.
    ! kdp     Vertical index for each gridpoint.  This index points into a
    !         vertical slice array whose vertical grid is given by sig.
    !         E.g.,   sig(kdp(i,k)) .le. sigdp(i,k) .lt. sig(kdp(i,k)+1).
    
    !  Local scalars:
    INTEGER :: i, iter, k
    
    !  Local arrays:
    REAL(dp) :: wmp(plon,plev)
    
    !!LK  !  External subroutines
    !!LK  EXTERNAL herzin, kdpfnd, vdplim
    
    
    !  Executable statements
    
    ! Loop over departure point iterates.
    
    DO iter = 1, iterdp
      
      ! Compute midpoint indices in half-index sigma-level arrays (use kdp
      ! as temporary storage).
      
      CALL kdpfnd(plevp1,pmap,sigh,sigmp,kdpmph,kdp)
      
      ! Interpolate sigma dot field to trajectory midpoints using Hermite
      ! cubic interpolant.
      
      CALL herzin(plevp1,1,wb,wst,wsb,sigh,dsigh,sigmp,kdp,wmp)
      
      ! Update estimate of trajectory midpoint.
      
      DO k = 1, plev
        DO i = 1, plon
          sigmp(i,k) = sig(k) - 0.5_dp*dt*wmp(i,k)
        END DO
      END DO

      ! Restrict vertical midpoints to be between the top and bottom half-
      ! index sigma levels.
      
      CALL vdplim(plevp1,sigh,sigmp)
    END DO

    ! Compute trajectory endpoints.
    
    DO k = 1, plev
      DO i = 1, plon
        sigdp(i,k) = sig(k) - dt*wmp(i,k)
      END DO
    END DO

    ! Restrict vertical departure points to be between the top and bottom
    ! full-index sigma levels.
    
    CALL vdplim(plev,sig,sigdp)
    
    ! Vertical indices for trajectory endpoints that point into full-index
    ! sigma level arrays.
    
    CALL kdpfnd(plev,pmap,sig,sigdp,kdpmpf,kdp)
    
  END SUBROUTINE vrtdep
  

#if (! defined __crayx1) || (defined __SX__) || (defined ES)
  SUBROUTINE whenfgt(n,x,incr,target,index,nn)

    !
    ! Authors:
    !
    ! U. Schulzweida, MPI, September 1997
    !
    
    !  Scalar arguments
    REAL(dp) :: target
    INTEGER :: incr, n, nn
    
    !  Array arguments
    REAL(dp) :: x(*)
    INTEGER :: INDEX(*)
    
    !  Local scalars:
    INTEGER :: ival, j
    
    !  Executable statements
    
    ival = 0
    DO j = 1, n, incr
      IF (x(j)>target) THEN
        ival = ival + 1
        INDEX(ival) = j
      END IF
    END DO
    
    nn = ival
    
  END SUBROUTINE whenfgt
#endif

  SUBROUTINE setup_overlap

    ! Set up communication scheme for the slt
    ! Current version does not work if second area of PEs is rotated !!
    ! Current version assumes static overlap in latitudinal direction !!
    ! This routine should be called after suslt in initialise

    ! local arrays
    INTEGER, ALLOCATABLE :: pe_a(:,:)
    INTEGER, ALLOCATABLE :: icountlats(:), ilat(:,:,:), znh(:,:,:)
    INTEGER, ALLOCATABLE :: ilat_own(:,:)
    LOGICAL, ALLOCATABLE :: lpe_tmp(:)

    INTEGER, ALLOCATABLE :: icount_jump(:,:)
    LOGICAL, ALLOCATABLE :: ljump(:)

    ! local scalars
    INTEGER :: row_pe, row_pe_tmp, lower_pe, sndpe_id, jj
    INTEGER :: counter_lats, count2, sndpelats, mypelats, receiver, sender
    INTEGER :: nprocx, nprocy, pid_lt, pid_ln
    INTEGER :: row_end, row_stride, upper_pe, ih_rcv, ih_snd
    INTEGER :: blks_tmp, col_pe, src_ln, src
    INTEGER :: icountlats_own
    LOGICAL :: lfound

    ! Be careful about the handling of pid_lt=1: 
    ! - towards equator: only sends
    ! - Pole extensions of pole PEs filled in extyv/extys and not in overlap. 
    ! - Pole PEs are the only ones which might send latitudes from the southern
    !   overlap (on the southern hemisphere) and from the northern overlap (on
    !   the northern hemisphere).
    !
    ! If owner turns out to be the same as receiver: do not send/receive, just
    ! copy data across equator from northern area to southern area or vice versa 
    ! (in routine overlap).

    nprocy = lc%nproca
    nprocx = lc%nprocb
    pid_lt = lc%set_a
    pid_ln = lc%set_b

    ALLOCATE(pe_snd_eq(nprocy))
    ALLOCATE(nsnd_eq(nprocy))
    ALLOCATE(snd_eq(nxpt+jintmx,2,nprocy))

    ALLOCATE(nrcv_eq(nprocy))
    ALLOCATE(rcv_eq(nxpt+jintmx,2,nprocy))

    ALLOCATE(pe_snd_po(nprocy))
    ALLOCATE(nsnd_po(nprocy))
    ALLOCATE(snd_po(nxpt+jintmx,2,nprocy))

    ALLOCATE(nh(nxpt+jintmx,2,nprocy))

    ALLOCATE(nrcv_po(nprocy))
    ALLOCATE(rcv_po(nxpt+jintmx,2,nprocy))
    ALLOCATE(rcv_po_own(nxpt+jintmx,2))

    ALLOCATE(blks_ea(nprocx))
    ALLOCATE(blks_we(nprocx))

    ! local arrays
    ALLOCATE(pe_a(nxpt+jintmx,nprocy))
    ALLOCATE(icountlats(nprocy))
    ALLOCATE(ilat(nxpt+jintmx,2,nprocy))
    ALLOCATE(ilat_own(nxpt+jintmx,2))
    ALLOCATE(znh(nxpt+jintmx,2,nprocy))
    ALLOCATE(lpe_tmp(nprocy))

    ALLOCATE(icount_jump(2,nprocy))
    ALLOCATE(ljump(2))

    IF (pid_lt == 1) THEN
       mypelats = lc%nglat/2 + nxpt +jintmx
    ELSE
       mypelats = lc%nglat/2
    END IF

    !===========================================================================

    ! Set up scheme for receives of lats that are sent towards the equator,
    ! i.e. receives for southern overlap in southern hemisphere and 
    ! northern overlap in northern hemisphere (except for pid_lt=1).
    ! Determine owner of needed overlap lats and local lat numbers of owner.

    ! determine pe_a globally

    pe_a(:,:) = -999

    DO row_pe = 2,nprocy      ! determine receive info for these PEs
       row_pe_tmp = row_pe
       counter_lats = 0
       DO jj = 1,nxpt+jintmx  ! loop over all lats of an overlap region (for
          ! one hemisphere)

          ! Now loop over all PEs south of current PE (north of current PE for 
          ! northern hemisphere) until the PE is found which owns the current
          ! overlap latitude.
          DO lower_pe = row_pe_tmp - 1, 1, -1
             sndpe_id = lc%mapmesh(1,lower_pe)  ! choose any column, e.g. 1
             !             isndpe_id = indx(sndpe_id,gc)
             IF (lower_pe == 1 ) THEN
                sndpelats = gc(sndpe_id)%nglat/2 + nxpt +jintmx
             ELSE
                sndpelats = gc(sndpe_id)%nglat/2
             END IF
             ! if counter of lats for lower_pe <= total #lats for lower_pe,
             ! i.e. if needed lat is owned by lower_pe, then ...
             IF (counter_lats + 1 <= sndpelats) THEN

                ! row_pe rcvs from lower_pe
                receiver = row_pe
                sender = lower_pe
                pe_a(jj,receiver) = sender
                ! lat was found on this lower_pe; increment counter_lats and
                ! start looking for next lat on this same lower_pe
                counter_lats = counter_lats + 1
                EXIT
             END IF
             ! if lat was not found on this lower_pe, set counter_lats=0 and
             ! loop to next lower PE
             counter_lats = 0
          END DO
          ! for next lat, start searching at last lower_pe (and not at row_pe)
          row_pe_tmp    = lower_pe + 1
       END DO
    END DO

    ! local send values
    pe_snd_eq(:) = -999
    npe_snd_eq = 0

    DO row_pe = 1,nprocy
       DO jj = 1,nxpt+jintmx
          IF (pe_a(jj,row_pe) == pid_lt) THEN  ! If I am a sender...
             ! number of PEs to send to (inludes myself in cases of copying)
             npe_snd_eq = npe_snd_eq + 1
             ! dest id
             pe_snd_eq(npe_snd_eq) = lc%mapmesh(pid_ln,row_pe)
             EXIT
          END IF
       END DO
    END DO

    icountlats(:) = 0
    ilat(:,:,:) = -999
    nsnd_eq(:) = 0
    snd_eq(:,:,:) = -999

    DO row_pe =1,nprocy
       DO jj=1,nxpt+jintmx
          IF (pe_a(jj,row_pe) == pid_lt) THEN
             ! icountlats counts the lats row_pe needs from me (per hemisph.)
             icountlats(row_pe) = icountlats(row_pe) + 1
             ! local lats (S)
             ilat(icountlats(row_pe),1,row_pe) = jstart + plato2 - &
                  &                                                icountlats(row_pe)
             ! local lats (N)        
             ilat(icountlats(row_pe),2,row_pe) = jstart + &
                  &                                                icountlats(row_pe) - 1 
          END IF
       END DO
    END DO

    nsnd_eq(:npe_snd_eq) = PACK(icountlats, icountlats /= 0)

    count2=0
    DO row_pe=1,nprocy
       IF (ilat(1,1,row_pe) /= -999) THEN
          count2 = count2 + 1
          snd_eq(:,:,count2) = ilat(:,:,row_pe)
       ENDIF
    ENDDO

    ! local receive values
    lpe_tmp(:) = .FALSE.
    npe_rcv_eq = 0

    row_pe = pid_lt
    DO jj = 1,nxpt+jintmx
       ! If I am a real receiver ...
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) /= pid_lt)) THEN
          ! number of PEs to receive from (does NOT include myself in cases
          ! of copying; these cases are handled in the send section only)
          lpe_tmp(pe_a(jj,row_pe)) = .TRUE.
       END IF
    END DO
    npe_rcv_eq = COUNT(lpe_tmp)

    icountlats(:) = 0
    ilat(:,:,:) = -999
    nrcv_eq(:) = 0
    rcv_eq(:,:,:) = -999

    row_pe = pid_lt
    DO jj = 1,nxpt+jintmx
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) /= pid_lt)) THEN
          ! icountlats counts the lats I need from pe_a(jj,row_pe) (per hemi.)
          sender = pe_a(jj,row_pe)
          icountlats(sender) = icountlats(sender) + 1
          ! local lats (S)
          ilat(icountlats(sender),1,sender) = nxpt + jintmx - jj + 1
          ! local lats (N)
          ilat(icountlats(sender),2,sender) = jstart + plato2 + jj - 1
       END IF
    END DO

    nrcv_eq(:) = icountlats(:)
    rcv_eq(:,:,:) = ilat(:,:,:)

    !============================================================================

    ! Set up scheme for receives of lats that are sent towards the poles and
    ! possibly across the equator,
    ! i.e. receives for northern overlap in southern hemisphere and 
    ! southern overlap in northern hemisphere.
    ! Determine owner of needed overlap lats and local lat numbers of owner.

    ! determine pe_a globally
    pe_a(:,:) = -999

    DO row_pe = 1,nprocy      ! determine receive info for these PEs
       row_pe_tmp = row_pe
       counter_lats = 0
       row_end = nprocy
       row_stride = 1
       upper_pe = row_pe_tmp + 1
       lfound = .FALSE.
       DO jj = 1,nxpt+jintmx  ! loop over all lats of an overlap region (for
          ! one hemisphere)

          ! Now loop over all PEs north of current PE (south of current PE for 
          ! northern hemisphere) until the PE is found which owns the current
          ! overlap latitude.
10        CONTINUE
          DO upper_pe = row_pe_tmp + 1, row_end, row_stride
             sndpe_id = lc%mapmesh(1,upper_pe)  ! choose any column, e.g. 1
             !             isndpe_id = indx(sndpe_id,gc)
             IF (upper_pe == 1 ) THEN
                sndpelats = gc(sndpe_id)%nglat/2 + nxpt +jintmx
             ELSE
                sndpelats = gc(sndpe_id)%nglat/2
             END IF
             ! if counter of lats for upper_pe <= total #lats for upper_pe,
             ! i.e. if needed lat is owned by upper_pe, then ...
             IF (counter_lats + 1 <= sndpelats) THEN

                ! row_pe rcvs from upper_pe
                receiver = row_pe
                sender = upper_pe
                pe_a(jj,receiver) = sender
                ! lat was found on this upper_pe; increment counter_lats and
                ! start looking for next lat on this same upper_pe
                counter_lats = counter_lats + 1
                lfound = .TRUE.
                EXIT
             END IF
             ! if lat was not found on this upper_pe, set counter_lats=0 and
             ! loop to next upper PE
             counter_lats = 0
          END DO
          ! for next lat, start searching at last upper_pe
          IF (lfound) THEN
             row_pe_tmp    = upper_pe - 1; lfound = .FALSE.; CYCLE
          END IF
          ! if not found on this hemisphere, jump to other hemisphere
          ! and repeat search for current jj
          IF (.NOT. lfound) THEN 
             row_stride = -1 
             row_end = 1
             row_pe_tmp = upper_pe - 2
             lfound = .FALSE.
             GO TO 10
          END IF
       END DO
    END DO

    ! local send values
    pe_snd_po(:) = -999
    npe_snd_po = 0

    DO row_pe = 1,nprocy
       DO jj = 1,nxpt+jintmx
          IF (pe_a(jj,row_pe) == pid_lt) THEN   ! If I am a sender...
             ! number of PEs to send to (inludes myself in cases of copying)
             npe_snd_po = npe_snd_po + 1
             ! dest index
             pe_snd_po(npe_snd_po) = lc%mapmesh(pid_ln,row_pe)
             EXIT
          END IF
       END DO
    END DO

    icountlats(:) = 0
    ilat(:,:,:) = -999
    nsnd_po(:) = 0
    snd_po(:,:,:) = -999

    nh(:,:,:) = 0
    znh(:,:,:) = 0

    ljump(:) =.FALSE.
    icount_jump(:,:) = 0

    DO row_pe =1,nprocy
       DO jj=1,nxpt+jintmx
          IF (pe_a(jj,row_pe) == pid_lt) THEN
             ! icountlats counts the lats row_pe needs from me (per hemisph.)
             icountlats(row_pe) = icountlats(row_pe) + 1
             !================================================
             ! my local lats for northern overlap in southern hem. of receiver
             ih_rcv = 1 ! sending to this hemsiphere of receiver
             IF (pid_lt > row_pe) THEN   ! if I am north of row_pe on southern 
                ! hem., I might have to send from 
                ! both my hemispheres
                ih_snd = 1 ! my hemisphere to send from
                IF (icountlats(row_pe) > mypelats) THEN
                   ih_snd = 2 ! jump to my other hem.
                   ljump(ih_rcv) = .TRUE.
                END IF
             ELSE  ! if I am south of row_pe on southern hem. (or same), I only
                ! have to send from my northern hemisphere
                ih_snd = 2
             END IF
             IF (.NOT.(ljump(ih_rcv))) THEN
                ilat(icountlats(row_pe),ih_rcv,row_pe) = jstart + &
                     &                                                       icountlats(row_pe) - 1 
             ELSE
                icount_jump(ih_rcv,row_pe) = icount_jump(ih_rcv,row_pe) + 1
                ilat(icountlats(row_pe),ih_rcv,row_pe) = jstart + &
                     &                                                icount_jump(ih_rcv,row_pe) - 1
             END IF
             znh(icountlats(row_pe),ih_rcv,row_pe) = ih_snd
             !==================================================
             ! my local lats for southern overlap in northern hem. of receiver
             ih_rcv = 2
             IF (pid_lt > row_pe) THEN ! if I am south of row_pe on northern..
                ih_snd = 2
                IF (icountlats(row_pe) > mypelats) THEN
                   ih_snd = 1 ! jump to my other hem.
                   ljump(ih_rcv) = .TRUE.
                END IF
             ELSE ! if I am north of row_pe on northern hem. (or same), I only
                ! have to send from my southern hem.
                ih_snd = 1
             END IF
             IF (.NOT.(ljump(ih_rcv))) THEN             
                ilat(icountlats(row_pe),ih_rcv,row_pe) = jstart + plato2 - &
                     &                                                        icountlats(row_pe)
             ELSE
                icount_jump(ih_rcv,row_pe) = icount_jump(ih_rcv,row_pe) + 1
                ilat(icountlats(row_pe),ih_rcv,row_pe) = jstart + plato2 - &
                     &                                                   icount_jump(ih_rcv,row_pe)
             END IF
             znh(icountlats(row_pe),ih_rcv,row_pe) = ih_snd
             !==================================================
          END IF
       END DO
    END DO

    nsnd_po(:npe_snd_po) = PACK(icountlats, icountlats /= 0)

    count2=0
    DO row_pe=1,nprocy
       IF (ilat(1,1,row_pe) /= -999) THEN
          count2=count2+1
          snd_po(:,:,count2) = ilat(:,:,row_pe)
          nh(:,:,count2)     = znh(:,:,row_pe)
       ENDIF
    ENDDO

    ! local receive values
    lpe_tmp(:) = .FALSE.
    npe_rcv_po = 0

    row_pe = pid_lt
    DO jj = 1,nxpt+jintmx
       ! If I am a receiver ...
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) /= pid_lt)) THEN
          ! number of PEs to receive from (does NOT include myself in cases
          ! of copying; these cases are handled in the send section only)
          lpe_tmp(pe_a(jj,row_pe)) = .TRUE.
       END IF
    END DO
    npe_rcv_po = COUNT(lpe_tmp)

    icountlats(:) = 0
    ilat(:,:,:) = -999
    nrcv_po(:) = 0
    rcv_po(:,:,:) = -999

    icountlats_own = 0
    ilat_own(:,:) = -999
    nrcv_po_own = 0
    rcv_po_own(:,:) = -999

    row_pe = pid_lt
    DO jj = 1,nxpt+jintmx
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) /= pid_lt)) THEN
          ! icountlats counts the lats I need from pe_a(jj,row_pe) (per hemi.)
          sender = pe_a(jj,row_pe)
          icountlats(sender) = icountlats(sender) + 1
          ! local lats (S)
          ilat(icountlats(sender),1,sender) = jstart + plato2 + jj - 1
          ! local lats (N)
          ilat(icountlats(sender),2,sender) = nxpt + jintmx - jj + 1
       END IF
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) == pid_lt)) THEN
          ! icountlats_own counts the lats I need from myself (per hemi.)
          icountlats_own = icountlats_own + 1
          ! local lats (S)
          ilat_own(icountlats_own,1) = jstart + plato2 + jj - 1
          ! local lats (N)
          ilat_own(icountlats_own,2) = nxpt + jintmx - jj + 1
       END IF
    END DO

    nrcv_po(:) = icountlats(:)
    rcv_po(:,:,:) = ilat(:,:,:)

    nrcv_po_own = icountlats_own
    rcv_po_own(:,:) = ilat_own(:,:)

    !============================================================================

    ! the major part of the setup for the longitudinal exchange must be 
    ! carried out dynamically, since the longitudinal overlap is dynamic 
    ! (the latitudinal overlap currently is static)

    ! blks_we(row_pe)/blks_ea(row_pe): 
    ! this many lons available from other PEs closer (than row_pe) to myself 
    ! in the west/east
    IF (nprocx > 1) THEN

       blks_we(:) = 0
       blks_tmp = 0
       col_pe = pid_ln
       DO src_ln = 1,nprocx - 1 ! this many PEs to the west of myself
          col_pe = col_pe - 1
          IF (col_pe < 1) col_pe = nprocx
          src = lc%mapmesh(col_pe,pid_lt)
          !          isrc = indx(src,gc)
          blks_we(col_pe) = blks_tmp
          blks_tmp = blks_tmp + gc(src)%nglon
       END DO

       blks_ea(:) = 0
       blks_tmp = 0
       col_pe = pid_ln
       DO src_ln = 1,nprocx - 1 ! this many PEs to the east of myself
          col_pe = col_pe + 1
          IF (col_pe > nprocx) col_pe = 1
          src = lc%mapmesh(col_pe,pid_lt)
          !          isrc = indx(src,gc)
          blks_ea(col_pe) = blks_tmp
          blks_tmp = blks_tmp + gc(src)%nglon
       END DO

    END IF

    !============================================================================


    ! deallocate local arrays
    DEALLOCATE(pe_a)
    DEALLOCATE(icountlats)
    DEALLOCATE(ilat)
    DEALLOCATE(ilat_own)
    DEALLOCATE(znh)
    DEALLOCATE(lpe_tmp)

    DEALLOCATE(icount_jump)
    DEALLOCATE(ljump)

  END SUBROUTINE setup_overlap

  SUBROUTINE overlap(ub,vb,fb)

    REAL(dp) :: ub(plond,plev,platd,2)
    REAL(dp) :: vb(plond,plev,platd,2)
    REAL(dp) :: fb(plond,plev,pcnst,platd,2)

    ! local arrays
    REAL(dp) :: buf_snd_eq(plond*plev*(nxpt+jintmx)*(2+pcnst)*2, npe_snd_eq)
    REAL(dp) :: buf_rcv_eq(plond*plev*(nxpt+jintmx)*(2+pcnst)*2, npe_rcv_eq)
    REAL(dp) :: buf_snd_po(plond*plev*(nxpt+jintmx)*(2+pcnst)*2, npe_snd_po)
    REAL(dp) :: buf_rcv_po(plond*plev*(nxpt+jintmx)*(2+pcnst)*2, npe_rcv_po)

    REAL(dp) :: buf_snd_ea(nxpt*plev*platd*(2+pcnst)*2,     lc%nprocb - 1)
    REAL(dp) :: buf_rcv_ea(nxpt*plev*platd*(2+pcnst)*2,     lc%nprocb - 1)
    REAL(dp) :: buf_snd_we((nxpt+1)*plev*platd*(2+pcnst)*2, lc%nprocb - 1)
    REAL(dp) :: buf_rcv_we((nxpt+1)*plev*platd*(2+pcnst)*2, lc%nprocb - 1)

    INTEGER :: tagtable_eq(lc%nproca), tagtable_po(lc%nproca)
    INTEGER :: tagtable_we(lc%nprocb), tagtable_ea(lc%nprocb)
    INTEGER :: nxpt_tmp(plev,platd,2)

    ! local scalars
    INTEGER :: nprocy, nprocx, jproc, imeslen, ih, ihrcv, n, lat, k, lon, m, i
    INTEGER :: tag_base_eq, tag_base_po, tagcount_eq, tagcount_po
    INTEGER :: tag_base_we, tag_base_ea, tagcount_we, tagcount_ea
    INTEGER :: tag, count, dest, idest, source, isource, src_lt, src_ln
    INTEGER :: pid_ln, pid_lt, col_pe, diff, iblkbeg, iblkend, total, elements

    !===========================================================================
    !  Update the overlap between processors.
    !  The arrays have been extended as if poles were present.
    !
    !  On entry:
    !    ub, vb and fb are extended arrays with pole extensions
    !    filled for processor zero and interior values
    !    set for all processors.
    !
    !  On exit:
    !    ub, vb and fb have been filled with the extended values
    !    from neighboring processors. The interior
    !    values are unchanged.
    !
    !=========================================================================

    nprocy = lc%nproca
    nprocx = lc%nprocb
    pid_ln = lc%set_b
    pid_lt = lc%set_a

    ! send towards equator

    tag_base_eq = 400+10+nprocx+2+nprocx   ! last change occurred in sltini...

    DO jproc = 1,npe_snd_eq
       imeslen = 0
       DO ih = 1,2
          DO n = 1,nsnd_eq(jproc) ! this many lats per my current hem.
             lat = snd_eq(n,ih,jproc)
             DO k = 1,plev
                DO lon = 1,plond
                   imeslen = imeslen + 1
                   buf_snd_eq(imeslen,jproc) = ub(lon,k,lat,ih)
                   imeslen = imeslen +1 
                   buf_snd_eq(imeslen,jproc) = vb(lon,k,lat,ih)
                END DO
             END DO
          END DO
       END DO

       DO ih = 1,2
          DO n = 1,nsnd_eq(jproc) ! this many lats per my current hem.
             lat = snd_eq(n,ih,jproc)
             DO m =1,pcnst
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      buf_snd_eq(imeslen,jproc) = fb(lon,k,m,lat,ih)
                   END DO
                END DO
             END DO
          END DO
       END DO

       tag = tag_base_eq + jproc
       count = imeslen
       dest = gc(pe_snd_eq(jproc))%pe
       CALL p_isend(buf_snd_eq(1,jproc),dest,tag,p_count=count)

    END DO

    ! send towards poles

    tag_base_po = tag_base_eq + nprocy

    DO jproc = 1,npe_snd_po
       imeslen = 0
       DO ihrcv = 1,2
          DO n = 1,nsnd_po(jproc) ! this many lats per hem. of receiver
             ! (note the difference to the send towards equator) 
             lat = snd_po(n,ihrcv,jproc)
             ih = nh(n,ihrcv,jproc)
             DO k = 1,plev
                DO lon = 1,plond
                   imeslen = imeslen + 1
                   buf_snd_po(imeslen,jproc) = ub(lon,k,lat,ih)
                   imeslen = imeslen +1 
                   buf_snd_po(imeslen,jproc) = vb(lon,k,lat,ih)
                END DO
             END DO
          END DO
       END DO

       DO ihrcv = 1,2
          DO n = 1,nsnd_po(jproc) ! this many lats per hem. of receiver
             ! (note the difference to the send towards equator) 
             lat = snd_po(n,ihrcv,jproc)
             ih = nh(n,ihrcv,jproc)
             DO m = 1,pcnst
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      buf_snd_po(imeslen,jproc) = fb(lon,k,m,lat,ih)
                   END DO
                END DO
             END DO
          END DO
       END DO

       tag   = tag_base_po + jproc 
       count = imeslen
       dest  = gc(pe_snd_po(jproc))%pe

       IF (dest /= lc%pe) THEN
          CALL p_isend(buf_snd_po(1,jproc),dest,tag,p_count=count)
       ELSE ! copy across equator
          imeslen = 0
          DO ih =1,2
             DO n = 1,nrcv_po_own
                lat = rcv_po_own(n,ih)
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      ub(lon,k,lat,ih) = buf_snd_po(imeslen,jproc)
                      imeslen = imeslen + 1
                      vb(lon,k,lat,ih) = buf_snd_po(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO

          DO ih = 1,2
             DO n = 1,nrcv_po_own
                lat = rcv_po_own(n,ih)
                DO m = 1,pcnst
                   DO k = 1,plev
                      DO lon = 1,plond
                         imeslen = imeslen + 1
                         fb(lon,k,m,lat,ih) = buf_snd_po(imeslen,jproc)
                      END DO
                   END DO
                END DO
             END DO
          END DO

       END IF
    END DO

    ! receive latitudes sent towards equator

    ! nprocy could be changed to npe_rcv_eq
    ! this would imply changes at other places...
    tagcount_eq = nprocy
    DO jproc = 1,nprocy
       tagtable_eq(jproc) = tag_base_eq + jproc
    END DO

    DO jproc = 1,npe_rcv_eq

       ! look for tags, not for sources; sources are not unique, could be 
       ! from pole sends as well.
       CALL p_probe(tagcount_eq,tagtable_eq,source,tag,count)
       CALL p_recv(buf_rcv_eq(1,jproc),source,tag,p_count=count)

       isource = indx(source,gc)
       src_lt = gc(isource)%set_a

       imeslen = 0
       DO ih = 1,2
          DO n = 1,nrcv_eq(src_lt)
             lat = rcv_eq(n,ih,src_lt)
             DO k = 1,plev
                DO lon = 1,plond
                   imeslen = imeslen + 1
                   ub(lon,k,lat,ih) = buf_rcv_eq(imeslen,jproc)
                   imeslen = imeslen +1 
                   vb(lon,k,lat,ih) = buf_rcv_eq(imeslen,jproc)
                END DO
             END DO
          END DO
       END DO

       DO ih = 1,2
          DO n = 1,nrcv_eq(src_lt)
             lat = rcv_eq(n,ih,src_lt)
             DO m = 1,pcnst
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      fb(lon,k,m,lat,ih) = buf_rcv_eq(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO
       END DO

    END DO

    ! receive latitudes sent towards poles 

    ! nprocy could be changed to npe_rcv_po
    ! this would imply changes at other places...
    tagcount_po = nprocy
    DO jproc = 1,nprocy
       tagtable_po(jproc) = tag_base_po + jproc
    END DO

    DO jproc = 1,npe_rcv_po

       ! look for tags, not for sources; sources are not unique, could be 
       ! from pole sends as well.
       CALL p_probe(tagcount_po,tagtable_po,source,tag,count)
       CALL p_recv(buf_rcv_po(1,jproc),source,tag,p_count=count)

       isource = indx(source,gc)
       src_lt = gc(isource)%set_a

       imeslen = 0
       DO ih = 1,2
          DO n = 1,nrcv_po(src_lt)
             lat = rcv_po(n,ih,src_lt)
             DO k = 1,plev
                DO lon = 1,plond
                   imeslen = imeslen + 1
                   ub(lon,k,lat,ih) = buf_rcv_po(imeslen,jproc)
                   imeslen = imeslen +1 
                   vb(lon,k,lat,ih) = buf_rcv_po(imeslen,jproc)
                END DO
             END DO
          END DO
       END DO

       DO ih = 1,2
          DO n = 1,nrcv_po(src_lt)
             lat = rcv_po(n,ih,src_lt)
             DO m = 1,pcnst
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      fb(lon,k,m,lat,ih) = buf_rcv_po(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO
       END DO

    END DO

    !===================================================

    IF (nprocx > 1) THEN

       ! send towards east

       tag_base_ea = tag_base_po + nprocy

       col_pe = pid_ln
       diff = 0

       col:   DO jproc = 1,nprocx-1  ! max # of PEs to send to

          col_pe = col_pe + 1
          IF (col_pe > nprocx) col_pe = 1
          dest = gc(lc%mapmesh(col_pe,pid_lt))%pe
          tag = tag_base_ea + jproc

          imeslen = 0
          DO ih = 1,2
             DO lat = 1,platd
                DO k = 1,plev
                   iblkbeg = 1
                   iblkend = MIN(plon,nxpt_a(k,lat,ih) - diff)
                   DO i = iblkbeg, iblkend
                      lon = istart + plon - i
                      imeslen = imeslen + 1
                      buf_snd_ea(imeslen,jproc) = ub(lon,k,lat,ih)
                      imeslen = imeslen +1 
                      buf_snd_ea(imeslen,jproc) = vb(lon,k,lat,ih)
                   END DO
                END DO
             END DO
          END DO

          IF (imeslen == 0) EXIT col

          DO ih = 1,2
             DO lat = 1,platd
                DO m = 1,pcnst
                   DO k = 1,plev
                      iblkbeg = 1
                      iblkend = MIN(plon,nxpt_a(k,lat,ih) - diff)
                      DO i = iblkbeg, iblkend
                         lon = istart + plon - i
                         imeslen = imeslen + 1
                         buf_snd_ea(imeslen,jproc) = fb(lon,k,m,lat,ih)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          count = imeslen
          CALL p_isend(buf_snd_ea(1,jproc),dest,tag,p_count=count)

          idest = indx(dest,gc)
          diff = diff + gc(idest)%nglon

       END DO col

       ! send towards west
       ! note the extra line sent towards the west (=to the eastern overlaps)

       tag_base_we = tag_base_ea + nprocx

       col_pe = pid_ln
       diff = 0

       col2:   DO jproc = 1,nprocx-1  ! max # of PEs to send to

          col_pe = col_pe - 1
          IF (col_pe < 1) col_pe = nprocx
          dest = gc(lc%mapmesh(col_pe,pid_lt))%pe
          tag = tag_base_we + jproc

          imeslen = 0
          DO ih = 1,2
             DO lat = 1,platd
                DO k = 1,plev
                   iblkbeg = 1
                   iblkend = MIN(plon,nxpt_a(k,lat,ih) + 1 - diff)
                   DO i = iblkbeg, iblkend
                      lon = istart - 1 + i
                      imeslen = imeslen + 1
                      buf_snd_we(imeslen,jproc) = ub(lon,k,lat,ih)
                      imeslen = imeslen +1 
                      buf_snd_we(imeslen,jproc) = vb(lon,k,lat,ih)
                   END DO
                END DO
             END DO
          END DO

          IF (imeslen == 0) EXIT col2

          DO ih = 1,2
             DO lat = 1,platd
                DO m = 1,pcnst
                   DO k = 1,plev
                      iblkbeg = 1
                      iblkend = MIN(plon,nxpt_a(k,lat,ih) + 1 - diff)
                      DO i = iblkbeg, iblkend
                         lon = istart - 1 + i
                         imeslen = imeslen + 1
                         buf_snd_we(imeslen,jproc) = fb(lon,k,m,lat,ih)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          count = imeslen
          CALL p_isend(buf_snd_we(1,jproc),dest,tag,p_count=count)

          idest = indx(dest,gc)
          diff = diff + gc(idest)%nglon

       END DO col2

       ! receive longitudes sent towards east (i.e. receive my western overlap)

       ! nprocx could be changed to actual number of PEs to receive from
       ! this would imply changes at other places...
       tagcount_ea = nprocx - 1
       DO jproc = 1,nprocx - 1
          tagtable_ea(jproc) = tag_base_ea + jproc
       END DO

       total = (2+pcnst)*SUM(nxpt_a)
       elements = 0

       DO jproc = 1,nprocx - 1 ! max. # of PEs to receive from

          CALL p_probe(tagcount_ea,tagtable_ea,source,tag,count)
          CALL p_recv(buf_rcv_ea(1,jproc),source,tag,p_count=count)

          isource = indx(source,gc)
          src_ln = gc(isource)%set_b

          imeslen = 0
          DO ih = 1,2
             DO lat = 1,platd
                DO k = 1,plev
                   iblkbeg = blks_we(src_ln) + 1
                   iblkend = iblkbeg + MIN(gc(isource)%nglon, nxpt_a(k,lat,ih) - &
                        &                                                          blks_we(src_ln)) - 1
                   DO i = iblkbeg, iblkend
                      lon = istart - i                
                      imeslen = imeslen + 1
                      ub(lon,k,lat,ih) = buf_rcv_ea(imeslen,jproc)
                      imeslen = imeslen +1 
                      vb(lon,k,lat,ih) = buf_rcv_ea(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO

          DO ih = 1,2
             DO lat = 1,platd
                DO m = 1,pcnst
                   DO k = 1,plev
                      iblkbeg = blks_we(src_ln) + 1
                      iblkend = iblkbeg + MIN(gc(isource)%nglon, nxpt_a(k,lat,ih) - &
                           &                                                          blks_we(src_ln)) - 1
                      DO i = iblkbeg, iblkend
                         lon = istart - i                                   
                         imeslen = imeslen + 1
                         fb(lon,k,m,lat,ih) = buf_rcv_ea(imeslen,jproc)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          elements = elements + imeslen
          IF (elements == total) EXIT 

       END DO

       ! receive longitudes sent towards west (i.e. receive my eastern overlap)

       ! nprocx could be changed to actual number of PEs to receive from
       ! this would imply changes at other places...
       tagcount_we = nprocx - 1
       DO jproc = 1,nprocx - 1
          tagtable_we(jproc) = tag_base_we + jproc
       END DO

       nxpt_tmp(:,:,:) = nxpt_a(:,:,:) + 1
       total = (2+pcnst)*SUM(nxpt_tmp)
       elements = 0

       DO jproc = 1,nprocx - 1 ! max. # of PEs to receive from

          CALL p_probe(tagcount_we,tagtable_we,source,tag,count)
          CALL p_recv(buf_rcv_we(1,jproc),source,tag,p_count=count)

          isource = indx(source,gc)
          src_ln = gc(isource)%set_b

          imeslen = 0
          DO ih = 1,2
             DO lat = 1,platd
                DO k = 1,plev
                   iblkbeg = blks_ea(src_ln) + 1
                   iblkend =iblkbeg + MIN(gc(isource)%nglon, nxpt_a(k,lat,ih) + 1&
                        &                                                        - blks_ea(src_ln)) - 1
                   DO i = iblkbeg, iblkend
                      lon = istart + plon - 1 + i            
                      imeslen = imeslen + 1
                      ub(lon,k,lat,ih) = buf_rcv_we(imeslen,jproc)
                      imeslen = imeslen +1 
                      vb(lon,k,lat,ih) = buf_rcv_we(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO

          DO ih = 1,2
             DO lat = 1,platd
                DO m = 1,pcnst
                   DO k = 1,plev
                      iblkbeg = blks_ea(src_ln) + 1
                      iblkend =iblkbeg + MIN(gc(isource)%nglon, nxpt_a(k,lat,ih) + 1&
                           &                                                        - blks_ea(src_ln)) - 1
                      DO i = iblkbeg, iblkend
                         lon = istart + plon - 1 + i                
                         imeslen = imeslen + 1
                         fb(lon,k,m,lat,ih) = buf_rcv_we(imeslen,jproc)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          elements = elements + imeslen
          IF (elements == total) EXIT 

       END DO

       !========================================================================

       ! end of the multiple longitudinal processors if
    END IF

    CALL p_wait

    !    CALL p_barrier

  END SUBROUTINE overlap

!!============================================================================

  SUBROUTINE init_mass_fixer

    ! Description:
    !
    ! Set up mass fixer
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, May 1995, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_advection, ONLY: jps
    USE mo_tracdef,   ONLY: trlist

    IMPLICIT NONE

    !  Local scalars: 

    INTEGER :: jt, ita

    !  Executable statements 

    ! Determine fixer type for transported constituents

    kftype(1:jps) = 1    ! fixer type for model base variables (like q, ql, qi)

    ita = 0
    DO jt = 1, trlist% ntrac
      if (trlist% ti(jt)% ntran == 0) cycle
      ita = ita + 1
      kftype(jps+ita) = trlist% ti(jt)% nfixtyp

      IF (kftype(jps+ita) < 0 .OR. kftype(jps+ita) > 2) THEN
        WRITE (nout, *) ' Wrong fixer type for variable no. ', jt
        WRITE (nout, *) '   type = ', kftype(jps+ita)
        CALL finish('init_mass_fixer','Run terminated.')
      END IF
    END DO

  END SUBROUTINE init_mass_fixer

#ifndef SLDIAG
  SUBROUTINE mass_fixer(qfcst, alpha, psm1cor, pscor, ztodt)
#else
  SUBROUTINE mass_fixer(qfcst, alpha, psm1cor, pscor, ztodt, hqfcst)
#endif
    
    ! Description:
    !
    ! Mass fixing of constituents and air mass, tendencies
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, March 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! T. Diehl, DKRZ, July 1999, parallel version
    ! 
    ! for more details see file AUTHORS
    !
    
    USE mo_advection,   ONLY: jps
    USE mo_tracdef,     ONLY: trlist
    USE mo_scan_buffer, ONLY: alps, qte, xlte, xite, xtte
    USE mo_memory_g1a,  ONLY: alpsm1
    
    !  Scalar arguments 
    REAL(dp), INTENT(in) :: pscor, psm1cor, ztodt
    
    !  Array arguments 
    REAL(dp), INTENT(in)    :: alpha(pcnst)
    REAL(dp), INTENT(inout) :: qfcst(plond,plev,pcnst,platd,2)
    
    !  Local scalars: 
    REAL(dp) :: zlpsc, zlpsm1c, ztmst
    INTEGER :: ilat, jlat, ita, jcen, jk, jl, jt, ih
    
    !  Local arrays: 
#ifdef SLDIAG
    REAL(dp) :: qf3m(plon,plev,pcnst), vqfm(plon,plev,pcnst)
    REAL(dp) :: hqfcst(plond,plev,pcnst,plat)
#endif
    
    !  Executable statements 
    
!    IF (lresume) THEN
!      ztmst = 0.5_dp*ztodt
!    ELSE
!      ztmst = ztodt
!    END IF
    ztmst = ztodt
   
    DO jlat  = 1, lc% nglat  ! local  index north -> south
      ilat = plat + 1 - jlat  ! continuous local index for hqfcst
      IF (ilat > plato2) THEN
        ih = 2
        jcen = js - 1 + ilat - plato2
      ELSE
        ih = 1
        jcen = js - 1 + ilat
      END IF
      
      !-- 1. Fix mass of semi lagrangian transported scalars
      
#ifdef SLDIAG
      CALL fixer(ztmst, alpha, fb(1,1,1,jcen,ih), qfcst(1,1,1,jcen,ih), &
                 hqfcst(1,1,1,ilat), vqfm, qf3m)
#else
      CALL fixer(alpha, fb(1,1,1,jcen,ih), qfcst(1,1,1,jcen,ih))
#endif
      
      !-- 2. Compute tendencies of semi lagrangian transport
        
      DO jk = 1, plev
!DIR$ IVDEP
!OCL NOVREC
        DO jl = 1, plon
          qte(jl,jk,jlat)                                            &
               = qte(jl,jk,jlat) + (qfcst(jl+nxpt,jk,1,jcen,ih)-     &
               fb(jl+nxpt,jk,1,jcen,ih))/ztmst
          xlte(jl,jk,jlat)                                            &
               = xlte(jl,jk,jlat) + (qfcst(jl+nxpt,jk,2,jcen,ih)-     &
               fb(jl+nxpt,jk,2,jcen,ih))/ztmst
          xite(jl,jk,jlat)                                            &
               = xite(jl,jk,jlat) + (qfcst(jl+nxpt,jk,2,jcen,ih)-     &
               fb(jl+nxpt,jk,2,jcen,ih))/ztmst
        END DO
      END DO
      
      ita = 0
      DO jt = 1, trlist% ntrac
        IF (trlist% ti(jt)% ntran /= 0) THEN
          ita = ita + 1
          DO jk = 1, plev
!DIR$ IVDEP
!OCL NOVREC
            DO jl = 1, plon
              xtte(jl,jk,jt,jlat) =                         &
                   xtte(jl,jk,jt,jlat) +                    &
                   (qfcst(jl+nxpt,jk,ita+jps,jcen,ih) -        &
                   fb(jl+nxpt,jk,ita+jps,jcen,ih))/ztmst
            END DO
          END DO
        END IF
      END DO
      
      !-- 3. Fix mass of air
      
      zlpsm1c = LOG(psm1cor)
      zlpsc = LOG(pscor)
      
      DO jl = 1, plon
        alpsm1(jl,jlat) = alpsm1(jl,jlat) + zlpsm1c
        alps  (jl,jlat) = alps  (jl,jlat) + zlpsc
      END DO
    END DO
    
  END SUBROUTINE mass_fixer

#ifndef SLDIAG
  SUBROUTINE fixer(alpha, qin, qout)
#else
  SUBROUTINE fixer(ztodt, alpha, qin, qout, hqfm, vqfm, qf3m)
#endif

    ! Description:
    !
    ! Computes  1) time tendency due to vertical advection.
    !           2) time tendency due to mass adjustment of the constituents.
    !
    ! Method:
    !
    ! Also, the constituent fields are updated based upon the fixer (mass
    ! adjustment) tendency.
    !
    ! qfix=alpha*F*q2*|q2 - q1|**Beta.
    ! Two options are available:
    !  1. kftype=1 : F=1.  and Beta=1.5
    !  2. kftype=2 : F=eta and Beta=1.
    !
    ! .On input
    ! lat     latitude index (values run from 1 to nlat; from southern-most
    !         to northern-most latitude).
    ! ztodt   length of time step (in seconds)
    ! alpha   array of fixer coefficients (one for each constituent) used
    !         in computing the mass adjustment tendency
    ! qin     constituent fields from the previous time step
    ! hqfm    horizontal time tendency  (optional)
    ! qout    set of FORECASTED constituent fields
    ! etamid  vert. coord. at full levels
    ! kftype  type of mass fixer
    !
    ! .On return.
    ! vqfm    vertical time tendency  (optional)
    ! qfm3    fixer time tendency     (optional)
    ! qout    set of FORECASTED AND FIXED constituent fields
    !
    ! .Local parameters: required.
    ! plev    Number of levels in global grid.
    ! plon    Number of longitudes in global grid.
    ! plond   number of longitudes in extended grid.
    ! pcnst   Number of constituents
        
    !  Scalar arguments 
#ifdef SLDIAG
    REAL(dp) :: ztodt
#endif
    
    !  Array arguments 
    REAL(dp) :: alpha(pcnst), &
                qin(plond,plev,pcnst), qout(plond,plev,pcnst) 
#ifdef SLDIAG
    REAL(dp) :: qf3m(plon,plev,pcnst), hqfm(plond,plev,pcnst), &
                vqfm(plon,plev,pcnst)
#endif

    !  Local scalars: 
    INTEGER :: i, jc, k
    
    !  Intrinsic functions 
    INTRINSIC ABS, MAX, SQRT
    
    !  Executable statements 
#ifdef SLDIAG
    DO jc = 1, pcnst
      DO k = 1, plev
        DO i = istart, istop
          vqfm(i-nxpt,k,jc) = (qout(i,k,jc)-qin(i,k,jc))/ztodt - &
               hqfm(i-nxpt,k,jc)
        END DO
      END DO
    END DO
#endif
    DO jc = 1, pcnst
          
      IF (kftype(jc) == 1) THEN
        
        DO k = 1, plev
          DO i = istart, istop
#ifdef SLDIAG
            qf3m(i-nxpt,k,jc) = alpha(jc)*qout(i,k,jc) &
                 *(SQRT(ABS(qout(i,k,jc)-qin(i,k,jc))))**3/ztodt
            qout(i,k,jc) = qout(i,k,jc) + ztodt*qf3m(i-nxpt,k,jc)
#else
            qout(i,k,jc) = qout(i,k,jc) + alpha(jc)*qout(i,k,jc) &
                 *(SQRT(ABS(qout(i,k,jc)-qin(i,k,jc))))**3
#endif
          END DO
        END DO
        
      ELSE IF (kftype(jc)==2) THEN
        
        DO k = 1, plev
          DO i = istart, istop
#ifdef SLDIAG
            qf3m(i-nxpt,k,jc) = alpha(jc)*qout(i,k,jc)*etamid(k)* &
                 ABS(qout(i,k,jc)-qin(i,k,jc))/ztodt
            qout(i,k,jc) = qout(i,k,jc)+ztodt*qf3m(i-nxpt,k,jc)
#else
            qout(i,k,jc) = qout(i,k,jc)+alpha(jc)*qout(i,k,jc)*etamid(k) &
                 *ABS(qout(i,k,jc)-qin(i,k,jc))
#endif
          END DO
        END DO
      END IF
      
      ! Fill up negative values
      
      DO k = 1, plev
        DO i = istart, istop
          qout(i,k,jc) = MAX(qout(i,k,jc), 0.0_dp)
        END DO
      END DO
    END DO
    
  END SUBROUTINE fixer

  SUBROUTINE zonal_mass_integral(cwava,w,q3,pdel,hw1lat,jlat)

    ! Description:
    !
    ! Calculate contribution of current latitude to mass of constituents
    ! being advected by slt.
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          P. Rasch, D. Williamson, August 1992
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    ! parallel version:  T. Diehl, DKRZ, July 1999
    !                    A. Rhodin, MPI, Sept 1999
    !
    ! for more details see file AUTHORS
    !



    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: cwava                ! normalization factor l/(g*plon)
    REAL(dp), INTENT (in) :: w                    ! gaussian weight this latitude
    INTEGER, INTENT (in) :: jlat ! local row index N->S

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: pdel(plond,plev)     ! pressure diff between interfaces
    REAL(dp), INTENT (in) :: q3(plond,plev,pcnst) ! constituents

    !  Array arguments with intent(Out):
    REAL(dp), INTENT (out) :: hw1lat(pcnst)       ! accumulator

    !  Local scalars:
    INTEGER :: m

    !  Executable statements

    ! longitude, level, constituent indices

    DO m = 1, pcnst
      hw1lat(m) = sum_zonal_sl(&
           q3(i1:plon+i1-1,:,m)*pdel(i1:plon+i1-1,:),jlat=jlat)
      !   more efficient but not identical to serial version:
      !   hw1lat(m) = sum_zonal (                                             &
      !     sum ( q3(i1:plon+i1-1,:,m)*pdel(i1:plon+i1-1,:), dim=2), jlat=jlat)
    END DO

    ! The 0.5 factor arises because gaussian weights sum to 2

    DO m = 1, pcnst
      hw1lat(m) = cwava*w*hw1lat(m)*0.5_dp
    END DO

  END SUBROUTINE zonal_mass_integral

  SUBROUTINE global_mass_integral(cwava,w,q1,q2,pdel,etamid,kftype,hwn,jlat)

    ! Description:
    !
    ! Compute contribution of current latitude to global integral of
    ! F*q2*|q2 - q1|**Beta.
    !
    ! Method:
    !
    ! Compute contribution of current latitude to global integral of
    ! F*q2*|q2 - q1|**Beta.
    ! This is a measure of the difference between the fields before and
    ! after the SLT "forecast".
    ! It is used in the "fixer" which enforces conservation in constituent
    ! fields transport via SLT.
    ! Two options are available:
    ! 1. kftype=1 : F=1.  and Beta=1.5
    ! 2. kftype=2 : F=eta and Beta=1.
    !
    ! Reference Rasch and Williamson, 1991
    !
    ! Authors:
    !
    ! Original version:  J. Olson
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          P. Rasch, D. Williamson, August 1992
    ! Modified:          U. Schlese, April 1995  (optional fixers)
    ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
    ! parallel version:  T. Diehl, DKRZ, July 1999
    !                    A. Rhodin, MPI, Sept 1999
    !
    ! for more details see file AUTHORS
    !

    !  Scalar arguments with intent(In):
    REAL(dp), INTENT (in) :: cwava, w
    INTEGER, INTENT (in) :: jlat ! local latitude index N->S

    !  Array arguments with intent(In):
    REAL(dp), INTENT (in) :: etamid(plev), pdel(plond,plev), q1(plond,plev,pcnst), &
         &      q2(plond,plev,pcnst)
    INTEGER :: kftype(pcnst)

    !  Array arguments with intent(InOut):
    REAL(dp), INTENT (inout) :: hwn(pcnst)

    ! cwava   l/(g*plon)
    ! w       Gaussian weight.
    ! q1      Untransported q-field.
    ! q2      Transported   q-field.
    ! pdel    array of pressure differences between layer interfaces
    !         (used for mass weighting ).
    ! hwn     Mass averaged constituent in units of kg/m**2.

    !  Local scalars:
    REAL(dp) :: hwava
    INTEGER :: k, m

    !  Local arrays:
    REAL(dp) :: zhwava (1:plon,plev)

    !  Intrinsic functions
    INTRINSIC ABS, SQRT

    !  Executable statements

    ! accumulator

    DO m = 1, pcnst

      IF (kftype(m)==1) THEN
        hwava = sum_zonal_sl(&
             (q2(i1:plon+i1-1,:,m)*(SQRT(ABS(q1(i1:plon+i1-1,:,m)-&
             q2(i1:plon+i1-1,:,m))))**3)*pdel(i1:plon+i1-1,:), jlat)

        !      more efficient but not identical to serial version:
        !      zhwava = 0.
        !      DO k = 1, plev
        !        zhwava(:) = zhwava(:) + &
        !          (q2(i1:plon+i1-1,k,m)*(SQRT(ABS(q1(i1:plon+i1-1,k,m)-&
        !           q2(i1:plon+i1-1,k,m))))**3)*pdel(i1:plon+i1-1,k)
        !      END DO
        !      hwava = sum_zonal ( zhwava, jlat )

      ELSE IF (kftype(m)==2) THEN
        zhwava = 0.0_dp
        DO k = 1, plev
          zhwava(:,k)=(q2(i1:plon+i1-1,k,m)*etamid(k)*ABS(q1(i1:plon+i1-1,k,m)-&
               q2(i1:plon+i1-1,k,m)))*pdel(i1:plon+i1-1,k)
        END DO
        hwava = sum_zonal_sl(zhwava, jlat )

        !      more efficient but not identical to serial version:
        !      zhwava = 0.
        !      DO k = 1, plev
        !        zhwava(:) = zhwava(:) + &
        !          (q2(i1:plon+i1-1,k,m)*etamid(k)*ABS(q1(i1:plon+i1-1,k,m)-&
        !           q2(i1:plon+i1-1,k,m)))*pdel(i1:plon+i1-1,k)
        !      END DO
        !      hwava = sum_zonal ( zhwava, jlat )

      END IF

      ! The 0.5 factor arises because gaussian weights sum to 2

      hwn(m) = hwn(m) + cwava*w*hwava*0.5_dp
    END DO

  END SUBROUTINE global_mass_integral

END MODULE mo_semi_lagrangian
