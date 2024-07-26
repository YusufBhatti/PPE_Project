!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE diag_water_balance(kbdim, kdeep, ntiles,  &
        zdtime, lmask, lpglac, fract, jllog, jlpoint,  &
        rain, snow, prunoff, pdrain, pevapl,  &
        psn, pws,pwl,pwl_sn, wsi,  &
        etrans, pevwsd, epot, zevwld, esnow,  &
        tsurf, redevap )
 
! Code converted using TO_F90 by Alan Miller
! Date: 2011-05-19  Time: 13:39:22

!     ******* Balancing of the water balance per time step and write Log output
!             in a Table - equivalent to XB2TAB from JSBACH_boxwsi.com

!     *** Vs. 1.0 -- Feb. 2008 -- Hag
!     ****REMO-Fields: 012 013 053 160 260 182 001 002 003 004 005 140
!     *** JSBACH-Version -> Fields: 76 67 64 63 260  44 59 55 57 77

!     *** Vs. 1.1 -- july 2009 -- Hag
!     *** Blizzard Version
!
!     *** Vs. 1.2 -- August 2009 -- Hag
!     *** Support for more grid cells on the same processor by output
!
!     ***  zdtime        Time step in seconds
!     ***  Rain and snow [kg/m**2]
!     ***  pwl           Water content [m] in skin reservoir
!     ***  pwl_sn        Snow content [m] on the canopy (is not included in skin reservoir!)
!     ***  pws           Soil water content [m]
!     ***  pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
!     ***  prunoff       Surface runoff [m water equivalent] at non-glacier points (accumul.)
!     ***  pdrain        Drainage at non-glacier points [m water equivalent] (accumul.)
!     ***  pevwsd        Evapotranspiration without sublimation and evaporation
!                        from interception reservoir [m]
!     ***  zevwld        Evaporation from the skin reservoir [m]
!     ***  esnow         Evaporation over snow [kg/m**2/s]
!     ***  soil%transpiration    Transpiration (inst.)' [kg/m**2/s]
!     ***  tsurf         Temperature of the surface
!     ***  redevap       DiagnostIc evaporation reduction obtained from 5 layer scheme [kg/m**2/s]
!
!     *** jllog          relativer Index der betreffenden Gitterbox zwischen 1 und kbdim
!     *** jlpoint        Absoluter Index aus kpoints of water balance test grid box --> Log-Filename
!
!     *** ngrid          Number of different grids for log output (ngridmax = 10)
!     *** mbox(:)        Array with different gridboxes indices (absolute)
!     *** mlu(:)         Array with Logical unit for the different grid boxes
!
USE mo_utils,             ONLY: average_tiles
USE mo_exception,         ONLY: message, message_text
USE mo_mpi,               ONLY: p_pe, p_io
USE mo_kind,              ONLY: dp
USE mo_jsbach_constants,  ONLY: Tmelt
!

INTEGER,  INTENT(IN)                    :: kbdim
INTEGER,  INTENT(IN)                    :: kdeep
INTEGER,  INTENT(IN)                    :: ntiles
REAL(dp), INTENT(IN)                    :: zdtime
LOGICAL,  INTENT(IN)                    :: lmask(kbdim,ntiles)
LOGICAL,  INTENT(IN)                    :: lpglac(kbdim,ntiles)
REAL(dp), INTENT(IN)                    :: fract(kbdim,ntiles)
INTEGER,  INTENT(IN)                    :: jllog
INTEGER,  INTENT(IN)                    :: jlpoint
REAL(dp), INTENT(IN)                    :: rain(kbdim)
REAL(dp), INTENT(IN)                    :: snow(kbdim)
REAL(dp), INTENT(IN)                    :: prunoff(kbdim, ntiles)
REAL(dp), INTENT(IN)                    :: pdrain(kbdim,ntiles)
REAL(dp), INTENT(IN)                    :: pevapl(kbdim,ntiles)
REAL(dp), INTENT(IN)                    :: psn(kbdim)
REAL(dp), INTENT(IN)                    :: pws(kbdim)
REAL(dp), INTENT(IN)                    :: pwl(kbdim,ntiles)
REAL(dp), INTENT(IN)                    :: pwl_sn(kbdim,ntiles)
REAL(dp), INTENT(IN)                    :: wsi(kbdim, kdeep)
REAL(dp), INTENT(IN)                    :: etrans(kbdim,ntiles)
REAL(dp), INTENT(IN)                    :: pevwsd(kbdim)
REAL(dp), INTENT(IN)                    :: epot(kbdim,ntiles)
REAL(dp), INTENT(IN)                    :: zevwld(kbdim)
REAL(dp), INTENT(IN)                    :: esnow(kbdim)
REAL(dp), INTENT(IN)                    :: tsurf(kbdim,ntiles)
REAL(dp), INTENT(IN)                    :: redevap(kbdim, ntiles)

INTEGER, PARAMETER :: narr=13
INTEGER, PARAMETER :: npl=5
INTEGER, PARAMETER :: ngridmax=10

CHARACTER (LEN=6)  :: c6
CHARACTER (LEN=80) :: dnout
REAL(dp) ::  fdat(narr+npl)       &  ! hydrological values for gridbox jlpoint
           , fevap(8)                ! Evap. fluxes

REAL(dp), SAVE :: fdatm1(ngridmax, narr+npl) ! ..M1 = last time step for grid cell jlpoint

REAL(dp) :: ufak

INTEGER  :: j, luout, jgrid
REAL(dp) :: zdum(kbdim)

INTEGER, SAVE :: ilauf = 0
INTEGER, SAVE :: ngrid = 0
INTEGER, SAVE :: mbox(ngridmax)
INTEGER, SAVE :: mlu(ngridmax)

  ufak = zdtime          ! kg/m**2/s == mm/s --> mm

  !     *** Determine if old or new grid box
  IF (ilauf == 0) THEN
    ngrid = 1
    jgrid = 1
    mbox(ngrid) = jlpoint
    mlu(ngrid)  = 220
  ELSE
    jgrid = 0
    DO j=1,ngrid
      IF (mbox(j) == jlpoint) jgrid = j
    END DO
    IF (jgrid == 0) THEN
      ngrid = ngrid+1
      mbox(ngrid) = jlpoint
      mlu(ngrid) = 218 + 2*ngrid
      jgrid = ngrid
      ilauf = 0
    END IF
  END IF
  luout = mlu(jgrid)

!     *** Unit conversion to mm
  CALL average_tiles(tsurf, lmask, fract, zdum)
  fdat( 1) = zdum(jllog) - Tmelt                  ! K --> deg C
  CALL average_tiles(pdrain, lmask, fract, zdum)
!cc      FDAT(2) = zdum(jllog) * UFAK
  fdat( 2) = zdum(jllog) * 1000._dp
  CALL average_tiles(prunoff, lmask, fract, zdum)
!cc      FDAT(3) = zdum(jllog) * UFAK                     ! 3 = Surface runoff
  fdat( 3) = zdum(jllog) * 1000._dp
  fdat( 4) = (rain(jllog) + snow(jllog) ) * ufak

  CALL average_tiles(pevapl, lmask, fract, zdum)
  fdat( 5) = zdum(jllog) * ufak                  ! * 1000/Rhoh2O = 1. --> mm
  fdat( 6) = psn(jllog) *1000._dp
  CALL average_tiles(pwl+pwl_sn, lmask, fract, zdum)
  fdat( 7) = zdum(jllog) * 1000._dp

  fdat( 8) = pws(jllog) *1000._dp
  fdat( 9) = wsi(jllog,1) *1000._dp
  fdat(10) = wsi(jllog,2) *1000._dp
  fdat(11) = wsi(jllog,3) *1000._dp
  fdat(12) = wsi(jllog,4) *1000._dp
  fdat(13) = wsi(jllog,5) *1000._dp

  CALL average_tiles(etrans, lmask, fract, zdum)
  fevap(1) = zdum(jllog) * ufak                    ! * 1000/Rhoh2O = 1. --> mm
  fevap(2) = pevwsd(jllog) * 1000._dp - fevap(1)   ! Bare Soil Evaporation
  fevap(3) = esnow(jllog) * ufak                   ! Snow evaporation
  fevap(4) = zevwld(jllog) * 1000._dp              ! Skin Reservoir evap.
  fevap(5) = sum(fevap(1:4))
  fevap(6) = fdat(5)
  CALL average_tiles(epot, lmask, fract, zdum)
  fevap(7) = zdum(jllog) * ufak                 ! * 1000/Rhoh2O = 1. --> mm
  CALL average_tiles(redevap, lmask, fract, zdum)
  fevap(8) = zdum(jllog) * ufak                 ! * 1000/Rhoh2O = 1. --> mm

!     *** Sum of 5 Layers and Water Balance Res.
  fdat(narr+1) = sum(fdat(narr-4:narr))

!     *** P-E-R = P - R - D + E
  fdat(narr+2) = fdat(4) - fdat(3) - fdat(2) + fdat(5)

  IF (ilauf == 0) THEN
    WRITE(message_text,*) 'jlpoint =', jlpoint,  &
       & '  jllog = ', jllog, ' bei kbdim = ', kbdim, p_pe, p_io
    CALL message('DIAG_WATER_BALANCE', message_text)

    WRITE(c6, '(I6.6)') jlpoint
    dnout='tab_xblock_' // c6 // '.txt'
    OPEN(luout, FILE=dnout, FORM='formatted')
    WRITE(luout, '(A8,1X, 4(A8,1X),A8,1X, 9(A9,1X),A9,1X,A8,1X,A9 )') 'Temp.',  &
       & 'Drain', 'Surf-Run', 'Prec.', 'Evap.', 'Snowpack', 'WSkin', 'WS buck',  &
       & 'WS 1', 'WS 2', 'WS 3', 'WS 4', 'WS 5',  &
       & 'Sum5WS', 'P-E-R', 'Del-Wges','Del-WSb', 'Del-WSI'

!       *** Nullify due to balancing of P-E-R (initialization value unknown)
    fdat(2:5) = 0._dp
    fdat(narr+2:narr+5) = 0._dp
    ilauf = ilauf + 1

!       *** Evap-File
    dnout = 'tab_evap_' // c6 // '.txt'
    OPEN(luout+1, FILE=dnout, FORM='formatted')
    WRITE(luout+1,'(8(A9,1X),A2,L1,L1,A3,L1,L1)') 'ETrans', 'EBsoil',  &
      &  'ESnow', 'ESkin', 'Sum 4E', 'Evap.', 'Epot', 'RedEvap',  &
      &  'L=',lmask(jllog,1), lmask(jllog,2), ' G=', lpglac(jllog,1), lpglac(jllog,2)

  ELSE
    fdat(narr+3) = fdat(6) - fdatm1(jgrid, 6) + fdat(7) - fdatm1(jgrid, 7) + &
       & fdat(narr+1) - fdatm1(jgrid, narr+1)     ! WSI, Wskin, and Snow changes
!!!     &             FDAT(8) - FDATM1(jgrid, 8)            ! WSb, Wskin, and Snow changes
    fdat(narr+4) = fdat(8) - fdatm1(jgrid, 8)           ! Bucket changes
    fdat(narr+5) = fdat(narr+1) - fdatm1(jgrid, narr+1) ! WSi changes
  END IF

!     *** Write Water Balance Table
  WRITE(luout, '( F8.3,1X,4(F8.5,1X),F9.5,1X,F8.5,1X,5(F9.5,1X),  &
     & 2(F9.4,1X), F9.5,1X,F9.5,1X,F8.4,1X,F9.4 )') (fdat(j),j=1,narr+npl)

!     *** Write Evaporation fluxes Table
  WRITE(luout+1, '(8(F9.6,1X))') (fevap(j),j=1,8)

  fdatm1(jgrid, :) = fdat(:)

  IF (ilauf <= 121) THEN
    ilauf = ilauf+1
    
    WRITE(message_text,*) ' ilauf = ', ilauf-1,  &
       & '  jlpoint =', jlpoint, '  jllog = ', jllog, ' CPU=', p_pe
    CALL message('DIAG_WATER_BALANCE', message_text)
  END IF

END SUBROUTINE diag_water_balance
