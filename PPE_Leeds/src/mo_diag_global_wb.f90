!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_diag_global_wb
IMPLICIT NONE

PRIVATE

PUBLIC diag_global_wb

CONTAINS

SUBROUTINE diag_global_wb(kbdim, ntiles, nsoil,  &
        zdtime, lmask, fract,  &
        rain, snow, prunoff, pdrain, pevapl,  &
        psn, pws, pwl, pwl_sn, wtot_m1,  &
        watbal, wsi )
 
! Code converted using TO_F90 by Alan Miller
! Date: 2011-05-19  Time: 13:39:06

!     ******* Calculates the water budget per time step in JSBACH and save it as
!             fluxes (P-E-Suru-Dr) - storage changes in field watbal

!     *** Vs. 1.0 -- Jan. 2010 -- Hag
!     *** Programming partly analogous to subroutine diag_water_balance.f
!
!     ***  zdtime        Time step in seconds
!     ***  Rain and snow [kg/m**2]
!     ***  pwl           Water content [m] in skin reservoir
!     ***  pwl_sn        Snow content [m] on the canopy (is not included in skin reservoir!)
!     ***  pws           Root zone soil moisture content (bucket for nsoil==1) [m]
!     ***  pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
!     ***  prunoff       Surface runoff [m water equivalent] at non-glacier points (accumul.)
!     ***  pdrain        Drainage at non-glacier points [m water equivalent] (accumul.)
!     ***  wtot_m1       Total water storage (WSI, Wskin, and Snow) of previous time step

USE mo_utils,             ONLY: average_tiles
USE mo_kind,              ONLY: dp

INTEGER,  INTENT(IN)                     :: kbdim
INTEGER,  INTENT(IN)                     :: ntiles
INTEGER,  INTENT(IN)                     :: nsoil ! StW: Inserted instead of i5layer and kdeep
REAL(dp), INTENT(IN)                     :: zdtime
LOGICAL,  INTENT(IN)                     :: lmask(kbdim,ntiles)
REAL(dp), INTENT(IN)                     :: fract(kbdim,ntiles)
REAL(dp), INTENT(IN)                     :: rain(kbdim)
REAL(dp), INTENT(IN)                     :: snow(kbdim)
REAL(dp), INTENT(IN)                     :: prunoff(kbdim, ntiles)
REAL(dp), INTENT(IN)                     :: pdrain(kbdim,ntiles)
REAL(dp), INTENT(IN)                     :: pevapl(kbdim,ntiles)
REAL(dp), INTENT(IN)                     :: psn(kbdim)
REAL(dp), INTENT(IN)                     :: pws(kbdim)
REAL(dp), INTENT(IN)                     :: pwl(kbdim,ntiles)
REAL(dp), INTENT(IN)                     :: pwl_sn(kbdim,ntiles)
REAL(dp), INTENT(IN)                     :: wtot_m1(kbdim)
REAL(dp), INTENT(OUT)                    :: watbal(kbdim)
REAL(dp), INTENT(IN), OPTIONAL           :: wsi(kbdim, nsoil)

INTEGER, PARAMETER :: narr=8
REAL(dp) :: fdat(narr, kbdim)          ! hydrological fields

REAL(dp) :: ufak
INTEGER  :: jl
REAL(dp) :: zdum(kbdim)

  ufak = zdtime          ! kg/m**2/s == mm/s --> mm

! *** Filling of variables and unit conversion to mm
  CALL average_tiles(pdrain, lmask, fract, zdum)
  fdat(2,:) = zdum(:) * 1000._dp                ! 2 = Drainage
  CALL average_tiles(prunoff, lmask, fract, zdum)
  fdat(3,:) = zdum(:) * 1000._dp                ! 3 = Surface runoff
  fdat(4,:) = (rain(:) + snow(:) ) * ufak       ! 4 = Total precipitation
  CALL average_tiles(pevapl, lmask, fract, zdum)
  fdat(5,:) = zdum(:) * ufak                    ! 5 = Evap., * 1000/Rhoh2O = 1. --> mm
  !
  fdat(6,:) = psn(:) *1000._dp                  ! 6 = snowpack
  CALL average_tiles(pwl+pwl_sn, lmask, fract, zdum)
  fdat(7,:) = zdum(:) * 1000._dp                ! 7 = skin reservoir

  IF (nsoil == 5) THEN
    DO jl=1, kbdim
      fdat(8,jl) = sum( wsi(jl,1:nsoil) ) *1000._dp  ! 8 = Total 5 layer soil moisture
    END DO
  ELSE
    fdat(8,:) = pws(:) *1000._dp                     ! 8 = Bucket soil moisture
  END IF

! *** P-E-R = P - R - D + E
  fdat(1,:) = fdat(4,:) - fdat(3,:) - fdat(2,:) + fdat(5,:)
  !
  zdum(:) = wtot_m1(:) * 1000._dp
  watbal(:)= fdat(1,:) - (fdat(6,:)+fdat(7,:)+fdat(8,:)-zdum(:))

! *** Unit watbal is mm ~ kg/m^2
! *** Here, this is a value valid for the time step

  RETURN

END SUBROUTINE diag_global_wb
!
END MODULE
