!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_climatology

!------------------------------------------------------------------------------
! Calculate instantaneous values from a monthly climatology
!
! This routine is based on 'update_lai', it just uses more general variable
! names
!------------------------------------------------------------------------------
  USE mo_interpo, ONLY: wgt1, wgt2, nmw1, nmw2
  USE mo_radiation_parameters, ONLY: nmonth
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PUBLIC :: get_from_climatology

  INTERFACE get_from_climatology
    MODULE PROCEDURE get_from_clim_tiles ! arrays with tiles dimension
    MODULE PROCEDURE get_from_clim       ! arrays without tiles dimension 
  END INTERFACE

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE get_from_clim_tiles(kdim, ntiles, clim, instant)
!------------------------------------------------------------------------------
  INTEGER,  INTENT(in)  :: kdim, ntiles
  REAL(dp), INTENT(in)  :: clim(kdim,ntiles,0:13)
  REAL(dp), INTENT(out) :: instant(kdim,ntiles)

  IF (nmonth == 0) THEN
     ! Interpolation in time
     instant(:,:) = wgt1*clim(:,:,nmw1) + wgt2*clim(:,:,nmw2)
  ELSE
     ! Perpetual month
     instant(:,:) = clim(:,:,nmonth)
  END IF

  END SUBROUTINE get_from_clim_tiles

!------------------------------------------------------------------------------
  SUBROUTINE get_from_clim (kdim, ntiles, clim, instant)
!------------------------------------------------------------------------------
  INTEGER,  INTENT(in)  :: kdim, ntiles
  REAL(dp), INTENT(in)  :: clim(kdim,0:13)
  REAL(dp), INTENT(out) :: instant(kdim,ntiles)

  REAL(dp)              :: inst_1d(kdim)


  IF (nmonth == 0) THEN
     ! Interpolation in time
     inst_1d(:) = wgt1*clim(:,nmw1) + wgt2*clim(:,nmw2)
  ELSE
     ! Perpetual month
     inst_1d(:) = clim(:,nmonth)
  END IF
  instant(:,:) = SPREAD(inst_1d,DIM=2,NCOPIES=ntiles)

  END SUBROUTINE get_from_clim

END MODULE mo_climatology
