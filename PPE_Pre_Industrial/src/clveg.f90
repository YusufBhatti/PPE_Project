#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE clveg(krow)

  ! Description:
  !
  ! Passes climate leaf area index and climate vegetation ratio 
  ! to atmosphere
  !
  ! Method:
  !
  ! This subroutine calculates the leaf area index and vegetation 
  ! ratio for each time step and updates vlt and vgrat.
  !
  ! *clveg* is called from *gpc*.
  !
  ! Authors: 
  !
  ! U. Schulzweida, MPI, July 1999
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_memory_g3b,    ONLY: vlt, vgrat
  USE mo_radiation_parameters,     ONLY: nmonth
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_interpo,       ONLY: wgt1, wgt2, nmw1, nmw2
  USE mo_clim,          ONLY: vltclim, vgratclim

  IMPLICIT NONE

  INTEGER :: krow

  !  Local scalars: 
  INTEGER :: im, jrow, jn, nproma


  !  Executable Statements 

  jrow  = krow

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

  ! Update leaf area index

  ! Annual cycle

  IF (nmonth == 0) THEN

    ! Interpolation in time

    DO jn = 1, nproma
      vlt(jn,jrow) = wgt1*vltclim(jn,jrow,nmw1) + wgt2*vltclim(jn,jrow,nmw2)
      vgrat(jn,jrow) = wgt1*vgratclim(jn,jrow,nmw1) + wgt2*vgratclim(jn,jrow,nmw2)
    END DO

  ELSE

    ! Perpetual month

    im = nmonth
    DO jn = 1, nproma
      vlt(jn,jrow) = vltclim(jn,jrow,im)
      vgrat(jn,jrow) = vgratclim(jn,jrow,im)
    END DO

  END IF

  RETURN
END SUBROUTINE clveg
