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
SUBROUTINE gpc(krow)

  ! Description:
  !
  ! Grid point computations.
  !
  ! Method:
  !
  ! This subroutine controls parts the computations in grid points, 
  ! that is the physical computations(*phys*) 
  ! and grid point contributions to the semi implicit scheme (*si1*).
  ! 
  !
  ! *gpc* is called from *scan1*.
  ! 
  ! Externals:
  !   *si1*       grid point contributions to the semi implicit scheme.
  !   *physc*     physical computations.
  !
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! F. Lunkeit, MI, June 1989, CLSST added
  ! U. Schlese, MPI, July 1989, add seaice computations
  ! U. Schlese, DKRZ, January 1995, initialization of soil temperatures
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ, and M. Esch, MPI July 1999, modifications for ECHAM5
  ! I. Kirchner, MPI Sepember 2000, nudging sst update
  ! I. Kirchner, MPI December 2000, time control
  ! U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
  ! U. Schlese, MPI, November 2008, "iaero" instead of "lso4"
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,       ONLY: lnudge, lcouple, lmlo
  USE mo_radiation_parameters,     ONLY: iaero
  USE mo_nudging_sst,   ONLY: NudgingSSTnew

  IMPLICIT NONE

  INTEGER  :: krow

  !  Executable statements 

!-- 1. Distribute climate values

  IF (.NOT. lcouple) THEN
     IF (lnudge) THEN
        CALL NudgingSSTnew(krow)
     ELSE
        IF (lmlo) THEN
          CALL ml_flux(krow)
        ELSE
          CALL clsst(krow)
        END IF
     END IF
  END IF

  CALL clveg(krow)

  IF(iaero == 4) CALL intaero(krow)

!-- 2. Parametrisation of diabatic processes

  CALL physc(krow)

END SUBROUTINE gpc
