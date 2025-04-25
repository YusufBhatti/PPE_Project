!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!!
!! These routines belong to the standalone Version of JSBACH
!!
SUBROUTINE jsbalone_iniphy

  ! Description:
  !
  ! Initialises physical constants of uncertain value.
  !
  ! Method:
  !
  ! This routine sets the values for the physical constants used
  ! in the parameterization routines (except for the radiation
  ! black-box) whenever these values are not well enough known to
  ! forbid any tuning or whenever they are subject to an arbitrary
  ! choice of the modeller. These constants will be in *comph2*.
  !
  ! *iniphy* is called from *physc*.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, December 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,         ONLY: dp
  USE mo_physc2,       ONLY: clam, ckap, cb, cc, cd, cchar, cfreec     &
                           , cgam, cvdifts, ctfreez, cz0ice, cqsncr    &
                           , csncri, cwlmax

  IMPLICIT NONE

  !  Executable statements 

!-- 1. Setting of constants

!-- 1.1 Constants for *vdiff*

  clam = 150._dp
  ckap = 0.4_dp
  cb = 5._dp
  cc = 5._dp
  cd = 5._dp
  cchar = 0.018_dp
  cfreec = 0.001_dp
  cgam = 1.25_dp
  cvdifts = 1.5_dp

!-- 1.2 Constants for *vdiff*, *clsst* and *atmice*

  ctfreez = 271.38_dp
  cz0ice = 0.001_dp

!-- 1.3 Constants for *vdiff* and *pysc*

  cqsncr = 0.95_dp
  csncri = 5.85036E-3_dp

!-- 1.7 Constants for *vdiff* and *surf*

  cwlmax = 2.E-4_dp

  RETURN

END SUBROUTINE jsbalone_iniphy
