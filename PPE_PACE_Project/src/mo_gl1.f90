!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_gl1

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

! global buffered variables

  REAL(dp), TARGET, ALLOCATABLE :: q(:,:)
  REAL(dp), TARGET, ALLOCATABLE :: x(:,:)
  REAL(dp), TARGET, ALLOCATABLE :: xt(:,:,:)

END MODULE mo_gl1
