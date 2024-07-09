!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_tmp_buffer

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! Matrix for Helmholtz-equation

  REAL(dp), ALLOCATABLE :: cn(:,:,:)   ! setdyn.f  (nlev,nlev,nkp1) -

CONTAINS
  !------------------------------------------------------------------------------
  SUBROUTINE cleanup_tmp_buffer
    !
    ! deallocate module variables
    !
    IF (ALLOCATED(cn)) DEALLOCATE (cn)
  END SUBROUTINE cleanup_tmp_buffer
  !------------------------------------------------------------------------------
END MODULE mo_tmp_buffer
