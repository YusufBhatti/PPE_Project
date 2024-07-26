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
SUBROUTINE intaero(krow)

  ! Description:
  !
  ! Passes aerosoles to atmosphere
  !
  ! Method:
  !
  ! This subroutine interpolates in time.
  !
  ! *intaero* is called from *gpc*.
  !
  ! Authors: 
  !
  ! U. Schlese, DKRZ, U. Schulzweida, M. Esch, MPI, Apr 2004, original source
  ! for more details see file AUTHORS
  !

  USE mo_memory_g3b,    ONLY: so4nat, so4all
  USE mo_so4,           ONLY: so4nat_io, so4all_io
  USE mo_control,       ONLY: nlev
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2

  IMPLICIT NONE

  INTEGER :: krow

  ! Local variables

  INTEGER :: jl      ! longitude loop index
  INTEGER :: jk      ! level loop index
  INTEGER :: jrow    ! local latitude index
  INTEGER :: nproma  ! number of longitudes on PE

  !  Executable statements

  jrow    = krow        ! local latitude index

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

!-- 1. Update temperatures and sea ice

  DO jk=1,nlev
    DO jl=1,nproma
!
!-- 1.1 Annual cycle
!
     so4nat(jl,jk,jrow)=wgt1*so4nat_io(jl,jk,jrow,nmw1)+wgt2*so4nat_io(jl,jk,jrow,nmw2)
     so4all(jl,jk,jrow)=wgt1*so4all_io(jl,jk,jrow,nmw1)+wgt2*so4all_io(jl,jk,jrow,nmw2)
!
    END DO
  END DO

  RETURN
END SUBROUTINE intaero
