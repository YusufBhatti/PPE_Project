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
SUBROUTINE presf(pf,kdimp,ph,klen)

  ! Description:
  !
  ! Compute full-level pressures from half-level values.
  !
  ! Method:
  !
  ! Full-level pressures are defined as the arithmetic
  ! average of the two adjoining half-level pressures.
  !
  ! *presf* is called from *physc*. 
  ! Parameters are:
  !    *pf*        *computed full-level pressures.
  !    *ph*        *half-level pressures.
  !    *kdimp*     *first dimension of 2-d arrays *pf* and *ph*
  !    *klen*      *number of points for which computation is
  !                 performed.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, May 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,    ONLY: wp
  USE mo_control, ONLY: nlev

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER ,INTENT(in) :: kdimp, klen

  !  Array arguments 
  REAL(wp),INTENT(in) :: ph(kdimp, *)
  REAL(wp),INTENT(inout):: pf(kdimp, *)

  !  Local scalars: 
  INTEGER :: jlev, jlon


  !  Executable statements 

!-- 1. Compute full-level pressure values

  DO jlev = 1, nlev
    DO jlon = 1, klen
      pf(jlon,jlev) = (ph(jlon,jlev)+ph(jlon,jlev+1))*.5_wp
    END DO
  END DO

  RETURN
END SUBROUTINE presf
