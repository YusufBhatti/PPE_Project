!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_sub_echam
!-----------------------------------------------------------------------
! This module holds the simple submodule provided by echam. It allows to
! request tracers via the namelist and provides routines for commonly
! used processes. All these routines are called from routines within
! call_submodels.f90
!
! Authors: 
! A.Rhodin MPI/DWD April 2002 original code
!------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  USE mo_kind,    ONLY: wp
  USE mo_tracdef, ONLY: trlist
  IMPLICIT NONE
  !----------------
  ! Public entities
  !----------------
  PUBLIC :: radionucl_sink       ! apply exponential decay
  PRIVATE

CONTAINS
!==============================================================================
  !==========
  ! Processes
  !==========
!------------------------------------------------------------------------------
  SUBROUTINE radionucl_sink (kproma, kbdim, klev, pxtm1, pxtte)
    !------------------------------------------------------------
    ! Description:
    !
    ! Calculates the decrease of tracer concentration for a given 
    ! exponential decay time
    !
    ! Method:
    !
    ! The mass mixing-ratio of tracers is multiplied with
    ! exp(time-step/decay-time).
    ! This routine could also be used for emission or sink
    ! above the surface.
    !
    ! *radionucl_sink* is called from *physc*
    !
    ! Authors:
    !
    ! J. Feichter, MI, August 1991, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    !------------------------------------------------------------
    USE mo_time_control,  ONLY: time_step_len
    !----------
    ! Arguments
    !----------
    INTEGER ,INTENT(in)    :: kproma, kbdim, klev             ! dimensions
    REAL(wp),INTENT(in)    :: pxtm1(kbdim,klev,trlist% ntrac) ! concentr. t-dt
    REAL(wp),INTENT(inout) :: pxtte(kbdim,klev,trlist% ntrac) ! tendency
    !----------------
    ! Local variables
    !----------------
    INTEGER :: jt                 ! tracer index
    REAL(wp):: zxtp1(kproma,klev) ! tracer concentration at t+dt
    REAL(wp):: ztmst              ! time step
    REAL(wp):: zqtmst             ! 1 / time step
    !----------------------
    ! Executable statements
    !----------------------
    ztmst  = time_step_len
    zqtmst = 1.0_wp/ztmst

    DO jt = 1, trlist% ntrac
       IF (trlist% ti(jt)% tdecay /= 0._wp) THEN
          zxtp1(:,:)           = pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt) * ztmst
          zxtp1(:,:)           = EXP (-ztmst/trlist% ti(jt)% tdecay) * zxtp1(:,:)
          pxtte(1:kproma,:,jt) = (zxtp1(:,:)-pxtm1(1:kproma,:,jt))*zqtmst
       END IF
    END DO

  END SUBROUTINE radionucl_sink
!==============================================================================
END MODULE mo_sub_echam
