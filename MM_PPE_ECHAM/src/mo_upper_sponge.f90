!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_upper_sponge

  !============================================================================
  !
  !- Description:
  !
  !  This module contains:
  !  1) the upper sponge routine
  !
  !   Note: the sponge parameters are part of the DYNCTL namelist
  ! 
  !============================================================================

  USE mo_kind, ONLY: dp
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: spdrag    ! upper sponge layer coefficient (sec)-1      setdyn
  PUBLIC :: enspodi   ! factor by which   upper sponge layer        setdyn
  !                     coefficient is increased from one
  !                     level to next level above.
  PUBLIC :: nlvspd1   ! last (uppermost) layer of upper sponge      setdyn
  PUBLIC :: nlvspd2   ! first (lowest) layer of upper sponge        setdyn
  PUBLIC :: uspnge    ! subroutine - relaxation on zonal waves      stepon

  REAL(dp) :: spdrag, enspodi 

  INTEGER  :: nlvspd1, nlvspd2

  
CONTAINS

  SUBROUTINE uspnge
    !
    !****   *uspnge* upper sponge for divergence and vorticity 
    !
    !   E. Manzini, MPI-HH, 1995, original version
    !   T. Diehl, DKRZ, July 1999, parallel version
    !   I. Kirchner, MPI, December, 2000, time control
    !
    !
    !   purpose
    !   ------
    !   
    !   sponge layer for vorticity and divergence:
    !   upper layer(s) linear relaxation on zonal waves 
    !
    !**   interface
    !     ---------
    !   *uspnge* is called from *stepon*
    !

    USE mo_decomposition, ONLY: lc => local_decomposition
    USE mo_control,       ONLY: nlev
    USE mo_memory_sp,     ONLY: sd, stp, svo
    USE mo_time_control,  ONLY: time_step_len

    REAL(dp) :: zlf(nlev), zs(nlev,2)
    REAL(dp) :: zspdrag, znul
    INTEGER :: jk, is, jlev, jr, snsp
    INTEGER :: mymsp(lc%snsp)

    snsp = lc%snsp
    mymsp = lc%mymsp

    ! 1.  compute linear relaxation
    !     ------- ------ ----------
    zspdrag = spdrag

    znul = 0.0_dp
    DO jk = 1,nlev
       zlf(jk)=znul
    END DO

    zlf(nlvspd2) = zspdrag

    DO jk = nlvspd2-1,nlvspd1,-1
       zspdrag = zspdrag*enspodi
       zlf(jk) = zspdrag
    END DO

    DO jk = 1,nlev
       zs(jk,1) = 1.0_dp/(1.0_dp+zlf(jk)*time_step_len)
       zs(jk,2) = 1.0_dp/(1.0_dp+zlf(jk)*time_step_len)
    END DO

    ! 3.  modify divergence and vorticity         
    !     ------ ---------- --- ---------         

    DO is = 1, snsp
       IF (mymsp(is) /= 0) THEN
          !          DO jlev = 1, nlev                                               
          sd (:,:,is)  = sd(:,:,is)*zs(:,:)
          svo(:,:,is) = svo(:,:,is)*zs(:,:)                   
          !          END DO

          DO jr = 1, 2
             DO jlev = 1, nlev
                stp(jlev,jr,is) = stp(jlev,jr,is)*zs(jlev,1)
             END DO
          END DO

       END IF
    END DO

  END SUBROUTINE uspnge
  
END MODULE mo_upper_sponge
