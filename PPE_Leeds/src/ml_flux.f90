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
SUBROUTINE ml_flux(krow)

  ! Description:
  !
  ! Passes flux corrections to mixed layer ocean
  !
  ! Method:
  !
  ! This subroutine interpolates the flux corrections ...
  ! preliminary....dummy routine
  !
  ! *ml_flux* is called from *gpc*.
  !
  ! Authors: 
  !
  ! U. Schlese, DKRZ, January 1993, original source (clsst)
  ! M. Esch, MPI, September 2002, rearranged for mixed layer ocean's
  !                               flux correction
  ! for more details see file AUTHORS
  ! 
  !

  USE mo_kind,          ONLY: dp
  USE mo_memory_g3b,    ONLY: slf, amlcorr, alake
  USE mo_sst,           ONLY: sst, aice, aflux
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2

  IMPLICIT NONE

  INTEGER :: krow

  ! Local variables

  INTEGER :: jl      ! longitude loop index
  INTEGER :: jrow    ! local latitude index
  INTEGER :: nproma  ! number of longitudes on PE
  REAL(dp):: zdmix, zcpwater, zrho_sea, zmixcap, zmonlen, zic

!  Executable statements
!
! 1. Set up constants
!
  zdmix=50._dp
  zcpwater=3994._dp
  zrho_sea=1025._dp
  zmixcap=zrho_sea*zcpwater*zdmix
  zmonlen=2592000._dp

  jrow    = krow        ! local latitude index

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

!-- 2. Update mixed layer correction


  DO jl=1,nproma

    IF(alake(jl,jrow).EQ.0._dp .AND. slf(jl,jrow).LT.1._dp) THEN
!
      IF(nmw2.LT.nmw1) THEN       ! first part of month
        amlcorr(jl,jrow)=wgt1*aflux(jl,jrow,nmw1)+wgt2*aflux(jl,jrow,nmw2) &
                    -zmixcap*(sst(jl,jrow,nmw1)-sst(jl,jrow,nmw2))         &
                    /zmonlen
      ELSE                        !second half of month
        amlcorr(jl,jrow)=wgt1*aflux(jl,jrow,nmw1)+wgt2*aflux(jl,jrow,nmw2) &
                    -zmixcap*(sst(jl,jrow,nmw2)-sst(jl,jrow,nmw1))         &
                    /zmonlen
      END IF
!
! no flux correction on climatological ice points
!
      zic=(wgt1*aice(jl,jrow,nmw1)+wgt2*aice(jl,jrow,nmw2))*0.01_dp
      zic=MAX(0._dp,MIN(1._dp,zic))
      IF(zic.GT.0.9_dp) amlcorr(jl,jrow)=0._dp
!
    ELSE                          ! land or lake
        amlcorr(jl,jrow)=0._dp
    END IF
  END DO

  RETURN
END SUBROUTINE ml_flux
