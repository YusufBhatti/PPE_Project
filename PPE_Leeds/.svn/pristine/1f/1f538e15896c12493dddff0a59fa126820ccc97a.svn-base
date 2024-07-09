!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE update_surfacetemp(klon, pcp,                             &
  &            pescoe, pfscoe, peqcoe, pfqcoe,                       &
  &            psold, pqsold, pdqsold,                               &
  &            pnetrad, pgrdfl,                                      &
  &            pcfh, pcair, pcsat, pfracsu, pgrdcap,                 &
  &            psnew)


  USE mo_time_control,          ONLY: time_step_len
  USE mo_physc2,                ONLY: cvdifts
  USE mo_jsbach_constants,      ONLY: LatentHeatVaporization, LatentHeatSublimation, &
                                      StefanBoltzmann, cemiss
  USE mo_kind,                  ONLY: dp
  IMPLICIT NONE

  INTEGER,  INTENT(in)    :: klon
  REAL(dp),     INTENT(in)    :: pcp(klon), pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
  REAL(dp),     INTENT(in)    :: psold(klon), pqsold(klon), pdqsold(klon)
  REAL(dp),     INTENT(in)    :: pnetrad(klon), pgrdfl(klon)
  REAL(dp),     INTENT(in)    :: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
  REAL(dp),     INTENT(in)    :: pgrdcap(klon)
  REAL(dp),     INTENT(out)   :: psnew(klon)
  REAL(dp) :: zcolin(klon), zcohfl(klon), zcoind(klon), zicp(klon), zca(klon), zcs(klon)
  REAL(dp) :: pdt, pemi, pboltz
  REAL(dp) :: platev, platsu
  REAL(dp) :: ztpfac1, ztmst

!-------------------------------------------------------------------------------------
! Constants

  ztpfac1 = cvdifts
  ztmst   = time_step_len
  pdt     = ztpfac1*ztmst                  ! zcons29 in 'old' vdiff
  pemi    = cemiss                         ! emissivity
  pboltz  = StefanBoltzmann
  platev  = LatentHeatVaporization
  platsu  = LatentHeatSublimation

!************************************************************************************
!
     zicp(:) = 1._dp / pcp(:)
!
     zca(:)    = platsu * pfracsu(:) +  platev * (pcair(:) - pfracsu(:))
     zcs(:)    = platsu * pfracsu(:) +  platev * (pcsat(:) - pfracsu(:))
!
     zcolin(:) = pgrdcap(:)*zicp(:) +                                                     &
                 pdt * (zicp(:) * 4._dp * pemi * pboltz * ((zicp(:) * psold(:))**3) -     &
                 pcfh(:) * (zca(:) * peqcoe(:) - zcs(:)) * zicp(:) * pdqsold(:))
!
     zcohfl(:) = -pdt * pcfh(:) * (pescoe(:) - 1._dp)
!
     zcoind(:) = pdt * (pnetrad(:) + pcfh(:) * pfscoe(:) +  pcfh(:) *                      &
                 ((zca(:) * peqcoe(:) - zcs(:)) * pqsold(:) + zca(:) * pfqcoe(:)) + pgrdfl(:))
!
    psnew(:)  = (zcolin(:) * psold(:) + zcoind(:)) / (zcolin(:) + zcohfl(:))

END SUBROUTINE update_surfacetemp
