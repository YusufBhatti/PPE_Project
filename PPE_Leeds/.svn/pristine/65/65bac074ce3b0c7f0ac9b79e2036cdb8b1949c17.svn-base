!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
  SUBROUTINE ml_ocean ( kproma                                         &
                  , pslf,       pslm                                   &
                  , lonorth                                            &
                  , pseaice,    psiced,    palake                      &
                  , ptsi,       ptsw                                   &
                  , pahflw,     pahfsw,    pfluxres                    &
                  , ptrflw,     psoflw,    pqres                       &
                  , pamlcorr,   pamlcorac, pamlheatac                  &
                  , pevapi,     psni,      pcvsi                       &
                  , pfri                      )  

!
!  ---------------------------------------------------------------------
!
  USE mo_kind,           ONLY: wp
  USE mo_physical_constants, ONLY: alf, rhoh2o, alice, alsno, rhosea,       &
                                   rhosno, cpsea, rhoice, dice, rhoilf
  USE mo_physc2,         ONLY: ctfreez
  USE mo_time_control,   ONLY: delta_time
  USE mo_control,        ONLY: lfractional_mask, lmlo_ice, rmlo_depth
!
! Arguments
!
  REAL(wp) ::                                                                &
          pseaice(kproma),     psiced(kproma),     palake(kproma)            &
         ,ptsi(kproma),        ptsw(kproma)                                  &
         ,pahflw(kproma),      pahfsw(kproma),     pfluxres(kproma)          &
         ,ptrflw(kproma),      psoflw(kproma),     pqres(kproma)             &
         ,pamlcorr(kproma),    pamlcorac(kproma),  pamlheatac(kproma)        &
         ,pevapi(kproma),      psni(kproma),       pcvsi(kproma)             &
         ,pfri(kproma),        pslf(kproma),       pslm(kproma)
!
  LOGICAL :: lonorth(kproma)      ! .true. for northern latitude

! Local variables

  REAL(wp) :: zdtime, zfw, zheat, zfbase, ziscond, zdtrilf, zfreez,          &
              zmixcap, zmcapdt, zmcaprilf

! Executable statements
!
! 1. Set up constants
!
  zdtime    = delta_time
  ziscond   = alice/alsno*rhoh2o/rhosno
  zdtrilf   = zdtime/rhoilf
  zfreez    = -dice/zdtrilf
  zmixcap   = rhosea*cpsea*rmlo_depth
  zmcapdt   = zdtime/zmixcap
  zmcaprilf = zmixcap/rhoilf
!
! 2. Mixed layer ocean temperature and ice thickness   
!
  DO jl=1,kproma
!
  IF ((lfractional_mask .AND. palake(jl) .EQ. 0.0_wp .AND. pslf(jl) .LT. 1.0_wp) .OR.  &
   (.NOT. lfractional_mask .AND. palake(jl) .LT. 0.5_wp .AND. pslm(jl) .LT. 0.5_wp)) THEN      ! no lake points
!
     IF (pseaice(jl) .LT. 0.5_wp .OR. .NOT. lmlo_ice) THEN      ! open water
                                                                ! only calculate if sea ice formation is
                                                                ! not turned off.
!
        zfluxw             = pahflw(jl)+pahfsw(jl)+ptrflw(jl)+psoflw(jl)-pamlcorr(jl)
!
!       Water temperature (ptsw)
!
        zts                = ptsw(jl)+zmcapdt*(zfluxw+pfluxres(jl))
        IF(zts.LT.ctfreez) THEN
          pamlcorr(jl)     = MIN(0.0_wp,pamlcorr(jl))
          zfluxw           = pahflw(jl)+pahfsw(jl)+ptrflw(jl)+psoflw(jl)-pamlcorr(jl)
          zts              = ptsw(jl)+zmcapdt*(zfluxw+pfluxres(jl))
        END IF
        ptsi(jl)           = ctfreez
        pfluxres(jl)       = 0._wp
        psiced(jl)         = 0._wp
        IF (zts.GE.ctfreez .OR. .NOT. lmlo_ice) THEN              ! open water (unchanged)
                                                                  ! only calculate if sea ice formation is
                                                                  ! not turned off.
           ptsw(jl)        = zts
        ELSE                                                      ! check ice formation
           ptsw(jl)        = ctfreez
           zfres           = (zts-ctfreez)/zmcapdt                ! < 0.
           IF (zfres.LE.zfreez) THEN                              ! ice formation
              psiced(jl)   = zmcaprilf*(ctfreez-zts)              ! > dice
              pseaice(jl)  = 1._wp
           ELSE
              pfluxres(jl) = zfres
           END IF
        END IF
!  ---------------------------------------------------------------------
     ELSE IF (psiced(jl) .GE. dice) THEN
!
        IF (lonorth(jl)) THEN
           zfbase       =  4.0_wp
        ELSE
           zfbase       = 15.0_wp
        END IF
!
!       Ice thickness (psiced)
!
        zconhflx        = alice*(ptsi(jl)-ctfreez)/(psiced(jl)+ziscond*psni(jl))
        zsubice         = (1._wp-pcvsi(jl))*pevapi(jl)*zdtime/rhoice
        pamlcorr(jl)    = MIN(0.0_wp,pamlcorr(jl))-zfbase
        zhi             = psiced(jl)-zdtrilf*(zconhflx+pqres(jl)+pfluxres(jl)-pamlcorr(jl))+zsubice
        ptsw(jl)        = ctfreez
        IF (zhi .GE. dice) THEN
           psiced(jl)   = zhi
           pseaice(jl)  = 1._wp
           pfluxres(jl) = 0._wp
        ELSE IF (zhi.LE.0._wp) THEN                     ! complete melting
           ptsw(jl)     = ctfreez-zhi/zmcaprilf         ! ptsw > ctfreez
           pfluxres(jl) = -rhoh2o*alf*psni(jl)/zdtime
           psiced(jl)   = 0._wp
           psni(jl)     = 0._wp
           pseaice(jl)  = 0._wp
        ELSE                                            ! incomplete melting
           psiced(jl)   = dice
           pseaice(jl)  = 1._wp
           pfluxres(jl) = (dice-zhi)/zdtrilf
        END IF
     END IF
     IF (lfractional_mask) THEN
        zfw                = 1._wp-pfri(jl)-pslf(jl)
     ELSE
        zfw                = 1._wp-pfri(jl)-pslm(jl)
     END IF
!
!    accumulate mixed layer variables
!
     zheat              = zmixcap*zfw*ptsw(jl)-rhoilf*pfri(jl)*psiced(jl)
     pamlheatac(jl)     = pamlheatac(jl)+zheat*zdtime
     pamlcorac(jl)      = pamlcorac(jl)+pamlcorr(jl)*zfw*zdtime
   END IF
  END DO
!  ---------------------------------------------------------------------
     RETURN
  END SUBROUTINE ml_ocean
