!>
!! @brief
!!       *cuflx*
!!        This routine does the final calculation of convective fluxes in the
!!        cloud layer and in the subcloud layer
!!        *cudtdq*
!!        This routine updates t and q tendencies, precipitation rates, does
!!        global diagnostics
!!        *cududv*
!!        This routine updates u and v tendencies, does global diagnostic of
!!        dissipation
!!
!! @author M. Tiedtke, ECMWF,    Dec 1989
!!
!! @par Revision History
!! - Taken from ECHAM6.3, unified for ICON/ECHAM by Monika Esch, MPI-M (2015-06)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cufluxdts

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: grav, alv, als, alf, tmelt, cpd, vtmpc2
#ifndef __ICON__
  USE mo_physical_constants, ONLY: rd
  USE mo_physc2,             ONLY: cevapcu
  USE mo_submodel,           ONLY: lanysubmodel, lham    ! ### explicit submodel dependency (zcucov)
  USE mo_submodel_interface, ONLY: cuflx_subm
  USE mo_tracdef,            ONLY : trlist
  USE mo_time_control,       ONLY : delta_time
  USE mo_vphysc,             ONLY: set_vphysc_var
  !++ for cfmip diagnostics
  USE mo_memory_cfdiag,      ONLY : locfdiag, mc, mcu, mcd, imc, smc, dmc
  USE mo_geoloc,             ONLY : gboxarea
  USE mo_cosp_offline,       ONLY: locospoffl, cospoffl_ccrain, cospoffl_ccsnow
  USE mo_memory_g3b,         ONLY : aprc_na
  !++-end-for cfmip diagnostics-------------------------------------------------------
#endif
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cuflx, cudtdq, cududv

!>>SF
#ifdef HAMMOZ
  INTEGER  :: ierr, ibc_detr_cond
#endif
!<<SF

CONTAINS
  !>
  !!
  !!++mgs: zcucov and zdpevap now 1d vectors
  !!
  SUBROUTINE cuflx(    kproma, kbdim, klev, klevp1,                                  &
    &        pqen,     pqsen,    ptenh,    pqenh,                                    &
    &        ktrac,                                                                  &
    &        ptime_step_len,                                                         &
#ifndef __ICON__
    &        krow,                                                                   &
    &        pxtte,    pxtu,     ptu,                                                &
    &        pmwc,     pmrateprecip,                                                 &
#else
    &        cevapcu,                                                                &
#endif
    &        pxtenh,   pmfuxt,   pmfdxt,                                             &
    &        paphp1,   pgeoh,                                                        &
    &        kcbot,    kctop,    kdtop,                                              &
    &        ktype,    lddraf,   ldcum,                                              &
    &        pmfu,     pmfd,     pmfus,    pmfds,                                    &
    &        pmfuq,    pmfdq,    pmful,                                              &
    &        pdmfup,   pdmfdp,   prfl,     prain,                                    &
    &        pcpcu,                                                                  &
    &        pten,     psfl,     pdpmel,   ktopm2                                    &
#ifdef HAMMOZ /* SF */
            ,pxtbound                                                                &
#endif
           )
    !
    INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktrac
    INTEGER, INTENT (OUT):: ktopm2
    REAL(wp),INTENT (IN) :: ptime_step_len
#ifdef __ICON__
    REAL(wp),INTENT (IN) :: cevapcu(klev)
#endif
    REAL(wp):: pqen(kbdim,klev),        pqsen(kbdim,klev),                           &
      &        ptenh(kbdim,klev),       pqenh(kbdim,klev),                           &
      &        paphp1(kbdim,klevp1),    pgeoh(kbdim,klev)
    REAL(wp):: pmfu(kbdim,klev),        pmfd(kbdim,klev),                            &
      &        pmfus(kbdim,klev),       pmfds(kbdim,klev),                           &
      &        pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                           &
      &        pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),                          &
      &        pmful(kbdim,klev),                                                    &
      &        prfl(kbdim),             prain(kbdim)
    REAL(wp):: pten(kbdim,klev),        pdpmel(kbdim,klev),                          &
      &        psfl(kbdim)
    REAL(wp):: pcpcu(kbdim,klev)
    INTEGER :: kcbot(kbdim),            kctop(kbdim),                                &
      &        kdtop(kbdim),            ktype(kbdim)
    LOGICAL :: lddraf(kbdim),           ldcum(kbdim)
    REAL(wp):: pxtenh(kbdim,klev,ktrac),                                             &
      &        pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)
!>>SF
#ifdef HAMMOZ
    REAL(wp), INTENT(inout) :: pxtbound(kbdim,ktrac) !SF boundary condition for xt_conv_massfix
#endif
!<<SF
    !
    INTEGER :: jl, jk, jt, ikb
    !++mgs
    REAL(wp):: zcons1, zcons2,         ztmelp2, zzp, zfac, zsnmlt                    &
      &      , zrfl, zrnew, zrmin, zrfln, zdrfl, zrsum
    REAL(wp):: zpsubcl(kbdim), zcucov(kbdim), zdpevap(kbdim)
    !--mgs
#ifndef __ICON__
    !
    REAL(wp) :: afac(kbdim) ! for cfmip diagnostics
    !
    !---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):---------
    !
    INTEGER, INTENT(IN) :: krow
    REAL(wp):: pxtte(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),                      &
      &        ptu(kbdim,klev)
    REAL(wp):: zmlwc(kbdim,klev),       zmiwc(kbdim,klev),                           &
      &        zmratepr(kbdim,klev),    zmrateps(kbdim,klev),                        &
      &        zfrain(kbdim,klev),      zfsnow(kbdim,klev),                          &
      &        zdpg(kbdim,klev),        zfevapr(kbdim,klev),                         &
      &        zfsubls(kbdim,klev),     zaclc(kbdim,klev),                           &
      &        zmsnowacl(kbdim,klev),   zrhou(kbdim,klev)
    REAL(wp):: pmwc(kbdim,klev),        pmrateprecip(kbdim,klev)
    REAL(wp):: ztc,         zzfac,      zwu
    REAL(wp):: zcaa, &   !(Constants for partitioning of cloud water 
      &        zcab      ! into liquid and solid part)
    REAL(wp):: zalpha    ! Fraction of cloud water in liquid phase
                         ! =>  (1-zalpha)      "   in solid  phase
    REAL(wp):: zpsubcl_sav(kbdim)
    !---End Included for scavenging---------------------------------------------------
#endif
    !
    !*             Specify constants
    !
    zcons1=cpd/(alf*grav*ptime_step_len)
    zcons2=1._wp/(grav*ptime_step_len)
    ztmelp2=tmelt+2._wp
#ifndef __ICON__
    !
    !---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
    !
    IF(lanysubmodel) THEN
      zcaa = 0.0059_wp   ! (Constants for partitioning of cloud water 
      zcab = 0.003102_wp !  into liquid and solid part)
      zmlwc(1:kproma,:)     = 0._wp
      zmiwc(1:kproma,:)     = 0._wp
      zmratepr(1:kproma,:)  = 0._wp
      zmrateps(1:kproma,:)  = 0._wp
      zfrain(1:kproma,:)    = 0._wp
      zfsnow(1:kproma,:)    = 0._wp
      zdpg(1:kproma,:)      = 0._wp
      zfevapr(1:kproma,:)   = 0._wp
      zfsubls(1:kproma,:)   = 0._wp
      zmsnowacl(1:kproma,:) = 0._wp
      zwu            = 2.0_wp
      !--- Set cloud cover to 1 below cloud top and bottom:
      !
      DO jk=1, klev
        DO jl=1, kproma
          IF (jk>=kctop(jl) .AND. jk<=kcbot(jl)) THEN
            zaclc(jl,jk) = 1._wp
          ELSE
            zaclc(jl,jk) = 0._wp
          END IF
        END DO
      END DO
    END IF
    !
    !---End Included for scavenging-----------------------------------------
#endif
    !
    !*    1.0          Determine final convection fluxes
    !                  ---------------------------------
    !
    !  itop=klev
    DO jl=1,kproma
      ! itop=MIN(itop,kctop(jl))
      IF(.NOT.ldcum(jl).OR.kdtop(jl).LT.kctop(jl)) lddraf(jl)=.FALSE.
      IF(.NOT.ldcum(jl)) ktype(jl)=0
    END DO
    ktopm2=1
    DO jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
      DO jl=1,kproma
        IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
          pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*                                     &
            &                  (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
          pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
          IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
            pmfds(jl,jk)=pmfds(jl,jk)-pmfd(jl,jk)*                                   &
              &                (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
            pmfdq(jl,jk)=pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
          ELSE
            pmfd(jl,jk)=0._wp
            pmfds(jl,jk)=0._wp
            pmfdq(jl,jk)=0._wp
            pdmfdp(jl,jk-1)=0._wp
          END IF
        END IF
      END DO
      !
      DO jt=1,ktrac
        DO jl=1,kproma
          IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
            pmfuxt(jl,jk,jt)=pmfuxt(jl,jk,jt)-pmfu(jl,jk)*pxtenh(jl,jk,jt)
            IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
              pmfdxt(jl,jk,jt)=pmfdxt(jl,jk,jt)-pmfd(jl,jk)*pxtenh(jl,jk,jt)
            ELSE
              pmfdxt(jl,jk,jt)=0._wp
            ENDIF
          ELSE
            pmfuxt(jl,jk,jt)=0._wp
            pmfdxt(jl,jk,jt)=0._wp
          ENDIF
        END DO
      END DO
      !
    END DO
    DO jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
      DO jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
          ikb=kcbot(jl)
          zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/(paphp1(jl,klevp1)-paphp1(jl,ikb)))
          zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
          pmfu(jl,jk)=pmfu(jl,ikb)*zzp
          pmfus(jl,jk)=pmfus(jl,ikb)*zzp
          pmfuq(jl,jk)=pmfuq(jl,ikb)*zzp
          pmful(jl,jk)=pmful(jl,ikb)*zzp
        END IF
      END DO
      !
      DO jt=1,ktrac
!DIR$ IVDEP
!OCL NOVREC
        DO jl=1,kproma
          IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
            ikb=kcbot(jl)
            zzp=(paphp1(jl,klevp1)-paphp1(jl,jk))/(paphp1(jl,klevp1)-paphp1(jl,ikb))
            zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
            pmfuxt(jl,jk,jt)=pmfuxt(jl,ikb,jt)*zzp
          ENDIF
        END DO
      END DO
    END DO
    !
    !*    2.  Calculate rain/snow fall rates
    !*        Calculate melting of snow
    !*        Calculate evaporation of precipitation
    !         --------------------------------------
    !
    DO jl=1,kproma
      prfl(jl)=0._wp
      psfl(jl)=0._wp
      prain(jl)=0._wp
    END DO
#ifndef __ICON__
    !++jsr interface for scavenging
    IF (lanysubmodel) THEN
      DO jk=ktopm2,klev
        !---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-----
        zdpg(1:kproma,jk)=(paphp1(1:kproma,jk+1)-paphp1(1:kproma,jk))/grav
        !---End Included for scavenging-----------------------------------------
      END DO
      DO jk=ktopm2,klev
        DO jl=1,kproma
          IF (ldcum(jl)) THEN
            ztc=ptu(jl,jk)-tmelt
            zzfac=MERGE(1._wp,0._wp,ztc<0._wp)
            zalpha=(1._wp-zzfac)+zzfac*(zcaa+(1._wp-zcaa)*exp(-zcab*ztc**2))
            zmlwc(jl,jk)=zalpha*pmwc(jl,jk)
            zmiwc(jl,jk)=(1._wp-zalpha)*pmwc(jl,jk)
            zmratepr(jl,jk)=zalpha*pmrateprecip(jl,jk)
            zmrateps(jl,jk)=(1._wp-zalpha)*pmrateprecip(jl,jk)
          END IF
        END DO
      END DO
    END IF
    !---End Included for scavenging-----------------------------------------
    !--jsr interface for scavenging
#endif
      
    DO jk=ktopm2,klev
      DO jl=1,kproma
        IF(ldcum(jl)) THEN
          prain(jl)=prain(jl)+pdmfup(jl,jk)
          IF(pten(jl,jk).GT.tmelt) THEN
            prfl(jl)=prfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
            IF(psfl(jl).GT.0._wp.AND.pten(jl,jk).GT.ztmelp2) THEN
              zfac=zcons1*(1._wp+vtmpc2*pqen(jl,jk))*(paphp1(jl,jk+1)-paphp1(jl,jk))
              zsnmlt=MIN(psfl(jl),zfac*(pten(jl,jk)-ztmelp2))
              pdpmel(jl,jk)=zsnmlt
              psfl(jl)=psfl(jl)-zsnmlt
              prfl(jl)=prfl(jl)+zsnmlt
            END IF
          ELSE
            psfl(jl)=psfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
          END IF
        END IF
      END DO
#ifndef __ICON__
      IF (lanysubmodel) THEN
      !---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
        zfrain(1:kproma,jk)=prfl(1:kproma)
        zfsnow(1:kproma,jk)=psfl(1:kproma)
      !---End Included for scavenging-----------------------------------------      
       END IF
#endif
    END DO
    DO jl=1,kproma
      prfl(jl)=MAX(prfl(jl),0._wp)
      psfl(jl)=MAX(psfl(jl),0._wp)
      zpsubcl(jl)=prfl(jl)+psfl(jl)
    END DO
#ifndef __ICON__
    !++jsr scavenging interface
    !++mgs ### new code
    !!IF (lanysubmodel) THEN
    !!   IF (lham) THEN
    !!      DO jk=ktopm2,klev
    !!         DO jl=1,kproma
    !!            IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_wp) THEN
    !!               zrhou(jl,jk)=(paphp1(jl,jk+1)+paphp1(jl,jk))* &
    !!                            0.5_wp/(ptu(jl,jk)*rd)
    !!               zcucov=pmfu(jl,jk)/(zwu*zrhou(jl,jk))
    !!               zrfl=zpsubcl(jl)
    !!               zrnew=(MAX(0._wp,SQRT(zrfl/zcucov)-                         &
    !!                      cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
    !!                      MAX(0._wp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
    !!               zrmin=zrfl-zcucov*MAX(0._wp,0.8_wp*pqsen(jl,jk)-pqen(jl,jk))&
    !!                     *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
    !!               zrnew=MAX(zrnew,zrmin)
    !!               zrfln=MAX(zrnew,0._wp)
    !!               zdrfl=MIN(0._wp,zrfln-zrfl)
    !!               pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
    !!!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
    !!               zdpevap=-zdrfl
    !!               zfevapr(jl,jk)=zdpevap*prfl(jl)/zpsubcl(jl)
    !!               zfsubls(jl,jk)=zdpevap*psfl(jl)/zpsubcl(jl)
    !!!---End Included for scavenging-----------------------------------------  
    !!               zpsubcl(jl)=zrfln
    !!            END IF
    !!         END DO
    !!      END DO
    !!   ELSE
    !!      DO jk=ktopm2,klev
    !!         DO jl=1,kproma
    !!            IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_wp) THEN
    !!               zrhou(jl,jk)=(paphp1(jl,jk+1)+paphp1(jl,jk))* &
    !!                            0.5_wp/(ptu(jl,jk)*rd)
    !!               zrfl=zpsubcl(jl)
    !!               zrnew=(MAX(0._wp,SQRT(zrfl/zcucov)-                         &
    !!                      cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
    !!                      MAX(0._wp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
    !!               zrmin=zrfl-zcucov*MAX(0._wp,0.8_wp*pqsen(jl,jk)-pqen(jl,jk))&
    !!                     *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
    !!               zrnew=MAX(zrnew,zrmin)
    !!               zrfln=MAX(zrnew,0._wp)
    !!               zdrfl=MIN(0._wp,zrfln-zrfl)
    !!               pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
    !!!---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
    !!               zdpevap=-zdrfl
    !!               zfevapr(jl,jk)=zdpevap*prfl(jl)/zpsubcl(jl)
    !!               zfsubls(jl,jk)=zdpevap*psfl(jl)/zpsubcl(jl)
    !!!---End Included for scavenging-----------------------------------------  
    !!               zpsubcl(jl)=zrfln
    !!            END IF
    !!         END DO
    !!      END DO
    !!   END IF
    !!ELSE
    !!  DO 240 jk=ktopm2,klev
    !!     DO 235 jl=1,kproma
    !!       IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_wp)  &
    !!                                                                  THEN
    !!           zrfl=zpsubcl(jl)
    !!           zrnew=(MAX(0._wp,SQRT(zrfl/zcucov)-                         &
    !!                        cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
    !!                        MAX(0._wp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
    !!           zrmin=zrfl-zcucov*MAX(0._wp,0.8_wp*pqsen(jl,jk)-pqen(jl,jk))&
    !!                        *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
    !!           zrnew=MAX(zrnew,zrmin)
    !!           zrfln=MAX(zrnew,0._wp)
    !!           zdrfl=MIN(0._wp,zrfln-zrfl)
    !!           pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
    !!           zpsubcl(jl)=zrfln
    !!       END IF
    !!235  END DO
    !!240 END DO
    !!END IF
    IF (lanysubmodel) THEN
      DO jk=ktopm2,klev
        zrhou(1:kproma,jk)=(paphp1(1:kproma,jk+1)+paphp1(1:kproma,jk))*              &
          &               0.5_wp/(ptu(1:kproma,jk)*rd)
      END DO
    END IF
#endif
    DO jk=ktopm2,klev
      zdpevap(1:kproma) = 0._wp
#ifndef __ICON__
      IF (lanysubmodel) zpsubcl_sav(1:kproma) = zpsubcl(1:kproma)
      ! ### value of zcucov depends on lham => explicit submodel dependence !
      IF (lham) THEN
        zcucov(1:kproma) = pmfu(1:kproma,jk)/(zwu*zrhou(1:kproma,jk))
      ELSE
        zcucov(1:kproma) = 0.05_wp
      END IF
#else
      zcucov(1:kproma) = 0.05_wp
#endif

      DO jl=1,kproma
        IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_wp) THEN
          zrfl=zpsubcl(jl)
          zrnew=(MAX(0._wp,SQRT(zrfl/zcucov(jl))-                                    &
            &    cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*                        &
            &    MAX(0._wp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov(jl)
          zrmin=zrfl-zcucov(jl)*MAX(0._wp,0.8_wp*pqsen(jl,jk)-pqen(jl,jk))           &
            &   *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
          zrnew=MAX(zrnew,zrmin)
          zrfln=MAX(zrnew,0._wp)
          zdrfl=MIN(0._wp,zrfln-zrfl)
          zdpevap(jl) = -zdrfl
          pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
          zpsubcl(jl)=zrfln
        END IF
      END DO

#ifndef __ICON__
      !---Included for scavenging in wetdep_interface (Philip Stier, 28/03/01):-------
      IF (lanysubmodel) THEN
        DO jl=1,kproma
          IF (zpsubcl_sav(jl) > 1.e-20_wp) THEN
            zfevapr(jl,jk)=zdpevap(jl)*prfl(jl)/zpsubcl_sav(jl)
            zfsubls(jl,jk)=zdpevap(jl)*psfl(jl)/zpsubcl_sav(jl)
          ELSE
            zfevapr(jl,jk) = 0._wp
            zfsubls(jl,jk) = 0._wp
          END IF
        END DO
      END IF
      !---End Included for scavenging-----------------------------------------
#endif
    END DO
    !--mgs ### new code
    !
#ifndef __ICON__
    IF (lanysubmodel) THEN
      CALL cuflx_subm(kbdim,  kproma,    klev,     ktopm2,   & ! dimensions
        &             krow,                                  & ! longitude (kproma-block) index
#ifdef HAMMOZ /* SF */
                      ktype, kcbot,    kctop,                & ! conv. cloud extent
#endif
        &             pxtenh, pxtu, zrhou,                   & ! tracers
        &             pmfu,   pmfuxt,                        & ! convective fluxes and corresp. mmr
        &             zmlwc,  zmiwc,     zmratepr, zmrateps, & ! cloud properties
        &             zfrain, zfsnow,    zfevapr,  zfsubls,  & !   "       "
        &             zaclc,  zmsnowacl,                     & !   "       "
        &             ptu,    zdpg,                          & ! thermodynamic quantities
#ifdef HAMMOZ /* SF */
!SF MGS for lightning
                      pten,   ptenh, paphp1,                 &             
#endif
        &             pxtte                                  &
#ifdef HAMMOZ /* SF */
                     ,pxtbound                               &
#endif
                   ) 
    END IF
#endif

    !!baustelle!! (?)
    DO jl=1,kproma
      zrsum=prfl(jl)+psfl(jl)
      zdpevap(jl)=zpsubcl(jl)-zrsum
      prfl(jl)=prfl(jl)+zdpevap(jl)*prfl(jl)*(1._wp/MAX(1.e-20_wp,zrsum))
      psfl(jl)=psfl(jl)+zdpevap(jl)*psfl(jl)*(1._wp/MAX(1.e-20_wp,zrsum))
    END DO
    !
#ifndef __ICON__
    IF ( locfdiag ) THEN
      DO jk=ktopm2, klev
        IF (lham) THEN
          zcucov(1:kproma) = pmfu(1:kproma,jk)/(zwu*zrhou(1:kproma,jk))
        ELSE
          zcucov(1:kproma) = 0.05_wp
        END IF
        afac(1:kproma)=zcucov(1:kproma)/gboxarea(1:kproma)
        imc(1:kproma,jk,krow)=afac(1:kproma)*(pmfu(1:kproma,jk)+pmfd(1:kproma,jk))
        DO jl=1,kproma
          IF ( ktype(jl).EQ.1 ) THEN
            dmc (jl,jk,krow)  = dmc(jl,jk,krow) + imc(jl,jk,krow) * delta_time
          END IF
          IF ( ktype(jl).EQ.2 ) THEN
            smc (jl,jk,krow)  = smc(jl,jk,krow) + imc(jl,jk,krow) * delta_time
          END IF
        END DO        
        mc (1:kproma,jk,krow)=mc(1:kproma,jk,krow)+imc(1:kproma,jk,krow)*delta_time
        mcu (1:kproma,jk,krow) = mcu(1:kproma,jk,krow) + afac(1:kproma) *            &
          &                     (pmfu(1:kproma,jk)) * delta_time
        mcd (1:kproma,jk,krow) = mcd(1:kproma,jk,krow) + afac(1:kproma) *            &
          &                     (pmfd(1:kproma,jk)) * delta_time
      END DO
    END IF ! locfdiag
    !
    IF ( locospoffl ) THEN
      DO jk=ktopm2,klev
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            IF(pten(jl,jk).GT.tmelt) THEN
              cospoffl_ccrain(jl,jk,krow) = pdmfup(jl,jk)+pdmfdp(jl,jk)
              cospoffl_ccsnow(jl,jk,krow) = 0._wp
            ELSE
              cospoffl_ccsnow(jl,jk,krow) = pdmfup(jl,jk)+pdmfdp(jl,jk)
              cospoffl_ccrain(jl,jk,krow) = 0._wp
            END IF
          END IF
        END DO
      END DO
    ENDIF !locospoffl
#endif
  END SUBROUTINE cuflx
  !>
  !!
  SUBROUTINE cudtdq(kproma, kbdim, klev, klevp1, ktopm2, ldcum, ktrac,               &
#ifndef __ICON__
    &               krow,                                                            &
    &               ptte,     pqte,                                                  &
    &               pxtte,                                                           &
    &               paprc,    paprs,                                                 &
#endif
    &               paphp1,   pten,                                                  &
    &               pmfuxt,   pmfdxt,                                                &
    &               pmfus,    pmfds,    pmfuq,    pmfdq,                             &
    &               pmful,    pdmfup,   pdmfdp,   plude,                             &
    &               pdpmel,   prfl,     psfl,                                        &
    &               pcpen,    palvsh,   pqtec,    pqude,                             &
    &               prsfc,    pssfc,                                                 &
    &               pch_con,  pcw_con,                                               &
#ifdef __ICON__
    &               pcon_dtrl, pcon_dtri, pcon_iqte,                                 &
    &               ptte_cnv, pqte_cnv, pxtte_cnv,                                   &
#endif
    &               pxtecl,   pxteci                                                )
    !

!>>SF
#ifdef HAMMOZ
  USE mo_boundary_condition, ONLY: bc_set, bc_find
  USE mo_time_control,       ONLY: lstart, lresume
  USE mo_param_switches,     ONLY: lcdnc_progn
#endif
!<<SF

    INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2, ktrac
#ifndef __ICON__
    INTEGER, INTENT (IN) :: krow
#else
    REAL(wp),INTENT(INOUT) :: pcon_dtrl(kbdim), pcon_dtri(kbdim)
    REAL(wp),INTENT(INOUT) :: pcon_iqte(kbdim) ! integrated qv tendency
    REAL(wp),INTENT(OUT) :: ptte_cnv(kbdim,klev)
    REAL(wp),INTENT(OUT) :: pqte_cnv(kbdim,klev), pxtte_cnv(kbdim,klev,ktrac)
#endif
    LOGICAL  llo1
    !
    REAL(wp) :: pten(kbdim,klev),        paphp1(kbdim,klevp1),                       &
#ifndef __ICON__
      &         ptte(kbdim,klev),        pqte(kbdim,klev),                           &
      &         pxtte(kbdim,klev,ktrac),                                             &
      &         paprc(kbdim),            paprs(kbdim),                               &
#endif
      &         prsfc(kbdim),            pssfc(kbdim)
    REAL(wp) :: pxtecl(kbdim,klev),      pxteci(kbdim,klev)
    REAL(wp) :: pmfus(kbdim,klev),       pmfds(kbdim,klev),                          &
      &         pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                          &
      &         pmful(kbdim,klev),       plude(kbdim,klev),                          &
      &         pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),                         &
      &         pqtec(kbdim,klev),       pqude(kbdim,klev),                          &
      &         prfl(kbdim)
    REAL(wp) :: pdpmel(kbdim,klev),      psfl(kbdim)
    REAL(wp) :: pcpen(kbdim,klev),       palvsh(kbdim,klev),                         &
      &         pch_con(kbdim),          pcw_con(kbdim)
    LOGICAL  :: ldcum(kbdim)
    REAL(wp) :: zmelt(kbdim), zcpten(kbdim,klev), zqten(kbdim,klev)
    REAL(wp) :: zsheat(kbdim)
    REAL(wp) :: pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)
    REAL(wp) :: zxtec(kbdim,klev)
    REAL(wp) :: zrcpm ! reciprocal value of specific heat of moist air
    INTEGER  :: jl, jk, jt
    REAL(wp) :: zalv, zdtdt, zdqdt, zdxtdt
#ifndef __ICON__
    REAL(wp) :: zdiagt
    !
    !----------------------------------------------------------------------
    !
    !*    1.0          Specify parameters
    !                  ------------------
    !
    zdiagt=delta_time
    !
#else
    ptte_cnv(:,:)       = 0._wp
    pqte_cnv(:,:)       = 0._wp
    pxtte_cnv(:,:,:)    = 0._wp
    pcon_dtrl(1:kproma) = 0._wp
    pcon_dtri(1:kproma) = 0._wp
    pcon_iqte(1:kproma) = 0._wp
#endif

!>>SF Missing initialization of zxtec
    zxtec(1:kproma,:) = 0._wp
!<<SF
    !----------------------------------------------------------------------
    !
    !*    2.0          Incrementation of t and q tendencies
    !                  ------------------------------------
    !
    DO jl=1,kproma
      zmelt(jl)=0._wp
      zsheat(jl)=0._wp
    END DO
    !
    zcpten(1:kproma,ktopm2:klev)=0._wp
    zqten (1:kproma,ktopm2:klev)=0._wp
    pch_con(1:kproma)=0._wp
    pcw_con(1:kproma)=0._wp
    !
    DO jk=ktopm2,klev
      !
      IF(jk.LT.klev) THEN
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            llo1=(pten(jl,jk)-tmelt).GT.0._wp
            zalv=MERGE(alv,als,llo1)
            zrcpm=1._wp/pcpen(jl,jk)
            zdtdt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*                      &
              &                 (pmfus(jl,jk+1)-pmfus(jl,jk)+                        &
              &                  pmfds(jl,jk+1)-pmfds(jl,jk)-                        &
              &                  alf*pdpmel(jl,jk)-                                  &
              &                  palvsh(jl,jk+1)*pmful(jl,jk+1)+                     &
              &                  palvsh(jl,jk)*pmful(jl,jk)+                         &
              &                  zalv*(plude(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
            zcpten(jl,jk)=zdtdt*pcpen(jl,jk)
            zdqdt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                            &
              &                 (pmfuq(jl,jk+1)-pmfuq(jl,jk)+                        &
              &                  pmfdq(jl,jk+1)-pmfdq(jl,jk)+                        &
              &                  pmful(jl,jk+1)-pmful(jl,jk)-                        &
              &                  plude(jl,jk)-                                       &
              &                  (pdmfup(jl,jk)+pdmfdp(jl,jk)))
            zxtec(jl,jk)=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*plude(jl,jk)
            zqten(jl,jk) = zdqdt + zxtec(jl,jk)
            pxteci(jl,jk)=MERGE(0.0_wp,zxtec(jl,jk),llo1)
            pxtecl(jl,jk)=MERGE(zxtec(jl,jk),0.0_wp,llo1)
            pqtec(jl,jk)=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*pqude(jl,jk)
            zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
            zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
#ifndef __ICON__
            ptte(jl,jk)=ptte(jl,jk)+zdtdt
            pqte(jl,jk)=pqte(jl,jk)+zdqdt
#else
            ptte_cnv(jl,jk)=zdtdt
            pqte_cnv(jl,jk)=zdqdt
#endif
          ENDIF
        END DO
        !
#ifndef __ICON__
        IF (trlist% anyconv /= 0) THEN
#endif
          DO jt=1,ktrac
#ifndef __ICON__
            IF (trlist% ti(jt)% nconv == 1) THEN
#endif
              DO jl=1,kproma
                IF(ldcum(jl)) THEN
                  zdxtdt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))                      &
                    &         *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt)                  &
                    &          +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))
#ifndef __ICON__
                  pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
#else
                  pxtte_cnv(jl,jk,jt)=zdxtdt
#endif
                ENDIF
              END DO
#ifndef __ICON__
            ENDIF
#endif
          END DO
#ifndef __ICON__
        ENDIF
#endif
      !
      ELSE
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            llo1=(pten(jl,jk)-tmelt).GT.0._wp
            zalv=MERGE(alv,als,llo1)
            zrcpm=1._wp/pcpen(jl,jk)
            zdtdt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*                     &
              &    (pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk)-                     &
              &     palvsh(jl,jk)*pmful(jl,jk)-                                      &
              &     zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk)+plude(jl,jk)))
            zcpten(jl,jk)=zdtdt*pcpen(jl,jk)
            zdqdt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                           &
              &       (pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)+                       &
              &       (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
            zxtec(jl,jk)=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*plude(jl,jk)
            zqten(jl,jk) = zdqdt + zxtec(jl,jk)
            pxteci(jl,jk)=MERGE(0.0_wp,zxtec(jl,jk),llo1)
            pxtecl(jl,jk)=MERGE(zxtec(jl,jk),0.0_wp,llo1)
            pqtec(jl,jk)=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*pqude(jl,jk)
            zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
            zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
#ifndef __ICON__
            ptte(jl,jk)=ptte(jl,jk)+zdtdt
            pqte(jl,jk)=pqte(jl,jk)+zdqdt
#else
            ptte_cnv(jl,jk)=zdtdt
            pqte_cnv(jl,jk)=zdqdt
#endif
          END IF
        END DO
        !
#ifndef __ICON__
        IF (trlist% anyconv /= 0) THEN
#endif
          DO jt=1,ktrac
#ifndef __ICON__
            IF (trlist% ti(jt)% nconv == 1) THEN
#endif
              DO jl=1,kproma
                IF(ldcum(jl)) THEN
                  zdxtdt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))                     &
                    &    *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
#ifndef __ICON__
                  pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
#else
                  pxtte_cnv(jl,jk,jt)=zdxtdt
#endif
                ENDIF
              END DO
#ifndef __ICON__
            END IF
#endif
          END DO
#ifndef __ICON__
        ENDIF
#endif
      !
      END IF
      !
    END DO

!>>SF #518: save detrained condensate as a new boundary condition to use in 2-moment stratiform cl-micro
#ifdef HAMMOZ
    IF (lcdnc_progn) THEN
        IF (lstart .OR. lresume) THEN
           CALL bc_find('Detrained condensate', ibc_detr_cond, ierr=ierr)
        ENDIF
        CALL bc_set(ibc_detr_cond, kbdim, krow, zxtec)
    ENDIF
#endif
!<<SF

    !
    !---------------------------------------------------------------------
    !
    !      3.          Update surface fields
    !                  ---------------------
    !
    DO jl=1,kproma
      prsfc(jl)=prfl(jl)
      pssfc(jl)=psfl(jl)
#ifndef __ICON__
    ! set not accumulated variables first
      aprc_na(jl,krow)=prfl(jl)+psfl(jl)
    ! set accumulated variables
      paprc(jl)=paprc(jl)+zdiagt*aprc_na(jl,krow)
      paprs(jl)=paprs(jl)+zdiagt*psfl(jl)
#endif
    END DO
    !
    ! column integral of convective heating and moistening
    DO jk=ktopm2,klev
      DO jl=1,kproma
        pch_con(jl)=pch_con(jl)+zcpten(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
        pcw_con(jl)=pcw_con(jl)+zqten(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
      END DO
    END DO
    DO jl=1,kproma
      pch_con(jl)=pch_con(jl)-(alv*prsfc(jl)+als*pssfc(jl))
      pcw_con(jl)=pcw_con(jl)+prsfc(jl)+pssfc(jl)
    END DO
    !
#ifdef __ICON__
    ! do we need to account for the surface, or for the top 2 layers?
    DO jk=ktopm2,klev
      DO jl=1,kproma
        ! water vapor
        pcon_iqte(jl)=pcon_iqte(jl)+                                                 &
          &   pqte_cnv(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
        ! detrained liquid water
        pcon_dtrl(jl)=pcon_dtrl(jl)+                                                 &
          &   pxtecl(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
        ! detrained ice
        pcon_dtri(jl)=pcon_dtri(jl)+                                                 &
          &   pxteci(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))/grav
      END DO
    END DO
#else
    ! calculate and store convective accumulated precipitation (mm)
    IF (lanysubmodel) THEN
      CALL set_vphysc_var (kproma, klev, krow, prflconv=prfl, psflconv=psfl)
    ENDIF
#endif
    !
  END SUBROUTINE cudtdq
  !>
  !!
  SUBROUTINE cududv(   kproma,   kbdim,    klev,     klevp1,                         &
    &        ktopm2,   ktype,    kcbot,    paphp1,   ldcum,                          &
    &        puen,     pven,                                                         &
#ifndef __ICON__
    &        pvom,     pvol,                                                         &
#else
    &        pvom_cnv, pvol_cnv,                                                     &
#endif
    &        puu,      pud,      pvu,      pvd,                                      &
    &        pmfu,     pmfd)
    !
    INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2
    REAL(wp):: puen(kbdim,klev),        pven(kbdim,klev),                            &
#ifndef __ICON__
      &        pvol(kbdim,klev),        pvom(kbdim,klev),                            &
#else
      &        pvom_cnv(kbdim,klev),    pvol_cnv(kbdim,klev),                        &
#endif
      &        paphp1(kbdim,klevp1)
    REAL(wp):: puu(kbdim,klev),         pud(kbdim,klev),                             &
      &        pvu(kbdim,klev),         pvd(kbdim,klev),                             &
      &        pmfu(kbdim,klev),        pmfd(kbdim,klev)
    INTEGER :: ktype(kbdim),            kcbot(kbdim)
    LOGICAL :: ldcum(kbdim)
    !
    REAL(wp):: zmfuu(kbdim,klev),       zmfdu(kbdim,klev),                           &
      &        zmfuv(kbdim,klev),       zmfdv(kbdim,klev)
    INTEGER :: jl, jk, ik, ikb
    REAL(wp):: zzp, zdudt, zdvdt

#ifdef __ICON__
    pvom_cnv(:,:) = 0._wp
    pvol_cnv(:,:) = 0._wp
#endif
    !
    !----------------------------------------------------------------------
    !
    !*    1.0          Calculate fluxes and update u and v tendencies
    !                  ----------------------------------------------
    !
    IF(ktopm2.EQ.1) THEN
      DO jk=2,klev
        ik=jk-1
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
            zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
            zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
            zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
        END DO
      END DO
      DO jl=1,kproma
        IF(ldcum(jl)) THEN
          zmfuu(jl,1)=zmfuu(jl,2)
          zmfuv(jl,1)=zmfuv(jl,2)
          zmfdu(jl,1)=zmfdu(jl,2)
          zmfdv(jl,1)=zmfdv(jl,2)
        END IF
      END DO
    ELSE
      DO jk=ktopm2,klev
        ik=jk-1
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
            zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
            zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
            zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
        END DO
      END DO
    END IF
    DO jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
      DO jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
          ikb=kcbot(jl)
          zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/(paphp1(jl,klevp1)-paphp1(jl,ikb)))
          zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
          zmfuu(jl,jk)=zmfuu(jl,ikb)*zzp
          zmfuv(jl,jk)=zmfuv(jl,ikb)*zzp
          zmfdu(jl,jk)=zmfdu(jl,ikb)*zzp
          zmfdv(jl,jk)=zmfdv(jl,ikb)*zzp
        END IF
      END DO
    END DO
    !
    DO jk=ktopm2,klev
      !
      IF(jk.LT.klev) THEN
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zdudt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                            &
              &         (zmfuu(jl,jk+1)-zmfuu(jl,jk)+zmfdu(jl,jk+1)-zmfdu(jl,jk))
            zdvdt=(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                            &
              &         (zmfuv(jl,jk+1)-zmfuv(jl,jk)+zmfdv(jl,jk+1)-zmfdv(jl,jk))
#ifndef __ICON__
            pvom(jl,jk)=pvom(jl,jk)+zdudt
            pvol(jl,jk)=pvol(jl,jk)+zdvdt
#else
            pvom_cnv(jl,jk)=zdudt
            pvol_cnv(jl,jk)=zdvdt
#endif
          END IF
        END DO
        !
      ELSE
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zdudt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                           &
              &            (zmfuu(jl,jk)+zmfdu(jl,jk))
            zdvdt=-(grav/(paphp1(jl,jk+1)-paphp1(jl,jk)))*                           &
              &            (zmfuv(jl,jk)+zmfdv(jl,jk))
#ifndef __ICON__
            pvom(jl,jk)=pvom(jl,jk)+zdudt
            pvol(jl,jk)=pvol(jl,jk)+zdvdt
#else
            pvom_cnv(jl,jk)=zdudt
            pvol_cnv(jl,jk)=zdvdt
#endif
          END IF
        END DO
      END IF
      !
    END DO
  END SUBROUTINE cududv

END MODULE mo_cufluxdts
