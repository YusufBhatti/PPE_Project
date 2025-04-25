#ifdef __xlC__
@PROCESS STRICT
#endif
#ifndef __xlC__
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#define SWDIV_NOCHK(a,b) ((a)/(b))
#endif

!>
!! @brief Master routine - provides interface for cumulus parameterization
!!
!! @remarks
!!     *cucall*
!!     This routine is called from *physc* (ECHAM) or *mo_echam_phy_main* (ICON)
!!     It provides input for *cumastr*, receives updated tendencies, precipitation
!!     *cumastr*
!!     This routine computes the physical tendencies of the prognostic variables
!!     t,q,u and v due to convective processes. Processes considered are: convective
!!     fluxes, formation of precipitation, evaporation of falling rain below cloud
!!     base, saturated cumulus downdrafts.
!!     *cumastr* is called from *cucall*
!!     It takes its input from the long-term storage t,q,u,v,phi and p and moisture
!!     tendencies. It returns its output to the same space
!!     1. modified tendencies of model variables
!!     2. rates of convective precipitation (for surface)
!!     Method: Parameterization is done using a massflux-scheme.
!!        (1) Define constants and parameters
!!        (2) Specify values (t,q,qs...) at half levels and
!!            initialize updraft- and downdraft-values in 'cuini'
!!        (3) Calculate cloud base in 'cubase'
!!            and specify cloud base massflux from pbl moisture budget
!!        (4) Do cloud ascent in 'cuasc' in absence of downdrafts
!!        (5) Do downdraft calculations:
!!              (A) Determine values at lfs in 'cudlfs'
!!              (B) Determine moist descent in 'cuddraf'
!!              (C) Recalculate cloud base massflux considering the 
!!                  effect of cumulus-downdrafts
!!        (6) Do final cloud ascent in 'cuasc'
!!        (7) Do final adjustments to convective fluxes in 'cuflx',
!!            do evaporation in subcloud layer
!!        (8) Calculate increments of t and q in 'cudtdq'
!!        (9) Calculate increments of u and v in 'cududv'
!!
!! @author M. Tiedtke, ECMWF,    1986/1987/1989
!!         M. Tiedtke, ECMWF,    Dec 1989
!!
!! @references.
!!       Tiedtke, 1989: Mon. Wea. Rev., 117, 1779-1800
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
MODULE mo_cumastr
  USE mo_kind,                 ONLY: wp
  USE mo_physical_constants,   ONLY: grav, alv, als, tmelt, vtmpc1, rd, cpd, cpv, vtmpc2
  USE mo_echam_convect_tables, ONLY: prepare_ua_index_spline,lookup_ua_spline, lookup_ubc
  USE mo_echam_conv_constants, ONLY: entrpen, entrscv, lmfdd, cmfdeps, lmfdudv, cmftau
  USE mo_cuinitialize,         ONLY: cuini, cubase
  USE mo_cuascent,             ONLY: cuasc
  USE mo_cudescent,            ONLY: cudlfs, cuddraf
  USE mo_cufluxdts,            ONLY: cuflx, cudtdq, cududv
#ifndef __ICON__
  USE mo_param_switches,       ONLY: lconvmassfix
  USE mo_tracer_processes,     ONLY: xt_conv_massfix
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cucall

CONTAINS 
  !>
  !!
  SUBROUTINE cucall(   kproma, kbdim, klev, klevp1, klevm1,            &! in
    &                  ktrac,                                          &! in
    &                  ptime_step_len,                                 &! in
    &                  ldland,                                         &! in
    &                  ptm1,     pum1,     pvm1,                       &! in
    &                  pqm1,     pxlm1,    pxim1,                      &! in
    &                  pxtm1,                                          &! in
    &                  pqte,     pxlte,    pxite,                      &! in
    &                  pverv,    pqhfla,   pgeo,                       &! in
    &                  papp1,    paphp1,   pthvsig,                    &! in
    &                  ptte,     pvom,     pvol,                       &! in
    &                  pxtte,                                          &! inout
    &                  pqtec,                                          &! inout
    &                  pch_con,  pcw_con,                              &! inout
    &                  prsfc,    pssfc,                                &! out
    &                  pxtecl,   pxteci,                               &! out (?)
    &                  ktype,    ictop,    ilab,                       &! out
#ifndef __ICON__
    &                  krow,                                           &! in
    &                  paprc,    paprs,                                &! out (?)
    &                  ptopmax                                        ) ! inout
#else
    &                  ptopmax,                                        &! inout
    &                  nmctop, cevapcu,                                &! in
    &                  pcd, pcv,                                       &! in
    &                  pqte_dyn, pqte_phy,                             &! in
    &                  pcon_dtrl, pcon_dtri, pcon_iqte,                &! inout
    &                  ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv,         &! out
    &                  pxtte_cnv                               )        ! out
#endif

    INTEGER, INTENT (IN) :: klev, klevm1, klevp1, kproma, kbdim, ktrac
    REAL(wp),INTENT (IN) :: ptime_step_len
#ifdef __ICON__
    REAL(wp),INTENT (IN) :: cevapcu(:)
    INTEGER, INTENT (IN) :: nmctop
    REAL(wp),INTENT (IN) :: pcd, pcv
#endif
    REAL(wp)::  ptm1(kbdim,klev),         pqm1(kbdim,klev),              &
      &         pum1(kbdim,klev),         pvm1(kbdim,klev),              &
      &         pqte(kbdim,klev),                                        &
      &         pverv(kbdim,klev),        pgeo(kbdim,klev),              &
      &         papp1(kbdim,klev),        paphp1(kbdim,klevp1),          &
      &         prsfc(kbdim),             pssfc(kbdim)
    REAL(wp)::  ptte(kbdim,klev),                                        &
      &         pvom(kbdim,klev),         pvol(kbdim,klev)
    REAL(wp)::  pxtte(kbdim,klev,ktrac)
#ifndef __ICON__
    INTEGER, INTENT (IN) :: krow
    REAL(wp)::  paprc(kbdim),             paprs(kbdim)
#else
    REAL(wp),INTENT(INOUT) :: pcon_dtrl(kbdim), pcon_dtri(kbdim)
    REAL(wp),INTENT(INOUT) :: pcon_iqte(kbdim)
    REAL(wp),INTENT(INOUT) :: ptte_cnv(kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pqte_cnv(kbdim,klev), pxtte_cnv(kbdim,klev,ktrac)
    REAL(wp),INTENT(IN)    :: pqte_dyn(kbdim,klev), pqte_phy(kbdim,klev)
#endif
    REAL(wp)::  pthvsig(kbdim)
    INTEGER ::  ktype(kbdim)
    REAL(wp)::  pqhfla(kbdim)
    REAL(wp)::  pch_con(kbdim),           pcw_con(kbdim)
    REAL(wp)::  ptopmax(kbdim)
    INTEGER ::  ilab(kbdim,klev)
    REAL(wp)::  pqtec(kbdim,klev)
    REAL(wp)::  pxtecl(kbdim,klev),       pxteci(kbdim,klev)
    REAL(wp)::  pxlm1(kbdim,klev),        pxim1(kbdim,klev),             &
      &         pxlte(kbdim,klev),        pxite(kbdim,klev)
    REAL(wp)::  ztp1(kbdim,klev),         zqp1(kbdim,klev),              &
      &         zxp1(kbdim,klev),         ztvp1(kbdim,klev),             &
      &         zup1(kbdim,klev),         zvp1(kbdim,klev),              &
      &         ztu(kbdim,klev),          zqu(kbdim,klev),               &
      &         zlu(kbdim,klev),          zlude(kbdim,klev),             &
      &         zqude(kbdim,klev),                                       &
      &         zcpq(kbdim,klev),         zcq(kbdim,klev),               &
      &         zmfu(kbdim,klev),         zmfd(kbdim,klev),              &
      &         zqsat(kbdim,klev),        zrain(kbdim),                  &
      &         za(kbdim),                ua(kbdim)
    INTEGER ::  itopec2(kbdim),           idx(kbdim)
    INTEGER ::  icbot(kbdim),             ictop(kbdim)
    REAL(wp)::  zxtp1(kbdim,klev,ktrac),  zxtu(kbdim,klev,ktrac),        &
      &         pxtm1(kbdim,klev,ktrac)
    REAL(wp)::  ztopmax(kbdim)
    LOGICAL ::  locum(kbdim),             ldland(kbdim)
    !
    !  Local scalars: 
    REAL(wp):: ztmst, zxlp1, zxip1
    INTEGER :: ilevmin, jk, jl, jt
#ifdef __ICON__
    REAL(wp)::  zqte_dyn_phy(kbdim,klev)
#endif
    !
    !  Executable statements 
    !
    !-----------------------------------------------------------------------
    !*    1.           Calculate t,q and qs at main levels
    !*                 -----------------------------------
    !
    ztmst=ptime_step_len
    DO jk=1,klev
!IBM* NOVECTOR
      DO jl=1,kproma
        ztp1(jl,jk)=ptm1(jl,jk)+ptte(jl,jk)*ztmst
        zqp1(jl,jk)=MAX(0._wp,pqm1(jl,jk)+pqte(jl,jk)*ztmst)
        zxlp1=pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
        zxip1=pxim1(jl,jk)+pxite(jl,jk)*ztmst
        zxp1(jl,jk)=MAX(0._wp,zxlp1+zxip1)
        ztvp1(jl,jk)=ztp1(jl,jk)*(1._wp+vtmpc1*zqp1(jl,jk)-zxp1(jl,jk))
        zup1(jl,jk)=pum1(jl,jk)+pvom(jl,jk)*ztmst
        zvp1(jl,jk)=pvm1(jl,jk)+pvol(jl,jk)*ztmst
      END DO

      CALL prepare_ua_index_spline('cucall',kproma,ztp1(1,jk),idx(1),za(1))
      CALL lookup_ua_spline(kproma,idx(1),za(1),ua(1))

!IBM* NOVECTOR
      DO jl=1,kproma
        zqsat(jl,jk)=ua(jl)/papp1(jl,jk)
        zqsat(jl,jk)=MIN(0.5_wp,zqsat(jl,jk))
        zqsat(jl,jk)=zqsat(jl,jk)/(1._wp-vtmpc1*zqsat(jl,jk))
#ifndef __ICON__
        zcpq(jl,jk)=cpd*(1.0_wp+vtmpc2*MAX(pqm1(jl,jk),0.0_wp))
#else
        zcpq(jl,jk)=cpd+(cpv-cpd)*MAX(pqm1(jl,jk),0.0_wp) ! cp of moist air for comp. of fluxes
        zcq (jl,jk)=pcd+(pcv-pcd)*MAX(pqm1(jl,jk),0.0_wp) ! cp or cv of moist air for temp tendency
#endif
      END DO

      DO jt=1,ktrac
        DO jl=1,kproma
          zxtp1(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst
        END DO
      END DO
    END DO
    DO jl=1,kproma
      zrain(jl)=0._wp
      locum(jl)=.FALSE.
    END DO
    !
    !-----------------------------------------------------------------------
    !
    !*    2.     Call 'cumastr'(master-routine for cumulus parameterization)
    !*           -----------------------------------------------------------
    !
#ifdef __ICON__
    ! Total moisture tendency due to transport and the parameterized processes
    ! computed before "cucall".
    zqte_dyn_phy(1:kproma,:) = pqte_dyn(1:kproma,:)+pqte_phy(1:kproma,:)
#endif

    CALL cumastr(kproma, kbdim, klev, klevp1, klevm1, ilab,           &
      &          ptime_step_len,                                      &
      &          ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
      &          ztvp1,    ktrac,    ldland,                          &
      &          zxtp1,    zxtu,                                      &
      &          pverv,    zqsat,    pqhfla,                          &
      &          paphp1,   pgeo,                                      &
#ifndef __ICON__
      &          pqte,                                                &
#else
      &          zqte_dyn_phy,                                        &
#endif
      &          prsfc,    pssfc,                                     &
      &          pqtec,    zqude,    zcpq,                            &
      &          pch_con,  pcw_con,                                   &
      &          locum,    ktype,    icbot,    ictop,                 &
      &          ztu,      zqu,      zlu,      zlude,                 &
      &          zmfu,     zmfd,     zrain,    pthvsig,               &
#ifndef __ICON__
      &          krow,     papp1,                                     &
      &          pxtte,                                               &
      &          ptte,     pvom,     pvol,                            &
      &          paprc,    paprs,                                     &
#else
      &          nmctop,   cevapcu,                                   &
      &          zcq,                                                 &
      &          pcon_dtrl, pcon_dtri, pcon_iqte,                     &
      &          ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv,pxtte_cnv,    &
#endif
      &          pxtecl,   pxteci                                     )
    !
    ! ------------------------------------------------------------------
    !
    !*     3.     Pressure altitude of convective cloud tops.
    !             -------- -------- -- ---------- ----- -----
    !
    ilevmin=klev-4
    !
    DO jl=1,kproma
      itopec2(jl)=klevp1
    END DO
    !
    DO jk=1,ilevmin
      DO jl=1,kproma
        IF(ilab(jl,jk).EQ.2 .AND. itopec2(jl).EQ.klevp1) THEN
          itopec2(jl)=jk
        END IF
      END DO
    END DO
    !
    ztopmax(1:kproma) = ptopmax(1:kproma)

    DO jl=1,kproma
      IF(itopec2(jl).EQ.1) THEN
        ptopmax(jl)=papp1(jl,1)
      ELSE IF(itopec2(jl).NE.klevp1) THEN
        ptopmax(jl)=paphp1(jl,itopec2(jl))
      ELSE
        ptopmax(jl)=99999._wp
      END IF
      ptopmax(jl)=MIN(ptopmax(jl),ztopmax(jl))
    END DO
    !
    !---------------------------------------------------------------------
    !
  END SUBROUTINE cucall
  !>
  !!
  SUBROUTINE cumastr(  kproma, kbdim, klev, klevp1, klevm1, ilab,         &
    &                  ptime_step_len,                                    &!in
    &        pten,     pqen,     pxen,     puen,     pven,                &
    &        ptven,    ktrac,    ldland,                                  &
    &        pxten,    pxtu,                                              &
    &        pverv,    pqsen,    pqhfla,                                  &
    &        paphp1,   pgeo,                                              &
    &        pqte,                                                        &
    &        prsfc,    pssfc,                                             &
    &        pqtec,    pqude,    pcpen,                                   &
    &        pch_con,  pcw_con,                                           &
    &        ldcum,    ktype,    kcbot,    kctop,                         &
    &        ptu,      pqu,      plu,      plude,                         &
    &        pmfu,     pmfd,     prain,    pthvsig,                       &
#ifndef __ICON__
    &        krow,    papp1,                                              &
    &        pxtte,                                                       &
    &        ptte,     pvom,     pvol,                                    &
    &        paprc,    paprs,                                             &
#else
    &        nmctop,   cevapcu,                                           &
    &        pcen,                                                        &
    &        pcon_dtrl, pcon_dtri, pcon_iqte,                             &
    &        ptte_cnv, pvom_cnv, pvol_cnv, pqte_cnv, pxtte_cnv,           &
#endif
    &        pxtecl,   pxteci                                             )
    !
    INTEGER, INTENT(IN)   :: kproma, kbdim, klev, klevp1, ktrac, klevm1
    REAL(wp),INTENT(IN)   :: ptime_step_len
    REAL(wp),INTENT(INOUT):: pch_con(kbdim), pcw_con(kbdim)
#ifndef __ICON__
    INTEGER, INTENT(IN)   :: krow
    REAL(wp),INTENT(IN)   :: papp1(kbdim,klev)
    REAL(wp),INTENT(IN)   :: ptte(kbdim,klev),                                       &
      &                      pvom(kbdim,klev),        pvol(kbdim,klev),              &
      &                      paprc(kbdim),            paprs(kbdim)
    REAL(wp),INTENT(INOUT):: pxtte(kbdim,klev,ktrac)
#else
    INTEGER, INTENT(IN)   :: nmctop
    REAL(wp),INTENT(IN)   :: cevapcu(:)
    REAL(wp),INTENT(INOUT):: pcon_dtrl(kbdim), pcon_dtri(kbdim)
    REAL(wp),INTENT(INOUT):: pcon_iqte(kbdim)
    REAL(wp),INTENT(INOUT):: ptte_cnv(kbdim,klev)
    REAL(wp),INTENT(INOUT):: pvom_cnv(kbdim,klev), pvol_cnv(kbdim,klev)
    REAL(wp),INTENT(INOUT):: pqte_cnv(kbdim,klev), pxtte_cnv(kbdim,klev,ktrac)
    REAL(wp),INTENT(IN)   :: pcen(kbdim,klev)
#endif
    REAL(wp),INTENT(IN)   :: pcpen(kbdim,klev), pqte(kbdim,klev)
    !
    REAL(wp):: pten(kbdim,klev),        pqen(kbdim,klev),                  &
      &        pxen(kbdim,klev),        ptven(kbdim,klev),                 &
      &        puen(kbdim,klev),        pven(kbdim,klev),                  &
      &        pqsen(kbdim,klev),       pgeo(kbdim,klev),                  &
      &        paphp1(kbdim,klevp1),                                       &
      &        pverv(kbdim,klev)                                            
    REAL(wp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
      &        plu(kbdim,klev),         plude(kbdim,klev),                 &
      &        pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
      &        prsfc(kbdim),            pssfc(kbdim),                      &
      &        prain(kbdim),            pqhfla(kbdim)
    REAL(wp):: pthvsig(kbdim)
    INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
      &        ktype(kbdim)
    REAL(wp):: pqtec(kbdim,klev),                                          &
      &        pxtecl(kbdim,klev),      pxteci(kbdim,klev),                &
      &        pqude(kbdim,klev)
    REAL(wp):: ztenh(kbdim,klev),       zqenh(kbdim,klev),                 &
      &        zxenh(kbdim,klev),       zalvsh(kbdim,klev),                &
      ! zalvsh: latent heat of vaporisation/sublimation defined at half levels
      &        zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                &
      &        ztd(kbdim,klev),         zqd(kbdim,klev),                   &
      &        zmfus(kbdim,klev),       zmfds(kbdim,klev),                 &
      &        zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                 &
      &        zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                &
      &        zmful(kbdim,klev),       zrfl(kbdim),                       &
      &        zuu(kbdim,klev),         zvu(kbdim,klev),                   &
      &        zud(kbdim,klev),         zvd(kbdim,klev)
    REAL(wp):: zcpcu(kbdim,klev)
    REAL(wp):: zentr(kbdim),            zhcbase(kbdim),                    &
      &        zmfub(kbdim),            zmfub1(kbdim),                     &
      &        zktype(kbdim),           zldcum(kbdim),                     &
      &        zcpcui(kbdim,klev),      zkcbot(kbdim),                     &
      &        zictop0(kbdim),          ztmp1(kbdim),                      &
      &        ztmp2(kbdim),            ztmp3(kbdim),                      &
      &        za(kbdim)
    REAL(wp):: zsfl(kbdim),             zdpmel(kbdim,klev)
    REAL(wp):: zcape(kbdim),            zheat(kbdim)
    REAL(wp):: ua(kbdim), dua(kbdim), ub(kbdim)
    REAL(wp):: zhmin(kbdim),            zihmin(kbdim)
    REAL(wp):: zhhatt(kbdim,klev),      zdqpbl(kbdim)
    REAL(wp):: zdqcv(kbdim)
    INTEGER :: ihmin(kbdim), ilo1(kbdim), ldidx(kbdim), loidx(kbdim)
    INTEGER :: ilab(kbdim,klev),        idtop(kbdim),                      &
      &        ictop0(kbdim),           ilwmin(kbdim)
    REAL(wp):: pxten(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),            &
      &        zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
      &        zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
    LOGICAL :: loddraf(kbdim),          ldland(kbdim)
    LOGICAL :: ldcum(kbdim)
    LOGICAL :: llo1
#ifndef __ICON__
    REAL(wp):: zmwc(kbdim,klev),        zmrateprecip(kbdim,klev)
#endif
    !
    INTEGER :: nl, jl, jk, ikb, jt, itopm2, locnt, ldcnt
    REAL(wp):: zcons2, zqumqe, zdqmin, zmfmax, zalvs, zalvdcp, zqalv       &
      &      , zhsat, zes, zcor, zqsat, zdqsdt, zgam, zzz, zhhat           &
      &      , zb, zbi, zroi, zdz, zdhdz, zdepth, zfac, zrh, zeps          &
      &      , zjk, zhelp, zlo1, zpaphp1i, ztenhi, zkctop, zpbmpt, za1, za2
    !
    !  INTRINSIC FUNCTIONS
    INTRINSIC MIN, MAX
    !
    !  Executable statements
    
#ifndef __ICON__
    IF (lconvmassfix) THEN
      CALL xt_conv_massfix(kproma,         kbdim,         klev,              &
        &                  klevp1,         ktrac,         krow,              &
        &                  papp1,          paphp1,        pxtte,             &
        &                  .TRUE.  )
    END IF
#endif
    !-----------------------------------------------------------------------
    !
    !     1.           Specify constants and parameters
    !                  --------------------------------
    !
    zcons2=1._wp/(grav*ptime_step_len)
    !
    !----------------------------------------------------------------------
    !
    !*    2.           Initialize values at vertical grid points in 'cuini'
    !                  ---------------------------------------------------
    !
    CALL cuini(kproma, kbdim, klev, klevp1, klevm1,                      &
      &        pten,     pqen,     pqsen,    pxen,     puen,     pven,   &
      &        ptven,    ktrac,                                          &
      &        pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt, &
      &        pverv,    pgeo,     paphp1,   zgeoh,                      &
      &        ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,           &
      &        ptu,      pqu,      ztd,      zqd,                        &
      &        zuu,      zvu,      zud,      zvd,                        &
      &        pmfu,     pmfd,     zmfus,    zmfds,                      &
      &        zmfuq,    zmfdq,    zdmfup,   zdmfdp,                     &
      &        pcpen,    zcpcu,    zalvsh,                               &
      &        zdpmel,   plu,      plude,    pqude,    ilab             )
    !
    !-----------------------------------------------------------------------
    !
    !*    3.0          Cloud base calculations
    !                  -----------------------
    !
    !*             (A) Determine cloud base values in 'cubase'
    !                  ---------------------------------------
    !
    CALL cubase(kproma, kbdim, klev, klevp1, klevm1,                     &
      &         ztenh,    zqenh,    zgeoh,    paphp1,    pthvsig,        &
      &         ptu,      pqu,      plu,                                 &
      &         puen,     pven,     zuu,      zvu,                       &
      &         zcpcu,                                                   &
      &         ldcum,    kcbot,    ilab)
    !
    !*             (B) Determine total moisture convergence and
    !*                 then decide on type of cumulus convection
    !                  -----------------------------------------
    !
    zkcbot(1:kproma) = REAL(kcbot(1:kproma),wp)

    jk=1
    DO jl=1,kproma
      zdqpbl(jl)=0.0_wp
      zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
      idtop(jl)=0
    END DO
    DO jk=2,klev
      zjk = REAL(jk,wp)
      DO jl=1,kproma
        zhelp      = paphp1(jl,jk+1)-paphp1(jl,jk)
        zdqcv(jl)  = zdqcv(jl)+pqte(jl,jk)*zhelp
        zdqpbl(jl) = zdqpbl(jl) + FSEL(zjk - zkcbot(jl),pqte(jl,jk),0._wp)*zhelp
#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,E18.10,L3)','zdqpbl',jk,jl,zdqpbl(jl),(FSEL(zjk - zkcbot(jl),1._wp,0._wp)==1._wp)
#endif
      END DO
    END DO
    !
    !*             (C) Determine moisture supply for boundary layer and determine 
    !*                 cloud base massflux ignoring the effects of downdrafts
    !*                 at this stage
    !                  ---------------------------------------------------------------
    !
    zldcum(1:kproma) = MERGE(1._wp,0._wp,ldcum(1:kproma))

!DIR$ IVDEP
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
    DO jl=1,kproma
      ikb=kcbot(jl)
      zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
      zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
      zlo1 = FSEL(-zdqpbl(jl),0._wp,1._wp)
      zlo1 = FSEL(zdqmin - zqumqe,0._wp,zlo1) * zldcum(jl)
      zmfub(jl)=FSEL(-zlo1,0.01_wp,(zdqpbl(jl)/(grav*MAX(zqumqe,zdqmin))))
      zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      zmfub(jl)=MIN(zmfub(jl),zmfmax)
      zldcum(jl) = zlo1
      zhelp = MAX(0._wp,-1.1_wp*pqhfla(jl)*grav)
      zktype(jl) = FSEL(zhelp - zdqcv(jl),2._wp,1._wp)
      zentr(jl)  = FSEL(zhelp - zdqcv(jl),entrscv,entrpen)
      ktype(jl) = INT(zktype(jl))
#ifdef __ibmdbg__
      PRINT '(A6,I3,I4,2 E18.10,I3,L3)','zmfub',ikb,jl,zmfub(jl),zentr(jl),ktype(jl),(zlo1==1._wp)
#endif
    END DO
    ldcum(1:kproma) = (zldcum(1:kproma).GT.0._wp)
    !
    !-----------------------------------------------------------------------
    !*    4.0          Determine cloud ascent for entraining plume
    !                  -------------------------------------------
    !
    !*   (A) Estimate cloud height for entrainment/detrainment calculations in
    !*       cuasc (max. possible cloud height for non-entraining plume,
    !*       following A.-S.,1974)
    !        -------------------------------------------------------------------------
    !
!DIR$ IVDEP
    DO jl=1,kproma
      ikb=kcbot(jl)
      zalvs=FSEL(tmelt-ptu(jl,ikb),als,alv)
      zhcbase(jl) = zcpcu(jl,ikb)*ptu(jl,ikb) + zgeoh(jl,ikb) + zalvs*pqu(jl,ikb)
      zictop0(jl) = zkcbot(jl)-1._wp
    END DO
    DO jk=klev,1,-1
      zcpcui(1:kproma,jk) = 1._wp/zcpcu(1:kproma,jk)
      IF (jk <= klevm1 .AND. jk >= 3) THEN
        ! mpuetz: too few instructions (FP dependencies)
        CALL prepare_ua_index_spline('cumastr',kproma,ztenh(1,jk),loidx(1),za(1))
        CALL lookup_ua_spline(kproma,loidx(1),za(1),ua(1),dua(1))
        CALL lookup_ubc('cumastr',kproma,ztenh(1,jk),ub(1))
        zjk = REAL(jk,wp)
!IBM* NOVECTOR
        DO jl=1,kproma
          ! mpuetz: move some of these into the previous loop
          zalvs=FSEL(tmelt - ztenh(jl,jk),als,alv)
          zalvdcp=zalvs*zcpcui(jl,jk)
          zqalv=1._wp/zalvs
          zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
          zpaphp1i = 1._wp/paphp1(jl,jk)
          zes=ua(jl)*zpaphp1i
          zes=MIN(0.5_wp,zes)
          zcor=1._wp/(1._wp-vtmpc1*zes)
          zqsat=zes*zcor
          zdqsdt=zpaphp1i*zcor**2*dua(jl)
          zgam=FSEL(zes - 0.4_wp,zqsat*zcor*ub(jl),zalvdcp*zdqsdt)
          zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
          zhhat=zhsat-((zzz+zgam*zzz)/(1._wp+zgam*zzz*zqalv))*             &
            &  MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._wp)
          zhhatt(jl,jk)=zhhat
          zlo1 = FSEL(zjk - zictop0(jl),0._wp,1._wp)
          zlo1 = FSEL(zhhat - zhcbase(jl),0._wp,zlo1)
          zictop0(jl) = FSEL(-zlo1,zictop0(jl),zjk)
#ifdef __ibmdbg__
          PRINT '(A6,I3,I4,3 E18.10,I3)','ictop',jk,jl,zes,zqsdt,zhhatt(jl,jk),INT(zictop0(jl))
#endif
        END DO
      END IF
    END DO
    ictop0(1:kproma) = INT(zictop0(1:kproma))
    !!
    !!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
    !!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
    !!
    !  DO jl=1,kproma
    !     ktype(jl)=MERGE(1,2,                                              &
    !                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_wp)
    !     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
    !  ENDDO
    !!
    !              Find lowest possible org. detrainment level
    !              -------------------------------------------
    !
    ldcnt = 0
    DO jl=1,kproma
      zhmin(jl)=0._wp
      zihmin(jl)=0._wp
      llo1=ldcum(jl).AND.ktype(jl).EQ.1
      IF(llo1) THEN
        zihmin(jl)=zkcbot(jl)
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
      ENDIF
    ENDDO
    !
    zb=25._wp
    zbi=1._wp/(zb*grav)
    DO jk=klev,1,-1
      ! mpuetz: compute the update criterion
      zjk = REAL(jk,wp)
!IBM* ASSERT(NODEPS)
      DO nl = 1,ldcnt
        jl = ldidx(nl)
        zlo1 = FSEL(zjk - zkcbot(jl),0._wp,1._wp)
        zlo1 = FSEL(zjk - zictop0(jl),zlo1,0._wp)
        zlo1 = FSEL(-ABS(zihmin(jl)-zkcbot(jl)),zlo1,0._wp)
        ilo1(nl) = INT(zlo1)
      END DO
      ! mpuetz: compute the indices of elements to be updated
      locnt = 0
      DO nl = 1,ldcnt
        IF (ilo1(nl).GT.0) THEN
          locnt = locnt + 1
          loidx(locnt) = ldidx(nl)
        END IF
      END DO
!IBM* ASSERT(NODEPS)
      DO nl=1,locnt
        jl = loidx(nl)
        zalvs = FSEL(tmelt - ztenh(jl,jk),als,alv)
      !  zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
        ikb   = kcbot(jl)
        zroi  = SWDIV_NOCHK(rd*ztenh(jl,jk)*(1._wp+vtmpc1*zqenh(jl,jk)),paphp1(jl,jk))
        zdz   = (paphp1(jl,jk)-paphp1(jl,jk-1))*zroi/grav
        za1 = (pcpen(jl,jk-1)*pten(jl,jk-1) - pcpen(jl,jk)*pten(jl,jk)               &
             + zalvs*(pqen(jl,jk-1) - pqen(jl,jk))+(pgeo(jl,jk-1)-pgeo(jl,jk)))*grav
        za2 = pgeo(jl,jk-1)-pgeo(jl,jk)
        zdhdz = SWDIV_NOCHK(za1, za2)
        zdepth    = zgeoh(jl,jk)-zgeoh(jl,ikb)
        ztmp1(nl) = zalvs
        ztmp2(nl) = zdz*zdhdz
        ztmp3(nl) = 1._wp+zdepth*zbi
      END DO
      ztmp3(1:locnt) = SQRT(ztmp3(1:locnt))
!IBM* ASSERT(NODEPS)
      DO nl=1,locnt
        jl = loidx(nl)
        zalvs     = ztmp1(nl)
        zfac      = ztmp3(nl)
        zdepth    = zgeoh(jl,jk)-zgeoh(jl,ikb)
        zhmin(jl) = zhmin(jl) + zfac*ztmp2(nl)
        zrh       =-zalvs*(zqsenh(jl,jk)-zqenh(jl,jk))*zfac
        zihmin(jl) = FSEL(zrh - zhmin(jl),zihmin(jl),zjk)
      !  IF(zhmin(jl).GT.zrh) ihmin(jl)=jk
      ENDDO
    ENDDO
!IBM* ASSERT(NODEPS)
    DO nl=1,ldcnt
      jl = ldidx(nl)
#ifdef __ibmdbg__
      PRINT '(A6,I4,I4,I4)','ihmin',jl,INT(zihmin(jl)),ictop0(jl)
#endif
      zihmin(jl) = FSEL(zihmin(jl)-zictop0(jl),zihmin(jl),zictop0(jl))
    !    IF(ihmin(jl).LT.ictop0(jl)) ihmin(jl)=ictop0(jl)
    ENDDO
    ihmin(1:kproma) = INT(zihmin(1:kproma))
    !
    !*         (B) Do ascent in 'cuasc' in absence of downdrafts
    !              ---------------------------------------------
    !
    CALL cuasc(kproma, kbdim, klev, klevp1, klevm1,                      &
      &        ztenh,    zqenh,    puen,     pven,                       &
      &        ktrac,                                                    &
      &        ptime_step_len,                                           &
      &        zxtenh,   pxten,    pxtu,     zmfuxt,                     &
      &        pten,     pqen,     pqsen,                                &
      &        pgeo,     zgeoh,    paphp1,   pthvsig,                    &
      &        pqte,               pverv,    ilwmin,                     &
      &        ldcum,    ldland,   ktype,    ilab,                       &
      &        ptu,      pqu,      plu,      zuu,      zvu,              &
      &        pmfu,     zmfub,    zentr,                                &
      &        zmfus,    zmfuq,                                          &
      &        zmful,    plude,    pqude,    zdmfup,                     &
      &        ihmin,    zhhatt,   zhcbase,  zqsenh,                     &
      &        pcpen,    zcpcu,                                          &
      &        kcbot,    kctop,    ictop0,                               &
#ifndef __ICON__
      &        zmwc,     zmrateprecip                      )
#else
      &        nmctop                                      )
#endif
    !
    !*     (C) Check cloud depth and change entrainment rate accordingly
    !          Calculate precipitation rate (for downdraft calculation)
    !          ---------------------------------------------------------
    !
!DIR$ IVDEP
    DO jl=1,kproma
      zpbmpt=paphp1(jl,kcbot(jl))-paphp1(jl,kctop(jl))
      IF(ldcum(jl).AND.ktype(jl).EQ.1.AND.zpbmpt.LT.2.e4_wp) ktype(jl)=2
      IF(ldcum(jl)) ictop0(jl)=kctop(jl)
      IF(ktype(jl).EQ.2) zentr(jl)=entrscv
      zrfl(jl)=zdmfup(jl,1)
    END DO
    DO jk=2,klev
      DO jl=1,kproma
        zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
      END DO
    END DO
    ! mpuetz: must recompute ldidx() since ktype could have changed
    ldcnt = 0
    DO jl=1,kproma
      llo1=ldcum(jl).AND.ktype(jl).EQ.1
      IF(llo1) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
      ENDIF
    ENDDO
    !
    !-----------------------------------------------------------------------
    !*    5.0          Cumulus downdraft calculations
    !                  ------------------------------
    !
    IF(lmfdd) THEN
      !
      !*             (A) Determine lfs in 'cudlfs'
      !                  -------------------------
      CALL cudlfs(kproma,   kbdim,    klev,     klevp1,                 &
        &         ztenh,    zqenh,    puen,     pven,                   &
        &         ktrac,                                                &
        &         zxtenh,   pxtu,     zxtd,     zmfdxt,                 &
        &         zgeoh,    paphp1,                                     &
        &         ptu,      pqu,      zuu,      zvu,                    &
        &         ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
        &         ztd,      zqd,      zud,      zvd,                    &
        &         pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
        &         zcpcu,                                                &
        &         idtop,    loddraf)
      !
      !*            (B)  Determine downdraft t,q and fluxes in 'cuddraf'
      !                  -----------------------------------------------
      CALL cuddraf(kproma,   kbdim,    klev,     klevp1,                &
        &          ztenh,    zqenh,    puen,     pven,                  &
        &          ktrac,                                               &
        &          zxtenh,   zxtd,     zmfdxt,                          &
        &          zgeoh,    paphp1,   zrfl,                            &
        &          ztd,      zqd,      zud,      zvd,                   &
        &          pmfd,     zmfds,    zmfdq,    zdmfdp,                &
        &          zcpcu,                                               &
        &          loddraf)
      !
    END IF
    !
    !*    5.5      Recalculate cloud base massflux from a cape closure for deep
    !*             convection (ktype=1) and by pbl equilibrium taking downdrafts
    !*             into account for shallow convection (ktype=2)
    !              -------------------------------------------------------------
    ! mpuetz: cuasc can modify kctop !!
    !
    zkcbot(1:kproma) = REAL(kcbot(1:kproma),wp)
    !
    DO jl=1,kproma
      zheat(jl)=0._wp
      zcape(jl)=0._wp
      zmfub1(jl)=zmfub(jl)
    ! zhelp(jl) =0._wp
    END DO
    !
    DO jk=1,klev
      zjk = REAL(jk,wp)
!IBM* ASSERT(NODEPS)
      DO nl = 1,ldcnt
        jl = ldidx(nl)
        zkctop = REAL(kctop(jl),wp)
        zlo1 = FSEL(zkcbot(jl) - zjk,1._wp,0._wp)* FSEL(zkctop - zjk,0._wp,1._wp)
        ilo1(nl) = INT(zlo1)
      END DO
      locnt = 0
      DO nl = 1,ldcnt
        jl = ldidx(nl)
        IF(jk.LE.kcbot(jl).AND.jk.GT.kctop(jl)) THEN
      !  IF (ilo1(nl).GT.0) THEN
          locnt = locnt + 1
          loidx(locnt) = jl
        END IF
      END DO
#ifdef __ibmdbg__
      PRINT '(A6,I3,2 I4)','cumas1',jk,ldcnt,locnt,INT(zkcbot(jl))
#endif
      ! mpuetz: there is reuse of zro, maybe we can calulate once and store
!IBM* ASSERT(NODEPS)
      DO nl=1,locnt
        jl = loidx(nl)
        zroi = SWDIV_NOCHK(rd*ztenh(jl,jk)*(1._wp+vtmpc1*zqenh(jl,jk)),paphp1(jl,jk))
        ztenhi = SWDIV_NOCHK(1._wp,ztenh(jl,jk))
        zdz    = (paphp1(jl,jk)-paphp1(jl,jk-1))*zroi/grav
        zheat(jl) = zheat(jl)+((pten(jl,jk-1)-pten(jl,jk) + grav*zdz*zcpcui(jl,jk))  &
          &         *ztenhi+ vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk)))                     &
          &         *(grav*(pmfu(jl,jk)+pmfd(jl,jk)))*zroi
        zcape(jl) = zcape(jl) + ( grav*(ptu(jl,jk)-ztenh(jl,jk))*ztenhi              &
                    +grav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))-grav*plu(jl,jk))*zdz
#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,5 E18.10,I4)','zheat',jk,jl,zheat(jl),zcape(jl),zroi,ztenhi,zdz,kcbot(jl)
        PRINT '(A6,I3,I4,6 E18.10)','zcape',jk,jl,pten(jl,jk-1)-pten(jl,jk),pqen(jl,jk-1)-pqen(jl,jk),pmfu(jl,jk)+pmfd(jl,jk),ptu(jl,jk),pqu(jl,jk),plu(jl,jk)
#endif
      ENDDO
    ENDDO
    !
    !  DO jl=1,kproma
    !     DO jk=2,klev
    !        zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*                         &
    !             (1._wp+vtmpc1*zqenh(jl,jk)))
    !        zdz=(paphp1(jl,jk)-paphp1(jl,jk-1))/(grav*zro)
    !        zhelp(jl)=zhelp(jl) +                               &
    !             (grav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)     &
    !             +grav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))      &
    !             -grav*plu(jl,jk) ) * zdz
    !     ENDDO
    !     if (zhelp(jl).lt.0._wp) zhelp(jl)=0._wp
    !  ENDDO
    DO jl=1,kproma
      IF(ldcum(jl).AND.ktype(jl).EQ.1) THEN
        ikb=kcbot(jl)
        zmfub1(jl) = SWDIV_NOCHK((zcape(jl)*zmfub(jl)),(zheat(jl)*cmftau))
        zmfub1(jl) = MAX(zmfub1(jl),0.001_wp)
        zmfmax     = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        zmfub1(jl) = MIN(zmfub1(jl),zmfmax)
      ENDIF
    ENDDO
    !
    !*      Recalculate convective fluxes due to effect of downdrafts on boundary
    !*      layer moisture budget for shallow convection (ktype=2)
    !       ---------------------------------------------------------------------
    ldcnt = 0
    DO jl=1,kproma
      llo1=ktype(jl).EQ.2
      IF(llo1) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
      ENDIF
    ENDDO
!DIR$ IVDEP
!IBM* ASSERT(NODEPS)
    DO nl=1,ldcnt
      jl = ldidx(nl)
      ikb=kcbot(jl)
      llo1=pmfd(jl,ikb).LT.0._wp.AND.loddraf(jl)
      zeps=MERGE(cmfdeps,0._wp,llo1)
      zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zeps*zqd(jl,ikb)-(1._wp-zeps)*zqenh(jl,ikb)
      zdqmin=MAX(0.01_wp*zqenh(jl,ikb),1.e-10_wp)
      zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      llo1=zdqpbl(jl).GT.0._wp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)                    &
        &                                          .AND.zmfub(jl).LT.zmfmax
      zmfub1(jl)=MERGE(zdqpbl(jl)/(grav*MAX(zqumqe,zdqmin)),zmfub(jl),llo1)
      zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                                         &
        &                          ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_wp*zmfub(jl))
#ifdef __ibmdbg__
      PRINT '(A6,I4,E18.10)','zmfub1',jl,zmfub1(jl)
#endif
    END DO
    ldcnt = 0
    DO jl=1,kproma
      IF (ldcum(jl)) THEN
        ldcnt = ldcnt + 1
        ldidx(ldcnt) = jl
      ENDIF
    ENDDO
    DO jk=1,klev
!IBM* ASSERT(NODEPS)
      DO nl=1,ldcnt
        jl = ldidx(nl)
        ztmp1(nl) = SWDIV_NOCHK(zmfub1(jl),MAX(zmfub(jl),1.e-10_wp))
        zfac      = ztmp1(nl)
        pmfd(jl,jk)   = pmfd(jl,jk)*zfac
        zmfds(jl,jk)  = zmfds(jl,jk)*zfac
        zmfdq(jl,jk)  = zmfdq(jl,jk)*zfac
        zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
#ifdef __ibmdbg__
        PRINT '(A6,I3,I4,E18.10)','zfac',jk,jl,ztmp1(nl)
#endif
      END DO
!IBM* unroll(4)
      DO jt=1,ktrac
!IBM* ASSERT(NODEPS)
        DO nl=1,ldcnt
          jl = ldidx(nl)
          zfac             = ztmp1(nl)
          zmfdxt(jl,jk,jt) = zmfdxt(jl,jk,jt)*zfac
        END DO
      END DO
      !
    END DO
    !
    !*       New values of cloud base mass flux
    !        ----------------------------------
    !
    DO nl=1,ldcnt
      jl = ldidx(nl)
      zmfub(jl) = zmfub1(jl)
    END DO
    !
    !-----------------------------------------------------------------------
    !*    6.0   Determine final cloud ascent for entraining plume
    !*          for penetrative convection (type=1),
    !*          for shallow to medium convection (type=2)
    !*          and for mid-level convection (type=3).
    !           --------------------------------------------------
    !
    CALL cuasc(kproma, kbdim, klev, klevp1, klevm1,                      &
      &        ztenh,    zqenh,    puen,     pven,                       &
      &        ktrac,                                                    &
      &        ptime_step_len,                                           &
      &        zxtenh,   pxten,    pxtu,     zmfuxt,                     &
      &        pten,     pqen,     pqsen,                                &
      &        pgeo,     zgeoh,    paphp1,   pthvsig,                    &
      &        pqte,               pverv,    ilwmin,                     &
      &        ldcum,    ldland,   ktype,    ilab,                       &
      &        ptu,      pqu,      plu,      zuu,      zvu,              &
      &        pmfu,     zmfub,    zentr,                                &
      &        zmfus,    zmfuq,                                          &
      &        zmful,    plude,    pqude,    zdmfup,                     &
      &        ihmin,    zhhatt,   zhcbase,  zqsenh,                     &
      &        pcpen,    zcpcu,                                          &
      &        kcbot,    kctop,    ictop0,                               &
#ifndef __ICON__
      &        zmwc,     zmrateprecip                     )
#else
      &        nmctop                                     )
#endif
    !-----------------------------------------------------------------------
    !*    7.0      Determine final convective fluxes in 'cuflx'
    !              --------------------------------------------
    !
    CALL cuflx(kproma,   kbdim,    klev,     klevp1,                     &
      &        pqen,     pqsen,    ztenh,    zqenh,                      &
      &        ktrac,                                                    &
      &        ptime_step_len,                                           &
#ifndef __ICON__
      &        krow,                                                     &
      &        pxtte,    pxtu,     ptu,                                  &
      &        zmwc,     zmrateprecip,                                   &
#else
      &        cevapcu,                                                  &
#endif
      &        zxtenh,   zmfuxt,   zmfdxt,                               &
      &        paphp1,   zgeoh,                                          &
      &        kcbot,    kctop,    idtop,                                &
      &        ktype,    loddraf,  ldcum,                                &
      &        pmfu,     pmfd,     zmfus,    zmfds,                      &
      &        zmfuq,    zmfdq,    zmful,                                &
      &        zdmfup,   zdmfdp,   zrfl,     prain,                      &
      &        zcpcu,                                                    &
      &        pten,     zsfl,     zdpmel,   itopm2                     )
    !-----------------------------------------------------------------------
    !
    !*    8.0      Update tendencies for t and q in subroutine'cudtdq'
    !              ---------------------------------------------------
    !
    CALL cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, ktrac,       &
#ifndef __ICON__
      &         krow,                                                    &
      &         ptte,     pqte,                                          &
      &         pxtte,                                                   &
      &         paprc,    paprs,                                         &
#endif
      &         paphp1,   pten,                                          &
      &         zmfuxt,   zmfdxt,                                        &
      &         zmfus,    zmfds,    zmfuq,    zmfdq,                     &
      &         zmful,    zdmfup,   zdmfdp,   plude,                     &
      &         zdpmel,   zrfl,     zsfl,                                &
#ifndef __ICON__
      &         pcpen,    zalvsh,   pqtec,    pqude,                     &
#else
      &         pcen,     zalvsh,   pqtec,    pqude,                     &
#endif
      &         prsfc,    pssfc,                                         &
      &         pch_con,  pcw_con,                                       &
#ifdef __ICON__
      &         pcon_dtrl, pcon_dtri, pcon_iqte,                         &
      &         ptte_cnv, pqte_cnv, pxtte_cnv,                           &
#endif
      &         pxtecl,   pxteci                                        )
    !
    !-----------------------------------------------------------------------
    !
    !*    9.0      Update tendencies for u and v in subroutine'cududv'
    !              ---------------------------------------------------
    !
    IF(lmfdudv) THEN
      CALL cududv(kproma,   kbdim,    klev,     klevp1,                 &
        &         itopm2,   ktype,    kcbot,    paphp1,   ldcum,        &
#ifndef __ICON__
        &         puen,     pven,     pvom,     pvol,                   &
#else
        &         puen,     pven,     pvom_cnv, pvol_cnv,               &
#endif
        &         zuu,      zud,      zvu,      zvd,                    &
        &         pmfu,     pmfd)
      !
    END IF
    !
#ifndef __ICON__
    IF (lconvmassfix) THEN
      CALL xt_conv_massfix(kproma,         kbdim,         klev,              &
        &                  klevp1,         ktrac,         krow,              &
        &                  papp1,          paphp1,        pxtte,             &
        &                  .FALSE.  )
    END IF
#endif
  END SUBROUTINE cumastr

END MODULE mo_cumastr
