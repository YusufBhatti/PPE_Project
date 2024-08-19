#ifdef __xlC__
@PROCESS HOT
#else
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE vdiff ( kproma, kbdim, ktdia, klev, klevm1, klevp1, ktrac   &
         , jrow                                                        &
!-----------------------------------------------------------------------
! - 3D from mo_memory_g1a
         , pxtm1                                                       &
! - 2D from mo_memory_g1a
         , pqm1,           ptm1,           pum1                        &
         , pvm1,           pxlm1,          pxim1                       &
! - 1D 
         , pahfl,          pahfs,          paz0                        &
         , pdew2,          pevap,          pforest                     &
         , ptemp2,         pt2max                                      &
         , pt2min,         pwind10w,       pvdis                       &
         , pu10,           pv10,           pustr                       &
         , pvstr,          pwimax,         pwind10                     &
!-----------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition  -----
         , pvgrat,         pvlt,           pcvw                        &
!---------- kai zhang 2009-07-01  for emission and dry deposition  -----
!-----------------------------------------------------------------------
         , ptsw,           ptsi                                        &
         , pocu,           pocv                                        &
         , paz0l,          paz0w,          paz0i                       &
         , pahfsl,         pahfsw,         pahfsi                      &
         , pahflw,         pevapw,         pevapi                      &
         , pahfslac,       pahfswac,       pahfsiac                    &
         , pahfllac,       pahflwac,       pahfliac                    &
         , pevaplac,       pevapwac,       pevapiac                    &
         , pustrl,         pustrw,         pustri                      &
         , pvstrl,         pvstrw,         pvstri                      &
         , palbedo,        palbedo_vis                                 &
         , palbedo_nir,    palsol                                      &
         , palbedo_vis_dir,                palbedo_nir_dir             &
         , palbedo_vis_dif,                palbedo_nir_dif             &
!--- Included PBL top level export (Jan Kazil 10/2008)------------------
         , ppbl                                                        &
!--- End included ------------------------------------------------------
! - 2D from mo_memory_g3b
         , ptke,           ptkem1,         ptkem                       &
         , paclc,          pemter                                      &
         , pthvvar,        pthvsig                                     &
         , psh_vdiff,      pev_vdiff                                   &
! - 2D within physics only
         , paphm1,         papm1,          pgeom1                      &
         , ptvm1                                                       &
! - 1D within physics only.
         , pqhfla,         pevapot                                     &
         , ptslnew,        pfrl,           lpland                      &
         , lpnorth,        pfrw,           pfri                        &
         , palake                                                      &
! - Tendencies
! - 3D
         , pxtte                                                       &
! - 2D
         , pvol,           pvom,           pqte                        &
         , ptte,           pxlte,          pxite                       &
         , psiced                                                      &
! - jsbach addons for whole surface calculations
         , prain,          psnow                                       &
         , jsswnir,        jsswdifnir                                  &
         , jsswvis,        jsswdifvis                                  &
         , jsswpar,        jsswdifpar                                  &
         , pi0,            ptrsol                                      &
         , ptrflw,         ptrfli,         pahfli                      &
         , psofll,         psoflw,         psofli                      &
         , ptrfllac,       ptrflwac,       ptrfliac                    &
         , psofllac,       psoflwac,       psofliac                    &
         , palsoi,         palsow,         zth_new                     &
         , palsobs,        palsom,         ptaus                       &
         , pfage,          psnifrac,       pbarefrac                   &
         , pameltdepth,    pameltfrac                                 &
         , pamlcorr,       pamlcorac,      pamlheatac                  &
         , psni,           pahfice,        pfluxres                    &
         , pqres,          pahfcon,        pahfres                     &
         , pzteffl4,       pztsnew,        ptsurf                      &
         , pradtemp_old                                                &
!--- Included for emissions in chem_bcond (Philip Stier 08/04/01)-------
!-----------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ------
         , paphp1                                                      & 
!---------- kai zhang 2009-07-01  for emission and dry deposition ------
!-----------------------------------------------------------------------
         , pseaice                    )
!--- End Included for for emissions ------------------------------------
                              
!
!**** *vdiff* - does the vertical exchange of u, v, t, q, xl, xi and xt
!               by turbulence.
!
!
!     Subject.
!
!       This routine computes the physical tendencies of the seven
!   prognostic variables u, v, t, q, xl, xi and xt due to the vertical
!   exchange by turbulent (= non-moist convective) processes.
!   These tendencies are obtained as the difference between
!   the result of an implicit time-step starting from values at t-1
!   and these t-1 (???) values.
!   all the diagnostic computations (exchange coefficients, ...) are
!   done from the t-1 values. As a by-product the roughness length
!   over sea is updated accordingly to the *charnock formula. heat and
!   moisture surface fluxes and their derivatives against ts and ws,
!   later to be used for soil processes treatment, are also
!   computed as well as a stability value to be used as a diagnostic
!   of the depth of the well mixed layer in convective computations.
!
!**   Interface.
!     ----------
!
!          *vdiff* is called from *physc*.
!
!      Arguments.
!      ----------
!
!  - 3d from mo_memory_g1a
!
!  pxtm1    : tracer variables (t-dt)
!
!  - 2d from mo_memory_g1a
!
!  pqm1     : humidity (t-dt)
!  ptm1     : temperature (t-dt)
!  pum1     : zonal wind (t-dt)
!  pvm1     : meridional wind (t-dt)
!  pxlm1    : cloud water (t-dt)
!  pxim1    : cloud ice (t-dt)

! - 2d from mo_memory_g3
!
!  paclc    : cloud cover
!
!  ptke     : turbulent kinetic energy at t+dt (unfiltered)
!  ptkem    :            "             at t    (unfiltered)
!  ptkem1   :            "             at t-dt   (filtered)
!  pthvvar  : variance of virtual potential temperature
!
!  - 1d from mo_memory_g3
!
!  pthvsig  : std dev of virtual pot temp at standard half level klevm1
!
!  ptsw     :             "       over water
!  ptsi     :             "       over ice
!
!  pocu     : ocean u-velocity
!  pocv     : ocean v-velocity
!
!  pahfs    : surface sensible heat flux (accumulated)
!  pahfsl   :             "              over land
!  pahfsw   :             "              over water
!  pahfsi   :             "              over ice
!
!  pahfl    : surface latent heat flux   (accumulated)
!  pahfll   :             "              over land
!  pahflw   :             "              over water
!  pahfli   :             "              over ice
!
!  pevap    : surface evaporation (accumulated)
!  pevapl   :             "        over land
!  pevapw   :             "        over water
!  pevapi   :             "        over ice
!
!  paz0     : roughness length
!  paz0l    :      "            over land
!  paz0w    :      "            over water
!  paz0i    :      "            over ice
!
!  pustr    : u-stress (accumulated)
!  pustrl   :     "     over land
!  pustrw   :     "     over sea
!  pustri   :     "     over ice
!
!  pvstr    : v-stress (accumulated)
!  pvstrl   :     "     over land
!  pvstrw   :     "     over water
!  pvstri   :     "     over ice
!
!  pdew2    : dew point temperature at 2 meter
!  peforest : forest coverage
!  ptemp2   : temperature at 2 meter
!  ptsm1    : surface temperature (t-dt)
!  pt2max   : maximum temp. at 2 m between output intervals
!  pt2min   : minimun temp. at 2 m between output intervals
!  pwind10w : 10m wind over water
!  pu10     : u-wind at 10 meter
!  pv10     : v-wind at 10 meter
!  pwind10  : wind speed at 10 meter (accumulated)
!  pwimax   : maximum windspeed at 10 m. between output intervals
!  pvdis    : boundary layer dissipation (accumulated)
!
! - 2d within physics only
!
!  paphm1   : half level pressure (t-dt)
!  papm1    : full level pressure (t-dt)
!  ptvm1    : virtual temperature at t-dt
!
! - 1d within physics only
!
!  pgeom1   : geopotential above surface (t-dt)
!  pqhfla   : moisture flux at the surface
!  pevapot  : potential evaporation
!  ktropo   : tropopause index
!  lpland   : land-sea flag
!
!        Tendencies
!
!  - 3d
!
!  pxtte    : tendencies of tracer variables
!
!  - 2d
!  pvol     : tendency of meridional wind
!  pvom     : tendency of zonal wind
!  pqte     : tendency of humidity
!  ptte     : tendency of temperature
!  pxlte    : tendency of cloud water
!  pxite    : tendency of cloud ice
!
!
!     Method.
!     -------
!
!        First an auxialiary variable cp(q)t+gz is created on which
!   the vertical diffusion process will work like on u,v and q. then
!   along the vertical and at the surface, exchange coefficients (with
!   the dimension of a pressure thickness) are computed for momentum
!   and for heat (sensible plus latent). the letters m and h are used
!   to distinguish them. the diffusioncoefficents depend on the
!   turbulent kinetic energy (tke) calculated by an additional
!   prognostic equation, which considers advection of tke.
!        In the second part of the routine the implicit linear
!   systems for u,v first and t,q second are solved by a *gaussian
!   elimination back-substitution method. for t and q the lower
!   boundary condition depends on the surface state.
!   for tke the lower boundary condition depends on the square of
!   the frictional velocity.
!   over land, two different regimes of evaporation prevail:
!   a stomatal resistance dependent one over the vegetated part
!   and a soil relative humidity dependent one over the
!   bare soil part of the grid mesh.
!   potential evaporation takes place over the sea, the snow
!   covered part and the liquid water covered part of the
!   grid mesh as well as in case of dew deposition.
!        Finally one returns to the variable temperature to compute
!   its tendency and the later is modified by the dissipation's effect
!   (one assumes no storage in the turbulent kinetic energy range) and
!   the effect of moisture diffusion on cp. z0 is updated and the
!   surface fluxes of t and q and their derivatives are prepared and
!   stored like the difference between the implicitely obtained
!   cp(q)t+gz and cp(q)t at the surface.
!
!
!     Reference.

!
!          See vertical diffusion's part of the model's documentation
!     for details about the mathematics of this routine.
!
!     Authors.
!
!     u. schlese     dkrz-hamburg  feb-93
!       modified     e. roeckner  - 1994
!
!     j.-p. schulz   mpi - 1997 : implementation of implicit
!                                 coupling between land surface
!                                 and atmosphere.
!     m. esch, mpi, june 1999, echam5-modifications
!
!
!     based  on  original ecmwf version by j.f. geleyn  - 1982
!                              modified by c.b. blondin - 1986
!                                          h. feichter  - 1991
!                                          s. brinkop   - 1992
!                                          m. claussen  - 1993
USE mo_kind,             ONLY: wp
USE mo_geoloc,           ONLY: coriol_2d, amu0_x
USE mo_param_switches,   ONLY: lvdiff
USE mo_physc2,           ONLY: clam, ckap, cb, cc, cvdifts
USE mo_physical_constants, ONLY: grav, vtmpc1, rd, tmelt, alv, als,    &
                                 cpd, vtmpc2
USE mo_tracdef,          ONLY: trlist
USE mo_echam_convect_tables,   ONLY: prepare_ua_index_spline, lookup_ua_spline
USE mo_exception,        ONLY: finish
USE mo_time_control,     ONLY: delta_time, lstart, time_step_len
USE mo_semi_impl,        ONLY: eps
USE mo_memory_g3b,       ONLY: tkevn, cdum, cduh, udif, vdif, ustarm
USE mo_surface,          ONLY: update_surface
USE mo_co2,              ONLY: ico2idx,            & ! Index of CO2 in tracer list
                               lco2_mixpbl,        & ! Mix CO2 in PBL?
                               co2m1,              & ! CO2 concentration (for lco2=false)
                               co2_flux_ocean,     & ! passed to update_surface
                               co2_flux_land,      & !          "
                               co2_flux,           & !          "
                               co2_flux_npp,       & !          "
                               co2_flux_soilresp,  & !          "
                               co2_flux_herbivory, & !          "
                               co2_flux_dynveg,    & !          "
                               co2_emission_lcc,   & !          "
                               co2_emission_harvest  !          "
!
USE mo_control,          ONLY: vct, nvclev, lnwp, lrce
USE mo_submodel,         ONLY: lanysubmodel
USE mo_submodel_interface, ONLY: vdiff_subm

!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------

!--- Included for emissions in xt_emiss (Philip Stier 08/04/01)----------
USE mo_vphysc,           ONLY: vphysc
USE mo_submodel,         ONLY: lanysubmodel
!--- End Included for emissions in xt_emiss -----------------------------

!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------

#ifdef _PROFILE
USE mo_profile,          ONLY: trace_start, trace_stop
#endif
!
IMPLICIT NONE

!
! Variables declared
!

  ! gauss grid description
  INTEGER, INTENT(in)                      :: jrow   ! sequential index

  INTEGER :: itop, itopp1, jk, jl, jt, klev, klevm1, klevp1
  INTEGER :: iblmax, iblmin, ktop
  INTEGER :: kproma, kbdim, ktdia, ktrac
  REAL(wp):: z2geomf, zalf, zalh2
  REAL(wp):: zb, zbet, zbuoy, zc
  REAL(wp):: zchneu
  REAL(wp):: zcons2, zcons3, zcons5, zcons6, zcons8, zcons9
  REAL(wp):: zcons13, zcons15, zcons18, zcons23, zcons25
  REAL(wp):: zcor, zcpd 
  REAL(wp):: zda1, zdisc, zdisl, zdisx, zdisxt, zdivv, zdivv1
  REAL(wp):: zdqdt, zdqtot, zds, zdtdt, zdudt, zdus1
  REAL(wp):: zdus2, zdvdt, zdximdt, zdxlmdt, zdz
  REAL(wp):: zepcor, zepdu2, zeps, zepshr
  REAL(wp):: zes, zfac, zfox
  REAL(wp):: zfux, zh1, zh2, zhexp
  REAL(wp):: zkap, zkappa, zktest, zlam, zm1, zm2
  REAL(wp):: zm4, zmix, zmult1, zmult2, zmult3, zmult4
  REAL(wp):: zmult5, zqddif, zqdp, zqtmit
  REAL(wp):: zrd, zrdrv, zri, zrvrd, zsdep1
  REAL(wp):: zsdep2, zsh, zshear, zshn, zsm, zsmn
  REAL(wp):: zteldif, ztest, ztkemin, ztkesq
  REAL(wp):: ztmst, ztpfac1, ztpfac2, ztpfac3, ztpfac4
  REAL(wp):: zucf, zusus1
  REAL(wp):: zz2geo, zzb, zztvm, zzzlam, zdtime
  REAL(wp):: ztvh, zdisv, zfav
  REAL(wp):: zdisthv, zfacthv, zthvprod, zthvdiss, zthvirdif

!     ------------------------------------------------------------------
!
  REAL(wp):: pxtm1(kbdim,klev,ktrac)
!
  REAL(wp)::                                                           &
       pqm1(kbdim,klev),     ptm1(kbdim,klev),     pum1(kbdim,klev)    &
      ,pvm1(kbdim,klev),     pxlm1(kbdim,klev),    pxim1(kbdim,klev)
  REAL(wp)::                                                           &
       pahfl(kbdim),         pahfs(kbdim),         paz0(kbdim)         &
      ,pdew2(kbdim),         pevap(kbdim),         pforest(kbdim)      &
      ,ptemp2(kbdim),        pt2max(kbdim)                             &
      ,pt2min(kbdim),        pwind10w(kbdim),      pvdis(kbdim)        &
      ,pu10(kbdim),          pv10(kbdim),          pustr(kbdim)        &
      ,pvstr(kbdim),         pwimax(kbdim),        pwind10(kbdim)      &
      ,ptsw(kbdim),          ptsi(kbdim)                               &
      ,pocu(kbdim),          pocv(kbdim)                               &
      ,paz0l(kbdim),         paz0w(kbdim),         paz0i(kbdim)        &
      ,pahfsl(kbdim),        pahfsw(kbdim),        pahfsi(kbdim)       &
      ,pahflw(kbdim),        pevapw(kbdim),        pevapi(kbdim)       &
      ,pahfslac(kbdim),      pahfswac(kbdim),      pahfsiac(kbdim)     &
      ,pahfllac(kbdim),      pahflwac(kbdim),      pahfliac(kbdim)     &
      ,pevaplac(kbdim),      pevapwac(kbdim),      pevapiac(kbdim)     &
      ,pustrl(kbdim),        pustrw(kbdim),        pustri(kbdim)       &
      ,pvstrl(kbdim),        pvstrw(kbdim),        pvstri(kbdim)       &
      ,palbedo(kbdim),       palbedo_vis(kbdim)                        &
      ,palbedo_nir(kbdim),   palsol(kbdim)                             &
      ,palbedo_vis_dir(kbdim),     palbedo_nir_dir(kbdim)              &
      ,palbedo_vis_dif(kbdim),     palbedo_nir_dif(kbdim)
  REAL(wp)::                                                           &
       ptke(kbdim,klev),     ptkem1(kbdim,klev),   ptkem(kbdim,klev)   &
      ,paclc(kbdim,klev),    pemter(kbdim,klevp1), psh_vdiff(kbdim)    &
      ,pthvvar(kbdim,klev),  pthvsig(kbdim),       pev_vdiff(kbdim)
  REAL(wp)::                                                           &
       paphm1(kbdim,klevp1), papm1(kbdim,klev),    pgeom1(kbdim,klev)  &
      ,ptvm1(kbdim,klev)
  REAL(wp)::                                                           &
       pqhfla(kbdim),        pevapot(kbdim)                            &
      ,ptslnew(kbdim)
  LOGICAL ::                 lpland(kbdim)
  LOGICAL ::                 lpnorth(kbdim)
!--- Included PBL top level export (Jan Kazil 10/2008)------------------
  REAL(wp) ::                                                          &
       ppbl(kbdim)
!--- End included ------------------------------------------------------
  REAL(wp)::                                                           &
       pfrl(kbdim),         pfrw(kbdim),           pfri(kbdim)         &
      ,palake(kbdim)
!
  REAL(wp)::                                                           &
       pxtte(kbdim,klev,ktrac)
!
  REAL(wp)::                                                           &
       pvol(kbdim,klev),     pvom(kbdim,klev),     pqte(kbdim,klev)    &
      ,ptte(kbdim,klev),     pxlte(kbdim,klev),    pxite(kbdim,klev)   &
      ,psiced(kbdim)                                                   &
! rain and snow total [kg/m2*s] over last time step 
      ,prain(kbdim),         psnow(kbdim)                              &
      ,jsswnir(kbdim),       jsswdifnir(kbdim)                         &
      ,jsswvis(kbdim),       jsswdifvis(kbdim)                         &
      ,jsswpar(kbdim),       jsswdifpar(kbdim)                         &
      ,pi0(kbdim),           ptrsol(kbdim,klevp1)                      &
      ,ptrflw(kbdim),        ptrfli(kbdim)                             &
      ,psofll(kbdim),        psoflw(kbdim),        psofli(kbdim)       &
      ,ptrfllac(kbdim),      ptrflwac(kbdim),      ptrfliac(kbdim)     &
      ,psofllac(kbdim),      psoflwac(kbdim),      psofliac(kbdim)     &
      ,palsoi(kbdim),        palsow(kbdim)                             &
      ,palsobs(kbdim),       palsom(kbdim),        ptaus(kbdim)        &
      ,pfage(kbdim),         psnifrac(kbdim),      pbarefrac(kbdim)    &
      ,pameltdepth(kbdim),   pameltfrac(kbdim)                         &
      ,pamlcorr(kbdim),      pamlcorac(kbdim),     pamlheatac(kbdim)   &
      ,pseaice_new(kbdim),   psiced_new(kbdim)

  REAL(wp)::                                                           &
       psni(kbdim),          pahfice(kbdim),       pfluxres(kbdim)     &
      ,pqres(kbdim),         pahfcon(kbdim),       pahfres(kbdim)      &
      ,pzti(kbdim)                                                     &
      ,pzteffl4(kbdim),      pztsnew(kbdim),       ptsurf(kbdim)       &
      ,ptsw_new(kbdim)                                                 &
      ,pradtemp_old(kbdim)
!
!--- Included for dry deposition in xt_drydep (Philip Stier 08/04/01)----
  REAL(wp)::                                                           &
       pseaice(kbdim)
!--- End Included for for dry deposition in xt_drydep--------------------

!     local variables
!
  REAL(wp):: zxt(kbdim,klev,ktrac)
  REAL(wp)::                                                           &
          zxtdif(kbdim,klev,ktrac), zxtems(kbdim,ktrac)
  REAL(wp):: zdxtdt(kbdim,klev), zdxtdtsum, zdxtdtmean
  REAL(wp)::                                                           &
          zcfm(kbdim,klev),    zdis(kbdim,klev)                        &
         ,zcfh(kbdim,klev),    zcptgz(kbdim,klev),   zebsm(kbdim,klev) &
         ,zudif(kbdim,klev),   zvdif(kbdim,klev)                       &
         ,ztcoe(kbdim)
  REAL(wp)::                                                           &
          ztdif(kbdim,klev)                                            &
         ,zqdif(kbdim,klev),   zebsh(kbdim,klev),    zvidis(kbdim)
  REAL(wp)::                                                           &
          zhdyn(kbdim),        zteta1(kbdim,klev)                      &
         ,zlteta1(kbdim,klev), zcpten(kbdim,klev),   zqten(kbdim,klev) &
         ,ztvir1(kbdim,klev),  zhh(kbdim,klevm1),    zqss(kbdim,klev)  &
         ,zxldif(kbdim,klev),  zxidif(kbdim,klev),   zedif(kbdim,klev) &
         ,ztkevn(kbdim,klev),  zx(kbdim,klev)
  REAL(wp)::                                                           &
          zqssm(kbdim,klevm1), ztmitte(kbdim,klevm1)                   &
         ,zqmit(kbdim,klevm1), ztvirmit(kbdim,klevm1)                  &
         ,zfaxen(kbdim,klevm1),zfaxe(kbdim,klev)                       &
         ,ztemit(kbdim,klevm1),zccover(kbdim,klevm1)                   &
         ,zcdum(kbdim,klev),   zlwcmit(kbdim,klevm1)                   &
         ,zcfv(kbdim,klev),    zebsv(kbdim,klev)                       &
         ,zcthv(kbdim,klev),   zebthv(kbdim,klev)                      &
         ,zthvdif(kbdim,klev), zthvvar(kbdim,klev)                     &
         ,ztthv(kbdim),        zua(kbdim)                              &
         ,zpapm1i(kbdim),      za(kbdim) 

  REAL(wp):: zghabl(kbdim)
  REAL(wp):: zph(klevp1), zp(klev), zh(klev)
  
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------

!++mgs: list cleaned !
  REAL(wp):: pvlt(kbdim)                 ! leaf area index
  REAL(wp):: pvgrat(kbdim)               ! vegetation ratio
  REAL(wp):: pcvw(kbdim)                 ! wet skin fraction
  REAL(wp):: paphp1(kbdim,klevp1)        ! air pressure at layer interface (t+dt)

  REAL(wp):: zsrfll(kbdim)               ! surface net solar radiation flux over land (W/m2)
  REAL(wp):: zvelo10m(kbdim)             ! 10m wind     !++mgs: renamed from v10m
  REAL(wp):: ztslm1(kbdim)               ! surface temperature (t-dt)
  REAL(wp):: zcvs(kbdim)                 ! snow cover fraction
  REAL(wp):: zcfnc(kbdim)               ! function of heat transfer coeff.
  REAL(wp):: zrib(kbdim)                 ! moist richardson number
  REAL(wp):: ztvl(kbdim)                 ! see mo_surface_land
  REAL(wp):: zcdnl(kbdim)                ! see mo_surface_land
!--mgs

!---------- kai zhang 2009-07-01  for emission and dry deposition ----------  
!---------------------------------------------------------------------------
 
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-03  for sat simulator               ----------
  REAL(wp):: pahfli(kbdim)
  REAL(wp):: pahfll(kbdim) 
!---------- kai zhang 2009-07-03  for sat simulator               ----------
!---------------------------------------------------------------------------
 
!
  INTEGER :: ihpbl(kbdim),  ihpblc(kbdim),   ihpbld(kbdim), idx(kbdim)
!
!     ------------------------------------------------------------------
!
!     THE FOLLOWING VARIABLES ARE NEEDED FOR THE SUBROUTINES UPDATE-SURFACE
!

  REAL(wp)::    ztdif_new(kbdim), zqdif_new(kbdim)
  REAL(wp)::    zqsurf_new(kbdim), ztvh_new(kbdim), zth_new(kbdim)
  REAL(wp)::    ztte_corr(kbdim)

  REAL(wp)::    zco2(kbdim)     ! CO2 concentration [kg/kg] in lowest level

!>>SF gf #78
#ifdef HAMMOZ
  REAL(wp) :: zcfml(kbdim),   zcfmw(kbdim),    zcfmi(kbdim) !< 
  REAL(wp) :: zcfncl(kbdim),  zcfncw(kbdim),   zcfnci(kbdim)
  REAL(wp) :: zril(kbdim),    zriw(kbdim),     zrii(kbdim)
  REAL(wp) ::                 ztvw(kbdim),     ztvi(kbdim)
  REAL(wp) :: zcdn(kbdim),    zcdnw(kbdim),    zcdni(kbdim)
#endif
!<<SF gf #78

!  Executable statements 
!
#ifdef _PROFILE
  CALL trace_start ('vdiff', 20)
#endif

!
!     ------------------------------------------------------------------
!
!*    PHYSICAL CONSTANTS.
!
!          *ZLAM* IS THE ASYMPTOTIC MIXING LENGTH FOR MOMENTUM EXCHANGE,
!     *ZKAP* IS THE VON KARMAN CONSTANT, *ZB*, *ZC* AND *ZD* ARE SOME
!     CONSTANTS FOR THE FORMULAE ABOUT STABILITY DEPENDENCY RESPECTIVELY
!     NEAR THE NEUTRAL CASE, IN THE UNSTABLE CASE AND IN THE STABLE
!     CASE AND *ZCHAR* IS THE CONSTANT OF THE *CHARNOCK FORMULA.
!     *ZQWSSAT* AND *ZQSNCR* ARE THE INVERSES OF CRITICAL VALUES FOR
!     SOIL WATER AND SNOW DEPTH THAT ARE USED IN THE COMPUTATION OF THE
!     EVAPOTRANSPIRATION'S EFFICIENCY.
!
  zlam=clam
  zkap=ckap
  zb=cb
  zc=cc
  zrvrd=vtmpc1+1._wp
  zrdrv=1._wp/zrvrd
  zcpd=cpd
  zrd=rd
  zkappa=zrd/zcpd
!
!*    SECURITY PARAMETERS.
!
  zepdu2=1.0_wp
  zepshr=1.e-5_wp
  zepcor=5.e-05_wp
  ztkemin=1.e-10_wp
!
!*    COMPUTATIONAL CONSTANTS.
!
  zdtime = delta_time
  ztmst  = time_step_len
  ztpfac1=cvdifts
  ztpfac2=1._wp/ztpfac1
  ztpfac3=1._wp-ztpfac2
  ztpfac4=1._wp+ztpfac3
  zzzlam=1._wp
  zcons2=0.5_wp*zkap/grav
  zcons3=zlam
  zcons5=3._wp*zb*zc*grav**2
  zcons6=1._wp/3._wp
  zcons8=2._wp*zb
  zcons9=3._wp*zb
  zcons13=1._wp/ztmst
  zcons15=1._wp/(grav*ztmst)
  zcons18=ztpfac1*ztmst*grav**2
  zcons25=zcons2/zcons3
  zchneu=.3_wp
  zh1= 2.22_wp
  zh2= 0.22_wp
  zm1= 1.24_wp
  zm2= 2.37_wp
  zm4= 3.69_wp
  zshn=zh1*zh2*SQRT(2._wp)
  zsmn=zshn*zm1*zm2/zm4
  zda1=1._wp/zsmn**3
!
  itop=1
  itopp1=itop+1


! store partially updated tracer mixing ratio for tendency limitation
  DO jt = 1, ktrac
    zxt(1:kproma,:,jt) = pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt)*time_step_len
  END DO

!
! ------------------------------------------------------------
!
  ! CO2 mixing in PBL.
  ! Compute lowest (klev-2, approx. 500m in L19/31) and highest
  !   (highest below 2km) level ktop.
  ! CO2 is ideally mixed between bottom and ktop.
  iblmin = klev-2

!-- half level pressure values, assuming 101320. Pa surface pressure

  DO jk=1,klevp1
    zph(jk)=vct(jk)+vct(jk+nvclev)*101320.0_wp
  END DO
!
! -- full level pressure
!
  DO jk = 1, klev
    zp(jk)=(zph(jk)+zph(jk+1))*0.5_wp
  END DO
!
  DO jk = 1, klev
    zh(jk)=(zph(klevp1)-zp(jk))/(grav*1.25_wp)
  END DO
!
! -- search for highest level below 2000m
!
  DO jk = 1, klev
    iblmax=jk
    IF(zh(jk).LT.2000.0_wp) EXIT
  END DO
      
!
!     ------------------------------------------------------------------
!
!*         2.     NEW THERMODYNAMIC VARIABLE AND BOUNDARY CONDITIONS.
!
!*         2.1     REPLACE T BY CP(Q)*T+GZ IN THE ATMOSPHERE.
!
#ifdef _PROFILE
  CALL trace_start ('vdiff_loop_1', 21)
#endif

  DO 212 jk=ktdia,klev
     CALL prepare_ua_index_spline('vdiff (1)',kproma,ptm1(1,jk),idx(1),za(1))
     CALL lookup_ua_spline(kproma,idx(1),za(1),zua(1))

     zpapm1i(1:kproma) = 1._wp/papm1(1:kproma,jk)
     zteta1(1:kproma,jk) = (100000._wp*zpapm1i(1:kproma))**zkappa

!IBM* NOVECTOR
     DO 211 jl=1,kproma
        zteta1(jl,jk)=ptm1(jl,jk)*zteta1(jl,jk)
        zx(jl,jk)=pxlm1(jl,jk)+pxim1(jl,jk)   ! total cloud water
        zcptgz(jl,jk)=pgeom1(jl,jk)+ptm1(jl,jk)                        &
                                   *cpd*(1._wp+vtmpc2*pqm1(jl,jk))
        ztvir1(jl,jk)=zteta1(jl,jk)*(1._wp+vtmpc1*pqm1(jl,jk)-zx(jl,jk))
!        lo=ptm1(jl,jk).GE.tmelt
!        zfaxe(jl,jk)=MERGE(alv,als,lo)
        zfaxe(jl,jk)=FSEL(ptm1(jl,jk)-tmelt,alv,als)
        zbet=zfaxe(jl,jk)/zcpd
        zusus1=zbet*zteta1(jl,jk)/ptm1(jl,jk)*zx(jl,jk)
        zlteta1(jl,jk)=zteta1(jl,jk)-zusus1
        zes=zua(jl)*zpapm1i(jl)
        zes=MIN(zes,0.5_wp)
        zqss(jl,jk)=zes/(1._wp-vtmpc1*zes)
211   END DO
212 END DO

  DO 214 jk=ktdia,klevm1
     DO 213 jl=1,kproma
        zhh(jl,jk)=(pgeom1(jl,jk)-pgeom1(jl,jk+1))
        zsdep1=(paphm1(jl,jk)-paphm1(jl,jk+1))                         &
              /(paphm1(jl,jk)-paphm1(jl,jk+2))
        zsdep2=(paphm1(jl,jk+1)-paphm1(jl,jk+2))                       &
              /(paphm1(jl,jk)  -paphm1(jl,jk+2))
        zqssm(jl,jk)=zsdep1*zqss(jl,jk)+zsdep2*zqss(jl,jk+1)
        ztmitte(jl,jk)=zsdep1*ptm1(jl,jk)+zsdep2*ptm1(jl,jk+1)
        ztvirmit(jl,jk)=zsdep1*ztvir1(jl,jk)+zsdep2*ztvir1(jl,jk+1)
        zfaxen(jl,jk)=zsdep1*zfaxe(jl,jk)+zsdep2*zfaxe(jl,jk+1)
        zlwcmit(jl,jk)=zsdep1*zx(jl,jk)+zsdep2*zx(jl,jk+1)
        zqmit(jl,jk)=zsdep1*pqm1(jl,jk)+zsdep2*pqm1(jl,jk+1)
        ztemit(jl,jk)=zsdep1*zteta1(jl,jk)+zsdep2*zteta1(jl,jk+1)
        zccover(jl,jk)=paclc(jl,jk)*zsdep1+paclc(jl,jk+1)*zsdep2
213  END DO
214 END DO

#ifdef _PROFILE
  CALL trace_stop ('vdiff_loop_1', 21)
  CALL trace_start ('vdiff_loop_2', 22)
#endif

!

  IF (lvdiff) THEN

!
!!------------------- MO_SURFACE Coupling-----------------
! INITIALISE udif and vdif which previously came from momentum
! transfer und is then used for wind stress calculations
! wind stress is now calculated before momentum transfer
!- kalle 010904
!-----------------------------------------------------------
     IF (lstart) THEN

        DO jl=1,kproma
           udif(jl,jrow)=ztpfac2*pum1(jl,klev)
           vdif(jl,jrow)=ztpfac2*pvm1(jl,klev)
           ustarm(jl,jrow)=1._wp
        END DO

     END IF
! JSBACH-end udif, vdif--------------------------------------

!
! Compute planetary boundary layer extension
!
! JSBACH note: in standard ECHAM5, ustarm is computed before zhdyn for
!     the current time step. Here, ustarm comes from the call to
!     mo_surface at the previous timestep.
!     But this should have only a minor effect on ihpbl and ghabl.
!
    DO jl = 1,kproma
       zcor=MAX(ABS(coriol_2d(jl,jrow)),zepcor)
       zhdyn(jl)=MIN(pgeom1(jl,1)/grav,zchneu*ustarm(jl,jrow)/zcor)
       ihpblc(jl)=klev
       ihpbld(jl)=klev
    END DO

     DO jk=klevm1,1,-1
        DO jl=1,kproma
           zds=zcptgz(jl,jk)-zcptgz(jl,klev)
           zdz=pgeom1(jl,jk)/grav-zhdyn(jl)
           ihpblc(jl)=MERGE(jk,ihpblc(jl),                             &
                            ihpblc(jl).EQ.klev.AND.zds.GT.0._wp)
           ihpbld(jl)=MERGE(jk,ihpbld(jl),                             &
                            ihpbld(jl).EQ.klev.AND.zdz.GE.0._wp)
     END DO
  END DO
!
     DO  jl=1,kproma
       IF (.NOT. lrce) THEN
        ihpbl(jl)=MIN(ihpblc(jl),ihpbld(jl))
       ELSE
        ihpbl(jl)=ihpblc(jl)         ! pbl extension in RCE only given by convective pbl height
       END IF
        zghabl(jl)=MIN(50000._wp,pgeom1(jl,ihpbl(jl)))
!--- Included PBL top level export (Jan Kazil 10/2008)------------------
        ppbl(jl) = REAL(ihpbl(jl),wp)
!--- End included ------------------------------------------------------
     END DO

!
!     ==================================================================
!
!*       3.5   Vertical loop: Computation of basic quantities:
!              wind shear, buoyancy, Ri-number, mixing length
!
#ifdef _PROFILE
     CALL trace_start ('vdiff_loop_6', 26)
#endif

     DO 372 jk=ktdia,klevm1
        DO 361 jl=1,kproma
           zqtmit=zlwcmit(jl,jk)+zqmit(jl,jk)
           zfux=zfaxen(jl,jk)/(zcpd*ztmitte(jl,jk))
           zfox=zfaxen(jl,jk)/(zrd*ztmitte(jl,jk))
           zmult1=1._wp+vtmpc1*zqtmit
           zmult2=zfux*zmult1-zrvrd
           zmult3=zrdrv*zfox*zqssm(jl,jk)                              &
                  /(1._wp+zrdrv*zfux*zfox*zqssm(jl,jk))
           zmult5=zmult1-zmult2*zmult3
           zmult4=zfux*zmult5-1._wp
           zdus1=zccover(jl,jk)*zmult5+(1._wp-zccover(jl,jk))*zmult1
           zdus2=zccover(jl,jk)*zmult4+(1._wp-zccover(jl,jk))*vtmpc1
           zteldif=(zlteta1(jl,jk)-zlteta1(jl,jk+1))/zhh(jl,jk)*grav
           zthvirdif=(ztvir1(jl,jk)-ztvir1(jl,jk+1))/zhh(jl,jk)*grav
           zdqtot=(pqm1(jl,jk)+zx(jl,jk))-(pqm1(jl,jk+1)+zx(jl,jk+1))
           zqddif=zdqtot/zhh(jl,jk)*grav
           zbuoy=(zteldif*zdus1+ztemit(jl,jk)*zdus2*zqddif)            &
                 *grav/ztvirmit(jl,jk)
           zdivv=(pum1(jl,jk)-pum1(jl,jk+1))**2
           zdivv1=(pvm1(jl,jk)-pvm1(jl,jk+1))**2
           zshear=(zdivv+zdivv1)*(grav/zhh(jl,jk))**2
           zri=zbuoy/MAX(zshear,zepshr)
!
!      ASYMPTOTIC MIXING LENGTH FOR MOMENTUM AND
!      HEAT (ZLAM) ABOVE THE PBL AS A FUNCTION OF HEIGHT
!      ACCORDING TO HOLTSLAG AND BOVILLE (1992), J. CLIMATE.
!
           zhexp=EXP(1._wp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
           zlam=zzzlam+(zcons3-zzzlam)*zhexp
           IF(jk.GE.ihpbl(jl)) THEN
              zcons23=zcons25
           ELSE
              zcons23=zcons2/zlam
           END IF
!
!     MIXING LENGTH (BLACKADAR) + STABILITY DEPENDENT FUNCTION
!
           z2geomf=pgeom1(jl,jk)+pgeom1(jl,jk+1)
           zz2geo=zcons2*z2geomf
           zmix=zz2geo/(1._wp+zcons23*z2geomf)
!
!      STABILITY FUNCTIONS (LOUIS, 1979)
!
           IF(zri.LT.0._wp) THEN
              zalh2=zmix*zmix
              zucf=1._wp/                                              &
                   (1._wp+zcons5*zalh2*SQRT(ABS(zri)*(((pgeom1(jl,jk)  &
                       /pgeom1(jl,jk+1))**zcons6-1._wp)/(pgeom1(jl,jk) &
                       -pgeom1(jl,jk+1)))**3/(pgeom1(jl,jk+1))))
              zsh=zshn*(1._wp-zcons9*zri*zucf)*zmix
              zsm=zsmn*(1._wp-zcons8*zri*zucf)*zmix
           ELSE
              zsh=zshn/(1._wp+zcons8*zri*SQRT(1._wp+zri))*zmix
              zsm=zsmn/(1._wp+zcons8*zri/SQRT(1._wp+zri))*zmix
           END IF
!
!       Dimensionless coefficients multiplied by pressure
!            thicknesses for momentum and heat exchange
!
           zzb=zshear*zsm-zbuoy*zsh
           zdisl=zda1*zmix/ztmst
           zktest=1._wp+(zzb*ztmst+SQRT(ptkem1(jl,jk))*2._wp)/zdisl
           IF (zktest.LE.1._wp) THEN
              ztkevn(jl,jk)=ztkemin
           ELSE
              ztkevn(jl,jk)=MAX(ztkemin,(zdisl*(SQRT(zktest)-1._wp))**2)
           END IF
           IF(lstart) THEN
              ptkem1(jl,jk)=ztkevn(jl,jk)
              ptkem(jl,jk)=ztkevn(jl,jk)
           END IF
           ztkesq=SQRT(MAX(ztkemin,ptkem1(jl,jk)))
           zztvm=(ptvm1(jl,jk)+ptvm1(jl,jk+1))*0.5_wp
           zalf=paphm1(jl,jk+1)/(zztvm*zhh(jl,jk)*zrd)
           zcfm(jl,jk)=zsm*ztkesq*zcons18*zalf
           zcfh(jl,jk)=zsh*ztkesq*zcons18*zalf
           zcfv(jl,jk)=0.5_wp*zcfh(jl,jk)
           zcdum(jl,jk)=zcfm(jl,jk)/ztkesq*SQRT(ztkevn(jl,jk))
           zcthv(jl,jk)=zcfh(jl,jk)/ztkesq*SQRT(ztkevn(jl,jk))
           zthvprod=2._wp*zsh*ztkesq*zthvirdif**2
           zthvdiss=pthvvar(jl,jk)*ztkesq/(zda1*zmix)
           zthvvar(jl,jk)=pthvvar(jl,jk)+(zthvprod-zthvdiss)*ztmst
           zthvvar(jl,jk)=MAX(ztkemin,zthvvar(jl,jk))
361     END DO
372  END DO
!
     DO jl=1,kproma
        zthvvar(jl,klev) = zthvvar(jl,klevm1)
     END DO
!
#ifdef _PROFILE
  CALL trace_stop ('vdiff_loop_6', 26)
  CALL trace_start ('vdiff_loop_7', 27)
#endif
!
!     ------------------------------------------------------------------
!
!*         5.     DIFFUSION IMPLICIT COMPUTATIONS FOR HEAT (S.+L.).
!
     DO 502 jk=1,klev
        DO 501 jl=1,kproma
           ztdif(jl,jk)=0._wp
           zqdif(jl,jk)=0._wp
           zxldif(jl,jk)=0._wp
           zxidif(jl,jk)=0._wp
501     END DO
502  END DO
!
!*         5.1     SETTING OF RIGHT HAND SIDES.
!
     DO 512 jk=itop,klev
        DO 511 jl=1,kproma
           ztdif(jl,jk)=ztpfac2*zcptgz(jl,jk)
           zqdif(jl,jk)=ztpfac2*pqm1(jl,jk)
           zxldif(jl,jk)=ztpfac2*pxlm1(jl,jk)
           zxidif(jl,jk)=ztpfac2*pxim1(jl,jk)
511     END DO
512  END DO
!
!*         5.2     TOP LAYER ELIMINATION.
!
     DO 521 jl=1,kproma
        zqdp=1._wp/(paphm1(jl,itopp1)-paphm1(jl,itop))
        zdisc=1._wp/(1._wp+zcfh(jl,itop)*zqdp)
        zebsh(jl,itop)=zdisc*(zcfh(jl,itop)*zqdp)
        zdisv=1._wp/(1._wp+zcfv(jl,itop)*zqdp)
        zebsv(jl,itop)=zdisv*(zcfv(jl,itop)*zqdp)
        ztdif(jl,itop)=zdisc*ztdif(jl,itop)
        zqdif(jl,itop)=zdisc*zqdif(jl,itop)
        zxldif(jl,itop)=zdisc*zxldif(jl,itop)
        zxidif(jl,itop)=zdisc*zxidif(jl,itop)
521  END DO
!
!*         5.3     ELIMINATION FOR MIDDLE LAYERS.
!
     DO 532 jk=itopp1,klevm1
        DO 531 jl=1,kproma
           zqdp=1._wp/(paphm1(jl,jk+1)-paphm1(jl,jk))
           zfac=zcfh(jl,jk-1)*zqdp
           zdisc=1._wp/(1._wp+zfac*(1._wp-zebsh(jl,jk-1))              &
                                                    +zcfh(jl,jk)*zqdp)
           zebsh(jl,jk)=zdisc*(zcfh(jl,jk)*zqdp)
           zfav=zcfv(jl,jk-1)*zqdp
           zdisv=1._wp/(1._wp+zfav*(1._wp-zebsv(jl,jk-1))                 &
                                                    +zcfv(jl,jk)*zqdp)
           zebsv(jl,jk)=zdisv*(zcfv(jl,jk)*zqdp)
           ztdif(jl,jk)=zdisc*(ztdif(jl,jk)+zfac*ztdif(jl,jk-1))
           zqdif(jl,jk)=zdisc*(zqdif(jl,jk)+zfac*zqdif(jl,jk-1))
           zxldif(jl,jk)=zdisc*(zxldif(jl,jk)+zfac*zxldif(jl,jk-1))
           zxidif(jl,jk)=zdisc*(zxidif(jl,jk)+zfac*zxidif(jl,jk-1))
531     END DO
532  END DO
!
!*         5.4     BOTTOM LAYER ELIMINATION.
!
     DO 541 jl=1,kproma
        zqdp=1._wp/(paphm1(jl,klevp1)-paphm1(jl,klev))
        zfac=zcfh(jl,klevm1)*zqdp
        zdisx=1._wp/(1._wp+zfac*(1._wp-zebsh(jl,klevm1)))
        zfav=zcfv(jl,klevm1)*zqdp
        zdisv=1._wp/(1._wp+zfav*(1._wp-zebsv(jl,klevm1)))
        zxldif(jl,klev)=zdisx*(zxldif(jl,klev)+zfac*zxldif(jl,klevm1))
        zxidif(jl,klev)=zdisx*(zxidif(jl,klev)+zfac*zxidif(jl,klevm1))
541  END DO

     IF (ico2idx > 0) THEN
        IF (trlist%ti(ico2idx)%nvdiff /= 1) &
             CALL finish('vdiff','Need %nvdiff = 1 for CO2 tracer')
        zco2(1:kproma) = pxtm1(1:kproma,klev,ico2idx)
     ELSE
        zco2(1:kproma) = co2m1(1:kproma,klev,jrow)
     END IF

     CALL update_surface(                 & !! logistic parameters
       kproma, jrow,                      & !! length of vector
       klev, klevp1, klevm1,              & !! lowest level number, lln plus 1, lln minus 1
       lpnorth(1:kproma),                 & !! true in northern hemisphere
!!---INPUT------------------!! Atmospheric conditions lowest level INPUT
       pxlm1(1:kproma,klev),              & !! liquid clouds (water)
       pxim1(1:kproma,klev),              & !! frozen clouds (ice)
       pgeom1(1:kproma,klev),             & !! geopotential above surface
       ptm1(1:kproma,klev),               & !! temperature (t-dt)
       pqm1(1:kproma,klev),               & !! humidity (t-dt)
       papm1(1:kproma,klev),              & !! full level pressure (means middle of the layer)
       paphm1(1:kproma,1:klevp1),         & !! half level pressures (bottom of the layer)
       pum1(1:kproma,klev),               & !! wind_u
       pvm1(1:kproma,klev),               & !! wind_v
       prain(1:kproma),                   & !! rain (convective and large scale) over last time step
       psnow(1:kproma),                   & !! snow (convective and large scale) over last time step
       pemter(1:kproma,klevp1),           & !! longwave net
       jsswvis(1:kproma),                 & !! net surface visible
       jsswdifvis(1:kproma),              & !! fraction of diffuse visible
       jsswnir(1:kproma),                 & !! net surface near infrared
       jsswdifnir(1:kproma),              & !! fraction of diffuse near infrared
       jsswpar(1:kproma),                 & !! downward surface PAR
       jsswdifpar(1:kproma),              & !! fraction of diffuse PAR
       amu0_x(1:kproma,jrow),             & !! solar zenith angle
       paclc(1:kproma,klev),              & !! cloud cover
       ptsw(1:kproma),                    & !! ocean_temp read from sst field or ocean model
       pocu(1:kproma),                    & !! ocean_u_velocity
       pocv(1:kproma),                    & !! ocean_v_velocity
       psiced(1:kproma),                  & !! seaice depth from clsst
       pameltdepth(1:kproma),             & !! melt pond depth
       pameltfrac(1:kproma),              & !! melt pond fraction
       pamlcorr(1:kproma),                & !! mixed layer ocean
       pamlcorac(1:kproma),               & !! mixed layer ocean
       pamlheatac(1:kproma),              & !! mixed layer ocean
       zcfh(1:kproma,1:klev),             & !! turb diffusion coeff for RM-scheme
       zebsh(1:kproma,1:klev),            & !! coeff from diffusion scheme for RM
       zqdif(1:kproma,1:klev),            & !! humid. coeff from diffusion scheme
       ztdif(1:kproma,1:klev),            & !! temp. coeff from diffusion scheme
       udif(1:kproma,jrow),              & !! wind speed coeff from momentum diffusion
       vdif(1:kproma,jrow),              & !! wind speed coeff ..for wind stress calc.
       zghabl(1:kproma),                  & !! Geopotential of PBL extension
       pi0(1:kproma),                     & !! solar incidence factor from radiation scheme
       ptrsol(1:kproma,klevp1),           & !! solar radiation
       zco2(1:kproma),                    & !! CO2 concentration in lowest level
!! - OUTPUT FROM SURFACE SCHEME ----------------------------------------
       palbedo(1:kproma),       palbedo_vis(1:kproma),                 &
       palbedo_nir(1:kproma),                                          &
       palbedo_vis_dir(1:kproma),     palbedo_nir_dir(1:kproma),       &
       palbedo_vis_dif(1:kproma),     palbedo_nir_dif(1:kproma),       &
       ptrflw(1:kproma),        ptrfli(1:kproma),                      &
       psofll(1:kproma),        psoflw(1:kproma),                      &
       psofli(1:kproma) ,       ptrfllac(1:kproma),                    &
       ptrflwac(1:kproma),      ptrfliac(1:kproma),                    &
       psofllac(1:kproma),      psoflwac(1:kproma),                    &
       psofliac(1:kproma),      palsol(1:kproma),                      &
       palsoi(1:kproma),        palsow(1:kproma),                      &
       palsobs(1:kproma),       palsom(1:kproma),                      &
       ptaus(1:kproma),         pfage(1:kproma),                       &
       psnifrac(1:kproma),      pbarefrac(1:kproma),                   &
       ustarm(1:kproma,jrow),   cdum(1:kproma,jrow), cduh(1:kproma,jrow), &
       tkevn(1:kproma,jrow),    ztdif_new(1:kproma),                   &
       zqdif_new(1:kproma),     ztvh_new(1:kproma),                    &
       zqsurf_new(1:kproma),    zth_new(1:kproma),                     &
       pwind10w(1:kproma),      pu10(1:kproma),                        &
       pv10(1:kproma),          pwimax(1:kproma),                      &
       pwind10(1:kproma),       pdew2(1:kproma),                       &
       ptemp2(1:kproma),        pt2max(1:kproma),                      &
       pt2min(1:kproma),        pevaplac(1:kproma),                    &
       pevapwac(1:kproma),      pevapiac(1:kproma),                    &
       pevap(1:kproma),         pahfllac(1:kproma),                    &
       pahflwac(1:kproma),      pahfliac(1:kproma),                    &
       pahfl(1:kproma),         pahfslac(1:kproma),                    &
       pahfswac(1:kproma),      pahfsiac(1:kproma),                    &
       pahfs(1:kproma),         pqhfla(1:kproma),                      &
       pevapw(1:kproma),        pevapi(1:kproma),  pahfsl(1:kproma),   &
       pahfsw(1:kproma),        pahfsi(1:kproma),                      &
       pevapot(1:kproma),       pahflw(1:kproma),                      &
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-03  for sat simulator               ----------
       pahfli(1:kproma),        pahfll(1:kproma),                      &
!---------- kai zhang 2009-07-03  for sat simulator               ----------
!---------------------------------------------------------------------------
       psni(1:kproma),          pahfice(1:kproma),                     &  !Note: psni is INOUT (comes already from ocean)
       pfluxres(1:kproma),      pqres(1:kproma),                       &
       pahfcon(1:kproma),       pahfres(1:kproma),                     &
       ptsi(1:kproma),          ptslnew(1:kproma),                     &
       pzti(1:kproma),          pzteffl4(1:kproma),                    &
       pztsnew(1:kproma),       ptsurf(1:kproma),                      &
       paz0w(1:kproma),         paz0i(1:kproma),                       &
       paz0l(1:kproma),         paz0(1:kproma),                        &
       pustrl(1:kproma),        pvstrl(1:kproma),                      &
       pustrw(1:kproma),        pvstrw(1:kproma),                      &
       pustri(1:kproma),        pvstri(1:kproma),                      &
       pustr(1:kproma),         pvstr(1:kproma),                       &
       ztte_corr(1:kproma),     ptsw_new(1:kproma),                    &
       pseaice(1:kproma),       pseaice_new(1:kproma),                 &
       psiced_new(1:kproma),    pradtemp_old(1:kproma),                &
#ifdef HAMMOZ /* SF */
!SF gf #78
       zcfml(1:kproma), zcfmw(1:kproma), zcfmi(1:kproma),              &
       zcfnc(1:kproma),                                                &
       zcfncl(1:kproma), zcfncw(1:kproma), zcfnci(1:kproma),           &
       zrib(1:kproma),                                                 &
       zril(1:kproma), zriw(1:kproma), zrii(1:kproma),                 &
       ztvl(1:kproma), ztvw(1:kproma), ztvi(1:kproma),                 &
       pfrl(1:kproma), pfrw(1:kproma), pfri(1:kproma),                 &
       palake(1:kproma),                                               &
       zcvs(1:kproma), zsrfll(1:kproma),                               &
       zcdn(1:kproma),                                                 &
       zcdnl(1:kproma), zcdnw(1:kproma), zcdni(1:kproma),              &
       zvelo10m(1:kproma),                                             & 
#else
       zcfnc(1:kproma), zrib(1:kproma), ztvl(1:kproma),                & 
       pfrl(1:kproma), pfrw(1:kproma), pfri(1:kproma),                 & 
       palake(1:kproma),                                               &
       zcvs(1:kproma), zsrfll(1:kproma), zcdnl(1:kproma),              & 
       zvelo10m(1:kproma),                                             & 
#endif
       co2_flux_ocean(1:kproma,jrow), co2_flux_land(1:kproma,jrow),         &
       co2_flux(1:kproma,jrow),                                             &
       co2_flux_npp(1:kproma,jrow), co2_flux_soilresp(1:kproma,jrow),       &
       co2_flux_herbivory(1:kproma,jrow), co2_flux_dynveg(1:kproma,jrow),   & 
       co2_emission_lcc(1:kproma,jrow), co2_emission_harvest(1:kproma,jrow))

  !initialize surface emission for tracers 
!!++mgs  ### WORKAROUND until ztslm1 is avaialble from update_surface as ztslp1
  ztslm1(1:kproma) = ptm1(1:kproma, klev)
!!--mgs END OF WORKAROUND

     IF (lanysubmodel) THEN
        vphysc%velo10m(1:kproma,jrow) = zvelo10m(1:kproma)
        vphysc%tsw(1:kproma,jrow) = ptsw(1:kproma)
     END IF

     IF(ktrac.GT.0) THEN
        DO jt=1,ktrac
           DO 3230 jl=1,kproma
              zxtems(jl,jt)=0._wp
3230       END DO
        END DO
     END IF !(ktrac.GT.0)

!>>SF gf #1 Some of these updates must be done before the call to vdiff_subm
!--------------------
! Hand it over ! begin
!--------------------

     DO jl=1,kproma
        zcdum(jl,klev)  = cdum(jl,jrow)
        zcfm(jl,klev)   = cdum(jl,jrow)
        zcthv(jl,klev)  = cduh(jl,jrow)
        ztkevn(jl,klev) = tkevn(jl,jrow)
        ztdif(jl,klev)  = ztdif_new(jl)
        zqdif(jl,klev)  = zqdif_new(jl)
        ptsw(jl)        = ptsw_new(jl)
        pseaice(jl)     = pseaice_new(jl)
        psiced(jl)      = psiced_new(jl)
     END DO
     IF(lstart) THEN
        DO jl=1,kproma
           ptkem1(jl,klev)=ztkevn(jl,klev)
           ptkem(jl,klev)=ztkevn(jl,klev)
        END DO
     END IF
!--------------------
! Hand it over ! end
!--------------------

!<<SF gf #1

     IF (lanysubmodel) THEN

        CALL vdiff_subm      (kproma, kbdim,  klev,   klevp1,        &
                              ktrac,  jrow,                          &
                              ptm1,   pum1,   pvm1,   pqm1,          &
                              papm1,  paphm1, paphp1, pgeom1, ztslm1,& 
                              pxtm1,  pseaice,pforest,               & 
                              pfrl,   pfrw,   pfri,   zcvs,   pcvw,  &
                              pvgrat, ptsw,   ptsi,                  &
                              pu10,   pv10,                          &
                              paz0,   paz0l,  paz0w,  paz0i,         & 
#ifdef HAMMOZ /* SF */
!SF gf #78
                              zcfml,  zcfmw,  zcfmi,                 &
                              zcfncl, zcfncw, zcfnci, zepdu2, zkap,  &
                              zril,   zriw,   zrii,   ztvir1,        &
                              ztvl,   ztvw,   ztvi,                  &
                              zsrfll, zcdnl,  zcdnw,  zcdni,         &
                              zqss,   pvlt,                          &
#else
                              zcfm,   zcfnc, zepdu2, zkap,          &
                              zrib,   ztvir1, ztvl,                  &
                              zsrfll, zcdnl,  zqss,   pvlt,          & 
#endif
                              lpland,                                & 
                              pxtte,  zxtems,                        & 
                              pxlm1,  pxim1                          &
#ifdef HAMMOZ /* SF */
!SF gf #161
                             ,ihpbl                                  &
#endif
                             )                                        
!--mgs
     END IF

     DO 505 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 504 jk=1,klev
           DO 503 jl=1,kproma
              zxtdif(jl,jk,jt)=0._wp
503        END DO
504     END DO

!!$        DO jk=ibl,klev
!!$           DO jl=1,kproma
!!$              zqdp = 1._wp/(paphm1(jl,klevp1)-paphm1(jl,ibl))
!!$              zxtdif(jl,jk,jt) = ztpfac2 * pxtm1(jl,jk,jt)             &
!!$                       + ztmst * g * zqdp * zxtems(jl,jt)
!!$           END DO
!!$        END DO
! original code:
        DO 516 jk=itop,klev
           DO 514 jl=1,kproma
              zxtdif(jl,jk,jt)=ztpfac2*pxtm1(jl,jk,jt)
514        END DO
516     END DO

        DO 526 jl=1,kproma
           zqdp=1._wp/(paphm1(jl,itopp1)-paphm1(jl,itop))
           zdisc=1._wp/(1._wp+zcfh(jl,itop)*zqdp)
           zxtdif(jl,itop,jt)=zdisc*zxtdif(jl,itop,jt)
526     END DO

        DO 536 jk=itopp1,klevm1
           DO 534 jl=1,kproma
              zqdp=1._wp/(paphm1(jl,jk+1)-paphm1(jl,jk))
              zfac=zcfh(jl,jk-1)*zqdp
              zdisc=1._wp/(1._wp+zfac*(1._wp-zebsh(jl,jk-1))           &
                                                    +zcfh(jl,jk)*zqdp)
              zxtdif(jl,jk,jt)=zdisc *                                 &
                          (zxtdif(jl,jk,jt) + zfac*zxtdif(jl,jk-1,jt))
534        END DO
536     END DO

        DO 543 jl=1,kproma
           zqdp=1._wp/(paphm1(jl,klevp1)-paphm1(jl,klev))
           zfac=zcfh(jl,klevm1)*zqdp
           zdisxt=1._wp/(1._wp+zfac*(1._wp-zebsh(jl,klevm1)))
           zxtdif(jl,klev,jt)=zdisxt * (zxtdif(jl,klev,jt)             &
!! the following line is NOT commented in the original code, so it's original code here!!!
                                 + ztmst * grav * zqdp * zxtems(jl,jt) &
                                    + zfac*zxtdif(jl,klevm1,jt)        &
                                       )
543     END DO

505  END DO


!
!*         5.5     BACK-SUBSTITUTION.
!
     DO 552 jk=klevm1,itop,-1
        DO 551 jl=1,kproma
           ztdif(jl,jk)=ztdif(jl,jk)+zebsh(jl,jk)*ztdif(jl,jk+1)
           zqdif(jl,jk)=zqdif(jl,jk)+zebsh(jl,jk)*zqdif(jl,jk+1)
           zxldif(jl,jk)=zxldif(jl,jk)+zebsh(jl,jk)*zxldif(jl,jk+1)
           zxidif(jl,jk)=zxidif(jl,jk)+zebsh(jl,jk)*zxidif(jl,jk+1)
551     END DO
552  END DO
!
     DO 558 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 556 jk=klevm1,itop,-1
           DO 554 jl=1,kproma
              zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt)+                       &
                               zebsh(jl,jk)*zxtdif(jl,jk+1,jt)
554        END DO
556     END DO
558  END DO
!     ==================================================================
!
!*       3.8  DIFFUSION IMPLICIT COMPUTATIONS FOR TKE
!                       AND VARIANCE OF VIRTUAL POTENTIAL TEMPEATURE
!
     DO 380 jk=ktdia,klev
        DO 381 jl=1,kproma
           zedif(jl,jk)=ztpfac2*ztkevn(jl,jk)
           zthvdif(jl,jk)=ztpfac2*zthvvar(jl,jk)
381     END DO
380  END DO
!
     DO 385 jl=1,kproma
        ztcoe(jl)=(zcdum(jl,itop)+zcdum(jl,itopp1))*0.5_wp
        ztthv(jl)=(zcthv(jl,itop)+zcthv(jl,itopp1))*0.5_wp
        zqdp=1._wp/(papm1(jl,itopp1)-papm1(jl,itop))
        zdisc=1._wp/(1._wp+(zcdum(jl,itop)+zcdum(jl,itopp1))           &
                                                         *0.5_wp*zqdp)
        zdisthv=1._wp/(1._wp+(zcthv(jl,itop)+zcthv(jl,itopp1))         &
                                                         *0.5_wp*zqdp)
        zebsm(jl,itop)=zdisc*(zcdum(jl,itop)+zcdum(jl,itopp1))         &
                                                         *0.5_wp*zqdp
        zebthv(jl,itop)=zdisthv*(zcthv(jl,itop)+zcthv(jl,itopp1))      &
                                                         *0.5_wp*zqdp
        zedif(jl,itop)=zdisc*zedif(jl,itop)
        zthvdif(jl,itop)=zdisthv*zthvdif(jl,itop)
385  END DO
!
     DO 386 jk=itopp1,klev-2
        DO 387 jl=1,kproma
           zqdp=1._wp/(papm1(jl,jk+1)-papm1(jl,jk))
           zfac=ztcoe(jl)*zqdp
           zfacthv=ztthv(jl)*zqdp
           ztcoe(jl)=(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5_wp
           ztthv(jl)=(zcthv(jl,jk+1)+zcthv(jl,jk))*0.5_wp
           zdisc=1._wp/(1._wp+zfac*(1._wp-zebsm(jl,jk-1))                 &
                     +(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5_wp*zqdp)
           zdisthv=1._wp/(1._wp+zfacthv*(1._wp-zebthv(jl,jk-1))           &
                     +(zcthv(jl,jk+1)+zcthv(jl,jk))*0.5_wp*zqdp)
           zebsm(jl,jk)=zdisc*(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5_wp*zqdp
           zebthv(jl,jk)=zdisthv*(zcthv(jl,jk+1)+zcthv(jl,jk))*0.5_wp*zqdp
           zedif(jl,jk)=zdisc*(zedif(jl,jk)+zfac*zedif(jl,jk-1))
           zthvdif(jl,jk)=zdisthv*(zthvdif(jl,jk)+zfacthv*zthvdif(jl,jk-1))
387     END DO
386  END DO
!
     DO 390 jl=1,kproma
        zqdp=1._wp/(papm1(jl,klev)-papm1(jl,klevm1))
        zfac=ztcoe(jl)*zqdp
        zfacthv=ztthv(jl)*zqdp
        ztcoe(jl)=(zcdum(jl,klev)+zcdum(jl,klevm1))*0.5_wp
        ztthv(jl)=(zcthv(jl,klev)+zcthv(jl,klevm1))*0.5_wp
        zdisc=1._wp/(1._wp+zfac*(1._wp-zebsm(jl,klev-2))               &
                  +(zcdum(jl,klev)+zcdum(jl,klevm1))*0.5_wp*zqdp)
        zdisthv=1._wp/(1._wp+zfacthv*(1._wp-zebthv(jl,klev-2))         &
                  +(zcthv(jl,klev)+zcthv(jl,klevm1))*0.5_wp*zqdp)
        zedif(jl,klevm1)=zdisc*((zcdum(jl,klev)+zcdum(jl,klevm1))      &
                    *0.5_wp*zqdp*zedif(jl,klev)+zedif(jl,klevm1)       &
                    +zfac*zedif(jl,klev-2))
        zthvdif(jl,klevm1)=zdisthv*((zcthv(jl,klev)+zcthv(jl,klevm1))  &
                    *0.5_wp*zqdp*zthvdif(jl,klev)+zthvdif(jl,klevm1)   &
                    +zfacthv*zthvdif(jl,klev-2))
        zthvdif(jl,klev) = zthvdif(jl,klevm1)
390  END DO
!
     DO 392 jk=klev-2,itop,-1
        DO 393 jl=1,kproma
           zedif(jl,jk)=zedif(jl,jk)+zebsm(jl,jk)*zedif(jl,jk+1)
           zthvdif(jl,jk)=zthvdif(jl,jk)+zebthv(jl,jk)*zthvdif(jl,jk+1)
393     END DO
392  END DO
!
!*    TIME INTEGRATION OF TURBULENT KINETIC ENERGY AND CHECK
!     NEW VARIANCE OF VIRTUAL POTENTIAL TEMPERATURE
!
     DO 394 jk=itop,klev
        ztest=0._wp
        DO 395 jl=1,kproma
           ptke(jl,jk)=zedif(jl,jk)+ztpfac3*ztkevn(jl,jk)
           ztest=ztest+MERGE(1._wp,0._wp,ptke(jl,jk)<0._wp)
           pthvvar(jl,jk)=zthvdif(jl,jk)+ztpfac3*zthvvar(jl,jk)
           pthvvar(jl,jk)=MAX(ztkemin,pthvvar(jl,jk))
395     END DO
        IF(ztest.NE.0._wp) CALL finish('vdiff','TKE IS NEGATIVE')
394  END DO
!
!*    TIME FILTER FOR TURBULENT KINETIC ENERGY
!
     IF(.NOT.lstart) THEN
       zeps=eps
     ELSE
       zeps=0._wp
     END IF

     IF (lnwp) zeps=0._wp

     DO 397 jk=ktdia,klev
       DO 396 jl=1,kproma
         ptkem1(jl,jk)=ptkem(jl,jk)                                    &
                   +zeps*(ptkem1(jl,jk)-2._wp*ptkem(jl,jk)+ptke(jl,jk))
         ptkem(jl,jk)=ptke(jl,jk)
396     END DO
397  END DO
!
!      STD DEV OF VIRTUAL POT TEMPERATURE AT STANDARD HALF LEVEL KLEV
!      (CORRESPONDING TO HALF LEVEL KLEV-1 FOR TKE AND PTHVVAR)
!
       DO jl=1,kproma
          pthvsig(jl) = SQRT(pthvvar(jl,klev-1))
       END DO
!
!     ------------------------------------------------------------------
!
!*         4.     DIFFUSION IMPLICIT COMPUTATIONS FOR MOMENTUM.
!
!*         4.1     SETTING OF RIGHT HAND SIDES.
!
     DO 412 jk=itop,klevm1
        DO 411 jl=1,kproma
           zudif(jl,jk)=ztpfac2*pum1(jl,jk)
           zvdif(jl,jk)=ztpfac2*pvm1(jl,jk)
411     END DO
412  END DO
         
        DO 413 jl=1,kproma
           zudif(jl,klev)=ztpfac2*(pum1(jl,klev)-                      &
                         pocu(jl)*(1._wp-pfrl(jl)))
           zvdif(jl,klev)=ztpfac2*(pvm1(jl,klev)-                      &
                         pocv(jl)*(1._wp-pfrl(jl)))
413     ENDDO
!
!*         4.2     TOP LAYER ELIMINATION.
!
     DO 421 jl=1,kproma
        zqdp=1._wp/(paphm1(jl,itopp1)-paphm1(jl,itop))
        zdisc=1._wp/(1._wp+zcfm(jl,itop)*zqdp)
        zebsm(jl,itop)=zdisc*(zcfm(jl,itop)*zqdp)
        zudif(jl,itop)=zdisc*zudif(jl,itop)
        zvdif(jl,itop)=zdisc*zvdif(jl,itop)
421  END DO
!
!*         4.3     ELIMINATION FOR MIDDLE LAYERS.
!
     DO 432 jk=itopp1,klevm1
        DO 431 jl=1,kproma
           zqdp=1._wp/(paphm1(jl,jk+1)-paphm1(jl,jk))
           zfac=zcfm(jl,jk-1)*zqdp
           zdisc=1._wp/(1._wp+zfac*(1._wp-zebsm(jl,jk-1))              &
                                                   +zcfm(jl,jk)*zqdp)
           zebsm(jl,jk)=zdisc*(zcfm(jl,jk)*zqdp)
           zudif(jl,jk)=zdisc*(zudif(jl,jk)+zfac*zudif(jl,jk-1))
           zvdif(jl,jk)=zdisc*(zvdif(jl,jk)+zfac*zvdif(jl,jk-1))
431     END DO
432  END DO
!
!*         4.4     BOTTOM LAYER ELIMINATION.
!
     DO 441 jl=1,kproma
        zqdp=1._wp/(paphm1(jl,klevp1)-paphm1(jl,klev))
        zfac=zcfm(jl,klevm1)*zqdp
        zdisc=1._wp/(1._wp+zfac*(1._wp-zebsm(jl,klevm1))               &
                                                  +zcfm(jl,klev)*zqdp)
        zudif(jl,klev)=zdisc*(zudif(jl,klev)+zfac*zudif(jl,klevm1))
        zvdif(jl,klev)=zdisc*(zvdif(jl,klev)+zfac*zvdif(jl,klevm1))
441  END DO
!
!*         4.5     BACK-SUBSTITUTION.
!
     DO 452 jk=klevm1,itop,-1
        DO 451 jl=1,kproma
           zudif(jl,jk)=zudif(jl,jk)+zebsm(jl,jk)*zudif(jl,jk+1)
           zvdif(jl,jk)=zvdif(jl,jk)+zebsm(jl,jk)*zvdif(jl,jk+1)
451     END DO
452  END DO
!
!*         4.6     INCREMENTATION OF U AND V TENDENCIES AND STORAGE OF
!*                 THE DISSIPATION.
!
     DO 461 jl=1,kproma
        zvidis(jl)=0._wp
461  END DO
!
     DO 471 jk=itop,klev
        DO 462 jl=1,kproma
           zdudt=(zudif(jl,jk)-ztpfac2*pum1(jl,jk))*zcons13
           pvom(jl,jk)=pvom(jl,jk)+zdudt
           zdvdt=(zvdif(jl,jk)-ztpfac2*pvm1(jl,jk))*zcons13
           pvol(jl,jk)=pvol(jl,jk)+zdvdt
           zdis(jl,jk)=0.5_wp*((ztpfac2*pum1(jl,jk)-zudif(jl,jk))      &
                              *(ztpfac4*pum1(jl,jk)+zudif(jl,jk))      &
                              +(ztpfac2*pvm1(jl,jk)-zvdif(jl,jk))      &
                              *(ztpfac4*pvm1(jl,jk)+zvdif(jl,jk)))
           zvidis(jl)=zvidis(jl)+                                      &
                      zdis(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
462     END DO
471  END DO
! HAND IT OVER (JSBACH - MO_SURFACE----------------------

     DO jl=1,kproma
        udif(jl,jrow) = zudif(jl,klev)
        vdif(jl,jrow) = zvdif(jl,klev)
! BLM dissipation -----------------------------------------
        pvdis(jl)=pvdis(jl)+zdtime*zcons15*zvidis(jl)

     END DO

!
!*         5.6     INCREMENTATION OF T AND Q TENDENCIES.
!
     DO 571 jk=itop,klev
        DO 561 jl=1,kproma
           zqdif(jl,jk)=zqdif(jl,jk)+ztpfac3*pqm1(jl,jk)
           zdqdt=(zqdif(jl,jk)-pqm1(jl,jk))*zcons13
           pqte(jl,jk)=pqte(jl,jk)+zdqdt
           ztdif(jl,jk)=ztdif(jl,jk)+ztpfac3*zcptgz(jl,jk)
           zdtdt=((ztdif(jl,jk)+zdis(jl,jk)-pgeom1(jl,jk))             &
                 /(cpd*(1._wp+vtmpc2*pqm1(jl,jk)))-ptm1(jl,jk))*zcons13
           zcpten(jl,jk) = (ztdif(jl,jk)-zcptgz(jl,jk))*zcons13
           ptte(jl,jk)=ptte(jl,jk)+zdtdt
           zqten(jl,jk)   = zdqdt
           zxldif(jl,jk)=zxldif(jl,jk)+ztpfac3*pxlm1(jl,jk)
           zxidif(jl,jk)=zxidif(jl,jk)+ztpfac3*pxim1(jl,jk)
           zdxlmdt=(zxldif(jl,jk)-pxlm1(jl,jk))*zcons13
           zdximdt=(zxidif(jl,jk)-pxim1(jl,jk))*zcons13
           pxlte(jl,jk)=pxlte(jl,jk)+zdxlmdt
           pxite(jl,jk)=pxite(jl,jk)+zdximdt
561     END DO
571  END DO
!
     ! Correction of tte for snow melt
     DO jl = 1,kproma
        ptte(jl,klev)=ptte(jl,klev)-ztte_corr(jl)
     END DO
!
     IF (trlist% anyvdiff /= 0) THEN
        DO jt=1,trlist% ntrac
           IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
!+++sschr: Bug #329
!          (initialize zdxtdt inside(!) the tracer loop)
           zdxtdt = 0._wp
!---sschr: Bug #329
           DO jk=itop,klev
              DO jl=1,kproma
                 zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt)+                    &
                                  ztpfac3*pxtm1(jl,jk,jt)
                 zdxtdt(jl,jk) = zxtdif(jl,jk,jt)-pxtm1(jl,jk,jt)
              END DO
           END DO
           IF (lco2_mixpbl .AND. jt == ico2idx) THEN    ! ideal mixing in the PBL for CO2
              DO jl=1,kproma
                 ktop = MIN(MAX(ihpbl(jl),iblmax),iblmin)
                 zdxtdtsum = 0._wp
                 DO jk=ktop,klev
                    zdxtdtsum = zdxtdtsum + zdxtdt(jl,jk) *            &
                                        (paphm1(jl,jk+1)-paphm1(jl,jk))
                 END DO
                 zdxtdtmean = zdxtdtsum /                              &
                                    (paphm1(jl,klevp1)-paphm1(jl,ktop))
                 DO jk=ktop,klev
!+++sschr: Issue #332
!          (avoid negative concentrations)
                   zdxtdt(jl,jk) = MAX(zdxtdtmean,-zxt(jl,jk,jt))
!---sschr: Issue #332
                 END DO
              END DO
           END IF
           pxtte(1:kproma,1:klev,jt) = pxtte(1:kproma,1:klev,jt) +     &
                                      zdxtdt(1:kproma,1:klev) * zcons13
        END DO
     END IF

#ifdef _PROFILE
  CALL trace_stop ('vdiff_loop_7', 27)
  CALL trace_start ('vdiff_loop_8', 28)
#endif
!
! back out moisture flux
!
    DO jl=1,kproma
       zdqtot=(pqm1(jl,klev)+zx(jl,klev))-zqsurf_new(jl)
    ENDDO !jl
    DO jk=itop+1,klevp1
       DO jl=1,kproma
          IF (jk<klevp1) THEN
             ztvh=(ptm1(jl,jk)*(1._wp+vtmpc1*pqm1(jl,jk)-zx(jl,jk))    &
                  +ptm1(jl,jk-1)*(1._wp+vtmpc1*pqm1(jl,jk-1)           &
                   -zx(jl,jk-1)))/2._wp
          ENDIF
       ENDDO !jl
    ENDDO !jk
    DO jk=itop,klev
       DO jl=1,kproma
          zhexp=EXP(1._wp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
          zlam=zzzlam+(zcons3-zzzlam)*zhexp
          IF(jk.GE.ihpbl(jl)) THEN
             zcons23=zcons25
          ELSE
             zcons23=zcons2/zlam
          END IF
          z2geomf=2._wp*pgeom1(jl,jk)
          zz2geo=zcons2*z2geomf
          zmix=zz2geo/(1._wp+zcons23*z2geomf)
          IF(jk.EQ.1) THEN
             ztkesq=SQRT(MAX(ztkemin,ptkem1(jl,1)))
          ELSE
             ztkesq=SQRT(MAX(ztkemin,0.5_wp*(ptkem1(jl,jk-1)           &
                                           +ptkem1(jl,jk))))
          END IF
       ENDDO !jl
    ENDDO !jk
!
!
! column heating and moistening due to vertical diffusion
! (compare locally to sensible heat flux and evaporation)
  psh_vdiff(1:kproma) = 0.0_wp
  pev_vdiff(1:kproma) = 0.0_wp
  DO jk = itop,klev
     DO jl = 1,kproma
        psh_vdiff(jl)=psh_vdiff(jl)+zcpten(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))/grav
        pev_vdiff(jl)=pev_vdiff(jl)+zqten(jl,jk) *(paphm1(jl,jk+1)-paphm1(jl,jk))/grav
     END DO
  END DO

!     ------------------------------------------------------------------
!
!*         6.     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
!
  ELSE
     DO  601 jl=1,kproma
        pevapot(jl) = 0.0_wp
        pqhfla(jl)  = 0.0_wp
        ppbl(jl)    = 0.0_wp
601  END DO
  END IF
!
!     ------------------------------------------------------------------
!
#ifdef _PROFILE
  CALL trace_stop ('vdiff', 20)
#endif
!
END SUBROUTINE vdiff
