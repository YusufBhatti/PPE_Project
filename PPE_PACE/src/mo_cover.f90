#ifdef __xlC__
@PROCESS HOT
@PROCESS XLF90(NOSIGNEDZERO)
#else
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#define SWDIV_NOCHK(a,b) ((a)/(b))
#endif

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module diagnoses cloud cover for current timestep
!!
!! @remarks
!!     In ECHAM4 all cloud cover calculations were performed
!!     in routine CLOUD.  The routine diagnosed the cloud cover
!!     from the previous timestep using m1 variables, but then
!!     used a "first guess" calculation of T and Q to give a
!!     cloud cover estimate for the next timestep.  This meant
!!     that the radiation scheme used cloud cover values that were
!!     from a different timestep to the temperature and water vapour
!!     values, and also that the T and Q values were anyway
!!     preliminary.  Finally, the cover calculation was performed
!!     twice each timestep, when one calculation suffices.
!!
!!     This scheme calculates cover diagnostically and is called
!!     once at the beginning of each timestep.  It uses the
!!     standard relative humidity calculation from Lohmann and
!!     Roeckner (96), or the method from the new prognostic
!!     scheme of Tompkins.  The choice of which scheme to use is
!!     controlled by the parameter switch ICOVER, which is set in
!!     namelist PHYSCTL along with lsurf etc... Note that even if
!!     icover.EQ.1 (RH scheme) you can't restart this model version
!!     from restart files saved from a different model version, since
!!     the two extra prognostic equations, pxvar and pxskew are still
!!     stored even though they are not actively used.  However, this means
!!     that once you have restart files from this version, you are able
!!     change icover at will.
!!
!!     In the new scheme the variable xskew is provided
!!     as outlined in the reference, this variable represents
!!     directly the Beta distribution shape parameter "q"
!!     The shape parameter "p" (zbetap) a tunable parameter and it is
!!     recommended that this be set to a low constant 1.5<p<2.0
!!     (This may be changed later to prognostic to allow negative skewness
!!     from downdraft detrainment, see ref. for details).
!!
!!     from xi,xl,q,xskew and zbetap, the Beta distribution is defined
!!     and cloud cover is diagnosable.  For the iteration, Ridders' method
!!     is used (see Numerical Recipes).
!!
!!     Attention:
!!     In the current version the advective tendencies of skewness
!!     and variance are set to zero.
!!
!! @references.
!!     Diagnostic CC scheme: Lohmann and Roeckner 96, Clim. Dyn.
!!     Prognostic CC scheme: Tompkins 2002, J. Atmos. Sci.
!!
!! @author A. Tompkins    MPI-Hamburg        2000
!!         K. Ketelesen   NEC,         April 2002
!!         L. Kornblueh   MPI-Hamburg, April 2002
!!
!! @par Revision History
!!    v2: first working version
!!    v5: lookup table added
!!    v8: zriddr and functions replaced for vectorization
!!    v9: optimizations, longer vector loop and less indirect addressing
!!       - introduction of additional arrays
!!       - change structure if "beta function scheme" IF block
!!       - scattered loops "ictit" and "ictdg" are collected over "kproma" and "klev"
!!       - intoducing additional arrays to hold data in "ictit" and "ictdg"
!!         addressing scheme
!! - Taken from ECHAM6.2, wrapped in module and modified for ICON
!!   by Monika Esch, MPI-M (2013-11)
!! - updated to ECHAM.6.3 and unified for use in ICON; removed Tompkins scheme
!!   by Monika Esch, MPI-M (2015-05)
!!
!
MODULE mo_cover

  USE mo_kind,                 ONLY : wp
  USE mo_physical_constants,   ONLY : vtmpc1, cpd, grav
  USE mo_echam_convect_tables, ONLY : prepare_ua_index_spline,lookup_ua_eor_uaw_spline
  USE mo_echam_cloud_params,   ONLY : jbmin, jbmax, csatsc, crt, crs, nex, nadd, cinv
#ifdef _PROFILE
  USE mo_profile,              ONLY : trace_start, trace_stop
#endif
  USE mo_p3_fields,            ONLY : iprog, ccclpwr, ccsupci, lallow_si
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cover, cover2

CONTAINS
  !>
  !!
  SUBROUTINE cover (         kproma,   kbdim, ktdia, klev, klevp1                  & !in
    &                      , ktype,    pfrw,     pfri                              & !in
    &                      , paphm1,   papm1,    pgeo                              & !in
    &                      , ptm1,     pqm1,     pxim1                             & !in
    &                      , paclc                                                 & !inout
    &                      , knvb,     printop                                     & !out
    &              )
    !---------------------------------------------------------------------------------
    !

    USE mo_echam_cloud_params,   ONLY: cthomi
    USE mo_cirrus,               ONLY: SCRHOM
    USE mo_physical_constants,   ONLY: tmelt

    INTEGER, INTENT(IN)    :: kbdim, klevp1, klev, kproma, ktdia
    INTEGER, INTENT(IN)    ::  &
      & ktype(kbdim)          !< type of convection
    REAL(wp),INTENT(IN)    ::  &
      & pfrw(kbdim)         ,&!< water mask
      & pfri(kbdim)           !< ice mask
    REAL(wp),INTENT(IN)    ::  &
      & paphm1(kbdim,klevp1),&!< pressure at half levels                   (n-1)
      & papm1(kbdim,klev)   ,&!< pressure at full levels                   (n-1)
      & pgeo(kbdim,klev)    ,&!<
      & pqm1(kbdim,klev)    ,&!< specific humidity                         (n-1)
      & ptm1(kbdim,klev)    ,&!< temperature                               (n-1)
      & pxim1(kbdim,klev)     !< cloud ice                                 (n-1)
    REAL(wp),INTENT(INOUT) ::  &
      & paclc(kbdim,klev)     !< cloud cover
    INTEGER, INTENT(OUT)   ::  &
      & knvb(kbdim)
    REAL(wp),INTENT(OUT)   ::  &
      & printop(kbdim)

    INTEGER :: jl, jk, jb, jbn
    INTEGER :: locnt, nl, ilev
    REAL(wp):: zdtdz, zcor, zrhc, zsat, zqr
    INTEGER :: itv1(kproma*klev), itv2(kproma*klev)

    !
    !   Temporary arrays
    !
    REAL(wp)   ::  zdtmin(kbdim), za(kbdim) 
    !
    !   Pointers and counters for iteration and diagnostic loop:
    !
    INTEGER :: nphase
    !
    !   variables required for zriddr iteration scheme:
    !
    REAL(wp) :: zjk, zgam
    REAL(wp) :: zqsm1(kproma*klev)
    REAL(wp) :: ua(kproma)

    LOGICAL :: lao, lao1

    REAL(wp) :: zpapm1i(kbdim,klev),     ztmp(kproma*klev)

    REAL(wp) :: zknvb(kbdim),            zphase(kbdim)

    INTEGER :: loidx(kproma*klev)
    REAL(wp) :: w, f, rh_koop, qsw, qsi, rh_wi, rh_upper, rh_lower, qs, rh

#ifdef _PROFILE
    CALL trace_start ('cover', 9)
#endif
    !
    !   Initialize variables
    !
    DO jl = 1,kproma
      zdtmin(jl) = -cinv * grav/cpd   ! fraction of dry adiabatic lapse rate
      zknvb(jl)  = 1.0_wp
      printop(jl)= 0.0_wp
    END DO
    !
    DO jk = ktdia,klev
      DO jl = 1,kproma
         zpapm1i(jl,jk) = SWDIV_NOCHK(1._wp,papm1(jl,jk))
      END DO
    END DO
    !
    !       1.3   Checking occurrence of low-level inversion
    !             (below 2000 m, sea points only, no convection)
    !  
    locnt = 0
    DO jl = 1,kproma
      IF (pfrw(jl).GT.0.5_wp.AND.pfri(jl).LT.1.e-12_wp.AND.ktype(jl).EQ.0) THEN
        locnt = locnt + 1
        loidx(locnt) = jl
      END IF
    END DO

    IF (locnt.GT.0) THEN
      DO jk = klev,jbmin,-1

!IBM* ASSERT(NODEPS)
!IBM* novector
        DO nl = 1,locnt
          jl = loidx(nl)
          ztmp(nl) = (ptm1(jl,jk-1)-ptm1(jl,jk))*grav/(pgeo(jl,jk-1)-pgeo(jl,jk))
        END DO

        zjk = REAL(jk,wp)
!IBM* ASSERT(NODEPS)
        DO nl = 1,locnt
          jl = loidx(nl)
          zdtdz       = MIN(0.0_wp, ztmp(nl))
          zknvb(jl)   = FSEL(zdtmin(jl)-zdtdz,zknvb(jl),zjk)
          zdtmin(jl)  = MAX(zdtdz,zdtmin(jl))
        END DO
      END DO
    END IF
    knvb(1:kproma) = INT(zknvb(1:kproma))
    !
    !       1.   Calculate the saturation mixing ratio
    !
    IF (ktdia < klev+1) THEN

      DO jk = ktdia,klev

        CALL prepare_ua_index_spline('cover (2)',kproma,ptm1(1,jk),itv1(1),          &
                                         za(1),pxim1(1,jk),nphase,zphase,itv2)
        CALL lookup_ua_eor_uaw_spline(kproma,itv1(1),za(1),nphase,itv2(1),ua(1))

!IBM* novector
        DO jl = 1,kproma
           ! RDnote: saturation mixing ratio analogously to cloud micro
           !zqsm1(jl) = saturation_mixing_ratio(ptm1(jl,jk), papm1(jl,jk))
           
           IF(iprog == 3) THEN
              zqsm1(jl) = MIN(ua(jl)*zpapm1i(jl,jk),0.5_wp)
              zcor      = 1._wp/(1._wp-vtmpc1*zqsm1(jl))
              zqsm1(jl) = zqsm1(jl)*zcor
           ELSE
              ! linear saturation mixing ratio between -35°C and 0°C (Morrison & Gettelman 2008)
              zqsm1(jl) = saturation_mixing_ratio(ptm1(jl,jk), papm1(jl,jk))
           ENDIF !iprog

          
        END DO
        !
        !       Threshold relative humidity, qsat and cloud cover
        !       This is from cloud, and is the original calculation for
        !       cloud cover, based on relative humidity
        !       (Lohmann and Roeckner, Clim. Dyn.  96)
        !
        DO jl = 1,kproma
        !
          zrhc=crt+(crs-crt)*EXP(1._wp-(paphm1(jl,klevp1)/papm1(jl,jk))**nex)
          zsat=1._wp
          jb=knvb(jl)
          jbn=jb+nadd                  ! mo_echam_cloud_params: nadd=0 except for T31
          lao=(jb.GE.jbmin .AND. jb.LE.jbmax)
          lao1=(jk.EQ.jb .OR. jk.EQ.jbn)
          ilev=klev
          IF (lao .AND. lao1) THEN
          !  ilev=klevp1-jb
            ilev=100
            printop(jl)=REAL(ilev,wp)
            zdtdz = (ptm1(jl,jb-1)-ptm1(jl,jb))*grav/(pgeo(jl,jb-1)-pgeo(jl,jb))
            zgam  = MAX(0.0_wp,-zdtdz*cpd/grav)
            zsat  = MIN(1.0_wp,csatsc+zgam)
         END IF
          ! Modified cloud cover scheme. Note that for warm clouds, it simplifies to 
          ! the original cloud cover scheme.

          ! cirrus flag
          w = MERGE(1._wp, 0._wp, ptm1(jl,jk) < cthomi .OR. lallow_si)

          ! mixed-phase ramp
          f = (ptm1(jl,jk) - tmelt)/(cthomi - tmelt)
          f = MIN(MAX(f, 0._wp), 1._wp)

          ! get supersaturation needed for homogeneous nucleation
          rh_koop = SCRHOM(ptm1(jl,jk))

          ! ice and liquid saturation humidities
          qsw = liquid_saturation(ptm1(jl,jk), papm1(jl,jk))
          qsi = ice_saturation(ptm1(jl,jk), papm1(jl,jk))

          ! relative humidity at liquid saturation
          rh_wi = qsw/qsi

          ! use minium of rh_iw and rh_koop as upper limit where b==1. Assert rh_upper >= 1.
          rh_upper = MAX(MIN(rh_koop, rh_wi), 1._wp)
          ! linearly adjust lower limit to 1 in cirrus regime
          rh_lower = f*1 + (1._wp - f)*zrhc

          ! define grid-box sat. spec. hum. w.r.t. liquid for T > 0 and w.r.t. ice for T < 0
          qs = MERGE(qsw, qsi, ptm1(jl,jk) > tmelt)
          
          ! ... and compute RH w.r.t. above saturation humidity (use ice in cirrus)
          rh = (pqm1(jl,jk)+w*pxim1(jl,jk))/(qs*zsat)
         
          ! Sundqvist formula
!>>DN bugfix: division by 0
!          paclc(jl,jk)=(rh-rh_lower)/(rh_upper-rh_lower)
          IF (rh_upper>rh_lower) THEN
             paclc(jl,jk)=(rh-rh_lower)/(rh_upper-rh_lower)
          ELSE
             paclc(jl,jk)=0._wp             
          END IF
!<<DN bugfix
          paclc(jl,jk)=MAX(MIN(paclc(jl,jk),1.0_wp),0.0_wp)
          paclc(jl,jk)=1._wp-SQRT(1._wp-paclc(jl,jk))
        END DO !jl
      END DO  !jk
    END IF

#ifdef _PROFILE
    CALL trace_stop ('cover', 9)
#endif
    !
    !
  END SUBROUTINE cover

SUBROUTINE cover2(&
              !--IN
              kproma, kbdim, ktdia, klev, &
              ptm1, papm1, paphm1, pqm1, pxim1, &
              !--OUT
              paclc)

  USE mo_echam_cloud_params,   ONLY: cthomi, cqtmin
  USE mo_cirrus,               ONLY: SCRHOM

  INTEGER, INTENT(IN) :: kproma, kbdim, ktdia, klev
  REAL(wp), DIMENSION(kbdim,klev), INTENT(IN) :: ptm1, papm1, pqm1, pxim1
  REAL(wp), DIMENSION(kbdim,klev+1), INTENT(IN) :: paphm1
  REAL(wp), DIMENSION(kbdim,klev), INTENT(OUT) :: paclc

  INTEGER :: jl, jk
  REAL(wp) :: qsw, qsi
  LOGICAL :: ll_cover

  DO jl = 1,kproma
     DO jk = ktdia,klev
        ! ice and liquid saturation humidities
        qsw = liquid_saturation(ptm1(jl,jk), papm1(jl,jk))
        qsi = ice_saturation(ptm1(jl,jk), papm1(jl,jk))

        ! cover if at or super saturation w.r.t. liquid (ice if ice is present)
        ll_cover = pqm1(jl,jk) >= qsw .OR. &
                  (pqm1(jl,jk) >= qsi .AND. pxim1(jl,jk) > cqtmin)

        ! set cover
        paclc(jl,jk) = MERGE(1._wp, 0._wp, ll_cover)
     END DO
  END DO

END SUBROUTINE cover2

FUNCTION liquid_saturation(t, p) RESULT(qs)

  USE mo_echam_convect_tables, ONLY: lookuperror, jptlucu1,    &
                                     jptlucu2, tlucua, tlucub, tlucuaw
  USE mo_echam_cloud_params,   ONLY: cthomi
  USE mo_physical_constants,   ONLY: vtmpc1, tmelt

  REAL(wp) :: t         !< temperature
  REAL(wp) :: p         !< pressure
  REAL(wp) :: qs        !< saturation mixing ratio

  REAL(wp) :: dum
  REAL(wp) :: frac
  INTEGER  :: idx       !< temperature index in lookup table

  dum = 1000._wp*t
  idx = NINT(dum)

  IF (idx<jptlucu1 .OR. idx>jptlucu2) CALL lookuperror ('mo_cover: saturation_mixing_ratio')

  idx = MAX(MIN(idx,jptlucu2),jptlucu1)

  dum = tlucuaw(idx)/p

  dum = MIN(dum, 0.5_wp)
  qs  = dum/(1._wp-vtmpc1*dum)

END FUNCTION liquid_saturation

FUNCTION ice_saturation(t, p) RESULT(qs)

  USE mo_echam_convect_tables, ONLY: lookuperror, jptlucu1,    &
                                     jptlucu2, tlucua, tlucub, tlucuaw
  USE mo_echam_cloud_params,   ONLY: cthomi
  USE mo_physical_constants,   ONLY: vtmpc1, tmelt

  REAL(wp) :: t         !< temperature
  REAL(wp) :: p         !< pressure
  REAL(wp) :: qs        !< saturation mixing ratio

  REAL(wp) :: dum
  REAL(wp) :: frac
  INTEGER  :: idx       !< temperature index in lookup table

  dum = 1000._wp*t
  idx = NINT(dum)

  IF (idx<jptlucu1 .OR. idx>jptlucu2) CALL lookuperror ('mo_cover: saturation_mixing_ratio')

  idx = MAX(MIN(idx,jptlucu2),jptlucu1)

  dum = tlucua(idx)/p

  dum = MIN(dum, 0.5_wp)
  qs  = dum/(1._wp-vtmpc1*dum)

END FUNCTION ice_saturation

FUNCTION saturation_mixing_ratio(t, p) RESULT(qs)

  USE mo_echam_convect_tables, ONLY: lookuperror, jptlucu1,    &
                                     jptlucu2, tlucua, tlucub, tlucuaw
  USE mo_echam_cloud_params,   ONLY: cthomi
  USE mo_physical_constants,   ONLY: vtmpc1, tmelt

  REAL(wp) :: t         !< temperature
  REAL(wp) :: p         !< pressure
  REAL(wp) :: qs        !< saturation mixing ratio

  REAL(wp) :: dum
  REAL(wp) :: qsi, qsw
  REAL(wp) :: frac
  INTEGER  :: idx       !< temperature index in lookup table

  dum = 1000._wp*t
  idx = NINT(dum)

  IF (idx<jptlucu1 .OR. idx>jptlucu2) CALL lookuperror ('mo_cover: saturation_mixing_ratio')

  idx = MAX(MIN(idx,jptlucu2),jptlucu1)
  
  frac = (t-tmelt)/(cthomi-tmelt)
  frac = MIN(MAX(0._wp, frac), 1._wp)
  frac = frac**ccclpwr

  qsi = tlucua(idx)/p
  qsi = MIN(qsi, 0.5_wp)
  qsi = qsi/(1._wp-vtmpc1*qsi)

  qsw = tlucuaw(idx)/p
  qsw = MIN(qsw, 0.5_wp)
  qsw = qsi/(1._wp-vtmpc1*qsw)

  qs = qsi*frac + qsw*(1._wp-frac)

END FUNCTION saturation_mixing_ratio

FUNCTION saturation_mixing_ratio2(t, p, xi) RESULT(qs)

  USE mo_echam_convect_tables, ONLY: lookuperror, jptlucu1,    &
                                     jptlucu2, tlucua, tlucub, tlucuaw
  USE mo_echam_cloud_params,   ONLY: cthomi, csecfrl
  USE mo_physical_constants,   ONLY: vtmpc1, tmelt

  REAL(wp) :: t         !< temperature
  REAL(wp) :: p         !< pressure
  REAL(wp) :: xi        !< ice mixing ratio
  REAL(wp) :: qs        !< saturation mixing ratio

  REAL(wp) :: dum
  REAL(wp) :: frac
  INTEGER  :: idx       !< temperature index in lookup table

  dum = 1000._wp*t
  idx = NINT(dum)

  IF (idx<jptlucu1 .OR. idx>jptlucu2) CALL lookuperror ('mo_cover: saturation_mixing_ratio_ice')

  idx = MAX(MIN(idx,jptlucu2),jptlucu1)

  IF(t < cthomi) THEN
     dum = tlucua(idx)/p
  ELSE
     dum = tlucuaw(idx)/p
  ENDIF

  dum = MIN(dum, 0.5_wp)
  qs  = dum/(1._wp-vtmpc1*dum)

END FUNCTION saturation_mixing_ratio2

FUNCTION saturation_mixing_ratio3(t, p) RESULT(qs)

  USE mo_echam_convect_tables, ONLY: lookuperror, jptlucu1, jptlucu2, tlucua, tlucuaw
  USE mo_physical_constants,   ONLY: vtmpc1, tmelt

  REAL(wp) :: t        !< temperature
  REAL(wp) :: p        !< pressure
  REAL(wp) :: qs       !< saturation mixing ratio
  REAL(wp) :: qsl      !< saturation specific humidity with respect to water
  REAL(wp) :: qsi      !< saturation specific humidity eith respect to ice

  REAL(wp) :: rto      !< mixing ratio of the two humidities (rto=0.: 100% qsi, rto=1.: 100% qsl)
  INTEGER  :: idx      !< temperature index in lookup table
  REAL(wp) :: dum

  dum = 1000._wp*t
  idx = NINT(dum)

  IF (idx<jptlucu1 .OR. idx>jptlucu2) CALL lookuperror ('mo_cover: saturation_mixing_ratio_ice')

  idx = MAX(MIN(idx,jptlucu2),jptlucu1)

  qsl = 1._wp/((p/tlucuaw(idx))-vtmpc1)
  qsi = 1._wp/((p/tlucua(idx))-vtmpc1)

  IF(t > tmelt) THEN
     rto = 1._wp
  ELSE IF(t < tmelt-35._wp) THEN
     rto = 0._wp
 ELSE
     rto = (t-tmelt)/35._wp+1._wp
  END IF

  qs = rto*qsl+(1._wp-rto)*qsi

END FUNCTION saturation_mixing_ratio3

END MODULE mo_cover
