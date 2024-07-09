!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_update_surf_up

PRIVATE

PUBLIC update_surf_up

CONTAINS

SUBROUTINE update_surf_up ( kbdim, nsoil, lperma &
     , pwl, pcvw, pwlmx &
     , pws, pwsmx &
     , ptsoil, pcvs          &
     , porostd &
     , pevapot, pevwsd                     &
     , prain                           &
     , lmask, lpglac                   &
     , pros_hd,    pdrain_hd                   &
     , pmlres &
     , zevwld &
!     Soil moisture layers, evap. fluxes over land
     , zdel, wsi, zwsfc           &
     , dzr, dzs, fksat , fmpot , bclapp, vpor     &
     , etrans, redevap                    &
     , vfc, vpwp, spsi                           &
     , jllog, ice                                &
     )
!
!     ------------------------------------------------------------------
! pwl           Water content [m] in skin reservoir
!               (vegetation and bare soil)
! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
!               of vegetation index and leaf area index)
! pcvw          Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
! pws           Soil water content [m]
! pwsmx         Water holding capacity [m] of the soil
! ptsoil        Temperature [K] in upper soil layer
! pcvs          Fractional snow cover (function of psn in *physc*)
! porostd       Subgrid standard diviation [m] used in runoff scheme
! etrans        Transpiration [kg/m**2/s]
! pevapot       Potential evaporation/sublimation [kg/m**2/s]
! pevwsd        Evapotranspiration without sublimation and evaporation from interception reservoir [kg/m**2]
!               BUT in this routine it is [m], such as zevwld.
! prain         Totalrainfall [kg/m**2/s]
! pros_hd       Runoff for HD-Model (does NOT include drainage) [m]
! pdrain_hd     Drainage for HD-Model [m]
! lpglac        Logical glacier mask
!
! The following local variables represent the respective fluxes
! integrated over one timestep (delta_time) and divided by the density of
! water (rhoh2o). Units are m water equivalent.
!
! zraind        Total rain [m] 
! zevwld        Evaporation from the skin reservoir [m]
! zros          Total runoff (including drainage) at non-glacier points
! zdrain        Drainage at non-glacier points
! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
! pmlres        Residual melt water available for infiltration into the
!               non-frozen soil after being intercepted by the
!               skin reservoir
!  
!  ******** multi-layer scheme
!     ***  kdeep = Number of soil layers
!     ***   zdel = Thickness of the soil layers [m]
!     ***  dzrsi = Rooted depth per soil layer (until rooting depth DZR) [m]
!     ***   dzsi = Depth that can be saturated with water per soil layer
!     ***          (until bedrock DZS) [m]
!     ***    dzr = Rooting depth in [m]
!     ***    dzs = Soil depth derived from textures in [m]
!     ***  fksat = Saturated hydraulic conductivity: unit [m/s]
!     ***  fmpot = Matrix potential [m]
!     *** bclapp = Exponent B in Clapp and Hornberger 
!     ***   vpor = Volumetric soil porosity [m/m]   following Cosby et al.
!     ***    vfc = Volumetric soil field capacity [m/m]   
!     ***   vpwp = Volumetric soil wilting point [m/m]   
!     ***   spsi = Soil Pore size distribution index    
!     ***  zwsat = MAXIMALE SPEICHERGEHALT DER SOIL LAYER I (POROSITAET) [M]
!     ***          WICHTIG FUER AUFNAHME DER INFILTRATION IN DEN BODENSCHICHTEN.
!     ***          IST DEFINIERT BIS ZUR TIEFE DES BODENS dzs, WO DER BEDROCK BEGINNT.
!     ***  zwsfc = Field capacity of soil layer I (analog zu ZWSAT) [m]
!     ***  zwpwp = Soil wilting point of soil layer I (analog to ZWSAT) [m]
!     ** aebsoil = Bare soil evaporation in [m], accumulated over TWODT/2
!     ** aetrans = Transpiration in [m] accumulated over one time step 
!     ***          (calculated from etrans which is in [kg/m**2/s])
!     *** redevap  = Diagnostic desired reduction of evapotranspiration
!     ***            due to limited storage in soil layers [kg/m**2/s]
!     *** wsi(I) = Soil moisture content in I. soil layer [m]
!     ***          Is updated with the fluxes calculated in SOILHYD
!     *** wsim1m(i) = Soil Moisture of previous time step in i. Soil layer [m] 
!     *** ZDRAIN = Drainage: Volumen in time step: unit [m]
!
!     *** jllog          relative Index of grid box for log output between 1 and kbdim
!
!       The rest is to be added later on ....
!
!     *SURF* - Updates land values of temperature, moisture and snow.
!              Calculate fluxes of total rain, total snow and evapo-
!              ration from the three reservoirs (SN, WS, WL)
!              convert fluxes (kg/m**2*s) to changes of water levels (m)
!              during timestep delta_time.
!
!     J.F.GELEYN     E.C.M.W.F.     08/06/82.
!     MODIFIED BY
!     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
!     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
!                 J.-P. SCHULZ   MPI - 1997 : Implementation of implicit coupling between 
!                                             land surface and atmosphere.
!     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
!     MODIFIED BY M. ESCH        MPI - APR  1999
!     MODIFIED BY E. ROECKNER    MPI - JAN  2001
!     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
!     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir 
!                                                for snow changed
!     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
!     MODIFIED BY ??? R. SCHNUR?        JSBACH Version
!
!     MODIFIED BY S. HAGEMANN    MPI - February  2009 
!     *** Implementation of 5 layer soil hydrology scheme
!
!     MODIFICATION
!
!     PURPOSE
!
!     INTERFACE.
!
!          *SURF* is called from *PHYSC*.
!
!     METHOD.
!
!     EXTERNALS.
!
!          NONE.
!
!     REFERENCE.
!
!          See soil processes' part of the model's documentation for
!     details about the mathematics of this routine.
!
  USE mo_param_switches,   ONLY: lsurf
  USE mo_jsbach_constants, ONLY: cvinter, rhoh2o, tmelt
  USE mo_time_control,     ONLY: delta_time
  USE mo_kind,             ONLY: dp
#ifdef STANDALONE
  USE mo_jsbach_comm_to_echam5mods, ONLY: nlat
#else
  USE mo_control,        ONLY: ngl
#endif

IMPLICIT NONE
 
  INTEGER,  INTENT(in)    :: kbdim
  INTEGER,  INTENT(in)    :: nsoil ! also used as switch for bucket (nsoil=1) vs. multi-layer (nsoil>1) calculation
  LOGICAL,  INTENT(in)    :: lperma
  LOGICAL,  INTENT(in)    :: lmask(kbdim)
  LOGICAL,  INTENT(in)    :: lpglac(kbdim)
  REAL(dp), INTENT(in)    :: pcvw(kbdim)
  REAL(dp), INTENT(in)    :: pwlmx(kbdim)
  REAL(dp), INTENT(in)    :: pwsmx(kbdim)
  REAL(dp), INTENT(in)    :: ptsoil(kbdim) 
  REAL(dp), INTENT(in)    :: pcvs(kbdim)
  REAL(dp), INTENT(in)    :: porostd(kbdim)
  REAL(dp), INTENT(in)    :: pevapot(kbdim)
  REAL(dp), INTENT(in)    :: prain(kbdim)
  REAL(dp), INTENT(in)    :: pmlres(kbdim)

  REAL(dp), INTENT(inout) :: pwl(kbdim)
  REAL(dp), INTENT(inout) :: pws(kbdim)
  REAL(dp), INTENT(inout) :: pevwsd(kbdim)

  REAL(dp), INTENT(inout) :: pros_hd(kbdim)
  REAL(dp), INTENT(inout) :: pdrain_hd(kbdim)
  REAL(dp), INTENT(out)   :: zevwld(kbdim) ! From optional part

  REAL(dp), INTENT(in), optional :: ice(kbdim,nsoil) ! Needed for permafrost model
!
! Soil moisture layers, evap. fluxes over land.
  REAL(dp), INTENT(in)   , optional :: zdel(nsoil)
  REAL(dp), INTENT(inout), optional :: wsi(kbdim,nsoil)
  REAL(dp), INTENT(inout), optional :: zwsfc(kbdim,nsoil)
  REAL(dp), INTENT(in)   , optional :: dzr(kbdim)
  REAL(dp), INTENT(in)   , optional :: dzs(kbdim)
  REAL(dp), INTENT(in)   , optional :: fksat(kbdim)
  REAL(dp), INTENT(in)   , optional :: fmpot(kbdim)
  REAL(dp), INTENT(in)   , optional :: bclapp(kbdim)
  REAL(dp), INTENT(in)   , optional :: vpor(kbdim)
  REAL(dp), INTENT(in)   , optional :: etrans(kbdim)
  REAL(dp), INTENT(out)  , optional :: redevap(kbdim) ! reduction of evapotranspiration from SOILCHANGE
  REAL(dp), INTENT(in)   , optional :: vfc(kbdim)
  REAL(dp), INTENT(in)   , optional :: vpwp(kbdim)
  REAL(dp), INTENT(in)   , optional :: spsi(kbdim)
  INTEGER,  INTENT(inout), optional :: jllog      ! gridbox number for logoutput in SOILHYD
!
!  local variables
!
  INTEGER :: jl,jk
  REAL(dp) ::                                                             &
       zraind (kbdim)      &
     , zros(kbdim),          zdrain(kbdim)

  REAL(dp) ::                                                         &
       zorvari, zorvars, zdrmin, zdrmax, zdrexp,                      &
       zmprcp, zwlp, zwdtr, zwslim, zconw2, zconw3,                   &
       zroeff, zbws, zb1, zbm, zconw1, zlyeps, zvol, zprfl,           &
       zlysic, zwsup
  REAL(dp) :: zdtime

  ! Variables necessary for evaporation fluxes and soil layers
  REAL(dp), PARAMETER :: zeps  = 1.e-10_dp
  REAL(dp), PARAMETER :: zfak1 = 1._dp           ! =1., Dummy for SOILCHANGE (not used)
  INTEGER,  PARAMETER :: kb1   = 1               ! = 1, Digit for 1st field value for transfer to subroutine
  INTEGER,  PARAMETER :: ilog  = 0               ! Switch for logoutput in SOILHYD
!!!!  PARAMETER (jllog=918, ilog=3)        ! No logoutput for ilog=0

  REAL(dp) ::                                                         &
       dzrsi(kbdim, nsoil), dzsi(kbdim, nsoil),                       &
       zwsat(kbdim,nsoil),                        &
       zwpwp(kbdim,nsoil),                                            &
       wsim1m(kbdim, nsoil), aebsoil(kbdim), aetrans(kbdim)

  REAL(dp) :: zinfil(kbdim)    ! Infiltration for transfer to SOILCHANGE - was previously a local scalar 
  INTEGER  :: isch             ! Switch to control SOILCHANGE

!
!      Parameters
!
  zdtime  = delta_time
  zorvari = 100._dp
#ifndef STANDALONE
  zorvars = 1000._dp*64._dp/ngl
#else
  zorvars = 1000._dp*64._dp/nlat
  IF (nlat == 2) zorvars = 1000._dp*64._dp/96  !vg: for site level jsbach: T63 value
#endif
  zdrmin  = 0.001_dp/(3600._dp*1000._dp)
  zdrmax  = 0.1_dp/(3600._dp*1000._dp)
  zdrexp  = 1.5_dp

!     ------------------------------------------------------------------
!
!*    1.     Convert water fluxes to [m water equivalent * timestep]

  zros   = 0._dp
  zdrain = 0._dp
  zraind = prain * zdtime / rhoh2o
  zevwld = 0._dp
  DO jl=1,kbdim
    IF (lmask(jl)) THEN
!!$      zros(jl) = 0._dp
!!$      zdrain(jl) = 0._dp
!!$      zraind(jl) = prain(jl) * zdtime/rhoh2o
      zevwld(jl) = (1._dp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o
    END IF
  ENDDO

! Initialization of working arrays for multi-layer soil scheme
  IF (nsoil > 1) THEN
    CALL SoilDef(kb1, kbdim, kbdim, nsoil, lperma,       &
                 dzr, dzs, vpor, vfc, vpwp, lmask, zdel, &
                 dzrsi, dzsi, zwsat, zwsfc, zwpwp, ice)

!   *** Unit conversion for Transpiration [kg/m**2/s] --> [m]
    aetrans(1:kbdim) = etrans(1:kbdim)  * zdtime / rhoh2o
    redevap(1:kbdim) = 0._dp

  ENDIF   

  IF (lsurf) THEN

!     ------------------------------------------------------------------
!*  4.     Water budget
!
!*  4.1    Skin reservoir (vegetation and bare soil)
!
    DO jl=1,kbdim
      IF (lmask(jl) .AND. .NOT. lpglac(jl)) THEN
!
!*      4.1.1  Interception of rain
!
        zmprcp     = MIN(zraind(jl)*cvinter,pwlmx(jl)-pwl(jl))
        zwlp       = pwl(jl)    + zmprcp
        zraind(jl) = zraind(jl) - zmprcp
!
!*      4.1.2  Evaporation or dew collection
!
        pwl(jl)    = MIN(MAX(0._dp,zwlp+zevwld(jl)),pwlmx(jl))
        pevwsd(jl) = pevwsd(jl)-(pwl(jl)-zwlp)
!HAG
!       Get sure that ESKIN comprises only amounts that really change the skin reservoir
        zevwld(jl) = pwl(jl)-zwlp ! This is supposedly the correct way
!        if (nsoil==5) zevwld(jl) = pwl(jl)-zwlp ! Use this for compatibility with older jsbach revisions
      ELSE
        pwl(jl) = 0._dp
      END IF
    ENDDO

!   Calculation of Bare Soil Evaporation
!   Unit of pevwsd is [m]
    IF (nsoil > 1) aebsoil(1:kbdim) = pevwsd(1:kbdim) - aetrans(1:kbdim)
!
!*  4.2    Soil reservoir
!
    DO jl=1,kbdim

      zinfil(jl) = 0._dp
      IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
        zwslim = 0.05_dp*pwsmx(jl)

        zroeff = MAX(0._dp, porostd(jl)-zorvari)   &
                     /(porostd(jl)+zorvars)
        zbws   = MAX(MIN(zroeff,0.5_dp),0.01_dp)
        zb1    = 1._dp+zbws
        zbm    = 1._dp/zb1
        zconw1 = pwsmx(jl)*zb1
        zlyeps = 0._dp
        zvol   = 0._dp

!*    4.2.1  Surface runoff, infiltration and evaporation from soil

        ! Positive Enoskin_nosnow is added to throughfall if the old bucket is used
        IF (nsoil == 1) THEN

          IF (pevwsd(jl) >= 0.0_dp) THEN
            zprfl = pmlres(jl)+zraind(jl)+pevwsd(jl)
          ELSE
            pws(jl) = pws(jl)    + pevwsd(jl)
            zprfl   = pmlres(jl) + zraind(jl)
          END IF

        ELSE

!         ******* Positive EBSoil and Etrans are added to throughfall separately
!         ******* Negative fluxes change soil moisture in SOILCHANGE.
          zprfl = pmlres(jl)+zraind(jl)
          IF (aebsoil(jl) >= 0.0_dp) zprfl = zprfl + aebsoil(jl)
          IF (aetrans(jl) >= 0.0_dp) zprfl = zprfl + aetrans(jl)

        ENDIF

        IF (ptsoil(jl) < tmelt) THEN
          zros(jl) = zprfl
        ELSE
          IF (zprfl > 0._dp .AND. pws(jl) > zwslim) THEN
            IF (pws(jl) > pwsmx(jl)) THEN
              zlyeps = pws(jl)-pwsmx(jl)
            ELSE
              zlyeps = 0._dp
            END IF
            zlysic     = MIN(1._dp,(pws(jl)-zlyeps)/pwsmx(jl))
            zvol       = (1._dp-zlysic)**zbm-zprfl/zconw1
            zros(jl)   = zprfl-(pwsmx(jl)-pws(jl))
            IF (zvol > 0._dp) zros(jl) = zros(jl)+pwsmx(jl)*zvol**zb1
            zros(jl)   = MAX(zros(jl),0._dp)
            zinfil(jl) = zprfl-zros(jl)
          ELSE
            zros(jl)   = 0._dp
            zinfil(jl) = zprfl
          END IF
          ! In multi-layer scheme, infiltration changes soil moisture in SOILCHANGE.
          IF (nsoil == 1) pws(jl) = pws(jl)+zinfil(jl)
        END IF
!
!*    4.2.2  Drainage and total runoff
!
      ENDIF
    ENDDO
!
!   ****** Drainage if the old bucket is used
    IF (nsoil == 1) THEN

      DO jl=1,kbdim

        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN

!         *** Term shifting moved from the beginning of the original soil hydrology loop
          zwdtr  = 0.90_dp*pwsmx(jl)
          zwslim = 0.05_dp*pwsmx(jl)
          zconw2 = pwsmx(jl)-zwdtr
          zconw3 = zdrmax-zdrmin

          IF (pws(jl) <= zwslim) THEN
            zdrain(jl) = 0._dp
          ELSE
            IF (ptsoil(jl) > tmelt) THEN
              zdrain(jl) = zdrmin*MAX(pws(jl)/pwsmx(jl),0.7_dp)
              IF (pws(jl) > zwdtr) THEN
                 zdrain(jl) = zdrain(jl)+zconw3*((pws(jl)-zwdtr)/zconw2)**zdrexp
              END IF
              zdrain(jl) = zdrain(jl)*zdtime
              zdrain(jl) = MIN(zdrain(jl),pws(jl)-zwslim)
              pws(jl)    = pws(jl)-zdrain(jl)
            ELSE
              zdrain(jl) = 0._dp
            END IF
          END IF
          zwsup    = MAX(pws(jl)-pwsmx(jl),0._dp)
          pws(jl)  = pws(jl)-zwsup
          zros(jl) = zros(jl)+zdrain(jl)+zwsup
        ELSE
           pws(jl) = 0._dp
        END IF
      ENDDO

    ELSE ! nsoil == 1
      !  Change of soil moisture in multi-layer scheme
      !  change belongs to one time step = zdtime

      ! Soil moisture of previous time step: wsim1m
      wsim1m(:,:) = wsi(:,:)

      ! Zero initialization for drainage field before SOILCHANGE
      zdrain(1:kbdim) = 0._dp

      ! Soilchange expects aetrans and aebsoil in [m]
      isch = 1     ! Changes by evaporation fluxes
                   ! Soilchange only uses negative evaporation fluxes as positive fluxes,
                   ! e.g. by dew formation, are part of the throughfall zprfl.
      CALL SoilChange(kb1, kbdim, kbdim, nsoil,     &
            isch, zdtime, zfak1, lmask,             &
            wsi , wsim1m,                           &
            dzr, dzrsi,    dzsi, zwsat, zwsfc,      &
            zinfil, zdrain, aetrans, aebsoil,       &
            redevap)
      ! Unit of zinfil is [m] but soilchange expects m/s
      ! For numerical reasons, half of the water is infiltrated 
      ! before SOILHYD, half if it after
      zinfil(1:kbdim) = zinfil(1:kbdim) / zdtime / 2.

      isch = 2     ! Changes by infiltration

      CALL SoilChange(kb1, kbdim, kbdim, nsoil,     &
            isch, zdtime, zfak1, lmask,             &
            wsi , wsim1m,                           &
            dzr, dzrsi,    dzsi, zwsat, zwsfc,      &
            zinfil, zdrain, aetrans, aebsoil,       &
            redevap)

      ! Drainage calculation [m] with multi-layer soil scheme
      CALL SoilHyd(kb1, kbdim, kbdim, nsoil, lperma, &
           zdtime, ilog, jllog,                      &
           wsi, dzsi, zwsat, zwsfc,                  &
           fksat, vpor, bclapp, fmpot,               &     
           zdrain, spsi, zwpwp, zdel, ice            &
           )

      ! Second half of the water is infiltrated after SOILHYD
      isch = 2     ! Changes by infiltration
      CALL SoilChange(kb1, kbdim, kbdim, nsoil,     &
            isch, zdtime, zfak1, lmask,             &
            wsi , wsim1m,                           &
            dzr, dzrsi,    dzsi, zwsat, zwsfc,      &
            zinfil, zdrain, aetrans, aebsoil,       &
            redevap)

      ! Does overflow occurs due to infiltration in the layers? --> added to drainage
      DO jk=1, nsoil
        DO jl=kb1, kbdim
          IF (dzsi(jl,1) > 0) THEN
            IF (wsi(jl,jk) > zwsfc(jl,jk)) THEN
              zdrain(jl) = zdrain(jl) + wsi(jl,jk) - zwsfc(jl,jk)
              wsi(jl,jk) = zwsfc(jl,jk)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      ! WS must be re-derived from wsi after SOILHYD
      pws(1:kbdim) = 0.0_dp
      DO jk=1, nsoil
        DO jl=kb1, kbdim
          IF (lmask(jl)) THEN
            IF (dzrsi(jl,jk) >= zdel(jk)) THEN
              pws(jl) = pws(jl) + wsi(jl,jk)
            ELSE IF (dzrsi(jl,jk) > 0._dp) THEN
              pws(jl) = pws(jl) + wsi(jl,jk) * dzrsi(jl,jk)/dzsi(jl,jk)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
          zros(jl) = zros(jl)+zdrain(jl)
        ENDIF
      ENDDO

    ENDIF ! end of nsoil if clause.
!
!*  4.2.3  Runoff and drainage [m]
!
    DO jl=1,kbdim
      IF (lmask(jl)) THEN
        pros_hd(jl)   = zros(jl)-zdrain(jl)
        pdrain_hd(jl) = zdrain(jl)
      END IF
    END DO
!
!     ------------------------------------------------------------------
  END IF

END SUBROUTINE update_surf_up

END MODULE
