!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE update_surf_down ( kbdim    &
     , pqm1, ptsl, pspeed10, ptm1      &
     , pwl, pwlmx                      &
     , psn, pcvs, psnc,       pgld     &
     , pgrndcapc                       &
     , pevapl,     pevapot, pevwsd     &
     , prain, psnow                    &
     , lmask, lpglac                   &
     , palac                           &
     , psnaccum, psnmel, pglmel, progl &
     , pmlres, tte_corr                &
     ) 
  !
  !     ------------------------------------------------------------------
  ! ptsl          Sfc temp [K] at timestep t+dt (unfiltered)
  ! pspeed10      Wind speed at 10m height [m/s] ... from 'vdiff'
  ! ptm1          Air temp [K] at timestep t-dt (filtered)
  ! pwl           Water content [m] in skin reservoir
  !               (vegetation and bare soil)
  ! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
  !               of vegetation index and leaf area index)
  ! psn           Snow depth [m water equivalent] at the ground
  ! pcvs          Fractional snow cover (function of psn in *physc*)
  ! psnc          Snow depth [m water equivalent] at the canopy
  ! pgld          Glacier depth (including snow) [m water equivalent]
  ! pgrndcapc     Heat capacity of the uppermost soil layer [j/m**2/K]
  ! pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
  ! pevapot       Potential evaporation/sublimation [kg/m**2/s]
  ! pevwsd        Evapotranspiration without sublimation and evaporation from interception reservoir [kg/m**2]
  ! prain         Totalrainfall [kg/m**2/s]
  ! psnow         Total snowfall [kg/m**2/s]
  ! lpglac        Logical glacier mask
  ! palac         Precipitation minus sublimation at glacier points [m water equivalent]
  ! psnaccum      Snow accumulation (budget) at non-glacier points [m water equivalent]
  ! psnmel        Snow/ice melt at land points [m water equivalent]
  ! pglmel        Snow/ice melt at glacier points [m water equivalent]
  ! progl         Glacier runoff (rain+snow/ice melt, but no calving) [m water equivalent]
  !
  ! The following local variables represent the respective fluxes
  ! integrated over one timestep (delta_time) and divided by the density of
  ! water (rhoh2o). Units are m water equivalent.
  !
  ! zraind        Total rain
  ! zsnowd        Total snowfall
  ! zevttd        Total evaporation
  ! zevsnd        Sublimation from snow
  ! zsncmelt      Snow melt in the canopy
  ! pmlres        Residual melt water available for infiltration into the
  !               non-frozen soil after being intercepted by the
  !               skin reservoir
  ! tte_corr      ???
  !
  !       Rest folgt spaeter ....
  !
  !     *SURF* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.
  !              CALCULATE FLUXES OF TOTAL RAIN, TOTAL SNOW AND EVAPO-
  !              RATION FROM THE THREE RESERVOIRS (SN, WS, WL)
  !              CONVERT FLUXES (KG/M**2*S) TO CHANGES OF WATER LEVELS (M)
  !              DURING TIMESTEP DELTA_TIME.
  !
  !     J.F.GELEYN     E.C.M.W.F.     08/06/82.
  !     MODIFIED BY
  !     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
  !     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
  !     J.-P. SCHULZ   MPI - 1997 : IMPLEMENTATION OF IMPLICIT
  !                                 COUPLING BETWEEN LAND SURFACE
  !                                 AND ATMOSPHERE.
  !     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
  !     MODIFIED BY M. ESCH        MPI - APR  1999
  !     MODIFIED BY E. ROECKNER    MPI - JAN  2001
  !     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
  !     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir 
  !                                                for snow changed
  !     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
  !
  !     MODIFICATION
  !
  !     PURPOSE
  !
  !     INTERFACE.
  !
  !          *SURF* IS CALLED FROM *PHYSC*.
  !
  !     METHOD.
  !
  !     EXTERNALS.
  !
  !          NONE.
  !
  !     REFERENCE.
  !
  !          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
  !     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
  !
  USE mo_param_switches,   ONLY: lsurf
  USE mo_physc2,           ONLY: cwlmax
  USE mo_jsbach_constants, ONLY: cvinter, tmelt, SpecificHeatDryAirConstPressure, gravity, rhoh2o, vtmpc2, LatentHeatFusion
  USE mo_time_control,   ONLY : lstart, delta_time
  USE mo_kind, ONLY: dp
  !
  IMPLICIT NONE

  INTEGER,  INTENT(in) :: kbdim
  LOGICAL,  INTENT(in) :: lmask(kbdim)
  LOGICAL,  INTENT(in) :: lpglac(kbdim)
  REAL(dp), INTENT(in) :: pqm1(kbdim)
  REAL(dp), INTENT(in) :: pspeed10(kbdim)
  REAL(dp), INTENT(in) :: ptm1(kbdim)
  REAL(dp), INTENT(in) :: pwlmx(kbdim)
  REAL(dp), INTENT(in) :: pcvs(kbdim)
  REAL(dp), INTENT(in) :: pgrndcapc(kbdim)
  REAL(dp), INTENT(in) :: pevapl(kbdim)
  REAL(dp), INTENT(in) :: pevapot(kbdim)
  REAL(dp), INTENT(in) :: prain(kbdim)
  REAL(dp), INTENT(in) :: psnow(kbdim)

  REAL(dp), INTENT(inout) :: ptsl(kbdim)
  REAL(dp), INTENT(inout) :: pwl(kbdim)
  REAL(dp), INTENT(inout) :: psn(kbdim)
  REAL(dp), INTENT(inout) :: psnc(kbdim)
  REAL(dp), INTENT(inout) :: pgld(kbdim)
  REAL(dp), INTENT(out)   :: psnaccum(kbdim)
  REAL(dp), INTENT(out)   :: psnmel(kbdim)
  REAL(dp), INTENT(out)   :: pglmel(kbdim)
  REAL(dp), INTENT(out)   :: progl(kbdim)

  REAL(dp), INTENT(out)   :: pevwsd(kbdim)
  REAL(dp), INTENT(out)   :: palac(kbdim)
  REAL(dp), INTENT(inout) :: tte_corr(kbdim)
  REAL(dp), INTENT(out)   :: pmlres(kbdim)
  !
  !  local arrays
  !
  INTEGER :: jl
  REAL(dp) ::                                                             &
       zraind (kbdim),       zsnowd(kbdim),         zevttd(kbdim)        &
       , zevsnd(kbdim),      zsncmelt(kbdim)     
  REAL(dp) :: zlfdcp(kbdim)
  !
  !  local scalars
  !
  REAL(dp) ::                                                             &
       zsmelt, zsnmlt,      &
       zmprcp,              &
       zc2, zc3, zsncp, zexpt, zexpw, zsncmax, zsncwind
  REAL(dp) :: zdtime, zrcp, zsncfac
  !
  !      Parameters
  !
  zdtime = delta_time
  zc2=1.87E5_dp
  zc3=1.56E5_dp
  zsncfac=rhoh2o*gravity/zdtime
  !
  !     ------------------------------------------------------------------
  !
  !*    1.     Convert water fluxes to [m water equivalent * timestep]
  !
  DO 110 jl=1,kbdim
     palac(jl)=0._dp
     progl(jl)=0._dp
     pmlres(jl)=0._dp
     psnaccum(jl)=0._dp
     psnmel(jl)=0._dp
     pglmel(jl)=0._dp
     zsncmelt(jl)=0._dp
     zraind(jl)=prain(jl)             *zdtime/rhoh2o
     zsnowd(jl)=psnow(jl)             *zdtime/rhoh2o
     IF(lmask(jl)) THEN
        zrcp=1._dp/(SpecificHeatDryAirConstPressure*(1._dp+vtmpc2*MAX(0.0_dp,pqm1(jl))))
        zlfdcp(jl)=LatentHeatFusion*zrcp
        zevttd(jl)=pevapl(jl)                        *zdtime/rhoh2o
        zevsnd(jl)=pcvs(jl)*pevapot(jl)              *zdtime/rhoh2o
        pevwsd(jl)=zevttd(jl)-zevsnd(jl)
     ELSE
        pevwsd(jl) = 0._dp ! Sain initialization for diagnostics
     END IF
110 END DO
  !
  IF(lsurf) THEN
     !
     !     ------------------------------------------------------------------
     !
     !*    2.     Budgets of snow (canopy and ground) and glaciers
     !
     !*    2.1    Snow changes in the canopy (interception of snowfall,
     !            sublimation, melting, unloading due to wind)
     !
     DO 210 jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
           psnaccum(jl)=zsnowd(jl)+zevsnd(jl)
           zsncmax=MAX(0.0_dp,pwlmx(jl)-cwlmax)
           zmprcp=MIN(zsnowd(jl)*cvinter,zsncmax-psnc(jl))
           zsncp=psnc(jl)+zmprcp
           zsnowd(jl)=zsnowd(jl)-zmprcp
           psnc(jl)=MIN(MAX(0._dp,zsncp+zevsnd(jl)),zsncmax)
           zevsnd(jl)=zevsnd(jl)-(psnc(jl)-zsncp)
           zexpt=MAX(0._dp,ptm1(jl)+3._dp-tmelt)*zdtime/zc2
           zexpw=pspeed10(jl)*zdtime/zc3
           zsncmelt(jl)=psnc(jl)*(1._dp-EXP(-zexpt))
           psnc(jl)=psnc(jl)-zsncmelt(jl)
           zsncwind=psnc(jl)*(1._dp-EXP(-zexpw))
           psnc(jl)=psnc(jl)-zsncwind
           zsnowd(jl)=zsnowd(jl)+zsncwind
           tte_corr(jl)=zsncmelt(jl)*zsncfac*zlfdcp(jl)
           !
           !   pwl(jl)=pwl(jl)+zsncmelt(jl) see section 2.5
           !
        ELSE
           psnc(jl)=0._dp
        END IF
210  END DO
     !
     !*    2.2    Snowfall and sublimation on land (excluding glaciers)
     !
     DO 220 jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
           psn(jl)=psn(jl)+zsnowd(jl)+zevsnd(jl)
           IF (psn(jl).LT.0._dp) THEN
              pevwsd(jl)=pevwsd(jl)+psn(jl)
              psn(jl)=0._dp
           END IF
        ELSE
           psn(jl)=0._dp
        END IF
220  END DO
     !
     !*    2.3    Snowfall and sublimation on glaciers and diagnostics
     !
     DO 230 jl=1,kbdim
        IF (lpglac(jl)) THEN
           pgld(jl)=pgld(jl)+zsnowd(jl)+zevsnd(jl)
           palac(jl)=zraind(jl)+zsnowd(jl)+zevttd(jl)
           progl(jl)=zraind(jl)
        END IF
230  END DO
     !
     !*    2.4    Snow and glacier melt
     !
     IF (.NOT. lstart) THEN
        DO 240 jl=1,kbdim
           IF (lmask(jl)) THEN
              IF (ptsl(jl).GT.tmelt) THEN
                 IF (lpglac(jl)) THEN
                    pglmel(jl)=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(LatentHeatFusion*rhoh2o)
                    pgld(jl)=pgld(jl)-pglmel(jl)
                    progl(jl)=progl(jl)+pglmel(jl)
                    ptsl(jl)=tmelt
                 ELSE IF (psn(jl).GT.0._dp) THEN
                    zsmelt=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(LatentHeatFusion*rhoh2o)
                    psnmel(jl)=MIN(psn(jl),zsmelt)
                    ptsl(jl)=ptsl(jl)-psnmel(jl)*LatentHeatFusion*rhoh2o/pgrndcapc(jl)
                    psn(jl)=MAX(psn(jl)-psnmel(jl),0._dp)
                 END IF
              END IF
           END IF
240     END DO
     END IF
     !
     !*    2.5    Snow budget and meltwater (glacier-free land only)
     !
     DO 250 jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
           pwl(jl)=pwl(jl)+zsncmelt(jl)                  ! Add melt water on canopy to skin reservoir
           zsnmlt=psnmel(jl)+MAX(0._dp,pwl(jl)-pwlmx(jl))   ! Excess water on canopy drips to ground and
           pwl(jl)=MIN(pwlmx(jl),pwl(jl))                ! skin reservoir on canopy is set to maximum
           pmlres(jl)=zsnmlt
!!$           psnmel(jl)=psnmel(jl)+zsncmelt(jl)
           psnaccum(jl)=psnaccum(jl)-psnmel(jl)-zsncmelt(jl)
        END IF
250  END DO
     !
     ! water fluxes
     ! Note: conversion from m water equivalent to kg/m^2s - multiply by rhoh2o / zdtime
!!$      psnmel(1:kbdim) = psnmel(1:kbdim) * rhoh2o / zdtime
!!$      psnaccum(1:kbdim) = psnaccum(1:kbdim) * rhoh2o / zdtime
!!$      progl(1:kbdim) = progl(1:kbdim) * rhoh2o / zdtime

  END IF
  !
  RETURN

END SUBROUTINE update_surf_down
