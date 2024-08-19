!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! mo_radiation_forcing: determine instantaneous radiative forcing by
!!                       calling the radiation calculation twice with different
!!                       input parameters
!!
!! @author Sebastian Rast, MPI Met, Hamburg, based on code by M.A.Thomas,
!!  S.J.Lorenz, and Ph.Stier
!!
!! $ID: n/a$
!!
!! @par Revision History
!! original source by J.S.Rast (2010-03-15) based on code by M.A.Thomas,
!!  S.J.Lorenz, and Ph.Stier
!!
MODULE mo_radiation_forcing
  USE mo_kind,                ONLY: wp
  USE mo_linked_list,         ONLY: t_stream
  USE mo_physical_constants,  ONLY: stbo,secperday
  USE mo_radiation_parameters,ONLY: psctm, lradforcing, cemiss
  USE mo_time_control,        ONLY: zdtime=>delta_time

  PUBLIC :: construct_forcing, prepare_forcing, calculate_forcing

  PRIVATE

  TYPE (t_stream), POINTER :: forcing

  !--- Auxiliary radiative fluxes:

  REAL(wp),  POINTER :: emter_for(:,:,:)
  REAL(wp),  POINTER :: emtef_for(:,:,:)
  REAL(wp),  POINTER :: trsol_for(:,:,:)
  REAL(wp),  POINTER :: trsof_for(:,:,:)

  !--- Forcing fields TOA, SUR:

  REAL(wp),  POINTER :: fsw_clear_top(:,:) !instantaneous sw forcing clear sky top of atmosphere
  REAL(wp),  POINTER :: fsw_total_top(:,:) !instantaneous sw forcing all sky top of atmosphere
  REAL(wp),  POINTER :: fsw_clear_sur(:,:) !instantaneous sw forcing clear sky surface
  REAL(wp),  POINTER :: fsw_total_sur(:,:) !instantaneous sw forcing all sky surface

  REAL(wp),  POINTER :: flw_clear_top(:,:) !instantaneous lw forcing clear sky top of atmosphere
  REAL(wp),  POINTER :: flw_total_top(:,:) !instantaneous lw forcing all sky top of atmosphere
  REAL(wp),  POINTER :: flw_clear_sur(:,:) !instantaneous lw forcing clear sky surface
  REAL(wp),  POINTER :: flw_total_sur(:,:) !instantaneous lw forcing all sky surface

  !--- Radiation flux forcing:
  REAL(wp),  POINTER :: d_aflx_sw(:,:,:)  !3d instantaneous sw forcing all sky
  REAL(wp),  POINTER :: d_aflx_lw(:,:,:)  !3d instantaneous lw forcing all sky
  REAL(wp),  POINTER :: d_aflx_swc(:,:,:) !3d instantaneous sw forcing clear sky
  REAL(wp),  POINTER :: d_aflx_lwc(:,:,:) !3d instantaneous lw forcing clear sky

  !--- Heating rate forcing:
  REAL(wp), POINTER  :: netht_lw(:,:,:)   !3d forcing of net lw heating rate (K/day) 
  REAL(wp), POINTER  :: netht_sw(:,:,:)   !3d forcing of net sw heating rate (K/day)

CONTAINS

  !-----------------------------------------------------------------------------
  !>
  !! construct_forcing: set up memory for forcing calculation
  !!
  !! @par Revision History
  !! original source by J.S.Rast (2010-04-15) based on code by M.A.Thomas,
  !!  S.J.Lorenz, and Ph.Stier
  !-----------------------------------------------------------------------------
  SUBROUTINE construct_forcing
    USE mo_memory_base,           ONLY: new_stream, add_stream_element, &
         default_stream_setting, &
         add_stream_reference, &
         default_output
    USE mo_linked_list,           ONLY: SURFACE, HYBRID, HYBRID_H

    IMPLICIT NONE

    LOGICAL :: lpost, lrerun


    !--- Create new stream:

    lrerun=lradforcing(1) .OR. lradforcing(2)
    lpost= lrerun .AND. default_output

    CALL new_stream (forcing ,'forcing', lpost=lpost, lrerun=lrerun, lcontnorest=.TRUE.)

    !--- 1) Add standard fields for post-processing:

    CALL add_stream_reference (forcing, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (forcing, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (forcing, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (forcing, 'gboxarea','geoloc',lpost=.TRUE.)

    !--- 2) Add auxiliary stream elements:

    CALL default_stream_setting (forcing, lrerun   = .TRUE.,   &
         laccu    = .FALSE.,  &
         lpost    = .FALSE.    )

    CALL add_stream_element(forcing,  'emter_for',  emter_for, leveltype=HYBRID_H, units='W m-2')
    CALL add_stream_element(forcing,  'emtef_for',  emtef_for, leveltype=HYBRID_H, units='W m-2')
    CALL add_stream_element(forcing,  'trsol_for',  trsol_for, leveltype=HYBRID_H, units='W m-2')
    CALL add_stream_element(forcing,  'trsof_for',  trsof_for, leveltype=HYBRID_H, units='W m-2')

    !--- 3) Add forcing stream elements:

    CALL default_stream_setting (forcing, lrerun    = .TRUE.,   &
         contnorest= .TRUE. ,  &
         laccu     = .TRUE.,   &
         lpost     = .TRUE.,   &
         table     = 200       )

    IF (lradforcing(1)) THEN

      !--- Forcing fields TOA, SUR:
      CALL default_stream_setting (forcing, leveltype = SURFACE )
      CALL add_stream_element(forcing,  'FSW_CLEAR_TOP',  fsw_clear_top,     units='W m-2',             &
           longname='Mean instantaneous sw forcing clear sky top of the atmosphere', &
           code=11                                                                   )
      CALL add_stream_element(forcing,  'FSW_TOTAL_TOP',  fsw_total_top,     units='W m-2',             &
           longname='Mean instantaneous sw forcing all sky top of the atmosphere',   &
           code=12                                                                   )
      CALL add_stream_element(forcing,  'FSW_CLEAR_SUR',  fsw_clear_sur,     units='W m-2',             &
           longname='Mean instantaneous sw forcing clear sky surface',               &
           code=13                                                                   )
      CALL add_stream_element(forcing,  'FSW_TOTAL_SUR',  fsw_total_sur,     units='W m-2',             &
           longname='Mean instantaneous sw forcing all sky surface',                 &
           code=14                  )
      !--- Radiation flux forcing:
      CALL default_stream_setting (forcing, leveltype = HYBRID_H)
      CALL add_stream_element (forcing, 'd_aflx_sw', d_aflx_sw, units='W m-2',                          &
           longname='Accumulated SW flux anomalies - all sky ',                     &
           code=15                                                                  )
      CALL add_stream_element (forcing, 'd_aflx_swc', d_aflx_swc, units='W m-2',                        &
           longname='Accumulated SW flux anomalies - clear',                        &
           code=16                                                                  )
      !--- Heating rate forcing:
      CALL default_stream_setting (forcing, leveltype = HYBRID)
      CALL add_stream_element (forcing, 'netht_sw', netht_sw, units='K/d',                              &
           longname='Net SW heating rate forcing',                                  &
           code=17                                                                  )

    END IF

    IF (lradforcing(2)) THEN

      !--- Forcing fields TOA, SUR:
      CALL default_stream_setting (forcing, leveltype = SURFACE )
      CALL add_stream_element(forcing,  'FLW_CLEAR_TOP',  flw_clear_top,     units='W m-2',             &
           longname='Mean instantaneous lw forcing clear sky top of the atmosphere', &
           code=21                                                                   )
      CALL add_stream_element(forcing,  'FLW_TOTAL_TOP',  flw_total_top,     units='W m-2',             &
           longname='Mean instantaneous lw forcing all sky top of the atmosphere',   &
           code=22                                                                   )
      CALL add_stream_element(forcing,  'FLW_CLEAR_SUR',  flw_clear_sur,     units='W m-2',             &
           longname='Mean instantaneous lw aerosol forcing clear sky surface',       &
           code=23                                                                   )
      CALL add_stream_element(forcing,  'FLW_TOTAL_SUR',  flw_total_sur,     units='W m-2',             &
           longname='Mean instantaneous lw forcing all sky surface',                 &
           code=24                                                                   )
      !--- Radiation flux forcing:
      CALL default_stream_setting (forcing, leveltype = HYBRID_H)
      CALL add_stream_element (forcing, 'd_aflx_lw', d_aflx_lw, units='W m-2',                          &
           longname='Accumulated LW flux anomalies - all sky ',                     &
           code=25                                                                  )
      CALL add_stream_element (forcing, 'd_aflx_lwc', d_aflx_lwc, units='W m-2',                        &
           longname='Accumulated LW flux anomalies - clear',                        &
           code=26                                                                  )
      !--- Heating rate forcing:
      CALL default_stream_setting (forcing, leveltype = HYBRID)
      CALL add_stream_element (forcing, 'netht_lw', netht_lw, units='K/d',                              &
           longname='Net LW heating rate forcing',                                  &
           code=27                                                                  )

    END IF
  END SUBROUTINE construct_forcing
  !-----------------------------------------------------------------------------
  !>
  !! prepare_forcing: calculate intermediate quantities for forcing calculation
  !!
  !! @par Revision History
  !! original source by J.S.Rast (2010-04-16) based on code by M.A.Thomas and
  !!  S.J.Lorenz
  !-----------------------------------------------------------------------------
  SUBROUTINE prepare_forcing(                                                   &
       & kproma           ,kbdim             ,klevp1           ,krow           ,&
       & pflx_dnlw        ,pflx_dnsw         ,pflx_dnlw_clr    ,pflx_dnsw_clr   )

    INTEGER, INTENT(in)     :: kproma, kbdim, klevp1, krow
    REAL(wp), INTENT(in)    :: &
         pflx_dnlw(kbdim,klevp1),    & !< Net dwnwrd LW flux [Wm2]
         pflx_dnsw(kbdim,klevp1),    & !< Net dwnwrd SW flux [Wm2] for forcing
         pflx_dnlw_clr(kbdim,klevp1),& !< Net dn LW flux (clear sky) [Wm2]
         pflx_dnsw_clr(kbdim,klevp1)   !< Net dn SW flux (clear sky) [Wm2]

    emter_for(1:kproma,1:klevp1,krow)=pflx_dnlw(1:kproma,1:klevp1)
    emtef_for(1:kproma,1:klevp1,krow)=pflx_dnlw_clr(1:kproma,1:klevp1)
    trsol_for(1:kproma,1:klevp1,krow)=pflx_dnsw(1:kproma,1:klevp1)
    trsof_for(1:kproma,1:klevp1,krow)=pflx_dnsw_clr(1:kproma,1:klevp1)

  END SUBROUTINE prepare_forcing
  !-----------------------------------------------------------------------------
  !>
  !! calculate_forcing: calculate radiative forcing 
  !!
  !! @par Revision History
  !! original source by J.S.Rast (2010-04-16) based on code by M.A.Thomas and
  !!  S.J.Lorenz
  !-----------------------------------------------------------------------------
  SUBROUTINE calculate_forcing( &
       &  kproma             ,kbdim               ,klevp1          &
       & ,krow               ,pi0                 ,pconvfact       &
       & ,pflxs              ,pflxs0              ,pflxt           &
       & ,pflxt0             ,pti                 ,pztsnew         )
    INTEGER, INTENT(in)          :: kproma, kbdim, klevp1, krow
    REAL(wp), INTENT(in)         :: &
         & pi0(kbdim),                 & !> solar irradiance at top of atmosphere
         & pconvfact(kbdim,klevp1-1),  & !> conversion factor radiation flux -> heating rate
         & pflxs(kbdim,klevp1),        & !> short wave net radiation flux per Watt irradiance (all sky)
         & pflxs0(kbdim,klevp1),       & !> short wave net radiation flux per Watt irradiance (clear sky)
         & pflxt(kbdim,klevp1),        & !> long wave net radiation flux per Watt irradiance (all sky)
         & pflxt0(kbdim,klevp1),       & !> long wave net radiation flux per Watt irradiance (clear sky)
         & pti(kbdim,klevp1),          & !> temperature at interfaces (atmosphere)
         & pztsnew(kbdim)                !> surface temperature

    INTEGER                     :: jk, jl, klev
    REAL(wp)                    :: zflxs_all_for(kbdim,klevp1), &
         zflxs_clear_for(kbdim,klevp1), &
         zflxt_all_for(kbdim,klevp1), &
         zflxt_clear_for(kbdim,klevp1)
    REAL(wp)                    :: zdtdt_sw, zdtdt_all_for_sw, &
         zdtdt_lw, zdtdt_all_for_lw

    klev=klevp1-1

    IF (lradforcing(1)) THEN
      !--- Radiation flux forcing:
      zflxs_all_for(1:kproma,1:klevp1)=SPREAD(pi0(1:kproma),2,klevp1)* &
           trsol_for(1:kproma,1:klevp1,krow)
      zflxs_clear_for(1:kproma,1:klevp1)=SPREAD(pi0(1:kproma),2,klevp1)* &
           trsof_for(1:kproma,1:klevp1,krow)
      !  forcing solar wave length bands:
      DO jk = 1, klev
        DO jl = 1, kproma
          zdtdt_sw=pconvfact(jl,jk)*(pflxs(jl,jk+1)-pflxs(jl,jk))
          zdtdt_all_for_sw=pconvfact(jl,jk)*(zflxs_all_for(jl,jk+1)-zflxs_all_for(jl,jk))
          netht_sw(jl,jk,krow) = netht_sw(jl,jk,krow)+secperday*(zdtdt_sw-zdtdt_all_for_sw)*zdtime
        ENDDO
      END DO
      DO jk = 1, klevp1
        d_aflx_sw(1:kproma,jk,krow) = d_aflx_sw(1:kproma,jk,krow) + &
             (pflxs(1:kproma,jk)-zflxs_all_for(1:kproma,jk))*zdtime
        d_aflx_swc(1:kproma,jk,krow) = d_aflx_swc(1:kproma,jk,krow) + &
             (pflxs0(1:kproma,jk)-zflxs_clear_for(1:kproma,jk))*zdtime
      END DO
      fsw_total_top(1:kproma,krow)=d_aflx_sw(1:kproma,1,krow)
      fsw_total_sur(1:kproma,krow)=d_aflx_sw(1:kproma,klevp1,krow)
      fsw_clear_top(1:kproma,krow)=d_aflx_swc(1:kproma,1,krow)
      fsw_clear_sur(1:kproma,krow)=d_aflx_swc(1:kproma,klevp1,krow)
    END IF

    IF (lradforcing(2)) THEN
      !--- Radiation flux forcing:
      zflxt_all_for(1:kproma,1:klev)=emter_for(1:kproma,1:klev,krow)
      ! in the following formulae, the order of calculation has to be as is
      ! if not, for an aerosol free atmosphere the forcing is not exactly 0
      ! because of numeric effects (see the corresponding formulae in radheat)
      zflxt_all_for(1:kproma,klevp1)=emter_for(1:kproma,klevp1,krow)+ &
           cemiss*stbo*pti(1:kproma,klevp1)**4
      zflxt_all_for(1:kproma,klevp1)=zflxt_all_for(1:kproma,klevp1)- &
           cemiss*stbo*pztsnew(1:kproma)**4
      zflxt_clear_for(1:kproma,1:klev)=emtef_for(1:kproma,1:klev,krow)
      zflxt_clear_for(1:kproma,klevp1)=emtef_for(1:kproma,klevp1,krow) + &
           cemiss*stbo*pti(1:kproma,klevp1)**4
      zflxt_clear_for(1:kproma,klevp1)=zflxt_clear_for(1:kproma,klevp1)- &
           cemiss*stbo*pztsnew(1:kproma)**4
      DO jk = 1,klev
        DO jl = 1, kproma
          zdtdt_lw=pconvfact(jl,jk)*(pflxt(jl,jk+1)-pflxt(jl,jk))
          zdtdt_all_for_lw=pconvfact(jl,jk)*(zflxt_all_for(jl,jk+1)-zflxt_all_for(jl,jk))
          netht_lw(jl,jk,krow) = netht_lw(jl,jk,krow)+secperday*(zdtdt_lw-zdtdt_all_for_lw)*zdtime
        END DO
      END DO
      DO jk = 1, klevp1
        d_aflx_lw(1:kproma,jk,krow) = d_aflx_lw(1:kproma,jk,krow) + &
             (pflxt(1:kproma,jk)-zflxt_all_for(1:kproma,jk))*zdtime
        d_aflx_lwc(1:kproma,jk,krow) = d_aflx_lwc(1:kproma,jk,krow) + &
             (pflxt0(1:kproma,jk)-zflxt_clear_for(1:kproma,jk))*zdtime
      END DO
      flw_total_top(1:kproma,krow) = d_aflx_lw(1:kproma,1,krow)
      flw_total_sur(1:kproma,krow) = d_aflx_lw(1:kproma,klevp1,krow)
      flw_clear_top(1:kproma,krow) = d_aflx_lwc(1:kproma,1,krow)
      flw_clear_sur(1:kproma,krow) = d_aflx_lwc(1:kproma,klevp1,krow)
    END IF
  END SUBROUTINE calculate_forcing
END MODULE mo_radiation_forcing
