!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! This module contains substantial contributions by the Max Planck Institute for Biogeochemistry, Jena.
!!
!! ---------------------------------------------------------------------------------------------
!!
!! This module belongs to the standalone version of JSBACH
!!
!+ Definition for offline driving secundary parameters (other than forcing [= weather])
!
Module mo_jsbalone

  !
  ! Description:
  !   Derives additional parameters for offline driving of jsbach
  !   like wind profile, stability params, Richtmeyr-Morton Coeffs
  !   ect.
  ! 
  !   updates to version 0.2 by SZ:
  !   - revised such that drag and zchl to conform to echam code
  !   - assuming explicit coupling to atmosphere for the calculation 
  !     of the Richtmeyr-Morton Coefficients
  !   - calculation of reference height should be based on vegetation
  !     height from lctlib, and not (as now) from roughness length
  !   obsolete:
  !   code needs cleaning from obsolete variables
  !   Some additional parameters use values from the last time step
  !   
  !
  ! Current Code Owner:  jsbach_admin
  ! 
  ! History:
  !
  ! Version       Date                Comment
  ! -------       ------              -------
  ! 0.1           2005/03/22          Org. Code admin
  ! 0.2           2011/06/03          updates to code to calculate drag and exchange 
  !                                   coefficients as coupled to Echam, SZ
  !
  ! Modules used:
  !
  USE mo_kind,              ONLY: dp
  USE mo_jsbach,            ONLY: options_type
  USE mo_jsbach_grid,       ONLY: grid_type, domain_type
  USE mo_linked_list,       ONLY: t_stream
  USE mo_atmosphere,        ONLY: sat_specific_humidity
 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Drive_type, init_driving, update_driving
  !
  !-----------------------------------------------------------------------------------------------------------------
  !
  TYPE Drive_type
     REAL(dp), POINTER, Dimension(:)  :: & !! Dimension (nland)
          Air_Temp_old , &                !! Air Temperature from last time step
          geopot, &
          wind, &
          temp_air, &
          qair, &
          eair, &
          precip_rain, &
          precip_snow, &
          co2, &
          lwdown, &
          vis_net, &
          nir_net, &
          czenith, &
          pressure, &
          etacoef, &
          eqacoef, &
          etbcoef, &
          eqbcoef, &
          cair, &
          csat, &
          cdrag, &
          z0h, &
          z0m, &
          albedo_vis, & 
          albedo_nir, &
          zchl, &
          evap_act, &
          evap_pot, &
          evap_act2pot, &
          soil_temp, &
          veg_height, &
! DIAG TYPE : upper elements have to be cleared !!
          geopot_acc, &
          wind_acc, &
          temp_air_acc, &
          qair_acc, &
          precip_rain_acc, &
          precip_snow_acc, &
          co2_acc, &
          lwdown_acc, &
          vis_net_acc, &
          nir_net_acc, &
          czenith_acc, &
          pressure_acc, &
          sensible_heat_acc, &
          evap_act_sum, &
          soil_temp_acc, &
          evap_pot_sum, &
          cdrag_acc, &
          ustar_acc, &
          zchl_acc, &
          latent_acc, &
          surf_runoff_hd, &
          drainage_hd, &
          glac_runoff_evap
  END TYPE Drive_type


  TYPE(t_stream),  POINTER, SAVE :: IO_driving      !! Memory stream for driving variables in offline mode

  !
  !-----------------------------------------------------------------------------------------------------------------
  !
  CONTAINS
  !
  !=================================================================================================================
  !
  SUBROUTINE driving_init_memory(g_nland, l_nland, Drv_Var, stream)

    USE mo_jsbach,        ONLY : missing_value
    USE mo_memory_base, ONLY : add =>add_stream_element
    USE mo_netCDF,      ONLY : max_dim_name
    
    INTEGER         , INTENT(in)             :: g_nland, l_nland
    TYPE(drive_type), INTENT(inout)          :: Drv_Var
    TYPE(t_stream)  , POINTER, OPTIONAL      :: stream

    INTEGER                     :: dim1p(1), dim1(1)
    CHARACTER(LEN=max_dim_name) :: dim1n(1)


    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    CALL add(stream, 'air_temp_old', Drv_Var%Air_Temp_old,dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'geopot'      , Drv_Var%geopot, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'wind'        , Drv_Var%wind, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'temp_air'    , Drv_Var%temp_air, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'qair'        , Drv_Var%qair, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'eair'        , Drv_Var%eair, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'precip_rain' , Drv_Var%precip_rain, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'precip_snow' , Drv_Var%precip_snow, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'co2'         , Drv_Var%co2, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'lwdown'      , Drv_Var%lwdown, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'vis_net'     , Drv_Var%vis_net, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'nir_net'     , Drv_Var%nir_net, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'czenith'     , Drv_Var%czenith, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'pressure'    , Drv_Var%pressure, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'etacoef'     , Drv_Var%etacoef, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'etbcoef'     , Drv_Var%etbcoef, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'eqacoef'     , Drv_Var%eqacoef, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'eqbcoef'     , Drv_Var%eqbcoef, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'cair'        , Drv_Var%cair , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'csat'        , Drv_Var%csat , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'cdrag'       , Drv_Var%cdrag , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'zchl'        , Drv_Var%zchl  , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'z0h'         , Drv_Var%z0h  , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'z0m'         , Drv_Var%z0m  , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'albedo_vis'  , Drv_Var%albedo_vis, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'albedo_nir'  , Drv_Var%albedo_nir, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'evap_act'    , Drv_Var%evap_act , dim1p, dim1, dimnames=dim1n, code=92, laccu=.false., lpost=.false., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(stream, 'evap_pot'    , Drv_Var%evap_pot , dim1p, dim1, dimnames=dim1n, code=93, laccu=.false., lpost=.false., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(stream, 'evap_act2pot' , Drv_Var%evap_act2pot , dim1p, dim1, dimnames=dim1n, code=94, laccu=.false., lpost=.true., &
         lmiss=.TRUE., missval=missing_value)
    CALL add(stream, 'soil_temp'   , Drv_Var%soil_temp , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'veg_height'   , Drv_Var%veg_height , dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'surf_runoff_hd', Drv_Var%surf_runoff_hd, dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'drainage_hd',    Drv_Var%drainage_hd,    dim1p, dim1, dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
    CALL add(stream, 'glac_runoff_evap', Drv_Var%glac_runoff_evap, dim1p,dim1,dimnames=dim1n, code=1, laccu=.false., lpost=.false.)
!
!   DIAG TYPE: Upper elements have to be cleaned ! 
!   
    CALL add(stream, 'geopot_acc'      , Drv_Var%geopot_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'wind_acc'        , Drv_Var%wind_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'temp_air_acc'    , Drv_Var%temp_air_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'qair_acc'        , Drv_Var%qair_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'precip_rain_acc' , Drv_Var%precip_rain_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'precip_snow_acc' , Drv_Var%precip_snow_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'co2_acc'         , Drv_Var%co2_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'lwdown_acc'      , Drv_Var%lwdown_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'vis_net_acc'     , Drv_Var%vis_net_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'nir_net_acc'     , Drv_Var%nir_net_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'czenith_acc'     , Drv_Var%czenith_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'pressure_acc'    , Drv_Var%pressure_acc, dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'evap_act_sum'    , Drv_Var%evap_act_sum, dim1p, dim1, dimnames=dim1n, code=95, laccu=.FALSE., &
       lmiss=.TRUE., missval=missing_value, lpost=.FALSE., contnorest=.TRUE.)
    CALL add(stream, 'evap_pot_sum'    , Drv_Var%evap_pot_sum, dim1p, dim1, dimnames=dim1n, code=96, laccu=.FALSE., &
       lmiss=.TRUE., missval=missing_value, lpost=.FALSE., contnorest=.TRUE.)
    CALL add(stream, 'cdrag_acc'       , Drv_Var%cdrag_acc , dim1p, dim1, dimnames=dim1n, code=1, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)
    CALL add(stream, 'ustar_acc'       , Drv_Var%ustar_acc , dim1p, dim1, dimnames=dim1n, code=200, laccu=.true., &
       lmiss=.TRUE., missval=missing_value)
    CALL add(stream, 'zchl_acc'        , Drv_Var%zchl_acc  , dim1p, dim1, dimnames=dim1n, code=2, laccu=.true., &
       lmiss=.TRUE., missval=missing_value, lpost=.false.)

  END SUBROUTINE driving_init_memory
  !
  !=================================================================================================================
  !
  SUBROUTINE init_driving(grid, domain, global_options, drv_var)         

    USE mo_linked_list, ONLY : LAND
    USE mo_memory_base, ONLY : new_stream, default_stream_setting

    TYPE(grid_type)   ,    INTENT(in)   :: grid
    TYPE(domain_type) ,    INTENT(in)   :: domain
    TYPE(options_type),    INTENT(in)   :: global_options
    TYPE(drive_type)  ,    INTENT(inout):: Drv_Var

    CALL new_stream(IO_driving, 'driving', filetype=global_options%FileType, ztype=global_options%FileZtype, &
                    lpost=.true., lrerun=.TRUE.)

    CALL default_stream_setting(IO_driving, lpost=.true., lrerun=.TRUE., repr=LAND)

    call driving_init_memory(grid%nland, domain%nland, Drv_Var, stream=IO_driving)

    Drv_Var%evap_act = 0._dp
    Drv_Var%evap_pot = 0._dp

    Drv_Var%geopot_acc = 0._dp
    Drv_Var%wind_acc = 0._dp
    Drv_Var%temp_air_acc = 0._dp
    Drv_Var%qair_acc = 0._dp
    Drv_Var%precip_rain_acc = 0._dp
    Drv_Var%precip_snow_acc = 0._dp
    Drv_Var%co2_acc = 0._dp
    Drv_Var%lwdown_acc = 0._dp
    Drv_Var%vis_net_acc = 0._dp
    Drv_Var%nir_net_acc = 0._dp
    Drv_Var%czenith_acc = 0._dp
    Drv_Var%pressure_acc = 0._dp
    Drv_Var%evap_act_sum = 0._dp
    Drv_Var%evap_pot_sum = 0._dp
    Drv_Var%cdrag_acc = 0._dp
    Drv_Var%ustar_acc = 0._dp
    Drv_Var%zchl_acc = 0._dp
    Drv_Var%evap_act2pot = 1._dp
    Drv_Var%z0m = 1._dp
    Drv_Var%z0h = 1._dp
    Drv_Var%albedo_vis = 0.1_dp
    Drv_Var%albedo_nir = 0.3_dp
    Drv_Var%cair = 0.5_dp
    Drv_Var%csat = 0.5_dp
    Drv_Var%veg_height = 10._dp
    Drv_Var%soil_temp = 273._dp
    Drv_Var%surf_runoff_hd = 0._dp
    Drv_Var%drainage_hd = 0._dp
    Drv_Var%glac_runoff_evap = 0._dp

  END SUBROUTINE init_driving
  !
  !=================================================================================================================
  !
  SUBROUTINE update_driving( &
       kdim, &                            
       kland, &                          
       drv_var, &                     
       HeightWind, HeightHumidity, &
       czenith, &
       wind, temp_air, qair,  &
       precip_rain, precip_snow, &
       lwdown, &
       rad_uv_down, &
       rad_par_down, &                       
       rad_nir_down, &                     
       pressure, &
       CO2_concentration)

    USE mo_time_control,     ONLY: delta_time, lstart, time_step_len
    USE mo_jsbach_constants, ONLY: SpecificHeatDryAirConstPressure, GasConstantDryAir, Gravity, vtmpc1
    USE mo_jsbach,           ONLY: new_day
    USE mo_physc2,           ONLY: cb, cc, cvdifts, ckap

    INTEGER,           INTENT(in)    :: kdim                          !! Length of vectors (if call from ECHAM5 this
    INTEGER,           INTENT(in)    :: kland                         !! Number of land points in vectors
    TYPE(drive_type),  INTENT(inout) :: drv_var                       !! driving variables
    REAL(dp),              INTENT(in)    :: HeightWind                !! Defines lowest layer height-> where measurements are taken
    REAL(dp),              INTENT(in)    :: HeightHumidity            !! Defines lowest layer height-> where measurements are taken
    REAL(dp),              INTENT(in)    :: czenith(kdim)             !! Cosine of solar zenith angle
    REAL(dp),              INTENT(in)    :: wind(kdim)                !! Lowest level wind speed [m/s]
    REAL(dp),              INTENT(in)    :: temp_air(kdim)            !! Lowest level air temperature [Kelvin]
    REAL(dp),              INTENT(in)    :: qair(kdim)                !! Lowest level specific humidity
    REAL(dp),              INTENT(in)    :: precip_rain(kdim)         !! Precipitation as rain [kg/(m^2 s)]
    REAL(dp),              INTENT(in)    :: precip_snow(kdim)         !! Precipitation as snow [kg/(m^2 s)]
    REAL(dp),              INTENT(in)    :: rad_uv_down(kdim)         !! UV reaching the surface [W/m^2]
    REAL(dp),              INTENT(in)    :: rad_par_down(kdim)        !! PAR reaching the surface [W/m^2]
    REAL(dp),              INTENT(in)    :: rad_nir_down(kdim)        !! solar radiation in the near infrared band reaching surface
                                                                      !!   [W/m^2]
    REAL(dp),              INTENT(in)    :: lwdown(kdim)              !! Downward longwave flux
    REAL(dp),              INTENT(in)    :: pressure(kdim)            !! Surface pressure
    REAL(dp),              INTENT(in)    :: CO2_concentration(kdim)   !! Atmospheric CO2 concentration [kg(CO2)/kg(air)]

    REAL(dp), DIMENSION(kland) :: zg, zchnl, zcfnchl, zcfhl, zdu2, zucfhl, zscfl, ztvd, ztvs, zril, zcons
    REAL(dp), DIMENSION(kland) :: qsat_surf, zcdnl, zcfml, zcfncl, zucfl, zcdn2m, zcdnr, zcfm2m, zustl, zustarl
    REAL(dp)                   :: zepz0o, zcons8, zcons9, zcons11, zcons12, zkappa

    !-------------------------------------------------------------------------------------------------
    ! Some constants from echam/mo_surface_land.f90
    !-------------------------------------------------------------------------------------------------
    zepz0o        = 2._dp
    zcons8        = 2._dp * cb
    zcons9        = 3._dp * cb
    zcons11       = 3._dp * cb * cc
    zcons12       = cvdifts * time_step_len * Gravity / GasConstantDryAir
    zkappa        = GasConstantDryAir / SpecificHeatDryAirConstPressure

    !------------------------------------------------------------------------------------
    ! Approximation of cdrag
    !------------------------------------------------------------------------------------
    ! squared wind shear ! minimum wind speed square from echam/mo_surface_land.f90:precalc_land
    zdu2(:) = MAX(wind(:)**2,1._dp)  
   
    ! virtual potential air temperature (see mo_surface_boundary.f90)
    ! according to Saucier, WJ Principles of Meteoroligical Analyses
    ! tv = t * (1 + 0.61 * q * t)    ! virtual temperature
    ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
    ! tvd = tair * (100000/p_pa)^zkappa * 1 + vtmpc1 * qair) ! virtual potential temperature
    ztvd(:) = ( temp_air(:) * ( 100000._dp/pressure(:))**zkappa ) * & 
              ( 1._dp + vtmpc1 * qair(:) )
    ! virtual potential surface temperature 
    qsat_surf(:) = sat_specific_humidity(Drv_var%soil_temp(:),pressure(:))
    ztvs = Drv_var%soil_temp(:) * ( 100000._dp/pressure(:))**zkappa * &
           ( 1._dp + vtmpc1 * ( Drv_var%csat(:) * qsat_surf(:) + ( 1._dp - Drv_var%cair(:) ) * qair(:)))

    ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
    ! If HeightWind is set, then the measurement height + an offset from the average vegetation height
    ! is used. Otherwise, the code defaults to the half-level of ECHAM's lowest atmospheric layer
    IF(HeightWind > 0._dp) THEN
       zg(:) = (0.75_dp * Drv_Var%veg_height + HeightWind) * Gravity
    ELSE
      zg(:) = ztvd(:) *  GasConstantDryAir * LOG(1._dp / ( 1._dp - 0.0025_dp ))
    ENDIF

    ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
    ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
    zril(:) = zg(:) * ( ztvd(:) - ztvs(:) ) / ( zdu2(:) * ztvd(:) )
    zril(:) = MAX(MIN(zril(:),5._dp),-5._dp)

    ! Neutral drag coefficient for momentum and heat
    zcdnl(:) = (ckap / LOG(1._dp + zg(:) / (Gravity * Drv_Var%z0m(:) )))**2
    zchnl(:) = (ckap / LOG(1._dp + zg(:) / (Gravity * Drv_Var%z0h(:) )))**2

    ! account for stable/unstable case: helper variables
    zscfl(:) = SQRT (  1._dp + ABS(zril(:)))
    zucfl(:) = 1._dp / (1._dp + zcons11 * zcdnl(:) * SQRT(ABS(zril(:)) * (1._dp  &
            + zg(:) / (Gravity * Drv_Var%z0m(:)))))
    zucfhl(:) = 1._dp / (1._dp + zcons11 * zchnl(:) * SQRT(ABS(zril(:)) * (1._dp  &
            + zg(:) / (Gravity * Drv_Var%z0h(:)))))

    ! ignoring cloud water correction (see mo_surface_land.f90)
    zcons(:) = zcons12 * pressure(:) / ( temp_air(:) * (1._dp + vtmpc1 * qair(:)))
    zcfncl(:)  = zcons(:) * SQRT(zdu2(:)) * zcdnl(:)
    zcfnchl(:)  = zcons(:) * SQRT(zdu2(:)) * zchnl(:)

    ! Stable / Unstable case
    WHERE ( zril(:) .GT. 0._dp )
       zcfml(:) = zcfncl(:) / (1._dp + zcons8 * zril(:) * zscfl(:))
       zcfhl(:) = zcfnchl(:) / (1._dp + zcons8 * zril(:) * zscfl(:))
    ELSEWHERE
       zcfml(:) = zcfncl(:) * (1._dp - zcons8 * zril(:) * zucfl(:))
       zcfhl(:) = zcfnchl(:) * (1._dp - zcons9 * zril(:) * zucfhl(:))
    ENDWHERE

    Drv_Var%cdrag(1:kland) = zcfhl(:)
    Drv_Var%zchl(1:kland)  = zcfhl(:) / zcfnchl(:) * zchnl(:)

    ! Computation of the PBL extension -> ustar
    WHERE(Drv_Var%z0m(:) .GT. zepz0o)
       zcdn2m(:)= ( ckap / LOG(1._dp + zg(:) / (Gravity * zepz0o)))**2
    ELSEWHERE
       zcdn2m(:)= zcdnl(:)
    END WHERE
    zcdnr(:)    = zcdn2m(:) / zcdnl(:)
    ! stable and unstable case
    WHERE(Drv_Var%z0m(:) .GT. zepz0o .AND. zril(:) .LT. 0._dp)
          zcfm2m(:)= zcfncl(:) * zcdnr(:) * (1._dp - zcons8 * zril(:) &
               / (1._dp + zcons11 * zcdn2m(:) * SQRT(ABS(zril(:)) &
               * (1._dp + zg(:) / (Gravity * zepz0o)))))
    ELSEWHERE
       zcfm2m(:)= zcfml(:) * zcdnr(:)
    END WHERE
    zustl(:)    = zcfm2m(:) * SQRT(zdu2(:))
    ! again omitting the cloud water correction term zx
    zustarl(:)  = SQRT(zustl(:) * temp_air(:) &
            * (1._dp + vtmpc1 * qair(:) ) &
            / (zcons12 * pressure(:)))

    !---------------------------------------------------------------------------------------------------------------
    ! Computation of Richtmeyr-Morton Coefficients
    ! This follows now Jan Polcher's explicit solution, i.e. atmospheric conditions at t+1 are assumed to be valid
    !---------------------------------------------------------------------------------------------------------------    
    Drv_Var%etacoef = 0.0_dp
    Drv_Var%etbcoef(1:kland) = SpecificHeatDryAirConstPressure * temp_air + HeightHumidity * Gravity 
    Drv_Var%eqacoef = 0.0_dp
    Drv_Var%eqbcoef(1:kland) = qair(1:kland)

    !-------------------------------------------------------------------------------------------------
    ! compute net radiation from downward radiation using albedo of previous time step 
    !-------------------------------------------------------------------------------------------------
    Drv_Var%vis_net = (rad_uv_down(:) + rad_par_down(:)) * (1._dp - Drv_Var%albedo_vis(:))
    Drv_Var%nir_net = rad_nir_down * (1._dp - Drv_Var%albedo_nir(:))

    !-------------------------------------------------------------------------------------------------
    ! Write fields to output stream
    !-------------------------------------------------------------------------------------------------
    Drv_Var%geopot(1:kland) = zg(:)
    Drv_Var%wind = wind
    Drv_Var%temp_air = temp_air
    Drv_Var%qair = qair 
    Drv_Var%precip_rain = precip_rain 
    Drv_Var%precip_snow = precip_snow 
    Drv_Var%co2 = CO2_concentration
    Drv_Var%lwdown = lwdown 
    Drv_Var%czenith = czenith
    Drv_Var%pressure = pressure

    Drv_Var%geopot_acc(1:kland) = Drv_Var%geopot_acc(1:kland) + delta_time*zg(:)
    Drv_Var%wind_acc = Drv_Var%wind_acc + delta_time*wind
    Drv_Var%temp_air_acc = Drv_Var%temp_air_acc + delta_time*temp_air
    Drv_Var%qair_acc = Drv_Var%qair_acc + delta_time*qair 
    Drv_Var%precip_rain_acc = Drv_Var%precip_rain_acc + delta_time*precip_rain 
    Drv_Var%precip_snow_acc = Drv_Var%precip_snow_acc + delta_time*precip_snow 
    Drv_Var%co2_acc = Drv_Var%co2_acc + delta_time*CO2_concentration
    Drv_Var%lwdown_acc = Drv_Var%lwdown_acc + delta_time*lwdown 
    Drv_Var%vis_net_acc= Drv_Var%vis_net_acc + delta_time*Drv_Var%vis_net
    Drv_Var%nir_net_acc = Drv_Var%nir_net_acc + delta_time*Drv_Var%nir_net
    Drv_Var%czenith_acc = Drv_Var%czenith_acc + delta_time*czenith
    Drv_Var%pressure_acc = Drv_Var%pressure_acc + delta_time*pressure
    Drv_Var%cdrag_acc(1:kland) = Drv_Var%cdrag_acc + delta_time*Drv_Var%cdrag(1:kland)
    Drv_Var%zchl_acc(1:kland) = Drv_Var%zchl_acc + delta_time*Drv_Var%zchl(1:kland)
    Drv_Var%ustar_acc(1:kland) = Drv_Var%ustar_acc + delta_time*zustarl(1:kland)

    ! SZ stuff needed for QAIR FORCING = NONE
    Drv_Var%evap_act_sum(1:kland) = Drv_Var%evap_act_sum(1:kland) + MIN(Drv_Var%evap_act(1:kland),-EPSILON(1._dp))
    Drv_Var%evap_pot_sum(1:kland) = Drv_Var%evap_pot_sum(1:kland) + MIN(Drv_Var%evap_pot(1:kland),-EPSILON(1._dp))
    IF (new_day .AND. .NOT. lstart ) THEN
       Drv_Var%evap_act2pot(1:kland) = Drv_Var%evap_act_sum / Drv_Var%evap_pot_sum
       Drv_Var%evap_act_sum(1:kland) = 0._dp
       Drv_Var%evap_pot_sum(1:kland) = 0._dp
    ENDIF

  END SUBROUTINE update_driving
  !
  !=================================================================================================================
  !
End Module mo_jsbalone
