!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!!++mgs: see lines 210 and 211
MODULE mo_surface_memory

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_decomposition,   ONLY: ldc=>local_decomposition
  USE mo_memory_base, ONLY: delete_stream, add => add_stream_element, &
                            default_stream_setting
  USE mo_surface_types, ONLY: land2atmos_type, ocean2atmos_type, ice2atmos_type, surface_type, grid_box_mean_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_surface     ! routine to construct the g3b table
  PUBLIC :: destruct_surface      ! routine to destruct  the g3b table
  PUBLIC :: surf                  ! the g3b table
  PUBLIC :: land, ocean, ice, box, surface

  TYPE(land2atmos_type) :: land
  TYPE(ocean2atmos_type) :: ocean
  TYPE(ice2atmos_type) :: ice
  TYPE(grid_box_mean_type) :: box
  TYPE(surface_type) :: surface

  ! declaration of predefined fields within this module 

  ! JSBACH variables for the interface call
  REAL(dp), POINTER, PUBLIC, DIMENSION(:,:) :: &
       jsswnir, jsswdifnir, jsswvis, jsswdifvis, jsswpar, jsswdifpar
  REAL(dp), POINTER, PUBLIC, DIMENSION(:,:) :: &
       jsswniracc, jsswvisacc, jsswparacc, jsswdifniracc, jsswdifvisacc, jsswdifparacc
  REAL(dp), POINTER, PUBLIC, DIMENSION(:,:) :: &
       jrsfl, jrsfc, jssfl, jssfc, ztrfli, zsofli

 TYPE (t_stream), POINTER     :: surf

CONTAINS
!-----------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE construct_surface

    INTEGER :: nproma, ngpblks

    nproma  = ldc%nproma
    ngpblks = ldc%ngpblks

    ALLOCATE( surface%is_land(nproma,ngpblks)  )
    ALLOCATE( surface%is_seaice(nproma,ngpblks)  )
    ALLOCATE( surface%is_ocean(nproma,ngpblks)   )
    ALLOCATE( surface%is_ocean_for_coup(nproma,ngpblks)   )
    ALLOCATE( surface%is_ice_for_coup(nproma,ngpblks)   )


    ! set default attributes for the surface stream
    CALL default_stream_setting (surf         &
                                ,lrerun=.TRUE.     &
                                ,lpost=.FALSE.      &
                                ,table=181 ,bits=16)

    ! Add fields to the surface stream.

    ! JSBACH variables for the interface call

    ! radiation

    CALL add (surf,'jsswnir'              ,jsswnir                     ,code=128,laccu=.FALSE.,lpost=.FALSE. , &
              longname='net surface NIR'                       ,units='W/m**2')
    CALL add (surf,'jsswdifnir'           ,jsswdifnir                  ,code=129,laccu=.FALSE.,lpost=.FALSE. , &
              longname='fraction of diffuse NIR'               ,units='')
    CALL add (surf,'jsswvis'              ,jsswvis                     ,code=130,laccu=.FALSE.,lpost=.FALSE. , &
              longname='net surface visible'                   ,units='W/m**2')
    CALL add (surf,'jsswdifvis'           ,jsswdifvis                  ,code=131,laccu=.FALSE.,lpost=.FALSE. , &
              longname='fraction of diffuse visible'           ,units='')
    CALL add (surf,'jsswpar'              ,jsswpar                     ,code=240,laccu=.FALSE.,lpost=.FALSE. , &
              longname='downward surface PAR'                  ,units='W/m**2')
    CALL add (surf,'jsswdifpar'           ,jsswdifpar                  ,code=241,laccu=.FALSE.,lpost=.FALSE. , &
              longname='fraction of diffuse PAR'               ,units='')
    CALL add (surf,'jsswniracc'           ,jsswniracc                  ,code=203,laccu=.TRUE. ,lpost=.TRUE. , &
              longname='net surface NIR accumulated'           ,units='W/m**2')
    CALL add (surf,'jsswvisacc'           ,jsswvisacc                  ,code=204,laccu=.TRUE. ,lpost=.TRUE. , &
              longname='net surface visible accumulated'       ,units='W/m**2')
    CALL add (surf,'jsswparacc'           ,jsswparacc                  ,code=244,laccu=.TRUE. ,lpost=.TRUE. , &
              longname='downward surface PAR accumulated'      ,units='W/m**2')
    CALL add (surf,'jsswdifniracc'        ,jsswdifniracc               ,code=205,laccu=.TRUE. ,lpost=.TRUE. , &
             longname='net surface diffuse NIR accumulated'    ,units='W/m**2')
    CALL add (surf,'jsswdifvisacc'        ,jsswdifvisacc               ,code=206,laccu=.TRUE. ,lpost=.TRUE. , &
              longname='net surface diffuse visible accumulated',units='W/m**2')
    CALL add (surf,'jsswdifparacc'        ,jsswdifparacc               ,code=246,laccu=.TRUE. ,lpost=.TRUE. , &
              longname='downward surface diffuse PAR accumulated',units='W/m**2')
    CALL add (surf,'longwave_down_acc'    ,box%longwave_down_acc       ,code=207,laccu=.TRUE. ,lpost=.TRUE. , &
             longname='downward longwave radiation'            ,units='W/m**2') 
    CALL add (surf,'ztrfli'               ,ztrfli                      ,code=137,laccu=.FALSE.,lpost=.FALSE., &
             longname='longwave flux ice'                      ,units='W/m**2')
    CALL add (surf,'zsofli'               ,zsofli                      ,code=138,laccu=.FALSE.,lpost=.FALSE., &
             longname='solar flux ice'                         ,units='W/m**2')

    ! precipitation

    CALL add (surf,'jrsfl'                ,jrsfl                       ,code=98 ,laccu=.FALSE.,lpost=.FALSE., &
             longname='large scale rain'                       ,units='kg/m**2s')
    CALL add (surf,'jrsfc'                ,jrsfc                       ,code=99 ,laccu=.FALSE.,lpost=.FALSE., &
             longname='conv scale rain'                        ,units='kg/m**2s')
    CALL add (surf,'jssfl'                ,jssfl                       ,code=100,laccu=.FALSE.,lpost=.FALSE., &
             longname='large scale snow'                       ,units='kg/m**2s')
    CALL add (surf,'jssfc'                ,jssfc                       ,code=101,laccu=.FALSE.,lpost=.FALSE., &
             longname='conv scale snow'                        ,units='kg/m**2s')

    ! visible and NIR albedo

    CALL add (surf,'alsoi_vis'            ,ice%albedo_vis              ,code=183,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of ice visible'                 ,units='')
    CALL add (surf,'alsoi_vis_dir'        ,ice%albedo_vis_dir                   ,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of ice visible direct'          ,units='')
    CALL add (surf,'alsoi_vis_dif'        ,ice%albedo_vis_dif                   ,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of ice visible diffuse'         ,units='')
    CALL add (surf,'alsoi_nir'            ,ice%albedo_nir              ,code=184,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of ice nir'                     ,units='')
    CALL add (surf,'alsoi_nir_dir'        ,ice%albedo_nir_dir                   ,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of ice nir direct'              ,units='')
    CALL add (surf,'alsoi_nir_dif'        ,ice%albedo_nir_dif                   ,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of ice nir diffuse'             ,units='')
    CALL add (surf,'alsow_vis'            ,ocean%albedo_vis            ,code=185,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of water visible'               ,units='')
    CALL add (surf,'alsow_vis_dir'        ,ocean%albedo_vis_dir                 ,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of water visible direct'        ,units='')
    CALL add (surf,'alsow_vis_dif'        ,ocean%albedo_vis_dif                 ,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of water visible diffuse'       ,units='')
    CALL add (surf,'alsow_nir'            ,ocean%albedo_nir            ,code=186,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of water nir'                   ,units='')
    CALL add (surf,'alsow_nir_dir'        ,ocean%albedo_nir_dir                 ,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of water nir direct'            ,units='')
    CALL add (surf,'alsow_nir_dif'        ,ocean%albedo_nir_dif                 ,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of water nir diffuse'           ,units='')
    CALL add (surf,'alsol_vis'            ,land%albedo_vis             ,code=187,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of land visible'                ,units='')
    CALL add (surf,'alsol_nir'            ,land%albedo_nir             ,code=188,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of land NIR'                    ,units='')
    CALL add (surf,'albedo_mean_vis'      ,box%albedo_vis              ,code=189,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface albedo visible mean'           ,units='')
    CALL add (surf,'albedo_mean_vis_dir'  ,box%albedo_vis_dir                   ,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface albedo visible mean'           ,units='')
    CALL add (surf,'albedo_mean_vis_dif'  ,box%albedo_vis_dif                   ,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface albedo visible mean'           ,units='')
    CALL add (surf,'albedo_mean_nir'      ,box%albedo_nir              ,code=190,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface albedo nir mean'               ,units='')
    CALL add (surf,'albedo_mean_nir_dir'  ,box%albedo_nir_dir                   ,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface albedo nir mean'               ,units='')
    CALL add (surf,'albedo_mean_nir_dif'  ,box%albedo_nir_dif                   ,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface albedo nir mean'               ,units='')
    
    ! temporary vars for mo_surface_X

    CALL add (surf,'atm_tot_cloud_water'  ,box%atm_tot_cloud_water     ,code=1  ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'atm_sat_spec_hum '    ,box%atm_sat_spec_hum        ,code=2  ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'atm_pot_temp'         ,box%atm_pot_temp            ,code=3  ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'atm_vir_pot_temp'     ,box%atm_vir_pot_temp        ,code=4  ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'atm_lat_heat_fact'    ,box%atm_lat_heat_fact       ,code=5  ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'atm_dry_stat_energy'  ,box%atm_dry_stat_energy     ,code=6  ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'atm_liq_wat_pot_temp' ,box%atm_liq_wat_pot_temp    ,code=7  ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'evaporation_inst'     ,box%evaporation_inst        ,code=8  ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zetnl'                ,land%zetnl                  ,code=9  ,laccu=.FALSE.,lpost=.FALSE.) !RM COEFF LAND
    CALL add (surf,'zeqnl'                ,land%zeqnl                  ,code=10 ,laccu=.FALSE.,lpost=.FALSE.) ! -"-
    CALL add (surf,'zftnl'                ,land%zftnl                  ,code=11 ,laccu=.FALSE.,lpost=.FALSE.) ! -"-
    CALL add (surf,'zfqnl'                ,land%zfqnl                  ,code=12 ,laccu=.FALSE.,lpost=.FALSE.) ! -"-
    CALL add (surf,'zetnw'                ,ocean%zetnw                 ,code=13 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zftnw'                ,ocean%zftnw                 ,code=14 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zeqnw'                ,ocean%zeqnw                 ,code=15 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zfqnw'                ,ocean%zfqnw                 ,code=16 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zetni'                ,ice%zetni                   ,code=17 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zftni'                ,ice%zftni                   ,code=18 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zeqni'                ,ice%zeqni                   ,code=19 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zfqni'                ,ice%zfqni                   ,code=20 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zqsw'                 ,ocean%zqsw                  ,code=21 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcptw'                ,ocean%zcptw                 ,code=22 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zriw'                 ,ocean%zriw                  ,code=23 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcfhw'                ,ocean%zcfhw                 ,code=24 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zchw'                 ,ocean%zchw                  ,code=25 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbnw'                 ,ocean%zbnw                  ,code=26 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbmw'                 ,ocean%zbmw                  ,code=27 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbhw'                 ,ocean%zbhw                  ,code=28 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zustarw'              ,ocean%zustarw               ,code=29 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ztkevw'               ,ocean%ztkevw                ,code=30 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcfmw'                ,ocean%zcfmw                 ,code=31 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'trfll'                ,land%trfll                  ,code=32 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'sofll'                ,land%sofll                  ,code=33 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'trflw'                ,ocean%trflw                 ,code=34 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'soflw'                ,ocean%soflw                 ,code=35 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'trfli'                ,ice%trfli                   ,code=36 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'sofli'                ,ice%sofli                   ,code=37 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ztklevl'              ,land%ztklevl                ,code=38 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zqklevl'              ,land%zqklevl                ,code=39 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'evaporation_inst_land',land%evaporation_inst       ,code=40 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'evaporation_pot_land' ,land%evaporation_pot        ,code=41 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'sens_heat_inst_land'  ,land%sensible_heat_flux_inst,code=42 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'latent_heat_inst_land',land%latent_heat_flux_inst  ,code=43 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'radiative_land_temp'  ,land%surface_temperature_rad,code=44 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'surface_temp_new'     ,land%surface_temperature_new,code=45 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'dry_static_energy_new',land%dry_static_energy_new  ,code=46 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'dewpoint_2_meter_land',land%dewpoint_2_meter       ,code=47 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcair_old'            ,land%zcair_old              ,code=48 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'land_fract'           ,box%land_fract              ,code=49 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ocean_fract'          ,box%ocean_fract             ,code=50 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'surface_qsat'         ,land%surface_qsat           ,code=51 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zdqsl'                ,land%zdqsl                  ,code=52 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zril'                 ,land%zril                   ,code=53 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zqsl'                 ,land%zqsl                   ,code=54 ,laccu=.FALSE.,lpost=.FALSE.)
!TS April 2010: CALL add (surf,'zcfncl'               ,land%zcfncl                 ,code=55 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcfnc'                ,box%zcfnc                   ,code=55 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zchl'                 ,land%zchl                   ,code=56 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcfhl'                ,land%zcfhl                  ,code=57 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbnl'                 ,land%zbnl                   ,code=58 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbhnl'                ,land%zbhnl                  ,code=59 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbml'                 ,land%zbml                   ,code=60 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbhl'                 ,land%zbhl                   ,code=61 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zustarl'              ,land%zustarl                ,code=62 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcfml'                ,land%zcfml                  ,code=63 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ztkevl'               ,land%ztkevl                 ,code=64 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zustl'                ,land%zustl                  ,code=65 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zspeedl'              ,land%wind_10_meter          ,code=66 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcair'                ,land%zcair                  ,code=67 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcsat'                ,land%zcsat                  ,code=68 ,laccu=.FALSE.,lpost=.FALSE.)
!    CALL add (surf,'z0h'                  ,land%z0h                    ,code=69 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcptl'                ,land%zcptl                  ,code=70 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcpq'                 ,land%zcpq                   ,code=71 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'surface_qsat_new'     ,land%surface_qsat_new       ,code=72 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'land_u_wind_10_meter' ,land%u_wind_10_meter        ,code=73 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'land_v_wind_10_meter' ,land%v_wind_10_meter        ,code=74 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'temp_2_meter'         ,land%temp_2_meter           ,code=75 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'momentum_exch_coeff'  ,box%momentum_ex_coef        ,code=76 ,laccu=.FALSE.,lpost=.FALSE.) 
    CALL add (surf,'zustw'                ,ocean%zustw                 ,code=78 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcvsi'                ,ice%zcvsi                   ,code=79 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'fluxres'              ,ocean%fluxres               ,code=80 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'sens_heat_flux_inst'  ,box%sensible_heat_flux_inst ,code=81 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'latent_heat_flux_inst',box%latent_heat_flux_inst   ,code=82 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'surf_vir_temp'        ,box%surf_vir_temp           ,code=83 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'surface_humidity'     ,box%surface_humidity        ,code=84 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ztklevw'              ,ocean%ztklevw               ,code=85 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zqklevw'              ,ocean%zqklevw               ,code=86 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'pevapw'               ,ocean%evaporation_inst      ,code=87 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'pahfsw'               ,ocean%sensible_heat_flux_inst,code=88 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'pahflw'               ,ocean%latent_heat_flux_inst ,code=89 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zt2w'                 ,ocean%temp_2_meter          ,code=90 ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'heat_exch_coeff'      ,box%heat_ex_coef            ,code=98 ,laccu=.FALSE.,lpost=.FALSE.) 
    CALL add (surf,'dewpoint_2_meter_oc'  ,ocean%dewpoint_2_meter      ,code=220,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ocean_u_wind_10_meter',ocean%u_wind_10_meter       ,code=221,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ocean_v_wind_10_meter',ocean%v_wind_10_meter       ,code=222,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ustar_mean'           ,box%ustar                   ,code=223,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'mean_tke_ztkevn'      ,box%ztkevn                  ,code=224,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'u_stress_mean'        ,box%u_stress                ,code=225,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'v_stress_mean'        ,box%v_stress                ,code=226,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'wind_10_meter_mean'   ,box%wind_10_meter           ,code=227,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zqsi'                 ,ice%zqsi                    ,code=228,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcpti'                ,ice%zcpti                   ,code=229,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zrii'                 ,ice%zrii                    ,code=230,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcfmi'                ,ice%zcfmi                   ,code=231,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zchi'                 ,ice%zchi                    ,code=232,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zcfhi'                ,ice%zcfhi                   ,code=233,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbni'                 ,ice%zbni                    ,code=234,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbmi'                 ,ice%zbmi                    ,code=235,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zbhi'                 ,ice%zbhi                    ,code=236,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zustari'              ,ice%zustari                 ,code=237,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ztkevi'               ,ice%ztkevi                  ,code=238,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zusti'                ,ice%zusti                   ,code=239,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ice_depth'            ,ice%ice_depth               ,code=240,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'SNOW_ON_ICE'          ,ice%SNOW_ON_ICE             ,code=241,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'surface_temperature'  ,box%surface_temperature     ,code=242,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'pcvi'                 ,box%seaice_fract            ,code=243,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'snow_cover_fract'     ,ice%snow_cover_fract        ,code=244,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ztklevi'              ,ice%ztklevi                 ,code=245,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zqklevi'              ,ice%zqklevi                 ,code=246,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'pevapi'               ,ice%evaporation_inst        ,code=247,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'pahfsi'               ,ice%sensible_heat_flux_inst ,code=248,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'pahfli'               ,ice%latent_heat_flux_inst   ,code=249,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'zt2i'                 ,ice%temp_2_meter            ,code=250,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'dewpoint_2_meter_ice' ,ice%dewpoint_2_meter        ,code=251,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'wind_10_meter_ice'    ,ice%wind_10_meter           ,code=252,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ice_u_wind_10_meter'  ,ice%u_wind_10_meter         ,code=253,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ice_v_wind_10_meter'  ,ice%v_wind_10_meter         ,code=254,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ameltdepth'           ,ice%ameltdepth                       ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'ameltfrac'            ,ice%ameltfrac                        ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'alsobs'               ,ice%alsobs                           ,laccu=.FALSE.,lpost=.FALSE.)
    CALL add (surf,'alsom'                ,ice%alsom                            ,laccu=.FALSE.,lpost=.FALSE.)

    ! old ECHAM5 output

    CALL add (surf,'trfliac'              ,ice%trfliac                 ,code=91 ,laccu=.TRUE. ,lpost=.FALSE., &
              longname='LW flux over ice'                      ,units='W/m**2')
    CALL add (surf,'trflwac'              ,ocean%trflwac               ,code=92 ,laccu=.TRUE. ,lpost=.FALSE., &
              longname='LW flux over water'                    ,units='W/m**2')
    CALL add (surf,'trfllac'              ,land%trfllac                ,code=93 ,laccu=.TRUE. ,lpost=.FALSE., &
              longname='LW flux over land'                     ,units='W/m**2')
    CALL add (surf,'sofliac'              ,ice%sofliac                 ,code=94 ,laccu=.TRUE. ,lpost=.FALSE., &
              longname='SW flux over ice'                      ,units='W/m**2')
    CALL add (surf,'soflwac'              ,ocean%soflwac               ,code=95 ,laccu=.TRUE. ,lpost=.FALSE., &
              longname='SW flux over water'                    ,units='W/m**2')
    CALL add (surf,'sofllac'              ,land%sofllac                ,code=96 ,laccu=.TRUE. ,lpost=.FALSE., &
              longname='SW flux over land'                     ,units='W/m**2')
    CALL add (surf,'friac'                ,box%frac_ice_cover_acc      ,code=97 ,laccu=.TRUE. ,lpost=.FALSE., &
              longname='ice cover (fraction of grid box)'      ,units='')
    CALL add (surf,'tsi'                  ,ice%surface_temperature     ,code=102,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface temperature of ice'            ,units='K')
    CALL add (surf,'tsw'                  ,ocean%surface_temperature   ,code=103,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface temperature of water'          ,units='K')
    CALL add (surf,'ustri'                ,ice%u_stress                ,code=104,laccu=.FALSE.,lpost=.FALSE., &
              longname='zonal wind stress over ice'            ,units='Pa')
    CALL add (surf,'vstri'                ,ice%v_stress                ,code=105,laccu=.FALSE.,lpost=.FALSE., &
              longname='meridional_wind_stress_over_ice'       ,units='Pa')
    CALL add (surf,'ustrw'                ,ocean%u_stress              ,code=106,laccu=.FALSE.,lpost=.FALSE., &
              longname='zonal wind stress over water'          ,units='Pa')
    CALL add (surf,'vstrw'                ,ocean%v_stress              ,code=107,laccu=.FALSE.,lpost=.FALSE., &
              longname='meridional wind stress over water'     ,units='Pa')
    CALL add (surf,'ustrl'                ,land%u_stress               ,code=108,laccu=.FALSE.,lpost=.FALSE., &
              longname='zonal wind stress over land'           ,units='Pa')
    CALL add (surf,'vstrl'                ,land%v_stress               ,code=109,laccu=.FALSE.,lpost=.FALSE., &
              longname='meridional wind stress over land'      ,units='Pa') 
    CALL add (surf,'ahfliac'              ,ice%latent_heat_flux_acc    ,code=110,laccu=.TRUE. ,lpost=.FALSE., &
              longname='latent heat flux over ice'             ,units='W/m**2')
    CALL add (surf,'ahflwac'              ,ocean%latent_heat_flux_acc  ,code=111,laccu=.TRUE. ,lpost=.FALSE., &
              longname='latent heat flux over water'           ,units='W/m**2')
    CALL add (surf,'ahfllac'              ,land%latent_heat_flux_acc   ,code=112,laccu=.TRUE. ,lpost=.FALSE., &
              longname='latent heat flux over land'            ,units='W/m**2')
    CALL add (surf,'evapiac'              ,ice%evaporation_acc         ,code=113,laccu=.TRUE. ,lpost=.FALSE., &
              longname='evaporation over ice'                  ,units='kg/m**2s')
    CALL add (surf,'evapwac'              ,ocean%evaporation_acc       ,code=114,laccu=.TRUE. ,lpost=.FALSE., &
              longname='evaporation over water'                ,units='kg/m**2s')
    CALL add (surf,'evaplac'              ,land%evaporation_acc        ,code=115,laccu=.TRUE. ,lpost=.FALSE., &
              longname='evaporation over land'                 ,units='kg/m**2s')
    CALL add (surf,'az0i'                 ,ice%roughness               ,code=116,laccu=.FALSE.,lpost=.FALSE., &
              longname='roughness length over ice'             ,units='m')
    CALL add (surf,'az0w'                 ,ocean%roughness             ,code=117,laccu=.FALSE.,lpost=.FALSE., &
              longname='roughness length over water'           ,units='m') 
    CALL add (surf,'az0lm'                ,land%roughness_momentum     ,code=118,laccu=.FALSE.,lpost=.FALSE., &
              longname='roughness length of momentum over land',units='m')
    CALL add (surf,'az0lh'                ,land%roughness_heat         ,code=212,laccu=.FALSE.,lpost=.TRUE., &
              longname='roughness length of heat over land'    ,units='m')
    CALL add (surf,'ahfsiac'              ,ice%sensible_flux_acc       ,code=119,laccu=.TRUE. ,lpost=.FALSE., &
              longname='sensible heat flux over ice'           ,units='W/m**2')
    CALL add (surf,'ahfswac'              ,ocean%sensible_flux_acc     ,code=120,laccu=.TRUE. ,lpost=.FALSE., &
              longname='sensible heat flux over water'         ,units='W/m**2')
    CALL add (surf,'ahfslac'              ,land%sensible_flux_acc      ,code=121,laccu=.TRUE. ,lpost=.FALSE., &
              longname='sensible heat flux over land'          ,units='W/m**2')
    CALL add (surf,'alsoi'                ,ice%albedo                  ,code=122,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of ice'                         ,units='')
    CALL add (surf,'alsow'                ,ocean%albedo                ,code=123,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of water'                       ,units='')
    CALL add (surf,'alsol'                ,land%albedo                 ,code=124,laccu=.FALSE.,lpost=.FALSE., &
              longname='albedo of land'                        ,units='')
    CALL add (surf,'ahfice'               ,ice%ahfice                  ,code=125,laccu=.FALSE.,lpost=.FALSE., &
              longname='conductive heat flux'                  ,units='W/m**2')
    CALL add (surf,'qres'                 ,ice%qres                    ,code=126,laccu=.FALSE.,lpost=.FALSE., &
              longname='residual heat flux for melting sea ice',units='W/m**2')
    CALL add (surf,'alake'                ,box%lake_fract              ,code=127,laccu=.FALSE.,lpost=.FALSE., &
              longname='lake fraction of grid box'             ,units='')
    CALL add (surf,'tslm1'                ,land%surface_temperature    ,code=139,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface temperature of land'           ,units='K')
    CALL add (surf,'ahfs'                 ,box%sensible_heat_flux_acc  ,code=146,laccu=.TRUE. ,lpost=.FALSE., &
              longname='sensible heat flux'                    ,units='W/m**2')
    CALL add (surf,'ahfl'                 ,box%latent_heat_flux_acc    ,code=147,laccu=.TRUE. ,lpost=.FALSE., &
              longname='latent heat flux'                      ,units='W/m**2')
    CALL add (surf,'wind10_water'         ,ocean%wind_10_meter         ,code=148,laccu=.FALSE.,lpost=.FALSE., &
              longname='10m windspeed over water'              ,units='m/s')
    CALL add (surf,'u10_mean'             ,box%u_wind_10_meter         ,code=165,laccu=.FALSE.,lpost=.FALSE., &
              longname='10m u-velocity'                        ,units='m/s')
    CALL add (surf,'v10_mean'             ,box%v_wind_10_meter         ,code=166,laccu=.FALSE.,lpost=.FALSE., &
              longname='10m v-velocity'                        ,units='m/s')
    CALL add (surf,'temp2_mean'           ,box%temp_2_meter            ,code=167,laccu=.FALSE.,lpost=.FALSE., &
              longname='2m temperature'                        ,units='K')
    CALL add (surf,'dewpoint_2_meter'     ,box%dewpoint_2_meter        ,code=168,laccu=.FALSE.,lpost=.FALSE., &
              longname='2m dew point temperature'              ,units='K')
    CALL add (surf,'tsurf'                ,box%surface_temperature_acc ,code=169,laccu=.TRUE. ,lpost=.FALSE., &
              longname='surface temperature'                   ,units='K')
    CALL add (surf,'wind10_mean_acc'      ,box%wind_10_meter_acc       ,code=171,laccu=.TRUE. ,lpost=.FALSE., &
              longname='10m windspeed'                         ,units='m/s')
    CALL add (surf,'az0'                  ,box%roughness               ,code=213,laccu=.FALSE.,lpost=.TRUE., &
              longname='roughness length'                      ,units='m')
    CALL add (surf,'albedo_mean'          ,box%albedo                  ,code=175,laccu=.FALSE.,lpost=.FALSE., &
              longname='surface_albedo'                        ,units='')
    CALL add (surf,'ustr'                 ,box%u_stress_acc            ,code=180,laccu=.TRUE. ,lpost=.FALSE., &
              longname='u-stress_mean'                         ,units='Pa')
    CALL add (surf,'vstr'                 ,box%v_stress_acc            ,code=181,laccu=.TRUE. ,lpost=.FALSE., &
              longname='v-stress_mean'                         ,units='Pa')
    CALL add (surf,'evap'                 ,box%evaporation_acc         ,code=182,laccu=.TRUE. ,lpost=.FALSE., &
              longname='evaporation_mean_acc'          ,bits=24,units='kg/m**2s')
    CALL add (surf,'landmask'             ,box%landpingo)
    CALL add (surf,'seamask'              ,box%oceanpingo)
    CALL add (surf,'icemask'              ,box%icepingo)
    CALL add (surf,'t2max_mean'           ,box%maxtemp_2_meter         ,code=201,reset=-99._dp,lpost=.FALSE., &
              longname='maximum 2m temperature'                ,units='K')
    CALL add (surf,'t2min_mean'           ,box%mintemp_2_meter         ,code=202,reset=999._dp,lpost=.FALSE., &
              longname='minimum 2m temperature'                ,units='K')
    CALL add (surf,'ahfcon'               ,ice%ahfcon                  ,code=208,laccu=.TRUE. ,lpost=.FALSE., &
              longname='conductive heat flux through ice'      ,units='W/m**2')
    CALL add (surf,'ahfres'               ,ice%melting                 ,code=209,laccu=.TRUE. ,lpost=.FALSE., &
              longname='melting of ice'                        ,units='W/m**2')
    CALL add (surf,'seaice'               ,box%seaice                  ,code=210,laccu=.FALSE.,lpost=.FALSE., &
              longname='ice cover (fraction of 1-SLM)'         ,units='')
    CALL add (surf,'siced'                ,ice%ice_depth               ,code=211,laccu=.FALSE.,lpost=.FALSE., &
              longname='ice depth'                             ,units='m')
    CALL add (surf,'sni'                  ,ice%snow_water_equivalent   ,code=214,laccu=.FALSE.,lpost=.FALSE., &
              longname='water equivalent of snow on ice'       ,units='m')
    CALL add (surf,'wimax'                ,box%wind_10_max             ,code=216,reset=-99._dp,lpost=.FALSE., &
              longname='maximum 10m-wind speed'                ,units='m/s')

  END SUBROUTINE construct_surface 

  SUBROUTINE destruct_surface

    CALL delete_stream (surf)

  END SUBROUTINE destruct_surface

END MODULE mo_surface_memory
