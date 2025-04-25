!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_g3b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add => add_stream_element, &
                            default_stream_setting,                   &
                            HYBRID_H, HYBRID

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g3b ! routine to construct the g3b table
  PUBLIC :: destruct_g3b  ! routine to destruct  the g3b table
  PUBLIC :: g3b           ! the g3b table

  ! declaration of predefined fields within this module

  REAL(dp), POINTER, PUBLIC :: geosp(:,:)
  REAL(dp), POINTER, PUBLIC :: tsl(:,:)
  REAL(dp), POINTER, PUBLIC :: ws(:,:)
  REAL(dp), POINTER, PUBLIC :: wl(:,:)
  REAL(dp), POINTER, PUBLIC :: sn(:,:)
  REAL(dp), POINTER, PUBLIC :: slm(:,:)
  REAL(dp), POINTER, PUBLIC :: az0(:,:)
  REAL(dp), POINTER, PUBLIC :: alb(:,:)
  REAL(dp), POINTER, PUBLIC :: forest(:,:)
  REAL(dp), POINTER, PUBLIC :: vgrat(:,:)
  REAL(dp), POINTER, PUBLIC :: vlt(:,:)
  REAL(dp), POINTER, PUBLIC :: wsmx(:,:)
  REAL(dp), POINTER, PUBLIC :: fao(:,:)
  REAL(dp), POINTER, PUBLIC :: aps(:,:)
  REAL(dp), POINTER, PUBLIC :: aprl(:,:)
  REAL(dp), POINTER, PUBLIC :: aprc(:,:)
  REAL(dp), POINTER, PUBLIC :: aprs(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrgw(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrgw(:,:)
  REAL(dp), POINTER, PUBLIC :: vdisgw(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcov(:,:)
  REAL(dp), POINTER, PUBLIC :: temp2(:,:)
  REAL(dp), POINTER, PUBLIC :: dew2(:,:)
  REAL(dp), POINTER, PUBLIC :: wind10(:,:)
  REAL(dp), POINTER, PUBLIC :: sswvis_frc(:,:)
  REAL(dp), POINTER, PUBLIC :: sswdifnir(:,:)
  REAL(dp), POINTER, PUBLIC :: sswdifvis(:,:)
  REAL(dp), POINTER, PUBLIC :: sswpar(:,:)
  REAL(dp), POINTER, PUBLIC :: sswdifpar(:,:)
  REAL(dp), POINTER, PUBLIC :: u10(:,:)
  REAL(dp), POINTER, PUBLIC :: v10(:,:)
  REAL(dp), POINTER, PUBLIC :: srads(:,:)
  REAL(dp), POINTER, PUBLIC :: trads(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0(:,:)
  REAL(dp), POINTER, PUBLIC :: trad0(:,:)
  REAL(dp), POINTER, PUBLIC :: vdis(:,:)
  REAL(dp), POINTER, PUBLIC :: ustr(:,:)
  REAL(dp), POINTER, PUBLIC :: vstr(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfs(:,:)
  REAL(dp), POINTER, PUBLIC :: evap(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfl(:,:)
  REAL(dp), POINTER, PUBLIC :: tslm(:,:)
  REAL(dp), POINTER, PUBLIC :: tslm1(:,:)
  REAL(dp), POINTER, PUBLIC :: emter(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsol(:,:,:)
  REAL(dp), POINTER, PUBLIC :: srad0u(:,:)
  REAL(dp), POINTER, PUBLIC :: sradsu(:,:)
  REAL(dp), POINTER, PUBLIC :: tradsu(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo_vis(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo_nir(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo_vis_dir(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo_vis_dif(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo_nir_dir(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo_nir_dif(:,:)
  REAL(dp), POINTER, PUBLIC :: tsurf(:,:)
  REAL(dp), POINTER, PUBLIC :: seaice(:,:)
  REAL(dp), POINTER, PUBLIC :: siced(:,:)
  REAL(dp), POINTER, PUBLIC :: relhum(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wind10w(:,:)
  REAL(dp), POINTER, PUBLIC :: glac(:,:)
  REAL(dp), POINTER, PUBLIC :: gld(:,:)
  REAL(dp), POINTER, PUBLIC :: aclc(:,:,:)
  REAL(dp), POINTER, PUBLIC :: aclcac(:,:,:)
  REAL(dp), POINTER, PUBLIC :: t2max(:,:)
  REAL(dp), POINTER, PUBLIC :: t2min(:,:)
  REAL(dp), POINTER, PUBLIC :: wimax(:,:)
  REAL(dp), POINTER, PUBLIC :: topmax(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcv(:,:)
  REAL(dp), POINTER, PUBLIC :: qvi(:,:)
  REAL(dp), POINTER, PUBLIC :: xlvi(:,:)
  REAL(dp), POINTER, PUBLIC :: xivi(:,:)
  REAL(dp), POINTER, PUBLIC :: runlnd(:,:)
  REAL(dp), POINTER, PUBLIC :: rgcgn(:,:)
  REAL(dp), POINTER, PUBLIC :: sodif(:,:)
  REAL(dp), POINTER, PUBLIC :: srafs(:,:)
  REAL(dp), POINTER, PUBLIC :: trafs(:,:)
  REAL(dp), POINTER, PUBLIC :: sraf0(:,:)
  REAL(dp), POINTER, PUBLIC :: traf0(:,:)
  REAL(dp), POINTER, PUBLIC :: emtef(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsof(:,:,:)
  REAL(dp), POINTER, PUBLIC :: emtef0(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsof0(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tke(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tkem(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tkem1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: srad0d(:,:)
  REAL(dp), POINTER, PUBLIC :: acdnc(:,:,:)
  REAL(dp), POINTER, PUBLIC :: alake(:,:)
  REAL(dp), POINTER, PUBLIC :: aprflux(:,:)
  REAL(dp), POINTER, PUBLIC :: acvtype(:,:)
  REAL(dp), POINTER, PUBLIC :: xtec(:,:,:)
  REAL(dp), POINTER, PUBLIC :: slf(:,:)
  REAL(dp), POINTER, PUBLIC :: snc(:,:)
  REAL(dp), POINTER, PUBLIC :: rtype(:,:)
  REAL(dp), POINTER, PUBLIC :: rintop(:,:)
  REAL(dp), POINTER, PUBLIC :: apmeb(:,:)
  REAL(dp), POINTER, PUBLIC :: apmebco(:,:)
  REAL(dp), POINTER, PUBLIC :: rain(:,:)
  REAL(dp), POINTER, PUBLIC :: qtnew(:,:)
  REAL(dp), POINTER, PUBLIC :: sicemin(:,:)
  REAL(dp), POINTER, PUBLIC :: sicepdw(:,:)
  REAL(dp), POINTER, PUBLIC :: sicepdi(:,:)
  REAL(dp), POINTER, PUBLIC :: sicepres(:,:)
  REAL(dp), POINTER, PUBLIC :: tsicepdi(:,:)
  REAL(dp), POINTER, PUBLIC :: ameltdepth(:,:)
  REAL(dp), POINTER, PUBLIC :: ameltfrac(:,:)
  REAL(dp), POINTER, PUBLIC :: abso4(:,:)
  REAL(dp), POINTER, PUBLIC :: so4nat(:,:,:)
  REAL(dp), POINTER, PUBLIC :: so4all(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ao3(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tropo(:,:)
  REAL(dp), POINTER, PUBLIC :: taus(:,:)
  !
  REAL(dp), POINTER, PUBLIC :: tpot(:,:,:)
  REAL(dp), POINTER, PUBLIC :: thvvar(:,:,:)
  REAL(dp), POINTER, PUBLIC :: thvsig(:,:)
  !
  !  variables for fractional surface coverage
  !
  REAL(dp), POINTER, PUBLIC :: tsi(:,:)
  REAL(dp), POINTER, PUBLIC :: tsw(:,:)
  REAL(dp), POINTER, PUBLIC :: sni(:,:)
  REAL(dp), POINTER, PUBLIC :: ustri(:,:)
  REAL(dp), POINTER, PUBLIC :: vstri(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrw(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrw(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrl(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrl(:,:)
  REAL(dp), POINTER, PUBLIC :: ahflw(:,:)
  REAL(dp), POINTER, PUBLIC :: az0i(:,:)
  REAL(dp), POINTER, PUBLIC :: az0w(:,:)
  REAL(dp), POINTER, PUBLIC :: az0l(:,:)
  REAL(dp), POINTER, PUBLIC :: alsoi(:,:)
  REAL(dp), POINTER, PUBLIC :: alsow(:,:)
  REAL(dp), POINTER, PUBLIC :: alsol(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfice(:,:)
  REAL(dp), POINTER, PUBLIC :: qres(:,:)
  REAL(dp), POINTER, PUBLIC :: qresmp(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfcon(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfres(:,:)
  REAL(dp), POINTER, PUBLIC :: fluxres(:,:)
  !
  REAL(dp), POINTER, PUBLIC :: ahfsiac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfswac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfslac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfliac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahflwac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfllac(:,:)
  REAL(dp), POINTER, PUBLIC :: evapiac(:,:)
  REAL(dp), POINTER, PUBLIC :: evapwac(:,:)
  REAL(dp), POINTER, PUBLIC :: evaplac(:,:)
  REAL(dp), POINTER, PUBLIC :: trfllac(:,:)
  REAL(dp), POINTER, PUBLIC :: trflwac(:,:)
  REAL(dp), POINTER, PUBLIC :: trfliac(:,:)
  REAL(dp), POINTER, PUBLIC :: sofllac(:,:)
  REAL(dp), POINTER, PUBLIC :: soflwac(:,:)
  REAL(dp), POINTER, PUBLIC :: sofliac(:,:)
  REAL(dp), POINTER, PUBLIC :: friac(:,:)
  REAL(dp), POINTER, PUBLIC :: sh_vdiff(:,:)
  REAL(dp), POINTER, PUBLIC :: ev_vdiff(:,:)
  REAL(dp), POINTER, PUBLIC :: ch_concloud(:,:)
  REAL(dp), POINTER, PUBLIC :: cw_concloud(:,:)
  !
  !  variables for ocean coupling only
  !
  REAL(dp), POINTER, PUBLIC :: awhea(:,:)
  REAL(dp), POINTER, PUBLIC :: awsol(:,:)
  REAL(dp), POINTER, PUBLIC :: awfre(:,:)
  REAL(dp), POINTER, PUBLIC :: resi(:,:)
  REAL(dp), POINTER, PUBLIC :: pmeb(:,:)
  REAL(dp), POINTER, PUBLIC :: awust(:,:)
  REAL(dp), POINTER, PUBLIC :: awvst(:,:)
  REAL(dp), POINTER, PUBLIC :: awsta(:,:)
  REAL(dp), POINTER, PUBLIC :: aicon(:,:)
  REAL(dp), POINTER, PUBLIC :: aiqre(:,:)
  REAL(dp), POINTER, PUBLIC :: aifre(:,:)
  REAL(dp), POINTER, PUBLIC :: aiust(:,:)
  REAL(dp), POINTER, PUBLIC :: aivst(:,:)
  REAL(dp), POINTER, PUBLIC :: ocu(:,:)
  REAL(dp), POINTER, PUBLIC :: ocv(:,:)
  !
  !  variables for sso parametrization
  !
  REAL(dp), POINTER, PUBLIC :: oromea(:,:)
  REAL(dp), POINTER, PUBLIC :: orostd(:,:)
  REAL(dp), POINTER, PUBLIC :: orosig(:,:)
  REAL(dp), POINTER, PUBLIC :: orogam(:,:)
  REAL(dp), POINTER, PUBLIC :: orothe(:,:)
  REAL(dp), POINTER, PUBLIC :: oropic(:,:)
  REAL(dp), POINTER, PUBLIC :: oroval(:,:)
  !
  !  variables for mixed layer ocean only
  !
  REAL(dp), POINTER, PUBLIC :: amlcorr(:,:)
  REAL(dp), POINTER, PUBLIC :: amlcorac(:,:)
  REAL(dp), POINTER, PUBLIC :: amlheatac(:,:)

  !
  !  variables for 200mb radiation
  !
  REAL(dp), POINTER, PUBLIC :: tradl(:,:)
  REAL(dp), POINTER, PUBLIC :: sradl(:,:)
  REAL(dp), POINTER, PUBLIC :: trafl(:,:)
  REAL(dp), POINTER, PUBLIC :: srafl(:,:)
  !
  ! temporary write out of variables from subr. albedos
  !
  REAL(dp), POINTER, PUBLIC ::  fage(:,:)
  REAL(dp), POINTER, PUBLIC ::  snifrac(:,:)
  REAL(dp), POINTER, PUBLIC ::  barefrac(:,:)
  REAL(dp), POINTER, PUBLIC ::  alsom(:,:)
  REAL(dp), POINTER, PUBLIC ::  alsobs(:,:)
  !
  ! aerosol optical parameters for volcanic forcing
  !
  REAL(dp), POINTER, PUBLIC :: ext_sw03(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ssa_sw03(:,:,:)
  REAL(dp), POINTER, PUBLIC :: asy_sw03(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ext_lw08(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ext_lw12(:,:,:)

  ! VARS for mo_surface JSBACH init

!!$  REAL(dp), POINTER, PUBLIC :: csat(:,:)
!!$  REAL(dp), POINTER, PUBLIC :: cair(:,:)
  REAL(dp), POINTER, PUBLIC :: tkevn(:,:)
  REAL(dp), POINTER, PUBLIC :: cdum(:,:)
  REAL(dp), POINTER, PUBLIC :: cduh(:,:)
  REAL(dp), POINTER, PUBLIC :: udif(:,:)
  REAL(dp), POINTER, PUBLIC :: vdif(:,:)
  REAL(dp), POINTER, PUBLIC :: radtemp(:,:)
  REAL(dp), POINTER, PUBLIC :: ustarm(:,:)

! variables for station diagnostics. This list contains variables that are
! accumulated in normal echam output but for the station ouput instantaneous 
! values are needed. They have the names of the accumulated counterparts with 
! a suffix _na (not accumulated).

  REAL(dp), POINTER, PUBLIC :: aclcov_na(:,:)
  REAL(dp), POINTER, PUBLIC :: aprl_na(:,:)
  REAL(dp), POINTER, PUBLIC :: xivi_na(:,:)
  REAL(dp), POINTER, PUBLIC :: xlvi_na(:,:)
  REAL(dp), POINTER, PUBLIC :: evap_na(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfl_na(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfs_na(:,:)
  REAL(dp), POINTER, PUBLIC :: tradsu_na(:,:)
  REAL(dp), POINTER, PUBLIC :: trad0_na(:,:)
  REAL(dp), POINTER, PUBLIC :: traf0_na(:,:)
  REAL(dp), POINTER, PUBLIC :: sradsu_na(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0u_na(:,:)
  REAL(dp), POINTER, PUBLIC :: aprc_na(:,:)
  REAL(dp), POINTER, PUBLIC :: aprs_na(:,:)
  REAL(dp), POINTER, PUBLIC :: qvi_na(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0d_na(:,:)
  REAL(dp), POINTER, PUBLIC :: wind10_na(:,:)
  REAL(dp), POINTER, PUBLIC :: srads_na(:,:)
  REAL(dp), POINTER, PUBLIC :: srafs_na(:,:)
  REAL(dp), POINTER, PUBLIC :: sraf0_na(:,:)
  REAL(dp), POINTER, PUBLIC :: ustr_na(:,:) 
  REAL(dp), POINTER, PUBLIC :: vstr_na(:,:)
  REAL(dp), POINTER, PUBLIC :: trads_na(:,:)
  REAL(dp), POINTER, PUBLIC :: trafs_na(:,:)
  REAL(dp), POINTER, PUBLIC :: tsurf_na(:,:)

  ! declaration of table with 2d- and 3d-field entries

  TYPE (t_stream), POINTER     :: g3b

CONTAINS

  SUBROUTINE construct_g3b

    USE mo_control,    ONLY: lmidatm, lcouple
    !!                      , l_volc

    ! set default attributes for the g3b stream

    CALL default_stream_setting (g3b               &
                                ,lrerun=.TRUE.     &
                                ,lpost=.TRUE.      &
                                ,table=128 ,bits=16)

    ! Add fields to the g3b stream.
    ! despite some 3-d fields (lpost=.FALSE.) these fields are written out by default

    !!
    !!  aerosol optical parameters from sbrt rad_int
    !!
    !IF (l_volc) THEN
    !  CALL add (g3b, 'ext_sw03', ext_sw03, code=31, lpost=.FALSE., leveltype=HYBRID, contnorest=.TRUE.)
    !  CALL add (g3b, 'ssa_sw03', ssa_sw03, code=32, lpost=.FALSE., leveltype=HYBRID, contnorest=.TRUE.)
    !  CALL add (g3b, 'asy_sw03', asy_sw03, code=33, lpost=.FALSE., leveltype=HYBRID, contnorest=.TRUE.)
    !  CALL add (g3b, 'ext_lw08', ext_lw08, code=34, lpost=.FALSE., leveltype=HYBRID, contnorest=.TRUE.)
    !  CALL add (g3b, 'ext_lw12', ext_lw12, code=35, lpost=.FALSE., leveltype=HYBRID, contnorest=.TRUE.)
    !ENDIF
    !
    ! fields for energy budget diagnostic
    !
    CALL add (g3b,'sh_vdiff',     sh_vdiff,     code=64,                                          &
                  longname='column heating due to vertical diffusion', units='W/m**2' )
    CALL add (g3b,'ev_vdiff',     ev_vdiff,     code=65,                                          &
                  longname='column moistening due to vertical diffusion', units='kg/m**2s' )
    CALL add (g3b,'ch_concloud',  ch_concloud,  code=66,                                          &
                  longname='convective heating', units='W/m**s' )
    CALL add (g3b,'cw_concloud',  cw_concloud,  code=67,                                          &
                  longname='convective moistening', units='kg/m**2s' )
    !
    !  temporary write out of fields from subr. albedos
    !
    CALL add (g3b,'fage',      fage      ,code= 68, longname='aging factor of snow on ice'                            )
    CALL add (g3b,'snifrac',   snifrac   ,code= 69, longname='fraction of ice covered with snow'                      )
    CALL add (g3b,'barefrac',  barefrac  ,code= 70, longname='bare ice fraction'                                      )
    CALL add (g3b,'alsom',     alsom     ,code= 71, longname='albedo of melt ponds'                                   )
    CALL add (g3b,'alsobs',    alsobs    ,code= 72, longname='albedo of bare ice and snow without ponds'              )

    CALL add (g3b,'sicepdw',   sicepdw   ,code= 73, longname='melt pond depth on sea-ice '           ,units='m'       )
    CALL add (g3b,'sicepdi',   sicepdi   ,code= 74, longname='ice thickness on melt pond '           ,units='m'       )
    CALL add (g3b,'tsicepdi',  tsicepdi  ,code= 75, longname='ice temperature on frozen melt pond '  ,units='K'       )
    CALL add (g3b,'sicepres',  sicepres  ,code= 76, longname='residual heat flux '                   ,units='W/m**2'  )
    CALL add (g3b,'ameltdepth',ameltdepth,code= 77, longname='total melt pond depth '                ,units='m'       )
    CALL add (g3b,'ameltfrac', ameltfrac ,code= 78, longname='fract area of melt ponds on sea-ice '                   )
    CALL add (g3b,'albedo_vis_dir', albedo_vis_dir ,code= 79    ,longname='surface albedo visible range direct'                    )
    CALL add (g3b,'albedo_nir_dir', albedo_nir_dir ,code= 80    ,longname='surface albedo NIR range direct'                        )
    CALL add (g3b,'albedo_vis_dif', albedo_vis_dif ,code= 81    ,longname='surface albedo visible range diffuse'                   )
    CALL add (g3b,'albedo_nir_dif', albedo_nir_dif ,code= 82    ,longname='surface albedo NIR range diffuse'                       )
    CALL add (g3b,'ocu',       ocu       ,code= 83,lpost=lcouple,longname='ocean eastward velocity'       ,bits=24,units='m/s'     )
    CALL add (g3b,'ocv',       ocv       ,code= 84,lpost=lcouple,longname='ocean northward velocity'      ,bits=24,units='m/s'     )
    !
    !  variables for 200mb radiation
    !
    CALL add (g3b,'tradl',     tradl     ,code= 85,laccu=.TRUE. ,longname='thermal radiation 200mb'               ,units='W/m**2'  )
    CALL add (g3b,'sradl',     sradl     ,code= 86,laccu=.TRUE. ,longname='solar radiation 200mb'                 ,units='W/m**2'  )
    CALL add (g3b,'trafl',     trafl     ,code= 87,laccu=.TRUE. ,longname='thermal radiation 200mb (clear sky)'   ,units='W/m**2'  )
    CALL add (g3b,'srafl',     srafl     ,code= 88,laccu=.TRUE. ,longname='solar radiation 200mb (clear sky)'     ,units='W/m**2'  )
    !
    !  variables for mixed layer ocean only
    !
    CALL add (g3b,'amlcorac',  amlcorac  ,code= 89,laccu=.TRUE. ,longname='mixed layer flux correction'           ,units='W/m**2'  )
    CALL add (g3b,'amlheatac', amlheatac ,code= 90,lpost=.FALSE.,longname='mixed layer heat content' ,laccu=.TRUE.,units='J/m**2'  )
!
    CALL add (g3b,'trfliac',   trfliac   ,code= 91,laccu=.TRUE. ,longname='LW flux over ice'                      ,units='W/m**2'  )
    CALL add (g3b,'trflwac',   trflwac   ,code= 92,laccu=.TRUE. ,longname='LW flux over water'                    ,units='W/m**2'  )
    CALL add (g3b,'trfllac',   trfllac   ,code= 93,laccu=.TRUE. ,longname='LW flux over land'                     ,units='W/m**2'  )
    CALL add (g3b,'sofliac',   sofliac   ,code= 94,laccu=.TRUE. ,longname='SW flux over ice'                      ,units='W/m**2'  )
    CALL add (g3b,'soflwac',   soflwac   ,code= 95,laccu=.TRUE. ,longname='SW flux over water'                    ,units='W/m**2'  )
    CALL add (g3b,'sofllac',   sofllac   ,code= 96,laccu=.TRUE. ,longname='SW flux over land'                     ,units='W/m**2'  )
    CALL add (g3b,'friac',     friac     ,code= 97,laccu=.TRUE. ,longname='ice cover (fraction of grid box)'                       )
    CALL add (g3b,'albedo_vis',albedo_vis,code=100              ,longname='surface albedo visible range'                           )
    CALL add (g3b,'albedo_nir',albedo_nir,code=101              ,longname='surface albedo NIR range'                               )
    CALL add (g3b,'tsi',       tsi       ,code=102              ,longname='surface temperature of ice'            ,units='K'       )
    CALL add (g3b,'tsw',       tsw       ,code=103              ,longname='surface temperature of water'          ,units='K'       )
    CALL add (g3b,'ustri',     ustri     ,code=104              ,longname='zonal      wind stress over ice'       ,units='Pa'      )
    CALL add (g3b,'vstri',     vstri     ,code=105              ,longname='meridional wind stress over ice'       ,units='Pa'      )
    CALL add (g3b,'ustrw',     ustrw     ,code=106              ,longname='zonal      wind stress over water'     ,units='Pa'      )
    CALL add (g3b,'vstrw',     vstrw     ,code=107              ,longname='meridional wind stress over water'     ,units='Pa'      )
    CALL add (g3b,'ustrl',     ustrl     ,code=108              ,longname='zonal      wind stress over land'      ,units='Pa'      )
    CALL add (g3b,'vstrl',     vstrl     ,code=109              ,longname='meridional wind stress over land'      ,units='Pa'      )
    CALL add (g3b,'ahfliac',   ahfliac   ,code=110,laccu=.TRUE. ,longname='latent heat flux over ice'             ,units='W/m**2'  )
    CALL add (g3b,'ahflwac',   ahflwac   ,code=111,laccu=.TRUE. ,longname='latent heat flux over water'           ,units='W/m**2'  )
    CALL add (g3b,'ahfllac',   ahfllac   ,code=112,laccu=.TRUE. ,longname='latent heat flux over land'            ,units='W/m**2'  )
    CALL add (g3b,'evapiac',   evapiac   ,code=113,laccu=.TRUE. ,longname='evaporation over ice'                  ,units='kg/m**2s')
    CALL add (g3b,'evapwac',   evapwac   ,code=114,laccu=.TRUE. ,longname='evaporation over water'                ,units='kg/m**2s')
    CALL add (g3b,'evaplac',   evaplac   ,code=115,laccu=.TRUE. ,longname='evaporation over land'                 ,units='kg/m**2s')
    CALL add (g3b,'az0i',      az0i      ,code=116              ,longname='roughness length over ice'             ,units='m'       )
    CALL add (g3b,'az0w',      az0w      ,code=117              ,longname='roughness length over water'           ,units='m'       )
    CALL add (g3b,'az0l',      az0l      ,code=118              ,longname='roughness length over land'            ,units='m'       )
    CALL add (g3b,'ahfsiac',   ahfsiac   ,code=119,laccu=.TRUE. ,longname='sensible heat flux over ice'           ,units='W/m**2'  )
    CALL add (g3b,'ahfswac',   ahfswac   ,code=120,laccu=.TRUE. ,longname='sensible heat flux over water'         ,units='W/m**2'  )
    CALL add (g3b,'ahfslac',   ahfslac   ,code=121,laccu=.TRUE. ,longname='sensible heat flux over land'          ,units='W/m**2'  )
    CALL add (g3b,'alsoi',     alsoi     ,code=122              ,longname='albedo of ice'                                          )
    CALL add (g3b,'alsow',     alsow     ,code=123              ,longname='albedo of water'                                        )
    CALL add (g3b,'alsol',     alsol     ,code=124              ,longname='albedo of land'                                         )
    CALL add (g3b,'ahfice',    ahfice    ,code=125              ,longname='conductive heat flux'                  ,units='W/m**2'  )
    CALL add (g3b,'qres',      qres      ,code=126              ,longname='residual heat flux for melting sea ice',units='W/m**2'  )
    CALL add (g3b,'alake',     alake     ,code=127,lpost=.FALSE.,longname='lake fraction of grid box'                              )
    CALL add (g3b,'rintop',    rintop    ,code=128,lpost=.FALSE.,longname='low level inversion      '                              )
    CALL add (g3b,'geosp',     geosp     ,code=129              ,longname='surface geopotential (orography)'     ,units='m**2/s**2')
    !               stp                        130                         temperature                                    K
    !                                          131                         u-velocity                                     m/s
    !                                          132                         v-velocity                                     m/s
    !                                          133                         specific humidity                              kg/kg
    CALL add (g3b,'aps',       aps       ,code=134              ,longname='surface pressure'                      ,units='Pa'      )
    !                                          135                         vertical velocity                              Pa/s
    CALL add (g3b,'acdnc',     acdnc     ,code=136,lpost=.FALSE.,longname='cloud droplet number concentration'    ,units='1/m**3'  )
    CALL add (g3b,'apmeb',     apmeb     ,code=137,laccu=.TRUE. ,longname='vertic integr tendenc of water',bits=24,units='kg/m**2s')
    !                          svo             138                         vorticity                                      1/s
    CALL add (g3b,'tslm1',     tslm1     ,code=139              ,longname='surface temperature of land'           ,units='K'       )
    CALL add (g3b,'ws',        ws        ,code=140              ,longname='soil wetness'                          ,units='m'       )
    CALL add (g3b,'sn',        sn        ,code=141              ,longname='snow depth'                            ,units='m'       )
    CALL add (g3b,'aprl',      aprl      ,code=142,laccu=.TRUE. ,longname='large scale precipitation'     ,bits=24,units='kg/m**2s')
    CALL add (g3b,'aprc',      aprc      ,code=143,laccu=.TRUE. ,longname='convective  precipitation'     ,bits=24,units='kg/m**2s')
    CALL add (g3b,'aprs',      aprs      ,code=144,laccu=.TRUE. ,longname='snow fall'                     ,bits=24,units='kg/m**2s')
    CALL add (g3b,'vdis',      vdis      ,code=145,laccu=.TRUE. ,longname='boundary layer dissipation'            ,units='W/m**2'  )
    CALL add (g3b,'ahfs',      ahfs      ,code=146,laccu=.TRUE. ,longname='sensible heat flux'                    ,units='W/m**2'  )
    CALL add (g3b,'ahfl',      ahfl      ,code=147,laccu=.TRUE. ,longname='latent heat flux'                      ,units='W/m**2'  )
    !                                          148                         streamfunction                                 m**2/s
    !                                          149                         velocity potential                             m**2/s
    CALL add (g3b,'xivi',      xivi      ,code=150,laccu=.TRUE. ,longname='vertically integrated cloud ice'       ,units='kg/m**2' )
    !                                          151                         mean sea level pressure                        Pa
    !                          stp(20)         152                         log surface pressure
    !                          xl              153                         cloud water                                    kg/kg
    !                          xi              154                         cloud ice                                      kg/kg
    !                          sd              155                         divergence                                     1/s
    !                                          156                         geopotential height                            gpm
    CALL add (g3b,'relhum',    relhum    ,code=157              ,longname='relative humidity'                                      )
    !                                          158                         tendency of surface pressure                   Pa/s
    CALL add (g3b,'wind10w',   wind10w   ,code=159,lpost=.FALSE.,longname='10m windspeed over water'              ,units='m/s'     )
    CALL add (g3b,'aclc',      aclc      ,code=162,lpost=.FALSE.,longname='cloud cover'                                            )
    CALL add (g3b,'aclcv',     aclcv     ,code=163,lpost=.FALSE.,longname='total cloud cover'                                      )
    CALL add (g3b,'aclcov',    aclcov    ,code=164,laccu=.TRUE. ,longname='total cloud cover'                                      )
    CALL add (g3b,'u10',       u10       ,code=165              ,longname='10m u-velocity'                        ,units='m/s'     )
    CALL add (g3b,'v10',       v10       ,code=166              ,longname='10m v-velocity'                        ,units='m/s'     )
    CALL add (g3b,'temp2',     temp2     ,code=167              ,longname='2m temperature'                        ,units='K'       )
    CALL add (g3b,'dew2',      dew2      ,code=168              ,longname='2m dew point temperature'              ,units='K'       )
    CALL add (g3b,'tsurf',     tsurf     ,code=169,laccu=.TRUE. ,longname='surface temperature'                   ,units='K'       )
    CALL add (g3b,'wind10',    wind10    ,code=171,laccu=.TRUE. ,longname='10m windspeed'                         ,units='m/s'     )
    CALL add (g3b,'slm',       slm       ,code=172              ,longname='land sea mask (1=land, 0=sea/lakes)'                    )
    CALL add (g3b,'az0',       az0       ,code=173,lpost=.FALSE.,longname='roughness length'                      ,units='m'       )
    CALL add (g3b,'alb',       alb       ,code=174,lpost=.FALSE.,longname='surface background albedo'                              )
    CALL add (g3b,'albedo',    albedo    ,code=175              ,longname='surface albedo'                                         )
    CALL add (g3b,'srads',     srads     ,code=176,laccu=.TRUE. ,longname='net surface solar radiation'           ,units='W/m**2'  )
    CALL add (g3b,'trads',     trads     ,code=177,laccu=.TRUE. ,longname='net surface thermal radiation'         ,units='W/m**2'  )
    CALL add (g3b,'srad0',     srad0     ,code=178,laccu=.TRUE. ,longname='net top solar radiation'               ,units='W/m**2'  )
    CALL add (g3b,'trad0',     trad0     ,code=179,laccu=.TRUE. ,longname='top thermal radiation (OLR)'           ,units='W/m**2'  )
    CALL add (g3b,'ustr',      ustr      ,code=180,laccu=.TRUE. ,longname='u-stress'                              ,units='Pa'      )
    CALL add (g3b,'vstr',      vstr      ,code=181,laccu=.TRUE. ,longname='v-stress'                              ,units='Pa'      )
    CALL add (g3b,'evap',      evap      ,code=182,laccu=.TRUE. ,longname='evaporation'                   ,bits=24,units='kg/m**2s')
    CALL add (g3b,'srad0d',    srad0d    ,code=184,laccu=.TRUE. ,longname='top incoming solar radiation'          ,units='W/m**2'  )
    CALL add (g3b,'srafs',     srafs     ,code=185,laccu=.TRUE. ,longname='net surface solar radiation (clear sky)',units='W/m**2' )
    CALL add (g3b,'trafs',     trafs     ,code=186,laccu=.TRUE. ,longname='net surface therm radiation (clear sky)',units='W/m**2' )
    CALL add (g3b,'sraf0',     sraf0     ,code=187,laccu=.TRUE. ,longname='net top solar radiation     (clear sky)',units='W/m**2' )
    CALL add (g3b,'traf0',     traf0     ,code=188,laccu=.TRUE. ,longname='net top thermal radiation   (clear sky)',units='W/m**2' )
    !                                          189                         surface solar cloud forcing                    W/m**2
    !                                          190                         surface thermal cloud forcing                  W/m**2
    !                                          191                         SW top cloud forcing (178-187)                 W/m**2
    !                                          192                         LW top cloud forcing (179-188)                 W/m**2
    CALL add (g3b,'wl',        wl        ,code=193              ,longname='skin reservoir content'                ,units='m'       )
    CALL add (g3b,'slf',       slf       ,code=194              ,longname='sea land fraction'                                      )
    CALL add (g3b,'ustrgw',    ustrgw    ,code=195,lpost=.FALSE.,longname='u-gravity wave stress'                 ,units='Pa'      )
    CALL add (g3b,'vstrgw',    vstrgw    ,code=196,lpost=.FALSE.,longname='v-gravity wave stress'                 ,units='Pa'      )
    CALL add (g3b,'vdisgw',    vdisgw    ,code=197              ,longname='gravity wave dissipation'              ,units='W/m**2'  )
    CALL add (g3b,'vgrat',     vgrat     ,code=198,lpost=.FALSE.,longname='vegetation ratio'                                       )
    CALL add (g3b,'orostd',    orostd    ,code=199,lpost=.FALSE.,longname='orographic standard deviation'         ,units='m'       )
    CALL add (g3b,'vlt',       vlt       ,code=200,lpost=.FALSE.,longname='leaf area index'                                        )
    CALL add (g3b,'t2max',     t2max     ,code=201              ,longname='maximum 2m temperature',reset=-99.0_dp ,units='K'       )
    CALL add (g3b,'t2min',     t2min     ,code=202              ,longname='minimum 2m temperature',reset=999.0_dp ,units='K'       )
    CALL add (g3b,'srad0u',    srad0u    ,code=203,laccu=.TRUE. ,longname='top solar radiation upward'            ,units='W/m**2'  )
    CALL add (g3b,'sradsu',    sradsu    ,code=204,laccu=.TRUE. ,longname='surface solar radiation upward'        ,units='W/m**2'  )
    CALL add (g3b,'tradsu',    tradsu    ,code=205,laccu=.TRUE. ,longname='surface thermal radiation upward'      ,units='W/m**2'  )
    CALL add (g3b,'ahfcon',    ahfcon    ,code=208,laccu=.TRUE. ,longname='conductive heat flux through ice'      ,units='W/m**2'  )
    CALL add (g3b,'ahfres',    ahfres    ,code=209,laccu=.TRUE. ,longname='melting of ice'                        ,units='W/m**2'  )
    CALL add (g3b,'seaice',    seaice    ,code=210              ,longname='ice cover (fraction of 1-SLM)'                          )
    CALL add (g3b,'siced',     siced     ,code=211              ,longname='ice depth'                             ,units='m'       )
    CALL add (g3b,'forest',    forest    ,code=212,lpost=.FALSE.,longname='forest fraction'                                        )
    CALL add (g3b,'gld',       gld       ,code=213              ,longname='glacier depth'                         ,units='m'       )
    CALL add (g3b,'sni',       sni       ,code=214              ,longname='water equivalent of snow on ice'       ,units='m'       )
    CALL add (g3b,'wimax',     wimax     ,code=216              ,longname='maximum 10m-wind speed' ,reset=-99.0_dp,units='m/s'     )
    CALL add (g3b,'topmax',    topmax    ,code=217            ,longname='max height of conv cloud tops',reset=99999.0_dp,units='Pa')
    CALL add (g3b,'runlnd',    runlnd    ,code=220,lpost=.FALSE.,longname='surface runoff not running into ocean' ,units='kg/m**2s')
    CALL add (g3b,'aclcac',    aclcac    ,code=223,laccu=.TRUE. ,longname='cloud cover'                                            )
    CALL add (g3b,'tke',       tke       ,code=224,lpost=.FALSE.,longname='turbulent kinetic energy'             ,units='m**2/s**2')
    CALL add (g3b,'tkem1',     tkem1     ,code=225,lpost=.FALSE.,longname='turbulent kinetic energy (t-1)'       ,units='m**2/s**2')
    CALL add (g3b,'fao',       fao       ,code=226,lpost=.FALSE.,longname='FAO data set (soil data flags 0-5)'                     )
    CALL add (g3b,'rgcgn',     rgcgn     ,code=227,lpost=.FALSE.,longname='volumetric heatcap of soil and landice',units='J/m**3*K')
    CALL add (g3b,'sodif',     sodif     ,code=228,lpost=.FALSE.,longname='diffusivity  of soil and land ice'      ,units='m**2/s' )
    CALL add (g3b,'wsmx',      wsmx      ,code=229              ,longname='field capacity of soil'                 ,units='m'      )
    CALL add (g3b,'qvi',       qvi       ,code=230,laccu=.TRUE. ,longname='vertically integrated water vapor'      ,units='kg/m**2')
    CALL add (g3b,'xlvi',      xlvi      ,code=231,laccu=.TRUE. ,longname='vertically integrated cloud water'      ,units='kg/m**2')
    CALL add (g3b,'glac',      glac      ,code=232              ,longname='fraction of land covered by glaciers'                   )
    CALL add (g3b,'snc',       snc       ,code=233              ,longname='snow depth at the canopy'               ,units='m'      )
    CALL add (g3b,'rtype',     rtype     ,code=234,lpost=.FALSE.,longname='type of convection 0-3 '                                )
    CALL add (g3b,'abso4',     abso4     ,code=235,laccu=.TRUE. ,longname='antropogenic sulfur burden'             ,units='kg/m**2')
    CALL add (g3b,'ao3',       ao3       ,code=236              ,longname='ipcc ozone'           ,leveltype=HYBRID ,units='kg/kg'  )
    CALL add (g3b,'tropo',     tropo     ,code=237              ,longname='WMO defined tropopause height'          ,units='Pa'     )
    CALL add (g3b,'thvsig',    thvsig    ,code=238              ,longname='stddev virt pot temp at halflev klevm1' ,units='K'      )
    CALL add (g3b,'tpot',      tpot      ,code=239              ,longname='potential temperature',leveltype=HYBRID ,units='K'      )
    !                                          259                         windspeed (sqrt(u**2+v**2))
    !                                          260                         total precipitation (142+143)
    !                                          261                         total top radiation (178+179)
    !                                          262                         total surface radiation (176+177)
    !                                          263                         net surface heat flux 
    !                                                                       (146+147+176+177-C*218-208*fice-209); C=3.345E5*fland
    !                                          264                         total surface water (142+143+182-160-221)

    ! Add fields not written to the output stream

    CALL add (g3b,'tsl',     tsl     ,lpost=.FALSE.)
    CALL add (g3b,'tslm',    tslm    ,lpost=.FALSE.)
    CALL add (g3b,'emter',   emter   ,lpost=.FALSE. ,leveltype=HYBRID_H)
    CALL add (g3b,'trsol',   trsol   ,lpost=.FALSE. ,leveltype=HYBRID_H)
    CALL add (g3b,'emtef0',  emtef0  ,lpost=.FALSE. ,leveltype=HYBRID_H)
    CALL add (g3b,'trsof0',  trsof0  ,lpost=.FALSE. ,leveltype=HYBRID_H)
    CALL add (g3b,'emtef',   emtef   ,lpost=.FALSE. ,klev=2)
    CALL add (g3b,'trsof',   trsof   ,lpost=.FALSE. ,klev=2)
    CALL add (g3b,'tkem',    tkem    ,lpost=.FALSE.)
    CALL add (g3b,'acvtype', acvtype ,lpost=.FALSE.)
    CALL add (g3b,'xtec',    xtec    ,lpost=.FALSE.)
    CALL add (g3b,'sswvis_frc',  sswvis_frc  ,lpost=.FALSE.)
    CALL add (g3b,'sswpar',  sswpar  ,lpost=.FALSE.  )
    CALL add (g3b,'sswdifnir',sswdifnir,lpost=.FALSE.)
    CALL add (g3b,'sswdifvis',sswdifvis,lpost=.FALSE.)
    CALL add (g3b,'sswdifpar',sswdifpar,lpost=.FALSE.)
    !
    ! variables for middle atmosphere only
    !
    IF (lmidatm) &
      CALL add (g3b,'aprflux',  aprflux ,lpost=.FALSE.)
    !
    !  variables for fractional surface coverage
    !
    CALL add (g3b,'ahflw',    ahflw    ,lpost=.FALSE.)
    CALL add (g3b,'fluxres',  fluxres  ,lpost=.FALSE.)
    !
    !  variables for mixed layer ocean only
    !
    CALL add (g3b,'amlcorr',  amlcorr  ,lpost=.FALSE.)
    !
    !  variables for ocean coupling only
    !
    CALL add (g3b,'apmebco',  apmebco  ,lpost=.FALSE.)
    CALL add (g3b,'rain',     rain     ,lpost=.FALSE.)
    !
    IF (lcouple) THEN
      CALL add (g3b,'awhea',    awhea   ,lpost=.FALSE.)
      CALL add (g3b,'awsol',    awsol   ,lpost=.FALSE.)
      CALL add (g3b,'resi',     resi    ,code=240, contnorest=.TRUE., lpost=.FALSE.)
      CALL add (g3b,'pmeb',     pmeb    ,code=241, contnorest=.TRUE., lpost=.FALSE.)
      CALL add (g3b,'awust',    awust   ,lpost=.FALSE.)
      CALL add (g3b,'awvst',    awvst   ,lpost=.FALSE.)
      CALL add (g3b,'awsta',    awsta   ,lpost=.FALSE.)
      CALL add (g3b,'aicon',    aicon   ,lpost=.FALSE.)
      CALL add (g3b,'aifre',    aifre   ,lpost=.FALSE.)
      CALL add (g3b,'aiust',    aiust   ,lpost=.FALSE.)
      CALL add (g3b,'aivst',    aivst   ,lpost=.FALSE.)
    END IF
    ! aiqre outside of lcouple because its used for meltpond also
    !
      CALL add (g3b,'aiqre',    aiqre   ,lpost=.FALSE.)

    !
    !  variables for sso parametrization
    !
    CALL add (g3b,'oromea', oromea,lpost=.FALSE.)
    CALL add (g3b,'orosig', orosig,lpost=.FALSE.)
    CALL add (g3b,'orogam', orogam,lpost=.FALSE.)
    CALL add (g3b,'orothe', orothe,lpost=.FALSE.)
    CALL add (g3b,'oropic', oropic,lpost=.FALSE.)
    CALL add (g3b,'oroval', oroval,lpost=.FALSE.)
    !
    ! JSBACH init vars from init field
    !
!!$    CALL add (g3b,'csat',csat ,lpost=.FALSE.)
!!$    CALL add (g3b,'cair',cair ,lpost=.FALSE.)
    CALL add (g3b,'tkevn',tkevn,lpost=.false.)
    CALL add (g3b,'cdum',cdum,lpost=.false.)
    CALL add (g3b,'cduh',cduh,lpost=.false.)
    CALL add (g3b,'udif',udif,lpost=.false.)
    CALL add (g3b,'vdif',vdif,lpost=.false.)
    CALL add (g3b,'radtemp',radtemp,lpost=.false.)
    CALL add (g3b,'ustarm',ustarm,lpost=.false.)
    CALL add (g3b,'qresmp',   qresmp   ,lpost=.FALSE.)
    CALL add (g3b,'qtnew',    qtnew    ,lpost=.FALSE.)
    CALL add (g3b,'so4nat',   so4nat   ,lpost=.FALSE.,leveltype=HYBRID,          &
         longname='natural sulfate'  ,units='kg/kg' )
    CALL add (g3b,'so4all',   so4all   ,lpost=.FALSE.,leveltype=HYBRID,          &
         longname='total sulfur'  ,units='kg/kg' )
    CALL add (g3b,'taus', taus ,lpost=.FALSE., longname='non-dimensional age of snow')
    CALL add (g3b,'sicemin',  sicemin,lpost=.FALSE.,longname='minimum annual ice thickness ' &
            ,units='m'         )
    CALL add (g3b,'thvvar',   thvvar    ,lpost=.FALSE.,leveltype=HYBRID,         &
         longname='variance of virt pot temperature' ,units='K**2' )

! Add stream elements for station diagnostic

  CALL add (g3b,'aclcov_na', aclcov_na, lpost=.FALSE.,                           &
            longname='instantaneous total cloud cover',                          &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'aprl_na', aprl_na, lpost=.FALSE.,                               &
            bits=24,             units='kg/m**2s',                               &
            longname='instantaneous large scale precipitation',                  &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'xivi_na', xivi_na, lpost=.FALSE.,                               &
            units='kg/m**2',                                                     &
            longname='instantaneous vertically integrated cloud ice',            &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'xlvi_na', xlvi_na, lpost=.FALSE.,                               &
            longname='instantaneous verically integrated cloud water',           &
            units='kg/m**2',                                                     &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'evap_na', evap_na, lpost=.FALSE.,                               &
            bits=24,units='kg/m**2s',                                            &
            longname='instantaneous surface evaporation',                        &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'ahfl_na', ahfl_na, lpost=.FALSE.,                               &
            longname='instantaneous latent heat flux',                           &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'ahfs_na', ahfs_na, lpost=.FALSE.,                               &
            longname='instantaneous sensible heat flux',                         &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'tradsu_na', tradsu_na, lpost=.FALSE.,                           &
            longname='instantaneous upward thermal rad. flux at surface',        &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'trad0_na', trad0_na, lpost=.FALSE.,                             &
            longname='instantaneous net thermal rad. flux at top of atm.',       &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'traf0_na', traf0_na, lpost=.FALSE.,                             &
            longname='instantaneous net thermal rad. flux clear sky at top of atm.', &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'sradsu_na', sradsu_na, lpost=.FALSE.,                           &
            longname='instantaneous upward solar rad. flux at surface',          &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'srad0u_na', srad0u_na, lpost=.FALSE.,                           &
            longname='instantaneous upward solar rad. flux at top of atm.',      &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'aprc_na', aprc_na, lpost=.FALSE.,                               &
            longname='instantaneous convective precipitation',                   &
            units='kg/m**2s',                                                    &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'aprs_na', aprs_na, lpost=.FALSE.,                               &
            longname='instantaneous snow fall',                                  &
            units='kg/m**2s',                                                    &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'qvi_na', qvi_na, lpost=.FALSE.,                                 &
            longname='instantaneous vertically integrated water vapour',         &
            units='kg/m**2',                                                     &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'srad0d_na', srad0d_na, lpost=.FALSE.,                           &
            longname='instantaneous incoming solar rad. flux at top of atm.',    &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'wind10_na', wind10_na, lpost=.FALSE.,                           &
            longname='instantaneous 10 metre wind',                              &
            units='m/s',                                                         &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'srads_na', srads_na, lpost=.FALSE.,                             &
            longname='instantaneous net solar rad. flux at surface',             &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'srafs_na', srafs_na, lpost=.FALSE.,                             &
            longname='instantaneous net solar rad. flux clear sky at surface',   &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'sraf0_na', sraf0_na, lpost=.FALSE.,                             &
            longname='instantaneous net solar rad. flux clear sky at top of atm.', &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'ustr_na', ustr_na, lpost=.FALSE.,                               &
            longname='instantaneous zonal wind stress',                          &
            units='Pa',                                                          &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'vstr_na', vstr_na, lpost=.FALSE.,                               &
            longname='instantaneous meridional wind stress',                     &
            units='Pa',                                                          &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'trads_na', trads_na, lpost=.FALSE.,                             &
            longname='instantaneous net thermal rad. flux at surface',           &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'trafs_na', trafs_na, lpost=.FALSE.,                             &
            longname='instantaneous net thermal rad. flux clear sky at surface', &
            units='W/m**2',                                                      &
            lrerun=.FALSE.                                                       )
  CALL add (g3b,'tsurf_na', tsurf_na, lpost=.FALSE.,                             &
            longname='instantaneous surface temperature',                        &
            units='K',                                                           &
            lrerun=.FALSE.                                                       )

  END SUBROUTINE construct_g3b

  SUBROUTINE destruct_g3b

    CALL delete_stream (g3b)

  END SUBROUTINE destruct_g3b

END MODULE mo_memory_g3b
