!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE MO_SURFACE_BOUNDARY

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  
  ! used by mo_surface

CONTAINS
  !
  !------------------------------------------------------------------------------------------------- 
  ! Some pre calculations for derived parameters for the lowest atmosphere level
  !
  SUBROUTINE atm_conditions ( kdim, &
       pxlm1   , pxim1,   &
       pgeom1  , ptm1,    &
       pqm1    , papm1,   &
       zx      , zcptgz,  &
       zteta1  , ztvir1,  &
       zfaxe   , zlteta1, &
       zqss               &
       )
    
    USE mo_physical_constants, ONLY: rd, vtmpc1, alv, als, tmelt, cpd, vtmpc2
  USE mo_echam_convect_tables,  ONLY : prepare_ua_index_spline, lookup_ua_spline
    
    INTEGER, INTENT(in)    :: kdim
    REAL(wp), INTENT(in)       :: pxlm1(kdim), pxim1(kdim), pgeom1(kdim)
    REAL(wp), INTENT(in)       :: ptm1(kdim), pqm1(kdim), papm1(kdim)
    REAL(wp), INTENT(out)      :: zx(kdim), zcptgz(kdim), zteta1(kdim), ztvir1(kdim), zfaxe(kdim)
    REAL(wp), INTENT(out)      :: zlteta1(kdim), zqss(kdim)
    
    INTEGER :: idx(kdim)
    REAL(wp)    :: zkappa, zes(kdim), za(kdim), ua(kdim)
    
    !--CONSTANTS
    !
    zkappa       = rd / cpd
    
    !-- MAIN CALCULATIONS
    !
    zx(:)        = pxlm1(:) + pxim1(:) 
    zcptgz(:)    = pgeom1(:) + ptm1(:) * cpd * (1._wp+ vtmpc2 * pqm1(:))
    zteta1(:)    = ptm1(:) * (100000._wp / papm1(:))**zkappa
    ztvir1(:)    = zteta1(:) * (1._wp+ vtmpc1 * pqm1(:) - zx(:))
    WHERE( ptm1(:) .GE. tmelt)
       zfaxe(:)  = alv
    ELSEWHERE
       zfaxe(:)  = als
    ENDWHERE
    zlteta1(:)   = zteta1(:) - zfaxe(:) / cpd * zteta1(:) / ptm1(:) * zx(:)
    CALL prepare_ua_index_spline('atm_conditions', kdim, ptm1(:), idx(:),za(:))
    CALL lookup_ua_spline(kdim, idx(:), za(:), ua(:))
    zes(:)       = MIN(ua(:) / papm1(:) ,0.5_wp)
    zqss(:)      = zes(:) / (1._wp- vtmpc1 * zes(:))
    
  END SUBROUTINE atm_conditions
  !----------------------------------------------------------------------------------------------
  !
  SUBROUTINE blend_zq_zt(     &
       kdim       , pfrl      &
       , pfrw     , pfri      &
       , ztklevl  , ztklevw   &
       , ztklevi  , zqklevl   &
       , zqklevw  , zqklevi   &    
       , ztdif    , zqdif     &
       , land     , ocean     &
       , ice       &
       )

    USE mo_physc2,         ONLY: cvdifts

    INTEGER,   INTENT(in)  :: kdim
    REAL(wp),      INTENT(in)  :: pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(wp),      INTENT(in)  :: ztklevl(kdim), ztklevw(kdim), ztklevi(kdim)
    REAL(wp),      INTENT(in)  :: zqklevl(kdim), zqklevw(kdim), zqklevi(kdim)
    REAL(wp),      INTENT(out) :: ztdif(kdim), zqdif(kdim)
    LOGICAL,   INTENT(in)  :: land(kdim), ocean(kdim), ice(kdim)

    REAL(wp)       :: ztpfac2
    REAL(wp)       :: zero(kdim)

    ztpfac2    = 1._wp / cvdifts
    zero(:)    = 0._wp
    
    ztdif(:)   = (pfrl(:) * MERGE(ztklevl(:),zero,land)  + pfrw(:) * MERGE(ztklevw(:),zero,ocean) &
         + pfri(:) * MERGE(ztklevi(:),zero,ice)) * ztpfac2
    zqdif(:)   = (pfrl(:) * MERGE(zqklevl(:),zero,land) + pfrw(:) * MERGE(zqklevw(:),zero,ocean) &
         + pfri(:) * MERGE(zqklevi(:),zero,ice)) * ztpfac2
    
  END SUBROUTINE blend_zq_zt
  !----------------------------------------------------------------------------------------------
  !
  SUBROUTINE average_ustar( &
       kdim      , klevp1   &
       , pfrl    , pfrw     &
       , pfri    , zustl    &
       , zustw   , zusti    &
       , ptm1    , pqm1     &
       , zx      , paphm1   &
       , zustarm , land     &
       , ocean   , ice      &
       )
        
    USE mo_physical_constants,      ONLY: grav, vtmpc1, rd
    USE mo_physc2,         ONLY: cvdifts
    USE mo_time_control,   ONLY: time_step_len
    
    INTEGER, INTENT(IN)    :: kdim, klevp1
    REAL(wp),    INTENT(IN)    :: pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(wp),    INTENT(IN)    :: zustl(kdim), zustw(kdim), zusti(kdim)
    REAL(wp),    INTENT(in)    :: ptm1(kdim), pqm1(kdim), zx(kdim)
    REAL(wp),    INTENT(in)    :: paphm1(kdim,klevp1)
    LOGICAL, INTENT(in)    :: land(kdim), ocean(kdim), ice(kdim)
    
    REAL(wp),    INTENT(OUT) :: zustarm(kdim)
    
    REAL(wp)  :: zcons12, ztmst, ztpfac1
    REAL(wp)  :: zust(kdim), zero(kdim)
    
    !   CONSTANTS
    ztmst   = time_step_len
    ztpfac1 = cvdifts 
    zcons12 = ztpfac1*ztmst*grav/rd   
    zero(:) = 0._wp
 
    !----------------------------------------------------------     
    !      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)
    !      For pbl-height calculation
    zust(:)        = pfrl(:) * MERGE(zustl(:),zero,land) + pfrw(:) * MERGE(zustw(:),zero,ocean) &
         + pfri(:) * MERGE(zusti(:),zero,ice)    
    zustarm(:)     = SQRT(zust(:) * ptm1(:) *  &
         ( 1._wp + vtmpc1 * pqm1(:) - zx(:)) / (zcons12 * paphm1(:,klevp1)))

  END SUBROUTINE average_ustar
  !
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE average_tvh_qsurf( &
       kdim       , ptslm1      & 
       , pqm1     , ptsw        &
       , ptsi     , pfrl        &
       , pfrw     , pfri        &
       , zcsat    , zcair       &
       , zqsw                   &
       , zqsi     , ztvh        &
       , zqsurf   , zqsl        &
       , land     , ocean       &
       , ice                    &
       )

    USE mo_physical_constants,      ONLY: vtmpc1

    INTEGER, INTENT(in)   :: kdim
    REAL(wp), INTENT(in)      :: ptslm1(kdim), pqm1(kdim), ptsw(kdim), ptsi(kdim)
    REAL(wp), INTENT(in)      :: pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(wp), INTENT(in)      :: zcsat(kdim), zcair(kdim), zqsw(kdim), zqsi(kdim), zqsl(kdim)
    REAL(wp), INTENT(out)     :: ztvh(kdim), zqsurf(kdim)
    LOGICAL, INTENT(in)   :: land(kdim), ocean(kdim), ice(kdim)

    REAL(wp) :: ztvlan(kdim), ztvsea(kdim), ztvice(kdim), zero(kdim)

    zero(:)        = 0._wp

    ztvlan(:)      = MERGE(ptslm1(:),zero,land)  * (1._wp + vtmpc1 *               &
                      MERGE(zcsat(:)*zqsl(:)+(1._wp-zcair(:))*pqm1(:),zero,land))
    ztvsea(:)      = MERGE(ptsw(:)  ,zero,ocean) * (1._wp + vtmpc1 * MERGE(zqsw(:) ,zero,ocean))
    ztvice(:)      = MERGE(ptsi(:)  ,zero,ice)   * (1._wp + vtmpc1 * MERGE(zqsi(:) ,zero,ice))
    ztvh(:)        = pfrl(:) * ztvlan(:) + pfrw(:) * ztvsea(:) + pfri(:) * ztvice(:)
    zqsurf(:)      = pfrl(:) * MERGE(zqsl(:),zero,land) * MERGE(zcsat(:),zero,land) &
         + pfrw(:) * MERGE(zqsw(:),zero,ocean) + pfri(:) * MERGE(zqsi(:),zero,ice)

  END SUBROUTINE average_tvh_qsurf
  !
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE surface_box_average(        &
       kdim      , land     &
       , ice     , ocean    &
       , box_avg            &
       , m_land  , m_ocean  &
       , m_ice   , f_land   &
       , f_ocean , f_ice    &
       )

    INTEGER, INTENT(in)  :: kdim
    REAL(wp),    INTENT(in)  :: land(kdim), ice(kdim), ocean(kdim)
    LOGICAL, INTENT(in)  :: m_land(kdim), m_ocean(kdim), m_ice(kdim)
    REAL(wp),    INTENT(in)  :: f_land(kdim), f_ocean(kdim), f_ice(kdim)
    REAL(wp),    INTENT(out) :: box_avg(kdim)
    REAL(wp)                 :: zero(kdim)
    
    zero(:)         = 0._wp
    box_avg(1:kdim) = f_land * MERGE(land,zero,m_land) + &
         f_ice * MERGE(ice,zero,m_ice)  + f_ocean * MERGE(ocean,zero,m_ocean)
    
  END SUBROUTINE surface_box_average
  !
  !-------------------------------------------------------------------------------------------------
  !
  SUBROUTINE longwave_down_rad(     &
       kdim          , pemter &
       ,ptslm1       , ptsw    &
       ,ptsi         , pfrl    &
       ,pfrw         , pfri    &
       ,ztrdown         &
       )

    USE mo_radiation_parameters,        ONLY: cemiss
    USE mo_physical_constants,          ONLY: stbo

    INTEGER, INTENT(in)    :: kdim
    REAL(wp), INTENT(in)       :: pemter(kdim), ptslm1(kdim), ptsw(kdim)
    REAL(wp), INTENT(in)       :: ptsi(kdim), pfrl(kdim), pfrw(kdim), pfri(kdim)
    REAL(wp), INTENT(out)      :: ztrdown(kdim)

    REAL(wp)    :: zteff4
    INTEGER :: jl
    
     DO jl = 1,kdim
        zteff4=pfrl(jl)*ptslm1(jl)**4                                  &
              +pfri(jl)*ptsi(jl)**4                                    &
              +pfrw(jl)*ptsw(jl)**4
        ztrdown(jl)=pemter(jl)+cemiss*stbo*zteff4
     END DO

   END SUBROUTINE longwave_down_rad
  !
  !-------------------------------------------------------------------------------------------------
  !
END MODULE MO_SURFACE_BOUNDARY
