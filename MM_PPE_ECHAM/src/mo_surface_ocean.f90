#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_surface_ocean
  
  USE mo_echam_convect_tables,  ONLY : lookup_ua_list_spline
  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  ! constants for albedo computation
  REAL(wp), PARAMETER :: calbsea = 0.07_wp    ! sea albedo (diffuse radiation)
  
CONTAINS
 
  !-------------------------------------------------------------------------------------------------
  !
  SUBROUTINE update_albedo_ocean(kdim,mask,                            &
       &             pmu0, pswvis, pswdifvis, pswnir, pswdifnir,       &
       &             palbedo_vis_dir, palbedo_vis_dif,                 &
       &             palbedo_nir_dir, palbedo_nir_dif,                 &
       &             palw1dir, palw2dir, palw1dif, palw2dif,           &
       &             palw1, palw2, palsow    )

  USE mo_control,         ONLY : lrce
  
    IMPLICIT NONE

    ! INPUT
    ! -----

    INTEGER,  INTENT(in)                  :: kdim        ! number of longitudes
    LOGICAL,  INTENT(in)                  :: mask(kdim)
    REAL(wp), INTENT(in), DIMENSION(kdim) :: pswvis      ! net surface visible
    REAL(wp), INTENT(in), DIMENSION(kdim) :: pswdifvis   ! fraction of diffuse visible
    REAL(wp), INTENT(in), DIMENSION(kdim) :: pswnir      ! net surface near infrared
    REAL(wp), INTENT(in), DIMENSION(kdim) :: pswdifnir   ! fraction of diffuse near infrared
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palbedo_vis_dir   ! grid-mean visible albedo (dir)
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palbedo_vis_dif   ! grid-mean visible albedo (dif)
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palbedo_nir_dir   ! grid-mean near infrared albedo (dir)
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palbedo_nir_dif   ! grid-mean near infrared albedo (dif)
    REAL(wp), INTENT(in), DIMENSION(kdim) :: pmu0        ! cos of solar zenith angle

    ! OUTPUT
    ! ------

    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palw1dir   ! water albedo (vis, dir)
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palw2dir   ! water albedo (nir, dir)
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palw1dif   ! water albedo (vis, dif)
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palw2dif   ! water albedo (nir, dif)
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palw1      ! albedo of water (vis, dir+dif)
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palw2      ! albedo of water (nir, dir+dif)
    REAL(wp), INTENT(inout), DIMENSION(kdim) :: palsow     ! water albedo (total)
    
    ! LOCAL
    ! -----

    REAL(wp) :: zalsow(kdim)    ! water albedo (total)

    REAL(wp) :: zeps        ! safety
    REAL(wp) :: zalw        ! function of solar zenith angle (for water albedo)
    REAL(wp) :: zvisp_net   ! net visible direct (parallel) radiation at surface
    REAL(wp) :: zvisd_net   ! net visible diffuse radiation at surface
    REAL(wp) :: znirp_net   ! net near-infrared direct (parallel) radiation at surface
    REAL(wp) :: znird_net   ! net near-infrared diffuse radiation at surface
    REAL(wp) :: zvisp_down  ! downward visible direct (parallel) radiation at surface
    REAL(wp) :: zvisd_down  ! downward visible diffuse radiation at surface
    REAL(wp) :: zvis_down   ! downward visible radiation at surface
    REAL(wp) :: znirp_down  ! downward near-infrared direct (parallel) radiation at surface
    REAL(wp) :: znird_down  ! downward near-infrared diffuse radiation at surface
    REAL(wp) :: znir_down   ! downward near-infrared radiation at surface
    REAL(wp) :: zrad_down   ! downward visible + near-infrared radiation at surface
    
    INTEGER :: jl    ! loop index

    zeps        = EPSILON(1.0_wp)
    zalsow      = palsow

    DO jl = 1, kdim

    IF(mask(jl)) THEN
       !
       !      Visible and near infrared downwelling solar radiation at the surface
       !      (separately for parallel and diffuse radiation)
       !
       zvisd_net     = pswvis(jl)*pswdifvis(jl)
       zvisp_net     = pswvis(jl)*(1.0_wp-pswdifvis(jl))
       znird_net     = pswnir(jl)*pswdifnir(jl)
       znirp_net     = pswnir(jl)*(1.0_wp-pswdifnir(jl))
       zvisp_down    = zvisp_net/(1.0_wp-palbedo_vis_dir(jl))
       zvisd_down    = zvisd_net/(1.0_wp-palbedo_vis_dif(jl))
       znirp_down    = znirp_net/(1.0_wp-palbedo_nir_dir(jl))
       znird_down    = znird_net/(1.0_wp-palbedo_nir_dif(jl))
       zvis_down     = zvisp_down + zvisd_down
       znir_down     = znirp_down + znird_down
       zrad_down     = zvis_down  + znir_down          
       !
       ! Albedo of sea water
       ! -------------------
       !
       IF (zrad_down .GT. 0.0_wp) THEN
          IF (.NOT. lrce) THEN
          !
          ! direct radiation
           zalw       = 0.026_wp/(pmu0(jl)**1.7_wp+0.065_wp)+0.015_wp                    &
                       *(pmu0(jl)-0.1_wp)*(pmu0(jl)-0.5_wp)*(pmu0(jl)-1.0_wp)
          ELSE
           zalw        = 0.07_wp                             !ocean surface albedo in RCE
          END IF
          palw1dir(jl)   = zalw+0.0082_wp
          palw2dir(jl)   = zalw-0.007_wp 
          !
          ! diffuse radiation
          palw1dif(jl)   = calbsea
          palw2dif(jl)   = calbsea
          !
          ! Mean albedo of sea water
          ! ------------------------
          ! dir/p + dif/d
          palw1(jl)      = (palw1dir(jl)*zvisp_down+palw1dif(jl)*zvisd_down)/(zvis_down+zeps)
          palw2(jl)      = (palw2dir(jl)*znirp_down+palw2dif(jl)*znird_down)/(znir_down+zeps)
          ! 1/vis + 2/nir
          palsow(jl)     = (palw1(jl)*zvis_down+palw2(jl)*znir_down)/(zrad_down+zeps)
       ELSE
          palw1dir(jl)   = zalsow(jl)
          palw2dir(jl)   = zalsow(jl)
          palw1dif(jl)   = calbsea
          palw2dif(jl)   = calbsea
          palw1(jl)      = zalsow(jl)
          palw2(jl)      = zalsow(jl)
          palsow(jl)     = zalsow(jl)
       END IF
    END IF   

    END DO

  END SUBROUTINE update_albedo_ocean
  !
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE update_z0_ocean(kdim, klevp1, &
       mask, zcfmw, zudif, zvdif, ptm1, pqm1, zx,&
       paphm1, paz0w&
       )
    
    USE mo_physical_constants,    ONLY: grav, rd, vtmpc1
    USE mo_physc2,       ONLY: cchar
    USE mo_time_control, ONLY: time_step_len

    INTEGER, INTENT(in)     :: kdim, klevp1
    LOGICAL, INTENT(in)     :: mask(kdim)
    REAL(wp), INTENT(in)        :: zcfmw(kdim), zudif(kdim), zvdif(kdim), ptm1(kdim)
    REAL(wp), INTENT(in)        :: pqm1(kdim), zx(kdim), paphm1(kdim, klevp1)
    REAL(wp), INTENT(out)       :: paz0w(kdim)

    REAL(wp) :: zepzzo, ztmst, zcons14

    ztmst  = time_step_len
    zcons14 = cchar * rd / (grav**2 * ztmst)
    zepzzo=1.5e-05_wp

    paz0w = 0.0005_wp

    WHERE(mask)
       
       paz0w(:)  = MAX(zepzzo, zcons14 * zcfmw(:) &
            * SQRT(zudif(:)**2 + zvdif(:)**2) * ptm1(:) &
            * (1._wp + vtmpc1 * pqm1(:) - zx(:)) / paphm1(:,klevp1))
    END WHERE

  END SUBROUTINE UPDATE_Z0_OCEAN
  !
  !-------------------------------------------------------------------------------------------------
  !
  SUBROUTINE update_stress_ocean(kdim, mask, zcfmw, zudif, zvdif &
       , pustrw, pvstrw)

    USE mo_time_control, ONLY: time_step_len
    USE mo_physical_constants,    ONLY: grav
    
    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(wp), INTENT(in)       :: zcfmw(kdim)
    REAL(wp), INTENT(in)       :: zudif(kdim), zvdif(kdim)
    REAL(wp), INTENT(out)      :: pustrw(kdim), pvstrw(kdim)
    
    REAL(wp) :: ztmst, zcons15
    
    ztmst   = time_step_len
    zcons15 = 1._wp / (grav * ztmst)

    pustrw = 0._wp
    pvstrw = 0._wp

    Where (mask)
       pustrw(:) = zcons15 * zcfmw(:) * zudif(:)
       pvstrw(:) = zcons15 * zcfmw(:) * zvdif(:)
    END Where
    
  END SUBROUTINE UPDATE_STRESS_OCEAN
  !
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE precalc_ocean( &
       kdim   , mask &
       , ptsw, paphm1 & 
       , pqm1, zx &
       , ptm1, zqss &
       , zteta1, ztvir1 &
       , zfaxe ,paclc &
       , zlteta1, pum1 &
       , pvm1                                                          &
       , pocu, pocv                                                    &
       , paz0w &
       , pgeom1, zghabl &
       , zqsw, zcptw &
       , zriw, zcfhw &
       , zchw, zbnw &
       , zbmw, zbhw &
       , zustarw, ztkevw &
       , zcfmw &
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
       , zcfncw, ztvw, zcdnw & 
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------
       , zustw &
       )

    USE mo_physical_constants,    ONLY: grav, rd, vtmpc1, cpd, vtmpc2
    USE mo_physc2,       ONLY: ckap, cc, cb, cvdifts, cfreec, cgam
    USE mo_time_control, ONLY: time_step_len

    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(wp),    INTENT(IN)    :: ptsw(kdim), paphm1(kdim), pqm1(kdim), zx(kdim)
    REAL(wp),    INTENT(in)    :: ptm1(kdim), zqss(kdim), zteta1(kdim)
    REAL(wp),    INTENT(IN)    :: ztvir1(kdim), zfaxe(kdim), paclc(kdim), zlteta1(kdim)
    REAL(wp),    INTENT(in)    :: pum1(kdim), pvm1(kdim), paz0w(kdim)
    REAL(wp),    INTENT(in)    :: pocu(kdim), pocv(kdim)
    REAL(wp),    INTENT(IN)    :: pgeom1(kdim), zghabl(kdim)
    
    REAL(wp),    INTENT(out)   :: zqsw(kdim), zcptw(kdim), zriw(kdim)
    REAL(wp),    INTENT(out)   :: zcfhw(kdim), zchw(kdim)
    REAL(wp),    INTENT(out)   :: zbnw(kdim), zbmw(kdim), zbhw(kdim), zustarw(kdim)
    REAL(wp),    INTENT(out)   :: ztkevw(kdim), zcfmw(kdim), zustw(kdim)
    
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
    REAL(wp),    INTENT(out)   :: zcfncw(kdim) 
    REAL(wp),    INTENT(out)   :: zcdnw(kdim) 
    REAL(wp),    INTENT(out)   :: ztvw(kdim) 
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------

    ! Achtung, (klevp1): paphm1 
    ! ACHTUNG, (klev): zteta1, zvir1, pqm1, zfaxe, paclc, zlteta1, pgeom1, pum1, pvm1
    
    ! Local Variables
    INTEGER  :: loidx(kdim), is, nl, jl
    REAL(wp)     :: zes(kdim), zqmitte(kdim), zqtmit(kdim), ztmit(kdim), zqsmit(kdim), zvirmitte(kdim)
    REAL(wp)     :: ztemitte(kdim), zfux(kdim), zfox(kdim), zmult1(kdim), zbuoy(kdim), zdu2oc(kdim)
    REAL(wp)     :: zmult2(kdim), zmult3(kdim), zmult4(kdim), zmult5(kdim), zdus1(kdim), zdus2(kdim)
    REAL(wp)     :: zteldif(kdim), z0h(kdim), zalo(kdim), zaloh(kdim), zcons(kdim), zdthv(kdim) 
    REAL(wp)     :: zucfw(kdim), zscfw(kdim), zcfnchw(kdim), zchnw(kdim), zcr(kdim)
    REAL(wp)     :: zcdn2m(kdim), zcdnr(kdim), zcfm2m(kdim), ztesw(kdim)
    REAL(wp)     :: zwstw(kdim), zconvs(kdim), zmonob(kdim), zstabf(kdim)
    REAL(wp)     :: ua(kdim)
    
    REAL(wp)     :: zkappa, zrvrd, zepdu2, zkap, zcons8, zcons11, zcons12, zepsec, zcons17, zepz0o
    REAL(wp)     :: zrdrv, zepsr, zcons6, zustf
    REAL(wp)     :: zsmn, zshn, zm1, zm2, zm4, zh1, zh2, zwstf
    
    
    ! Constants
    zepsec      = 1.e-2_wp
    zkappa      = rd / cpd
    zrvrd       = vtmpc1 + 1._wp
    zepdu2      = 1.0_wp
    zkap        = ckap
    zcons8      = 2._wp * cb
    zepz0o      = 2._wp
    zcons11     = 3._wp * cb * cc
    zcons12     = cvdifts * time_step_len * grav / rd
    zcons17     = 1._wp / zkap**2
    zrdrv       = 1._wp / zrvrd
    zcons6        = 1._wp / 3._wp
    zh1           = 2.22_wp
    zh2           = 0.22_wp
    zm1           = 1.24_wp
    zm2           = 2.37_wp
    zm4           = 3.69_wp
    zshn          = zh1 * zh2 * SQRT(2._wp)
    zsmn          = zshn * zm1 * zm2 / zm4
    zustf         = 1._wp / zsmn**2
    zwstf         = 0.2_wp
    zepsr         = 1.e-10_wp
 
    zqsw =  0._wp
    zcptw =  0._wp
    zriw =  0._wp
    zcfhw =  0._wp
    zchw =  0._wp
    zbnw =  0._wp
    zbmw =  0._wp
    zbhw =  0._wp
    zustarw =  0._wp
    ztkevw =  0._wp
    zcfmw =  0._wp
    zustw =  0._wp
    loidx = 0
   
!grazia   #78
    ztvw(:)   = 0._wp
    zcdnw(:)  = 0._wp
    zcfncw(:) = 0._wp
!grazia #78

    !-----------------------------------------------------------------------
    !      2.2   surface humidity and virtual temperature
    !
    is = 0
    DO jl=1,kdim
      IF(mask(jl)) THEN
        is = is + 1
        loidx(is) = jl
      ENDIF
    ENDDO
      CALL lookup_ua_list_spline('precalc_ocean', kdim, is, loidx(1), ptsw(1), ua(1))

    DO nl=1,is
      jl = loidx(nl)
      zes(jl)       = ua(nl) / paphm1(jl)
      zqsw(jl)      = zes(jl) / (1._wp- vtmpc1 * zes(jl))
      zcptw(jl)     = ptsw(jl) * cpd * (1._wp+ vtmpc2 * zqsw(jl))
      ztesw(jl)     = ptsw(jl) * (1.e5_wp / paphm1(jl))**zkappa
      ztvw(jl)      = ztesw(jl) * (1._wp + vtmpc1 * zqsw(jl))
    ENDDO
    
    !     ------------------------------------------------------------------
    !        3.     COMPUTATION OF THE EXCHANGE COEFFICIENTS.
    !        3.1       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
    !                  RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
    !                  AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
    !                  COMMON PART OF THE DRAG COEFFICIENTS.
    !
    WHERE(mask)
    !
    ! correction for water and ice points
    !
       zdu2oc(:)    = MAX(zepdu2,(pum1(:)-pocu(:))**2                  &
                                +(pvm1(:)-pocv(:))**2)
       zcons(:)     = zcons12 * paphm1(:) / (ptm1(:) * (1._wp + vtmpc1 * pqm1(:) - zx(:)))
       zqmitte(:)   = (pqm1(:) + zqsw(:)) / 2._wp
       zqtmit(:)    = zx(:) * 0.5_wp + zqmitte(:)
       ztmit(:)     = (ptm1(:) + ptsw(:)) / 2._wp
       zqsmit(:)    = (zqss(:) + zqsw(:)) / 2._wp
       ztemitte(:)  = (zteta1(:) + ztesw(:)) / 2._wp
       zvirmitte(:) = (ztvir1(:) + ztvw(:)) / 2._wp
       zfux(:)      = zfaxe(:) / (cpd * ztmit(:))
       zfox(:)      = zfaxe(:) / (rd * ztmit(:))
       zmult1(:)    = 1._wp + vtmpc1 * zqtmit(:)
       zmult2(:)    = zfux(:) * zmult1 - zrvrd
       zmult3(:)    = zrdrv * zfox(:) * zqsmit(:) / (1._wp + zrdrv * zfox(:) * zfux(:) * zqsmit(:))
       zmult5(:)    = zmult1(:) - zmult2(:) * zmult3(:)
       zmult4(:)    = zfux(:) * zmult5(:) - 1._wp
       zdus1(:)     = paclc(:) * zmult5(:) + ( 1._wp - paclc(:)) * zmult1(:)
       zdus2(:)     = paclc(:) * zmult4(:) + ( 1._wp - paclc(:)) * vtmpc1
       zteldif(:)   = zlteta1(:) - ztesw(:)
       zbuoy(:)     = zdus1(:) * zteldif(:) + zdus2(:) * ztemitte(:) * ((pqm1(:) + zx(:)) - zqsw(:))
       zriw(:)      = pgeom1(:) * zbuoy(:) / (zvirmitte(:) * zdu2oc(:))
       z0h(:)       = paz0w(:) * EXP(2._wp-86.276_wp * paz0w(:)**0.375_wp)
       zalo(:)      = LOG(1._wp + pgeom1(:) / (grav * paz0w(:)))
       zaloh(:)     = LOG(1._wp + pgeom1(:) / (grav * z0h(:)))
       zcdnw(:)     = (zkap / zalo(:))**2
       zchnw(:)     = zkap**2 / (zalo(:) * zaloh(:))
       zucfw(:)     = 1._wp / (1._wp + zcons11 * zcdnw(:) * SQRT(ABS(zriw(:)) &
            * (1._wp + pgeom1(:) / (grav * paz0w(:)))))
       zscfw(:)     = SQRT(1._wp + ABS(zriw(:)))
       zcfncw(:)    = zcons(:) * SQRT(zdu2oc(:)) * zcdnw(:)
       zcfnchw(:)   = zcons(:) * SQRT(zdu2oc(:)) * zchnw(:)
       zdthv(:)     = MAX(0._wp,(ztvw(:) - ztvir1(:)))
       zwstw(:)     = zdthv(:) * SQRT(zdu2oc(:)) / zvirmitte(:)
       zcr(:)       = (cfreec / (zchnw(:) * SQRT(zdu2oc(:)))) * ABS(zbuoy(:))**(1._wp/3._wp)
    ENDWHERE
    
    !----------------------------------------------------------------------------------------------
    !     3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
    !          BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
    !
    WHERE(mask) 
       WHERE(zriw(:).GT.0._wp)
          zcfmw(:)  = zcfncw(:) / (1._wp + zcons8 * zriw(:) / zscfw(:))
          zcfhw(:)  = zcfnchw(:) / (1._wp + zcons8 * zriw(:) * zscfw(:))
          zchw(:)   = zcfhw(:) / zcfnchw(:) * zchnw(:)
       ELSEWHERE
          zcfmw(:)  = zcfncw(:) *(1._wp - zcons8 * zriw(:) * zucfw(:))
          zcfhw(:)  = zcfnchw(:) * (1._wp + zcr(:)**cgam)**(1._wp/cgam)
          zchw(:)   = zcfhw(:) / zcfnchw(:) * zchnw(:)
       ENDWHERE
    ENDWHERE
    
    !-----------------------------------------------------------------------------------------------
    !     interpolation functions for diagnostics
    !
    WHERE(mask)
       zbnw(:)       = zkap / SQRT(zcdnw(:))
       zbmw(:)       = MAX(zepsec, SQRT(zcfmw(:) * zcdnw(:) * zcons17 / zcfncw(:)))
       zbhw(:)       = MAX(zepsec, zchw(:) / zbmw(:) * zcons17)
       zbmw(:)       = 1._wp / zbmw(:)
       zbhw(:)       = 1._wp / zbhw(:)
    ENDWHERE
    
    !-----------------------------------------------------------------------------------------------
    !*       3.4       COMPUTATION OF THE PBL EXTENSION.
    !
    WHERE(mask)
       WHERE(paz0w(:).GT.zepz0o)
          zcdn2m(:)   = (zkap / LOG(1._wp + pgeom1(:) / (grav * zepz0o)))**2
       ELSEWHERE
          zcdn2m(:)   = zcdnw(:)
       ENDWHERE
       
       zcdnr(:)       = zcdn2m(:) / zcdnw(:)
       
       WHERE(paz0w(:).GT.zepz0o.AND.zriw(:).LT.0._wp)
          zcfm2m(:)   = zcfncw(:) * zcdnr(:) * (1._wp - zcons8 *zriw(:) &
               / (1._wp+ zcons11 * zcdn2m(:) * SQRT(ABS(zriw(:)) &
               * (1._wp+ pgeom1(:) / (grav * zepz0o)))))
       ELSEWHERE
          zcfm2m(:)   = zcfmw(:)*zcdnr(:)
       ENDWHERE
       zustw(:)    = zcfm2m(:) * SQRT(zdu2oc(:))
       zustarw(:)  = SQRT(zustw(:) * ptm1(:) &
            * (1._wp + vtmpc1 * pqm1(:) - zx(:)) &
            / (zcons12 * paphm1(:)))
    ENDWHERE
    !----------------------------------------------------------------------------------------------
    !      CONVECTIVE VELOCITY SCALE, MONIN-OBUKHOV LENGTH AND
    !      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)
    WHERE(mask)
       WHERE(zwstw(:) .GT. zepsr)
          zconvs(:) = (zwstw(:) * zchw(:) * zghabl(:))**zcons6
          zmonob(:) = (zustarw(:)**3) / (zkap * grav * zwstw(:) * zchw(:))
          zstabf(:) = (pgeom1(:) / (grav * zmonob(:)))**(zcons6 * 2._wp)
          zstabf(:) = MIN(zustf*3._wp, zstabf(:))
       ELSEWHERE
          zconvs=0._wp
          zstabf=0._wp
       ENDWHERE
       ztkevw(:)    = (zustf + zstabf(:)) * (zustarw(:)**2) + zwstf * (zconvs(:)**2)
    ENDWHERE
    
    
  END SUBROUTINE precalc_ocean
  
  !------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE richtmyer_ocean(kdim, mask, &
       klev, klevp1, klevm1, &
       paphm1, zcfh, &
       zebsh, zqdif, &
       ztdif, zcfhw, &
       zetnw, zftnw, &
       zeqnw, zfqnw &
       )

    USE mo_physc2,       ONLY: cvdifts

    INTEGER,   INTENT(IN)  :: kdim, klev, klevp1, klevm1
    LOGICAL,   INTENT(IN)  :: mask(kdim)
    REAL(wp),      INTENT(IN)  :: zcfhw(kdim)
    REAL(wp),      INTENT(IN)  :: paphm1(kdim,klevp1), zcfh(kdim,klev), ztdif(kdim,klev)
    REAL(wp),      INTENT(IN)  :: zqdif(kdim,klev), zebsh(kdim,klev)
    REAL(wp),      INTENT(OUT) :: zetnw(kdim), zftnw(kdim), zeqnw(kdim), zfqnw(kdim)
    ! klevm1   : zebsh
    ! klevm1   : zcfh
    
    REAL(wp)   :: zdiscw(kdim), zdisqw(kdim)
    REAL(wp)   :: zqdp(kdim), zfac(kdim)
    REAL(wp)   :: ztpfac1
    
    !------------------------
    ! CONSTANTS
    ztpfac1   = cvdifts
    
    zetnw = 0._wp
    zftnw = 0._wp
    zeqnw = 0._wp
    zfqnw = 0._wp

    WHERE(mask)
       zqdp(:)       = 1._wp / (paphm1(:,klevp1) - paphm1(:,klev))
       zfac(:)       = zcfh(:,klevm1) * zqdp(:)

       !*  CALCULATION OF THE EN AND FN COEFFICIENTS OF THE RICHTMYER-
       !*  MORTON-SCHEME CONCERNING THE EQUATION:
       !
       !*  XN = EN * XS + FN
       !
       !*  WITH XN = S_ATM  OR  XN = QATM : ATM. VALUE OF S OR Q
       !*  AND  XS = SSURF  OR  XS = QSAT : SURFACE VALUE OF S OR SAT. SPEC.
       !*                                   HUM. AT THE SURFACE
       !
       zdiscw(:)    = 1._wp / (1._wp + zfac(:) * (1._wp - zebsh(:,klevm1)) + zcfhw(:) * zqdp(:))
       zdisqw(:)    = zdiscw(:)
       zetnw(:)     = zdiscw(:) * zcfhw(:) * zqdp(:)
       zftnw(:)     = zdiscw(:) * (ztdif(:,klev) + zfac(:) * ztdif(:,klevm1)) * ztpfac1
       zeqnw(:)     = zdisqw(:) * zcfhw(:) * zqdp(:)
       zfqnw(:)     = zdisqw(:) * (zqdif(:,klev) + zfac(:) * zqdif(:,klevm1)) * ztpfac1
       
    ENDWHERE
    
  END SUBROUTINE richtmyer_ocean
  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE UPDATE_OCEAN (kdim, mask &
       , zetnw, zftnw &
       , zeqnw, zfqnw &
       , zcptw, zqsw &
       , ztklevw, zqklevw &
       )

    INTEGER, INTENT(IN) :: kdim
    LOGICAL, INTENT(IN) :: mask(kdim)
    REAL(wp), INTENT(IN)    :: zetnw(kdim), zftnw(kdim), zeqnw(kdim), zfqnw(kdim)
    REAL(wp), INTENT(IN)    :: zcptw(kdim), zqsw(kdim)
    REAL(wp), INTENT(OUT)   :: ztklevw(kdim), zqklevw(kdim)
    
    !
    !*  CALCULATION OF SKLEV AND QKLEV USING THE NEW SURFACE VALUES
    !*  ZSNEW AND ZQSNEW WHICH WERE CALCULATED IN SUBROUTINE SURFTEMP
    !
    
    ztklevw = 0._wp
    zqklevw = 0._wp

    WHERE(MASK)
       ztklevw(:)   = zetnw(:) * zcptw(:) + zftnw(:)
       zqklevw(:)   = zeqnw(:) * zqsw(:) + zfqnw(:)
    ENDWHERE
  END SUBROUTINE UPDATE_OCEAN
  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE postproc_ocean( &
       kdim,                 &
       klevp1,    mask,      &
       zcfhw,     zqsw,      & 
       zqdif,     pqm1,      &
       pgeom1,    ztdif,     &
       zcptgz,    zcptw,     &
       ptsw,      zbnw,      &
       zbhw,      zriw,      &
       ptm1,      papm1,     &
       zx,        pum1,      &
       pvm1,                 &
       pocu,      pocv,      & 
       pevapw,               &
       pahfsw,    pahflw,    &
       pwind10w,  zt2w,      &
       paphm1,    zbmw,      &
       zdew2w,    zu10w,     &
       zv10w      &
       )

    USE mo_physical_constants,    ONLY: grav, rd, vtmpc1, alv, tmelt, vtmpc2, cpd, &
                                        c3les, c4les, c3ies, c4ies, c2es
    USE mo_physc2,       ONLY: cvdifts
    USE mo_time_control, ONLY: time_step_len

    INTEGER, INTENT(in)  :: kdim, klevp1
    LOGICAL, INTENT(in)  :: mask(kdim)
    REAL(wp), INTENT(in)     :: zcfhw(kdim), zqsw(kdim)
    REAL(wp), INTENT(in)     :: zqdif(kdim), pqm1(kdim), pgeom1(kdim) 
    REAL(wp), INTENT(in)     :: ztdif(kdim), zcptgz(kdim)
    REAL(wp), INTENT(in)     :: zcptw(kdim), ptsw(kdim), zbnw(kdim), zbhw(kdim), zriw(kdim)
    REAL(wp), INTENT(in)     :: ptm1(kdim), papm1(kdim), zx(kdim)
    REAL(wp), INTENT(in)     :: pum1(kdim), pvm1(kdim), paphm1(kdim,klevp1)
    REAL(wp), INTENT(in)     :: pocu(kdim), pocv(kdim)
    REAL(wp), INTENT(in)     :: zbmw(kdim)
    REAL(wp), INTENT(out)    :: pevapw(kdim), pahfsw(kdim), pahflw(kdim), pwind10w(kdim)
    REAL(wp), INTENT(out)    :: zt2w(kdim), zu10w(kdim), zv10w(kdim), zdew2w(kdim)
    
    INTEGER  :: loidx(kdim), is, nl, jl
    REAL(wp)     :: ztmst, zcons15, ztpfac1, ztpfac2, zcons16, zhuv, zhtq, zephum
    REAL(wp)     :: zcoefw(kdim), zqhflw(kdim), zthflw(kdim)
    REAL(wp)     :: zrat(kdim), zcbn(kdim), zcbs(kdim), zcbu(kdim), zmerge(kdim), zred(kdim)
    REAL(wp)     :: zh2m(kdim), zqs1(kdim), zrh2m(kdim), zcvm3(kdim), zcvm4(kdim)
    REAL(wp)     :: zaph2m(kdim), zqs2(kdim), zq2m(kdim), zfrac(kdim)
    REAL(wp)     :: zmerge1(kdim)
    REAL(wp)     :: ua(kdim)
    
    
    !CONSTANTS
    ztmst         = time_step_len
    ztpfac1       = cvdifts
    ztpfac2       = 1._wp / ztpfac1
    !ztpfac3       = 1._wp - ztpfac2
    zcons15       = 1._wp / (grav * ztmst)
    zcons16       = cpd * vtmpc2
    zhuv          = 10._wp * grav
    zhtq          = 2._wp * grav
    zephum        = 5.e-2_wp
   
    pevapw = 0._wp
    pahfsw = 0._wp
    pahflw = 0._wp
    pwind10w = 0._wp
    zt2w = 0._wp
    zdew2w = 0._wp
    zu10w = 0._wp
    zv10w = 0._wp
    loidx = 0

    WHERE(mask)
       
       !*         5.8     Surface fluxes of heat and moisture
       !
       !*    Moisture fluxes
       zcoefw(:)     = zcons15 * zcfhw(:) * ztpfac2
       zqhflw(:)     = zcoefw(:) * (zqdif(:) - zqsw(:))
       !*    Sensible heat fluxes
       zthflw(:)     = zcoefw(:) * (ztdif(:) - zcptw(:))
       ! 
       pevapw(:)     = zqhflw(:)
       pahfsw(:)     = zthflw(:)
       !     Latent heat fluxes
       pahflw(:)     = alv * zqhflw(:)
    ENDWHERE
       !
       !     Compute new t2m, t2m_max t2m_min
       !
    WHERE(mask)
       zrat(:)       = zhtq / pgeom1(:)
       zcbn(:)       = LOG(1._wp + (EXP (zbnw(:)) - 1._wp) * zrat(:) )
       zcbs(:)       = -(zbnw(:) - zbhw(:)) * zrat(:)
       zcbu(:)       = -LOG(1._wp + (EXP (zbnw(:) - zbhw(:)) - 1._wp) * zrat(:))
       WHERE(zriw(:) .GT. 0._wp)
          zmerge(:)  = zcbs(:)
       ELSEWHERE
          zmerge(:)  = zcbu(:)
       ENDWHERE
       zred(:)       = (zcbn(:) + zmerge(:)) / zbhw(:)
       zh2m(:)       = zcptw(:) + zred(:) * (zcptgz(:) - zcptw(:))
       zt2w(:)       = (zh2m(:) - zhtq ) / (cpd * (1._wp + vtmpc2 * pqm1(:)))    
    ENDWHERE
       !
       !           5.96   2M DEW POINT
       !
    is = 0
    DO jl=1,kdim
      IF(mask(jl)) THEN
        is = is + 1
        loidx(is) = jl
      ENDIF
    ENDDO

    CALL lookup_ua_list_spline('postproc_ocean(1)', kdim, is, loidx(1), ptm1(1), ua(1))

    DO nl=1,is
      jl = loidx(nl)
      zqs1(jl)      = ua(nl) / papm1(jl)
      zqs1(jl)      = zqs1(jl) / (1._wp- vtmpc1 * zqs1(jl))
      zrh2m(jl)     = MAX(zephum, pqm1(jl) / zqs1(jl))
    ENDDO

    WHERE(mask)
       WHERE(zt2w(:) .GT. tmelt)
          zcvm3(:)   = c3les
          zcvm4(:)   = c4les
       ELSEWHERE
          zcvm3(:)   = c3ies
          zcvm4(:)   = c4ies
       ENDWHERE
       zaph2m(:)     = paphm1(:,klevp1) * &
            (1._wp - zhtq / ( rd * zt2w(:) * (1._wp + vtmpc1 * pqm1(:) - zx(:))))
    ENDWHERE

    CALL lookup_ua_list_spline('postproc_ocean(2)', kdim, is, loidx(1), zt2w(1), ua(1))

    DO nl=1,is
      jl = loidx(nl)
      zqs2(jl)      = ua(nl) / zaph2m(jl)
      zqs2(jl)      = zqs2(jl) / (1._wp- vtmpc1 * zqs2(jl))
      zq2m(jl)      = zrh2m(jl) * zqs2(jl)
      zfrac(jl)     = LOG(zaph2m(jl) * zq2m(jl) / (c2es * (1._wp + vtmpc1 * zq2m(jl)))) / zcvm3(jl)
      zdew2w(jl)    = MIN(zt2w(jl), (tmelt - zfrac(jl) * zcvm4(jl)) / (1._wp - zfrac(jl)))
    ENDDO

       !
       !*          5.97   10M WIND COMPONENTS, MAX 10M WINDSPEED
       !
    WHERE(mask)
       zrat(:)       = zhuv / pgeom1(:)
       zcbn(:)       = LOG(1._wp + (EXP (zbnw(:)) - 1._wp) * zrat(:) )

       zcbs(:)       = -(zbnw(:) - zbmw(:)) * zrat(:)
       zcbu(:)       = -LOG(1._wp + (EXP (zbnw(:) - zbmw(:)) - 1._wp) * zrat(:))
       WHERE(zriw(:) .GT. 0._wp)
          zmerge1(:) =zcbs(:)
       ELSEWHERE
          zmerge1(:) =zcbu(:)
       ENDWHERE
       zred(:)       = (zcbn(:) + zmerge1(:)) / zbmw(:)
       zu10w(:)      = zred(:) * pum1(:)
       zv10w(:)      = zred(:) * pvm1(:)
       pwind10w(:)   = zred(:)*SQRT((pum1(:)-pocu(:))**2               &
                                   +(pvm1(:)-pocv(:))**2)
    ENDWHERE
    
  END SUBROUTINE postproc_ocean
  !---------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE s_lake(          &
         kdim     , pseaice   &
       , psiced   , palake    &
       , ptsi     , ptsw      &
       , pahflw   , pahfsw    &
       , pqres    , pfluxres  &
       , ptrflw               &
       , psoflw   , pahfli    &
       , psni     , pcvsi     &
       )
    
    USE mo_physical_constants, ONLY: tmelt, alice, snicecond, dice, rhoilf,     &
                                 hcapmix, hcaprilf, rilfhcap, tfreez, dmix
    USE mo_column,         ONLY: mld, nfor_ts
    USE mo_control,        ONLY: lcolumn, lfractional_mask
    USE mo_time_control,   ONLY: delta_time

    INTEGER,       INTENT(in)   :: kdim
    REAL(wp),      INTENT(in)   :: palake(kdim)
    REAL(wp),      INTENT(inout):: psiced(kdim), pseaice(kdim)
    REAL(wp),      INTENT(inout):: ptsi(kdim),   ptsw(kdim)  
    REAL(wp),      INTENT(in)   :: pahflw(kdim), pahfsw(kdim)
    REAL(wp),      INTENT(in)   :: pqres(kdim)          ! ice melt
    REAL(wp),      INTENT(inout):: pfluxres(kdim)       ! flux residuum
    REAL(wp),      INTENT(in)   :: ptrflw(kdim), psoflw(kdim)  
    REAL(wp),      INTENT(in)   :: pahfli(kdim), psni(kdim)    
    REAL(wp),      INTENT(in)   :: pcvsi(kdim)

    REAL(wp)   :: zfluxw(kdim), zts(kdim), zfres(kdim), zconhflx(kdim)
    REAL(wp)   :: zsubice(kdim), zhi(kdim) 
    REAL(wp)   :: zdtime, zhcapdt, zdthcap, zdtrilf, zrilfdt

    ! CONSTANTS

    zdtime            = delta_time
    zdtrilf           = zdtime/rhoilf    
    zrilfdt           = rhoilf/zdtime   
    zdthcap           = zdtime/hcapmix 
    zhcapdt           = hcapmix/zdtime
    IF(lcolumn .AND. nfor_ts(1) == 0) THEN  ! Use mixed layer depth from namelist
       zdthcap = zdthcap*dmix/mld
       zhcapdt = zhcapdt*mld/dmix
    END IF    


    WHERE ((lfractional_mask .AND. palake(:).GT.0.0_wp) .OR.          &
           (.NOT. lfractional_mask .AND. palake(:).GE.0.5_wp))   
       ! lake points
       WHERE (pseaice(:) .LT. 0.5_wp)     ! open water; for lakes ‘pseaice‘ is either 0. or 1.
          
          zfluxw(:)   = pahflw(:) + pahfsw(:) + ptrflw(:) + psoflw(:)
          
          !--------       Lake temperature (ptsw)         ------------------------------------
          
          zts(:)            = ptsw(:) + zdthcap * (zfluxw(:) + pfluxres(:))
          ptsi(:)           = tmelt
          pfluxres(:)       = 0._wp
          psiced(:)         = 0._wp
          WHERE (zts(:) .GE. tmelt)                                   ! open water (unchanged)
             ptsw(:)        = zts(:)
          ELSEWHERE                                                   ! check ice formation
             ptsw(:)        = tmelt
             zfres(:)       = (zts(:) - tmelt) * zhcapdt              ! < 0.
             WHERE (zts(:) .LE. tmelt - tfreez)                       ! ice formation
                psiced(:)   = hcaprilf * (tmelt - zts(:))             ! >= dice
                pseaice(:)  = 1._wp
             ELSEWHERE
                pfluxres(:)    = zfres(:)
             END WHERE
          END WHERE
       ELSEWHERE (psiced(:) .GE. dice) 
         
          !--------       Ice thickness (psiced)        --------------------------------------
          
          zconhflx(:)       = alice * (ptsi(:) - tmelt) / (psiced(:) + snicecond * psni(:))
          zsubice(:)        = (1._wp - pcvsi(:)) * pahfli(:)
          zhi(:)            = psiced(:)-zdtrilf*(zconhflx(:)+pqres(:)+pfluxres(:)+zsubice(:))
          ptsw(:)           = tmelt
          WHERE (zhi(:) .GE. dice)
             psiced(:)      = zhi(:)
             pseaice(:)     = 1._wp
             pfluxres(:)    = 0._wp
          ELSEWHERE (zhi(:) .LE. 0._wp)                                ! complete melting
             ptsw(:)        = tmelt - zhi(:) * rilfhcap                ! ptsw > tmelt
             psiced(:)      = 0._wp
             pseaice(:)     = 0._wp
             pfluxres(:)    = 0._wp
          ELSEWHERE                                                    ! incomplete melting
             psiced(:)      = dice
             pseaice(:)     = 1._wp
             pfluxres(:)    = (dice - zhi(:)) * zrilfdt                ! > 0
          END WHERE
       END WHERE  ! open water or ice
    END WHERE     ! palake
    
  END SUBROUTINE s_lake

  SUBROUTINE ocean_rad(           &
      kdim         , mask         &
     ,ztrdown      , ptsi         &
     ,pi0          , ptrsol       &
     ,palsoi       , p_albedo     &
     ,psofli       , ptrfli       &
     )

  USE mo_radiation_parameters,         ONLY: cemiss
  USE mo_physical_constants,           ONLY: stbo

  INTEGER, INTENT(in)        :: kdim
  LOGICAL, INTENT(in)        :: mask(kdim)
  REAL(wp),    INTENT(in)    :: ptsi(kdim), pi0(kdim), ptrsol(kdim)
  REAL(wp),    INTENT(in)    :: palsoi(kdim), ztrdown(kdim), p_albedo(kdim)
  REAL(wp),    INTENT(out)   :: psofli(kdim), ptrfli(kdim)

  REAL(wp)                   :: zflxs(kdim)

  psofli = 0._wp
  ptrfli = 0._wp

  WHERE(mask)    
       
     ptrfli(:)  = ztrdown(:) - cemiss * stbo * ptsi(:)**4
     zflxs(:)   = pi0(:) * ptrsol(:)     
     psofli(:)  = (1._wp - palsoi(:)) * zflxs(:) / (1._wp - p_albedo(:))
      
  END WHERE

END SUBROUTINE ocean_rad
  
END MODULE mo_surface_ocean
