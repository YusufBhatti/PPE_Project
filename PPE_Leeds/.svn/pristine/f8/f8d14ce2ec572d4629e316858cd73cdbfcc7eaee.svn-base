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
MODULE MO_SURFACE_LAND
  
  USE mo_echam_convect_tables,  ONLY : lookup_ua_list_spline
!   USE mo_surface_memory,       ONLY: stest0, stest1, stest2, stest3, stest4, stest5,&
!       stest6, stest7, stest8, stest9
   USE mo_kind, ONLY: wp

  IMPLICIT NONE
  
CONTAINS
  !
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_stress_land(  &
       kdim          , mask       &
       , zcfml       , zudif      &
       , zvdif       , pustrl     &
       , pvstrl                   &
       )

    USE mo_time_control, ONLY: time_step_len
    USE mo_physical_constants,    ONLY: grav
    
    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(wp), INTENT(in)       :: zcfml(kdim)
    REAL(wp), INTENT(in)       :: zudif(kdim), zvdif(kdim)
    REAL(wp), INTENT(out)      :: pustrl(kdim), pvstrl(kdim)

    REAL(wp) :: ztmst, zcons15

    pustrl = 0._wp
    pvstrl = 0._wp

    ztmst   = time_step_len
    zcons15 = 1._wp / (grav * ztmst)

    WHERE(mask)
       pustrl(:) = zcons15 * zcfml(:) * zudif(:)
       pvstrl(:) = zcons15 * zcfml(:) * zvdif(:)
    ENDWHERE

  END SUBROUTINE UPDATE_STRESS_LAND
  !
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE precalc_land(  &
       kdim       , mask    &
       , ptslm1   , paphm1  &
       , pum1     , pvm1    &
       , pqm1     , zx      &
       , zqss     , zteta1  &
       , ztvir1   , zfaxe   &
       , paclc    , zlteta1 &
       , pgeom1   , paz0lh, paz0lm &
       , ptm1     , zghabl  &              
       , zdqsl    , zril    &
       , zqsl     , zcfncl  &
       , zchl     , zcfhl   &
       , zbnl     , zbhnl   &
       , zbml     , zbhl    &
       , zustarl  , ztkevl  &
       , zcfml              &
       , zustl    , zcptl   &
       , zcpq     , zcair   &
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
       , ztvl     , zcdnl   & 
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------
       , zcsat)
    
    USE mo_physical_constants,    ONLY: grav, rd, vtmpc1, cpd, vtmpc2
    USE mo_physc2,       ONLY: ckap, cc, cb, cvdifts
    USE mo_time_control, ONLY: time_step_len

    INTEGER, INTENT(in)    :: kdim
    LOGICAL, INTENT(in)    :: mask(kdim)
    REAL(wp),    INTENT(IN)    :: ptslm1(kdim), paphm1(kdim), pum1(kdim), pvm1(kdim)
    REAL(wp),    INTENT(in)    :: pqm1(kdim), zx(kdim), zqss(kdim)
    REAL(wp),    INTENT(in)    :: zteta1(kdim), ztvir1(kdim)
    REAL(wp),    INTENT(in)    :: zfaxe(kdim), zlteta1(kdim), zghabl(kdim)
    REAL(wp),    INTENT(in)    :: pgeom1(kdim), paz0lh(kdim), paz0lm(kdim), paclc(kdim), ptm1(kdim)
    REAL(wp),    INTENT(in)    :: zcsat(kdim), zcair(kdim)
    
    REAL(wp),    INTENT(out)   :: zdqsl(kdim), zril(kdim)
    REAL(wp),    INTENT(out)   :: zqsl(kdim) 
    REAL(wp),    INTENT(out)   :: zcfncl(kdim), zchl(kdim), zcfhl(kdim)
    REAL(wp),    INTENT(out)   :: zbnl(kdim), zbhnl(kdim), zbml(kdim), zbhl(kdim), zustarl(kdim)
    REAL(wp),    INTENT(out)   :: ztkevl(kdim), zustl(kdim), zcptl(kdim), zcpq(kdim)
    
    REAL(wp),    INTENT(out)   :: zcfml(kdim)
    
!---------------------------------------------------------------------------
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
    REAL(wp),    INTENT(out)   :: ztvl(kdim) 
    REAL(wp),    INTENT(out)   :: zcdnl(kdim) 
!---------- kai zhang 2009-07-01  for emission and dry deposition ----------
!---------------------------------------------------------------------------

    INTEGER   :: loidx(kdim), is, nl, jl
    REAL(wp)      :: ztesl(kdim)
    REAL(wp)      :: zdu2(kdim), zqmitte(kdim), ztmit(kdim), zqsmit(kdim)
    REAL(wp)      :: ztemitte(kdim)
    REAL(wp)      :: zvirmitte(kdim), zfux(kdim), zfox(kdim), zmult1(kdim),zmult2(kdim), zmult3(kdim)
    REAL(wp)      :: zmult4(kdim), zmult5(kdim), zdus1(kdim), zdus2(kdim), zteldif(kdim), zqddif(kdim)
    REAL(wp)      :: zbuoy(kdim)
    REAL(wp)      :: zchnl(kdim), zucfl(kdim), zucfhl(kdim), zscfl(kdim), zcons(kdim), zcfnchl(kdim)
    REAL(wp)      :: zdthv(kdim), zcdn2m(kdim), zcdnr(kdim), zcfm2m(kdim), zes(kdim), zqtmit(kdim)
    REAL(wp)      :: zwstl(kdim), zconvs(kdim), zmonob(kdim), zstabf(kdim)
    REAL(wp)      :: ua(kdim), dua(kdim)
    
    
    REAL(wp)      :: zcpd, zrd, zkappa, zepdu2, zrvrd, zrdrv, zkap, zcons11, zcons12, ztmst
    REAL(wp)      :: ztpfac1, zcons8, zepsec, zcons17, zepz0o, zcons9, zepsr, zcons6, zustf
    REAL(wp)      :: zsmn, zshn, zm1, zm2, zm4, zh1, zh2, zwstf, zpaphm1, zcor

    ! Constants
    ztmst         = time_step_len
    ztpfac1       = cvdifts
    zepsec        = 1.e-2_wp
    zepsr         = 1.e-10_wp
    
    zcpd          = cpd
    zrd           = rd
    zkappa        = zrd/zcpd
    zepdu2        = 1.0_wp
    zrvrd         = vtmpc1+1._wp
    zrdrv         = 1._wp/ zrvrd
    zkap          = ckap
    zcons11       = 3._wp * cb * cc
    zcons12       = ztpfac1 * ztmst * grav / rd
    zcons8        = 2._wp * cb
    zcons9        = 3._wp * cb
    zcons17       = 1._wp / zkap**2
    zepz0o        = 2._wp
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

    zdqsl = 0._wp
    zril = 0._wp
    zqsl = 0._wp
    zcfncl = 0._wp
    zchl = 0._wp
    zcfhl = 0._wp
    zbnl = 0._wp
    zbhnl = 0._wp
    zbml = 0._wp
    zbhl = 0._wp
    zustarl = 0._wp    
    ztkevl = 0._wp
    zustl = 0._wp
    zcptl = 0._wp
    zcpq = 0._wp
    zcfml = 0._wp
!++mgs    
    ztvl(:) = 0._wp
    zcdnl(:) = 0._wp
!--mgs
    loidx = 0


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

    CALL lookup_ua_list_spline('precalc_land', kdim, is, loidx(1), ptslm1(1), ua(1), dua(1))

    DO nl=1,is
      jl = loidx(nl)
      zpaphm1       = 1._wp / paphm1(jl)
      zes(jl)       = ua(nl) * zpaphm1
      zcor          = 1._wp / (1._wp - vtmpc1 * zes(jl))
      zqsl(jl)      = zes(jl) * zcor
      zdqsl(jl)     = zpaphm1*zcor**2*dua(nl)
      ztesl(jl)     = ptslm1(jl) * (1.e5_wp/ paphm1(jl))**zkappa
      ztvl(jl)      = ztesl(jl) * (1._wp + vtmpc1 * &
                      (zcsat(jl) * zqsl(jl) + (1._wp - zcair(jl)) * pqm1(jl)))
    ENDDO

    !     ------------------------------------------------------------------
    !        3.     COMPUTATION OF THE EXCHANGE COEFFICIENTS.
    !        3.1       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
    !                  RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
    !                  AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
    !                  COMMON PART OF THE DRAG COEFFICIENTS.
    !

    WHERE(mask)
       zdu2(:)     = MAX(zepdu2, pum1(:)**2 + pvm1(:)**2)
       zqmitte(:)  = (pqm1(:) + zcsat(:) * zqsl(:) + (1._wp - zcair(:)) * pqm1(:)) / 2._wp
       zqtmit(:)   = zx(:) * 0.5_wp + zqmitte(:)
       ztmit(:)    = (ptm1(:) + ptslm1(:)) / 2._wp


       zqsmit(:)   = (zqss(:) + zqsl(:)) / 2._wp
       ztemitte(:) = (zteta1(:) + ztesl(:)) / 2._wp
       zvirmitte(:)= (ztvir1(:) + ztvl(:)) / 2._wp
       zfux(:)     = zfaxe(:) / (zcpd * ztmit(:))
       zfox(:)     = zfaxe(:) / (zrd * ztmit(:))
       zmult1(:)   = 1._wp + vtmpc1 * zqtmit(:)
       zmult2(:)   = zfux(:) * zmult1(:) - zrvrd
       zmult3(:)   = zrdrv * zfox(:) * zqsmit(:) / (1._wp + zrdrv * zfox(:) * zfux(:) * zqsmit(:))
       zmult5(:)   = zmult1(:) - zmult2(:) * zmult3(:)
       zmult4(:)   = zfux(:) * zmult5(:) - 1._wp
       zdus1(:)    = paclc(:) * zmult5(:) + (1._wp - paclc(:)) * zmult1(:)
       zdus2(:)    = paclc(:) * zmult4(:) + (1._wp - paclc(:)) * vtmpc1
       zteldif(:)  = zlteta1(:) - ztesl(:)
       zqddif(:)   = zcair(:) * pqm1(:) + zx(:) - zcsat(:) * zqsl(:)


       zbuoy(:)    = zdus1(:) * zteldif(:) + zdus2(:) * ztemitte(:) * zqddif(:)
       zril(:)     = pgeom1(:) * zbuoy(:) / (zvirmitte(:) * zdu2(:))

       zcdnl(:)    = (zkap / LOG(1._wp + pgeom1(:) / (grav * paz0lm(:))))**2

  
       zchnl(:)    = (zkap / LOG(1._wp + pgeom1(:) / (grav * paz0lh(:))))**2


       zucfl(:)    = 1._wp / (1._wp + zcons11 * zcdnl(:) * SQRT(ABS(zril(:)) * ( 1._wp &
            + pgeom1(:) / (grav * paz0lm(:)))))
       zucfhl(:)   = 1._wp / (1._wp + zcons11 * zchnl(:) * SQRT(ABS(zril(:)) * (1._wp  &
            + pgeom1(:) / (grav * paz0lh))))
       zscfl(:)    = SQRT(1._wp + ABS(zril(:)))
       zcons(:)    = zcons12 * paphm1(:) / &
            (ptm1(:) * (1._wp + vtmpc1 * pqm1(:) - zx(:)))
       zcfncl(:)   = zcons(:) * SQRT(zdu2(:)) * zcdnl(:)
       zcfnchl(:)  = zcons(:) * SQRT(zdu2(:)) * zchnl(:)
       zdthv(:)    = MAX(0._wp,(ztvl(:) - ztvir1(:)))
       zwstl(:)    = zdthv(:) * SQRT(zdu2(:)) / zvirmitte(:)

    ENDWHERE
    !----------------------------------------------------------------------------------------------
    !     3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
    !          BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
    !

    WHERE(mask)
       WHERE(zril(:).GT.0._wp)
          zcfml(:) = zcfncl(:) / (1._wp + zcons8 * zril(:) / zscfl(:))
          zcfhl(:) = zcfnchl(:) / (1._wp + zcons8 * zril(:) * zscfl(:))
          zchl(:)  = zcfhl(:) / zcfnchl(:) * zchnl(:)
       ELSEWHERE
          zcfml(:) = zcfncl(:) * (1._wp - zcons8 * zril(:) * zucfl(:))
          zcfhl(:) = zcfnchl(:) * (1._wp - zcons9 * zril(:) * zucfhl(:))
          zchl(:)  = zcfhl(:) / zcfnchl(:) * zchnl(:)
       ENDWHERE
    ENDWHERE

    !-----------------------------------------------------------------------------------------------
    !     interpolation functions for diagnostics
    !

   WHERE(mask)
       zbnl(:)     = zkap / SQRT(zcdnl(:))
       zbhnl(:)    = zkap / SQRT(zchnl(:))
       zbml(:)     = MAX(zepsec, SQRT(zcfml(:) * zcdnl(:) * zcons17 / zcfncl(:)))
       zbhl(:)     = MAX(zepsec, zchl(:) / zbml(:) * zcons17)
       zbml(:)     = 1._wp / zbml(:)
       zbhl(:)     = 1._wp / zbhl(:)
    ENDWHERE

    !-----------------------------------------------------------------------------------------------
    !*       3.4       COMPUTATION OF THE PBL EXTENSION.
    !

  WHERE(mask)
       WHERE(paz0lm(:) .GT. zepz0o)
          zcdn2m(:)= (zkap / LOG(1._wp+ pgeom1(:) / (grav * zepz0o)))**2
       ELSEWHERE
          zcdn2m(:)= zcdnl(:)
       ENDWHERE
       zcdnr(:)    = zcdn2m(:) / zcdnl(:)
       WHERE(paz0lm(:) .GT. zepz0o .AND. zril(:) .LT. 0._wp)
          zcfm2m(:)= zcfncl(:) * zcdnr(:) * (1._wp - zcons8 * zril(:) &
               / (1._wp + zcons11 * zcdn2m(:) * SQRT(ABS(zril(:)) &
               * (1._wp + pgeom1(:) / (grav * zepz0o)))))
       ELSEWHERE
          zcfm2m(:)= zcfml(:) * zcdnr(:)
       ENDWHERE
       zustl(:)    = zcfm2m(:) * SQRT(zdu2(:))
       zustarl(:)  = SQRT(zustl(:) * ptm1(:) &
            * (1._wp + vtmpc1 * pqm1(:) - zx(:)) &
            / (zcons12 * paphm1(:)))
    ENDWHERE

    !----------------------------------------------------------------------------------------------
    !      CONVECTIVE VELOCITY SCALE, MONIN-OBUKHOV LENGTH AND
    !      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)

    WHERE(mask)
       WHERE(zwstl(:) .GT. zepsr)
          zconvs(:)= (zwstl(:) * zchl(:) * zghabl(:))**zcons6
          zmonob(:)= (zustarl(:)**3) / (zkap * grav * zwstl(:) * zchl(:))
          zstabf(:)= (pgeom1(:) / (grav * zmonob(:)))**(zcons6 * 2._wp)
          zstabf(:)= MIN(zustf * 3._wp, zstabf(:))
       ELSEWHERE
          zconvs(:)= 0._wp
          zstabf(:)= 0._wp
       ENDWHERE
       ztkevl(:)   = (zustf + zstabf(:)) * (zustarl(:)**2) + zwstf * (zconvs(:)**2)
    ENDWHERE

    WHERE(mask)
       zcptl(:)    = ptslm1(:) * cpd * (1._wp + vtmpc2 * &
            ( zcsat(:) * zqsl(:) + (1._wp-zcair(:)) * pqm1(:)))
       zcpq(:)     = zcptl(:) / ptslm1(:)
    ENDWHERE

  END SUBROUTINE PRECALC_LAND
  !------------------------------------------------------------------------------------------------------
  !
  SUBROUTINE richtmyer_land( kdim, mask, &
       klev, klevp1, klevm1 , &
       paphm1, zcfh, &
       zebsh, zqdif, &
       ztdif, &
       zcfhl, zcair, &
       zcsat, &
       zetnl, zftnl, &
       zeqnl, zfqnl &
       )

    USE mo_physc2,       ONLY: cvdifts

    INTEGER,   INTENT(IN)  :: kdim, klev, klevp1, klevm1
    LOGICAL,   INTENT(IN)  :: mask(kdim)
    REAL(wp),      INTENT(IN)  :: zcfhl(kdim), zcair(kdim), zcsat(kdim)
    REAL(wp),      INTENT(IN)  :: paphm1(kdim,klevp1), zcfh(kdim,klev), ztdif(kdim,klev)
    REAL(wp),      INTENT(IN)  :: zqdif(kdim,klev), zebsh(kdim,klev)
    REAL(wp),      INTENT(OUT) :: zetnl(kdim), zftnl(kdim), zeqnl(kdim), zfqnl(kdim)
    
    REAL(wp)   :: zdiscl(kdim), zdisql(kdim)
    REAL(wp)   :: zqdp(kdim), zfac(kdim)
    REAL(wp)   :: ztpfac1

    !------------------------
    ! CONSTANTS
    ztpfac1   = cvdifts
    zetnl = 0._wp
    zftnl = 0._wp
    zeqnl = 0._wp
    zfqnl = 0._wp

    WHERE(mask)
       zqdp(:)       = 1._wp / (paphm1(:,klevp1) - paphm1(:,klev))
       zfac(:)       = zcfh(:,klevm1) * zqdp(:)
       !-------------------------------------------------------------------
       !*  CALCULATION OF THE EN AND FN COEFFICIENTS OF THE RICHTMYER-
       !*  MORTON-SCHEME CONCERNING THE EQUATION:
       !
       !*  XN = EN * XS + FN
       !
       !*  WITH XN = S_ATM  OR  XN = QATM : ATM. VALUE OF S OR Q
       !*  AND  XS = SSURF  OR  XS = QSAT : SURFACE VALUE OF S OR SAT. SPEC.
       !*                                   HUM. AT THE SURFACE
       !
       
       zdiscl(:)    = 1._wp / (1._wp + zfac(:) * (1._wp - zebsh(:,klevm1)) + zcfhl(:) * zqdp(:))
       zdisql(:)    = 1._wp / (1._wp + zfac(:) * (1._wp - zebsh(:,klevm1)) + zcair(:) * zcfhl(:) * zqdp(:))
       zetnl(:)     = zdiscl(:) * zcfhl(:) * zqdp(:)
       zftnl(:)     = zdiscl(:) * (ztdif(:,klev) + zfac(:) * ztdif(:,klevm1)) * ztpfac1
       zeqnl(:)     = zdisql(:) * zcsat(:) * zcfhl(:) * zqdp(:)
       zfqnl(:)     = zdisql(:) * (zqdif(:,klev) + zfac(:) * zqdif(:,klevm1)) * ztpfac1

    ENDWHERE    

  END SUBROUTINE richtmyer_land
  !-----------------------------------------------------------------------------------
  !  
  SUBROUTINE UPDATE_LAND ( &
       kdim ,mask          &
       , zetnl, zftnl      &
       , zeqnl, zfqnl      &
       , zcptlnew, zqslnew &
       , ztklevl, zqklevl  &
       )
    
    INTEGER, INTENT(IN) :: kdim  
    LOGICAL, INTENT(IN) :: mask(kdim)
    REAL(wp), INTENT(IN)    :: zetnl(kdim), zftnl(kdim), zeqnl(kdim), zfqnl(kdim)
    REAL(wp), INTENT(IN)    :: zcptlnew(kdim), zqslnew(kdim)
    REAL(wp), INTENT(OUT)   :: ztklevl(kdim), zqklevl(kdim)
    
    !
    !*  CALCULATION OF SKLEV AND QKLEV USING THE NEW SURFACE VALUES
    !*  ZSNEW AND ZQSNEW WHICH WERE CALCULATED IN SUBROUTINE SURFTEMP

    ztklevl = 0._wp
    zqklevl = 0._wp

    WHERE(MASK)
       
       ztklevl(:) = zetnl(:) * zcptlnew(:) + zftnl(:)
       zqklevl(:) = zeqnl(:) * zqslnew(:)  + zfqnl(:)
       
    ENDWHERE

  END SUBROUTINE UPDATE_LAND

  SUBROUTINE UPDATE_ALBEDO_LAND ( &
       kdim ,mask                 &
       , albedo_vis, albedo_nir   &
       , sw_vis, sw_nir           &
       , albedo                  )
    
    INTEGER, INTENT(IN)     :: kdim  
    LOGICAL, INTENT(IN)     :: mask(kdim)
    REAL(wp), INTENT(IN)    :: albedo_vis(kdim), albedo_nir(kdim)
    REAL(wp), INTENT(IN)    :: sw_vis(kdim), sw_nir(kdim)
    REAL(wp), INTENT(OUT)   :: albedo(kdim)

    ! local variables
    REAL(wp)  ::  frac_down_vis(kdim)

    albedo(:) = 0._wp

    !* CALCULATION OF LAND ALBEDO BY WEIGHTING VIS & NIR LAND ALBEDO WITH VIS & NIR INSOLATION
    !* AT GRID POINTS WITH SOLAR INSOLATION (DAY)
    WHERE(MASK(:) .AND. (sw_vis(:) + sw_nir(:)) > 1.0e-09_wp)
       frac_down_vis(:) = (sw_vis(:) / (1.0_wp - albedo_vis(:))) /  &
                          ((sw_vis(:) / (1.0_wp - albedo_vis(:))) + &
                          (sw_nir(:) / (1.0_wp - albedo_nir(:))))

       albedo(:) = frac_down_vis(:) * albedo_vis(:) + (1.0_wp - frac_down_vis(:)) * albedo_nir(:)
    !* CALCULATION OF LAND ALBEDO BY AVERAGING VIS & NIR LAND ALBEDO
    !* At GRID POINTS WITHOUT SOLAR INSOLATION (NIGHT)
    ELSEWHERE(MASK(:))
       albedo(:) = (albedo_vis(:) + albedo_nir(:)) / 2.0_wp
    END WHERE

  END SUBROUTINE UPDATE_ALBEDO_LAND

  !-----------------------------------------------------------------------------------
  !
  SUBROUTINE postproc_land( &
    kdim      , mask        &
    , pgeom1  , zbnl        &
    , zbml    , pum1        &
    , pvm1    , zril        &
    , zspeedl , zbhnl       &
    , zbhl    , zcptl       &
    , zcptgz  , pqm1        &
    , zt2l    , ptm1        &
    , papm1   , paphm1      &
    , klevp1  , zx          &
    , zdew2l  , zu10l       &
    , zv10l &
    )

    USE mo_physical_constants,    ONLY: grav, rd, tmelt, vtmpc1, vtmpc2, cpd, c3les, &
                                        c4les, c3ies, c4ies, c2es

    INTEGER,  INTENT(in)   :: kdim, klevp1
    LOGICAL,  INTENT(in)   :: mask(kdim)
    REAL(wp),     INTENT(in)   :: pgeom1(kdim), zbnl(kdim), zbml(kdim), zril(kdim)
    REAL(wp),     INTENT(in)   :: pum1(kdim), pvm1(kdim)
    REAL(wp),     INTENT(out)  :: zspeedl(kdim)
    REAL(wp),     INTENT(in)   :: zbhnl(kdim), zbhl(kdim), zcptl(kdim)
    REAL(wp),     INTENT(in)   :: zcptgz(kdim), pqm1(kdim)
    REAL(wp),     INTENT(out)  :: zt2l(kdim)
    REAL(wp),     INTENT(in)   :: ptm1(kdim), papm1(kdim), paphm1(kdim,klevp1), zx(kdim)
    REAL(wp),     INTENT(out)  :: zdew2l(kdim), zu10l(kdim), zv10l(kdim)

    REAL(wp)     :: zhtq, zephum, zhuv
    REAL(wp)     :: zrat(kdim), zcbn(kdim), zcbs(kdim), zmerge1(kdim), zred(kdim)
    REAL(wp)     :: zmerge(kdim), zcbu(kdim)
    INTEGER  :: loidx(kdim), is, nl, jl
    REAL(wp)     :: zqs1(kdim), zrh2m(kdim), zcvm3(kdim), zcvm4(kdim), zqs2(kdim)
    REAL(wp)     :: zaph2m(kdim), zq2m(kdim), zfrac(kdim), zh2m(kdim)
    REAL(wp)     :: ua(kdim)

    zspeedl = 0._wp
    zt2l = 0._wp
    zdew2l = 0._wp
    zu10l = 0._wp
    zv10l = 0._wp
    loidx = 0

    zhtq          = 2._wp * grav
    zephum        = 5.e-2_wp
    zhuv          = 10._wp * grav
    !--------------------------------------------------------------------------------
    ! 10M WIND COMPONENTS, MAX 10M WINDSPEED
    WHERE(mask)
       zrat(:)        = zhuv / pgeom1(:)
       zcbn(:)        = LOG(1._wp+ (EXP (zbnl(:)) - 1._wp) * zrat(:))
       zcbs(:)        = -(zbnl(:) - zbml(:)) * zrat(:)
       zcbu(:)        = -LOG(1._wp +(EXP (zbnl(:) - zbml(:)) - 1._wp) * zrat(:))
       WHERE(zril(:).GT.0._wp)
          zmerge1(:)  = zcbs(:)
       ELSEWHERE
          zmerge1(:)  = zcbu(:)
       ENDWHERE
       zred(:)        = (zcbn(:) + zmerge1(:)) / zbml(:)
       zu10l(:)       = zred(:) * pum1(:)
       zv10l(:)       = zred(:) * pvm1(:)
       zspeedl(:)     = SQRT(zu10l(:)**2 + zv10l(:)**2)
    ENDWHERE
    !--------------------------------------------------------------------------------
    ! 2M temperature
    WHERE(mask)    
       zrat(:)        = zhtq / pgeom1(:)
       zcbn(:)        = LOG(1._wp+ (EXP (zbnl(:)) - 1._wp) * zrat(:))

       zcbn(:) = LOG( 1._wp + (EXP ( zbhnl(:)) - 1._wp) * zrat(:))
       zcbs(:) = -(zbhnl(:) - zbhl(:)) * zrat(:)
       zcbu(:) = -LOG(1._wp + (EXP ( zbhnl(:) - zbhl(:)) -1._wp) * zrat(:))
       WHERE(zril(:) .GT. 0._wp)
          zmerge(:)   = zcbs(:)
       ELSEWHERE
          zmerge(:)   = zcbu(:)
       ENDWHERE
       zred(:)    = (zcbn(:) + zmerge(:)) / zbhl(:)
       zh2m(:)    = zcptl(:) + zred(:) * ( zcptgz(:) - zcptl(:))
       zt2l(:)    = (zh2m(:) - zhtq) / (cpd * (1._wp + vtmpc2 * pqm1(:)))
    ENDWHERE
    !--------------------------------------------------------------------------------
    ! 2M dew point
    !
    is = 0
    DO jl=1,kdim
      IF(mask(jl)) THEN
        is = is + 1
        loidx(is) = jl
      ENDIF
    ENDDO

    CALL lookup_ua_list_spline('postproc_land(1)', kdim, is, loidx(1), ptm1(1), ua(1))

    DO nl=1,is
      jl = loidx(nl)
      zqs1(jl)      = ua(nl) / papm1(jl)
      zqs1(jl)      = zqs1(jl) / (1._wp- vtmpc1 * zqs1(jl))
      zrh2m(jl)     = MAX(zephum, pqm1(jl) / zqs1(jl))
    ENDDO

    WHERE(mask)
       WHERE(zt2l(:) .GT. tmelt)
          zcvm3(:)    = c3les
          zcvm4(:)    = c4les
       ELSEWHERE
          zcvm3(:)    = c3ies
          zcvm4(:)    = c4ies
       ENDWHERE
       zaph2m(:)      = paphm1(:,klevp1) * &
            (1._wp - zhtq / (rd * zt2l(:) * (1._wp + vtmpc1 * pqm1(:) - zx(:))))
    ENDWHERE

    CALL lookup_ua_list_spline('postproc_land(2)', kdim, is, loidx(1), zt2l(1), ua(1))

    DO nl=1,is
      jl = loidx(nl)
      zqs2(jl)      = ua(nl) / zaph2m(jl)
      zqs2(jl)      = zqs2(jl) / (1._wp- vtmpc1 * zqs2(jl))
      zq2m(jl)      = zrh2m(jl) * zqs2(jl)
      zfrac(jl)     = LOG(zaph2m(jl) * zq2m(jl) / (c2es * (1._wp + vtmpc1 * zq2m(jl)))) / zcvm3(jl)
      zdew2l(jl)    = MIN(zt2l(jl), (tmelt - zfrac(jl) * zcvm4(jl)) / (1._wp - zfrac(jl)))
    ENDDO

  END SUBROUTINE postproc_land

  SUBROUTINE land_rad(        &
     kdim          , mask         &
     ,ztrdown      , ptslm1       &
     ,pi0          , ptrsol       &
     ,palsoi       , palbedo      &
     ,psofli       , ptrfli       &
     ,ptslnew      , zteffl4      &
     )

  USE mo_radiation_parameters,         ONLY: cemiss
  USE mo_physical_constants,           ONLY: stbo

  INTEGER, INTENT(in)    :: kdim
  LOGICAL, INTENT(in)    :: mask(kdim)
  REAL(wp),    INTENT(in)    :: ptslm1(kdim), pi0(kdim), ptrsol(kdim), ptslnew(kdim)
  REAL(wp),    INTENT(in)    :: palsoi(kdim), ztrdown(kdim), palbedo(kdim)
  REAL(wp),    INTENT(out)   :: psofli(kdim), ptrfli(kdim), zteffl4(kdim)


  REAL(wp) :: zflxs(kdim)

  psofli  = 0._wp
  ptrfli  = 0._wp
  zteffl4 = 0._wp
  WHERE(mask)    
     zteffl4(:) = ptslm1(:)**3*(4._wp*ptslnew(:)-3._wp*ptslm1(:))
   
     ptrfli(:)  = ztrdown(:) - cemiss * stbo * zteffl4(:)
     zflxs(:)   = pi0(:) * ptrsol(:)     
     psofli(:)  = (1._wp - palsoi(:)) * zflxs(:) / (1._wp - palbedo(:))
      
  ENDWHERE

END SUBROUTINE land_rad

END MODULE MO_SURFACE_LAND
