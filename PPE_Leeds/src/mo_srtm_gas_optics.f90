!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_srtm_gas_optics

  USE mo_kind,       ONLY : wp
  USE mo_rrtm_params, ONLY: nbndsw, ngptsw
  USE mo_srtm_setup,  ONLY: ngb, ngc, ngs

  IMPLICIT NONE

  LOGICAL :: initialized = .FALSE.
  REAL(wp), PARAMETER :: oneminus = 1.0_wp - 1.0e-6_wp
  REAL(wp) ::                &
       rayl_24_js(9,ngptsw), &
       xka( 1175,2,ngptsw) , &
       xk_mnr(   2,ngptsw) , &
       ref_self(10,ngptsw) , &
       ref_frgn( 4,ngptsw) , &
       ref_sflx( 9,ngptsw) , &
       s_wght0(2,nbndsw)   , &
       s_wght1(2,nbndsw)   , &
       s_wght2(2,nbndsw)   , &
       rayl(ngptsw)
  INTEGER  :: igas(4,2,nbndsw), nsp(2,nbndsw), ioff(2,nbndsw), layreffr(nbndsw)
  INTEGER  :: laysolfr_ref(nbndsw), laysolfr

  INTEGER, PARAMETER :: nogas = 0
  INTEGER, PARAMETER :: ih2o  = 1
  INTEGER, PARAMETER :: ico2  = 2
  INTEGER, PARAMETER :: ich4  = 3
  INTEGER, PARAMETER :: io2   = 4
  INTEGER, PARAMETER :: io3   = 5

CONTAINS

  SUBROUTINE gpt_taumol( &
       & klev          ,ig            , &
       & jp            ,jt            ,jt1           ,jk_trop       , &
       & i_self        ,i_frgn        ,gas_abundance ,col_mol       , &
       & x00           ,x01           ,x10           ,x11           , &
       & x_self        ,y_self        ,x_frgn        ,y_frgn        , &
       & sflx_zen      ,taug          ,taur          )

    INTEGER, INTENT (IN) ::  &
         klev               ,&
         ig                 ,&
         jp(:)              ,&
         jt(:)              ,&
         jt1(:)             ,&
         jk_trop            ,&
         i_self(:)          ,&
         i_frgn(:)     

    REAL(wp), INTENT (IN) :: &
         col_mol(:)         ,&
         gas_abundance(:,:) ,&
         x00(:)             ,&
         x01(:)             ,&
         x10(:)             ,&
         x11(:)             ,&
         x_self(:)          ,&
         y_self(:)          ,&
         x_frgn(:)          ,&
         y_frgn(:)         

    REAL(wp), INTENT (INOUT):: &
         sflx_zen         ,&
         taug(:)          ,&
         taur(:)    

    INTEGER, PARAMETER :: noffset(2) = (/1,61/)
    INTEGER  :: jk, js, ii, ib, itrop, jk1, jk2, laysolfr
    INTEGER  :: i000, i001, i010, i011, i100, i101, i110, i111
    REAL(wp) :: fs, gas_mjr, specmult, gas_key, fs_wght
    REAL(wp) :: x000, x001, x010, x011, x100, x101, x110, x111

    IF (.NOT.initialized) CALL setup_taumol(klev)

    js = 0
    ib = ngb(ig)-15
    ! IF (laysolfr_ref(ib) == 0) laysolfr = jk_trop
    laysolfr = jk_trop
    DO jk = 1,klev
      itrop = MAX(0,MIN(1,jk-jk_trop)) + 1    ! itrop = 1 below 100 hPa, 2 above

      gas_key = 0.0_wp
      IF (igas(1,itrop,ib) > 0) gas_key = gas_abundance(jk,igas(1, itrop, ib)) 
      gas_mjr = gas_key

      fs = 0.0_wp
      ii = (jp(jk)*5 - noffset(itrop) )*nsp(itrop,ib) + 1 
      IF (igas(2,itrop,ib) > 0) THEN
        gas_mjr = gas_key + s_wght1(itrop,ib)*gas_abundance(jk,igas(2,itrop,ib))
        specmult = s_wght2(itrop,ib) * MIN(oneminus,gas_key/gas_mjr)
        fs   = specmult - AINT(specmult)
        ii   = ii + INT(specmult)
      END IF

      i000 = ii + (jt(jk) - 5)*nsp(itrop,ib)
      i001 = ii + (jt1(jk))*nsp(itrop,ib)
      i100 = i000+1
      i101 = i001+1
      i110 = i100+ioff(itrop,ib)
      i010 = i000+ioff(itrop,ib)
      i111 = i101+ioff(itrop,ib)
      i011 = i001+ioff(itrop,ib)

      x100 = x00(jk) * fs
      x000 = x00(jk) - x100
      x101 = x01(jk) * fs
      x001 = x01(jk) - x101
      x110 = x10(jk) * fs
      x010 = x10(jk) - x110
      x111 = x11(jk) * fs
      x011 = x11(jk) - x111

      taug(jk) = 0.0_wp

      IF (igas(1,itrop,ib) > 0) taug(jk) = taug(jk) +                   &
           & gas_mjr *                                                        &
           &  ( x000 * xka(i000,itrop,ig) + x001 * xka(i001,itrop,ig)         &
           &  + x010 * xka(i010,itrop,ig) + x011 * xka(i011,itrop,ig) )

      IF (igas(2,itrop,ib) > 0) taug(jk) = taug(jk) +                   &
           & gas_mjr*(                                                        &
           &  + x100 * xka(i100,itrop,ig) + x101 * xka(i101,itrop,ig)         &
           &  + x110 * xka(i110,itrop,ig) + x111 * xka(i111,itrop,ig) )

      IF (igas(3,itrop,ib) > 0) taug(jk) = taug(jk) +                   &
           &   gas_abundance(jk,igas(3, itrop, ib))  *                        &
           &   (x_frgn(jk) * (ref_frgn(i_frgn(jk),ig) + y_frgn(jk) *          &
           &        (ref_frgn(i_frgn(jk)+1,ig) - ref_frgn(i_frgn(jk),ig) )))

      IF (igas(3,itrop,ib) > 0 .AND. x_self(jk) > EPSILON(1.0_wp))            &
           & taug(jk) = taug(jk) + gas_abundance(jk,igas(3, itrop, ib)) &
           &   *(x_self(jk) * (ref_self(i_self(jk),ig) + y_self(jk)   *       &
           &        (ref_self(i_self(jk)+1,ig) - ref_self(i_self(jk),ig))) )

      IF (igas(4,itrop,ib) > 0) taug(jk) = taug(jk) +                   &
           gas_abundance(jk,igas(4, itrop, ib)) * xk_mnr(itrop,ig) 

      taur(jk) = col_mol(jk)*rayl(ig)

      IF (ib == 9 .AND. itrop == 1) THEN
        js = 1 + INT(specmult)
        taur(jk) = col_mol(jk) * (rayl_24_js(js,ig) + &
             fs * (rayl_24_js(js+1,ig) - rayl_24_js(js,ig)))
        js = 0
      END IF

      jk1 = jk+1-itrop
      jk2 = jk1+1
      IF (jp(jk1) < layreffr(ib) .AND. jp(jk2) >= layreffr(ib))               &
           laysolfr = MIN(jk+1,jk_trop)
      IF (laysolfr <=  jk_trop .and. laysolfr_ref(ib) > jk_trop)              &
           laysolfr = jk

      IF (jk == laysolfr .AND. layreffr(ib) >= 0) THEN
        js      = 1  + INT(specmult)
        fs_wght = fs
      END IF
    ENDDO
    !
    ! if an effective layer is specified for the solar flux then replace 
    ! the default solar flux at zenith with one interpolated using a value
    ! at an effective layer -- an RRTM special feature ;-)
    !
    sflx_zen = ref_sflx(1,ig)
    IF (js > 0) sflx_zen =                                             &
         ref_sflx(js,ig) + fs_wght*(ref_sflx(js+1,ig) -ref_sflx(js,ig)) 

  END SUBROUTINE gpt_taumol

  SUBROUTINE setup_taumol(klev)

    USE mo_srtm_setup,  ONLY: nspa, nspb
    USE rrsw_kg16,      ONLY: absa16=>absa, absb16=>absb, forref16=>forref,  &
         & selfref16=>selfref, sflx_ref16 =>sfluxref, rayl16=>rayl, ng16
    USE rrsw_kg17,      ONLY: absa17=>absa, absb17=>absb, forref17=>forref,  &
         & selfref17=>selfref, sflx_ref17 =>sfluxref, rayl17=>rayl, ng17
    USE rrsw_kg18,      ONLY: absa18=>absa, absb18=>absb, forref18=>forref,  &
         & selfref18=>selfref, sflx_ref18 =>sfluxref, rayl18=>rayl, ng18
    USE rrsw_kg19,      ONLY: absa19=>absa, absb19=>absb, forref19=>forref,  &
         & selfref19=>selfref, sflx_ref19 =>sfluxref, rayl19=>rayl, ng19
    USE rrsw_kg20,      ONLY: absa20=>absa, absb20=>absb, forref20=>forref,  &
         & selfref20=>selfref, sflx_ref20 =>sfluxref, rayl20=>rayl, ng20,    &
         & abs_mnr20=>absch4
    USE rrsw_kg21,      ONLY: absa21=>absa, absb21=>absb, forref21=>forref,  &
         & selfref21=>selfref, sflx_ref21 =>sfluxref, rayl21=>rayl, ng21
    USE rrsw_kg22,      ONLY: absa22=>absa, absb22=>absb, forref22=>forref,  &
         & selfref22=>selfref, sflx_ref22 =>sfluxref, rayl22=>rayl, ng22
    USE rrsw_kg23,      ONLY: absa23=>absa, forref23=>forref,                &
         & selfref23=>selfref, sflx_ref23 =>sfluxref, rayl23=>rayl, ng23
    USE rrsw_kg24,      ONLY: absa24=>absa, absb24=>absb, forref24=>forref,  &
         & selfref24=>selfref, sflx_ref24=>sfluxref, ng24, abs_m24a=>abso3a, &
         & abs_m24b=>abso3b, rayl24a=>rayla, rayl24b=>raylb
    USE rrsw_kg25,      ONLY: absa25=>absa, sflx_ref25 =>sfluxref, ng25,     &
         & rayl25=>rayl, abs_m25a=>abso3a,abs_m25b=>abso3b
    USE rrsw_kg26,      ONLY: sflx_ref26=>sfluxref, rayl26=>rayl, ng26
    USE rrsw_kg27,      ONLY: absa27=>absa, absb27=>absb,                    &
         sflx_ref27=>sfluxref, rayl27=>rayl, ng27
    USE rrsw_kg28,      ONLY: absa28=>absa, absb28=>absb,                    &
         sflx_ref28=>sfluxref, rayl28=>rayl, ng28
    USE rrsw_kg29,      ONLY: absa29=>absa, absb29=>absb, forref29=>forref,  &
         selfref29=>selfref, sflx_ref29=>sfluxref, ng29, rayl29=>rayl,       &
         abs_m29a=>absco2,abs_m29b=>absh2o

    INTEGER, INTENT (IN) :: klev

    INTEGER :: ig, jg, ibeg, iband, ib
    !
    !  Prepare band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
    !
    ibeg  = 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/9,1/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o,ich4,ih2o,nogas,     &
         &                               ich4,nogas,nogas,nogas/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/252.131_wp, 0.0_wp/)
    s_wght2(:,ib) = (/8.0_wp, 8.0_wp/)
    layreffr(ib)  = -18
    laysolfr_ref(ib)  = klev

    DO ig=ibeg,ibeg+ng16-1
      jg = ig - ibeg + 1
      rayl(ig)      = rayl16
      xka(1:585,1,ig)   = absa16(:,jg)
      xka(1:235,2,ig)   = absb16(:,jg)
      ref_self(1:10,ig) = selfref16(:,jg)
      ref_frgn(1:3 ,ig) = forref16(:,jg)
      ref_sflx(1,ig)    = sflx_ref16(jg)
    END DO
    !
    !  Prepare band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/9,5/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o,ico2,ih2o,nogas,  &
         &                              ih2o,ico2,ih2o,nogas/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.364641_wp, 0.364641_wp/)
    s_wght2(:,ib) = (/8.0_wp, 4.0_wp/)
    layreffr(ib)  = 30
    laysolfr_ref(ib)  = klev

    DO ig=ibeg,ibeg+ng17-1
      jg = ig - ibeg + 1
      rayl(ig)          = rayl17
      xka(1:585,1,ig)   = absa17(:,jg)
      xka(1:1175,2,ig)  = absb17(:,jg)
      ref_self(1:10,ig) = selfref17(:,jg)
      ref_frgn(1:4 ,ig) = forref17(:,jg)
      ref_sflx(1:5,ig)  = sflx_ref17(jg,:)
    END DO
    !
    !  Prepare band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)    = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)   = (/9,1/)
    igas(:,:,ib) = RESHAPE ( SOURCE = (/ih2o,ich4,ih2o,nogas,     &
         &                              ich4,nogas,nogas,nogas/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/38.9589_wp, 0.0_wp/)
    s_wght2(:,ib) = (/8.0_wp, 0.0_wp/)
    layreffr(ib) = 6
    laysolfr_ref(ib)  = 0

    DO ig=ibeg,ibeg+ng18-1
      jg = ig - ibeg + 1
      rayl(ig)          = rayl18
      xka(1:585,1,ig)   = absa18(:,jg)
      xka(1:235,2,ig)   = absb18(:,jg)
      ref_self(1:10,ig) = selfref18(:,jg)
      ref_frgn(1:3 ,ig) = forref18(:,jg)
      ref_sflx(1:9,ig)  = sflx_ref18(jg,:)
    END DO
    !
    !  Prepare band 19: 4650-5150 cm-1 (low - h2o,co2; high - co2)
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/9,1/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o,ico2,ih2o,nogas,    &
         &                              ico2,nogas,nogas,nogas/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/5.49281_wp, 0.0_wp/)
    s_wght2(:,ib) = (/8.0_wp, 0.0_wp/)
    layreffr(ib)  = 3
    laysolfr_ref(ib)  = 0

    DO ig=ibeg,ibeg+ng19-1
      jg = ig - ibeg + 1
      rayl(ig)          = rayl19
      xka(1:585,1,ig)   = absa19(:,jg)
      xka(1:235,2,ig)   = absb19(:,jg)
      ref_self(1:10,ig) = selfref19(:,jg)
      ref_frgn(1:3 ,ig) = forref19(:,jg)
      ref_sflx(1:9,ig)  = sflx_ref19(jg,:)
    END DO
    !
    !  Prepare band 20: 5150-6150 cm-1 (low - h2o; high - h2o)
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/1,1/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o,nogas,ih2o,ich4,  &
         &                              ih2o,nogas,ih2o,ich4/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -3
    laysolfr_ref(ib)  = 0

    DO ig=ibeg,ibeg+ng20-1
      jg = ig - ibeg + 1
      rayl(ig)          = rayl20
      xk_mnr(1,ig)      = abs_mnr20(jg)
      xk_mnr(2,ig)      = abs_mnr20(jg)
      xka(1:65,1,ig)    = absa20(:,jg)
      xka(1:235,2,ig)   = absb20(:,jg)
      ref_self(1:10,ig) = selfref20(:,jg)
      ref_frgn(1:4 ,ig) = forref20(:,jg)
      ref_sflx(1,ig)    = sflx_ref20(jg)
    END DO
    !
    !  Prepare band 21: 6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/9,5/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o,ico2,ih2o,nogas, &
         &                              ih2o,ico2,ih2o,nogas/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0045321_wp, 0.0045321_wp/)
    s_wght2(:,ib) = (/8.0_wp, 4.0_wp/)
    layreffr(ib)  = 8
    laysolfr_ref(ib)  = 0

    DO ig=ibeg,ibeg+ng21-1
      jg = ig - ibeg + 1
      rayl(ig)          = rayl21
      xka(1:585,1,ig)   = absa21(:,jg)
      xka(1:1175,2,ig)  = absb21(:,jg)
      ref_self(1:10,ig) = selfref21(:,jg)
      ref_frgn(1:4 ,ig) = forref21(:,jg)
      ref_sflx(1:9,ig)  = sflx_ref21(jg,:)
    END DO
    !
    ! Prepare band 22: 6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
    ! 
    ! In this band the ratio of the ratio of total O2 band intensity (lines 
    ! and Mate continuum) to O2 band intensity (line only) is 1.6 and is used
    ! to adjust the optical depths since the k's include only lines.  This is
    ! done by multiplying swght1 by 1.6 and xka(:,2,:) by 1.6 to account for 
    ! this difference in below and above 100 hPa respectively.  Also note that
    ! the minor gas does not have a g-point dependent absorption because the
    ! o2 abosrption is a continuum effect
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/9,1/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o,io2,ih2o,io2,    &
         &                              io2,nogas,nogas,io2/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.022708_wp, 0.0_wp/)*1.6_wp
    s_wght2(:,ib) = (/8.0_wp, 0.0_wp/)
    layreffr(ib)  = 2
    laysolfr_ref(ib)  = 0

    DO ig=ibeg,ibeg+ng22-1
      jg = ig - ibeg + 1
      rayl(ig)          = rayl22
      xk_mnr(1,ig)      = 4.35e-4_wp/(350.0_wp*2.0_wp) 
      xk_mnr(2,ig)      = 4.35e-4_wp/(350.0_wp*2.0_wp) 
      xka(1:585,1,ig)   = absa22(:,jg)
      xka(1:235,2,ig)   = absb22(:,jg)*1.6_wp
      ref_self(1:10,ig) = selfref22(:,jg)
      ref_frgn(1:3 ,ig) = forref22(:,jg)
      ref_sflx(1:9,ig)  = sflx_ref22(jg,:)
    END DO
    !
    !  Prepare band 23:  8050-12850 cm-1 (low - h2o; high - nothing)
    ! 
    !  Average Giver et al. correction factor for this band is 1.029
    !  and is implemented by multiplying the abosrption coefficients
    ! 
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)    = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)   = (/1,1/)
    igas(:,:,ib) = RESHAPE ( SOURCE = (/ih2o,nogas,ih2o,nogas,     &
         &                              nogas,nogas,nogas,nogas/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -6
    laysolfr_ref(ib)  = 0

    DO ig=ibeg,ibeg+ng23-1
      jg = ig - ibeg + 1
      rayl(ig)          = rayl23(jg)
      xka(1:65,1,ig)    = absa23(:,jg) * 1.029_wp
      ref_self(1:10,ig) = selfref23(:,jg)
      ref_frgn(1:3 ,ig) = forref23(:,jg)
      ref_sflx(1,ig)    = sflx_ref23(jg)
    END DO
    !
    !  Prepare band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)
    ! 
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)    = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)   = (/9,1/)
    igas(:,:,ib) = RESHAPE ( SOURCE = (/ih2o,io2,ih2o,io3,     &
         &                              io2,nogas,nogas,io3/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.124692_wp, 0.0_wp/)
    s_wght2(:,ib) = (/8.0_wp, 0.0_wp/)
    layreffr(ib)  = 1
    laysolfr_ref(ib)  = 0

    DO ig=ibeg,ibeg+ng24-1
      jg = ig - ibeg + 1
      rayl_24_js(:,ig)  = rayl24a(jg,:)
      rayl(ig)          = rayl24b(jg)
      xk_mnr(1,ig)      = abs_m24a(jg)
      xk_mnr(2,ig)      = abs_m24b(jg)
      xka(1:585,1,ig)   = absa24(:,jg)
      xka(1:235,2,ig)   = absb24(:,jg)
      ref_self(1:10,ig) = selfref24(:,jg)
      ref_frgn(1:3 ,ig) = forref24(:,jg)
      ref_sflx(1:9,ig)  = sflx_ref24(jg,:)
    END DO
    !
    !  Prepare band 25: 16000-22650 cm-1 (low - h2o; high - nothing)
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)    = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)   = (/1,1/)
    igas(:,:,ib) = RESHAPE ( SOURCE = (/ih2o,nogas,nogas,io3,  &
         &                              nogas,nogas,nogas,io3/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -2
    laysolfr_ref(ib)  = 0

    DO ig=ibeg,ibeg+ng25-1
      jg = ig - ibeg + 1
      rayl(ig)        = rayl25(jg)
      xk_mnr(1,ig)    = abs_m25a(jg)
      xk_mnr(2,ig)    = abs_m25b(jg)
      xka(1:65,1,ig)  = absa25(:,jg)
      ref_sflx(1,ig)  = sflx_ref25(jg)
    END DO
    !
    !  Prepare band 26: 22650-29000 cm-1 (low - nothing; high - nothing)
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/1,1/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/nogas,nogas,nogas,nogas,  &
         &                              nogas,nogas,nogas,nogas/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -1
    laysolfr_ref(ib)  = 0

    DO ig=ibeg,ibeg+ng26-1
      jg = ig - ibeg + 1
      rayl(ig)        = rayl26(jg)
      ref_sflx(1,ig)  = sflx_ref26(jg)
    END DO
    !
    !  Prepare band 27:  29000-38000 cm-1 (low - o3; high - o3)
    ! 
    ! Kurucz solar source function
    ! The values in sfluxref were obtained using the "low resolution"
    ! version of the Kurucz solar source function.  For unknown reasons,
    ! the total irradiance in this band differs from the corresponding
    ! total in the "high-resolution" version of the Kurucz function.
    ! Therefore, these values are scaled below by the factor SCALEKUR.
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/1,1/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/io3,nogas,nogas,nogas,  &
         &                              io3,nogas,nogas,nogas/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -32
    laysolfr_ref(ib)  = klev

    DO ig=ibeg,ibeg+ng27-1
      jg = ig - ibeg + 1
      rayl(ig)        = rayl27(jg)
      xka(1:65,1,ig)  = absa27(:,jg)
      xka(1:235,2,ig) = absb27(:,jg)
      ref_sflx(1,ig)  = sflx_ref27(jg) * 50.15_wp/48.37_wp
    END DO
    !
    !  Prepare band 28:   38000-50000 cm-1 (low - o3,o2; high - o3,o2)
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/9,5/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/io3,io2,nogas,nogas,  &
         &                              io3,io2,nogas,nogas/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/6.67029e-07_wp, 6.67029e-07_wp/)
    s_wght2(:,ib) = (/8.0_wp, 4.0_wp/)
    layreffr(ib)  = 58
    laysolfr_ref(ib)  = klev

    DO ig=ibeg,ibeg+ng28-1
      jg = ig - ibeg + 1
      rayl(ig)          = rayl28
      xka(1:585,1,ig)   = absa28(:,jg)
      xka(1:1175,2,ig)  = absb28(:,jg)
      ref_sflx(1:5,ig)  = sflx_ref28(jg,:)
    END DO
    !
    !  Prepare band 29:   820-2600 cm-1 (low - h2o; high - co2)
    !
    ibeg  = ngs(ib) + 1
    iband = ngb(ibeg)
    ib = iband - 15
    nsp(:,ib)     = (/nspa(iband),nspb(iband)/)
    ioff(:,ib)    = (/1,1/)
    igas(:,:,ib)  = RESHAPE ( SOURCE = (/ih2o,nogas,ih2o,ico2,  &
         &                             ico2,nogas,nogas,ih2o/), &
         &                  SHAPE = (/4,2/) )
    s_wght1(:,ib) = (/0.0_wp, 0.0_wp/)
    s_wght2(:,ib) = (/0.0_wp, 0.0_wp/)
    layreffr(ib)  = -49
    laysolfr_ref(ib)  = klev

    DO ig=ibeg,ibeg+ng29-1
      jg = ig - ibeg + 1
      rayl(ig)          = rayl29
      xka(1:65,1,ig)    = absa29(:,jg)
      xka(1:235,2,ig)   = absb29(:,jg)
      xk_mnr(1,ig)      = abs_m29a(jg)
      xk_mnr(2,ig)      = abs_m29b(jg)
      ref_self(1:10,ig) = selfref29(:,jg)
      ref_frgn(1:4 ,ig) = forref29(:,jg)
      ref_sflx(1,ig)    = sflx_ref29(jg)
    END DO

    initialized = .TRUE.
  END SUBROUTINE setup_taumol

END MODULE mo_srtm_gas_optics
