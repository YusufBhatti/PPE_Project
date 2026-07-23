!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide coefficients for use in longwave solver
!!
!! @remarks
!!   For a given atmosphere, calculate the indices and fractions related to the 
!!   pressure and temperature interpolations. Also calculate the values of the 
!!   integrated Planck functions for each band at the level and layer temperatures.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2010-08)
!!         Robert Pincus, U. Colorado, visiting MPI-M, Hamburg (2011-07)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard) 
!!   code previously provided by AER and copyrighed by them.  The authors of the
!!   original AER code are: Eli J. Mlawer, Jennifer S. Delamere, Michael J. 
!!   Iacono and Shepard A. Clough with acknowledgments to Steven J. Taubman, 
!!   Karen Cady-Pereira, Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert 
!!   Bergstrom. The rewrites were designed to better interface with the structure
!!   of the ICON family of models and elements of the ICON programming standard.
!!
!! @par Copyright
!!   The AER copyright
!!   on the original code is as follows: Copyright 2002-2009, Atmospheric and
!!   Environmental Research, Inc. (AER). This software may be used, copied, or
!!   redistributed as long as it is not sold and this copyright notice is
!!   reproduced on each copy made.  This model is provided as is without any
!!   express or implied warranties. (http://www.rtweb.aer.com/)               
!! 
!
MODULE mo_rrtm_coeffs

  USE mo_kind, ONLY : wp
  USE mo_rrtm_params, ONLY : preflog, tref
  USE rrlw_planck,    ONLY : chi_mls
  IMPLICIT NONE
  
  REAL(wp), PARAMETER :: stpfac  = 296._wp/1013._wp

CONTAINS
! --------------------------------------------------------------------------------------------
  SUBROUTINE lrtm_coeffs( kproma, kbdim                              , &
       & klev           ,play          ,tlay          ,coldry        , &
       & wkl            ,wbroad        ,laytrop       ,jp            , &
       & jt             ,jt1                                         , &
       &                 colh2o        ,colco2        ,colo3         , &
       & coln2o         ,colco         ,colch4        ,colo2         , &
       & colbrd         ,fac00         ,fac01         ,fac10         , &
       & fac11          ,rat_h2oco2    ,rat_h2oco2_1  ,rat_h2oo3     , &
       & rat_h2oo3_1    ,rat_h2on2o    ,rat_h2on2o_1  ,rat_h2och4    , &
       & rat_h2och4_1   ,rat_n2oco2    ,rat_n2oco2_1  ,rat_o3co2     , &
       & rat_o3co2_1    ,selffac       ,selffrac      ,indself       , &
       & forfac         ,forfrac       ,indfor        ,minorfrac     , &
       & scaleminor     ,scaleminorn2  ,indminor)

    INTEGER, INTENT(in) ::  &
         kproma,            & ! number of columns
         kbdim,             & ! maximum number of column as first dim is declared in calling (sub)prog.
         klev                 ! total number of layers
    REAL(wp), INTENT(in) :: &
         play(kbdim,klev)          , & ! layer pressures (mb) 
         tlay(kbdim,klev)          , & ! layer temperatures (K)
         coldry(kbdim,klev)        , & ! dry air column density (mol/cm2)
         wbroad(kbdim,klev)        , & ! broadening gas column density (mol/cm2)
         wkl(:,:,:)             !< molecular amounts (mol/cm-2) (mxmol,klev)
    !
    ! Output Dimensions kproma, klev unless otherwise specified
    !
    INTEGER, INTENT(out) :: &
         laytrop(kbdim)         , & !< tropopause layer index
         jp(kbdim,klev)            , & ! 
         jt(kbdim,klev)            , & !
         jt1(kbdim,klev)           , & !
         indself(kbdim,klev)       , & !
         indfor(kbdim,klev)        , & !
         indminor(kbdim,klev)          !
    REAL(wp), INTENT(out) :: &
         colh2o(kbdim,klev)         , & !< column amount (h2o)
         colco2(kbdim,klev)         , & !< column amount (co2)
         colo3(kbdim,klev)          , & !< column amount (o3)
         coln2o(kbdim,klev)         , & !< column amount (n2o)
         colco(kbdim,klev)          , & !< column amount (co)
         colch4(kbdim,klev)         , & !< column amount (ch4)
         colo2(kbdim,klev)          , & !< column amount (o2)
         colbrd(kbdim,klev)         , & !< column amount (broadening gases)
         selffac(kbdim,klev)        , & !<
         selffrac(kbdim,klev)       , & !<
         forfac(kbdim,klev)         , & !<
         forfrac(kbdim,klev)        , & !<
         fac00(kbdim,klev)          , & !<
         fac01(kbdim,klev)          , &
         fac10(kbdim,klev)          , & 
         fac11(kbdim,klev)          , &
         minorfrac(kbdim,klev)      , &
         scaleminor(kbdim,klev)     , &
         scaleminorn2(kbdim,klev)   , &
         rat_h2oco2(kbdim,klev)     , &
         rat_h2oco2_1(kbdim,klev)   , &
         rat_h2oo3(kbdim,klev)      , &
         rat_h2oo3_1(kbdim,klev)    , & 
         rat_h2on2o(kbdim,klev)     , &
         rat_h2on2o_1(kbdim,klev)   , &
         rat_h2och4(kbdim,klev)     , &
         rat_h2och4_1(kbdim,klev)   , &
         rat_n2oco2(kbdim,klev)     , &
         rat_n2oco2_1(kbdim,klev)   , &
         rat_o3co2(kbdim,klev)      , &
         rat_o3co2_1(kbdim,klev)
    
    INTEGER  :: jk
    REAL(wp) :: colmol(kbdim,klev), factor(kbdim,klev) 

    ! ------------------------------------------------
    CALL srtm_coeffs( kproma, kbdim,  klev                          , &
       & play          ,tlay          ,coldry        ,wkl           , &
       & laytrop       ,jp            ,jt            ,jt1           , &
       & colch4        ,colco2        ,                               &
       & colh2o        ,colmol        ,coln2o        ,colo2         , &
       & colo3         ,fac00         ,fac01         ,fac10         , &
       & fac11         ,selffac       ,selffrac      ,indself       , &
       & forfac        ,forfrac       ,indfor)

    
    colbrd(1:kproma,1:klev) = 1.e-20_wp * wbroad(1:kproma,1:klev)
    colco (1:kproma,1:klev) = MERGE(1.e-20_wp * wkl(1:kproma,5,1:klev),  &
                                    1.e-32_wp * coldry(1:kproma,1:klev), &
                                    wkl(1:kproma,5,1:klev) > 0._wp)
    
    !
    ! Water vapor continuum broadening factors are used differently in LW and SW? 
    !
     forfac(1:kproma,1:klev) =  forfac(1:kproma,1:klev) * colh2o(1:kproma,1:klev)
    selffac(1:kproma,1:klev) = selffac(1:kproma,1:klev) * colh2o(1:kproma,1:klev)
    
    !
    !  Setup reference ratio to be used in calculation of binary species parameter.
    !
    DO jk = 1, klev
      rat_h2oco2  (1:kproma, jk) = chi_mls(1,jp(1:kproma, jk)  )/chi_mls(2,jp(1:kproma, jk)  )
      rat_h2oco2_1(1:kproma, jk) = chi_mls(1,jp(1:kproma, jk)+1)/chi_mls(2,jp(1:kproma, jk)+1)
      !
      ! Needed only in lower atmos (plog > 4.56_wp) 
      !
      rat_h2oo3   (1:kproma, jk) = chi_mls(1,jp(1:kproma, jk)  )/chi_mls(3,jp(1:kproma, jk)  )
      rat_h2oo3_1 (1:kproma, jk) = chi_mls(1,jp(1:kproma, jk)+1)/chi_mls(3,jp(1:kproma, jk)+1)
      rat_h2on2o  (1:kproma, jk) = chi_mls(1,jp(1:kproma, jk)  )/chi_mls(4,jp(1:kproma, jk)  )
      rat_h2on2o_1(1:kproma, jk) = chi_mls(1,jp(1:kproma, jk)+1)/chi_mls(4,jp(1:kproma, jk)+1)
      rat_h2och4  (1:kproma, jk) = chi_mls(1,jp(1:kproma, jk)  )/chi_mls(6,jp(1:kproma, jk)  )
      rat_h2och4_1(1:kproma, jk) = chi_mls(1,jp(1:kproma, jk)+1)/chi_mls(6,jp(1:kproma, jk)+1)
      rat_n2oco2  (1:kproma, jk) = chi_mls(4,jp(1:kproma, jk)  )/chi_mls(2,jp(1:kproma, jk)  )
      rat_n2oco2_1(1:kproma, jk) = chi_mls(4,jp(1:kproma, jk)+1)/chi_mls(2,jp(1:kproma, jk)+1)
      !
      ! Needed only in upper atmos (plog <= 4.56_wp) 
      !
      rat_o3co2   (1:kproma, jk) = chi_mls(3,jp(1:kproma, jk)  )/chi_mls(2,jp(1:kproma, jk)  )
      rat_o3co2_1 (1:kproma, jk) = chi_mls(3,jp(1:kproma, jk)+1)/chi_mls(2,jp(1:kproma, jk)+1)         
    END DO 

    !
    !  Set up factors needed to separately include the minor gases
    !  in the calculation of absorption coefficient
    !
    scaleminor  (1:kproma,1:klev) = play(1:kproma,1:klev)/tlay(1:kproma,1:klev)
    scaleminorn2(1:kproma,1:klev) = scaleminor(1:kproma,1:klev) * &
         &                             (wbroad(1:kproma,1:klev)/(coldry(1:kproma,1:klev)+wkl(1:kproma,1,1:klev)))
    factor(1:kproma,1:klev) = (tlay(1:kproma,1:klev)-180.8_wp)/7.2_wp
    indminor(1:kproma,1:klev)     = MIN(18, MAX(1, INT(factor(1:kproma,1:klev))))
    minorfrac(1:kproma,1:klev)    = (tlay(1:kproma,1:klev)-180.8_wp)/7.2_wp - FLOAT(indminor(1:kproma,1:klev))

  END SUBROUTINE lrtm_coeffs    

! --------------------------------------------------------------------------------------------

  SUBROUTINE srtm_coeffs(kproma       , kbdim        ,klev          , &
       & play          ,tlay          ,coldry        ,wkl           , &
       & laytrop       ,jp            ,jt            ,jt1           , &
       & colch4        ,colco2        ,                               &
       & colh2o        ,colmol        ,coln2o        ,colo2         , &
       & colo3         ,fac00         ,fac01         ,fac10         , &
       & fac11         ,selffac       ,selffrac      ,indself       , &
       & forfac        ,forfrac       ,indfor)

    INTEGER, INTENT(in) ::  &
         kproma,            & ! number of columns
         kbdim,             & ! maximum number of col. as declared in calling (sub)programs
         klev              ! total number of layers
    REAL(wp), INTENT(in) :: &
         play(kbdim,klev)         , & ! layer pressures (mb) 
         tlay(kbdim,klev)         , & ! layer temperatures (K)
         coldry(kbdim,klev)        , & ! dry air column density (mol/cm2)
         wkl(:,:,:)             !< molecular amounts (mol/cm-2) (mxmol,klev)
    !
    ! Output Dimensions kproma, klev unless otherwise specified
    !
    INTEGER, INTENT(out) :: &
         laytrop(kbdim)         , & !< tropopause layer index
         jp(kbdim,klev)            , & ! 
         jt(kbdim,klev)            , & !
         jt1(kbdim,klev)           , & !
         indself(kbdim,klev)       , & !
         indfor(kbdim,klev)        
    REAL(wp), INTENT(out) :: &
         colh2o(kbdim,klev)         , & !< column amount (h2o)
         colco2(kbdim,klev)         , & !< column amount (co2)
         colo3(kbdim,klev)          , & !< column amount (o3)
         coln2o(kbdim,klev)         , & !< column amount (n2o)
         colch4(kbdim,klev)         , & !< column amount (ch4)
         colo2(kbdim,klev)          , & !< column amount (o2)
         colmol(kbdim,klev)         , &
         selffac(kbdim,klev)        , & !<
         selffrac(kbdim,klev)       , & !<
         forfac(kbdim,klev)         , & !<
         forfrac(kbdim,klev)        , & !<
         fac00(kbdim,klev)          , & !<
         fac01(kbdim,klev)          , &
         fac10(kbdim,klev)          , & 
         fac11(kbdim,klev)          

    INTEGER :: jp1(kbdim,klev), jk
    REAL(wp) :: &
         plog  (kbdim,klev), fp      (kbdim,klev), &
         ft    (kbdim,klev), ft1     (kbdim,klev), &
         water (kbdim,klev), scalefac(kbdim,klev), &
         compfp(kbdim,klev), factor  (kbdim,klev) 
    ! -------------------------------------------------------------------------
    !
    !  Find the two reference pressures on either side of the
    !  layer pressure.  Store them in JP and JP1.  Store in FP the
    !  fraction of the difference (in ln(pressure)) between these
    !  two values that the layer pressure lies.
    !
    plog(1:kproma,1:klev) = LOG(play(1:kproma,1:klev))
    jp (1:kproma,1:klev)  = MIN(58,MAX(1,INT(36._wp - 5*(plog(1:kproma,1:klev)+0.04_wp))))
    jp1(1:kproma,1:klev)  = jp(1:kproma,1:klev) + 1
    do jk = 1, klev
      fp(1:kproma,jk)   = 5._wp *(preflog(jp(1:kproma,jk)) - plog(1:kproma,jk))
    end do 
    !
    !  Determine, for each reference pressure (JP and JP1), which
    !  reference temperature (these are different for each  
    !  reference pressure) is nearest the layer temperature but does
    !  not exceed it.  Store these indices in JT and JT1, resp.
    !  Store in FT (resp. FT1) the fraction of the way between JT
    !  (JT1) and the next highest reference temperature that the 
    !  layer temperature falls.
    !
    DO jk = 1, klev
      jt(1:kproma,jk)                                           &
                   = MIN(4,MAX(1,INT(3._wp + (tlay(1:kproma,jk) &
                                              - tref(jp (1:kproma,jk)))/15._wp)))
      jt1(1:kproma,jk)                                           &
                   = MIN(4,MAX(1,INT(3._wp + (tlay(1:kproma,jk) & 
                                              - tref(jp1(1:kproma,jk)))/15._wp)))
    END DO 
    DO jk = 1, klev
      ft(1:kproma,jk)     = ((tlay(1:kproma,jk)-tref(jp (1:kproma,jk)))/15._wp) &
                            - float(jt (1:kproma,jk)-3)
      ft1(1:kproma,jk)    = ((tlay(1:kproma,jk)-tref(jp1(1:kproma,jk)))/15._wp) &
                            - float(jt1(1:kproma,jk)-3)
    END DO 
    water(1:kproma,1:klev)    = wkl(1:kproma,1,1:klev)/coldry(1:kproma,1:klev)
    scalefac(1:kproma,1:klev) = play(1:kproma,1:klev) * stpfac / tlay(1:kproma,1:klev)

    !
    !  We have now isolated the layer ln pressure and temperature,
    !  between two reference pressures and two reference temperatures 
    !  (for each reference pressure).  We multiply the pressure 
    !  fraction FP with the appropriate temperature fractions to get 
    !  the factors that will be needed for the interpolation that yields
    !  the optical depths (performed in routines TAUGBn for band n).`
    !
    compfp(1:kproma,1:klev) = 1. - fp(1:kproma,1:klev)
    fac10(1:kproma,1:klev) = compfp(1:kproma,1:klev) * ft(1:kproma,1:klev)
    fac00(1:kproma,1:klev) = compfp(1:kproma,1:klev) * (1._wp - ft(1:kproma,1:klev))
    fac11(1:kproma,1:klev) = fp(1:kproma,1:klev) * ft1(1:kproma,1:klev)
    fac01(1:kproma,1:klev) = fp(1:kproma,1:klev) * (1._wp - ft1(1:kproma,1:klev))
    
    ! Tropopause defined in terms of pressure (~100 hPa)
    !   We're looking for the first layer (counted from the bottom) at which the pressure reaches
    !   or falls below this value
    ! 
    laytrop(1:kproma) = COUNT(plog(1:kproma,1:klev) > 4.56_wp, DIM = 2) 

    !
    !  Calculate needed column amounts.
    !    Only a few ratios are used in the upper atmosphere but masking may be less efficient
    !
    colh2o(1:kproma,1:klev) = MERGE(1.e-20_wp * wkl(1:kproma,1,1:klev),  &
                                    1.e-32_wp * coldry(1:kproma,1:klev), &
                                    wkl(1:kproma,1,1:klev) > 0._wp)
    colco2(1:kproma,1:klev) = MERGE(1.e-20_wp * wkl(1:kproma,2,1:klev),  &
                                    1.e-32_wp * coldry(1:kproma,1:klev), &
                                    wkl(1:kproma,2,1:klev) > 0._wp)
    colo3 (1:kproma,1:klev) = MERGE(1.e-20_wp * wkl(1:kproma,3,1:klev),  &
                                    1.e-32_wp * coldry(1:kproma,1:klev), &
                                    wkl(1:kproma,3,1:klev) > 0._wp)
    coln2o(1:kproma,1:klev) = MERGE(1.e-20_wp * wkl(1:kproma,4,1:klev),  &
                                    1.e-32_wp * coldry(1:kproma,1:klev), &
                                    wkl(1:kproma,4,1:klev) > 0._wp)
    colch4(1:kproma,1:klev) = MERGE(1.e-20_wp * wkl(1:kproma,6,1:klev),  &
                                    1.e-32_wp * coldry(1:kproma,1:klev), &
                                    wkl(1:kproma,6,1:klev) > 0._wp)
    colo2 (1:kproma,1:klev) = MERGE(1.e-20_wp * wkl(1:kproma,7,1:klev),  &
                                    1.e-32_wp * coldry(1:kproma,1:klev), &
                                    wkl(1:kproma,7,1:klev) > 0._wp)
    colmol(1:kproma,1:klev) = 1.e-20_wp * coldry(1:kproma,1:klev) + colh2o(1:kproma,1:klev)


    ! ------------------------------------------
    ! Interpolation coefficients 
    !
    forfac(1:kproma,1:klev) = scalefac(1:kproma,1:klev) / (1._wp+water(1:kproma,1:klev))
    !
    !  Set up factors needed to separately include the water vapor
    !  self-continuum in the calculation of absorption coefficient.
    !
    selffac(1:kproma,1:klev)  = water(1:kproma,1:klev) * forfac(1:kproma,1:klev)
    
    !
    !  If the pressure is less than ~100mb, perform a different set of species
    !  interpolations.
    !
    factor(1:kproma,1:klev) = (332.0_wp-tlay(1:kproma,1:klev))/36.0_wp
    indfor(1:kproma,1:klev) =                 &
      MERGE(3,                                &
            MIN(2, MAX(1, INT(factor(1:kproma,1:klev)))), & 
            plog(1:kproma,1:klev) <= 4.56_wp)
            
    forfrac(1:kproma,1:klev) =                                  & 
      MERGE((tlay(1:kproma,1:klev)-188.0_wp)/36.0_wp - 1.0_wp, &
            factor(1:kproma,1:klev) - FLOAT(indfor(1:kproma,1:klev)),      &
            plog(1:kproma,1:klev) <= 4.56_wp)

    ! In RRTMG code, this calculation is done only in the lower atmosphere (plog > 4.56) 
    ! 
    factor(1:kproma,1:klev)  = (tlay(1:kproma,1:klev)-188.0_wp)/7.2_wp
    indself (1:kproma,1:klev)  = MIN(9, MAX(1, INT(factor(1:kproma,1:klev))-7))
    selffrac(1:kproma,1:klev) = factor(1:kproma,1:klev) - float(indself(1:kproma,1:klev) + 7)
  END SUBROUTINE srtm_coeffs
END MODULE mo_rrtm_coeffs

