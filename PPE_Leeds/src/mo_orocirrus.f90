!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! @brief mo_orocirrus computes the updraft velocity and the additional 
!! cirrus cloud coverage due to gravity waves.              
!!
!! This routine computes the updraft velocity (pvervx_2d) over mountains due to
!! the gravity waves. The newly calculated updraft velocity is used 
!! for the calculation of the homogeneously frozen ice crystals.
!! Consequently also the cloud coverage is updated.
!! Attention: 
!! You can use this module only for Bernd's cirrus cloud scheme 
!! (nic_cirrus=2)
!!
!! @par References
!!     Lott und Miller 97, Eq.(11)                                (updraft velocity)
!!     Kaercher and Lohmann, 2002, JGR.        (freezing threshold saturation ratio)    
!!     Murphy and Koop, 2005, G.J.R.Met.Soc.                 (water pressure of ice)
!!     Joos et al.,(2008),JGR., 113, D18205, doi:10.1029/2007JD009605
!!     Joos et al.,(2010),JGR., 115, D19129, doi:10.1029/2010JD013824
!! @author 
!! <ol>
!! <li> H.Joos      ETH Zuerich  2008
!! </ol>
!! @par Revision History
!! <ol>
!!<li>G.Frontoso    ETH Zuerich  2013             : implemented into echam6-ham2.1
!!
!! @par This module is used by
!! cloud_cdnc_icnc
!!
!! @par Responsible coder
!! grazia.frontoso@env.ethz.ch
!!
!! @par Copyright
!! 2009 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ECHAM is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!! violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!! copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!! an according license agreement with MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_orocirrus

#ifdef HAMMOZ

  USE mo_kind,           ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: orocirrus_w, orocirrus_cc

  INTEGER,  ALLOCATABLE :: ktest(:)        !< flag to indicate active points (oro scheme)

  REAL(dp), ALLOCATABLE :: pw_gwd_kpr(:,:) !< max updraft velocity over mountains
  REAL(dp), ALLOCATABLE :: pampl_gwd(:,:)  !< amplitude of the gravity waves
  REAL(dp), ALLOCATABLE :: pl_z(:)         !< hortizontal wavelength

  CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE orocirrus_w(kproma, kbdim, klev, ktdia, paphm1, papm1,   &
                         pgeom1, ptm1, pqm1, pum1, pvm1, pstd, psig,  &
                         pmea, pgam, pthe, ppic, pval, pg_rcp,        &
                         prho_rcp, pvervel, pqsi_2d, pvervx_2d)

  !---------------------------------------------------------------------------
  !>
  !! @brief Calculate the updraft velocity taking into account the gravity waves
  !!

    USE mo_ssodrag
    USE mo_ssortns,        ONLY: orosetup, gwstress
    USE mo_physical_constants, ONLY: grav,cpd,rd,rv
    USE mo_math_constants,     ONLY: pi
 
    IMPLICIT NONE

  INTEGER, INTENT(in) :: kproma        !< block size
  INTEGER, INTENT(in) :: kbdim         !< max. block size
  INTEGER, INTENT(in) :: klev          !< number of vertical levels
  INTEGER, INTENT(in) :: ktdia         !< highest vertical level for "diagnostics" (use??)-careful with this if /=1!

  REAL(dp), INTENT(in) :: pg_rcp       !< reciprocal of gravity acceleration

  REAL(dp), INTENT(in) :: pstd(kbdim)  !< SSO standard deviation (m)
  REAL(dp), INTENT(in) :: psig(kbdim)  !< SSO slope
  REAL(dp), INTENT(in) :: pmea(kbdim)  !< mean Orography (m)
  REAL(dp), INTENT(in) :: pgam(kbdim)  !< SSO Anisotropy
  REAL(dp), INTENT(in) :: pthe(kbdim)  !< SSO angle
  REAL(dp), INTENT(in) :: ppic(kbdim)  !< SSO peaks elevation (m)
  REAL(dp), INTENT(in) :: pval(kbdim)  !< SSO valleys elevation (m) 

  REAL(dp), INTENT(in) :: ptm1(kbdim,klev)     !< temperature (t-dt)
  REAL(dp), INTENT(in) :: pum1(kbdim,klev)     !< zonal wind (t-dt)
  REAL(dp), INTENT(in) :: pvm1(kbdim,klev)     !< meridional wind (t-dt)
  REAL(dp), INTENT(in) :: papm1(kbdim,klev)    !< full level pressure (t-dt)
  REAL(dp), INTENT(in) :: pgeom1(kbdim,klev)   !< geopotential above surface (t-dt)
  REAL(dp), INTENT(in) :: pqm1(kbdim,klev)     !< specific humidity   
  REAL(dp), INTENT(in) :: prho_rcp(kbdim,klev) !< inverse air density
  REAL(dp), INTENT(in) :: pvervel(kbdim,klev)  !< large scale vertical velocity [Pa s-1]
  REAL(dp), INTENT(in) :: pqsi_2d(kbdim,klev)  !< saturation specific humidity w.r.t. ice [kg/kg]

  REAL(dp), INTENT(in) :: paphm1(kbdim,klev+1) !< half level pressure (t-dt)

  REAL(dp), INTENT(inout) :: pvervx_2d(kbdim,klev) !< updraft velocity !SF inout is IMPORTANT!!!

!   Local scalars/arrays

    INTEGER  :: ji,jl,jk 
    INTEGER  :: kgwd           !< total points where orocirrus is active
    INTEGER  :: kdx(kbdim)     !< physical location of an active point
    INTEGER  :: ikcrith(kbdim) !< security value for top of low level flow
    INTEGER  :: ikcrit(kbdim)  !< critical level
    INTEGER  :: ikenvh(kbdim)  !< top of blocked flow layer
    INTEGER  :: iknu(kbdim)    !< layer that sees mountain peaks
    INTEGER  :: iknu2(kbdim)   !< layer that sees mountain peaks above mountain mean
    INTEGER  :: icrit(kbdim)  

    REAL(dp)  :: zsqr,zalfa,zdz2n,zalpha,zriw,zdel,zb,zvt1,zvt2
    REAL(dp)  :: zdelp, zdelpt, zeps

    REAL(dp)  :: zdmod(kbdim)               !< Norme of pdi
    REAL(dp)  :: zulow(kbdim), zvlow(kbdim) !< Low-level wind

    REAL(dp)  :: zd1(kbdim),zd2(kbdim),                &
                 znu(kbdim),znorm(kbdim),zoro(kbdim),  &
                 zdz2 (kbdim,klev), zdep_b(kbdim)

    REAL(dp)  :: zzdep(kbdim,klev), zw_gwd(kbdim,klev),&
                 zvpf(kbdim,klev)

    REAL(dp)  :: zstab(kbdim,klev+1) !< Brunt-Vaisala freq. at 1/2 layers
    REAL(dp)  :: zri(kbdim,klev+1)   !< Background Richardson Number, Wind shear measured along GW stress
    REAL(dp)  :: zrho_h(kbdim,klev+1)!< Air density [kg/m3] at 1/2 layers
    REAL(dp)  :: zvph(kbdim,klev+1)  !< Wind in  plan of GW stress at 1/2 layers
    REAL(dp)  :: zpsi(kbdim,klev+1)  !< Angle between low level wind and SS0 main axis
    REAL(dp)  :: ztau(kbdim,klev+1)  !< Surface stress due to gravity waves
    REAL(dp)  :: ptau(kbdim,klev+1)  !< Surface stress due to gravity waves
    
  REAL(dp) :: zi,zstep_long,tempgw,wsin,zice,Scr,Scr0,rhi,rhival,lambda

  REAL(dp) :: zsusati(kbdim,klev),zwnucl(kbdim,klev)

!< Kaercher and Lohmann, 2002, JGR

  REAL(dp), PARAMETER :: zscr1=2.583_dp
  REAL(dp), PARAMETER :: zscr2=207.83_dp
  REAL(dp), PARAMETER :: zscr3=2.418_dp
  REAL(dp), PARAMETER :: zscr4=245.68_dp

!< Murphy and Koop, 2005, G.J.R.Met.Soc.

  REAL(dp), PARAMETER :: zzice1=9.550426_dp
  REAL(dp), PARAMETER :: zzice2=5723.265_dp
  REAL(dp), PARAMETER :: zzice3=3.53068_dp
  REAL(dp), PARAMETER :: zzice4=0.00728332_dp

! Executable statements

  IF (.NOT. ALLOCATED(pampl_gwd))   ALLOCATE(pampl_gwd(kbdim,klev))
  IF (.NOT. ALLOCATED(pw_gwd_kpr))  ALLOCATE(pw_gwd_kpr(kbdim,klev))
  IF (.NOT. ALLOCATED(pl_z))        ALLOCATE(pl_z(kbdim))
  IF (.NOT. ALLOCATED(ktest))       ALLOCATE(ktest(kbdim))

!---Initialization

 pampl_gwd(1:kproma,:) =0._dp
 pw_gwd_kpr(1:kproma,:)=0._dp
 pl_z(1:kproma )       =0._dp
 ktest(1:kproma )      =0._dp

!---Computational constants:

  zeps=EPSILON(1._dp)
  zstep_long  = 1000._dp

  CALL sugwd(klev)

  !
  !*         1.    orographic gravity wave drag
  !                -----------------------------

  !  SELECTION  POINTS WHERE THE SCHEME IS ACTIVE

  kdx(:) = 0
  kgwd=0
  DO jl=1,kproma
     ktest(jl)=0
     IF (((ppic(jl)-pmea(jl)) > gpicmea).AND.(pstd(jl) > gstd)) THEN
        ktest(jl)=1
        kgwd=kgwd+1
        kdx(kgwd)=jl
     ENDIF
  ENDDO


  !
  !*         2.    call orosetup
  !                --------------

  CALL orosetup(kproma, kbdim, klev, kgwd, kdx,                       &
                ikcrit, ikcrith, icrit, ikenvh, iknu, iknu2,          &
                paphm1, papm1, pum1, pvm1, ptm1, pgeom1,              &
                zrho_h, zri, zstab, ztau, zvph, zpsi, zzdep,          & 
                zulow, zvlow, pthe, pgam, pmea, ppic, pval,           &
                znu, zd1, zd2, zdmod)


!
  !*         3.    calculate znorm
  !                ------------------------

!IBM* ASSERT(NODEPS)
  DO 2110 ji=1,kgwd
     jl = kdx(ji)
     znorm(jl)=MAX(SQRT(zulow(jl)**2+zvlow(jl)**2),gvsec)
2110 END DO

  !
  !*         4.    calculate zvpf
  !                ------------------------
  !
  !  ************ projet flow in plane of lowlevel stress *************
  !  ************ Find critical levels...                 *************
  !

  DO 213 jk=1,klev
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
     DO 212 ji=1,kgwd
        jl = kdx(ji)
        zvt1       = zulow(jl)*pum1(jl,jk)+zvlow(jl)*pvm1(jl,jk)
        zvt2       =-zvlow(jl)*pum1(jl,jk)+zulow(jl)*pvm1(jl,jk)
        zvpf(jl,jk)=(zvt1*zd1(jl)+zvt2*zd2(jl))/(znorm(jl)*zdmod(jl))
212  END DO
213 END DO

  !
  !*         5.    calculate zdep_b
  !                --------------

   DO 253 ji=1,kgwd
          jl = kdx(ji)
          zdep_b(jl)=MIN(SQRT(zzdep(jl,ikenvh(jl))),1._dp)
           IF (ikenvh(jl)== klev .OR. ikenvh(jl)== klev+1) THEN
              zdep_b(jl)=1._dp
           ENDIF
253  END DO


  !
  !*         6.    call gwstress
  !                ------------------------

  CALL gwstress(kproma, kbdim, klev, kgwd, kdx,       &
                ikenvh, zrho_h, zstab, zvph,          &
                pstd, psig, ppic, pval, ptau, zdmod)

 
  !
  !*         7.    calculate pampl_gwd 
  !                ------------------------

!CDIR NODEP
!DIR$ CONCURRENT
  DO 400 ji=1,kgwd
     jl=kdx(ji)
     zoro(jl)=psig(jl)*zdmod(jl)/4._dp/pstd(jl)
     ztau(jl,klev+1)=ptau(jl,klev+1)
     ztau(jl,ikcrith(jl))=grahilo*ptau(jl,klev+1)
400 END DO

  !
  DO 430 jk=klev+1,2,-1
     !
     !
     !*         4.1    constant shear stress until top of the
     !                 low-level breaking/trapped layer

     !
!DIR$ CONCURRENT
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
     DO 411 ji=1,kgwd
        jl=kdx(ji)
        IF(jk > ikcrith(jl)) THEN
           zdelp=paphm1(jl,jk)-paphm1(jl,klev+1)
           zdelpt=paphm1(jl,ikcrith(jl))-paphm1(jl,klev+1)
           ptau(jl,jk)=ztau(jl,klev+1)+zdelp/zdelpt*                    &
                (ztau(jl,ikcrith(jl))-ztau(jl,klev+1))
        ELSE
           ptau(jl,jk)=ztau(jl,ikcrith(jl))
        ENDIF
411  END DO

430 END DO

  !
  !*         4.15   constant shear stress until the top of the
  !                 low level flow layer.

  !
  !
  !*         4.2    wave displacement at next level.
  !
  !

  !
  !*         4.4    wave richardson number, new wave displacement
  !*                and stress:  breaking evaluation and critical
  !                 level
  !

  DO 440 jk=klev,2,-1
!DIR$ CONCURRENT
     DO 441 ji=1,kgwd
        jl=kdx(ji)
        znorm(jl)=zrho_h(jl,jk)*SQRT(zstab(jl,jk))*zvph(jl,jk)
        zdz2(jl,jk)=ptau(jl,jk)/MAX(znorm(jl),gssec)/zoro(jl)
        pampl_gwd(jl,jk)=zdz2(jl,jk)
441  END DO

!CDIR NODEP
!DIR$ CONCURRENT
     DO 442 ji=1,kgwd
        jl=kdx(ji)
        IF(jk < ikcrith(jl)) THEN
           IF((ptau(jl,jk+1) < gtsec).OR.(jk <= icrit(jl))) THEN
              ptau(jl,jk)=0.0_dp
           ELSE
              zsqr=SQRT(zri(jl,jk))
              zalfa=SQRT(zstab(jl,jk)*zdz2(jl,jk))/zvph(jl,jk)
              zriw=zri(jl,jk)*(1._dp-zalfa)/(1+zalfa*zsqr)**2
              IF(zriw < grcrit) THEN
                 zdel=4._dp/zsqr/grcrit+1._dp/grcrit**2+4._dp/grcrit
                 zb=1._dp/grcrit+2._dp/zsqr
                 zalpha=0.5_dp*(-zb+SQRT(zdel))
                 zdz2n=(zvph(jl,jk)*zalpha)**2/zstab(jl,jk)
                 pampl_gwd(jl,jk)=(zvph(jl,jk)*zalpha)**2/zstab(jl,jk)
                 ptau(jl,jk)=znorm(jl)*zdz2n*zoro(jl)
              ENDIF

              ptau(jl,jk)=MIN(ptau(jl,jk),ptau(jl,jk+1))

           ENDIF
        ENDIF

442  END DO
440 END DO

  !
  !*         8.    calculate pl_z and zw_gwd
  !                ------------------------

    DO 524 jk=1,klev
     DO 523 ji=1,kgwd
            jl = kdx(ji)
            pl_z(jl) = 2._dp*MAX(ABS(pstd(jl)/psig(jl)*  &
                       COS(zpsi(jl,klev+1))),            &
                       ABS(pstd(jl)/(pgam(jl)*psig(jl))* &
                       SIN(zpsi(jl,klev+1))))*zdep_b(jl)


        IF (ABS(pl_z(jl)) > zeps) THEN
           zw_gwd(jl,jk) = 2._dp*pi/(2.0_dp*pl_z(jl))* &
                           SQRT(MAX(pampl_gwd(jl,jk),0._dp))*MAX(0._dp,zvpf(jl,jk))
        ELSE
           zw_gwd(jl,jk)=0._dp
        END IF

523  END DO
524 END DO

  !
  !*         9.    calculate pw_gwd_kpr
  !                --------------------

  DO jk=1,klev
     DO ji=1,kgwd
        jl = kdx(ji)
        pw_gwd_kpr(jl,jk)=zw_gwd(jl,jk)
     END DO
  END DO

  !
  !*         10.    calculate zsusati
  !                ------------------------

    DO jk = klev, ktdia, -1
       DO jl = 1, kproma
          zsusati(jl,jk)= MAX(pqm1(jl,jk)/pqsi_2d(jl,jk)-1._dp,0._dp)
       ENDDO
    ENDDO

  !
  !*         11.    calculate pvervx_2d
  !                ----------------

! Reduction of the maximum vertical velocity occuring in the gravity waves
! See Joos et al. 2010, section 2.2 and Eq. 12

    DO jk = klev, ktdia, -1
     DO ji=1,kgwd
            jl = kdx(ji)
             IF (pl_z(jl) < 100._dp) THEN
                zwnucl(jl,jk)=0._dp
             ELSE
                lambda=2._dp*pl_z(jl)
                Scr0=(zscr1-ptm1(jl,jk)/zscr2)*100._dp
                zwnucl(jl,jk)=0._dp
                rhival=zsusati(jl,jk)*100._dp+100._dp
                zi=1._dp
                
                IF (rhival > Scr0) THEN
                   zwnucl(jl,jk)=0._dp
                ELSE
                  DO WHILE (rhival <= Scr0 .AND. zi < lambda)
                      tempgw=MAX(100._dp,ptm1(jl,jk)-((grav/cpd*SQRT(pampl_gwd(jl,jk))) &
                             *SIN(3._dp*pi/2._dp+2._dp*pi*zi/lambda) &
                             +(grav/cpd*SQRT(pampl_gwd(jl,jk)))))
                      wsin=pw_gwd_kpr(jl,jk)*SIN(2._dp*pi*zi/lambda)
                      zice=EXP(zzice1-(zzice2/tempgw)+ &
                           (zzice3*LOG(tempgw))-(zzice4*tempgw))
                      rhi=(100._dp*papm1(jl,jk)*pqm1(jl,jk))/((rd/rv)*zice)
                      Scr=(zscr3-tempgw/zscr4)*100._dp
                      rhival=rhi
                      Scr0=Scr
                      zwnucl(jl,jk)=wsin                         !in m/s
                      zi=zi+zstep_long
                   END DO

                END IF ! rhival
             END IF ! pl_z

             zwnucl(jl,jk)=MAX(0.001_dp,zwnucl(jl,jk)) ! for safety reasons (Hanna)


             pvervx_2d(jl,jk) = -100.0_dp*pvervel(jl,jk)*prho_rcp(jl,jk)*pg_rcp &
                                 +(zwnucl(jl,jk)*100.0_dp)

      ENDDO ! jl loop
     ENDDO ! jk loop


  END  SUBROUTINE orocirrus_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE orocirrus_cc(kproma, kbdim, klev, ktdia, krow,      &
                          pqinucl, psusatix_2d, paclc)

  !---------------------------------------------------------------------------
  !>
  !! @brief Calculate the additional cloud cover due to the modified updraft velocity
  !!

    USE mo_geoloc,             ONLY: sqcst_2d  !<cos of latitude
    USE mo_echam_cloud_params, ONLY: ccwmin

    IMPLICIT NONE

  INTEGER, INTENT(in)  :: kproma       !< block size
  INTEGER, INTENT(in)  :: kbdim        !< max. block size
  INTEGER, INTENT(in)  :: klev         !< number of vertical levels
  INTEGER, INTENT(in)  :: ktdia        !< highest vertical level for "diagnostics" (use??)-careful with this if /=1!
  INTEGER, INTENT(in)  :: krow         !< geographic block number

  REAL(dp), INTENT(in) :: pqinucl(kbdim,klev)     !< mixing ratio of newly nucleated IC [kg/kg]
  REAL(dp), INTENT(in) :: psusatix_2d(kbdim,klev) !< supersat. with respect to ice

  REAL(dp), INTENT(inout) :: paclc(kbdim,klev)     !< cloud cover

! Local scalars/arrays

  INTEGER ::  jl,jk

  REAL(dp) :: occ(kbdim,klev),zgridbl(kbdim) 

! Executable statements

  !
  !*         1.    calculate additional cloud cover (occ)
  !                ------------------------
 
        zgridbl(1:kproma)=0._dp
        occ(1:kproma,:)=0._dp

  DO jk =ktdia,klev
     DO jl= 1,kproma
        IF (ktest(jl) == 1 .AND. pw_gwd_kpr(jl,jk) > 0._dp .AND. pqinucl(jl,jk) > ccwmin) THEN
           zgridbl(jl) = MAX(sqcst_2d(jl,krow)*111325._dp*1.8_dp,0.1_dp)
           occ(jl,jk) = MIN((2.0_dp*pl_z(jl)/zgridbl(jl)),1._dp)
           IF (psusatix_2d(jl,jk) > 0._dp) THEN
              occ(jl,jk) = 1._dp
           ENDIF
        ENDIF
        paclc(jl,jk)=paclc(jl,jk)+occ(jl,jk)
        paclc(jl,jk)=MIN(paclc(jl,jk),1._dp)
     ENDDO
  ENDDO

  END  SUBROUTINE orocirrus_cc
#endif
END MODULE mo_orocirrus
