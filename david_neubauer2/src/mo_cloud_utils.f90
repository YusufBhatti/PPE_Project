MODULE mo_cloud_utils
!-----------------------------------------------------------------------
! *mo_cloud_utils* contains utility routines and constants related to the 2-m cloud microphysics scheme
!
! Authors:
! --------
! Sylvaine Ferrachat, ETHZ   2008, 2009
!-----------------------------------------------------------------------
!

 USE mo_kind,               ONLY: dp
 USE mo_math_constants,     ONLY: pi
 USE mo_physical_constants, ONLY: grav, rhoh2o
 USE mo_echam_cloud_params, ONLY: cqtmin, crhoi, cn0s, crhosno !SF <--- BEWARE of crhosno!! It's different from 
                                  ! rhosno in mo_physical_constants.
 USE mo_param_switches,    ONLY: lconv     !csld #455

 IMPLICIT NONE

!----------------
! Public entities
!----------------

 REAL(dp), PARAMETER :: epsec = 1.0e-12_dp
 REAL(dp), PARAMETER :: xsec  = 1._dp - epsec
 REAL(dp), PARAMETER :: qsec  = 1._dp - cqtmin 
 REAL(dp), PARAMETER :: eps   = EPSILON(1.0_dp)
 REAL(dp), PARAMETER :: cri   = 10.E-6_dp  ! to estimate the number of produced  
                                           ! cloud droplets from ice melting in  
                                           ! case of licnc=.FALSE. [m]=> 10 um
 REAL(dp), PARAMETER :: mi = 4._dp/3._dp*cri**3*pi*crhoi ! assumed mass of ice crystals with 
                                                         ! corresponding volume mean radius cri
 REAL(dp), PARAMETER :: ri_vol_mean_1 = 2.166E-9_dp ! vol mean ice crystal radius, range border 1
 REAL(dp), PARAMETER :: ri_vol_mean_2 = 4.264E-8_dp ! vol mean ice crystal radius, range border 2
 REAL(dp), PARAMETER :: alfased_1 = 63292.4_dp ! for ice crystal fall velocity 
 REAL(dp), PARAMETER :: alfased_2 = 8.78_dp    ! for ice crystal fall velocity
 REAL(dp), PARAMETER :: alfased_3 = 329.75_dp  ! for ice crystal fall velocity
 REAL(dp), PARAMETER :: betased_1 = 0.5727_dp  ! for ice crystal fall velocity 
 REAL(dp), PARAMETER :: betased_2 = 0.0954_dp  ! for ice crystal fall velocity
 REAL(dp), PARAMETER :: betased_3 = 0.3091_dp  ! for ice crystal fall velocity
 !>>SF #475
 REAL(dp), PARAMETER :: cdnc_min_upper = 40.e6_dp       !upper value for min CDNC (used when ldyn_cdnc_min is TRUE)
 REAL(dp), PARAMETER :: cdnc_min_lower =  1.e6_dp       !lower value for min CDNC (used when ldyn_cdnc_min is TRUE)
 REAL(dp), PARAMETER :: rcd_vol_max = 19.e-6_dp         !maximum droplet volume radius [m] 
 !SF The above has been updated on 2016.04.04 (David Neubauer, pure atm run, HAM-M7, AR&G activation) 
 !<<SF #475
 REAL(dp), PARAMETER :: icemin  = 10._dp
 REAL(dp), PARAMETER :: icemax  = 1.e7_dp       ! maximum ice crystal number concentration
 REAL(dp), PARAMETER :: sigmaw  = 0.28_dp 
 REAL(dp), PARAMETER :: disp    = EXP(0.5_dp * sigmaw**2)
 REAL(dp), PARAMETER :: dw0     = 10.e-6_dp*disp       ! dispersion parameter in mu
 REAL(dp), PARAMETER :: cdi     = 0.6_dp
 REAL(dp), PARAMETER :: mw0     = 4.19e-12_dp
 REAL(dp), PARAMETER :: mi0     = 1.e-12_dp
 REAL(dp), PARAMETER :: mi0_rcp = 1.e12_dp
 REAL(dp), PARAMETER :: ka      = 0.024_dp
 REAL(dp), PARAMETER :: kb      = 1.38e-23_dp
 REAL(dp), PARAMETER :: alpha   = 0.5_dp        ! deposition coefficent alpha
 REAL(dp), PARAMETER :: xmw     = 2.992e-26_dp  ! mass of a h2o molecule [kg]
 REAL(dp), PARAMETER :: fall    = 3._dp
 REAL(dp), PARAMETER :: rhoice  = 925._dp
 REAL(dp), PARAMETER :: conv_effr2mvr = 0.9_dp  ! conversion ice crystal effective radius to mean volume radius
 REAL(dp), PARAMETER :: clc_min = 0.01_dp      ! minimum cloud cover below which cloud cover should be considered 0.
 REAL(dp), PARAMETER :: exm1_1    = 2.47_dp-1.0_dp
 REAL(dp), PARAMETER :: exp_1     = -1._dp / exm1_1
 REAL(dp), PARAMETER :: exm1_2    = 4.7_dp-1.0_dp
 REAL(dp), PARAMETER :: exp_2     = -1._dp / exm1_2
 REAL(dp), PARAMETER :: pirho     = pi*rhoh2o
 REAL(dp), PARAMETER :: pirho_rcp = 1._dp / pirho
 REAL(dp), PARAMETER :: cap       = 2._dp / pi
 REAL(dp), PARAMETER :: cons4     = 1._dp / ( pi*crhosno*cn0s )**0.8125_dp
 REAL(dp), PARAMETER :: cons5     = 1._dp / ( pi*crhosno*cn0s )**0.875_dp
 REAL(dp), PARAMETER :: fact_coll_eff = 0.09_dp ! Factor for the temperature-dep collection efficiency of snow by 
                                                ! cold hydrometeors
                                                ! Seifert & Beheng, Meteorol Atmos Phys 92, 45–66 (2006), eq. 67
                                                ! (#471)
 REAL(dp), PARAMETER :: fact_tke = 0.7_dp !SF #345

!>>SF #176  constants for ice crystal mass to effective radius relationship, from Pruppacher & Klett 1997 !     (tables 2.2a; 2.4a)
!     m = zfact_PK x (2r_i)**zpow_PK (m in g; r_i in cm)

!>>SF obsolete, former parametrization
 ! plate P1a, from table 2.4a (size range: 0.3-1.5 mm)
 !REAL(dp), PARAMETER :: fact_PK = 3.76e-2_dp
 !REAL(dp), PARAMETER :: pow_PK  = 3.31_dp
!<<SF obsolete

 ! plate P1a, from table 2.2a (size range: 10-3000 micrometers)
 REAL(dp), PARAMETER :: fact_PK = 8.253e-3_dp ! 9.17 x 0.9, because m=rho_c x V_c, where rho_c=0.9g.cm-3
 REAL(dp), PARAMETER :: pow_PK  = 2.475_dp
!<<SF

 CONTAINS

 SUBROUTINE get_util_var(kproma, kbdim, ktdia, klev, klevp1, &
                         paphm1, pgeo, papm1, ptm1,          &
                         pgeoh, pdp, pdpg, pdz, paaa, pviscos)

 !Get several utility variables:
 !    geopotential at half levels (pgeoh)
 !    pressure- and height-differences (pdp, pdz)
 !    air density correction for computing ice crystal fall velocity (paaa)
 !    dynamic viscosity of air (pviscos)
 !

   INTEGER, INTENT(IN) :: kproma, kbdim, ktdia, klev, klevp1

   REAL(dp), INTENT(IN)  :: paphm1(kbdim,klevp1), pgeo(kbdim,klev), papm1(kbdim,klev), ptm1(kbdim,klev)
   REAL(dp), INTENT(OUT) :: pgeoh(kbdim,klevp1),  pdp(kbdim,klev), pdpg(kbdim,klev), pdz(kbdim,klev), &
                            paaa(kbdim,klev), pviscos(kbdim,klev)

   REAL(dp) :: g_rcp

   g_rcp = 1._dp / grav

   !Geopotential at half levels:
   pgeoh(1:kproma,ktdia+1:klev) = 0.5_dp * (pgeo(1:kproma,ktdia+1:klev) + pgeo(1:kproma,ktdia:klev-1))
   pgeoh(1:kproma,ktdia)        = pgeo(1:kproma,ktdia) + (pgeo(1:kproma,ktdia)-pgeoh(1:kproma,ktdia+1))
   pgeoh(1:kproma,klevp1)       = 0.0_dp

   !Pressure differences:
   pdp(1:kproma,ktdia:klev)           = paphm1(1:kproma,ktdia+1:klevp1) - paphm1(1:kproma,ktdia:klev)
   pdpg(1:kproma,ktdia:klev)          = g_rcp * pdp(1:kproma,ktdia:klev)

   !Height differences:
   pdz(1:kproma,ktdia:klev)           = g_rcp * (pgeoh(1:kproma,ktdia:klev) - pgeoh(1:kproma,ktdia+1:klevp1))

   !Air density correction:
   paaa(1:kproma,:) = ( (papm1(1:kproma,:)/30000._dp)**(-0.178_dp) ) * ( (ptm1(1:kproma,:)/233.0_dp)**(-0.394_dp) )

   !Dynamic viscosity of air:
   pviscos(1:kproma,:) = (1.512_dp + 0.0052_dp*(ptm1(1:kproma,:)-233.15_dp))*1.e-5_dp

 END SUBROUTINE get_util_var

 SUBROUTINE get_cloud_bounds(kproma, kbdim, ktdia, klev, paclc, ktop, kbas, kcl_minustop, kcl_minusbas)

 ! *get_cloud_bounds* flags the top, base and intermediate levels for each cloud  
 !
   INTEGER, INTENT(IN)  :: kproma, kbdim, ktdia, klev
   INTEGER, INTENT(OUT) :: ktop(kbdim,klev), &      !flag for cloud tops 
                           kbas(kbdim,klev), &      !flag for cloud bases
                           kcl_minustop(kbdim,klev), & !flag for all cloud levels excepted their top
                          kcl_minusbas(kbdim,klev) !flag for all cloud levels excepted their base

   REAL(dp), INTENT(IN) :: paclc(kbdim,klev)    !cloud cover

   !Local variables:

   INTEGER  :: jk, jl, jm, jnumb, jtop, jbas
   INTEGER  :: iindex(kbdim,klev), &  !index of levels
               iclnb(kbdim),       &  !number of clouds per column
               iclbounds(2,klev/2+1)  !bounds infos per cloud

   REAL(dp) :: zaclcm(kbdim,klev), zaclcp(kbdim,klev)

   LOGICAL  :: ll(kbdim,klev), llm(kbdim,klev), llp(kbdim,klev), lltop(kbdim,klev), llbas(kbdim,klev)

   !Initialization:
   zaclcm(1:kproma,:)       = 0._dp
   zaclcp(1:kproma,:)       = 0._dp
   iclnb(1:kproma)          = 0
   kcl_minustop(1:kproma,:) = 0
   kcl_minusbas(1:kproma,:) = 0

   ll(1:kproma,:)    = .false.
   llm(1:kproma,:)   = .false.
   llp(1:kproma,:)   = .false.
   lltop(1:kproma,:) = .false.
   llbas(1:kproma,:) = .false.
 
   DO jk=ktdia,klev
      iindex(1:kproma,jk)=jk
   ENDDO 

   !Duplicates paclc at level-1 and level+1
   zaclcm(1:kproma,ktdia+1:klev) = paclc(1:kproma,ktdia:klev-1)
   zaclcp(1:kproma,ktdia:klev-1) = paclc(1:kproma,ktdia+1:klev)

   !Sets logical switches
   ll(1:kproma,:)  = (paclc(1:kproma,:)  >= epsec)   
   llm(1:kproma,:) = (zaclcm(1:kproma,:) <  epsec)  
   llp(1:kproma,:) = (zaclcp(1:kproma,:) <  epsec)

   lltop(1:kproma,:) = (ll(1:kproma,:) .AND. llm(1:kproma,:)) !true if cloud top detected
   llbas(1:kproma,:) = (ll(1:kproma,:) .AND. llp(1:kproma,:)) !true if cloud base detected

   !Sets itop and ibas
   ktop(1:kproma,:) = MERGE(iindex(1:kproma,:),0,lltop(1:kproma,:))
   kbas(1:kproma,:) = MERGE(iindex(1:kproma,:),0,llbas(1:kproma,:))

   !Resets the logical switches
   lltop(1:kproma,:) = .false.
   llbas(1:kproma,:) = .false.

   lltop(1:kproma,:) = (ktop(1:kproma,:) > 0)
   llbas(1:kproma,:) = (kbas(1:kproma,:) > 0)
 
   !Counts the number of clouds per column
   iclnb(1:kproma) = COUNT(lltop(1:kproma,:),2)

   DO jl=1,kproma
      jnumb = iclnb(jl)

      !Sets the bounds in a compact array
      iclbounds(1:2,:) = 0

      iclbounds(1,1:jnumb) = PACK(ktop(jl,:),lltop(jl,:))   !cloud tops
      iclbounds(2,1:jnumb) = PACK(kbas(jl,:),llbas(jl,:))   !cloud bases

      !Flags cloud levels excepted their base (or top -this later is not yet needed-)
      DO jm=1,jnumb
         jtop=iclbounds(1,jm)
         jbas=iclbounds(2,jm)
         kcl_minusbas(jl,jtop:jbas-1)=jbas
         kcl_minustop(jl,jtop+1:jbas)=jtop
      ENDDO
   ENDDO

 END SUBROUTINE get_cloud_bounds
!---------------------------------------------------------------------------------------------

 SUBROUTINE init_cloud_micro_2m

 ! Creates boundary conditions inherited from echam 
 ! (currently: quantities from convection)

   USE mo_boundary_condition,       ONLY: bc_nml, bc_define
   USE mo_external_field_processor, ONLY: EF_MODULE

   INTEGER :: ibc_cvcbot, ibc_wcape, ibc_tconv, ibc_detr_cond

   TYPE(bc_nml) :: bc_cvcbot, bc_wcape, bc_tconv, bc_detr_cond
    
   IF (lconv) THEN !csld #455 
      bc_cvcbot%ef_type = EF_MODULE
      ibc_cvcbot = bc_define('Convective cloud base index', &
                              bc_cvcbot, 2, .TRUE.)

      bc_wcape%ef_type = EF_MODULE
      ibc_wcape = bc_define('CAPE contrib. to conv. vertical velocity', &
                             bc_wcape, 2, .TRUE.)

      !>>SF #368
      bc_tconv%ef_type = EF_MODULE
      ibc_tconv = bc_define('Temperature in convective scheme', &
                            bc_tconv, 3, .TRUE.)
      !<<SF #368

      !>>SF #518
      bc_detr_cond%ef_type = EF_MODULE
      ibc_detr_cond = bc_define('Detrained condensate', &
                                bc_detr_cond, 3, .TRUE.)
      !<<SF #518
   ENDIF !csld #455 

 END SUBROUTINE init_cloud_micro_2m
  
END MODULE mo_cloud_utils
