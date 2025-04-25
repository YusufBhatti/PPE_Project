!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_chemistry.f90
!!
!! \brief
!! mo_ham_chemistry includes aerosol sulfur chemistry subroutines
!!
!! \author Martin Schultz (FZJ)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZJ) - general idea and code frame (2009-09)
!!   -#  K. Zhang  (MPI-Met) - rewrite the variable definition part
!!                             implementation in ECHAM6 and testing (2009-07)
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!!   - Feichter et al., Simulation of the tropospheric sulfur cycle in a global climate model,
!!     Atmospheric Enivronment 30, 1693-1706, 1996
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  MODULE mo_ham_chemistry
  
  PRIVATE
   
  PUBLIC :: ham_wet_chemistry   
  PUBLIC :: ham_gas_chemistry  
  
  CONTAINS
 
!--- Liquid-phase chemistry for sulfate (and other aerosol compounds)
!    This routine was originally part of xt_chemistry

!Note: the routine has been cleaned to avoid unnecessary variables etc., but some may remain
!++mgs: added pxtp1 argument - needed outside ham_wet_chemistry for MOZ

SUBROUTINE ham_wet_chemistry(kproma,  kbdim,  klev,      &
           ktop,       krow,       pmlwc,      pmiwc,    &
           paclc,      ptm1,       ptte,                 &
           pxtm1,      pxtte,      pxtp1,                &
           pxtp10,     pxtp1c,                           &
           paphp1,     papp1,      prhop1,     plfrac,   &
           pdpg     )


  !      J. Feichter     University of Hamburg           30/06/1992
  !      H.-S. Bauer     MPI (conversion to f90)         15/06/2000
  !      P. Stier        MPI-MET (changes, updates, 
  !                               modularisation   )      2001-2003
  !      P. Stier        MPI-MET (changes, updates, 
  !                               modularisation   )      2001-2003
  !      L.Pozzoli, J.S.Rast MPI-MET(coupled chemistry)  07/2005
  !      M. Schultz,     FZ Juelich, Dec 2008, extracted wet chemistry part from xt_chemistry
  !
  !      Purpose
  !      ---------
  !      Computes the homogeneous and heterogeneous oxidation 
  !      of the species for the sulfur cycle chemistry as described
  !      in Feichter et al. (1996). 
  !
  !      Alternatively, if HAMMOZ is activated, the routine uses the results from the MOZART
  !      sulfur chemistry.
  !
  USE mo_kind,               ONLY: dp
  USE mo_time_control,       ONLY: time_step_len, delta_time
  USE mo_physical_constants, ONLY: avo, argas
  USE mo_boundary_condition, ONLY: bc_apply
  USE mo_species,            ONLY: speclist
  USE mo_tracdef,            ONLY: ntrac
  USE mo_ham,                ONLY: ibc_o3, ibc_h2o2, mw_s, mw_so2,&
                                   mw_so4, nham_subm,             &
                                   HAM_M7, HAM_SALSA,HAM_BULK
  !TB added salsa specific variables
  USE mo_ham_salsa_trac,     ONLY: idt_so2_salsa => idt_so2, idt_ms4, idt_n
  USE mo_ham_m7_trac,        ONLY: idt_so2_m7 => idt_so2,         &
                                   idt_ms4as,  idt_ms4cs,         &
                                   idt_nas,    idt_ncs
  USE mo_ham_species,        ONLY: id_o3, id_h2o2, id_so2
  USE mo_ham_streams,        ONLY: d_prod_so4_liq,                &
                                   d_prod_ms4as,                  &
                                   d_prod_ms4cs
  USE mo_exception,          ONLY: finish
  USE mo_ham_salsactl,       ONLY: in1a, in2b, fn2b, in2a, fn2a !TB for SALSA
  
  IMPLICIT NONE

  !----------------------------------------------------
  ! Arguments
  INTEGER, INTENT(in)     :: kbdim, klev, kproma, ktop, krow
  REAL(dp), INTENT(in)    :: pmlwc(kbdim,klev),       &
                             pmiwc(kbdim,klev),       &
                             paclc(kbdim,klev),       &
                             ptm1(kbdim,klev),        &
                             ptte(kbdim,klev),        &
                             prhop1(kbdim,klev),      &
                             papp1(kbdim,klev),       &
                             paphp1(kbdim,klev+1),    &
                             pxtm1(kbdim,klev,ntrac), &
                             pxtp1(kbdim,klev,ntrac), &
                             pdpg(kbdim,klev)

  REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ntrac), &
                             plfrac(kbdim,klev),      &
                             pxtp1c(kbdim,klev,ntrac),&
                             pxtp10(kbdim,klev,ntrac)

  ! Local variables
  INTEGER , PARAMETER :: niter = 5  !gf cleanup (see #103)

  INTEGER :: jk, jl, jt, jn
  INTEGER :: idt_so2 !TD / SF: dummy index functional for both M7 and SALSA

!SF needed for selecting only relevant tracers
  INTEGER :: jnchem !total number of relevant tracers
  !TB chamged jchem to be allocatable
  INTEGER,ALLOCATABLE     :: jchem(:)
  INTEGER                 :: jtrac
!SF

  REAL(dp):: zrkh2o2(kbdim,klev),                               &
             zso4(kbdim,klev),         zdso4(kbdim,klev),       &
             ztp1(kbdim,klev)

  REAL(dp)::                           zzh2o2(kbdim,klev),      &  ! HAMMOZ uses h2o2 and o3
   &         zzo3(kbdim,klev)
  REAL(dp):: zxtte(kbdim,klev,ntrac)

  LOGICAL, PARAMETER :: lstrat=.TRUE.  ! Flag indicating call of
                                       ! xt_wetdep from the
                                       ! stratiform cloud scheme.

  REAL(dp):: zmin,zmolgh2o2, zmolgair, zmolgw,                   &
          ze1k, ze1h, ze3k, ze3h, zq298, zrgas,                  &
          zhpbase, znamair, zlwcmin, zlwcl, zlwcv,               &
          zhp, zqtp1, zrk, zrke, zh_so2, zpfac, zp_so2, zf_so2,  &
          zh_h2o2, zp_h2o2, zf_h2o2, x, zc, ze1, ze2,            &
          ze3, zfac1, zrkfac, zza, za21, za22, zph_o3, zf_o3,    &
          zdt, zh2o2m, zso2m, zso4m, zq, zso2mh, zdso2h,         &
          zso2l, zso4l, zzb, zzp, zzq, zzp2, zqhp,               &
          za2, zheneff, zrko3, zso2mo, zdso2o, zdso2tot,         &
          ztmst, zqtmst,                                         &
          zfraction

  REAL(dp) :: zhenry_so2(2)

  !--- Constants:

  REAL(dp), PARAMETER :: zso4massc=3.25E-15_dp ! Mass of coarse mode so4 particle
                                        ! as assumed for the in-cloud
                                        ! formation of sulfate aerosol in the
                                        ! absence of other particles. [kg]
                                        ! 3.25E-15 kg is equivalent to a
                                        ! mono-disperse radius of 0.75 um.

  !--- 0) Initialisations: -----------------------------------------

!csld commented and called in mo_submdel_interface.f90
!!$  !--- Mass conserving correction of negative tracer values:
!!$
!!$ CALL xt_borrow(kproma, kbdim,  klev,  ilevp1, ntrac, &
!!$                 papp1,  paphp1,                      &
!!$                 pxtm1,  pxtte                        )
!end csld

  !--- Constants:
  ztmst=time_step_len
  zqtmst=1._dp/ztmst
  zdt=ztmst/REAL(niter,dp) !gf cleanup (see #103)

  !--- Initialisation

  plfrac(:,:)     = 0._dp
  zxtte(:,:,:)    = 0._dp             ! Preset initial tracer tendency
  zmin            = 1.E-20_dp         
  zdso4(:,:)      = 0._dp             ! Sulfate produced by heterogeneous
  zhenry_so2(:)   = speclist(id_so2)%henry(:) ! Henry's law constant and activation energy

  !--- Molecular weights in g:

!!mgs(S)!!  zmolgs=32.064_dp
!!mgs(S)!!  zmolgso2=64.06_dp          ! HAMMOZ
  zmolgh2o2=34.01474_dp
  zmolgair=28.84_dp
  zmolgw=18.015_dp
  zhpbase=2.5e-06_dp
  ze1k=1.1e-02_dp
  ze1h=2300._dp
  ze3k=1.2e-02_dp
  ze3h=2010._dp
  zq298=1._dp/298._dp
  zrgas=8.2e-02_dp
  znamair=1.e-03_dp*avo/zmolgair
  zlwcmin=1.e-07_dp                    ! Liquid water content threshold

  !SF-- Creates an array of indices for the relevant tracers
  !
  !     this allows to avoid looping over all tracers, and also avoids side-effects on other ones
  !     with operations that should cancel each other, but don't really, numerically.

  !TB moved jnchem here so it can depend on the aerosolmodel
  SELECT CASE (nham_subm)
      !CASE(HAM_BULK)

      CASE(HAM_M7)
         !TB
         !fixed jnchem for M7 for sulfate produced in clouds
         jnchem=5
         ALLOCATE(jchem(jnchem))
         idt_so2  = idt_so2_m7

         jchem(1) = idt_so2
         jchem(2) = idt_ms4cs
         jchem(3) = idt_ms4as
         jchem(4) = idt_ncs
         jchem(5) = idt_nas

  CASE(HAM_SALSA)
         !TB 
         !set jnchem dynamically depending on number of bins
         jnchem=1+fn2a*2
         ALLOCATE(jchem(jnchem))
         idt_so2  = idt_so2_salsa

         jchem(1) = idt_so2
         jchem((1+in1a):(1+fn2a)) = idt_ms4(in1a:fn2a)
         jchem(2+fn2a:1+fn2a*2) = idt_n(in1a:fn2a)

  END SELECT

  !
  !--- Calculate temperature at t=t+1:

  ztp1(1:kproma,:) = ptm1(1:kproma,:)+ptte(1:kproma,:)*time_step_len

  !--- Get oxidant concentrations
  !--- and convert from mass mixing ratio to concentration [molecule cm-3]
!Note: this takes care of offline oxidants as well as online (HAMMOZ) fields!

  IF (ibc_o3 .gt. 0) THEN 
     CALL bc_apply(ibc_o3,   kproma, krow, zzo3)
  ELSE
     CALL finish('ham_wet_chemistry','Boundary condition for O3 not found') !gf #264
  ENDIF

  IF (ibc_h2o2 .gt. 0) THEN 
     CALL bc_apply(ibc_h2o2, kproma, krow, zzh2o2)
  ELSE
     CALL finish('ham_wet_chemistry','Boundary condition for H2O2 not found') !gf #264
  ENDIF

  DO jk=1,klev
    DO jl=1,kproma
      zc=(papp1(jl,jk)/(argas*ptm1(jl,jk)))*1.E-6_dp*avo*zmolgair
      !  (          mol m-3              )
      !  (          mol cm-3                   )
      !  (          molecules cm-3                 )
      zzh2o2(jl,jk) = zzh2o2(jl,jk) * zc / speclist(id_h2o2)%moleweight
      zzo3(jl,jk)   = zzo3(jl,jk)   * zc / speclist(id_o3)%moleweight
    END DO
  END DO
  !
  !

  !--- Calculate total sulfate in the accumulation and coarse mode:
  !TB / SF: adapted for handling SALSA
  SELECT CASE(nham_subm)
      !CASE(HAM_BULK)
         
      CASE(HAM_M7)
         zso4(1:kproma,:)=pxtp1(1:kproma,:,idt_ms4as)+ pxtp1(1:kproma,:,idt_ms4cs)
         zso4(1:kproma,:)=MAX(zso4(1:kproma,:),0._dp)
         
      CASE(HAM_SALSA)
         zso4(1:kproma,:)=sum(pxtp1(1:kproma,:,idt_ms4(in2a:fn2a)),3)
         zso4(1:kproma,:)=MAX(zso4(1:kproma,:),0._dp)

  END SELECT

  !
  !--- 1) Calculate the reaction-rates for SO2-H2O2:
  !
  DO jk=ktop,klev
     DO jl=1,kproma
!gf see #103        IF(pmlwc(jl,jk).GT.zmin) THEN
        IF(pmlwc(jl,jk).GT.zlwcmin) THEN
           zlwcl=pmlwc(jl,jk)*prhop1(jl,jk)*1.e-06_dp
           zlwcv=pmlwc(jl,jk)*prhop1(jl,jk)*1.e-03_dp
!!mgs(S)!!           zhp=zhpbase+zso4(jl,jk)*1000._dp/(pmlwc(jl,jk)*zmolgs)
           zhp=zhpbase+zso4(jl,jk)*1000._dp/(pmlwc(jl,jk)*mw_so4)
           zqtp1=1._dp/ptm1(jl,jk)-zq298
           zrk=8.e+04_dp * exp(-3650._dp*zqtp1)/(0.1_dp+zhp)
           zrke=zrk/(zlwcl*avo)
           !
           zh_so2=zhenry_so2(1) * exp(zhenry_so2(2)*zqtp1)
           zpfac=zrgas*zlwcv*ptm1(jl,jk)
           zp_so2=zh_so2*zpfac
           zf_so2=zp_so2/(1._dp+zp_so2)
           !
           zh_h2o2=9.7e+04_dp * exp(6600._dp*zqtp1)
           zp_h2o2=zh_h2o2*zpfac
           zf_h2o2=zp_h2o2/(1._dp+zp_h2o2)
           !
           zrkh2o2(jl,jk)=zrke*zf_so2*zf_h2o2
        ELSE
           zrkh2o2(jl,jk)=0._dp
        END IF
     END DO
  END DO
  !
  !--- 2) Heterogenous chemistry:
  !

  !OCL NOVREC,NOALIAS

  DO jk=ktop,klev
     DO jl=1,kproma

!gf see #103       IF(pxtp1(jl,jk,idt_so2).GT.zmin.AND.pmlwc(jl,jk).GT.zmin) THEN
        IF(pxtp1(jl,jk,idt_so2).GT.zmin.AND.pmlwc(jl,jk).GT.zlwcmin) THEN
           x=prhop1(jl,jk)
           !
           zqtp1=1._dp/ptm1(jl,jk)-zq298
           ze1=ze1k * exp(ze1h*zqtp1)
           ze2=zhenry_so2(1) * exp(zhenry_so2(2)*zqtp1)
           ze3=ze3k * exp(ze3h*zqtp1)

           !--- Calculate the liquid water content in [l/cm**3]:

           zlwcl=pmlwc(jl,jk)*prhop1(jl,jk)*1.e-06_dp

           !--- Calculate conversion factor for [molec/cm**3(air)] to [mole/l(aq)]:

           zfac1=1._dp/(zlwcl*avo)

           !--- Calculate the liquid water volume fraction [m**3/m**3]:

           zlwcv=pmlwc(jl,jk)*prhop1(jl,jk)*1.e-03_dp

           !--- Calculate conversion factor R_gas * T * V(aq)/V(air):

           zrkfac=zrgas*ptm1(jl,jk)*zlwcv

           !--- Multiplication with the Henry coefficient for SO2:

           zza=ze2*zrkfac

           za21=4.39e+11_dp*EXP(-4131._dp/ptm1(jl,jk))    ! k51 in Feichter et al.
           za22=2.56e+03_dp*EXP(-926._dp/ptm1(jl,jk))     ! k52 in Feichter et al.
           zph_o3=ze1*zrkfac
           zf_o3=zph_o3/(1._dp+zph_o3)
           !
           zh2o2m=zzh2o2(jl,jk)
           zso2m=pxtp1(jl,jk,idt_so2)*xtoc(x,mw_so2)
           zso4m=zso4(jl,jk)*xtoc(x,mw_so4)

!CDIR UNROLL=5
           DO jn=1,niter
              zq=zrkh2o2(jl,jk)*zh2o2m
              zso2mh=zso2m*EXP(-zq*zdt)
              !
              zdso2h=zso2m-zso2mh
              zh2o2m=zh2o2m-zdso2h
              zh2o2m=MAX(0._dp,zh2o2m)
              !
              zso4m=zso4m+zdso2h
              !   calculate the ph value
              zso2l=zso2mh*zfac1
              zso4l=zso4m*zfac1
              zzb=zhpbase+zso4l
              zzp=(zza*ze3-zzb-zza*zzb)/(1._dp+zza)
              zzq=-zza*ze3*(zzb+zso2l)/(1._dp+zza)
              zzp=0.5_dp*zzp
              zzp2=zzp*zzp
              zhp=-zzp+SQRT(zzp2-zzq)
              zqhp=1._dp/zhp
              !
              !--- Calculate the reaction rate for SO2-O3:
              !
              za2=(za21+za22*zqhp)*zfac1

              !--- Calculate effective Henry coefficient neglecting SO3-- :

              zheneff=1._dp+ze3*zqhp
              zp_so2=zza*zheneff

              !--- Calculate the liquid fraction of total SO2, i.e. SO2(aq)/SO2(tot):

              zf_so2=zp_so2/(1._dp+zp_so2)

              zrko3=za2*zf_o3*zf_so2
              !
              zq=zzo3(jl,jk)*zrko3
              zso2mo=zso2mh*EXP(-zq*zdt)
              zdso2o=zso2mh-zso2mo
              zso4m=zso4m+zdso2o
              zso2m=zso2mo
           END DO

           zdso2tot=pxtp1(jl,jk,idt_so2)-zso2m*ctox(x,mw_so2)
           zdso2tot=MIN(zdso2tot,pxtp1(jl,jk,idt_so2))
           pxtp1c(jl,jk,idt_so2)=pxtp1(jl,jk,idt_so2)-zdso2tot
           !--- zdso2tot is the total amount of SO2 oxidised by
           !    heterogeneous chemistry to sulfate, i.e. d(SO4):
!gf           zdso4(jl,jk)=zdso2tot
           zdso4(jl,jk)=zdso2tot*(mw_so4/mw_so2)

           !--- Liquid fraction of the total SO2:

           plfrac(jl,jk)=zf_so2

        END IF
     END DO
  END DO

  !---   Sulfate produced within cloud droplets is distributed
  !      to accumulation and coarse mode sulfuric acid mass:

  DO jk=ktop,klev
     DO jl=1,kproma
        SELECT CASE(nham_subm)
            CASE(HAM_M7)

               !--- Accumulation and Coarse Mode particles present:
               
               IF ((pxtp1(jl,jk,idt_nas) >= 1.E-3_dp) .AND. (pxtp1(jl,jk,idt_ncs) >= 1.E-3_dp)) THEN
                  
                  zfraction=pxtp1(jl,jk,idt_nas)/(pxtp1(jl,jk,idt_nas)+pxtp1(jl,jk,idt_ncs))
                  
                  pxtp1c(jl,jk,idt_ms4as)=pxtp1(jl,jk,idt_ms4as) + zdso4(jl,jk)*zfraction
                  
                  pxtp1c(jl,jk,idt_ms4cs)=pxtp1(jl,jk,idt_ms4cs) + zdso4(jl,jk)*(1-zfraction)
                  
                  !--- Accumulation but no Coarse Mode particles present:
                  
               ELSE IF ((pxtp1(jl,jk,idt_nas) >= 1.E-3_dp) .AND. (pxtp1(jl,jk,idt_ncs) < 1.E-3_dp)) THEN
                  
                  pxtp1c(jl,jk,idt_ms4as)=pxtp1(jl,jk,idt_ms4as) + zdso4(jl,jk)
                  
                  !--- No Accumulation but Coarse Mode particles present:
                  
               ELSE IF ((pxtp1(jl,jk,idt_nas) < 1.E-3_dp) .AND. (pxtp1(jl,jk,idt_ncs) >= 1.E-3_dp)) THEN
                  
                  pxtp1c(jl,jk,idt_ms4cs)=pxtp1(jl,jk,idt_ms4cs) + zdso4(jl,jk)
                  
                  !--- Neither Accumulation Mode nor Coarse Mode particles present.
                  !
                  !    Assumption: Cloud water sulfate is assumed to exits as
                  !                coarse mode particles as this resembles most
                  !                the observed "droplet mode"
                  !                (See Seinfeld & Pandis, p. 821, 1998)
                  
               ELSE IF ((pxtp1(jl,jk,idt_nas) < 1.E-3_dp) .AND. (pxtp1(jl,jk,idt_ncs) < 1.E-3_dp)) THEN
                  
                  pxtp1c(jl,jk,idt_ms4cs)=pxtp1(jl,jk,idt_ms4cs) + zdso4(jl,jk)
                  
                  !--- Increase the coarse mode soluble particle number tendency:
                  
    !gf see #103           pxtp1c(jl,jk,idt_ncs)=pxtp1c(jl,jk,idt_ncs)+(zdso4(jl,jk)/zso4massc)
                  pxtp1c(jl,jk,idt_ncs)=pxtp1c(jl,jk,idt_ncs)+(zdso4(jl,jk)/(zso4massc*mw_so4/mw_s))
                  
               END IF

            CASE (HAM_SALSA)
    
               IF ( sum(pxtp1(jl,jk,idt_n(in2a:fn2a))) >= 1.E-3_dp .AND. sum(pxtp1(jl,jk,idt_ms4(in2a:fn2a))) >= 1.E-30_dp ) THEN
                  
                  pxtp1c( jl,jk,idt_ms4(in2a:fn2a) ) = pxtp1( jl,jk, idt_ms4(in2a:fn2a) ) + &
                       zdso4(jl,jk) *pxtp1( jl,jk,idt_ms4(in2a:fn2a) ) / &
                       sum( pxtp1(jl,jk,idt_ms4(in2a:fn2a)) )
               END IF

        END SELECT

     END DO
  END DO

tracer : DO jt=1, jnchem !SF relevant tracers only

     jtrac = jchem(jt) 

     DO jk=ktop,klev
        DO jl=1,kproma
           zxtte(jl,jk,jtrac)=((pxtp10(jl,jk,jtrac)*(1._dp-paclc(jl,jk))+ &
                             pxtp1c(jl,jk,jtrac)*paclc(jl,jk))-pxtm1(jl,jk,jtrac))*zqtmst - &
                           pxtte(jl,jk,jtrac)
           !
           !--- Security check:
           !
           !    Total tendency is not allowed to exeed 
           !    threshold that produces negative values.
           !
           !    pxtp1 = pxtm1 + (pxtte+zxtte)*dt
           !
           !<=> 0     < pxtm1 + (pxtte+zxtte)*dt
           !
           !               pxtm1  
           ! => zxtte > - ------- - pxtte
           !                dt

           !SF note: be aware that it has side-effects on tracers that
           !         ie idt_so2, idt_ms4cs, idt_ms4as, idt_ncs
           !         this is not 

           zxtte(jl,jk,jtrac)=MAX(-pxtm1(jl,jk,jtrac)*zqtmst-pxtte(jl,jk,jtrac) , zxtte(jl,jk,jtrac))

           !--- Change of total tendencies:

           pxtte(jl,jk,jtrac)=pxtte(jl,jk,jtrac)+zxtte(jl,jk,jtrac)

        END DO
     END DO
  END DO tracer

  !--- Store sulfate production for diagnostics:
  SELECT CASE(nham_subm)
      CASE(HAM_M7)
         
         DO jk=1, klev
            
            d_prod_ms4as(1:kproma,krow)=d_prod_ms4as(1:kproma,krow) + &
                                        zxtte(1:kproma,jk,idt_ms4as)*pdpg(1:kproma,jk)*delta_time
            
            d_prod_ms4cs(1:kproma,krow)=d_prod_ms4cs(1:kproma,krow) + &
                                        zxtte(1:kproma,jk,idt_ms4cs)*pdpg(1:kproma,jk)*delta_time
            
            d_prod_so4_liq(1:kproma,jk,krow)=d_prod_so4_liq(1:kproma,jk,krow) + &
                                             (zxtte(1:kproma,jk,idt_ms4as)+zxtte(1:kproma,jk,idt_ms4cs)) * &
                                             pdpg(1:kproma,jk)*delta_time
    
         END DO

      CASE(HAM_SALSA)

         DO jk=1, klev
            
    !        d_prod_ms4(1:kproma,krow)=d_prod_ms4(1:kproma,krow) + &
    !             sum(zxtte(1:kproma,jk,idt_ms4(:))*pdpg(1:kproma,jk)*delta_time,3)
            
            d_prod_so4_liq(1:kproma,jk,krow)=d_prod_so4_liq(1:kproma,jk,krow) + &
                 sum(zxtte(1:kproma,jk,idt_ms4(:)),2) * pdpg(1:kproma,jk)*delta_time
            
         END DO

  END SELECT
END SUBROUTINE ham_wet_chemistry

!--- Gas-phase chemistry for SO2 and DMS oxidation
!    This routine was originally part of xt_chemistry

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ham_gas_chemistry(kbdim, klev, kproma, krow,        &
                             paph, pap, pt, pq, pxtm1, pxtte   )

  ! This routine is called from *chemi_physc_2* in *mo_submodel_interface*
  ! if lham is set and lhmzoxi is not set.
  !
  ! Method: compute gas-phase chemistry for HAM aerosol module
  !
  ! Authors:
  !      J. Feichter,    University of Hamburg           30/06/1992, original code 
  !      H.-S. Bauer     MPI (conversion to f90)         15/06/2000
  !      P. Stier        MPI-MET (changes, updates,
  !                               modularisation   )      2001-2003
  !      P. Stier        MPI-MET (changes, updates,
  !                               modularisation   )      2001-2003
  !      L.Pozzoli, J.S.Rast MPI-MET(coupled chemistry)  07/2005
  !      M. Schultz,     FZ-Juelich, Dec 2008,    extracted from xt_chemistry
  !      D. O'Donnell, FMI, July 2011: bugfixes, cleanup, some optimisation

  USE mo_kind,                ONLY: dp
  USE mo_math_constants,      ONLY: pi
  USE mo_physical_constants,  ONLY: grav, avo, argas, rd, vtmpc1
  USE mo_time_control,        ONLY: time_step_len, delta_time, current_date
  USE mo_time_conversion,     ONLY: TC_get
  USE mo_geoloc,              ONLY: rdayl_x, philon_2d
  USE mo_boundary_condition,  ONLY: bc_apply
  USE mo_species,             ONLY: speclist
  USE mo_ham_species,         ONLY: id_oh, id_o3, id_no2, id_so2, id_dms, id_so4g, &
                                    id_no3 !gf #146
  USE mo_tracdef,             ONLY: ntrac
  USE mo_ham,                 ONLY: ibc_oh, ibc_o3, ibc_no2, mw_dms, &
                                    ibc_no3 !gf #146
  USE mo_ham_streams,         ONLY: daylength,          &
                                    d_prod_so2_dms_oh,  &
                                    d_prod_so4_dms_oh,  &
                                    d_prod_so2_dms_no3, &
                                    d_prod_so4_so2_oh,  & 
                                    d_prod_h2so4
  USE mo_exception,          ONLY: finish

  IMPLICIT NONE

  !--- Arguments
  INTEGER, INTENT(in)     :: kbdim, klev, kproma, krow
  REAL(dp), INTENT(in)    :: pap(kbdim,klev),      &          ! full level pressure
                             paph(kbdim,klev+1),   &          ! half level pressure
!!mgs!!                             prhop1(kbdim,klev),   &
                             pt(kbdim,klev),       &          ! temperature
                             pq(kbdim,klev),       &          ! specific humidity
                             pxtm1(kbdim,klev,ntrac)          ! tracer concentrations
  REAL(dp), INTENT(inout) :: pxtte(kbdim,klev,ntrac)          ! tracer concentration tendencies
                             

!>>gf see #103
  !---local data
  !   parameters

  !   Reaction rate data
  !   SO2 + OH + M
  REAL(dp), PARAMETER :: zk2i=2.0e-12_dp        ! High pressure limit
  REAL(dp), PARAMETER :: zk2=4.0e-31_dp         ! Low pressure limit
  REAL(dp), PARAMETER :: zk2f=0.45_dp           ! Broadening factor

  !   DMS-NO3:
  REAL(dp), PARAMETER :: zk3=1.9e-13_dp

  REAL(dp), PARAMETER :: zkn2o5aq=0.1e-04_dp

  !   Miscellaneous:
  REAL(dp), PARAMETER :: molgair=28.84_dp      ! Molecular weight of dry air  
  REAL(dp), PARAMETER :: znamair=1.e-03_dp*avo/molgair   ! Molecules air cm-3
!<<gf

  !--- Local variables
  INTEGER    :: jk, jl, &
                iday, isecond

  INTEGER    :: idt_dms, idt_so2, idt_so4 ! Tracer indices !gf see #103

  REAL(dp)   :: zrdayl(kbdim), &
                zdayfac(kbdim)            ! Weighting factor for diurnal OH cycle when monthly
                                          ! mean OH fields are used. (zdayfac=1/(relative daylength))
  REAL(dp)   :: zzoh(kbdim,klev),       &
                zzo3(kbdim,klev),       &
                zzno2(kbdim,klev),      &
!gf see #146
                zzno3(kbdim,klev)
!gf
  !>>dod changed these to scalars 
  REAL(dp)   :: zdso2,                  & ! Loss of SO2, SO2+OH
                zso4so2oh,              & ! Production of SO4, SO2+OH
                zdms2so4,               & ! Production of SO4, DMS+OH
                zdms2so2oh,             & ! Production of SO2, DMS+OH
                zdms2so2no3               ! Production of SO2, DMS+NO3
  !<<dod
  REAL(dp)   :: zdpg(kbdim,klev)

  REAL(dp)   :: slt    (kbdim)

  !>>dod bugfix negative dayfac values
  LOGICAL :: ldaylight(kbdim)
  !<<dod

  !>>dod bugfix(#49)
  REAL(dp) :: zconv_so2_so4, zconv_dms_so2, zconv_dms_so4  ! Factors to convert mass(SO2) -> mass(SO4), etc.
  !<<dod

  REAL(dp)   :: zt, ztk1, ztk2, ztk12, ztk3, ztk23b,     & ! temperature and rates
                zkno2o3, zkno2no3, zkn2o5,               &
                zeqn2o5,                                 & ! N2O5 equilibrium
                zc, zhil,                                &  
                ztmst, zqtmst, zqt, zqt3, zzq,           & ! inverse time step etc.
                zo2, zdms, zno3,                         & ! species mixing ratios
                ztoso2, zexp, zrx1, zrx2, zeps             ! auxillary variabes 

  !>>dod changed to scalar and back
  REAL(dp) :: zrho(kbdim,klev)              !! ++mgs
  REAL(dp) :: zmair(kbdim,klev)
  !<<dod
  REAL(dp) :: zdms0(kbdim,klev), zso20(kbdim,klev)
  REAL(dp) :: ztmp1(kbdim) !SF #458


  !--- Initialisation
  zeps=EPSILON(1._dp)
  ztmst=time_step_len
  zqtmst=1._dp/ztmst

!>>gf see #103
  zno3 = 0._dp
  zdms = 0._dp

  !---tracer indices
  idt_dms = speclist(id_dms)%idt
  idt_so2 = speclist(id_so2)%idt
  idt_so4 = speclist(id_so4g)%idt

  !---ratios of molecular weights: S compounds, for mass conversions in chemistry
  zconv_so2_so4 = speclist(id_so4g)%moleweight / speclist(id_so2)%moleweight
  zconv_dms_so2 = speclist(id_so2)%moleweight / speclist(id_dms)%moleweight
  zconv_dms_so4 = speclist(id_so4g)%moleweight / speclist(id_dms)%moleweight
!<<gf

  !--- Current UT:
  CALL TC_get(current_date,iday,isecond)

  !>> dod loop to calculate air density and air mass/unit area deleted, merged into main loop instead <<dod

  zrdayl(1:kproma) = rdayl_x(1:kproma,krow)

  ! Solar local time (as fraction of the day):
  slt(1:kproma) = mod(1.15741E-5_dp*(isecond + philon_2d(1:kproma,krow)*240._dp),1._dp)

  !>>SF #458 (replacing WHERE statements)

  !dod correction of negative zdayfac values:
  !      check if solar local time is within half a day of local noon 
  ldaylight(1:kproma) = (daylength(1:kproma,krow) > zeps) .AND. &
                        (slt(1:kproma) > (0.5_dp-0.5_dp*daylength(1:kproma,krow)) ) .AND. &
                        (slt(1:kproma) < (0.5_dp+0.5_dp*daylength(1:kproma,krow)) )

  ztmp1(1:kproma) = MERGE(daylength(1:kproma,krow), 1._dp, ldaylight(1:kproma)) !SF 1. is a dummy val.

  ! First we account for daylength, so that the integral of the OH
  ! concentration over a month gives the monthly mean (this introduces an
  ! artificial variation of the OH concentration over the month):
  ztmp1(1:kproma) = 1._dp / ztmp1(1:kproma)
  
  ! Then we account for the diurnal cycle so that the OH concentration
  ! follows a cosine peak between sunrise and sunset:
  ztmp1(1:kproma) = ztmp1(1:kproma)*0.5_dp*pi*cos(pi*(slt(1:kproma)-0.5_dp)*ztmp1(1:kproma))
  
  zdayfac(1:kproma) = MERGE(ztmp1(1:kproma), 0._dp, ldaylight(1:kproma))
  !<<SF #458 (replacing WHERE statements)

  !--- Get oxidant concentrations and convert to [molecules cm-3]

  IF (ibc_oh .gt. 0) THEN 
     CALL bc_apply(ibc_oh,   kproma, krow, zzoh)
  ELSE
     CALL finish('ham_wet_chemistry','Boundary condition for OH not found') !gf #264
  ENDIF

  IF (ibc_o3 .gt. 0) THEN 
     CALL bc_apply(ibc_o3,   kproma, krow, zzo3)
  ELSE
     CALL finish('ham_wet_chemistry','Boundary condition for O3 not found') !gf #264
  ENDIF

  IF (ibc_no2 .gt. 0) THEN 
     CALL bc_apply(ibc_no2,   kproma, krow, zzno2)
  ELSE
     CALL finish('ham_wet_chemistry','Boundary condition for NO2 not found') !gf #264
  ENDIF

  IF (ibc_no3 .gt. 0) THEN 
     CALL bc_apply(ibc_no3,   kproma, krow, zzno3)
  ELSE
     CALL finish('ham_wet_chemistry','Boundary condition for NO3 not found') !gf #264
  ENDIF

  DO jk=1,klev
    DO jl=1,kproma
      zc=(pap(jl,jk)/(argas*pt(jl,jk)))*1.E-6_dp*avo*molgair
      !  (          mol m-3              )
      !  (          mol cm-3                   )
      !  (          molecules cm-3                 )
      zzoh(jl,jk)  = zzoh(jl,jk) * zc / speclist(id_oh)%moleweight
      zzo3(jl,jk)  = zzo3(jl,jk) * zc / speclist(id_o3)%moleweight
      zzno2(jl,jk) = zzno2(jl,jk)* zc / speclist(id_no2)%moleweight

!gf see #146
      zzno3(jl,jk) = zzno3(jl,jk)* zc / speclist(id_no3)%moleweight
!gf

    END DO
  END DO

!>>gf see #103
  !---vectorizable calculations
  DO jk=1,klev
     DO jl=1,kproma
        !---air density [kg m-3]
        zrho(jl,jk) = pap(jl,jk)/(pt(jl,jk)*rd*(1._dp+vtmpc1*pq(jl,jk)))

        !---air mass per unit area (converts mmr -> kg m-2)
        zdpg(jl,jk) = (paph(jl,jk+1) - paph(jl,jk))/grav

        !---tracer concentrations at t=t+dt before chemistry for DMS and SO2
        zdms0(jl,jk) = pxtm1(jl,jk,idt_dms) + pxtte(jl,jk,idt_dms)*ztmst
        zso20(jl,jk) = pxtm1(jl,jk,idt_so2) + pxtte(jl,jk,idt_so2)*ztmst

        ! molecular density of air
        zmair(jl,jk)=zrho(jl,jk)*znamair  
     END DO
  END DO
!<<gf

  !++mgs: kept original numbering for better recognition
  !--- 4) Gas phase chemistry:
  !
  !OCL NOVREC,NOALIAS
   
  DO jk=1,klev
     DO jl=1,kproma
        IF(NINT(zrdayl(jl)).EQ.1) THEN  !--- 4.1) Day-time chemistry (oxidation with OH)

           ztk2=zk2*(pt(jl,jk)/300._dp)**(-3.3_dp)                ! k0 in Feichter et al.
           zhil=ztk2*zmair(jl,jk)/zk2i                            ! see Appendix A in Feichter et al.
           zexp=LOG10(zhil)
           zexp=1._dp/(1._dp+zexp*zexp)
           ztk23b=ztk2*zmair(jl,jk)/(1._dp+zhil)*zk2f**zexp       ! reaction rate SO2+OH+M 

           zdso2=zso20(jl,jk)*zzoh(jl,jk)*ztk23b*zdayfac(jl)      ! change in SO2
           zdso2=MIN(zdso2,zso20(jl,jk)*zqtmst)                   ! limit change to available SO2
           zdso2=MAX(zdso2,0._dp)                                 ! probably unnecessary

           pxtte(jl,jk,idt_so2)=pxtte(jl,jk,idt_so2)-zdso2        ! update tendencies for SO2
           !>>dod bugfix(#49)
           zso4so2oh = zdso2*zconv_so2_so4
           pxtte(jl,jk,idt_so4)=pxtte(jl,jk,idt_so4)+zso4so2oh    ! update tendencies for SO4
           !<<dod

           zt=pt(jl,jk)
           !   zo2 = 21% of air density [molec cm-3]
           zo2=0.21_dp*zmair(jl,jk)
           !   H abstraction
           ztk1=9.6e-12_dp*EXP(-234._dp/zt)
           !   OH addition
           ztk2=1.7e-42_dp*EXP(7810._dp/zt)*zo2/(1._dp+5.5e-31_dp*EXP(7460._dp/zt)*zo2)
           ztk12=ztk1+ztk2

           zdms=zdms0(jl,jk)*zzoh(jl,jk)*zdayfac(jl)*ztk12       !>>dod deleted double conversion of DMS <<dod
           zdms=MIN(zdms,zdms0(jl,jk)*zqtmst)

           pxtte(jl,jk,idt_dms)=pxtte(jl,jk,idt_dms)-zdms
              !
              !--- ztoso2 is the fraction of DMS oxidized to SO2:
              !
           ztoso2=(ztk1+0.75_dp*ztk2)/(ztk1+ztk2)

           zdms2so2oh = zdms*ztoso2*zconv_dms_so2       !>>dod bugfix(#49) <<dod

           pxtte(jl,jk,idt_so2)=pxtte(jl,jk,idt_so2)+zdms2so2oh

           !--- (1-ztoso2) is converted to MSA and assumed to occur as sulfate:

           zdms2so4 = zdms*(1._dp-ztoso2)*zconv_dms_so4

           pxtte(jl,jk,idt_so4)=pxtte(jl,jk,idt_so4)+zdms2so4

           !--- Store production of sulfate and SO2 for diagnostics:
           d_prod_so4_dms_oh(jl,krow)=d_prod_so4_dms_oh(jl,krow) + &
                                      zdms2so4*zdpg(jl,jk)*delta_time

           d_prod_so2_dms_oh(jl,krow)=d_prod_so2_dms_oh(jl,krow) + &
                                      zdms2so2oh*zdpg(jl,jk)*delta_time

           d_prod_so4_so2_oh(jl,krow)=d_prod_so4_so2_oh(jl,krow) + &
                                      zso4so2oh*zdpg(jl,jk)*delta_time

           d_prod_h2so4(jl,jk,krow)=d_prod_h2so4(jl,jk,krow) + &
                                    (zso4so2oh+zdms2so4)*zdpg(jl,jk)*delta_time

        ELSE   !--- 4.2) Night-time chemistry 

              ztk3=zk3*EXP(520._dp/pt(jl,jk))

!gf see #146 - Values from input file
              zno3=zzno3(jl,jk)
!gf
              zdms=zdms0(jl,jk)*zno3*ztk3
              zdms=MIN(zdms,zdms0(jl,jk)*zqtmst)
              pxtte(jl,jk,idt_dms)=pxtte(jl,jk,idt_dms)-zdms
              zdms2so2no3=zdms*zconv_dms_so2                                    !>>dod bugfix(#49) <<dod
              pxtte(jl,jk,idt_so2)=pxtte(jl,jk,idt_so2)+zdms2so2no3  

   !--- Store producion of sulfate and SO2 for diagnostics:

              d_prod_so2_dms_no3(jl,krow)=d_prod_so2_dms_no3(jl,krow) + &
                                          zdms2so2no3*zdpg(jl,jk)*delta_time

        ENDIF

     END DO
  END DO


  END SUBROUTINE ham_gas_chemistry
  

  FUNCTION xtoc(x,y)

    !
    !     Function for changing the units
    !     from mass-mixing ratio to molecules per cm**3 and vice versa
    !
    !   x = density of air, y = mol weight in gramm

    USE mo_kind,     ONLY: dp
    IMPLICIT NONE

    ! Function return value
    REAL(dp) :: xtoc

    ! Scalar arguments
    REAL(dp) :: x, y

    xtoc=x*6.022e+20_dp/y

  END FUNCTION xtoc

  FUNCTION ctox(x,y)

    !
    !     Function for changing the units
    !     from molecules per cm**3 to mass-mixing ratio.
    !
    !   x = density of air, y = mol weight in gramm

    USE mo_kind,     ONLY: dp
    IMPLICIT NONE

    ! Function return value
    REAL(dp) :: ctox

    ! Scalar arguments
    REAL(dp) :: x, y

    ctox=y/(6.022e+20_dp*x)

  END FUNCTION ctox
  
  
  END MODULE mo_ham_chemistry
  
