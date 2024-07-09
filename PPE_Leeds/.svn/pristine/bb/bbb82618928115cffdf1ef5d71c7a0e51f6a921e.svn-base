!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_drydep.f90
!!
!! \brief
!! This module handles all input terms required for the calculations
!! of the surface aerosol dry deposition in HAM.
!!
!! \author Martin Schultz (MPIfM)
!! \author Hans-Stefan Bauer (MPIfM)
!!
!! \responsible_coder
!! Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz and H.-S. Bauer (MPIfM) - original code (2000-07)
!!   -# L. Ganzeveld and A. Rhodin () - (2001-10)
!!   -# P. Stier (MPIfM) - (2002-2006)
!!   -# M. Schultz (FZJ) and S. Ferrachat (ETHZ) - new module structure (2009-11)
!!   -# G. Frontoso (C2SM) - distinguish between land / water / ice to
!!                           account for non-linearity (2012-02)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!
!! \limitations
!! None
!!
!! \details
!! This module handles all input terms required for the calculations
!! of the surface aerosol dry deposition in HAM. In particular, 
!! the routine ham_vdaer calculates the dry deposition velocities for
!! HAM aerosol tracers.
!! Most of parameters are vegetation and soil data derived from satellite data
!! a high-resolution geograhical databases. For more details see the
!! routine where the actual reading occurs. The data are monthly mean
!! values, whenever there is an annual cycle in the data. In this 
!! module, there is also the initialisation of the surface resistances
!! used for the dry deposition calculations.
!! Individual subroutines for calculation of vd by mode and integrated
!!
!! \bibliographic_references
!! None
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

MODULE mo_ham_drydep

  USE mo_kind,              ONLY: dp
  USE mo_tracdef,           ONLY: ntrac
  USE mo_ham_streams,       ONLY: rwet, densaer
 
  IMPLICIT NONE
  !----------------
  ! Public entities
  !----------------
  PRIVATE

  PUBLIC :: ham_vdaer       ! 'big leaf' aerosol dry deposition scheme
  PUBLIC :: ham_vd_presc    ! prescribed drydep velocities for HAM tracers

  LOGICAL, PARAMETER :: lwhitecap  = .TRUE.  ! consider role of whitecaps
  LOGICAL, PARAMETER :: ldr_rh     = .FALSE. ! consider relative hum. effect
  LOGICAL, PARAMETER :: lvd_bymode = .TRUE.  ! calculate drydep velocities
                                             ! individually by mode or as integrated
                                             ! quantity.



CONTAINS

SUBROUTINE ham_vdaer( kproma,   kbdim,   klev,    krow,     loland,  pvgrat, &
                      pcvs,     pcvw,    pfri,    pcvbs,    pfrw,    pum1,   &  
                      pvm1, pustarl, pustarw, pustari, pustveg, pustslsn,    &
                      pu10, pv10,                                            &
                      paz0w, paz0i, prahwat, prahice, prahveg, prahslsn,     &
                      ptslm1, prh, pxtm1, pvd                                )

  !=====================================================================
  !    Program to calculate the dry deposition velocity for
  !    aerosols as a function of the particle radius considering the
  !    impaction and difussion. Sedimentation is also calculated but not
  !    included in the calculated deposition velocity since this process
  !    is considered in a separate routine
  !=====================================================================

  !   Author:
  !   -------
  !   Laurens Ganzeveld, MPI Mainz                           04/2001
  !
  !   Modifications:
  !   --------------
  !   Philip Stier,      MPI Hamburg (adaption to ECHAM/M7,
  !                                   non-integrated scheme) 12/2001
  !   Sylvaine Ferrachart, Martin Schultz: new module structure  11/2009
  !
  !   Methodology:
  !
  !     This program calculates the integrated deposition velocity from
  !     the mass size distribution for a log normal aerosol distribution
  !     defined by the radius and the log sigma. The model calculates the
  !     Vd over land and over sea considering the diffusion,impaction and
  !     sedimentation. Over sea the effect op particle growth due to the
  !     large relative humidity is accounted for and the effect of
  !     bubble bursting is also considered. The bubble bursting enhances
  !     the dry deposition since it causes the breakdown of the laminar
  !     boundary layer and the scavenging of the particles by the sea spray.
  !     This model version does not contain yet a parameterization which
  !     specifically considers the deposition to vegetated surfaces. For
  !     these surfaces, a surface resistance as a function of canopy
  !     structure should be incorporated. This might possible be implemented
  !     in the future (Laurens Ganzeveld, 29-01-2002)
  !
  !     The model requires as input parameters:
  !     -------------------------------------------------------------------
  !     - UM, windspeed at reference height
  !     - PAZ0M, surface roughness, over sea this term is calculated from UM
  !     - PTSLM1, Air or surface temperature
  !     - LOLAND, land-sea mask
  !     - N, aerosol number [cm -3]
  !     - R, aerosol radius [m]
  !     - LSIGMA, LOG sigma of the log normal distribution
  !     - RHOA_AEROS, the density of the aerosol
  !                   replaced by: densaer (Philip Stier)
  !========================================================================
  !     LG- the parameter names in M7 for the parameters R, LSIGMA, and RHOA 
  !         are pm6rp, sigma (are similar for all the modes of M7, 
  !         LG 23-01-2002) and densaer 
  ! =======================================================================      

  USE mo_ham_m7ctl,          ONLY: cmr2mmr, sigmaln
  USE mo_ham,                ONLY: nclass,nham_subm,HAM_BULK, HAM_M7, HAM_SALSA
  USE mo_tracdef,            ONLY: trlist, AEROSOLNUMBER, AEROSOLMASS  
  USE mo_math_constants,     ONLY: pi
  USE mo_physical_constants, ONLY: grav

  IMPLICIT NONE

  !--- input parameters from echam, the parameter pevap has still to be
  !    checked more carefully, if this is a used constant or actually the
  !    echam evaporation

  INTEGER, INTENT(in)     :: kproma, kbdim, klev, krow

  LOGICAL, INTENT(in)     :: loland(kbdim)

  REAL(dp), INTENT(in)    :: pvgrat(kbdim),     & 
                             pcvs(kbdim),       &   
                             pcvw(kbdim),       &     
                             pfri(kbdim),       &   
                             pcvbs(kbdim),      &    
                             pfrw(kbdim),       &   
                             pum1(kbdim,klev),  &
                             pvm1(kbdim,klev),  &
                             pustarl(kbdim),    &
                             pustarw(kbdim),    &
                             pustari(kbdim),    &
                             pustveg(kbdim),    &
                             pustslsn(kbdim),   & 
                             pu10(kbdim),       &     
                             pv10(kbdim),       &   
                             paz0w(kbdim),      &
                             paz0i(kbdim),      &
                             prahwat(kbdim),    &
                             prahice(kbdim),    &
                             prahveg(kbdim),    &  
                             prahslsn(kbdim),   & 
                             ptslm1(kbdim),     &   
                             prh(kbdim)

!gf  REAL(dp), INTENT(in)    :: pevap(kbdim)

  REAL(dp), INTENT(in)    :: pxtm1(kbdim,klev,ntrac)

  REAL(dp), INTENT(inout) :: pvd(kbdim,ntrac)


  !--- stream variable references
  REAL(dp),            POINTER :: rwet_p(:,:,:)
  REAL(dp),            POINTER :: densaer_p(:,:,:)

  !--- Local Variables:

  INTEGER :: ji, jl, jt, jclass, jtype

  INTEGER, PARAMETER :: NNUMBER=1, MASS=2

  !    Aerosol properties and parameters relevant to the calculations:

  REAL(dp)    :: rx1,rx2,um10,s,sc,cunning,vb_veg,vb_slsn,                            &
                 vim_veg, vim_slsn, st_veg, st_slsn, st_wat, st_ice, rew, vb_sea,     & 
                 vb_ice, vim_sea, vim_ice, vkdaccsea,                                 &
                 eff, rdrop, zdrop, qdrop, vkd_veg, vkd_slsn, vkd_wat, vkd_ice,       &
                 vkc_veg, vkc_slsn, vkc_wat, vkc_ice, vdpart_veg, vdpart_slsn,        &
                 vdpart_wat, vdpart_ice, dc, relax,                                   &
                 eps,phi,alpha1,alpharat,vk1,vk2,zm,zmvd,zn,zdm,zdn,zdlnr

  ! mz_lg_20031014+ added for some modifications of the vegetation surface
  !     resistance based on the paper by Gallagher et al., JGR, 2002. The
  !     paper shows some more details about the parameterization by 
  !     Slinn, Atmos. Environment, Vol 16, 1785-1794, 1982
  REAL(dp)    :: vin_veg, zrebound
  ! mz_lg_20031014-

  REAL(dp)    :: zrint
  REAL(dp)    :: zvdrydep(kbdim)            ! Aerosol Dry Deposition Velocity defined for each mode
  REAL(dp)    :: zr(kbdim)                  ! Aerosol Radius of the respective tracer at klev
  REAL(dp)    :: znmr(kbdim)                ! Aerosol Number Mixing ratio at t=t+1

  REAL(dp)    :: zevap(kbdim) !gf

  !--- Auxiliary fields:

  REAL(dp)    :: um(kbdim),     alpha(kbdim), beta(kbdim), &
             alphae(kbdim), bubble(kbdim)

  LOGICAL :: ldrydep

  INTEGER, PARAMETER  :: nint=199           ! number of integration intervals
  REAL(dp)            :: zintstep(nint)     ! integration intervals

  REAL(dp), PARAMETER :: crmin=0.01E-6_dp   ! smallest radius for which a dry deposition 
                                            ! velocity is calculated [cm]

  !--- Assign values to used constants:

  REAL(dp), PARAMETER :: fln10=2.302585_dp
  REAL(dp), PARAMETER :: w2pi=2.506638_dp
  REAL(dp), PARAMETER :: g=grav*1.e2_dp      ! cm s-2 at sea level (SF #369: changed to model-wide constant)
  REAL(dp), PARAMETER :: dynvisc=1.789e-4_dp ! g cm-1 s-1
  REAL(dp), PARAMETER :: cl=0.066*1e-4_dp    ! mean free path [cm] (particle size also in cm)
  REAL(dp), PARAMETER :: bc= 1.38e-16_dp     ! boltzman constant [g cm-2 s-1 K-1] (1.38e-23 J deg-1)
  REAL(dp), PARAMETER :: kappa=1._dp         ! shapefactor
  REAL(dp), PARAMETER :: visc=0.15_dp        ! molecular viscocity [cm2 s-1]
  REAL(dp), PARAMETER :: daccm=0._dp         ! factor which corrects for evaporation (see paper slinn)
  REAL(dp), PARAMETER :: vkar=0.40_dp        ! von karman constant

  ! mz_lg_20031014+ added for the Vd aerosol over vegetation
  ! mz_lg_20040602+ modified
  REAL(dp), PARAMETER :: zAS    = 10.E-6_dp*1.E2_dp  ! um -> CM, see paper Gallagher and Slinn, 1982
  ! here the smallest collector size is set at 10 um
  ! Those are the values for Slinn's 82 model (see Table 1)
  REAL(dp), PARAMETER :: zAL    = 1.E-3_dp*1.E2_dp   ! mm -> CM, see paper Gallagher and Slinn, 1982
  ! here the largest collector size is set at 1 mm  
  ! mz_lg_20040602-+ modified
  ! mz_lg_20031014-

  !@@@ Check with Laurens (see comments above)

  zevap(1:kproma)=0._dp

  !--- Integration stepsizes [m], 0-10 and 10-100,100-1000 um
  !    the first interval is divided into 1000 steps of 0.1 um
  !    whereas for the interval >10 um 90 steps of 1 um are selected
  !    and for the radius > 100 um 90 steps of 10 um are selected

  DATA zintstep /100*.01E-6_dp,90*1.E-6_dp,9*10.E-6_dp/

  !--- 1) Calculate correction factors:

  DO jl=1,kproma

     ! LG- calculation of some parameters required for the deposition calculations

     um(jl)=MAX(0.001_dp,SQRT(pum1(jl,klev)**2+pvm1(jl,klev)**2))
     um10=MAX(0.001_dp,SQRT(pu10(jl)**2+pv10(jl)**2))

     ! bubble bursting effect,see Hummelshoj, equation 10
     ! relationship by Wu (1988), note that Hummelshoj has not
     ! considered the cunningham factor which yields a different
     ! vb curve, with smaller values for small particles

     IF (lwhitecap) THEN
        alpha(jl)=MAX(1.e-10_dp,1.7e-6_dp*um10**3.75_dp)     ! old 10 m windspeed !!
        eff=0.5_dp
        rdrop=0.005_dp        ! cm
        zdrop=10.0_dp         ! cm
        qdrop=5._dp*(100._dp*alpha(jl))
        bubble(jl)=((100._dp*pustarw(jl))**2)/(100._dp*um(jl))+eff* &
                   (2._dp*pi*rdrop**2)*(2._dp*zdrop)*(qdrop/alpha(jl))
     ELSE
        alpha(jl)=0._dp
        bubble(jl)=0._dp
     ENDIF

     !--- Correction of particle radius for particle growth close to the
     !    surface according to Fitzgerald, 1975, the relative humidity over
     !    the ocean is restricted to 98% (0.98) due to the salinity

     s=MIN(0.98_dp,prh(jl))
     eps=0.6_dp
     beta(jl)=EXP((0.00077_dp*s)/(1.009_dp-s))
     phi=1.058_dp-((0.0155_dp*(s-0.97_dp))/(1.02_dp-s**1.4_dp))
     alpha1=1.2_dp*EXP((0.066_dp*s)/(phi-s))
     vk1=10.2_dp-23.7_dp*s+14.5_dp*s**2
     vk2=-6.7_dp+15.5_dp*s-9.2_dp*s**2
     alpharat=1._dp-vk1*(1._dp-eps)-vk2*(1._dp-eps**2)
     alphae(jl)=alpharat*alpha1

     !--- Over land no correction for large humidity close to the surface:

     IF (loland(jl).OR..NOT.ldr_rh) THEN
        alphae(jl)=1._dp
        beta(jl)=1._dp
     ENDIF

  END DO! jl=1, kproma

  !--- 2)Calculate dry deposition velocity for each internally mixed mode:

  IF (lvd_bymode) THEN
    CALL ham_vdaer_bymode
  ELSE
    CALL ham_vdaer_integrated
  END IF

           

  CONTAINS

!---------------------------------------------------------------------------------------------
!  calculation by aerosol mode

  SUBROUTINE ham_vdaer_bymode

     DO jclass=1, nclass

       rwet_p    => rwet(jclass)%ptr
       densaer_p => densaer(jclass)%ptr

       DO jtype=1, 2       ! number (= 1) and mass (= 2)

         !--- Calculations are done in [cm] to get results in [cm s-1]:

         SELECT CASE(jtype)
            CASE(1)
                ldrydep=ANY(trlist%ti(:)%ndrydep==2 .AND. trlist%ti(:)%mode==jclass    &
                           .AND. trlist%ti(:)%nphase==AEROSOLNUMBER)
                zr(1:kproma)=rwet_p(1:kproma,klev,krow)*1.E2_dp
            CASE(2)
                ldrydep=ANY(trlist%ti(:)%ndrydep==2 .AND. trlist%ti(:)%mode==jclass    &
                           .AND. trlist%ti(:)%nphase==AEROSOLMASS)
            SELECT CASE(nham_subm)
               CASE(HAM_M7)
                   zr(1:kproma)=rwet_p(1:kproma,klev,krow)*1.E2_dp*cmr2mmr(jclass)
               CASE(HAM_SALSA) 
                   !Deposition size is the same for mass and number in SALSA
                   zr(1:kproma)=rwet_p(1:kproma,klev,krow)*1.E2_dp
            END SELECT
         END SELECT

         IF(ldrydep) THEN

           DO jl=1, kproma

             !--- Do calculations only for particles larger than 0.0001 um :

             IF (zr(jl) > crmin) THEN

                !--- Cunningham factor:

                cunning=1._dp+(cl/(alphae(jl)*zr(jl))**beta(jl))                 &
                              *(2.514_dp+0.800_dp*EXP(-0.55_dp*(zr(jl))/cl))

                !--- Diffusivity:

                dc=(bc*ptslm1(jl)*cunning)/(3._dp*pi*dynvisc*(alphae(jl)*zr(jl))**beta(jl))

                ! Relaxation:

                ! mz_lg_20021311, modified based on discussions with
                ! Philip Stier [kg m-3 => g cm-3]

                relax=(densaer_p(jl,klev,krow)*1.E-3_dp*(((alphae(jl)*zr(jl))**beta(jl))**2)* &
                          cunning)/(18._dp*dynvisc*kappa)

                ! Calculation of schmidt and stokes number

                sc      =visc/dc
                st_veg  =(relax*(100._dp*pustveg(jl))**2)/visc

                ! mz_lg_20031014+, modified calculation of the Stokes number
                !     over vegetated surfaces, see paper by Gallagher et al., 
                !     JGR 2002. zAL is a characteristic radius for the
                !     largest collectors comprising the surface
                ! ham_ps_20031201+ included minumum value for Stokes numbers:
                st_veg  =MAX((relax*(100._dp*pustveg(jl))**2)/(g*zAL),1.E-1_dp)
                ! mz_lg_20031014-
                st_slsn =MAX((relax*(100._dp*pustslsn(jl))**2)/visc,1.E-1_dp)
                st_wat  =MAX((relax*(100._dp*pustarw(jl))**2)/visc,1.E-1_dp)
                st_ice  =MAX((relax*(100._dp*pustari(jl))**2)/visc,1.E-1_dp)
                ! ham_ps_20031201-

                !--- Calculation of the dry deposition velocity
                !    See paper slinn and slinn, 1980, vd is related to d**2/3
                !    over land, whereas over sea there is accounted for slip
                !    vb represents the contribution in vd of the brownian diffusion
                !    and vi represents the impaction.
                !
                !--- Over land, the vegetation and snow and bare soil fractions
                !    are considered:

                vb_veg   =(1._dp/vkar)*((pustveg(jl)/um(jl))**2)*100._dp*um(jl)*(sc**(-2._dp/3._dp))
                vb_slsn  =(1._dp/vkar)*((pustslsn(jl)/um(jl))**2)*100._dp*um(jl)*(sc**(-2._dp/3._dp))
                vb_ice  =(1._dp/vkar)*((pustari(jl)/um(jl))**2)*100._dp*um(jl)*(sc**(-2._dp/3._dp))
                vim_veg  =(1._dp/vkar)*((pustveg(jl)/um(jl))**2)*100._dp*um(jl)     &
                          *(10._dp**(-3._dp/st_veg))

                ! mz_lg_20031014+, modified calculation of the impaction over
                !     vegetated surfaces, see paper by Gallagher et al., JGR 2002. 
                !     We have applied here the parameterization by Slinn [1982] 
                !     over vegetated surfaces
                vim_veg   =(1._dp/vkar)*((pustveg(jl)/um(jl))**2)*100._dp*um(jl)* &
                           (st_veg**2/(1._dp+st_veg**2))
                ! mz_lg_20031014-

                vim_slsn  =(1._dp/vkar)*((pustslsn(jl)/um(jl))**2)*100._dp*um(jl)       &
                           *(10._dp**(-3._dp/st_slsn))
                vim_ice  =(1._dp/vkar)*((pustari(jl)/um(jl))**2)*100._dp*um(jl)         &
                           *(10._dp**(-3._dp/st_ice))

                ! LG- the term evap has not been defined yet (30-01-2002) and the
                !     term daccm is set to zero anyhow. It still must be checked if
                !     this term should be included and how it relates to the factors
                !     that correct for particle growth close to the surface, according
                !     to Fitzjarald, 1975 (see ALPHAE and BETA)

                ! mz_lg_20031014+, modified calculation of the surface resistance over 
                !     vegetated surfaces, see paper by Gallagher et al., JGR 2002. 
                !     The calculation includes the interception collection efficiency 
                !     vim and a rebound correction factor R
                !vkd_veg  =pevap(jl)*daccm+(vb_veg+vim_veg)
                vin_veg  =(1._dp/vkar)*((pustveg(jl)/um(jl))**2)*100._dp*um(jl)*      &
                          (1._dp/2._dp)*(zr(jl)/zAS)**2          ! equation 18
                zrebound =EXP(-st_veg**0.5_dp)                ! equation 19
                vkd_veg  =zevap(jl)*daccm+zrebound*(vb_veg+vim_veg+vin_veg)
                ! mz_lg_20031014-

                vkd_slsn =zevap(jl)*daccm+(vb_slsn+vim_slsn)
                vkd_ice  =zevap(jl)*daccm+(vb_ice+vim_ice)

                !--- Over sea:
                !    Brownian diffusion for rough elements, see Hummelshoj
                !    re is the reynolds stress:

                rew       =(100._dp*pustarw(jl)*100._dp*paz0w(jl))/visc
                vb_sea    =(1._dp/3._dp)*100._dp*pustarw(jl)*((sc**(-0.5_dp))*rew**(-0.5_dp))
                vim_sea   =100._dp*pustarw(jl)*10._dp**(-3._dp/st_wat)
                vkdaccsea =vb_sea+vim_sea
                vkd_wat   =(1._dp-alpha(jl))*vkdaccsea+alpha(jl)*(bubble(jl))

                ! Slinn and Slinn parameterization, calculation of vd:
                ! ====================================================
                ! LG- without considering the role of sedimentation !
                !     This process is being calculated in a separate
                !     subroutine in echam5!
                ! ====================================================

                vkc_veg     =(100._dp/prahveg(jl))
                vkc_slsn    =(100._dp/prahslsn(jl))
                vkc_wat     =(100._dp/prahwat(jl))
                vkc_ice     =(100._dp/prahice(jl))

                !@@@ Included security checks for small values of vxx_yyy:

                IF(vkc_veg > 1.E-5_dp .AND. vkd_veg > 1.E-5_dp) THEN
                   vdpart_veg  =1._dp/((1._dp/vkc_veg)+(1._dp/vkd_veg))
                ELSE
                   vdpart_veg  = 0._dp
                END IF
                IF(vkc_slsn > 1.E-5_dp .AND. vkd_slsn > 1.E-5_dp) THEN 
                   vdpart_slsn = 1._dp/((1._dp/vkc_slsn)+(1._dp/vkd_slsn))
                ELSE
                   vdpart_slsn = 0._dp
                END IF
                IF(vkc_wat > 1.E-5_dp .AND. vkd_wat > 1.E-5_dp) THEN
                   vdpart_wat  =1._dp/((1._dp/vkc_wat)+(1._dp/vkd_wat))
                ELSE
                   vdpart_wat = 0._dp
                END IF
                IF(vkc_ice > 1.E-5_dp .AND. vkd_ice > 1.E-5_dp) THEN
                   vdpart_ice  =1._dp/((1._dp/vkc_ice)+(1._dp/vkd_ice))
                ELSE
                   vdpart_ice = 0._dp
                END IF

                !--- Calculate the dry deposition velocity weighted according to the surface types:
                !   pcvs:   snow fraction
                !   pcvbs:  bare soil fraction
                !   pvgrat: vegetation ratio
                !   pcvw:   wet skin fraction
                !   pfri:   (sea) ice fraction
                !   pfrw:   open water fraction

                zvdrydep(jl) = pcvs(jl)*vdpart_slsn                                       & 
                             + pcvbs(jl)*vdpart_slsn                                      &
                             + (1._dp-pcvs(jl))*(1._dp-pcvw(jl))*pvgrat(jl)*vdpart_veg    &
                             + (1._dp-pcvs(jl))*pcvw(jl)*vdpart_veg                       &
                             + pfri(jl)*vdpart_ice                                        &
                             + pfrw(jl)*vdpart_wat    

             ELSE

                zvdrydep(jl)=0._dp

             END IF !zr(jl) > crmin

          END DO ! jl=1, kproma

          !--- save drydep velocity
          SELECT CASE(jtype)
          CASE(1)
            DO jt=1,ntrac
              IF (trlist%ti(jt)%ndrydep>0  .AND.     &
                  trlist%ti(jt)%mode==jclass .AND.     &
                  trlist%ti(jt)%nphase==AEROSOLNUMBER)THEN

                 pvd(1:kproma,jt) = zvdrydep(1:kproma)*1.E-2_dp  ! conversion from cm s-1 to m s-1
      
              END IF
            END DO
          CASE(2)
            DO jt=1,ntrac
              IF (trlist%ti(jt)%ndrydep>0  .AND.   &
                  trlist%ti(jt)%mode==jclass .AND.   &
                  trlist%ti(jt)%nphase==AEROSOLMASS)THEN

                pvd(1:kproma,jt) = zvdrydep(1:kproma)*1.E-2_dp  ! conversion from cm s-1 to m s-1

              END IF
            END DO
          END SELECT
      
        END IF ! ldrydep

      END DO ! jtype=1, 2

    END DO ! jclass=1, nclass

  END SUBROUTINE ham_vdaer_bymode


!---------------------------------------------------------------------------------------------
!  integrated calculation

  SUBROUTINE ham_vdaer_integrated

     DO jl=1, kproma

       DO jtype=1, 2       ! number (= 1) and mass (= 2)

         DO jclass=1, nclass

           rwet_p    => rwet(jclass)%ptr
           densaer_p => densaer(jclass)%ptr

           ! mz_lg_20030521+, small change in criteria to determine if calculations
           !     should be bypassed, see also the same criteria for lvd_integrated=.false.

           IF ( ANY(trlist%ti(:)%mode==jclass.AND.trlist%ti(:)%ndrydep>0) .AND. &
                rwet_p(jl,klev,krow) > crmin) THEN

             !--- Initialisations:

             DO jt=1, ntrac
               IF (trlist%ti(jt)%mode==jclass .AND. trlist%ti(jt)%nphase==AEROSOLNUMBER ) THEN
                 znmr(jl)=pxtm1(jl,klev,jt)
               END IF
             END DO

             SELECT CASE(jtype)
             CASE(1)
               zr(jl)=rwet_p(jl,klev,krow)
             CASE(2)
               zr(jl)=rwet_p(jl,klev,krow)*cmr2mmr(jclass)
             END SELECT

             rx1  = 0.0_dp
             zmvd = 0.0_dp
             zm   = 0.0_dp

             DO ji=1,nint-1

               ! Integration size interval is rx2-rx1:

               rx1=rx1+zintstep(ji)
               rx2=rx1+zintstep(ji+1)

               !--- Calculate the number of particles in the intervall r+dr
               !    for a log-normal distribution with the moments:
               !    N=znmr, sigma, CountMeanRadius=pm6rp
               !       dn=n(ln(r))dln(r) [m-3]
               !    where:
               !       r = (rx1+rx2)/2
               !    and
               !       n(ln(r)) = N/(sqrt(2 PI)*ln(sigma)) * exp(-(ln(r)-ln(CMR))/(2ln(sigma)**2)

               zrint     = (rx1+rx2)/2._dp

               !--- Particle Numbers:

               zdlnr    = LOG(rx2)-LOG(rx1)

               zn       = znmr(jl)/(2._dp*pi*sigmaln(jclass)) * &
                          EXP(-(LOG(zrint)-LOG(zr(jl)))/(2._dp*sigmaln(jclass)**2))
               zdn      = zn * zdlnr

               !--- Average particle mass:
               !    in intervall (zrint-dr,zrint+dr) [kg(tracer)/kg(air)]

               zdm  = 4._dp/3._dp*pi*zrint**3 * densaer_p(jl,klev,krow) * zdn

               !--- Cunningham factor:
               !    (Philip Stier, 11-2001, Separation into dry and wet not needed
               !     as the ambient radius pm6rp incorporates the particle growth)
               ! =================================================================
               ! LG- It still must be checked if it is indeed allright to
               !     calculate only the parameters for the aerosol modal
               !     parameters given by M7 to include the effect of particle
               !     growth. The parameterization of the effect of particle growth
               !     that was included in the original scheme of the aerosol
               !     dry deposition module corrects for the growth close to the
               !     surface for a large relative humidity in quasi-laminar layer,
               !     which is actual a sub-grid scale effect occuring in the
               !     surface layer of echam. We have included again the calculation
               !     of the parameters that are used to introduce this effect,
               !     ALPHAE and BETA and use of this terms can easily be considered
               !     or ignored, for the latter option setting the calculated values
               !     to 1 (this is default done for over-land grids)
               !
               !     Laurens Ganzeveld, 30-01-2002
               ! =================================================================

               !--- Convert radii to cm to get the results in cm s-1:

               cunning=1._dp+(cl/(alphae(jl)*zrint*1.E2_dp)**beta(jl))       &
                             *(2.514_dp+0.800_dp*EXP(-0.55_dp*(zrint*1.E2_dp)/cl))

               !--- Diffusivity:

               dc=(bc*ptslm1(jl)*cunning)/(3._dp*pi*dynvisc*(alphae(jl)*zrint*1.E2_dp)**beta(jl))

               ! Relaxation:

               ! mz_lg_20021311, modified based on discussions with
               ! Philip Stier [kg m-3 => g cm-3]

               relax=(densaer_p(jl,klev,krow)*1.E-3_dp*(((alphae(jl)*zrint)**beta(jl))**2)* &
                     cunning)/(18._dp*dynvisc*kappa)

               ! Calculation of schmidt and stokes number

               sc      =visc/dc
               st_veg  =(relax*(100._dp*pustveg(jl))**2)/visc

               ! mz_lg_20031014+, modified calculation of the Stokes number
               !     over vegetated surfaces, see paper by Gallagher et al., 
               !     JGR 2002. zAL is a characteristic radius for the
               !     largest collectors comprising the surface
               st_veg  =(relax*(100._dp*pustveg(jl))**2)/(g*zAL)
               ! mz_lg_20031014-

               st_slsn =(relax*(100._dp*pustslsn(jl))**2)/visc
               st_wat  =(relax*(100._dp*pustarw(jl))**2)/visc
               st_ice  =(relax*(100._dp*pustari(jl))**2)/visc

               !--- Calculation of the dry deposition velocity
               !    See paper slinn and slinn, 1980, vd is related to d**2/3
               !    over land, whereas over sea there is accounted for slip
               !    vb represents the contribution in vd of the brownian diffusion
               !    and vi represents the impaction.
               !
               !--- Over land, the vegetation and snow and bare soil fractions
               !    are considered:

               vb_veg   =(1._dp/vkar)*((pustveg(jl)/um(jl))**2)*100._dp*um(jl)*(sc**(-2._dp/3._dp))
               vb_slsn  =(1._dp/vkar)*((pustslsn(jl)/um(jl))**2)*100._dp*um(jl)*(sc**(-2._dp/3._dp))
               vb_ice  =(1._dp/vkar)*((pustari(jl)/um(jl))**2)*100._dp*um(jl)*(sc**(-2._dp/3._dp))

               ! mz_lg_20031014+, modified calculation of the impaction over
               !     vegetated surfaces, see paper by Gallagher et al., JGR 2002. 
               !     We have applied here the parameterization by Slinn [1982] 
               !     over vegetated surfaces
               !     (Philip Stier 30/10/2003, corrected um to um(jl))
               !vim_veg   =(1./vkar)*((pustveg(jl)/um(jl))**2)*100.*um(jl)*(10.**(-3./st_veg))

               vim_veg   =(1._dp/vkar)*((pustveg(jl)/um(jl))**2)*100._dp*um(jl)      &
                          *(st_veg**2/(1._dp+st_veg**2))
               ! mz_lg_20031014-

               vim_slsn  =(1._dp/vkar)*((pustslsn(jl)/um(jl))**2)*100._dp*um(jl)     &
                          *(10._dp**(-3._dp/st_slsn))
               vim_ice  =(1._dp/vkar)*((pustari(jl)/um(jl))**2)*100._dp*um(jl)     &
                          *(10._dp**(-3._dp/st_ice))

               ! LG- the term evap has not been defined yet (30-01-2002) and the
               !     term daccm is set to zero anyhow. It still must be checked if
               !     this term should be included and how it relates to the factors
               !     that correct for particle growth close to the surface, according
               !     to Fitzjarald, 1975 (see ALPHAE and BETA)


               ! mz_lg_20031014+, modified calculation of the surface resistance over 
               !     vegetated surfaces, see paper by Gallagher et al., JGR 2002. 
               !     The calculation includes the interception collection efficiency 
               !     vim and a rebound correction factor R
               !vkd_veg  =pevap(jl)*daccm+(vb_veg+vim_veg)

               vin_veg  =(1._dp/vkar)*((pustveg(jl)/um(jl))**2)*100._dp*um(jl)*      &
                    (1._dp/2._dp)*(zrint/zAS)**2           ! equation 18
               zrebound =EXP(-st_veg**0.5_dp)            ! equation 19
               vkd_veg  =zevap(jl)*daccm+zrebound*(vb_veg+vim_veg+vin_veg)
               ! mz_lg_20031014-

               vkd_slsn =zevap(jl)*daccm+(vb_slsn+vim_slsn)

               !--- Over sea:
               !    Brownian diffusion for rough elements, see Hummelshoj
               !    re is the reynolds stress:

               rew        =(100._dp*pustarw(jl)*100._dp*paz0w(jl))/visc
               vb_sea     =(1._dp/3._dp)*100._dp*pustarw(jl)*((sc**(-0.5_dp))*rew**(-0.5_dp))
               vim_sea     =100._dp*pustarw(jl)*10._dp**(-3._dp/st_wat)
               vkdaccsea =vb_sea+vim_sea
               vkd_wat   =(1._dp-alpha(jl))*vkdaccsea+alpha(jl)*(bubble(jl))

               ! Slinn and Slinn parameterization, calculation of vd:
               ! ====================================================
               ! LG- without considering the role of sedimentation !
               !     This process is being calculated in a separate
               !     subroutine in echam5!
               ! ====================================================

               vkc_veg     =(100._dp/prahveg(jl))
               vkc_slsn    =(100._dp/prahslsn(jl))
               vkc_wat     =(100._dp/prahwat(jl))
               vkc_ice     =(100._dp/prahice(jl))
               vdpart_veg  =1._dp/((1._dp/vkc_veg)+(1._dp/vkd_veg))
               vdpart_slsn =1._dp/((1._dp/vkc_slsn)+(1._dp/vkd_slsn))
               vdpart_wat  =1._dp/((1._dp/vkc_wat)+(1._dp/vkd_wat))
               vdpart_ice  =1._dp/((1._dp/vkc_ice)+(1._dp/vkd_ice))

               !--- Calculate sum(dm*vd(dm)):
               !   pcvs:   snow fraction
               !   pcvbs:  bare soil fraction
               !   pvgrat: vegetation ratio
               !   pcvw:   wet skin fraction
               !   pfri:   (sea) ice fraction
               !   pfrw:   open water fraction


               zmvd = zmvd + zdm*( &
                      pcvs(jl)*vdpart_slsn+                                    & ! snow fraction
                      pcvbs(jl)*vdpart_slsn+                                   & ! bare soil fraction
                      (1._dp-pcvs(jl))*(1._dp-pcvw(jl))*pvgrat(jl)*vdpart_veg+ & ! vegetation fraction
                      (1._dp-pcvs(jl))*pcvw(jl)*vdpart_veg+                    & ! wet skin fraction
                      pfri(jl)*vdpart_ice+                                     & ! sea ice fraction
                      pfrw(jl)*vdpart_wat                                      ) ! water fraction

               !--- Calculate m=sum(dm) [kg(tracer)/kg(air)]:

               zm   = zm   + zdm

             END DO !ji=1,nint-1

             !--- Calculate integrated dry deposition velocity [cm s-1]:
             !    vd = sum(dm*vd) / m

             zvdrydep(jl)= zmvd/zm

           ELSE

             zvdrydep(jl)=0._dp
  
           END IF

         END DO ! jl=1, kproma

         !--- Assign the calculated deposition velocity to all tracers
         !    in the respecive mode and re-convert it to [m s-1]:
        
         ! this can cause the same fields to be written many times if drydep diag is BYTRACER
         ! needs improvement

         SELECT CASE(jtype)
         CASE(1)
           DO jt=1,ntrac
             IF (trlist%ti(jt)%ndrydep>0  .AND.     &
                 trlist%ti(jt)%mode==jclass .AND.     &
                 trlist%ti(jt)%nphase==AEROSOLNUMBER)THEN

               pvd(1:kproma,jt) = zvdrydep(1:kproma)*1.E-2_dp   ! conversion cm s-1 to m s-1 
             END IF
           END DO
         CASE(2)
           DO jt=1,ntrac
             IF (trlist%ti(jt)%ndrydep>0  .AND.   &
                 trlist%ti(jt)%mode==jclass .AND.   &
                 trlist%ti(jt)%nphase==AEROSOLMASS)THEN

               pvd(1:kproma,jt) = zvdrydep(1:kproma)*1.E-2_dp   ! conversion cm s-1 to m s-1 
             END IF
           END DO
         END SELECT

       END DO ! jclass=1, nclass

     END DO ! jtype=1, 2

  END SUBROUTINE ham_vdaer_integrated


END SUBROUTINE ham_vdaer

!------------------------------------------------------------------------------

SUBROUTINE ham_vd_presc(kproma, kbdim,  klev,    krow,  loland,   &
                        paphp1, pcvs,   pforest, pfri,  ptsi,     &
                        pcvw,   ptslm1, pws,     pwsmx, pdensair, &
                        pvd                                        )
  
  ! Purpose:
  ! ---------
  ! This routine prescribes dry deposition velocities for
  ! aerosols and gases in dependency of surface type, etc.
  ! ...
  ! 
  ! Authors:
  ! ----------
  ! Hans Feichter, MPI-MET
  ! Philip Stier,  MPI-MET
  !
  ! Method:
  ! -------
  ! Currently rudimentary!  
  !
  ! Interface:
  ! ----------
  ! *xt_vdrydep_presc* is called from *xt_drydep*
  
  USE mo_kind,              ONLY: dp
  USE mo_control,           ONLY: nlev
  USE mo_time_control,      ONLY: time_step_len
  USE mo_physical_constants, ONLY: tmelt, grav
  USE mo_tracdef,           ONLY: trlist, ntrac, AEROSOLMASS, AEROSOLNUMBER, GAS
  USE mo_physc2,            ONLY: ctfreez, & ! sea water freezing temp.
                                 csncri     ! m water equivalent critical snow depth
  USE mo_vphysc,            ONLY: vphysc
!!  USE mo_submodel_diag,     ONLY: get_diag_pointer
  
  
  IMPLICIT NONE  

  !--- Parameters:
  !--- Local variables:
  !
  ! zvdrd   dry deposition velocity in m/s
  ! zvdrd(jl,1)  for so2 gas 
  ! zvdrd(jl,2)  for aerosols
  ! zvdcoarse    for dust and seasalt aerosols (coarse mode)

  ! Parameters:

  INTEGER, INTENT(in)     :: kproma, kbdim, klev, krow

  REAL(dp), INTENT(in)    :: pfri(kbdim),         &
                             pcvs(kbdim),         &
                             pcvw(kbdim),         &
                             pforest(kbdim),      &
                             ptsi(kbdim),         &
                             ptslm1(kbdim),       &
                             pws(kbdim),          &
                             pwsmx(kbdim),        &
                             pdensair(kbdim),     &
                             paphp1(kbdim,klev+1)

  REAL(dp), INTENT(inout) :: pvd(kbdim,ntrac)

  ! Local scalars:

  INTEGER :: jl, jt

  REAL(dp):: zvd2ice,        zvd4ice,        zvd2nof,       zvd4nof,    &
             zvwc2,          zvw02,          zvwc4,         zvw04

  ! Local arrays:

  REAL(dp):: zmaxvdry(kbdim)
  REAL(dp):: zvdrd(kbdim,2)

  LOGICAL :: loland(kbdim)

  ! Constants:

  REAL(dp), PARAMETER :: zvdcoarse = 0.01_dp, & ! V_dep coarse mode [m/s]
                         zvdaitken = 0.01_dp    ! V_dep aitken mode [m/s]


  !--- 1) Calculation of the dry deposition velocity: -----------------------

  !--- Coefficients for zvdrd = function of soil moisture:
  !
  zvwc2=(0.8e-2_dp - 0.2e-2_dp)/(1._dp - 0.9_dp)
  zvw02=zvwc2-0.8e-2_dp
  zvwc4=(0.2e-2_dp - 0.025e-2_dp)/(1._dp - 0.9_dp)
  zvw04=zvwc4-0.2e-2_dp

  !--- Maximal deposition velocity is maximal vertical grid velocity:

  DO jl=1,kproma
    zmaxvdry(jl)=(paphp1(jl,nlev+1)-paphp1(jl,nlev))/(grav*pdensair(jl)*time_step_len)
  END DO   
             

  DO jl=1,kproma

    !--- 1.1) Over sea:

    IF (.NOT.loland(jl)) THEN
      !         - sea ice -
      !           - melting/not melting seaice-
      IF (ptsi(jl).GE.(ctfreez-0.01_dp)) THEN
        zvd2ice=0.8e-2_dp                             ! SO2
        zvd4ice=0.2e-2_dp                             ! others
      ELSE
        zvd2ice=0.1e-2_dp
        zvd4ice=0.025e-2_dp
      END IF
      zvdrd(jl,1)=(1._dp-pfri(jl))*1.0e-2_dp+pfri(jl)*zvd2ice
      zvdrd(jl,2)=(1._dp-pfri(jl))*0.2e-2_dp+pfri(jl)*zvd4ice
    ELSE

      !--- 1.2) Over land:
      !        - non-forest areas -
      !         -  snow/no snow -
      IF (pcvs(jl).GT.csncri) THEN
        !  - melting/not melting snow -
        IF (vphysc%smelt(jl,krow).GT.0._dp) THEN
          zvd2nof=0.8e-2_dp
          zvd4nof=0.2e-2_dp
        ELSE
          zvd2nof=0.1e-2_dp
          zvd4nof=0.025e-2_dp
        END IF
      ELSE
        !  -  frozen/not frozen soil -
        IF (ptslm1(jl).LE.tmelt) THEN
          zvd2nof=0.2e-2_dp
          zvd4nof=0.025e-2_dp
        ELSE
          !  - wet/dry -
          !  - completely wet -
          IF (pcvw(jl).GE.0.01_dp .OR. pws(jl).EQ.pwsmx(jl)) THEN
            zvd2nof=0.8e-2_dp
            zvd4nof=0.2e-2_dp
          ELSE
            !   - dry -
            IF ((pws(jl)/pwsmx(jl)).LT.0.9_dp) THEN
              zvd2nof=0.2e-2_dp
              zvd4nof=0.025e-2_dp
            ELSE
            !  - partly wet -
              zvd2nof=zvwc2*(pws(jl)/pwsmx(jl))-zvw02
              zvd4nof=zvwc4*(pws(jl)/pwsmx(jl))-zvw04
            END IF
          END IF
        END IF
      END IF
      zvdrd(jl,1)=pforest(jl)*0.8e-2_dp+(1._dp-pforest(jl))*zvd2nof
      zvdrd(jl,2)=pforest(jl)*0.2e-2_dp+(1._dp-pforest(jl))*zvd4nof
    END IF
    !
    !---  Ask zvdrd for maximum:
    !
    zvdrd(jl,1)=MIN(zvdrd(jl,1),zmaxvdry(jl))
    zvdrd(jl,2)=MIN(zvdrd(jl,2),zmaxvdry(jl))

  END DO

  DO jt=1, ntrac

    IF (trlist%ti(jt)%ndrydep/=1) CYCLE   ! only set vd for tracers that request prescribed values

    ! HAM gas-phase tracers
    IF (trlist%ti(jt)%basename=='SO2' .OR. &
        (trlist%ti(jt)%basename=='SO4'.AND. trlist%ti(jt)%nphase==GAS) ) THEN

      pvd(1:kproma,jt) = zvdrd(1:kproma,1)

    ! HAM aerosol tracers
    ELSE IF(trlist%ti(jt)%nphase==AEROSOLMASS .OR. trlist%ti(jt)%nphase==AEROSOLNUMBER) THEN

      SELECT CASE (trlist%ti(jt)%mode)
      CASE(1)
        pvd(1:kproma,jt) = 0._dp
      CASE(2)
        pvd(1:kproma,jt) = zvdaitken
      CASE(3)  
        pvd(1:kproma,jt) = zvdrd(1:kproma,2)
      CASE(4) 
        pvd(1:kproma,jt) = zvdcoarse
      CASE(5) 
        pvd(1:kproma,jt) = zvdaitken
      CASE(6) 
        pvd(1:kproma,jt) = zvdrd(1:kproma,2)
      CASE(7) 
        pvd(1:kproma,jt) = zvdcoarse
      END SELECT
           
    END IF
        
  END DO

END SUBROUTINE ham_vd_presc

END MODULE mo_ham_drydep
