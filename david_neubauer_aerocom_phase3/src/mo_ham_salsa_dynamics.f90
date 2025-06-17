!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_salsa_dynamics.f90
!!
!! \brief
!! mo_ham_salsa_dynamics calculates coagulation and condensation for SALSA
!!
!! \author Harri Kokkola (FMI)
!!
!! \responsible_coder
!! Harri Kokkola, harri.kokkola@fmi.fi
!!
!! \revision_history
!!   -# Hannele Korhonen (FMI) - original code 2005
!!   -# Harri Kokkola (FMI) - 2006-2014
!!   -# Tommi Bergman (FMI) 2012
!!   -# Matti Niskanen(FMI) 2012
!!   -# Anton Laakso  (FMI) 2013
!!   -# Thomas Kuehn  (UEF) 2015
!! 
!! \limitations
!! None
!!
!! \details
!! None
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

MODULE mo_ham_salsa_dynamics


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! coagulation: Calculates particle loss and change in size distribution
!! due to (Brownian) coagulation. Semi-implicit, non-iterative method:
!! Volume concentrations of the smaller colliding particles added to the
!! bin of the larger colliding particles. Start from first bin and use the 
!! updated number and volume for calculation of following bins.
!! 
!! @author see module info 
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! salsa
!!
!! @par Responsible coder
!! harri.kokkola@fmi.fi
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE coagulation(kproma, kbdim,  klev,        &
                         pnaero, pvols,  pdwet,       &
                         pcore,  ptstep, ptemp, ppres)

    USE mo_ham_salsactl, ONLY:       &
         in1a,                       & ! size bin indices
         in2a, fn2a,                 &
         in2b, fn2b,                 &
         dpmid                         ! mid dry diameter of size bins [m]
    !USE mo_aero_mem_salsa, ONLY :    &
    !     d_cond_so4,                 & ! diagnostic variable for condensated mass of so4
    !     d_nuc_so4                     !diagnostic variable for nucleated mass of so4
    USE mo_ham_salsa_init, only: coagc

    USE mo_kind, ONLY : dp

    USE mo_math_constants, ONLY:    &
         pi_6

    USE mo_physical_constants, ONLY: &
         avo

    ! >> thk: VBS
    USE mo_ham, ONLY:&         
         sizeclass,                  & ! list of bin properties
         subm_naerospec_nowat,       &
         nsoa                          ! switch to turn on VBS

    USe mo_ham_vbsctl, ONLY:    &
         vbs_ngroup,            & ! number of VBS groups
         vbs_set,               & ! list of VBS grouops
         laqsoa,                & ! switch to turn on wet SOA
         aqsoa_ngroup,          & ! number of wet SOA groups 
         aqsoa_set                ! list of wet SOA grouops
    ! << thk

    IMPLICIT NONE


    !-- Input and output variables -------------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid kproma 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical klev 

    REAL(dp), INTENT(INOUT) ::      &     
         pnaero(kbdim,klev,fn2b),   & ! particle concentration [#/m3]
         ! >> thk: VBS
         pvols(kbdim,klev,fn2b,subm_naerospec_nowat)! total volume concentrations of each
                                      !   chem. compound in a size bin [fxm]
         ! << thk
         
    REAL(dp), INTENT(IN) ::         &
         pdwet(kbdim,klev,fn2b),    & ! particle wet radius [m]
         pcore(kbdim,klev,fn2b),    & ! particle dry volume [fxm]
         ptstep,                    & ! time step [s]
         ptemp(kbdim,klev),         &
         ppres(kbdim,klev)
    !-- Local variables ------------------------
    INTEGER ::                      &
         ii, jj, kk, ll, mm, nn,    & ! loop indices 
         index_2a, index_2b           ! corresponding bin in regime 2a/2b

    REAL(dp) ::                     &
         zntemp(kbdim,klev),        & ! variable for saving pnaero(fn2b) temporarily
         zcc(fn2b,fn2b), & ! updated coagulation coefficients [m3/s]
         zminusterm,                & ! coagulation loss in a bin [1/s] 
         ! >> thk: VBS
         zplusterm(subm_naerospec_nowat)         ! coagulation gain in a bin [fxm/s]
                                      ! (for each chemical compound)
         ! << thk

    REAL(dp) :: &
         zmpart(fn2b)   ! approximate mass of particles [kg]
    REAL(dp) :: &
         temppi,pressi,pdmm,pdnn

    ! >> thk: VBS
    ! for easier code readability
    INTEGER :: salsa_oc = 2
    ! << thk

    !-----------------------------------------------------------------------------
    !-- 1) Coagulation to coarse mode calculated in a simplified way: ------------
    !      CoagSink ~ Dp in continuum regime, thus we calculate
    !      'effective' number concentration of coarse particles

    zntemp = pnaero(:,:,fn2b)

    !-- 2) Updating coagulation coefficients -------------------------------------

    !-- particle mass; density of 1500 kg/m3 assumed [kg] 
    zmpart = pi_6*dpmid(1:fn2b)**3*1500.


     DO jj = 1,klev      ! vertical grid
        DO ii = 1,kproma ! horizontal kproma in the slab
           temppi=ptemp(ii,jj)
           pressi=ppres(ii,jj)

        
           DO mm = 1,fn2b         ! smaller colliding particle
              DO nn = mm,fn2b            ! larger colliding particle 
              pdmm=pdwet(ii,jj,mm)
              pdnn=pdwet(ii,jj,nn)
              zcc(mm,nn) = coagc(dpmid(mm),dpmid(nn),zmpart(mm),zmpart(nn),temppi,pressi)&
                     *dpmid(mm)*pdnn/(dpmid(nn)*pdmm)
              zcc(nn,mm)=zcc(mm,nn)
             END DO
          END DO
    
    ! Original version:
    ! loops over
    !DO nn = 1,fn2b            ! larger colliding particle
    !   DO mm = 1,fn2b         ! smaller colliding particle
    !      DO jj = 1,klev      ! vertical grid
    !         DO ii = 1,kproma ! horizontal kproma in the slab

               !zcc(ii,jj,mm,nn) = coagtable(ii,jj,mm,nn)*dpmid(mm)*pdwet(ii,jj,nn)/(dpmid(nn)*pdwet(ii,jj,mm))
    !            zcc(ii,jj,mm,nn) = coagc(dpmid(mm),dpmid(nn),zmpart(mm),zmpart(nn),ptemp(ii,jj),ppres(ii,jj))*dpmid(mm)*pdwet(ii,jj,nn) &
    !                 /(dpmid(nn)*pdwet(ii,jj,mm))

    !         END DO
    !      END DO
    !   END DO
    !END DO

    !-- 3) New particle and volume concentrations after coagulation -------------

    ! loops over
    DO kk = in1a,fn2b      ! bins that we want to update
       !DO jj = 1,klev      ! vertical grid
          !DO ii = 1,kproma ! horizontal kproma in the slab

             !-- 3.1) Bin in regime 1 -----------------------------------

             IF (kk < in2a) THEN

                !-- Particles lost by coagulation with larger particles
                !zminusterm = sum(zcc(ii,jj,kk,kk+1:fn2b)*pnaero(ii,jj,kk+1:fn2b))
                zminusterm = sum(zcc(kk,kk+1:fn2b)*pnaero(ii,jj,kk+1:fn2b))

                !-- Particle volume gained by coagulation with smaller particles

                zplusterm = 0._dp

                IF(kk > in1a) THEN

                   DO ll = in1a, kk-1
                      
                      ! >> thk: VBS
                      zplusterm(:) = zplusterm(:) + zcc(ll,kk)*pvols(ii,jj,ll,:)                      
                      ! << thk
                   END DO

                END IF

                !-- Volume and number concentrations after coagulation update [fxm]
                !pvols(ii,jj,kk,1:2) = pvols(ii,jj,kk,1:2)/(1._dp + ptstep*(zminusterm - 1./pcore(ii,jj,kk)*zplusterm(1:2)))
                ! >> thk: VBS
                pvols(ii,jj,kk,:) = (pvols(ii,jj,kk,:)+ptstep*zplusterm(:)*pnaero(ii,jj,kk))/ &
                                    (1._dp + ptstep*zminusterm) 
                ! << thk
                !pnaero(ii,jj,kk)    = pnaero(ii,jj,kk)/(1. + ptstep*zminusterm  + &
                !                      0.5_dp*ptstep*zcc(ii,jj,kk,kk)*pnaero(ii,jj,kk))
                pnaero(ii,jj,kk)    = pnaero(ii,jj,kk)/(1. + ptstep*zminusterm  + &
                                      0.5_dp*ptstep*zcc(kk,kk)*pnaero(ii,jj,kk))
             ELSE 

                !-- 3.2) Bin in regime 2a -----------------------------------

                IF (kk < in2b) THEN

                   zplusterm = 0._dp

                   !-- Find corresponding size bin in subregime 2b
                   index_2b = kk - in2a + in2b

                   !-- Particles lost to larger particles in regimes 2a, 2b and 3
                   zminusterm = 0._dp

                   DO ll = kk+1, fn2a

                      !zminusterm = zminusterm + zcc(ii,jj,kk,ll)*pnaero(ii,jj,ll) ! 2a
                      zminusterm = zminusterm + zcc(kk,ll)*pnaero(ii,jj,ll) ! 2a
                   END DO

                   DO ll = index_2b+1, fn2b

                      !zminusterm = zminusterm + zcc(ii,jj,kk,ll)*pnaero(ii,jj,ll) ! 2b 
                      zminusterm = zminusterm + zcc(kk,ll)*pnaero(ii,jj,ll) ! 2b + 3

                   END DO

                   !-- Particle volume gained from smaller particles in regimes 1, 2a and 2b

                   ! sulphate

                   DO ll = in1a, kk-1

                      !zplusterm(1:2) = zplusterm(1:2) + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,1:2) ! 1 + 2a
                      ! >> thk: VBS
                      zplusterm(:) = zplusterm(:) + zcc(ll,kk)*pvols(ii,jj,ll,:) ! 1 + 2a
                      ! << thk

                   END DO
                   
                   DO ll = in2a, kk-1

                      !zplusterm(3:5) = zplusterm(3:5) + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,3:5)
                      ! >> thk: VBS
                      zplusterm(:) = zplusterm(:) + zcc(ll,kk)*pvols(ii,jj,ll,:)
                      ! << thk
                   
                   END DO
                      
                   DO ll = in2b, index_2b

                      !zplusterm = zplusterm + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,:) ! 2b
                      zplusterm = zplusterm + zcc(ll,kk)*pvols(ii,jj,ll,:) ! 2b

                   END DO
                   
                   !-- 3.2) Bin in regime 2b -----------------------------------

                ELSE

                   !-- Find corresponding size bin in subregime 2a
                   index_2a = kk - in2b + in2a

                   !-- Particles lost to larger particles in regimes 2b and 3, 
                   !  as well as to 'same sized' and larger particles in regime 2a             

                   zminusterm = 0._dp
                   zplusterm  = 0._dp

                   DO ll = kk+1, fn2b

                      !zminusterm = zminusterm + zcc(ii,jj,kk,ll)*pnaero(ii,jj,ll) ! 2b + 3
                      zminusterm = zminusterm + zcc(kk,ll)*pnaero(ii,jj,ll) ! 2b + 3

                   END DO

                   DO ll = index_2a, fn2a

                      !zminusterm = zminusterm + zcc(ii,jj,kk,ll)*pnaero(ii,jj,ll)
                      zminusterm = zminusterm + zcc(kk,ll)*pnaero(ii,jj,ll)

                   END DO

                   !-- Particle volume gained from smaller particles in regimes 1, 2a and 2b
                   !  (except not corresponding bin in regime 2a even if it is slightly smaller)

                   ! sulphate

                   DO ll = in1a, index_2a-1

                      !zplusterm(1:2) = zplusterm(1:2) + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,1:2)
                      ! >> thk: VBS
                      zplusterm(:) = zplusterm(:) + zcc(ll,kk)*pvols(ii,jj,ll,:)
                      ! << thk
                   END DO


                   DO ll = in2a, index_2a-1

                      !zplusterm(3:5) = zplusterm(3:5) + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,3:5)
                      ! >> thk: VBS
                      zplusterm(:) = zplusterm(:) + zcc(ll,kk)*pvols(ii,jj,ll,:)
                      ! << thk

                   END DO

                   DO ll = in2b, kk-1

                      !zplusterm = zplusterm + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,1:5)
                      ! >> thk: VBS
                      zplusterm = zplusterm + zcc(ll,kk)*pvols(ii,jj,ll,:)
                      ! << thk
                   END DO

                END IF

                
                !-- Volume and number concentrations after coagulation update [fxm]
                !pvols(ii,jj,kk,:) = pvols(ii,jj,kk,:)/ &
                !                    (1._dp + ptstep*(zminusterm - 1./pcore(ii,jj,kk)*zplusterm))
                pvols(ii,jj,kk,:) = (pvols(ii,jj,kk,:)+ptstep*zplusterm*pnaero(ii,jj,kk))/ &
                                    (1._dp + ptstep*zminusterm)

                !pnaero(ii,jj,kk) = pnaero(ii,jj,kk)/(1. + ptstep*zminusterm  + &
                !                   0.5_dp*ptstep*zcc(ii,jj,kk,kk)*pnaero(ii,jj,kk))
                pnaero(ii,jj,kk) = pnaero(ii,jj,kk)/(1. + ptstep*zminusterm  + &
                                   0.5_dp*ptstep*zcc(kk,kk)*pnaero(ii,jj,kk))

             END IF
          END DO
       END DO
    END DO

    ! fxm: here we approximate that the sea salt regime 2b particles have
    ! gained by coagulation can be treated as sulphate
    pvols(1:kproma,:,in2b:fn2b,1) = pvols(1:kproma,:,in2b:fn2b,1) + pvols(1:kproma,:,in2b:fn2b,4)
    pvols(1:kproma,:,in2b:fn2b,4) = 0._dp

    ! >> thk: VBS
    IF (nsoa == 2) THEN
       ! Regime 2b does not contain semi-volatiles, which leaves us two options
       ! 1) convert everything into OC
       ! 2) move the accumulated semi-volatiles back to the gas phase
       ! we go for option one first, as it is much faster to implement
       DO kk = 1, fn2b
          IF (.NOT. (sizeclass(kk)%lsoainclass)) THEN
             DO ll = 1,vbs_ngroup
! >> thk: bugfix
!                IF (ll == salsa_oc) CYCLE
                IF (vbs_set(ll)%id_vols == salsa_oc) CYCLE
                pvols(1:kproma,:,kk,salsa_oc) = pvols(1:kproma,:,kk,salsa_oc)&
                     + pvols(1:kproma,:,kk,vbs_set(ll)%id_vols)
                pvols(1:kproma,:,kk,vbs_set(ll)%id_vols)=0._dp
! << thk
             END DO ! ll
          END IF ! .NOT. lsoainclass
       END DO ! kk
    END IF ! nsoa == 2
    IF (laqsoa) THEN 
       ! Regime 2b does not contain semi-volatiles, which leaves us two options
       ! 1) convert everythin into OC
       ! 2) move the accumulated semi-volatiles back to the gas phase
       ! we go for option one first, as it is much faster to implement
       DO kk = 1, fn2b
          IF (.NOT. (sizeclass(kk)%lsoainclass)) THEN
                DO ll = 1,aqsoa_ngroup
                   pvols(1:kproma,:,kk,salsa_oc) = pvols(1:kproma,:,kk,salsa_oc)&
                        + pvols(1:kproma,:,kk,aqsoa_set(ll)%id_aqsoa)
                   pvols(1:kproma,:,kk,aqsoa_set(ll)%id_aqsoa)=0._dp
                END DO ! ll             
          END IF ! .NOT. lsoainclass
       END DO ! kk
    END IF ! laqsoa


    pvols(1:kproma,:,:,:) = MAX(pvols(1:kproma,:,:,:), 0._dp)

  END SUBROUTINE coagulation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! coagulation: Calculates the increase in particle volume and 
!!  decrease in gas phase concentrations due to condensation 
!!  of sulphuric acid and two organic compounds (non-volatile
!!  and semivolatile)
!! 
!! @author see module info 
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! salsa
!!
!! @par Responsible coder
!! harri.kokkola@fmi.fi
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE condensation(kproma, kbdim,  klev,   krow,              &
                          pnaero, pvols,  pdwet,   pc_gas,           &
                          prh,    ptemp,  ppres,  ptstep, ppbl       &
                          )

    USE mo_species,                 ONLY: speclist

    USE mo_ham_species,             ONLY: id_so4, &
                                          id_ss,  &
                                          id_oc,  &
                                          id_wat
    USE mo_ham_salsa_nucleation

    USE mo_physical_constants, ONLY:&
         avo,                       &
         argas,                     &
         ak,                        &
         p0sl_bg

    USE mo_math_constants, ONLY:    &
         pi,                        &
         pi_6

    USE mo_ham_salsactl,   ONLY :   &
         in1a, in2a,                & ! size bin indices
         fn2b,                      & 
         nbin,                      & ! number of size bins in each regime
         d_sa,                      & ! diameter of H2SO4 molecule [m]
         
         epsoc,                     & ! soluble fraction of organics (scaled to sulphate)
         massacc,                   & ! mass accomodation coefficients in each bin
         dpmid,                     & ! mid dry diameter of each bin [m]
         n3,                        & ! number of molecules in one 3 nm particle [1]
         nsnucl                       ! nucleation

    USE mo_kind, ONLY : dp
    !USE mo_aero_mem_salsa, only: d_nuc_so4,d_cond_so4
    USE mo_time_control,       ONLY: delta_time,time_step_len

    ! >> thk: volatility basis set (VBS)
    USE mo_ham, ONLY:      &
         subm_ngasspec,             & ! number of gas species in HAM
         subm_naerospec_nowat,      & ! number of non-water aerosol species in HAM
         nsoa                         ! switch to turn VBS on and off

    USE mo_ham_vbsctl, ONLY:        &
         laqsoa

    USE mo_ham_vbs_partition, ONLY: &
         vbs_condensation_salsa, &  ! compute condensation of organics for SALSA
         aqsoa_condensation_salsa

    USE mo_ham_subm_species,   ONLY: isubm_ocnv, isubm_so4g
    ! << thk

   IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid kproma 
         kbdim,                     & ! dimension for arrays 
         klev,                      & ! number of vertical klev 
         krow

    REAL(dp), INTENT(IN) ::        &  
         pdwet(kbdim,klev,fn2b),   & ! wet diameter of particles in each bin [m]
         prh(kbdim,klev),          & ! relative humidity [0-1]
         ptemp(kbdim,klev),        & ! ambient temperature [K]
         ppres(kbdim,klev),        & ! ambient pressure [Pa]
         ptstep                     ! timestep [s]


    INTEGER :: ppbl(kbdim)           ! Planetary boundary layer top level

    ! >> thk: VBS
    REAL(dp), INTENT(INOUT) ::     &
         pnaero(kbdim,klev,fn2b),  & ! number concentration of particles in each bin [#/m3]
         pvols(kbdim,klev,fn2b,subm_naerospec_nowat),& ! volume concentration of each chem. compound in regimes 1&2 [fxm]
         pc_gas(kbdim,klev,subm_ngasspec)    ! gas phase tracers
    ! << thk


    !-- Local variables ----------------------
    INTEGER :: ii, jj                ! loop indices

    REAL(dp) ::                    &
         zvisc,                    & ! viscosity of air [kg/(m s)]
         zdfvap,                   & ! air diffusion coefficient [m2/s]
         zmfp,                     & ! mean free path of condensing vapour [m]
         zcs_tot,                  & ! total condensation sink [1/s]
         zcs_su,                   & ! condensation sink for sulfate [1/s]
         zcs_ocnv,                 & ! condensation sink for nonvolatile organics [1/s]
                                     ! vapour concentration after time step [#/m3]
         zcvap_new1,               & ! sulphuric acid
         zcvap_new2,               & ! nonvolatile organics
         zcvap_new3,               & ! semivolatile organics
                                     ! change in vapour concentration [#/m3]
         zdvap1,                   & ! sulphuric acid
         zdvap2,                   & ! nonvolatile organics
         zdvap3,                   & ! semivolatile organics
         
         zdfpart(in1a+1),          & ! particle diffusion coefficient
         zknud(fn2b),              & ! particle Knudsen number
         zbeta(fn2b),              & ! transitional correction factor
         zcolrate(fn2b),           & ! collision rate of molecules to particles [1/s]
         zcolrate_ocnv(fn2b),      & ! collision rate of organic molecules to particles [1/s]
         zdvolsa(fn2b),            & ! change of sulphate volume in each bin [fxm]
         zdvoloc(fn2b),            & !    - " - organics 
         zj3n3(kbdim,klev,2),      & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
         zn_vs_c,                  & ! ratio of nucleation of all mass transfer in the smallest bin
         zxsa(kbdim,klev),         & ! ratio of sulphuric acid and organic vapor in 3nm particles 
         zxocnv(kbdim,klev),       & !
         zmvsu,                    & ! molar volume of sulfate
         zmvoc                       ! molar volume of organics

    real(dp)::ztmst,zqtmst

    !-- Initializations
    zmvsu = (speclist(id_so4)%moleweight/1000.)/avo/speclist(id_so4)%density
    zmvoc = (speclist(id_oc)%moleweight/1000.)/avo/speclist(id_oc)%density
    ztmst = time_step_len
    zqtmst = 1.0_dp/ztmst
    zxsa = 0._dp
    zj3n3 = 0._dp
    
    !------------------------------------------------------------------------------
    
    IF (nsnucl > 0) CALL nucleation(kproma, kbdim,  klev,   krow,                     &
                                 pnaero, pdwet,                                    &
                                 ptemp,  prh,    ppres,                            &
                                 ptstep, pc_gas,                                   &
                                 zj3n3, zxsa, zxocnv, ppbl                         &
                                 )
    zdvolsa=0._dp
    zn_vs_c=0._dp
    DO jj = 1,klev
       DO ii = 1,kproma

          zdvoloc = 0._dp

          !-- 1) Properties of air and condensing gases --------------------

          zvisc  = (7.44523e-3_dp*ptemp(ii,jj)**1.5_dp)/(5093._dp*(ptemp(ii,jj)+110.4_dp))! viscosity of air [kg/(m s)] 
          zdfvap = 5.1111e-10_dp*ptemp(ii,jj)**1.75_dp*p0sl_bg/ppres(ii,jj)                ! diffusion coefficient [m2/s]
          zmfp   = 3._dp*zdfvap*sqrt(pi*(speclist(id_so4)%moleweight/1000.)/ &
               (8._dp*argas*ptemp(ii,jj)))                      ! mean free path [m]


          !-- 2) Transition regime correction factor for particles ---------
          !  
          !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
          !  Topics in current aerosol research, Pergamon.  
          !
          !  Size of condensing molecule considered only for 
          !  nucleation mode (3 - 20 nm) 
          !

          !-- particle Knudsen number
          zknud(in1a:in1a+1) = 2.*zmfp/(pdwet(ii,jj,in1a:in1a+1)+d_sa) 
          zknud(in1a+2:fn2b) = 2.*zmfp/pdwet(ii,jj,in1a+2:fn2b)

          !-- transitional correction factor
          zbeta = (zknud + 1.)/(0.377_dp*zknud+1._dp+4._dp/ &
                  (3._dp*massacc)*(zknud+zknud**2))  

          !-- 3) Collision rate of molecules to particles -------------------
          !
          !  Particle diffusion coefficient considered only for 
          !  nucleation mode (3 - 20 nm) 
          !

          !-- particle diffusion coefficient [m2/s]
          zdfpart = ak*ptemp(ii,jj)*zbeta(in1a:in1a+1)/ &    
                    (3._dp*pi*zvisc*pdwet(ii,jj,in1a:in1a+1))  

          !-- collision rate [1/s]
          zcolrate(in1a:in1a+1) = 2._dp*pi*(pdwet(ii,jj,in1a:in1a+1)+d_sa)* & 
                                  (zdfvap+zdfpart)*zbeta(in1a:in1a+1)* &
                                  pnaero(ii,jj,in1a:in1a+1)
         
          zcolrate(in1a+2:fn2b) = 2._dp*pi*pdwet(ii,jj,in1a+2:fn2b)*zdfvap* &
                                  zbeta(in1a+2:fn2b)*pnaero(ii,jj,in1a+2:fn2b)

          !-- 4) Condensation sink [1/s] -------------------------------------

          zcs_tot = sum(zcolrate)             ! total sink    
          
          !-- 5) Changes in gas-phase concentrations and particle volume -----
          !
          !  Sulphuric acid and non-volatile organic compound
          !  condense onto all sized particles. Semivolatile
          !  organic compound condenses only onto particles
          !  in regimes 2 and 3.
          !
          !  Regime 3 particles act only as a sink for vapours
          !  and their size and composition are not updated.
          !  (except for soluble fraction of subregime 3c) 
          !
          

          
          !--- 5.1) Organic vapours ------------------------


          ! >> thk: in the case VBS is on, call own condensation routine for
          !         organics
          IF (nsoa == 2) THEN
             
             call vbs_condensation_salsa(&
                      ptemp(ii,jj),                           & ! ambient conditions
                      ptstep,                                 & ! model parameters
                      pc_gas(ii,jj,:),                        & ! gas phase concentrations
                      pnaero(ii,jj,:), pvols(ii,jj,:,:),      & ! aerosol properties
                      pdwet(ii,jj,:), zcolrate                & ! quantities from condensation code
                      )

!             !-- Change of number concentration in the smallest bin caused by nucleation
!             !   Jacobson (2005), equation (16.75)
!             ! If zxocnv = 0, then the chosen nucleation mechanism does not take into account
!             ! the nonvolatile organic vapors and thus the pnaero does not have to be updated.
!             IF (zxocnv(ii,jj) > 0._dp) THEN 
!                pnaero(ii,jj,in1a) = pnaero(ii,jj,in1a) + &
!                     zn_vs_c * zdvoloc(in1a)/zmvoc/(n3*zxocnv(ii,jj))
!             END IF

             IF (laqsoa) THEN 
                
                call aqsoa_condensation_salsa(&
                     ptemp(ii,jj),                                  & ! ambient conditions
                     pc_gas(ii,jj,:),                                 & ! gas phase concentrations
                     pnaero(ii,jj,:),                                 & ! 
                     pvols(ii,jj,:,:),                                  & ! aerosol properties
                     pdwet(ii,jj,:)                                    &
                     ) 
             ENDIF
          ELSE
          
          !---- 5.1.1) Non-volatile organic compound: condenses onto all bins 

             ! >> thk: VBS
             IF(pc_gas(ii,jj,isubm_ocnv) > 1.e-10_dp .and. zcs_tot > 1.e-30_dp) THEN
             !<< thk

                zn_vs_c = 0._dp

                ! >> thk: VBS
                IF(zj3n3(ii,jj,2) > 1._dp) zn_vs_c = (zj3n3(ii,jj,2))/(zj3n3(ii,jj,2) + &
                     pc_gas(ii,jj,isubm_ocnv) * zcolrate(in1a))
                ! << thk

                !   collision rate in the smallest bin, including nucleation and condensation
                !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
                !   equation (16.73)
                zcolrate_ocnv = zcolrate
                ! >> thk: VBS
                zcolrate_ocnv(in1a) = zcolrate_ocnv(in1a) + zj3n3(ii,jj,2)/pc_gas(ii,jj,isubm_ocnv)

                zcs_ocnv = zcs_tot + zj3n3(ii,jj,2)/pc_gas(ii,jj,isubm_ocnv) ! total sink for organic vapor
                ! << thk

                ! >> thk: VBS
                zcvap_new2 = pc_gas(ii,jj,isubm_ocnv)/(1.+ptstep*zcs_ocnv) ! new gas phase concentration [#/m3]
                zdvap2 = pc_gas(ii,jj,isubm_ocnv) - zcvap_new2             ! change in gas concentration [#/m3]
                pc_gas(ii,jj,isubm_ocnv) = zcvap_new2                      ! updating vapour concentration [#/m3]
                ! << thk

                zdvoloc = zcolrate_ocnv(in1a:fn2b)/zcs_ocnv*zmvoc*zdvap2 ! volume change of particles 
                !  [m3(OC)/m3(air)]

                pvols(ii,jj,in1a:fn2b,2) = pvols(ii,jj,in1a:fn2b,2) + & !-- change of volume
                     zdvoloc                      !   due to condensation in 1a-2b

                !-- Change of number concentration in the smallest bin caused by nucleation
                !   Jacobson (2005), equation (16.75)
                ! If zxocnv = 0, then the chosen nucleation mechanism does not take into account
                ! the nonvolatile organic vapors and thus the pnaero does not have to be updated.
                IF (zxocnv(ii,jj) > 0._dp) THEN 
                   pnaero(ii,jj,in1a) = pnaero(ii,jj,in1a) + &
                        zn_vs_c * zdvoloc(in1a)/zmvoc/(n3*zxocnv(ii,jj))
                END IF

             END IF

             ! >> thk: VBS
             ! old semivolatile organics code removed 
             ! << thk
          END IF


          ! << thk

          !--- 5.2) Sulphuric acid -------------------------
          !
          ! nucleh2so4=nucleh2so4+pcsa(ii,jj)

          ! >> thk: VBS
          IF(pc_gas(ii,jj,isubm_so4g) > 1.e-10_dp .and. zcs_tot > 1.e-30_dp) THEN
          ! << thk
             !-- Ratio of mass transfer betwwn nucleation and condensation

             zn_vs_c = 0._dp

             ! >> thk: VBS
             IF(zj3n3(ii,jj,1) > 1._dp) zn_vs_c = (zj3n3(ii,jj,1)) / &
                                              (zj3n3(ii,jj,1) +  &
                                              pc_gas(ii,jj,isubm_so4g) * zcolrate(in1a))
             ! << thk
             !   collision rate in the smallest bin, including nucleation and condensation
             !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
             !   equation (16.73)
             ! >> thk: VBS
             zcolrate(in1a) = zcolrate(in1a) + zj3n3(ii,jj,1) / pc_gas(ii,jj,isubm_so4g)
             zcs_su = zcs_tot + zj3n3(ii,jj,1) / pc_gas(ii,jj,isubm_so4g) ! total sink for sulfate

             zcvap_new1 = pc_gas(ii,jj,isubm_so4g) /(1.+ptstep*zcs_su)! new gas phase concentration [#/m3]
             zdvap1 = pc_gas(ii,jj,isubm_so4g) - zcvap_new1           ! change in gas concentration [#/m3]
             pc_gas(ii,jj,isubm_so4g) = zcvap_new1                    ! updating vapour concentration [#/m3]
             ! << thk

             zdvolsa = zcolrate(in1a:fn2b)/zcs_su*zmvsu*zdvap1     ! volume change of particles
             ! [m3(SO4)/m3(air)] by condensation

             !-- Change of volume concentration of sulphate in aerosol [fxm]
             pvols(ii,jj,in1a:fn2b,1) = pvols(ii,jj,in1a:fn2b,1) + zdvolsa

             !-- Change of number concentration in the smallest bin caused by nucleation
             !   Jacobson (2005), equation (16.75)
             IF (zxsa(ii,jj) > 0._dp) THEN
                pnaero(ii,jj,in1a) = pnaero(ii,jj,in1a) +          &
                     zn_vs_c * zdvolsa(in1a)/zmvsu/(n3*zxsa(ii,jj))
             END IF

             ! Diagnostic for nucleation of so4

             !d_nuc_so4(ii,krow)=d_nuc_so4(ii,krow)+zn_vs_c*zdvolsa(in1a)*pdpg(ii,jj)*speclist(id_so4)%density*zqtmst*delta_time/3.0_dp
             ! Diagnostic for condensation of so4
             !d_cond_so4(ii,krow)=d_cond_so4(ii,krow)+(sum(zdvolsa)-zn_vs_c*zdvolsa(in1a))*zqtmst*delta_time*speclist(id_so4)%density*pdpg(ii,jj)/3.0_dp
             
          END IF

       END DO
    END DO


  END SUBROUTINE condensation


END MODULE mo_ham_salsa_dynamics
