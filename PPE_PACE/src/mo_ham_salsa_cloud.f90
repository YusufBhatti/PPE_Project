!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_salsa.f90
!!
!! \brief
!! Implementation of Abdul-Razzak & Ghan activation scheme for SALSA.
!! 
!!
!! \author Harri Kokkola (FMI)
!!
!! \responsible_coder
!! Thomas Kuehn -- thomas.h.kuhn@uef.fi
!!
!! \revision_history
!! T. Anttila (FMI)     2007
!! H. Kokkola (FMI)     2007
!! A.-I. Partanen (FMI) 2007
!! T. Kuehn (UEF)       2015
!!
!! \limitations
!! None
!!
!! \details
!! Purpose: Calculates the number of activated cloud 
!! droplets according to parameterizations by:
!!
!! Abdul-Razzak et al: "A parameterization of aerosol activation - 
!!                      3. Sectional representation"
!!                      J. Geophys. Res. 107, 10.1029/2001JD000483, 2002. 
!!                      [Part 3]
!!
!! Abdul Razzak et al: "A parameterization of aerosol activation - 
!!                      1. Single aerosol type"
!!                      J. Geophys. Res. 103, 6123-6130, 1998. 
!!                      [Part 1]
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

MODULE mo_ham_salsa_cloud

  PRIVATE

  PUBLIC :: salsa_abdul_razzak_ghan

CONTAINS 

  SUBROUTINE cloud_activation(&
       kproma,   kbdim, klev,  krow,  ktdia,    &
       pnaero,   pvols, ptm1,  papm1, pqm1,     &    
       pcdncact, pw,     pwpdf, psc,   psmax    &
       )

    USE mo_kind,        ONLY : dp

    USE mo_ham, ONLY:                           &
         aerocomp, sizeclass, nclass,           &
         subm_naerospec_nowat,                  & 
         subm_aerospec

    USE mo_species, ONLY :                      &
         speclist

    USE mo_ham_species, ONLY:                   &
         id_so4,                                &
         id_ss,                                 &
         id_oc,                                 &
         id_wat

    USE mo_ham_salsactl, ONLY :                 &
         slim,                                  &
         surfw0,                                & ! surface tension 
                                                  ! of water [J/m2]
         nbin,                                  & ! number of size bins 
                                                  ! in subranges
         nlim,                                  & ! lowest possible particle conc. in a bin [#/m3]
         in1a,in2a,in2b,fn1a,fn2a, fn2b,        & ! size regime bin indices
         vlolim, vhilim,                        & ! low and high volume
                                                  ! limits for bins
         vratiohi, vratiolo                       ! ratio of volume in the bin (high and low) 
                                                  ! boundary and the middle of the bin

    USE mo_activ,        ONLY: nw

    USE mo_echam_cloud_params, ONLY: cthomi ! minimum T for mixed clouds

    USE mo_physical_constants, ONLY: &
         argas,                      &
         grav,                       &
         tmelt,                      &
         p0sl_bg,                    &
         cpd,                        &
         mair

    USE mo_math_constants, ONLY:     &
         pi

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::                      &
             kproma,                            & ! number of horiz. grid points 
             kbdim,                             & ! dimension for arrays 
             klev,                              & ! number of vertical levels 
             krow,                              & ! block number
             ktdia                                ! highest vertical level for "diagnostics"

    REAL(dp), INTENT(IN)::                      &
             pnaero(kbdim,klev,fn2b),           & ! number concentration  [#/m3]
             papm1(kbdim,klev),                 & ! atmospheric pressure at grid point [Pa]
             pqm1(kbdim,klev),                  & ! specific humidity
             ptm1(kbdim,klev),                  & ! temperature [K]
             pvols(kbdim,klev,fn2b,subm_naerospec_nowat), & ! total volume concentrations of each
                                                  ! chem. compound in a size bin
             pw(kbdim,klev, nw),                & ! mean or bins of updraft velocity (>0.0) [m/s]
             pwpdf(kbdim,klev,nw)!,              & ! PDF of updraft velocity [s m-1]
  
    REAL(dp), INTENT(OUT) ::                    &
             pcdncact(kbdim,klev),              & ! number of cloud droplets
             psc(kbdim,klev,nclass),            & ! critical supersaturation [% 0-1]
             psmax(kbdim,klev,nw)                 ! maximum supersaturation [% 0-1]


    !-- local variables --------------
    REAL(dp) ::                         &
             zeps,                      & ! small number
             sil,                       & ! critical supersaturation 
                                          !     at the upper bound of the bin 
             siu,                       & !  "  at the lower bound of the bin
             zaa,                        & ! curvature (Kelvin) effect [m]
             zbb,                        & ! solute (Raoult) effect [m3/mol]
             zkoehl,                     & ! argument of exp for critical supersaturation [mol]
             zns(fn2b),                  & ! number of moles of solute
             znshi,                      & !             " at the upper bound of the bin
             znslo,                      & !             " at the lower bound of the bin
             zseff,                     & ! effective supersaturation
             zka,                       & ! thermal conductivity
             zdv,                       & ! diffusion coefficient
             zgc,                        & ! growth coefficient
             zalpha,                     & ! see Abdul-Razzak and Ghan, part 3
             zgamma,                     & ! see Abdul-Razzak and Ghan, part 3
             zevap,                      & ! latent heat of evaporation
             zps,                       & ! saturation vapor pressure of water [Pa]
             za1,                       & ! helper variable in calculation of zps [unitless]
             zkhi,                       & ! see Abdul-Razzak and Ghan, part 3
             ztheta,                     & ! see Abdul-Razzak and Ghan, part 3
             zfrac(fn2b),                & ! fraction of activated droplets in a bin
             zntot,                      & ! total number conc of particles [#/m3]
             zsum1,                      & ! sum in denominator of Eq (8) in A-R&G p. 3
             zsum2,                      & ! sum in denominator of Eq (9) in A-R&g p.3
             zpdfsum                       ! for pdf normalization

    REAL(dp), DIMENSION(kbdim, klev, fn2b) :: &
         zvols_sol,   & ! volume mixing ratio of soluable substances
         zvols_insol, & ! volume mixing ratio of insoluable substances
         zvols_tot      ! total volume concentration


    INTEGER :: ii, jj, jclass, js, jw         ! loop indices

    zeps=EPSILON(1._dp)

    zbb = 6._dp*(speclist(id_wat)%moleweight/1000.)/ &
         (pi*speclist(id_wat)%density)    ! Raoult effect [m3/mol]
                                          ! NOTE!
                                          ! bb must be multiplied
                                          ! by the number of moles of
                                          ! solute

    pcdncact(1:kproma,:) = 0._dp
    psc(1:kproma,:,:)    = 0._dp
    psmax(1:kproma,:,:)  = 0._dp

    ! Keeping track of what is soluble and what is insoluble
    zvols_sol(1:kproma,:,:)  =0.0_dp
    zvols_insol(1:kproma,:,:)=0.0_dp
    zvols_tot(1:kproma,:,:)  =0.0_dp

    DO jj = ktdia, klev
       DO ii = 1, kproma
          DO jclass = in2a, fn2b
             zvols_sol(ii,jj,jclass)   = pvols(ii,jj,jclass,1) + pvols(ii,jj,jclass,2) + pvols(ii,jj,jclass,4)
             zvols_insol(ii,jj,jclass) = pvols(ii,jj,jclass,3) + pvols(ii,jj,jclass,5)
             zvols_tot(ii,jj,jclass)   = zvols_sol(ii,jj,jclass)+ zvols_insol(ii,jj,jclass)
          END DO
       END DO
    END DO

    DO jj = ktdia,klev    ! vertical grid
       DO ii = 1,kproma   ! horizontal grid

          zaa = 4._dp*(speclist(id_wat)%moleweight/1000.)*surfw0/ &
               (argas*speclist(id_wat)%density*ptm1(ii,jj)) ! Kelvin effect [m]
          zkoehl  = 4._dp*zaa**3/(27._dp*zbb)
          zns(:) = 0._dp

          IF( & ! are the minimal conditions fulfilled?
               nw > 1 .OR. pw(ii,jj,1) > zeps .AND. &
               pqm1(ii,jj)>zeps               .AND. &
               ptm1(ii,jj)>cthomi                   &
               ) THEN  

           
             zaa = 4._dp*(speclist(id_wat)%moleweight/1000.)*surfw0/ &
                  (argas*speclist(id_wat)%density*ptm1(ii,jj)) ! Kelvin effect [m]
             zkoehl  = 4._dp*zaa**3/(27._dp*zbb)

             zntot  = 0._dp
             zsum1  = 0._dp

             !-- subrange 1a -- neglected

             !-- subrange 2a
             DO jclass = in2a, fn2a
                IF ( pnaero(ii,jj,jclass) > nlim .AND. &
                     zvols_sol(ii,jj,jclass) > 1e-20_dp*zvols_insol(ii,jj,jclass) .AND. &
                     zvols_sol(ii,jj,jclass) > 1e-3_dp*vlolim(jclass)*pnaero(ii,jj,jclass)&
                   ) THEN

                   !-- number of moles of solute in one particle [mol]
                   zns(jclass) = (&
                        3._dp*pvols(ii,jj,jclass,1) * speclist(id_so4)%density/  &
                        (speclist(id_so4)%moleweight/1000.)                      &
                        + 1._dp*pvols(ii,jj,jclass,2) * speclist(id_oc)%density/ &
                        (speclist(id_oc)%moleweight/1000.)                       &
                        + 2._dp*pvols(ii,jj,jclass,4) * speclist(id_ss)%density/ &
                        (speclist(id_ss)%moleweight/1000.)                       &
                        )/pnaero(ii,jj,jclass)  

                   !-- critical supersaturation, Kohler equation
                   psc(ii,jj,jclass) = exp(sqrt(zkoehl/zns(jclass))) - 1._dp

                   !-- sums in equation (8), part 3
                   zntot = zntot + pnaero(ii,jj,jclass)
                   zsum1 = zsum1 + pnaero(ii,jj,jclass)/psc(ii,jj,jclass)**(2._dp/3._dp)
                   
                END IF

             END DO

             !-- subrange 2b
             DO jclass = in2b, fn2b
                
                IF ( pnaero(ii,jj,jclass) > nlim .AND. &
                     zvols_sol(ii,jj,jclass) > 1e-20_dp*zvols_insol(ii,jj,jclass) .AND. &
                     zvols_sol(ii,jj,jclass) > 1e-3_dp*vlolim(jclass)*pnaero(ii,jj,jclass)&
                     ) THEN

                   !-- number of moles of solute in one particle [mol]
                   zns(jclass) = (&
                          3._dp*pvols(ii,jj,jclass,1) * speclist(id_so4)%density/ &
                          (speclist(id_so4)%moleweight/1000.)                     &
                          + 1._dp*pvols(ii,jj,jclass,2) * speclist(id_oc)%density/&
                          (speclist(id_oc)%moleweight/1000.)                      &
                          )/pnaero(ii,jj,jclass)  
                
                   !-- critical supersaturation, Kohler equation
                   psc(ii,jj,jclass) = exp(sqrt(zkoehl/zns(jclass))) - 1._dp

                   !-- sums in equation (8), part 3
                   zntot = zntot + pnaero(ii,jj,jclass)
                   zsum1 = zsum1 + pnaero(ii,jj,jclass)/psc(ii,jj,jclass)**(2._dp/3._dp)
                   
                END IF

             END DO

             psmax(ii,jj,:) = 1._dp

             IF(zntot < nlim .OR. zsum1 < nlim) CYCLE

             !-- latent heat of evaporation [J/kg]
             zevap  = 2.501e6_dp-2370._dp*(ptm1(ii,jj)-273.15_dp)

             !-- saturation vapor pressure of water [Pa]
             !   Seinfeld & Pandis (1.10)
             za1    = 1._dp-(373.15_dp/ptm1(ii,jj))
             zps    = p0sl_bg*exp(13.3185_dp*za1-1.976_dp*za1**2-0.6445_dp*za1**3-0.1299_dp*za1**4)     

             !-- part 1, eq (11)
             zalpha = grav*(speclist(id_wat)%moleweight/1000.)*zevap/(cpd*argas*ptm1(ii,jj)**2)-                            &
                  grav*mair/(argas*ptm1(ii,jj))

             !-- part 1, eq (12)
             zgamma = argas*ptm1(ii,jj)/(zps*(speclist(id_wat)%moleweight/1000.)) &
                  + (speclist(id_wat)%moleweight/1000.)*zevap**2/(cpd*papm1(ii,jj)*mair*ptm1(ii,jj))

             !-- diffusivity [m2/s], Seinfeld and Pandis (15.65)
             !  Eq (17.61) in second edition
             zdv= 1.e-4_dp * (0.211_dp*p0sl_bg/papm1(ii,jj)) * ((ptm1(ii,jj)/tmelt)**1.94_dp)

             !-- thermal conductivity [J/(m s K)], Seinfeld and Pandis (15.75)
             ! Eq (17.71) in second edition
             zka= 1.e-3_dp * (4.39_dp + 0.071_dp * ptm1(ii,jj))

             !-- growth coefficient, part 1, eq (16)
             !-- (note: here uncorrected diffusivities and conductivities are used
             !    based on personal communication with H. Abdul-Razzak, 2007)
             zgc = 1._dp/(speclist(id_wat)%density*argas*ptm1(ii,jj)/                  &
                  (zps*zdv*(speclist(id_wat)%moleweight/1000.)) +                      &
                  zevap*speclist(id_wat)%density/(zka*ptm1(ii,jj)) *                   &
                  (zevap*(speclist(id_wat)%moleweight/1000.)/(ptm1(ii,jj)*argas)-1._dp))

             !-- effective critical supersaturation: part 3, eq (8)
             zseff = (zntot/zsum1)**(3._dp/2._dp)

             DO jw = 1,nw    ! updraft bins

                !-- part 3, equation (5)
                ztheta = ((zalpha*pw(ii,jj,jw)/zgc)**(3._dp/2._dp))/&
                     (2._dp*pi*speclist(id_wat)%density*zgamma*zntot)

                !-- part 3, equation (6)
                zkhi = (2._dp/3._dp)*zaa*SQRT(zalpha*pw(ii,jj,jw)/zgc)

                !-- maximum supersaturation of the air parcel: part 3, equation (9)
                zsum2 = 0.5_dp*(zkhi/ztheta)**(3._dp/2._dp)     &
                     + ((zseff**2)/(ztheta+3._dp*zkhi))**(3._dp/4._dp)

                IF (zsum2 > zeps) THEN
                   psmax(ii,jj,jw) = zseff / SQRT(zsum2)
                END IF

             END DO !jw

          ELSE
             psmax(ii,jj,:)=0._dp
          END IF

          ! diagnostics
          DO jw=1, nw

             zfrac(:) = 0._dp

             DO jclass = in2a, fn2b

                IF (pnaero(ii,jj,jclass) > nlim .AND. zns(jclass) > 0._dp) THEN

                   !-- moles of solute in particle at the upper bound of the bin
                   znshi = zns(jclass)*pnaero(ii,jj,jclass)*vhilim(jclass)/zvols_tot(ii,jj,jclass)

                   !-- critical supersaturation
                   sil = exp(sqrt(zkoehl/znshi)) - 1._dp

                   IF(psmax(ii,jj,jw) < sil) CYCLE
             
                   !-- moles of solute at the lower bound of the bin:
                   znslo = zns(jclass)*pnaero(ii,jj,jclass)*vlolim(jclass)/zvols_tot(ii,jj,jclass)
             
                   !-- critical supersaturation
                   siu = exp(sqrt(zkoehl/znslo)) - 1._dp
             
                   !-- fraction of activated in a bin, eq (13), part 3
                   zfrac(jclass) = min(1._dp,log(psmax(ii,jj,jw)/sil)/log(siu/sil))
             
                END IF
          
             END DO ! jclass

             !-- number of cloud droplets in a grid cell
             pcdncact(ii,jj) = pcdncact(ii,jj) + sum(zfrac(in2a:fn2b)*pnaero(ii,jj,in2a:fn2b))*pwpdf(ii,jj,jw)
             
          END DO ! jw
          
          ! according to ham_activ_abdulrazzak_ghan,
          ! normalizing with the integral over the pdf (why is pwpdf not already normalized?)
          pcdncact(ii,jj) = pcdncact(ii,jj)/sum(pwpdf(ii,jj,:))

       END DO ! ii
       
    END DO ! jj

  END SUBROUTINE cloud_activation

  SUBROUTINE salsa_abdul_razzak_ghan(&
       kproma,   kbdim,   klev,  krow,  ktdia, &
       pcdncact, pesw,    prho,                &
       pxtm1,    ptm1,    papm1, pqm1,         &
       pw,       pwpdf, &!  pa,    pb,    prdry, &
       !pnact,    pfracn,  
       psc,&!   prc,   
       psmax  &
    )

    USE mo_species, ONLY :                      &
         speclist

    USE mo_ham_species, ONLY:                   &
         id_so4,                                &
         id_ss,                                 &
         id_oc,                                 &
         id_du,                                 &
         id_bc

    USE mo_kind,        ONLY : dp

    USE mo_tracdef,     ONLY: ntrac

    USE mo_ham_salsactl, ONLY :                      &
         in1a, in2a, in2b, fn1a, fn2a, fn2b,         & ! size regime bin indices
         iso4b, iocb, ibcb, idub, issb                 ! tracer look-up tables

    USE mo_ham, ONLY:                                &
         aerocomp, sizeclass, nclass, nsoa, subm_naerospec_nowat

    USE mo_activ,        ONLY: nw

    USE mo_ham_vbsctl, ONLY: &
         t_vbs_group,        & ! the VBS data structure
         vbs_ngroup,         & ! number of VBS bins
         vbs_set,            & ! VBS
         laqsoa,             & ! whether or not aqsoa is on
         t_aq_soa,           & ! the aqsoa data structure
         aqsoa_ngroup,       & ! same as number of VBS bins?
         aqsoa_set             ! aqsoa set

    USE mo_species,      ONLY: speclist
    !--- Arguments:

    IMPLICIT NONE

    INTEGER, INTENT(in)   :: kproma, kbdim, klev, krow, ktdia

    REAL(dp), INTENT(out) :: &
         pcdncact(kbdim,klev),       & ! number of activated particles
         psc(kbdim,klev,nclass),     & ! critical supersaturation [% 0-1]
         psmax(kbdim,klev,nw)          ! maximum supersaturation

    REAL(dp), INTENT(in)  :: &
         ptm1(kbdim,klev),           & ! temperature
         papm1(kbdim,klev),          & ! pressure 
         prho(kbdim,klev),           & ! air density
         pqm1(kbdim,klev),           & ! specific humidity
         pesw(kbdim,klev),           & ! saturation water vapour pressure
         pw(kbdim,klev,nw),          & ! mean or bins of updraft velocity (>0.0) [m/s]
         pwpdf(kbdim,klev,nw),       & ! pdf of updraft velocity [s/m]
         pxtm1(kbdim,klev,ntrac)       ! tracer mixing ratios at t-dt


    !-- local variables ----------
    REAL(dp) ::                        &
             znaero(kbdim,klev,fn2b),  & ! number concentration  [#/m3]
             zvols(kbdim,klev,fn2b,subm_naerospec_nowat) ! total volume concentrations of each
                                               ! chem. compound in a size bin [m3/m3]

    TYPE(t_vbs_group), POINTER :: zgroup      ! local copy for easier reading
    TYPE(t_aq_soa), POINTER :: zgroupaq       ! local copy for easier reading

    INTEGER :: ii, jt, jn                     ! for indexing
    INTEGER :: jl, jd, jv


    !--- 3.2 Particle mass:
    ! calculating zvols from pxtm1:
  
    zvols(1:kproma,:,:,:) = 0.0_dp

    DO ii = in1a, fn1a
       !--- Sulfate volume
       jt = aerocomp(iso4b(ii))%idt
       zvols(1:kproma,:,ii,1) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_so4)%density

       !--- Organic carbon
       jt = aerocomp(iocb(ii))%idt
       zvols(1:kproma,:,ii,2) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_oc)%density
    END DO

    DO ii = in2a, fn2a
       !--- Sulfate volume
       jt = aerocomp(iso4b(ii))%idt
       zvols(1:kproma,:,ii,1) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_so4)%density

       !--- Organic carbon
       jt = aerocomp(iocb(ii))%idt
       zvols(1:kproma,:,ii,2) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_oc)%density

       !--- Black carbon
       jt = aerocomp(ibcb(ii))%idt
       zvols(1:kproma,:,ii,3) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_bc)%density

       !--- Sea salt
       jt = aerocomp(issb(ii))%idt
       zvols(1:kproma,:,ii,4) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_ss)%density

       !--- Mineral dust
       jt = aerocomp(idub(ii))%idt
       zvols(1:kproma,:,ii,5) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_du)%density
    END DO

    DO ii = in2b, fn2b
       !--- Sulfate volume
       jt = aerocomp(iso4b(ii))%idt
       zvols(1:kproma,:,ii,1) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_so4)%density

       !--- Organic carbon
       jt = aerocomp(iocb(ii))%idt
       zvols(1:kproma,:,ii,2) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_oc)%density

       !--- Black carbon
       jt = aerocomp(ibcb(ii))%idt
       zvols(1:kproma,:,ii,3) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_bc)%density

       !--- Mineral dust
       jt = aerocomp(idub(ii))%idt
       zvols(1:kproma,:,ii,5) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(id_du)%density
    END DO
    
    ! >> thk: converting VBS groups
    ! REMARK: vbs is not entirely implemented yet in A-R&G, therefore we temporarily
    ! lump all VBS species into OC (id in pvols is 2) for now:
    IF (nsoa == 2) THEN
       DO ii = 1,nclass
          IF ( sizeclass(ii)%lsoainclass ) THEN
             DO jv = 1, vbs_ngroup
                zgroup => vbs_set(jv)
                ! >> thk:
                IF (zgroup%id_vols == 2) CYCLE
                ! << thk
                jt = aerocomp(zgroup%idx(ii))%idt
!                zvols(1:kproma, :, ii, zgroup%id_vols) = &
!                     prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(zgroup%spid)%density
                zvols(1:kproma, :, ii, 2) =  zvols(1:kproma, :, ii, 2)+&
                     prho(1:kproma,:)*pxtm1(1:kproma,:,jt)/speclist(zgroup%spid)%density
             END DO
          END IF
       END DO
       IF (laqsoa) THEN    
          DO ii = 1,nclass
             IF ( sizeclass(ii)%lsoainclass ) THEN
                DO jv = 1, aqsoa_ngroup
                   zgroupaq => aqsoa_set(jv)
                   ! >> thk:
                   IF (zgroupaq%id_aqsoa == 2) CYCLE
                   ! << thk
                   jt = aerocomp(zgroupaq%idx(ii))%idt
!                   zvols(1:kproma, :, ii, zgroupaq%id_aqsoa) = &
!                        prho(1:kproma,:)*pxtm1(1:kproma, :, jt)/speclist(zgroupaq%spid)%density
                   zvols(1:kproma, :, ii, 2) =  zvols(1:kproma, :, ii, 2)+&
                        prho(1:kproma,:)*pxtm1(1:kproma, :, jt)/speclist(zgroupaq%spid)%density
                END DO
             END IF
          END DO
       END IF
    END IF
    ! << thk

    !--- 3.3) Particle numbers:
    ! calculating znaero from pxtm1:

    DO jn=1, nclass
       jt = sizeclass(jn)%idt_no
       znaero(1:kproma,:,jn) = prho(1:kproma,:)*pxtm1(1:kproma,:,jt)
    END DO

    CALL cloud_activation(&
         kproma,   kbdim, klev,  krow,  ktdia,    &
         znaero,   zvols, ptm1,  papm1, pqm1,     &
         pcdncact, pw,    pwpdf, psc,  psmax      &
         )

  END SUBROUTINE salsa_abdul_razzak_ghan

END MODULE mo_ham_salsa_cloud
