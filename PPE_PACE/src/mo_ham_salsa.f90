!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_salsa.f90
!!
!! \brief
!! Main routine for SALSA microphysics.
!! 
!!
!! \author Harri Kokkola (FMI)
!!
!! \responsible_coder
!! Harri Kokkola, harri.kokkola@fmi.fi
!!
!! \revision_history
!!   -# H. Kokkola (FMI) - original code (2006-2014)
!!
!! \limitations
!! None
!!
!! \details
!! The main routine of aerosol microphysics module SALSA. Routine call 
!! subroutines for microphysical processes sequentially. After microphysical
!! processing, particles are remapped to size distribution.
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

MODULE mo_ham_salsa

  PRIVATE

  ! -- subroutines
  PUBLIC :: salsa

CONTAINS

  SUBROUTINE salsa(&
       kproma,   kbdim,   klev,    krow,                &
       ppres,    prh,     ptemp,   ptstep,              &
       pc_gas,   paerml,  pnaero,                       &
       pm6rp,    pm6dry,  prhop,   pww,    ppbl         &
    )

    USE mo_ham_salsa_properties
    USE mo_ham_salsa_dynamics
    USE mo_ham_salsa_nucleation
    USE mo_ham_salsa_update
    USE mo_ham_salsa_init

    USE mo_physical_constants, ONLY: avo

    USE mo_math_constants,     ONLY: pi_6

    USE mo_ham,                ONLY: naerocomp, lscond, lscoag

    USE mo_ham_salsactl,       ONLY : &
         in1a,  in2a,               & ! size section and composition indices
         in2b,  fn1a,               &
         fn2a,  fn2b,               &
         dpmid,                     &
         nlim,                      &
         recalc                       ! logical switch for recalculation of zdwet

    USE mo_ham_salsactl, ONLY: &
         iso4b,&
         iocb,&
         ibcb,&
         idub,&
         issb
    USE mo_ham_salsa_trac, ONLY:    &
         idt_ms4,                   &
         idt_moc,                   &
         idt_mbc,                   &
         idt_mss,                   &
         idt_mdu

    USE mo_kind, ONLY : dp

    ! >> thk: VBS
    USE mo_ham_vbsctl, ONLY: &
         t_vbs_group,        & ! the VBS data structure
         vbs_ngroup,         & ! number of VBS bins
         vbs_set,            & ! VBS
         laqsoa,             & ! Wet SOA cotrol 
         t_aq_soa,           & ! the aqsoa data structure
         aqsoa_ngroup,       & ! same as number of VBS bins?
         aqsoa_set             ! aqsoa set

    USE mo_ham, ONLY: &
         nsoa,               & ! SOA scheme control switch
         subm_ngasspec,      & ! number of gas species in HAM
         subm_naerospec_nowat,&! number of non-water aerosol species
         nclass,             & ! number of bins
         sizeclass             ! aerosol size classes

    USE mo_species, ONLY: &
         speclist              ! list of chemical species in ham

    USE mo_ham_species, ONLY: &
         id_so4, id_oc, id_ss, id_bc, id_du, id_wat

    USE mo_ham_subm_species, ONLY: &
         isubm_so4g,         & ! pc_gas id for sulfuric acid
         isubm_ocnv            ! pc_gas id for semi-volatile OC
    ! << thk

    IMPLICIT NONE

    !-- Input parameters and variables --------------------------------------------------
    INTEGER, INTENT(in) ::          &
         kproma,                    & ! number of horiz. grid points in a slab (='kproma')
         kbdim,                     & ! dimension for arrays (='kbdim')
         klev,                      & ! number of vertical levels (='klev')
         krow                         ! local latitude index


    REAL(dp), INTENT(in) ::            &
         ppres(kbdim,klev),            & ! atmospheric pressure at each grid point [Pa]
         prh  (kbdim,klev),            & ! relative humidity at each grid point [0-1]
         ptemp(kbdim,klev),            & ! temperature at each grid point [K]
         ptstep                          ! time step [s]


    !-- Input variables that are changed within --------------------------------------
    REAL(dp), INTENT(inout) ::         & ! gas phase concentrations at each grid point [#/m3]
         paerml(kbdim,klev,naerocomp), & ! aerosol mass for individual compounds 
                                         ! [molec. cm-3 for sulfate and ug m-3 for bc, oc, ss, and dust] 
         pc_gas(kbdim,klev,subm_ngasspec),    & ! gas phase tracers
                                ! aerosol distribution properties at each grid point
         pnaero(kbdim,klev,fn2b)   ,& ! number concentration in each size bin [#/m3]
         pm6rp (kbdim,klev,fn2b)   ,& ! mean mode actual radius (wet radius for soluble modes
                                ! and dry radius for insoluble modes) [cm]
         pm6dry(kbdim,klev,fn2a)   ,& ! dry radius [cm]
         pww   (kbdim,klev,fn2b)   ,& ! aerosol water content for each mode [kg(water) m-3(air)]
         prhop (kbdim,klev,fn2b)   ,& ! mean mode particle density [g cm-3]
         ppbl(kbdim)
    

    INTEGER :: zpbl(kbdim)            ! Planetary boundary layer top level

    !-- Output variables -----------------------------------------------------------------

    !   in each bin

    !-- Local variables ------------------------------------------------------------------

    LOGICAL::moving_center
    INTEGER :: ii, jj, kk, nn

    REAL(dp) ::                       &
         zcore   (kbdim,klev,fn2b),   & ! volume of the core in one size bin
         zdwet   (kbdim,klev,fn2b),   &
         zvq     (kbdim,klev,fn2b),   &
         zddry,                       &
         zmvsu 

    REAL(dp) :: &
         zvols (kbdim,klev,fn2b,subm_naerospec_nowat)   ! volume concentrations

    ! >> thk: VBS
    TYPE(t_vbs_group), POINTER :: zgroup      ! local copy for easier reading
    TYPE(t_aq_soa), POINTER :: zgroupaq       ! local copy for easier reading
    INTEGER :: idt_group                      ! for easier indexing
    INTEGER :: jv                             ! index for VBS
    ! << thk

    zmvsu = (speclist(id_so4)%moleweight/1000.)/avo/speclist(id_so4)%density

    zpbl(1:kproma) = int(ppbl(1:kproma))

    moving_center=.false.

    zcore = 0._dp
    zddry = 0._dp
    zvols(1:kproma,:,:,:)=0.0_dp

    !>> convert mass concentrations to volume concentrations for SALSA

    DO ii = in1a, fn1a
       !--- Sulfate volume
       zvols(1:kproma,:,ii,1) = paerml(1:kproma,:,iso4b(ii)) * 1.e6_dp *zmvsu
       !--- Organic carbon
       zvols(1:kproma,:,ii,2) = paerml(1:kproma,:,iocb(ii))/speclist(id_oc)%density * 1.e-9_dp
    END DO

    DO ii = in2a, fn2a
       !--- Sulfate volume
       zvols(1:kproma,:,ii,1) = paerml(1:kproma,:,iso4b(ii)) * 1.e6_dp *zmvsu
       !--- Organic carbon
       zvols(1:kproma,:,ii,2) = paerml(1:kproma,:,iocb(ii))/speclist(id_oc)%density * 1.e-9_dp
       !--- Black carbon
       zvols(1:kproma,:,ii,3) = paerml(1:kproma,:,ibcb(ii))/speclist(id_bc)%density * 1.e-9_dp
       !--- Sea salt
       zvols(1:kproma,:,ii,4) = paerml(1:kproma,:,issb(ii))/speclist(id_ss)%density * 1.e-9_dp
       !--- Mineral dust
       zvols(1:kproma,:,ii,5) = paerml(1:kproma,:,idub(ii))/speclist(id_du)%density * 1.e-9_dp
    END DO

    DO ii = in2b, fn2b
       !--- Sulfate volume
       zvols(1:kproma,:,ii,1) = paerml(1:kproma,:,iso4b(ii)) * 1.e6_dp *zmvsu
       !--- Organic carbon
       zvols(1:kproma,:,ii,2) = paerml(1:kproma,:,iocb(ii))/speclist(id_oc)%density * 1.e-9_dp
       !--- Black carbon
       zvols(1:kproma,:,ii,3) = paerml(1:kproma,:,ibcb(ii))/speclist(id_bc)%density * 1.e-9_dp
       !--- Mineral dust
       zvols(1:kproma,:,ii,5) = paerml(1:kproma,:,idub(ii))/speclist(id_du)%density * 1.e-9_dp
    END DO

    ! >> thk: converting VBS groups
    IF (nsoa == 2) THEN
       DO ii = 1,nclass
          IF ( sizeclass(ii)%lsoainclass ) THEN
             DO jv = 1, vbs_ngroup
                zgroup => vbs_set(jv)
                zvols(1:kproma, :, ii, zgroup%id_vols) = &
                     paerml(1:kproma, :, zgroup%idx(ii))/speclist(zgroup%spid)%density*1e-9_dp
             END DO
          END IF
       END DO
    END IF
    if (laqsoa) THEN    
       DO ii = 1,nclass
          IF ( sizeclass(ii)%lsoainclass ) THEN
             DO jv = 1, aqsoa_ngroup
                zgroupaq => aqsoa_set(jv)
                zvols(1:kproma, :, ii, zgroupaq%id_aqsoa) = &
                     paerml(1:kproma, :, zgroupaq%idx(ii))/speclist(zgroupaq%spid)%density*1e-9_dp
             END DO
          END IF
       END DO
    END IF
    ! << thk




    DO kk = in1a,fn1a      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid
             IF(pnaero(ii,jj,kk) < nlim) CYCLE
             zddry = (sum(zvols(ii,jj,kk,1:2))/pnaero(ii,jj,kk)/pi_6)**(1._dp/3._dp)
             IF(zddry < 1.e-10_dp) THEN
                pc_gas(ii,jj,isubm_so4g)   = pc_gas(ii,jj,isubm_so4g) + &
                     zvols(ii,jj,kk,1) * speclist(id_so4)%density /     &
                     (speclist(id_so4)%moleweight/1000.) * avo
                pc_gas(ii,jj,isubm_ocnv)   = pc_gas(ii,jj,isubm_ocnv) + &
                     zvols(ii,jj,kk,2) * speclist(id_oc)%density /      &
                     (speclist(id_oc)%moleweight/1000.) * avo
                pnaero(ii,jj,kk)  = 0.0_dp
                zvols(ii,jj,kk,:) = 0.0_dp
             END IF
          END DO
       END DO
    END DO

    DO kk = in2a,fn2a      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid
             IF(pnaero(ii,jj,kk) < nlim) CYCLE
             zddry = (sum(zvols(ii,jj,kk,:))/pnaero(ii,jj,kk)/pi_6)**(1._dp/3._dp)
             IF(zddry < 1.e-10_dp) THEN
                pc_gas(ii,jj,isubm_so4g) = pc_gas(ii,jj,isubm_so4g) + &
                     zvols(ii,jj,kk,1) * speclist(id_so4)%density /   &
                     (speclist(id_so4)%moleweight/1000.) * avo
                pc_gas(ii,jj,isubm_ocnv) = pc_gas(ii,jj,isubm_ocnv) + &
                     zvols(ii,jj,kk,2) * speclist(id_oc)%density /    &
                     (speclist(id_oc)%moleweight/1000.) * avo
                pnaero(ii,jj,kk)  = 0.0_dp
                zvols(ii,jj,kk,:) = 0.0_dp
             END IF
          END DO
       END DO
    END DO

    DO kk = in2b,fn2b      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid
             IF(pnaero(ii,jj,kk) < nlim) CYCLE
             zddry = (sum(zvols(ii,jj,kk,:))/pnaero(ii,jj,kk)/pi_6)**(1._dp/3._dp)
             IF(zddry < 1.e-10_dp) THEN
                pc_gas(ii,jj,isubm_so4g) = pc_gas(ii,jj,isubm_so4g) + &
                     zvols(ii,jj,kk,1) * speclist(id_so4)%density /   &
                     (speclist(id_so4)%moleweight/1000.) * avo
                pc_gas(ii,jj,isubm_ocnv) = pc_gas(ii,jj,isubm_ocnv) + &
                     zvols(ii,jj,kk,2) * speclist(id_oc)%density /    &
                     (speclist(id_oc)%moleweight/1000.) * avo
                pnaero(ii,jj,kk)  = 0.0_dp
                zvols(ii,jj,kk,:) = 0.0_dp
             END IF
          END DO
       END DO
    END DO

    CALL equilibration(kproma, kbdim, klev,        &
         pnaero, zvols, prh, ptemp,  &
         zcore,  zdwet) 

    IF (lscoag) THEN

       !--------------------------------------------------------------------------------
       !
       !  Calculate coagulation coefficients for particles:
       !  different set for each vertical (pressure) level
       !
       !  The values are calculated for bin mid-diameters and scaled each
       !  time step according to actual particle wet size
       !
       CALL coagulation(kproma, kbdim,  klev,        &
            pnaero, zvols,  zdwet,       &
            zcore,  ptstep, ptemp, ppres)

    END IF


    !- For more accurate condensation sink, wet diameter must be recalculated
    !  NOTE: this is much more time consuming
    IF (recalc) CALL equilibration(kproma, kbdim, klev,       &
         pnaero, zvols, prh, ptemp, &
         zcore,  zdwet) 

    IF (lscond) CALL condensation(&
         kproma,   kbdim,   klev,    krow,               &
         pnaero,   zvols,   zdwet,                       &
         pc_gas,                                         &
         prh,      ptemp,   ppres,   ptstep, zpbl        &
    )

    DO nn = 1, 2

       CALL distr_update(kproma, kbdim, klev, &
                         pnaero, zvols)

       DO kk = in1a,fn1a      ! size bin
          DO jj = 1,klev      ! vertical grid
             DO ii = 1,kproma ! horizontal grid
                IF(pnaero(ii,jj,kk) < nlim) CYCLE
                zddry = (sum(zvols(ii,jj,kk,1:2))/pnaero(ii,jj,kk)/pi_6)**(1._dp/3._dp)
                IF(zddry < 1.e-10_dp) THEN
                   pc_gas(ii,jj,isubm_so4g) = pc_gas(ii,jj,isubm_so4g) + &
                        zvols(ii,jj,kk,1) * speclist(id_so4)%density /   &
                        (speclist(id_so4)%moleweight/1000.) * avo
                   pc_gas(ii,jj,isubm_ocnv) = pc_gas(ii,jj,isubm_ocnv) + &
                        zvols(ii,jj,kk,2) * speclist(id_oc)%density /    &
                        (speclist(id_oc)%moleweight/1000.) * avo
                   pnaero(ii,jj,kk)  = 0.0_dp
                   zvols(ii,jj,kk,:) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

       DO kk = in2a,fn2a      ! size bin
          DO jj = 1,klev      ! vertical grid
             DO ii = 1,kproma ! horizontal grid
                IF(pnaero(ii,jj,kk) < nlim) CYCLE
                zddry = (sum(zvols(ii,jj,kk,:))/pnaero(ii,jj,kk)/pi_6)**(1._dp/3._dp)
                IF(zddry < 1.e-10_dp) THEN
                   pc_gas(ii,jj,isubm_so4g) = pc_gas(ii,jj,isubm_so4g) + &
                        zvols(ii,jj,kk,1) * speclist(id_so4)%density /   &
                        (speclist(id_so4)%moleweight/1000.) * avo
                   pc_gas(ii,jj,isubm_ocnv) = pc_gas(ii,jj,isubm_ocnv) + &
                        zvols(ii,jj,kk,2) * speclist(id_oc)%density /    &
                        (speclist(id_oc)%moleweight/1000.) * avo
                   pnaero(ii,jj,kk)    = 0.0_dp
                   zvols(ii,jj,kk,:) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

       DO kk = in2b,fn2b      ! size bin
          DO jj = 1,klev      ! vertical grid
             DO ii = 1,kproma ! horizontal grid
                IF(pnaero(ii,jj,kk) < nlim) CYCLE
                zddry = (sum(zvols(ii,jj,kk,:))/pnaero(ii,jj,kk)/pi_6)**(1._dp/3._dp)
                IF(zddry < 1.e-10_dp) THEN
                   pc_gas(ii,jj,isubm_so4g) = pc_gas(ii,jj,isubm_so4g) + &
                        zvols(ii,jj,kk,1) * speclist(id_so4)%density /   &
                        (speclist(id_so4)%moleweight/1000.) * avo
                   pc_gas(ii,jj,isubm_ocnv) = pc_gas(ii,jj,isubm_ocnv) + &
                        zvols(ii,jj,kk,2) * speclist(id_oc)%density /    &
                        (speclist(id_oc)%moleweight/1000.) * avo
                   pnaero(ii,jj,kk)  = 0.0_dp
                   zvols(ii,jj,kk,:) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

       CALL equilibration(kproma, kbdim, klev,        &
            pnaero, zvols, prh, ptemp,  &
            zcore,  zdwet) 

       pm6dry(1:kproma,:,1:fn2a) = (zcore(1:kproma,:,1:fn2a)/pi_6)**(1._dp/3._dp)/2._dp * 100._dp
       pm6rp(1:kproma,:,1:fn2b)  = zdwet(1:kproma,:,1:fn2b)/2._dp * 100._dp
       pww(1:kproma,:,1:fn2b)    = (pi_6*zdwet(1:kproma,:,1:fn2b)**3 -                            &
                                   zcore(1:kproma,:,1:fn2b))*speclist(id_wat)%density*pnaero(1:kproma,:,1:fn2b)

    END DO

    !---------------------------------------------------------------


    DO jj = 1,klev      ! vertical grid
       DO ii = 1,kproma ! horizontal kproma in the slab
          DO kk = in1a, fn1a
             zvq(ii,jj,kk) = (zvols(ii,jj,kk,1) + zvols(ii,jj,kk,2) + &
                  pww(ii,jj,kk)/ speclist(id_wat)%density)
             
             ! check if enough aerosol mass present
             
             IF(zvq(ii,jj,kk) > 1.e-30_dp) THEN
                prhop(ii,jj,kk) = (zvols(ii,jj,kk,1)*speclist(id_so4)%density            &
                     + zvols(ii,jj,kk,2)*speclist(id_oc)%density + pww(ii,jj,kk)) /     &
                     zvq(ii,jj,kk)                                    &
                     ! conversion from kg/m3 to g/cm3
                     /1000._dp
             ELSE
                
                ! if not enough mass, assume density of water in g/cm3 to avoid NaN
                prhop(ii,jj,kk) = speclist(id_wat)%density/1000._dp

             END IF

          END DO

          DO kk = in2a, fn2a
             zvq(ii,jj,kk) = (zvols(ii,jj,kk,1) +                     &
                              zvols(ii,jj,kk,2) +                     &
                              zvols(ii,jj,kk,3) +                     &
                              zvols(ii,jj,kk,4) +                     &
                              zvols(ii,jj,kk,5) +                     &
                              pww(ii,jj,kk)/speclist(id_wat)%density)

             IF(zvq(ii,jj,kk) > 1.e-30_dp) THEN 

                prhop(ii,jj,kk) = (zvols(ii,jj,kk,1)*speclist(id_so4)%density +          &
                                   zvols(ii,jj,kk,2)*speclist(id_oc)%density +          &
                                   zvols(ii,jj,kk,3)*speclist(id_bc)%density +          &
                                   zvols(ii,jj,kk,4)*speclist(id_ss)%density +          &
                                   zvols(ii,jj,kk,5)*speclist(id_du)%density +          &
                                   pww(ii,jj,kk))/                    &
                                   zvq(ii,jj,kk) / 1000._dp

             ELSE

                prhop(ii,jj,kk) = speclist(id_wat)%density/1000._dp

             END IF

          END DO

          DO kk = in2b, fn2b 

             zvq(ii,jj,kk) = (zvols(ii,jj,kk,1) +                      &
                              zvols(ii,jj,kk,2) +                      &
                              zvols(ii,jj,kk,3) +                      &
                              zvols(ii,jj,kk,5) +                      &
                              pww(ii,jj,kk) / speclist(id_wat)%density)

             IF(zvq(ii,jj,kk) > 1.e-30_dp) THEN
                prhop(ii,jj,kk) = (zvols(ii,jj,kk,1)*speclist(id_so4)%density +           &
                                   zvols(ii,jj,kk,2)*speclist(id_oc)%density +           &
                                   zvols(ii,jj,kk,3)*speclist(id_bc)%density +           &
                                   zvols(ii,jj,kk,5)*speclist(id_du)%density +           &
                                   pww(ii,jj,kk))/                     &
                                   zvq(ii,jj,kk)/1000._dp
             ELSE
                prhop(ii,jj,kk) = speclist(id_bc)%density/1000._dp
             END IF

          END DO
       END DO
    END DO

    !>> convert volume concentrations to mass concentrations for HAM

    DO ii = in1a, fn1a
       !--- Sulfate volume
       paerml(1:kproma,:,iso4b(ii)) = zvols(1:kproma,:,ii,1) / (1.e6_dp * zmvsu)
       !--- Organic carbon
       paerml(1:kproma,:,iocb(ii)) = zvols(1:kproma,:,ii,2) * speclist(id_oc)%density * 1.e9_dp
    END DO

    DO ii = in2a, fn2a
       !--- Sulfate volume
       paerml(1:kproma,:,iso4b(ii)) = zvols(1:kproma,:,ii,1)  / (1.e6_dp * zmvsu)
       !--- Organic carbon
       paerml(1:kproma,:,iocb(ii)) = zvols(1:kproma,:,ii,2) * speclist(id_oc)%density * 1.e9_dp
       !--- Black carbon
       paerml(1:kproma,:,ibcb(ii)) = zvols(1:kproma,:,ii,3) * speclist(id_bc)%density * 1.e9_dp
       !--- Sea salt
       paerml(1:kproma,:,issb(ii)) = zvols(1:kproma,:,ii,4) * speclist(id_ss)%density * 1.e9_dp
       !--- Mineral dust
       paerml(1:kproma,:,idub(ii)) = zvols(1:kproma,:,ii,5) * speclist(id_du)%density * 1.e9_dp
    END DO

    DO ii = in2b, fn2b
       !--- Sulfate volume
       paerml(1:kproma,:,iso4b(ii)) = zvols(1:kproma,:,ii,1)  / (1.e6_dp * zmvsu)
       !--- Organic carbon
       paerml(1:kproma,:,iocb(ii)) = zvols(1:kproma,:,ii,2) * speclist(id_oc)%density * 1.e9_dp
       !--- Black carbon
       paerml(1:kproma,:,ibcb(ii)) = zvols(1:kproma,:,ii,3) * speclist(id_bc)%density * 1.e9_dp
       !--- Mineral dust
       paerml(1:kproma,:,idub(ii)) = zvols(1:kproma,:,ii,5) * speclist(id_du)%density * 1.e9_dp
    END DO

    ! >> thk: converting VBS groups
    IF (nsoa == 2) THEN
       DO ii = 1,nclass
          IF ( sizeclass(ii)%lsoainclass ) THEN
             DO jv = 1, vbs_ngroup
                zgroup => vbs_set(jv)
                paerml(1:kproma, :, zgroup%idx(ii))= &
                     zvols(1:kproma, :, ii, zgroup%id_vols)*speclist(zgroup%spid)%density*1e9_dp
             END DO
          END IF
       END DO
    END IF
    if (laqsoa) THEN    
       DO ii = 1,nclass
          IF ( sizeclass(ii)%lsoainclass ) THEN
             DO jv = 1, aqsoa_ngroup
                zgroupaq => aqsoa_set(jv)
                paerml(1:kproma, :, zgroupaq%idx(ii)) = &
                     zvols(1:kproma, :, ii, zgroupaq%id_aqsoa)*speclist(zgroupaq%spid)%density*1e9_dp
             END DO
          END IF
       END DO
    END IF
    ! << thk


  END SUBROUTINE salsa

END MODULE mo_ham_salsa
