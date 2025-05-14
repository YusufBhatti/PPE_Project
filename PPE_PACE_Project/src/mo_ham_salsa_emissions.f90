!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_salsa_emissions.f90
!!
!! \brief
!! Module for HAM-SALSA specific emissions that cannot be dealt with the standard
!! processing scheme.
!!
!! Originally for M7:
!! \author M.G. Schultz and S. Schroeder, (FZ Juelich)
!! \author D. O'Donnell, ETH-Z (ETH Zurich)
!! SALSA:
!! \author A.Laakso, FMI
!!
!! \responsible_coder
!! Anton Laakso, anton.laakso@fmi.fi
!!
!! \revision_history
!!   -# M.G. Schultz and S. Schroeder (FZ Juelich) - original code (2010-02-11)
!!   -# D. O'Donnell (ETH-Zurich) - added AEROCOM emissions (2010-04-21)
!!   -# A.Laakso (FMI) - modified module for SALSA (2013-06)
!!
!! \limitations
!! None
!!
!! \details
!!  Routine ham_salsa_emissions is called from emi_interface and must return
!!  lprocessed=.TRUE. when HAM decided to do something for this sector and 
!!  species. In this case, ham_salsa_emissions must set the boundary conditions
!!  for emission mass fluxes.
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

MODULE mo_ham_salsa_emissions

  USE mo_kind,             ONLY: dp
  USE mo_ham,              ONLY: mw_so2, mw_so4, mw_s
  USE mo_ham_salsactl,     ONLY: in2a,fn2a,in2b, fn2b,in1a,fn1a
  USE mo_ham_dust,         ONLY: ntrace

  IMPLICIT NONE

  PRIVATE

  ! public variables  (see declaration below)

  ! subprograms
  PUBLIC                       :: ham_salsa_init_emissions
  PUBLIC                       :: ham_salsa_emissions
  PUBLIC                       :: ham_salsa_init_dust_emissions

  ! sector indices 
  ! IPCC emissions have separate sectors for various fossil and fire emissions
  ! for AEROCOM emissions use fossil and fire
  !INTEGER           :: ibc_dust =0, ibc_seasalt = 0 !alaak moved to the mo_ham
  !INTEGER           :: idsec_dust
  INTEGER           :: idsec_seasalt
  INTEGER           :: idsec_fire, idsec_ffire, idsec_gfire, idsec_awb
  INTEGER           :: idsec_fossil, idsec_dom, idsec_ene, idsec_ind, idsec_tra, idsec_wst
  INTEGER           :: idsec_ships
  !INTEGER           :: idsec_biogenic            ! prescribed SOA or SOA precursors
  INTEGER           :: idsec_biofuel             ! aerocom sector

!----------------------------------------------------------------------

  ! indices for weighting factors to distribute mass flux among modes
  INTEGER :: idx_ms4(fn2b)    
  INTEGER :: idx_moc(fn2b)  
  INTEGER :: idx_mbc(fn2b)    
  INTEGER :: idx_mss(fn2a)    
  INTEGER :: idx_mdu(fn2b)    
  INTEGER :: idx_ns4(fn2b)    
  INTEGER :: idx_noc(fn2b)   
  INTEGER :: idx_nbc(fn2b)    
  INTEGER :: idx_nss(fn2a)   
  INTEGER :: idx_ndu(fn2b)  
  INTEGER :: idx_ocnv

  !weighting matrix to distribute dust emissions from bgc scheme to salsa bins
  REAL(dp)      ::  memidust(ntrace,fn2b)

!-------------------------------------------------------------------------

  ! Module parameters to distribute emission mass flux across modes and to convert
  ! mass flux in number flux
  REAL(dp), PARAMETER:: zbb_wsoc_perc  = 0.65_dp,       &! Biom. Burn. Percentage of Water Soluble OC (WS
                                                         ! (M.O. Andreae; Talk: Smoke and Climate)
                        zbg_wsoc_perc  = 0.65_dp,       &! Assume same Percentage of WSOC for biogenic OC
                        zom2oc         = 1.4_dp,        &! Mass ratio organic species to organic carbon
                                                         ! (Seinfeld and Pandis, 1998, p709;
                                                         !  Ferek et al., JGR, 1998)
                        zfacso2        = 0.975_dp,      &! factor to scale primary SO4 emissions
                                                         ! AEROCOM assumption 2.5 % of the SO2 emissions
                                                         ! in the from of SO4
                        zso2ts         = mw_s/mw_so2,   &! conversion factor SO2 to S
                        zso2tso4       = mw_so4/mw_so2   ! conversion mass SO2 to mass SO4

  CONTAINS

  !! ---------------------------------------------------------------------------
  !! subroutine to initialize HAM specific emissions: must set appropriate local sector
  !! indices idsec_xyz

  SUBROUTINE ham_salsa_init_emissions(nsectors)

  USE mo_emi_matrix,             ONLY: em_get_sector_info
  USE mo_util_string,            ONLY: toupper
  USE mo_species,                ONLY: speclist, spec_idt, spec_ntrac
  USE mo_ham,                    ONLY: sizeclass,idsec_dust, idsec_biogenic
  USE mo_ham_species,            ONLY: id_bc, id_oc, id_so4, id_du, id_ss

  USE mo_ham_salsa_trac,         ONLY: idt_ms4,   idt_moc,   idt_mbc,   idt_mss,      & ! SALSA indices
                                       idt_mdu,   idt_n, idt_mwa, idt_ocnv
 
  USE mo_ham_dust,               ONLY: bgc_dust_initialize

  INTEGER, INTENT(in)      :: nsectors    ! number of sectors defined in emi matrix

  INTEGER           :: i, nvars, jt, nvars0, jx
  CHARACTER(LEN=64) :: secname

  !-- initialize
  idsec_dust      = -1
  idsec_seasalt   = -1
  idsec_fire      = -1    ! fire emissions all in one (AEROCOM)
  idsec_ffire     = -1    ! forest fire emissions (IPCC)
  idsec_gfire     = -1    ! grass fire emissions (IPCC)
  idsec_awb       = -1    ! agricultural waste burning (AEROCOM and IPCC)
  idsec_fossil    = -1    ! fossil fuel emissions all in one (AEROCOM)
  idsec_dom       = -1    ! domestic (fossil and biofuel) emissions (IPCC)
  idsec_ene       = -1    ! energy sector emissions (IPCC)
  idsec_ind       = -1    ! industry sector emissions (IPCC)
  idsec_tra       = -1    ! traffic emissions (IPCC)
  idsec_wst       = -1    ! waste handling emissions (IPCC)
  idsec_ships     = -1    ! international ship traffic emissions (IPCC)
  idsec_biogenic  = -1
  idsec_biofuel   = -1    ! biofuel emissions (AEROCOM)

  DO i = 1, nsectors
    CALL em_get_sector_info(i, secname, nvars)
    IF ((TRIM(toupper(secname)) == 'DUST') .AND. (nvars > 0))     idsec_dust = i
    IF ((TRIM(toupper(secname)) == 'SEASALT') .AND. (nvars > 0))  idsec_seasalt = i
    IF ((TRIM(toupper(secname)) == 'FIRE') .AND. (nvars > 0))     idsec_fire = i
    IF ((TRIM(toupper(secname)) == 'FFIRE') .AND. (nvars > 0))    idsec_ffire = i
    IF ((TRIM(toupper(secname)) == 'GFIRE') .AND. (nvars > 0))    idsec_gfire = i
    IF ((TRIM(toupper(secname)) == 'AWB') .AND. (nvars > 0))      idsec_awb = i
    IF ((TRIM(toupper(secname)) == 'FOSSIL') .AND. (nvars > 0))   idsec_fossil = i
    IF ((TRIM(toupper(secname)) == 'DOM') .AND. (nvars > 0))      idsec_dom = i
    IF ((TRIM(toupper(secname)) == 'ENE') .AND. (nvars > 0))      idsec_ene = i
    IF ((TRIM(toupper(secname)) == 'IND') .AND. (nvars > 0))      idsec_ind = i
    IF ((TRIM(toupper(secname)) == 'TRA') .AND. (nvars > 0))      idsec_tra = i
    IF ((TRIM(toupper(secname)) == 'WST') .AND. (nvars > 0))      idsec_wst = i
    IF ((TRIM(toupper(secname)) == 'SHIPS') .AND. (nvars > 0))    idsec_ships = i
    IF ((TRIM(toupper(secname)) == 'BIOGENIC') .AND. (nvars > 0)) idsec_biogenic = i
    IF ((TRIM(toupper(secname)) == 'BIOFUEL')  .AND. (nvars > 0)) idsec_biofuel = i
  END DO

  idx_ndu(:)=-1
  idx_mdu(:)=-1
  idx_nss(:)=-1
  idx_mss(:)=-1
  idx_nbc(:)=-1
  idx_mbc(:)=-1
  idx_noc(:)=-1
  idx_moc(:)=-1
  idx_ns4(:)=-1
  idx_ms4(:)=-1
  idx_ocnv=1
  
  !alaak: spec_idt points specific species' specific bin place in the tracers list 
  IF (idsec_dust > 0) THEN
  !-- call submodel specific initialisations
    CALL bgc_dust_initialize !alaak: oisko se ok?   
   CALL ham_salsa_init_dust_emissions
!spec_ntrac(id_du) is set to 2 subregion
!output of number tracers is done by BYNUMMODE
    nvars0 = spec_ntrac(id_du)+1
    nvars = spec_ntrac(id_du)+(fn2b-fn1a)
    spec_ntrac(id_du) = nvars
    !in2a -> fn2b
    DO i = nvars0, nvars
    spec_idt(id_du, i) = sizeclass(in2a+(i-nvars0))%idt_no
    idx_ndu(in2a+(i-nvars0))=i
    END DO
    DO jt=1,nvars0-1
        DO jx=in2a,fn2b 
              IF (spec_idt(id_du, jt) == idt_mdu(jx)) idx_mdu(jx) = jt
        END DO
    END DO
  END IF

  !siisolt
  IF (idsec_seasalt > 0) THEN
    nvars0 = spec_ntrac(id_ss)+1
    nvars = spec_ntrac(id_ss)+(fn2a-fn1a)
    !in2a -> fn2a
    DO i = nvars0, nvars
        spec_idt(id_ss, i) = sizeclass(in2a+(i-nvars0))%idt_no
        idx_nss(in2a+(i-nvars0))=i
    END DO
    spec_ntrac(id_ss) = nvars
    DO jt=1,nvars0-1
        DO jx=in2a,fn2a 
              IF (spec_idt(id_ss, jt) == idt_mss(jx)) idx_mss(jx) = jt
        END DO
    END DO
  END IF

  !-- expand spec_idt list to include aerosol number tracers and find appropriate aerosol mass tracers
  !   for the distribution of emission among aerosol modes
  !-- black carbon
  nvars0 = spec_ntrac(id_bc)+1
  nvars = spec_ntrac(id_bc)+(fn2b-fn1a)
  DO i = nvars0, nvars
    spec_idt(id_bc, i) = sizeclass(in2a+(i-nvars0))%idt_no
    idx_nbc(in2a+(i-nvars0))=i
  END DO
  spec_ntrac(id_bc) = nvars
   DO jt=1,nvars0-1
        DO jx=in1a,fn2b 
              IF (spec_idt(id_bc, jt) == idt_mbc(jx)) idx_mbc(jx) = jt
        END DO
   END DO

  !-- organic carbon
  nvars0 = spec_ntrac(id_oc)+1
  nvars = spec_ntrac(id_oc)+(fn2b-in1a+1)
  DO i = nvars0, nvars
    spec_idt(id_oc, i) = sizeclass(in1a+(i-nvars0))%idt_no
    idx_noc(in1a+(i-nvars0))=i
  END DO
  spec_ntrac(id_oc) = nvars
    DO jt=1,nvars0-1
        DO jx=in1a,fn2b 
              IF (spec_idt(id_oc, jt) == idt_moc(jx)) idx_moc(jx) = jt
            IF (spec_idt(id_oc, jt) == idt_ocnv) idx_ocnv = jt 
        END DO
    END DO

  !-- sulfur
  nvars0 = spec_ntrac(id_so4)+1
  nvars = spec_ntrac(id_so4)+(fn2b-in1a+1)
    !in2a -> fn2a !emissions go to soluble bins
  DO i = nvars0, nvars
    spec_idt(id_so4, i) = sizeclass(in1a+(i-nvars0))%idt_no
    idx_ns4(in1a+(i-nvars0))=i
  END DO
  spec_ntrac(id_so4) = nvars

    DO jt=1,nvars0-1
        DO jx=in1a,fn2b 
              IF (spec_idt(id_so4, jt) == idt_ms4(jx)) idx_ms4(jx) = jt
        END DO
    END DO

  END SUBROUTINE ham_salsa_init_emissions

  SUBROUTINE ham_salsa_init_dust_emissions
  !Calculates weighting matrix for distributing dust emissions from bins of MPI BGC dust     
  !emissionc scheme to SALSA bins. Necessary do only in the beginning of the run.
  USE mo_ham_dust,      ONLY: bgc_dust_calc_emis, flux_6h, Dmin, Dmax,nbinit,Dstep
  USE mo_ham_salsactl,  ONLY: vhilim, vlolim

  USE mo_math_constants,ONLY: pi_6

  IMPLICIT NONE
  INTEGER        ::  nclassdu, ii, kk,nn, ntracemax
  REAL(dp)        ::  dbmin(ntrace), dbmax(ntrace), dpk(ntrace), vlobgc(ntrace), vhibgc(ntrace)
  REAL(dp)        ::  rdp, dlast
  
    nclassdu=ntrace*nbinit
    nn=1     
    rdp  =Dmin  ! minimum particules diameter (cm)
    dlast=Dmin
    dpk  (:)=0._dp
    dbmin(:)=0._dp
    dbmax(:)=0._dp
    memidust(:,:) = 0._dp

    DO kk=1,nclassdu                              ! assign fluxes to bins
       IF (mod(kk,nbinit).eq.0) THEN
          dbmax(nn)=rdp*0.01_dp*0.5_dp      ! calculate bin minimum/maximum radius in m
          dbmin(nn)=dlast*0.01_dp*0.5_dp     
          dpk(nn)=sqrt(dbmax(nn)*dbmin(nn))
          nn=nn+1
          dlast=rdp
       ENDIF
       rdp = rdp * exp(Dstep)
    ENDDO !kk      
    
    DO nn=1,ntrace
        vlobgc(nn)=pi_6*dbmin(nn)**3_dp
        vhibgc(nn)=pi_6*dbmax(nn)**3_dp
    ENDDO
    !DO ii=in2b,fn2b
   !     dpsmin(ii) = ((vlolim(ii))/(pi_6))**(1._dp/3._dp)
   !     dpsmax(ii) = ((vhilim(ii))/(pi_6))**(1._dp/3._dp)
!    END DO
    ntracemax=1
       DO nn = 1, ntrace  
        IF (vhibgc(nn)<=vhilim(fn2b)) THEN
            ntracemax=ntracemax+1
        ENDIF
    ENDDO

    
    DO nn=1,ntracemax
        DO ii=in2b,fn2b
        !///      |      |        |  salsa
         !///xxx|x|    |    |    |    bgc (emissions)
        IF (vhibgc(nn)<=vlolim(in2b)) THEN
            memidust(nn,in2b) = 1.0_dp
         ENDIF
                        
        !///       |      |        | salsa   
        !///   | |x    |    |    |   bgc 
        IF (vlobgc(nn)<=vlolim(in2b) .and. vhibgc(nn)>=vlolim(in2b)) THEN
            memidust(nn,in2b) = (vlolim(in2b)-vlobgc(nn))/(vhibgc(nn)-vlobgc(nn)) 
        ENDIF

        !///      |      |        |  salsa
        !///   | |    |    |xxxx|    bgc 
        IF (vlobgc(nn)>=vlolim(ii) .and. vhibgc(nn)<=vhilim(ii)) THEN
            memidust(nn,ii) = 1.0_dp
        ENDIF

        !///      |       |   |        salsa
        !///   | |    | |  xxx    |  bgc   
         IF (vlobgc(nn)<=vlolim(ii) .and. vhibgc(nn)>=vhilim(ii)) THEN
            memidust(nn,ii) = (vhilim(ii)-vlolim(ii))/(vhibgc(nn)-vlobgc(nn)) 
        ENDIF
        
        !///      |       |    |     salsa
        !///   | |    | |x       |   bgc
         IF (vlobgc(nn)>=vlolim(ii) .and. vlobgc(nn)<=vhilim(ii) .and. vhibgc(nn)>=vhilim(ii)) THEN
            memidust(nn,ii) = (vhilim(ii)-vlobgc(nn))/(vhibgc(nn)-vlobgc(nn)) 
        ENDIF

        !///      |       |    |   | salsa
        !///   | |    | |       x|   bgc        
          IF (vhibgc(nn)>=vlolim(ii) .and. vlolim(ii)>=vlobgc(nn) .and. vhibgc(nn)<=vhilim(ii)) THEN
            memidust(nn,ii) = (vhibgc(nn)-vlolim(ii))/(vhibgc(nn)-vlobgc(nn)) 
        ENDIF
            
        ENDDO
    ENDDO
                                          
END SUBROUTINE ham_salsa_init_dust_emissions                                            

SUBROUTINE ham_salsa_emissions(kproma, kbdim, klev, krow, ktrac, ksec, kspec,     &
                              ibc_extra, pfactor)

  !    --A. Laakso 5-2013
  ! Calculates and returns weighting factors for specific species and sector

  USE mo_time_control,       ONLY: delta_time
  USE mo_submodel_diag,      ONLY: t_diag_list, get_diag_pointer
  USE mo_submodel_streams,   ONLY: emi_lpost
  USE mo_species,            ONLY: nmaxtrspec  !maximum number of tracers per species
  USE mo_ham,                ONLY: aerocomp, nsoa, sizeclass, nseasalt, idsec_dust, idsec_biogenic, &
                                   sigma_fine, sigma_coarse !SF #320
  USE mo_ham_salsactl,       ONLY: locgas
  ! for dust and seasalt emissions
  USE mo_ham_species,        ONLY: id_du, id_ss
  USE mo_ham_salsactl, ONLY:  &
         dpmid,     &
         vhilim,    &
         vlolim
  ! for carbon emissions
  USE mo_species,            ONLY: speclist
  USE mo_ham_species,        ONLY: id_bc, id_oc
  ! for sulfur emissions
  USE mo_ham_species,        ONLY: id_so2, id_so4
  !>>dod
  USE mo_ham_soa,            ONLY: soaprop 
  !<<dod
  !>>dod (redmine #44) import seasalt emission schemes from HAM2
  USE mo_ham_m7_emi_seasalt, ONLY:seasalt_emissions_monahan, seasalt_emissions_lsce,    &
                                  seasalt_emissions_mh,      seasalt_emissions_guelle,  &
                                  seasalt_emissions_gong,    seasalt_emissions_long,    &
                                  seasalt_emissions_gong_SST
  !<<dod
  USE mo_boundary_condition, ONLY: bc_set

  !alaak
  USE mo_ham,   ONLY: ibc_dust, ibc_seasalt

  USE mo_math_constants, ONLY: pi, pi_6

  ! arguments

  INTEGER,  INTENT(in)    :: kproma                   ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                    ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: klev                     ! number of levels
  INTEGER,  INTENT(in)    :: krow                     ! geographic block number
  INTEGER,  INTENT(in)    :: ktrac                    ! number of tracers
  INTEGER,  INTENT(in)    :: ksec                     ! emission sector number (match with idsec_xyz)
  INTEGER,  INTENT(in)    :: kspec                    ! species number
  INTEGER,  INTENT(inout) :: ibc_extra(nmaxtrspec)    ! some sectors (dust, seasalt) use more than one bc
  REAL(dp), INTENT(inout) :: pfactor(nmaxtrspec)      ! weighting factor for size distribution

  ! local variables

  REAL(dp)         :: zmassf_bins(fn2b,kbdim)         ! mass flux of dust for salsa bins
  REAL(dp)         :: zm2n_bins(fn2b)                 ! mass to number conversion factor for bins
  REAL(dp)         :: zmassfss_bins(fn2b,kbdim)       ! mass flux of seasalt for salsa bins
  REAL(dp)         :: znumfss_bins(fn2b,kbdim)        ! number flux of seasalt for salsa bins

  REAL(dp)         :: zmasstemp(3,kbdim), znumtemp(3,kbdim) 

  INTEGER          :: ierr

  INTEGER, PARAMETER :: nmod = 3

  REAL(dp) ::                &
       deltadp,              & ! bin width [m]
       ntot(nmod),           & ! total number concentration of a mode [m-3]
       dpg(nmod),            & ! number median diameter of a mode [m]
       dpg_m(nmod),          & ! mass median diameter of a mode [m]
       sigmag(nmod)            ! standard deviation of a mode

  INTEGER :: ii
  !-- Initialisation

  ibc_extra(:) = 0

  !-- Test for sector-species combinations where HAM should take action

  !---------------------------------------------------------------------------

  !--- Dust emissions
  IF (ksec == idsec_dust .AND. kspec == id_du) THEN
     ! note: id_du = 0 should not occur here. Let it crash if it does ...
     ! calculate dust emissions and map onto HAM modes
     CALL ham_salsa_dust_emissions(kproma, kbdim, krow, zmassf_bins,zm2n_bins)

     DO ii=in2b,fn2b
        ibc_extra(idx_mdu(ii)) = ibc_dust + (ii-in2b)
        ibc_extra(idx_ndu(ii)) = ibc_dust + (ii-in2b)
     ENDDO

     ! next statement for efficiency reasons
     ibc_extra(1) = 0
     DO ii=in2b,fn2b
        CALL bc_set(ibc_dust+(ii-in2b), kproma, krow, zmassf_bins(ii,:))
     ENDDO

     pfactor(:) =0._dp                  ! reset weighting factors
     !dust emissions are weighted already
     DO ii=in2b,fn2b
        pfactor(idx_mdu(ii)) = 1.0_dp
        pfactor(idx_ndu(ii)) = zm2n_bins(ii)
     ENDDO

  END IF

  !--------------------------------------------------------------------------
  !--- Seasalt emissions
  IF (ksec == idsec_seasalt .AND. kspec == id_ss) THEN
     zmassfss_bins(:,1:kproma) = 0._dp
     znumfss_bins(:,1:kproma) = 0._dp
     ibc_extra(1) = 0

     !---emissions from M7 modes to SALSA bins: 
     !M7 seasalt emissions
     SELECT CASE(nseasalt)

     CASE (1)
        CALL seasalt_emissions_monahan(kproma, kbdim, krow, zmasstemp(2,:), zmasstemp(3,:), znumtemp(2,:), znumtemp(3,:))

     CASE (2)
        CALL seasalt_emissions_lsce(kproma, kbdim, krow, zmasstemp(2,:), zmasstemp(3,:), znumtemp(2,:), znumtemp(3,:))
        !>>dod (redmine #44) import of seasalt schemes from HAM2   
        !   CASE(3)      
        !...not yet implemented

     CASE(4)
        CALL seasalt_emissions_mh(kproma, kbdim, krow, zmasstemp(2,:), zmasstemp(3,:), znumtemp(2,:), znumtemp(3,:))

     CASE(5)
        CALL seasalt_emissions_guelle(kproma, kbdim, krow, zmasstemp(2,:), zmasstemp(3,:), znumtemp(2,:), znumtemp(3,:))

     CASE(6)
        CALL seasalt_emissions_gong(kproma, kbdim, krow, zmasstemp(2,:), zmasstemp(3,:), znumtemp(2,:), znumtemp(3,:))

     CASE(7)
        CALL seasalt_emissions_long(kproma, kbdim, krow, zmasstemp(2,:), zmasstemp(3,:), znumtemp(2,:), znumtemp(3,:))

     CASE(8)
        CALL seasalt_emissions_gong_SST(kproma, kbdim, krow, zmasstemp(2,:), zmasstemp(3,:), znumtemp(2,:), znumtemp(3,:))

     END SELECT

     dpg    = (/ 2._dp*.022e-6_dp, 2._dp*.13e-6_dp, 2._dp*.75e-6_dp /)

     sigmag = (/ sigma_fine,          sigma_fine,         sigma_coarse           /)

     dpg_m  = dpg * exp(3._dp * (log(sigmag)**2))    

     DO ii  = in2a, fn2a

        deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
             pi_6**(1._dp/3._dp)

        zmassfss_bins(ii,1:kproma) =zmasstemp(2,1:kproma)*deltadp/         &
             (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(2)))*                    &
             exp(-log(dpmid(ii)/dpg_m(2))**2/                              &
             (2._dp*log(sigmag(2))**2))+zmasstemp(3,1:kproma)*deltadp/     &
             (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(3)))*                    &
             exp(-log(dpmid(ii)/dpg_m(3))**2/                              &
             (2._dp*log(sigmag(3))**2))
        znumfss_bins(ii,1:kproma) = znumtemp(2,1:kproma)*deltadp/          &
             (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(2)))*                    &
             exp(-log(dpmid(ii)/dpg(2))**2/                                &
             (2._dp*log(sigmag(2))**2)) +znumtemp(3,1:kproma)*deltadp/     &
             (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(3)))*                    &
             exp(-log(dpmid(ii)/dpg(3))**2/                                &
             (2._dp*log(sigmag(3))**2)) 

     ENDDO
     !-------------------------------

     DO ii=in2a,fn2a
        ibc_extra(idx_mss(ii)) = ibc_seasalt + (ii-in2a)
        ibc_extra(idx_nss(ii)) = ibc_seasalt + (ii-in2a) + fn2a - fn1a
     ENDDO
     DO ii=in2a,fn2a
        CALL bc_set(ibc_seasalt+(ii-in2a), kproma, krow, zmassfss_bins(ii,:))
        CALL bc_set(ibc_seasalt+(ii-in2a + (fn2a-fn1a)), kproma, krow, znumfss_bins(ii,:))
     ENDDO

     pfactor(:) = 0._dp                  ! reset weighting factors
     pfactor(idx_mss(in2a:fn2a)) = 1._dp
     pfactor(idx_nss(in2a:fn2a)) = 1._dp
  END IF
  !---------------------------------------------------------------------------
  dpg    = (/ 2._dp*.03e-6_dp, 2._dp*0.075e-6_dp , 2._dp*0.5e-6_dp /)
  sigmag = (/ sigma_fine,         sigma_fine           , sigma_coarse /)
  dpg_m  = dpg * exp(3._dp * (log(sigmag)**2))

  !--- Black carbon emissions
  IF (kspec == id_bc) THEN
     pfactor(:) = 0._dp                  ! reset weighting factors
     IF (ksec == idsec_awb .OR. ksec == idsec_fossil .OR. ksec == idsec_dom      .OR. &
          ksec == idsec_ene .OR. ksec == idsec_ind    .OR. ksec == idsec_tra      .OR. &
          ksec == idsec_wst .OR. ksec == idsec_ships  .OR. ksec == idsec_biogenic .OR. &
          ksec == idsec_biofuel) THEN
        !only insoluble
        DO ii  = in2b, fn2b
           !-- width of the size bin
           deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
                pi_6**(1._dp/3._dp)

           !-- Black carbon and insoluble organic emissions according to AEROCOM
           !mfrac_bc
           pfactor(idx_mbc(ii)) = deltadp/                                            &
                (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(1)))*          &
                exp(-log(dpmid(ii)/dpg_m(1))**2/                    &
                (2._dp*log(sigmag(1))**2))
        END DO
     ELSE !(_fire, _ffire, _gfire
        DO ii  = in2a, fn2a
           deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
                pi_6**(1._dp/3._dp)
           !mfrac_bb
           pfactor(idx_mbc(ii)) = deltadp/                                        &
                (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(2)))*                   &
                exp(-log(dpmid(ii)/dpg_m(2))**2/                             &
                (2._dp*log(sigmag(2))**2))
        END DO
     END IF
     pfactor(idx_mbc(in2a:fn2b))=pfactor(idx_mbc(in2a:fn2b))*1._dp/sum(pfactor(idx_mbc(in2a:fn2b)))
     pfactor(idx_nbc(in2a:fn2b))=pfactor(idx_mbc(in2a:fn2b))/(pi_6*dpmid(in2a:fn2b)**3*speclist(id_bc)%density)
  END IF

  !--- Organic carbon emissions
  IF (kspec == id_oc) THEN
     pfactor(:) = 0._dp                  ! reset weighting factors

     !Fossil fuel: all organic carbon mass goes into insoluble aitken mode
     IF (ksec == idsec_fossil .OR. ksec == idsec_ene .OR. ksec == idsec_ind .OR.   &
          ksec == idsec_tra .OR. ksec == idsec_wst .OR. ksec == idsec_ships) THEN
        !insoluble
        !  
        DO ii  = in2b, fn2b
           !-- width of the size bin
           deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
                pi_6**(1._dp/3._dp)     
           !mfrac_ff
           pfactor(idx_moc(ii)) = deltadp/                                        &
                (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(1)))*                   &
                exp(-log(dpmid(ii)/dpg_m(1))**2/                             &
                (2._dp*log(sigmag(1))**2))
        END DO
        pfactor(idx_moc(:))=zom2oc*pfactor(idx_moc(:))*1._dp/sum(pfactor(idx_moc(:)))
        pfactor(idx_noc(:))=pfactor(idx_moc(:))/(pi_6*dpmid(:)**3*speclist(id_oc)%density)
     END IF

     !Biomass:organic carbon mass goes into soluble and insoluble aitken mode
     IF (ksec == idsec_fire .OR. ksec == idsec_ffire .OR. ksec == idsec_gfire     &
          .OR. ksec == idsec_awb .OR. ksec == idsec_dom .OR. ksec ==  idsec_biofuel) THEN
        !soluble
        DO ii  = in1a, fn2a
           deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
                pi_6**(1._dp/3._dp)
           !mfrac_bb
           pfactor(idx_moc(ii)) = deltadp/                                        &
                (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(1)))*                   &
                exp(-log(dpmid(ii)/dpg_m(1))**2/                             &
                (2._dp*log(sigmag(1))**2))*(zbg_wsoc_perc)
        END DO
        !insoluble
        DO ii  = in2b, fn2b
           deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
                pi_6**(1._dp/3._dp)
           !mfrac_bb
           pfactor(idx_moc(ii)) = deltadp/                                        &
                (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(1)))*                   &
                exp(-log(dpmid(ii)/dpg_m(1))**2/                             &
                (2._dp*log(sigmag(1))**2))*(1._dp-zbg_wsoc_perc)
        END DO
        !pfactor(idx_moc(:))=zom2oc*pfactor(idx_moc(:))*1._dp/sum(pfactor(idx_moc(:)))
        pfactor(idx_moc(in1a:fn2a))=zom2oc*pfactor(idx_moc(in1a:fn2a))*1._dp/sum(pfactor(idx_moc(in1a:fn2a)))*(zbg_wsoc_perc)
        pfactor(idx_moc(in2b:fn2b))=zom2oc*pfactor(idx_moc(in2b:fn2b))*1._dp/sum(pfactor(idx_moc(in2b:fn2b)))*(1._dp-zbg_wsoc_perc)
        pfactor(idx_noc(:))=pfactor(idx_moc(:))/(pi_6*dpmid(:)**3*speclist(id_oc)%density)
     END IF

     !--- biogenic
     !--- if SOA scheme from D. O'Donnell is active, the biogenic sector does not emit aerosols
     !    (only precursors)
     IF (ksec == idsec_biogenic .AND. nsoa /= 1) THEN        !>>dod soa

        !pfactor(idx_mocki) = (1._dp-zbg_wsoc_perc)               ! insoluble fraction
        !pfactor(idx_mocks) = 0.5*zbg_wsoc_perc                   ! soluble fraction (aitken)
        !pfactor(idx_mocas) = 0.5*zbg_wsoc_perc                   ! soluble fraction (accumulation)
        !pfactor(idx_nocki) = (1._dp-zbg_wsoc_perc)*zm2n_ocki_bg  ! number for insoluble aitken

        DO ii  = in1a, fn2a
           deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
                pi_6**(1._dp/3._dp)
           !mfrac_bb
           pfactor(idx_moc(ii)) = SUM(deltadp/                                        &
                (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(1:2)))*                   &
                exp(-log(dpmid(ii)/dpg_m(1:2))**2/                             &
                (2._dp*log(sigmag(1:2))**2)))!*(zbg_wsoc_perc)
        END DO

        !insoluble
        DO ii  = in2b, fn2b
           deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
                pi_6**(1._dp/3._dp)
           !mfrac_bb
           pfactor(idx_moc(ii)) = deltadp/                                        &
                (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(1)))*                   &
                exp(-log(dpmid(ii)/dpg_m(1))**2/                             &
                (2._dp*log(sigmag(1))**2))!*(1._dp-zbg_wsoc_perc)
        END DO

        IF(locgas) THEN
           pfactor(idx_ocnv)=zbg_wsoc_perc
        ELSE
           pfactor(idx_moc(in1a:fn2a)) = pfactor(idx_moc(in1a:fn2a))*1._dp &
                / sum(pfactor(idx_moc(in1a:fn2a)))*(zbg_wsoc_perc)
        END IF

        pfactor(idx_moc(in2b:fn2b)) = pfactor(idx_moc(in2b:fn2b))*1._dp &
             / sum(pfactor(idx_moc(in2b:fn2b)))*(1._dp-zbg_wsoc_perc)
        ! note that soluble number concentration is not altered because it is assumed that
        ! soluble biogenic particles immediately condense on other particles 
        pfactor(idx_noc(:))=0 
        pfactor(idx_noc(in2b:fn2b))=pfactor(idx_moc(in2b:fn2b))/(pi_6*dpmid(in2b:fn2b)**3*speclist(id_oc)%density)
     END IF

  END IF

  !--- Sulfur emissions
  ! ### NOTE: lmoz=.true. will require the sulfur tracer in units of kg(SO2) rather than kg(S) !! ###
  ! Note: emi_matrix only contains SO2 as species, but a fraction (1-zso2frac) of SO2 
  ! emissions is assumed to go directly into the (SO4) aerosol phase
  IF (kspec == id_so2) THEN
     pfactor(:) = 0._dp
     pfactor(1) = zfacso2
  END IF
  IF (kspec == id_so4) THEN
     IF (ksec == idsec_ene .OR. ksec == idsec_ships) THEN ! .OR. ksec == idsec_ind?
        DO ii  = in1a, fn2a
           deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
                pi_6**(1._dp/3._dp)
           pfactor(idx_ms4(ii)) = sum(deltadp/                                           &
                (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(2:3)))*                    &
                exp(-log(dpmid(ii)/dpg_m(2:3))**2/                              &
                (2._dp*log(sigmag(2:3))**2)))
        END DO
     ELSE 
        DO ii  = in1a, fn2a
           deltadp = (vhilim(ii)**(1._dp/3._dp)-vlolim(ii)**(1._dp/3._dp))/   &
                pi_6**(1._dp/3._dp)
           pfactor(idx_ms4(ii)) = sum(deltadp/                                       &
                (dpmid(ii)*sqrt(2._dp*pi)*log(sigmag(1:2)))*                  &
                exp(-log(dpmid(ii)/dpg_m(1:2))**2/                            &
                (2._dp*log(sigmag(1:2))**2)))
        END DO
     END IF

     pfactor(idx_ms4(:))=zso2tso4*pfactor(idx_ms4(:))*1._dp/sum(pfactor(idx_ms4(:)))*(1._dp-zfacso2)
     pfactor(idx_ns4(:))=pfactor(idx_ms4(:))/(pi_6*dpmid(:)**3*speclist(id_so4)%density)

  END IF

END SUBROUTINE ham_salsa_emissions

SUBROUTINE ham_salsa_dust_emissions(kproma, kbdim, krow, pmassf_bins,     &
                                   pm2n_bins)

  ! Description:
  ! ------------
  ! *ham_salsa_dust_emissions* is the interface between the MPI BGC dust emission scheme
  ! and the HAM/SALSA aerosol module. Dust number flux is calculated from the mass flux and
  ! both fluxes are stored as ECHAM5 tracers.
  !
  ! Author:
  ! -------
  ! M. Werner, MPI BGC, November 2002
  !
  ! Interface:
  ! ----------
  ! *ham_salsa_dust_emissions* is called from *ham_salsa_emissions* which is called from
  ! *emi_interface*.
  !
  ! Adaptation to ECHAM6-HAMMOZ code structure by Martin Schultz, FZ Juelich (2010-02-11)
  ! A. Laakso, FMI, 2013 for SALSA
  !
  ! Emissions from BGC is distributed to the salsa bins by using distribution matrix 
  ! memidust (ham_salsa_init_dust_emissions)
  !

  USE mo_memory_g3b,     ONLY: slm, glac
  USE mo_species,        ONLY: speclist
  USE mo_ham_species,    ONLY: id_du
  USE mo_exception,      ONLY: finish
  USE mo_ham_dust,       ONLY: bgc_dust_calc_emis, flux_6h
  USE mo_ham_salsactl,   ONLY: dpmid
  USE mo_math_constants, ONLY: pi_6

  IMPLICIT NONE

  INTEGER,  INTENT(in)    :: kproma                   ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                    ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: krow                     ! geographic block number
  REAL(dp), INTENT(out)   :: pmassf_bins(fn2b,kbdim)         ! mass flux of dust for salsa bins

  REAL(dp), INTENT(out)   :: pm2n_bins(fn2b)                  ! mass to number conversion factor for bins


  !--- Local Variables:

  !INTEGER,  PARAMETER :: min_bin = 1                   ! index boundaries for BGC_Dust bin numbers
  !INTEGER,  PARAMETER :: max_bin = 5                   ! dust emissions are calculated for size range 0,1 um -> 220 um

  REAL(dp)            :: densdust                     ! Particle density in kg/m3
  INTEGER             :: nn, ii, kk

  !--- Get density of dust from the trlist:
  densdust=speclist(id_du)%density

  !--- 1) Calculate BGC dust emissions: result is in g m-2 s-1
  CALL bgc_dust_calc_emis(kproma, krow)

  ! Sum over BGC_Dust internal tracer classes for HAM modal scheme:
  !pmassf_ai(:) = 0._dp
  !pmassf_ci(:) = 0._dp
  pmassf_bins(:,:) = 0._dp

  DO nn=1, ntrace
    DO ii = in2b, fn2b
    pmassf_bins(ii,1:kproma) = pmassf_bins(ii,1:kproma) + flux_6h(1:kproma,nn,krow)*memidust(nn,ii)
     ENDDO 
  ENDDO 

  ! Make sure mass flux is positive
  !WHERE (pmassf_bins.le.0._dp) pmassf_bins(:) = 0._dp

  ! Mask out glacier and ocean values on ECHAM grid
  !WHERE (glac(1:kproma,krow) > 0.5_dp .OR. slm(1:kproma,krow) < 0.5_dp)
  !  pmassf_bins(:,1:kproma) = 0._dp
  !END WHERE
  DO kk = 1,kproma
    IF (glac(kk,krow) > 0.5_dp .OR. slm(kk,krow) < 0.5_dp) THEN
        pmassf_bins(:,kk) = 0._dp
    ENDIF
  ENDDO

  ! Reduce flux over glacier grids according to glacier fraction
  !WHERE (glac(1:kproma,krow) > 0.1_dp)
  !  pmassf_bins(:,1:kproma) = pmassf_bins(:,1:kproma) * (1._dp-glac(1:kproma,krow))**3  
  !END WHERE
  DO kk = 1,kproma
    IF (glac(kk,krow) > 0.1_dp) THEN
        pmassf_bins(:,kk) = pmassf_bins(:,kk) * (1._dp-glac(kk,krow))**3  
    ENDIF
  ENDDO

  ! Reduce flux over not-all-land grids according to land fraction
  !WHERE (slm(1:kproma,krow) < 0.99_dp)
  !  pmassf_bins(:,1:kproma)= pmassf_bins(:,1:kproma) * min(slm(1:kproma,krow),0.7_dp)**4
  !END WHERE
  DO kk = 1,kproma
    IF (slm(kk,krow) < 0.99_dp) THEN
        pmassf_bins(:,kk) = pmassf_bins(:,kk)* min(slm(kk,krow),0.7_dp)**4
    ENDIF
  ENDDO

  ! Convert to kg m-2 s-1
  pmassf_bins(:,1:kproma)= MAX(pmassf_bins(:,1:kproma), 0._dp) / 1000._dp 

  !--- Calculate conversion factor from mass flux to number flux for
  !    log-normal distributions:
  pm2n_bins(:)=(pi_6*(dpmid(:))**3*densdust)

END SUBROUTINE ham_salsa_dust_emissions

END MODULE mo_ham_salsa_emissions
