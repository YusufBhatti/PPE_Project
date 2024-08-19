!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_freezing.f90
!!
!! \brief
!! This module contains all subroutines necessary to handle HAM freezing calculations
!!
!! \author Sylvaine Ferrachat (ETH Zurich)
!!
!! \responsible_coder
!! Sylvaine Ferrachat, sylvaine.ferrachat@env.ethz.ch
!!
!! \revision_history
!!   -# Sylvaine Ferrachat (ETH Zurich) - original code (2010-03)
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

MODULE mo_ham_freezing

  USE mo_kind,          ONLY: dp
  USE mo_tracdef,       ONLY: ntrac
  USE mo_time_control,  ONLY : delta_time, time_step_len

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: HAM_IN_setup

  REAL(dp), PARAMETER :: zeps   = EPSILON(1.0_dp)

 CONTAINS

  SUBROUTINE ham_IN_setup(kproma, kbdim, klev, krow,                                  &
                          prho, pxtm1, pxtte, pcdncact,                               &
                          prwetki, prwetai, prwetci,                                  &
                          pfracdusol, pfracduai, pfracduci, pfracbcsol, pfracbcinsol, &
                          pascs, papnx, paprx, papsigx, ld_het                        )

    USE mo_ham_m7ctl,      ONLY : crdiv,                    &  
                                  iaiti,iacci,icoai,        &     
                                  iaccs, inucs    
    USE mo_ham,            ONLY : lhetfreeze, sizeclass, nclass, nham_subm, HAM_M7, HAM_SALSA
    USE mo_activ,          ONLY: nfrzmod
    USE mo_ham_streams,    ONLY: rwet, ndusol_strat, nbcsol_strat,       &
                                 nduinsolai, nduinsolci, nbcinsol, &
                                 naerinsol
    USE mo_ham_salsactl,   ONLY: in2b, in2a, in1a
    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(in)  :: prho(kbdim,klev)     ! air density
    REAL(dp), INTENT(IN)  :: pxtm1(kbdim,klev,ntrac)
    REAL(dp), INTENT(IN)  :: pxtte(kbdim,klev,ntrac) 
    REAL(dp), INTENT(IN)  :: pcdncact(kbdim,klev)  ! number of activated particules
    !-- for mixed-phase freezing:
    REAL(dp), INTENT(out) :: prwetki(kbdim,klev)  ! wet radius, aitken insoluble mode
    REAL(dp), INTENT(out) :: prwetai(kbdim,klev)  ! wet radius, accumulation insoluble mode
    REAL(dp), INTENT(out) :: prwetci(kbdim,klev)  ! wet radius, coarse insoluble mode
    REAL(dp), INTENT(out) :: pfracdusol(kbdim,klev)   ! total fraction of dust particules in all soluble modes
    REAL(dp), INTENT(out) :: pfracduai(kbdim,klev)    ! fraction of dust particules in the accum. soluble mode 
    REAL(dp), INTENT(out) :: pfracduci(kbdim,klev)    ! fraction of dust particules in the coarse soluble mode 
    REAL(dp), INTENT(out) :: pfracbcsol(kbdim,klev)   ! total fraction of BC particules in all soluble modes 
    REAL(dp), INTENT(out) :: pfracbcinsol(kbdim,klev) ! total fraction of BC particules in all insoluble modes 
    !-- for cirrus freezing:
    REAL(dp), INTENT(out) :: pascs(kbdim,klev)    ! soluble aerosol number conc.
    REAL(dp), INTENT(out) :: papnx(kbdim,klev,nfrzmod)   ! aerosol number available for freezing [1/cm3] 
    REAL(dp), INTENT(out) :: paprx(kbdim,klev,nfrzmod)   ! radius of aerosols avail. for freezing  [cm] 
    REAL(dp), INTENT(out) :: papsigx(kbdim,klev,nfrzmod) ! std. dev. of aerosols available for freezing
    LOGICAL, INTENT(out)  :: ld_het  !switch to set heterogeneous freezing below 235K (cirrus scheme)

    INTEGER  :: jclass, jtn, idx
    REAL(dp) :: zdtime
    REAL(dp) :: ztmst
    REAL(dp) :: zlim_min

    zdtime = delta_time
    ztmst  = time_step_len

    !--- Mixed-phase freezing setup:

     !--- Wet radii:
     SELECT CASE(nham_subm)
         CASE(HAM_M7)
             prwetki(1:kproma,:) = rwet(iaiti)%ptr(1:kproma,:,krow)
             prwetai(1:kproma,:) = rwet(iacci)%ptr(1:kproma,:,krow)
             prwetci(1:kproma,:) = rwet(icoai)%ptr(1:kproma,:,krow)
         CASE(HAM_SALSA)
             prwetki(1:kproma,:) = rwet(in2b)%ptr(1:kproma,:,krow)
             prwetai(1:kproma,:) = rwet(in2b+2)%ptr(1:kproma,:,krow)
             prwetci(1:kproma,:) = rwet(in2b+4)%ptr(1:kproma,:,krow)
     END SELECT

     !--- Get specific aerosol number concentrations for diagnostics and freezing
     CALL get_aerofreez_nc(kproma, kbdim, klev, krow, prho, pxtm1)
 
     !--- Various dust and BC fractions:
     pfracdusol(1:kproma,:)   = MIN(ndusol_strat(1:kproma,:,krow)/(pcdncact(1:kproma,:)      +zeps), 1._dp)
     pfracduai(1:kproma,:)    = MIN(nduinsolai(1:kproma,:,krow)  /(naerinsol(1:kproma,:,krow)+zeps), 1._dp)
     pfracduci(1:kproma,:)    = MIN(nduinsolci(1:kproma,:,krow)  /(naerinsol(1:kproma,:,krow)+zeps), 1._dp)
     pfracbcsol(1:kproma,:)   = MIN(nbcsol_strat(1:kproma,:,krow)/(pcdncact(1:kproma,:)      +zeps), 1._dp)
     pfracbcinsol(1:kproma,:) = MIN(nbcinsol(1:kproma,:,krow)    /(naerinsol(1:kproma,:,krow)+zeps), 1._dp)
 
    !--- Cirrus freezing setup:

     !--- Soluble aerosol number concentration:
     pascs(1:kproma,:) = 0._dp

     DO jclass=1,nclass
        jtn = sizeclass(jclass)%idt_no !mode number tracer identity
        IF (sizeclass(jclass)%lsoluble .AND. jclass .NE. inucs) THEN !soluble modes except nucleation
            pascs(1:kproma,:) = pascs(1:kproma,:) + pxtm1(1:kproma,:,jtn) + pxtte(1:kproma,:,jtn)*ztmst
        ENDIF
     ENDDO

     pascs(1:kproma,:) = MAX(pascs(1:kproma,:), 10.E6_dp)

     !-- Number, radius and std. dev. of aerosols available for cirrus freezing
     !   note: the number will be later reduced by the actual icnc in the cloud routine

!SF note: the following is relying on the fact that nfrzmod=1, ie one uses a monodisperse
!         aerosol distribution for the cirrus freezing calculations (see mo_ham_init.f90).
!         It has been shown that choosing monodisperse distrib is saving a lot of computing time,
!         while not significantly degrading the results, as compared to a more physical
!         3-compounds distribution. This latter possibility has been therefore 'de-implemented' from the code
!         while the parameter nfrzmod is still present. 
!         If this is needed for some testings, this will have to be re-implemented here.

     SELECT CASE(nham_subm)
         CASE(HAM_M7)
             idx      = iaccs
             zlim_min = crdiv(3)
         CASE(HAM_SALSA)
             idx      = in2a+2
             zlim_min = 5.e-6_dp !SF <-- please avoid hardcoding!!
     END SELECT 
     
     ld_het = lhetfreeze

     IF (ld_het) THEN
        papnx(1:kproma,:,1) = ndusol_strat(1:kproma,:,krow)
        
        paprx(1:kproma,:,1) = 100._dp*rwet(idx)%ptr(1:kproma,:,krow) ![cm]
        paprx(1:kproma,:,1) = MAX(paprx(1:kproma,:,1), zlim_min)
        
        papsigx(1:kproma,:,1) = 1.0_dp
        
     ELSE !ld_het=.false.
        
        papnx(1:kproma,:,1) = prho(1:kproma,:) * pascs(1:kproma,:) ![1/m3]
        
        paprx(1:kproma,:,1) = 100._dp * rwet(idx)%ptr(1:kproma,:,krow) ![cm]
        paprx(1:kproma,:,1) = MAX(paprx(1:kproma,:,1), zlim_min)
        
        papsigx(1:kproma,:,1) = 1._dp
        
     ENDIF

  END SUBROUTINE HAM_IN_setup

  !---------------------------------------------------------------------------
  !>
  !! @brief Helper routine for get_aerofreez_nc
  !! 
  !! @remarks computes the ratio of mass (resp. volume)
  !!                     of a given aerosol species in a given mode
  !!                     over the whole mass (resp. vol) of this mode
  !!

  SUBROUTINE aero_massvolratio(kproma, kbdim,  klev,                &
                               thismod, thisspec, calc_type, pxtm1, &
                               paero_massvolratio)
 
     USE mo_ham,                 ONLY: naerocomp, aerocomp
     USE mo_exception,           ONLY: finish
     !>>dod(S) removed id_so2...etc <<dod
  
     INTEGER, INTENT(IN) :: kproma, kbdim, klev, thismod, thisspec
 
     REAL(dp), INTENT(in) :: pxtm1(kbdim, klev, ntrac)
     REAL(dp), INTENT(out) :: paero_massvolratio(kbdim,klev)
     CHARACTER(LEN=4), INTENT(IN) :: calc_type

     INTEGER  :: jn, jclass, jspec, jt
     REAL(dp) :: zratio(kbdim,klev) 
     REAL(dp) :: zdenom(kbdim,klev), znum(kbdim,klev)
     REAL(dp) :: zdens_rcp
 
 
     LOGICAL :: ll(kbdim,klev)
 
     !--1. Main Calc:
     zdenom(1:kproma,:) = 0._dp
 
     DO jn = 1,naerocomp    !loop over all mode-species
 
        jclass = aerocomp(jn)%iclass
        jspec  = aerocomp(jn)%spid
        jt     = aerocomp(jn)%idt
 
        IF ( jclass == thismod ) THEN  !select relevant mode

           !>>dod(S) removed sulfate mass conversion <<dod

           !---set densities or a dummy variable in case of mass ratio calc
           SELECT CASE (TRIM(ADJUSTL(calc_type)))
             CASE ('volu')  !volume ratio calculation
                zdens_rcp = 1000._dp/aerocomp(jn)%species%density !SF 1000. is for the unit conversion from 
                                                                  ! kg.m-3 to g.cm-3
             CASE ('mass') !mass ratio calculation
                zdens_rcp = 1._dp
             CASE DEFAULT
                CALL finish('aero_massvolratio','Wrong calc_type argument',1)
           END SELECT 
 
           !---Add up the different contribs to the denominator 
           !DN #296: put back this line *outside the IF statement*  as it should:
           zdenom (1:kproma,:) = zdenom(1:kproma,:) + zdens_rcp*pxtm1(1:kproma,:,jt)

           !---Computes the numerator
           IF ( jspec == thisspec ) THEN  !select relevant species
              znum(1:kproma,:) = zdens_rcp*pxtm1(1:kproma,:,jt)
           ENDIF 
 
        ENDIF !end if thismod
  
     ENDDO !end loop over all mode-species
 
     !--2. Final calculation:
     ll(1:kproma,:)     = (zdenom(1:kproma,:) > zeps)
     zdenom(1:kproma,:) = MERGE(zdenom(1:kproma,:), 1._dp, ll(1:kproma,:)) !SF to prevent div by ~zero
                                                                           !in  this latter case, the ratio
                                                                           !must be 0 anyway (see right below)
 
     zratio(1:kproma,:) = znum(1:kproma,:) / zdenom(1:kproma,:)
     zratio(1:kproma,:) = MERGE(zratio(1:kproma,:), 0._dp, ll(1:kproma,:)) 
 
     ll(1:kproma,:) = (zratio(1:kproma,:) < zeps)
     IF (ANY(ll(1:kproma,:))) zratio(1:kproma,:)  = MERGE(0._dp, zratio(1:kproma,:), ll(1:kproma,:))
 
     paero_massvolratio (1:kproma,:) = zratio(1:kproma,:)

     IF(kproma < kbdim) paero_massvolratio(kproma+1:kbdim,:)=0._dp
 
   END SUBROUTINE aero_massvolratio
 
   !---------------------------------------------------------------------------
   !>
   !! @brief Helper routine for get_aerofreez_nc
   !! 
   !! @remarks computes the number concentration of a given species
   !!                 by using the surface weighting method,
   !!                 provided its mass or vol ratio as calculated by aero_massvolratio
   !!                 and the total number concentration of the corresponding set
   !!                 of aerosols from which it is extracted 
   !!                 (in general this set is the tot aerosol content of the corresponding mode,
   !!                  but it can also be a subset of that, e.g. the CCNs) 
   !!
 
   SUBROUTINE aero_nc_surfw(kproma, kbdim,  klev, pvolr, pnum, &
                            paero_nc_surfw)
 
     INTEGER, INTENT(IN)   :: kproma, kbdim, klev
     REAL(dp),INTENT(IN)   :: pvolr(kbdim,klev)
     REAL(dp), INTENT(IN)  :: pnum(kbdim,klev)
     REAL(dp), INTENT(OUT) :: paero_nc_surfw(kbdim,klev)
 
     paero_nc_surfw(1:kproma,:) = pvolr(1:kproma,:)**(2._dp/3._dp) &
                                * pnum(1:kproma,:) 

   END SUBROUTINE aero_nc_surfw
 
   !---------------------------------------------------------------------------
   !>
   !! @brief Get specific aerosol number concentrations
   !! 
   !! @remarks loop over all needed aerosol species and
   !!          computes specific number concentrations requested
   !!          in the freezing calculations (both contact and immersion)
   !!
  
  SUBROUTINE get_aerofreez_nc(kproma, kbdim,  klev, krow, prho, pxtm1)
  
  USE mo_species,        ONLY: speclist
  USE mo_ham_species,    ONLY: id_bc, id_du
  USE mo_ham_m7ctl,      ONLY: iaiti, iacci, icoai, iaits, iaccs, icoas
  USE mo_ham,            ONLY: aerocomp, sizeclass, nclass, &
                               nham_subm, HAM_M7, HAM_SALSA
  USE mo_ham_streams,    ONLY: nact_strat, nact_conv,                    &
                               nbcsol_ait, nbcsol_acc,                  & 
                               nbcsol_diag, nbcinsol_diag, nduinsol_diag,&
                               ndusol_strat, ndusol_cv,                  &
                               nbcsol_strat, nbcsol_cv,                  &
                               nduinsolai, nduinsolci, nbcinsol,         &
                               naerinsol
  USE mo_ham_salsa_trac, ONLY: idt_mdu,idt_mbc
  USE mo_ham_salsactl,   ONLY: in1a, in2b, fn2b, in2a, fn2a, idub, ibcb

  IMPLICIT none
 
  INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
  REAL(dp), INTENT(IN)  :: prho(kbdim,klev), pxtm1(kbdim,klev,ntrac)

  INTEGER :: naero_frz 

  !SF remark: Declan O'Donnell had initially removed the 'allocatable' prop for iaero_frz and znaero_frz,
  !           as this seemed to create pbs with OpenMP (in rev 410 of echam6-hammoz/branches/declan)
  !           Tommi Bergman needed to re-introduce the 'allocatable' prop, in order to handle
  !           both M7 and SALSA.
  !           I suspect that the original OpenMP pb was more related to the fact that these arrays
  !           were allocated at each timestep before Declan removed this, rather than to the 
  !           sole ability for them to be allocated at all.
  !           For this reason, I keep the 'allocatable' re-introduction, but I add some security to
  !           avoid multiple allocations.
  !           I have currently no way to test if this is acceptable under OpenMP conditions. If it is
  !           shown that the current implementation is still causing troubles with OMP, then one will
  !           have to set the 2 relevant arrays to the max size (salsa), which sounds a bit silly for now.

  INTEGER, ALLOCATABLE :: iaero_frz(:) 
  INTEGER              :: jf, jclass, jspec, jtn 
 
  REAL(dp) :: zdtime
  REAL(dp), ALLOCATABLE :: znaero_frz(:,:,:)
  REAL(dp)              :: zratio(kbdim,klev), znumb(kbdim,klev)

  INTEGER :: i,j
  LOGICAL :: llsol

  zdtime = delta_time

  !--- Set the needed species-modes in a compact array

  SELECT CASE(nham_subm)
      CASE(HAM_M7)

         naero_frz = 8

         IF (.NOT. ALLOCATED(iaero_frz))  ALLOCATE(iaero_frz(naero_frz))
         IF (.NOT. ALLOCATED(znaero_frz)) ALLOCATE(znaero_frz(kbdim,klev,naero_frz))
    
         iaero_frz(1) = speclist(id_bc)%iaerocomp(iaiti)    !! im7table(iaiti,id_bc)
         iaero_frz(2) = speclist(id_du)%iaerocomp(iacci)    !! im7table(iacci,id_du)
         iaero_frz(3) = speclist(id_du)%iaerocomp(icoai)    !! im7table(icoai,id_du)
         iaero_frz(4) = speclist(id_bc)%iaerocomp(iaits)    !! im7table(iaits,id_bc)
         iaero_frz(5) = speclist(id_bc)%iaerocomp(iaccs)    !! im7table(iaccs,id_bc)
         iaero_frz(6) = speclist(id_bc)%iaerocomp(icoas)    !! im7table(icoas,id_bc)
         iaero_frz(7) = speclist(id_du)%iaerocomp(iaccs)    !! im7table(iaccs,id_du)
         iaero_frz(8) = speclist(id_du)%iaerocomp(icoas)    !! im7table(icoas,id_du)

      CASE(HAM_SALSA) 
           
         naero_frz = count(ibcb(:)>0)+count(idub(:)>0)!n_insoluble+n_soluble

         IF (.NOT. ALLOCATED(iaero_frz))  ALLOCATE(iaero_frz(naero_frz))
         IF (.NOT. ALLOCATED(znaero_frz)) ALLOCATE(znaero_frz(kbdim,klev,naero_frz))

         j=0
         DO i =in1a,fn2b
            IF(ibcb(i)>0) THEN
               j=j+1
               iaero_frz(j) =ibcb(i)!speclist(id_bc)%iaerocomp(i) 
            END IF
         END DO
         DO i =in1a,fn2b
            IF(ibcb(i)>0) THEN
               j=j+1
               iaero_frz(j) =idub(i)!speclist(id_bc)%iaerocomp(i) 
            END IF
         END DO

  END SELECT
  
  !--- Loop over the relevant species-modes for freezing:
  DO jf=1,naero_frz
 
     jclass = aerocomp(iaero_frz(jf))%iclass !mode identity
     jspec  = aerocomp(iaero_frz(jf))%spid   !species identity
     jtn    = sizeclass(jclass)%idt_no   !mode number tracer identity
     llsol  = sizeclass(jclass)%lsoluble !soluble property
 
     IF (llsol) THEN !soluble mode, immersion freezing: 
                     !  need to use volume ratio
                     !  need to use 'the activated fraction of the total mode number concentration'
 
        CALL aero_massvolratio(kproma, kbdim, klev, jclass, jspec, 'volu', pxtm1,zratio)
        znumb(1:kproma,:)  = nact_strat(jclass)%ptr(1:kproma,:,krow)  !SF note: prho is already contained in it  
 
     ELSE            !insoluble mode, contact freezing:
                     !   need to use mass ratio
                     !   need to use the total mode number concentration
 
        CALL aero_massvolratio(kproma, kbdim, klev, jclass, jspec, 'mass', pxtm1, zratio)
        znumb(1:kproma,:)  = pxtm1(1:kproma,:,jtn) * prho(1:kproma,:) 
 
     ENDIF !SF llsol

     !>>dod correct array shape
     CALL aero_nc_surfw(kproma, kbdim, klev, zratio, znumb, znaero_frz(:,:,jf))
 
  ENDDO !end loop over species-modes for freezing
 
  !--- Remap number concentrations to the diag arrays:
  SELECT CASE(nham_subm)
      CASE (HAM_M7) 
         nbcinsol(1:kproma,:,krow)    = znaero_frz(1:kproma,:,1)
         nduinsolai(1:kproma,:,krow)  = znaero_frz(1:kproma,:,2)
         nduinsolci(1:kproma,:,krow)  = znaero_frz(1:kproma,:,3)
         
         nbcsol_strat(1:kproma,:,krow) = znaero_frz(1:kproma,:,5) & !DN #295 removed the contribution from the aitken mode
                                       + znaero_frz(1:kproma,:,6)
     
         ndusol_strat(1:kproma,:,krow) = znaero_frz(1:kproma,:,7) &
                                       + znaero_frz(1:kproma,:,8)
         
         !For diagnostics:
         nbcsol_ait(1:kproma,:,krow)    = nbcsol_ait(1:kproma,:,krow)  + znaero_frz(1:kproma,:,4)*zdtime
         !>>dod bugfix
         nbcsol_acc(1:kproma,:,krow)    = nbcsol_acc(1:kproma,:,krow)  + znaero_frz(1:kproma,:,5)*zdtime
         !<<dod

      CASE(HAM_SALSA)
         nbcinsol(1:kproma,:,krow)     = sum(znaero_frz(1:kproma,:,8:14),3)
         nduinsolai(1:kproma,:,krow)   = sum(znaero_frz(1:kproma,:,22:25),3)
         nduinsolci(1:kproma,:,krow)   = sum(znaero_frz(1:kproma,:,26:28),3)
         
         nbcsol_strat(1:kproma,:,krow) = sum(znaero_frz(1:kproma,:,4:7),3) !HK: Check if these are correct indices
         
         ndusol_strat(1:kproma,:,krow) = sum(znaero_frz(1:kproma,:,15:21),3) 
     
  END SELECT

  nbcsol_diag(1:kproma,:,krow)   = nbcsol_diag(1:kproma,:,krow)   + nbcsol_strat(1:kproma,:,krow)*zdtime
  nbcinsol_diag(1:kproma,:,krow) = nbcinsol_diag(1:kproma,:,krow) + nbcinsol(1:kproma,:,krow)*zdtime
  nduinsol_diag(1:kproma,:,krow) = nduinsol_diag(1:kproma,:,krow)                                     &
                                 + (nduinsolai(1:kproma,:,krow) + nduinsolci(1:kproma,:,krow))*zdtime
 
  !--- Get the total number of aerosols in insoluble modes:
  naerinsol(1:kproma,:,krow) = 0._dp
  DO jclass=1,nclass
     jtn = sizeclass(jclass)%idt_no !mode number tracer identity
     IF (.NOT. sizeclass(jclass)%lsoluble) THEN !insoluble modes
        naerinsol(1:kproma,:,krow) = naerinsol(1:kproma,:,krow) + pxtm1(1:kproma,:,jtn)
     ENDIF
  ENDDO
  naerinsol(1:kproma,:,krow) = naerinsol(1:kproma,:,krow)*prho(1:kproma,:)
 
  END SUBROUTINE get_aerofreez_nc
 
END MODULE mo_ham_freezing

