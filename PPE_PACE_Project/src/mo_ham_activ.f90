!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_activ.f90
!!
!! \brief
!! This code handles all ham-dependent activation processes
!! 
!!
!! \author Philip Stier (MPI-Met)
!! \author Sylvaine Ferrachat (ETHZ)
!!
!! \responsible_coder
!! Sylvaine Ferrachat, sylvaine.ferrachat@env.ethz.ch
!!
!! \revision_history
!!   -# P. Stier (MPI-Met)  - original code 
!!   -# S. Ferrachat (ETHZ) - new code structure (separate HAM-dependent pieces 
!!                            from the general activation schemes) - (2010-03)
!!
!! \limitations
!! None
!!
!! \details
!! Implementation of :
!!   - the Abdul-Razzak & Ghan scheme which calculates the number of activated aerosol 
!!     particles from the aerosol size-distribution, composition and ambient supersaturation 
!!     (see subroutine ham_activ_abdulrazzak_ghan)
!!   - a preparatory routine for Lin & Leaitch activation scheme (HAM-specific), which computes
!!     the fractional mass and number of each mode larger than the cutoff of the instrument and add them up
!!
!! \bibliographic_references
!!    - Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
!!    - Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
!!    - Pruppbacher and Klett, Kluewer Ac. Pub., 1997.
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

MODULE mo_ham_activ

  USE mo_kind,          ONLY: dp
  USE mo_ham,           ONLY: sizeclass,nclass

  IMPLICIT NONE

  PUBLIC ham_activ_abdulrazzak_ghan, ham_activ_diag_abdulrazzak_ghan_strat, ham_activ_diag_abdulrazzak_ghan_conv
  PUBLIC ham_avail_activ_lin_leaitch, ham_activ_diag_lin_leaitch
  PUBLIC ham_activ_koehler_ab

  PRIVATE

  !--- Subroutines:

CONTAINS

  SUBROUTINE ham_activ_abdulrazzak_ghan(kproma,   kbdim,   klev,  krow,  ktdia, &
                                        pcdncact, pesw,    prho,                &
                                        pxtm1,    ptm1,    papm1, pqm1,         &
                                        pw,       pwpdf,   pa,    pb,           &
                                        prdry,    pnact,   pfracn,              &
                                        psc,      prc,     psmax  )


    ! *ham_activ_abdulrazzak_ghan* calculates the number of activated aerosol 
    !              particles from the aerosol size-distribution,
    !              composition and ambient supersaturation
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-MET, Caltech, University of Oxford  2002-2009
    !
    ! Method:
    ! -------
    ! The calculation of the activation can be reduced to 4 tasks:
    ! 
    ! 0)   Calculation of Koehler A/B coefficients 
    !      (now done in ham_activ_koehler_ab)
    ! I)   Calculate the maximum supersaturation
    ! II)  Calculate the corresponding radius of activation
    !      for each mode
    ! III) Calculate the number of particles that are larger
    !      than the radius of activation for each mode.
    ! 
    ! III) Calculation of the number of activated particles:
    !      See the routine ham_m7_logtail below.
    !
    ! The calculations are now performed separately for 
    ! stratiform and convective updraft velocities.
    !
    ! References:
    ! -----------
    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997.

    !>>dod soa
    USE mo_ham_m7ctl,   ONLY: sigmaln, &
    !>>UP #810
                              l1ccn, n1ccn
    !<<UP
    !<<dod
    USE mo_ham_tools,   ONLY: ham_m7_logtail
    USE mo_tracdef,     ONLY: trlist, ntrac, AEROSOLMASS
    USE mo_ham,         ONLY: naerocomp, aerocomp, sizeclass, nclass, HAM_M7, naeroclass
    USE mo_math_constants,     ONLY: pi
    USE mo_physical_constants, ONLY: rhoh2o, argas, rv, cpd, grav, alv, amw, amd, tmelt
    USE mo_echam_cloud_params, ONLY: cthomi ! minimum T for mixed clouds
    !>>SF
    USE mo_conv,        ONLY: cdncact_cv
    !<<SF
    !>>DN #364
    USE mo_ham_streams,  ONLY: frac
    USE mo_activ,        ONLY: na, nw
    !<<DN #364

    IMPLICIT NONE

    !--- Arguments:

    INTEGER, INTENT(in)   :: kproma, kbdim, klev, krow, ktdia

    REAL(dp), INTENT(out) :: pcdncact(kbdim,klev),       & ! number of activated particles
                             pnact(kbdim,klev,nclass),   & ! number of activated particles per mode [m-3]
                             pfracn(kbdim,klev,nclass),  & ! fraction of activated particles per mode
                             psc(kbdim,klev,nclass),     & ! critical supersaturation [% 0-1]
                             prc(kbdim,klev,nclass,nw),  & ! critical radius of activation per mode [m]
                             psmax(kbdim,klev,nw)          ! maximum supersaturation [% 0-1]

    REAL(dp), INTENT(in)  :: ptm1(kbdim,klev),           & ! temperature
                             papm1(kbdim,klev),          & ! pressure 
                             prho(kbdim,klev),           & ! air density
                             pqm1(kbdim,klev),           & ! specific humidity
                             pesw(kbdim,klev),           & ! saturation water vapour pressure
                             pw(kbdim,klev,nw),          & ! mean or bins of updraft velocity (>0.0) [m s-1]
                             pwpdf(kbdim,klev,nw),       & ! PDF of updraft velocity [s m-1]
                             pxtm1(kbdim,klev,ntrac),    &
                             pa(kbdim,klev,nclass),      & ! curvature parameter A of the Koehler equation
                             pb(kbdim,klev,nclass),      & ! hygroscopicity parameter B of the Koehler equation
                             prdry(kbdim,klev,nclass)      ! dry radius for each mode

    !--- Local variables:

    INTEGER :: jclass,      jk,            &
               jt,          jl,            &
               jw

    REAL(dp):: zkv,                        &
               zeps,      zalpha,          &
               zgamma,                     &
               zgrowth,                    &
               zdif,      zxv,             &
               zk,        zka

    REAL(dp):: zamw,                       & ! molecular weight of water [kg mol-1]
               zamd                          ! molecular weight of dry air [kg mol-1]

    REAL(dp):: zf(nclass), zg(nclass)

    REAL(dp):: zfracn(kbdim,klev)            ! fraction of activated aerosol numbers for current mode and w bin

    REAL(dp):: zfracn_top(kbdim,klev,nclass), & ! w * dw weighted activated aerosol fraction
               zfracn_bot(kbdim,klev,nclass)    ! weighting factor  (sum of w * dw)

    REAL(dp):: zn(kbdim,klev,nclass),      & ! aerosol number concentration for each mode [m-3]
               zsm(kbdim,klev,nclass)        ! critical supersaturation for activating particles
                                             !    with the mode number median radius

    REAL(dp) :: zeta(nw), zxi(nw), zsum(nw)

    REAL(dp), PARAMETER :: zsten = 75.0E-3_dp ! surface tension of H2O [J m-2] 
                                              !   neglecting salts and temperature
                                               !   (also tried P&K 5.12 - erroneous!)

    LOGICAL, PARAMETER :: ll_numb = .TRUE.  ! switch between number/mass in logtail calculation !SF


    !--- 0) Initializations:

    zsm(1:kproma,:,:)        = 0._dp
    zfracn_top(1:kproma,:,:) = 0._dp
    zfracn_bot(1:kproma,:,:) = 0._dp
    prc(1:kproma,:,:,:)      = 1._dp ! [m] initialized with 1m, only changed if activation occurs
    pnact(1:kproma,:,:)      = 0._dp
    pfracn(1:kproma,:,:)     = 0._dp
    pcdncact(1:kproma,:)     = 0._dp

    zeps=EPSILON(1._dp)

    !--- Conversions to SI units [g mol-1 to kg mol-1]:
    
    zamw=amw*1.E-3_dp
    zamd=amd*1.E-3_dp

    !--- Number per unit volume for each mode:
    DO jclass=1, nclass
       jt = sizeclass(jclass)%idt_no
       !>>dod #377
       IF (sizeclass(jclass)%lactivation) THEN
          zn(1:kproma,:,jclass)=pxtm1(1:kproma,:,jt)*prho(1:kproma,:)
       END IF
       !<<dod
    END DO

    ! (7):
    zf(:)=0.5_dp*EXP(2.5_dp*sigmaln(:)**2._dp)

    ! (8):
    zg(:)=1._dp+0.25_dp*sigmaln(:)

    !--- 1) Calculation of Koehler A/B coefficients: 
    !       Now done in ham_activ_koehler_ab once so that they can be used in
    !       convective and stratiform activation 

    !--- 2) Calculate maximum supersaturation:

    !--- 2.1) Abbdul-Razzak and Ghan (2000):
    !         (Equations numbers from this paper unless otherwise quoted)

    DO jk=ktdia, klev
       DO jl=1, kproma

          !--- Water vapour pressure:

          IF( (nw>1 .OR. pw(jl,jk,1)>zeps) .AND. &
              pqm1(jl,jk)>zeps             .AND. &
              ptm1(jl,jk)>cthomi                 ) THEN  

             !--- Abdul-Razzak et al. (1998) (Eq. 11):

             zalpha=(grav*zamw*alv)/(cpd*argas*ptm1(jl,jk)**2) - &
                    (grav*zamd)/(argas*ptm1(jl,jk))

             zgamma=(argas*ptm1(jl,jk))/(pesw(jl,jk)*zamw) +  &
                    (zamw*alv**2)/(cpd*papm1(jl,jk)*zamd*ptm1(jl,jk))

             !--- Diffusivity of water vapour in air (P&K, 13.3) [m2 s-1]:

             zdif=0.211_dp * (ptm1(jl,jk)/tmelt)**1.94_dp * (101325._dp/papm1(jl,jk)) *1.E-4_dp

             !--- Thermal conductivity zk (P&K, 13.18) [cal cm-1 s-1 K-1]:

             ! Mole fraction of water:

             zxv=pqm1(jl,jk)*(zamd/zamw)

             zka=(5.69_dp+0.017_dp*(ptm1(jl,jk)-273.15_dp))*1.E-5_dp

             zkv=(3.78_dp+0.020_dp*(ptm1(jl,jk)-273.15_dp))*1.E-5_dp

             ! Moist air, convert to [J m-1 s-1 K-1]:

             zk =zka*(1._dp-(1.17_dp-1.02_dp*zkv/zka)*zxv) * 4.1868_dp*1.E2_dp

             !--- Abdul-Razzak et al. (1998) (Eq. 16):

             zgrowth=1._dp/                                                   &
                       ( (rhoh2o*argas*ptm1(jl,jk))/(pesw(jl,jk)*zdif*zamw) + &
                       (alv*rhoh2o)/(zk*ptm1(jl,jk)) * ((alv*zamw)/(ptm1(jl,jk)*argas) -1._dp) )

             !--- Summation for equation (6):

             zsum(:)=0._dp

             DO jclass=1, nclass
                !>>dod #377
                IF (sizeclass(jclass)%lactivation) THEN
                   IF (zn(jl,jk,jclass)    > zeps     .AND. &
                        prdry(jl,jk,jclass) > 1.E-9_dp .AND. &
                        pa(jl,jk,jclass)    > zeps     .AND. &
                        pb(jl,jk,jclass)    > zeps           ) THEN

                      ! (9):

                      zsm(jl,jk,jclass)=2._dp/SQRT(pb(jl,jk,jclass)) * &
                                   (pa(jl,jk,jclass)/(3._dp*prdry(jl,jk,jclass)))**1.5_dp

                      ! (10):

                      zxi(:)=2._dp*pa(jl,jk,jclass)/3._dp * SQRT(zalpha*pw(jl,jk,:)/zgrowth)

                      ! (11):
                      
                      zeta(:)=((zalpha*pw(jl,jk,:)/zgrowth)**1.5_dp) / &
                             (2._dp*pi*rhoh2o*zgamma*zn(jl,jk,jclass))
                      
                      ! (6):

                      WHERE (pw(jl,jk,:)>zeps)
                         zsum(:)=zsum(:) + ( 1._dp/zsm(jl,jk,jclass)**2                  &
                                          * ( zf(jclass)*(zxi(:)/zeta(:))**1.5_dp     &
                                              + zg(jclass)*( zsm(jl,jk,jclass)**2._dp &
                                                             / (zeta(:)+3._dp*zxi(:)) )**0.75_dp ) )
                      END WHERE

                   ENDIF
                END IF
                !<<dod
             END DO ! jclass

             WHERE (zsum(:) > zeps)
                psmax(jl,jk,:)=1._dp/SQRT(zsum(:))
             ELSEWHERE
                psmax(jl,jk,:)=0._dp
             END WHERE

          ELSE
             psmax(jl,jk,:)=0._dp
          END IF

       END DO ! jl
    END DO ! jk
    
    !--- Diagnostics:
    DO jclass=1, nclass
       !>>dod #377
       IF (sizeclass(jclass)%lactivation) THEN
          psc(1:kproma,:,jclass) = zsm(1:kproma,:,jclass)
       END IF
       !<<dod
    END DO

    !--- 3) Calculate activation:

    DO jw=1, nw
       DO jclass=1, nclass

          IF (sizeclass(jclass)%lactivation) THEN !>>dod<< #377

             WHERE (psmax(1:kproma,ktdia:klev,jw)>zeps      .AND. &
                    zsm(1:kproma,ktdia:klev,jclass)>zeps    .AND. &
                    zn(1:kproma,ktdia:klev,jclass)>zeps     .AND. &
                    prdry(1:kproma,ktdia:klev,jclass)>1.E-9_dp       )

                prc(1:kproma,ktdia:klev,jclass,jw)                            &
                    = prdry(1:kproma,ktdia:klev,jclass)                       &
                        * ( zsm(1:kproma,ktdia:klev,jclass)                   &
                            / psmax(1:kproma,ktdia:klev,jw) )**(2._dp/3._dp)

             END WHERE

             !--- 3.2) Calculate the fractional number of each mode
             !         larger than the mode critical radius:
             CALL ham_m7_logtail(kproma,    kbdim,  klev,   krow, jclass, &
                                 ll_numb,   prdry(:,:,jclass),            &
                                 prc(:,:,jclass,jw), zfracn(:,:))

             !--- 3.3) Sum up the total number of activated particles, integrating over updraft PDF [m-3]:
             ! The weighting here should be correct (up to discretisation error)
             ! provided that the w bins are equally spaced and cover the whole
             ! range of integration.
             zfracn_top(1:kproma,:,jclass) = zfracn_top(1:kproma,:,jclass) &
                                           + zfracn(1:kproma,:)*pwpdf(1:kproma,:,jw)
             zfracn_bot(1:kproma,:,jclass) = zfracn_bot(1:kproma,:,jclass) &
                                           + pwpdf(1:kproma,:,jw)

          END IF ! lactivation !>>dod<<

       END DO ! jclass
    END DO ! jw

    DO jclass=1, nclass
       IF (sizeclass(jclass)%lactivation) THEN !>>dod<< #377
          !>>UP #810
          IF (l1ccn) THEN
                  IF (jclass == n1ccn) THEN
                          pfracn(1:kproma,:,jclass) = zfracn_top(1:kproma,:,jclass)/zfracn_bot(1:kproma,:,jclass)
                  ELSE
                          pfracn(1:kproma,:,jclass) = 0.0_dp
                  END IF
          ELSE
                  pfracn(1:kproma,:,jclass) = zfracn_top(1:kproma,:,jclass)/zfracn_bot(1:kproma,:,jclass)
          END IF
          !<<UP #810
          pnact(1:kproma,:,jclass)  = pfracn(1:kproma,:,jclass) * zn(1:kproma,:,jclass)
          pcdncact(1:kproma,:)      = pcdncact(1:kproma,:) + pnact(1:kproma,:,jclass)
       END IF
    END DO

  END SUBROUTINE ham_activ_abdulrazzak_ghan

  SUBROUTINE ham_activ_diag_abdulrazzak_ghan_strat(kproma, kbdim, klev,       &
                                                   krow,   pnact, pfracn,     &
                                                   prc,    psmax )

    USE mo_activ,       ONLY: na, nw, swat_max_strat
    USE mo_ham_streams, ONLY: frac, nact_strat, rc_strat
    USE mo_param_switches, ONLY: nactivpdf

    INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(IN) :: pnact(kbdim,klev,nclass), & ! number of activated particles per mode [m-3]
                            pfracn(kbdim,klev,nclass),& ! fraction of activated particles per mode
                            prc(kbdim,klev,nclass,nw),& ! critical radius of activation per mode and w bin [m]
                            psmax(kbdim,klev,nw)        ! maximum supersaturation per w bin [% 0-1]

    INTEGER :: jclass, jw

    IF (nactivpdf <= 0) THEN
      DO jw=1, nw
        swat_max_strat(jw)%ptr(1:kproma,:,krow) = psmax(1:kproma,:,jw)

        DO jclass=1, nclass
           !>>dod #377
           IF (sizeclass(jclass)%lactivation) THEN
              WHERE(prc(1:kproma,:,jclass,jw) /= 1._dp)
                 rc_strat(jclass,jw)%ptr(1:kproma,:,krow)=prc(1:kproma,:,jclass,jw)
              ENDWHERE
           END IF
           !>>dod
        END DO
      END DO
    END IF

    na(1:kproma,:,krow) = 0._dp

    DO jclass=1, nclass
       !>>dod #377
       IF (sizeclass(jclass)%lactivation) THEN
          nact_strat(jclass)%ptr(1:kproma,:,krow)=pnact(1:kproma,:,jclass)
          frac(jclass)%ptr(1:kproma,:,krow)=pfracn(1:kproma,:,jclass)
          na(1:kproma,:,krow)=na(1:kproma,:,krow)+pnact(1:kproma,:,jclass)
       END IF
       !>>dod
    END DO

  END SUBROUTINE ham_activ_diag_abdulrazzak_ghan_strat

  SUBROUTINE ham_activ_diag_abdulrazzak_ghan_conv(kproma, kbdim, klev,       &
                                                  krow,   pnact, prc,  psmax )

    USE mo_activ,       ONLY: nw, swat_max_conv
    USE mo_ham_streams, ONLY: nact_conv, rc_conv
    USE mo_param_switches, ONLY: nactivpdf

    INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(IN) :: pnact(kbdim,klev,nclass), & ! number of activated particles per mode [m-3]
                            prc(kbdim,klev,nclass,nw),& ! critical radius of activation per mode and w bin [m]
                            psmax(kbdim,klev,nw)        ! maximum supersaturation per w bin [% 0-1]

    INTEGER :: jclass, jw

    IF (nactivpdf <= 0) THEN
      DO jw=1,nw
         swat_max_conv(jw)%ptr(1:kproma,:,krow) = psmax(1:kproma,:,jw)
         DO jclass=1, nclass
            !>>dod #377
            IF (sizeclass(jclass)%lactivation) THEN
               WHERE(prc(1:kproma,:,jclass,jw) /= 1._dp)
                  rc_conv(jclass,jw)%ptr(1:kproma,:,krow)=prc(1:kproma,:,jclass,jw)
               ENDWHERE
            END IF
            !>>dod
         END DO
      END DO
    END IF

    DO jclass=1, nclass
       !>>dod #377
       IF (sizeclass(jclass)%lactivation) THEN
          nact_conv(jclass)%ptr(1:kproma,:,krow)=pnact(1:kproma,:,jclass)
       END IF
       !>>dod
    END DO

  END SUBROUTINE ham_activ_diag_abdulrazzak_ghan_conv


  SUBROUTINE ham_activ_koehler_ab(kproma,   kbdim,   klev,  krow,  ktdia, &
                                  pxtm1,    ptm1,  pa,    pb              )

    ! *ham_activ_koehler_ab* calculates the Koehler A and B coefficients
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford, 2013
    !
    ! References:
    ! -----------
    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997.

    USE mo_physical_constants,   ONLY: rhoh2o, argas, amw
    USE mo_ham,         ONLY: naerocomp, aerocomp, sizeclass
    USE mo_tracdef,     ONLY: ntrac
!davidn
    USE mo_param_switches,     ONLY: lepsfix
    USE mo_ham_streams,        ONLY: densaer
!davidn

    IMPLICIT NONE

    !--- Arguments:

    INTEGER :: kproma, kbdim, klev, krow, ktdia

    REAL(dp), INTENT(IN)  :: ptm1(kbdim,klev),        & ! temperature
                             pxtm1(kbdim,klev,ntrac)

    REAL(dp), INTENT(OUT) :: pa(kbdim,klev,nclass),   & ! curvature parameter A of the Koehler equation
                             pb(kbdim,klev,nclass)      ! hygroscopicity parameter B of the Koehler equation

    !--- Local variables:

    INTEGER :: jclass, jt, jl, jk

    REAL(dp):: zmoleweight,                &
               znion,     zosm,            &
               zrhoaer,   zeps

    REAL(dp):: zamw                          ! molecular weight of water [kg mol-1]

    REAL(dp):: zmassfrac(kbdim,klev)

    REAL(dp):: zsumtop(kbdim,klev,nclass), & ! temporary summation field
               zsumbot(kbdim,klev,nclass)    ! temporary summation field

    REAL(dp):: zmasssum(kbdim,klev,nclass)

    REAL(dp), PARAMETER :: zsten = 75.0E-3_dp ! surface tension of H2O [J m-2] 

    REAL(dp) :: zafac                        ! for calculation of curvature term in Koehler equations

    INTEGER :: jn, jspec
    INTEGER :: nion
!davidn
    REAL(dp),            POINTER :: densaer_p(:,:,:)
!davidn

    zmasssum(1:kproma,:,:) = 0._dp
    zmassfrac(1:kproma,:)  = 0._dp !davidn
    pa(1:kproma,:,:)       = 0._dp
    pb(1:kproma,:,:)       = 0._dp

    zsumtop(1:kproma,:,:)  = 0._dp
    zsumbot(1:kproma,:,:)  = 0._dp

    zeps=EPSILON(1._dp)

    !--- Conversions to SI units [g mol-1 to kg mol-1]:
    
    zamw=amw*1.E-3_dp

    !---for calculation of curvature parameter A (equation (5) in first paper)
    zafac = 2._dp*zsten*zamw / (rhoh2o*argas)

    !--- Sum mass mixing ratio for each mode:
    DO jn = 1,naerocomp
       jclass = aerocomp(jn)%iclass
       jspec = aerocomp(jn)%spid
       jt = aerocomp(jn)%idt
       !>>dod #377
       IF (sizeclass(jclass)%lactivation) THEN
          zmasssum(1:kproma,:,jclass)=zmasssum(1:kproma,:,jclass)+pxtm1(1:kproma,:,jt)
       END IF
       !<<dod
    END DO

    !--- 1) Calculate properties for each aerosol mode:

    !--- 1.1) Calculate the auxiliary parameters A and B of the Koehler equation:

       !--- 1) Calculate weighted properties:
       !       (Abdul-Razzak & Ghan, 2000)

    DO jn=1,naerocomp
       zmoleweight = aerocomp(jn)%species%moleweight*1.E-3_dp         ! [kg mol-1]
       nion = aerocomp(jn)%species%nion
       znion = REAL(nion,dp)
       zosm = aerocomp(jn)%species%osm
       zrhoaer = aerocomp(jn)%species%density
       jt = aerocomp(jn)%idt
       jclass = aerocomp(jn)%iclass
       
!davidn 
       densaer_p => densaer(jclass)%ptr
!davidn 
       IF (nion > 0 .AND. sizeclass(jclass)%lactivation) THEN      !>>dod<< #377
!davidn
          IF(lepsfix)THEN !#680
             WHERE(zmasssum(1:kproma,:,jclass)>zeps)           
                zmassfrac(1:kproma,:)=zmassfrac(1:kproma,:)+pxtm1(1:kproma,:,jt)/zmasssum(1:kproma,:,jclass)
             
                zsumtop(1:kproma,:,jclass)=zsumtop(1:kproma,:,jclass)+pxtm1(1:kproma,:,jt)*znion*zosm/zmoleweight
                zsumbot(1:kproma,:,jclass)=zsumbot(1:kproma,:,jclass)+pxtm1(1:kproma,:,jt)/densaer_p(1:kproma,:,krow)                   
             END WHERE
          ELSE
!davidn
          WHERE(zmasssum(1:kproma,:,jclass)>zeps)

             zmassfrac(1:kproma,:)=pxtm1(1:kproma,:,jt)/zmasssum(1:kproma,:,jclass)
             
             zsumtop(1:kproma,:,jclass)=zsumtop(1:kproma,:,jclass)+pxtm1(1:kproma,:,jt)*znion*zosm*zmassfrac(1:kproma,:)/zmoleweight
             zsumbot(1:kproma,:,jclass)=zsumbot(1:kproma,:,jclass)+pxtm1(1:kproma,:,jt)/zrhoaer
                   
          END WHERE
          !davidn
          END IF
          !davidn

       END IF! nion>0

    END DO !naerocomp

    DO jclass=1,nclass
       !>>dod #377
       IF(sizeclass(jclass)%lactivation) THEN
!davidn
         IF(lepsfix)THEN
          WHERE (zsumbot(1:kproma,:,jclass)>zeps)

             !--- 1.1.1) Hygroscopicity parameter B (Eq. 4) [1]:

             pb(1:kproma,:,jclass)=(zamw*zsumtop(1:kproma,:,jclass))/(rhoh2o*zsumbot(1:kproma,:,jclass))&
                  *zmassfrac(1:kproma,:)

             !--- 1.1.2) Calculate the curvature parameter A [m]:

             pa(1:kproma,:,jclass)= zafac/ptm1(1:kproma,:)

          END WHERE
          ELSE
!davidn
          WHERE (zsumbot(1:kproma,:,jclass)>zeps)

             !--- 1.1.1) Hygroscopicity parameter B (Eq. 4) [1]:

             pb(1:kproma,:,jclass)=(zamw*zsumtop(1:kproma,:,jclass))/(rhoh2o*zsumbot(1:kproma,:,jclass))

             !--- 1.1.2) Calculate the curvature parameter A [m]:

             pa(1:kproma,:,jclass)= zafac/ptm1(1:kproma,:)

          END WHERE
!davidn
          END IF
!davidn
       END IF
       !<<dod
    END DO !jclass=1, nclass

  END SUBROUTINE ham_activ_koehler_ab
    

  !---------------------------------------------------------------------------
  !>
  !! @brief Computes available particules for activation
  !! 
  !! @remarks Preparatory routine for Lin&Leaitch activation scheme (HAM-specific)
  !! Basically it computes the fractional mass and number of each mode
  !! larger than the cutoff of the instrument and adds them up
  !! Derived from the former aero_activ_lin_leaitch subroutine

  SUBROUTINE ham_avail_activ_lin_leaitch(kproma, kbdim, klev, krow, &
                                       prho, pxtm1)

    USE mo_activ,        ONLY: na
    USE mo_conv,         ONLY: na_cv
    USE mo_ham_tools,    ONLY: ham_m7_logtail
    USE mo_tracdef,      ONLY: ntrac
    USE mo_ham_streams,  ONLY: frac, rwet

    !SF note: lwetrad is true in this routine, ie wet radius is used
    USE mo_ham,          ONLY: nham_subm, HAM_M7, HAM_SALSA, sizeclass !>>dod<< #377

    INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(IN) :: prho(kbdim,klev)        ! air density
    REAL(dp), INTENT(IN) :: pxtm1(kbdim,klev,ntrac) ! tracer mmr

    REAL(dp), PARAMETER :: crcut=0.03*1E-6_dp ! Assumed lower cut-off of the
                                              ! aerosol size distribution [m]

    !--- Ulrike: included for activation in convective clouds
    REAL(dp), PARAMETER :: crcut_cv=0.02*1E-6_dp ! Assumed lower cut-off of the
                                                 ! aerosol size distribution in convective clouds [m]

    REAL(dp) :: cfracn(nclass)

    LOGICAL, PARAMETER :: lcut=.TRUE. ! explicit calculation of cut-off crcut or
                                      ! usage of prescribed values cfracn

    LOGICAL, PARAMETER :: ll_numb = .TRUE. ! switch between number/mass in logtail calculation !SF

    INTEGER  :: jclass, it

    REAL(dp) :: zr(kbdim,klev,nclass)
    REAL(dp) :: zfracn(kbdim,klev,nclass)
    REAL(dp) :: zfracn_cv(kbdim,klev,nclass)

    SELECT CASE(nham_subm)
        CASE(HAM_M7)
            cfracn(:) = (/1.0_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp/)
        CASE(HAM_SALSA)
            cfracn(:) = 1.0_dp
    END SELECT
    !>>dod bugfix
    na(1:kproma,:,krow) = 0._dp
    na_cv(1:kproma,:,krow) = 0._dp
    !<<dod

    !>>dod redmine #377: calculations and diagnotics restricted to relevant modes only. 
    !---  Calculate the fractional number of each mode
    !     larger than the cutoff of the instrument:
    SELECT CASE(nham_subm)
        CASE(HAM_M7)
            IF (lcut) THEN
     
               DO jclass=1, nclass !SF the ham_m7_logtail calculation is now done mode per mode for better efficiency
                  IF (sizeclass(jclass)%lactivation) THEN
                     IF (cfracn(jclass) > EPSILON(1._dp)) THEN !SF #279: only performs this calculation when relevant
                        !SF stratiform:
                        zr(1:kproma,:,jclass)=crcut
        
                        CALL ham_m7_logtail(kproma, kbdim,  klev,  krow, jclass, &
                                         ll_numb, rwet(jclass)%ptr(:,:,krow), &
                                         zr(:,:,jclass), zfracn(:,:,jclass) )
        
                        !SF convective:
                        zr(1:kproma,:,jclass)=crcut_cv
        
                        CALL ham_m7_logtail(kproma, kbdim,  klev,  krow, jclass, &
                                         ll_numb, rwet(jclass)%ptr(:,:,krow), &
                                         zr(:,:,jclass), zfracn_cv(:,:,jclass) )
                     ELSE !SF #279: ensures that zfracn and zfracn_cv are set to 0. when size class 
                          !         is not relevant for activation
                        zfracn(1:kproma,:,jclass)    = 0._dp
                        zfracn_cv(1:kproma,:,jclass) = 0._dp
                     ENDIF
                  END IF
               END DO

            ELSE
     
               DO jclass=1, nclass
                  IF (sizeclass(jclass)%lactivation) THEN
                     zfracn(1:kproma,:,jclass)    = cfracn(jclass)
                     zfracn_cv(1:kproma,:,jclass) = cfracn(jclass) !SF was missing in previous version. Mistake??
                  END IF
               END DO
     
            END IF

        CASE(HAM_SALSA)
            DO jclass=1, nclass
               zfracn(1:kproma,:,jclass)    = cfracn(jclass)
               zfracn_cv(1:kproma,:,jclass) = cfracn(jclass) !SF was missing in previous version. Mistake??
            END DO
            ! If lcut for SALSA cut out the particles in 1a r=~25nm
            ! Approximately the same as with M7 having crcut=30nm
            IF (lcut) THEN
               zfracn(1:kproma,:,1:3)=0.0_dp
               ! For convective clouds.. cut out 1a r=~25nm, (M7 crcut=20nm)
               zfracn_cv(1:kproma,:,1:3)=0.0_dp
            END IF
    END SELECT
    !--- Sum up aerosol number concentrations and convert from [kg-1] to [m-3]:
    DO jclass=1, nclass
       !>>dod #377
       IF(sizeclass(jclass)%lactivation) THEN

          !>>dod soa
          it = sizeclass(jclass)%idt_no
          !<<dod
          na(1:kproma,:,krow) = na(1:kproma,:,krow)                   &
                              + pxtm1(1:kproma,:,it)*prho(1:kproma,:) &
                              *zfracn(1:kproma,:,jclass)*cfracn(jclass)

          !--- Ulrike: included for NA from convection ---
          na_cv(1:kproma,:,krow) = na_cv(1:kproma,:,krow)                   &
                                 + pxtm1(1:kproma,:,it)*prho(1:kproma,:)    &
                                 *zfracn_cv(1:kproma,:,jclass)*cfracn(jclass)
          !--- end included

          frac(jclass)%ptr(1:kproma,:,krow)=zfracn(1:kproma,:,jclass)
          
       END IF
       !<<dod
    END DO

  END SUBROUTINE ham_avail_activ_lin_leaitch

  SUBROUTINE ham_activ_diag_lin_leaitch(kproma, kbdim, klev, krow, prho, pxtm1, pcdncact, pcdnc_min, &
                                        ll_actccn, ll_actccn_wobase, ll_bas, zrwet)

    USE mo_activ,        ONLY: na, &
                        !>>UP
                        ccn_strat_cdncmin, ccn_strat_cdncmin_ctr, &
                        ccn_strat_actccn, ccn_strat_actccn_ctr, &
                        ccn_strat_noactccn, ccn_strat_noactccn_ctr, &
                        ccn_strat_actccnwobase, ccn_strat_actccnwobase_ctr, &
                        ccn_strat_nosmall40, ccn_strat_nosmall40_ctr, &
                        ccn_strat_nosmall1, ccn_strat_nosmall1_ctr, &
                        ccn_strat_so4ks, &
                        ccn_strat_max, ccn_strat_max_ctr
    USE mo_ham_m7_trac, ONLY: idt_ms4ks, idt_nks, idt_mocks, idt_mbcks, &
                              idt_ms4as, idt_nas, idt_mocas, idt_mbcas, &
                              idt_mduas, idt_mssas
    USE mo_ham_m7ctl,   ONLY: iaits, iaccs
    USE mo_activ,       ONLY: rwet_KS_diag, rwet_KS_diag_cloudbase, &
                              rwet_KS_diag_cloudbase_ctr
    USE mo_param_switches, ONLY: lccnclimdiags
    !<<UP
    USE mo_conv,         ONLY: cdncact_cv, na_cv
    USE mo_ham_streams,  ONLY: frac, nact_strat, nact_conv, &
                               ccn_strat, ccn_conv !SF #333
    USE mo_tracdef,      ONLY: ntrac
    !>>SF #333
    USE mo_ham,          ONLY: lccnclim_diag
    USE mo_time_control, ONLY: delta_time
    !<<SF #333
    USE mo_ham,         ONLY: sizeclass !>>dod<< #377

    INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(IN) :: prho(kbdim,klev)        ! air density
    REAL(dp), INTENT(IN) :: pxtm1(kbdim,klev,ntrac) ! tracer mmr
    REAL(dp), INTENT(IN) :: pcdncact(kbdim,klev) ! number of activated particles
    REAL(dp), INTENT(IN) :: pcdnc_min(kbdim,klev) ! minimal CDNC
    LOGICAL, INTENT(IN) :: ll_actccn(kbdim,klev) ! whether these CCN are
    ! potentially used for CD activation
    LOGICAL, INTENT(IN)  :: ll_actccn_wobase(kbdim,klev) ! logical for whether
    ! these CCN will be potentially used for CD activation, not including the cloud base check
    !>>UP
    INTEGER, INTENT(IN)  :: ll_bas(kbdim,klev) ! logical for cloud base
    REAL(dp), INTENT(IN) :: zrwet(kbdim,klev,krow) ! wet number mean radius, Aitken soluble
    !<<UP

    INTEGER  :: jclass, it, jccn, jw
    LOGICAL  :: ll1(kbdim,klev)
    REAL(dp) :: zeps
    REAL(dp) :: ztmp1(kbdim,klev), ztmp2(kbdim,klev)

    zeps = EPSILON(1._dp)

    DO jclass=1, nclass
       !>>dod #377
       IF(sizeclass(jclass)%lactivation) THEN
          !>>dod soa
          it = sizeclass(jclass)%idt_no
          !<<dod
          !>>SF #458 (replacing where statements)
          ll1(1:kproma,:) = (na(1:kproma,:,krow) > zeps)
          ztmp1(1:kproma,:) = MERGE(na(1:kproma,:,krow), 1._dp, ll1(1:kproma,:)) !SF 1._dp is a dummy val.
          ztmp2(1:kproma,:) = pxtm1(1:kproma,:,it) * prho(1:kproma,:) &
                            * frac(jclass)%ptr(1:kproma,:,krow) / ztmp1(1:kproma,:)
                            !UP note: fraction of particles in that mode (pxtm1) * fraction of activated particles in that mode (frac) / na
                            !-> activated particles in that mode per total
                            !activated particles
                            !UP note: I don't get the units, though, because
                            !pxtm1 is mass...

          nact_strat(jclass)%ptr(1:kproma,:,krow) = MERGE( &
                                                          pcdncact(1:kproma,:) * ztmp2(1:kproma,:), &
                                                          !UP note: * CDNC
                                                          !If the two scale the
                                                          !same, this gives the
                                                          !CDNC per mode
                                                          nact_strat(jclass)%ptr(1:kproma,:,krow), &
                                                          ll1(1:kproma,:))

          nact_conv(jclass)%ptr(1:kproma,:,krow) = MERGE( &
                                                        cdncact_cv(1:kproma,:,krow) * ztmp2(1:kproma,:), &
                                                        nact_conv(jclass)%ptr(1:kproma,:,krow), &
                                                        ll1(1:kproma,:))
          !<<SF #458 (replacing where statements)
       END IF
       !<<dod
    END DO

    !>>SF #333
    IF (lccnclim_diag) THEN
       ccn_strat(1:kproma,:,krow) = ccn_strat(1:kproma,:,krow) + na(1:kproma,:,krow)    * delta_time
       ccn_conv(1:kproma,:,krow)  = ccn_conv(1:kproma,:,krow)  + na_cv(1:kproma,:,krow) * delta_time
       !>>UP
       ccn_strat_cdncmin(1:kproma,:,krow) = ccn_strat_cdncmin(1:kproma,:,krow) + & 
                                            MERGE(na(1:kproma,:,krow) * delta_time, 0._dp, na(1:kproma,:,krow) > pcdnc_min)
       ccn_strat_cdncmin_ctr(1:kproma,:,krow) = ccn_strat_cdncmin_ctr(1:kproma,:,krow) + &
                                            MERGE(1._dp*delta_time, 0._dp, na(1:kproma,:,krow) > pcdnc_min)
       ccn_strat_actccn(1:kproma,:,krow) = ccn_strat_actccn(1:kproma,:,krow) + & 
                                            MERGE(na(1:kproma,:,krow) * delta_time, 0._dp, ll_actccn(1:kproma,:))
       ccn_strat_actccn_ctr(1:kproma,:,krow) = ccn_strat_actccn_ctr(1:kproma,:,krow) + &
                                            MERGE(1._dp*delta_time, 0._dp, ll_actccn(1:kproma,:))
       ccn_strat_noactccn(1:kproma,:,krow) = ccn_strat_noactccn(1:kproma,:,krow) + & 
                                            MERGE(0._dp, na(1:kproma,:,krow) * delta_time, ll_actccn(1:kproma,:))
       ccn_strat_noactccn_ctr(1:kproma,:,krow) = ccn_strat_noactccn_ctr(1:kproma,:,krow) + &
                                            MERGE(0._dp, 1._dp*delta_time, ll_actccn(1:kproma,:))
       ccn_strat_actccnwobase(1:kproma,:,krow) = ccn_strat_actccnwobase(1:kproma,:,krow) + & 
                                            MERGE(na(1:kproma,:,krow) * delta_time, 0._dp, ll_actccn_wobase(1:kproma,:))
       ccn_strat_actccnwobase_ctr(1:kproma,:,krow) = ccn_strat_actccnwobase_ctr(1:kproma,:,krow) + &
                                            MERGE(1._dp*delta_time, 0._dp, ll_actccn_wobase(1:kproma,:))
       ccn_strat_nosmall40(1:kproma,:,krow) = ccn_strat_nosmall40(1:kproma,:,krow) + &
                                            MERGE(na(1:kproma,:,krow) * delta_time, 0._dp, na(1:kproma,:,krow) > 40.0E6_dp)
       ccn_strat_nosmall40_ctr(1:kproma,:,krow) = ccn_strat_nosmall40_ctr(1:kproma,:,krow) + &
                                            MERGE(1._dp*delta_time, 0._dp, na(1:kproma,:,krow) > 40.0E6_dp)
       ccn_strat_nosmall1(1:kproma,:,krow) = ccn_strat_nosmall1(1:kproma,:,krow) + &
                                            MERGE(na(1:kproma,:,krow) * delta_time, 0._dp, na(1:kproma,:,krow) > 1.0E6_dp)
       ccn_strat_nosmall1_ctr(1:kproma,:,krow) = ccn_strat_nosmall1_ctr(1:kproma,:,krow) + &
                                            MERGE(1._dp*delta_time, 0._dp, na(1:kproma,:,krow) > 1.0E6_dp)
       ccn_strat_so4ks(1:kproma,:,krow) = ccn_strat_so4ks(1:kproma,:,krow) + &
                                          MERGE( &
                                          ! Mass of SO4 in KS mode 
                                          pxtm1(1:kproma,:,idt_ms4ks) / &
                                          ! Mass of the KS mode 
                                          (pxtm1(1:kproma,:,idt_ms4ks) + pxtm1(1:kproma,:,idt_mbcks) + pxtm1(1:kproma,:,idt_mocks)) &
                                          ! Number concentration in the KS mode 
                                          * pxtm1(1:kproma,:,idt_nks) &
                                          ! CCN fraction of the KS mode 
                                          * frac(iaits)%ptr(1:kproma,:,krow) &
                                          ! Conversion and time
                                          * prho(1:kproma,:) &
                                          * delta_time, &
                                          ! MERGE 
                                          0._dp, &
                                          ! MERGE 
                                          (pxtm1(1:kproma,:,idt_ms4ks) + pxtm1(1:kproma,:,idt_mbcks) + &
                                          pxtm1(1:kproma,:,idt_mocks)) > 0._dp)
       ccn_strat_max(1:kproma,:,krow) = MAX(na(1:kproma,:,krow), ccn_strat_max(1:kproma,:,krow))
       ccn_strat_max_ctr(1:kproma,:,krow) = 1._dp
                                        

       !<<UP
    ENDIF
    !<<SF #333

    !>>UP
    IF (lccnclimdiags) THEN
    ! Diagnose wet radius at cloud base and outside
    rwet_KS_diag(1:kproma,:,krow) = rwet_KS_diag(1:kproma,:,krow) + zrwet(1:kproma,:,krow) * delta_time
    rwet_KS_diag_cloudbase(1:kproma,:,krow) = rwet_KS_diag_cloudbase(1:kproma,:,krow) + &
                MERGE(zrwet(1:kproma,:,krow) * delta_time, 0._dp, ll_bas(1:kproma,:) > 0._dp)
    rwet_KS_diag_cloudbase_ctr(1:kproma,:,krow) = rwet_KS_diag_cloudbase_ctr(1:kproma,:,krow) + &
                MERGE(1._dp * delta_time, 0._dp, ll_bas(1:kproma,:) > 0._dp)
    ENDIF
    !<<UP
 
  END SUBROUTINE ham_activ_diag_lin_leaitch

END MODULE mo_ham_activ
