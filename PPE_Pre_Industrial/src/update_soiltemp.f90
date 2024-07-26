!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE update_soiltemp (klon, nsoil, cdel, cmid          &
                         , pts, psn                          &
                         , psodif, prgcgn                    &
                         , pgrndc, pgrndd                    &
                         , ptsoil                            &
                         , pgrndcapc, pgrndhflx              &
                         , lmask, ldglac)
!
!   AUTHOR:  FREDERIC HOURDIN     30/01/92
!
!            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
!            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
!
!            J.-P. SCHULZ   MPI - OCTOBER 1997 :
!               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
!               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
!               ECHAM4 GCM.
!            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
!            U.Schlese DKRZ - February 2000  new soil temperatures
!            L Kornblueh, MPI, January 2003, removed MERGE
!            
!            Adapted to JSBACH by Thomas Raddatz, Mai 2004
!
!   OBJECTIVE:  COMPUTATION OF:
!               THE GROUND TEMPERATURE EVOLUTION
!               THE GROUND SPECIFIC HEAT "CAPCAL"
!               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
!
!
!   METHOD:  IMPLICIT TIME INTEGRATION
!
!   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
!           T(K+1) = C(K) + D(K)*T(K)  (1)
!   THE COEFFICIENTS C (=pgrndc) AND D (=pgrndd) ARE COMPUTED AT THE
!   T-DT TIME-STEP.
!   ROUTINE STRUCTURE:
!   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
!   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
!     PROFILE FOR THE T+DT TIME-STEP
!   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
!     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
!            FDIFF = A + B TS(T+DT)
!     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
!            WITH F0 = A + B (TS(T))
!                 CAPCAL = B*DT
!
!     ------------------------------------------------------------------
!
!   DECLARATIONS:
!
USE mo_jsbach_constants   , ONLY: RhoH2O, snow_density
USE mo_time_control       , ONLY: lstart, delta_time
USE mo_kind               , ONLY: dp

IMPLICIT NONE
!
!-----------------------------------------------------------------------
!  ARGUMENTS
!
  INTEGER, Intent(in)  ::  klon                 !! length of the vector
  INTEGER, Intent(in)  ::  nsoil                !! number of soil layers 
  REAL(dp), Intent(in)     ::  cdel(nsoil)          !! soil layer thickness [m]
  REAL(dp), Intent(in)     ::  cmid(nsoil)          !! depth of mids of soil layers [m]
  REAL(dp), Intent(in)     ::  pts(klon)            !! surface temperature at top of soil [K]
  REAL(dp), Intent(in)     ::  psn(klon)            !! equivalent snow depth [m water]
  REAL(dp), Intent(in)     ::  psodif(klon)         !! soil temperature diffusivity [m^2/s]
  REAL(dp), Intent(in)     ::  prgcgn(klon)         !! soil heat capacity [J/m^3K]
  REAL(dp), Intent(inout)  ::  pgrndc(klon,nsoil)   !!
  REAL(dp), Intent(inout)  ::  pgrndd(klon,nsoil)   !!
  REAL(dp), Intent(inout)  ::  ptsoil(klon,nsoil)   !! soil temperature [K]
  REAL(dp), Intent(out)    ::  pgrndcapc(klon)      !!
  REAL(dp), Intent(out)    ::  pgrndhflx(klon)      !! ground heat flux
  LOGICAL, Intent(in)  ::  lmask(klon)
  LOGICAL, Intent(in)  ::  ldglac(klon)         !! glacier mask
!
!     ------------------------------------------------------------------
!
!  local Variables
!
  INTEGER :: jk
  REAL(dp) :: zso_cond(klon), zso_capa(klon)
  REAL(dp) :: z1(klon)
  REAL(dp) :: zd1(nsoil)
  REAL(dp) :: zdz1(klon,nsoil),   zdz2(klon,nsoil)
  REAL(dp) :: zkappa(klon,nsoil), zcapa(klon,nsoil)
  REAL(dp) :: zsnow_h(klon), zx1(klon), zx2(klon)
  REAL(dp) :: zrici, zdifiz, zsn_cond, zsn_capa
!
!     ------------------------------------------------------------------
!
!*    1.  SPECIFYING THE DEPTHS OF THE TEMPERATURE LEVELS.
!
!*    1.1 SOME CONSTANTS.
!
  zrici = 2.09e+06_dp                                 !! volumetric heat capacity of ice [j/m**3/k]
  zdifiz = 12.e-07_dp                                 !! temperature diffusivity of ice  [m**2/s]
  zsn_cond = 0.31_dp                                  !! snow thermal conductivity [j/s/m/k]
  zsn_capa = 634500.0_dp                              !! snow  heat capacity   [j/m**3/k]
!
!*    1.2 COMPUTING SOME USEFUL CONSTANTS.
!
  DO jk = 1,nsoil-1
     zd1(jk) = 1._dp / (cmid(jk+1) - cmid(jk))
  END DO
!
!*    1.3 COMPUTE OF THE SOIL THERMAL CONDUCTIVITY [J/S/M/K] FROM
!*        THE SOIL TEMPERATURE DIFFUSIVITY [M**2/S].
!
  WHERE (lmask(:) .AND. ldglac(:))
    zso_capa(:) = zrici
    zso_cond(:) = zso_capa(:) * zdifiz
  ELSEWHERE (lmask(:))
    zso_capa(:) = prgcgn(:)
    zso_cond(:) = zso_capa(:) * psodif(:)
  END WHERE
!
!*    1.4 PRE-SET THERMAL CONDUCTIVITY AT ALL LEVELS.
!
  DO jk = 1,nsoil
     WHERE (lmask)
        zkappa(:,jk) = zso_cond(:)
        zcapa(:,jk)  = zso_capa(:)
     END WHERE
  END DO
!
!   --------------------------------------------------------------
!   COMPUTATION OF THE GROUND TEMPERATURES USING THE CGRD AND DGRD
!   COEFFICIENTS COMPUTED AT THE PREVIOUS TIME-STEP
!   --------------------------------------------------------------
!
  IF(.NOT.lstart) THEN
!
!   Upper layer
!
    ptsoil(:,1) = pts(:)
!
!   Deeper layers
!
    DO jk = 1,nsoil-1
      WHERE(lmask) ptsoil(:,jk+1) = pgrndc(:,jk) + pgrndd(:,jk) * ptsoil(:,jk)
    END DO
  END IF
!
!   ---------------------------------------------------------------
!   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
!   ---------------------------------------------------------------
!
  WHERE (lmask) 
     zsnow_h(:) = psn(:) * RhoH2O / snow_density
!
!*       Special treatment for first layer
!
     WHERE ( zsnow_h(:) > cmid(2) )
        zcapa(:,1) = zsn_capa
        zkappa(:,1) = zsn_cond
     ELSEWHERE( zsnow_h(:) > 0.0_dp .AND. zsnow_h(:) <= cmid(2) )
        zx1 = zsnow_h(:) / cmid(2)
        zx2 = ( cmid(2) - zsnow_h(:)) / cmid(2)
        zcapa(:,1) = zx1 * zsn_capa + zx2 * zso_capa(:)
        zkappa(:,1) = 1._dp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
     ELSEWHERE
        zcapa(:,1) = zso_capa(:)
        zkappa(:,1) = zso_cond(:)
     ENDWHERE
  END WHERE
!
  DO jk = 2, nsoil - 2
    WHERE (lmask)
       WHERE ( zsnow_h(:) > cmid(jk+1) )
          zcapa(:,jk) = zsn_capa
          zkappa(:,jk) = zsn_cond
       ELSEWHERE ( zsnow_h(:) > cmid(jk) .AND. zsnow_h(:) <= cmid(jk+1) )
          zx1 = (zsnow_h(:) - cmid(jk)) * zd1(jk)
          zx2 = ( cmid(jk+1) - zsnow_h(:)) * zd1(jk)
          zcapa(:,jk) = zx1 * zsn_capa + zx2 * zso_capa(:)
          zkappa(:,jk) = 1._dp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
       ELSEWHERE
          zcapa(:,jk) = zso_capa(:)
          zkappa(:,jk) = zso_cond(:)
       ENDWHERE
    END WHERE
  END DO
!
  DO jk=1,nsoil
    WHERE (lmask) zdz2(:,jk) = zcapa(:,jk) * cdel(jk) / delta_time
  END DO
!
  DO jk=1,nsoil-1
    WHERE (lmask) zdz1(:,jk) = zd1(jk) * zkappa(:,jk)
  END DO
!
  WHERE (lmask)
     z1(:) = zdz2(:,nsoil) + zdz1(:,nsoil-1)
     pgrndc(:,nsoil-1) = zdz2(:,nsoil) * ptsoil(:,nsoil) / z1(:)
     pgrndd(:,nsoil-1) = zdz1(:,nsoil-1) / z1(:)
  END WHERE
!
  DO jk=nsoil-1,2,-1
     WHERE (lmask)
        z1(:) = 1._dp / (zdz2(:,jk) + zdz1(:,jk-1) + zdz1(:,jk) * (1._dp - pgrndd(:,jk)))
        pgrndc(:,jk-1) = (ptsoil(:,jk) * zdz2(:,jk) + zdz1(:,jk) * pgrndc(:,jk)) * z1(:)
        pgrndd(:,jk-1) = zdz1(:,jk-1) * z1(:)
     END WHERE
  END DO
!
!   ---------------------------------------------------------
!   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
!   CALORIFIC CAPACITY OF THE GROUND:
!   ---------------------------------------------------------
!
  WHERE (lmask)
     pgrndhflx(:) = zdz1(:,1) * (pgrndc(:,1) + (pgrndd(:,1) - 1._dp) * ptsoil(:,1))
     pgrndcapc(:) = (zdz2(:,1) * delta_time + delta_time * (1._dp - pgrndd(:,1)) * zdz1(:,1))
  ELSEWHERE
     pgrndhflx = 0._dp
     pgrndcapc = 0._dp
  END WHERE

END SUBROUTINE update_soiltemp

SUBROUTINE update_soiltemp_5lyr (klon, nsoil, zsoil, smid    &
                         , pts, psn                          &
                         , pgrndc, pgrndd                    &
                         , ptsoil                            &
                         , pgrndcapc, pgrndhflx              &
                         , lmask, ldglac                     &
                         , snow_c, snow_d, snow_t            &
                         , org_c, org_d, org_t               &
                         , ice, c, k, heatcap, heatcond      &
                         , dzs, FMPOT, BCLAPP, VPOR, VFC, WSI&
                         , thaw_depth, k_snow, c_snow        &
                         , nsnwlyr, snw_density, psn_old )
!
!   AUTHOR:  FREDERIC HOURDIN     30/01/92
!
!            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
!            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
!
!            J.-P. SCHULZ   MPI - OCTOBER 1997 :
!               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
!               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
!               ECHAM4 GCM.
!            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
!            U.Schlese DKRZ - February 2000  new soil temperatures
!            L Kornblueh, MPI, January 2003, removed MERGE
!
!            Adapted to JSBACH by Thomas Raddatz, Mai 2004
!
!   OBJECTIVE:  COMPUTATION OF:
!               THE GROUND TEMPERATURE EVOLUTION
!               THE GROUND SPECIFIC HEAT "CAPCAL"
!               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
!
!
!   METHOD:  IMPLICIT TIME INTEGRATION
!
!   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
!           T(K+1) = C(K) + D(K)*T(K)  (1)
!   THE COEFFICIENTS C (=pgrndc) AND D (=pgrndd) ARE COMPUTED AT THE
!   T-DT TIME-STEP.
!   ROUTINE STRUCTURE:
!   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
!   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
!     PROFILE FOR THE T+DT TIME-STEP
!   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
!     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
!            FDIFF = A + B TS(T+DT)
!     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
!            WITH F0 = A + B (TS(T))
!                 CAPCAL = B*DT
!
!
!     ------------------------------------------------------------------
!
!     MODIFIED BY ALTUG EKICI
!
!     INCLUDED:
!     -PHASE CHANGE 
!     -HEAT TRANSFER PARAMETERS AS A FUNCTION OF SOIL WATER AND ICE
!     -SUPERCOOLED WATER EQUATION
!     -THAW DEPTH CALCULATION
!     -MULTI LAYERED SNOW SCHEME and KEEPING 5 SOIL LAYERS AT ALL TIMES 
!     -ORGANIC LAYER ABOVE THE SOIL
!     
!     ------------------------------------------------------------------
!
!   DECLARATIONS:
!
USE mo_time_control       , ONLY: lstart, delta_time
USE mo_kind               , ONLY: dp
USE mo_physical_constants , ONLY: rhoh2o, rhoice, alf, cpliq, cpice, cpsno, alice, alsno, tmelt, grav
USE mo_jsbach_constants   , ONLY: snow_density
USE mo_land_surface       , ONLY: fract_small

IMPLICIT NONE
!
!-----------------------------------------------------------------------
!  ARGUMENTS
!
INTEGER,  INTENT(in)    ::  klon                 !! length of the vector
INTEGER,  INTENT(in)    ::  nsoil                !! number of soil layers
REAL(dp), INTENT(in)    ::  zsoil(nsoil)         !! Thickness of individual soil layers [m]
REAL(dp), INTENT(in)    ::  smid(nsoil)          !! Midpoint depths of soil layers below surface [m]
REAL(dp), INTENT(in)    ::  pts(klon)            !! surface temperature at top of soil [K]
REAL(dp), INTENT(in)    ::  psn(klon)            !! equivalent snow depth [m water]
REAL(dp), INTENT(inout) ::  pgrndc(klon,nsoil)   !!
REAL(dp), INTENT(inout) ::  pgrndd(klon,nsoil)   !!
REAL(dp), INTENT(inout) ::  ptsoil(klon,nsoil)   !! soil temperature [K]
REAL(dp), INTENT(out)   ::  pgrndcapc(klon)      !! ground heat capacity
REAL(dp), INTENT(out)   ::  pgrndhflx(klon)      !! ground heat flux
LOGICAL,  INTENT(in)    ::  lmask(klon)          !! land mask
LOGICAL,  INTENT(in)    ::  ldglac(klon)         !! glacier mask
REAL(dp), INTENT(inout) ::  snow_c(klon,nsoil)
REAL(dp), INTENT(inout) ::  snow_d(klon,nsoil)
REAL(dp), INTENT(inout) ::  snow_t(klon,nsoil)   !! snow temperature [K]
REAL(dp), INTENT(inout) ::  org_c(klon)
REAL(dp), INTENT(inout) ::  org_d(klon)
REAL(dp), INTENT(inout) ::  org_t(klon)          !! org layer temperature [K]
REAL(dp), INTENT(in)    ::  c(klon)              !! soil mineral volumetric heat capacity [J/m3/K]
REAL(dp), INTENT(in)    ::  k(klon)              !! soil mineral heat conductivity [W/m/K]
REAL(dp), INTENT(out)   ::  heatcap(klon,nsoil)  !! modified soil heat capacity [J/m3/K]
REAL(dp), INTENT(out)   ::  heatcond(klon,nsoil) !! modified soil heat conductivity [W/m/K]
REAL(dp), INTENT(inout) ::  ice(klon,nsoil)      !! soil ice content [m]
REAL(dp), INTENT(inout) ::  WSI(klon,nsoil)      !! soil water content [m]
REAL(dp), INTENT(in)    ::  dzs(klon)            !! soil depth until bedrock [m]
REAL(dp), INTENT(in)    ::  FMPOT(klon)          !! saturated moisture potential [m]
REAL(dp), INTENT(in)    ::  BCLAPP(klon)         !! Exponent b of Clapp and Hornberger
REAL(dp), INTENT(in)    ::  VPOR(klon)           !! Volumetric soil porosity [m/m]
REAL(dp), INTENT(in)    ::  VFC(klon)            !! Volumetric soil field capacity [m/m]
REAL(dp), INTENT(out)   ::  thaw_depth(klon)     !! Thaw depth [m]
REAL(dp), INTENT(out)   ::  k_snow(klon)         !! Snow heat conductivity [W/m/K]
REAL(dp), INTENT(out)   ::  c_snow(klon)         !! Snow heat capacity [J/m3/K]
REAL(dp), INTENT(inout) ::  nsnwlyr(klon)        !! Number of snow layers in last timestep
REAL(dp), INTENT(inout) ::  snw_density(klon)    !! Flexible snow density [kg/m3]
REAL(dp), INTENT(inout) ::  psn_old(klon)        !! snow water equivalent of previous timestep [m]
!
!     ------------------------------------------------------------------
!
!  local Variables
!
INTEGER  :: jl,jk
INTEGER  :: nsnow,n_org,nlayers,nlayers_tot,extra,delta_snw
LOGICAL  :: lsnow,lheatcap,lheatcond,lfreeze,lsupercool,lorglayer,ldynsnow,ldynmoss
REAL(dp) :: org_cond,org_capa,z_org(1),snwmindens,snwmaxdens,pack_density,newsnow
REAL(dp) :: melt,frozen,KE,ksat,sat,kdry,rho_bulk,rho_soil,alwat,k_bedrock,c_bedrock
REAL(dp) :: v_wsi, v_ice
REAL(dp) :: z1(klon),zsnow_h(klon)
REAL(dp) :: zso_cond(klon,nsoil), zso_capa(klon,nsoil),  wsi_max(klon,nsoil)
REAL(dp) :: acc_depth(0:nsoil), sat_depth(nsoil)
REAL(dp), POINTER :: zlayers(:)
REAL(dp), POINTER :: zd1(:),cdel(:),cmid(:)
REAL(dp), POINTER :: zdz1(:),zdz2(:),zkappa(:),zcapa(:),c_coef(:),d_coef(:),lyr_t(:)

! rho_bulk:   soil bulk density [kg/m3]
! kdry    :   soil dry heat conductivity
! KE      :   Kersten number
! ksat    :   soil saturated heat conductivity
! sat     :   saturation degree 
! melt    :   amount of melted ice
! frozen  :   amount of frozen water
! v_wsi   :   volumetric fraction of soil water with respect to air-free part of soil
! v_ice   :   volumetric fraction of soil ice   with respect to air-free part of soil

! SWITCHES

! important note: when you remove organic layer, update this line in mo_soil:
! "soil%soil_temperature(kidx0:kidx1,1,itile) = surface_temperature_avg(:)"
  lsnow      = .TRUE.  !switch for the snow scheme
  ldynsnow   = .FALSE. !switch for dynamic snow parameters
  lfreeze    = .TRUE.  !switch for phase change
  lheatcap   = .TRUE.  !switch for heat capacity
  lheatcond  = .TRUE.  !switch for heat conductivity
  lsupercool = .TRUE.  !switch for supercooled water
  lorglayer  = .TRUE.  !switch for organic layer
  ldynmoss   = .FALSE. !switch for dynamic organic layer parameters

  ! CONSTANTS
  snwmindens = 50._dp     !! minimum snow density [kg/m3]
  snwmaxdens = 300._dp    !! maximum snow density [kg/m3]

  alwat     = 0.57_dp     !! thermal conductivity of liquid water [W/K/m]    (Hillel 1982)
  k_bedrock = 2._dp       !! thermal conductivity of bedrock [W/m/K]         (Tarnocai 2004)
  c_bedrock = 2000000._dp !! volumetric heat capacity of bedrock [J/m3/K]    (Bonan 2002)
  rho_soil  = 2700._dp    !! typical soil partical density [kg/m3]           (Hillel 1980, Johansen 1977)
  org_cond  = 0.25_dp     !! organic layer heat conductivity [W/m/K]         (Beringer etal 2001)
  org_capa  = 2500000._dp !! organic layer volumetric heat capacity [J/m3/K] (Beringer etal 2001)
  n_org     = 1           !! number of organic layers (actually 1 is the only supported value)
  z_org     = (/0.1_dp/)  !! organic layer depths [m]

  ! Initializations
  k_snow    (:) = 0._dp
  c_snow    (:) = 0._dp
  thaw_depth(:) = 0._dp

  ! Allocate layer data and set fixed soil depths
  IF (.NOT. lorglayer) n_org = 0
  nlayers_tot = 5+n_org+nsoil
  ALLOCATE(zlayers(nlayers_tot),cmid (nlayers_tot),zd1   (nlayers_tot),zdz1 (nlayers_tot),zdz2(nlayers_tot),    &
           zkappa (nlayers_tot),zcapa(nlayers_tot),c_coef(nlayers_tot),d_coef(nlayers_tot),lyr_t(nlayers_tot+1))
  zlayers(1:5) = 0.05_dp
  zlayers(6+n_org:nlayers_tot) = zsoil(:)
  IF (lorglayer) zlayers(6) = z_org(1) ! Should be moved inside the grid-box loop if thickness of organic layer becomes dynamic

  ! GRIDBOXES
  DO jl = 1,klon  !for each gridbox
    IF (lmask(jl)) THEN  !only for land points

      IF (lorglayer) THEN
        cdel => zlayers(6:nlayers_tot)
        nlayers = nsoil+n_org
      ELSE
        cdel => zlayers(6+n_org:nlayers_tot)
        nlayers = nsoil
      ENDIF

      IF (lsnow) THEN !snow switch is on

        IF (ldynsnow) THEN !dynamic snow density
          !SNOW density and depth calculations
          IF (nsnwlyr(jl) == 0) then    !there was no snow before so use minimum density
            snw_density(jl) = snwmindens
          ELSE  ! calculation from Verseghy,1991
            pack_density   = (snw_density(jl)-snwmaxdens)*EXP(-0.01_dp*delta_time/3600._dp)+snwmaxdens
            IF (psn(jl) > psn_old(jl)) THEN  !new snowfall
              newsnow = psn(jl) - psn_old(jl)
              snw_density(jl) = (pack_density*psn_old(jl) + snwmindens*newsnow) / psn(jl)
            ELSE
              snw_density(jl) = pack_density
            END IF
          END IF
        ELSE
          snw_density(jl) = snow_density
        END IF

        zsnow_h(jl) = psn(jl) * rhoh2o / snw_density(jl)    ! calculating snow depth

        ! Distribute snow on snow_layers
        ! StW: This algorithm copies the original scheme of Altug, but does in general not conserve show depth - may need revision
        ! The snow depths are mapped as follows (layers from top to bottom):
        !        zsnow_h <= 0.02 m => no snow layers
        ! 0.02 < zsnow_h <= 0.05 m => 1 5 cm snow layer
        ! 0.05 < zsnow_h <= 0.10 m => 2 5 cm snow layers
        ! 0.10 < zsnow_h <= 0.15 m => 3 5 cm snow layers
        ! 0.15 < zsnow_h <= 0.25 m => 4 5 cm snow layers
        ! 0.25 < zsnow_h           => 4 5 cm and 1 (zsnow_h - 0.2 m) snow layers
        IF (zsnow_h(jl) > 0.25_dp) THEN
          zlayers(5) = zsnow_h(jl) - .2_dp
          nlayers    = nlayers_tot
          cdel => zlayers(1:nlayers_tot)
        ELSE
          zlayers(5) = 0.05_dp
          nlayers = MIN(4,INT(CEILING(zsnow_h(jl)/0.05_dp)))
          IF (zsnow_h(jl) <= 0.02_dp) nlayers = 0
          cdel => zlayers(6-nlayers:nlayers_tot)
          nlayers = nlayers + nsoil + n_org
        ENDIF
      END IF

      lyr_t(nlayers+1) = 0._dp
      lyr_t(nlayers-nsoil+1:nlayers) = ptsoil(jl,:)

      !calculating mid layer depths
      cmid(1) = cdel(1)/2._dp
      DO jk = 2, nlayers
        cmid(jk) = (cdel(jk)+cdel(jk-1))/2._dp+cmid(jk-1)
        zd1(jk-1) = 1._dp / (cmid(jk) - cmid(jk-1))
      END DO

      !calculating water saturatable depths for each layer
      acc_depth(0) = 0._dp
      acc_depth(1) = zsoil(1)
      sat_depth(1) = zsoil(1)
      DO jk = 2, nsoil
        acc_depth(jk) = acc_depth(jk-1) + zsoil(jk)

        IF ((acc_depth(jk) > dzs(jl)) .AND. (dzs(jl) > acc_depth(jk-1))) THEN
          sat_depth(jk) = dzs(jl) - acc_depth(jk-1)
        ELSE
          sat_depth(jk) = zsoil(jk)
        END IF
      END DO

      extra = nlayers-nsoil  !number of layers above soil
      nsnow = extra-n_org    !number of snow layers
      delta_snw = nsnow - nsnwlyr(jl)  !difference in number of snow layers from last timestep

      IF (nsnow > 0) THEN
        !****************************
        !******  SNOWLAYERS  ********
        !****************************
        IF (ldynsnow) THEN
           IF (nsnwlyr(jl) == 0) then    !there was no snow before
              k_snow(jl) = 0.1_dp
              c_snow(jl) = 526500._dp
           ELSE
              !k_snow(jl) = (2.576_dp*(10._dp**(-6._dp))*(snw_density(jl)**2._dp))+0.074_dp  ! (Verseghy, 1991)
              k_snow(jl) = 2.9_dp*10**(-6._dp)*(snw_density(jl)**2._dp)  !(Goodrich,1981))
              c_snow(jl) = cpice * snw_density(jl)
           END IF
        ELSE
           k_snow(jl) = alsno
           c_snow(jl) = cpsno*snw_density(jl)
        END IF

        DO jk = 1, nsnow  !for all the snow layers
          zkappa(jk) = k_snow(jl)
          zcapa(jk)  = c_snow(jl)
          IF ((delta_snw > 0) .AND. (jk <= delta_snw)) THEN !for the new snow layers 
            IF (nsnwlyr(jl) == 0) THEN !there was no snow before, so using first soil layer  values
              c_coef(jk) = pgrndc(jl,1)
              d_coef(jk) = pgrndd(jl,1)
              lyr_t(jk)  = ptsoil(jl,1)
            ELSE  !use top snow layer values for the new layers
              c_coef(jk) = snow_c(jl,1)
              d_coef(jk) = snow_d(jl,1)
              lyr_t(jk)  = snow_t(jl,1)
            END IF
          ELSE !there has been snow melt, remove the layer at the top, use the bottom layers
            c_coef(jk) = snow_c(jl,jk-delta_snw)
            d_coef(jk) = snow_d(jl,jk-delta_snw)
            lyr_t (jk) = snow_t(jl,jk-delta_snw)
          END IF

          IF (.NOT. lstart) THEN
            IF (jk == 1) THEN  ! Surface layer
              lyr_t(jk) = pts(jl)
            ELSE               ! Deeper layers
              lyr_t(jk) = c_coef(jk-1) + d_coef(jk-1) * lyr_t(jk-1)
            END IF
          END IF
        END DO
      END IF

      IF (lorglayer) THEN  !for the organic layers
        !******************************
        !******  ORGANICLAYERS ********
        !******************************
        DO jk = 1, n_org  !for all the organic layers

          IF (ldynmoss) THEN !dynamic organic layer parameters
            zcapa(jk+nsnow) = (1.0_dp-0.8_dp)*org_capa+(WSI(jl,1)/sat_depth(1))*cpliq*rhoh2o+&
                              (ice(jl,1)/sat_depth(1))*cpice*rhoice

            kdry = 0.05_dp !value from O'Donnel et al 2009
            sat = MAX(0._dp,MIN(1._dp,((WSI(jl,1)+ice(jl,1))/sat_depth(1))/0.8_dp))
            IF (sat > 0.1_dp) THEN
              KE = LOG10(sat)+1._dp
            ELSE
              KE = 0.0_dp !log10(0.1_dp)+1._dp
            END IF
            v_wsi = wsi(jl,1)/ ((1._dp-0.8_dp+VFC(jl))*sat_depth(1))
            v_ice = ice(jl,1)/ ((1._dp-0.8_dp+VFC(jl))*sat_depth(1))
            ksat            = (org_cond**(1._dp-(v_wsi+v_ice)))*(alwat**(v_wsi))*(alice**(v_ice))

            zkappa(jk+nsnow) = ksat*KE + (1._dp-KE)*kdry
          ELSE
            zkappa(jk+nsnow) = org_cond
            zcapa(jk+nsnow)  = org_capa
          END IF

          c_coef(jk+nsnow) = org_c(jl)
          d_coef(jk+nsnow) = org_d(jl)
          lyr_t(jk+nsnow)  = org_t(jl)

          IF (.NOT. lstart) THEN
            IF ((nsnow == 0) .AND. (jk == 1)) THEN  ! Surface layer
              lyr_t(jk) = pts(jl)
            ELSE               ! Deeper layers
              lyr_t(jk+nsnow) = c_coef(jk+nsnow-1) + d_coef(jk+nsnow-1) * lyr_t(jk+nsnow-1)
            END IF
          END IF
        END DO
      END IF

      !****************************
      !******  SOILLAYERS  ********
      !****************************
      DO jk = 1, nsoil
        IF (ldglac(jl)) THEN  !just for glaciers
          zso_capa(jl,jk)  = cpice*rhoice
          zso_cond(jl,jk)  = alice
        ELSE !just for glacier-free land points

          !soil volumetric heat capacity calculation from deVries,2963
          IF (lheatcap) THEN !if the heat capacity switch is on
            IF (jk > 1 .AND. acc_depth(jk-1) > dzs(jl)) THEN !bedrock
              zso_capa(jl,jk) = c_bedrock
            ELSE
              zso_capa(jl,jk) = (1._dp-VPOR(jl))*c(jl)+(WSI(jl,jk)/sat_depth(jk))*cpliq*rhoh2o+&
                                (ice(jl,jk)/sat_depth(jk))*cpice*rhoice
              IF (sat_depth(jk) < zsoil(jk)) THEN !special condition for semi-saturatable layer (soil and bedrock mixture)
                 zso_capa(jl,jk) = (zso_capa(jl,jk)*sat_depth(jk)+c_bedrock*(zsoil(jk)-sat_depth(jk)))/zsoil(jk)
              END IF
            END IF
          ELSE
            zso_capa(jl,jk) = c(jl)
          END IF

          !soil heat conductivity calculation from Johansen,1977
          rho_bulk        = rho_soil*(1._dp-VPOR(jl))
          kdry            = (0.135_dp*rho_bulk+64.7_dp)/(rho_soil-0.947_dp*rho_bulk)
          IF (lheatcond) THEN !if the heat conductivity switch is on
            IF (jk > 1 .AND. acc_depth(jk-1) < dzs(jl)) THEN !bedrock
               zso_cond(jl,jk) = k_bedrock
            ELSE
              sat = ((WSI(jl,jk)+ice(jl,jk))/sat_depth(jk))/VPOR(jl)
              IF (lyr_t(jk+extra) > 0.0_dp) THEN !unfrozen case
                IF(sat > 0.1_dp) THEN
                  KE = LOG10(sat)+1
                ELSE
                  KE = 0.0_dp !log10(0.1)+1
                END IF
              ELSE !frozen case
                KE = sat
              END IF
              v_wsi = wsi(jl,jk)/ ((1._dp-VPOR(jl)+VFC(jl))*sat_depth(jk))
              v_ice = ice(jl,jk)/ ((1._dp-VPOR(jl)+VFC(jl))*sat_depth(jk))
              ksat  = (k(jl)**(1._dp-(v_wsi+v_ice)))*(alwat**(v_wsi))*(alice**(v_ice))

              !actual soil layer heat conductivity using dry and saturated states as weighted by Kersten number
              zso_cond(jl,jk) = ksat*KE + kdry*(1._dp-KE)

              IF (sat_depth(jk) < zsoil(jk)) THEN !special condition for semi-saturatable layer (soil and bedrock mixture)
                zso_cond(jl,jk) = (zso_cond(jl,jk)*sat_depth(jk)+k_bedrock*(zsoil(jk)-sat_depth(jk)))/zsoil(jk)
              END IF
            END IF
          ELSE
            zso_cond(jl,jk) = kdry
          END IF
        END IF

        zkappa(jk+extra) = zso_cond(jl,jk)
        zcapa (jk+extra) = zso_capa(jl,jk)
        c_coef(jk+extra) = pgrndc  (jl,jk)
        d_coef(jk+extra) = pgrndd  (jl,jk)
        lyr_t (jk+extra) = ptsoil  (jl,jk)

        IF (.NOT. lstart) THEN

          IF ((extra == 0) .AND. (jk == 1)) THEN  ! Surface layer
            lyr_t(jk) = pts(jl)
          ELSE               ! Deeper layers
            lyr_t(jk+extra) = c_coef(jk+extra-1) + d_coef(jk+extra-1) * lyr_t(jk+extra-1)
          END IF

          !phase change calculations in soil layers
          IF (lfreeze) THEN !if the freeze/thaw switch on
            IF (lyr_t(jk+extra) < tmelt) THEN
              !***** FREEZING *****
              IF (lsupercool) THEN
              !SUPERCOOLED WATER equation (Niu&Yang,2006)
                wsi_max(jl,jk)= sat_depth(jk)*VPOR(jl)*(((alf*(tmelt-lyr_t(jk+extra)))&
                                /MAX(fract_small,grav*lyr_t(jk+extra)*FMPOT(jl)))**(-1._dp/MAX(1._dp,BCLAPP(jl))))
              ELSE
                wsi_max(jl,jk) = 0._dp
              END IF
              IF (WSI(jl,jk) > wsi_max(jl,jk)) THEN
                frozen = MIN((WSI(jl,jk)-wsi_max(jl,jk)),((sat_depth(jk)*zcapa(jk+extra)*&
                         (tmelt-lyr_t(jk+extra)))/(alf*rhoh2o)))
                WSI(jl,jk)   = WSI(jl,jk) - frozen
                ice(jl,jk)   = ice(jl,jk) + frozen
                lyr_t(jk+extra) = lyr_t(jk+extra) + ((frozen*alf*rhoh2o)/(zcapa(jk+extra)*sat_depth(jk)))
              END IF
            ELSE IF ((lyr_t(jk+extra) > tmelt) .AND. (ice(jl,jk) > 0._dp)) THEN
              !****** MELTING ******
              melt   = MIN(ice(jl,jk),((sat_depth(jk)*zcapa(jk+extra)*(lyr_t(jk+extra)-tmelt))/(alf*rhoice)))
              WSI(  jl,jk)    = WSI(jl,jk) + melt
              ice  (jl,jk)    = ice(jl,jk) - melt
              lyr_t(jk+extra) = lyr_t(jk+extra) - ((melt*alf*rhoice)/(zcapa(jk+extra)*sat_depth(jk)))
            END IF
          END IF !end of freeze/thaw switch if
        END IF !end of lstart if
      END DO !end of loop over layers

      !***** THAWING DEPTH *******
      IF (.NOT. lstart) THEN
        IF (lyr_t(extra+1) <= tmelt) THEN  ! fully-frozen soil grid, assign zero
          thaw_depth(jl) = 0.0_dp
        ELSE
          thaw_depth(jl) = smid(nsoil) !first assign the maximum depth in case its fully thawed
          DO jk = 1,nsoil-1
            IF ((lyr_t(jk+extra) > tmelt) .AND. (lyr_t(jk+extra+1) <= tmelt)) THEN
              thaw_depth(jl) = smid(jk) + (smid(jk+1)-smid(jk)) * (lyr_t(jk+extra)-tmelt)             &
                                                                / (lyr_t(jk+extra)-lyr_t(jk+extra+1))
            END IF
          END DO  ! layers
        END IF  ! lyr_t < tmelt
      END IF  ! lstart
      
!
!   ---------------------------------------------------------------
!   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
!   ---------------------------------------------------------------
!
      DO jk=1,nlayers
        zdz2(jk) = zcapa(jk) * cdel(jk) / delta_time
      END DO

      DO jk=1,nlayers-1
        zdz1(jk) = zd1(jk) * zkappa(jk)
      END DO

      z1(jl) = zdz2(nlayers) + zdz1(nlayers-1)
      c_coef(nlayers-1) = zdz2(nlayers) * lyr_t(nlayers) / z1(jl)
      d_coef(nlayers-1) = zdz1(nlayers-1) / z1(jl)

      DO jk=nlayers-1,2,-1
        z1(jl) = 1._dp / (zdz2(jk) + zdz1(jk-1) + zdz1(jk) * (1._dp - d_coef(jk)))
        c_coef(jk-1) = (lyr_t(jk) * zdz2(jk) + zdz1(jk) * c_coef(jk)) * z1(jl)
        d_coef(jk-1) = zdz1(jk-1) * z1(jl)
      END DO
!
!   ---------------------------------------------------------
!   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
!   CALORIFIC CAPACITY OF THE GROUND:
!   ---------------------------------------------------------
!
      pgrndhflx(jl) = zdz1(1) * (c_coef(1) + (d_coef(1) - 1._dp) * lyr_t(1))
      pgrndcapc(jl) = (zdz2(1) * delta_time + delta_time * (1._dp - d_coef(1)) * zdz1(1))

      !  assigning values for output
      IF (lsnow) THEN
        psn_old(jl) = psn(jl) !keep the amount of snow for the next timestep
        nsnwlyr(jl) = nsnow  !keep the number of snow layers for the next timestep
        IF(nsnow > 0) THEN
          DO jk = 1, nsnow     !for the snow layers
            snow_c(jl,jk) = c_coef(jk)
            snow_d(jl,jk) = d_coef(jk)
            snow_t(jl,jk) = lyr_t(jk)
          END DO
        END IF
      END IF
      IF (lorglayer) THEN   !for the organic layers
        DO jk = 1, n_org
          org_c(jl) = c_coef(jk+nsnow)
          org_d(jl) = d_coef(jk+nsnow)
          org_t(jl) = lyr_t(jk+nsnow)
        END DO
      END IF
      DO jk = 1, nsoil        !for the soil layers
        pgrndc(jl,jk)  = c_coef(jk+extra)
        pgrndd(jl,jk)  = d_coef(jk+extra)
        ptsoil(jl,jk)  = lyr_t(jk+extra)
        heatcap(jl,jk) = zso_capa(jl,jk)
        heatcond(jl,jk)= zso_cond(jl,jk)
      END DO
    ELSE !not a land point
      pgrndhflx(jl)   = 0._dp
      pgrndcapc(jl)   = 0._dp
      heatcap  (jl,:) = 0._dp
      heatcond (jl,:) = 0._dp
    END IF

  END DO !end of loop over grid boxes
  DEALLOCATE (zlayers,cmid,zd1,zdz1,zdz2,zkappa,zcapa,c_coef,d_coef,lyr_t)

END SUBROUTINE update_soiltemp_5lyr

