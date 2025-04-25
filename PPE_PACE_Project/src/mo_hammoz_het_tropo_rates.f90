!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_hammoz_het_tropo_rates.f90
!!
!! \brief
!!  This module includes tropospheric heterogeneous chemistry subroutines
!!
!! \author M. Schultz   (FZ Juelich)
!! \author S. Stadtler  (Uni Bonn)
!! \author S. Schroeder (FZ Juelich)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2015-02-20)
!!
!! \limitations
!! only runs with also HAM module running
!!
!! \details
!! None
!!
!! \bibliographic_references
!! Evans, 2005
!! Liao, 2005 (dust)
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_hammoz_het_tropo_rates

  USE mo_kind,               ONLY: dp
  USE mo_parameters
  USE mo_time_control,       ONLY: time_step_len,delta_time
  USE mo_math_constants,     ONLY: pi 
  USE mo_physical_constants, ONLY: grav, avo, argas, tmelt, rhoh2o, &
                                   zmolgair => amd, zmolgs => ams, zmolgso4 => amso4
  USE mo_tracdef,            ONLY: ntrac, trlist, AEROSOLMASS
  USE mo_tracer,             ONLY: get_tracer
  USE mo_ham,                ONLY: nmaxclass, nclass, nham_subm, HAM_M7, HAM_SALSA
  USE mo_ham_m7ctl,          ONLY: cmr2ras,   &     
                                   inucs, iaits, iaccs,  &
                                   icoas, iacci, icoai
  USE mo_boundary_condition, ONLY: bc_nml


  IMPLICIT NONE

  PRIVATE

  ! public variables  (see declaration below)

  ! subprograms
  PUBLIC :: ratecon_sftropo_inti
  PUBLIC :: ratecon_sftropo
  PUBLIC :: bc_aerosol
  PUBLIC :: bc_sad_cloud
  PUBLIC :: NAERO_SPEC
  PUBLIC :: aero_spec

  !-- name list fine control of boundary conditions
  TYPE(bc_nml), SAVE          :: bc_sad_cloud, &! soll ich hier ueberhaupt ein neues definieren? Eigentlich unnoetig 
                                 bc_aerosol     ! variable names fixed
                                                ! - only allowed to change ef_type, ef_geometry 
                                                ! normally, in HAMMOZ, all fields will be taken from HAM (ef_type=EF_MODULE)

  ! internal tracer indices (don't confuse with idt_xxx!)
  INTEGER, PARAMETER          :: idx_no3 = 1,     &
                                 idx_n2o5 = 2,    &
                                 idx_ho2 = 3,     &
                                 idx_hno3 = 4,    &
                                 idx_o3 = 5,      &
                                 idx_no2 = 6

  INTEGER, PARAMETER          :: MAXTRAC = 6          ! maximum number of gas-phase tracers
                                                      ! (must be consistent with idx_xxx below)
  INTEGER, PARAMETER          :: NAERO_SPEC = 5       ! maximum number of aerosol components
  INTEGER, PARAMETER          :: idx_aero_so4 = 1, &
                                 idx_aero_bc  = 2, &
                                 idx_aero_oc  = 3, &
                                 idx_aero_ss  = 4, &
                                 idx_aero_du  = 5
  CHARACTER(LEN=3), PARAMETER :: aero_spec(NAERO_SPEC) = (/'SO4', 'BC ','OC ','SS ','DU '/)

  ! boundary condition indices
  INTEGER, PUBLIC          ::   ibc_rwet(nmaxclass) = -1,    &    ! wet radii per mode
                                ibc_rdry(nmaxclass) = -1,    &    ! dry radii per mode
                                ibc_numden_modes(nmaxclass) = -1, &
                                ibc_aero_comp_modes(nmaxclass,NAERO_SPEC) = -1,&
                                ibc_radl = -1,&                  ! mean volume radius
                                ibc_cdnc = -1                    ! cloud droplet number concentration
 

  ! tracer indices (these are the normal idt_xxx values!)
  INTEGER, SAVE              :: idt(MAXTRAC)
  ! molecular weights
  REAL(dp)                   :: mweight(MAXTRAC)
  REAL(dp)                   :: pmpma(MAXTRAC)
  REAL(dp)                   :: const, dg_const1, dg_const2


  CONTAINS

  !!----------------------------------------------------------------------------
  SUBROUTINE het_tropo_m7_boundary_condition_inti()
  USE mo_tracer,                   ONLY: get_tracer
  USE mo_ham,                      ONLY: sizeclass
  USE mo_boundary_condition,       ONLY: bc_define

  ! initialisation routine

  !--- local variables
  INTEGER                    :: imod,jt,ierr,idx
  LOGICAL                    :: flag_aero_comp_modes(NAERO_SPEC,nclass)

  bc_aerosol%ef_actual_unit = 'm'
  DO imod = 1, nclass
      bc_aerosol%ef_varname = "rwet_"//TRIM(sizeclass(imod)%shortname)
      ibc_rwet(imod) = bc_define('aerosol wet radius: '//TRIM(sizeclass(imod)%classname), bc_aerosol, 3, .TRUE.)
  ENDDO

  ibc_aero_comp_modes = -1
                                     !SO4      BC        OC      SS       DU
  flag_aero_comp_modes = reshape ( (/ .TRUE. , .FALSE., .FALSE., .FALSE., .FALSE. , & !NS
                                      .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE. , & !KS
                                      .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE.  , & !AS
                                      .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE.  , & !CS
                                      .FALSE., .TRUE. , .TRUE. , .FALSE., .FALSE. , & !KI
                                      .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.  , & !AI
                                      .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.  /), (/NAERO_SPEC, nclass/) )  !CI
  bc_aerosol%ef_actual_unit = 'kg kg-1'
  DO imod = 1, nclass
    bc_aerosol%ef_varname = "NUM_"//TRIM(sizeclass(imod)%shortname)
    ibc_numden_modes(imod)=bc_define('number mixing ratio: ' &
                                     //TRIM(sizeclass(imod)%classname), bc_aerosol, 3, .TRUE.)
    DO jt = 1, NAERO_SPEC
      bc_aerosol%ef_varname = TRIM(aero_spec(jt))//'_'//TRIM(sizeclass(imod)%shortname)
      IF (flag_aero_comp_modes(jt,imod)) THEN
        ibc_aero_comp_modes(imod,jt)=bc_define(TRIM(aero_spec(jt))//' mass mixing ratio: ' &
                                               //TRIM(sizeclass(imod)%classname), bc_aerosol, 3, .TRUE.)
      ENDIF
    ENDDO
  ENDDO



  END SUBROUTINE het_tropo_m7_boundary_condition_inti

  !!----------------------------------------------------------------------------
  SUBROUTINE het_tropo_salsa_boundary_condition_inti()
  USE mo_tracer,                   ONLY: get_tracer
  USE mo_ham,                      ONLY: sizeclass
  USE mo_boundary_condition,       ONLY: bc_define

  ! initialisation routine

  !--- local variables
  INTEGER                    :: imod,jt,ierr,idx
  LOGICAL                    :: flag_aero_comp_modes(NAERO_SPEC,nclass)

  bc_aerosol%ef_actual_unit = 'm'
  DO imod = 1, nclass
      bc_aerosol%ef_varname = "rwet_"//TRIM(sizeclass(imod)%shortname)
      ibc_rwet(imod) = bc_define('aerosol wet radius: '//TRIM(sizeclass(imod)%classname), bc_aerosol, 3, .TRUE.)
  ENDDO

  ibc_aero_comp_modes = -1
                                     ! SO4      BC        OC      SS       DU
  flag_aero_comp_modes = reshape ( (/ .TRUE. , .FALSE., .TRUE. , .FALSE., .FALSE. , & !1a1
                                      .TRUE. , .FALSE., .TRUE. , .FALSE., .FALSE. , & !1a2
                                      .TRUE. , .FALSE., .TRUE. , .FALSE., .FALSE. , & !1a3
                                      .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE.  , & !2a1
                                      .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE.  , & !2a2
                                      .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE.  , & !2a3
                                      .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE.  , & !2a4
                                      .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE.  , & !2a5
                                      .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE.  , & !2a6
                                      .TRUE. , .TRUE. , .TRUE. , .TRUE. , .TRUE.  , & !2a7
                                      .TRUE. , .TRUE. , .TRUE. , .FALSE., .TRUE.  , & !2b1
                                      .TRUE. , .TRUE. , .TRUE. , .FALSE., .TRUE.  , & !2b2
                                      .TRUE. , .TRUE. , .TRUE. , .FALSE., .TRUE.  , & !2b3
                                      .TRUE. , .TRUE. , .TRUE. , .FALSE., .TRUE.  , & !2b4
                                      .TRUE. , .TRUE. , .TRUE. , .FALSE., .TRUE.  , & !2b5
                                      .TRUE. , .TRUE. , .TRUE. , .FALSE., .TRUE.  , & !2b6
                                      .TRUE. , .TRUE. , .TRUE. , .FALSE., .TRUE./), (/NAERO_SPEC, nclass/) )  !2b7
  bc_aerosol%ef_actual_unit = 'kg kg-1'
  DO imod = 1, nclass
    bc_aerosol%ef_varname = "NUM_"//TRIM(sizeclass(imod)%shortname)
    ibc_numden_modes(imod)=bc_define('number mixing ratio: ' &
                                     //TRIM(sizeclass(imod)%classname), bc_aerosol, 3, .TRUE.)
    DO jt = 1, NAERO_SPEC
      bc_aerosol%ef_varname = TRIM(aero_spec(jt))//'_'//TRIM(sizeclass(imod)%shortname)
      IF (flag_aero_comp_modes(jt,imod)) THEN
        ibc_aero_comp_modes(imod,jt)=bc_define(TRIM(aero_spec(jt))//' massmixing ratio: ' &
                                               //TRIM(sizeclass(imod)%classname), bc_aerosol, 3, .TRUE.)
      ENDIF
    ENDDO
  ENDDO

  END SUBROUTINE het_tropo_salsa_boundary_condition_inti

  !------------------------------------------------------------------------------------------------
  SUBROUTINE ratecon_sftropo_inti()

  ! define boundary conditions for aerosol properties
  ! identify (gas-phase) tracers (NO3, N2O5) and set constants for diffusion

  USE mo_external_field_processor, ONLY: EF_FILE, EF_MODULE
  USE mo_exception,                ONLY: finish

  !--- local variables
  INTEGER      :: ii, ierr, idum, inml, iunit, jt

  !
  ! 1) define boundary conditions
  !

  IF ((bc_aerosol%ef_type == EF_MODULE) .OR.  (bc_aerosol%ef_type == EF_FILE)) THEN
    SELECT CASE(nham_subm)

       CASE(HAM_M7) 
         CALL het_tropo_m7_boundary_condition_inti()

       CASE(HAM_SALSA)
         CALL het_tropo_salsa_boundary_condition_inti()

    END SELECT 
  ELSE
    CALL finish("ratecon_sftropo_inti","only works for aerosol boundary condition of type MODULE or FILE")
  ENDIF

  !
  ! 2) identify (gas-phase) tracers
  !

  ! Don't use argument ierr=... ==> finish is called if tracers are not found
  CALL get_tracer('NO3', idx=idt(idx_no3),ierr=ierr)
  IF (ierr .ne. 0) CALL get_tracer('NO3_moz', idx=idt(idx_no3))
  CALL get_tracer('N2O5',idx=idt(idx_n2o5))
  CALL get_tracer('HO2',idx=idt(idx_ho2))
  CALL get_tracer('HNO3',idx=idt(idx_hno3))
  CALL get_tracer('O3',idx=idt(idx_o3)) 
  CALL get_tracer('NO2', idx=idt(idx_no2),ierr=ierr)

  !
  ! 3) constants for diffusion
  !

  ! mmspd= (8RT/pi*PM)^1/2; PM=[g/mol]; mmspd=[cm/s]
  const=1.e2_dp*(8._dp*argas/pi)**0.5_dp
  dg_const1=0.0185_dp*1.e20_dp/avo 
  ![3/(8*avo*dq^2)]; A typical value of dg is 4.5*1.e-10 [m] (M.Z.Jacobson,pag.456)
  dg_const2=(argas*zmolgair*1.e-03_dp/2._dp*2._dp*pi)  ![R*PMair/2pi]

  END SUBROUTINE ratecon_sftropo_inti

  !!----------------------------------------------------------------------------
  SUBROUTINE ratecon_sftropo(temp, prho, h2ovmr, pmid, m, lat, kproma, &
                          zk_het_no3, zk_het_no2, zk_het_n2o5, &
                          zk_het_ho2, zk_het_hno3, zk_het_o3)  
  USE mo_exception,     ONLY : message_text, finish
  USE mo_moz_mods,      ONLY : plonl, plev, pcnstm1
  USE mo_moz_diag,      ONLY : dpk_no3,dpk_no2, dpk_n2o5, dpk_ho2, dpk_hno3, dpk_o3

  ! wrapper routine for calculation of khet(s)
  ! wrapper builded for different possibilities of microphysics (lm7, lsalsa, ...)

  !--- arguments
  REAL(dp), INTENT(IN) :: temp(plonl,plev), &           ! temperature (K)
                          prho(plonl, plev), &          ! air density
                          m(plonl,plev), &              ! total atm density (1/cm^3)
                          h2ovmr(plonl,plev), &         ! water vapor (mol/mol)
                          pmid(plonl,plev)              ! midpoint pressure (Pa)
  INTEGER, INTENT(IN)  :: lat, kproma                   ! latitude index, kproma
  REAL(dp),INTENT(OUT) :: zk_het_no3(plonl,plev), &
                          zk_het_no2(plonl,plev), &
                          zk_het_n2o5(plonl,plev),&
                          zk_het_ho2(plonl,plev), &
                          zk_het_hno3(plonl,plev),&
                          zk_het_o3(plonl,plev)
 
  !--- local variables
  REAL(dp)             :: mmspd(plonl,plev), &          ! mean molecular speed
                          dg(plonl,plev), &             ! gas-phase molecular diffusion coefficient [cm2 s-1]
                          khet(plonl,plev), &           ! local variable for calculation of 
                          gamma(plonl,plev)             ! heterogeneous reaction rate coefficients
  INTEGER              :: jt


  DO jt = 1,MAXTRAC
    CALL calc_dg(plonl, plev, temp, prho, jt, mmspd, dg)

    SELECT CASE(nham_subm)

      CASE(HAM_M7)  
        CALL calc_khet_m7(temp, jt, mmspd, dg, lat, kproma, prho, khet) 

      CASE(HAM_SALSA)
        CALL calc_khet_salsa(temp, jt, mmspd, dg, lat, kproma, prho, khet)

    END SELECT

      IF (jt == idx_no3) THEN
        zk_het_no3 = khet
        dpk_no3(1:kproma,:,lat) = khet(1:kproma,:)
      ELSE IF (jt == idx_n2o5) THEN
        zk_het_n2o5 = khet
        dpk_n2o5(1:kproma,:,lat) = khet(1:kproma,:)
      ELSE IF (jt == idx_ho2) THEN
        zk_het_ho2 = khet
        dpk_ho2(1:kproma,:,lat) = khet(1:kproma,:)
      ELSE IF (jt == idx_hno3) THEN
        zk_het_hno3 = khet
        dpk_hno3(1:kproma,:,lat) = khet(1:kproma,:)
      ELSE IF (jt == idx_o3) THEN
        zk_het_o3 = khet
        dpk_o3(1:kproma,:,lat) = khet(1:kproma,:)
      ELSE IF (jt == idx_no2) THEN
        zk_het_no2 = khet
        dpk_no2(1:kproma,:,lat) = khet(1:kproma,:)
      ELSE
        WRITE(message_text,*) "ERROR jt =", jt
        CALL finish('ratecon_sftropo', message_text)
      END IF
  ENDDO
  END SUBROUTINE ratecon_sftropo

  !!----------------------------------------------------------------------------
  SUBROUTINE calc_dg(plonl, plev, ptm1, prhop1, idx, mmspd, dg)

  ! calculation of diffusion coefficient
  
  !--- arguments
  INTEGER, INTENT(in)         :: plonl,plev, idx          ! idx is the tracer index eg. idx_h2o2
  REAL(dp), INTENT(in)        :: ptm1(plonl,plev), prhop1(plonl,plev)! temperature and density
  REAL(dp), INTENT(out)       :: mmspd(plonl,plev), &     ! mean molecular speed
                                 dg(plonl,plev)           ! gas-phase molecular diffusion coefficient [cm2 s-1]

  ! obtain tracer moleweights
  mweight(idx) = trlist%ti(idt(idx))%moleweight
  pmpma(idx) = 1._dp + zmolgair/mweight(idx)

  ! Dg = [cm^2/s]; Fundamentals of atmospheric modelling, M.Z.Jacobson,pag.456.
  mmspd(1:plonl,1:plev) = const * (1.e03_dp*ptm1(1:plonl,1:plev)/mweight(idx))**0.5_dp
  dg(1:plonl,1:plev)    = 1.e04_dp*dg_const1/prhop1(1:plonl,1:plev) &
                          *(dg_const2*ptm1(1:plonl,1:plev)*pmpma(idx))**0.5_dp

  END SUBROUTINE calc_dg

  !!----------------------------------------------------------------------------
  SUBROUTINE calc_gamma(plonl, plev, kproma, temp, radius, mmspd, idx, lat, zmass, gamma)

  USE mo_exception,     ONLY : message_text, finish
  USE mo_memory_g3b,    ONLY: relhum
  !USE mo_debugs,        ONLY: ddf01,ddf02,ddf03,ddf04,ddf05

  ! determine reaction probabilities

  !--- arguments
  INTEGER, INTENT(IN)  :: idx, plonl, plev, lat, kproma 
  REAL(dp), INTENT(IN) :: temp(plonl,plev), zmass(plonl,plev,NAERO_SPEC), mmspd(plonl,plev), radius(plonl,plev)
  REAL(dp), INTENT(OUT):: gamma(plonl,plev)

  !--- local variables
  REAL(dp), PARAMETER  :: c1= 2.79e-04_dp, &       ! N2O5 parameterization
                          c2= 1.30e-04_dp, &
                          c3=-3.43e-06_dp, &
                          c4= 7.52e-08_dp, &
                          c_HO2_gas = 1.0e8_dp, &     ! Assumed gas phase concentration of HO2 [molec cm-3]
                          R_HO2  = 0.082057_dp        ! Universal gas constant [atmLmol-1k-1]

  REAL(dp)             :: beta, alpha, gamma_n2o5(plonl,plev,NAERO_SPEC),gamma_ho2(plonl,plev,NAERO_SPEC), &
                          gamma_hno3(plonl,plev,NAERO_SPEC),gamma_o3(plonl,plev,NAERO_SPEC), ztotal_mass(plonl,plev), relhum_dust, &
                          var, r_p(plonl,plev), pH, H_HO2_eff, k_eff
  INTEGER              :: jk, jl, jt, ibc


  ztotal_mass = 0.0_dp
  DO jt = 1,NAERO_SPEC
    DO jk=1,plev
      DO jl=1,kproma
        IF (zmass(jl,jk,jt) > 0.0_dp) ztotal_mass(jl,jk) = ztotal_mass(jl,jk) + zmass(jl,jk,jt)
      ENDDO
    ENDDO
  ENDDO
 
  IF(idx == idx_no3) THEN
    gamma(:,:) = 0.001_dp 
  ELSE IF(idx == idx_no2) THEN
    gamma(:,:) = 1.e-4_dp  
  ELSE IF (idx == idx_n2o5) THEN

    ! Gammas calculated like in Evans 2005
    ! Sulfate
  
    ! [Kane,2001;Hallquist,2003]
    DO jk=1,plev
      DO jl=1,kproma ! plonl
        IF (temp(jl,jk) >= 282._dp) THEN
          beta=-4.e-02_dp*(temp(jl,jk)-294._dp)
        ELSE
          beta=0.48_dp
        ENDIF
          gamma_n2o5(jl,jk,idx_aero_so4)=(10._dp**beta)*(c1+c2*relhum(jl,jk,lat)*100&
                                         +c3*(100*relhum(jl,jk,lat)**2)+c4*(100*relhum(jl,jk,lat)**3))
      ENDDO
    ENDDO

    ! Organic Carbon
    ! [Thornton,2003]
    DO jk=1,plev
      DO jl=1,kproma
        IF (relhum(jl,jk,lat)*100 >= 57._dp) THEN
          gamma_n2o5(jl,jk,idx_aero_oc)= 0.03_dp
        ELSE
          gamma_n2o5(jl,jk,idx_aero_oc)= 5.2e-04_dp*relhum(jl,jk,lat)*100
        ENDIF
      ENDDO
    ENDDO

    ! Black Carbon
    ! [Sander,2003]
    gamma_n2o5(:,:,idx_aero_bc)= 0.005_dp


    ! Sea Salt
    ! [Atkinson,2004]
    DO jk=1,plev
      DO jl=1,kproma
        IF (relhum(jl,jk,lat)*100 >= 50._dp) THEN
          gamma_n2o5(jl,jk,idx_aero_ss)= 0.03_dp
        ELSE
          gamma_n2o5(jl,jk,idx_aero_ss)= 0.005_dp
        ENDIF
      ENDDO
    ENDDO

    ! Dust
    ! [Bauer,2004]

    DO jk=1,plev
      DO jl=1,kproma
        relhum_dust = min(70.0_dp,max(30.0_dp,relhum(jl,jk,lat)*100)) ! only valid in this range of relhum
        gamma_n2o5(jl,jk,idx_aero_du) = 4.25e-4_dp*relhum_dust - 9.75e-3_dp
      ENDDO
    ENDDO

    ! Calculation of a mean gamma_n2o5 per mode over all aerosol components weighted
    ! by the relative contribution of each to the total mass.
   
    gamma = 0._dp 
    DO jt = 1, NAERO_SPEC
      DO jk=1,plev
        DO jl=1,kproma
          IF (ztotal_mass(jl,jk) > 0.0_dp .AND. zmass(jl,jk,jt) > 0.0_dp) THEN
            gamma(jl,jk) = gamma(jl,jk) + zmass(jl,jk,jt)/ztotal_mass(jl,jk)*gamma_n2o5(jl,jk,jt)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  ELSE IF (idx == idx_ho2) THEN
!-------------------------------HO2--------------------------------------!
    !gamma(:,:) = 0.2_dp 

!!    ! Gammas calculated like in Macintyre and Evans 2011
!!    beta  = 0.0_dp
!!    alpha = 0.0_dp
!!    ! Sulfate
!!    DO jk=1,plev
!!      DO jl=1,kproma ! plonl
!!        IF (temp(jl,jk) >= 282._dp) THEN
!!          IF (relhum(jl,jk,lat)*100 >= 35._dp) THEN
!!            IF (temp(jl,jk) > 0.0_dp) alpha = 5.14545e-4_dp*exp(1560._dp/temp(jl,jk))
!!            beta  = -26.1818_dp*exp(-0.078_dp*relhum(jl,jk,lat)*100) + 1.74545_dp
!!            gamma_ho2(jl,jk,idx_aero_so4) = max(0._dp,alpha*beta)
!!          ELSE
!!            gamma_ho2(jl,jk,idx_aero_so4) = 0.01_dp
!!          ENDIF
!!        ELSE
!!          IF (relhum(jl,jk,lat)*100 >= 35._dp) THEN
!!            alpha = 5.14545e-4_dp*exp(1560._dp/282._dp)
!!            beta  = -26.1818_dp*exp(-0.078_dp*relhum(jl,jk,lat)*100) + 1.74545_dp
!!            gamma_ho2(jl,jk,idx_aero_so4) = max(0._dp,alpha*beta)
!!          ELSE
!!            gamma_ho2(jl,jk,idx_aero_so4) = 0.01_dp
!!          ENDIF
!!        ENDIF
!!      ENDDO
!!    ENDDO
!!
!!    ! Organic Carbon
!!    gamma_ho2(:,:,idx_aero_oc) = 0.025_dp
!!
!!    ! Black Carbon
!!    gamma_ho2(:,:,idx_aero_bc)= 0.01_dp
!!
!!
!!    ! Sea Salt
!!    gamma_ho2(:,:,idx_aero_ss)= 0.05_dp
!!
!!    ! Dust
!!    DO jk=1,plev
!!      DO jl=1,kproma
!!        IF (relhum(jl,jk,lat)*100 >= 50._dp) THEN
!!          gamma_ho2(jl,jk,idx_aero_du) = 0.1_dp
!!        ELSE
!!          gamma_ho2(jl,jk,idx_aero_du) = 0.05_dp
!!        END IF
!!      ENDDO
!!    ENDDO
 
    gamma = 0.0_dp
    H_HO2_eff = 0.0_dp
    pH = 1.e-5_dp ! Assume a pH value of 5
    k_eff = (8.6e-5_dp + (2.1e-5_dp/pH)*1.e8_dp)/(1. + (2.1e-5_dp/pH)) ![1/(M s)]

    ! reaction probability on cloud droplets Thornton et. al 2008
    !        1       1                3 N_A mmspd
    !     ------ = ----- + ----------------------------------------
    !      gamma   alpha   8000 (H_eff*R_g*T)^2 k_eff [HO2]_gas r_p
    
    !       1        1            mmspd
    !<=>  ----- =  ----- + const -------
    !     gamma    alpha           r_p

    ! calculate constat out of constant values
    
    ! alpha_HO2 = 1 --> 1/alpha = 1 so the variable alpha will not be used in this case 
     
    r_p(1:kproma,:) = radius(1:kproma,:)*1.0e6_dp ! mean volume radius in micro meters
    DO jk=1,plev
      DO jl=1,kproma
        H_HO2_eff = 9.4e-6_dp*exp(5.92e3_dp/temp(jl,jk))*(1._dp + 2.1e-5_dp/pH) 
        var = 3.e5_dp*avo/(8000.0_dp*((H_HO2_eff*R_HO2*temp(jl,jk))**2.0_dp)*k_eff*c_HO2_gas) ! Multiplicaton by 10**5 due to unit conversion
        IF (r_p(jl,jk) > 0.0_dp) THEN
          gamma(jl,jk) = (1.0_dp + var*mmspd(jl,jk)/r_p(jl,jk))**(-1.0_dp)
        IF (gamma(jl,jk) < 0.0_dp) gamma(jl,jk) = 0.0_dp
        END IF
      ENDDO
    ENDDO
    ! Sulfate
    gamma_ho2(:,:,idx_aero_so4) = gamma(:,:)

    ! Organic Carbon
    gamma_ho2(:,:,idx_aero_oc) = gamma(:,:)

    ! Black Carbon
    gamma_ho2(:,:,idx_aero_bc)= gamma(:,:) 

    ! Sea Salt
    gamma_ho2(:,:,idx_aero_ss)= gamma(:,:)

    ! Dust
    gamma_ho2(:,:,idx_aero_du) = 0.1_dp

    ! Calculation of a mean gamma_ho2 per mode over all aerosol components weighted
    ! by the relative contribution of each to the total mass.
   
    gamma = 0._dp 
    DO jt = 1, NAERO_SPEC
      DO jk=1,plev
        DO jl=1,kproma
          !IF (ztotal_mass(jl,jk) > 0) THEN
          IF (ztotal_mass(jl,jk) > 0.0_dp .AND. zmass(jl,jk,jt) > 0.0_dp) THEN
            gamma(jl,jk) = gamma(jl,jk) + zmass(jl,jk,jt)/ztotal_mass(jl,jk)*gamma_ho2(jl,jk,jt)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
   
!-------------------------------HNO3-------------------------------------!
  ELSE IF (idx == idx_hno3) THEN
    ! Sulfate
    gamma_hno3(:,:,idx_aero_so4) = 0.0_dp

    ! Organic Carbon
    gamma_hno3(:,:,idx_aero_oc) = 0.0_dp

    ! Black Carbon
    gamma_hno3(:,:,idx_aero_bc)= 0.0_dp

    ! Sea Salt
    ! [Davies 1998]
    gamma_hno3(:,:,idx_aero_ss)= 0.01_dp

    ! Dust
    ! [Hodzic 2006]
    gamma_hno3(:,:,idx_aero_du) = 0.1_dp
 
    ! Calculation of a mean gamma_hno3 per mode over all aerosol components weighted
    ! by the relative contribution of each to the total mass.
   
    gamma = 0._dp 
    DO jt = 1, NAERO_SPEC
      DO jk=1,plev
        DO jl=1,kproma
          !IF (ztotal_mass(jl,jk) > 0) THEN
          IF (ztotal_mass(jl,jk) > 0.0_dp .AND. zmass(jl,jk,jt) > 0.0_dp) THEN
            gamma(jl,jk) = gamma(jl,jk) + zmass(jl,jk,jt)/ztotal_mass(jl,jk)*gamma_hno3(jl,jk,jt)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
   
!---------------------------------O3-------------------------------------!
  ELSE IF (idx == idx_o3) THEN
    ! Sulfate
    gamma_o3(:,:,idx_aero_so4) = 0.0_dp

    ! Organic Carbon
    gamma_o3(:,:,idx_aero_oc) = 0.0_dp

    ! Black Carbon
    gamma_o3(:,:,idx_aero_bc)= 0.0_dp

    ! Sea Salt
    gamma_o3(:,:,idx_aero_ss)= 0.0_dp

    ! Dust
    gamma_o3(:,:,idx_aero_du) = 1.0e-6_dp
 
    ! Calculation of a mean gamma_o3 per mode over all aerosol components weighted
    ! by the relative contribution of each to the total mass.
   
    gamma = 0._dp 
    DO jt = 1, NAERO_SPEC
      DO jk=1,plev
        DO jl=1,kproma
          IF (ztotal_mass(jl,jk) > 0.0_dp .AND. zmass(jl,jk,jt) > 0.0_dp) THEN
            gamma(jl,jk) = gamma(jl,jk) + zmass(jl,jk,jt)/ztotal_mass(jl,jk)*gamma_o3(jl,jk,jt)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
   
  ELSE
    WRITE(message_text,*) "ERROR idx =",idx 
    CALL finish('calc_gamma', message_text)
  END IF


  END SUBROUTINE calc_gamma


  !!----------------------------------------------------------------------------

  SUBROUTINE calc_gamma_cloud(plonl, plev, kproma, temp, idx, mmspd, zmvolr, lat, gamma)

  USE mo_exception,     ONLY : message_text, finish
  USE mo_memory_g3b,    ONLY: relhum
  !USE mo_debugs,        ONLY: ddf01,ddf02,ddf03,ddf04,ddf05

  ! determine reaction probabilities

  !--- arguments
  INTEGER, INTENT(IN)  :: idx, plonl, plev, lat, kproma
  REAL(dp), INTENT(IN) :: mmspd(plonl,plev), zmvolr(plonl,plev),temp(plonl,plev) 
  REAL(dp), INTENT(OUT):: gamma(plonl,plev)

  !--- local variables --- Values taken from Whalley et. al 2015 ---
  REAL(dp), PARAMETER  :: c_HO2_gas = 1.0e8_dp, &     ! Assumed gas phase concentration of HO2 [molec cm-3]
                          R_HO2  = 0.082057_dp, &     ! Universal gas constant [atmLmol-1k-1]
                          alpha_HO2 = 1.0_dp          ! Accommodation coefficient HO2

  REAL(dp)            :: var, r_p(plonl,plev), pH, H_HO2_eff, k_eff              
  INTEGER              :: jk, jl

  ! Initialization
  gamma = 0.0_dp
  H_HO2_eff = 0.0_dp
  pH = 1.e-5_dp ! Assume a pH value of 5
  k_eff = (8.6e-5_dp + (2.1e-5_dp/pH)*1.e8_dp)/(1. + (2.1e-5_dp/pH)) ![1/(M s)]

    ! reaction probability on cloud droplets Thornton et. al 2008
    !        1       1                3 N_A mmspd
    !     ------ = ----- + ----------------------------------------
    !      gamma   alpha   8000 (H_eff*R_g*T)^2 k_eff [HO2]_gas r_p
    
    !       1        1            mmspd
    !<=>  ----- =  ----- + const -------
    !     gamma    alpha           r_p

    ! calculate constat out of constant values
    
    ! alpha_HO2 = 1 --> 1/alpha = 1 so the variable alpha will not be used in this case 
     
    r_p(1:kproma,:) = zmvolr(1:kproma,:)*1.0e6_dp ! mean volume radius in micro meters
    DO jk=1,plev
      DO jl=1,kproma
        H_HO2_eff = 9.4e-6_dp*exp(5.92e3_dp/temp(jl,jk))*(1._dp + 2.1e-5_dp/pH) 
        var = 3.e5_dp*avo/(8000.0_dp*((H_HO2_eff*R_HO2*temp(jl,jk))**2.0_dp)*k_eff*c_HO2_gas) ! Multiplicaton by 10**5 due to unit conversion
        IF (r_p(jl,jk) > 0.0_dp) THEN
          gamma(jl,jk) = (1.0_dp + var*mmspd(jl,jk)/r_p(jl,jk))**(-1.0_dp)
        IF (gamma(jl,jk) < 0.0_dp) gamma(jl,jk) = 0.0_dp
        END IF
      ENDDO
    ENDDO
   
  END SUBROUTINE calc_gamma_cloud


  !!----------------------------------------------------------------------------
  SUBROUTINE calc_khet_m7(temp, idx, mmspd, dg,  &
                         lat, kproma, prhop1, khet)
  
  USE mo_exception,          ONLY: message, message_text, em_info  
  USE mo_math_constants,     ONLY: pi
  USE mo_species,            ONLY: speclist, aero_idx
  USE mo_ham_species,        ONLY: ham_naerospec
  USE mo_ham_m7ctl,          ONLY: cmr2ras
  USE mo_ham,                ONLY: sizeclass
  USE mo_moz_mods,           ONLY: plonl, plev, pcnstm1
  USE mo_moz_diag,           ONLY: dpsa_aerosol, dpsa_aerosol_wet, dpsa_cloud, &
                                   dpgamma_n2o5, dpgamma_n2o5_a, &
                                   dpgamma_ho2, dpgamma_ho2_a, dpgamma_ho2_cloud, dpgamma_hno3, &
                                   dpgamma_hno3_a, dpgamma_o3, dpgamma_o3_a, dpgamma_ho2_tot
                                    
  USE mo_boundary_condition, ONLY: bc_apply, bc_find
  USE mo_vphysc,             ONLY: vphysc

  ! calculation of khet(s) for m7 microphysics

  !--- arguments
  INTEGER, INTENT(IN)         :: lat,idx, kproma
  REAL(dp),INTENT(IN)         :: mmspd(plonl,plev), dg(plonl,plev), &
                                 prhop1(plonl,plev), temp(plonl,plev) 
  REAL(dp),INTENT(OUT)        :: khet(plonl,plev)


  !--- local variables
  INTEGER              :: jclass, jk, jl, ibc, ierr, jt
  REAL(dp)             :: gamma(plonl,plev),            &
                          radius(plonl,plev),           & ! dry or wet radius [m]
                          zmvolr(plonl,plev),           & ! cloud volume mean radius [m]
                          za,                           & ! dry or wet surface area density
                          za_tot(plonl,plev),           & ! dry or wet surface total area denstiy
                          za_cloud_tot(plonl,plev),     & ! cloud droplet surface total area denstiy
                          khet_cloud(plonl,plev),       & ! reaction rate coefficient HO2 on cloud [s-1]
                          khet_aero(plonl,plev),        & ! reaction rate coefficient on aerosol [s-1]
                          gamma_tot(plonl,plev),        & ! Sum over reactionprobabilities for single species
                          numden(plonl,plev),           & ! number density [1/kg]
                          zcdnc(plonl,plev),            & ! cloud droplet number concentration [1/m3]
                          zclw(plonl,plev),             & ! cloud liquid water [kg/kg]
                          zaclc(plonl,plev),            & ! cloud fraction [1]
                          zmass(plonl,plev,NAERO_SPEC), & ! [kg/kg] 
                          taudg(plonl,plev),            & ! only variables to calculate khet, but
                          tauim(plonl,plev)               ! named as characteristic times they can be interpreted
                                                          ! taudg = diffusion life time 
                                                          ! tauim = interfaical mass transport life time

  LOGICAL              :: lhydrolysis                     ! switch for hydrolysis heterogeneous reactions
                                                          ! only happening on aqueous aerosol
  ! Initialization
  za_tot = 0.0_dp
  za_cloud_tot = 0.0_dp
  gamma_tot = 0.0_dp
  khet = 0._dp
  khet_aero = 0._dp
  khet_cloud = 0._dp
  zcdnc(:,:) = 0.0_dp 
  zclw(:,:)  = 0.0_dp  
  zaclc(:,:) = 0.0_dp
  lhydrolysis = .FALSE.
 
  ! Calculate khet
  !                A 
  ! khet = --------------------
  !        r/dg + 4/gamma mmspd
  ! Chang 2011

  IF (idx == idx_no3)  lhydrolysis = .TRUE.
  IF (idx == idx_no2)  lhydrolysis = .TRUE.
  IF (idx == idx_ho2)  lhydrolysis = .TRUE.
  IF (idx == idx_hno3) lhydrolysis = .TRUE.

  DO jclass=1, nclass
    zmass = 0.0_dp
    DO jt=1,NAERO_SPEC
      ibc = ibc_aero_comp_modes(jclass,jt)
      IF (ibc > 0) CALL bc_apply(ibc,kproma,lat,zmass(:,:,jt))
    END DO
    CALL bc_apply(ibc_numden_modes(jclass),kproma,lat,numden(:,:))
    IF(lhydrolysis .AND. .NOT. sizeclass(jclass)%lsoluble) CONTINUE ! Falls lhydrolysis=False, macht er alle jclass
    CALL bc_apply(ibc_rwet(jclass),kproma,lat,radius)
    CALL calc_gamma(plonl, plev, kproma, temp, radius, mmspd, idx, lat, zmass, gamma)
    DO jk=1, plev
      DO jl=1,kproma !plonl
        za= 4.e-02_dp*pi*numden(jl,jk)*prhop1(jl,jk)*(radius(jl,jk) & ! [numden] = [#/kg]  
                       *cmr2ras(jclass))**2._dp ![cm2/cm3]
        IF (za < 0._dp) za = 0._dp
        za_tot(jl,jk) = za_tot(jl,jk) + za ! just for diagnostics
        IF(gamma(jl,jk) > 1.e-12_dp) THEN
          taudg = radius(jl,jk)/dg(jl,jk)
          tauim = 4/(gamma(jl,jk)*mmspd(jl,jk))
          khet_aero(jl,jk) = khet_aero(jl,jk) + za/(taudg(jl,jk) + tauim(jl,jk))
          gamma_tot(jl,jk) = gamma_tot(jl,jk) + gamma(jl,jk)*za ! just for diagnostics
        ENDIF
      ENDDO
    ENDDO
  END DO

 
  IF (idx == idx_no3 ) khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_no2 ) khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_n2o5) khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_hno3) khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_o3)   khet(1:kproma,:) = khet_aero(1:kproma,:)
  !IF (idx == idx_ho2)  khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_ho2) THEN
    zcdnc(1:kproma,:) = vphysc%cdnc(1:kproma,:,lat) 
    zclw(1:kproma,:)  = vphysc%clw(1:kproma,:,lat)  
    zaclc(1:kproma,:) = vphysc%aclc(1:kproma,:,lat) 

    ! Calculation mean volume radius as in mo_cloud_micro_2m.f90
    zmvolr(1:kproma,:) = (0.75_dp*MAX(zclw(1:kproma,:),0._dp) * prhop1(1:kproma,:) &      ! zclw = zxlb
            / (pi*rhoh2o*zcdnc(1:kproma,:)))**(1._dp/3._dp)
    
    CALL calc_gamma_cloud(plonl, plev, kproma, temp, idx, mmspd, zmvolr, lat, gamma)
    za = 0.0_dp
    DO jk=1, plev
      DO jl=1,kproma !plonl
        za = 4.e-02_dp*pi*zmvolr(jl,jk)**2._dp*zcdnc(jl,jk)*0.01_dp    ![cm2/cm3]
        IF (za < 0._dp) za = 0._dp ! set to 0, because zcdnc or zmvolr can be < 0
        IF (zmvolr(jl,jk) < 0._dp) zmvolr(jl,jk) = 0._dp
        za_cloud_tot(jl,jk) = za_cloud_tot(jl,jk) + za ! just for diagnostics
        IF(gamma(jl,jk) > 1.e-12_dp) THEN
          taudg = zmvolr(jl,jk)/dg(jl,jk)
          tauim = 4/(gamma(jl,jk)*mmspd(jl,jk))
          khet_cloud(jl,jk) = za/(taudg(jl,jk) + tauim(jl,jk))
        END IF
      ENDDO
    ENDDO
                                                         
    DO jk=1, plev
      DO jl=1,kproma !plonl 
        khet(jl,jk) = (1.0_dp-zaclc(jl,jk))*khet_aero(jl,jk) + zaclc(jl,jk)*khet_cloud(jl,jk)
      END DO
    END DO
  END IF 


  ! diagnostics (see ldetail2 in mo_moz_diag.f90)
  IF (idx == idx_n2o5) THEN
    dpsa_aerosol(1:kproma,:,lat) = za_tot(1:kproma,:)
    dpgamma_n2o5(1:kproma,:,lat) = 0._dp
    dpgamma_n2o5_a(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_tot(jl,jk) > 0._dp) THEN
          dpgamma_n2o5(jl,jk,lat) = gamma_tot(jl,jk)/za_tot(jl,jk) ! gamma weighted with modes and surfaces
          dpgamma_n2o5_a(jl,jk,lat) = 4.0_dp*khet(jl,jk)/(za_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO

  ELSE IF (idx == idx_ho2) THEN
    dpsa_aerosol_wet(1:kproma,:,lat) = za_tot(1:kproma,:)
    dpsa_cloud(1:kproma,:,lat) = za_cloud_tot(1:kproma,:)
    dpsa_cloud(1:kproma,:,lat) = 0.0_dp
    dpgamma_ho2(1:kproma,:,lat) = 0._dp
    dpgamma_ho2_a(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_tot(jl,jk) > 0._dp) THEN
          dpgamma_ho2(jl,jk,lat) = gamma_tot(jl,jk)/za_tot(jl,jk) ! gamma weighted with modes and surfaces
          dpgamma_ho2_a(jl,jk,lat) = 4.0_dp*khet_aero(jl,jk)/(za_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO

    dpgamma_ho2_cloud(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_cloud_tot(jl,jk) > 0._dp) THEN
          dpgamma_ho2_cloud(jl,jk,lat) = 4.0_dp*khet_cloud(jl,jk)/(za_cloud_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO

    dpgamma_ho2_tot(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_cloud_tot(jl,jk) > 0._dp .and. za_tot(jl,jk) > 0._dp) THEN
          dpgamma_ho2_tot(jl,jk,lat) = 4.0_dp*khet(jl,jk)/((za_cloud_tot(jl,jk)+za_tot(jl,jk))*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO

  ELSE IF (idx == idx_hno3) THEN
    dpgamma_hno3(1:kproma,:,lat) = 0._dp
    dpgamma_hno3_a(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_tot(jl,jk) > 0._dp) THEN
          dpgamma_hno3(jl,jk,lat) = gamma_tot(jl,jk)/za_tot(jl,jk) ! gamma weighted with modes and surfaces
          dpgamma_hno3_a(jl,jk,lat) = 4.0_dp*khet(jl,jk)/(za_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO

  ELSE IF (idx == idx_o3) THEN
    dpgamma_o3(1:kproma,:,lat) = 0._dp
    dpgamma_o3_a(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_tot(jl,jk) > 0._dp) THEN
          dpgamma_o3(jl,jk,lat) = gamma_tot(jl,jk)/za_tot(jl,jk) ! gamma weighted with modes and surfaces
          dpgamma_o3_a(jl,jk,lat) = 4.0_dp*khet(jl,jk)/(za_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  END SUBROUTINE calc_khet_m7
!----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_khet_salsa(temp, idx, mmspd, dg,  &
                         lat, kproma, prhop1, khet)

  USE mo_exception,          ONLY: message, message_text, em_info
  USE mo_math_constants,     ONLY: pi
  USE mo_species,            ONLY: speclist, aero_idx
  USE mo_ham_species,        ONLY: ham_naerospec
  USE mo_ham,                ONLY: sizeclass, nsol
  USE mo_moz_mods,           ONLY: plonl, plev, pcnstm1
  USE mo_moz_diag,           ONLY: dpsa_aerosol, dpsa_aerosol_wet, dpsa_cloud, &
                                   dpgamma_n2o5, dpgamma_n2o5_a, &
                                   dpgamma_ho2, dpgamma_ho2_a, dpgamma_ho2_cloud, dpgamma_hno3, &
                                   dpgamma_hno3_a, dpgamma_o3, dpgamma_o3_a, dpgamma_ho2_tot

  USE mo_boundary_condition, ONLY: bc_apply, bc_find
  USE mo_vphysc,             ONLY: vphysc


  ! calculation of khet(s) for salsa microphysics

  !--- arguments
  INTEGER, INTENT(IN)         :: lat,idx, kproma
  REAL(dp),INTENT(IN)         :: mmspd(plonl,plev), dg(plonl,plev), &
                                 prhop1(plonl,plev), temp(plonl,plev)
  REAL(dp),INTENT(OUT)        :: khet(plonl,plev)


  !--- local variables
  INTEGER              :: jclass, jk, jl, ibc, ierr, jt, ilastmod
  REAL(dp)             :: gamma(plonl,plev),            &
                          radius(plonl,plev),           & ! dry or wet radius [m]
                          zmvolr(plonl,plev),           & ! cloud volume mean radius [m]
                          za,                           & ! dry or wet surface area density
                          za_tot(plonl,plev),           & ! dry or wet surface total area denstiy
                          za_cloud_tot(plonl,plev),     & ! cloud droplet surface total area denstiy
                          khet_cloud(plonl,plev),       & ! reaction rate coefficient HO2 on cloud [s-1]
                          khet_aero(plonl,plev),        & ! reaction rate coefficient on aerosol [s-1]
                          gamma_tot(plonl,plev),        & ! Sum over reactionprobabilities for single species
                          numden(plonl,plev),           & ! number density [1/kg]
                          zcdnc(plonl,plev),            & ! cloud droplet number concentration [1/m3]
                          zclw(plonl,plev),             & ! cloud liquid water [kg/kg]
                          zaclc(plonl,plev),            & ! cloud fraction [1]
                          zmass(plonl,plev,NAERO_SPEC), & ! [kg/kg]
                          taudg(plonl,plev),            & ! only variables to calculate khet, but
                          tauim(plonl,plev)               ! named as characteristic times they can be interpreted
                                                          ! taudg = diffusion
                                                          ! life time
                                                          ! tauim = interfaical
                                                          ! mass transport life
                                                          ! time
  LOGICAL              :: lhydrolysis                     ! switch for hydrolysis heterogeneous reactions
                                                          ! only happening on aqueous aerosol

  ! Initialization
  za_tot = 0.0_dp
  za_cloud_tot = 0.0_dp
  gamma_tot = 0.0_dp
  khet = 0._dp
  khet_aero = 0._dp
  khet_cloud = 0._dp
  zcdnc(:,:) = 0.0_dp
  zclw(:,:)  = 0.0_dp
  zaclc(:,:) = 0.0_dp
  lhydrolysis = .FALSE.


  ! Calculate khet
  !                A
  ! khet = --------------------
  !        r/dg + 4/gamma mmspd
  ! Chang 2011


  IF (idx == idx_no3)  lhydrolysis = .TRUE.
  IF (idx == idx_no2)  lhydrolysis = .TRUE.
  IF (idx == idx_ho2)  lhydrolysis = .TRUE.
  IF (idx == idx_hno3) lhydrolysis = .TRUE.

  DO jclass=1, nclass 
    zmass = 0.0_dp
    DO jt=1,NAERO_SPEC
      ibc = ibc_aero_comp_modes(jclass,jt)
      IF (ibc > 0) CALL bc_apply(ibc,kproma,lat,zmass(:,:,jt))
    END DO
    CALL bc_apply(ibc_numden_modes(jclass),kproma,lat,numden(:,:))
!!    IF(jclass .le. nsol) THEN
    IF(lhydrolysis .AND. .NOT. sizeclass(jclass)%lsoluble) CONTINUE ! Falls lhydrolysis=False, macht er alle jclass
      CALL bc_apply(ibc_rwet(jclass),kproma,lat,radius)
!!    ENDIF
    CALL calc_gamma(plonl, plev, kproma, temp, radius, mmspd, idx, lat, zmass, gamma)
    DO jk=1, plev
      DO jl=1,kproma !plonl
        za= 4.e-02_dp*pi*numden(jl,jk)*prhop1(jl,jk)*radius(jl,jk)**2._dp ! [numden] = [#/kg]
        IF (za < 0._dp) za = 0._dp
        za_tot(jl,jk) = za_tot(jl,jk) + za ! just for diagnostics
        IF(gamma(jl,jk) > 1.e-12_dp .AND. radius(jl,jk) > 1.e-9_dp) THEN
          taudg = radius(jl,jk)/dg(jl,jk)
          tauim = 4/(gamma(jl,jk)*mmspd(jl,jk))
          khet_aero(jl,jk) = khet_aero(jl,jk) + za/(taudg(jl,jk) + tauim(jl,jk))
          gamma_tot(jl,jk) = gamma_tot(jl,jk) + gamma(jl,jk)*za ! just for diagnostics
        ENDIF
      ENDDO
    ENDDO
  END DO


  IF (idx == idx_no3 ) khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_no2 ) khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_n2o5) khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_hno3) khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_o3)   khet(1:kproma,:) = khet_aero(1:kproma,:)
  !IF (idx == idx_ho2)  khet(1:kproma,:) = khet_aero(1:kproma,:)
  IF (idx == idx_ho2) THEN
    zcdnc(1:kproma,:) = vphysc%cdnc(1:kproma,:,lat)
    zclw(1:kproma,:)  = vphysc%clw(1:kproma,:,lat)
    zaclc(1:kproma,:) = vphysc%aclc(1:kproma,:,lat)

    ! Calculation mean volume radius as in mo_cloud_micro_2m.f90
    zmvolr(1:kproma,:) = (0.75_dp*MAX(zclw(1:kproma,:),0._dp) * prhop1(1:kproma,:) &      ! zclw = zxlb
                         / (pi*rhoh2o*zcdnc(1:kproma,:)))**(1._dp/3._dp)

    CALL calc_gamma_cloud(plonl, plev, kproma, temp, idx, mmspd, zmvolr, lat, gamma)
    za = 0.0_dp
    DO jk=1, plev
      DO jl=1,kproma !plonl
        za = 4.e-02_dp*pi*zmvolr(jl,jk)**2._dp*zcdnc(jl,jk)*0.01_dp ![cm2/cm3]
        IF (za < 0._dp) za = 0._dp ! set to 0, because zcdnc or zmvolr can be < 0
        IF (zmvolr(jl,jk) < 0._dp) zmvolr(jl,jk) = 0._dp
        za_cloud_tot(jl,jk) = za_cloud_tot(jl,jk) + za ! just for diagnostics
        IF(gamma(jl,jk) > 1.e-12_dp) THEN
          taudg = zmvolr(jl,jk)/dg(jl,jk)
          tauim = 4/(gamma(jl,jk)*mmspd(jl,jk))
          khet_cloud(jl,jk) = za/(taudg(jl,jk) + tauim(jl,jk))
        END IF
      ENDDO
    ENDDO

    DO jk=1, plev
      DO jl=1,kproma !plonl
        khet(jl,jk) = (1.0_dp-zaclc(jl,jk))*khet_aero(jl,jk) + zaclc(jl,jk)*khet_cloud(jl,jk)
      END DO
    END DO
  END IF


  ! diagnostics (see ldetail2 in mo_moz_diag.f90)
  IF (idx == idx_n2o5) THEN
    dpsa_aerosol(1:kproma,:,lat) = za_tot(1:kproma,:)
    dpgamma_n2o5(1:kproma,:,lat) = 0._dp
    dpgamma_n2o5_a(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_tot(jl,jk) > 0._dp) THEN
          dpgamma_n2o5(jl,jk,lat) = gamma_tot(jl,jk)/za_tot(jl,jk) ! gamma weighted with modes and surfaces
          dpgamma_n2o5_a(jl,jk,lat) = 4.0_dp*khet(jl,jk)/(za_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO
  ELSE IF (idx == idx_ho2) THEN
    dpsa_aerosol_wet(1:kproma,:,lat) = za_tot(1:kproma,:)
    dpsa_cloud(1:kproma,:,lat) = za_cloud_tot(1:kproma,:)
    dpsa_cloud(1:kproma,:,lat) = 0.0_dp
    dpgamma_ho2(1:kproma,:,lat) = 0._dp
    dpgamma_ho2_a(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_tot(jl,jk) > 0._dp) THEN
          dpgamma_ho2(jl,jk,lat) = gamma_tot(jl,jk)/za_tot(jl,jk) ! gamma weighted with modes and surfaces
          dpgamma_ho2_a(jl,jk,lat) = 4.0_dp*khet_aero(jl,jk)/(za_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO

    dpgamma_ho2_cloud(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_cloud_tot(jl,jk) > 0._dp) THEN
          dpgamma_ho2_cloud(jl,jk,lat) = 4.0_dp*khet_cloud(jl,jk)/(za_cloud_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO

    dpgamma_ho2_tot(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_cloud_tot(jl,jk) > 0._dp .and. za_tot(jl,jk) > 0._dp) THEN
          dpgamma_ho2_tot(jl,jk,lat) = 4.0_dp*khet(jl,jk)/((za_cloud_tot(jl,jk)+za_tot(jl,jk))*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO

  ELSE IF (idx == idx_hno3) THEN
    dpgamma_hno3(1:kproma,:,lat) = 0._dp
    dpgamma_hno3_a(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_tot(jl,jk) > 0._dp) THEN
          dpgamma_hno3(jl,jk,lat) = gamma_tot(jl,jk)/za_tot(jl,jk) ! gamma weighted with modes and surfaces
          dpgamma_hno3_a(jl,jk,lat) = 4.0_dp*khet(jl,jk)/(za_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO

  ELSE IF (idx == idx_o3) THEN
    dpgamma_o3(1:kproma,:,lat) = 0._dp
    dpgamma_o3_a(1:kproma,:,lat) = 0._dp
    DO jk=1, plev
      DO jl=1,kproma
        IF (za_tot(jl,jk) > 0._dp) THEN
          dpgamma_o3(jl,jk,lat) = gamma_tot(jl,jk)/za_tot(jl,jk) ! gamma weighted with modes and surfaces
          dpgamma_o3_a(jl,jk,lat) = 4.0_dp*khet(jl,jk)/(za_tot(jl,jk)*mmspd(jl,jk)) ! gamma calculated over k
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  END SUBROUTINE calc_khet_salsa

END MODULE mo_hammoz_het_tropo_rates
