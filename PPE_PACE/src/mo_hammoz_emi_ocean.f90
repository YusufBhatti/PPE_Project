!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_emi_ocean.f90
!!
!! \brief
!! Module for interactive ocean emisisons (DMS, NH3)
!!
!! \author M. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2010-03-01)
!!
!! \limitations
!! None
!!
!! \details
!! Based on code by Silvia Kloster and Philip Stier
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
MODULE mo_hammoz_emi_ocean

  USE mo_kind,             ONLY: dp

  IMPLICIT NONE

  PRIVATE

  ! public variables  (see declaration below)

  ! subprograms
  PUBLIC            :: init_emi_ocean
  PUBLIC            :: init_emi_ocean_stream
  PUBLIC            :: oceani_emissions

  ! boundary condition index for DMS ocean concentrations
  INTEGER           :: ibc_dms_sea

  ! boundary condition indices for emission mass fluxes
  INTEGER           :: ibc_emi_dms, ibc_emi_nh3   ! ibc_emi_co2

  ! diagnostic pointers
  REAL(dp), POINTER :: vp_dms(:,:), vp_nh3(:,:), vp_co2(:,:)


  CONTAINS

  !! ---------------------------------------------------------------------------
  !! initialize interactive ocean emissions for DMS, NH3 (and CO2 ?)
  !! init_emi_ocean is only called when an OCEANI sector is defined in the emi_matrix
  !! Note: should actually check if ef_type is defined as EF_MODULE to allow for file alternative!

  SUBROUTINE init_emi_ocean

  USE mo_exception,            ONLY: message, message_text, em_info
  USE mo_boundary_condition,   ONLY: bc_find, bc_define, bc_nml, BC_REPLACE, BC_BOTTOM
  USE mo_external_field_processor, ONLY: EF_FILE, EF_LONLAT, EF_IGNOREYEAR

  IMPLICIT NONE

  CHARACTER(LEN=64)            :: bc_name
  TYPE(bc_nml)                 :: bc_struc
  INTEGER                      :: ierr

  ! -- Initialisation
  ibc_emi_dms = -1
  ibc_emi_nh3 = -1
  ibc_dms_sea = -1

  !-- locate boundary condition for emission mass flux of DMS and NH3
  !   bc_define was called from init_emissions in mo_emi_interface
  bc_name = 'OCEANI emissions of DMS'
  CALL bc_find(bc_name, ibc_emi_dms, ierr)

  bc_name = 'OCEANI emissions of NH3'
  CALL bc_find(bc_name, ibc_emi_nh3, ierr)

!!  bc_name = 'OCEANI emissions of CO2'
!!  CALL bc_find(bc_name, ibc_emi_co2)

  !-- define boundary condition for ocean DMS concentrations
  IF (ibc_emi_dms > 0) THEN
    bc_struc%ef_type      = EF_FILE
    bc_struc%ef_template  = 'conc_aerocom_DMS_sea.nc'
    bc_struc%ef_varname   = 'DMS_sea'
    bc_struc%ef_geometry  = EF_LONLAT
    bc_struc%ef_timedef   = EF_IGNOREYEAR
    bc_struc%bc_mode      = BC_REPLACE
    bc_struc%bc_domain    = BC_BOTTOM
    ibc_dms_sea = bc_define('DMS seawater concentrations', bc_struc, 2, .TRUE.)
  END IF

  END SUBROUTINE init_emi_ocean


  !! ---------------------------------------------------------------------------
  !! subroutine to define additional diagnostics for the emi diagnostic stream
  SUBROUTINE init_emi_ocean_stream(nsectors, emi_stream, ldiagdetail)

  ! add field to emis stream (init_emi_stream in mo_emi_interface must have been called!)
  ! in principle one could use ldiagdetail to turn off detailed diagnostics for this sector ...

  USE mo_linked_list,   ONLY: SURFACE
  USE mo_memory_base,   ONLY: t_stream, add_stream_element, default_stream_setting, &
                              AUTO, SURFACE

  INTEGER, INTENT(in)      :: nsectors
  TYPE(t_stream), POINTER  :: emi_stream
  LOGICAL, INTENT(inout)   :: ldiagdetail(nsectors)

  ! -- Diagnostics for DMS (and NH3) emissions (piston velocity)
  IF (ibc_emi_dms > 0 .OR. ibc_emi_nh3 > 0) THEN
    CALL default_stream_setting (emi_stream, lrerun    = .FALSE.,     &
                                 laccu     = .FALSE.,     &
                                 lpost     = .TRUE. ,     &
                                 leveltype = SURFACE,     &
                                 table     = 199,         &
                                 code      = AUTO         )
  END IF

  IF (ibc_emi_dms > 0) THEN
    CALL add_stream_element (emi_stream, 'vp_DMS', vp_dms,     &
                             longname='piston velocity for DMS emissions',   &
                             units='m s-1'  )
  END IF
  IF (ibc_emi_nh3 > 0) THEN
    CALL add_stream_element (emi_stream, 'vp_NH3', vp_nh3,     &
                             longname='piston velocity for NH3 emissions',   &
                             units='m s-1'  )
  END IF
!!  IF (ibc_emi_co2 > 0) THEN
!!  CALL add_stream_element (emi_stream, 'vp_co2', vp_co2,     &
!!                           longname='piston velocity for CO2 emissions',   &
!!                           units='m s-1'  )
!!  END IF

  END SUBROUTINE init_emi_ocean_stream

  !! ---------------------------------------------------------------------------
  !! subroutine to handle interactive ocean emissions
  !! call piston velocity and call dms_emisisons routine and store result in boundary condition
  SUBROUTINE oceani_emissions(kproma, kbdim, krow, npist)

  USE mo_boundary_condition,   ONLY: bc_set
  ! arguments

  INTEGER,  INTENT(in)    :: kproma                   ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                    ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: krow                     ! geographic block number
  INTEGER,  INTENT(in)    :: npist                    ! choice of method for computing piston velocity

  ! local variables
  REAL(dp)          :: zmassf1(kbdim)

  ! compute piston velocity
  CALL piston_velocity(kproma, kbdim, krow, npist)
 
  ! calculate DMS emissions
  IF (ibc_emi_dms > 0) THEN
    CALL dms_emissions(kproma, kbdim, krow, zmassf1) 
    CALL bc_set(ibc_emi_dms, kproma, krow, zmassf1(1:kproma))
  END IF

  END SUBROUTINE oceani_emissions

  !! ---------------------------------------------------------------------------
  ! DMS (NH3 and CO2) emissions: first compute piston velocity, then emission mass flux

  SUBROUTINE piston_velocity(kproma, kbdim, krow, npist)

  ! *piston_velocity* calculation of the piston velocity
  !                   for the tracers DMS, CO2 and NH3
  !
  ! Authors:
  ! --------
  ! S. Kloster MPI-MET 28/10/02
  !

  USE mo_memory_g3b,         ONLY: slm
  USE mo_vphysc,             ONLY: vphysc
  USE mo_physical_constants, ONLY: tmelt

  IMPLICIT NONE

  INTEGER,  INTENT(in) :: kproma               ! geographic block number of locations
  INTEGER,  INTENT(in) :: kbdim                ! geographic block maximum number of locations
  INTEGER,  INTENT(in) :: krow                 ! geographic block number
  INTEGER,  INTENT(in) :: npist                ! mode for calculating the piston velocity

  !--- Local Variables:
  INTEGER                :: jl
  REAL(dp)               :: zsst, zzspeed
  REAL(dp)               :: zschmidt_dms, zschmidt_co2, zschmidt_nh3
  REAL(dp)               :: zkw
  REAL(dp)               :: zvp_dms(kbdim), zvp_nh3(kbdim)
  LOGICAL                :: loocean(kbdim)

  loocean(1:kproma)  = ( slm(1:kproma,krow) < 1.e-2_dp )
  zvp_dms(:) = 0._dp
  zvp_nh3(:) = 0._dp

  DO jl=1,kproma

    !--- 1.) Schmidt number

    zsst=vphysc%tsw(jl,krow)-tmelt

    !--- 1.1.)  DMS Schmidt number after Andreae

    !--- Limit SST to valid range as the polynomial gets negative
    !    for T > 35 C:

    zsst=MIN( 35.0_dp , zsst )

    zschmidt_dms=3652.047271_dp-246.99_dp*zsst+8.536397_dp*zsst*zsst     &
                 -0.124397_dp*zsst*zsst*zsst

    !--- 1.2.) CO2 Schmidt number after Jaehne et al. (1987b)
!      zschmidt_co2=2073.1     -125.62*zsst+3.6276  *zsst*zsst     &
!                   -0.043219*zsst*zsst*zsst

    !--- 1.2.) NH3 Schmidt number
    zschmidt_nh3=2073.1     -125.62*zsst+3.6276  *zsst*zsst     &
                 -0.043219*zsst*zsst*zsst

    !--- 2.)  piston velocity in [cm/hr]

    IF (loocean(jl)) THEN   ! ocean
      zzspeed=vphysc%velo10m(jl,krow)

      IF (npist == 1) THEN ! Calculate piston velocity (Liss&Merlivat, 1986)
        IF (zzspeed.GT.3.6_dp.AND.zzspeed.LE.13._dp) THEN
          zkw=2.85_dp*zzspeed-9.65_dp
          zvp_dms(jl)=zkw*(zschmidt_dms/600._dp)**(-0.5_dp)
!         zvp_co2(jl)=zkw*(zschmidt_co2/600.)**(-0.5)
          zvp_nh3(jl)=zkw*(zschmidt_nh3/600.)**(-0.5)
        ELSE IF(zzspeed.LE.3.6_dp) THEN
          zkw=0.17_dp*zzspeed
          zvp_dms(jl)=zkw*(zschmidt_dms/600._dp)**(-2._dp/3._dp)
!         zvp_co2(jl)=zkw*(zschmidt_co2/600.)**(-2./3.)
          zvp_nh3(jl)=zkw*(zschmidt_nh3/600.)**(-2./3.)
        ELSE
          zkw=5.9_dp*zzspeed-49.3_dp
          zvp_dms(jl)=zkw*(zschmidt_dms/600._dp)**(-0.5_dp)
!         zvp_co2(jl)=zkw*(zschmidt_co2/600.)**(-0.5)
          zvp_nh3(jl)=zkw*(zschmidt_nh3/600.)**(-0.5)
        END IF

      ELSE IF (npist == 2) THEN ! Calculate piston velocity (Wanninkhof, 1992)
        zkw=0.31_dp*zzspeed*zzspeed
        zvp_dms(jl)=zkw*(zschmidt_dms/660._dp)**(-0.5_dp)
!       zvp_co2(jl)=zkw*(zschmidt_co2/660.)**(-0.5)
        zvp_nh3(jl)=zkw*(zschmidt_nh3/660.)**(-0.5)

      ELSE IF (npist == 3) THEN  ! Calculate piston velocity (Nightingale, 2000)
        zkw=(0.222_dp*zzspeed*zzspeed+0.333_dp*zzspeed)
        zvp_dms(jl)=zkw*(zschmidt_dms/660._dp)**(-0.5_dp)
!       zvp_co2(jl)=zkw*(zschmidt_co2/660.)**(-0.5)
        zvp_nh3(jl)=zkw*(zschmidt_nh3/660.)**(-0.5)
      END IF
    ELSE                    ! land
      zkw=0._dp
    END IF
  END DO

  !-- store piston velocity as diagnostics
  IF (ibc_emi_dms > 0)  vp_dms(1:kproma, krow) = zvp_dms(1:kproma)
  IF (ibc_emi_nh3 > 0)  vp_nh3(1:kproma, krow) = zvp_nh3(1:kproma)

  END SUBROUTINE piston_velocity


  SUBROUTINE dms_emissions(kproma, kbdim, krow, pmassf)

  ! *dms_emissions* calculated DMS emission fluxes from the
  !                 precalculated piston velocities
  !
  !    Authors:
  !    --------
  !    Johann Feichter, MPI-Met
  !    Philip Stier,    MPI-Met        06/2002
  !    Silvia Kloster,  MPI-Met        10/2002

  USE mo_memory_g3b,         ONLY: slm, seaice
  USE mo_boundary_condition, ONLY: bc_apply
  USE mo_ham,                ONLY: mw_dms

  IMPLICIT NONE

  INTEGER, INTENT(in)    :: kproma                    ! geographic block number of locations
  INTEGER, INTENT(in)    :: kbdim                     ! geographic block maximum number of locations
  INTEGER, INTENT(in)    :: krow                      ! geographic block number
  REAL(dp), INTENT(out)  :: pmassf(kbdim)             ! emission mass flux

  !--- Local Variables:

  REAL(dp)               :: zdmscon(kbdim)
  REAL(dp)               :: zseaice(kbdim)
  LOGICAL                :: loocean(kbdim)            ! ocean mask

  !--- Initialisation
  loocean(1:kproma)  = ( slm(1:kproma,krow) < 1.e-2_dp )
  zseaice(1:kproma)  = seaice(1:kproma, krow)
  pmassf(:) = 0._dp

  !--- Get seawater DMS concentration from boundary condition
  CALL bc_apply(ibc_dms_sea, kproma, krow, zdmscon)

  !--- Mask out land points !SF #458 (replacement of WHERE statements)
  zdmscon(1:kproma) = MERGE(zdmscon(1:kproma), 0._dp, loocean(1:kproma))

  !--- Calculate DMS emissions
  !    conversion of piston velocity from cm/hr to m/s and nmol/L to mmr
  pmassf(1:kproma) = (1._dp-zseaice(1:kproma))                    &
                     * zdmscon(1:kproma) * vp_dms(1:kproma,krow)  &
                     * mw_dms*1.e-11_dp/3600._dp
!!mgs(S)!!                     * 32.064e-11_dp/3600._dp

  END SUBROUTINE dms_emissions

END MODULE mo_hammoz_emi_ocean
