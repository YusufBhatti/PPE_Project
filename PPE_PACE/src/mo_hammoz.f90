!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz.f90
!!
!! \brief
!! contains routines which couple the HAM and MOZ models
!!
!! \author M. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2010-04-20)
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
MODULE mo_hammoz

  IMPLICIT NONE

  PRIVATE

  SAVE

  PUBLIC  :: hammoz_initialize 
  PUBLIC  :: hammoz_set_oxidants
  PUBLIC  :: hammoz_set_aerosol_bc

! Variable declarations
  INTEGER :: idt_oh, idt_o3, idt_h2o2, idt_no2, idt_no3


  CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! hammoz_initialize: perform cross-initialidsation of HAM and MOZ submodels
!! This includes initialisation of aircraft module
!!
!! @author see module info
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Notes
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE hammoz_initialize

  USE mo_exception,          ONLY: message, em_warn
  USE mo_submodel,           ONLY: lhmzoxi
  USE mo_ham,                ONLY: ibc_oh, ibc_o3, ibc_h2o2, ibc_no2, ibc_no3
  USE mo_tracer,             ONLY: get_tracer
  USE mo_external_field_processor, ONLY: EF_MODULE
  USE mo_boundary_condition, ONLY: bc_modify
  USE mo_util_string,        ONLY: toupper

  CHARACTER(LEN=16)  :: cmodule


  cmodule = ''

  !-- check boundary conditions for HAM oxidant fields
  !   and initialize if needed
  !   HAM defines these BCs as EF_FILE. HAMMOZ coupling requires them as EF_MODULE
  IF (lhmzoxi) THEN
    !-- OH
    CALL get_tracer('OH', idx=idt_oh, modulename=cmodule)
    IF (TRIM(toupper(cmodule)) == 'MOZ') THEN
      CALL bc_modify(ibc_oh, ef_type=EF_MODULE)
    ELSE
      CALL message('hammoz_initialize', 'Cannot redefine boundary condition for HAM as EF_MODULE!'// &
                   ' No MOZ tracer OH found!', level=em_warn)
    END IF
    !-- O3
    CALL get_tracer('O3', idx=idt_o3, modulename=cmodule)
    IF (TRIM(toupper(cmodule)) == 'MOZ') THEN
      CALL bc_modify(ibc_o3, ef_type=EF_MODULE)
    ELSE
      CALL message('hammoz_initialize', 'Cannot redefine boundary condition for HAM as EF_MODULE!'// &
                   ' No MOZ tracer O3 found!', level=em_warn)
    END IF
    !-- H2O2
    CALL get_tracer('H2O2', idx=idt_h2o2, modulename=cmodule)
    IF (TRIM(toupper(cmodule)) == 'MOZ') THEN
      CALL bc_modify(ibc_h2o2, ef_type=EF_MODULE)
    ELSE
      CALL message('hammoz_initialize', 'Cannot redefine boundary condition for HAM as EF_MODULE!'// &
                   ' No MOZ tracer H2O2 found!', level=em_warn)
    END IF
    !-- NO2
    CALL get_tracer('NO2', idx=idt_no2,  modulename=cmodule)
    IF (TRIM(toupper(cmodule)) == 'MOZ') THEN
      CALL bc_modify(ibc_no2, ef_type=EF_MODULE)
    ELSE
      CALL message('hammoz_initialize', 'Cannot redefine boundary condition for HAM as EF_MODULE!'// &
                   ' No MOZ tracer NO2 found!', level=em_warn)
    END IF
    !-- NO3 (SOA)
    idt_no3 = -1
    IF (ibc_no3 > 0) THEN
      CALL get_tracer('NO3', idx=idt_no3, modulename=cmodule)
      IF (TRIM(toupper(cmodule)) == 'MOZ') THEN
        CALL bc_modify(ibc_no3, ef_type=EF_MODULE)
      ELSE
        CALL message('hammoz_initialize', 'Cannot redefine boundary condition for HAM as EF_MODULE!'// &
                     ' No MOZ tracer NO3 found!', level=em_warn)
      END IF
    END IF
  ELSE       ! lhmzoxi
    idt_oh   = -1
    idt_o3   = -1
    idt_h2o2 = -1
    idt_no2  = -1
    idt_no3  = -1
  END IF     ! lhmzoxi


END SUBROUTINE hammoz_initialize

!!@brief: set MOZ aerosols as HAM boundary condition values
!!bc_set does nothing if ibc_x%bc_ef%ef_type != EF_MODULE
SUBROUTINE hammoz_set_aerosol_bc(kproma, kbdim, klev, krow, ktrac, pxt)
  USE mo_kind,                     ONLY: dp
  USE mo_ham,                      ONLY: sizeclass, aerocomp, naerocomp, nclass, nsol
  USE mo_ham_streams,              ONLY: rdry, rwet
  USE mo_moz,                      ONLY: ltrophet
  USE mo_tracdef,                  ONLY: trlist
  USE mo_boundary_condition,       ONLY: bc_set
  USE mo_hammoz_het_tropo_rates,   ONLY: ibc_rwet, ibc_rdry, &
                                         ibc_numden_modes, &
                                         ibc_aero_comp_modes, &
                                         NAERO_SPEC, aero_spec


  !--- arguments
  INTEGER, INTENT(in)    :: kbdim, kproma, klev, krow, ktrac
  REAL(dp), INTENT(in)   :: pxt(kbdim, klev, ktrac)

  !--- local variables
  REAL(dp), POINTER :: rdry_p(:,:,:)
  REAL(dp), POINTER :: rwet_p(:,:,:)
  INTEGER           :: jclass, jspec, jtrac, ihelp


  !--- boundary conditions ibc_rwet, ibc_rdry, ... are only defined 
  !--- if ltrophet = .true.
  IF (ltrophet) THEN
      DO jclass = 1, nclass
        rwet_p     => rwet(jclass)%ptr
        CALL bc_set(ibc_rwet(jclass),kproma,krow,rwet_p(1:kproma,:,krow))
        CALL bc_set(ibc_numden_modes(jclass),kproma,krow,pxt(1:kproma,:,sizeclass(jclass)%idt_no))
      END DO

      ! preliminary code (depends on variables/constants of HAM -- ltrophet is now only able to run
      ! if lham = .true.
      ! after potential changes considering constants of HAM(M7/SALSA) the code of
      ! mo_hammoz_het_tropo_rates.f90 will change
      ! for now:
      ! attention: NAERO_SPEC, aero_spec (mo_hammoz_het_tropo_rates.f90) <-> naerocomp, aerocomp (mo_ham.f90)!!!
      !            (This is not equal -- especially ordering ...)
      ! ### string comparison should be removed in final version

      DO jtrac = 1, naerocomp
        ihelp = 0
        DO jspec = 1, NAERO_SPEC
          IF (TRIM(aerocomp(jtrac)%species%shortname) == TRIM(aero_spec(jspec))) ihelp = jspec
        END DO
        IF (ihelp > 0) THEN

          ! bc set will be done for class,spec,aerocomp:", &
          ! sizeclass(aerocomp(jtrac)%iclass)%shortname,aero_spec(ihelp),trlist%ti(aerocomp(jtrac)%idt)%fullname

          jclass = aerocomp(jtrac)%iclass
          CALL bc_set(ibc_aero_comp_modes(jclass,ihelp),kproma,krow,pxt(1:kproma,:,aerocomp(jtrac)%idt)) 
        END IF
      END DO

    !!Scarlet END IF
  END IF

END SUBROUTINE hammoz_set_aerosol_bc

!!@brief: set MOZ oxidant concentrations as HAM boundary condition values
!!bc_set does nothing if ibc_x%bc_ef%ef_type != EF_MODULE
SUBROUTINE hammoz_set_oxidants(kproma, kbdim, klev, krow, ktrac, pxtm1, pxtte)

  USE mo_kind,               ONLY: dp
  USE mo_time_control,       ONLY: time_step_len
  USE mo_boundary_condition, ONLY: bc_set
  USE mo_ham,                ONLY: ibc_oh, ibc_o3, ibc_h2o2, ibc_no2, ibc_no3

  INTEGER, INTENT(in)    :: kbdim, kproma, klev, krow, ktrac
  REAL(dp), INTENT(in)   :: pxtm1(kbdim, klev, ktrac),     &
                            pxtte(kbdim, klev, ktrac)

  REAL(dp)           :: zxtp1(kbdim, klev)

  IF (idt_oh > 0) THEN
    zxtp1(1:kproma,:) = pxtm1(1:kproma,:,idt_oh)+pxtte(1:kproma,:,idt_oh)*time_step_len
    CALL bc_set(ibc_oh, kproma, krow, zxtp1(1:kproma,:))
  END IF
  IF (idt_o3 > 0) THEN
    zxtp1(1:kproma,:) = pxtm1(1:kproma,:,idt_o3)+pxtte(1:kproma,:,idt_o3)*time_step_len
    CALL bc_set(ibc_o3, kproma, krow, zxtp1(1:kproma,:))
  END IF
  IF (idt_h2o2 > 0) THEN
    zxtp1(1:kproma,:) = pxtm1(1:kproma,:,idt_h2o2)+pxtte(1:kproma,:,idt_h2o2)*time_step_len
    CALL bc_set(ibc_h2o2, kproma, krow, zxtp1(1:kproma,:))
  END IF
  IF (idt_no2 > 0) THEN
    zxtp1(1:kproma,:) = pxtm1(1:kproma,:,idt_no2)+pxtte(1:kproma,:,idt_no2)*time_step_len
    CALL bc_set(ibc_no2, kproma, krow, zxtp1(1:kproma,:))
  END IF
  IF (idt_no3 > 0) THEN
    zxtp1(1:kproma,:) = pxtm1(1:kproma,:,idt_no3)+pxtte(1:kproma,:,idt_no3)*time_step_len
    CALL bc_set(ibc_no3, kproma, krow, zxtp1(1:kproma,:))
  END IF

END SUBROUTINE hammoz_set_oxidants

END MODULE mo_hammoz


