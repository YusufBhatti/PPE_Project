!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_moz_lbc_ubc.f90
!!
!! \brief
!! This module defines and applies lower and upper species boundary conditions for MOZART
!!
!! \author Martin Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2011-05-19)
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


MODULE mo_moz_lbc_ubc

  USE mo_kind,               ONLY: dp
  USE mo_moz_mods,           ONLY: pcnstm1                   ! number of MOZ tracers

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: moz_lbc_ubc_initialize, moz_lbc_ubc

  !-- private module variables
  INTEGER          :: idt_lbc(pcnstm1),         & ! tracer indices for boundary conditions
                      idt_ubc(pcnstm1)

  INTEGER          :: nlbc, nubc                  ! number of lower and upper boundary conditions
  INTEGER          :: ibc_lbc, ibc_ubc            ! first BC index 



CONTAINS

! ------------------------------------------------------------------------------
! private routine to define lower or upper BC for one species
!
SUBROUTINE define_species_bc(specname, title, bc_struc, ifirstbc, bccount, idt_moz)

  USE mo_exception,          ONLY: message, em_error, em_info
  USE mo_boundary_condition, ONLY: bc_nml, bc_define
  USE mo_tracdef,            ONLY: trlist
  USE mo_physical_constants, ONLY: amd
  USE mo_moz_util,           ONLY: get_spc_ndx
  USE mo_moz_mods,           ONLY: adv_mass

  CHARACTER(LEN=*), INTENT(in)     :: specname        ! species name
  CHARACTER(LEN=*), INTENT(in)     :: title           ! boundary condition title
  TYPE(bc_nml), INTENT(in)         :: bc_struc        ! namelist structure for BC
  INTEGER, INTENT(inout)           :: ifirstbc        ! index of BC if not already > 0
  INTEGER, INTENT(inout)           :: bccount         ! number of BCs defined
  INTEGER, INTENT(out)             :: idt_moz         ! tracer index (will be <= 0 if tracer not found)

  TYPE(bc_nml)         :: mybc_struc
  INTEGER              :: bcindex

  ! find out if species exists in MOZART scheme
  idt_moz = get_spc_ndx(TRIM(specname))
  IF (idt_moz > 0) THEN
    mybc_struc = bc_struc
    mybc_struc%ef_varname = TRIM(specname)
    mybc_struc%ef_factor=adv_mass(idt_moz)/amd
    bcindex = bc_define(TRIM(specname)//' '//TRIM(title),mybc_struc, 2, .TRUE.)
    IF (ifirstbc <= 0) ifirstbc = bcindex
    bccount = bccount + 1
  ELSE
    CALL message('define_species_bc', TRIM(title)//' requested for '//TRIM(specname)//    &
                 ' but tracer not defined!', level=em_error)
  END IF

END SUBROUTINE define_species_bc
! ------------------------------------------------------------------------------
!
SUBROUTINE moz_lbc_ubc_initialize (lbc_species, ubc_species, ndeflbc, ndefubc, bc_lbc, bc_ubc)

  USE mo_exception,          ONLY: message, message_text, em_warn, em_info
  USE mo_boundary_condition, ONLY: bc_nml, BC_RELAX, BC_REPLACE
  USE mo_util_string,        ONLY: tolower

  CHARACTER(LEN=*), INTENT(in)     :: lbc_species(:),  &     ! namelist input for LBC species
                                      ubc_species(:)         ! dto. for UBC species
  INTEGER, INTENT(in)              :: ndeflbc, ndefubc       ! number of namelist entries
  TYPE(bc_nml), INTENT(in)         :: bc_lbc, bc_ubc         ! namelist definitions (type, geometry,  etc.)

  INTEGER          :: ii, idt
  TYPE(bc_nml)     :: bc_ch4

  !-- initialize
  idt_lbc(:) = 0
  idt_ubc(:) = 0

  !-----------------------------------------------
  !-- define lower boundary conditions
  nlbc = 0
  ibc_lbc = 0
has_lbc:     &
  IF (ndeflbc > 0) THEN
    DO ii = 1, ndeflbc
      IF (ii == 1 .AND. tolower(TRIM(lbc_species(ii))) == 'default') THEN
        ! evaluate default list
        CALL define_species_bc('CO2', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        ! original MOZ3.5 has CH4 and H2 - ECHAM-MOZ has emissions for these
!+++sschr
        bc_ch4 = bc_lbc
        bc_ch4%bc_mode = BC_RELAX
        bc_ch4%bc_relaxtime = 3._dp * 86400._dp 
        CALL define_species_bc('CH4', 'lower boundary', bc_ch4, ibc_lbc, nlbc, idt)
!---sschr
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('N2O', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CFC11', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CFC12', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CFC113', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CFC114', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CFC115', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CCL4', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CH3CCL3', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('HCFC22', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('HCFC141B', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('HCFC142B', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CH3CL', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CH3BR', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CH2BR2', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CHBR3', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('H1202', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('H2402', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CF2CLBR', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('CF3BR', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
        CALL define_species_bc('SF6', 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
      ELSE
        ! evaluate species names from namelist
        CALL define_species_bc(lbc_species(ii), 'lower boundary', bc_lbc, ibc_lbc, nlbc, idt)
        IF (idt > 0) idt_lbc(nlbc) = idt
      END IF
    END DO
    WRITE (message_text, '(i0,a)') nlbc, ' lower boundary conditions defined.'
    CALL message('moz_lbc_initialize', message_text, level=em_info)
  ELSE
    CALL message('moz_lbc_initialize', 'No lower species boundary conditions defined!', &
                 level=em_warn)
  END IF has_lbc


  !-----------------------------------------------
  !-- define upper boundary conditions
  nubc = 0
  ibc_ubc = 0
has_ubc:     &
  IF (ndefubc > 0) THEN
    DO ii = 1, ndefubc
      ! (no defaults for UBC)
      CALL define_species_bc(ubc_species(ii), 'upper boundary', bc_ubc, ibc_ubc, nubc, idt)
      IF (idt > 0) idt_ubc(nubc) = idt
    END DO
    WRITE (message_text,'(i0,a)') nubc, ' upper boundary conditions defined.'
    CALL message('moz_ubc_initialize', message_text, level=em_info)
  ELSE
    CALL message('moz_ubc_initialize', 'No upper species boundary conditions defined!',   &
                 level=em_info)
  END IF has_ubc

END SUBROUTINE moz_lbc_ubc_initialize


! ------------------------------------------------------------------------------
!

SUBROUTINE moz_lbc_ubc (kproma, krow, pmmr)
!!! kproma <-> plonl ?!
  USE mo_boundary_condition, ONLY: bc_apply
  USE mo_moz_mods,           ONLY: plonl, plev

  INTEGER,  INTENT(in)    :: kproma, krow
  REAL(dp), INTENT(inout) :: pmmr(plonl,plev,pcnstm1)

  integer :: i

  DO i = 1, nlbc
    CALL bc_apply(ibc_lbc+i-1, kproma, krow, pmmr(:,:,idt_lbc(i)))
  END DO
END SUBROUTINE moz_lbc_ubc

END MODULE mo_moz_lbc_ubc
