!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_m7_trac.f90
!!
!! \brief
!! mo_ham_m7_trac contains routines to requests tracers for ECHAM/HAM and 
!! prescribes their physical and chemical properties.
!! It controls the aerosol physics by providing the necessary switches.
!!
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# P. Stier (MPI-Met) - original code (2001)
!!   -# D. O'Donnell (MPI-Met) - code generalization and changes for soa (2009-02-xx)
!!   -# K. Zhang (MPI-Met) - adaption for new species list and tracer defination (2009-08-11) 
!!   -# M.G. Schultz (FZ Juelich) - cleanup and adaptation to new structure (2009-11-20)
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

!!### Questionable whether this module is needed with new emissions scheme.
!! idt_ seem to be used primarily for identiyfing tracers for emissions.

MODULE mo_ham_m7_trac

  ! Parameters:
  ! -----------
  ! User defined flags: density   density                    [kg m-3]
  !                     osm       osmotic coefficient        [???]
  !                     nion      number of ions the tracer 
  !                               dissolves into             [1]

  USE mo_kind,          ONLY: dp
  USE mo_tracdef,       ONLY: ntrac,                & ! number of tracers
                              OFF, ON,              & ! ON/OFF index
                              GAS,                  & ! phase indicators
                              AEROSOLMASS,          & !
                              AEROSOLNUMBER,        & !
                              SOLUBLE,              & ! soluble indicator
                              INSOLUBLE, &
                              itrprog, itrdiag, itrpresc
  USE mo_species,       ONLY: speclist
  USE mo_physical_constants, ONLY: rhoh2o 
  USE mo_ham_species,   ONLY: id_dms, id_so2, id_so4g, id_oh, id_h2o2, id_o3, &
                              id_no2, id_so4, id_bc, id_oc, id_ss, id_du, id_wat 
  
  IMPLICIT NONE

  !--- Public entities:

  PUBLIC :: idt_dms,   idt_so2,   idt_so4,   idt_ocnv,     &
            idt_ms4ns, idt_ms4ks, idt_ms4as, idt_ms4cs,    &
            idt_mbcki, idt_mbcks, idt_mbcas, idt_mbccs,    &
            idt_mocki, idt_mocks, idt_mocas, idt_moccs,    &
            idt_mssas, idt_msscs,                          &
            idt_mduai, idt_mduas, idt_mduci, idt_mducs,    &
            idt_nns,   idt_nki,   idt_nks,   idt_nai,      &
            idt_nas,   idt_nci,   idt_ncs,                 &
            idt_cdnc_ham,  idt_icnc_ham,                   &
            idt_mwans, idt_mwaks, idt_mwaas, idt_mwacs
 
  PUBLIC:: ham_m7_set_idt
  PUBLIC:: ham_get_class_flag           ! for tracer diagnostics

  !--- Module variables:
  !
  !    Tracer indices:
  !
  !    Legend: iABBCD
  !
  !            A:  m  = particle mass mixing ratio, n number mixing ratio
  !            BB: s4 = sulfate, bc/oc = black/organic carbon, du = dust, ss = seasalt
  !            C:  n  = nucleation , k = Aitken, a = accumulation, c = coarse mode
  !            D:  i  = insoluble,  s = soluble

  INTEGER :: idt_dms    ! mass mixing ratio dms
  INTEGER :: idt_so2    ! mass mixing ratio so2
  INTEGER :: idt_so4    ! mass mixing ratio so4
  INTEGER :: idt_ocnv   ! mass mixing ratio nonvolatile organic

  INTEGER :: idt_ms4ns  ! mass mixing ratio sulfate        nuclea. soluble
  INTEGER :: idt_ms4ks  ! mass mixing ratio sulfate        aitken  soluble
  INTEGER :: idt_ms4as  ! mass mixing ratio sulfate        accum.  soluble
  INTEGER :: idt_ms4cs  ! mass mixing ratio sulfate        coarse  soluble
  INTEGER :: idt_mbcki  ! mass mixing ratio black carbon   aitken  insoluble
  INTEGER :: idt_mbcks  ! mass mixing ratio black carbon   aitken  soluble
  INTEGER :: idt_mbcas  ! mass mixing ratio black carbon   accum.  soluble
  INTEGER :: idt_mbccs  ! mass mixing ratio black carbon   coarse  soluble
  INTEGER :: idt_mocki  ! mass mixing ratio organic carbon aitken  insoluble
  INTEGER :: idt_mocks  ! mass mixing ratio organic carbon aitken  soluble
  INTEGER :: idt_mocas  ! mass mixing ratio organic carbon accum.  soluble
  INTEGER :: idt_moccs  ! mass mixing ratio organic carbon coarse  soluble
  INTEGER :: idt_mssas  ! mass mixing ratio seasalt        accum.  soluble
  INTEGER :: idt_msscs  ! mass mixing ratio seasalt        coarse  soluble
  INTEGER :: idt_mduai  ! mass mixing ratio dust           accum.  insoluble
  INTEGER :: idt_mduas  ! mass mixing ratio dust           accum.  soluble
  INTEGER :: idt_mduci  ! mass mixing ratio dust           coarse  insoluble
  INTEGER :: idt_mducs  ! mass mixing ratio dust           coarse  soluble
  INTEGER :: idt_mwans  ! mass mixing ratio aerosol water  nuclea. soluble
  INTEGER :: idt_mwaks  ! mass mixing ratio aerosol water  aitken  soluble
  INTEGER :: idt_mwaas  ! mass mixing ratio aerosol water  accum.  soluble
  INTEGER :: idt_mwacs  ! mass mixing ratio aerosol water  coarse  soluble


  INTEGER :: idt_nns    ! number mixing ratio              nuclea. soluble
  INTEGER :: idt_nki    ! number mixing ratio              aitken  insoluble
  INTEGER :: idt_nks    ! number mixing ratio              aitken  soluble
  INTEGER :: idt_nai    ! number mixing ratio              accum.  insoluble
  INTEGER :: idt_nas    ! number mixing ratio              accum.  soluble
  INTEGER :: idt_nci    ! number mixing ratio              coarse  insoluble
  INTEGER :: idt_ncs    ! number mixing ratio              coarse  soluble

  INTEGER :: idt_cdnc_ham   ! cloud droplet number concentration
  INTEGER :: idt_icnc_ham   ! ice   cristal number concentration

  
  CONTAINS




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! Define HAM tracers 
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! to_be_filled
!!
!! @par Externals:
!! <ol>
!! <li>none
!! </ol>
!!
!! @par Notes
!! 
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 

  SUBROUTINE ham_m7_set_idt

  USE mo_ham,           ONLY: aerocomp, aerowater, sizeclass
  USE mo_ham_m7ctl,     ONLY: inucs,  iaits,  iaccs,  icoas,   &
                              iaiti,  iacci,  icoai,           &
                              iso4ns, iso4ks, iso4as, iso4cs,  &
                              ibcks,  ibcas,  ibccs,  ibcki,   &
                              iocks,  iocas,  ioccs,  iocki,   &
                              issas,  isscs,                   &
                              iduas,  iducs,  iduai,  iduci

 
       idt_dms  = speclist(id_dms)%idt
       idt_so2  = speclist(id_so2)%idt
       idt_so4  = speclist(id_so4g)%idt
  
       idt_ms4ns = aerocomp(iso4ns)%idt
       idt_ms4ks = aerocomp(iso4ks)%idt
       idt_ms4as = aerocomp(iso4as)%idt
       idt_ms4cs = aerocomp(iso4cs)%idt

       idt_mbcks = aerocomp(ibcks)%idt
       idt_mbcas = aerocomp(ibcas)%idt
       idt_mbccs = aerocomp(ibccs)%idt
       idt_mbcki = aerocomp(ibcki)%idt

       idt_mocks = aerocomp(iocks)%idt
       idt_mocas = aerocomp(iocas)%idt
       idt_moccs = aerocomp(ioccs)%idt
       idt_mocki = aerocomp(iocki)%idt

       idt_mssas = aerocomp(issas)%idt
       idt_msscs = aerocomp(isscs)%idt

       idt_mduas = aerocomp(iduas)%idt
       idt_mducs = aerocomp(iducs)%idt
       idt_mduai = aerocomp(iduai)%idt
       idt_mduci = aerocomp(iduci)%idt

       idt_nns = sizeclass(inucs)%idt_no
       idt_nks = sizeclass(iaits)%idt_no
       idt_nas = sizeclass(iaccs)%idt_no
       idt_ncs = sizeclass(icoas)%idt_no
       idt_nki = sizeclass(iaiti)%idt_no
       idt_nai = sizeclass(iacci)%idt_no
       idt_nci = sizeclass(icoai)%idt_no

       idt_mwans = aerowater(inucs)%idt
       idt_mwaks = aerowater(iaits)%idt
       idt_mwaas = aerowater(iaccs)%idt
       idt_mwacs = aerowater(icoas)%idt

END SUBROUTINE ham_m7_set_idt 

!@brief: set a list of flag values and mode names for th eindividual aerosol modes
!
! The flag values depend on the optional property flags ldrydep, ...
!
! @author: Martin Schultz, FZ Juelich (2010-04-16)
!
SUBROUTINE ham_get_class_flag(nclass, classflag, classname, classnumname,        &
                             ldrydep, lwetdep, lsedi)         !>>dod added lsedi <<dod

  USE mo_ham,             ONLY: sizeclass, aero_nclass=>nclass
  USE mo_tracdef,         ONLY: ln
  !>>dod
  USE mo_exception,       ONLY:finish
  !<<dod

  INTEGER,           INTENT(out) :: nclass
  LOGICAL,           INTENT(out) :: classflag(aero_nclass)
  CHARACTER(len=ln), INTENT(out) :: classname(aero_nclass)
  CHARACTER(len=ln), INTENT(out) :: classnumname(aero_nclass) !SF #299
  LOGICAL, OPTIONAL, INTENT(in)  :: ldrydep, lwetdep  ! set flag true only for modes that are deposited
  LOGICAL, OPTIONAL, INTENT(in)  :: lsedi             !>>dod

  INTEGER        :: jclass

  !>>dod
  IF (PRESENT(lwetdep) .AND. PRESENT(ldrydep)) CALL finish('mo_ham_m7_trac', 'ham_get_class_flag received multiple requests')
  IF (PRESENT(lwetdep) .AND. PRESENT(lsedi))   CALL finish('mo_ham_m7_trac', 'ham_get_class_flag received multiple requests')
  IF (PRESENT(ldrydep) .AND. PRESENT(lsedi))   CALL finish('mo_ham_m7_trac', 'ham_get_class_flag received multiple requests')
  !
  ! define number of values returned
  nclass = aero_nclass
  ! define values
  DO jclass = 1,nclass
    classname(jclass) = sizeclass(jclass)%shortname
    !>>SF #299: added mode number name
    classnumname(jclass) = 'NUM_'//TRIM(sizeclass(jclass)%shortname)
    !<<SF
  END DO
  ! define flag values
  classflag(:) = .FALSE.

  !>>dod: drydep is true for all modes
  IF (PRESENT(ldrydep)) THEN     
     IF (ldrydep) classflag(:) = .TRUE.
  END IF
  !<<dod
  
  !>>dod: so is wetdep...
  IF (PRESENT(lwetdep)) THEN
     IF (lwetdep) classflag(:) = .TRUE.
  END IF

  !>>dod
  IF (PRESENT(lsedi)) THEN
    IF (lsedi) THEN
      DO jclass = 1,nclass
        classflag(jclass) = sizeclass(jclass)%lsed
      END DO
    END IF
  END IF
  !<<dod

END SUBROUTINE ham_get_class_flag

END MODULE mo_ham_m7_trac
