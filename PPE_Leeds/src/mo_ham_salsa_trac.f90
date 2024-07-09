!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_salsa_trac.f90
!!
!! \brief
!! mo_ham_salsa_trac contains routines to requests tracers for ECHAM/HAM and 
!! prescribes their physical and chemical properties.
!! It controls the aerosol physics by providing the necessary switches.
!!
!! \author Tommi Bergman (FMI)
!!
!! \responsible_coder
!! Tommi Bergman tommi.bergman@iki.fi
!!
!! \revision_history
!!   -# P. Stier (MPI-Met) - original code (2001)
!!   -# D. O'Donnell (MPI-Met) - code generalization and changes for soa (2009-02-xx)
!!   -# K. Zhang (MPI-Met) - adaption for new species list and tracer defination (2009-08-11) 
!!   -# M.G. Schultz (FZ Juelich) - cleanup and adaptation to new structure (2009-11-20)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!   -# A. Laakso (FMI) - Tracers for SALSA
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


MODULE mo_ham_salsa_trac

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
  USE mo_ham_species,   ONLY: id_dms, id_so2, id_so4g, id_ocnv, id_oh, id_h2o2, id_o3, &
                              id_no2, id_so4, id_bc, id_oc, id_ss, id_du, id_wat
  USE mo_ham_salsactl,  ONLY: in1a, fn1a,           &
                              in2a, fn2a,           &
                              in2b, fn2b, nbins
  
  IMPLICIT NONE

  !--- Public entities:

  PUBLIC :: idt_dms,   idt_so2,   idt_so4,   idt_ocnv,     &
            idt_cdnc_ham,  idt_icnc_ham,                   &
            !idt_mwans, idt_mwaks, idt_mwaas, idt_mwacs     &
        idt_ms4,   idt_moc,   idt_mbc,   idt_mss,      & ! SALSA indices
            idt_mdu,   idt_mws,   idt_n, idt_mwa
 
  PUBLIC:: ham_salsa_set_idt
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

  INTEGER :: idt_cdnc_ham   ! cloud droplet number concentration
  INTEGER :: idt_icnc_ham   ! ice   cristal number concentration

  INTEGER :: idt_ms4(fn2b)    ! mass mixing ratio sulfate
  INTEGER :: idt_moc(fn2b)    ! mass mixing ratio organic carbon
  INTEGER :: idt_mbc(fn2b)    ! mass mixing ratio black carbon
  INTEGER :: idt_mss(fn2a)    ! mass mixing ratio seasalt
  INTEGER :: idt_mdu(fn2b)    ! mass mixing ratio mineral dust
  INTEGER :: idt_mws(fn2b)    ! mass mixing ratio water soluble
  INTEGER :: idt_n(fn2b)      ! number mixing ratio
  INTEGER :: idt_mwa(nbins)   ! mass mixing ratio aerosol water

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
  !! ham_initialize
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

  SUBROUTINE ham_salsa_set_idt

  USE mo_ham,           ONLY: aerocomp, aerowater, sizeclass
  USE mo_ham_salsactl,  ONLY:  iso4b, ibcb,iocb,issb, idub,&
                               in1a, in2a, in2b,  &
                               fn1a, fn2a, fn2b

   INTEGER :: i

   idt_dms  = speclist(id_dms)%idt
   idt_so2  = speclist(id_so2)%idt
   idt_so4  = speclist(id_so4g)%idt
   idt_ocnv = speclist(id_ocnv)%idt
   idt_ms4(:)=0
   idt_moc(:)=0
   idt_n(:)=0
   idt_mss(:)=0
   idt_mbc(:)=0
   idt_mdu(:)=0
   DO i = in1a,fn2b
       idt_ms4(i) = aerocomp(iso4b(i))%idt
       idt_moc(i) = aerocomp(iocb(i))%idt
       idt_n(i) = sizeclass(i)%idt_no
   END DO

   DO i = in2a,fn2b
       idt_mbc(i) = aerocomp(ibcb(i))%idt
       idt_mdu(i) = aerocomp(idub(i))%idt
   END DO

   DO i = in2a,fn2a
        idt_mss(i) = aerocomp(issb(i))%idt
   END DO

   DO i = in1a,fn2a
        idt_mwa(i) = aerowater(i)%idt
   END DO

END SUBROUTINE ham_salsa_set_idt 

!@brief: set a list of flag values and mode names for th eindividual aerosol modes
!
! The flag values depend on the optional property flags ldrydep, ...
!
! @author: Martin Schultz, FZ Juelich (2010-04-16)
! 
!
SUBROUTINE ham_get_class_flag(nclass, modflag, modname,         &
                             ldrydep, lwetdep, lsedi)         !>>dod added lsedi <<dod

  USE mo_ham,             ONLY: sizeclass, aero_nmod=>nclass
  USE mo_tracdef,         ONLY: ln
  !>>dod
  USE mo_exception,       ONLY:finish
  !<<dod

  INTEGER,           INTENT(out) :: nclass
  LOGICAL,           INTENT(out) :: modflag(aero_nmod)
  CHARACTER(len=ln), INTENT(out) :: modname(aero_nmod)
  LOGICAL, OPTIONAL, INTENT(in)  :: ldrydep, lwetdep  ! set flag true only for modes that are deposited
  LOGICAL, OPTIONAL, INTENT(in)  :: lsedi             !>>dod

  INTEGER        :: jclass

  !>>dod
  IF (PRESENT(lwetdep) .AND. PRESENT(ldrydep)) CALL finish('mo_ham_salsa_trac', 'ham_get_bin_flag received multiple requests')
  IF (PRESENT(lwetdep) .AND. PRESENT(lsedi))   CALL finish('mo_ham_salsa_trac', 'ham_get_bin_flag received multiple requests')
  IF (PRESENT(ldrydep) .AND. PRESENT(lsedi))   CALL finish('mo_ham_salsa_trac', 'ham_get_bin_flag received multiple requests')
  !
  ! define number of values returned
  nclass = aero_nmod
  ! define values
  DO jclass = 1,nclass
    modname(jclass) = sizeclass(jclass)%shortname
  END DO
  ! define flag values
  modflag(:) = .FALSE.

  !>>dod: drydep is true for all modes
  IF (PRESENT(ldrydep)) THEN     
     IF (ldrydep) modflag(:) = .TRUE.
  END IF
  !<<dod
  
  !>>dod: so is wetdep...
  IF (PRESENT(lwetdep)) THEN
     IF (lwetdep) modflag(:) = .TRUE.
  END IF

  !>>dod
  IF (PRESENT(lsedi)) THEN
    IF (lsedi) THEN
      DO jclass = 1,nclass
        modflag(jclass) = sizeclass(jclass)%lsed
      END DO
    END IF
  END IF
  !<<dod

END SUBROUTINE ham_get_class_flag

END MODULE mo_ham_salsa_trac
