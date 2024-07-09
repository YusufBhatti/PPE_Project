!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mo_species assigns identities and properties to species.
!! Former ECHAM versions had physical and chemical information stored with tracers which
!! meant duplication of fields when more than one tracer was defined for a given species.
!! Examples are aerosol modes or tagged gas-phase tracers.
!! With mo_species, these species characteristics are defined once and the tracerinfo
!! is trimmed down to the essential quantities that vary with tracer. 
!! Some duplication is unavoidable in order to keep flexibility.
!! 
!! @author 
!! <ol> 
!! <li>Declan O'Donnell (MPI-Met)
!! <li>M. Schultz (FZ-Juelich) 
!! </ol>
!!
!! $Id: 1423$
!!
!! @par Revision History
!! <ol> 
!! <li>Declan O'Donnell (MPI-Met) - original version - (2008-xx-xx) 
!! <li>M. Schultz (FZ-Juelich) - generalisation and only one species list for everything  - (2009-06-xx) 
!! <li>K. Zhang (MPI-Met) - merge Declan's and Martin's version - (2009-08-11) 
!! <li>M. Schultz (FZ-Juelich) - further cleanup (2009-09-21)
!! </ol>
!!
!! @par This module is used by
!! to_be_added
!!  
!! @par Notes
!! 
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! ### nphase is additive and lumps all phases of a chemical species
!! ### into one definition of a species. 
!! ### nphase:  1 = gas, 2 = aerosol, (4 = liquid)   (3 = gas + aerosol, etc.)
!! ### must replace IF (speclist()%nphase == GAS) by IF (IAND(speclist()%nphase,GAS))
!! ### must think about diagnostics: BYSPECIES = BYSPECIES_AND_PHASE ??


MODULE mo_species

  !--- inherited types, data and functions
  
  USE mo_kind,         ONLY: dp
  USE mo_tracdef,      ONLY: GAS, AEROSOL, GAS_OR_AEROSOL, jptrac
  USE mo_exception,    ONLY: message, finish, em_error, em_warn, em_none
 

  IMPLICIT NONE

  PRIVATE

  !--- public member functions

  PUBLIC :: init_splist
  PUBLIC :: new_species
  PUBLIC :: query_species
  PUBLIC :: get_nspec         ! return number of species defined, optionally keyed by submodel name
  PUBLIC :: printspec
 
  !--- Parameters and integer declarations
  
  INTEGER, PARAMETER, PUBLIC :: nmaxspec   = jptrac ! maximum number of species in the model
  INTEGER, PARAMETER, PUBLIC :: nmaxtrspec =  42 ! maximum number of tracers per species

  !   Tracer types
  
  INTEGER, PARAMETER, PUBLIC :: itrnone    = 0   ! No ECHAM tracer
  INTEGER, PARAMETER, PUBLIC :: itrpresc   = 1   ! Prescribed tracer, read from file
  INTEGER, PARAMETER, PUBLIC :: itrdiag    = 2   ! Diagnostic tracer, no transport
  INTEGER, PARAMETER, PUBLIC :: itrprog    = 3   ! Prognostic tracer, transport by ECHAM
  
  !   Number of species 
  INTEGER,            PUBLIC :: nspec            ! number of all  species defined
  INTEGER,            PUBLIC :: naerospec        ! number of species in aerosol phase

  !   Species-Tracer relationship

  INTEGER,            PUBLIC :: spec_idt(nmaxspec, nmaxtrspec)
  INTEGER,            PUBLIC :: spec_ntrac(nmaxspec)

  !   Basic 'species' type common for gas and aerosol.

  TYPE, PUBLIC :: t_species
   
     ! --- basic information --- 
     CHARACTER(LEN=64) :: longname            ! Name of the species
     CHARACTER(LEN=32) :: shortname           ! Shorthand name  
     CHARACTER(LEN=15) :: units               ! Units for 'outside world', usually kg kg-1 
     REAL(dp)          :: moleweight          ! Molecular weight [g mol-1] 
     ! --- transport --- 
     INTEGER           :: itrtype             ! Tracer transport mechanism
     ! --- physical and chemical parameters ---  
     INTEGER           :: nphase              ! phase of specie 
                                              ! && gas &&
     REAL(dp)          :: tdecay              ! Exponential decay rate (simple tracers)
                                              ! e-folding time in seconds
                                              ! will always be applied if /= 0 !!
     REAL(dp)          :: henry(2)            ! Henry's law coefficient and activation energy
     REAL(dp)          :: dryreac             ! Dry reactivity if species is subject to dry deposition 
                                              ! && aerosol &&
     REAL(dp)          :: density             ! Density of solid or liquid phase [g cm-3]
     LOGICAL           :: lwatsol             ! Species dissolves in water
                                              ! NBB THIS IS NOT THE SAME AS M7 'soluble', which means
                                              ! that an M7 mode can take up water. This flag means
                                              ! that the species can physically dissolve in water.
                                              ! Should be set = TRUE for weakly soluble substances
                                              ! and FALSE only for e.g. dust, black carbon
     LOGICAL           :: lelectrolyte        ! Species dissociates into ions in aqueous solution
     INTEGER           :: nion                ! Number of ions if electrolyte
     REAL(dp)          :: osm                 ! Osmotic coefficient 
     REAL(dp)          :: kappa               ! Kappa-Koehler coefficient for the species.  
                                              ! See: Petters and Kreidenweis
                                              ! ACP 7, 1961-1971, 2007 
     ! --- process instructions ---  
     LOGICAL           :: lburden             ! Calculate column burden for this species     
     LOGICAL           :: lrad                ! Corresponding tracer should interact with radiation
     LOGICAL           :: lemis               ! Corresponding tracer should be emitted
     LOGICAL           :: ldrydep             ! Corresponding tracer should be dry deposited
     LOGICAL           :: lwetdep             ! Corresponding tracer should be scavenged  
     INTEGER           :: nbudg               ! Budget calculations
     INTEGER           :: iaerorad            ! Index to aerosol radiation table
     ! --- accounting ---  
     LOGICAL           :: ltrreq              ! Tracer request
     INTEGER           :: idt                 ! Tracer index for gas-phase tracer
     INTEGER, ALLOCATABLE :: iaerocomp(:)     ! Index to aerocomp array. Dimension nmod.
     ! --- in which submodel --- 
     CHARACTER(LEN=15) :: tsubmname           ! the submodel that the tracer belongs to 
  END TYPE t_species

  ! derived variables  
  
  TYPE(t_species), TARGET, PUBLIC :: speclist(nmaxspec)     ! aerosol and gas-phase species lists
                                                            ! (analogy to trlist in mo_tracdef)
  INTEGER, PUBLIC                 :: aero_idx(nmaxspec)     ! indices of aerosol species

  CONTAINS
   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! Set defaults in the species lists
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! initialize before calling init_trlist 
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


  SUBROUTINE init_splist 

  USE mo_kind,         ONLY: dp
  USE mo_physical_constants,    ONLY: amd       ! dry air mass

  nspec     = 0
  naerospec = 0
   
  speclist(:)% longname     = ''
  speclist(:)% shortname    = ''
  speclist(:)% units        = ''
  speclist(:)% moleweight   = amd
  speclist(:)% idt          = -1
  speclist(:)% tdecay       = 0._dp 
  speclist(:)% itrtype      = itrprog
  speclist(:)% nbudg        = 0 
  speclist(:)% lburden      = .TRUE.
  speclist(:)% henry(1)     = 0._dp
  speclist(:)% henry(2)     = 0._dp
  speclist(:)% dryreac      = 0._dp
  speclist(:)% density      = 0._dp
  speclist(:)% lwatsol      = .FALSE.
  speclist(:)% lelectrolyte = .FALSE.
  speclist(:)% nion         = 0
  speclist(:)% osm          = 0._dp
  speclist(:)% kappa        = 0._dp
  speclist(:)% lrad         = .FALSE.
  speclist(:)% lemis        = .FALSE.
  speclist(:)% ldrydep      = .FALSE.
  speclist(:)% lwetdep      = .FALSE. 
  speclist(:)% iaerorad     = 0 
  speclist(:)% ltrreq       = .FALSE.
  speclist(:)% tsubmname    = '' 
  
  aero_idx(:) = 0

  spec_idt(:,:) = -1
  spec_ntrac(:) = 0

  END SUBROUTINE init_splist 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! new_species defines the physical and chemical properties of a chemical species.
!! Normally, each species will also be defined as tracer (new_tracer routine in mo_tracer).
!! The same species can be used in the definition of various tracers (for example
!! different aerosol modes) 
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! <ol>
!! <li>ham_define_tracers
!! <li>moz_define_tracers
!! </ol> 
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
 
  SUBROUTINE new_species(nphase,                                  &
                         longname,                                &
                         shortname,                               &
                         units,                                   &
                         mw,                                      & 
                         tsubmname,                               & 
                         itrtype,                                 &
                         tdecay,                                  &!gas
                         henry,                                   &
                         dryreac,                                 &
                         density,                                 &!aerosol 
                         lwatsol,                                 &
                         lelectrolyte,                            &
                         nion,                                    &
                         osm,                                     &
                         kappa,                                   &
                         lemis,                                   &
                         ldrydep,                                 &
                         lwetdep,                                 & 
                         nbudg,                                   &  
                         lrad,                                    &
                         iaerorad,                                & 
                         lburden,                                 & 
                         ltrreq,                                  & 
                         idx)



  USE mo_exception,      ONLY: message, message_text
 
  
  INTEGER,          INTENT(IN)   :: nphase               ! GAS or AEROSOL
  CHARACTER(LEN=*), INTENT(IN)   :: longname             ! full name (e.g. for diagnostics)
  CHARACTER(LEN=*), INTENT(IN)   :: shortname            ! short name (e.g. variable name in output files) 
  CHARACTER(LEN=*), INTENT(IN)   :: units                ! concentration units, eg 'kg kg-1', 'VMR'
  REAL(dp),         INTENT(IN)   :: mw                   ! Molecular weight (required for conversions)
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tsubmname    ! name of module defining the species
  INTEGER,  INTENT(IN), OPTIONAL :: itrtype              ! apply tracer transport?   
  INTEGER,  INTENT(IN), OPTIONAL :: nbudg                ! Calculate species budget   
  REAL(dp), INTENT(IN), OPTIONAL :: tdecay               ! e-folding time in seconds
                                                         ! will always be applied if /= 0 !!
  REAL(dp), INTENT(IN), OPTIONAL :: henry(2)             ! Henry coefficient
  REAL(dp), INTENT(IN), OPTIONAL :: dryreac              ! dry reactivity coefficient 
  REAL(dp), INTENT(IN), OPTIONAL :: density              ! (aerosol) density
  LOGICAL,  INTENT(IN), OPTIONAL :: lwatsol              ! flag for water soluble (SOA scheme)
  LOGICAL,  INTENT(IN), OPTIONAL :: lelectrolyte         ! flag for dissociation in water
  INTEGER,  INTENT(IN), OPTIONAL :: nion                 ! number of ions formed
  REAL(dp), INTENT(IN), OPTIONAL :: osm                  ! osmotic coefficient
  REAL(dp), INTENT(IN), OPTIONAL :: kappa                ! Kappa-Koehler coefficient
  LOGICAL,  INTENT(IN), OPTIONAL :: lemis                ! flag for emission
  LOGICAL,  INTENT(IN), OPTIONAL :: ldrydep              ! flag for dry deposition
  LOGICAL,  INTENT(IN), OPTIONAL :: lwetdep              ! flag for wet deposition
  LOGICAL,  INTENT(IN), OPTIONAL :: lrad                 ! flag for radiation interaction
  INTEGER,  INTENT(IN), OPTIONAL :: iaerorad             ! index to aerosol radiation table
  LOGICAL,  INTENT(IN), OPTIONAL :: lburden              ! Calculate column burden
  LOGICAL,  INTENT(IN), OPTIONAL :: ltrreq               ! if or not request tracer 
  INTEGER,  INTENT(OUT),OPTIONAL :: idx                  ! index to species list 
  

  !--- Local data
  INTEGER              :: i
  CHARACTER(len=16)    :: cphasenam

  !--- executable procedure
  ! Default return value
  IF (PRESENT(idx))   idx = 0

  ! Check validity of nphase
  SELECT CASE(nphase)
  CASE(GAS)
    cphasenam = 'Gas-phase'
  CASE(AEROSOL)
    cphasenam = 'Aerosol'
  CASE(GAS_OR_AEROSOL)
    cphasenam = 'Gas or Aerosol'
  CASE DEFAULT
    CALL message('new_species', 'Invalid value for nphase. Only GAS or AEROSOL allowed.', level=em_error)
    RETURN
  END SELECT
  
  ! Check if a species of this name has already been defined
  CALL query_species(longname=longname, nphase=nphase, index=i)
!!mgs!! write(0,*) ' new_species: query longname='//longname,' nphase=',nphase,', index=',i
  
  IF (i > 0) THEN
     CALL message('new_species', TRIM(cphasenam)//' species '//TRIM(longname)//' already defined.', &
                  level=em_warn)
     CALL message('', 'Will use existing properties.', level=em_none)
     IF (PRESENT(idx)) idx = i    ! return index
     RETURN
  END IF
 
!!bug fix mgs: added shortname= !  2010-02-01 
  CALL query_species(shortname=shortname, nphase=nphase, index=i)
!!mgs!! write(0,*) ' new_species: query shortname='//shortname,' nphase=',nphase,', index=',i
  
  IF (i > 0) THEN
     CALL message('new_species', TRIM(cphasenam)//' species '//TRIM(shortname)//' already defined.', &
                  level=em_warn)
     CALL message('', 'Will use existing properties.', level=em_none)
     IF (PRESENT(idx)) idx = i    ! return index
     RETURN
  END IF
 
  ! Check if table is full
  IF (nspec >= nmaxspec) THEN
    WRITE (message_text,*) 'Species list full. nmaxspec = ',nmaxspec
    CALL finish('new_species', message_text)
  END IF
     
  ! Increment number of species instances and store data
  nspec    = nspec + 1
  i        = nspec
     
  ! set species properties 
  speclist(i)% nphase       = nphase
  speclist(i)% longname     = longname
  speclist(i)% shortname    = shortname
  speclist(i)% units        = units
  speclist(i)% moleweight   = mw
  IF (PRESENT(tsubmname)) speclist(i)% tsubmname    = tsubmname 
  IF (PRESENT(tdecay))    speclist(i)% tdecay       = tdecay
  IF (PRESENT(nbudg))     speclist(i)% nbudg        = nbudg
  IF (PRESENT(lburden))   speclist(i)% lburden      = lburden
  IF (PRESENT(lrad))      speclist(i)% lrad         = lrad
  IF (PRESENT(lemis))     speclist(i)% lemis        = lemis
  IF (PRESENT(ldrydep))   speclist(i)% ldrydep      = ldrydep
  IF (PRESENT(lwetdep))   speclist(i)% lwetdep      = lwetdep
  IF (PRESENT(ltrreq))    speclist(i)% ltrreq       = ltrreq
  IF (PRESENT(itrtype))   speclist(i)% itrtype      = itrtype
 
!!mgs!! write(0,*) '### new_species: defining new species: index=',i,', shortname=',shortname,', tsubmname=',tsubmname 
  ! for gas only
  IF (PRESENT(dryreac))      speclist(i)% dryreac      = dryreac    ! ## rename to drygamma ???
  IF (PRESENT(henry))        speclist(i)% henry        = henry
  ! for aerosol only 
  IF (PRESENT(density))      speclist(i)% density      = density
  IF (PRESENT(lwatsol))      speclist(i)% lwatsol      = lwatsol
  IF (PRESENT(lelectrolyte)) speclist(i)% lelectrolyte = lelectrolyte
  IF (PRESENT(nion))         speclist(i)% nion         = nion
  IF (PRESENT(osm))          speclist(i)% osm          = osm
  IF (PRESENT(kappa))        speclist(i)% kappa        = kappa
  IF (PRESENT(iaerorad))     speclist(i)% iaerorad     = iaerorad
    
  ! capture species definition errors
  IF (IAND(nphase, AEROSOL) == 0) THEN

     ! capture errors 
     IF (PRESENT(density))      CALL message('new_species',  &
                                  'Property "density" undefined for gas-phase species '//TRIM(shortname),   &
                                  level=em_error)
     IF (PRESENT(lwatsol))      CALL message('new_species',  &
                                  'Property "lwatsol" undefined for gas-phase species '//TRIM(shortname),   &
                                  level=em_error)
     IF (PRESENT(lelectrolyte)) CALL message('new_species',  &
                                'Property "lelectrolyte" undefined for gas-phase species '//TRIM(shortname), &
                                  level=em_error)
     IF (PRESENT(nion))         CALL message('new_species',  &
                                  'Property "nion" undefined for gas-phase species '//TRIM(shortname),   &
                                  level=em_error)
     IF (PRESENT(osm))          CALL message('new_species',  &
                                  'Property "osm" undefined for gas-phase species '//TRIM(shortname),   &
                                  level=em_error)
     IF (PRESENT(kappa))        CALL message('new_species',  &
                                  'Property "kappa" undefined for gas-phase species '//TRIM(shortname),   &
                                  level=em_error)
     IF (PRESENT(iaerorad))     CALL message('new_species',  &
                                  'Property "iaerorad" undefined for gas-phase species '//TRIM(shortname),   &
                                  level=em_error)

  ELSE IF (IAND(nphase, GAS) == 0) THEN

     ! capture errors 
     IF (PRESENT(dryreac))      CALL message('new_species',  &
                                  'Property "dryreac" undefined for aerosol species '//TRIM(shortname),   &
                                  level=em_error)
     IF (PRESENT(henry))        CALL message('new_species',  &
                                  'Property "henry" undefined for aerosol species '//TRIM(shortname),   &
                                  level=em_error)

  END IF

  ! test validity of aerosol species
  IF (IAND(nphase, AEROSOL) /= 0) THEN
     
     CALL check_aerosol_species(i)      
     ! add aerosol species to aero_idx list
     naerospec = naerospec + 1
     aero_idx(naerospec) = i

  END IF

! CALL message('new_species', speclist(i)%shortname)

  ! Note: no warning if phase is neither gas or aerosol. Nothing will happen: quiet ignorance...
  ! idx will be 0 -- this can be used as indicator for an error

  IF (PRESENT(idx)) idx = i    ! return index

 
  END SUBROUTINE new_species 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! Check the metadata for consistency
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! <ol>
!! <li>new_species 
!! </ol> 
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

  SUBROUTINE check_aerosol_species (idx)
  
 
  INTEGER,  INTENT(IN) :: idx                  ! species index

  ! lelectrolyte & nion
  
  IF (speclist(idx)% lelectrolyte) THEN
     IF (speclist(idx)% nion < 2) THEN
        CALL message('new_species', TRIM(speclist(idx)% shortname) //   &
                     ': lelectrolyte=T requires nion > 1.', level=em_error)
     END IF
  ELSE              ! electrolyte=F so nion cannot be accepted unless 0 or 1
     IF (speclist(idx)% nion > 1) THEN
        CALL message('new_species', TRIM(speclist(idx)% shortname) //   &
                     ': nion > 1 requires lelectrolyte=T.', level=em_error)
     END IF
  END IF

  ! kappa & lwatsol
  
  IF (ABS(speclist(idx)% kappa) > 1.e-20_dp) THEN
     IF (.NOT. speclist(idx)% lwatsol) THEN
        CALL message('new_species', TRIM(speclist(idx)% shortname) //   &
                     ': kappa /= 0 requires lwatsol=T.', level=em_error)
     END IF
  END IF

  ! lrad & iaerorad
  
  IF (speclist(idx)% lrad) THEN 
 
!SF Note: I think the following is a pretty much reasonable requirement, even if outside of a
!         HAM-related context. I therefore think this should stay as is.
    IF (speclist(idx)% iaerorad < 0) THEN
       CALL message('new_species', TRIM(speclist(idx)% shortname) //   &
                       ': iradkey requires >= 0', level=em_error) 
    END IF
  
!SF Note: in the process of fixing #360, I wanted to clean up the following, but the simple fix
!         introduces a circular dependency (importing naeroradspec from mo_ham_rad_data, enclosed into
!         #ifdef HAMMOZ / #endif statements).
!         --> I leave that off for the moment...
    IF (speclist(idx)% iaerorad > 7) THEN 
       CALL message('new_species', TRIM(speclist(idx)% shortname) //   &
                       ': iradkey requires <= 7', level=em_error) 
    END IF
  
  END IF

  END SUBROUTINE check_aerosol_species 
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  
!! get properties of a specific species identified by
!! either its longname or shortname. Query by index is ambivalent due to split
!! in gas-phase and aerosol species lists. If you want to know species properties
!! for a given index in speclist, you can get the result directly as 
!! speclist(idx)% <property>.
!! This routine will return ierr/=0 if the species is undefined.
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! <ol>
!! <li>new_species 
!! </ol> 
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

  SUBROUTINE query_species (longname, shortname,            &      ! search criteria
                           index,                           &      ! query results ...
                           nphase, units,                   &
                           mw,                              &
                           tsubmname,                       & 
                           tdecay,                          &
                           henry, dryreac,                  &
                           density,                         &
                           lwatsol, lelectrolyte,           &
                           nion, osm, kappa,                &
                           lrad,                            &
                           lemis, ldrydep, lwetdep,         &
                           iaerorad,                        &
                           itrtype,                         &
                           idt,                             &
                           ierr)

  
  CHARACTER(len=*) ,INTENT(in)  ,OPTIONAL :: longname     ! name of species
  CHARACTER(len=*) ,INTENT(in)  ,OPTIONAL :: shortname    ! subname of species
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nphase       ! GAS or AEROSOL

  INTEGER          ,INTENT(out) ,OPTIONAL :: index        ! index of species as query result
  CHARACTER(len=*) ,INTENT(out) ,OPTIONAL :: units        ! physical units of tracer
  REAL(dp)         ,INTENT(out) ,OPTIONAL :: mw           ! molecular weight
  CHARACTER(len=*) ,INTENT(out) ,OPTIONAL :: tsubmname    ! tsubmname 
  REAL(dp)         ,INTENT(out) ,OPTIONAL :: tdecay       ! e-folding time in seconds
  REAL(dp)         ,INTENT(out) ,OPTIONAL :: henry(2)     ! Henry coefficient
  REAL(dp)         ,INTENT(out) ,OPTIONAL :: dryreac      ! dry reactivity coefficient
  REAL(dp)         ,INTENT(out) ,OPTIONAL :: density      ! (aerosol) density
  LOGICAL          ,INTENT(out) ,OPTIONAL :: lwatsol      ! flag for water soluble (SOA scheme)
  LOGICAL          ,INTENT(out) ,OPTIONAL :: lelectrolyte ! flag for dissociation in water
  INTEGER          ,INTENT(out) ,OPTIONAL :: nion         ! number of ions formed
  REAL(dp)         ,INTENT(out) ,OPTIONAL :: osm          ! osmotic coefficient
  REAL(dp)         ,INTENT(out) ,OPTIONAL :: kappa        ! Kappa-Koehler coefficient
  LOGICAL          ,INTENT(out) ,OPTIONAL :: lrad         ! flag for radiation interaction
  LOGICAL          ,INTENT(out) ,OPTIONAL :: lemis        ! flag for emissions
  LOGICAL          ,INTENT(out) ,OPTIONAL :: ldrydep      ! flag for dry deposition
  LOGICAL          ,INTENT(out) ,OPTIONAL :: lwetdep      ! flag for wet deposition
  INTEGER          ,INTENT(out) ,OPTIONAL :: iaerorad     ! index to aerosol radiation table
  INTEGER          ,INTENT(out) ,OPTIONAL :: itrtype      ! tracer type (transport characteristics)
  INTEGER          ,INTENT(out) ,OPTIONAL :: idt          ! associated tracer index
  INTEGER          ,INTENT(out) ,OPTIONAL :: ierr         ! error return value

  ! --- local variables
  INTEGER             :: i
  INTEGER             :: ispec   ! local species id 
  INTEGER             :: iphase
  INTEGER, PARAMETER  :: OK = 0
  INTEGER, PARAMETER  :: UNDEF = 1

  ! -- set defaults
  ispec = 0
  IF (PRESENT(ierr))         ierr         = UNDEF       ! default error: species not found
  IF (PRESENT(index))        index        = 0
  IF (PRESENT(units))        units        = ''
  IF (PRESENT(mw))           mw           = 0._dp
  IF (PRESENT(tsubmname))    tsubmname    = '' 
  IF (PRESENT(tdecay))       tdecay       = 0._dp
  IF (PRESENT(henry))        henry(:)     = 0._dp
  IF (PRESENT(dryreac))      dryreac      = 0._dp
  IF (PRESENT(density))      density      = 0._dp
  IF (PRESENT(lwatsol))      lwatsol      = .FALSE.
  IF (PRESENT(lelectrolyte)) lelectrolyte = .FALSE.
  IF (PRESENT(nion))         nion         = 0
  IF (PRESENT(osm))          osm          = 0._dp
  IF (PRESENT(kappa))        kappa        = 0._dp
  IF (PRESENT(lrad))         lrad         = .FALSE.
  IF (PRESENT(lemis))        lemis        = .FALSE.
  IF (PRESENT(ldrydep))      ldrydep      = .FALSE.
  IF (PRESENT(lwetdep))      lwetdep      = .FALSE.
  IF (PRESENT(itrtype))      itrtype  = 0
  IF (PRESENT(idt))          idt          = 0
  IF (PRESENT(iaerorad))     iaerorad     = 0

  ! -- test if phase if given
  iphase = -1    ! undefined
  IF (PRESENT(nphase)) THEN
    iphase = nphase
    IF (IAND(iphase, GAS+AEROSOL) == 0) CALL message('query_species', &
                'Invalid value for nphase!', level=em_error)
  END IF

  ! -- search for species: always search for longname first
  IF (PRESENT(longname)) THEN
!!mgs!!    IF (iphase == 1 .OR. iphase == 2 .OR. iphase == 4) THEN
!!mgs!!      DO i=1, nspec
!!mgs!!        IF (TRIM(speclist(i)%longname) == TRIM(longname) .AND. IAND(speclist(i)%nphase,iphase) /= 0) THEN 
!!mgs!!          ispec = i
!!mgs!!          EXIT
!!mgs!!        END IF
!!mgs!!      END DO
!!mgs!!    ELSE
      DO i=1, nspec
        IF (TRIM(speclist(i)%longname) == TRIM(longname)) THEN 
          ispec = i
          EXIT
        END IF
      END DO
!!mgs!!    END IF
  END IF

  ! -- if no species found or longname not given look for shortname
  IF (ispec == 0 .AND. PRESENT(shortname)) THEN
!!mgs!!    IF (iphase == 1 .OR. iphase == 2 .OR. iphase == 4) THEN
!!mgs!!      DO i=1, nspec
!!mgs!!        IF (TRIM(speclist(i)%shortname) == TRIM(shortname) .AND. IAND(speclist(i)%nphase,iphase) /= 0) THEN 
!!mgs!!          ispec = i
!!mgs!!          EXIT 
!!mgs!!        END IF
!!mgs!!      END DO
!!mgs!!    ELSE
      DO i=1, nspec
        IF (TRIM(speclist(i)%shortname) == TRIM(shortname)) THEN 
          ispec = i
          EXIT 
        END IF
      END DO
!!mgs!!    END IF
  END IF

  ! --- Collect species properties
  IF (ispec > 0) THEN
    IF (PRESENT(ierr))         ierr         = OK  ! error status: OK
    IF (PRESENT(index))        index        = ispec
    IF (PRESENT(units))        units        = speclist(ispec)% units
    IF (PRESENT(mw))           mw           = speclist(ispec)% moleweight 
    IF (PRESENT(tsubmname))    tsubmname    = speclist(ispec)% tsubmname  
    IF (PRESENT(tdecay))       tdecay       = speclist(ispec)% tdecay
    IF (PRESENT(henry))        henry(:)     = speclist(ispec)% henry
    IF (PRESENT(dryreac))      dryreac      = speclist(ispec)% dryreac
    IF (PRESENT(density))      density      = speclist(ispec)% density
    IF (PRESENT(lwatsol))      lwatsol      = speclist(ispec)% lwatsol
    IF (PRESENT(lelectrolyte)) lelectrolyte = speclist(ispec)% lelectrolyte
    IF (PRESENT(nion))         nion         = speclist(ispec)% nion
    IF (PRESENT(osm))          osm          = speclist(ispec)% osm
    IF (PRESENT(kappa))        kappa        = speclist(ispec)% kappa
    IF (PRESENT(lrad))         lrad         = speclist(ispec)% lrad
    IF (PRESENT(lemis))        lemis        = speclist(ispec)% lemis   
    IF (PRESENT(ldrydep))      ldrydep      = speclist(ispec)% ldrydep
    IF (PRESENT(lwetdep))      lwetdep      = speclist(ispec)% lwetdep
    IF (PRESENT(itrtype))      itrtype      = speclist(ispec)% itrtype
    IF (PRESENT(idt))          idt          = speclist(ispec)% idt
    IF (PRESENT(iaerorad))     iaerorad     = speclist(ispec)% iaerorad
  END IF

  END SUBROUTINE query_species


!! get_nspec: return number of specers defined, optionally keyed by submodel name

  FUNCTION get_nspec (tsubmname)  RESULT (ispec)

    INTEGER                                  :: ispec
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: tsubmname

    INTEGER            :: jt

    ispec = 0

    IF (PRESENT(tsubmname)) THEN
      DO jt = 1,nspec
        IF (speclist(jt)%tsubmname == tsubmname) ispec = ispec + 1
      END DO
    ELSE
      ispec = nspec       ! return number of all specers
    END IF

  END FUNCTION

  SUBROUTINE printspec
  !! print species definitions
  !! and fill spec_idt matrix

    USE mo_exception,      ONLY: message, message_text, em_info, em_param
    USE mo_tracdef,        ONLY: ntrac, trlist
    USE mo_util_string,    ONLY: separator

    INTEGER               :: jt, ispec, i
    CHARACTER(len=11)    :: cstring

    CALL message('', 'Species definitions:', level=em_info) 
    CALL message('','',level=em_param)
    WRITE(message_text,'(a)') 'Id shortname         phase       submodel         '//&
                              'mol.wght unit         itrtype idt longname'
    CALL message('', message_text, level=em_param)
    DO jt = 1, nspec

      SELECT CASE(speclist(jt)%nphase)
      CASE (GAS)
        cstring = 'gas        '
      CASE (AEROSOL)
        cstring = 'aerosol    '
      CASE (GAS_OR_AEROSOL)
        cstring = 'gas|aerosol'
      END SELECT
      
!SF Note: in the line below, all string fields are meant to be printed with their original declared length
!         (which should make sense!), except for the shortname property (otherwise these lines are
!         extremely long and unreadable). Therefore it may be that the shortname property is truncated in the
!         printed lines (> 18 characters...)
      WRITE(message_text,'(i0,t4,a,t22,2(a,1x),f7.3,3x,a,1x,i3,1x,i3,2x,a)') &
                               jt,TRIM(speclist(jt)%shortname), cstring, speclist(jt)%tsubmname, &
                               speclist(jt)%moleweight, speclist(jt)%units, speclist(jt)%itrtype, &
                               speclist(jt)%idt, speclist(jt)%longname

      CALL message('', message_text, level=em_param)

    END DO

    CALL message('','',level=em_param)
    WRITE(message_text,'(a)') 'Notes:'
    CALL message('', message_text, level=em_param)
    WRITE(message_text,'(a)') '1. itrtype meaning: '
    CALL message('', message_text, level=em_param)
    WRITE(message_text,'(a,i0,a)') ' --> ',itrnone,' = No ECHAM tracer, '
    CALL message('', message_text, level=em_param)
    WRITE(message_text,'(a,i0,a)') ' --> ',itrpresc,' = Prescribed tracer, read from file'
    CALL message('', message_text, level=em_param)
    WRITE(message_text,'(a,i0,a)') ' --> ',itrdiag,' = Diagnostic tracer, no transport'
    CALL message('', message_text, level=em_param)
    WRITE(message_text,'(a,i0,a)') ' --> ',itrprog,' = Prognostic tracer, transport by ECHAM'
    CALL message('', message_text, level=em_param)
    CALL message('','',level=em_param)
    WRITE(message_text,'(a)') '2. idt is the tracer index for gas-phase tracers only'
    CALL message('', message_text, level=em_param)

    !! fill spec_idt matrix
    DO jt = 1, ntrac
      ispec = trlist%ti(jt)%spid 
      IF (ispec > 0) THEN
        spec_ntrac(ispec) = spec_ntrac(ispec) + 1
        spec_idt(ispec, spec_ntrac(ispec)) = jt
      END IF
    END DO

    CALL message('','',level=em_param)
    CALL message('', 'Species-tracer matrix:', level=em_param)
    CALL message('', 'Species index : Tracer index/indices', level=em_param)
    !>>dod removed use of list-directed I/O to variable
    message_text=''
    DO jt = 1, nspec
       WRITE(message_text, '(i0,t15,a1)') jt,':'
       DO i=1,spec_ntrac(jt)
          WRITE(cstring,'(i0,1x)') spec_idt(jt,i)
          WRITE(message_text,'(a,1x,a)') TRIM(message_text), TRIM(cstring)
       END DO
      CALL message('', message_text, level=em_param)
    END DO
    !<<dod
    CALL message('','',level=em_param)
    CALL message('',separator)

  END SUBROUTINE printspec

END MODULE mo_species

