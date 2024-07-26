!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief mo_submodel_diag: provide easy access to tracer-dependent process diagnostics
!!
!! @remarks
!! This module defines the t_diag and t_diag_list structure and provides routines to manage these.
!! Example use: see mo_hammoz_drydep, mo_hammoz_emi, mo_hammoz_sedimentation and mo_hammoz_wetdep
!!
!! @author Martin Schultz, FZ Juelich (2009-09-22)
!!
!! $ID: n/a$
!!
!! @par Origin
!! Created after a concept from Declan O'Donnell, MPI.
!!
!! @par This module is used by
!! various submodels
!!
!! @par Notes:
!! The present implementation uses a loop 1,ndiag to find the diagnostics associated with a specific
!! tracer. May be more efficient to define a map_index array (1:ntrac) which can be populated the first time
!! when get_diag_* is called.
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!! @par Copyright
!! 2009 by FZ-Juelich and MPI-M
!! This software is provided for non-commercial use only.
!!
!
MODULE mo_submodel_diag

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  USE mo_tracdef,       ONLY: ln           ! character length


  IMPLICIT NONE

  PRIVATE

  SAVE

  PUBLIC  :: t_diag, t_diag_list                  ! diagnostic structures
  PUBLIC  :: new_diag_list                        ! subroutine for creating a new diagnostic list
  PUBLIC  :: new_diag                             ! subroutine for defining diagnostic entries
  PUBLIC  :: new_diag_element                     ! subroutine for defining a single diagnostic entry
  PUBLIC  :: get_diag_pointer                     ! subroutine to retrieve diagnostics
  PUBLIC  :: get_diag_index                       ! return index of a specific diag in diag_list
  PUBLIC  :: get_diag_stream                      ! return stream pointer for a diagnostics
  PUBLIC  :: vmem2d, vmem3d
  PUBLIC  :: OFF, ON, BYTRACER, BYSPECIES, BYMODE, BYNUMMODE, BYUSER  ! key types

! Interface
  INTERFACE get_diag_pointer
    MODULE PROCEDURE get_diag_pointer2d
    MODULE PROCEDURE get_diag_pointer3d
  END INTERFACE 

! Variable declarations

  TYPE :: t_diag                         ! 2-D or 3-D output field plus meta data
     INTEGER            :: key_type      ! Output type: 
                                         !   1 = by tracer
                                         !   2 = by species
                                         !   3 = by aerosol mode or bin (mass)
                                         !   4 = by aerosol mode or bin (number)
                                         !   5 = by anything else
     INTEGER            :: key           ! index to tracer, species, mode etc.
     INTEGER            :: ndims         ! number of dimensions for this diag
     CHARACTER (LEN=ln) :: name          ! name of the diagnostic field
     REAL(dp), POINTER  :: fld2d(:,:)
     REAL(dp), POINTER  :: fld3d(:,:,:)
  END TYPE t_diag


  TYPE t_diag_list                       ! list of t_diag elements
     INTEGER            :: ndims         ! number of dimensions fo rthe entire list (2 or 3)
     INTEGER            :: ndiag         ! number of diagnostic entries
     INTEGER            :: nmaxdiag      ! max. number of allowed diag entries (automatic!)
     INTEGER            :: nmaxkey(5)    ! maximum allowed value for each key_type
     CHARACTER (len=ln) :: tsubmname     ! name of submodel
     CHARACTER (len=ln) :: diag_list_name
     CHARACTER (len=64) :: longname      ! default variable longname template
     CHARACTER (len=64) :: units         ! default variable units
     LOGICAL            :: lpost         ! include in output
     TYPE(t_stream), POINTER    :: tstream       ! stream reference
     TYPE(t_diag), ALLOCATABLE  :: diag(:)
  END TYPE

  !! pointer lists for diagnostics not keyed on tracers
  TYPE :: vmem2d
     REAL(dp), POINTER  :: ptr(:,:)
  END TYPE vmem2d

  TYPE :: vmem3d
     REAL(dp), POINTER  :: ptr(:,:,:)
  END TYPE vmem3d

  INTEGER, PARAMETER    :: OFF       = 0,  &
                           ON        = 1,  &
                           BYTRACER  = 1,  &
                           BYSPECIES = 2,  &
                           BYMODE    = 3,  &
                           BYNUMMODE = 4,  &
                           BYUSER    = 5

  CONTAINS


!
!! @brief new_diag_list: define a new diagnostics type
!!
!! @remarks
!! This subroutine initializes a variable of type t_diag_list and associates the list
!! with an output stream. The user can either provide a pointer to an existing stream 
!! or provide a name for a new stream which will then be created.
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!

  SUBROUTINE new_diag_list(diag_list, tstream, diagname, tsubmname,       &
                           longname, units, ndims, nmaxkey, lpost)

  USE mo_exception,     ONLY: finish, message, message_text, em_info, em_error
  USE mo_tracdef,       ONLY: ntrac       ! for default list size
  USE mo_species,       ONLY: nspec       ! for default list size

  ! -- parameters
  TYPE(t_diag_list),           INTENT(inout)   :: diag_list
  TYPE(t_stream),    POINTER,  INTENT(inout)   :: tstream    ! pointer reference to existing stream
  CHARACTER(len=*),  OPTIONAL, INTENT(in)      :: diagname   ! name of the diagnostic type
  CHARACTER(len=*),  OPTIONAL, INTENT(in)      :: tsubmname  ! name of the submodel managing this diag_list
  CHARACTER(len=*),  OPTIONAL, INTENT(in)      :: longname   ! descriptive name of diagnostics: "of XY" wil be added
  CHARACTER(len=*),  OPTIONAL, INTENT(in)      :: units      ! physical unit for output
  INTEGER,           OPTIONAL, INTENT(in)      :: ndims      ! dimensionality of diagnostics (2 or 3-D)
  INTEGER,           OPTIONAL, INTENT(in)      :: nmaxkey(5) ! maximum number of key values for each of
                                                             ! the five key types
  LOGICAL,           OPTIONAL, INTENT(in)      :: lpost      ! include in output (default=true)

  ! -- local variables
  CHARACTER(len=ln)                :: cdiagname, csubmname
  CHARACTER(len=64)                :: clongname, cunits
  INTEGER                          :: idims, imaxkey(5)
  INTEGER                          :: imaxdiag

  ! -- executable code

  ! 1) -- set defaults
  cdiagname   = 'unknown'
  csubmname   = 'unknown'
  clongname   = ''
  cunits      = ''
  idims       = 2          ! set 2-D as default
  imaxkey     = (/ ntrac, nspec, 0, 0, 0 /)

  ! 1) -- check input arguments and set local variables
  
  IF (PRESENT(diagname))   cdiagname   = diagname
  IF (PRESENT(tsubmname))  csubmname   = tsubmname 
  IF (PRESENT(longname))   clongname   = longname
  IF (PRESENT(units))      cunits      = units
  IF (PRESENT(ndims))      idims       = ndims
  IF (PRESENT(nmaxkey))    imaxkey     = nmaxkey

  ! -- error checks
  IF (idims /= 2 .AND. idims /= 3) THEN
    CALL message('new_diag_list', TRIM(cdiagname)//': Invalid value of ndims. Only 2 or 3 allowed!', &
                 level=em_error)
    idims = 2      ! this should get the code through the initialisation after which it stops
  END IF

  imaxdiag = MAXVAL(imaxkey)     ! maximum allowed number of diagnostics
  !SFnote: I don't quite understand why 'MAXVAL(imaxkey)' and not 'SUM(imaxkey)' (in the hypothesis
  !         of allowing several diag types at the same time)
  !        Since #299, I've modified the diagnostics to not allow BYMODE diags if BYTRACER or BYSPECIES
  !        are chosen. This way, using MAXVAL(imaxkey) is logical again, but this will
  !        will have to be modified if one want to implement multiple diag types at the same time.
  IF (imaxdiag > ntrac) THEN
    write(message_text,*) 'cdiagname=',TRIM(cdiagname),                    &
                          ': Invalid value in nmaxval (=',imaxdiag,')! ',  &
                          'Only up to ntrac (=',ntrac,') allowed!'
    CALL message('new_diag_list', message_text, level=em_error)
  END IF

  ! 2) -- associate diagnostics with existing stream or create new stream
  IF (ASSOCIATED(tstream)) THEN 
    !! ++mgs: error check necessary if stream is actually allocated !?
    CALL message('new_diag_list', 'Defining diagnostics type '//TRIM(cdiagname)//' using stream ' &
                 //TRIM(tstream%name), level=em_info)
  ELSE
    CALL finish('new_diag_list', 'Routine requires a valid stream reference!')
  END IF

  ! 3) -- set t_diag_list structure fields
  diag_list%ndims          = idims
  diag_list%ndiag          = 0
  diag_list%nmaxdiag       = imaxdiag
  diag_list%nmaxkey        = imaxkey
  diag_list%tsubmname      = csubmname
  diag_list%longname       = clongname
  diag_list%units          = cunits
  diag_list%diag_list_name = cdiagname
  diag_list%tstream        => tstream
  diag_list%lpost          = .true.
  IF (PRESENT(lpost)) THEN
    diag_list%lpost        = lpost
  END IF
  
  ! 4) -- allocate diag field

  IF (ALLOCATED(diag_list%diag)) CALL finish('new_diag_list', 'diag field already allocated!')
  ALLOCATE (diag_list%diag(imaxdiag))

  END SUBROUTINE new_diag_list    



!
!! @brief new_diag_element: this subroutine adds one tracer-, species- or mode-specific diagnostics 
!! to the diag_list
!!
!! @remarks
!! This subroutine increases the number of diagnostics defined in the diag_list and populates the
!! respective diag structure fields. A stream element is added to the stream referenced in the
!! diag_list.
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!

  SUBROUTINE new_diag_element(diag_list, diagname, ikey_type, ikey,     &
                              longname, units, lpost)

  USE mo_exception,     ONLY: message, message_text, em_info, em_error
  USE mo_memory_base,   ONLY: add_stream_element

  ! -- parameters
  TYPE(t_diag_list),          INTENT(inout)   :: diag_list
  CHARACTER(len=*),           INTENT(in)      :: diagname   ! name of the diagnostic variable
  INTEGER,                    INTENT(in)      :: ikey_type  ! key type (1 to 4)
  INTEGER,                    INTENT(in)      :: ikey       ! tracer-, species- or mode-number
  CHARACTER(len=*), OPTIONAL, INTENT(in)      :: longname   ! descriptive name of diagnostics
  CHARACTER(len=*), OPTIONAL, INTENT(in)      :: units      ! physical unit for output
  LOGICAL,          OPTIONAL, INTENT(in)      :: lpost      ! add diagnostics to file output

  ! -- local variables
  INTEGER               :: i, jt
  LOGICAL               :: lok, lopost
  CHARACTER(len=64)     :: clongname, cunits
  CHARACTER(len=32)     :: cdebug, cdebug2

  ! 1) -- error checks
  lok = .true. 
  IF (ikey_type < 1 .OR. ikey_type > 5) THEN
    CALL message('new_diag_element', diag_list%diag_list_name//': Invalid value of key_type (1..5)!', &
                 level=em_error)
    lok = .false.
  ELSE IF (ikey < 1 .OR. ikey > diag_list%nmaxkey(ikey_type)) THEN
    write (message_text, *) diag_list%diag_list_name//': Invalid key value ',ikey,  &
                            '. For key_type ',ikey_type, &
                            ' the maximum allowed value is ',diag_list%nmaxkey(ikey_type)
    CALL message('new_diag_element', message_text, level=em_error)
    lok = .false.
  END IF
  IF (diag_list%ndiag == diag_list%nmaxdiag) THEN
    write (message_text, *) 'Diagnostics list '//diag_list%diag_list_name        &
                            //' full. Maximum allowed number of diag entries is ',diag_list%nmaxdiag
    CALL message('new_diag_element', message_text, level=em_error)
    CALL message('', 'Check the nmaxkey values in the call to new_diag_element_list.', level=em_info)
    lok = .false.
  END IF

  !>>SF check for duplicate diags
  DO jt=1,diag_list%ndiag
     IF (TRIM(diag_list%diag(jt)%name) == TRIM(diagname)) THEN
        lok = .FALSE.
        EXIT
     ENDIF
  ENDDO
  !<<SF

  IF (.NOT. lok) RETURN

  !-- optional arguments
  IF (PRESENT(longname)) THEN
    clongname = TRIM(longname)
  ELSE
    clongname = TRIM(diag_list%longname)//' of '//TRIM(diagname)
  END IF
  IF (PRESENT(units)) THEN
    cunits = units
  ELSE
    cunits = diag_list%units
  END IF

  ! 2) -- increase number of diagnostics and populate diag element
  diag_list%ndiag = diag_list%ndiag + 1
  i = diag_list%ndiag
  diag_list%diag(i)%name     = diagname
  diag_list%diag(i)%key_type = ikey_type
  diag_list%diag(i)%key      = ikey

  ! 3) -- add stream element and associate pointer
  lopost = diag_list%lpost              ! set default
  IF (PRESENT(lpost))  lopost = lpost   ! overwrite if given explicitly
!### debug
    IF (ikey_type == BYTRACER)  cdebug = 'BYTRACER'
    IF (ikey_type == BYSPECIES) cdebug = 'BYSPECIES'
    IF (ikey_type == BYMODE)    cdebug = 'BYMODE'
    IF (ikey_type == BYNUMMODE) cdebug = 'BYNUMMODE'
    IF (ikey_type == BYUSER)    cdebug = 'BYUSER'
    cdebug2 = ', lpost=.true. (default)'
    IF (PRESENT(lpost)) THEN
      IF (lpost) THEN
         WRITE(cdebug2,*) ', lpost=.true.'
      ELSE
         WRITE(cdebug2,*) ', lpost=.false.'
      END IF
    END IF
  IF (diag_list%ndims == 2) THEN
    WRITE(message_text, '(A,A24,1x,4(A))') 'defining new 2D diagnostics ', TRIM(diagname), &
                        TRIM(cdebug), ': units=', TRIM(cunits), TRIM(cdebug2)
    CALL message('', message_text, level=em_info)
    CALL add_stream_element(diag_list%tstream, diagname, diag_list%diag(i)%fld2d, units=cunits, &
                            longname=clongname, lpost=lopost )
    diag_list%diag(i)%ndims = 2
  ELSE
    WRITE(message_text, '(A,A24,1x,4(A))') 'defining new 3D diagnostics ', TRIM(diagname), &
                        TRIM(cdebug), ': units=', TRIM(cunits), TRIM(cdebug2)
    CALL message('', message_text, level=em_info)
    CALL add_stream_element(diag_list%tstream, diagname, diag_list%diag(i)%fld3d, units=cunits, &
                            longname=clongname, lpost=lopost )
    diag_list%diag(i)%ndims = 3
  END IF

  END SUBROUTINE new_diag_element



!
!! @brief new_diag: this subroutine populates a diag_list with a given array of tracer-, species- 
!! or mode-specific diagnostics 
!!
!! @remarks
!! This subroutine loops over a given valflag array and adds individual diagnostic elements to
!! the diag_list. Longname and units will overwrite the diag_list defaults if provided.
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!

  SUBROUTINE new_diag(diag_list, nval, valflag, valname, ikey_type,     &
                      longname, units, lpost)

  USE mo_exception,     ONLY: message, em_error

  ! -- parameters
  TYPE(t_diag_list),          INTENT(inout)   :: diag_list
  INTEGER,                    INTENT(in)      :: nval          ! number of diag elements to be created
  LOGICAL,                    INTENT(in)      :: valflag(nval) ! decides whether to create diag element
  CHARACTER(len=*),           INTENT(in)      :: valname(nval) ! name of the diagnostic type
  INTEGER,                    INTENT(in)      :: ikey_type     ! key type (1 to 4)
  CHARACTER(len=*), OPTIONAL, INTENT(in)      :: longname      ! descriptive name of diagnostics
  CHARACTER(len=*), OPTIONAL, INTENT(in)      :: units         ! physical unit for output
  LOGICAL,          OPTIONAL, INTENT(in)      :: lpost         ! add diagnostics to file output

  ! -- local variables
  INTEGER               :: jt
  LOGICAL               :: lok, lopost
  CHARACTER(len=64)     :: clongname, cunits

  ! 1) -- error checks
  lok = .true. 
  IF (nval <= 0) THEN
    CALL message('new_diag', diag_list%diag_list_name//': No flag and name values provided!', &
                 level=em_error)
    lok = .false.
  END IF
  IF (ikey_type < 1 .OR. ikey_type > 5) THEN
    CALL message('new_diag', diag_list%diag_list_name//': Invalid value of key_type (1..5)!', &
                 level=em_error)
    lok = .false.
  END IF
  !! Abort processing in case of errors
  IF (.NOT. lok) RETURN

  ! 2) -- handle optional arguments
  IF (PRESENT(longname)) THEN
    clongname = longname
  ELSE
    clongname = diag_list%longname
  END IF
  IF (PRESENT(units)) THEN
    cunits = units
  ELSE
    cunits = diag_list%units
  END IF
  IF (PRESENT(lpost)) THEN
    lopost = lpost
  ELSE
    lopost = diag_list%lpost
  END IF

  ! 3) -- add diag elements for each valname that has valflag=.true.
  DO jt = 1,nval
    IF (valflag(jt)) CALL new_diag_element(diag_list,                                              &
                                           TRIM(diag_list%diag_list_name)//'_'//TRIM(valname(jt)), &
                                           ikey_type, jt,                                          &
                                           longname=TRIM(clongname)//' of '//TRIM(valname(jt)),    & 
                                           units=TRIM(cunits), lpost=lopost)
  END DO

  END SUBROUTINE new_diag


!
!! @brief get_diag_index: return the index number in the diag_list that is associated with a given
!!  tracer id.
!!
!! @remarks
!! This function will always return the first occurence of a matching diagnostics. The user is
!! responsible to ensure that new_diag entries don't lead to duplicate matches.
!! If any new_diag defines a key_type of 4, then ivalue must be given as argument to this function.
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!

  FUNCTION get_diag_index(diag_list, ldebug, itrac, ivalue)    RESULT (index)

  USE mo_exception,        ONLY: finish, message, message_text, em_debug
  USE mo_tracdef,          ONLY: trlist, AEROSOLMASS, AEROSOLNUMBER

  INTEGER                         :: index
  TYPE(t_diag_list), INTENT(in)   :: diag_list
  LOGICAL,           INTENT(in)   :: ldebug       ! print debugging info
  INTEGER,           INTENT(in)   :: itrac        ! tracer index in trlist
  INTEGER, OPTIONAL, INTENT(in)   :: ivalue       ! key value to look for if key_type==BYUSER

  ! -- local variables
  INTEGER            :: i, iphase, ispec, iclass
  LOGICAL            :: lfound

  ! 1) -- collect information
  index  = -1
  lfound = .false.
  iphase = trlist%ti(itrac)%nphase
  ispec  = trlist%ti(itrac)%spid
  iclass  = trlist%ti(itrac)%mode

  IF (ldebug) THEN
    WRITE(message_text,*) 'itrac, iphase, ispec, iclass = ',itrac, iphase, ispec, iclass
    CALL message('get_diag_index', message_text, level=em_debug)
  END IF

  ! 2) -- find matching diag element
  DO i = 1, diag_list%ndiag
    SELECT CASE(diag_list%diag(i)%key_type) 
    CASE(1)       ! BYTRACER : index = tracer number
      IF (diag_list%diag(i)%key == itrac) THEN
        lfound = .true.
        index = i
        IF (ldebug) THEN
          WRITE(message_text,*) 'itrac: lfound=.true. for diag = ',i,'. diag_list%diag(i)%key_type = ', &
                                diag_list%diag(i)%key_type
          CALL message('', message_text, level=em_debug)
        END IF
      END IF
    CASE(2)       ! BYSPECIES : index = species number
!### ispec == 0 shouldn't occur...
      IF (ispec > 0 .AND. diag_list%diag(i)%key == ispec) THEN
        lfound = .true.
        IF (ldebug) THEN
          WRITE(message_text,*) 'ispec: lfound=.true. for diag = ',i,'. diag_list%diag(i)%key_type = ', &
                                diag_list%diag(i)%key_type
          CALL message('', message_text, level=em_debug)
        END IF
        index = i
      END IF
    CASE(3)       ! BYMODE : index = mode
      IF (iclass > 0 .AND. diag_list%diag(i)%key == iclass .AND. iphase == AEROSOLMASS) THEN
        lfound = .true.
        index = i
        IF (ldebug) THEN
          WRITE(message_text,*) 'iclass: lfound=.true. for diag = ',i,'. diag_list%diag(i)%key_type = ', &
                                diag_list%diag(i)%key_type
          CALL message('', message_text, level=em_debug)
        END IF
      END IF
    CASE(4)       ! BYNUMMODE : index = mode
      IF (iclass > 0 .AND. diag_list%diag(i)%key == iclass .AND. iphase == AEROSOLNUMBER) THEN
        lfound = .true.
        index = i
        IF (ldebug) THEN
          WRITE(message_text,*) 'iclass (NUM): lfound=.true. for diag = ',i,'. diag_list%diag(i)%key_type = ', &
                                diag_list%diag(i)%key_type
          CALL message('', message_text, level=em_debug)
        END IF
      END IF
    CASE(5)       ! BYUSER : index = something else, user-defined
      IF (.NOT. PRESENT(ivalue)) CALL finish('get_diag_index', 'ivalue required but not given!')
      IF (diag_list%diag(i)%key == ivalue) THEN
        lfound = .true.
        index = i
      END IF
    END SELECT
    IF (lfound) EXIT          ! first occurence of diag takes precedence! (usually BYTRACER)
  END DO
  ! IF (ldebug .AND. .NOT. lfound) THEN
    ! WRITE(message_text,*) 'No diag found for tracer itrac = ',itrac
    ! CALL message('', message_text, level=em_debug)
  ! END IF

  END FUNCTION get_diag_index

!
!! @brief get_diag_pointer: return the 2-D or 3-D pointer of the diagnostics associated with the
!! tracer number itrac or the given ivalue for key_type=4.
!!
!! @remarks
!! This function first calls get_diag_index and then returns the 2-D or 3-D pointer depending
!! on the shape of the ptr argument 
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!

  SUBROUTINE get_diag_pointer2d(diag_list, ptr2d, itrac, ivalue, ierr, ldebug)

  USE mo_exception,        ONLY: finish, message_text

  TYPE(t_diag_list), INTENT(in)   :: diag_list
  INTEGER,           INTENT(in)   :: itrac       ! tracer index in trlist
  INTEGER, OPTIONAL, INTENT(in)   :: ivalue      ! key value to look for if key_type==BYUSER
  INTEGER, OPTIONAL, INTENT(out)  :: ierr        ! error value (code will abort in case of error if not present)
  LOGICAL, OPTIONAL, INTENT(in)   :: ldebug
  REAL(dp), POINTER, INTENT(out)  :: ptr2d(:,:)

  INTEGER, PARAMETER :: OK=0, LIST_EMPTY=1, NOT_FOUND=2
  INTEGER            :: ind, myierr
  LOGICAL            :: lodebug


  NULLIFY (ptr2d)
  myierr = OK

  lodebug = .FALSE.
  IF (PRESENT(ldebug)) THEN
    lodebug = ldebug
  END IF

  ! 1) -- check if diag_list is empty
  IF (diag_list%ndiag == 0) THEN
    myierr = LIST_EMPTY
  ELSE

    ! 2) -- find diag index
    IF (PRESENT(ivalue)) THEN
      ind = get_diag_index(diag_list, lodebug, itrac, ivalue) 
    ELSE
      ind = get_diag_index(diag_list, lodebug, itrac) 
    END IF

    IF (ind <= 0) THEN
      IF (PRESENT(ivalue)) THEN
        WRITE (message_text, *) 'Cannot find index for diag '//TRIM(diag_list%diag_list_name), &
                                ': itrac, ivalue = ', itrac, ivalue
      ELSE 
        WRITE (message_text, *) 'Cannot find index for diag '//TRIM(diag_list%diag_list_name), &
                                ': itrac = ', itrac
      END IF
      myierr = NOT_FOUND
      IF (.NOT. PRESENT(ierr)) CALL finish('get_diag_pointer', message_text)
    END IF
  END IF

  IF (myierr == OK) ptr2d => diag_list%diag(ind)%fld2d
  IF (PRESENT(ierr)) ierr = myierr

  END SUBROUTINE get_diag_pointer2d   

  ! --

  SUBROUTINE get_diag_pointer3d(diag_list, ptr3d, itrac, ivalue, ierr, ldebug)

  USE mo_exception,        ONLY: finish, message_text

  TYPE(t_diag_list), INTENT(in)   :: diag_list
  INTEGER,           INTENT(in)   :: itrac       ! tracer index in trlist
  INTEGER, OPTIONAL, INTENT(in)   :: ivalue      ! key value to look for if key_type==4
  INTEGER, OPTIONAL, INTENT(out)  :: ierr        ! error value (code will abort in case of error if not present)
  LOGICAL, OPTIONAL, INTENT(in)   :: ldebug
  REAL(dp), POINTER, INTENT(out)  :: ptr3d(:,:,:)

  INTEGER, PARAMETER :: OK=0, LIST_EMPTY=1, NOT_FOUND=2
  INTEGER            :: ind, myierr
  LOGICAL            :: lodebug


  NULLIFY (ptr3d)
  myierr = OK

  lodebug = .FALSE.
  IF (PRESENT(ldebug)) THEN
    lodebug = ldebug
  END IF

  ! 1) -- check if diag_list is empty
  IF (diag_list%ndiag == 0) THEN
    myierr = LIST_EMPTY
  ELSE

    ! 2) -- find diag index
    IF (PRESENT(ivalue)) THEN
      ind = get_diag_index(diag_list, lodebug, itrac, ivalue)
    ELSE
      ind = get_diag_index(diag_list, lodebug, itrac)
    END IF

    IF (ind <= 0) THEN
      IF (PRESENT(ivalue)) THEN
        WRITE (message_text, *) 'Cannot find index for diag '//diag_list%diag_list_name, &
                                ': itrac, ivalue = ', itrac, ivalue
      ELSE 
        WRITE (message_text, *) 'Cannot find index for diag '//diag_list%diag_list_name, &
                                ': itrac = ', itrac
      END IF
      myierr = NOT_FOUND
      IF (.NOT. PRESENT(ierr)) CALL finish('get_diag_pointer', message_text)
    END IF
  END IF

  IF (myierr == OK) ptr3d => diag_list%diag(ind)%fld3d
  IF (PRESENT(ierr)) ierr = myierr

  END SUBROUTINE get_diag_pointer3d   



  SUBROUTINE get_diag_stream(diag_list, tstream)

  USE mo_exception,        ONLY: finish
  USE mo_linked_list,      ONLY: t_stream

  TYPE(t_diag_list),       INTENT(in)   :: diag_list
  TYPE(t_stream), POINTER, INTENT(out)  :: tstream

  NULLIFY (tstream)
  IF (ASSOCIATED(diag_list%tstream)) THEN
    tstream => diag_list%tstream
  ELSE
    CALL finish('get_diag_stream', 'tstream pointer in '//diag_list%diag_list_name//' not associated!')
  END IF

  END SUBROUTINE get_diag_stream


END MODULE mo_submodel_diag

