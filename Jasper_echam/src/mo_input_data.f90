!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! This module contains procedures for mo_input dealing with data lists 
! Never directly use this module from outside the mo_input framework - only via mo_input.f90
MODULE mo_input_data
USE mo_input_strings
USE mo_input_types
IMPLICIT NONE

PUBLIC :: InputDataNew, InputDataDone, InputDataUnuse, InputDataDisp, InputDataGetRefGlobal, InputDataGetRef
PUBLIC :: InputDataGetSpecified, InputDataIntrpAlloc, ReplaceAuxData1D, GetNextUserDataNumber

CONTAINS

! ----------------------------------------------------
!  Creating, adding, setting, modifying, displaying:
!     Data buffers
! ----------------------------------------------------

! Creates a new buffer element with default values
  SUBROUTINE InputDataNew(buffer)
  TYPE (input_data_list), POINTER :: buffer

  TYPE (input_data_list), POINTER :: new

    ALLOCATE(new)
    NULLIFY(new%next,new%prev,new%dta)
    new%time(:) = -HUGE(i8)
    new%flags   =  0

    buffer => new

  END SUBROUTINE InputDataNew

! Returns pointer and evt. number of the global buffer associated with the same data as the given buffer
  FUNCTION InputDataGetRefGlobal(ref,no)
  TYPE (input_data_list), POINTER :: InputDataGetRefGlobal
  TYPE (input_data_list), POINTER :: ref
  INTEGER, OPTIONAL,  INTENT(out) :: no

  TYPE (input_data_list), POINTER :: curr
  INTEGER :: i

    i = 1
    curr => AllInputData
    DO WHILE (ASSOCIATED(curr))
      IF (ASSOCIATED(curr%dta,ref%dta)) THEN
        IF (PRESENT(no)) no = i
        InputDataGetRefGlobal => curr
        RETURN
      ENDIF
      i = i + 1
      curr => curr%next
    ENDDO
    IF (PRESENT(no)) no = -1
    NULLIFY(InputDataGetRefGlobal)

  END FUNCTION InputDataGetRefGlobal

! Returns pointer to the data of an available buffer with at least the requested size (may need some more advanced search criteria)
! The returned buffer may be larger along the last (slowest varying) dimension
  FUNCTION InputDataGetRef(dims,lShape,ref,no)
  REAL (dp),            POINTER  :: InputDataGetRef(:,:,:,:,:,:,:)
  INTEGER, OPTIONAL, INTENT(in)  :: dims(:)
  LOGICAL, OPTIONAL, INTENT(in)  :: lShape
  TYPE (input_data_list), OPTIONAL, POINTER :: ref
  INTEGER, OPTIONAL, INTENT(out) :: no

  TYPE (input_data_list), POINTER :: curr, prev
  INTEGER :: nDims, Sz(7), i

    NULLIFY(prev)
    nDims = SIZE(dims)
    i = 1
    curr => AllInputData
    DO WHILE (ASSOCIATED(curr))
      IF (ASSOCIATED(curr%dta)) THEN
        Sz = SHAPE(curr%dta)
      ELSE
        Sz(:) = 0
      ENDIF
      IF (ALL(Sz(1:nDims-1)==dims(1:nDims-1)) .AND. (Sz(nDims) >= dims(nDims)) .AND. (IAND(curr%flags,INPUT_DATA_IN_USE)==0)) THEN
        IF (PRESENT(ref)) ref => curr
        curr%flags = curr%flags + 1 ! Warning: No overflow check!!! 
        InputDataGetRef => curr%dta
        RETURN
      ENDIF
      IF (PRESENT(lShape)) THEN
        IF (.NOT. lShape) THEN
          IF ((SIZE(curr%dta)==PRODUCT(dims)) .AND. (IAND(curr%flags,INPUT_DATA_IN_USE)==0)) THEN
            IF (PRESENT(ref)) ref => curr
            IF (PRESENT(no))  no  =  i
            curr%flags = curr%flags + 1 ! Warning: No overflow check!!! 
            InputDataGetRef => curr%dta
            RETURN
          ENDIF
        ENDIF
      ENDIF
      i = i + 1
      prev => curr
      curr => curr%next
    ENDDO

    CALL InputDataNew(curr)
    IF (ASSOCIATED(prev)) THEN
      prev%next => curr
    ELSE
      AllInputData => curr
    ENDIF
    IF (ALL(dims>0)) THEN
      Sz(:) = 1
      Sz(1:nDims) = dims
      ALLOCATE(curr%dta(Sz(1),Sz(2),Sz(3),Sz(4),Sz(5),Sz(6),Sz(7)))
      curr%flags = IOR(curr%flags,INPUT_DATA_ALLOC)
    ENDIF
    IF (PRESENT(ref)) ref => curr
    IF (PRESENT(no))  no  =  i
    curr%flags = IOR(curr%flags,1)
    InputDataGetRef => curr%dta

  END FUNCTION InputDataGetRef

! Finds a buffer with matching type and id
  FUNCTION InputDataGetSpecified(list,typ,id,ref)
  REAL (dp),              POINTER :: InputDataGetSpecified(:,:,:,:,:,:,:)
  TYPE (input_data_list), POINTER :: list
  INTEGER,             INTENT(in) :: typ
  INTEGER, OPTIONAL,   INTENT(in) :: id
  TYPE (input_data_list), OPTIONAL, POINTER :: ref

  TYPE (input_data_list), POINTER :: curr
  LOGICAL :: lFound

    curr => list
    DO WHILE (ASSOCIATED(curr))
      lFound = IAND(curr%flags,INPUT_DATA_TYPE_MASK) == typ
      IF (lFound .AND. PRESENT(id)) THEN
        IF ((IAND(curr%flags,INPUT_DATA_NO_MASK) /= id) .AND. (id /= -1)) lFound = .FALSE.
      ENDIF
      IF (lFound) THEN
        InputDataGetSpecified => curr%dta
        IF (PRESENT(ref)) ref => curr
        RETURN
      ELSE
        curr => curr%next
      ENDIF
    ENDDO
    NULLIFY(InputDataGetSpecified)
    IF (PRESENT(ref)) NULLIFY(ref)

  END FUNCTION InputDataGetSpecified
  
! Allocates and links buffers circularly for interpolation
  FUNCTION InputDataIntrpAlloc(nIntrp,dims)
  TYPE (input_data_list), POINTER :: InputDataIntrpAlloc
  INTEGER :: nIntrp
  INTEGER :: dims(7)

  TYPE (input_data_list), POINTER :: curr, new, head
  INTEGER :: i

    NULLIFY(curr)
    DO i=1,nIntrp
      ALLOCATE(new)
      IF (ASSOCIATED(curr)) THEN
        curr%next => new
        new%prev  => curr
      ELSE
        head => new
      ENDIF
      ALLOCATE(new%dta(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
      new%time(:) = -HUGE(i8)
      new%flags   =  0
      curr => new
    ENDDO
    head%prev => new
    new%next  => head

    InputDataIntrpAlloc => head

  END FUNCTION InputDataIntrpAlloc

! Display buffers in a list
  SUBROUTINE InputDataDisp(buf,cnt,unt,msg)
  TYPE (input_data_list), OPTIONAL, POINTER :: buf
  INTEGER, OPTIONAL, INTENT(in) :: cnt
  INTEGER, OPTIONAL, INTENT(in) :: unt
  CHARACTER (len=*), INTENT(in), OPTIONAL :: msg

  TYPE (input_data_list), POINTER :: curr
  CHARACTER (len=256) :: tmpstr
  INTEGER :: un, i, cr

    un = GetLogUnit(); IF (PRESENT(unt)) un = unt
    cr = -1;           IF (PRESENT(cnt)) cr = cnt

    IF (PRESENT(buf)) THEN
      curr => buf
    ELSE
      curr => AllInputData
    ENDIF
    i = 1
    IF (PRESENT(msg)) WRITE(un,'(a)') msg
    DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))

      WRITE(un,'(a,i4)') 'Data buffer:',i
      IF (ASSOCIATED(curr%dta)) THEN
        WRITE(tmpstr,'(7i6)') SHAPE(curr%dta)
      ELSE
        tmpstr = ' <undefined>'
      ENDIF
      WRITE(un,'(2a)') '  Data array size : ',TRIM(tmpstr)

      IF (ASSOCIATED(curr%prev)) THEN
        tmpstr = 'Yes'
      ELSE
        tmpstr = 'No'
      ENDIF
      WRITE(un,'(2a)') '  Data back link  : ',TRIM(tmpstr)

      WRITE(un,'(2a)') '  Data time       : ',DateString(curr%time)
      WRITE(message_text,*) '_#using_objects__:_',IAND(curr%flags,INPUT_DATA_NREF_MASK)
      CALL RmSpc(message_text)
      WRITE(un,*) TRIM(message_text)
      WRITE(message_text,*) '_Id_within_type__:_',IAND(curr%flags,INPUT_DATA_NO_MASK)/INPUT_DATA_NO_STEP
      CALL RmSpc(message_text)
      WRITE(un,*) TRIM(message_text)
      CALL DispFlagsSingleBit(curr%flags,INPUT_DATA_ALLOC,flg_data_n,flg_data,Txt=message_text)
      WRITE(un,'(3a,Z8.8,a)')   '  Data flags      : ',TRIM(message_text),' (0x',curr%flags,')'

      ! Next buffer
      i = i + 1
      curr => curr%next
      cr = cr - 1
    ENDDO

  END SUBROUTINE InputDataDisp

! Deletes a single buffer from a list. Buffer is only actually removed if its usage count gets 0
  SUBROUTINE InputDataUnuse(list,buf,dta)
  TYPE (input_data_list), POINTER :: list
  TYPE (input_data_list), POINTER, OPTIONAL :: buf
  REAL(dp),               POINTER, OPTIONAL :: dta(:,:,:,:,:,:,:)

  TYPE (input_data_list), POINTER :: prev, curr
  REAL(dp),               POINTER :: Tst(:,:,:,:,:,:,:)

    IF (PRESENT(dta)) THEN
      Tst => dta
    ELSE
      IF (PRESENT(buf)) THEN
        Tst => buf%dta
      ELSE
        CALL local_error('InputDataUnuse','Either buf or Dta must be specified')
      ENDIF
    ENDIF

    NULLIFY(prev)
    curr => list
    DO WHILE (ASSOCIATED(curr))
      IF (ASSOCIATED(curr%dta,Tst)) THEN
        IF (IAND(curr%flags,INPUT_DATA_NREF_MASK) /= 0) curr%flags = curr%flags - 1
        IF (IAND(curr%flags,INPUT_DATA_NREF_MASK) == 0) THEN
          IF (ASSOCIATED(prev)) THEN
            prev%next => curr%next
          ELSE
            list => curr%next
          ENDIF
          IF (IAND(curr%flags,INPUT_DATA_ALLOC) /= 0) DEALLOCATE(curr%dta)
          DEALLOCATE(curr)
        ENDIF
        NULLIFY(curr)
      ELSE
        prev => curr
        curr => curr%next
      ENDIF
    ENDDO

  END SUBROUTINE InputDataUnuse

! Deletes a list of buffers, eventually keeping the data (for shared usage)
  SUBROUTINE InputDataDone(list,lDel)
  TYPE (input_data_list), POINTER :: list
  LOGICAL, OPTIONAL,   INTENT(in) :: lDel

  TYPE (input_data_list), POINTER :: curr, prev
  LOGICAL :: Del

    Del = .TRUE.
    IF (PRESENT(lDel)) Del = lDel

    NULLIFY(prev)
    curr => list
    DO WHILE (ASSOCIATED(curr))
      IF (ASSOCIATED(curr%dta) .AND. Del) DEALLOCATE(curr%dta)
      prev => curr
      curr => curr%next
      DEALLOCATE(prev)
    ENDDO

  END SUBROUTINE InputDataDone

! Adds or replaces a specified 1-dimensional auxillary data set of a variable with a new one
  SUBROUTINE ReplaceAuxData1D(var,typ,id,dta)
  TYPE (input_var_list), POINTER :: var
  INTEGER,            INTENT(in) :: typ
  INTEGER,            INTENT(in) :: id
  REAL(dp)                       :: dta(:)

  REAL(dp),               POINTER :: aux(:,:,:,:,:,:,:)
  TYPE (input_data_list), POINTER :: aux_data

    aux => InputDataGetSpecified(var%aux_data,typ,id,aux_data)
    IF (ASSOCIATED(aux_data)) THEN
      IF (SIZE(aux,1) /= SIZE(dta)) DEALLOCATE(aux_data%dta) ! Just replace data - no reason to rebuild header
    ELSE
      CALL InputDataNew(aux_data)
      aux_data%flags = typ+id
      aux_data%next => var%aux_data
      var%aux_data  => aux_data
    ENDIF
    IF (.NOT. ASSOCIATED(aux_data%dta)) ALLOCATE(aux_data%dta(SIZE(dta,1),1,1,1,1,1,1))
    aux_data%dta(:,1,1,1,1,1,1) = dta

  END SUBROUTINE ReplaceAuxData1D

! Associate user data with variable
  INTEGER FUNCTION GetNextUserDataNumber(var)
  TYPE (input_var_list), POINTER :: var

  REAL(dp), POINTER :: tst(:,:,:,:,:,:,:)
  INTEGER :: i
  
    NULLIFY(tst)
    i = 0
    DO WHILE ((i==0) .OR. ASSOCIATED(tst))
      tst => InputDataGetSpecified(var%aux_data,INPUT_DATA_USER,i*INPUT_DATA_NO_STEP)
      IF (ASSOCIATED(tst)) THEN
        i = i + 1
      ELSE
        EXIT
      ENDIF
    ENDDO
    GetNextUserDataNumber = i

  END FUNCTION GetNextUserDataNumber

END MODULE mo_input_data
