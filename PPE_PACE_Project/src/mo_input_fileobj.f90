!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! ----------------------------------------------------
!  Creating, adding, setting, modifying, displaying:
!    File objects 
! (holding either file dimensions or file variables)
! ----------------------------------------------------
MODULE mo_input_fileobj
USE mo_input_strings
USE mo_input_types
IMPLICIT NONE

PUBLIC InputFileObjNew, InputFileObjsDone, InputFileObjGetRef, InputFileObjExpandVar, InputFileObjsCmp, InputFileObjDisp

CONTAINS

! New file_obj element
  SUBROUTINE InputFileObjNew(Obj)
  TYPE (input_fileobj_list), POINTER :: Obj

  TYPE (input_fileobj_list), POINTER :: new

    ALLOCATE(new)
    NULLIFY(new%next,new%points)
    new%name_obj    = ''
    new%id          = -1
    new%size_global =  0
    new%flags       =  0
    new%fids_dim(:) =  0

    Obj => new

  END SUBROUTINE InputFileObjNew

! Deallocate a chain of file objects 
  SUBROUTINE InputFileObjsDone(Objs)
  TYPE (input_fileobj_list), POINTER :: Objs

  TYPE (input_fileobj_list), POINTER :: curr, next 

    curr => Objs
    DO WHILE (ASSOCIATED(curr))
      next => curr%next
      IF (ASSOCIATED(curr%points)) DEALLOCATE(curr%points)
      DEALLOCATE(curr)
      curr => next
    END DO

  END SUBROUTINE InputFileObjsDone

! Returns reference to the file object with the appropriate name
  FUNCTION InputFileObjGetRef(objs,name_obj,name_alt)
  TYPE (input_fileobj_list),      POINTER :: InputFileObjGetRef
  TYPE (input_fileobj_list),      POINTER :: objs
  CHARACTER (len=*),           INTENT(in) :: name_obj
  CHARACTER (len=*), OPTIONAL, INTENT(in) :: name_alt

  TYPE (input_fileobj_list), POINTER :: curr

    curr => objs
    DO WHILE (ASSOCIATED(curr))
      IF (TRIM(name_obj) == TRIM(curr%name_obj)) THEN
        InputFileObjGetRef => curr
        RETURN
      ENDIF
      IF (PRESENT(name_alt)) THEN
        IF (TRIM(name_alt) == TRIM(curr%name_obj)) THEN
          InputFileObjGetRef => curr
          RETURN
        ENDIF
      ENDIF
      curr => curr%next
    ENDDO
    NULLIFY(InputFileObjGetRef)

  END FUNCTION InputFileObjGetRef

! Replaces a variable name containing a wildcard with the list of matching variable names found in the file
! No checks are done for string length overflows
  SUBROUTINE InputFileObjExpandVar(fvars,name_var,Pos,lExtra)
  TYPE (input_fileobj_list), POINTER :: fvars
  CHARACTER(len=128), INTENT(inout)  :: name_var
  INTEGER,            INTENT(in)     :: Pos      ! Position of the wildcard to replace
  LOGICAL,            INTENT(in)     :: lExtra   ! Are extra underscores allowed in the expansion?

  TYPE (input_fileobj_list), POINTER :: curr
  CHARACTER(len=128) :: Res, Search
  INTEGER :: Before, After, i

    ! Split the name string and prepare the exact string to search for
    Before = INDEX(name_var(1:Pos),',',BACK=.TRUE.)
    After  = Pos+1
    Search = name_var(Before+1:Pos-1)
    IF (Search(1:1) == '+') THEN
      i = INDEX(name_var,',')
      IF (i==0) CALL local_error('InputFileObjExpandVar','Illegal use of "+" in variable name')
      i = INDEX(name_var(1:i),'_',BACK=.TRUE.)
      Search = name_var(1:i) // TRIM(Search)
    ENDIF
    i = LEN_TRIM(Search)

    ! Make list of matching names
    Res    = ''
    curr => fvars
    DO WHILE (ASSOCIATED(curr))
      IF (curr%name_obj(1:i) == TRIM(Search)) THEN
        IF (lExtra .OR. (INDEX(curr%name_obj(i+1:),'_') == 0)) THEN
          Res = TRIM(Res)//','//TRIM(curr%name_obj)
        ENDIF
      ENDIF
      curr => curr%next
    ENDDO

    ! Join the resultant string with the parts not changes
    IF (LEN_TRIM(Res) > 0) THEN
      name_var = name_var(1:Before) // TRIM(Res(2:)) // TRIM(name_var(After:))
    ELSE
      IF (Before > 0) THEN
        name_var = name_var(1:Before-1) // TRIM(name_var(After:)) ! Remove additional comma
      ELSE
        IF (name_var(After:After)==',') THEN
          name_var = name_var(2:)
        ELSE
          name_var = ''
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE InputFileObjExpandVar

! Compare two chains of file objects
  LOGICAL FUNCTION InputFileObjsCmp(O1, O2)
  TYPE (input_fileobj_list), POINTER :: O1, O2

  TYPE (input_fileobj_list), POINTER :: c1, c2

    InputFileObjsCmp = .FALSE.
    IF (.NOT. (ASSOCIATED(O1) .AND. ASSOCIATED(O2))) RETURN ! Comparison with at least one empty list
    c1 => O1
    c2 => O2
    DO WHILE (ASSOCIATED(c1) .AND. ASSOCIATED(c2))
      IF (TRIM(c1%name_obj) /= TRIM(c2%name_obj)) RETURN
      IF (c1%flags       /= c2%flags            ) RETURN
      IF (c1%id          /= c2%id               ) RETURN
      IF ((c1%size_global /= c2%size_global) .AND. (IAND(c1%flags,INPUT_FILEDIM_DIR_UNLIM)==0)) RETURN
      IF (ASSOCIATED(c1%points)) THEN
        IF (.NOT. ASSOCIATED(c2%points))        RETURN
        IF (SIZE(c1%points) /= SIZE(c2%points)) RETURN
        IF (ANY(c1%points-c2%points /= 0))      RETURN
      ENDIF
      IF (ANY(c1%fids_dim-c2%fids_dim /= 0))      RETURN
      c1 => c1%next
      c2 => c2%next
    ENDDO

    InputFileObjsCmp = .NOT. (ASSOCIATED(c1) .OR. ASSOCIATED(c2)) ! Both lists must have been processed to the end

  END FUNCTION InputFileObjsCmp

! Display buffers in a list
  SUBROUTINE InputFileObjDisp(obj,cnt,unt,msg,list)
  TYPE (input_fileobj_list), POINTER :: obj
  INTEGER, OPTIONAL,      INTENT(in) :: cnt
  INTEGER, OPTIONAL,      INTENT(in) :: unt
  CHARACTER (len=*), INTENT(in), OPTIONAL :: msg
  LOGICAL, OPTIONAL,      INTENT(in) :: list

  TYPE (input_fileobj_list), POINTER :: curr
  CHARACTER (len=128) :: tmpstr
  INTEGER :: un, cr

    un = GetLogUnit(); IF (PRESENT(unt)) un = unt
    cr = -1;           IF (PRESENT(cnt)) cr = cnt

    curr => obj
    IF (PRESENT(list)) THEN
      IF (list) THEN
        message_text = ''
        DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
          message_text = TRIM(message_text)//','//TRIM(curr%name_obj)
          curr => curr%next
          cr = cr - 1
        ENDDO
        IF (PRESENT(msg)) THEN
          WRITE(un,'(2a)') msg,TRIM(message_text(2:))
        ELSE
          WRITE(un,'(a)') TRIM(message_text(2:))
        ENDIF
        RETURN
      ENDIF
    ENDIF
    IF (PRESENT(msg)) WRITE(un,'(a)') msg
    DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))

      WRITE(un,'(2a)')   '  ',TRIM(curr%name_obj)
      WRITE(message_text,*) curr%id
      CALL RmSpc(message_text)
      WRITE(un,'(2a)')   '    File id    : ',TRIM(message_text)
      IF (IAND(obj%flags,INPUT_FILEOBJ_VAR)/=0) THEN
        CALL DispFlagsSingleBit(curr%flags,INPUT_FILEVAR_HAS_SCALE,flg_filevar_n,flg_fileobj_var,Txt=tmpstr)
        WRITE(message_text,*) '  Dimensions : ',curr%size_global
        WRITE(message_text,*) TRIM(message_text),' (fids:',curr%fids_dim(1:curr%size_global),',',TRIM(tmpstr),')'
        CALL RmDblSpc(message_text(13:))
        WRITE(un,'(a)') TRIM(message_text)
      ELSEIF (IAND(obj%flags,INPUT_FILEOBJ_DIM)/=0) THEN
        WRITE(message_text,*) '   Size       : ',curr%size_global,'(', &
        TRIM(GetStringIndexed(flg_fileobj_dir,IAND(curr%flags,INPUT_FILEDIM_DIR_MASK)/INPUT_FILEDIM_DIR_ERROR+1,',')),')'
        CALL RmDblSpc(message_text(16:))
        WRITE(un,'(a)') TRIM(message_text)
        IF (ASSOCIATED(curr%points)) THEN
          WRITE(message_text,*) '   Data range : ',MINVAL(curr%points),' - ',MAXVAL(curr%points)
          CALL RmDblSpc(message_text(5:))
          WRITE(un,'(a)') TRIM(message_text)
        ENDIF
      ELSE
        WRITE(un,'(a)') '  File object is of undetermined type'
      ENDIF

      ! Next obj 
      curr => curr%next
      cr = cr - 1
    ENDDO

  END SUBROUTINE InputFileObjDisp

END MODULE mo_input_fileobj

