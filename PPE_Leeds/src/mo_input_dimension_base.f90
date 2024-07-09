!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! This module contains procedures for mo_input dealing with basic list operation on dimensions 
! Never directly use this module from outside the mo_input framework - only via mo_input.f90
MODULE mo_input_dimension_base
USE mo_input_strings
USE mo_input_arrays
USE mo_input_types
USE mo_input_fileobj
IMPLICIT NONE

PUBLIC InputDimNewCtl, InputDimNew, InputDimDone, InputDimListElmRemove, InputDimCopyList
PUBLIC InputDimGetRef, InputDimGetRefEx, InputDimDisp, InputDimInq

CONTAINS

! ----------------------------------------------------
!  Creating, adding, setting, modifying, displaying:
!     Dimensions
! ----------------------------------------------------

! New input_dim with default values
  SUBROUTINE InputDimNewCtl(new)
  TYPE (input_dim_list), POINTER :: new

    ALLOCATE(new)
    NULLIFY(new%src_dims,new%next,new%dim_data)
    new%lo            = -1
    new%hi            = -1
    new%flags         =  0
    new%fid           =  0
    new%file_dim_name = ''

  END SUBROUTINE InputDimNewCtl

! New input_dim_list element with default values
  SUBROUTINE InputDimNew(new)
  TYPE (input_dim_list), POINTER :: new

    CALL InputDimNewCtl(new)
    ALLOCATE(new%dim_data)
    NULLIFY(new%dim_data%mask,new%dim_data%local)

    new%dim_data%nsrc        =  1
    new%dim_data%name_dim    = ''
    new%dim_data%name_alt    = ''
    new%dim_data%sub_model   =  1
    new%dim_data%nDst        =  0
    new%dim_data%size_global = -1
    new%dim_data%size_local  = -1
    new%dim_data%local_lo    = -1
    new%dim_data%local_hi    = -1
    new%dim_data%chunk_lo    = -1
    new%dim_data%chunk_hi    = -1
    new%dim_data%flags       =  0
    new%dim_data%origin      = UNLIKELY_RVAL

  END SUBROUTINE InputDimNew

! Destroys a list of dimensions, eventually only headers
  SUBROUTINE InputDimDone(dims,ldel,lnonmodel)
  TYPE (input_dim_list), POINTER :: dims
  LOGICAL, OPTIONAL,  INTENT(in) :: ldel, lnonmodel

  TYPE (input_dim_list), POINTER :: curr, next
  LOGICAL :: del, nonmodel

    del = .FALSE.
    IF (PRESENT(ldel)) del = ldel
    nonmodel = .FALSE.
    IF (PRESENT(lnonmodel)) nonmodel = lnonmodel
    curr => dims
    DO WHILE (ASSOCIATED(curr))
      next => curr%next
      IF (ASSOCIATED(curr%dim_data)) THEN
        IF (del .OR. (nonmodel .AND. (IAND(curr%dim_data%flags,INPUT_DIM_NON_MODEL) /= 0)) .OR. &
        (IAND(curr%flags,INPUT_DIM_LOCAL_DATA) /= 0)) THEN
          IF (ASSOCIATED(curr%dim_data%local)) DEALLOCATE(curr%dim_data%local)
          ! The mask is a pointer into the list of global masks, thus the destruction is a mask and not a dim. issue
          DEALLOCATE(curr%dim_data)
        ENDIF
      ENDIF
      DEALLOCATE(curr)
      curr => next
    ENDDO

  END SUBROUTINE InputDimDone

! Removes a single dimension header identified by its number in the list from a list
  FUNCTION InputDimListElmRemove(list,no)
  TYPE (input_dim_list), POINTER :: InputDimListElmRemove
  TYPE (input_dim_list), POINTER :: list
  INTEGER,            INTENT(in) :: no

  TYPE (input_dim_list), POINTER :: curr, prev
  INTEGER                        :: i

    NULLIFY(prev)
    curr => list
    DO i=2,no
      prev => curr
      curr => curr%next
    ENDDO

    IF (ASSOCIATED(prev)) THEN
      prev%next => curr%next
      prev => list
    ELSE
      prev => curr%next
    ENDIF
    DEALLOCATE(curr)

    InputDimListElmRemove => prev

  END FUNCTION InputDimListElmRemove

! Copies an entire list of dimensions (only headers)
  FUNCTION InputDimCopyList(dims)
  TYPE (input_dim_list), POINTER :: InputDimCopyList
  TYPE (input_dim_list), POINTER :: dims

  TYPE (input_dim_list), POINTER :: head, prev, new, curr

    NULLIFY(head,prev)
    curr => dims
    DO WHILE (ASSOCIATED(curr))

      CALL InputDimNewCtl(new)
      IF (ASSOCIATED(prev)) THEN
        prev%next => new
      ELSE
        head => new
      ENDIF
      prev => new

      new%src_dims      => curr%src_dims
      new%dim_data      => curr%dim_data
      new%lo            =  curr%lo
      new%hi            =  curr%hi
      new%flags         =  IAND(curr%flags,INPUT_ALLFLAGS - INPUT_DIM_LOCAL_DATA)
      new%fid           =  curr%fid
      new%file_dim_name =  TRIM(curr%file_dim_name)

      curr => curr%next
    ENDDO

    InputDimCopyList => head

  END FUNCTION InputDimCopyList

! Returns handle to a dimension which matches the name and sub-model
  FUNCTION InputDimGetRef(dim_name,alt_name,sub_model,search,no)
  TYPE (input_dim_list),              POINTER :: InputDimGetRef
  CHARACTER(len=*),                INTENT(in) :: dim_name
  CHARACTER(len=*),      OPTIONAL, INTENT(in) :: alt_name
  INTEGER,               OPTIONAL, INTENT(in) :: sub_model
  INTEGER,               OPTIONAL, INTENT(out):: no
  TYPE (input_dim_list), OPTIONAL,    POINTER :: search 

  TYPE (input_dim_list), POINTER :: curr
  INTEGER :: sm, i

    sm = CurrSubModel
    IF (PRESENT(sub_model)) sm = sub_model
    i = 1
    curr => AllInputDims
    IF (PRESENT(search)) curr => search
    DO WHILE (ASSOCIATED(curr))
      IF ((curr%dim_data%sub_model < 0) .OR. (sm < 0) .OR. (curr%dim_data%sub_model == sm)) THEN
        IF ((TRIM(dim_name)==TRIM(curr%dim_data%name_dim)) .OR. &
            (TRIM(dim_name)==TRIM(curr%dim_data%name_alt)) .OR. &
            (TRIM(dim_name)==TRIM(curr%file_dim_name))) THEN
          InputDimGetRef => curr
          IF (PRESENT(no)) no = i
          RETURN
        ENDIF
        IF (PRESENT(alt_name)) THEN
          IF (LEN_TRIM(alt_name) > 0) THEN
            IF ((TRIM(alt_name)==TRIM(curr%dim_data%name_dim)) .OR. &
                (TRIM(alt_name)==TRIM(curr%dim_data%name_alt)) .OR. &
                (TRIM(alt_name)==TRIM(curr%file_dim_name))) THEN
              InputDimGetRef => curr
              IF (PRESENT(no)) no = i
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      i = i + 1
      curr => curr%next
    ENDDO

    NULLIFY(InputDimGetRef)

  END FUNCTION InputDimGetRef

! Returns the relation of a dimension to dimensions in a dimension list (INPUT_DIM_FND_xxx flags):
! Through recursion, complex relations like: dimension is part of an equivalent dimension which should be packed are recognised
  RECURSIVE INTEGER FUNCTION InputDimGetRefEx(se_dim,search,ref,no,it,fdims) RESULT(ret)
  TYPE (input_dim_list),              POINTER :: se_dim 
  INTEGER,               OPTIONAL, INTENT(out):: no
  TYPE (input_dim_list), OPTIONAL,    POINTER :: search 
  TYPE (input_dim_list), OPTIONAL,    POINTER :: ref 
  INTEGER,               OPTIONAL, INTENT(in ):: it
  TYPE (input_fileobj_list), OPTIONAL,POINTER :: fdims

  TYPE (input_eqdim_list),            POINTER :: curr_eq
  TYPE (input_dim_list),              POINTER :: curr_dim, ref_dim
  TYPE (input_fileobj_list),          POINTER :: file_dim
  INTEGER                                     :: se_no, res, i, j, k, kmax, iter
  LOGICAL                                     :: lFound

    iter = 10 ! Arbitrary value > 2
    IF (PRESENT(it)) iter = it
    ret = INPUT_DIM_FND_NOT

    ! Set scatter and reverse flags
    lFound = .FALSE.
    IF (ASSOCIATED(se_dim%dim_data%local)) THEN
      IF (se_dim%dim_data%local%npts /= se_dim%dim_data%local%nvalid) lFound = .TRUE.
    ENDIF
    res = INPUT_DIM_FND_NOT
    IF ((n_pe > 1) .AND. (se_dim%dim_data%chunk_lo > 0) .AND. (se_dim%dim_data%local_lo > 0)) THEN
      IF ((se_dim%dim_data%chunk_lo /= 1)                           .OR. &
          (se_dim%dim_data%chunk_hi /= se_dim%dim_data%size_global) .OR. &
          (se_dim%dim_data%chunk_lo /= se_dim%dim_data%local_lo)    .OR. &
          (se_dim%dim_data%chunk_hi /= se_dim%dim_data%local_hi)    .OR. &
          lFound) res=INPUT_DIM_FND_SCATTER
    ENDIF
    IF (IAND(se_dim%dim_data%flags,INPUT_DIM_REVERSE) /= IAND(se_dim%flags,INPUT_DIM_REVERSE)) res = res + INPUT_DIM_FND_REVERSE

    ! Determine data swap of cyclic dimensions
    IF ((IAND(se_dim%dim_data%flags,INPUT_DIM_CYCLIC) /= 0) .AND. PRESENT(fdims) .AND. &
        (se_dim%dim_data%origin /= REAL(UNLIKELY_VAL,dp)/1000._dp)) THEN
      file_dim => InputFileObjGetRef(fdims,se_dim%dim_data%name_dim,se_dim%file_dim_name)
      IF (.NOT. ASSOCIATED(file_dim) .AND. (LEN_TRIM(se_dim%dim_data%name_alt) > 0)) &
        file_dim => InputFileObjGetRef(fdims,se_dim%dim_data%name_alt)
      IF (ASSOCIATED(file_dim)) THEN
        IF (ASSOCIATED(file_dim%points)) THEN
          IF (se_dim%dim_data%origin /= file_dim%points(1)) res = IOR(res,INPUT_DIM_FND_SWAP)
        ENDIF
      ENDIF
    ENDIF

    ! Is dimension in the list?
    ref_dim => InputDimGetRef(se_dim%dim_data%name_dim,se_dim%dim_data%name_alt,se_dim%dim_data%sub_model,search,se_no)
    IF (.NOT. ASSOCIATED(ref_dim) .AND. (LEN_TRIM(se_dim%file_dim_name)>0)) & 
      ref_dim => InputDimGetRef(se_dim%file_dim_name,sub_model=se_dim%dim_data%sub_model,search=search,no=se_no)
    IF (ASSOCIATED(ref_dim)) THEN
      IF (PRESENT(no))  no  =  se_no
      IF (PRESENT(ref)) ref => ref_dim
      ret = res + INPUT_DIM_FND_NORMAL
      RETURN
    ENDIF

    IF (iter == 2) RETURN ! Stop recursion here if requested

    ! Is dimension packed and all sources are in the list?
    IF (IAND(se_dim%dim_data%flags,INPUT_DIM_PACKED) /= 0) THEN
      curr_dim => se_dim%src_dims
      DO WHILE (ASSOCIATED(curr_dim))
        i = InputDimGetRefEx(curr_dim,search,ref_dim,se_no,it=iter+1)
        lFound = i > 0
        IF (lFound) THEN          
          IF (ASSOCIATED(curr_dim,se_dim%src_dims)) THEN ! Save first for output
            IF (PRESENT(ref)) ref => ref_dim
            IF (PRESENT(no))  no  = se_no
          ENDIF
          curr_dim => curr_dim%next
        ELSE
          NULLIFY(curr_dim) ! At least one source dimension not found -> this is not an option
        ENDIF
      ENDDO
      IF (lFound) THEN
        ! The dimension may be distributed by doesn't here matter since the dimension is not directly used
        ret = res + INPUT_DIM_FND_SOURCE
        IF (IAND(curr_dim%dim_data%flags,INPUT_DIM_DISTRIB) /= 0) ret = IOR(ret,INPUT_DIM_FND_DISTRIB)
        RETURN
      ENDIF
      res = IAND(res,INPUT_ALLFLAGS - INPUT_DIM_FND_DISTRIB)
    ENDIF

    ! Is dimension a part of the source chain of a packed dimension in the searched list?
    IF (PRESENT(search)) THEN
      curr_dim => search
    ELSE
      curr_dim => AllInputDims
    ENDIF
    DO WHILE (ASSOCIATED(curr_dim))
      IF (IAND(curr_dim%dim_data%flags,INPUT_DIM_PACKED) /= 0) THEN
        i = InputDimGetRefEx(se_dim,curr_dim%src_dims,ref_dim,se_no,it=iter+1)
        IF (ASSOCIATED(ref_dim)) THEN
          IF (PRESENT(ref)) ref => curr_dim
          IF (PRESENT(no))  no  =  se_no
          se_dim%lo = ref_dim%lo
          se_dim%hi = ref_dim%hi
          ! The dimension may be distributed by doesn't here matter since the dimension is not directly used
          ret = res + INPUT_DIM_FND_PACKED
          IF (IAND(curr_dim%dim_data%flags,INPUT_DIM_DISTRIB) /= 0) ret = IOR(ret,INPUT_DIM_FND_DISTRIB)
          RETURN
        ENDIF
      ENDIF
      curr_dim => curr_dim%next
    ENDDO

    ! Is dimension a part of the source chain of equivalent dimensions in the searched list?
    kmax = 2
    IF (ColDimAssociated) kmax=1
    curr_eq => AllInputEqDims
    i = 1
    DO k=1,kmax
      DO WHILE (ASSOCIATED(curr_eq))
        ! Test if dimension is part of the source of an equivalent declaration
        ref_dim => InputDimGetRef(se_dim%dim_data%name_dim,se_dim%dim_data%name_alt,se_dim%dim_data%sub_model, &
                                  curr_eq%src_dims,se_no)
        IF (ASSOCIATED(ref_dim)) THEN
          lFound = .TRUE.
          IF (PRESENT(search)) THEN
            curr_dim => curr_eq%dst_dims
            DO WHILE (ASSOCIATED(curr_dim))
              j = InputDimGetRefEx(curr_dim,search,it=iter+1,ref=ref_dim)
              IF (j <= 0) THEN
                lFound = .FALSE.
                EXIT
              ENDIF
              curr_dim => curr_dim%next
            ENDDO
          ENDIF
          IF (lFound) THEN
            IF (PRESENT(ref)) NULLIFY(ref)
            IF (PRESENT(no))  no = i
            ret = res + INPUT_DIM_FND_EQ_SRC
            IF (IAND(j,INPUT_DIM_FND_DISTRIB+INPUT_DIM_FND_SCATTER) /= 0) ret = IOR(ret,INPUT_DIM_FND_DISTRIB)
            RETURN
          ENDIF
        ENDIF
        i = i + 1
        curr_eq => curr_eq%next
      ENDDO
      curr_eq => AllInputColDimDistrib
    ENDDO

    ! None of it!
    IF (PRESENT(ref)) NULLIFY(ref)
    IF (PRESENT(no))  no = -1
    Ret = INPUT_DIM_FND_NOT

  END FUNCTION InputDimGetRefEx

! Display dimensions from a list
  SUBROUTINE InputDimDisp(dims,cnt,unt,msg,list)
  TYPE (input_dim_list), OPTIONAL, POINTER :: dims
  INTEGER, OPTIONAL, INTENT(in) :: cnt
  INTEGER, OPTIONAL, INTENT(in) :: unt
  CHARACTER (len=*), INTENT(in), OPTIONAL :: msg
  LOGICAL, OPTIONAL, INTENT(in) :: list

  TYPE (input_dim_list),  POINTER :: curr, curr_src
  TYPE (input_mask_list), POINTER :: mask
  LOGICAL, POINTER :: lmask(:)
  CHARACTER (LEN=256) :: tmpstr
  LOGICAL :: lCont
  INTEGER :: un, i, j, k, mn, mx, cr

    un = GetLogUnit(); IF (PRESENT(unt)) un = unt
    cr = -1;           IF (PRESENT(cnt)) cr = cnt

    IF (PRESENT(dims)) THEN
      curr => dims
    ELSE
      curr => AllInputDims
    ENDIF
    IF (PRESENT(list)) THEN
      IF (list) THEN
        tmpstr = ''
        DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
          tmpstr = TRIM(tmpstr)//','//TRIM(curr%dim_data%name_dim)
          curr => curr%next
          cr = cr - 1
        ENDDO
        IF (PRESENT(msg)) THEN
          WRITE(un,'(3a)') msg,TRIM(tmpstr(2:))
        ELSE
          WRITE(un,'(2a)') TRIM(tmpstr(2:))
        ENDIF
        RETURN
      ENDIF
    ENDIF
    IF (PRESENT(msg)) WRITE(un,'(a)') msg
    DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
      i = 0
      IF (curr%dim_data%name_alt /= '') i = i + 1
      IF (curr%file_dim_name     /= '') i = i + 2
      SELECT CASE (i)
        CASE (0)
          tmpstr = TRIM(curr%dim_data%name_dim)
        CASE (1)
          tmpstr = TRIM(curr%dim_data%name_dim)//'%or%'//TRIM(curr%dim_data%name_alt)//'%(A)'
        CASE (2)
          tmpstr = TRIM(curr%dim_data%name_dim)//'%or%'//TRIM(curr%file_dim_name)//'%(F)'
        CASE (3)
          tmpstr = TRIM(curr%dim_data%name_dim)//',%'//TRIM(curr%dim_data%name_alt)//'%(A)%or%'//TRIM(curr%file_dim_name)//'%(F)'
      END SELECT
      WRITE(tmpstr,*) TRIM(tmpstr),'%(',curr%dim_data%sub_model,')'
      CALL RmSpc(tmpstr,'%')
      WRITE(un,'(a)') TRIM(tmpstr)
      WRITE(tmpstr,*) curr%dim_data%size_global
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)') '  Global size       :',TRIM(tmpstr)
      IF (ASSOCIATED(curr%dim_data%local)) THEN
        WRITE(tmpstr,*) curr%dim_data%size_local,'_(Mask)'
      ELSE
        WRITE(tmpstr,*) curr%dim_data%size_local,'_(',curr%dim_data%local_lo,'-',curr%dim_data%local_hi,')'
      ENDIF
      CALL RmSpc(tmpstr)
      WRITE(un,'(2a)') '  Local  size       : ',TRIM(tmpstr)
      IF (liodata) THEN
        IF ((curr%lo == -1) .AND. (curr%hi == -1)) THEN
          WRITE(tmpstr,*) curr%dim_data%chunk_lo,'-',curr%dim_data%chunk_hi
        ELSE
          WRITE(tmpstr,*) curr%lo,'-',curr%hi,'_(sub-dim)'
        ENDIF
        CALL RmSpc(tmpstr)
        WRITE(un,'(2a)') '  Read chunk        : ',TRIM(tmpstr)
        IF (lioserver .AND. ASSOCIATED(curr%dim_data%local)) THEN
          ALLOCATE(lmask(curr%dim_data%chunk_hi - curr%dim_data%chunk_lo + 1))
          j = 1
          tmpstr = ''
          mask => curr%dim_data%local
          DO WHILE (ASSOCIATED(mask))
            lCont = j < n_pe
            DO WHILE(lCont)
              IF (src_pe(j) /= my_pe) THEN
                j = j + 1
              ELSE
                lCont = .FALSE.
              ENDIF
              IF (lCont) lCont = j < n_pe
            ENDDO
            CALL MaskInt2Logical(lmask,mask%mask)
            mn = -1
            mx = -1
            DO k=1,SIZE(lmask)
              IF (lmask(k)) THEN
                IF (mn == -1) mn = k
                mx = k
              ENDIF
            ENDDO
            WRITE(tmpstr,*) TRIM(tmpstr),j-1,'_(',COUNT(lmask),'_pts_in_range:_',mn,'-',mx,'),_'
            IF (LEN_TRIM(tmpstr) > MAX_LINE_LEN) THEN
              CALL RmSpc(tmpstr)
              WRITE(un,'(2a)') '  Distribute to     : ',tmpstr(1:LEN_TRIM(tmpstr)-1)
              tmpstr = ''
            ENDIF
            mask => mask%next
            j = j + 1
          ENDDO
          CALL RmSpc(tmpstr)
          IF (TRIM(tmpstr) /= '') WRITE(un,'(2a)') '  Distribute to     : ',tmpstr(1:LEN_TRIM(tmpstr)-1)
        ENDIF
      ELSE
        WRITE(tmpstr,*) src_pe(my_pe+1)
        CALL RmDblSpc(tmpstr)
        WRITE(un,'(2a)') '  Recieve from PE   :',TRIM(tmpstr)
      ENDIF
      i = 0
      mask => curr%dim_data%mask
      DO WHILE (ASSOCIATED(mask))
        i = i + 1
        mask => mask%next
      ENDDO
      j = 0
      mask => curr%dim_data%local
      DO WHILE (ASSOCIATED(mask))
        j = j + 1
        mask => mask%next
      ENDDO
      WRITE(tmpstr,*) i,'/',j
      CALL RmSpc(tmpstr)
      WRITE(un,'(2a)') '  Alloc. mask/local : ',TRIM(tmpstr)
      IF (ASSOCIATED(curr%src_dims)) THEN
        tmpstr = ''
        curr_src => curr%src_dims
        DO WHILE (ASSOCIATED(curr_src)) 
          tmpstr = TRIM(tmpstr)//'_'//TRIM(curr_src%dim_data%name_dim)//','
          curr_src => curr_src%next
        ENDDO
        WRITE(message_text,*) tmpstr(2:LEN_TRIM(tmpstr)-1),'_(',curr%dim_data%nsrc,'_points)'
        CALL RmSpc(message_text)
        WRITE(un,'(2a)') '  Source dimensions : ',TRIM(message_text)
      ENDIF
      i = IOR(curr%flags,curr%dim_data%flags)
      CALL DispFlagsSingleBit(i,INPUT_DIM_UNLIM,flg_dim_n,flg_dim,Txt=tmpstr)
      WRITE(un,'(3a,Z8.8,a)') '  Flags             : ',TRIM(tmpstr),' (0x',i,')'
      IF (IAND(curr%dim_data%flags,INPUT_DIM_CYCLIC) /= 0) &
        WRITE(un,'(a,g10.4)')   '  Base value(origin):',curr%dim_data%origin
      WRITE(tmpstr,*) curr%fid
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)') '  File id if avail  :',TRIM(tmpstr)

      curr => curr%next
      cr = cr - 1
    ENDDO

  END SUBROUTINE InputDimDisp

! Inquires information about a dimension
  SUBROUTINE InputDimInq(name_dim,ref,fid,alt_name,global,local,local_lo,local_hi,chunk_lo,chunk_hi)
  CHARACTER (len=*),     OPTIONAL, INTENT(in)  :: name_dim
  TYPE (input_dim_list), OPTIONAL, POINTER     :: ref
  INTEGER,               OPTIONAL, INTENT(out) :: fid
  CHARACTER (len=64),    OPTIONAL, INTENT(out) :: alt_name
  INTEGER,               OPTIONAL, INTENT(out) :: global
  INTEGER,               OPTIONAL, INTENT(out) :: local, local_lo, local_hi
  INTEGER,               OPTIONAL, INTENT(out) :: chunk_lo, chunk_hi
 
  TYPE (input_dim_list), POINTER :: curr

    IF (PRESENT(name_dim)) THEN
      curr => InputDimGetRef(name_dim)
      IF (PRESENT(ref)) THEN
        IF (ASSOCIATED(ref) .AND. .NOT. ASSOCIATED(ref,curr)) &
          CALL local_error('InputDimInq','Name and reference refere to different dimensions')
        ref => curr
      ENDIF
    ELSE
      IF (.NOT. PRESENT(ref))    CALL local_error('InputDimInq','Neither name nor reference specified')
      IF (.NOT. ASSOCIATED(ref)) CALL local_error('InputDimInq','Neither name nor reference specified')
      curr => ref
    ENDIF

    IF (PRESENT(alt_name)) THEN
      IF (LEN_TRIM(curr%file_dim_name) > 0) THEN
        alt_name = TRIM(curr%file_dim_name)
      ELSE
        alt_name = TRIM(curr%dim_data%name_alt)
      ENDIF
    ENDIF

    IF (PRESENT(fid))      fid      = curr%fid
    IF (PRESENT(global))   global   = curr%dim_data%size_global
    IF (PRESENT(local))    local    = curr%dim_data%size_local
    IF (PRESENT(local_lo)) local_lo = curr%dim_data%local_lo
    IF (PRESENT(local_hi)) local_hi = curr%dim_data%local_hi
    IF (PRESENT(chunk_lo)) chunk_lo = curr%dim_data%chunk_lo
    IF (PRESENT(chunk_hi)) chunk_hi = curr%dim_data%chunk_hi

  END SUBROUTINE InputDimInq
  
END MODULE mo_input_dimension_base
