!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! This module contains procedures for mo_input dealing with dimensions and equivalent dimensions 
! Never directly use this module from outside the mo_input framework - only via mo_input.f90
MODULE mo_input_dimension
USE mo_input_strings
USE mo_input_arrays
USE mo_input_types
USE mo_input_base
USE mo_input_dimension_base
USE mo_input_mask
IMPLICIT NONE

INTERFACE InputDimPackedAdd
  MODULE PROCEDURE InputDimPackedAddnd
  MODULE PROCEDURE InputDimPackedAdd1d
  MODULE PROCEDURE InputDimPackedAdd2d
  MODULE PROCEDURE InputDimPackedAdd3d
  MODULE PROCEDURE InputDimPackedAdd4d
  MODULE PROCEDURE InputDimPackedAdd5d
  MODULE PROCEDURE InputDimPackedAdd6d
  MODULE PROCEDURE InputDimPackedAdd7d
END INTERFACE
INTERFACE InputDimLocalSetMulti
  MODULE PROCEDURE InputDimLocalSetMultiMask
  MODULE PROCEDURE InputDimLocalSetMultiInt
  MODULE PROCEDURE InputDimLocalSetMulti1d
  MODULE PROCEDURE InputDimLocalSetMulti2d
  MODULE PROCEDURE InputDimLocalSetMulti3d
  MODULE PROCEDURE InputDimLocalSetMulti4d
  MODULE PROCEDURE InputDimLocalSetMulti5d
  MODULE PROCEDURE InputDimLocalSetMulti6d
  MODULE PROCEDURE InputDimLocalSetMulti7d
END INTERFACE

PUBLIC :: InputDimCopy, InputDimListFromString, InputDimAdd, InputDimPackedAdd, InputDimPackedAddMask, InputDimLocalSet
PUBLIC :: InputDimLocalSetMulti, InputDimSync, InputDimEquivalence, InputEqDimDisp

CONTAINS

! New input_dim_list element which is a copy of an existing one and mark it as a local copy
  SUBROUTINE InputDimCopy(src,new,add,ldata)
  TYPE (input_dim_list), POINTER :: src, new
  CHARACTER (len=*), OPTIONAL, INTENT(in) :: add
  LOGICAL,           OPTIONAL, INTENT(in) :: ldata

  CHARACTER (len=32) :: aadd
  LOGICAL :: copydata

    aadd = ''
    IF (PRESENT(add)) aadd = TRIM(add)
    copydata = .TRUE.
    IF (PRESENT(ldata)) copydata = ldata
    IF (copydata) THEN
      CALL InputDimNew(new)

      new%dim_data%name_dim    = TRIM(src%dim_data%name_dim)//TRIM(aadd)
      new%dim_data%name_alt    = TRIM(src%dim_data%name_alt)
      new%dim_data%sub_model   = src%dim_data%sub_model
      new%dim_data%flags       = src%dim_data%flags + INPUT_DIM_LOCAL_DATA
      new%dim_data%size_global = src%dim_data%size_global
      new%dim_data%size_local  = src%dim_data%size_local
      new%dim_data%chunk_lo    = src%dim_data%chunk_lo
      new%dim_data%chunk_hi    = src%dim_data%chunk_hi
      new%dim_data%local_lo    = src%dim_data%local_lo
      new%dim_data%local_hi    = src%dim_data%local_hi
      new%dim_data%nDst        = src%dim_data%nDst
      new%dim_data%nSrc        = src%dim_data%nSrc
      new%dim_data%origin      = src%dim_data%origin
      IF (ASSOCIATED(src%dim_data%local)) THEN
        new%dim_data%local => src%dim_data%local ! Is it sufficient to reference?
      ENDIF
      IF (ASSOCIATED(src%dim_data%local)) CALL InputMaskCopy(new%dim_data%local,src%dim_data%local)
    ELSE
      CALL InputDimNewCtl(new)
      new%dim_data => src%dim_data
    ENDIF

    new%file_dim_name =  TRIM(src%file_dim_name)
    new%lo            =  src%lo
    new%hi            =  src%hi
    new%flags         =  src%flags
    new%src_dims      => src%src_dims
    new%fid           =  src%fid

  END SUBROUTINE InputDimCopy

! Build list of dimensions from a comma separated string of names of existing global dimension
  FUNCTION InputDimListFromString(dims,flag,txt,sub_model,sz,add)
  TYPE (input_dim_list), POINTER :: InputDimListFromString
  CHARACTER (len=*), INTENT(in ) :: dims
  INTEGER,           INTENT(in ) :: flag
  CHARACTER (len=*), INTENT(in ) :: txt
  INTEGER, OPTIONAL, INTENT(in ) :: sub_model
  INTEGER, OPTIONAL, INTENT(out) :: sz
  CHARACTER (len=*), INTENT(in ), OPTIONAL :: add

  TYPE (input_eqdim_list), POINTER :: tst
  TYPE (input_dim_list),   POINTER :: head, prev, curr, search
  INTEGER cSz, st, en, le, cnt, nfnd
  LOGICAL lnext

    NULLIFY(head,curr, prev)
    cnt = 0
    cSz = 1
    le  = LEN_TRIM(dims)
    st  = 1
    DO WHILE (st < le)
      en = INDEX(dims(st:),',') + st - 1
      IF (en < st) en = le + 1
      NULLIFY(search)
      IF (PRESENT(add)) THEN
        IF (LEN_TRIM(add) > 0) search => InputDimGetRef(dims(st:en-1)//TRIM(add),sub_model=sub_model)
        IF (ASSOCIATED(search)) cnt = cnt + 1
      ENDIF
      IF (.NOT. ASSOCIATED(search)) THEN
        search => InputDimGetRef(dims(st:en-1),sub_model=sub_model)
        IF (.NOT. ASSOCIATED(search)) CALL local_error(Txt,'Reference to unknown dimension: '//dims(st:en-1))
      ENDIF
      CALL InputDimCopy(search,curr,ldata=.FALSE.)
      curr%flags = IAND(curr%flags,INPUT_ALLFLAGS-INPUT_DIM_TYPE) + Flag
      IF (search%dim_data%size_local > 0) THEN
        IF (cSz > 0) cSz = cSz * search%dim_data%size_local
      ELSE
        IF (search%dim_data%size_global < 0) cSz = -1
        IF (cSz > 0) cSz = cSz * search%dim_data%size_global
      ENDIF
      IF (ASSOCIATED(prev)) THEN
        prev%next => curr
      ELSE
        head => curr
      ENDIF
      prev => curr
      st = en + 1
    ENDDO
    IF (cnt > 0) THEN ! At least one dimension is found with the extension "_iloc", indicating a part of collectively distrib dim.
      ! Check if the present _iloc dimensions exactly matches the dimensions of a collective distribution
      nfnd = -1
      tst => AllInputColDimDistrib
      DO WHILE (ASSOCIATED(tst))
        nfnd =  0
        curr => head
        lnext =  ASSOCIATED(curr)
        DO WHILE (lnext)
          le = LEN_TRIM(curr%dim_data%name_dim)
          lnext = le < 5
          IF (.NOT. lnext) lnext = curr%dim_data%name_dim(le-4:le)/='_iloc'
          IF (.NOT. lnext) THEN
            prev => InputDimGetRef(TRIM(curr%dim_data%name_dim),search=tst%dst_dims,sub_model=sub_model)
            IF (ASSOCIATED(prev)) nfnd = nfnd + 1
          ENDIF
          curr => curr%next
          lnext = ASSOCIATED(curr)
        ENDDO
        IF (nfnd==cnt) THEN
          NULLIFY(tst)
        ELSE
          tst => tst%next
        ENDIF
      ENDDO
      IF (nfnd /= cnt) CALL local_error('InputDimListFromString', & 
        'Dimension list seems to involve a collective dimension distribution, but no matching is found')
    ENDIF

    IF (PRESENT(sz)) sz = cSz
    InputDimListFromString => head

  END FUNCTION InputDimListFromString

! Add a dimension
  SUBROUTINE InputDimAdd(dim_name,lRecSep,lCyclic,Extent,alt_name,ldec,dim_ref,origin,sub_model, &
                         local_lo,local_hi,llocal,ilocal,intr)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  LOGICAL, OPTIONAL,   INTENT(IN) :: lRecSep, lCyclic
  INTEGER, OPTIONAL,   INTENT(IN) :: Extent, local_lo, local_hi, sub_model
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec, intr
  REAL(dp),OPTIONAL,   INTENT(IN) :: origin
  LOGICAL, OPTIONAL,   INTENT(IN), POINTER :: llocal(:)
  INTEGER, OPTIONAL,   INTENT(IN), POINTER :: ilocal(:)
  CHARACTER (LEN=*),      OPTIONAL, INTENT(IN) :: alt_name
  TYPE (input_dim_list),  OPTIONAL,    POINTER :: dim_ref

  TYPE (input_dim_list),  POINTER :: new, prev
  TYPE (input_file_list), POINTER :: curr_file
  TYPE (input_mask_list), POINTER :: local_mask
  INTEGER, POINTER :: mask(:)
  LOGICAL :: ldummy(1)

    IF (.NOT. moInputInitialized) CALL local_error('InputDimAdd','Please call InputInit before adding dimensions')
    IF (.NOT. PRESENT(intr) .AND. lInLoop) CALL local_error('InputDimAdd','Cannot add dimensions after entering the main loop')

    ! Test if dimension already exist
    new => InputDimGetRef(dim_name,alt_name,sub_model)
    IF (ASSOCIATED(new)) THEN
      IF (PRESENT(dim_ref)) dim_ref => new
      RETURN
    ENDIF

    ! Create new dimension and add it to the end of the global list
    CALL InputDimNew(new)
    prev => AllInputDims
    IF (ASSOCIATED(prev)) THEN
      DO WHILE (ASSOCIATED(prev%next))
        prev => prev%next
      ENDDO
      prev%next => new
    ELSE
      AllInputDims => new
    ENDIF

    new%dim_data%name_dim = dim_name
    IF (PRESENT(lRecSep)) THEN 
      IF (lRecSep) new%dim_data%flags = IOR(new%dim_data%flags,INPUT_DIM_UNLIM)
    ENDIF

    IF (PRESENT(alt_name)) new%dim_data%name_alt = alt_name
    IF (PRESENT(lCyclic)) THEN
      IF (lCyclic) THEN
        new%dim_data%flags = IOR(new%dim_data%flags,INPUT_DIM_CYCLIC)
        IF (PRESENT(origin)) new%dim_data%origin = origin
      ENDIF
    ENDIF

    ! Get dimension size
    new%dim_data%size_global = -1
    IF (PRESENT(Extent)) THEN
      new%dim_data%size_global = Extent
      new%dim_data%size_local  = Extent
      new%dim_data%local_lo    = 1
      new%dim_data%local_hi    = Extent
      new%dim_data%chunk_lo    = 1
      new%dim_data%chunk_hi    = Extent
    ENDIF
    NULLIFY(local_mask)
    IF (PRESENT(llocal)) THEN
      IF (SIZE(llocal)/=new%dim_data%size_global) CALL local_error('InputDimAdd','Size of logical mask must match dimension size')
      mask => InputMaskCreate(llocal,SIZE(llocal),Ref=local_mask,Flags=INPUT_MASK_TYPE_PE+INPUT_MASK_GLOBAL)
      mask(1) = mask(1) ! To avoid compiler warnings
    ELSE
      IF (PRESENT(ilocal)) THEN
        ldummy(:) = .FALSE.
        mask => InputMaskCreate(ldummy,1,Ref=local_mask,mask=ilocal,Flags=INPUT_MASK_TYPE_PE+INPUT_MASK_GLOBAL)
        mask(1) = mask(1) ! To avoid compiler warnings
      ENDIF
    ENDIF
    IF ((PRESENT(local_lo) .AND. PRESENT(local_hi)) .OR. ASSOCIATED(local_mask)) &
      CALL InputDimLocalSet(local_lo,local_hi,local_mask,dim_ref=new)
    IF (PRESENT(ldec)) THEN
      IF (ldec) THEN
        new%dim_data%flags = IOR(new%dim_data%flags,INPUT_DIM_DEC)
        new%flags          = IOR(new%flags,INPUT_DIM_DEC)
      ELSE
        new%dim_data%flags = IAND(new%dim_data%flags,INPUT_ALLFLAGS - INPUT_DIM_DEC)
        new%flags          = IAND(new%flags,INPUT_ALLFLAGS - INPUT_DIM_DEC)
      ENDIF
    ENDIF
    new%flags = IOR(new%flags,INPUT_DIM_GLOBAL)

    IF (PRESENT(dim_ref)) dim_ref => new

    new%dim_data%sub_model = CurrSubModel
    IF (PRESENT(sub_model)) new%dim_data%sub_model = sub_model

    ! Force all input files to recheck their dimensions
    curr_file => AllInputFiles
    DO WHILE (ASSOCIATED(curr_file))
      curr_file%stat = IAND(curr_file%stat,INPUT_ALLFLAGS - INPUT_FILE_CHECK_DIM_SIZE - INPUT_FILE_CHECK_VAR_DIMS)
      curr_file => curr_file%next
    ENDDO

  END SUBROUTINE InputDimAdd

! Adds a mask to a packed dimension
  SUBROUTINE InputDimPackedAddMask(pdim,mask,lmask,lmask_sz,imask,Extent)
  TYPE (input_dim_list), POINTER  :: pdim
  LOGICAL, OPTIONAL               :: lmask(*)
  INTEGER, OPTIONAL,  INTENT(IN)  :: lmask_sz
  TYPE (input_mask_list), POINTER, OPTIONAL :: mask
  INTEGER, OPTIONAL,     POINTER  :: imask(:)
  INTEGER, OPTIONAL,  INTENT(IN)  :: Extent

  TYPE (input_mask_list), POINTER :: new
  INTEGER, POINTER :: NewMask(:)
  LOGICAL :: lDummy(1)

    IF (.NOT. ((PRESENT(lmask) .AND. PRESENT(lmask_sz)) .OR. PRESENT(imask) .OR. PRESENT(mask))) &
      CALL local_error('InputDimPackedAddMask','Either mask or imask must be present')
    IF (PRESENT(Extent)) THEN
      pdim%dim_data%size_global = Extent
    ELSE
      IF (PRESENT(lmask)) THEN 
        pdim%dim_data%size_global = COUNT(lmask(1:lmask_sz))
      ELSEIF (PRESENT(imask)) THEN
        pdim%dim_data%size_global = InputMaskGetTotal(dta=imask)
      ELSE
        pdim%dim_data%size_global = InputMaskGetTotal(mask)
      ENDIF
    ENDIF
    pdim%dim_data%size_local  =  pdim%dim_data%size_global
    pdim%dim_data%local_lo    =  1
    pdim%dim_data%local_hi    =  pdim%dim_data%size_global
    pdim%dim_data%chunk_lo    =  1
    pdim%dim_data%chunk_hi    =  pdim%dim_data%size_global
    IF (PRESENT(lmask)) THEN
      NewMask => InputMaskCreate(lmask,lmask_sz,pdim%dim_data%mask,Flags=INPUT_MASK_TYPE_MODEL+INPUT_MASK_GLOBAL)
    ELSEIF (PRESENT(imask)) THEN
      lDummy(:) = .FALSE.
      NewMask => InputMaskCreate(lDummy,1,pdim%dim_data%mask,mask=imask,Flags=INPUT_MASK_TYPE_MODEL+INPUT_MASK_GLOBAL)
    ELSE ! TODO: Not very nice - should check for presence in global list before adding
      CALL InputMaskCopy(new,mask)
      new%flags = INPUT_MASK_TYPE_MODEL+INPUT_MASK_GLOBAL
      new%next => AllInputMasks
      AllInputMasks => new%next
      NewMask => mask%mask
    ENDIF
    NewMask(1) = NewMask(1) ! To avoid compiler warnings

  END SUBROUTINE InputDimPackedAddMask

! Add a recognised packed dimension. Source dimensions (separated by commas) must already have been registered
  SUBROUTINE InputDimPackedAdd_gen(dim_name,src_dims,mask,lmask,lmask_sz,imask,alt_name,Extent,ldec,dim_ref,sub_model, &
                                   local_lo,local_hi,llocal,ilocal)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  CHARACTER (LEN=*),   INTENT(IN) :: src_dims
  TYPE (input_mask_list), POINTER, OPTIONAL :: mask
  LOGICAL, OPTIONAL,   INTENT(IN) :: lmask(*)
  INTEGER, OPTIONAL,   INTENT(IN) :: lmask_sz
  INTEGER, OPTIONAL,      POINTER :: imask(:)
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec
  CHARACTER (LEN=*),     OPTIONAL, INTENT(IN) :: alt_name
  INTEGER,               OPTIONAL, INTENT(IN) :: Extent, sub_model, local_lo, local_hi
  TYPE (input_dim_list), OPTIONAL,    POINTER :: dim_ref
  LOGICAL, OPTIONAL,      POINTER :: llocal(:)
  INTEGER, OPTIONAL,      POINTER :: ilocal(:)

  TYPE (input_mask_list), POINTER :: local_mask 
  TYPE (input_dim_list),  POINTER :: new, prev
  LOGICAL :: ldummy(1)
  INTEGER, POINTER :: mask_l(:)

    IF (lInLoop) CALL local_error('InputDimPackedAdd','Cannot add dimensions (packed or normal) after entering the main loop')

    ! Test if dimesion already exist
    new => InputDimGetRef(dim_name,alt_name,sub_model)
    IF (ASSOCIATED(new)) RETURN

    ! Create new packed dimension and add it to the end of the global list
    CALL InputDimNew(new)
    prev => AllInputDims
    IF (ASSOCIATED(prev)) THEN
      DO WHILE (ASSOCIATED(prev%next))
        prev => prev%next
      ENDDO
      prev%next => new
    ELSE
      AllInputDims => new
    ENDIF

    new%dim_data%name_dim = dim_name
    new%dim_data%flags    = IOR(new%dim_data%flags,INPUT_DIM_PACKED)
    new%src_dims => InputDimListFromString(src_dims,INPUT_DIM_SOURCE,'InputDimPackedAdd',sub_model,Sz=new%dim_data%nsrc)
    IF (PRESENT(lmask)) THEN
      CALL InputDimPackedAddMask(new,lmask=lmask,lmask_sz=lmask_sz,Extent=Extent)
    ELSEIF (PRESENT(imask)) THEN
      CALL InputDimPackedAddMask(new,imask=imask,Extent=Extent)
    ELSEIF (PRESENT(mask)) THEN
      CALL InputDimPackedAddMask(new,mask=mask,Extent=Extent)
    ELSE
      IF (PRESENT(Extent)) THEN
        new%dim_data%size_global = Extent
        new%dim_data%size_local  = Extent
        new%dim_data%local_lo    =      1
        new%dim_data%local_hi    = Extent
        new%dim_data%chunk_lo    =      1
        new%dim_data%chunk_hi    = Extent
      ENDIF
    ENDIF
    NULLIFY(local_mask)
    IF (PRESENT(llocal)) THEN
      IF (SIZE(llocal) /= new%dim_data%size_global) &
        CALL local_error('InputDimPackedAdd','Size of logical mask must match dimension size')
      mask_l => InputMaskCreate(llocal,SIZE(llocal),Ref=local_mask,Flags=INPUT_MASK_TYPE_PE+INPUT_MASK_GLOBAL)
      mask_l(1) = mask_l(1) ! To avoid compiler warnings
    ELSE
      IF (PRESENT(ilocal)) THEN
        ldummy(:) = .FALSE.
        mask_l => InputMaskCreate(ldummy,1,Ref=local_mask,mask=ilocal,Flags=INPUT_MASK_TYPE_PE+INPUT_MASK_GLOBAL)
        mask_l(1) = mask_l(1) ! To avoid compiler warnings
      ENDIF
    ENDIF
    IF ((PRESENT(local_lo) .AND. PRESENT(local_hi)) .OR. ASSOCIATED(local_mask)) &
      CALL InputDimLocalSet(local_lo,local_hi,local_mask,dim_ref=new)
    IF (PRESENT(alt_name))   new%dim_data%name_alt = alt_name
    IF (PRESENT(ldec)) THEN
      IF (ldec) THEN
        new%dim_data%flags = IOR(new%dim_data%flags,INPUT_DIM_DEC)
        new%flags          = IOR(new%flags,INPUT_DIM_DEC)
      ELSE
        new%dim_data%flags = IAND(new%dim_data%flags,INPUT_ALLFLAGS - INPUT_DIM_DEC)
        new%flags          = IAND(new%flags,INPUT_ALLFLAGS - INPUT_DIM_DEC)
      ENDIF
    ENDIF
    new%flags = IOR(new%flags,INPUT_DIM_GLOBAL)

    IF (PRESENT(dim_ref)) dim_ref => new

  END SUBROUTINE InputDimPackedAdd_gen

  SUBROUTINE InputDimPackedAddnd(dim_name,src_dims,imask,alt_name,Extent,ldec,dim_ref,sub_model,local_lo,local_hi,llocal,ilocal)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  CHARACTER (LEN=*),   INTENT(IN) :: src_dims
  INTEGER,                POINTER :: imask(:),ilocal(:)
  LOGICAL, OPTIONAL,      POINTER :: llocal(:)
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec
  CHARACTER (LEN=*),     OPTIONAL, INTENT(IN) :: alt_name
  INTEGER,               OPTIONAL, INTENT(IN) :: Extent, sub_model, local_lo, local_hi
  TYPE (input_dim_list), OPTIONAL,    POINTER :: dim_ref
    CALL InputDimPackedAdd_gen(dim_name,src_dims,imask=imask,alt_name=alt_name,Extent=Extent,ldec=ldec,dim_ref=dim_ref, &
                               sub_model=sub_model,local_lo=local_lo,local_hi=local_hi,llocal=llocal,ilocal=ilocal)
  END SUBROUTINE InputDimPackedAddnd

  SUBROUTINE InputDimPackedAdd1d(dim_name,src_dims,lmask,alt_name,Extent,ldec,dim_ref,sub_model,local_lo,local_hi,llocal,ilocal)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  CHARACTER (LEN=*),   INTENT(IN) :: src_dims
  LOGICAL,             INTENT(IN) :: lmask(:)
  LOGICAL, OPTIONAL,      POINTER :: llocal(:)
  INTEGER, OPTIONAL,      POINTER :: ilocal(:)
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec
  CHARACTER (LEN=*),     OPTIONAL, INTENT(IN) :: alt_name
  INTEGER,               OPTIONAL, INTENT(IN) :: Extent, sub_model, local_lo, local_hi
  TYPE (input_dim_list), OPTIONAL,    POINTER :: dim_ref
    CALL InputDimPackedAdd_gen(dim_name,src_dims,lmask=lmask,lmask_sz=SIZE(lmask),alt_name=alt_name,Extent=Extent,ldec=ldec, &
                               dim_ref=dim_ref,sub_model=sub_model,local_lo=local_lo,local_hi=local_hi,llocal=llocal,ilocal=ilocal)
  END SUBROUTINE InputDimPackedAdd1d

  SUBROUTINE InputDimPackedAdd2d(dim_name,src_dims,lmask,alt_name,Extent,ldec,dim_ref,sub_model,local_lo,local_hi,llocal,ilocal)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  CHARACTER (LEN=*),   INTENT(IN) :: src_dims
  LOGICAL,             INTENT(IN) :: lmask(:,:)
  LOGICAL, OPTIONAL,      POINTER :: llocal(:)
  INTEGER, OPTIONAL,      POINTER :: ilocal(:)
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec
  CHARACTER (LEN=*),     OPTIONAL, INTENT(IN) :: alt_name
  INTEGER,               OPTIONAL, INTENT(IN) :: Extent, sub_model, local_lo, local_hi
  TYPE (input_dim_list), OPTIONAL,    POINTER :: dim_ref
    CALL InputDimPackedAdd_gen(dim_name,src_dims,lmask=lmask,lmask_sz=SIZE(lmask),alt_name=alt_name,Extent=Extent,ldec=ldec, &
                               dim_ref=dim_ref,sub_model=sub_model,local_lo=local_lo,local_hi=local_hi,llocal=llocal,ilocal=ilocal)
  END SUBROUTINE InputDimPackedAdd2d

  SUBROUTINE InputDimPackedAdd3d(dim_name,src_dims,lmask,alt_name,Extent,ldec,dim_ref,sub_model,local_lo,local_hi,llocal,ilocal)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  CHARACTER (LEN=*),   INTENT(IN) :: src_dims
  LOGICAL,             INTENT(IN) :: lmask(:,:,:)
  LOGICAL, OPTIONAL,      POINTER :: llocal(:)
  INTEGER, OPTIONAL,      POINTER :: ilocal(:)
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec
  CHARACTER (LEN=*),     OPTIONAL, INTENT(IN) :: alt_name
  INTEGER,               OPTIONAL, INTENT(IN) :: Extent, sub_model, local_lo, local_hi
  TYPE (input_dim_list), OPTIONAL,    POINTER :: dim_ref
    CALL InputDimPackedAdd_gen(dim_name,src_dims,lmask=lmask,lmask_sz=SIZE(lmask),alt_name=alt_name,Extent=Extent,ldec=ldec, &
                               dim_ref=dim_ref,sub_model=sub_model,local_lo=local_lo,local_hi=local_hi,llocal=llocal,ilocal=ilocal)
  END SUBROUTINE InputDimPackedAdd3d

  SUBROUTINE InputDimPackedAdd4d(dim_name,src_dims,lmask,alt_name,Extent,ldec,dim_ref,sub_model,local_lo,local_hi,llocal,ilocal)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  CHARACTER (LEN=*),   INTENT(IN) :: src_dims
  LOGICAL,             INTENT(IN) :: lmask(:,:,:,:)
  LOGICAL, OPTIONAL,      POINTER :: llocal(:)
  INTEGER, OPTIONAL,      POINTER :: ilocal(:)
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec
  CHARACTER (LEN=*),     OPTIONAL, INTENT(IN) :: alt_name
  INTEGER,               OPTIONAL, INTENT(IN) :: Extent, sub_model, local_lo, local_hi
  TYPE (input_dim_list), OPTIONAL,    POINTER :: dim_ref
    CALL InputDimPackedAdd_gen(dim_name,src_dims,lmask=lmask,lmask_sz=SIZE(lmask),alt_name=alt_name,Extent=Extent,ldec=ldec, &
                               dim_ref=dim_ref,sub_model=sub_model,local_lo=local_lo,local_hi=local_hi,llocal=llocal,ilocal=ilocal)
  END SUBROUTINE InputDimPackedAdd4d

  SUBROUTINE InputDimPackedAdd5d(dim_name,src_dims,lmask,alt_name,Extent,ldec,dim_ref,sub_model,local_lo,local_hi,llocal,ilocal)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  CHARACTER (LEN=*),   INTENT(IN) :: src_dims
  LOGICAL,             INTENT(IN) :: lmask(:,:,:,:,:)
  LOGICAL, OPTIONAL,      POINTER :: llocal(:)
  INTEGER, OPTIONAL,      POINTER :: ilocal(:)
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec
  CHARACTER (LEN=*),     OPTIONAL, INTENT(IN) :: alt_name
  INTEGER,               OPTIONAL, INTENT(IN) :: Extent, sub_model, local_lo, local_hi
  TYPE (input_dim_list), OPTIONAL,    POINTER :: dim_ref
    CALL InputDimPackedAdd_gen(dim_name,src_dims,lmask=lmask,lmask_sz=SIZE(lmask),alt_name=alt_name,Extent=Extent,ldec=ldec, &
                               dim_ref=dim_ref,sub_model=sub_model,local_lo=local_lo,local_hi=local_hi,llocal=llocal,ilocal=ilocal)
  END SUBROUTINE InputDimPackedAdd5d

  SUBROUTINE InputDimPackedAdd6d(dim_name,src_dims,lmask,alt_name,Extent,ldec,dim_ref,sub_model,local_lo,local_hi,llocal,ilocal)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  CHARACTER (LEN=*),   INTENT(IN) :: src_dims
  LOGICAL,             INTENT(IN) :: lmask(:,:,:,:,:,:)
  LOGICAL, OPTIONAL,      POINTER :: llocal(:)
  INTEGER, OPTIONAL,      POINTER :: ilocal(:)
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec
  CHARACTER (LEN=*),     OPTIONAL, INTENT(IN) :: alt_name
  INTEGER,               OPTIONAL, INTENT(IN) :: Extent, sub_model, local_lo, local_hi
  TYPE (input_dim_list), OPTIONAL,    POINTER :: dim_ref
    CALL InputDimPackedAdd_gen(dim_name,src_dims,lmask=lmask,lmask_sz=SIZE(lmask),alt_name=alt_name,Extent=Extent,ldec=ldec, &
                               dim_ref=dim_ref,sub_model=sub_model,local_lo=local_lo,local_hi=local_hi,llocal=llocal,ilocal=ilocal)
  END SUBROUTINE InputDimPackedAdd6d

  SUBROUTINE InputDimPackedAdd7d(dim_name,src_dims,lmask,alt_name,Extent,ldec,dim_ref,sub_model,local_lo,local_hi,llocal,ilocal)
  CHARACTER (LEN=*),   INTENT(IN) :: dim_name
  CHARACTER (LEN=*),   INTENT(IN) :: src_dims
  LOGICAL,             INTENT(IN) :: lmask(:,:,:,:,:,:,:)
  LOGICAL, OPTIONAL,      POINTER :: llocal(:)
  INTEGER, OPTIONAL,      POINTER :: ilocal(:)
  LOGICAL, OPTIONAL,   INTENT(IN) :: ldec
  CHARACTER (LEN=*),     OPTIONAL, INTENT(IN) :: alt_name
  INTEGER,               OPTIONAL, INTENT(IN) :: Extent, sub_model, local_lo, local_hi
  TYPE (input_dim_list), OPTIONAL,    POINTER :: dim_ref
    CALL InputDimPackedAdd_gen(dim_name,src_dims,lmask=lmask,lmask_sz=SIZE(lmask),alt_name=alt_name,Extent=Extent,ldec=ldec, &
                               dim_ref=dim_ref,sub_model=sub_model,local_lo=local_lo,local_hi=local_hi,llocal=llocal,ilocal=ilocal)
  END SUBROUTINE InputDimPackedAdd7d

! Set the local size of a dimension given by either reference or name
  SUBROUTINE InputDimLocalSet(Lo,Hi,mask,lmask,imask,dim_name,dim_ref,sub_model)
  INTEGER,                OPTIONAL, INTENT(in) :: Lo, Hi
  TYPE (input_mask_list), OPTIONAL, POINTER    :: mask
  LOGICAL,                OPTIONAL, INTENT(in) :: lmask(:)
  INTEGER,                OPTIONAL, POINTER    :: imask(:)
  CHARACTER (len=*),      OPTIONAL, INTENT(in) :: dim_name
  TYPE (input_dim_list),  OPTIONAL, POINTER    :: dim_ref
  INTEGER,                OPTIONAL, INTENT(in) :: sub_model

  TYPE (input_mask_list), POINTER :: msk
  INTEGER :: tst_model, i, Extent, dum(1)
  INTEGER, POINTER :: im(:)
  REAL(dp), POINTER :: packed(:), unpacked(:)
  LOGICAL :: lSpec, lHasMask, ldummy(1)
  LOGICAL, POINTER :: lmsk(:)
  TYPE (input_dim_list), POINTER :: curr_dim, search_dim
  TYPE (input_var_list), POINTER :: curr_var

    tst_model = CurrSubModel
    IF (PRESENT(sub_model)) tst_model = sub_model

    ! Update size of dimension itself
    IF (PRESENT(dim_ref)) THEN
      curr_dim => dim_ref
    ELSE
      curr_dim => InputDimGetRef(dim_name,sub_model = tst_model)
      IF (.NOT. ASSOCIATED(curr_dim)) CALL local_error('InputDimLocalSet','Attempt to set local size of unknown dimension')
    ENDIF

    lHasMask = PRESENT(mask)
    IF (lHasMask) msk => mask

    IF (PRESENT(lmask)) THEN
      lHasMask = .TRUE.
      CALL MaskLogical2Int(im,lmask,size(lmask))
    ELSEIF (PRESENT(imask)) THEN
      im => imask
    ENDIF
    IF (PRESENT(lmask) .OR. PRESENT(imask)) THEN
      im => InputMaskCreate(ldummy,1,msk,mask=im)
    ENDIF
    IF (PRESENT(Lo) .AND. PRESENT(Hi)) THEN
      ALLOCATE(lmsk(curr_dim%dim_data%size_global))
      lmsk = .FALSE.
      lmsk(Lo:Hi) = .TRUE.
      im => InputMaskCreate(lmsk,SIZE(lmsk),msk)
      lHasMask = .TRUE.
      DEALLOCATE(lmsk)
    ENDIF

    lSpec = .FALSE.
    IF (lHasMask) THEN
      lSpec = .TRUE.
      Extent = InputMaskGetValid(msk)
      ALLOCATE(packed(Extent),unpacked(InputMaskGetTotal(msk)))
      DO i=1,Extent
        packed(i) = REAL(i,dp)
      ENDDO
      CALL MaskUnpack(packed,unpacked,msk%mask,Fill=0._dp)
      dum = MINLOC(unpacked,mask = unpacked > 0._dp)
      curr_dim%dim_data%local_lo = dum(1)
      dum = MAXLOC(unpacked)
      curr_dim%dim_data%local_hi = dum(1)
      curr_dim%dim_data%size_local = Extent
      DEALLOCATE(packed,unpacked)
      curr_dim%dim_data%local => msk
    ENDIF
    IF (PRESENT(Lo) .AND. PRESENT(Hi) .AND. .NOT. lSpec) THEN
      lSpec = .TRUE.
      Extent = Hi - Lo + 1
      IF ((curr_dim%dim_data%size_local >= 0) .AND. (Extent /= curr_dim%dim_data%size_local) .AND. &
          (curr_dim%dim_data%size_local /= curr_dim%dim_data%size_global)) THEN
        WRITE(message_text,'(3a,i3,a,i3)') 'Local size of ',TRIM(dim_name),' ',Extent, &
                                           ' is inconsistent with previous definition: ', curr_dim%dim_data%size_local
        CALL local_error('InputDimLocalSet',message_text)
      ENDIF
      curr_dim%dim_data%size_local = Extent
      curr_dim%dim_data%local_lo   = Lo
      curr_dim%dim_data%local_hi   = Hi
    ENDIF
    IF (.NOT. lSpec) CALL local_error('InputDimLocalSet','Either boundaries or mask must be specified')
    IF (curr_dim%dim_data%size_local /= curr_dim%dim_data%size_global) &
      curr_dim%dim_data%flags = IOR(curr_dim%dim_data%flags,INPUT_DIM_DISTRIB)

    ! Update size of associated variables
    curr_var => AllInputVars
    DO WHILE (ASSOCIATED(curr_var))
      IF ((curr_var%sub_model == tst_model) .AND. (IAND(curr_var%flags,INPUT_VAR_GLOBAL) == 0))  THEN
        search_dim => curr_var%dims
        i = 1
        DO WHILE (ASSOCIATED(search_dim))
          IF (ASSOCIATED(search_dim,curr_dim)) THEN
            curr_var%edims(i) = Extent
            IF (ASSOCIATED(curr_var%dta)) THEN
              DEALLOCATE(curr_var%dta)
              ALLOCATE(curr_var%dta(curr_var%edims(1),curr_var%edims(2),curr_var%edims(3),curr_var%edims(4),curr_var%edims(5), &
                                    curr_var%edims(6),curr_var%edims(7)))
              ! Needs manual reassociation of model var. Impossible to do since mo_input cannot not contain a proper reference
            ENDIF
            NULLIFY(search_dim)
          ELSE
            i = i + 1
            search_dim => search_dim%next
          ENDIF
        ENDDO
      ENDIF
      curr_var => curr_var%next
    ENDDO

  END SUBROUTINE InputDimLocalSet

! Creates a collective PE-distribution of multiple dimensions
  SUBROUTINE InputDimLocalSetMultiMask(dim_names,mask,sub_model,ndim)
  CHARACTER (len=*), INTENT(in)   :: dim_names
  TYPE (input_mask_list), POINTER :: mask
  INTEGER, OPTIONAL, INTENT(in)   :: sub_model
  INTEGER, OPTIONAL, INTENT(in)   :: ndim

  TYPE (input_dim_list), POINTER :: curr_dim, undef_dim
  CHARACTER (len=128) :: new_name, sub_dim
  INTEGER :: cnt, st, en, sz, sz_old, sz_all, i

    NULLIFY(undef_dim)
    new_name = TRIM(dim_names)
    sz_all = 1
    sz_old = 1
    cnt    = 0
    st     = 1
    DO WHILE (st<=LEN_TRIM(new_name))
      en = INDEX(new_name,',')
      IF (en<=0) THEN
        en = LEN_TRIM(new_name)+1
      ELSE
        new_name(en:en) = '_'
      ENDIF
      curr_dim => InputDimGetRef(new_name(st:en-1),sub_model=sub_model)
      IF (.NOT. ASSOCIATED(curr_dim)) CALL local_error('InputDimLocalSetMulti','Unknown source dimension: '//new_name(st:en-1))
      sz_all = sz_all * curr_dim%dim_data%size_global
      sub_dim = new_name(st:en-1)//'_iloc'
      CALL InputDimAdd(TRIM(sub_dim),dim_ref=curr_dim)
      IF (curr_dim%dim_data%size_global /= -1) THEN ! If this is an already defined dim with specified extent,accept previous size
        sz_old = sz_old * curr_dim%dim_data%size_global
      ELSE
        undef_dim => curr_dim
        curr_dim%dim_data%chunk_lo    =  1 ! Dim is defined to have size 1, since size of indiv. dims can't be determined
        curr_dim%dim_data%local_lo    =  1
        curr_dim%dim_data%chunk_hi    =  1
        curr_dim%dim_data%local_hi    =  1
        curr_dim%dim_data%size_local  =  1
        curr_dim%dim_data%size_global =  1
      ENDIF
      st  = en  + 1
      cnt = cnt + 1
    ENDDO

    ! Last undefined temporary dimension is set to be the only with a non-1 extent 
    ! (total extent is known, but not necessarily that of the individual dims)
    IF (ASSOCIATED(undef_dim)) THEN
      sz = INT(CEILING(REAL(sz_all,dp)/REAL(sz_old,dp)))
      undef_dim%dim_data%chunk_lo    =  1
      undef_dim%dim_data%local_lo    =  1
      undef_dim%dim_data%chunk_hi    = sz
      undef_dim%dim_data%local_hi    = mask%nValid/sz_old
      undef_dim%dim_data%size_local  = mask%nValid/sz_old
      undef_dim%dim_data%size_global = sz
    ELSE
      IF (sz_old < sz_all) CALL local_error('InputDimLocalSetMulti','All local dimensions defined, but total size is too small')
    ENDIF
    IF (PRESENT(ndim)) THEN
      IF (cnt /= ndim) CALL local_error('InputDimLocalSetMulti','Provided mask does not fit number of collective dimensions')
    ENDIF
    IF (mask%npts/=sz_all) CALL local_error('InputDimLocalSetMulti','Size mismatch between provided mask and specified dimensions')
    CALL InputDimAdd(TRIM(new_name),Extent=sz_all)
    sub_dim = ''
    st = 1
    DO i=1,cnt
      en = INDEX(dim_names(st:),',') + st - 1
      IF (en < st) en = LEN_TRIM(dim_names) + 1
      sub_dim = TRIM(sub_dim)//dim_names(st:en-1)//'_iloc,'
      st = en + 1
    ENDDO
    CALL InputDimEquivalence(TRIM(new_name),TRIM(dim_names),sub_model,AllInputColDimDistrib,NoAdd=.TRUE.)
    CALL InputDimLocalSet(dim_name=TRIM(new_name),mask=mask,sub_model=sub_model)
    CALL InputDimEquivalence(sub_dim(1:LEN_TRIM(sub_dim)-1),TRIM(new_name),sub_model,AllInputColDimDistrib)

  END SUBROUTINE InputDimLocalSetMultiMask

  SUBROUTINE InputDimLocalSetMultiInt(dim_names,imask,sub_model)
  CHARACTER (len=*), INTENT(in) :: dim_names
  INTEGER,              POINTER :: imask(:)
  INTEGER, OPTIONAL, INTENT(in) :: sub_model

  TYPE (input_mask_list), POINTER :: mask
  INTEGER, POINTER :: imask_ref(:)
  LOGICAL :: ldummy(1)

    imask_ref => InputMaskCreate(ldummy,1,mask,mask=imask,Flags=INPUT_MASK_TYPE_DISTRIB)
    imask_ref(1) = imask_ref(1) ! Avoid compiler warnings
    CALL InputDimLocalSetMultiMask(dim_names,mask,sub_model)

  END SUBROUTINE InputDimLocalSetMultiInt

  SUBROUTINE InputDimLocalSetMulti1d(dim_names,lmask,sub_model)
  CHARACTER (len=*), INTENT(in) :: dim_names
  LOGICAL,              POINTER :: lmask(:)
  INTEGER, OPTIONAL, INTENT(in) :: sub_model

  TYPE (input_mask_list), POINTER :: mask
  INTEGER, POINTER :: imask(:)

    imask => InputMaskCreate(lmask,SIZE(lmask),mask,Flags=INPUT_MASK_TYPE_DISTRIB)
    imask(1) = imask(1) ! Avoid compiler warnings
    CALL InputDimLocalSetMultiMask(dim_names,mask,sub_model,1)

  END SUBROUTINE InputDimLocalSetMulti1d

  SUBROUTINE InputDimLocalSetMulti2d(dim_names,lmask,sub_model)
  CHARACTER (len=*), INTENT(in) :: dim_names
  LOGICAL,              POINTER :: lmask(:,:)
  INTEGER, OPTIONAL, INTENT(in) :: sub_model

  TYPE (input_mask_list), POINTER :: mask
  INTEGER, POINTER :: imask(:)

    imask => InputMaskCreate(lmask,SIZE(lmask),mask,Flags=INPUT_MASK_TYPE_DISTRIB)
    imask(1) = imask(1) ! Avoid compiler warnings
    CALL InputDimLocalSetMultiMask(dim_names,mask,sub_model,2)

  END SUBROUTINE InputDimLocalSetMulti2d

  SUBROUTINE InputDimLocalSetMulti3d(dim_names,lmask,sub_model)
  CHARACTER (len=*), INTENT(in) :: dim_names
  LOGICAL,              POINTER :: lmask(:,:,:)
  INTEGER, OPTIONAL, INTENT(in) :: sub_model

  TYPE (input_mask_list), POINTER :: mask
  INTEGER, POINTER :: imask(:)

    imask => InputMaskCreate(lmask,SIZE(lmask),mask,Flags=INPUT_MASK_TYPE_DISTRIB)
    imask(1) = imask(1) ! Avoid compiler warnings
    CALL InputDimLocalSetMultiMask(dim_names,mask,sub_model,3)

  END SUBROUTINE InputDimLocalSetMulti3d

  SUBROUTINE InputDimLocalSetMulti4d(dim_names,lmask,sub_model)
  CHARACTER (len=*), INTENT(in) :: dim_names
  LOGICAL,              POINTER :: lmask(:,:,:,:)
  INTEGER, OPTIONAL, INTENT(in) :: sub_model

  TYPE (input_mask_list), POINTER :: mask
  INTEGER, POINTER :: imask(:)

    imask => InputMaskCreate(lmask,SIZE(lmask),mask,Flags=INPUT_MASK_TYPE_DISTRIB)
    imask(1) = imask(1) ! Avoid compiler warnings
    CALL InputDimLocalSetMultiMask(dim_names,mask,sub_model,4)

  END SUBROUTINE InputDimLocalSetMulti4d

  SUBROUTINE InputDimLocalSetMulti5d(dim_names,lmask,sub_model)
  CHARACTER (len=*), INTENT(in) :: dim_names
  LOGICAL,              POINTER :: lmask(:,:,:,:,:)
  INTEGER, OPTIONAL, INTENT(in) :: sub_model

  TYPE (input_mask_list), POINTER :: mask
  INTEGER, POINTER :: imask(:)

    imask => InputMaskCreate(lmask,SIZE(lmask),mask,Flags=INPUT_MASK_TYPE_DISTRIB)
    imask(1) = imask(1) ! Avoid compiler warnings
    CALL InputDimLocalSetMultiMask(dim_names,mask,sub_model,5)

  END SUBROUTINE InputDimLocalSetMulti5d

  SUBROUTINE InputDimLocalSetMulti6d(dim_names,lmask,sub_model)
  CHARACTER (len=*), INTENT(in) :: dim_names
  LOGICAL,              POINTER :: lmask(:,:,:,:,:,:)
  INTEGER, OPTIONAL, INTENT(in) :: sub_model

  TYPE (input_mask_list), POINTER :: mask
  INTEGER, POINTER :: imask(:)

    imask => InputMaskCreate(lmask,SIZE(lmask),mask,Flags=INPUT_MASK_TYPE_DISTRIB)
    imask(1) = imask(1) ! Avoid compiler warnings
    CALL InputDimLocalSetMultiMask(dim_names,mask,sub_model,6)

  END SUBROUTINE InputDimLocalSetMulti6d

  SUBROUTINE InputDimLocalSetMulti7d(dim_names,lmask,sub_model)
  CHARACTER (len=*), INTENT(in) :: dim_names
  LOGICAL,              POINTER :: lmask(:,:,:,:,:,:,:)
  INTEGER, OPTIONAL, INTENT(in) :: sub_model

  TYPE (input_mask_list), POINTER :: mask
  INTEGER, POINTER :: imask(:)

    imask => InputMaskCreate(lmask,SIZE(lmask),mask,Flags=INPUT_MASK_TYPE_DISTRIB)
    imask(1) = imask(1) ! Avoid compiler warnings
    CALL InputDimLocalSetMultiMask(dim_names,mask,sub_model,7)

  END SUBROUTINE InputDimLocalSetMulti7d

! Updates links between reading and recieving PE's
  SUBROUTINE InputDimSync()
  TYPE (input_dim_list),   POINTER :: curr_dim, src_dim
  TYPE (input_mask_list),  POINTER :: first, prev, curr, new
  TYPE (input_eqdim_list), POINTER :: last
  LOGICAL,  POINTER :: mask(:), local(:)
  INTEGER,  POINTER :: dst(:,:), dummy(:)
  REAL(dp), POINTER :: rmask(:), rlocal(:)
  INTEGER :: i, j, n, Sz(7), Hi(7), Lo(7), SrcLo(7), SrcHi(7), mn, mx
  LOGICAL :: lPar, lEquiv
  CHARACTER (len=128) :: str

    ! Since this subroutine is called only once at the end of the initialization, 
    ! linking of AllInputEqDims and AllInputColDimDistrib can be performed here
    last => AllInputEqDims
    IF (ASSOCIATED(last)) THEN
      DO WHILE (ASSOCIATED(last%next))
        last => last%next
      ENDDO
      last%next => AllInputColDimDistrib
    ELSE
      AllInputEqDims => AllInputColDimDistrib
    ENDIF
    ColDimAssociated = .TRUE.

    IF (n_pe <= 1) RETURN ! Not a parallel run => nothing to syncronize

    n = MAX(1,COUNT(src_pe==my_pe))
    ! Distribute local dimension information from all PEs to their respective data pusher PE
    curr_dim => AllInputDims
    DO WHILE (ASSOCIATED(curr_dim))
      ! Syncronize dimension settings if dimension belongs to current sub-model and is distributed across PEs
      IF ((curr_dim%dim_data%sub_model==CurrSubModel) .AND. (curr_dim%dim_data%size_local/=curr_dim%dim_data%size_global)) THEN
        ALLOCATE(dst(3,n))
        ! IOData-PEs recieve information from all PEs which it provides with data
        NULLIFY(first,prev)
        IF (liodata) THEN
          ALLOCATE(mask(curr_dim%dim_data%size_global))
          j = 0
          DO i=1,n_pe
            IF (src_pe(i) == my_pe) THEN
              j = j + 1
              IF (my_pe == i-1) THEN ! Don't communicate with yourself via MPI
                dst(:,j) = (/curr_dim%dim_data%local_lo,curr_dim%dim_data%local_hi,0/)
                IF (ASSOCIATED(curr_dim%dim_data%local)) dst(3,j) = SIZE(curr_dim%dim_data%local%mask)
              ELSE
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
                CALL p_recv(dst(:,j),i-1,1,3)
#endif
              ENDIF
              IF (dst(3,j) > 0) THEN ! If mask exist, use it...
                CALL InputMaskNew(new)
                ALLOCATE(new%mask(dst(3,j)))
                IF (my_pe == i-1) THEN
                  new%mask(:) = curr_dim%dim_data%local%mask(:)
                ELSE
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
                  CALL p_recv(new%mask,i-1,2,dst(3,j))
#endif
                ENDIF
              ELSE ! Otherwise construct manually...
                mask(:) = .FALSE.
                mask(dst(1,j):dst(2,j)) = .TRUE.
                dummy => InputMaskCreate(mask,SIZE(mask),new,Flags=INPUT_MASK_GLOBAL+INPUT_MASK_TYPE_PE)
              ENDIF
              IF (ASSOCIATED(prev)) THEN
                prev%next => new
              ELSE
                first => new
              ENDIF
              prev => new
            ENDIF
          ENDDO
        ELSE ! Non-IOData-PEs just send their own information
          dst(:,1) = (/curr_dim%dim_data%local_lo,curr_dim%dim_data%local_hi,0/)
          IF (ASSOCIATED(curr_dim%dim_data%local)) dst(3,1) = SIZE(curr_dim%dim_data%local%mask)
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
          CALL p_send(dst(:,1),src_pe(my_pe+1),1,3)
          IF (ASSOCIATED(curr_dim%dim_data%local)) &
            CALL p_send(curr_dim%dim_data%local%mask,src_pe(my_pe+1),2,SIZE(curr_dim%dim_data%local%mask))
#endif
        ENDIF
        DEALLOCATE(dst)
        ! IOData-PEs now have all necessary information, now determine chunk size and create relevant masks wrt. the chunk
        IF (liodata) THEN
          mask(:) = .FALSE.
          curr => first
          DO WHILE (ASSOCIATED(curr))
            CALL MaskInt2Logical(mask,curr%mask,lFalse=.FALSE.) ! Mark only .true. entries at each pass
            curr => curr%next
          ENDDO
          curr_dim%dim_data%chunk_lo = -1
          curr_dim%dim_data%chunk_hi = -1
          DO i=1,curr_dim%dim_data%size_global
            IF (mask(i)) THEN
              IF (curr_dim%dim_data%chunk_lo == -1) curr_dim%dim_data%chunk_lo = i
              curr_dim%dim_data%chunk_hi = i
            ENDIF
          ENDDO
          IF (ASSOCIATED(curr_dim%dim_data%local)) DEALLOCATE(curr_dim%dim_data%local) ! Potential source of a memory leak!
          ! This is a packed, distributed dimension - minimize the chunk of the source dimensions to be read and the packing mask
          IF (IAND(curr_dim%dim_data%flags,INPUT_DIM_PACKED) /= 0) THEN
            i     = 1
            Sz(:) = 1
            SrcLo(:) = 1
            SrcHi(:) = 1
            lPar  = .FALSE.
            src_dim => curr_dim%src_dims
            DO WHILE (ASSOCIATED(src_dim))
              Sz(i) = src_dim%dim_data%size_global
              ! The src dims are defined and thus processed before the packed, therefore the chunk size of src dims is reliable
              SrcLo(i) = src_dim%dim_data%chunk_lo
              SrcHi(i) = src_dim%dim_data%chunk_hi
              IF ((SrcLo(i) /= 1) .OR. (SrcHi(i) /= Sz(i))) lPar = .TRUE.
              i = i + 1
              src_dim => src_dim%next
            ENDDO
            new => curr_dim%dim_data%mask
            CALL InputMaskReduce(Sz,new,curr_dim%dim_data%chunk_lo,curr_dim%dim_data%chunk_hi,Lo,Hi, &
                                 INPUT_MASK_CHUNK+INPUT_MASK_TYPE_MODEL)
            ! If any of the source dimensions are distributed, build also a chunk mask for unpacking
            IF (lPar) THEN
              curr => curr_dim%dim_data%mask
              CALL InputMaskReduce(Sz,curr,SrcLo,SrcHi,mn,mx,INPUT_MASK_CHUNK+INPUT_MASK_TYPE_MODEL)
              new%next => curr ! Link as a second mask
              curr_dim%lo = mn
              curr_dim%hi = mx
            ENDIF
            ! Assign new chunk-packing mask(s) to dimension
            DEALLOCATE(curr_dim%dim_data%mask)
            curr_dim%dim_data%mask => new
            ! Assign alternative chunk to the source dimensions
            i = 1
            src_dim => curr_dim%src_dims
            DO WHILE (ASSOCIATED(src_dim))
              src_dim%lo = Lo(i)
              src_dim%hi = Hi(i)
              i = i + 1
              src_dim => src_dim%next
            ENDDO
          ENDIF
          ! Make a chain of distributing masks relative to the read chunk
          NULLIFY(curr_dim%dim_data%local,prev)
          curr => first
          DO WHILE (ASSOCIATED(curr))
            CALL MaskInt2Logical(mask,curr%mask)
            dummy => InputMaskCreate(mask(curr_dim%dim_data%chunk_lo:curr_dim%dim_data%chunk_hi), &
                                     curr_dim%dim_data%chunk_hi-curr_dim%dim_data%chunk_lo+1,new, &
                                     Flags=INPUT_MASK_CHUNK+INPUT_MASK_TYPE_PE)
            dummy(1) = dummy(1) ! To avoid compiler warnings
            IF (ASSOCIATED(prev)) THEN
              prev%next => new
            ELSE
              curr_dim%dim_data%local => new
            ENDIF
            prev => new
            curr => curr%next
          ENDDO
          ! Destroy local list of global PE masks
          CALL InputMasksDone(first,.FALSE.) ! Another potential source of a memory leak!
          NULLIFY(first)
          DEALLOCATE(mask)
        ENDIF
      ENDIF ! liodata

      ! Create local versions of packed dimensions with distributed source dimensions
      IF (IAND(curr_dim%dim_data%flags,INPUT_DIM_PACKED) /= 0) THEN
        lEquiv  = .FALSE.
        Lo(1:2) = 1
        Sz      = 1
        i       = 1
        n       = 1
        mn      = 1
        str     = ''
        src_dim => curr_dim%src_dims
        DO WHILE (ASSOCIATED(src_dim))
          ! Build name of evt. collectively distributed dimensions, build of same dimensions as the packed one
          mx = LEN_TRIM(src_dim%dim_data%name_dim)
          IF (mx > 4) THEN
            IF (src_dim%dim_data%name_dim(mx-4:mx)=='_iloc') THEN
              str = TRIM(str)//'_'//src_dim%dim_data%name_dim(1:mx-5)
              mn  = mn + 1
            ENDIF
          ENDIF
          Sz(i) = src_dim%dim_data%size_global
          IF (Sz(i) /= src_dim%dim_data%size_local) THEN
            n = n + 1
            first => src_dim%dim_data%local
          ELSE
            Lo(n) = Lo(n) * src_dim%dim_data%size_global
          ENDIF
          i = i + 1
          src_dim => src_dim%next
        ENDDO
        IF (n>2) CALL local_error('InputDimSync','Sorry - present version only allows packed dimensions ('// &
                                  TRIM(curr_dim%dim_data%name_dim)//') to have at most one distributed source dimension')
        IF (mn==i) THEN ! All source dimensions were '_iloc' dimensions => very likely to be a col. distrib. dim.
          src_dim => InputDimGetRef(TRIM(str(2:)))
          IF (.NOT. ASSOCIATED(src_dim) .AND. .NOT. ASSOCIATED(first)) CALL local_error('InputDimSync','Local packing of col'// &
            'lectively distributed dimensions is limited to the case that the source are the same as the distributed dimesions')
          first => src_dim%dim_data%local
          lEquiv = .TRUE.
        ENDIF
        ! For liodata PEs, the relevant mask may not be the first one
        IF (liodata .AND. ASSOCIATED(first)) THEN
          DO j=0,my_pe-1
            IF (src_pe(j+1) == my_pe) first => first%next
          ENDDO
        ENDIF
        ! Now create the local mask if necessary
        IF (ASSOCIATED(first)) THEN
          IF (lEquiv) THEN
            n = first%nValid
          ELSE
            n = first%nValid*Lo(1)*Lo(2)
          ENDIF
          i = PRODUCT(Sz)
          ALLOCATE(mask(i),rmask(i),local(n),rlocal(n))
          CALL MaskInt2Logical(mask,curr_dim%dim_data%mask%mask)
          rmask(:) = 0._dp
          WHERE(mask)
            rmask = 1._dp
          END WHERE
          IF (lEquiv) THEN
            CALL MaskPack(rmask,rlocal,first%mask)
          ELSE
            CALL MaskPack(rmask,rlocal,first%mask,Lo(1),Lo(2))
          ENDIF
          local(:) = rlocal(:)>0.5_dp 
          dummy => InputMaskCreate(local,n,curr,Flags=INPUT_MASK_LOCAL+INPUT_MASK_TYPE_MODEL)
          ! and add it to the head of the packed dimension mask list
          curr_dim%dim_data%size_local = curr%nValid
          curr_dim%dim_data%local_lo   = 1           ! Correct setting would require all PEs having
          curr_dim%dim_data%local_hi   = curr%nValid !   points before this one to communicate
          NULLIFY(new)
          CALL InputMaskCopy(new,curr)
          new%next => curr_dim%dim_data%mask
          curr_dim%dim_data%mask => new
          DEALLOCATE(mask,rmask,local,rlocal)
        ENDIF
      ENDIF
      curr_dim => curr_dim%next
    ENDDO

    ! Count total number of data pusher PEs
    ALLOCATE(mask(n_pe))
    mask(:) = .FALSE.
    DO i=1,n_pe
      mask(src_pe(i)+1) = .TRUE.
    ENDDO
    n_iope = COUNT(mask)
    DEALLOCATE(mask)

  END SUBROUTINE InputDimSync

! Declare list of dimensions to be equivalent
  SUBROUTINE InputDimEquivalence(dst_dims,src_dims,sub_model,list,NoAdd)
  CHARACTER (len=*), INTENT(in) :: dst_dims
  CHARACTER (len=*), INTENT(in) :: src_dims
  INTEGER, OPTIONAL, INTENT(in) :: sub_model
  TYPE (input_eqdim_list), OPTIONAL, POINTER :: list ! List to add the new equivalent dimension statement to
  LOGICAL, OPTIONAL, INTENT(in) :: NoAdd

  TYPE (input_eqdim_list), POINTER :: new
  TYPE (input_dim_list),   POINTER :: dcurr, scurr
  CHARACTER (len=5) :: add
  INTEGER :: SSz, DSz
  INTEGER, POINTER :: mask(:), dummy(:)
  LOGICAL :: ldummy(1)

    ! Create new element and put it into the global list
    ALLOCATE(new)
    NULLIFY(new%next,new%dst_dims,new%src_dims,new%mask)
    IF (PRESENT(list)) THEN
      new%next => list
      list => new
    ELSE
      new%next => AllInputEqDims ! Swaps order, but that doesn't matter
      AllInputEqDims => new
    ENDIF

    add = '_iloc'
    IF (PRESENT(NoAdd)) THEN
      IF (NoAdd) add = ''
    ENDIF
    new%dst_dims => InputDimListFromString(dst_dims,INPUT_DIM_EQ_DST,'InputDimEquivalence',sub_model,DSz)
    new%src_dims => InputDimListFromString(src_dims,INPUT_DIM_EQ_SRC,'InputDimEquivalence',sub_model,SSz,add=add)

    IF ((DSz > 0) .AND. (SSz > 0) .AND. (DSz < SSz)) CALL local_error('InputDimEquivalence', &
      'Destination dimensions ('//dst_dims//') have fewer points than source dimensions ('//src_dims//')')
    dcurr => new%dst_dims
    DO WHILE (ASSOCIATED(dcurr))
      scurr => new%src_dims
      DO WHILE (ASSOCIATED(scurr))
        IF (ASSOCIATED(dcurr%dim_data,scurr%dim_data)) CALL local_error('InputDimEquivalence','Dimension: '// &
          TRIM(dcurr%dim_data%name_dim)//' appears as both source and destination')
        scurr => scurr%next
      ENDDO
      dcurr => dcurr%next
    ENDDO
    ALLOCATE(mask(2))
    mask = (/SSz,DSz-SSz/)
    dummy => InputMaskCreate(ldummy,1,new%mask,mask=mask,Flags=INPUT_MASK_TYPE_EQUIV)
    dummy(1) = dummy(1) ! To avoid compiler warnings

  END SUBROUTINE InputDimEquivalence

! Display dimensions from a list
  SUBROUTINE InputEqDimDisp(eqdims,cnt,unt,msg)
  TYPE (input_eqdim_list), OPTIONAL, POINTER :: eqdims
  INTEGER, OPTIONAL, INTENT(in) :: cnt
  INTEGER, OPTIONAL, INTENT(in) :: unt
  CHARACTER (len=*), INTENT(in), OPTIONAL :: msg

  TYPE (input_eqdim_list),  POINTER :: curr
  TYPE (input_dim_list),    POINTER :: src
  CHARACTER (LEN=256) :: tmpstr
  INTEGER :: un, i, cr
  LOGICAL :: lAll

    un = GetLogUnit(); IF (PRESENT(unt)) un = unt
    cr = -1;           IF (PRESENT(cnt)) cr = cnt

    lAll = .FALSE.
    IF (PRESENT(eqdims)) THEN
      curr => eqdims
    ELSE
      curr => AllInputEqDims
      IF (ASSOCIATED(curr)) THEN
        lAll = .NOT. ColDimAssociated
      ELSE
        curr => AllInputColDimDistrib
      ENDIF
    ENDIF
    IF (PRESENT(msg)) WRITE(un,'(a)') msg
    DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
      tmpstr = ''
      src => curr%src_dims
      DO WHILE (ASSOCIATED(src))
        tmpstr = TRIM(tmpstr)//','//TRIM(src%dim_data%name_dim)
        src => src%next
      ENDDO
      tmpstr = TRIM(tmpstr(2:))//' =>;'
      src => curr%dst_dims
      DO WHILE (ASSOCIATED(src))
        tmpstr = TRIM(tmpstr)//TRIM(src%dim_data%name_dim)//','
        src => src%next
      ENDDO
      i = INDEX(tmpstr,';')
      tmpstr(i:i) = ' '
      WRITE(un,*) tmpstr(:LEN_TRIM(tmpstr)-1)
      curr => curr%next
      IF (lAll .AND. .NOT. ASSOCIATED(curr)) THEN
        curr => AllInputColDimDistrib ! If displaying all dims - display also those needed for collective dim distribution
        lAll = .FALSE.
      ENDIF
      cr = cr - 1
    ENDDO

  END SUBROUTINE InputEqDimDisp

END MODULE mo_input_dimension
