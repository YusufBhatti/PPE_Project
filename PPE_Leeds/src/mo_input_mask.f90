!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! This module contains procedures for mo_input dealing with masks
! Never directly use this module from outside the mo_input framework - only via mo_input.f90
MODULE mo_input_mask
USE mo_input_strings
USE mo_input_arrays
USE mo_input_types
USE mo_input_dimension_base

INTERFACE InputMaskReduce
  MODULE PROCEDURE InputMaskReduce_Single
  MODULE PROCEDURE InputMaskReduce_Multi
END INTERFACE

PUBLIC InputMaskNew, InputMaskCopy, InputMaskDisp, InputMasksDone, InputMaskCreate, InputMaskGetTotal, InputMaskGetValid
PUBLIC InputMaskCreateScatter, InputMaskReduce

CONTAINS

! ----------------------------------------------------
!  Creating, adding, setting, modifying, displaying:
!     Masks
! ----------------------------------------------------

! New input_mask_list element with default element
  SUBROUTINE InputMaskNew(new)
  TYPE (input_mask_list), POINTER :: new

    ALLOCATE(new)
    NULLIFY(new%mask,new%next)
    new%npts   = 0
    new%nvalid = 0
    new%flags  = 0

  END SUBROUTINE InputMaskNew

! Creates a copy of a mask, using the same mask data
  SUBROUTINE InputMaskCopy(dst,src)
  TYPE (input_mask_list), POINTER :: dst, src

    IF (.NOT. ASSOCIATED(dst)) CALL InputMaskNew(dst)
    dst%npts   =  src%npts
    dst%nvalid =  src%nvalid
    dst%flags  =  src%flags
    dst%mask   => src%mask
    NULLIFY(dst%next)

  END SUBROUTINE InputMaskCopy

! Display masks from a list
  SUBROUTINE InputMaskDisp(masks,cnt,unt,msg)
  TYPE (input_mask_list), OPTIONAL, POINTER :: masks
  INTEGER, OPTIONAL, INTENT(in) :: cnt
  INTEGER, OPTIONAL, INTENT(in) :: unt
  CHARACTER (len=*), INTENT(in), OPTIONAL :: msg

  CHARACTER (len=256) :: tmpstr
  TYPE (input_mask_list), POINTER :: curr
  INTEGER :: un, Sz, cr

    un = GetLogUnit(); IF (PRESENT(unt)) un = unt
    cr = -1;           IF (PRESENT(cnt)) cr = cnt

    IF (PRESENT(masks)) THEN
      curr => masks
    ELSE
      curr => AllInputMasks
    ENDIF
    IF (PRESENT(msg)) WRITE(un,*) msg
    IF (ASSOCIATED(curr)) WRITE(un,*) 'Mask size, Src size, Dst size (flags)'
    DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
      CALL DispFlagsSingleBit(curr%flags,INPUT_MASK_CONT,flg_mask_n,flg_mask,Txt=message_text,lEmpty=.TRUE.)
      tmpstr = TRIM(GetStringIndexed(mask_type,IAND(curr%flags,INPUT_MASK_TYPE_MASK)  /INPUT_MASK_TYPE_DISTRIB+1,    ' '))//' '// &
               TRIM(GetStringIndexed(flg_mask, IAND(curr%flags,INPUT_MASK_DOMAIN_MASK)/INPUT_MASK_CHUNK+flg_mask_n+1,' '))//' '// &
               TRIM(message_text)
      Sz = -1
      IF (ASSOCIATED(curr%mask)) Sz=SIZE(curr%mask)
      WRITE(message_text,*) Sz,curr%npts,curr%nvalid,' (',TRIM(tmpstr),')'
      CALL RmDblSpc(message_text)
      WRITE(un,'(a)') TRIM(message_text)

      curr => curr%next
      cr = cr - 1
    ENDDO

  END SUBROUTINE InputMaskDisp

! Destroys a list of mask, eventually the associated data
  SUBROUTINE InputMasksDone(masks,lDel)
  TYPE (input_mask_list), POINTER :: masks
  LOGICAL, OPTIONAL,   INTENT(in) :: lDel

  TYPE (input_mask_list), POINTER :: curr, next
  LOGICAL :: Del

    Del = .TRUE.
    IF (PRESENT(lDel)) Del = lDel

    curr => masks
    DO WHILE (ASSOCIATED(curr))
      next => curr%next
      IF (ASSOCIATED(curr%mask) .AND. Del) DEALLOCATE(curr%mask)
      DEALLOCATE(curr)
      curr => next
    ENDDO

  END SUBROUTINE InputMasksDone

! Creates new integer mask from a logical. If the resulting mask is not previously known, it is added to the
! list of known masks, otherwise the just created mask is thrown away and the old one is returned
  FUNCTION InputMaskCreate(Src,Sz,Ref,lOld,mask,Flags)
  INTEGER, POINTER :: InputMaskCreate(:)
  LOGICAL, DIMENSION(*) :: Src
  INTEGER, INTENT(in)   :: Sz
  INTEGER, OPTIONAL, POINTER :: mask(:)
  TYPE (input_mask_list), OPTIONAL, POINTER :: ref
  LOGICAL, OPTIONAL, INTENT(OUT) :: lOld
  INTEGER, OPTIONAL, INTENT(IN)  :: Flags

  TYPE (input_mask_list), POINTER :: curr, new
  INTEGER, POINTER :: Tmp(:)
  INTEGER :: dp
  LOGICAL :: lCont

    IF (PRESENT(mask)) THEN
      Tmp => mask
    ELSE
      CALL MaskLogical2Int(Tmp,Src,Sz)
    ENDIF
    dp = SIZE(Tmp)

    ! Check if this mask is already in the list
    NULLIFY(new)
    curr => AllInputMasks
    DO WHILE (ASSOCIATED(curr))
      lCont = SIZE(curr%mask)==dp
      IF (lCont) lCont = ALL(curr%mask(:)-Tmp(:)==0)
      IF (lCont) THEN
        IF (.NOT. PRESENT(mask)) DEALLOCATE(Tmp)
        IF (PRESENT(ref)) THEN
          CALL InputMaskCopy(new,curr)
          ref  => new
        ENDIF
        IF (PRESENT(lOld)) lOld = .TRUE.
        InputMaskCreate => curr%mask
        RETURN
      ENDIF
      curr => curr%next
    ENDDO

    ! This is a new mask
    CALL InputMaskNew(curr)
    curr%next     => AllInputMasks ! Swaps order, but order doesn't matter
    AllInputMasks => curr

    IF (PRESENT(mask)) THEN
      ALLOCATE(curr%mask(dp))
      curr%mask(:) = Tmp(:)
    ELSE
      curr%mask => Tmp
    ENDIF

    ! Necessary mask statistics
    IF (PRESENT(Flags)) curr%flags = Flags
    curr%npts   = MaskGetTotal(curr%mask)
    curr%nvalid = MaskGetValid(curr%mask)
    IF (MaskIsCont(curr%mask,dp)) curr%flags = IOR(curr%flags,INPUT_MASK_CONT + dp)

    IF (PRESENT(ref)) THEN
      CALL InputMaskCopy(new,curr)
      ref => new
    ENDIF
    IF (PRESENT(lOld)) lOld = .FALSE.
    InputMaskCreate => curr%mask

  END FUNCTION InputMaskCreate

! Number of points for an unpacked field
  INTEGER FUNCTION InputMaskGetTotal(mask,dta)
  TYPE (input_mask_list), POINTER, OPTIONAL :: mask
  INTEGER,             INTENT(in), OPTIONAL :: dta(:)

    IF (PRESENT(mask)) THEN
      InputMaskGetTotal = mask%npts
    ELSE
      InputMaskGetTotal = MaskGetTotal(dta)
    ENDIF

  END FUNCTION InputMaskGetTotal

! Number of points for a packed field
  INTEGER FUNCTION InputMaskGetValid(mask,dta)
  TYPE (input_mask_list), POINTER, OPTIONAL :: mask
  INTEGER,             INTENT(in), OPTIONAL :: dta(:)

    IF (PRESENT(mask)) THEN
      InputMaskGetValid = mask%nValid
    ELSE
      InputMaskGetValid = MaskGetValid(dta)
    ENDIF

  END FUNCTION InputMaskGetValid

! Internal subroutine of InputMaskCreateScatter for processing a single scattered dimension
  RECURSIVE SUBROUTINE MaskCreateScatter(dims,masks,nmask,ndim,Tst,Lo)
  TYPE (input_dim_list),  POINTER :: dims  ! INTENT(in)
  TYPE (input_mask_list), POINTER :: masks ! INTENT(out)
  INTEGER,            INTENT(out) :: nmask
  INTEGER,            INTENT(in)  :: ndim
  INTEGER,            INTENT(in)  :: Lo(7)
  LOGICAL,                POINTER :: Tst(:,:,:,:,:,:,:)

  TYPE (input_mask_list), POINTER :: mask, prev, first, curr
  TYPE (input_dim_list),  POINTER :: curr_dim
  LOGICAL,                POINTER :: bak(:,:,:,:,:,:,:), lmask(:)
  INTEGER,                POINTER :: imask(:)
  INTEGER                         :: Sz(7), n, nm, mul, rep

    nmask = 0
    Sz(:) = SHAPE(Tst)
    n = ndim
    curr_dim => dims
    DO WHILE (ASSOCIATED(curr_dim))
      IF (curr_dim%dim_data%size_local /= curr_dim%dim_data%size_global) THEN ! This a dimension to scatter
        ! Make a backup of the test array
        ALLOCATE(bak(Sz(1),Sz(2),Sz(3),Sz(4),Sz(5),Sz(6),Sz(7)))
        bak(:,:,:,:,:,:,:) = Tst(:,:,:,:,:,:,:)
        ! Loop over destination PEs
        NULLIFY(first,prev)
        nm = 0
        curr => curr_dim%dim_data%local
        DO WHILE (ASSOCIATED(curr))
          curr%flags = IAND(curr%flags,INPUT_ALLFLAGS - INPUT_MASK_PROCESSED)
          curr => curr%next
        ENDDO
        curr => curr_dim%dim_data%local
        DO WHILE (ASSOCIATED(curr))
          IF (IAND(curr%flags,INPUT_MASK_PROCESSED) == 0) THEN
            ! Set parts of Tst which do not belong on PE to false, leave the rest as they were
            Tst(:,:,:,:,:,:,:) = bak(:,:,:,:,:,:,:)
            mul = PRODUCT(Sz(1:n-1))
            rep = PRODUCT(Sz(n+1:7))
            imask => curr%mask
            IF (curr%npts /= Sz(n)) THEN ! Mask referes to a chunk which is not in use here (e.g. equiv. dims)
              ALLOCATE(lmask(Sz(n)))     ! Remap mask to cover global size of dimension
              lmask(:) = .FALSE.
              CALL MaskInt2Logical(lmask(Lo(n):),imask)
              CALL MaskLogical2Int(imask,lmask,Sz(n))
            ENDIF
            CALL MaskInt2Logical(Tst,imask,lTrue=.FALSE.,mul=mul,rep=rep) 
            CALL MaskCreateScatter(curr_dim%next,mask,nm,ndim-1,Tst,Lo) ! Set Tst according to further distributed dimensions
            IF (curr%npts /= Sz(n)) DEALLOCATE(lmask,imask)
            nmask = nmask + nm
            IF (ASSOCIATED(prev)) THEN
              prev%next => mask
            ELSE
              first     => mask
            ENDIF
            prev => mask
            DO WHILE (ASSOCIATED(prev%next))
              prev => prev%next
            ENDDO
            ! Mark all PEs with this mask as been processed
            mask => curr_dim%dim_data%local
            DO WHILE (ASSOCIATED(mask))
              ! Since all mask data here are in the global list, it is sufficient to test for associations
              IF (ASSOCIATED(mask%mask,curr%mask)) curr%flags = IOR(curr%flags,INPUT_MASK_PROCESSED)
              mask => mask%next
            ENDDO
          ENDIF
          ! Escape if all necessary masks have been created by the recursive call. Is this always the right thing to do?
          IF (nm <= 1) THEN
            curr => curr%next
          ELSE
            NULLIFY(curr)
          ENDIF
        ENDDO
        masks => first
        DEALLOCATE(bak)
        RETURN ! The rest is beeing taken care of in the recursive call
      ELSE
        curr_dim => curr_dim%next
        n = n - 1
      ENDIF
    ENDDO

    ! No more distributed dimensions, now actually create the mask
    imask => InputMaskCreate(Tst,PRODUCT(Sz),mask,Flags=INPUT_MASK_TYPE_DISTRIB+INPUT_MASK_CHUNK)
    IF (.NOT. ASSOCIATED(imask)) CALL local_error('MaskCreateScatter','Could not create mask')
    nmask =  1
    masks => mask

  END SUBROUTINE MaskCreateScatter

! Creates masks for distributing dimensions with MPI. In case of multiple scattered dimensions, a PE mesh is assumed
  SUBROUTINE InputMaskCreateScatter(dims,masks,nmask)
  TYPE (input_dim_list),  POINTER :: dims  ! INTENT(in)
  TYPE (input_mask_list), POINTER :: masks ! INTENT(out)
  INTEGER,            INTENT(out) :: nmask

  TYPE (input_dim_list),  POINTER :: curr_dim, rev_dim, first
  TYPE (input_mask_list), POINTER :: mask
  LOGICAL,                POINTER :: lmask(:,:,:,:,:,:,:)
  INTEGER :: Sz(7), Lo(7), n

    NULLIFY(masks)
    nmask = 0
    IF (.NOT. liodata .OR. (n_pe == 1)) RETURN ! Scatter masks are only of interest to IOData-PEs during parallel runs

    ! Find size of data read on this PE and make a copy of the dimensions in reverse order
    NULLIFY(first)
    Sz(:) = 1
    Lo(:) = 1
    n = 0
    curr_dim => dims
    DO WHILE (ASSOCIATED(curr_dim))
      n = n + 1
      IF (n > 7) CALL local_error('InputMaskCreateScatter','Too many dimensions')
      IF ((IAND(curr_dim%flags,INPUT_DIM_USEASGLOBAL) /= 0) .AND. (curr_dim%lo > 0)) THEN
        Sz(n) = curr_dim%hi - curr_dim%lo + 1
        Lo(n) = curr_dim%dim_data%chunk_lo
      ELSE
        Sz(n) = curr_dim%dim_data%chunk_hi - curr_dim%dim_data%chunk_lo + 1
      ENDIF
      CALL InputDimNewCtl(rev_dim)
      rev_dim%flags    =  IAND(curr_dim%flags,INPUT_ALLFLAGS - INPUT_DIM_TYPE) + INPUT_DIM_MASK
      rev_dim%dim_data => curr_dim%dim_data
      rev_dim%next     => first
      first            => rev_dim
      curr_dim         => curr_dim%next
    ENDDO

    ! Create base array
    ALLOCATE(lmask(Sz(1),Sz(2),Sz(3),Sz(4),Sz(5),Sz(6),Sz(7)))
    lmask(:,:,:,:,:,:,:) = .TRUE.
    CALL MaskCreateScatter(first,mask,nmask,n,lmask,lo)
    CALL InputDimDone(first,.FALSE.,.FALSE.)
    DEALLOCATE(lmask)
    masks => mask

  END SUBROUTINE InputMaskCreateScatter

  SUBROUTINE InputMaskReduce_Multi(Sz,mask,mn,mx,lo,hi,Flags)
  INTEGER,            INTENT(in ) :: Sz(7)        ! Shape of a source array for the mask
  TYPE (input_mask_list), POINTER :: mask         ! Mask to be reduced. Will on output hold the reduced mask. Original preserved
  INTEGER,            INTENT(in ) :: mn(7), mx(7) ! Minimum and maximum allowed points in each dimension
  INTEGER, OPTIONAL,  INTENT(out) :: lo, hi       ! Minimum and maximum of index within range
  INTEGER, OPTIONAL,  INTENT(in ) :: Flags        ! Flags for mask in global list

  TYPE (input_mask_list), POINTER :: res
  REAL(dp), POINTER :: One2n(:), FullField(:,:,:,:,:,:,:)
  INTEGER, POINTER :: dummy(:)
  INTEGER :: i, n

    n = PRODUCT(Sz)
    ALLOCATE(One2n(n),FullField(Sz(1),Sz(2),Sz(3),Sz(4),Sz(5),Sz(6),Sz(7)))
    DO i=1,n
      One2n(i) = REAL(i,dp)
    ENDDO

    CALL MaskUnpack(One2n,FullField,mask%mask,Fill=0._dp)
    IF (PRESENT(lo) .OR. PRESENT(hi)) THEN
      FullField(1:mn(1)-1,:,:,:,:,:,:) = 0._dp
      FullField(:,1:mn(2)-1,:,:,:,:,:) = 0._dp
      FullField(:,:,1:mn(3)-1,:,:,:,:) = 0._dp
      FullField(:,:,:,1:mn(4)-1,:,:,:) = 0._dp
      FullField(:,:,:,:,1:mn(5)-1,:,:) = 0._dp
      FullField(:,:,:,:,:,1:mn(6)-1,:) = 0._dp
      FullField(:,:,:,:,:,:,1:mn(7)-1) = 0._dp
      FullField(mx(1)+1:Sz(1),:,:,:,:,:,:) = 0._dp
      FullField(:,mx(2)+1:Sz(2),:,:,:,:,:) = 0._dp
      FullField(:,:,mx(3)+1:Sz(3),:,:,:,:) = 0._dp
      FullField(:,:,:,mx(4)+1:Sz(4),:,:,:) = 0._dp
      FullField(:,:,:,:,mx(5)+1:Sz(5),:,:) = 0._dp
      FullField(:,:,:,:,:,mx(6)+1:Sz(6),:) = 0._dp
      FullField(:,:,:,:,:,:,mx(7)+1:Sz(7)) = 0._dp
    ENDIF
    IF (PRESENT(lo)) lo = INT(MINVAL(FullField,mask=FullField>0._dp))
    IF (PRESENT(hi)) hi = INT(MAXVAL(FullField))
    dummy => InputMaskCreate(FullField(mn(1):mx(1),mn(2):mx(2),mn(3):mx(3),mn(4):mx(4),mn(5):mx(5),mn(6):mx(6),mn(7):mx(7)) &
                               /=0._dp,PRODUCT(mx(:)-mn(:)+1),Ref=res,Flags=Flags)
    dummy(1) = dummy(1) ! Avoid compiler warning
    DEALLOCATE(One2n,FullField)
    mask => res

  END SUBROUTINE InputMaskReduce_Multi

! Reduce mask to other coordinate frame (i.e. global => chunk, global => local or chunk => local)
  SUBROUTINE InputMaskReduce_Single(Sz,mask,mn,mx,lo,hi,Flags)
  INTEGER,            INTENT(in ) :: Sz(7)        ! Shape of a source array for the mask
  TYPE (input_mask_list), POINTER :: mask         ! Mask to be reduced. Will on output hold the reduced mask. Original preserved
  INTEGER,            INTENT(in ) :: mn, mx       ! Minimum and maximum allowed points in output mask
  INTEGER, OPTIONAL,  INTENT(out) :: lo(7), hi(7) ! Minimum and maximum of each dimension which hold data
  INTEGER, OPTIONAL,  INTENT(in ) :: Flags        ! Flags for mask in global list

  TYPE (input_mask_list), POINTER :: res
  REAL(dp), POINTER :: One2n(:), FullField(:,:,:,:,:,:,:)
  INTEGER, POINTER :: dummy(:)
  INTEGER :: i, n, l(7), h(7)

    n = PRODUCT(Sz)
    ALLOCATE(One2n(n),FullField(Sz(1),Sz(2),Sz(3),Sz(4),Sz(5),Sz(6),Sz(7)))
    DO i=1,n
      One2n(i) = REAL(i,dp)
    ENDDO

    CALL MaskUnpack(One2n,FullField,mask%mask,Fill=0._dp)
    WHERE (FullField(:,:,:,:,:,:,:) < REAL(mn,dp) .OR. FullField(:,:,:,:,:,:,:) > REAL(mx,dp))
      FullField(:,:,:,:,:,:,:) = 0._dp
    END WHERE
    l(:) = 1
    h(:) = Sz(:)
    DO WHILE (l(1)<Sz(1) .AND. ALL(FullField(l(1),:,:,:,:,:,:)==0._dp))
      l(1) = l(1) + 1
    ENDDO
    DO WHILE (l(2)<Sz(2) .AND. ALL(FullField(:,l(2),:,:,:,:,:)==0._dp))
      l(2) = l(2) + 1
    ENDDO
    DO WHILE (l(3)<Sz(3) .AND. ALL(FullField(:,:,l(3),:,:,:,:)==0._dp))
      l(3) = l(3) + 1
    ENDDO
    DO WHILE (l(4)<Sz(4) .AND. ALL(FullField(:,:,:,l(4),:,:,:)==0._dp))
      l(4) = l(4) + 1
    ENDDO
    DO WHILE (l(5)<Sz(5) .AND. ALL(FullField(:,:,:,:,l(5),:,:)==0._dp))
      l(5) = l(5) + 1
    ENDDO
    DO WHILE (l(6)<Sz(6) .AND. ALL(FullField(:,:,:,:,:,l(6),:)==0._dp))
      l(6) = l(6) + 1
    ENDDO
    DO WHILE (l(7)<Sz(7) .AND. ALL(FullField(:,:,:,:,:,:,l(7))==0._dp))
      l(7) = l(7) + 1
    ENDDO
    DO WHILE (h(1) >  1  .AND. ALL(FullField(h(1),:,:,:,:,:,:)==0._dp))
      h(1) = h(1) - 1
    ENDDO
    DO WHILE (h(2) >  1  .AND. ALL(FullField(:,h(2),:,:,:,:,:)==0._dp))
      h(2) = h(2) - 1
    ENDDO
    DO WHILE (h(3) >  1  .AND. ALL(FullField(:,:,h(3),:,:,:,:)==0._dp))
      h(3) = h(3) - 1
    ENDDO
    DO WHILE (h(4) >  1  .AND. ALL(FullField(:,:,:,h(4),:,:,:)==0._dp))
      h(4) = h(4) - 1
    ENDDO
    DO WHILE (h(5) >  1  .AND. ALL(FullField(:,:,:,:,h(5),:,:)==0._dp))
      h(5) = h(5) - 1
    ENDDO
    DO WHILE (h(6) >  1  .AND. ALL(FullField(:,:,:,:,:,h(6),:)==0._dp))
      h(6) = h(6) - 1
    ENDDO
    DO WHILE (h(7) >  1  .AND. ALL(FullField(:,:,:,:,:,:,h(7))==0._dp))
      h(7) = h(7) - 1
    ENDDO
    IF (PRESENT(lo) .AND. PRESENT(hi)) THEN
      lo(:) = l(:)
      hi(:) = h(:)
    ENDIF
    dummy => InputMaskCreate(FullField(l(1):h(1),l(2):h(2),l(3):h(3),l(4):h(4),l(5):h(5),l(6):h(6),l(7):h(7))/=0._dp, &
                             PRODUCT(h(:)-l(:)+1),Ref=res,Flags=Flags)
    dummy(1) = dummy(1) ! Avoid compiler warning
    DEALLOCATE(One2n,FullField)
    mask => res

  END SUBROUTINE InputMaskReduce_Single

END MODULE mo_input_mask
