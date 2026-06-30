!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! This module contains procedures for mo_input to performs data transformation on the read data 
! Never directly use this module from outside the mo_input framework - only via mo_input.f90
MODULE mo_input_proc
USE mo_input_strings
USE mo_input_arrays
USE mo_input_types
USE mo_input_data
USE mo_input_base
USE mo_input_mask
USE mo_input_netcdf
IMPLICIT NONE

PUBLIC :: GroupReadMsg
PUBLIC :: DataProcNoAction, DataProcMessage, DataProcRead, DataProcReadMulti, DataProcPermute, DataProcReverse, DataProcSwap
PUBLIC :: DataProcPack, DataProcUnpack, DataProcEquiv, DataProcCopy, DataProcSum, DataProcAvg, DataProcNormalize
PUBLIC :: DataProcFillBad, DataProcFillOutOfRange, DataProcRescale, DataProcTimeDiv, DataProcSpread, DataProcSubset
PUBLIC :: DataProcInterpolateTimeLinear, DataProcBroadcast, DataProcScatter
PUBLIC :: DataProcAdd, DataProcSub, DataProcSubr, DataProcMul, DataProcDiv, DataProcDivr
PUBLIC :: DataProcFlux, DataProcNudge, DataProcNudgeFld
#ifdef HAVE_F2003
PUBLIC :: InputVarGetNextBuffer, InputVarGetNextMask
#else
PUBLIC :: InputDoDataActions
#endif

CONTAINS

! -------------------------------------------
! Processing of data from file
! Mainly "plug-in" routines designed to be called via function pointers
! -------------------------------------------

! Writes the message, that a group has been read
  SUBROUTINE GroupReadMsg(ActGroup)
  TYPE (input_group_list), POINTER :: ActGroup

  TYPE (input_var_list), POINTER :: curr
  INTEGER :: add
  CHARACTER (len=2) :: sep

    IF (ldbg .AND. (     (IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_UPDATE_GROUP) &
                    .OR. (IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_READ_GROUP))  &
             .AND. (IAND(ActGroup%flags,INPUT_GROUP_MSG_DONE) == 0)) THEN
      curr => AllInputVars
      DO WHILE (ASSOCIATED(curr))
        IF (ASSOCIATED(curr%group,ActGroup)) EXIT
        curr => curr%next
      ENDDO
      IF (.NOT. ASSOCIATED(curr)) RETURN
      IF (curr%nToRead > 0) THEN
        IF (ActGroup%curr_rec >= 0) ActGroup%curr_rec = -1
      ELSE
        IF (ActGroup%curr_rec < 0) THEN
          WRITE(message_text,'(5a)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Read: ', &
                                     TRIM(ActGroup%name_group)
        ELSE
          WRITE(message_text,'(5a)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Read: Valid: ', &
            TRIM(DateString(ActGroup%time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': ',TRIM(ActGroup%name_group)
        ENDIF
        add = 1
        sep = ' ('
        IF (ASSOCIATED(ActGroup%group_file) .AND. (IAND(ActGroup%flags,INPUT_GROUP_MULTI_FILE)==0)) THEN
          IF (LEN_TRIM(ActGroup%group_file%name_file) > 0) THEN
            WRITE(message_text(LEN_TRIM(message_text)+1:),'(3a)') ' (',TRIM(ActGroup%group_file%name_file),')'
            add = 0
            sep = ', '
          ENDIF
        ENDIF
        IF (ActGroup%curr_rec /= UNLIKELY_VAL) THEN
          IF (ActGroup%curr_rec >= 0) THEN
            WRITE(message_text(LEN_TRIM(message_text)+add:),'(a,i10,a)') sep,ActGroup%curr_rec,')'
          ELSE
            WRITE(message_text(LEN_TRIM(message_text)+add:),'(a,i10,a)') sep,-ActGroup%curr_rec,' records)'
            ActGroup%curr_rec = 0
          ENDIF
        ENDIF
        CALL RmDblSpc(message_text)
        CALL local_message('',TRIM(message_text))
        ActGroup%flags = IOR(ActGroup%flags,INPUT_GROUP_MSG_DONE) ! Prevent sucessive cycles from producing output
      ENDIF
    ENDIF

  END SUBROUTINE GroupReadMsg

! Get the actual buffer in a variable's buffer list and progress buffer pointer
  FUNCTION InputVarGetNextBuffer(var,buf,lProg)
  TYPE (input_data_list), POINTER :: InputVarGetNextBuffer
  TYPE (input_var_list),  POINTER :: var
  TYPE (input_data_list), POINTER :: buf
  LOGICAL, INTENT(in)             :: lProg

    IF (lProg) buf => buf%next
    IF (ASSOCIATED(buf%dta)) THEN
      InputVarGetNextBuffer => buf
    ELSE
      InputVarGetNextBuffer => var%Intrp
    ENDIF

  END FUNCTION InputVarGetNextBuffer

! Get the next mask in a variable's mask list and progress mask pointer
  FUNCTION InputVarGetNextMask(var,mask,lProg)
  TYPE (input_mask_list), POINTER :: InputVarGetNextMask
  TYPE (input_var_list),  POINTER :: var
  TYPE (input_mask_list), POINTER :: mask
  LOGICAL, INTENT(in)            :: lProg

    InputVarGetNextMask => mask
    IF (lProg) mask => mask%next
    RETURN
    var=>var ! Eliminate compiler warnings

  END FUNCTION InputVarGetNextMask

! Dummy function doing nothing to the data
  SUBROUTINE DataProcNoAction(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),            POINTER :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = SIZE(var%saction) + 1 ! Skip over all actions
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'NoAction'

  END SUBROUTINE DataProcNoAction

! Add debug message for non-file variables
  SUBROUTINE DataProcMessage(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Msg'
    IF ((IAND(level,INPUT_PROC_DATA_ACTION) /= 0) .AND. (IAND(mo_debug,INPUT_DBG_MASK) >= INPUT_DBG_READ_VAR) .AND. ldbg) THEN
      IF (ASSOCIATED(var%group) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) THEN
        var%group%time     = var%time
        var%group%curr_rec = 0
      ENDIF
      IF (ASSOCIATED(CurrGroup) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) CALL GroupReadMsg(CurrGroup)
      CurrGroup => var%group
      WRITE(message_text,'(5a)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Updt; ',TRIM(var%name_var)
    ENDIF

  END SUBROUTINE DataProcMessage

! Read data from single file
  SUBROUTINE DataProcRead(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data
  INTEGER :: nDim, cp

    IF (IAND(level,INPUT_PROC_PROGRESS) /= 0) THEN
      cp   = curr%acp
      nDim = var%saction(cp+2)
      curr%acp = curr%acp + 3 + 2*nDim
    ENDIF
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Read'
    IF (IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      IF (ASSOCIATED(var%group) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) THEN
        var%group%time     = var%time
        IF (var%group%curr_rec >= 0) THEN
          var%group%curr_rec = var%parent_file%curr_rec
        ELSE
          var%group%curr_rec = var%group%curr_rec - 1
        ENDIF
      ENDIF
      IF (ASSOCIATED(CurrGroup) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) CALL GroupReadMsg(CurrGroup)
      CurrGroup => var%group
      act_data => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
#ifdef HAVE_F2003
      IF (.NOT. var%parent_file%ft%FcnRead(var%parent_file,var,var%saction(cp+1),var%saction(cp+2),var%saction(cp+3:cp+2+2*nDim), &
                          1,var%parent_file%curr_rec,act_data%dta)) CALL local_error('DataProcRead',TRIM(message_text(2:)))
#else
      IF (.NOT. NCFileRead(var%parent_file,var,var%saction(cp+1),var%saction(cp+2),var%saction(cp+3:cp+2+2*nDim), &
                          1,var%parent_file%curr_rec,act_data%dta)) CALL local_error('DataProcRead',TRIM(message_text(2:)))
#endif
      IF (( (IAND(mo_debug,INPUT_DBG_MASK) >= INPUT_DBG_READ_VAR) .OR. &
           ((IAND(mo_debug,INPUT_DBG_READ_GROUP) /= 0) .AND. .NOT. ASSOCIATED(var%group))) .AND. ldbg) THEN
        message_text = TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK)))//': Read; '//TRIM(var%name_var)//' ('// &
                       TRIM(var%parent_file%name_file)
        IF (IAND(var%stat,INPUT_VAR_HAS_UNLIM) /= 0) &
          WRITE(message_text(LEN_TRIM(message_text)+1:),'(a,i5)') ',',var%parent_file%curr_rec
        cp = LEN_TRIM(message_text)+1
        message_text(cp:cp) = ')'
        CALL RmDblSpc(message_text)
      ENDIF
      var%parent_file%nRead_dbg = var%parent_file%nRead_dbg + 1
    ENDIF

  END SUBROUTINE DataProcRead

! Read data from multiple files 
  SUBROUTINE DataProcReadMulti(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_file_list), POINTER :: curr_file
  TYPE (input_data_list), POINTER :: act_data
  CHARACTER (len=256) :: tmp_txt
  INTEGER :: i, nDim, ofs, mul, tmp, cp

    IF (IAND(level,INPUT_PROC_PROGRESS) /= 0) THEN
      cp   = curr%acp
      tmp  = var%saction(cp+1)
      nDim = var%saction(cp+tmp+4)
      curr%acp = curr%acp + 5 + 2*nDim + tmp
    ENDIF
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Multi_read'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      IF (ASSOCIATED(var%group) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) THEN
        var%group%time     = var%time
        IF (var%group%curr_rec >= 0) THEN
          var%group%curr_rec = var%parent_file%curr_rec
        ELSE
          var%group%curr_rec = var%group%curr_rec - 1
        ENDIF
      ENDIF
      IF (ASSOCIATED(CurrGroup) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) CALL GroupReadMsg(CurrGroup)
      CurrGroup => var%group
      act_data => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      tmp_txt = ''
      curr_file => var%parent_file
      ofs  = cp + tmp + nDim + 4
      mul  = PRODUCT(var%saction(ofs+1:ofs+nDim-1)) ! Last dimension is artificial, exclude from product
      ofs  = 1
      DO i = 1,tmp
#ifdef HAVE_F2003
        IF (.NOT. curr_file%ft%FcnRead(curr_file,var,var%saction(cp+1+i),var%saction(cp+tmp+4)-1, &
                                       (/var%saction(cp+tmp     +5:cp+tmp+  nDim+3),   & ! Last dim is artificial, rm from read
                                         var%saction(cp+tmp+nDim+5:cp+tmp+2*nDim+3)/), &
                                     ofs,var%parent_file%curr_rec,act_data%dta))       &
#else
        IF (.NOT. NCFileRead(curr_file,var,var%saction(cp+1+i),var%saction(cp+tmp+4)-1, &
                                        (/var%saction(cp+tmp     +5:cp+tmp+  nDim+3),   & ! Last dim is artificial, rm from read
                                          var%saction(cp+tmp+nDim+5:cp+tmp+2*nDim+3)/), &
                                      ofs,var%parent_file%curr_rec,act_data%dta))       &
#endif
          CALL local_error('DataProcReadMulti',message_text)
        ofs = ofs + mul
        IF (IAND(var%stat,INPUT_VAR_MULTI_SPREAD) /= 0) THEN
          IF ((IAND(mo_debug,INPUT_DBG_MASK)>=INPUT_DBG_READ_VAR) .AND. ldbg) &
            tmp_txt = TRIM(tmp_txt)//','//TRIM(curr_file%name_file)
          curr_file => curr_file%next
        ENDIF
      ENDDO
      IF (( (IAND(mo_debug,INPUT_DBG_MASK) >= INPUT_DBG_READ_VAR) .OR. &
           ((IAND(mo_debug,INPUT_DBG_READ_GROUP) /= 0) .AND. .NOT. ASSOCIATED(var%group))) .AND. ldbg) THEN
        IF (LEN_TRIM(tmp_txt) == 0) tmp_txt = ','//TRIM(var%parent_file%name_file)
        message_text = TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK)))//': Read; '//TRIM(var%name_var)//' ('// &
                       TRIM(tmp_txt(2:))
        IF (IAND(var%stat,INPUT_VAR_HAS_UNLIM) /= 0) &
          WRITE(message_text(LEN_TRIM(message_text)+1:),'(a,i5)') ',',var%parent_file%curr_rec
        cp = LEN_TRIM(message_text)+1
        message_text(cp:cp) = ')'
        CALL RmDblSpc(message_text)
      ENDIF
      var%parent_file%nRead_dbg = var%parent_file%nRead_dbg + 1
    ENDIF

  END SUBROUTINE DataProcReadMulti

! Permute array dimensions
  SUBROUTINE DataProcPermute(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  INTEGER :: nDim, cp

    IF (IAND(level,INPUT_PROC_PROGRESS) /= 0) THEN
      cp = curr%acp
      nDim = var%saction(cp+1)
      curr%acp = curr%acp + 2 + 2*nDim
    ENDIF
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Permute'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      IF (.NOT. ArrayPermute(act_data%dta,next_data%dta,var%saction(cp+2:cp+1+nDim),var%saction(cp+2+nDim:cp+1+2*nDim))) &
        CALL local_error('DataProcPermute','Permuting error')
    ENDIF

  END SUBROUTINE DataProcPermute

! Reverse array data across a dimension
  SUBROUTINE DataProcReverse(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  INTEGER :: nDim, cp

    IF (IAND(level,INPUT_PROC_PROGRESS) /= 0) THEN
      cp = curr%acp
      nDim = var%saction(cp+1)
      curr%acp = curr%acp + 3 + nDim
    ENDIF
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Reverse'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      CALL ArrayReverse(act_data%dta,var%saction(cp+2:cp+nDim+1),var%saction(cp+nDim+2),next_data%dta)
    ENDIF

  END SUBROUTINE DataProcReverse

! Swap two data blocks
  SUBROUTINE DataProcSwap(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  INTEGER :: cp

    IF (IAND(level,INPUT_PROC_PROGRESS) /= 0) THEN
      cp = curr%acp
      curr%acp = curr%acp + 4
    ENDIF
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Swap'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      CALL ArraySwap(act_data%dta,var%saction(cp+1),var%saction(cp+2),var%saction(cp+3),next_data%dta)
    ENDIF

  END SUBROUTINE DataProcSwap

! Pack data according to mask
  SUBROUTINE DataProcPack(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  TYPE (input_mask_list), POINTER :: act_mask

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 3
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Pack'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      act_mask  => InputVarGetNextMask(var,curr%mask,.TRUE.)
      CALL MaskPack(act_data%dta,next_data%dta,act_mask%mask,var%saction(curr%acp-2),var%saction(curr%acp-1))
    ENDIF

  END SUBROUTINE DataProcPack

! Unpack data according to mask
  SUBROUTINE DataProcUnpack(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  TYPE (input_mask_list), POINTER :: act_mask
  REAL(dp), POINTER :: def(:,:,:,:,:,:,:)
  LOGICAL :: lDone

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 3
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Unpack'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      act_mask  => InputVarGetNextMask(var,curr%mask,.TRUE.)
      def => InputDataGetSpecified(var%aux_data,INPUT_DATA_DEFAULT,0)
      lDone = .FALSE.
      IF (ASSOCIATED(def)) THEN
        IF (SIZE(def) == 1) THEN
          CALL MaskUnpack(act_data%dta,next_data%dta,act_mask%mask,var%saction(curr%acp-2),var%saction(curr%acp-1), &
                          def(1,1,1,1,1,1,1))
          lDone = .TRUE.
        ENDIF
      ENDIF
      IF (.NOT. lDone) &
        CALL MaskUnpack(act_data%dta,next_data%dta,act_mask%mask,var%saction(curr%acp-2),var%saction(curr%acp-1))
    ENDIF

  END SUBROUTINE DataProcUnpack

! Map to equivalent dimension (actually the same as unpack at this stage)
  SUBROUTINE SetData(buf,st,sz,val)
  REAL(dp), INTENT(inout) :: buf(*)
  INTEGER,  INTENT(in)    :: st,sz
  REAL(dp), INTENT(in)    :: val

    buf(st:st+sz-1) = val

  END SUBROUTINE SetData

  SUBROUTINE DataProcEquiv(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  TYPE (input_mask_list), POINTER :: act_mask
  INTEGER rem, sz

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 3
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Equivalence'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_mask  => InputVarGetNextMask(var,curr%mask,.TRUE.)
      IF (var%saction(curr%acp-3) < 0) THEN ! Nothing to do, except perhaps to fill non-defined areas of the buffer 
        rem = act_mask%mask(SIZE(act_mask%mask)) ! Dirty direct way to get the number of values to fill
        IF (rem > 0) THEN
          act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
          sz = SIZE(act_data%dta)
          rem = rem * var%saction(curr%acp-2)
          CALL SetData(act_data%dta,sz-rem+1,rem,0._dp)
        ENDIF
        RETURN ! Nothing to do - the equivalence has in reality been performed by the last action
      ENDIF
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      CALL MaskUnpack(act_data%dta,next_data%dta,act_mask%mask,var%saction(curr%acp-2),var%saction(curr%acp-1),0._dp)
    ENDIF

  END SUBROUTINE DataProcEquiv

! Copy to an array of identical shape
  SUBROUTINE DataProcCopy(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Copy'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      IF (var%saction(curr%acp-1) < 0) RETURN ! Nothing to do - the copy has in reality been performed by the last action
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      next_data%dta(:,:,:,:,:,:,:) = act_data%dta(:,:,:,:,:,:,:)
    ENDIF

  END SUBROUTINE DataProcCopy

! Weights individual entries along a dimension with individual weights
  SUBROUTINE DataProcWeight_int(dta,weights,dimen,sz)
  REAL(dp), POINTER    :: dta(:,:,:,:,:,:,:)
  REAL(dp), INTENT(in) ::weights(:)
  INTEGER, INTENT(in)  :: dimen, sz

  INTEGER :: i

    IF (SIZE(weights) /= sz) CALL local_error('DataProcWeight_int','Internal error - weight-dimension size mismatch')
    DO i=1,sz
      SELECT CASE (dimen)
        CASE (1)
          dta(i,:,:,:,:,:,:) = dta(i,:,:,:,:,:,:) * weights(i)
        CASE (2)
          dta(:,i,:,:,:,:,:) = dta(:,i,:,:,:,:,:) * weights(i)
        CASE (3)
          dta(:,:,i,:,:,:,:) = dta(:,:,i,:,:,:,:) * weights(i)
        CASE (4)
          dta(:,:,:,i,:,:,:) = dta(:,:,:,i,:,:,:) * weights(i)
        CASE (5)
          dta(:,:,:,:,i,:,:) = dta(:,:,:,:,i,:,:) * weights(i)
        CASE (6)
          dta(:,:,:,:,:,i,:) = dta(:,:,:,:,:,i,:) * weights(i)
        CASE (7)
          dta(:,:,:,:,:,:,i) = dta(:,:,:,:,:,:,i) * weights(i)
      END SELECT
    ENDDO

  END SUBROUTINE DataProcWeight_int

! Sum data over an array dimension
  SUBROUTINE DataProcSum(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  REAL(dp), POINTER :: weights(:,:,:,:,:,:,:)

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 3
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Sum'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      weights   => InputDataGetSpecified(var%aux_data,INPUT_DATA_WEIGHTS,0)
      IF (ASSOCIATED(weights)) & ! Weight data with provided weights
        CALL DataProcWeight_int(act_data%dta,weights(:,1,1,1,1,1,1),var%saction(curr%acp-2),var%saction(curr%acp-1))
      next_data%dta(:,:,:,:,:,:,1) = SUM(act_data%dta,DIM=var%saction(curr%acp-2))
    ENDIF

  END SUBROUTINE DataProcSum


! Average data over an array dimension
  SUBROUTINE DataProcAvg(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  REAL(dp), POINTER :: weights(:,:,:,:,:,:,:)

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 3
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Avg'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      weights   => InputDataGetSpecified(var%aux_data,INPUT_DATA_WEIGHTS,0)
      IF (ASSOCIATED(weights)) & ! Weight data with provided weights
        CALL DataProcWeight_int(act_data%dta,weights(:,1,1,1,1,1,1),var%saction(curr%acp-2),var%saction(curr%acp-1))
      next_data%dta(:,:,:,:,:,:,1) = SUM(act_data%dta,DIM=var%saction(curr%acp-2)) / REAL(var%saction(curr%acp-1),dp)
    ENDIF

  END SUBROUTINE DataProcAvg

! Normalize eventually weighted data over an array dimension
  SUBROUTINE DataProcNormalize(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  REAL(dp), POINTER :: weights(:,:,:,:,:,:,:)
  REAL (dp) :: rmul

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 3
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Normalize'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      weights   => InputDataGetSpecified(var%aux_data,INPUT_DATA_WEIGHTS,0)
      IF (ASSOCIATED(weights)) THEN ! Weight data with provided weights
        CALL DataProcWeight_int(act_data%dta,weights(:,1,1,1,1,1,1),var%saction(curr%acp-2),var%saction(curr%acp-1))
        rmul = 1._dp / SUM(weights)
      ELSE
        rmul = 1._dp / REAL(var%saction(curr%acp-1),dp)
      ENDIF
      next_data%dta(:,:,:,:,:,:,1) = SUM(act_data%dta,DIM=var%saction(curr%acp-2)) * rmul
    ENDIF

  END SUBROUTINE DataProcNormalize

! Replace bad values from file with specified value
  SUBROUTINE DataProcFillBad(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data
  REAL(dp), POINTER :: fill_val(:,:,:,:,:,:,:)
  CHARACTER (len=256) :: tmpstr

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Fill_bad'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      fill_val  => InputDataGetSpecified(var%aux_data,INPUT_DATA_FILL,0)
      IF (ASSOCIATED(fill_val)) THEN
        IF ((IAND(var%flags,INPUT_VAR_BAD_VAL_WARN)/=0) .OR. (IAND(var%flags,INPUT_VAR_BAD_VAL_MASK)==INPUT_VAR_BAD_VAL_ERR)) THEN
          IF (ANY(act_data%dta(:,:,:,:,:,:,:) == fill_val(2,1,1,1,1,1,1))) THEN
            SELECT CASE (IAND(var%flags,INPUT_VAR_BAD_VAL_MASK))
              CASE (INPUT_VAR_BAD_VAL_ERR)
                CALL local_error('DataProcFillBad',TRIM(var%name_var)//' contains invalid fill values')
              CASE (INPUT_VAR_BAD_VAL_WARN, INPUT_VAR_BAD_VAL_WARN_REPLACE)
                tmpstr = message_text
                IF (IAND(var%flags,INPUT_VAR_BAD_VAL_REPLACE) /= 0) THEN
                  CALL local_message('WARNING',TRIM(var%name_var)//' contained fill values which have been replaced')
                ELSE
                  CALL local_message('WARNING',TRIM(var%name_var)//' contains (bad) fill values')
                ENDIF
                message_text = tmpstr
            END SELECT
            var%flags = IAND(var%flags,INPUT_ALLFLAGS - INPUT_VAR_BAD_VAL_WARN) + INPUT_VAR_BAD_VAL_IGNORE
          ENDIF
        ENDIF
        IF (IAND(var%flags,INPUT_VAR_BAD_VAL_REPLACE) /= 0) THEN ! Utilizes that WARN_REPLACE and REPLACE have same bit set
          WHERE (act_data%dta(:,:,:,:,:,:,:) == fill_val(2,1,1,1,1,1,1))
            act_data%dta(:,:,:,:,:,:,:)       = fill_val(1,1,1,1,1,1,1)
          ENDWHERE
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE DataProcFillBad

! Replace out-of-valid-range values from file with specified value
  SUBROUTINE DataProcFillOutOfRange(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data
  REAL(dp), POINTER :: fill_val(:,:,:,:,:,:,:)
  CHARACTER (len=256) :: tmpstr

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Fill_range'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      fill_val  => InputDataGetSpecified(var%aux_data,INPUT_DATA_RANGE,0)
      IF (ASSOCIATED(fill_val)) THEN
        IF ((IAND(var%flags,INPUT_VAR_BAD_VAL_WARN)/=0) .OR. (IAND(var%flags,INPUT_VAR_BAD_VAL_MASK)==INPUT_VAR_BAD_VAL_ERR)) THEN
          IF (ANY(     (act_data%dta(:,:,:,:,:,:,:) < fill_val(2,1,1,1,1,1,1))        &
                  .OR. (act_data%dta(:,:,:,:,:,:,:) > fill_val(3,1,1,1,1,1,1)))) THEN
            SELECT CASE (IAND(var%flags,INPUT_VAR_BAD_VAL_MASK))
              CASE (INPUT_VAR_BAD_VAL_ERR)
                CALL local_error('DataProcFillBad',TRIM(var%name_var)//' contains out-of-range values')
              CASE (INPUT_VAR_BAD_VAL_WARN, INPUT_VAR_BAD_VAL_WARN_REPLACE)
                tmpstr = message_text
                IF (IAND(var%flags,INPUT_VAR_BAD_VAL_REPLACE) /= 0) THEN
                  CALL local_message('WARNING',TRIM(var%name_var)//' contained out-of-range values which have been replaced')
                ELSE
                  CALL local_message('WARNING',TRIM(var%name_var)//' contains out-of-range values')
                ENDIF
                message_text = tmpstr
            END SELECT
            var%flags = IAND(var%flags,INPUT_ALLFLAGS - INPUT_VAR_BAD_VAL_WARN) + INPUT_VAR_BAD_VAL_IGNORE
          ENDIF
        ENDIF
        IF (IAND(var%flags,INPUT_VAR_BAD_VAL_REPLACE) /= 0) THEN ! Utilizes that WARN_REPLACE and REPLACE have same bit set
          WHERE ((act_data%dta(:,:,:,:,:,:,:) < fill_val(2,1,1,1,1,1,1)) .OR. &
                 (act_data%dta(:,:,:,:,:,:,:) > fill_val(3,1,1,1,1,1,1)))
            act_data%dta(:,:,:,:,:,:,:)       = fill_val(1,1,1,1,1,1,1)
          ENDWHERE
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE DataProcFillOutOfRange

! Rescale data according to multiplicative scale and additional constant 
  SUBROUTINE DataProcRescale(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data
  REAL(dp), POINTER :: rescale(:,:,:,:,:,:,:)

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 2
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Rescale'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      rescale   => InputDataGetSpecified(var%aux_data,INPUT_DATA_RESCALE,var%saction(curr%acp-1))
      IF (ASSOCIATED(rescale)) &
        act_data%dta(:,:,:,:,:,:,:) = act_data%dta(:,:,:,:,:,:,:) * rescale(1,1,1,1,1,1,1) + rescale(2,1,1,1,1,1,1)
    ENDIF

  END SUBROUTINE DataProcRescale

! Rescale data according to multiplicative scale and additional constant 
  SUBROUTINE DataProcTimeDiv(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data
  REAL (dp) :: mul

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'TimeDiv'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      IF (IAND(var%flags,INPUT_VAR_DIV_YEAR) /= 0) THEN
        mul = 1._dp / REAL(CalendarYearLen(model_time(1)),dp)
      ELSE
        mul = 1._dp / REAL(CalendarMonthLen(model_time(1),model_time(2)),dp)
      ENDIF
      act_data%dta(:,:,:,:,:,:,:) = act_data%dta(:,:,:,:,:,:,:) * mul
    ENDIF

  END SUBROUTINE DataProcTimeDiv

! Spread data over one or more array dimension(s)
  SUBROUTINE DataProcSpread(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  INTEGER :: nDim, tmp, cp

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) THEN
      cp   = curr%acp
      nDim = var%saction(cp+1)
      tmp  = var%saction(cp+nDim+2)
      curr%acp  = curr%acp + 3 + nDim + 2*tmp
    ENDIF
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Spread'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      IF (.NOT. ArrayCopy(act_data%dta,next_data%dta,var%saction(cp+2:cp+1+nDim), &
                          CopyDims=var%saction(cp+3+nDim    :cp+2+nDim+  tmp),    &
                          nCopy   =var%saction(cp+3+nDim+tmp:cp+2+nDim+2*tmp)))   &
        CALL local_error('DataProcSpread','Spreading error')
    ENDIF

  END SUBROUTINE DataProcSpread

! Setting only a single dimension of destination data (TODO: Untested)
  SUBROUTINE DataProcSubset(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 3
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Subset'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      SELECT CASE (var%saction(curr%acp-2))
        CASE (1)
          next_data%dta(var%saction(curr%acp-1),:,:,:,:,:,:) = act_data%dta(:,:,:,:,:,:,1)
        CASE (2)
          next_data%dta(:,var%saction(curr%acp-1),:,:,:,:,:) = act_data%dta(:,:,:,:,:,:,1)
        CASE (3)
          next_data%dta(:,:,var%saction(curr%acp-1),:,:,:,:) = act_data%dta(:,:,:,:,:,:,1)
        CASE (4)
          next_data%dta(:,:,:,var%saction(curr%acp-1),:,:,:) = act_data%dta(:,:,:,:,:,:,1)
        CASE (5)
          next_data%dta(:,:,:,:,var%saction(curr%acp-1),:,:) = act_data%dta(:,:,:,:,:,:,1)
        CASE (6)
          next_data%dta(:,:,:,:,:,var%saction(curr%acp-1),:) = act_data%dta(:,:,:,:,:,:,1)
        CASE (7)
          next_data%dta(:,:,:,:,:,:,var%saction(curr%acp-1)) = act_data%dta(:,:,:,:,:,:,1)
      END SELECT
    ENDIF

  END SUBROUTINE DataProcSubset

! Interpolate data in time
  SUBROUTINE DataProcInterpolateTimeLinear(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: next_data
  INTEGER (i8) :: dt, st
  REAL (dp) :: frac

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Time_interpolate_linear'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      IF (var%intrp%time(1) == -HUGE(i8)) THEN
        IF (var%intrp%prev%time(1) == -HUGE(i8)) CALL local_error('DataProcInterpolateTimeLinear', &
                                                   'No valid data to interpolate for variable '//TRIM(var%name_var))
        next_data%dta(:,:,:,:,:,:,:) = var%intrp%prev%dta(:,:,:,:,:,:,:)
        RETURN
      ENDIF
      IF (var%intrp%prev%time(1) == -HUGE(i8)) THEN
        next_data%dta(:,:,:,:,:,:,:) = var%intrp%dta(:,:,:,:,:,:,:)
        RETURN
      ENDIF
      dt = CalendarTimeDiff(var%intrp%prev%time,var%intrp%time)
      st = CalendarTimeDiff(var%intrp%prev%time,model_time)
      frac = REAL(st,dp)/REAL(dt,dp)
      IF (ASSOCIATED(var%aux_data)) THEN
        IF (IAND(var%aux_data%flags,INPUT_DATA_TYPE_MASK) == INPUT_DATA_FUTURE_FRAC) & ! future_frac are always first in list
          st = CalendarTimeDiff(var%intrp%prev%time,var%aux_data%time)
      ENDIF
      frac = REAL(st,dp)/REAL(dt,dp)
      next_data%dta(:,:,:,:,:,:,:) = var%intrp%dta     (:,:,:,:,:,:,:) *          frac  &
                                   + var%intrp%prev%dta(:,:,:,:,:,:,:) * (1._dp - frac)
    ENDIF

  END SUBROUTINE DataProcInterpolateTimeLinear

! Broadcast data to all destination PE's
  SUBROUTINE DataProcBroadcast(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) THEN
      IF (liodata) THEN
        msg = 'Broadcast'
      ELSE
        msg = 'Recieve'
      ENDIF
    ENDIF
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      IF (ASSOCIATED(var%group) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) THEN
        var%group%time     = var%time
        var%group%curr_rec = var%parent_file%curr_rec
      ENDIF
      IF (ASSOCIATED(CurrGroup) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) CALL GroupReadMsg(CurrGroup)
      CurrGroup => var%group
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      CALL InputVarBroadcast(act_data)
#endif
      IF ((IAND(mo_debug,INPUT_DBG_MASK) >= INPUT_DBG_READ_VAR) .AND. ldbg .AND. .NOT. liodata) &
        WRITE(message_text,'(3a)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Recv; ',TRIM(var%name_var)
    ENDIF

  END SUBROUTINE DataProcBroadcast

! Scatter data to all destination PE's
  SUBROUTINE DataProcScatter(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data, src_data, dst_data, buf_data
  TYPE (input_mask_list), POINTER :: act_mask
  REAL (dp),              POINTER :: rtmp(:,:,:,:,:,:,:)

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) THEN
      IF (liodata) THEN
        msg = 'Scatter'
      ELSE
        msg = 'Recieve'
      ENDIF
    ENDIF
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      IF (ASSOCIATED(var%group) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) THEN
        var%group%time     = var%time
        var%group%curr_rec = var%parent_file%curr_rec
      ENDIF
      IF (ASSOCIATED(CurrGroup) .AND. .NOT. ASSOCIATED(CurrGroup,var%group)) CALL GroupReadMsg(CurrGroup)
      CurrGroup => var%group
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      NULLIFY(src_data,buf_data,next_data)
      IF (liodata) THEN
        next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
        src_data => act_data
        rtmp => act_data%dta
        IF (IAND(var%stat,INPUT_VAR_CONT_SEND) /= 0) THEN
          dst_data => next_data
        ELSE
          buf_data => next_data
          dst_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
        ENDIF
      ELSE
        dst_data => act_data
        rtmp => act_data%dta ! Not used, but it seems to be illegal to transfere a nullified pointer to a dimension(*)
      ENDIF
      act_mask  => InputVarGetNextMask(var,curr%mask,.FALSE.) ! InputVarScatter takes care of progessing the mask pointer
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      curr%mask => InputVarScatter(act_mask,src_data,rtmp,dst_data,buf_data)
#endif
      IF ((IAND(mo_debug,INPUT_DBG_MASK) >= INPUT_DBG_READ_VAR) .AND. ldbg .AND. .NOT. liodata) &
        WRITE(message_text,'(3a)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Recv; ',TRIM(var%name_var)
    ENDIF

  END SUBROUTINE DataProcScatter

! Adding src and dst data to dst
  SUBROUTINE DataProcAdd(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Add'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      next_data%dta(:,:,:,:,:,:,:) = next_data%dta(:,:,:,:,:,:,:) + act_data%dta(:,:,:,:,:,:,:)
    ENDIF

  END SUBROUTINE DataProcAdd

! Subtract src from dst data to dst
  SUBROUTINE DataProcSub(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Sub'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      next_data%dta(:,:,:,:,:,:,:) = next_data%dta(:,:,:,:,:,:,:) - act_data%dta(:,:,:,:,:,:,:)
    ENDIF

  END SUBROUTINE DataProcSub

! Subtracting dst from src data to dst
  SUBROUTINE DataProcSubr(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Subr'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      next_data%dta(:,:,:,:,:,:,:) = act_data%dta(:,:,:,:,:,:,:) - next_data%dta(:,:,:,:,:,:,:)
    ENDIF

  END SUBROUTINE DataProcSubr

! Multiplying src and dst data to dst
  SUBROUTINE DataProcMul(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Mul'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      next_data%dta(:,:,:,:,:,:,:) = next_data%dta(:,:,:,:,:,:,:) * act_data%dta(:,:,:,:,:,:,:)
    ENDIF

  END SUBROUTINE DataProcMul

! Dividing dst with src data to dst
  SUBROUTINE DataProcDiv(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Div'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      next_data%dta(:,:,:,:,:,:,:) = next_data%dta(:,:,:,:,:,:,:) / act_data%dta(:,:,:,:,:,:,:)
    ENDIF

  END SUBROUTINE DataProcDiv

! Dividing src with dst data to dst
  SUBROUTINE DataProcDivr(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Divr'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      next_data%dta(:,:,:,:,:,:,:) = act_data%dta(:,:,:,:,:,:,:) / next_data%dta(:,:,:,:,:,:,:)
    ENDIF

  END SUBROUTINE DataProcDivr

! Adding src flux and dst data to dst
  SUBROUTINE DataProcFlux(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Flux'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      next_data%dta(:,:,:,:,:,:,:) = next_data%dta(:,:,:,:,:,:,:) + act_data%dta(:,:,:,:,:,:,:) * model_dt
    ENDIF

  END SUBROUTINE DataProcFlux

! Nudge src and dst data to dst
  SUBROUTINE DataProcNudge(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  REAL(dp), POINTER :: weight(:,:,:,:,:,:,:)

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'Nudge'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      weight => InputDataGetSpecified(var%aux_data,INPUT_DATA_NUDGE_SCA,0)
      IF (ASSOCIATED(weight)) &
        next_data%dta(:,:,:,:,:,:,:) = next_data%dta(:,:,:,:,:,:,:) * (1._dp - weight(1,1,1,1,1,1,1)) &
                                     +  act_data%dta(:,:,:,:,:,:,:) *          weight(1,1,1,1,1,1,1)
    ENDIF

  END SUBROUTINE DataProcNudge

! Nudge src and dst data to dst with individual weights
  SUBROUTINE DataProcNudgeFld_private(modd,readd,weights,len) ! Extra subroutine necessary to get rid of array shapes
  REAL(dp), DIMENSION(*), INTENT(inout) :: modd
  REAL(dp), DIMENSION(*), INTENT(in)    :: readd, weights
  INTEGER,                INTENT(in)    :: len

    modd(1:len) = modd(1:len) * (1._dp - weights(1:len)) + readd(1:len) * weights(1:len)

  END SUBROUTINE DataProcNudgeFld_private

  SUBROUTINE DataProcNudgeFld(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),        POINTER     :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

  TYPE (input_data_list), POINTER :: act_data, next_data
  REAL(dp), POINTER :: weight(:,:,:,:,:,:,:)

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = curr%acp + 1
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'NudgeFld'
    IF ( IAND(level,INPUT_PROC_DATA_ACTION) /= 0) THEN
      act_data  => InputVarGetNextBuffer(var,curr%buf,.FALSE.)
      next_data => InputVarGetNextBuffer(var,curr%buf,.TRUE.)
      weight => InputDataGetSpecified(var%aux_data,INPUT_DATA_NUDGE_FLD,0)
      IF (.NOT. ASSOCIATED(weight)) CALL local_error('DataProcNudgeFld','No weighting data found: '//TRIM(var%name_var))
      IF (SIZE(next_data%dta) /= SIZE(weight)) &
        CALL local_error('DataProcNudgeFld','Field size mismatch between field and weights ('//TRIM(var%name_var)//')')
      CALL DataProcNudgeFld_private(next_data%dta,act_data%dta,weight,SIZE(next_data%dta))
    ENDIF

  END SUBROUTINE DataProcNudgeFld

#ifndef HAVE_F2003
  SUBROUTINE InputDoDataActions(var,curr,lterm)
  TYPE (input_var_list),  POINTER :: var
  TYPE (input_curr),INTENT(inout) :: curr
  LOGICAL,          INTENT(in   ) :: lterm

  LOGICAL :: lCont

    lCont = ASSOCIATED(var%saction)
    DO WHILE (lCont)
      SELECT CASE (ABS(var%saction(curr%acp)))
        CASE (INPUT_VAR_DO_NOACTION);          CALL DataProcNoAction             (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_TOREAD);            CALL DataProcRead                 (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_MULTI_READ);        CALL DataProcReadMulti            (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_PERMUTE);           CALL DataProcPermute              (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_PACK);              CALL DataProcPack                 (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_UNPACK);            CALL DataProcUnpack               (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_EQUIV);             CALL DataProcEquiv                (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COPY);              CALL DataProcCopy                 (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_REVERSE);           CALL DataProcReverse              (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_SWAP);              CALL DataProcSwap                 (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_SUM);               CALL DataProcSum                  (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_AVG);               CALL DataProcAvg                  (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_NORM);              CALL DataProcNormalize            (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_FILL_BAD);          CALL DataProcFillBad              (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_FILL_RANGE);        CALL DataProcFillOutofRange       (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_RESCALE);           CALL DataProcRescale              (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_SPREAD);            CALL DataProcSpread               (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_BROADCAST);         CALL DataProcBroadcast            (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_SCATTER);           CALL DataProcScatter              (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_SUBSET);            CALL DataProcSubSet               (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_INTERPOLATE);       CALL DataProcInterpolateTimeLinear(INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COMBINE_ADD);       CALL DataProcAdd                  (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COMBINE_SUB);       CALL DataProcSub                  (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COMBINE_SUBR);      CALL DataProcSubr                 (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COMBINE_MUL);       CALL DataProcMul                  (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COMBINE_DIV);       CALL DataProcDiv                  (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COMBINE_DIVR);      CALL DataProcDivr                 (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COMBINE_ADD_FLUX);  CALL DataProcFlux                 (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COMBINE_NUDGE);     CALL DataProcNudge                (INPUT_PROC_DO_DATA_ACTION,var,curr)
        CASE (INPUT_VAR_DO_COMBINE_NUDGE_FLD); CALL DataProcNudgeFld             (INPUT_PROC_DO_DATA_ACTION,var,curr)
      END SELECT
      lCont = (curr%acp <= SIZE(var%saction))
      IF (lCont .AND. lterm) lCont = (var%saction(curr%acp) /= INPUT_VAR_DO_INTERPOLATE)
    ENDDO

  END SUBROUTINE InputDoDataActions
#endif

END MODULE mo_input_proc
