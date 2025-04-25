!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! HISTORY:
! V.1.00: 06-2011 - Stiig Wilkenskjeld, MPI-M - first attempt
! V.2.00: 04-2012 - Stiig Wilkenskjeld, MPI-M - Support for absolute time, user interpolation, improved interface and more... 
! V.2.10: 09-2012 - Stiig Wilkenskjeld, MPI-M - Minor modifications of time control
! V.2.20: 09-2012 - Stiig Wilkenskjeld, MPI-M - Support for Julian calendar
! V.3.00: 09-2013 - Stiig Wilkenskjeld, MPI-M - New internal data structures, <=7D vars., timestep unit,
!                                               interface support for multiple calendars, separating NetCDF-lib calls, 
!                                               namelist control, different ways of combining model and read data, 
!                                               chuncked reading, handling of mismatching dimensions
! V.3.09: 12-2013 - Stiig Wilkenskjeld, MPI-M - Equivalent dimensions, bug fixes
! V.3.10: 01-2014 - Stiig Wilkenskjeld, MPI-M - Joined distributed dimensions, restructured real type auxillary data,
!                                               removal of unnecessary data copying, prep of support for multiple file formats 
! V.3.20: 01-2014 - Stiig Wilkenskjeld, MPI-M - Exchanged concept "initial file" for "initial variable", cyclic dimensions
! V.3.21: 04-2014 - Stiig Wilkenskjeld, MPI-M - Bug fixes, mainly rewriting InputSetVarStartTime
! V.3.30: 04-2014 - Stiig Wilkenskjeld, MPI-M - Variable groups, new valid_time and debug options, extended error treatment
! V.3.40: 07-2014 - Stiig Wilkenskjeld, MPI-M - Support for handling NetCDF-fill-values and flexible handling of initial variables
!                                               with time variable. Extensions of debug routines.
! V.3.50: 08-2014 - Stiig Wilkenskjeld, MPI-M - Selectable treatment of fill/out-of-range data, dependent variables, inquire fcns,
!                                               set dt of var to dt of file.
! V.3.60: 09-2014 - Stiig Wilkenskjeld, MPI-M - Change of file structure, variables witout model data, returning 4D-var-data, 
!                                               group default "parameters", scalar variables, reading of meta-data, 'model' vars.
!
! Planned changes: 
!  - Chunked reading should also work together with equivalent dimensions, collective dim dist and data swapped cyclic dimensions.
!  - Build a memory mechanism for mapping combinations of variable and file dims to a process chain for performance reasons.

MODULE mo_input

USE mo_input_strings
USE mo_input_calendar
USE mo_input_arrays
USE mo_input_types
USE mo_input_opt
USE mo_input_base
USE mo_input_fileobj
USE mo_input_dimension_base
USE mo_input_mask
USE mo_input_dimension
USE mo_input_misc
USE mo_input_netcdf
USE mo_input_proc
IMPLICIT NONE

! Make all relevant entities from mo_input_types publically available
#ifdef HAVE_F2003
PUBLIC input_fcn_list, input_file_type_fcns_list
#endif
PUBLIC :: input_data_list, input_mask_list, input_curr, input_dim_list, input_dim_data, input_eqdim_list
PUBLIC :: input_fileobj_list, input_file_list, input_group_list, input_var_list, input_opt_list
PUBLIC :: INPUT_ACTION_OBTAIN, INPUT_ACTION_AFTER_READ, INPUT_ACTION_AFTER_CONDENSE, INPUT_ACTION_AFTER_BROADCAST
PUBLIC :: INPUT_ACTION_AFTER_RESCALE, INPUT_ACTION_AFTER_REVERSE, INPUT_ACTION_AFTER_PACK, INPUT_ACTION_AFTER_INTERPOLATE
PUBLIC :: INPUT_PROC_DATA_ACTION, INPUT_PROC_PROGRESS, INPUT_PROC_DO_DATA_ACTION, INPUT_PROC_NAME
PUBLIC :: INPUT_PROC_MOD_DIM_FORWARD, INPUT_PROC_MOD_DIM_BACKWARD
PUBLIC :: INPUT_DATA_NREF_MASK, INPUT_DATA_NO_STEP, INPUT_DATA_NO_MASK, INPUT_DATA_ALLOC, INPUT_DATA_GLOBAL, INPUT_DATA_IN_USE
PUBLIC :: INPUT_DATA_DEFAULT, INPUT_DATA_FILL, INPUT_DATA_RANGE, INPUT_DATA_RESCALE, INPUT_DATA_NUDGE_SCA, INPUT_DATA_NUDGE_FLD 
PUBLIC :: INPUT_DATA_WEIGHTS, INPUT_DATA_FUTURE_FRAC, INPUT_DATA_USER
PUBLIC :: INPUT_MASK_FIRST_VALID, INPUT_MASK_TYPE_MODEL, INPUT_MASK_TYPE_DISTRIB, INPUT_MASK_TYPE_PE, INPUT_MASK_TYPE_EQUIV
PUBLIC :: INPUT_MASK_GLOBAL, INPUT_MASK_CHUNK, INPUT_MASK_LOCAL, INPUT_MASK_CONT, INPUT_MASK_PROCESSED
PUBLIC :: INPUT_FILEOBJ_UNDET, INPUT_FILEOBJ_DIM, INPUT_FILEOBJ_VAR
PUBLIC :: INPUT_FILEDIM_DIR_UNDET, INPUT_FILEDIM_DIR_ERROR, INPUT_FILEDIM_DIR_UNLIM, INPUT_FILEDIM_DIR_DEC, INPUT_FILEDIM_DIR_INC
PUBLIC :: INPUT_FILEVAR_HAS_SCALE, INPUT_FILEVAR_HAS_BAD, INPUT_FILEVAR_HAS_RANGE, INPUT_FILEVAR_MASK
PUBLIC :: INPUT_DIM_UNLIM, INPUT_DIM_PACKED, INPUT_DIM_EXIST, INPUT_DIM_REVERSE, INPUT_DIM_CYCLIC, INPUT_DIM_COLLAPSED
PUBLIC :: INPUT_DIM_PARALLEL, INPUT_DIM_DISTRIB, INPUT_DIM_INC, INPUT_DIM_DEC
PUBLIC :: INPUT_DIM_NON_MODEL, INPUT_DIM_READ, INPUT_DIM_LOCAL_DATA 
PUBLIC :: INPUT_DIM_EXTRA, INPUT_DIM_FILE, INPUT_DIM_VAR, INPUT_DIM_GLOBAL, INPUT_DIM_SOURCE, INPUT_DIM_CHUNKED
PUBLIC :: INPUT_DIM_EQ_SRC, INPUT_DIM_EQ_DST, INPUT_DIM_USEASGLOBAL, INPUT_DIM_TYPE
PUBLIC :: INPUT_DIM_FND_NOT, INPUT_DIM_FND_NORMAL, INPUT_DIM_FND_SOURCE, INPUT_DIM_FND_PACKED, INPUT_DIM_FND_MASK
PUBLIC :: INPUT_DIM_FND_REVERSE, INPUT_DIM_FND_SWAP, INPUT_DIM_FND_SCATTER, INPUT_DIM_FND_DISTRIB 
PUBLIC :: INPUT_DIM_FND_EQ_SRC, INPUT_DIM_FND_EQ_DST
PUBLIC :: INPUT_FILE_INITIAL, INPUT_FILE_TIME_ABSOLUTE, INPUT_FILE_TIME_FROM_NAME, INPUT_FILE_NOCYCLE
PUBLIC :: INPUT_FILE_NAME_YEAR_SIGN, INPUT_FILE_NAME_YEAR_LONG, INPUT_FILE_NAME_YEAR, INPUT_FILE_NAME_MONTH, INPUT_FILE_NAME_DAY
PUBLIC :: INPUT_FILE_NAME_HOUR, INPUT_FILE_NAME_MINUTE, INPUT_FILE_NAME_SECOND, INPUT_FILE_MULTI
PUBLIC :: INPUT_FILE_TIME_ERR, INPUT_FILE_TIME_WARN, INPUT_FILE_TIME_IGNORE
PUBLIC :: INPUT_FILE_OPEN, INPUT_FILE_MULTI_OPENED, INPUT_FILE_DIMS_GOT, INPUT_FILE_VARS_GOT, INPUT_FILE_TIME_GOT
PUBLIC :: INPUT_FILE_FIRST_CYCLE, INPUT_FILE_END_REACHED, INPUT_FILE_CYCLED
PUBLIC :: INPUT_FILE_READ, INPUT_FILE_ENDDATE_GOT, INPUT_FILE_EXTERN_APPLIED
PUBLIC :: INPUT_FILE_CHECK_DIM_SIZE, INPUT_FILE_CHECK_NRECS, INPUT_FILE_CHECK_TIME_STEP, INPUT_FILE_CHECK_TIME_CYCLE
PUBLIC :: INPUT_FILE_CHECK_VAR, INPUT_FILE_CHECK_VAR_DIMS, INPUT_FILE_CHECK_INITIAL, INPUT_FILE_CHECK_MASK
PUBLIC :: INPUT_FILE_REPEAT_NONE, INPUT_FILE_REPEAT_ALL, INPUT_FILE_REPEAT_YEAR, INPUT_FILE_REPEAT_DAY, INPUT_FILE_REPEAT_RECS
PUBLIC :: INPUT_FILE_REPEAT_SECS, INPUT_FILE_REPEAT_MONTHS
PUBLIC :: INPUT_GROUP_MULTI_FILE, INPUT_GROUP_MSG_DONE
PUBLIC :: INPUT_VAR_MISS_DEFAULT, INPUT_VAR_MISS_IGNORE, INPUT_VAR_MISS_WARN, INPUT_VAR_MISS_WARN_IGNORE, INPUT_VAR_MISS_ERR
PUBLIC :: INPUT_VAR_MISS_ERR_IGNORE, INPUT_VAR_MISS_ERR_WARN
PUBLIC :: INPUT_VAR_BAD_VAL_DEFAULT, INPUT_VAR_BAD_VAL_REPLACE, INPUT_VAR_BAD_VAL_WARN, INPUT_VAR_BAD_VAL_WARN_REPLACE
PUBLIC :: INPUT_VAR_BAD_VAL_ERR, INPUT_VAR_BAD_VAL_IGNORE
PUBLIC :: INPUT_VAR_INTERPOL_USER, INPUT_VAR_GLOBAL, INPUT_VAR_AUTORESET_OFF
PUBLIC :: INPUT_VAR_VALID_ACTUAL, INPUT_VAR_VALID_BEFORE, INPUT_VAR_VALID_START, INPUT_VAR_VALID_MIDPOINT, INPUT_VAR_VALID_END 
PUBLIC :: INPUT_VAR_VALID_MIDDAY, INPUT_VAR_VALID_MIDMONTH, INPUT_VAR_VALID_MIDYEAR
PUBLIC :: INPUT_VAR_DIV_MONTH, INPUT_VAR_DIV_YEAR
PUBLIC :: INPUT_VAR_ACTION_REPLACE, INPUT_VAR_ACTION_ADD, INPUT_VAR_ACTION_SUB, INPUT_VAR_ACTION_SUBR, INPUT_VAR_ACTION_MUL
PUBLIC :: INPUT_VAR_ACTION_DIV, INPUT_VAR_ACTION_DIVR, INPUT_VAR_ACTION_ADD_FLUX, INPUT_VAR_ACTION_NUDGE,INPUT_VAR_ACTION_NUDGE_FLD
PUBLIC :: INPUT_VAR_DIM_MIS_ERR, INPUT_VAR_DIM_MIS_WARN, INPUT_VAR_DIM_MIS_VAR_MISS, INPUT_VAR_DIM_MIS_RESOLVE
PUBLIC :: INPUT_VAR_DIM_MIS_SPREAD, INPUT_VAR_DIM_MIS_SUBSET, INPUT_VAR_DIM_MIS_SUM, INPUT_VAR_DIM_MIS_AVG, INPUT_VAR_DIM_MIS_NORM
PUBLIC :: INPUT_VAR_INITIALIZED, INPUT_VAR_EXIST, INPUT_VAR_PROC_IN_FILE, INPUT_VAR_UPDATED, INPUT_VAR_READ, INPUT_VAR_HAS_UNLIM
PUBLIC :: INPUT_VAR_MULTI_SPREAD, INPUT_VAR_MULTI, INPUT_VAR_CONT_SEND, INPUT_VAR_BEFORE_FIRST, INPUT_VAR_DT_SET 
PUBLIC :: INPUT_VAR_DEPENDENT, INPUT_VAR_FILE_TIME_STEP
PUBLIC :: INPUT_VAR_DO_NOACTION, INPUT_VAR_DO_MESSAGE, INPUT_VAR_DO_TOREAD, INPUT_VAR_DO_MULTI_READ
PUBLIC :: INPUT_VAR_DO_SUM, INPUT_VAR_DO_AVG, INPUT_VAR_DO_NORM, INPUT_VAR_DO_SPREAD, INPUT_VAR_DO_SUBSET
PUBLIC :: INPUT_VAR_DO_REVERSE, INPUT_VAR_DO_SWAP, INPUT_VAR_DO_PERMUTE, INPUT_VAR_DO_PACK, INPUT_VAR_DO_UNPACK, INPUT_VAR_DO_EQUIV
PUBLIC :: INPUT_VAR_DO_FILL_BAD, INPUT_VAR_DO_FILL_RANGE, INPUT_VAR_DO_RESCALE, INPUT_VAR_DO_DIVTIME, INPUT_VAR_DO_INTERPOLATE
PUBLIC :: INPUT_VAR_DO_BROADCAST, INPUT_VAR_DO_SCATTER, INPUT_VAR_DO_COPY
PUBLIC :: INPUT_VAR_DO_COMBINE_ADD, INPUT_VAR_DO_COMBINE_SUB, INPUT_VAR_DO_COMBINE_SUBR, INPUT_VAR_DO_COMBINE_MUL
PUBLIC :: INPUT_VAR_DO_COMBINE_DIV, INPUT_VAR_DO_COMBINE_DIVR, INPUT_VAR_DO_COMBINE_ADD_FLUX 
PUBLIC :: INPUT_VAR_DO_COMBINE_NUDGE, INPUT_VAR_DO_COMBINE_NUDGE_FLD, INPUT_VAR_DO_USER_ACTION
PUBLIC :: INPUT_DBG_TIME_YMDHMS, INPUT_DBG_TIME_YMDHM, INPUT_DBG_TIME_YMDH, INPUT_DBG_TIME_YMD, INPUT_DBG_TIME_YM, INPUT_DBG_TIME_Y
PUBLIC :: INPUT_DBG_AVG, INPUT_DBG_EXTERN, INPUT_DBG_DBG_FLUSH, INPUT_DBG_DBG_IO, INPUT_DBG_DBG_IODATA, INPUT_DBG_DBG_ALL
PUBLIC :: INPUT_DBG_QUIET, INPUT_DBG_WARN, INPUT_DBG_OPENCLOSE, INPUT_DBG_READ_SUMMARY, INPUT_DBG_READ_GROUP, INPUT_DBG_READ_FILE
PUBLIC :: INPUT_DBG_READ_VAR, INPUT_DBG_UPDATE_SUMMARY, INPUT_DBG_UPDATE_GROUP, INPUT_DBG_UPDATE_VAR, INPUT_DBG_DEFAULT
PUBLIC :: my_pe

! Public procedures from mo_input_data
PUBLIC :: InputDataNew, InputDataDone, InputDataUnuse, InputDataDisp, InputDataGetRefGlobal, InputDataGetRef
PUBLIC :: InputDataGetSpecified, InputDataIntrpAlloc

! Public procedures from mo_input_debug
PUBLIC :: InputIsDebug, InputDebugLevelSet

! Public procedures from mo_input_base
PUBLIC :: GetRecTime, GetTimeScaleFromStr
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
PUBLIC :: InputParSetSourcePE, InputVarBroadcast, InputVarScatter
#endif

! Public procedures from mo_input_fileobj
PUBLIC :: InputFileObjNew, InputFileObjsDone, InputFileObjGetRef, InputFileObjExpandVar, InputFileObjsCmp, InputFileObjDisp

! Public procedures from mo_input_opt
PUBLIC :: InputOptNew, InputOptDisp, InputOptImport

! Public procedures from mo_input_dimension_base
PUBLIC InputDimNewCtl, InputDimNew, InputDimDone, InputDimListElmRemove, InputDimCopyList
PUBLIC InputDimGetRef, InputDimGetRefEx, InputDimDisp, InputDimInq

! Public procedures from mo_input_mask
PUBLIC InputMaskNew, InputMaskCopy, InputMaskDisp, InputMasksDone, InputMaskCreate, InputMaskGetTotal, InputMaskGetValid
PUBLIC InputMaskCreateScatter, InputMaskReduce_multi, InputMaskReduce_single

! Public procedures from mo_input_dimension
PUBLIC :: InputDimCopy, InputDimListFromString, InputDimAdd, InputDimPackedAdd, InputDimPackedAddMask, InputDimLocalSet
PUBLIC :: InputDimLocalSetMulti, InputDimSync, InputDimEquivalence, InputEqDimDisp

! Public procedures from mo_input_misc
PUBLIC :: InputSubModelGetNew, InputSubModelGetCurr, InputSubModelSelect
PUBLIC :: InputFileSplitName, BuildGenericName, GetFileName
PUBLIC :: InputTimeGet, InputTimeSet, InputTimeString
PUBLIC :: InputVarInterpolBoxCar, InputVarInterpolLinear, InputVarGetInterpolData

#ifdef HAVE_F2003
PUBLIC :: InputFcnNew, InputFcnAdd, InputFcnDone
#endif
! Public procedures from mo_input_proc
PUBLIC :: GroupReadMsg
PUBLIC :: DataProcNoAction, DataProcMessage, DataProcRead, DataProcReadMulti, DataProcPermute, DataProcReverse, DataProcSwap
PUBLIC :: DataProcPack, DataProcUnpack, DataProcEquiv, DataProcCopy, DataProcSum, DataProcAvg, DataProcNormalize
PUBLIC :: DataProcFillBad, DataProcFillOutOfRange, DataProcRescale, DataProcTimeDiv, DataProcSpread, DataProcSubset
PUBLIC :: DataProcInterpolateTimeLinear, DataProcBroadcast, DataProcScatter
PUBLIC :: DataProcAdd, DataProcSub, DataProcSubr, DataProcMul, DataProcDiv, DataProcDivr
PUBLIC :: DataProcFlux, DataProcNudge, DataProcNudgeFld
#ifdef HAVE_F2003
PUBLIC :: InputVarGetNextBuffer, InputVarGetNextMask
#endif

PRIVATE
INTERFACE InputVarAdd
  MODULE PROCEDURE InputVarNDAdd
  MODULE PROCEDURE InputVar0dAdd
  MODULE PROCEDURE InputVar1dAdd
  MODULE PROCEDURE InputVar2dAdd
  MODULE PROCEDURE InputVar3dAdd
  MODULE PROCEDURE InputVar4dAdd
  MODULE PROCEDURE InputVar5dAdd
  MODULE PROCEDURE InputVar6dAdd
  MODULE PROCEDURE InputVar7dAdd
END INTERFACE
INTERFACE InputVarDataPush
  MODULE PROCEDURE InputVarDataPush0d
  MODULE PROCEDURE InputVarDataPush1d
  MODULE PROCEDURE InputVarDataPush2d
  MODULE PROCEDURE InputVarDataPush3d
  MODULE PROCEDURE InputVarDataPush4d
  MODULE PROCEDURE InputVarDataPush5d
  MODULE PROCEDURE InputVarDataPush6d
  MODULE PROCEDURE InputVarDataPush7d
END INTERFACE
INTERFACE InputVarUserDataSet
  MODULE PROCEDURE InputVarUserDataSet0d
  MODULE PROCEDURE InputVarUserDataSet1d
  MODULE PROCEDURE InputVarUserDataSet7d
END INTERFACE

! ---------------------------------------------------------------------
! Public subroutines and functions
PUBLIC InputInit
PUBLIC InputClose
PUBLIC InputGetData
PUBLIC InputUpdate
PUBLIC InputFileAdd
PUBLIC InputFileGetRef
PUBLIC InputFileRemove
PUBLIC InputFileDisp
PUBLIC InputFileInq
PUBLIC InputFileOpen
PUBLIC InputFileClose
PUBLIC InputVarNew
PUBLIC InputVarAdd
PUBLIC InputVarNewDummy
PUBLIC InputVarGetRef
PUBLIC InputVarDisp
PUBLIC InputVarInq
PUBLIC InputVarDataPush
PUBLIC InputVarUpdated
PUBLIC InputVarLastUpdate
PUBLIC InputVarLastRead
PUBLIC InputVarIsInitialized
PUBLIC InputVarIsPresent
PUBLIC InputVarIsUpdated
PUBLIC InputVarIsRead
PUBLIC InputVarWillBeUpdated
PUBLIC InputVarWillBeRead
PUBLIC InputVarIsExternal
PUBLIC InputVarUserDataGet
PUBLIC InputVarUserDataSet
PUBLIC InputVarGroupNew
PUBLIC InputVarGroupGetRef
PUBLIC InputVarGroupDisp
PUBLIC InputVarGroupStart
PUBLIC InputVarGroupEnd
PUBLIC InputAttrGet
PUBLIC InputVarGetDimForOutput
#if defined(INPUT_IN_ECHAM) || defined(INPUT_IN_ICON)
PUBLIC InputVar2Output
#endif
#ifdef HAVE_F2003
PUBLIC InputAddDataAction
PUBLIC InputVarNewBuffer
#endif

CONTAINS

! ----------------------------------------------------
!  Creating, adding, setting, modifying, displaying:
!     Files
! ----------------------------------------------------

! Create a default input_file_list element
  SUBROUTINE InputFileNew(nfile)
  TYPE (input_file_list), POINTER :: nfile

    ALLOCATE(nfile)
    NULLIFY(nfile%next,nfile%dims,nfile%fvars,nfile%fdims)
#ifdef HAVE_F2003
    NULLIFY(nfile%ft)
#endif
    nfile%name_pre        = ''
    nfile%name_suf        = ''
    nfile%name_file       = ''
    nfile%time_unit       = ''
    nfile%sub_model       =  1
    nfile%flags           =  0
    nfile%stat            =  0
    nfile%nUsed           =  0
    nfile%nVar            =  0
    nfile%nRead           =  0
    nfile%nRead_dbg       =  0
    nfile%fid             =  0
    nfile%nRec            = -1
    nfile%curr_rec        =  0
    nfile%init_rec        =  UNLIKELY_VAL
    nfile%rec_step        =  0
    nfile%dt_file         =  UNLIKELY_VAL
    nfile%time_scale(:)   =  0
    nfile%time_offset(:)  =  0
    nfile%time_fid(:)     =  0
    nfile%time_open(:)    = -HUGE(i8)
    nfile%time_file(:)    = -HUGE(i8)
    nfile%time_start(:)   = -HUGE(i8)
    nfile%time_enddata(:) = -HUGE(i8)
    nfile%time_cycle(1:4) = -HUGE(i8)
    nfile%time_cycle(5:7) = -1

  END SUBROUTINE InputFileNew

! Returns handle to a file which matches the name and sub-model
  FUNCTION InputFileGetRef(name_pre,name_suf,flags,sub_model)
  TYPE (input_file_list), POINTER :: InputFileGetRef
  CHARACTER(len=*),    INTENT(in) :: name_pre, name_suf
  INTEGER,             INTENT(in) :: flags
  INTEGER,   OPTIONAL, INTENT(in) :: sub_model

  CHARACTER (len=64) :: pre, suf 
  TYPE (input_file_list), POINTER :: curr
  INTEGER :: sm, tmpflags

    sm = CurrSubModel
    IF (PRESENT(sub_model)) sm = sub_model

    pre = name_pre
    suf = name_suf
    tmpflags = flags
    IF ((suf=='') .AND. (flags==0)) CALL InputFileSplitName(name_pre,pre,suf,tmpflags)

    curr => AllInputFiles
    DO WHILE (ASSOCIATED(curr))
      IF ((curr%sub_model < 0) .OR. (sm < 0) .OR. (curr%sub_model == sm)) THEN
        IF ((TRIM(pre)==TRIM(curr%name_pre)) .AND. &
            (TRIM(suf)==TRIM(curr%name_suf)) .AND. &
            (IAND(tmpflags,INPUT_FILE_MULTI) == IAND(curr%flags,INPUT_FILE_MULTI))) THEN
          InputFileGetRef => curr
          RETURN
        ENDIF
      ENDIF
      curr => curr%next
    ENDDO

    NULLIFY(InputFileGetRef)

  END FUNCTION InputFileGetRef

! Remove an input file from the list of files to be processed
  SUBROUTINE InputFileRemove(IFile)
  TYPE (input_file_list), POINTER :: IFile

  TYPE (input_file_list), POINTER :: curr_file, prev_file
  TYPE (input_var_list),  POINTER :: curr_var,  prev_var

    NULLIFY(prev_file)
    curr_file => AllInputFiles
    DO WHILE (ASSOCIATED(curr_file))
      IF (ASSOCIATED(curr_file,IFile)) THEN
        IF (IAND(IFile%stat,INPUT_FILE_OPEN)/=0) CALL InputFileClose(IFile)
        IF (ASSOCIATED(prev_file)) THEN
          prev_file%next => curr_file%next
        ELSE
          AllInputFiles => curr_file%next
        ENDIF
        ! Remove all variables associated with this file
        NULLIFY(prev_var)
        curr_var => AllInputVars
        DO WHILE (ASSOCIATED(curr_var))
          IF (ASSOCIATED(curr_var%parent_file,curr_file)) THEN
            IF (ASSOCIATED(prev_var)) THEN
              prev_var%next => curr_var%next
            ELSE
              AllInputVars => curr_var%next
            ENDIF
            NULLIFY(curr_var) ! TODO - here some memory freeing may still be missing
          ELSE
            prev_var => curr_var
            curr_var => curr_var%next
          ENDIF
        ENDDO
        NULLIFY(curr_file)
      ELSE
        prev_file => curr_file
        curr_file => curr_file%next
      ENDIF
    ENDDO

  END SUBROUTINE InputFileRemove

! Finds an instance of a multi-file if possible (Which takes the values: 0: Any (default), 1: First, 2: Last)
  FUNCTION InputFileMultiInstanceGet(IFile,Which,BaseTime,MaxIt)
  INTEGER :: InputFileMultiInstanceGet(4)
  TYPE (input_file_list), POINTER :: IFile
  INTEGER, OPTIONAL, INTENT(in)   :: Which
  INTEGER, OPTIONAL, INTENT(in)   :: BaseTime(4)
  INTEGER, OPTIONAL, INTENT(in)   :: MaxIt

  TYPE (input_var_list),  POINTER :: curr, grpvar
  CHARACTER (len=128) :: file_name
  INTEGER :: time_plus(4), time_minus(4), nIt, It
  REAL(dp), POINTER :: tst(:,:,:,:,:,:,:)
  LOGICAL :: lExist, lPlus

    IF (lio) THEN
      time_plus    = exp_start
      time_plus(1) = time_plus(1) + IFile%time_offset(1)
      IF (PRESENT(BaseTime)) THEN
        IF (ALL(BaseTime > -HUGE(It))) time_plus = BaseTime
      ENDIF
      time_minus   = InputFileTimePrevGet(IFile,time_plus)

      IF (PRESENT(MaxIt)) THEN
        nIt = MaxIt
      ELSE
        IF (IAND(IFile%flags,INPUT_FILE_NAME_SECOND+INPUT_FILE_NAME_MINUTE+INPUT_FILE_NAME_HOUR) /= 0) nIt = 100000000
        IF (IAND(IFile%flags,INPUT_FILE_NAME_DAY)   /= 0) nIt = 43800000
        IF (IAND(IFile%flags,INPUT_FILE_NAME_MONTH) /= 0) nIt = 120000
        IF (IAND(IFile%flags,INPUT_FILE_NAME_YEAR)  /= 0) nIt = 10000
      ENDIF
      IF (IAND(IFile%flags,INPUT_FILE_MULTI) == 0) nIt = 2
      It = 0
      lExist = .FALSE.
      DO WHILE (.NOT. lExist .AND. (It < nIt))
        lPlus = .TRUE.
        file_name = GetFileName(IFile,time_plus(1),time_plus(2),time_plus(3),time_plus(4))
        INQUIRE(file=TRIM(file_name),exist=lExist)
        IF (.NOT. lExist) THEN
          file_name = GetFileName(IFile,time_minus(1),time_minus(2),time_minus(3),time_minus(4))
          INQUIRE(file=TRIM(file_name),exist=lExist)
          lPlus = .FALSE.
        ENDIF
        IF (.NOT. lExist) THEN
          time_plus  = InputFileTimeNextGet(IFile,time_plus)
          time_minus = InputFileTimePrevGet(IFile,time_minus)
        ENDIF
        It = It + 1
      ENDDO
      IF (lPlus .AND. It < nIt) time_minus = time_plus

      IF (PRESENT(Which) .AND. (IAND(IFile%flags,INPUT_FILE_MULTI) /= 0) .AND. It < nIt) THEN
        SELECT CASE (Which)
          CASE (0) ! Just some instance - nothing to do 
          CASE (1) ! Find first instance
            IF ((.NOT. lPlus .AND. (It>1)) .OR. PRESENT(BaseTime)) THEN
              DO WHILE (lExist)
                time_plus = InputFileTimePrevGet(IFile,time_minus)
                file_name = GetFileName(IFile,time_plus(1),time_plus(2),time_plus(3),time_plus(4))
                INQUIRE(file=TRIM(file_name),exist=lExist)
                IF (lExist) time_minus = time_plus
              ENDDO
            ENDIF
          CASE (2) ! Find last instance
            IF ((lPlus .AND. (It>1)) .OR. PRESENT(BaseTime)) THEN
              DO WHILE (lExist)
                time_plus = InputFileTimeNextGet(IFile,time_minus)
                file_name = GetFileName(IFile,time_plus(1),time_plus(2),time_plus(3),time_plus(4))
                INQUIRE(file=TRIM(file_name),exist=lExist)
                IF (lExist) time_minus = time_plus
              ENDDO
            ENDIF
        END SELECT
      ENDIF
    ENDIF

#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
    CALL p_bcast(nIt,io_pe)
    CALL p_bcast(It,io_pe)
#endif

    ! No instance found?
    IF (It==nIt) THEN
      message_text = ''
      curr => AllInputVars
      DO WHILE (ASSOCIATED(curr))
        IF (ASSOCIATED(curr%parent_file,IFile)) THEN
          tst => InputDataGetSpecified(curr%aux_data,INPUT_DATA_DEFAULT,0)
          IF (lio .AND. (.NOT. ASSOCIATED(tst) .OR. (IAND(curr%flags,INPUT_VAR_MISS_MASK) >= INPUT_VAR_MISS_ERR))) &
            CALL local_error('InputFileOpen','No instance of multi-file: '//TRIM(BuildGenericName(IFile))//        &
                                             ' found and defaults are not available for all file variables')
          ! Set variable to default value, detach variable from non-existing file and issue warning if requested
          file_name = BuildGenericName(IFile)
          NULLIFY(curr%parent_file)
          IF (ASSOCIATED(curr%dta)) THEN
            IF (SIZE(tst) == 1) THEN
              curr%dta(:,:,:,:,:,:,:) = tst(1,1,1,1,1,1,1)
            ELSE
              IF (.NOT. ArrayCopy(tst,curr%dta,(/SIZE(tst)/),SHAPE(curr%dta))) &
                CALL local_error('InputFileOpen','Spread error')
            ENDIF
          ENDIF
          curr%time = HUGE(It)
          IF (lio .AND. (IAND(curr%flags,INPUT_VAR_MISS_MASK) > INPUT_VAR_MISS_IGNORE)) &
            message_text = TRIM(message_text)//', '//TRIM(curr%name_var)
          ! Reassociate group-file of variable group if necessary and possible
          IF (ASSOCIATED(curr%group)) THEN
            IF (ASSOCIATED(curr%group%group_file,IFile)) THEN
              NULLIFY(curr%group%group_file)
              grpvar => AllInputVars
              DO WHILE (ASSOCIATED(grpvar))
                IF (ASSOCIATED(curr%group,grpvar%group) .AND. ASSOCIATED(grpvar%parent_file) .AND. &
                    .NOT. ASSOCIATED(grpvar%parent_file,IFile)) curr%group%group_file => grpvar%parent_file
                grpvar => grpvar%next
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        curr => curr%next
      ENDDO
      IF (LEN_TRIM(message_text) > 0) CALL local_message('WARNING','No instance of multi-file: '//                 &
        TRIM(file_name)//' found. Using default value(s) for variable(s): '//TRIM(message_text(3:)))
      time_minus(:) = -HUGE(It)
    ENDIF

#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
    CALL p_bcast(time_minus,io_pe)
#endif

    InputFileMultiInstanceGet = time_minus

  END FUNCTION InputFileMultiInstanceGet

! Display files in a list
  SUBROUTINE InputFileDisp(fil,cnt,unt,msg,list)
  TYPE (input_file_list), OPTIONAL, POINTER :: fil
  INTEGER, OPTIONAL, INTENT(in) :: cnt
  INTEGER, OPTIONAL, INTENT(in) :: unt
  CHARACTER (len=*), INTENT(in), OPTIONAL :: msg
  LOGICAL, OPTIONAL, INTENT(in) :: list

  TYPE (input_file_list), POINTER :: curr_file
  TYPE (input_var_list),  POINTER :: curr_var
  TYPE (input_dim_list),  POINTER :: curr_dim
  CHARACTER (len=256) :: tmpstr
  INTEGER :: un, cr

    un = GetLogUnit(); IF (PRESENT(unt)) un = unt
    cr = -1;           IF (PRESENT(cnt)) cr = cnt

    IF (PRESENT(fil)) THEN
      curr_file => fil
    ELSE
      curr_file => AllInputFiles
    ENDIF
    IF (PRESENT(list)) THEN
      IF (list) THEN
        tmpstr = ''
        DO WHILE (ASSOCIATED(curr_file) .AND. (cr /= 0))
          tmpstr = TRIM(tmpstr)//','//TRIM(BuildGenericName(curr_file))
          curr_file => curr_file%next
          cr = cr - 1
        ENDDO
        IF (PRESENT(msg)) THEN
          WRITE(un,'(2a)') msg,TRIM(tmpstr(2:))
        ELSE
          WRITE(un,'(a)') TRIM(tmpstr(2:))
        ENDIF
        RETURN
      ENDIF
    ENDIF
    IF (PRESENT(msg)) WRITE(un,'(a)') msg
    DO WHILE (ASSOCIATED(curr_file) .AND. (cr /= 0))

      ! File name
      tmpstr = BuildGenericName(curr_file)
      WRITE(tmpstr,'(2a,i2)') TRIM(tmpstr),'%(',curr_file%sub_model
      IF ( curr_file%name_file /= '' .AND. (IAND(curr_file%flags,INPUT_FILE_MULTI) /= 0)) &
        tmpstr = TRIM(tmpstr)//', current: '//TRIM(curr_file%name_file)
      CALL RmSpc(tmpstr,'%')
      WRITE(un,'(2a)') TRIM(tmpstr),')'

      ! File type
#ifdef HAVE_F2003
      WRITE(un,'(2a)') '  File type                    : ',TRIM(curr_file%ft%FileType)
#endif

      ! #Declared
      IF (curr_file%nUsed > 1) THEN
        WRITE(tmpstr,*) curr_file%nUsed
        CALL RmDblSpc(tmpstr)
        WRITE(unt,'(2a)') '  #Declarations:',TRIM(tmpstr)
      ENDIF

      ! Dimensions
      tmpstr = ''
      curr_dim => curr_file%dims
      DO WHILE (ASSOCIATED(curr_dim))
        tmpstr = TRIM(tmpstr)//', '//TRIM(curr_dim%dim_data%name_dim)
        IF (curr_dim%file_dim_name /= '') THEN
          tmpstr = TRIM(tmpstr)//' (alias: '//TRIM(curr_dim%file_dim_name)//')'
        ELSE
          IF (curr_dim%dim_data%name_alt /= '') tmpstr = TRIM(tmpstr)//' ('//TRIM(curr_dim%dim_data%name_alt)//')'
        ENDIF
        curr_dim => curr_dim%next
      ENDDO
      WRITE(un,'(2a)') '  Dimensions: ',TRIM(tmpstr(3:))

      ! Variables
      IF (curr_file%nVar > 0) THEN
        WRITE(tmpstr,'(a,i2,a)') '__Variables_(',curr_file%nVar,'):'
        CALL RmSpc(tmpstr)
        WRITE(un,'(a)') TRIM(tmpstr)
        curr_var => AllInputVars
        DO WHILE (ASSOCIATED(curr_var))
          IF (ASSOCIATED(curr_var%parent_file,curr_file)) THEN
            tmpstr = '    '//TRIM(curr_var%name_var)
            IF (LEN_TRIM(curr_var%name_alt) > 0) tmpstr = TRIM(tmpstr)//' ('//TRIM(curr_var%name_alt)//')'
            IF ((IAND(curr_file%stat,INPUT_FILE_CHECK_VAR) /= 0) .AND. (IAND(curr_var%stat,INPUT_VAR_EXIST) == 0)) &
              tmpstr = TRIM(tmpstr)//' - not present.'
            WRITE(un,'(a)') TRIM(tmpstr)
          ENDIF
          curr_var => curr_var%next
        ENDDO
      ENDIF

      ! Read status
      IF (IAND(curr_file%stat,INPUT_FILE_OPEN) /= 0) THEN
        WRITE(tmpstr,*) curr_file%nRead
        CALL RmDblSpc(tmpstr)
        WRITE(un,'(2a)') '  Currently read variables ... :',TRIM(tmpstr)
        WRITE(tmpstr,*) curr_file%curr_rec,'   of',curr_file%nRec
        CALL RmDblSpc(tmpstr)
        WRITE(un,'(2a)') '  ... from record              :',TRIM(tmpstr)
      ENDIF
      IF ((curr_file%init_rec /= UNLIKELY_VAL) .AND. (curr_file%init_rec > 0)) THEN
        WRITE(tmpstr,*) curr_file%init_rec
      ELSE
        WRITE(tmpstr,*) curr_file%nRec
      ENDIF
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)') '  Initial var. from record     :',TRIM(tmpstr)

      ! Times
      WRITE(un,'(2a)')     '  File time step               : ',TimeStepString(curr_file%dt_file,UNLIKELY_VAL)
      WRITE(message_text,*) curr_file%rec_step
      CALL RmDblSpc(message_text)
      WRITE(un,'(2a)')     '  File record step             :',TRIM(message_text)
      IF (IAND(curr_file%flags,INPUT_FILE_MULTI) /= 0) THEN
        WRITE(un,'(a,a19)')  '  Current file name time       : ',DateString(curr_file%time_open)
        WRITE(un,'(a,a19)')  '  Next    file name time       : ',DateString(curr_file%time_file)
      ENDIF
      WRITE(un,'(a,a19)')  '  First available data time    : ',DateString(curr_file%time_start(1:4))
      WRITE(un,'(a,a19)')  '  First unavailable data time  : ',DateString(curr_file%time_enddata(1:4))
      WRITE(un,'(a,a19)')  '  Data cycle start time        : ',DateString(curr_file%time_cycle(1:4))
      WRITE(tmpstr,*) MAX(-1,curr_file%time_start(5)), MAX(-1,curr_file%time_enddata(5)), curr_file%time_cycle(5)
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)')   '  Record # of start, end, cycle:',TRIM(tmpstr)
      WRITE(tmpstr,*) MAX(-1,curr_file%time_enddata(6)),MAX(-1,curr_file%time_cycle(6))
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)')   '  Total #rec and #recs repeated:',TRIM(tmpstr)
      WRITE(tmpstr,*) curr_file%time_offset(4)
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)')   '  File record offset           :',TRIM(tmpstr)
      WRITE(tmpstr,*) curr_file%time_offset(1)
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)')   '  Year offset (model->data)    :',TRIM(tmpstr)
      WRITE(tmpstr,*) curr_file%time_offset(2)
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)')   '  Time check day correction    :',TRIM(tmpstr)
      WRITE(tmpstr,*) curr_file%time_offset(3)
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)')   '  Time shift (sec)             :',TRIM(tmpstr)

      ! NetCDF internals
      WRITE(tmpstr,*) curr_file%fid
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)')   '  File identification number   :',TRIM(tmpstr)
      WRITE(tmpstr,*) curr_file%time_fid(1)
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)')   '  Record separating dimension  :',TRIM(tmpstr)

      ! Flags
      WRITE(un,'(2(a,Z8.8))') '  File flags (flags,status)    : 0x',curr_file%flags,', 0x',curr_file%stat
      CALL DispFlagsSingleBit(curr_file%stat,INPUT_FILE_OPEN,flg_file_stat_n,flg_file_stat,Txt=message_text)
      WRITE(un,'(2a)')     '  File status                  : ',TRIM(message_text)
      CALL DispFlagsSingleBit(curr_file%flags,INPUT_FILE_INITIAL,flg_file_time_n,flg_file_time,Txt=message_text)
      WRITE(un,'(2a)')     '  File time treatment          : ',TRIM(message_text)
      WRITE(un,'(2a)')     '  File time mismatch action    : ',TRIM(GetStringIndexed(opt_time_err, &
                                                               IAND(curr_file%flags,INPUT_FILE_TIME_MASK)/INPUT_FILE_TIME_WARN+1))
      CALL DispFlagsSingleBit(curr_file%stat,INPUT_FILE_CHECK_DIM_SIZE,flg_file_check_n,flg_file_check,Txt=message_text)
      WRITE(un,'(2a)')     '  File checks done             : ',TRIM(message_text)

      ! Next file
      curr_file => curr_file%next
      cr = cr - 1
    ENDDO

  END SUBROUTINE InputFileDisp

! Returns various information about the file
  SUBROUTINE InputFileInq(name_file,ref,name_var,var_ref,fid,name_curr,nVar,nRec,dt_file, &
                          time_curr,time_next,time_start,time_end,time_cycle,time_rec_curr,time_rec_prev)
  CHARACTER (len=*),      OPTIONAL, INTENT(in)  :: name_file,name_var
  TYPE (input_file_list), OPTIONAL, POINTER     :: ref
  TYPE (input_var_list),  OPTIONAL, POINTER     :: var_ref
  INTEGER,                OPTIONAL, INTENT(out) :: fid
  CHARACTER (len=128),    OPTIONAL, INTENT(out) :: name_curr
  INTEGER,                OPTIONAL, INTENT(out) :: nRec, nVar
  INTEGER,                OPTIONAL, INTENT(out) :: dt_file
  INTEGER,                OPTIONAL, INTENT(out) :: time_curr(4),     time_next(4), time_cycle(4), time_start(4), time_end(4)
  INTEGER,                OPTIONAL, INTENT(out) :: time_rec_curr(4), time_rec_prev(4)

  TYPE (input_file_list), POINTER :: curr
  TYPE (input_var_list),  POINTER :: var
  INTEGER :: tmm(4)

    NULLIFY(curr,var)
    IF (PRESENT(ref)) THEN
      IF (ASSOCIATED(ref)) curr => ref
    ENDIF
    IF (PRESENT(name_file)) curr => InputFileGetRef(TRIM(name_file),'',0)
    IF (PRESENT(var_ref))  var => var_ref
    IF (PRESENT(name_var)) var => InputVarGetRef(name_var)
    IF (ASSOCIATED(var)) THEN
      IF (ASSOCIATED(var%parent_file)) THEN
        IF (ASSOCIATED(curr) .AND. .NOT. ASSOCIATED(var%parent_file,curr)) &
          CALL local_error('InputFileInq','File and variable reference refere to different files')
        curr => var%parent_file
      ENDIF
    ENDIF
    IF (PRESENT(ref)) THEN
      IF (ASSOCIATED(ref)) THEN
        IF (ASSOCIATED(curr) .AND. .NOT. ASSOCIATED(curr,ref)) &
          CALL local_error('InputFileInq','Name/variable reference and file reference refere to different files')
        curr => ref
      ELSE
        IF (ASSOCIATED(curr)) ref => curr
      ENDIF
    ENDIF
    IF (.NOT. ASSOCIATED(curr)) CALL local_error('InputFileInq','No method to identify requested file')

    IF (PRESENT(name_curr))  name_curr  = TRIM(curr%name_file)
    IF (PRESENT(fid))        fid        = curr%fid
    IF (PRESENT(nRec))       nRec       = curr%nRec
    IF (PRESENT(nVar))       nVar       = curr%nVar
    IF (PRESENT(dt_file))    dt_file    = curr%dt_file
    IF (PRESENT(time_curr))  time_curr  = curr%time_open   (1:4)
    IF (PRESENT(time_next))  time_next  = curr%time_file   (1:4)
    IF (PRESENT(time_start)) time_start = curr%time_start  (1:4)
    IF (PRESENT(time_cycle)) time_cycle = curr%time_cycle  (1:4)
    IF (PRESENT(time_end))   time_end   = curr%time_enddata(1:4)
    IF (PRESENT(time_rec_curr)) THEN
#ifdef HAVE_F2003
      IF (lio) tmm = IFile%ft%FcnGetRecTime(curr,1,.FALSE.) ! Time of first record in found file
#else
      IF (lio) tmm = NCFileGetRecTime(curr,1,.FALSE.)
#endif
      CALL CalendarTimeAdd(tmm,INT(curr%dt_file,i8)*INT((curr%curr_rec-1),i8))
      time_rec_curr = tmm
    ENDIF
    IF (PRESENT(time_rec_prev)) THEN
#ifdef HAVE_F2003
      IF (lio) tmm = IFile%ft%FcnGetRecTime(curr,1,.FALSE.) ! Time of first record in found file
#else
      IF (lio) tmm = NCFileGetRecTime(curr,1,.FALSE.)
#endif
      CALL CalendarTimeAdd(tmm,INT(curr%dt_file,i8)*INT(MAX(0,curr%curr_rec-2),i8))
      time_rec_prev = tmm
    ENDIF

  END SUBROUTINE InputFileInq
    
! Set parameters of existing file based on external settings
  SUBROUTINE InputFileApplyExternalSettings(IFile,opt,cnt)
  TYPE (input_file_list), POINTER :: IFile
  TYPE (input_opt_list),  POINTER :: opt
  INTEGER,             INTENT(in) :: cnt

  TYPE (input_opt_list),  POINTER :: curr_opt
  TYPE (input_dim_list),  POINTER :: curr_dim, prev_dim
  INTEGER :: i, j, Le, St, En, Split
  LOGICAL :: lSet

    curr_opt => opt
    DO j=1,cnt

      ! Get labstime setting
      IF (TRIM(curr_opt%labstime) /= '') THEN
        i = GetStringIndex(curr_opt%labstime,opt_logical,ErrMsg='labstime')
        IF (i==2) IFile%flags = IOR (IFile%flags,                 INPUT_FILE_TIME_ABSOLUTE)
        IF (i==1) IFile%flags = IAND(IFile%flags,INPUT_ALLFLAGS - INPUT_FILE_TIME_ABSOLUTE)
      ENDIF

      ! Get name time setting
      IF (TRIM(curr_opt%lnametime) /= '') THEN
        i = GetStringIndex(curr_opt%lnametime,opt_logical,ErrMsg='lnametime')
        IF (i==2) IFile%flags = IOR (IFile%flags,                 INPUT_FILE_TIME_FROM_NAME)
        IF (i==1) IFile%flags = IAND(IFile%flags,INPUT_ALLFLAGS - INPUT_FILE_TIME_FROM_NAME)
      ENDIF

      ! Get time error handling
      IF (TRIM(curr_opt%action_time_err) /= '') THEN
        i = GetStringIndex(curr_opt%action_time_err,opt_time_err,ErrMsg='action_time_err')
        IF (i > 0) IFile%flags = IOR(IAND(IFile%flags,INPUT_ALLFLAGS-INPUT_FILE_TIME_MASK),(i-1)*INPUT_FILE_TIME_WARN)
      ENDIF

      ! Set time parameters
      IF ((curr_opt%dt_file /= UNLIKELY_VAL) .OR. (TRIM(curr_opt%dt_file_unit) /= '')) THEN
        IFile%dt_file = InterpretDtUnit(curr_opt%dt_file,curr_opt%dt_file_unit)
        IF (IFile%dt_file == 0) IFile%flags = IOR(IFile%flags,INPUT_FILE_INITIAL)
      ENDIF

      ! Start and end of available data
      IF (TRIM(curr_opt%data_start) /= '') THEN 
        CALL GetDateFromString(curr_opt%data_start,IFile%time_start(1:4),'data_start')
        IF (IFile%time_start(2) == -HUGE(i)) IFile%time_start(2) = 1
        IF (IFile%time_start(3) == -HUGE(i)) IFile%time_start(3) = 1
        IFile%time_start(4) = 0
      ENDIF
      IF (TRIM(curr_opt%data_end) /= '') THEN 
        CALL GetDateFromString(curr_opt%data_end,IFile%time_enddata(1:4),'data_end')
        IF (IFile%time_enddata(2)==-HUGE(i)) IFile%time_enddata(2) = 12
        IF (IFile%time_enddata(3)==-HUGE(i)) IFile%time_enddata(3) = CalendarMonthLen(IFile%time_enddata(1),IFile%time_enddata(2))
        IFile%time_enddata(4) = 0
      ENDIF

      ! Offset times and initial variable record
      IF (curr_opt%offset_year /= UNLIKELY_VAL) IFile%time_offset(1) = curr_opt%offset_year
      IF (curr_opt%offset_day  /= UNLIKELY_VAL) IFile%time_offset(2) = curr_opt%offset_day
      IF (curr_opt%offset_sec  /= UNLIKELY_VAL) IFile%time_offset(3) = curr_opt%offset_sec
      IF (curr_opt%offset_rec  /= UNLIKELY_VAL) IFile%time_offset(4) = curr_opt%offset_rec
      IF (curr_opt%init_rec    /= UNLIKELY_VAL) IFile%init_rec       = curr_opt%init_rec

      ! Prepare cycle time
      IF (TRIM(curr_opt%action_cycle) /= '') THEN
        i = GetStringIndex(TRIM(curr_opt%action_cycle),opt_cycle,ErrMsg='action_cycle')
        IF (i>0) IFile%time_cycle(1) = HUGE(i8)-i
        IF (i>2) IFile%time_cycle(4) = curr_opt%cycle_length
        IF (i==1) THEN
          IFile%flags = IOR (IFile%flags,                 INPUT_FILE_NOCYCLE)
        ELSE
          IFile%flags = IAND(IFile%flags,INPUT_ALLFLAGS - INPUT_FILE_NOCYCLE)
        ENDIF
      ENDIF

      ! Set renamed dimension names
      NULLIFY(prev_dim)
      Le = LEN_TRIM(curr_opt%dim_rename)
      St = 1
      DO WHILE (St < Le)
        En = INDEX(curr_opt%dim_rename(St:),',') + St - 1
        IF (En < St) En = Le + 1
        Split = INDEX(curr_opt%dim_rename(St:En-1),'=>') + St - 1
        IF (Split < St) CALL local_error('InputFileSetParams','Dimension rename format error: '//curr_opt%dim_rename(St:En-1))
        Split = Split + 2
        DO WHILE (curr_opt%dim_rename(Split:Split)==' ')
          Split = Split + 1
        ENDDO
        lSet = .FALSE.
        curr_dim => IFile%dims
        DO WHILE (ASSOCIATED(curr_dim))
          IF (TRIM(curr_dim%dim_data%name_dim) == TRIM(curr_opt%dim_rename(Split:En-1))) THEN
            curr_dim%file_dim_name = TRIM(curr_opt%dim_rename(St:Split-3))
            lSet = .TRUE.
            NULLIFY(curr_dim)
          ELSE
            curr_dim => curr_dim%next
          ENDIF
        ENDDO
        IF (.NOT. lSet) THEN
          CALL InputDimNew(curr_dim)
          IF (ASSOCIATED(prev_dim)) THEN
            prev_dim%next => curr_dim
          ELSE
            IFile%dims    => curr_dim
          ENDIF
          prev_dim => curr_dim
          curr_dim%file_dim_name      = TRIM(curr_opt%dim_rename(St:Split-3))
          curr_dim%dim_data%name_dim  = TRIM(curr_opt%dim_rename(Split:En-1))
          curr_dim%dim_data%sub_model = IFile%sub_model
          curr_dim%flags              = IOR(curr_dim%flags,INPUT_DIM_FILE)
        ENDIF
        St = En + 1
        DO WHILE ((St <= Le) .AND. (curr_opt%dim_rename(St:St)==' '))
          St = St + 1
        ENDDO
      ENDDO

      curr_opt => curr_opt%next
    ENDDO

    IFile%stat = IOR(IFile%stat,INPUT_FILE_EXTERN_APPLIED)

  END SUBROUTINE InputFileApplyExternalSettings

! Finds external settings applying to given file
  SUBROUTINE InputFileGetExternalSettings(IFile,opt,cnt)
  TYPE (input_file_list), POINTER :: IFile
  TYPE (input_opt_list),  POINTER :: opt
  INTEGER,            INTENT(out) :: cnt

  TYPE (input_opt_list),  POINTER :: curr, prev, last, prevopt
  TYPE (input_file_list), POINTER :: curr_file
  TYPE (input_var_list),  POINTER :: curr_var
  CHARACTER (len=128) :: GenName, OName, NName
  LOGICAL :: lCheck
  INTEGER :: st, en, nst, flags

    cnt = 0
    NULLIFY(opt,prev,last,prevopt)
    curr => AllInputOpts
    DO WHILE (ASSOCIATED(curr))
      lCheck = LEN_TRIM(curr%file_name) > 0
      IF (lCheck) THEN
        GenName = BuildGenericName(IFile)
        IF ((curr%sub_model/=IFile%sub_model) .AND. (curr%sub_model>0) .AND. (IFile%sub_model>0)) lCheck = .FALSE.
        IF (lCheck) THEN
          st = INDEX(curr%file_name,'=')
          IF ((st > 0) .AND. (curr%file_name(st+1:st+1)=='>')) THEN ! Renaming of a file is part of the specification
            st = st + 2
            en = INDEX(curr%file_name(st+2:),',') + st - 1
            IF (en < st) en = LEN_TRIM(curr%file_name) + 1
            IF (TRIM(GenName) == curr%file_name(st:en-1)) THEN  ! This is the file to be renamed
              nst = INDEX(curr%file_name(1:st),',',.TRUE.) + 1 ! curr%file_name(nst:st-2) is the new file name
              OName = TRIM(GenName)
              NName = curr%file_name(nst:st-3)
              curr%file_name = curr%file_name(1:st-3)//TRIM(curr%file_name(en:))
              ! Replace name in file specifications
              curr_file => AllInputFiles
              DO WHILE (ASSOCIATED(curr_file))
                GenName = BuildGenericName(curr_file)
                IF (TRIM(GenName) == TRIM(OName)) THEN
                  CALL InputFileSplitName(NName,curr_file%name_pre,curr_file%name_suf,flags)
                  curr_file%flags = IOR(IAND(curr_file%flags,INPUT_ALLFLAGS-INPUT_FILE_MULTI),flags)
                ENDIF
                curr_file => curr_file%next
              ENDDO
              ! Replace names in external setting specifications
              curr => AllInputOpts
              DO WHILE (ASSOCIATED(curr))
                IF (GetStringIndex(TRIM(OName),TRIM(curr%file_name),Sep=',',lCaseS=.TRUE.,St=st) > 0) THEN
                  en = INDEX(curr%file_name(st:),',') + st - 1
                  IF (en < st) en = LEN_TRIM(curr%file_name) + 1
                  curr%file_name = curr%file_name(1:st-1)//TRIM(NName)//TRIM(curr%file_name(en:))
                ENDIF
                curr => curr%next
              ENDDO
              ! Now test all external settings again
              curr => AllInputOpts
              CYCLE
            ENDIF
          ENDIF
          IF (GetStringIndex(TRIM(GenName),TRIM(curr%file_name),Sep=',',lCaseS=.TRUE.) <= 0) lCheck = .FALSE.
        ENDIF
      ENDIF
      ! Check if any file variables belongs to a group for which settings are specified
      IF (.NOT. lCheck) THEN
        curr_var => AllInputVars
        DO WHILE (ASSOCIATED(curr_var))
          IF (ASSOCIATED(curr_var%parent_file,IFile) .AND. ASSOCIATED(curr_var%group)) &
            lCheck = lCheck .OR. (GetStringIndex(curr_var%group%name_group,curr%var_name,Sep=',',lCaseS=.TRUE.) > 0)
          curr_var => curr_var%next
        ENDDO
      ENDIF
      IF (lCheck) THEN
        cnt = cnt + 1
        ! Swap elements, so that all regarding this file are contiguous
        IF (ASSOCIATED(last)) THEN
          IF (IAND(curr%flags,INPUT_OPT_NAMELIST) == 0) THEN ! Put specification objects first, namelist specs last
            prev%next => curr%next
            curr%next => opt
            IF (ASSOCIATED(prevopt)) prevopt%next => curr
            opt       => curr
            curr      => prev%next
          ELSE
            IF (.NOT. ASSOCIATED(last%next,curr)) THEN ! Are they already ordered?
              prev%next => curr%next ! Remove curr from its current position
              curr%next => last%next ! Forward link of new position
              last%next => curr      ! Link to new position from last found (insert curr)
              last      => curr      ! Last found
              curr      => prev%next ! Next options to test
            ELSE
              prev => curr      ! Dazzling pointer to be able to take out 
              curr => curr%next ! Already in order -> just go to next element
            ENDIF
          ENDIF
        ELSE ! No previous match
          prevopt => prev      ! Pointer to last uninteresting opt ahead of the interesting ones
          opt     => curr      ! First match should be returned
          last    => curr      ! Last in list of matching
          prev    => curr      ! Dazzling pointer to be able to take out 
          curr    => curr%next ! Go to next element
        ENDIF
      ELSE ! curr is not of interest
        prev => curr        ! Dazzling pointer to be able to take out 
        curr => curr%next   ! Go to next element
      ENDIF
    ENDDO

  END SUBROUTINE InputFileGetExternalSettings

! Add a new input file to the global list
  FUNCTION InputFileAdd(pre,suf,lyear,lmonth,lday,lhour,lmin,lsec,dt_file,dt_unit,year_offset,   &
  start_year,start_month,start_day,enddata_year,enddata_month,enddata_day, &
  labstime,time_error,cycle_act,day_offset,time_shift,sub_model,SetFlags)
  TYPE (input_file_list),        POINTER :: InputFileAdd 
  CHARACTER(len=*),           INTENT(in) :: pre
  CHARACTER(len=*), OPTIONAL, INTENT(in) :: suf
  LOGICAL, OPTIONAL,          INTENT(in) :: lyear, lmonth, lday, lhour, lmin, lsec
  INTEGER, OPTIONAL,          INTENT(in) :: dt_file, time_error, cycle_act, year_offset, day_offset, time_shift
  INTEGER, OPTIONAL,          INTENT(in) :: start_year,start_month,start_day,enddata_year,enddata_month,enddata_day
  INTEGER, OPTIONAL,          INTENT(in) :: sub_model, SetFlags
  LOGICAL, OPTIONAL,          INTENT(in) :: labstime
  CHARACTER(len=*), OPTIONAL, INTENT(in) :: dt_unit

#ifdef HAVE_F2003
  TYPE (input_file_type_fcns_list), POINTER :: ft
  CHARACTER (len=10) :: ext
  INTEGER :: i
#endif
  TYPE (input_file_list), POINTER :: new
  INTEGER :: flags
  LOGICAL :: lOk
  CHARACTER (len=64) :: name_pre, name_suf

    IF (.NOT. moInputInitialized) CALL local_error('InputFileAdd','Please call InputInit before adding files')
    IF (lInLoop) CALL local_error('InputFileAdd','Cannot add files after entering the main loop')

    ! Is this a dummy declaration?
    IF (TRIM(pre) == '') THEN
      IF (.NOT. ASSOCIATED(DummyFile)) CALL InputFileNew(DummyFile)
      InputFileAdd => DummyFile
      RETURN
    ENDIF

    ! Get effective name
    name_pre = pre
    name_suf = ''
    IF (PRESENT(suf)) name_suf = suf
    flags = 0

    IF (PRESENT(lyear))  THEN
      IF (lyear) flags  = IOR(flags,INPUT_FILE_NAME_YEAR) ! This option do not support 6 digit and/or signed years
    ENDIF
    IF (PRESENT(lmonth)) THEN
      IF (lmonth) flags = IOR(flags,INPUT_FILE_NAME_MONTH)
    ENDIF
    IF (PRESENT(lday))   THEN
      IF (lday) flags   = IOR(flags,INPUT_FILE_NAME_DAY)
    ENDIF
    IF (PRESENT(lhour))  THEN
      IF (lhour) flags  = IOR(flags,INPUT_FILE_NAME_HOUR)
    ENDIF
    IF (PRESENT(lmin))   THEN
      IF (lmin) flags   = IOR(flags,INPUT_FILE_NAME_MINUTE)
    ENDIF
    IF (PRESENT(lsec))   THEN
      IF (lsec) flags   = IOR(flags,INPUT_FILE_NAME_SECOND)
    ENDIF

    ! Extract date positions from file name
    IF ((flags==0) .AND. (LEN_TRIM(name_suf)==0)) CALL InputFileSplitName(pre,name_pre,name_suf,flags)

    ! Test if file already assigned
    new => InputFileGetRef(name_pre,name_suf,flags,sub_model)
    IF (ASSOCIATED(new)) THEN
      InputFileAdd => new
      RETURN
    ENDIF

    ! Create new file and add it to the end of the global list
    CALL InputFileNew(new)
    IF (ASSOCIATED(LastFile)) THEN
      LastFile%next => new
    ELSE
      AllInputFiles => new
    ENDIF
    LastFile => new

    new%name_pre = name_pre
    new%name_suf = name_suf
    new%flags    = flags
    new%nUsed    = 1

    ! Set time parameters
    IF (PRESENT(dt_file)) THEN
      IF (dt_file /= UNLIKELY_VAL) THEN
        IF (PRESENT(dt_unit)) THEN
          new%dt_file = InterpretDtUnit(dt_file,dt_unit)
        ELSE
          new%dt_file = dt_file
        ENDIF
        new%stat = IOR(new%stat,INPUT_FILE_TIME_GOT)
        IF (new%dt_file == 0) new%flags = IOR(new%flags,INPUT_FILE_INITIAL)
      ENDIF
    ENDIF
    IF (PRESENT(start_year   )) new%time_start(1:4)  = (/start_year,1,1,0/)
    IF (PRESENT(start_month  )) new%time_start(2)  = start_month
    IF (PRESENT(start_day    )) new%time_start(3)  = start_day
    IF (PRESENT(enddata_year )) new%time_enddata(1:4) = (/enddata_year,12,31,0/) ! hmmmmmmm
    IF (PRESENT(enddata_month)) new%time_enddata(2) = enddata_month
    IF (PRESENT(enddata_day  )) new%time_enddata(3) = enddata_day
    IF (PRESENT(year_offset  )) new%time_offset(1) = year_offset
    IF (PRESENT(day_offset   )) new%time_offset(2) = day_offset
    IF (PRESENT(time_shift   )) new%time_offset(3) = time_shift

    new%time_file = model_time
    new%time_file(1) = new%time_file(1) + new%time_offset(1)

    ! Prepare (and if possible: determine) start-of-repeat time
    IF (PRESENT(cycle_act)) THEN
      SELECT CASE (cycle_act)
        CASE (-2100000000:-1) ! Repeat (negative) number of months
          new%time_cycle(4) = -cycle_act
          new%time_cycle(1) = INPUT_FILE_REPEAT_MONTHS
        CASE (0:1799) ! Repeat number of records
          new%time_cycle(4) = cycle_act
          new%time_cycle(1) = INPUT_FILE_REPEAT_RECS
        CASE (1800:2100000000) ! Repeat number of seconds
          new%time_cycle(4) = cycle_act
          new%time_cycle(1) = INPUT_FILE_REPEAT_SECS
        CASE (INPUT_FILE_REPEAT_NONE) ! Don't repeat anything
          new%flags = IOR(new%flags,INPUT_FILE_NOCYCLE)
        CASE (INPUT_FILE_REPEAT_ALL) ! Repeat entire forcing
          new%time_cycle(1) = INPUT_FILE_REPEAT_ALL
        CASE (INPUT_FILE_REPEAT_YEAR) ! Repeat last year
          new%time_cycle(1) = INPUT_FILE_REPEAT_YEAR
        CASE (INPUT_FILE_REPEAT_DAY) ! Repeat last day
          new%time_cycle(1) = INPUT_FILE_REPEAT_DAY
      END SELECT
    ELSE
      IF (PRESENT(enddata_year)) THEN
        new%time_cycle(1:4) = new%time_start(1:4)
        new%time_cycle(5:6) = 0
      ELSE
        new%flags = IOR(new%flags,INPUT_FILE_NOCYCLE)
      ENDIF
    ENDIF

    IF (IAND(new%flags,INPUT_FILE_INITIAL)==0) THEN
      IF (PRESENT(labstime)) THEN
        lOK = .TRUE.
        IF (PRESENT(SetFlags)) THEN
          IF (IAND(SetFlags,1) == 0) lOK = .FALSE.
        ENDIF
        IF (lOK) THEN
          IF (labstime) THEN
            new%flags = IOR(new%flags,INPUT_FILE_TIME_ABSOLUTE)
          ELSE
            new%flags = IAND(new%flags,INPUT_ALLFLAGS-INPUT_FILE_TIME_ABSOLUTE)
          ENDIF
        ENDIF
      ELSE
        new%flags = IOR(new%flags,INPUT_FILE_TIME_ABSOLUTE)
      ENDIF
    ENDIF

    IF (PRESENT(time_error)) THEN
      new%flags = IOR(IAND(new%flags,INPUT_ALLFLAGS-INPUT_FILE_TIME_MASK),IAND(time_error,INPUT_FILE_TIME_MASK))
    ELSE
      IF (.NOT. PRESENT(SetFlags)) new%flags = IOR(new%flags,INPUT_FILE_TIME_ERR) ! Default : Time mismatch is an error
    ENDIF

    new%sub_model = CurrSubModel
    IF (PRESENT(sub_model)) new%sub_model = sub_model

    ! Determine file type if necessary
#ifdef HAVE_F2003
    IF (LEN_TRIM(name_suf) > 0) THEN
      i = INDEX(name_suf,'.',.TRUE.)
      IF (i>0) ext = TRIM(name_suf(i+1:))
    ELSE
      i = INDEX(name_pre,'.',.TRUE.)
      IF (i>0) ext = TRIM(name_pre(i+1:))
    ENDIF
    IF (i <= 0) CALL local_error('InputFileAdd','Could not determine file type')
    ft => AllInputFileTypes
    DO WHILE (ASSOCIATED(ft))
      IF (TRIM(ext) == TRIM(ft%Suffix)) THEN
        new%ft => ft
        NULLIFY(ft)
      ELSE
        ft => ft%next
      ENDIF
    ENDDO
    IF (.NOT. ASSOCIATED(new%ft)) CALL local_error('InputFileAdd','Unknown file type: '//TRIM(ext))
#endif

    InputFileAdd => new

  END FUNCTION InputFileAdd

! Find the model unlimited dimension in a file
  FUNCTION InputGetUnlimitedDim(IFile)
  TYPE (input_dim_list),  POINTER :: InputGetUnlimitedDim
  TYPE (input_file_list), POINTER :: IFile

  TYPE (input_dim_list),     POINTER :: curr_dim, file_dim
  TYPE (input_fileobj_list), POINTER :: obj

    ! Can we be sure that dimension and variable information is always available?
    curr_dim => AllInputDims
    DO WHILE (ASSOCIATED(curr_dim))
      IF (IAND(curr_dim%dim_data%flags,INPUT_DIM_UNLIM) /= 0) THEN
        IF (ANY(IFile%time_fid==0)) THEN
          IF (IAND(IFile%stat,INPUT_FILE_DIMS_GOT) == 0) &
            CALL InputFileOpen(IFile,IFile%time_start(1),IFile%time_start(2), &
                               IFile%time_start(3),check=INPUT_FILE_CHECK_DIM_SIZE)
          file_dim => InputDimGetRef(curr_dim%dim_data%name_dim,curr_dim%dim_data%name_alt,search=IFile%dims)
          IF (ASSOCIATED(file_dim)) THEN
            IFile%time_fid(1) = file_dim%fid
            obj => InputFileObjGetRef(IFile%fvars,file_dim%dim_data%name_dim,file_dim%file_dim_name)
            IF (ASSOCIATED(obj)) THEN
              IFile%time_fid(2) = obj%id
            ELSE
              IFile%time_fid(2) = -1
            ENDIF
          ELSE
            IFile%time_fid(1) = -1
          ENDIF
        ENDIF
        IF (IFile%time_fid(1) > 0) THEN
          InputGetUnlimitedDim => curr_dim
          RETURN
        ENDIF
      ENDIF
      curr_dim => curr_dim%next
    ENDDO

    NULLIFY(InputGetUnlimitedDim)

  END FUNCTION InputGetUnlimitedDim

! ----------------------------------------------------
!  Creating, adding, setting, modifying, displaying:
!     Variable groups 
! ----------------------------------------------------

! Creates a new group element with default values
  SUBROUTINE InputVarGroupNew(group)
  TYPE (input_group_list), POINTER :: group

  TYPE (input_group_list), POINTER :: new

    ALLOCATE(new)
    NULLIFY(new%next,new%group_file)
    new%flags      =  0
    new%dt_group   = UNLIKELY_VAL
    new%nIntrp     = UNLIKELY_VAL
    new%nFuture    = UNLIKELY_VAL
    new%nAtOnce    = UNLIKELY_VAL
    new%name_group = ''

    group => new

  END SUBROUTINE InputVarGroupNew

! Searches for a specific group and eventually create it if it doesn't exist
  FUNCTION InputVarGroupGetRef(group,lCreate,lOld)
  TYPE (input_group_list), POINTER :: InputVarGroupGetRef
  CHARACTER (len=*),    INTENT(in) :: group
  LOGICAL, OPTIONAL,    INTENT(in) :: lCreate
  LOGICAL, OPTIONAL,   INTENT(out) :: lOld

  TYPE (input_group_list), POINTER :: curr

    curr => AllInputGroups
    DO WHILE (ASSOCIATED(curr))
      IF (TRIM(group) == TRIM(curr%name_group)) THEN
        IF (PRESENT(lOld)) lOld = .TRUE.
        InputVarGroupGetRef => curr
        RETURN
      ENDIF
      curr => curr%next
    ENDDO

    NULLIFY(InputVarGroupGetRef)
    IF (PRESENT(lCreate)) THEN
      IF (lCreate) THEN
        CALL InputVarGroupNew(curr)
        curr%name_group = TRIM(group)
        IF (PRESENT(lOld)) lOld = .FALSE.
        InputVarGroupGetRef => curr
      ENDIF
    ENDIF

  END FUNCTION InputVarGroupGetRef

! Displays one or all groups in a list
  SUBROUTINE InputVarGroupDisp(group,cnt,unt,msg,list)
  TYPE (input_group_list), OPTIONAL, POINTER :: group
  INTEGER, OPTIONAL, INTENT(in) :: cnt
  INTEGER, OPTIONAL, INTENT(in) :: unt
  CHARACTER (len=*), INTENT(in), OPTIONAL :: msg
  LOGICAL, OPTIONAL, INTENT(in) :: list

  TYPE (input_group_list), POINTER :: curr
  TYPE (input_var_list),   POINTER :: var
  INTEGER :: un, cr

    un = GetLogUnit(); IF (PRESENT(unt)) un = unt
    cr = -1;           IF (PRESENT(cnt)) cr = cnt

    IF (PRESENT(group)) THEN
      curr => group
    ELSE
      curr => AllInputGroups
    ENDIF
    IF (PRESENT(list)) THEN
      IF (list) THEN
        message_text = ''
        DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
          message_text = TRIM(message_text)//', '//TRIM(curr%name_group)
          curr => curr%next
          cr = cr - 1
        ENDDO
        IF (PRESENT(msg)) THEN
          WRITE(un,'(2a)') TRIM(msg),TRIM(message_text(3:))
        ELSE
          WRITE(un,'(a)') TRIM(message_text(3:))
        ENDIF
      ENDIF
    ENDIF
    IF (PRESENT(msg)) WRITE(un,'(a)') msg
    DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
      WRITE(un,'(a)') TRIM(curr%name_group)
      message_text = ''
      var => AllInputVars
      DO WHILE (ASSOCIATED(var))
        IF (ASSOCIATED(var%group,curr)) message_text = TRIM(message_text)//', '//TRIM(var%name_var)
        var => var%next
      ENDDO
      WRITE(un,'(2a)') '  Group variables    : ',TRIM(message_text(3:))
      WRITE(un,'(2a)') '  Group time step    : ',TRIM(TimeStepString(curr%dt_group,UNLIKELY_VAL))
      WRITE(un,'(a,i2)') '  Group #Interpol    : ',curr%nIntrp
      WRITE(un,'(a,i2)') '  Group #Future steps: ',curr%nFuture
      WRITE(un,'(a,i2)') '  Group #Read at once: ',curr%nAtOnce
      WRITE(un,'(2a)') '  Group time validity: ',TRIM(GetStringIndexed(opt_time_val, &
                                                      IAND(curr%flags,INPUT_VAR_VALID_MASK)/INPUT_VAR_VALID_BEFORE+1))
      CALL DispFlagsSingleBit(curr%flags,INPUT_GROUP_MULTI_FILE,flg_group_n,opt_group,Txt=message_text)
      WRITE(un,'(2a)') '  Group properties   : ',TRIM(message_text)
      WRITE(un,'(a,Z8.8)') '  Group flags        : 0x',curr%flags
      IF (ASSOCIATED(curr%group_file)) &
        WRITE(un,'(2a)') '  Group file template: ',TRIM(BuildGenericName(curr%group_file))
      var => InputVarGetRef(TRIM(curr%name_group),list=AllInputGrpVars)
      IF (ASSOCIATED(var)) THEN
        WRITE(un,'(a)')  '  Group defaults in variable:'
        CALL InputVarDisp(var,1)
      ENDIF

      curr => curr%next
      cr = cr - 1
    ENDDO
 
  END SUBROUTINE InputVarGroupDisp

! Set the current group
  SUBROUTINE InputVarGroupStart(group,group_ref,IFile,file_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time, &
                                depend,depend_ref,lauto,nIntrp,nAtOnce,nFuture,subset_index,miss_val,weight,weights,sub_model)
  CHARACTER (len=*), INTENT(in) :: group
  TYPE (input_group_list), POINTER, OPTIONAL :: group_ref
  TYPE (input_file_list), OPTIONAL, POINTER    :: IFile   ! File to add variable to
  CHARACTER (LEN=*),      OPTIONAL, INTENT(in) :: file_name  ! Name of file to add variable to
  CHARACTER (LEN=*),      OPTIONAL, INTENT(in) :: depend     ! Comma separated list of names of vars on which this var depends
  TYPE (input_var_list),  OPTIONAL, POINTER    :: depend_ref ! Reference to variable on which this variable depends
  CHARACTER (len=*), OPTIONAL, INTENT(in) :: dt_unit      ! Unit of dt_update
  INTEGER,           OPTIONAL, INTENT(in) :: dt_update    ! Time between updates
  INTEGER,           OPTIONAL, INTENT(in) :: miss         ! What to do when var is missing in file : See structure flag desc 
  INTEGER,           OPTIONAL, INTENT(in) :: mis_dim      ! What to do when dimensions are mismatching between file and variable 
  INTEGER,           OPTIONAL, INTENT(in) :: data_action  ! How to combine file and model data 
  LOGICAL,           OPTIONAL, INTENT(in) :: lauto        ! Switch off auto-reset of READ and UPDATE flags a every InputUpdate call
  INTEGER,           OPTIONAL, INTENT(in) :: valid_time   ! When in interval are data valid
  INTEGER,           OPTIONAL, INTENT(in) :: nIntrp       ! Number of time steps to store at once for (user) interpolation
  INTEGER,           OPTIONAL, INTENT(in) :: nAtOnce      ! Number of time steps to read at once
  INTEGER,           OPTIONAL, INTENT(in) :: nFuture      ! Number of future time steps needed for (user) interpolation
  INTEGER,           OPTIONAL, INTENT(in) :: subset_index ! Index of missing dimension to set for partial reading
  INTEGER,           OPTIONAL, INTENT(in) :: sub_model    ! Sub model to which the variable should belong
  REAL(dp),          OPTIONAL, INTENT(in) :: miss_val     ! Value to replace bad/out-of-range values with and bad value
  REAL(dp),          OPTIONAL, INTENT(in) :: weight       ! Default value to use if var is missing and data weight for combinations
  REAL(dp),          OPTIONAL, INTENT(in) :: weights(:)   ! Weighting of sub-dimensions on contraction

  TYPE (input_group_list), POINTER :: curr
  TYPE (input_var_list),   POINTER :: grp_var
  INTEGER :: dummy(1)
  LOGICAL :: lOld

    curr => InputVarGroupGetRef(group,.TRUE.,lOld)
    IF (.NOT. lOld) THEN
      curr%next => AllInputGroups
      AllInputGroups => curr
    ENDIF

    IF (PRESENT(IFile)    .OR. PRESENT(file_name)   .OR. PRESENT(dt_update)  .OR. PRESENT(dt_unit) .OR. PRESENT(miss)         .OR.&
        PRESENT(mis_dim)  .OR. PRESENT(data_action) .OR. PRESENT(valid_time) .OR. PRESENT(depend)  .OR. PRESENT(depend_ref)   .OR.&
        PRESENT(lauto)    .OR. PRESENT(nIntrp)      .OR. PRESENT(nAtOnce)    .OR. PRESENT(nFuture) .OR. PRESENT(subset_index) .OR.&
        PRESENT(miss_val) .OR. PRESENT(weight)      .OR. PRESENT(weights)    .OR. PRESENT(sub_model)) THEN
      NULLIFY(CurrGroup) ! Don't add group to group variable
      grp_var => InputVarAdd_gen('',-2,dummy,lOld,var_name=group,IFile=IFile,file_name=file_name,dt_update=dt_update,           &
                                 dt_unit=dt_unit,miss=miss,mis_dim=mis_dim,data_action=data_action,valid_time=valid_time,       &
                                 depend=depend,depend_ref=depend_ref,lauto=lauto,nIntrp=nIntrp,nAtOnce=nAtOnce,nFuture=nFuture, &
                                 subset_index=subset_index,miss_val=miss_val,weight=weight,weights=weights,sub_model=sub_model)
      ! Add variable remove it from main list again and add to group-var list if its a new variable
      IF (.NOT. lOld) THEN
        AllInputVars    => AllInputVars%next
        grp_var%next    => AllInputGrpVars
        AllInputGrpVars => grp_var
        grp_var%stat    =  0 ! InputVarAdd_gen may have used %stat as should be done for normal variables, but not for group vars
      ENDIF
      ! Set variable flags
      grp_var%flags = IOR(grp_var%flags,INPUT_VAR_GROUP)
      IF (PRESENT(IFile)     .OR. PRESENT(file_name))  grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_IFILE)
      IF (PRESENT(dt_update) .OR. PRESENT(dt_unit))    grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_DT)
      IF (PRESENT(depend)    .OR. PRESENT(depend_ref)) grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_DEPEND)
      IF (PRESENT(miss))         grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_MISS_ACT)
      IF (PRESENT(mis_dim))      grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_MIS_DIM)
      IF (PRESENT(data_action))  grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_DATA_ACT)
      IF (PRESENT(valid_time))   grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_VALID_TIME)
      IF (PRESENT(lauto))        grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_LAUTO)
      IF (PRESENT(nIntrp))       grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_NINTRP)
      IF (PRESENT(nAtOnce))      grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_NATONCE)
      IF (PRESENT(nFuture))      grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_NFUTURE)
      IF (PRESENT(subset_index)) grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_INDEX)
      IF (PRESENT(miss_val))     grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_MISS_VAL)
      IF (PRESENT(weight))       grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_WEIGHT)
      IF (PRESENT(weights))      grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_WEIGHTS)
      IF (PRESENT(sub_model))    grp_var%stat = IOR(grp_var%stat,INPUT_GRPVAR_SET_SUB_MODEL)
    ENDIF

    CurrGroup => curr
    IF (PRESENT(group_ref)) group_ref => curr

  END SUBROUTINE InputVarGroupStart

! Reset current group
  SUBROUTINE InputVarGroupEnd()

    NULLIFY(CurrGroup)

  END SUBROUTINE InputVarGroupEnd

! ----------------------------------------------------
!  Creating, adding, setting, modifying, displaying:
!     Variables
! ----------------------------------------------------

! Creates a new input_var_list element with default values
  SUBROUTINE InputVarNew(var)
  TYPE (input_var_list), POINTER :: var

  TYPE (input_var_list), POINTER :: new 

    ALLOCATE(new)
    NULLIFY(new%next,new%dta,new%depend,new%masks,new%masks_intrp,new%Intrp,new%parent_file,new%saction, &
            new%dims,new%dims_at_read,new%buffers,new%buffers_intrp,new%aux_data,new%group)
    new%name_var      = ''
    new%name_alt      = ''
    new%sub_model     =  1
    new%ndim          =  0
    new%nIntrp        =  0
    new%nAtOnce       =  1
    new%nFuture       =  1
    new%edims(:)      = -1
    new%nupdate       =  0
    new%next_update   =  0
    new%nToRead       =  0
    new%fid           =  0
    new%flags         =  0
    new%stat          =  0
    new%nRead         =  0
    new%nLastRead     = -1
    new%nLastUpdate   = -1
    new%subset_index  = -1
    new%dt_update     =  0
    new%acp_intrp     =  0
    new%time(:)       = -HUGE(i8)
#ifdef HAVE_F2003
    NULLIFY(new%DataFcns,new%TimeIntrpFcn,new%AddActionFcn,new%ProcInterpol)
#endif
    var => new

  END SUBROUTINE InputVarNew

! Returns a default input_var_list element
  FUNCTION InputVarNewDummy()
  TYPE (input_var_list), POINTER :: InputVarNewDummy

    IF (.NOT. ASSOCIATED(DummyVar)) CALL InputVarNew(DummyVar)
    InputVarNewDummy => DummyVar

  END FUNCTION InputVarNewDummy

! Returns handle to a variable which matches the variable name and sub-model
  FUNCTION InputVarGetRef(var_name,alt_name,sub_model,list)
  TYPE (input_var_list),           POINTER    :: InputVarGetRef
  CHARACTER(len=*),                INTENT(in) :: var_name
  CHARACTER(len=*),      OPTIONAL, INTENT(in) :: alt_name
  INTEGER,               OPTIONAL, INTENT(in) :: sub_model
  TYPE (input_var_list), OPTIONAL, POINTER    :: list

  TYPE (input_var_list), POINTER :: curr
  INTEGER :: sm

    sm = CurrSubModel
    IF (PRESENT(sub_model)) sm = sub_model

    curr => AllInputVars
    IF (PRESENT(list)) curr => list
    DO WHILE (ASSOCIATED(curr))
      IF ((curr%sub_model < 0) .OR. (sm < 0) .OR. (curr%sub_model == sm)) THEN
        IF ((TRIM(var_name)==TRIM(curr%name_var)) .OR. &
            (TRIM(var_name)==TRIM(curr%name_alt))) THEN
          InputVarGetRef => curr
          RETURN
        ENDIF
        IF (PRESENT(alt_name)) THEN
          IF ((TRIM(alt_name)==TRIM(curr%name_var)) .OR. &
              (TRIM(alt_name)==TRIM(curr%name_alt))) THEN
            InputVarGetRef => curr
            RETURN
          ENDIF
        ENDIF
      ENDIF
    curr => curr%next
    ENDDO

    NULLIFY(InputVarGetRef)

  END FUNCTION InputVarGetRef

! Display variables in a list
  SUBROUTINE InputVarDisp(var,cnt,unt,msg,list)
  TYPE (input_var_list), OPTIONAL, POINTER :: var
  INTEGER, OPTIONAL, INTENT(in) :: cnt
  INTEGER, OPTIONAL, INTENT(in) :: unt
  CHARACTER (len=*), INTENT(in), OPTIONAL :: msg
  LOGICAL, OPTIONAL, INTENT(in) :: list

  TYPE (input_var_list),  POINTER :: curr
  TYPE (input_dim_list),  POINTER :: curr_dim, ref_dim
  TYPE (input_data_list), POINTER :: curr_buf, glob_buf
  TYPE (input_mask_list), POINTER :: curr_mask
  TYPE (input_file_list), POINTER :: curr_file
  CHARACTER (len=256) :: tmpstr
  INTEGER :: un, i, j, st, en, tmp(7), cr
  LOGICAL :: lCont
#ifdef HAVE_F2003
  TYPE (input_curr)               :: curr_proc
  TYPE (input_fcn_list),  POINTER :: fcn
#endif

    un = GetLogUnit(); IF (PRESENT(unt)) un = unt
    cr = -1;           IF (PRESENT(cnt)) cr = cnt

    IF (PRESENT(var)) THEN
      curr => var
    ELSE
      curr => AllInputVars
    ENDIF
    IF (PRESENT(list)) THEN
      IF (list) THEN
        tmpstr = ''
        DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
          tmpstr = TRIM(tmpstr)//','//TRIM(curr%name_var)
          curr => curr%next
          cr = cr - 1
        ENDDO
        IF (PRESENT(msg)) THEN
          WRITE(un,'(2a)') msg,TRIM(tmpstr(2:))
        ELSE
          WRITE(un,'(a)') TRIM(tmpstr(2:))
        ENDIF
        RETURN
      ENDIF
    ENDIF
    IF (PRESENT(msg)) WRITE(un,'(a)') msg
    DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
      IF (LEN_TRIM(curr%name_alt) /= 0) THEN
        WRITE(message_text,*) TRIM(curr%name_var),'%(',TRIM(curr%name_alt),',',curr%sub_model,')'
      ELSE
        WRITE(message_text,*) TRIM(curr%name_var),'%(',curr%sub_model,')'
      ENDIF
      CALL RmSpc(message_text,'%')
      WRITE(un,'(a)') TRIM(message_text)
      IF (ASSOCIATED(curr%group))  WRITE(un,'(2a)') '  Variable group                : ',TRIM(curr%group%name_group)
      IF (ASSOCIATED(curr%depend)) WRITE(un,'(2a)') '  Dependent on variable         : ',TRIM(curr%depend%name_var)

      message_text = ''
      curr_dim => curr%dims
      DO WHILE (ASSOCIATED(curr_dim))
        message_text = TRIM(message_text)//','//TRIM(curr_dim%dim_data%name_dim)
        curr_dim => curr_dim%next
      ENDDO
      WRITE(tmpstr,*) curr%ndim,' (',TRIM(message_text(2:)),')'
      CALL RmDblSpc(tmpstr)
      WRITE(un,'(2a)') '  Dimensions                    :',TRIM(tmpstr)
      IF (liodata .AND. ASSOCIATED(curr%parent_file)) THEN
        i  = 0
        IF ((IAND(curr%flags,INPUT_VAR_GROUP) == 0) .AND. (IAND(curr%stat,INPUT_VAR_MULTI_SPREAD) /= 0)) THEN
          message_text = ''
          curr_file => curr%parent_file
          st = 1
          DO WHILE (st < LEN_TRIM(curr%name_alt))
            en = INDEX(curr%name_alt(st:),',') + st - 1
            IF (en < st) en = LEN_TRIM(curr%name_alt) + 1
            message_text = TRIM(message_text)//','//TRIM(curr_file%name_file)
            curr_file => curr_file%next
            i  = i  + 1
            st = en + 1
          ENDDO
        ELSE
          message_text = ','//TRIM(curr%parent_file%name_file)
        ENDIF
        WRITE(un,'(2a)') '  Source file(s)                : ',TRIM(message_text(2:))
        IF (IAND(curr%flags,INPUT_VAR_GROUP) == 0) THEN
          IF ((i > 0) .AND. ASSOCIATED(curr%saction)) THEN
            WRITE(message_text,*) ' Variable file identifications :',curr%saction(3:2+i) ! Multiread is ALWAYS the first action
          ELSE
            WRITE(message_text,*) ' Variable file identification  :',curr%fid
          ENDIF
          CALL RmDblSpc(message_text(32:))
          WRITE(un,'(2a)') TRIM(message_text)
          message_text = ''
          curr_dim => curr%dims_at_read
          DO WHILE (ASSOCIATED(curr_dim))
            lCont = .TRUE.
            NULLIFY(ref_dim)
            IF (IAND(InputDimGetRefEx(curr_dim,curr%dims,ref=ref_dim,no=i),INPUT_DIM_FND_MASK) == INPUT_DIM_FND_PACKED) THEN
              ref_dim => ref_dim%src_dims
              DO j=1,i-1
                ref_dim => ref_dim%next
              ENDDO
              IF (ref_dim%lo > 0) THEN
                lCont = .FALSE.
                WRITE(message_text,*) TRIM(message_text),'%',TRIM(curr_dim%dim_data%name_dim),&
                                      '%(',ref_dim%lo,'-',ref_dim%hi,'),'
              ENDIF
            ENDIF          
            IF (lCont) THEN
              WRITE(message_text,*) TRIM(message_text),'%',TRIM(curr_dim%dim_data%name_dim),&
                                    '%(',curr_dim%dim_data%chunk_lo,'-',curr_dim%dim_data%chunk_hi,'),'
            ENDIF
            curr_dim => curr_dim%next
          ENDDO        
          CALL RmSpc(message_text,'%')
          WRITE(un,'(2a)') '  Reading      size             :',message_text(1:LEN_TRIM(message_text)-1)
        ENDIF
      ENDIF
      IF (IAND(curr%flags,INPUT_VAR_GROUP) == 0) THEN
        WRITE(message_text,*) curr%edims(1:curr%ndim)
        CALL RmDblSpc(message_text)
        WRITE(un,'(2a)') '  Expected mem size             :',TRIM(message_text)
      ENDIF

      IF (ASSOCIATED(curr%dta)) THEN
        WRITE(message_text,*) SHAPE(curr%dta)
        CALL RmDblSpc(message_text)
        WRITE(un,'(2a)') '  Model data size               :',TRIM(message_text)
      ENDIF
      IF (ASSOCIATED(curr%Intrp)) THEN
        WRITE(un,'(a,4L2)') '  Data associated (Model, Intrp1, Intrp2, Masks):',&
        ASSOCIATED(curr%dta),ASSOCIATED(curr%Intrp%next%dta),ASSOCIATED(curr%Intrp%dta), &
        ASSOCIATED(curr%masks)
      ELSE
        WRITE(un,'(a,2L2)') '  Data associated (Model, Masks):',ASSOCIATED(curr%dta),ASSOCIATED(curr%masks)
      ENDIF
      IF (ASSOCIATED(curr%masks)) THEN
        curr_mask => curr%masks
        i = 0
        message_text = '___Masks_(size_(valid_points))___:_'
        DO WHILE (ASSOCIATED(curr_mask))
          i = i + 1
          WRITE(message_text,*) TRIM(message_text(2:)),SIZE(curr_mask%mask),'_(',curr_mask%nvalid,'),_'
          IF (LEN_TRIM(message_text) > 170) THEN
            CALL RmSpc(message_text)
            WRITE(un,'(a)') message_text(1:LEN_TRIM(message_text)-1)
            message_text = '_________________________________:_'
          ENDIF
          curr_mask => curr_mask%next
        ENDDO
        CALL RmSpc(message_text)
        IF (LEN_TRIM(message_text) > 35) WRITE(un,'(a)') message_text(1:LEN_TRIM(message_text)-1)
      ENDIF
      WRITE(un,'(2a)') '  Variable time                 : ',DateString(curr%time)
      WRITE(un,'(a,i2)') '  #Interpol. fields             :',curr%nIntrp
      WRITE(un,'(a,i2)') '  #File read at once            :',curr%nAtOnce
      WRITE(un,'(a,i2)') '  #Future fields                :',curr%nFuture
      IF (ASSOCIATED(curr%Intrp)) THEN
        IF (ASSOCIATED(curr%Intrp%dta)) THEN
          WRITE(message_text,*) SHAPE(curr%Intrp%dta)
          CALL RmDblSpc(message_text)
          WRITE(un,'(2a)') '  Intrp data size               :',TRIM(message_text)
        ENDIF
        curr_buf => curr%Intrp%next
        message_text = ''
        DO i = 1,curr%nIntrp
          message_text = TRIM(message_text)//', '//DateString(curr_buf%time)
          curr_buf => curr_buf%next
        ENDDO
        CALL RmDblSpc(message_text)
        WRITE(un,'(3a)') '  Time of interpolation fields  : ',TRIM(message_text(3:))
      ENDIF
      message_text(1:11) = TimeStepString(curr%dt_update,-HUGE(i))
      WRITE(message_text(12:),*) '_(',curr%nupdate(1)
      CALL RmSpc(message_text)
      WRITE(un,'(3a)') '  Time between updates          : ',TRIM(message_text),' time steps)'
      IF (IAND(curr%flags,INPUT_VAR_GROUP) == 0) THEN
        WRITE(message_text,*) curr%next_update
        CALL RmDblSpc(message_text)
        WRITE(un,'(2a)') '  Next update in (time steps)   :',TRIM(message_text)
        WRITE(message_text,*) curr%nLastUpdate
        CALL RmDblSpc(message_text)
        WRITE(un,'(2a)') '  Steps since last update       :',TRIM(message_text)
        WRITE(message_text,*) curr%nLastRead
        CALL RmDblSpc(message_text)
        WRITE(un,'(2a)') '  Steps since last read         :',TRIM(message_text)
      ENDIF
      WRITE(un,'(2a)') '  Relative time                 : ',TRIM(GetStringIndexed(opt_time_val, &
                                                            IAND(curr%flags,INPUT_VAR_VALID_MASK)/INPUT_VAR_VALID_BEFORE+1))
      WRITE(un,'(2(a,Z8.8))') '  Flags (flags,status)          : 0x',curr%flags,', 0x',curr%stat
      IF (IAND(curr%flags,INPUT_VAR_GROUP) == 0) THEN
        CALL DispFlagsSingleBit(curr%stat,INPUT_VAR_INITIALIZED,flg_var_stat_n,flg_var_stat,Txt=message_text)
        WRITE(un,'(2a)') '  Variable status               : ',TRIM(message_text)
      ENDIF
      message_text = ''
      IF (IAND(curr%flags,INPUT_VAR_DIV_MASK) /= 0) message_text = ' '//TRIM(GetStringIndexed(opt_timediv, &
                                                                  IAND(curr%flags,INPUT_VAR_DIV_MASK)/INPUT_VAR_DIV_MONTH))//','
      WRITE(un,'(2a)') '  Model/read data merging       :',TRIM(message_text)//' '//TRIM(GetStringIndexed(opt_data, &
                                                            IAND(curr%flags,INPUT_VAR_ACTION_MASK)/INPUT_VAR_ACTION_ADD+1))
      WRITE(un,'(2a)') '  Missing action                : ',TRIM(GetStringIndexed(opt_var_miss, &
                                                            IAND(curr%flags,INPUT_VAR_MISS_MASK)/INPUT_VAR_MISS_IGNORE))
      WRITE(un,'(2a)') '  Dimension mismatch action     : ',TRIM(GetStringIndexed(opt_dim_mis, &
                                                     IAND(curr%flags,INPUT_VAR_DIM_MIS_MASK)/INPUT_VAR_DIM_MIS_WARN   +1,sep=','))
      WRITE(un,'(2a)') '  Out of range data action      : ',TRIM(GetStringIndexed(opt_bad_val, &
                                                     IAND(curr%flags,INPUT_VAR_BAD_VAL_MASK)/INPUT_VAR_BAD_VAL_REPLACE+1,sep=' '))
      WRITE(message_text,*) curr%subset_index
      CALL RmDblSpc(message_text)
      WRITE(un,'(2a)') '  Index for subset reading      :',TRIM(message_text)
      IF (ASSOCIATED(curr%saction)) THEN
        message_text = ''
#ifdef HAVE_F2003
        curr_proc%acp  =  1
        curr_proc%buf  => curr_buf
        curr_proc%mask => curr_mask
        fcn => curr%DataFcns
        lCont = ASSOCIATED(fcn)
        DO WHILE (lCont)
          st = curr_proc%acp
          CALL fcn%fcn(INPUT_PROC_NAME+INPUT_PROC_PROGRESS,curr,curr_proc,tmpstr(1:32))
          message_text = TRIM(message_text) // ' '// TRIM(tmpstr(1:32))
          fcn => fcn%next
          lCont = ASSOCIATED(fcn)
          IF (lCont) THEN
            ! Data to be pushed manually add a NoAction to prevent further immediate processing - still display
            IF ((st==1) .AND. (curr%saction(1)==INPUT_VAR_DO_NOACTION)) curr_proc%acp = 2
            lCont = (curr_proc%acp <= SIZE(curr%saction))
          ENDIF
        ENDDO
        j = curr_proc%acp
#else
        j = 0
        i = 1
        DO WHILE (i<=SIZE(curr%saction))
          IF (.NOT. liodata .AND. ((curr%saction(i)==INPUT_VAR_DO_BROADCAST) .OR. (curr%saction(i)==INPUT_VAR_DO_SCATTER))) THEN
            message_text = TRIM(message_text)//' recieve'
          ELSE
            message_text = TRIM(message_text)//' '//TRIM(GetStringIndexed(flg_var_spatial,curr%saction(i)+1))
          ENDIF
          SELECT CASE (ABS(curr%saction(i)))
          CASE (INPUT_VAR_DO_NOACTION)
            j = i
            i = SIZE(curr%saction)+1
          CASE (INPUT_VAR_DO_TOREAD)
            i = i + 2 * curr%saction(i+2) + 3
          CASE (INPUT_VAR_DO_MULTI_READ)
            i = i + 2 + curr%saction(i+1)
          CASE (INPUT_VAR_DO_PERMUTE)
            i = i + 2 * curr%saction(i+1) + 2
          CASE (INPUT_VAR_DO_REVERSE)
            i = i + 3 + curr%saction(i+1)
          CASE (INPUT_VAR_DO_SPREAD)
            i = i + 3 + curr%saction(i+1) + 2*curr%saction(i+2+curr%saction(i+1))
          CASE (INPUT_VAR_DO_USER_ACTION)
            i = i + curr%saction(i+1)
          CASE (INPUT_VAR_DO_SWAP)
            i = i + 4
          CASE (INPUT_VAR_DO_SUM,INPUT_VAR_DO_AVG,INPUT_VAR_DO_NORM,INPUT_VAR_DO_UNPACK,INPUT_VAR_DO_PACK,INPUT_VAR_DO_EQUIV)
            i = i + 3
          CASE (INPUT_VAR_DO_RESCALE)
            i = i + 2
          CASE (INPUT_VAR_DO_SUBSET,INPUT_VAR_DO_INTERPOLATE,INPUT_VAR_DO_BROADCAST,INPUT_VAR_DO_SCATTER,INPUT_VAR_DO_COPY, &
                INPUT_VAR_DO_COMBINE_ADD, INPUT_VAR_DO_COMBINE_SUB,  INPUT_VAR_DO_COMBINE_SUBR, INPUT_VAR_DO_COMBINE_MUL,   &
                INPUT_VAR_DO_COMBINE_DIV, INPUT_VAR_DO_COMBINE_DIVR, INPUT_VAR_DO_COMBINE_ADD_FLUX,                         &
                INPUT_VAR_DO_COMBINE_NUDGE, INPUT_VAR_DO_COMBINE_NUDGE_FLD, INPUT_VAR_DO_FILL_BAD, INPUT_VAR_DO_FILL_RANGE)
            i = i + 1
          END SELECT
        ENDDO
        IF (j==0) j = i
#endif
        WRITE(un,'(2a)') '  Process chain                 :',TRIM(message_text)
        WRITE(message_text,*) curr%saction(1:j-1)
        CALL RmDblSpc(message_text)
        WRITE(un,'(2a)') '  Process chain (details)       :',TRIM(message_text)
      ENDIF
      IF (IAND(curr%flags,INPUT_VAR_GROUP) == 0) THEN
        message_text = 'Select case'
#ifdef HAVE_F2003
        IF (ASSOCIATED(curr%DataFcns)) message_text = 'Procedure pointers'
#endif
        WRITE(un,'(2a)') '  Processing method             : ',TRIM(message_text)
#ifdef HAVE_F2003
        WRITE(un,'(a,l1)') '  Add action function associated: ',ASSOCIATED(curr%AddActionFcn)
#endif
      ENDIF
      CALL DispFlagsSingleBit(curr%flags,INPUT_VAR_INTERPOL_USER,flg_var_type_n,flg_var_type,Txt=message_text)
      WRITE(un,'(2a)') '  Variable type                 : ',TRIM(message_text)

      IF (IAND(curr%flags,INPUT_VAR_GROUP) == 0) THEN
        WRITE(message_text,*) curr%nRead
        CALL RmDblSpc(message_text)
        WRITE(un,'(2a)') '  Records read                  :',TRIM(message_text)
      ELSE
        CALL DispFlagsSingleBit(curr%stat,INPUT_GRPVAR_SET_IFILE,flg_grpvar_stat_n,flg_grpvar_stat,Txt=message_text)
        WRITE(un,'(2a)') '  Group parameters set          : ',TRIM(message_text)
      ENDIF
      IF (ASSOCIATED(curr%aux_data)) THEN
        WRITE(un,'(a)') '  Auxillary real data           : '
        curr_buf => curr%aux_data
        DO WHILE (ASSOCIATED(curr_buf))
          i = IAND(curr_buf%flags,INPUT_DATA_TYPE_MASK)
          CALL DispFlagsSingleBit(curr_buf%flags,INPUT_DATA_ALLOC,flg_data_n,flg_data,Txt=message_text)
          st = LEN_TRIM(message_text)
          IF ((i==INPUT_DATA_USER) .OR. (i==INPUT_DATA_RESCALE) .OR. (i==INPUT_DATA_FILL) .OR. (i==INPUT_DATA_RANGE)) THEN
            WRITE(message_text(st+1:),*) IAND(curr_buf%flags,INPUT_DATA_NO_MASK)/INPUT_DATA_NO_STEP
            CALL RmDblSpc(message_text)
          ENDIF
          message_text(29:29) = ':'
          IF (ASSOCIATED(curr_buf%dta)) THEN
            tmp = SHAPE(curr_buf%dta)
            IF (ANY(tmp(2:6) > 1)) THEN
              WRITE(tmpstr,'(a,7i4,3(a,G10.4))')  'Size: ',                                tmp,' min: ',MINVAL(curr_buf%dta), &
                                                 ' mean: ',SUM(curr_buf%dta)/SIZE(curr_buf%dta),'max: ',MAXVAL(curr_buf%dta)
              CALL RmDblSpc(tmpstr)
            ELSE
              tmpstr = ''
              IF (tmp(1) > 10) THEN
                WRITE(tmpstr,*) 'Length: ',tmp(1),' :'
                CALL RmDblSpc(tmpstr)
              ENDIF
              DO i=1,MIN(tmp(1),10)
                st = LEN_TRIM(tmpstr)
                WRITE(tmpstr(st+2:),'(G10.4)') curr_buf%dta(i,1,1,1,1,1,1)
              ENDDO
              CALL RmDblSpc(tmpstr)
              IF (tmp(1) > 10) tmpstr = TRIM(tmpstr)//' ...'
            ENDIF
          ELSE
            tmpstr = '<Unassociated buffer>'
          ENDIF
          WRITE(un,'(4a)') '    ',TRIM(message_text),' ',TRIM(tmpstr)
          curr_buf => curr_buf%next
        ENDDO
      ENDIF
      IF (ASSOCIATED(curr%buffers)) THEN
        message_text = '  Tmp buffers (-ve = non shared): '
        curr_buf => curr%buffers
        DO WHILE (ASSOCIATED(curr_buf))
          glob_buf => InputDataGetRefGlobal(ref=curr_buf,no=i)
          IF (ASSOCIATED(glob_buf)) i = i ! Just to get rid of a silly compiler warning
          IF (ASSOCIATED(curr_buf%dta)) THEN
            WRITE(message_text(35:),'(a,i3,a,i6,6i5)') '(',i,')',SHAPE(curr_buf%dta)
          ELSE
            WRITE(message_text(35:),'(a,i6,6i5,a,i2,a)') '( -1)',curr%edims,' (',curr%nIntrp,' x )'
          ENDIF
          WRITE(un,'(a)') TRIM(message_text)
          message_text = ''
          curr_buf => curr_buf%next
        ENDDO
      ENDIF

      curr => curr%next
      cr = cr - 1
    ENDDO

  END SUBROUTINE InputVarDisp

! Returns various information on a variable
  SUBROUTINE InputVarInq(name_var,ref,fid,name_file,name_depend,name_group,name_alt,dt_update,nTimes,nAtOnce,nFuture, &
                         state_update,state_read)
  CHARACTER (len=*),     OPTIONAL, INTENT(in)  :: name_var
  TYPE (input_var_list), OPTIONAL, POINTER     :: ref
  INTEGER,               OPTIONAL, INTENT(out) :: fid
  CHARACTER (len=64),    OPTIONAL, INTENT(out) :: name_file, name_depend, name_group, name_alt
  INTEGER,               OPTIONAL, INTENT(out) :: dt_update, state_update, state_read
  INTEGER,               OPTIONAL, INTENT(out) :: nTimes, nAtOnce, nFuture

  TYPE (input_var_list), POINTER :: curr
  INTEGER :: stat, tme(4)

    IF (PRESENT(name_var)) THEN
      curr => InputVarGetRef(name_var)
      IF (PRESENT(ref)) THEN
        IF (ASSOCIATED(ref) .AND. .NOT. ASSOCIATED(ref,curr)) &
          CALL local_error('InputVarInq','Name and reference refere to different dimensions')
        ref => curr
      ENDIF
    ELSE
      IF (.NOT. PRESENT(ref))    CALL local_error('InputVarInq','Neither name nor reference specified')
      IF (.NOT. ASSOCIATED(ref)) CALL local_error('InputVarInq','Neither name nor reference specified')
      curr => ref
    ENDIF

    IF (PRESENT(name_depend)) THEN
      IF (ASSOCIATED(curr%depend)) THEN
        name_depend = TRIM(curr%depend%name_var)
      ELSE
        name_depend = ''
      ENDIF
    ENDIF

    IF (PRESENT(name_file)) THEN
      IF (ASSOCIATED(curr%parent_file)) THEN
        name_file = TRIM(BuildGenericName(curr%parent_file))
      ELSE
        name_file = ''
      ENDIF
    ENDIF

    IF (PRESENT(name_group)) THEN
      IF (ASSOCIATED(curr%group)) THEN
        name_group = TRIM(curr%group%name_group)
      ELSE
        name_group = ''
      ENDIF
    ENDIF

    IF (PRESENT(fid))       fid       = curr%fid
    IF (PRESENT(dt_update)) dt_update = curr%dt_update
    IF (PRESENT(nTimes))    nTimes    = curr%nIntrp
    IF (PRESENT(nFuture))   nFuture   = curr%nFuture
    IF (PRESENT(nAtOnce))   nAtOnce   = curr%nAtOnce
    IF (PRESENT(name_alt))  name_alt  = curr%name_alt

    IF (PRESENT(state_update)) THEN
      stat = 0
      IF (IAND(curr%stat,INPUT_VAR_UPDATED) /= 0) stat = 1
      IF (curr%next_update <= 0)           stat = stat + 2
      state_update = stat
    ENDIF
    IF (PRESENT(state_read)) THEN
      stat = 0
      IF (IAND(curr%stat,INPUT_VAR_READ) /= 0) stat = 1
      IF (curr%next_update <= 0) THEN
        tme = model_time
        CALL CalendarTimeAdd(tme,INT(model_dt,i8))
        IF (CalendarTimeDiff(curr%time,tme) >= 0) stat = stat + 2
      ENDIF
      state_read = stat
    ENDIF

  END SUBROUTINE InputVarInq

! Extracts parameters from external settings and apply them to a variable 
  SUBROUTINE InputVarApplyExternalSettings(var,opt,cnt)
  TYPE (input_var_list), POINTER :: var
  TYPE (input_opt_list), POINTER :: opt
  INTEGER,            INTENT(in) :: cnt

  TYPE (input_opt_list),  POINTER :: curr_opt, file_def
  TYPE (input_file_list), POINTER :: new_file, prev_file
  REAL(dp),               POINTER :: aux(:), tst(:,:,:,:,:,:,:)
  REAL(dp)                        :: add, mul, tmp
  CHARACTER (len=128) :: dt_file_unit
  INTEGER :: i, j, dt_file, st, en
  LOGICAL :: lCont

    NULLIFY(file_def)
    curr_opt => opt
    DO j=1,cnt
      ! Set those variable properties which are specified
      IF (LEN_TRIM(curr_opt%weights) > 0) THEN
        aux => String2Vec(curr_opt%weights)
        CALL ReplaceAuxData1D(var,INPUT_DATA_WEIGHTS,0,aux)
        DEALLOCATE(aux)
      ENDIF
      IF (LEN_TRIM(curr_opt%var_file_name) > 0) THEN
        var%name_alt = TRIM(curr_opt%var_file_name)
        IF ((INDEX(var%name_alt,',') > 0) .OR. (INDEX(var%name_alt,'*') > 0) .OR. (INDEX(var%name_alt,'?') > 0)) &
          var%stat = IOR(var%stat,INPUT_VAR_MULTI) ! Multiple names given => multi-variable
      ENDIF
      mul = 1._dp
      add = 0._dp
      IF (curr_opt%mul           /= UNLIKELY_RVAL) mul = curr_opt%mul
      IF (curr_opt%add           /= UNLIKELY_RVAL) add = curr_opt%add
      IF ((mul/=1._dp) .OR. (add /= 0._dp)) CALL ReplaceAuxData1D(var,INPUT_DATA_RESCALE,INPUT_DATA_NO_STEP,(/mul,add/))
      IF (curr_opt%def_val       /= UNLIKELY_RVAL) THEN
        IF (IAND(var%flags,INPUT_VAR_MISS_MASK) == INPUT_VAR_MISS_DEFAULT) &
          var%flags = IOR(var%flags,INPUT_VAR_MISS_WARN)
        CALL ReplaceAuxData1D(var,INPUT_DATA_DEFAULT,0,(/curr_opt%def_val/))
      ENDIF
      IF (LEN_TRIM(curr_opt%def_vec) > 0) THEN
        IF (IAND(var%flags,INPUT_VAR_MISS_MASK) == INPUT_VAR_MISS_DEFAULT) &
          var%flags = IOR(var%flags,INPUT_VAR_MISS_WARN)
        aux => String2Vec(curr_opt%def_vec)
        CALL ReplaceAuxData1D(var,INPUT_DATA_DEFAULT,0,aux)
        DEALLOCATE(aux)
      ENDIF
      IF (curr_opt%fill_val /= UNLIKELY_RVAL) THEN
        tmp = UNLIKELY_RVAL
        tst => InputDataGetSpecified(var%aux_data,INPUT_DATA_FILL,INPUT_DATA_NO_STEP)
        IF (ASSOCIATED(tst)) tmp = tst(2,1,1,1,1,1,1)
        CALL ReplaceAuxData1D(var,INPUT_DATA_FILL,INPUT_DATA_NO_STEP,(/curr_opt%fill_val,tmp/))
      ENDIF
      IF (curr_opt%missing_val /= UNLIKELY_RVAL) THEN
        tmp = UNLIKELY_RVAL
        tst => InputDataGetSpecified(var%aux_data,INPUT_DATA_FILL,INPUT_DATA_NO_STEP)
        IF (ASSOCIATED(tst)) THEN
          tmp = tst(1,1,1,1,1,1,1)
        ELSE
          tst => InputDataGetSpecified(var%aux_data,INPUT_DATA_DEFAULT,0)
          IF (ASSOCIATED(tst)) THEN
            IF (SIZE(tst) == 1) tmp = tst(1,1,1,1,1,1,1)
          ENDIF
        ENDIF
        CALL ReplaceAuxData1D(var,INPUT_DATA_FILL,INPUT_DATA_NO_STEP,(/tmp,curr_opt%missing_val/))
      ENDIF
      IF (ANY(curr_opt%valid_range /= UNLIKELY_RVAL)) CALL ReplaceAuxData1D(var,INPUT_DATA_RANGE,INPUT_DATA_NO_STEP, &
                                                        (/UNLIKELY_RVAL,curr_opt%valid_range(1),curr_opt%valid_range(2)/))
      IF (curr_opt%file_weight   /= UNLIKELY_RVAL) &
        CALL ReplaceAuxData1D(var,INPUT_DATA_NUDGE_SCA,0,(/curr_opt%file_weight/))
      IF (curr_opt%subset_index  /=      UNLIKELY_VAL                ) var%subset_index = curr_opt%subset_index
      IF ((curr_opt%init_rec  /= UNLIKELY_VAL) .AND. ASSOCIATED(var%parent_file)) var%parent_file%init_rec = curr_opt%init_rec
      IF ((curr_opt%dt_update /= UNLIKELY_VAL) .OR. (LEN_TRIM(curr_opt%dt_update_unit) > 0)) THEN
        var%dt_update = InterpretDtUnit(curr_opt%dt_update,curr_opt%dt_update_unit,INPUT_VAR_FILE_TIME_STEP)
        IF (var%dt_update /= INPUT_VAR_FILE_TIME_STEP) var%stat = IOR(var%stat,INPUT_VAR_DT_SET)
      ENDIF
      IF (LEN_TRIM(curr_opt%action_var_miss) > 0) THEN
        i = GetStringIndex(curr_opt%action_var_miss,opt_var_miss,ErrMsg='action_var_miss')
        IF (i > 0) var%flags = IOR(IAND(var%flags,INPUT_ALLFLAGS-INPUT_VAR_MISS_MASK),i*INPUT_VAR_MISS_IGNORE)
      ENDIF
      IF (LEN_TRIM(curr_opt%action_dim_mis) > 0) &
        var%flags = IOR(IAND(var%flags,INPUT_ALLFLAGS - INPUT_VAR_DIM_MIS_MASK - INPUT_VAR_DIM_MIS_WARN), &
                    GetFlagsSingleBit(curr_opt%action_dim_mis,INPUT_VAR_DIM_MIS_WARN,opt_dim_mis,.TRUE.,'action_dim_mis',ofs=1))
      IF (LEN_TRIM(curr_opt%action_bad_val) > 0) THEN
        i = GetStringIndex(curr_opt%action_bad_val,opt_bad_val,ErrMsg='action_bad_val')-1
        IF (i > 0) var%flags = IOR(IAND(var%flags,INPUT_ALLFLAGS-INPUT_VAR_BAD_VAL_MASK),i*INPUT_VAR_BAD_VAL_REPLACE)
      ENDIF
      st = 1
      lCont = .TRUE.
      DO WHILE (lCont)
        en = INDEX(curr_opt%action_data(st:),',') + st - 1
        IF (en <  st) THEN
          en = LEN_TRIM(curr_opt%action_data) + 1
          lCont = .FALSE.
        ENDIF
        i = GetStringIndex(curr_opt%action_data(st:en-1),opt_data)-1
        IF (i >= 0) THEN
          var%flags = IOR(IAND(var%flags,INPUT_ALLFLAGS-INPUT_VAR_ACTION_MASK),i*INPUT_VAR_ACTION_ADD)
        ELSE
          i = GetStringIndex(curr_opt%action_data(st:en-1),opt_timediv,ErrMsg='action_data')
          var%flags = IOR(IAND(var%flags,INPUT_ALLFLAGS-INPUT_VAR_DIV_MASK),i*INPUT_VAR_DIV_MONTH)
        ENDIF
        st = en + 1
      ENDDO
      i = GetStringIndex(curr_opt%valid_time,opt_time_val,ErrMsg='valid_time   ')
      IF (i > 0) var%flags = IOR(IAND(var%flags,INPUT_ALLFLAGS-INPUT_VAR_VALID_MASK),(i-1)*INPUT_VAR_VALID_BEFORE)

      IF (LEN_TRIM(curr_opt%file_name) > 0) file_def => curr_opt ! New source file specified

      IF (LEN_TRIM(curr_opt%depend) > 0) CALL SetVarDependencies(var,curr_opt%depend)

      curr_opt => curr_opt%next
    ENDDO

    ! Reset "model" flag for the variable if settings have been applied (mo_input should take care of the variable update)
    IF (Cnt > 0) var%flags = IAND(var%flags,INPUT_ALLFLAGS - INPUT_VAR_MODEL)

    ! Associate new file(s) with variable if requested
    IF (ASSOCIATED(file_def)) THEN
      IF (TRIM(file_def%file_name) == '<model>') THEN ! No file, data will come from somewhere else in the model
        NULLIFY(var%parent_file)
      ELSE
        NULLIFY(prev_file)
        j  = 0
        st = 1
        DO WHILE (st <= LEN_TRIM(file_def%file_name))
          en = INDEX(file_def%file_name(st:),',') + st - 1
          IF (en < st) THEN
            en = LEN_TRIM(file_def%file_name) + 1
          ELSE
            IF (IAND(var%stat,INPUT_VAR_MULTI) == 0) CALL local_error('InputVarApplyExternalSettings', &
                                                                        'Multiple file names given for single-named variable')
            var%stat = IOR(var%stat,INPUT_VAR_MULTI_SPREAD) ! Variable is a superposition of vars from different files
          ENDIF          
          new_file => InputFileGetRef(file_def%file_name(st:en-1),'',0,sub_model=var%sub_model)
          IF (ASSOCIATED(new_file)) THEN
            IF (j==1) CALL local_error('InputVarApplyExternalSettings','File setting mismatch for multi-file variable')
            IF ((st>1) .AND. ASSOCIATED(prev_file)) THEN
               IF (.NOT. ASSOCIATED(new_file,prev_file%next)) CALL local_error('InputVarApplyExternalSettings', &
                                                    'File order mismatch for multi-file variable - previous file(s) were new')
            ENDIF
            j = 2
          ELSE ! Create new file, applying only basic settings, the rest will be set later
            IF (j==2) CALL local_error('InputVarApplyExternalSettings', &
                                       'File setting mismatch for multi-file variable - previous file(s) already defined')
            j = 1
            dt_file      = UNLIKELY_VAL
            dt_file_unit = ''
            IF (ASSOCIATED(var%parent_file)) THEN
              dt_file = var%parent_file%dt_file
              dt_file_unit = TRIM(var%parent_file%time_unit)
            ENDIF
            IF (file_def%dt_file /= UNLIKELY_VAL) dt_file = file_def%dt_file
            IF (LEN_TRIM(file_def%dt_file_unit) /= 0) dt_file_unit = TRIM(file_def%dt_file_unit)
            new_file => InputFileAdd(TRIM(file_def%file_name(st:en-1)), & 
            dt_file   = dt_file,                           &
            dt_unit   = TRIM(dt_file_unit),                &
            sub_model = var%sub_model)
          ENDIF
          IF (st>1) THEN
            IF (.NOT. ASSOCIATED(new_file,prev_file%next)) CALL local_error('InputVarApplyExternalSettings', &
                                                                            'File order mismatch for multi-file variable')
          ENDIF
          prev_file => new_file
          IF (st==1) THEN
            var%parent_file => new_file
          ELSE
            new_file%nvar = new_file%nvar + 1
          ENDIF
          st = en + 1
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE InputVarApplyExternalSettings

! Condenses dimension information from the chain in the variables edims
  LOGICAL FUNCTION InputVarGetDims(var)
  TYPE (input_var_list), POINTER :: var

  TYPE (input_dim_list), POINTER :: curr
  INTEGER :: i
  LOGICAL :: lGlobal

    lGlobal = .TRUE.
    curr => var%dims
    DO i=1,var%ndim
      IF (.NOT. ASSOCIATED(curr)) CALL local_error('InputVarAdd','#dimensions given for '//TRIM(var%name_var)// &
        ' is smaller than the #dimensions of its data pointer')
      var%edims(i) = curr%dim_data%size_local
      IF (curr%dim_data%size_local /= curr%dim_data%size_global) lGlobal = .FALSE.
      curr => curr%next
    ENDDO
      IF (ASSOCIATED(curr)) CALL local_error('InputVarAdd','#dimensions given for '//TRIM(var%name_var)// &
        ' is larger than the #dimensions of its data pointer')
    IF (lGlobal) var%flags = IOR(var%flags,INPUT_VAR_GLOBAL)

    InputVarGetDims = ALL(var%edims >= 0)

  END FUNCTION InputVarGetDims

! Finds external settings applying to given variable
  SUBROUTINE InputVarGetExternalSettings(var,opt,cnt)
  TYPE (input_var_list), POINTER :: var
  TYPE (input_opt_list), POINTER :: opt
  INTEGER,           INTENT(out) :: cnt

  TYPE (input_opt_list), POINTER :: curr, prev, last, prevopt
  LOGICAL :: lUse

    cnt = 0
    NULLIFY(opt,prev,last,prevopt)
    curr => AllInputOpts
    DO WHILE (ASSOCIATED(curr))
      lUse = (curr%sub_model==var%sub_model) .OR. (curr%sub_model<0) .OR. (var%sub_model<0)
      IF (lUse) THEN
        lUse = GetStringIndex(var%name_var,curr%var_name,Sep=',',lCaseS=.TRUE.) > 0
        IF (.NOT. lUse .AND. ASSOCIATED(var%group)) &
          lUse = GetStringIndex(var%group%name_group,curr%var_name,Sep=',',lCaseS=.TRUE.) > 0
      ENDIF
      IF (lUse) THEN
        cnt = cnt + 1
        ! Swap elements, so that all regarding this var are contiguous
        IF (ASSOCIATED(last)) THEN
          IF (IAND(curr%flags,INPUT_OPT_NAMELIST) == 0) THEN ! Put specification objects first, namelist specs last
            prev%next => curr%next
            curr%next => opt
            IF (ASSOCIATED(prevopt)) prevopt%next => curr
            opt       => curr
            curr      => prev%next
          ELSE
            IF (.NOT. ASSOCIATED(last%next,curr)) THEN ! Are they already ordered?
              prev%next => curr%next ! Remove curr from its current position
              curr%next => last%next ! Forward link of new position
              last%next => curr      ! Link to new position from last found (insert curr)
              last      => curr      ! Last found
              curr      => prev%next ! Next options to test
            ELSE
              prev => curr      ! Dazzling pointer to be able to take out 
              curr => curr%next ! Already in order -> just go to next element
            ENDIF
          ENDIF
        ELSE ! No previous match
          prevopt => prev      ! Pointer to last non-interesting opt ahead of the interesting ones
          opt     => curr      ! First match should be returned
          last    => curr      ! Last in list of matching
          prev    => curr      ! Dazzling pointer to be able to take out 
          curr    => curr%next ! Go to next element
        ENDIF
      ELSE ! curr is not of interest
        prev => curr      ! Dazzling pointer to be able to take out 
        curr => curr%next ! Go to next element
      ENDIF
    ENDDO

  END SUBROUTINE InputVarGetExternalSettings

! Recursively test yet unresolved dependencies
  RECURSIVE FUNCTION GetRemainingDependencies(depend) RESULT (remain)
  CHARACTER (len=256)            :: remain
  CHARACTER (len=*),  INTENT(in) :: depend

  TYPE (input_var_list), POINTER :: curr
  CHARACTER (len=256) :: curr_dep
  INTEGER Le, St, En, Fn, L2

    curr_dep = TRIM(depend)
    Le = LEN_TRIM(curr_dep)
    St = 1
    En = INDEX(curr_dep,',')
    IF (En < 1) En = Le + 1
    DO WHILE (St <= Le)
      curr => InputVarGetRef(curr_dep(St:En-1))
      IF (.NOT. ASSOCIATED(curr)) CALL local_error('InputVarAdd','Dependency involves a non-existing variable: '// &
                                                                 TRIM(curr_dep(St:En-1)))
      IF (ASSOCIATED(curr%depend)) THEN
        Fn = INDEX(TRIM(curr_dep),TRIM(curr%depend%name_var))
        IF (Fn >= 1) THEN
          L2 = LEN_TRIM(curr%depend%name_var)
          IF ((Fn+L2 > En) .OR. (curr_dep(Fn+L2:Fn+L2) == ',')) THEN ! This is one of the requested dependencies
            curr_dep = curr_dep(1:Fn-2)//curr_dep(Fn+L2:Le)
            curr_dep = GetRemainingDependencies(curr_dep)
            Le = LEN_TRIM(curr_dep)
            St = 1
          ELSE
            St = En + 1
          ENDIF
        ELSE
          St = En + 1
        ENDIF
      ELSE
        St = En + 1
      ENDIF
      En = INDEX(curr_dep(St:),',') + St - 1
      IF (En < St) En = Le + 1
    ENDDO
    remain = TRIM(curr_dep)

  END FUNCTION GetRemainingDependencies

! Set the dependencies of a variable
  SUBROUTINE SetVarDependencies(var,depend,depend_ref)
  TYPE (input_var_list),           POINTER    :: var
  CHARACTER (len=*),     OPTIONAL, INTENT(in) :: depend
  TYPE (input_var_list), OPTIONAL, POINTER    :: depend_ref

  TYPE (input_var_list), POINTER :: curr
  CHARACTER (len=256) :: remain
  INTEGER :: Le, St, En

    IF (PRESENT(depend)) THEN
      St = INDEX(depend,',')
      IF (St < 1) THEN ! Only one dependency - just assign it
        curr => InputVarGetRef(depend)
        IF (.NOT. ASSOCIATED(curr)) CALL local_error('InputVarAdd',TRIM(var%name_var)// &
                                                     ' is requested to depend on non-existing variable: '//TRIM(depend))
        var%depend => curr
      ELSE ! Multiple dependencies given
        ! Reset status of all dependency chain
        curr => AllInputVars
        DO WHILE (ASSOCIATED(curr))
          curr%stat = IAND(curr%stat,INPUT_ALLFLAGS - INPUT_VAR_DEPENDENT)
          curr => curr%next
        ENDDO
        remain = GetRemainingDependencies(depend) ! Resolve which are already in the dependency chain of each other
        ! Now distribute the remaining ones
        Le = LEN_TRIM(remain)
        St = 1
        En = INDEX(remain,',')
        IF (En < 1) En = Le + 1
        curr => var
        DO WHILE (St <= Le)
          curr%depend => InputVarGetRef(remain(St:En-1))
          DO WHILE (ASSOCIATED(curr%depend) .AND. (IAND(curr%stat,INPUT_VAR_DEPENDENT)==0))
            curr%stat = IOR(curr%stat,INPUT_VAR_DEPENDENT)
            curr => curr%depend
            DO WHILE(ASSOCIATED(curr%depend))
              curr => curr%depend
            ENDDO
            St = En + 1
            En = INDEX(remain(St:),',') + St - 1
            IF (En < St) En = Le + 1
          ENDDO
        ENDDO
      ENDIF
    ELSE
      IF (.NOT. PRESENT(depend_ref)) CALL local_error('InputVarAdd',TRIM(var%name_var)//': No method to resolve dependencies')
      ! This is easy - just assign the variable of dependency
      var%depend => depend_ref
    ENDIF

    ! Reset status of all dependency chain
    curr => AllInputVars
    DO WHILE (ASSOCIATED(curr))
      curr%stat = IAND(curr%stat,INPUT_ALLFLAGS - INPUT_VAR_DEPENDENT)
      curr => curr%next
    ENDDO

    ! Check for recursive dependencies - is this actually necessary? Since all dependence variable must be defined?!
    remain = ''
    curr => var
    DO WHILE (ASSOCIATED(curr))
      IF (IAND(curr%stat,INPUT_VAR_DEPENDENT) /= 0) &
        CALL local_error('InputVarAdd',TRIM(var%name_var)//' is part of a recursive dependency chain: '//TRIM(remain(3:)))
      curr%stat = IOR(curr%stat,INPUT_VAR_DEPENDENT)
      remain = TRIM(remain)//'=>'//TRIM(curr%name_var)
      curr => curr%depend
    ENDDO

  END SUBROUTINE SetVarDependencies

! Make sure that any variable will be processed after the variable(s) on which it depends
  SUBROUTINE SortDependencies
  TYPE (input_var_list), POINTER :: curr, dep, prev
  LOGICAL :: lCont

    NULLIFY(dep,prev)
    curr => AllInputVars
    DO WHILE (ASSOCIATED(curr))
      IF (ASSOCIATED(curr%depend)) THEN
        dep  => curr%next
        prev => curr
        lCont = .TRUE.
        DO WHILE (ASSOCIATED(dep))
          IF (ASSOCIATED(curr%depend,dep)) THEN ! Variable before its dependencies - we have to reorder
            prev%next => dep%next
            dep%next  => AllInputVars
            AllInputVars => dep
            NULLIFY(dep)
            curr => AllInputVars ! Start all over again
            lCont = .FALSE.
          ELSE
            prev => dep
            dep  => dep%next
          ENDIF
        ENDDO
        IF (lCont) curr => curr%next
      ELSE
        curr => curr%next
      ENDIF
    ENDDO

  END SUBROUTINE SortDependencies

! Removes the dependency on a specified variable from all other variables
  SUBROUTINE RemoveDependency(var)
  TYPE (input_var_list), POINTER :: var

  TYPE (input_var_list), POINTER :: curr

    curr => AllInputVars
    DO WHILE (ASSOCIATED(curr))
      IF (ASSOCIATED(curr%depend,var)) NULLIFY(curr%depend)
      curr => curr%next
    ENDDO

  END SUBROUTINE RemoveDependency

! Internal function to add (or find existing) variable independent of the number of dimensions. Returns true if already found
  FUNCTION InputVarAdd_gen(src_dims,ndim,dims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time,&
                           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val, &
                           miss_val,valid_range,weight,weights,lmodel,sub_model)
  TYPE (input_var_list), POINTER    :: InputVarAdd_gen  ! New variable
  CHARACTER (LEN=*),  INTENT(in)    :: src_dims         ! Comma separated list of names of variable dimensions (fastest first)
  INTEGER,            INTENT(in)    :: ndim             ! # dimensions
  INTEGER,            INTENT(inout) :: dims(:)          ! Size of dimensions
  LOGICAL,            INTENT(out)   :: lOld             ! Variable already in list?
  TYPE (input_file_list), OPTIONAL, POINTER    :: IFile ! File to add variable to
  CHARACTER (LEN=*),      OPTIONAL, INTENT(in) :: file_name  ! Name of file to add variable to
  CHARACTER (LEN=*),      OPTIONAL, INTENT(in) :: depend     ! Comma separated list of names of vars on which this var depends
  TYPE (input_var_list),  OPTIONAL, POINTER    :: depend_ref ! Reference to variable on which this variable depends
  INTEGER, OPTIONAL,  INTENT(in)    :: dt_update        ! Time between updates
  INTEGER, OPTIONAL,  INTENT(in)    :: miss             ! What to do when var is missing in file : See structure flag description 
  INTEGER, OPTIONAL,  INTENT(in)    :: mis_dim          ! What to do when dimensions are mismatching between file and variable 
  INTEGER, OPTIONAL,  INTENT(in)    :: data_action      ! How to combine file and model data 
  LOGICAL, OPTIONAL,  INTENT(in)    :: lauto            ! Switch off auto-reset of READ and UPDATE flags a every InputUpdate call
  INTEGER, OPTIONAL,  INTENT(in)    :: valid_time       ! When in interval are data valid
  INTEGER, OPTIONAL,  INTENT(in)    :: nIntrp           ! Number of time steps to store at once for (user) interpolation
  INTEGER, OPTIONAL,  INTENT(in)    :: nAtOnce          ! Number of time steps to read at once
  INTEGER, OPTIONAL,  INTENT(in)    :: nFuture          ! Number of future time steps needed for (user) interpolation
  INTEGER, OPTIONAL,  INTENT(in)    :: subset_index     ! Index of missing dimension to set for partial reading
  INTEGER, OPTIONAL,  INTENT(in)    :: sub_model        ! Sub model to which the variable should belong
  LOGICAL, OPTIONAL,  INTENT(in)    :: lmodel           ! Normally the model will take care of updating this variable
  REAL(dp),OPTIONAL,  INTENT(in)    :: mul,add          ! Scaling to correct for different units (d_mod = d_file*mul + add)
  REAL(dp),OPTIONAL,  INTENT(in)    :: fill_val,miss_val! Value to replace bad/out-of-range values with and bad value
  REAL(dp),OPTIONAL,  INTENT(in)    :: def_val,weight   ! Default value to use if var is missing and data weight for combinations
  REAL(dp),OPTIONAL,  INTENT(in)    :: valid_range(2)   ! Valid range for data
  REAL(dp),OPTIONAL,  INTENT(in)    :: weights(:)       ! Weighting of sub-dimensions on contraction
  REAL(dp),OPTIONAL,  INTENT(in)    :: def_vec(:)       ! Vector default value
  CHARACTER (len=*), OPTIONAL, INTENT(in) :: var_name   ! Name of variable
  CHARACTER (len=*), OPTIONAL, INTENT(in) :: alt_name   ! Alternative name of variable (may be overwritten externally)
  CHARACTER (len=*), OPTIONAL, INTENT(in) :: dt_unit    ! Unit of dt_update

  TYPE (input_var_list),  POINTER :: new_var, grp_var
  TYPE (input_opt_list),  POINTER :: opt
  TYPE (input_data_list), POINTER :: curr_data
  INTEGER :: dt, cnt, chk, grpset
  LOGICAL :: lGetDims, lOpen, lDtSet, lNormVar
  REAL(dp) :: mu, ad, fill(2)
  REAL(dp), POINTER :: tmp1(:), tmp7(:,:,:,:,:,:,:)

    IF (.NOT. moInputInitialized) CALL local_error('InputVarAdd','Please call InputInit before adding variables')
    IF (lInLoop) CALL local_error('InputVarAdd','Cannot add variables after entering the main loop')

    ! First check if variable already exist (we may reuse group-parameter variables)
    lNormVar = ndim >= -1
    IF (lNormVar) THEN
      new_var => InputVarGetRef(var_name,alt_name,sub_model)
      IF (ASSOCIATED(new_var)) THEN
        lOld = .TRUE.
        InputVarAdd_gen => new_var
        RETURN
      ENDIF
    ELSE
      new_var => InputVarGetRef(var_name,sub_model=sub_model,list=AllInputGrpVars)
    ENDIF
    lOld = ASSOCIATED(new_var)

    ! Create new variable and add it to list of variables
    IF (.NOT. ASSOCIATED(new_var)) THEN
      CALL InputVarNew(new_var)
      new_var%next => AllInputVars ! This swaps the order of the elements, but since we in all cases should
      AllInputVars => new_var      ! process them all, this doesn't really matter.
    ENDIF
    lOpen = .FALSE.
    lDtSet = .FALSE.

    ! Process parameters specified for the group
    grpset = 0
    IF (ASSOCIATED(CurrGroup)) THEN
      grp_var => InputVarGetRef(TRIM(CurrGroup%name_group),list=AllInputGrpVars)
      IF (ASSOCIATED(grp_var)) THEN
        grpset = grp_var%stat
        IF (IAND(grpset,INPUT_GRPVAR_SET_IFILE)      /= 0) new_var%parent_file  => grp_var%parent_file
        IF (IAND(grpset,INPUT_GRPVAR_SET_DT)         /= 0) new_var%dt_update    =  grp_var%dt_update
        IF (IAND(grpset,INPUT_GRPVAR_SET_DEPEND)     /= 0) new_var%depend       => grp_var%depend
        IF (IAND(grpset,INPUT_GRPVAR_SET_NINTRP)     /= 0) new_var%nIntrp       =  grp_var%nIntrp
        IF (IAND(grpset,INPUT_GRPVAR_SET_NATONCE)    /= 0) new_var%nAtOnce      =  grp_var%nAtOnce
        IF (IAND(grpset,INPUT_GRPVAR_SET_NFUTURE)    /= 0) new_var%nFuture      =  grp_var%nFuture
        IF (IAND(grpset,INPUT_GRPVAR_SET_INDEX)      /= 0) new_var%subset_index =  grp_var%subset_index
        IF (IAND(grpset,INPUT_GRPVAR_SET_SUB_MODEL)  /= 0) new_var%sub_model    =  grp_var%sub_model
        IF (IAND(grpset,INPUT_GRPVAR_SET_MISS_ACT)   /= 0) &
          new_var%flags=IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_MISS_MASK)    ,IAND(grp_var%flags,INPUT_VAR_MISS_MASK))
        IF (IAND(grpset,INPUT_GRPVAR_SET_MIS_DIM)    /= 0) &
          new_var%flags=IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_DIM_MIS_MASK) ,IAND(grp_var%flags,INPUT_VAR_DIM_MIS_MASK))
        IF (IAND(grpset,INPUT_GRPVAR_SET_DATA_ACT)   /= 0) &
          new_var%flags=IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_ACTION_MASK)  ,IAND(grp_var%flags,INPUT_VAR_ACTION_MASK))
        IF (IAND(grpset,INPUT_GRPVAR_SET_VALID_TIME) /= 0) &
          new_var%flags=IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_VALID_MASK)   ,IAND(grp_var%flags,INPUT_VAR_VALID_MASK))
        IF (IAND(grpset,INPUT_GRPVAR_SET_LAUTO)      /= 0) &
          new_var%flags=IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_AUTORESET_OFF),IAND(grp_var%flags,INPUT_VAR_AUTORESET_OFF))
        IF (IAND(grpset,INPUT_GRPVAR_SET_MISS_VAL)   /= 0) THEN
          tmp7 => InputDataGetSpecified(grp_var%aux_data,INPUT_DATA_FILL,INPUT_DATA_NO_STEP)
          tmp1 => tmp7(:,1,1,1,1,1,1)
          CALL ReplaceAuxData1D(new_var,INPUT_DATA_FILL,INPUT_DATA_NO_STEP,tmp1)
        ENDIF
        IF (IAND(grpset,INPUT_GRPVAR_SET_WEIGHT)     /= 0) THEN
          tmp7 => InputDataGetSpecified(grp_var%aux_data,INPUT_DATA_NUDGE_SCA,0)
          tmp1 => tmp7(:,1,1,1,1,1,1)
          CALL ReplaceAuxData1D(new_var,INPUT_DATA_NUDGE_SCA,0,tmp1)
        ENDIF
        IF (IAND(grpset,INPUT_GRPVAR_SET_WEIGHTS)    /= 0) THEN
          tmp7 => InputDataGetSpecified(grp_var%aux_data,INPUT_DATA_WEIGHTS,0)
          tmp1 => tmp7(:,1,1,1,1,1,1)
          CALL ReplaceAuxData1D(new_var,INPUT_DATA_WEIGHTS,0,tmp1)
        ENDIF
      ENDIF
    ENDIF

    ! Set variable properties which are independent of the data source
    IF (PRESENT(alt_name)) new_var%name_alt = alt_name
    IF (PRESENT(var_name)) new_var%name_var = var_name
    IF (lNormVar) THEN
      IF (ndim >= 0) THEN
        new_var%ndim = ndim
      ELSE
        new_var%ndim = CountSubStrings(src_dims)
      ENDIF
    ENDIF
    new_var%sub_model = CurrSubModel
    IF (PRESENT(sub_model)) new_var%sub_model = sub_model
    IF (PRESENT(miss)) THEN
      new_var%flags = IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_MISS_MASK),IAND(miss,INPUT_VAR_MISS_MASK))
    ELSE
      IF (IAND(grpset,INPUT_GRPVAR_SET_MISS_ACT) == 0) THEN
        IF (PRESENT(def_val) .OR. PRESENT(def_vec)) THEN
          new_var%flags = IOR(new_var%flags,INPUT_VAR_MISS_WARN)
        ELSE
          new_var%flags = IOR(new_var%flags,INPUT_VAR_MISS_ERR)
        ENDIF
      ENDIF
    ENDIF
    IF (PRESENT(lauto)) THEN
      IF (lauto) THEN
        new_var%flags = IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_AUTORESET_OFF),INPUT_VAR_AUTORESET_OFF)
      ELSE
        new_var%flags = IAND(new_var%flags,INPUT_ALLFLAGS - INPUT_VAR_AUTORESET_OFF)
      ENDIF
    ENDIF
    IF (PRESENT(mis_dim)) THEN
      new_var%flags = IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_DIM_MIS_MASK),IAND(mis_dim,INPUT_VAR_DIM_MIS_MASK))
    ELSE
      IF (IAND(grpset,INPUT_GRPVAR_SET_MIS_DIM) == 0) new_var%flags = IAND(new_var%flags,INPUT_ALLFLAGS - INPUT_VAR_DIM_MIS_MASK)
    ENDIF
    IF (PRESENT(data_action)) THEN
      new_var%flags = IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_ACTION_MASK),IAND(data_action,INPUT_VAR_ACTION_MASK))
    ELSE
      IF (IAND(grpset,INPUT_GRPVAR_SET_DATA_ACT) == 0) new_var%flags = IAND(new_var%flags,INPUT_ALLFLAGS - INPUT_VAR_ACTION_MASK)
    ENDIF
    IF (PRESENT(subset_index)) new_var%subset_index = subset_index
    new_var%flags = IAND(new_var%flags,INPUT_ALLFLAGS - INPUT_VAR_MODEL)
    IF (PRESENT(lmodel)) THEN
      IF (lmodel) new_var%flags = IOR(new_var%flags,INPUT_VAR_MODEL)
    ENDIF
    fill(:) = UNLIKELY_RVAL
    mu      = 1._dp
    ad      = 0._dp
    IF (PRESENT(mul))    mu = mul ! E.g. kg -> g
    IF (PRESENT(add))    ad = add ! E.g. degC -> K
    IF ((mu/=1._dp) .OR. (ad/=0._dp)) &
                              CALL ReplaceAuxData1D(new_var,INPUT_DATA_RESCALE,INPUT_DATA_NO_STEP,(/mu,ad/))
    IF (PRESENT(weight))      CALL ReplaceAuxData1D(new_var,INPUT_DATA_NUDGE_SCA,               0,(/weight /))
    IF (PRESENT(weights))     CALL ReplaceAuxData1D(new_var,INPUT_DATA_WEIGHTS,                 0,  weights  )
    IF (PRESENT(def_val))     CALL ReplaceAuxData1D(new_var,INPUT_DATA_DEFAULT,                 0,(/def_val/))
    IF (PRESENT(def_vec))     CALL ReplaceAuxData1D(new_var,INPUT_DATA_DEFAULT,                 0,  def_vec  )
    IF (PRESENT(valid_range)) CALL ReplaceAuxData1D(new_var,INPUT_DATA_RANGE,  INPUT_DATA_NO_STEP,(/UNLIKELY_RVAL,               &
                                                                                                  valid_range(1),valid_range(2)/))
    IF (PRESENT(fill_val)) fill(1) = fill_val
    IF (PRESENT(miss_val)) fill(2) = miss_val
    IF (PRESENT(miss_val) .OR. PRESENT(fill_val)) &
      CALL ReplaceAuxData1D(new_var,INPUT_DATA_FILL,INPUT_DATA_NO_STEP,fill)

    IF (PRESENT(dt_update) .OR. PRESENT(dt_unit)) THEN
      dt = InterpretDtUnit(dt_update,dt_unit,INPUT_VAR_FILE_TIME_STEP)
      IF ((dt < model_dt) .AND. (dt > 0)) THEN
        CALL local_message('InputVarAdd','Warning: Too small update time step set to model time step.')
        dt = model_dt
      ENDIF
      IF (dt > 0) THEN ! dt given in seconds -> possible to check and set fixed update interval
        IF (dt / model_dt /= INT(REAL(dt,dp)/REAL(model_dt,dp))) THEN
          CALL local_error('InputVarAdd','Not an integer number of model time steps between updates');
        ENDIF
        new_var%nupdate = dt / model_dt
      ENDIF
      new_var%dt_update = dt
      IF (dt /= INPUT_VAR_FILE_TIME_STEP) THEN
        new_var%stat = IOR(new_var%stat,INPUT_VAR_DT_SET)
        lDtSet = .TRUE.
      ENDIF
    ENDIF
    IF (PRESENT(valid_time)) THEN
      new_var%flags = IOR(IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_VALID_MASK),IAND(valid_time,INPUT_VAR_VALID_MASK))
    ELSE ! Default is INPUT_VAR_VALID_ACTUAL
      IF (IAND(grpset,INPUT_GRPVAR_SET_VALID_TIME) == 0) new_var%flags = IAND(new_var%flags,INPUT_ALLFLAGS-INPUT_VAR_VALID_MASK)
    ENDIF

    ! Time settings for correct interpolation fields
    IF (IAND(grpset,INPUT_GRPVAR_SET_NATONCE) == 0) new_var%nAtOnce   = 1
    IF (PRESENT(nIntrp)) THEN
      new_var%nIntrp    =  nIntrp
      new_var%flags     =  IOR(new_var%flags,INPUT_VAR_INTERPOL_USER)
      IF (IAND(grpset,INPUT_GRPVAR_SET_NFUTURE) == 0) new_var%nFuture = NINT(REAL(nIntrp,dp)/2._dp+0.01_dp)
    ELSE
      IF (IAND(grpset,INPUT_GRPVAR_SET_NFUTURE) == 0) new_var%nFuture = 1
    ENDIF
    IF (PRESENT(nAtOnce)) THEN
      IF (.NOT. PRESENT(nIntrp)) CALL local_error('InputVarAdd','Variable: '//TRIM(new_var%name_var)// &
                                                  'Specifying nAtOnce requires specified nIntrp')
      new_var%nAtOnce = nAtOnce
    ENDIF
    IF (PRESENT(nFuture)) THEN
      IF (.NOT. PRESENT(nIntrp)) CALL local_error('InputVarAdd','Variable: '//TRIM(new_var%name_var)// & 
                                                  'Specifying nFuture requires specified nIntrp')
      new_var%nFuture = nFuture
    ENDIF

    ! Associate dimensions
    IF (lNormVar) THEN
      new_var%dims => InputDimListFromString(src_dims,INPUT_DIM_VAR,'InputVarAdd: '//TRIM(new_var%name_var), &
                                             new_var%sub_model,add='_iloc')

      ! Get sizes of dimensions (dim_lo, dim_hi, edims)
      new_var%edims(:) = dims(:)               ! Get any explicitly specified dimension sizes
      new_var%edims(new_var%ndim+1:7) = 1      ! For variables without model data and thus only indirectly known #dims
      lGetDims = ANY(dims(1:new_var%ndim) < 0) ! Check - do we know every size?
      IF (lGetDims) lGetDims = .NOT. InputVarGetDims(new_var) ! No - try to get it from the list of dimensions

      ! Associate group to get all settings (accept for now, consistency checks are done below)
      IF (ASSOCIATED(CurrGroup)) new_var%group => CurrGroup
    ELSE
      lGetDims = .FALSE.
    ENDIF 

    IF (PRESENT(depend) .OR. PRESENT(depend_ref)) CALL SetVarDependencies(new_var,depend,depend_ref)

    ! Include externally specified settings
    IF (PRESENT(IFile)) new_var%parent_file => IFile ! Provide link to file to obtain "defaults". File may be overwritten
    IF (PRESENT(file_name)) new_var%parent_file => InputFileAdd(file_name)
    CALL InputVarGetExternalSettings(new_var,opt,cnt)
    CALL InputVarApplyExternalSettings(new_var,opt,cnt)
    IF ((IAND(new_var%stat,INPUT_VAR_DT_SET) /= 0) .OR. (IAND(grpset,INPUT_GRPVAR_SET_DT) /= 0)) &
      lDtSet = new_var%dt_update /= INPUT_VAR_FILE_TIME_STEP
    IF (IAND(new_var%flags,INPUT_VAR_MISS_MASK) == INPUT_VAR_MISS_DEFAULT) &
      new_var%flags = IOR(new_var%flags,INPUT_VAR_MISS_ERR)

    IF ((IAND(new_var%flags,INPUT_VAR_DIM_MIS_SUBSET) /= 0) .AND. (new_var%subset_index < 1)) &
      CALL local_error('InputVarAdd','Variable: '//TRIM(new_var%name_var)//                   &
                       ': Set subset option selected without specifying a valid subset')

    IF (ASSOCIATED(new_var%parent_file)) THEN
      ! Update file status and force checks of variables and variable dimensions
      IF (lNormVar) new_var%parent_file%nVar = new_var%parent_file%nVar + 1
      new_var%parent_file%stat = IAND(new_var%parent_file%stat,INPUT_ALLFLAGS-INPUT_FILE_CHECK_VAR_DIMS-INPUT_FILE_CHECK_VAR)

      ! Crosschecking with file parameters
      IF (.NOT. ASSOCIATED(new_var%parent_file,DummyFile)) THEN
        IF (.NOT. lDtSet .AND. (new_var%dt_update /= INPUT_VAR_FILE_TIME_STEP)) THEN
          IF (IAND(new_var%parent_file%flags,INPUT_FILE_INITIAL) == 0) THEN
            WRITE(message_text,*) 'Warning: No requested update time interval given for ',TRIM(new_var%name_var), &
                                  '. Assuming model time step: ',model_dt
            IF (IAND(mo_debug,INPUT_DBG_MASK) /= 0) CALL local_message('InputVarAdd',message_text)
            new_var%dt_update = model_dt
          ELSE
            WRITE(message_text,*) 'Warning: No requested update time interval given for ',TRIM(new_var%name_var), &
                                  '. Its associated file is an initial file, and variable is thus assumed to be initial.'
            IF (IAND(mo_debug,INPUT_DBG_MASK) /= 0) CALL local_message('InputVarAdd',message_text)
            new_var%dt_update = 0
          ENDIF
          new_var%stat = IOR(new_var%stat,INPUT_VAR_DT_SET)
          lDtSet = .TRUE. ! In the case the variable is associated with an initial file, new%dt_update must stay 0
        ELSE
          IF ((new_var%dt_update /= 0) .AND. (new_var%dt_update /= INPUT_VAR_FILE_TIME_STEP) .AND. &
              (IAND(new_var%parent_file%flags,INPUT_FILE_INITIAL) /= 0))         CALL local_error( &
            'InputVarAdd','Attept to assign non-initial variable '//TRIM(new_var%name_var)//' to an initial file.')
        ENDIF
      ENDIF

      IF (ASSOCIATED(new_var%parent_file,DummyFile) .AND. (IAND(new_var%flags,INPUT_VAR_VALID_MASK)/=INPUT_VAR_VALID_ACTUAL)) THEN
        CALL local_message('InputVarAdd','Variable: '//TRIM(new_var%name_var)// &
                                         ' is not associated with a file and must be treated as being valid now!')
        new_var%flags = IAND(new_var%flags,INPUT_ALLFLAGS - INPUT_VAR_VALID_MASK)
      ENDIF

      ! Get dimensions sizes if possible if they are not already available
      chk = 0
      IF (lGetDims)     chk =       INPUT_FILE_CHECK_DIM_SIZE
      IF (.NOT. lDtSet) chk = chk + INPUT_FILE_CHECK_TIME_STEP + INPUT_FILE_CHECK_VAR + INPUT_FILE_CHECK_DIM_SIZE
      IF (lGetDims .OR. .NOT. lDtSet) THEN
        IF (ASSOCIATED(new_var%parent_file,DummyFile)) CALL local_error('InputVarAdd',TRIM(var_name)// &
           ': No file associated from which dimension sizes and/or time parameters can be obtained.')
        IF (IAND(new_var%stat,INPUT_VAR_MULTI) /= 0) CALL local_error('InputVarAdd','Variable: '//TRIM(new_var%name_var)// &
          'Dimension sizes must be determined before adding multi-file variables')
        CALL InputFileOpen(new_var%parent_file,new_var%parent_file%time_start(1), &
                                               new_var%parent_file%time_start(2), &
                                               new_var%parent_file%time_start(3), &
                                               check=chk)
        IF (lNormVar) lGetDims = .NOT. InputVarGetDims(new_var) ! Try again to get it from the list of dimensions
        IF (.NOT. lDtSet) THEN
          IF (ASSOCIATED(new_var%parent_file)) THEN
            IF (new_var%parent_file%dt_file /= UNLIKELY_VAL) new_var%dt_update = new_var%parent_file%dt_file
          ELSE
            new_var%dt_update = 0 ! Default provided and file not found
            new_var%stat = IOR(new_var%stat,INPUT_VAR_INITIALIZED)
          ENDIF
          lDtSet = .TRUE.
          new_var%stat = IOR(new_var%stat,INPUT_VAR_DT_SET)
        ENDIF
        lOpen = .TRUE.
      ENDIF
    ENDIF

    IF (.NOT. lDtSet) THEN ! Not assigned to a file and time step not specified
      new_var%dt_update = model_dt
      WRITE(message_text,*) 'Warning: No requested update time interval given for ',TRIM(new_var%name_var), &
                            '. Assuming model time step: ',model_dt
      IF ((IAND(mo_debug,INPUT_DBG_MASK) /= 0) .AND. (IAND(new_var%flags,INPUT_VAR_MODEL)==0)) &
        CALL local_message('InputVarAdd',message_text)
    ENDIF

    IF (lGetDims) CALL local_error('InputVarAdd','Could not obtain size of one or more dimensions of variable '// &
                                                  TRIM(new_var%name_var))
    dims(1:new_var%ndim) = new_var%edims(1:new_var%ndim)
    IF (new_var%ndim >= 0) THEN
      ALLOCATE(new_var%dta(new_var%edims(1),new_var%edims(2),new_var%edims(3),new_var%edims(4), &
                           new_var%edims(5),new_var%edims(6),new_var%edims(7)))

      ! Fill data with default value if provided
      tmp7 => InputDataGetSpecified(new_var%aux_data,INPUT_DATA_DEFAULT,0)
      IF (ASSOCIATED(tmp7)) THEN
        IF (SIZE(tmp7) == 1) new_var%dta(:,:,:,:,:,:,:) = tmp7(1,1,1,1,1,1,1)
      ENDIF

      ! Associate last buffer with real variable
      IF (lOpen) THEN
        curr_data => new_var%buffers
        IF (ASSOCIATED(curr_data)) THEN
          DO WHILE (ASSOCIATED(curr_data%next))
            curr_data => curr_data%next
          ENDDO
          curr_data%dta => new_var%dta
        ENDIF
      ENDIF
    ENDIF

    ! Determine if variable should be added to the current group
    IF (ASSOCIATED(CurrGroup)) THEN
      IF (CurrGroup%dt_group == UNLIKELY_VAL) THEN ! No variable in group - add this one
        CurrGroup%dt_group   =  new_var%dt_update
        CurrGroup%nFuture    =  new_var%nFuture
        CurrGroup%nIntrp     =  new_var%nIntrp
        CurrGroup%nAtOnce    =  new_var%nAtOnce
        CurrGroup%flags      =  IAND(new_var%flags,INPUT_VAR_VALID_MASK)
        CurrGroup%group_file => new_var%parent_file
        new_var%group        => CurrGroup
      ELSE
        IF ((CurrGroup%dt_group == new_var%dt_update) .AND. (CurrGroup%nFuture == new_var%nFuture) .AND. &
            (CurrGroup%nIntrp   == new_var%nIntrp)    .AND. (CurrGroup%nAtOnce == new_var%nAtOnce) .AND. &
            (IAND(CurrGroup%flags,INPUT_VAR_VALID_MASK) ==  IAND(new_var%flags,INPUT_VAR_VALID_MASK))) THEN 
          new_var%group => CurrGroup
          IF (.NOT. ASSOCIATED(new_var%parent_file,CurrGroup%group_file) .AND. ASSOCIATED(CurrGroup%group_file)) THEN
            IF (ASSOCIATED(new_var%parent_file)) THEN
              IF (IAND(new_var%parent_file%flags,INPUT_FILE_INITIAL) /= IAND(CurrGroup%group_file%flags,INPUT_FILE_INITIAL)) THEN
                IF (IAND(mo_debug,INPUT_DBG_MASK) /= 0) &
                  CALL local_message('WARNING',TRIM(new_var%name_var)//' could not be added to group '// &
                                     TRIM(CurrGroup%name_group)//' due to mismatching file time parameters.')
                NULLIFY(new_var%group)
              ENDIF
            ENDIF
            CurrGroup%flags = IOR(CurrGroup%flags,INPUT_GROUP_MULTI_FILE)
          ENDIF
        ELSE
          IF (IAND(mo_debug,INPUT_DBG_MASK) /= 0) &
            CALL local_message('WARNING',TRIM(new_var%name_var)//' could not be added to group '// &
                               TRIM(CurrGroup%name_group)//' due to mismatching variable time parameters.')
        ENDIF
      ENDIF
    ENDIF

    InputVarAdd_gen => new_var

  END FUNCTION InputVarAdd_gen

! Adds a new variable without model data to the variable list
  SUBROUTINE InputVarNDAdd(src_dims,var_name,IFile,dt_update,dt_unit,file_name,miss,mis_dim,data_action,dims,                    &
                           valid_time,lauto,var_ref,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,weight,weights, &
                           def_vec,fill_val,miss_val,valid_range,depend,depend_ref,lmodel,sub_model)
  CHARACTER (len=*),                INTENT(IN)    :: src_dims
  INTEGER,                OPTIONAL, INTENT(in)    :: dt_update, miss, mis_dim, data_action, nIntrp, nAtOnce, nFuture
  INTEGER,                OPTIONAL, INTENT(in)    :: sub_model, valid_time, subset_index
  LOGICAL,                OPTIONAL, INTENT(in)    :: lauto, lmodel
  INTEGER,                OPTIONAL, INTENT(inout) :: dims(:)
  REAL(dp),               OPTIONAL, INTENT(in)    :: mul,add,def_val,fill_val,miss_val,valid_range(2),weight,weights(:),def_vec(:)
  TYPE (input_var_list),  OPTIONAL, POINTER       :: var_ref, depend_ref
  CHARACTER (len=*),      OPTIONAL, INTENT(in)    :: var_name, alt_name, dt_unit, file_name, depend
  TYPE (input_file_list), OPTIONAL, POINTER       :: IFile

  TYPE (input_var_list), POINTER :: new
  INTEGER :: VDims(7)
  LOGICAL :: lOld

    VDims(:) = -1
    IF (PRESENT(dims)) VDims(1:SIZE(dims)) = dims(:)

    new => InputVarAdd_gen(src_dims,-1,VDims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time,&
           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val,                &
           miss_val,valid_range,weight,weights,lmodel,sub_model)
    IF (PRESENT(var_ref)) var_ref => new

  END SUBROUTINE InputVarNDAdd

! Adds a new scalar variable to the variable list
  SUBROUTINE InputVar0dAdd(ptr,var_name,IFile,dt_update,dt_unit,file_name,miss,mis_dim,data_action,lcopy,ldealloc,               &
                           valid_time,lauto,var_ref,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,weight,weights, &
                           def_vec,fill_val,miss_val,valid_range,depend,depend_ref,lmodel,sub_model,p4)
  REAL (dp),                        POINTER       :: ptr
  INTEGER,                OPTIONAL, INTENT(in)    :: dt_update, miss, mis_dim, data_action, nIntrp, nAtOnce, nFuture
  INTEGER,                OPTIONAL, INTENT(in)    :: sub_model, valid_time, subset_index
  LOGICAL,                OPTIONAL, INTENT(in)    :: lcopy, ldealloc, lmodel, lauto
  REAL(dp),               OPTIONAL, INTENT(in)    :: mul,add,def_val,fill_val,miss_val,valid_range(2),weight,weights(:),def_vec(:)
  TYPE (input_var_list),  OPTIONAL, POINTER       :: var_ref, depend_ref
  CHARACTER (len=*),      OPTIONAL, INTENT(in)    :: var_name, alt_name, dt_unit, file_name, depend
  TYPE (input_file_list), OPTIONAL, POINTER       :: IFile
  REAL(dp),               OPTIONAL, POINTER       :: P4(:,:,:,:)

  TYPE (input_var_list), POINTER :: new
  INTEGER :: VDims(7)
  LOGICAL :: lOld

    VDims = (/1,1,1,1,1,1,1/)

    new => InputVarAdd_gen('',0,VDims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time, &
           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val,          &
           miss_val,valid_range,weight,weights,lmodel,sub_model)
    IF (.NOT. lOld .AND. ASSOCIATED(ptr)) THEN
      IF (PRESENT(lcopy)) THEN
        IF (lcopy) new%dta(1,1,1,1,1,1,1) = ptr
      ENDIF
      IF (PRESENT(ldealloc)) THEN 
        IF (ldealloc) DEALLOCATE(ptr)
      ENDIF
    ENDIF
    ptr => new%dta(1,1,1,1,1,1,1)

    IF (PRESENT(var_ref)) var_ref => new
    IF (PRESENT(P4)) P4 => new%dta(:,:,:,:,1,1,1)

  END SUBROUTINE InputVar0dAdd

! Adds a new 1d variable to the variable list
  SUBROUTINE InputVar1dAdd(src_dims,ptr,var_name,IFile,dt_update,dt_unit,file_name,miss,mis_dim,data_action,dims,lcopy,ldealloc, &
                           valid_time,lauto,var_ref,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,weight,weights, &
                           def_vec,fill_val,miss_val,valid_range,depend,depend_ref,lmodel,sub_model,p4)
  REAL (dp),                        POINTER       :: ptr(:)
  CHARACTER (len=*),                INTENT(IN)    :: src_dims
  INTEGER,                OPTIONAL, INTENT(in)    :: dt_update, miss, mis_dim, data_action, nIntrp, nAtOnce, nFuture
  INTEGER,                OPTIONAL, INTENT(in)    :: sub_model, valid_time, subset_index
  LOGICAL,                OPTIONAL, INTENT(in)    :: lcopy, ldealloc, lmodel, lauto
  INTEGER,                OPTIONAL, INTENT(inout) :: dims(:)
  REAL(dp),               OPTIONAL, INTENT(in)    :: mul,add,def_val,fill_val,miss_val,valid_range(2),weight,weights(:),def_vec(:)
  TYPE (input_var_list),  OPTIONAL, POINTER       :: var_ref, depend_ref
  CHARACTER (len=*),      OPTIONAL, INTENT(in)    :: var_name, alt_name, dt_unit, file_name, depend
  TYPE (input_file_list), OPTIONAL, POINTER       :: IFile
  REAL(dp),               OPTIONAL, POINTER       :: P4(:,:,:,:)

  TYPE (input_var_list), POINTER :: new
  INTEGER :: VDims(7)
  LOGICAL :: lOld

    VDims = (/-1,1,1,1,1,1,1/)
    IF (PRESENT(dims))   VDims(1) = dims(1)
    IF (ASSOCIATED(ptr)) VDims(1) = SIZE(ptr,1)

    new => InputVarAdd_gen(src_dims,1,VDims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time, &
           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val,                &
           miss_val,valid_range,weight,weights,lmodel,sub_model)
    IF (.NOT. lOld .AND. ASSOCIATED(ptr)) THEN
      IF (PRESENT(lcopy)) THEN
        IF (lcopy) new%dta(:,1,1,1,1,1,1) = ptr(:)
      ENDIF
      IF (PRESENT(ldealloc)) THEN 
        IF (ldealloc) DEALLOCATE(ptr)
      ENDIF
    ENDIF
    ptr => new%dta(:,1,1,1,1,1,1)

    IF (PRESENT(var_ref)) var_ref => new
    IF (PRESENT(P4)) P4 => new%dta(:,:,:,:,1,1,1)

  END SUBROUTINE InputVar1dAdd

! Adds a new 2d variable to the variable list
  SUBROUTINE InputVar2dAdd(src_dims,ptr,var_name,IFile,dt_update,dt_unit,file_name,miss,mis_dim,data_action,dims,lcopy,ldealloc, &
                           valid_time,lauto,var_ref,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,weight,weights, &
                           def_vec,fill_val,miss_val,valid_range,depend,depend_ref,lmodel,sub_model,P4)
  REAL (dp),                        POINTER       :: ptr(:,:)
  CHARACTER (len=*),                INTENT(IN)    :: src_dims
  INTEGER,                OPTIONAL, INTENT(in)    :: dt_update, miss, mis_dim, data_action, nIntrp, nAtOnce, nFuture
  INTEGER,                OPTIONAL, INTENT(in)    :: sub_model, valid_time, subset_index
  LOGICAL,                OPTIONAL, INTENT(in)    :: lcopy, ldealloc, lmodel, lauto
  INTEGER,                OPTIONAL, INTENT(inout) :: dims(:)
  REAL(dp),               OPTIONAL, INTENT(in)    :: mul,add,def_val,fill_val,miss_val,valid_range(2),weight,weights(:),def_vec(:)
  TYPE (input_var_list),  OPTIONAL, POINTER       :: var_ref, depend_ref
  CHARACTER (len=*),      OPTIONAL, INTENT(in)    :: var_name, alt_name, dt_unit, file_name, depend
  TYPE (input_file_list), OPTIONAL, POINTER       :: IFile
  REAL(dp),               OPTIONAL, POINTER       :: P4(:,:,:,:)

  TYPE (input_var_list), POINTER :: new
  INTEGER :: VDims(7)
  LOGICAL :: lOld

    VDims = (/-1,-1,1,1,1,1,1/)
    IF (PRESENT(dims))   VDims(1:2) = dims(1:2)
    IF (ASSOCIATED(ptr)) VDims(1:2) = SHAPE(ptr)

    new => InputVarAdd_gen(src_dims,2,VDims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time, &
           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val,                &
           miss_val,valid_range,weight,weights,lmodel,sub_model)
    IF (.NOT. lOld .AND. ASSOCIATED(ptr)) THEN
      IF (PRESENT(lcopy)) THEN
        IF (lcopy) new%dta(:,:,1,1,1,1,1) = ptr(:,:)
      ENDIF
      IF (PRESENT(ldealloc)) THEN
        IF (ldealloc) DEALLOCATE(ptr)
      ENDIF
    ENDIF
    ptr => new%dta(:,:,1,1,1,1,1)

    IF (PRESENT(var_ref)) var_ref => new
    IF (PRESENT(P4)) P4 => new%dta(:,:,:,:,1,1,1)

  END SUBROUTINE InputVar2dAdd

! Adds a new 3d variable to the variable list
  SUBROUTINE InputVar3dAdd(src_dims,ptr,var_name,IFile,dt_update,dt_unit,file_name,miss,mis_dim,data_action,dims,lcopy,ldealloc, &
                           valid_time,lauto,var_ref,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,weight,weights, &
                           def_vec,fill_val,miss_val,valid_range,depend,depend_ref,lmodel,sub_model,P4)
  REAL (dp),                        POINTER       :: ptr(:,:,:)
  CHARACTER (len=*),                INTENT(IN)    :: src_dims
  INTEGER,                OPTIONAL, INTENT(in)    :: dt_update, miss, mis_dim, data_action, nIntrp, nAtOnce, nFuture
  INTEGER,                OPTIONAL, INTENT(in)    :: sub_model, valid_time, subset_index
  LOGICAL,                OPTIONAL, INTENT(in)    :: lcopy, ldealloc, lmodel, lauto
  INTEGER,                OPTIONAL, INTENT(inout) :: dims(:)
  REAL(dp),               OPTIONAL, INTENT(in)    :: mul,add,def_val,fill_val,miss_val,valid_range(2),weight,weights(:),def_vec(:)
  TYPE (input_var_list),  OPTIONAL, POINTER       :: var_ref, depend_ref
  CHARACTER (len=*),      OPTIONAL, INTENT(in)    :: var_name, alt_name, dt_unit, file_name, depend
  TYPE (input_file_list), OPTIONAL, POINTER       :: IFile
  REAL(dp),               OPTIONAL, POINTER       :: P4(:,:,:,:)

  TYPE (input_var_list), POINTER :: new
  INTEGER :: VDims(7)
  LOGICAL :: lOld

    VDims = (/-1,-1,-1,1,1,1,1/)
    IF (PRESENT(dims))   VDims(1:3) = dims(1:3)
    IF (ASSOCIATED(ptr)) VDims(1:3) = SHAPE(ptr)

    new => InputVarAdd_gen(src_dims,3,VDims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time, &
           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val,                &
           miss_val,valid_range,weight,weights,lmodel,sub_model)
    IF (.NOT. lOld .AND. ASSOCIATED(ptr)) THEN
      IF (PRESENT(lcopy)) THEN 
        IF (lcopy) new%dta(:,:,:,1,1,1,1) = ptr(:,:,:)
      ENDIF
      IF (PRESENT(ldealloc)) THEN
        IF (ldealloc) DEALLOCATE(ptr)
      ENDIF
    ENDIF
    ptr => new%dta(:,:,:,1,1,1,1)

    IF (PRESENT(var_ref)) var_ref => new
    IF (PRESENT(P4)) P4 => new%dta(:,:,:,:,1,1,1)

  END SUBROUTINE InputVar3dAdd

! Adds a new 4d variable to the variable list
  SUBROUTINE InputVar4dAdd(src_dims,ptr,var_name,IFile,dt_update,dt_unit,file_name,miss,mis_dim,data_action,dims,lcopy,ldealloc, &
                           valid_time,lauto,var_ref,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,weight,weights, &
                           def_vec,fill_val,miss_val,valid_range,depend,depend_ref,lmodel,sub_model,P4)
  REAL (dp),                        POINTER       :: ptr(:,:,:,:)
  CHARACTER (len=*),                INTENT(IN)    :: src_dims
  INTEGER,                OPTIONAL, INTENT(in)    :: dt_update, miss, mis_dim, data_action, nIntrp, nAtOnce, nFuture
  INTEGER,                OPTIONAL, INTENT(in)    :: sub_model, valid_time, subset_index
  LOGICAL,                OPTIONAL, INTENT(in)    :: lcopy, ldealloc, lmodel, lauto
  INTEGER,                OPTIONAL, INTENT(inout) :: dims(:)
  REAL(dp),               OPTIONAL, INTENT(in)    :: mul,add,def_val,fill_val,miss_val,valid_range(2),weight,weights(:),def_vec(:)
  TYPE (input_var_list),  OPTIONAL, POINTER       :: var_ref, depend_ref
  CHARACTER (len=*),      OPTIONAL, INTENT(in)    :: var_name, alt_name, dt_unit, file_name, depend
  TYPE (input_file_list), OPTIONAL, POINTER       :: IFile
  REAL(dp),               OPTIONAL, POINTER       :: P4(:,:,:,:)

  TYPE (input_var_list), POINTER :: new
  INTEGER :: VDims(7)
  LOGICAL :: lOld

    VDims = (/-1,-1,-1,-1,1,1,1/)
    IF (PRESENT(dims))   VDims(1:4) = dims(1:4)
    IF (ASSOCIATED(ptr)) VDims(1:4) = SHAPE(ptr)

    new => InputVarAdd_gen(src_dims,4,VDims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time, &
           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val,                &
           miss_val,valid_range,weight,weights,lmodel,sub_model)
    IF (.NOT. lOld .AND. ASSOCIATED(ptr)) THEN
      IF (PRESENT(lcopy)) THEN
        IF (lcopy) new%dta(:,:,:,:,1,1,1) = ptr(:,:,:,:)
      ENDIF
      IF (PRESENT(ldealloc)) THEN
        IF (ldealloc) DEALLOCATE(ptr)
      ENDIF
    ENDIF
    ptr => new%dta(:,:,:,:,1,1,1)

    IF (PRESENT(var_ref)) var_ref => new
    IF (PRESENT(P4)) P4 => new%dta(:,:,:,:,1,1,1)

  END SUBROUTINE InputVar4dAdd

! Adds a new 5d variable to the variable list
  SUBROUTINE InputVar5dAdd(src_dims,ptr,var_name,IFile,dt_update,dt_unit,file_name,miss,mis_dim,data_action,dims,lcopy,ldealloc, &
                           valid_time,lauto,var_ref,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,weight,weights, &
                           def_vec,fill_val,miss_val,valid_range,depend,depend_ref,lmodel,sub_model)
  REAL (dp),                        POINTER       :: ptr(:,:,:,:,:)
  CHARACTER (len=*),                INTENT(IN)    :: src_dims
  INTEGER,                OPTIONAL, INTENT(in)    :: dt_update, miss, mis_dim, data_action, nIntrp, nAtOnce, nFuture
  INTEGER,                OPTIONAL, INTENT(in)    :: sub_model, valid_time, subset_index
  LOGICAL,                OPTIONAL, INTENT(in)    :: lcopy, ldealloc, lmodel, lauto
  INTEGER,                OPTIONAL, INTENT(inout) :: dims(:)
  REAL(dp),               OPTIONAL, INTENT(in)    :: mul,add,def_val,fill_val,miss_val,valid_range(2),weight,weights(:),def_vec(:)
  TYPE (input_var_list),  OPTIONAL, POINTER       :: var_ref, depend_ref
  CHARACTER (len=*),      OPTIONAL, INTENT(in)    :: var_name, alt_name, dt_unit, file_name, depend
  TYPE (input_file_list), OPTIONAL, POINTER       :: IFile

  TYPE (input_var_list), POINTER :: new
  INTEGER :: VDims(7)
  LOGICAL :: lOld

    VDims = (/-1,-1,-1,-1,-1,1,1/)
    IF (PRESENT(dims))   VDims(1:5) = dims(1:5)
    IF (ASSOCIATED(ptr)) VDims(1:5) = SHAPE(ptr)

    new => InputVarAdd_gen(src_dims,5,VDims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time, &
           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val,                &
           miss_val,valid_range,weight,weights,lmodel,sub_model)
    IF (.NOT. lOld .AND. ASSOCIATED(ptr)) THEN
      IF (PRESENT(lcopy)) THEN
        IF (lcopy) new%dta(:,:,:,:,:,1,1) = ptr(:,:,:,:,:)
      ENDIF
      IF (PRESENT(ldealloc)) THEN
        IF (ldealloc) DEALLOCATE(ptr)
      ENDIF
    ENDIF
    ptr => new%dta(:,:,:,:,:,1,1)

    IF (PRESENT(var_ref)) var_ref => new

  END SUBROUTINE InputVar5dAdd

! Adds a new 6d variable to the variable list
  SUBROUTINE InputVar6dAdd(src_dims,ptr,var_name,IFile,dt_update,dt_unit,file_name,miss,mis_dim,data_action,dims,lcopy,ldealloc, &
                           valid_time,lauto,var_ref,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,weight,weights, &
                           def_vec,fill_val,miss_val,valid_range,depend,depend_ref,lmodel,sub_model)
  REAL (dp),                        POINTER       :: ptr(:,:,:,:,:,:)
  CHARACTER (len=*),                INTENT(IN)    :: src_dims
  INTEGER,                OPTIONAL, INTENT(in)    :: dt_update, miss, mis_dim, data_action, nIntrp, nAtOnce, nFuture
  INTEGER,                OPTIONAL, INTENT(in)    :: sub_model, valid_time, subset_index
  LOGICAL,                OPTIONAL, INTENT(in)    :: lcopy, ldealloc, lmodel, lauto
  INTEGER,                OPTIONAL, INTENT(inout) :: dims(:)
  REAL(dp),               OPTIONAL, INTENT(in)    :: mul,add,def_val,fill_val,miss_val,valid_range(2),weight,weights(:),def_vec(:)
  TYPE (input_var_list),  OPTIONAL, POINTER       :: var_ref, depend_ref
  CHARACTER (len=*),      OPTIONAL, INTENT(in)    :: var_name, alt_name, dt_unit, file_name, depend
  TYPE (input_file_list), OPTIONAL, POINTER       :: IFile

  TYPE (input_var_list), POINTER :: new
  INTEGER :: VDims(7)
  LOGICAL :: lOld

    VDims = (/-1,-1,-1,-1,-1,-1,1/)
    IF (PRESENT(dims))   VDims(1:6) = dims(1:6)
    IF (ASSOCIATED(ptr)) VDims(1:6) = SHAPE(ptr)

    new => InputVarAdd_gen(src_dims,6,VDims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time, &
           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val,                &
           miss_val,valid_range,weight,weights,lmodel,sub_model)
    IF (.NOT. lOld .AND. ASSOCIATED(ptr)) THEN
      IF (PRESENT(lcopy)) THEN
        IF (lcopy) new%dta(:,:,:,:,:,:,1) = ptr(:,:,:,:,:,:)
      ENDIF
      IF (PRESENT(ldealloc)) THEN
        IF (ldealloc) DEALLOCATE(ptr)
      ENDIF
    ENDIF
    ptr => new%dta(:,:,:,:,:,:,1)

    IF (PRESENT(var_ref)) var_ref => new

  END SUBROUTINE InputVar6dAdd

! Adds a new 7d variable to the variable list
  SUBROUTINE InputVar7dAdd(src_dims,ptr,var_name,IFile,dt_update,dt_unit,file_name,miss,mis_dim,data_action,dims,lcopy,ldealloc, &
                           valid_time,lauto,var_ref,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,weight,weights, &
                           def_vec,fill_val,miss_val,valid_range,depend,depend_ref,lmodel,sub_model)
  REAL (dp),                        POINTER       :: ptr(:,:,:,:,:,:,:)
  CHARACTER (len=*),                INTENT(IN)    :: src_dims
  INTEGER,                OPTIONAL, INTENT(in)    :: dt_update, miss, mis_dim, data_action, nIntrp, nAtOnce, nFuture
  INTEGER,                OPTIONAL, INTENT(in)    :: sub_model, valid_time, subset_index
  LOGICAL,                OPTIONAL, INTENT(in)    :: lcopy, ldealloc, lmodel, lauto
  INTEGER,                OPTIONAL, INTENT(inout) :: dims(:)
  REAL(dp),               OPTIONAL, INTENT(in)    :: mul,add,def_val,fill_val,miss_val,valid_range(2),weight,weights(:),def_vec(:)
  TYPE (input_var_list),  OPTIONAL, POINTER       :: var_ref, depend_ref
  CHARACTER (len=*),      OPTIONAL, INTENT(in)    :: var_name, alt_name, dt_unit, file_name, depend 
  TYPE (input_file_list), OPTIONAL, POINTER       :: IFile

  TYPE (input_var_list), POINTER :: new
  INTEGER :: VDims(7)
  LOGICAL :: lOld

    VDims = (/-1,-1,-1,-1,-1,-1,-1/)
    IF (PRESENT(dims))   VDims(1:7) = dims(1:7)
    IF (ASSOCIATED(ptr)) VDims(1:7) = SHAPE(ptr)

    new => InputVarAdd_gen(src_dims,7,VDims,lOld,IFile,file_name,var_name,dt_update,dt_unit,miss,mis_dim,data_action,valid_time, &
           depend,depend_ref,lauto,alt_name,nIntrp,nAtOnce,nFuture,subset_index,mul,add,def_val,def_vec,fill_val,                &
           miss_val,valid_range,weight,weights,lmodel,sub_model)
    IF (.NOT. lOld .AND. ASSOCIATED(ptr)) THEN
      IF (PRESENT(lcopy)) THEN
        IF (lcopy) new%dta(:,:,:,:,:,:,:) = ptr(:,:,:,:,:,:,:)
      ENDIF
      IF (PRESENT(ldealloc)) THEN
        IF (ldealloc) DEALLOCATE(ptr)
      ENDIF
    ENDIF
    ptr => new%dta(:,:,:,:,:,:,:)

    IF (PRESENT(var_ref)) var_ref => new

  END SUBROUTINE InputVar7dAdd

! Create a new data buffer and link it to the end of a variables data chain
  SUBROUTINE InputVarNewBuffer(var,Sz,bufs,global,lShape,dta)
  TYPE (input_var_list),  POINTER :: var                     ! Variable to add buffer to
  INTEGER,             INTENT(in) :: Sz(7)                   ! Size of new buffer
  TYPE (input_curr),INTENT(inout) :: bufs                    ! Structure with this and last buffer to be updated
  TYPE (input_curr),INTENT(inout), OPTIONAL :: global        ! Previous buffers in global list
  LOGICAL,  OPTIONAL              :: lShape                  ! Does shape of buffer matter?
  REAL(dp), OPTIONAL,     POINTER :: dta(:,:,:,:,:,:,:)      ! Insert these data as buffer

  TYPE (input_data_list), POINTER :: new, ref, rprev
  REAL(dp),               POINTER :: new_dta(:,:,:,:,:,:,:)
  LOGICAL                         :: lDone

    lDone = .FALSE.
    IF (PRESENT(dta)) THEN
      IF (ASSOCIATED(dta)) THEN
        new_dta => dta
        NULLIFY(ref)
        lDone = .TRUE.
      ENDIF
    ENDIF
    IF (.NOT. lDone) THEN
      new_dta => InputDataGetRef(Sz,lShape,ref)
      ref%flags = IOR(ref%flags,INPUT_DATA_IN_USE)
    ENDIF

    IF (PRESENT(global)) THEN
      IF (ASSOCIATED(global%buf)) global%buf%flags = IAND(global%buf%flags,INPUT_ALLFLAGS - INPUT_DATA_IN_USE) ! Allow reuse
      global%buf_pre => global%buf
      global%buf     => ref
    ENDIF

    CALL InputDataNew(new)
    IF (ASSOCIATED(bufs%buf)) THEN
      bufs%buf_pre => bufs%buf
      rprev => bufs%buf
    ELSE
      rprev => var%buffers
      IF (ASSOCIATED(rprev)) THEN
        DO WHILE (ASSOCIATED(rprev%next))
          rprev => rprev%next
        ENDDO
      ENDIF
    ENDIF
    IF (ASSOCIATED(rprev)) THEN
      rprev%next => new
      IF (.NOT. PRESENT(global)) THEN
        ref => AllInputData
        DO WHILE (ASSOCIATED(ref))
          IF (ASSOCIATED(ref%dta,rprev%dta)) THEN
            ref%flags = IAND(ref%flags,INPUT_ALLFLAGS - INPUT_DATA_IN_USE)
            NULLIFY(ref)
          ELSE
            ref => ref%next
          ENDIF
        ENDDO
      ENDIF
    ELSE
      var%buffers => new
    ENDIF
    bufs%buf => new

    new%dta => new_dta

  END SUBROUTINE InputVarNewBuffer

! Push data into a variables buffer system. Individual subroutines for 0-7d data + a main one
  SUBROUTINE InputVarDataPush_int(var,buf,tim)
  TYPE (input_var_list), POINTER  :: var
  TYPE (input_data_list), POINTER :: buf
  INTEGER,  OPTIONAL, INTENT(in)  :: tim(4)

  TYPE (input_curr) :: curr
  INTEGER :: tt(4)
#ifdef HAVE_F2003
  TYPE (input_fcn_list), POINTER :: fcn
#endif
 
    IF (PRESENT(tim)) THEN
      tt = tim
    ELSE
      tt = model_time
      CALL CalendarTimeSub(tt,INT(model_dt,i8))
    ENDIF
    IF (ASSOCIATED(var%buffers%dta)) THEN
      var%time = tt
    ELSE
      var%Intrp%next%time = tt
    ENDIF
    IF (ASSOCIATED(var%saction)) THEN
      curr%acp  =  2
      curr%buf  => buf
      curr%mask => var%masks
#ifdef HAVE_F2003
      fcn => var%DataFcns%next ! First is a NoAction to be skipped
      DO WHILE (ASSOCIATED(fcn))
        CALL fcn%fcn(INPUT_PROC_DO_DATA_ACTION,var,curr)
        fcn => fcn%next
      ENDDO
#else
      CALL InputDoDataActions(var,curr,.FALSE.)
#endif

    ENDIF

  END SUBROUTINE InputVarDataPush_int

  SUBROUTINE InputVarDataPush0d(var,d0,tim)
  TYPE (input_var_list), POINTER  :: var
  REAL(dp),           INTENT(in)  :: d0
  INTEGER,  OPTIONAL, INTENT(in)  :: tim(4)
  TYPE (input_data_list), POINTER :: buf
    buf => InputVarGetNextBuffer(var,var%buffers,.FALSE.)
    buf%dta(:,:,:,:,:,:,:) = d0
    CALL InputVarDataPush_int(var,buf,tim)
  END SUBROUTINE InputVarDataPush0d
  SUBROUTINE InputVarDataPush1d(var,d1,tim)
  TYPE (input_var_list), POINTER  :: var
  REAL(dp),           INTENT(in)  :: d1(:)
  INTEGER,  OPTIONAL, INTENT(in)  :: tim(4)
  TYPE (input_data_list), POINTER :: buf
    buf => InputVarGetNextBuffer(var,var%buffers,.FALSE.)
    buf%dta(:,1,1,1,1,1,1) = d1(:)
    CALL InputVarDataPush_int(var,buf,tim)
  END SUBROUTINE InputVarDataPush1d
  SUBROUTINE InputVarDataPush2d(var,d2,tim)
  TYPE (input_var_list), POINTER  :: var
  REAL(dp),           INTENT(in)  :: d2(:,:)
  INTEGER,  OPTIONAL, INTENT(in)  :: tim(4)
  TYPE (input_data_list), POINTER :: buf
    buf => InputVarGetNextBuffer(var,var%buffers,.FALSE.)
    buf%dta(:,:,1,1,1,1,1) = d2(:,:)
    CALL InputVarDataPush_int(var,buf,tim)
  END SUBROUTINE InputVarDataPush2d
  SUBROUTINE InputVarDataPush3d(var,d3,tim)
  TYPE (input_var_list), POINTER  :: var
  REAL(dp),           INTENT(in)  :: d3(:,:,:)
  INTEGER,  OPTIONAL, INTENT(in)  :: tim(4)
  TYPE (input_data_list), POINTER :: buf
    buf => InputVarGetNextBuffer(var,var%buffers,.FALSE.)
    buf%dta(:,:,:,1,1,1,1) = d3(:,:,:)
    CALL InputVarDataPush_int(var,buf,tim)
  END SUBROUTINE InputVarDataPush3d
  SUBROUTINE InputVarDataPush4d(var,d4,tim)
  TYPE (input_var_list), POINTER  :: var
  REAL(dp),           INTENT(in)  :: d4(:,:,:,:)
  INTEGER,  OPTIONAL, INTENT(in)  :: tim(4)
  TYPE (input_data_list), POINTER :: buf
    buf => InputVarGetNextBuffer(var,var%buffers,.FALSE.)
    buf%dta(:,:,:,:,1,1,1) = d4(:,:,:,:)
    CALL InputVarDataPush_int(var,buf,tim)
  END SUBROUTINE InputVarDataPush4d
  SUBROUTINE InputVarDataPush5d(var,d5,tim)
  TYPE (input_var_list), POINTER  :: var
  REAL(dp),           INTENT(in)  :: d5(:,:,:,:,:)
  INTEGER,  OPTIONAL, INTENT(in)  :: tim(4)
  TYPE (input_data_list), POINTER :: buf
    buf => InputVarGetNextBuffer(var,var%buffers,.FALSE.)
    buf%dta(:,:,:,:,:,1,1) = d5(:,:,:,:,:)
    CALL InputVarDataPush_int(var,buf,tim)
  END SUBROUTINE InputVarDataPush5d
  SUBROUTINE InputVarDataPush6d(var,d6,tim)
  TYPE (input_var_list), POINTER  :: var
  REAL(dp),           INTENT(in)  :: d6(:,:,:,:,:,:)
  INTEGER,  OPTIONAL, INTENT(in)  :: tim(4)
  TYPE (input_data_list), POINTER :: buf
    buf => InputVarGetNextBuffer(var,var%buffers,.FALSE.)
    buf%dta(:,:,:,:,:,:,1) = d6(:,:,:,:,:,:)
    CALL InputVarDataPush_int(var,buf,tim)
  END SUBROUTINE InputVarDataPush6d
  SUBROUTINE InputVarDataPush7d(var,d7,tim)
  TYPE (input_var_list), POINTER  :: var
  REAL(dp),           INTENT(in)  :: d7(:,:,:,:,:,:,:)
  INTEGER,  OPTIONAL, INTENT(in)  :: tim(4)
  TYPE (input_data_list), POINTER :: buf
    buf => InputVarGetNextBuffer(var,var%buffers,.FALSE.)
    buf%dta(:,:,:,:,:,:,:) = d7(:,:,:,:,:,:,:)
    CALL InputVarDataPush_int(var,buf,tim)
  END SUBROUTINE InputVarDataPush7d

  SUBROUTINE InputVarGetDimForOutput(var,ldims,gdims,dimnames)
  TYPE (input_var_list), POINTER :: var
  INTEGER,           INTENT(out) :: ldims(4), gdims(4)
  CHARACTER(len=64), INTENT(out) :: dimnames(4)

  TYPE (input_dim_list), POINTER :: curr_dim
  INTEGER :: i

    ldims(:) = 1
    gdims(:) = 1
    dimnames(:) = ''
 
    i = 1
    curr_dim => var%dims
    DO WHILE (ASSOCIATED(curr_dim))
      ldims(i) = curr_dim%dim_data%size_local
      gdims(i) = curr_dim%dim_data%size_global
      dimnames(i) = TRIM(curr_dim%dim_data%name_dim) ! This may require more work for ECHAM variable declared on equivalent dims
      i = i + 1
      curr_dim => curr_dim%next
    ENDDO
    
  END SUBROUTINE InputVarGetDimForOutput

#if defined(INPUT_IN_ECHAM) || defined(INPUT_IN_ICON)
! Add variable(s) to an output stream
  SUBROUTINE InputVar2Output(OStream,IFile,var,start_code,codes)
  USE mo_linked_list, ONLY : t_stream
  USE mo_memory_base, ONLY : add_stream_element
  USE mo_netcdf,      ONLY : max_dim_name

  TYPE (t_stream),          POINTER           :: OStream
  TYPE (input_file_list),   POINTER, OPTIONAL :: IFile 
  TYPE (input_var_list) ,   POINTER, OPTIONAL :: var
  INTEGER, INTENT(in)              , OPTIONAL :: start_code
  INTEGER, INTENT(in), DIMENSION(:), OPTIONAL :: codes

  TYPE (input_file_list), POINTER :: curr_file
  TYPE (input_var_list),  POINTER :: curr_var
  TYPE (input_dim_list),  POINTER :: curr_dim
  REAL(dp), POINTER :: dum1(:), dum2(:,:), dum3(:,:,:), dum4(:,:,:,:)
  INTEGER :: ldims(7), gdims(7)
  CHARACTER (max_dim_name) :: dimn(7)

  INTEGER :: code, iCurr, i
  LOGICAL :: lCont

    IF (PRESENT(IFile)) THEN
      curr_file => IFile
    ELSE
      curr_file => AllInputFiles
      lCont = ASSOCIATED(curr_file)
      IF (lCont) lCont = IAND(curr_file%flags,INPUT_FILE_INITIAL) /= 0
      DO WHILE (lCont)
        curr_file => curr_file%next
        lCont = ASSOCIATED(curr_file)
        IF (lCont) lCont = IAND(curr_file%flags,INPUT_FILE_INITIAL) /= 0
      ENDDO
    ENDIF
    IF (.NOT. ASSOCIATED(curr_file)) RETURN

    IF (PRESENT(var)) THEN
      IF (PRESENT(IFile) .AND. .NOT. ASSOCIATED(var%parent_file,IFile)) &
      CALL local_error('InputAdd2Output','Variable specified for output not within specificed file')
      curr_var => var
    ELSE
      curr_var => AllInputVars
      DO WHILE (ASSOCIATED(curr_var) .AND. .NOT. ASSOCIATED(curr_var%parent_file,curr_file))
        curr_var => curr_var%next
      ENDDO
      IF (.NOT. ASSOCIATED(curr_var)) RETURN
    ENDIF

    iCurr = 1
    DO WHILE (ASSOCIATED(curr_var))

      ! Determine grib code
      IF (PRESENT(codes)) THEN
        IF (SIZE(codes) >= iCurr) THEN
          code = codes(iCurr)
        ELSE
          code = iCurr - SIZE(codes) + codes(SIZE(codes))
        ENDIF
      ELSE
        IF (PRESENT(start_code)) THEN
          code = iCurr + start_code - 1
        ELSE
          code = iCurr
        ENDIF
      ENDIF

      ! Extract dimension information
      ldims = curr_var%edims
      curr_dim => AllInputEqDims%src_dims !curr_var%dims ! Dirty hack to give echam the dimensions it want's!
      DO i=1,curr_var%ndim
        gdims(i) = curr_dim%dim_data%size_global
        dimn(i) = curr_dim%dim_data%name_dim
        curr_dim => curr_dim%next
      ENDDO

      ! Insert variable in output stream
#ifdef INPUT_IN_ECHAM
      dum4 => curr_var%dta(:,:,:,:,1,1,1)
      SELECT CASE (curr_var%ndim)
        CASE (1)
          CALL add_stream_element(OStream,TRIM(curr_var%name_var),dum1,p4=dum4,code=code,lrerun=.FALSE., &
          ldims=ldims(1:1),gdims=gdims(1:1),dimnames=dimn(1:1),lpost=.TRUE.)
        CASE (2)
          CALL add_stream_element(OStream,TRIM(curr_var%name_var),dum2,p4=dum4,code=code,lrerun=.FALSE., &
          ldims=ldims(1:2),gdims=gdims(1:2),dimnames=dimn(1:2),lpost=.TRUE.)
        CASE (3)
          CALL add_stream_element(OStream,TRIM(curr_var%name_var),dum3,p4=dum4,code=code,lrerun=.FALSE., &
          ldims=ldims(1:3),gdims=gdims(1:3),dimnames=dimn(1:3),lpost=.TRUE.)
        CASE (4)
          CALL add_stream_element(OStream,TRIM(curr_var%name_var),dum4,p4=dum4,code=code,lrerun=.FALSE., &
          ldims=ldims(1:4),gdims=gdims(1:4),dimnames=dimn(1:4),lpost=.TRUE.)
        CASE DEFAULT
          WRITE(message_text,'(3a,i1,a)') 'Warning: Variables to add to output must have 1-4 dimensions. ', &
                                    TRIM(curr_var%name_var),' has ',curr_var%ndim,' dimensions and will not be added to the output'
          IF (IAND(mo_debug,INPUT_DBG_MASK) /= 0) CALL local_message('InputAdd2Output',message_text)
      END SELECT
#endif      
      ! Next variable
      IF (PRESENT(var)) THEN
        NULLIFY(curr_var)
      ELSE
        curr_var => curr_var%next
        DO WHILE (ASSOCIATED(curr_var) .AND. .NOT. ASSOCIATED(curr_var%parent_file,curr_file))
        curr_var => curr_var%next
        ENDDO
        IF (.NOT. ASSOCIATED(curr_var) .AND. .NOT. PRESENT(IFile)) THEN
          curr_file => curr_file%next
          lCont = ASSOCIATED(curr_file)
          IF (lCont) lCont = IAND(curr_file%flags,INPUT_FILE_INITIAL) /= 0
          DO WHILE (lCont)
            curr_file => curr_file%next
            lCont = ASSOCIATED(curr_file)
            IF (lCont) lCont = IAND(curr_file%flags,INPUT_FILE_INITIAL) /= 0
          ENDDO
          IF (ASSOCIATED(curr_file)) THEN
            curr_var => AllInputVars
            DO WHILE (ASSOCIATED(curr_var) .AND. .NOT. ASSOCIATED(curr_var%parent_file,curr_file))
            curr_var  => curr_var%next
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      iCurr = iCurr + 1

    ENDDO

  END SUBROUTINE InputVar2Output
#endif

! Return if a variable is marked as initialized
  LOGICAL FUNCTION InputVarIsInitialized(var)
  TYPE (input_var_list), POINTER :: var

    InputVarIsInitialized = .FALSE. 
    IF (ASSOCIATED(var)) InputVarIsInitialized = IAND(var%stat,INPUT_VAR_INITIALIZED) /= 0

  END FUNCTION InputVarIsInitialized

! Return if a variable is marked as present in currently opened file
  LOGICAL FUNCTION InputVarIsPresent(var)
  TYPE (input_var_list), POINTER :: var

    InputVarIsPresent = .FALSE.
    IF (ASSOCIATED(var)) InputVarIsPresent = IAND(var%stat,INPUT_VAR_EXIST) /= 0

  END FUNCTION InputVarIsPresent

! Return if a variable has been updated in this time step or not
  LOGICAL FUNCTION InputVarIsUpdated(var)
  TYPE (input_var_list), POINTER :: var

    InputVarIsUpdated = .FALSE.
    IF (ASSOCIATED(var)) InputVarIsUpdated = IAND(var%stat,INPUT_VAR_UPDATED) /= 0

  END FUNCTION InputVarIsUpdated

! Return if a variable has been read from file in this time step or not
  LOGICAL FUNCTION InputVarIsRead(var)
  TYPE (input_var_list), POINTER :: var

    InputVarIsRead = .FALSE.
    IF (ASSOCIATED(var)) InputVarIsRead = IAND(var%stat,INPUT_VAR_READ) /= 0

  END FUNCTION InputVarIsRead

! Return if a variable will be updated in the next call to InputUpdate or not
  LOGICAL FUNCTION InputVarWillBeUpdated(var)
  TYPE (input_var_list), POINTER :: var

    InputVarWillBeUpdated = .FALSE.
    IF (ASSOCIATED(var)) InputVarWillBeUpdated = var%next_update <= 0 

  END FUNCTION InputVarWillBeUpdated

! Return if a variable will be read during the next call to InputUpdate or not
  LOGICAL FUNCTION InputVarWillBeRead(var)
  TYPE (input_var_list), POINTER :: var

  INTEGER :: tme(4)

    IF (ASSOCIATED(var)) THEN
      tme = model_time
      CALL CalendarTimeAdd(tme,INT(model_dt,i8))
      InputVarWillBeRead = InputVarWillBeUpdated(var) .AND. (CalendarTimeDiff(var%time,tme) >= 0) 
    ELSE
      InputVarWillBeRead = .FALSE.
    ENDIF

  END FUNCTION InputVarWillBeRead

! Return if a variable is read from an external file
  LOGICAL FUNCTION InputVarIsExternal(var)
  TYPE (input_var_list), POINTER :: var

    InputVarIsExternal = .FALSE.

    IF (.NOT. ASSOCIATED(var)) RETURN
    IF (.NOT. ASSOCIATED(var%parent_file)) RETURN
    IF (ASSOCIATED(var%parent_file,DummyFile)) RETURN
    InputVarIsExternal = .TRUE. 

  END FUNCTION InputVarIsExternal

! Return number of InputUpdate calls since last read of variable
  INTEGER FUNCTION InputVarLastRead(var)
  TYPE (input_var_list), POINTER :: var

    InputVarLastRead = HUGE(mo_debug)
    IF (ASSOCIATED(var)) InputVarLastRead = var%nLastRead

  END FUNCTION InputVarLastRead

! Return number of InputUpdate calls since last update of variable
  INTEGER FUNCTION InputVarLastUpdate(var)
  TYPE (input_var_list), POINTER :: var

    InputVarLastUpdate = HUGE(mo_debug)
    IF (ASSOCIATED(var)) InputVarLastUpdate = var%nLastUpdate

  END FUNCTION InputVarLastUpdate

! Resets read and update flags of a variable
  SUBROUTINE InputVarUpdated(var)
  TYPE (input_var_list), POINTER :: var

    IF (ASSOCIATED(var)) var%stat = IAND(var%stat,INPUT_ALLFLAGS - INPUT_VAR_READ - INPUT_VAR_UPDATED)

  END SUBROUTINE InputVarUpdated

  INTEGER FUNCTION InputVarUserDataSet0D(var,dta,no)
  TYPE (input_var_list), POINTER :: var
  REAL (dp),         INTENT(in)  :: dta
  INTEGER, OPTIONAL, INTENT(in)  :: no

  REAL(dp), POINTER :: org_dta(:,:,:,:,:,:,:)
  TYPE (input_data_list), POINTER :: new
  INTEGER :: i, id

    i = -1
    IF (PRESENT(no)) THEN
      IF (no==-1) THEN
        org_dta => InputDataGetSpecified(var%aux_data,INPUT_DATA_NUDGE_FLD,0,new)
        id = INPUT_DATA_NUDGE_FLD
      ELSE
        org_dta => InputDataGetSpecified(var%aux_data,INPUT_DATA_USER,no*INPUT_DATA_NO_STEP,new)
        id = INPUT_DATA_USER + i*INPUT_DATA_NO_STEP
      ENDIF
      IF (ASSOCIATED(org_dta)) THEN
        new%dta = dta
        InputVarUserDataSet0D = no
        RETURN
      ENDIF
    ENDIF
    i = GetNextUserDataNumber(var)
    CALL InputDataNew(new)
    new%flags    =  id
    new%next     => var%aux_data
    var%aux_data => new
    ALLOCATE(new%dta(1,1,1,1,1,1,1))
    new%dta      =  dta
    InputVarUserDataSet0D = i

  END FUNCTION InputVarUserDataSet0D

  INTEGER FUNCTION InputVarUserDataSet1D(var,dta,no)
  TYPE (input_var_list), POINTER :: var
  REAL (dp), TARGET, INTENT(in)  :: dta(:)
  INTEGER, OPTIONAL, INTENT(in)  :: no

  REAL(dp), POINTER :: org_dta(:,:,:,:,:,:,:)
  TYPE (input_data_list), POINTER :: new
  INTEGER :: i, id

    NULLIFY(new)
    i = -1
    IF (PRESENT(no)) THEN
      IF (no==-1) THEN
        org_dta => InputDataGetSpecified(var%aux_data,INPUT_DATA_NUDGE_FLD,0,new)
        id = INPUT_DATA_NUDGE_FLD
      ELSE
        org_dta => InputDataGetSpecified(var%aux_data,INPUT_DATA_USER,no*INPUT_DATA_NO_STEP,new)
        id = INPUT_DATA_USER + i*INPUT_DATA_NO_STEP
      ENDIF
      IF (ASSOCIATED(org_dta)) THEN
        i = no
        IF ((SIZE(org_dta,1) /= SIZE(dta)) .OR. (SIZE(org_dta,2) /= 1) .OR. &
            (SIZE(org_dta,3) /= 1)         .OR. (SIZE(org_dta,4) /= 1) .OR. &
            (SIZE(org_dta,5) /= 1)         .OR. (SIZE(org_dta,6) /= 1) .OR. &
            (SIZE(org_dta,7) /= 1)) THEN
          DEALLOCATE(new%dta)
        ELSE
          new%dta(:,1,1,1,1,1,1) = dta
          InputVarUserDataSet1D = i
          RETURN
        ENDIF
        i = no
      ENDIF
    ENDIF
    IF (.NOT. ASSOCIATED(new)) THEN
      i = GetNextUserDataNumber(var)
      CALL InputDataNew(new)
      new%flags    =  id
      new%next     => var%aux_data
      var%aux_data => new
    ENDIF
    ALLOCATE(new%dta(SIZE(dta),1,1,1,1,1,1))
    new%dta(:,1,1,1,1,1,1) = dta
    InputVarUserDataSet1D = i

  END FUNCTION InputVarUserDataSet1D

  INTEGER FUNCTION InputVarUserDataSet7D(var,dta,no)
  TYPE (input_var_list), POINTER :: var
  REAL (dp), TARGET, INTENT(in)  :: dta(:,:,:,:,:,:,:)
  INTEGER, OPTIONAL, INTENT(in)  :: no

  REAL(dp), POINTER :: org_dta(:,:,:,:,:,:,:)
  TYPE (input_data_list), POINTER :: new
  INTEGER :: i, id

    i = -1
    IF (PRESENT(no)) THEN
      IF (no==-1) THEN
        org_dta => InputDataGetSpecified(var%aux_data,INPUT_DATA_NUDGE_FLD,0,new)
        id = INPUT_DATA_NUDGE_FLD
      ELSE
        org_dta => InputDataGetSpecified(var%aux_data,INPUT_DATA_USER,no*INPUT_DATA_NO_STEP,new)
        id = INPUT_DATA_USER + i*INPUT_DATA_NO_STEP
      ENDIF
      IF (ASSOCIATED(org_dta)) THEN
        new%dta => dta
        InputVarUserDataSet7D = no
        RETURN
      ENDIF
    ENDIF
    i = GetNextUserDataNumber(var)
    CALL InputDataNew(new)
    new%flags    =  id
    new%next     => var%aux_data
    var%aux_data => new
    new%dta      => dta
    InputVarUserDataSet7D = i

  END FUNCTION InputVarUserDataSet7D

! Return user data of a variable
  FUNCTION InputVarUserDataGet(var,no)
  REAL (dp), POINTER :: InputVarUserDataGet(:,:,:,:,:,:,:)
  TYPE (input_var_list), POINTER :: var
  INTEGER, OPTIONAL, INTENT(in)  :: no

  INTEGER :: i

    i = -1
    IF (PRESENT(no)) i = no
    InputVarUserDataGet => InputDataGetSpecified(var%aux_data,INPUT_DATA_USER,i*INPUT_DATA_NO_STEP)

  END FUNCTION InputVarUserDataGet

! Return time between the read fields for interpolation
  INTEGER FUNCTION InputVarGetInterpolStep(var)
  TYPE (input_var_list), POINTER :: var

    InputVarGetInterpolStep = 0
    IF (ASSOCIATED(var%Intrp)) &
      InputVarGetInterpolStep = INT(CalendarTimeDiff(var%Intrp%prev%time,var%Intrp%time))

  END FUNCTION InputVarGetInterpolStep

! Returns the string value of a named file or variable attribute
  FUNCTION InputAttrGet(attr,IFile,var,stat)
  CHARACTER(len=256)           :: InputAttrGet
  CHARACTER(len=*), INTENT(in) :: Attr
  TYPE (input_file_list), POINTER, OPTIONAL :: IFile
  TYPE (input_var_list),  POINTER, OPTIONAL :: var
  INTEGER,          INTENT(out),   OPTIONAL :: stat

  TYPE (input_file_list), POINTER :: ActFile
  CHARACTER(len=256) :: tmp
  INTEGER :: itmp

    InputAttrGet = ''
    IF (PRESENT(stat)) stat = 3
    NULLIFY(ActFile)
    IF (PRESENT(var)) THEN
      IF (ASSOCIATED(var)) ActFile => var%parent_file
    ENDIF
    IF (.NOT. ASSOCIATED(ActFile) .AND. PRESENT(IFile)) THEN
      ActFile => IFile
    ENDIF
    IF (ASSOCIATED(ActFile)) THEN
      IF (IAND(ActFile%stat,INPUT_FILE_OPEN+INPUT_FILE_CHECK_VAR) /= INPUT_FILE_OPEN+INPUT_FILE_CHECK_VAR) &
        CALL InputFileOpen(ActFile,ActFile%time_file(1),ActFile%time_file(2),ActFile%time_file(3),check=INPUT_FILE_CHECK_VAR)
      IF (lio) THEN
#ifdef HAVE_F2003
        tmp = TRIM(ActFile%ft%FcnGetAttr(attr,ActFile,var,itmp))
#else
        tmp = TRIM(NCFileGetAttr(attr,ActFile,var,itmp))
#endif
      ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      CALL p_bcast(tmp,io_pe)
      CALL p_bcast(itmp,io_pe)
#endif
      IF (PRESENT(stat)) stat = itmp
      InputAttrGet = TRIM(tmp)
    ENDIF

  END FUNCTION InputAttrGet

! ---------------------------------------------------------
!  Check contents of files against previously defined data
! ---------------------------------------------------------

! Check or set sizes of all model dimensions found in file. .FALSE. is returned if any mismatches are found,
! message_text contains info
  RECURSIVE LOGICAL FUNCTION InputCheckDimSizes(IFile,sub_model)
  TYPE (input_file_list),    POINTER :: IFile
  INTEGER,                INTENT(in) :: sub_model

  TYPE (input_fileobj_list), POINTER :: curr_fdim
  TYPE (input_dim_list),     POINTER :: curr_dim, new_dim, prev_dim, file_dim
  CHARACTER (len=64)                 :: TmpStr, DimName
  LOGICAL :: lCont
  INTEGER :: Sz, i, tst(3)

    InputCheckDimSizes = .TRUE.
    message_text = ''
    lCont = .TRUE.
    IF (lio) THEN
      curr_dim => AllInputDims
      DO WHILE (ASSOCIATED(curr_dim))
        curr_dim%dim_data%flags = IAND(curr_dim%dim_data%flags,INPUT_ALLFLAGS - INPUT_DIM_EXIST - INPUT_DIM_FILE)
        curr_dim => curr_dim%next
      ENDDO 
      NULLIFY(prev_dim)
      file_dim  => IFile%dims
      curr_fdim => IFile%fdims
      DO WHILE (ASSOCIATED(curr_fdim))
        TmpStr = ''
        curr_dim  => InputDimGetRef(TRIM(curr_fdim%name_obj),sub_model=sub_model,search=file_dim) ! Check for dim remappings
        IF (.NOT. ASSOCIATED(curr_dim)) THEN
          curr_dim  => InputDimGetRef(TRIM(curr_fdim%name_obj),sub_model=sub_model)
        ELSE
          TmpStr = curr_dim%file_dim_name
          curr_dim  => InputDimGetRef(TRIM(curr_dim%dim_data%name_dim),sub_model=sub_model) ! Model version of remaped dim
        ENDIF
        IF (ASSOCIATED(curr_dim)) THEN
          Sz = curr_fdim%size_global
          IF ((curr_dim%dim_data%size_global < 0) .OR. (IAND(curr_dim%dim_data%flags,INPUT_DIM_UNLIM) /= 0)) THEN
            curr_dim%dim_data%size_global = Sz
            IF (curr_dim%dim_data%size_local < 0) curr_dim%dim_data%size_local = Sz
            IF (curr_dim%dim_data%local_lo   < 0) curr_dim%dim_data%local_lo   =  1
            IF (curr_dim%dim_data%local_hi   < 0) curr_dim%dim_data%local_hi   = Sz
            IF ((curr_dim%dim_data%chunk_lo  < 0) .AND. lio) curr_dim%dim_data%chunk_lo   =  1
            IF ((curr_dim%dim_data%chunk_hi  < 0) .AND. lio) curr_dim%dim_data%chunk_hi   = Sz
          ELSE
            IF (IAND(curr_dim%dim_data%flags,INPUT_DIM_EXIST) == 0) THEN
              IF ((curr_dim%dim_data%size_global /= Sz) .AND. (Sz /= 1)) THEN ! DimSize 1 is still allowed,check further with vars
                WRITE(message_text,*) TRIM(message_text),',_',TRIM(curr_dim%dim_data%name_dim),'=',Sz, &
                                      '_(expected:_',curr_dim%dim_data%size_global,')'
                lCont = .FALSE.
              ENDIF
            ELSE
              NULLIFY(curr_dim) ! This dimension was already determined and this is an unlucky aliasing of dim. names
            ENDIF
          ENDIF
        ENDIF
        IF (ASSOCIATED(curr_dim)) THEN
!          IF (IAND(curr_dim%dim_data%flags,INPUT_DIM_UNLIM) /= 0) IFile%nRec = Sz
          IF (ASSOCIATED(curr_fdim%points) .AND. (IAND(curr_dim%dim_data%flags,INPUT_DIM_CYCLIC) /= 0)) THEN
            IF (curr_dim%dim_data%origin == UNLIKELY_RVAL) curr_dim%dim_data%origin = curr_fdim%points(1)
          ENDIF
          CALL InputDimNewCtl(new_dim)
          IF (ASSOCIATED(prev_dim)) THEN
            prev_dim%next => new_dim
          ELSE
            IFile%dims    => new_dim
          ENDIF
          prev_dim => new_dim
          new_dim%file_dim_name = TmpStr 
          new_dim%dim_data => curr_dim%dim_data
          new_dim%src_dims => curr_dim%src_dims
          new_dim%lo       =  curr_dim%lo
          new_dim%hi       =  curr_dim%hi
          IF (IAND(curr_fdim%flags,INPUT_FILEDIM_DIR_MASK) == INPUT_FILEDIM_DIR_DEC) new_dim%flags = INPUT_DIM_REVERSE
          new_dim%flags = IOR(new_dim%flags,INPUT_DIM_EXIST+INPUT_DIM_FILE)
          new_dim%dim_data%flags = IOR(new_dim%dim_data%flags,INPUT_DIM_EXIST+INPUT_DIM_FILE)
          new_dim%fid = curr_fdim%id
        ELSE ! Make a list of file dimensions which are not part of the model - may be somehow aggregated
          CALL InputDimNew(new_dim)
          IF (ASSOCIATED(prev_dim)) THEN
            prev_dim%next => new_dim
          ELSE
            IFile%dims    => new_dim
          ENDIF
          prev_dim => new_dim
          new_dim%dim_data%name_dim    = TRIM(curr_fdim%name_obj)
          new_dim%dim_data%size_global = curr_fdim%size_global
          new_dim%dim_data%size_local  = curr_fdim%size_global
          new_dim%dim_data%local_lo    =                     1
          new_dim%dim_data%local_hi    = curr_fdim%size_global
          new_dim%dim_data%chunk_lo    =                     1
          new_dim%dim_data%chunk_hi    = curr_fdim%size_global
          new_dim%dim_data%flags       = INPUT_DIM_EXIST + INPUT_DIM_NON_MODEL + INPUT_DIM_FILE
          new_dim%fid = curr_fdim%id
          IF (IAND(curr_fdim%flags,INPUT_FILEDIM_DIR_MASK)==INPUT_FILEDIM_DIR_DEC) &
            new_dim%dim_data%flags = new_dim%dim_data%flags + INPUT_DIM_REVERSE
        ENDIF
        curr_fdim => curr_fdim%next
      ENDDO
      IF (ASSOCIATED(file_dim)) CALL InputDimDone(file_dim,.FALSE.,.TRUE.)
    ENDIF
    IF (.NOT. lCont) THEN
      CALL RmSpc(message_text)
      message_text = 'Dimension(s) size mismatch:'//TRIM(message_text(2:))
    ENDIF

#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      CALL p_bcast(lCont,io_pe)
#endif

    InputCheckDimSizes = lCont

    ! Distribute results to other processors
    IF (n_pe >= 0) THEN
      IF (lCont) THEN
        IF (lio) THEN
          curr_dim => IFile%dims
          Sz = 0
          DO WHILE(ASSOCIATED(curr_dim))
            Sz = Sz + 1
            curr_dim => curr_dim%next
          ENDDO
        ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
        CALL p_bcast(Sz,io_pe)
#endif
        NULLIFY(prev_dim)
        curr_dim => IFile%dims
        DO i=1,Sz
          IF (.NOT. lio) THEN
            CALL InputDimNewCtl(new_dim)
            IF (ASSOCIATED(prev_dim)) THEN
              prev_dim%next => new_dim
            ELSE
              IFile%dims    => new_dim
            ENDIF
            prev_dim => new_dim
          ELSE
            DimName = TRIM(curr_dim%dim_data%name_dim)
            new_dim => curr_dim
          ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
          CALL p_bcast(DimName,io_pe)
#endif
          IF (lio) THEN
            curr_dim => curr_dim%next
          ELSE
            curr_dim  => InputDimGetRef(DimName,sub_model=sub_model)
            IF (.NOT. ASSOCIATED(curr_dim)) CALL InputDimAdd(DimName,dim_ref=curr_dim,intr=.true.) ! Declare unkown file dim.glob.
            new_dim%dim_data => curr_dim%dim_data
            new_dim%src_dims => curr_dim%src_dims
            new_dim%lo       =  curr_dim%lo
            new_dim%hi       =  curr_dim%hi
            new_dim%flags = IOR(IAND(curr_dim%flags,INPUT_ALLFLAGS - INPUT_DIM_TYPE),INPUT_DIM_FILE)
          ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
          CALL p_bcast(new_dim%file_dim_name       ,io_pe)
          CALL p_bcast(new_dim%flags               ,io_pe)
          CALL p_bcast(new_dim%fid                 ,io_pe)
          CALL p_bcast(new_dim%dim_data%flags      ,io_pe)
          CALL p_bcast(new_dim%dim_data%size_global,io_pe)
          CALL p_bcast(new_dim%dim_data%origin     ,io_pe)
          tst = (/new_dim%dim_data%chunk_lo,new_dim%dim_data%chunk_hi,0/)
          CALL p_bcast(tst,io_pe)
          IF (.NOT. lio) THEN
            IF (new_dim%dim_data%chunk_lo   < 0) new_dim%dim_data%chunk_lo   = tst(1)
            IF (new_dim%dim_data%chunk_hi   < 0) new_dim%dim_data%chunk_hi   = tst(2)
          ENDIF
          tst = (/new_dim%dim_data%size_local,new_dim%dim_data%local_lo,new_dim%dim_data%local_hi/)
          CALL p_bcast(tst ,io_pe)
          IF (.NOT. lio) THEN
            IF (new_dim%dim_data%size_local < 0) new_dim%dim_data%size_local = tst(1)
            IF (new_dim%dim_data%local_lo   < 0) new_dim%dim_data%local_lo   = tst(2)
            IF (new_dim%dim_data%local_hi   < 0) new_dim%dim_data%local_hi   = tst(3)
          ENDIF
#endif
        ENDDO
      ENDIF  
    ENDIF

  END FUNCTION InputCheckDimSizes

! Adds a new action to the action list
  SUBROUTINE InputAddDataAction(var,act,curr &
#ifdef HAVE_F2003
  ,NewAction &
#endif
  )
  TYPE (input_var_list),   POINTER :: var
  INTEGER,              INTENT(in) :: act(:)
  TYPE (input_curr), INTENT(inout) :: curr
#ifdef HAVE_F2003
  TYPE (input_fcn_list),   POINTER, OPTIONAL :: NewAction
#endif

  INTEGER, POINTER :: actions(:)
#ifdef HAVE_F2003
  TYPE (input_fcn_list), POINTER :: act_fcn

    act_fcn => InputFcnAdd(var%DataFcns)
    IF (.NOT. ASSOCIATED(var%DataFcns)) var%DataFcns => act_fcn
    IF (PRESENT(NewAction)) NewAction => act_fcn
    SELECT CASE (ABS(act(1)))
      CASE (INPUT_VAR_DO_NOACTION);          act_fcn%fcn => DataProcNoAction
      CASE (INPUT_VAR_DO_MESSAGE) ;          act_fcn%fcn => DataProcMessage
      CASE (INPUT_VAR_DO_TOREAD);            act_fcn%fcn => DataProcRead
      CASE (INPUT_VAR_DO_MULTI_READ);        act_fcn%fcn => DataProcReadMulti
      CASE (INPUT_VAR_DO_PERMUTE);           act_fcn%fcn => DataProcPermute
      CASE (INPUT_VAR_DO_PACK);              act_fcn%fcn => DataProcPack
      CASE (INPUT_VAR_DO_UNPACK);            act_fcn%fcn => DataProcUnpack
      CASE (INPUT_VAR_DO_EQUIV);             act_fcn%fcn => DataProcEquiv
      CASE (INPUT_VAR_DO_COPY);              act_fcn%fcn => DataProcCopy
      CASE (INPUT_VAR_DO_REVERSE);           act_fcn%fcn => DataProcReverse
      CASE (INPUT_VAR_DO_SWAP);              act_fcn%fcn => DataProcSwap
      CASE (INPUT_VAR_DO_SUM);               act_fcn%fcn => DataProcSum
      CASE (INPUT_VAR_DO_AVG);               act_fcn%fcn => DataProcAvg
      CASE (INPUT_VAR_DO_NORM);              act_fcn%fcn => DataProcNormalize
      CASE (INPUT_VAR_DO_FILL_BAD);          act_fcn%fcn => DataProcFillBad
      CASE (INPUT_VAR_DO_FILL_RANGE);        act_fcn%fcn => DataProcFillOutOfRange
      CASE (INPUT_VAR_DO_RESCALE);           act_fcn%fcn => DataProcRescale
      CASE (INPUT_VAR_DO_DIVTIME);           act_fcn%fcn => DataProcTimeDiv
      CASE (INPUT_VAR_DO_SPREAD);            act_fcn%fcn => DataProcSpread
      CASE (INPUT_VAR_DO_BROADCAST);         act_fcn%fcn => DataProcBroadcast
      CASE (INPUT_VAR_DO_SCATTER);           act_fcn%fcn => DataProcScatter
      CASE (INPUT_VAR_DO_SUBSET);            act_fcn%fcn => DataProcSubSet
      CASE (INPUT_VAR_DO_INTERPOLATE);       act_fcn%fcn => DataProcInterpolateTimeLinear
      CASE (INPUT_VAR_DO_COMBINE_ADD);       act_fcn%fcn => DataProcAdd
      CASE (INPUT_VAR_DO_COMBINE_SUB);       act_fcn%fcn => DataProcSub
      CASE (INPUT_VAR_DO_COMBINE_SUBR);      act_fcn%fcn => DataProcSubr
      CASE (INPUT_VAR_DO_COMBINE_MUL);       act_fcn%fcn => DataProcMul
      CASE (INPUT_VAR_DO_COMBINE_DIV);       act_fcn%fcn => DataProcDiv
      CASE (INPUT_VAR_DO_COMBINE_DIVR);      act_fcn%fcn => DataProcDivr
      CASE (INPUT_VAR_DO_COMBINE_ADD_FLUX);  act_fcn%fcn => DataProcFlux
      CASE (INPUT_VAR_DO_COMBINE_NUDGE);     act_fcn%fcn => DataProcNudge
      CASE (INPUT_VAR_DO_COMBINE_NUDGE_FLD); act_fcn%fcn => DataProcNudgeFld
    END SELECT
#endif

    IF (.NOT. ASSOCIATED(var%saction)) THEN
      ALLOCATE(var%saction(100))
      var%saction(:) = 0
    ENDIF
    IF (SIZE(var%saction) < curr%acp + SIZE(act) - 1) THEN
      ALLOCATE(actions(SIZE(var%saction)+100))
      actions(1:curr%acp-1) = var%saction(1:curr%acp-1)
      actions(curr%acp:SIZE(actions)) = 0
      DEALLOCATE(var%saction)
      var%saction => actions
    ENDIF

    var%saction(curr%acp:curr%acp+SIZE(act)-1) = act
    curr%acp_pre = curr%acp
    curr%acp = curr%acp + SIZE(act)

  END SUBROUTINE InputAddDataAction

  SUBROUTINE AddUserAction(Level,var,Sz,fdim,ids,curr_proc,global)
  TYPE (input_var_list),  POINTER :: var
  TYPE (input_dim_list),  POINTER :: fdim
  TYPE (input_curr),INTENT(inout) :: curr_proc, global
  INTEGER, INTENT(in)             :: Level
  INTEGER, INTENT(inout)          :: Sz(:), ids(:)

#ifdef HAVE_F2003
  TYPE (input_data_list), POINTER :: curr
  TYPE (input_var_list),  POINTER :: tmp
  INTEGER :: bak

    IF (.NOT. ASSOCIATED(var%AddActionFcn)) RETURN
    bak = curr_proc%acp
    tmp => var ! For some reason at least nagfor do not accept calls to dummy argument procedure references
    CALL tmp%AddActionFcn(Level,var,Sz,fdim,curr_proc,ids)
    ! Update data buffer positions if action has been added
    IF (bak /= curr_proc%acp) THEN
      curr_proc%acp_pre = bak
      curr => var%buffers
      IF (ASSOCIATED(curr)) THEN
        DO WHILE (ASSOCIATED(curr%next))
          curr => curr%next
        ENDDO
      ENDIF
      curr_proc%buf => curr
      IF (ASSOCIATED(curr)) THEN
        curr => AllInputData
        DO WHILE (ASSOCIATED(curr))
          IF (ASSOCIATED(curr%dta,curr_proc%buf%dta)) THEN
            global%buf_pre => global%buf
            global%buf     => curr
            NULLIFY(curr)
          ELSE
            curr => curr%next
          ENDIF
        ENDDO
      ENDIF
    ENDIF
#endif

  END SUBROUTINE AddUserAction

! Perform necessary dimension permutations to join certain dimensions of a variable
  SUBROUTINE Permute2Collect(var,nDim,DimIds,Sz,fdims,dim_type,dim_order,lAddAction,curr_proc,global,mul,rep)
  TYPE (input_var_list),   POINTER :: var
  INTEGER,           INTENT(in   ) :: nDim
  INTEGER,           INTENT(inout) :: DimIds(7)
  INTEGER,           INTENT(inout) :: Sz(7)
  TYPE (input_dim_list),   POINTER :: fdims
  INTEGER,           INTENT(in   ) :: dim_type
  TYPE (input_dim_list),   POINTER :: dim_order
  LOGICAL,           INTENT(in   ) :: lAddAction
  TYPE (input_curr), INTENT(inout) :: curr_proc, global
  INTEGER,           INTENT(out  ) :: mul, rep

  TYPE (input_dim_list),   POINTER :: org, curr, prev, head, tmp_dim, tt
  INTEGER :: Order(7), SzOrg(7), tmpaction(16)
  INTEGER :: i, j, k, n, add, tmp, t2
  LOGICAL :: lCont, lReorder

    ! Preparations
    NULLIFY(org)
    DO i=1,nDim
      Order(i) = i
    ENDDO
    SzOrg(:) = Sz(:)

    IF (dim_type /= 0) THEN ! If dim_type "not given", the permutations should be done on all dims according to dim_order 
      ! Count expected number of dimensions to collect
      n = 0
      curr => dim_order
      DO WHILE (ASSOCIATED(curr))
        n = n + 1
        curr => curr%next
      ENDDO
      IF (COUNT(IAND(DimIds(1:nDim),INPUT_DIM_FND_MASK) == dim_type) /= n) CALL local_error('Permute2Collect', &
        'Inconsistent dimension mapping')

      ! --- Collect interesting dimensions
      ! Find first relevant dimension
      NULLIFY(org)
      head => fdims
      DO i=1,nDim
        IF (IAND(DimIds(i),INPUT_DIM_FND_MASK) == dim_type) EXIT
        org  => head       ! org => last dimension before the interesting ones
        head => head%next
      ENDDO
      curr => head%next
      prev => head
      mul = PRODUCT(Sz(1:i-1))
      ! Find any other relevant dimension
      add = 1
      j = i+1
      DO WHILE (j<=nDim)
        IF ((IAND(DimIds(j),INPUT_DIM_FND_MASK) == dim_type) .AND. (j-i /= add)) THEN
          ! If necessary: Swap dimensions in list together with Order, DimIds and Sz
          tmp_dim   => head%next
          head%next => curr
          prev%next => curr%next
          curr%next => tmp_dim
          head      => curr
          curr      => prev%next
          tmp = DimIds(j)
          t2  = Sz(j)
          DO k=j,i+add+1,-1
            Order (k) = Order (k-1)
            DimIds(k) = DimIds(k-1)
            Sz    (k) = Sz    (k-1)
          ENDDO
          Order (i+add) = j
          DimIds(i+add) = tmp
          Sz    (i+add) = t2
          add = add + 1
          j   = add + i
        ELSE ! This is an uninteresting dimension
          prev => curr
          curr => curr%next
          j = j + 1
        ENDIF
      ENDDO
      rep = PRODUCT(Sz(i+add+1:nDim))
    ELSE
      rep = 1
      mul = 1
    ENDIF

    ! --- Now secure that all interesting dimensions are in the right order
    lReorder = .FALSE.
    IF (ASSOCIATED(org)) THEN
      head => org%next
    ELSE
      head => fdims
      i = 1
    ENDIF
    ! head => First interesting dimension
    curr => dim_order             ! Correct dimension at this position in dim_order list
    DO WHILE (ASSOCIATED(curr))
      j = i
      NULLIFY(prev)
      tmp_dim => head
      lCont = ASSOCIATED(tmp_dim)
      IF (lCont) lCont = .NOT. ASSOCIATED(tmp_dim%dim_data,curr%dim_data)
      DO WHILE (lCont)
        j = j + 1
        prev    => tmp_dim
        tmp_dim => tmp_dim%next
        lCont = ASSOCIATED(tmp_dim)
        IF (lCont) lCont = .NOT. ASSOCIATED(tmp_dim%dim_data,curr%dim_data)
      ENDDO
      ! tmp_dim => Correct dimension to be placed at this position in fdims list
      IF (ASSOCIATED(tmp_dim)) THEN ! curr dimension found in dim_order list
        IF (.NOT. ASSOCIATED(head%dim_data,tmp_dim%dim_data)) THEN ! This is not the position where it was found
          lReorder = .TRUE.
          ! Swap Order, DimIds and Sz
          t2  = Sz(j)
          DO k=j,i+1,-1
            Order (k) = Order (k-1)
            Sz    (k) = Sz    (k-1)
          ENDDO
          Order (i) = j
          Sz    (i) = t2
        ELSE
          head => head%next
        ENDIF
      ELSE
        ! Dimensions missing in the "head" chain are considered to be at the right place
      ENDIF
      i = i + 1
      curr => curr%next
    ENDDO
    ! Permute the dimensions themselves
    IF (lReorder) THEN
      NULLIFY(prev)
      head => fdims
      i = 1
      DO WHILE (Order(i)==i)
        prev => head
        head => head%next
        i = i + 1
      ENDDO
      org => head ! First dimension out of order
      DO k=i,nDim
        tt => org
        DO j=i+1,Order(k)
          tt => tt%next
        ENDDO
        CALL InputDimCopy(tt,tmp_dim,ldata=.FALSE.)
        IF (ASSOCIATED(prev)) THEN
          prev%next => tmp_dim
        ELSE
          fdims => tmp_dim
        ENDIF
        prev => tmp_dim
      ENDDO
      DO WHILE (ASSOCIATED(org))
        prev => org%next
        DEALLOCATE(org)
        org => prev
      ENDDO
    ENDIF

    ! Add permutation if needed
    IF (lAddAction .AND. ANY(Order(2:nDim)-Order(1:nDim-1) /= 1)) THEN
      tmpaction(1     :2       ) = (/INPUT_VAR_DO_PERMUTE,nDim/)
      tmpaction(3     :2  +nDim) = SzOrg(1:nDim)
      tmpaction(3+nDim:2+2*nDim) = Order(1:nDim)
      CALL InputAddDataAction(var,tmpaction(1:2+2*nDim),curr_proc)
      CALL InputVarNewBuffer(var,Sz,curr_proc,global)
    ENDIF

  END SUBROUTINE Permute2Collect

! Add a mask to the current list
  SUBROUTINE AddVarMask(var,mask,curr)
  TYPE (input_var_list),    POINTER :: var
  TYPE (input_mask_list),  POINTER :: mask
  TYPE (input_curr), INTENT(inout) :: curr

    IF (ASSOCIATED(curr%mask)) THEN
      curr%mask%next => mask
    ELSE
      var%masks => mask
    ENDIF
    curr%mask => mask

  END SUBROUTINE AddVarMask

! Check file variable dimensions against those of the model and construct data processing chain
  RECURSIVE LOGICAL FUNCTION InputCheckVarDims(IFile)
  TYPE (input_file_list), POINTER :: IFile

  TYPE (input_var_list),     POINTER :: curr_var
  TYPE (input_eqdim_list),   POINTER :: eq_dim
  TYPE (input_dim_list),     POINTER :: curr_dim, tmp_dim, new_dim, fdim, prev_dim, src_dim
  TYPE (input_data_list),    POINTER :: buf
  TYPE (input_mask_list),    POINTER :: curr_mask, tmp_mask
  TYPE (input_file_list),    POINTER :: curr_file
  TYPE (input_fileobj_list), POINTER :: file_var, file_dim
#ifdef HAVE_F2003
  TYPE (input_fcn_list),     POINTER :: NewAction
#endif
  TYPE (input_curr)                  :: curr_proc, global
  REAL(dp),                  POINTER :: buf_dta(:,:,:,:,:,:,:)
  INTEGER,                   POINTER :: actions(:)
  INTEGER  :: nDim, i, j, k, m, n, st, en, ThisDim, DimIds(8), tmp_act(32), Sz(7), SzCheck, nIt, MaxLen, mul, rep
  INTEGER  :: chunk_lo(7), chunk_hi(7), dumar(1)
  REAL(dp) :: tst, mn, mx, add
  CHARACTER (len=64) :: name_dim, name_var
  CHARACTER (len=128) :: msg
  LOGICAL :: lCont, lThisPE, lLocal, lScatter, lTmp, lInit, lBroadCast, lActionDone

    ! Determine variables of this file
    curr_var => AllInputVars
    DO WHILE (ASSOCIATED(curr_var))
      lInit = .FALSE.
      DimIds(1:7) = INPUT_DIM_FND_NORMAL
      NULLIFY(global%buf,global%buf_pre,fdim,curr_mask,curr_proc%buf,curr_proc%buf_pre,curr_proc%mask)
      curr_proc%acp     = 1
      curr_proc%acp_pre = 1
      IF (.NOT. ASSOCIATED(curr_var%saction)) &
        CALL AddUserAction(INPUT_ACTION_OBTAIN,curr_var,curr_var%edims,curr_var%dims,DimIds,curr_proc,global)

      IF (ASSOCIATED(curr_var%parent_file,IFile) .AND. (curr_proc%acp==1) .AND. & ! If obtaining fcn is provided or
        IAND(curr_var%stat,INPUT_VAR_PROC_IN_FILE) == 0) THEN                     ! processing chain already made: Don't go on

        ! Build dimension chain for the variable as it is in the file
        CALL InputDimDone(curr_var%dims_at_read,.FALSE.,.FALSE.)
        NULLIFY(curr_var%dims_at_read,curr_dim,tmp_dim)
        DimIds(1:7) = curr_var%saction(1:7)
        IF (ASSOCIATED(curr_var%saction)) DEALLOCATE(curr_var%saction)
        i = 1
        DO WHILE ((i<=7) .AND. (DimIds(i) /= 0))
          tmp_dim => IFile%dims
          lCont = ASSOCIATED(tmp_dim)
          DO WHILE (lCont)
            IF (tmp_dim%fid == DimIds(i)) THEN
              lCont = .FALSE.
            ELSE
              tmp_dim => tmp_dim%next
              lCont = ASSOCIATED(tmp_dim)
            ENDIF
          ENDDO
          IF (.NOT. ASSOCIATED(tmp_dim)) CALL local_error('InputCheckVarDims','Internal error1 - dimension disappeared')
          IF (IAND(tmp_dim%dim_data%flags,INPUT_DIM_UNLIM) /= 0) THEN
            curr_var%stat = IOR(curr_var%stat,INPUT_VAR_HAS_UNLIM)
          ELSE
            CALL InputDimCopy(tmp_dim,new_dim,ldata=.FALSE.)
            IF (ASSOCIATED(curr_dim)) THEN
              curr_dim%next => new_dim
            ELSE
              curr_var%dims_at_read => new_dim
            ENDIF
            curr_dim => new_dim
            new_dim%flags    =  IAND(tmp_dim%flags,INPUT_ALLFLAGS - INPUT_DIM_TYPE) + INPUT_DIM_VAR + INPUT_DIM_READ
          ENDIF
          i = i + 1
        ENDDO

        ! If variable is not available, assume it has got the right shape
        IF (IAND(curr_var%stat,INPUT_VAR_EXIST) == 0) curr_var%dims_at_read => InputDimCopyList(curr_var%dims)

        ! Check if variable is a multi-variable and add appropriate dimension if necessary
        IF (IAND(curr_var%stat,INPUT_VAR_MULTI) /= 0) THEN
          ! Replace wildcards in variable names with names of matching variables in file
          st = 1
          DO WHILE (st > 0)
            lTmp = .FALSE.
            st = INDEX(curr_var%name_alt,'*')
            IF (st > 0) THEN
              lTmp = .TRUE.
            ELSE
              st = INDEX(curr_var%name_alt,'?')
            ENDIF
            IF (st > 0) CALL InputFileObjExpandVar(IFile%fvars,curr_var%name_alt,st,lTmp)
          ENDDO
          ! Count variables to form new dimension
          MaxLen = INDEX(curr_var%name_alt,'=>')
          IF (MaxLen == 0) MaxLen = LEN_TRIM(curr_var%name_alt) + 1
          st = 1
          n  = 0
          DO WHILE (st <= MaxLen)
            en = INDEX(curr_var%name_alt(st:MaxLen),',') + st - 1
            IF (en < st) en = MaxLen
            n = n + 1
            st = en + 1
          ENDDO
          ! Map extra dimension from multi-var to existing model dimension...
          st = INDEX(curr_var%name_alt,'=>')
          IF (st > 0) THEN
            name_dim = TRIM(curr_var%name_alt(st+2:))
            new_dim => InputDimGetRef(TRIM(name_dim))
            IF (.NOT. ASSOCIATED(new_dim)) CALL local_error('InputCheckVarDims',            &
              'Attempt to map extra dimension of multi-variable "'//TRIM(curr_var%name_var) &
              //'" to nonexisting dimension "'//TRIM(name_dim)//'"')
            IF (new_dim%dim_data%size_global /= n) THEN
              WRITE(message_text,*) 'Multi-variable "',TRIM(curr_var%name_var),'" has ',n,                 &
                ' members. Size mismatch when mapping to dimension "',TRIM(name_dim),'", which has size ', &
                new_dim%dim_data%size_global
              CALL RmDblSpc(message_text)
              CALL local_error('InputCheckVarDims',TRIM(message_text)//'.')
            ENDIF
          ELSE
          ! ... or find or create an internal extra dimension
            WRITE(name_dim,*) n
            st = 1
            DO WHILE (name_dim(st:st)==' ')
              st = st + 1
            ENDDO
            name_dim = 'internal_dim'//TRIM(name_dim(st:))
            new_dim => InputDimGetRef(TRIM(name_dim))
            IF (.NOT. ASSOCIATED(new_dim)) THEN
              CALL InputDimAdd(TRIM(name_dim),Extent=n,dim_ref=new_dim,intr=.true.)
              new_dim%flags = IOR(new_dim%flags,INPUT_DIM_EXTRA)
            ENDIF
          ENDIF
          ! Add new/extra dimension to variable dimensions
          CALL InputDimNewCtl(tmp_dim)
          tmp_dim%dim_data => new_dim%dim_data
          tmp_dim%src_dims => new_dim%src_dims
          tmp_dim%flags    =  IAND(new_dim%flags,INPUT_ALLFLAGS - INPUT_DIM_TYPE) + INPUT_DIM_VAR
          IF (ASSOCIATED(curr_dim)) THEN
            curr_dim%next => tmp_dim
          ELSE
            curr_var%dims_at_read => tmp_dim
          ENDIF
        ENDIF

        ! Create a working copy of the dimension list, which can step by step be modified to approach the expected model dims
        fdim => InputDimCopyList(curr_var%dims_at_read)

        ! Determine if variable should be scattered or broadcast (needed if any PE does not read its own data)
        lScatter = .FALSE.
        IF (n_pe>1) THEN
          DO i=1,n_pe
            IF (src_pe(i) /= i-1) lScatter = .TRUE.
          ENDDO
        ENDIF

        ! Other preparations
        CALL InputMasksDone(curr_var%masks,  .FALSE.)
        CALL InputDataDone (curr_var%buffers,.FALSE.)
        NULLIFY(curr_var%masks,curr_var%masks_intrp,curr_var%buffers,curr_var%buffers_intrp)
#ifdef HAVE_F2003
        CALL InputFcnDone(curr_var%DataFcns)
        NULLIFY(curr_var%DataFcns,curr_var%ProcInterpol)
#endif
        curr_var%acp_intrp = 0
        lLocal   = (n_pe <= 1)
        lThisPE  = liodata
        ! Find any non-model dimension and dimensions which are not expected to be associated with this variable
        DimIds(:) = -1   ! Use as flag
        msg = ''
        i = 1
        curr_dim => fdim
        DO WHILE (ASSOCIATED(curr_dim))
          DimIds(i) = InputDimGetRefEx(curr_dim,curr_var%dims,new_dim,ThisDim,fdims=curr_var%parent_file%fdims)
          IF (DimIds(i) == INPUT_DIM_FND_NOT) msg = TRIM(msg)//', '//TRIM(curr_dim%dim_data%name_dim)
          curr_dim => curr_dim%next
          i = i + 1
        ENDDO
        ! Errors and warnings from non-expected dimensions (sub-read, sum or average)
        i = COUNT(DimIds(:)==INPUT_DIM_FND_NOT)

        ! In this list INPUT_VAR_DIM_MIS_NORM = .._SUM + .._AVG - therefore not mentioned
        IF (((IAND(curr_var%flags,INPUT_VAR_DIM_MIS_VAR_MISS+INPUT_VAR_DIM_MIS_SUM+INPUT_VAR_DIM_MIS_AVG+ &
                                  INPUT_VAR_DIM_MIS_SUBSET+INPUT_VAR_DIM_MIS_RESOLVE) == 0) .AND.         &
            (i>0)) .OR. ((IAND(curr_var%flags,INPUT_VAR_DIM_MIS_SUBSET) /= 0) .AND. (i>1))) THEN
          CALL MessageDims(curr_var,fdim)
          CALL local_error('InputCheckVarDims',TRIM(curr_var%name_var)//': Too many unexpected dimensions: '// &
                                               TRIM(msg(3:)))
        ENDIF
        IF ((i>0) .AND. (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_WARN) /= 0) .AND. ldbg) THEN
          name_var = ''
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_VAR_MISS) == INPUT_VAR_DIM_MIS_VAR_MISS) &
            name_var = TRIM(name_var)//', by treating variable as missing'
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_RESOLVE)  == INPUT_VAR_DIM_MIS_RESOLVE ) &
            name_var = TRIM(name_var)//', user action'
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_SUBSET)   == INPUT_VAR_DIM_MIS_SUBSET  ) &
            name_var = TRIM(name_var)//', subsetting'
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_SUM)      == INPUT_VAR_DIM_MIS_SUM     ) &
            name_var = TRIM(name_var)//', cumulating'
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_AVG)      == INPUT_VAR_DIM_MIS_AVG     ) &
            name_var = TRIM(name_var)//', averaging'
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_NORM)     == INPUT_VAR_DIM_MIS_NORM    ) &
            name_var = TRIM(name_var)//', normalized averaging'
          CALL local_message(TRIM(curr_var%name_var),'Unexpected dimension(s) ('//TRIM(msg(3:))// &
                             ') found. Resolved by '//TRIM(name_var(3:))//'.')
        ENDIF
        IF ((i>0) .AND. (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_MASK) == INPUT_VAR_DIM_MIS_VAR_MISS)) THEN
          IF ( (IAND(curr_var%flags,INPUT_VAR_MISS_MASK) == INPUT_VAR_MISS_ERR) .OR. &
              ((IAND(curr_var%flags,INPUT_VAR_MISS_MASK) >  INPUT_VAR_MISS_WARN_IGNORE) .AND. (curr_var%nRead == 0))) THEN
            CALL MessageDims(curr_var,fdim)
            CALL local_error('InputCheckVarDims','Variable: '//TRIM(curr_var%name_var)//': Too many unexpected dimensions: ' &
                                                 //TRIM(msg(3:)))
          ENDIF
          buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_DEFAULT,0)
          IF (ASSOCIATED(buf_dta)) THEN
            ! Set variable to default value, detach variable from non-existing file and issue warning if requested
            curr_var%parent_file%nVar = curr_var%parent_file%nVar - 1
            curr_var%next_update = HUGE(i)
            NULLIFY(curr_var%parent_file,curr_var%group)
            IF (ASSOCIATED(curr_var%saction)) DEALLOCATE(curr_var%saction)
            IF (ASSOCIATED(curr_var%dta)) THEN
              IF (SIZE(buf_dta) == 1) THEN
                curr_var%dta(:,:,:,:,:,:,:) = buf_dta(1,1,1,1,1,1,1)
              ELSE
                IF (.NOT. ArrayCopy(buf_dta,curr_var%dta,(/SIZE(buf_dta)/),SHAPE(curr_var%dta))) &
                  CALL local_error('InputCheckVarDims','Spread error')
              ENDIF
            ENDIF
            IF (IAND(curr_var%flags,INPUT_VAR_MISS_MASK) >  INPUT_VAR_MISS_IGNORE) THEN
              IF (SIZE(buf_dta) == 1) THEN
                WRITE(message_text,'(2a,G10.5)') TRIM(curr_var%name_var), &
                  ': No method to map file dimensions to expected model dimensions. Using default value: ',buf_dta(1,1,1,1,1,1,1)
              ELSE
                WRITE(message_text,'(2a)') TRIM(curr_var%name_var), & 
                  ' No method to map file dimensions to expeced model dimensions. Using values from default vector.'
              ENDIF
              CALL local_message('WARNING',TRIM(message_text))
            ENDIF
            CYCLE
          ELSE
            CALL MessageDims(curr_var,fdim)
            CALL local_error('InputCheckVarDims',TRIM(curr_var%name_var)//': Too many unexpected dimensions ('// &
                                                 TRIM(message_text(3:))//') and no default value has been provided')
          ENDIF
        ENDIF

        ! Start action chain
        nDim = COUNT(DimIds(:) > -1) ! Current number of dimensions

        ! Get file size of variable and set tmp_act for reading
        tmp_act(1:3) = (/INPUT_VAR_DO_TOREAD,curr_var%fid,nDim/)
        j = 1
        curr_dim => fdim
        DO WHILE (ASSOCIATED(curr_dim))        
          IF (IAND(curr_var%flags,INPUT_VAR_GLOBAL) /= 0) THEN
            chunk_lo(j) = 1
            chunk_hi(j) = curr_dim%dim_data%size_global
          ELSE
            IF (((IAND(DimIds(j),INPUT_DIM_FND_SOURCE) /= 0) .OR. (IAND(DimIds(j),INPUT_DIM_FND_EQ_SRC) /= 0)) &
                .AND. (curr_dim%lo > 0)) THEN
              chunk_lo(j) = curr_dim%lo ! TODO: Missing something for reversed dimension?
              chunk_hi(j) = curr_dim%hi 
            ELSE
              chunk_lo(j) = curr_dim%dim_data%chunk_lo
              chunk_hi(j) = curr_dim%dim_data%chunk_hi
            ENDIF
          ENDIF
          IF ((IAND(curr_var%flags,INPUT_VAR_DIM_MIS_SUBSET) /= 0) .AND. (DimIds(j)==INPUT_DIM_FND_NOT)) THEN
            chunk_lo(j) = curr_var%subset_index
            chunk_hi(j) = chunk_lo(j)
          ENDIF
          IF (IAND(DimIds(j),INPUT_DIM_FND_REVERSE) == 1) THEN ! Dimension reversed
            chunk_lo(j) = curr_dim%dim_data%size_global - chunk_hi(j) + 1 ! Is this right
            chunk_hi(j) = chunk_hi(j) - i + 1
          ENDIF
          tmp_act(j     +3) = chunk_lo(j)
          tmp_act(j+nDim+3) = chunk_hi(j) - chunk_lo(j) + 1
          j = j + 1
          curr_dim => curr_dim%next
        ENDDO
        Sz(1:nDim)   = chunk_hi(1:nDim) - chunk_lo(1:nDim) + 1!tmp_act(nDim+4:2*nDim+3)
        Sz(nDim+1:7)       = 1
        chunk_lo(nDim+1:7) = 1
        chunk_hi(nDim+1:7) = 1

        ! Add multi option to the read-action if necessary
        j = 0
        IF (lThisPE .AND. (IAND(curr_var%stat,INPUT_VAR_MULTI) /= 0)) THEN
          ! Find position of last underscore in first name for eventual later name substitution
          m = 1
          k = INDEX(curr_var%name_alt,',')
          IF (k > 0) m = INDEX(curr_var%name_alt(1:k),'_',BACK=.TRUE.)
          ! Prepare multi-read request
          tmp_act(3+n:2*nDim+n+5) = tmp_act(1:2*nDim+3) ! Save read-request while making space for multi-request (n from above)
          tmp_act(1) = INPUT_VAR_DO_MULTI_READ
          tmp_act(2) = n
          j = n + 2
          ! Process names of individual variable names
          MaxLen = INDEX(curr_var%name_alt,'=>')
          IF (MaxLen == 0) MaxLen = LEN_TRIM(curr_var%name_alt) + 1
          st = 1
          curr_file => IFile
          DO k=1,n
            en = INDEX(curr_var%name_alt(st:),',') + st - 1
            IF (en < st) en = MaxLen
            IF (curr_var%name_alt(st:st) == '+') THEN
              name_var = curr_var%name_alt(1:m) // curr_var%name_alt(st+1:en-1)
            ELSE
              name_var = curr_var%name_alt(st:en-1)
            ENDIF
            file_var => InputFileObjGetRef(curr_file%fvars,name_var)
            IF (.NOT. ASSOCIATED(file_var)) &
              CALL local_error('InputCheckVarDims','Internal error - part of multi-variable disappeared')
            tmp_act(2+k) = file_var%id
            IF (IAND(curr_var%stat,INPUT_VAR_MULTI_SPREAD) /= 0) curr_file => curr_file%next
            st = en + 1
          ENDDO
        ENDIF

        ! Read on io-pes
        IF (lThisPE) THEN
          CALL InputAddDataAction(curr_var,tmp_act(1:3+2*nDim+j),curr_proc)
          CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global)
        ENDIF

        ! Throw away sub-read dimension
        IF ((i>0) .AND. (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_SUBSET) /= 0)) THEN
          j = 1
          DO WHILE (DimIds(j) /= INPUT_DIM_FND_NOT)
            j = j + 1
          ENDDO
          fdim => InputDimListElmRemove(fdim,j)
          DimIds(1:nDim-1) = PACK(DimIds(1:nDim),mask=DimIds(1:nDim) /= INPUT_DIM_FND_NOT)
          DimIds(nDim) = -1
          nDim = nDim - 1
          i = 0
        ENDIF

        CALL AddUserAction(INPUT_ACTION_AFTER_READ,curr_var,Sz,fdim,DimIds,curr_proc,global)

        ! Sum/average over any unexpected dimension(s)
        IF ((i>0) .AND. (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_NORM) /= 0)) THEN
          SELECT CASE (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_NORM))
            CASE (INPUT_VAR_DIM_MIS_SUM)
              tmp_act(1) = INPUT_VAR_DO_SUM
            CASE (INPUT_VAR_DIM_MIS_AVG)
              tmp_act(1) = INPUT_VAR_DO_AVG
            CASE (INPUT_VAR_DIM_MIS_NORM)
              tmp_act(1) = INPUT_VAR_DO_NORM
          END SELECT
          DO j=1,i
            k = 1
            DO WHILE (DimIds(k) /= INPUT_DIM_FND_NOT)
              k = k + 1
            ENDDO
            tmp_act(2) = k
            tmp_act(3) = Sz(k)
            buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_WEIGHTS,0)
            IF (ASSOCIATED(buf_dta)) THEN
              IF (SIZE(buf_dta) /= Sz(k)) THEN
                WRITE(message_text,*) 'Variable: ',TRIM(curr_var%name_var),': Length of weight vector (',SIZE(buf_dta) &
                                     ,') does not match size of extra dimension (',Sz(k),')'
                CALL local_error('InputCheckVarDims',TRIM(message_text))
              ENDIF
            ENDIF
            ! Size of new needed buffer
            Sz(k:6) = Sz(k+1:7) 
            Sz(7) = 1
            IF (lThisPE) THEN ! Still process only on io-pes
              CALL InputAddDataAction(curr_var,tmp_act(1:3),curr_proc)
              CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global)
            ENDIF
            fdim => InputDimListElmRemove(fdim,k)
            DimIds(k:7) = DimIds(k+1:8)
            nDim = nDim - 1
          ENDDO
        ENDIF

        CALL AddUserAction(INPUT_ACTION_AFTER_CONDENSE,curr_var,Sz,fdim,DimIds,curr_proc,global)

        ! Broadcast if necessary (i.e. variable contains no scattered dimensions and only a subset of pes are reading)
        lBroadCast = .FALSE.
        IF (lScatter .AND. (ALL(IAND(DimIds(1:nDim),INPUT_DIM_FND_DISTRIB+INPUT_DIM_FND_SCATTER)==0))) THEN
          tmp_act(1) = INPUT_VAR_DO_BROADCAST
          CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
          IF (.NOT. liodata) CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global) ! New buffer needed only on recv pes
          lThisPE = .TRUE. ! After a broadcast, all pes should do their work
          lBroadCast = .TRUE.
        ENDIF

        IF (lThisPE) THEN
          ! Add rescaling accounting for time
          IF (IAND(curr_var%flags,INPUT_VAR_DIV_MASK) /= 0) THEN
            tmp_act(1) = INPUT_VAR_DO_DIVTIME
            CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
          ENDIF

          ! Mark for rescaling if requested
          DO i=0,1
            buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_RESCALE,i*INPUT_DATA_NO_STEP)
            IF (ASSOCIATED(buf_dta)) THEN
              tmp_act(1) = INPUT_VAR_DO_RESCALE
              tmp_act(2) = i*INPUT_DATA_NO_STEP
              CALL InputAddDataAction(curr_var,tmp_act(1:2),curr_proc)
            ENDIF
          ENDDO
        ENDIF

        CALL AddUserAction(INPUT_ACTION_AFTER_RESCALE,curr_var,Sz,fdim,DimIds,curr_proc,global)

        nIt = 0 ! Limit number of iterations so that unresolved dim. problems results in  an error instead of an infinte loop
        DO WHILE (ANY(ABS(DimIds(:)) /= INPUT_DIM_FND_NORMAL) .AND. (nIt < 10)) ! "not defined"(-1) and "normal"(1) same abs value
          nIt = nIt + 1
          lActionDone = .FALSE.

          ! Swap data if necessary
          IF (lThisPE) THEN
            curr_dim => fdim
            DO i=1,ndim
              IF (IAND(DimIds(i),INPUT_DIM_FND_SWAP) /= 0) THEN
                file_dim => InputFileObjGetRef(curr_var%parent_file%fdims,curr_dim%dim_data%name_dim,curr_dim%file_dim_name)
                IF (.NOT. ASSOCIATED(file_dim)) CALL local_error('InputCheckVarDims','Internal error - file dim disappeared')
                tst = curr_dim%dim_data%origin
                IF (IAND(curr_dim%dim_data%flags,INPUT_DIM_REVERSE) /= 0) THEN ! TODO: untested case swap and reverse together!
                  mn  = file_dim%points(SIZE(file_dim%points))
                  mx  = file_dim%points(1)
                  add = mx - mn + mx - file_dim%points(2)
                ELSE
                  mx  = file_dim%points(SIZE(file_dim%points))
                  mn  = file_dim%points(1)
                  add = mx - mn + file_dim%points(2) - mn
                ENDIF
                DO WHILE (tst < mn)
                  tst = tst + add
                ENDDO
                DO WHILE (tst > mx)
                  tst = tst - add
                ENDDO
                dumar = MINLOC(ABS(file_dim%points-tst))
                IF (file_dim%points(dumar(1)) /= tst) CALL local_error('InputCheckVarDims','Values of dimension '// &
                                                                       TRIM(curr_dim%dim_data%name_dim)//' do not match')
                mul = PRODUCT(Sz(  1: i-1))
                rep = PRODUCT(Sz(i+1:ndim))
                tmp_act(1) = INPUT_VAR_DO_SWAP
                tmp_act(2) = mul*(dumar(1)-1)
                tmp_act(3) = mul*(Sz(i)-dumar(1)+1)
                tmp_act(4) = rep
                CALL InputAddDataAction(curr_var,tmp_act(1:4),curr_proc)
                CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global) ! Local recv-buffer allocated on every PE
                DimIds(i) = IAND(DimIds(i),INPUT_ALLFLAGS-INPUT_DIM_FND_SWAP)
              ENDIF
              curr_dim => curr_dim%next
            ENDDO

            ! Reverse dimensions which need it
            DO i=1,nDim
              IF (IAND(DimIds(i),INPUT_DIM_FND_REVERSE) /= 0) THEN
                tmp_act(1) = INPUT_VAR_DO_REVERSE
                tmp_act(2) = nDim
                tmp_act(3:2+nDim) = Sz(1:nDim)
                tmp_act(3+nDim) = i
                IF (lThisPE) THEN
                  CALL InputAddDataAction(curr_var,tmp_act(1:3+nDim),curr_proc)
                  CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global)
                ENDIF
                DimIds(i) = IAND(DimIds(i),INPUT_ALLFLAGS - INPUT_DIM_FND_REVERSE)
              ENDIF
            ENDDO

            CALL AddUserAction(INPUT_ACTION_AFTER_REVERSE,curr_var,Sz,fdim,DimIds,curr_proc,global)

          ENDIF ! lThisPE

          ! Distribute to non-io-pes if necessary (lScatter is false if each PE is an IO-PE
          IF (lScatter .AND. .NOT. lLocal .AND. ANY(IAND(DimIds(1:nDim),INPUT_DIM_FND_SCATTER) /= 0)) THEN 
            tmp_act(1) = INPUT_VAR_DO_SCATTER
            IF (n_iope < n_pe) THEN ! Scattering is only really needed if not every PE is doing it's own reading
              CALL InputMaskCreateScatter(fdim,curr_mask,k)
              CALL AddVarMask(curr_var,curr_mask,curr_proc)

              ! Forward to last mask and determine size of temporary send buffer if needed
              k = 0
              lCont = .TRUE. ! Here in the meaning "continuous data"
              curr_proc%mask => curr_mask
              IF (ASSOCIATED(curr_proc%mask)) THEN
                IF (IAND(curr_proc%mask%flags,INPUT_MASK_CONT) == 0) THEN
                  k = k + curr_proc%mask%nvalid
                  lCont = .FALSE.
                ENDIF
                DO WHILE (ASSOCIATED(curr_proc%mask%next))
                  curr_proc%mask => curr_proc%mask%next
                  IF (IAND(curr_proc%mask%flags,INPUT_MASK_CONT) == 0) THEN
                    k = k + curr_proc%mask%nvalid
                    lCont = .FALSE.
                  ENDIF
                ENDDO
              ENDIF
              IF (lCont) curr_var%stat = IOR(curr_var%stat,INPUT_VAR_CONT_SEND)
              ! This now should be done on every PE
              CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
              ! IO-server may need an extra buffer to reorder the data before sending. Could be a bit smaller, but it's easier so
              IF (lThisPE .AND. .NOT. lCont) CALL InputVarNewBuffer(curr_var,(/k,1,1,1,1,1,1/),curr_proc,global)
              ! Adjust working size and flags
            ENDIF
            curr_dim => fdim
            DO i=1,ndim
              IF (IAND(DimIds(i),INPUT_DIM_FND_SCATTER) /= 0) Sz(i) = curr_dim%dim_data%size_local
              DimIds(i) = IAND(DimIds(i),INPUT_ALLFLAGS - INPUT_DIM_FND_SCATTER - INPUT_DIM_FND_DISTRIB)
              curr_dim => curr_dim%next
            ENDDO
            IF (n_iope < n_pe) CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global) ! Local recv-buffer allocated on every PE
            lThisPE = .TRUE. ! From now on, everything should be performed on every PE
            lLocal  = .TRUE.
          ENDIF

          CALL AddUserAction(INPUT_ACTION_AFTER_BROADCAST,curr_var,Sz,fdim,DimIds,curr_proc,global)

          ! Model requires a dimension which is a pack of multiple non-distributed dimensions in the file?
          k = COUNT(IAND(DimIds(1:nDim),INPUT_DIM_FND_MASK)==INPUT_DIM_FND_PACKED)
          IF (k>0) THEN
            ! Find the correct packed dimension in the global list
            src_dim => AllInputDims
            m = -1
            DO WHILE (m /= k)
              lCont = ASSOCIATED(src_dim)
              DO WHILE (lCont)
                IF (ASSOCIATED(src_dim%src_dims)) THEN
                  lCont = .FALSE.
                ELSE
                  src_dim => src_dim%next
                  lCont = ASSOCIATED(src_dim)
                ENDIF
              ENDDO
              IF (.NOT. ASSOCIATED(src_dim)) CALL local_error('InputCheckVarDims','Internal error2 - dimension disappeared')
              m = 0
              ! Count number of source dimensions
              tmp_dim => src_dim%src_dims
              DO WHILE (ASSOCIATED(tmp_dim))
                tmp_dim => tmp_dim%next
                m = m + 1
              ENDDO
              IF (m == k) THEN ! This may be the right dimension one - test the individual dimensions
                i = 1
                curr_dim => fdim
                DO WHILE (ASSOCIATED(curr_dim))
                  IF (IAND(DimIds(i),INPUT_DIM_FND_MASK) == INPUT_DIM_FND_PACKED) THEN
                    tmp_dim => InputDimGetRef(TRIM(curr_dim%dim_data%name_dim),search=src_dim%src_dims)
                    IF (.NOT. ASSOCIATED(tmp_dim)) THEN
                      m = -1 ! Force iteration for another packed dimension in the global list
                      src_dim => src_dim%next
                      NULLIFY(curr_dim)
                    ELSE
                      curr_dim => curr_dim%next
                      i = i + 1
                    ENDIF
                  ELSE
                    curr_dim => curr_dim%next
                  ENDIF
                ENDDO
              ELSE
                src_dim => src_dim%next
              ENDIF
            ENDDO

            ! Collect all source dimensions of the packed dimension
            CALL Permute2Collect(curr_var,nDim,DimIds,Sz,fdim,INPUT_DIM_FND_PACKED,src_dim%src_dims, &
                                 lThisPE,curr_proc,global,mul,rep)

            ! Now prepare packing
            tmp_act(1:3) = (/INPUT_VAR_DO_PACK,mul,rep/)
            NULLIFY(prev_dim)
            curr_dim => fdim
            j = 1
            DO WHILE ((IAND(DimIds(j),INPUT_DIM_FND_MASK) /= INPUT_DIM_FND_PACKED) .AND. (DimIds(j) > INPUT_DIM_FND_NOT))
              prev_dim => curr_dim
              curr_dim => curr_dim%next
              j = j + 1
            ENDDO
            i = InputDimGetRefEx(curr_dim,search=curr_var%dims,ref=tmp_dim,fdims=curr_var%parent_file%fdims)
            ! Associate mask
            buf_dta => InputDataGetRef((/1,1,1,1,1,1,1/),ref=buf) ! Get a size 1 buffer to transfere a single scalar
            IF (lThisPE) THEN
              tmp_mask => tmp_dim%dim_data%mask
              lCont = ASSOCIATED(tmp_mask)
              m = PRODUCT(Sz(j:nDim))/rep
              IF (lLocal) THEN
              n = tmp_dim%dim_data%size_local
              ELSE
                n = tmp_dim%dim_data%chunk_hi-tmp_dim%dim_data%chunk_lo+1
              ENDIF
              DO WHILE (lCont)
                IF ((tmp_mask%nPts == m) .AND. (tmp_mask%nvalid == n)) THEN
                  lCont = .FALSE.
                ELSE
                  tmp_mask => tmp_mask%next
                  lCont = ASSOCIATED(tmp_mask)
                ENDIF
              ENDDO
              IF (.NOT. ASSOCIATED(tmp_mask)) CALL local_error('InputCheckVarDims','Variable: '// &
                                                TRIM(curr_var%name_var)//': No mask consistent with requested packing was found')
              NULLIFY(curr_mask)
              CALL InputMaskCopy(curr_mask,tmp_mask)
              CALL AddVarMask(curr_var,curr_mask,curr_proc)
              IF (lLocal) THEN
                n = curr_mask%nValid
              ELSE
                buf_dta(1,1,1,1,1,1,1) = REAL(curr_mask%nValid,dp)
              ENDIF
            ENDIF
            IF (.NOT. lLocal) THEN
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
              CALL InputVarBroadcast(buf)        
#endif
              n = INT(buf_dta(1,1,1,1,1,1,1))
            ENDIF
            ! Replace dimensions in dimension list
            k = j
            DO WHILE (IAND(DimIds(k),INPUT_DIM_FND_MASK) == INPUT_DIM_FND_PACKED)
              prev_dim => curr_dim
              curr_dim => curr_dim%next
              DEALLOCATE(prev_dim)
              k = k + 1
            ENDDO
            NULLIFY(new_dim)
            CALL InputDimCopy(tmp_dim,new_dim,ldata=.FALSE.)
            new_dim%next => curr_dim
            IF (ASSOCIATED(prev_dim)) THEN
              prev_dim%next => new_dim
            ELSE
              fdim => new_dim
            ENDIF
            new_dim%flags = IAND(tmp_dim%flags,INPUT_ALLFLAGS - INPUT_DIM_TYPE) + INPUT_DIM_VAR
            ! Set flags of new packed dimension and update array size
            DimIds(j) = InputDimGetRefEx(new_dim,curr_var%dims,fdims=curr_var%parent_file%fdims)
            DimIds(j+1:8-k+j) = DimIds(k:7)
            DimIds(9-k+j:8) = -1
            IF (lLocal) DimIds(j) = IAND(DimIds(j),INPUT_ALLFLAGS - INPUT_DIM_FND_DISTRIB - INPUT_DIM_FND_SCATTER)
            Sz(j) = n
            Sz(j+1:8-k+j) = Sz(k:7)
            Sz(9-k+j:7) = 1
            IF (lThisPE) THEN
              CALL InputAddDataAction(curr_var,tmp_act(1:3),curr_proc)
              CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global)
            ENDIF
            nDim = nDim - k + j + 1
            lActionDone = .TRUE.
          ENDIF

          ! Does the variable contain a packed dimension which is requested unpacked in the model?
          IF (.NOT. lActionDone .AND. (ANY(IAND(DimIds(1:nDim),INPUT_DIM_FND_MASK) == INPUT_DIM_FND_SOURCE))) THEN
            ! Prepare action
            tmp_act(1) = INPUT_VAR_DO_UNPACK
            j = 1
            DO WHILE (IAND(DimIds(j),INPUT_DIM_FND_MASK) /= INPUT_DIM_FND_SOURCE)
              j = j + 1
            ENDDO
            tmp_act(2) = PRODUCT(Sz(1:j-1))
            tmp_act(3) = PRODUCT(Sz(j+1:nDim))
            NULLIFY(prev_dim)
            curr_dim => fdim
            DO i=2,j
              prev_dim => curr_dim
              curr_dim => curr_dim%next
            ENDDO
            i = 0
            k = 1
            tmp_dim => curr_dim%src_dims
            DO WHILE (ASSOCIATED(tmp_dim))
              i = i + 1
              IF (lLocal) THEN
                k = k * tmp_dim%dim_data%size_local
              ELSE
                k = k * (tmp_dim%dim_data%chunk_hi - tmp_dim%dim_data%chunk_lo + 1)
              ENDIF
              tmp_dim => tmp_dim%next
            ENDDO
            ! Add mask
            tmp_mask => curr_dim%dim_data%mask
            IF ((curr_dim%lo > 0) .AND. ASSOCIATED(tmp_mask%next)) tmp_mask => tmp_mask%next
            NULLIFY(curr_mask)
            CALL InputMaskCopy(curr_mask,tmp_mask)
            CALL AddVarMask(curr_var,curr_mask,curr_proc)
            IF (lLocal) THEN
              SzCheck = curr_dim%dim_data%size_local
            ELSE
              SzCheck = chunk_hi(j) - chunk_lo(j) + 1
            ENDIF
            IF ((curr_mask%nValid /= SzCheck) .OR. (curr_mask%nPts /= k)) &
              CALL local_error('InputCheckVarDims','Variable: '//TRIM(curr_var%name_var)//': Dimension mask size mismatch')
            ! Replace dimensions in list
            DimIds(j+i:8) = DimIds(j+1:9-i)
            DimIds(j:j+i-1) = INPUT_DIM_FND_NORMAL
            Sz(j+i-1:7) = Sz(j:8-i)
            tmp_dim => curr_dim%src_dims
            DO k=1,i
              NULLIFY(new_dim)
              CALL InputDimCopy(tmp_dim,new_dim,ldata=.FALSE.)
              IF (ASSOCIATED(prev_dim)) THEN
                prev_dim%next => new_dim
              ELSE
                fdim => new_dim
              ENDIF
              prev_dim => new_dim
              new_dim%flags = IAND(tmp_dim%flags,INPUT_ALLFLAGS - INPUT_DIM_TYPE) + INPUT_DIM_VAR
              DimIds(j+k-1) = InputDimGetRefEx(new_dim,curr_var%dims,fdims=curr_var%parent_file%fdims)
              IF (lLocal) DimIds(j+k-1) = IAND(DimIds(j+k-1),INPUT_ALLFLAGS - INPUT_DIM_FND_DISTRIB - INPUT_DIM_FND_SCATTER)
              IF (lScatter) THEN
                Sz  (j+k-1) = new_dim%dim_data%chunk_hi - new_dim%dim_data%chunk_lo + 1
              ELSE
                Sz  (j+k-1) = new_dim%dim_data%local_hi - new_dim%dim_data%local_lo + 1
              ENDIF
              tmp_dim => tmp_dim%next
            ENDDO
            new_dim%next => curr_dim%next ! Link following dimensions
            DEALLOCATE(curr_dim)
            IF (lThisPE) THEN
              CALL InputAddDataAction(curr_var,tmp_act(1:3),curr_proc)
              CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global)
            ENDIF
            nDim = nDim + i - 1
            lActionDone = .TRUE.
          ENDIF

          ! Does the variable contain equivalent dimensions which is requested as their equivalents in the model?
          IF (.NOT. lActionDone .AND. (ANY(IAND(DimIds(1:nDim),INPUT_DIM_FND_MASK) == INPUT_DIM_FND_EQ_SRC))) THEN

            ! First perform eventual necessary permutation to bring source dims of equiv together in the right order
            eq_dim => AllInputEqDims
            DO WHILE (ASSOCIATED(eq_dim))
              lCont = .FALSE.
              tmp_dim => eq_dim%src_dims
              DO WHILE (ASSOCIATED(tmp_dim))
                IF (InputDimGetRefEx(tmp_dim,fdim,it=2) == INPUT_DIM_FND_NOT) THEN
                  NULLIFY(tmp_dim)
                  lCont = .TRUE.
                ELSE
                  tmp_dim => tmp_dim%next
                ENDIF
              ENDDO
              IF (.NOT. lCont) EXIT
              eq_dim => eq_dim%next
            ENDDO
            IF (.NOT. ASSOCIATED(eq_dim)) CALL local_error('InputCheckVarDims','Equivalent dimension definition disappeared')
            CALL Permute2Collect(curr_var,nDim,DimIds,Sz,fdim,INPUT_DIM_FND_EQ_SRC,eq_dim%src_dims, &
                                 lThisPE,curr_proc,global,mul,rep)
            ! Now replace the source dimension(s) with the destination one(s)
            i = 1
            j = 0
            k = 0
            NULLIFY(prev_dim)
            curr_dim => fdim
            DO WHILE (.NOT. ASSOCIATED(curr_dim%dim_data,eq_dim%src_dims%dim_data))
              i = i + 1
              prev_dim => curr_dim
              curr_dim => curr_dim%next
            ENDDO
            tmp_dim => curr_dim%next
            ! Copy the destination dimensions
            src_dim => eq_dim%dst_dims
            DO WHILE (ASSOCIATED(src_dim))
              j = j + 1
              CALL InputDimCopy(src_dim,new_dim,ldata=.FALSE.)
              new_dim%flags = IAND(new_dim%flags,INPUT_ALLFLAGS-INPUT_DIM_TYPE) + INPUT_DIM_VAR + INPUT_DIM_USEASGLOBAL
              new_dim%lo    = 1 ! Equivalent dimensions presently cannot use chunked reading!
              new_dim%hi    = new_dim%dim_data%size_global
              IF (ASSOCIATED(prev_dim)) THEN
                prev_dim%next => new_dim
              ELSE
                fdim => new_dim
              ENDIF
              prev_dim => new_dim
              src_dim  => src_dim%next
            ENDDO
            ! Delete the source dimensions
            src_dim => eq_dim%src_dims
            DO WHILE (ASSOCIATED(src_dim))
              k = k + 1
              new_dim => curr_dim%next
              DEALLOCATE(curr_dim)
              curr_dim => new_dim
              src_dim => src_dim%next
            ENDDO
            prev_dim%next => curr_dim ! Link eventually following dimensions

            ! Rebuild dimension types and sizes (easier to do a complete rebuild than to push data in a priori unknown direction)
            DimIds(1:7) = -1
            Sz(1:7)     =  1
            src_dim => fdim
            nDim = 0
            DO WHILE (ASSOCIATED(src_dim))
              nDim = nDim + 1
              DimIds(nDim) = InputDimGetRefEx(src_dim,curr_var%dims,fdims=curr_var%parent_file%fdims)
              IF (lLocal) DimIds(nDim) = IAND(DimIds(nDim),INPUT_ALLFLAGS - INPUT_DIM_FND_DISTRIB - INPUT_DIM_FND_SCATTER)
              IF (.NOT. lLocal .OR. (src_dim%dim_data%local_lo < 0)) THEN
                Sz(nDim) = src_dim%dim_data%size_global ! chunk_hi-src_dim%dim_data%chunk_lo + 1 !Chunked reading of eq. dim not impl
              ELSE
                Sz(nDim) = src_dim%dim_data%local_hi - src_dim%dim_data%local_lo + 1
              ENDIF
              src_dim => src_dim%next
            ENDDO

            ! Prepare the equivalence operation
            IF (lThisPE) THEN
              ! Add mask
              NULLIFY(curr_mask)
              CALL InputMaskCopy(curr_mask,eq_dim%mask)
              CALL AddVarMask(curr_var,curr_mask,curr_proc)
              j = curr_proc%acp_pre
              tmp_act(1:3) = (/INPUT_VAR_DO_EQUIV,mul,rep/)
              CALL InputAddDataAction(curr_var,tmp_act(1:3),curr_proc)
              ! If possible, remove last buffer to save data copying (time and space)
              IF ((rep==1) .OR. (curr_mask%mask(SIZE(curr_mask%mask))==0)) THEN ! No unused internal space
                dumar = MINLOC(ABS(INPUT_VAR_1D_SET - ABS(curr_var%saction(j))))
                IF ((INPUT_VAR_1D_SET(dumar(1)) == ABS(curr_var%saction(j))) .AND. &
                    (.NOT. liodata .OR. (curr_var%saction(j) /= INPUT_VAR_DO_BROADCAST))) THEN
                  buf_dta => curr_proc%buf%dta
                  CALL InputDataUnuse(curr_var%buffers,curr_proc%buf)     ! Remove buffer from variable's list
                  curr_proc%buf => curr_proc%buf_pre                      ! Attach next buffer to the one before
                  global%buf    => global%buf_pre                         ! Attach next global buffer to the one before
                  CALL InputDataUnuse(AllInputData,dta=buf_dta)           ! Remove buffer from global list
                  curr_var%saction(curr_proc%acp-3) = -INPUT_VAR_DO_EQUIV ! Mark that no data copy should actually take place
                ENDIF
              ENDIF
              CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global)
            ENDIF
          ENDIF

        ENDDO ! Any dimension not "normal"

        IF (nIt == 5) THEN
          lCont = .TRUE.
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_VAR_MISS) /= 0) THEN
            buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_DEFAULT,0)
            IF (ASSOCIATED(buf_dta)   .AND. ((IAND(curr_var%flags,INPUT_VAR_MISS_MASK) < INPUT_VAR_MISS_ERR) .OR. &
                ((curr_var%nRead > 0) .AND.  (IAND(curr_var%flags,INPUT_VAR_MISS_MASK) > INPUT_VAR_MISS_ERR)))) THEN
              lCont = .FALSE.
              ! Set variable to default value, detach variable from non-existing file and issue warning if requested
              curr_var%next_update = HUGE(i)
              curr_var%parent_file%nVar = curr_var%parent_file%nVar - 1
              NULLIFY(curr_var%parent_file,curr_var%group)
              DEALLOCATE(curr_var%saction)
              IF (ASSOCIATED(curr_var%dta)) THEN
                IF (SIZE(buf_dta) == 1) THEN
                  curr_var%dta(:,:,:,:,:,:,:) = buf_dta(1,1,1,1,1,1,1)
                ELSE
                  IF (.NOT. ArrayCopy(buf_dta,curr_var%dta,(/SIZE(buf_dta)/),SHAPE(curr_var%dta))) &
                    CALL local_error('InputCheckVarDims','Spread error')
                ENDIF
                IF (IAND(curr_var%flags,INPUT_VAR_MISS_MASK) >  INPUT_VAR_MISS_IGNORE) THEN
                  IF (SIZE(buf_dta) == 1) THEN 
                    WRITE(message_text,'(3a,G10.5)') 'No method to maps file dimensions of ',TRIM(curr_var%name_var), &
                      ' to the expected model dimensions. Using default value: ',buf_dta(1,1,1,1,1,1,1)
                  ELSE
                    WRITE(message_text,'(3a)') 'No method to maps file dimensions of ',TRIM(curr_var%name_var), &
                      ' to the expected model dimensions. Using values from default-value vector.'
                  ENDIF
                  CALL local_message('WARNING',TRIM(message_text))
                ENDIF
              ENDIF
              CYCLE
            ENDIF
          ENDIF
          IF (lCont) THEN
            CALL MessageDims(curr_var,fdim,'Dimensions%after%maximum%processing%of%dimensions%found%in%file')
            CALL local_error('InputCheckVardims','No method to correct dimensions of '//TRIM(curr_var%name_var))
          ENDIF
        ENDIF

        ! Broadcast if necessary (i.e. variable contains no scattered dimensions and only a subset of pes are reading)
        IF ((n_pe > 1) .AND. (n_iope < n_pe) .AND. .NOT. lLocal .AND. .NOT. lBroadcast) THEN
          tmp_act(1) = INPUT_VAR_DO_BROADCAST
          CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
          IF (.NOT. liodata) CALL InputVarNewBuffer(curr_var,Sz,curr_proc,global) ! New buffer needed only on recv pes
          lThisPE = .TRUE. ! After a broadcast, all pes should do their work
        ENDIF

        CALL AddUserAction(INPUT_ACTION_AFTER_PACK,curr_var,Sz,fdim,DimIds,curr_proc,global)

        ! Put dimensions in the same order as expected by the model - still missing dimensions may be in between
        CALL Permute2Collect(curr_var,nDim,DimIds,Sz,fdim,0,curr_var%dims,lThisPE,curr_proc,global,mul,rep)

        ! Test for missing dimensions
        i = 24
        tmp_dim => curr_var%dims
        DO WHILE (ASSOCIATED(tmp_dim))
          tmp_act(i) = InputDimGetRefEx(tmp_dim,fdim)
          i = i + 1
          tmp_dim => tmp_dim%next
        ENDDO
        j = COUNT(tmp_act(24:i-1)==INPUT_DIM_FND_NOT) ! # undefined dimensions
        IF (j > 0) THEN ! Model expects dimensions not in the file variable
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_SPREAD+INPUT_VAR_DIM_MIS_SUBSET) == 0) THEN
            IF ((IAND(curr_var%flags,INPUT_VAR_DIM_MIS_WARN) /= 0) .AND. lDbg) &
              CALL MessageDims(curr_var,fdim,'Dimensions%after%maximum%processing%of%dimensions%found%in%file')
            message_text = ''
            tmp_dim => curr_var%dims
            DO j=24,i-1
              IF (tmp_act(j)==INPUT_DIM_FND_NOT) message_text = TRIM(message_text)//','//TRIM(tmp_dim%dim_data%name_dim)
              tmp_dim => tmp_dim%next
            ENDDO
            CALL local_error('InputCheckVarDims','Variable: '//TRIM(curr_var%name_var)// &
                             ': Expected dimension(s) missing in file variable: '//TRIM(message_text(2:)))
          ENDIF
          IF ((j > 1) .AND. (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_SUBSET) /= 0)) THEN
            CALL MessageDims(curr_var,fdim,'Dimensions%after%maximum%processing%of%dimensions%found%in%file')
            CALL local_error('InputCheckVarDims',TRIM(curr_var%name_var)//': Can only subset one dimension')
          ENDIF
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_WARN) /= 0) THEN
            IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_SUBSET) /= 0) THEN
              WRITE(message_text,*) 'setting only index: ',curr_var%subset_index
            ELSE
              message_text = 'spreading data across dimension(s)'
            ENDIF
            IF (ldbg .AND. (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_WARN) /= 0)) THEN
              msg = ''
              curr_dim => curr_var%dims
              DO k=24,i-1
                IF (tmp_act(k) == INPUT_DIM_FND_NOT) &
                  msg = TRIM(msg)//','//TRIM(curr_dim%dim_data%name_dim)
                curr_dim => curr_dim%next
              ENDDO
              CALL local_message('WARNING',TRIM(curr_var%name_var)//': Requested dimension(s) ('//TRIM(msg(2:))// &
                                           ') not found in file. Resolved by '//TRIM(message_text))
            ENDIF
          ENDIF
          k = 1
          DO WHILE (tmp_act(k+23) /= INPUT_DIM_FND_NOT)
            k = k + 1
          ENDDO
          lTmp = .TRUE.
          IF (IAND(curr_var%flags,INPUT_VAR_DIM_MIS_SUBSET) /= 0) THEN ! Subset variable
            tmp_act(1:3) = (/INPUT_VAR_DO_SUBSET,k,curr_var%subset_index/)
            CALL InputAddDataAction(curr_var,tmp_act(1:3),curr_proc)
          ELSE ! Spread variable
            tmp_act(1:2     ) = (/INPUT_VAR_DO_SPREAD,nDim/)
            tmp_act(3:2+nDim) = Sz(1:nDim)
            tmp_act(  3+nDim) = j
            lTmp = .FALSE.
            i = 0
            DO WHILE (k<=nDim+j)
              tmp_act(4  +i+nDim) = k
              tmp_act(4+j+i+nDim) = curr_var%edims(k)
              IF (curr_var%edims(k) > 1) lTmp = .TRUE.
              k = k + 1
              DO WHILE ((k<=nDim+j) .AND. (tmp_act(k+23) /= INPUT_DIM_FND_NOT))
                k = k + 1
              ENDDO
              i = i + 1
            ENDDO
            IF (lTmp) CALL InputAddDataAction(curr_var,tmp_act(1:3+2*j+nDim),curr_proc)
          ENDIF
          IF (lTmp) CALL InputVarNewBuffer(curr_var,curr_var%edims,curr_proc,global) ! This the last action in the action chain
        ELSE
          lCont = .FALSE.
        ENDIF

        ! If missing dimension issues were resolved the previous reorder may have screwed up - thus redo
        IF (lCont) CALL Permute2Collect(curr_var,nDim,DimIds,Sz,fdim,0,curr_var%dims,lThisPE,curr_proc,global,mul,rep)

        CALL AddUserAction(INPUT_ACTION_AFTER_DIM_CORRECT,curr_var,Sz,fdim,DimIds,curr_proc,global)

        ! Replace fill values from file if possible
        IF (IAND(curr_var%flags,INPUT_VAR_BAD_VAL_MASK) /= INPUT_VAR_BAD_VAL_IGNORE) THEN
          lCont = .FALSE.
          tst = UNLIKELY_RVAL
          buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_FILL,INPUT_DATA_NO_STEP)
          IF (.NOT. ASSOCIATED(buf_dta)) &
            buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_FILL,0)
          IF (ASSOCIATED(buf_dta)) THEN
            IF (buf_dta(2,1,1,1,1,1,1) /= UNLIKELY_RVAL) THEN
              tmp_act(1) = INPUT_VAR_DO_FILL_BAD
              CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
              lCont = .TRUE.
            ENDIF
            tst = buf_dta(1,1,1,1,1,1,1)
          ENDIF
          buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_RANGE,INPUT_DATA_NO_STEP)
          IF (.NOT. ASSOCIATED(buf_dta)) &
            buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_RANGE,0)
          IF (ASSOCIATED(buf_dta)) THEN
            IF (buf_dta(2,1,1,1,1,1,1) /= UNLIKELY_RVAL) THEN
              tmp_act(1) = INPUT_VAR_DO_FILL_RANGE
              CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
              lCont = .TRUE.
            ENDIF
            IF (tst == UNLIKELY_RVAL) tst = buf_dta(1,1,1,1,1,1,1)
          ENDIF
          IF ((tst == UNLIKELY_RVAL) .AND. (IAND(curr_var%flags,INPUT_VAR_BAD_VAL_REPLACE) /= 0)) &
            curr_var%flags = IAND(curr_var%flags,INPUT_ALLFLAGS - INPUT_VAR_BAD_VAL_REPLACE)
          IF (IAND(curr_var%flags,INPUT_VAR_BAD_VAL_MASK) == INPUT_VAR_BAD_VAL_DEFAULT) THEN
            IF (lCont) THEN
              IF (tst /= UNLIKELY_RVAL) THEN
                curr_var%flags = IOR(IAND(curr_var%flags,INPUT_ALLFLAGS - INPUT_VAR_BAD_VAL_MASK), INPUT_VAR_BAD_VAL_WARN_REPLACE)
              ELSE
                curr_var%flags = IOR(IAND(curr_var%flags,INPUT_ALLFLAGS - INPUT_VAR_BAD_VAL_MASK), INPUT_VAR_BAD_VAL_ERR)
              ENDIF
            ELSE
              curr_var%flags = IOR(IAND(curr_var%flags,INPUT_ALLFLAGS - INPUT_VAR_BAD_VAL_MASK), INPUT_VAR_BAD_VAL_ERR)
            ENDIF
          ENDIF
        ENDIF

        ! Add interpolation step if required
        lCont = curr_var%dt_update /= 0
        IF (lCont) THEN
          IF (INT(curr_var%dt_update,i8)*INT(IFile%dt_file,i8) > 0) THEN ! 32 bit integers may overflow
            lCont = ABS(curr_var%dt_update) < ABS(IFile%dt_file)
          ELSE
            lCont = .TRUE.
          ENDIF
        ENDIF
        IF (lCont .OR. (curr_var%nIntrp > 1)) THEN
          NULLIFY(curr_proc%buf%dta) ! Empty data in buffer means: Use actual interpolation buffer
          curr_var%buffers_intrp => curr_proc%buf
          curr_var%masks_intrp   => curr_mask
          curr_var%acp_intrp     =  curr_proc%acp
          CALL InputVarNewBuffer(curr_var,curr_var%edims,curr_proc,global)
          tmp_act(1) = INPUT_VAR_DO_INTERPOLATE
#ifdef HAVE_F2003
          CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc,NewAction)
          IF (ASSOCIATED(curr_var%TimeIntrpFcn)) NewAction%fcn => curr_var%TimeIntrpFcn
          curr_var%ProcInterpol => NewAction
#else
          CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
#endif
        ENDIF

        ! Add extra buffer for subsetting if necessary - TODO: Cannot be right to put it here!!!
        IF (tmp_act(1) == INPUT_VAR_DO_SUBSET) CALL InputVarNewBuffer(curr_var,curr_var%edims,curr_proc,global)

        CALL AddUserAction(INPUT_ACTION_AFTER_INTERPOLATE,curr_var,Sz,fdim,DimIds,curr_proc,global)

        curr_var%stat = IOR(curr_var%stat,INPUT_VAR_INITIALIZED)
        lInit = .TRUE.
      ELSE
        ! If this is a non-initialized non-file variable, make buffers and primitive time settings for it
        IF (.NOT. ASSOCIATED(curr_var%parent_file) .AND. IAND(curr_var%stat,INPUT_VAR_INITIALIZED) == 0) THEN
          ! Make sure that update messages are correctly created
          tmp_act(1) = INPUT_VAR_DO_MESSAGE
          CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
          CALL InputVarNewBuffer(curr_var,curr_var%edims,curr_proc,global)
          IF (curr_var%time(1) == -HUGE(i)) THEN
            curr_var%time  = model_time
            CALL CalendarTimeSub(curr_var%time,INT(curr_var%dt_update,i8))
          ENDIF
          IF (curr_var%nIntrp < 1) THEN
            curr_var%nIntrp  = 1
            curr_var%nAtOnce = 1
            curr_var%nFuture = 1
          ENDIF
          ! Add actions which can also be done for manual data providing
          lTmp = .FALSE.
          DO i=0,1
            buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_RESCALE,i*INPUT_DATA_NO_STEP)
            lTmp = lTmp .OR. ASSOCIATED(buf_dta)
          ENDDO
          ! For data to be manually pushed, add a NoAction to the processing chain to prevent the normal processing,
          !   which would operate on still unavailable data
          IF ((curr_proc%acp==1) .AND. (lTmp .OR. (IAND(curr_var%flags,INPUT_VAR_ACTION_MASK) /= INPUT_VAR_ACTION_REPLACE))) THEN
            tmp_act(1:1) = INPUT_VAR_DO_NOACTION
            CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
          ENDIF
          ! Add rescaling
          IF (lTmp) THEN
            DO i=0,1
              buf_dta => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_RESCALE,i*INPUT_DATA_NO_STEP)
              IF (ASSOCIATED(buf_dta)) THEN
                tmp_act(1:2) = (/INPUT_VAR_DO_RESCALE,i*INPUT_DATA_NO_STEP/)
                CALL InputAddDataAction(curr_var,tmp_act(1:2),curr_proc)
              ENDIF
            ENDDO
          ENDIF
          curr_var%stat = IOR(curr_var%stat,INPUT_VAR_INITIALIZED+INPUT_VAR_EXIST+INPUT_VAR_PROC_IN_FILE)
          lInit = .TRUE.
        ENDIF
      ENDIF ! Has parent file - else

      IF ((curr_proc%acp > 1) .AND. ASSOCIATED(curr_var%dta)) THEN ! Data are obtained in some way and should be passed to model

        ! Add model_data read_data combination operator if needed
        SELECT CASE (IAND(curr_var%flags,INPUT_VAR_ACTION_MASK))
          CASE (INPUT_VAR_ACTION_REPLACE);   tmp_act(1) = INPUT_VAR_DO_NOACTION
          CASE (INPUT_VAR_ACTION_ADD);       tmp_act(1) = INPUT_VAR_DO_COMBINE_ADD
          CASE (INPUT_VAR_ACTION_SUB);       tmp_act(1) = INPUT_VAR_DO_COMBINE_SUB
          CASE (INPUT_VAR_ACTION_SUBR);      tmp_act(1) = INPUT_VAR_DO_COMBINE_SUBR
          CASE (INPUT_VAR_ACTION_MUL);       tmp_act(1) = INPUT_VAR_DO_COMBINE_MUL
          CASE (INPUT_VAR_ACTION_DIV);       tmp_act(1) = INPUT_VAR_DO_COMBINE_DIV
          CASE (INPUT_VAR_ACTION_DIVR);      tmp_act(1) = INPUT_VAR_DO_COMBINE_DIVR
          CASE (INPUT_VAR_ACTION_ADD_FLUX);  tmp_act(1) = INPUT_VAR_DO_COMBINE_ADD_FLUX
          CASE (INPUT_VAR_ACTION_NUDGE);     tmp_act(1) = INPUT_VAR_DO_COMBINE_NUDGE
          CASE (INPUT_VAR_ACTION_NUDGE_FLD); tmp_act(1) = INPUT_VAR_DO_COMBINE_NUDGE_FLD
        END SELECT
        IF (tmp_act(1) /= INPUT_VAR_DO_NOACTION) THEN
          CALL InputAddDataAction(curr_var,tmp_act(1:1),curr_proc)
          CALL InputVarNewBuffer(curr_var,curr_var%edims,curr_proc,global)
        ENDIF

        ! Clean up
        IF (ASSOCIATED(fdim)) CALL InputDimDone(fdim,.FALSE.,.FALSE.)
        IF (ASSOCIATED(global%buf)) global%buf%flags = IAND(global%buf%flags,INPUT_ALLFLAGS-INPUT_DATA_IN_USE)

        ! Now resize the action list to its actual size
        ALLOCATE(actions(curr_proc%acp-1))
        actions(:) = curr_var%saction(1:curr_proc%acp-1)
        DEALLOCATE(curr_var%saction)
        curr_var%saction => actions

      ENDIF ! Data obtained

      ! Last buffer is the model data themselves if available
      IF (lInit .AND. ASSOCIATED(curr_var%dta)) THEN
        CALL InputDataUnuse(AllInputData,dta=curr_proc%buf%dta) ! Remove last buffer which is not actually needed
        curr_proc%buf%dta => curr_var%dta
      ENDIF

      curr_var => curr_var%next
    ENDDO ! Variable loop

    InputCheckVarDims = .TRUE.

  END FUNCTION InputCheckVarDims

! Checks variable presence and dimensions in file
  RECURSIVE LOGICAL FUNCTION InputCheckVars(IFile,lOrderChg)
  TYPE (input_file_list), POINTER :: IFile
  LOGICAL, OPTIONAL,  INTENT(out) :: lOrderChg

  INTEGER, PARAMETER :: WarnMap(0:12) = (/1,1,2,2,4,4,4,1,2,1,4,1,2/)

  TYPE (input_var_list),     POINTER :: curr_var, last_var, first_var, prev_var, prev_var_ass, back_var, ins_after
  TYPE (input_fileobj_list), POINTER :: curr_obj
  TYPE (input_file_list),    POINTER :: curr_file
  INTEGER :: i, st, en, MaxLen, shap(7)
  LOGICAL :: lCont, lTest, lOC
  LOGICAL, SAVE :: lRecurse = .TRUE.
  REAL(dp), POINTER :: def_value(:,:,:,:,:,:,:)
  REAL(dp) :: fill_model(5), tmp
  CHARACTER(len=128) :: str

    InputCheckVars = .FALSE.
    message_text = ''
    curr_var => AllInputVars
    DO WHILE (ASSOCIATED(curr_var))
      IF (ASSOCIATED(curr_var%parent_file,IFile)) THEN
        ! Allocate saction as temporary buffer for file dimension ids and variable id
        IF (lRecurse) THEN ! Only the first time
          IF (ASSOCIATED(curr_var%saction)) DEALLOCATE(curr_var%saction)
          ALLOCATE(curr_var%saction(7))
          curr_var%saction(:) = INPUT_VAR_DO_NOACTION ! Prevent error when displaying variable in this state
        ENDIF

        lCont = .FALSE.
        lTest = .TRUE.
        curr_var%stat = IAND(curr_var%stat,INPUT_ALLFLAGS - INPUT_VAR_EXIST - INPUT_VAR_HAS_UNLIM - INPUT_VAR_PROC_IN_FILE)

        ! Secure that extra files for multi-file variables are open, and that the necessary variables are present in the files
        ! Notice, that presently no checks are done on the dimensionality of the files
        IF ((IAND(curr_var%stat,INPUT_VAR_MULTI) /= 0) .AND. lRecurse) THEN
          lRecurse = .FALSE. ! Prevent infinite recursions between InputFileOpen and this function
          lTest    = .FALSE. ! Switch off main test, since sub-test is performed
          ! Replace variable names including wildcards with matching names from file
          IF (ASSOCIATED(IFile%fvars)) THEN
            str = TRIM(curr_var%name_alt)
            st = 1
            DO WHILE (st > 0)
              lTest = .FALSE.
              st = INDEX(curr_var%name_alt,'*')
              IF (st > 0) THEN
                lTest = .TRUE.
              ELSE
                st = INDEX(curr_var%name_alt,'?')
              ENDIF
              IF (st > 0) CALL InputFileObjExpandVar(IFile%fvars,curr_var%name_alt,st,lTest)
            ENDDO
            IF (LEN_TRIM(curr_var%name_alt)==0) &
              CALL local_message(TRIM(curr_var%name_var),'No matches found for variable name: '//TRIM(str))
          ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
          CALL p_bcast(curr_var%name_alt,io_pe)
#endif
          MaxLen = INDEX(curr_var%name_alt,'=>')
          IF (MaxLen == 0) MaxLen = LEN_TRIM(curr_var%name_alt) + 1
          curr_file => IFile
          st = 1
          DO WHILE (st < MaxLen)
            en = INDEX(curr_var%name_alt(st:),',') + st - 1
            IF (en < st) en = MaxLen
            IF ((IAND(curr_file%stat,INPUT_FILE_OPEN) == 0) .OR. (                      &
                 IAND(curr_file%stat,INPUT_FILE_CHECK_DIM_SIZE+INPUT_FILE_CHECK_VAR) /= &
                                     INPUT_FILE_CHECK_DIM_SIZE+INPUT_FILE_CHECK_VAR))   & 
              CALL InputFileOpen(curr_file,check=INPUT_FILE_CHECK_DIM_SIZE+INPUT_FILE_CHECK_VAR)
            curr_obj => InputFileObjGetRef(curr_file%fvars,curr_var%name_alt(st:en-1))
            IF (ASSOCIATED(curr_obj)) THEN
              IF (st==1) curr_var%saction(:) = curr_obj%fids_dim(:) ! Get dimensions from first part of sub-variable
            ELSE
              lCont = .TRUE.
            ENDIF
            IF (IAND(curr_var%stat,INPUT_VAR_MULTI_SPREAD) /= 0) curr_file => curr_file%next
            st = en + 1
          ENDDO
          IF (.NOT. lCont) curr_var%stat = IOR(curr_var%stat,INPUT_VAR_EXIST) ! Exist only if all sub-vars exist
          lRecurse = .TRUE.
        ENDIF

        ! Check presence of variable
        fill_model(:) = UNLIKELY_RVAL
        IF (liodata .AND. .NOT. lCont .AND. lTest .AND. lRecurse) THEN
          lCont = .TRUE.
          curr_obj => InputFileObjGetRef(IFile%fvars,curr_var%name_var,curr_var%name_alt)
          IF (ASSOCIATED(curr_obj)) THEN
            curr_var%stat = IOR(curr_var%stat,INPUT_VAR_EXIST)
            curr_var%fid   = curr_obj%id
            curr_var%saction(:) = curr_obj%fids_dim(:)
            ! Special processing related to NetCDF attributes
            IF (ASSOCIATED(curr_obj%points)) THEN
              IF (IAND(curr_obj%flags,INPUT_FILEVAR_HAS_SCALE) /= 0) & ! Scaling given
                CALL ReplaceAuxData1D(curr_var,INPUT_DATA_RESCALE,0,curr_obj%points(1:2))
              IF (IAND(curr_var%flags,INPUT_VAR_BAD_VAL_MASK) /= INPUT_VAR_BAD_VAL_IGNORE) THEN
                def_value => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_FILL,INPUT_DATA_NO_STEP)
                IF (.NOT. ASSOCIATED(def_value)) &
                  def_value => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_FILL,0)
                IF (.NOT. ASSOCIATED(def_value)) THEN
                  def_value => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_DEFAULT,0)
                  IF (ASSOCIATED(def_value)) THEN
                    IF (SIZE(def_value) > 1) THEN
                      CALL local_message('WARNING','InputFileOpen: '//TRIM(curr_var%name_var)// &
                                                   ': Cannot use vector-default as fill value')
                      NULLIFY(def_value)
                    ENDIF
                  ENDIF
                ENDIF
                tmp = UNLIKELY_RVAL
                IF (ASSOCIATED(def_value)) tmp = def_value(1,1,1,1,1,1,1)
                IF (IAND(curr_obj%flags,INPUT_FILEVAR_HAS_BAD) /= 0) & ! Bad value given
                  fill_model(1:2) = (/tmp,curr_obj%points(3)/) !Needed on every PE -split
                IF (IAND(curr_obj%flags,INPUT_FILEVAR_HAS_RANGE) /= 0) & ! Valid range given
                  fill_model(3:5) = (/tmp,curr_obj%points(4:5)/) !Needed on every PE -split
              ENDIF
            ENDIF
            lCont = .FALSE.
          ENDIF
        ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
        CALL p_bcast(curr_var%flags,io_pe)
        CALL p_bcast(curr_var%stat, io_pe)
        CALL p_bcast(lCont,         io_pe)
        CALL p_bcast(curr_var%fid,  io_pe)
        CALL p_bcast(fill_model,    io_pe)
        IF (lRecurse) CALL p_bcast(curr_var%saction,io_pe)
#endif
        IF (lCont) THEN ! Variable not found
          i = IAND(curr_var%flags,INPUT_VAR_MISS_MASK)
          IF  (curr_var%nRead > 0) i = i + 6
          IF  (WarnMap(i) >  1) WRITE(message_text,'(3a)') 'Variable not found in ',TRIM(curr_var%parent_file%name_file)
          IF  (WarnMap(i) == 4) CALL local_error('InputFileOpen','Error: '//TRIM(curr_var%name_var)//': '//TRIM(message_text))
          IF ((WarnMap(i) <  4) .AND. (curr_var%nIntrp > 0)) &
            curr_var%Intrp%next%dta(:,:,:,:,:,:,:) = curr_var%Intrp%dta(:,:,:,:,:,:,:)
          def_value => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_DEFAULT,0)
          IF (ASSOCIATED(def_value) .AND. ASSOCIATED(curr_var%dta)) THEN
            IF (SIZE(def_value) == 1) THEN
              curr_var%dta(:,:,:,:,:,:,:) = def_value(1,1,1,1,1,1,1)
              WRITE(message_text(LEN_TRIM(message_text)+1:),'(a,g10.5)') '. Using default value: ',def_value(1,1,1,1,1,1,1)
            ELSE
              IF ((SIZE(def_value) > 1) .AND. ALL(SHAPE(curr_var%dta) /= SIZE(def_value))) THEN
                shap = SHAPE(curr_var%dta)
                i = 7
                DO WHILE (shap(i) == 1)
                  i = i - 1
                ENDDO
                WRITE(message_text,*) ': Vector default length ',SIZE(def_value),' does not match size of any dimension: ', &
                                      shap(1:i)
                CALL RmDblSpc(message_text)
                CALL local_error('InputFileOpen',TRIM(curr_var%name_Var)//TRIM(message_text))
              ENDIF
              IF (.NOT. ArrayCopy(def_value,curr_var%dta,(/SIZE(def_value)/),SHAPE(curr_var%dta))) &
                    CALL local_error('InputFileOpen','Spread error')
              WRITE(message_text(LEN_TRIM(message_text)+1:),'(a)') '. Using values from vector-defaults.'
            ENDIF
          ENDIF
          IF ((WarnMap(i) > 1) .AND. (IAND(mo_debug,INPUT_DBG_MASK) /= 0)) &
            CALL local_message('','WARNING: '//TRIM(curr_var%name_var)//': '//TRIM(message_text))
        ELSE
          def_value => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_FILL,INPUT_DATA_NO_STEP)
          IF (ASSOCIATED(def_value)) THEN
            IF (def_value(1,1,1,1,1,1,1) /= UNLIKELY_RVAL) fill_model(1) = def_value(1,1,1,1,1,1,1)
            IF (def_value(2,1,1,1,1,1,1) /= UNLIKELY_RVAL) fill_model(2) = def_value(2,1,1,1,1,1,1)
          ENDIF
          def_value => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_RANGE,INPUT_DATA_NO_STEP)
          IF (ASSOCIATED(def_value)) THEN
            IF (def_value(1,1,1,1,1,1,1) /= UNLIKELY_RVAL) fill_model(3) = def_value(1,1,1,1,1,1,1)
            IF (def_value(2,1,1,1,1,1,1) /= UNLIKELY_RVAL) fill_model(4) = def_value(2,1,1,1,1,1,1)
            IF (def_value(3,1,1,1,1,1,1) /= UNLIKELY_RVAL) fill_model(5) = def_value(3,1,1,1,1,1,1)
          ENDIF
          IF (fill_model(2) /= UNLIKELY_RVAL) & ! Fill value now on every PE
            CALL ReplaceAuxData1D(curr_var,INPUT_DATA_FILL,0,fill_model(1:2))
          IF (fill_model(4) /= UNLIKELY_RVAL) & ! Fill value now on every PE
            CALL ReplaceAuxData1D(curr_var,INPUT_DATA_RANGE,0,fill_model(3:5))
        ENDIF
      ENDIF
      curr_var => curr_var%next
    ENDDO

    ! Reorder variables so that the processing order is the same as the data order in file - saves time for strided reading
    lOC = .FALSE.
    NULLIFY(first_var,last_var,prev_var,prev_var_ass,back_var,ins_after)
    curr_var => AllInputVars
    DO WHILE (ASSOCIATED(curr_var))
      lCont = .TRUE.
      IF (ASSOCIATED(curr_var%parent_file,IFile)) THEN ! curr_var belongs to file
        IF (ASSOCIATED(first_var)) THEN                ! File already have associated var(s)
          IF (curr_var%fid <= first_var%fid) THEN      ! curr_var should be before the first in file
            lOC = .TRUE.
            ins_after => prev_var_ass                  ! Ins curr_var after the var before the first in this file (may be undef)
            first_var => curr_var                      ! curr_var is now the first in file
          ELSE
            IF (curr_var%fid > last_var%fid) THEN      ! curr_var should be after the last in file
              ins_after => last_var                    ! insert curr_var after last in file
              last_var  => curr_var                    ! curr_var is now the last in file
            ELSE                                       ! curr_var is somewhere in the middle of the file
              ins_after => first_var                   ! test from the first in file
              lCont = ASSOCIATED(ins_after%next)
              DO WHILE (lCont)
                lCont = ASSOCIATED(ins_after%next)     ! For security,should never return false, since last is take care of above
                IF (lCont) lCont = curr_var%fid > ins_after%next%fid ! Test until next should be after curr_var
                IF (lCont) ins_after => ins_after%next ! More tests to come, go to next var to test
              ENDDO
            ENDIF
          ENDIF
          lCont = .FALSE.                              ! Normal progress to next variable is not the default
          IF (ASSOCIATED(ins_after)) THEN              ! Not the first variable in system
            IF (ASSOCIATED(ins_after%next,curr_var)) lCont = .TRUE. ! curr_var happens to be at the right position, just progress
          ELSE
            back_var     => AllInputVars               ! Save first variable in system for later
            AllInputVars => curr_var                   ! curr_var should be the first in the list 
          ENDIF
          IF (.NOT. lCont) THEN                        ! curr_var should be moved            
            lOC = .TRUE.
            prev_var%next => curr_var%next             ! curr_var is taken out of the list
            IF (ASSOCIATED(ins_after)) THEN            ! curr_var is not the first in system
              curr_var%next  => ins_after%next         ! Secure forward link
              ins_after%next => curr_var               ! Insert curr_var
            ELSE
              curr_var%next => back_var                ! Secure forward link
            ENDIF
            curr_var => prev_var%next                  ! This is the next to process
          ENDIF
        ELSE ! curr_var is the first found variable belonging to this file - that is curr_var is in place
          prev_var_ass => prev_var ! Last variable before those belonging to this file
          first_var    => curr_var ! first_var points to first found variable in this file
          last_var     => curr_var ! last_var  points to last  found variable in this file
          lCont = .TRUE.           ! Go to next variable the normal way
        ENDIF
      ENDIF
      IF (lCont) THEN ! Normal list progressing
        prev_var => curr_var
        curr_var => curr_var%next
      ENDIF
    ENDDO
    IF (PRESENT(lOrderChg)) lOrderChg = lOC

    InputCheckVars = .TRUE.

  END FUNCTION InputCheckVars

! Checks dimensions, variables and time settings for a file
  RECURSIVE SUBROUTINE InputFileCheck(IFile,rcheck,lOrderChg)
  TYPE (input_file_list), POINTER :: IFile
  INTEGER,            INTENT(in)  :: rcheck
  LOGICAL, OPTIONAL,  INTENT(out) :: lOrderChg

  TYPE (input_fileobj_list), POINTER :: new_obj
  TYPE (input_var_list), POINTER :: curr_var
  TYPE (input_opt_list), POINTER :: opt
  TYPE (input_dim_list), POINTER :: udim
  INTEGER :: cnt, tt(4)
  LOGICAL :: lCont

    ! Check dimensions
    IF ((IAND(rcheck,INPUT_FILE_CHECK_DIM_SIZE) /= 0) .AND. (IAND(IFile%stat,INPUT_FILE_CHECK_DIM_SIZE) == 0)) THEN
      IF (liodata .AND. (IAND(IFile%stat,INPUT_FILE_DIMS_GOT)==0)) THEN
#ifdef HAVE_F2003
        new_obj => IFile%ft%FcnGetDimensions(IFile,31)   ! Get dimensions from file
#else
        new_obj => NCFileGetDimensions(IFile,31)         ! Get dimensions from file
#endif
        IF (InputFileObjsCmp(IFile%fdims,new_obj)) THEN ! If file dimensions are the same, we can skip some checks
          CALL InputFileObjsDone(new_obj)
          IFile%stat = IOR(IFile%stat,INPUT_FILE_CHECK_DIM_SIZE)
        ELSE
          IF (ASSOCIATED(IFile%fdims)) CALL InputFileObjsDone(IFile%fdims)
          IFile%fdims => new_obj
          CALL InputFileGetExternalSettings(IFile,opt,cnt)
          CALL InputFileApplyExternalSettings(IFile,opt,cnt)
        ENDIF
      ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      IF (n_pe >= 0) CALL p_bcast(IFile%flags,io_pe)
      IF (n_pe >= 0) CALL p_bcast(IFile%stat,io_pe)
#endif
      IF (IAND(IFile%stat,INPUT_FILE_CHECK_DIM_SIZE) == 0) THEN
        IF (.NOT. InputCheckDimSizes(IFile,IFile%sub_model)) & ! Check or set model dimension sizes agains file dim. sizes
          CALL local_error('InputFileOpen',message_text)
        IFile%stat = IOR(IFile%stat,INPUT_FILE_CHECK_DIM_SIZE)
      ENDIF
      IF (IAND(IFile%stat,rcheck) == rcheck) RETURN  ! All requested checks have been performed
    ENDIF

    ! Check if required variables are present
    IF (((IAND(rcheck,INPUT_FILE_CHECK_VAR     ) /= 0) .AND. (IAND(IFile%stat,INPUT_FILE_CHECK_VAR     ) == 0)) .OR.  &
        ((IAND(rcheck,INPUT_FILE_CHECK_VAR_DIMS) /= 0) .AND. (IAND(IFile%stat,INPUT_FILE_CHECK_VAR_DIMS) == 0))) THEN
      IF (liodata .AND. (IAND(IFile%stat,INPUT_FILE_VARS_GOT) == 0)) THEN
#ifdef HAVE_F2003
        new_obj => IFile%ft%FcnGetVariables(IFile)
#else        
        new_obj => NCFileGetVariables(IFile)
#endif
        IF (InputFileObjsCmp(IFile%fvars,new_obj)) THEN ! If file variables are the same, we can skip some checks
          IFile%stat = IOR(IFile%stat,INPUT_FILE_CHECK_VAR)
          CALL InputFileObjsDone(new_obj)
        ELSE
          IF (ASSOCIATED(IFile%fdims)) CALL InputFileObjsDone(IFile%fvars)
          IFile%fvars => new_obj
        ENDIF
      ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      IF (n_pe >= 0) THEN
        CALL p_bcast(IFile%flags,io_pe)
        CALL p_bcast(IFile%stat,io_pe)
      ENDIF
#endif

      IF ((IAND(IFile%stat,INPUT_FILE_CHECK_VAR) == 0) .OR. (IAND(IFile%stat,INPUT_FILE_CHECK_VAR_DIMS) == 0)) THEN
        IF (.NOT. InputCheckVars(IFile,lOrderChg)) &
          CALL local_error('InputFileOpen','One or more required variables not found in '//TRIM(IFile%name_file))
        IF (PRESENT(lOrderChg)) THEN
          IF (lOrderChg) CALL SortDependencies() ! If variable order changed, make sure that dependencies are properly resolved
        ENDIF
        IFile%stat = IOR(IFile%stat,INPUT_FILE_CHECK_VAR)
      ENDIF
      IF (IAND(IFile%stat,rcheck) == rcheck) RETURN  ! All requested checks have been performed
    ENDIF

    ! Number of records along the unlimited dimension axis in file
    IF (lio) THEN
      lCont = .FALSE.
      udim => InputGetUnlimitedDim(IFile)
      IF (ASSOCIATED(udim)) THEN
#ifdef HAVE_F2003
        IFile%nRec = IFile%ft%FcnGetNRec(IFile,udim)
#else
        IFile%nRec = NCFileGetNRec(IFile,udim)
#endif
      ELSE
        IF (IAND(IFile%flags,INPUT_FILE_INITIAL) == 0) THEN
          lCont = .TRUE. ! Flags are to be changed on every PE
          curr_var => AllInputVars
          DO WHILE (ASSOCIATED(curr_var))
            IF (ASSOCIATED(curr_var%parent_file,IFile) .AND.                                    &
              (curr_var%dt_update /= 0) .AND. (curr_var%dt_update /= INPUT_VAR_FILE_TIME_STEP)) &
                CALL local_error('InputFileOpen','No unlimited dimension present in '//TRIM(IFile%name_file))
            curr_var => curr_var%next
          ENDDO
        ENDIF
      ENDIF
    ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
    CALL p_bcast(IFile%nRec, io_pe)
    CALL p_bcast(lCont,      io_pe)
#endif

    IF (lCont) THEN
      IF (IAND(mo_debug,INPUT_DBG_MASK) /= 0) CALL local_message('WARNING',TRIM(BuildGenericName(IFile))// &
                                                ' contains no unlimited dimension. Assuming it is an initial file')
      IFile%flags = IOR(IFile%flags,INPUT_FILE_INITIAL) ! No unlimited dimension found - must be an initial file
      curr_var => AllInputVars
      DO WHILE (ASSOCIATED(curr_var))
        IF (ASSOCIATED(curr_var%parent_file,IFile) .AND. (curr_var%dt_update == INPUT_VAR_FILE_TIME_STEP)) &
          curr_var%dt_update = 0
        curr_var => curr_var%next
      ENDDO
    ENDIF

    ! Obtain time step in file and loop end (only relevant for non-initial files - i.e. files without non-initial variables)
    lCont = ((IAND(rcheck,INPUT_FILE_CHECK_TIME_STEP) /= 0) .AND. (IAND(IFile%flags,INPUT_FILE_INITIAL) == 0)) .AND. &
             ((IFile%dt_file == 0) .OR. (IFile%dt_file == UNLIKELY_VAL))
    IF (lCont) THEN
      lCont = .FALSE.
      curr_var => AllInputVars
      DO WHILE (ASSOCIATED(curr_var))
        IF (ASSOCIATED(curr_var%parent_file,IFile) .AND. (curr_var%dt_update /= 0)) THEN 
          lCont = .TRUE.
          NULLIFY(curr_var)
        ELSE
          curr_var => curr_var%next
        ENDIF
      ENDDO
    ENDIF

    IF (lCont) THEN

      ! Time step (space between successive entries on the unlimited dimension axis) of data in file
      IFile%dt_file = InputFileTimeStepGet(IFile,.TRUE.)
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      IF (n_pe >= 0) CALL p_bcast(IFile%dt_file,io_pe)
#endif
      ! If loop end is within this file do as if the file was appropriately shorter
      IF (IFile%time_enddata(5) > 0) THEN
        tt = IFile%time_open(1:4)
        CALL CalendarTimeAdd(tt,INT(IFile%dt_file,i8)*INT(IFile%nRec-1,i8))
        IF (CalendarTimeCmp(tt,IFile%time_enddata(1:4)) <= 0) IFile%nRec = IFile%time_enddata(5) - 1
      ENDIF

    ENDIF

    ! Check variable dimensions
    IF ((IAND(rcheck,INPUT_FILE_CHECK_VAR_DIMS) /= 0) .AND. (IAND(IFile%stat,INPUT_FILE_CHECK_VAR_DIMS) == 0)) THEN
      IF (.NOT. InputCheckVarDims(IFile)) &
        CALL local_error('InputFileOpen','Variable dimension check failed in '//TRIM(IFile%name_file))
      IFile%stat = IOR(IFile%stat,INPUT_FILE_CHECK_VAR_DIMS)
    ENDIF

  END SUBROUTINE InputFileCheck

! ---------------------------------------------------
!  Finding time of single files of a multifile
! ---------------------------------------------------

! Get time step of file (file must be open)
  INTEGER FUNCTION InputFileTimeStepGet(IFile,lKeep)
  TYPE (input_file_list), POINTER :: IFile
  LOGICAL, INTENT(in)             :: lKeep ! Reopen this exact file even if it was necessary to access another

  TYPE (input_dim_list),  POINTER :: udim
  INTEGER :: tmp_time(4), file_time(4), curr_time(4)
  INTEGER :: res
  INTEGER (KIND=i8) :: res64
  LOGICAL :: lCont, lSuc, lExist
  CHARACTER (len=128) :: FileName

    ! Make sure, the number of records in file is available
    lExist = .TRUE.
    IF (lio .AND. (IFile%nRec <= 0)) THEN ! Not initialized yet
      lExist = .FALSE.
      udim => InputGetUnlimitedDim(IFile)
      IF (ASSOCIATED(udim)) THEN
        lExist = .TRUE.
#ifdef HAVE_F2003
        IFile%nRec = IFile%ft%FcnGetNRec(IFile,udim)
#else
        IFile%nRec = NCFileGetNRec(IFile,udim)
#endif
      ELSE 
        lExist = .FALSE.
        IFile%nRec = 1
      ENDIF
      IF (IFile%nRec < 0) CALL local_error('InputFileTimeStepGet',&
                                  'Unlimited dimension '//TRIM(udim%dim_data%name_dim)//' not present in '//TRIM(IFile%name_file))
    ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
    CALL p_bcast(IFile%nRec,io_pe)
    CALL p_bcast(lExist,io_pe)
#endif

    ! Time step can only be determined if at least two time steps are available (may be in different files)
    IF (lExist) THEN
      IF (IFile%nRec < 2) THEN
        IF (IAND(IFile%flags,INPUT_FILE_MULTI) == 0) &
          CALL local_error('InputFileTimeStepGet',TRIM(IFile%name_file)// &
                                                  ': At least two records must be present in file for autodetermining time step')
#ifdef HAVE_F2003
        IF (lio) tmp_time  = IFile%ft%FcnGetRecTime(IFile,1,.FALSE.)
#else
        IF (lio) tmp_time  = NCFileGetRecTime(IFile,1,.FALSE.)
#endif
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
        CALL p_bcast(tmp_time,io_pe)
#endif
        curr_time = IFile%time_open
        ! Find neighbouring file: First look for files after the current - if unavailable, look before
        lSuc  = .TRUE.
        lCont = .TRUE.
        DO WHILE (lCont)
          IF (lSuc) THEN
            curr_time = InputFileTimeNextGet(IFile,curr_time)
            IF (CalendarTimeCmp(curr_time,IFile%time_enddata(1:4)) < 0) THEN
              lSuc = .FALSE.
              curr_time = IFile%time_open
            ENDIF
          ENDIF
          IF (.NOT. lSuc) THEN
            curr_time = InputFileTimePrevGet(IFile,curr_time)
            IF (CalendarTimeCmp(curr_time,IFile%time_start(1:4)) > 0) &
              CALL local_error('InputFileTimeStepGet','No other file available - cannot determine time step')
          ENDIF
          FileName = GetFileName(IFile,curr_time(1),curr_time(2),curr_time(3),curr_time(4))
          INQUIRE(file=FileName,exist=lExist)
          IF (lExist) lCont = .FALSE.                  
        ENDDO
        file_time = IFile%time_open ! Save to be able to re-open
        CALL InputFileOpen(IFile,curr_time(1),curr_time(2),curr_time(3),Check = 0)
#ifdef HAVE_F2003
        IF (lio) curr_time = IFile%ft%FcnGetRecTime(IFile,1,.FALSE.)
#else
        IF (lio) curr_time = NCFileGetRecTime(IFile,1,.FALSE.)
#endif
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
        CALL p_bcast(curr_time,io_pe)
#endif
        IF (lKeep) CALL InputFileOpen(IFile,file_time(1),file_time(2),file_time(3),Check = 0)
        IF (lSuc) THEN
          res64 = CalendarTimeDiff(tmp_time,curr_time)
        ELSE
          res64 = CalendarTimeDiff(curr_time,tmp_time)
        ENDIF
        lExist = .TRUE.
      ELSE
#ifdef HAVE_F2003
        IF (lio) tmp_time = IFile%ft%FcnGetRecTime(IFile,1,.TRUE.)
#else
        IF (lio) tmp_time = NCFileGetRecTime(IFile,1,.TRUE.)
#endif
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
        CALL p_bcast(tmp_time,io_pe)
#endif
        res64 = INT(tmp_time(4),i8)
      ENDIF
      IF (res64 >= GREGORIAN_MIN_MONTH) THEN
        res64= -INT(REAL(res64,dp) / REAL(GREGORIAN_MEAN_MONTH,dp) + 0.5_dp,i8)
      ENDIF
      res = INT(res64)
      IF (res == 0) CALL local_error('InputFileTimeStepGet','File: '//TRIM(IFile%name_file)//': Time step could not be determined')
    ELSE
      InputFileTimeStepGet = 0
    ENDIF
    IFile%stat = IOR(IFile%stat,INPUT_FILE_TIME_GOT)

    InputFileTimeStepGet = res

  END FUNCTION InputFileTimeStepGet

! --------------------------------------------
!  Reading and updating of relevant variables
! --------------------------------------------

! Determines start, end and cycle time for a file
  LOGICAL FUNCTION InputSetFileDataLimits(IFile,lOrderChg)
  TYPE (input_file_list), POINTER :: IFile
  LOGICAL,            INTENT(out) :: lOrderChg

  TYPE (input_opt_list),  POINTER :: opt
  INTEGER(i8) :: dt
  INTEGER :: tmf_curr(4), tmf_test(4)
  INTEGER :: Tmp, Sel, am1, am2
  LOGICAL :: lCont, lOC, lFirstLast

    InputSetFileDataLimits = .FALSE.
    lOrderChg = .FALSE.
    IF (IAND(IFile%stat,INPUT_FILE_EXTERN_APPLIED) == 0) THEN
      CALL InputFileGetExternalSettings(IFile,opt,Tmp)
      CALL InputFileApplyExternalSettings(IFile,opt,Tmp)
    ENDIF

    IF (ANY(IFile%time_start > -HUGE(Tmp))) THEN
      tmf_curr = InputFileTimeCurrGet(IFile,IFile%time_start(1:4))
    ELSE
      tmf_curr = InputFileMultiInstanceGet(IFile,1)
    ENDIF
    IF (ANY(tmf_curr==-HUGE(Tmp))) RETURN

    IF (IAND(IFile%flags,INPUT_FILE_TIME_ABSOLUTE) /= 0) THEN

      ! Find first real available time and record/file of first user available
      CALL InputFileOpen(IFile,tmf_curr(1),tmf_curr(2),tmf_curr(3),lOrderChg=lOC)
      IF (IAND(IFile%stat,INPUT_FILE_OPEN) == 0) RETURN ! Probably no instance found
      lOrderChg = lOrderChg .OR. lOC
#ifdef HAVE_F2003
      IF (lio) tmf_test = IFile%ft%FcnGetRecTime(IFile,1,.FALSE.) ! Time of first record in found file
#else
      IF (lio) tmf_test = NCFileGetRecTime(IFile,1,.FALSE.)
#endif
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      IF (n_pe >= 0) CALL p_bcast(tmf_test,io_pe)
#endif
      dt = CalendarTimeDiffSec(tmf_test,tmf_curr) ! Time difference between time given by file name and first record in file
      IF (ANY(IFile%time_start > -HUGE(i8))) THEN ! time_start set?
        IF (IAND(IFile%flags,INPUT_FILE_MULTI) /= 0) THEN
          IF (CalendarTimeCmp(IFile%time_start(1:4),tmf_test) < 0) &     ! Record time before first time to use?
            tmf_curr = InputFileTimeCurrGet(IFile,IFile%time_start(1:4)) ! File time of first time to use
        ENDIF
        CALL CalendarTimeAddSec(tmf_curr,dt) ! tmf_curr is now the record time of first rec in the relevant file
        lCont = CalendarTimeCmp(IFile%time_start(1:4),tmf_curr) < 0 ! Still before first valid time?
        IFile%time_start(5) = 1
        DO WHILE (lCont) ! Loop until first record to use is reached
          IFile%time_start(5) = IFile%time_start(5) + 1
          CALL CalendarTimeAdd(tmf_curr,INT(IFile%dt_file,i8))
          lCont = CalendarTimeCmp(IFile%time_start(1:4),tmf_curr) < 0
        ENDDO
      ELSE ! time_start not set
        IF (IAND(IFile%flags,INPUT_FILE_MULTI) /= 0) THEN
          CALL CalendarTimeAddSec(tmf_curr,dt) ! Exact time of first record
        ELSE ! tmf_curr is not reliable for non-multi-files, since it depends on the model time
          tmf_curr = tmf_test
        ENDIF
        IFile%time_start(5) = 1 ! Valid from first record
      ENDIF
      IFile%time_start(1:4) = tmf_curr       
      CALL CalendarTimeSub(tmf_curr,INT(IFile%dt_file,i8)*INT(IFile%time_start(5)-1,i8)) ! Adjust time for time_enddata calc.
    ENDIF

    ! Determine hypothetical time of first non-existing data if not already done and cyclicity is requested
    IF ((IAND(IFile%stat,INPUT_FILE_ENDDATE_GOT) == 0)) THEN ! Only needed for cyclic data!?
      IF (ANY(IFile%time_enddata(1:4) > -HUGE(i8))) THEN ! time_enddata set?
        tmf_curr = InputFileTimeCurrGet(IFile,IFile%time_enddata(1:4))
      ELSE
        tmf_curr = InputFileMultiInstanceGet(IFile,2,IFile%time_start(1:4))
      ENDIF

      IF (ANY(tmf_curr==-HUGE(Tmp))) RETURN

      ! Now we known the name of the last file: Open it and determine time of first non-existing record
      tmf_test = IFile%time_open
      CALL InputFileOpen(IFile,tmf_curr(1),tmf_curr(2),tmf_curr(3),lOrderChg=lOC)
      lOrderChg = lOrderChg .OR. lOC
#ifdef HAVE_F2003
      IF (lio) tmf_curr = IFile%ft%FcnGetRecTime(IFile,1,.FALSE.)
#else          
      IF (lio) tmf_curr = NCFileGetRecTime(IFile,1,.FALSE.)
#endif
      lFirstLast = CalendarTimeCmp(tmf_test,IFile%time_open) == 0
      IF (lFirstLast) THEN
        IFile%name_pre = TRIM(IFile%name_file)
        IFile%name_suf = ''
        IFile%flags = IAND(IFile%flags,INPUT_ALLFLAGS - INPUT_FILE_MULTI) ! First and last file are the same - don't treat as multi
      ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      IF (n_pe >= 0) CALL p_bcast(tmf_curr,io_pe)
#endif
      IF (ANY(IFile%time_enddata(1:4) > -HUGE(i8))) THEN ! time_enddata set?
        CALL CalendarTimeAdd(IFile%time_enddata(1:4),INT(IFile%dt_file,i8)) ! User set last valid time - we want first invalid
        IFile%time_enddata(5) = 0
        IF (lFirstLast) IFile%time_enddata(5) = IFile%time_start(5)
        tmf_curr = IFile%time_start(1:4)
        lCont = CalendarTimeCmp(tmf_curr,IFile%time_enddata(1:4)) >= 0
        DO WHILE (lCont) ! Loop as long as record is before time_enddata
          IFile%time_enddata(5) = IFile%time_enddata(5) + 1
          CALL CalendarTimeAdd(tmf_curr,INT(IFile%dt_file,i8))
          lCont = CalendarTimeCmp(tmf_curr,IFile%time_enddata(1:4)) >= 0
        ENDDO
        CALL CalendarTimeSub(tmf_curr,INT(IFile%dt_file,i8))
      ELSE ! time_enddata not set
        CALL CalendarTimeAdd(tmf_curr,INT(IFile%dt_file,i8)*INT(IFile%nRec,i8)) ! Assume all data in file should be used
        IF (IAND(IFile%flags,INPUT_FILE_MULTI) /= 0) THEN
          IFile%time_enddata(5) = 1 ! First non-existing rec is the first in the first non-existing file
        ELSE
          IFile%time_enddata(5) = IFile%nRec + 1 ! First non-existing rec is the one after the end of the file
        ENDIF
      ENDIF
      IFile%time_enddata(1:4) = tmf_curr
      IFile%stat = IOR(IFile%stat,INPUT_FILE_ENDDATE_GOT)

      ! Find total number of available records
      IF (IFile%dt_file < 0) THEN
        am1 = IFile%time_start  (1)*12 + IFile%time_start  (2) 
        am2 = IFile%time_enddata(1)*12 + IFile%time_enddata(2)
        IFile%time_enddata(6) = (am1 - am2) / IFile%dt_file
      ELSE
        dt = CalendarTimeDiffSec(IFile%time_enddata(1:4),IFile%time_start(1:4))
        IFile%time_enddata(6) = INT(dt / INT(IFile%dt_file,i8))
      ENDIF
    ENDIF ! Determine time_enddata

    ! Determine cycle time if requested
    IF (lio) THEN
      IF (IFile%time_cycle(2) == -HUGE(i8)) THEN   ! time_cycle not yet set?
        IF (IFile%time_cycle(1) /= -HUGE(i8)) THEN ! time_cycle holds unprocessed info to construct the cycle 
          IF (IFile%time_cycle(1) == INPUT_FILE_REPEAT_NONE) THEN
            IFile%flags = IOR(IFile%flags,INPUT_FILE_NOCYCLE)
          ELSE
            IFile%flags = IAND(IFile%flags,INPUT_ALLFLAGS - INPUT_FILE_NOCYCLE)
          ENDIF
        ENDIF
        IF (IAND(IFile%flags,INPUT_FILE_NOCYCLE) == 0) THEN ! Cyclic data used?
          Sel = IFile%time_cycle(1) ! Type of cycle length
          Tmp = IFile%time_cycle(4) ! Cycle length for those types needing it
          IF (Tmp<0) Tmp = 1        ! Map dummy defaults to a sensible value
          ! Get time step of file (updated within the IFile structure)
          IF ((IFile%dt_file == 0) .OR. (IFile%dt_file == UNLIKELY_VAL)) &
            CALL InputFileOpen(IFile,IFile%time_start(1),IFile%time_start(2),IFile%time_start(3),lOrderChg=lOC)
          lOrderChg = lOrderChg .OR. lOC
          SELECT CASE (Sel)
            CASE (INPUT_FILE_REPEAT_ALL)
              IFile%time_cycle(1:5) = IFile%time_start(1:5)
              IFile%time_cycle(6)   = IFile%time_enddata(6)
            CASE (INPUT_FILE_REPEAT_YEAR)
              IFile%time_cycle(1:4) = IFile%time_enddata(1:4)
              IFile%time_cycle(1)   = IFile%time_cycle(1) - Tmp
            CASE (INPUT_FILE_REPEAT_DAY)
              IFile%time_cycle(1:4) = IFile%time_enddata(1:4)
              CALL CalendarTimeSub(IFile%time_cycle(1:4),INT(86400,i8)*INT(Tmp,i8))
            CASE (INPUT_FILE_REPEAT_SECS)
              IFile%time_cycle(1:4) = IFile%time_enddata(1:4)
              CALL CalendarTimeSub(IFile%time_cycle(1:4),INT(Tmp,i8))
            CASE (INPUT_FILE_REPEAT_RECS)
              IFile%time_cycle(1:4) = IFile%time_enddata(1:4)
              CALL CalendarTimeSub(IFile%time_cycle,INT(IFile%dt_file,i8)*INT(Tmp,i8))
              IFile%time_cycle(6)   = Tmp
            CASE (INPUT_FILE_REPEAT_MONTHS)
              IFile%time_cycle(1:4) = IFile%time_enddata(1:4)
              CALL CalendarTimeSub(IFile%time_cycle(1:4),INT(-1,i8)*INT(Tmp,i8))
          END SELECT
          ! Determine cycle length in records if not already done
          IF (IFile%time_cycle(6) < 0) THEN
            IF (IFile%dt_file < 0) THEN
              am1 = IFile%time_cycle  (1)*12 + IFile%time_cycle  (2) 
              am2 = IFile%time_enddata(1)*12 + IFile%time_enddata(2)
              IFile%time_cycle(6) = (am1 - am2) / IFile%dt_file
            ELSE
              dt = CalendarTimeDiffSec(IFile%time_enddata(1:4),IFile%time_cycle(1:4))
              IFile%time_cycle(6) = INT(dt / INT(IFile%dt_file,i8))
            ENDIF
          ENDIF

          ! Determine start record of cycle time
          IF (Sel /= INPUT_FILE_REPEAT_ALL) THEN
            tmf_curr = InputFileTimeCurrGet(IFile,IFile%time_cycle(1:4)) ! Time of file holding the 
            CALL CalendarTimeAddSec(tmf_curr,dt) ! Adjust from file name to rec time.
            IFile%time_cycle(5) = 0
            DO WHILE (CalendarTimeCmp(tmf_curr,IFile%time_cycle(1:4)) >= 0) ! Loop until a record after time_cycle found
              IFile%time_cycle(5) = IFile%time_cycle(5) + 1
              CALL CalendarTimeAdd(tmf_curr,INT(IFile%dt_file,i8))
            ENDDO
          ENDIF
        ELSE ! no cyclic data used
          IFile%time_cycle(1:4) = -HUGE(i8)
          IFile%time_cycle(5:7) = -1
        ENDIF
      ENDIF
    ENDIF

    InputSetFileDataLimits = .TRUE.

  END FUNCTION InputSetFileDataLimits

! Determine start and cycle times of all variables read at regular intervals
  SUBROUTINE InputSetVarStartTime()
  TYPE (input_var_list),  POINTER :: var
  TYPE (input_file_list), POINTER :: IFile
  TYPE (input_data_list), POINTER :: dta
  INTEGER :: tmm_valid(4), tmm_test(4), tmm_var(4)                                      ! Time variables holding model times
  INTEGER :: tmf_valid(4), tmf_test(4), tmf_var(4), tmf_mod(4), tmf_rec(4), tmf_file(4) ! Time variables holding file  times
  INTEGER :: Tmp, nRec, nSub
  INTEGER(i8) :: dt
  LOGICAL  :: lOC, lOrderChg, lIntrp, lEndReached, lCont

    var => AllInputVars
    DO WHILE (ASSOCIATED(var))

      lCont = .FALSE.
      ! Tests to do if a file is associated with the variable
      IF ((var%sub_model == CurrSubModel) .AND. ASSOCIATED(var%parent_file)) THEN

        lCont = .TRUE.
        lOrderChg = .FALSE.

        ! Get time limits for file if not already done
        IFile => var%parent_file
        IF ((IAND(IFile%stat,INPUT_FILE_EXTERN_APPLIED+INPUT_FILE_ENDDATE_GOT) &
                          /= INPUT_FILE_EXTERN_APPLIED+INPUT_FILE_ENDDATE_GOT) &
            .OR. (ANY(IFile%time_start == -HUGE(i8)))                          & 
            .OR. (      (IFile%time_cycle(1) /= -HUGE(i8))                     &
                  .AND. (IFile%time_cycle(2) == -HUGE(i8)))) THEN
          IF (.NOT. InputSetFileDataLimits(IFile,lOrderChg)) lCont = .FALSE.
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
          IF (n_pe >= 0) THEN
            CALL p_bcast(IFile%time_start,  io_pe)
            CALL p_bcast(IFile%time_cycle,  io_pe)
            CALL p_bcast(IFile%time_enddata,io_pe)
          ENDIF
#endif
        ENDIF
      ENDIF

      IF (lCont) THEN

        ! Prepare interpolations if requested and not already done
        dt = INT(IFile%dt_file,i8)
        IFile%rec_step = 1
        lIntrp = .FALSE.
        IF (INT(var%dt_update,i8)*INT(IFile%dt_file,i8) > 0) THEN ! dt_file and dt_update have same unit (32bit ints may overflow)
          lIntrp = ABS(var%dt_update) < ABS(IFile%dt_file)
          IF (ABS(var%dt_update) > ABS(IFile%dt_file)) dt = INT(var%dt_update,i8)
          IFile%rec_step = MAX(1,ABS(var%dt_update)/ABS(IFile%dt_file))
        ELSE
          IF (var%dt_update < 0) CALL local_error('InputSetVarStartTime',TRIM(var%name_var)// &
                                   ': Sorry, update time in months and file time step in seconds currently not supported')
          lIntrp = .TRUE.
        ENDIF

        IF (lIntrp) THEN ! Should interpolation be used?
          IF (var%nIntrp == 0) var%nIntrp = 2
          tmp = 0
          IF (var%dt_update < 0)                       tmp = tmp + 1
          IF (IFile%dt_file < 0)                       tmp = tmp + 2
          IF (ABS(IFile%dt_file) < ABS(var%dt_update)) tmp = tmp + 4
          SELECT CASE (tmp)
            CASE (0) ! Both given in seconds, file > update
              var%nupdate = var%dt_update / model_dt
            CASE (3) ! Both given in months, file > update. nupdate must be determined from month to month
            CASE (1,5) ! update in mth, file in sec
              CALL local_error('InputSetVarStartTime', &
                               'Sorry. Cannot handle update time steps in months while file time steps is in seconds')
            CASE (2,6) ! file in mth, update in sec
              IF (var%dt_update >= GREGORIAN_MIN_MONTH) CALL local_error('InputSetVarStartTime', &
                                                       'Relation between file time step and variable update step undeterminable')
              var%nupdate = var%dt_update / model_dt
            CASE (4,7) ! Both in seconds or both in month, update > file
              IF (     ((var%dt_update < 0) .AND. (MOD(IFile%dt_file,var%dt_update)) /= 0)      &
                  .OR. ((var%dt_update > 0) .AND. (MOD(var%dt_update,IFile%dt_file)) /= 0)) THEN
                WRITE(message_text,*) TRIM(var%name_var),': When update time step (',var%dt_update, &
                  ') is larger than file time step (',IFile%dt_file,'), it must be an integer multiplum of the file time step'
                CALL local_error('InputSetVarStartTime',TRIM(message_text(2:)))
              ENDIF
              tmp = var%dt_update / IFile%dt_file
              IF ((IFile%rec_step > 0) .AND. (IFile%rec_step /= tmp)) CALL local_error('InputSetVarStartTime', &
                                                    'File record progressing inconsistent with requirements of previous variable')
              IFile%rec_step = MAX(1,tmp)
          END SELECT
        ELSE ! No interpolation used
          IF (var%nIntrp == 0) var%nIntrp = 1
        ENDIF

        ! Allocate data buffers for interpolation data if needed
        IF ((var%nIntrp > 1) .AND. .NOT. ASSOCIATED(var%Intrp)) &
          var%intrp => InputDataIntrpAlloc(var%nIntrp,var%edims)
#ifdef HAVE_F2003
        ! If an interpolation function is explicitly given, it is assumed to match the nIntrp, nFuture, nAtOnce 
        ! settings and can be called automatically
        IF ((IAND(var%flags,INPUT_VAR_INTERPOL_USER) /= 0) .AND. ASSOCIATED(var%TimeIntrpFcn)) &
          var%flags = IAND(var%flags,INPUT_ALLFLAGS - INPUT_VAR_INTERPOL_USER)
#endif

        ! --- Determine file(s) and record(s) to read ---
        IF (IAND(IFile%flags,INPUT_FILE_TIME_ABSOLUTE) /= 0) THEN

          ! File time corresponding to model time and global record number
          tmf_mod    = model_time
          tmf_mod(1) = tmf_mod(1) + IFile%time_offset(1)
          ! Latest record time before or at this time
          nRec       = 1 + IFile%time_offset(4) ! Records are counted from first valid
          tmf_valid  = IFile%time_start(1:4)
          tmf_rec    = tmf_valid
          tmf_test   = ModVarTime(IFile,var%flags,tmf_valid)
          IF (CalendarTimeCmp(tmf_test,tmf_mod) >= 0) THEN
            nRec = nRec - IFile%rec_step
            DO WHILE (CalendarTimeCmp(tmf_test,tmf_mod) >= 0)
              tmf_rec  = tmf_valid
              CALL CalendarTimeAdd(tmf_valid,dt)
              tmf_test = ModVarTime(IFile,var%flags,tmf_valid)
              nRec = nRec + IFile%rec_step
            ENDDO
          ELSE
            DO WHILE (CalendarTimeCmp(tmf_test,tmf_mod) < 0)
              CALL CalendarTimeSub(tmf_rec,1_i8)
              CALL CalendarTimeSub(tmf_rec,dt)
              CALL CalendarTimeAdd(tmf_rec,1_i8)
              tmf_test = ModVarTime(IFile,var%flags,tmf_rec)
              nRec = nRec - IFile%rec_step
            ENDDO
          ENDIF
          tmf_valid = ModVarTime(IFile,var%flags,tmf_rec)

          ! Rewind to get first record needed
          IF ((var%nIntrp > 1) .AND. (IAND(var%stat,INPUT_VAR_BEFORE_FIRST) == 0)) THEN
            CALL CalendarTimeSub(tmf_rec,1_i8)
            CALL CalendarTimeSub(tmf_rec,INT(var%nIntrp-var%nFuture-1,i8)*dt)
            CALL CalendarTimeAdd(tmf_rec,1_i8)
            nRec = nRec - (var%nIntrp-var%nFuture-1)*IFile%rec_step
          ENDIF
          tmm_var = tmf_rec
          tmm_valid = tmf_rec

          ! Data before first valid data?
          IF (nRec < 1) THEN
            IF (IAND(IFile%flags,INPUT_FILE_NOCYCLE) == 0) THEN
              IF ((-nRec >= IFile%time_cycle(6)) .AND. (CalendarTimeCmp(IFile%time_cycle(1:4),IFile%time_start(1:4)) /= 0)) &
                CALL local_error('InputSetVarStartTime',TRIM(var%name_var)//                                                &
                  ': Requested data is more than one cycle before first available data but should not be entirely cycled.')
              IFile%stat = IOR(IFile%stat,INPUT_FILE_FIRST_CYCLE)
              nSub    = (-nRec / IFile%time_cycle(6)) * IFile%time_cycle(6)
              nRec    = nRec + nSub
              nSub    = IFile%time_cycle(6) + nRec ! - 1
              nRec    = IFile%time_enddata(6) - IFile%time_cycle(6) + nSub ! + 1
              tmf_rec = IFile%time_cycle(1:4)
              CALL CalendarTimeAdd(tmf_rec,INT(nSub,i8)*dt)
            ELSE
              var%stat  = IOR(var%stat,INPUT_VAR_BEFORE_FIRST)
              nRec      = IFile%time_start(5)
              tmf_rec   = IFile%time_start(1:4)
              tmm_valid = tmf_rec
            ENDIF
          ENDIF

          ! Data after last valid data?
          lEndReached = .FALSE.
          IF (nRec > IFile%time_enddata(6)) THEN
            IF (IAND(IFile%flags,INPUT_FILE_NOCYCLE) == 0) THEN
              nSub = ((nRec - IFile%time_enddata(6))/IFile%time_cycle(6) + 1) * IFile%time_cycle(6)
              nRec = MODULO(nRec-nSub-1,IFile%time_cycle(6)) + 1
              IF (IFile%dt_file < 0) CALL CalendarTimeSub(tmf_rec,1_i8)
              CALL CalendarTimeSub(tmf_rec,INT(nSub,i8)*dt)
              IF (IFile%dt_file < 0) CALL CalendarTimeAdd(tmf_rec,1_i8)
              IF (CalendarTimeCmp(tmf_rec,IFile%time_start) > 0) CALL CalendarTimeAdd(tmf_rec,INT(IFile%time_cycle(6),i8)*dt)
            ELSE
              nRec = IFile%time_enddata(6)
              tmf_rec  = IFile%time_enddata(1:4)
              lEndReached = .TRUE.
            ENDIF
          ENDIF

          ! Model time for which the previous  - actually never read - record of the variable is valid
          tmm_valid(1) = tmm_valid(1) - IFile%time_offset(1)
          IF (IAND(var%flags,INPUT_VAR_VALID_MASK) /= INPUT_VAR_VALID_ACTUAL) &
            tmm_valid = ModVarTime(IFile,INPUT_VAR_VALID_START,tmm_valid)
          nSub = var%nIntrp-var%nFuture
          IF (IAND(var%stat,INPUT_VAR_BEFORE_FIRST) == 0) THEN
            CALL CalendarTimeSub(tmm_valid,1_i8)
            CALL CalendarTimeSub(tmm_valid,INT(nSub,i8)*dt)
            CALL CalendarTimeAdd(tmm_valid,1_i8)
          ENDIF
          tmm_test = tmm_valid
          IF (IAND(var%stat,INPUT_VAR_BEFORE_FIRST) == 0) CALL CalendarTimeAdd(tmm_test,dt)
          tmm_test = ModVarTime(IFile,var%flags,tmm_test)
          IF ((CalendarTimeCmp(tmm_test,model_time) < 0)) THEN
            CALL CalendarTimeSub(tmm_valid,1_i8)
            CALL CalendarTimeSub(tmm_valid,dt) ! Hack needed in unpredictable cases
            CALL CalendarTimeAdd(tmm_valid,1_i8)
          ENDIF

          ! Get the file and file record corresponding to variable time
          tmf_file = InputFileTimeCurrGet(IFile,tmf_rec)
          CALL InputFileOpen(IFile,tmf_file(1),tmf_file(2),tmf_file(3),lOrderChg=lOC) ! Right file to avoid dest. of act rec #
          lOrderChg = lOrderChg .OR. lOC
#ifdef HAVE_F2003
          IF (lio) tmf_var = IFile%ft%FcnGetRecTime(IFile,1,.FALSE.)
#else
          IF (lio) tmf_var = NCFileGetRecTime(IFile,1,.FALSE.)
#endif
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
          CALL p_bcast(tmf_var,io_pe)
#endif

          IF (lEndReached) THEN
            IF (IFile%nRec < 1) Tmp = InputFileTimeStepGet(IFile,.FALSE.)
            IFile%curr_rec = IFile%nRec
            CALL CalendarTimeSub(tmm_valid,dt)
          ELSE
            IFile%curr_rec = 1
            tmf_var = ModVarTime(IFile,INPUT_VAR_VALID_START,tmf_var)
            tmf_rec = ModVarTime(IFile,INPUT_VAR_VALID_START,tmf_rec)
            DO WHILE (CalendarTimeCmp(tmf_var,tmf_rec) > 0)
              IFile%curr_rec = IFile%curr_rec + IFile%rec_step
              CALL CalendarTimeAdd(tmf_var,dt)
            ENDDO
          ENDIF
          IF (IFile%curr_rec > IFile%nRec) THEN
            WRITE(message_text,*) IFile%curr_rec,' of ',IFile%nRec
            CALL RmDblSpc(message_text)
            CALL local_error('InputSetVarStartTime',TRIM(var%name_var)//': Internal error in time system. '// &
              'Trying to access non-existing record'//TRIM(message_text)//' in file: '//TRIM(IFile%name_file))
          ENDIF

          ! Future fraction if current time does not match an update step
          IF (var%nIntrp > 1) THEN
            tmm_test = exp_start
            tmm_var  = tmm_test
            DO WHILE (CalendarTimeCmp(tmm_test,model_time) > 0)
              tmm_var = tmm_test
              CALL CalendarTimeAdd(tmm_test,INT(var%dt_update,i8))
            ENDDO
            IF (CalendarTimeCmp(tmm_test,model_time) /= 0) THEN
              CALL InputDataNew(dta)
              ALLOCATE(dta%dta(1,1,1,1,1,1,1))
              dta%flags    =  INPUT_DATA_FUTURE_FRAC
              dta%time     =  tmm_var
              dta%next     => var%aux_data
              dta%dta      =  0._dp
              var%aux_data => dta
            ENDIF
          ENDIF
          
        ELSE ! No absolute time
          IFile%curr_rec = 1
          tmf_valid    = model_time
          tmf_valid(1) = tmf_valid(1) + IFile%time_offset(1) 
          tmf_file = InputFileTimeCurrGet(IFile,tmf_valid)
          CALL InputFileOpen(IFile,tmf_file(1),tmf_file(2),tmf_file(3),lOrderChg=lOC) ! Right file to avoid dest. of act rec #
          tmm_valid = model_time
          CALL CalendarTimeSub(tmm_valid,1_i8)
          CALL CalendarTimeSub(tmm_valid,INT(var%nIntrp-var%nFuture,i8)*dt)
          CALL CalendarTimeAdd(tmm_valid,1_i8)
        ENDIF

        var%time = ModVarTime(IFile,var%flags,tmm_valid)

        ! Time until next update
        tmm_test = exp_start
        DO WHILE (CalendarTimeCmp(tmm_test,model_time) >= 0)
          CALL CalendarTimeAdd(tmm_test,INT(var%dt_update,i8))
        ENDDO
        dt = CalendarTimeDiff(model_time,tmm_test)
        var%nupdate(2) = INT(dt/INT(model_dt,i8))

      ENDIF ! parent_file exist and variable belongs to current sub-model

      var => var%next
      IF (lOrderChg) var => AllInputVars ! If the variable order have changed, redo all to be sure that all are processed
    ENDDO

  END SUBROUTINE InputSetVarStartTime

! Forwards multi-file to next record
  SUBROUTINE InputFileNextRec(IFile)
  TYPE (input_file_list), POINTER :: IFile

  INTEGER :: time_ret(5)
  LOGICAL :: lClose

    IFile%nRead = 0
    IFile%curr_rec = IFile%curr_rec + IFile%rec_step
    IF ((IFile%curr_rec > IFile%nRec) .OR. ((IFile%time_enddata(5) > 0) .AND. (IFile%curr_rec >= IFile%time_enddata(5)))) THEN
!      time_ret(1:4) = IFile%time_open
!      time_ret(1) = time_ret(1) + IFile%time_offset(1)
#ifdef HAVE_F2003
      IF (lio) time_ret(1:4) = IFile%ft%FcnGetRecTime(IFile,1,.FALSE.) ! Time of first record in found file
#else
      IF (lio) time_ret(1:4) = NCFileGetRecTime(IFile,1,.FALSE.)
#endif
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
      CALL p_bcast(time_ret,io_pe)
#endif
      CALL CalendarTimeAdd(time_ret(1:4),INT(IFile%dt_file,i8)*INT(IFile%curr_rec-1,i8))
      IFile%curr_rec = IFile%curr_rec - IFile%nRec
      IF (CalendarTimeCmp(time_ret(1:4),IFile%time_enddata(1:4)) <= 0) THEN ! Time to cycle?
        lClose = IAND(IFile%flags,INPUT_FILE_MULTI) /= 0
        IF (IAND(IFile%flags,INPUT_FILE_NOCYCLE) == 0) THEN
          IF (IAND(IFile%stat,INPUT_FILE_FIRST_CYCLE) /= 0) THEN
            IFile%stat = IAND(IFile%stat,INPUT_ALLFLAGS - INPUT_FILE_FIRST_CYCLE)
            time_ret = IFile%time_start(1:5)
          ELSE
            time_ret = IFile%time_cycle(1:5)
          ENDIF
          IF (CalendarTimeCmp(IFile%time_open,time_ret(1:4)) /= 0) THEN
            IF (lClose) CALL InputFileClose(IFile) ! Close file to secure open of new one at next read
            IFile%time_file = time_ret(1:4)
          ENDIF
          IFile%curr_rec = time_ret(5)
        ELSE ! No more data are available
          CALL InputFileClose(IFile)
          IFile%stat = IOR(IFile%stat,INPUT_FILE_END_REACHED)
        ENDIF
      ELSE ! Nop - just close - that corresponds to progressing to next file
        CALL InputFileClose(IFile)
      ENDIF
    ENDIF

  END SUBROUTINE InputFileNextRec

! Open a file for reading of variables and check consistency
  RECURSIVE SUBROUTINE InputFileOpen(IFile,year,month,day,dt_file,check,lOrderChg)
  TYPE (input_file_list), POINTER :: IFile
  INTEGER, OPTIONAL  , INTENT(IN) :: year, month, day
  INTEGER, OPTIONAL  , INTENT(IN) :: dt_file, check
  LOGICAL, OPTIONAL  ,INTENT(OUT) :: lOrderChg

  TYPE (input_file_list), POINTER :: file_save
  TYPE (input_var_list), POINTER :: curr_var
  TYPE (input_opt_list), POINTER :: opt
  REAL(dp), POINTER :: tst(:,:,:,:,:,:,:)
  LOGICAL :: lClose
  CHARACTER (len=128) :: name_new
  INTEGER :: rcheck, tim(4), Tmp

    IF (PRESENT(lOrderChg)) lOrderChg = .FALSE.
    rcheck = INPUT_FILE_CHECK_MASK ! Default: Do all checks
    IF (PRESENT(check)) rcheck = IAND(check,INPUT_FILE_CHECK_MASK)

    IF (IAND(IFile%stat,INPUT_FILE_EXTERN_APPLIED) == 0) THEN
      CALL InputFileGetExternalSettings(IFile,opt,Tmp)
      CALL InputFileApplyExternalSettings(IFile,opt,Tmp)
    ENDIF

    tim(4) = 0
    IF (PRESENT(year) .AND. PRESENT(month) .AND. PRESENT(day)) THEN
      tim(1) = year
      tim(2) = month
      tim(3) = day
    ELSE
      tim(1:3) = IFile%time_file(1:3)
    ENDIF
    WHERE (INT(tim,i8)+INT(HUGE(rcheck),i8) == 0)
      tim = model_time
    ENDWHERE
    tim = InputFileMultiInstanceGet(IFile,BaseTime=tim)
    IF (ANY(tim==-HUGE(Tmp))) RETURN

    lClose = .FALSE.
    IFile%stat = IAND(IFile%stat,INPUT_ALLFLAGS-INPUT_FILE_CYCLED)
    IF (liodata) THEN
      name_new = GetFileName(IFile,tim(1),tim(2),tim(3))

      ! Open file
      IF (IAND(IFile%stat,INPUT_FILE_OPEN) /= 0) THEN
        IF (ALL(tim(1:3)-IFile%time_open(1:3)==0) .OR. (IAND(IFile%flags,INPUT_FILE_MULTI) == 0)) THEN
          IF (IAND(IFile%stat,INPUT_FILE_CHECK_MASK) == rcheck) RETURN ! This file already opened and req. checks are done
          ! File open, but one or more checks are missing        
        ELSE
          lClose = .TRUE.
        ENDIF
      ENDIF
    ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
    CALL p_bcast(lClose,io_pe)
#endif
    IF (lClose) CALL InputFileClose(IFile) ! Another instance of a multi-file was open - close to get the right one
    IF (liodata) THEN
      IFile%name_file = TRIM(name_new)
#ifdef HAVE_F2003
      IFile%fid = IFile%ft%FcnOpen(IFile%name_file)
#else
      IFile%fid = NCFileOpen(IFile%name_file)
#endif
      IF (IFile%fid < 0) THEN
        IF (IAND(IFile%stat,INPUT_FILE_MULTI_OPENED) == 0) THEN
          ! If all variables in file have default values, it's ok if the file is not found.
          file_save => IFile ! We need to save the pointer for proper output after decoupling since IFile and curr_var%parent_file
                             ! may actually be the same object
          message_text = ''
          curr_var => AllInputVars
          DO WHILE (ASSOCIATED(curr_var))
            IF (ASSOCIATED(curr_var%parent_file,IFile)) THEN
              tst => InputDataGetSpecified(curr_var%aux_data,INPUT_DATA_DEFAULT,0)
              IF (.NOT. ASSOCIATED(tst) .OR. (IAND(curr_var%flags,INPUT_VAR_MISS_MASK) >= INPUT_VAR_MISS_ERR)) &
                CALL local_error('InputFileOpen','File: '//TRIM(IFile%name_file)// &
                                                 ' does not exist and defaults are not available for all file variables')
              ! Set variable to default value, detach variable from non-existing file and issue warning if requested
              NULLIFY(curr_var%parent_file)
              IF (ASSOCIATED(curr_var%dta)) THEN
                IF (SIZE(tst) == 1) THEN
                  curr_var%dta(:,:,:,:,:,:,:) = tst(1,1,1,1,1,1,1)
                ELSE
                  IF (.NOT. ArrayCopy(tst,curr_var%dta,(/SIZE(tst)/),SHAPE(curr_var%dta))) &
                    CALL local_error('InputFileOpen','Spread error')
                ENDIF
              ENDIF
              IF (IAND(curr_var%flags,INPUT_VAR_MISS_MASK) > INPUT_VAR_MISS_IGNORE) &
                message_text = TRIM(message_text)//', '//TRIM(curr_var%name_var)
            ENDIF
            curr_var => curr_var%next
          ENDDO
          IF (LEN_TRIM(message_text) > 0) CALL local_message('WARNING','File '//TRIM(file_save%name_file)// &
            ' does not exist. Using default value(s) for variable(s): '//TRIM(message_text(3:)))
          RETURN ! Do as if the file was opened
        ELSE
          IF (IFile%time_enddata(2) < 0) THEN
            WRITE(message_text,*) 'WARNING: Expected file: ',TRIM(IFile%name_file),' not found.'
            CALL local_message('InputFileOpen',message_text)
            IF (IAND(IFile%flags,INPUT_FILE_NOCYCLE) /= 0) THEN
              WRITE(message_text,*) 'WARNING: Keeping all variables from this file constant from now on'
            ELSE
              WRITE(message_text,*) 'WARNING: (Unintended?) cycling of multi-file performed'
            ENDIF
            IF (ldbg .AND. (IAND(mo_debug,INPUT_DBG_MASK) /= 0)) CALL local_message('InputFileOpen',message_text)
          ENDIF
        ENDIF
      ELSE
        IF (IAND(IFile%flags,INPUT_FILE_MULTI) /= 0) IFile%stat = IOR(IFile%stat,INPUT_FILE_MULTI_OPENED)
      ENDIF
    ENDIF !liodata
    IF ((IAND(mo_debug,INPUT_DBG_MASK)>=INPUT_DBG_OPENCLOSE) .AND. ldbg .AND. liodata &
        .AND. (IAND(IFile%stat,INPUT_FILE_OPEN)==0)) THEN
      WRITE(message_text,'(3a)') &
        TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Open: ',TRIM(IFile%name_file)
      CALL local_message('',message_text)
    ENDIF
    IFile%stat = IOR(IFile%stat,INPUT_FILE_OPEN)
    IF (IAND(IFile%flags,INPUT_FILE_MULTI) /= 0) IFile%stat = IOR(IFile%stat,INPUT_FILE_MULTI_OPENED)
    IFile%time_open = tim
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
    CALL p_bcast(IFile%flags    ,io_pe)
    CALL p_bcast(IFile%stat     ,io_pe)
    CALL p_bcast(IFile%time_open,io_pe)
    CALL p_bcast(IFile%name_file,io_pe)
#endif

    ! Set time step if needed and provided
    IF ((IFile%dt_file == 0) .OR. (IFile%dt_file == UNLIKELY_VAL)) THEN
      IF (PRESENT(dt_file)) THEN
        IFile%dt_file = dt_file
        IFile%stat = IOR(IFile%stat,INPUT_FILE_TIME_GOT)
      ENDIF
    ENDIF

    ! Perform consistency checks on file contents and prepare data handling chains
    CALL InputFileCheck(IFile,rcheck,lOrderChg)

    ! Progress file "pointer" to next file to open
    IF ((IAND(IFile%flags,INPUT_FILE_MULTI) /= 0) .AND. (IFile%dt_file /= UNLIKELY_VAL)) THEN
      IFile%time_file = IFile%time_open
      CALL CalendarTimeAdd(IFile%time_file,INT(IFile%dt_file*IFile%nRec,i8))
      IF (CalendarTimeCmp(IFile%time_file,IFile%time_enddata(1:4)) < 0) THEN
        IF (IAND(IFile%flags,INPUT_FILE_NOCYCLE) == 0) THEN
          IF (ABS(IFile%time_cycle(1)) <= HUGE(i8)-10) THEN
            IFile%time_file = InputFileTimeCurrGet(IFile,IFile%time_cycle(1:4))
          ELSE
            IFile%time_file = -HUGE(i8)
          ENDIF
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE InputFileOpen

! Closes an input file
  SUBROUTINE InputFileClose(IFile)
  TYPE (input_file_list), POINTER :: IFile

    IF (IAND(IFile%stat,INPUT_FILE_OPEN) /= 0) THEN
      IF (liodata) THEN
        IF ((IAND(mo_debug,INPUT_DBG_MASK) >= INPUT_DBG_OPENCLOSE) .AND. ldbg) &
          CALL local_message('',TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK)))//': Done: '//TRIM(IFile%name_file))
#ifdef HAVE_F2003
        IF (.NOT. IFile%ft%FcnClose(IFile%fid)) CALL local_error('InputFileClose',message_text)
#else        
        IF (.NOT. NCFileClose(IFile%fid)) CALL local_error('InputFileClose',message_text)
#endif
        IFile%fid = -1 ! Assume this is an invalid nc-identifier
      ENDIF
      IFile%stat = IAND(IFile%stat,INPUT_ALLFLAGS- &
        (INPUT_FILE_OPEN+INPUT_FILE_CHECK_MASK+INPUT_FILE_DIMS_GOT+INPUT_FILE_VARS_GOT+INPUT_FILE_TIME_GOT))
      IFile%name_file = ''
    ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
    CALL p_bcast(IFile%flags,io_pe)
    CALL p_bcast(IFile%stat,io_pe)
#endif

  END SUBROUTINE InputFileClose

! Closes all initital files and removes them and their associated variables from the list of files to process, same for init vars.
  SUBROUTINE InputDoneInit()
  TYPE (input_file_list), POINTER :: curr_file, prev_file
  TYPE (input_var_list),  POINTER :: curr_var,  prev_var
  LOGICAL lKeep

    ! Remove all variables read from initial files/obtained only in initialization from processing list
    NULLIFY(prev_var)
    curr_var => AllInputVars
    DO WHILE (ASSOCIATED(curr_var))
      IF (curr_var%sub_model == CurrSubModel) THEN
        lKeep = (curr_var%dt_update /= 0) .AND. (IAND(curr_var%flags,INPUT_VAR_MODEL) == 0)
        IF (.NOT. lKeep .AND. ASSOCIATED(curr_var%parent_file)) &
          curr_var%parent_file%nVar = curr_var%parent_file%nVar - 1
        IF (lKeep) THEN
          prev_var => curr_var
        ELSE
          CALL RemoveDependency(curr_var)
          IF (ASSOCIATED(prev_var)) THEN
            prev_var%next => curr_var%next
          ELSE
            AllInputVars => curr_var%next
          ENDIF
        ENDIF
      ENDIF
      curr_var => curr_var%next
    ENDDO

    ! Remove all initial files from processing list
    NULLIFY(prev_file)
    curr_file => AllInputFiles
    DO WHILE (ASSOCIATED(curr_file))
      IF ((curr_file%sub_model == CurrSubModel) .AND. &
          ((IAND(curr_file%flags,INPUT_FILE_INITIAL) /= 0) .OR. (curr_file%nVar <= 0)))THEN
        CALL InputFileClose(curr_file)
        IF (ASSOCIATED(prev_file)) THEN
          prev_file%next => curr_file%next
        ELSE
          AllInputFiles  => curr_file%next
        ENDIF
        curr_file => curr_file%next
      ELSE
        prev_file => curr_file
        curr_file => curr_file%next
      ENDIF
    ENDDO

  END SUBROUTINE InputDoneInit

! Closes all inputs
  SUBROUTINE InputClose()
  TYPE (input_file_list), POINTER :: curr

    IF (.NOT. moInputInitialized) RETURN
    curr => AllInputFiles
    DO WHILE (ASSOCIATED(curr))
      CALL InputFileClose(curr)
      curr => curr%next
    ENDDO
    IF (LogOpened) CLOSE(GetLogUnit())
    LogOpened = .FALSE.

  END SUBROUTINE InputClose

! ----------------------------------------------------------
! Main routines to obtain data from files
! ----------------------------------------------------------

! Read a new record from a file
  SUBROUTINE InputFileGetNewData(IFile,var)
  TYPE (input_file_list), POINTER :: IFile
  TYPE (input_var_list) , POINTER :: var

  TYPE (input_curr)               :: curr
#ifdef HAVE_F2003
  TYPE (input_fcn_list),  POINTER :: curr_fcn
#endif
  REAL(dp), POINTER :: dta(:,:,:,:,:,:,:)
  INTEGER :: tmp_time(4), mul, i
  INTEGER(i8) :: dt
  LOGICAL :: lCont, lMsg

    ! Close file if all data have been read or recycle single file
    IF (ASSOCIATED(IFile)) THEN
      IF ((IFile%curr_rec > IFile%nRec) .AND. (IAND(IFile%flags,INPUT_FILE_INITIAL) == 0) .AND. (var%dt_update /= 0)) THEN
        IF (IAND(IFile%flags,INPUT_FILE_MULTI) /= 0) THEN
          CALL InputFileClose(IFile)
        ELSE
          IFile%nRead = 0
        ENDIF
      ENDIF

      ! Open new file if necessary
      IF (IAND(IFile%stat,INPUT_FILE_OPEN+INPUT_FILE_END_REACHED) == 0) &
        CALL InputFileOpen(IFile)
      IF (IFile%curr_rec == 0) IFile%curr_rec = MAX(1,IFile%nRec) ! In some cases nRec is not set at this stage
      ! Test if file match with eventual group file
      IF (ASSOCIATED(var%group)) THEN
        IF (.NOT. ASSOCIATED(var%group%group_file,IFile)) THEN
          IF ((var%group%group_file%dt_file  /= IFile%dt_file)  .OR. &
              (var%group%group_file%curr_rec /= IFile%curr_rec) .OR. &
              (var%group%group_file%nRec     /= IFile%nRec)     .OR. &
              (IAND(var%group%group_file%flags,INPUT_FILE_MULTI) /= IAND(IFile%flags,INPUT_FILE_MULTI))) THEN
            IF (IAND(mo_debug,INPUT_DBG_MASK) /= 0) &
              CALL local_message('WARNING',TRIM(var%name_var)//' excluded from group '//TRIM(var%group%name_group)// &
                                 ' due to mismatch of file time parameters.')
            NULLIFY(var%group)
          ENDIF
        ENDIF
      ENDIF
    ENDIF

    IF (IAND(var%stat,INPUT_VAR_EXIST) /= 0) THEN ! Variable exist (in file)
      ! Update time for non-initial fields
      lCont = var%dt_update /= 0
      dt = INT(var%dt_update,i8)
      IF (ASSOCIATED(IFile) .AND. lCont) dt = INT(IFile%dt_file,i8)*INT(IFile%rec_step,i8)

      IF (lCont) THEN
        IF (ASSOCIATED(var%parent_file) .AND. (IAND(var%stat,INPUT_VAR_HAS_UNLIM)==0)) CALL local_error('InputFileGetNewData', &
          'Attempt to change time of variable '//TRIM(var%name_var)//' which does not contain a time dimension')
        ! Adding the time double and determining the half makes sure to process months of different length correctly
        tmp_time = var%time
        SELECT CASE (IAND(var%flags,INPUT_VAR_VALID_MASK))
          CASE (INPUT_VAR_VALID_ACTUAL,INPUT_VAR_VALID_START,INPUT_VAR_VALID_END,INPUT_VAR_VALID_MIDDAY)
            mul = 1
          CASE (INPUT_VAR_VALID_BEFORE,INPUT_VAR_VALID_MIDPOINT,INPUT_VAR_VALID_MIDMONTH,INPUT_VAR_VALID_MIDYEAR) 
            mul = 2
        END SELECT
        CALL CalendarTimeAdd(tmp_time,dt*INT(mul,i8))
        dt = CalendarTimeDiff(var%time,tmp_time)/INT(mul,i8)
        tmp_time = var%time
        CALL CalendarTimeAdd(tmp_time,dt)
        ! Turn the interpolation wheel to make room for new data
        IF (ASSOCIATED(var%Intrp)) THEN
          var%Intrp => var%Intrp%next
          var%Intrp%time = tmp_time
        ENDIF
        var%time = tmp_time
      ELSE
        IF (IAND(var%stat,INPUT_VAR_HAS_UNLIM) /= 0) THEN ! Use right record for initial variables
          IFile%curr_rec = IFile%init_rec
          IF ((IFile%curr_rec == UNLIKELY_VAL) .OR. (IFile%curr_rec < 1) .OR. (IFile%curr_rec > IFile%nRec)) &
            IFile%curr_rec = IFile%nRec
        ELSE
          IFile%curr_rec = 1 ! For safety when the variable do not contain the time dimension
        ENDIF
      ENDIF

      ! Perform data action chain
      curr%mask => var%masks
      curr%buf  => var%buffers
      curr%acp  =  1
#ifdef HAVE_F2003
      curr_fcn => var%DataFcns
      lCont = ASSOCIATED(curr_fcn)
      DO WHILE (lCont)
        CALL curr_fcn%fcn(INPUT_PROC_DO_DATA_ACTION,var,curr)
        curr_fcn => curr_fcn%next
        lCont = ASSOCIATED(curr_fcn)
        IF (lCont) lCont = (curr%acp <= SIZE(var%saction))
        IF (lCont) lCont = (var%saction(curr%acp) /= INPUT_VAR_DO_INTERPOLATE)
      ENDDO
#else
      CALL InputDoDataActions(var,curr,.TRUE.)
#endif
      lMsg = ldbg
      IF (lMsg) THEN
        lMsg = .FALSE.
        i = IAND(mo_debug,INPUT_DBG_MASK)
        IF ( (i == INPUT_DBG_READ_VAR) .OR. (i == INPUT_DBG_UPDATE_VAR) .OR. (i == INPUT_DBG_UPDATE_SUMMARY) .OR. &
            .NOT. ASSOCIATED(var%group) .AND.                                & 
            ((i == INPUT_DBG_READ_GROUP) .OR. (i == INPUT_DBG_UPDATE_GROUP))) lMsg = .TRUE.
      ENDIF
      IF (lMsg .AND. IAND(mo_debug,INPUT_DBG_AVG) /= 0) THEN
        IF (ASSOCIATED(var%intrp)) THEN
          dta => var%intrp%dta
        ELSE
          dta => var%dta
        ENDIF
        IF (SIZE(dta)>0) THEN
          WRITE(message_text(LEN_TRIM(message_text)+1:),'(a,g10.5)') ', avg = ',SUM(dta)/SIZE(dta)
        ELSE
          message_text = TRIM(message_text)//', avg = NaN'
        ENDIF
      ENDIF
      IF (lMsg) THEN
        lCont = .FALSE.
        i = INDEX(message_text,';')
        IF (i>0) message_text(i:i) = ':' ! Non-file vars don't have any text until here!
        IF (ASSOCIATED(IFile)) THEN
          IF (IAND(IFile%flags,INPUT_FILE_INITIAL) == 0) lCont = .TRUE.
        ENDIF
        IF (lCont) THEN
          IF (var%dt_update /= 0) THEN
            CALL local_message('',message_text(1:i+1)//'Valid: '//TRIM(DateString(var%time,IAND(mo_debug,INPUT_DBG_TIME_MASK)))//&
                                  ': '//TRIM(message_text(i+2:)))
          ELSE
            CALL local_message('',TRIM(message_text))
          ENDIF
        ELSE
          CALL local_message('',TRIM(message_text))
        ENDIF
      ENDIF
      message_text = '' ! Text otherwise accumulates on non-ldbg PEs

      ! Update variable status
      var%nRead = var%nRead + 1
      var%stat = IOR(var%stat,INPUT_VAR_READ)
      var%nLastRead = 0
      IF (.NOT. ASSOCIATED(var%Intrp)) &
        var%stat = IOR(var%stat,INPUT_VAR_UPDATED) ! The only way to update non-interpolated vars are to read them

    ENDIF ! Var exist

    ! File status update
    IF (ASSOCIATED(IFile)) THEN
      IFile%stat = IOR(IFile%stat,INPUT_FILE_READ)
      IF (var%dt_update /= 0) IFile%nRead = IFile%nRead + 1 ! Initial variables don't count for further processing
      ! Close file if all data read and not single file cycled
      IF ((IFile%nRead >= IFile%nVar) .AND. (IAND(IFile%flags,INPUT_FILE_INITIAL) == 0)) CALL InputFileNextRec(IFile)
    ENDIF

  END SUBROUTINE InputFileGetNewData

! Updates linear interpolation of all variables in a list when necessary
  SUBROUTINE InputFileUpdateVars(lFirstStep)
  LOGICAL, INTENT(IN) :: lFirstStep

  TYPE (input_file_list), POINTER :: curr_file
  TYPE (input_var_list) , POINTER :: var
  TYPE (input_data_list), POINTER :: dta
  TYPE (input_curr)               :: curr
  CHARACTER (len=256)             :: msg
  INTEGER :: tt(4), nRead, i, nFile, nUpdate
  INTEGER(i8) :: dt1, dt2
  LOGICAL lCont, lFirst, lTmp, lMsg
#ifdef HAVE_F2003
  TYPE (input_fcn_list),  POINTER :: curr_fcn
#endif

    ! Update file flags and debug info
    curr_file => AllInputFiles
    DO WHILE (ASSOCIATED(curr_file))
      IF (curr_file%sub_model == CurrSubModel) THEN
        curr_file%stat      =  IAND(curr_file%stat,INPUT_ALLFLAGS - INPUT_FILE_READ - INPUT_FILE_CYCLED)
        curr_file%nRead_dbg =  0
      ENDIF
      curr_file => curr_file%next
    ENDDO
    CurrGroup => AllInputGroups
    DO WHILE (ASSOCIATED(CurrGroup))
      CurrGroup%flags = IAND(CurrGroup%flags,INPUT_ALLFLAGS - INPUT_GROUP_MSG_DONE)
      CurrGroup => CurrGroup%next
    ENDDO
    NULLIFY(CurrGroup)

    ! Output group messages?
    lMsg = ldbg
    IF (lMsg) THEN
      i = IAND(mo_debug,INPUT_DBG_MASK)
      IF ((i /= INPUT_DBG_READ_GROUP) .AND. (i /= INPUT_DBG_UPDATE_GROUP)) lMsg = .FALSE.
    ENDIF

    ! Main processing
    nRead  = -1 
    lFirst = .TRUE.
    lCont  = .TRUE.
    DO WHILE (lCont) ! Loop as long as more time steps of any variable should be read
      lCont = .FALSE.
      var => AllInputVars
      DO WHILE (ASSOCIATED(var))
        IF (var%sub_model == CurrSubModel) THEN
          ! Reading of variable
          !   Should variable be read? (No reason to read if no update should be performed)
          IF (lFirst) THEN
            IF (IAND(var%flags,INPUT_VAR_AUTORESET_OFF) == 0) &
              var%stat = IAND(var%stat,INPUT_ALLFLAGS - INPUT_VAR_READ - INPUT_VAR_UPDATED) ! Reset read and update status
            var%nLastRead   = var%nLastRead   + 1
            var%nLastUpdate = var%nLastUpdate + 1
            var%nToRead = 0
            IF (var%time(1)==-HUGE(i)) THEN
              var%time = model_time
              CALL CalendarTimeSub(var%time,INT(var%parent_file%dt_file,i8))
            ENDIF
            i = CalendarTimeCmp(var%time,model_time)
            IF ((i <= 0) .AND. lFirstStep .AND. (IAND(var%stat,INPUT_VAR_BEFORE_FIRST) /= 0)) i = 1
            IF ((i >  0) .AND. (var%next_update <= 0)) THEN
              var%nToRead = var%nAtOnce
              IF (lFirstStep) THEN
                var%nToRead = var%nIntrp ! Fill all buffers from the start
                IF (IAND(var%stat,INPUT_VAR_BEFORE_FIRST) /= 0) THEN
                  IF (IAND(var%flags,INPUT_VAR_INTERPOL_USER) == 0) var%nToRead = 1 ! 1st data record is used until it is "expired"
                  var%stat = IAND(var%stat,INPUT_ALLFLAGS - INPUT_VAR_BEFORE_FIRST)
                ENDIF
                tt = var%time
                IF (ASSOCIATED(var%parent_file)) THEN
                  CALL CalendarTimeAdd(tt,INT(var%parent_file%dt_file,i8))
                ELSE
                  CALL CalendarTimeAdd(tt,INT(var%dt_update,i8))
                ENDIF
                IF (CalendarTimeDiff(tt,model_time) == 0) var%nToRead = 1
              ENDIF
            ENDIF
          ENDIF
          IF (var%nToRead > 0) THEN
            lTmp = .TRUE.
            IF (ASSOCIATED(var%parent_file)) THEN
              IF (IAND(var%parent_file%stat,INPUT_FILE_END_REACHED) == 0) THEN
                lCont = .TRUE.
              ELSE
                lTmp = .FALSE.
              ENDIF
            ELSE
              lCont = .TRUE.
            ENDIF
            IF (lTmp) THEN
              CALL InputFileGetNewData(var%parent_file,var)
              var%nLastUpdate = 0
            ENDIF
            var%nToRead = var%nToRead - 1
          ENDIF
          IF (var%nToRead > 0) lCont = .TRUE.
        ENDIF
        var => var%next
      ENDDO
      lFirst = .FALSE.
      IF (ASSOCIATED(CurrGroup)) CALL GroupReadMsg(CurrGroup)
      NULLIFY(CurrGroup)
    ENDDO

    ! Update of model variable
    IF (IAND(mo_debug,INPUT_DBG_MASK) /= INPUT_DBG_UPDATE_GROUP) lMsg = .FALSE.
    NULLIFY(CurrGroup)
    nUpdate = 0
    var => AllInputVars
    DO WHILE (ASSOCIATED(var))
      IF (var%sub_model == CurrSubModel) THEN
        ! Should variable be updated?
        IF (var%next_update <= 0) THEN

          ! All data read? If yes: keep at same level from now on
          IF (ASSOCIATED(var%parent_file)) THEN
            IF ((IAND(var%parent_file%stat,INPUT_FILE_END_REACHED) /= 0) .AND. (CalendarTimeCmp(var%time,model_time)) >= 0) THEN
              var%next_update = HUGE(i)
              var%nupdate     = HUGE(i)
              var%dt_update   = HUGE(i)
            ENDIF
          ENDIF

          IF (ASSOCIATED(var%intrp) .AND. ((var%nRead > 1) .OR. lFirstStep)) THEN
            IF (IAND(var%flags,INPUT_VAR_INTERPOL_USER) == 0) THEN

              ! Perform data action chain from the point of interpolation
              curr%buf  => var%buffers_intrp
              curr%mask => var%masks_intrp
              curr%acp  =  var%acp_intrp
#ifdef HAVE_F2003
              curr_fcn  => var%ProcInterpol
              DO WHILE (ASSOCIATED(curr_fcn))
                CALL curr_fcn%fcn(INPUT_PROC_DO_DATA_ACTION,var,curr)
                curr_fcn => curr_fcn%next
              ENDDO
#else
              CALL InputDoDataActions(var,curr,.FALSE.)
#endif
            ENDIF
            IF (lMsg .AND. ASSOCIATED(CurrGroup) .AND. .NOT. ASSOCIATED(var%group,CurrGroup)) THEN
              i = INDEX(message_text,': Interpolation time: ')
              IF (i > 0) THEN
                msg = TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK)))//': Data: '// &
                      TRIM(CurrGroup%name_group)//TRIM(message_text(i:))
              ELSE 
                WRITE(msg,'(4a,f7.4)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Data: ', &
                                       TRIM(CurrGroup%name_group),': Future frac = ',REAL(dt2,dp)/REAL(dt1,dp)
              ENDIF
              CALL local_message('',TRIM(msg))
            ENDIF
            CurrGroup => var%group ! Adapt to new group
            nUpdate = nUpdate + 1
            var%stat = IOR(var%stat,INPUT_VAR_UPDATED) ! Indicate that user should redo interpolation or that it has been done
            var%nLastUpdate = 0

            IF (ldbg .AND. ((IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_UPDATE_VAR) .OR. &
                (IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_UPDATE_GROUP))) THEN
              lTmp = IAND(var%flags,INPUT_VAR_INTERPOL_USER) == 0
              WRITE(message_text,'(3a)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Data: ', &
                                                TRIM(var%name_var)
              IF (lTmp .AND. (var%intrp%time(1) /= -HUGE(i8)) .AND. (var%intrp%prev%time(1) /= -HUGE(i8))) THEN
                dt1 = CalendarTimeDiff(var%intrp%prev%time,var%intrp%Time)
                dt2 = CalendarTimeDiff(var%intrp%prev%time,model_time)
                IF (ASSOCIATED(var%aux_data)) THEN
                  IF (IAND(var%aux_data%flags,INPUT_DATA_TYPE_MASK) == INPUT_DATA_FUTURE_FRAC) THEN
                    dt2 = CalendarTimeDiff(var%intrp%prev%time,var%aux_data%time)
                    WRITE(message_text(LEN_TRIM(message_text)+1:),'(2a)') &
                      ': Interpolation time: ',TRIM(InputTimeString(time=var%aux_data%time))
                  ENDIF
                ENDIF
                WRITE(message_text(LEN_TRIM(message_text)+1:),'(a,f7.4)') ': Future frac = ',REAL(dt2,dp)/REAL(dt1,dp)
              ENDIF
              IF (.NOT. ASSOCIATED(var%group) .OR. (IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_UPDATE_VAR)) THEN
                IF (IAND(mo_debug,INPUT_DBG_AVG) /= 0) THEN
                  IF (SIZE(var%dta) > 0) THEN
                    WRITE(message_text(LEN_TRIM(message_text)+1:),'(a,g10.5)') ', avg = ',SUM(var%dta)/SIZE(var%dta)
                  ELSE
                    message_text = TRIM(message_text)//', avg = NaN'
                  ENDIF
                ENDIF
                CALL RmDblSpc(message_text)
                CALL local_message('',TRIM(message_text))
              ENDIF
            ENDIF
          ENDIF

          ! Remove preset future fraction for future interpolations
          IF (ASSOCIATED(var%aux_data)) THEN
            IF (IAND(var%aux_data%flags,INPUT_DATA_TYPE_MASK) == INPUT_DATA_FUTURE_FRAC) THEN
              dta => var%aux_data
              var%aux_data => var%aux_data%next
              DEALLOCATE(dta) ! TODO: Also dta%dta should be deallocated?
            ENDIF
          ENDIF

          ! Set flags and prepare for next update
          IF ((var%dt_update < 0) .OR. (var%nupdate(1) < -1)) THEN
            IF (var%nupdate(2) == 0) THEN
              tt = model_time
              CALL CalendarTimeAdd(tt,INT(var%dt_update,i8))
              i = INT(CalendarTimeDiff(model_time,tt) / model_dt)
              IF (var%nupdate(1) < -1) THEN
                var%nupdate(1) = var%nupdate(1) + i
              ELSE
                var%nupdate(1) = i
              ENDIF
            ELSE
              var%nupdate(1) = var%nupdate(2)
              var%nupdate(2) = 0
            ENDIF
          ELSE
            var%nupdate(1) = var%nupdate(2)
            var%nupdate(2) = var%nupdate(3)
          ENDIF

          var%next_update = var%nupdate(1)

        ENDIF ! Interpolation update

        ! When to do next update
        var%next_update = var%next_update - 1

      ENDIF

      ! Process next variable
      var => var%next

    ENDDO ! ASSOCIATED(var)
    IF (lMsg .AND. ASSOCIATED(CurrGroup)) THEN
      i = INDEX(message_text,': Interpolation time: ')
      IF (i > 0) THEN
        msg = TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK)))//': Data: '// &
              TRIM(CurrGroup%name_group)//TRIM(message_text(i:))
      ELSE 
        WRITE(msg,'(4a,f7.4)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Data: ', &
                               TRIM(CurrGroup%name_group),': Future frac = ',REAL(dt2,dp)/REAL(dt1,dp)
      ENDIF
      CALL local_message('',TRIM(msg))
    ENDIF

    ! Summarized debug output
    IF (ldbg .AND. ((IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_READ_SUMMARY)   .OR. &
                    (IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_READ_FILE   ))) THEN
      nFile = 0
      nRead = 0
      curr_file => AllInputFiles
      DO WHILE (ASSOCIATED(curr_file))
        IF (curr_file%sub_model == CurrSubModel) THEN
          nRead = nRead + curr_file%nRead_dbg
          IF (curr_file%nRead_dbg > 0) nFile = nFile + 1
          IF ((curr_file%nRead_dbg /= 0) .AND. (IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_READ_FILE)) THEN
            WRITE(message_text,'(2a,i10,3a)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Read:', &
                                              curr_file%nRead_dbg,' variable(s) x record(s) (',TRIM(curr_file%name_file),')'
            CALL RmDblSpc(message_text)
            CALL local_message('',TRIM(message_text))
          ENDIF
        ENDIF
        curr_file => curr_file%next
      ENDDO
      IF (nRead > 0 .AND. IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_READ_SUMMARY) THEN
        WRITE(message_text,'(2a,2(i10,a))') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Read:', &
                                         nRead,' variable(s) x record(s) from ',nFile,' file(s)'
        CALL RmDblSpc(message_text)
        CALL local_message('',TRIM(message_text))
      ENDIF
    ENDIF

    IF (ldbg .AND. (IAND(mo_debug,INPUT_DBG_MASK) == INPUT_DBG_UPDATE_SUMMARY) .AND. (nUpdate > 0)) THEN
      WRITE(message_text,'(2a,i10,a)') TRIM(DateString(model_time,IAND(mo_debug,INPUT_DBG_TIME_MASK))),': Data: ', &
                                        nUpdate,' variable(s)'
      CALL RmDblSpc(message_text)
      CALL local_message('',TRIM(message_text))
    ENDIF

  END SUBROUTINE InputFileUpdateVars

! Reads initial variables which has not yet been read
  SUBROUTINE InputGetData(IFile,var,recno)
  TYPE (input_file_list), POINTER, OPTIONAL :: IFile
  TYPE (input_var_list),  POINTER, OPTIONAL :: var
  INTEGER,             INTENT(in), OPTIONAL :: recno

  TYPE (input_var_list), POINTER :: curr_var
  TYPE (input_opt_list), POINTER :: opt
  INTEGER :: cnt
  LOGICAL :: lCont

    IF (.NOT. moInputInitialized) CALL local_error('InputGetData','Please call InputInit before attempting to access data')
    ! Make sure that dependencies are properly resolved
    CALL SortDependencies()

    ! Go through variable list and read what should be read/obtained
    IF (PRESENT(var)) THEN
      curr_var => var
    ELSE
      curr_var => AllInputVars
    ENDIF
    DO WHILE(ASSOCIATED(curr_var))
      CALL InputVarGetExternalSettings(curr_var,opt,cnt)
      CALL InputVarApplyExternalSettings(curr_var,opt,cnt)
      lCont = .FALSE.
      IF (curr_var%sub_model == CurrSubModel) THEN
        IF (ASSOCIATED(curr_var%parent_file)) THEN
          IF ((curr_var%dt_update == 0) .AND. (IAND(curr_var%stat,INPUT_VAR_READ) == 0)) THEN ! Process only initial variables
            lCont = .TRUE.
            IF (PRESENT(IFile)) THEN 
              IF (.NOT. ASSOCIATED(curr_var%parent_file,IFile)) lCont = .FALSE.
            ENDIF
            IF (lCont) THEN
              lCont = .FALSE.
              IF ((IAND(curr_var%parent_file%stat,INPUT_FILE_OPEN) == 0) .OR.                            &
                   IAND(curr_var%parent_file%stat,INPUT_FILE_CHECK_INITIAL) /= INPUT_FILE_CHECK_INITIAL) &
                CALL InputFileOpen(curr_var%parent_file,check=INPUT_FILE_CHECK_INITIAL,lOrderChg=lCont)
              IF (PRESENT(recno) .AND. (curr_var%parent_file%init_rec == UNLIKELY_VAL)) curr_var%parent_file%init_rec = recno
              CALL InputFileGetNewData(curr_var%parent_file,curr_var)
              curr_var%stat = IOR(curr_var%stat,INPUT_VAR_READ)
            ENDIF
          ENDIF
        ELSE
          IF (IAND(curr_var%stat,INPUT_VAR_INITIALIZED) == 0) &
            lCont = InputCheckVarDims(curr_var%parent_file) ! File should be NULL
          IF (curr_var%dt_update == 0) &
            CALL InputFileGetNewData(curr_var%parent_file,curr_var) ! File should be NULL
        ENDIF
      ENDIF
      IF (lCont) THEN
        curr_var => AllInputVars
      ELSE
        curr_var => curr_var%next
      ENDIF
      IF (PRESENT(var)) NULLIFY(curr_var)
    ENDDO

  END SUBROUTINE InputGetData

! Update all inputs by closing, opening, reading and interpolating all necessary input data as requested
  SUBROUTINE InputUpdate
  LOGICAL, SAVE :: lFirst = .TRUE.

    IF (.NOT. moInputInitialized) RETURN
    IF (lFirst) THEN
      lInLoop = .TRUE.
      CALL InputGetData()         ! Secure to read all init-vars. Specifically vars which are init by file-association.
      CALL InputDoneInit()        ! Don't waste time scanning through initial files and their variables
      CALL InputSetVarStartTime() ! Determine time and relevant files for the first read of each variable
    ENDIF

    CALL InputFileUpdateVars(lFirst)

    ! Advance "model" time
    CALL CalendarTimeAdd(model_time,INT(model_dt,i8))

    lFirst = .FALSE.

  END SUBROUTINE InputUpdate

! ------------------------------------
!  Initializing functions/subroutines
! ------------------------------------

! Initialize module variables
  SUBROUTINE InputInit(dt_model,year,month,day,SOD,startyear,startmonth,startday,startSOD,Debug,sDbg,calendar,dt_unit,npe,pe,iope, &
                       lallread)
  INTEGER, INTENT(in) :: dt_model, year
  INTEGER, OPTIONAL, INTENT(in) :: month, day, SOD, startyear, startmonth, startday, startSOD, calendar, npe, pe, iope, Debug
  LOGICAL, OPTIONAL, INTENT(in) :: lallread
  CHARACTER (len=*), INTENT(in), OPTIONAL :: dt_unit,sDbg

  INTEGER :: i, mod_time(4), exp_time(4)
  INTEGER, POINTER :: tmp_arr(:,:)

    NULLIFY(AllInputFiles,AllInputVars,AllInputGrpVars,AllInputDims,AllInputEqDims,AllInputColDimDistrib,AllInputData, &
            AllInputOpts,AllInputMasks,AllInputGroups,CurrGroup,LastFile,DummyFile,DummyVar,src_pe)

    ! Parallel setup
    my_pe     = -1     ! Default is "not a parallel run"
    n_pe      =  1
    n_iope    =  1
    io_pe     =  0
    lio       = .TRUE.
    liodata   = .TRUE.
    lioserver = .FALSE.
    IF (PRESENT(npe)) THEN
      IF (PRESENT(iope)) io_pe = iope
      IF (.NOT.PRESENT(pe)) CALL local_error('InputInit','Local PE number not specified for parallel run')
      my_pe = pe
      n_pe  = npe
      lio   = my_pe==io_pe ! For non-parallel runs, master io should always be performed from this (i.e. the only) pe
      ALLOCATE(src_pe(n_pe))
      ! Default is that io_pe does all the reading and distributes
      src_pe(:) = io_pe 
      lioserver = lio
      liodata   = lio
      ldbg      = lio
      IF (PRESENT(lallread)) THEN
        IF (lallread) THEN
          DO i=1,n_pe
            src_pe(i) = i-1
          ENDDO
          lioserver = .FALSE.
          liodata   = .TRUE.
        ENDIF
      ENDIF
    ENDIF

    ! Set debug level
    mo_debug = INPUT_DBG_DEFAULT
    IF (PRESENT(Debug)) CALL InputDebugLevelSet(Debug,Flag=.TRUE.)
    IF (PRESENT(sDbg)) THEN 
      IF (LEN_TRIM(sDbg) > 0) CALL InputDebugLevelSet(sDebug=sDbg,Flag=.TRUE.)
    ENDIF
    IF (IAND(mo_debug,INPUT_DBG_DBG_IODATA+INPUT_DBG_DBG_ALL) /= 0) THEN
      WRITE(message_text,*) pe
      CALL local_message('InputInit','This file contains log output from PE:'//TRIM(message_text))
    ENDIF
 
    ! Time setup
    CALL InputTimeSet(year,month,day,SOD,dt_model,dt_unit,mod_time)
    IF (ANY(model_time /= UNLIKELY_VAL) .AND. ANY(model_time - mod_time /= 0)) CALL local_error('InputInit', &
      'Attempt to initialize sub-models with different run starting times')
    model_time = mod_time
    exp_time = mod_time
    IF (PRESENT(startyear))  exp_time = (/startyear,1,1,0/) 
    IF (PRESENT(startmonth)) exp_time(2) = startmonth
    IF (PRESENT(startday))   exp_time(3) = startday
    IF (PRESENT(startSOD))   exp_time(4) = startSOD
    IF (ANY(exp_start /= UNLIKELY_VAL) .AND. ANY(exp_start - exp_time /= 0)) CALL local_error('InputInit', &
      'Attempt to initialize sub-models with different experiment starting times')
    exp_start = exp_time
    ALLOCATE(tmp_arr(CurrSubModel,7)) ! Year, month, day, SOD, dt, mo_debug, log_unt
    tmp_arr(:,1:4) = SPREAD(model_time,DIM=1,NCOPIES=CurrSubModel)
    IF (ASSOCIATED(submodel_settings)) THEN
      tmp_arr(1:CurrSubModel-1,5:7) = submodel_settings(:,5:7)
      DEALLOCATE(submodel_settings)
    ENDIF
    tmp_arr(CurrSubModel,5:7) = (/model_dt,mo_debug,GetLogUnit()/)
    submodel_settings => tmp_arr

    ! Set calendar if given
    IF (PRESENT(calendar)) THEN
      IF (.NOT. CalendarSet(calendar)) CALL local_error('InitInput','Bad calendar specification')
    ELSE
      IF (.NOT. CalendarSet(INPUT_CALENDAR_GREGORIAN)) CALL local_error('InitInput','Internal bad calendar specification')
    ENDIF
    IF (lio) THEN
      i = CalendarGet(message_text)
      CALL local_message('InputInit','Calendar selected: '//TRIM(message_text))
    ENDIF

    ! File type support
#ifdef HAVE_F2003
    CALL NCRegister()       ! Register functions to handle NetCDF files
#endif

    moInputInitialized = .TRUE.

  END SUBROUTINE InputInit

END MODULE mo_input
