!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! This module contains procedures for various mo_input data objects and for time and file name management 
! Never directly use this module from outside the mo_input framework - only via mo_input.f90
MODULE mo_input_misc
USE mo_input_strings
USE mo_input_calendar
USE mo_input_types
IMPLICIT NONE

PUBLIC :: InputSubModelGetNew, InputSubModelGetCurr, InputSubModelSelect
PUBLIC :: InputFileSplitName, BuildGenericName, GetFileName
PUBLIC :: ModVarTime
PUBLIC :: InputTimeGet, InputTimeSet, InputTimeString
PUBLIC :: InputVarInterpolBoxCar, InputVarInterpolLinear, InputVarGetInterpolData
#ifdef HAVE_F2003
PUBLIC :: InputFcnNew, InputFcnAdd, InputFcnDone
#endif

CONTAINS

! ----------------------------------------------------
!  Creating, adding, setting, modifying, displaying:
!     Procedural references (Fortran 2003 only)
! ----------------------------------------------------

#ifdef HAVE_F2003
! Creates a new procedure ref element with default values
  SUBROUTINE InputFcnNew(fcn)
  TYPE (input_fcn_list), POINTER :: fcn

  TYPE (input_fcn_list), POINTER :: new

    ALLOCATE(new)
    NULLIFY(new%next,new%fcn)

    fcn => new

  END SUBROUTINE InputFcnNew

! Creates a new procedure ref element and adds it to a list
  FUNCTION InputFcnAdd(list)
  TYPE (input_fcn_list), POINTER :: InputFcnAdd
  TYPE (input_fcn_list), POINTER :: list

  TYPE (input_fcn_list), POINTER :: new, curr

    CALL InputFcnNew(new)
    InputFcnAdd => new
    IF (ASSOCIATED(list)) THEN
      curr => list
      DO WHILE (ASSOCIATED(curr%next))
        curr => curr%next
      ENDDO
      curr%next => new
    ENDIF

  END FUNCTION InputFcnAdd

! Destroys a list of procedure refs
  SUBROUTINE InputFcnDone(fcn)
  TYPE (input_fcn_list), POINTER :: fcn

  TYPE (input_fcn_list), POINTER :: curr, next

    curr => fcn
    DO WHILE (ASSOCIATED(curr))
      next => curr%next
      DEALLOCATE(curr)
      curr => next
    ENDDO

  END SUBROUTINE InputFcnDone
#endif

! --------------------------
!  Getting/setting submodel
! --------------------------

! Create new submodel and set this as the current
  INTEGER FUNCTION InputSubModelGetNew()
    NextSubModel = NextSubModel + 1
    CurrSubModel = NextSubModel
    InputSubModelGetNew = NextSubModel
  END FUNCTION InputSubModelGetNew

! Get number of current submodel
  INTEGER FUNCTION InputSubModelGetCurr()
    InputSubModelGetCurr = CurrSubModel
  END FUNCTION InputSubModelGetCurr

! Set current submodel
  SUBROUTINE InputSubModelSelect(SubModel)
  INTEGER, INTENT(in) :: SubModel

    IF ((SubModel < 1) .OR. (SubModel > NextSubModel)) CALL local_error('InputSubModelSelect','Unknown sub-model')
    IF (ASSOCIATED(submodel_settings)) THEN
      ! Save settings for the previous current sub-model
      submodel_settings(CurrSubModel,1:4) = model_time
      submodel_settings(CurrSubModel,6:7) = (/mo_debug,GetLogUnit()/) ! Ind 5, dt_model can only be specified using InputInit
      ! Extract settings for new actual sub-model
      model_time = submodel_settings(SubModel,1:4)
      model_dt   = submodel_settings(SubModel,5)
      mo_debug   = submodel_settings(SubModel,6)
      CALL SetLogUnit(submodel_settings(SubModel,7))
    ENDIF
    CurrSubModel = SubModel

  END SUBROUTINE InputSubModelSelect

! --------------------------
!  Working with file names
! --------------------------

! Split one-string filename info prefix, suffix and time flags - in case of comma sep multi-name, only the first is split
  SUBROUTINE InputFileSplitName(Str,Pre,Suf,Flags)
  CHARACTER (len=* ), INTENT(in ) :: Str
  CHARACTER (len=64), INTENT(out) :: Pre, Suf
  INTEGER,            INTENT(out) :: Flags

  CHARACTER (len=18) :: TimeAbb = '+-yYmMdDhHiIsS>>'
  INTEGER :: pos, cfp, ind, en

    en = INDEX(Str,',')
    IF (en <= 0) en = LEN_TRIM(Str)+1
    Flags = 0
    pos = INDEX(Str(1:en-1),'<')
    IF (pos < 1) THEN
      Pre = TRIM(Str(1:en-1))
      Suf = ''
      RETURN
    ENDIF
    Pre = Str(1:pos-1)
    pos = pos + 1
    cfp = 0
    ind = 1
    DO WHILE ((ind > 0) .AND. (ind < 8))
      ind = (INDEX(TimeAbb,Str(pos:pos))+1)/2
      IF (ind == 0) THEN
        WRITE(message_text,*) 'Invalid time macro in file name: ',TRIM(Str(1:en-1));
        CALL local_error('InputFileSplitName',message_text)
      ENDIF
      IF (ind <= cfp) THEN
        WRITE(message_text,*) 'Time macros must follow in order of decreasing intervals in file name: ',TRIM(Str(1:en-1));
        CALL local_error('InputFileSplitName',message_text)
      ENDIF
      IF (ind < 8) THEN
        Flags = IOR(Flags,ISHFT(INPUT_FILE_NAME_YEAR_SIGN,ind-1))
        IF (ind==1) pos = pos + 1 ! For signed years, skip over the y's instead of only + or -
        cfp = 1
        DO WHILE (Str(pos+1:pos+1)==Str(pos:pos))
          cfp = cfp + 1
          pos = pos + 1
        ENDDO
        IF (ind==2 .AND. cfp==4) Flags = IOR(Flags,INPUT_FILE_NAME_YEAR_SIGN) ! INPUT_FN_YEAR = IFN_YEAR_SIGN + IFN_YEAR_LONG
        pos = pos + 1
        cfp = ind
      ENDIF
    ENDDO
    Suf = Str(pos+1:en-1)

  END SUBROUTINE InputFileSplitName

! Builds the generic name (with time macros) of a file 
  FUNCTION BuildGenericName(IFile)
  CHARACTER (len=128) :: BuildGenericName
  TYPE (input_file_list), POINTER :: IFile

  CHARACTER (len=128) :: tmpstr

    tmpstr = TRIM(IFile%name_pre)
    IF (IAND(IFile%flags,INPUT_FILE_MULTI) /= 0) THEN
      tmpstr = TRIM(tmpstr)//'<'
      IF (IAND(IFile%flags,INPUT_FILE_NAME_YEAR) == INPUT_FILE_NAME_YEAR_SIGN) tmpstr = TRIM(tmpstr)//'+yyyyyy'
      IF (IAND(IFile%flags,INPUT_FILE_NAME_YEAR) == INPUT_FILE_NAME_YEAR_LONG) tmpstr = TRIM(tmpstr)//'yyyyyy'
      IF (IAND(IFile%flags,INPUT_FILE_NAME_YEAR) == INPUT_FILE_NAME_YEAR     ) tmpstr = TRIM(tmpstr)//'yyyy'
      IF (IAND(IFile%flags,INPUT_FILE_NAME_MONTH ) /= 0) tmpstr = TRIM(tmpstr)//'mm'
      IF (IAND(IFile%flags,INPUT_FILE_NAME_DAY   ) /= 0) tmpstr = TRIM(tmpstr)//'dd'
      IF (IAND(IFile%flags,INPUT_FILE_NAME_HOUR  ) /= 0) tmpstr = TRIM(tmpstr)//'hh'
      IF (IAND(IFile%flags,INPUT_FILE_NAME_MINUTE) /= 0) tmpstr = TRIM(tmpstr)//'ii'
      IF (IAND(IFile%flags,INPUT_FILE_NAME_SECOND) /= 0) tmpstr = TRIM(tmpstr)//'ss'
      tmpstr = TRIM(tmpstr)//'>'//TRIM(IFile%name_suf)
    ENDIF
    BuildGenericName = TRIM(tmpstr)

  END FUNCTION BuildGenericName

! Construct filename with time macros
  FUNCTION GetFileName(IFile,year,month,day,sod)
  CHARACTER (LEN=128) :: GetFileName
  TYPE (input_file_list), POINTER :: IFile
  INTEGER, INTENT(IN) :: year,month,day
  INTEGER, INTENT(in), OPTIONAL :: sod

  CHARACTER (LEN=7)   :: YearAdd, MonthAdd, DayAdd, HourAdd, MinAdd, SecAdd

    YearAdd  = ''
    MonthAdd = ''
    DayAdd   = ''
    HourAdd  = ''
    MinAdd   = ''
    SecAdd   = ''
    IF (IAND(IFile%flags,INPUT_FILE_NAME_YEAR) == INPUT_FILE_NAME_YEAR_SIGN) THEN
      YearAdd(1:1) = '+'
      IF (year < 0) YearAdd(1:1) = '-'
      WRITE(YearAdd(2:),'(I6.6)') ABS(year)
    ENDIF
    IF (IAND(IFile%flags,INPUT_FILE_NAME_YEAR) == INPUT_FILE_NAME_YEAR_LONG) WRITE(YearAdd ,'(I6.6)') year
    IF (IAND(IFile%flags,INPUT_FILE_NAME_YEAR) == INPUT_FILE_NAME_YEAR     ) WRITE(YearAdd ,'(I4.4)') year
    IF (IAND(IFile%flags,INPUT_FILE_NAME_MONTH) /= 0) WRITE(MonthAdd,'(I2.2)') month
    IF (IAND(IFile%flags,INPUT_FILE_NAME_DAY  ) /= 0) WRITE(DayAdd  ,'(I2.2)') day
    IF (PRESENT(sod)) THEN
      IF (IAND(IFile%flags,INPUT_FILE_NAME_HOUR  ) /= 0) WRITE(HourAdd,'(I2.2)') sod/3600
      IF (IAND(IFile%flags,INPUT_FILE_NAME_MINUTE) /= 0) WRITE(MinAdd ,'(I2.2)') mod(sod/60,60)
      IF (IAND(IFile%flags,INPUT_FILE_NAME_SECOND) /= 0) WRITE(SecAdd ,'(I2.2)') mod(sod,60)
    ENDIF

    GetFileName = TRIM(IFile%name_pre)//TRIM(YearAdd)//TRIM(MonthAdd)//TRIM(DayAdd)// &
                  TRIM(HourAdd)//TRIM(MinAdd)//TRIM(SecAdd)//TRIM(IFile%name_suf)

  END FUNCTION GetFileName

! --------------------------
!    Working with times
! --------------------------

! Modifies a variable time according to the time validity
  RECURSIVE FUNCTION ModVarTime(IFile,sel,tme,tp) RESULT(tio)
  INTEGER :: tio(4)
  TYPE (input_file_list), POINTER :: IFile
  INTEGER,             INTENT(in) :: sel
  INTEGER,             INTENT(in) :: tme(4)
  INTEGER, OPTIONAL,   INTENT(in) :: tp

  INTEGER :: vl
  INTEGER :: tmp1(4), tmp2(4)
  INTEGER(i8) :: dt

    vl = IAND(sel,INPUT_VAR_VALID_MASK)
    IF (PRESENT(tp)) vl = tp
    tmp1 = tme
    SELECT CASE (vl)
      CASE (INPUT_VAR_VALID_START)
        IF (IFile%dt_file > 0) THEN
          IF ((IFile%dt_file > 86400) .OR. ((86400/IFile%dt_file)*IFile%dt_file /= 86400)) CALL local_error('ModVarTime', &
            'Sorry. To round to start of an interval in seconds, the a day must be an integer #intervals')
          tmp1(4) = (tmp1(4)/IFile%dt_file) * IFile%dt_file
        ELSE
          tmp1(3) = 1
          tmp1(4) = 0
          IF (IFile%dt_file <= -12) THEN
            IF (MOD(-IFile%dt_file,12) /= 0) CALL local_error('ModVarTime', &
              'Sorry. To round to start of an interval > one year, the interval must be an integer number of years')
            tmp1(2) = 1
            tmp1(1) = (tmp1(1)/(IFile%dt_file/12)) * (IFile%dt_file/12)
          ELSE
            SELECT CASE (-IFile%dt_file)
              CASE (1)
              CASE (2,3,4,6)
                tmp1(2) = (tmp1(2)/IFile%dt_file) * IFile%dt_file
              CASE (5,7,8,9,10,11)
                CALL local_error('ModVarTime', &
                                 'Sorry. To round to start of an interval in months, one year must be an integer #intervals')
            END SELECT
          ENDIF
        ENDIF
        tio = tmp1
      CASE (INPUT_VAR_VALID_BEFORE)
        CALL CalendarTimeSub(tmp1,1_i8) ! Secure length of previous interval is used
        CALL CalendarTimeSub(tmp1,INT(IFile%dt_file,i8))
        CALL CalendarTimeAdd(tmp1,1_i8)
        tmp2 = ModVarTime(IFile,sel,tmp1,INPUT_VAR_VALID_MIDPOINT)
        tio = tmp2
      CASE (INPUT_VAR_VALID_ACTUAL)
        tio = tme
      CASE (INPUT_VAR_VALID_MIDPOINT)
        IF (IFile%dt_file > 0) THEN
          CALL CalendarTimeAdd(tmp1,INT(IFile%dt_file/2,i8))
        ELSE
          tmp2 = tmp1
          CALL CalendarTimeAdd(tmp2,INT(IFile%dt_file,i8))
          dt = CalendarTimeDiff(tmp1,tmp2)
          CALL CalendarTimeAdd(tmp1,dt/2)
        ENDIF
        tio = tmp1
      CASE (INPUT_VAR_VALID_END)
        CALL CalendarTimeAdd(tmp1,INT(IFile%dt_file,i8))
        tmp2 = ModVarTime(IFile,sel,tmp1,INPUT_VAR_VALID_START)
        CALL CalendarTimeSub(tmp2,INT(model_dt,i8))
        tio = tmp2
      CASE (INPUT_VAR_VALID_MIDDAY)
        tmp1(4) = 43200
        tio = tmp1
      CASE (INPUT_VAR_VALID_MIDMONTH)
        tmp1(3:4) = (/1,0/)
        tmp2 = tmp1
        CALL CalendarTimeAdd(tmp2,-1_i8)
        dt = CalendarTimeDiff(tmp1,tmp2)/2_i8
        CALL CalendarTimeAdd(tmp1,dt)
        tio = tmp1
      CASE (INPUT_VAR_VALID_MIDYEAR)
        tmp1(2:4) = (/1,1,0/)
        tmp2 = tmp1
        CALL CalendarTimeAdd(tmp2,-12_i8)
        dt = CalendarTimeDiff(tmp1,tmp2)/2_i8
        CALL CalendarTimeAdd(tmp1,dt)
        tio = tmp1
    END SELECT

  END FUNCTION ModVarTime

! Return the current time in the module Notice: always one time step ahead of the main model)
  FUNCTION InputTimeGet(lAdjust)
  INTEGER :: InputTimeGet(4)
  LOGICAL, OPTIONAL, INTENT(in) :: lAdjust

  INTEGER :: res(4)

    res = model_time
    IF (PRESENT(lAdjust)) THEN
      IF (lAdjust) CALL CalendarTimeSub(res,INT(model_dt,i8))
    ELSE
      CALL CalendarTimeSub(res,INT(model_dt,i8))
    ENDIF
    InputTimeGet = res

  END FUNCTION InputTimeGet

! Set module time
  SUBROUTINE InputTimeSet(year,month,day,SOD,dt_model,dt_unit,mod_time)
  INTEGER,                     INTENT(in)  :: year
  INTEGER,           OPTIONAL, INTENT(in)  :: month,day,SOD,dt_model
  CHARACTER (len=*), OPTIONAL, INTENT(in)  :: dt_unit
  INTEGER,           OPTIONAL, INTENT(out) :: mod_time(4)

  INTEGER :: mod_tim(4)
    IF (lInLoop) CALL local_error('InputTimeSet','Cannot change time after entering main loop')

    IF (PRESENT(dt_model) .OR. PRESENT(dt_unit)) model_dt = InterpretDtUnit(dt_model,dt_unit)

    mod_tim  = (/year,1,1,0/)
    IF (PRESENT(month)) mod_tim(2) = month
    IF (PRESENT(day  )) mod_tim(3) = day
    IF (PRESENT(SOD  )) mod_tim(4) = SOD
    IF (PRESENT(mod_time)) mod_time = mod_tim

  END SUBROUTINE InputTimeSet

! Return time (default: adjusted current input time) as a string in the selected input debug message format 
  FUNCTION InputTimeString(lAdjust,Time)
  CHARACTER (len=19) :: InputTimeString
  LOGICAL, OPTIONAL, INTENT(in) :: lAdjust
  INTEGER, OPTIONAL, INTENT(in) :: Time(4)

  INTEGER :: Tim(4)

    IF (PRESENT(Time)) THEN
      Tim = Time
    ELSE
      Tim = InputTimeGet(lAdjust)
    ENDIF
    InputTimeString = DateString(Tim,IAND(mo_debug,INPUT_DBG_TIME_MASK))

  END FUNCTION InputTimeString

! ---------------------------------
!  Build-in interpolation routines 
! ---------------------------------

! Return fields and times necessary for user interpolation
  SUBROUTINE InputVarGetInterpolData(Var,Fields,t,t1,t2,t3,t4,dt1,dt2,dt3,dt4,d,d1,d2,d3,d4)
  TYPE (input_var_list), POINTER :: Var
  INTEGER, OPTIONAL, INTENT(IN),  TARGET :: Fields(:)
  INTEGER, OPTIONAL, INTENT(OUT) :: t(4)
  INTEGER, OPTIONAL, INTENT(OUT) :: t1(4), t2(4), t3(4), t4(4)
  INTEGER, OPTIONAL, INTENT(OUT), TARGET :: dt1, dt2, dt3, dt4
  REAL(dp),OPTIONAL, POINTER     :: d(:,:,:,:,:,:,:), d1(:,:,:,:,:,:,:), d2(:,:,:,:,:,:,:), d3(:,:,:,:,:,:,:), d4(:,:,:,:,:,:,:)

  TYPE (input_data_list), POINTER :: curr
  INTEGER, POINTER :: EffFields(:)
  INTEGER, POINTER :: pdt
  INTEGER :: i, j, Dst(3), iCurr, iMax, Tme(4)

    iMax = MIN(4,Var%nIntrp)
    IF (PRESENT(Fields)) THEN
      EffFields => Fields
      iMax = MIN(size(Fields),iMax)
    ELSE
      ALLOCATE(EffFields(iMax))
      DO i=1,iMax
        EffFields(i) = i
      ENDDO
    ENDIF

    IF (PRESENT(dt1) .OR. PRESENT(dt2) .OR. PRESENT(dt3) .OR. PRESENT(dt4) .OR. PRESENT(t)) THEN
      Tme = model_time
      CALL CalendarTimeSub(Tme,INT(model_dt,i8))
      IF (PRESENT(t)) t = Tme
    ENDIF
   
    IF (.NOT. ASSOCIATED(Var%Intrp)) CALL local_error('InputVarGetInterpolData', &
      'Attempt to get data from non-interpolated variable: '//TRIM(Var%name_var)) 
    curr => Var%Intrp%prev
    iCurr = 1
    DO i=1,iMax
      Dst(1) = iCurr-EffFields(i)
      Dst(2) = Dst(1) - var%nIntrp
      IF (ABS(Dst(2)) < ABS(Dst(1))) Dst(1) = Dst(2)
      Dst(2) = EffFields(i)
      Dst(3) = EffFields(i) - var%nIntrp
      IF (Dst(2) < ABS(Dst(3))) Dst(2) = Dst(3)
      IF (ABS(Dst(2)) < ABS(Dst(1))) THEN
        curr => Var%Intrp
        Dst(1) = Dst(2)
      ENDIF
      IF (Dst(1) < 0) THEN
        DO j=1,-Dst(1)
          curr => curr%prev
        ENDDO
      ELSE
        DO j=1,Dst(1)
          curr => curr%next
        ENDDO
      ENDIF
      iCurr = EffFields(i)
      NULLIFY(pdt)
      SELECT CASE (i)
        CASE (1)
          IF (PRESENT( t1)) t1 = curr%time
          IF (PRESENT(dt1)) pdt => dt1
          IF (PRESENT( d1)) d1 => curr%dta
        CASE (2)
          IF (PRESENT( t2)) t2 = curr%time
          IF (PRESENT(dt2)) pdt => dt2
          IF (PRESENT( d2)) d2 => curr%dta
        CASE (3)
          IF (PRESENT( t3)) t3 = curr%time
          IF (PRESENT(dt3)) pdt => dt2
          IF (PRESENT( d3)) d3 => curr%dta
        CASE (4)
          IF (PRESENT( t4)) t4 = curr%time
          IF (PRESENT(dt4)) pdt => dt4
          IF (PRESENT( d4)) d4 => curr%dta
      END SELECT
      IF (ASSOCIATED(pdt)) THEN
        IF (curr%time(2) > 0) THEN
          IF (CalendarTimeCmp(curr%time,tme) < 0) THEN
            pdt =  CalendarTimeDiff(tme,curr%time)
          ELSE
            pdt = -CalendarTimeDiff(curr%time,tme)
          ENDIF
        ELSE
          pdt = 0
        ENDIF
      ENDIF
    ENDDO

    IF (PRESENT(d)) d => var%dta

  END SUBROUTINE InputVarGetInterpolData

! Make a boxcar "interpolation" by assuming that the last field before or at current time is valid
  SUBROUTINE InputVarInterpolBoxCar(var)
  TYPE (input_var_list), POINTER :: var

  REAL(dp), POINTER :: d(:,:,:,:,:,:,:), d1(:,:,:,:,:,:,:)
  INTEGER :: dt1, Curr

    Curr = var%nIntrp - var%nFuture ! Maximum field number which can be in the past
    DO WHILE (Curr > 0)
      CALL InputVarGetInterpolData(var,Fields=(/Curr/),dt1=dt1,d=d,d1=d1)
      IF (dt1 <= 0) THEN
        d(:,:,:,:,:,:,:) = d1(:,:,:,:,:,:,:)
        RETURN
      ENDIF
      Curr = Curr - 1
    ENDDO
    CALL local_error('InputVarInterpolBoxCar','No past data found')

  END SUBROUTINE InputVarInterpolBoxCar

! Make a linear interpolation between the last field before current time and the first one after
  SUBROUTINE InputVarInterpolLinear(var)
  TYPE (input_var_list), POINTER :: var

  REAL(dp), POINTER :: d(:,:,:,:,:,:,:), d1(:,:,:,:,:,:,:), d2(:,:,:,:,:,:,:)
  REAL(dp) :: frac
  INTEGER :: dt1, dt2, Curr

    Curr = 1
    DO WHILE (Curr < var%nIntrp-1)
      CALL InputVarGetInterpolData(var,Fields=(/Curr,Curr+1/),dt1=dt1,dt2=dt2,d=d,d1=d1,d2=d2)
      IF (((dt1 <= 0) .AND. (dt2 >  0)) .OR. &
          ((dt1 <  0) .AND. (dt2 >= 0))) THEN ! One field in the past, the other in the future
        frac = -REAL(dt1,dp) / ( REAL(dt2,dp) - REAL(dt1,dp) )
        d(:,:,:,:,:,:,:) = d1(:,:,:,:,:,:,:) * (1.0_dp - frac) + d2(:,:,:,:,:,:,:) * frac
        RETURN
      ENDIF
      Curr = Curr + 1
    ENDDO
    CALL local_error('InputVarInterpolLinear','No past or no future data found')

  END SUBROUTINE InputVarInterpolLinear

END MODULE mo_input_misc
