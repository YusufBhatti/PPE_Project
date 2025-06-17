!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_input_strings
#ifdef INPUT_IN_ECHAM
  USE mo_exception,   ONLY: finish, message, message_text
  USE mo_util_string, ONLY: tolower
#endif
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
  USE mo_kind,        ONLY: dp
IMPLICIT NONE
#else
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
#endif

! String functions for mo_input
! V.1.00, 2013-03-05, Stiig Wilkenskjeld - MPI-M

PUBLIC local_message
PUBLIC local_error
PUBLIC CountSubStrings
PUBLIC GetStringIndex
PUBLIC GetStringIndexed
PUBLIC InterpretDtUnit
PUBLIC TimeStepString
PUBLIC GetDateFromString
PUBLIC DateString
PUBLIC DispFlagsSingleBit
PUBLIC GetFlagsSingleBit
PUBLIC SetLogUnit
PUBLIC GetLogUnit
PUBLIC RmSpc
PUBLIC RmDblSpc
PUBLIC String2Vec
#ifndef INPUT_IN_ECHAM
  PUBLIC tolower
  CHARACTER (len=256), PUBLIC :: message_text
#endif

INTEGER, PRIVATE :: local_log_unit = -1
LOGICAL, PRIVATE :: lstdout = .TRUE.
LOGICAL, PRIVATE :: lflush  = .FALSE.

CONTAINS

! Set the unit for the PE-local log file
  SUBROUTINE SetLogUnit(unt,stdout,flush)
  INTEGER, INTENT(in) :: unt
  LOGICAL, INTENT(in), OPTIONAL :: stdout
  LOGICAL, INTENT(in), OPTIONAL :: flush
  
    local_log_unit = unt
    IF (PRESENT(stdout)) lstdout = stdout
    IF (PRESENT(flush )) lflush  = flush 

  END SUBROUTINE SetLogUnit

! Set the unit for the PE-local log file
  INTEGER FUNCTION GetLogUnit()
  
    IF (local_log_unit < 0) THEN
      GetLogUnit = 6
    ELSE
      GetLogUnit = local_log_unit
    ENDIF

  END FUNCTION GetLogUnit

! Local error: Here just write message and stop. Collected here for easy implementation into different projects
  SUBROUTINE local_message(subprg,text)
  CHARACTER (len=*), INTENT(in) :: subprg, text

    IF ((local_log_unit > 0) .AND. (local_log_unit /= 6)) THEN
      IF (subprg == '') THEN
        WRITE(local_log_unit,'(a)') TRIM(text)
      ELSE
        WRITE(local_log_unit,'(a,": ",a)') TRIM(subprg), TRIM(text) ! Local test version
      ENDIF
    ENDIF
#ifndef NO_FLUSH_AVAIL
    IF (lflush) flush(local_log_unit)
#endif
    IF (lstdout) THEN 
#ifdef INPUT_IN_ECHAM
      CALL message(subprg,text) ! cbalone/jsbach version
#else
      IF (subprg == '') THEN
        WRITE(*,'(a)') TRIM(text)
      ELSE
        WRITE(*,'(a,": ",a)') TRIM(subprg), TRIM(text) ! Local test version
      ENDIF
#endif
    ENDIF

  END SUBROUTINE local_message

! Local error: Here just write message and stop. Collected here for easy implementation into different projects
  SUBROUTINE local_error(subprg,text)
  CHARACTER (len=*), INTENT(in) :: subprg, text

    IF ((local_log_unit > 0) .AND. (local_log_unit /= 6)) THEN
      IF (subprg == '') THEN
        WRITE(local_log_unit,'(a)') TRIM(text)
      ELSE
        WRITE(local_log_unit,'(a,": ",a)') TRIM(subprg), TRIM(text) ! Local test version
      ENDIF
    ENDIF
    IF (lstdout) THEN
#ifdef INPUT_IN_ECHAM
      CALL finish(subprg,text) ! cbalone/jsbach version
#else
      WRITE(*,'(a,": ",a)') subprg, text ! Local test version
      STOP 1
#endif
    ENDIF

  END SUBROUTINE local_error

! Count number of substrings in a string separated by character (default ',')
  INTEGER FUNCTION CountSubStrings(str,sep)
  CHARACTER (len=*), INTENT(IN)           :: str
  CHARACTER (len=1), INTENT(IN), OPTIONAL :: sep
  
  CHARACTER (len=1) :: rsep
  INTEGER           :: n,curr,ind,le
  
    rsep = ','
    IF (PRESENT(sep)) rsep = sep
  
    n = 0
    ind = 1
    curr = 1
    le = LEN_TRIM(str)
    DO WHILE (curr<le)
      n   = n + 1
      ind = INDEX(str(curr+1:),rsep) + curr - 1
      IF (ind < curr) ind = le + 1
      curr = ind + 1
    ENDDO

    CountSubStrings = n  
  
  END FUNCTION CountSubStrings

#ifndef INPUT_IN_ECHAM
! Conversion to lowercase (adopted from echam: mo_util_string.f90)
  FUNCTION tolower(str)
  CHARACTER(LEN=*)            ,INTENT(in) :: str
  CHARACTER(LEN=LEN_TRIM(str))            :: tolower

  INTEGER            :: i
  INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

    DO i=1,LEN_TRIM(str)
      IF (ICHAR(str(i:i)) >= ICHAR('A') .AND. &
          ICHAR(str(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(str(i:i)) + idel )
      ELSE
        tolower(i:i) = str(i:i)
      END IF
    END DO

  END FUNCTION tolower
#endif

! Returns the index of a string in a list of strings
  INTEGER FUNCTION GetStringIndex(Search,List,lIgnoreS,Sep,ErrMsg,lCaseS,St)
  CHARACTER (len=*), INTENT(in) :: Search, List
  LOGICAL, OPTIONAL, INTENT(in) :: lIgnoreS, lCaseS
  CHARACTER (len=1), INTENT(in), OPTIONAL :: Sep
  CHARACTER (len=*), INTENT(in), OPTIONAL :: ErrMsg
  INTEGER,          INTENT(out), OPTIONAL :: St

  CHARACTER (len=1) :: LSep
  CHARACTER (len=len_trim(Search)) :: Str
  INTEGER           :: Start, iEnd, Res, Tmp
  LOGICAL           :: lCont

    GetStringIndex = -1

    IF (Search == ' ') THEN
      GetStringIndex = 0
      RETURN
    ENDIF

    ! Defaults for optional parameters
    LSep = " "
    IF (PRESENT(Sep)) LSep = Sep

    iEnd = LEN_TRIM(Search)
    IF (PRESENT(lIgnoreS)) THEN
      IF ((lIgnoreS) .AND. (Search(iEnd:iEnd)=='s')) iEnd = iEnd-1
    ENDIF

    ! Compare strings in lower case
    IF (PRESENT(lCaseS)) THEN
      IF (lCaseS) THEN
        Str = TRIM(Search)
      ELSE
        Str = tolower(Search)
      ENDIF
    ELSE
      Str = tolower(Search)
    ENDIF

    ! Find string position which is surrounded by LSep's
    Start = 1
    lCont  = .TRUE.
    DO WHILE (lCont)
      Res = INDEX(List(Start:),Str(1:iEnd)) + Start - 1
      IF (Res < Start) THEN ! String not found
        IF (PRESENT(ErrMsg)) CALL local_error('GetStringIndex','Invalid string for '//ErrMsg//': '//TRIM(Search)// &
                                              '. Valid options are: '//TRIM(List))
        RETURN
      ENDIF
      IF (Res == 1) THEN 
        lCont = .FALSE.
      ELSE
        IF (List(Res-1:Res-1)==LSep) lCont = .FALSE.
      ENDIF
      IF (.NOT. lCont) THEN
        Tmp = Res + iEnd
        IF (iEnd < LEN_TRIM(List) .AND. (iEnd < LEN_TRIM(Search))) THEN
          IF (List(Tmp:Tmp+1) /= 's'//LSep) lCont = .TRUE.
        ELSE
          IF (Tmp <= LEN_TRIM(List)) THEN
            IF (List(Tmp:Tmp) /= LSep) lCont = .TRUE.
          ENDIF
        ENDIF
      ENDIF
      Start = Res + 1
    ENDDO
    IF (PRESENT(St)) St = Start - 1

    ! Count # LSep's before found occurrence
    Tmp = 1
    Start = 1
    DO WHILE (Start < Res)
      iEnd = INDEX(List(Start:),LSep) + Start - 1
      IF (iEnd < Start) iEnd = LEN_TRIM(List) + 1
      Start = iEnd + 1
      Tmp = Tmp + 1
    ENDDO

    GetStringIndex = Tmp

  END FUNCTION GetStringIndex

! Return the substring with the given number in a list of strings
  FUNCTION GetStringIndexed(Strs,Ind,Sep)
  CHARACTER (len=16)             :: GetStringIndexed
  CHARACTER (len=* ), INTENT(in) :: Strs
  INTEGER,            INTENT(in) :: Ind
  CHARACTER (len=1 ), INTENT(in), OPTIONAL :: Sep

  CHARACTER (len=1) :: RSep
  INTEGER :: i, Curr, Next

    RSep = ' '
    IF (PRESENT(Sep)) RSep = Sep

    GetStringIndexed = ''
    Curr = 1
    DO i=1,Ind
      Next = INDEX(Strs(Curr:),RSep) + Curr - 1
      IF (Next < Curr) THEN
        IF (i < Ind) RETURN
        Next = LEN_TRIM(Strs) + 1
      ENDIF
      IF (i < Ind) Curr = Next + 1
    ENDDO

    GetStringIndexed = Strs(Curr:Next-1)

  END FUNCTION GetStringIndexed

! Return a time step in model format given an integer time step and eventual a unit for the time step
  INTEGER FUNCTION InterpretDtUnit(dt,Units,def_val)
  INTEGER,          INTENT(in), OPTIONAL :: dt
  CHARACTER(len=*), INTENT(in), OPTIONAL :: Units
  INTEGER,          INTENT(in), OPTIONAL :: def_val

  CHARACTER (len=*), PARAMETER :: ShortUnits = 'sSiIhHdDmMyYeEcClL'
  CHARACTER (len=*), PARAMETER :: UnitStrs   = 'second minute hour day month year decade centuries millenia file'
  INTEGER,           PARAMETER :: UnitMul(9) = (/1,    60,    3600,86400, -1, -12,  -120,    -1200,  -12000/)

  CHARACTER (len=4) :: SFmt
  INTEGER :: Unt, i, St, Tmp

    InterpretDtUnit = 0
    Tmp = 1
    IF (PRESENT(dt)) THEN
      Tmp = dt
      InterpretDtUnit = dt
    ENDIF
    IF (.NOT. PRESENT(Units)) RETURN
    IF (LEN_TRIM(Units)==0) RETURN

    St = 1
    IF (INDEX('1234567890',Units(1:1)) > 0) THEN
      DO WHILE (INDEX('1234567890',Units(St:St)) > 0)
        St = St + 1
      ENDDO
      IF (St > 10) CALL local_error('InterpretDtUnit','Too large numeric time constant as string')
      WRITE(SFmt,'(a,i1,a)') '(i',St-1,')'
      READ(Units(1:St-1),SFmt) i
      DO WHILE ((Units(St:St)==' ') .OR. (ICHAR(Units(St:St))==9))
        St = St + 1
      ENDDO
      Tmp = i
    ENDIF

    IF (LEN_TRIM(Units(St:))==1) THEN
      Unt = (INDEX(ShortUnits,Units(St:St))+1)/2
      IF (Unt == 0) THEN
        InterpretDtUnit = 0
        RETURN
      ENDIF
    ELSE
      i = LEN_TRIM(Units(St:))
      IF ((Units(i+St-1:i+St-1) == 's') .OR. (Units(i+St-1:i+St-1) == 'S')) i = i - 1
      Unt = GetStringIndex(Units(St:i+St-1),UnitStrs,lIgnoreS=.TRUE.,ErrMsg='time unit')
      IF (Unt < 0) THEN
        InterpretDtUnit = 0
        RETURN
      ENDIF
      IF ((Unt == 10) .AND. (PRESENT(def_val))) THEN
        InterpretDtUnit = def_val
        RETURN
      ENDIF
    ENDIF

    InterpretDtUnit = Tmp * UnitMul(Unt)

  END FUNCTION InterpretDtUnit

! Return a string with a time interval "humanly scaled"
  FUNCTION TimeStepString(dt,bad_val)
  CHARACTER (len=11)  :: TimeStepString
  INTEGER, INTENT(in) :: dt
  INTEGER, INTENT(in), OPTIONAL :: bad_val

  INTEGER :: Tmp, i
  CHARACTER (len=50) :: Str

    Tmp = dt
    IF (PRESENT(bad_val)) THEN
      IF (dt == bad_val) Tmp = HUGE(dt)
    ENDIF
    IF (ABS(Tmp)==HUGE(Tmp)) THEN
      TimeStepString = '<undefined>'
      RETURN
    ENDIF

    Str = ''
    IF (Tmp < 0) THEN ! Time step in months
      IF (Tmp < -11) THEN
        i = 1
        IF ((Tmp/12)*12 /= Tmp) THEN
          Str(1:1) = '~'
          i = 2
        ENDIF
        WRITE(Str(i:),*) -Tmp/12,'years'
      ELSE
        WRITE(Str,*) -Tmp,'months'
      ENDIF
      TimeStepString = TRIM(Str)
      RETURN
    ELSE      
      i = 1
      IF (Tmp >= 86400) THEN
        WRITE(Str(i:),*) Tmp/86400,'d'
        Tmp = MODULO(Tmp,86400)
        i   = LEN_TRIM(Str) + 1
      ENDIF
      IF (Tmp >= 3600) THEN
        WRITE(Str(i:),*) Tmp/3600,'h'
        Tmp = MODULO(Tmp,3600)
        i   = LEN_TRIM(Str) + 1
      ENDIF
      IF (Tmp >= 60) THEN
        WRITE(Str(i:),*) Tmp/60,'m'
        Tmp = MODULO(Tmp,60)
        i   = LEN_TRIM(Str) + 1
      ENDIF
      IF ((Tmp > 0) .OR. (i==1)) WRITE(Str(i:),*) Tmp,'s'
    ENDIF
    CALL RmSpc(Str)
    TimeStepString = TRIM(Str)

  END FUNCTION TimeStepString

! Convert a 4-11 digit string to a date
  SUBROUTINE GetDateFromString(Str,Date,ErrMsg)
  CHARACTER (len=*), INTENT(in   ) :: Str
  INTEGER,           INTENT(inout) :: Date(3)
  CHARACTER (len=*), INTENT(in   ) :: ErrMsg

  INTEGER i, st;

    st = 1
    IF (Str(1:1)=='-' .OR. Str(1:1)=='+') st = 4

    i = LEN_TRIM(Str(st:))
    SELECT CASE (i)
      CASE (0,4,6,8)
        ! These are the allowed ones
      CASE DEFAULT
        CALL local_error('GetDateFromString','Invalid date format for '//TRIM(ErrMsg)//': '//TRIM(Str))
    END SELECT
    IF (i >= 4) READ(Str(st  :st+3),'(i4)') Date(1)
    IF (i >= 6) READ(Str(st+4:st+5),'(i2)') Date(2)
    IF (i >= 8) READ(Str(st+6:st+7),'(i2)') Date(3)
    IF (st > 1) THEN
      READ(Str(2:3),'(i2)') i
      Date(1) = Date(1) + i * 10000
      IF (Str(1:1)=='-') Date(1) = -Date(1)
    ENDIF

  END SUBROUTINE GetDateFromString

! Display names of set single bit flags occupying a number of consecutive bits
  SUBROUTINE DispFlagsSingleBit(Flags,FlagLow,Cnt,StrList,PreText,Txt,Sep,lEmpty)
  INTEGER,                     INTENT(in)  :: Flags, FlagLow, Cnt
  CHARACTER (len=*),           INTENT(in)  :: StrList
  CHARACTER (len=*), OPTIONAL, INTENT(in)  :: PreText
  CHARACTER (len=*), OPTIONAL, INTENT(out) :: Txt
  CHARACTER (len=1), OPTIONAL, INTENT(in)  :: Sep
  LOGICAL,           OPTIONAL, INTENT(in)  :: lEmpty

  INTEGER CurrList, NextList, CurrFlag, i
  CHARACTER (len=1) :: RSep
  LOGICAL Empty

    RSep = " "
    IF (PRESENT(Sep)) RSep = Sep
    message_text = ""
    CurrList = 1
    CurrFlag = FlagLow
    DO i=1,Cnt
      NextList = INDEX(StrList(CurrList:),RSep) + CurrList - 1
      IF (NextList < CurrList) NextList = LEN_TRIM(StrList) + 1
      IF (IAND(Flags,CurrFlag) /= 0) message_text = TRIM(message_text)//" "//StrList(CurrList:NextList-1)
      IF (i<Cnt) CurrFlag = CurrFlag * 2
      CurrList = NextList + 1
    ENDDO

    Empty = .FALSE.
    IF (PRESENT(lEmpty)) Empty = lEmpty
    IF (.NOT. Empty .AND. LEN_TRIM(message_text) == 0) message_text = " <none>"

    IF (PRESENT(Txt)) THEN
      Txt = message_text(2:)
    ELSE
      IF (PRESENT(PreText)) THEN
        message_text = TRIM(PreText)//TRIM(message_text(2:))
        CALL local_message('',message_text)
      ELSE
        CALL local_message('DispFlagsSingleBit',message_text(2:))
      ENDIF
    ENDIF

  END SUBROUTINE DispFlagsSingleBit

! Get multi-flags from string
  INTEGER FUNCTION GetFlagsSingleBit(Txt,FlagLow,StrList,lLin,ErrMsg,Sep,Ofs)
  CHARACTER (len=*),           INTENT(in) :: Txt
  INTEGER,                     INTENT(in) :: FlagLow
  CHARACTER (len=*),           INTENT(in) :: StrList
  LOGICAL,           OPTIONAL, INTENT(in) :: lLin
  CHARACTER (len=*), OPTIONAL, INTENT(in) :: ErrMsg
  CHARACTER (len=1), OPTIONAL, INTENT(in) :: Sep
  INTEGER,           OPTIONAL, INTENT(in) :: Ofs

  INTEGER :: i, Res, st, en, iofs
  LOGICAL :: lCont
  CHARACTER (len=1) :: RSep

    RSep = ","
    IF (PRESENT(Sep)) RSep = Sep
    iofs = 1
    IF (PRESENT(Ofs)) iofs = Ofs
    Res = 0

    IF (LEN_TRIM(Txt) > 0) THEN
      st = 1
      lCont = .TRUE.
      DO WHILE (lCont)
        en = INDEX(Txt(st:),',') + st - 1
        IF (en < st) THEN
          en = LEN_TRIM(Txt) + 1
          lCont = .FALSE.
        ENDIF
        i = GetStringIndex(Txt(st:en-1),StrList,ErrMsg=ErrMsg,sep=RSep)
        IF (PRESENT(lLin)) THEN
          IF (lLin) THEN
            IF ((i-iofs >= 0) .AND. (i > 0)) Res = IOR(Res,FlagLow*(i-iofs))
            i = 0
          ENDIF
        ENDIF
        IF ((i-iofs >= 0) .AND. (i > 0)) Res = IOR(Res,ISHFT(FlagLow,(i-iofs)))
        st = en + 1
      ENDDO
    ENDIF

    GetFlagsSingleBit = Res

  END FUNCTION GetFlagsSingleBit

! Writes a date in mo_input_calendar time format as string. Special output for year = +-HUGE(i4)
  FUNCTION DateString(Date,Flag)
  CHARACTER (len=19)            :: DateString
  INTEGER, INTENT(in)           :: Date(4)
  INTEGER, INTENT(in), OPTIONAL :: Flag

  CHARACTER (len=16), SAVE :: fmts = '(i4.4,5(a,i2.2))'
  INTEGER, SAVE :: add = 0
  INTEGER, PARAMETER :: Pos(0:6) = (/19,16,13,10,7,4,0/)
  INTEGER f
  CHARACTER (len=22) :: tmp

    f = 0
    IF (PRESENT(Flag)) f = Flag
    IF ((f < 0) .OR. (f > 5)) f = 0

    IF (ABS(Date(1)) /= HUGE(Date(1))) THEN
      IF (Date(1)>9999) THEN
        IF (fmts(3:3) == '4') THEN
          fmts(3:5) = '6.6'
          add = 2
        ELSE
          add = 3
        ENDIF
        f = MAX(1,f)
      ENDIF
      IF (Date(1)<0) THEN
        fmts(3:5) = '7.6'
        f = MAX(1,f)
        IF (f < 6) add = 3
      ENDIF
      WRITE(tmp,fmts) Date(1),'-',Date(2),'-',Date(3),' ', &
                                    Date(4)/3600,':',MODULO(Date(4),3600)/60,':',MODULO(Date(4),60)
      DateString = TRIM(tmp(1:Pos(f)+add))
      IF (f==2) DateString = TRIM(DateString(1:19))//'h'
    ELSE 
      IF (DATE(1)<0) THEN
        SELECT CASE (f)
          CASE (0)!       1234567890123456789
            DateString = "   <undefined time>"
          CASE (1)
            DateString = "<undefined time>"
          CASE (2)
            DateString = "  <undefined>"
          CASE (3)
            DateString = "   <undef>"
          CASE (4)
            DateString = "<undef>"
          CASE (5)
            DateString = "----"
        END SELECT
      ELSE
        SELECT CASE (f)
          CASE (0)!       1234567890123456789
            DateString = "  <end of all time>"
          CASE (1)
            DateString = "    <time's end>"
          CASE (2)
            DateString = " <time's end>"
          CASE (3)
            DateString = '<infinite>'
          CASE (4)
            DateString = '  <inf>'
          CASE (5)
            DateString = ' inf'
        END SELECT
      ENDIF
    ENDIF

  END FUNCTION DateString

! Removes unwanted spaces from string. Spaces will be substituted for given special character (default "_")
  SUBROUTINE RmSpc(Str,Spc)
  CHARACTER (len=*), INTENT(inout)        :: Str
  CHARACTER (len=1), INTENT(in), OPTIONAL :: Spc

  INTEGER :: i, j
  CHARACTER (len=1) :: Sub

    Sub = '_'
    IF (PRESENT(Spc)) Sub = Spc

    j = 1
    DO i=1,LEN_TRIM(Str)
      IF (Str(i:i) /= ' ') THEN
        Str(j:j) = Str(i:i)
        j = j + 1
      ENDIF
    ENDDO
    Str(j:LEN_TRIM(Str)) = ' '
    DO i=1,LEN_TRIM(Str)
      IF (Str(i:i) == Sub) Str(i:i) = ' '
    ENDDO

  END SUBROUTINE RmSpc

! Removes duplicated spaces from string.
  SUBROUTINE RmDblSpc(Str)
  CHARACTER (len=*), INTENT(inout) :: Str

  INTEGER :: i, j

    j = 2
    DO i=2,LEN_TRIM(Str)
      IF ((Str(i:i) /= ' ') .OR. (Str(j-1:j-1) /= ' ')) THEN
        Str(j:j) = Str(i:i)
        j = j + 1
      ENDIF
    ENDDO
    Str(j:LEN_TRIM(Str)) = ' '

  END SUBROUTINE RmDblSpc

! Converts a string with comma separated real values to a numeric vector of appropriate length
  FUNCTION String2Vec(Str)
  REAL (dp), POINTER :: String2Vec(:)
  CHARACTER (len=*), INTENT(in) :: Str

  REAL (dp), POINTER :: Tmp(:)
  INTEGER :: i, st, en
  LOGICAL :: lCont

    i = 0
    st = 1
    lCont = .TRUE.
    DO WHILE (lCont)
      en = INDEX(Str(st:),',') + st - 1
      IF (en < st) THEN
        en = LEN_TRIM(Str) + 1
        lCont = .FALSE.
      ENDIF
      i = i + 1
      st = en + 1
    ENDDO
    ALLOCATE(Tmp(i))
    i  = 1
    st = 1
    lCont = .TRUE.
    DO WHILE (lCont)
      en = INDEX(Str(st:),',') + st - 1
      IF (en < st) THEN
        en = LEN_TRIM(Str) + 1
        lCont = .FALSE.
      ENDIF
      READ(Str(st:en-1),'(G50.48)') Tmp(i) ! Is it really necessary with such a crap format desc?
      i  = i  + 1
      st = en + 1
    ENDDO

    String2Vec => Tmp

  END FUNCTION String2Vec

END MODULE mo_input_strings
