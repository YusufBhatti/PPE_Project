!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_input_calendar
USE mo_kind, ONLY: i8, dp
IMPLICIT NONE

! V.1.00, 2013-03-05, Stiig Wilkenskjeld MPI-M, extracted from mo_input

PUBLIC CalendarYearLen
PUBLIC CalendarMonthLen
PUBLIC CalendarTimeDiff
PUBLIC CalendarTimeDiffSec
PUBLIC CalendarTimeAdd
PUBLIC CalendarTimeAddSec
PUBLIC CalendarTimeSub
PUBLIC CalendarTimeSubSec
PUBLIC CalendarTimeCmp
PUBLIC CalendarSet

! Constant minimum and average length of a month in the Gregorian calendar (s)
INTEGER, PARAMETER :: GREGORIAN_MIN_MONTH  =   2419200
INTEGER, PARAMETER :: GREGORIAN_MEAN_MONTH =   2629746
PUBLIC GREGORIAN_MIN_MONTH, GREGORIAN_MEAN_MONTH
! Calendar enumeration
INTEGER, PARAMETER :: INPUT_CALENDAR_UNDEFINED =     0
INTEGER, PARAMETER :: INPUT_CALENDAR_STANDARD  =     1
INTEGER, PARAMETER :: INPUT_CALENDAR_GREGORIAN =     1 ! Actually proleptic-Gregorian
INTEGER, PARAMETER :: INPUT_CALENDAR_JULIAN    =     2
INTEGER, PARAMETER :: INPUT_CALENDAR_DAY360    =     3
INTEGER, PARAMETER :: INPUT_CALENDAR_DAY365    =     4
INTEGER, PARAMETER :: INPUT_CALENDAR_DAY366    =     5
CHARACTER (len=*), PARAMETER :: CalendarNames = 'undefined Gregorian Julian 360-Day 365-Day 366-Day'
PUBLIC INPUT_CALENDAR_UNDEFINED, INPUT_CALENDAR_STANDARD, INPUT_CALENDAR_GREGORIAN, INPUT_CALENDAR_JULIAN
PUBLIC INPUT_CALENDAR_DAY360, INPUT_CALENDAR_DAY365, INPUT_CALENDAR_DAY366

INTEGER, PARAMETER, PRIVATE :: nYearDays(3:5) = (/360,365,366/)
INTEGER, PARAMETER, PRIVATE :: nFebDays (4:5) = (/     28, 29/)
INTEGER,            PRIVATE :: calendar  = INPUT_CALENDAR_UNDEFINED
INTEGER,            PRIVATE :: FixedYearLen = 0

CONTAINS

! ------------------------------------
!  Calendar functions and subroutines
! ------------------------------------

! Set current calendar
  LOGICAL FUNCTION CalendarSet(cal)
  INTEGER, INTENT(in) :: cal

    IF ((cal<1) .OR. (cal>INPUT_CALENDAR_DAY366)) THEN
      CalendarSet = .FALSE.
    ELSE
      CalendarSet = .TRUE.
      calendar = cal
      IF (cal > INPUT_CALENDAR_JULIAN) FixedYearLen = nYearDays(cal)
    ENDIF

  END FUNCTION CalendarSet

! Get the current calendar
  INTEGER FUNCTION CalendarGet(CName)
  USE mo_input_strings
  CHARACTER (len=*), OPTIONAL, INTENT(out) :: CName

    CalendarGet = calendar
    IF (PRESENT(CName)) CName = GetStringIndexed(CalendarNames,calendar+1,' ')

  END FUNCTION CalendarGet

! Returns the length of a year in the (proleptic) Gregorian or the (proleptic) Julian calendar dependent on the calendar setting
  INTEGER FUNCTION CalendarYearLen(Year)
  INTEGER, INTENT(IN) :: Year

    IF (calendar > INPUT_CALENDAR_JULIAN) THEN
      CalendarYearLen = FixedYearLen
      RETURN
    ENDIF
    IF ((MOD(Year,4)==0) .AND. ((calendar==INPUT_CALENDAR_JULIAN) .OR. (MOD(Year,100) /= 0) .OR. (MOD(Year,400) == 0))) THEN
      CalendarYearLen = 366
    ELSE
      CalendarYearLen = 365
    ENDIF

  END FUNCTION CalendarYearLen

! Returns the length of a month in the (proleptic) Gregorian or the (proleptic) Julian calendar dependent on the ljulian switch
  INTEGER FUNCTION CalendarMonthLen(Year,Month)
  INTEGER, INTENT(IN) :: Year, Month

    CalendarMonthLen = 0
    IF (calendar==INPUT_CALENDAR_DAY360) THEN
      CalendarMonthLen = 30
      RETURN
    ENDIF
    SELECT CASE (Month)
      CASE (1,3,5,7,8,10,12)
        CalendarMonthLen = 31
      CASE (4,6,9,11)
        CalendarMonthLen = 30
      CASE (2)
        IF (calendar > INPUT_CALENDAR_JULIAN) THEN
          CalendarMonthLen = nFebDays(calendar)
          RETURN
        ENDIF
        IF ((MOD(Year,4)==0) .AND. ((calendar==INPUT_CALENDAR_JULIAN) .OR. (MOD(Year,100)/=0) .OR. (MOD(Year,400)==0))) THEN
          CalendarMonthLen = 29
        ELSE
          CalendarMonthLen = 28
        ENDIF
    END SELECT

  END FUNCTION CalendarMonthLen

! Time difference between two times (t2 must be >= t1)
  INTEGER(i8) FUNCTION CalendarTimeDiff(t1,t2)
  INTEGER, INTENT(IN) :: t1(4), t2(4)

  INTEGER tsum, i

    tsum = 0
    IF ((t1(2) /= t2(2)) .OR. (t1(1) /= t2(1))) THEN
      ! Uncomment four lines below if i8 has been mapped to a 32-bit integer
      !  IF ((t2(1) - t1(1) > 67) .AND. (REAL(HUGE(tsum),dp) < 3e9)) THEN
      !    WRITE(*,*) t1,t2
      !    CALL local_error('CalendarTimeDiff','Time interval too large for integer type - please change to 64 bit integers')
      !  ENDIF
      IF (calendar > INPUT_CALENDAR_JULIAN) THEN
        tsum = FixedYearLen * (t2(1)-t1(1)-2)
      ELSE
        DO i=t1(1)+1,t2(1)-1
          tsum = tsum + CalendarYearLen(i)
        ENDDO
      ENDIF
      IF (t1(1) /= t2(1)) THEN
        DO i=t1(2),12
          tsum = tsum + CalendarMonthLen(t1(1),i)
        ENDDO
        DO i=1,t2(2)-1
          tsum = tsum + CalendarMonthLen(t2(1),i)
        ENDDO
      ELSE
        DO i=t1(2),t2(2)-1
          tsum = tsum + CalendarMonthLen(t2(1),i)
        ENDDO
      ENDIF
    ENDIF
    CalendarTimeDiff = INT((tsum + t2(3) - t1(3)),i8) * 86400 + t2(4) - t1(4)

  END FUNCTION CalendarTimeDiff

! Signed time difference between two times (t1-t2) in seconds
  INTEGER(i8) FUNCTION CalendarTimeDiffSec(t1,t2)
  INTEGER, INTENT(IN) :: t1(4), t2(4)

  INTEGER :: t3(4), t4(4)
  INTEGER (i8) :: dt
  LOGICAL :: sgn

    IF (CalendarTimeCmp(t1,t2) < 0) THEN
      sgn = .TRUE.
      t3  = t2
      t4  = t1
    ELSE
      sgn = .FALSE.
      t3  = t1
      t4  = t2
    ENDIF
    dt = CalendarTimeDiff(t3,t4)
    IF (sgn) dt = -dt
    CalendarTimeDiffSec = dt

  END FUNCTION CalendarTimeDiffSec

! Add time and time difference, +ve dt is in seconds, -ve in month
  SUBROUTINE CalendarTimeAdd(t,dt)
  INTEGER,   INTENT(INOUT) :: t(4)
  INTEGER(i8), INTENT(IN) :: dt

  INTEGER(i8) :: tmp, tmp2

    IF (dt < 0) THEN
      IF (MOD(-dt,12_i8)==0) THEN
        t(1) = t(1) - INT(dt/12_i8)
        RETURN
      ENDIF
      tmp  = t(1)
      tmp2 = t(2)
      tmp2 = tmp2 - dt
      t(1) = INT(tmp  + (tmp2-1) / 12)
      t(2) = INT(MOD(tmp2-1,INT(12,i8))) + 1
      tmp = CalendarMonthLen(t(1),t(2))
      IF (t(3) > tmp) t(3) = tmp
      RETURN
    ENDIF
    IF (dt <= GREGORIAN_MIN_MONTH) THEN ! Cannot step more than one month ahead
      t(4) = t(4) + dt
      t(3) = t(3) + t(4)/86400 ! Integer arithmetics secures correct result
      t(4) = MOD(t(4),86400)
      tmp = CalendarMonthLen(t(1),t(2))
      IF (t(3) > tmp) THEN
        t(3) = t(3) - tmp
        t(2) = t(2) + 1
        IF (t(2) > 12) THEN
          t(1) = t(1) + 1
          t(2) = 1
        ENDIF
      ENDIF
    ELSE
      ! TODO: This can be optimized for fixed length years
      tmp = dt + CalendarTimeDiff((/t(1),1,1,0/),t)
      tmp2 = CalendarYearLen(t(1)) * 86400
      DO WHILE (tmp - tmp2 >= 0)
        tmp = tmp - tmp2
        t(1) = t(1) + 1
        tmp2 = CalendarYearLen(t(1)) * 86400
      ENDDO
      t(2) = 1
      tmp2 = CalendarMonthLen(t(1),t(2)) * 86400
      DO WHILE (tmp - tmp2 >= 0)
        tmp = tmp - tmp2
        t(2) = t(2) + 1
        tmp2 = CalendarMonthLen(t(1),t(2)) * 86400
      ENDDO
      t(3) =     tmp/86400 + 1
      t(4) = MOD(tmp,86400_i8)
    ENDIF

  END SUBROUTINE CalendarTimeAdd

! Substracts a +ve time period from a time
  SUBROUTINE CalendarTimeSub(t,dt)
  INTEGER,   INTENT(INOUT) :: t(4)
  INTEGER(i8), INTENT(IN) :: dt

  INTEGER(i8) :: tmp, tmp2, tmp3, tmp4
  REAL(dp) :: mf

    tmp = dt
    IF (dt < 0) THEN
      IF (MOD(-dt,12_i8) == 0) THEN
        t(1) = t(1) + INT(dt/12_i8)
        RETURN
      ENDIF
      tmp3 = CalendarMonthLen(t(1),t(2))
      mf = REAL(t(3),dp)/REAL(tmp3,dp)
      tmp2 = CalendarMonthLen(t(1),t(2))-t(3)
      t(1) = t(1) + tmp/12
      tmp = MOD(-tmp,12_i8)
      t(2) = t(2) - tmp
      IF (t(2) <= 0) THEN
        t(1) = t(1) - 1
        t(2) = t(2) + 12
      ENDIF
      tmp = CalendarMonthLen(t(1),t(2))
      t(3) = INT(mf*REAL(tmp,dp)+.5_dp)
    ELSE
      IF (dt > GREGORIAN_MIN_MONTH) THEN ! Step more than one month backwards
        ! TODO: This can be optimized for fixed length years
        tmp2 = CalendarYearLen(t(1)-1)
        DO WHILE (tmp2*86400 < tmp)
          t(1) = t(1) - 1 
          tmp  = tmp  - tmp2 * 86400
          tmp2 = CalendarYearLen(t(1)-1)
        ENDDO
        tmp3 = t(2) - 1
        tmp4 = t(1)
        IF (tmp3 == 0) THEN
          tmp4 = tmp4 - 1
          tmp3 = 12
        ENDIF
        tmp2 = CalendarMonthLen(INT(tmp4),INT(tmp3))
        DO WHILE (tmp2*86400 < tmp)
          t(2) = t(2) - 1
          IF (t(2) == 0) THEN
            t(2) = 12
            t(1) = t(1) - 1
          ENDIF
          tmp = tmp - tmp2 * 86400
          tmp4 = t(1)
          tmp3 = t(2) - 1
          IF (tmp3 == 0) THEN
            tmp4 = tmp4 - 1
            tmp3 = 12
          ENDIF
          tmp2 = CalendarMonthLen(INT(tmp4),INT(tmp3))
        ENDDO
      ENDIF
      t(4) = t(4) - tmp
      IF (t(4) < 0) THEN
        t(3) = t(3) + (t(4)+1)/86400 - 1 ! t(4) will be < 86400 and probably < 0
        IF (t(4) < 0) t(4) = t(4) + 86400*CEILING(-REAL(t(4),dp)/86400)
        t(4) = MOD(t(4),86400)
        IF (t(3) <= 0) THEN
          t(2) = t(2) - 1
          IF (t(2)==0) THEN
            t(2) = 12
           t(1) = t(1) - 1
          ENDIF
          t(3) = t(3) + CalendarMonthLen(t(1),t(2))
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE CalendarTimeSub

! Add time and (signed) time difference in seconds
  SUBROUTINE CalendarTimeAddSec(t,dt)
  INTEGER,   INTENT(INOUT) :: t(4)
  INTEGER(i8), INTENT(IN) :: dt

    IF (dt < 0) THEN
      CALL CalendarTimeSub(t,-dt)
    ELSE
      CALL CalendarTimeAdd(t,dt)
    ENDIF

  END SUBROUTINE CalendarTimeAddSec

! Subtract time and (signed) time difference in seconds
  SUBROUTINE CalendarTimeSubSec(t,dt)
  INTEGER,   INTENT(INOUT) :: t(4)
  INTEGER(i8), INTENT(IN) :: dt

    IF (dt < 0) THEN
      CALL CalendarTimeAdd(t,-dt)
    ELSE
      CALL CalendarTimeSub(t,dt)
    ENDIF

  END SUBROUTINE CalendarTimeSubSec

! Returns which of two time are largest (-ve for t1 > t2)
  INTEGER FUNCTION CalendarTimeCmp(t1,t2)
  INTEGER, INTENT(IN) :: t1(4), t2(4)

  INTEGER i

    DO i=1,4
      IF (t1(i) /= t2(i)) THEN
        CalendarTimeCmp = t2(i) - t1(i)
        RETURN
      ENDIF
    ENDDO
    CalendarTimeCmp = 0

  END FUNCTION CalendarTimeCmp

END MODULE mo_input_calendar
