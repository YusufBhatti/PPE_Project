!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! Helper module for mo_input containing simple procedures which are called at rather low level
! Never use this file directly from outside the mo_input framework - only via mo_input.f90
MODULE mo_input_base
USE mo_input_calendar
USE mo_input_strings
USE mo_input_types
USE mo_input_debug
USE mo_input_mask
IMPLICIT NONE

PUBLIC :: GetRecTime, GetTimeScaleFromStr
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
PUBLIC :: InputParSetSourcePE, InputVarBroadcast, InputVarScatter
#endif

CONTAINS
! ---------------------------------------------------------------------
!  Processing of time strings to convert numeric date/time info to 
! model time. These functions are constructed to process typical 
! NetCDF time information, but do no contain any NetCDF library calls.
! ---------------------------------------------------------------------

! Convert a specific time given in file units to model time
  FUNCTION GetRecTime(IFile,Time)
  INTEGER :: GetRecTime(4)
  TYPE (input_file_list), POINTER :: IFile
  REAL(dp),            INTENT(in) :: Time

  CHARACTER (len=*), PARAMETER :: SpecStrs   = 'sMhdmY'
  INTEGER,           PARAMETER :: SpecLen(6) = (/2,2,2,2,2,4/)
  INTEGER,           PARAMETER :: UnitMul(6) = (/1,60,3600,86400,-1,-12/)

  INTEGER  :: Ofs(6), CurrPos, iLen, i, j
  REAL(dp) :: Tme
  CHARACTER(len=64) :: Units

    IF (IFile%time_scale(2) < 0) THEN ! CDO "absolute" time
      Units = ' '//IFile%time_unit(-IFile%time_scale(2):)
      Ofs = 0
      Ofs(2:3) = 1
      i = 0
      iLen = LEN_TRIM(Units)
      j = INDEX(Units,'.')
      IF ((j <= iLen) .AND. (j>0)) THEN
        CurrPos = INDEX(Units,'.') - 1
      ELSE
        CurrPos = iLen
      ENDIF
      DO WHILE ((CurrPos > 1) .AND.Units(CurrPos:CurrPos) /= ' ')
        j = INDEX(SpecStrs,Units(CurrPos:CurrPos))
        CurrPos = CurrPos - 1
        IF (Units(CurrPos:CurrPos) == '%') CurrPos = CurrPos - 1
        Ofs(7-j) = MOD(FLOOR(Time/(10._dp**i)),NINT(10._dp**SpecLen(j)))
        i = i + SpecLen(j)
      ENDDO
      ! Process fraction if given (fraction of a second is ignored)
      i = INDEX(Units,'.')
      IF ((i < iLen) .AND. (i > 0)) THEN
        IF (Units(i-1:i-1) /= 's') THEN
          Tme = MOD(Time,1._dp)
          j = INDEX(SpecStrs,Units(i-1:i-1))
          SELECT CASE (j)
            CASE (2,3,4) ! fraction of {minute, hour, day}
              Ofs(4) = Ofs(4)+NINT(Tme * UnitMul(j))
            CASE (5) ! fraction of month
              CALL CalendarTimeAdd(Ofs(1:4),NINT(Tme*CalendarMonthLen(Ofs(1),Ofs(2))*86400,i8))
            CASE (6) ! fraction of year
              CALL CalendarTimeAdd(Ofs(1:4),NINT(Tme*CalendarYearLen(Ofs(1))*86400,i8))
          END SELECT
        ENDIF
      ENDIF
    ELSE ! NetCDF std time (CDO "relative" time)
      Ofs(1:4) = IFile%time_scale(1:4)
      IF (Time < 0._dp) THEN
        CALL CalendarTimeSub(Ofs(1:4),NINT(-Time*IFile%time_scale(5),i8))
      ELSE
        CALL CalendarTimeAdd(Ofs(1:4),NINT(Time*IFile%time_scale(5),i8))
      ENDIF
    ENDIF
    CALL CalendarTimeAddSec(Ofs(1:4),INT(IFile%time_offset(2)*86400+IFile%time_offset(3),i8))
    GetRecTime = Ofs(1:4)

  END FUNCTION GetRecTime

! Interpret time description string in the formats known by the CDOs
  SUBROUTINE GetTimeScaleFromStr(IFile,Units)
  TYPE(input_file_list), POINTER :: IFile
  CHARACTER(len=*),   INTENT(in) :: Units

  INTEGER :: UnitMul, CurrPos, NextPos, iLen, i, j, Ofs(6)
  LOGICAL :: lOfs

    IFile%time_unit = Units

    ! Determine unit of time
    i       = INDEX(Units,'_') + 1 ! Ignore until first '_' for compatibility with udunits units like 'Gregorian_year'
    NextPos = INDEX(Units,' ')
    IF ((i <= 1) .OR. (i >= NextPos)) i = 1
    UnitMul = InterpretDtUnit(1,Units(i:NextPos-1))
    IF (UnitMul == 0) RETURN

    ! Determining type of time unit
    CurrPos = 1
    lOfs = .TRUE.
    iLen = LEN_TRIM(Units)
    IF (CurrPos <= iLen) THEN
      CurrPos = NextPos+1
      NextPos = INDEX(Units(CurrPos:),' ') + NextPos
      IF (Units(CurrPos:NextPos) == 'as') lOfs = .FALSE.
    ENDIF

    IF (lOfs) THEN ! Offset of time (std. NetCDF time)
      IFile%time_scale(5) = UnitMul
      IF (CurrPos <= iLen) THEN ! Classical NetCDF format time
        CurrPos = NextPos + 1
        i = 1
        Ofs(:) = 0
        DO WHILE ((CurrPos <= iLen) .AND. (i <= 6))
          j = ICHAR(Units(CurrPos:CurrPos))
          IF ((j < 48) .OR. (j > 57)) THEN
            i = i + 1
          ELSE
            Ofs(i) = 10*Ofs(i) + j - 48
          ENDIF
          CurrPos = CurrPos + 1
        ENDDO
        IFile%time_scale(1:3) = Ofs(1:3)
        IFile%time_scale(4)   = Ofs(4)*3600 + Ofs(5)*60 + Ofs(6)
      ELSE
        IF (UnitMul==-12) THEN
          IFile%time_scale(2:3) = 1 ! Account for non-standard but common unit "year(s)"
        ELSE
          CALL local_error('GetTimeScaleFromStr','Cannot interpret unit of time dimension: '//TRIM(Units))
        ENDIF
      ENDIF
    ELSE ! Interpreter string (CDO "Absolute" time)
      ! For now just save the position where the important format information is
      ! Save negative value to distinguish from normal NetCDF time (CDO "relative" time)
      IFile%time_scale(2) = -NextPos-1
    ENDIF

  END SUBROUTINE GetTimeScaleFromStr

! -------------------------
!  Working with file times 
! -------------------------

! Get the start time of file fitting to a given time
  FUNCTION InputFileTimeCurrGet(IFile,ftime,def)
  INTEGER :: InputFileTimeCurrGet(4)
  TYPE (input_file_list), POINTER :: IFile
  INTEGER, INTENT(IN) :: ftime(4)
  INTEGER, INTENT(IN), OPTIONAL :: def(4)

  INTEGER :: res(4)

    res = ftime
    res(4) = 0
    IF (IAND(IFile%flags,INPUT_FILE_NAME_MONTH) == 0) & ! File name is constructed using month number
    res(2) = 1
    IF (IAND(IFile%flags,INPUT_FILE_NAME_DAY) == 0) & ! File name is constructed using day number
    res(3) = 1
    IF (IAND(IFile%flags,INPUT_FILE_MULTI) == 0) THEN
      IF (PRESENT(def)) THEN
        res = def
      ELSE
        res = exp_start
      ENDIF
    ENDIF

    InputFileTimeCurrGet = res

  END FUNCTION InputFileTimeCurrGet

! Time to expect the previous file name to refere to
  FUNCTION InputFileTimePrevGet(IFile,ftime)
  INTEGER :: InputFileTimePrevGet(4)
  TYPE (input_file_list), POINTER :: IFile
  INTEGER, INTENT(IN) :: ftime(4)

  INTEGER :: res(4)

    res = ftime
    res(4) = 0
    IF (IAND(IFile%flags,INPUT_FILE_NAME_DAY) /= 0) THEN ! File name is constructed using day number
      CALL CalendarTimeSub(res,86400_i8)
    ELSEIF (IAND(IFile%flags,INPUT_FILE_NAME_MONTH) /= 0) THEN ! File name is constructed using month number
      CALL CalendarTimeSub(res,-1_i8)
    ELSEIF (IAND(IFile%flags,INPUT_FILE_NAME_YEAR) /= 0) THEN ! File name is constructed using year num (4 or 6 digit (evt +sign))
      res(1) = res(1) - 1
    ENDIF

    InputFileTimePrevGet = res

  END FUNCTION InputFileTimePrevGet

! Time to expect the next file name to refere to
  FUNCTION InputFileTimeNextGet(IFile,ftime)
  INTEGER :: InputFileTimeNextGet(4)
  TYPE (input_file_list), POINTER :: IFile
  INTEGER, INTENT(IN) :: ftime(4)

  INTEGER :: res(4)

    res = ftime
    res(4) = 0
    IF (IAND(IFile%flags,INPUT_FILE_NAME_DAY) /= 0) THEN ! File name is constructed using day number
      res(3) = res(3) + 1
      IF (res(3) > CalendarMonthLen(res(1),res(2))) THEN
        res(3) = 1
        res(2) = res(2) + 1
        IF (res(2) > 12) THEN
          res(1) = res(1) + 1
          res(2) = 1
        ENDIF
      ENDIF
      ELSEIF (IAND(IFile%flags,INPUT_FILE_NAME_MONTH) /= 0) THEN ! File name is constructed using month number
      res(2) = res(2) + 1
      IF (res(2) > 12) THEN
        res(1) = res(1) + 1
        res(2) = 1
      ENDIF
      ELSEIF (IAND(IFile%flags,INPUT_FILE_NAME_YEAR) /= 0) THEN ! File name is constructed using year num (4 or 6 digit (evt +sign))
      res(1) = res(1) + 1
    ENDIF

    InputFileTimeNextGet = res

  END FUNCTION InputFileTimeNextGet

! ---------------------------------------------------------------------
!  Parallel functions. Could advantageously be outsourced to elsewhere
! ---------------------------------------------------------------------

#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
! Sets which PE should provide this PE with data
  SUBROUTINE InputParSetSourcePE(SrcPe,comm)
#if !defined(INPUT_IN_ECHAM) || !defined(NOMPI)
  USE mpi ! Unfortunately mo_mpi doesn't define a wrapper for MPI_AllGather, thus we need to include it directly from mpi
#endif
  INTEGER,           INTENT(in) :: SrcPe
  INTEGER, OPTIONAL, INTENT(in) :: comm

  LOGICAL, SAVE :: lIgnore = .FALSE.
  INTEGER :: icomm, ierr


    IF ((n_pe == 1) .OR. .NOT. ASSOCIATED(src_pe) .OR. lIgnore) RETURN ! Not parallel run or inproper initialization has been done
#if !defined(INPUT_IN_ECHAM) || !defined(NOMPI)
    icomm = MPI_COMM_WORLD
    IF (PRESENT(comm)) icomm = comm

    CALL MPI_AllGather(SrcPe,1,MPI_INT,src_pe,1,MPI_INT,icomm,ierr)
    IF (ierr /= MPI_SUCCESS) CALL local_error('InputParSetSourcePE','MPI_AllGather error')
#else
    src_pe(1) = SrcPe
#endif

    ierr = COUNT(src_pe==my_pe)
    liodata   = .FALSE.
    lioserver = .FALSE.
    IF (ierr > 0) liodata   = .TRUE.
    IF (ierr > 1) lioserver = .TRUE.
    ldbg = InputIsDebug() 

    IF (lio .AND. .NOT. liodata) CALL local_error('InputParSetSourcePE','The master IO PE must be an IOData-PE as well')
    IF (liodata .AND. (src_pe(my_pe+1) /= my_pe)) CALL local_error('InputParSetSourcePE','Any IOData-PE must be its own source PE')

    lIgnore = .TRUE.

  END SUBROUTINE InputParSetSourcePE

! Broadcasts a variable from each lioserver-pe to all its destination pes
  SUBROUTINE InputVarBroadcast(buf)
  USE mo_mpi, ONLY: p_recv, p_isend, p_wait
  TYPE (input_data_list), POINTER :: buf

  REAL(dp), POINTER :: buf_data(:,:,:,:)
  INTEGER :: i

    IF (liodata .AND. .NOT. lioserver) RETURN ! PE should only provide itself with data - no action required

    buf_data => buf%dta(:,:,:,:,1,1,1) ! Since the size is not checked by subroutines p_*, this is sufficient also for 5-7D data
    IF (liodata) THEN ! Send from this PE
      DO i=1,n_pe
        IF (src_pe(i) == my_pe) THEN
          IF (my_pe /= i-1) THEN
            CALL p_isend(buf_data,i-1,2,SIZE(buf%dta))
          ENDIF
        ENDIF
      ENDDO
      CALL p_wait()
    ELSE ! Recieve to this PE
      CALL p_recv(buf_data,src_pe(my_pe+1),2,SIZE(buf%dta))
    ENDIF

  END SUBROUTINE InputVarBroadcast

! Scatters a variable subject to masks
  FUNCTION InputVarScatter(mask,src,srcdta,dst,buf)
  USE mo_mpi, ONLY: p_recv, p_send ! p_isend, p_wait ! Should also work with non-blocking sends, but data recieved gets wrong!
  TYPE (input_mask_list), POINTER :: InputVarScatter
  TYPE (input_mask_list), POINTER :: mask
  TYPE (input_data_list), POINTER :: src, dst, buf
  REAL(dp),         INTENT(inout) :: srcdta(*)

  TYPE (input_mask_list), POINTER :: my_mask, curr_mask
  REAL(dp), POINTER :: src_data(:,:,:,:), dst_data(:,:,:,:)  
  INTEGER :: i, ofs

    dst_data => dst%dta(:,:,:,:,1,1,1) ! Since the size is not checked by subroutines p_*, this is sufficient also for 5-7D data
    curr_mask => mask
    IF (liodata) THEN ! Data read on this PE
      IF (ASSOCIATED(buf)) THEN
        src_data => buf%dta(:,:,:,:,1,1,1)        
        buf%dta=0._dp
      ELSE
        src_data => src%dta(:,:,:,:,1,1,1)
      ENDIF
      ofs = 1
      NULLIFY(my_mask)
      DO i=0,n_pe-1
        IF (src_pe(i+1) == my_pe) THEN
          IF (my_pe == i) THEN ! These are the data to "send" to the io-pe itself, save ref for data copy while waiting for sends
            my_mask => curr_mask
          ELSE
            IF (IAND(curr_mask%flags,INPUT_MASK_CONT) == 0) THEN ! Non-continuous data needs extra buffering
              CALL MaskPack(src%dta,buf%dta(ofs:,:,:,:,:,:,:),curr_mask%mask)
              CALL p_send(src_data(ofs:ofs+curr_mask%nvalid-1,:,:,:),i,1,curr_mask%nvalid)
              ofs = ofs + curr_mask%nvalid ! Wenn using isends, we're not allowed to reuse buffer before send is completed
            ELSE ! Continuous data can be send directly
              CALL p_send(srcdta(IAND(curr_mask%flags,INPUT_MASK_FIRST_VALID)),i,1,curr_mask%nvalid)
            ENDIF
          ENDIF
          curr_mask => curr_mask%next
        ENDIF
      ENDDO
      CALL MaskPack(src%dta,dst%dta,my_mask%mask)
!      CALL p_wait()
    ELSE ! Recieve to this PE
      CALL p_recv(dst_data,src_pe(my_pe+1),1,SIZE(dst%dta))
    ENDIF

    InputVarScatter => curr_mask

  END FUNCTION InputVarScatter
#endif

END MODULE mo_input_base
