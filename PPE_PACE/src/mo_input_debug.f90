!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! Helper module for mo_input containing procedures dealing with the debug output
! Never use this file directly from outside the mo_input framework - only via mo_input.f90
MODULE mo_input_debug
USE mo_input_strings
USE mo_input_types
IMPLICIT NONE

PUBLIC InputIsDebug, InputDebugLevelSet, MessageDims

CONTAINS

! Tests if this PE should produce debug output
  LOGICAL FUNCTION InputIsDebug()

  LOGICAL ldbg

    SELECT CASE (IAND(mo_debug,INPUT_DBG_DBG_IODATA+INPUT_DBG_DBG_ALL))
      CASE (0)
        ldbg = lio
      CASE (INPUT_DBG_DBG_IODATA)
        ldbg = liodata
      CASE (INPUT_DBG_DBG_ALL)
        ldbg = .TRUE.
    END SELECT
    InputIsDebug = ldbg

  END FUNCTION InputIsDebug

! Set which information is posted to the log
  SUBROUTINE InputDebugLevelSet(Debug,sDebug,Flag)
  INTEGER,          INTENT(in), OPTIONAL :: Debug
  CHARACTER(len=*), INTENT(in), OPTIONAL :: sDebug
  LOGICAL,          INTENT(in), OPTIONAL :: Flag

  CHARACTER (len= 10) :: str
  CHARACTER (len=128) :: fname
  INTEGER :: Dbg, St, En, iLen, i, j, k, m, Sub, Unt

    Dbg = INPUT_DBG_DEFAULT 
    IF (PRESENT(Debug)) Dbg = Debug

    ! Process flags when passed as string
    IF (PRESENT(sDebug)) THEN
      Dbg = 0
      iLen = LEN_TRIM(sDebug)
      St = 1
      DO WHILE (St <= iLen)
        En = INDEX(sDebug(St:),',') + St - 1
        IF (En < St) En = iLen+1
        IF (En - St > 10) THEN
          str = tolower(sDebug(St:St+9))
          IF (str(1:10)=='input_dbg_') St = St + 10
        ENDIF
        ! Disentangle log file setup for parallel runs
        Sub = St
        DO WHILE ((Sub < En-1) .AND. (sDebug(Sub:Sub) < '0' .OR. sDebug(Sub:Sub) > '9'))
          Sub = Sub + 1
        ENDDO
        IF (Sub < En-1) THEN
          i = Sub
          Unt = 0
          DO WHILE ((Sub < En-1) .AND. (sDebug(i:i) >= '0') .AND. (sDebug(i:i) <= '9'))
            Unt = Unt * 10 + ICHAR(sDebug(i:i)) - ICHAR('0')
            i = i + 1
          ENDDO
          IF (sDebug(i:i)=='+') THEN
            i = i + 1
            Unt = Unt + my_pe
          ENDIF
          fname = ''
          IF ((i < En-1) .AND. (sDebug(i:i) /= ',')) THEN
            j = INDEX(sDebug(i:En-1),'<') + i - 1
            IF (j >= i) THEN
              k = INDEX(sDebug(j:En-1),'>') + j - 1
              IF (k < j) CALL local_error('InputDebugLevelSet','Parameter format error: '//sDebug(St:En-1))
              fname = sDebug(i:j-1)
              IF (sDebug(j+1:j+2)=='PE') THEN
                m = my_pe
              ELSE
                m = Unt
              ENDIF
              IF ((sDebug(k-1:k-1) >= '0') .AND. (sDebug(k-1:k-1) <= '9')) THEN
                str = '(a,i'//sDebug(k-1:k-1)//'.'//sDebug(k-1:k-1)//')'
              ELSE
                str = '(a,i5)'
              ENDIF
              WRITE(fname,TRIM(str)) sDebug(i:j-1),m
              IF (k+1 < En) WRITE(fname(LEN_TRIM(fname)+1:),*) sDebug(k+1:En-1)
              CALL RmSpc(fname,';')
            ELSE
              fname = sDebug(i:En-1)
            ENDIF
          ENDIF
          Sub = Sub - 1
        ENDIF
        ! Now get the main option
        i = GetStringIndex(sDebug(St:Sub),flg_dbg,.FALSE.,' ')
        SELECT CASE (i)
          CASE (-1: 0)  ! Bad option, just ignore it
            CALL local_message('InputDebugLevelSet','Option: '//sDebug(St:En-1)//' unknown. Ignored')
          CASE ( 1: 6)  ! Time output format option
            Dbg = IAND(Dbg,INPUT_ALLFLAGS-INPUT_DBG_TIME_MASK) + i-1
          CASE ( 7:11)  ! Avg, Extern and IOPe/All flags - independent bit flags
            Dbg = IOR(Dbg,ISHFT(INPUT_DBG_AVG,i-7))
            IF ((i>9) .AND. (Sub /= En)) THEN
              IF (LEN_TRIM(fname) > 0) THEN
                OPEN(unit=Unt,file=TRIM(fname),action='write',form='formatted')
                LogOpened = .TRUE.
              ENDIF
              CALL SetLogUnit(Unt,Flush=IAND(Dbg,INPUT_DBG_DBG_FLUSH) /= 0)
            ENDIF
          CASE (12:21) ! Output level flags
            Dbg = IAND(Dbg,INPUT_ALLFLAGS-INPUT_DBG_MASK) + INPUT_DBG_WARN * (i-12)
        END SELECT
        St = En + 1
      ENDDO
    ENDIF

    IF (IAND(Dbg,INPUT_DBG_MASK) == 0) THEN
      mo_debug = IOR(IAND(mo_debug,INPUT_ALLFLAGS - INPUT_DBG_MASK),Dbg)
    ELSE
      IF (IAND(Dbg,INPUT_DBG_TIME_MASK) == 0) THEN
        mo_debug = IOR(IAND(mo_debug,INPUT_ALLFLAGS - INPUT_DBG_TIME_MASK),Dbg)
      ELSE
        mo_debug = IAND(Dbg,INPUT_DBG_TIME_MASK+INPUT_DBG_MASK+INPUT_DBG_AVG+INPUT_DBG_EXTERN)
      ENDIF
    ENDIF
    IF (PRESENT(Flag)) THEN ! Set all regardless of masks
      IF (Flag) mo_debug = IAND(Dbg, &
        INPUT_DBG_TIME_MASK+INPUT_DBG_MASK+INPUT_DBG_AVG+INPUT_DBG_EXTERN+INPUT_DBG_DBG_IODATA+INPUT_DBG_DBG_ALL)
    ENDIF
    IF (IAND(Dbg,INPUT_DBG_TIME_MASK) == INPUT_DBG_TIME_MASK) &
      mo_debug = IAND(mo_debug,INPUT_ALLFLAGS - INPUT_DBG_TIME_MASK) + INPUT_DBG_TIME_YMDH
    IF ((IAND(mo_debug,INPUT_DBG_MASK) < INPUT_DBG_READ_VAR) .AND. (IAND(mo_debug,INPUT_DBG_READ_GROUP)==0)) & 
      mo_debug = IAND(mo_debug,INPUT_ALLFLAGS - INPUT_DBG_AVG)

    ldbg = InputIsDebug()

  END SUBROUTINE InputDebugLevelSet

! Write messages about expected and present dimentions of a variable
  SUBROUTINE MessageDims(var,dims,alt_msg)
  TYPE (input_var_list), POINTER :: var
  TYPE (input_dim_list), POINTER :: dims
  CHARACTER (len=*), OPTIONAL, INTENT(in) :: alt_msg

  TYPE (input_dim_list), POINTER :: curr
  CHARACTER (len=256) :: my_msg

    IF (ldbg) THEN ! Output only from debug PEs
      my_msg = TRIM(var%name_var)//':%Expected_dimensions:%'
      curr => var%dims
      DO WHILE (ASSOCIATED(curr))
        WRITE(my_msg(LEN_TRIM(my_msg)+1:),'(2a)') TRIM(curr%dim_data%name_dim),','
        curr => curr%next
      ENDDO
      IF (PRESENT(alt_msg)) THEN
        WRITE(my_msg(LEN_TRIM(my_msg):),'(3a)') '.%',TRIM(alt_msg),':%'
      ELSE
        WRITE(my_msg(LEN_TRIM(my_msg):),'(a)') '.%Dimensions%found%in%file:%'
      ENDIF
      curr => dims
      DO WHILE (ASSOCIATED(curr))
        WRITE(my_msg(LEN_TRIM(my_msg)+1:),'(2a)') TRIM(curr%dim_data%name_dim),','
        curr => curr%next
      ENDDO
      my_msg(LEN_TRIM(my_msg):LEN_TRIM(my_msg)) = ' '
      CALL RmSpc(my_msg,'%')
      CALL local_message('WARNING',TRIM(my_msg))
    ENDIF

  END SUBROUTINE MessageDims

END MODULE mo_input_debug
