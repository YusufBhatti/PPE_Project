MODULE mo_util_string

CONTAINS

  FUNCTION tolower (upper)
    !-----------------------------------
    ! Conversion: Uppercase -> Lowercase
    !-----------------------------------
    CHARACTER(len=*)              ,INTENT(in) :: upper
    CHARACTER(len=LEN_TRIM(upper))            :: tolower
    
    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')
    
    DO i=1,LEN_TRIM(upper)
      IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
           ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
      ELSE
        tolower(i:i) = upper(i:i)
      END IF
    END DO
    
  END FUNCTION tolower

END MODULE mo_util_string

