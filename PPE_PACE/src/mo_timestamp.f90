!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_timestamp

  USE mo_exception, ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_timestamp

CONTAINS

  SUBROUTINE write_timestamp(application)
    
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: application

    CHARACTER(len= 8) :: date
    CHARACTER(len=10) :: time
    CHARACTER(len=10) :: iso_date
    CHARACTER(len=12) :: new_time

    CALL date_and_time(date, time)

    ! reformat date to ISO-8601

    WRITE(iso_date,'(a4,a1,a2,a1,a2)') &
         date(1:4), '-', date(5:6), '-', date(7:8)

    ! reformat time string

    WRITE (new_time,'(a2,a1,a2,a1,a2,a1,a3)') &
         time(1:2), ':', time(3:4), ':', time(5:6), '.', time(8:10)

    ! write out ...

    IF (PRESENT(application)) THEN
      WRITE (message_text,'(5(a))') TRIM(application), ': ', iso_date, ' ',new_time
    ELSE
      WRITE (message_text,'(3(a))') iso_date, ' ',new_time
    ENDIF

    CALL message('',message_text)

  END SUBROUTINE write_timestamp

END MODULE mo_timestamp
