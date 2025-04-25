!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_util_logging

  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_CHAR, C_NULL_CHAR

  IMPLICIT NONE

  PRIVATE

  INTERFACE 
    FUNCTION open_network_logging(ip_address, port) RESULT(iret) BIND(C,NAME='open_network_logging')
      IMPORT :: C_INT, C_CHAR
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR) :: ip_address(*)
      INTEGER(C_INT), VALUE :: port
    END FUNCTION open_network_logging
  END INTERFACE

  INTERFACE
    FUNCTION send_network_log_message(message) RESULT(iret) BIND(C,NAME='send_network_log_message')
      IMPORT :: C_INT, C_CHAR
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR) :: message(*)
    END FUNCTION send_network_log_message
  END INTERFACE

  INTERFACE
    SUBROUTINE close_network_logging() BIND(C,NAME='close_network_logging')
    END SUBROUTINE close_network_logging
  END INTERFACE

  LOGICAL :: connected = .FALSE.

  PUBLIC :: open_network_log
  PUBLIC :: send_network_log
  PUBLIC :: close_network_log

CONTAINS

  SUBROUTINE open_network_log(hostname, port)
    CHARACTER(len=*,kind=C_CHAR), INTENT(in) :: hostname
    INTEGER,          INTENT(in) :: port
    INTEGER :: iret
    iret = open_network_logging(TRIM(hostname)//C_NULL_CHAR, port)
    IF (iret == 0) THEN
      connected = .TRUE.
    ENDIF
  END SUBROUTINE open_network_log

  SUBROUTINE send_network_log(message)
    CHARACTER(len=*,kind=C_CHAR), INTENT(in) :: message
    INTEGER :: iret
    IF (connected) THEN
      iret = send_network_log_message(TRIM(message)//C_NULL_CHAR)
    ENDIF
  END SUBROUTINE send_network_log

  SUBROUTINE close_network_log
    IF (connected) THEN
      CALL close_network_logging
    ENDIF
  END SUBROUTINE close_network_log

END MODULE mo_util_logging
