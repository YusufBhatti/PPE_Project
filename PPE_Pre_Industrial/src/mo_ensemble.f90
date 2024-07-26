#if (__INTEL_COMPILER < 1200) || defined (__SUNPRO_F95) || defined (__PGI) || (__IBMC__ < 1210)
#define _FORTRAN_STANDARD 2003
#else
#define _FORTRAN_STANDARD 2008         
#endif

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_ensemble

  IMPLICIT NONE

  PRIVATE

  ! ___________________________________________________________________________
  ! 
  ! The number range is defined by:
  ! WMO, Manual of Codes, Vol. I.2, FM 92-XIV GRIB, Code table 4.6.
  !
  ! This range cannot be changed as this information needs to be 
  ! stored in one octet (8-bit).
  !
  ! unperturbed high-resolution control forecast
  INTEGER, PARAMETER, PUBLIC :: HR_CONTROL        =   0 
  ! unperturbed low-resolution control forecast
  INTEGER, PARAMETER, PUBLIC :: LR_CONTROL        =   1
  ! negatively perturbed forecast
  INTEGER, PARAMETER, PUBLIC :: NEGATIV_PERTURBED =   2
  ! positively perturbed forecast
  INTEGER, PARAMETER, PUBLIC :: POSITIV_PERTURBED =   3
  ! multi-model forecast
  INTEGER, PARAMETER, PUBLIC :: MULTI_MODEL       =   4
  !   5 - 191 Reserved
  ! 192 - 254 Reserved for local use
  ! missing
  INTEGER, PARAMETER, PUBLIC :: MISSING           = 255
  !___________________________________________________________________________

  INTEGER, PARAMETER :: UNDEFINED = -1

  INTEGER :: ensemble_member  = UNDEFINED
  INTEGER :: ensemble_size    = UNDEFINED
  INTEGER :: forecast_type    = MISSING

  CHARACTER(len=8) :: ensemble_suffix = ''

  PUBLIC :: init_ensemble
  PUBLIC :: finalize_ensemble
  PUBLIC :: get_ensemble_member
  PUBLIC :: get_ensemble_size
  PUBLIC :: get_forecast_type
  PUBLIC :: get_ensemble_suffix
  PUBLIC :: is_ensemble

CONTAINS
  
  SUBROUTINE init_ensemble(finish, broadcast, root_rank)
    INTERFACE
      SUBROUTINE finish(routine, message, exit_no)
        CHARACTER(len=*), INTENT(in) :: routine
        CHARACTER(len=*), INTENT(in), OPTIONAL :: message
        INTEGER,          INTENT(in), OPTIONAL :: exit_no
      END SUBROUTINE finish
    END INTERFACE

    OPTIONAL :: broadcast
    INTERFACE
      ! broadcast: wrapper function for a parallel broadcast
      SUBROUTINE broadcast(buffer, root, comm)
        ! buffer: assumed-size array for sending arbitrary integer arrays
        integer, intent(inout) :: buffer(:) 
        ! root: the rank from which the data are distributed to all others
        !       optional, if not given assume rank 0  
        integer, intent(in) :: root
        ! comm: a communicator for the broadcast
        integer, intent(in), optional :: comm
      END SUBROUTINE broadcast
    END INTERFACE

    INTEGER,  INTENT(in), OPTIONAL :: root_rank

    INTEGER :: broadcast_buffer(3)
    INTEGER :: funit

    LOGICAL :: lexist

    NAMELIST /ensctl/ ensemble_member, ensemble_size, forecast_type
    
    ! check if ensemble is already initialized
    IF (is_ensemble()) RETURN
    
    lexist = .FALSE.
    INQUIRE(FILE='namelist.ens', EXIST=lexist)

    IF (lexist) THEN
      funit = UNDEFINED
#if _FORTRAN_STANDARD < 2008
      funit=newunit(finish)
      OPEN(UNIT=funit, FILE='namelist.ens', FORM='formatted', STATUS='old')
#else
      OPEN(NEWUNIT=funit, FILE='namelist.ens', FORM='formatted', STATUS='old')
#endif
      READ(funit, ensctl)
      CLOSE(funit)
      IF (ensemble_member < 1 .OR. ensemble_member > ensemble_size) THEN
        CALL finish('','Ensemble member is not in valid range.')
      ENDIF
      IF (forecast_type < 0 .OR. forecast_type > 255) THEN
        CALL finish('','Forecast type is not defined.')
      ENDIF
      IF (PRESENT(broadcast)) THEN
        broadcast_buffer(1) = ensemble_member
        broadcast_buffer(2) = ensemble_size
        broadcast_buffer(3) = forecast_type 
        CALL broadcast(broadcast_buffer, root_rank)
        ensemble_member = broadcast_buffer(1)
        ensemble_size   = broadcast_buffer(2)
        forecast_type   = broadcast_buffer(3) 
      ENDIF
      WRITE(ensemble_suffix,'(a4,i4.4)') '.ens', ensemble_member
    ENDIF
     
  END SUBROUTINE init_ensemble

  SUBROUTINE finalize_ensemble()
    ensemble_member = UNDEFINED
    ensemble_size   = UNDEFINED
    forecast_type   = MISSING
    ensemble_suffix = ''
  END SUBROUTINE finalize_ensemble

  INTEGER FUNCTION get_ensemble_member()
    get_ensemble_member = ensemble_member
  END FUNCTION get_ensemble_member

  INTEGER FUNCTION get_ensemble_size()
    get_ensemble_size = ensemble_size
  END FUNCTION get_ensemble_size

  INTEGER FUNCTION get_forecast_type()
    get_forecast_type = forecast_type
  END FUNCTION get_forecast_type

  FUNCTION get_ensemble_suffix() RESULT(filename_suffix)
    CHARACTER(len=8) :: filename_suffix
    filename_suffix = ensemble_suffix
  END FUNCTION get_ensemble_suffix
  
  LOGICAL FUNCTION is_ensemble()
    is_ensemble = (ensemble_member > 0)
  END FUNCTION is_ensemble

#if _FORTRAN_STANDARD < 2008
  INTEGER FUNCTION newunit(finish)
    INTERFACE
      SUBROUTINE finish(routine, message, exit_no)
        CHARACTER(len=*), INTENT(in) :: routine
        CHARACTER(len=*), INTENT(in), OPTIONAL :: message
        INTEGER,          INTENT(in), OPTIONAL :: exit_no
      END SUBROUTINE finish
    END INTERFACE
    INTEGER :: iunit
    LOGICAL :: lopen
    lopen = .FALSE.
    DO iunit = 10, 1000
      INQUIRE(UNIT=iunit, OPENED=lopen)
      IF (.NOT. lopen) THEN
        newunit = iunit
        RETURN
      ENDIF
    ENDDO
    CALL finish('','No free unit available for namelist reading.')
  END FUNCTION newunit
#endif

END MODULE mo_ensemble
