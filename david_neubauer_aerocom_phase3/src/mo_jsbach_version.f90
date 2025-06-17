#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! ----------------------------------------------------------------------------------------------
!!
!! Handle JSBACH version
!!
!! Define and print JSBACH version number and label a simulation.
!!
!! @author Reiner Schnur, MPI Met, Hamburg
!!
!! @par Revision History
!! Initial version derived from ECHAM6 mo_version.f90 by Reiner Schnur, MPI Met (2011-02-10)
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_jsbach_version

  USE mo_exception,         ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: jsbach_version, jsbach_init_version, jsbach_label_run

  CHARACTER(len=  4), PARAMETER :: jsbach_version = '3.10' !< JSBACH version number

  CHARACTER(len=256) :: executable                         !< Name of executable
  CHARACTER(len= 80) :: label(4)                           !< Space for storing model information

CONTAINS

  !>
  !! Define JSBACH version and labels
  !!
  !! Define JSBACH version number and get some information about a simulation
  !! at model initialization. 
  !!
  !! @par Revision History
  !! Initial version derived from ECHAM6 mo_version::init_version by Reiner Schnur, MPI Met (2011-02-10)
  !!
  SUBROUTINE jsbach_init_version(standalone)

    USE mo_netcdf,       ONLY: global_att, put_att
    USE mo_util_vcs,     ONLY: util_repository_url, util_branch_name, util_revision_key

    LOGICAL, INTENT(in), OPTIONAL :: standalone !< true: offline JSBACH, false: ECHAM/JSBACH

    LOGICAL             :: l_standalone !< true: offline JSBACH, false: ECHAM/JSBACH
    CHARACTER (len=  8) :: ydate        !< Date of simulation
    CHARACTER (len= 10) :: ytime        !< Time of simulation
    CHARACTER (len=256) :: os_name      !< Name of OS for simulation
    CHARACTER (len=256) :: user_name    !< User that runs simulation
    CHARACTER (len=256) :: host_name    !< Host name for simulation
    CHARACTER (len=256) :: repository = ''
    CHARACTER (len=256) :: branch     = ''
    CHARACTER (len=256) :: revision   = ''

    INTEGER :: i_lena, i_lenb, i_lenc
    INTEGER :: i
    INTEGER :: nlen

    !  External subroutines
    EXTERNAL :: util_os_system, util_user_name, util_node_name

    !  Executable statements

    l_standalone = .TRUE.
    IF (PRESENT(standalone)) l_standalone = standalone

    IF (l_standalone) THEN
      ! Print version

      nlen = 256
      CALL util_repository_url(repository, nlen)
      nlen = 256
      CALL util_branch_name(branch, nlen)
      nlen = 256
      CALL util_revision_key(revision, nlen)

      CALL message ('','')
      CALL message ('','')
      message_text = '==========================================================='
      CALL message ('',TRIM(message_text))        
      CALL message ('','')
  
      message_text = '  JSBACH - Version '//jsbach_version
      CALL message ('',TRIM(message_text))        
      message_text = '  Copyright by Max-Planck-Institute for Meteorology, 2010'
      CALL message ('',TRIM(message_text))        
      WRITE(message_text,'(a,a)') 'Repository: ', TRIM(repository)
      CALL message('',message_text)
      WRITE(message_text,'(a,a)') 'Branch    : ', TRIM(branch)
      CALL message('',message_text)
      WRITE(message_text,'(a,a)') 'Revision  : ', TRIM(revision)
      CALL message('',message_text)
  
      CALL message ('','')
      message_text = '==========================================================='
      CALL message ('',TRIM(message_text))        
      CALL message ('','')
      CALL message ('','')
    END IF


    IF (l_standalone) THEN
      os_name   = 'unknown'
      user_name = 'unknown'
      host_name = 'unknown'

      CALL util_os_system (os_name,   i_lena)
      CALL util_user_name (user_name, i_lenb)
      CALL util_node_name (host_name, i_lenc)

      CALL DATE_AND_TIME(ydate, ytime)
    END IF

    DO i = 1, SIZE(label)
      label(i) = ' '
    ENDDO

    WRITE (label(1), '(a)') ' Surface model version: JSBACH '//jsbach_version
    IF (l_standalone) THEN
      WRITE (label(2), '(a)') ' Date - '//ydate(1:8)//' Time - '//ytime(1:6)
      WRITE (label(3), '(a)') ' '//user_name(1:i_lenb)//' on '//host_name(1:i_lenc)
      WRITE (label(4), '(a)') ' '//os_name(1:i_lena)
    END IF

    ! set global attributes to be written to NetCDF file
    CALL put_att (global_att,'jsbach_version',jsbach_version)
    IF (l_standalone) THEN
      CALL put_att (global_att,'date_time',ydate(1:8)//' '//ytime(1:6))
      CALL put_att (global_att,'operating_system',os_name(1:i_lena))
      CALL put_att (global_att,'user_name',user_name(1:i_lenb))
      CALL put_att (global_att,'host_name',host_name(1:i_lenc))
   END IF

  END SUBROUTINE jsbach_init_version

  !>
  !! Label a simulation
  !!
  !! Write out details of a forecast run after the setup is complete, just before
  !! computing the first timestep
  !!
  !! @par Revision History
  !! Initial version derived from ECHAM6 mo_version::label_run by Reiner Schnur, MPI Met (2011-02-10)
  !!
  SUBROUTINE jsbach_label_run(standalone)

    USE mo_mpi,               ONLY: p_parallel_io, p_io, p_bcast
    USE mo_util_string,       ONLY: separator
 
    LOGICAL, INTENT(in), OPTIONAL :: standalone !< true: offline JSBACH, false: ECHAM/JSBACH

    LOGICAL :: l_standalone !< true: offline JSBACH, false: ECHAM/JSBACH
    INTEGER :: i

    ! Executable statements

    l_standalone = .TRUE.
    IF (PRESENT(standalone)) l_standalone = standalone

    IF (l_standalone) THEN
      ! Get model name from environment
       IF (p_parallel_io) THEN
         CALL get_command_argument(0, executable)
      END IF

      CALL p_bcast (executable, p_io)

      CALL message ('','')
      CALL message('',separator)
      CALL message ('',' Model: '//TRIM(executable))
      CALL message('',separator)
    END IF

    ! Type of run
    IF (l_standalone) CALL message('',separator)
    DO i = 1, SIZE(label)
      IF (label(i) /= ' ') THEN
        WRITE (message_text,'(a)') label(i)
        CALL message ('',message_text)
      END IF
    ENDDO

    IF (l_standalone) THEN
      CALL message('',separator)
      CALL message('','')
    END IF

  END SUBROUTINE jsbach_label_run

END MODULE mo_jsbach_version
