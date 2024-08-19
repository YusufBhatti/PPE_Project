!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_diags.f90
!!
!! \brief
!! Main module for AeroCom diagnostics
!!
!! \author D. Neubauer (ETH Zurich)
!!  adapted from AerChemMIP diagnostics from S. Ferrachat
!!
!! \responsible_coder
!! D. Neubauer, david.neubauer@env.ethz.ch
!!
!! \revision_history
!!   -# D. Neubauer (ETH Zurich) - original code (2018-08-14)
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_hammoz_aerocom_diags

  USE mo_kind, ONLY: dp
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_aerocom
  PUBLIC :: init_aerocom_streams
  PUBLIC :: update_aerocom_diags

! Variable declarations
  LOGICAL, PUBLIC :: lTraj = .FALSE.    !< switch for 'Traj' diags
  LOGICAL, PUBLIC :: lHEmon = .FALSE.   !< switch for 'HEmon' diags
  LOGICAL, PUBLIC :: lHEaci = .FALSE.   !< switch for 'HEaci' diags
  LOGICAL, PUBLIC :: lHEaci_activ = .TRUE.   !< switch for 'HEaci' activation diags
  LOGICAL, PUBLIC :: lHEpro = .FALSE.   !< switch for 'HEpro' diags
  LOGICAL, PUBLIC :: lMMPPE = .FALSE.   !< switch for 'MMPPE' diags
  LOGICAL, PUBLIC :: lAP3M = .FALSE.   !< switch for 'AP3 monthly' diags
  LOGICAL, PUBLIC :: lAP3D = .FALSE.   !< switch for 'AP3 daily' diags

  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! init_aerocom: initialization routine for AeroCom diags
!!
!! @author see module info
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Notes
!!
!! @par Responsible coder
!! david.neubauer@env.ethz.ch
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE init_aerocom

    USE mo_mpi,       ONLY: p_io, p_parallel_io, p_parallel, p_bcast
    USE mo_namelist,  ONLY: open_nml, position_nml, POSITIONED, &
                            LENGTH_ERROR, READ_ERROR
    USE mo_submodel,  ONLY: print_value, laoa
    USE mo_exception, ONLY: finish, message, em_error, em_warn
    USE mo_control,   ONLY: ltdiag
    USE mo_hammoz_aerocom_HEaci, ONLY: init_HEaci
    USE mo_hammoz_aerocom_AP3M,  ONLY: init_AP3M
    
    INTEGER :: ierr, inml, iunit

    INCLUDE 'hammoz_aerocomctl.inc'

    !-- Handle aerocom namelist
    IF (p_parallel_io) THEN
       inml = open_nml('namelist.echam')
       iunit = position_nml ('HAMMOZ_AEROCOMCTL', inml, status=ierr)
       SELECT CASE (ierr)
       CASE (POSITIONED)
          READ (iunit, hammoz_aerocomctl)
       CASE (LENGTH_ERROR)
          CALL finish ('init_aerocom', &
                       'length error in namelist hammoz_aerocomctl')
       CASE (READ_ERROR)
          CALL finish ('init_aerocom', &
                       'read error in namelist.echam')
       END SELECT
    END IF

    IF (p_parallel) THEN
        CALL p_bcast (lTraj, p_io)
        CALL p_bcast (lHEmon, p_io)
        CALL p_bcast (lHEaci, p_io)
        CALL p_bcast (lHEaci_activ, p_io)
        CALL p_bcast (lHEpro, p_io)
        CALL p_bcast (lMMPPE, p_io)
        CALL p_bcast (lAP3M, p_io)
        CALL p_bcast (lAP3D, p_io)
    END IF

    !-- Security for dependencies

    !-- Print status
    IF (p_parallel_io) THEN
        CALL print_value('init_aerocom, lTraj', lTraj)
        CALL print_value('init_aerocom, lHEmon', lHEmon)
        CALL print_value('init_aerocom, lHEaci', lHEaci)
        CALL print_value('init_aerocom, lHEaci_activ', lHEaci)
        CALL print_value('init_aerocom, lHEpro', lHEpro)
        CALL print_value('init_aerocom, lMMPPE', lMMPPE)
        CALL print_value('init_aerocom, lAP3M', lAP3M)
        CALL print_value('init_aerocom, lAP3D', lAP3D)
    ENDIF
    IF (lHEaci) THEN
       CALL init_HEaci
    END IF
    IF (lAP3M) THEN
       CALL init_AP3M
    END IF

END SUBROUTINE init_aerocom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! init_aerocom_streams: initialize AeroCom diags streams
!!
!! @author see module info
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! init_subm
!!
!! @par Notes
!!
!! @par Responsible coder
!! david.neubauer@env.ethz.ch
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE init_aerocom_streams

    USE mo_hammoz_aerocom_Traj,  ONLY: construct_Trajg_stream
    USE mo_hammoz_aerocom_Traj,  ONLY: construct_Traj_stream
    USE mo_hammoz_aerocom_HEmon, ONLY: construct_HEmon_stream
    USE mo_hammoz_aerocom_HEaci, ONLY: construct_HEaci_stream
    USE mo_hammoz_aerocom_HEpro, ONLY: construct_HEpro_stream
    USE mo_hammoz_aerocom_MMPPE, ONLY: construct_MMPPE_stream
    USE mo_hammoz_aerocom_AP3M,  ONLY: construct_AP3M_stream
    USE mo_hammoz_aerocom_AP3D,  ONLY: construct_AP3D_stream
                                     

    IF (lTraj)         CALL construct_Trajg_stream
    IF (lTraj)         CALL construct_Traj_stream
    IF (lHEmon)        CALL construct_HEmon_stream
    IF (lHEaci)        CALL construct_HEaci_stream

    IF (lHEpro)        CALL construct_HEpro_stream
    IF (lMMPPE)        CALL construct_MMPPE_stream
    IF (lAP3M)         CALL construct_AP3M_stream
    IF (lAP3D)         CALL construct_AP3D_stream

END SUBROUTINE init_aerocom_streams
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! update_aerocom_diags: driver for all AeroCom diags
!!
!! @author see module info
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! physc
!!
!! @par Notes
!!
!! @par Responsible coder
!! david.neubauer@env.ethz.ch
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE update_aerocom_diags(kproma, kbdim, klev, krow, pi0)

    USE mo_hammoz_aerocom_Traj,         ONLY: update_Traj_diags
    USE mo_hammoz_aerocom_HEmon,        ONLY: update_HEmon_diags
    USE mo_hammoz_aerocom_HEaci,        ONLY: update_HEaci_diags
    USE mo_hammoz_aerocom_HEpro,        ONLY: update_HEpro_diags
    USE mo_hammoz_aerocom_MMPPE,        ONLY: update_MMPPE_diags
    USE mo_hammoz_aerocom_AP3M,         ONLY: update_AP3M_diags
    USE mo_hammoz_aerocom_AP3D,         ONLY: update_AP3D_diags
      
    USE mo_tracdef, ONLY: ntrac

    INTEGER, INTENT(in)  :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(IN )  :: pi0(kbdim)

    IF (lTraj)  CALL update_Traj_diags(kproma, kbdim, klev, krow)
    IF (lHEaci) CALL update_HEaci_diags(kproma, kbdim, klev, krow, pi0)
    IF (lHEmon) CALL update_HEmon_diags(kproma, kbdim, klev, krow)
    IF (lHEpro) CALL update_HEpro_diags(kproma, kbdim, klev, krow)
    IF (lMMPPE) CALL update_MMPPE_diags(kproma, kbdim, klev, krow)
    IF (lAP3M)  CALL update_AP3M_diags(kproma, kbdim, klev, krow)
    IF (lAP3D)  CALL update_AP3D_diags(kproma, kbdim, klev, krow)

    lHEaci_activ = .TRUE.

END SUBROUTINE update_aerocom_diags

END MODULE mo_hammoz_aerocom_diags
