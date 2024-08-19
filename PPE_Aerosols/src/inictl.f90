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
SUBROUTINE inictl

  ! Description:
  !
  ! Preset constants in mo_control.
  !
  ! Method:
  !
  ! Calculate space for memory manager
  !
  ! *inictl* is called from *initialize*.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, August 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! H.-S. Bauer, MPI, Jul 1998, changed
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! M. Esch, MPI, July 1999, remove nudging, nmi switches
  ! M. Esch, MPI, July 1999, modifications for ECHAM5
  ! U. Schlese, DKRZ, December 1999, modifications for coupling
  ! I. Kirchner, MPI, October 2000, revision, time control
  ! L. Kornblueh, MPI, January 2001, revision of time control
  ! I. Kirchner, MPI, March 2001, revision
  ! A. Rhodin, MPI, June 2001, g3x,g4x fields removed
  ! L. Kornblueh, MPI, October 2001, added missing broadcast of nsub in runctl
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! M. Esch, MPI, September 2002, switch for mixed layer ocean
  ! L. Kornblueh, MPI, April 2003, switch for port test added
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! S. Lorenz, MPI, November 2007, added switch for volcanic forcing
  ! U. Schlese, MPI, November 2008, setting of "iaero" removed, "lso4" removed
  ! D. Popke, MPI, Nov 2013, Implementation of RCE (lrce)
  !
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_mpi,           ONLY: p_io, p_parallel, p_parallel_io, p_bcast
  USE mo_control,       ONLY: lamip, lcouple, ldebug, lmidatm, lmlo,          &
                              lnmi, lnudge, lnwp, ltdiag, lhd,                &
                              lcolumn, ldiagamip, lprint_m0, ldebugio,        &
                              ldebugmem, ltimer, ltctest, nproma,             &
                              ldailysst, lmeltpond, lcouple_co2,              &
                              ldebugs, l_volc, lindependent_read,             &
                              lcollective_write, lfractional_mask,            &
                              lforcererun,                                    & !SF
                              earth_angular_velocity, lrce,                   &
                              rmlo_depth, lmlo_ice
  USE mo_port_test,     ONLY: lport
  USE mo_namelist,      ONLY: open_nml, position_nml, close_nml,              &
                              POSITIONED, MISSING     
  USE mo_time_base,     ONLY: set_calendar_type, CYL360
  USE mo_time_control,  ONLY: delta_time, no_days, no_steps, no_cycles,       &
                              dt_start, dt_stop, dt_resume,                   &
                              l_orbvsop87,                                    &
                              putdata, putrerun, putcheckpoint,               &
                              putocean, getocean,                             &
                              trigjob, nsub, subflag,                         &
                              p_bcast_event, NSUB_MAX, trigfiles,             &
                              lresume, lresume_cmd, ldebugev, lfirst_cycle
  USE mo_advection,     ONLY: iadvec
  USE mo_hdiff,         ONLY: ndiahdf
  USE mo_filename,      ONLY: out_datapath, out_expname,                      &
                              out_ztype,                                      &
                              out_filetype, trac_filetype, rerun_filetype
  USE mo_memory_base,   ONLY: default_output

#ifdef __prism
  USE mo_couple,        ONLY: lcouple_parallel, ldebugcpl
#endif

  IMPLICIT NONE

  ! Local scalars: 

  LOGICAL :: ly360   = .FALSE.

  INTEGER ::  i, inml, iunit

  INTEGER :: ierr  ! error return value from position_nml

#include "runctl.inc"

  ! Executable statements 

  ! 1. Preset constants

  ! 2. Read namelist runctl

  IF (p_parallel_io) THEN
     inml = open_nml ('namelist.echam')
     iunit = position_nml ('RUNCTL', inml, status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       trac_filetype = 0
       READ (iunit, runctl)
       IF(trac_filetype == 0) trac_filetype = out_filetype
     CASE(MISSING)
       CALL message('inictl', 'Cannot find namelist RUNCTL in namelist.echam!')
     END SELECT
     CALL close_nml(inml)
  ENDIF
  IF (p_parallel) THEN
     CALL p_bcast (ltimer, p_io)
     CALL p_bcast (ldebugio, p_io)
     CALL p_bcast (ldebugmem, p_io)
     CALL p_bcast (ldebugev, p_io)
     CALL p_bcast (ltctest, p_io)
     CALL p_bcast (ldebugs, p_io)
     CALL p_bcast (lindependent_read, p_io)
     CALL p_bcast (lcollective_write, p_io)

     CALL p_bcast (lresume, p_io)

     CALL p_bcast (out_datapath, p_io)
     CALL p_bcast (out_expname,  p_io)
     CALL p_bcast (rerun_filetype, p_io)
     CALL p_bcast (out_filetype, p_io)
     CALL p_bcast (trac_filetype, p_io)

     CALL p_bcast (out_ztype, p_io)
     CALL p_bcast (default_output, p_io)

     CALL p_bcast (lprint_m0, p_io)
     CALL p_bcast (ldebug, p_io)
     CALL p_bcast (ldailysst, p_io)
     CALL p_bcast (lamip, p_io)
     CALL p_bcast (lcolumn, p_io)
     CALL p_bcast (ldiagamip, p_io)
     CALL p_bcast (lnwp, p_io)
     CALL p_bcast (lnudge, p_io)
     CALL p_bcast (lmidatm, p_io)
     CALL p_bcast (lmlo, p_io)
     CALL p_bcast (lmeltpond, p_io)
     CALL p_bcast (iadvec, p_io)
     CALL p_bcast (lnmi, p_io)
     CALL p_bcast (ltdiag, p_io)
     CALL p_bcast (lport, p_io)
     CALL p_bcast (nproma, p_io)
     CALL p_bcast (lcouple, p_io)
     CALL p_bcast (lcouple_co2, p_io)
     CALL p_bcast (lfractional_mask, p_io)

     CALL p_bcast (l_orbvsop87, p_io)
     CALL p_bcast (ly360, p_io)
     CALL p_bcast (delta_time, p_io)
        
     CALL p_bcast (dt_start, p_io)
     CALL p_bcast (dt_resume, p_io)
     CALL p_bcast (dt_stop, p_io)

     CALL p_bcast (no_days, p_io)
     CALL p_bcast (no_cycles, p_io)
     CALL p_bcast (no_steps, p_io)

     CALL p_bcast (nsub, p_io)

     CALL p_bcast_event(putdata, p_io)
     CALL p_bcast_event(trigfiles, p_io)
     CALL p_bcast_event(putrerun, p_io)
     CALL p_bcast_event(putcheckpoint, p_io)

     nsub = MIN(MAX(nsub,0),NSUB_MAX)
     CALL p_bcast(subflag, p_io)
     DO i=1,nsub
        CALL p_bcast_event(trigjob(i), p_io)
     END DO
     CALL p_bcast_event(putocean, p_io)
     CALL p_bcast_event(getocean, p_io)

     CALL p_bcast (lhd, p_io)

     CALL p_bcast (l_volc, p_io)

     CALL p_bcast (ndiahdf, p_io)

     CALL p_bcast(earth_angular_velocity, p_io)
     CALL p_bcast(lrce, p_io)

     CALL p_bcast (lmlo_ice, p_io)
     CALL p_bcast (rmlo_depth, p_io)

!>>SF: switchable internal reruns
     CALL p_bcast (lforcererun, p_io)
!<<SF 

#ifdef __prism
     CALL p_bcast (lcouple_parallel, p_io)
     CALL p_bcast (ldebugcpl, p_io)
#endif
  ENDIF

  ! if lresume is set via command line, ignore namelist entry

  IF (lresume_cmd) THEN
    CALL message('inictl', &
         ' lresume taken from CLI - ignore namelist entry.')
    lresume = .FALSE.
  ENDIF

  ! reset lresume for the second rerun cycle during an initial run
  IF (.NOT. lfirst_cycle) lresume = .TRUE.

  ! check consistency of orbit and year length

  IF (p_parallel_io) THEN
    IF (l_orbvsop87 .AND. ly360) THEN
      CALL finish('inictl', &
           ' ly360=.TRUE. cannot run with real orbit (l_orbvsop87=.TRUE.).')
    ENDIF
  ENDIF
  IF (ly360) THEN
    CALL set_calendar_type(CYL360)
  ENDIF

  IF ( nproma < 0 ) nproma = 0  ! default is lon/lat ordering

!lk  ! correct number of rerun cycles
!lk  no_cycles = MAX(no_cycles,1)
  IF (no_cycles < 1) THEN
    no_cycles = 0
    CALL message('inictl', &
         ' Debug mode - DT_STOP is finishing this specific run.')
  END IF

  ! correct maximal number of subjobs
  nsub = MIN(MAX(nsub,0),NSUB_MAX)

  ! if NWP mode we always have a restart
  IF (lnwp) lresume = .TRUE.

  IF (lamip .AND. ldailysst) THEN
    CALL finish('inictl', &
           ' ldailysst=.TRUE. does not work with lamip=.TRUE.')
  END IF

  IF (lrce) THEN
    WRITE (message_text, '(a,f12.8)') &
         'earth_angular_velocity --> Solid Earth rotation velocity =', earth_angular_velocity
    CALL message('inictl',message_text)

    WRITE (message_text, '(a,L5)') &
         'lrce --> switch for radiative-convective eq. radiation =', lrce
    CALL message('inictl',message_text)
  END IF

END SUBROUTINE inictl
