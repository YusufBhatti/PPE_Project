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
MODULE mo_version

  IMPLICIT NONE

  !---------------------------------------------------------------------------
  !
  ! Authors:
  !
  ! J. K. Gibson, ECMWF, May 1980, original source
  ! L. Kornblueh, MPI, November 2008, split for version and labeling and 
  !                                   io-units
  ! L. Kornblueh, MPI, Februrary 2012, remove spitfire
  !
  !---------------------------------------------------------------------------
  !
  ! cycle number of original ECMWF model (pre IFS)
  
  INTEGER, PARAMETER :: ncycle = 36 

  ! model version number (6.0.00 => 6.0 * 10) 

  INTEGER :: nversion

  ! basic information

  CHARACTER (len=256) :: executable
  CHARACTER (len=256) :: version_tag
  CHARACTER (len=256) :: os_name
  CHARACTER (len=256) :: user_name
  CHARACTER (len=256) :: host_name


  ! labels contain space for storing model information

  CHARACTER (len= 80) :: label(8)

CONTAINS

  SUBROUTINE init_version

    USE mo_netcdf,         ONLY: global_att, put_att
    USE mo_jsbach_version, ONLY: jsbach_init_version
    
    !  Local scalars: 
    REAL :: zvers
    
    CHARACTER (len=  3)   :: yovers
    CHARACTER (len=  6)   :: version_tag

    CHARACTER (len=  8)   :: ydate
    CHARACTER (len= 10)   :: ytime

    CHARACTER (len=256)   :: ytmp
    
    INTEGER :: nlena, nlenb, nlenc
    
    INTEGER :: i
    
    !  External subroutines
    EXTERNAL :: util_os_system, util_user_name, util_node_name
    
    !  Executable statements 
    
    version_tag = '6.3.02'
    yovers    = '6.3'
    READ (yovers,'(f8.1)') zvers
    nversion = NINT(10*zvers)
    
    os_name   = ''
    user_name = ''
    host_name = ''
    
    ytmp = ''
    CALL util_os_system (ytmp, nlena)
    os_name = ytmp(1:nlena)

    ytmp = ''
    CALL util_user_name (ytmp, nlenb)
    user_name = ytmp(1:nlenb)

    ytmp = ''
    CALL util_node_name (ytmp, nlenc)
    host_name = ytmp(1:nlenc)    

    DO i = 1, SIZE(label)
      label(i) = ' '
    ENDDO
    
    CALL DATE_AND_TIME(ydate, ytime)
    
    WRITE (label(1), '(a80)') ' Atmospheric model version '//version_tag
    WRITE (label(2), '(a80)') ' Library 18-Feb-2016'
    WRITE (label(3), '(a80)') ' Lin & Rood ADVECTION is default'
    WRITE (label(4), '(a80)') ' Modified ECMWF physics'
    WRITE (label(5), '(a80)') ' Using PSrad/RRTMG radiation'
    WRITE (label(6), '(a80)') ' Date - '//ydate(1:8)//' Time - '//ytime(1:6)
    WRITE (label(7), '(a80)') ' '//TRIM(user_name)//' on '//TRIM(host_name) 
    WRITE (label(8), '(a80)') ' '//TRIM(os_name)
    
    ! set global attributes to be written to NetCDF file
    
    CALL put_att (global_att,'echam_version',version_tag)
    CALL put_att (global_att,'advection','Lin & Rood')
    CALL put_att (global_att,'physics','Modified ECMWF physics')
    CALL put_att (global_att,'radiation','Using PSrad/RRTMG radiation')
    CALL put_att (global_att,'date_time',ydate(1:8)//' '//ytime(1:6))
    CALL put_att (global_att,'operating_system',TRIM(os_name))
    CALL put_att (global_att,'user_name',TRIM(user_name))
    CALL put_att (global_att,'host_name',TRIM(host_name))
    
    CALL jsbach_init_version(standalone=.FALSE.)

  END SUBROUTINE init_version

  SUBROUTINE label_run

    ! Description:
    !
    ! Label a forecast run.
    !
    ! Method:
    !
    ! Write out details of a forecast run after the set-up is
    ! complete, just before computing the first timestep.
    !
    ! Various items are printed from modules, and the forecast
    ! namelists are written.
    !
    ! Authors:
    !
    ! J. K. Gibson, ECMWF, February 1983, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! U. Schlese, DKRZ, Aug 1999, modifications for ECHAM5 (SPITFIRE)
    ! L. Kornblueh, MPI, November 2008, modifications for adding HAM and JSBACH
    ! 

    USE mo_control,           ONLY: nlev, ngl, nlon, lmlo, lmidatm, vct,  &
                                    lamip, nvclev, lhd, lcouple,          &
                                    ldailysst, lmeltpond
    USE mo_semi_impl,         ONLY: betadt, betazq, apr, tr, eps
    USE mo_param_switches,    ONLY: lphys, lgwdrag, lsurf, lcond,         &
                                    lvdiff, lconv, lice, lrad
    USE mo_tracdef,           ONLY: ntrac
    USE mo_time_control,      ONLY: delta_time, lstart, lresume
    USE mo_mpi,               ONLY: p_parallel_io, p_io, p_bcast
    USE mo_exception,         ONLY: message, message_text
    USE mo_util_string,       ONLY: separator
    USE mo_advection,         ONLY: iadvec, semi_lagrangian, tpcore
    USE mo_jsbach_version,    ONLY: jsbach_label_run

    INTEGER :: i, istart, ilast

    !  Executable statements

    ! get model name from environment

    IF (p_parallel_io) THEN

      CALL get_command_argument(0, executable)

    END  IF

    CALL p_bcast (executable, p_io)

    CALL message ('','')
    CALL message('',separator)
    CALL message ('',' Model: '//TRIM(executable)) 
    CALL message('',separator)

    ! type of run

    CALL message('',separator)
    DO i = 1, SIZE(label)
      WRITE (message_text,'(a)') label(i)
      CALL message ('',message_text)
    ENDDO

    ! Print copyright for advection scheme

    SELECT CASE (iadvec)
    CASE (semi_lagrangian)
      CALL message('',separator)
      CALL message ('',' The semi Lagrangian transport scheme is based on the')
      CALL message ('',' NCAR Community Climate Model (CCM2)')
      CALL message ('',' Version 2.1.2 [02/07/94]/, Copyright (C) 1993')
      CALL message ('',' University Corporation for Atmospheric Research')
      CALL message ('',' All Rights Reserved.')
      CALL message('',separator)
    CASE (tpcore)
      CALL message('',separator)
      CALL message ('',' TransPort of NASA Goddard Chemistry Transport Model')
      CALL message ('',' using a Flux Form Semi-Lagrangian (FFSL) scheme')
      CALL message ('',' by Shian-Jiann Lin et al., NASA - GSFC, 2001')
      CALL message('',separator)
    END SELECT

    IF(lhd) THEN
      CALL message('',separator)
      CALL message ('',' Running Version 2.0 of Hydrological Discharge Model   ')
      CALL message('',separator)
    END IF
    
    IF (lresume) THEN
       CALL message ('',' Restarted run (from history files)')
    ELSE IF (lstart) THEN
       CALL message ('',' Initial run')
    END IF
    
    CALL message('',separator)
    
    CALL message ('',' General runtime parameter: ')
    WRITE (message_text,'(a,i7)')   '   number of vertical levels.                         (nlev) = ', nlev
    CALL message ('',message_text)
    WRITE (message_text,'(a,i7)')   '   number of gaussian latitudes.                       (ngl) = ', ngl
    CALL message ('',message_text)
    WRITE (message_text,'(a,i7)')   '   max number of points on each latitude line         (nlon) = ', nlon
    CALL message ('',message_text)
    CALL message ('','')
    WRITE (message_text,'(a,f6.1)') '   integration time stepping                  (2*delta_time) = ', 2*delta_time
    CALL message ('',message_text)
    WRITE (message_text,'(a,f7.3)') '   time filtering coefficient                          (eps) = ', eps             
    CALL message ('',message_text)
    WRITE (message_text,'(a,f4.1)') '   explicit scheme for d, t, alps (= 0.0)           (betadt) = ', betadt
    CALL message ('',message_text)
    CALL message ('','   semi implicit scheme (= 1.0)                                ')
    WRITE (message_text,'(a,f4.1)') '   explicit scheme for vo, q (= 0.0)                (betazq) = ', betazq
    CALL message ('',message_text)
    CALL message ('','   semi implicit scheme (= 1.0)                                ')
    WRITE (message_text,'(a,e9.3)') '   reference surface pressure for semi-implicit scheme (apr) = ', apr
    CALL message ('',message_text)
    WRITE (message_text,'(a,e9.3)') '   reference temperature for semi-implicit scheme       (tr) = ', tr
    CALL message ('',message_text)
    
    CALL message ('',' Physics switches: ')
    WRITE (message_text,'(a,l2)') '   physics                  (lphys)    = ', lphys
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   radiation                (lrad)     = ', lrad
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   gravity wave drag        (lgwdrag)  = ', lgwdrag
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   surface exchanges        (lsurf)    = ', lsurf
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   large scale condensation (lcond)    = ', lcond
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   vertical diffusion       (lvdiff)   = ', lvdiff
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   convection               (lconv)    = ', lconv
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   surface ice              (lice)     = ', lice
    
    CALL message ('',' Run control switches: ')
    WRITE (message_text,'(a,l2)') '   middle atmosphere        (lmidatm)  = ', lmidatm
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   mixed layer ocean        (lmlo)     = ', lmlo
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   meltpond calculation     (lmeltpond)= ', lmeltpond
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   full ocean coupling      (lcouple)  = ', lcouple
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   AMIP run                 (lamip)    = ', lamip
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   daily SST and SIC        (ldailysst)= ', ldailysst
    CALL message ('',message_text)
    WRITE (message_text,'(a,l2)') '   hydr. discharge model    (lhd)      = ', lhd
    CALL message ('',message_text)
    
    CALL message ('','')

    CALL message ('',' Vertical coordinate table (VCT)')
    CALL message ('',' Parameter A:')
    DO i = 1, nvclev, 10
       istart = i
       ilast  = i+9
       IF (ilast > nvclev) ilast = nvclev
       WRITE (message_text, '(10f7.0)') vct(istart:ilast)
       CALL message ('',message_text)
    ENDDO
    CALL message ('',' Parameter B:')
    DO i = nvclev+1, 2*nvclev, 10
       istart = i
       ilast  = i+9
       IF (ilast > 2*nvclev) ilast = 2*nvclev
       WRITE (message_text, '(10f7.4)') vct(istart:ilast)
       CALL message ('',message_text)
    ENDDO
    CALL message('',separator)
    
    ! Print some tracer information
    
    CALL message ('','ECHAM6 - transport of specific humidity, ')
    WRITE (message_text,'(a,i0,a)') 'cloud water and ice, and ', ntrac, ' trace gase(s).'
    CALL message ('',message_text)
    CALL message ('',separator)

    CALL jsbach_label_run(standalone=.FALSE.)
    
  END SUBROUTINE label_run
  
END MODULE mo_version
