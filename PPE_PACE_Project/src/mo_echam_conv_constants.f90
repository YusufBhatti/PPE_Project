!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!>
!! Constants used in the cumulus convection parameterization scheme
!! of the ECHAM physics package
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2010-07)
!! Second vers. by M. Esch, MPI (2015-06)
!! - updated to ECHAM-6.3 and unified for use in ICON
!!
MODULE mo_echam_conv_constants

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_conv_constants'

  PUBLIC :: cuparam

  PUBLIC :: entrpen, entrscv, entrmid, entrdd 
  PUBLIC :: centrmax, cmfctop, cminbuoy, cmaxbuoy, cbfac
  PUBLIC :: cmfcmax, cmfcmin, cmfdeps, cprcon, cmftau
  PUBLIC :: lmfmid, lmfdd, lmfdudv
#ifndef __ICON__
  PUBLIC :: nmctop, lmfpen, lmfscv
#endif

  REAL(wp) :: entrpen      !< entrainment rate for penetrative convection
  REAL(wp) :: entrscv      !< entrainment rate for shallow convection
  REAL(wp) :: entrmid      !< entrainment rate for midlevel convection
  REAL(wp) :: entrdd       !< entrainment rate for cumulus downdrafts
  REAL(wp) :: centrmax     !<
  REAL(wp) :: cmfctop      !< relat. cloud massflux at level above nonbuoyanc
  REAL(wp) :: cminbuoy     !< minimum excess buoyancy
  REAL(wp) :: cmaxbuoy     !< maximum excess buoyancy
  REAL(wp) :: cbfac        !< factor for std dev of virtual pot temp
  REAL(wp) :: cmfcmax      !< maximum massflux value allowed for
  REAL(wp) :: cmfcmin      !< minimum massflux value (for safety)
  REAL(wp) :: cmfdeps      !< fractional massflux for downdrafts at lfs
  REAL(wp) :: cprcon       !< coefficients for determining conversion
                           !<  from cloud water to rain
  REAL(wp) :: cmftau       !< characteristic adjustment time scale (s)

  LOGICAL :: lmfmid        !< true if midlevel    convection is switched on
  LOGICAL :: lmfdd         !< true if cumulus downdraft      is switched on
  LOGICAL :: lmfdudv       !< true if cumulus friction       is switched on

#ifndef __ICON__
  INTEGER :: nmctop        !< max. level for cloud base of mid level conv.
  LOGICAL :: lmfpen        !< true if penetrative convection is switched on
  LOGICAL :: lmfscv        !< true if shallow     convection is switched on
#endif

CONTAINS
  !>
  !! Description:
  !!
  !! Defines disposable parameters for massflux scheme
  !!
  !! Method:
  !!
  !! This routine is called from *iniphy*
  !!
  !! Authors:
  !!
  !! M. Tiedtke, ECMWF, February 1989, original source
  !! L. Kornblueh, MPI, May 1998, f90 rewrite
  !! U. Schulzweida, MPI, May 1998, f90 rewrite
  !! A. Rhodin, MPI, Jan 1999, subroutine cuparam -> module mo_cumulus_flux
  !! M. Esch, MPI, July 1999, modifications for ECHAM5
  !! U. Schlese, MPI, August 2000, mid level cloud base *nmctop*
  !! 
  !! for more details see file AUTHORS
  !! 
  SUBROUTINE cuparam

   USE mo_control,      ONLY: nlev, nlevp1, nvclev, vct, nn, lmidatm
   USE mo_exception,    ONLY: finish, message, message_text, em_param
   !>>SF
   USE mo_param_switches, ONLY: ncd_activ, & !SF ncd_activ replaces former ncdnc
                                lcdnc_progn !SF
   !<<SF

  IMPLICIT NONE

! local variables
  REAL(wp)    :: za, zb, zph(nlevp1), zp(nlev)
  INTEGER :: jk

  !  Executable Statements 

!-- 1. Specify parameters for massflux-scheme

  entrpen  = 1.0E-4_wp !
  entrscv  = 3.0E-3_wp 
  cminbuoy = 0.2_wp
  cmaxbuoy = 1.0_wp
  cbfac    = 1.0_wp
  entrmid  = 1.0E-4_wp ! Average entrainment rate for midlevel convection
  entrdd   = 2.0E-4_wp ! Average entrainment rate for downdrafts
  centrmax = 3.E-4_wp !
  cmfcmax  = 1.0_wp ! Maximum massflux value allowed for updrafts etc
  cmfcmin  = 1.E-10_wp ! Minimum massflux value (for safety)
  cmfdeps  = 0.3_wp ! Fractional massflux for downdrafts at lfs

#ifdef __ICON__
  cmfctop = 0.2_wp
  cprcon  = 2.5E-4_wp
  cmftau  = 7200._wp
#else

!
  IF (nn == 31) THEN
    cmfctop = 0.2_wp
    cprcon  = 2.5E-4_wp
  ELSE IF (nn == 63) THEN
    cmfctop = 0.2_wp
    cprcon  = 2.5E-4_wp
  ELSE IF (nn == 127) THEN
    cmfctop = 0.23_wp
    cprcon  = 1.5E-4_wp
  ELSE IF (nn == 255) THEN
    cmfctop = 0.23_wp
    cprcon  = 1.5E-4_wp
  ELSE
    CALL finish ('mo_echam_conv_constants', 'Truncation not supported.')
  ENDIF
  cmftau  = MIN(3._wp*3600._wp,7200._wp*63._wp/nn)

!>>SF Modifications for double-moment cloud microphysics scheme
  IF (lcdnc_progn) THEN
     IF (nn == 63) THEN
        IF (nlev == 31) THEN
           SELECT CASE (ncd_activ)
              CASE(1) ! LL activation
                 !SF: updated on 2015.02.25 (David Neubauer / Katty Huang, pure atm run, HAM-M7, LL activation)
                 entrscv = 3.0E-3_wp
                 entrpen = 1.E-4_wp
                 cmfctop = 0.2_wp
                 cprcon  = 9.E-04_wp
              CASE(2) ! AR&G activation
                 !SF: updated on 2017.01.27 (David Neubauer, pure atm run, HAM-M7) 
                 entrscv = 3.0E-3_wp
                 entrpen = 2.E-4_wp
                 cmfctop = 0.2_wp
                 cprcon  = 9.E-04_wp
           END SELECT
        ELSE IF (nlev == 47) THEN
           SELECT CASE (ncd_activ)
              CASE(1) ! LL activation
                 !SF: updated on 2015.02.19 (David Neubauer, pure atm run, HAM-M7, LL activation)
                 entrscv = 3.0E-3_wp
                 entrpen = 1.0E-4_wp
                 cmfctop = 0.2_wp
                 cprcon  = 2.5E-04_wp
              CASE(2) ! AR&G activation
                 !SF: updated on 2017.01.27 (David Neubauer, pure atm run, HAM-M7) 
                 entrscv = 3.0E-3_wp
                 entrpen = 2.E-4_wp
                 cmfctop = 0.2_wp
                 cprcon  = 9.E-04_wp
           END SELECT
        ENDIF 
     ENDIF
  ENDIF
!<<SF

!>>SF
  WRITE(message_text,'(a,e25.15)') 'Entrainment rate for shallow convection: entrscv= ',entrscv
  CALL message('mo_cumulus_flux', message_text, level=em_param)
  WRITE(message_text,'(a,e25.15)') 'Entrainment rate for penetrative convection: entrpen= ',entrpen
  CALL message('mo_cumulus_flux', message_text, level=em_param)
  WRITE(message_text,'(a,e25.15)') 'relat. cloud massflux at level above nonbuoyanc: cmfctop=',cmfctop
  CALL message('mo_cumulus_flux', message_text, level=em_param)
  WRITE(message_text,'(a,e25.15)') 'autoconversion for convect. clouds: cprcon=',cprcon
  CALL message('mo_cumulus_flux', message_text, level=em_param)
!<<SF

! Determine highest level *nmctop* for cloud base of midlevel convection
! assuming nmctop=9 (300 hPa) for the standard 19 level model

!-- half level pressure values, assuming 101320. Pa surface pressure 

  DO jk=1,nlevp1
    za=vct(jk)
    zb=vct(jk+nvclev)
    zph(jk)=za+zb*101320.0_wp
  END DO
!
! -- full level pressure
!
  DO jk = 1, nlev
    zp(jk)=(zph(jk)+zph(jk+1))*0.5_wp
  END DO
!
! -- search for 300 hPa level
!
  DO jk = 1, nlev
    nmctop=jk
    IF(zp(jk).GE.30000.0_wp) EXIT
  END DO
!
  WRITE (message_text,'(a,i0)') &
       'max. level for cloud base of mid level convection: nmctop = ', nmctop
  CALL message('',message_text)
#endif
  
END SUBROUTINE cuparam

END MODULE mo_echam_conv_constants
