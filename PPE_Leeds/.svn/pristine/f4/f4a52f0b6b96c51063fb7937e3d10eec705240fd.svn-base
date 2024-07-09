!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide scaling factors for daily varying photolysis rates.
!!
!! @remarks
!!   This module contains routines two routines to read (read_dailyetf) and set (set_dailyetf)
!!   the photolysis rate scaling factors (i.e. daily varying extra terrestrial solar flux
!!   at a variety of spaectral bands). read_dailyetf is called once a year to read a full
!!   etfphot_YYYY.nc file into module variables. At the start of a new day set_dailyetf reads
!!   the previous, current and next day of etfphot, etfphot_ms93 and rsf_sclfac into module
!!   variables. The interpolation onto the current time of day is done in mo_moz_j* using
!!   ECHAM interpolation routines. Routines to generate the etfphot_YYYY.nc input files are
!!   available from swahl@geomar.de
!!
!! @author J. Kieser, MPI, July 2008
!!
!!         H. Schmidt, MPI, July 2008
!!
!!         S. Wahl, GEOMAR, June 2015
!!              code copied from HAMMONIA (ECHAM5.4) and adjusted for ECHAM6.3-HAMMOZ
!!
!! $ID: n/a$
!!
!
MODULE mo_moz_dailyetf

USE mo_kind,               ONLY: dp
USE mo_exception,          ONLY: finish, message, message_text


IMPLICIT NONE

PRIVATE

INTEGER, PARAMETER :: dimlength=368
INTEGER, PARAMETER :: numwl_etf=33 ! Hauke had 34
INTEGER, PARAMETER :: numwl_etf_ms=4
INTEGER, PARAMETER :: numwl_rsf=67 ! Hauke had 122

PUBLIC:: read_dailyetf, set_dailyetf

REAL(dp),PUBLIC,POINTER,save:: etfphot_year(:,:)
REAL(dp),PUBLIC,POINTER,save:: etfphot_ms93_year(:,:)
REAL(dp),PUBLIC,POINTER,save:: rsf_sclfac_year(:,:)
INTEGER ,PUBLIC,POINTER,save:: etftime_year(:)

REAL(dp),DIMENSION(numwl_etf,3),PUBLIC        :: etfphot_daily      !requested by jshort
REAL(dp),DIMENSION(numwl_etf_ms,3),PUBLIC     :: etfphot_ms93_daily   !requested by jshort
REAL(dp),DIMENSION(numwl_rsf,3),PUBLIC        :: rsf_sclfac_daily      !requested by jlong

CONTAINS

SUBROUTINE read_dailyetf(iyear)  ! called by src/stepon and  src/control

  USE mo_kind,         ONLY: dp
  USE mo_mpi,          ONLY: p_parallel_io,&
                           p_parallel,   &
                           p_io,         &
                           p_bcast,      &
                           p_pe
  USE mo_io,           ONLY: io_inq_dimid, io_inq_dimlen

  USE mo_read_netcdf77,ONLY: read_var_nf77_1d, read_var_nf77_2d

  !----variables-declaration------------
  INTEGER, INTENT(in) :: iyear

  ! netcdf input variables
  CHARACTER(len=256) :: time_nc
  CHARACTER(len=256) :: filename_nc
  CHARACTER(len=4)  :: cyr
  CHARACTER(len=256) :: varname_v1_nc, varname_v2_nc, varname_v3_nc, varname_time_nc
  CHARACTER(len=256) :: dim_wl_etfphot, dim_wl_ms93, dim_wl_rsf

  REAL(dp), ALLOCATABLE :: varptr_time_nc(:)

  REAL(dp), ALLOCATABLE :: varptr_1d_v1_nc(:)
  REAL(dp), ALLOCATABLE :: varptr_1d_v2_nc(:)

  REAL(dp), ALLOCATABLE :: varptr_2d_v1_nc(:,:)
  REAL(dp), ALLOCATABLE :: varptr_2d_v2_nc(:,:)
  REAL(dp), ALLOCATABLE :: varptr_2d_v3_nc(:,:)

  INTEGER    :: ierr

  ! read data

  time_nc="time"
  dim_wl_etfphot="wl_etfphot"
  dim_wl_ms93="wl_etfphot_ms93"
  dim_wl_rsf="wl_rsf_sclfac"

  varname_time_nc="date"
  varname_v1_nc="etfphot"
  varname_v2_nc="etfphot_ms93"
  varname_v3_nc="rsf_sclfac"

  WRITE(cyr,'(i4.4)') iyear
  filename_nc='etfphot_'//TRIM(ADJUSTL(cyr))//'.nc'

  CALL message('read_etfdaily','Read etfphot* and rsf_sclfac from file '//filename_nc)

  ALLOCATE (varptr_time_nc(dimlength))
  ALLOCATE (varptr_2d_v1_nc(numwl_etf,dimlength))
  ALLOCATE (varptr_2d_v2_nc(numwl_etf_ms,dimlength))
  ALLOCATE (varptr_2d_v3_nc(numwl_rsf,dimlength))

  IF (p_parallel_io) THEN
    CALL read_var_nf77_1d(filename_nc,time_nc,varname_time_nc,varptr_time_nc,&
       ierr)
    CALL read_var_nf77_2d(filename_nc,dim_wl_etfphot,time_nc,varname_v1_nc,varptr_2d_v1_nc,&
       ierr)
    CALL read_var_nf77_2d(filename_nc,dim_wl_ms93,time_nc,varname_v2_nc,varptr_2d_v2_nc,&
       ierr)
    CALL read_var_nf77_2d(filename_nc,dim_wl_rsf,time_nc,varname_v3_nc,varptr_2d_v3_nc,&
       ierr)
    IF (ierr /= 0) CALL finish ('read_dailyetf ', 'error reading UV netcdf file')
  ENDIF

  IF (p_parallel) THEN
    CALL p_bcast(varptr_time_nc,p_io)
    CALL p_bcast(varptr_2d_v1_nc,p_io)
    CALL p_bcast(varptr_2d_v2_nc,p_io)
    CALL p_bcast(varptr_2d_v3_nc,p_io)
  ENDIF

  ALLOCATE(etfphot_year(numwl_etf,dimlength))
  ALLOCATE(etfphot_ms93_year(numwl_etf_ms,dimlength))
  ALLOCATE(rsf_sclfac_year(numwl_rsf,dimlength))
  ALLOCATE(etftime_year(dimlength))

  etfphot_year(:,:)=varptr_2d_v1_nc(:,:)
  etfphot_ms93_year(:,:)=varptr_2d_v2_nc(:,:)
  rsf_sclfac_year(:,:)=varptr_2d_v3_nc(:,:)
  etftime_year(:)=varptr_time_nc(:)

  DEALLOCATE (varptr_2d_v1_nc, varptr_2d_v2_nc, varptr_2d_v3_nc, varptr_time_nc)


END SUBROUTINE read_dailyetf

SUBROUTINE set_dailyetf  ! called by src/stepon and  src/control

  USE mo_time_control,  ONLY: next_date, get_date_components

  INTEGER :: inextmonth, inextyear, inextday, r_c_date, lt


  CALL get_date_components (next_date,  month=inextmonth, year=inextyear,          &
                              day=inextday )

  r_c_date = inextyear*10000 + inextmonth*100 + inextday

  DO lt=1,dimlength

    IF (etftime_year(lt) .EQ. r_c_date) THEN

      etfphot_daily(:,:)=etfphot_year(:,lt-1:lt+1)
      etfphot_ms93_daily(:,:)=etfphot_ms93_year(:,lt-1:lt+1)
      rsf_sclfac_daily(:,:)=rsf_sclfac_year(:,lt-1:lt+1)
      write(message_text,'(a,i8)') 'set_dailyetf: etfphot(:,2) for ',r_c_date
      CALL message('', message_text)
      write(message_text,'(1p,10g15.7)') etfphot_daily(1:10,2)
      CALL message('', message_text)
      write(message_text,'(1p,10g15.7)') etfphot_daily(11:20,2)
      CALL message('', message_text)
      write(message_text,'(1p,10g15.7)') etfphot_daily(21:30,2)
      CALL message('', message_text)
      write(message_text,'(1p,3g15.7)') etfphot_daily(31:33,2)
      CALL message('', message_text)
      write(message_text,'(a,i8)') 'set_dailyetf: etfphot_ms93(:,2) for ',r_c_date
      CALL message('', message_text)
      write(message_text,'(1p,5g15.7)') etfphot_ms93_daily(:,2)
      CALL message('', message_text)
      write(message_text,'(a,i8)') 'set_dailyetf: rsf_sclfac(1:20,2) for ',r_c_date
      CALL message('', message_text)
      write(message_text,'(1p,10g15.7)') rsf_sclfac_daily(1:10,2)
      CALL message('', message_text)
      write(message_text,'(1p,10g15.7)') rsf_sclfac_daily(11:20,2)
      CALL message('', message_text)
      EXIT

    ENDIF

  ENDDO


END SUBROUTINE set_dailyetf

END MODULE mo_moz_dailyetf
