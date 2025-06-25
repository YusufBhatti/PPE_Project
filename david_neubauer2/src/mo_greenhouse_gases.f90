!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_greenhouse_gases

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_ghg
  PUBLIC :: interpolate_ghg, interpolate_co2
  PUBLIC :: cleanup_greenhouse_gases

  PUBLIC :: ghg_no_cfc
  PUBLIC :: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcvmr

  INTEGER, PARAMETER :: ghg_no_cfc = 2
  CHARACTER(len=*), PARAMETER :: ghg_cfc_names(ghg_no_cfc) = (/ &
       "CFC_11", "CFC_12" /)

  INTEGER :: ghg_no_years

  REAL(dp), ALLOCATABLE :: ghg_years(:)
  REAL(dp), ALLOCATABLE :: ghg_co2(:)
  REAL(dp), ALLOCATABLE :: ghg_ch4(:)
  REAL(dp), ALLOCATABLE :: ghg_n2o(:)
  REAL(dp), ALLOCATABLE :: ghg_cfc(:,:)

  REAL(dp) :: ghg_base_year

  REAL(dp) :: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr
  REAL(dp) :: ghg_cfcvmr(ghg_no_cfc)

CONTAINS

  SUBROUTINE init_ghg(ighg)

    ! Description:
    !
    ! Read time series of greenhouse gases
    !
    ! Method:
    !
    ! Time series of various greenhouse gases are read from
    ! file greenhouse_gases.nc (CO2, CH4, N2O, and CFC's).
    ! 
    ! Authors:
    !
    ! U. Schlese, DKRZ, June 1995, original source
    ! L. Kornblueh, MPI, November 2001, changed to read netCDF input,
    !             packed in a module, f90 rewrite, and parallelization        
    !
    ! for more details see file AUTHORS
    !

    USE mo_netCDF,    ONLY: FILE_INFO, &
                            IO_inq_dimid, IO_inq_dimlen, &
                            IO_inq_varid, IO_get_var_double
    USE mo_io,        ONLY: IO_open, IO_READ, IO_close
    USE mo_mpi,       ONLY: p_parallel_io, p_io, p_bcast

    INTEGER, INTENT (IN) :: ighg

    INTEGER :: nvarid, ndimid
    INTEGER :: i

    TYPE (FILE_INFO) :: ghg

    IF (p_parallel_io) THEN 
      ghg%opened=.FALSE.
    ! select scenario
      IF(ighg == 1) CALL IO_open('greenhouse_gases.nc', ghg, IO_READ)

      CALL IO_inq_dimid (ghg%file_id, 'time', ndimid) 
      CALL IO_inq_dimlen (ghg%file_id, ndimid, ghg_no_years) 
    ENDIF

    CALL p_bcast(ghg_no_years, p_io)

    ALLOCATE (ghg_years(ghg_no_years))
    ALLOCATE (ghg_co2(ghg_no_years))
    ALLOCATE (ghg_ch4(ghg_no_years))
    ALLOCATE (ghg_n2o(ghg_no_years))
    ALLOCATE (ghg_cfc(ghg_no_years,ghg_no_cfc))

    IF (p_parallel_io) THEN 
      CALL IO_inq_varid(ghg%file_id, 'time', nvarid)
      CALL IO_get_var_double (ghg%file_id, nvarid, ghg_years)      

      CALL IO_inq_varid (ghg%file_id, 'CO2', nvarid)
      CALL IO_get_var_double (ghg%file_id, nvarid, ghg_co2)

      CALL IO_inq_varid (ghg%file_id, 'CH4', nvarid)
      CALL IO_get_var_double (ghg%file_id, nvarid, ghg_ch4)

      CALL IO_inq_varid (ghg%file_id, 'N2O', nvarid)
      CALL IO_get_var_double (ghg%file_id, nvarid, ghg_n2o)

      DO i = 1, ghg_no_cfc
        CALL IO_inq_varid (ghg%file_id, TRIM(ghg_cfc_names(i)), nvarid)
        CALL IO_get_var_double (ghg%file_id, nvarid, ghg_cfc(:,i))
      ENDDO

      CALL IO_close(ghg)
    END IF

    CALL p_bcast(ghg_years, p_io)
    CALL p_bcast(ghg_co2, p_io)
    CALL p_bcast(ghg_ch4, p_io)
    CALL p_bcast(ghg_n2o, p_io)
    CALL p_bcast(ghg_cfc, p_io)

    ghg_base_year = ghg_years(1)

  END SUBROUTINE init_ghg

  SUBROUTINE interpolate_ghg
    
    USE mo_exception,       ONLY: message, message_text
    USE mo_physical_constants, ONLY: amd, amco2, amch4, amn2o
    USE mo_time_conversion, ONLY: time_native, &
                                  TC_get, TC_convert, &
                                  year_len, day_in_year, &
                                  print_date
    USE mo_time_control,    ONLY: radiation_date, NDAYLEN
    USE mo_radiation_parameters, ONLY: fco2

    TYPE(time_native) :: ghg_date 

    REAL(dp) :: zsecref, zsecnow
    REAL(dp) :: zw1, zw2
    REAL(dp) :: zco2int, zch4int, zn2oint
    REAL(dp) :: zcfc(ghg_no_cfc)

    INTEGER :: iyear, iyearm, iyearp
    INTEGER :: yr, mo, dy, hr, mn, se
    INTEGER :: day, second
    
    CHARACTER(len=32)  :: cdate, cformat
    CHARACTER(len=256) :: ccfc

    ! Description:
    !
    ! Read time series of greenhouse gases
    !
    ! Method:
    !
    ! 1. Interpolation in time.
    ! 2. Convert from volume mixing ratio to mass mixing ratio of 
    !    CO2, CH4, and N2O - not for CFC's!
    ! 
    ! Authors:
    !
    ! U. Schlese, DKRZ, June 1995, original source
    ! L. Kornblueh, MPI, November 2001, changed to read netCDF input,
    !             packed in a module, f90 rewrite, and parallelization        
    ! M. Esch,    MPI, MAY 2004, modified for scenarios
    ! M. Esch,    MPI, DEC 2009, modified for CMIP5
    ! J.S. Rast,  MPI, AUG 2010, modified interpolation to time step of radiation
    !
    ! for more details see file AUTHORS
    !

    ! 1. Interpolation in time.

    CALL TC_convert(radiation_date, ghg_date) 
    CALL TC_get (ghg_date, yr, mo, dy, hr, mn, se)
    CALL TC_get (radiation_date, day, second)

    zsecref = year_len(ghg_date)*REAL(NDAYLEN,dp)
    zsecnow = (day_in_year(radiation_date)-1)*REAL(NDAYLEN,dp)+second 

    iyear = INT(yr)-ghg_base_year+1   ! set right index to access in ghg fields
    iyearm = iyear-1
    iyearp = iyear+1

    IF(mo <= 6) THEN     ! First half of year
      zw1 = zsecnow/zsecref+0.5_dp
      zw2 = 1.0_dp-zw1

      zco2int = zw1*ghg_co2(iyear)+zw2*ghg_co2(iyearm)
      zch4int = zw1*ghg_ch4(iyear)+zw2*ghg_ch4(iyearm)
      zn2oint = zw1*ghg_n2o(iyear)+zw2*ghg_n2o(iyearm)
      zcfc(:) = zw1*ghg_cfc(iyear,:)+zw2*ghg_cfc(iyearm,:)
    ELSE                 ! Second half of year
      zw2= zsecnow/zsecref-0.5_dp
      zw1= 1.0_dp-zw2

      zco2int = zw1*ghg_co2(iyear)+zw2*ghg_co2(iyearp)
      zch4int = zw1*ghg_ch4(iyear)+zw2*ghg_ch4(iyearp)
      zn2oint = zw1*ghg_n2o(iyear)+zw2*ghg_n2o(iyearp)
      zcfc(:) = zw1*ghg_cfc(iyear,:)+zw2*ghg_cfc(iyearp,:)
    END IF

    IF (ABS(fco2-1._dp) > EPSILON(1._dp)) zco2int = fco2*zco2int

    ! Write greenhouse gas time series
    
    !IF (p_parallel_io) THEN
    !  WRITE(8) get_time_step(),zco2int,zch4int,zn2oint,(zcfc(jc),jc=1,2)
    !END IF

    CALL print_date(ghg_date, mess=cdate)
    WRITE(cformat,'(a,i0,a)') '(a,', ghg_no_cfc, 'f7.2)'
    WRITE(ccfc,cformat) ' CFC = ', zcfc(1:ghg_no_cfc)
    WRITE (message_text,'(a,a,a,f9.4,a,f8.2,a,f8.3,a)') &
         'Greenhouse gas concentrations ', TRIM(cdate), ' CO2 = ', zco2int, &
         ' CH4 = ', zch4int,' N2O = ', zn2oint, TRIM(ccfc)
    CALL message('', TRIM(message_text))

    ! consistency check

    ! 2. Convert to mass mixing ratio.

    ghg_co2mmr    = zco2int*1.0e-06_dp*amco2/amd 
    ghg_ch4mmr    = zch4int*1.0e-09_dp*amch4/amd
    ghg_n2ommr    = zn2oint*1.0e-09_dp*amn2o/amd

    ! Scale CFCs only, keep the volume mixing ratio 

    ghg_cfcvmr(:) = zcfc(:)*1.0e-12_dp

  END SUBROUTINE interpolate_ghg

  SUBROUTINE interpolate_co2
    
    USE mo_physical_constants, ONLY: amd, amco2 
    USE mo_time_conversion, ONLY: time_native, &
                                  TC_get, TC_convert, &
                                  year_len, day_in_year
    USE mo_time_control,    ONLY: current_date, NDAYLEN
    USE mo_co2,             ONLY: co2m1

    TYPE(time_native) :: ghg_date 

    REAL(dp) :: zsecref, zsecnow
    REAL(dp) :: zw1, zw2
    REAL(dp) :: zco2int

    INTEGER :: iyear, iyearm, iyearp
    INTEGER :: yr, mo, dy, hr, mn, se
    INTEGER :: day, second
    
    ! Description:
    !
    ! Get CO2 concentration for current time step
    !
    ! Method:
    !
    ! 1. Interpolation in time.
    ! 2. Convert from volume mixing ratio to mass mixing ratio of 
    !    CO2
    ! 
    ! Authors:
    !
    ! R. Schnur,  MPI, NOV 2010, copy of interpolate_ghg, but for current time step
    !                            and CO2 only (for use in mo_co2)
    !
    ! for more details see file AUTHORS
    !

    ! 1. Interpolation in time.

    CALL TC_convert(current_date, ghg_date) 
    CALL TC_get (ghg_date, yr, mo, dy, hr, mn, se)
    CALL TC_get (current_date, day, second)

    zsecref = year_len(ghg_date)*REAL(NDAYLEN,dp)
    zsecnow = (day_in_year(current_date)-1)*REAL(NDAYLEN,dp)+second 

    iyear = INT(yr)-ghg_base_year+1   ! set right index to access in ghg fields
    iyearm = iyear-1
    iyearp = iyear+1

    IF(mo <= 6) THEN     ! First half of year
      zw1 = zsecnow/zsecref+0.5_dp
      zw2 = 1.0_dp-zw1

      zco2int = zw1*ghg_co2(iyear)+zw2*ghg_co2(iyearm)
    ELSE                 ! Second half of year
      zw2= zsecnow/zsecref-0.5_dp
      zw1= 1.0_dp-zw2

      zco2int = zw1*ghg_co2(iyear)+zw2*ghg_co2(iyearp)
    END IF

    ! 2. Convert to mass mixing ratio.

    co2m1 = zco2int*1.0e-06_dp*amco2/amd 

  END SUBROUTINE interpolate_co2

  SUBROUTINE cleanup_greenhouse_gases

    IF (ALLOCATED(ghg_years))   DEALLOCATE(ghg_years)
    IF (ALLOCATED(ghg_co2))     DEALLOCATE(ghg_co2)
    IF (ALLOCATED(ghg_ch4))     DEALLOCATE(ghg_ch4)
    IF (ALLOCATED(ghg_n2o))     DEALLOCATE(ghg_n2o)
    IF (ALLOCATED(ghg_cfc))     DEALLOCATE(ghg_cfc)

  END SUBROUTINE cleanup_greenhouse_gases

END MODULE mo_greenhouse_gases
