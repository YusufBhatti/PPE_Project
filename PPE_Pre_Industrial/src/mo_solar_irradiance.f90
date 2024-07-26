!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! mo_solar_irradiance: spectrally resolved time varying solar irradiance
!!
!! Read spectrally resolved solar irradiance yearly, apply primitive time
!! interpolation to monthly mean values and apply
!!
!! @author Sebastian Rast, MPI Met, Hamburg
!!
!! $ID: n/a$
!!
!! @par Revision History
!! original source by J.S.Rast (2010-03-22)
!!
MODULE mo_solar_irradiance
  USE mo_kind,                 ONLY: dp
  USE mo_exception,            ONLY: finish

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC        :: init_solar_irradiance, cleanup_solar_irradiance, &
                   get_solar_irradiance, get_solar_irradiance_m,    &
                   set_solar_irradiance, set_solar_irradiance_m

  REAL(dp), ALLOCATABLE           :: ssi_m(:,:) !< ssi: spectrally
            !< resolved time dependent ssi, ssi_m: for interpolation 
            !< in the case of l_trigrad=.true., may be different years because
            !< radiation_date is in the future compared to current_date
  REAL(dp), ALLOCATABLE           :: tsi(:), tsi_m(:) !< total solar irradiance
            !< compare with ssi
  LOGICAL                         :: ltsi_set=.false., ltsi_set_m=.false.
CONTAINS
!-----------------------------------------------------------------------------
!>
!! su_solar_irradiance: set up memory for time varying solar radiation
!!
!! @par Revision History
!! original source by J.S.Rast (2010-03-23)
!-----------------------------------------------------------------------------
SUBROUTINE init_solar_irradiance(nb_sw)
  INTEGER, INTENT(in)               :: nb_sw
  INTEGER, PARAMETER                :: nmonths=12
  
  ALLOCATE(tsi(0:nmonths+1))
  ALLOCATE(ssi_m(nb_sw,0:nmonths+1))
  ALLOCATE(tsi_m(0:nmonths+1))
END SUBROUTINE init_solar_irradiance
!-----------------------------------------------------------------------------
!>
!! get_solar_irradiance: yearly reading of spectrally dependent
!! solar irradiance
!!
!! @par Revision History
!! original source by J.S.Rast (2010-03-24)
!-----------------------------------------------------------------------------
SUBROUTINE get_solar_irradiance(currentdate, nextdate)

  !USES
  USE mo_time_conversion,    ONLY: time_days
  USE mo_time_control,       ONLY: get_date_components, next_date

  !INPUT PARAMETERS
  TYPE(time_days), INTENT(in), OPTIONAL :: currentdate, nextdate

  !LOCAL VARIABLES
  LOGICAL                       :: lnewyear
  INTEGER                       :: icurrentyear, inextyear
  INTEGER                       :: iyrm1, iyr, iyrp1
  CHARACTER(len=20)             :: cfname_base,cyr
  CHARACTER(len=25)             :: cfname
  INTEGER                       :: yr, mo, dy, hr, mn, se

  IF (PRESENT(currentdate) .AND. PRESENT(nextdate)) THEN
    CALL get_date_components (currentdate, year=icurrentyear)
    CALL get_date_components (nextdate, year=inextyear)
    lnewyear=icurrentyear/=inextyear
    iyrm1 = inextyear-1
    iyr   = inextyear
    iyrp1 = inextyear+1
  ELSE
    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)
    iyrm1 = yr-1
    iyr   = yr
    iyrp1 = yr+1
  ENDIF

  IF (lnewyear .OR. .NOT. ltsi_set) THEN 
! irradiance file
     cfname_base='swflux'
     WRITE(cyr,'(i4.4)') iyrm1
     cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
     CALL read_solar_irradiance(tsi(0:0),    'TSI',         12,             &
                                12,          'time',        cfname          )
     WRITE(cyr,'(i4.4)') iyr
     cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
     CALL read_solar_irradiance(tsi(1:12),   'TSI',          1,             &
                                12,          'time',         cfname         )
     WRITE(cyr,'(i4.4)') iyrp1
     cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
     CALL read_solar_irradiance(tsi(13:13),  'TSI',          1,             &
                                1,           'time',         cfname         )
     ltsi_set=.true.
  END IF
END SUBROUTINE get_solar_irradiance
!-----------------------------------------------------------------------------
!>
!! get_solar_irradiance_m: yearly reading of spectrally dependent
!! solar irradiance for l_trigrad=.true. (radiation time step)
!!
!! @par Revision History
!! original source by J.S.Rast (2010-03-24)
!-----------------------------------------------------------------------------
SUBROUTINE get_solar_irradiance_m(rad_step_1, rad_step_2, nb_sw)

  !USES
  USE mo_time_conversion,    ONLY: time_days
  USE mo_time_control,       ONLY: get_date_components

  !INPUT PARAMETERS
  TYPE(time_days), INTENT(in)   :: rad_step_1, rad_step_2
  INTEGER, INTENT(in)           :: nb_sw

  !LOCAL VARIABLES
  LOGICAL                       :: lnewyear
  INTEGER                       :: icurrentyear, inextyear
  INTEGER                       :: iyrm1, iyr, iyrp1
  CHARACTER(len=20)             :: cfname_base,cyr
  CHARACTER(len=25)             :: cfname
  CALL get_date_components (rad_step_1, year=icurrentyear)
  CALL get_date_components (rad_step_2, year=inextyear)
  lnewyear=icurrentyear/=inextyear

  IF (lnewyear .OR. .NOT. ltsi_set_m) THEN 
     iyrm1=inextyear-1
     iyr=inextyear
     iyrp1=inextyear+1
! irradiance file
     cfname_base='swflux'
     WRITE(cyr,'(i4.4)') iyrm1
     cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
     CALL read_solar_irradiance_m(tsi_m(0:0),    'TSI',       ssi_m(:,0:0), &
                                  'SSI',         nb_sw,       'numwl',      &
                                  12,            12,          'time',       &
                                  cfname                                    )
     WRITE(cyr,'(i4.4)') iyr
     cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
     CALL read_solar_irradiance_m(tsi_m(1:12),   'TSI',       ssi_m(:,1:12),&
                                  'SSI',          nb_sw,      'numwl',      &
                                  1,              12,         'time',       &
                                  cfname                                    )
     WRITE(cyr,'(i4.4)') iyrp1
     cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
     CALL read_solar_irradiance_m(tsi_m(13:13),   'TSI',      ssi_m(:,13:13),&
                                  'SSI',          nb_sw,      'numwl',       &
                                  1,              1,          'time',        &
                                  cfname                                     )
     ltsi_set_m=.true.
  END IF
END SUBROUTINE get_solar_irradiance_m
!-----------------------------------------------------------------------------
!>
!! set_solar_irradiance: interpolate total solar irradiance 
!! (for l_trigrad=.false.)
!! solar irradiance
!!
!! @par Revision History
!! original source by J.S.Rast (2010-03-29)
!-----------------------------------------------------------------------------
SUBROUTINE set_solar_irradiance(ptsi_i)

  !USES
  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2

  !INPUT PARAMETERS
  REAL(dp), INTENT(out)   :: ptsi_i !time interpolated total solar irradiance

  ptsi_i=wgt1*tsi(nmw1)+wgt2*tsi(nmw2)
END SUBROUTINE set_solar_irradiance
!-----------------------------------------------------------------------------
!>
!! set_solar_irradiance: interpolate and set spectrally dependent 
!! solar irradiance
!!
!! @par Revision History
!! original source by J.S.Rast (2010-03-29)
!-----------------------------------------------------------------------------
SUBROUTINE set_solar_irradiance_m(ptsi_i, pssi_i, nb_sw)

  !USES
  USE mo_interpo,       ONLY: nmw1_m, nmw2_m, wgt1_m, wgt2_m

  !INPUT PARAMETERS
  INTEGER, INTENT(in)     :: nb_sw    !number of solar wave length bands
  REAL(dp), INTENT(out)   :: ptsi_i !time interpolated total solar irradiance
  REAL(dp), INTENT(out)   :: pssi_i(nb_sw) !time interpolated solar irradiance
  
  ptsi_i=wgt1_m*tsi_m(nmw1_m)+wgt2_m*tsi_m(nmw2_m)
  pssi_i(1:nb_sw)=wgt1_m*ssi_m(1:nb_sw,nmw1_m)+wgt2_m*ssi_m(1:nb_sw,nmw2_m)
END SUBROUTINE set_solar_irradiance_m
!-----------------------------------------------------------------------------
!>
!! solar_irradiance_read: handle netcdf-routines for reading
!!
!! @par Revision History
!! original source by J.S.Rast (2010-03-24)
!-----------------------------------------------------------------------------
SUBROUTINE read_solar_irradiance(ptsi,       cptsi,             imnthb,     &
                                 imnthe,     ctime,             cfname      )
  
  !USES
  USE mo_mpi,                      ONLY: p_parallel_io, p_io, p_bcast
  USE mo_read_netcdf77,            ONLY: read_var_hs_nf77_0d

  !INPUT PARAMETERS
  INTEGER, INTENT(in)           :: imnthb, imnthe ! begin and end month
  REAL(dp), INTENT(out)         :: ptsi(imnthb:imnthe) ! total solar irradiance
  CHARACTER(len=*), INTENT(in)  :: cptsi,ctime,cfname 
                                   ! names of variables in netcdf file

  !LOCAL VARIABLES
  LOGICAL                       :: lex
  INTEGER                       :: j, ierr
  
  IF (p_parallel_io) THEN
     INQUIRE (file=TRIM(cfname), exist=lex)
     IF (.NOT. lex) THEN
        CALL finish('read_solar_irradiance','file '//TRIM(cfname)// &
                    ' does not exist')
     END IF
     DO j=imnthb,imnthe
        CALL read_var_hs_nf77_0d (TRIM(cfname), TRIM(ctime), &
                                  j, TRIM(cptsi), ptsi(j), ierr)
     END DO
  END IF
  CALL p_bcast(ptsi(imnthb:imnthe), p_io)
END SUBROUTINE read_solar_irradiance
!-----------------------------------------------------------------------------
!>
!! read_solar_irradiance_m: handle netcdf-routines for reading
!!
!! @par Revision History
!! original source by J.S.Rast (2010-03-24)
!-----------------------------------------------------------------------------
SUBROUTINE read_solar_irradiance_m(ptsi,       cptsi,             pssi,       &
                                   cpssi,      nwl,               cnwl,       &
                                   imnthb,     imnthe,            ctime,      &
                                   cfname                                     )
  
  !USES
  USE mo_mpi,                      ONLY: p_parallel_io, p_io, p_bcast
  USE mo_read_netcdf77,            ONLY: read_diml_nf77, read_var_hs_nf77_0d, &
                                         read_var_hs_nf77_1d

  !INPUT PARAMETERS
  INTEGER, INTENT(in)           :: nwl ! number of wave length bands
  INTEGER, INTENT(in)           :: imnthb, imnthe ! begin and end month
  REAL(dp), INTENT(out)         :: ptsi(imnthb:imnthe) ! total solar irradiance
  REAL(dp), INTENT(out)         :: pssi(nwl,imnthb:imnthe) 
                                   ! spectrally resolved solar irradiance
  CHARACTER(len=*), INTENT(in)  :: cptsi,cpssi,cnwl,ctime,cfname 
                                   ! names of variables in netcdf file

  !LOCAL VARIABLES
  LOGICAL                       :: lex
  INTEGER                       :: inwl, j, ierr
  
  IF (p_parallel_io) THEN
     INQUIRE (file=TRIM(cfname), exist=lex)
     IF (.NOT. lex) THEN
        CALL finish('read_solar_irradiance_m','file '//TRIM(cfname)// &
                    ' does not exist')
     END IF
     DO j=imnthb,imnthe
        CALL read_var_hs_nf77_0d (TRIM(cfname), TRIM(ctime), &
                                  j, TRIM(cptsi), ptsi(j), ierr)
     END DO
     inwl=read_diml_nf77(TRIM(cfname),TRIM(cnwl))
     IF (inwl /= nwl) THEN
        CALL finish('aero_read_opt_volc_sw', &
             'incompatible number of wavelengths in file '//TRIM(cfname))
     END IF
     DO j=imnthb,imnthe
        CALL read_var_hs_nf77_1d(TRIM(cfname),TRIM(cnwl),TRIM(ctime),j, &
             TRIM(cpssi),pssi(:,j),ierr)        
     END DO
  END IF
  CALL p_bcast(ptsi(imnthb:imnthe), p_io)
  CALL p_bcast(pssi(:,imnthb:imnthe), p_io)
END SUBROUTINE read_solar_irradiance_m
!-----------------------------------------------------------------------------
!>
!! cleanup_solar_irradiance: deallocate allocated memory
!!
!! @par Revision History
!! original source by J.S.Rast (2011-01-31)
!-----------------------------------------------------------------------------
SUBROUTINE cleanup_solar_irradiance

  IF (ALLOCATED(tsi))   DEALLOCATE(tsi)
  IF (ALLOCATED(ssi_m)) DEALLOCATE(ssi_m)
  IF (ALLOCATED(tsi_m)) DEALLOCATE(tsi_m)

  ltsi_set   = .FALSE.
  ltsi_set_m = .FALSE.

END SUBROUTINE cleanup_solar_irradiance
END MODULE mo_solar_irradiance
