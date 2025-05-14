!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! Module that reads and processes external fields in NetCDF format
!!
!! concepts of the routines:
!! see also: http://hammoz.icg.fz-juelich.de/data/BoundaryConditions
!!
!! @author S. Schroeder, FZ-Juelich
!!
!! $Id: 1423$
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
!! @par Copyright
!! 2009 by MPI-M and FZJ
!! This software is provided for non-commercial use only.
!!
MODULE mo_external_field_processor
  USE mo_kind,            ONLY: dp
  USE mo_time_conversion, ONLY: time_days

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ef_get_next_timestep, ef_get_first_timestep
  PUBLIC :: ef_read_single, ef_read_1d, ef_read_2d, ef_read_3d

  PUBLIC :: EF_INACTIVE, EF_VALUE, EF_FILE, EF_MODULE
  PUBLIC :: EF_UNDEFINED, EF_3D, EF_LONLAT, EF_LATLEV, EF_LEV, EF_LAT, EF_SINGLE
  PUBLIC :: EF_TIMERESOLVED, EF_IGNOREYEAR, EF_CONSTANT
  PUBLIC :: EF_NOINTER, EF_LINEAR, EF_CUBIC

! flags of ef_type

  INTEGER, PARAMETER :: EF_INACTIVE = 0
  INTEGER, PARAMETER :: EF_VALUE    = 1
  INTEGER, PARAMETER :: EF_FILE     = 2
  INTEGER, PARAMETER :: EF_MODULE   = 3

! flags of ef_geometry

  INTEGER, PARAMETER :: EF_UNDEFINED = 0
  INTEGER, PARAMETER :: EF_3D        = 1
  INTEGER, PARAMETER :: EF_LONLAT    = 2
  INTEGER, PARAMETER :: EF_LATLEV    = 3
  INTEGER, PARAMETER :: EF_LEV       = 4
  INTEGER, PARAMETER :: EF_LAT       = 5
  INTEGER, PARAMETER :: EF_SINGLE    = 6

! flags of ef_timedef

  INTEGER, PARAMETER :: EF_TIMERESOLVED = 1
  INTEGER, PARAMETER :: EF_IGNOREYEAR   = 2
  INTEGER, PARAMETER :: EF_CONSTANT     = 3

! flags of ef_interpolate

  INTEGER, PARAMETER :: EF_NOINTER = 0
  INTEGER, PARAMETER :: EF_LINEAR  = 1
  INTEGER, PARAMETER :: EF_CUBIC   = 2

  TYPE, PUBLIC :: external_field   ! public only to mo_boundary_condition_processor
                                   ! type of external field
    INTEGER                      :: ef_type = EF_INACTIVE, &    ! EF_MODULE, EF_FILE, EF_VALUE, EF_INACTIVE
                                                                ! information for type='file'
                                    ef_nzval, &                 ! number of levels in input file (not to be set by user!)
                                    ef_geometry=-1, &           ! EF_3D, EF_LONLAT, EF_LATLEV, EF_LEV, EF_LAT, EF_SINGLE
                                    ef_timedef=1, &             ! EF_TIMERESOLVED, EF_IGNOREYEAR, EF_CONSTANT
                                    ef_timeindex = 0, &         ! time record in ef file
                                    ef_timeindex_prior = 0, &   ! (for time interpolation:
                                                                ! prior time step might be from a previous file)
                                    ef_interpolate = EF_NOINTER ! no interpolation in time
    LOGICAL                      :: ef_ltpredict=.false.        ! is the next timestep predictable from the actual one?
                                                                ! (important for the last timestep in a file)
    REAL(dp)                     :: ef_timeoffset = 0.0_dp, &   ! offset value for file time
                                                                ! (unit: unit of time record in input file)
                                    ef_timetol = 2.0_dp, &      ! ++mgs tolerance value for time comparison (later user-defined?)
                                    ef_value = 0._dp, &         ! information for type='value'
                                    ef_factor = 1.0_dp          ! scaling factor
    REAL(dp), ALLOCATABLE        :: ef_zvalueh(:)               ! half levels either in height[m]    (descending order)
                                                                !             or     in pressure[Pa] (ascending order)
                                                                ! (not to be set by user!)
    CHARACTER(LEN=512)           :: &
                                    ef_file= '', &              ! filename
                                    ef_template= '', &          ! template for filename
                                    ef_varname= ''              ! variable name (user must define it by namelist!)
    CHARACTER(LEN=30)            :: ef_actual_unit=''
    TYPE(time_days), ALLOCATABLE :: ef_times(:)                 ! not to be set by user!
    TYPE(time_days), ALLOCATABLE :: ef_times_prior(:)           ! (for time interpolation:
                                                                ! prior time step might be from a previous file)
  END TYPE external_field

  REAL(dp), ALLOCATABLE  :: value1d(:)
  REAL(dp), ALLOCATABLE  :: value2d(:,:)
  REAL(dp), ALLOCATABLE  :: value3d(:,:,:)
  PUBLIC :: value1d, value2d, value3d    ! to let these arrays be safely allocated at run time

  ! subprograms

  CONTAINS

!-----------------------------------------------------------------------
!>
!! return one variable's description from a list of variables by index
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
! FUNCTION ef_get_element_from_list(efield_varlist,index) RESULT (efield_element)
! TYPE(ext_field_nml), INTENT(in) :: efield_varlist
! INTEGER,             INTENT(in) :: index
! 
! TYPE(external_field) :: efield_element
!
! efield_element%ef_type        = efield_varlist%ef_type
! efield_element%ef_template    = efield_varlist%ef_template
! efield_element%ef_varname     = efield_varlist%ef_varlist(index)
! efield_element%ef_geometry    = efield_varlist%ef_geometry
! efield_element%ef_timedef     = efield_varlist%ef_timedef
! efield_element%ef_ltpredict   = efield_varlist%ef_ltpredict
! efield_element%ef_timeoffset  = efield_varlist%ef_timeoffset
! efield_element%ef_timeindex   = efield_varlist%ef_timeindex
! efield_element%ef_interpolate = efield_varlist%ef_interpolate
! efield_element%ef_value       = efield_varlist%ef_value
! efield_element%ef_factor      = efield_varlist%ef_factor
! efield_element%ef_actual_unit = efield_varlist%ef_actual_unit
!
! END FUNCTION ef_get_element_from_list

!>
!! consistency check of an external field description
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  SUBROUTINE ef_check (efield, ncid)
  ! should be CALLed by p_io
  ! values should be broadcasted to other PEs after the CALL

  USE mo_netcdf
  USE mo_exception, ONLY: message, em_error

  TYPE(external_field), INTENT(inout) :: efield
  INTEGER, INTENT(in)                 :: ncid

  character                           :: cinter
  INTEGER                             :: dimid, dimlen, varid, status

  IF (trim(efield%ef_varname) .EQ. '') THEN
    CALL message('ef_check', 'ef_varname has to be defined!',level=em_error)
  ELSE
    status = nf_inq_varid( ncid, efield%ef_varname, varid )
    IF (status /= nf_noerr) THEN
      CALL message('ef_check', trim(efield%ef_varname)//' does not exist in file '// &
                                         TRIM(efield%ef_file)//'!',level=em_error)
    ELSE
      IF (.not. ef_check_var(efield,ncid)) CALL message('ef_check', &
          'wrong dimensionality for '//TRIM(efield%ef_varname)//' in '//TRIM(efield%ef_file),level=em_error)
    ENDIF
  ENDIF
  SELECT CASE (efield%ef_geometry)
     CASE (EF_3D, EF_LONLAT, EF_LATLEV, EF_LEV, EF_LAT, EF_SINGLE)
       IF (.NOT. ef_check_geometry(efield,ncid)) CALL message('ef_check', &
              '  missing dimension or wrong dimension length found in '//TRIM(efield%ef_file),level=em_error)
     CASE DEFAULT
       CALL message('ef_check', 'not a defined external_field_geometry given!',level=em_error)
  END SELECT
  SELECT CASE (efield%ef_timedef)
     CASE (EF_TIMERESOLVED, EF_IGNOREYEAR, EF_CONSTANT)
       IF (efield%ef_timedef .eq. EF_CONSTANT) THEN

  ! check whether ef_timeindex is set

         IF (efield%ef_timeindex .eq. 0) THEN

  ! better handled in other routines if set!

           efield%ef_timeindex = 1
           CALL nf_check(nf_inq_dimid  (ncid, "time", dimid))
           CALL nf_check(nf_inq_dimlen (ncid, dimid, dimlen))
           IF (dimlen /= 1) CALL message('ef_check', 'no timeindex (of multiple timesteps) '// &
                                           'for input file defined but timedef "constant" is given!',level=em_error)
         ENDIF
       ENDIF
     CASE DEFAULT
       CALL message('ef_check', 'not a defined external_field_timedef given!',level=em_error)
  END SELECT
  SELECT CASE (efield%ef_interpolate)
     CASE (0, 1, 2)
     CASE DEFAULT
       write(cinter,'(i1)') efield%ef_interpolate
       CALL message('ef_check', cinter//' is not a defined external_field_interpolation!',level=em_error)
  END SELECT
  IF (.NOT. ef_check_units(efield)) CALL message('ef_check', 'unit mismatch in field definition and file!',level=em_error)
  END SUBROUTINE ef_check

!>
!! check whether the variable of an external field corresponds to the given dimensions
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
!++mgs (14.06.2011): ef_get_geometry recursive to allow for automatic determination of geometry
  RECURSIVE FUNCTION ef_get_geometry(efield,ncid,jgeom) RESULT (found)
!--mgs

  USE mo_netcdf
!++mgs
  USE mo_exception,                ONLY: message, em_info, em_error
!--mgs

  TYPE(external_field), INTENT(inout) :: efield
  INTEGER, INTENT(in)                 :: ncid
  INTEGER, POINTER, INTENT(in)        :: jgeom

  INTEGER, PARAMETER                  :: MAXDIM=4
  INTEGER                             :: ndim
  CHARACTER*4                         :: dims(MAXDIM)
  INTEGER                             :: varid
  LOGICAL                             :: found
  CHARACTER*25                        :: varname
!++mgs
  CHARACTER(len=9), PARAMETER         :: geometry_name(6) = (/ 'EF_3D    ', &
                                                               'EF_LONLAT', &
                                                               'EF_LATLEV', &
                                                               'EF_LEV   ', &
                                                               'EF_LAT   ', &
                                                               'EF_SINGLE' /)
!--mgs

  found = .FALSE.

  SELECT CASE (efield%ef_geometry)
     CASE (EF_3D)
        ndim=3
        dims(1) = 'lon'
        dims(2) = 'lat'
        dims(3) = 'mlev'
     CASE (EF_LONLAT)
        ndim=2
        dims(1) = 'lon'
        dims(2) = 'lat'
     CASE (EF_LATLEV)
        ndim=2
        dims(1) = 'lat'
        dims(2) = 'mlev'
     CASE (EF_LEV)
        ndim=1
        dims(1) = 'mlev'
     CASE (EF_LAT)
        ndim=1
        dims(1) = 'lat'
     CASE (EF_SINGLE)
        ndim=0
!!++mgs (14.06.2011): automatically determine geometry
     CASE (EF_UNDEFINED)
        !! try out various options by calling ef_get_geometry recursively
        ndim=-1
        efield%ef_geometry = EF_3D
        found = ef_get_geometry(efield,ncid,jgeom)
        IF (.NOT. found) THEN
          efield%ef_geometry = EF_LONLAT
          found = ef_get_geometry(efield,ncid,jgeom)
        END IF
        IF (.NOT. found) THEN
          efield%ef_geometry = EF_LATLEV
          found = ef_get_geometry(efield,ncid,jgeom)
        END IF
        IF (.NOT. found) THEN
          efield%ef_geometry = EF_LEV
          found = ef_get_geometry(efield,ncid,jgeom)
        END IF
        IF (.NOT. found) THEN
          efield%ef_geometry = EF_LAT
          found = ef_get_geometry(efield,ncid,jgeom)
        END IF
        IF (.NOT. found) THEN
          efield%ef_geometry = EF_SINGLE
          found = ef_get_geometry(efield,ncid,jgeom)
        END IF
        IF (found) THEN
          jgeom = efield%ef_geometry
          CALL message('ef_get_geometry', 'Determined geometry for '//TRIM(efield%ef_varname)//  &
                       ' in '//TRIM(efield%ef_file)//' as '//                           &
                       TRIM(geometry_name(efield%ef_geometry)), level=em_info)
        ELSE
          CALL message('ef_get_geometry', 'Could not determine geometry for '//TRIM(efield%ef_varname)//  &
                       ' in '//TRIM(efield%ef_file), level=em_error)
        END IF
!!--mgs
  END SELECT
  IF ((efield%ef_timedef /= EF_CONSTANT) .AND. (ndim /= -1)) THEN
    ndim = ndim + 1
    dims(ndim) = 'time'
  ENDIF

! the given variable must be of right dimensionality (see above)

!!++mgs (14.06.2011): automatically determine geometry
  IF (.NOT. found) THEN
    CALL nf_check( nf_inq_varid( ncid, efield%ef_varname, varid ),efield%ef_file)
    found = ef_check_vardims(ncid,varid, ndim, dims)
    IF (found) jgeom = efield%ef_geometry
  END IF
!!--mgs

  END FUNCTION ef_get_geometry

  FUNCTION ef_check_var (efield, ncid) RESULT(lfound)
  TYPE(external_field), INTENT(inout) :: efield
  INTEGER, INTENT(in)                 :: ncid
  
  INTEGER, TARGET  :: igeom
  INTEGER, POINTER :: pgeom
  LOGICAL :: lfound

  igeom = EF_UNDEFINED
  pgeom => igeom
  lfound = ef_get_geometry(efield,ncid,pgeom)
  IF (lfound) THEN
    efield%ef_geometry = igeom
  ENDIF
  
  END FUNCTION ef_check_var

!>
!! check whether the variable of an external field corresponds to the given unit
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  FUNCTION ef_check_units(efield) RESULT (found)

  USE mo_netcdf

  TYPE(external_field), INTENT(in) :: efield
  LOGICAL                          :: found

  found = .true.
!!!baustelle!!!
  END FUNCTION ef_check_units

!>
!! check whether all dimensions of the external field (here: file)
!! match the given geometry (are all mandatory dimensions given and
!! are they of the same size as needed for the current run?)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  FUNCTION ef_check_geometry(efield, io_file_id) RESULT (lok)

  USE mo_exception, ONLY: message, em_info, em_error, em_debug
  USE mo_control,   ONLY  : nlev, ngl, nlon
  USE mo_gaussgrid, ONLY: philat, philon
  USE mo_netcdf

  TYPE(external_field), INTENT(inout) :: efield
  INTEGER, INTENT(in)                 :: io_file_id

  INTEGER, PARAMETER               :: MAXDIM=4
  INTEGER                          :: ndim, i, dimmodell(MAXDIM), io_dimid, io_var_id, iret
  REAL(dp)                         :: zlatvals(ngl), zlonvals(nlon), zepslat, zepslon
  REAL(dp)                         :: zlatmax, zlatmin, zlonmax, zlonmin
  CHARACTER*4                      :: dimnames(MAXDIM)
  LOGICAL                          :: lok, lcheckgrid

  lok = .TRUE.
  lcheckgrid = .FALSE.
  SELECT CASE (efield%ef_geometry)
     CASE (EF_3D)
        ndim=3
        dimnames(1) = 'mlev'
        dimmodell(1) = nlev
        dimnames(2) = 'lat'
        dimmodell(2) = ngl
        dimnames(3) = 'lon'
        dimmodell(3) = nlon
        lcheckgrid = .TRUE.
     CASE (EF_LONLAT)
        ndim=2
        dimnames(1) = 'lat'
        dimmodell(1) = ngl
        dimnames(2) = 'lon'
        dimmodell(2) = nlon
        lcheckgrid = .TRUE.
     CASE (EF_LATLEV)
        ndim=2
        dimnames(1) = 'mlev'
        dimmodell(1) = nlev
        dimnames(2) = 'lat'
        dimmodell(2) = ngl
     CASE (EF_LEV)
        ndim=1
        dimnames(1) = 'mlev'
        dimmodell(1) = nlev
     CASE (EF_LAT)
        ndim=1
        dimnames(1) = 'lat'
        dimmodell(1) = ngl
     CASE (EF_SINGLE)
        ndim=0
  END SELECT
  DO i = 1, ndim
    iret = nf_inq_dimid  (io_file_id, dimnames(i), io_var_id)
    IF ( (iret /= 0 ) .AND. ( dimnames(i) .eq. 'mlev' ) ) CALL nf_check(nf_inq_dimid  (io_file_id, 'lev', io_var_id))
    CALL nf_check(nf_inq_dimlen (io_file_id, io_var_id, io_dimid))
    !!### mgs: bug fix 09/02/2010
    IF (io_dimid /= dimmodell(i)) THEN
!sschr: The following might be true in case of e.g. aircraft data
      IF (dimnames(i) .eq. 'mlev') THEN
         !SF WRITE(message_text, '(a,i0,a,i0)') 'ef_check_geometry: dimlen(lev) = ',io_dimid,', model levels = ',dimmodell(i)
         !SF CALL message('ef_check_geometry', message_text, level=em_debug)
      ELSE
        lok = .FALSE.
      END IF
    END IF
  ENDDO

  !SF WRITE(message_text, '(a,l2,l2)') "lok,lcheckgrid=",lok,lcheckgrid
  !SF CALL message('ef_check_geometry', message_text, level=em_debug)

  IF (lok .AND. lcheckgrid) THEN
    CALL nf_check (nf_inq_varid(io_file_id, 'lat', io_var_id ),TRIM(efield%ef_file))
    CALL nf_check (nf_get_vara_double(io_file_id,io_var_id,(/1/),(/ ngl /),zlatvals),TRIM(efield%ef_file))
    zepslat = abs((zlatvals(2)-zlatvals(1)) * 0.05_dp)
    zlatmax = maxval(zlatvals(:)-philat(:))
    zlatmin = minval(zlatvals(:)-philat(:))
    lok = (max(abs(zlatmax),abs(zlatmin)) < zepslat)
    IF (.not. lok) then
      IF (((zlatvals(2)-zlatvals(1)) > 0.0_dp) .AND. ((philat(2)-philat(1)) < 0.0_dp)) THEN 
        CALL message('ef_check_geometry', 'Input latitudes are ascending while model latitudes are descending!',level=em_error)
      ELSE
        IF (((zlatvals(2)-zlatvals(1)) < 0.0_dp) .AND. ((philat(2)-philat(1)) > 0.0_dp)) THEN 
          CALL message('ef_check_geometry', 'Input latitudes are descending while model latitudes are ascending!',level=em_error)
        ELSE
          CALL message('ef_check_geometry', 'Input latitudes do not match model latitudes!',level=em_error)
        ENDIF
      ENDIF
    ELSE
      CALL nf_check (nf_inq_varid(io_file_id, 'lon',io_var_id ),TRIM(efield%ef_file))
      CALL nf_check (nf_get_vara_double(io_file_id,io_var_id,(/1/),(/ nlon /),zlonvals),TRIM(efield%ef_file))
      zepslon = abs((zlonvals(2)-zlonvals(1)) * 0.05)
      zlonmax = maxval(zlonvals(:)-philon(:))
      zlonmin = minval(zlonvals(:)-philon(:))
      lok = (max(abs(zlonmax),abs(zlonmin)) < zepslon)
      IF (.not. lok) then
        IF (((zlonvals(2)-zlonvals(1)) > 0.0_dp) .AND. ((philon(2)-philon(1)) < 0.0_dp)) THEN 
          CALL message('ef_check_geometry', 'Input longitudes are ascending while model longitudes are descending!',level=em_error)
        ELSE
          IF (((zlonvals(2)-zlonvals(1)) < 0.0_dp) .AND. ((philon(2)-philon(1)) > 0.0_dp)) THEN 
            CALL message('ef_check_geometry', 'Input longitudes are descending while model longitudes are ascending!', &
                         level=em_error)
          ELSE
            CALL message('ef_check_geometry', 'Input longitudes do not match model longitudes!',level=em_error)
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  
  END FUNCTION ef_check_geometry

!>
!! check whether the variable's dimensions of an external opened
!! NetCDF file match the given dimensions
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  FUNCTION ef_check_vardims(ncid,varid, ndim, dims) RESULT (found)
  USE mo_netcdf

  INTEGER, INTENT(IN)                   :: ncid,varid, ndim
  CHARACTER*(*), INTENT(IN)             :: dims(*)

  INTEGER                               :: i, j, vardim
  LOGICAL                               :: found, found2
  INTEGER, DIMENSION(NF_MAX_VAR_DIMS)   :: vardimids
  CHARACTER*10                          :: vardimname

  found = .FALSE.
  CALL nf_check(nf_inq_varndims(ncid, varid, vardim))
  IF (vardim .eq. ndim) THEN
    found = .TRUE.
    CALL nf_check(nf_inq_vardimid(ncid, varid, vardimids))
    DO i = 1, ndim
      CALL nf_check(nf_inq_dimname(ncid, vardimids(i), vardimname))
      found2=.FALSE.
      DO j = 1, ndim
        IF (dims(j) == vardimname) found2 = .TRUE.
      ENDDO
      IF ((.NOT. found2) .AND. (vardimname == 'lev')) THEN
        vardimname ='mlev'
        DO j = 1, ndim
          IF (dims(j) == vardimname) found2 = .TRUE.
        ENDDO
      ENDIF
      found = found .AND. found2
    ENDDO
  ENDIF
  END FUNCTION ef_check_vardims

!>
!! convert intimes (at the moment: from "days since" and month indices) to format
!! time_days used by the model run
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  SUBROUTINE ef_convert_times(dimlen, intimes, units, timedef, seasonal,     &
                              outtimes, lmonthly, newyear)
  USE mo_exception,       ONLY: message, em_error
  USE mo_time_conversion, ONLY: TC_set, TC_get, TC_convert, add_date, &
                                time_native
  USE mo_time_control,    ONLY: current_date, write_date, get_year_len

  INTEGER, INTENT(in)            :: dimlen
  REAL(dp), INTENT(inout)        :: intimes(:)
  CHARACTER (len = *), INTENT(in):: units
  INTEGER, INTENT(in)            :: timedef
  LOGICAL, INTENT(in)            :: seasonal
  TYPE(time_days), INTENT(inout) :: outtimes(:)
  LOGICAL, INTENT(out)           :: lmonthly
  INTEGER, INTENT(in)            :: newyear

  INTEGER           :: nstart, i
  INTEGER           :: year, month, day, hour, minute, second
  INTEGER           :: nyear,nmonth,nday,nhour,nminute,nsecond
  INTEGER           :: cyear,cmonth,cday,chour,cminute,csecond
  INTEGER           :: add_day, add_second
  CHARACTER(LEN=128):: clunits       ! local copy with minimum length
  LOGICAL           :: lsince
  TYPE(time_native) :: dummy_date

  ! 1) --- identify time units 
  second = 0
  lsince   = .FALSE.
  lmonthly = .FALSE.
  clunits = TRIM(units)//'                             '
  IF (clunits(1:9) == "day since") THEN
    lsince = .TRUE.
    nstart = 10 
    CALL nf_getreftime(units(nstart:), year, month, day, hour, minute)
  ELSE IF (clunits(1:10) == "days since") THEN
    lsince = .TRUE.
    nstart = 11
    CALL nf_getreftime(units(nstart:), year, month, day, hour, minute)
  ELSE IF (clunits(1:11) == "month since") THEN
    lsince = .TRUE.
    nstart = 12
    CALL nf_getreftime(units(nstart:), year, month, day, hour, minute)
    lmonthly = .TRUE.
  ELSE IF (clunits(1:12) == "months since") THEN
    lsince = .TRUE.
    nstart = 13
    CALL nf_getreftime(units(nstart:), year, month, day, hour, minute)
    lmonthly = .TRUE.
  ELSEIF (clunits(1:5) == "month") THEN
    lmonthly = .TRUE.
  ELSE
    CALL message('ef_convert_times', 'time format not supported!',level=em_error)
  ENDIF
    
  ! 2) --- convert date to standard format
  IF (lsince) THEN
    DO i = 1, dimlen
!     IF ((timedef == EF_IGNOREYEAR) .AND. (.NOT. seasonal)) THEN
!     IF (timedef == EF_IGNOREYEAR) &
!       year = newyear
!     ENDIF
      CALL TC_set(year,month,day,hour,minute,second,dummy_date)
!++mgs
      IF (lmonthly) THEN
!### PRELIMINARY: wrong results, just to make things work!!!!
!### WARNING: this means that ef_read cannot find the correct time step!!!
!sschr  add_day=30*int(intimes(i))
!! write(0,*) '### lmonthly=true: adding ',add_day,' days...'
!sschr  add_second=nint((30._dp*intimes(i)-add_day)*86400.0_dp)
!sschr  CALL add_date(add_day,add_second,dummy_date)
        nyear = year + (month+int(intimes(i))) / 12
        nmonth = mod(month+int(intimes(i)),12)
        if (nmonth == 0) nmonth = 12
        CALL TC_set(nyear,nmonth,day,hour,minute,second,dummy_date)
      ELSE
        add_day=int(intimes(i))
        add_second=nint((intimes(i)-add_day)*86400.0_dp)
        CALL add_date(add_day,add_second,dummy_date)
      END IF
!--mgs
!     IF ((timedef == EF_IGNOREYEAR) .AND. (seasonal)) THEN
      IF (timedef == EF_IGNOREYEAR) THEN
        CALL TC_get(dummy_date,cyear,cmonth,cday,chour,cminute,csecond)
        IF (get_year_len(cyear) == get_year_len(newyear)) THEN
          cyear = newyear
          CALL TC_set(cyear,cmonth,cday,chour,cminute,csecond,dummy_date)
        ELSE
          cyear = newyear
          IF ((cmonth == 2) .AND. (cday == 29)) THEN
!           quick fix: set date to February 28, 23:59:59
!           (because the mean between the previous and the next input timestep
!           might also be on February 29th, if for example hourly files are used)
            CALL TC_set(cyear,2,28,23,59,59,dummy_date)
          ELSE
            CALL TC_set(cyear,cmonth,cday,chour,cminute,csecond,dummy_date)
          ENDIF
        ENDIF
      ENDIF
      CALL TC_convert(dummy_date, outtimes(i))
    END DO
  ELSE
    DO i = 1, dimlen
      CALL TC_convert(current_date,dummy_date)
      CALL TC_get(dummy_date,year,month,day,hour,minute,second)
      IF (timedef == EF_IGNOREYEAR) year = newyear
      nmonth=int(intimes(i))
      CALL TC_set(year,nmonth,1,0,0,0,dummy_date)
      CALL TC_convert(dummy_date, outtimes(i))
    END DO
  END IF

  END SUBROUTINE ef_convert_times

!++mgs: extract reference date from "days|months since ..." string
! only the string part containing the date must be passed!

  SUBROUTINE nf_getreftime(cdate, year, month, day, hour, minute)

  CHARACTER(len=*), INTENT(in)   :: cdate
  INTEGER, INTENT(out)           :: year, month, day, hour, minute

  INTEGER   :: nlen, nstart, ndx

  nlen = LEN_TRIM(cdate)

  nstart = 1
  ndx = INDEX(cdate(nstart:nlen),'-')
  read(cdate(nstart:nstart+ndx-2),*) year
  nstart = nstart+ndx
  ndx = INDEX(cdate(nstart:nlen),'-') 
  read(cdate(nstart:nstart+ndx-2),*) month
  nstart = nstart+ndx
  ndx = INDEX(cdate(nstart:nlen),' ') 
  IF (ndx > 0) THEN
    read(cdate(nstart:nstart+ndx-2),*) day
    nstart = nstart+ndx
    ndx = INDEX(cdate(nstart:nlen),':') 
    read(cdate(nstart:nstart+ndx-2),*) hour
    nstart = nstart+ndx
    read(cdate(nstart:nlen),'(i2)') minute
  ELSE
    read(cdate(nstart:nlen),*) day
    hour = 0
    minute = 0
  ENDIF
  END SUBROUTINE nf_getreftime
!--mgs


!>
!! convert intimes to format time_days used by the model run
!! (using udunits: all formats of intimes can be processed)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
! SUBROUTINE ef_convert_times_withUdUnits(dimlen, intimes, units, timedef, outtimes)
! USE udunits
! USE mo_time_conversion, ONLY: TC_set, TC_get, TC_convert, &
!                               time_native
! USE mo_time_control,    ONLY: current_date

! INTEGER, INTENT(in)            :: dimlen
! REAL(dp), INTENT(inout)        :: intimes(:)
! CHARACTER (len = *), INTENT(in):: units
! INTEGER, INTENT(in)            :: timedef
! TYPE(time_days), INTENT(inout) :: outtimes(:)

! INTEGER           :: i, status
! INTEGER           :: year, month, day, hour, minute, isecond
! REAL              :: second
! INTEGER           :: cyear,cmonth,cday,chour,cminute,csecond
! INTEGER*8         :: timecenters_unit
! TYPE(time_native) :: dummy_date

! status = utopen("/home2/icg2/ich212/udunits-1.12.9/etc/udunits.dat")
! timecenters_unit = utmake()
! status = utdec(units, timecenters_unit)
! DO i = 1, dimlen
!     
!   status = utcaltime(intimes(i), timecenters_unit, year, month, day, &
!                        hour, minute, second)
!   isecond = nint(second)
!   if (isecond .eq. 60) then
!     intimes(i) = intimes(i) + 5.78e-6
!     status = utcaltime(intimes(i), timecenters_unit, year, month, day, &
!                        hour, minute, second)
!     isecond = nint(second)
!   endif
!   IF (timedef == EF_IGNOREYEAR) THEN
!     CALL TC_convert(current_date,dummy_date)
!     CALL TC_get(dummy_date,cyear,cmonth,cday,chour,cminute,csecond)
!     CALL TC_set(cyear,month,day,hour,minute,isecond,dummy_date)
!   ELSE
!     CALL TC_set(year,month,day,hour,minute,isecond,dummy_date)
!   ENDIF
!   CALL TC_convert(dummy_date, outtimes(i))
! END DO
! status = utfree(timecenters_unit)
! status = utcls()
! END SUBROUTINE ef_convert_times_withUdUnits

!++mgs
!>
!! check time values from EF_FILE for plausibility
!!
!! Description?
!!
!! ... assuming the values have been converted into days since ...
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  SUBROUTINE ef_check_time(efield, time, lmonthly)

  USE mo_exception,       ONLY: message, message_text, em_error, em_warn, em_info

  TYPE(external_field), INTENT(in) :: efield
  REAL(dp),             INTENT(in) :: time(:)
  LOGICAL,              INTENT(in) :: lmonthly

  INTEGER    :: ntime, i
  REAL(dp)   :: testval
  LOGICAL    :: lerr

  ntime = SIZE(time)

  ! -- case 1: file has only one time step
  IF (ntime == 1) THEN
    IF (time(1) < 0._dp .OR. time(1) > 15000._dp) THEN
      WRITE(message_text, *) 'Suspicious time value (', time(1), '). File has only one time value.'
      CALL message('ef_check_time', message_text, level=em_warn)
    END IF
  END IF

  ! -- case 2: monthly data with monthly units
  IF (lmonthly) THEN
    IF (ntime == 1) THEN
      CALL message('ef_check_time', 'Only one monthly time value in file.', level=em_warn)
    ELSE IF (ntime /= 12) THEN
      WRITE(message_text, *) ntime, ' time values with monthly time unit in file.'
      CALL message('ef_check_time', message_text, level=em_warn)
    END IF
    DO i=1,ntime-1
      testval = time(i+1)-time(i)
      IF (testval < 0._dp .OR. ABS(testval) /= 1._dp) THEN
        WRITE(message_text, *) 'Suspicious time values at i=', i, i+1,  &
                               ' (', time(i),',', time(i+1), ').'
        CALL message('ef_check_time', message_text, level=em_warn)
      END IF
    END DO

  ! -- case 3: presumably monthly data (daily units)
  ELSE IF (ntime == 12) THEN
!   CALL message('ef_check_time', '12 time values detected. Assuming these are monthly data...', &
!                level=em_info)
    DO i=1,ntime-1
      testval = time(i+1)-time(i)
      IF (testval < 0._dp .OR. ABS(testval) < 28._dp .OR. ABS(testval) > 31._dp) THEN
        WRITE(message_text, *) 'Suspicious time interval at i=', i, i+1,  &
                               ' (', time(i),',', time(i+1), ').'
        CALL message('ef_check_time', message_text, level=em_warn)
      END IF
    END DO
    IF (time(1) < 0._dp) THEN
      WRITE(message_text, *) 'Suspicious time value <0.'
      CALL message('ef_check_time', message_text, level=em_warn)
    END IF
    IF (time(12) > 366._dp) THEN
      WRITE(message_text, *) 'Time value(12) > 366.'
      CALL message('ef_check_time', message_text, level=em_info)
    END IF
  END IF

  ! -- case 3: presumably daily/annual data ...
 

  END SUBROUTINE ef_check_time
!--mgs

  SUBROUTINE ef_get_actual_filename (efield, changed)
  USE mo_time_conversion, ONLY: time_native, TC_convert, TC_get
  USE mo_time_control,    ONLY: current_date, previous_date

  TYPE(external_field), INTENT(inout) :: efield
  LOGICAL, INTENT(out)                :: changed
  CHARACTER(LEN=512)                  :: efile

  TYPE(time_native)            :: date
  INTEGER                      :: pyear, cyear, month, day, hour, minute, second

  changed = .FALSE.
  efile = expand_template(efield%ef_template, efield%ef_varname, .TRUE.)
  IF (efile /= efield%ef_file) THEN 
    efield%ef_file = efile
    changed = .TRUE.
  ENDIF
  IF (efield%ef_timedef == EF_IGNOREYEAR) THEN
! now compare year of last step with actual step
! previous_date is initialized in init_manager of mo_time_control; therefore it IS
! actually set to a valid value at the very first step of the model
    CALL TC_convert(current_date,date)
    CALL TC_get(date,cyear, month, day, hour, minute, second)
    CALL TC_convert(previous_date,date)
    CALL TC_get(date,pyear, month, day, hour, minute, second)
! ef_times have to be recalculated (and adjusted to the current year)!
    IF (cyear /= pyear) changed = .TRUE.
  ENDIF
  END SUBROUTINE ef_get_actual_filename

  SUBROUTINE ef_get_filename (efield, changed, ioffset, newyear)

! ioffset = -1 : get previous filename
! ioffset =  0 : get actual   filename
! ioffset =  1 : get next     filename

  USE mo_time_conversion, ONLY: time_native, TC_convert, TC_get, add_date
  USE mo_time_control,    ONLY: current_date, &
                                get_month_len, get_year_len
  USE mo_control,         ONLY: nn, nlev
  USE mo_filename,        ONLY: str_filter

  TYPE(external_field), INTENT(inout) :: efield
  LOGICAL, INTENT(out)                :: changed
  INTEGER, INTENT(in)                 :: ioffset
  INTEGER, INTENT(out)                :: newyear

  CHARACTER(LEN=512)                  :: efile
  CHARACTER(LEN=16)                   :: cnn, cnlev
  TYPE(time_native)                   :: date
  INTEGER                             :: year, month, day, hour, minute, second, &
                                         ndx, idays, nform

  CALL TC_convert(current_date,date)
  CALL TC_get(date,year,month,day,hour,minute,second)
  newyear = year
  changed = .FALSE.
  if (ioffset == 0) THEN
    CALL ef_get_actual_filename(efield, changed)
  ELSE
    write(cnn,'(a1,i0)') 'T',nn
    write(cnlev,'(a1,i0)') 'L',nlev
    ndx = INDEX(efield%ef_template,'%S')
    IF (ndx > 0) THEN
      CALL add_date(0,ioffset,date)
    ELSE
      ndx = INDEX(efield%ef_template,'%I')
      IF (ndx > 0) THEN
        CALL add_date(0,ioffset*60,date)
      ELSE
        ndx = INDEX(efield%ef_template,'%H')
        IF (ndx > 0) THEN
          CALL add_date(0,ioffset*3600,date)
        ELSE
          ndx = INDEX(efield%ef_template,'%D')
          IF (ndx > 0) THEN
            CALL add_date(ioffset,0,date)
          ELSE
            ndx = INDEX(efield%ef_template,'%M')
            IF (ndx > 0) THEN
              IF (ioffset == -1) THEN
                month = month - 1
                IF (month < 1) THEN
                  month = 12
                  year = year -1
                ENDIF
              ENDIF
              idays = get_month_len(year,month)
              CALL add_date(ioffset*idays,0,date)
            ELSE
              ndx = INDEX(efield%ef_template,'%Y')
              IF (ndx > 0) THEN
                IF (ioffset == -1) year = year - 1
                IF (ioffset == 1)  year = year + 1
                newyear = year
                idays = NINT(get_year_len(year))
                CALL add_date(ioffset*idays,0,date)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    CALL TC_get(date,year,month,day,hour,minute,second)
    efile = str_filter(efield%ef_template,year, month, day, hour, minute, second,nn,trim(efield%ef_varname),cnn,cnlev)
    IF (efile /= efield%ef_file) THEN 
      efield%ef_file = efile
      changed = .TRUE.
    ENDIF
  ENDIF
  END SUBROUTINE ef_get_filename

!>
!! get timestep in input file corresponding to actual model time
!!
!! Description?                 
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  SUBROUTINE ef_get_next_timestep(efield)
  USE mo_netcdf
  USE mo_exception,       ONLY: message, em_error, number_of_errors
  USE mo_time_control,    ONLY: stop_date, current_date
  USE mo_time_conversion, ONLY: add_date, time_native, TC_convert, TC_get

  TYPE(external_field), INTENT(inout) :: efield

  LOGICAL                      :: lchanged, lexist, lmonthly, lastindex
  INTEGER                      :: newyear
  INTEGER                      :: year, month, day, hour, minute, second
  TYPE(time_native)            :: cdate

  lexist = .TRUE.
  lastindex = (efield%ef_timeindex == size(efield%ef_times))
  CALL ef_get_filename (efield, lchanged, 0, newyear)
  IF (efield%ef_timedef /= EF_CONSTANT) THEN
    efield%ef_timeindex = 0
    CALL ef_get_filetimes(efield,.TRUE.,lexist, lmonthly,newyear)
    efield%ef_timeindex = ef_get_next_timeindex(efield%ef_timeindex, efield%ef_times, efield%ef_timetol, lmonthly, &
                                                efield%ef_varname)
    IF (efield%ef_timeindex == 0) THEN
! look for matching date in next file
      CALL ef_get_filename (efield, lchanged, 1, newyear)
      IF (lastindex .AND. (.NOT. lchanged) .AND. (efield%ef_timedef == EF_IGNOREYEAR)) THEN
        CALL TC_convert(current_date,cdate)
        CALL TC_get(cdate,year, month, day, hour, minute, second)
        newyear=year+1
      ENDIF
      CALL ef_get_filetimes(efield,.FALSE.,lexist, lmonthly, newyear)
      IF (lexist) THEN
        efield%ef_timeindex = ef_get_next_timeindex(0, efield%ef_times, efield%ef_timetol, lmonthly, &
                                                    efield%ef_varname)
      ELSE
! keep data until the end of the model run
        efield%ef_timeindex = 1
        ALLOCATE(efield%ef_times(1))
        efield%ef_times(1) = stop_date
        CALL add_date(1, 0, efield%ef_times(1))
      ENDIF
    ENDIF
  ENDIF
  IF (efield%ef_timeindex == 0) THEN
    IF (.NOT. lexist) THEN
      CALL message('ef_get_next_timestep','Data file '//TRIM(efield%ef_file)//' does not exist!', level=em_error)
    ELSE
      CALL message('ef_get_next_timestep','No data found in file '//TRIM(efield%ef_file)//'!', level=em_error)
    ENDIF
  ENDIF
  END SUBROUTINE ef_get_next_timestep

  SUBROUTINE ef_get_first_timestep(efield)
  USE mo_netcdf
  USE mo_exception,       ONLY: message, em_error, number_of_errors
  USE mo_time_control,    ONLY: stop_date
  USE mo_time_conversion, ONLY: add_date
  USE mo_submodel,        ONLY: emi_scenario

  TYPE(external_field), INTENT(inout) :: efield

  LOGICAL            :: lchanged, lexist, lmonthly
  INTEGER            :: newyear, ncid, varid, ndx_r, str_len
  CHARACTER(len=256) :: str2

  lexist=.FALSE.

! after the next command a possibly existing scenario string (%R) is
! permanently replaced in efield%ef_template!
  ndx_r = INDEX(efield%ef_template,'%R')
  IF (ndx_r > 0) THEN
    str_len  = LEN_TRIM(TRIM(efield%ef_template))
    str2 = efield%ef_template(1:ndx_r-1) // TRIM(emi_scenario)// efield%ef_template(ndx_r+3:str_len)
    efield%ef_template = TRIM(str2)
  ENDIF
  IF (efield%ef_timedef == EF_CONSTANT) THEN
! check whether data file exists
    CALL ef_get_filename (efield, lchanged, 0, newyear)
    CALL ef_open (efield, ncid, .FALSE.)
    IF (ncid <= 0) THEN
      CALL message('ef_get_first_timestep','Data file '//TRIM(efield%ef_file)//' does not exist!', level=em_error)
      efield%ef_timeindex = 0
    ELSE
! read unit of field once (lfirst!)
      CALL nf_check(nf_inq_varid( ncid, efield%ef_varname, varid ),efield%ef_file)
      efield%ef_actual_unit=''
      CALL nf_check(nf_get_att_text(ncid, varid, 'units',efield%ef_actual_unit),efield%ef_file)
! keep data until the end of the model run
! check whether data step exists IS STILL MISSING!
      CALL nf_check( nf_close( ncid ) )
      ALLOCATE(efield%ef_times(efield%ef_timeindex))
      efield%ef_times(efield%ef_timeindex) = stop_date
      CALL add_date(1, 0, efield%ef_times(efield%ef_timeindex))
      lexist=.TRUE.
    ENDIF
  ELSE
    efield%ef_timeindex = 0
    CALL ef_get_filename (efield, lchanged, 0, newyear)
    CALL ef_get_filetimes(efield,.TRUE.,lexist, lmonthly, newyear)
    IF (lexist) efield%ef_timeindex = ef_get_first_timeindex(efield%ef_times, efield%ef_timetol, lmonthly, &
                                                 efield%ef_varname)
    IF (efield%ef_timeindex == 0) THEN
! look for matching date in "previous" file
      CALL ef_get_filename (efield, lchanged, -1, newyear)
      IF (.NOT. lchanged) newyear = newyear - 1
      CALL ef_get_filetimes(efield,.TRUE.,lexist, lmonthly, newyear)
      IF (lexist) efield%ef_timeindex = ef_get_first_timeindex(efield%ef_times, efield%ef_timetol, lmonthly, &
                                                   efield%ef_varname)
    ENDIF
! read unit of field once (lfirst!)
    CALL ef_open (efield, ncid, .FALSE.)
    CALL nf_check(nf_inq_varid( ncid, efield%ef_varname, varid ),efield%ef_file)
    efield%ef_actual_unit=''
    CALL nf_check(nf_get_att_text(ncid, varid, 'units',efield%ef_actual_unit),efield%ef_file)
    CALL nf_check( nf_close( ncid ) )
  ENDIF
  IF (efield%ef_timeindex == 0) THEN
    IF (.NOT. lexist) THEN
      CALL message('ef_get_first_timestep','Data file '//TRIM(efield%ef_file)//' does not exist!', level=em_error)
    ELSE
      CALL message('ef_get_first_timestep','No data found in file '//TRIM(efield%ef_file)//'!', level=em_error)
    ENDIF
  ENDIF
  END SUBROUTINE ef_get_first_timestep

!>
!!
!!
!>
  SUBROUTINE ef_get_filetimes (efield, labort, lexist, lmonthly, newyear)
  USE mo_netcdf
  USE mo_exception,       ONLY: message, em_error, em_warn, number_of_errors

  TYPE(external_field), INTENT(inout) :: efield
  LOGICAL, INTENT(in)                 :: labort
  LOGICAL, INTENT(out)                :: lexist
  LOGICAL, INTENT(out)                :: lmonthly
  INTEGER, INTENT(in)                 :: newyear

  INTEGER               :: ncid, dimlen, dimid, varid
  REAL(dp), ALLOCATABLE :: filetimes(:)
  CHARACTER (len = 80)  :: units

  lexist = .FALSE.
  IF (ALLOCATED(efield%ef_times)) DEALLOCATE(efield%ef_times)
  CALL ef_open (efield, ncid, .TRUE.)
  IF (ncid > 0) THEN 
    lexist = .TRUE.
    CALL nf_check(nf_inq_dimid  (ncid, "time", dimid))
    CALL nf_check(nf_inq_dimlen (ncid, dimid, dimlen))
    ALLOCATE(filetimes(dimlen))
    ALLOCATE(efield%ef_times(dimlen))
    CALL nf_check(nf_inq_varid( ncid, "time", varid ),efield%ef_file)
    CALL nf_check (nf_get_var_double(ncid,varid,filetimes),efield%ef_file)
    filetimes = filetimes + efield%ef_timeoffset
    units=''
    CALL nf_check(nf_get_att_text(ncid, varid, 'units', units),efield%ef_file)
    CALL nf_check( nf_close( ncid ) )
  ELSE 
    IF (labort) THEN
      CALL message('ef_get_filetimes', "file "//trim(efield%ef_file)//" does not exist!", level=em_error)
    ELSE
      CALL message('ef_get_filetimes', "file "//trim(efield%ef_file)//" does not exist!", level=em_warn)
      CALL message('ef_get_filetimes', "old data will be used until the end of the model run!", level=em_warn)
    ENDIF
    !### mgs TEMPORARY CODE ###
    RETURN
  ENDIF

! Flag "seasonal" not yet implemented!!! (just test both cases)

! CALL ef_convert_times(dimlen,filetimes,trim(units),efield%ef_timedef,.TRUE.,efield%ef_times, lmonthly)
  CALL ef_convert_times(dimlen,filetimes,trim(units),efield%ef_timedef,.FALSE.,efield%ef_times, lmonthly,newyear)
  CALL ef_check_time(efield, filetimes, lmonthly)     !++mgs
  DEALLOCATE(filetimes)

  END SUBROUTINE ef_get_filetimes
!>
!!
!!
!>
  FUNCTION ef_get_next_timeindex(nstart, times, timetol, lmonthly, vname) RESULT (nind)
  USE mo_time_conversion, ONLY: time_days, &
                                OPERATOR(==), OPERATOR(<), OPERATOR(>), time_native, &
                                TC_convert, TC_get
  USE mo_time_control,    ONLY: previous_date, current_date, write_date, stop_date

  INTEGER, INTENT(in)            :: nstart
  TYPE(time_days), INTENT(in)    :: times(:)
  REAL(dp), INTENT(in)           :: timetol
  LOGICAL, INTENT(in)            :: lmonthly
  CHARACTER (len = *), INTENT(in):: vname

  INTEGER             :: i, nind, dimlen
  INTEGER             :: idayc, idayf, idummy
  LOGICAL             :: lfound
  CHARACTER (len = 4) :: cind

  nind = 0 
  lfound = .FALSE. 
  i = nstart + 1
  dimlen = SIZE(times)
  DO WHILE ((.NOT. lfound) .AND. (i <= dimlen))
    IF (current_date < times(i)) THEN
      nind = i
      lfound = .TRUE.
      write(cind,'(i4)') nind
      CALL write_date(times(nind),'.........next reading for '//TRIM(vname)//' will be done at index '//&
                      TRIM(cind)//' and file date: ')
    ENDIF
!++mgs : test tolerance
!sschr : only for ef_get_first_timeindex! (?)
!   CALL TC_get(current_date, day=idayc, second=idummy)
!   CALL TC_get(times(i), day=idayf, second=idummy)
!   IF (lmonthly .AND. (nind == 0) .AND. (ABS(idayc-idayf) < timetol)) THEN
!     nind = i
!     lfound = .TRUE.
!     write(cind,'(i4)') nind
!     CALL write_date(times(nind),'.........next reading for '//TRIM(vname)//'will be done at index '//&
!                     TRIM(cind)//' and file date: ')
!   ENDIF
!--mgs
    i = i + 1
  ENDDO
  END FUNCTION ef_get_next_timeindex

  FUNCTION ef_get_first_timeindex(times, timetol, lmonthly, vname) RESULT (nind)
  USE mo_time_conversion, ONLY: time_days, &
                                OPERATOR(==), OPERATOR(<), OPERATOR(>), time_native, &
                                TC_convert, TC_get
  USE mo_time_control,    ONLY: previous_date, current_date, write_date, stop_date

  TYPE(time_days), INTENT(in)    :: times(:)
  REAL(dp), INTENT(in)           :: timetol
  LOGICAL, INTENT(in)            :: lmonthly
  CHARACTER (len = *), INTENT(in):: vname

  INTEGER           :: i, nind, dimlen
  INTEGER           :: idayc, idayf, idummy
  LOGICAL           :: lfound
  CHARACTER (len=4) :: cind

  nind = 0  
  lfound = .FALSE. 
  i = 1
  dimlen = SIZE(times)
  DO WHILE ((.NOT. lfound) .AND. (i <= dimlen))
    IF ((current_date < times(i)) .OR. (current_date == times(i))) THEN 
      IF (current_date == times(i)) THEN 
        nind = i
      ELSE 
        nind = i - 1
      ENDIF
      lfound = .TRUE.
    ENDIF
!++mgs : test tolerance
    CALL TC_get(current_date, day=idayc, second=idummy)
    CALL TC_get(times(i), day=idayf, second=idummy)
    IF (lmonthly .AND. (nind == 0) .AND. (ABS(idayc-idayf) < timetol)) THEN 
      nind = i
      lfound = .TRUE.
    ENDIF
!--mgs
    i= i + 1
  ENDDO
  IF ((i > dimlen) .AND. (current_date > times(dimlen))) nind = dimlen
  IF (nind /= 0) THEN
    write(cind,'(i4)') nind
    CALL write_date(times(nind),'.........first reading for '//TRIM(vname)//' is done at index '//&
                    TRIM(cind)//' and file date: ')
  ENDIF
  END FUNCTION ef_get_first_timeindex

!>
! expand template of filename corresponding to actual model time
!! and actual processed variable
!! 
!! uses str_filter from mo_filename
!! (see there to find out about token replacements)
!! 
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  FUNCTION expand_template (template,spec,ltime) RESULT (efile)

  USE mo_control,         ONLY: nn, nlev
  USE mo_time_control,    ONLY: current_date
  USE mo_time_conversion, ONLY: TC_convert, TC_get, time_native
  USE mo_filename,        ONLY: str_filter

  CHARACTER(LEN=512)             :: efile
  CHARACTER(LEN=512), INTENT(IN) :: template
  CHARACTER(LEN=512), INTENT(IN) :: spec
  LOGICAL,            INTENT(IN) :: ltime
  
  INTEGER                   :: ndx, year, month, day, hour, minute, second, nform
  TYPE(time_native)         :: cdate
  CHARACTER(LEN=16)         :: cnn, cnlev

  ndx = INDEX(template,'%')
  IF (ndx > 0) THEN
    write(cnn,'(a1,i0)') 'T',nn
    write(cnlev,'(a1,i0)') 'L',nlev
    IF (ltime) THEN
      CALL TC_convert(current_date,cdate)
      CALL TC_get(cdate,year, month, day, hour, minute, second)
      efile = str_filter(template,year, month, day, hour, minute, second,nn,trim(spec),cnn,cnlev)
    ELSE
      efile = str_filter(template,0, 0, 0, 0, 0, 0,nn,trim(spec),cnn,cnlev)
    ENDIF
  ELSE
    efile=template
  ENDIF
  END FUNCTION expand_template

!>
!! open NetCDF file corresponding to actual model time
!! and actual processed variable, when a certain template
!! for the filename is given
!!
!! Description?
!! 
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!  
  SUBROUTINE ef_open(efield,io_file_id,lwith_check)

  USE mo_netcdf
  USE mo_exception, ONLY: message, em_error

  TYPE(external_field), INTENT(inout) :: efield
  INTEGER,              INTENT(out)   :: io_file_id
  LOGICAL, INTENT(in)                 :: lwith_check

  LOGICAL                           :: lexist
  
  IF (trim(efield%ef_template) /= '') THEN
    INQUIRE (file=efield%ef_file, exist=lexist)
    IF (lexist) THEN
      CALL nf_check( nf_open( TRIM(efield%ef_file), NF_NOWRITE, io_file_id ), TRIM(efield%ef_file))
      IF (lwith_check) CALL ef_check( efield, io_file_id )
    ELSE
      io_file_id = -1
    ENDIF
  ELSE
    CALL message('ef_open', 'no template name defined!', level=em_error)
  ENDIF
  END SUBROUTINE ef_open

  SUBROUTINE ef_get_full_and_half_level_information (efield, nzval, clevname, lreverse, lispressure, llevel)
  USE mo_netcdf
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, em_error, em_warn
  USE mo_control,       ONLY  : nlev

  TYPE(external_field), INTENT(inout) :: efield
  INTEGER, INTENT(IN)                 :: nzval
  CHARACTER(len=5),INTENT(IN)         :: clevname
  LOGICAL, INTENT(OUT)                :: lreverse, lispressure, llevel

  INTEGER               :: iret, io_file_id, io_var_id, ilev
  REAL(dp), ALLOCATABLE :: zzvalues(:), zzvaltemp(:)
  CHARACTER(len=80)     :: units
  CHARACTER(LEN=128)    :: clunits       ! local copy with minimum length

! get values of levels (altitude or pressure values)

  lreverse = .FALSE.
  llevel = .FALSE.
  ALLOCATE(zzvalues(nzval))
  CALL ef_open (efield, io_file_id, .FALSE.)
  iret = nf_inq_varid( io_file_id, TRIM(clevname), io_var_id )
  IF (iret /= 0 ) CALL nf_check( nf_inq_varid( io_file_id, "lev", io_var_id ),efield%ef_file)
  CALL nf_check (nf_get_var_double(io_file_id,io_var_id,zzvalues),efield%ef_file)
  units=' '
  CALL nf_check(nf_get_att_text(io_file_id, io_var_id, 'units', units),efield%ef_file)

! how to know that km are read in?!
! ==> units attribute for dimension lev MUST be set in file!!!

  clunits = TRIM(units)//'                             '
  IF (clunits(1:2) == "km") THEN
! convert from km to m
    zzvalues = zzvalues* 1000._dp
  ELSE
    IF (clunits(1:3) == "hPa") THEN
! convert from hPa to Pa 
      zzvalues = zzvalues* 100._dp
    ENDIF
  ENDIF

! Check order of given level values
! and resort, if needed ==> we want descending altitude levels (ascending pressure levels)

  IF ((clunits(1:2) == "km") .or. (clunits(1:2) == "m ")) THEN
    lispressure = .FALSE.
    IF (zzvalues(1) < zzvalues(2)) THEN
! Resort altitude array (report reversion of altitude array to calling routine as values array has also to be resorted!)
      lreverse = .TRUE.
      ALLOCATE(zzvaltemp(nzval))
      zzvaltemp = zzvalues
      DO ilev = 1, nzval
        zzvalues(ilev) = zzvaltemp(nzval - ilev + 1) 
      ENDDO
      DEALLOCATE (zzvaltemp)
    ENDIF
        
! calculate half levels (in the arithmetic middle of the full levels
    IF (ALLOCATED(efield%ef_zvalueh)) DEALLOCATE(efield%ef_zvalueh)
    ALLOCATE(efield%ef_zvalueh(nzval+1))
    efield%ef_zvalueh(nzval+1) = 0._dp
    DO ilev=nzval,2,-1
      efield%ef_zvalueh(ilev) = (zzvalues(ilev-1) + zzvalues(ilev))/2._dp
    ENDDO
    efield%ef_zvalueh(1) = efield%ef_zvalueh(2) + (zzvalues(1) - efield%ef_zvalueh(2))*2._dp
    DEALLOCATE(zzvalues)
  ELSE
    IF ((clunits(1:3) == "hPa") .or. (clunits(1:2) == "Pa")) THEN
      lispressure = .TRUE.
      IF (zzvalues(2) < zzvalues(1)) THEN 
! Resort altitude array (report reversion of altitude array to calling routine as values array has also to be resorted!)
        lreverse = .TRUE.
        ALLOCATE(zzvaltemp(nzval))
        zzvaltemp = zzvalues
        DO ilev = 1, nzval
          zzvalues(ilev) = zzvaltemp(nzval - ilev + 1) 
        ENDDO
        DEALLOCATE (zzvaltemp)
      ENDIF
             
! calculate half levels (in the arithmetic middle of the full levels
      IF (ALLOCATED(efield%ef_zvalueh)) DEALLOCATE(efield%ef_zvalueh)
      ALLOCATE(efield%ef_zvalueh(nzval+1))
      efield%ef_zvalueh(1) = 0._dp
      DO ilev=2,nzval
        efield%ef_zvalueh(ilev) = (zzvalues(ilev-1) + zzvalues(ilev))/2._dp
      ENDDO
      efield%ef_zvalueh(nzval+1) = efield%ef_zvalueh(nzval) + (zzvalues(nzval) - efield%ef_zvalueh(nzval))*2._dp
      DEALLOCATE(zzvalues)
    ELSE
      IF ((clunits(1:5) == "level") .or. (clunits(1:1) == "1")) THEN
        CALL message('ef_get_full_and_half_level_information', &
                     'Vertical coordinate system given in model levels. This is deprecated!',level=em_warn)
        IF (nzval /= nlev) THEN
          CALL message('ef_get_full_and_half_level_information', "number of model levels of input file don't fit!",level=em_error)
        ELSE
          llevel = .TRUE.
        ENDIF
      ELSE
        CALL message('ef_get_full_and_half_level_information', 'vertical coordinate system unknown!',level=em_error)
      ENDIF
    ENDIF
  ENDIF
  CALL nf_check(nf_close(io_file_id))
  END SUBROUTINE ef_get_full_and_half_level_information

  SUBROUTINE ef_get_levname_and_size(filename, clevname,nzval)
  USE mo_read_netcdf77, ONLY: read_diml_nf77

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  CHARACTER(LEN=5), INTENT(OUT) :: clevname
  INTEGER, INTENT(OUT)          :: nzval

  clevname = "lev"
  nzval=read_diml_nf77(filename,TRIM(clevname))
  IF (nzval < 0) THEN
    clevname = "mlev"
    nzval=read_diml_nf77(filename,TRIM(clevname))
  ENDIF
  END SUBROUTINE ef_get_levname_and_size

!>  
!! if needed: read new input data (single value version)
!!
!! Description?                 
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!  
  SUBROUTINE ef_read_single(efield, value)
  ! should be CALLed by p_io
  ! values should be broadcasted to other PEs after the CALL

  USE mo_read_netcdf77, ONLY: read_var_hs_nf77_0d
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, em_info

  TYPE(external_field), INTENT(inout) :: efield
  REAL(dp),             INTENT(out)   :: value
  INTEGER                             :: ierr

  call read_var_hs_nf77_0d(efield%ef_file, "time", efield%ef_timeindex, trim(efield%ef_varname), value, ierr)
  value=value*efield%ef_factor
  CALL message('ef_read_single', 'new single value read for variable '//TRIM(efield%ef_varname) &
               //' from file '//TRIM(efield%ef_file), level=em_info)

  END SUBROUTINE ef_read_single

!>
!! if needed: read new input data (1D version)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  SUBROUTINE ef_read_1d(efield, lispressure, llevel)
  ! should be CALLed by p_io
  ! values should be broadcasted to other PEs after the CALL

  USE mo_read_netcdf77, ONLY: read_var_hs_nf77_1d
  USE mo_control,       ONLY: ngl
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, em_info

  TYPE(external_field), INTENT(inout) :: efield
  LOGICAL, OPTIONAL,    INTENT(OUT)   :: lispressure, llevel

  INTEGER               :: nzval, ilev, ierr
  REAL(dp), ALLOCATABLE :: valuestemp(:)
  CHARACTER(len=5)      :: clevname
  LOGICAL               :: lreverse, lpres, llev

  IF (efield%ef_geometry .eq. EF_LEV) THEN
    CALL ef_get_levname_and_size(efield%ef_file,clevname,nzval)
    efield%ef_nzval=nzval
    ALLOCATE(value1d(nzval))
    call read_var_hs_nf77_1d(efield%ef_file, TRIM(clevname), "time", efield%ef_timeindex, trim(efield%ef_varname), value1d, ierr)
    CALL ef_get_full_and_half_level_information(efield,nzval,clevname,lreverse,lpres,llev)
    IF (lreverse) THEN
      ALLOCATE(valuestemp(nzval))
      valuestemp = value1d
      DO ilev = 1, nzval
        value1d(ilev) = valuestemp(nzval - ilev + 1) 
      ENDDO
      DEALLOCATE (valuestemp)
    ENDIF
    IF (PRESENT(lispressure)) lispressure = lpres
    IF (PRESENT(llevel)) llevel = llev
  ELSE
    ALLOCATE(value1d(ngl))
    call read_var_hs_nf77_1d(efield%ef_file, "lat", "time", efield%ef_timeindex, trim(efield%ef_varname), value1d, ierr)
  ENDIF
  value1d(:)=value1d(:)*efield%ef_factor
  CALL message('ef_read_1d', 'new 1d-field read for variable '//TRIM(efield%ef_varname) &
               //' from file '//TRIM(efield%ef_file), level=em_info)

  END SUBROUTINE ef_read_1d

!>
!! if needed: read new input data (2D version)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  SUBROUTINE ef_read_2d(efield, lispressure, llevel)
  ! should be CALLed by p_io
  ! values should be broadcasted to other PEs after the CALL

  USE mo_control,       ONLY: nlon, ngl
  USE mo_read_netcdf77, ONLY: read_var_hs_nf77_2d
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, em_info

  TYPE(external_field), INTENT(inout) :: efield
  LOGICAL, OPTIONAL,    INTENT(OUT)   :: lispressure, llevel

  INTEGER               :: ierr, nzval, ilev
  REAL(dp), ALLOCATABLE :: valuestemp(:,:)
  CHARACTER(len=5)      :: clevname
  LOGICAL               :: lreverse, lpres, llev
  
  IF (efield%ef_geometry .eq. EF_LONLAT) THEN
    ALLOCATE(value2d(nlon,ngl))
    call read_var_hs_nf77_2d(efield%ef_file, "lon", "lat", "time", efield%ef_timeindex, &
                             trim(efield%ef_varname), value2d, ierr)
  ELSE
    CALL ef_get_levname_and_size(efield%ef_file,clevname,nzval)
    efield%ef_nzval=nzval
    ALLOCATE(value2d(ngl,nzval))
    call read_var_hs_nf77_2d(efield%ef_file, "lat", TRIM(clevname), "time", efield%ef_timeindex, &
                             trim(efield%ef_varname), value2d, ierr)
    CALL ef_get_full_and_half_level_information(efield,nzval,clevname,lreverse,lpres,llev)
    IF (lreverse) THEN
      ALLOCATE(valuestemp(ngl,nzval))
      valuestemp = value2d
      DO ilev = 1, nzval
        value2d(:,ilev) = valuestemp(:,nzval - ilev + 1) 
      ENDDO
      DEALLOCATE (valuestemp)
    ENDIF
    IF (PRESENT(lispressure)) lispressure = lpres
    IF (PRESENT(llevel)) llevel = llev
  ENDIF
  value2d(:,:)=value2d(:,:)*efield%ef_factor
  CALL message('ef_read_2d', 'new 2d-field read for variable '//TRIM(efield%ef_varname) &
               //' from file '//TRIM(efield%ef_file), level=em_info)

  END SUBROUTINE ef_read_2d

!>
!! if needed: read new input data (3D version)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-08)
!!
  SUBROUTINE ef_read_3d(efield, lispressure, llevel)
  USE mo_control,       ONLY: ngl, nlon
  USE mo_read_netcdf77, ONLY: read_var_hs_nf77_3d
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, em_info

  TYPE(external_field), INTENT(INOUT) :: efield
  LOGICAL, OPTIONAL,    INTENT(OUT)   :: lispressure, llevel

  INTEGER               :: io_file_id, io_var_id, ndimid, nstep, nzval, iret, ilev, ierr
  REAL(dp), ALLOCATABLE :: valuestemp(:,:,:)
  CHARACTER(len=80)     :: units
  CHARACTER(len=5)      :: clevname
  CHARACTER(LEN=128)    :: clunits       ! local copy with minimum length
  LOGICAL               :: lreverse, lpres, llev

  CALL ef_get_levname_and_size(efield%ef_file,clevname,nzval)
  efield%ef_nzval=nzval
  ALLOCATE(value3d(nlon,ngl,nzval))
  call read_var_hs_nf77_3d(efield%ef_file, "lon", "lat", TRIM(clevname), "time", efield%ef_timeindex, &
                           trim(efield%ef_varname), value3d, ierr)
  value3d(:,:,:)=value3d(:,:,:)*efield%ef_factor
  CALL ef_get_full_and_half_level_information(efield,nzval,clevname,lreverse,lpres,llev)
  IF (lreverse) THEN
    ALLOCATE(valuestemp(nlon,ngl,nzval))
    valuestemp = value3d
    DO ilev = 1, nzval
      value3d(:,:,ilev) = valuestemp(:,:,nzval - ilev + 1) 
    ENDDO
    DEALLOCATE (valuestemp)
  ENDIF
  IF (PRESENT(lispressure)) lispressure = lpres
  IF (PRESENT(llevel)) llevel = llev
  CALL message('ef_read_3d', 'new 3d-field read for variable '//TRIM(efield%ef_varname) &
               //' from file '//TRIM(efield%ef_file), level=em_info)

  END SUBROUTINE ef_read_3d

END MODULE mo_external_field_processor 
