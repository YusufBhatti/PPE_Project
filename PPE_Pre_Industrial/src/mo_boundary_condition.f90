!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! handling of boundary conditions
!!
!! concepts of the routines:
!! see document: "ECHAM6 Boundary condition scheme" by M. G. Schultz, S. Schroeder, et al. - August 2009
!! see also: http://hammoz.icg.fz-juelich.de/data/BoundaryConditions
!!
!! @author S. Schroeder, FZ-Juelich
!!
!! $Id: 1423$
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
!! @par Copyright
!! 2009 by MPI-M and FZJ
!! This software is provided for non-commercial use only.
!!
MODULE mo_boundary_condition
  USE mo_kind,                     ONLY: dp
  USE mo_external_field_processor, ONLY: external_field, EF_INACTIVE, &
                                         EF_UNDEFINED, EF_TIMERESOLVED, &
                                         EF_NOINTER
  USE mo_memory_base,              ONLY: t_stream
#ifdef __PGI  
  USE mo_decomposition,            ONLY: gl_dc => global_decomposition, lc_dc => local_decomposition
#endif

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '1.0'

! flags of bc_domain

  INTEGER, PARAMETER :: BC_EVERYWHERE = 0
  INTEGER, PARAMETER :: BC_BOTTOM     = 1
  INTEGER, PARAMETER :: BC_TOP        = 2
  INTEGER, PARAMETER :: BC_LEVEL      = 3
  INTEGER, PARAMETER :: BC_ALTITUDE   = 4
  INTEGER, PARAMETER :: BC_PRESSURE   = 5

! flags of bc_mode

  INTEGER, PARAMETER :: BC_REPLACE  = 1
  INTEGER, PARAMETER :: BC_ADD      = 2
  INTEGER, PARAMETER :: BC_RELAX    = 3
  INTEGER, PARAMETER :: BC_SPECIAL  = 4

! flags of bc_vertint

  INTEGER, PARAMETER :: BC_VERTICAL_NONE                   = 0
  INTEGER, PARAMETER :: BC_VERTICAL_INTERPOLATION          = 1
  INTEGER, PARAMETER :: BC_VERTICAL_WEIGHTED_INTERPOLATION = 2
  INTEGER, PARAMETER :: BC_VERTICAL_INTEGRATION            = 3

  TYPE, PUBLIC :: bc_nml
  ! external field properties:
    ! source of the boundary condition
    INTEGER               :: ef_type = EF_INACTIVE   
                             ! possible values:
                             ! EF_INACTIVE = 0
                             ! EF_VALUE = 1 
                             ! EF_FILE = 2
                             ! EF_MODULE = 3
    ! information for type=ef_file
    ! template for filename. Example: bc_ch4.%Y4.%T0.nc
    CHARACTER(LEN=512)    :: ef_template = ''
    ! variable name in file
    CHARACTER(LEN=512)    :: ef_varname= ''
    ! geometry of the file variable
    INTEGER               :: ef_geometry = EF_UNDEFINED
                             ! possible values:
                             ! EF_UNDEFINED = 0
                             ! EF_3D        = 1
                             ! EF_LONLAT    = 2
                             ! EF_LATLEV    = 3
                             ! EF_LEV       = 4
                             ! EF_LAT       = 5
                             ! EF_SINGLE    = 6
    ! definition of the time values in the file
    INTEGER               :: ef_timedef = EF_TIMERESOLVED  
                             ! possible values:
                             ! EF_TIMERESOLVED = 1
                             ! EF_IGNOREYEAR   = 2
                             ! EF_CONSTANT     = 3
    ! offset to time values (e.g. to shift from mid-month to 1st day)
    REAL(dp)              :: ef_timeoffset = 0.0_dp  ! (time unit)
    ! fixed record number to select a specific time in the file
    INTEGER               :: ef_timeindex = 0
    ! time interpolation
    INTEGER               :: ef_interpolate = EF_NOINTER
                             ! possible values:
                             ! EF_NOINTER = 0
                             ! EF_LINEAR  = 1
                             ! EF_CUBIC   = 2
    ! scaling factor
    REAL(dp)              :: ef_factor = 1.0_dp 
    ! actual unit string in netcdf file
    CHARACTER(LEN=30)     :: ef_actual_unit = ''
    ! information for type=ef_value
    ! globally uniform boundary condition value
    REAL(dp)              :: ef_value       = 0.0_dp
  ! define application of boundary condition
    ! (vertical) domain where field should be applied
    INTEGER              :: bc_domain       = BC_EVERYWHERE
                             ! possible values:
                             ! BC_EVERYWHERE = 0
                             ! BC_BOTTOM     = 1
                             ! BC_TOP        = 2
                             ! BC_LEVEL      = 3
                             ! BC_ALTITUDE   = 4
                             ! BC_PRESSURE   = 5
    ! minimum and maximum model level where field shall be applied
    INTEGER              :: bc_minlev       = -1
    INTEGER              :: bc_maxlev       = 10000
    ! mode of application
    INTEGER              :: bc_mode         = BC_REPLACE
                             ! possible values:
                             ! BC_REPLACE  = 1
                             ! BC_ADD      = 2
                             ! BC_RELAX    = 3
                             ! BC_SPECIAL  = 4
    ! relaxation time for mode = bc_relax
    REAL(dp)             :: bc_relaxtime    = 0._dp
  END TYPE bc_nml

  PUBLIC :: bc_list_read
  PUBLIC :: bc_define
  PUBLIC :: bc_set
  PUBLIC :: bc_apply
  PUBLIC :: bc_query
  PUBLIC :: bc_modify      ! ++mgs 2010-02-25
  PUBLIC :: bc_find
  PUBLIC :: p_bcast_bc

! the following three lines are not supported via ICON style rules...

  PUBLIC :: BC_EVERYWHERE, BC_BOTTOM, BC_TOP, BC_LEVEL, BC_ALTITUDE, BC_PRESSURE                                     ! flags of bc_domain
  PUBLIC :: BC_REPLACE, BC_ADD, BC_RELAX, BC_SPECIAL                                                                 ! flags of bc_mode
  PUBLIC :: BC_VERTICAL_NONE, BC_VERTICAL_INTERPOLATION, BC_VERTICAL_WEIGHTED_INTERPOLATION, BC_VERTICAL_INTEGRATION ! flags of bc_vertint

  INTEGER, PARAMETER             :: MAXNBC = 800 ! Scarlet: for heterogeneous chemistry coupling with SALSA

  TYPE :: boundary_condition
    CHARACTER(LEN=128)   :: bc_name                        = ''             ! just for logging interests
    TYPE(external_field) :: bc_ef
    INTEGER              :: bc_domain                      = BC_EVERYWHERE
    INTEGER              :: bc_minlev                      = -1
    INTEGER              :: bc_maxlev                      = -1
    INTEGER              :: bc_mode                        = BC_REPLACE
    INTEGER              :: bc_ndim                        = -1
    INTEGER              :: bc_vertint                     = BC_VERTICAL_NONE
    REAL(dp)             :: bc_relaxtime                   = 0._dp
    REAL(dp), POINTER    :: bc_values_pointer(:,:,:)       => NULL()
    REAL(dp), POINTER    :: bc_values_pointer_prior(:,:,:) => NULL()
    REAL(dp), POINTER    :: bc_values_pointer_later(:,:,:) => NULL()
    LOGICAL              :: bc_ldefined                    = .false.
    LOGICAL              :: bc_lfirst                      = .TRUE.
  END TYPE boundary_condition

  TYPE(boundary_condition), SAVE :: bc_list(MAXNBC) 

  TYPE (t_stream), POINTER, SAVE :: bc_stream
  TYPE (t_stream), POINTER, SAVE :: bc_stream_prior
  TYPE (t_stream), POINTER, SAVE :: bc_stream_later

  REAL(dp), POINTER              :: tmp_glob_values2d(:,:)
  REAL(dp), POINTER              :: tmp_glob_values3d(:,:,:)

  INTEGER, SAVE                  :: nbc = 0

  LOGICAL, SAVE                  :: lfirst = .TRUE.

  INTERFACE bc_set
    MODULE PROCEDURE bc_set1d
    MODULE PROCEDURE bc_set2d
  END INTERFACE

  INTERFACE bc_apply
    MODULE PROCEDURE bc_apply1d
    MODULE PROCEDURE bc_apply2d
  END INTERFACE

  ! subprograms

  CONTAINS

!-----------------------------------------------------------------------
!>
!! expand newly read external field (0d = single) to boundary condition dimension needed
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE bc_expand0d(jbc,value)
  USE mo_control,                  ONLY: nlon, nlev, ngl

  INTEGER,  INTENT(in) :: jbc
  REAL(dp), INTENT(in) :: value

  SELECT CASE(bc_list(jbc)%bc_domain)
  CASE (BC_BOTTOM, BC_TOP)
    ALLOCATE(tmp_glob_values2d(nlon,ngl))
    tmp_glob_values2d(:,:) = value
  CASE (BC_LEVEL)
    ALLOCATE(tmp_glob_values3d(nlon,nlev,ngl))
    tmp_glob_values3d(:,:,:) = value
  END SELECT
  END SUBROUTINE bc_expand0d

!>
!! expand newly read external field to boundary condition dimension needed
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE bc_expand1d(jbc)
  USE mo_control,                  ONLY: nlon, nlev, ngl
  USE mo_external_field_processor, ONLY: EF_LAT, EF_LEV, value1d
  USE mo_exception,                ONLY: message, em_error

  INTEGER,              INTENT(in) :: jbc

  INTEGER  :: ilat, ilev

  SELECT CASE(bc_list(jbc)%bc_ef%ef_geometry)
  CASE (EF_LAT)
    IF (bc_list(jbc)%bc_domain == BC_LEVEL) THEN
      CALL message ('bc_expand', TRIM(bc_list(jbc)%bc_name)//' bc_domain=BC_LEVEL && ef_geometry=EF_LAT not possible!', &
                    level=em_error)
    ELSE
      ALLOCATE(tmp_glob_values2d(nlon,ngl))
      DO ilat = 1, ngl
        tmp_glob_values2d(:,ilat) = value1d(ilat)
      END DO
    ENDIF
  CASE (EF_LEV)
    ALLOCATE(tmp_glob_values3d(nlon,nlev,ngl))
    DO ilev = 1, nlev
      tmp_glob_values3d(:,ilev,:) = value1d(ilev)
    END DO
  END SELECT
  END SUBROUTINE bc_expand1d

!>
!! expand newly read external field to boundary condition dimension needed
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE bc_expand2d(jbc)
  USE mo_control,                  ONLY: nlon, nlev, ngl
  USE mo_external_field_processor, ONLY: EF_LONLAT, EF_LATLEV, value2d

  INTEGER,              INTENT(in) :: jbc

  INTEGER  :: ilev, ilon

  IF (bc_list(jbc)%bc_ndim == 2) THEN
    ALLOCATE(tmp_glob_values2d(nlon,ngl))
  ELSE
    ALLOCATE(tmp_glob_values3d(nlon,bc_list(jbc)%bc_ef%ef_nzval,ngl))
  ENDIF
  SELECT CASE(bc_list(jbc)%bc_domain)
  CASE (BC_BOTTOM)
    IF (bc_list(jbc)%bc_ndim == 2) THEN
      tmp_glob_values2d(:,:) = value2d(:,:)
    ELSE
      tmp_glob_values3d(:,1,:) = value2d(:,:)
    ENDIF
  CASE (BC_TOP)
    IF (bc_list(jbc)%bc_ndim == 2) THEN
      tmp_glob_values2d(:,:) = value2d(:,:)
    ELSE
      tmp_glob_values3d(:,nlev,:) = value2d(:,:)
    ENDIF
  CASE (BC_LEVEL)
    SELECT CASE (bc_list(jbc)%bc_ef%ef_geometry)
    CASE (EF_LONLAT)
      DO ilev = bc_list(jbc)%bc_minlev, bc_list(jbc)%bc_maxlev
        tmp_glob_values3d(:,ilev,:) = value2d(:,:)
      END DO
    CASE (EF_LATLEV)
      DO ilon = 1, nlon
        DO ilev = bc_list(jbc)%bc_minlev, bc_list(jbc)%bc_maxlev
          tmp_glob_values3d(ilon,ilev,:) = value2d(:,ilev)
        END DO
      END DO
    END SELECT
  CASE (BC_EVERYWHERE)
    IF (bc_list(jbc)%bc_ndim == 2) THEN
      tmp_glob_values2d(:,:) = value2d(:,:)
    ELSE
      SELECT CASE (bc_list(jbc)%bc_ef%ef_geometry)
      CASE (EF_LONLAT)
        DO ilev = 1, nlev
          tmp_glob_values3d(:,ilev,:) = value2d(:,:)
        END DO
      CASE (EF_LATLEV)
        DO ilon = 1, nlon
          tmp_glob_values3d(ilon,:,:) = TRANSPOSE(value2d(:,:))
        END DO
      END SELECT
    ENDIF
  END SELECT
  END SUBROUTINE bc_expand2d

!>
!! expand newly read external field to boundary condition dimension needed
!!
!! Description?
!!  
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE bc_expand3d(jbc)
  USE mo_control,                  ONLY: nlon, ngl
  USE mo_external_field_processor, ONLY: value3d

  INTEGER,              INTENT(in) :: jbc

  INTEGER :: ilev, nzval

  nzval = bc_list(jbc)%bc_ef%ef_nzval
  ALLOCATE(tmp_glob_values3d(nlon,nzval,ngl))
  DO ilev = 1,nzval
    tmp_glob_values3d(:,ilev,:) = value3d(:,:,ilev)
  END DO
  END SUBROUTINE bc_expand3d

  FUNCTION bc_timefac(jbc) RESULT(timefac)
  USE mo_time_control,             ONLY: current_date
  USE mo_time_conversion,          ONLY: time_days, TC_get

  INTEGER,  INTENT(in)  :: jbc

  INTEGER               :: npdays, ncdays, nndays, npseconds, ncseconds, nnseconds
  REAL(dp)              :: timefac, nom_diff, denom_diff
  TYPE(time_days)       :: ztime

  ztime = bc_list(jbc)%bc_ef%ef_times_prior(bc_list(jbc)%bc_ef%ef_timeindex_prior)
  CALL TC_get(ztime,npdays,npseconds)
  ztime  = bc_list(jbc)%bc_ef%ef_times(bc_list(jbc)%bc_ef%ef_timeindex)
  CALL TC_get(ztime,nndays,nnseconds)
  CALL TC_get(current_date,ncdays,ncseconds)

  denom_diff   = REAL(nndays-npdays,dp) + REAL(nnseconds-npseconds,dp)/86400._dp
  nom_diff     = REAL(ncdays-npdays,dp) + REAL(ncseconds-npseconds,dp)/86400._dp

  timefac   = nom_diff / denom_diff
  END FUNCTION bc_timefac

  SUBROUTINE bc_in_geom(jbc)
  USE mo_control,                  ONLY: nlev
  USE mo_external_field_processor, ONLY: EF_SINGLE, EF_LAT, EF_LEV, EF_LONLAT, EF_LATLEV, EF_3D, &
                                         ef_read_single, ef_read_1d, ef_read_2d, ef_read_3D, &
                                         value1d, value2d, value3d

  INTEGER, INTENT(in)    :: jbc
  REAL(dp)               :: value
  LOGICAL                :: lispressure, llevel

  SELECT CASE(bc_list(jbc)%bc_ef%ef_geometry)
  CASE (EF_SINGLE)
    CALL ef_read_single (bc_list(jbc)%bc_ef, value)
    CALL bc_expand0d(jbc,value)
  CASE (EF_LAT)
    CALL ef_read_1d (bc_list(jbc)%bc_ef)
    CALL bc_expand1d(jbc)
    DEALLOCATE(value1d)
  CASE (EF_LEV)
    CALL ef_read_1d (bc_list(jbc)%bc_ef, lispressure, llevel)
! this is the first time to know about bc_domain (BC_ALTITUDE || BC_PRESSURE)
    IF (llevel) THEN
      bc_list(jbc)%bc_domain = BC_LEVEL
      bc_list(jbc)%bc_minlev = 1
      bc_list(jbc)%bc_maxlev = nlev
    ELSE
      IF (lispressure) THEN
        bc_list(jbc)%bc_domain = BC_PRESSURE
      ELSE
        bc_list(jbc)%bc_domain = BC_ALTITUDE
      ENDIF
    ENDIF
    CALL bc_expand1d(jbc)
    DEALLOCATE(value1d)
  CASE (EF_LONLAT)
    CALL ef_read_2d (bc_list(jbc)%bc_ef)
    CALL bc_expand2d(jbc)
    DEALLOCATE(value2d)
  CASE (EF_LATLEV)
    CALL ef_read_2d (bc_list(jbc)%bc_ef, lispressure, llevel)
! this is the first time to know about bc_domain (BC_ALTITUDE || BC_PRESSURE)
    IF (llevel) THEN
      bc_list(jbc)%bc_domain = BC_LEVEL
      bc_list(jbc)%bc_minlev = 1
      bc_list(jbc)%bc_maxlev = nlev
    ELSE
      IF (lispressure) THEN
        bc_list(jbc)%bc_domain = BC_PRESSURE
      ELSE
        bc_list(jbc)%bc_domain = BC_ALTITUDE
      ENDIF
    ENDIF
    CALL bc_expand2d(jbc)
    DEALLOCATE(value2d)
  CASE (EF_3D)
    CALL ef_read_3d (bc_list(jbc)%bc_ef, lispressure, llevel)
! this is the first time to know about bc_domain (BC_ALTITUDE || BC_PRESSURE)
    IF (llevel) THEN
      bc_list(jbc)%bc_domain = BC_LEVEL
      bc_list(jbc)%bc_minlev = 1
      bc_list(jbc)%bc_maxlev = nlev
    ELSE
      IF (bc_list(jbc)%bc_domain == BC_EVERYWHERE) THEN
        bc_list(jbc)%bc_minlev = 1
        bc_list(jbc)%bc_maxlev = nlev
! if bc_domain = BC_PRESSURE/BC_ALTITUDE, bc_minlev/bc_maxlev have to be calculated
! (not yet ready)!
      ENDIF
      IF (lispressure) THEN
        bc_list(jbc)%bc_domain = BC_PRESSURE
      ELSE
        bc_list(jbc)%bc_domain = BC_ALTITUDE
      ENDIF
    ENDIF
    CALL bc_expand3d(jbc)
    DEALLOCATE(value3d)
  END SELECT
  END SUBROUTINE bc_in_geom

  SUBROUTINE bc_store_in_list(jbc)
#ifndef __PGI  
  USE mo_decomposition,            ONLY: gl_dc => global_decomposition, lc_dc => local_decomposition
#endif
  USE mo_mpi,                      ONLY: p_io, p_parallel, p_parallel_io, p_bcast, p_pe
  USE mo_transpose,                ONLY: scatter_gp
  USE mo_memory_base,              ONLY: add_stream_element

  INTEGER, INTENT(in)    :: jbc

  REAL(dp), ALLOCATABLE  :: tmp_loc_values2d(:,:)
  REAL(dp), ALLOCATABLE  :: tmp_loc_values3d(:,:,:)
  CHARACTER(len=5)       :: cbc = "bc000"

  IF (bc_list(jbc)%bc_ndim == 2) THEN

! from mo_memory_base:
!     info%ndim == 3:
!       info% dima(1) = ldc% nproma
!       info% dima(3) = ldc% ngpblks

    bc_list(jbc)%bc_ef%ef_nzval = 1
    IF (bc_list(jbc)%bc_lfirst) THEN 
      write(cbc(3:5),'(i3.3)') jbc
      CALL add_stream_element (bc_stream, cbc, &
                               bc_list(jbc)%bc_values_pointer, &
                               longname=TRIM(bc_list(jbc)%bc_name), &
                               klev=1)
      IF (bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) THEN
        CALL add_stream_element (bc_stream_prior, cbc, &
                                 bc_list(jbc)%bc_values_pointer_prior, &
                                 longname=TRIM(bc_list(jbc)%bc_name), &
                                 klev=1)
        CALL add_stream_element (bc_stream_later, cbc, &
                                 bc_list(jbc)%bc_values_pointer_later, &
                                 longname=TRIM(bc_list(jbc)%bc_name), &
                                 klev=1)
      ENDIF
      bc_list(jbc)%bc_lfirst = .FALSE.
    ENDIF
    ALLOCATE(tmp_loc_values2d(lc_dc%nproma,lc_dc%ngpblks))
    CALL scatter_gp(tmp_glob_values2d,tmp_loc_values2d,gl_dc)
    IF (bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) THEN
      bc_list(jbc)%bc_values_pointer_prior(1:lc_dc%nproma,1,1:lc_dc%ngpblks) = &
              bc_list(jbc)%bc_values_pointer(1:lc_dc%nproma,1,1:lc_dc%ngpblks)
    ENDIF
    bc_list(jbc)%bc_values_pointer(1:lc_dc%nproma,1,1:lc_dc%ngpblks) = tmp_loc_values2d
    IF (bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) THEN
      IF (p_parallel_io) THEN
        DEALLOCATE(tmp_glob_values2d)
        CALL bc_in_geom(jbc)
      ENDIF
      CALL scatter_gp(tmp_glob_values2d,tmp_loc_values2d,gl_dc)
      bc_list(jbc)%bc_values_pointer_later(1:lc_dc%nproma,1,1:lc_dc%ngpblks) = tmp_loc_values2d
    ENDIF
    DEALLOCATE(tmp_loc_values2d)
    IF (p_parallel_io) DEALLOCATE(tmp_glob_values2d)
  ELSE 
    IF (p_parallel) call p_bcast(bc_list(jbc)%bc_ef%ef_nzval, p_io)
    IF (bc_list(jbc)%bc_lfirst) THEN
      write(cbc(3:5),'(i3.3)') jbc
      CALL add_stream_element (bc_stream, cbc, &
                               bc_list(jbc)%bc_values_pointer, &
                               longname=TRIM(bc_list(jbc)%bc_name), &
                               klev=bc_list(jbc)%bc_ef%ef_nzval)
      IF (bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) THEN
        CALL add_stream_element (bc_stream_prior, cbc, &
                                 bc_list(jbc)%bc_values_pointer_prior, &
                                 longname=TRIM(bc_list(jbc)%bc_name), &
                                 klev=bc_list(jbc)%bc_ef%ef_nzval)
        CALL add_stream_element (bc_stream_later, cbc, &
                                 bc_list(jbc)%bc_values_pointer_later, &
                                 longname=TRIM(bc_list(jbc)%bc_name), &
                                 klev=bc_list(jbc)%bc_ef%ef_nzval)
      ENDIF
      bc_list(jbc  )%bc_lfirst = .FALSE.   
    ENDIF
    ALLOCATE(tmp_loc_values3d(lc_dc%nproma,bc_list(jbc)%bc_ef%ef_nzval,lc_dc%ngpblks))
    CALL scatter_gp(tmp_glob_values3d,tmp_loc_values3d,gl_dc)
    IF (bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) THEN
      bc_list(jbc)%bc_values_pointer_prior(1:lc_dc%nproma,1:bc_list(jbc)%bc_ef%ef_nzval,1:lc_dc%ngpblks) = &
         bc_list(jbc)%bc_values_pointer(1:lc_dc%nproma,1:bc_list(jbc)%bc_ef%ef_nzval,1:lc_dc%ngpblks)
    ENDIF
    bc_list(jbc)%bc_values_pointer(1:lc_dc%nproma,1:bc_list(jbc)%bc_ef%ef_nzval,1:lc_dc%ngpblks) = &
        tmp_loc_values3d
    IF (bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) THEN
      IF (p_parallel_io) THEN
        DEALLOCATE(tmp_glob_values3d)
        CALL bc_in_geom(jbc)
      ENDIF
      CALL scatter_gp(tmp_glob_values3d,tmp_loc_values3d,gl_dc)
      bc_list(jbc)%bc_values_pointer_later(1:lc_dc%nproma,1:bc_list(jbc)%bc_ef%ef_nzval,1:lc_dc%ngpblks) = &
         tmp_loc_values3d
    ENDIF
    DEALLOCATE(tmp_loc_values3d)
    IF (p_parallel_io) DEALLOCATE(tmp_glob_values3d)
    IF (p_parallel) THEN
      IF (((bc_list(jbc)%bc_domain == BC_ALTITUDE) .or. (bc_list(jbc)%bc_domain == BC_PRESSURE)) &
        .AND. (p_pe /= p_io)) THEN
        IF (ALLOCATED(bc_list(jbc)%bc_ef%ef_zvalueh)) DEALLOCATE(bc_list(jbc)%bc_ef%ef_zvalueh)
        ALLOCATE(bc_list(jbc)%bc_ef%ef_zvalueh(bc_list(jbc)%bc_ef%ef_nzval+1))
      ENDIF
      IF ((bc_list(jbc)%bc_domain == BC_ALTITUDE) .or. (bc_list(jbc)%bc_domain == BC_PRESSURE)) &
        CALL p_bcast(bc_list(jbc)%bc_ef%ef_zvalueh,p_io)
    ENDIF
  ENDIF
  END SUBROUTINE bc_store_in_list

!>
!! update all boundary conditions (if needed: read new ones from external fields)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE bc_list_read 
  USE mo_mpi,                      ONLY: p_io, p_parallel, p_parallel_io, p_bcast
  USE mo_time_conversion,          ONLY: p_bcast_time_days
  USE mo_external_field_processor, ONLY: EF_FILE, EF_VALUE, EF_CONSTANT, ef_get_next_timestep, ef_get_first_timestep
  USE mo_time_control,             ONLY: current_date
  USE mo_time_conversion,          ONLY: time_days, OPERATOR(==), OPERATOR(>)
  USE mo_exception,                ONLY: finish, number_of_errors
  USE mo_memory_base,              ONLY: new_stream, add_stream_element
  USE mo_control,                  ONLY: nlev

  INTEGER                :: jbc, ntimes, ntimes_prior, i
  TYPE(time_days)        :: znexttime
  LOGICAL                :: lread
  CHARACTER(len=5)       :: cbc = "bc000"

  IF (lfirst) THEN
    CALL new_stream (bc_stream,'bc',lpost=.FALSE.,lrerun=.FALSE.)
    CALL new_stream (bc_stream_prior,'bcp',lpost=.FALSE.,lrerun=.FALSE.)
    CALL new_stream (bc_stream_later,'bcl',lpost=.FALSE.,lrerun=.FALSE.)

! define individual bcs AFTER first reading
! ==> if input file consists of different number of altitude levels
! than model levels, the user (and therefore the program) might not
! know about this number at this point

! bcs can't be defined by ef_varname, as some bcs actually don't have one
! (f. ex.: "DUST emissions of DU", "SEASALT emissions of SS")
! or ef_varname isn't unique (f. ex.: "grassfire", "forestfire")
! ==> use longname attribute to identify bcs
! efield%ef_nlevels must be introduced by user!
! bc_list must know first about number of levels (in order to allocate memory) BEFORE reading

  ENDIF
  DO jbc = 1, nbc
    ! cycle if ef_type is not ef_file (but allocate memory in stream!)
    IF (bc_list(jbc)%bc_ef%ef_type /= EF_FILE) THEN
      IF (lfirst .AND. (bc_list(jbc)%bc_ef%ef_type /= EF_VALUE)) THEN
        write(cbc(3:5),'(i3.3)') jbc
        IF (bc_list(jbc)%bc_ndim == 2) THEN ! time interpolation is only done for type EF_FILE (ef_interpolate)
          CALL add_stream_element (bc_stream, cbc, &
                                   bc_list(jbc)%bc_values_pointer, &
                                   longname=TRIM(bc_list(jbc)%bc_name), &
                                   klev=1)
        ELSE
          CALL add_stream_element (bc_stream, cbc, &
                                   bc_list(jbc)%bc_values_pointer, &
                                   longname=TRIM(bc_list(jbc)%bc_name))
          bc_list(jbc)%bc_ef%ef_nzval = nlev
        ENDIF
      ENDIF
      CYCLE
    ENDIF

    lread = .FALSE.
    IF (p_parallel_io) THEN

! store to local temporary arrays -- data is finally stored to values2d/values3d in bc_expand

      IF (lfirst) THEN
! in the special case that the data for the first model timestep is in the "previous" file,
! the routine must know about lfirst
! handling of first and next timestep is slightly different -- call special routine for first
! timestep (only once called per bc), so ef_get_next_timestep (called frequently during model run)
! hasn't to deal with lfirst
        CALL ef_get_first_timestep(bc_list(jbc)%bc_ef)
      ENDIF
      IF (bc_list(jbc)%bc_ef%ef_timeindex /= 0) THEN
        znexttime = bc_list(jbc)%bc_ef%ef_times(bc_list(jbc)%bc_ef%ef_timeindex)

        IF (lfirst .OR. ((current_date > znexttime) .OR. (current_date == znexttime))) THEN
          lread = .TRUE.
          CALL bc_in_geom(jbc)
          IF (bc_list(jbc)%bc_ef%ef_timedef /= EF_CONSTANT) THEN
            IF (bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) THEN
              bc_list(jbc)%bc_ef%ef_timeindex_prior = bc_list(jbc)%bc_ef%ef_timeindex
              ntimes_prior = size(bc_list(jbc)%bc_ef%ef_times)
              IF (ALLOCATED(bc_list(jbc)%bc_ef%ef_times_prior)) DEALLOCATE(bc_list(jbc)%bc_ef%ef_times_prior)
              ALLOCATE(bc_list(jbc)%bc_ef%ef_times_prior(ntimes_prior))
              bc_list(jbc)%bc_ef%ef_times_prior = bc_list(jbc)%bc_ef%ef_times
            ENDIF
            CALL ef_get_next_timestep(bc_list(jbc)%bc_ef)
            ntimes = size(bc_list(jbc)%bc_ef%ef_times)
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    IF (p_parallel) THEN
      call p_bcast(lread, p_io)
      IF (lfirst .OR. lread) THEN
        call p_bcast(bc_list(jbc)%bc_ef%ef_actual_unit, p_io)
        call p_bcast(bc_list(jbc)%bc_vertint, p_io)
        call p_bcast(bc_list(jbc)%bc_domain, p_io)
        call p_bcast(bc_list(jbc)%bc_ef%ef_timeindex, p_io)
        call p_bcast(bc_list(jbc)%bc_ef%ef_timeindex_prior, p_io)
        ! allocation of ef_times and ef_times_prior has to be done on ALL PEs (but after
        ! reading data from file via PE_io!)
        IF (bc_list(jbc)%bc_ef%ef_timedef /= EF_CONSTANT) THEN
          IF (bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) THEN
            call p_bcast(ntimes, p_io)
            IF (.not. p_parallel_io) THEN
              IF (ALLOCATED(bc_list(jbc)%bc_ef%ef_times)) DEALLOCATE(bc_list(jbc)%bc_ef%ef_times)
              ALLOCATE(bc_list(jbc)%bc_ef%ef_times(ntimes))
            ENDIF
            DO i = 1, ntimes
              call p_bcast_time_days(bc_list(jbc)%bc_ef%ef_times(i), p_io)
            ENDDO
            call p_bcast(ntimes_prior, p_io)
            IF (.not. p_parallel_io) THEN
              IF (ALLOCATED(bc_list(jbc)%bc_ef%ef_times_prior)) DEALLOCATE(bc_list(jbc)%bc_ef%ef_times_prior)
              ALLOCATE(bc_list(jbc)%bc_ef%ef_times_prior(ntimes_prior))
            ENDIF
            DO i = 1, ntimes_prior
              call p_bcast_time_days(bc_list(jbc)%bc_ef%ef_times_prior(i), p_io)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    IF (lread) CALL bc_store_in_list(jbc)
    bc_list(jbc)%bc_ldefined = .true.  !! ignore read errors, because program will stop then anyway...
  END DO
  !! ### mgs TEMPORARY CODE ###
  IF (number_of_errors > 0) CALL finish('bc_list_read', 'Abort run because errors were encountered.')

  lfirst = .FALSE.
  END SUBROUTINE bc_list_read

  FUNCTION  bc_get_actual_unit(jbc) RESULT(strunit)
  INTEGER :: jbc

  CHARACTER(LEN=30) :: strunit

  strunit=bc_list(jbc)%bc_ef%ef_actual_unit
  END FUNCTION bc_get_actual_unit 

  SUBROUTINE bc_vert_integration(jbc,kproma,jrow,bc_values,values)
! distribute f. ex. number of molecules (number doesn't change in in vertical
! column through integration)
  USE mo_kind,      ONLY: dp
  USE mo_physical_constants, ONLY: grav
  USE mo_vphysc,    ONLY: vphysc
  USE mo_control,   ONLY: nlev

  INTEGER,               INTENT(IN)    :: jbc,kproma,jrow
  REAL(dp), ALLOCATABLE, INTENT(IN)    :: bc_values(:,:,:)
  REAL(dp), ALLOCATABLE, INTENT(INOUT) :: values(:,:)

  INTEGER  :: jk, jl, jz, jincr, jzstart, jzstop, jkstart, jkstop
  REAL(dp) :: zgi
  REAL(dp) :: zmha(kproma,nlev+1)
  REAL(dp) :: zfrac, zfdelta, zfdeltai

! initialize with 0

  values(:,:) = 0._dp

  IF (bc_list(jbc)%bc_domain == BC_ALTITUDE) THEN

! Inverse gravitational acceleration
    zgi = 1._dp/grav
! calculate height of full and half levels in m
    zmha(1:kproma,:) =  vphysc%geohm1(1:kproma,:,jrow) * zgi
    jincr   = -1
! zvalueh (= half levels in height[m] (descending order) of input) has one element more
! than zvalue (which means, it is dimensioned with ef_nzval + 1)
    jzstart = bc_list(jbc)%bc_ef%ef_nzval
    jzstop  = 1
! zmha (= half levels in height[m] (descending order) at actual geolocation) has one element more
! (which means, it is dimensioned with nlev + 1)
    jkstart = bc_list(jbc)%bc_maxlev
! as model has "no" height boundaries geohm1 has an "endless" value at index 1
! ==> therefore stop at index 2 
    jkstop  = bc_list(jbc)%bc_minlev 
    IF (jkstop == 1) jkstop=2
  ELSE
    zmha(1:kproma,:) =  vphysc%aphm1(1:kproma,:,jrow)
    jincr  = 1
    jzstart = 2
    jzstop  = bc_list(jbc)%bc_ef%ef_nzval + 1
    jkstart = bc_list(jbc)%bc_minlev + 1
    jkstop  = bc_list(jbc)%bc_maxlev
  ENDIF

! Loop over latitudes

  DO jl=1,kproma

! Loop over emission file levels        

    DO jz = jzstart, jzstop, jincr
! Calculate increment of file's half levels

      zfdelta = bc_list(jbc)%bc_ef%ef_zvalueh(jz) - bc_list(jbc)%bc_ef%ef_zvalueh(jz-jincr)
      zfdeltai = 1._dp / zfdelta

! Loop over model levels

      DO jk = jkstart, jkstop, jincr

! fraction of model level overlapping with file level 

        zfrac = MAX(0._dp, (MIN(bc_list(jbc)%bc_ef%ef_zvalueh(jz),zmha(jl,jk)) - &
                          MAX(bc_list(jbc)%bc_ef%ef_zvalueh(jz-jincr),zmha(jl,jk-jincr))) * zfdeltai)

! determine fraction of model level coinciding with file level

        zfrac = MIN(zfrac, 1._dp)
        values(jl,jk) = values(jl,jk) + zfrac * bc_values(jl,jz,jrow)

      ENDDO
    ENDDO
  ENDDO
  END SUBROUTINE bc_vert_integration


  SUBROUTINE bc_vert_weighted_interpolation(jbc,kproma,jrow,bc_values,values)
  USE mo_kind,      ONLY: dp
  USE mo_physical_constants, ONLY: grav
  USE mo_vphysc,    ONLY: vphysc
  USE mo_control,   ONLY: nlev

  INTEGER,               INTENT(IN)    :: jbc, kproma, jrow
  REAL(dp), ALLOCATABLE, INTENT(IN)    :: bc_values(:,:,:)
  REAL(dp), ALLOCATABLE, INTENT(INOUT) :: values(:,:)

  INTEGER  :: jk, jl, jz, jincr, jzstart, jzstop, jkstart, jkstop
  REAL(dp) :: zgi
  REAL(dp) :: zmha(kproma,nlev+1)
  REAL(dp) :: zfrac, zmdelta, zmdeltai

! initialize with 0

  values(:,:) = 0._dp

  IF (bc_list(jbc)%bc_domain == BC_ALTITUDE) THEN

! Inverse gravitational acceleration
    zgi = 1._dp/grav
! calculate height of full and half levels in m
    zmha(1:kproma,:) =  vphysc%geohm1(1:kproma,:,jrow) * zgi
    jincr   = -1
! zvalueh (= half levels in height[m] (descending order) of input) has one element more
! than zvalue (which means, it is dimensioned with ef_nzval + 1)
    jzstart = bc_list(jbc)%bc_ef%ef_nzval
    jzstop  = 1
! zmha (= half levels in height[m] (descending order) at actual geolocation) has one element more
! (which means, it is dimensioned with nlev + 1)
    jkstart = bc_list(jbc)%bc_maxlev
! as model has "no" height boundaries geohm1 has an "endless" value at index 1
! ==> therefore stop at index 2 
    jkstop  = bc_list(jbc)%bc_minlev 
    IF (jkstop == 1) jkstop=2
  ELSE
    zmha(1:kproma,:) =  vphysc%aphm1(1:kproma,:,jrow)
    jincr  = 1
    jzstart = 2
    jzstop  = bc_list(jbc)%bc_ef%ef_nzval + 1
    jkstart = bc_list(jbc)%bc_minlev + 1
    jkstop  = bc_list(jbc)%bc_maxlev
  ENDIF

! Loop over latitudes

  DO jl=1,kproma

! Loop over model levels

    DO jk = jkstart, jkstop, jincr

! Calculate increment of model's half levels

      zmdelta = zmha(jl,jk) - zmha(jl,jk-jincr)
      zmdeltai = 1._dp / zmdelta

! Loop over emission file levels        

      DO jz = jzstart, jzstop, jincr

! fraction of file level overlapping with model level 

        zfrac = MAX(0._dp, (MIN(zmha(jl,jk),bc_list(jbc)%bc_ef%ef_zvalueh(jz)) - &
                          MAX(zmha(jl,jk-jincr),bc_list(jbc)%bc_ef%ef_zvalueh(jz-jincr))) * zmdeltai)

! determine fraction of file level coinciding with model level

        zfrac = MIN(zfrac, 1._dp)
        values(jl,jk) = values(jl,jk) + zfrac * bc_values(jl,jz,jrow)

      ENDDO
    ENDDO
  ENDDO
  END SUBROUTINE bc_vert_weighted_interpolation

  SUBROUTINE bc_vert_interpolation(jbc,kproma,jrow,bc_values,values)
  USE mo_kind,      ONLY: dp
  USE mo_physical_constants, ONLY: grav
  USE mo_vphysc,    ONLY: vphysc
  USE mo_control,   ONLY: nlev

  INTEGER,               INTENT(IN)    :: jbc, kproma, jrow
  REAL(dp), ALLOCATABLE, INTENT(IN)    :: bc_values(:,:,:)
  REAL(dp), ALLOCATABLE, INTENT(INOUT) :: values(:,:)

  INTEGER  :: jk, jl, jz, jincr, jzstart, jzstop, jkstart, jkstop
  REAL(dp) :: zgi
  REAL(dp) :: zmha(kproma,nlev+1)
  REAL(dp) :: denom_diff, nom_diff, vertintfac

  values(:,:) = 0._dp

  IF (bc_list(jbc)%bc_domain == BC_ALTITUDE) THEN
! Inverse gravitational acceleration
    zgi = 1._dp/grav
! calculate height of full and half levels in m
    zmha(1:kproma,:) =  vphysc%geohm1(1:kproma,:,jrow) * zgi
! zvalueh (= half levels in height[m] (descending order) of input) has one element more
! than zvalue (which means, it is dimensioned with ef_nzval + 1)
    jzstart = bc_list(jbc)%bc_ef%ef_nzval
    jzstop = 1
! zmha (= half levels in height[m] (descending order) at actual geolocation) has one element more
! (which means, it is dimensioned with nlev + 1)
    jkstart = bc_list(jbc)%bc_maxlev
! as model has "no" height boundaries geohm1 has an "endless" value at index 1
! ==> therefore stop at index 2
    jkstop  = bc_list(jbc)%bc_minlev
    IF (jkstop == 1) jkstop=2

! Loop over latitudes

    DO jl=1,kproma

! Index for input levels

      jz = jzstart

! Don't extrapolate --> set upper and lower values to the values of upper and
! lower values found in input file!
! Loop over model levels

      DO jk = jkstart, jkstop, -1
  
        IF (zmha(jl,jk) .ge. bc_list(jbc)%bc_ef%ef_zvalueh(jzstart)) THEN
           values(jl,jk) = bc_values(jl,jzstart,jrow)
        ELSE IF (zmha(jl,jk) .le. bc_list(jbc)%bc_ef%ef_zvalueh(jzstop)) THEN
           values(jl,jk) = bc_values(jl,jzstop,jrow)
        ELSE
          denom_diff   = bc_list(jbc)%bc_ef%ef_zvalueh(jz) - bc_list(jbc)%bc_ef%ef_zvalueh(jz-1)
          nom_diff     = bc_list(jbc)%bc_ef%ef_zvalueh(jz) - zmha(jl,jk)

          IF (nom_diff .ge. 0.0_dp) THEN
            vertintfac   = nom_diff / denom_diff
            values(jl,jk) = bc_values(jl,jz,jrow) + &
                            vertintfac * &
                             (bc_values(jl,jz-1,jrow) &
                            - bc_values(jl,jz,jrow))
          ELSE
            jz = MAX(jz - 1, jzstop)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ELSE
    zmha(1:kproma,:) =  vphysc%aphm1(1:kproma,:,jrow)
    jzstart = 2
    jzstop  = bc_list(jbc)%bc_ef%ef_nzval + 1
    jkstart = bc_list(jbc)%bc_minlev + 1
    jkstop  = bc_list(jbc)%bc_maxlev

! Loop over latitudes

    DO jl=1,kproma

! Index for input levels

      jz = jzstart

! Don't extrapolate --> set upper and lower values to the values of upper and
! lower values found in input file!
! Loop over model levels

      DO jk = jkstart, jkstop

        IF (zmha(jl,jk) .le. bc_list(jbc)%bc_ef%ef_zvalueh(jzstart)) THEN
           values(jl,jk) = bc_values(jl,jzstart,jrow)
        ELSE IF (zmha(jl,jk) .ge. bc_list(jbc)%bc_ef%ef_zvalueh(jzstop)) THEN
           values(jl,jk) = bc_values(jl,jzstop,jrow)
        ELSE
          denom_diff   = bc_list(jbc)%bc_ef%ef_zvalueh(jz) - bc_list(jbc)%bc_ef%ef_zvalueh(jz-1)
          nom_diff     = bc_list(jbc)%bc_ef%ef_zvalueh(jz) - zmha(jl,jk)

          IF (nom_diff .ge. 0.0_dp) THEN
            vertintfac   = nom_diff / denom_diff
            values(jl,jk) = bc_values(jl,jz,jrow) + &
                            vertintfac * &
                             (bc_values(jl,jz-1,jrow) &
                            - bc_values(jl,jz,jrow))
          ELSE
            jz = MIN(jz + 1, jzstop)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  END SUBROUTINE bc_vert_interpolation
!>
!! interpolate boundary conditions vertically to the levels needed
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE p_bcast_bc (bc_struc, p_source, comm)
  USE mo_mpi, ONLY: p_bcast, p_all_comm

  TYPE(bc_nml),      INTENT(INOUT) :: bc_struc
  INTEGER,           INTENT(in)    :: p_source
  INTEGER, OPTIONAL, INTENT(in)    :: comm

  INTEGER :: p_comm

  IF (PRESENT(comm)) THEN
     p_comm = comm
  ELSE
     p_comm = p_all_comm
  ENDIF

  CALL p_bcast(bc_struc%ef_type,    p_source, p_comm)
  CALL p_bcast(bc_struc%ef_template, p_source, p_comm)
  CALL p_bcast(bc_struc%ef_varname,     p_source, p_comm)
  CALL p_bcast(bc_struc%ef_geometry,     p_source, p_comm)
  CALL p_bcast(bc_struc%ef_timedef,      p_source, p_comm)
  CALL p_bcast(bc_struc%ef_timeoffset,   p_source, p_comm)
  CALL p_bcast(bc_struc%ef_timeindex,    p_source, p_comm)
  CALL p_bcast(bc_struc%ef_interpolate,  p_source, p_comm)
  CALL p_bcast(bc_struc%ef_value,        p_source, p_comm)
  CALL p_bcast(bc_struc%ef_factor,       p_source, p_comm)
  CALL p_bcast(bc_struc%ef_actual_unit,  p_source, p_comm)
  CALL p_bcast(bc_struc%bc_domain,       p_source, p_comm)
  CALL p_bcast(bc_struc%bc_mode,         p_source, p_comm)
  CALL p_bcast(bc_struc%bc_maxlev,       p_source, p_comm)
  CALL p_bcast(bc_struc%bc_minlev,       p_source, p_comm)

  END SUBROUTINE p_bcast_bc

  SUBROUTINE bc_printstat (bc, out)
  USE mo_io_units,                 ONLY: nerr
  USE mo_external_field_processor, ONLY: EF_INACTIVE, EF_VALUE, EF_FILE, EF_MODULE, &
                                         EF_UNDEFINED, EF_3D, EF_LONLAT, EF_LAT, EF_LEV, EF_LATLEV, EF_SINGLE, &
                                         EF_TIMERESOLVED, EF_IGNOREYEAR, EF_CONSTANT, &
                                         EF_LINEAR, EF_CUBIC
  USE mo_exception,                ONLY: message, message_text, em_info, em_param

  !! print value of all bc fields ... ###
  !! private routine: called from bc_define if lverbose=.true. (debugging help)

  TYPE(boundary_condition), INTENT(IN) :: bc
  INTEGER, INTENT(in), OPTIONAL        :: out ! choice of output unit. Default: error and log file

  INTEGER           :: iout
  CHARACTER(LEN=24) :: cc, cc2

  IF (PRESENT(out)) THEN
    iout = out 
  ELSE
    iout = nerr
  END IF

  CALL message('','')
  WRITE (message_text, '(a)') 'Boundary condition:'
  CALL message('',message_text,iout, level=em_info)
  SELECT CASE (bc%bc_mode)
    CASE (BC_REPLACE)
      cc = ' BC_REPLACE'
    CASE (BC_ADD)
      cc = ' BC_ADD'
    CASE (BC_RELAX)
      cc = ' BC_RELAX'
    CASE (BC_SPECIAL)
      cc = ' BC_SPECIAL'
  END SELECT
  WRITE (message_text, '(3a,i1,2a)') 'name: ',TRIM(bc%bc_name),', ndims: ',bc%bc_ndim,', mode: ',TRIM(cc)
  CALL message('',message_text,iout,level=em_param)

  SELECT CASE (bc%bc_domain)
    CASE (BC_EVERYWHERE)
      cc = ' BC_EVERYWHERE'
    CASE (BC_BOTTOM)
      cc = ' BC_BOTTOM'
    CASE (BC_TOP)
      cc = ' BC_TOP'
    CASE (BC_LEVEL)
      cc = ' BC_LEVEL'
    CASE (BC_ALTITUDE)
      cc = ' BC_ALTITUDE'
    CASE (BC_PRESSURE)
      cc = ' BC_PRESSURE'
  END SELECT
  WRITE (message_text, '(a,L1,3a,2i4)') 'ldefined: ',bc%bc_ldefined, ', domain: ',TRIM(cc),      &
                ' minlev/maxlev : ',bc%bc_minlev,bc%bc_maxlev
  CALL message('',message_text,iout,level=em_param)
  IF (bc%bc_mode == BC_RELAX) THEN
    WRITE (message_text, '(a,E11.4)') 'relaxtime:', bc%bc_relaxtime
    CALL message('',message_text,iout,level=em_param)
  ENDIF
  SELECT CASE (bc%bc_ef%ef_type)
    CASE (EF_INACTIVE) 
      cc = ' EF_INACTIVE'
    CASE (EF_VALUE)    
      cc = ' EF_VALUE'
    CASE (EF_FILE)     
      cc = ' EF_FILE'
    CASE (EF_MODULE)   
      cc = ' EF_MODULE'
  END SELECT
  WRITE (message_text, '(2a)') 'external field type:', cc
  CALL message('',message_text,iout,level=em_param)

  IF (bc%bc_ef%ef_type == EF_FILE) THEN
    WRITE (message_text, '(2a)') 'ef_template: ', TRIM(bc%bc_ef%ef_template)
    CALL message('',message_text,iout,level=em_param)
    WRITE (message_text, '(5a,e11.4)') 'ef_varname:  ', TRIM(bc%bc_ef%ef_varname),   &
                  ', units: ', TRIM(bc%bc_ef%ef_actual_unit),            &
                  ', conversion factor : ',bc%bc_ef%ef_factor
    CALL message('',message_text,iout,level=em_param)
    SELECT CASE (bc%bc_ef%ef_geometry)
      CASE (EF_UNDEFINED)
        cc = ' EF_UNDEFINED'
      CASE (EF_3D)
        cc = ' EF_3D'
      CASE (EF_LONLAT)
        cc = ' EF_LONLAT'
      CASE (EF_LATLEV)
        cc = ' EF_LATLEV'
      CASE (EF_LEV)
        cc = ' EF_LEV'
      CASE (EF_LAT)
        cc = ' EF_LAT'
      CASE (EF_SINGLE)
        cc = ' EF_SINGLE'
    END SELECT
    WRITE (message_text, '(2a)') 'ef_geometry: ', cc
    CALL message('',message_text,iout,level=em_param)
    SELECT CASE (bc%bc_ef%ef_timedef)
      CASE (EF_TIMERESOLVED)
        cc = ' EF_TIMERESOLVED'
      CASE (EF_IGNOREYEAR)
        cc = ' EF_IGNOREYEAR'
      CASE (EF_CONSTANT)
        cc = ' EF_CONSTANT'
    END SELECT
    SELECT CASE (bc%bc_ef%ef_interpolate)
      CASE (EF_NOINTER)
        cc2 = ' EF_NOINTER'
      CASE (EF_LINEAR)
        cc2 = ' EF_LINEAR'
      CASE (EF_CUBIC)
        cc2 = ' EF_CUBIC'
    END SELECT
    WRITE (message_text, '(3a,e11.4,a,i4,2a)') 'ef_timedef: ', TRIM(cc), ', ef_timeoffset: ',bc%bc_ef%ef_timeoffset,  &
                  ', ef_timeindex: ', bc%bc_ef%ef_timeindex, ', time interpolation: ', TRIM(cc2)
    CALL message('',message_text,iout,level=em_param)
  ELSE IF (bc%bc_ef%ef_type == EF_VALUE) THEN
    WRITE (message_text, '(a,e11.4)') 'constant BC value: ', bc%bc_ef%ef_value
    CALL message('',message_text,iout,level=em_param)
  END IF
  CALL message('','',level=em_param)

  END SUBROUTINE bc_printstat

!++mgs 2010-03-01
  SUBROUTINE bc_find (name, ibc, ierr)

  USE mo_exception,           ONLY: finish, message, message_text, em_info, em_error
  USE mo_util_string,         ONLY: tolower

  CHARACTER(len=*), INTENT(in)   :: name
  INTEGER, INTENT(out)           :: ibc   ! index in bc_list
  INTEGER, OPTIONAL, INTENT(out) :: ierr  ! error status

  INTEGER     :: i

  !-- Initialize and check validity of arguments
  ibc = 0
  IF (PRESENT(ierr)) THEN
    ierr = 0
  END IF
  IF (TRIM(name) == '') CALL finish('bc_find', 'Invalid (empty) name string!')

  DO i=1, nbc
    IF (TRIM(tolower(name)) == TRIM(tolower(bc_list(i)%bc_name))) ibc = i
  END DO

  IF (ibc > 0) THEN
    WRITE(message_text,'(a,i4)') 'Located boundary condition for "'//TRIM(name)//'" as index ', ibc
    CALL message('bc_find', message_text, level=em_info)
  ELSE
    IF (PRESENT(ierr)) THEN
      ierr = 1
    ELSE
      CALL message('bc_find', 'Cannot find index for boundary condition '//TRIM(name),   &
                   level=em_error)
    END IF
  END IF

  END SUBROUTINE bc_find
!--mgs

!++mgs
  SUBROUTINE bc_query (ibc, name, ef_type, ef_actual_unit)

  USE mo_exception,         ONLY: finish

  INTEGER, INTENT(in)                      :: ibc   ! index in bc_list
  CHARACTER(len=*), INTENT(out), OPTIONAL  :: name
  INTEGER,          INTENT(out), OPTIONAL  :: ef_type
  CHARACTER(len=*), INTENT(out), OPTIONAL  :: ef_actual_unit

  !-- check validity of arguments
  IF (ibc < 0 .OR. ibc > nbc) CALL finish('bc_query', 'Invalid index to bc_list!')

  IF (PRESENT(name)) THEN
    name = bc_list(ibc)%bc_name
  END IF
  IF (PRESENT(ef_type)) THEN
    ef_type = bc_list(ibc)%bc_ef%ef_type
  END IF
  IF (PRESENT(ef_actual_unit)) THEN
    ef_actual_unit = bc_list(ibc)%bc_ef%ef_actual_unit
  END IF

  END SUBROUTINE bc_query
!--mgs

!++mgs  2010-02-25
  SUBROUTINE bc_modify (ibc, name, ef_type, bc_domain, bc_ndims, bc_vertint, ef_actual_unit)

  USE mo_exception,         ONLY: finish

  INTEGER, INTENT(in)                      :: ibc   ! index in bc_list
  CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
  INTEGER,          INTENT(in), OPTIONAL   :: ef_type
  INTEGER,          INTENT(in), OPTIONAL   :: bc_domain
  INTEGER,          INTENT(in), OPTIONAL   :: bc_ndims
  INTEGER,          INTENT(in), OPTIONAL   :: bc_vertint
  CHARACTER(len=*), INTENT(in), OPTIONAL   :: ef_actual_unit

  !-- check validity of arguments
  IF (ibc < 0 .OR. ibc > nbc) CALL finish('bc_query', 'Invalid index to bc_list!')

  IF (PRESENT(name)) THEN
    bc_list(ibc)%bc_name = TRIM(name)
  END IF
  IF (PRESENT(ef_type)) THEN
    bc_list(ibc)%bc_ef%ef_type = ef_type     ! add error checks (?)
  END IF
  IF (PRESENT(bc_domain)) THEN
    bc_list(ibc)%bc_domain = bc_domain       ! add error checks (?)
  END IF
  IF (PRESENT(bc_ndims)) THEN
    bc_list(ibc)%bc_ndim = bc_ndims         ! add error checks (?)
    ! may need to do more here to ensure consistency...
  END IF
  IF (PRESENT(bc_vertint)) THEN
    bc_list(ibc)%bc_vertint = bc_vertint
  END IF
  IF (PRESENT(ef_actual_unit)) THEN
    bc_list(ibc)%bc_ef%ef_actual_unit = ef_actual_unit
  END IF

  END SUBROUTINE bc_modify
!--mgs

!>
!! function to define a new boundary condition and return its index in the
!! boundary condition list
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  FUNCTION bc_define(bc_name,bc_struc,bc_ndim,lverbose) RESULT(ibc)
  USE mo_external_field_processor, ONLY: EF_LAT, EF_LEV, EF_LATLEV, &
                                         EF_INACTIVE, EF_VALUE, EF_FILE, EF_MODULE
  USE mo_exception,                ONLY: finish, message, message_text, em_error, em_warn
  USE mo_control,                  ONLY: nlev
  USE mo_mpi,                      ONLY: p_parallel_io

  CHARACTER(LEN=*),  INTENT(in) :: bc_name
  TYPE(bc_nml),      INTENT(in) :: bc_struc
  INTEGER, INTENT(in)           :: bc_ndim
  LOGICAL, OPTIONAL, INTENT(in) :: lverbose

  CHARACTER(LEN=24)            :: cc
  INTEGER                      :: ibc

  ibc = -1

! not to be used in this way bc_name(CHARACTER(128)), message_text(CHARACTER(132))
! WRITE (message_text, '(a)') 'Boundary condition ' // TRIM(bc_name) // ' of type'
  SELECT CASE (bc_struc%ef_type)
    CASE (EF_INACTIVE) 
      cc = ' EF_INACTIVE'
    CASE (EF_VALUE)    
      cc = ' EF_VALUE'
    CASE (EF_FILE)     
      cc = ' EF_FILE'
    CASE (EF_MODULE)   
      cc = ' EF_MODULE'
    CASE DEFAULT       
      cc = ' INVALID'
      CALL message('bc_define', 'Invalid value of EF_TYPE for boundary condition ' &
                   // TRIM(bc_name) //'!', level=em_error)
      RETURN
  END SELECT
! not to be used in this way bc_name(CHARACTER(128)), message_text(CHARACTER(132))
! CALL message('bc_define', TRIM(message_text) // cc, level=em_info)
  IF (bc_struc%ef_type == EF_INACTIVE) RETURN  ! nothing to be done

  nbc = nbc + 1
  ! Error checks
  IF (nbc > MAXNBC) THEN
    WRITE (message_text, '(a,i4,a)') 'Too many boundary conditions (MAXNBC=', MAXNBC, ')'
    CALL finish('bc_define', message_text)
  END IF
  ibc = nbc
  ! bc_mode
  IF (bc_struc%bc_mode < BC_REPLACE .OR. bc_struc%bc_mode > BC_SPECIAL) THEN
    WRITE (message_text, '(a,i2)') 'Invalid value for bc_mode : ', bc_struc%bc_mode
    CALL message('bc_define', message_text, level=em_error)
  END IF
  ! bc_domain
  IF (bc_struc%bc_domain < BC_EVERYWHERE .OR. bc_struc%bc_domain > BC_PRESSURE) THEN
    WRITE (message_text, '(a,i2)') 'Invalid value for bc_domain : ', bc_struc%bc_domain
    CALL message('bc_define', message_text, level=em_error)
  END IF
  ! bc_relaxtime
  IF (bc_struc%bc_relaxtime < 0._dp) THEN
    WRITE (message_text, '(a,e11.4)') 'Invalid value for bc_relaxtime : ', bc_struc%bc_relaxtime
    CALL message('bc_define', message_text, level=em_error)
  END IF
  IF ((bc_struc%bc_mode == BC_RELAX) .AND. (bc_struc%bc_relaxtime == 0._dp)) &
    CALL message('bc_define', 'value for bc_relaxtime must be set if bc_mode=BC_RELAX!', level=em_error)

  ! settings for all ef_type values
  bc_list(nbc)%bc_name              = bc_name
  bc_list(nbc)%bc_ef%ef_type        = bc_struc%ef_type
  bc_list(nbc)%bc_domain            = bc_struc%bc_domain
  bc_list(nbc)%bc_minlev            = MAX(bc_struc%bc_minlev, 1)  ! i.e. also 1 if "undefined"
  bc_list(nbc)%bc_maxlev            = MIN(bc_struc%bc_maxlev, nlev)
  IF (bc_list(nbc)%bc_maxlev < 1) bc_list(nbc)%bc_maxlev = nlev   ! any value < 1 could otherwise cause segfault
  bc_list(nbc)%bc_mode              = bc_struc%bc_mode  
  bc_list(nbc)%bc_relaxtime         = bc_struc%bc_relaxtime
  bc_list(nbc)%bc_ef%ef_nzval       = nlev
  ! ef_actual_unit has to be reported by ALL boundary conditions (also by modules!)
  bc_list(nbc)%bc_ef%ef_actual_unit = bc_struc%ef_actual_unit

  SELECT CASE (bc_struc%ef_type)
    CASE (EF_VALUE)    
      bc_list(nbc)%bc_ef%ef_value       = bc_struc%ef_value

! 2014/05/16 sschr: fix for bug #346
! the above statement sets the value of this boundary condition already at definition time
      bc_list(nbc)%bc_ldefined = .true.

      ! note: scale factor always = 1.0 (=default) !
      IF (bc_struc%ef_factor /= 1.0_dp) CALL message ('bc_define',  &
          ' ef_factor is set though ef_type is EF_VALUE (the factor is ignored!!)', level=em_warn)

    CASE (EF_FILE)     
      bc_list(nbc)%bc_ef%ef_template    = bc_struc%ef_template
      bc_list(nbc)%bc_ef%ef_varname     = bc_struc%ef_varname
      bc_list(nbc)%bc_ef%ef_geometry    = bc_struc%ef_geometry
      bc_list(nbc)%bc_ef%ef_timedef     = bc_struc%ef_timedef
      bc_list(nbc)%bc_ef%ef_timeoffset  = bc_struc%ef_timeoffset
      bc_list(nbc)%bc_ef%ef_timeindex   = bc_struc%ef_timeindex
      bc_list(nbc)%bc_ef%ef_interpolate = bc_struc%ef_interpolate
      bc_list(nbc)%bc_ef%ef_factor      = bc_struc%ef_factor
      IF ( ( ((bc_struc%bc_domain == BC_BOTTOM) .OR. (bc_struc%bc_domain == BC_TOP)) &
           .AND. ((bc_struc%ef_geometry == EF_LEV) .OR. (bc_struc%ef_geometry == EF_LATLEV)) ) &
         .OR. ( (bc_struc%bc_domain == BC_LEVEL) .AND. (bc_struc%ef_geometry == EF_LAT) ) ) &
        CALL message ('bc_define', ' bc_domain and ef_geometry do not fit!', level=em_error)

    CASE (EF_MODULE)   
      ! nothing to be done here. Module must call bc_set to set BC value. Reserve memory in bc_list ==> bc_ndim must be set!
  END SELECT

  bc_list(nbc)%bc_ndim = bc_ndim

  IF (PRESENT(lverbose)) THEN
    IF (lverbose .AND. p_parallel_io) CALL bc_printstat(bc_list(nbc))
  ENDIF

  END FUNCTION bc_define

!>
!! set boundary condition from ECHAM module (2D)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE bc_set1d(jbc,kproma,jrow,values)
  USE mo_external_field_processor, ONLY: EF_MODULE

  INTEGER,  INTENT(in) :: jbc, kproma, jrow
  REAL(dp), INTENT(in) :: values(:)

  ! overwrite bc value only if ef_type=EF_MODULE
  IF (bc_list(jbc)%bc_ef%ef_type /= EF_MODULE) RETURN

  ! check SHAPE of incoming array
  ! careful with nproma ...
  ! in the end do
  bc_list(jbc)%bc_values_pointer(1:kproma,1,jrow) = values(1:kproma)
  bc_list(jbc)%bc_ldefined = .true.

  END SUBROUTINE bc_set1d

!>
!! set boundary condition in ECHAM module (3D)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE bc_set2d(jbc,kproma, jrow,values)
  USE mo_external_field_processor, ONLY: EF_MODULE

  INTEGER,  INTENT(in) :: jbc, kproma, jrow
  REAL(dp), INTENT(in) :: values(:,:)

  ! overwrite bc value only if ef_type=EF_MODULE
  IF (bc_list(jbc)%bc_ef%ef_type /= EF_MODULE) RETURN

  ! check SHAPE of incoming array
  ! careful with nproma ...
  ! in the end do
  bc_list(jbc)%bc_values_pointer(1:kproma,:,jrow) = values(1:kproma,:)
  bc_list(jbc)%bc_ldefined = .true.

  END SUBROUTINE bc_set2d


!>
!! apply boundary condition to user field (1D - vertical column)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE bc_apply1d(jbc,kproma,jrow,values)
  USE mo_time_control,             ONLY: delta_time
  USE mo_exception,                ONLY: finish
  USE mo_external_field_processor, ONLY: EF_VALUE, EF_FILE
#ifndef __PGI  
  USE mo_decomposition,            ONLY: lc_dc => local_decomposition
#endif

  INTEGER,  INTENT(in)    :: jbc, kproma, jrow
  REAL(dp), INTENT(inout) :: values(:)

  REAL(dp)             :: factor
  REAL(dp),ALLOCATABLE :: bc_values(:,:,:) 

  IF ( .NOT. bc_list(jbc)%bc_ldefined ) &
    CALL finish('bc_apply', TRIM(bc_list(jbc)%bc_name)//' has no values defined yet!')

  ALLOCATE(bc_values(1:lc_dc%nproma,1,1:lc_dc%ngpblks))
  IF (bc_list(jbc)%bc_ef%ef_type == EF_VALUE) THEN
    bc_values = bc_list(jbc)%bc_ef%ef_value
  ELSE
    IF ((bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) &
         .AND. (bc_list(jbc)%bc_ef%ef_type == EF_FILE)) THEN
      bc_values = (1._dp - bc_timefac(jbc)) * bc_list(jbc)%bc_values_pointer &
                + bc_timefac(jbc) * bc_list(jbc)%bc_values_pointer_later
    ELSE
      bc_values = bc_list(jbc)%bc_values_pointer
    ENDIF
  ENDIF
!!!!  - following de-activated because BC_BOTTOM or BC_TOP should be valid conditions here
!!   IF ( .NOT. (bc_list(jbc)%bc_domain == BC_EVERYWHERE) ) &
!!     CALL finish('bc_apply', TRIM(bc_list(jbc)%bc_name)//' bc_domain indicates 3d field, but 2d field passed!')
  SELECT CASE(bc_list(jbc)%bc_mode)
  CASE (BC_REPLACE)
    values(1:kproma) = bc_values(1:kproma,1,jrow)
  CASE (BC_ADD)
    values(1:kproma) = values(1:kproma) + bc_values(1:kproma,1,jrow)
  CASE (BC_RELAX)
    factor = exp(-(delta_time/bc_list(jbc)%bc_relaxtime))
    values(1:kproma) = (1.0_dp - factor) * values(1:kproma) + factor * bc_values(1:kproma,1,jrow)
  CASE (BC_SPECIAL)
! user should pass a temporary field (not the real boundary condition array)
    values(1:kproma) = bc_values(1:kproma,1,jrow)
  END SELECT
  DEALLOCATE(bc_values)
  END SUBROUTINE bc_apply1d

!>
!! apply boundary condition to user field (2D)
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2009-08-25)
!!
  SUBROUTINE bc_apply2d(jbc,kproma,jrow,values,minlev,maxlev)
  USE mo_time_control,             ONLY: delta_time
  USE mo_exception,                ONLY: finish
  USE mo_control,                  ONLY: nlev
  USE mo_external_field_processor, ONLY: EF_VALUE, EF_FILE
#ifndef __PGI  
  USE mo_decomposition,            ONLY: lc_dc => local_decomposition
#endif

  INTEGER,  INTENT(in)    :: jbc, kproma, jrow
  REAL(dp), INTENT(inout) :: values(:,:)
  INTEGER,  OPTIONAL,  INTENT(out) :: minlev, maxlev  ! return minlev and maxlev from structure for SPECIAL mode

  INTEGER               :: k, kmin, kmax
  REAL(dp)              :: factor
  REAL(dp), ALLOCATABLE :: bc_values(:,:,:) 
  REAL(dp), ALLOCATABLE :: ztmpvalues(:,:)

  IF ( .NOT. bc_list(jbc)%bc_ldefined ) &
    CALL finish('bc_apply', TRIM(bc_list(jbc)%bc_name)//' has no values defined yet!')

  ALLOCATE(bc_values(1:lc_dc%nproma,bc_list(jbc)%bc_ef%ef_nzval,1:lc_dc%ngpblks))
  IF (bc_list(jbc)%bc_ef%ef_type == EF_VALUE) THEN
    bc_values = bc_list(jbc)%bc_ef%ef_value
  ELSE
    IF ((bc_list(jbc)%bc_ef%ef_interpolate /= EF_NOINTER) &
         .AND. (bc_list(jbc)%bc_ef%ef_type == EF_FILE)) THEN
      bc_values = (1._dp - bc_timefac(jbc)) * bc_list(jbc)%bc_values_pointer &
                + bc_timefac(jbc) * bc_list(jbc)%bc_values_pointer_later
    ELSE
      bc_values = bc_list(jbc)%bc_values_pointer
    ENDIF
  ENDIF

  SELECT CASE(bc_list(jbc)%bc_domain)
  CASE (BC_EVERYWHERE)
  ! if BC_EVERYWHERE is the domain of a field given in altitude/pressure levels,
  ! then BC_EVERYWHERE has already been converted (by bc_in_geom) to BC_ALTITUDE/BC_PRESSURE
  ! with bc_minlev/bc_maxlev=1/nlev
    kmin = 1
    kmax = nlev
  CASE (BC_BOTTOM)
    kmin = nlev
    kmax = nlev
  CASE (BC_TOP)
    kmin = 1
    kmax = 1
  CASE (BC_LEVEL)
    kmin =bc_list(jbc)%bc_minlev
    kmax =bc_list(jbc)%bc_maxlev
  CASE (BC_ALTITUDE, BC_PRESSURE)
    ALLOCATE(ztmpvalues(SIZE(values,1),SIZE(values,2)))
    SELECT CASE(bc_list(jbc)%bc_vertint)
    CASE (BC_VERTICAL_NONE)
    CASE (BC_VERTICAL_INTERPOLATION)    
      call bc_vert_interpolation(jbc,kproma,jrow,bc_values,ztmpvalues)
    CASE (BC_VERTICAL_WEIGHTED_INTERPOLATION)    
      call bc_vert_weighted_interpolation(jbc,kproma,jrow,bc_values,ztmpvalues)
    CASE (BC_VERTICAL_INTEGRATION)    
      call bc_vert_integration(jbc,kproma,jrow,bc_values,ztmpvalues)
    END SELECT
    ! not yet ready:
    ! bc_minlev/bc_maxlev have to be calculated into indices (by the routine that reads the application domain)
    ! at this point bc_minlev/bc_maxlev are always 1/nlev
    kmin =bc_list(jbc)%bc_minlev
    kmax =bc_list(jbc)%bc_maxlev
  END SELECT

  IF (bc_list(jbc)%bc_ndim == 2) THEN
! now apply a 1d field to the 2d argument field
    SELECT CASE(bc_list(jbc)%bc_mode)
    CASE (BC_REPLACE)
      DO k=kmin,kmax
        values(1:kproma,k) = bc_values(1:kproma,1,jrow)
      END DO
    CASE (BC_ADD)
      DO k=kmin,kmax
        values(1:kproma,k) = values(1:kproma,k) + bc_values(1:kproma,1,jrow)
      END DO
    CASE (BC_RELAX)
      factor = exp(-(delta_time/bc_list(jbc)%bc_relaxtime))
      DO k=kmin,kmax
        values(1:kproma,k) = (1.0_dp - factor) * values(1:kproma,k) + factor * bc_values(1:kproma,1,jrow)
      END DO
    CASE (BC_SPECIAL)
! user passed a temporary 2d field, but would get a 1d field back
      CALL finish('bc_apply', TRIM(bc_list(jbc)%bc_name)//' a 3d field was passed, but bc_list contains a 2d field!')
    END SELECT
  ELSE
! now apply a 2d field to the 2d argument field
    SELECT CASE(bc_list(jbc)%bc_mode)
    CASE (BC_REPLACE)
      IF (ALLOCATED(ztmpvalues)) THEN
        DO k=kmin,kmax
          values(1:kproma,k) = ztmpvalues(1:kproma,k)
        END DO
        DEALLOCATE(ztmpvalues)
      ELSE
        DO k=kmin,kmax
          values(1:kproma,k) = bc_values(1:kproma,k,jrow)
        END DO
      ENDIF
    CASE (BC_ADD)
      IF (ALLOCATED(ztmpvalues)) THEN 
        DO k=kmin,kmax
          values(1:kproma,k) = values(1:kproma,k) + ztmpvalues(1:kproma,k)
        END DO
        DEALLOCATE(ztmpvalues)
      ELSE
        DO k=kmin,kmax
          values(1:kproma,k) = values(1:kproma,k) + bc_values(1:kproma,k,jrow)
        END DO
      ENDIF
    CASE (BC_RELAX)
      factor = exp(-(delta_time/bc_list(jbc)%bc_relaxtime))
      IF (ALLOCATED(ztmpvalues)) THEN 
        DO k=kmin,kmax
          values(1:kproma,k) = (1.0_dp - factor) * values(1:kproma,k) + factor * ztmpvalues(1:kproma,k)
        END DO
        DEALLOCATE(ztmpvalues)
      ELSE 
        DO k=kmin,kmax
          values(1:kproma,k) = (1.0_dp - factor) * values(1:kproma,k) + factor * bc_values(1:kproma,k,jrow)
        END DO
      ENDIF  
    CASE (BC_SPECIAL)
      IF (ALLOCATED(ztmpvalues)) THEN
        DO k=kmin,kmax
! user should pass a temporary field (not the real boundary condition array)
          values(1:kproma,k) = ztmpvalues(1:kproma,k)
        END DO
        DEALLOCATE(ztmpvalues)
      ELSE
        DO k=kmin,kmax
! user should pass a temporary field (not the real boundary condition array)
          values(1:kproma,k) = bc_values(1:kproma,k,jrow)
        END DO
      ENDIF
    END SELECT
  ENDIF
  !! return kmin and kmax as minlev, maxlev if requested
  IF (PRESENT(minlev)) minlev = kmin
  IF (PRESENT(maxlev)) maxlev = kmax
  DEALLOCATE(bc_values)
  END SUBROUTINE bc_apply2d

END MODULE mo_boundary_condition
