!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_global_diag.f90
!!
!! \brief
!! Evaluate global diagnostics (after row loop in scan1) for various parts of the HAMMOZ model
!!
!! \author Martin Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (XXXX-XX-XX)
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

MODULE mo_hammoz_global_diag

  !-- include use statements for diag pointers here so that they are available to
  !   initialisation and diagnostics routine
  USE mo_kind,                ONLY: dp
  USE mo_submodel_diag,       ONLY: t_diag

  IMPLICIT NONE

  PRIVATE

  ! -- public module routines
  PUBLIC  :: init_hammoz_gdiag
  PUBLIC  :: hammoz_gdiag

  ! -- interface to provate routine
  INTERFACE define_variable
    MODULE PROCEDURE define_variable_2d
    MODULE PROCEDURE define_variable_3d
  END INTERFACE

  ! -- private module variables
  CHARACTER(len=80), SAVE   :: cfilename                   ! name of netcdf output file
  INTEGER, SAVE             :: nfvar_time                  ! variable id of time in netcdf file
  LOGICAL, SAVE             :: ldebug                      ! write diagnostics at each time step
  INTEGER, SAVE             :: navg                        ! number of time steps for averaging
                                                           ! (daily output)
  INTEGER, PARAMETER        :: maxvars = 100
  INTEGER, PARAMETER        :: GSUM = 1,              &    ! global sum
                               GAVG = 2,              &    ! global mean
                               GWAVG = 3,             &    ! global weighted average (surface area)
                               GTAU = 4                    ! global lifetime average (inverse)
  INTEGER                   :: ndiag                       ! number of diagnostics defined
  TYPE(t_diag)              :: dps(maxvars)                ! list of diagnostic pointers
                                                           ! use key_type for processing info (see GSUM etc.)
                                                           ! use key for ncdf variable id
  REAL(dp)                  :: dbuf(maxvars)               ! buffer for global results
  REAL(dp)                  :: scalefac(maxvars)           ! scaling factor (unit conversion etc.)
  REAL(dp)                  :: garea                       ! sum(gboxarea_2d)


  CONTAINS


  SUBROUTINE init_hammoz_gdiag

    USE mo_exception,           ONLY: message, message_text, em_info
    USE mo_filename,            ONLY: out_expname
    USE mo_netcdf,              ONLY: nf_check, nf_create, nf_close, nf_clobber,      &
                                      nf_double, nf_real, nf_def_dim, nf_def_var,     &
                                      nf_unlimited, nf_put_att_text, nf_global,       &
                                      nf_enddef
    USE mo_time_control,        ONLY: time_step_len, current_date, get_date_components
    USE mo_moz_diag,            ONLY: dpalb, dpohconc, dpratech4, dpratemcl
    USE mo_moz_lightning,       ONLY: dpff, dpemino
    USE mo_geoloc,              ONLY: gboxarea_2d
    USE mo_global_op,           ONLY: sum_global

    INTEGER              :: year, month, day, hour, minute, second
    INTEGER              :: nfid, nfdim_time
    CHARACTER(len=9)     :: cdate
    CHARACTER(len=27)    :: crefdate

    ndiag = 0

    ldebug = .TRUE.             ! detailed diagnostics and output every time step
    IF (ldebug) THEN
      navg = 1
      CALL message('init_hammoz_gdiag', 'DEBUGGING MODE: Output every time step!', level=em_info)
    ELSE
      navg = NINT(86400._dp/time_step_len)
      write (message_text,'(a,i0,a)') 'Averaging over ',navg,     &
                             ' time steps for daily output of global diagnostics.'
      CALL message('init_hammoz_gdiag', message_text, level=em_info)
    END IF

    !-- construct filename for netcdf file
    ! (we cannot use streams, because there are no scalar stream variables possible)
    CALL get_date_components(current_date, year = year, month = month, day = day,   &
                             hour = hour, minute = minute)
    write(cdate,'(i4.4,i2.2,".",i2.2)') year, month, day
    cfilename = TRIM(out_expname)//'_'//cdate//'_hammoz_global.nc'
    IF (ldebug) CALL message('', 'Output to file '//TRIM(cfilename), level=em_info)

    !-- create netcdf file
    CALL nf_check( nf_create(TRIM(cfilename), nf_clobber, nfid), fname=TRIM(cfilename) )
    !-- define dimension
    CALL nf_check( nf_def_dim(nfid, 'time', nf_unlimited, nfdim_time) )
    !-- define variables
!!### ADD YOUR GLOBAL DIAGNOSTIC VARIABLES HERE !! ###
    CALL nf_check( nf_def_var(nfid, 'time', nf_double, 1, (/nfdim_time/), nfvar_time) )
    write(crefdate, '(A,i4.4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2)')       &
          'days since ', year, month, day, hour, minute 
    CALL nf_check( nf_put_att_text(nfid, nfvar_time, 'units', 27, crefdate) )
    CALL define_variable( dpohconc, 'OH_conc', nfid, nfdim_time, GWAVG, 1._dp )
    CALL define_variable( dpratech4, 'rate_CH4', nfid, nfdim_time, GWAVG, 1._dp )
    CALL define_variable( dpratemcl, 'rate_CH3CCl3', nfid, nfdim_time, GWAVG, 1._dp )
    CALL define_variable( dpff,     'flash_frequency', nfid, nfdim_time, GSUM, 1._dp )
    CALL define_variable( dpemino,  'emi_NO_lghtng', nfid, nfdim_time, GSUM, 1._dp )
    !-- assign global attributes
    CALL nf_check(nf_put_att_text(nfid, nf_global, 'Conventions', 6, 'CF-1.5'))
    CALL nf_check(nf_put_att_text(nfid, nf_global, 'title', 18, 'Global Diagnostics'))
    CALL nf_check(nf_put_att_text(nfid, nf_global, 'source', 13, 'ECHAM6-HAMMOZ'))

    !-- end definition mode and close netcdf file
    CALL nf_check( nf_enddef(nfid) )
    CALL nf_check( nf_close(nfid) )

    !-- set buffer to zero
    dbuf(:) = 0._dp

    !-- compute global gridbox area
    garea = sum_global(gboxarea_2d)

  END SUBROUTINE init_hammoz_gdiag


  SUBROUTINE define_variable_2d( dpointer, vname, nfid, nfdim, nproc, factor )

    ! defines netcdf variables and associated the diagnostic pointers

    USE mo_netcdf,           ONLY: nf_check, nf_def_var, nf_real
    USE mo_exception,        ONLY: message, em_error

    REAL(dp), POINTER             :: dpointer(:,:)      ! reference to stream element
    CHARACTER(len=*), INTENT(in)  :: vname              ! variable name
    INTEGER, INTENT(in)           :: nfid               ! netcdf file ID
    INTEGER, INTENT(in)           :: nfdim              ! id of time dimension in netcdf file
    INTEGER, INTENT(in)           :: nproc              ! type of processing required
    REAL(dp), INTENT(in)          :: factor             ! scale factor for output
    INTEGER                       :: nftmp              ! temporary variable id in netcdf file

    IF (associated(dpointer)) THEN
      CALL nf_check( nf_def_var(nfid, TRIM(vname), nf_real, 1, (/nfdim/), nftmp) )
      ndiag = ndiag + 1
      dps(ndiag)%ndims = 2
      dps(ndiag)%key   = nftmp
      dps(ndiag)%key_type = nproc
      dps(ndiag)%fld2d => dpointer
      scalefac(ndiag) = factor

      IF (nproc == GTAU) CALL message('init_hammoz_gdiag',       &
                                      'Lifetime diagnostics not defined for 2D field '//TRIM(vname), &
                                      level=em_error)
    END IF

  END SUBROUTINE define_variable_2d

  SUBROUTINE define_variable_3d( dpointer, vname, nfid, nfdim, nproc, factor )

    ! defines netcdf variables and associated the diagnostic pointers

    USE mo_netcdf,           ONLY: nf_check, nf_def_var, nf_real

    REAL(dp), POINTER             :: dpointer(:,:,:)    ! reference to stream element
    CHARACTER(len=*), INTENT(in)  :: vname              ! variable name
    INTEGER, INTENT(in)           :: nfid               ! netcdf file ID
    INTEGER, INTENT(in)           :: nfdim              ! id of time dimension in netcdf file
    INTEGER, INTENT(in)           :: nproc              ! type of processing required
    REAL(dp), INTENT(in)          :: factor             ! scale factor for output
    INTEGER                       :: nftmp              ! temporary variable id in netcdf file

    IF (associated(dpointer)) THEN
      CALL nf_check( nf_def_var(nfid, TRIM(vname), nf_real, 1, (/nfdim/), nftmp) )
      ndiag = ndiag + 1
      dps(ndiag)%ndims = 3
      dps(ndiag)%key   = nftmp
      dps(ndiag)%key_type = nproc
      dps(ndiag)%fld3d => dpointer
      scalefac(ndiag) = factor
    END IF

  END SUBROUTINE define_variable_3d


  SUBROUTINE hammoz_gdiag

    USE mo_netcdf,              ONLY: nf_check, nf_open, nf_write, nf_put_vara_double,  &
                                      nf_close
    USE mo_control,             ONLY: ngl, nlon, nlev
    USE mo_time_control,        ONLY: time_step_len, current_date, get_date_components
    USE mo_geoloc,              ONLY: gboxarea_2d
    USE mo_global_op,           ONLY: sum_global
!##DEBUG##
    USE mo_mpi,                 ONLY: p_parallel_io

!   INTEGER              :: year, month, day, hour, minute, second
    INTEGER              :: nfid
    INTEGER              :: i, jk
    INTEGER, SAVE        :: nstep = 0
    INTEGER, SAVE        :: noutstep = 0
    INTEGER, PARAMETER   :: vcount = 1
    REAL(dp), SAVE       :: curdate = 0._dp         ! can be improved by providing real dates ###
    REAL(dp)             :: zbuf

!   CALL get_date_components(current_date, year = year, month = month, day = day,   &
!                            hour = hour, minute = minute, second = second)
    !-- advance time stamp and time step
    curdate = curdate + time_step_len
    nstep = nstep + 1
    !-- evaluate diagnostics
    DO i=1, ndiag
      zbuf = 0._dp
      SELECT CASE ( dps(i)%key_type )
        CASE (GSUM)
          IF ( dps(i)%ndims == 2 ) THEN
            dbuf(i) = dbuf(i) + sum_global(dps(i)%fld2d)
          ELSE
            DO jk=1, nlev
              zbuf = zbuf + sum_global(dps(i)%fld3d(:,jk,:))
            END DO
            dbuf(i) = dbuf(i) + zbuf
          END IF

        CASE (GAVG)
          IF ( dps(i)%ndims == 2 ) THEN
            dbuf(i) = dbuf(i) + sum_global(dps(i)%fld2d)/real(ngl*nlon,dp)
          ELSE
            DO jk=1, nlev
              zbuf = zbuf + sum_global(dps(i)%fld3d(:,jk,:))
            END DO
            dbuf(i) = dbuf(i) + zbuf/real(ngl*nlon*nlev,dp)
          END IF

        CASE (GWAVG)
          IF ( dps(i)%ndims == 2 ) THEN
            zbuf = sum_global(dps(i)%fld2d*gboxarea_2d)/garea
            dbuf(i) = dbuf(i) + sum_global(dps(i)%fld2d)/real(ngl*nlon,dp)
          ELSE
            zbuf = 0._dp
            ! ### WARNING: no vertical weight yet !! ###
            DO jk=1, nlev
              zbuf = zbuf + sum_global(dps(i)%fld3d(:,jk,:)*gboxarea_2d)
            END DO
            dbuf(i) = dbuf(i) + zbuf/garea/real(nlev,dp)
          END IF
    
        CASE (GTAU)
          ! (only 3D)
          zbuf = 0._dp
          ! ### WARNING: no vertical weight yet !! ###
          DO jk=1, nlev
            zbuf = zbuf + sum_global(gboxarea_2d/(dps(i)%fld3d(:,jk,:)))
!! if (p_parallel_io) write(0,*) '##DEBUG## jk, zbuf, min/max(term) = ',  &
!! jk,zbuf,minval(dps(i)%fld3d(:,jk,:)),&
!! maxval(dps(i)%fld3d(:,jk,:)),minval(gboxarea_2d),maxval(gboxarea_2d)
          END DO
          dbuf(i) = dbuf(i) + garea*real(nlev,dp)/zbuf
if (p_parallel_io) write(0,*) '##DEBUG## i, dbuf(i), garea, real(nlev,dp) =',i, dbuf(i), garea, real(nlev,dp)
      END SELECT
    END DO

    !-- decide if output needs to be written
    IF (MOD(nstep, navg) == 0) THEN
      CALL nf_check( nf_open(TRIM(cfilename), nf_write, nfid) )
      !-- write time stamp
      noutstep = noutstep + 1
      CALL nf_check( nf_put_vara_double(nfid, nfvar_time, noutstep, vcount, curdate/real(navg,dp)) )
      DO i=1, ndiag
      !-- write buffer value average
        CALL nf_check( nf_put_vara_double(nfid, dps(i)%key, noutstep, vcount, &
                                          scalefac(i)*dbuf(i)/real(navg,dp)) )
      END DO
      !-- close file
      CALL nf_check( nf_close(nfid) )
      !-- reset buffer
      dbuf(:) = 0._dp 
    END IF

  END SUBROUTINE hammoz_gdiag


END MODULE mo_hammoz_global_diag
