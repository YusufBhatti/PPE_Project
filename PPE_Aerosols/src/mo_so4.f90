!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_so4

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_so4nat, read_so4all, cleanup_so4
  PUBLIC :: so4nat_io
  PUBLIC :: so4all_io

  !=======================================================================

  ! unit number of sulfate file
  INTEGER, PARAMETER :: niso4nat=18

  REAL(dp), ALLOCATABLE  :: so4nat_io(:,:,:,:) 
  REAL(dp), ALLOCATABLE  :: so4all_io(:,:,:,:) 

  !=======================================================================

CONTAINS

  SUBROUTINE read_so4nat

    USE mo_control,       ONLY: ngl, nlon, nlev
    USE mo_io_units,        ONLY: nout, nerr
    USE mo_mpi,           ONLY: p_pe, p_io
    USE mo_exception,     ONLY: finish
    USE mo_io,            ONLY: io_open_unit, io_close, io_read, &
                                io_var_id, io_file_id, ini_field
    USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
                                io_inq_varid, io_get_vara_double
    USE mo_decomposition, ONLY: lc => local_decomposition, &
                                gl_dc => global_decomposition
    USE mo_transpose,     ONLY: scatter_gp

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:,:)
    REAL(dp), POINTER :: gl_so4(:,:,:,:)

    INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
    INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
    INTEGER               :: io_nlev  ! number of levels in NetCDF file
    INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
    INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read

    INTEGER               :: jk ,i     ! loop index

    ! Read sulfate file
    ! ===============

    IF (p_pe==p_io) THEN

       WRITE(nout,'(/,A,I2)') ' Read sulfate climatology from unit ', niso4nat

       CALL io_open_unit (niso4nat, ini_field, io_read)
       io_file_id = ini_field%file_id

       ! Check resolution
       CALL io_inq_dimid  (io_file_id, 'lat', io_var_id)
       CALL io_inq_dimlen (io_file_id, io_var_id, io_ngl)
       CALL io_inq_dimid  (io_file_id, 'lon', io_var_id)
       CALL io_inq_dimlen (io_file_id, io_var_id, io_nlon)
       CALL io_inq_dimid  (io_file_id, 'lev', io_var_id)
       CALL io_inq_dimlen (io_file_id, io_var_id, io_nlev)

       IF (io_ngl/=ngl) THEN
          WRITE(nerr,*) 'read_so4: unexpected resolution ',io_nlon,io_ngl
          WRITE(nerr,*) 'expected number of latitudes = ',ngl
          WRITE(nerr,*) 'number of latitudes of sulfate = ',io_ngl
          CALL finish ('read_so4','unexpected resolution')
       END IF

       IF (io_nlev/=nlev) THEN
          WRITE(nerr,*) 'read_so4: unexpected number of levels ',io_nlev
          WRITE(nerr,*) 'expected number of levels = ',nlev
          CALL finish ('read_so4','unexpected resolution')
       END IF

    END IF

    !     Allocate memory for sst per PE

    IF ( .NOT. ALLOCATED(so4nat_io)) ALLOCATE(so4nat_io(lc%nproma,nlev,lc%ngpblks,0:13))

    IF (p_pe == p_io) THEN

      !     Allocate memory for so4 global fields

      ALLOCATE (zin(nlon,nlev,ngl,0:13))

      ! read sulfate
      CALL io_inq_varid (io_file_id, 'sulfate', io_var_id)
      DO jk = 1, nlev
        io_start(:) = (/       1,   1, jk,  1 /)
        io_count(:) = (/ io_nlon, ngl,  1, 12 /)
        ! for level jk: read io_nlon longitudes, ngl latitudes and 12 months
        CALL io_get_vara_double(io_file_id,io_var_id,io_start,io_count, &
             &                  zin(1:io_nlon,jk,:,1:12))
      END DO

      CALL io_close (ini_field)

      ! copy December to month 0
      zin(:,:,:,0)  = zin(:,:,:,12)

      ! copy January to month 13
      zin(:,:,:,13) = zin(:,:,:,1)

    ENDIF

    NULLIFY (gl_so4)
    DO i = 0, 13
      IF (p_pe == p_io) gl_so4 => zin(:,:,:,i:i)
      CALL scatter_gp(gl_so4, so4nat_io(:,:,:,i:i), gl_dc)
    END DO

    IF (p_pe == p_io) THEN
      DEALLOCATE (zin)
    END IF

  END SUBROUTINE read_so4nat


  SUBROUTINE read_so4all

    USE mo_control,       ONLY: ngl, nlon, nlev
    USE mo_io_units,        ONLY: nout, nerr
    USE mo_mpi,           ONLY: p_pe, p_io
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io,            ONLY: io_open, io_close, io_read, &
                                io_var_id, sstnc0, sstnc1, sstnc2
    USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
                                io_inq_varid, io_get_vara_double
    USE mo_decomposition, ONLY: lc => local_decomposition, &
                                gl_dc => global_decomposition
    USE mo_transpose,     ONLY: scatter_gp

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:,:)
    REAL(dp), POINTER :: gl_so4(:,:,:,:)

    INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
    INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
    INTEGER               :: io_nlev  ! number of levels in NetCDF file
    INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
    INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read

    INTEGER :: IO_file_id0, IO_file_id1, IO_file_id2
    INTEGER               :: jk ,i     ! loop index
    CHARACTER (8) :: fn0, fn1, fn2
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex0, lex1, lex2

    ! Read sulfate file
    ! ===============

    IF (p_pe==p_io) THEN

       CALL set_years(ihy0, ihy1, ihy2)

       WRITE (fn0, '("aero",i4)') ihy0
       WRITE (fn1, '("aero",i4)') ihy1
       WRITE (fn2, '("aero",i4)') ihy2

       WRITE (nout, '(/)')

       WRITE(message_text,*) 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), ' fn2: ',TRIM(fn2)
       CALL message('read_so4all',message_text)

       INQUIRE (file=fn0, exist=lex0)
       INQUIRE (file=fn1, exist=lex1)
       INQUIRE (file=fn2, exist=lex2)

       IF ( .NOT. lex0 ) THEN
         WRITE (message_text,*) 'Could not open file <',fn0,'>'
         CALL message('',message_text)
         CALL finish ('read_so4all', 'run terminated.')
       END IF

       IF ( .NOT. lex1 ) THEN
         WRITE (message_text,*) 'Could not open file <',fn1,'>'
         CALL message('',message_text)
         CALL finish ('read_so4all', 'run terminated.')
       END IF

       IF ( .NOT. lex2 ) THEN
         WRITE (message_text,*) 'Could not open file <',fn2,'>'
         CALL message('',message_text)
         CALL finish ('read_so4all', 'run terminated.')
       END IF
        
       CALL IO_open (fn0, sstnc0, IO_READ)
       CALL IO_open (fn1, sstnc1, IO_READ)
       CALL IO_open (fn2, sstnc2, IO_READ)

       io_file_id0 = sstnc0%file_id
       io_file_id1 = sstnc1%file_id
       io_file_id2 = sstnc2%file_id

       ! Check resolution
       CALL io_inq_dimid  (io_file_id1, 'lat', io_var_id)
       CALL io_inq_dimlen (io_file_id1, io_var_id, io_ngl)
       CALL io_inq_dimid  (io_file_id1, 'lon', io_var_id)
       CALL io_inq_dimlen (io_file_id1, io_var_id, io_nlon)
       CALL io_inq_dimid  (io_file_id1, 'mlev', io_var_id)
       CALL io_inq_dimlen (io_file_id1, io_var_id, io_nlev)

       IF (io_ngl/=ngl) THEN
          WRITE(nerr,*) 'read_so4: unexpected resolution ',io_nlon,io_ngl
          WRITE(nerr,*) 'expected number of latitudes = ',ngl
          WRITE(nerr,*) 'number of latitudes of sulfate = ',io_ngl
          CALL finish ('read_so4','unexpected resolution')
       END IF

       IF (io_nlev/=nlev) THEN
          WRITE(nerr,*) 'read_so4: unexpected number of levels ',io_nlev
          WRITE(nerr,*) 'expected number of levels = ',nlev
          CALL finish ('read_so4','unexpected resolution')
       END IF

    END IF

    !     Allocate memory for sst per PE

    IF ( .NOT. ALLOCATED(so4all_io) ) ALLOCATE(so4all_io(lc%nproma,nlev,lc%ngpblks,0:13))

    IF (p_pe == p_io) THEN

      !     Allocate memory for so4 global fields
       
      ALLOCATE (zin(nlon,nlev,ngl,0:13))

      ! read sulfate 12 month
      CALL io_inq_varid (io_file_id1, 'sulfate', io_var_id)
      DO jk = 1, nlev
        io_start(:) = (/       1,   1, jk,  1 /)
        io_count(:) = (/ io_nlon, ngl,  1, 12 /)
        ! for level jk: read io_nlon longitudes, ngl latitudes and 12 months
        CALL io_get_vara_double(io_file_id1,io_var_id,io_start,io_count, &
             &                  zin(1:io_nlon,jk,:,1:12))
      END DO

      ! read sulfate december of last year
      CALL io_inq_varid (io_file_id0, 'sulfate', io_var_id)
      DO jk = 1, nlev
        io_start(:) = (/       1,   1, jk, 12 /)
        io_count(:) = (/ io_nlon, ngl,  1,  1 /)
        ! for level jk: read io_nlon longitudes, ngl latitudes and 1 months
        CALL io_get_vara_double(io_file_id0,io_var_id,io_start,io_count, &
             &                  zin(1:io_nlon,jk,:,0))
      END DO

      ! read sulfate january of next year
      CALL io_inq_varid (io_file_id2, 'sulfate', io_var_id)
      DO jk = 1, nlev
        io_start(:) = (/       1,   1, jk,  1 /)
        io_count(:) = (/ io_nlon, ngl,  1,  1 /)
        ! for level jk: read io_nlon longitudes, ngl latitudes and 1 months
        CALL io_get_vara_double(io_file_id2,io_var_id,io_start,io_count, &
             &                  zin(1:io_nlon,jk,:,13))
      END DO

      CALL io_close (sstnc0)
      CALL io_close (sstnc1)
      CALL io_close (sstnc2)

    ENDIF

    NULLIFY (gl_so4)
    DO i = 0, 13
      IF (p_pe == p_io) gl_so4 => zin(:,:,:,i:i)
      CALL scatter_gp(gl_so4, so4all_io(:,:,:,i:i), gl_dc)
    END DO

    IF (p_pe == p_io) THEN
      DEALLOCATE (zin)
    END IF

  END SUBROUTINE read_so4all

  SUBROUTINE set_years(y1,y2,y3)

    USE mo_time_control,    ONLY: next_date, get_date_components

    INTEGER, INTENT(out) :: y1, y2, y3
    INTEGER :: yr, mo, dy, hr, mn, se

    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)

    y1 = yr - 1
    y2 = yr
    y3 = yr + 1

  END SUBROUTINE set_years

!------------------------------------------------------------------------------
  SUBROUTINE cleanup_so4
    !----------------------------
    ! deallocate module variables
    !----------------------------
    IF (ALLOCATED(so4nat_io)) DEALLOCATE (so4nat_io)
    IF (ALLOCATED(so4all_io)) DEALLOCATE (so4all_io)

  END SUBROUTINE cleanup_so4
!------------------------------------------------------------------------------
END MODULE mo_so4
