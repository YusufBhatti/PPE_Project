!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_sst

  USE mo_kind,            ONLY: dp
  USE mo_time_control,    ONLY: next_date, get_date_components
  USE mo_netCDF,          ONLY: IO_info_print
  USE mo_decomposition,   ONLY: lc => local_decomposition, global_decomposition

  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: dailysst(:,:,:)  ! (nlon,ngl,3) in global coordinates
  REAL(dp), ALLOCATABLE :: dailyice(:,:,:)  ! (nlon,ngl,3) in global coordinates
  REAL(dp), ALLOCATABLE :: sst(:,:,:)  ! (nlon,ngl,0:13) in global coordinates
  REAL(dp), ALLOCATABLE :: aice(:,:,:)  ! (nlon,ngl,0:13) in global coordinates
  REAL(dp), ALLOCATABLE :: aflux(:,:,:) ! (nlon,ngl,0:13) in global coordinates

CONTAINS

  SUBROUTINE readsst

    ! U. Schlese, DKRZ,  May 1993, original version
    ! U. Schulzweida, MPI, May 1999, netCDF version
    ! L. Kornblueh, MPI, November 2001, cleanup for parallel environment
    ! U. Schulzweida, MPI, May 2002, blocking (nproma)

    USE mo_control,       ONLY: lamip, lmlo, nist
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
   
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_sst(:,:,:)

    CHARACTER (12) :: fn0, fn1, fn2
    INTEGER       :: i, iy
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex, lex0, lex1, lex2
    INTEGER       :: start(3), COUNT(3), nvarid

    CALL set_years(ihy0, ihy1, ihy2)
    iy = ihy1

    IF (iy < 0) THEN
       WRITE (fn0, '("sst",i0)') ihy0
       WRITE (fn1, '("sst",i0)') ihy1
       WRITE (fn2, '("sst",i0)') ihy2
    ELSE IF (iy < 100) THEN
       WRITE (fn0, '("sst",i2.2)') ihy0
       WRITE (fn1, '("sst",i2.2)') ihy1
       IF(iy/= 99) THEN
          WRITE (fn2, '("sst",i2.2)') ihy2
       ELSE
          WRITE (fn2, '("sst",i3)') ihy2
       ENDIF
    ELSE IF (iy< 1000) THEN
       IF (iy/= 100) THEN
          WRITE (fn0, '("sst",i3)') ihy0
       ELSE
          WRITE (fn0, '("sst",i2.2)') ihy0
       ENDIF
       WRITE (fn1, '("sst",i3)') ihy1
       IF(iy/= 999) THEN
          WRITE (fn2, '("sst",i3)') ihy2
       ELSE
          WRITE (fn2, '("sst",i4)') ihy2
       ENDIF
    ELSE
       IF(iy/= 1000) THEN
          WRITE (fn0, '("sst",i4)') ihy0
       ELSE
          WRITE (fn0, '("sst",i3)') ihy0
       ENDIF
       WRITE (fn1, '("sst",i4)') ihy1
       WRITE (fn2, '("sst",i4)') ihy2
    ENDIF

    WRITE(message_text,'(7a,i4)') 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), &
                         ' fn2: ',TRIM(fn2),' nist: ',nist
    CALL message('readsst',message_text)

    ! Amip-type:

    IF (p_parallel_io) THEN
      IF(lamip .AND. .NOT. lmlo) THEN
        CALL message('','This is an AMIP run (lamip = .true.).')
        INQUIRE (file=fn0, exist=lex0)
        INQUIRE (file=fn1, exist=lex1)
        INQUIRE (file=fn2, exist=lex2)
        IF (lex1) THEN
          CALL IO_open (fn1, sstnc1, IO_READ)
          WRITE (message_text,*) 'Reading sst from files ',&
               fn0, ', ',fn1,', ',fn2
          CALL message('',message_text)
          IF(lex0) THEN
            CALL IO_open (fn0, sstnc0, IO_READ)
          ELSE
            WRITE (message_text,*) 'Could not open file <',fn0,'>'
            CALL message('',message_text)
            CALL finish ('readsst', 'run terminated.')
          ENDIF
          IF(lex2) THEN
            CALL IO_open (fn2, sstnc2, IO_READ)
          ELSE
            WRITE (message_text,*) 'Could not open file <',fn2,'>'
            CALL message('',message_text)
            CALL finish ('readsst', 'run terminated.')
          ENDIF
        ELSE
          WRITE (message_text,*) 'Could not open file <',fn1,'>'
          CALL message('',message_text)
          CALL finish ('readsst', 'run terminated.')
        ENDIF
      ELSE
        CALL message('','This is no AMIP run (lamip = .false.).')
        INQUIRE (nist, exist=lex)
        WRITE(message_text,*) 'lex: ', lex
        CALL message('readsst',message_text)
        IF (lex) THEN
          sstnc1%format = NETCDF
          CALL IO_open_unit (nist, sstnc1, IO_READ)
          ! has to be fixed...
          !          CALL IO_read_header(sstnc1)
          !          CALL IO_info_print(sstnc1)
        ELSE
          CALL finish ('readsst', 'Could not open sst file')
        ENDIF
      ENDIF
    ENDIF
    
    !     Allocate memory for sst per PE

    IF (.NOT. ALLOCATED(sst)) ALLOCATE (sst(lc%nproma, lc%ngpblks,0:13))

    !     Read sst-file
    IF (p_parallel_io) THEN

      !     Allocate memory for sst global fields
       
      ALLOCATE (zin(lc%nlon,lc%nlat,0:13))

      CALL IO_INQ_VARID (sstnc1%file_id, 'sst', nvarid)
      CALL IO_GET_VAR_DOUBLE (sstnc1%file_id, nvarid, zin(:,:,1:12))

      IF(.NOT.lamip .OR. lmlo) THEN
        zin(:,:,0)  = zin(:,:,12)
        zin(:,:,13) = zin(:,:,1)
      ELSE 
        CALL IO_INQ_VARID (sstnc0%file_id, 'sst', nvarid)
        COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 12 /)
        CALL IO_GET_VARA_DOUBLE (sstnc0%file_id,nvarid,start,count, &
             zin(1,1,0))

        CALL IO_INQ_VARID (sstnc2%file_id, 'sst', nvarid)
        COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 1 /)
        CALL IO_GET_VARA_DOUBLE (sstnc2%file_id,nvarid,start,count, &
             zin(1,1,13))
      END IF
    END IF

    NULLIFY (gl_sst)
    DO i = 0, 13
      IF (p_pe == p_io) gl_sst => zin(:,:,i:i)
      CALL scatter_gp (gl_sst, sst(:,:,i:i), global_decomposition)
    END DO

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)

       !    Close file(s)

       CALL IO_close(sstnc1)

       IF(lamip) THEN
         CALL IO_close(sstnc0)
         CALL IO_close(sstnc2)
       ENDIF
     ENDIF

  END SUBROUTINE readsst

  SUBROUTINE readice

    ! U. Schlese, DKRZ,  May 1993, original version
    ! L. Kornblueh, MPI, November 2001, cleanup for parallel environment

    USE mo_control,       ONLY: lamip, nice
    USE mo_exception,     ONLY: finish
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
    
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_ice(:,:,:)

    CHARACTER (12) :: fn0, fn1, fn2
    INTEGER       :: i, iy
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex, lex0, lex1, lex2
    INTEGER       :: start(3), COUNT(3), nvarid

    CALL set_years(ihy0, ihy1, ihy2)
    iy = ihy1

    IF (iy < 0) THEN
       WRITE (fn0, '("ice",i0)') ihy0
       WRITE (fn1, '("ice",i0)') ihy1
       WRITE (fn2, '("ice",i0)') ihy2
    ELSE IF (iy < 100) THEN
       WRITE (fn0, '("ice",i2.2)') ihy0
       WRITE (fn1, '("ice",i2.2)') ihy1
       IF(iy/= 99) THEN
          WRITE (fn2, '("ice",i2.2)') ihy2
       ELSE
          WRITE (fn2, '("ice",i3)') ihy2
       ENDIF
    ELSE IF (iy< 1000) THEN
       IF (iy/= 100) THEN
          WRITE (fn0, '("ice",i3)') ihy0
       ELSE
          WRITE (fn0, '("ice",i2.2)') ihy0
       ENDIF
       WRITE (fn1, '("ice",i3)') ihy1
       IF(iy/= 999) THEN
          WRITE (fn2, '("ice",i3)') ihy2
       ELSE
          WRITE (fn2, '("ice",i4)') ihy2
       ENDIF
    ELSE
       IF(iy/= 1000) THEN
          WRITE (fn0, '("ice",i4)') ihy0
       ELSE
          WRITE (fn0, '("ice",i3)') ihy0
       ENDIF
       WRITE (fn1, '("ice",i4)') ihy1
       WRITE (fn2, '("ice",i4)') ihy2
    ENDIF

    WRITE(message_text,*) 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), &
                          ' fn2: ',TRIM(fn2),' nice: ',nice
    CALL message('readice',message_text)

    ! Amip-type:

    IF (p_parallel_io) THEN
      IF(lamip) THEN
        CALL message('','This is an AMIP run (lamip = .true.).')
        INQUIRE (file=fn0, exist=lex0)
        INQUIRE (file=fn1, exist=lex1)
        INQUIRE (file=fn2, exist=lex2)
        IF (lex1) THEN
          CALL IO_open (fn1, icenc1, IO_READ)
          WRITE (message_text,*) 'Reading ice from files ',fn0, ', ', &
               fn1,', ',fn2
          CALL message('',message_text)
          IF(lex0) THEN
            CALL IO_open (fn0, icenc0, IO_READ)
          ELSE
            WRITE (message_text,*) 'Could not open file <',fn0,'>'
            CALL message('',message_text)
            CALL finish ('readice', 'run terminated.')
          ENDIF
          IF(lex2) THEN
            CALL IO_open (fn2, icenc2, IO_READ)
          ELSE
            WRITE (message_text,*) 'Could not open file <',fn2,'>'
            CALL message('',message_text)
            CALL finish ('readice', 'run terminated.')
          ENDIF
        ELSE
          WRITE (message_text,*) 'Could not open file <',fn1,'>'
          CALL message('',message_text)
          CALL finish ('readice', 'run terminated.')
        ENDIF
      ELSE
        CALL message('','This is no AMIP run (lamip = .false.).')
        
        INQUIRE (nice, exist=lex)
        WRITE(message_text,*) 'lex: ', lex
        CALL message('readice',message_text)
        IF (lex) THEN
          icenc1%format = NETCDF
          CALL IO_open_unit (nice, icenc1, IO_READ)
          ! has to be fixed...
          !          CALL IO_read_header(icenc1)
          !          CALL IO_info_print(icenc1)
        ELSE
          CALL finish ('readice', 'Could not open ice file')
        ENDIF
      ENDIF
    END IF

    !     Allocate memory for ice per PE

    IF (.NOT. ALLOCATED(aice)) ALLOCATE (aice(lc%nproma, lc%ngpblks,0:13))

    !     Read ice-file
    IF (p_parallel_io) THEN

      !     Allocate memory for ice global fields
       
      ALLOCATE (zin(lc%nlon,lc%nlat,0:13))
      
      CALL IO_INQ_VARID (icenc1%file_id, 'sic', nvarid)
      CALL IO_GET_VAR_DOUBLE (icenc1%file_id, nvarid, zin(:,:,1:12))
      
      IF(.NOT.lamip) THEN
        zin(:,:,0)  = zin(:,:,12)
        zin(:,:,13) = zin(:,:,1)
      ELSE 
        CALL IO_INQ_VARID (icenc0%file_id, 'sic', nvarid)
        COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 12 /)
        CALL IO_GET_VARA_DOUBLE (icenc0%file_id,nvarid,start,count, &
             zin(1,1,0))
        
        CALL IO_INQ_VARID (icenc2%file_id, 'sic', nvarid)
        COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 1 /)
        CALL IO_GET_VARA_DOUBLE (icenc2%file_id,nvarid,start,count, &
             zin(1,1,13))
      END IF
    END IF

    NULLIFY (gl_ice)
    DO i = 0, 13
      IF (p_parallel_io) gl_ice => zin(:,:,i:i)
      CALL scatter_gp (gl_ice, aice(:,:,i:i), global_decomposition)
    END DO
    
    IF (p_parallel_io) THEN
      DEALLOCATE (zin)

      !    Close file(s)
      
      CALL IO_close(icenc1)
      
      IF(lamip) THEN
        CALL IO_close(icenc0)
        CALL IO_close(icenc2)
      ENDIF
    ENDIF

  END SUBROUTINE readice

  SUBROUTINE readflux

    ! U. Schlese, DKRZ,  May 1993, original version (readsst)
    ! M. Esch,    MPI,   Sep 2002, modified for flux correction

    USE mo_control,       ONLY: nflu
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
   
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_aflux(:,:,:)

    CHARACTER (7) :: fn0
    INTEGER       :: i, iy
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex
    INTEGER       :: nvarid

    CALL set_years(ihy0, ihy1, ihy2)
    iy = ihy1

    IF (iy < 100) THEN
       WRITE (fn0, '("flu",i2.2)') ihy0
    ELSE IF (iy< 1000) THEN
       IF (iy/= 100) THEN
          WRITE (fn0, '("flu",i3)') ihy0
       ELSE
          WRITE (fn0, '("flu",i2.2)') ihy0
       ENDIF
    ELSE
       IF(iy/= 1000) THEN
          WRITE (fn0, '("flu",i4)') ihy0
       ELSE
          WRITE (fn0, '("flu",i3)') ihy0
       ENDIF
    ENDIF

    WRITE(message_text,*) 'fn0: ', TRIM(fn0),' nflu: ',nflu
    CALL message('readflux',message_text)

    ! 

    IF (p_parallel_io) THEN
        CALL message('','This is an MLO run (lmlo = .true.).')
        INQUIRE (nflu, exist=lex)
        WRITE(message_text,*) 'lex: ', lex
        CALL message('readflux',message_text)
        IF (lex) THEN
          flunc1%format = NETCDF
          CALL IO_open_unit (nflu, flunc1, IO_READ)
          ! has to be fixed...
          !          CALL IO_read_header(flunc1)
          !          CALL IO_info_print(flunc1)
        ELSE
          CALL finish ('readflux', 'Could not open flux file')
        ENDIF
    ENDIF
    
    !     Allocate memory for aflux per PE

    IF (.NOT. ALLOCATED(aflux)) ALLOCATE (aflux(lc%nproma, lc%ngpblks,0:13))

    !     Read flux-file
    IF (p_parallel_io) THEN

      !     Allocate memory for flux global fields
       
      ALLOCATE (zin(lc%nlon,lc%nlat,0:13))

      CALL IO_INQ_VARID (flunc1%file_id, 'aflux', nvarid)
      CALL IO_GET_VAR_DOUBLE (flunc1%file_id, nvarid, zin(:,:,1:12))

        zin(:,:,0)  = zin(:,:,12)
        zin(:,:,13) = zin(:,:,1)
    END IF

    NULLIFY (gl_aflux)
    DO i = 0, 13
      IF (p_pe == p_io) gl_aflux => zin(:,:,i:i)
      CALL scatter_gp (gl_aflux, aflux(:,:,i:i), global_decomposition)
    END DO

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)

       !    Close file(s)

       CALL IO_close(flunc1)

     ENDIF

  END SUBROUTINE readflux


  SUBROUTINE readdailysst

    ! U. Schulzweida, MPI, March 2007

    USE mo_control,       ONLY: nist
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
   
    REAL(dp), ALLOCATABLE :: timevals(:)
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_sst(:,:,:)

    INTEGER :: i
    LOGICAL :: lex
    INTEGER :: start(3), COUNT(3), nvarid, ndimid, nts, tsID
    INTEGER :: yr, mo, dy, hr, mn, se, ydate

    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)

    ydate = yr*10000+mo*100+dy

    IF (p_parallel_io) THEN
      WRITE(message_text,*) ' date: ', ydate, ' nist: ',nist
      CALL message('readdailysst',message_text)
    ENDIF

    IF (p_parallel_io) THEN
      INQUIRE (nist, exist=lex)
      IF (lex) THEN
        sstnc1%format = NETCDF
        CALL IO_open_unit (nist, sstnc1, IO_READ)
      ELSE
        CALL finish ('readdailysst', 'Could not open sst file')
      ENDIF
    ENDIF
    
    !     Allocate memory for sst per PE

    IF (.NOT. ALLOCATED(dailysst)) ALLOCATE (dailysst(lc%nproma, lc%ngpblks,3))

    !     Read sst-file
    IF (p_parallel_io) THEN

      CALL IO_INQ_DIMID (sstnc1%file_id, 'time', ndimid)

      CALL IO_INQ_DIMLEN (sstnc1%file_id, ndimid, nts)

      CALL IO_INQ_VARID (sstnc1%file_id, 'time', nvarid)

      IF ( nts .lt. 3 ) CALL finish('readdailysst', 'To few time steps')

      ALLOCATE (timevals(nts))

      CALL IO_GET_VAR_DOUBLE(sstnc1%file_id, nvarid, timevals)

      DO tsID = 2, nts-1
        IF ( ydate .eq. INT(timevals(tsID)) ) exit;
      END DO

      DEALLOCATE (timevals)

      IF ( tsID .ge. nts ) CALL finish('readdailysst', 'Date not found')

      !     Allocate memory for sst global fields
       
      ALLOCATE (zin(lc%nlon,lc%nlat,3))

      CALL IO_INQ_VARID (sstnc1%file_id, 'sst', nvarid)
      COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
      start(:) = (/ 1, 1, tsID-1 /)
      CALL IO_GET_VARA_DOUBLE (sstnc1%file_id,nvarid,start,count, zin(1,1,1))
      start(:) = (/ 1, 1, tsID /)
      CALL IO_GET_VARA_DOUBLE (sstnc1%file_id,nvarid,start,count, zin(1,1,2))
      start(:) = (/ 1, 1, tsID+1 /)
      CALL IO_GET_VARA_DOUBLE (sstnc1%file_id,nvarid,start,count, zin(1,1,3))
    END IF

    NULLIFY (gl_sst)
    DO i = 1, 3
      IF (p_pe == p_io) gl_sst => zin(:,:,i:i)
      CALL scatter_gp (gl_sst, dailysst(:,:,i:i), global_decomposition)
    END DO

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)

       !    Close file(s)

       CALL IO_close(sstnc1)
     ENDIF

  END SUBROUTINE readdailysst

  SUBROUTINE readdailyice

    ! U. Schulzweida, MPI, March 2007

    USE mo_control,       ONLY: nice
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
   
    REAL(dp), ALLOCATABLE :: timevals(:)
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_ice(:,:,:)

    INTEGER :: i
    LOGICAL :: lex
    INTEGER :: start(3), COUNT(3), nvarid, ndimid, nts, tsID
    INTEGER :: yr, mo, dy, hr, mn, se, ydate

    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)

    ydate = yr*10000+mo*100+dy

    IF (p_parallel_io) THEN
      WRITE(message_text,*) ' date: ', ydate, ' nice: ',nice
      CALL message('readdailyice',message_text)
    ENDIF

    IF (p_parallel_io) THEN
      CALL message('','This is no AMIP run (lamip = .false.).')
      INQUIRE (nice, exist=lex)
      IF (lex) THEN
        icenc1%format = NETCDF
        CALL IO_open_unit (nice, icenc1, IO_READ)
      ELSE
        CALL finish ('readdailyice', 'Could not open ice file')
      ENDIF
    ENDIF
    
    !     Allocate memory for ice per PE

    IF (.NOT. ALLOCATED(dailyice)) ALLOCATE (dailyice(lc%nproma, lc%ngpblks,3))

    !     Read ice-file
    IF (p_parallel_io) THEN

      CALL IO_INQ_DIMID (icenc1%file_id, 'time', ndimid)

      CALL IO_INQ_DIMLEN (icenc1%file_id, ndimid, nts)

      CALL IO_INQ_VARID (icenc1%file_id, 'time', nvarid)

      IF ( nts .lt. 3 ) CALL finish('readdailyice', 'To few time steps')

      ALLOCATE (timevals(nts))

      CALL IO_GET_VAR_DOUBLE(icenc1%file_id, nvarid, timevals)

      DO tsID = 2, nts-1
        IF ( ydate .eq. INT(timevals(tsID)) ) exit;
      END DO

      DEALLOCATE (timevals)

      IF ( tsID .ge. nts ) CALL finish('readdailyice', 'Date not found')

      !     Allocate memory for ice global fields
       
      ALLOCATE (zin(lc%nlon,lc%nlat,3))

      CALL IO_INQ_VARID (icenc1%file_id, 'sic', nvarid)
      COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
      start(:) = (/ 1, 1, tsID-1 /)
      CALL IO_GET_VARA_DOUBLE (icenc1%file_id,nvarid,start,count, zin(1,1,1))
      start(:) = (/ 1, 1, tsID /)
      CALL IO_GET_VARA_DOUBLE (icenc1%file_id,nvarid,start,count, zin(1,1,2))
      start(:) = (/ 1, 1, tsID+1 /)
      CALL IO_GET_VARA_DOUBLE (icenc1%file_id,nvarid,start,count, zin(1,1,3))
    END IF

    NULLIFY (gl_ice)
    DO i = 1, 3
      IF (p_pe == p_io) gl_ice => zin(:,:,i:i)
      CALL scatter_gp (gl_ice, dailyice(:,:,i:i), global_decomposition)
    END DO

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)

       !    Close file(s)

       CALL IO_close(icenc1)
     ENDIF

  END SUBROUTINE readdailyice


  SUBROUTINE set_years(y1,y2,y3)
    INTEGER, INTENT(out) :: y1, y2, y3
    INTEGER :: yr, mo, dy, hr, mn, se

    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)

    y1 = yr - 1
    y2 = yr
    y3 = yr + 1

  END SUBROUTINE set_years
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_sst
    !
    ! deallocate module variables
    !
    IF (ALLOCATED(dailysst)) DEALLOCATE (dailysst)
    IF (ALLOCATED(dailyice)) DEALLOCATE (dailyice)
    IF (ALLOCATED(sst))   DEALLOCATE (sst)
    IF (ALLOCATED(aice))  DEALLOCATE (aice)
    IF (ALLOCATED(aflux)) DEALLOCATE (aflux)
  END SUBROUTINE cleanup_sst
!------------------------------------------------------------------------------
END MODULE mo_sst
