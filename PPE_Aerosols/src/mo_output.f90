!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_output

  ! Description:
  !
  ! This module contains the output routines based on CDI version 1.4.3.
  !
  ! For a detailed description of the contained subroutines look on their
  ! header.
  !
  ! Most important new features:
  ! - Use of a code table in section 1 which allows the use of 128 tables
  !   each containg 128 different variables.
  ! - Uses the sub center entry (232 for ECHAM). Center is still 98 for
  !   ECMWF
  ! - The century parameter is correctly used. Now, in 360 day mode,
  !   25599 years can be simulated before an overflow in the date occures.
  ! - Uses correct dates in 365/366 day mode. Allows for clean forecast and
  !   nudging runs.
  ! - The output is now real standard !!!
  !
  ! Authors:
  !
  ! L. Kornblueh,   MPI, December 1998
  ! L. Kornblueh,   MPI, April    1999, added NWP forecast mode
  ! U. Schulzweida, MPI, April    2000, EMOS compatible
  ! I. Kirchner,    MPI, December 2000, time control
  ! A. Rhodin,      MPI, Mai      2001, condensed to one output routine
  ! A. Rhodin,  DWD/MPI, October  2001, NetCDF calls, parallel GRIB encoding
  ! L. Kornblueh,   MPI, October  2001, changed subcenter to 232
  !                                     (ECMWF assigned)
  ! A. Rhodin,  DWD/MPI, February 2001, bug fixes for parallel mode
  ! U. Schulzweida, MPI, May      2002, blocking (nproma)
  ! R. Johanni, IPP Garching, May-2002, parallel nudging
  ! I. Kirchner,    MPI, August   2002, lpout flag, add backup_output_streams
  ! U. Schulzweida, MPI, February 2003, change codegb5 interface to gribex
  ! A. Rhodin,      DWD, March    2003, bug fix for no_cycles > 1
  ! L. Kornblueh    MPI, April    2003, global communicator specified
  ! U. Schulzweida  MPI, November 2007, changed output to CDI library
  ! U. Schulzweida  MPI, March    2009, update for jsbach
  ! U. Schulzweida  MPI, January  2010, added support for seconds
  ! U. Schulzweida  MPI, February 2010, use complex packing for spherical harmonics
  ! L. Kornblueh    MPI, March    2010, new gather calls
  ! S. K. Cheedela  MPI, January  2010, support for output in single column
  ! D. Goll         MPI, February 2013, code table for yasso soil carbon

  USE mo_kind,                      ONLY: dp
  USE mo_exception,                 ONLY: finish, message, message_text
  USE mo_util_string,               ONLY: separator            ! format string (----)
  USE mo_time_control,              ONLY: ev_putdata, l_putdata, &
                                          get_interval_seconds, next_date, start_date, &
                                          get_date_components, get_forecast_hours,     &
                                          write_date
  USE mo_netcdf,                    ONLY: io_dim, IO_dim_ids, global_att,              &
                                          BELOWSUR, SURFACE, HYBRID, HYBRID_H,         &
                                          TILES, SOILLEV, ROOTZONES, CANOPY, IO_GET_VARINDX, & 
                                          COSP_LEVELS !ZK #474
  USE mo_decomposition,             ONLY: ld => local_decomposition,  &
                                          gd => global_decomposition
  USE mo_control,                   ONLY: nn, lnwp, lcolumn, ldebugio, nprocio
  USE mo_time_conversion,           ONLY: TC_get
  USE mo_time_control,              ONLY: get_forecast_hours,  &
                                          start_date,          & ! reference for time axis
                                          get_date_components, & ! split date into components
                                          next_date
  USE mo_linked_list,               ONLY: list_element,       & ! output stream entry
                                          t_stream,           & ! output stream data type
                                          memory_info,        & ! stream entry data type
                                          ZIP, SZIP
  USE mo_memory_base,               ONLY: maxstr, ostreams, nstreams,     & ! output streams
                                          nofiles, cfilenames,            & ! subjob interface
                                          GRIB, NETCDF, NETCDF2, NETCDF4, & ! allowed file types
                                          GAUSSIAN, SPECTRAL, LAND          ! grid representations
  USE mo_mpi,                       ONLY: p_pe, p_io, p_bcast
  USE mo_transpose,                 ONLY: reorder
  USE mo_tr_gather,                 ONLY: gather_field, gather_spectral
  USE mo_filename,                  ONLY: compose_filenames,     &
                                          standard_output_file,  &
                                          find_next_free_unit,   &
                                          out_expname, out_filetype, out_ztype
  USE mo_tracdef,                   ONLY: ntrac, & ! number of tracers actually defined
                                          trlist   ! tracer info variable

  USE mo_math_constants,            ONLY: pi
  USE mo_control,                   ONLY: nvclev, nlon, ngl, nhgl, nn, lnmi, &
                                          lnudge, vct
  USE mo_gaussgrid,                 ONLY: gl_gmu
  USE mo_version,                   ONLY: nversion

  USE mo_jsbach_comm_to_echam5mods, ONLY: mask, domain_mask, nland, domain_nland
  USE mo_column,                    ONLY: dlat,dlon
#ifdef HAVE_CDIPIO
  USE mo_echam_yaxt,                ONLY: gp => gp_gdeco, sp => sp_sdeco
  USE yaxt, ONLY: xt_idxlist_get_num_indices, xt_is_null
#endif

  IMPLICIT NONE

  INTEGER :: gaussianID, spectralID, taxisIDa, taxisIDr
  INTEGER :: belowsurID, surfaceID, hybridID, hybrid_hID
  INTEGER :: tilesID, soillevID, rootzonesID, canopyID
  INTEGER :: cosp_levelsID !ZK #474
  INTEGER :: instID, modelID
  INTEGER :: local_tableID, nudging_tableID, tracer_tableID, &
             land_tableID, ldiag_tableID, ldiag2_tableID, veg_tableID, &
             nitrogen_tableID, yasso_tableID, chem_tableID, dist_tableID

  INTEGER ::  ksec1(43)

  INTEGER, PARAMETER :: local_table     = 128 !  local code table
  INTEGER, PARAMETER :: nudging_table   = 129 !  nudging code table
  INTEGER, PARAMETER :: tracer_table    = 131 !  tracer code table
  INTEGER, PARAMETER :: land_table      = 180 !  land surface code table
  INTEGER, PARAMETER :: ldiag_table     = 180 !  code table for land diagnostics
  INTEGER, PARAMETER :: veg_table       = 181 !  land vegetation code table
  INTEGER, PARAMETER :: nitrogen_table  = 182 !  land nitrogen variabel table
  INTEGER, PARAMETER :: yasso_table     = 184 !  yasso carbon code table
  INTEGER, PARAMETER :: ldiag2_table    = 185 !  land nitrogen variabel table
  INTEGER, PARAMETER :: dist_table      = 186 !  disturbance code table
  INTEGER, PARAMETER :: chem_table      = 199 !  chemie code table
  INTEGER, PARAMETER :: center_id       = 252 !  identification of centre
  INTEGER            :: model_id              !  model identification
  INTEGER, PARAMETER :: grid_type       = 255 !  grid definition
  INTEGER, PARAMETER :: nflag           = 128 !  flag(see code table 1)
  INTEGER            :: code_parameter        !  data field code
                                              !  (see code table 2 )

  INCLUDE 'cdi.inc'
#ifdef HAVE_CDIPIO
  INCLUDE 'cdipio.inc'
#endif


  ! reference time of data

  INTEGER, PUBLIC    :: year                  !  year of century
  INTEGER, PUBLIC    :: month                 !  month
  INTEGER, PUBLIC    :: day                   !  day
  INTEGER, PUBLIC    :: hour                  !  hour
  INTEGER, PUBLIC    :: minute                !  minute
  INTEGER, PUBLIC    :: second                !  second

  INTEGER            :: time_unit       =   0 ! unit of time range
                                              ! (see code table 4)
  INTEGER            :: time_p1         =   0 ! time range 1
  INTEGER            :: time_p2         =   0 ! time range 2
  INTEGER            :: range_flag      =  10 ! time range flag
                                              ! (see code table 5)
  INTEGER, PUBLIC    :: century               ! century
  INTEGER, PARAMETER :: subcenter_id    =   1 ! subcenter

  INTEGER, PUBLIC    :: forecast_hours  =   0 ! number of forecast hours in
                                              ! NWP mode

  LOGICAL :: vlistsInitialised(maxstr) = .FALSE.
  LOGICAL :: p_call_cdi

  TYPE manage_file
   CHARACTER(LEN=8)   :: suffix   = ' '
   CHARACTER(LEN=128) :: filename = ' '
   INTEGER            :: filetype = -1
   INTEGER            :: ztype    = -1
   INTEGER            :: vlistID  = CDI_UNDEFID
   INTEGER            :: fileID   = CDI_UNDEFID
  END type manage_file

  TYPE(manage_file), ALLOCATABLE :: post(:)
  INTEGER                     :: npl ! actual number of output files

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE set_output_time

    IF (lnwp) THEN
       CALL get_date_components(start_date,year,month,day,hour,minute,second)
       forecast_hours = get_forecast_hours()
    ELSE
       CALL get_date_components(next_date,year,month,day,hour,minute,second)
       forecast_hours = 0
    END IF
    century = year/100+1
    year    = MOD(year,100)

  END SUBROUTINE set_output_time

!------------------------------------------------------------------------------
  SUBROUTINE close_output_streams

    !-------------------------------------------------------------------
    ! Loop over all the output files and close the associated files, set
    ! unit numbers to CDI_UNDEFID
    !-------------------------------------------------------------------
    CHARACTER(*), PARAMETER :: &
      label = 'Debug:I/O:mo_output::close_output_streams'
    INTEGER :: i

    IF (ldebugio) CALL message(label, '')
    nofiles = npl
    !---------------------------
    ! loop over all output files
    !---------------------------
    DO i = 1, npl
      IF (p_call_cdi) THEN
        CALL streamClose(post(i)%fileID)
      ENDIF
      post(i)%fileID = CDI_UNDEFID

      ! store subjob command
      cfilenames(i) = 'OSTREAM' // &
           TRIM(post(i)%suffix) // '=' &
           // TRIM(post(i)%filename)

    END DO

    ostreams(1:nstreams)%fileID = CDI_UNDEFID

  END SUBROUTINE close_output_streams

!------------------------------------------------------------------------------
  SUBROUTINE backup_output_streams

    !---------------------------------------------------------------------
    ! make a copy of all output streams
    !---------------------------------------------------------------------
    INTEGER            :: i, ilenc
    CHARACTER(len=512) :: ycopy

    ! the command is machine dependent, now only for linux testet
    CHARACTER(len=512), PARAMETER :: ycopy_cmd = 'cp'     ! backup command

    INTEGER, EXTERNAL :: util_system

    IF (p_pe == p_io) THEN

      outfile_loop: DO i = 1, npl

        IF (post(i)%fileID /= CDI_UNDEFID) THEN
          ycopy = TRIM(ycopy_cmd) // ' ' &
               // TRIM(post(i)%filename) // ' ' &
               // TRIM(post(i)%filename) // '.bak'
          ilenc = MIN(LEN_TRIM(ycopy),512)
          WRITE(message_text,*) 'backup: ',TRIM(post(i)%filename),' <',TRIM(ycopy),'>'
          CALL message('backup_output_streams',message_text)
          IF (util_system(ycopy(:ilenc)) /= 0) &
               CALL message('backup_output_streams','copy failed')
!!!               CALL finish('backup_output_streams','copy failed')
        END IF

      END DO outfile_loop

    END IF

  END SUBROUTINE backup_output_streams

!------------------------------------------------------------------------------
  SUBROUTINE open_output_streams

  ! Loop over all the output streams and open the associated files. Set
  ! unit numbers (file IDs) for all streams associated with a file.

    INTEGER                                  :: i ,j, ia  ! loop indices
    CHARACTER (LEN(standard_output_file)+4)  :: base  !
    CHARACTER (LEN(standard_output_file)+16) :: file  !
    INTEGER                                  :: used(255) ! used codes
    INTEGER                                  :: iunit ! codefile unit
    INTEGER                                  :: status
    INTEGER                                  :: bcast_buffer(3)

    IF (ldebugio) THEN
      WRITE (message_text, '(i0,a)') p_pe, ': START open_output_streams'
      IF (ldebugio) CALL message('Debug:I/O:mo_output: ', TRIM(message_text))
    ENDIF

    ! Derive base of filename
    ! filename is: standard_output_file[+forecast_hours][suffix][.nc]
    ! output streams with the same 'suffix' use the same file

    CALL compose_filenames
    IF (lnwp) THEN
      CALL set_output_time
      IF (forecast_hours > 744) THEN
        CALL message ('','NWP mode makes no sense for this time range.')
        CALL message ('','Please change to the climate mode.')
        CALL finish('open_output_streams','Run terminated.')
      END IF
    ENDIF
    base = standard_output_file
    !------------------------------
    ! print output streams (header)
    !------------------------------
    CALL message('',separator)
    CALL message('','')
    CALL message('','Open output streams:')
    CALL message('','')
    CALL message('','base filename '//base(1:index(base,'/',.TRUE.)-1))
    CALL message('','')
    WRITE(message_text,'(t5,a,t38,a,t54,a,t64,a,t70,a5)') &
         'file', 'stream', 'fileID', 'lpost', 'lrest'
    CALL message('',message_text)
    CALL message('','')
    !--------------------------------------------------------------------
    ! Loop over all the output streams and open the associated files. Set
    ! file IDs for all streams associated with a file.
    !--------------------------------------------------------------------
    outfile_loop: DO i = 1, npl
      !------------------------------------------------
      ! Derive filename
      ! Skip if column model is running and GRIB output
      !------------------------------------------------
      SELECT CASE (post(i)%filetype)
      CASE (GRIB)
        IF (lcolumn) CYCLE
        file = TRIM(base)//post(i)%suffix
        !---------------------
        ! open code table file
        !---------------------
        IF (p_pe == p_io) THEN
          iunit = find_next_free_unit (80, 100)
          OPEN (iunit,file=TRIM(file)//'.codes')
        ENDIF
      CASE (NETCDF,NETCDF2,NETCDF4)
        file = TRIM(base)//TRIM(post(i)%suffix)//'.nc'
      CASE default
        CALL finish('open_output_streams','unknown filetype')
      END SELECT

      ! Open output files
      IF (ldebugio) CALL message ('open_output_streams', 'Open '//TRIM(file))
      IF (p_call_cdi) THEN
        SELECT CASE (post(i)%filetype)
        CASE (GRIB)
          post(i)%fileID = streamOpenWrite(file, FILETYPE_GRB)
          IF ( post(i)%ztype == SZIP ) THEN
            CALL streamDefCompType(post(i)%fileID, COMPRESS_SZIP)
            CALL streamDefCompLevel(post(i)%fileID, 0)
          END IF
        CASE (NETCDF)
          post(i)%fileID = streamOpenWrite(file, FILETYPE_NC)
        CASE (NETCDF2)
          post(i)%fileID = streamOpenWrite(file, FILETYPE_NC2)
        CASE (NETCDF4)
          post(i)%fileID = streamOpenWrite(file, FILETYPE_NC4)
          IF ( post(i)%ztype == ZIP ) THEN
            CALL streamDefCompType(post(i)%fileID, COMPRESS_ZIP)
            CALL streamDefCompLevel(post(i)%fileID, 1)
          END IF
        CASE default
          CALL finish('open_output_streams','unknown filetype')
        END SELECT

        IF ( post(i)%fileID == CDI_UNDEFID ) THEN
          WRITE(message_text,'(a)') cdiStringError(post(i)%fileID)
          CALL message('',message_text)
          CALL finish ('open_output_streams', 'Open failed on '//TRIM(file))
        END IF
        post(i)%filename = file(1:16)

        !--------------------------------------------------------------
        ! due to performance reasons switch off fill mode (netCDF only)
        !--------------------------------------------------------------
        ! CALL nf(nf_set_fill(post(i)%fileID, nf_nofill, old_mode))

        !------------------------------------------
        ! define variable list associated with file
        !------------------------------------------
        CALL streamDefVlist(post(i)%fileID, post(i)%vlistID)

      ENDIF ! p_call_cdi

      IF (nprocio == 0) THEN
         CALL p_bcast (post(i)%fileID, p_io)
      ENDIF

      !----------------------------------------------------------------
      ! loop over all output streams and set file IDs of all associated
      ! output streams
      !----------------------------------------------------------------
      used = 0 ! initialize check list for test_codes
      stream_loop: DO j = 1, nstreams
        !--------------------------------------------------------
        ! Skip if stream does not correspond to the current file,
        ! or is not marked for output.
        !--------------------------------------------------------
        IF (.NOT. ostreams(j)%lpost .OR. &
            ostreams(j)%post_suf /= post(i)%suffix)   CYCLE
        IF (ostreams(j)%filetype /= post(i)%filetype)   &
          CALL finish('mo_output:open_output_streams',              &
                      'different file types for same file')
        !----------------------------
        ! check for valid grib codes
        ! print code table
        !----------------------------
        CALL test_codes

        ! set file ID and initialize timestep
        ostreams(j)%fileID = post(i)%fileID
        ostreams(j)%timestep = 0

        IF (p_pe==p_io) THEN
          WRITE(message_text,'(t5,a32,t38,a,t54,i9,t66,l1,t72,l1)') &
            file(index(file,'/',.TRUE.)+1:), &
            ostreams(j)%name, ostreams(j)%fileID, &
            ostreams(j)%lpost, ostreams(j)%lrerun
          CALL message('',message_text)
        ENDIF

        ostreams(j)%filename = file(1:16)

      END DO stream_loop


      IF (p_pe==p_io) THEN
        !----------------------
        ! close code-table file
        !----------------------
        IF (post(i)%filetype == GRIB) THEN
          IF (SUM(used) > 0) THEN
            CLOSE (iunit)
          ELSE
            CLOSE (iunit, status='DELETE')
          ENDIF
        ENDIF
      ENDIF
    END DO outfile_loop

    CALL message('','')

    IF (ldebugio) THEN
      WRITE (message_text, '(i0,a)') p_pe, ': END open_output_streams'
      IF (ldebugio) CALL message('Debug:I/O:mo_output: ', TRIM(message_text))
    ENDIF

  CONTAINS
    !
    ! check grib codes, write code table
    !
    SUBROUTINE test_codes
      TYPE (list_element) ,POINTER :: next
      INTEGER :: newcode, minl(1)
      !
      ! loop over stream entries
      !
      next => ostreams(j)%first_list_element
      DO
        IF (.NOT.ASSOCIATED(next)) EXIT
        !
        ! check if fields can be written
        !
        IF (next%field%info%lpost) THEN
          !
          ! check for NETCDF output
          !
          IF (ostreams(j)%filetype == NETCDF .OR.    &
              ostreams(j)%filetype == NETCDF2 .OR.   &
              ostreams(j)%filetype == NETCDF4) THEN
            SELECT CASE (next%field%info%repr)
            CASE (GAUSSIAN)
              SELECT CASE (next%field%info%ndim)
              CASE (2,3,4)
              CASE default
                next%field%info%lpost = .FALSE.
              END SELECT
            CASE (LAND)
              SELECT CASE (next%field%info%ndim)
              CASE (1,2,3)
              CASE default
                next%field%info%lpost = .FALSE.
              END SELECT
            CASE (SPECTRAL)
              IF (lcolumn) next%field%info%lpost = .FALSE.
            CASE default
              next%field%info%lpost = .FALSE.
            END SELECT
          ENDIF
          !
          ! check for GRIB output
          !
          IF (ostreams(j)%filetype == GRIB) THEN
            !
            ! check dimensions
            !
            SELECT CASE (next%field%info%repr)
            CASE (GAUSSIAN)
              SELECT CASE (next%field%info%ndim)
              CASE (2,3)
              CASE default
                next%field%info%lpost = .FALSE.
              END SELECT
            CASE (LAND)
              SELECT CASE (next%field%info%ndim)
              CASE (1,2)
              CASE default
                next%field%info%lpost = .FALSE.
              END SELECT
            END SELECT
            !
            ! check codes for GRIB output
            !
            IF (next%field%info%gribcode == 0 &
            .OR.next%field%info%gribcode >255 ) THEN
              !
              ! no valid code
              !
              CALL message('','No GRIB code for '//TRIM(next%field%info%name))
            ELSE IF (next%field%info%gribcode < 0) THEN
              !
              ! automatic code determination
              !
              minl    = MINLOC(used)
!>>SF #412 (fix for GRIB output corruption)
              !orig newcode = minl(1)
              newcode = minl(1) + SUM(used(minl(1):))
!<<SF #412
              IF (used(newcode)==0) THEN
                next%field%info%gribcode = newcode
              ELSE
                  CALL message('','No GRIB code for '//TRIM(next%field%info%name))
                  next%field%info%gribcode = 0
              ENDIF
            ENDIF
            !
            ! check for codes used twice
            !
            newcode = next%field%info%gribcode
            IF (newcode > 0) THEN
              IF (used(newcode)/=0) THEN
                CALL message('','GRIB code used twice for '//TRIM(next%field%info%name))
                next%field%info%gribcode = 0
              ELSE
                used (newcode) = used (newcode) + 1
              ENDIF
            ENDIF
            !
            ! write code table
            !
            IF (p_pe==p_io) THEN  ! must be 'i' here!
              IF (next%field%info%gribcode /= 0) THEN
                WRITE(iunit,'(i4,i4,1x,a,1x,f7.2,1x,f7.2,1x,a," [",a,"]")') &
                  next%field%info%gribcode,      &
                  next%field%info%klev,          &
                  next%field%info%name,          &
                  0., 1.,                     &
                  TRIM(next%field%info%longname),&
                  TRIM(next%field%info%units)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        next => next%next_list_element
      END DO
    END SUBROUTINE test_codes

  END SUBROUTINE open_output_streams
  !------------------------------------------------------------------------------
  SUBROUTINE addStreamToVlist(stream, vlistID)

    TYPE (t_stream) ,INTENT(in) :: stream ! output stream
    INTEGER,INTENT(inout) :: vlistID

    TYPE (list_element) ,POINTER    :: le
    TYPE (list_element) ,TARGET     :: first
    TYPE (memory_info)  ,POINTER    :: info

    INTEGER                :: varID, gridID, zaxisID
    INTEGER                :: n
    INTEGER                :: prec          ! float prec. (default)
    CHARACTER(len=5)       :: axis
    CHARACTER(len=32)      :: grid_type

    IF (ldebugio) THEN
      WRITE (message_text, '(a,i0,a)') 'pe ', p_pe, 'START addStreamToVlist'
      IF (ldebugio) CALL message('Debug:I/O:mo_output: ', TRIM(message_text))
    ENDIF

    first%next_list_element => stream%first_list_element
    !---------------------------------------
    ! 1st loop, define additional dimensions
    !---------------------------------------
    le => first
    DO ! loop over elements in linked list
      le => le%next_list_element
      IF (.NOT.ASSOCIATED(le)) EXIT
      info => le%field%info
      IF (.NOT. info%lpost)      CYCLE ! skip if lpost flag not set
      IF (.NOT.(info%repr == GAUSSIAN .OR. info%repr == LAND .OR. &
                info%repr == SPECTRAL)) CYCLE
      !-----------------------------------------------------------
      ! Only a standard 2D field (nlon,ngl) and a
      ! standard 3D field (nlon,...,ngl)
      ! and spectral representations are allowed for netcdf/GRIB
      ! output. For all other fields, info%lpost is set to .FALSE.
      ! surface pressure is written in any case.
      !-----------------------------------------------------------
      n = info%ndim
      IF (info%repr == SPECTRAL .AND. info%gdim(n) == 1) n=n-1
      !----------------------------------------
      ! write message for grid types not written
      !-----------------------------------------
      IF (.NOT.info%lpost) THEN
        WRITE(message_text,'(a,a,a,4i0)') &
             '   ',TRIM(info%name), ' is non-standard: info%gdim(:) = ' &
             ,info%gdim(1), info%gdim(2), info%gdim(3), info%gdim(4)
        CALL message('addStreamToVlist',message_text)
        CYCLE
      ENDIF
    END DO
    !------------------------------------------
    ! 2nd loop, define variables and attributes
    !------------------------------------------
    le => first
    DO ! loop over elements in linked list
      le => le%next_list_element
      IF (.NOT.ASSOCIATED(le)) EXIT
      info => le%field%info
      IF (.NOT. info%lpost)      CYCLE ! skip if lpost flag not set

      n = info%ndim
      SELECT CASE (info%repr)
      CASE (GAUSSIAN)
        info%gridID = gaussianID
        grid_type = 'gaussian'
        IF ((n==4)) THEN
          CALL finish('addStreamToVlist', '5-dimensional arrays unsupported!')
        ELSE IF ((n==3)) THEN
          !------------------------------
          ! 3d data, transpose dimensions
          !------------------------------
          axis       = 'tzyx'
        ELSE
          !-------------
          ! regular data
          !-------------
          axis      = 'tyx'
        ENDIF
      CASE (LAND)
        info%gridID = gaussianID
        grid_type = 'gaussian'
        IF ((n==4)) THEN
          CALL finish('addStreamToVlist', '5-dimensional arrays unsupported!')
        ENDIF
      CASE (SPECTRAL)
        info%gridID = spectralID
        IF (info%gdim(1) == 1) THEN
          n=2
          axis      = 't--'
        ELSE
          axis      = 'tz--'
        ENDIF
        grid_type = 'spectral, triangular truncation'
      CASE default
        CYCLE
      END SELECT

      IF ( info%levelindx == SURFACE ) THEN
        info%zaxisID = surfaceID
      ELSE IF ( info%levelindx == HYBRID ) THEN
        info%zaxisID = hybridID
      ELSE IF ( info%levelindx == HYBRID_H ) THEN
        info%zaxisID = hybrid_hID
      ELSE IF ( info%levelindx == BELOWSUR ) THEN
        info%zaxisID = belowsurID
      ELSE IF ( info%levelindx == TILES ) THEN
        info%zaxisID = tilesID
      ELSE IF ( info%levelindx == SOILLEV ) THEN
        info%zaxisID = soillevID
      ELSE IF ( info%levelindx == ROOTZONES ) THEN
        info%zaxisID = rootzonesID
      ELSE IF ( info%levelindx == CANOPY ) THEN
        info%zaxisID = canopyID
!>>ZK #474
      ELSE IF ( info%levelindx == COSP_LEVELS ) THEN
        info%zaxisID = cosp_levelsID 
!<<ZK #474
      ELSE
        WRITE(message_text,*) 'info%levelindx = ', info%levelindx
        CALL message('',message_text)
        CALL finish('addStreamToVlist', 'unsupported info%levelindx!')
      END IF

      gridID  = info%gridID
      zaxisID = info%zaxisID

      IF ( gridID  == -1 ) THEN
        WRITE(message_text,*) 'GRID undefined for ', info%name
        CALL message('',message_text)
        CALL finish('addStreamToVlist', 'GRID definition missing')
      END IF
      IF ( zaxisID == -1 ) THEN
        WRITE(message_text,*) 'ZAXIS undefined for ', info%name
        CALL message('',message_text)
        CALL finish('addStreamToVlist', 'ZAXIS definition missing')
      END IF

      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE)
      info%IO_var_stid = varID ! store for later use

      prec = DATATYPE_FLT32
      IF ( stream%filetype == GRIB ) THEN
        prec =  info%gribbits
      ELSE
        IF (info%gribbits > 32) prec = DATATYPE_FLT64
      END IF

      CALL vlistDefVarDatatype(vlistID, varID, prec)
      CALL vlistDefVarName(vlistID, varID, info%name)

      IF (info%longname /= '')  CALL vlistDefVarLongname(vlistID, varID, info%longname)
      IF (info%units /= '')     CALL vlistDefVarUnits(vlistID, varID, info%units)
      IF (info%gribcode > 0)    CALL vlistDefVarCode(vlistID, varID, info%gribcode)
      IF (info%lmiss)           CALL vlistDefVarMissval(vlistID, varID, info%missval)

      IF (info%gribtable > 0) THEN
        SELECT CASE (info%gribtable)
        CASE(local_table)
          CALL vlistDefVarTable(vlistID, varID, local_tableID)
        CASE(nudging_table)
          CALL vlistDefVarTable(vlistID, varID, nudging_tableID)
        CASE(tracer_table)
          CALL vlistDefVarTable(vlistID, varID, tracer_tableID)
        CASE(land_table)
          CALL vlistDefVarTable(vlistID, varID, land_tableID)
        CASE(ldiag2_table)
          CALL vlistDefVarTable(vlistID, varID, ldiag2_tableID)
        CASE(nitrogen_table)
          CALL vlistDefVarTable(vlistID, varID, nitrogen_tableID)
        CASE(yasso_table)
          CALL vlistDefVarTable(vlistID, varID, yasso_tableID)
        CASE(dist_table)
           CALL vlistDefVarTable(vlistID, varID, dist_tableID)
        CASE(veg_table)
          CALL vlistDefVarTable(vlistID, varID, veg_tableID)
        CASE(chem_table)
          CALL vlistDefVarTable(vlistID, varID, chem_tableID)
        END SELECT
      END IF

      CALL vlistDefVarInstitut(vlistID, varID, instID)
      CALL vlistDefVarModel(vlistID, varID, modelID)

      !------------------------------------------
      ! print tracer attributes into netcdf file:
      !------------------------------------------
      IF (info%tracidx > 0) CALL write_tracer_header(info%tracidx, vlistID, varID)

    END DO


    IF (ldebugio) THEN
      WRITE (message_text, '(a,i0,a)') 'pe ', p_pe, 'END addStreamToVlist'
      IF (ldebugio) CALL message('Debug:I/O:mo_output: ', TRIM(message_text))
    ENDIF

  END SUBROUTINE addStreamToVlist

  SUBROUTINE write_tracer_header(tracidx, vlistID, varID)

    INTEGER ,INTENT(in) :: tracidx             ! tracer index
    INTEGER ,INTENT(in) :: vlistID, varID
    INTEGER             :: status

    status = vlistDefAttInt(vlistID, varID, 'index',      DATATYPE_INT32, 1, tracidx)
    status = vlistDefAttFlt(vlistID, varID, 'molar_mass', DATATYPE_FLT64, 1, trlist%ti(tracidx)%moleweight)
!++mgs
!!    status = vlistDefAttFlt(vlistID, varID, 'Henry',      DATATYPE_FLT64, 1, trlist%ti(tracidx)%henry)
!!    status = vlistDefAttFlt(vlistID, varID, 'dryreac',    DATATYPE_FLT64, 1, trlist%ti(tracidx)%dryreac)
!--mgs
    status = vlistDefAttInt(vlistID, varID, 'ndrydep',    DATATYPE_INT32, 1, trlist%ti(tracidx)%ndrydep)
    status = vlistDefAttInt(vlistID, varID, 'ntran',      DATATYPE_INT32, 1, trlist%ti(tracidx)%ntran)
    status = vlistDefAttInt(vlistID, varID, 'nvdiff',     DATATYPE_INT32, 1, trlist%ti(tracidx)%nvdiff)
    status = vlistDefAttInt(vlistID, varID, 'nconv',      DATATYPE_INT32, 1, trlist%ti(tracidx)%nconv)
    status = vlistDefAttInt(vlistID, varID, 'nwetdep',    DATATYPE_INT32, 1, trlist%ti(tracidx)%nwetdep)
    status = vlistDefAttInt(vlistID, varID, 'nsoluble',   DATATYPE_INT32, 1, trlist%ti(tracidx)%nsoluble)
    status = vlistDefAttInt(vlistID, varID, 'nphase',     DATATYPE_INT32, 1, trlist%ti(tracidx)%nphase)
    status = vlistDefAttInt(vlistID, varID, 'mode',       DATATYPE_INT32, 1, trlist%ti(tracidx)%mode)

  END SUBROUTINE write_tracer_header

!------------------------------------------------------------------------------
  SUBROUTINE init_output

    TYPE (io_dim) ,POINTER :: pdim
    INTEGER :: year, month, day, hour, minute, second ! date/time variables
    INTEGER :: idate, itime
    CHARACTER (32)  :: ymodel

    INTEGER :: i, j, ia ! loop indices
    INTEGER :: status

    IF (ldebugio) THEN
      WRITE (message_text, '(a,i0,a)') 'pe ', p_pe, 'START init_output'
      IF (ldebugio) CALL message('Debug:I/O:mo_output: ', TRIM(message_text))
    ENDIF

    IF (nprocio == 0 .AND. p_pe /= p_io) THEN
      p_call_cdi = .FALSE.
    ELSE
      p_call_cdi = .TRUE.
    ENDIF

    ! get all output file suffixes
    ALLOCATE(post(nstreams))
    CALL build_list_of_output_files(ostreams, post)

    ! set model version
    model_id = nversion + 10
    WRITE(ymodel,'(a,f3.1)') 'ECHAM', REAL(nversion,dp)/10.0_dp
    CALL setup_tables

    ! set fixed values in GRIB block 1, common to all GRIB sets

    ! to set standard output to analysis mode lnwp must be .false. and
    ! lanalysis .true.. Later with working adiabatic NMI this must be
    ! changed to range_flag = 1. range_flag = 0 is default for output
    ! of the OI or 3DVAR analysis, this are given external. In case a
    ! 4DVAR is running inside ECHAM this value must be adjusted to
    ! range_flag = 0, as it is now.
    ! For usual climate mode runs range_flag is set to 10, which means
    ! The time given in the GRIB header in section 1 is the valid time
    ! and time_p1 and time_p2 do not mean anything.

    IF (lnmi .AND. .NOT.lnudge) range_flag = 0
    ksec1(18) = range_flag

    ! create horizontal grids
    CALL setup_horizontal_grids

    ! create vertical axes
    IF (nvclev > 127) CALL message('GRIB', 'VCT not defined for more than 126 levels!')
    CALL setup_zaxes

    ! create time axes
    CALL setup_taxes

    ! create lists of variables
    CALL setup_vlists

    ! add streams to lists of variables (i.e. define variables)
    DO i = 1, nstreams
      IF(.NOT. ostreams(i)%lpost) CYCLE
      CALL addStreamToVlist(ostreams(i), ostreams(i)%vlistID)
    ENDDO

    IF (ldebugio) THEN
      WRITE (message_text, '(a,i0,a)') 'pe ', p_pe, 'END init_output'
      IF (ldebugio) CALL message('Debug:I/O:mo_output: ', TRIM(message_text))
    ENDIF

  CONTAINS
    SUBROUTINE build_list_of_output_files(ostreams, post)

      TYPE(t_stream),    INTENT(INOUT)  :: ostreams(:)
      TYPE(manage_file), INTENT(OUT) :: post(:)

      INTEGER :: i, j

      ostreams%first = .FALSE.

      npl = 0
      stream_loop: DO i = 1, nstreams
        IF (ostreams(i)%lpost) THEN ! account for output streams only
          suffix_loop: DO j = 1, npl
            IF (ostreams(i)%post_suf == post(j)%suffix) THEN
              ! Found a match in list of suffixes so go to the next stream
              CYCLE stream_loop
            ENDIF
          ENDDO suffix_loop

          ! No match found so add suffix to the list
          ostreams(i)%first = .TRUE.
          npl = npl + 1
          post(npl)%suffix   = ostreams(i)%post_suf
          post(npl)%filetype = ostreams(i)%filetype
          post(npl)%ztype    = ostreams(i)%ztype
          post(npl)%vlistID  = CDI_UNDEFID
          post(npl)%fileID   = CDI_UNDEFID
        ENDIF
      ENDDO stream_loop

    END SUBROUTINE build_list_of_output_files

    SUBROUTINE setup_tables()

      instID  = institutDef(center_id, subcenter_id, "MPIMET", "Max-Planck-Institute for Meteorology")
      modelID = modelDef(instID,  model_id, ymodel)

      ! define tables
      local_tableID   = tableDef(modelID, local_table,   "echam6")
      nudging_tableID = tableDef(modelID, nudging_table, "echam6")
      tracer_tableID  = tableDef(modelID, tracer_table,  "echam6")
      land_tableID    = tableDef(modelID, land_table,    "echam6")
      ldiag_tableID   = tableDef(modelID, ldiag_table,   "echam6")
      ldiag2_tableID  = tableDef(modelID, ldiag2_table,  "echam6")
      veg_tableID     = tableDef(modelID, veg_table,     "echam6")
      nitrogen_tableID= tableDef(modelID, nitrogen_table,"echam6")
      yasso_tableID   = tableDef(modelID, yasso_table,   "echam6")
      dist_tableID     = tableDef(modelID, dist_table,   "echam6")
      chem_tableID    = tableDef(modelID, chem_table,    "echam6")

    END SUBROUTINE setup_tables

    SUBROUTINE setup_horizontal_grids()
      REAL(dp), ALLOCATABLE :: xvals(:), yvals(:)

      IF(lcolumn) THEN
        gaussianID = gridCreate(GRID_GAUSSIAN, 1)
        CALL gridDefXsize(gaussianID, 1)
        CALL gridDefYsize(gaussianID, 1)
        ALLOCATE(xvals(1))
        ALLOCATE(yvals(1))
        xvals = REAL(dlon, dp)
        yvals = REAL(dlat, dp)
        CALL gridDefXvals(gaussianID,xvals)
        CALL gridDefYvals(gaussianID,yvals)
        DEALLOCATE(xvals)
        DEALLOCATE(yvals)
      ELSE
        gaussianID = gridCreate(GRID_GAUSSIAN, nlon*ngl)
        CALL gridDefXsize(gaussianID, nlon)
        CALL gridDefYsize(gaussianID, ngl)
        ALLOCATE(xvals(nlon))
        ALLOCATE(yvals(ngl))
        DO i = 1, nlon
          xvals(i) = (i-1)*360.0_dp/REAL(nlon,dp)
        END DO
        DO i = 1, nhgl
          yvals(i) = ASIN(gl_gmu(i))*180.0_dp/pi
          yvals(ngl-i+1) = -yvals(i)
        END DO
        CALL gridDefXvals(gaussianID, xvals)
        CALL gridDefYvals(gaussianID, yvals)
        DEALLOCATE(xvals)
        DEALLOCATE(yvals)
      END IF

      spectralID = gridCreate(GRID_SPECTRAL, (nn+1)*(nn+2))
      CALL gridDefTrunc(spectralID, nn)
      CALL gridDefComplexPacking(spectralID, 1)
    END SUBROUTINE setup_horizontal_grids

    SUBROUTINE setup_zaxes()
      REAL(dp), ALLOCATABLE :: levels(:)

      surfaceID  = zaxisCreate(ZAXIS_SURFACE, 1)
#ifndef STANDALONE
      hybridID   = zaxisCreate(ZAXIS_HYBRID, nvclev-1)
      hybrid_hID = zaxisCreate(ZAXIS_HYBRID_HALF, nvclev)
#endif

      ALLOCATE(levels(1))
      levels(1) = 0.0_dp
      CALL zaxisDefLevels(surfaceID, levels)
      DEALLOCATE(levels)

      belowsurID = -1
      IF ( IO_get_varindx("belowsurface") >= 0 ) THEN
        pdim => IO_dim_ids(BELOWSUR)
        belowsurID = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, pdim%dim_len)
        CALL zaxisDefName(belowsurID, pdim%dim_name)
        IF (ASSOCIATED (pdim%value)) CALL zaxisDefLevels(belowsurID, pdim%value)
      END IF

      tilesID = -1
      IF ( IO_get_varindx("tiles") >= 0 ) THEN
        pdim => IO_dim_ids(TILES)
        tilesID = zaxisCreate(ZAXIS_GENERIC, pdim%dim_len)
        CALL zaxisDefName(tilesID, pdim%dim_name)
        CALL zaxisDefLongname(tilesID, pdim%longname)
        CALL zaxisDefUnits(tilesID, pdim%units)
        CALL zaxisDefLtype(tilesID, pdim%levtyp)
        IF (ASSOCIATED (pdim%value)) CALL zaxisDefLevels(tilesID, pdim%value)
      END IF

      soillevID = -1
      IF ( IO_get_varindx("soil_layer") >= 0 ) THEN
        pdim => IO_dim_ids(SOILLEV)
        soillevID = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, pdim%dim_len)
        CALL zaxisDefName(soillevID, pdim%dim_name)
        CALL zaxisDefLongname(soillevID, pdim%longname)
        CALL zaxisDefUnits(soillevID, pdim%units)
        CALL zaxisDefLtype(soillevID, pdim%levtyp)
        IF (ASSOCIATED (pdim%value)) CALL zaxisDefLevels(soillevID, pdim%value)
      END IF

      rootzonesID = -1
      IF ( IO_get_varindx("root_zone") >= 0 ) THEN
        pdim => IO_dim_ids(ROOTZONES)
        rootzonesID = zaxisCreate(ZAXIS_GENERIC, pdim%dim_len)
        CALL zaxisDefName(rootzonesID, pdim%dim_name)
        CALL zaxisDefLongname(rootzonesID, pdim%longname)
        CALL zaxisDefUnits(rootzonesID, pdim%units)
        CALL zaxisDefLtype(rootzonesID, pdim%levtyp)
        IF (ASSOCIATED (pdim%value)) CALL zaxisDefLevels(rootzonesID, pdim%value)
      END IF

      canopyID = -1
      IF ( IO_get_varindx("canopy_layer") >= 0 ) THEN
        pdim => IO_dim_ids(CANOPY)
        canopyID = zaxisCreate(ZAXIS_GENERIC, pdim%dim_len)
        CALL zaxisDefName(canopyID, pdim%dim_name)
        CALL zaxisDefLongname(canopyID, pdim%longname)
        CALL zaxisDefUnits(canopyID, pdim%units)
        CALL zaxisDefLtype(canopyID, pdim%levtyp)
        IF (ASSOCIATED (pdim%value)) CALL zaxisDefLevels(canopyID, pdim%value)
      END IF

!>>ZK #474
    cosp_levelsID = -1
    IF ( IO_get_varindx("cosplev") >= 0 ) THEN
      pdim => IO_dim_ids(COSP_LEVELS)
      cosp_levelsID = zaxisCreate(ZAXIS_ALTITUDE, pdim%dim_len)
      CALL zaxisDefName(cosp_levelsID, pdim%dim_name)
      CALL zaxisDefLongname(cosp_levelsID, pdim%longname)
      CALL zaxisDefUnits(cosp_levelsID, pdim%units)
      CALL zaxisDefLtype(cosp_levelsID, pdim%levtyp)
      IF (ASSOCIATED (pdim%value)) CALL zaxisDefLevels(cosp_levelsID, pdim%value)
    END IF
!<<ZK #474
#ifndef STANDALONE
      ALLOCATE(levels(nvclev))
      DO i = 1, nvclev
        levels(i) = REAL(i,dp)
      END DO
      CALL zaxisDefLevels(hybridID, levels)
      CALL zaxisDefLevels(hybrid_hID, levels)
      DEALLOCATE(levels)

      CALL zaxisDefVct(hybridID,   2*nvclev, vct(1:2*nvclev))
      CALL zaxisDefVct(hybrid_hID, 2*nvclev, vct(1:2*nvclev))
#endif
    END SUBROUTINE setup_zaxes

    SUBROUTINE setup_taxes()
      taxisIDa = taxisCreate(TAXIS_ABSOLUTE)
      taxisIDr = taxisCreate(TAXIS_RELATIVE)

      CALL taxisDefCalendar(taxisIDa, CALENDAR_PROLEPTIC)
      CALL taxisDefCalendar(taxisIDr, CALENDAR_PROLEPTIC)

      CALL get_date_components(start_date,year,month,day,hour,minute,second)
      idate = ABS(year)*10000+month*100+day
      IF ( year < 0 ) idate = -idate
      itime = hour*10000+minute*100+second

      CALL taxisDefRdate(taxisIDr, idate)
      CALL taxisDefRtime(taxisIDr, itime)
      CALL taxisDefTunit(taxisIDr, TUNIT_DAY)

    END SUBROUTINE setup_taxes

    SUBROUTINE setup_vlists()

      ! create lists of variables
      suffix_loop: DO i = 1, npl
        stream_loop: DO j = 1, nstreams
          IF (ostreams(j)%post_suf == post(i)%suffix) THEN
            IF (post(i)%vlistID == CDI_UNDEFID) THEN
              ostreams(j)%vlistID = vlistCreate()
              post(i)%vlistID = ostreams(j)%vlistID
              post(i)%filetype = ostreams(j)%filetype
            ELSE
              ostreams(j)%vlistID = post(i)%vlistID
              IF (post(i)%filetype /= ostreams(j)%filetype) THEN
                CALL finish('mo_output::init_output','wrong mix of filetypes...')
              ENDIF
            ENDIF
          ENDIF
        ENDDO stream_loop
      ENDDO suffix_loop

      ! define time axes for lists of variables
      DO i = 1, npl
        IF (post(i)%filetype == GRIB) THEN
          CALL vlistDefTaxis(post(i)%vlistID, taxisIDa)
        ELSE
          CALL vlistDefTaxis(post(i)%vlistID, taxisIDr)
        ENDIF
      ENDDO

      ! define global attributes for lists of variables
      DO i = 1, npl
        status = vlistDefAttTxt(post(i)%vlistID, CDI_GLOBAL, 'title', &
             LEN(TRIM(out_expname)), TRIM(out_expname))
        DO ia = 1, SIZE(global_att)
          IF (global_att(ia)%name /= '') &
            status = vlistDefAttTxt(post(i)%vlistID, CDI_GLOBAL, &
                 global_att(ia)%name, LEN(TRIM(global_att(ia)%text)),&
                 TRIM(global_att(ia)%text) )
        END DO
        status = vlistDefAttInt(post(i)%vlistID, CDI_GLOBAL, &
             'truncation', DATATYPE_INT32, 1, nn)
      ENDDO

    END SUBROUTINE setup_vlists

  END SUBROUTINE init_output

!------------------------------------------------------------------------------
  SUBROUTINE out_streams
    !-------------------------------------------------
    ! loop over all output files
    ! write if flag lpost is set and
    !       if output time (flag L_PUTDATA) is reached
    !-------------------------------------------------

    LOGICAL :: time_written, write_info, data_written
    INTEGER :: idate, itime, ts
    INTEGER :: i, j

    write_info = .TRUE.

    !-----------------------------------------------
    ! loop over all output files
    ! pick up first stream associated with each file
    !-----------------------------------------------
    DO i = 1, npl
      time_written = .FALSE.
      !-----------------------------------------------
      ! loop over all streams associated with the file
      !-----------------------------------------------
      DO j = i, nstreams
        IF (ostreams(j)%filetype == post(i)%filetype .AND. &
            ostreams(j)%fileID   == post(i)%fileID   ) THEN
          !---------------------------
          ! check condition for output
          !---------------------------
          IF(l_putdata (ostreams(j)%post_idx) &
             .AND.      ostreams(j)%lpost     &
             .AND.      ostreams(j)%lpout     ) THEN
            !
            IF (write_info) THEN
              SELECT CASE (out_filetype)
              CASE (NETCDF)
                CALL write_date(next_date,'Write netCDF output for : ')
              CASE (NETCDF2)
                CALL write_date(next_date,'Write netCDF2 output for : ')
              CASE (NETCDF4)
                IF ( out_ztype == ZIP ) THEN
                  CALL write_date(next_date,'Write netCDF4/ZIP output for : ')
                ELSE
                  CALL write_date(next_date,'Write netCDF4 output for : ')
                END IF
              CASE (GRIB)
                IF ( out_ztype == SZIP ) THEN
                  CALL write_date(next_date,'Write GRIB1/SZIP output for : ')
                ELSE
                  CALL write_date(next_date,'Write GRIB1 output for : ')
                END IF
              END SELECT
              write_info = .FALSE.
            ENDIF

            !------------------------------------
            ! increase time slice in output files
            !------------------------------------
            IF (.NOT. time_written) THEN
              CALL write_time_to_stream(ostreams(j), idate, itime)
              time_written = .TRUE.
            ENDIF
          ENDIF ! stream scheduled and enabled for output
        END IF ! stream associated with current file
      END DO ! streams
    END DO ! files

    ! Time step have to be defined for ALL output files
    ! BEFORE data write
    data_written = .FALSE.
    DO i = 1, npl
      !-----------------------------------------------
      ! loop over all streams associated with the file
      !-----------------------------------------------
      DO j = i, nstreams
        IF (ostreams(j)%filetype == post(i)%filetype .AND. &
            ostreams(j)%fileID   == post(i)%fileID   ) THEN
          !---------------------------
          ! check condition for output
          !---------------------------
          IF(l_putdata (ostreams(j)%post_idx) &
             .AND.      ostreams(j)%lpost     &
             .AND.      ostreams(j)%lpout     ) THEN
            !
            !----------------
            ! write variables
            !----------------
            CALL out_stream (ostreams(j))
            data_written = .TRUE.
            IF ( nprocio == 0 .AND. p_call_cdi ) &
                  CALL streamSync(ostreams(j)%fileID)

          ENDIF ! stream scheduled and enabled for output
        END IF ! stream associated with current file
      END DO ! streams
    END DO ! files

#ifdef HAVE_CDIPIO
    IF (nprocio > 0 .AND. data_written) THEN
      !------------------------------------------------
      ! communicate model time to collecting I/O PE
      ! expose MPI memory windows to collecting I/O PEs
      !------------------------------------------------
      CALL pioWriteTimestep
    ENDIF
#endif

    DO i = 1, nstreams
      IF (l_putdata(ostreams(i)%post_idx)) THEN
        ! Do reset for all (output and not output) streams
        CALL reset_stream(ostreams(i))
        ! increment streams timestep
        ostreams(i)%timestep = ostreams(i)%timestep + 1
      END IF
    END DO ! streams


  end SUBROUTINE out_streams
!------------------------------------------------------------------------------
  !>
  !! Control postprocessing of output fields.
  !!
  !! Generic routine for all streams
  !! If a stream is not output (lpost=.FALSE. or lpout=.FALSE.) its variables that
  !! have laccu=.TRUE. should still be reset to the reset value because these
  !! accumulated variables might be used elsewhere.
  !!
  SUBROUTINE reset_stream (stream)
    !
    ! Argument
    !
    TYPE (t_stream) ,INTENT(in) :: stream
    !
    ! variables of derived type used in linked list
    !
    TYPE (list_element) ,POINTER :: element

    IF(ldebugio) THEN
      CALL message('Debug:I/O:mo_output::reset_stream: stream', TRIM(stream%name))
    END IF

    element => stream%first_list_element
    DO
      IF(.NOT.ASSOCIATED(element)) EXIT

      !-------------------------------------------------
      ! reset field if accumulation or reset flag is set
      !-------------------------------------------------
      IF (element%field%info%laccu .OR. element%field%info%reset /= 0._dp) THEN

        IF(ldebugio) THEN
          CALL message('Debug:I/O:mo_output::reset_stream: element', &
                       TRIM(element%field%info%name))
        END IF

        element%field%ptr = element%field%info%reset
      ENDIF

      element => element%next_list_element
    END DO

  END SUBROUTINE reset_stream
!------------------------------------------------------------------------------
  SUBROUTINE out_stream (stream)
  !
  ! Description:
  !
  ! Control postprocessing of output fields.
  ! Generic routine for all streams
  !
    !
    ! Argument
    !
    TYPE (t_stream) ,INTENT(IN) :: stream
    !
    !  Local scalars:
    !
    REAL(dp)          :: interval_factor, output_factor
    INTEGER           :: gridtype, interval_seconds
    CHARACTER(len=16) :: yname
    INTEGER           :: level_idx ! index to level definition table
    !
    ! variables of derived type used in linked list
    !
    TYPE (memory_info)  ,POINTER :: info
    TYPE (list_element) ,POINTER :: element
    TYPE (list_element) ,TARGET  :: start
    !
    !  Local arrays:
    !
    REAL(dp),POINTER     :: ptr4d (:,:,:,:) ! field distributed over processors
    REAL(dp),POINTER     :: z4d   (:,:,:,:) ! field gathered on I/O processor
                                            ! or locally reordered
    REAL(dp),ALLOCATABLE :: ptr4dx(:,:,:,:) ! for LAND: domain (lon, tiles, lat, 1)
    INTEGER              :: dimx(4)         ! for LAND: dimensions of ptr4dx
    INTEGER              :: nlandgc         !number of land grid cells (local or global)
    INTEGER              :: i

    IF(ldebugio) THEN
      CALL message('Debug:I/O:mo_output::out_stream: stream', TRIM(stream%name))
    END IF
    !
    ! Definition of grib blocks
    !
    CALL set_output_time

    IF (lnwp) THEN
      ksec1(15) = time_unit
      ksec1(16) = forecast_hours
      ksec1(18) = range_flag
    END IF

    ! Get accumulation interval
    interval_seconds = get_interval_seconds(ev_putdata(stream%post_idx))
    interval_factor = 1._dp
    IF(interval_seconds > 0) interval_factor = 1._dp/REAL(interval_seconds, dp)

    IF(ldebugio) THEN
      WRITE (message_text, '(i0)') interval_seconds
      CALL message('Debug:I/O:mo_output::out_stream: interval', &
                   TRIM(message_text))
      WRITE (message_text, *) interval_factor
      CALL message('Debug:I/O:mo_output::out_stream: 1/interval', &
                   TRIM(message_text))
    END IF

    !
    ! Loop over all fields in linked list
    !
    element => start
    element%next_list_element => stream%first_list_element
    DO
      element => element%next_list_element
      IF(.NOT.ASSOCIATED(element)) EXIT
      !-----------------------------------------------------
      ! retrieve information from actual linked list element
      !-----------------------------------------------------
      info           => element%field%info
      ptr4d          => element%field%ptr(:,:,:,:)
      yname          = info%name(1:16)
      code_parameter = info%gribcode
      gridtype       = info%repr
      level_idx      = info%levelindx
      !------------------
      ! skip this field ?
      !------------------
      IF (.NOT. info%lpost                                    .OR. &
           yname == ' '                                        .OR. &
           (stream%filetype == GRIB .AND. code_parameter <= 0) .OR. &
           (level_idx <1 .OR. level_idx > SIZE(IO_dim_ids))         &
         ) THEN
         ! Nothing else has to be done
         CYCLE
      END IF

      IF(ldebugio) THEN
        CALL message('Debug:I/O:mo_output::out_stream: element', &
                     TRIM(info%name) // ' ' // MERGE('accu', '    ', info%laccu))
      END IF
      !------------------------------------------
      ! rescale field if accumulation flag is set
      !------------------------------------------
      output_factor = MERGE(interval_factor, 1._dp, info%laccu)

      !----------------------------------------------------
      ! Allocate temporary global array on output processor
      ! Gather field from other processors
      !----------------------------------------------------
      NULLIFY(z4d)
      IF (.NOT.lcolumn) THEN
        SELECT CASE (gridtype)
!===============================================================================
        CASE (GAUSSIAN)
          IF (nprocio == 0) THEN  ! root IO
            IF (p_pe==p_io) ALLOCATE (z4d (info%gdim(1),info%gdim(2), &
                 &                         info%gdim(3),info%gdim(4)) )
            CALL gather_field(z4d, info%gdim, ptr4d*output_factor)
          ELSE ! parallel IO
            ALLOCATE (z4d (info%dim(1), info%dim(2), info%dim(3), info%dim(4)))
            ! FIXME: reodering (should be eliminated later)
            SELECT CASE (info%ndim)
              CASE (2)
                CALL reorder (z4d(:,:,1,1), ptr4d(:,:,1,1)*output_factor)
              CASE (3)
                CALL reorder (z4d(:,:,:,1), ptr4d(:,:,:,1)*output_factor)
              CASE (4)
                CALL reorder (z4d(:,:,:,:), ptr4d(:,:,:,:)*output_factor)
              CASE DEFAULT
                CALL finish('out_stream','unsupported number of dimensions')
            END SELECT
          ENDIF
!===============================================================================
        CASE (LAND)
          IF (info%ndim > 2) &
             CALL finish('out_stream','Only 2 dimensions in LAND streams allowed')

          ALLOCATE(ptr4dx(SIZE(domain_mask,1),info%gdim(2),SIZE(domain_mask,2),info%gdim(3)))
          DO i = 1, info%gdim(2) ! tiles or soil levels
            IF (info%lmiss) THEN
              ptr4dx(:,i,:,1) = UNPACK(ptr4d(:,i,1,1)*output_factor, MASK=domain_mask, &
                   FIELD=info%missval)
            ELSE
              ptr4dx(:,i,:,1) = UNPACK(ptr4d(:,i,1,1)*output_factor, MASK=domain_mask, &
                   FIELD=0._dp)
            END IF
          ENDDO

          IF (nprocio == 0) THEN ! root IO
            IF (p_pe == p_io) &
                ALLOCATE(z4d(SIZE(mask,1),info%gdim(2),SIZE(mask,2),info%gdim(3)))
            dimx = (/ SIZE(mask,1), info%gdim(2), &
                      SIZE(mask,2), info%gdim(3) /)
            CALL gather_field(z4d, dimx, ptr4dx)
          ELSE ! parallel IO
            ALLOCATE(z4d(ld%nglon,info%gdim(2),ld%nglat,info%gdim(3)))
            CALL reorder (z4d(:,:,:,1), ptr4dx(:,:,:,1))
          ENDIF
          DEALLOCATE(ptr4dx)
!===============================================================================
        CASE (SPECTRAL)
          IF (nprocio == 0) THEN
            IF (p_pe==p_io) ALLOCATE (z4d(info%gdim(1),info%gdim(2), &
                                           info%gdim(3),info%gdim(4)) )
            CALL gather_spectral(z4d, info%gdim, ptr4d*output_factor)
          ELSE
            ALLOCATE (z4d (info%dim(1),info%dim(2), info%dim(3),info%dim(4)) )
            ! FIXME: check reordering
            z4d(:,:,:,:) = ptr4d(:,:,:,:)*output_factor
          ENDIF
        CASE default
          CALL finish('out_stream','unknown grid type')
        END SELECT
      ENDIF ! .not. lcolumn
      !----------------------
      ! Write data
      !----------------------
      nlandgc = MERGE(nland, domain_nland, nprocio == 0)
      IF (p_call_cdi) THEN
        IF (lcolumn) THEN
          CALL write_var (info, stream, ptr4d(:,:,:,:)*output_factor, nlandgc)
        ELSE
          CALL write_var (info, stream, z4d(:,:,:,:), nlandgc)
        ENDIF
      ENDIF

      !-----------------------------------
      ! Deallocate temporary global arrays
      !-----------------------------------
      IF (ASSOCIATED(z4d)) DEALLOCATE (z4d)
    END DO

  END SUBROUTINE out_stream

!------------------------------------------------------------------------------
  SUBROUTINE cleanup_output
    !
    ! Deallocate module variables
    !
    INTEGER :: i

    DO i=1,nstreams
      IF (ostreams(i)%vlistID /= CDI_UNDEFID) THEN
        CALL vlistDestroy(ostreams(i)%vlistID)
      ENDIF

      WHERE (ostreams%filetype == ostreams(i)%filetype  &
             .AND. ostreams%fileID   == ostreams(i)%fileID)   &
             ostreams%vlistID = CDI_UNDEFID
    ENDDO

    CALL gridDestroy(gaussianID)
    CALL gridDestroy(spectralID)

    CALL zaxisDestroy(surfaceID)
    CALL zaxisDestroy(hybridID)
    CALL zaxisDestroy(hybrid_hID)
    IF ( belowsurID  >= 0 ) CALL zaxisDestroy(belowsurID)
    IF ( tilesID     >= 0 ) CALL zaxisDestroy(tilesID)
    IF ( soillevID   >= 0 ) CALL zaxisDestroy(soillevID)
    IF ( rootzonesID >= 0 ) CALL zaxisDestroy(rootzonesID)
    IF ( canopyID    >= 0 ) CALL zaxisDestroy(canopyID)

    DEALLOCATE(post)

  END SUBROUTINE cleanup_output

!------------------------------------------------------------------------------
  SUBROUTINE write_var (info, stream, xzy, nland)

    !>>dod
    USE mo_control,   ONLY: ldebugio
    USE mo_exception, ONLY: message, message_text, em_debug
    !<<dod

    TYPE (memory_info) ,INTENT(in) :: info         ! field description
    TYPE (t_stream)    ,INTENT(in) :: stream       ! output stream description
    REAL(dp)           ,INTENT(in) :: xzy(:,:,:,:) ! output field
    INTEGER            ,INTENT(in) :: nland        ! number of land points

    REAL(dp)    ,ALLOCATABLE :: xyz(:,:,:,:) ! transposed field
    INTEGER :: jy, jz                    ! indices used for transposition
    INTEGER :: fileID                    ! File ID
    INTEGER :: varID                     ! Variable ID
    INTEGER :: n                         ! rank of field to write
    INTEGER :: nmiss

    fileID  = stream%fileID
    n     = info%ndim
    nmiss = 0
    !-----------------------------------------------------------
    ! variable ID may be invalid for no_cycle > 0 (rerun cycle).
    !  request ID from NetCDF file in this case.
    !-----------------------------------------------------------
    varID   = info%IO_var_stid

    !----------------------------------------
    ! write 3D,2D Gaussian or spectral fields
    !----------------------------------------
    SELECT CASE (info%repr)
    CASE (GAUSSIAN, LAND)

      ! 4d array
      IF ( SIZE(xzy,4) .GT. 1 ) THEN
        !-----------------------------------------------
        ! The array xzy is sorted xzy(lon,lev,?,lat) but
        ! we require (lon,lat,lev,?).
        !-----------------------------------------------

        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,4),SIZE(xzy,2),SIZE(xzy,3)))
        FORALL (jy=1:SIZE(xzy,4),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,:)=xzy(:,jz,:,jy) ! switch lat and lev
        END FORALL

        CALL finish('4-dimensional arrays unsupported for output!')

        DEALLOCATE (xyz)

      ! 3d array
      ELSE IF ( SIZE(xzy,3) .GT. 1 ) THEN
        !-------------------------------------------------------------------
        ! The array xzy is sorted xzy(lon,lev,lat) - ONLY FOR single netcdf
        ! output. Parallel netcdf data hasn't been 'gather'd, so is still
        ! scrambled and needs to be 'reorder'd and unpacked -  but the
        ! COARDS convention for netcdf requires (lon,lat,lev). Applies
        ! to /all/ fields here.
        !-------------------------------------------------------------------

        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,3),SIZE(xzy,2),1))
        FORALL (jy=1:SIZE(xzy,3),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,1)=xzy(:,jz,jy,1) ! switch lat and lev
        END FORALL

        IF ( info%repr == LAND ) nmiss = (SIZE(xzy,1)*SIZE(xzy,3) - nland) * SIZE(xzy,2)

        !>>dod
        IF (ldebugio) THEN
           message_text = 'Write variable: stream '//TRIM(stream%name)//' variable '//TRIM(info%name)
           CALL message('mo_output', message_text, level=em_debug)
        END IF
        !<<dod

        IF (nprocio == 0) THEN
          CALL StreamWriteVar(fileID, varID, xyz(:,:,:,1), nmiss)
        ELSE
#ifdef HAVE_CDIPIO
          !PRINT *, 'pe ', p_pe, ', varid=', varid, 'size(xyz) = ', size(xyz), &
          !     'xt_idxlist_get_num_indices = ', xt_idxlist_get_num_indices(gp(SIZE(xyz,3))%idxlist), &
          !     'dim=', info%dim, SIZE(xyz,3)
          CALL streamWriteVarPart(fileID, varID, xyz(:,:,:,1), nmiss, &
          gp(SIZE(xyz,3))%idxlist)
#endif
        ENDIF
        DEALLOCATE (xyz)

      ! 2d array
      ELSE
        IF ( info%repr == LAND ) nmiss = SIZE(xzy,1)*SIZE(xzy,2) - nland

        !>>dod
        IF (ldebugio) THEN
           message_text = 'Write variable: stream '//TRIM(stream%name)//' variable '//TRIM(info%name)
           CALL message('mo_output', message_text, level=em_debug)
        END IF
        !<<dod

        IF (nprocio == 0) THEN
          CALL streamWriteVar(fileID, varID, xzy(:,:,1,1), nmiss)
        ELSE
#ifdef HAVE_CDIPIO
          CALL streamWriteVarPart(fileID, varID, xzy(:,:,1,1), nmiss, &
          gp(1)%idxlist)
#endif
        ENDIF
      END IF

!===============================================================================
    CASE (SPECTRAL)
      IF (info%gdim(1) == 1) n = n - 1  ! number of non-singleton dimensions

      SELECT CASE (n)
      CASE (3)

        ALLOCATE (xyz (SIZE(xzy,2),SIZE(xzy,3),SIZE(xzy,1),1))
        FORALL (jz=1:SIZE(xzy,1))
          xyz(:,:,jz,1)=xzy(jz,:,:,1) ! switch lat and lev
        END FORALL

        IF (nprocio == 0) THEN
          CALL streamWriteVar(fileID, varID, xyz(:,:,:,1), 0)
        ELSE
#ifdef HAVE_CDIPIO
          CALL streamWriteVarPart(fileID, varID, xyz(:,:,:,1), nmiss, &
          sp(SIZE(xyz,3))%idxlist)
#endif
        ENDIF
        DEALLOCATE (xyz)

      CASE (2)

        IF (nprocio == 0) THEN
          CALL streamWriteVar(fileID, varID, xzy(1,:,:,1), 0)
        ELSE
#ifdef HAVE_CDIPIO
          CALL streamWriteVarPart(fileID, varID, xzy(1,:,:,1), nmiss, &
          sp(1)%idxlist)
#endif
        ENDIF

      END SELECT

    END SELECT

  END SUBROUTINE write_var

!------------------------------------------------------------------------------
  SUBROUTINE write_time_to_stream (stream, idate, itime)

    TYPE (t_stream), INTENT(IN)  :: stream    ! output stream description
    INTEGER,         INTENT(OUT) :: idate, itime

    INTEGER :: start_day, start_sec, present_day, present_sec
    INTEGER :: year, month, day, hour, minute, second ! date/time variables
    INTEGER :: iret

    !-----------------------------------------
    ! convert time_days format into 2 integers
    !-----------------------------------------
    CALL TC_get(start_date,    start_day,   start_sec)
    CALL TC_get(next_date,   present_day, present_sec)

    !--------------------------------
    ! write alternative time yyyymmdd
    !--------------------------------
    CALL get_date_components(next_date,year,month,day,hour,minute,second)

    idate = cdiEncodeDate(year, month, day)
    itime = cdiEncodeTime(hour, minute, second)

    CALL taxisDefVdate(taxisIDa, idate)
    CALL taxisDefVtime(taxisIDa, itime)

    CALL taxisDefVdate(taxisIDr, idate)
    CALL taxisDefVtime(taxisIDr, itime)


    IF (ldebugio) THEN
      WRITE (message_text, '(i0,a,i0)') p_pe, ': define timestep ', &
             stream%timestep
      IF (ldebugio) CALL message('Debug:I/O:mo_output: ', &
                    TRIM(message_text))
    ENDIF
    IF (p_call_cdi) THEN
      iret = streamDefTimestep(stream%fileID, stream%timestep)
    ENDIF

  END SUBROUTINE write_time_to_stream

END MODULE mo_output