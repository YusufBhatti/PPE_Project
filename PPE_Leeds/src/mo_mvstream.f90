!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_mvstream

  !--------------------------------------------------------------------------
  !
  ! Module that allows to open an extra redoubled stream for some tracers
  ! or elements of any other stream (3d or 2d grid point fields) 
  ! in order to calculate for example the mean values over a certain interval
  ! -----------------------------------------
  !
  ! Authors:
  !  J.S. Rast, MPI-Met, Hamburg, December 2003, original source
  !  P. Stier,  MPI-Met, Hamburg,          2004, modifications
  !  J.S. Rast, MPI-Met, Hamburg, June     2004, optimisation
  !--------------------------------------------------------------------------

  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_kind,                  ONLY: wp
  USE mo_linked_list,           ONLY: t_stream, memory_info, &
                                      GRIDPOINT, SPECTRAL, LAND, UNKNOWN
  USE mo_mpi,                   ONLY: p_parallel_io
  USE mo_memory_base,           ONLY: maxstr, get_stream, remove_stream_element, &
                                      default_output
  USE mo_memory_gl,             ONLY: xt
  USE mo_namelist,              ONLY: open_nml, position_nml, close_nml, &
                                      POSITIONED
  USE mo_time_control,          ONLY: delta_time, ev_putdata, putdata
  USE mo_time_event,            ONLY: io_time_event, &
                                      TIME_INC_MONTHS, TRIG_FIRST, &
                                      events_nonequal
  USE mo_tracdef,               ONLY: trlist
  USE mo_util_string,           ONLY: int2string
  USE mo_string_utls,           ONLY: st1_in_st2_idx
  USE mo_decomposition,         ONLY: ldc => local_decomposition

  IMPLICIT NONE

  PRIVATE

  PUBLIC                    :: init_mvstream, &
                               mvstream_update_cache, &
                               mvstream_accumulate

  TYPE(t_stream) :: dummy_stream  ! Dummy for reading defaults and sizes
  TYPE(memory_info) :: dummy_info ! Dummy for reading defaults and sizes

  INTEGER, PARAMETER        :: nam_max=500

  TYPE(io_time_event), PARAMETER :: &
       default_interval = io_time_event(1, TIME_INC_MONTHS, TRIG_FIRST, 0)

  ! @todo Check necessity of workaround for GNU Fortran bug
  ! in array member initialization of dummy_meta
  CHARACTER(*), PARAMETER   :: empty_name = REPEAT(' ', LEN(dummy_info%name))

  TYPE metastream
    CHARACTER(LEN(dummy_stream%name)) :: target = ''
    TYPE(io_time_event) :: interval = io_time_event(0, '', '', 0)
    CHARACTER(LEN(dummy_stream%post_suf)) :: filetag = ''
    CHARACTER(LEN(dummy_stream%name)) :: source = ''
    CHARACTER(LEN(dummy_info%name)) :: meannam(nam_max) = empty_name
    CHARACTER(LEN(dummy_info%name)) :: sqrmeannam(nam_max) = empty_name
  END TYPE metastream

  TYPE(metastream), SAVE :: dummy_meta ! Dummy for reading defaults and sizes

  TYPE stream_cache
    TYPE(t_stream), POINTER :: stream => NULL()
  END TYPE stream_cache

  LOGICAL                   :: lvarmean(nam_max, maxstr), &
                               lstddev(nam_max, maxstr), &
                               lgridpoint(nam_max, maxstr)
  INTEGER                   :: ivarrank(nam_max, maxstr)
  CHARACTER(LEN(dummy_stream%name)) :: mstreamname(maxstr)
  CHARACTER(LEN(dummy_info%name)) :: mvarname(nam_max, maxstr), &
                               ystddev
  INTEGER                   :: nmstream,nmvar(maxstr),nlist_elements(maxstr)
  TYPE(io_time_event), SAVE :: putmean_stream(maxstr)
  TYPE(metastream), SAVE    :: metastream_data(maxstr)
  TYPE(stream_cache), SAVE  :: stream_cache_data(maxstr)

  TYPE vptr
    REAL(wp), POINTER      :: ptr3d(:,:,:)
    REAL(wp), POINTER      :: ptr3d_s(:,:,:)
    REAL(wp), POINTER      :: ptr2d(:,:)
    REAL(wp), POINTER      :: ptr2d_s(:,:)
    REAL(wp), POINTER      :: ptr1d(:,:)
    REAL(wp), POINTER      :: ptr1d_s(:,:)
    LOGICAL                :: source_accu = .FALSE.
  END TYPE vptr

  TYPE tptr
    REAL(wp), POINTER      :: ptt3d(:,:,:)
    REAL(wp), POINTER      :: ptt2d(:,:)
    REAL(wp), POINTER      :: ptt1d(:,:)
  END TYPE tptr

  TYPE(vptr), SAVE          :: varptr(nam_max, maxstr)

  TYPE(tptr), SAVE          :: tarptr(nam_max, maxstr)

CONTAINS

  SUBROUTINE init_mvstream
    !------------------------------------------------------------------------
    ! initialize mean value stream
    !
    ! Authors:
    ! J.S. Rast, MPI-Met, Hamburg, December 2003, original source
    !------------------------------------------------------------------------
    USE mo_tracdef,             ONLY: trlist, ln, ll
    USE mo_tracer,              ONLY: get_tracer
    USE mo_util_string,         ONLY: separator
    USE mo_mpi,                 ONLY: p_bcast, p_parallel, p_io
    USE mo_time_control,        ONLY: p_bcast_event
    USE mo_linked_list,         ONLY: list_element
    USE mo_memory_base,         ONLY: new_stream, &
                                      add_stream_element, &
                                      add_stream_reference, &
                                      default_stream_setting, &
                                      get_stream_element_info, &
                                      get_stream_element, AUTO

    LOGICAL                        :: lmean(nam_max), lstd_all, lmean_all
    CHARACTER(len=256)             :: ierrchar
    INTEGER                        :: ierr, kstream, knam, kdx, nunit, iunit
    INTEGER                        :: irank(nam_max), istd(nam_max), ilist_cnt
    TYPE(memory_info)              :: se_info
    TYPE(t_stream), POINTER        :: mv, stream
    TYPE (list_element) ,POINTER   :: this_list_element
    LOGICAL :: mvctl_exists
    TYPE(metastream) :: meta

    ! Namelist variables - make sure that this is consistent with include below
    TYPE(io_time_event), SAVE :: putmean
    CHARACTER(LEN(se_info%name)) :: meannam(nam_max)
    INTEGER :: stddev(nam_max)

    INCLUDE 'mvctl.inc'

    IF (p_parallel_io) THEN

      ! Read simplified configuration if available
      ! Side effects: updates nmstream and mstreamname
      CALL read_metastream_data

      ! Check given list of source streams against list of stream names

      DO kstream = 1, nmstream
        ! if number of streams greater than maxstr, error message by new_stream

        ! Check if source streams exists
        ! and if number of variables does not exceed buffer size.
        IF ( mstreamname(kstream) /= 'tracer' ) THEN

          ! Default handling of source streams

          ! proof that stream exists
          CALL get_stream ( mv, TRIM(mstreamname(kstream)) )
          IF (.NOT.ASSOCIATED(mv)) THEN
            CALL finish ('mo_mvstream   :', &
                 'stream '//TRIM(mstreamname(kstream))//' does not exist')
          END IF
          ! proof that stream contains <= nam_max variables
          IF (mv%list_elements > nam_max) THEN 
            WRITE(ierrchar,*) mv%list_elements
            CALL finish ('mo_mvstream   , init_mvstream ', &
                 'number of stream elements of stream '// &
                 TRIM(mstreamname(kstream))// &
                 'exceeds '//TRIM(ierrchar)// &
                 ', may not get all mean value variables. '// &
                 'Please, increase nam_max in mo_mvstream   .')
          END IF

        ELSE ! mstreamname(kstream) /= 'tracer'

          ! Special handling of tracer stream

          ! proof that number of tracers <= nam_max
          IF (trlist%ntrac > nam_max) THEN
            WRITE(ierrchar,*) trlist%ntrac
            CALL finish ('mo_mvstream   , init_mvstream ', &
                 'number of tracer exceeds '//TRIM(ierrchar)// &
                 ', may not get all mean value tracers. '// &
                 'Please, increase nam_max in mo_mvstream   .')
          END IF

        END IF ! ELSE mstreamname(kstream) /= 'tracer'

      END DO ! kstream = 1, nmstream

      IF (nmstream > 0) THEN
        CALL message('','')
        CALL message('', separator)
        CALL message('','mean values calculated for following variables:')
        CALL message('','')
      ENDIF
      mvarname(:,:)      = ''
      lgridpoint = .FALSE.
      DO kstream = 1, nmstream

        meta = metastream_data(kstream)

        ! set default values
        putmean            = meta%interval
        meannam(:)         = meta%meannam
        lmean(:)           = .FALSE.
        irank(:)           = 0
        stddev(:)          = 0
        istd(:)            = 0
        nmvar(kstream)     = 0
        lmean_all          = .FALSE.
        lstd_all           = .FALSE.

        ! If old style mvctl namelist file exists, read config data.

        INQUIRE(file=TRIM(mstreamname(kstream))//'.nml', exist=mvctl_exists)
        IF(mvctl_exists) THEN

          ! read time intervals and variable names for each stream
          nunit = open_nml(TRIM(mstreamname(kstream))//'.nml')
          iunit = position_nml ('MVCTL', nunit, status=ierr)
          SELECT CASE (ierr)
          CASE (POSITIONED)
            READ(iunit, mvctl)
          END SELECT
          CALL close_nml (nunit)

          ! all variables in stream (deprecated convention)
          IF (TRIM(meannam(1)) == 'all' .OR. TRIM(meannam(1)) == 'ALL' ) THEN
            meannam(1) = ''
          END IF

          ! New namelist overrides old style
          IF(meta%interval%unit /= dummy_meta%interval%unit) THEN
            CALL message('Warning', "any 'putmean' settings from '" // &
                 TRIM(mstreamname(kstream))//'.nml' // &
                 "' will be overridden by 'mvstreamctl'")
            putmean = meta%interval
          END IF
          IF(meta%meannam(1) /= dummy_meta%meannam(1)) THEN
            CALL message('Warning', "any 'meannam' and 'stddev' settings from '" // &
                 TRIM(mstreamname(kstream))//'.nml' // &
                 "' will be overridden by 'mvstreamctl'")
            meannam = meta%meannam
            stddev = 0 ! Ignored, as meannam lists might be inconsistent
          END IF

        END IF ! mvctl_exists

        ! Set dependent defaults
        IF(putmean%unit == dummy_meta%interval%unit) THEN
          IF(default_output) THEN
            putmean = default_interval
          ELSE
            putmean = putdata
          END IF
        END IF

        ! save putmean interval
        putmean_stream(kstream) = putmean

        ! empty list of mean variables takes all
        IF(meannam(1) == '' .OR. meannam(1) == '*') THEN
          meannam(1) = ''
          lmean_all = .TRUE.
        END IF

        ! mean of square for all variables in tracer stream
        IF (stddev(1) == -1 .OR. meta%sqrmeannam(1) == '*') THEN
          meta%sqrmeannam(1) = ''
          lstd_all=.TRUE.
        END IF

        ! count variables for each stream and assure that variables exist
        IF (TRIM(mstreamname(kstream)) == 'tracer') THEN

          ! Special handling for tracer stream.
          !
          ! Tracer stream is supposed to exist,
          ! and tracer variables do not need to be checked for representation.
          IF (lmean_all) THEN

            ! select all tracers
            mvarname(1:trlist%ntrac,kstream)=trlist%ti(1:trlist%ntrac)%fullname
            lmean(1:trlist%ntrac)  = .TRUE.
            nmvar(kstream)         = trlist%ntrac
            IF (lstd_all) THEN
              istd(1:trlist%ntrac)= 1
            ELSE
              istd(1:trlist%ntrac)= stddev(1:trlist%ntrac)
            END IF

          ELSE ! lmean_all

            ! individual tracers
            ilist_cnt=0
            DO knam = 1, trlist%ntrac
              IF (meannam(knam) /= '') THEN
                ilist_cnt=ilist_cnt+1
                CALL get_tracer (meannam(knam), idx=kdx, ierr = ierr)
                IF (ierr == 1)THEN
                  CALL finish ('mo_mvstream, init_mvstream','tracer ' &
                       //TRIM(meannam(knam))//' in stream '//TRIM(mstreamname(kstream))// &
                       ' is not defined')
                ELSE
                  lmean(kdx) = .TRUE.
                  mvarname(kdx,kstream) = meannam(knam)
                  nmvar(kstream) = nmvar(kstream) + 1
                  IF (lstd_all) THEN
                    istd(kdx) = 1
                  ELSE
                    istd(kdx) = stddev(nmvar(kstream))
                  END IF  !lstd_all
                END IF  ! tracer defined
              ELSE
                ! Empty name signals end of list
                EXIT
              END IF ! meannam(knam) /= ''
            END DO  ! loop over tracers

          END IF  ! ELSE lmean_all

        ELSE 

          ! General stream, not the tracer stream
          !
          ! Need to check whether stream name actually exists,
          ! for supported representation (gridpoint, spectral),
          ! and for accumulation variables (laccu)
          CALL get_stream ( mv, TRIM(mstreamname(kstream)) )
          IF (lmean_all) THEN

            ! select all variables
            IF (lstd_all) THEN
              istd(1:mv%list_elements) = 1
            ELSE
              istd(1:mv%list_elements) = stddev(1:mv%list_elements)
            END IF

            this_list_element => mv%first_list_element
            DO knam = 1, mv%list_elements
              meannam(knam) = TRIM(this_list_element%field%info%name)
              CALL get_stream_element_info (mv, TRIM(meannam(knam)), se_info)
              ! allow accumulation only for 3d, 2d or 1d
              ! gaussian grid point, spectral, or land fields
              IF (se_info%ndim == UNKNOWN) THEN
                CALL finish ('mo_mvstream, init_mvstream','stream element ' &
                     //TRIM(meannam(knam))//' in stream '//TRIM(mstreamname(kstream))// &
                     ' is not defined')
              END IF
              IF ( ANY(se_info%repr == (/ GRIDPOINT, SPECTRAL, LAND /)) ) THEN
                IF ( ANY(se_info%ndim == (/ 1, 2, 3 /)) ) THEN
                  ! variable type is supported
                  lmean(knam) = .TRUE.
                  mvarname(knam,kstream) = meannam(knam)
                  nmvar(kstream) = nmvar(kstream) + 1
                  irank(knam) = se_info%ndim
                  lgridpoint(knam, kstream) = se_info%repr == GRIDPOINT
                ELSE
                  istd(knam) = 0
                  CALL message ('mo_mvstream, init_mvstream','variable ' &
                       //meannam(knam)//' in stream '//TRIM(mstreamname(kstream))// &
                       ': neither 1d, 2d, nor 3d => skipped')
                END IF ! ndim
              ELSE
                istd(knam) = 0
                CALL message ('mo_mvstream, init_mvstream','variable ' &
                     //meannam(knam)//' in stream '//TRIM(mstreamname(kstream))// &
                     ': neither gridpoint, spectral, nor land => skipped')
              END IF ! repr
              this_list_element => this_list_element%next_list_element
            END DO  ! loop over all stream elements

          ELSE ! lmean_all

            ! individual variables
            ilist_cnt=0
            DO knam = 1, mv%list_elements
              IF (meannam(knam) /= '') THEN
                ilist_cnt=ilist_cnt+1
                CALL get_stream_element_info (mv, TRIM(meannam(knam)), se_info)
                ! allow accumulation only for 3d, 2d or 1d
                ! gaussian grid point, spectral, or land fields
                IF (se_info%ndim == UNKNOWN) THEN
                  CALL finish ('mo_mvstream, init_mvstream','stream element ' &
                       //TRIM(meannam(knam))//' in stream '//TRIM(mstreamname(kstream))// &
                       ' is not defined')
                END IF
                IF ( ANY(se_info%repr == (/ GRIDPOINT, SPECTRAL, LAND /)) ) THEN
                  IF ( ANY(se_info%ndim == (/ 1, 2, 3 /)) ) THEN
                    ! variable type is supported
                    lmean(knam) = .TRUE.
                    mvarname(knam,kstream) = meannam(knam)
                    nmvar(kstream) = nmvar(kstream) + 1
                    irank(knam) = se_info%ndim
                    lgridpoint(knam, kstream) = se_info%repr == GRIDPOINT
                    IF (lstd_all) THEN
                      istd(knam) = 1
                    ELSE
                      istd(knam) = stddev(ilist_cnt)
                    END IF
                  ELSE
                    CALL message ('mo_mvstream, init_mvstream','variable ' &
                         //TRIM(meannam(knam))//' in stream '//TRIM(mstreamname(kstream))// &
                         ': neither 1d, 2d, nor 3d => skipped')
                  END IF ! ndim
                ELSE  ! not a supported grid variable or variable has accumulation flag
                  CALL message ('mo_mvstream, init_mvstream','variable ' &
                       //TRIM(meannam(knam))//' in stream '//TRIM(mstreamname(kstream))// &
                       ': neither gridpoint, spectral, nor land => skipped')
                END IF ! repr
              END IF ! meannam(knam) /= ''
            END DO  ! loop over stream elements

          END IF ! ELSE lmean_all

        END IF  ! tracer stream or other stream

        ! Update stddev array according to sqrmeannam spec
        sqrmeannam_entries: DO knam = 1, nam_max

          IF(meta%sqrmeannam(knam) == '') EXIT

          kdx = st1_in_st2_idx(meta%sqrmeannam(knam), &
                               mvarname(:nmvar(kstream), kstream))
          IF(kdx > 0) THEN
            istd(kdx) = 1
          ELSE
            CALL message('mo_mvstream, init_mvstream','square mean for '// &
                 TRIM(meta%sqrmeannam(knam))//' in stream '// &
                 TRIM(mstreamname(kstream))// &
                 ': base variable skipped or non-existent => skipped')
          END IF

        END DO sqrmeannam_entries

        ! table of variables that require mean value calculation
        IF (nmvar(kstream) > 0) THEN
          CALL message('','')
          CALL message('','            stream '//TRIM(mstreamname(kstream))//':')
          message_text = ''
          DO knam = 1, nmvar(kstream)
            message_text = message_text//TRIM(meannam(knam))
          ENDDO
          CALL message('',message_text)
          message_text = ''
          DO knam = 1, nmvar(kstream)
            WRITE (ystddev,'(1x,i6)') stddev(knam)
            message_text = message_text//TRIM(ystddev)
          ENDDO
          CALL message('',message_text)
          CALL message('',separator)
        ENDIF
        lvarmean(:, kstream) = lmean(:)
        ivarrank(:, kstream) = irank(:)
        lstddev(:, kstream) = istd(:) == 1
      END DO ! kstream
      !!debug
      !     WRITE(0,*)
      !     WRITE(0,*) '******************** begin debug output mo_mvstream ***********************'
      !     WRITE(0,*) 'mo_mvstream   , init_mvstream:'
      !     WRITE(0,*) 'number of streams, nmstream=',nmstream
      !     WRITE(0,*) 'names of streams, mstreamname=', &
      !          (TRIM(mstreamname(kstream)),' ',kstream=1,nmstream)
      !     WRITE(0,*) 'number of variables in each stream:',nmvar(1:nmstream)
      !     WRITE(0,*) 'variables in each stream:'
      !     WRITE(0,*) ((TRIM(mvarname(knam,kstream)),' ',knam=1,nam_max),kstream=1,nmstream)
      !     WRITE(0,*) 'mean value flag in each stream', &
      !          ((lvarmean(knam,kstream),' ',knam=1,nam_max),kstream=1,nmstream)
      !     WRITE(0,*) 'variable rank in each stream', &
      !          ((ivarrank(knam,kstream),' ',knam=1,nam_max),kstream=1,nmstream)
      !     WRITE(0,*) 'stddev flag in each stream', &
      !          ((lstddev(knam,kstream),' ',knam=1,nam_max),kstream=1,nmstream)
      !     WRITE(0,*) 'put mean interval for each stream', putmean_stream(1:nmstream)
      !     WRITE(0,*) '******************** end debug output mo_mvstream *************************'
      !     WRITE(0,*)
      !!debug.
    END IF ! (p_parallel_io)
    ! send information on all processors
    IF (p_parallel) THEN
      CALL p_bcast(nmstream, p_io)
      CALL p_bcast(mstreamname, p_io)
      CALL p_bcast(nmvar, p_io)
      CALL p_bcast(metastream_data%target, p_io)
      CALL p_bcast(metastream_data%filetag, p_io)
      DO kstream = 1, nmstream
        CALL p_bcast(mvarname(:, kstream), p_io)
        CALL p_bcast(lvarmean(:, kstream), p_io)
        CALL p_bcast(ivarrank(:, kstream), p_io)
        CALL p_bcast(lstddev(:, kstream), p_io)
        CALL p_bcast(lgridpoint(:, kstream), p_io)
        CALL p_bcast_event(putmean_stream(kstream), p_io)
      END DO
    END IF
    ! open new streams
    DO kstream = 1, nmstream

      IF (nmvar(kstream) /= 0) THEN

        CALL new_stream (mv, metastream_data(kstream)%target, &
             post_suf='_'//metastream_data(kstream)%filetag, &
             interval=putmean_stream(kstream), &
             lpost=.TRUE., lcontnorest=.TRUE.)

        ! Store stream pointer for later use
        stream_cache_data(kstream)%stream => mv

        CALL default_stream_setting (mv, lrerun = .TRUE., contnorest = .TRUE., &
             lpost = .TRUE., laccu = .TRUE., code = AUTO)
        IF (TRIM(mstreamname(kstream)) == 'tracer') THEN
          DO knam = 1, trlist%ntrac
            IF (lvarmean(knam,kstream)) THEN
              CALL add_stream_element (mv, trlist%ti(knam)%fullname, varptr(knam,kstream)%ptr3d, &
                   units = trlist%ti(knam)%units, longname = trlist%ti(knam)%longname)
            END IF  ! lvarmean is true
            IF (lstddev(knam,kstream)) THEN
              IF (LEN_TRIM(trlist%ti(knam)%fullname) > ln - 2) THEN
                WRITE(ierrchar,*) ln
                CALL finish('mo_mvstream   , init_mvstream', 'tracer name '// &
                     TRIM(trlist%ti(knam)%fullname)//' exceeds '//ierrchar// &
                     '-2 alphanumeric characters, cannot define mean of squares tracer name.')
              END IF
              IF (LEN_TRIM(trlist%ti(knam)%units) > ln - 5) THEN
                WRITE(ierrchar,*) ln
                CALL finish('mo_mvstream   , init_mvstream', 'unit name '// &
                     TRIM(trlist%ti(knam)%units)//' of tracer '//TRIM(trlist%ti(knam)%fullname)// &
                     ' exceeds '//ierrchar// &
                     '-5 alphanumeric characters, cannot define mean of squares variable units.')
              END IF
              IF (LEN_TRIM(trlist%ti(knam)%longname) > ll - 10) THEN
                WRITE(ierrchar,*) ll
                CALL message('mo_mvstream, init_mvstream', 'longname '// &
                     TRIM(trlist%ti(knam)%longname)//' of tracer ' &
                     //TRIM(trlist%ti(knam)%fullname)//' exceeds '//ierrchar// &
                     '-10 alphanumeric characters, cannot define longname of square.')
              END IF
              CALL add_stream_element (mv, TRIM(trlist%ti(knam)%fullname)//'_s', varptr(knam,kstream)%ptr3d_s, &
                   units = '('//TRIM(trlist%ti(knam)%units)//')**2', &
                   longname = 'square of '//TRIM(trlist%ti(knam)%longname))
            END IF  ! calculation of mean of square
          END DO  ! tracer list
        ELSE  ! general stream
          CALL get_stream (stream, TRIM(mstreamname(kstream)))
          nlist_elements(kstream)=stream%list_elements
          DO knam = 1, nlist_elements(kstream)
            IF (lvarmean(knam, kstream)) THEN
              CALL get_stream_element_info (stream, TRIM(mvarname(knam, kstream)), &
                   se_info)

              IF (se_info%laccu) THEN
                ! If value is already accu'd, just reference it.
                ! This reference might later be deleted if the final output
                ! intervals do not match (see mvstream_update_cache).
                ! Remove internal accumulation marker,
                ! but set marker for external accumulation.
                ! Any stddev computation requested for these fields is disabled.
                lvarmean(knam, kstream) = .FALSE.
                varptr(knam, kstream)%source_accu = .TRUE.
                IF(lstddev(knam, kstream)) THEN
                  CALL message('mo_mvstream   , init_mvstream', 'field ' // &
                       TRIM(se_info%name) // ' is accumulated in source ' // &
                       'stream; disabling stddev')
                  lstddev(knam, kstream) = .FALSE.
                END IF
                CALL add_stream_reference(mv, se_info%name, mstreamname(kstream))
              ELSE IF (se_info%ndim == 1) THEN
                CALL add_stream_element (mv, TRIM(se_info%name), varptr(knam,kstream)%ptr1d, &
                     units = se_info%units, longname = se_info%longname, &
                     ldims = se_info%dim, gdims = se_info%gdim, &
                     klev = se_info%klev, repr = se_info%repr, &
                     table = se_info%gribtable, code = se_info%gribcode, &
                     bits = se_info%gribbits, leveltype = se_info%levelindx, &
                     lmiss = se_info%lmiss, missval = se_info%missval, &
                     lpost = se_info%lpost) ! lrerun must be true
                CALL get_stream_element (stream, TRIM(se_info%name), tarptr(knam,kstream)%ptt1d)
              ELSE IF (se_info%ndim == 2) THEN
                CALL add_stream_element (mv, TRIM(se_info%name), varptr(knam,kstream)%ptr2d, &
                     units = se_info%units, longname = se_info%longname, &
                     ldims = se_info%dim, gdims = se_info%gdim, &
                     klev = se_info%klev, repr = se_info%repr, &
                     table = se_info%gribtable, code = se_info%gribcode, &
                     bits = se_info%gribbits, leveltype = se_info%levelindx, &
                     lmiss = se_info%lmiss, missval = se_info%missval, &
                     lpost = se_info%lpost) ! lrerun must be true
                CALL get_stream_element (stream, TRIM(se_info%name), tarptr(knam,kstream)%ptt2d)
              ELSE IF (se_info%ndim == 3) THEN
                CALL add_stream_element (mv, TRIM(se_info%name), varptr(knam,kstream)%ptr3d, &
                     units = se_info%units, longname = se_info%longname, &
                     ldims = se_info%dim, gdims = se_info%gdim, &
                     klev = se_info%klev, repr = se_info%repr, &
                     table = se_info%gribtable, code = se_info%gribcode, &
                     bits = se_info%gribbits, leveltype = se_info%levelindx, &
                     lmiss = se_info%lmiss, missval = se_info%missval, &
                     lpost = se_info%lpost) ! lrerun must be true
                CALL get_stream_element (stream, TRIM(se_info%name), tarptr(knam,kstream)%ptt3d)
              END IF ! se_info%laccu | se_info%ndim in [1,2,3]
            END IF ! lvarmean
            IF (lstddev(knam,kstream)) THEN
              CALL get_stream_element_info (stream, TRIM(mvarname(knam, kstream)), &
                   se_info)
              IF(LEN_TRIM(se_info%name) > 30) THEN
                CALL finish('mo_mvstream   , init_mvstream', 'variable name '// &
                     TRIM(se_info%name)//' exceeds 30' // &
                     ' alphanumeric characters, cannot define mean of squares variable name.')
              END IF
              IF (LEN_TRIM(se_info%units) > 59) THEN
                CALL finish('mo_mvstream   , init_mvstream', 'unit name '// &
                     TRIM(se_info%units)//' of variable '//TRIM(se_info%name)// &
                     ' exceeds 59' // &
                     'alphanumeric characters, cannot define mean of squares variable units.')
              END IF
              IF (LEN_TRIM(se_info%longname) > 118) THEN
                CALL finish('mo_mvstream   , init_mvstream', 'longname '// &
                     TRIM(se_info%longname)//' of variable '//TRIM(se_info%name)// &
                     ' exceeds 118' // &
                     'alphanumeric characters, cannot define longname of square.')
              END IF
              SELECT CASE (se_info%ndim)
              CASE (1)
                CALL add_stream_element (mv, TRIM(se_info%name)//'_s', varptr(knam,kstream)%ptr1d_s, &
                     units = '('//TRIM(se_info%units)//')**2', &
                     longname = 'square of '//TRIM(se_info%longname), &
                     ldims = se_info%dim, gdims = se_info%gdim, &
                     klev = se_info%klev, repr = se_info%repr, &
                     table = se_info%gribtable, code = se_info%gribcode, &
                     bits = se_info%gribbits, leveltype = se_info%levelindx, &
                     lmiss = se_info%lmiss, missval = se_info%missval, &
                     lpost = se_info%lpost) ! lrerun must be true
              CASE (2)
                CALL add_stream_element (mv, TRIM(se_info%name)//'_s', varptr(knam,kstream)%ptr2d_s, &
                     units = '('//TRIM(se_info%units)//')**2', &
                     longname = 'square of '//TRIM(se_info%longname), &
                     ldims = se_info%dim, gdims = se_info%gdim, &
                     klev = se_info%klev, repr = se_info%repr, &
                     table = se_info%gribtable, code = se_info%gribcode, &
                     bits = se_info%gribbits, leveltype = se_info%levelindx, &
                     lmiss = se_info%lmiss, missval = se_info%missval, &
                     lpost = se_info%lpost) ! lrerun must be true
              CASE (3)
                CALL add_stream_element (mv, TRIM(se_info%name)//'_s', varptr(knam,kstream)%ptr3d_s, &
                     units = '('//TRIM(se_info%units)//')**2', &
                     longname = 'square of '//TRIM(se_info%longname), &
                     ldims = se_info%dim, gdims = se_info%gdim, &
                     klev = se_info%klev, repr = se_info%repr, &
                     table = se_info%gribtable, code = se_info%gribcode, &
                     bits = se_info%gribbits, leveltype = se_info%levelindx, &
                     lmiss = se_info%lmiss, missval = se_info%missval, &
                     lpost = se_info%lpost) ! lrerun must be true
              END SELECT
            END IF ! lstddev
          END DO ! loop over stream elements
        ENDIF ! tracer stream
      END IF ! there are elements in stream for which calculation is required
    END DO  ! kstream
  END SUBROUTINE init_mvstream

  !
  ! Read in data from metastream namelist groups and store for later use.
  !
  SUBROUTINE read_metastream_data()

    CHARACTER(*), PARAMETER :: func = 'mo_mvstream::read_metastream_data'
    INTEGER :: namelist_id, handle, state, idx, target_idx
    CHARACTER(100) :: target_pattern
    TYPE(metastream) :: buffer

    ! Namelist variables - make sure that this is consistent with include below
    CHARACTER(100) :: target
    TYPE(io_time_event) :: interval
    CHARACTER(100) :: filetag
    CHARACTER(LEN(dummy_stream%name)) :: source(maxstr)
    CHARACTER(LEN(dummy_info%name)) :: meannam(nam_max)
    CHARACTER(LEN(dummy_info%name)) :: sqrmeannam(nam_max)
    CHARACTER(LEN(dummy_stream%name)) :: m_stream_name(maxstr)

    ! Namelist definition - make sure that this is consistent with list above
    INCLUDE 'mvstreamctl.inc'

    ! Loop over all 'mvstreamctl' namelist entries

    ! Initialize global variables for gathering results
    nmstream = 0
    mstreamname = ''

    namelist_id = open_nml('namelist.echam')
    metastream_entries: DO

      ! Read all metastream definitions

      ! Seek next entry, bail out if none is available
      handle = position_nml('mvstreamctl', namelist_id, rewind=.FALSE., status=state)
      IF(state /= positioned) EXIT

      ! Set default values
      target = dummy_meta%target
      interval = dummy_meta%interval
      filetag = dummy_meta%filetag
      source = dummy_meta%source
      meannam = dummy_meta%meannam
      sqrmeannam = dummy_meta%sqrmeannam
      ! @deprecated
      ! Read list of source streams.
      ! This usually only overrides some leading elements of m_stream_name,
      ! so make sure the others are initialized correctly.
      m_stream_name = dummy_meta%source

      ! Read values from namelist file
      READ(handle, mvstreamctl)

      ! Check input data

      IF(m_stream_name(1) /= dummy_meta%source) THEN
        IF(source(1) /= dummy_meta%source) THEN
          ! Do not mix old and new field names
          CALL finish(func, "using both 'm_stream_name' and 'source' is not allowed")
        END IF
        ! Put out fields' values into new structures
        source = m_stream_name
      ELSE IF(source(1) == dummy_meta%source) THEN
        ! Must have at least one source stream
        CALL finish(func, "value for 'source' is missing")
      END IF

      ! Provide default target. Save target pattern for later use
      IF(target == dummy_meta%target) target = '*m'
      target_idx = INDEX(target, '*')
      target_pattern = target

      ! For more than one source stream we cannot use a single target
      ! but must use a target pattern (*)
      ! Note: source(1) /= dummy_meta%source was already checked above
      IF(source(2) /= dummy_meta%source .AND. target_idx == 0) &
           THEN
        CALL finish(func, "must use '*' pattern for 'target' with more than one source stream")
      END IF

      ! Go through source list until first empty element is encountered
      DO idx = 1, maxstr
        IF (source(idx) == dummy_meta%source) EXIT

        ! Generate target name if pattern was given
        IF(target_idx > 0) THEN
          target = target_pattern(:target_idx-1) // TRIM(source(idx)) // target_pattern(target_idx+1:)
        END IF

        CALL check_str_len(func, 'target', target, LEN(dummy_stream%name))
        CALL check_str_len(func, 'filetag', filetag, LEN(dummy_stream%post_suf)-1)
        IF(filetag == dummy_meta%filetag) THEN
          CALL check_str_len(func, 'target', target, LEN(dummy_stream%post_suf)-1, &
               'because it is used as default output file tag')
        END IF

        ! Store data for import into mvstream structure

        IF(nmstream >= maxstr) &
             CALL finish(func, "maximum number of streams reached")
        nmstream = nmstream + 1
        mstreamname(nmstream) = source(idx)

        buffer = metastream(target, interval, filetag, source(idx), meannam, sqrmeannam)

        ! Set dependent defaults
        IF(buffer%filetag == dummy_meta%filetag) buffer%filetag = TRIM(buffer%target)
        ! Defaulting for interval is done in init_mvstream,
        ! after reading of legacy namelist mvctl from separate file.

        metastream_data(nmstream) = buffer

      END DO ! idx = 1, maxstr

    END DO metastream_entries

  END SUBROUTINE read_metastream_data

  !>
  !! Update cache for accumulation info.
  !!
  !! Currently this only checks streams for changed output intervals.
  !!
  !! The other relevant setting would be field accumulation (laccu).
  !! In general, laccu might be changed by customization,
  !! but this would also propose source code changes for doing accu'n.
  !! Thus, the programmer is given responsibility to set laccu consistently,
  !! and changing laccu for source streams via namelist is not supported.
  !!
  !! Changing intervals might lead to incompatible averaging intervals.
  !! This is also checked here.
  !!
  !! Needs to be called after stream customization is finished.
  !!
  SUBROUTINE mvstream_update_cache

    INTEGER :: stream_idx, element_idx
    TYPE(t_stream), POINTER :: source, stream

    ! Retrieve stream info for all registered streams.
    DO stream_idx = 1, nmstream

      ! Get source and cached stream
      CALL get_stream(source, mstreamname(stream_idx))
      stream => stream_cache_data(stream_idx)%stream
      IF(.NOT. ASSOCIATED(stream)) THEN
        CALL message( 'mo_mvstream, update_cache', &
             'skipping empty stream ' // TRIM(mstreamname(stream_idx)))
        CYCLE
      END IF
        
      ! Consistency check

      ! Disable shadow accumulation for incompatible intervals
      IF( events_nonequal(ev_putdata(stream%post_idx), &
           ev_putdata(source%post_idx)) ) &
           THEN
        DO element_idx = 1, nlist_elements(stream_idx)
          IF(varptr(element_idx,stream_idx)%source_accu) THEN

            ! Remove element from stream
            ! and make sure it is never touched during accumulation.
            ! Yes, and issue a warning, just in case.

            CALL remove_stream_element(stream, mvarname(element_idx, stream_idx))

            varptr(element_idx,stream_idx)%source_accu = .FALSE.

            CALL message( 'mo_mvstream, update_cache', &
                 'incompatible intervals: variable ' // &
                 TRIM(mvarname(element_idx, stream_idx)) // &
                 ' removed from stream ' // TRIM(stream%name) )

          END IF! source_accu
        END DO! element_idx
      END IF! events_nonequal

    END DO! stream_idx

  END SUBROUTINE mvstream_update_cache

  !>
  !! Accumulate values for all streams.
  !!
  SUBROUTINE mvstream_accumulate
    CALL mvstream_accumulate_tracer
    CALL mvstream_accumulate_others
  END SUBROUTINE mvstream_accumulate

  !>
  !! Accumulate values for tracer stream.
  !!
  SUBROUTINE mvstream_accumulate_tracer
    INTEGER :: kstream, kvar

    DO kstream = 1, nmstream
      IF (mstreamname(kstream) == 'tracer') THEN
        ! Check elements from tracer stream.
        DO kvar = 1, trlist%ntrac
          IF (lvarmean(kvar, kstream)) THEN
            ! Only variables that are scheduled for averaging.
            varptr(kvar,kstream)%ptr3d(:, :, :ldc%ngpblks-1) = &
                 varptr(kvar,kstream)%ptr3d(:, :, :ldc%ngpblks-1) + &
                 xt(:, :, kvar, :ldc%ngpblks-1) * delta_time
            ! For row ngpblks ignore possibly undefined values
            varptr(kvar,kstream)%ptr3d(:ldc%npromz, :, ldc%ngpblks) = &
                 varptr(kvar,kstream)%ptr3d(:ldc%npromz, :, ldc%ngpblks) + &
                 xt(:ldc%npromz, :, kvar, ldc%ngpblks) * delta_time
            IF (lstddev(kvar,kstream)) THEN
              ! Add square of current value to standard deviation.
              varptr(kvar,kstream)%ptr3d_s(:, :, :ldc%ngpblks-1) = &
                   varptr(kvar,kstream)%ptr3d_s(:, :, :ldc%ngpblks-1) + &
                   xt(:, :, kvar, :ldc%ngpblks-1) * &
                   xt(:, :, kvar, :ldc%ngpblks-1) * delta_time
              ! For row ngpblks ignore possibly undefined values
              varptr(kvar,kstream)%ptr3d_s(:ldc%npromz, :, ldc%ngpblks) = &
                   varptr(kvar,kstream)%ptr3d_s(:ldc%npromz, :, ldc%ngpblks) + &
                   xt(:ldc%npromz, :, kvar, ldc%ngpblks) * &
                   xt(:ldc%npromz, :, kvar, ldc%ngpblks) * delta_time
            END IF! lstddev
          END IF! lvarmean
        END DO! kvar
      END IF! mstreamname == 'tracer'
    END DO! kstream
  END SUBROUTINE mvstream_accumulate_tracer

  !>
  !! Accumulate values for all source streams besides tracer.
  !!
  SUBROUTINE mvstream_accumulate_others
    INTEGER :: kstream, kvar

    DO kstream = 1, nmstream
      IF (mstreamname(kstream) /= 'tracer') THEN
        ! Check elements from all non-tracer streams.
        DO kvar = 1, nlist_elements(kstream)
          IF(lvarmean(kvar, kstream)) THEN
            ! Only variables that are scheduled for averaging.
            SELECT CASE (ivarrank(kvar,kstream))
            CASE (1)
              ! For one-dimensional variables, add current value.
              ! Cannot be gridpoint variable, so no npromz handling
              varptr(kvar,kstream)%ptr1d = &
                   varptr(kvar,kstream)%ptr1d + &
                   tarptr(kvar,kstream)%ptt1d * delta_time
              IF (lstddev(kvar,kstream)) THEN
                ! Add square of current value to standard deviation.
                varptr(kvar,kstream)%ptr1d_s = &
                     varptr(kvar,kstream)%ptr1d_s + &
                     tarptr(kvar,kstream)%ptt1d * tarptr(kvar,kstream)%ptt1d * &
                     delta_time
              END IF! lstddev
            CASE (2)
              ! For two-dimensional variables, add current value.
              IF (lgridpoint(kvar,kstream)) THEN
                varptr(kvar,kstream)%ptr2d(:, :ldc%ngpblks-1) = &
                     varptr(kvar,kstream)%ptr2d(:, :ldc%ngpblks-1) + &
                     tarptr(kvar,kstream)%ptt2d(:, :ldc%ngpblks-1) * delta_time
                ! For row ngpblks ignore possibly undefined values
                varptr(kvar,kstream)%ptr2d(:ldc%npromz, ldc%ngpblks) = &
                     varptr(kvar,kstream)%ptr2d(:ldc%npromz, ldc%ngpblks) + &
                     tarptr(kvar,kstream)%ptt2d(:ldc%npromz, ldc%ngpblks) * delta_time
              ELSE
                varptr(kvar,kstream)%ptr2d = &
                     varptr(kvar,kstream)%ptr2d + &
                     tarptr(kvar,kstream)%ptt2d * delta_time
              END IF
              IF (lstddev(kvar,kstream)) THEN
                ! Add square of current value to standard deviation.
                IF (lgridpoint(kvar,kstream)) THEN
                  varptr(kvar,kstream)%ptr2d_s(:, :ldc%ngpblks-1) = &
                       varptr(kvar,kstream)%ptr2d_s(:, :ldc%ngpblks-1) + &
                       tarptr(kvar,kstream)%ptt2d(:, :ldc%ngpblks-1) * &
                       tarptr(kvar,kstream)%ptt2d(:, :ldc%ngpblks-1) * &
                       delta_time
                  ! For row ngpblks ignore possibly undefined values
                  varptr(kvar,kstream)%ptr2d_s(:ldc%npromz, ldc%ngpblks) = &
                       varptr(kvar,kstream)%ptr2d_s(:ldc%npromz, ldc%ngpblks) + &
                       tarptr(kvar,kstream)%ptt2d(:ldc%npromz, ldc%ngpblks) * &
                       tarptr(kvar,kstream)%ptt2d(:ldc%npromz, ldc%ngpblks) * &
                       delta_time
                ELSE
                  varptr(kvar,kstream)%ptr2d_s = &
                       varptr(kvar,kstream)%ptr2d_s + &
                       tarptr(kvar,kstream)%ptt2d * tarptr(kvar,kstream)%ptt2d * &
                       delta_time
                END IF
              END IF! lstddev
            CASE (3) ! ivarrank
              ! For three-dimensional variables, add current value.
              IF (lgridpoint(kvar,kstream)) THEN
                varptr(kvar,kstream)%ptr3d(:, :, :ldc%ngpblks-1) = &
                     varptr(kvar,kstream)%ptr3d(:, :, :ldc%ngpblks-1) + &
                     tarptr(kvar,kstream)%ptt3d(:, :, :ldc%ngpblks-1) * delta_time
                ! For row ngpblks ignore possibly undefined values
                varptr(kvar,kstream)%ptr3d(:ldc%npromz, :, ldc%ngpblks) = &
                     varptr(kvar,kstream)%ptr3d(:ldc%npromz, :, ldc%ngpblks) + &
                     tarptr(kvar,kstream)%ptt3d(:ldc%npromz, :, ldc%ngpblks) * delta_time
              ELSE
                varptr(kvar,kstream)%ptr3d = &
                     varptr(kvar,kstream)%ptr3d + &
                     tarptr(kvar,kstream)%ptt3d * delta_time
              END IF
              IF (lstddev(kvar,kstream)) THEN
                ! Add square of current value to standard deviation.
                IF (lgridpoint(kvar,kstream)) THEN
                  varptr(kvar,kstream)%ptr3d_s(:, :, :ldc%ngpblks-1) = &
                       varptr(kvar,kstream)%ptr3d_s(:, :, :ldc%ngpblks-1) + &
                       tarptr(kvar,kstream)%ptt3d(:, :, :ldc%ngpblks-1) * &
                       tarptr(kvar,kstream)%ptt3d(:, :, :ldc%ngpblks-1) * &
                       delta_time
                  ! For row ngpblks ignore possibly undefined values
                  varptr(kvar,kstream)%ptr3d_s(:ldc%npromz, :, ldc%ngpblks) = &
                       varptr(kvar,kstream)%ptr3d_s(:ldc%npromz, :, ldc%ngpblks) + &
                       tarptr(kvar,kstream)%ptt3d(:ldc%npromz, :, ldc%ngpblks) * &
                       tarptr(kvar,kstream)%ptt3d(:ldc%npromz, :, ldc%ngpblks) * &
                       delta_time
                ELSE
                  varptr(kvar,kstream)%ptr3d_s = &
                       varptr(kvar,kstream)%ptr3d_s + &
                       tarptr(kvar,kstream)%ptt3d * tarptr(kvar,kstream)%ptt3d * &
                       delta_time
                END IF
              END IF! lstddev
            END SELECT! ivarrank
          END IF! lvarmean
        END DO! kvar
      END IF! mstreamname /= 'tracer'
    END DO! kstream
  END SUBROUTINE mvstream_accumulate_others

  !
  ! Utility routines
  !

  SUBROUTINE check_str_len(section, key, string, max_length, message)
    CHARACTER(*), INTENT(in) :: section
    CHARACTER(*), INTENT(in) :: key
    CHARACTER(*), INTENT(in) :: string
    INTEGER, INTENT(in) :: max_length
    CHARACTER(*), INTENT(in), OPTIONAL :: message

    CHARACTER(10) :: max_length_str

    IF(LEN_TRIM(string) > max_length) THEN
      max_length_str = int2string(max_length)
      IF(PRESENT(message)) THEN
        CALL finish(section, "value for '" // key // "' may not exceed " // &
             TRIM(max_length_str) // " characters " // message)
      ELSE
        CALL finish(section, "value for '" // key // "' may not exceed " // &
             TRIM(max_length_str) // " characters")
      END IF
    END IF
  END SUBROUTINE check_str_len

  ! @todo Remove this if no longer needed
  !  subroutine check_str_val(section, key, string, allowed, message)
  !    character(*), intent(in) :: section
  !    character(*), intent(in) :: key
  !    character(*), intent(in) :: string
  !    character(*), intent(in) :: allowed(:)
  !    character(*), intent(in), optional :: message
  !
  !    if(all(allowed /= string)) then
  !        if(present(message)) then
  !            call finish(section, "value for '" // key // "' is not supported " // &
  !                        message)
  !        else
  !            call finish(section, "value for '" // key // "' is not supported")
  !        end if
  !    end if
  !  end subroutine

END MODULE mo_mvstream
