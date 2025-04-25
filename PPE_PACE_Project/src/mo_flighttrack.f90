!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham.f90
!!
!! \brief
!! Flight-track simulator. Interpolates model variables onto a specified flight track.
!!
!! \author Zak Kipling, Department of Physics, University of Oxford
!!
!! \responsible_coder
!! Zak Kipling, zak.kipling@physics.ox.ac.uk
!!
!! \revision_history
!!   -# Zak Kipling (Uni Oxford) - ported from ECHAM5.5-HAM2 (2013)
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

MODULE mo_flighttrack

  USE mo_kind,            ONLY: dp, wp
  USE mo_linked_list,     ONLY: t_stream
  USE mo_time_conversion, ONLY: time_days
  USE mo_tracdef,         ONLY: ln

  IMPLICIT NONE
  PRIVATE

  ! Fixed parameters
  INTEGER, PARAMETER :: max_tracers = 50  ! Maximum number of tracers to sample (old method)

  ! Type definitions
  TYPE t_track_point
    TYPE(time_days) :: datetime ! Date and time
    REAL(wp) :: lat             ! Latitude (degN)
    REAL(wp) :: lon             ! Latitude (degE)
    REAL(wp) :: p               ! Pressure level (Pa)
  END TYPE

  TYPE t_element_map
    TYPE(t_stream), POINTER :: stream
    CHARACTER(len=64)       :: name
    CHARACTER(len=64)       :: units
    REAL(dp),       POINTER :: src(:,:,:)
    INTEGER                 :: ktra
  END TYPE t_element_map

  ! Switches
  LOGICAL :: lseconds     = .FALSE.  ! Use seconds or only minutes

  ! I/O configuration
  INTEGER            :: track_max_points = 0     ! Maximum # points in flight track
  INTEGER            :: max_stream_elements = 50 ! Maximum number of stream elements to sample
  CHARACTER(LEN=ln)  :: tracer_names(max_tracers) = '' ! Names of tracers to sample (old method)

  CHARACTER(LEN=100) :: track_filename = ''      ! Name of text file containing track
  CHARACTER(LEN=100) :: output_filename = ''     ! Name of text file containing output

  ! Flight track data
  INTEGER :: n_track_points = 0                           ! Actual # points in track
  TYPE(t_track_point), ALLOCATABLE, SAVE :: track_data(:) ! Track specification

  ! Output specification
  INTEGER :: n_stream_elements = 0          ! Actual # stream elements to sample
  TYPE(t_element_map), ALLOCATABLE, SAVE, TARGET :: element_map(:)

  ! Runtime data
  INTEGER :: nout = 0                    ! Unit number for output file
  INTEGER :: track_pos = 0               ! Index into track_data of next point expected
  CHARACTER(LEN=100) :: out_format  = '' ! Format for output lines


  PUBLIC :: flighttrack_initialize
  PUBLIC :: flighttrack_init_stream_elements
  PUBLIC :: flighttrack_init_output
  PUBLIC :: flighttrack_diag
  PUBLIC :: flighttrack_finalize_output

  PRIVATE :: fracgeoindex
  PRIVATE :: fracplevindex
  PRIVATE :: flighttrack_interpolate
  PRIVATE :: interp_axis

CONTAINS

  SUBROUTINE flighttrack_initialize
  
    USE mo_exception,       ONLY: finish, message, message_text
    USE mo_filename,        ONLY: find_next_free_unit
    USE mo_mpi,             ONLY: p_parallel, p_parallel_io, p_bcast, p_io
    USE mo_namelist,        ONLY: close_nml, open_nml, position_nml, POSITIONED
    USE mo_time_conversion, ONLY: tc_convert, tc_get, tc_set, time_native, operator(<)

    IMPLICIT NONE

    INCLUDE 'flighttrackctl.inc'

    INTEGER :: ierr, inml, iunit, itrack, ipt
    TYPE(t_track_point) :: pt
    INTEGER :: year, month, day, hour, minute, second, p_int
    INTEGER, ALLOCATABLE :: all_days(:), all_seconds(:)
    TYPE(time_native) :: datetime

    IF (p_parallel_io) THEN
      ! Read namelist
      inml = open_nml('namelist.echam')
      iunit = position_nml ('FLIGHTTRACKCTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
         READ (iunit, flighttrackctl)
      CASE DEFAULT
         WRITE(message_text,'(a,i0)') 'Namelist flighttrackctl not correctly read! ierr = ', ierr
         CALL finish('flighttrack_initialize', message_text)
      END SELECT
      CALL close_nml(inml)
    END IF

    IF (p_parallel) THEN
      CALL p_bcast(lseconds,         p_io)
      CALL p_bcast(track_max_points, p_io)
      CALL p_bcast(track_filename,   p_io)
      CALL p_bcast(output_filename,  p_io)
      CALL p_bcast(tracer_names,     p_io)
    END IF

    ALLOCATE(track_data(track_max_points))

    IF (p_parallel_io) THEN
      ! Read flight track
      itrack = find_next_free_unit(20, 99)
      OPEN(itrack, file=track_filename, action='READ', &
           status='OLD', form='FORMATTED')
      READ(itrack,*) ! Skip header line
      ipt = 0
      DO WHILE(.TRUE.)
        IF (lseconds) THEN
          READ(itrack, fmt='(i4,2(i2),i3,2(i2),f10.4,f9.4,i5)',     &
               iostat=ierr) year, month, day, hour, minute, second, &
                            pt%lon, pt%lat, p_int
        ELSE
          READ(itrack, fmt='(i4,2(i2),i3,i2,f10.4,f9.4,i5)', &
               iostat=ierr) year, month, day, hour, minute,  &
                            pt%lon, pt%lat, p_int
          second = 0
        END IF
        SELECT CASE(ierr)
          CASE(-1)
            EXIT ! End of file
          CASE(0)
            ! Successful read
            CALL tc_set(year, month, day, hour, minute, second, datetime)
            CALL tc_convert(datetime, pt%datetime)
            IF (ipt > 0) THEN
              ! Not .AND. because Fortran doesn't do short-circuit logic.
              IF (pt%datetime < track_data(ipt)%datetime) THEN
                CALL finish('flighttrack_initialize', &
                            'Flight track times not monotonically increasing')
              END IF
            END IF
            pt%lon = MODULO(pt%lon, 360._wp) ! Ensure longitude in [0,360)
            pt%p = p_int * 100 ! Convert hPa to Pa
            WRITE (message_text,                                           &
                   fmt='(A,i7,1x,i4.4,2(i2.2),1x,3(i2.2),f8.3,f9.3,f9.1)') &
                  'Track pt', ipt+1, year, month, day,                     &
                  hour, minute, second, pt%lon, pt%lat, pt%p
            CALL message('flighttrack_initialize', message_text)
            IF (pt%lat < -90._wp .OR. pt%lat > 90._wp) THEN
              CALL finish('flighttrack_initialize', 'Track latitude out of range')
            END IF
            IF (ipt < track_max_points) THEN
              ipt = ipt + 1
              track_data(ipt) = pt
            ELSE
              CALL finish('flighttrack_initialize', &
                          'Too many points in flight track')
            END IF
          CASE DEFAULT
            CALL finish('flighttrack_initialize', 'Error reading track file')
        END SELECT
      END DO
      CLOSE(itrack)
      n_track_points = ipt
      WRITE (message_text, fmt='(A,i7,A)') 'Read', n_track_points, &
                                  ' flight track points'
      CALL message('flighttrack_initialize', message_text)
    ENDIF ! p_parallel_io
    track_pos = 1

    IF (p_parallel) THEN
      CALL p_bcast(n_track_points,         p_io)
      CALL p_bcast(track_data(:)%lat,      p_io)
      CALL p_bcast(track_data(:)%lon,      p_io)
      CALL p_bcast(track_data(:)%p,        p_io)

      ! Can't p_bcast TYPE(time_days) directly...
      ALLOCATE(all_days(track_max_points))
      ALLOCATE(all_seconds(track_max_points))

      IF (p_parallel_io) THEN
        DO ipt=1,n_track_points
          CALL tc_get(track_data(ipt)%datetime,        &
                      all_days(ipt), all_seconds(ipt))
        END DO
      END IF

      CALL p_bcast(all_days,    p_io)
      CALL p_bcast(all_seconds, p_io)

      IF (.NOT. p_parallel_io) THEN
        DO ipt=1,n_track_points
          CALL tc_set(all_days(ipt), all_seconds(ipt), &
                      track_data(ipt)%datetime)
        END DO
      END IF

      DEALLOCATE(all_seconds)
      DEALLOCATE(all_days)
    END IF

  END SUBROUTINE flighttrack_initialize
   

  SUBROUTINE flighttrack_init_stream_elements
  
    USE mo_exception,   ONLY: finish, message, message_text
    USE mo_linked_list, ONLY: GAUSSIAN, HYBRID, UNKNOWN, &
                              memory_info
    USE mo_memory_base, ONLY: new_stream, get_stream,  &
                              add_stream_element,      &
                              default_stream_setting,  &
                              get_stream_element_info, & 
                              get_stream_element,      &
                              nstreams, ostreams, AUTO
    USE mo_memory_g1a,  ONLY: xtm1
    USE mo_mpi,         ONLY: p_bcast, p_parallel, p_parallel_io, p_io
    USE mo_namelist,    ONLY: open_nml, position_nml, close_nml, &
                              POSITIONED
    USE mo_time_control,ONLY: p_bcast_event
    USE mo_tracdef,     ONLY: trlist
    USE mo_tracer,      ONLY: get_tracer
    !USE mo_util_string, ONLY: tolower

    IMPLICIT NONE

    CHARACTER(len=64)            :: stream, element
    INTEGER                      :: i, ierr, inml, iunit, kstream
    TYPE(memory_info)            :: info
    TYPE(t_stream),      POINTER :: strm
    TYPE(t_element_map), POINTER :: el

    INCLUDE 'set_flighttrack_element.inc'

    n_stream_elements = 0
    IF (max_stream_elements> 0) THEN
      ALLOCATE(element_map(max_stream_elements))

      ! Identify tracers to output (old method)
      DO i=1,max_tracers
        IF (LEN_TRIM(tracer_names(i)) > 0) THEN
          n_stream_elements = n_stream_elements + 1
          el => element_map(n_stream_elements)
          CALL get_tracer(TRIM(tracer_names(i)), idx=el%ktra, ierr=ierr)
          IF (ierr == 0) THEN
            WRITE (message_text, fmt='(A,i4,3(A),i4,A,A)')                  &
                  'Found tracer', n_stream_elements, ' (',                  &
                  TRIM(trlist%ti(el%ktra)%fullname), ') at index', el%ktra, &
                  ' with units ', TRIM(trlist%ti(el%ktra)%units)
            CALL message('flighttrack_init_stream_elements', message_text)
            el%name=trlist%ti(el%ktra)%fullname
            el%units=trlist%ti(el%ktra)%units
            NULLIFY(el%src)
          ELSE
            WRITE (message_text, fmt='(3(A))') 'Tracer ',       &
                                      TRIM(tracer_names(i)),    &
                                      ' not found'
            CALL finish('flighttrack_init_stream_elements', message_text)
          END IF
        END IF
      END DO

      inml = open_nml('namelist.echam')
      DO
        IF (p_parallel_io) THEN
          stream = ''
          element = ''
          iunit = position_nml('SET_FLIGHTTRACK_ELEMENT', inml, &
                               REWIND=.false., STATUS=ierr)
          IF (ierr == POSITIONED) THEN
            READ (iunit, set_flighttrack_element)
            !stream = tolower(stream)
          END IF
        ENDIF

        IF (p_parallel) THEN
          CALL p_bcast(ierr,    p_io)
          CALL p_bcast(stream,  p_io)
          CALL p_bcast(element, p_io)
        ENDIF

        IF (ierr /= POSITIONED) EXIT

        IF (LEN_TRIM(stream) == 0 .OR. TRIM(stream) == 'tracer') THEN
          el => element_map(n_stream_elements+1)

          CALL get_tracer(TRIM(element), idx=el%ktra, ierr=ierr)
          IF (ierr == 0) THEN
            WRITE (message_text, fmt='(A,i4,3(A),i4,A,A)')                  &
                  'Found tracer', n_stream_elements, ' (',                  &
                  TRIM(trlist%ti(el%ktra)%fullname), ') at index', el%ktra, &
                  ' with units ', TRIM(trlist%ti(el%ktra)%units)
            CALL message('flighttrack_init_stream_elements', message_text)

            n_stream_elements = n_stream_elements + 1

            el%name=trlist%ti(el%ktra)%fullname
            el%units=trlist%ti(el%ktra)%units
            NULLIFY(el%src)
            CYCLE
          ELSE IF (LEN_TRIM(stream) /= 0) THEN
            WRITE (message_text, fmt='(3(A))') 'Tracer ',         &
                                      TRIM(element), ' not found'
            CALL finish('flighttrack_init_stream_elements', message_text)
          END IF
        END IF

        NULLIFY(strm)
        DO i = 1, nstreams
          IF (LEN_TRIM(stream) == 0 .OR. &
              TRIM(stream) == TRIM(ostreams(i)%name)) THEN
            strm => ostreams(i)
            CALL get_stream_element_info(strm, element, info)
            IF (LEN_TRIM(info%name) /= 0) EXIT
          ENDIF
        END DO
        IF (.NOT. ASSOCIATED(strm)) THEN
          CALL finish ('flighttrack_init_stream_elements','cannot find stream' &
                       //TRIM(stream))
        END IF

        IF (LEN_TRIM(info%name) == 0) THEN
          CALL finish ('flighttrack_init_stream_elements','cannot find ' &
                       //TRIM(element)//' in stream '//TRIM(stream))
        END IF

        IF (info%ndim == UNKNOWN) THEN
          CALL finish('flighttrack_init_stream_elements','stream element ' &
                      //TRIM(element)//' in stream '//TRIM(stream)         &
                      //' is not defined')
        ELSE IF (info%ndim /= 3) THEN
          WRITE (*,*) 'ndim', info%ndim
          CALL finish('flighttrack_init_stream_elements','stream element ' &
                      //TRIM(element)//' in stream '//TRIM(stream)         &
                      //' is not three-dimensional')
        ELSE IF (info%repr /= GAUSSIAN) THEN
          CALL finish('flighttrack_init_stream_elements','stream element ' &
                      //TRIM(element)//' in stream '//TRIM(stream)         &
                      //': not gaussian')
        ELSE IF (info%levelindx /= HYBRID) THEN
          CALL finish('flighttrack_init_stream_elements','stream element ' &
                      //TRIM(element)//' in stream '//TRIM(stream)         &
                      //': not on model levels')
        ELSE IF (info%laccu) THEN
          CALL finish('flighttrack_init_stream_elements','stream element ' &
                      //TRIM(element)//' in stream '//TRIM(stream)         &
                      //' is an accumulated field')
        END IF

        n_stream_elements = n_stream_elements + 1
        el => element_map(n_stream_elements)

        el%stream => strm
        el%name   = info%name
        el%units  = info%units
        el%ktra   = 0

        CALL get_stream_element(strm, TRIM(element), el%src)

        WRITE (message_text, fmt='(A,i4,6(A))')          &
              'Found element', n_stream_elements,        &
              ' (', TRIM(element), ' in stream ',        &
              TRIM(stream), ') with units ', info%units

        CALL message('flighttrack_init_stream_elements', message_text)

      END DO  ! loop over stream elements
      CALL close_nml(inml)

    END IF ! max_stream_elements > 0

  END SUBROUTINE flighttrack_init_stream_elements
   

  SUBROUTINE flighttrack_init_output
  
    USE mo_filename,  ONLY: find_next_free_unit
    USE mo_mpi,       ONLY: p_parallel_io

    IMPLICIT NONE

    CHARACTER(LEN=2*64+10) :: head(n_stream_elements)
    CHARACTER(LEN=100) :: fmt
    INTEGER :: max_head_len
    INTEGER :: i
    INTEGER :: ierr
    TYPE(t_element_map), POINTER :: el

    IF (p_parallel_io) THEN
      ! Set up header format
      max_head_len = 12 ! To allow space for numeric output
      DO i=1,n_stream_elements
        el => element_map(i)
        IF (LEN_TRIM(el%units) > 0) THEN
          WRITE (head(i), fmt='(4A)') TRIM(el%name), ' (', TRIM(el%units), ')'
        ELSE
          head(i) = el%name
        END IF
        max_head_len = MAX(max_head_len, LEN_TRIM(head(i)))
      END DO
      WRITE (fmt, fmt='(A,i3,A,i2.2,A)') '(2(1x,a8),6(1x,a10),', &
          n_stream_elements, '(1x,a', max_head_len, '))'

      ! Append to existing output file or create new one
      nout = find_next_free_unit(20, 99)
      OPEN(nout, file=output_filename, form='FORMATTED', action='WRITE', &
                 status='OLD', position='APPEND', iostat=ierr)
      IF (ierr /= 0) THEN
        ! No existing output file -- create new one and write header
        OPEN(nout, file=output_filename, form='FORMATTED', action='WRITE', &
                   status='NEW')
        WRITE (nout,fmt)                                        &
                   'Time', 'YYYYMMDD','lon(deg)', 'lat(deg)',   &
                   'p(Pa)', 'h(m)', 'T(K)', 'q', head
      END IF

      ! Set up output format
      WRITE (out_format, fmt='(a,i3,a,i2.2,a)')                                         &
            '(1x,i8,1x,i4.4,2(i2.2),1x,2(f10.4,1x),f10.1,1x,f10.3,1x,f10.4,1x,es10.4,', &
            n_stream_elements, '(1x,es', max_head_len, '.5))'
    END IF ! p_parallel_io

  END SUBROUTINE flighttrack_init_output


  SUBROUTINE flighttrack_diag
  
    USE mo_ham_tools,          ONLY: geoindex
    USE mo_physical_constants, ONLY: grav, rd, vtmpc1
    USE mo_control,            ONLY: nlon, nlev, ngl
    USE mo_decomposition,      ONLY: gd => global_decomposition, &
                                     ld => local_decomposition
    USE mo_exception,          ONLY: finish, message, message_text
    USE mo_memory_g1a,         ONLY: alpsm1, tm1, qm1, xim1, xlm1, xtm1
    USE mo_memory_g3b,         ONLY: oromea
    USE mo_mpi,                ONLY: p_parallel_io, p_bcast, p_io, p_pe
    USE mo_scan_buffer,        ONLY: alnpr, alpha
    USE mo_time_control,       ONLY: delta_time, previous_date
    USE mo_time_conversion,    ONLY: add_date, day_in_year, tc_convert, tc_get, &
                                     time_native, operator(>)
    USE mo_transpose,          ONLY: gather_gp3


    IMPLICIT NONE

    ! Local arrays
    REAL(dp), POINTER :: local_ap(:,:,:)
    REAL(dp), POINTER :: local_geoh(:,:,:)
    REAL(dp), POINTER :: local_aph_tmp(:,:)
    REAL(dp), POINTER :: local_ztv_tmp(:,:)
    REAL(dp), POINTER :: local_zgeo_tmp(:)
    REAL(dp), POINTER :: local_tmp(:,:,:)

    ! Gathered arrays
    REAL(dp), POINTER :: gathered_ap(:,:,:)
    REAL(dp), POINTER :: gathered_output(:,:,:,:)
    REAL(dp), POINTER :: gathered_tmp(:,:,:)

    INTEGER :: ifield, irow
    REAL(dp) :: ilon, ilev, ilat ! fractional lon/lev/lat indices
    INTEGER :: time_val ! time value for output (minute/second in year)
    TYPE(time_days) :: cutoff_date 
    TYPE(time_native) :: dt_native
    INTEGER :: year, month, day, hour, minute, second
    INTEGER :: nproma_tmp
    LOGICAL :: l_output_any

    TYPE(t_element_map), POINTER :: el

    EXTERNAL :: geopot, pres, presf

    cutoff_date = previous_date
    CALL add_date(0, -FLOOR(delta_time), cutoff_date)
    l_output_any = .FALSE.
    DO WHILE(track_pos <= n_track_points)
      IF (track_data(track_pos)%datetime > previous_date) THEN
        ! We haven't reached this point yet
        EXIT
      ELSE IF (track_data(track_pos)%datetime > cutoff_date) THEN
        ! Output this point now
        ! Allocate gathered fields on I/O processor
        IF (p_parallel_io) THEN
          IF (.NOT. l_output_any) THEN
            ALLOCATE(gathered_ap(nlon,nlev,ngl))
            ALLOCATE(gathered_output(ld%nlon,ld%nlev,ld%nlat,n_stream_elements+3))
          END IF
        ELSE
          NULLIFY(gathered_ap)
          NULLIFY(gathered_output)
          NULLIFY(gathered_tmp)
        END IF

        IF (.NOT. l_output_any) THEN
          ! Calculate local pressure and geopotential fields
          ! These only appear to exist on a per-block basis within
          ! physc, not as full per-domain fields, so we replicate
          ! the calculations here
          ALLOCATE(local_ap(ld%nproma,ld%nlev,ld%ngpblks))
          ALLOCATE(local_geoh(ld%nproma,ld%nlev,ld%ngpblks))
          ALLOCATE(local_aph_tmp(ld%nproma,ld%nlev+1))
          ALLOCATE(local_ztv_tmp(ld%nproma,ld%nlev))
          ALLOCATE(local_zgeo_tmp(ld%nproma))

          DO irow=1,ld%ngpblks
            IF (irow == ld%ngpblks) THEN
              nproma_tmp = ld%npromz
            ELSE
              nproma_tmp = ld%nproma
            END IF
            CALL pres(local_aph_tmp, ld%nproma, EXP(alpsm1(:,irow)), nproma_tmp)
            CALL presf(local_ap(:,:,irow), ld%nproma, local_aph_tmp, nproma_tmp)
            local_ztv_tmp(1:nproma_tmp,:) = tm1(1:nproma_tmp,:,irow)*(1._dp+vtmpc1*qm1(1:nproma_tmp,:,irow) &
                                                                      -(xlm1(1:nproma_tmp,:,irow)           &
                                                                        +xim1(1:nproma_tmp,:,irow)))
            local_zgeo_tmp(1:nproma_tmp) = 0._dp
            CALL geopot(local_geoh(:,:,irow), local_ztv_tmp(:,:), alnpr(:,:,irow), alpha(:,:,irow), &
                        local_zgeo_tmp(:), ld%nproma, nproma_tmp)
            ! convert geopotential to (geometric surface height above msl) + (geopotential height above surface)
            local_geoh(1:nproma_tmp,:,irow) = local_geoh(1:nproma_tmp,:,irow) / grav + SPREAD(oromea(1:nproma_tmp,irow),2,nlev)
          END DO

          ! Gather pressure field
          CALL gather_gp3(gathered_ap, local_ap, gd)

          ! Gather geopotential height field
          IF (p_parallel_io) THEN
            gathered_tmp => gathered_output(:,:,:,1) 
          END IF
          CALL gather_gp3(gathered_tmp, local_geoh, gd)

          ! Gather temperature field
          IF (p_parallel_io) THEN
            gathered_tmp => gathered_output(:,:,:,2) 
          END IF
          CALL gather_gp3(gathered_tmp, tm1, gd)

          ! Gather humidity field
          IF (p_parallel_io) THEN
            gathered_tmp => gathered_output(:,:,:,3) 
          END IF
          CALL gather_gp3(gathered_tmp, qm1, gd)

          ! Gather output fields
          DO ifield=1,n_stream_elements
            el => element_map(ifield)
            IF (ASSOCIATED(el%src)) THEN
              local_tmp => el%src
            ELSE IF (el%ktra > 0) THEN
              local_tmp => xtm1(:,:,el%ktra,:)
            ELSE
              WRITE (message_text, fmt='(a,i4,3A)') &
                    'no data for element', ifield, ' (', el%name, ')'
              CALL finish('flighttrack_diag', message_text)
            END IF

            ! It seems we need to pass a named pointer here, rather than just
            ! referring directly to the array section...
            IF (p_parallel_io) THEN
              gathered_tmp => gathered_output(:,:,:,ifield+3) 
            END IF

            CALL gather_gp3(gathered_tmp, local_tmp, gd)
          END DO

          ! Deallocate local fields after gathering
          DEALLOCATE(local_zgeo_tmp)
          DEALLOCATE(local_ztv_tmp)
          DEALLOCATE(local_aph_tmp)
          DEALLOCATE(local_geoh)
          DEALLOCATE(local_ap)

          ! We only need to prepare the gathered arrays once at each model
          ! timestep, regardless of how many points are being output.
          l_output_any = .TRUE.
        END IF

        IF (p_parallel_io) THEN
          ! Interpolate and output
          CALL fracgeoindex(track_data(track_pos)%lon, &
                            track_data(track_pos)%lat, &
                            ilon, ilat)
          IF (ilon /= -999._dp .AND. ilat /= -999._dp) THEN
            CALL fracplevindex(flighttrack_interpolate_horiz(gathered_ap(:,:,:), &
                                                             ilon, ilat),        &
                               track_data(track_pos)%p, ilev)
            IF (ilev /= -999._dp) THEN
              time_val = (day_in_year(track_data(track_pos)%datetime)-1)
              CALL tc_convert(track_data(track_pos)%datetime, dt_native)
              CALL tc_get(dt_native, year, month, day, hour, minute, second)

              write (message_text,                                                  &
                     fmt='(a,i4.4,2(i2.2),1x,3(i2.2),2(f10.4,1x),f10.1,a,3(f8.3))') &
                    'sampling at ', year, month, day, hour, minute, second,         &
                    track_data(track_pos)%lon, track_data(track_pos)%lat,           &
                    track_data(track_pos)%p, ': indices', ilon, ilev, ilat
              call message('flighttrack_diag', message_text)

              time_val = (time_val * 24 + hour) * 60 + minute
              IF (lseconds) THEN
                time_val = time_val * 60 + second
              END IF
              WRITE(nout, out_format) time_val, year, month, day,     &
                                      track_data(track_pos)%lon,      &
                                      track_data(track_pos)%lat,      &
                                      track_data(track_pos)%p,        &
                                      flighttrack_interpolate(        &
                                        gathered_output(:,:,:,:),     &
                                        ilon,ilev,ilat                &
                                      )
            END IF
          END IF
        END IF ! p_parallel_io
     !ELSE we're already too late for this point
      END IF
      track_pos = track_pos + 1
    END DO ! track_pos

    IF (l_output_any .AND. p_parallel_io) THEN
      ! Deallocate gathered fields
      DEALLOCATE(gathered_output)
      DEALLOCATE(gathered_ap)
    END IF

  END SUBROUTINE flighttrack_diag


  SUBROUTINE flighttrack_finalize_output
  
    USE mo_mpi,       ONLY: p_parallel_io

    IMPLICIT NONE

    IF (p_parallel_io) THEN
      ! Close output file
      CLOSE(nout)
    END IF

  END SUBROUTINE flighttrack_finalize_output

  SUBROUTINE fracgeoindex (plon, plat, klon, klat)

    USE mo_exception, ONLY: finish, message
    USE mo_control,   ONLY: ngl, nlon
    USE mo_gaussgrid, ONLY: philat, philon
    IMPLICIT NONE

    REAL(dp),INTENT(IN)  :: plon, plat   ! Coordinates in degrees
    REAL(dp),INTENT(OUT) :: klon, klat   ! Corresponding indices

    INTEGER             :: jlon, jlat
      

    IF(plon<0._dp.OR.plon>360._dp.OR.plat<-90._dp.OR.plat>90._dp) THEN
      CALL finish('geoindex:', 'Coordinates out of range')
    END IF

    klon=-999._dp
    klat=-999._dp

    DO jlon = 1, nlon
      IF (plon <= philon(jlon)) THEN
        IF (jlon > 1) THEN
          klon = jlon - (philon(jlon) - plon) / (philon(jlon) - philon(jlon-1))
        ELSE
          klon = jlon - (philon(jlon) - plon) / (philon(jlon) - (philon(nlon)-360._dp))
        END IF
        EXIT
      END IF
    END DO
    IF (klon==-999._dp) THEN
      klon = nlon + (plon - philon(nlon)) / (philon(1)+360._dp - philon(nlon))
    END IF

    DO jlat = 1, ngl
      IF (plat >= philat(jlat)) THEN
        IF (jlat > 1) THEN
          klat = jlat - (plat - philat(jlat)) / (philat(jlat-1) - philat(jlat))
        ELSE
          CALL message('fracgeoindex:', 'Beyond northern edge of model')
        END IF
        EXIT
      END IF
    END DO
    IF (klat==-999._dp) THEN
      CALL message('fracgeoindex:', 'Beyond southern edge of model')
    END IF

  END SUBROUTINE fracgeoindex

  SUBROUTINE fracplevindex (colp, p, klev)

    USE mo_exception, ONLY: message
    USE mo_control,   ONLY: nlev

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: colp(:) ! full-level pressures in column
    REAL(dp), INTENT(IN) :: p       ! pressure to find
    REAL(dp), INTENT(OUT) :: klev   ! (output) fractional level in log-pressure

    INTEGER :: jlev

    klev = -999._dp

    DO jlev = 1, nlev
      IF (p <= colp(jlev)) THEN
        IF (jlev > 1) THEN
          klev = jlev - (LOG(colp(jlev)) - LOG(p)) / (LOG(colp(jlev)) - LOG(colp(jlev-1)))
        ELSE
          CALL message('fracplevindex:', 'Beyond top of model')
        END IF
        EXIT
      END IF
    END DO

    IF (klev==-999._dp) THEN
      CALL message('fracgeoindex:', 'Beyond bottom of model')
    END IF

  END SUBROUTINE fracplevindex

  FUNCTION flighttrack_interpolate (fields, ilon, ilev, ilat)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: fields(:,:,:,:)
    REAL(dp), INTENT(IN) :: ilon, ilev, ilat

    REAL(dp) :: flighttrack_interpolate(SIZE(fields,4))

    INTEGER :: lon1, lon2
    INTEGER :: lev1, lev2
    INTEGER :: lat1, lat2

    REAL(dp) :: lonw, levw, latw

    CALL interp_axis(ilon, lon1, lon2, lonw)
    IF (lon1 == 0) lon1 = UBOUND(fields,1)
    IF (lon2 == UBOUND(fields,1)+1) lon2 = 1

    CALL interp_axis(ilev, lev1, lev2, levw)
    CALL interp_axis(ilat, lat1, lat2, latw)

    flighttrack_interpolate =                            &
          (1._dp-latw) * (                               &
            (1._dp-levw) * (                             &
              (1._dp-lonw) * fields(lon1, lev1, lat1, :) &
              +      lonw  * fields(lon2, lev1, lat1, :) &
            )                                            &
            +      levw  * (                             &
              (1._dp-lonw) * fields(lon1, lev2, lat1, :) &
              +      lonw  * fields(lon2, lev2, lat1, :) &
            )                                            &
          )                                              &
          +      latw  * (                               &
            (1._dp-levw) * (                             &
              (1._dp-lonw) * fields(lon1, lev1, lat2, :) &
              +      lonw  * fields(lon2, lev1, lat2, :) &
            )                                            &
            +      levw  * (                             &
              (1._dp-lonw) * fields(lon1, lev2, lat2, :) &
              +      lonw  * fields(lon2, lev2, lat2, :) &
            )                                            &
          )

  END FUNCTION flighttrack_interpolate

  FUNCTION flighttrack_interpolate_horiz(field, ilon, ilat)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: field(:,:,:)
    REAL(dp), INTENT(IN) :: ilon, ilat

    REAL(dp) :: flighttrack_interpolate_horiz(SIZE(field,2))

    INTEGER :: lon1, lon2
    INTEGER :: lat1, lat2

    REAL(dp) :: lonw, latw

    CALL interp_axis(ilon, lon1, lon2, lonw)
    IF (lon1 == 0) lon1 = UBOUND(field,1)
    IF (lon2 == UBOUND(field,1)+1) lon2 = 1

    CALL interp_axis(ilat, lat1, lat2, latw)

    flighttrack_interpolate_horiz =              &
          (1._dp-latw) * (                       &
            (1._dp-lonw) * field(lon1, :, lat1) &
            +      lonw  * field(lon2, :, lat1) &
          )                                      &
          +      latw  * (                       &
            (1._dp-lonw) * field(lon1, :, lat2) &
            +      lonw  * field(lon2, :, lat2) &
          )

  END FUNCTION flighttrack_interpolate_horiz

  SUBROUTINE interp_axis (pos, first, second, ratio)
      
    IMPLICIT NONE
    REAL(dp), INTENT(IN)  :: pos
    INTEGER,  INTENT(OUT) :: first, second
    REAL(dp), INTENT(OUT) :: ratio

    first = FLOOR(pos)
    second = CEILING(pos)
    IF (first /= second) THEN
      ratio = pos - first
    ELSE
      ratio = 0._dp
    END IF

  END SUBROUTINE interp_axis

END MODULE mo_flighttrack
