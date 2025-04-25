!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
#define DEBUG 1
!>
!! @brief Handles MPI transfer of derived data types
!!
!! @author S. Saarinen, CSC, June 2009, prepared basic version
!! @author L. Kornblueh, D. Kleberg, MPI, July 2009, implemented into 
!!                       echam, expand for level/stream based balancing,
!!                       added network (TCP) logging of I/O server, and
!!                       changed to minimal invasive interface to model
!!                       for generic version for SCALES I/O server library                   
!!
!------------------------------------------------------------------------------
MODULE mo_io_server

  USE mo_kind,             ONLY: dp
  USE mo_exception,        ONLY: message_text, message
  USE mo_util_logging,     ONLY: open_network_log, send_network_log, close_network_log
  USE mo_mpi,              ONLY: p_stop, p_send, p_recv, p_is_iope, p_global_comm,       &
                                 p_pe, p_nprocs, p_num_iopes, p_iopes, p_any_source,     &
                                 p_ioservers, p_parallel_io, p_gio,                      &
                                 p_io_server_info, p_io_info_pair
  USE mo_mpi_types,        ONLY: p_send, p_recv,                                         &
                                 p_define_memory_info_type, p_define_decomposed_type,    &
                                 p_define_stream_type,                                   &
                                 p_free_memory_info_type, p_free_decomposed_type,        &
                                 p_free_stream_type
  USE mo_memory_base,      ONLY: nstreams, ostreams
  USE mo_linked_list,      ONLY: list_element, memory_info, add_list_element,            &
                                 default_stream
  USE mo_decomposition,    ONLY: pe_decomposed, global_decomposition
  USE mo_decompose_io,     ONLY: decompose_io, io_field_balance
  USE mo_util_buffer_pool, ONLY: init_buffer_pool, setup_buffer_pool,                    &
                                 allocate_buffer_pool, buffer_pool_entry,                &
                                 copy_chunk_to_buffer_pool  
  USE mo_output,           ONLY: set_output_time, reset_stream,                          & 
                                 century, year, month, day, hour, minute, second,        &
                                 forecast_hours  
  USE mo_time_control,     ONLY: l_putdata, next_date, write_date


  IMPLICIT NONE

  PUBLIC
 
  INTEGER, PARAMETER :: p_io_server_case_tag          = 2009
  INTEGER, PARAMETER :: p_io_server_nstream_tag       = 2010
  INTEGER, PARAMETER :: p_io_server_ostream_tag       = 2011
  INTEGER, PARAMETER :: p_io_server_decomposition_tag = 2012
  INTEGER, PARAMETER :: p_io_server_nfields_tag       = 2013
  INTEGER, PARAMETER :: p_io_server_field_tag         = 2014

  INTEGER, PARAMETER :: p_io_server_shutdown       =  0
  INTEGER, PARAMETER :: p_io_server_init           =  1
  INTEGER, PARAMETER :: p_io_server_open_streams   =  2
  INTEGER, PARAMETER :: p_io_server_close_streams  =  3
  INTEGER, PARAMETER :: p_io_server_write_streams  =  4
  INTEGER, PARAMETER :: p_io_server_time_set       =  5
  
  INTEGER, PARAMETER :: p_io_server_ncases = 5

  CHARACTER(len=*), PARAMETER :: p_io_server_cases(0:p_io_server_ncases) = (/ &
       'Shutdown server(s)     ', &
       'Initialize             ', &
       'Open Output Streams    ', &
       'Close Output Streams   ', &
       'Write to Output Streams', &
       'Time set               '  /)

  ! keep handlers for RDMA buffers on compute PEs

  INTEGER, ALLOCATABLE :: pd(:)
  
CONTAINS

  !===========================================================================

  SUBROUTINE init_io_server
    CALL p_define_memory_info_type
    CALL p_define_decomposed_type
    CALL p_define_stream_type
  END SUBROUTINE init_io_server

  !===========================================================================

  SUBROUTINE cleanup_io_server
    CALL p_free_memory_info_type  
    CALL p_free_decomposed_type
    CALL p_free_stream_type
  END SUBROUTINE cleanup_io_server

  !===========================================================================

  SUBROUTINE io_server

    LOGICAL  :: lloop, lmaster_iope, lslave_iope, lerror
    INTEGER  :: icase, nslaves, njobs, nfields, iosv, i
    REAL(dp) :: t1, t2, dt
    REAL(dp) :: times(p_io_server_ncases)
    INTEGER  :: ncalls(p_io_server_ncases)

    TYPE(memory_info)           :: recv_buffer    
    TYPE(list_element), POINTER :: new_list_element => NULL()
!#ifdef DEBUG
!    TYPE(list_element), POINTER :: list_entry
!    TYPE(list_element), TARGET  :: first
!    TYPE(memory_info),  POINTER :: info
!#endif
    CHARACTER(len=80) :: logtext

    REAL(dp), EXTERNAL:: util_walltime

    IF (.NOT. p_is_iope()) RETURN

    ! I/O-server from  1 to p_num_iopes

    iosv = p_pe - p_nprocs + 1 

!    CALL open_network_log('creus.mpi.zmaw.de', 19876)

    WRITE(logtext,'(a,i4,a,i4)') 'I/O server: ', iosv, ' global rank: ', p_pe
    CALL send_network_log(TRIM(logtext))

    lmaster_iope = (p_pe == p_nprocs)
    lslave_iope  = (p_pe  > p_nprocs)
    nslaves      = p_num_iopes - 1
    
    times(:)  = 0._dp
    ncalls(:) = 0
    
    njobs = 0
    dt    = 0._dp

    lerror = .FALSE.
    lloop  = .TRUE.

    DO WHILE (lloop)

!      WRITE(logtext,'(a,i4,a)') 'I/O server: ', iosv, ' enter work polling loop.'
      CALL send_network_log(TRIM(logtext))

      CALL p_recv(icase, p_any_source, p_io_server_case_tag, comm=p_global_comm)
      
      IF (icase >= 0 .AND. icase <= p_io_server_ncases) THEN
        WRITE(message_text,'(a,i4,a,:,i12,1x,a)')                  &
             'I/O_server: ', iosv, ' received case:', icase, &
             TRIM(p_io_server_cases(icase))
!        CALL message('io_server', message_text)
        CALL send_network_log(TRIM(message_text))

      ELSE
        WRITE(message_text,'(a,i4,a,:,i12,1x,a)')                  &
             'I/O_server: ', iosv, ' received case:', icase, &
             'Unknown case ...'
!        CALL message('io_server', message_text)  
        CALL send_network_log(TRIM(message_text))      
      ENDIF
      
      t1 = util_walltime()

      SELECT CASE (icase)

      !=======================================================================
      CASE (p_io_server_shutdown)

        lloop = .FALSE.

      !=======================================================================
      CASE (p_io_server_init)

        WRITE(logtext,'(a,i4,a)') 'I/O server: ', iosv, ' initialization ...'
        CALL send_network_log(TRIM(logtext))

        !-- 1. recv no of stream headers 
        CALL p_recv(nstreams,  p_any_source, p_io_server_nstream_tag, comm=p_global_comm)
!#ifdef DEBUG
!        WRITE(logtext,'(a,i4,1x,a,i4)') 'I/O server: ', iosv, ' no of streams: ', nstreams
!        CALL send_network_log(TRIM(logtext))
!#endif        
        !-- 2. recv actual stream header 
        CALL p_recv(ostreams, p_any_source, p_io_server_ostream_tag, count=nstreams, comm=p_global_comm)        
        DO i = 1, nstreams
          ostreams(i) = default_stream
        ENDDO
!#ifdef DEBUG
!        DO i = 1, nstreams
!          WRITE(logtext,'(a,i4,1x,a,a)') 'I/O server: ', iosv, ' stream: ', TRIM(ostreams(i)%name)
!          CALL send_network_log(TRIM(logtext))
!        ENDDO
!#endif        
        !-- 3. recv decomposition
        ALLOCATE(global_decomposition(1:p_nprocs))
        CALL p_recv(global_decomposition,  p_any_source, p_io_server_decomposition_tag,count=p_nprocs,comm=p_global_comm)
        
        !-- 4. send no of field headers per I/O server
        CALL p_recv(nfields, p_any_source, p_io_server_nfields_tag, comm=p_global_comm)
!#ifdef DEBUG
!        WRITE(logtext,'(a,i4,1x,a,i8)') 'I/O server: ', iosv, ' no of fields: ', nfields
!        CALL send_network_log(TRIM(logtext))
!#endif        
        !-- 5. send field headers per I/O server        
        DO i = 1, nfields
          CALL p_recv(recv_buffer, p_any_source, p_io_server_field_tag, comm=p_global_comm)          
!#ifdef DEBUG
!          WRITE(logtext,'(a,i4,1x,a,i6,a,i4)') &
!               'I/O server: ', iosv, ' received field no ', i, ' for stream ', recv_buffer%IO_stream_indx
!          CALL send_network_log(TRIM(logtext))
!#endif
          CALL add_list_element(ostreams(recv_buffer%IO_stream_indx), new_list_element)
          new_list_element%field%info = recv_buffer
          new_list_element%field%ptr => NULL()
        ENDDO
!#ifdef DEBUG        
!        DO i = 1, nstreams
!          first%next_list_element => ostreams(i)%first_list_element
!          list_entry => first
!          cycle_list: DO
!            list_entry => list_entry%next_list_element
!            IF (.NOT.ASSOCIATED(list_entry)) EXIT cycle_list
!            info => list_entry%field%info
!            WRITE(logtext,'(a,i4,1x,a,a)') &
!                 'I/O server: ', iosv, ' field: ', info%name
!            CALL send_network_log(TRIM(logtext))
!          END DO cycle_list
!        END DO
!#endif
        WRITE(logtext,'(a,i4,a)') 'I/O server: ', iosv, ' initialization finished.'
        CALL send_network_log(TRIM(logtext))
        
      !=======================================================================
      CASE (p_io_server_open_streams)
!        CALL open_output_streams
        WRITE(logtext,'(a,i4,a)') 'I/O server: ', iosv, ' open output streams'
        CALL send_network_log(TRIM(logtext))
      !=======================================================================
      CASE (p_io_server_close_streams)
!        CALL close_output_streams
        WRITE(logtext,'(a,i4,a)') 'I/O server: ', iosv, ' close output streams'
        CALL send_network_log(TRIM(logtext))
      !=======================================================================
      CASE (p_io_server_write_streams)
!        CALL out_streams
        WRITE(logtext,'(a,i4,a)') 'I/O server: ', iosv, ' out streams'
        CALL send_network_log(TRIM(logtext))
      !=======================================================================
      CASE (p_io_server_time_set)
! for sure not the right call ;-)
!        CALL time_set
        WRITE(logtext,'(a,i4,a)') 'I/O server: ', iosv, ' output time setting'
        CALL send_network_log(TRIM(logtext))
      CASE default
        lloop  = .FALSE.
        lerror = .TRUE.
      END SELECT

      t2 = util_walltime()
      dt = dt + (t2-t1)
      
      IF (lloop) THEN
        njobs = njobs + 1
        WRITE(message_text,'(a,i4,a,:,i4,1x,a,:,1x,a,f13.9,a)')   &
             'I/O server: ', iosv, ' finished case:', icase, &
             TRIM(p_io_server_cases(icase)), ', time = ', t2-t1, ' s'
!        CALL message('io_server', message_text)  
        CALL send_network_log(TRIM(message_text))      
        IF (icase >= 1 .AND. icase <= p_io_server_ncases) THEN
          times(icase) = times(icase) + (t2-t1)
          ncalls(icase) = ncalls(icase) + 1
        ENDIF
      ENDIF
    ENDDO
    
    IF (lerror) THEN
      WRITE(message_text,'(a,i4,a,i5,1x,f10.3,a)')               &
           'I/O server: ', iosv,                                 &
           ' finished with error(s): # of calls, total time = ', &
           njobs, dt, ' s'
!        CALL message('io_server', message_text)        
        CALL send_network_log(TRIM(message_text))
    ELSE
      WRITE(message_text,'(a,i4,a,i5,1x,f10.3,a)')               &
           'I/O server: ', iosv,                                 &
           ' finished ok: # of calls, total time = ',            &
           njobs, dt, ' s'
!      CALL message('io_server', message_text)    
      CALL send_network_log(TRIM(message_text))    
      DO icase=1,p_io_server_ncases
        WRITE(message_text,'(a,i4,a,i1,a,a30,a,f10.3,a,i5,a)')   &
             'I/O server: ', iosv,                               &
             '(', icase,')', TRIM(p_io_server_cases(icase)),     &
             ': ', times(icase), ' s ', ncalls(icase), ' calls'
      ENDDO
    ENDIF
    
!    WRITE(logtext,'(a,i4,a)') 'I/O server: ', iosv, ' finish ...'
!    CALL send_network_log(TRIM(logtext))

!    CALL close_network_log 

  END SUBROUTINE io_server

  !===========================================================================

  SUBROUTINE io_server_start
    
    TYPE (list_element) ,POINTER    :: list_entry
    TYPE (list_element) ,TARGET     :: first
    TYPE (memory_info)  ,POINTER    :: info

#ifdef DEBUG
    CHARACTER(len=80) :: logtext
#endif
    INTEGER :: i, j, iope, icount, size1, size2, ndim, istat

    IF (p_num_iopes == 0 .OR. p_is_iope()) RETURN

    IF (p_ioservers) CALL io_wakeup(p_io_server_init)
    
    !-- decompose I/O

    CALL decompose_io(p_num_iopes)

    IF (p_pe ==  p_io_info_pair(p_pe)%comp_source .AND. p_pe < p_num_iopes) THEN

      iope = p_io_info_pair(p_pe)%io_target
!#ifdef DEBUG
!      WRITE(logtext,'(a,i4,1x,a,i4)') 'Compute server (source): ', p_pe, ' I/O server (target): ', iope
!      CALL send_network_log(TRIM(logtext))
!#endif        

      !-- 1. send no of stream headers 
      CALL p_send(nstreams, iope, p_io_server_nstream_tag, comm=p_global_comm)

      !-- 2. send actual stream header 
      CALL p_send(ostreams, iope, p_io_server_ostream_tag,count=nstreams, comm=p_global_comm)

      !-- 3. send decomposition
      CALL p_send(global_decomposition, iope, p_io_server_decomposition_tag, count=p_nprocs, comm=p_global_comm)

      !-- 4. send no of field headers per I/O server
      CALL p_send(SUM(io_field_balance), iope, p_io_server_nfields_tag, comm=p_global_comm)
     
      !-- 5. send all field headers to all I/O server
      icount = 0
      DO i = 1, nstreams
        first%next_list_element => ostreams(i)%first_list_element
        list_entry => first
        cycle_list: DO
          list_entry => list_entry%next_list_element
          IF (.NOT.ASSOCIATED(list_entry)) EXIT cycle_list
          info => list_entry%field%info
!#ifdef DEBUG
!          WRITE(logtext,'(a,i4,1x,a,i4,a,a,a,i5)') &
!               'Compute server: ', p_pe, 'stream ', i, ': ', info%name, ' count = ', icount
!          CALL send_network_log(TRIM(logtext))
!#endif
          CALL p_send(info, iope, p_io_server_field_tag, comm=p_global_comm)          
          icount = icount+1
        END DO cycle_list
      END DO

    ENDIF

    !-- initialize RDMA buffer 
      
    IF (.NOT. ALLOCATED(pd)) THEN
      ALLOCATE(pd(p_num_iopes))
    ENDIF

    DO iope = 1, p_num_iopes
      WRITE(logtext,'(a,i3,a,i3,a,i4)') &
           'Compute node: ', p_pe, &
           'handler: ', iope, &
           ' pool entries: ', io_field_balance(iope-1)
      CALL send_network_log(TRIM(logtext))
      pd(iope) = init_buffer_pool(io_field_balance(iope-1))
    ENDDO
 
    DO i = 1, nstreams
      first%next_list_element => ostreams(i)%first_list_element
      list_entry => first
      cycle_buffer_list: DO
        list_entry => list_entry%next_list_element
        IF (.NOT.ASSOCIATED(list_entry)) EXIT cycle_buffer_list
        info => list_entry%field%info
        !
        !-- check arguments:
        !   1) size 
        size1 = 1
        DO j = 1, 4
          size1 = size1*list_entry%field%info%dima(j)
        ENDDO
        size2 = SIZE(list_entry%field%ptr)
        IF (size1 /= size2) THEN
          WRITE(logtext,'(a,6i7)') 'Field size specification', &
               size1, size2, list_entry%field%info%dima(1:4)
          CALL send_network_log(TRIM(logtext))
        ENDIF
        !   2) dimensions
        ndim = list_entry%field%info%ndim 
        !   3) target I/O PE
        iope = list_entry%field%info%IO_comm_indx+1
        !   4) register
        istat = buffer_pool_entry(pd(iope), size1, ndim, list_entry%field%ptr) 
        !
        !--
        icount = icount+1
      END DO cycle_buffer_list
    END DO
    
    DO iope = 1, p_num_iopes
      istat = allocate_buffer_pool(pd(iope), 1);
      istat = setup_buffer_pool(pd(iope))
    ENDDO
    
  END SUBROUTINE io_server_start

  !===========================================================================

  SUBROUTINE io_server_stop

    IF (p_num_iopes == 0 .OR. p_is_iope()) RETURN

    IF (p_ioservers) CALL io_wakeup(p_io_server_shutdown)

  END SUBROUTINE io_server_stop

  !===========================================================================

  SUBROUTINE io_server_open_streams

    IF (p_num_iopes == 0 .OR. p_is_iope()) RETURN

    IF (p_ioservers) CALL io_wakeup(p_io_server_open_streams)
    
  END SUBROUTINE io_server_open_streams

  !===========================================================================

  SUBROUTINE io_server_write_streams

    TYPE (list_element) ,POINTER    :: list_entry
    TYPE (list_element) ,TARGET     :: first

    REAL(dp) :: t1, t2

    INTEGER, ALLOCATABLE :: icount(:)
    INTEGER :: i, iope, istat

    REAL(dp), EXTERNAL:: util_walltime

    LOGICAL :: write_info

    IF (p_num_iopes == 0 .OR. p_is_iope()) RETURN

    t1 = util_walltime()

    write_info = .FALSE.

    IF (p_ioservers) CALL io_wakeup(p_io_server_write_streams)

    ALLOCATE(icount(p_num_iopes))
    icount(:) = 0
    DO i = 1, nstreams
      IF(l_putdata(ostreams(i)%post_idx)) THEN
        write_info = .TRUE.
        first%next_list_element => ostreams(i)%first_list_element
        list_entry => first
        cycle_list: DO
          list_entry => list_entry%next_list_element
          IF (.NOT.ASSOCIATED(list_entry)) EXIT cycle_list
          iope = list_entry%field%info%IO_comm_indx+1
          istat = copy_chunk_to_buffer_pool(pd(iope), icount(iope), list_entry%field%ptr)
          icount(iope) = icount(iope)+1
        ENDDO cycle_list
        CALL reset_stream(ostreams(i))
      ENDIF
    ENDDO

    t2 = util_walltime()
    
    IF (write_info) THEN
      CALL write_date(next_date,'Copy output data for : ')
    ENDIF

  END SUBROUTINE io_server_write_streams

  !===========================================================================

  SUBROUTINE io_server_close_streams

    IF (p_num_iopes == 0 .OR. p_is_iope()) RETURN

    IF (p_ioservers) CALL io_wakeup(p_io_server_close_streams)

  END SUBROUTINE io_server_close_streams

  !===========================================================================

  SUBROUTINE io_wakeup(icase, iope)

    INTEGER,           INTENT(in) :: icase
    INTEGER, OPTIONAL, INTENT(in) :: iope

    INTEGER :: jp

    IF (p_num_iopes == 0 .OR. p_is_iope()) RETURN

    IF (icase >= 0 .AND. icase <= p_io_server_ncases) THEN
      WRITE(message_text,'(i5,a,i12,:,1x,a)')             &
           p_pe, ': io_wakeup() : sending case=', icase,  &
           TRIM(p_io_server_cases(icase))
    ELSE
      WRITE(message_text,'(i5,a,i12,:,1x,a)')            &
           p_pe, ': io_wakeup() : sending case=', icase, &
           '(unknown case ...)'
    ENDIF

    IF (PRESENT(iope)) THEN
      CALL p_send(icase, iope, p_io_server_case_tag, comm=p_global_comm)
    ELSE
      DO jp = 0, p_num_iopes-1
        IF (p_pe == p_io_info_pair(jp)%comp_source) &
             CALL p_send(icase, p_io_info_pair(jp)%io_target, p_io_server_case_tag, comm=p_global_comm)
      ENDDO
    ENDIF
    
  END SUBROUTINE io_wakeup
  
  !===========================================================================

  SUBROUTINE io_shutdown

    IF (p_num_iopes == 0 .OR. p_is_iope()) RETURN
    CALL io_wakeup(p_io_server_shutdown, -1)

  END SUBROUTINE io_shutdown

  !===========================================================================

  SUBROUTINE io_set_output_time

    INTEGER :: buffer(8)  ! century, year, month, day, hour, minute, second, forecat_hours
    INTEGER :: jp

    IF (p_is_iope()) THEN
      CALL p_recv(buffer, p_any_source, p_io_server_time_set, comm=p_global_comm)
      century        = buffer(1)
      year           = buffer(2)
      month          = buffer(3)
      day            = buffer(4)
      hour           = buffer(5)
      minute         = buffer(6)
      second         = buffer(7)
      forecast_hours = buffer(8)
    ELSE
      CALL set_output_time
      buffer(:) = (/ century, year, month, day, hour, minute, second, forecast_hours /)
      DO jp = p_nprocs, p_nprocs+p_num_iopes-1
        CALL p_send(buffer, p_iopes(jp), p_io_server_time_set, comm=p_global_comm)
      ENDDO
    ENDIF
  END SUBROUTINE io_set_output_time
  
  !===========================================================================

END MODULE mo_io_server

