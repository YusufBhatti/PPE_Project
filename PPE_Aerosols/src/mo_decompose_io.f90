!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_decompose_io

  USE mo_kind,        ONLY: dp
  
  USE mo_memory_base, ONLY: ostreams,         & ! output streams
                            nstreams            ! # elem in ostreams 
  USE mo_exception,   ONLY: message_text, message, finish
  USE mo_linked_list, ONLY: list_element,     & ! linked list element type
                            memory_info         ! stream entry data type

  IMPLICIT NONE
 
  PRIVATE

  PUBLIC :: decompose_io

  PUBLIC :: io_load_balance, io_field_balance

  INTEGER, ALLOCATABLE :: io_load_balance(:)
  INTEGER, ALLOCATABLE :: io_field_balance(:)

  LOGICAL :: do_decompose_io = .true. 

CONTAINS

  !---------------------------------------------------------

  SUBROUTINE loadbalancing(comm_io_size, fields_nlev, fields_mapping)

    INTEGER, INTENT(in)    :: comm_io_size
    INTEGER, INTENT(in)    :: fields_nlev(:)
    INTEGER, INTENT(inout) :: fields_mapping(:)
    
    INTEGER, PARAMETER :: initvalue = 99999
    INTEGER, PARAMETER :: noindex = -1
    INTEGER, PARAMETER :: min = -99
    
    INTEGER, ALLOCATABLE :: fields_sorted(:,:)

    INTEGER :: fields_n, i, j
    INTEGER :: load(comm_io_size)
    INTEGER :: currCapacity, nextCapacity
    INTEGER :: currIndex, nextIndex
    REAL(dp) :: mean, capacityLeft

    fields_n       = SIZE(fields_nlev)
    fields_mapping = initvalue
    load           = 0

    ALLOCATE(fields_sorted(fields_n, 2))
    
    DO i = 1, fields_n
      fields_sorted(i, 1) = fields_nlev(i)
      fields_sorted(i, 2) = i
    END DO
    
    CALL iQsort(fields_sorted)

    mean = REAL(SUM(fields_nlev),dp)/REAL(comm_io_size,dp)

    prbl: DO i = 1, fields_n

      fit: DO j = 1, comm_io_size
        capacityLeft = ABS(mean - REAL(load(j),dp) - REAL(fields_sorted(i,1),dp))
        IF (capacityLeft <= 0.0_dp) THEN
          fields_mapping(fields_sorted(i,2)) = j-1
          load(j) = load(j)+fields_sorted(i,1)
          EXIT fit
        END IF
      END DO fit

      IF (fields_mapping(fields_sorted(i,2)) /= INITVALUE) CYCLE prbl
      
      currCapacity = MIN
      currIndex    = NOINDEX
      
      leastload: DO j = 1,comm_io_size
        nextIndex    = j
        nextCapacity = INT(mean) - load(j)
        IF (nextCapacity > currCapacity) THEN
          currIndex    = nextIndex
          currCapacity = nextCapacity
        END IF
      END DO leastload
      
      fields_mapping(fields_sorted(i,2)) = currIndex-1
      load(currIndex) = load(currIndex)+fields_sorted(i,1)
    END DO prbl
    
    DEALLOCATE(fields_sorted)

  CONTAINS

    RECURSIVE SUBROUTINE iQsort(A)
      
      INTEGER, INTENT(inout) :: A(:,:)
      INTEGER :: pivot(1,2), temp(1,2)
      INTEGER :: n, pivotIndex, currIndex
      
      IF (SIZE(A,2) /=2 ) &
           CALL finish ('iQsort(A)','dimension(A,2) has to be 2.')
      
      n = SIZE(A,1)
      IF (n <= 1) RETURN
      
      pivot(1,:) = A(n,:)
      pivotIndex = n
      currIndex = 1
      
      DO WHILE (currIndex < pivotIndex)
        IF (A(currIndex,1) < pivot(1,1)) THEN
          temp = pivot
          A(pivotIndex,:)   = A(currIndex,:)
          A(currIndex,:)    = A(pivotIndex-1,:)
          A(pivotIndex-1,:) = temp(1,:)
          pivotIndex = pivotIndex-1
        ELSE
          currIndex = currIndex + 1
        END IF
      END DO
      
      IF (pivotIndex /= 1) CALL iQsort(A(1:pivotIndex-1,:))
      IF (pivotIndex /= n) CALL iQsort(A(pivotIndex+1:n,:))
      
    END SUBROUTINE iQsort
    
  END SUBROUTINE loadbalancing

!---------------------------------------------------------

  SUBROUTINE decompose_io(iopes)

    INTEGER, INTENT(in) :: iopes

    TYPE (list_element) ,POINTER :: le  
    TYPE (list_element) ,TARGET  :: first 
    TYPE (memory_info)  ,POINTER :: info

    INTEGER, ALLOCATABLE :: fields_nlev(:), fields_mapping(:)

    INTEGER :: i, fields_n, index

    CHARACTER(len=16) :: cformat

    IF (iopes < 1) RETURN

    IF (.NOT. do_decompose_io) RETURN

    fields_n = 0
    index = 0

    ALLOCATE(io_load_balance(0:iopes-1))
    ALLOCATE(io_field_balance(0:iopes-1))

    io_load_balance(:) = 0
    io_field_balance(:) = 0

    DO i = 1, nstreams
      fields_n = fields_n + ostreams(i)%list_elements
    END DO

    ALLOCATE(fields_nlev(fields_n), fields_mapping(fields_n))

    DO i = 1, nstreams 
      first%next_list_element => ostreams(i)%first_list_element
      le => first
      listcycle1: DO
        le => le%next_list_element
        IF (.NOT.ASSOCIATED(le)) EXIT listcycle1
        info => le%field%info
        index = index + 1
        fields_nlev(index) = info%klev       
      END DO listcycle1
    END DO

    CALL loadbalancing(iopes, fields_nlev, fields_mapping)

    index = 0
    
    DO i = 1, nstreams 
      first%next_list_element => ostreams(i)%first_list_element
      le => first
      listcycle2: DO
        le => le%next_list_element
        IF (.NOT.ASSOCIATED(le)) EXIT listcycle2
        info => le%field%info
        index = index + 1
        info%IO_comm_indx = fields_mapping(index)
        info%IO_stream_indx = i
        io_load_balance(info%IO_comm_indx) = &
             io_load_balance(info%IO_comm_indx) + info%klev  
        io_field_balance(info%IO_comm_indx) = &
             io_field_balance(info%IO_comm_indx) + 1
      END DO listcycle2
    END DO

#ifdef DEBUG
    CALL message ('','')
    WRITE (message_text,'(a,i5)') 'I/O PEs:   ', iopes
    CALL message ('',message_text)
    WRITE (message_text,'(a,i5)') 'I/O fields:', fields_n
    CALL message ('',message_text)
    WRITE (message_text,'(a)') 'I/O field cost:'
    CALL message ('',message_text)
    IF (fields_n <= 12) THEN
      WRITE(cformat,'(a,i0,a)') '(9x,', fields_n, 'i5)'
      WRITE (message_text,TRIM(cformat)) fields_nlev(1:fields_n)
      CALL message ('',message_text)
    ELSE
      WRITE (message_text,'(9x,12i5)') fields_nlev(1:12)
      CALL message ('',message_text)
      DO i = 13, fields_n, 12
        IF (fields_n <= i+11) THEN
          WRITE(cformat,'(a,i0,a)') '(9x,', fields_n, 'i5)'
          WRITE (message_text,TRIM(cformat)) fields_nlev(i:fields_n)
          CALL message ('',message_text)
        ELSE
          WRITE (message_text,'(9x,12i5)') fields_nlev(i:i+11)
          CALL message ('',message_text)
        ENDIF
      ENDDO
    ENDIF
    CALL message ('','')
    WRITE (message_text,'(a)') 'I/O field mapping:'
    CALL message ('',message_text)
    IF (fields_n <= 12) THEN
      WRITE(cformat,'(a,i0,a)') '(9x,', fields_n, 'i5)'
      WRITE (message_text,TRIM(cformat)) fields_mapping(1:fields_n)
      CALL message ('',message_text)
    ELSE
      WRITE (message_text,'(9x,12i5)') fields_mapping(1:12)
      CALL message ('',message_text)
      DO i = 13, fields_n, 12
        IF (fields_n <= i+11) THEN
          WRITE(cformat,'(a,i0,a)') '(9x,', fields_n, 'i5)'
          WRITE (message_text,TRIM(cformat)) fields_mapping(i:fields_n)
          CALL message ('',message_text)
        ELSE
          WRITE (message_text,'(9x,12i5)') fields_mapping(i:i+11)
          CALL message ('',message_text)
        ENDIF
      ENDDO
    ENDIF

    CALL message ('','')
    WRITE (message_text,'(a)') 'I/O load balance:'
    CALL message ('',message_text)
    IF (iopes <= 6) THEN
      WRITE(cformat,'(a,i0,a)') '(9x,', iopes, 'i9)'
      WRITE (message_text,TRIM(cformat)) io_load_balance(0:iopes-1)
      CALL message ('',message_text)
    ELSE
      WRITE (message_text,'(9x,6i9)') io_load_balance(0:5)
      CALL message ('',message_text)
      DO i = 6, iopes, 6
        IF (iopes <= i+5) THEN
          WRITE(cformat,'(a,i0,a)') '(9x,', iopes, 'i9)'
          WRITE (message_text,TRIM(cformat)) io_load_balance(i:iopes-1)
          CALL message ('',message_text)
        ELSE
          WRITE (message_text,'(9x,6i9)') io_load_balance(i:i+5)
          CALL message ('',message_text)
        ENDIF
      ENDDO
    ENDIF
    CALL message ('','')
#endif

    DEALLOCATE(fields_nlev, fields_mapping)

    do_decompose_io = .false.

  END SUBROUTINE decompose_io

!---------------------------------------------------------

END MODULE mo_decompose_io
