!>
!! Configuration of data processing streams (mvstreams)
!!
!! The term 'mvstream' may refer to Mean Value Stream or Meta Value Stream.
!!
NAMELIST /mvstreamctl/ &
  target, &       ! name of new stream ['*m']
  interval, &   ! time interval over which reduction is effected [1, 'months', 'last', 0]
  filetag, &    ! name tag for output file ['<target>']
  source, &     ! name of stream for which data reduction is requested
  meannam, &    ! list of variables for which means are calculated [all applicable source variables]
  sqrmeannam, & ! list of variables for which means of the squared values are calculated [none]
  ! Deprecated fields
  m_stream_name ! names of streams for the elements of which processing may be desired
