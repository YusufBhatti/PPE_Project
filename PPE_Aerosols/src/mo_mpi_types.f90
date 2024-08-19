!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Handles MPI transfer of derived data types
!!
!! @author D. Kleberg, MPI, July 2009, prepared basic version
!! @author L. Kornblueh, MPI, July 2009, added more datatypes
!!
!! @note
!! For using Fortran derived datatypes as MPI datatypes the datatype
!! should have the SEQUENCE attribute. 
!! 
!! @par
!! To assure the proper send/recv of vectors of MPI datatypes the last
!! element in a Fortran derived datatype must be one of the prefedined 
!! Fortran and MPI datatypes. Therefore we include an INTEGER :: last
!! in all respective datatypes.
!!
!
MODULE mo_mpi_types

#ifndef NOMPI 
  USE mpi
#endif

#ifdef DEBUG
  USE mo_mpi,           ONLY: p_real_dp, p_int_i8, p_pe, p_global_comm 
  USE mo_exception,     ONLY: message, message_text, finish
#else
  USE mo_mpi,           ONLY: p_real_dp, p_int_i8, p_global_comm 
#endif
  USE mo_linked_list,   ONLY: memory_info, t_stream
  USE mo_decomposition, ONLY: pe_decomposed
  
  IMPLICIT NONE

  PRIVATE  
  
  INTEGER :: p_memory_info_type
  INTEGER :: p_decomposed_type
  INTEGER :: p_stream_type

  LOGICAL :: p_memory_info_defined  = .FALSE.
  LOGICAL :: p_decomposed_defined   = .FALSE.
  LOGICAL :: p_stream_defined       = .FALSE.


  PUBLIC :: p_define_memory_info_type
  PUBLIC :: p_define_decomposed_type
  PUBLIC :: p_define_stream_type

  PUBLIC :: p_free_memory_info_type  
  PUBLIC :: p_free_decomposed_type
  PUBLIC :: p_free_stream_type

  PUBLIC :: p_send
  PUBLIC :: p_recv

  INTERFACE p_send
    MODULE PROCEDURE p_send_memory_info
    MODULE PROCEDURE p_send_memory_info_1d
    MODULE PROCEDURE p_send_decomposed
    MODULE PROCEDURE p_send_decomposed_1d
    MODULE PROCEDURE p_send_stream
    MODULE PROCEDURE p_send_stream_1d
  END INTERFACE

  INTERFACE p_recv
    MODULE PROCEDURE p_recv_memory_info
    MODULE PROCEDURE p_recv_memory_info_1d
    MODULE PROCEDURE p_recv_decomposed
    MODULE PROCEDURE p_recv_decomposed_1d
    MODULE PROCEDURE p_recv_stream
    MODULE PROCEDURE p_recv_stream_1d
  END INTERFACE

  INTEGER :: p_error = 0

CONTAINS
  
  !---------------------------------------------------------------------------
  
  SUBROUTINE p_define_memory_info_type

#ifndef NOMPI 
    INTEGER, PARAMETER :: count = 36

    INTEGER                   :: blockcounts (count)
    INTEGER                   ::    mpitypes (count)
    INTEGER(MPI_ADDRESS_KIND) ::       displ (count)

    INTEGER(MPI_ADDRESS_KIND) :: base

    TYPE(memory_info) :: dummy

    CALL MPI_GET_ADDRESS (dummy%name          , displ( 1), p_error)          
    CALL MPI_GET_ADDRESS (dummy%units         , displ( 2), p_error)
    CALL MPI_GET_ADDRESS (dummy%longname      , displ( 3), p_error) 
    CALL MPI_GET_ADDRESS (dummy%gdim          , displ( 4), p_error) 
    CALL MPI_GET_ADDRESS (dummy%dim           , displ( 5), p_error)
    CALL MPI_GET_ADDRESS (dummy%dima          , displ( 6), p_error)
    CALL MPI_GET_ADDRESS (dummy%lreg          , displ( 7), p_error)
    CALL MPI_GET_ADDRESS (dummy%ndim          , displ( 8), p_error)
    CALL MPI_GET_ADDRESS (dummy%klev          , displ( 9), p_error)
    CALL MPI_GET_ADDRESS (dummy%iklev         , displ(10), p_error)
    CALL MPI_GET_ADDRESS (dummy%alloc         , displ(11), p_error)
    CALL MPI_GET_ADDRESS (dummy%repr          , displ(12), p_error)
    CALL MPI_GET_ADDRESS (dummy%lpost         , displ(13), p_error)
    CALL MPI_GET_ADDRESS (dummy%laccu         , displ(14), p_error)
    CALL MPI_GET_ADDRESS (dummy%lmiss         , displ(15), p_error)
    CALL MPI_GET_ADDRESS (dummy%missval       , displ(16), p_error)
    CALL MPI_GET_ADDRESS (dummy%reset         , displ(17), p_error)
    CALL MPI_GET_ADDRESS (dummy%lrerun        , displ(18), p_error)
    CALL MPI_GET_ADDRESS (dummy%contnorest    , displ(19), p_error)
    CALL MPI_GET_ADDRESS (dummy%restart_read  , displ(20), p_error)
    CALL MPI_GET_ADDRESS (dummy%restart_pref  , displ(21), p_error)
    CALL MPI_GET_ADDRESS (dummy%gribtable     , displ(22), p_error)
    CALL MPI_GET_ADDRESS (dummy%gribcode      , displ(23), p_error)
    CALL MPI_GET_ADDRESS (dummy%gribbits      , displ(24), p_error)
    CALL MPI_GET_ADDRESS (dummy%levelindx     , displ(25), p_error)
    CALL MPI_GET_ADDRESS (dummy%tracidx       , displ(26), p_error)
    CALL MPI_GET_ADDRESS (dummy%IO_var_indx   , displ(27), p_error)
    CALL MPI_GET_ADDRESS (dummy%IO_var_id     , displ(28), p_error)
    CALL MPI_GET_ADDRESS (dummy%IO_var_stid   , displ(29), p_error)
    CALL MPI_GET_ADDRESS (dummy%IO_name       , displ(30), p_error)
    CALL MPI_GET_ADDRESS (dummy%IO_unit       , displ(31), p_error)
    CALL MPI_GET_ADDRESS (dummy%gridID        , displ(32), p_error)
    CALL MPI_GET_ADDRESS (dummy%zaxisID       , displ(33), p_error)
    CALL MPI_GET_ADDRESS (dummy%IO_comm_indx  , displ(34), p_error)
    CALL MPI_GET_ADDRESS (dummy%IO_stream_indx, displ(35), p_error)
    CALL MPI_GET_ADDRESS (dummy%last          , displ(36), p_error)

    base = displ(1)
    displ(:) = displ(:) - base

    mpitypes( 1) = MPI_CHARACTER
    mpitypes( 2) = MPI_CHARACTER
    mpitypes( 3) = MPI_CHARACTER
    mpitypes( 4) = MPI_INTEGER
    mpitypes( 5) = MPI_INTEGER
    mpitypes( 6) = MPI_INTEGER
    mpitypes( 7) = MPI_LOGICAL
    mpitypes( 8) = MPI_INTEGER
    mpitypes( 9) = MPI_INTEGER
    mpitypes(10) = MPI_INTEGER
    mpitypes(11) = MPI_LOGICAL
    mpitypes(12) = MPI_INTEGER
    mpitypes(13) = MPI_LOGICAL
    mpitypes(14) = MPI_LOGICAL
    mpitypes(15) = MPI_LOGICAL
    mpitypes(16) = p_real_dp
    mpitypes(17) = p_real_dp
    mpitypes(18) = MPI_LOGICAL
    mpitypes(19) = MPI_LOGICAL
    mpitypes(20) = MPI_LOGICAL
    mpitypes(21) = MPI_CHARACTER
    mpitypes(22) = MPI_INTEGER
    mpitypes(23) = MPI_INTEGER
    mpitypes(24) = MPI_INTEGER
    mpitypes(25) = MPI_INTEGER
    mpitypes(26) = MPI_INTEGER
    mpitypes(27) = MPI_INTEGER
    mpitypes(28) = MPI_INTEGER
    mpitypes(29) = MPI_INTEGER
    mpitypes(30) = MPI_CHARACTER
    mpitypes(31) = MPI_CHARACTER
    mpitypes(32) = MPI_INTEGER
    mpitypes(33) = MPI_INTEGER
    mpitypes(34) = MPI_INTEGER
    mpitypes(35) = MPI_INTEGER
    mpitypes(36) = MPI_INTEGER

    blockcounts(:) = 1

    blockcounts( 1) = LEN(dummy%name) 
    blockcounts( 2) = LEN(dummy%units)
    blockcounts( 3) = LEN(dummy%longname)
    blockcounts( 4) = SIZE(dummy%gdim)
    blockcounts( 5) = SIZE(dummy%dim)
    blockcounts( 6) = SIZE(dummy%dima)
    blockcounts(27) = SIZE(dummy%IO_var_indx)
    blockcounts(30) = LEN(dummy%IO_name)
    blockcounts(31) = LEN(dummy%IO_unit)

    CALL MPI_TYPE_CREATE_STRUCT(count, blockcounts, displ, mpitypes, p_memory_info_type, p_error)
    CALL MPI_TYPE_COMMIT(p_memory_info_type, p_error)

    p_memory_info_defined = .TRUE. 
#endif

  END SUBROUTINE p_define_memory_info_type

  !---------------------------------------------------------------------------

  SUBROUTINE p_free_memory_info_type

#ifndef NOMPI 
    IF (p_memory_info_defined) CALL MPI_TYPE_FREE(p_memory_info_type, p_error)

    p_memory_info_defined = .FALSE.
#endif

  END SUBROUTINE p_free_memory_info_type

  !---------------------------------------------------------------------------

  SUBROUTINE p_define_decomposed_type

#ifndef NOMPI 
    INTEGER, PARAMETER :: count = 55

    INTEGER                   :: blockcounts (count)
    INTEGER                   ::    mpitypes (count)
    INTEGER(MPI_ADDRESS_KIND) ::       displ (count)

    INTEGER(MPI_ADDRESS_KIND) :: base

    TYPE(pe_decomposed) :: dummy

    CALL MPI_GET_ADDRESS (dummy%nlon    , displ( 1), p_error)
    CALL MPI_GET_ADDRESS (dummy%nlat    , displ( 2), p_error)
    CALL MPI_GET_ADDRESS (dummy%npts    , displ( 3), p_error)
    CALL MPI_GET_ADDRESS (dummy%nlev    , displ( 4), p_error)
    CALL MPI_GET_ADDRESS (dummy%nm      , displ( 5), p_error)   
   
    CALL MPI_GET_ADDRESS (dummy%d_nprocs, displ( 6), p_error)
    CALL MPI_GET_ADDRESS (dummy%spe     , displ( 7), p_error)
    CALL MPI_GET_ADDRESS (dummy%epe     , displ( 8), p_error)
    CALL MPI_GET_ADDRESS (dummy%nprocb  , displ( 9), p_error)
    CALL MPI_GET_ADDRESS (dummy%nproca  , displ(10), p_error)     
    
    CALL MPI_GET_ADDRESS (dummy%pe      , displ(11), p_error)
    CALL MPI_GET_ADDRESS (dummy%set_b   , displ(12), p_error)
    CALL MPI_GET_ADDRESS (dummy%set_a   , displ(13), p_error)
    CALL MPI_GET_ADDRESS (dummy%col_1d  , displ(14), p_error)       
   
    CALL MPI_GET_ADDRESS (dummy%nglat   , displ(15), p_error)
    CALL MPI_GET_ADDRESS (dummy%nglatmax, displ(16), p_error)
    CALL MPI_GET_ADDRESS (dummy%nglon   , displ(17), p_error)
    CALL MPI_GET_ADDRESS (dummy%nglonmax, displ(18), p_error)  
    CALL MPI_GET_ADDRESS (dummy%nglh    , displ(19), p_error)
    CALL MPI_GET_ADDRESS (dummy%glats   , displ(20), p_error)
    CALL MPI_GET_ADDRESS (dummy%glate   , displ(21), p_error)
    CALL MPI_GET_ADDRESS (dummy%glons   , displ(22), p_error)
    CALL MPI_GET_ADDRESS (dummy%glone   , displ(23), p_error)
  
    CALL MPI_GET_ADDRESS (dummy%ngpts   , displ(24), p_error)
    CALL MPI_GET_ADDRESS (dummy%gptss   , displ(25), p_error)
    CALL MPI_GET_ADDRESS (dummy%gptse   , displ(26), p_error) 
  
    CALL MPI_GET_ADDRESS (dummy%ngpblks , displ(27), p_error)
    CALL MPI_GET_ADDRESS (dummy%nproma  , displ(28), p_error)
    CALL MPI_GET_ADDRESS (dummy%npromz  , displ(29), p_error)
    CALL MPI_GET_ADDRESS (dummy%lreg    , displ(30), p_error)
  
    CALL MPI_GET_ADDRESS (dummy%lfused  , displ(31), p_error)
    CALL MPI_GET_ADDRESS (dummy%nflat   , displ(32), p_error)
    CALL MPI_GET_ADDRESS (dummy%nflev   , displ(33), p_error)
    CALL MPI_GET_ADDRESS (dummy%nflevp1 , displ(34), p_error)
    CALL MPI_GET_ADDRESS (dummy%flats   , displ(35), p_error)
    CALL MPI_GET_ADDRESS (dummy%flate   , displ(36), p_error)
    CALL MPI_GET_ADDRESS (dummy%flevs   , displ(37), p_error)
    CALL MPI_GET_ADDRESS (dummy%fleve   , displ(38), p_error)
   
    CALL MPI_GET_ADDRESS (dummy%nlm     , displ(39), p_error)
    CALL MPI_GET_ADDRESS (dummy%lnsp    , displ(40), p_error)
    CALL MPI_GET_ADDRESS (dummy%nlnm0   , displ(41), p_error)
  
    CALL MPI_GET_ADDRESS (dummy%nllev   , displ(42), p_error)
    CALL MPI_GET_ADDRESS (dummy%nllevp1 , displ(43), p_error)
    CALL MPI_GET_ADDRESS (dummy%llevs   , displ(44), p_error)
    CALL MPI_GET_ADDRESS (dummy%lleve   , displ(45), p_error)
   
    CALL MPI_GET_ADDRESS (dummy%snsp    , displ(46), p_error)
    CALL MPI_GET_ADDRESS (dummy%snsp2   , displ(47), p_error)
    CALL MPI_GET_ADDRESS (dummy%ssps    , displ(48), p_error)
    CALL MPI_GET_ADDRESS (dummy%sspe    , displ(49), p_error)
  
    CALL MPI_GET_ADDRESS (dummy%lfirstc , displ(50), p_error)
    CALL MPI_GET_ADDRESS (dummy%ifirstc , displ(51), p_error)
   
    CALL MPI_GET_ADDRESS (dummy%nns     , displ(52), p_error)
  
    CALL MPI_GET_ADDRESS (dummy%nsm     , displ(53), p_error)
    CALL MPI_GET_ADDRESS (dummy%nsnm0   , displ(54), p_error)
   
    CALL MPI_GET_ADDRESS (dummy%last    , displ(55), p_error)

    base = displ(1)
    displ(:) = displ(:) - base

    mpitypes(1)  = MPI_INTEGER
    mpitypes(2)  = MPI_INTEGER
    mpitypes(3)  = MPI_INTEGER
    mpitypes(4)  = MPI_INTEGER
    mpitypes(5)  = MPI_INTEGER  
  
    mpitypes(6)  = MPI_INTEGER
    mpitypes(7)  = MPI_INTEGER
    mpitypes(8)  = MPI_INTEGER
    mpitypes(9)  = MPI_INTEGER
    mpitypes(10) = MPI_INTEGER 

    mpitypes(11) = MPI_INTEGER
    mpitypes(12) = MPI_INTEGER
    mpitypes(13) = MPI_INTEGER
    mpitypes(14) = MPI_LOGICAL  

    mpitypes(15) = MPI_INTEGER
    mpitypes(16) = MPI_INTEGER
    mpitypes(17) = MPI_INTEGER
    mpitypes(18) = MPI_INTEGER  
    mpitypes(19) = MPI_INTEGER
    mpitypes(20) = MPI_INTEGER
    mpitypes(21) = MPI_INTEGER
    mpitypes(22) = MPI_INTEGER
    mpitypes(23) = MPI_INTEGER
  
    mpitypes(24) = MPI_INTEGER
    mpitypes(25) = MPI_INTEGER
    mpitypes(26) = MPI_INTEGER
 
    mpitypes(27) = MPI_INTEGER
    mpitypes(28) = MPI_INTEGER
    mpitypes(29) = MPI_INTEGER
    mpitypes(30) = MPI_LOGICAL

    mpitypes(31) = MPI_LOGICAL
    mpitypes(32) = MPI_INTEGER
    mpitypes(33) = MPI_INTEGER
    mpitypes(34) = MPI_INTEGER
    mpitypes(35) = MPI_INTEGER
    mpitypes(36) = MPI_INTEGER
    mpitypes(37) = MPI_INTEGER
    mpitypes(38) = MPI_INTEGER

    mpitypes(39) = MPI_INTEGER
    mpitypes(40) = MPI_INTEGER
    mpitypes(41) = MPI_INTEGER

    mpitypes(42) = MPI_INTEGER
    mpitypes(43) = MPI_INTEGER
    mpitypes(44) = MPI_INTEGER
    mpitypes(45) = MPI_INTEGER

    mpitypes(46) = MPI_INTEGER
    mpitypes(47) = MPI_INTEGER
    mpitypes(48) = MPI_INTEGER
    mpitypes(49) = MPI_INTEGER

    mpitypes(50) = MPI_LOGICAL
    mpitypes(51) = MPI_INTEGER

    mpitypes(52) = MPI_INTEGER

    mpitypes(53) = MPI_INTEGER
    mpitypes(54) = MPI_INTEGER
    
    mpitypes(55) = MPI_INTEGER

    blockcounts(:) = 1
    blockcounts(19) = SIZE(dummy%nglh)
    blockcounts(20) = SIZE(dummy%glats)
    blockcounts(21) = SIZE(dummy%glate)
    blockcounts(22) = SIZE(dummy%glons)
    blockcounts(23) = SIZE(dummy%glone)
    blockcounts(35) = SIZE(dummy%flats)
    blockcounts(36) = SIZE(dummy%flate)

    CALL MPI_TYPE_CREATE_STRUCT(count, blockcounts, displ, mpitypes, p_decomposed_type, p_error)
    CALL MPI_TYPE_COMMIT(p_decomposed_type, p_error)

    p_decomposed_defined = .TRUE. 
#endif

  END SUBROUTINE p_define_decomposed_type

  !---------------------------------------------------------------------------

  SUBROUTINE p_free_decomposed_type

#ifndef NOMPI 
    IF (p_decomposed_defined) CALL MPI_TYPE_FREE(p_decomposed_type, p_error)

    p_decomposed_defined = .FALSE.
#endif
    
  END SUBROUTINE p_free_decomposed_type

  !---------------------------------------------------------------------------

  SUBROUTINE p_define_stream_type

#ifndef NOMPI 
    INTEGER, PARAMETER :: count = 19

    INTEGER                   :: blockcounts(count)
    INTEGER                   ::    mpitypes(count)
    INTEGER(MPI_ADDRESS_KIND) ::       displ(count)

    INTEGER(MPI_ADDRESS_KIND) :: base

    TYPE(t_stream) :: dummy

    CALL MPI_GET_ADDRESS(dummy%name         , displ( 1) , p_error)
    CALL MPI_GET_ADDRESS(dummy%memory_used  , displ( 2) , p_error)
    CALL MPI_GET_ADDRESS(dummy%list_elements, displ( 3) , p_error)
    CALL MPI_GET_ADDRESS(dummy%lpost        , displ( 4) , p_error)
    CALL MPI_GET_ADDRESS(dummy%lpout        , displ( 5) , p_error)
    CALL MPI_GET_ADDRESS(dummy%lrerun       , displ( 6) , p_error)
    CALL MPI_GET_ADDRESS(dummy%linit        , displ( 7) , p_error)
    CALL MPI_GET_ADDRESS(dummy%filename     , displ( 8) , p_error)
    CALL MPI_GET_ADDRESS(dummy%post_suf     , displ( 9) , p_error)
    CALL MPI_GET_ADDRESS(dummy%rest_suf     , displ(10) , p_error)
    CALL MPI_GET_ADDRESS(dummy%init_suf     , displ(11) , p_error)
    CALL MPI_GET_ADDRESS(dummy%first        , displ(12) , p_error)
    CALL MPI_GET_ADDRESS(dummy%post_idx     , displ(13) , p_error)
    CALL MPI_GET_ADDRESS(dummy%filetype     , displ(14) , p_error)
    CALL MPI_GET_ADDRESS(dummy%ztype        , displ(15) , p_error)
    CALL MPI_GET_ADDRESS(dummy%fileID       , displ(16) , p_error)
    CALL MPI_GET_ADDRESS(dummy%vlistID      , displ(17) , p_error)
    CALL MPI_GET_ADDRESS(dummy%timestep     , displ(18) , p_error)
    CALL MPI_GET_ADDRESS(dummy%last         , displ(19) , p_error)
    
    base = displ(1)
    displ(:) = displ(:) - base
    
    mpitypes( 1) = MPI_CHARACTER
    mpitypes( 2) = p_int_i8
    mpitypes( 3) = MPI_INTEGER
    mpitypes( 4) = MPI_LOGICAL
    mpitypes( 5) = MPI_LOGICAL
    mpitypes( 6) = MPI_LOGICAL
    mpitypes( 7) = MPI_LOGICAL
    mpitypes( 8) = MPI_CHARACTER
    mpitypes( 9) = MPI_CHARACTER
    mpitypes(10) = MPI_CHARACTER
    mpitypes(11) = MPI_CHARACTER
    mpitypes(12) = MPI_LOGICAL
    mpitypes(13) = MPI_INTEGER
    mpitypes(14) = MPI_INTEGER
    mpitypes(15) = MPI_INTEGER
    mpitypes(16) = MPI_INTEGER
    mpitypes(17) = MPI_INTEGER
    mpitypes(18) = MPI_INTEGER
    mpitypes(19) = MPI_INTEGER
    
    blockcounts(:) = 1
    blockcounts( 1) = LEN(dummy%name)
    blockcounts( 8) = LEN(dummy%filename)
    blockcounts( 9) = LEN(dummy%post_suf)
    blockcounts(10) = LEN(dummy%rest_suf)
    blockcounts(11) = LEN(dummy%init_suf)
    
    CALL MPI_TYPE_CREATE_STRUCT(count, blockcounts, displ, mpitypes, p_stream_type, p_error)
    CALL MPI_TYPE_COMMIT(p_stream_type, p_error)

    p_stream_defined = .true.
#endif

  END SUBROUTINE p_define_stream_type

  !---------------------------------------------------------------------------

  SUBROUTINE p_free_stream_type

#ifndef NOMPI 
    IF (p_stream_defined) CALL MPI_TYPE_FREE(p_stream_type, p_error)
    
    p_stream_defined = .FALSE.
#endif    

  END SUBROUTINE p_free_stream_type

  !---------------------------------------------------------------------------

  SUBROUTINE p_send_memory_info (buffer, dest, tag, count, comm)
      
    TYPE(memory_info), INTENT(in) :: buffer
    INTEGER          , INTENT(in) :: dest, tag
    INTEGER, OPTIONAL, INTENT(in) :: count, comm
    
#ifndef NOMPI
    INTEGER :: local_comm
      
    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_SEND (buffer, count, p_memory_info_type, dest, tag, &
           local_comm, p_error)
    ELSE
      CALL MPI_SEND (buffer, 1, p_memory_info_type, dest, tag, &
           local_comm, p_error)
    ENDIF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
           ' to ', dest, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)
    END IF
#endif
#endif   

  END SUBROUTINE p_send_memory_info

  !---------------------------------------------------------------------------

  SUBROUTINE p_send_memory_info_1d (buffer, dest, tag, count, comm)
      
    TYPE(memory_info),     INTENT(in) :: buffer(:)
    INTEGER,               INTENT(in) :: dest, tag
    INTEGER, OPTIONAL,     INTENT(in) :: count, comm

#ifndef NOMPI
    INTEGER :: local_comm
    
    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
    
    IF (PRESENT(count)) THEN
      CALL MPI_SEND (buffer, count, p_memory_info_type, dest, tag, &
           local_comm, p_error)
    ELSE
      CALL MPI_SEND (buffer, SIZE(buffer), p_memory_info_type, dest, tag, &
           local_comm, p_error)
    END IF
    
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
           ' to ', dest, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)
    END IF
#endif
#endif
    
  END SUBROUTINE p_send_memory_info_1d

  !---------------------------------------------------------------------------

  SUBROUTINE p_send_decomposed (buffer, dest, tag, count, comm)
      
    TYPE(pe_decomposed), INTENT(in) :: buffer
    INTEGER            , INTENT(in) :: dest, tag
    INTEGER, OPTIONAL  , INTENT(in) :: count, comm
    
#ifndef NOMPI
    INTEGER :: local_comm
      
    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_SEND (buffer, count, p_decomposed_type, dest, tag, &
           local_comm, p_error)
    ELSE
      CALL MPI_SEND (buffer, 1, p_decomposed_type, dest, tag, &
           local_comm, p_error)
    END IF
      
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
           ' to ', dest, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)
    END IF
#endif
#endif   

  END SUBROUTINE p_send_decomposed

  !---------------------------------------------------------------------------

  SUBROUTINE p_send_decomposed_1d (buffer, dest, tag, count, comm)
      
    TYPE(pe_decomposed), INTENT(in) :: buffer(:)
    INTEGER            , INTENT(in) :: dest, tag
    INTEGER, OPTIONAL  , INTENT(in) :: count, comm
    
#ifndef NOMPI
    INTEGER :: local_comm
      
    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_SEND (buffer, count, p_decomposed_type, dest, tag, &
           local_comm, p_error)
    ELSE
      CALL MPI_SEND (buffer, 1, p_decomposed_type, dest, tag, &
           local_comm, p_error)
    END IF
      
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
           ' to ', dest, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)
    END IF
#endif
#endif   

  END SUBROUTINE p_send_decomposed_1d

  !---------------------------------------------------------------------------

  SUBROUTINE p_send_stream (buffer, dest, tag, count, comm)
      
    TYPE(t_stream),      INTENT(in) :: buffer
    INTEGER            , INTENT(in) :: dest, tag
    INTEGER, OPTIONAL  , INTENT(in) :: count, comm
    
#ifndef NOMPI
    INTEGER :: local_comm
      
    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_SEND (buffer, count, p_stream_type, dest, tag, &
           local_comm, p_error)
    ELSE
      CALL MPI_SEND (buffer, 1, p_stream_type, dest, tag, &
           local_comm, p_error)
    END IF
      
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
           ' to ', dest, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)
    END IF
#endif
#endif   

  END SUBROUTINE p_send_stream

  !---------------------------------------------------------------------------

  SUBROUTINE p_send_stream_1d (buffer, dest, tag, count, comm)
      
    TYPE(t_stream),      INTENT(in) :: buffer(:)
    INTEGER            , INTENT(in) :: dest, tag
    INTEGER, OPTIONAL  , INTENT(in) :: count, comm
    
#ifndef NOMPI
    INTEGER :: local_comm

    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_SEND (buffer, count, p_stream_type, dest, tag, &
           local_comm, p_error)
    ELSE
      CALL MPI_SEND (buffer, SIZE(buffer), p_stream_type, dest, tag, &
           local_comm, p_error)
    END IF
      
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
           ' to ', dest, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)
    END IF
#endif
#endif   

  END SUBROUTINE p_send_stream_1d

  !---------------------------------------------------------------------------

  SUBROUTINE p_recv_memory_info (buffer, source, tag, count, comm)

    TYPE(memory_info), INTENT(out) :: buffer
    INTEGER          , INTENT(in)  :: source, tag
    INTEGER, OPTIONAL, INTENT(in)  :: count, comm 
    
#ifndef NOMPI
    INTEGER :: status(MPI_STATUS_SIZE)

    INTEGER :: local_comm

    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_RECV (buffer, count, p_memory_info_type, source, tag, &
           local_comm, status, p_error)
    ELSE

      CALL MPI_RECV (buffer, 1, p_memory_info_type, source, tag, &
           local_comm, status, p_error)
    END IF
      
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
           ' from ', source, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)      
    END IF
#endif
#endif
    
  END SUBROUTINE p_recv_memory_info

  !---------------------------------------------------------------------------

  SUBROUTINE p_recv_memory_info_1d (buffer, source, tag, count, comm)
      
    TYPE(memory_info),    INTENT(out) :: buffer(:)
    INTEGER,          INTENT(in)  :: source, tag
    INTEGER, OPTIONAL, INTENT(in) :: count, comm
    
#ifndef NOMPI
    INTEGER :: status(MPI_STATUS_SIZE)

    INTEGER :: local_comm
    
    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
    
    IF (PRESENT(count)) THEN
      CALL MPI_RECV (buffer, count, p_memory_info_type, source, tag, &
           local_comm, status, p_error)
    ELSE
      CALL MPI_RECV (buffer, SIZE(buffer), p_memory_info_type, source, tag, &
           local_comm, status, p_error)
    END IF
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
           ' from ', source, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)      
    END IF
#endif
#endif
    
  END SUBROUTINE p_recv_memory_info_1d

  !---------------------------------------------------------------------------

  SUBROUTINE p_recv_decomposed (buffer, source, tag, count, comm)

    TYPE(pe_decomposed), INTENT(out) :: buffer
    INTEGER            , INTENT(in)  :: source, tag
    INTEGER, OPTIONAL  , INTENT(in)  :: count, comm 
    
#ifndef NOMPI
    INTEGER :: status(MPI_STATUS_SIZE)

    INTEGER :: local_comm

    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_RECV (buffer, count, p_decomposed_type, source, tag, &
           local_comm, status, p_error)
    ELSE

      CALL MPI_RECV (buffer, 1, p_decomposed_type, source, tag, &
           local_comm, status, p_error)
    END IF
      
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
           ' from ', source, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)      
    END IF
#endif
#endif
    
  END SUBROUTINE p_recv_decomposed

  !---------------------------------------------------------------------------

  SUBROUTINE p_recv_decomposed_1d (buffer, source, tag, count, comm)

    TYPE(pe_decomposed), INTENT(out) :: buffer(:)
    INTEGER            , INTENT(in)  :: source, tag
    INTEGER, OPTIONAL  , INTENT(in)  :: count, comm 
    
#ifndef NOMPI
    INTEGER :: status(MPI_STATUS_SIZE)

    INTEGER :: local_comm

    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_RECV (buffer, count, p_decomposed_type, source, tag, &
           local_comm, status, p_error)
    ELSE

      CALL MPI_RECV (buffer, SIZE(buffer), p_decomposed_type, source, tag, &
           local_comm, status, p_error)
    END IF
      
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
           ' from ', source, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)      
    END IF
#endif
#endif
    
  END SUBROUTINE p_recv_decomposed_1d

  !---------------------------------------------------------------------------

  SUBROUTINE p_recv_stream (buffer, source, tag, count, comm)

    TYPE(t_stream),      INTENT(out) :: buffer
    INTEGER            , INTENT(in)  :: source, tag
    INTEGER, OPTIONAL  , INTENT(in)  :: count, comm 
    
#ifndef NOMPI
    INTEGER :: status(MPI_STATUS_SIZE)

    INTEGER :: local_comm

    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_RECV (buffer, count, p_stream_type, source, tag, &
           local_comm, status, p_error)
    ELSE

      CALL MPI_RECV (buffer, 1, p_stream_type, source, tag, &
           local_comm, status, p_error)
    END IF
      
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
           ' from ', source, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)      
    END IF
#endif
#endif
    
  END SUBROUTINE p_recv_stream

  !---------------------------------------------------------------------------

  SUBROUTINE p_recv_stream_1d (buffer, source, tag, count, comm)

    TYPE(t_stream),      INTENT(out) :: buffer(:)
    INTEGER            , INTENT(in)  :: source, tag
    INTEGER, OPTIONAL  , INTENT(in)  :: count, comm 
    
#ifndef NOMPI
    INTEGER :: status(MPI_STATUS_SIZE)

    INTEGER :: local_comm

    IF (PRESENT(comm)) THEN
      local_comm = comm
    ELSE
      local_comm = p_global_comm
    ENDIF
      
    IF (PRESENT(count)) THEN
      CALL MPI_RECV (buffer, count, p_stream_type, source, tag, &
           local_comm, status, p_error)
    ELSE

      CALL MPI_RECV (buffer, SIZE(buffer), p_stream_type, source, tag, &
           local_comm, status, p_error)
    END IF
      
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (message_text,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
           ' from ', source, ' for tag ', tag, ' failed.'
      CALL message('',message_text)
      WRITE (message_text,'(a,i4)') ' p_error = ', p_error
      CALL finish('',message_text)      
    END IF
#endif
#endif
    
  END SUBROUTINE p_recv_stream_1d

  !---------------------------------------------------------------------------

END MODULE mo_mpi_types
