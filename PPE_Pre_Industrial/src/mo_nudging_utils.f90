!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_nudging_utils
!BOP
  ! !MODULE: mo_nudging_utils (layer 3)

  ! !DESCRIPTION: 
  ! contains utility functions used for nudging

  ! !REVISION HISTORY: 
  ! Ingo Kirchner, MPI Hamburg, April-2001
  ! R. Johanni, IPP Garching, May-2002, parallel version
  ! Ingo Kirchner, MPI Hamburg, Aug 2002, revision

!BOX
  IMPLICIT NONE

  PRIVATE
!EOX
  ! !PUBLIC MEMBER FUNCTIONS:

  PUBLIC :: Nudg_Correl
  PUBLIC :: cpbread

#ifdef LITTLE_ENDIAN
  PUBLIC :: swap64
  INTERFACE swap64
     MODULE PROCEDURE r_swap64
     MODULE PROCEDURE i_swap64
     MODULE PROCEDURE c_swap64
  END INTERFACE
#endif
!EOP

  INTEGER, PARAMETER, PUBLIC :: HEAD_LEN = 8   ! number of header words
  INTEGER, PARAMETER, PUBLIC :: WORD_LEN = 8   ! number of bytes in words, used for Cray binary files


CONTAINS

  !======================================================================
!BOP
  ! !IROUTINE:  Nudg_Correl
  ! !INTERFACE:

  SUBROUTINE Nudg_Correl

    ! !DESCRIPTION: 
    ! correlations between the nudging term and the tendency diagnostic terms

    ! !USES:
    USE mo_kind,           ONLY: dp
    USE mo_nudging_buffer, ONLY: sdtac, svtac, sttac, sstac, nio_index
    USE mo_nudging_constants, ONLY: ndunit
    USE mo_control,        ONLY: nlev
    USE mo_diag_tendency,  ONLY: &
         ndvor, nddiv, ndtem, ndprs, &! number of diagnostic terms
         pdvor, pddiv, pdtem, pdprs, &! accumulated terms of equations
         dio_index
    USE mo_mpi,            ONLY: p_parallel
    USE mo_exception,      ONLY: finish, message
    USE mo_time_conversion,ONLY: print_date, TC_PRN_NATIVE
    USE mo_time_control,   ONLY: next_date, get_interval_seconds, ev_putdata
    USE mo_spectral,       ONLY: corrsp
    ! !INPUT PARAMETERS: 

!EOP

    INTEGER       :: jk, jl, isec1, isec2
    REAL(kind=dp) :: cc_tem(ndtem), cc_vor(ndvor), cc_div(nddiv), cc_prs(ndprs)
    CHARACTER(len=256) :: mess
    REAL(dp), POINTER :: p3(:,:,:), p2a(:,:), p2b(:,:)

    IF (p_parallel) CALL finish ('Nudg_Correl', &
         'ltdiag must NOT be set in parallel run!!!!!!')

    isec1 = get_interval_seconds(ev_putdata(nio_index))
    isec2 = get_interval_seconds(ev_putdata(dio_index))

    CALL print_date(next_date,TC_PRN_NATIVE,mess=mess)
    WRITE(ndunit,'(a,a)') '###  ',TRIM(mess)

    IF (isec1 /= isec2) THEN
      CALL message('Nudg_Correl',&
           'accumulation interval mismatch (Nudging, TDiag) - NO correlation performed!')
      WRITE(ndunit,'(a)') '### accumulation interval mismatch -- no correlation performed'
      RETURN
    END IF

    ! divergence equation
    WRITE(ndunit,'(a,20(1x,a7))') 'DIV_    ',&
         'dynAPCG', 'VerAdv ', 'VerDiff', 'GWDrag ', &
         'CuCall ', 'TimFilt', 'SemiImp', 'HorDiff', &
         'SumTend'
    
    DO jk=1,NDDIV
      p3 => pddiv(:,:,:,jk)
      cc_div(jk) = corrsp(sdtac,p3)
    END DO
    WRITE(ndunit,'(a,20(1x,f7.3))') 'DIV_ALL ',cc_div(:)
    
    DO jl=1,nlev
      DO jk=1,NDDIV
        p2a => sdtac(jl,:,:)
        p2b => pddiv(jl,:,:,jk)
        cc_div(jk) = corrsp(p2a,p2b)
      END DO
      WRITE(ndunit,'(a,i3,20(1x,f7.3))') 'DIV_L',jl,cc_div(:)
    END DO
  
    ! vorticity equation
    WRITE(ndunit,'(a,20(1x,a7))') 'VOR_    ',&
         'AdPrCor', 'VerAdv ', 'VerDiff', 'GWDrag ', &
         'CuCall ', 'TimFilt', 'SemiImp', 'HorDiff', &
         'SumTend'
    DO jk=1,NDVOR
      p3 => pdvor(:,:,:,jk)
      cc_vor(jk) = corrsp(svtac,p3)
    END DO
    WRITE(ndunit,'(a,20(1x,f7.3))') 'VOR_ALL ',cc_vor(:)
    
    DO jl=1,nlev
      DO jk=1,NDVOR
        p2a => svtac(jl,:,:)
        p2b => pdvor(jl,:,:,jk)
        cc_vor(jk) = corrsp(p2a,p2b)
      END DO
      WRITE(ndunit,'(a,i3,20(1x,f7.3))') 'VOR_L',jl,cc_vor(:)
    END DO

    ! temperature equation
    WRITE(ndunit,'(a,20(1x,a7))') 'TEM_    ',&
         'HorAdv ', 'VerAdv ', 'EnConv ', 'Radheat', &
         'VerDiff', 'GWDrag ', 'CuCall ', 'Convect', &
         'TimFilt', 'SemiImp', 'HorDiff', 'RadLong', &
         'RadSola', 'SumTend'
    DO jk=1,NDTEM
      p3 => pdtem(:,:,:,jk)
      cc_tem(jk) = corrsp(sttac,p3)
    END DO
    WRITE(ndunit,'(a,20(1x,f7.3))') 'TEM_ALL ',cc_tem(:)
    
    DO jl=1,nlev
      DO jk=1,NDTEM
        p2a => sttac(jl,:,:)
        p2b => pdtem(jl,:,:,jk)
        cc_tem(jk) = corrsp(p2a,p2b)
      END DO
      WRITE(ndunit,'(a,i3,20(1x,f7.3))') 'TEM_L',jl,cc_tem(:)
    END DO
    
    ! pressure equation
    WRITE(ndunit,'(a,20(1x,a7))') 'PRS_    ',&
         'ConvInt', 'TimFilt', 'SemiImp', 'SumTend'
    p2a => sstac(1,:,:)
    DO jk=1,NDPRS
      p2b => pdprs(:,:,jk)
      cc_prs(jk) = corrsp(p2a,p2b)
    END DO
    WRITE(ndunit,'(a,i3,20(1x,f7.3))') 'PRS_L',0,cc_prs(:)
    
  END SUBROUTINE Nudg_Correl

  !======================================================================

#ifdef LITTLE_ENDIAN

  SUBROUTINE r_swap64(rfield,idx)
    USE mo_kind, ONLY : dp
    REAL(kind=dp) :: rfield(:)
    INTEGER :: idx
    CALL swap64_main(rfield=rfield,idx=idx)
  END SUBROUTINE r_swap64

  SUBROUTINE i_swap64(ifield,idx)
    INTEGER :: ifield(:)
    INTEGER :: idx
    CALL swap64_main(ifield=ifield,idx=idx)
  END SUBROUTINE i_swap64

  SUBROUTINE c_swap64(cfield,idx)
    CHARACTER(len=8) :: cfield(:)
    INTEGER :: idx
    CALL swap64_main(cfield=cfield,idx=idx)
  END SUBROUTINE c_swap64

!BOP
  ! !IROUTINE:  swap64
  ! ! !INTERFACE:

  SUBROUTINE swap64_main(rfield,ifield,cfield,idx)

    ! !DESCRIPTION: 
    ! swap 8-byte sequences of real, integer or 8-byte words
    ! from big to little endian

    ! !USES:
    USE mo_kind, ONLY : dp

    ! !INPUT/OUTPUT PARAMETERS:
    REAL(kind=dp), OPTIONAL, INTENT(inout)    :: rfield(:)
    INTEGER, OPTIONAL, INTENT(inout)          :: ifield(:)
    CHARACTER(len=WORD_LEN), OPTIONAL, INTENT(inout) :: cfield(:)

    ! !INPUT PARAMETERS: 
    INTEGER                                   :: idx  ! no of elements
!EOP

    INTEGER :: i, j
    CHARACTER(len=WORD_LEN),ALLOCATABLE :: ctr(:)
    CHARACTER(len=1)                    :: tmp(8)

    ALLOCATE(ctr(idx))

    IF (PRESENT(rfield)) THEN
       ctr(1:idx) = TRANSFER(rfield(1:idx),ctr)

    ELSE IF(PRESENT(ifield)) THEN
       ctr(1:idx) = TRANSFER(ifield(1:idx),ctr)

    ELSE IF(PRESENT(cfield)) THEN
       ctr(1:idx) = cfield(1:idx)

    END IF

    ! switch bytes
    DO i=1,idx
       DO j=1,WORD_LEN
          tmp(j) = ctr(i)(j:j)
       END DO
       DO j=1,WORD_LEN
          ctr(i)(j:j) = tmp(9-j)
       END DO
    END DO

    IF (PRESENT(rfield)) THEN
       rfield(1:idx) = TRANSFER(ctr(1:idx),rfield)

    ELSE IF(PRESENT(ifield)) THEN
       ifield(1:idx) = TRANSFER(ctr(1:idx),ifield)

    ELSE IF(PRESENT(cfield)) THEN
       cfield(1:idx) = ctr(1:idx)

    END IF

    DEALLOCATE(ctr)

  END SUBROUTINE swap64_main

#endif

  !======================================================================
!BOP
  ! !IROUTINE:  cpbread
  ! !INTERFACE:

  SUBROUTINE cpbread(unit,cfield,nbytes,ierr)

    ! !DESCRIPTION: 
    ! read and swap service file header

    ! !USES:
    USE mo_kind, ONLY: dp

    ! !INPUT PARAMETERS: 
    INTEGER, INTENT(in)             :: unit       ! file ID
    INTEGER, INTENT(in)             :: nbytes     ! no of bytes to read

    ! !INPUT/OUTPUT PARAMETERS: 
    CHARACTER(len=WORD_LEN), INTENT(inout) :: cfield(:)  ! service header

    ! !OUTPUT PARAMETERS:
    INTEGER, INTENT(out)            :: ierr       ! return code
!EOP

    REAL(kind=dp), ALLOCATABLE :: field(:)

    EXTERNAL pbread

    INTEGER :: no

    no = INT(nbytes/WORD_LEN)
    ALLOCATE(field(no))

    CALL pbread(unit,field(1),nbytes,ierr)
    IF (ierr >= 0) THEN
       cfield = TRANSFER(field(1:no),cfield)
    END IF

#ifdef LITTLE_ENDIAN
    CALL swap64(cfield,no)
#endif

    DEALLOCATE(field)

  END SUBROUTINE cpbread

  !======================================================================

END MODULE mo_nudging_utils
