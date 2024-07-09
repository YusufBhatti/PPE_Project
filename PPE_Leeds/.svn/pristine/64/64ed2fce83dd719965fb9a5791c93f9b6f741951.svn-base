!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!!
!! These routines belong to the standalone Version of JSBACH
!!
  SUBROUTINE jsbalone_init_decomp(nlon, nlat, mask, nproca, nprocb, npedim)

    USE mo_decomposition, ONLY : decompose,local_decomposition, global_decomposition
    USE mo_mpi,           ONLY : p_nprocs,p_parallel_io, p_parallel, p_io, p_pe
    USE mo_jsbach,        ONLY : debug
    USE mo_transpose,     ONLY : indx
    USE mo_exception,     ONLY : message, finish
    USE mo_io_units,      ONLY : nout
    USE mo_decomposition, ONLY : print_decomposition

    INTEGER, INTENT(in) :: nlon, nlat
    INTEGER, INTENT(in) :: nproca, nprocb, npedim
    LOGICAL, INTENT(in) :: mask(nlon,nlat)

    INTEGER :: p, i
    INTEGER :: debug_parallel = -1   ! Parallel debug flag: -1= no debugging, 0,1=debugging

    IF (debug) CALL message('jsbalone_init_decomp','Entering ...')

#ifdef NOMPI
    IF (nproca*nprocb > 1) THEN
       CALL message('', 'Parallel execution selected with wrong compiler options')
       CALL message('', 'Please recompile without -DNOMPI')
       CALL finish('jsbalone_init_decomp', 'Program aborted')
    ENDIF
#endif

    ! nprocb has to be 1, otherwise the consitency checks in
    ! mo_decomposition will fail (nprocb - must be < nlev).
    IF (nprocb > 1) THEN
       CALL finish('jsbalone_init_decomp', 'nprocb has to be 1')
    ENDIF

    IF (p_nprocs == nproca*nprocb) THEN
       debug_parallel = -1
    ELSE IF (p_nprocs == nproca*nprocb+1) THEN
       debug_parallel = 0
    ELSE
       CALL finish('jsbalone_init_decomp', 'Number of runtime PEs doesn''t fit nproca*nprocb(+1)')
    ENDIF

    IF (p_parallel .AND. p_parallel_io) THEN
       WRITE(nout, '(/,a,i4,a,i3,a,i3)') ' Total number of PEs: ', p_nprocs, &
                                         ' Set A: ', nproca, ' Set B: ', nprocb
    ENDIF

    ! Derive decomposition
    ALLOCATE(global_decomposition(1:p_nprocs))

    ! The truncation is set to nproca, otherwise the consitency checks in
    ! mo_decomposition will fail 
    CALL decompose(global_decomposition, npedim, nproca, nprocb, nlat, nlon, nlev=1, &
                   nm=nproca, nn=nproca, nk=nproca, nadvec=0, mask=mask, &
                   debug=debug_parallel)

    ! Keep index values, not id values
    DO p=1,p_nprocs
       global_decomposition(p)%mapmesh(:,:) = &
            indx(global_decomposition(p)%mapmesh(:,:),global_decomposition)
    ENDDO

    ! Copy global decomposition table entry to local decomposition
    DO i=1,p_nprocs
       IF (global_decomposition(i)% pe == p_pe) THEN
          local_decomposition = global_decomposition(i)
       END IF
    END DO

    IF (p_parallel_io) THEN
       WRITE(nout,*) '---------------------------------------------'
       WRITE(nout,*) ' Blocking information:'
       WRITE(nout,*) '   lreg  nproma ngpblks  npromz'
       WRITE(nout,'(L5,3I8)') local_decomposition% lreg, &
            local_decomposition% nproma, &
            local_decomposition% ngpblks, &
            local_decomposition% npromz
       WRITE(nout,*) '---------------------------------------------'
    ENDIF
    
   
    ! decomposition printout
    IF (debug .AND. p_pe == p_io .AND. p_parallel) THEN
       DO i = 1, p_nprocs
          WRITE (nout,'(78("-"))')
          CALL print_decomposition(global_decomposition(i))
       END DO
       WRITE (nout,'(78("-"),/)')
    END IF

    IF (debug) CALL message('jsbalone_init_decomp','... returning')

  END SUBROUTINE jsbalone_init_decomp

