!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!> I/O servers do not participate in the coupling. This routine combines calls
!> to collective OASIS3-MCT interface routines and is called by all I/O 
!> servers after communicator split
SUBROUTINE pio_uncouple
    USE mo_mpi,         ONLY: p_comm_null
    USE mo_io_units,    ONLY: nerr
#ifdef __prism
#ifdef __oa3mct
    USE mo_couple_wrap, ONLY: oasis_set_couplcomm, oasis_def_partition, &
                              oasis_def_var, oasis_enddef, oasis_abort, &
                              oasis_ok, oasis_in, oasis_out, oasis_real, &
                              namdstfld, namsrcfld

    IMPLICIT NONE 

    INTEGER :: ierror, retid, partid, jfld
    INTEGER, PARAMETER :: ig_paral(3) = (/0, 0, 0/), &
                          var_nodims(2) = (/2, 1/), &
                          var_actual_shape(4) = (/1,1,1,1/)

    ! set coupling communicator to p_comm_null for I/O servers since they are
    ! not involved in the coupling
    CALL oasis_set_couplcomm(p_comm_null, ierror)

    ! oasis_def_partition is collective and must be called by all processes
    CALL oasis_def_partition (partid, ig_paral, ierror)

    ! oasis_def_var is collective and must be called by all processes
    !   define ports of incoming fields
    DO jfld = 1,size(namdstfld)
       CALL oasis_def_var (retid, TRIM(namdstfld(jfld)), partid, var_nodims, &
                           oasis_in, var_actual_shape, oasis_real, ierror)
       IF (ierror /= oasis_ok) &
          CALL oasis_abort(0, 'pio_uncouple', &
              'oasis_def_var failed for model echam6 (I/O servers)')
    ENDDO

    ! oasis_enddef is collective and must be called by all processes
    CALL oasis_enddef(ierror)
    IF (ierror /= oasis_ok) &
       CALL oasis_abort(0, 'pio_uncouple', &
           'oasis_enddef failed for model echam6 (I/O servers)')
#else
    CALL oasis_abort(0, 'pio_uncouple', &
        'coupler oasis3 unhandled yet, use nprocio=0')
#endif
#endif

END SUBROUTINE pio_uncouple



