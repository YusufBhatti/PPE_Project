!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_g2a

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add_stream_element
  USE mo_netCDF,      ONLY: max_dim_name

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g2a ! construct the g2a table
  PUBLIC :: destruct_g2a  ! destruct  the g2a table

  PUBLIC :: g2a           ! the g2a table

  PUBLIC :: add_stream_element

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: um1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: vm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: dudlm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: dvdlm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: dtlm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: dtmm1(:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: g2a

CONTAINS

  SUBROUTINE construct_g2a (lnlon, lnlev, lngl, nlon, nlev, ngl)

    INTEGER, INTENT (in) :: lnlon, lnlev, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ngl

    INTEGER :: dim1(3), dim1p(3)
    CHARACTER (max_dim_name) :: dim1n(3)

    ! construct the g2a table
    !
    ! all information specific to this table is set in this subroutine


    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    dim1p = (/ lnlon,  lnlev, lngl  /)
    dim1  = (/  nlon,   nlev,  ngl  /)
    dim1n = (/ "lon",  "lev", "lat" /)

    CALL add_stream_element (g2a, 'um1' , um1, dim1p, dim1, dimnames=dim1n, &
                                               units='m/s', lrerun=.TRUE.)
    CALL add_stream_element (g2a, 'vm1' , vm1, dim1p, dim1, dimnames=dim1n, &
                                               units='m/s', lrerun=.TRUE.)
    CALL add_stream_element (g2a, 'dudlm1' , dudlm1, dim1p, dim1, &
                                            dimnames=dim1n, lrerun=.TRUE.)
    CALL add_stream_element (g2a, 'dvdlm1' , dvdlm1, dim1p, dim1, &
                                            dimnames=dim1n, lrerun=.TRUE.)
    CALL add_stream_element (g2a, 'dtlm1',dtlm1,dim1p,dim1, dimnames=dim1n, &
                                                            lrerun=.TRUE.)
    CALL add_stream_element (g2a, 'dtmm1',dtmm1,dim1p,dim1, dimnames=dim1n, &
                                                            lrerun=.TRUE.)

  END SUBROUTINE construct_g2a

  SUBROUTINE destruct_g2a

    CALL delete_stream (g2a)

  END SUBROUTINE destruct_g2a

END MODULE mo_memory_g2a
