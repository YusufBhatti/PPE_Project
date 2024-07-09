!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_g2b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add_stream_element

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g2b ! construct the g2b table
  PUBLIC :: destruct_g2b  ! destruct  the g2b table

  PUBLIC :: g2b           ! the g2b table

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: uf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: vf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: dudlf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: dvdlf(:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: g2b

CONTAINS

  SUBROUTINE construct_g2b (lnlon, lnlev, lngl, nlon, nlev, ngl)

    INTEGER, INTENT (in) :: lnlon, lnlev, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ngl

    INTEGER :: dim1(3), dim1p(3)

    ! construct the g2b table
    !
    ! all information specific to this table is set in this subroutine

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    dim1p = (/ lnlon, lnlev, lngl /)
    dim1  = (/  nlon,  nlev,  ngl /)

    CALL add_stream_element (g2b, 'uf', uf, dim1p, dim1)
    CALL add_stream_element (g2b, 'vf', vf, dim1p, dim1)
    CALL add_stream_element (g2b, 'dudlf', dudlf, dim1p, dim1)
    CALL add_stream_element (g2b, 'dvdlf', dvdlf, dim1p, dim1)

  END SUBROUTINE construct_g2b

  SUBROUTINE destruct_g2b

    CALL delete_stream (g2b)

  END SUBROUTINE destruct_g2b

END MODULE mo_memory_g2b
