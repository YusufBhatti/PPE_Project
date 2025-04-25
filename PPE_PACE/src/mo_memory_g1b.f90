!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_g1b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add_stream_element

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g1b ! construct the g1b table
  PUBLIC :: destruct_g1b  ! destruct  the g1b table

  PUBLIC :: g1b           ! the g1b table

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: vof(:,:,:)
  REAL(dp), POINTER, PUBLIC :: df(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: alpsf(:,:)
  REAL(dp), POINTER, PUBLIC :: qf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xlf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xif(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xtf(:,:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: g1b

CONTAINS

  SUBROUTINE construct_g1b (lnlon, lnlev, lntrac, lngl, &
                             nlon,  nlev,  ntrac,  ngl)

    INTEGER, INTENT (in) :: lnlon, lnlev, lntrac, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ntrac,  ngl

    INTEGER :: dim1(3), dim1p(3)
    INTEGER :: dim2(2), dim2p(2)
    INTEGER :: dim3(4), dim3p(4)

    ! construct the g1b table
    !
    ! all information specific to this table is set in this subroutine

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    dim1p = (/ lnlon, lnlev, lngl  /)
    dim1  = (/  nlon,  nlev,  ngl  /)

    dim2p = (/ lnlon, lngl  /)
    dim2  = (/  nlon,  ngl  /)

    dim3p = (/ lnlon, lnlev,  lntrac, lngl /)
    dim3  = (/  nlon,  nlev,   ntrac,  ngl /)

    CALL add_stream_element (g1b, 'vof',   vof,   dim1p, dim1)
    CALL add_stream_element (g1b, 'df',    df,    dim1p, dim1)
    CALL add_stream_element (g1b, 'tf',    tf,    dim1p, dim1)
    CALL add_stream_element (g1b, 'alpsf', alpsf, dim2p, dim2)
    CALL add_stream_element (g1b, 'qf',    qf,    dim1p, dim1)
    CALL add_stream_element (g1b, 'xlf',   xlf,   dim1p, dim1)
    CALL add_stream_element (g1b, 'xif',   xif,   dim1p, dim1)
    CALL add_stream_element (g1b, 'xtf',   xtf,   dim3p, dim3)

  END SUBROUTINE construct_g1b

  SUBROUTINE destruct_g1b

    CALL delete_stream (g1b)

  END SUBROUTINE destruct_g1b

END MODULE mo_memory_g1b
