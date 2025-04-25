!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_f
  !
  ! declaration of predefined fields within this module 
  !
  ! used in inverse Legendre transformation
  !

  ! Modules used

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: add_stream_element, delete_stream, &
                            default_stream_setting, FOURIER
  USE mo_netCDF,      ONLY: max_dim_name
  IMPLICIT NONE

  ! public entities

  PRIVATE
  PUBLIC :: construct_f ! construct the f table
  PUBLIC :: destruct_f  ! destruct  the f table
  PUBLIC :: f           ! the f table

  ! pointers for *f1* space.

  REAL(dp), POINTER, PUBLIC  :: fsvo(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: favo(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fsu(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fau(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fsv(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fav(:,:,:,:)

  ! pointers for *f3* space.

  REAL(dp), POINTER, PUBLIC  :: fsd(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fad(:,:,:,:)

  ! pointers for *f4* space.

  REAL(dp), POINTER, PUBLIC  :: fstp(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fatp(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fstpm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fatpm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fsu0(:,:)
  REAL(dp), POINTER, PUBLIC  :: fau0(:,:)
  REAL(dp), POINTER, PUBLIC  :: fsdu0(:,:)
  REAL(dp), POINTER, PUBLIC  :: fadu0(:,:)

  ! used in direct Legendre transformation
  !
  ! pointers for *f5* space.

  REAL(dp), POINTER, PUBLIC :: fszl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fazl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fszm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fazm(:,:,:,:)

  ! pointers for *f7* space.

  REAL(dp), POINTER, PUBLIC :: fsdl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fadl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fsdm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fadm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fsr(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: far(:,:,:,:)

  ! pointers for *f8* space.

  REAL(dp), POINTER, PUBLIC :: fstp1(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fatp1(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fsul(:,:)
  REAL(dp), POINTER, PUBLIC :: faul(:,:)

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: f

CONTAINS

  SUBROUTINE construct_f (lnlev, lnlevp1, lnmp1, lnhgl, nlev, nmp1, nhgl)

    INTEGER, INTENT (in) :: lnlev, lnlevp1, lnmp1, lnhgl
    INTEGER, INTENT (in) ::  nlev,           nmp1,  nhgl
    INTEGER :: dim1(4), dim1p(4)
    INTEGER :: dim2(4), dim2p(4)
    INTEGER :: dim3(2), dim3p(2)
    CHARACTER (max_dim_name) :: dim1n(4), dim2n(4), dim3n(2)

    ! construct the f table
    !
    ! all information specific to this table is set in this subroutine
    !
    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields
    !
    ! assign pointers

    dim1p = (/ lnlev,     2,       lnmp1,    lnhgl    /)
    dim1  = (/  nlev,     2,        nmp1,     nhgl    /)
    dim1n = (/ "lev    ","complex","nmp1   ","nhgl   "/)

    dim2p = (/ lnlevp1,   2,       lnmp1,    lnhgl    /)
    dim2  = (/  nlev+1,   2,        nmp1,     nhgl    /)
    dim2n = (/ "ilev   ","complex","nmp1   ","nhgl   "/)

    dim3p = (/ lnlev, lnhgl /)
    dim3  = (/  nlev,  nhgl /)
    dim3n = (/ "lev ","nhgl"/)

    ! Arrays used by inverse transform

    CALL default_stream_setting (f, repr   = FOURIER, &
                                    lpost  = .FALSE., &
                                    lrerun = .TRUE.)

    CALL add_stream_element (f, 'fsvo',  fsvo,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'favo',  favo,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsu',   fsu,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fau',   fau,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsv',   fsv,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fav',   fav,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsd',   fsd,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fad',   fad,   dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fstp',  fstp,  dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fatp',  fatp,  dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fstpm', fstpm, dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fatpm', fatpm, dim2p, dim2, dimnames=dim2n)

    CALL add_stream_element (f, 'fsu0',  fsu0,  dim3p, dim3, dimnames=dim3n)
    CALL add_stream_element (f, 'fau0',  fau0,  dim3p, dim3, dimnames=dim3n)
    CALL add_stream_element (f, 'fsdu0', fsdu0, dim3p, dim3, dimnames=dim3n)
    CALL add_stream_element (f, 'fadu0', fadu0, dim3p, dim3, dimnames=dim3n)

    ! Arrays used by direct transform (not yet in memory buffer)

    CALL default_stream_setting (f, lrerun= .FALSE.)

    CALL add_stream_element (f, 'fsdl',  fsdl, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fadm',  fadm, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsr',   fsr,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fszl',  fszl, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fazm',  fazm, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fstp1', fstp1,dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'fsul',  fsul, dim3p, dim3, dimnames=dim3n)


    CALL add_stream_element (f, 'fadl', fadl, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fsdm', fsdm, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'far',  far,  dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fazl', fazl, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fszm', fszm, dim1p, dim1, dimnames=dim1n)
    CALL add_stream_element (f, 'fatp1',fatp1,dim2p, dim2, dimnames=dim2n)
    CALL add_stream_element (f, 'faul', faul, dim3p, dim3, dimnames=dim3n)

  END SUBROUTINE construct_f

  SUBROUTINE destruct_f

    CALL delete_stream (f)

  END SUBROUTINE destruct_f

END MODULE mo_memory_f
