!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_g1a

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add_stream_element, &
                            default_stream_setting
  USE mo_netCDF,      ONLY: max_dim_name
  USE mo_tracdef,     ONLY: trlist
  USE mo_memory_gl,   ONLY: tracer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g1a ! construct the g1a table
  PUBLIC :: destruct_g1a  ! destruct  the g1a table

  PUBLIC :: g1a           ! the g1a table

  PUBLIC :: add_stream_element

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: vom1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: dm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: alpsm1(:,:)
  REAL(dp), POINTER, PUBLIC :: dalpslm1(:,:)
  REAL(dp), POINTER, PUBLIC :: dalpsmm1(:,:)
  REAL(dp), POINTER, PUBLIC :: qm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xlm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xim1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xtm1(:,:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: g1a

  ! private storage for tracer fields

  REAL(dp), POINTER :: pxtm1 (:,:,:,:,:)

CONTAINS

  SUBROUTINE construct_g1a (lnlon, lnlev, lntrac, lngl, &
                             nlon,  nlev,  ntrac,  ngl)

    INTEGER, INTENT (in) :: lnlon, lnlev, lntrac, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ntrac,  ngl

    INTEGER                      :: dim1(3), dim1p(3)
    INTEGER                      :: dim2(2), dim2p(2)
    INTEGER                      :: dim3(4), dim3p(4)
    CHARACTER (max_dim_name)     :: dim1n(3), dim2n(2), dim3n(4)
    INTEGER                      :: i
    REAL(dp), POINTER            :: p3(:,:,:), p4(:,:,:,:)

    ! construct the g1a table
    !
    ! all information specific to this table is set in this subroutine

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    dim1p = (/ lnlon,  lnlev, lngl  /)
    dim1  = (/  nlon,   nlev,  ngl  /)
    dim1n = (/  "lon ","lev ","lat "/)

    dim2p = (/ lnlon, lngl  /)
    dim2  = (/  nlon,  ngl  /)
    dim2n = (/  "lon","lat" /)

    dim3p = (/ lnlon,    lnlev,   lntrac,  lngl    /)
    dim3  = (/  nlon,     nlev,    ntrac,   ngl    /)
    dim3n = (/  "lon   ","lev   ","ntrac ","lat   "/)

    CALL default_stream_setting (g1a, lrerun=.TRUE.)

    CALL add_stream_element (g1a, 'vom1',    vom1,    dim1p,dim1,dimnames=dim1n)
    CALL add_stream_element (g1a, 'dm1',     dm1,     dim1p,dim1,dimnames=dim1n)
    CALL add_stream_element (g1a, 'tm1',     tm1,     dim1p,dim1,dimnames=dim1n, units='K')
    CALL add_stream_element (g1a, 'alpsm1',  alpsm1,  dim2p,dim2,dimnames=dim2n)
    CALL add_stream_element (g1a, 'dalpslm1',dalpslm1,dim2p,dim2,dimnames=dim2n)
    CALL add_stream_element (g1a, 'dalpsmm1',dalpsmm1,dim2p,dim2,dimnames=dim2n)
    CALL add_stream_element (g1a, 'qm1',     qm1,     dim1p,dim1,dimnames=dim1n)
    CALL add_stream_element (g1a, 'xlm1',    xlm1,    dim1p,dim1,dimnames=dim1n)
    CALL add_stream_element (g1a, 'xim1',    xim1,    dim1p,dim1,dimnames=dim1n)
    !
    ! Special handling for tracers
    !

    IF (ntrac > 0) THEN
      !
      ! Allocate a 5d-array (with dummy index 5) for tracers. This array
      ! is referenced by the 4-d array 'XTM1' and by the 3-d arrays of
      ! individual tracers.
      !
      ALLOCATE (pxtm1(lnlon, lnlev, lntrac, lngl, 1))
      !
      ! Set meta information on XTM1 array.
      ! Set restart flag, obtain reference to memory info entry.
      !
      p4 => pxtm1(:,:,:,:,1)
      CALL add_stream_element (g1a, 'xtm1', xtm1, dim3p, dim3, &
             dimnames   = dim3n,                               &
             lrerun     = trlist% oldrestart,                  &
             contnorest = .TRUE.,                              &
             mem_info   = trlist% mixtm1,                      &
             p4         = p4)
      !
      ! provide additional meta-information for individual tracers.
      !
      DO i = 1, lntrac
        p4 => pxtm1(:,:,i,:,:)
        CALL add_stream_element (tracer, TRIM(trlist% ti(i)% fullname)//'_m1',&
                                     p3, dim1p, dim1,                         &
                        dimnames   = dim1n,                                   &
                        units      = trlist% ti(i)% units,                    &
                        lrerun     = trlist% ti(i)% nrerun==1,                &
                        contnorest = .TRUE.,                                  &
                        mem_info   = trlist% mi(i)% xtm1,                     &
                        lpost      = .FALSE.,                                 &
                        p4         = p4)
      END DO
    ELSE
      NULLIFY (pxtm1)
      CALL add_stream_element (g1a, 'xtm1', xtm1, dim3p, dim3, &
                        mem_info = trlist% mixtm1,             &
                        lrerun   = .FALSE.)
    END IF


  END SUBROUTINE construct_g1a

  SUBROUTINE destruct_g1a

    CALL delete_stream (g1a)

    IF(ASSOCIATED (pxtm1)) DEALLOCATE (pxtm1)

  END SUBROUTINE destruct_g1a

END MODULE mo_memory_g1a
