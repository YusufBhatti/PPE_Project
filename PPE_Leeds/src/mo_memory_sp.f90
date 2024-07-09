!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_sp

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream ,default_stream_setting  &
                           ,add_stream_element ,get_stream_element &
                           ,SPECTRAL ,HYBRID ,SURFACE ,UNKNOWN

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_sp      ! construct the sp table
  PUBLIC :: destruct_sp       ! destruct  the sp table

  PUBLIC :: sp                ! the sp table

  PUBLIC :: svo, sd, stp, su0 ! predefined fields

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER  :: svo(:,:,:) ! vorticity
  REAL(dp), POINTER  :: sd (:,:,:) ! divergence
  REAL(dp), POINTER  :: stp(:,:,:) ! surface pressure, temperature
  REAL(dp), POINTER  :: su0(:,:)   ! 

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: sp

  ! 4D pointer stp4 to be shared by stp, sp, t
  ! 3D pointer for dummy entries sp, t used for output

  REAL(dp), POINTER  :: stp4(:,:,:,:)
  REAL(dp), POINTER  :: st4 (:,:,:,:)
  REAL(dp), POINTER  :: sp4 (:,:,:,:)
  REAL(dp), POINTER  :: st3 (:,:,:)
  REAL(dp), POINTER  :: sp3 (:,:,:)

CONTAINS

  SUBROUTINE construct_sp (lnlev, lnnp1, lnsp, nlev, nnp1, nsp)

    INTEGER, INTENT (in) :: lnlev, lnnp1, lnsp
    INTEGER, INTENT (in) :: nlev, nnp1, nsp

    ! construct the sp table
    !
    ! all information specific to this table is set in this subroutine

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    CALL default_stream_setting (sp &
                                ,table     = 128      &
                                ,bits      = 16       &
                                ,repr      = SPECTRAL &
                                ,leveltype = HYBRID   &
                                ,lpost     = .TRUE.)

    CALL add_stream_element (sp ,'svo' ,svo ,(/lnlev,  2,lnsp/) ,(/nlev,  2,nsp/) ,code=138 ,longname='vorticity' ,units='1/s')
    CALL add_stream_element (sp ,'sd'  ,sd  ,(/lnlev,  2,lnsp/) ,(/nlev,  2,nsp/) ,code=155 ,longname='divergence',units='1/s')
    CALL add_stream_element (sp ,'stp' ,stp ,(/lnlev+1,2,lnsp/) ,(/nlev+1,2,nsp/) ,lpost=.FALSE. )
    CALL add_stream_element (sp ,'su0' ,su0 ,(/lnlev,lnnp1/)    ,(/nlev,nnp1/)    ,lpost=.FALSE. &
                                            ,leveltype=UNKNOWN)

    ! for output create dummy entries sp, t (surface pressure, temperature)

    CALL get_stream_element (sp,'stp',stp4)
    st4 => stp4 (1:lnlev        ,:,:,:)
    sp4 => stp4 (lnlev+1:lnlev+1,:,:,:)

    CALL add_stream_element (sp ,'st' ,st3 ,(/lnlev, 2, lnsp/) ,(/nlev, 2, nsp/) ,p4=st4 &
                            ,code=130 ,longname='temperature',units='K')
    CALL add_stream_element (sp ,'lsp',sp3 ,(/1,     2, lnsp/) ,(/1,    2, nsp/) ,p4=sp4 &
                            ,code=152 ,longname='log surface pressure' ,leveltype=SURFACE)

  END SUBROUTINE construct_sp

  SUBROUTINE destruct_sp

    CALL delete_stream (sp)

  END SUBROUTINE destruct_sp

END MODULE mo_memory_sp
