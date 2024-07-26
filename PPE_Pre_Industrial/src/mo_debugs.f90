!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_debugs
  !-----------------------------------------------------------------------
  ! Module containing variables in order to debug mozech
  !
  ! Authors:
  ! J.S. Rast, MPI, September 2003, original source
  !-----------------------------------------------------------------------

  USE mo_kind,                    ONLY: dp
  USE mo_memory_base,             ONLY: new_stream, add_stream_element, &
                                        default_stream_setting, AUTO, t_stream
  USE mo_time_event,              ONLY: io_time_event
  USE mo_time_control,            ONLY: p_bcast_event
  USE mo_netcdf,                  ONLY: HYBRID_H

  IMPLICIT NONE

  PRIVATE
  ! debugs stream
  PUBLIC :: init_debugs
  PUBLIC :: ddf01, ddf02, ddf03, ddf04, ddf05, ddf06, &
            ddf07, ddf08, ddf09, ddf10, ddf11, ddf12, &
            ddf13, ddf14
  PUBLIC :: ddfh01, ddfh02, ddfh03, ddfh04, ddfh05, ddfh06
  PUBLIC :: zdf01, zdf02, zdf03, zdf04, zdf05, zdf06, &
            zdf07, zdf08, zdf09, zdf10, zdf11, zdf12, &
            zdf13, zdf14
  PUBLIC :: nddf, nddfh, nzdf, pvddf, pvddfh, pvzdf

  TYPE t_pvddf
    REAL(dp), POINTER :: v(:,:,:)
  END TYPE t_pvddf
  TYPE t_pvzdf
    REAL(dp), POINTER :: v(:,:)
  END TYPE t_pvzdf
    
  INTEGER           :: nddf=0, nddfh=0, nzdf=0
  REAL(dp), POINTER :: ddf01(:,:,:),ddf02(:,:,:),ddf03(:,:,:),ddf04(:,:,:), &
                       ddf05(:,:,:),ddf06(:,:,:),ddf07(:,:,:),ddf08(:,:,:), &
                       ddf09(:,:,:),ddf10(:,:,:),ddf11(:,:,:),ddf12(:,:,:), &
                       ddf13(:,:,:),ddf14(:,:,:)
  REAL(dp), POINTER :: ddfh01(:,:,:),ddfh02(:,:,:),ddfh03(:,:,:), &
                       ddfh04(:,:,:),ddfh05(:,:,:),ddfh06(:,:,:)
  REAL(dp), POINTER :: zdf01(:,:),  zdf02(:,:),  zdf03(:,:),  zdf04(:,:), &
                       zdf05(:,:),  zdf06(:,:),  zdf07(:,:),  zdf08(:,:), &
                       zdf09(:,:),  zdf10(:,:),  zdf11(:,:),  zdf12(:,:), &
                       zdf13(:,:),  zdf14(:,:)
  TYPE(io_time_event), SAVE     :: putdebug_stream
  TYPE(t_pvddf), POINTER  :: pvddf(:), pvddfh(:)
  TYPE(t_pvzdf), POINTER  :: pvzdf(:)

CONTAINS

  SUBROUTINE init_debugs

    USE mo_mpi,                   ONLY: p_parallel_io, p_barrier, p_bcast, p_io
    USE mo_namelist,              ONLY: open_nml, position_nml, POSITIONED
    USE mo_exception,             ONLY: finish

    CHARACTER(LEN=32) :: ichar
    INTEGER           :: i, j, ierr, inml, iunit, ilength, nlength
    TYPE (t_stream), POINTER      :: debugs

    INCLUDE 'debugsctl.inc'

    ! set default output interval
    putdebug_stream%counter      = 6
    putdebug_stream%unit         = 'hours'
    putdebug_stream%adjustment   = 'first'
    putdebug_stream%offset       = 0
    IF (p_parallel_io) THEN
      inml = open_nml('namelist.echam')
      iunit = position_nml ('DEBUGSCTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
        READ(iunit, debugsctl)
      END SELECT
    END IF
    CALL p_barrier
    CALL p_bcast (nddf, p_io)
    CALL p_bcast (nddfh, p_io)
    CALL p_bcast (nzdf, p_io)
    CALL p_bcast_event (putdebug_stream, p_io)
    CALL new_stream (debugs,'debugs',lrerun=.false.,interval=putdebug_stream)
    CALL default_stream_setting (debugs,lrerun=.false.,contnorest=.true., &
         laccu=.false.,lpost=.true.,table=199, &
         code=AUTO)
    CALL add_stream_element (debugs,'ddf01',ddf01)
    CALL add_stream_element (debugs,'ddf02',ddf02)
    CALL add_stream_element (debugs,'ddf03',ddf03)
    CALL add_stream_element (debugs,'ddf04',ddf04)
    CALL add_stream_element (debugs,'ddf05',ddf05)
    CALL add_stream_element (debugs,'ddf06',ddf06)
    CALL add_stream_element (debugs,'ddf07',ddf07)
    CALL add_stream_element (debugs,'ddf08',ddf08)
    CALL add_stream_element (debugs,'ddf09',ddf09)
    CALL add_stream_element (debugs,'ddf10',ddf10)
    CALL add_stream_element (debugs,'ddf11',ddf11)
    CALL add_stream_element (debugs,'ddf12',ddf12)
    CALL add_stream_element (debugs,'ddf13',ddf13)
    CALL add_stream_element (debugs,'ddf14',ddf14)

    CALL add_stream_element (debugs,'ddfh01',ddfh01,leveltype=HYBRID_H)
    CALL add_stream_element (debugs,'ddfh02',ddfh02,leveltype=HYBRID_H)
    CALL add_stream_element (debugs,'ddfh03',ddfh03,leveltype=HYBRID_H)
    CALL add_stream_element (debugs,'ddfh04',ddfh04,leveltype=HYBRID_H)
    CALL add_stream_element (debugs,'ddfh05',ddfh05,leveltype=HYBRID_H)
    CALL add_stream_element (debugs,'ddfh06',ddfh06,leveltype=HYBRID_H)

    CALL add_stream_element (debugs,'zdf01',zdf01)
    CALL add_stream_element (debugs,'zdf02',zdf02)
    CALL add_stream_element (debugs,'zdf03',zdf03)
    CALL add_stream_element (debugs,'zdf04',zdf04)
    CALL add_stream_element (debugs,'zdf05',zdf05)
    CALL add_stream_element (debugs,'zdf06',zdf06)
    CALL add_stream_element (debugs,'zdf07',zdf07)
    CALL add_stream_element (debugs,'zdf08',zdf08)
    CALL add_stream_element (debugs,'zdf09',zdf09)
    CALL add_stream_element (debugs,'zdf10',zdf10)
    CALL add_stream_element (debugs,'zdf11',zdf11)
    CALL add_stream_element (debugs,'zdf12',zdf12)
    CALL add_stream_element (debugs,'zdf13',zdf13)
    CALL add_stream_element (debugs,'zdf14',zdf14)

! allocate nddf 3-d variables on full levels (layer centres) and store pointers
    ALLOCATE(pvddf(nddf))
    WRITE(ichar,*) nddf
    nlength=LEN_TRIM(ADJUSTL(ichar))
    IF ( nlength > 27 ) THEN
       CALL finish('mo_debugs, init_debugs', 'too many 3d debugs variables')
    END IF
    DO i=1,nddf
      WRITE(ichar,*) i
      ilength=LEN_TRIM(ADJUSTL(ichar))
      DO j=ilength+1,nlength
         ichar='0'//TRIM(ADJUSTL(ichar))
      END DO
      CALL add_stream_element (debugs,'vddf'//TRIM(ADJUSTL(ichar)),pvddf(i)%v)
    END DO

! allocate nddfh 3-d variables on half levels (layer interfaces) and store pointers
    ALLOCATE(pvddfh(nddfh))
    WRITE(ichar,*) nddfh
    nlength=LEN_TRIM(ADJUSTL(ichar))
    IF ( nlength > 26 ) THEN
       CALL finish('mo_debugs, init_debugs', 'too many 3d debugs variables on half levels')
    END IF
    DO i=1,nddfh
      WRITE(ichar,*) i
      ilength=LEN_TRIM(ADJUSTL(ichar))
      DO j=ilength+1,nlength
         ichar='0'//TRIM(ADJUSTL(ichar))
      END DO
      CALL add_stream_element (debugs,'vddfh'//TRIM(ADJUSTL(ichar)),pvddfh(i)%v, &
                               leveltype=HYBRID_H)
    END DO

! allocate nzdf 2-d variables
    ALLOCATE(pvzdf(nzdf))
    WRITE(ichar,*) nzdf
    nlength=LEN_TRIM(ADJUSTL(ichar))
    IF ( nlength > 27 ) THEN
       CALL finish('mo_debugs, init_debugs', 'too many 2d debugs variables')
    END IF
    DO i=1,nzdf
      WRITE(ichar,*) i
      ilength=LEN_TRIM(ADJUSTL(ichar))
      DO j=ilength+1,nlength
         ichar='0'//TRIM(ADJUSTL(ichar))
      END DO
      CALL add_stream_element (debugs,'vzdf'//TRIM(ADJUSTL(ichar)),pvzdf(i)%v)
    END DO
  END SUBROUTINE init_debugs
END MODULE mo_debugs
