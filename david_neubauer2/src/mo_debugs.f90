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
            ddf13, ddf14, ddf15, ddf16, ddf17, ddf18, &
            ddf19, ddf20, ddf21
  PUBLIC :: ddfh01, ddfh02, ddfh03, ddfh04, ddfh05, ddfh06
  PUBLIC :: zdf01, zdf02, zdf03, zdf04, zdf05, zdf06, &
            zdf07, zdf08, zdf09, zdf10, zdf11, zdf12, &
            zdf13, zdf14, zdf15, zdf16, zdf17, zdf18, &
            zdf19, zdf20, zdf21, zdf22, zdf23, zdf24, &
            zdf25, zdf26, zdf27, zdf28, zdf29, zdf30, &
            zdf31, zdf32, zdf33, zdf34, zdf35, zdf36, &
            zdf37, zdf38, zdf39, zdf40, zdf41, zdf42, &
            zdf43, zdf44, zdf45, zdf46, zdf47, zdf48, &
            zdf49, zdf50, zdf51, zdf52, zdf53, zdf54, &
            zdf55, zdf56, zdf57, zdf58, zdf59, zdf60, &
            zdf61, zdf62, zdf63, zdf64
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
                       ddf13(:,:,:),ddf14(:,:,:),ddf15(:,:,:),ddf16(:,:,:), &
                       ddf17(:,:,:),ddf18(:,:,:),ddf19(:,:,:),ddf20(:,:,:), &
                       ddf21(:,:,:)
  REAL(dp), POINTER :: ddfh01(:,:,:),ddfh02(:,:,:),ddfh03(:,:,:), &
                       ddfh04(:,:,:),ddfh05(:,:,:),ddfh06(:,:,:)
  REAL(dp), POINTER :: zdf01(:,:),  zdf02(:,:),  zdf03(:,:),  zdf04(:,:), &
                       zdf05(:,:),  zdf06(:,:),  zdf07(:,:),  zdf08(:,:), &
                       zdf09(:,:),  zdf10(:,:),  zdf11(:,:),  zdf12(:,:), &
                       zdf13(:,:),  zdf14(:,:),  zdf15(:,:),  zdf16(:,:), &
                       zdf17(:,:),  zdf18(:,:),  zdf19(:,:),  zdf20(:,:), &
                       zdf21(:,:),  zdf22(:,:),  zdf23(:,:),  zdf24(:,:), &
                       zdf25(:,:),  zdf26(:,:),  zdf27(:,:),  zdf28(:,:), &
                       zdf29(:,:),  zdf30(:,:),  zdf31(:,:),  zdf32(:,:), &
                       zdf33(:,:),  zdf34(:,:),  zdf35(:,:),  zdf36(:,:), &
                       zdf37(:,:),  zdf38(:,:),  zdf39(:,:),  zdf40(:,:), &
                       zdf41(:,:),  zdf42(:,:),  zdf43(:,:),  zdf44(:,:), &
                       zdf45(:,:),  zdf46(:,:),  zdf47(:,:),  zdf48(:,:), &
                       zdf49(:,:),  zdf50(:,:),  zdf51(:,:),  zdf52(:,:), &
                       zdf53(:,:),  zdf54(:,:),  zdf55(:,:),  zdf56(:,:), &
                       zdf57(:,:),  zdf58(:,:),  zdf59(:,:),  zdf60(:,:), &
                       zdf61(:,:),  zdf62(:,:),  zdf63(:,:),  zdf64(:,:)
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
         laccu=.true.,lpost=.true.,table=199, &
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
    CALL add_stream_element (debugs,'ddf15',ddf15)
    CALL add_stream_element (debugs,'ddf16',ddf16)
    CALL add_stream_element (debugs,'ddf17',ddf17)
    CALL add_stream_element (debugs,'ddf18',ddf18)
    CALL add_stream_element (debugs,'ddf19',ddf19)
    CALL add_stream_element (debugs,'ddf20',ddf20)
    CALL add_stream_element (debugs,'ddf21',ddf21)

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
    CALL add_stream_element (debugs,'zdf15',zdf15)
    CALL add_stream_element (debugs,'zdf16',zdf16)
    CALL add_stream_element (debugs,'zdf17',zdf17)
    CALL add_stream_element (debugs,'zdf18',zdf18)
    CALL add_stream_element (debugs,'zdf19',zdf19)
    CALL add_stream_element (debugs,'zdf20',zdf20)
    CALL add_stream_element (debugs,'zdf21',zdf21)
    CALL add_stream_element (debugs,'zdf22',zdf22)
    CALL add_stream_element (debugs,'zdf23',zdf23)
    CALL add_stream_element (debugs,'zdf24',zdf24)
    CALL add_stream_element (debugs,'zdf25',zdf25)
    CALL add_stream_element (debugs,'zdf26',zdf26)
    CALL add_stream_element (debugs,'zdf27',zdf27)
    CALL add_stream_element (debugs,'zdf28',zdf28)
    CALL add_stream_element (debugs,'zdf29',zdf29)
    CALL add_stream_element (debugs,'zdf30',zdf30)
    CALL add_stream_element (debugs,'zdf31',zdf31)
    CALL add_stream_element (debugs,'zdf32',zdf32)
    CALL add_stream_element (debugs,'zdf33',zdf33)
    CALL add_stream_element (debugs,'zdf34',zdf34)
    CALL add_stream_element (debugs,'zdf35',zdf35)
    CALL add_stream_element (debugs,'zdf36',zdf36)
    CALL add_stream_element (debugs,'zdf37',zdf37)
    CALL add_stream_element (debugs,'zdf38',zdf38)
    CALL add_stream_element (debugs,'zdf39',zdf39)
    CALL add_stream_element (debugs,'zdf40',zdf40)
    CALL add_stream_element (debugs,'zdf41',zdf41)
    CALL add_stream_element (debugs,'zdf42',zdf42)
    CALL add_stream_element (debugs,'zdf43',zdf43)
    CALL add_stream_element (debugs,'zdf44',zdf44)
    CALL add_stream_element (debugs,'zdf45',zdf45)
    CALL add_stream_element (debugs,'zdf46',zdf46)
    CALL add_stream_element (debugs,'zdf47',zdf47)
    CALL add_stream_element (debugs,'zdf48',zdf48)
    CALL add_stream_element (debugs,'zdf49',zdf49)
    CALL add_stream_element (debugs,'zdf50',zdf50)
    CALL add_stream_element (debugs,'zdf51',zdf51)
    CALL add_stream_element (debugs,'zdf52',zdf52)
    CALL add_stream_element (debugs,'zdf53',zdf53)
    CALL add_stream_element (debugs,'zdf54',zdf54)
    CALL add_stream_element (debugs,'zdf55',zdf55)
    CALL add_stream_element (debugs,'zdf56',zdf56)
    CALL add_stream_element (debugs,'zdf57',zdf57)
    CALL add_stream_element (debugs,'zdf58',zdf58)
    CALL add_stream_element (debugs,'zdf59',zdf59)
    CALL add_stream_element (debugs,'zdf60',zdf60)
    CALL add_stream_element (debugs,'zdf61',zdf61)
    CALL add_stream_element (debugs,'zdf62',zdf62)
    CALL add_stream_element (debugs,'zdf63',zdf63)
    CALL add_stream_element (debugs,'zdf64',zdf64)

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
