#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE ionwp

  ! Description:
  ! Read initial data for a NWP restart
  !
  ! Method:
  ! An abstraction layer is used to access netCDF data files and
  ! retrieve data.
  ! Further information is contained in a note on the IO package
  !
  ! Authors:
  ! D. Klocke, MPI, March 2010, modified from ioinitial/iorestart
  !

  USE mo_kind,            ONLY: dp
  USE mo_io,              ONLY: ini_spec, io_file_id           &
                              , io_read, io_get_var_double    &
                              , io_open_unit, io_inq_varid, io_read_streams
  USE mo_mpi,             ONLY: p_io, p_pe
  USE mo_control,         ONLY: nlevp1, lvctch, nisp, nlev
  USE mo_memory_sp,       ONLY: sd, sp, stp, su0, svo
  USE mo_memory_gl,       ONLY: xt
  USE mo_memory_g1a,      ONLY: xtm1
  USE mo_memory_base,     ONLY: memory_info, get_stream_element_info
  USE mo_tracer_processes,ONLY: xt_initialize
  USE mo_decomposition,   ONLY: global_decomposition
  USE mo_transpose,       ONLY: scatter_sp
  USE mo_filename,        ONLY: NETCDF
  USE mo_jsbach_interface,ONLY: jsbach_restart
  USE mo_hyb,             ONLY: apsurf
  IMPLICIT NONE

  !  Local scalars: 

  ! number of codes read from surface initialfile

  INTEGER :: nsvoid, nsdid, nstpid, jlev

  REAL(dp), POINTER :: zin(:,:,:), zsu0(:,:), zin1(:,:,:), zin2(:,:,:)

  TYPE (memory_info) :: info

  ! Executable Statements

  ! Initial file information already read in initialize (CALL IO_init)

  ! skip if column model runs with changed hybrid levels

  ! read from restart file, same as in iorestart

  CALL IO_read_streams

  CALL xt_initialize (xt, xtm1)

  CALL jsbach_restart

  IF (.NOT.lvctch) THEN

  IF (p_pe == p_io) THEN

    ini_spec%format = NETCDF
    CALL IO_open_unit(nisp, ini_spec, IO_READ)

    IO_file_id = ini_spec%file_id

    ! 1. Process spectral files

    CALL IO_INQ_VARID (IO_file_id, 'SVO', nsvoid)
    CALL get_stream_element_info (sp, 'svo', info)
    ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
    CALL IO_GET_VAR_DOUBLE (IO_file_id, nsvoid, zin)

  END IF

  CALL scatter_sp (zin, svo, global_decomposition)

  ! 1.1 derive su0 from svo, saves one transpose operation and is 
  !     fast enough

  IF (p_pe == p_io) THEN

    CALL get_stream_element_info (sp, 'su0', info)
    ALLOCATE (zsu0(info%gdim(1), info%gdim(2)))
    CALL init_su0 (zin, zsu0)      

  END IF

  CALL scatter_sp (zsu0, su0, global_decomposition)

  ! finish setup of svo and calculation of su0  

  IF (p_pe == p_io) DEALLOCATE (zin, zsu0)

  IF (p_pe == p_io) THEN

    CALL IO_INQ_VARID (IO_file_id, 'SD', nsdid)
    CALL get_stream_element_info (sp, 'sd', info)
    ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
    CALL IO_GET_VAR_DOUBLE (IO_file_id, nsdid, zin)

  END IF

  CALL scatter_sp (zin, sd, global_decomposition)

  IF (p_pe == p_io) DEALLOCATE (zin)


  ! read ST from nwp restart file (linked to unit.23)

  IF (p_pe == p_io) THEN

    CALL IO_INQ_VARID (IO_file_id, 'ST', nstpid)
    CALL get_stream_element_info (sp, 'st', info)
    ALLOCATE (zin1(info%gdim(1), info%gdim(2), info%gdim(3)))
    CALL IO_GET_VAR_DOUBLE (IO_file_id, nstpid, zin1)

  END IF

  ! read LSP from nwp restart file (linked to unit.23)

  IF (p_pe == p_io) THEN

    CALL IO_INQ_VARID (IO_file_id, 'LSP', nstpid)
    CALL get_stream_element_info (sp, 'lsp', info) 
    ALLOCATE (zin2(info%gdim(1), info%gdim(2), info%gdim(3)))
    CALL IO_GET_VAR_DOUBLE (IO_file_id, nstpid, zin2)

  END IF

  ! combine ST and LSP to STP, as read for initial runs

  IF (p_pe == p_io) THEN
    CALL get_stream_element_info (sp, 'stp', info) 
    ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
    DO jlev = 1, nlev
    zin(jlev,:,:)=zin1(jlev,:,:)
    ENDDO
    zin(nlevp1,:,:) = zin2(1,:,:)
    zin(nlevp1,1,1) = LOG(apsurf-286._dp) ! 286 Pa is difference between 
  END IF

  IF (p_pe == p_io) DEALLOCATE (zin1)
  IF (p_pe == p_io) DEALLOCATE (zin2)

  CALL scatter_sp (zin, stp, global_decomposition)

  IF (p_pe == p_io) DEALLOCATE (zin)

  ENDIF

  RETURN

CONTAINS

  SUBROUTINE init_su0 (psvo, psu0)

    ! Description:
    !
    ! Compute initial spectral components for the zonal mean wind used
    ! in the linearization of the vorticity and humidity equations.
    !
    ! Method:
    !
    ! This subroutine computes initial spectral components
    ! of the mean zonal wind used in the semi-implicit treatment of
    ! vorticity and humidity equations from the vorticity zonal
    ! spectral components.
    !
    ! *inisu0* is called from *ioinitial* / *ionwp*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, February 1983, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_kind,       ONLY: dp
    USE mo_control,    ONLY: nlev, nn, nnp1
    USE mo_physical_constants,  ONLY: earth_radius

    IMPLICIT NONE

    REAL(dp), POINTER :: psvo(:,:,:), psu0(:,:)

    !  Local scalars: 
    REAL(dp) :: zeps1, zeps2, zn
    INTEGER :: jlev, jn

    !  Intrinsic functions 
    INTRINSIC SQRT


    !  Executable statements 

    !-- 1. Set up *su0*

    DO jlev = 1, nlev

      zeps2 = earth_radius/SQRT(3._dp)
      psu0(jlev,1) = zeps2*psvo(jlev,1,2)

      DO jn = 2, nn
        zeps1 = zeps2
        zn = 4.0_dp*jn*jn-1.0_dp
        zeps2 = earth_radius/SQRT(zn)
        psu0(jlev,jn) = -zeps1*psvo(jlev,1,jn-1)+zeps2*psvo(jlev,1,jn+1)
      END DO

      zeps1 = zeps2

      psu0(jlev,nnp1) = -zeps1*psvo(jlev,1,nn)
    END DO

  END SUBROUTINE init_su0

END SUBROUTINE ionwp
