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
SUBROUTINE ioinitial

  ! Description:
  !
  ! Read initial data
  !
  ! Method:
  !
  ! An abstraction layer is used to access netCDF data files and
  ! retrieve data.   
  ! 
  ! Further information is contained in a note on the IO package
  !
  ! Authors:
  !
  ! U. Schulzweida, MPI, May 1999, original version
  ! L. Kornblueh, MPI, April 2010, update scatter routines calls and 
  !       add support for reading of multiple processes  
  !

  USE mo_kind,             ONLY: dp
  USE mo_mpi,              ONLY: p_io, p_pe, p_bcast
  USE mo_tr_gather,        ONLY: gather_field
  USE mo_tr_scatter,       ONLY: scatter_field, scatter_spectral
  USE mo_exception,        ONLY: message, message_text
  USE mo_control,          ONLY: lvctch, lcolumn, lindependent_read,      &
                                 ldebugio, lfractional_mask,              &
                                 nigp, nisp,                              & 
                                 nlev, nlevp1, ngl, nlon, nn, nnp1 
  USE mo_memory_base,      ONLY: get_stream_element, memory_info,         &
                                 get_stream_element_info
  USE mo_memory_sp,        ONLY: sd, sp, stp, su0, svo
  USE mo_memory_gl,        ONLY: gl, q, xl, xi, xt
  USE mo_memory_g1a,       ONLY: xtm1
  USE mo_memory_g3a
  USE mo_memory_g3b
  USE mo_util_string,      ONLY: toupper
  USE mo_filename,         ONLY: NETCDF  
  USE mo_io,               ONLY: ini_spec, ini_surf, io_file_id,          &
                                 io_read, io_var_id, io_get_var_double,   &
                                 io_open_unit, io_inq_varid, io_close
  USE mo_tracer_processes, ONLY: xt_initialize
  USE mo_hyb,              ONLY: apsurf
  USE mo_physical_constants, ONLY: tmelt, earth_radius
  USE mo_decomposition,    ONLY: ldc => local_decomposition
  USE mo_io,               ONLY: slm_glob
  USE mo_gaussgrid,        ONLY: gl_budw

  IMPLICIT NONE

  !  Local scalars: 

  ! number of codes read from surface initialfile
  INTEGER, PARAMETER :: nrec_surf = 12
  CHARACTER(len=8)   :: cname, csurf(nrec_surf)

  INTEGER :: nsvoid, nsdid, nstpid, nqid

  INTEGER :: irec

  REAL(dp), POINTER :: zin(:,:,:), zsu0(:,:), zptr(:,:,:)

  TYPE (memory_info) :: info

  LOGICAL :: l_read = .FALSE.

  !  Executable Statements

  ! Initial file information already read in initialize (CALL IO_init)

  ! skip if column model runs with changed hybrid levels

  IF (.NOT.lvctch) THEN

    IF (p_pe == p_io .OR. lindependent_read) THEN
      l_read = .TRUE.
    ENDIF

    ! 1. process upper air data

    IF (l_read) THEN
      ini_spec%format = NETCDF
      CALL IO_open_unit(nisp, ini_spec, IO_READ)
      IO_file_id = ini_spec%file_id
    ENDIF

    CALL get_stream_element_info (sp, 'svo', info)
    IF (l_read) THEN
      CALL IO_INQ_VARID (IO_file_id, 'SVO', nsvoid)
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nsvoid, zin)
    END IF
    
    CALL scatter_spectral(zin, svo, lindependent_read)

    ! 1.1 derive su0 from svo, saves one transpose operation and is 
    !     fast enough

    CALL get_stream_element_info (sp, 'su0', info)
    IF (l_read) THEN
      ALLOCATE (zsu0(info%gdim(1), info%gdim(2)))
      CALL init_su0 (zin, zsu0)      
    END IF

    CALL scatter_spectral(zsu0, su0, lindependent_read) 

    ! finish setup of svo and calculation of su0  

    IF (l_read) DEALLOCATE (zin, zsu0)

    CALL get_stream_element_info (sp, 'sd', info)
    IF (l_read) THEN
      CALL IO_INQ_VARID (IO_file_id, 'SD', nsdid)
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nsdid, zin)
    END IF

    CALL scatter_spectral(zin, sd, lindependent_read)

    IF (l_read) DEALLOCATE (zin)

    CALL get_stream_element_info (sp, 'stp', info)
    IF (l_read) THEN
      CALL IO_INQ_VARID (IO_file_id, 'STP', nstpid)
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nstpid, zin)
      
      ! Set global mean of surface pressure for initial field: stp(nlevp1,1,1)
      ! 286 Pa is difference between spectral and gridpoint mean. apsurf sets 
      ! the grid point mean
      zin(nlevp1,1,1) = LOG(apsurf-286.0_dp) 
    END IF
    
    CALL scatter_spectral(zin, stp, lindependent_read)

    IF (l_read) DEALLOCATE (zin)

    ! 2. Read grid point data, advected by SL

    CALL get_stream_element_info (gl, 'q', info)
    IF (l_read) THEN
      CALL IO_INQ_VARID (IO_file_id, 'Q', nqid)
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nqid, zin)
    END IF
    
    CALL scatter_field(zin, info%gdim, q, lindependent_read)

    IF (l_read) DEALLOCATE (zin)

    IF (l_read) THEN
      CALL IO_close(ini_spec)
    ENDIF
    
  ELSE          ! column model runs with changed hybrid levels
    q(:,:,:) = 0.0_dp
  ENDIF

  ! 2.1 cloud water and cloud ice set initial to zero
  
  xl(:,:,:) = 0.0_dp
  xi(:,:,:) = 0.0_dp

  ! 2.2 initialize optional tracer fields

  CALL xt_initialize(xt, xtm1)

  ! 3. Prepare grid point surface fields

  IF (l_read) THEN
    ini_surf%format = NETCDF
    CALL IO_open_unit(nigp, ini_surf, IO_READ)
    IO_file_id = ini_surf%file_id
  ENDIF

  ! Codes read from surface initialfile (unit:24)

  csurf( 1) = 'geosp'    ! Surface geopotential
  csurf( 2) = 'slm'      ! Land sea mask (1/0)
  csurf( 3) = 'glac'     ! Glacier mask
  csurf( 4) = 'alake'    ! Lake mask
  csurf( 5) = 'oromea'   ! Mean orography (m)
  csurf( 6) = 'orostd'   ! Orographic standard deviation (m)
  csurf( 7) = 'orosig'   ! Orographic slope
  csurf( 8) = 'orogam'   ! Orographic anisotropy
  csurf( 9) = 'orothe'   ! Orographic angle
  csurf(10) = 'oropic'   ! Orographic peak elevation (m)
  csurf(11) = 'oroval'   ! Orographic valley elevation (m)
  csurf(12) = 'slf'      ! fractional land sea mask

  DO irec = 1, nrec_surf

    cname = csurf(irec)

    IF (l_read) THEN
      IF (ldebugio) THEN
        CALL message('IO_initial','Read '//TRIM(csurf(irec)))
      ENDIF
      CALL IO_INQ_VARID (IO_file_id, toupper(cname), IO_var_id)
    ENDIF

    CALL get_stream_element_info (g3b, cname, info)

    IF (l_read) THEN
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, zin)
    END IF

    CALL get_stream_element (g3b, cname, zptr)

    CALL scatter_field(zin, info%gdim, zptr, lindependent_read)

    IF (l_read)  DEALLOCATE (zin)

  END DO

  IF (l_read) CALL IO_close (ini_surf)

  CALL init_g3()

  CALL message('','')
  CALL message('',' File setup done and arrays allocated.')
  CALL message('','')

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
    ! inisu0 is called from ioinitial.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, February 1983, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    IMPLICIT NONE

    REAL(dp), POINTER :: psvo(:,:,:), psu0(:,:)

    !  Local scalars: 
    REAL(dp) :: zeps1, zeps2, zn
    INTEGER :: jlev, jn

    !  Executable statements 

    !-- 1. Set up su0

    DO jlev = 1, nlev

      zeps2 = earth_radius/SQRT(3.0_dp)
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

  SUBROUTINE init_g3 ()

    !
    ! init_g3a - initialize parameterisation scheme data.
    !
    ! J. K. Gibson, ECMWF, April 1983
    !
    ! Purpose: to prepare the g3a work buffer from an netCDF initial file.
    !
    ! Method: Initial values are set for appropriate variables.
    !

    IMPLICIT NONE
    !  Local array

    REAL(dp), POINTER :: zslf(:,:) => NULL()
    REAL(dp)          :: zslf_sum(ngl)

    !  Local scalars: 

    INTEGER :: jl, jg, jlat

    !  Executable Statements

    ! Make non fractional land sea mask comprising all grid cells with
    ! land fraction (overwrites slm from initial file, necessary for JSBACH)

    IF (lfractional_mask) THEN
       WHERE (slf(:,:) > 0._dp)
         slm(:,:) = 1._dp
       ELSEWHERE
         slm(:,:) = 0._dp
       ENDWHERE
    END IF

    ! global land fraction for water budget correction
    ! used for ocean coupling only, may be skipped in column mode

    IF (.NOT. lcolumn) THEN 

      IF (p_pe == p_io) THEN
        ALLOCATE (zslf(nlon,ngl))
      END IF

      CALL gather_field (zslf, slf)

      IF (p_pe == p_io) THEN
        DO jlat = 1, ngl
          zslf_sum(jlat) = SUM(zslf(1:nlon,jlat))*gl_budw(jlat)
        END DO

        slm_glob = SUM(zslf_sum)

        DEALLOCATE (zslf)

      END IF

      CALL p_bcast (slm_glob,  p_io)

    ELSE
      slm_glob = 0.0_dp    ! set to some value to be written to the rerun file
    END IF

    ! Initialize *g3a* variables not read

    aprlm(:,:)     = 0.0_dp
    aprcm(:,:)     = 0.0_dp
    aprsm(:,:)     = 0.0_dp
    sradsm(:,:)    = 0.0_dp
    tradsm(:,:)    = 0.0_dp
    srad0m(:,:)    = 0.0_dp
    trad0m(:,:)    = 0.0_dp
    vdism(:,:)     = 0.0_dp
    ustrm(:,:)     = 0.0_dp
    vstrm(:,:)     = 0.0_dp
    ahfsm(:,:)     = 0.0_dp
    evapm(:,:)     = 0.0_dp
    ahflm(:,:)     = 0.0_dp
    wind10m(:,:)   = 0.0_dp
    ustrgwm(:,:)   = 0.0_dp
    vstrgwm(:,:)   = 0.0_dp
    vdisgwm(:,:)   = 0.0_dp
    temp2m(:,:)    = 0.0_dp
    dew2m(:,:)     = 0.0_dp
    u10m(:,:)      = 0.0_dp
    v10m(:,:)      = 0.0_dp
    tsurfm(:,:)    = 0.0_dp
    srad0um(:,:)   = 0.0_dp
    tradsum(:,:)   = 0.0_dp
    sradsum(:,:)   = 0.0_dp
    t2maxm(:,:)    = 0.0_dp
    t2minm(:,:)    = 999.0_dp
    wimaxm(:,:)    = 0.0_dp
    topmaxm(:,:)   = 99999.0_dp
    aclcvm(:,:)    = 0.0_dp
    aclcovm(:,:)   = 0.0_dp
    qvim(:,:)      = 0.0_dp
    xlvim(:,:)     = 0.0_dp
    xivim(:,:)     = 0.0_dp
    wl(:,:)        = 0.0_dp
    siced(:,:)     = 0.0_dp
    sni(:,:)       = 0.0_dp
    gld(:,:)       = 0.0_dp
    rain(:,:)      = 0.0_dp
    acvtype(:,:)   = 0.0_dp
    xtec(:,:,:)    = 0.0_dp
    snc(:,:)       = 0.0_dp
    sswvism(:,:)   = 0.0_dp
    sswparm(:,:)   = 0.0_dp
    sswdifnirm(:,:)= 0.0_dp
    sswdifvism(:,:)= 0.0_dp
    sswdifparm(:,:)= 0.0_dp

    amlcorr(:,:)   = 0.0_dp
    amlcorac(:,:)  = 0.0_dp
    amlheatac(:,:) = 0.0_dp
    
 
    srad0dm(:,:)   = 0.0_dp

    emterm(:,:,:)  = 0.0_dp
    trsolm(:,:,:)  = 0.0_dp

    aclcm(:,:,:)   = 0.0_dp
    aclcacm(:,:,:) = 0.0_dp
    tkem(:,:,:)    = 1.0e-4_dp
    tkem1(:,:,:)   = tkem(:,:,:)
    thvvar(:,:,:)  = 1.0e-4_dp
    thvsig(:,:)    = 1.0e-2_dp

    ! Set variables for fractional surface coverage

    ahfswac(:,:)   = 0.0_dp
    ahfsiac(:,:)   = 0.0_dp
    ahfslac(:,:)   = 0.0_dp
    ahflwac(:,:)   = 0.0_dp
    ahfliac(:,:)   = 0.0_dp
    ahfllac(:,:)   = 0.0_dp
    evapwac(:,:)   = 0.0_dp
    evapiac(:,:)   = 0.0_dp
    evaplac(:,:)   = 0.0_dp
    trfllac(:,:)   = 0.0_dp
    trflwac(:,:)   = 0.0_dp
    trfliac(:,:)   = 0.0_dp
    sofllac(:,:)   = 0.0_dp
    soflwac(:,:)   = 0.0_dp
    sofliac(:,:)   = 0.0_dp
    friac(:,:)     = 0.0_dp
    ustrw(:,:)     = 0.0_dp
    ustri(:,:)     = 0.0_dp
    ustrl(:,:)     = 0.0_dp
    vstrw(:,:)     = 0.0_dp
    vstri(:,:)     = 0.0_dp
    vstrl(:,:)     = 0.0_dp
    alsom(:,:)     = alb(:,:)
    alsobs(:,:)    = alb(:,:)
! set albedo (remove init_surface)
    alsol(:,:)     = 0.2_dp
    alsow(:,:)     = 0.07_dp       ! calbsea
    alsoi(:,:)     = 0.55_dp       ! update_albedo_ice
    albedo_vis(:,:)     = 0.0_dp
    albedo_vis_dir(:,:) = 0.0_dp
    albedo_vis_dif(:,:) = 0.0_dp
    albedo_nir(:,:)     = 0.0_dp
    albedo_nir_dir(:,:) = 0.0_dp
    albedo_nir_dif(:,:) = 0.0_dp
! end set albedo
    ahfice(:,:)    = 0.0_dp   
    qres(:,:)      = 0.0_dp
    ahfcon(:,:)    = 0.0_dp
    ahfres(:,:)    = 0.0_dp
    seaice(:,:)    = 0.0_dp
    tsi(:,:)       = tmelt  ! dummy setting
    tsw(:,:)       = tmelt  ! dummy setting 
    wind10w(:,:)   = 0.0_dp
    rtype(:,:)     = 0.0_dp
    rintop(:,:)    = 0.0_dp
    apmeb(:,:)     = 0.0_dp
    apmebco(:,:)   = 0.0_dp
    qtnew(:,:)     = 0.0_dp
    taus(:,:)      = 0.0_dp
    fage(:,:)      = 0.0_dp
    snifrac(:,:)   = 0.0_dp
    barefrac(:,:)  = 0.0_dp
    !
    !   new variables for melt ponds on sea-ice
    !
    sicepdw(:,:)    = 0.0_dp
    sicepdi(:,:)    = 0.0_dp
    sicepres(:,:)   = 0.0_dp
    sicemin(:,:)    = 1.0_dp
    tsicepdi(:,:)   = tmelt
    ameltdepth(:,:) = 0.0_dp
    ameltfrac(:,:)  = 0.0_dp
    qresmp(:,:)     = 0.0_dp
    aiqre(:,:)      = 0.0_dp
    !
    ! Initialize tropopause height to 200 hPA
    !
    tropo(:,:)      = 20000.0_dp
    !
    ! variables for exchange with ocean model
    !
    ocu(:,:)        = 0.0_dp
    ocv(:,:)        = 0.0_dp
    ! 
    ! sulfate aerosols
    !
    abso4(:,:)      = 0.0_dp
    so4nat(:,:,:)   = 0.0_dp
    so4all(:,:,:)   = 0.0_dp
    !
    ! Setting of array of variable soil characteristics to be use
    ! in *surf*
    ! Input: FAO soils interpolated from 0.5 degree resolution
    !        to model resolution (simple average).
    ! and setting of array of variable available water storage capacity
    ! to be used in *surf*
    ! Input: Patterson data interpolated from 0.5 degree resolution
    !        to model resolution (simple average).
    !
    DO jg = 1, ldc%ngpblks
      DO jl = 1, ldc%nproma
        IF (NINT(faom(jl,jg)) == 1) THEN
          rgcgnm(jl,jg) = 1.93e+06_dp
          sodifm(jl,jg) = 8.7e-7_dp
        ELSE IF (NINT(faom(jl,jg)) == 2) THEN
          rgcgnm(jl,jg) = 2.10e+06_dp
          sodifm(jl,jg) = 8.0e-7_dp
        ELSE IF (NINT(faom(jl,jg)) == 3) THEN
          rgcgnm(jl,jg) = 2.25e+06_dp
          sodifm(jl,jg) = 7.4e-7_dp
        ELSE IF (NINT(faom(jl,jg)) == 4) THEN
          rgcgnm(jl,jg) = 2.36e+06_dp
          sodifm(jl,jg) = 7.1e-7_dp
        ELSE IF (NINT(faom(jl,jg)) == 5) THEN
          rgcgnm(jl,jg) = 2.48e+06_dp
          sodifm(jl,jg) = 6.7e-7_dp
        ELSE
          IF (NINT(faom(jl,jg)) == 0) THEN
            rgcgnm(jl,jg) = 2.25e+06_dp
            sodifm(jl,jg) = 7.4e-7_dp
          ELSE
            WRITE (message_text,'(a,i0,a,i0,a,f3.0)') &
                 'faom(', jl, ',', jg, ') = ', faom(jl,jg)
            CALL message('',message_text)
          END IF
        END IF
        ws(jl,jg) = MIN(ws(jl,jg),wsmx(jl,jg))
      END DO
    END DO

  END SUBROUTINE init_g3

END SUBROUTINE ioinitial
