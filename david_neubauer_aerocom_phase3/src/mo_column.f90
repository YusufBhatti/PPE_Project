!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_column
  !> Main driver module for ECHAM Single Column Model
  !>
  !>
  !> Authors:
  !>
  !> A. Rhodin,      MPI, November  1999, original source
  !> I. Kirchner,    MPI, December  2000, time control update
  !> A. Rhodin,      DWD, March     2003, nproma blocking
  !> L. Kornblueh,   MPI, August    2009, bug fix with respect to pole filter
  !> S. K. Cheedela, MPI, January   2010, rewritten
  !>
  USE mo_kind,          ONLY: wp
  USE mo_mpi,           ONLY: p_pe, p_io, p_parallel_io                                
  USE mo_io_units,      ONLY: nout
  USE mo_exception,     ONLY: finish,message
  USE mo_time_control,  ONLY: lstart, delta_time
  USE mo_control,       ONLY: vct, nvclev, nlevp1, nlev
  USE mo_iocolumn

  USE mo_read_netcdf77

  IMPLICIT NONE

  PRIVATE 
   
  PUBLIC :: setup_column
  PUBLIC :: init_column
  PUBLIC :: get_col_states
  PUBLIC :: get_col_div
  PUBLIC :: get_col_omega
  PUBLIC :: get_col_dyns
  PUBLIC :: get_col_trans
  PUBLIC :: update_column
  PUBLIC :: column_diagnostics_ls
  PUBLIC :: column_diagnostics_final
  PUBLIC :: resetcolumn
  PUBLIC :: nfor_div
  PUBLIC :: nfor_omega
  PUBLIC :: nfor_ts
  ! Forcing flags for states and tendencies
  !> variable(tendencyswitch,relaxswitch,extrapolation)
  !> Default is free run 
  !> states with tendencies
  INTEGER  :: nfor_t(3)   =  (/0,0,0/)
  INTEGER  :: nfor_ps(3)  =  (/0,0,0/)
  INTEGER  :: nfor_uv(3)  =  (/0,0,0/)
  INTEGER  :: nfor_q(3)   =  (/0,0,0/)

  !> states with no tendencies
  INTEGER  :: nfor_uvgeo(2)  =  (/0,0/)
  INTEGER  :: nfor_omega(2)  =  (/0,0/)
  INTEGER  :: nfor_div(2)    =  (/0,0/)
  INTEGER  :: nfor_ts(2)  =  (/0,0/)
  INTEGER  :: nfor_lhf(2) =  (/0,0/)
  INTEGER  :: nfor_shf(2) =  (/0,0/)
    
  !> additional tracers
  INTEGER  :: nfor_xi(3)  =  (/0,0,0/)
  INTEGER  :: nfor_xl(3)  =  (/0,0,0/)
  !INTEGER  :: nfor_xt(3)  =  (/0,0,0/)     !> currently not yet implemented

  INTEGER,PUBLIC   :: lat_1d(1)  = -1       !> nearest latitude longitude index 
                                            !> in current resolution
  INTEGER,PUBLIC   :: lon_1d(1)  = -1
  REAL(wp),PUBLIC  :: dlat(1)    = -100._wp !> latitude longitude in degrees
  REAL(wp),PUBLIC  :: dlon(1)    = -1._wp         
  REAL(wp),PUBLIC  :: mld        = 10._wp
  REAL(wp)         :: zqcst
  REAL(wp)         :: zcorio
  !> geostrophic winds 
  !REAL(wp)     :: ug      =  0._wp
  !REAL(wp)     :: vg      =  0._wp
  REAL(wp),PUBLIC,ALLOCATABLE :: sst_1d(:,:)
  REAL(wp),PUBLIC,ALLOCATABLE :: sst_1d_input(:)
  INTEGER       :: tendency_optn = 0
  
  LOGICAL,PUBLIC       :: ml_input =  .FALSE.
  INTEGER              :: ierr_ff

  INCLUDE 'columnctl.inc' 

CONTAINS
 
  SUBROUTINE setup_column(lcolumn)
  USE mo_advection,       ONLY:iadvec, no_advection
  USE mo_mpi,             ONLY:p_nprocs
  USE mo_namelist,        ONLY:open_nml, position_nml, POSITIONED
  USE mo_netcdf,          ONLY:io_init_dims
  LOGICAL ,INTENT(in)      :: lcolumn    ! .True. for Single column model
  INTEGER                  :: ncol ,ierr, inml, iunit 
 
  !> Read namelist  
  IF (p_parallel_io) THEN
    inml = open_nml('namelist.echam')
    iunit = position_nml ('COLUMNCTL', inml, status=ierr)
    SELECT CASE (ierr)
    CASE (POSITIONED)
      READ (iunit, columnctl)
    END SELECT
  ENDIF
 
  CALL init_forcing

  !> Only for column model
  IF(lcolumn) THEN
    IF(p_pe == p_io) THEN

    iadvec  = no_advection  !> switch off advection
    
    WRITE(nout,'(a)') REPEAT('*',72)
    WRITE(nout,*) 'You are using ECHAM Single Column Model'
    WRITE(nout,'(a)') REPEAT('-',72)
    
    IF(p_nprocs>1) CALL finish('','Single column model cannot use more than 1 processor')
   
    CALL get_col_location(dlat,dlon) !> read lat,lon from forcing file in degrees
   
    IF (dlat(1)>90_wp .OR. dlat(1)<-90_wp .OR. dlon(1)>360_wp .OR. dlon(1)<0_wp) THEN
     CALL finish('latitude(deg N) and longitude(deg E)','Out of range in the forcing file')
    END IF

    WRITE(nout,*) 'latitude(deg N) and logitude(deg E) for single column model are', &
                  & REAL(dlat),REAL(dlon)
 
    CALL cal_latlon_index(dlat,dlon,lat_1d,lon_1d) !> compute nearest lat,lon index on
                                                   !> current grid 

    ncol = SIZE(lat_1d)*SIZE(lon_1d) !> Check size and limits of lat and lon

    IF(lat_1d(1) < 0 .OR. lon_1d(1) < 0 .OR. ncol >1)   &
                  & CALL finish('setup_column','Check your grid, wrong index')

    WRITE(nout,*) 'latitude and longitude indices of nearest grid point are',lat_1d,lon_1d

    CALL get_col_vct(nlev, nlevp1, nvclev, vct) !> always read and change vertical 
                                                !> levels from forcingfile
    !> initialize dimensions
    CALL io_init_dims

    IF(nfor_omega(1)>0 .AND. nfor_div(1)>0 )            &
                  & CALL finish('only one of nfor_omega,nfor_div','can be greater than 1')

    IF(nfor_omega(1)>0 .OR. nfor_div(1)>0 ) tendency_optn=1 

    IF((nfor_uv(2).GT.0) .AND.(nfor_uvgeo(1).GT.0))     &
                  & CALL finish('only one of nfor_uv and nfor_uvgeo','can be greater than 1')
    ENDIF
  ENDIF
    
  END SUBROUTINE setup_column

  SUBROUTINE init_column
  !> Purpose is to initialize column model 
  !> called after ioinitialize and before init_surface
  !> initialize initial atm fields --essential T,Q,U,V ---optional  Xl,Xi
  USE mo_memory_gl,   ONLY:q,xl,xi
  USE mo_memory_g3b,  ONLY:g3b,slm,slf,alake
  USE mo_memory_base, ONLY:get_stream_element
  ! local 
  CHARACTER (len=10)  :: cname , csurf(16)
  INTEGER             :: ivar
  REAL(wp) , POINTER  :: zptr(:,:)
 
  IF (p_pe == p_io) THEN
    !> Atmospheric States
    
    CALL read_initial3('t',scm_t,'essential')
    CALL read_initial3('u',scm_u,'essential')
    CALL read_initial3('v',scm_v,'essential')
    CALL read_initial3('q',q,'essential')
    CALL read_initial2('aps',scm_aps,'essential') 
    CALL read_initial3('xl',xl,'optional')
    CALL read_initial3('xi',xi,'optional')
    !CALL read_initial3('Xt',xt,'optional') !> yet to implement for tracers
 
    !> Optional surface fields 
    !> Surface fields from g3b
    csurf( 1) = 'geosp'    !> Surface geopotential
    csurf( 2) = 'ws'       !> Surface soil wetness
    csurf( 3) = 'sn'       !> Snow depth
    csurf( 4) = 'slm'      !> Land sea mask (1/0)
    csurf( 5) = 'alb'      !> Albedo
    csurf( 6) = 'forest'   !> Vegetation type
    csurf( 7) = 'wsmx'     !> Field capacity of soil
    csurf( 8) = 'fao'      !> Fao data set
    csurf( 9) = 'glac'     !> Glacier mask
    csurf(10) = 'oromea'   !> Mean orography (m)
    csurf(11) = 'orostd'   !> Orographic standard deviation (m)
    csurf(12) = 'orosig'   !> Orographic slope
    csurf(13) = 'orogam'   !> Orographic anisotropy
    csurf(14) = 'orothe'   !> Orographic angle
    csurf(15) = 'oropic'   !> Orographic peak elevation (m)
    csurf(16) = 'oroval'   !> Orographic valley elevation (m)
   
    DO ivar=1,size(csurf) 
      cname = csurf(ivar) 
      CALL get_stream_element (g3b, cname, zptr )
      CALL read_initial2(cname,zptr,'optional')
    END DO 
    
    !> Surface fields from jsbach surface memory
    !> Pending list
   
    ! Set default surface parameters mainly for ocean , sice unclear
    IF (nfor_ts(1)>0 .AND. slm(1,1) <1 )  THEN !> force ocean with surface temp 
      alake = 0._wp 
      slf = slm 
      ALLOCATE(sst_1d(1,1))
    ELSE !> interactive ocean ,surface temperature calculated by lake model
      alake = 1._wp -slm 
      slf = slm  
    END IF
   
    IF (ml_input) THEN
      ALLOCATE(sst_1d_input(4))
    ENDIF
 
    CALL message('Single Column Model','was sucessfully initialized')
    WRITE(nout,'(a)') REPEAT('-',72)
  
  ENDIF 
 
  END SUBROUTINE init_column
 
  SUBROUTINE get_col_states
  USE mo_memory_gl,   ONLY: q,xl,xi
  USE mo_scan_buffer, ONLY: u,v,t,alps,ugeo,vgeo
  USE mo_gaussgrid,   ONLY: gl_sqcst,gl_coriol
  EXTERNAL pres
 
  IF (lstart) THEN !> initialization before dyn
    zqcst=gl_sqcst(lat_1d(1))
    zcorio=gl_coriol(lat_1d(1))
    t= scm_t
    u= zqcst*scm_u
    v= zqcst*scm_v
    alps = LOG(scm_aps)
    CALL pres(pressure(:,:,1),1,EXP(alps(:,1)),1) 
  END IF
    
  IF(nfor_t(2)>0)  CALL get_ext_state("t",t,nfor_t(2),nfor_t(3))
  IF(nfor_uv(2)>0) CALL get_ext_state("u",u,nfor_uv(2),nfor_uv(3))
  IF(nfor_uv(2)>0) CALL get_ext_state("v",v,nfor_uv(2),nfor_uv(3))
 
  IF(nfor_uvgeo(1)>0) CALL get_ext_state("u",ugeo,nfor_uvgeo(1),nfor_uvgeo(2))
  IF(nfor_uvgeo(1)>0) CALL get_ext_state("v",vgeo,nfor_uvgeo(1),nfor_uvgeo(2))
 
  IF(nfor_q(2)>0)  CALL get_ext_state("q",q,nfor_q(2),nfor_q(3))
  IF(nfor_xi(2)>0) CALL get_ext_state("xi",xi,nfor_xi(2),nfor_xi(3))
  IF(nfor_xl(2)>0) CALL get_ext_state("xl",xl,nfor_xl(2),nfor_xl(3))
    
  IF(nfor_ps(2)>0) THEN     !> need to call before dyn.f90 
    CALL get_ext_state("ps",scm_aps,nfor_ps(2),nfor_ps(3))
    alps = LOG(scm_aps)
  END IF
 
  !> Surface parameters 
 
  ! sets the surface-temperature and the sea-surface-temperature to the values of the 
  !  surface temperature 'ts' from forcingfile, if ml_input is set to true  
  IF(lstart) THEN
    IF(ml_input)      CALL read_var_nf77_1d(forcingfile,"time","ts",sst_1d_input,ierr_ff)
  ENDIF
    
  IF(nfor_ts(1)>0)  CALL get_ext_state("ts",sst_1d,nfor_ts(1),nfor_ts(2)) 
    
  !IF(nfor_lhf(1)>0)  CALL get_ext_state("lhf",lhf,nfor_lhf(1),nfor_lhf(2))
  !IF(nfor_shf(2)>0)  CALL get_ext_state("shf",shf,nfor_shf(2),nfor_shf(2)) 
    
  END SUBROUTINE get_col_states
 
  SUBROUTINE get_col_div(d)
  REAL(wp), INTENT(inout) :: d(:,:,:)
  
  CALL get_ext_state("div",scm_div,nfor_div(1),nfor_div(2))
  d = scm_div 
 
  END SUBROUTINE get_col_div
   
  SUBROUTINE get_col_omega(zsdiv,d,delp)
  !> prescribe pressure velocity in SCM works only in pressure coordinates or
  !> sigma??? 
  REAL(wp), INTENT(inout) :: d(:,:,:),zsdiv(:,:)
  REAL(wp), INTENT(in)    :: delp(:,:)
  INTEGER                 :: jk,jkp 
   
  CALL get_ext_state("omega",scm_omega,nfor_omega(1),nfor_omega(2))
  zsdiv(1,:) = 0.0_wp
   
  DO jk =1,nlev-1
    jkp =jk+1
    zsdiv(1,jkp)= (scm_omega(1,jk,1) + scm_omega(1,jkp,1))/2.0
  END DO
 
  zsdiv = -1*zsdiv 
  
  DO jk = 1,nlev
    jkp =jk+1   
    d(1,jk,1)    =(zsdiv(1,jkp)-zsdiv(1,jk))/delp(1,jk)
  END DO
 
  scm_div   = d !> divergence for output

  END SUBROUTINE get_col_omega

  SUBROUTINE get_col_dyns
  USE mo_scan_buffer, ONLY: tte, vom, vol
  !> change of names for ddt_t to _adv or ls required
  IF(nfor_t(1)>0) CALL get_ext_tend("ddt_t",tte,nfor_t(3),tendency_optn)
  IF(nfor_uv(1)>0) CALL get_ext_tend("ddt_u",vom,nfor_uv(3),tendency_optn)
  IF(nfor_uv(1)>0) CALL get_ext_tend("ddt_v",vol,nfor_uv(3),tendency_optn)
  
  END SUBROUTINE get_col_dyns

  SUBROUTINE get_col_trans
  ! tracers are seperated for future use
  USE mo_scan_buffer, ONLY: qte, xlte, xite 
  USE mo_memory_gl,   ONLY: q, xi, xl 
  ! Maintain tracers positive definite; force 0 if negative
  ! done after non monotonic vertical advection
  q = MAX(q,0.0_wp)
  xi= MAX(xi,0.0_wp)
  xl= MAX(xl,0.0_wp)

  IF(nfor_q(1)>0) CALL get_ext_tend("ddt_q",qte,nfor_q(3),tendency_optn)
  IF(nfor_xl(1)>0) CALL get_ext_tend("ddt_ql",xlte,nfor_xl(3),tendency_optn)
  IF(nfor_xi(1)>0) CALL get_ext_tend("ddt_qi",xite,nfor_xi(3),tendency_optn)

  END SUBROUTINE get_col_trans

  SUBROUTINE update_column
  USE mo_scan_buffer, ONLY:u, v, t, vom, vol, tte , ugeo,vgeo
  USE mo_memory_g2a,  ONLY: um1, vm1
  USE mo_memory_g1a,  ONLY: tm1

  IF( nfor_uvgeo(1) .EQ. 0) THEN  ! changed to 2*dt and added coriolis , fxp 01/2012
    u   = um1  + 2*delta_time*vom 
    v   = vm1  + 2*delta_time*vol
    t   = tm1  + 2*delta_time*tte
    ! alps = alpsm1 + 2*delta_time*alpste
  ELSE 
    u   = um1  + 2*delta_time*(vom-zcorio*(vgeo*zqcst-vm1))
    v   = vm1  + 2*delta_time*(vol+zcorio*(ugeo*zqcst-um1))
    t   = tm1  + 2*delta_time*tte
  END IF 
  
  CALL column_diagnostics_final
 
  END SUBROUTINE update_column

  SUBROUTINE column_diagnostics_ls
  USE mo_scan_buffer, ONLY: xlte, xite, qte, vom, vol, tte

  ddt_u_adv  = vom
  ddt_v_adv  = vol
  ddt_t_adv  = tte
  ddt_q_adv  = qte
  ddt_xl_adv = xlte
  ddt_xi_adv = xite

  END SUBROUTINE column_diagnostics_ls

  SUBROUTINE column_diagnostics_final 
  ! added zqcst conversion, fxp 7-11
  !
  USE mo_scan_buffer, ONLY:u, v, t, qte, xlte, xite, &
                      vom, vol, tte, vervel
  scm_u     = u/zqcst
  scm_v     = v/zqcst
  scm_omega = vervel   ! computed vertical velocity
  scm_t     = t
  !Total tendency from physics
  ddt_u_phy = vom-ddt_u_adv
  ddt_v_phy = vol-ddt_v_adv
  ddt_t_phy = tte-ddt_t_adv
  ddt_q_phy = qte-ddt_q_adv
  ddt_xl_phy= xlte-ddt_xl_adv
  ddt_xi_phy= xite-ddt_xi_adv

  END SUBROUTINE column_diagnostics_final

  SUBROUTINE resetcolumn

  CALL destruct_scm

  IF ( ALLOCATED( sst_1d )) DEALLOCATE (sst_1d)
  IF ( ALLOCATED( sst_1d_input)) DEALLOCATE (sst_1d_input)

  END SUBROUTINE resetcolumn
 
  SUBROUTINE cal_latlon_index(dlat,dlon,lat_1d,lon_1d)
  USE mo_gaussgrid, ONLY:inigau,philat,philon
  REAL(wp) , INTENT(in)  :: dlat(1),dlon(1)
  INTEGER  , INTENT(out) :: lat_1d(1),lon_1d(1)
    
  IF (dlat(1)>90_wp .OR. dlat(1)<-90_wp .OR. dlon(1)>360_wp .OR. dlon(1)<0_wp ) THEN 
    CALL finish('init_column','check your latitude(deg N) and longitude(deg E) range')
  END IF
    ! very crude approx only first decimal accurate
    ! dlatlen= 180.0/(ngl-1) !-1 because of non cyclic 90/=-90
    ! dlonlen= 360.0/nlon    
    ! lat_1d = NINT((90.0-dlat(1))/dlatlen)+1 !
    ! lon_1d = NINT(dlon(1)/dlonlen)+1        

    ! below accurate way
  CALL inigau ! check for changing horizontal resolution...??
    
  lat_1d = MINLOC(ABS(philat-dlat(1)))
  lon_1d = MINLOC(ABS(philon-dlon(1)))
 
  END SUBROUTINE cal_latlon_index
  
END MODULE mo_column
