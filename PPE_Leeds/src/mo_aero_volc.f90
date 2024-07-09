!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_aero_volc
!-------------------------------------------------------------------------
!
!    Sebastian Rast, MPI Met, Hamburg, February 2010
!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: mo_aero_volc
!
! !DESCRIPTION:
! Introduction of optical properties of stratospheric volcanic aerosols from 
! transient data set
!
! !REVISION HISTORY:
! original source by J.S.Rast (2010-02-18)
!
! !USES:
  USE mo_kind,                 ONLY: wp
  USE mo_decomposition,        ONLY: ldc=>local_decomposition
  USE mo_rrtm_params,          ONLY: nbndlw
  USE mo_exception,            ONLY: finish

  IMPLICIT NONE

  PRIVATE
  PUBLIC                           :: su_aero_volc, read_aero_volc, &
                                      add_aop_volc, cleanup_aero_volc

! !LOCAL VARIABLES

  REAL(wp), ALLOCATABLE            :: aod_v_s(:,:,:),   & ! volcanic AOD solar
                                      ssa_v_s(:,:,:,:), & ! volcanic SSA solar
                                      asy_v_s(:,:,:,:), & ! volcanic ASY solar
                                      ext_v_s(:,:,:,:), & ! volcanic EXT solar
                                      aod_v_t(:,:,:),   & ! volcanic AOD therm
                                      ssa_v_t(:,:,:,:), & ! volcanic SSA therm
                                      ext_v_t(:,:,:,:), & ! volcanic EXT therm
                                      p_lim_clim(:)    ! limit lev press. data
  INTEGER, PARAMETER               :: lev_clim=40
  LOGICAL                          :: laero_set = .FALSE.

CONTAINS
!EOP
!-------------------------------------------------------------------------
!BOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: su_aero_volc
!
! !SUBROUTINE INTERFACE:
SUBROUTINE su_aero_volc(nb_sw)
!
! !DESCRIPTION:
! set up memory for optical aerosol parameters for new aerosol climatology
! compiled by S. Kinne, called by setrad
!
! !REVISION HISTORY:
! original source by J.S. Rast (2010-02-18)
!
! !USES:
! !INPUT PARAMETERS
  INTEGER, INTENT(in)            :: nb_sw

! !LOCAL VARIABLES

  INTEGER, PARAMETER             :: nmonths=12

! allocate memory for optical properties
  ALLOCATE(aod_v_s(nb_sw,ldc%nlat,0:nmonths+1))
  ALLOCATE(ext_v_s(nb_sw,lev_clim,ldc%nlat,0:nmonths+1))
  ALLOCATE(ssa_v_s(nb_sw,lev_clim,ldc%nlat,0:nmonths+1))
  ALLOCATE(asy_v_s(nb_sw,lev_clim,ldc%nlat,0:nmonths+1))
  ALLOCATE(aod_v_t(nbndlw,ldc%nlat,0:nmonths+1))
  ALLOCATE(ext_v_t(nbndlw,lev_clim,ldc%nlat,0:nmonths+1))
  ALLOCATE(ssa_v_t(nbndlw,lev_clim,ldc%nlat,0:nmonths+1))
  ALLOCATE(p_lim_clim(lev_clim+1))
END SUBROUTINE su_aero_volc
!EOP
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_aero_volc
!
!SUBROUTINE interface
SUBROUTINE read_aero_volc(rad_step_1, rad_step_2, nb_sw)
!
! !DESCRIPTION:
! This routine reads the aerosol data on a yearly basis 
!
! !REVISION HISTORY:
! original source by J.S. Rast (2010-02-18)
! 
! !USES:

  USE mo_time_conversion,    ONLY: time_days
  USE mo_time_control,       ONLY: get_date_components
  
  !INPUT PARMETES
  TYPE(time_days), INTENT(in)   :: rad_step_1, rad_step_2
  INTEGER, INTENT(in)           :: nb_sw

  !LOCAL VARIABLES
  INTEGER                       :: icurrentyear, inextyear
  INTEGER                       :: zyrm1, zyr, zyrp1
  LOGICAL                       :: lnewyear
  CHARACTER(len=20)             :: cfname_base,cyr
  CHARACTER(len=25)             :: cfname

  CALL get_date_components (rad_step_1, year=icurrentyear)
  CALL get_date_components (rad_step_2, year=inextyear)
  lnewyear=icurrentyear/=inextyear

  IF (.NOT.lnewyear .AND. laero_set) return
  zyrm1=inextyear-1
  zyr=inextyear
  zyrp1=inextyear+1
! (stratospheric) volcanic aerosol, solar radiation
  cfname_base='strat_aerosol_sw'
  WRITE(cyr,*) zyrm1
  cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
  CALL aero_read_opt_volc_sw ( &
    ldc%nlat,             nb_sw,         lev_clim,                      &
    aod_v_s(:,:,0:0),     'aod',         ssa_v_s(:,:,:,0:0),    'ssa',  &
    asy_v_s(:,:,:,0:0),   'asy',         ext_v_s(:,:,:,0:0),    'ext',  &
    p_lim_clim,           'p_lim',       12,                    12,     &
    cfname                                                              )
  WRITE(cyr,*) zyr
  cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
  CALL aero_read_opt_volc_sw ( &
    ldc%nlat,             nb_sw,         lev_clim,                      &
    aod_v_s(:,:,1:12),    'aod',         ssa_v_s(:,:,:,1:12),   'ssa',  &
    asy_v_s(:,:,:,1:12),  'asy',         ext_v_s(:,:,:,1:12),   'ext',  &
    p_lim_clim,           'p_lim',       1,                     12,     &
    cfname                                                              )
  WRITE(cyr,*) zyrp1
  cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
  CALL aero_read_opt_volc_sw ( &
    ldc%nlat,             nb_sw,         lev_clim,                       &
    aod_v_s(:,:,13:13),   'aod',         ssa_v_s(:,:,:,13:13),   'ssa',  &
    asy_v_s(:,:,:,13:13), 'asy',         ext_v_s(:,:,:,13:13),   'ext',  &
    p_lim_clim,           'p_lim',       1,                      1,      &
    cfname                                                               )
! (stratospheric) volcanic aerosol, IR radiation
  cfname_base='strat_aerosol_ir'
  WRITE(cyr,*) zyrm1
  cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
  CALL aero_read_opt_volc_ir ( &
    ldc%nlat,             nbndlw,         lev_clim,                     &
    aod_v_t(:,:,0:0),     'aod',         ssa_v_t(:,:,:,0:0),    'ssa',  &
                                         ext_v_t(:,:,:,0:0),    'ext',  &
    p_lim_clim,           'p_lim',       12,                    12,     &
    cfname                                                              )
  WRITE(cyr,*) zyr
  cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
  CALL aero_read_opt_volc_ir ( &
    ldc%nlat,             nbndlw,         lev_clim,                     &
    aod_v_t(:,:,1:12),    'aod',         ssa_v_t(:,:,:,1:12),   'ssa',  &
                                         ext_v_t(:,:,:,1:12),   'ext',  &
    p_lim_clim,           'p_lim',       1,                     12,     &
    cfname                                                              )
  WRITE(cyr,*) zyrp1
  cfname=TRIM(cfname_base)//'_'//TRIM(ADJUSTL(cyr))//'.nc'
  CALL aero_read_opt_volc_ir ( &
    ldc%nlat,             nbndlw,         lev_clim,                      &
    aod_v_t(:,:,13:13),   'aod',         ssa_v_t(:,:,:,13:13),   'ssa',  &
                                         ext_v_t(:,:,:,13:13),   'ext',  &
    p_lim_clim,           'p_lim',       1,                      1,      &
    cfname                                                               )
  laero_set = .TRUE.
END SUBROUTINE read_aero_volc
!EOP
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_aop_volc
!
!SUBROUTINE interface
SUBROUTINE add_aop_volc ( &
          & kproma,                 kbdim,              klev,             &
          & krow,                   nb_lw,              nb_sw,            &
          & paer_tau_lw_vr,         paer_tau_sw_vr,     paer_piz_sw_vr,   &
          & paer_cg_sw_vr,          ppd_hl,             pp_fl,            &
          & tk_fl)
!
! !DESCRIPTION:
! add aerosol optical properties for all wave length bands (solar and IR)
! in the case of volcanic aerosols of Stenchikov
! The height profile is taken into account.
!
! !REVISION HISTORY:
! original source by J.S. Rast (2010-02-19)
! 
! !USES:
  USE mo_interpo,       ONLY: nmw1_m, nmw2_m, wgt1_m, wgt2_m
  USE mo_physical_constants,     ONLY: rgrav, rd
  USE mo_geoloc,        ONLY: ilat

! !INPUT PARAMETERS
  INTEGER,INTENT(in)  :: kproma, &! actual block length
                         kbdim,  &! maximum block length
                         krow,   &! block index
                         klev,   &! number of vertical levels
                         nb_lw,  &! number of wave length bands (far IR)
                         nb_sw    ! number of wave length bands (solar)
  REAL(wp),INTENT(in) :: ppd_hl(kbdim,klev)  ,& ! layer pressure thickness 
                         pp_fl(kbdim,klev)   ,& ! pressure at "full levels"
                         tk_fl(kbdim,klev)      ! temperature at "full lev."
! !OUTPUT PARAMETERS
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_lw):: &
   paer_tau_lw_vr      !aerosol optical depth (far IR)
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_sw):: &
   paer_tau_sw_vr,   & !aerosol optical depth (solar), sum_i(tau_i)
   paer_piz_sw_vr,   & !weighted sum of single scattering albedos, 
                       !sum_i(tau_i*omega_i)
   paer_cg_sw_vr       !weighted sum of asymmetry factors, 
                       !sum_i(tau_i*omega_i*g_i)

! !LOCAL VARIABLES
  
  INTEGER                     :: jl,jk,jki,jwl,idx_lat,idx_lev
  REAL(wp), DIMENSION(kbdim,klev)   :: zdeltag    ! layer thickness [m]
  REAL(wp), DIMENSION(kbdim,klev,nb_sw)  :: zext_s, zomg_s, zasy_s
  REAL(wp), DIMENSION(kbdim,nb_sw)       :: zaod_s, zext_s_int, zfact_s
  REAL(wp), DIMENSION(kbdim,klev,nb_lw)  :: zext_t, zomg_t
  REAL(wp), DIMENSION(kbdim,nb_lw)       :: zaod_t, zext_t_int, zfact_t 
  INTEGER, DIMENSION(kbdim,klev)    :: kindex ! index field
  REAL(wp), PARAMETER               :: rdog=rd*rgrav

! 1. calculate for each echam gridbox the index of the data set layer 
!     in which p_mid_echam is located and geometrical height of layers
  CALL pressure_index(kproma,        kbdim,         klev,              &
                      pp_fl,         lev_clim,      p_lim_clim,        &
                      kindex)
  zdeltag(1:kproma,1:klev)= &
       & ppd_hl(1:kproma,1:klev)* &
       & tk_fl(1:kproma,1:klev)/pp_fl(1:kproma,1:klev)*rdog
! 2. Solar radiation
! 2.1 interpolate optical properties solar radiation
  DO jwl=1,nb_sw
     DO jk=1,klev
        DO jl=1,kproma
           idx_lat=ilat(jl,krow)
           idx_lev=kindex(jl,jk)
           zext_s(jl,jk,jwl)=wgt1_m*ext_v_s(jwl,idx_lev,idx_lat,nmw1_m)+ &
                             wgt2_m*ext_v_s(jwl,idx_lev,idx_lat,nmw2_m)
           zomg_s(jl,jk,jwl)=wgt1_m*ssa_v_s(jwl,idx_lev,idx_lat,nmw1_m)+ &
                             wgt2_m*ssa_v_s(jwl,idx_lev,idx_lat,nmw2_m)
           zasy_s(jl,jk,jwl)=wgt1_m*asy_v_s(jwl,idx_lev,idx_lat,nmw1_m)+ &
                             wgt2_m*asy_v_s(jwl,idx_lev,idx_lat,nmw2_m)
        END DO
     END DO
  END DO
  DO jwl=1,nb_sw
     DO jl=1,kproma
        idx_lat=ilat(jl,krow)
        zaod_s(jl,jwl)=wgt1_m*aod_v_s(jwl,idx_lat,nmw1_m)+ &
                       wgt2_m*aod_v_s(jwl,idx_lat,nmw2_m)
     END DO
  END DO
! 2.2 normalize zext to the correct total optical depth
!     the normalization factor generally depends on the wavelength if
!     the ratios of the extinction at different wavelengths are not 
!     independent of the height level. Generally, the aerosol composition
!     depends on height, this leads to different ratios of the extinction
!     between two given wavelengths at different heights.
  zext_s_int(1:kproma,1:nb_sw)=0._wp
  DO jwl=1,nb_sw
     DO jk=1,klev
        zext_s_int(1:kproma,jwl)=zext_s_int(1:kproma,jwl) + &
          zext_s(1:kproma,jk,jwl)*zdeltag(1:kproma,jk)
     END DO
  END DO
  WHERE (zext_s_int(1:kproma,1:nb_sw) > 0._wp) 
     zfact_s(1:kproma,1:nb_sw)=zaod_s(1:kproma,1:nb_sw)/ &
                              zext_s_int(1:kproma,1:nb_sw)
  ELSEWHERE
     zfact_s(1:kproma,1:nb_sw)=1._wp
  END WHERE
  DO jwl=1,nb_sw
     DO jk=1,klev
        zext_s(1:kproma,jk,jwl)=zext_s(1:kproma,jk,jwl)* &
             zdeltag(1:kproma,jk)*zfact_s(1:kproma,jwl)
     END DO
  END DO
! 2.3 add optical parameters to the optical parameters of aerosols
!     inverse height profile
  DO jk=1,klev
     jki=klev-jk+1
     WHERE (zext_s(1:kproma,jki,1:nb_sw)>0._wp) 
     paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*&
       paer_piz_sw_vr(1:kproma,jk,1:nb_sw)*paer_cg_sw_vr(1:kproma,jk,1:nb_sw)+&
       zext_s(1:kproma,jki,1:nb_sw)*zomg_s(1:kproma,jki,1:nb_sw)*&
       zasy_s(1:kproma,jki,1:nb_sw)
     paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*&
       paer_piz_sw_vr(1:kproma,jk,1:nb_sw)+&
       zext_s(1:kproma,jki,1:nb_sw)*zomg_s(1:kproma,jki,1:nb_sw)
     paer_tau_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)+&
       zext_s(1:kproma,jki,1:nb_sw)
     paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=paer_piz_sw_vr(1:kproma,jk,1:nb_sw)/&
          paer_tau_sw_vr(1:kproma,jk,1:nb_sw)
     paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=paer_cg_sw_vr(1:kproma,jk,1:nb_sw)/&
          (paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*paer_piz_sw_vr(1:kproma,jk,1:nb_sw))
     END WHERE
  END DO  
! 3. far infrared
! 2.1 interpolate optical properties solar radiation
  DO jwl=1,nb_lw
     DO jk=1,klev
        DO jl=1,kproma
           idx_lat=ilat(jl,krow)
           idx_lev=kindex(jl,jk)
           zext_t(jl,jk,jwl)=wgt1_m*ext_v_t(jwl,idx_lev,idx_lat,nmw1_m)+ &
                             wgt2_m*ext_v_t(jwl,idx_lev,idx_lat,nmw2_m)
           zomg_t(jl,jk,jwl)=wgt1_m*ssa_v_t(jwl,idx_lev,idx_lat,nmw1_m)+ &
                             wgt2_m*ssa_v_t(jwl,idx_lev,idx_lat,nmw2_m)
        END DO
     END DO
  END DO
  DO jwl=1,nb_lw
     DO jl=1,kproma
        idx_lat=ilat(jl,krow)
        zaod_t(jl,jwl)=wgt1_m*aod_v_t(jwl,idx_lat,nmw1_m)+ &
                       wgt2_m*aod_v_t(jwl,idx_lat,nmw2_m)
     END DO
  END DO
! 2.2 normalize zext to the correct total optical depth
!     the normalization factor generally depends on the wavelength if
!     the ratios of the extinction at different wavelengths are not 
!     independent of the height level. Generally, the aerosol composition
!     depends on height, this leads to different ratios of the extinction
!     between two given wavelengths at different heights.
  zext_t_int(1:kproma,1:nb_lw)=0._wp
  DO jwl=1,nb_lw
     DO jk=1,klev
        zext_t_int(1:kproma,jwl)=zext_t_int(1:kproma,jwl) + &
          zext_t(1:kproma,jk,jwl)*zdeltag(1:kproma,jk)
     END DO
  END DO
  WHERE (zext_t_int(1:kproma,1:nb_lw) > 0._wp) 
     zfact_t(1:kproma,1:nb_lw)=zaod_t(1:kproma,1:nb_lw)/ &
                              zext_t_int(1:kproma,1:nb_lw)
  ELSEWHERE
     zfact_t(1:kproma,1:nb_lw)=1._wp
  END WHERE
  DO jwl=1,nb_lw
     DO jk=1,klev
        zext_t(1:kproma,jk,jwl)=zext_t(1:kproma,jk,jwl)* &
             zdeltag(1:kproma,jk)*zfact_t(1:kproma,jwl)
     END DO
  END DO
! 2.3 add optical parameters to the optical parameters of aerosols
!     inverse height profile
  DO jk=1,klev
     jki=klev-jk+1
     paer_tau_lw_vr(1:kproma,jk,1:nb_lw)=paer_tau_lw_vr(1:kproma,jk,1:nb_lw)+ &
          zext_t(1:kproma,jki,1:nb_lw)*(1._wp-zomg_t(1:kproma,jki,1:nb_lw))
  END DO  
END SUBROUTINE add_aop_volc
!EOP
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aero_read_opt_volc_sw
!
! !SUBROUTINE INTERFACE:
SUBROUTINE aero_read_opt_volc_sw ( & 
  nlat,            nwl,                lev_clim,                         &
  aod,             caod,               ssa,               cssa,          &
  asy,             casy,               ext,               cext,          &
  p_lim,           cp_lim,             imnthb,            imnthe,        &
  cfname                                                                 )
!
! !DESCRIPTION:
! read optical aerosol parameters from file containg
! aod, ssa, asy, ext (altitude dependent extinction),
! lev_clim (number of levels), p_mid, p_lim
!
! !REVISION HISTORY:
! original source by J.S. Rast (2010-01-08)
!
! !USES:
  USE mo_mpi,                 ONLY: p_parallel_io, p_io, p_bcast
  USE mo_read_netcdf77,       ONLY: read_var_hs_nf77_3d, &
                                    read_var_hs_nf77_2d, &
                                    read_diml_nf77, &
                                    read_var_nf77_1d
                                    

! !INPUT PARAMETERS
  INTEGER, INTENT(in)            :: imnthb, imnthe !begin and end month 2b read
  INTEGER, INTENT(in)            :: nlat, nwl, lev_clim
  CHARACTER(len=*), INTENT(in)   :: cfname
  CHARACTER(len=*), INTENT(in)   :: caod, cssa, casy, cext, cp_lim
  REAL(wp), INTENT(out)          :: aod(nwl,nlat,imnthb:imnthe), &
                                    ext(nwl,lev_clim,nlat,imnthb:imnthe), &
                                    ssa(nwl,lev_clim,nlat,imnthb:imnthe), &
                                    asy(nwl,lev_clim,nlat,imnthb:imnthe)
  REAL(wp), INTENT(out)          :: p_lim(lev_clim+1)

! !LOCAL VARIABLES
  LOGICAL                        :: lex
  INTEGER                        :: j,ierr,inwl,ilev_clim

  IF (p_parallel_io) THEN
     INQUIRE (file=TRIM(cfname), exist=lex)
     IF (.NOT. lex) THEN
        CALL finish('aero_read_opt_volc_sw','file '//TRIM(cfname)// &
                    ' does not exist')
     END IF
     inwl=read_diml_nf77(TRIM(cfname),'nwl')
     IF (inwl /= nwl) THEN
        CALL finish('aero_read_opt_volc_sw', &
                   'incompatible number of wavelengths in file '//TRIM(cfname))
     END IF
     ilev_clim=read_diml_nf77(TRIM(cfname),'mlev')
     IF (ilev_clim /= lev_clim) THEN
        CALL finish('aero_read_opt_volc_sw', &
                    'incompatible number of levels in file '//TRIM(cfname))
     END IF
     CALL read_var_nf77_1d(TRIM(cfname), 'ilev', TRIM(cp_lim), p_lim, ierr)
  END IF
  CALL p_bcast(p_lim,p_io)
  DO j=imnthb,imnthe
     IF (p_parallel_io) THEN
        CALL read_var_hs_nf77_3d(TRIM(cfname),'nwl','mlev','lat','time',j, &
             TRIM(cext),ext(:,:,:,j),ierr)
        CALL read_var_hs_nf77_3d(TRIM(cfname),'nwl','mlev','lat','time',j, &
             TRIM(cssa),ssa(:,:,:,j),ierr)
        CALL read_var_hs_nf77_3d(TRIM(cfname),'nwl','mlev','lat','time',j, &
             TRIM(casy),asy(:,:,:,j),ierr)
        CALL read_var_hs_nf77_2d(TRIM(cfname),'nwl','lat','time',j, &
             TRIM(caod),aod(:,:,j),ierr)
     END IF
  ENDDO
  CALL p_bcast(ext(:,:,:,imnthb:imnthe), p_io)
  CALL p_bcast(ssa(:,:,:,imnthb:imnthe), p_io)
  CALL p_bcast(asy(:,:,:,imnthb:imnthe), p_io)
  CALL p_bcast(aod(:,:,imnthb:imnthe), p_io)
END SUBROUTINE aero_read_opt_volc_sw
!END SUBROUTINE 
!EOP
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aero_read_opt_volc_ir
!
! !SUBROUTINE INTERFACE:
SUBROUTINE aero_read_opt_volc_ir ( & 
  nlat,            nwl,                lev_clim,                         &
  aod,             caod,               ssa,               cssa,          &
                                       ext,               cext,          &
  p_lim,           cp_lim,             imnthb,            imnthe,        &
  cfname                                                                 )
!
! !DESCRIPTION:
! read optical aerosol parameters from file containg
! aod, ssa, ext (altitude dependent extinction),
! lev_clim (number of levels), p_mid, p_lim
!
! !REVISION HISTORY:
! original source by J.S. Rast (2010-02-23)
!
! !USES:
  USE mo_mpi,                 ONLY: p_parallel_io, p_io, p_bcast
  USE mo_read_netcdf77,       ONLY: read_var_hs_nf77_3d, &
                                    read_var_hs_nf77_2d, &
                                    read_diml_nf77, &
                                    read_var_nf77_1d
                                    

! !INPUT PARAMETERS
  INTEGER, INTENT(in)            :: imnthb, imnthe !begin and end month 2b read
  INTEGER, INTENT(in)            :: nlat, nwl, lev_clim
  CHARACTER(len=*), INTENT(in)   :: cfname
  CHARACTER(len=*), INTENT(in)   :: caod, cssa, cext, cp_lim
  REAL(wp), INTENT(out)          :: aod(nwl,nlat,imnthb:imnthe), &
                                    ext(nwl,lev_clim,nlat,imnthb:imnthe), &
                                    ssa(nwl,lev_clim,nlat,imnthb:imnthe)
  REAL(wp), INTENT(out)          :: p_lim(lev_clim+1)

! !LOCAL VARIABLES
  LOGICAL                        :: lex
  INTEGER                        :: j,ierr,inwl,ilev_clim

  IF (p_parallel_io) THEN
     INQUIRE (file=TRIM(cfname), exist=lex)
     IF (.NOT. lex) THEN
        CALL finish('aero_read_opt_volc_ir','file '//TRIM(cfname)// &
                    ' does not exist')
     END IF
     inwl=read_diml_nf77(TRIM(cfname),'nwl')
     IF (inwl /= nwl) THEN
        CALL finish('aero_read_opt_volc_ir', &
                   'incompatible number of wavelengths in file '//TRIM(cfname))
     END IF
     ilev_clim=read_diml_nf77(TRIM(cfname),'mlev')
     IF (ilev_clim /= lev_clim) THEN
        CALL finish('aero_read_opt_volc_ir', &
                    'incompatible number of levels in file '//TRIM(cfname))
     END IF
     CALL read_var_nf77_1d(TRIM(cfname), 'ilev', TRIM(cp_lim), p_lim, ierr)
  END IF
  CALL p_bcast(p_lim,p_io)
  DO j=imnthb,imnthe
     IF (p_parallel_io) THEN
        CALL read_var_hs_nf77_3d(TRIM(cfname),'nwl','mlev','lat','time',j, &
             TRIM(cext),ext(:,:,:,j),ierr)
        CALL read_var_hs_nf77_3d(TRIM(cfname),'nwl','mlev','lat','time',j, &
             TRIM(cssa),ssa(:,:,:,j),ierr)
        CALL read_var_hs_nf77_2d(TRIM(cfname),'nwl','lat','time',j, &
             TRIM(caod),aod(:,:,j),ierr)
     END IF
  ENDDO
  CALL p_bcast(ext(:,:,:,imnthb:imnthe), p_io)
  CALL p_bcast(ssa(:,:,:,imnthb:imnthe), p_io)
  CALL p_bcast(aod(:,:,imnthb:imnthe), p_io)
END SUBROUTINE aero_read_opt_volc_ir
!END SUBROUTINE 
!EOP
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pressure index
!
SUBROUTINE pressure_index(kproma,        kbdim,         klev,              &
                          pp_mid,        klevels,       pp_bound,          &
                          kindex)
!
! !DESCRIPTION:
! find index of pressure layer in which mid level pressures of echam
! are located
!
! !REVISION HISTORY:
! original source by J.S. Rast (2010-02-17)
!
! !USES:
  USE mo_kind,                ONLY: wp
  IMPLICIT NONE

! !INPUT PARAMETERS
  INTEGER, INTENT(in)    :: kbdim, kproma, klev
  INTEGER, INTENT(in)    :: klevels !number of layers for indices are searched
  REAL(wp), INTENT(in)   :: pp_mid(kbdim,klev), & !echam midlevel pressures
                            pp_bound(klevels+1) !pressure at layer 
                                    !bounds of reference pressures
  INTEGER, INTENT(out)   :: kindex(kbdim,klev) !layer indices for echam press.

! !LOCAL VARIABLES

  LOGICAL                :: lp(kproma), lrepeat
  INTEGER                :: jk,il,kidx(kbdim)
  
  kidx(1:kproma)=2
  DO jk=1,klev
10   CONTINUE
     lrepeat=.FALSE.
     DO il=1,kproma
        lp(il)=pp_mid(il,jk).GT.pp_bound(kidx(il)).AND.kidx(il).LE.klevels
     END DO
     DO il=1,kproma
        IF (lp(il)) THEN
           kidx(il)=kidx(il)+1
           lrepeat=.TRUE.
        END IF
     END DO
     IF (lrepeat) THEN
        GOTO 10
     ELSE
        kindex(1:kproma,jk)=kidx(1:kproma)-1
     END IF
  END DO
END SUBROUTINE pressure_index

!END SUBROUTINE 
!EOP
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cleanup_aero_volc
!
SUBROUTINE cleanup_aero_volc
!
! !DESCRIPTION: deallocate allocated fields
! 
!
! !REVISION HISTORY: original source J.S.Rast, MPI-Met, 2001-01-31
! 
!
! !USES:
! !LOCAL VARIABLES

  IF(ALLOCATED(aod_v_s)) DEALLOCATE(aod_v_s)
  IF(ALLOCATED(ext_v_s)) DEALLOCATE(ext_v_s)
  IF(ALLOCATED(ssa_v_s)) DEALLOCATE(ssa_v_s)
  IF(ALLOCATED(asy_v_s)) DEALLOCATE(asy_v_s)
  IF(ALLOCATED(aod_v_t)) DEALLOCATE(aod_v_t)
  IF(ALLOCATED(ext_v_t)) DEALLOCATE(ext_v_t)
  IF(ALLOCATED(ssa_v_t)) DEALLOCATE(ssa_v_t)
  IF(ALLOCATED(p_lim_clim)) DEALLOCATE(p_lim_clim)

  laero_set = .FALSE.

END SUBROUTINE cleanup_aero_volc
!EOP
!-------------------------------------------------------------------------
!EOC
!-------------------------------------------------------------------------
END MODULE mo_aero_volc
