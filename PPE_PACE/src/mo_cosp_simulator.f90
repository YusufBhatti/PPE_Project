 ! Purpose:
 ! --------
 ! Model output for comparisons with CloudSat, CALIPSO, MODIS, ISCCP, MISR and RTTOV satellite data
 !
 ! The COSP simulator is (c) copyright by the British Crown / Met Office 2008
 ! See http://cfmip.metoffice.com/cosp/cosp.v0.3/  and
 ! Refer to Met_Office_licencse.text for details
 !
 ! Version v0.1.beta CICCS implemented by Christine Nam, July    2008
 ! Version v0.3.beta COSP - Lidar implemented by C.Nam, October 2008
 ! Version v0.3.beta COSP - Radar implemented by C.Nam, January 2009
 ! Version v1.0 COSP - Lidar updated by C.Nam, April 2009
 ! Version v1.1 COSP - Lidar updated by C.Nam, May 2009
 ! Version v1.1 COSP - Radar updated by C.Nam, October 2009
 ! Version v1.2.1 COSP - Lidar & Radar update implemented by C.Nam, Febuary 2010
 !
 ! M.Salzmann (cms): Initial implementation into ECHAM6
 !	- Update Lidar and ISCCP simulator to v1.3
 !      - De-implemented radar, threw away cosp types & convective cloud 
 !	- Include time control, openMP compatible, some bug fixes.
 !
 ! Version v1.4 COSP - Version control implementation C.Nam, March 2014
 !		     - Include: Radar, Lidar, MODIS, ISCCP
 !                   - Uses 'TYPES' for version control
 !                   - MISR, RRTOV still not implemented yet. 
 !

!>>KH change name for backward compatibility with echam-ham with prev version of cosp
!MODULE mo_cosp_echam
MODULE mo_cosp_simulator
!<<KH
  USE mo_kind,         ONLY: wp
  USE mo_linked_list,  ONLY: t_stream, HYBRID, SURFACE, BELOWSUR, &
                              ABOVESUR10, NETCDF, GRIB
  USE mo_memory_base,  ONLY: new_stream, add_stream_element, &
                              default_stream_setting, add_stream_reference, AUTO
  USE mo_cosp_v1p4_cosp
  USE mo_cosp_types
  USE mo_tr_omp_decomposition, ONLY: omp_get_thread_num

  IMPLICIT NONE

  INTEGER, SAVE      :: cospo_index ! event index of the COSP stream

  TYPE, PUBLIC :: vmem3d
     REAL(wp), POINTER :: ptr(:,:,:)
  END TYPE vmem3d
  
  TYPE, PUBLIC :: vmem2d
     REAL(wp), POINTER :: ptr2(:,:)
  END TYPE vmem2d

  PUBLIC:: cosp_initialize
  PUBLIC:: construct_stream_cosp
!>>KH change name for compatibility with echam-ham with prev version of cosp
  !PUBLIC:: echam_cospsimulator
  PUBLIC:: call_cospsimulator
!<<KH

  LOGICAL :: locosp = .FALSE.     ! Doing satellite output?
  LOGICAL, SAVE :: Lradar_sim, LcfadDbze94, Llidar_sim, LcfadLidarsr532, Lisccp_sim, &
                   Lmodis_sim, Lmisr_sim, Lrttov_sim, Lstats
  LOGICAL :: l_fixed_reff ! use fixed effective radii for testing only
  LOGICAL :: extra_output
  LOGICAL :: use_netcdf
  LOGICAL :: use_vgrid
  LOGICAL :: csat_vgrid
  INTEGER :: offl2dout

  !----------------------------------------------------------------------------------
  ! Inputs: Following values from cosp_input_nl.txt
  !----------------------------------------------------------------------------------
  INTEGER :: Ncolumns          ! number of sub-columns used for each profile (namelist input)
  INTEGER, PARAMETER :: cosp_overlap = 3       ! Overlap Type: 1=max, 2=rand, 3=max/rand 
  INTEGER, PARAMETER :: Nlr = 40  ! Number of levels in statistical outputs (only used if USE_VGRID=.true.)

  !----------------------------------------------------------------------------------
  !Inputs related to radar simulations
  !----------------------------------------------------------------------------------
  REAL(wp),PARAMETER :: RADAR_FREQ=94.0_wp 	! CloudSat radar frequency (GHz)
  REAL(wp),PARAMETER :: k2=-1_wp           	! |K|^2, -1=use frequency dependent default
  INTEGER, PARAMETER :: SURFACE_RADAR=0 	! surface=1, spaceborne=0
  INTEGER, PARAMETER :: use_mie_tables=0	! use a precomputed lookup table? yes=1,no=0
  INTEGER, PARAMETER :: use_gas_abs=1   	! include gaseous absorption? yes=1,no=0
  INTEGER, PARAMETER :: do_ray=0        	! calculate/output Rayleigh refl=1, not=0
  INTEGER, PARAMETER :: melt_lay=0      	! melting layer model off=0, on=1
  LOGICAL, PARAMETER :: use_reff = .TRUE. 	! .true.: use effective radii from model  
                                          	! .false.: defaults
  LOGICAL, PARAMETER :: use_precipitation_fluxes=.true.  ! True if precipitation fluxes are input to the algorithm 
  !----------------------------------------------------------------------------------
  ! Inputs related to lidar simulations
  !----------------------------------------------------------------------------------
  INTEGER, PARAMETER :: Nprmts_max_hydro=12 ! Max number of parameters for hydrometeor size distributions
  INTEGER, PARAMETER :: Naero=1             ! Number of aerosol species (Not used)
  INTEGER, PARAMETER :: Nprmts_max_aero=1   ! Max number of parameters for aerosol size distributions (Not used)
  INTEGER, PARAMETER :: lidar_ice_type=0    ! Ice particle shape in lidar calculations (0=ice-spheres ; 1=ice-non-spherical)
                                            ! (0=ice-spheres ; 1=ice-non-spherical)

  !----------------------------------------------------------------------------------
  ! Inputs related to ISCCP simulator
  !----------------------------------------------------------------------------------
  INTEGER, PARAMETER :: isccp_top_height = 1 ! 1 = adjust cloud top pressur using IR brightness 
                                             ! temperature and visible optical depth
  INTEGER, PARAMETER :: isccp_top_height_direction = 2 ! 2= find the lowest pressure level with interpolated  
                                             ! temperature equal to the radiance determined cloud-top temperature
  REAL(wp), PARAMETER :: isccp_emsfc_lw = 0.996_wp  ! Surface emissivity at 10.5 micron (fraction)
                                           !cms ... pretty close, check sensitivity ...!
  !----------------------------------------------------------------------------------
  !Inputs related to RTTOV inputs  : CNam
  !----------------------------------------------------------------------------------
!  INTEGER, PARAMETER :: Platform=1    ! satellite platform
!  INTEGER, PARAMETER :: Satellite=15  ! satellite
!  INTEGER, PARAMETER :: Instrument=0  ! instrument
!  INTEGER, PARAMETER :: Nchannels=8   ! Number of channels to be computed
!  INTEGER, PARAMETER :: Channels=1,3,5,6,8,10,11,13        ! Channel numbers (please be sure that you supply Nchannels)
!  REAL(wp), PARAMETER ::  Surfem=0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp  ! Surface emissivity (please be sure that you supply Nchannels)
!  REAL(wp), PARAMETER ::  ZenAng=50.0_wp ! Satellite Zenith Angle
!  REAL(wp), PARAMETER ::  CO2=5.241e-04_wp ! Mixing ratios of trace gases
!  REAL(wp), PARAMETER ::  CH4=9.139e-07_wp
!  REAL(wp), PARAMETER ::  N2O=4.665e-07_wp
!  REAL(wp), PARAMETER ::  CO=2.098e-07_wp

  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------
  !3D
  TYPE (vmem3d), PRIVATE, POINTER :: cosp_cfadze(:)
                                             ! COSP CFADSR for the 15 dbze bins of radar reflectivity
  TYPE (vmem3d), PRIVATE, POINTER :: cosp_cfadsr(:)
                                             ! COSP CFADSR for the 15 srbins of lidar scattering ratio
  TYPE (vmem2d), PRIVATE, POINTER :: cosp_parasol_refl(:)   
                                             ! PARASOL reflectance for 5 bins of solar zenith angle
  TYPE (vmem2d), PRIVATE, POINTER :: isccp_cldtypes(:)   
                                             ! ISCCP cloud fraction for the 49 ISCCP cloud types
  TYPE (vmem2d), PRIVATE, POINTER :: modis_cldtypes(:)   
                                             ! MODIS cloud fraction for the 49 ISCCP cloud types
  TYPE (vmem2d), PRIVATE, POINTER :: misr_cldtypes(:)   
                                             ! MISR cloud fraction for the 49 ISCCP cloud types

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_f3d ! cms: 3-d cloud fraction as seen in radiation

  !radar
  REAL(wp), PRIVATE, POINTER :: cosp_radar_lidar_tcc(:,:)		! COSP LidarRadar
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_only_freq_cloud(:,:,:)  	! COSP Lidar_only_freq_cloud


  !lidar
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_cld(:,:,:)          ! COSP Cloud frequency of occurance as seen by CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_cldice(:,:,:)       ! COSP Cloud freq ice as seen by CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_cldliq(:,:,:)       ! COSP Cloud freq liq as seen by CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_cldun(:,:,:)        ! COSP Cloud freq undef as seen by CALIPSO

  REAL(wp), PRIVATE, POINTER :: cosp_lidar_cldtmp(:,:,:)       ! COSP Cloud temp as seen by CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_cldtmpice(:,:,:)    ! COSP Cloud temp ice as seen by CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_cldtmpliq(:,:,:)    ! COSP Cloud temp liq as seen by CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_cldtmpun(:,:,:)     ! COSP Cloud temp undef as seen by CALIPSO

  REAL(wp), PRIVATE, POINTER :: cosp_lidar_lowcloud(:,:)      ! COSP Low-level cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_midcloud(:,:)      ! COSP Mid-level cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_highcloud(:,:)     ! COSP High-level cloud cover from CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_totalcloud(:,:)    ! COSP Total cloud cover from CALIPSO

  REAL(wp), PRIVATE, POINTER :: cosp_lidar_lowcloudice(:,:)   ! COSP Low-level ice cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_midcloudice(:,:)   ! COSP Mid-level ice cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_highcloudice(:,:)  ! COSP High-level ice cloud cover from CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_totalcloudice(:,:) ! COSP Total ice cloud cover from CALIPSO

  REAL(wp), PRIVATE, POINTER :: cosp_lidar_lowcloudliq(:,:)   ! COSP Low-level liq cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_midcloudliq(:,:)   ! COSP Mid-level liq cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_highcloudliq(:,:)  ! COSP High-level liq cloud cover from CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_totalcloudliq(:,:) ! COSP Total liq cloud cover from CALIPSO

  REAL(wp), PRIVATE, POINTER :: cosp_lidar_lowcloudun(:,:)    ! COSP Low-level undef cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_midcloudun(:,:)    ! COSP Mid-level undef cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_highcloudun(:,:)   ! COSP High-level undef cloud cover from CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_totalcloudun(:,:)  ! COSP Total undef cloud cover from CALIPSO


  !isccp
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_isccp_totalcloud    ! ISCCP total cloud cover 
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_isccp_meanptop      ! mean ISCCP cloud top pressure
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_isccp_meanalbedocld ! mean ISCCP cloud albedo
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_isccp_meantaucld    ! mean ISCCP cloud optical thickness

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cisccp_cldtau3d  ! Cloud optical thickness
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cisccp_cldemi3d  ! Cloud emissivity @ 10.5 µm


  !modis
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_clt   ! modis total cloud cover 
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_clw   ! modis liq cloud cover 
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_cli   ! modis ice cloud cover 
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_clh   ! modis high cloud cover 
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_clm   ! modis mid cloud cover 
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_cll   ! modis low cloud cover 

  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_taut  ! mean modis total optical thickness
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_taul  ! mean modis liq optical thickness
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_taui  ! mean modis ice optical thickness

  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_tautlog  ! Log mean modis total optical thickness
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_taullog  ! Log mean modis liq optical thickness
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_tauilog  ! Log mean modis ice optical thickness

  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_reffclw  ! mean modis liquid particle size
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_reffcli  ! mean modis ice particle size

  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_pct      ! mean modis cloud top pressure
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_lwp      ! mean modis liquid water path
  REAL(wp), PRIVATE, POINTER, DIMENSION(:,:) :: cosp_modis_iwp      ! mean modis ice water path

  !hydro
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_lsrain	! Flux large scale cloud rain [kg m^-2 s^-1]
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_lssnow	! Flux large scale cloud snow [kg m^-2 s^-1]
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_lsgrpl	! Flux large scale graupel [kg m^-2 s^-1]
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_ccrain	! Flux convective cloud rain [kg m^-2 s^-1]
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_ccsnow	! Flux convective cloud snow [kg m^-2 s^-1]

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_reffl     ! Liquid water droplet effective radius [um]
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_reffi     ! Ice crystal effective radius [um]

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: xi_cosp     
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: xl_cosp     
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: tm1_cosp     

  !sunlit
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:) :: cosp_sunlit    ! sunlit grid points
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:) :: cosp_sunlit_av ! average sunlit fraction (needed for correctly 
                                                              !  averaging output)
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:) :: cosp_freq ! frequency of calls 

  !N_miss
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_radar   ! radar sim: number of missing points

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_lidar      ! lidar sim: number of missing points
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_lidarice   ! lidar sim: n. missing points phase
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_lidarliq   ! lidar sim: n. missing points phase
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_lidarun    ! lidar sim: n. missing points phase
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_lidartmp      ! lidar sim: n. missing points tmp
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_lidartmpice   ! lidar sim: n. missing points tmp
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_lidartmpliq   ! lidar sim: n. missing points tmp
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_lidartmpun    ! lidar sim: n. missing points tmp

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_lidar_mid  ! lidar sim: number of missing points midcld
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_lidar_midice ! lidar sim: n. missing points ice 
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_lidar_midliq ! lidar sim: n. missing points liq
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_lidar_midun  ! lidar sim: n. missing points undef
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_lidar_low ! lidar sim: n. missing points lowcld
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_lidar_lowice ! lidar sim: n. missing points ice
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_lidar_lowliq ! lidar sim: n. missing points liq
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_lidar_lowun  ! lidar sim: n. missing points undef

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: N_miss_modis   ! modis sim: number of missing points

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_modis_clt ! modis sim: n. missing points totalcld
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_modis_clw ! modis sim: n. missing points liqcld
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_modis_cli ! modis sim: n. missing points icecld
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_modis_clh ! modis sim: n. missing points highcld
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_modis_clm ! modis sim: n. missing points midcld
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:)   :: N_miss_modis_cll ! modis sim: n. missing points lowcld


  REAL(wp), SAVE :: od_cosp ! fraction of day calls to satellite simulator 
                            ! (currently equal to radiation time step). 

  TYPE(cosp_config)  :: cfg             ! Configuration options
  TYPE(cosp_gridbox) :: gbx             ! Gridbox information. Input for COSP
  TYPE(cosp_subgrid) :: sgx             ! Subgrid outputs
  TYPE(cosp_sgradar) :: sgradar         ! Output from radar simulator
  TYPE(cosp_sglidar) :: sglidar         ! Output from lidar simulator
  TYPE(cosp_isccpstats)   :: isccp   	! Output from ISCCP simulator !CNam: renamed TYPE
  TYPE(cosp_modisstats)   :: modis   	! Output from MODIS simulator !CNam: renamed TYPE
  TYPE(cosp_misrstats)    :: misr    	! Output from MISR simulator  !CNam: renamed TYPE
!CNam:  TYPE(cosp_rttovstats):: rttov   	! Output from RTTOV   !CNam: renamed TYPE
  TYPE(cosp_vgrid)   :: vgrid   	! Information on vertical grid of stats
  TYPE(cosp_radarstats) :: stradar      ! Summary statistics from radar simulator
  TYPE(cosp_lidarstats) :: stlidar      ! Summary statistics from lidar simulator

  TYPE (t_stream), PUBLIC, POINTER :: scosp


CONTAINS
  
!--------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------

  SUBROUTINE cosp_initialize

   USE mo_mpi,         ONLY: p_parallel, p_parallel_io, p_bcast, p_io
   USE mo_namelist,    ONLY: open_nml, position_nml, POSITIONED
   USE mo_time_control,ONLY: trigrad
   USE mo_exception,   ONLY: finish
   USE mo_kind,        ONLY: wp

   IMPLICIT NONE 

   INTEGER :: ierr, inml, iunit

   NAMELIST /cospctl/ locosp, Lradar_sim, LcfadDbze94, Llidar_sim, LcfadLidarsr532, Lisccp_sim, &
                       Lmodis_sim, Lmisr_sim, Lrttov_sim, &
                       l_fixed_reff, Ncolumns, extra_output, use_netcdf, use_vgrid, csat_vgrid, &
                       offl2dout


    Lradar_sim  = .false.
    LcfadDbze94 = .false.
    Llidar_sim  = .false.
    LcfadLidarsr532 = .false.
    Lisccp_sim  = .false.
    Lmodis_sim  = .false.
    Lmisr_sim   = .false.
    Lrttov_sim  = .false.

    if (Lmodis_sim) Lisccp_sim = .true.

    l_fixed_reff = .false.
    extra_output = .false.
    use_netcdf  = .true.
    use_vgrid = .false.
    csat_vgrid = .false. ! CloudSat vertical grid (if .true. then USE_VGRID needs be .true.)
    Ncolumns = 12 !CNam: 25 CMS:12
    offl2dout=-1
    
      !--- Local variables:
     IF (p_parallel_io) THEN
       inml = open_nml ('namelist.echam')
       iunit = position_nml ('COSPCTL', inml, status=ierr)
       SELECT CASE (ierr)
       CASE (POSITIONED)
         READ (iunit, cospctl)          
       END SELECT
     END IF

    IF (offl2dout .GT. 0) THEN
      extra_output = .true.
    END IF

      !--- 2) Broadcast over processors:
     IF (p_parallel) THEN
         CALL p_bcast (locosp, p_io)
         CALL p_bcast (Lradar_sim, p_io)
         CALL p_bcast (LcfadDbze94, p_io)
         CALL p_bcast (Llidar_sim, p_io)
         CALL p_bcast (LcfadLidarsr532, p_io)
         CALL p_bcast (Lisccp_sim, p_io)
         CALL p_bcast (Lmisr_sim, p_io)
         CALL p_bcast (Lmodis_sim, p_io)
         CALL p_bcast (Lrttov_sim, p_io)
         CALL p_bcast (l_fixed_reff, p_io)
         CALL p_bcast (use_vgrid, p_io)
         CALL p_bcast (csat_vgrid, p_io)
         CALL p_bcast (Ncolumns, p_io)
         CALL p_bcast (extra_output, p_io)
         CALL p_bcast (offl2dout, p_io)
      END IF ! p_parallel

   cfg%Lstats = .false.
   IF ((Lradar_sim).or.(Llidar_sim).OR.(Lisccp_sim)) cfg%Lstats = .TRUE.
   
   ! Copy final flags to cfg structure
   cfg%Lradar_sim = Lradar_sim
       cfg%LcfadDbze94 = LcfadDbze94
   cfg%Llidar_sim = Llidar_sim
       cfg%LcfadLidarsr532 = LcfadLidarsr532  
   cfg%Lisccp_sim = Lisccp_sim
   cfg%Lmodis_sim = Lmodis_sim
   cfg%Lmisr_sim  = Lmisr_sim
!   cfg%Lrttov_sim = Lrttov_sim 

   vgrid%use_vgrid = use_vgrid
   vgrid%csat_vgrid = csat_vgrid


   IF (.NOT. locosp) RETURN
    ! timestep for satellite simulator call
   IF ( trigrad%unit .EQ. 'hours' ) THEN
      od_cosp =  3600._wp * REAL(trigrad%counter, wp) / 86400._wp
   ELSE
!>>KH changed name of module for compatibility with echam-ham
     !CALL finish('mo_cosp_echam ','expecting different trigrad unit')
     CALL finish('mo_cosp_simulator ','expecting different trigrad unit')
!<<KH
   END IF
 
END SUBROUTINE cosp_initialize

!-----------------------------------------------------------------------------------------------------------------

  SUBROUTINE construct_stream_cosp

    USE mo_time_event, ONLY: io_time_event
    USE mo_exception,  ONLY: finish

    IMPLICIT NONE

    INTEGER :: i
    CHARACTER*40 :: name
    !
    ! Local
    ! 

!--- 1) Construct the cosp stream: ------------------------------------------------------------------------------

! do not change output frequency
    IF (use_netcdf) THEN
     IF (offl2dout .LT. 0 ) THEN !default
      CALL new_stream (scosp,'cosp',filetype=NETCDF, interval=io_time_event(1,'days','last',0) )
     ELSE
      CALL new_stream (scosp,'cosp',filetype=NETCDF, interval=io_time_event(offl2dout,'minutes','last',0) )
     END IF 
    ELSE
!>>KH changed name of module for compatibility with echam-ham
      !CALL finish('mo_cosp_echam ','GRIB OUTPUT NOT YET IMPLEMENTED')
      CALL finish('mo_cosp_simulator ','GRIB OUTPUT NOT YET IMPLEMENTED')
!<<KH
      CALL new_stream (scosp,'cosp',filetype=GRIB, interval=io_time_event(1,'days','last',0) )
    END IF
    cospo_index  = scosp%post_idx

    !--- Add standard fields for post-processing:

    CALL add_stream_reference (scosp, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (scosp, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (scosp, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (scosp, 'gboxarea','geoloc',lpost=.TRUE.)

    IF (offl2dout .GT. 0 ) THEN
        CALL add_stream_reference (scosp, 'u10'   ,'g3b'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'v10'   ,'g3b'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'slm'   ,'g3b'   ,lpost=.TRUE.)

        CALL add_stream_reference (scosp, 'aclc'   ,'g3b'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'ao3'    ,'g3b'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'xl'     ,'gl'    ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'xi'     ,'gl'    ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'q'      ,'gl'    ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'tm1'    ,'g1a'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'relhum' ,'g3b'   ,lpost=.TRUE.)

    CALL add_stream_element (scosp, 'xi_cosp', xi_cosp,  laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=.TRUE., &
         longname='xi_cosp', units='')

    CALL add_stream_element (scosp, 'xl_cosp', xl_cosp,  laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=.TRUE., &
         longname='xl_cosp', units='')

    CALL add_stream_element (scosp, 'tm1_cosp', tm1_cosp,  laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=.TRUE., &
         longname='tm1_cosp', units='')


    END IF

!--- 2) Add stream elements: ------------------------------------------------------------------------------------

    CALL default_stream_setting (scosp, units     = '',          &
                                         lrerun    = .FALSE. ,    &
                                         lpost     = .TRUE. ,    &
                                         laccu     = .FALSE. ,   &
                                         reset     = 1.e-50_wp,  &
                                         leveltype = SURFACE ,   &
                                         contnorest = .TRUE. )

! lpost set to TRUE for debugging 
   CALL add_stream_element (scosp, 'lsrain', cosp_lsrain, &
                             leveltype = HYBRID, lpost=.TRUE., &
                             longname='Large-scale rain', units='kg kg-1')
    CALL add_stream_element (scosp, 'lssnow', cosp_lssnow, &
                             leveltype = HYBRID, lpost=.TRUE., &
                             longname='Large-scale snow', units='kg kg-1')
    CALL add_stream_element (scosp, 'ccrain', cosp_ccrain, &
                             leveltype = HYBRID, lpost=.TRUE., &
                             longname='Convective rain', units='kg kg-1')
    CALL add_stream_element (scosp, 'ccsnow', cosp_ccsnow, &
                             leveltype = HYBRID, lpost=.TRUE., &
                             longname='Convective snow', units='kg kg-1')

    CALL add_stream_element (scosp, 'reffl', cosp_reffl,  laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=extra_output, &
         longname='Liquid water droplet effective radius', units='um')
    CALL add_stream_element (scosp, 'reffi', cosp_reffi,  laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=extra_output, &
         longname='Ice crystal effective radius', units='um')

    CALL add_stream_element (scosp, 'cosp_sunlit', cosp_sunlit, &
                              longname='COSP sunlit fraction')
    CALL add_stream_element (scosp, 'cosp_sunlit_av', cosp_sunlit_av, &
                              longname='COSP average sunlit fraction')
    CALL add_stream_element (scosp, 'cosp_freq', cosp_freq, &
                              longname='COSP frequency of calls')
    CALL add_stream_element (scosp, 'cosp_f3d', cosp_f3d, &
                              leveltype = HYBRID, &
                              longname='COSP cloud fraction',lpost=extra_output)

	!------------------------------------------------------------------------------
     IF (cfg%Lradar_sim) THEN
        CALL add_stream_element (scosp, 'N_miss_radar', N_miss_radar, &
                            leveltype = HYBRID, &
                            longname='COSP radar number of missing values')

        CALL add_stream_element (scosp, 'radar_lidar_tcc', cosp_radar_lidar_tcc, &
                             leveltype = HYBRID,&
                             longname='calipso_and_cloudsat_total_cloud_fraction', units='fraction')
	
	! 2d mean quantities
	CALL add_stream_element (scosp, 'lidar_only_freq_cloud', cosp_lidar_only_freq_cloud, &
     	 		     leveltype = HYBRID,&
                             longname='Cld Freq. seen by CALIPSO not CloudSat', units='fraction')

        IF (cfg%LcfadDbze94) THEN
	ALLOCATE(cosp_cfadze(15))
	    DO i=1,15
                 WRITE(name, '(a7,i2.2)') 'cfad_ze',i
	       CALL add_stream_element (scosp, name, cosp_cfadze(i)%ptr, &
	            leveltype = HYBRID, longname='COSP CFAD of Radar Reflect. 94GHz')
	    ENDDO
        ENDIF !cfg%LcfadDbze94

     ENDIF !cfg%Lradar_sim

	!------------------------------------------------------------------------------

     IF (cfg%Llidar_sim) THEN
        ! CNAM: Exclude 1D stream: srbval

	!3D (cloud) fields
        CALL add_stream_element (scosp, 'clcalipso', cosp_lidar_cld, &
                             leveltype = HYBRID,&
                             longname='3D Lidar Cloud Fraction (532nm)', units='fraction')
        CALL add_stream_element (scosp, 'clcalipsoice', cosp_lidar_cldice, &
                             leveltype = HYBRID,&
                             longname='3D Lidar Cloud Ice', units='fraction')
        CALL add_stream_element (scosp, 'clcalipsoliq', cosp_lidar_cldliq, &
                             leveltype = HYBRID,&
                             longname='3D Lidar Cloud Liquid', units='fraction')
        CALL add_stream_element (scosp, 'clcalipsoun', cosp_lidar_cldun, &
                             leveltype = HYBRID,&
                             longname='3D Lidar Cloud Undef', units='fraction')

      IF (vgrid%use_vgrid) THEN
        CALL add_stream_element (scosp, 'clcalipsotmp', cosp_lidar_cldtmp, &
                             leveltype = HYBRID,&
                             longname='3D Lidar CloudTemp', units='fraction')
        CALL add_stream_element (scosp, 'clcalipsotmpice', cosp_lidar_cldtmpice, &
                             leveltype = HYBRID,&
                             longname='3D Lidar CloudTemp Ice', units='fraction')
        CALL add_stream_element (scosp, 'clcalipsotmpliq', cosp_lidar_cldtmpliq, &
                             leveltype = HYBRID,&
                             longname='3D Lidar CloudTemp Liquid', units='fraction')
        CALL add_stream_element (scosp, 'clcalipsotmpun', cosp_lidar_cldtmpun, &
                             leveltype = HYBRID,&
                             longname='3D Lidar CloudTemp Undef', units='fraction')
      ENDIF

	!high, mid, low-cloud cover
        CALL add_stream_element (scosp, 'cllcalipso', cosp_lidar_lowcloud, &
                             longname='CALIPSO Low-level cloud fraction', units='fraction')
        CALL add_stream_element (scosp, 'clmcalipso', cosp_lidar_midcloud, &
                             longname='CALIPSO Mid-level cloud fraction', units='fraction')
        CALL add_stream_element (scosp, 'clhcalipso', cosp_lidar_highcloud, &
                             longname='CALIPSO High-level cloud fraction', units='fraction')
        CALL add_stream_element (scosp, 'cltcalipso', cosp_lidar_totalcloud, &
                             longname='CALIPSO Total cloud fraction', units='fraction')
	!HML ice clouds
        CALL add_stream_element (scosp, 'cllcalipsoice', cosp_lidar_lowcloudice, &
                             longname='CALIPSO Low-level cloud ice fraction', units='fraction')
        CALL add_stream_element (scosp, 'clmcalipsoice', cosp_lidar_midcloudice, &
                             longname='CALIPSO Mid-level cloud ice fraction', units='fraction')
        CALL add_stream_element (scosp, 'clhcalipsoice', cosp_lidar_highcloudice, &
                             longname='CALIPSO High-level cloud ice fraction', units='fraction')
        CALL add_stream_element (scosp, 'cltcalipsoice', cosp_lidar_totalcloudice, &
                             longname='CALIPSO Total cloud ice fraction', units='fraction')
	!HML liquid clouds
        CALL add_stream_element (scosp, 'cllcalipsoliq', cosp_lidar_lowcloudliq, &
                             longname='CALIPSO Low-level cloud liq fraction', units='fraction')
        CALL add_stream_element (scosp, 'clmcalipsoliq', cosp_lidar_midcloudliq, &
                             longname='CALIPSO Mid-level cloud liq fraction', units='fraction')
        CALL add_stream_element (scosp, 'clhcalipsoliq', cosp_lidar_highcloudliq, &
                             longname='CALIPSO High-level cloud liq fraction', units='fraction')
        CALL add_stream_element (scosp, 'cltcalipsoliq', cosp_lidar_totalcloudliq, &
                             longname='CALIPSO Total cloud liq fraction', units='fraction')

	!HML undef clouds
        CALL add_stream_element (scosp, 'cllcalipsoun', cosp_lidar_lowcloudun, &
                             longname='CALIPSO Low-level cloud un fraction', units='fraction')
        CALL add_stream_element (scosp, 'clmcalipsoun', cosp_lidar_midcloudun, &
                             longname='CALIPSO Mid-level cloud un fraction', units='fraction')
        CALL add_stream_element (scosp, 'clhcalipsoun', cosp_lidar_highcloudun, &
                             longname='CALIPSO High-level cloud un fraction', units='fraction')
        CALL add_stream_element (scosp, 'cltcalipsoun', cosp_lidar_totalcloudun, &
                             longname='CALIPSO Total cloud un fraction', units='fraction')

        IF (cfg%LcfadLidarsr532) THEN
           ALLOCATE(cosp_cfadsr(15))
              DO i=1,15
                WRITE(name, '(a7,i2.2)') 'cfad_sr',i
                CALL add_stream_element (scosp, name, cosp_cfadsr(i)%ptr, &
                  leveltype = HYBRID, longname='COSP CFAD Lidar Scattering Ratio 532nm')
              ENDDO
        END IF !cfg%LcfadLidarsr532

        ALLOCATE(cosp_parasol_refl(PARASOL_NREFL))
        DO i=1, PARASOL_NREFL
           WRITE(name, '(a11,i1)') 'parasolRefl',i
           CALL add_stream_element (scosp, name, cosp_parasol_refl(i)%ptr2, &
                            longname='Parasol like mono-directional reflectance', units='fraction')
        ENDDO

        CALL add_stream_element (scosp, 'N_miss_lidar_mid', N_miss_lidar_mid, &
                            longname='COSP lidar number of missing values mid level cloud')
        CALL add_stream_element (scosp, 'N_miss_lidar_midice', N_miss_lidar_midice, &
                            longname='COSP lidar number of missing values mid ice level cloud')
        CALL add_stream_element (scosp, 'N_miss_lidar_midliq', N_miss_lidar_midliq, &
                            longname='COSP lidar number of missing values mid liq level cloud')
        CALL add_stream_element (scosp, 'N_miss_lidar_midun', N_miss_lidar_midun, &
                            longname='COSP lidar number of missing values mid undef level cloud')

        CALL add_stream_element (scosp, 'N_miss_lidar_low', N_miss_lidar_low, &
                            longname='COSP lidar number of missing values low level cloud')
        CALL add_stream_element (scosp, 'N_miss_lidar_lowice', N_miss_lidar_lowice, &
                            longname='COSP lidar number of missing values low ice level cloud')
        CALL add_stream_element (scosp, 'N_miss_lidar_lowliq', N_miss_lidar_lowliq, &
                            longname='COSP lidar number of missing values low liq level cloud')
        CALL add_stream_element (scosp, 'N_miss_lidar_lowun', N_miss_lidar_lowun, &
                            longname='COSP lidar number of missing values low undef level cloud')

        CALL add_stream_element (scosp, 'N_miss_lidar', N_miss_lidar, &
                            leveltype = HYBRID, &
                            longname='COSP lidar number of missing values')
        CALL add_stream_element (scosp, 'N_miss_lidarice', N_miss_lidarice, &
                            leveltype = HYBRID, &
                            longname='COSP lidar ice number of missing values')
        CALL add_stream_element (scosp, 'N_miss_lidarliq', N_miss_lidarliq, &
                            leveltype = HYBRID, &
                            longname='COSP lidar liq number of missing values')
        CALL add_stream_element (scosp, 'N_miss_lidarun', N_miss_lidarun, &
                            leveltype = HYBRID, &
                            longname='COSP lidar undef number of missing values')

      IF (vgrid%use_vgrid) THEN
        CALL add_stream_element (scosp, 'N_miss_lidartmp', N_miss_lidartmp, &
                            leveltype = HYBRID, &
                            longname='COSP lidar n. missing values')
        CALL add_stream_element (scosp, 'N_miss_lidartmpice', N_miss_lidartmpice, &
                            leveltype = HYBRID, &
                            longname='COSP lidar ice n. missing values')
        CALL add_stream_element (scosp, 'N_miss_lidartmpliq', N_miss_lidartmpliq, &
                            leveltype = HYBRID, &
                            longname='COSP lidar temp liq n. missing values')
        CALL add_stream_element (scosp, 'N_miss_lidartmpun', N_miss_lidartmpun, &
                            leveltype = HYBRID, &
                            longname='COSP lidar temp undef n. missing values')
     ENDIF	

     END IF !cfg%Llidar_sim

	!------------------------------------------------------------------------------
        ! Following cms for ECHAM6

     IF (cfg%Lisccp_sim) THEN

        CALL add_stream_element (scosp, 'cltisccp', cosp_isccp_totalcloud, &
                           longname='ISCCP total cloud area fraction', units='fraction')
        CALL add_stream_element (scosp, 'pctisccp', cosp_isccp_meanptop, &
                           longname='ISCCP cloud top pressure', units='Pa')
        CALL add_stream_element (scosp, 'albisccp', cosp_isccp_meanalbedocld, &
                           longname='ISCCP mean cloud albedo', units='fraction')
        CALL add_stream_element (scosp, 'tauisccp', cosp_isccp_meantaucld, &
                           longname='ISCCP mean cloud tau')
        CALL add_stream_element (scosp, 'cisccp_tau3d', cisccp_cldtau3d, &
                           leveltype = HYBRID, lpost=extra_output, &
                           longname='COSP-ISCCP cloud optical thickness')
        CALL add_stream_element (scosp, 'cisccp_emi3d', cisccp_cldemi3d, &
                           leveltype = HYBRID, lpost=extra_output, &
                           longname='COSP-ISCCP cloud optical thickness')

        ALLOCATE(isccp_cldtypes(49))
        DO i=1,49
           WRITE(name, '(a7,i2.2)') 'isccp_cldtype',i
           CALL add_stream_element (scosp, name, isccp_cldtypes(i)%ptr2, &
                longname='Fraction covered by ISCCP cloud type')
        ENDDO

     END IF!cfg%Lisccp_sim

	!------------------------------------------------------------------------------
     IF (cfg%Lmodis_sim) THEN

      CALL add_stream_element (scosp, 'cltmodis', cosp_modis_clt, &
                           longname='MODIS total cloud area fraction', units='fraction')
      CALL add_stream_element (scosp, 'clwmodis', cosp_modis_clw, &
                           longname='MODIS liquid cloud area fraction', units='fraction')
      CALL add_stream_element (scosp, 'climodis', cosp_modis_cli, &
                           longname='MODIS ice cloud area fraction', units='fraction')

      CALL add_stream_element (scosp, 'clhmodis', cosp_modis_clh, &
                           longname='MODIS high cloud area fraction', units='fraction')
      CALL add_stream_element (scosp, 'clmmodis', cosp_modis_clm, &
                           longname='MODIS mid cloud area fraction', units='fraction')
      CALL add_stream_element (scosp, 'cllmodis', cosp_modis_cll, &
                           longname='MODIS low cloud area fraction', units='fraction')

      CALL add_stream_element (scosp, 'tautmodis', cosp_modis_taut, &
                           longname='MODIS Optical_Thickness_Total_Mean', units='1')
      CALL add_stream_element (scosp, 'tauwmodis', cosp_modis_taul, &
                           longname='MODIS Optical_Thickness_Liq_Mean', units='1')
      CALL add_stream_element (scosp, 'tauimodis', cosp_modis_taui, &
                           longname='MODIS Optical_Thickness_Ice_Mean', units='1')

      CALL add_stream_element (scosp, 'tautlogmodis', cosp_modis_tautlog, &
                           longname='MODIS Optical_Thickness_Total_LogMean', units='1')
      CALL add_stream_element (scosp, 'tauwlogmodis', cosp_modis_taullog, &
                           longname='MODIS Optical_Thickness_Liq_LogMean', units='1')
      CALL add_stream_element (scosp, 'tauilogmodis', cosp_modis_tauilog, &
                           longname='MODIS Optical_Thickness_Ice_LogMean', units='1')

      CALL add_stream_element (scosp, 'reffclwmodis', cosp_modis_reffclw, &
                           longname='MODIS Cloud_Particle_Size_Liq_Mean', units='m')
      CALL add_stream_element (scosp, 'reffclimodis', cosp_modis_reffcli, &
                           longname='MODIS Cloud_Particle_Size_Ice_Mean', units='m')

      CALL add_stream_element (scosp, 'pctmodis', cosp_modis_pct, &
                           longname='MODIS Cloud_Top_Pressure_Total_Mean', units='Pa')
      CALL add_stream_element (scosp, 'lwpmodis', cosp_modis_lwp, &
                           longname='MODIS Liquid_Water_Path_Mean', units='kg/m^2')
      CALL add_stream_element (scosp, 'iwpmodis', cosp_modis_iwp, &
                           longname='MODIS Ice_Water_Path_Mean', units='kg/m^2')


      CALL add_stream_element (scosp, 'N_miss_modis_total', N_miss_modis_clt, &
                            longname='COSP modis number of missing values total level cloud')
      CALL add_stream_element (scosp, 'N_miss_modis_liq', N_miss_modis_clw, &
                            longname='COSP modis number of missing values liq level cloud')
      CALL add_stream_element (scosp, 'N_miss_modis_ice', N_miss_modis_cli, &
                            longname='COSP modis number of missing values ice level cloud')
      CALL add_stream_element (scosp, 'N_miss_modis_high', N_miss_modis_clh, &
                            longname='COSP modis number of missing values high level cloud')
      CALL add_stream_element (scosp, 'N_miss_modis_mid', N_miss_modis_clm, &
                            longname='COSP modis number of missing values mid level cloud')
      CALL add_stream_element (scosp, 'N_miss_modis_low', N_miss_modis_cll, &
                            longname='COSP modis number of missing values low level cloud')

        ALLOCATE(modis_cldtypes(49))
        DO i=1,49
           WRITE(name, '(a7,i2.2)') 'modis_cldtype',i
           CALL add_stream_element (scosp, name, modis_cldtypes(i)%ptr2, &
                longname='Fraction covered by MODIS cloud type')
        ENDDO

     END IF!cfg%Lmodis_sim

	!------------------------------------------------------------------------------

     IF (cfg%Lmisr_sim) THEN

        print*,'----- Needs to be debugged -----' !CNam

!       DO i=1,49
!           WRITE(name, '(a7,i2.2)') 'misr_cldtype',i
!           CALL add_stream_element (scosp, name, misr_cldtypes(i)%ptr2, &
!                longname='Fraction covered by MISR cloud type')
!        ENDDO

     END IF!cfg%Lmisr_sim

	!------------------------------------------------------------------------------
  END SUBROUTINE construct_stream_cosp

!--------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
!>>KH change name for compatibility with echam-ham with prev version of cosp
  !SUBROUTINE echam_cospsimulator(             &
  SUBROUTINE call_cospsimulator(             &
!<<KH
               kproma,        klev,                   jrow,               &
               p,             ph,                     pgeo,               &
               pgeospm,       pslm,                   ptm1,               &
               qm,            rh,                     xlm1,               &
               xim1,          zfrl,                   zfrw,               &
               zfri,          tslm1,                  tsw,                &
               tsi                                                           )

    USE mo_kind,             ONLY: wp
    USE mo_physical_constants,ONLY: grav        ! Grativation acceleration
    USE mo_cosp_constants,   ONLY: R_UNDEF
    USE mo_cosp_types	      
    USE mo_time_control,     ONLY: l_trigrad, l_putdata
    USE mo_exception,        ONLY: message, message_text
    USE mo_random_numbers,   ONLY: set_seed_random
    USE mo_cloudtypes,       ONLY: nclusters, cslabel, isccp_ct, isccp_freq, isccp_sunlit
    USE mo_histogram,        ONLY: cosp_ct_lf, cosp_ct_cc, cosp_ct_lf_time, cosp_ct_cc_time
    IMPLICIT NONE

    !
    ! Input to COSP-Simulator set here
    !
    INTEGER                              :: kproma, klev, jrow ! Dimensions
    REAL(wp), DIMENSION(kproma, klev)    :: height       ! Height of model levels [m]
    REAL(wp), DIMENSION(kproma, klev+1)  :: height_half  ! Height @ layer interfaces  [m]

    REAL(wp), DIMENSION(kproma, klev)    :: p            ! Pressure @ full levels [Pa]
    REAL(wp), DIMENSION(kproma, klev+1)  :: ph           ! Pressure @ layer interfaces [Pa]
    REAL(wp), DIMENSION(kproma, klev)    :: pgeo         ! geopotential height above surface
    REAL(wp), DIMENSION(kproma)  :: pgeospm ! surface geopotential height
    REAL(wp), DIMENSION(kproma)  :: pgeospm_l ! surface geopotential height
    REAL(wp), DIMENSION(kproma)  :: pslm

    REAL(wp), DIMENSION(kproma)  :: zfrl, zfrw, zfri, tslm1, tsw, tsi
   
    REAL(wp), DIMENSION(kproma, klev)    :: ptm1         ! Temperature at model level [K]
    REAL(wp), DIMENSION(kproma, klev)    :: qm           ! Water vapour mixing ratio [kg kg-1]

    REAL(wp), DIMENSION(kproma, klev)    :: rh		 ! Relative humidity [%]

    REAL(wp), DIMENSION(kproma, klev)    :: cca          ! Convective cloud fraction [0-1]
    REAL(wp), DIMENSION(kproma, klev)    :: tca		 ! Total Cloud Fraction at eac hlevel, from TOA to SURFACE [0-1]

    REAL(wp), DIMENSION(kproma, klev)    :: xlm1         ! Mixing ratio large scale cloud liquid [kg kg-1]
    REAL(wp), DIMENSION(kproma, klev)    :: xim1         ! Mixing ratio large scale cloud ice [kg kg-1]

    REAL(wp), DIMENSION(:,:,:), POINTER  :: cs,cz
    REAL(wp), DIMENSION(:,:), POINTER    :: ct,cm,cr

    REAL(wp), DIMENSION(kproma) :: dum

!gbx
    REAL(wp) :: mr_lsliq(kproma,klev)
    REAL(wp) :: mr_lsice(kproma,klev)
    REAL(wp) :: mr_ccliq(kproma,klev)
    REAL(wp) :: mr_ccice(kproma,klev)
    REAL(wp) :: fl_lsrain(kproma,klev)
    REAL(wp) :: fl_lssnow(kproma,klev)
    REAL(wp) :: fl_lsgrpl(kproma,klev)
    REAL(wp) :: fl_ccrain(kproma,klev)
    REAL(wp) :: fl_ccsnow(kproma,klev)

    REAL(wp) :: Reff(kproma,klev,N_hydro)

    INTEGER :: k,j,i, itau, ipres, jk, jkb, jl, Nlr, nc

!---------------------------------------------------------------------------------

    IF (l_putdata(cospo_index) .AND. jrow .eq.1 ) THEN
       WRITE(message_text,*) 'output step', jrow
       CALL message('cosp ',message_text) 
    END IF

    !only call at radiation step
    IF ( .NOT.  l_trigrad ) THEN
       RETURN
    END IF

!---------------------------------------------------------------------------------

    ! Height of model levels [m]: ECHAM5 Top of atmosphere is 1, surface is klev (i.e 19)

   DO jl = 1,kproma
     IF (pslm(jl) .GT. 0._wp) THEN
        pgeospm_l(jl) = pgeospm( jl)
     ELSE
        pgeospm_l(jl) =  0._wp
     END IF
   END DO

   DO jk = 1,klev
      DO jl = 1,kproma
          height(jl,jk) = (pgeo(jl,jk)  +  pgeospm_l( jl) )/grav
      END DO
  END DO

    !cms++
    DO jk = 2,klev
       DO jl = 1,kproma
          height_half(jl,jk) = 0.5_wp*(height(jl,jk)+height(jl,jk-1))
       END DO
    END DO
    DO jl = 1,kproma
       height_half(jl,1) =  height(jl,1) + ( height(jl,1) - height_half(jl,2) ) 
       height_half(jl,klev+1) =   pgeospm_l( jl ) /grav 
    END DO
    !cms--

    cca(1:kproma,1:klev) = 0.0_wp  !CNam: no convective cloud fraction

!---------------------------------------------------------------------------------    
    ! Create gbx structure (Original): For values see above (CNam)
    !------------------------------------------------------------------------------    
    ! npoints: kproma
    ! nlevels: klev

    CALL construct_cosp_gridbox(radar_freq,surface_radar,use_mie_tables,use_gas_abs,do_ray,melt_lay,k2, &
         kproma,klev,Ncolumns, nclusters, &
         N_HYDRO,Nprmts_max_hydro,Naero,Nprmts_max_aero, &
         lidar_ice_type,isccp_top_height,isccp_top_height_direction,cosp_overlap,isccp_emsfc_lw, &
         use_precipitation_fluxes,use_reff, &
!CNam:         Platform,Satellite,Instrument,Nchannels,ZenAng, &
!CNam:         channels(1:Nchannels),surfem(1:Nchannels),co2,ch4,n2o,co,&
         gbx)

    !------------------------------------------------------------------------------    
    ! Code to populate input structure
    !------------------------------------------------------------------------------    
    gbx%p(1:kproma,1:klev) = p(1:kproma,klev:1:-1)
    gbx%ph(1:kproma,1:klev) = ph(1:kproma,klev+1:2:-1)
    gbx%zlev(1:kproma,1:klev) = height(1:kproma,klev:1:-1)
    gbx%zlev_half(1:kproma,1:klev) = height_half(1:kproma,klev+1:2:-1)        
    gbx%T(1:kproma,1:klev) = ptm1(1:kproma,klev:1:-1)
    gbx%q(1:kproma,1:klev) = rh(1:kproma,klev:1:-1)
    gbx%cca(1:kproma,1:klev) = cca(1:kproma,klev:1:-1)
    gbx%tca(1:kproma,1:klev) = cosp_f3d(1:kproma,klev:1:-1,jrow) !cms: tca(1:kproma,klev:1:-1) 
    gbx%psfc(1:kproma) = ph(1:kproma,klev)
    gbx%land(1:kproma) = pslm(1:kproma)        
    !gbx%mr_ozone(1:kproma,1:klev) = CNam: moved directly to rttov
    !gbx%u_wind(1:kproma,1:klev)  = CNam: moved directly to rttov
    !gbx%v_wind(1:kproma,1:klev)  = CNam: moved directly to rttov
    gbx%sunlit(1:kproma)   = cosp_sunlit(1:kproma,jrow)

    ! allocate labels for formation pathways
    gbx%labels(1:kproma,1:klev) = cslabel(1:kproma,klev:1:-1,jrow)


! hydrometeors  
! CNam: Note hydrometeors needs to be filled in (gbx%mr_hydro) with the 
!      appropriate precipitation mixing ratios in the driver program
!      when running COSP in Cloud Resolving Mode (Ncolumns=1) ICON
!      Refer to readme(Sec 5.3.11 & 6.3)
    DO jk=1,klev
       jkb = klev+1-jk
       DO jl=1,kproma
             IF ( cosp_f3d(jl,jkb,jrow) .GT. 0._wp ) THEN
               gbx%mr_hydro(jl,jk,I_LSCLIQ) = xlm1(jl,jkb) 
               gbx%mr_hydro(jl,jk,I_LSCICE) = xim1(jl,jkb) 
	     ELSE
               gbx%mr_hydro(jl,jk,I_LSCLIQ) = 0.0_wp	
               gbx%mr_hydro(jl,jk,I_LSCICE) = 0.0_wp	
             END IF
	     gbx%mr_hydro(jl,jk,I_CVCLIQ) = 0.0_wp	
	     gbx%mr_hydro(jl,jk,I_CVCICE) = 0.0_wp	
       END DO
    END DO

  	gbx%rain_ls(1:kproma,1:klev) = cosp_lsrain(1:kproma,klev:1:-1,jrow)  
  	gbx%snow_ls(1:kproma,1:klev) = cosp_lssnow(1:kproma,klev:1:-1,jrow)  
  	gbx%grpl_ls(1:kproma,1:klev) = 0.0_wp ! CNam: graupel not explcit in ECHAM  
  	gbx%rain_cv(1:kproma,1:klev) = cosp_ccrain(1:kproma,klev:1:-1,jrow)  
    	gbx%snow_cv(1:kproma,1:klev) = cosp_ccsnow(1:kproma,klev:1:-1,jrow) 

    IF ( l_fixed_reff ) THEN
       ! Constant values for cloud droplet effective radius of liquid and ice in the test
       gbx%Reff(1:kproma,1:klev,I_LSCLIQ)=10.0e-6_wp
       gbx%Reff(1:kproma,1:klev,I_LSCICE)= 40.0e-6_wp
    ELSE
       gbx%Reff(1:kproma,1:klev,I_LSCLIQ)=cosp_reffl(1:kproma,klev:1:-1,jrow)*1.e-6_wp
       gbx%Reff(1:kproma,1:klev,I_LSCICE)=cosp_reffi(1:kproma,klev:1:-1,jrow)*1.e-6_wp 
       gbx%Reff(1:kproma,1:klev,I_CVCLIQ)=cosp_reffl(1:kproma,klev:1:-1,jrow)*1.e-6_wp 
       gbx%Reff(1:kproma,1:klev,I_CVCICE)=cosp_reffi(1:kproma,klev:1:-1,jrow)*1.e-6_wp 
    END IF


    ! ISCCP simulator
    ! following cms as in ECHAM6
    IF (cfg%Lisccp_sim) THEN
       gbx%skt(1:kproma) = zfrl(1:kproma)*tslm1(1:kproma)              &
             +zfri(1:kproma)*tsi(1:kproma)                &
             +zfrw(1:kproma)*tsw(1:kproma)
    
       gbx%sh(1:kproma,1:klev) = qm(1:kproma,klev:1:-1)   

       DO jk=1,klev
          jkb = klev+1-jk
          DO jl=1,kproma
             IF ( cosp_f3d(jl,jkb,jrow) .GT. 0._wp ) THEN
                !cms: cisccp_cldtau3d is already in-cloud, (and the all-sky comment
                ! in mo_srtm.f90 close to line 1111 is not correct)
                gbx%dtau_s(jl,jk) = cisccp_cldtau3d(jl,jkb,jrow)
                gbx%dem_s(jl,jk) = cisccp_cldemi3d(jl,jkb,jrow)
             ELSE
                gbx%dtau_s(jl,jk) = 0._wp
                gbx%dem_s(jl,jk) = 0._wp
             END IF
          END DO
       END DO

    END IF


! extra diagnostic for offl test only
   IF  ( offl2dout .GT. 0) THEN
     DO jk=1,klev
       DO jl=1,kproma
          tm1_cosp(jl,jk,jrow)= ptm1 (jl,jk)
          xi_cosp(jl,jk,jrow) = xim1 (jl,jk)
          xl_cosp(jl,jk,jrow) = xlm1 (jl,jk)
       END DO
      END DO
   END IF

!---------------------------------------------------------------------------------    
    !cms: Initialize global random number generator used inside COSP
    !   This is not thread or parallel safe but affects only COSP diagnostics
    !
        CALL set_seed_random(INT((gbx%T(1,1:klev)-INT(gbx%T(1,1:klev)))*10000000))

    !++++++++++ Define new vertical grid ++++++++++
        call construct_cosp_vgrid(gbx,Nlr,use_vgrid,csat_vgrid,vgrid)

    !++++++++++ Allocate memory ++++++++++
        call construct_cosp_subgrid(kproma, Ncolumns, klev, sgx)
        call construct_cosp_sgradar(cfg,kproma,Ncolumns,klev,N_HYDRO,sgradar)
        call construct_cosp_radarstats(cfg,kproma,Ncolumns,klev,N_HYDRO,stradar)
        call construct_cosp_sglidar(cfg,kproma,Ncolumns,klev,N_HYDRO,PARASOL_NREFL,sglidar)
        call construct_cosp_lidarstats(cfg,kproma,Ncolumns,klev,N_HYDRO,PARASOL_NREFL,stlidar)
        call construct_cosp_isccpstats(cfg,kproma,Ncolumns,klev,nclusters,isccp)
        call construct_cosp_modisstats(cfg,kproma,modis)
        call construct_cosp_misrstats(cfg,kproma,misr)
!CNam:        call construct_cosp_rttov(cfg,kproma,Nchannels,rttov)

    IF (cfg%Llidar_sim) THEN
        sglidar%temp_tot(:,:) = gbx%T(:,:) !CNam: assigned sglidar%temp_tot to gbx%T.
    ENDIF

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Call simulator
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!CNam: #ifdef RTTOV
!CNam:  call cosp(cosp_overlap,kproma,klev,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
!CNam: #else
        call cosp(cosp_overlap,kproma,klev,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
!CNam: #endif

!---------------------------------------------------------------------------------    

    DO j=1,kproma
       cosp_sunlit_av(j,jrow) = cosp_sunlit_av(j,jrow)+ cosp_sunlit(j,jrow)*od_cosp
       cosp_freq(j,jrow) = cosp_freq(j,jrow)+ od_cosp
    END DO


    IF (cfg%Llidar_sim) THEN
       ! For 3D cloud fields
       DO k=1,klev
          jkb = klev+1-k
          IF ( jkb .LE.  klev ) THEN
             DO j=1,kproma
                IF (stlidar%lidarcld(j,jkb) .NE. R_UNDEF) THEN
                   cosp_lidar_cld(j,k,jrow)= cosp_lidar_cld(j,k,jrow)&
			 + stlidar%lidarcld(j,jkb)*od_cosp 
                ELSE
                   N_miss_lidar(j,k,jrow) =  N_miss_lidar(j,k,jrow) +  1._wp*od_cosp
                END IF
		!phase
                IF (stlidar%lidarcldphase(j,jkb,1) .NE. R_UNDEF) THEN
                   cosp_lidar_cldice(j,k,jrow)= cosp_lidar_cldice(j,k,jrow)&
			 + stlidar%lidarcldphase(j,jkb,1)*od_cosp 
                ELSE
                   N_miss_lidarice(j,k,jrow) =  N_miss_lidarice(j,k,jrow) +  1._wp*od_cosp
                END IF
                IF (stlidar%lidarcldphase(j,jkb,2) .NE. R_UNDEF) THEN
                   cosp_lidar_cldliq(j,k,jrow)= cosp_lidar_cldliq(j,k,jrow)&
			 + stlidar%lidarcldphase(j,jkb,2)*od_cosp 
                ELSE
                   N_miss_lidarliq(j,k,jrow) =  N_miss_lidarliq(j,k,jrow) +  1._wp*od_cosp
                END IF
                IF (stlidar%lidarcldphase(j,jkb,3) .NE. R_UNDEF) THEN
                   cosp_lidar_cldun(j,k,jrow)= cosp_lidar_cldun(j,k,jrow)&
			 + stlidar%lidarcldphase(j,jkb,3)*od_cosp 
                ELSE
                   N_miss_lidarun(j,k,jrow) =  N_miss_lidarun(j,k,jrow) +  1._wp*od_cosp
                END IF
             END DO !kproma
          END IF
       ENDDO

       ! currently the 40 level output grid does not work in ECHAM6.3-HAM2.3
       ! Therefore use histogram fields here
       DO k=1,40
          jkb = 40+1-k
          ! cloud fraction
          dum(1:kproma) = MERGE(1._wp, 0._wp, &
                                stlidar%lidarcldtmp(1:kproma,jkb,1) .NE. R_UNDEF)
          cosp_ct_cc(k)%v(1:kproma,jrow) = cosp_ct_cc(k)%v(1:kproma,jrow) &
               + od_cosp*dum(1:kproma)*stlidar%lidarcldtmp(1:kproma,jkb,1)
          cosp_ct_cc_time(k)%v(1:kproma,jrow) = cosp_ct_cc_time(k)%v(1:kproma,jrow) &
               + od_cosp*dum(1:kproma)

          ! phase fraction
          dum(1:kproma) = MERGE(1._wp, 0._wp, &
                                stlidar%lidarcldtmp(1:kproma,jkb,5) .NE. R_UNDEF)
          cosp_ct_lf(k)%v(1:kproma,jrow) = cosp_ct_lf(k)%v(1:kproma,jrow) &
               + od_cosp*dum(1:kproma)*stlidar%lidarcldtmp(1:kproma,jkb,5)
          cosp_ct_lf_time(k)%v(1:kproma,jrow) = cosp_ct_lf_time(k)%v(1:kproma,jrow) &
               + od_cosp*dum(1:kproma)
       ENDDO

       !CNam: Note 3D lidarcldtmp has 40 levels hardcoded in cosp_lmd_ipsl_stats and cosp_constants.
       !RD: NOTE that there is a bug in the if statements, this should be inversed vertically aswell!
       IF (vgrid%use_vgrid) THEN
       DO k=1,40
             jkb = 40+1-k !RD: FIX if statements
             DO j=1,kproma
                IF (stlidar%lidarcldtmp(j,jkb,1) .NE. R_UNDEF) THEN
                   cosp_lidar_cldtmp(j,k,jrow)= cosp_lidar_cldtmp(j,k,jrow)&
			 + stlidar%lidarcldtmp(j,jkb,1)*od_cosp 
                ELSE
                   N_miss_lidartmp(j,k,jrow) =  N_miss_lidartmp(j,k,jrow) +  1._wp*od_cosp
                END IF
                IF (stlidar%lidarcldtmp(j,jkb,2) .NE. R_UNDEF) THEN
                   cosp_lidar_cldtmpice(j,k,jrow)= cosp_lidar_cldtmpice(j,k,jrow)&
			 + stlidar%lidarcldtmp(j,jkb,2)*od_cosp 
                ELSE
                   N_miss_lidartmpice(j,k,jrow) =  N_miss_lidartmpice(j,k,jrow) +  1._wp*od_cosp
                END IF
                IF (stlidar%lidarcldtmp(j,jkb,3) .NE. R_UNDEF) THEN
                   cosp_lidar_cldtmpliq(j,k,jrow)= cosp_lidar_cldtmpliq(j,k,jrow)&
			 + stlidar%lidarcldtmp(j,jkb,3)*od_cosp 
                ELSE
                   N_miss_lidartmpliq(j,k,jrow) =  N_miss_lidartmpliq(j,k,jrow) +  1._wp*od_cosp
                END IF
                IF (stlidar%lidarcldtmp(j,jkb,4) .NE. R_UNDEF) THEN
                   cosp_lidar_cldtmpun(j,k,jrow)= cosp_lidar_cldtmpun(j,k,jrow)&
			 + stlidar%lidarcldtmp(j,jkb,4)*od_cosp 
                ELSE
                   N_miss_lidartmpun(j,k,jrow) =  N_miss_lidartmpun(j,k,jrow) +  1._wp*od_cosp
                END IF
             END DO !kproma
       ENDDO
       ENDIF

       ! For 'cosp_lidar_parasol_refl'
       DO k=1,PARASOL_NREFL
           ct => cosp_parasol_refl(k)%ptr2
          DO j=1,kproma
             ct(j,jrow) = ct(j,jrow) + stlidar%parasolrefl(j,k)*od_cosp 
          END DO
       END DO

       !For 'cosp_lidar_lowcloud'
       DO j=1,kproma
          IF ( stlidar%cldlayer(j,1) .NE. R_UNDEF  ) THEN
		cosp_lidar_lowcloud(j,jrow)= cosp_lidar_lowcloud(j,jrow) + &
			 stlidar%cldlayer(j,1)*od_cosp 
          ELSE
             N_miss_lidar_low(j,jrow)= N_miss_lidar_low(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( stlidar%cldlayerphase(j,1,1) .NE. R_UNDEF  ) THEN
		cosp_lidar_lowcloudice(j,jrow)= cosp_lidar_lowcloudice(j,jrow) + &
			 stlidar%cldlayerphase(j,1,1)*od_cosp 
          ELSE
             N_miss_lidar_lowice(j,jrow)= N_miss_lidar_lowice(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( stlidar%cldlayerphase(j,1,2) .NE. R_UNDEF  ) THEN
		cosp_lidar_lowcloudliq(j,jrow)= cosp_lidar_lowcloudliq(j,jrow) + &
			 stlidar%cldlayerphase(j,1,2)*od_cosp 
          ELSE
             N_miss_lidar_lowliq(j,jrow)= N_miss_lidar_lowliq(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( stlidar%cldlayerphase(j,1,3) .NE. R_UNDEF  ) THEN
		cosp_lidar_lowcloudun(j,jrow)= cosp_lidar_lowcloudun(j,jrow) + &
			 stlidar%cldlayerphase(j,1,3)*od_cosp 
          ELSE
             N_miss_lidar_lowun(j,jrow)= N_miss_lidar_lowun(j,jrow) + 1._wp*od_cosp
          END IF
       END DO


       !For 'cosp_lidar_midcloud'
       DO j=1,kproma
          IF ( stlidar%cldlayer(j,2) .NE. R_UNDEF  ) THEN
		cosp_lidar_midcloud(j,jrow)= cosp_lidar_midcloud(j,jrow) + &
			 stlidar%cldlayer(j,2)*od_cosp 
          ELSE
             N_miss_lidar_mid(j,jrow)= N_miss_lidar_mid(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( stlidar%cldlayerphase(j,2,1) .NE. R_UNDEF  ) THEN
		cosp_lidar_midcloudice(j,jrow)= cosp_lidar_midcloudice(j,jrow) + &
			 stlidar%cldlayerphase(j,2,1)*od_cosp 
          ELSE
             N_miss_lidar_midice(j,jrow)= N_miss_lidar_midice(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( stlidar%cldlayerphase(j,2,2) .NE. R_UNDEF  ) THEN
		cosp_lidar_midcloudliq(j,jrow)= cosp_lidar_midcloudliq(j,jrow) + &
			 stlidar%cldlayerphase(j,2,2)*od_cosp 
          ELSE
             N_miss_lidar_midliq(j,jrow)= N_miss_lidar_midliq(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( stlidar%cldlayerphase(j,2,3) .NE. R_UNDEF  ) THEN
		cosp_lidar_midcloudun(j,jrow)= cosp_lidar_midcloudun(j,jrow) + &
			 stlidar%cldlayerphase(j,2,3)*od_cosp 
          ELSE
             N_miss_lidar_midun(j,jrow)= N_miss_lidar_midun(j,jrow) + 1._wp*od_cosp
          END IF
       END DO

       !For 'cosp_lidar_highcloud'
       DO j=1,kproma
		cosp_lidar_highcloud(j,jrow)= cosp_lidar_highcloud(j,jrow) + &
			stlidar%cldlayer(j,3)*od_cosp 		
		cosp_lidar_highcloudice(j,jrow)= cosp_lidar_highcloudice(j,jrow) + &
			 stlidar%cldlayerphase(j,3,1)*od_cosp 
		cosp_lidar_highcloudliq(j,jrow)= cosp_lidar_highcloudliq(j,jrow) + &
			 stlidar%cldlayerphase(j,3,2)*od_cosp 
		cosp_lidar_highcloudun(j,jrow)= cosp_lidar_highcloudun(j,jrow) + &
			 stlidar%cldlayerphase(j,3,3)*od_cosp 
       ENDDO

       !For 'cosp_lidar_totalcloud'
       DO j=1,kproma
		cosp_lidar_totalcloud(j,jrow)= cosp_lidar_totalcloud(j,jrow) + &
			stlidar%cldlayer(j,4)*od_cosp
		cosp_lidar_totalcloudice(j,jrow)= cosp_lidar_totalcloudice(j,jrow) + &
			 stlidar%cldlayerphase(j,4,1)*od_cosp 
		cosp_lidar_totalcloudliq(j,jrow)= cosp_lidar_totalcloudliq(j,jrow) + &
			 stlidar%cldlayerphase(j,4,2)*od_cosp 
		cosp_lidar_totalcloudun(j,jrow)= cosp_lidar_totalcloudun(j,jrow) + &
			 stlidar%cldlayerphase(j,4,3)*od_cosp 
       ENDDO

       !For CFAD of Lidar scattering ratio
     IF (cfg%LcfadLidarsr532) THEN
       DO i=1,15
          DO k=1,klev
             jkb = klev+1-k
             IF ( jkb .LE. klev ) THEN
                DO j=1,kproma
                   cs => cosp_cfadsr(i)%ptr
                   cs(j,k,jrow)= cs(j,k,jrow) + stlidar%cfad_sr(j,i,jkb) * od_cosp 
                END DO
             END IF
          ENDDO
       ENDDO
     END IF

    END IF !Llidar_sim
   
	!---------------------------------------------------------------------------------    

   IF (cfg%Lradar_sim) THEN
     ! For 'cosp_radar_lidar_tcc'
       DO j=1,kproma
		cosp_radar_lidar_tcc(j,jrow)= cosp_radar_lidar_tcc(j,jrow) + &
			stradar%radar_lidar_tcc(j)*od_cosp 
       ENDDO

     ! For 'cosp_lidar_only_freq_cloud'
       DO k=1,klev
          jkb = klev+1-k
          IF ( jkb .LE.  klev ) THEN
             DO j=1,kproma
                IF (stradar%lidar_only_freq_cloud(j,jkb) .NE.  R_UNDEF  ) THEN
                  cosp_lidar_only_freq_cloud(j,k,jrow)= cosp_lidar_only_freq_cloud(j,k,jrow) + &
			stradar%lidar_only_freq_cloud(j,jkb)*od_cosp
                ELSE
                  N_miss_radar(j,k,jrow) =  N_miss_radar(j,k,jrow) + 1._wp*od_cosp
                ENDIF
             ENDDO
         ENDIF
       ENDDO

     ! For CFAD of Radar reflectivity
     IF (cfg%LcfadDbze94) THEN
       DO i = 1,15
          DO k=1,klev
             jkb = klev+1-k
             IF ( jkb .LE.  klev ) THEN
                DO j=1,kproma
                   cz => cosp_cfadze(i)%ptr
                   cz(j,k,jrow)= cz(j,k,jrow) + stradar%cfad_ze(j,i,jkb)*od_cosp
                END DO
             END IF
          ENDDO
       ENDDO
     END IF

    END IF !Lradar_sim

	!---------------------------------------------------------------------------------    

    IF (cfg%Lisccp_sim) THEN

       DO j=1,kproma
          cosp_isccp_totalcloud(j,jrow) = cosp_isccp_totalcloud(j,jrow)+ isccp%totalcldarea(j) * od_cosp
       ENDDO

       ! weight by totalcloud for time mean
       DO j=1,kproma
          cosp_isccp_meanptop(j,jrow) = cosp_isccp_meanptop(j,jrow) + &
		isccp%totalcldarea(j) * isccp%meanptop(j) * od_cosp
       ENDDO

       ! weight by totalcloud for time mean
       DO j=1,kproma
          cosp_isccp_meanalbedocld(j,jrow) = cosp_isccp_meanalbedocld(j,jrow) + &
               isccp%totalcldarea(j) * isccp%meanalbedocld(j) * od_cosp
       ENDDO

       ! weight by totalcloud for time mean
       DO j=1,kproma
          cosp_isccp_meantaucld(j,jrow) = cosp_isccp_meantaucld(j,jrow) + &
               isccp%totalcldarea(j) * isccp%meantaucld(j) * od_cosp
       ENDDO

       i=0
       k=0
       DO itau=1,7
          DO ipres=1,7
             i=i+1
             ct => isccp_cldtypes(i)%ptr2
             ct(1:kproma,jrow) =  ct(1:kproma,jrow) + isccp%fq_isccp(1:kproma,itau,ipres) * od_cosp
             ! include convective contribution
             DO nc=1,nclusters+1
                k=k+1
                isccp_ct(k)%v(1:kproma,jrow) = isccp_ct(k)%v(1:kproma,jrow) &
                                             + isccp%fq_isccp_ct(1:kproma,nc,itau,ipres) * od_cosp
             ENDDO !nc
          ENDDO
       ENDDO
       isccp_freq(1:kproma,jrow) = isccp_freq(1:kproma,jrow) + od_cosp
       isccp_sunlit(1:kproma,jrow) = isccp_sunlit(1:kproma,jrow) + cosp_sunlit(1:kproma,jrow)*od_cosp

    END IF !Lisccp_sim

	!---------------------------------------------------------------------------------    

    IF (cfg%Lmodis_sim) THEN

       !total, liq, ice cloud
       DO j=1,kproma
          IF ( modis%Cloud_Fraction_Total_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_clt(j,jrow) = cosp_modis_clt(j,jrow)+ modis%Cloud_Fraction_Total_Mean(j) * od_cosp
          ELSE
             N_miss_modis_clt(j,jrow)= N_miss_modis_clt(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( modis%Cloud_Fraction_Water_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_clw(j,jrow) = cosp_modis_clw(j,jrow)+ modis%Cloud_Fraction_Water_Mean(j) * od_cosp
          ELSE
             N_miss_modis_clw(j,jrow)= N_miss_modis_clw(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( modis%Cloud_Fraction_Ice_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_cli(j,jrow) = cosp_modis_cli(j,jrow)+ modis%Cloud_Fraction_Ice_Mean(j) * od_cosp
          ELSE
             N_miss_modis_cli(j,jrow)= N_miss_modis_cli(j,jrow) + 1._wp*od_cosp
          END IF
       ENDDO

       !high, mid, low cloud
       DO j=1,kproma
          IF ( modis%Cloud_Fraction_High_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_clh(j,jrow) = cosp_modis_clh(j,jrow)+ modis%Cloud_Fraction_High_Mean(j) * od_cosp
          ELSE
             N_miss_modis_clh(j,jrow)= N_miss_modis_clh(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( modis%Cloud_Fraction_Mid_Mean(j) .NE. R_UNDEF  ) THEN
               cosp_modis_clm(j,jrow) = cosp_modis_clm(j,jrow)+ modis%Cloud_Fraction_Mid_Mean(j) * od_cosp
          ELSE
             N_miss_modis_clm(j,jrow)= N_miss_modis_clm(j,jrow) + 1._wp*od_cosp
          END IF
          IF ( modis%Cloud_Fraction_Low_Mean(j) .NE. R_UNDEF  ) THEN
               cosp_modis_cll(j,jrow) = cosp_modis_cll(j,jrow)+ modis%Cloud_Fraction_Low_Mean(j) * od_cosp
          ELSE
             N_miss_modis_cll(j,jrow)= N_miss_modis_cll(j,jrow) + 1._wp*od_cosp
          END IF
       ENDDO

       !tau !CNam: weight by totalcloud for time mean following cms
       DO j=1,kproma 
          IF ( modis%Optical_Thickness_Total_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_taut(j,jrow) = cosp_modis_taut(j,jrow)+ modis%Cloud_Fraction_Total_Mean(j) &
		 * modis%Optical_Thickness_Total_Mean(j) * od_cosp
          END IF
          IF ( modis%Optical_Thickness_Water_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_taul(j,jrow) = cosp_modis_taul(j,jrow)+ modis%Cloud_Fraction_Water_Mean(j) &
		* modis%Optical_Thickness_Water_Mean(j) * od_cosp
          END IF
          IF ( modis%Optical_Thickness_Ice_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_taui(j,jrow) = cosp_modis_taui(j,jrow)+ modis%Cloud_Fraction_Ice_Mean(j)  &
		* modis%Optical_Thickness_Ice_Mean(j) * od_cosp
          END IF
       ENDDO

       !tau log !CNam: weight by totalcloud for time mean following cms
       DO j=1,kproma 
          IF ( modis%Optical_Thickness_Total_LogMean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_tautlog(j,jrow) = cosp_modis_tautlog(j,jrow)+ modis%Cloud_Fraction_Total_Mean(j) &
		 * modis%Optical_Thickness_Total_LogMean(j) * od_cosp
          END IF
          IF ( modis%Optical_Thickness_Water_LogMean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_taullog(j,jrow) = cosp_modis_taullog(j,jrow)+ modis%Cloud_Fraction_Water_Mean(j) &
		 * modis%Optical_Thickness_Water_LogMean(j) * od_cosp
          END IF
          IF ( modis%Optical_Thickness_Ice_LogMean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_tauilog(j,jrow) = cosp_modis_tauilog(j,jrow)+ modis%Cloud_Fraction_Ice_Mean(j) &
		 * modis%Optical_Thickness_Ice_LogMean(j) * od_cosp
          END IF
       ENDDO

       !reff !CNam: weight by totalcloud for time mean following cms
       DO j=1,kproma 
          IF ( modis%Cloud_Particle_Size_Water_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_reffclw(j,jrow) = cosp_modis_reffclw(j,jrow)+ modis%Cloud_Fraction_Water_Mean(j)  &
		* modis%Cloud_Particle_Size_Water_Mean(j) * od_cosp
          END IF
          IF ( modis%Cloud_Particle_Size_Ice_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_reffcli(j,jrow) = cosp_modis_reffcli(j,jrow)+ modis%Cloud_Fraction_Ice_Mean(j)  &
		* modis%Cloud_Particle_Size_Ice_Mean(j) * od_cosp
          END IF
       ENDDO

       !pct, lwp, iwp !CNam: weight by totalcloud for time mean following cms
       DO j=1,kproma 
          IF ( modis%Cloud_Top_Pressure_Total_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_pct(j,jrow) = cosp_modis_pct(j,jrow)+ modis%Cloud_Fraction_Total_Mean(j)  &
		*  modis%Cloud_Top_Pressure_Total_Mean(j) * od_cosp
          END IF
          IF ( modis%Liquid_Water_Path_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_lwp(j,jrow) = cosp_modis_lwp(j,jrow)+ modis%Cloud_Fraction_Water_Mean(j)  &
		*  modis%Liquid_Water_Path_Mean(j) * od_cosp
          END IF
          IF ( modis%Ice_Water_Path_Mean(j) .NE. R_UNDEF  ) THEN
             cosp_modis_iwp(j,jrow) = cosp_modis_iwp(j,jrow)+ modis%Cloud_Fraction_Ice_Mean(j)  &
		*  modis%Ice_Water_Path_Mean(j) * od_cosp
          END IF
       ENDDO

       !clmodis
       i=0
       DO itau=1,7
          DO ipres=1,7
             i=i+1
             cm => modis_cldtypes(i)%ptr2
             cm(1:kproma,jrow) =  cm(1:kproma,jrow) + &
		modis%Optical_Thickness_vs_Cloud_Top_Pressure(1:kproma,itau,ipres)*100 * od_cosp
          ENDDO
       ENDDO

    END IF !Lmodis_sim

	!---------------------------------------------------------------------------------    

    IF (cfg%Lmisr_sim) THEN
       !clmisr
       i=0
       DO itau=1,7
          DO ipres=1,7
             i=i+1
             cr => misr_cldtypes(i)%ptr2
             cr(1:kproma,jrow) =  cr(1:kproma,jrow) + misr%fq_misr(1:kproma,itau,ipres)*100 * od_cosp
          ENDDO
       ENDDO

    END IF !Lmisr_sim

	!---------------------------------------------------------------------------------    
    
   !  Deallocate dynamic memory
        call free_cosp_gridbox(gbx)
        call free_cosp_subgrid(sgx)
        call free_cosp_sgradar(sgradar)
        call free_cosp_radarstats(stradar)
        call free_cosp_sglidar(sglidar)
        call free_cosp_lidarstats(stlidar)
        call free_cosp_isccpstats(isccp)
        call free_cosp_misrstats(misr)
        call free_cosp_modisstats(modis)
!CNam:        call free_cosp_rttovstats(rttov)

!>>KH change name for compatibility with echam-ham with prev version of cosp
  !END SUBROUTINE echam_cospsimulator
  END SUBROUTINE call_cospsimulator
!<<KH
!--------------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------------------

!>>KH changed name of module for compatibility with echam-ham
!END MODULE mo_cosp_echam
END MODULE mo_cosp_simulator
!<<KH
