! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!CNam: Did NOT include *_CPSECTION or 'PRINT SUBROUTINES'


MODULE MO_COSP_TYPES
    USE mo_kind,		ONLY: wp
    USE mo_cosp_constants
    USE mo_cosp_utils

    use mo_cosp_radar_simulator_types, only: class_param, nd, mt_nd, dmax, dmin

    IMPLICIT NONE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------- DERIVED TYPES ----------------------------    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Configuration choices (simulators, variables)
  TYPE COSP_CONFIG
     logical :: Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim,Lstats,Lwrite_output, &
                Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
                LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp,Lcllcalipso, &
                Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp,Lcltisccp, &
                Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
                Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun, &
                Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun, &
	        Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun, &
                Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun, &
                Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun, &
                Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun, &
                Lfracout,LlidarBetaMol532,Ltbrttov, &
		Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis, &
		Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
                Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
                Liwpmodis,Lclmodis

     character(len=32) :: out_list(N_OUT_LIST)
  END TYPE COSP_CONFIG
  
  ! Outputs from RTTOV
  TYPE COSP_RTTOV
     ! Dimensions
     integer :: Npoints   ! Number of gridpoints
     integer :: Nchan     ! Number of channels
     
     ! Brightness temperatures (Npoints,Nchan)
     real(wp),pointer :: tbs(:,:)
     
  END TYPE COSP_RTTOV
  
  ! Outputs from MISR simulator
  TYPE cosp_misrstats
     ! Dimensions
     integer :: Npoints   ! Number of gridpoints
     integer :: Ntau      ! Number of tau intervals
     integer :: Nlevels   ! Number of cth levels

     ! --- (npoints,ntau,nlevels)
     !  the fraction of the model grid box covered by each of the MISR cloud types
     real(wp),pointer :: fq_MISR(:,:,:)  
     
     ! --- (npoints)
     real(wp),pointer :: MISR_meanztop(:), MISR_cldarea(:)
     ! --- (npoints,nlevels)
     real(wp),pointer :: MISR_dist_model_layertops(:,:)
  END TYPE cosp_misrstats

  ! Outputs from ISCCP simulator
  TYPE COSP_ISCCPSTATS
     ! Dimensions
     integer :: Npoints   ! Number of gridpoints
     integer :: Ncolumns  ! Number of columns
     integer :: Nlevels   ! Number of levels
     integer :: nclusters ! Number formation pathway cloud types

    
     ! --- (npoints,tau=7,pressure=7)
     !  the fraction of the model grid box covered by each of the 49 ISCCP D level cloud types
     real(wp),pointer :: fq_isccp(:,:,:)
     !  ISCCP cloud types per formation pathway (npoints,nclusters,tau=7,pressure=7)
     real(wp),pointer :: fq_isccp_ct(:,:,:,:)
     
     ! --- (npoints) ---
     ! The fraction of model grid box columns with cloud somewhere in them.
     ! This should equal the sum over all entries of fq_isccp
     real(wp),pointer :: totalcldarea(:)
     ! mean all-sky 10.5 micron brightness temperature
     real(wp),pointer ::  meantb(:)
     ! mean clear-sky 10.5 micron brightness temperature
     real(wp),pointer ::  meantbclr(:)
     
     ! The following three means are averages over the cloudy areas only.  If no
     ! clouds are in grid box all three quantities should equal zero.
     
     !  mean cloud top pressure (mb) - linear averaging in cloud top pressure.
     real(wp),pointer :: meanptop(:)
     !  mean optical thickness linear averaging in albedo performed.
     real(wp),pointer :: meantaucld(:)
     ! mean cloud albedo. linear averaging in albedo performed 
     real(wp),pointer :: meanalbedocld(:)  
     
     !--- (npoints,ncol) ---
     !  optical thickness in each column     
     real(wp),pointer :: boxtau(:,:)
     !  cloud top pressure (mb) in each column
     real(wp),pointer :: boxptop(:,:)        
  END TYPE COSP_ISCCPSTATS
  
  ! Summary statistics from radar
  TYPE COSP_VGRID
    logical :: use_vgrid ! Logical flag that indicates change of grid
    logical :: csat_vgrid ! Flag for Cloudsat grid
    integer :: Npoints   ! Number of sampled points
    integer :: Ncolumns  ! Number of subgrid columns
    integer :: Nlevels   ! Number of model levels
    integer :: Nlvgrid   ! Number of levels of new grid
    ! Array with dimensions (Nlvgrid)
    real(wp), dimension(:), pointer :: z,zl,zu ! Height and lower and upper boundaries of new levels
    ! Array with dimensions (Nlevels)
    real(wp), dimension(:), pointer :: mz,mzl,mzu ! Height and lower and upper boundaries of model levels
  END TYPE COSP_VGRID
  
  ! Output data from lidar code
  TYPE COSP_SGLIDAR
    ! Dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Ncolumns  ! Number of columns
    integer :: Nlevels   ! Number of levels
    integer :: Nhydro    ! Number of hydrometeors    
    integer :: Nrefl     ! Number of parasol reflectances
    ! Arrays with dimensions (Npoints,Nlevels)
    real(wp),dimension(:,:),pointer :: beta_mol   ! Molecular backscatter
    real(wp),dimension(:,:),pointer :: temp_tot
    ! Arrays with dimensions (Npoints,Ncolumns,Nlevels)
    real(wp),dimension(:,:,:),pointer :: betaperp_tot   ! Total backscattered signal
    real(wp),dimension(:,:,:),pointer :: beta_tot   ! Total backscattered signal
    real(wp),dimension(:,:,:),pointer :: tau_tot    ! Optical thickness integrated from top to level z
    ! Arrays with dimensions (Npoints,Ncolumns,Nrefl)
    real(wp),dimension(:,:,:),pointer :: refl       ! parasol reflectances
  END TYPE COSP_SGLIDAR
  
  ! Output data from radar code
  TYPE COSP_SGRADAR
    ! Dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Ncolumns  ! Number of columns
    integer :: Nlevels   ! Number of levels
    integer :: Nhydro    ! Number of hydrometeors
    ! output vertical levels: spaceborne radar -> from TOA to SURFACE
    ! Arrays with dimensions (Npoints,Nlevels)
    real(wp),dimension(:,:),pointer :: att_gas ! 2-way attenuation by gases [dBZ]
    ! Arrays with dimensions (Npoints,Ncolumns,Nlevels)
    real(wp),dimension(:,:,:),pointer :: Ze_tot ! Effective reflectivity factor [dBZ]
  END TYPE COSP_SGRADAR

  
  ! Summary statistics from radar
  TYPE COSP_RADARSTATS
    integer :: Npoints  ! Number of sampled points
    integer :: Ncolumns ! Number of subgrid columns
    integer :: Nlevels  ! Number of model levels
    integer :: Nhydro   ! Number of hydrometeors
    ! Array with dimensions (Npoints,dBZe_bins,Nlevels)
    real(wp), dimension(:,:,:),pointer :: cfad_ze ! Ze CFAD
    ! Array with dimensions (Npoints)
    real(wp),dimension(:),pointer :: radar_lidar_tcc ! Radar&lidar total cloud amount, grid-box scale
    ! Arrays with dimensions (Npoints,Nlevels)
    real(wp), dimension(:,:),pointer :: lidar_only_freq_cloud
  END TYPE COSP_RADARSTATS

  ! Summary statistics from lidar
  TYPE COSP_LIDARSTATS
    integer :: Npoints  ! Number of sampled points
    integer :: Ncolumns ! Number of subgrid columns
    integer :: Nlevels  ! Number of model levels
    integer :: Nhydro   ! Number of hydrometeors
    integer :: Nrefl    ! Number of parasol reflectances
    
    ! Arrays with dimensions (SR_BINS)
    real(wp), dimension(:),pointer :: srbval ! SR bins in cfad_sr
    ! Arrays with dimensions (Npoints,SR_BINS,Nlevels)
    real(wp), dimension(:,:,:),pointer :: cfad_sr   ! CFAD of scattering ratio
    ! Arrays with dimensions (Npoints,Nlevels)
    real(wp), dimension(:,:),pointer :: lidarcld    ! 3D "lidar" cloud fraction 
    ! Arrays with dimensions (Npoints,LIDAR_NCAT)
    real(wp), dimension(:,:),pointer :: cldlayer      ! low, mid, high-level, total lidar cloud cover
    ! Arrays with dimensions (Npoints,Nlevels,Nphase)
    real(wp), dimension(:,:,:),pointer :: lidarcldphase    ! 3D "lidar" phase cloud fraction 
    ! Arrays with dimensions (Npoints,LIDAR_NCAT,Nphase)
    real(wp), dimension(:,:,:),pointer :: cldlayerphase      ! low, mid, high-level lidar phase cloud cover
    ! Arrays with dimensions (Npoints,Ntemps,Nphase)
    real(wp), dimension(:,:,:),pointer :: lidarcldtmp    ! 3D "lidar" phase cloud temperature
    ! Arrays with dimensions (Npoints,PARASOL_NREFL)
    real(wp), dimension(:,:),pointer :: parasolrefl   ! mean parasol reflectance

  END TYPE COSP_LIDARSTATS

    
  ! Input data for simulator. Subgrid scale.
  ! Input data from SURFACE to TOA
  TYPE COSP_SUBGRID
    ! Dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Ncolumns  ! Number of columns
    integer :: Nlevels   ! Number of levels
    integer :: Nhydro    ! Number of hydrometeors
    
    real(wp),dimension(:,:,:),pointer :: prec_frac  ! Subgrid precip array. Dimensions (Npoints,Ncolumns,Nlevels)
    real(wp),dimension(:,:,:),pointer :: frac_out  ! Subgrid cloud array. Dimensions (Npoints,Ncolumns,Nlevels)
  END TYPE COSP_SUBGRID

  ! Input data for simulator at Subgrid scale.
  ! Used on a reduced number of points
  TYPE COSP_SGHYDRO
    ! Dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Ncolumns  ! Number of columns
    integer :: Nlevels   ! Number of levels
    integer :: Nhydro    ! Number of hydrometeors
    real(wp),dimension(:,:,:,:),pointer :: mr_hydro ! Mixing ratio of each hydrometeor 
                                                ! (Npoints,Ncolumns,Nlevels,Nhydro) [kg/kg]
    real(wp),dimension(:,:,:,:),pointer :: Reff     ! Effective Radius of each hydrometeor
                                                ! (Reff==0 means use default size)   
                                                ! (Npoints,Ncolumns,Nlevels,Nhydro) [m]
    real(wp),dimension(:,:,:,:),pointer :: Np       ! Total # concentration each hydrometeor 
                                                ! (Optional, ignored if Reff > 0).
                                                ! (Npoints,Ncolumns,Nlevels,Nhydro) [#/kg]
                                                ! Np = Ntot / rho_a  = [#/m^3] / [kg/m^3) 
                                                ! added by Roj with Quickbeam V3
  END TYPE COSP_SGHYDRO
  
  ! Input data for simulator. Gridbox scale.
  TYPE COSP_GRIDBOX
    ! Scalars and dimensions
    integer :: Npoints   ! Number of gridpoints
    integer :: Nlevels   ! Number of levels
    integer :: Ncolumns  ! Number of columns
    integer :: nclusters ! Number of columns
    integer :: Nhydro    ! Number of hydrometeors
    integer :: Nprmts_max_hydro    ! Max number of parameters for hydrometeor size distributions
    integer :: Naero    ! Number of aerosol species
    integer :: Nprmts_max_aero    ! Max number of parameters for aerosol size distributions
!CNam    integer :: Npoints_it   ! Max number of gridpoints to be processed in one iteration
    
    ! Time [days]
    double precision :: time
    double precision :: time_bnds(2)
    
    ! Radar ancillary info
    real(wp) :: radar_freq, & ! Radar frequency [GHz]
            k2 ! |K|^2, -1=use frequency dependent default
    integer :: surface_radar, & ! surface=1, spaceborne=0
           use_mie_tables, & ! use a precomputed loopup table? yes=1,no=0
           use_gas_abs, & ! include gaseous absorption? yes=1,no=0
           do_ray, & ! calculate/output Rayleigh refl=1, not=0
           melt_lay ! melting layer model off=0, on=1
 
    ! structures used by radar simulator that need to be set only ONCE per radar configuration (e.g. freq, pointing direction) ... added by roj Feb 2008
    type(class_param) ::  hp    ! structure used by radar simulator to store Ze and N scaling constants and other information
    integer :: nsizes       ! number of discrete drop sizes (um) used to represent the distribution
    
    ! Lidar
    integer :: lidar_ice_type !ice particle shape hypothesis in lidar calculations 
                              !(ice_type=0 for spheres, ice_type=1 for non spherical particles)
    
    ! Radar
    logical ::  use_precipitation_fluxes  ! True if precipitation fluxes are input to the algorithm 
    logical ::  use_reff          ! True if Reff is to be used by radar (memory not allocated
    
    
    ! Geolocation (Npoints)
    real(wp),dimension(:),pointer :: toffset   ! Time offset of esch point from the value in time
    real(wp),dimension(:),pointer :: longitude ! longitude [degrees East]
    real(wp),dimension(:),pointer :: latitude  ! latitude [deg North]
    ! Gridbox information (Npoints,Nlevels)
    real(wp),dimension(:,:),pointer :: zlev ! Height of model levels [m]
    real(wp),dimension(:,:),pointer :: zlev_half ! Height at half model levels [m] (Bottom of model layer)
    real(wp),dimension(:,:),pointer :: dlev ! Depth of model levels  [m]
    real(wp),dimension(:,:),pointer :: p  ! Pressure at full model levels [Pa]
    real(wp),dimension(:,:),pointer :: ph ! Pressure at half model levels [Pa]
    real(wp),dimension(:,:),pointer :: T ! Temperature at model levels [K]
    real(wp),dimension(:,:),pointer :: q  ! Relative humidity to water (%)
    real(wp),dimension(:,:),pointer :: sh ! Specific humidity to water [kg/kg]
    real(wp),dimension(:,:),pointer :: dtau_s ! mean 0.67 micron optical depth of stratiform
                                          !  clouds in each model level
                                          !  NOTE:  this the cloud optical depth of only the
                                          !  cloudy part of the grid box, it is not weighted
                                          !  with the 0 cloud optical depth of the clear
                                          !         part of the grid box
    real(wp),dimension(:,:),pointer :: dtau_c !  mean 0.67 micron optical depth of convective
                                          !  clouds in each model level.  Same note applies as in dtau_s.
    real(wp),dimension(:,:),pointer :: dem_s  !  10.5 micron longwave emissivity of stratiform
                                          !  clouds in each model level.  Same note applies as in dtau_s.
    real(wp),dimension(:,:),pointer :: dem_c  !  10.5 micron longwave emissivity of convective
                                          !  clouds in each model level.  Same note applies as in dtau_s.
    real(wp),dimension(:,:),pointer :: mr_ozone !  Ozone mass mixing ratio [kg/kg]
    real(wp),dimension(:,:),pointer :: labels ! cloud formation labels

    ! Point information (Npoints)
    real(wp),dimension(:),pointer :: land !Landmask [0 - Ocean, 1 - Land]
    real(wp),dimension(:),pointer :: psfc !Surface pressure [Pa]
    real(wp),dimension(:),pointer :: sunlit ! (npoints) 1 for day points, 0 for nightime
    real(wp),dimension(:),pointer :: skt  ! Skin temperature (K)
    real(wp),dimension(:),pointer :: u_wind  ! eastward wind [m s-1]
    real(wp),dimension(:),pointer :: v_wind  ! northward wind [m s-1]

    ! TOTAL and CONV cloud fraction for SCOPS
    real(wp),dimension(:,:),pointer :: tca ! Total cloud fraction
    real(wp),dimension(:,:),pointer :: cca ! Convective cloud fraction
    ! Precipitation fluxes on model levels
    real(wp),dimension(:,:),pointer :: rain_ls ! large-scale precipitation flux of rain [kg/m2.s]
    real(wp),dimension(:,:),pointer :: rain_cv ! convective precipitation flux of rain [kg/m2.s]
    real(wp),dimension(:,:),pointer :: snow_ls ! large-scale precipitation flux of snow [kg/m2.s]
    real(wp),dimension(:,:),pointer :: snow_cv ! convective precipitation flux of snow [kg/m2.s]
    real(wp),dimension(:,:),pointer :: grpl_ls ! large-scale precipitation flux of graupel [kg/m2.s]
    ! Hydrometeors concentration and distribution parameters
!     real(wp),dimension(:,:,:),pointer :: fr_hydro ! Fraction of the gridbox occupied by each hydrometeor (Npoints,Nlevels,Nhydro)
    real(wp),dimension(:,:,:),pointer :: mr_hydro ! Mixing ratio of each hydrometeor (Npoints,Nlevels,Nhydro) [kg/kg]
    real(wp),dimension(:,:),pointer   :: dist_prmts_hydro !Distributional parameters for hydrometeors (Nprmts_max_hydro,Nhydro)

    ! Effective radius [m]. (Npoints,Nlevels,Nhydro) -- OPTIONAL, value of 0 mean use fixed default
    real(wp),dimension(:,:,:),pointer :: Reff

    ! Total Number Concentration [#/kg]. (Npoints,Nlevels,Nhydro) -- OPTIONAL, value of 0 mean use fixed default
    real(wp),dimension(:,:,:),pointer :: Np ! added by Roj with Quickbeam V3
 
    ! Aerosols concentration and distribution parameters
    real(wp),dimension(:,:,:),pointer :: conc_aero ! Aerosol concentration for each species (Npoints,Nlevels,Naero)
    integer,dimension(:),pointer :: dist_type_aero ! Particle size distribution type for each aerosol species (Naero)
    real(wp),dimension(:,:,:,:),pointer :: dist_prmts_aero ! Distributional parameters for aerosols 
                                                       ! (Npoints,Nlevels,Nprmts_max_aero,Naero)
    ! ISCCP simulator inputs
    integer :: isccp_top_height !  1 = adjust top height using both a computed
                                !  infrared brightness temperature and the visible
                                !  optical depth to adjust cloud top pressure. Note
                                !  that this calculation is most appropriate to compare
                                !  to ISCCP data during sunlit hours.
                                !  2 = do not adjust top height, that is cloud top
                                !  pressure is the actual cloud top pressure
                                !  in the model
                                !  3 = adjust top height using only the computed
                                !  infrared brightness temperature. Note that this
                                !  calculation is most appropriate to compare to ISCCP
                                !  IR only algortihm (i.e. you can compare to nighttime
                                !  ISCCP data with this option)
    integer :: isccp_top_height_direction ! direction for finding atmosphere pressure level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                 ! with interpolated temperature equal to the radiance 
                                 ! determined cloud-top temperature
                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                 ! 1 = default setting, and matches all versions of 
                                 ! ISCCP simulator with versions numbers 3.5.1 and lower
                                 ! 2 = experimental setting  
    integer :: cosp_overlap !  overlap type (1=max, 2=rand, 3=max/rand)
    real(wp) :: isccp_emsfc_lw      ! 10.5 micron emissivity of surface (fraction)
  
    ! RTTOV inputs/options
    integer :: plat      ! satellite platform
    integer :: sat       ! satellite
    integer :: inst      ! instrument
    integer :: Nchan     ! Number of channels to be computed
    integer, dimension(:), pointer :: Ichan   ! Channel numbers
    real(wp), dimension(:), pointer :: Surfem  ! Surface emissivity
    real(wp) :: ZenAng ! Satellite Zenith Angles
    real(wp) :: co2,ch4,n2o,co ! Mixing ratios of trace gases

  END TYPE COSP_GRIDBOX
 
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_RTTOV -------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!CNam:  SUBROUTINE CONSTRUCT_COSP_RTTOV(cfg,Npoints,Nchan,x)
!    type(cosp_config),intent(in) :: cfg ! Configuration options
!    integer,intent(in) :: Npoints  ! Number of sampled points
!    integer,intent(in) :: Nchan ! Number of channels
!    type(cosp_rttov),intent(out) :: x
!    ! Local variables
!    integer :: i,j
!    
!    ! Allocate minumum storage if simulator not used
!    if (cfg%Lrttov_sim) then
!      i = Npoints
!      j = Nchan
!    else
!      i = 1
!      j = 1
!    endif
!    x%Npoints  = i
!    x%Nchan    = j
!      
!    ! --- Allocate arrays ---
!    allocate(x%tbs(i, j))
!    ! --- Initialise to zero ---
!    x%tbs     = 0.0_wp
!  END SUBROUTINE CONSTRUCT_COSP_RTTOV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_RTTOV ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  SUBROUTINE FREE_COSP_RTTOV(x)
!    type(cosp_rttov),intent(inout) :: x
!    
!    ! --- Deallocate arrays ---
!    deallocate(x%tbs)
!  END SUBROUTINE FREE_COSP_RTTOV
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_MISRSTATS ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_MISRSTATS(cfg,Npoints,x) !CNam: renamed to misrstats
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints   ! Number of gridpoints
    type(cosp_misrstats),intent(out) :: x
    ! Local variables
    integer :: i,j,k
    
   
    ! Allocate minumum storage if simulator not used
    if (cfg%Lmisr_sim) then
      i = Npoints
      j = 7
      k = MISR_N_CTH
    else
      i = 1
      j = 1
      k = 1
    endif
    
    ! Dimensions
    x%Npoints = i
    x%Ntau    = j
    x%Nlevels = k
    
    ! allocate space for MISR simulator outputs ...
    allocate(x%fq_MISR(i,j,k), x%MISR_meanztop(i),x%MISR_cldarea(i), x%MISR_dist_model_layertops(i,k))
    x%fq_MISR = 0.0_wp
    x%MISR_meanztop = 0.0_wp
    x%MISR_cldarea = 0.0_wp
    x%MISR_dist_model_layertops = 0.0_wp
    
  END SUBROUTINE CONSTRUCT_COSP_MISRSTATS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_MISRSTATS ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_MISRSTATS(x)
    type(cosp_misrstats),intent(inout) :: x
    deallocate(x%fq_MISR, x%MISR_meanztop,x%MISR_cldarea, x%MISR_dist_model_layertops)
    
  END SUBROUTINE FREE_COSP_MISRSTATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_ISCCPSTATS ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_ISCCPSTATS(cfg,Npoints,Ncolumns,Nlevels,nclusters,x) !CNam: renamed to isccpstats
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    integer,intent(in) :: nclusters ! Number formation pathway clusters
    type(cosp_isccpstats),intent(out) :: x
    ! Local variables
    integer :: i,j,k,nc
    
    ! Allocate minumum storage if simulator not used
    if (cfg%Lisccp_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
      nc = nclusters
    else
      i = 1
      j = 1
      k = 1
      nc = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    x%nclusters  = nc
    
    ! --- Allocate arrays ---
    allocate(x%fq_isccp(i,7,7), x%fq_isccp_ct(i,nc+1,7,7),x%totalcldarea(i), &
         x%meanptop(i), x%meantaucld(i), &
         x%meantb(i), x%meantbclr(i), &
         x%boxtau(i,j), x%boxptop(i,j), &
         x%meanalbedocld(i))
    ! --- Initialise to zero ---
    x%fq_isccp     = 0.0_wp
    x%fq_isccp_ct  = 0.0_wp
    x%totalcldarea = 0.0_wp
    x%meanptop     = 0.0_wp
    x%meantaucld   = 0.0_wp
    x%meantb       = 0.0_wp
    x%meantbclr    = 0.0_wp
    x%boxtau       = 0.0_wp
    x%boxptop      = 0.0_wp
    x%meanalbedocld= 0.0_wp
  END SUBROUTINE CONSTRUCT_COSP_ISCCPSTATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_ISCCPSTATS -----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_ISCCPSTATS(x)
    type(cosp_isccpstats),intent(inout) :: x
    
    deallocate(x%fq_isccp, x%fq_isccp_ct, x%totalcldarea, &
         x%meanptop, x%meantaucld, x%meantb, x%meantbclr, &
         x%boxtau, x%boxptop, x%meanalbedocld)
  END SUBROUTINE FREE_COSP_ISCCPSTATS
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_VGRID ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_VGRID(gbx,Nlvgrid,use_vgrid,cloudsat,x)
    type(cosp_gridbox),intent(in) :: gbx ! Gridbox information
    integer,intent(in) :: Nlvgrid  ! Number of new levels    
    logical,intent(in) :: use_vgrid! Logical flag that controls the output on a different grid
    logical,intent(in) :: cloudsat ! TRUE if a CloudSat like grid (480m) is requested
    type(cosp_vgrid),intent(out) :: x
    
    ! Local variables
    integer :: i
    real(wp):: zstep
    
    x%use_vgrid  = use_vgrid
    x%csat_vgrid = cloudsat
    
    ! Dimensions
    x%Npoints  = gbx%Npoints
    x%Ncolumns = gbx%Ncolumns
    x%Nlevels  = gbx%Nlevels
    
    ! --- Allocate arrays ---
    if (use_vgrid) then
      x%Nlvgrid = Nlvgrid
    else 
      x%Nlvgrid = gbx%Nlevels
    endif
    allocate(x%z(x%Nlvgrid),x%zl(x%Nlvgrid),x%zu(x%Nlvgrid))
    allocate(x%mz(x%Nlevels),x%mzl(x%Nlevels),x%mzu(x%Nlevels))
    
    ! --- Model vertical levels ---
    ! Use height levels of first model gridbox
    x%mz  = gbx%zlev(1,:)
    x%mzl = gbx%zlev_half(1,:)
    x%mzu(1:x%Nlevels-1) = gbx%zlev_half(1,2:x%Nlevels)
    x%mzu(x%Nlevels) = gbx%zlev(1,x%Nlevels) + (gbx%zlev(1,x%Nlevels) - x%mzl(x%Nlevels))
    
    if (use_vgrid) then
      ! --- Initialise to zero ---
      x%z  = 0.0_wp
      x%zl = 0.0_wp
      x%zu = 0.0_wp
      if (cloudsat) then ! --- CloudSat grid requested ---
         zstep = 480.0_wp
      else
         ! Other grid requested. Constant vertical spacing with top at 20 km
         zstep = 20000.0_wp/x%Nlvgrid
      endif
      do i=1,x%Nlvgrid
         x%zl(i) = (i-1)*zstep
         x%zu(i) = i*zstep
      enddo
      x%z = (x%zl + x%zu)/2.0
    else
      x%z  = x%mz
      x%zl = x%mzl
      x%zu = x%mzu
    endif
    
  END SUBROUTINE CONSTRUCT_COSP_VGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_VGRID ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_VGRID(x)
    type(cosp_vgrid),intent(inout) :: x

    deallocate(x%z, x%zl, x%zu, x%mz, x%mzl, x%mzu)
  END SUBROUTINE FREE_COSP_VGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_SGLIDAR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_SGLIDAR(cfg,Npoints,Ncolumns,Nlevels,Nhydro,Nrefl,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    integer,intent(in) :: Nhydro   ! Number of hydrometeors
    integer,intent(in) :: Nrefl    ! Number of parasol reflectances ! parasol
    type(cosp_sglidar),intent(out) :: x
    ! Local variables
    integer :: i,j,k,l,m
    
    ! Allocate minumum storage if simulator not used
    if (cfg%Llidar_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
      l = Nhydro
      m = Nrefl
    else
      i = 1
      j = 1
      k = 1
      l = 1
      m = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    x%Nhydro   = l
    x%Nrefl    = m
    
    ! --- Allocate arrays ---
    allocate(x%beta_mol(i,k), x%beta_tot(i,j,k), &
             x%tau_tot(i,j,k),x%refl(i,j,m), &
             x%temp_tot(i,k),x%betaperp_tot(i,j,k))
    ! --- Initialise to zero ---
    x%beta_mol   = 0.0_wp
    x%beta_tot   = 0.0_wp
    x%tau_tot    = 0.0_wp
    x%refl       = 0.0_wp ! parasol
    x%temp_tot   	= 0.0_wp
    x%betaperp_tot 	= 0.0_wp	
  END SUBROUTINE CONSTRUCT_COSP_SGLIDAR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_SGLIDAR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_SGLIDAR(x)
    type(cosp_sglidar),intent(inout) :: x

    deallocate(x%beta_mol, x%beta_tot, x%tau_tot, x%refl, &
               x%temp_tot, x%betaperp_tot)

  END SUBROUTINE FREE_COSP_SGLIDAR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_SGRADAR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_SGRADAR(cfg,Npoints,Ncolumns,Nlevels,Nhydro,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    integer,intent(in) :: Nhydro   ! Number of hydrometeors
    type(cosp_sgradar),intent(out) :: x
    ! Local variables
    integer :: i,j,k,l
    
    if (cfg%Lradar_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
      l = Nhydro
    else ! Allocate minumum storage if simulator not used
      i = 1
      j = 1
      k = 1
      l = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    x%Nhydro   = l
    
    ! --- Allocate arrays ---
    allocate(x%att_gas(i,k), x%Ze_tot(i,j,k))
    ! --- Initialise to zero ---
    x%att_gas   = 0.0_wp
    x%Ze_tot    = 0.0_wp
    ! The following line give a compilation error on the Met Office NEC
!     call zero_real(x%Z_hydro, x%att_hydro)
!     f90: error(666): cosp_types.f90, line nnn:
!                                        Actual argument corresponding to dummy
!                                        argument of ELEMENTAL subroutine
!                                        "zero_real" with INTENET(OUT) attribute
!                                        is not array.
  END SUBROUTINE CONSTRUCT_COSP_SGRADAR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_SGRADAR ----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_SGRADAR(x)
    type(cosp_sgradar),intent(inout) :: x

    deallocate(x%att_gas, x%Ze_tot)
  END SUBROUTINE FREE_COSP_SGRADAR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------- SUBROUTINE CONSTRUCT_COSP_RADARSTATS ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_RADARSTATS(cfg,Npoints,Ncolumns,Nlevels,Nhydro,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    integer,intent(in) :: Nhydro   ! Number of hydrometeors
    type(cosp_radarstats),intent(out) :: x    
    ! Local variables
    integer :: i,j,k,l
    
    ! Allocate minumum storage if simulator not used
    if (cfg%Lradar_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
      l = Nhydro
    else
      i = 1
      j = 1
      k = 1
      l = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    x%Nhydro   = l
    
    ! --- Allocate arrays ---
    allocate(x%cfad_ze(i,DBZE_BINS,k),x%lidar_only_freq_cloud(i,k))
    allocate(x%radar_lidar_tcc(i))
    ! --- Initialise to zero ---
    x%cfad_ze = 0.0_wp
    x%lidar_only_freq_cloud = 0.0_wp
    x%radar_lidar_tcc = 0.0_wp
  END SUBROUTINE CONSTRUCT_COSP_RADARSTATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_RADARSTATS -------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_RADARSTATS(x)
    type(cosp_radarstats),intent(inout) :: x

    deallocate(x%cfad_ze,x%lidar_only_freq_cloud,x%radar_lidar_tcc)
  END SUBROUTINE FREE_COSP_RADARSTATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------- SUBROUTINE CONSTRUCT_COSP_LIDARSTATS ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_LIDARSTATS(cfg,Npoints,Ncolumns,Nlevels,Nhydro,Nrefl,x)
    type(cosp_config),intent(in) :: cfg ! Configuration options
    integer,intent(in) :: Npoints  ! Number of sampled points
    integer,intent(in) :: Ncolumns ! Number of subgrid columns
    integer,intent(in) :: Nlevels  ! Number of model levels
    integer,intent(in) :: Nhydro   ! Number of hydrometeors
    integer,intent(in) :: Nrefl    ! Number of parasol reflectance
    type(cosp_lidarstats),intent(out) :: x
    ! Local variables
    integer :: i,j,k,l,m
    
    ! Allocate minumum storage if simulator not used
    if (cfg%Llidar_sim) then
      i = Npoints
      j = Ncolumns
      k = Nlevels
      l = Nhydro
      m = Nrefl
    else
      i = 1
      j = 1
      k = 1
      l = 1
      m = 1
    endif
    
    ! Dimensions
    x%Npoints  = i
    x%Ncolumns = j
    x%Nlevels  = k
    x%Nhydro   = l
    x%Nrefl    = m
    
    ! --- Allocate arrays ---
    allocate(x%srbval(SR_BINS),x%cfad_sr(i,SR_BINS,k), &
             x%lidarcld(i,k), x%cldlayer(i,LIDAR_NCAT), x%parasolrefl(i,m))
    allocate(x%lidarcldphase(i,k,6),x%lidarcldtmp(i,LIDAR_NTEMP,5),&
             x%cldlayerphase(i,LIDAR_NCAT,6))
    ! --- Initialise to zero ---
    x%srbval    = 0.0_wp
    x%cfad_sr   = 0.0_wp
    x%lidarcld  = 0.0_wp
    x%cldlayer  = 0.0_wp
    x%parasolrefl  = 0.0_wp
    x%lidarcldphase  = 0.0_wp
    x%cldlayerphase  = 0.0_wp
    x%lidarcldtmp  = 0.0_wp

   END SUBROUTINE CONSTRUCT_COSP_LIDARSTATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------ SUBROUTINE FREE_COSP_LIDARSTATS -------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_LIDARSTATS(x)
    type(cosp_lidarstats),intent(inout) :: x

    deallocate(x%srbval, x%cfad_sr, x%lidarcld, x%cldlayer, x%parasolrefl)
    deallocate(x%cldlayerphase, x%lidarcldtmp, x%lidarcldphase)
  END SUBROUTINE FREE_COSP_LIDARSTATS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_SUBGRID ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_SUBGRID(Npoints,Ncolumns,Nlevels,y)
    integer,intent(in) :: Npoints, & ! Number of gridpoints
                                        Ncolumns, & ! Number of columns
                                        Nlevels   ! Number of levels
    type(cosp_subgrid),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels

    ! --- Allocate arrays ---
    allocate(y%frac_out(Npoints,Ncolumns,Nlevels))
    if (Ncolumns > 1) then
      allocate(y%prec_frac(Npoints,Ncolumns,Nlevels))
    else ! CRM mode, not needed
      allocate(y%prec_frac(1,1,1))
    endif
    ! --- Initialise to zero ---
    y%prec_frac = 0.0_wp
    y%frac_out  = 0.0_wp
    ! The following line gives a compilation error on the Met Office NEC
!     call zero_real(y%mr_hydro)
!     f90: error(666): cosp_types.f90, line nnn:
!                                        Actual argument corresponding to dummy
!                                        argument of ELEMENTAL subroutine
!                                        "zero_real" with INTENET(OUT) attribute
!                                        is not array.

  END SUBROUTINE CONSTRUCT_COSP_SUBGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_SUBGRID -----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_SUBGRID(y)
    type(cosp_subgrid),intent(inout) :: y
    
    ! --- Deallocate arrays ---
    deallocate(y%prec_frac, y%frac_out)
        
  END SUBROUTINE FREE_COSP_SUBGRID

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_SGHYDRO -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_SGHYDRO(Npoints,Ncolumns,Nlevels,Nhydro,y)
    integer,intent(in) :: Npoints, & ! Number of gridpoints
                                        Ncolumns, & ! Number of columns
                                        Nhydro, & ! Number of hydrometeors
                                        Nlevels   ! Number of levels
    type(cosp_sghydro),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    y%Nhydro   = Nhydro

    ! --- Allocate arrays ---
    allocate(y%mr_hydro(Npoints,Ncolumns,Nlevels,Nhydro), &
             y%Reff(Npoints,Ncolumns,Nlevels,Nhydro), &
             y%Np(Npoints,Ncolumns,Nlevels,Nhydro)) ! added by roj with Quickbeam V3
             
    ! --- Initialise to zero ---
    y%mr_hydro = 0.0_wp
    y%Reff     = 0.0_wp
    y%Np       = 0.0_wp                    ! added by roj with Quickbeam V3

  END SUBROUTINE CONSTRUCT_COSP_SGHYDRO

 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_SGHYDRO -----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_SGHYDRO(y)
    type(cosp_sghydro),intent(inout) :: y
    
    ! --- Deallocate arrays ---
    deallocate(y%mr_hydro, y%Reff, y%Np)        ! added by Roj with Quickbeam V3
        
  END SUBROUTINE FREE_COSP_SGHYDRO
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_COSP_GRIDBOX ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_GRIDBOX(radar_freq,surface_radar,use_mie_tables,use_gas_abs,do_ray,melt_lay,k2, &
                                   Npoints,Nlevels,Ncolumns,nclusters,Nhydro,Nprmts_max_hydro,Naero,Nprmts_max_aero, &
                                   lidar_ice_type,isccp_top_height,isccp_top_height_direction,cosp_overlap,isccp_emsfc_lw, &
                                   use_precipitation_fluxes,use_reff, &
                                   ! RTTOV inputs
!CNam:RTTOV                                   Plat,Sat,Inst,Nchan,ZenAng,Ichan,SurfEm,co2,ch4,n2o,co,&
                                   y,load_LUT)

!CNam    double precision,intent(in) :: time ! Time since start of run [days] 
!CNam    double precision,intent(in) :: time_bnds(2) ! Time boundaries

    real(wp),intent(in) :: radar_freq, & ! Radar frequency [GHz]
				k2        ! |K|^2, -1=use frequency dependent default
    integer, intent(in) :: &
        surface_radar, &  ! surface=1,spaceborne=0
        use_mie_tables, & ! use a precomputed lookup table? yes=1,no=0,2=use first column everywhere
        use_gas_abs, &    ! include gaseous absorption? yes=1,no=0
        do_ray, &         ! calculate/output Rayleigh refl=1, not=0
        melt_lay          ! melting layer model off=0, on=1
    integer,intent(in) :: Npoints   ! Number of gridpoints
    integer,intent(in) :: Nlevels   ! Number of levels
    integer,intent(in) :: Ncolumns  ! Number of columns
    integer,intent(in) :: nclusters  ! Number of formation pathway clusters
    integer,intent(in) :: Nhydro    ! Number of hydrometeors
    integer,intent(in) :: Nprmts_max_hydro    ! Max number of parameters for hydrometeor size distributions
    integer,intent(in) :: Naero    ! Number of aerosol species
    integer,intent(in) :: Nprmts_max_aero    ! Max number of parameters for aerosol size distributions
!CNam    integer,intent(in) :: Npoints_it   ! Number of gridpoints processed in one iteration
    integer,intent(in) :: lidar_ice_type ! Ice particle shape in lidar calculations (0=ice-spheres ; 1=ice-non-spherical)
    integer,intent(in) :: isccp_top_height
    integer,intent(in) :: isccp_top_height_direction
    integer,intent(in) :: cosp_overlap
    real(wp),intent(in):: isccp_emsfc_lw
    logical,intent(in) :: use_precipitation_fluxes,use_reff
!    integer,intent(in) :: Plat    !CNam: RTTOV
!    integer,intent(in) :: Sat
!    integer,intent(in) :: Inst
!    integer,intent(in) :: Nchan
!    integer,intent(in) :: Ichan(Nchan)
!    real(wp), intent(in):: SurfEm(Nchan)
!    real(wp), intent(in):: ZenAng
!    real(wp), intent(in):: co2,ch4,n2o,co
    type(cosp_gridbox), intent(out) :: y
    logical,intent(in),optional :: load_LUT


    ! local variables
    character*240 :: LUT_file_name
    logical :: local_load_LUT

    if (present(load_LUT)) then
      local_load_LUT = load_LUT
    else
      local_load_LUT = RADAR_SIM_LOAD_scale_LUTs_flag
    endif

    ! Dimensions and scalars
    y%radar_freq       = radar_freq
    y%surface_radar    = surface_radar
    y%use_mie_tables   = use_mie_tables
    y%use_gas_abs      = use_gas_abs
    y%do_ray           = do_ray
    y%melt_lay         = melt_lay
    y%k2               = k2
    y%Npoints          = Npoints
    y%Nlevels          = Nlevels
    y%nclusters        = nclusters
    y%Ncolumns         = Ncolumns
    y%Nhydro           = Nhydro
    y%Nprmts_max_hydro = Nprmts_max_hydro
    y%Naero            = Naero
    y%Nprmts_max_aero  = Nprmts_max_aero
!CNam    y%Npoints_it       = Npoints_it
    y%lidar_ice_type   = lidar_ice_type
    y%isccp_top_height = isccp_top_height
    y%isccp_top_height_direction = isccp_top_height_direction
    y%cosp_overlap    = cosp_overlap
    y%isccp_emsfc_lw   = isccp_emsfc_lw
    y%use_precipitation_fluxes = use_precipitation_fluxes
    y%use_reff = use_reff
    
!CNam    y%time = time
!CNam    y%time_bnds = time_bnds
    
    ! RTTOV parameters !CNam: RTTOV
!    y%Plat   = Plat 
!    y%Sat    = Sat  
!    y%Inst   = Inst
!    y%Nchan  = Nchan
!    y%ZenAng = ZenAng
!    y%co2    = co2
!    y%ch4    = ch4
!    y%n2o    = n2o
!    y%co     = co

    ! --- Allocate arrays ---
    ! Gridbox information (Npoints,Nlevels)
    allocate(y%zlev(Npoints,Nlevels), y%zlev_half(Npoints,Nlevels), y%dlev(Npoints,Nlevels), &
             y%p(Npoints,Nlevels), y%ph(Npoints,Nlevels+1), y%T(Npoints,Nlevels), &    !CNam:orig y%ph(Npoints,Nlevels)
             y%q(Npoints,Nlevels), y%sh(Npoints,Nlevels), y%labels(Npoints,Nlevels), &
             y%dtau_s(Npoints,Nlevels), y%dtau_c(Npoints,Nlevels), &
             y%dem_s(Npoints,Nlevels), y%dem_c(Npoints,Nlevels), &
             y%tca(Npoints,Nlevels), y%cca(Npoints,Nlevels), &
             y%rain_ls(Npoints,Nlevels), y%rain_cv(Npoints,Nlevels), y%grpl_ls(Npoints,Nlevels), &
             y%snow_ls(Npoints,Nlevels), y%snow_cv(Npoints,Nlevels),y%mr_ozone(Npoints,Nlevels))
             
             
    ! Surface information and geolocation (Npoints)
    allocate(y%toffset(Npoints), y%longitude(Npoints),y%latitude(Npoints),y%psfc(Npoints), y%land(Npoints), &
             y%sunlit(Npoints),y%skt(Npoints),y%u_wind(Npoints),y%v_wind(Npoints))
    ! Hydrometeors concentration and distribution parameters
    allocate(y%mr_hydro(Npoints,Nlevels,Nhydro), &
             y%dist_prmts_hydro(Nprmts_max_hydro,Nhydro), &
             y%Reff(Npoints,Nlevels,Nhydro), &
             y%Np(Npoints,Nlevels,Nhydro))      ! added by Roj with Quickbeam V3
    ! Aerosols concentration and distribution parameters
    allocate(y%conc_aero(Npoints,Nlevels,Naero), y%dist_type_aero(Naero), &
             y%dist_prmts_aero(Npoints,Nlevels,Nprmts_max_aero,Naero))
    
    ! RTTOV channels and sfc. emissivity
!    allocate(y%ichan(Nchan),y%surfem(Nchan)) !CNam: RTTOV
    
    ! RTTOV parameters
!    y%ichan   =  ichan	 !CNam: RTTOV
!    y%surfem  =  surfem !CNam: RTTOV
    
    ! --- Initialise to zero ---
    y%zlev      = 0.0_wp
    y%zlev_half = 0.0_wp
    y%dlev      = 0.0_wp
    y%p         = 0.0_wp
    y%ph        = 0.0_wp
    y%T         = 0.0_wp
    y%q         = 0.0_wp
    y%sh        = 0.0_wp
    y%dtau_s    = 0.0_wp
    y%dtau_c    = 0.0_wp
    y%dem_s     = 0.0_wp
    y%dem_c     = 0.0_wp
    y%tca       = 0.0_wp
    y%cca       = 0.0_wp
    y%rain_ls   = 0.0_wp
    y%rain_cv   = 0.0_wp
    y%grpl_ls   = 0.0_wp
    y%snow_ls   = 0.0_wp
    y%snow_cv   = 0.0_wp
    y%Reff      = 0.0_wp
    y%Np        = 0.0_wp ! added by Roj with Quickbeam V3
    y%mr_ozone  = 0.0_wp
    y%u_wind    = 0.0_wp
    y%v_wind    = 0.0_wp

    
    ! (Npoints)
    y%toffset = 0.0_wp
    y%longitude = 0.0_wp
    y%latitude = 0.0_wp
    y%psfc = 0.0_wp
    y%land = 0.0_wp
    y%sunlit = 0.0_wp
    y%skt = 0.0_wp
    ! (Npoints,Nlevels,Nhydro)
!     y%fr_hydro = 0.0_wp
    y%mr_hydro = 0.0_wp
    ! Others
    y%dist_prmts_hydro = 0.0_wp ! (Nprmts_max_hydro,Nhydro)
    y%conc_aero        = 0.0_wp ! (Npoints,Nlevels,Naero)
    y%dist_type_aero   = 0   ! (Naero)
    y%dist_prmts_aero  = 0.0_wp ! (Npoints,Nlevels,Nprmts_max_aero,Naero)


    ! NOTE: This location use to contain initialization of some radar simulator variables
    ! this initialization (including use of the variable "dist_prmts_hydro" - now obselete) 
    ! has been unified in the quickbeam v3 subroutine "cosp_radar_simulator_init".   Roj, June 2010

    ! --- Initialize the distributional parameters for hydrometeors in radar simulator

    !CNam: write(*,*) 'RADAR_SIM microphysics scheme is set to: ', &
    !            trim(RADAR_SIM_MICROPHYSICS_SCHEME_NAME)


    if(y%Nhydro.ne.N_HYDRO) then

        write(*,*) 'Number of hydrometeor input to subroutine', &
               ' CONSTRUCT_COSP_GRIDBOX does not match value', &
               ' specified in cosp_constants.f90!'
        write(*,*) 
    endif

    ! NOTE: SAVE_scale_LUTs_flag is hard codded as .false. here 
    ! so that radar simulator will NOT update LUT each time it 
    ! is called, but rather will update when "Free_COSP_GRIDBOX" is called!
    ! Roj, June 2010

    LUT_file_name = trim(RADAR_SIM_LUT_DIRECTORY) // &
                trim(RADAR_SIM_MICROPHYSICS_SCHEME_NAME)

    call cosp_radar_simulator_init(radar_freq,k2, &
                      use_gas_abs,do_ray,R_UNDEF, &
                      y%Nhydro, &
                      HCLASS_TYPE,HCLASS_PHASE, &
                      HCLASS_DMIN,HCLASS_DMAX, &
                      HCLASS_APM,HCLASS_BPM,HCLASS_RHO, &
                      HCLASS_P1,HCLASS_P2,HCLASS_P3, &
                      local_load_LUT,    &
                      .false., &
                      LUT_file_name, &
                      y%hp)

END SUBROUTINE CONSTRUCT_COSP_GRIDBOX


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE FREE_COSP_GRIDBOX -----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_GRIDBOX(y,dglobal,save_LUT)

    use mo_cosp_scale_LUTs_io

    type(cosp_gridbox),intent(inout) :: y
    logical,intent(in),optional :: dglobal
    logical,intent(in),optional :: save_LUT

    logical :: local_save_LUT

    if (present(save_LUT)) then
      local_save_LUT = save_LUT
    else
      local_save_LUT = RADAR_SIM_UPDATE_scale_LUTs_flag
    endif

    ! save any updates to radar simulator LUT
    if (local_save_LUT) call save_scale_LUTs(y%hp)

    deallocate(y%zlev, y%zlev_half, y%dlev, y%p, y%ph, y%T, y%q, &
               y%sh, y%dtau_s, y%dtau_c, y%dem_s, y%dem_c, y%labels, &
               y%toffset, y%longitude,y%latitude,y%psfc, y%land, y%tca, y%cca, &
               y%mr_hydro, y%dist_prmts_hydro, &
               y%conc_aero, y%dist_type_aero, y%dist_prmts_aero, &
               y%rain_ls, y%rain_cv, y%snow_ls, y%snow_cv, y%grpl_ls, &
               y%sunlit, y%skt, y%Reff,y%Np, &
!CNam: rttov               y%ichan,y%surfem, &
               y%mr_ozone,y%u_wind,y%v_wind)

  END SUBROUTINE FREE_COSP_GRIDBOX

END MODULE MO_COSP_TYPES
