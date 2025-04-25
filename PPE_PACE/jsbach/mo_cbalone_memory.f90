!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cbalone_memory
  USE mo_kind,             ONLY: dp
  USE mo_exception,        ONLY: finish
  USE mo_jsbach_grid,      ONLY: grid_type
  USE mo_cbal_bethy,       ONLY: cbalance_type

  implicit none

  private

  public :: cbal_offline_type
  public :: nbal_offline_type
  public :: dynveg_clim_type
  public :: initializePools       !! reads initial carbon pools from filepack
  public :: initializeYassoPools  !! reads initial pools for Yasso
  public :: initializeNPools      !! reads initial nitrogen pools from filepack
  public :: initSurface           !! initializes the structure "surface" from file
  public :: initSoil              !! initializes the structure "soil" from file
  public :: initCbalance          !! initializes the structure "cbalance" (except carbon pools)
  public :: initNbalance          !! initializes the structure "nbalance" (except nitrogen pools)
  public :: check_err

  TYPE cbal_offline_type
     real(dp),pointer :: LAI_yDayMean(:,:,:)            !! mean value of LAI yesterday (from  LAI_sum())
     real(dp),pointer :: LAI_previousDayMean(:,:)       !! mean value of LAI the day before yesterday
     real(dp),pointer :: NPP_yDayMean(:,:,:)            !! mean value of potential NPP-Rate yesterday, i.e. before N-limitation 
                                                        !!    (from NPP_sum()) [mol(CO2)/(m^2(canopy) s)]
     real(dp),pointer :: topSoilTemp_yDayMean(:,:,:)    !! mean value of upper layer soil temperature yesterday 
                                                        !!    (from topSoilTemp_sum()) [K]
     real(dp),pointer :: alpha_yDayMean(:,:,:)          !! mean value of water stress coefficient alpha yesterday (from alpha_sum())
     real(dp),pointer :: pseudo_temp_yDay(:,:,:)        !! 15 day mean air temperature for yasso
     real(dp),pointer :: pseudo_precip_yDay(:,:,:)      !! 15 day mean precipitation for yasso

     real(dp),pointer :: box_Cpools_total(:)            !! Sum of all carbon pools [mol(C)/m^2(grid box)]
     real(dp),pointer :: box_test_Ccons(:)              !! Difference between all C before and after a time step: should be zero if
                                                        !!    C is conserved [mol(C)/m^2(grid box)]
     real(dp),pointer :: box_CO2_flux_2_atm(:)          !! Output period mean CO2 flux to atm (kg/(m^2(grid box)s))
     real(dp),pointer :: avg_Cpool_green(:,:)           !! As Cpool_green() but time averaged and relative to grid box area 
                                                        !!     [mol(C)/m^2(grid box)] 
     real(dp),pointer :: avg_Cpool_reserve(:,:)         !! As Cpool_reserve() but time averaged and relative to grid box area 
                                                        !!    [mol(C)/m^2(grid box)]    
     real(dp),pointer :: avg_Cpool_woods(:,:)           !! As Cpool_woods() but time averaged and relative to grid box area 
                                                        !!    [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_Cpool_litter_green_ag(:,:) !! As Cpool_litter_green() but time averaged and relative to grid box area 
                                                        !!    [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_Cpool_litter_green_bg(:,:) !! As Cpool_litter_green_bg() but time averaged and relative to grid box
                                                        !!    area [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_Cpool_litter_wood_ag(:,:)  !! As Cpool_litter_wood_ag() but time averaged and relative to grid box area
                                                        !!    [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_Cpool_litter_wood_bg(:,:)  !! As Cpool_litter_wood_bg() but time averaged and relative to grid box area
                                                        !!    [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_Cpool_slow(:,:)            !! As Cpool_slow() but time averaged and relative to grid box area 
                                                        !!    [mol(C)/m^2(grid box)
                                                        !! Yasso above ground litter pools  (time averages)
                                                        !! Size class 1: green litter
     real(dp),pointer :: avg_YCpool_acid_ag1(:,:)        !!     for acid soluble litter     [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_water_ag1(:,:)       !!     for water soluble litter    [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_ethanol_ag1(:,:)     !!     for ethanol soluble litter  [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_nonsoluble_ag1(:,:)  !!     for non-soluble litter      [mol(C)/m^2(grid box)]
                                                        !! Yasso below ground litter pools  (time averages)
     real(dp),pointer :: avg_YCpool_acid_bg1(:,:)        !!     for acid soluble litter     [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_water_bg1(:,:)       !!     for water soluble litter    [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_ethanol_bg1(:,:)     !!     for ethanol soluble litter  [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_nonsoluble_bg1(:,:)  !!     for non-soluble litter      [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_humus_1(:,:)          !!     for slow carbon compartment [mol(C)/m²(grid box)]
                                                        !! Size class 2: woody litter
     real(dp),pointer :: avg_YCpool_acid_ag2(:,:)        !!     for acid soluble litter     [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_water_ag2(:,:)       !!     for water soluble litter    [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_ethanol_ag2(:,:)     !!     for ethanol soluble litter  [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_nonsoluble_ag2(:,:)  !!     for non-soluble litter      [mol(C)/m^2(grid box)]
                                                        !! Yasso below ground litter pools  (time averages)
     real(dp),pointer :: avg_YCpool_acid_bg2(:,:)        !!     for acid soluble litter     [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_water_bg2(:,:)       !!     for water soluble litter    [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_ethanol_bg2(:,:)     !!     for ethanol soluble litter  [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_nonsoluble_bg2(:,:)  !!     for non-soluble litter      [mol(C)/m^2(grid box)]
     real(dp),pointer :: avg_YCpool_humus_2(:,:)          !!     for slow carbon compartment [mol(C)/m²(grid box)]

     real(dp),pointer :: avg_soil_respiration(:,:)    !! Average soil respiration rate [mol(CO2)/m^2(grid box) s]
     real(dp),pointer :: avg_soil_respiration_pot(:,:)   !!Average soil respiration rate without N limitation !dk new
     real(dp),pointer :: avg_NPP_yDayMean(:,:)        !! Average NPP rate [mol(CO2)/m^2(grid box) s]
     real(dp),pointer :: avg_NPP_flux_correction(:,:) !! Average flux correction for NPP [mol(CO2)/m^2(grid box) s]
     real(dp),pointer :: avg_excess_NPP(:,:)          !! That part of NPP that could not be stored in the carbon pools of the living
                                                      !!   plant but had to be dropped into the fast pool [mol(CO2)/m^2(grid box) s]
     real(dp),pointer :: avg_root_exudates(:,:)
     real(dp),pointer :: avg_Cflux_herbivory(:,:)
     real(dp),pointer :: avg_Cflux_herbivory_LG(:,:)
     real(dp),pointer :: avg_Cflux_herbivory_2_atm(:,:)
     real(dp),pointer :: avg_box_NEP(:)               !! Net ecosystem productivity [kg(CO2)/m^2(grid box) s]
     real(dp),pointer :: NPP_act_yDayMean(:,:)        !! mean value of actual NPP-Rate yesterday , i.e. after N-limitation 
                                                      !!    [mol(CO2)/(m^2(canopy) s)]
     real(dp),pointer :: avg_NPP_act_yDayMean(:,:)       
     real(dp),pointer :: testCconserv(:,:)            !! Array for testing of Carbon conservation. Should be zero if Carbon is 
                                                      !!    conserved.
     ! Anthropogenic C pools
     real(dp),pointer :: avg_Cpool_onSite              (:)
     real(dp),pointer :: avg_Cpool_paper               (:)
     real(dp),pointer :: avg_Cpool_paper_harvest       (:)
     real(dp),pointer :: avg_Cpool_construction        (:)
     real(dp),pointer :: avg_Cpool_construction_harvest(:)
  end TYPE cbal_offline_type

  TYPE Nbal_offline_type
     real(dp),pointer :: avg_Npool_green(:,:)           !! As Npool_green() but time averaged and relative to grid box area 
                                                        !!    [mol(N)/m^2(grid box)] 
     real(dp),pointer :: avg_Npool_mobile(:,:)          !! As Npool_mobile() but time averaged and relative to grid box area 
                                                        !!    [mol(N)/m^2(grid box)]    
     real(dp),pointer :: avg_Npool_woods(:,:)           !! As Npool_woods() but time averaged and relative to grid box area 
                                                        !!    [mol(N)/m^2(grid box)]
     real(dp),pointer :: avg_Npool_litter_green_ag(:,:) !! As Npool_litter_green() but time averaged and relative to grid box area
                                                        !!    [mol(N)/m^2(grid box)]
     real(dp),pointer :: avg_Npool_litter_green_bg(:,:) !! As Npool_litter_green_bg() but time averaged and relative to grid box 
                                                        !!    area [mol(N)/m^2(grid box)]
     real(dp),pointer :: avg_Npool_litter_wood_ag(:,:)  !! As Npool_litter_wood_ag() but time averaged and relative to grid box area
                                                        !!    [mol(N)/m^2(grid box)]
     real(dp),pointer :: avg_Npool_litter_wood_bg(:,:)  !! As Npool_litter_wood_bg() but time averaged and relative to grid box area
                                                        !!    [mol(N)/m^2(grid box)]
     real(dp),pointer :: avg_Npool_slow(:,:)            !! As Npool_slow() but time averaged and relative to grid box area 
                                                        !!    [mol(N)/m^2(grid box)
     real(dp),pointer :: avg_SMINN_pool(:,:)

     real(dp),pointer :: box_Npools_total(:)            !! Sum of all nitrogen pools [mol(N)/m^2(grid box)]
     real(dp),pointer :: box_test_Ncons(:)              !! Difference between all N before and after a time step: should be zero if
                                                        !!    N is conserved [mol(N)/m^2(grid box)]
     real(dp),pointer :: avg_Nplant_demand(:,:)         !! As Nplant_demand() but time averaged and relative to grid box area 
                                                        !!    [mol(N)/m^2(grid box)] 
     real(dp),pointer :: avg_Nsoil_demand(:,:)          !! As Nsoil_demand() but time averaged and relative to grid box area 
                                                        !!    [mol(N)/m^2(grid box)] 
     real(dp),pointer :: avg_Ntotal_demand(:,:)         !! As Ntotal_demand() but time averaged and relative to grid box area 
                                                        !!    [mol(N)/m^2(grid box)] 
     real(dp),pointer :: avg_SMINN_herbivory(:,:)       !! SMINN gain from herbivory faeces and dungs [mol(N)/m^2(grid box)]      
     real(dp),pointer :: avg_nfix_to_sminn(:,:)  
     real(dp),pointer :: avg_sminn_to_denit(:,:)
     real(dp),pointer :: avg_sminn_leach(:,:)
     real(dp),pointer :: avg_nfert_to_sminn(:,:)        !! constatnt rate of N fertiliser application rate 
     real(dp),pointer :: avg_N2O_emissions_ExternalN(:,:)  !! N2O emissions [mol(N)/m^2(grid box)] due to deposition & fixation
     real(dp),pointer :: avg_N2O_emissions_mineraliz(:,:)  !! N2O emissions [mol(N)/m^2(grid box)] due to N -mineralization
     real(dp),pointer :: avg_NetEcosyst_N_flux(:,:)     !! Total balance of N gains and losses (positive for ecosystem gain) 
                                                        !!    [mol(N)/m^2(canopy)s]
     real(dp),pointer :: avg_N2O_emissions_nfert(:,:)   !! N2O emissions [mol(N)/m^2(grid box)] due to N-fertilizers use
     real(dp),pointer :: avg_N2O_emissions_grazing(:,:) !! N2O emissions [mol(N)/m^2(grid box)] due to N input by herbivores
     REAL(dp),POINTER :: minNflux_litter_green_ag(:,:)
     REAL(dp),POINTER :: minNflux_litter_green_bg(:,:)
     REAL(dp),POINTER :: minNflux_litter_wood_ag(:,:)
     REAL(dp),POINTER :: minNflux_litter_wood_bg(:,:)
     REAL(dp),POINTER :: minNflux_slow(:,:)  
     REAL(dp),POINTER :: avg_minNflux_litter_green_ag(:,:)
     REAL(dp),POINTER :: avg_minNflux_litter_green_bg(:,:)
     REAL(dp),POINTER :: avg_minNflux_litter_wood_ag(:,:)
     REAL(dp),POINTER :: avg_minNflux_litter_wood_bg(:,:)
     REAL(dp),POINTER :: avg_minNflux_slow(:,:)

     real(dp),pointer :: test_Ngreen(:,:)
     real(dp),pointer :: test_Nwoods(:,:)
     real(dp),pointer :: test_Nlitter_green_ag(:,:)
     real(dp),pointer :: test_Nlitter_green_bg(:,:)
     real(dp),pointer :: test_Nlitter_wood_ag(:,:)
     real(dp),pointer :: test_Nlitter_wood_bg(:,:)
     real(dp),pointer :: test_Nslow(:,:)
     real(dp),pointer :: testNconserv(:,:)     !! Array for testing of Nitrogen conservation. Should be zero if Carbon is conserved.
     real(dp),pointer :: Runoff_yDayMean(:,:,:)     !! mean value of runoff yesterday (from runoff_sum())
  END TYPE Nbal_offline_type

  TYPE dynveg_clim_type                            !! input for dynamic vegetation (dimensions including time)
     REAL(dp), POINTER :: ave_npp5(:,:,:) 
     REAL(dp), POINTER :: bio_exist(:,:,:)
     REAL(dp), POINTER :: rel_hum_air(:,:)
     REAL(dp), POINTER :: max_wind10(:,:)
     REAL(dp), POINTER :: prev_day_max_wind10(:,:)
     REAL(dp), POINTER :: prev_day_mean_wind10(:,:)
     REAL(dp), POINTER :: prev_day_temp_min(:,:)
     REAL(dp), POINTER :: prev_day_temp_max(:,:)
!     REAL(dp), POINTER :: dew_point_temp(:,:)
     REAL(dp), POINTER :: prev_day_precip_mean(:,:)
     REAL(dp), POINTER :: prev_day_mean_vol_moist(:,:)
  END TYPE dynveg_clim_type

  !! Parameters

  logical,parameter  ::  debug = .false.  !! If this is set to .true. debug information is printed out

contains

  ! --- initializePools() --------------------------------------------------------------------------------
  !
  ! Initializes carbon pools from a file (any output file from Cbalone can be taken as initialization file)
  ! If carbon pools are read in from a file, they must be in LatLon format, packed format is currently not
  ! supported.
  !
  SUBROUTINE initializePools(grid, cbalance, cbalance_diag, surface, filename)
    USE mo_cbal_cpools,  ONLY: printCpoolParameters
    USE mo_land_surface, ONLY: land_surface_type

    include 'netcdf.inc'
    type(grid_type),        intent(in)   :: grid
    TYPE(cbal_offline_type),  INTENT(inout) :: cbalance_diag
    TYPE(cbalance_type),      INTENT(inout) :: cbalance
    type(land_surface_type),intent(inout):: surface
    character(len=*),OPTIONAL,intent(in) :: filename

    !! --- internal variables

    real(dp),allocatable,dimension(:,:,:,:)  :: read_Cpools
    ! error status return
    integer  ::  iret
    ! netCDF id
    integer  ::  ncid
    ! dimension
    integer  ::  ndim,nvar,nattr,nunlimdim
    integer,allocatable,dimension(:)  ::  ndimlen
    character(len=50),allocatable,dimension(:)  ::  fdimnam
    ! variables
    integer ::  nin(100)
    integer,allocatable ,dimension(:)  ::  nvartyp
    integer,allocatable,dimension(:)   ::  nvardim
    integer,allocatable,dimension(:,:) ::  nvardimid
    integer,allocatable,dimension(:)   ::  nvaratt
    character(len=128),allocatable,dimension(:)  ::  fvarnam
    integer :: noTimeSteps
    logical :: time_dimension
    ! indices
    integer  ::  i,ii,itile

    !! Printout parameters 

    call printCpoolParameters

    !! --- initialize the pools
    
    cbalance_diag%LAI_previousDayMean(:,:) = 0.0_dp
    cbalance%Cpool_litter_green_ag   (:,:) = 0.0_dp
    cbalance%Cpool_litter_green_bg   (:,:) = 0.0_dp
    cbalance%Cpool_litter_wood_ag    (:,:) = 0.0_dp
    cbalance%Cpool_litter_wood_bg    (:,:) = 0.0_dp
    cbalance%Cpool_slow              (:,:) = 0.0_dp
    if( .not. present(filename) ) then            !! Initialize from scratch
       cbalance%Cpool_green             (:,:) = 0.0_dp
       cbalance%Cpool_woods             (:,:) = 0.0_dp
       cbalance%Cpool_reserve           (:,:) = 0.0_dp
       IF (ASSOCIATED(cbalance%Cpool_onSite              )) cbalance%Cpool_onSite              (:) = 0._dp
       IF (ASSOCIATED(cbalance%Cpool_paper               )) cbalance%Cpool_paper               (:) = 0._dp
       IF (ASSOCIATED(cbalance%Cpool_construction        )) cbalance%Cpool_construction        (:) = 0._dp
       IF (ASSOCIATED(cbalance%Cpool_paper_harvest       )) cbalance%Cpool_paper_harvest       (:) = 0._dp
       IF (ASSOCIATED(cbalance%Cpool_construction_harvest)) cbalance%Cpool_construction_harvest(:) = 0._dp
    else                                          !! Initialize from file
       !! open the file
       write (*,*) 'initializePools(): read Cpools from: ',filename
       iret = nf_open(filename,nf_nowrite,ncid)
       call check_err(iret)
       ! check what is in the input-file
       iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
       call check_err(iret)
       if (debug) write (*,*) 'initializePools(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim
       ! get the dimension name and length and find out whether there is a time dimension
       allocate (ndimlen(ndim))
       allocate (fdimnam(ndim))
       time_dimension=.false.
       do i=1,ndim
          iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
          call check_err(iret)
          if (debug) write (*,*) 'initializePools(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
          if (TRIM(fdimnam(i)) == 'time') time_dimension=.true.
       end do
       ! get variable names, types and shapes
       allocate (fvarnam(nvar))
       allocate (nvartyp(nvar))
       allocate (nvardim(nvar))
       allocate (nvardimid(nvar,100))
       allocate (nvaratt(nvar))
       do i=1,nvar
          iret = nf_inq_var(ncid,i,fvarnam(i),nvartyp(i),nvardim(i),nin,nvaratt(i))
          call check_err(iret)
          do ii=1,nvardim(i)
             nvardimid(i,ii)=nin(ii)
          end do
          if (debug) write (*,*) 'initializePools(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
               i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
       end do
       ! get data
       do i=1,nvar
          if (nvardim(i) >= 3) then  ! carbon pools have at least 3 dimensions (lon, lat, tiles)
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                call finish ('initializePools', 'dimensions of '//TRIM(fvarnam(i))//' wrong')
             endif
             if (time_dimension) then
                noTimeSteps = ndimlen(nvardimid(i,nvardim(i)))
             else
                noTimeSteps = 1
             end if
             allocate (read_Cpools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),noTimeSteps))
          endif

          if (fvarnam(i) == 'Cpool_green') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%Cpool_green(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read Cpool_green"
          else if(fvarnam(i) == 'Cpool_reserve') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%Cpool_reserve(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read Cpool_reserve"
          else if(fvarnam(i) == 'Cpool_woods') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%Cpool_woods(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read Cpool_woods"
          else if(fvarnam(i) == 'Cpool_litter_green_ag') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%Cpool_litter_green_ag(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read Cpool_litter_green_ag"

          else if(fvarnam(i) == 'Cpool_litter_green_bg') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%Cpool_litter_green_bg(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read Cpool_litter_green_bg"

          else if(fvarnam(i) == 'Cpool_litter_wood_ag') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%Cpool_litter_wood_ag(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read Cpool_litter_wood_ag"

          else if(fvarnam(i) == 'Cpool_litter_wood_bg') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%Cpool_litter_wood_bg(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read Cpool_litter_wood_bg"

          else if(fvarnam(i) == 'Cpool_slow') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%Cpool_slow(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read Cpool_slow"
          else if(fvarnam(i) == 'LAI_previousDayMean') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance_diag%LAI_previousDayMean(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read LAI_previousDayMean"
          else if (fvarnam(i) == 'cover_fract') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                surface%cover_fract(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read cover_fract"
          else if (fvarnam(i) == 'cover_fract_pot') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                surface%cover_fract_pot(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializePools():   read cover_fract_pot"                      
          else if (fvarnam(i) == 'LCC_C_onSite') then
             if (associated(cbalance%Cpool_onSite)) then
               iret = nf_get_var_double(ncid,i,read_Cpools)
               call check_err(iret)
               cbalance%Cpool_onSite(:) = PACK(read_Cpools(:,:,noTimeSteps,1),MASK=grid%mask)
               write(*,*) "initializePools():   read LCC_C_onSite"
             endif
          else if (fvarnam(i) == 'LCC_C_paper') then
             if (associated(cbalance%Cpool_paper)) then
               iret = nf_get_var_double(ncid,i,read_Cpools)
               call check_err(iret)
               cbalance%Cpool_paper(:) = PACK(read_Cpools(:,:,noTimeSteps,1),MASK=grid%mask)
               write(*,*) "initializePools():   read LCC_C_paper"
             endif
          else if (fvarnam(i) == 'LCC_C_construction') then
             if (associated(cbalance%Cpool_construction)) then
               iret = nf_get_var_double(ncid,i,read_Cpools)
               call check_err(iret)
               cbalance%Cpool_construction(:) = PACK(read_Cpools(:,:,noTimeSteps,1),MASK=grid%mask)
               write(*,*) "initializePools():   read LCC_C_construction"
             endif
          else if (fvarnam(i) == 'LCC_C_paper_harvest') then
             if (associated(cbalance%Cpool_paper_harvest)) then
               iret = nf_get_var_double(ncid,i,read_Cpools)
               call check_err(iret)
               cbalance%Cpool_paper_harvest(:) = PACK(read_Cpools(:,:,noTimeSteps,1),MASK=grid%mask)
               write(*,*) "initializePools():   read LCC_C_paper_harvest"
             endif
          else if (fvarnam(i) == 'LCC_C_construction_harvest') then
             if (associated(cbalance%Cpool_construction_harvest)) then
               iret = nf_get_var_double(ncid,i,read_Cpools)
               call check_err(iret)
               cbalance%Cpool_construction_harvest(:) = PACK(read_Cpools(:,:,noTimeSteps,1),MASK=grid%mask)
               write(*,*) "initializePools():   read LCC_C_construction_harvest"
             endif
          else if (fvarnam(i) == 'frac_litter_wood_new') then
             if (associated(cbalance%frac_litter_wood_new)) then
               iret = nf_get_var_double(ncid,i,read_Cpools)
               call check_err(iret)
               do itile = 1,surface%ntiles
                  cbalance%frac_litter_wood_new(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
               end do
               write(*,*) "initializePools():   read frac_litter_wood_new"
             endif
          end if
          if (allocated(read_Cpools)) deallocate(read_Cpools)
       end do
       write(*,*) "initializePools(): Carbon pools read in. Took values from last time step ",noTimeSteps," found in file."
       
       ! close the file
       iret = nf_close(NCID)
       call check_err(iret)
       ! deallocate
       deallocate(fdimnam,ndimlen)
       deallocate(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
    end if
    
  end subroutine initializePools


  ! --- initialize Nitrogen Pools() --------------------------------------------------------------------------------
  !
  ! Initializes 9 nitrogen pools from a file (any output file from Nbalone can be taken as initialization file)
  ! If nitrogen pools are read in from a file, they must be in LatLon format, packed format is currently not
  ! supported.
  !
  subroutine initializeNPools(grid, ntiles, nbalance, filename, cbalance)
    USE mo_cbal_cpools,     ONLY: printNpoolParameters
    USE mo_cbal_bethy,      ONLY: nbalance_type
    USE mo_cbal_parameters, ONLY: cn_green, cn_woods, cn_litter_green, cn_litter_wood, cn_slow

    include 'netcdf.inc'
    TYPE(grid_type),                   INTENT(in)    :: grid
    INTEGER,                           INTENT(in)    :: ntiles
    TYPE(nbalance_type),               INTENT(inout) :: nbalance
    character(len=*),        OPTIONAL, INTENT(in)    :: filename
    TYPE(cbalance_type),     OPTIONAL, INTENT(inout) :: cbalance

    !! --- internal variables

    real,allocatable,dimension(:,:,:,:)  :: read_Npools
    ! error status return
    integer  ::  iret
    ! netCDF id
    integer  ::  ncid
    ! dimension
    integer  ::  ndim,nvar,nattr,nunlimdim
    integer,allocatable,dimension(:)  ::  ndimlen
    character(len=50),allocatable,dimension(:)  ::  fdimnam
    ! variables
    integer ::  nin(100)
    integer,allocatable ,dimension(:)  ::  nvartyp
    integer,allocatable,dimension(:)   ::  nvardim
    integer,allocatable,dimension(:,:) ::  nvardimid
    integer,allocatable,dimension(:)   ::  nvaratt
    character(len=128),allocatable,dimension(:)  ::  fvarnam
    integer :: noTimeSteps
    ! indices
    integer  ::  i,ii,itile

    !! Printout parameters 

    call printNpoolParameters

    !! --- initialize the nitrogen  pools
    
    if( .not. present(filename) ) then            !! Initialize from scratch
       nbalance%Npool_green(:,:)           = 0.0_dp
       nbalance%Npool_woods(:,:)           = 0.0_dp
       nbalance%Npool_mobile(:,:)          = 0.0_dp
       nbalance%Npool_litter_green_ag(:,:) = 0.0_dp
       nbalance%Npool_litter_green_bg(:,:) = 0.0_dp
       nbalance%Npool_litter_wood_ag(:,:)  = 0.0_dp
       nbalance%Npool_litter_wood_bg(:,:)  = 0.0_dp
       nbalance%Npool_slow(:,:)            = 0.0_dp
       nbalance%SMINN_pool(:,:)            = 0.5_dp    !! set to 0.5 at begining 
                                                          
       IF (PRESENT(cbalance)) THEN
          !! Initialize nitrogen pools consistently with carbon pools          ! brp

          nbalance%Npool_green(:,:)           = cbalance%Cpool_green(:,:)/cn_green
          nbalance%Npool_woods(:,:)           = cbalance%Cpool_woods(:,:)/cn_woods
          nbalance%Npool_litter_green_ag(:,:) = cbalance%Cpool_litter_green_ag (:,:)/cn_litter_green
          nbalance%Npool_litter_green_bg(:,:) = cbalance%Cpool_litter_green_bg (:,:)/cn_litter_green
          nbalance%Npool_litter_wood_ag(:,:)  = cbalance%Cpool_litter_wood_ag(:,:)/cn_litter_wood
          nbalance%Npool_litter_wood_bg(:,:)  = cbalance%Cpool_litter_wood_bg(:,:)/cn_litter_wood
          nbalance%Npool_slow(:,:)            = cbalance%Cpool_slow (:,:)/cn_slow

          nbalance%Npool_mobile(:,:)          = nbalance%Npool_green(:,:)        !! use N-pool green for initialization
          nbalance%SMINN_pool(:,:)            = nbalance%Npool_litter_green_bg(:,:) +0.5_dp !! use N-pool fast for initialization

          WRITE(*,*) "initializeNPools(): Nitrogen pools initialized consistently via C/N-ratios from C-pools"
      END IF
    else                                         
       !! open the file
       write (*,*) 'initializeNPools(): read Npools from: ',filename
       iret = nf_open(filename,nf_nowrite,ncid)
       call check_err(iret)
       ! check what is in the input-file
       iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
       call check_err(iret)
       if (debug) write (*,*) 'initializeNPools(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim
       ! get the dimension name and length
       allocate (ndimlen(ndim))
       allocate (fdimnam(ndim))
       do i=1,ndim
          iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
          call check_err(iret)
          if (debug) write (*,*) 'initializeNPools(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
       end do
       ! get variable names, types and shapes
       allocate (fvarnam(nvar))
       allocate (nvartyp(nvar))
       allocate (nvardim(nvar))
       allocate (nvardimid(nvar,100))
       allocate (nvaratt(nvar))
       do i=1,nvar
          iret = nf_inq_var(ncid,i,fvarnam(i),nvartyp(i),nvardim(i),nin,nvaratt(i))
          call check_err(iret)
          do ii=1,nvardim(i)
             nvardimid(i,ii)=nin(ii)
          end do
          if (debug) write (*,*) 'initializeNPools(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
               i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
       end do
       ! get data
       noTimeSteps = 0
       do i=1,nvar
          if (fvarnam(i) == 'Npool_green') then
             allocate (read_Npools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  CALL finish('initializeNPools','dimensions of Npool_green wrong')
             iret = nf_get_var_real(ncid,i,read_Npools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   call finish('initializeNPools','inconsistent length of time series for nitrogen pools in '//trim(filename))
                end if
             end if
             do itile = 1,ntiles
                nbalance%Npool_green(:,itile) = PACK(read_Npools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Npools)
!! check if Npool_mobile intilization file needed : Switch case On/Off  
          else if(fvarnam(i) == 'Npool_mobile') then
             allocate (read_Npools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  call finish('initializeNPools','dimensions Npool_mobile wrong')
             iret = nf_get_var_real(ncid,i,read_Npools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   call finish('initializeNPools', 'inconsistent length of time series for nitrogen pools in '//trim(filename))
                end if
             end if
             do itile = 1,ntiles
                nbalance%Npool_mobile(:,itile) = PACK(read_Npools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Npools)

!!$             nbalance%Npool_mobile(:,:) = 0.0_dp      !! set to Npool_mobile initialization at 0.5 in transient offline(1992)
                                                      !! or read from prev file ?? 

          else if(fvarnam(i) == 'Npool_woods') then
             allocate (read_Npools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  call finish('initializeNPools', 'dimensions of Npool_woods wrong')
             iret = nf_get_var_real(ncid,i,read_Npools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   call finish('initializeNPools', 'inconsistent length of time series for nitrogen pools in '//trim(filename))
                end if
             end if
             do itile = 1,ntiles
                nbalance%Npool_woods(:,itile) = PACK(read_Npools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Npools)
          else if(fvarnam(i) == 'Npool_litter_green_ag') then
             allocate (read_Npools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  call finish('initializeNPools', 'dimensions Npool_litter_green_ag wrong')
             iret = nf_get_var_real(ncid,i,read_Npools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   call finish('initializeNPools', 'inconsistent length of time series for nitrogen pools in '//trim(filename))
                end if
             end if
             do itile = 1,ntiles
                nbalance%Npool_litter_green_ag(:,itile) = PACK(read_Npools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Npools)

          else if(fvarnam(i) == 'Npool_litter_green_bg') then
             allocate (read_Npools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  call finish('initializeNPools', 'dimensions of Npool_litter_green_bg wrong')
             iret = nf_get_var_real(ncid,i,read_Npools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   call finish('initializeNPools', 'inconsistent length of time series for nitrogen pools in '//trim(filename))
                end if
             end if
             do itile = 1,ntiles
                nbalance%Npool_litter_green_bg(:,itile) = PACK(read_Npools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Npools)

          else if(fvarnam(i) == 'Npool_litter_wood_ag') then
             allocate (read_Npools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  call finish('initializeNPools', 'dimensions Npool_litter_wood_ag wrong')
             iret = nf_get_var_real(ncid,i,read_Npools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   call finish('initializeNPools', 'inconsistent length of time series for nitrogen pools in '//trim(filename))
                end if
             end if
             do itile = 1,ntiles
                nbalance%Npool_litter_wood_ag(:,itile) = PACK(read_Npools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Npools)

          else if(fvarnam(i) == 'Npool_litter_wood_bg') then
             allocate (read_Npools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  call finish('initializeNPools', 'dimensions Npool_litter_wood_bg wrong')
             iret = nf_get_var_real(ncid,i,read_Npools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   call finish('initializeNPools', 'inconsistent length of time series for nitrogen pools in '//trim(filename))
                end if
             end if
             do itile = 1,ntiles
                nbalance%Npool_litter_wood_bg(:,itile) = PACK(read_Npools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Npools)
          else if(fvarnam(i) == 'Npool_slow') then
             allocate (read_Npools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  call finish('initializeNPools', 'dimensions Npool_slow wrong')
             iret = nf_get_var_real(ncid,i,read_Npools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   call finish('initializeNPools', 'inconsistent length of time series for nitrogen pools in '//trim(filename))
                end if
             end if
             do itile = 1,ntiles
                nbalance%Npool_slow(:,itile) = PACK(read_Npools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Npools)

!! check if SMINN_pool intilization file needed : Switch case On/Off  
          else if(fvarnam(i) == 'SMINN_pool') then
             allocate (read_Npools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),ndimlen(nvardimid(i,4))))
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) &
                  call finish('initializeNPools', 'dimensions SMINN pool wrong')
             iret = nf_get_var_real(ncid,i,read_Npools)
             call check_err(iret)
             if(noTimeSteps ==0) then
                noTimeSteps = ndimlen(nvardimid(i,4))
             else
                if(noTimeSteps /= ndimlen(nvardimid(i,4))) then
                   call finish('initializeNPools', 'inconsistent length of time series for nitrogen pools in '//trim(filename))
                end if
             end if
             do itile = 1,ntiles
                nbalance%SMINN_pool(:,itile) = PACK(read_Npools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             deallocate(read_Npools)

        end if
       end do
       write(*,*) "initializeNPools(): Nitrogen pools read in. Took values from last time step ",noTimeSteps," found in file."

       ! close the file
       iret = nf_close(NCID)
       call check_err(iret)
       ! deallocate
       deallocate(fdimnam,ndimlen)
       deallocate(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
    end if
    
  end subroutine initializeNPools


  ! --- initSurface() -----------------------------------------------------------------------------------

  subroutine initSurface(filename, grid, lctlib, surface)

    USE mo_jsbach_lctlib,         ONLY: lctlib_type
    USE mo_land_surface,          ONLY: land_surface_type, fract_small

    include 'netcdf.inc'

    character(len=*),             intent(in)    :: filename
    TYPE(grid_type),              INTENT(in)    :: grid
    TYPE(lctlib_type),            INTENT(in)    :: lctlib
    TYPE(land_surface_type),      INTENT(inout) :: surface

    !! locals

    integer :: iret,ncid,ndim,nvar,nattr,nunlimdim
    integer,allocatable            :: ndimlen(:)
    character(len=128),allocatable :: fdimnam(:),fvarnam(:)
    integer, allocatable           :: nvartyp(:),nvardim(:),nvaratt(:)
    integer,allocatable            :: nvardimid(:,:)
    integer,allocatable            :: coverType(:,:,:)
    real(dp),allocatable           :: coverFract(:,:,:)
    real(dp),allocatable           :: vegRatioMax(:,:)
    real(dp),allocatable           :: naturalVeg(:,:,:)
    integer :: i, ii, itile, ilct
    integer ::  nin(100)
    
    !! open the file
    write (*,*) 'initSurface(): read surface data from: ',trim(filename)
    iret = nf_open(filename,nf_nowrite,ncid)
    call check_err(iret)
    ! check what is in the input-file
    iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
    call check_err(iret)
    if (debug) write (*,*) 'initSurface(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim
    ! get the dimension name and length
    allocate (ndimlen(ndim))
    allocate (fdimnam(ndim))
    do i=1,ndim
       iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
       call check_err(iret)
       if (debug) write (*,*) 'initSurface(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
    end do
    ! get variable names, types and shapes
    allocate (fvarnam(nvar))
    allocate (nvartyp(nvar))
    allocate (nvardim(nvar))
    allocate (nvardimid(nvar,100))
    allocate (nvaratt(nvar))
    do i=1,nvar
       iret = nf_inq_var(ncid,i,fvarnam(i),nvartyp(i),nvardim(i),nin,nvaratt(i))
       call check_err(iret)
       do ii=1,nvardim(i)
          nvardimid(i,ii)=nin(ii)
       end do
       if (debug) write (*,*) 'initSurface(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
            i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
    end do
    ! get ntiles
    surface%ntiles = 0
    do i=1,ndim
       if (fdimnam(i) == 'ntiles') surface%ntiles = ndimlen(i)
    end do
    if (surface%ntiles == 0) call finish('initSurface', 'ntiles not found') 
    if(debug) write(*,*) "initSurface(): ntiles=",surface%ntiles

    ! get cover type data
    do i=1,nvar
       if (fvarnam(i) == 'cover_type') then
          allocate (coverType(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
               ndimlen(nvardimid(i,3))))
          if (debug) write (*,*) 'initSurface(): cover_type dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)),ndimlen(nvardimid(i,3))
          iret = nf_get_var_int(ncid,i,coverType)
          call check_err(iret)
       elseif (fvarnam(i) == 'cover_fract') then
          allocate (coverFract(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
               ndimlen(nvardimid(i,3))))
          if (debug) write (*,*) 'initSurface(): cover_fract dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)),ndimlen(nvardimid(i,3))
          iret = nf_get_var_double(ncid,i,coverFract)
          call check_err(iret)
       elseif (fvarnam(i) == 'natural_veg') then
          allocate (naturalVeg(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
               ndimlen(nvardimid(i,3))))
          if (debug) write (*,*) 'initSurface(): natural_veg dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)),ndimlen(nvardimid(i,3))
          iret = nf_get_var_double(ncid,i,naturalVeg)
          call check_err(iret)
       elseif (fvarnam(i) == 'veg_ratio_max') then
          allocate (vegRatioMax(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))))
          if (debug) write (*,*) 'initSurface(): veg_ratio_max dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))
          iret = nf_get_var_double(ncid,i,vegRatioMax)
          call check_err(iret)
       end if
    end do

    ! pack cover type

    if(.not. allocated(coverType)) then
       call finish('initSurface', 'Variable cover_type not found in input file')
    else
       allocate (surface%cover_type(grid%nland,surface%ntiles))
       do itile=1,surface%ntiles
          surface%cover_type(:,itile) = PACK(coverType(:,:,itile),mask=grid%mask)
       end do
       deallocate(coverType)
    end if
    
    ! pack cover fractions

    if(.not. allocated(coverFract)) then
       call finish('initSurface', 'Variable cover_fract not found in input file')
    else
       allocate (surface%cover_fract(grid%nland,surface%ntiles))
       do itile=1,surface%ntiles
          surface%cover_fract(:,itile) = PACK(coverFract(:,:,itile),mask=grid%mask)
       end do
       deallocate(coverFract)
    end if

    ! pack natural vegetation

    if(.not. allocated(naturalVeg)) then
       call finish('initSurface', 'Variable natural_veg not found in input file')
    else
       allocate (surface%cover_fract_pot(grid%nland,surface%ntiles))
       do itile=1,surface%ntiles
          surface%cover_fract_pot(:,itile) = PACK(naturalVeg(:,:,itile),mask=grid%mask)
       end do
       deallocate(naturalVeg)
    end if

    ! pack surface fractions

    if(.not. allocated(vegRatioMax)) then
       call finish('initSurface', 'Variable veg_ratio_max not found in input file')
    else
       allocate (surface%veg_ratio_max(grid%nland))
       surface%veg_ratio_max(:) = PACK(vegRatioMax(:,:),mask=grid%mask)
       deallocate(vegRatioMax)
    end if

    ! Define rock fraction (used with dynamic vegetation)
    ALLOCATE (surface%rock_fract(grid%nland))
    surface%rock_fract(:) = 0._dp

    ! Define logical masks for convenience

    ALLOCATE(surface%is_bare_soil(grid%nland, surface%ntiles), &
             surface%is_vegetation(grid%nland, surface%ntiles), &
             surface%is_C4vegetation(grid%nland, surface%ntiles), &
             surface%is_naturalVeg(grid%nland, surface%ntiles), &
             surface%is_forest(grid%nland, surface%ntiles), &
             surface%is_grass(grid%nland, surface%ntiles), &
             surface%is_pasture(grid%nland, surface%ntiles), &
             surface%is_crop(grid%nland, surface%ntiles), &
             surface%is_lake(grid%nland, surface%ntiles), &
             surface%is_glacier(grid%nland, surface%ntiles), &
             surface%is_present(grid%nland, surface%ntiles) )

    surface%is_bare_soil    = .FALSE.
    surface%is_vegetation   = .FALSE.
    surface%is_C4vegetation = .FALSE.
    surface%is_naturalVeg   = .FALSE.
    surface%is_forest       = .FALSE.
    surface%is_grass        = .FALSE.
    surface%is_pasture      = .FALSE.
    surface%is_crop         = .FALSE.
    surface%is_lake         = .FALSE.
    surface%is_glacier      = .FALSE.
    surface%is_present      = .TRUE.

    DO ilct=1,lctlib%nlct
       WHERE (surface%cover_type == ilct .AND. surface%cover_fract > 0._dp)
          surface%is_vegetation   = lctlib%NaturalVegFlag(ilct) .OR. lctlib%CropFlag(ilct) .OR. lctlib%PastureFlag(ilct)
          surface%is_C4vegetation = lctlib%C4flag(ilct) 
          surface%is_naturalVeg   = lctlib%NaturalVegFlag(ilct)
          surface%is_forest       = lctlib%ForestFlag    (ilct)
          surface%is_grass        = lctlib%GrassFlag     (ilct)   
          surface%is_pasture      = lctlib%PastureFlag   (ilct)
          surface%is_crop         = lctlib%CropFlag      (ilct)
          surface%is_lake         = lctlib%LakeFlag      (ilct)
          surface%is_glacier      = lctlib%GlacierFlag   (ilct)
          surface%is_bare_soil    = lctlib%BareSoilFlag  (ilct)
       END WHERE
    END DO

    ! In the current model version no fractional glacier cells are allowed.
    IF (ANY(surface%is_glacier(:,:) .AND. &
         surface%cover_fract(:,:) > fract_small + EPSILON(1._dp) .AND. &
         surface%cover_fract(:,:) < 1._dp - EPSILON(1._dp))) THEN
       CALL finish('initSurface','Fractional glacier cells are not allowed in the current model version')
    END IF

    ! final deallocations

    deallocate (ndimlen)
    deallocate (fdimnam)
    deallocate (fvarnam)
    deallocate (nvartyp)
    deallocate (nvardim)
    deallocate (nvardimid)
    deallocate (nvaratt)

  end subroutine initSurface

  ! --- initSoil() -----------------------------------------------------------------------------------

  subroutine initSoil(filename, grid, soil_param)

    USE mo_soil,       ONLY: soil_param_type

    include 'netcdf.inc'

    character(len=*),             intent(in)    :: filename
    TYPE(grid_type),              intent(in)    :: grid
    TYPE(soil_param_type),        intent(inout) :: soil_param

    !! locals

    integer :: iret,ncid,ndim,nvar,nattr,nunlimdim
    integer,allocatable            :: ndimlen(:)
    character(len=128),allocatable :: fdimnam(:),fvarnam(:)
    integer, allocatable           :: nvartyp(:),nvardim(:),nvaratt(:)
    integer,allocatable            :: nvardimid(:,:)
    real(dp),allocatable           :: max_moisture(:,:)
    integer :: i, ii
    integer ::  nin(100)
    
    !! open the file
    write (*,*) 'initSoil(): read soil data from: ',trim(filename)
    iret = nf_open(filename,nf_nowrite,ncid)
    call check_err(iret)
    ! check what is in the input-file
    iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
    call check_err(iret)
    if (debug) write (*,*) 'initSoil(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim
    ! get the dimension name and length
    allocate (ndimlen(ndim))
    allocate (fdimnam(ndim))
    do i=1,ndim
       iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
       call check_err(iret)
       if (debug) write (*,*) 'initSoil(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
    end do
    ! get variable names, types and shapes
    allocate (fvarnam(nvar))
    allocate (nvartyp(nvar))
    allocate (nvardim(nvar))
    allocate (nvardimid(nvar,100))
    allocate (nvaratt(nvar))
    do i=1,nvar
       iret = nf_inq_var(ncid,i,fvarnam(i),nvartyp(i),nvardim(i),nin,nvaratt(i))
       call check_err(iret)
       do ii=1,nvardim(i)
          nvardimid(i,ii)=nin(ii)
       end do
       if (debug) write (*,*) 'initSoil(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
            i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
    end do

    ! get cover type data
    do i=1,nvar
       if (fvarnam(i) == 'maxmoist') then      
          allocate (max_moisture(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))))
          if (debug) write (*,*) 'initSoil(): maxmoist dimensions: ', &
               ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2))
          iret = nf_get_var_double(ncid,i,max_moisture)
          call check_err(iret)
       end if
    end do

    ! pack max_moisture

    if(.not. allocated(max_moisture)) then
       call finish('initSoil', 'Variable maxmoist not found in input file')
    else
       allocate (soil_param%maxMoisture(grid%nland))
       soil_param%maxMoisture(:) = PACK(max_moisture(:,:),mask=grid%mask)
       deallocate(max_moisture)
    end if

    ! final deallocations

    deallocate (ndimlen)
    deallocate (fdimnam)
    deallocate (fvarnam)
    deallocate (nvartyp)
    deallocate (nvardim)
    deallocate (nvardimid)
    deallocate (nvaratt)

  end subroutine initSoil

  ! --- initCbalance()  -----------------------------------------------------------------------------------

  SUBROUTINE initCbalance(grid, ntiles, cbalance, cbalance_diag, UseLandcoverTransitions, lcc_scheme, with_yasso)

    TYPE(grid_type),         INTENT(in)    :: grid
    INTEGER,                 INTENT(in)    :: ntiles
    TYPE(cbalance_type),     INTENT(inout) :: cbalance
    TYPE(cbal_offline_type), INTENT(inout) :: cbalance_diag
    LOGICAL,                 INTENT(in)    :: UseLandcoverTransitions
    INTEGER,                 INTENT(in)    :: lcc_scheme
    LOGICAL,                 INTENT(in)    :: with_yasso
    
    !! --- Allocate cpools

    allocate (cbalance%Cpool_green(grid%nland,ntiles))
    allocate (cbalance%Cpool_reserve(grid%nland,ntiles))
    allocate (cbalance%Cpool_woods(grid%nland,ntiles))
    allocate (cbalance%Cpool_litter_green_ag(grid%nland,ntiles))
    allocate (cbalance%Cpool_litter_green_bg(grid%nland,ntiles))
    allocate (cbalance%Cpool_litter_wood_ag(grid%nland,ntiles))
    allocate (cbalance%Cpool_litter_wood_bg(grid%nland,ntiles))
    allocate (cbalance%Cpool_slow(grid%nland,ntiles))
    IF (with_yasso) THEN
      allocate (cbalance%YCpool_acid_ag1(grid%nland,ntiles))
      allocate (cbalance%YCpool_water_ag1(grid%nland,ntiles))
      allocate (cbalance%YCpool_ethanol_ag1(grid%nland,ntiles))
      allocate (cbalance%YCpool_nonsoluble_ag1(grid%nland,ntiles))
      allocate (cbalance%YCpool_acid_bg1(grid%nland,ntiles))
      allocate (cbalance%YCpool_water_bg1(grid%nland,ntiles))
      allocate (cbalance%YCpool_ethanol_bg1(grid%nland,ntiles))
      allocate (cbalance%YCpool_nonsoluble_bg1(grid%nland,ntiles))
      allocate (cbalance%YCpool_humus_1(grid%nland,ntiles))
      allocate (cbalance%YCpool_acid_ag2(grid%nland,ntiles))
      allocate (cbalance%YCpool_water_ag2(grid%nland,ntiles))
      allocate (cbalance%YCpool_ethanol_ag2(grid%nland,ntiles))
      allocate (cbalance%YCpool_nonsoluble_ag2(grid%nland,ntiles))
      allocate (cbalance%YCpool_acid_bg2(grid%nland,ntiles))
      allocate (cbalance%YCpool_water_bg2(grid%nland,ntiles))
      allocate (cbalance%YCpool_ethanol_bg2(grid%nland,ntiles))
      allocate (cbalance%YCpool_nonsoluble_bg2(grid%nland,ntiles))
      allocate (cbalance%YCpool_humus_2(grid%nland,ntiles))
    else ! In this case, the variables below are not use, but we need to be able to pass them through interfaces using indices
      cbalance%YCpool_acid_ag1       => cbalance%Cpool_green(:,:)
      cbalance%YCpool_water_ag1      => cbalance%Cpool_green(:,:)
      cbalance%YCpool_ethanol_ag1    => cbalance%Cpool_green(:,:)
      cbalance%YCpool_nonsoluble_ag1 => cbalance%Cpool_green(:,:)
      cbalance%YCpool_acid_bg1       => cbalance%Cpool_green(:,:)
      cbalance%YCpool_water_bg1      => cbalance%Cpool_green(:,:)
      cbalance%YCpool_ethanol_bg1    => cbalance%Cpool_green(:,:)
      cbalance%YCpool_nonsoluble_bg1 => cbalance%Cpool_green(:,:)
      cbalance%YCpool_humus_1        => cbalance%Cpool_green(:,:)
      cbalance%YCpool_acid_ag2       => cbalance%Cpool_green(:,:)
      cbalance%YCpool_water_ag2      => cbalance%Cpool_green(:,:)
      cbalance%YCpool_ethanol_ag2    => cbalance%Cpool_green(:,:)
      cbalance%YCpool_nonsoluble_ag2 => cbalance%Cpool_green(:,:)
      cbalance%YCpool_acid_bg2       => cbalance%Cpool_green(:,:)
      cbalance%YCpool_water_bg2      => cbalance%Cpool_green(:,:)
      cbalance%YCpool_ethanol_bg2    => cbalance%Cpool_green(:,:)
      cbalance%YCpool_nonsoluble_bg2 => cbalance%Cpool_green(:,:)
      cbalance%YCpool_humus_2        => cbalance%Cpool_green(:,:)
    endif
    nullify(cbalance%Cpool_onSite,        &
            cbalance%Cpool_paper,         &
            cbalance%Cpool_construction,  &
            cbalance%Cpool_paper_harvest, &
            cbalance%Cpool_construction_harvest)
    if (lcc_scheme==2) then
       allocate(cbalance%Cpool_onSite      (grid%nland), &
                cbalance%Cpool_paper       (grid%nland), &
                cbalance%Cpool_construction(grid%nland))
       cbalance%Cpool_onSite      (:) = 0._dp
       cbalance%Cpool_paper       (:) = 0._dp
       cbalance%Cpool_construction(:) = 0._dp
       if (UseLandcoverTransitions) then
          allocate(cbalance%Cpool_paper_harvest       (grid%nland), &
                   cbalance%Cpool_construction_harvest(grid%nland))
          cbalance%Cpool_paper_harvest       (:) = 0._dp
          cbalance%Cpool_construction_harvest(:) = 0._dp
       endif
    endif

    allocate (cbalance%soil_respiration(grid%nland,ntiles))
    allocate (cbalance%soil_respiration_pot(grid%nland,ntiles))
    allocate (cbalance%NPP_flux_correction(grid%nland,ntiles))
    allocate (cbalance%excess_NPP(grid%nland,ntiles))
    allocate (cbalance%root_exudates(grid%nland,ntiles))
    allocate (cbalance%Cflux_herbivory(grid%nland,ntiles))
    allocate (cbalance%Cflux_herbivory_LG(grid%nland,ntiles))
    allocate (cbalance%Cflux_herbivory_2_atm(grid%nland,ntiles))

    allocate (cbalance_diag%LAI_yDayMean(grid%nland,ntiles,31))
    allocate (cbalance_diag%LAI_previousDayMean(grid%nland,ntiles))
    allocate (cbalance_diag%NPP_yDayMean(grid%nland,ntiles,31))
    allocate (cbalance_diag%topSoilTemp_yDayMean(grid%nland,ntiles,31))
    allocate (cbalance_diag%alpha_yDayMean(grid%nland,ntiles,31))

    allocate (cbalance_diag%pseudo_temp_yDay(grid%nland,ntiles,31))
    allocate (cbalance_diag%pseudo_precip_yDay(grid%nland,ntiles,31))

    allocate (cbalance_diag%box_Cpools_total(grid%nland))
    cbalance_diag%box_Cpools_total(:)=0.0_dp
    allocate (cbalance_diag%box_test_Ccons(grid%nland))
    cbalance_diag%box_test_Ccons(:)=0.0_dp
    allocate (cbalance_diag%box_CO2_flux_2_atm(grid%nLand))
    cbalance_diag%box_CO2_flux_2_atm(:) = 0.0_dp
    allocate (cbalance_diag%avg_Cpool_green(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cpool_reserve(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cpool_woods(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cpool_litter_green_ag(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cpool_litter_green_bg(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cpool_litter_wood_ag(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cpool_litter_wood_bg(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cpool_slow(grid%nland,ntiles))

! Average pools of Yasso, THT 21.12.2007, DSG
    if (with_yasso) then
      allocate (cbalance_diag%avg_YCpool_acid_ag1(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_water_ag1(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_ethanol_ag1(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_nonsoluble_ag1(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_acid_bg1(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_water_bg1(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_ethanol_bg1(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_nonsoluble_bg1(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_humus_1(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_acid_ag2(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_water_ag2(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_ethanol_ag2(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_nonsoluble_ag2(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_acid_bg2(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_water_bg2(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_ethanol_bg2(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_nonsoluble_bg2(grid%nland,ntiles))
      allocate (cbalance_diag%avg_YCpool_humus_2(grid%nland,ntiles))
    endif
! end THT

    allocate (cbalance_diag%avg_soil_respiration(grid%nland,ntiles))
    allocate (cbalance_diag%avg_soil_respiration_pot(grid%nland,ntiles))
    allocate (cbalance_diag%avg_NPP_yDayMean(grid%nland,ntiles))
    allocate (cbalance_diag%avg_NPP_flux_correction(grid%nland,ntiles))
    allocate (cbalance_diag%avg_excess_NPP(grid%nland,ntiles))
    allocate (cbalance_diag%avg_root_exudates(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cflux_herbivory(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cflux_herbivory_LG(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cflux_herbivory_2_atm(grid%nland,ntiles))
    allocate (cbalance_diag%avg_box_NEP(grid%nland))
    allocate (cbalance%NPP_act_yDayMean(grid%nland,ntiles))
    cbalance%NPP_act_yDayMean(:,:)=0.0_dp
    allocate (cbalance%frac_litter_wood_new(grid%nland,ntiles))
    cbalance%frac_litter_wood_new(:,:)=0.0
    allocate (cbalance_diag%avg_NPP_act_yDayMean(grid%nland,ntiles))
    cbalance_diag%avg_NPP_act_yDayMean(:,:)=0.0_dp
    allocate (cbalance_diag%testCconserv(grid%nland,ntiles))
    allocate (cbalance_diag%avg_Cpool_onSite              (grid%nland), &
              cbalance_diag%avg_Cpool_paper               (grid%nland), &
              cbalance_diag%avg_Cpool_construction        (grid%nland), &
              cbalance_diag%avg_Cpool_paper_harvest       (grid%nland), &
              cbalance_diag%avg_Cpool_construction_harvest(grid%nland))

  end subroutine initCbalance


! --- initNbalance()  -----------------------------------------------------------------------------------

  subroutine initNbalance(grid, ntiles, nbalance, nbalance_offline)

    USE mo_cbal_bethy, ONLY: nbalance_type

    TYPE(grid_type),          INTENT(in) :: grid
    INTEGER,                  INTENT(in) :: ntiles
    TYPE(nbalance_type),      INTENT(out):: nbalance
    TYPE(nbal_offline_type),  INTENT(out):: nbalance_offline

    !! --- Allocate npools

    ALLOCATE (nbalance%Npool_green(grid%nland,ntiles))
    ALLOCATE (nbalance%Npool_mobile(grid%nland,ntiles))
    ALLOCATE (nbalance%Npool_woods(grid%nland,ntiles))
    ALLOCATE (nbalance%Npool_litter_green_ag(grid%nland,ntiles))
    ALLOCATE (nbalance%Npool_litter_green_bg(grid%nland,ntiles))
    ALLOCATE (nbalance%Npool_litter_wood_ag(grid%nland,ntiles))
    ALLOCATE (nbalance%Npool_litter_wood_bg(grid%nland,ntiles))
    ALLOCATE (nbalance%Npool_slow(grid%nland,ntiles))
    ALLOCATE (nbalance%SMINN_pool(grid%nland,ntiles))
    ALLOCATE (nbalance%redFact_Nlimitation(grid%nland,ntiles))
    nbalance%redFact_Nlimitation(:,:)=1.0_dp

    ALLOCATE (nbalance%minNflux_litter_green_ag(grid%nland,ntiles))
    nbalance%minNflux_litter_green_ag(:,:)=0.0_dp
    ALLOCATE (nbalance%minNflux_litter_green_bg(grid%nland,ntiles)) 
    nbalance%minNflux_litter_green_bg(:,:)=0.0_dp
    ALLOCATE (nbalance%minNflux_litter_wood_ag(grid%nland,ntiles))
    nbalance%minNflux_litter_wood_ag=0.0_dp
    ALLOCATE (nbalance%minNflux_litter_wood_bg(grid%nland,ntiles))
    nbalance%minNflux_litter_wood_bg=0.0_dp
    ALLOCATE (nbalance%minNflux_slow(grid%nland,ntiles))
    nbalance%minNflux_slow=0.0_dp
    ALLOCATE (nbalance%sum_N_pools(grid%nland,ntiles)) 
    nbalance%sum_N_pools=0.0_dp

    ALLOCATE (nbalance%nfix_to_sminn(grid%nland,ntiles))
    nbalance%nfix_to_sminn(:,:)=0.0_dp
    ALLOCATE (nbalance%ndep_to_sminn(grid%nland,ntiles))
    nbalance%ndep_to_sminn(:,:)=0.0_dp
    ALLOCATE (nbalance%nfert_to_sminn(grid%nland,ntiles))
    nbalance%nfert_to_sminn(:,:)=0.0_dp
    ALLOCATE (nbalance%nfert_forc_2d(grid%nland,ntiles))          ! not yet tested for cbalone
    nbalance%nfert_forc_2d(:,:)=0.0_dp
    ALLOCATE (nbalance%ndep_forc(grid%nland))
    nbalance%ndep_forc(:)=0.0_dp
    ALLOCATE (nbalance%nfert_forc(grid%nland))
    nbalance%nfert_forc(:)=0.0_dp                                  ! for usage of read_nitrogen_fert/dep by cbalone
    ALLOCATE (nbalance%sminn_to_denit(grid%nland,ntiles))
    nbalance%sminn_to_denit(:,:)=0.0_dp
    ALLOCATE (nbalance%sminn_leach(grid%nland,ntiles))
    nbalance%sminn_leach(:,:)=0.0_dp
    ALLOCATE (nbalance%N2O_emissions_depfix(grid%nland,ntiles))
    nbalance%N2O_emissions_depfix(:,:)=0.0_dp
    ALLOCATE (nbalance%N2O_emissions_mineraliz(grid%nland,ntiles))
    nbalance%N2O_emissions_mineraliz(:,:)=0.0_dp
    ALLOCATE (nbalance%N2O_emissions_slow(grid%nland,ntiles))
    nbalance%N2O_emissions_slow(:,:)=0.0_dp
    ALLOCATE (nbalance%N2O_emissions_nfert(grid%nland,ntiles))
    nbalance%N2O_emissions_nfert(:,:)=0.0_dp
    ALLOCATE (nbalance%N2O_emissions_grazing(grid%nland,ntiles))
    nbalance%N2O_emissions_grazing(:,:)=0.0_dp

    ALLOCATE (nbalance%Nplant_demand(grid%nland,ntiles))
    nbalance%Nplant_demand(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_Nplant_demand(grid%nland,ntiles))
    nbalance_offline%avg_Nplant_demand(:,:)=0.0_dp
    ALLOCATE (nbalance%Nsoil_demand(grid%nland,ntiles))
    nbalance%Nsoil_demand(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_Nsoil_demand(grid%nland,ntiles))
    nbalance_offline%avg_Nsoil_demand(:,:)=0.0_dp
    ALLOCATE (nbalance%Ntotal_demand(grid%nland,ntiles))
    nbalance%Ntotal_demand(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_Ntotal_demand(grid%nland,ntiles))
    nbalance_offline%avg_Ntotal_demand(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_minNflux_litter_green_ag(grid%nland,ntiles))
    nbalance_offline%avg_minNflux_litter_green_ag(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_minNflux_litter_green_bg(grid%nland,ntiles))
    nbalance_offline%avg_minNflux_litter_green_bg(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_minNflux_litter_wood_ag(grid%nland,ntiles))
    nbalance_offline%avg_minNflux_litter_wood_ag(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_minNflux_litter_wood_bg(grid%nland,ntiles))
    nbalance_offline%avg_minNflux_litter_wood_bg(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_minNflux_slow(grid%nland,ntiles))
    nbalance_offline%avg_minNflux_slow(:,:)=0.0_dp

    ALLOCATE (nbalance%SMINN_herbivory(grid%nland,ntiles))
    nbalance%SMINN_herbivory(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_SMINN_herbivory(grid%nland,ntiles))
    nbalance_offline%avg_SMINN_herbivory(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_nfix_to_sminn(grid%nland,ntiles))
    nbalance_offline%avg_nfix_to_sminn(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_sminn_to_denit(grid%nland,ntiles))
    nbalance_offline%avg_sminn_to_denit(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_sminn_leach(grid%nland,ntiles))
    nbalance_offline%avg_sminn_leach(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_nfert_to_sminn(grid%nland,ntiles))
    nbalance_offline%avg_nfert_to_sminn(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_N2O_emissions_ExternalN(grid%nland,ntiles))
    nbalance_offline%avg_N2O_emissions_ExternalN(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_N2O_emissions_mineraliz(grid%nland,ntiles))
    nbalance_offline%avg_N2O_emissions_mineraliz(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_N2O_emissions_nfert(grid%nland,ntiles))
    nbalance_offline%avg_N2O_emissions_nfert(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_N2O_emissions_grazing(grid%nland,ntiles))
    nbalance_offline%avg_N2O_emissions_grazing(:,:)=0.0_dp
    ALLOCATE (nbalance%NetEcosyst_N_flux(grid%nland,ntiles))
    nbalance%NetEcosyst_N_flux(:,:)=0.0_dp
    ALLOCATE (nbalance_offline%avg_NetEcosyst_N_flux(grid%nland,ntiles))
    nbalance_offline%avg_NetEcosyst_N_flux(:,:)=0.0_dp
    ALLOCATE (nbalance%NPP_run_mean(grid%nland,ntiles))
    nbalance%NPP_run_mean(:,:)=0.0_dp

    ALLOCATE (nbalance_offline%testNconserv(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%test_Ngreen(grid%nland,ntiles)) 
    ALLOCATE (nbalance_offline%test_Nwoods(grid%nland,ntiles)) 
    ALLOCATE (nbalance_offline%test_Nlitter_green_ag(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%test_Nlitter_green_bg(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%test_Nlitter_wood_ag(grid%nland,ntiles)) 
    ALLOCATE (nbalance_offline%test_Nlitter_wood_bg(grid%nland,ntiles)) 
    ALLOCATE (nbalance_offline%test_Nslow(grid%nland,ntiles)) 

    ALLOCATE (nbalance_offline%box_Npools_total(grid%nland))
    ALLOCATE (nbalance_offline%box_test_Ncons(grid%nland))

    ALLOCATE (nbalance_offline%avg_Npool_green(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%avg_Npool_mobile(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%avg_Npool_woods(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%avg_Npool_litter_green_ag(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%avg_Npool_litter_green_bg(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%avg_Npool_litter_wood_ag(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%avg_Npool_litter_wood_bg(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%avg_Npool_slow(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%avg_SMINN_pool(grid%nland,ntiles))
    ALLOCATE (nbalance_offline%Runoff_yDayMean(grid%nland,ntiles,31))

  end subroutine initNbalance

  ! --- check_err --------------------------------------------------------------------------------------------

  subroutine check_err(IRET,text)
    include 'netcdf.inc'
    integer,intent(in)          :: iret
    character(len=*),optional,intent(in) :: text
    if (iret /= NF_NOERR) then
       if(present(text)) write(*,*) 'check_err(): '//trim(text)
       call finish('check_err', 'netcdf says: '//TRIM(nf_strerror(iret)))
    end if
  end subroutine check_err

  ! --- initializeYPools() --------------------------------------------------------------------------------
  !
  ! Initializes yasso carbon pools from a file (any output file from simulation using yasso can be taken as initialization file)
  ! If carbon pools are read in from a file, they must be in LatLon format, packed format is currently not
  ! supported.
  !
  SUBROUTINE initializeYassoPools(grid, cbalance, surface, filename)
    USE mo_cbal_cpools,  ONLY: printCpoolParameters
    USE mo_land_surface, ONLY: land_surface_type

    include 'netcdf.inc'
    type(grid_type),        intent(in)   :: grid
    TYPE(cbalance_type),      INTENT(inout) :: cbalance
    type(land_surface_type),intent(inout):: surface
    character(len=*),OPTIONAL,intent(in) :: filename

    !! --- internal variables

    real(dp),allocatable,dimension(:,:,:,:)  :: read_Cpools
    ! error status return
    integer  ::  iret
    ! netCDF id
    integer  ::  ncid
    ! dimension
    integer  ::  ndim,nvar,nattr,nunlimdim
    integer,allocatable,dimension(:)  ::  ndimlen
    character(len=50),allocatable,dimension(:)  ::  fdimnam
    ! variables
    integer ::  nin(100)
    integer,allocatable ,dimension(:)  ::  nvartyp
    integer,allocatable,dimension(:)   ::  nvardim
    integer,allocatable,dimension(:,:) ::  nvardimid
    integer,allocatable,dimension(:)   ::  nvaratt
    character(len=128),allocatable,dimension(:)  ::  fvarnam
    integer :: noTimeSteps
    logical :: time_dimension
    ! indices
    integer  ::  i,ii,itile

    !! Printout parameters 

    call printCpoolParameters

    !! --- initialize the pools
    
    !! Currently cbalance litter and slow are calculated when yasso is used
    !! but given out.
    !! therefore they must set to zero as they may not be in the cpool file
      cbalance%Cpool_litter_green_ag(:,:)     = 0.0_dp 
      cbalance%Cpool_litter_wood_ag(:,:)      = 0.0_dp 
      cbalance%Cpool_litter_green_bg(:,:)     = 0.0_dp 
      cbalance%Cpool_litter_wood_bg(:,:)      = 0.0_dp 
      cbalance%Cpool_slow(:,:)                = 0.0_dp 
     
    if( .not. present(filename) ) then            !! Initialize from scratch
       cbalance%YCpool_acid_ag1(:,:)           = 0.0_dp
       cbalance%YCpool_water_ag1(:,:)          = 0.0_dp
       cbalance%YCpool_ethanol_ag1(:,:)        = 0.0_dp
       cbalance%YCpool_nonsoluble_ag1(:,:)     = 0.0_dp
       cbalance%YCpool_acid_bg1(:,:)           = 0.0_dp
       cbalance%YCpool_water_bg1(:,:)          = 0.0_dp
       cbalance%YCpool_ethanol_bg1(:,:)        = 0.0_dp
       cbalance%YCpool_nonsoluble_bg1(:,:)     = 0.0_dp
       cbalance%YCpool_humus_1(:,:)             = 0.0_dp
       cbalance%YCpool_acid_ag2(:,:)           = 0.0_dp
       cbalance%YCpool_water_ag2(:,:)          = 0.0_dp
       cbalance%YCpool_ethanol_ag2(:,:)        = 0.0_dp
       cbalance%YCpool_nonsoluble_ag2(:,:)     = 0.0_dp
       cbalance%YCpool_acid_bg2(:,:)           = 0.0_dp
       cbalance%YCpool_water_bg2(:,:)          = 0.0_dp
       cbalance%YCpool_ethanol_bg2(:,:)        = 0.0_dp
       cbalance%YCpool_nonsoluble_bg2(:,:)     = 0.0_dp
       cbalance%YCpool_humus_2(:,:)             = 0.0_dp
    else                                          !! Initialize from file
       !! open the file
       write (*,*) 'initializeYPools(): read YCpools from: ',filename
       iret = nf_open(filename,nf_nowrite,ncid)
       call check_err(iret)
       ! check what is in the input-file
       iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
       call check_err(iret)
       if (debug) write (*,*) 'initializeYPools(): ndim,nvar,nattr,nunlimdim',ndim,nvar,nattr,nunlimdim
       ! get the dimension name and length and find out whether there is a time dimension
       allocate (ndimlen(ndim))
       allocate (fdimnam(ndim))
       time_dimension=.false.
       do i=1,ndim
          iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
          call check_err(iret)
          if (debug) write (*,*) 'initializeYPools(): i,fdimnam,ndimlen',i,trim(fdimnam(i)),ndimlen(i)
          if (TRIM(fdimnam(i)) == 'time') time_dimension=.true.
       end do
       ! get variable names, types and shapes
       allocate (fvarnam(nvar))
       allocate (nvartyp(nvar))
       allocate (nvardim(nvar))
       allocate (nvardimid(nvar,100))
       allocate (nvaratt(nvar))
       do i=1,nvar
          iret = nf_inq_var(ncid,i,fvarnam(i),nvartyp(i),nvardim(i),nin,nvaratt(i))
          call check_err(iret)
          do ii=1,nvardim(i)
             nvardimid(i,ii)=nin(ii)
          end do
          if (debug) write (*,*) 'initializeYPools(): i,fvarnam,nvartyp,nvardim,nvardimid,nvaratt',&
               i,trim(fvarnam(i)),nvartyp(i),nvardim(i),nvardimid(i,1),nvaratt(i)
       end do
       ! get data
       do i=1,nvar
          if (nvardim(i) >= 3) then  ! carbon pools have at least 3 dimensions (lon, lat, tiles)
             if (ndimlen(nvardimid(i,1)) /= grid%nlon .or. ndimlen(nvardimid(i,2)) /= grid%nlat) then
                call finish ('initializeYPools', 'dimensions of '//TRIM(fvarnam(i))//' wrong')
             endif
             if (time_dimension) then
                noTimeSteps = ndimlen(nvardimid(i,nvardim(i)))
             else
                noTimeSteps = 1
             end if
             allocate (read_Cpools(ndimlen(nvardimid(i,1)),ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)),noTimeSteps))
          endif

          if (fvarnam(i) == 'YCpool_acid_ag1') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_acid_ag1(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_acid_ag1"

          else if(fvarnam(i) == 'YCpool_acid_bg1') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_acid_bg1(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_acid_bg1"

          else if(fvarnam(i) == 'YCpool_water_ag1') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_water_ag1(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_water_ag1"

          else if(fvarnam(i) == 'YCpool_water_bg1') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_water_bg1(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_water_bg1"

          else if(fvarnam(i) == 'YCpool_ethanol_ag1') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_ethanol_ag1(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_ethanol_ag1"

          else if(fvarnam(i) == 'YCpool_ethanol_bg1') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_ethanol_bg1(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_ethanol_bg1"

          else if(fvarnam(i) == 'YCpool_nonsoluble_ag1') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_nonsoluble_ag1(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_nonsoluble_ag1"

          else if(fvarnam(i) == 'YCpool_nonsoluble_bg1') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_nonsoluble_bg1(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_nonsoluble_bg1"

          else if(fvarnam(i) == 'YCpool_humus_1') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_humus_1(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_humus_1"

          else if (fvarnam(i) == 'YCpool_acid_ag2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_acid_ag2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_acid_ag2"

          else if(fvarnam(i) == 'YCpool_acid_bg2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_acid_bg2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_acid_bg2"

          else if(fvarnam(i) == 'YCpool_water_ag2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_water_ag2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_water_ag2"

          else if(fvarnam(i) == 'YCpool_water_ag2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_water_ag2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_water_ag2"

          else if(fvarnam(i) == 'YCpool_water_bg2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_water_bg2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_water_bg2"

          else if(fvarnam(i) == 'YCpool_ethanol_ag2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_ethanol_ag2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_ethanol_ag2"

          else if(fvarnam(i) == 'YCpool_ethanol_bg2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_ethanol_bg2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_ethanol_bg2"

          else if(fvarnam(i) == 'YCpool_nonsoluble_ag2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_nonsoluble_ag2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_nonsoluble_ag2"

          else if(fvarnam(i) == 'YCpool_nonsoluble_bg2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_nonsoluble_bg2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_nonsoluble_bg2"

          else if(fvarnam(i) == 'YCpool_humus_2') then
             iret = nf_get_var_double(ncid,i,read_Cpools)
             call check_err(iret)
             do itile = 1,surface%ntiles
                cbalance%YCpool_humus_2(:,itile) = PACK(read_Cpools(:,:,itile,noTimeSteps),MASK=grid%mask)
             end do
             write(*,*) "initializeYPools():   read YCpool_humus_2"

          end if !DSG
          if (allocated(read_Cpools)) deallocate(read_Cpools)
       end do
       write(*,*) "initializePools(): Carbon pools read in. Took values from last time step ",noTimeSteps," found in file."
       
       ! close the file
       iret = nf_close(NCID)
       call check_err(iret)
       ! deallocate
       deallocate(fdimnam,ndimlen)
       deallocate(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
    end if
    
  END SUBROUTINE initializeYassoPools

END MODULE mo_cbalone_memory
