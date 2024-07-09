!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! This module has been contributed by the Max Planck Institute for Biogeochemistry, Jena.
!!
MODULE mo_phenology_ccdas

! 
! This module contains the phenology routine from W. Knorr as implemented in CCDAS. It has been changed to use in JSBACH.
!
!
! Author: G.J. Schuermann, MPI-BGC, Jena,  2010
!
! Remarks:
!    -- The main part of the module has been published by W.Knorr et al. , 2010, Journal of Geophysical Research
!
! (0 GJS)  I have to think about the calc_transpiration_pot and the problem with
! division by zero. At the moment this is a work around that might lead to
! errenous results.
! (1 GJS) Have to pay attention with initialization when dealing with restart.(But I do not think we have to deal with restarts)... 
! (2 GJS) Try to use bethy%canopy_conductance_unlimited for canopy_conductance

  USE mo_kind,        ONLY: dp
  USE mo_exception,   ONLY: finish
  USE mo_jsbach_grid, ONLY: kstart, kend, nidx
  USE mo_time_control,    ONLY: get_date_components, current_date, previous_date,delta_time
  USE mo_jsbach_constants, ONLY: Tmelt, Emissivity, StefanBoltzmann, molarMassDryAir_kg, &
               UniversalGasConstant, SpecificHeatDryAirConstPressure ,vonKarman, Eps
  USE mo_helper,     ONLY: errf, maxx, mins, minx
  USE mo_jsbach_lctlib, ONLY: lctlib_type


  IMPLICIT NONE  


  ! === BEGIN OF PUBLIC PART =======================================================================================================

  PUBLIC :: init_phenology_ccdas       ! Allocates and initializes memory, sets standard parameters and reads in phenology-type data
  PUBLIC :: update_phenology_ccdas     ! Recomputes the LAIs for all phenology types at all grid points.
  PUBLIC ::  reset_phenology_ccdas

  ! --- phenology type coding --- NEVER NEVER CHANGE THIS !!!! from original
  ! routine
  INTEGER, PARAMETER, PUBLIC :: no_of_phenotypes = 5 !! Number of phenology types (not including "unvegetated")
  INTEGER, PARAMETER, PUBLIC :: unvegetated = 0   !! no vegetation present
  INTEGER, PARAMETER, PUBLIC :: evergreen   = 1
  INTEGER, PARAMETER, PUBLIC :: summergreen = 2
  INTEGER, PARAMETER, PUBLIC :: raingreen   = 3
  INTEGER, PARAMETER, PUBLIC :: grasses     = 4
  INTEGER, PARAMETER, PUBLIC :: crops       = 5


  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE ! Make ALL following objects private

  ! === Private declarations =======================================================================================================
  
  ! Output variables
  REAL(dp), SAVE, POINTER :: laim(:,:)   !! : water_limited lai memory
  REAL(dp), SAVE, POINTER :: tmpm(:,:)   !! : temperature memory
  REAL(dp), SAVE, POINTER :: temp_last(:) !! Previous day temperature
  REAL(dp), SAVE, POINTER :: lai_last(:,:) !! Previous day lai
  REAL(dp), SAVE, POINTER :: root_moisture_mean(:,:)  !! daily averaged root-moisture
  REAL(dp), SAVE, POINTER :: root_moisture_sum(:,:)
  REAL(dp), SAVE, POINTER :: evap_pot_last(:,:) !! Prevous day pot-evap
  REAL(dp), SAVE, POINTER :: temp_sum(:)
  REAL(dp), SAVE, POINTER :: lai_sum(:,:)
  REAL(dp), SAVE, POINTER :: evap_pot_sum(:,:)

  ! Some other declaratoins
  INTEGER, SAVE         :: ntiles  !! number of tiles
  REAL(dp)              :: laimmult !! exponent of water-limited calculation
  REAL(dp)              :: tmpmmult !! exponent of temperature memory-calculatoin
  LOGICAL, SAVE         :: first_ts_of_day = .FALSE. !! true during the first time step of a new day
  INTEGER, SAVE        :: time_steps_per_day
  
  ! Here we have the parameters that will be optimized in CCDAS
  !vg REAL(dp),PUBLIC,ALLOCATABLE              :: ptauw(:) !! Water stress time scale
  REAL(dp),PUBLIC                            :: plaimax !!Maximum LAI 
  !vg REAL(dp),PUBLIC                          :: plgr !! This is xi - linear growth rate of LAI shortly after bud burst.
  !vg REAL(dp),PUBLIC,ALLOCATABLE              :: ptphen(:),ptphenr(:) !! is T_phi and T_r in W.Knorr
  !vg REAL(dp),PUBLIC,ALLOCATABLE              :: pdphen(:),pdphenr(:) !! tc tr of W.Knorr
  !vg REAL(dp),PUBLIC,ALLOCATABLE              :: pkl(:) !=1/tau_l (longelivety) is k_l in W.Knorr
 

CONTAINS

  ! --- init_phenology() -----------------------------------------------------------------------------------------------------------
  !
  ! This routine initializes the phenology module. It has to be called before the first time step.
  !
  SUBROUTINE init_phenology_ccdas(g_nland,l_nland,zntiles,isRestart,fileformat,stream,nlct)
    USE mo_linked_list, ONLY: t_stream, LAND, TILES
    USE mo_memory_base, ONLY: new_stream,default_stream_setting, &
                              add =>add_stream_element
    USE mo_time_event,  ONLY: io_time_event
    USE mo_output,        ONLY: land_table
    USE mo_netCDF,      ONLY: max_dim_name

    INTEGER, INTENT(in)                   :: g_nland       ! number of global landpoints
    INTEGER, INTENT(in)                   :: l_nland       ! number of landpoints in domain
    INTEGER, INTENT(in)                   :: zntiles       ! The number of plant functional types to be used by the phenology module
    LOGICAL, INTENT(in)                   :: isRestart
    INTEGER, INTENT(in)                   :: fileformat    ! output file format (grib/netcdf)
    INTEGER, INTENT(in)                   :: nlct          ! Number of Land classes
    TYPE(t_stream), POINTER, OPTIONAL     :: stream

  !local declarations
    TYPE(t_stream),POINTER      :: phenStream
    INTEGER                     :: dim2p(2), dim2(2)
    CHARACTER(LEN=max_dim_name) :: dim2n(2)

   !In CCDAS-mode we check here, wheter nlct is consistent with what we assumed in mo_init_ccdas.
#ifdef _CCDAS
   IF (nlct .ne. SIZE(pkl)) CALL finish("init_phenology","nlct from lclibt does not correspond to assumptin in init_ccdas")
#endif
    
 
    ntiles = zntiles
    laimmult = exp (-1._dp/30._dp)
    tmpmmult = exp (-1._dp/30._dp)


    IF( mod(86400,INT(delta_time)) /= 0) THEN ! For computing the mean day temperature and others it is assumed that each day ..
                                              ! is computed with a fixed number of time steps. Therefore the program is  ..
                                              ! stopped wheen day length is not an integer multiple of the time step.
       CALL finish("init_phenology","ERROR: Day length is not an integer multiple of the time step!")
    ELSE
       time_steps_per_day = 86400/int(delta_time)
    END IF

        ! --- define phenology stream 

    IF (PRESENT(stream)) THEN
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'pheno', filetype=fileformat, interval=io_time_event(1,'day','first',0))
          ! Set default stream options
          CALL default_stream_setting(stream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       phenStream => stream
    ELSE
       ! Add new stream
       CALL new_stream(phenStream, 'pheno', filetype=fileformat, interval=io_time_event(1,'day','first',0))
       ! Set default stream options
       CALL default_stream_setting(phenStream, table=land_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    dim2p = (/ l_nland, zntiles /)
    dim2  = (/ g_nland, zntiles /)
    dim2n(1) = 'landpoint'
    dim2n(2) = 'tiles'

    ! code numbers of this routine are 185, 187, 194-197 and 205-208 using the GRIB veg_table
    CALL add(phenStream,'laim',                 laim,        longname='Water limited LAI memory', &
             units='',                         ldims=dim2p, gdims=dim2, dimnames=dim2n, code=194,  &
             lpost=.TRUE.)
    laim=0._dp !plaimax !Initialization
    CALL add(phenStream,'tmpm',                 tmpm,        longname='Temperature memory', &
             units='K',                         ldims=dim2p, gdims=dim2, dimnames=dim2n, code=195,  &
             lpost=.TRUE.)
    tmpm=273._dp !Initialization
    CALL add(phenStream,'temp_last',          temp_last, longname='Air Temperature of the Previous Day', &
             units='K', code=185 )
    temp_last=0._dp
    CALL add(phenStream,'temp_sum',          temp_sum, longname='Sum of Air Temperature since Midnight', &
             units='K', code=187, lpost=.FALSE.)
     temp_sum=0._dp
    CALL add(phenStream,'previous_day_lai',  lai_last, longname='LAI of the Previous Day', &
             units='', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=205 )
    lai_last=0._dp  !vg: plaimax not yet available
    CALL add(phenStream,'root_moist_daily_mean',root_moisture_mean, longname='Daily average of plant available water', &
             units='m', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=196 )
    root_moisture_mean=0._dp
    CALL add(phenStream,'root_moist_daily_sum',root_moisture_sum, longname='Sum of plant available water since Midnight', &
             units='m', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=197,lpost=.FALSE. )
    root_moisture_sum=0._dp
    CALL add(phenStream,'lai_sum',  lai_sum, longname='Sum of LAI since Midnight', &
             units='', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=207, lpost=.FALSE.)
    lai_sum=0_dp
    CALL add(phenStream,'previous_day_evap_pot',  evap_pot_last, longname='Potential Evapotranspiration of the Pervious Day', &
             units='?', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=206 )
    evap_pot_last=0._dp
    CALL add(phenStream,'evap_pot_sum',  evap_pot_sum, longname='Sum of Potential Evapotranspiration since Midnight', &
             units='?', ldims=dim2p, gdims=dim2, dimnames=dim2n, code=208, lpost=.TRUE.)
    evap_pot_sum=0._dp


     ! If this is a restart run we can exit now since the model state variables are read from restart file
    IF (isRestart) THEN
       RETURN
    ENDIF
  END SUBROUTINE init_phenology_ccdas


  ! --- update_phenology() ---------------------------------------------------------------------------------------------------------
  ! 
  ! This is the main routine of the phenology module. 
  !
  ! Remark: This routine has to be called every time step. In case the whole grid is not processed in a single call, it may be 
  ! called several times during one time step --- in that case the interface variables do not cover all landpoints associated 
  ! with the particular processor, but only only those of a single block of the processor domain. Therefore, the routine has to 
  ! be called with the appropriate section of the processor fields. Example: Let  procArray(1:domain%nland) be a field covering all 
  ! landpoints a particular processor handles (e.g. the LAI). The block currently processed by the processor is 
  ! procArray(kstart:kend) --- this is what has to be handed over to update_phenology().
  !
  SUBROUTINE update_phenology_ccdas (kidx,lctlib,lai,soil_rel_moist,pheno_type_of_tile,pft_type,tmp_air, &
                                 sinlat,coslat,rad_longwave_down,rad_shortwave_net,rad_shortwave_down, &
                                surface_temp,apar,par,pressure,wind10,veg_height,gc0,laimax_logrop, &
                                 spec_humidity)

  INTEGER,INTENT(IN)               :: kidx                        !the number of elemtns in the current call
  TYPE(lctlib_type),INTENT(IN)     :: lctlib
  REAL(dp),INTENT(INOUT),TARGET    :: lai(:,:)                    !lai to be newly computed
  REAL(dp),INTENT(IN)              :: soil_rel_moist(:,:)         !! relative soil moisture within the root zone
  INTEGER, INTENT(IN)              :: pheno_type_of_tile(:,:)     !! The phenology types of the tiles of grid points (including 
                                                                  !!    "unvegetated");  Dimension: (1:nidx,1:ntiles)
  INTEGER, INTENT(IN)              :: pft_type(:,:)               !! PFT-type
  REAL(dp),INTENT(IN)              :: tmp_air(:)                  !! Air temperature at time-step (for ECHAM in the lowest layer)
  REAL(dp),INTENT(IN)              :: sinlat(:),coslat(:)
  REAL(dp),INTENT(IN)              :: rad_longwave_down(:),rad_shortwave_net(:),rad_shortwave_down(:) !! (kidx)
  REAL(dp),INTENT(IN)              :: surface_temp(:,:)              !! Surface temperature (kidx,ntiles)
  REAL(dp),INTENT(IN)              :: apar(:,:),par(:,:)
  REAL(dp),INTENT(IN)              :: pressure(:)                  !!surface air pressure [Pa]
  REAL(dp),INTENT(IN)              ::wind10(:)
  REAL(dp),INTENT(IN)              :: veg_height(:,:)  !!vegetation height
  REAL(dp),INTENT(IN)              :: gc0(:,:)  ! unstressed canopy conductance
  REAL(dp),INTENT(IN)              :: laimax_logrop(:,:) !Laimax from definition file
  REAL(dp),INTENT(IN)              :: spec_humidity(:)

  ! local variables

  INTEGER               :: ni,nk,k,kidx0
  REAL(dp)              :: laimax                       ! holds the water-limited LAI      
  REAL(dp)              :: root_moisture(kidx,ntiles)                !The plant-available soil moisture
  REAL(dp)              :: r,lailim,ft,fd
  REAL(dp)              :: daylen                       !length of the day.
  INTEGER           :: day_of_month = -1                     !! day in the current month (in the sense of date)
  INTEGER,SAVE      :: day_of_month_at_prev_ts = -1          !! day_of_month at previous time step
  REAL(dp)          :: transpiration_pot(kidx,ntiles)                     !! Potential transpiration
  REAL(dp)          :: fapar(kidx,ntiles)


  ! --- set block range indices
  kidx0 = kstart    ! first value of kindex() (the index of the first element of a block in the processor domain)
    CALL get_date_components(current_date,DAY=day_of_month)
    CALL get_date_components(previous_date, DAY=day_of_month_at_prev_ts)
    first_ts_of_day = ( day_of_month /= day_of_month_at_prev_ts )

    !we need fpar
    WHERE (par .gt. 0._dp)
       fapar=apar/par
    ELSEWHERE
       fapar=0._dp
    END WHERE
    !! We now have to calculate the potential transpiration rate per time-step.
    !! pressure*spec_humidity/(Eps-spec_humidity*Eps+spec_humidity) gives vapor_pressure that is not passed to the interface 
    !! and hence needs to be recalculated from spec_humidity. Calculation is inverse of the one in mo_jsbalone_forcing.f90:2106
    CALL calc_transpiration_pot(kidx,transpiration_pot,tmp_air,lai,rad_longwave_down,rad_shortwave_net, &
                                   rad_shortwave_down,surface_temp,fapar, &
                                   pressure,pressure*spec_humidity/(Eps-spec_humidity*Eps+spec_humidity), &
                                   wind10,veg_height,gc0,pheno_type_of_tile)

    root_moisture(:,:) = soil_rel_moist(1:kidx,:)

    CALL update_day(tmp_air,transpiration_pot,lai,root_moisture)
   IF (.not. first_ts_of_day) THEN
     ! phenology is only calculated once a day. 
     RETURN
   ENDIF
   DO nk=1,kidx
     k=kidx0+nk-1
     CALL calc_daylen(sinlat(nk),coslat(nk),daylen)
     DO ni=1,ntiles
         ! Vegetation types as used by W. Knorr are not given explicitely.
         ! Instead phenology types are used, but parameters are landclass
         ! dependent (is called here pft_type) 
         SELECT CASE(pheno_type_of_tile(nk,ni))
            CASE(unvegetated)
                lai(nk,ni)=0.
            CASE(raingreen)
            ! PFT 1,2,3,7 (14lctlib)  warm-evergreen and warm-deciduous corresponds to
            ! raingreen in JSBACH
                plaimax=laimax_logrop(nk,ni)
                laimax=root_moisture_mean(k,ni)*1000._dp*lai_last(k,ni) / lctlib%knorr_Tau_w(pft_type(nk,ni))/&
                       maxx(evap_pot_last(k,ni),1e-3_dp,1e-2_dp)
                laimax = mins (laimax, plaimax, 0.99_dp)
                laim(k,ni) = laimax * (1._dp - laimmult) + laim(k,ni) * laimmult
                r = lctlib%knorr_leaf_growth_rate(pft_type(nk,ni))
                lailim=laim(k,ni)
                lai(nk,ni)=lailim-(lailim-lai(nk,ni))*exp(-r)
            CASE(evergreen,summergreen)
            ! Here comes the switch for PFT's Numbers 4,5,6,8,11 (14lctlib)  of W.Knorr 
            ! (cold-evergreen and cold-deciduous. cooresponds to evergreen and
            ! summergreen with the exceptoin of tundra (nr. 11) which is phenotype grasses -
                plaimax=laimax_logrop(nk,ni)
              tmpm(k,ni)=temp_last(k)*(1_dp-tmpmmult) + tmpm(k,ni)*tmpmmult
              if (lctlib%knorr_T_phi(pft_type(nk,ni)).le.0) call finish("update_phenology()",&
                "We never should end here, wrong assignement of p_phi to PFT (GJS)")
              ft=errf((tmpm(k,ni) - Tmelt - lctlib%knorr_T_phi(pft_type(nk,ni)))/lctlib%knorr_T_r(pft_type(nk,ni)))
              if (lctlib%knorr_Day_r(pft_type(nk,ni)).le.0) call finish("update_phenology()",&
                "We never should end here, wrong assignement of t_r to PFT (GJS)")
              fd=errf((daylen - lctlib%knorr_Day_c(pft_type(nk,ni)))/lctlib%knorr_Day_r(pft_type(nk,ni)))
              r=ft*fd*lctlib%knorr_leaf_growth_rate(pft_type(nk,ni))+(1._dp-ft*fd)*lctlib%knorr_k_l(pft_type(nk,ni))+1e-9_dp
              lailim = maxx(ft*fd*lctlib%knorr_leaf_growth_rate(pft_type(nk,ni))*plaimax/r,1e-9_dp,1e-3_dp)
              lai(nk,ni)=lailim-(lailim-lai(nk,ni))*exp(-r)
            CASE(crops,grasses)
            ! And all the rest
                plaimax=laimax_logrop(nk,ni)
              tmpm(k,ni)=temp_last(k)*(1_dp-tmpmmult) + tmpm(k,ni)*tmpmmult
              if (lctlib%knorr_T_phi(pft_type(nk,ni)).le.0) call finish("update_phenology()",&
               "We never should end here, wrong assignement of p_phi to PFT (GJS)")
              ft=errf((tmpm(k,ni) - Tmelt - lctlib%knorr_T_phi(pft_type(nk,ni)))/lctlib%knorr_T_r(pft_type(nk,ni)))
                laimax=root_moisture_mean(k,ni)*1000._dp*lai_last(k,ni) / lctlib%knorr_Tau_w(pft_type(nk,ni))/&
                maxx(evap_pot_last(k,ni),1e-3_dp,1e-2_dp)
                laimax = mins (laimax, plaimax, 0.99_dp)
                laim(k,ni) = laimax * (1._dp - laimmult) + laim(k,ni) * laimmult
                r=ft*lctlib%knorr_leaf_growth_rate(pft_type(nk,ni))+(1._dp-ft)*lctlib%knorr_k_l(pft_type(nk,ni))+1e-9_dp
                lailim=maxx(ft*lctlib%knorr_leaf_growth_rate(pft_type(nk,ni))*laim(k,ni)/r,1e-9_dp,1e-3_dp)
                lai(nk,ni)=lailim-(lailim-lai(nk,ni)) * exp(-r)
            CASE default
                   ! We should never end here, This seems to be a not
                   ! implemented phenology type
                   CALL finish("update_phenology()","There is a phenology type that has not been considered")
            END SELECT
    ENDDO
   ENDDO

   RETURN
  END SUBROUTINE update_phenology_ccdas

  SUBROUTINE calc_daylen(sinlat,coslat,daylen)
  USE mo_jsbach_constants, ONLY: pi
  USE mo_orbit,            ONLY: inquire_declination

  ! Calculates daylength 
  REAL(dp),INTENT(IN)::sinlat,coslat
  REAL(dp),INTENT(OUT)::daylen
  !local variables
  REAL(dp) ::spds,cpds,arg
  REAL(dp)::declination


  CALL inquire_declination(declination)

  spds=sinlat*SIN(declination)
  cpds=coslat*COS(declination)
  arg=-spds/cpds
  IF (arg > 1._dp) THEN
    daylen=0._dp
  ELSE IF(arg < -1._dp) THEN
    daylen=1._dp
  ELSE
    daylen=ACOS(arg)/pi
  ENDIF
    daylen=daylen*24._dp

  RETURN
  END SUBROUTINE calc_daylen

  SUBROUTINE update_day(air_temp,evap_pot,lai,root_moisture)
  ! This is to get daily values for the daily time-step.
  REAL(dp),INTENT(IN)::air_temp(1:nidx)
  REAL(dp),INTENT(IN)::evap_pot(1:nidx,1:ntiles)
  REAL(dp),INTENT(IN)::lai(1:nidx,1:ntiles)
  REAL(dp),INTENT(IN)::root_moisture(1:nidx,1:ntiles)
  INTEGER :: kidx0,kidx1
  kidx0=kstart
  kidx1=kend
  IF (first_ts_of_day) THEN  ! Starting with a new day and preparing daily variables for Phenology
   temp_last(kidx0:kidx1)=temp_sum(kidx0:kidx1)/time_steps_per_day
   temp_sum(kidx0:kidx1)=air_temp(1:nidx)
   lai_last(kidx0:kidx1,1:ntiles)=lai_sum(kidx0:kidx1,1:ntiles)/time_steps_per_day
   lai_sum(kidx0:kidx1,1:ntiles)=lai(1:nidx,1:ntiles)
   evap_pot_last(kidx0:kidx1,1:ntiles)=evap_pot_sum(kidx0:kidx1,1:ntiles)   !I just need the sum
   evap_pot_sum(kidx0:kidx1,1:ntiles)=evap_pot(1:nidx,1:ntiles)*24_dp/time_steps_per_day*60_dp*60_dp
   root_moisture_mean(kidx0:kidx1,1:ntiles)=root_moisture_sum(kidx0:kidx1,1:ntiles)/time_steps_per_day
   root_moisture_sum(kidx0:kidx1,1:ntiles)=root_moisture(1:nidx,1:ntiles)
  ELSE ! during a day just summing up
   temp_sum(kidx0:kidx1)=temp_sum(kidx0:kidx1)+air_temp(1:nidx)
   lai_sum(kidx0:kidx1,1:ntiles)=lai_sum(kidx0:kidx1,1:ntiles)+lai(1:nidx,1:ntiles)
   evap_pot_sum(kidx0:kidx1,1:ntiles)=evap_pot_sum(kidx0:kidx1,1:ntiles)+evap_pot(1:nidx,1:ntiles)*24_dp/&
   time_steps_per_day*60_dp*60_dp
            ! the factors at the end is to bring a unit per second onto a "time-step sum"
   root_moisture_sum(kidx0:kidx1,1:ntiles)=root_moisture_sum(kidx0:kidx1,1:ntiles)+root_moisture(1:nidx,1:ntiles)
  ENDIF
  RETURN
  END SUBROUTINE update_day

  SUBROUTINE calc_transpiration_pot(kidx,transpiration_pot,temp_air,lai,rad_longwave_down,rad_shortwave_net,  &
                                      rad_shortwave_down,surface_temp,fapar,pressure,vapor_pressure,wind10,veg_height,gc0,&
                                      pheno_type_of_tile)
  !!!GJS most parts of this routine have been adapted from surface subroutine as used in CCDAS from W.Knorr  
  !!!GJS description in  W.Knorr 1997 ("Satellite Remote Sensing and Modelling of the Global CO2 Exchange of Land Vegetation:
  !!!GJS                                A Synthesis Study"

    USE mo_bethy_constants,       ONLY: FcMax, FcMin, LaiLimit

  INTEGER,INTENT(IN)  :: kidx
  REAL(dp),INTENT(OUT):: transpiration_pot(:,:)  !Potential transpiration  (kidx,ntiles)
  REAL(dp),INTENT(IN) :: temp_air(:)              !Air temperature in [K] (kidx)
  REAL(dp),INTENT(IN) :: lai(:,:)                ! Leaf area index (kidx,ntiles)
  REAL(dp),INTENT(IN) :: rad_longwave_down(:),rad_shortwave_net(:),rad_shortwave_down(:) !(kidx)
  REAL(dp),INTENT(IN) :: surface_temp(:,:)        !surface temperature (kidx,ntiles)
  REAL(dp),INTENT(IN) :: fapar(:,:) ,pressure(:),vapor_pressure(:),wind10(:)
  REAL(dp),INTENT(IN) :: veg_height(:,:)
  REAL(dp),INTENT(IN) :: gc0(:,:) !unstressed canopy conductance
  INTEGER, INTENT(IN) :: pheno_type_of_tile(:,:)     !! The phenology types of the tiles of grid points (including "unvegetated")

  !local
  REAL(dp)  :: slope(kidx)   !slope of water vapour pressure curven
  REAL(dp)  :: temp_degree(kidx)      !Temperature in degree centigrade
  REAL(dp)  :: sw_snow(kidx)          !Switch to account for Snow/NoSnow
  REAL(dp)  :: a0(kidx),b0(kidx)
  REAL(dp)  :: sat_press(kidx)
  REAL(dp)  :: fc_veg(kidx,ntiles)       ! fractional_vegetation_cover
  INTEGER   :: i,j
  REAL(dp)  :: transm(kidx,ntiles)      !transmissivity of vegetation
  REAL(dp)  :: rad_veg_net(kidx,ntiles)            !vegetation radiation
  REAL(dp)  :: rad_longwave_net(kidx,ntiles)            !vegetation radiation
  REAL(dp)  :: soil_heat_flux(kidx,ntiles)
  REAL(dp)  :: veg_absorb(kidx,ntiles)        !absorbivity of vegetation
  REAL(dp)  :: alv,as0        ! constants to calculate veg_absorb
  REAL(dp)  :: rho(kidx)      !density of air
  REAL(dp)  :: vap_deficit(kidx)
  REAL(dp)  :: href,rz0 ! reference height for wind speed. 
  REAL(dp)  :: ga(kidx,ntiles)
  REAL(dp)  :: rlam0,rlam1,rlam2,lambda(kidx) !lambda and constants to calculate it
  REAL(dp)  :: gamma_psy(kidx)  !psychrometric constant


  alv=0.15_dp
  as0=0.05_dp

  transpiration_pot=0._dp

  temp_degree = temp_air - Tmelt

  !! vapor saturatoin pressure and slope of the curve
    sw_snow=SIGN(0.5_dp,temp_degree)+0.5_dp
    a0 = 17.269_dp * sw_snow+ 22.33_dp * (1._dp-sw_snow)
    b0 = 237.3_dp * sw_snow + 271.15_dp * (1._dp-sw_snow)
    sat_press = 610.78_dp*EXP(a0*temp_degree/(b0+temp_degree))
    slope = a0*b0*sat_press/(b0+temp_degree)**2
  !! Radiation balance of vegetation
     !!First we need fractional vegetation cover
      DO i=1,kidx
       DO j=1,ntiles
         fc_veg(i,j) = lai(i,j)/LaiLimit * FcMax
         fc_veg(i,j) = minx(fc_veg(i,j),FcMax,1e-3_dp)
         fc_veg(i,j) = maxx(fc_veg(i,j), FcMin,FcMin*10._dp)
       ENDDO
      ENDDO
     !!transmissivity
     transm =  (1._dp - fc_veg) +  fc_veg * EXP( -0.5_dp * lai / fc_veg )
     !! Net longwave radiation and soil heat flux
     DO i=1,ntiles
           !!!GJS surface_temp seems to be the same for all tiles  
             rad_longwave_net(:,i)=rad_longwave_down(:) - Emissivity * StefanBoltzmann * surface_temp(:,i)**4
             soil_heat_flux(:,i) = 0.036_dp*(rad_shortwave_net(:) + rad_longwave_net(:,i))
     ENDDO
     !! absorbed shortwave radiation
     veg_absorb=(1._dp-alv-as0)*fapar
     DO i=1,ntiles
        rad_veg_net(:,i)=(1._dp-transm(:,i))*(rad_longwave_net(:,i)-soil_heat_flux(:,i))+rad_shortwave_down(:)*veg_absorb(:,i)
     ENDDO

     !!air density
     rho=molarMassDryAir_kg/UniversalGasConstant/temp_air*pressure
     !!Water vapor deficit
     vap_deficit=(sat_press-vapor_pressure     )
     !! Aerodynamic conductancce
     href=10._dp
     rz0=0.1_dp
     DO i=1,ntiles
      !!!GJS An altarnative would be to use a global minum veg_heigth..
       WHERE (pheno_type_of_tile(:,i) .ne. unvegetated)
        ga(:,i) = vonKarman**2._dp*wind10(:)/log(href/(rz0*veg_height(:,i))+1._dp)**2._dp
       ELSEWHERE
        ga(:,i) = 0._dp
       END WHERE
     ENDDO
     !!lambda
     rlam0=3.1512E6_dp
     rlam1=-2.38E3_dp
     rlam2=3.834E6_dp
     lambda=( rlam0 + rlam1 * temp_air ) * sw_snow + rlam2 * ( 1._dp - sw_snow )
     !!Psychrometric constant
     gamma_psy=pressure*SpecificHeatDryAirConstPressure/eps/lambda

     !! And the final calculation to get potential transpiration
     DO i=1,ntiles
         !WHERE(gc0(:,i) .ne. 0._dp)
         WHERE (pheno_type_of_tile(:,i) .ne. unvegetated)
          transpiration_pot(:,i)=(slope(:)*rad_veg_net(:,i)+rho(:)*SpecificHeatDryAirConstPressure*vap_deficit(:)*ga(:,i))/ &
                                 !( slope(:)+gamma_psy(:)*(1._dp+ga(:,i)/max(gc0(:,i),0.0000001_dp)))/ &
                                 ( slope(:)+gamma_psy(:)*(1._dp+ga(:,i)/gc0(:,i)))/ &
                                 lambda(:)
         ELSEWHERE
          transpiration_pot(:,i)=0._dp
         END WHERE
     ENDDO
    RETURN
  END SUBROUTINE calc_transpiration_pot
  
  SUBROUTINE reset_phenology_ccdas()
   !This routine resets the values of variables initialized in init to the init-values
   ! (GJS) not entirely sure how to deal with phenStream (and streams in general)  
   ! (GJS) I am pretty sure that lai in update has not to be reseted, because this should be done elsewhere
   laim=0._dp !plaimax 
   tmpm=273._dp
   temp_last=0._dp
   lai_last=plaimax
   root_moisture_mean=0._dp
   root_moisture_sum=0._dp
   evap_pot_last=0._dp
   temp_sum=0._dp
   lai_sum=0_dp
   evap_pot_sum=0._dp
   first_ts_of_day=.FALSE.

   RETURN
  END SUBROUTINE


END MODULE mo_phenology_ccdas
