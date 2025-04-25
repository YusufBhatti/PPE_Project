!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
module mo_cbal_cpools
    USE mo_kind,             ONLY: dp
    USE mo_cbal_parameters,  ONLY: cn_green, cn_woods, cn_litter_green, cn_litter_wood, cn_slow
    USE mo_land_surface,     ONLY: fract_small
    USE mo_jsbach_constants, ONLY: Tmelt

    implicit none

    public :: update_Cpools
    public :: N_process
    public :: relocate_CarbonAndNitrogen
    public :: C_relocation_from_LUtransitions
    public :: C_relocation_from_harvest
    public :: relocate_carbon_desert
    public :: relocate_carbon_fire
    public :: relocate_carbon_damage
    public :: printCpoolParameters
    public :: printNpoolParameters
    private :: yasso

    REAL(dp), Parameter,public  :: g_C_perMol            = 12.011_dp ! self explaining, or ? 
    REAL(dp), Parameter,public  :: g_N_perMol            = 14.0067_dp

    REAL(dp), parameter,public  :: frac_wood_aboveGround  = 0.7_dp  ! Fraction of Carbon above ground in wood pool (for separation
                                                                    !    of woody litter into above and below ground litter pools)
    REAL(dp), parameter,public  :: frac_green_aboveGround = 0.5_dp  ! Fraction of Carbon above ground in green pool (for separation
                                                                    !    of green litter into above and below ground litter pools)
   

    ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART =============================================================================

    private

    REAL(dp), parameter :: Q10                    = 1.8_dp       ! Empirical parameter for temperature dependence of heterotrophic
                                                                 !    respiration rate
    REAL(dp), parameter :: kappa                  = 1.0_dp       ! Empirical parameter for soil moisture dependence of heterotrophic
                                                                 !    respiration rate
    REAL(dp), parameter :: referenceTemp_Q10      = Tmelt            ! Reference temperature in the Q10 formula [Kelvin]
    REAL(dp), parameter :: tau_Cpool_reserve      = 1.0_dp *365._dp  ! time constant by which reserve pool is depreciated [days]

    REAL(dp), parameter :: tau_Cpool_slow         = 100._dp*365._dp  ! time constant by which slow pool is depreciated  [days]
    REAL(dp), parameter :: greenC2leafC           = 4.0_dp           ! ratio of carbon green pool (leaves, fine roots, starches,
                                                                     !    sugars) to leaf carbon
    REAL(dp), parameter :: frac_C_litter_wood2atmos = 0.2_dp     ! fraction of heterotrophic loss of woody litter emitted to the
                                                                 !    atmosphere
    REAL(dp), parameter :: frac_C_faeces2_LG      = 0.3_dp       ! fraction of C from faces put into the litter green pool
    REAL(dp), parameter :: alpha_critical         = 0.35_dp      ! critical value of relative soil moisture below which
                                                                 !    heterotrophic respiration almost stops
    REAL(dp), parameter :: alpha_min              = 0.1_dp       ! effective lowest value for soil moisture entering the
                                                                 !    heterotrophic respiration
    REAL(dp), parameter :: N2O_rate_nitrification = 0.001_dp     ! JSBACH assumes lower bound of estimation in the literture as: 
                                                                 !   (0.1%) due to nitrification 
    REAL(dp), parameter :: N2O_rate_denitrification = 0.00125_dp ! 0.0025_dp  ! & (0.2%) due to denitrification 
                                                                 ! actually 0.2% of denitrified N, not of N inputs
                                                                 ! assumption: 50% loss of NO3 by leaching
    REAL(dp), parameter :: N2O_rate               = 0.005_dp     ! JSBACH assumes 0.5-1%
                                                                 !    N2O emissions rate per day (<0.1%-0.2% nitrification ); 
                                                                 !    0.2-4.7% due to denit in Xu Ri, 2008 * therein
    REAL(dp), parameter :: sminn_NH4_fraction     = 0.4_dp       ! soil mineral NH4 is 40% of SMINN pool and rest in the form of NO3
                                                                 ! Ref: Xu Ri, 2008

    REAL(dp), parameter :: sec_per_day            = 86400._dp    ! seconds per day
    REAL(dp), parameter :: days_per_year          = 365.25_dp
    REAL(dp), parameter :: sec_per_year           = days_per_year*sec_per_day

contains

  ! --- printCpoolParameters() -----------------------------------------------------------------------------------------------
  !
  ! Prints out the parameters used for the carbon pool model
  !
  SUBROUTINE printCpoolParameters
    USE mo_exception,  ONLY: message, message_text

    WRITE (message_text,*) "=== CpoolParameters ========================================="
    CALL message("",message_text)
    WRITE (message_text,*) "                     Q10=", Q10
    CALL message("",message_text)
    WRITE (message_text,*) "                   kappa=", kappa
    CALL message("",message_text)
    WRITE (message_text,*) "       referenceTemp_Q10=", referenceTemp_Q10, " Kelvin"
    CALL message("",message_text)
    WRITE (message_text,*) "       tau_Cpool_reserve=", tau_Cpool_reserve/365., " years"
    CALL message("",message_text)
    WRITE (message_text,*) "          tau_Cpool_slow=", tau_Cpool_slow/365., " years"
    CALL message("",message_text)
    WRITE (message_text,*) "            greenC2leafC=", greenC2leafC
    CALL message("",message_text)
    WRITE (message_text,*) "frac_C_litter_wood2atmos=", frac_C_litter_wood2atmos
    CALL message("",message_text)
    WRITE (message_text,*) "  frac_green_aboveGround=", frac_green_aboveGround
    CALL message("",message_text)
    WRITE (message_text,*) "   frac_wood_aboveGround=", frac_wood_aboveGround
    CALL message("",message_text)
    WRITE (message_text,*) "       frac_C_faeces2_LG=", frac_C_faeces2_LG
    CALL message("",message_text)
    WRITE (message_text,*) "          alpha_critical=", alpha_critical
    CALL message("",message_text)
    WRITE (message_text,*) "               alpha_min=", alpha_min
    CALL message("",message_text)
    WRITE (message_text,*) "============================================================="
    CALL message("",message_text)

  END SUBROUTINE printCpoolParameters

 ! --- printNpoolParameters() ------------------------------------------------------------------------------------------------------
  !
  ! Prints out the parameters used for the nitrogen pool model
  !
  SUBROUTINE printNpoolParameters
    USE mo_exception,       ONLY: message, message_text

    WRITE (message_text,*) "=== NpoolParameters ========================================="
    CALL message("",message_text)
    WRITE (message_text,*) "        cn_green=", cn_green
    CALL message("",message_text)
    WRITE (message_text,*) "        cn_woods=", cn_woods
    CALL message("",message_text)
    WRITE (message_text,*) " cn_litter_green=", cn_litter_green
    CALL message("",message_text)
    WRITE (message_text,*) "  cn_litter_wood=", cn_litter_wood
    CALL message("",message_text)
    WRITE (message_text,*) "         cn_slow=", cn_slow
    CALL message("",message_text)
    WRITE (message_text,*) "============================================================="

  END SUBROUTINE printNpoolParameters 


  ! --- update_Cpools() -----------------------------------------------------------------------------------------------------------
  !
  ! Updates carbon and nitrogen pools and computes soil respiration rate and Net Ecosystem Product [mol(C)/m^2 s].
  ! Has to be called once each day for each covertype, where carbon and nitrogen pools exists. 
  !

  elemental subroutine update_Cpools(LAI, LAI_previous, NPPrate, topSoilTemp, alpha,                          &
                                          frac_npp_2_woodPool,frac_npp_2_reservePool, frac_npp_2_exudates,    &
                                          frac_green_2_herbivory,                                             &
                                          tau_Cpool_litter_green,tau_Cpool_litter_wood,tau_Cpool_woods,       &
                                          LAI_shed_constant,                                                  &
                                          frac_C_litter_green2atmos,Max_C_content_woods,                      &
                                          specific_leaf_area_C, reserveC2leafC,                               &
                                          max_LAI, is_vegetation, Cpool_green, Cpool_woods, Cpool_reserve,    &
                                          Cpool_litter_green_ag,Cpool_litter_green_bg,                        &
                                          Cpool_litter_wood_ag,Cpool_litter_wood_bg,                          & 
                                          Cpool_slow,                                                         &
                                          soilResp_rate, soilResp_rate_pot,                                   &
                                          NPP_flux_correction, excess_NPP,root_exudates,                      &
                                          Cflx_litterTotal,                                                   &
                                          Cflx_herbivory,Cflx_herbivory_LG, Cflx_herbivory_2_atm,             &
                                          NPP_act, frac_litter_wood_new,                                      &
                                          ! Yasso variables
                                          temp2_30d, precip_30d,                                              &
                                          YCpool_acid_ag1, YCpool_water_ag1, YCpool_ethanol_ag1,              &
                                          YCpool_nonsoluble_ag1,                                              & 
                                          YCpool_acid_bg1, YCpool_water_bg1, YCpool_ethanol_bg1,              &
                                          YCpool_nonsoluble_bg1, YCpool_humus_1,                              & 
                                          YCpool_acid_ag2, YCpool_water_ag2, YCpool_ethanol_ag2,              &
                                          YCpool_nonsoluble_ag2,                                              & 
                                          YCpool_acid_bg2, YCpool_water_bg2, YCpool_ethanol_bg2,              &
                                          YCpool_nonsoluble_bg2, YCpool_humus_2,                              & 
                                          LeafLit_coef1, LeafLit_coef2, LeafLit_coef3, LeafLit_coef4,         &
                                          LeafLit_coef5, WoodLit_coef1, WoodLit_coef2,                        &
                                          WoodLit_coef3, WoodLit_coef4, WoodLit_coef5,                        &
                                          WoodLitterSize,                                                     &
                                          ! Nitrogen variables
                                          redFact_Nlimit,                                                     &
                                          Npool_green, Npool_woods, Npool_mobile,                             &
                                          Npool_litter_green_ag,Npool_litter_green_bg,                        &
                                          Npool_litter_wood_ag,Npool_litter_wood_bg,                          &
                                          Npool_slow, SMINN_pool,                                             &
                                          minNflux_litter_green_ag,minNflux_litter_green_bg,                  &
                                          minNflux_litter_wood_ag,minNflux_litter_wood_bg,                    &
                                          minNflux_slow,Nplant_demand,Nsoil_demand,                           &
                                          Ntotal_demand,SMINN_herbivory,N2O_emissions_mineraliz,              &      
                                          N2O_emissions_slow,                                                 &
                                          N2O_emissions_grazing                                               &
                                           ) 

    real(dp),intent(in)    :: LAI                    !! Yesterdays mean LAI
    real(dp),intent(in)    :: LAI_previous           !! The day before yesterdays mean LAI
    real(dp),intent(in)    :: NPPrate                !! Yesterdays mean NPPrate [mol(C)/m^2 s]
    real(dp),intent(in)    :: topSoilTemp            !! Yesterdays mean temperature of upper soil layer [degree Kelvin]
    real(dp),intent(in)    :: alpha                  !! Yesterdays mean water stress factor (between 0 and 1)
    real(dp),intent(in)    :: frac_npp_2_woodPool    !! Fraction of NPP to be put maximally into the green pool
    real(dp),intent(in)    :: frac_npp_2_reservePool !! Fraction of NPP to be put into the optimally into the reserve pool
    real(dp),intent(in)    :: frac_npp_2_exudates    !! Fraction of NPP to be put into the optimally into the root exudates
    real(dp),intent(in)    :: frac_green_2_herbivory
    real(dp),intent(in)    :: tau_Cpool_litter_green !! Time constant by which green litter pool is depreciated  [days]
    real(dp),intent(in)    :: tau_Cpool_litter_wood  !! Time constant by which woody litter pool is depreciated  [days]
    real(dp),intent(in)    :: tau_Cpool_woods        !! Time constant by which woods Pool is depreciated  [years]
    real(dp),intent(in)    :: LAI_shed_constant      !! Leaf shedding at a constant rate for evergreens etc.
    real(dp),intent(in)    :: frac_C_litter_green2atmos !! Fraction of heterotrophic loss of green litter pools emitted to the 
                                                        !!    atmosphere
    real(dp),intent(in)    :: Max_C_content_woods    !! Maximum carbon content of wood pool (from lctLib)
    real(dp),intent(in)    :: specific_leaf_area_C   !! Specific leaf area (from lctLib) [m^2(leaf)/mol(Carbon)]
    real(dp),intent(in)    :: reserveC2leafC         !! Ratio of max. carbon in reserve pool to max. carbon of leaves
    real(dp),intent(in)    :: max_LAI                !! Maximum value of LAI
    logical, intent(in)    :: is_vegetation          !! logical vegetation mask
    real(dp),intent(inout) :: Cpool_green            !! Green carbon pool: on input last value; updated on output 
                                                     !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods            !! Wood carbon pool: on input last value; updated on output
                                                     !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve          !! Reserve carbon pool: on input last value; updated on output
                                                     !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_green_ag  !! Above ground green litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_green_bg  !! Below ground green litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood_ag   !! Above ground woody litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood_bg   !! Below ground woody litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_slow             !! Slow soil carbon pool: on input last value; updated on output
                                                     !!    [mol(C)/m^2(canopy)]
    real(dp),intent(out)   :: soilResp_rate          !! Soil (=heterotrophic) respiration rate  [mol(C)/m^2 s] Note: this is a loss,
                                                     !!    therefore negative!
    real(dp),intent(out)   :: soilResp_rate_pot      !! soil respiration without N limitation     
    real(dp),intent(out)   :: NPP_flux_correction    !! Amount by which the NPP rate entering the routine has to be corrected. This
                                                     !!    correction arises either because otherwise the reserve pool would get
                                                     !!    negative (positive correction), or the wood pool would exceed its
                                                     !!    maximum value (negative correction). [mol(C)/m^2 s]
    real(dp),intent(out)   :: excess_NPP             !! Part of NPP that could not be stored in one of the plant carbon pools 
                                                     !!    (green, wood, reserve), but had to be  thrown away (a posteriori 
                                                     !!    reduction of NPP) [mol(C)/m^2 s]
    real(dp),intent(out)   :: root_exudates          !! Value of NPP_2_rootExudates had to be dropped into the litter_green_bg pool.
                                                     !!    [mol(C)/m^2 s]
    real(dp),intent(out)   :: Cflx_litterTotal       !! Total carbon flux from vegetation to litter (green+woody, ag+bg)
    real(dp),intent(out)   :: Cflx_herbivory         !!
    real(dp),intent(out)   :: Cflx_herbivory_LG      !!
    real(dp),intent(out)   :: Cflx_herbivory_2_atm   !!
    real(dp),intent(out)   :: NPP_act                !! Actual NPP after N-limitation and excess carbon drop [mol(C)/m^2 s]

    ! added for thonicke fire algorithm (spitfire) needed to correctly split the carbon pools into the different fuel classes
    real(dp),optional,intent(inout) :: frac_litter_wood_new   !! new fraction in above ground wood litter pool  

    ! Meteorology for Yasso
    real(dp),optional,intent(in)    :: temp2_30d              !! 30 day mean temperature
    real(dp),optional,intent(in)    :: precip_30d             !! 30 day mean precipitation

    ! Yasso pools
    ! Size class 1; green litter
    !   Aboveground  
    real(dp),optional,intent(inout) :: YCpool_acid_ag1        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_ag1       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_ag1     !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_ag1  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground  
    real(dp),optional,intent(inout) :: YCpool_acid_bg1        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_bg1       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_bg1     !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_bg1  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_humus_1         !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2; woody litter
    !   Aboveground  
    real(dp),optional,intent(inout) :: YCpool_acid_ag2        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_ag2       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_ag2     !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_ag2  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground  
    real(dp),optional,intent(inout) :: YCpool_acid_bg2        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_bg2       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_bg2     !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_bg2  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_humus_2         !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Parameters for Yasso
       ! Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    real(dp),optional,intent(in)    :: LeafLit_coef1         !! fraction going to the acid soluble pools 
    real(dp),optional,intent(in)    :: LeafLit_coef2         !! fraction going to the water soluble pools
    real(dp),optional,intent(in)    :: LeafLit_coef3         !! fraction going to the ethanol soluble pools
    real(dp),optional,intent(in)    :: LeafLit_coef4         !! fraction going to the non soluble pools
    real(dp),optional,intent(in)    :: LeafLit_coef5         !! fraction going to the humus pool
       ! Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    real(dp),optional,intent(in)    :: WoodLit_coef1         !! fraction going to the acid soluble pools 
    real(dp),optional,intent(in)    :: WoodLit_coef2         !! fraction going to the water soluble pools
    real(dp),optional,intent(in)    :: WoodLit_coef3         !! fraction going to the ethanol soluble pools
    real(dp),optional,intent(in)    :: WoodLit_coef4         !! fraction going to the non soluble pools
    real(dp),optional,intent(in)    :: WoodLit_coef5         !! fraction going to the humus pool
    real(dp),optional,intent(in)    :: WoodLitterSize        !! size of coarse debris


    !! NOTE: all optional arguments are needed with runs with nitrogen
    real(dp),optional,intent(out)   :: redFact_Nlimit         !! Reduction factor for NPP due to N-limitation [1]
    real(dp),optional,intent(inout) :: Npool_green            !! Green N-pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_woods            !! Woods N-pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_mobile           !! Mobile plant N-pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_litter_green_ag  !! Above ground green litter N-pool: updated on output
                                                              !!    [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_litter_green_bg  !! Below ground green litter N-pool: updated on output
                                                              !!    [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_litter_wood_ag   !! Above ground woody litter N-pool: updated on output
                                                              !!    [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_litter_wood_bg   !! Below ground woody litter N-pool: updated on output
                                                              !!    [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_slow             !! Slow soil organic N-pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: SMINN_pool             !! Soil mineral nitrogen pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(out)   :: minNflux_litter_green_ag !! Uptake of mineral N from above ground green litter pool
                                                                !!    [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: minNflux_litter_green_bg !! Uptake of mineral N from below ground green litter pool
                                                                !!    [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: minNflux_litter_wood_ag  !! Uptake of mineral N from above ground wood litter pool
                                                                !!    [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: minNflux_litter_wood_bg  !! Uptake of mineral N from below ground wood litter pool
                                                                !!    [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: minNflux_slow          !! Uptake of mineral N from slow soil pool  [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: Nplant_demand          !! Uptake of N for plant from soil [mol(N)/m^2 s] 
    real(dp),optional,intent(out)   :: Nsoil_demand           !! Uptake of N for soil microbes [mol(N)/m^2 s] 
    real(dp),optional,intent(out)   :: Ntotal_demand          !! Total N demand for both plant and soil [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: SMINN_herbivory        !! SMINN gained from herbivory  [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: N2O_emissions_mineraliz  !! Nitros oxide emission [mol(N)/m^2(canopy)s] due to mineralization
    real(dp),optional,intent(out)   :: N2O_emissions_slow     !! N2O emission [mol(N)/m^2(canopy)s] due to mineralization from
                                                              !!      slow pool
    real(dp),optional,intent(out)   :: N2O_emissions_grazing  !! N2O emission [mol(N)/m^2(canopy)s] due to N input by herbivores


    ! locals

    real(dp) :: d_Cpool_slow
    real(dp) :: d_Cpool_litter_green_ag, d_Cpool_litter_green_bg
    real(dp) :: d_Cpool_litter_wood_ag, d_Cpool_litter_wood_bg
    real(dp) :: Cpool_green_max,Cpool_reserve_optimal
    real(dp) :: NPP_2_greenPool,NPP_2_woodPool,NPP_2_reservePool,NPP_2_rootExudates,Green_2_herbivory
    real(dp) :: C_2_litter_greenPool,C_2_litter_woodPool
    real(dp) :: excess_carbon
    real(dp) :: leaf_shedding
    real(dp) :: alpha_mod

    real(dp) :: Cpool_green_pot, Cpool_woods_pot         !! values of potential carbon pools i.e. before accounting for N-limitation
    real(dp) :: Cpool_slow_pot, Cpool_reserve_pot
    real(dp) :: Cpool_litter_green_ag_pot, Cpool_litter_green_bg_pot
    real(dp) :: Cpool_litter_wood_ag_pot, Cpool_litter_wood_bg_pot

    real(dp) :: Cflx_NPP_2_green_pot,Cflx_NPP_2_wood_pot,Cflx_NPP_2_reserve_pot 
    real(dp) :: Cflx_litterGreen_ag_2_atmos_pot, Cflx_litterGreen_bg_2_atmos_pot
    real(dp) :: Cflx_litterGreen_ag_2_slow_pot, Cflx_litterGreen_bg_2_slow_pot
    real(dp) :: Cflx_litter_wood_ag_2_atmos_pot,Cflx_litter_wood_bg_2_atmos_pot
    real(dp) :: Cflx_litter_wood_ag_2_slow_pot,Cflx_litter_wood_bg_2_slow_pot
    real(dp) :: Cflx_slow_2_atmos
    real(dp) :: Cflx_faeces_2_LG, Cflx_faeces_2_atm
    real(dp) :: litter_wood_new

    ! local variables for Yasso
    real(dp), dimension(18)  :: Yasso_io_pools
    real(dp), dimension(2)  :: Weather
    real(dp), dimension(10) :: Yasso_out
    real(dp), dimension(5)  :: LeafLit_coefV
    real(dp), dimension(5)  :: WoodLit_coefV   
    real(dp)                :: leafLitter      ! amount of leafLitter 
    real(dp)                :: sizeZero

    real(dp) :: Nflux_NPP_2_green_pot,Nflux_NPP_2_wood_pot 
    real(dp) :: Nflux_litterGreen_ag_2_slow_pot, Nflux_litterGreen_bg_2_slow_pot
    real(dp) :: Nflux_litter_wood_ag_2_slow_pot,Nflux_litter_wood_bg_2_slow_pot
    real(dp) :: minNflux_litter_green_ag_pot,minNflux_litter_green_bg_pot
    real(dp) :: minNflux_litter_wood_ag_pot,minNflux_litter_wood_bg_pot
    real(dp) :: N2O_emissions_litter_pot 

    real(dp) :: minN_plant_demand,minN_soil_demand,total_minN_demand
    real(dp) :: decompRate_litter_green

    real(dp) :: NPP_2_green_via_mobileN, NPP_2_wood_via_mobileN
    real(dp) :: Nmobile_2_green, Nmobile_2_wood, Nmobile_2_green_pot, Nmobile_2_wood_pot
    real(dp) :: redFact_Nmobile, redFact_Nlimitation

    LOGICAL  :: with_Nitrogen        !! .FALSE.: only Carbon pools are updated; .TRUE.: in addition also Nitrogen pools are updated 
    LOGICAL  :: with_yasso           !! .FALSE.: cbalance model for litter and soil carbon is used; .TRUE.: yasso model is used 
    REAL(dp), PARAMETER :: N2O_ef_grazing = 0.0125_dp               !analogue to fertilizer ef


    !! Determine whether call is with/without nitrogen
    IF (PRESENT(Npool_green)) THEN
       with_Nitrogen = .true.
    ELSE
       with_Nitrogen = .false.
    END IF

    !! Determine whether call is with/without yasso model
    IF (PRESENT(YCpool_acid_ag1)) THEN
       with_yasso = .true.
    ELSE
       with_yasso = .false.
    END IF



    ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    soilResp_rate            = 0.0_dp
    soilResp_rate_pot        = 0.0_dp
    NPP_act                  = 0.0_dp
    redFact_Nlimitation      = 1.0_dp

    IF (with_Nitrogen) THEN
       redFact_Nlimit           = 1.0_dp
       minNflux_litter_green_ag = 0.0_dp
       minNflux_litter_green_bg = 0.0_dp
       minNflux_litter_wood_ag  = 0.0_dp
       minNflux_litter_wood_bg  = 0.0_dp
       minNflux_slow            = 0.0_dp
       Nplant_demand            = 0.0_dp
       Nsoil_demand             = 0.0_dp
       Ntotal_demand            = 0.0_dp
       SMINN_herbivory          = 0.0_dp
       N2O_emissions_mineraliz  = 0.0_dp
       N2O_emissions_slow       = 0.0_dp
       N2O_emissions_grazing    = 0.0_dp
    END IF
    IF (PRESENT(frac_litter_wood_new)) THEN
       litter_wood_new = 0.0_dp
    END IF

    ! Preparations
    !-------------------------------------------

    alpha_mod = max(alpha_min,(alpha-alpha_critical)/(1._dp-alpha_critical)) !! modified soil moist. such that heterotrophic
                                                                             !!     respiration is zero below alpha_critical

    ! Initializations
    NPP_flux_correction   = 0.0_dp
    Cpool_green_max       = 0.0_dp
    Cpool_reserve_optimal = 0.0_dp
    excess_NPP            = 0.0_dp
    root_exudates         = 0.0_dp
    Cflx_litterTotal      = 0.0_dp
    Cflx_herbivory        = 0.0_dp
    Cflx_herbivory_LG     = 0.0_dp
    Cflx_herbivory_2_atm  = 0.0_dp

    IF (with_yasso) THEN
       Yasso_io_pools(1)=YCpool_acid_ag1
       Yasso_io_pools(2)=YCpool_water_ag1
       Yasso_io_pools(3)=YCpool_ethanol_ag1
       Yasso_io_pools(4)=YCpool_nonsoluble_ag1
       Yasso_io_pools(5)=YCpool_acid_bg1
       Yasso_io_pools(6)=YCpool_water_bg1
       Yasso_io_pools(7)=YCpool_ethanol_bg1
       Yasso_io_pools(8)=YCpool_nonsoluble_bg1
       Yasso_io_pools(9)=YCpool_humus_1
       Yasso_io_pools(10)=YCpool_acid_ag2
       Yasso_io_pools(11)=YCpool_water_ag2
       Yasso_io_pools(12)=YCpool_ethanol_ag2
       Yasso_io_pools(13)=YCpool_nonsoluble_ag2
       Yasso_io_pools(14)=YCpool_acid_bg2
       Yasso_io_pools(15)=YCpool_water_bg2
       Yasso_io_pools(16)=YCpool_ethanol_bg2
       Yasso_io_pools(17)=YCpool_nonsoluble_bg2
       Yasso_io_pools(18)=YCpool_humus_2

       Weather(1) = temp2_30d
       Weather(2) = precip_30d

       LeafLit_coefV(1) = LeafLit_coef1
       LeafLit_coefV(2) = LeafLit_coef2
       LeafLit_coefV(3) = LeafLit_coef3
       LeafLit_coefV(4) = LeafLit_coef4
       LeafLit_coefV(5) = LeafLit_coef5
       WoodLit_coefV(1) = WoodLit_coef1
       WoodLit_coefV(2) = WoodLit_coef2
       WoodLit_coefV(3) = WoodLit_coef3
       WoodLit_coefV(4) = WoodLit_coef4
       WoodLit_coefV(5) = WoodLit_coef5

       Yasso_out(1:10) = 0.0              ! array to store the Yasso pools/fluxes 
       leafLitter = 0.0                   ! amount of leafLitter 
       sizeZero = 0.0
    END IF

    if (is_vegetation) then      !! only for covertypes with positive specific leaf area computations are meaningful

       Cpool_green_max = greenC2leafC * LAI / specific_leaf_area_C
       Cpool_reserve_optimal = reserveC2leafC * max_LAI / specific_leaf_area_C

       !! ============== NON-N-LIMITED CARBON AND NITROGEN ALLOCATION =========================================================
       !!
       !! Perform those C- and N-allocation steps that are not restricted by N-availabilty
       !! 1. Leaf shedding
       !! 2. Wood shedding
       !! 3. Depletion of the reserve pool
       !! 4. Decomposition of slow soil pool (whether this is really independent of N-availability is debatable!!)
       !!
       !! 
       !! ----------------------------------------------------------------------------
       !! 1. Leaf shedding: Transfer of C from green pool to leaf litter litter pool

!!$ tr
!!$       IF (LAI >= LAI_previous) THEN                              !! If LAI is increasing or constant the leaf shedding is given
!!$          leaf_shedding = LAI * LAI_shed_constant                 !!    by a constant loss of leaves for evergreens, raingreens
!!$                                                                  !!    and grasses.
!!$       ELSE                                                       !! Otherwise
!!$          leaf_shedding = LAI_previous - LAI                      !!    the leaf shedding is given by the decrease in LAI.
!!$       END IF
!!$ tr
       leaf_shedding = MAX(LAI * LAI_shed_constant,LAI_previous - LAI)  !! Leaf shedding is assured at a minimum constant rate 
                                                                        !! which maybe exceeded by the decrease in LAI.       

       C_2_litter_greenPool =  &                                     !! Maximally the whole carbon from the green pool is transfered
          min(greenC2leafC * leaf_shedding / specific_leaf_area_C, & !!    by leaf and root shedding to the green litter pool
          Cpool_green)                                               
       Cpool_green = Cpool_green - C_2_litter_greenPool              !! This removes carbon from the green pool
       excess_carbon = max(0.0_dp,Cpool_green - Cpool_green_max)     !! If by the reduction of LAI the maximum value of the green
                                                                     !!    pool is still smaller than the current value, there is 
       C_2_litter_greenPool = C_2_litter_greenPool + excess_carbon   !!    excess C that has also to be shedded to the green litter
       Cpool_green = Cpool_green - excess_carbon                     !!    pool and also substracted from the green pool
       Cpool_litter_green_ag = Cpool_litter_green_ag + frac_green_aboveGround*C_2_litter_greenPool         !! Move shedded Carbon to
       Cpool_litter_green_bg = Cpool_litter_green_bg + (1._dp-frac_green_aboveGround)*C_2_litter_greenPool !!    green litter pools

       IF (with_Nitrogen) THEN
          Npool_green        = Npool_green - C_2_litter_greenPool/cn_green   !! Shed N according to C/N-ratio of green litter
          Npool_litter_green_ag = &
              Npool_litter_green_ag +         frac_green_aboveGround*C_2_litter_greenPool /cn_litter_green   
          Npool_litter_green_bg = &
              Npool_litter_green_bg + (1._dp-frac_green_aboveGround)*C_2_litter_greenPool /cn_litter_green

          !! retranslocation the surplus N from shedded green parts to mobile N   ("leaf-N retranslocation")             
          Npool_mobile = Npool_mobile + (1._dp/cn_green - 1._dp/cn_litter_green) * C_2_litter_greenPool  !! makes CN constant
                                                                                                         !!    b/w green-litterpool
       END IF

       Cflx_litterTotal=C_2_litter_greenPool

       !! ---------------------------------------------------------------------------------
       !! Herbivory loss from Cpool_green & Npool_green w.r.t PFTs in lctlib

       Green_2_herbivory  = Cpool_green * frac_green_2_herbivory             !! Grazing flux indirect comput.from NPP losses 
       Cpool_green        = Cpool_green - Green_2_herbivory                  !! Loss due to grazing

       Cflx_faeces_2_LG  = Green_2_herbivory * frac_C_faeces2_LG             !! Cflux_faeces_2-LG to be stored in litter green pool
                                                                             !!      (not considered in Cflx_litterTotal) 
       Cflx_faeces_2_atm = Green_2_herbivory * (1.0_dp - frac_C_faeces2_LG)  !! Cflux_faeces_2_atm to be lost to atmosphere 
                                                                             !!      (dk: not added to CO2?)
       Cpool_litter_green_ag = Cpool_litter_green_ag + Cflx_faeces_2_LG      !! Cflux_faeces_2_LG put into LG pool  

       IF (with_Nitrogen) THEN
          ! new: corrections for use of cn_litter_green with Cflx_2_atm
          Npool_green = max(0.0_dp,Npool_green - Green_2_herbivory/cn_green) 
          Npool_litter_green_ag = Npool_litter_green_ag + Cflx_faeces_2_LG /cn_litter_green
          !! assumed that there is no N lossess to atm and infact all goes to the animals 
          !! faeces N according to intended  C/N-ration of green litter 
          !! animal faeces are assigned to green_litter_ag
          !! animals minus CO2/CH4 flux to atmosphere: urine, assigned to sminn_pool 

          !! retranslocation the surplus N from grazed animals faeces (urine) to SMINN pool 
          SMINN_pool = SMINN_pool + (1._dp/cn_green - 1._dp/cn_litter_green) * Cflx_faeces_2_LG    &
               + 1._dp/cn_green * Cflx_faeces_2_atm     
          !!cn of herbivory = cn_green as there is no N loss assumed

          !! N2O emissions from mineral N released by grazing; first assumption: EF similar to fertilizer application
          N2O_emissions_grazing = N2O_ef_grazing * ( (1._dp/cn_green - 1._dp/cn_litter_green) * Cflx_faeces_2_LG  &
                                                 + 1._dp/cn_green * Cflx_faeces_2_atm )   

          !!DIAGNOSTIC OUTPUTS: SMINN gain from herbivory faeces and dungs 
          SMINN_herbivory = SMINN_herbivory + (1._dp/cn_green - 1._dp/cn_litter_green) * Cflx_faeces_2_LG    &
               + 1._dp/cn_green * Cflx_faeces_2_atm  - N2O_emissions_grazing 
          SMINN_herbivory = SMINN_herbivory / sec_per_day
       END IF

       !! ---------------------------------------------------------------------------------
       !! 2. Wood shedding: Transfer of C from wood pool to wood litter pools
       !! 
                                                            !! Assuming that forests continously die at the inverse lifetime of trees
       C_2_litter_woodPool = Cpool_woods / tau_Cpool_woods !*365._dp     !! .. the shedded wood (MAX FUNCTION FOR NONWOODY PFTS)
       Cpool_litter_wood_ag = Cpool_litter_wood_ag &                    !! .. is put partly into the above ground woody litter pool
            + frac_wood_aboveGround*C_2_litter_woodPool
       if (present(frac_litter_wood_new)) &
          litter_wood_new = frac_wood_aboveGround * C_2_litter_woodPool  
       Cpool_litter_wood_bg = Cpool_litter_wood_bg &                    !! .. and the rest into the below ground woody litter pool
            + (1.0_dp-frac_wood_aboveGround)*C_2_litter_woodPool
       Cpool_woods = Cpool_woods -  C_2_litter_woodPool                 !! .. and then all is substracted from the wood pool.

       Cflx_litterTotal = Cflx_litterTotal + C_2_litter_woodPool

       if(with_Nitrogen) then 
          !! Handle the associated organic N-fluxes
          Npool_woods = Npool_woods - C_2_litter_woodPool / cn_woods            !! N removed from wood pool..
          Npool_litter_wood_ag = Npool_litter_wood_ag &                         !! .. is partly put to the ag woody N litter pool..
               + frac_wood_aboveGround*C_2_litter_woodPool /cn_litter_wood
          Npool_litter_wood_bg = Npool_litter_wood_bg &                         !! .. and the rest to the bg woody N litter pool..
               + (1.0_dp-frac_wood_aboveGround)*C_2_litter_woodPool /cn_litter_wood
          SMINN_pool = SMINN_pool &                                             !! .. while the soil mineral N-pool gains the 
               + (1._dp/cn_woods - 1._dp/cn_litter_wood) * C_2_litter_woodPool  !!    N-difference
          !! ASSUMPTION: cn_woods <= cn_litter_wood
       end if

       !! ----------------------------------------------------------------------------
       !! 3. Depletion of reserve pool: Transfer of C from reserve pool to green litter litter pool
       !!    Some organic carbon of the reserve pool is always lost with mortality of plants. 
       !!    Note that the reserve pool contains no Nitrogen (starches and sugar are free of N).
       !!    Note also that the green litter pool has no fixed C/N-ratio (therefore no compensation
       !!    flux is needed)

       C_2_litter_greenPool = Cpool_reserve / tau_Cpool_reserve   
       Cpool_litter_green_ag = Cpool_litter_green_ag + frac_green_aboveGround*C_2_litter_greenPool
       Cpool_litter_green_bg = Cpool_litter_green_bg + (1._dp-frac_green_aboveGround)*C_2_litter_greenPool
       Cpool_reserve = Cpool_reserve - C_2_litter_greenPool 

       Cflx_litterTotal = Cflx_litterTotal + C_2_litter_greenPool

       !! ----------------------------------------------------------------------------
       !! 4. Decomposition of slow soil pool (i) C-flux (ii) Mineral N-flux
       !!
       d_Cpool_slow = &
            -MIN(alpha_mod**kappa * q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  / & !! Slow soil respiration
            tau_Cpool_slow,1._dp)  * Cpool_slow                                          !! .. according to Q10-model
       Cpool_slow = Cpool_slow + d_Cpool_slow                 !! Remove respired Carbon (note the negative sign!)

       Cflx_slow_2_atmos            = -d_Cpool_slow          !! Remember the decomposition flux to the atmosphere (CO2-release)

       if(with_Nitrogen) then
          minNflux_slow = Cflx_slow_2_atmos/cn_slow             !! Mineral N-flux from slow pool to soil mineral N-pool
          !! part of the minNflux_slow flux emits as N2O gas by considering implicit NH4 & NO3 availability in soil       
          N2O_emissions_slow = (minNflux_slow * sminn_NH4_fraction * N2O_rate_nitrification              & 
                             + minNflux_slow * (1._dp - sminn_NH4_fraction) * N2O_rate_denitrification * alpha) &
                             * MAX(0._dp,MIN(topSoilTemp/30._dp,1._dp))
                                                                    ! alpha for soil moisture regulation of denitrification
                                                                    !!  N2O emissions due to nitrification  & denitrification 
                                                                    !!  a small fraction of N loss (N2O) to atm.
                                                                    !! .. but NO is ignored here     
!          N2O_emissions_slow = minNflux_slow * sminn_NH4_fraction * N2O_rate_nitrification              & 
!                             + minNflux_slow * (1._dp - sminn_NH4_fraction) * N2O_rate_denitrification 
!                                                                   !!  N2O emissions due to nitrification  & denitrification 
                                                                    !!  a small fraction of N loss (N2O) to atm.
                                                                    !! .. but NO is ignored here
                                                                    !! dk: all fractions are constants, no need to seperate between
                                                                    !!     Nh4 and NO3
          minNflux_slow = minNflux_slow - N2O_emissions_slow        !! Remained N-flux after N2O loss to atm.
          SMINN_pool = SMINN_pool + minNflux_slow                !! Put released mineral N to SMINN_pool

          !! Handle the associated organic N-flux
          Npool_slow = max(0.0_dp,Npool_slow + d_Cpool_slow / cn_slow)  !! check
       end if

       !!
       !! ============== START OF POTENTIAL PLANT-CARBON ALLOCATION TO ESTIMATE PLANT-N-DEMAND =================================
       !!                       (exception: case NPP<0 handles actual allocation)

       Cpool_green_pot   = Cpool_green
       Cpool_woods_pot   = Cpool_woods
       Cpool_reserve_pot = Cpool_reserve

       !! Preliminary determination of distribution of NPP to the various pools (may be corrected afterwards according to 
       !!    avaialable Nitrogen)

       if(NPPrate <= 0.0_dp) then                                      !! In case of negative or zero NPP
          NPP_2_reservePool =  NPPrate  * sec_per_day                  !! .. the plant looses C and  
          Cpool_reserve     = Cpool_reserve + NPP_2_reservePool        !! .. we try to take it from the reserve pool 
          if(Cpool_reserve < 0.0_dp) then                              !! But if thereby the reserve pool gets negative 
             NPP_2_reservePool = NPP_2_reservePool - Cpool_reserve     !! .. we absorb negative NPP only up to zero reseve pool and 
             NPP_flux_correction = -Cpool_reserve / sec_per_day        !! .. give the deficit as a flux correction back to the
             Cpool_reserve = 0.0_dp                                    !! .. calling routine. Hence all Carbon from the reserve pool
          end if                                                       !! .. is used up.

          !! For correct handling of nitrogen limitation below (which handles only the case of positive NPP) set NPP to zero:
          NPP_2_reservePool = 0.0_dp
          NPP_2_greenPool   = 0.0_dp   !! All other transfer rates are zero
          NPP_2_woodPool    = 0.0_dp
          NPP_2_rootExudates= 0.0_dp
          Cpool_reserve_pot = Cpool_reserve      !! remember modified reserve pool


       else  !! NPP is positive                                                            

          NPP_2_woodPool    = frac_npp_2_woodPool * NPPrate * sec_per_day    !! NPP it is distributed according to predefined 
                                                                             !!    relative fraction
          NPP_2_reservePool = frac_npp_2_reservePool * NPPrate * sec_per_day
          NPP_2_rootExudates= frac_npp_2_exudates * NPPrate * sec_per_day
          NPP_2_greenPool   = (1.0_dp - frac_npp_2_woodPool - frac_npp_2_reservePool                 & 
               - frac_npp_2_exudates) * NPPrate * sec_per_day

          !! Growth of Wood Pool (structural carbon of living plants)

          Cpool_woods_pot = Cpool_woods_pot + NPP_2_woodPool                 !! Then it is attempted to put all NPP share into it
          excess_carbon = max(0.0_dp,Cpool_woods_pot-Max_C_content_woods)    !! .. but thereby the pool may get too large
          Cpool_woods_pot = Cpool_woods_pot - excess_carbon                  !! .. so that the excess carbon has once more to be 
                                                                             !! .. subtracted
          NPP_2_greenPool = NPP_2_greenPool + excess_carbon                  !! .. and instead made available to the green pool
          NPP_2_woodPool = NPP_2_woodPool - excess_carbon                    !! .. and the actual amount of NPP put to the wood pool
                                                                             !! .. is less


          !! Growth of Reserve Pool (sugars and starches) and determination how much will enter the Green Pool
          if(Cpool_reserve_pot < Cpool_reserve_optimal) then                 !! .. If the reserve pool is smaller than optimal,
             Cpool_reserve_pot = Cpool_reserve_pot +  NPP_2_reservePool      !! .... it is filled by the available NPP
             excess_carbon =                                         &       !! .... Thereby it may happen that it gets larger than
                  max(0.0_dp,Cpool_reserve_pot - Cpool_reserve_optimal)      !! .... optimal so that there is excess carbon
             Cpool_reserve_pot = Cpool_reserve_pot - excess_carbon           !! .... that needs not be taken up
             NPP_2_greenPool = NPP_2_greenPool + excess_carbon               !! .... but can better be used to increase the green
             NPP_2_reservePool = NPP_2_reservePool - excess_carbon           !! .... pool and the actual amount of NPP put to the
                                                                             !! .... reserve pool is less.
          else                                                               !! .... Otherwise (reserve pool is larger than optimal)
             NPP_2_greenPool = NPP_2_greenPool + NPP_2_reservePool           !! .... all NPP is left for the green pool
             NPP_2_reservePool = 0.0_dp                                      !! .... so that nothing is stored in the reserve pool.
          end if

          !! Growth of Green Pool (leaves and fine roots): (in case of too much NPP, try to put it into the reserve pool).

          Cpool_green_pot = Cpool_green_pot + NPP_2_greenPool                !! Green pool is filled by the available NPP.
          excess_carbon = max(0.0_dp,Cpool_green_pot - Cpool_green_max)      !! .. Thereby it may get larger than appropriate for 
                                                                             !!    current LAI.
          NPP_2_greenPool = NPP_2_greenPool - excess_carbon                  !! .. Hence the actual amount of NPP put to the green
                                                                             !!    pool is less.
          Cpool_green_pot = Cpool_green_pot - excess_carbon                  !! .. This excess carbon needs not be taken up but
          if(Cpool_reserve_pot < Cpool_reserve_optimal) then                 !! .... if the reserve pool is smaller than optimal
             Cpool_reserve_pot = Cpool_reserve_pot + excess_carbon           !! .... it is tried to put the carbon there,
             NPP_2_reservePool = NPP_2_reservePool + excess_carbon           !! .... which means that additional NPP is put to the 
                                                                             !!      reserve pool.
             excess_carbon =                                          &      !! .... Thereby it may happen that the reserve pool
                  max(0.0_dp,Cpool_reserve_pot - Cpool_reserve_optimal)      !!      increases beyond the optimal value.
             Cpool_reserve_pot = Cpool_reserve_pot - excess_carbon           !! .... In that case the excess carbon is once more
                                                                             !!      removed from the reserve pool
             NPP_2_reservePool = NPP_2_reservePool - excess_carbon           !! .... so that the actual amount of NPP put into the
                                                                             !!      reserve pool is less.
          end if                                                             !!

          Cpool_litter_green_bg = Cpool_litter_green_bg + NPP_2_rootExudates
          ! THT 30.10.2012: This taken into account in Yasso too.. Put NPP_2_rootExudates to Yasso-water-pool. Added 20.11.2012
          ! Yasso_io_pools(2) = Yasso_io_pools(2) + NPP_2_rootExudates ! DSG: exudates goes into yasso

          excess_NPP = excess_carbon / sec_per_day

          root_exudates = NPP_2_rootExudates / sec_per_day

       end if !! NPP > 0 end


       !! ============== START OF POTENTIAL LITTER-CARBON ALLOCATION ==============================================

       Cpool_litter_green_ag_pot   = Cpool_litter_green_ag
       Cpool_litter_green_bg_pot   = Cpool_litter_green_bg
       Cpool_litter_wood_ag_pot    = Cpool_litter_wood_ag
       Cpool_litter_wood_bg_pot    = Cpool_litter_wood_bg

       ! Update Green Litter Pools (litter from dead leaves, fruits, debris from bark and fine roots)

       decompRate_litter_green =  &                                     !! decomposition rate of green litter according to Q10 model
            -MIN(alpha_mod**kappa * q10**((topSoilTemp - referenceTemp_Q10)/10._dp)/ tau_Cpool_litter_green,1._dp)

       d_Cpool_litter_green_ag      = decompRate_litter_green * Cpool_litter_green_ag_pot     !! Decomposition of green plant litter
       Cpool_litter_green_ag_pot    = Cpool_litter_green_ag_pot + d_Cpool_litter_green_ag
       d_Cpool_litter_green_bg      = decompRate_litter_green * Cpool_litter_green_bg_pot     !! Decomposition of green plant litter
       Cpool_litter_green_bg_pot    = Cpool_litter_green_bg_pot + d_Cpool_litter_green_bg

       ! Update Wood Litter Pools (litter from all dead woody material, above and below ground)

       d_Cpool_litter_wood_ag = -MIN(q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  / &   !! Decomposition of woody litter
            tau_Cpool_litter_wood,1._dp) * Cpool_litter_wood_ag_pot                         !! .. according to Q10-model
       d_Cpool_litter_wood_bg = -MIN(q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  / &   !! Decomposition of woody litter
            tau_Cpool_litter_wood,1._dp) * Cpool_litter_wood_bg_pot                         !! .. according to Q10-model
       !! .. without soil moisture dependence
       Cpool_litter_wood_ag_pot = Cpool_litter_wood_ag_pot + d_Cpool_litter_wood_ag
       Cpool_litter_wood_bg_pot = Cpool_litter_wood_bg_pot + d_Cpool_litter_wood_bg


       !! ============== START OF POTENTIAL ALLOCATION OF SLOW SOIL POOL ==============================================


       ! Update Slow Soil Pool (uptake of C from litter pools, decomposition has already been handled above)
       Cpool_slow_pot = Cpool_slow  
       Cpool_slow_pot = Cpool_slow_pot                                                           &
            - (1._dp - frac_C_litter_wood2atmos) * (d_Cpool_litter_wood_ag + d_Cpool_litter_wood_bg)  &
            - (1._dp - frac_C_litter_green2atmos) * (d_Cpool_litter_green_ag + d_Cpool_litter_green_bg)           



       !! ============== use nitrogen from N mobile pool to satisfy (at least partly) for potential carbon fluxes ===
       ! Remarks: Nmobile_2_green + Nmobile_2_wood <= Npool_mobile
       ! Npool_mobile accounts minN_plant_demand by trasfering N into green and wood pool 

       Cflx_NPP_2_green_pot          = NPP_2_greenPool                                       !! potential growth of green carbon 
       Cflx_NPP_2_wood_pot           = NPP_2_woodPool                                        !! potential growth of wood carbon 

       NPP_2_green_via_mobileN = 0.0_dp        
       NPP_2_wood_via_mobileN = 0.0_dp 

       IF (with_Nitrogen) THEN

          Nmobile_2_green_pot = Cflx_NPP_2_green_pot/cn_green
          Nmobile_2_wood_pot  = Cflx_NPP_2_wood_pot /cn_woods 

          !... CN mode 
          if (NPPrate > 0.0_dp) then                                   !! (Npool_mobile> 0) & avoid Npool_mobile negative 
             if(Nmobile_2_green_pot + Nmobile_2_wood_pot <= Npool_mobile) then
                Nmobile_2_green = Nmobile_2_green_pot
                Nmobile_2_wood  = Nmobile_2_wood_pot
                NPP_2_green_via_mobileN= Cflx_NPP_2_green_pot         
                NPP_2_wood_via_mobileN = Cflx_NPP_2_wood_pot    
             else
                redFact_Nmobile = Npool_mobile/(Nmobile_2_green_pot + Nmobile_2_wood_pot)  !! avoid Npool_mobile < 0; NPPrate > 0
                Nmobile_2_green = redFact_Nmobile * Nmobile_2_green_pot
                Nmobile_2_wood  = redFact_Nmobile * Nmobile_2_wood_pot
                NPP_2_green_via_mobileN= Nmobile_2_green * cn_green
                NPP_2_wood_via_mobileN = Nmobile_2_wood  * cn_woods
             end if

             !! update C-N green and wood pools 
             Cpool_green = Cpool_green + NPP_2_green_via_mobileN
             Cpool_woods = Cpool_woods + NPP_2_wood_via_mobileN
             Npool_green = Npool_green + Nmobile_2_green
             Npool_woods = Npool_woods + Nmobile_2_wood

             Npool_mobile = max(0.0_dp, Npool_mobile - Nmobile_2_green - Nmobile_2_wood)           !!  update Npool_mobile
          end if
       END IF

       Cflx_NPP_2_green_pot = Cflx_NPP_2_green_pot - NPP_2_green_via_mobileN
       Cflx_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot  - NPP_2_wood_via_mobileN

       Nflux_NPP_2_green_pot = Cflx_NPP_2_green_pot/cn_green
       Nflux_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot/cn_woods


       !! ========= COMPUTE MINERAL-PLANT N DEMAND =================

       minN_plant_demand = Nflux_NPP_2_green_pot + Nflux_NPP_2_wood_pot 


       !! ============== COLLECT POTENTIAL C-FLUXES ===============================================================

       Cflx_NPP_2_reserve_pot        = NPP_2_reservePool                                     !! potential growth of reserve carbon
       Cflx_litterGreen_ag_2_atmos_pot = -frac_C_litter_green2atmos * d_Cpool_litter_green_ag !! ag green litter decomp. -> atm.
       Cflx_litterGreen_bg_2_atmos_pot = -frac_C_litter_green2atmos * d_Cpool_litter_green_bg !! bg green litter decomp. -> atm.
       Cflx_litterGreen_ag_2_slow_pot  = &
           - (1._dp -frac_C_litter_green2atmos) * d_Cpool_litter_green_ag                    !! ag green litter decomp. -> slow pool
       Cflx_litterGreen_bg_2_slow_pot  = &
            -(1._dp -frac_C_litter_green2atmos) * d_Cpool_litter_green_bg                    !! bg green litter decomp. -> slow pool
       Cflx_litter_wood_ag_2_atmos_pot  = -frac_C_litter_wood2atmos * d_Cpool_litter_wood_ag !! ag-wood litter decomp.  -> atm.
       Cflx_litter_wood_bg_2_atmos_pot  = -frac_C_litter_wood2atmos * d_Cpool_litter_wood_bg !! bg-wood litter decomp.  -> atm.
       Cflx_litter_wood_ag_2_slow_pot   = &
            -(1._dp - frac_C_litter_wood2atmos) * d_Cpool_litter_wood_ag                     !! ag-wood litter decomp.  -> slow pool
       Cflx_litter_wood_bg_2_slow_pot   = &
            -(1._dp - frac_C_litter_wood2atmos) * d_Cpool_litter_wood_bg                     !! bg_wood litter decomp.  -> slow pool

       IF (.NOT. with_Nitrogen) THEN

          redFact_Nlimitation = 1.0_dp   !! This assures absence of nitrogen limitation in this mixed C-N code

       ELSE

          !! ============== DERIVE ASSOCIATED POTENTIAL FLUXES OF ORGANIC N ==============================================
          
          Nflux_litterGreen_ag_2_slow_pot = Cflx_litterGreen_ag_2_slow_pot / cn_slow !! N-content determined by N-sink 
          Nflux_litterGreen_bg_2_slow_pot = Cflx_litterGreen_bg_2_slow_pot / cn_slow !! N-content determined by N-sink 
          Nflux_litter_wood_ag_2_slow_pot  = Cflx_litter_wood_ag_2_slow_pot  / cn_slow !! N-content determined by N-sink 
          Nflux_litter_wood_bg_2_slow_pot  = Cflx_litter_wood_bg_2_slow_pot  / cn_slow !! N-content determined by N-sink 

          !! ============== DETERMINE POTENTIAL FLUXES OF MINERAL N ==============================================
          !
          ! Remark: Mineral-N-fluxes are considerd positive when they increase the SMINN-pool
          !
          ! How to derive the mineral N-fluxes:
          ! -----------------------------------
          !
          ! Mineral fluxes are determined as a kind of residues: assuming different but fixed C/N-ratios for the organic pools these
          ! could not be maintained without compensating Nitrogen fluxes: Therefore the strategy is to compute the mineral N-fluxes
          ! such that they assure N-mass conservation. 
          !
          ! Derivation of mineral N-flux (also called "immobilization flux"):
          !
          ! Consider a pool of organic carbon with N-content N and C-content C. Assume an influx of (organic) C called c1 associated
          ! with a N-influx called n1. Assume further two Carbon outfluxes, an organic one called c2 associated with a N-outflux n2,
          ! and a mineral Carbon outflux called c0.
          ! In addition there is a mineral N outflux m. This mineral flux m is coming from the mineral N-pool that must compensate
          ! for potential C/N-deviations. The aim is to determine m:
          !
          ! The change in C and N is then given by
          !
          !      dC     
          ! (1)  -- = c1 - c2 -c0
          !      dt 
          !
          !      dN   
          ! (2)  -- = n1 - n2 - m
          !      dt   
          !
          ! It is further assumed that the N-to-C ratios nc11 and nc2 for in- and outfluxes, respectively, is fixed:
          !
          ! (3)  n1 = nc1 * c1,   n2 = nc2 * c2
          !
          ! Combining the last 2 equations and solving for the compensating mineral flux gives:
          !
          !                                dN
          ! (4)  m = nc1 * c1 - nc2 * c2 - --
          !                                dt
          !
          ! If the N/C ratio of the considered pool shall be kept constant under the prescribed fluxes c1,c2,n1,n2, i.e. 
          !
          !           N   dN      
          ! (5) nc := - = -- 
          !           C   dC
          !
          ! so that
          !
          !      dN        dC
          ! (6)  -- = nc * --
          !      dt        dt
          !
          ! equations (4), (6) and (1) give
          !
          ! (7)  m = nc1 * c1 - nc2 * c2 - nc*(c1 - c2 -c0) = (nc1 - nc)*c1 + (nc - nc2)*c2 + nc*c0.
          !
          ! If in addition it is known that the process by which the organic C-pool is depleted is determined by the rate equation
          !
          !      dC        
          ! (8)  -- = - r * C
          !      dt        
          !
          ! an alternative equation to (7) can be derived: First, using (5) in (8) shows that N is depleted by the same rate:
          !
          !      dN        
          ! (9)  -- = - r * N
          !      dt        
          !
          ! and entering this into (4) gives the alternative result
          !
          ! (10) m = nc1 * c1 - nc2 * c2 + r*N.
          !
          ! This latter form is especially useful when C and N vary independently so that the N/C-ratio changes in time. (Note that
          ! this is no contradiction to (5), because other processes than the depletion (8) may change N/C-ratio, e.g. influxes of
          ! C and N with different N/C-ratio.)
          !
          ! ---------------------------------------------------------------
          ! Now Compute the potential mineral N fluxes (Note: fluxes into SMINN-pool are positive) 

          minNflux_litter_green_ag_pot =    &      !! ag green litter -> SMINN (use eq. (10), because green litter C/N is not fixed)
               - Nflux_litterGreen_ag_2_slow_pot - decompRate_litter_green * Npool_litter_green_ag

          minNflux_litter_green_bg_pot =    &      !! ag green litter -> SMINN (use eq. (10), because green litter C/N is not fixed)
               - Nflux_litterGreen_bg_2_slow_pot - decompRate_litter_green * Npool_litter_green_bg

          minNflux_litter_wood_ag_pot =    &  !! above ground woody litter -> SMINN (use eq. (7), because woody litter C/N is fixed)
               (1._dp/cn_litter_wood - 1._dp/cn_slow) * Cflx_litter_wood_ag_2_slow_pot &
               + Cflx_litter_wood_ag_2_atmos_pot/cn_litter_wood 

          minNflux_litter_wood_bg_pot =    &  !! below ground woody litter -> SMINN (use eq. (7), because woody litter C/N is fixed)
               (1._dp/cn_litter_wood - 1._dp/cn_slow) * Cflx_litter_wood_bg_2_slow_pot &
               + Cflx_litter_wood_bg_2_atmos_pot/cn_litter_wood

          !! ============== COMPUTE N2O emissions==============================================
          ! --- Nitros oxide emissions(N2O)--N2O_Emissions() due to N mineralization ----
          !! N2O emissions rate per day (<0.1%-0.2% nitrification ); 0.2-4.7% due to denit in Xu Ri, 2008 * therein

          !! minNflux_litter can be negativ, as cn_slow much smaller as cn_litter_wood and release of N in decomp not always 
          !!       sufficient  
          IF ( (minNflux_litter_green_ag_pot + minNflux_litter_green_bg_pot + minNflux_litter_wood_ag_pot &
                   + minNflux_litter_wood_bg_pot) .gt. 0._dp ) THEN
            !!if N is released in litter decomposition

            ! dk: new calculation of N2O from litter:
            !  above ground: treated as fertilizer: same emission factor
            !  below ground: decreasing SOC with soil depth represented by assuming that 50% of SOC are within top soil layers
            !                relevant for N2O production (usually 30 cm); 
            !  plus: moisture dependency of denitrification: alpha
            N2O_emissions_litter_pot = (minNflux_litter_green_ag_pot + 0.5_dp * minNflux_litter_green_bg_pot    &
                                      + minNflux_litter_wood_ag_pot + 0.5_dp * minNflux_litter_wood_bg_pot)     &
                                      * sminn_NH4_fraction * N2O_rate_nitrification * MAX(0._dp,MIN(topSoilTemp/30._dp,1._dp))

            N2O_emissions_litter_pot = N2O_emissions_litter_pot                                                 &
                                     + (minNflux_litter_green_ag_pot + 0.5_dp * minNflux_litter_green_bg_pot    &
                                     + minNflux_litter_wood_ag_pot + 0.5_dp * minNflux_litter_wood_bg_pot)      &
                                     * (1._dp - sminn_NH4_fraction)* N2O_rate_denitrification * alpha           &
                                     * MAX(0._dp,MIN(topSoilTemp/30._dp,1._dp)) !!  N2O emissions due to denitrification
          ELSE
            N2O_emissions_litter_pot = 0._dp
          END IF

          !!
          !! ============== COMPUTE MINERAL-N DEMAND ==============================================
          !!

          minN_soil_demand = -( &   !! "-" because of sign convention for fluxes to SMINN-pool. Note: Slow soil pool only releases N
                                 minNflux_litter_green_ag_pot + minNflux_litter_green_bg_pot  &  
                               + minNflux_litter_wood_ag_pot  + minNflux_litter_wood_bg_pot   &
                              )
          ! demand only if sum(minNflux_litter) negativ, then minN_soil_demand positiv
          minN_soil_demand = minN_soil_demand - N2O_emissions_litter_pot       !! Remained minN_soil_demand after N2O loss to atm.
                                                                               !!dk: N2O_em_litter_pot negativ!!
          !!
          !! ============== DETERMINE REDUCTION FACTOR FROM NITROGEN LIMITATION ========================================
          !!

          total_minN_demand = minN_plant_demand + max(0.0_dp,minN_soil_demand)  &
                             + N2O_emissions_litter_pot
          ! minN_soil_demand is negativ for sminn release = no soil demand
          ! N2O release is a function of N availability
          ! N2O emissions do not determine soil N demand, however, N2O 
          ! will be subtracted from sminn, hence considered in total demand

          if( total_minN_demand  <= SMINN_pool + epsilon(1.0_dp) ) then     !! enough mineral N available
             redFact_Nlimitation = 1.0_dp
             Nplant_demand = max(0.0_dp,minN_plant_demand) 
             Nsoil_demand  = max(0.0_dp,minN_soil_demand)            ! no negativ minN_soil_demand in output
             Ntotal_demand = max(0.0_dp,total_minN_demand)                                 
          else
             redFact_Nlimitation = SMINN_pool / total_minN_demand
             Nplant_demand = max(0.0_dp,minN_plant_demand * redFact_Nlimitation )
             Nsoil_demand  = max(0.0_dp,minN_soil_demand  * redFact_Nlimitation )
             Ntotal_demand = max(0.0_dp,total_minN_demand * redFact_Nlimitation )
          end if
       END IF

       !! ============== UPDATE ALL C-POOLS =======================================================================

       Cpool_green       = Cpool_green       + redFact_Nlimitation * Cflx_NPP_2_green_pot
       Cpool_woods       = Cpool_woods       + redFact_Nlimitation * Cflx_NPP_2_wood_pot
       Cpool_reserve     = Cpool_reserve     +  Cflx_NPP_2_reserve_pot        !!dk removed Nlimit, as no dependent on N

       Cpool_litter_green_ag = Cpool_litter_green_ag       &
            - redFact_Nlimitation * (Cflx_litterGreen_ag_2_atmos_pot + Cflx_litterGreen_ag_2_slow_pot)

       Cpool_litter_green_bg = Cpool_litter_green_bg       &
            - redFact_Nlimitation * (Cflx_litterGreen_bg_2_atmos_pot + Cflx_litterGreen_bg_2_slow_pot)

       ! fraction of new wood litter, assume old and new litter are respired in the same way, respiration does not change the fraction
       if (present(frac_litter_wood_new)) then
          if ((Cpool_litter_wood_ag >  EPSILON(1._dp)) .AND. (litter_wood_new < Cpool_litter_wood_ag)) then
             frac_litter_wood_new = litter_wood_new / Cpool_litter_wood_ag
          else
             frac_litter_wood_new = 1._dp    ! frac_litter_wood_new=1 means fuel fractions will be reset to original values
          endif
       endif
    
       Cpool_litter_wood_ag = Cpool_litter_wood_ag       &
            - redFact_Nlimitation * (Cflx_litter_wood_ag_2_atmos_pot + Cflx_litter_wood_ag_2_slow_pot)

       Cpool_litter_wood_bg = Cpool_litter_wood_bg       &
            - redFact_Nlimitation * (Cflx_litter_wood_bg_2_atmos_pot + Cflx_litter_wood_bg_2_slow_pot)

       Cpool_slow        = Cpool_slow     &
            + redFact_Nlimitation * (  Cflx_litterGreen_ag_2_slow_pot + Cflx_litterGreen_bg_2_slow_pot &
                                     + Cflx_litter_wood_ag_2_slow_pot  + Cflx_litter_wood_bg_2_slow_pot  )

       IF (with_Nitrogen) THEN

          !! ============== UPDATE ORGANIC-N-POOLS =======================================================================
          !!
          !! N-Pools should be related to respective C-pools simply by division with approriate C/N-ratio (except for green litter
          !! pool). But use here more complicated update following the logic of the above allocation scheme to allow for checking
          !! of N-conservation!
          !!

          Npool_green          = Npool_green       + redFact_Nlimitation * Nflux_NPP_2_green_pot

          Npool_woods          = Npool_woods       + redFact_Nlimitation * Nflux_NPP_2_wood_pot

          Npool_litter_green_ag   = Npool_litter_green_ag                              &
               - redFact_Nlimitation * (Nflux_litterGreen_ag_2_slow_pot + minNflux_litter_green_ag_pot)

          Npool_litter_green_bg   = Npool_litter_green_bg                              &
               - redFact_Nlimitation * (Nflux_litterGreen_bg_2_slow_pot + minNflux_litter_green_bg_pot)

          Npool_litter_wood_ag = Npool_litter_wood_ag                              &
               - redFact_Nlimitation * (Nflux_litter_wood_ag_2_slow_pot + minNflux_litter_wood_ag_pot)

          Npool_litter_wood_bg = Npool_litter_wood_bg                              &
               - redFact_Nlimitation * (Nflux_litter_wood_bg_2_slow_pot + minNflux_litter_wood_bg_pot)

          Npool_slow           = Npool_slow       &
              + redFact_Nlimitation *  ( Nflux_litterGreen_ag_2_slow_pot + Nflux_litterGreen_bg_2_slow_pot &
                                       + Nflux_litter_wood_ag_2_slow_pot + Nflux_litter_wood_bg_2_slow_pot ) 

          !! ============== UPDATE MINERAL-N-POOL =======================================================================

          SMINN_pool        = max(SMINN_pool - redFact_Nlimitation * total_minN_demand,0.0_dp)

          !!... Actual N2O_emissions_mineraliz() based on the N-availability
          N2O_emissions_mineraliz =  N2O_emissions_slow     &                           !! Net N2O emissions  
                                   + redFact_Nlimitation * N2O_emissions_litter_pot
       end IF

       !! ============== COMPUTE ACTUAL NPP (is needed outside the routine e.g. to recompute transpiration through stomata) =======

       if(NPPrate > 0.0_dp) then 
          NPP_act = (NPP_2_green_via_mobileN + NPP_2_wood_via_mobileN + NPP_2_rootExudates +               &
               redFact_Nlimitation * (Cflx_NPP_2_green_pot + Cflx_NPP_2_wood_pot) + Cflx_NPP_2_reserve_pot &
               )/sec_per_day
       else
          NPP_act = NPPrate + NPP_flux_correction  !dk: NPPrate negatave, NPP_flux_correction could not be removed from reserve pool
       endif


       !! ============== COMPUTE SOIL RESPIRATION ===================================================================

       IF (.NOT. with_yasso) THEN
          soilResp_rate = -(                                                                              &
               redFact_Nlimitation * (  Cflx_litterGreen_ag_2_atmos_pot + Cflx_litterGreen_bg_2_atmos_pot &
                                      + Cflx_litter_wood_ag_2_atmos_pot + Cflx_litter_wood_bg_2_atmos_pot &
                                    ) + Cflx_slow_2_atmos                                                 &
                           ) / sec_per_day

          soilResp_rate_pot = -( (  Cflx_litterGreen_ag_2_atmos_pot + Cflx_litterGreen_bg_2_atmos_pot    &
                                   + Cflx_litter_wood_ag_2_atmos_pot +  Cflx_litter_wood_bg_2_atmos_pot  &
                                  ) + Cflx_slow_2_atmos                                                  &
                        ) / sec_per_day
  
       ELSE 
       
          leafLitter    = Cflx_litterTotal - C_2_litter_woodPool + Cflx_faeces_2_LG 
          soilResp_rate = 0.
            
          !! call yasso for leaf litter
          CALL yasso(Yasso_io_pools(1:9), Weather, leafLitter, LeafLit_coefV,               & 
                           sizeZero, Yasso_out, frac_green_aboveGround, NPP_2_rootExudates)

          YCpool_acid_ag1         =Yasso_out(1)
          YCpool_water_ag1        =Yasso_out(2)
          YCpool_ethanol_ag1      =Yasso_out(3)
          YCpool_nonsoluble_ag1   =Yasso_out(4)
          YCpool_acid_bg1         =Yasso_out(5)
          YCpool_water_bg1        =Yasso_out(6)
          YCpool_ethanol_bg1      =Yasso_out(7)
          YCpool_nonsoluble_bg1   =Yasso_out(8)
          YCpool_humus_1          =Yasso_out(9)


          !! update respiration with respiration from leaf litter decomposition
          soilResp_rate = Yasso_out(10) / sec_per_day 

          !! call yasso for woody litter
          CALL yasso(Yasso_io_pools(10:18), Weather, C_2_litter_woodPool, WoodLit_coefV,       & 
                           WoodLitterSize, Yasso_out, frac_wood_aboveGround, 0.0_dp) ! no exudates
  
          YCpool_acid_ag2         =Yasso_out(1)
          YCpool_water_ag2        =Yasso_out(2)
          YCpool_ethanol_ag2      =Yasso_out(3)
          YCpool_nonsoluble_ag2   =Yasso_out(4)
          YCpool_acid_bg2         =Yasso_out(5)
          YCpool_water_bg2        =Yasso_out(6)
          YCpool_ethanol_bg2      =Yasso_out(7)
          YCpool_nonsoluble_bg2   =Yasso_out(8)
          YCpool_humus_2          =Yasso_out(9)
       
          !! Yasso total respiration
          ! add respiration from wood litter decomposition
           soilResp_rate = soilResp_rate  + Yasso_out(10) / sec_per_day 

       END IF

       !! ============== COMPUTE Herbivory RESPIRATION =====
        
        Cflx_herbivory_2_atm = -(Cflx_faeces_2_atm) /sec_per_day

        !! ============== Litter and herbivory diagnostics (flux from day-1 to s-1)
        Cflx_litterTotal = Cflx_litterTotal /sec_per_day
        Cflx_herbivory    = Green_2_herbivory / sec_per_day
        Cflx_herbivory_LG = Cflx_faeces_2_LG / sec_per_day   ! DSG: goes into yasso


        IF (with_Nitrogen) THEN

           redFact_Nlimit = redFact_Nlimitation

           !! ============== OTHER DIAGNOSTIC OUTPUTS ===================================================================

          !! redFact_Nlimit was missing first, but minNflux is reduced with N limitation
           minNflux_litter_green_ag   = redFact_Nlimitation * minNflux_litter_green_ag_pot / sec_per_day
           minNflux_litter_green_bg   = redFact_Nlimitation * minNflux_litter_green_bg_pot / sec_per_day
           minNflux_litter_wood_ag    = redFact_Nlimitation * minNflux_litter_wood_ag_pot  / sec_per_day
           minNflux_litter_wood_bg    = redFact_Nlimitation * minNflux_litter_wood_bg_pot  / sec_per_day
           minNflux_slow              = minNflux_slow / sec_per_day                                      ! no N limited
           Nplant_demand = Nplant_demand / sec_per_day                                                   ! demands already N limited
           Nsoil_demand  = Nsoil_demand  / sec_per_day
           Ntotal_demand = Ntotal_demand / sec_per_day
           N2O_emissions_mineraliz =  N2O_emissions_mineraliz / sec_per_day                              ! already N limited 
                                                                                               ! only N2O_emissions_litter_pot	   
           N2O_emissions_slow =  N2O_emissions_slow / sec_per_day                  ! N2O lost from slow to atmosphere; not N limited
           N2O_emissions_grazing = N2O_emissions_grazing / sec_per_day

        end IF

    end if

  end subroutine update_Cpools


  ! --- N_process() ---------------------------------------------------------------------------------------------------  
  !!
  !! ====Biological N Fixation, Deposition, Denitrification(N-forms gases loss) & Leaching===========
  !DESCRIPTION: On the radiation time step, update the nitrogen fixation rate as a function of annual total NPP (60PgC/Y) & ET. 
  !This rate gets updated once per year.! All N fixation goes to the soil mineral N pool. 
  !Forcing Nfix to global constant = 0.4 gN/m2/yr; BNF= 105 TgN/y; Range by Cleveland, 1999 is 100-290 TgN/y. 
  !BNF_NPP = 0.001666*NPP_annual & BNF_ET = 0.008 * (ET -40); unit: 40 cm/y or 400/365*86400 mm/s

  !BNF Calibration: For 100 Tg N/Y; -----------------------------------
  !                                     100 Tg   100 * e12
  !                             factor  = ------ = ---------- = 0.00166
  !                                     60 Pg    60  * e15 
  !converts to BNF in mole ->12(14); calibrated as 12 (to gN)->100Tg N/y =0.00166 -> more in NCAR CLM model codes 
  !->BNF_NPP = 1.8 *(1.0- exp(-0.003 * NPP_annual))/12 

  !                   OR
  ! BNF_NPP = 0.005*NPP (CLM-if NPP is 1PgC); BNF_NPP = 0.0005*NPP (TRIPLEX)
  ! BNF_NPP = 0.005*NPP/365   per day

  !!..Fixation based on ET -> Century Model, Schimel, 1996; Xu-Ri, 2008 
  !->BNF_ET=0.008 * (ET - 400)/(365 * 86400) 

  !....N Deposition (wet + Dry)----------
  ! Global const = 0.5 g/m2/y = 1.5*10-8 g/m2/s (Ref -> Class model, 2007 & Dickison, 2002 
  !->Ndepo= (0.4/(365* 86400))/14  !! g - mole Nitrogen 
  !Reference -> Frank Dentener, 2006 & Galloway et al., 2004 // Inputs for 1860, 1990, 2000, 2030, 2050 
  !Data downloaded from ftp://daac.ornl.gov/data/global_climate/global_N_deposition_maps
  !Global Maps of Atmospheric Nitrogen Deposition, 1860, 1993, and 2050 -> TM3 model: 3-D Chemistry trasport 
  !->Data used: Dentener, 2006 GBC V20 --1860, 2000(S1), 2030(S4-SRES A2) (by Personal communicattion)

  !....Denitrification (gaseous loss of Nitrogen: NO, N2O, N2)--- 
  !!-> Denit = 1.15*SMINN**0.57   -> Ref: Century model, Thomas,2002 per day
  !!-> Denit = dnp*SMINN          -> Ref: CLM model dnp = 0.01 per day
  
  !---Leaching descriptions-----------------
  !! L = (sf*SMINN)*runoff/tot_liquid_water -> Ref:CLM model and also in Century (Thomas, 2002); sf = 0.1
  !! L = f_leach * (SMINN_pool * 0.1_dp)    -> Ref:Xu Ri, 2008 
  !! unit conversion -> 1m water = 1000 kg water & 1mm water = 1 kg water 
  !!->CLM model
  !! L_online = (SMINN_pool * sf) * runoff_acc/moisture  !!.. multiplied with 0.001 m->kg H2O--(63 code)
  !! L_offline = SMINN * sf * R/M                        !!-> maxmoist in m --runoff from Echam output 160 code
  !! L = SMINN * sf * R/(alpha*max_moisture)
  
  ! --- Nitros oxide emissions(N2O)--N2O_Emissions() ----
  !! ====N2O Emissions from Natural Ecosystems Soils & Agricultural Ecosystems Soils
  !! N2O_emissions = 1% of N Deposition per day -->> IPCC methodology; Skiba et al, 1998; Bloemerts, M., 2009.
  !! & from BNF, the emission factor is 1.25% in IPCC methodology
  !! BNF converts N2 to NH4+ which undergone nitrification & denitrification process

  elemental pure subroutine N_process (NPP_act_yDayMean, nfix_to_sminn, ndep_to_sminn,                &
                                       ndep_forc, nfert_forc, nfert_to_sminn, sminn_to_denit, sminn_leach, &
                                       SMINN_pool, NetEcosyst_N_flux, runoff, alpha,   &
                                       max_moisture, NPP_run_mean, N2O_emissions_depfix, is_vegetation, &
                                       N2O_emissions_nfert) 

    real(dp),intent(in)    :: max_moisture              !! Depth of soil bucket in [m]
    real(dp),intent(out)   :: nfix_to_sminn             !! SMINN increase from biological nitrogen fixation [mol(N)/m^2(canopy)s]
                                                        !!    based on NPP
    real(dp),intent(out)   :: ndep_to_sminn             !! SMINN increase from N-deposition [mol(N)/m^2(canopy) s]
    real(dp),intent(in)    :: ndep_forc                 !! Deposition of N [mg(N)/m^2(grid box) year] based on atm forcing.
    real(dp),intent(in)    :: nfert_forc                !! Fertilizer of N [g(N)/m^2(grid box) year] based on external forcing.
    real(dp),intent(out)   :: nfert_to_sminn            !! SMINN increase from N-fertilizer [mol(N)/m^2(canopy) s] only in crops   
    real(dp),intent(out)   :: sminn_to_denit            !! SMINN loss by dentrification [mol(N)/m^2(canopy)s]
    real(dp),intent(out)   :: sminn_leach               !! SMINN loss by leaching [mol(N)/m^2(canopy)s] 
    real(dp),intent(in)    :: NPP_act_yDayMean          !! Mean NPP at the current day (this routine is called once a day) 
                                                        !!    [mol(C)/m^2 s]
    real(dp),intent(inout) :: SMINN_pool                !! mineral N in soils [mol(N)/m^2(canopy)]
    real(dp),intent(in)    :: runoff                    !! Yesterdays mean runoff [m/s] 
    real(dp),intent(in)    :: alpha                     !! Yesterdays mean water stress factor (between 0 and 1)
    real(dp),intent(out)   :: NetEcosyst_N_flux         !! Total balance of N gains and losses (positive for ecosystem gain) 
                                                        !!    [mol(N)/m^2(canopy)s]
    real(dp),intent(inout) :: NPP_run_mean              !! exponential running mean of NPP used for nitrogen fixation estimates
    real(dp),intent(out)   :: N2O_emissions_depfix      !! Nitros oxide emission [mol(N)/m^2(canopy)s]
    logical, intent(in)    :: is_vegetation             !!
    real(dp),intent(out)   :: N2O_emissions_nfert       !! Nitros oxide emission [mol(N)/m^2(canopy)s] from N fertilizer use    

    !! parameters

    REAL(dp),parameter  :: delta_time = 1.0_dp            !! time step (=calling interval) of this routine [days]
    REAL(dp),parameter  :: decay_time = 365.0_dp          !! Decay time of memory loss for exponential running mean of NPP [days]
    REAL(dp),parameter  :: sf  = 0.0005_dp                !! fraction of soil mineral N in soluble form < 0.05%
    REAL(dp),parameter  :: dnp = 0.0002_dp                !! denitrification rate per day (< .02%)
    REAL(dp),parameter  :: N2O_ef_Ndepo = 0.01_dp         !! N2O emissions Factor (.01); 1% in IPCC methodology
    REAL(dp),parameter  :: N2O_ef_Nfix  = 0.0125_dp       !! N2O emissions Factor (.0125); 1.25% in IPCC methodology
    REAL(dp),parameter  :: N2O_ef_Nfert = 0.0125_dp       !! N2O emissions Factor (.025); 2.5% in Davidson E.A (2009)

    !! locals

    REAL(dp) :: F_NPP_act               !! Exponential weighting factor for computing NPP_run_mean [1]
    REAL(dp) :: N_NPP_act               !! Normalization for computing NPP_run_mean [1] 
    REAL(dp) :: f_leach                 !! Fraction of leach
    REAL(dp) :: Nitrous_oxide_emis_Ndepo !! Nitros oxide emission by deposition

    nfix_to_sminn           = 0.0_dp
    ndep_to_sminn           = 0.0_dp
    sminn_to_denit          = 0.0_dp
    sminn_leach             = 0.0_dp
    NetEcosyst_N_flux       = 0.0_dp
    N2O_emissions_depfix    = 0.0_dp
    nfert_to_sminn          = 0.0_dp
    N2O_emissions_nfert     = 0.0_dp    

    !! -- Computation of an exponentially weighted running mean of NPP
    !!
    !!                        NPP(t+1)         delta_T
    !!    NPP_runmean(t+1) =  -------- + exp(- ------- ) * NPP_runmean(t)
    !!                          Norm             tau
    !!
    !! where the normalization Norm is given by
    !!
    !!                                delta_T
    !!    Norm = SUM(k=0,infty) exp(- ------- k) = (1 - exp(-delta_T/tau))^(-1)
    !!                                  tau

    F_NPP_act = exp(-delta_time/decay_time)
    N_NPP_act = 1._dp / (1._dp - F_NPP_act) 
    NPP_run_mean = NPP_act_yDayMean/N_NPP_act + F_NPP_act * NPP_run_mean

    !! --- Biological nitrogen fixation  based on NPP 

    nfix_to_sminn = &                               !! 0.7 is calibrated to match global n2fix (present day):e.g  see Cleveland 1999
         (0.7_dp * (1.0_dp - exp(-0.003_dp * NPP_run_mean * sec_per_year)))/sec_per_year   &
        *(g_N_perMol/g_C_perMol)                    !! Nitrogen fixation in [mol(N)/m^2 s]
    nfix_to_sminn = max(0.0_dp,nfix_to_sminn)       !! To prevent n-fixation for negative NPP 

    !! -- N Deposition (wet + Dry) ---------- 
    !! Distribute N depostion into all vegetation types; note: N-depostion in [mg(N)/m^2(grid box) year] 

    !! Distribute N depostion into all vegetation types; note: N-depostion in [mg(N)/m^2(grid box) year]
    !!                 (1e-3_dp/14 = conversion fact. mg(N) to mol(N); no cover_fract: homogeneous flux at all points in grid)
    IF (is_vegetation) THEN
      ndep_to_sminn = (ndep_forc/sec_per_year) * 1e-3_dp/g_N_perMol
    END IF

    !! -- Denitrification--------
                      
    sminn_to_denit = SMINN_pool * dnp/sec_per_day * alpha  !!dk: dnp depends on T and H2O; future plan: scale with N2O

    !! -- Nitros oxide emissions(N2O)-------
    Nitrous_oxide_emis_Ndepo = ndep_to_sminn * N2O_ef_Ndepo           !! N2O emissions due to N deposition
    N2O_emissions_depfix     = Nitrous_oxide_emis_Ndepo

    !! -- Leaching ----------------

    f_leach = min(1.0_dp, runoff*sec_per_day/max(1e-13_dp, alpha*max_moisture)) !! fraction of water removed during one day 
                                                                                !! (leaching limited to total bucket content)
    sminn_leach = f_leach * SMINN_pool * sf/sec_per_day !! Mol(N) leached per second (only soluble fraction  can be leached)
    !!-- Distribute N fertiliser into all crops types; note: N-fertilizer read in [g (N)/m^2(grid box) year]
    nfert_to_sminn = (nfert_forc * 0.5_dp/sec_per_year)/g_N_perMol     !!.. gm/14 = conversion fact. g(N) to mol(N)
                                                                       !!..& half lost as harvested N
    N2O_emissions_nfert = (nfert_forc /sec_per_year)/g_N_perMol * N2O_ef_Nfert   !! 2.5% of total Nfert use in crops  


    !--------------update SMINN_pool----

    SMINN_pool = SMINN_pool + nfix_to_sminn * sec_per_day           !! N-Fixation        
    SMINN_pool = SMINN_pool + ndep_to_sminn * sec_per_day           !! N-Deposition
    SMINN_pool = SMINN_pool + nfert_to_sminn * sec_per_day          !! N-Fertilizer application to crops      !dk new   
    SMINN_pool = SMINN_pool - sminn_to_denit* sec_per_day           !! N-Denitrification
    SMINN_pool = SMINN_pool - sminn_leach   * sec_per_day           !! N-Leaching
    SMINN_pool = SMINN_pool - N2O_emissions_depfix * sec_per_day    !! update SMINN_pool
    SMINN_pool = SMINN_pool - N2O_emissions_nfert * sec_per_day     !! update SMINN_pool

 !--------------Assessment of N inputs & lossess--- NetEcosyst_N_flux
 
    NetEcosyst_N_flux = nfix_to_sminn + ndep_to_sminn + nfert_to_sminn - (sminn_to_denit + sminn_leach) &
                       - N2O_emissions_depfix  - N2O_emissions_nfert 

  end subroutine N_process


  ! --- relocate_carbonAndNitrogen() -----------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of shifts in the cover fractions of the tiles for the carbon (C) and nitrogen (N) pools.
  ! This change in vegetation composition is either due to human land-cover change or due to vegetation dynamics.
  ! Concerning the C and N stored on those partitions of shrinking tiles that are lost to extending tiles during the time interval
  ! considered: In the case of landcover change the C and N from living plant pools (Cpool_green, Cpool_reserve and Cpool_woods plus
  ! associated N-pools) is relocated, partly by a release to the atmosphere (e.g. by fire stubbing) and the remaining C and N into
  ! soil pools.
  ! ATTENTION: In the case of vegetation dynamics it is assumed that the C and N from living plant pools is already lost to other
  ! reservoirs (atmosphere, litter, soil) by those processes that removed the vegetation (general carbon flow reflecting minimum
  ! mortality, fire, wind break) -- therefore in the case of vegetation dynamics this routine has to be called in combination with
  ! other routines that handle those destructive processes in order to assure C and N conservation. 
  !
  ! ATTENTION: 
  ! 1. If this routine is called with N-pools, then not only carbon is relocated, but also nitrogen.
  ! 2. To handle C and N relocation for landcover change, the necessary arrays must be present in the call.
  !    In absence of N-pools ("carbon only mode") only carbon is reshuffled.  
  !

  subroutine relocate_CarbonAndNitrogen(lctlib, surface, cf_old, cf_new, veg_fract_correction,      &
                                        ! Carbon
                                        Cpool_green, Cpool_woods, Cpool_reserve,                    &
                                        ! soil carbon pools of the old cbalance scheme
                                        Cpool_litter_green_ag, Cpool_litter_green_bg,               &
                                        Cpool_litter_wood_ag, Cpool_litter_wood_bg,                 &
                                        Cpool_slow,                                                 &
                                        ! Yasso litter and soil carbon pools
                                        YCpool_acid_ag1,                                             &
                                        YCpool_water_ag1,                                            &
                                        YCpool_ethanol_ag1,                                          &
                                        YCpool_nonsoluble_ag1,                                       &
                                        YCpool_acid_bg1,                                             &
                                        YCpool_water_bg1,                                            &
                                        YCpool_ethanol_bg1,                                          &
                                        YCpool_nonsoluble_bg1,                                       &
                                        YCpool_humus_1,                                              &
                                        YCpool_acid_ag2,                                             &
                                        YCpool_water_ag2,                                            &
                                        YCpool_ethanol_ag2,                                          &
                                        YCpool_nonsoluble_ag2,                                       &
                                        YCpool_acid_bg2,                                             &
                                        YCpool_water_bg2,                                            &
                                        YCpool_ethanol_bg2,                                          &
                                        YCpool_nonsoluble_bg2,                                       &
                                        YCpool_humus_2,                                              &
                                        LeafLit_coef,                                               &
                                        WoodLit_coef,                                               &
                                        ! LCC scheme 1: standard JSBACH
                                        C_2_atmos,                                                  &
                                        C_2_litterGreenPools,                                       &
                                        C_2_litterWoodPool_ag,C_2_litterWoodPool_bg,                &
                                        ! Nitrogen
                                        Npool_green, Npool_woods, Npool_mobile,                     &
                                        Npool_litter_green_ag, Npool_litter_green_bg,               &
                                        Npool_litter_wood_ag, Npool_litter_wood_bg,                 &
                                        Npool_slow,                                                 &
                                        SMINN_pool,                                                 &
                                        ! nitrogen variables for LCC scheme 1
                                        Nitrogen_2_atmos, Nitrogen_2_litterGreenPools,              &
                                        Nitrogen_2_litterWoodPool_ag, Nitrogen_2_litterWoodPool_bg, &
                                        Nitrogen_2_SMINNpool,                                       &
                                        ! LCC scheme 2: fluxes after Houghton
                                        Cpool_onSite,Cpool_paper,Cpool_construction,                &
                                        C_onSite_2_atmos, C_paper_2_atmos, C_construction_2_atmos,  &
                                        C_2_onSite, C_2_paper, C_2_construction,                    &
                                        lcc_scheme)
    
    USE mo_exception,        ONLY: finish
    USE mo_cbal_parameters,  ONLY: frac_wood_2_atmos, frac_green_2_atmos, frac_reserve_2_atmos, frac_mobile_2_atmos
    USE mo_jsbach_lctlib,    ONLY: lctlib_type 
    USE mo_land_surface,     ONLY: land_surface_type

    type(lctlib_type),intent(in) :: lctlib               !! PFT-specific constants
    type(land_surface_type),intent(in) :: surface        !! Access to cover_type
    real(dp),intent(in)    :: cf_old(:,:)                !! Cover fractions before landcover change or vegetation dynamics [] 
    real(dp),intent(in)    :: cf_new(:,:)                !! Cover fraction after landcover change or vegetation dynamics []
    real(dp),intent(in)    :: veg_fract_correction(:,:)  !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts for
                                                         !!    sparseness of vegetation)
    real(dp),intent(inout) :: Cpool_green(:,:)           !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)           !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)         !! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    
    ! Cbalance litter & soil pools:
    real(dp), optional, intent(inout) :: Cpool_litter_green_ag(:,:) !! above ground green litter carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: Cpool_litter_green_bg(:,:) !! below ground green litter carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: Cpool_litter_wood_ag(:,:)  !! above ground wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: Cpool_litter_wood_bg(:,:)  !! below ground wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: Cpool_slow(:,:)            !! slow soil carbon pool [mol(C)/m^2(canopy)]
    
    ! Yasso litter & soil pools:
    ! Size class 1; green litter
    real(dp), optional, intent(inout) :: YCpool_acid_ag1(:,:)       !! above ground acid soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_water_ag1(:,:)      !! above ground water soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_ethanol_ag1(:,:)    !! above ground ethanol soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_nonsoluble_ag1(:,:) !! above ground non carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_acid_bg1(:,:)       !! below ground acid soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_water_bg1(:,:)      !! below ground water soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_ethanol_bg1(:,:)    !! below ground ethanol soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_nonsoluble_bg1(:,:) !! below ground non soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_humus_1(:,:)         !! humus carbon pool [mol(C)/m^(canopy)2]
    ! Size class 2; woody litter
    real(dp), optional, intent(inout) :: YCpool_acid_ag2(:,:)       !! above ground acid soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_water_ag2(:,:)      !! above ground water soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_ethanol_ag2(:,:)    !! above ground ethanol soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_nonsoluble_ag2(:,:) !! above ground non carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_acid_bg2(:,:)       !! below ground acid soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_water_bg2(:,:)      !! below ground water soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_ethanol_bg2(:,:)    !! below ground ethanol soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_nonsoluble_bg2(:,:) !! below ground non soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_humus_2(:,:)         !! humus carbon pool [mol(C)/m^(canopy)2]

    ! Variables needed with standard jsbach LCC scheme
    real(dp),optional,intent(out) :: C_2_atmos(:)                 !! Amount of carbon directly emitted to atmosphere in this
                                                                  !!    timestep [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: C_2_litterGreenPools(:)      !! Amount of carbon relocated by land cover change from green
                                                                  !!    and reserve pool to below and above ground green litter
    real(dp),optional,intent(out) :: C_2_litterWoodPool_ag(:)     !! Amount of carbon relocated by land cover change from wood
                                                                  !!    pool to above ground woody litter pool
                                                                  !!    [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: C_2_litterWoodPool_bg(:)     !! Amount of carbon relocated by land cover change from wood
                                                                  !!    pool to below ground woody litter pool
                                                                  !!    [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(inout) ::  Cpool_onSite(:)           !! Amount of carbon remains in onSite anthro annual pool 
                                                                  !!  [mol(C)/m^2(vegetated area)] 
    real(dp),optional,intent(inout) ::  Cpool_paper(:)            !! Amount of carbon remains in anthro decadal pool
    real(dp),optional,intent(inout) ::  Cpool_construction(:)     !! Amount of carbon remains in anthro centinnial pool
    real(dp),optional,intent(out)   ::  C_onSite_2_atmos(:)       !! Carbon flux from anthro annual pool to atmos 
                                                                  !!  [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_paper_2_atmos(:)              !! Carbon flux from green anthro annual pool to atmos 
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_construction_2_atmos(:)       !! Carbon flux from woody anthro annual pool to atmos
                                                                        !!   [mol(C)/m^2(vegetated area)/day]
    real(dp),optional,intent(out)   ::  C_2_onSite(:)                   !! Carbon flux to onsite [mol(C)/m^2(vegetated area)/day]
    real(dp),optional,intent(out)   ::  C_2_construction(:)             !! Carbon flux from woody anthro annual pool to paper
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_2_paper(:)                    !! Carbon flux from woody anthro annual pool to
                                                                        !!   construction [mol(C)/m^2(vegetated area)/day]

    real(dp),intent(inout),optional :: Npool_green(:,:)           !! Green nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_woods(:,:)           !! Wood nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_mobile(:,:)          !! Plant mobile nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_green_ag(:,:) !! Above ground green litter nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional :: Npool_litter_green_bg(:,:) !! Below ground green litter nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional :: Npool_litter_wood_ag(:,:)  !! Wood litter nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional :: Npool_litter_wood_bg(:,:)  !! Wood litter nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional :: Npool_slow(:,:)            !! Slow soil nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: SMINN_pool(:,:)            !! Soil mineral nitrogen pool [mol(N)/m^2(canopy)]

    real(dp),optional,intent(out) :: Nitrogen_2_atmos(:)             !! Amount of nitrogen directly emitted to atmosphere in this 
                                                                     !! .. timestep [mol(N)/m^2(vegetated area)]
                                                                     !! .. to green litter pool [mol(N)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: Nitrogen_2_litterGreenPools(:)  !! Amount of nitrogen relocated by land cover change from green
                                                                     !!    and reserve pool to above and below ground green litter
                                                                     !!    pools [mol(N)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: Nitrogen_2_litterWoodPool_ag(:) !! Amount of nitrogen relocated by land cover change from wood
                                                                     !!    pool to above ground woody litter pool 
                                                                     !!    [mol(N)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: Nitrogen_2_litterWoodPool_bg(:) !! Amount of nitrogen relocated by land cover change from wood
                                                                     !!    pool to below ground woody litter pool
                                                                     !!    [mol(N)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: Nitrogen_2_SMINNpool(:)         !! Amount of nitrogen relocated by land cover change from wood
                                                                     !!    pool to soil mineral N pool (this is the surplus N that 
                                                                     !!    the woody litter N-pool cannot take up)   

    REAL(dp),dimension(5),intent(in),optional   :: LeafLit_coef      !! fractions to seperate fresh litter to non woody yasso pools
    REAL(dp),dimension(5),intent(in),optional   :: WoodLit_coef      !! fractions to seperate fresh litter to woody yasso pools
    
    integer,optional, intent(in) :: lcc_scheme
    
    !! locals
    !!   yasso variables 
    real(dp) :: Chumus_1_2_humus_1(1:size(cf_old,DIM=1))                 !! Total carbon relocated from humus pools of shrinking 
                                                                     !! ..  to humus pools of extending tiles
    real(dp) :: Cacid_ag1_2_acid_ag1(1:size(cf_old,DIM=1))             !! Total carbon relocated from above ground acid pools of
                                                                     !! ..  shrinking to above ground acid pools of extending tiles
    real(dp) :: Cacid_bg1_2_acid_bg1(1:size(cf_old,DIM=1))             !! Total carbon relocated from below ground acid pools of
                                                                     !! ..  shrinking to below ground acid pools of extending tiles
    real(dp) :: Cwater_ag1_2_water_ag1(1:size(cf_old,DIM=1))           !! Total carbon relocated from above ground water pools of
                                                                     !! ..  shrinking to above ground water pools of extending tiles
    real(dp) :: Cwater_bg1_2_water_bg1(1:size(cf_old,DIM=1))           !! Total carbon relocated from below ground water pools of
                                                                     !! ..  shrinking to below ground water pools of extending tiles
    real(dp) :: Cethanol_ag1_2_ethanol_ag1(1:size(cf_old,DIM=1))       !! Total carbon relocated from above ground ethanol pools of
                                                                     !! ..  shrinking to above ground ethanol pools of ext. tiles
    real(dp) :: Cethanol_bg1_2_ethanol_bg1(1:size(cf_old,DIM=1))       !! Total carbon relocated from below ground ethanol pools of
                                                                     !! ..  shrinking to below ground ethanol pools of ext. tiles
    real(dp) :: Cnonsoluble_ag1_2_nonsoluble_ag1(1:size(cf_old,DIM=1)) !! Total carbon relocated from above ground nonsoluble pools
                                                                       !! of shrinking to above ground nonsoluble pools of ext. tiles
    real(dp) :: Cnonsoluble_bg1_2_nonsoluble_bg1(1:size(cf_old,DIM=1)) !! Total carbon relocated from below ground nonsoluble pools
                                                                       !! of shrinking to below ground nonsoluble pools of ext. tiles
    real(dp) :: Chumus_2_2_humus_2(1:size(cf_old,DIM=1))                 !! Total carbon relocated from humus pools of shrinking 
                                                                     !! ..  to humus pools of extending tiles
    real(dp) :: Cacid_ag2_2_acid_ag2(1:size(cf_old,DIM=1))             !! Total carbon relocated from above ground acid pools of
                                                                     !! ..  shrinking to above ground acid pools of extending tiles
    real(dp) :: Cacid_bg2_2_acid_bg2(1:size(cf_old,DIM=1))             !! Total carbon relocated from below ground acid pools of
                                                                     !! ..  shrinking to below ground acid pools of extending tiles
    real(dp) :: Cwater_ag2_2_water_ag2(1:size(cf_old,DIM=1))           !! Total carbon relocated from above ground water pools of
                                                                     !! ..  shrinking to above ground water pools of extending tiles
    real(dp) :: Cwater_bg2_2_water_bg2(1:size(cf_old,DIM=1))           !! Total carbon relocated from below ground water pools of
                                                                     !! ..  shrinking to below ground water pools of extending tiles
    real(dp) :: Cethanol_ag2_2_ethanol_ag2(1:size(cf_old,DIM=1))       !! Total carbon relocated from above ground ethanol pools of
                                                                     !! ..  shrinking to above ground ethanol pools of ext. tiles
    real(dp) :: Cethanol_bg2_2_ethanol_bg2(1:size(cf_old,DIM=1))       !! Total carbon relocated from below ground ethanol pools of
                                                                     !! ..  shrinking to below ground ethanol pools of ext. tiles
    real(dp) :: Cnonsoluble_ag2_2_nonsoluble_ag2(1:size(cf_old,DIM=1)) !! Total carbon relocated from above ground nonsoluble pools
                                                                       !! of shrinking to above ground nonsoluble pools of ext. tiles
    real(dp) :: Cnonsoluble_bg2_2_nonsoluble_bg2(1:size(cf_old,DIM=1)) !! Total carbon relocated from below ground nonsoluble pools
                                                                       !! of shrinking to below ground nonsoluble pools of ext. tiles





    real(dp) :: C_2_litterGreenPool_ag(1:size(cf_old,DIM=1))    !! Carbon relocated by land cover change from green and reserve
                                                                !!    pool to above ground green litter
    real(dp) :: C_2_litterGreenPool_bg(1:size(cf_old,DIM=1))    !! Carbon relocated by land cover change from green and reserve
                                                                !!    pool to below ground green litter
    real(dp) :: Nitrogen_2_litterGreenPool_ag(1:size(cf_old,DIM=1))  !! Nitrogen relocated by land cover change from green and 
                                                                     !!    reserve pool to above ground green litter
    real(dp) :: Nitrogen_2_litterGreenPool_bg(1:size(cf_old,DIM=1))  !! Nitrogen relocated by land cover change from green and
                                                                     !!    reserve pool to below ground green litter
    real(dp) :: Cslow_2_slow(1:size(cf_old,DIM=1))                   !! Total carbon relocated from slow pools of shrinking 
                                                                     !! ..  to slow pools of extending tiles
    real(dp) :: Nslow_2_slow(1:size(cf_old,DIM=1))                   !! Total nitrogen relocated from slow pools of shrinking 
                                                                     !! ..  to slow pools of extending tiles
    real(dp) :: ClitGreen_ag_2_litGreen_ag(1:size(cf_old,DIM=1))  !! Total carbon relocated from above ground green litter pools of 
                                                                  !! shrinking to above ground green litter pool of extending tiles
    real(dp) :: ClitGreen_bg_2_litGreen_bg(1:size(cf_old,DIM=1))  !! Total carbon relocated from below ground green litter pools of 
                                                                  !! shrinking to below ground green litter pool of extending tiles
    
    real(dp) :: NlitGreen_ag_2_litGreen_ag(1:size(cf_old,DIM=1))  !! Total N relocated from above ground green litter pools of
                                                                  !! shrinking to above ground green litter pools of extending tiles
    real(dp) :: NlitGreen_bg_2_litGreen_bg(1:size(cf_old,DIM=1))  !! Total N relocated from below ground green litter pools of
                                                                  !! shrinking to below ground green litter pools of extending tiles
    real(dp) :: ClitterWood_ag_2_litterWood_ag(1:size(cf_old,DIM=1))    !! Total C relocated from ag wood litter pools of shrinking
                                                                        !! ..  to ag wood litter pools of extending tiles
    real(dp) :: ClitterWood_bg_2_litterWood_bg(1:size(cf_old,DIM=1))    !! Total C relocated from bg wood litter pools of shrinking
                                                                        !! ..  to bg wood litter pools of extending tiles
    real(dp) :: NlitterWood_ag_2_litterWood_ag(1:size(cf_old,DIM=1))    !! Total N relocated from ag wood litter pools of shrinking
                                                                        !! ..  to ag wood litter pools of extending tiles
    real(dp) :: NlitterWood_bg_2_litterWood_bg(1:size(cf_old,DIM=1))    !! Total N relocated from bg wood litter pools of shrinking
                                                                        !! ..  to bg wood litter pools of extending tiles
    real(dp) :: SMINN_2_SMINN(1:size(cf_old,DIM=1))                     !! Total nitrogen relocated from soil mineral N pools of
                                                                        !! ..  shrinking to SMINN pools of extending tile

    logical  :: tiles_extend(1:size(cf_old,DIM=1),1:size(cf_old,DIM=2)) !! Logical mask indicating those tiles that are extending
                                                                        !! by landcover change
    real(dp) :: cf_rescale(1:size(cf_old,DIM=1),1:size(cf_old,DIM=2))   !! Scale factor depending on changing landcover
    real(dp) :: shrinkRatio (1:size(cf_old,DIM=1),1:size(cf_old,DIM=2)) !! Ratio old and new cover fractions of extending tiles
    real(dp) :: sum_cf_new_old(1:size(cf_old,DIM=1))                    !! Sum of cover fractions added to old cover fractions of
                                                                        !! extending tiles
    logical  :: land_cover_change                                       !! true:  C relocation due to land cover change
                                                                        !! false: C relocation due to (natural) vegetation dynamics 
    integer  :: ntiles, i, nidx
    logical  :: nitrogenMode !! Flag: whether nitrogen shall be handled: only if N-pools are present in the call
    logical  :: with_yasso   !! Flag: whether yasso litter and soil pools shall be handled

    land_cover_change = .false.
    NitrogenMode      = .false.
    with_yasso        = .false.

    !! Avoid wrong use of the routine!
    !! This routine handles either land cover change or vegetation dynamics:
    !! For land cover change pass frac_wood_2_atmos, frac_green_2_atmos, frac_reserve_2_atmos, carbon_2_atmos,
    !! carbon_2_litter_green_bgSoilPool and carbon_2_slowSoilPool to this routine.

    !! Check for landcover change in carbon mode
    IF (PRESENT(C_2_atmos)             .OR. PRESENT(C_2_litterGreenPools) .OR. &
        PRESENT(C_2_litterWoodPool_ag) .OR. PRESENT(C_2_litterWoodPool_bg)) THEN
        IF (.NOT. (PRESENT(C_2_atmos) .AND. PRESENT(C_2_litterGreenPools) .AND. &
             PRESENT(C_2_litterWoodPool_ag) .AND. PRESENT(C_2_litterWoodPool_bg)) &
             ) CALL finish('relocate_CarbonAndNitrogen()','at least one variable missing to handle land cover change')
        land_cover_change=.TRUE.
    END IF
 
    !! Check for nitrogen mode
    if(present(Npool_green)           .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(Npool_slow)            .or. &
       present(SMINN_pool) ) then
       nitrogenMode=.true.
       if (.not. (present(Npool_green)           .and. present(Npool_woods)           .and. &
                  present(Npool_mobile)          .and. present(Npool_litter_green_ag) .and. &
                  present(Npool_litter_wood_ag)  .and. present(Npool_litter_wood_bg)  .and. &
                  present(Npool_litter_green_bg) .and. present(Npool_slow)            .and. &
                  present(SMINN_pool)) &
          ) call finish('relocate_CarbonAndNitrogen()','at least one variable missing to handle nitrogen')
    end if

    !!check for .not. yasso mode
    if (present(Cpool_litter_green_ag) .or. present(Cpool_litter_green_bg) .or. &
        present(Cpool_litter_wood_ag)  .or. present(Cpool_litter_wood_bg)  .or. &
        present(Cpool_slow) ) then
       if (.not.(present(Cpool_litter_green_ag) .and. present(Cpool_litter_green_bg) .and. &
          present(Cpool_litter_wood_ag)  .and. present(Cpool_litter_wood_bg)         .and. &
          present(Cpool_slow) ) ) &
          call finish('relocate_CarbonAndNitrogen()','at least one variable missing to handle cbalance')
     end if

    !! Check if yasso pools shall be handled
    if(present(YCpool_acid_ag1)        .or. present(YCpool_acid_bg1)            .or. &
       present(YCpool_water_ag1)       .or. present(YCpool_water_bg1)           .or. &
       present(YCpool_ethanol_ag1)     .or. present(YCpool_ethanol_bg1)         .or. &
       present(YCpool_nonsoluble_ag1)  .or. present(YCpool_nonsoluble_bg1)      .or. &
       present(YCpool_humus_1)         .or.                                          &
       present(YCpool_acid_ag2)        .or. present(YCpool_acid_bg2)            .or. &
       present(YCpool_water_ag2)       .or. present(YCpool_water_bg2)           .or. &
       present(YCpool_ethanol_ag2)     .or. present(YCpool_ethanol_bg2)         .or. &
       present(YCpool_nonsoluble_ag2)  .or. present(YCpool_nonsoluble_bg2)      .or. &
       present(YCpool_humus_2) ) then
       with_yasso=.true.
       if (.not. (present(YCpool_acid_ag1)        .and. present(YCpool_acid_bg1)        .and. &
                  present(YCpool_water_ag1)       .and. present(YCpool_water_bg1)       .and. &
                  present(YCpool_ethanol_ag1)     .and. present(YCpool_ethanol_bg1)     .and. &
                  present(YCpool_nonsoluble_ag1)  .and. present(YCpool_nonsoluble_bg1)  .and. &
                  present(YCpool_humus_1)                                               .and. &
                  present(YCpool_acid_ag2)        .and. present(YCpool_acid_bg2)        .and. &
                  present(YCpool_water_ag2)       .and. present(YCpool_water_bg2)       .and. &
                  present(YCpool_ethanol_ag2)     .and. present(YCpool_ethanol_bg2)     .and. &
                  present(YCpool_nonsoluble_ag2)  .and. present(YCpool_nonsoluble_bg2)  .and. &
                  present(YCpool_humus_2)) &
          ) call finish('relocate_CarbonAndNitrogen()','at least one variable missing to handle yasso pools')
        if(present(Npool_green)) call finish('relocate_CarbonAndNitrogen()','yasso in combination with nitrogen: no!')
    end if
   

    !! preparations

    ntiles = size(cf_old,DIM=2)
    nidx = size(cf_old,DIM=1)

    !! determine extending and shrinking tiles

    WHERE (cf_new(:,:) > cf_old(:,:))  
       tiles_extend(:,:) = .TRUE.
       cf_rescale(:,:) = 0._dp
    ELSEWHERE 
       tiles_extend(:,:) = .FALSE.
       cf_rescale(:,:) = (cf_old(:,:) - cf_new(:,:)) * veg_fract_correction(:,:)
    END WHERE

    !! determine amount of carbon released from living plant pools to atmosphere, litter (or anthropogenic pools) and soil

    IF (land_cover_change) THEN
       
       IF (lcc_scheme==1) THEN ! Standard JSBACH scheme

          !! C to be relocate from carbon living pools of shrinking tiles to litter pools of expanding tiles
          !! ... notice cf_rescale is now 0 for expanding tiles -> no need for usage of mask
          C_2_litterWoodPool_ag(:) = SUM(cf_rescale(:,:) *  (1._dp - frac_wood_2_atmos )           &
                                         *        frac_wood_aboveGround  * Cpool_woods(:,:),DIM=2)
          C_2_litterWoodPool_bg(:) = SUM(cf_rescale(:,:) *  (1._dp - frac_wood_2_atmos )           &
                                         * (1._dp-frac_wood_aboveGround) * Cpool_woods(:,:),DIM=2)
          C_2_litterGreenPools(:)  = SUM(cf_rescale(:,:) * (  (1._dp - frac_green_2_atmos  ) * Cpool_green  (:,:)   &
                                                            + (1._dp - frac_reserve_2_atmos) * Cpool_reserve(:,:))  &
                                        ,DIM=2)
          C_2_litterGreenPool_ag(:) = frac_green_aboveGround * C_2_litterGreenPools(:)        
          C_2_litterGreenPool_bg(:) = C_2_litterGreenPools(:) - C_2_litterGreenPool_ag(:)
         
          !! Living C from shrinking tiles to atmosphere
          C_2_atmos(:) = SUM((  frac_green_2_atmos   * Cpool_green  (:,:) &
                              + frac_wood_2_atmos    * Cpool_woods  (:,:) &
                              + frac_reserve_2_atmos * Cpool_reserve(:,:) &
                             ) * cf_rescale(:,:),DIM=2)
          
          IF (nitrogenMode) THEN 
             
             !! N to be relocated from carbon living pools of shrinking tiles to litter pools of expanding tiles
             !! so that the C/N-ratio is appropriate for the respective litter pool
             Nitrogen_2_litterWoodPool_ag(:) = SUM(cf_rescale(:,:) * (1.0_dp - frac_wood_2_atmos )                           &
                                                   *        frac_wood_aboveGround  *  Cpool_woods(:,:)/cn_litter_wood,DIM=2)
             Nitrogen_2_litterWoodPool_bg(:) = SUM(cf_rescale(:,:) * (1.0_dp - frac_wood_2_atmos )                           &
                                                   * (1._dp-frac_wood_aboveGround) *  Cpool_woods(:,:)/cn_litter_wood,DIM=2)
             !! ... Note: reserve pool contains no N!
             Nitrogen_2_litterGreenPools(:)  = SUM(cf_rescale(:,:) * (1.0_dp - frac_green_2_atmos) * Npool_green(:,:),DIM=2)
             
             Nitrogen_2_litterGreenPool_ag(:) = frac_green_aboveGround * Nitrogen_2_litterGreenPools(:)
             Nitrogen_2_litterGreenPool_bg(:) = Nitrogen_2_litterGreenPools(:) -Nitrogen_2_litterGreenPool_ag(:) 

             !! When stubbing wood to woody litter the (assume: cn_litter_wood >= cn_woods) surplus N not allocated in the
             !! woody litter pool and the remaining N from the plant mobile N pool are put into the soil mineral N pool
             Nitrogen_2_SMINNpool(:) = SUM((    (1._dp - frac_wood_2_atmos) * Cpool_woods(:,:)   & 
                                              * (1._dp/cn_woods - 1._dp/cn_litter_wood )         &
                                            + (1.0_dp - frac_mobile_2_atmos) * Npool_mobile(:,:) &
                                           ) * cf_rescale(:,:),DIM=2)
            
             !! Living N from shrinking tiles to atmosphere
             Nitrogen_2_atmos(:) = SUM(  (  frac_green_2_atmos  * Npool_green(:,:)    &
                                          + frac_wood_2_atmos   * Npool_woods(:,:)    &
                                          + frac_mobile_2_atmos * Npool_mobile(:,:) ) &
                                       * cf_rescale(:,:),DIM=2)
             
          END IF
         
       ELSE 
       !! Anthropogenic C-reloc (lcc_scheme==2). details: T Kato et al. 2009 based on grand slam protocol by Houghton et al. (1983)

          ! Update anthro pools
          IF (.NOT. with_yasso) THEN
             CALL C_loss_and_update_anthro_pools(nidx,ntiles,lcc_scheme,lctlib,surface,       &
                                                 Cpool_green,Cpool_woods,Cpool_reserve,       &
                 Cpool_litter_green_ag        =  Cpool_litter_green_ag,                       &
                 Cpool_litter_wood_ag         =  Cpool_litter_wood_ag,                        &
                 Cpool_onSite                 =  Cpool_onSite     ,                           &
                 Cpool_paper                  =  Cpool_paper,                                 &
                 Cpool_construction           =  Cpool_construction,                          &
                 scale_fac                    =  cf_rescale,                                  &
                 C_2_onSite                   =  C_2_onSite,                                  &
                 C_2_paper                    =  C_2_paper,                                   &
                 C_2_construction             =  C_2_construction,                            &
                 C_onSite_2_atmos             =  C_onSite_2_atmos,                            &
                 C_paper_2_atmos              =  C_paper_2_atmos,                             &
                 C_construction_2_atmos       =  C_construction_2_atmos,                      &
                 C_2_atmos                    =  C_2_atmos)

             C_2_litterWoodPool_bg (:) = SUM(cf_rescale(:,:) * Cpool_woods(:,:), DIM=2) * (1._dp - frac_wood_aboveGround) 

             C_2_litterGreenPool_bg(:) = (1._dp - frac_green_aboveGround) * &
                                         SUM((Cpool_green(:,:) + Cpool_reserve(:,:)) * cf_rescale(:,:), DIM=2)
          ELSE
             CALL C_loss_and_update_anthro_pools(nidx,ntiles,lcc_scheme,lctlib,surface,       &
                                                 Cpool_green,Cpool_woods,Cpool_reserve,       &
                 YCpool_acid_ag1              =  YCpool_acid_ag1,                             &
                 YCpool_water_ag1             =  YCpool_water_ag1,                            &
                 YCpool_ethanol_ag1           =  YCpool_ethanol_ag1,                          &
                 YCpool_nonsoluble_ag1        =  YCpool_nonsoluble_ag1,                       &
                 YCpool_acid_ag2              =  YCpool_acid_ag2,                             &
                 YCpool_water_ag2             =  YCpool_water_ag2,                            &
                 YCpool_ethanol_ag2           =  YCpool_ethanol_ag2,                          &
                 YCpool_nonsoluble_ag2        =  YCpool_nonsoluble_ag2,                       &
                 Cpool_onSite                 =  Cpool_onSite,                                &
                 Cpool_paper                  =  Cpool_paper,                                 &
                 Cpool_construction           =  Cpool_construction,                          &
                 scale_fac                    =  cf_rescale,                                  &
                 C_2_onSite                   =  C_2_onSite,                                  &
                 C_2_paper                    =  C_2_paper,                                   &
                 C_2_construction             =  C_2_construction,                            &
                 C_onSite_2_atmos             =  C_onSite_2_atmos,                            &
                 C_paper_2_atmos              =  C_paper_2_atmos,                             &
                 C_construction_2_atmos       =  C_construction_2_atmos,                      &
                 C_2_atmos                    =  C_2_atmos)
          END IF
          C_2_litterWoodPool_bg (:) = SUM(cf_rescale(:,:) * Cpool_woods(:,:), DIM=2) * (1._dp - frac_wood_aboveGround)
          C_2_litterGreenPool_bg(:) = (1._dp - frac_green_aboveGround) * &
                                      SUM((Cpool_green(:,:) + Cpool_reserve(:,:)) * cf_rescale(:,:), DIM=2)
       ENDIF
    
    ELSE
!! C-relocation from cover fraction change from the dynamic vegetation
!! CHR: IN THIS PART CONSIDERATION OF NITROGEN IS STILL MISSING!!!

       !! rescale living plant pools of shrinking tiles 
       
       WHERE(.NOT. tiles_extend(:,:) .AND. cf_new(:,:) >= fract_small)
          shrinkRatio(:,:) = cf_old(:,:) / cf_new(:,:)
       ELSEWHERE
          shrinkRatio(:,:) = 1._dp
       END WHERE
       
       Cpool_green  (:,:) = Cpool_green  (:,:) * shrinkRatio(:,:)
       Cpool_reserve(:,:) = Cpool_reserve(:,:) * shrinkRatio(:,:)
       Cpool_woods  (:,:) = Cpool_woods  (:,:) * shrinkRatio(:,:)

       if (nitrogenMode) then
         Npool_green  (:,:) = Npool_green  (:,:) * shrinkRatio(:,:)
         Npool_mobile (:,:) = Npool_mobile (:,:) * shrinkRatio(:,:)
         Npool_woods  (:,:) = Npool_woods  (:,:) * shrinkRatio(:,:)
       end if
       
    END IF ! Land cover change/dynveg 

    !! C and N to be relocated from slow pool of shrinking to extending tiles (all schemes and reasons for cover frac changes)
    IF (.NOT. with_yasso) THEN
       Cslow_2_slow(:) = SUM(cf_rescale(:,:) * Cpool_slow(:,:),DIM=2)
    ELSE
       Chumus_1_2_humus_1(:) = SUM(cf_rescale(:,:) * YCpool_humus_1(:,:),DIM=2)
       Chumus_2_2_humus_2(:) = SUM(cf_rescale(:,:) * YCpool_humus_2(:,:),DIM=2)
    END IF

    IF (nitrogenMode) THEN
      Nslow_2_slow(:) = SUM(cf_rescale(:,:) * Npool_slow(:,:),DIM=2)
    ENDIF

    IF (lcc_scheme==1 .OR. .NOT. land_cover_change) THEN !! dynveg or lcc standard JSBACH
       
       !! C to be relocated from litter pool of shrinking to extending tiles
       IF (.NOT. with_yasso) THEN
          ClitGreen_ag_2_litGreen_ag    (:) = SUM(cf_rescale(:,:) * Cpool_litter_green_ag(:,:),DIM=2)
          ClitterWood_ag_2_litterWood_ag(:) = SUM(cf_rescale(:,:) * Cpool_litter_wood_ag (:,:),DIM=2)
       ELSE
          Cacid_ag1_2_acid_ag1            (:) = SUM(cf_rescale(:,:) * YCpool_acid_ag1      (:,:),DIM=2)
          Cwater_ag1_2_water_ag1          (:) = SUM(cf_rescale(:,:) * YCpool_water_ag1     (:,:),DIM=2)
          Cethanol_ag1_2_ethanol_ag1      (:) = SUM(cf_rescale(:,:) * YCpool_ethanol_ag1   (:,:),DIM=2)
          Cnonsoluble_ag1_2_nonsoluble_ag1(:) = SUM(cf_rescale(:,:) * YCpool_nonsoluble_ag1(:,:),DIM=2)
          Cacid_bg1_2_acid_bg1            (:) = SUM(cf_rescale(:,:) * YCpool_acid_bg1      (:,:),DIM=2)
          Cwater_bg1_2_water_bg1          (:) = SUM(cf_rescale(:,:) * YCpool_water_bg1     (:,:),DIM=2)
          Cethanol_bg1_2_ethanol_bg1      (:) = SUM(cf_rescale(:,:) * YCpool_ethanol_bg1   (:,:),DIM=2)
          Cnonsoluble_bg1_2_nonsoluble_bg1(:) = SUM(cf_rescale(:,:) * YCpool_nonsoluble_bg1(:,:),DIM=2)   
          Cacid_ag2_2_acid_ag2            (:) = SUM(cf_rescale(:,:) * YCpool_acid_ag2      (:,:),DIM=2)
          Cwater_ag2_2_water_ag2          (:) = SUM(cf_rescale(:,:) * YCpool_water_ag2     (:,:),DIM=2)
          Cethanol_ag2_2_ethanol_ag2      (:) = SUM(cf_rescale(:,:) * YCpool_ethanol_ag2   (:,:),DIM=2)
          Cnonsoluble_ag2_2_nonsoluble_ag2(:) = SUM(cf_rescale(:,:) * YCpool_nonsoluble_ag2(:,:),DIM=2)
          Cacid_bg2_2_acid_bg2            (:) = SUM(cf_rescale(:,:) * YCpool_acid_bg2      (:,:),DIM=2)
          Cwater_bg2_2_water_bg2          (:) = SUM(cf_rescale(:,:) * YCpool_water_bg2     (:,:),DIM=2)
          Cethanol_bg2_2_ethanol_bg2      (:) = SUM(cf_rescale(:,:) * YCpool_ethanol_bg2   (:,:),DIM=2)
          Cnonsoluble_bg2_2_nonsoluble_bg2(:) = SUM(cf_rescale(:,:) * YCpool_nonsoluble_bg2(:,:),DIM=2)   
       END IF

       IF (nitrogenMode) THEN

          !! N to be relocated from litter pool of shrinking to extending tiles
          NlitGreen_ag_2_litGreen_ag    (:) = SUM(cf_rescale(:,:) * Npool_litter_green_ag(:,:),DIM=2)
          NlitGreen_bg_2_litGreen_bg    (:) = SUM(cf_rescale(:,:) * Npool_litter_green_bg(:,:),DIM=2)
          NlitterWood_ag_2_litterWood_ag(:) = SUM(cf_rescale(:,:) * Npool_litter_wood_ag (:,:),DIM=2)  
          NlitterWood_bg_2_litterWood_bg(:) = SUM(cf_rescale(:,:) * Npool_litter_wood_bg (:,:),DIM=2)
          !! N to be relocated from soil mineral N pool of shrinking to extending tiles 
          SMINN_2_SMINN                 (:) = SUM(cf_rescale(:,:) * SMINN_pool           (:,:),DIM=2)
          
       END IF
       
    END IF

    IF (.NOT. with_yasso) THEN
       ClitGreen_bg_2_litGreen_bg    (:) = SUM(cf_rescale(:,:) * Cpool_litter_green_bg(:,:),DIM=2)
       ClitterWood_bg_2_litterWood_bg(:) = SUM(cf_rescale(:,:) * Cpool_litter_wood_bg (:,:),DIM=2)
    ENDIF

    !! Lower down the C and N density [mol/m^2(canopy)] for extending tiles (all schemes and reasons for cover frac changes)
    WHERE(tiles_extend(:,:) .AND. cf_new(:,:) >= fract_small)
       shrinkRatio(:,:) = cf_old(:,:) / cf_new(:,:)
    ELSEWHERE
       shrinkRatio(:,:) = 1._dp
    END WHERE

    Cpool_green(:,:)           = Cpool_green(:,:)           * shrinkRatio(:,:) 
    Cpool_reserve(:,:)         = Cpool_reserve(:,:)         * shrinkRatio(:,:)
    Cpool_woods(:,:)           = Cpool_woods(:,:)           * shrinkRatio(:,:)
    
    IF (.NOT. with_yasso) THEN
      Cpool_litter_green_ag(:,:) = Cpool_litter_green_ag(:,:) * shrinkRatio(:,:)
      Cpool_litter_green_bg(:,:) = Cpool_litter_green_bg(:,:) * shrinkRatio(:,:)
      Cpool_litter_wood_ag(:,:)  = Cpool_litter_wood_ag(:,:)  * shrinkRatio(:,:)
      Cpool_litter_wood_bg(:,:)  = Cpool_litter_wood_bg(:,:)  * shrinkRatio(:,:)
      Cpool_slow(:,:)            = Cpool_slow(:,:)            * shrinkRatio(:,:)
    ELSE  
      YCpool_acid_ag1      (:,:) = YCpool_acid_ag1      (:,:)   * shrinkRatio(:,:) 
      YCpool_water_ag1     (:,:) = YCpool_water_ag1     (:,:)   * shrinkRatio(:,:) 
      YCpool_ethanol_ag1   (:,:) = YCpool_ethanol_ag1   (:,:)   * shrinkRatio(:,:) 
      YCpool_nonsoluble_ag1(:,:) = YCpool_nonsoluble_ag1(:,:)   * shrinkRatio(:,:) 
      YCpool_acid_bg1      (:,:) = YCpool_acid_bg1      (:,:)   * shrinkRatio(:,:) 
      YCpool_water_bg1     (:,:) = YCpool_water_bg1     (:,:)   * shrinkRatio(:,:) 
      YCpool_ethanol_bg1   (:,:) = YCpool_ethanol_bg1   (:,:)   * shrinkRatio(:,:) 
      YCpool_nonsoluble_bg1(:,:) = YCpool_nonsoluble_bg1(:,:)   * shrinkRatio(:,:) 
      YCpool_humus_1       (:,:) = YCpool_humus_1       (:,:)   * shrinkRatio(:,:) 
      YCpool_acid_ag2      (:,:) = YCpool_acid_ag2      (:,:)   * shrinkRatio(:,:) 
      YCpool_water_ag2     (:,:) = YCpool_water_ag2     (:,:)   * shrinkRatio(:,:) 
      YCpool_ethanol_ag2   (:,:) = YCpool_ethanol_ag2   (:,:)   * shrinkRatio(:,:) 
      YCpool_nonsoluble_ag2(:,:) = YCpool_nonsoluble_ag2(:,:)   * shrinkRatio(:,:) 
      YCpool_acid_bg2      (:,:) = YCpool_acid_bg2      (:,:)   * shrinkRatio(:,:) 
      YCpool_water_bg2     (:,:) = YCpool_water_bg2     (:,:)   * shrinkRatio(:,:) 
      YCpool_ethanol_bg2   (:,:) = YCpool_ethanol_bg2   (:,:)   * shrinkRatio(:,:) 
      YCpool_nonsoluble_bg2(:,:) = YCpool_nonsoluble_bg2(:,:)   * shrinkRatio(:,:) 
      YCpool_humus_2       (:,:) = YCpool_humus_2       (:,:)   * shrinkRatio(:,:) 
    END IF
     
    IF (nitrogenMode) THEN 
       Npool_green(:,:)           = Npool_green(:,:)           * shrinkRatio(:,:) 
       Npool_woods(:,:)           = Npool_woods(:,:)           * shrinkRatio(:,:)
       Npool_mobile(:,:)          = Npool_mobile(:,:)          * shrinkRatio(:,:)
       Npool_litter_green_ag(:,:) = Npool_litter_green_ag(:,:) * shrinkRatio(:,:)
       Npool_litter_green_bg(:,:) = Npool_litter_green_bg(:,:) * shrinkRatio(:,:)
       Npool_litter_wood_ag(:,:)  = Npool_litter_wood_ag(:,:)  * shrinkRatio(:,:)
       Npool_litter_wood_bg(:,:)  = Npool_litter_wood_bg(:,:)  * shrinkRatio(:,:)
       Npool_slow(:,:)            = Npool_slow(:,:)            * shrinkRatio(:,:)
       SMINN_pool(:,:)            = SMINN_pool(:,:)            * shrinkRatio(:,:)
    END IF

    !! sum of cover fractions added to cover fractions of extending tiles
    sum_cf_new_old(:) = 0._dp
    DO i = 1,ntiles
       WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))
          sum_cf_new_old(:) = sum_cf_new_old(:) + (cf_new(:,i) - cf_old(:,i)) * veg_fract_correction(:,i)
       END WHERE
    END DO

    !! New temporary usage of cf_rescale for fewer computations and thus faster execution
    DO i = 1,ntiles
       WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))
          cf_rescale(:,i) = (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
       ELSEWHERE
          cf_rescale(:,i) = 0._dp
       END WHERE
    END DO

    !! Distribute freed carbon and nitrogen to ground pools of extending tiles
    DO i = 1,ntiles

       !! Soil pools are treated equally, no matter the source for the cf-change
       IF (.NOT. with_yasso) THEN
          Cpool_slow(:,i) = Cpool_slow(:,i) + Cslow_2_slow(:) * cf_rescale(:,i)
       ELSE                   
          YCpool_humus_1(:,i) = YCpool_humus_1(:,i) + Chumus_1_2_humus_1(:) * cf_rescale(:,i)
          YCpool_humus_2(:,i) = YCpool_humus_2(:,i) + Chumus_2_2_humus_2(:) * cf_rescale(:,i)
       END IF
       IF (nitrogenMode) THEN
          Npool_slow(:,i) = Npool_slow(:,i) + Nslow_2_slow(:) * cf_rescale(:,i)
          SMINN_pool(:,i) = SMINN_pool(:,i) + SMINN_2_SMINN(:) * cf_rescale(:,i)
       ENDIF

       IF (land_cover_change) THEN

          IF (nitrogenMode) THEN
             SMINN_pool(:,i) = SMINN_pool(:,i) + Nitrogen_2_SMINNpool(:) * cf_rescale(:,i)
          ENDIF

          IF (lcc_scheme==1) THEN
             IF (.NOT. with_yasso) THEN
                WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))
                   Cpool_litter_green_ag(:,i) =  &
                   Cpool_litter_green_ag(:,i) + (C_2_litterGreenPool_ag(:) + ClitGreen_ag_2_litGreen_ag    (:)) * cf_rescale(:,i)
                   Cpool_litter_wood_ag (:,i) =  &
                   Cpool_litter_wood_ag (:,i) + (C_2_litterWoodPool_ag (:) + ClitterWood_ag_2_litterWood_ag(:)) * cf_rescale(:,i)
                END WHERE
             ELSE
                WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))

                   YCpool_acid_ag1(:,i) =                                                       & 
                        YCpool_acid_ag1(:,i) + (C_2_litterGreenPool_ag(:) *LeafLit_coef(1)      & ! Input leaf litter
                        + Cacid_ag1_2_acid_ag1(:)) *cf_rescale(:,i)                                ! Input from other tile

                   YCpool_acid_ag2(:,i) =                                                       & 
                        YCpool_acid_ag2(:,i) + (C_2_litterWoodPool_ag(:) * WoodLit_coef(1)      & ! Input woody litter
                        + Cacid_ag2_2_acid_ag2(:)) *cf_rescale(:,i)                                ! Input from other tile

                   YCpool_acid_bg1(:,i) =                                                       & 
                        YCpool_acid_bg1(:,i) + (C_2_litterGreenPool_bg(:) *LeafLit_coef(1)      & ! Input leaf litter
                        + Cacid_bg1_2_acid_bg1(:)) *cf_rescale(:,i)                                ! Input from other tile

                   YCpool_acid_bg2(:,i) =                                                       & 
                        YCpool_acid_bg2(:,i) + (C_2_litterWoodPool_bg(:) * WoodLit_coef(1)      & ! Input woody litter
                        + Cacid_bg2_2_acid_bg2(:)) *cf_rescale(:,i)                                ! Input from other tile


                   YCpool_water_ag1(:,i) =                                                      & 
                        YCpool_water_ag1(:,i)  + (C_2_litterGreenPool_ag(:) *LeafLit_coef(2)     &
                        + Cwater_ag1_2_water_ag1(:)) *cf_rescale(:,i)  

                   YCpool_water_ag2(:,i) =                                                      & 
                        YCpool_water_ag2(:,i)  + ( C_2_litterWoodPool_ag(:) * WoodLit_coef(2)   & 
                        + Cwater_ag2_2_water_ag2(:)) *cf_rescale(:,i)  

                   YCpool_water_bg1(:,i) =                                                      & 
                        YCpool_water_bg1(:,i)  + (C_2_litterGreenPool_bg(:) *LeafLit_coef(2)     &
                        + Cwater_bg1_2_water_bg1(:)) *cf_rescale(:,i)  

                   YCpool_water_bg2(:,i) =                                                      & 
                        YCpool_water_bg2(:,i)  + ( C_2_litterWoodPool_bg(:) * WoodLit_coef(2)   & 
                        + Cwater_bg2_2_water_bg2(:)) *cf_rescale(:,i)  


                   YCpool_ethanol_ag1(:,i) =                                                      & 
                        YCpool_ethanol_ag1(:,i)  + (C_2_litterGreenPool_ag(:) *LeafLit_coef(3)     &
                        + Cethanol_ag1_2_ethanol_ag1(:)) *cf_rescale(:,i)  

                   YCpool_ethanol_ag2(:,i) =                                                      & 
                        YCpool_ethanol_ag2(:,i)  + ( C_2_litterWoodPool_ag(:) * WoodLit_coef(3)   & 
                        + Cethanol_ag2_2_ethanol_ag2(:)) *cf_rescale(:,i)  

                   YCpool_ethanol_bg1(:,i) =                                                      & 
                        YCpool_ethanol_bg1(:,i)  + (C_2_litterGreenPool_bg(:) *LeafLit_coef(3)     &
                        + Cethanol_bg1_2_ethanol_bg1(:)) *cf_rescale(:,i)  

                   YCpool_ethanol_bg2(:,i) =                                                      & 
                        YCpool_ethanol_bg2(:,i)  + ( C_2_litterWoodPool_bg(:) * WoodLit_coef(3)   & 
                        + Cethanol_bg2_2_ethanol_bg2(:)) *cf_rescale(:,i)  

                   YCpool_nonsoluble_ag1(:,i) =                                                      & 
                        YCpool_nonsoluble_ag1(:,i)  + (C_2_litterGreenPool_ag(:) *LeafLit_coef(4)    &
                        + Cnonsoluble_ag1_2_nonsoluble_ag1(:)) *cf_rescale(:,i)  

                   YCpool_nonsoluble_ag2(:,i) =                                                      & 
                        YCpool_nonsoluble_ag2(:,i)  + ( C_2_litterWoodPool_ag(:) * WoodLit_coef(4)   & 
                        + Cnonsoluble_ag2_2_nonsoluble_ag2(:)) *cf_rescale(:,i)  

                   YCpool_nonsoluble_bg1(:,i) =                                                      & 
                        YCpool_nonsoluble_bg1(:,i)  + (C_2_litterGreenPool_bg(:) *LeafLit_coef(4)    &
                        + Cnonsoluble_bg1_2_nonsoluble_bg1(:)) *cf_rescale(:,i)  

                   YCpool_nonsoluble_bg2(:,i) =                                                      & 
                        YCpool_nonsoluble_bg2(:,i)  + ( C_2_litterWoodPool_bg(:) * WoodLit_coef(4)   & 
                        + Cnonsoluble_bg2_2_nonsoluble_bg2(:)) *cf_rescale(:,i)  

                   YCpool_humus_1(:,i) =                                                              & 
                        YCpool_humus_1(:,i) &
                                            + ((C_2_litterGreenPool_ag(:) + C_2_litterGreenPool_bg(:) &
                                                      )*LeafLit_coef(5)) *cf_rescale(:,i)  
                   YCpool_humus_2(:,i) =                                                              & 
                        YCpool_humus_2(:,i) & 
                                            + ((C_2_litterWoodPool_ag(:)  + C_2_litterWoodPool_bg(:)  &
                                                      )*WoodLit_coef(5)) *cf_rescale(:,i)  

                END WHERE
             END IF

             IF (nitrogenMode) THEN 
                WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))
                   Npool_litter_green_ag(:,i)= &
                   Npool_litter_green_ag(:,i)+(Nitrogen_2_litterGreenPool_ag(:)+NlitGreen_ag_2_litGreen_ag    (:))*cf_rescale(:,i)
                   Npool_litter_green_bg(:,i)= &
                   Npool_litter_green_bg(:,i)+(Nitrogen_2_litterGreenPool_bg(:)+NlitGreen_bg_2_litGreen_bg    (:))*cf_rescale(:,i)
                   Npool_litter_wood_ag (:,i)= &
                   Npool_litter_wood_ag (:,i)+(Nitrogen_2_litterWoodPool_ag (:)+NlitterWood_ag_2_litterWood_ag(:))*cf_rescale(:,i)
                   Npool_litter_wood_bg (:,i)= &
                   Npool_litter_wood_bg (:,i)+(Nitrogen_2_litterWoodPool_bg (:)+NlitterWood_bg_2_litterWood_bg(:))*cf_rescale(:,i)
                END WHERE
             END IF
          ELSE
             !! When lcc and anthropogenic pools are used, the freed carbon goes into the anthropogenic pools instead (see above)
          END IF
          IF (.NOT. with_yasso) THEN
             WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp)) ! StW: Will C_2_litterXPool_bg always be 0 for shrinking tiles?
                Cpool_litter_green_bg(:,i) =  &
                Cpool_litter_green_bg(:,i) + (C_2_litterGreenPool_bg(:) + ClitGreen_bg_2_litGreen_bg    (:)) * cf_rescale(:,i)
                Cpool_litter_wood_bg (:,i) =  &
                Cpool_litter_wood_bg (:,i) + (C_2_litterWoodPool_bg (:) + ClitterWood_bg_2_litterWood_bg(:)) * cf_rescale(:,i)
             END WHERE
          ENDIF
       ELSE
          !! cover frac changes from dynveg
          IF (.NOT. with_yasso) THEN
             WHERE((cf_new(:,i) - cf_old(:,i)) > EPSILON(1._dp) .AND. sum_cf_new_old(:) > EPSILON(1._dp))
                Cpool_litter_green_ag(:,i) = Cpool_litter_green_ag(:,i) + ClitGreen_ag_2_litGreen_ag    (:) * cf_rescale(:,i)
                Cpool_litter_green_bg(:,i) = Cpool_litter_green_bg(:,i) + ClitGreen_bg_2_litGreen_bg    (:) * cf_rescale(:,i)
                Cpool_litter_wood_ag (:,i) = Cpool_litter_wood_ag (:,i) + ClitterWood_ag_2_litterWood_ag(:) * cf_rescale(:,i)
                Cpool_litter_wood_bg (:,i) = Cpool_litter_wood_bg (:,i) + ClitterWood_bg_2_litterWood_bg(:) * cf_rescale(:,i)
             END WHERE
          ELSE
             WHERE((cf_new(:,i) - cf_old(:,i)) > EPSILON(1._dp) .AND. sum_cf_new_old(:) > EPSILON(1._dp))
                YCpool_acid_ag1(:,i)       = YCpool_acid_ag1(:,i)       + Cacid_ag1_2_acid_ag1                (:) * cf_rescale(:,i)
                YCpool_acid_bg1(:,i)       = YCpool_acid_bg1(:,i)       + Cacid_bg1_2_acid_bg1                (:) * cf_rescale(:,i)
                YCpool_water_ag1(:,i)      = YCpool_water_ag1(:,i)      + Cwater_ag1_2_water_ag1              (:) * cf_rescale(:,i)
                YCpool_water_bg1(:,i)      = YCpool_water_bg1(:,i)      + Cwater_bg1_2_water_bg1              (:) * cf_rescale(:,i)
                YCpool_ethanol_ag1(:,i)    = YCpool_ethanol_ag1(:,i)    + Cethanol_ag1_2_ethanol_ag1          (:) * cf_rescale(:,i)
                YCpool_ethanol_bg1(:,i)    = YCpool_ethanol_bg1(:,i)    + Cethanol_bg1_2_ethanol_bg1          (:) * cf_rescale(:,i)
                YCpool_nonsoluble_ag1(:,i) = YCpool_nonsoluble_ag1(:,i) + Cnonsoluble_ag1_2_nonsoluble_ag1    (:) * cf_rescale(:,i)
                YCpool_nonsoluble_bg1(:,i) = YCpool_nonsoluble_bg1(:,i) + Cnonsoluble_bg1_2_nonsoluble_bg1    (:) * cf_rescale(:,i)
                YCpool_acid_ag2(:,i)       = YCpool_acid_ag2(:,i)       + Cacid_ag2_2_acid_ag2                (:) * cf_rescale(:,i)
                YCpool_acid_bg2(:,i)       = YCpool_acid_bg2(:,i)       + Cacid_bg2_2_acid_bg2                (:) * cf_rescale(:,i)
                YCpool_water_ag2(:,i)      = YCpool_water_ag2(:,i)      + Cwater_ag2_2_water_ag2              (:) * cf_rescale(:,i)
                YCpool_water_bg2(:,i)      = YCpool_water_bg2(:,i)      + Cwater_bg2_2_water_bg2              (:) * cf_rescale(:,i)
                YCpool_ethanol_ag2(:,i)    = YCpool_ethanol_ag2(:,i)    + Cethanol_ag2_2_ethanol_ag2          (:) * cf_rescale(:,i)
                YCpool_ethanol_bg2(:,i)    = YCpool_ethanol_bg2(:,i)    + Cethanol_bg2_2_ethanol_bg2          (:) * cf_rescale(:,i)
                YCpool_nonsoluble_ag2(:,i) = YCpool_nonsoluble_ag2(:,i) + Cnonsoluble_ag2_2_nonsoluble_ag2    (:) * cf_rescale(:,i)
                YCpool_nonsoluble_bg2(:,i) = YCpool_nonsoluble_bg2(:,i) + Cnonsoluble_bg2_2_nonsoluble_bg2    (:) * cf_rescale(:,i)
             END WHERE
          END IF
          IF (nitrogenMode) THEN 
             WHERE((cf_new(:,i) - cf_old(:,i)) > EPSILON(1._dp) .AND. sum_cf_new_old(:) > EPSILON(1._dp))
                Npool_litter_green_ag(:,i) = Npool_litter_green_ag(:,i) + NlitGreen_ag_2_litGreen_ag    (:) * cf_rescale(:,i)
                Npool_litter_green_bg(:,i) = Npool_litter_green_bg(:,i) + NlitGreen_bg_2_litGreen_bg    (:) * cf_rescale(:,i)
                Npool_litter_wood_ag (:,i) = Npool_litter_wood_ag (:,i) + NlitterWood_ag_2_litterWood_ag(:) * cf_rescale(:,i)
                Npool_litter_wood_bg (:,i) = Npool_litter_wood_bg (:,i) + NlitterWood_bg_2_litterWood_bg(:) * cf_rescale(:,i)
             END WHERE
          END IF

       END IF
    END DO

  END SUBROUTINE relocate_CarbonAndNitrogen

  ! --- relocate_carbon_desert() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of desert dynamics for the carbon pools. More precisely:
  ! It is assumed that no carbon fluxes have to be considered as these are already taken into account by the other routines in this
  ! module. So, the amount of carbon per unit area is only rescaled to conserve the carbon mass.
  ! The only exception is the green pool, which should never exceed its prescribed maximum carbon content.
  ! In this (very exceptional) case the excessing carbon is transferred to the green litter pool.

  subroutine relocate_carbon_desert(nidx, ntiles, is_present, is_glacier, veg_ratio_max, veg_ratio_max_old, &
                                         lai, sla, Cpool_green, Cpool_woods, Cpool_reserve,             &
                                         Cpool_litter_green_ag, Cpool_litter_green_bg,                  &
                                         Cpool_litter_wood_ag,Cpool_litter_wood_bg, Cpool_slow,         &
                                         !Nitrogen pools
                                         Npool_green, Npool_woods, Npool_mobile, Npool_litter_green_ag, &
                                         Npool_litter_green_bg, Npool_litter_wood_ag, Npool_litter_wood_bg, &
                                         Npool_slow, SMINN_pool,                                        &
                                         ! Yasso litter and soil carbon pools
                                         YCpool_acid_ag1,                                                &
                                         YCpool_water_ag1,                                               &
                                         YCpool_ethanol_ag1,                                             &
                                         YCpool_nonsoluble_ag1,                                          &
                                         YCpool_acid_bg1,                                                &
                                         YCpool_water_bg1,                                               &
                                         YCpool_ethanol_bg1,                                             &
                                         YCpool_nonsoluble_bg1,                                          &
                                         YCpool_humus_1,                                                 &
                                         YCpool_acid_ag2,                                                &
                                         YCpool_water_ag2,                                               &
                                         YCpool_ethanol_ag2,                                             &
                                         YCpool_nonsoluble_ag2,                                          &
                                         YCpool_acid_bg2,                                                &
                                         YCpool_water_bg2,                                               &
                                         YCpool_ethanol_bg2,                                             &
                                         YCpool_nonsoluble_bg2,                                          &
                                         YCpool_humus_2,                                                 &
                                         LeafLitcoef)
                                       
    USE mo_exception,        ONLY: finish

    integer,intent(in)     :: nidx                        !! Vector length
    integer,intent(in)     :: ntiles                      !! Number of tiles
    logical,intent(in)     :: is_present(:,:)             !! Logical mask for grid points treated by jsbach
    logical,intent(in)     :: is_glacier(:,:)             !! Glacier mask
    real(dp),intent(in)    :: veg_ratio_max(:)            !! Fractional cover of vegetated area
    real(dp),intent(in)    :: veg_ratio_max_old(:)        !! Fractional cover of vegetated area of last year
    real(dp),intent(in)    :: lai(:,:)                    !! leaf area index [m^2(leaf)/m^2(canopy)]
    real(dp),intent(in)    :: sla(:,:)                    !! specific leaf area [m^2(leaf)/mol(C)]
    real(dp),intent(inout) :: Cpool_green(:,:)            !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)            !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)          !! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_green_ag(:,:)    !! Above ground green litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_green_bg(:,:)    !! Below ground green litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_wood_ag(:,:)     !! Above ground woody litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_wood_bg(:,:)     !! Below ground woody litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_slow(:,:)               !! Value of slow soil carbon pool [mol(C)/m^2(canopy)]
    !   YASSO 
    ! SIZE CLASS 1 : green litter
    !   above ground C pools
    real(dp),intent(inout),optional  :: YCpool_acid_ag1(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_water_ag1(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_ethanol_ag1(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_nonsoluble_ag1(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout),optional  :: YCpool_acid_bg1(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_water_bg1(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_ethanol_bg1(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_nonsoluble_bg1(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_humus_1(:,:)            !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 : woody litter
    !   above ground C pools
    real(dp),intent(inout),optional  :: YCpool_acid_ag2(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_water_ag2(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_ethanol_ag2(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_nonsoluble_ag2(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout),optional  :: YCpool_acid_bg2(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_water_bg2(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_ethanol_bg2(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_nonsoluble_bg2(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_humus_2(:,:)            !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]


    !   parameters 
    real(dp),intent(in)   ,optional  :: LeafLitcoef(:,:,:)           !! Factor to spread non woody litterfall into yasso pools [ ]
    real(dp),intent(inout),optional :: Npool_green(:,:)              !! Green nitrogen pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_woods(:,:)              !! Wood nitrogen pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_mobile(:,:)             !! Mobile nitrogen pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_green_ag(:,:)    !! Above ground green litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_green_bg(:,:)    !! Below ground green litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_wood_ag(:,:)     !! Above ground woody litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_wood_bg(:,:)     !! Below ground woody litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_slow(:,:)               !! Slow soil nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: SMINN_pool(:,:)               !! Mineral soil N [mol(N)/m^2(canopy)]

    !! locals

    real(dp)  ::  veg_ratio_old_new(nidx)          ! ratio of old fraction of vegetated land to new fraction of vegetated land
    real(dp)  ::  cpool_green_excess(nidx,ntiles),npool_green_excess(nidx,ntiles)
    logical   ::  nitrogenMode
    logical  :: with_yasso
    nitrogenMode = .false.
    with_yasso  = .FALSE.    

    if(present(Npool_green)           .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(Npool_slow)            .or. &
       present(SMINN_pool) ) then
       nitrogenMode=.true.
    end if

    !! Check if yasso or cbalance litter and soil pools should be handled

    IF (PRESENT(YCpool_acid_ag1)        .OR. PRESENT(YCpool_acid_bg1)            .OR. &
        PRESENT(YCpool_water_ag1)       .OR. PRESENT(YCpool_water_bg1)           .OR. &
        PRESENT(YCpool_ethanol_ag1)     .OR. PRESENT(YCpool_ethanol_bg1)         .OR. &
        PRESENT(YCpool_nonsoluble_ag1)  .OR. PRESENT(YCpool_nonsoluble_bg1)      .OR. &
        PRESENT(YCpool_humus_1)         .OR.                                          &
        PRESENT(YCpool_acid_ag2)        .OR. PRESENT(YCpool_acid_bg2)            .OR. &
        PRESENT(YCpool_water_ag2)       .OR. PRESENT(YCpool_water_bg2)           .OR. &
        PRESENT(YCpool_ethanol_ag2)     .OR. PRESENT(YCpool_ethanol_bg2)         .OR. &
        PRESENT(YCpool_nonsoluble_ag2)  .OR. PRESENT(YCpool_nonsoluble_bg2)      .OR. &
        PRESENT(YCpool_humus_2)         .OR. PRESENT(LeafLitcoef)  ) THEN
        with_yasso=.TRUE.
        IF (.NOT. (PRESENT(YCpool_acid_ag1)        .AND. PRESENT(YCpool_acid_bg1)        .AND. &
                   PRESENT(YCpool_water_ag1)       .AND. PRESENT(YCpool_water_bg1)       .AND. &
                   PRESENT(YCpool_ethanol_ag1)     .AND. PRESENT(YCpool_ethanol_bg1)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag1)  .AND. PRESENT(YCpool_nonsoluble_bg1)  .AND. &
                   PRESENT(YCpool_humus_1)         .AND.                                       &
                   PRESENT(YCpool_acid_ag2)        .AND. PRESENT(YCpool_acid_bg2)        .AND. &
                   PRESENT(YCpool_water_ag2)       .AND. PRESENT(YCpool_water_bg2)       .AND. &
                   PRESENT(YCpool_ethanol_ag2)     .AND. PRESENT(YCpool_ethanol_bg2)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag2)  .AND. PRESENT(YCpool_nonsoluble_bg2)  .AND. &
                   PRESENT(YCpool_humus_2)         .AND. PRESENT(LeafLitcoef)  )             &
            ) THEN 
            CALL finish('relocate_carbon_desert()','at least one variable missing to handle yasso pools')
       END IF
    END IF

    WHERE (veg_ratio_max(:) > 0.5_dp*fract_small)
       veg_ratio_old_new(:) = veg_ratio_max_old(:) / veg_ratio_max(:)
    ELSEWHERE
       veg_ratio_old_new(:) = 1._dp                ! never happens on non-glacier land points
    ENDWHERE

    WHERE (.NOT. is_glacier(:,:) .AND. is_present(:,:))
       cpool_green_excess(:,:) = MAX(0._dp,cpool_green(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) &
                                           - greenC2leafC * lai(:,:) / sla(:,:))
       cpool_green(:,:) = MIN(cpool_green(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2), &
                              greenC2leafC * lai(:,:) / sla(:,:))
       cpool_reserve(:,:) = cpool_reserve(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
       cpool_woods(:,:) = cpool_woods(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
    END WHERE
    IF (.NOT. with_yasso) THEN
       WHERE (.NOT. is_glacier(:,:) .AND. is_present(:,:))
          cpool_litter_green_ag(:,:) = cpool_litter_green_ag(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) &
                                      + cpool_green_excess(:,:)
          cpool_litter_green_bg(:,:) = cpool_litter_green_bg(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          cpool_litter_wood_ag(:,:) = cpool_litter_wood_ag(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          cpool_litter_wood_bg(:,:) = cpool_litter_wood_bg(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          cpool_slow(:,:) = cpool_slow(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
       END WHERE
    ELSE
       WHERE (.NOT. is_glacier(:,:) .AND. is_present(:,:))
          YCpool_acid_ag1       (:,:) = YCpool_acid_ag1       (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) &
                                      + cpool_green_excess(:,:)   * LeafLitcoef(:,:,1)
          YCpool_water_ag1      (:,:) = YCpool_water_ag1      (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) &
                                      + cpool_green_excess(:,:)   * LeafLitcoef(:,:,2)
          YCpool_ethanol_ag1    (:,:) = YCpool_ethanol_ag1    (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) &
                                      + cpool_green_excess(:,:)   * LeafLitcoef(:,:,3)
          YCpool_nonsoluble_ag1 (:,:) = YCpool_nonsoluble_ag1 (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) &
                                      + cpool_green_excess(:,:)   * LeafLitcoef(:,:,4)
          YCpool_humus_1        (:,:) = YCpool_humus_1        (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) &
                                      + cpool_green_excess(:,:)   * LeafLitcoef(:,:,5)
          YCpool_acid_bg1       (:,:) = YCpool_acid_bg1       (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_water_bg1      (:,:) = YCpool_water_bg1      (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_ethanol_bg1    (:,:) = YCpool_ethanol_bg1    (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_nonsoluble_bg1 (:,:) = YCpool_nonsoluble_bg1 (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)

          YCpool_acid_ag2       (:,:) = YCpool_acid_ag2       (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_water_ag2      (:,:) = YCpool_water_ag2      (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_ethanol_ag2    (:,:) = YCpool_ethanol_ag2    (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_nonsoluble_ag2 (:,:) = YCpool_nonsoluble_ag2 (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_humus_2        (:,:) = YCpool_humus_2        (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_acid_bg2       (:,:) = YCpool_acid_bg2       (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_water_bg2      (:,:) = YCpool_water_bg2      (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_ethanol_bg2    (:,:) = YCpool_ethanol_bg2    (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_nonsoluble_bg2 (:,:) = YCpool_nonsoluble_bg2 (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
       END WHERE
    END IF

    
    if (nitrogenMode) then
      WHERE (.NOT. is_glacier(:,:) .AND. is_present(:,:))
      !dk needs to be changed for potential future variable C/N values in JSBACH
         npool_green_excess(:,:) = cpool_green_excess(:,:)  / cn_green    
         npool_green(:,:) = cpool_green(:,:) / cn_green
         npool_mobile(:,:) = npool_mobile(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
         npool_woods(:,:) = cpool_woods(:,:) / cn_woods
         npool_litter_green_ag(:,:) = npool_litter_green_ag(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) &
                                    + npool_green_excess(:,:)
         npool_litter_green_bg(:,:) = npool_litter_green_bg(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
         npool_litter_wood_ag(:,:) = cpool_litter_wood_ag(:,:) / cn_litter_wood 
         npool_litter_wood_bg(:,:) = cpool_litter_wood_bg(:,:) / cn_litter_wood
         npool_slow(:,:) = cpool_slow(:,:) / cn_slow 
         sminn_pool(:,:) = sminn_pool(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) 
      END WHERE
    end if           

  end subroutine relocate_carbon_desert


  ! --- relocate_carbon_fire() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of fire for the carbon pools. More precisely:
  ! It is assumed that for the burned area the carbon from the above ground litter pools (Cpool_litter_green_ag,
  ! Cpool_litter_wood_ag) is completely released to the atmosphere and the carbon from the living plant pools
  ! (Cpool_green, Cpool_reserve, Cpool_woods) is partly released to the atmosphere and partly relocated
  ! into the litter pools.
  !
  subroutine relocate_carbon_fire(nidx, ntiles, with_yasso, frac_burned, cf,                    &
                                       veg_fract_correction, frac_wood_2_atmos_fire,            &
                                       Cpool_green, Cpool_woods, Cpool_reserve,                 &
                                       carbon_2_GreenLitterPools,                               &
                                       carbon_2_WoodLitterPools, carbon_2_atmos,                &
                                       Cpool_litter_green_ag,Cpool_litter_green_bg,             &
                                       Cpool_litter_wood_ag, Cpool_litter_wood_bg,              &
                                       ! Yasso litter and soil carbon pools
                                       YCpool_acid_ag1,      YCpool_acid_ag2,                   &
                                       YCpool_water_ag1,     YCpool_water_ag2,                  &
                                       YCpool_ethanol_ag1,   YCpool_ethanol_ag2,                &
                                       YCpool_nonsoluble_ag1,YCpool_nonsoluble_ag2,             &
                                       YCpool_acid_bg1,      YCpool_acid_bg2,                   &
                                       YCpool_water_bg1,     YCpool_water_bg2,                  &
                                       YCpool_ethanol_bg1,   YCpool_ethanol_bg2,                &
                                       YCpool_nonsoluble_bg1,YCpool_nonsoluble_bg2,             &
                                       YCpool_humus_1,       YCpool_humus_2,                    &
                                       LeafLitcoef,          WoodLitcoef,                       &
                                       !Nitrogen Pools and fluxes
                                       Npool_green, Npool_woods, Npool_mobile, SMINN_pool,      &
                                       Npool_litter_green_ag,Npool_litter_green_bg,             &
                                       Npool_litter_wood_ag,Npool_litter_wood_bg,               &
                                       nitrogen_2_GreenLitterPools, nitrogen_2_WoodLitterPools, &
                                       nitrogen_2_atmos, nitrogen_2_sminn)
!!
!! CHR 09-11-03: DIESE ROUTINE MUSS AUS DYNVEG HERAUSGEZOGEN WERDEN UM AUCH AGRARFLAECHEN IN DAS FEUER MIT EINZUBEZIEHEN
!!  

    USE mo_exception,        ONLY: finish

    integer,intent(in)     :: nidx                        !! Vector length
    integer,intent(in)     :: ntiles                      !! Number of tiles
    logical,intent(in)     :: with_yasso
    real(dp),intent(in)    :: frac_burned(:,:)            !! Fraction of the vegetated area of each tile burned till the last call
                                                          !!    of this routine
    real(dp),intent(in)    :: cf(:,:)                     !! Cover fractions
    real(dp),intent(in)    :: veg_fract_correction(:,:)   !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts for
                                                          !!    sparseness of vegetation)
    real(dp),intent(in)    :: frac_wood_2_atmos_fire      !! Fraction of above ground wood immediately emitted to atmosphere by fire
    real(dp),intent(inout) :: Cpool_green(:,:)            !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)            !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)          !! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: carbon_2_GreenLitterPools(:)!! Amount of carbon relocated by wind damage
                                                          !! .. to the green litter pools [mol(C)/m^2(vegetated area)]
    real(dp),intent(inout) :: carbon_2_WoodLitterPools(:) !! Amount of carbon relocated by wind damage and fire
                                                          !! .. to the wood litter pools [mol(C)/m^2(vegetated area)]
    real(dp),intent(out)   :: carbon_2_atmos(:)           !! Amount of carbon immediately released by fire
    real(dp),intent(inout) :: Cpool_litter_green_ag(:,:)    !! Above ground green litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_green_bg(:,:)    !! Below ground green litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood_ag(:,:)     !! Wood litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood_bg(:,:)     !! Wood litter carbon pool [mol(C)/m^2(canopy)]

    !   YASSO 
    ! SIZE CLASS 1 ; GREEN
    !   above ground C pools
    real(dp),intent(inout) :: YCpool_acid_ag1(:,:)             !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_water_ag1(:,:)            !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_ethanol_ag1(:,:)          !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_nonsoluble_ag1(:,:)       !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout) :: YCpool_acid_bg1(:,:)             !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_water_bg1(:,:)            !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_ethanol_bg1(:,:)          !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_nonsoluble_bg1(:,:)       !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_humus_1(:,:)              !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 ; WOOD
    !   above ground C pools
    real(dp),intent(inout) :: YCpool_acid_ag2(:,:)             !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_water_ag2(:,:)            !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_ethanol_ag2(:,:)          !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_nonsoluble_ag2(:,:)       !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout) :: YCpool_acid_bg2(:,:)             !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_water_bg2(:,:)            !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_ethanol_bg2(:,:)          !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_nonsoluble_bg2(:,:)       !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_humus_2(:,:)              !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    !   parameters 
    real(dp),intent(in)    :: LeafLitcoef(:,:,:)              !! Factor to spread non woody litterfall into yasso pools [ ]
    real(dp),intent(in)    :: WoodLitcoef(:,:,:)              !! Factor to spread woody litterfall into yasso pools [ ]

    ! Nitrogen
    real(dp),intent(inout),optional :: Npool_green(:,:)              !! Green nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_woods(:,:)              !! Wood nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_mobile(:,:)             !! Mobile nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: SMINN_pool(:,:)               !! Mineral nitrogen pool [mol(N)/m^2(canopy)]    
    real(dp),intent(inout),optional :: Npool_litter_green_ag(:,:)    !! Above ground green litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_green_bg(:,:)    !! Below ground green litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_wood_ag(:,:)     !! Wood litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_wood_bg(:,:)     !! Wood litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: nitrogen_2_GreenLitterPools(:)!! Amount of nitrogen relocated by fire to the green litter
                                                                     !!    pools [mol(N)/m^2(vegetated area)]
    real(dp),intent(inout),optional :: nitrogen_2_WoodLitterPools(:) !! Amount of nitrogen relocated by fire and wind damage 
                                                                     !!    to the wood litter pools [mol(N)/m^2(vegetated area)]
    real(dp),intent(out),  optional :: nitrogen_2_atmos(:)           !! Amount of nitrogen immediately released by fire and wind 
                                                                     !!    damage to the atmosphere [mol(N)/m^2(vegetated area)] 
    real(dp),intent(inout),optional :: nitrogen_2_sminn(:)           !! Amount of nitrogen relocated by fire and wind damage to the
                                                                     !!    soil mineral N pool [mol(N)/m^2(vegetated area)]

    !! locals

    real(dp)  ::  carbon_2_GreenLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  carbon_2_WoodLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  carbon_2_atmos_tiled(nidx,ntiles)
    real(dp)  ::  N_2_GreenLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  N_2_WoodLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  N_2_atmos_tiled(nidx,ntiles)
    real(dp)  ::  N_2_sminn_tiled(nidx,ntiles)
    logical   ::  nitrogenMode

    nitrogenMode=.false.
    if(present(Npool_green)           .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(SMINN_pool) ) then
       nitrogenMode=.true.
       if (.not. (present(Npool_green)           .and. present(Npool_woods)           .and. &
                  present(Npool_mobile)          .and. present(Npool_litter_green_ag) .and. &
                  present(Npool_litter_wood_ag)  .and. present(Npool_litter_wood_bg)  .and. &
                  present(Npool_litter_green_bg) .and. present(SMINN_pool) ) ) &
                  CALL finish('relocate_carbon_fire()','at least one variable missing to handle nitrogen dynamics')
    end if    

    !! preparations
    carbon_2_GreenLitterPools_tiled(:,:) = 0._dp
    carbon_2_WoodLitterPools_tiled(:,:) = 0._dp
    N_2_GreenLitterPools_tiled(:,:) = 0._dp
    N_2_WoodLitterPools_tiled(:,:) = 0._dp
    N_2_atmos_tiled(:,:) = 0._dp   
    N_2_sminn_tiled(:,:)   = 0._dp

    !! diagnose amount of carbon released from the green and reserve pool to the green litter pools
    !! diagnose amount of carbon released from the wood pool to the woody litter pools
    !! determine amount of carbon released from the wood pool, living tissue pools, and litter pools to the atmosphere

    carbon_2_GreenLitterPools_tiled(:,:) = (Cpool_green(:,:) + Cpool_reserve(:,:)) * (1._dp - frac_green_aboveGround) * &
                                           cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
    carbon_2_WoodLitterPools_tiled(:,:) = Cpool_woods(:,:) * (1._dp - frac_wood_2_atmos_fire * frac_wood_aboveGround) * &
                                          cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)


    IF (.NOT. with_yasso) THEN
       carbon_2_atmos_tiled(:,:) = ((Cpool_green(:,:) + Cpool_reserve(:,:)) * frac_green_aboveGround + &
                                   Cpool_litter_green_ag(:,:) + Cpool_litter_wood_ag(:,:) + &
                                   Cpool_woods(:,:) * frac_wood_2_atmos_fire * frac_wood_aboveGround) * &
                                   cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
    ELSE
       carbon_2_atmos_tiled(:,:) = ((Cpool_green(:,:) + Cpool_reserve(:,:)) * frac_green_aboveGround      + &
                                   YCpool_water_ag1(:,:) + YCpool_acid_ag1(:,:) + YCPool_ethanol_ag1(:,:) + &
                                   YCpool_nonsoluble_ag1(:,:) +                                             &
                                   YCpool_water_ag2(:,:) + YCpool_acid_ag2(:,:) + YCPool_ethanol_ag2(:,:) + &
                                   YCpool_nonsoluble_ag2(:,:) +                                             &
                                   Cpool_woods(:,:) * frac_wood_2_atmos_fire * frac_wood_aboveGround) *     &
                                   cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
    END IF
    
    if (nitrogenMode) then
      !Npool_mobile above ground to atmosphere
      !Npool_mobile below ground to sminn_pool
      !Npool_litter_ag to atmosphere
      !SMINN pool not affected by fire, only through N released by other pools
        N_2_atmos_tiled(:,:) = ((Npool_green(:,:) + Npool_mobile(:,:)) * frac_green_aboveGround + &
                                   Npool_litter_green_ag(:,:) + Npool_litter_wood_ag(:,:) + &                
                                   Npool_woods(:,:) * frac_wood_2_atmos_fire * frac_wood_aboveGround) * &    
                                   cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:) 
        N_2_GreenLitterPools_tiled(:,:) = Npool_green(:,:) * (1._dp - frac_green_aboveGround) * &
                                              cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
        N_2_WoodLitterPools_tiled(:,:) = Cpool_woods(:,:) * (1._dp - frac_wood_2_atmos_fire * frac_wood_aboveGround) * &
                                             1._dp/cn_litter_wood * &
                                             cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
        N_2_sminn_tiled(:,:) = ( Npool_mobile(:,:)* (1._dp - frac_green_aboveGround) + &
                                  Cpool_woods(:,:)  * (1._dp - frac_wood_2_atmos_fire * frac_wood_aboveGround) * &
                                  (1._dp/cn_woods - 1._dp/cn_litter_wood) )                    * &
                                    cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
    end if

    !! sum carbon fluxes for all tiles

    carbon_2_GreenLitterPools(:) = carbon_2_GreenLitterPools(:) + SUM(carbon_2_GreenLitterPools_tiled(:,:),DIM=2)
    carbon_2_WoodLitterPools(:) = carbon_2_WoodLitterPools(:) + SUM(carbon_2_WoodLitterPools_tiled(:,:),DIM=2)
    carbon_2_atmos(:) = SUM(carbon_2_atmos_tiled(:,:),DIM=2)
    
    if (nitrogenMode) then
      nitrogen_2_GreenLitterPools(:) = nitrogen_2_GreenLitterPools(:) + SUM(N_2_GreenLitterPools_tiled(:,:),DIM=2)
      nitrogen_2_WoodLitterPools(:) = nitrogen_2_WoodLitterPools(:) + SUM(N_2_WoodLitterPools_tiled(:,:),DIM=2)
      nitrogen_2_atmos(:) = SUM(N_2_atmos_tiled(:,:),DIM=2)
      nitrogen_2_sminn(:) = SUM(N_2_sminn_tiled(:,:),DIM=2)
    end if 

    !! lower down the carbon density of above ground litter pools
    !! transfer carbon from the living tissue pools (green and reserve) to the below ground green litter pool
    !! transfer carbon from the wood pool to the above ground woody litter pool
    !! transfer carbon from the wood pool to the below ground woody litter pool
    IF (.NOT. with_yasso) THEN
       Cpool_litter_green_ag(:,:) = Cpool_litter_green_ag(:,:) * (1._dp - frac_burned(:,:))
       Cpool_litter_wood_ag (:,:) = Cpool_litter_wood_ag (:,:) * (1._dp - frac_burned(:,:))
       Cpool_litter_green_bg(:,:) = Cpool_litter_green_bg(:,:) + &
            (Cpool_green(:,:) + Cpool_reserve(:,:)) * (1._dp - frac_green_aboveGround) * frac_burned(:,:)
       Cpool_litter_wood_ag(:,:) = Cpool_litter_wood_ag(:,:) + &
            Cpool_woods(:,:) * ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround) * frac_burned(:,:)
       Cpool_litter_wood_bg(:,:) = Cpool_litter_wood_bg(:,:) + &
            Cpool_woods(:,:) * (1._dp - frac_wood_aboveGround) * frac_burned(:,:)
       if (nitrogenMode) then
          Npool_litter_green_ag(:,:) = Npool_litter_green_ag(:,:) * (1._dp - frac_burned(:,:))
          Npool_litter_wood_ag (:,:) = Npool_litter_wood_ag (:,:) * (1._dp - frac_burned(:,:))
          Npool_litter_green_bg(:,:) = Npool_litter_green_bg(:,:) + &
                                Npool_green(:,:) * (1._dp - frac_green_aboveGround) * frac_burned(:,:)
          Npool_litter_wood_ag(:,:) = Npool_litter_wood_ag(:,:) + &
                                      Cpool_woods(:,:) * (1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround * &
                                      1._dp/cn_litter_wood * frac_burned(:,:)
          Npool_litter_wood_bg(:,:) = Npool_litter_wood_bg(:,:) + &
                                      Cpool_woods(:,:) * (1._dp - frac_wood_aboveGround) * 1._dp/cn_litter_wood * &
                                      frac_burned(:,:)
          SMINN_pool(:,:) = SMINN_pool(:,:) + &
                            ( Npool_mobile(:,:)* (1._dp - frac_green_aboveGround) + &
                             Cpool_woods(:,:)  * (1._dp - frac_wood_2_atmos_fire * frac_wood_aboveGround) * &
                             (1._dp/cn_woods - 1._dp/cn_litter_wood) ) * frac_burned(:,:)                        
       end if 

    ELSE

       ! Lower down density of aboveground litter pools
       YCpool_acid_ag1      (:,:) = YCpool_acid_ag1       (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_water_ag1     (:,:) = YCpool_water_ag1      (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_ethanol_ag1   (:,:) = YCpool_ethanol_ag1    (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_nonsoluble_ag1(:,:) = YCpool_nonsoluble_ag1 (:,:) * (1._dp - frac_burned(:,:))   ! lower down density

       YCpool_acid_ag2      (:,:) = YCpool_acid_ag2       (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_water_ag2     (:,:) = YCpool_water_ag2      (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_ethanol_ag2   (:,:) = YCpool_ethanol_ag2    (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_nonsoluble_ag2(:,:) = YCpool_nonsoluble_ag2 (:,:) * (1._dp - frac_burned(:,:))   ! lower down density

       ! Add green and reserve vegetation carbon to belowground green litter pools
       YCpool_acid_bg1      (:,:) = YCpool_acid_bg1       (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,1)
       YCpool_water_bg1     (:,:) = YCpool_water_bg1      (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,2)
       YCpool_ethanol_bg1   (:,:) = YCpool_ethanol_bg1    (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,3)
       YCpool_nonsoluble_bg1(:,:) = YCpool_nonsoluble_bg1 (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,4)
       YCpool_humus_1       (:,:) = YCpool_humus_1        (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,5)

       ! Add the remaings of burned wood to the aboveground wood pools:
       YCpool_acid_ag2      (:,:) = YCpool_acid_ag2      (:,:)                                    & 
                                   +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround)  &
                                    * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,1)
       YCpool_water_ag2     (:,:) = YCpool_water_ag2                                              & 
                                   +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround)  &
                                    * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,2)
       YCpool_ethanol_ag2   (:,:) = YCpool_ethanol_ag2   (:,:)                                    & 
                                   +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround)  &
                                    * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,3)
       YCpool_nonsoluble_ag2(:,:) = YCpool_nonsoluble_ag2(:,:)                                    & 
                                    +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround) &
                                     * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,4)
       YCpool_humus_2       (:,:) = YCpool_humus_2       (:,:)                                    &
                                    +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround) &
                                     * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,5)

       ! Add the remaing of burned wood to the belowground wood pools:
       YCpool_acid_bg2      (:,:) = &
       YCpool_acid_bg2      (:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,1)
       YCpool_water_bg2     (:,:) = &
       YCpool_water_bg2     (:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,2)
       YCpool_ethanol_bg2   (:,:) = &
       YCpool_ethanol_bg2   (:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,3)
       YCpool_nonsoluble_bg2(:,:) = &
       YCpool_nonsoluble_bg2(:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,4)
       YCpool_humus_2       (:,:) = &
       YCpool_humus_2       (:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,5)
    END IF

    !! lower down the carbon density of living plant pools
    Cpool_green(:,:) = Cpool_green(:,:) * (1._dp - frac_burned(:,:))
    Cpool_woods(:,:) = Cpool_woods(:,:) * (1._dp - frac_burned(:,:))
    Cpool_reserve(:,:) = Cpool_reserve(:,:) * (1._dp - frac_burned(:,:))
     
    if (nitrogenMode) then
       Npool_green(:,:)  = Npool_green(:,:) * (1._dp - frac_burned(:,:))
       Npool_woods(:,:)  = Npool_woods(:,:) * (1._dp - frac_burned(:,:))
       Npool_mobile(:,:) = Npool_mobile(:,:) * (1._dp - frac_burned(:,:))
    end if

  end subroutine relocate_carbon_fire

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of damages to the vegetation (e.g. wind break) for the carbon pools. More precisely:
  ! It is assumed that for the damaged area the carbon from the living plant pools (Cpool_green, Cpool_reserve and Cpool_woods)
  ! is partly relocated into the litter pools.
  !
  subroutine relocate_carbon_damage(nidx, ntiles, with_yasso, frac_damaged, cf,             &
                                         veg_fract_correction,                              &
                                         Cpool_green, Cpool_woods, Cpool_reserve,           &
                                         carbon_2_GreenLitterPools,                         &
                                         carbon_2_WoodLitterPools,                          &            
                                         Cpool_litter_green_ag,Cpool_litter_green_bg,       &
                                         Cpool_litter_wood_ag, Cpool_litter_wood_bg,        &
                                         ! Yasso pools
                                         YCpool_acid_ag1,      YCpool_acid_ag2,             &
                                         YCpool_water_ag1,     YCpool_water_ag2,            &
                                         YCpool_ethanol_ag1,   YCpool_ethanol_ag2,          &
                                         YCpool_nonsoluble_ag1,YCpool_nonsoluble_ag2,       &
                                         YCpool_acid_bg1,      YCpool_acid_bg2,             &
                                         YCpool_water_bg1,     YCpool_water_bg2,            &
                                         YCpool_ethanol_bg1,   YCpool_ethanol_bg2,          &
                                         YCpool_nonsoluble_bg1,YCpool_nonsoluble_bg2,       &
                                         YCpool_humus_1,       YCpool_humus_2,              &
                                         LeafLitcoef,          WoodLitcoef,                 &
                                         !Nitrogen pools and fluxes
                                         Npool_green, Npool_woods, Npool_mobile, SMINN_pool,&
                                         Npool_litter_green_ag,Npool_litter_green_bg,       &
                                         Npool_litter_wood_ag, Npool_litter_wood_bg,        &
                                         nitrogen_2_GreenLitterPools,                       &
                                         nitrogen_2_WoodLitterPools,                        &
                                         nitrogen_2_sminn)

    USE mo_exception,        ONLY: finish

    integer,intent(in)     :: nidx                        !! Vector length
    integer,intent(in)     :: ntiles                      !! Number of tiles
    logical,intent(in)     :: with_yasso
    real(dp),intent(in)    :: frac_damaged(:,:)           !! Fraction of the vegetated area of each tile damaged till the last call
                                                          !!    of this routine
    real(dp),intent(in)    :: cf(:,:)                     !! Cover fractions
    real(dp),intent(in)    :: veg_fract_correction(:,:)   !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts for
                                                          !!    sparseness of vegetation)
    real(dp),intent(inout) :: Cpool_green(:,:)            !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)            !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)          !! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(out)   :: carbon_2_GreenLitterPools(:)!! Amount of carbon relocated by wind damage
                                                          !! .. to the green litter pools [mol(C)/m^2(vegetated area)]
    real(dp),intent(out)   :: carbon_2_WoodLitterPools(:) !! Amount of carbon relocated by wind damage and fire
    
    real(dp),intent(inout) :: Cpool_litter_green_ag(:,:)  !! Value of above ground green litter carbon pool 
                                                          !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_green_bg(:,:)  !! Value of below ground green litter carbon pool 
                                                          !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood_ag(:,:)   !! Value of wood litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood_bg(:,:)   !! Value of wood litter carbon pool [mol(C)/m^2(canopy)]
    !   YASSO 
    ! Size class 1: green
    !   above ground C pools
    real(dp),intent(inout) :: YCpool_acid_ag1(:,:)        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_water_ag1(:,:)       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_ethanol_ag1(:,:)     !! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(dp),intent(inout) :: YCpool_nonsoluble_ag1(:,:)  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout) :: YCpool_acid_bg1(:,:)        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_water_bg1(:,:)       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_ethanol_bg1(:,:)     !! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(dp),intent(inout) :: YCpool_nonsoluble_bg1(:,:)  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_humus_1(:,:)         !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2: wood 
    !   above ground C pools
    real(dp),intent(inout) :: YCpool_acid_ag2(:,:)        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_water_ag2(:,:)       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_ethanol_ag2(:,:)     !! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(dp),intent(inout) :: YCpool_nonsoluble_ag2(:,:)  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout) :: YCpool_acid_bg2(:,:)        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_water_bg2(:,:)       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_ethanol_bg2(:,:)     !! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(dp),intent(inout) :: YCpool_nonsoluble_bg2(:,:)  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: YCpool_humus_2(:,:)         !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    !   parameters
    real(dp),intent(in)    :: LeafLitcoef(:,:,:)          !! Factor to spread non woody litter into yasso pools [ ]
    real(dp),intent(in)    :: WoodLitcoef(:,:,:)          !! Factor to spread woody litterfall into yasso pools [ ]

    real(dp),intent(inout),optional  :: Npool_green(:,:)                !! Value of green nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional  :: Npool_woods(:,:)                !! Value of wood nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional  :: Npool_mobile(:,:)               !! Value of mobile nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional  :: SMINN_pool(:,:)                 !! Value of mineral nitrogen pool [mol(N)/m^(canopy)2]   
    real(dp),intent(inout),optional  :: Npool_litter_green_ag(:,:)      !! Value of above ground green litter nitrogen pool
                                                                        !!    [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional  :: Npool_litter_green_bg(:,:)      !! Value of below ground green litter nitrogen pool
                                                                        !!    [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional  :: Npool_litter_wood_ag(:,:)       !! Value of above ground wood litter nitrogen pool
                                                                        !!    [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional  :: Npool_litter_wood_bg(:,:)       !! Value of below wood litter nitrogen pool
                                                                        !!    [mol(N)/m^(canopy)2]
    real(dp),intent(out),optional    :: nitrogen_2_GreenLitterPools(:)  !! Amount of nitrogen relocated by wind damage
                                                                        !! .. to the green litter pools [mol(N)/m^2(vegetated area)]
    real(dp),intent(out),optional    :: nitrogen_2_WoodLitterPools(:)   !! Amount of nitrogen relocated by wind damage
                                                                        !! .. to the wood litter pools [mol(N)/m^2(vegetated area)]
    real(dp),intent(out),optional    :: nitrogen_2_sminn(:)             !! Amount of nitrogen relocated by wind damage
                                                                        !! ... to the mineral N pool [mol(N)/m2(vegetated area) ]
    !! locals

    real(dp)  ::  carbon_2_WoodLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  carbon_2_GreenLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  N_2_WoodLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  N_2_GreenLitterPools_tiled(nidx,ntiles) 
    real(dp)  ::  N_2_SMINN_tiled(nidx,ntiles)
    real(dp)  ::  Cpool_temp(nidx,ntiles)
    logical   ::  nitrogenMode
    
    nitrogenMode=.false.
    if (present(Npool_green)          .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(SMINN_pool) ) then
       nitrogenMode = .true.
       if (.not. (present(Npool_green)           .and. present(Npool_woods)           .and. &
                  present(Npool_mobile)          .and. present(Npool_litter_green_ag) .and. &
                  present(Npool_litter_wood_ag)  .and. present(Npool_litter_wood_bg)  .and. &
                  present(Npool_litter_green_bg) .and. present(SMINN_pool) ) )&
            CALL finish('relocate_carbon_damage()','at least one variable missing to handle nitrogen dynamics')
    end if    

    !! preparations

    carbon_2_WoodLitterPools_tiled(:,:) = 0._dp
    carbon_2_greenLitterPools_tiled(:,:) = 0._dp
    N_2_WoodLitterPools_tiled(:,:) = 0._dp
    N_2_greenLitterPools_tiled(:,:) = 0._dp  
    N_2_SMINN_tiled(:,:) = 0._dp  

    !! transfer carbon from the wood pool to the wood litter pools nd
    !! transfer carbon from the green and reserve pool to the green litter pool
    !! determine amount of carbon relocated

    carbon_2_WoodLitterPools_tiled(:,:) = Cpool_woods(:,:) &
                                          * cf(:,:) * veg_fract_correction(:,:) * frac_damaged(:,:)
    carbon_2_GreenLitterPools_tiled(:,:) = (Cpool_green(:,:) + Cpool_reserve(:,:)) &
                                          * cf(:,:) * veg_fract_correction(:,:) * frac_damaged(:,:)
    
    if (nitrogenMode) then
!dk Npool_mobile assigned to litter green pool    
       Npool_litter_wood_bg(:,:) = Npool_litter_wood_bg(:,:) &
                               + Cpool_woods(:,:) * (1._dp - frac_wood_aboveGround)* 1._dp/cn_litter_wood  &
                               * frac_damaged(:,:)
       Npool_litter_wood_ag(:,:) = Npool_litter_wood_ag(:,:) &
                               + Cpool_woods(:,:) * frac_wood_aboveGround * 1._dp/cn_litter_wood  &
                               * frac_damaged(:,:)
       Npool_litter_green_bg(:,:) = Npool_litter_green_bg(:,:) &
                                + (Npool_green(:,:) + Npool_mobile(:,:)) * (1._dp - frac_green_aboveGround) &
                                * frac_damaged(:,:)
       Npool_litter_green_ag(:,:) = Npool_litter_green_ag(:,:) &
                                + (Npool_green(:,:) + Npool_mobile(:,:)) * frac_green_aboveGround &
                                * frac_damaged(:,:)
       SMINN_pool(:,:) = SMINN_pool(:,:) &
                               + Cpool_woods(:,:) * (1._dp/cn_woods - 1._dp/cn_litter_wood)     &
                               * frac_damaged(:,:)
       N_2_WoodLitterPools_tiled(:,:) = Cpool_woods(:,:) * 1._dp/cn_litter_wood &
                                           * cf(:,:) * veg_fract_correction(:,:) * frac_damaged(:,:)
       N_2_GreenLitterPools_tiled(:,:) = (Npool_green(:,:) + Npool_mobile(:,:)) &
                                           * cf(:,:) * veg_fract_correction(:,:) * frac_damaged(:,:)
       N_2_sminn_tiled(:,:) = Cpool_woods(:,:) * (1._dp/cn_woods - 1._dp/cn_litter_wood) * frac_damaged(:,:)
    end if

    IF (.NOT. with_yasso) THEN
       Cpool_litter_wood_bg(:,:) = Cpool_litter_wood_bg(:,:) &
                                 + Cpool_woods(:,:) * (1._dp - frac_wood_aboveGround)*frac_damaged(:,:)
       Cpool_litter_wood_ag(:,:) = Cpool_litter_wood_ag(:,:) &
                                 + Cpool_woods(:,:) * frac_wood_aboveGround * frac_damaged(:,:)
       Cpool_litter_green_bg(:,:) = Cpool_litter_green_bg(:,:) &
                                  + (Cpool_green(:,:) + Cpool_reserve(:,:)) * (1._dp - frac_green_aboveGround) &
                                  * frac_damaged(:,:)
       Cpool_litter_green_ag(:,:) = Cpool_litter_green_ag(:,:) &
                                  + (Cpool_green(:,:) + Cpool_reserve(:,:)) * frac_green_aboveGround &
                                  * frac_damaged(:,:)
    ELSE

       ! Add aboveground damaged green and reserve
       Cpool_temp(:,:) = Cpool_green(:,:) + Cpool_reserve(:,:)
       YCpool_acid_ag1      (:,:) = &
       YCpool_acid_ag1      (:,:) + frac_green_aboveGround * frac_damaged(:,:) * LeafLitcoef(:,:,1) * Cpool_temp(:,:)
       YCpool_water_ag1     (:,:) = &
       YCpool_water_ag1     (:,:) + frac_green_aboveGround * frac_damaged(:,:) * LeafLitcoef(:,:,2) * Cpool_temp(:,:)
       YCpool_ethanol_ag1   (:,:) = &
       YCpool_ethanol_ag1   (:,:) + frac_green_aboveGround * frac_damaged(:,:) * LeafLitcoef(:,:,3) * Cpool_temp(:,:)
       YCpool_nonsoluble_ag1(:,:) = &
       YCpool_nonsoluble_ag1(:,:) + frac_green_aboveGround * frac_damaged(:,:) * LeafLitcoef(:,:,4) * Cpool_temp(:,:)

       ! Add belowground damaged green and reserve
       YCpool_acid_bg1      (:,:) = &
       YCpool_acid_bg1      (:,:) + (1._dp - frac_green_aboveGround) * frac_damaged(:,:) * LeafLitcoef(:,:,1) * Cpool_temp(:,:)
       YCpool_water_bg1     (:,:) = &
       YCpool_water_bg1     (:,:) + (1._dp - frac_green_aboveGround) * frac_damaged(:,:) * LeafLitcoef(:,:,2) * Cpool_temp(:,:)
       YCpool_ethanol_bg1   (:,:) = &
       YCpool_ethanol_bg1   (:,:) + (1._dp - frac_green_aboveGround) * frac_damaged(:,:) * LeafLitcoef(:,:,3) * Cpool_temp(:,:)
       YCpool_nonsoluble_bg1(:,:) = &
       YCpool_nonsoluble_bg1(:,:) + (1._dp - frac_green_aboveGround) * frac_damaged(:,:) * LeafLitcoef(:,:,4) * Cpool_temp(:,:)
       YCpool_humus_1       (:,:) = &
       YCpool_humus_1       (:,:) +                                    frac_damaged(:,:) * LeafLitcoef(:,:,5) * Cpool_temp(:,:)

       ! Add aboveground damaged wood
       YCpool_acid_ag2      (:,:) = &
       YCpool_acid_ag2      (:,:) + frac_wood_aboveGround  * frac_damaged(:,:) * WoodLitcoef(:,:,1) * Cpool_woods(:,:)
       YCpool_water_ag2     (:,:) = &
       YCpool_water_ag2     (:,:) + frac_wood_aboveGround  * frac_damaged(:,:) * WoodLitcoef(:,:,2) * Cpool_woods(:,:)
       YCpool_ethanol_ag2   (:,:) = &
       YCpool_ethanol_ag2   (:,:) + frac_wood_aboveGround  * frac_damaged(:,:) * WoodLitcoef(:,:,3) * Cpool_woods(:,:)
       YCpool_nonsoluble_ag2(:,:) = &
       YCpool_nonsoluble_ag2(:,:) + frac_wood_aboveGround  * frac_damaged(:,:) * WoodLitcoef(:,:,4) * Cpool_woods(:,:)

       ! Add belowground damaged wood
       YCpool_acid_bg2      (:,:) = &
       YCpool_acid_bg2      (:,:) + (1._dp -frac_wood_aboveGround)  * frac_damaged(:,:) * WoodLitcoef(:,:,1) * Cpool_woods(:,:)
       YCpool_water_bg2     (:,:) = &
       YCpool_water_bg2     (:,:) + (1._dp -frac_wood_aboveGround)  * frac_damaged(:,:) * WoodLitcoef(:,:,2) * Cpool_woods(:,:)
       YCpool_ethanol_bg2   (:,:) = &
       YCpool_ethanol_bg2   (:,:) + (1._dp -frac_wood_aboveGround)  * frac_damaged(:,:) * WoodLitcoef(:,:,3) * Cpool_woods(:,:)
       YCpool_nonsoluble_bg2(:,:) = &
       YCpool_nonsoluble_bg2(:,:) + (1._dp -frac_wood_aboveGround)  * frac_damaged(:,:) * WoodLitcoef(:,:,4) * Cpool_woods(:,:)
       YCpool_humus_2       (:,:) = &
       YCpool_humus_2       (:,:) +                                   frac_damaged(:,:) * WoodLitcoef(:,:,5) * Cpool_woods(:,:)

    END IF

    !! sum carbon fluxes for all tiles

    carbon_2_WoodLitterPools(:) = SUM(carbon_2_WoodLitterPools_tiled(:,:),DIM=2)
    carbon_2_GreenLitterPools(:) = SUM(carbon_2_GreenLitterPools_tiled(:,:),DIM=2)
    if (nitrogenMode) then
      nitrogen_2_WoodLitterPools(:) = SUM(N_2_WoodLitterPools_tiled(:,:),DIM=2)
      nitrogen_2_GreenLitterPools(:) = SUM(N_2_GreenLitterPools_tiled(:,:),DIM=2) 
      nitrogen_2_sminn(:) = SUM(N_2_sminn_tiled(:,:),DIM=2)   
    end if

    !! lower down the carbon density of living plant pools
    Cpool_green(:,:) = Cpool_green(:,:) * (1._dp - frac_damaged(:,:))
    Cpool_woods(:,:) = Cpool_woods(:,:) * (1._dp - frac_damaged(:,:))
    Cpool_reserve(:,:) = Cpool_reserve(:,:) * (1._dp - frac_damaged(:,:))
    
    if (nitrogenMode) then
      Npool_green(:,:) = Npool_green(:,:) * (1._dp - frac_damaged(:,:))
      Npool_woods(:,:) = Npool_woods(:,:) * (1._dp - frac_damaged(:,:))
      Npool_mobile(:,:) = Npool_mobile(:,:) * (1._dp - frac_damaged(:,:))        
    end if

  end subroutine relocate_carbon_damage

  ! --- C_relocation_from_LUtransitions() -----------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of landcover transitions, i.e. here the landcover changes is given in form of a transition
  ! matrix. Otherwise the logic of this routine is quite similar to the routine relocate_carbonAndNitrogen() from above.
  !

  SUBROUTINE C_relocation_from_LUtransitions(lctlib,surface,nidx,ntiles,is_vegetation,                &
                                             cf_old, cf_new, Tile_TransMtrx, veg_fract_correction,    &
                                             frac_wood_2_atmos, frac_green_2_atmos,                   &
                                             frac_reserve_2_atmos,                                    &
                                             C_2_atmos,                                               &
                                             Cpool_green, Cpool_woods, Cpool_reserve,                 &
                                             Cpool_litter_green_ag,Cpool_litter_green_bg,             &
                                             Cpool_litter_wood_ag,Cpool_litter_wood_bg,               &
                                             Cpool_slow,                                              &
                                             !for N option
                                             Npool_green,Npool_woods,Npool_mobile,Npool_litter_green_ag, &
                                             Npool_litter_green_bg,Npool_litter_wood_ag,Npool_litter_wood_bg, &
                                             Npool_slow,SMINN_pool,Nitrogen_2_atmos,frac_mobile_2_atmos, &    
                                             !yasso
                                             YCpool_acid_ag1, YCpool_acid_bg1,                          &
                                             YCpool_water_ag1, YCpool_water_bg1,                        &
                                             YCpool_ethanol_ag1, YCpool_ethanol_bg1,                    &
                                             YCpool_nonsoluble_ag1, YCpool_nonsoluble_bg1,              &
                                             YCpool_humus_1,                                            &
                                             YCpool_acid_ag2, YCpool_acid_bg2,                          &
                                             YCpool_water_ag2, YCpool_water_bg2,                        &
                                             YCpool_ethanol_ag2, YCpool_ethanol_bg2,                    &
                                             YCpool_nonsoluble_ag2, YCpool_nonsoluble_bg2,              &
                                             YCpool_humus_2,                                            &
                                             LeafLitcoef, WoodLitcoef,                                &
                                             ! LCC fluxes after Houghton
                                             C_2_litterGreenPools,                                    &
                                             C_2_litterWoodPool_ag,C_2_litterWoodPool_bg,             &         
                                             Cpool_onSite,Cpool_paper,Cpool_construction,             &
                                             C_onSite_2_atmos,C_paper_2_atmos,C_construction_2_atmos, &
                                             C_2_onSite,C_2_paper,C_2_construction,                   &                             
                                             lcc_scheme)   

    USE mo_jsbach_lctlib,   ONLY: lctlib_type
    USE mo_land_surface,    ONLY: land_surface_type
    USE mo_exception,       ONLY: finish

    type(lctlib_type),       intent(in) :: lctlib           !! PFT-specific constants
    type(land_surface_type), intent(in) :: surface
    integer,  intent(in)   :: lcc_scheme                    !! scheme for land cover change 
    integer, intent(in)    :: nidx                          !! Number of gridpoints in vectors
    integer, intent(in)    :: ntiles                        !! Number of tiles
    logical, intent(in)    :: is_vegetation(:,:)            !! Logical mask indcating vegetated tiles
    real(dp),intent(in)    :: cf_old(:,:)                   !! Cover fractions before landcover change or vegetation dynamics [] 
    real(dp),intent(in)    :: cf_new(:,:)                   !! Cover fraction after landcover change or vegetation dynamics []
    real(dp),intent(in)    :: Tile_TransMtrx(:,:,:)         !! Transition matrix between all tiles describing the landcover change
                                                            !!    (noOfgridPoints x noOfTiles x noOfTiles)
    real(dp),intent(in)    :: veg_fract_correction(:,:)     !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts for
                                                            !!    sparseness of vegetation)
    real(dp),intent(in)    :: frac_wood_2_atmos             !! Fraction of wood pool immediately emitted to atmosphere (e.g. by 
                                                            !!    burning vegetation down)
    real(dp),intent(in)    :: frac_green_2_atmos            !! Fraction of green pool immediately emitted to atmosphere (e.g. by 
                                                            !!    burning vegetation down)
    real(dp),intent(in)    :: frac_reserve_2_atmos          !! Fraction of reserve pool immediately emitted to atmosphere (e.g. by
                                                            !!    burning vegetation down)
    real(dp),intent(out)   :: C_2_atmos(:)                  !! Amount of carbon directly emitted to atmosphere in this timestep 
                                                            !!    [mol(C)/m^2(vegetated area)]
    real(dp),intent(inout) :: Cpool_green(:,:)              !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)              !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)            !! Value of reserve carbon pool [mol(C)/m^(canopy)2]

    ! litter and soil C pools of old cbalance scheme
    real(dp),optional,intent(inout) :: Cpool_litter_green_ag(:,:)      !! above ground green litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),optional,intent(inout) :: Cpool_litter_green_bg(:,:)      !! below ground green litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),optional,intent(inout) :: Cpool_litter_wood_ag(:,:)       !! above ground wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),optional,intent(inout) :: Cpool_litter_wood_bg(:,:)       !! below ground wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),optional,intent(inout) :: Cpool_slow(:,:)                 !! slow soil carbon pool [mol(C)/m^2(canopy)]

    ! Variables only needed for nitrogen option
    real(dp),intent(inout),optional :: Npool_green(:,:)
    real(dp),intent(inout),optional :: Npool_woods(:,:)
    real(dp),intent(inout),optional :: Npool_mobile(:,:)
    real(dp),intent(inout),optional :: Npool_litter_green_ag(:,:)
    real(dp),intent(inout),optional :: Npool_litter_green_bg(:,:)
    real(dp),intent(inout),optional :: Npool_litter_wood_ag(:,:)
    real(dp),intent(inout),optional :: Npool_litter_wood_bg(:,:)
    real(dp),intent(inout),optional :: Npool_slow(:,:)
    real(dp),intent(inout),optional :: SMINN_pool(:,:)
    real(dp),intent(out),optional :: Nitrogen_2_atmos(:)
    real(dp),intent(in),optional  :: frac_mobile_2_atmos

    ! variables only needed with yasso
    ! size class 1: green
    !   above ground C pools  
    real(dp),optional,intent(inout) :: YCpool_acid_ag1(:,:)             !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_ag1(:,:)            !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_ag1(:,:)          !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_ag1(:,:)       !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),optional,intent(inout) :: YCpool_acid_bg1(:,:)             !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_bg1(:,:)            !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_bg1(:,:)          !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_bg1(:,:)       !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_humus_1(:,:)               !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! size class 2: wood 
    !   above ground C pools  
    real(dp),optional,intent(inout) :: YCpool_acid_ag2(:,:)             !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_ag2(:,:)            !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_ag2(:,:)          !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_ag2(:,:)       !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),optional,intent(inout) :: YCpool_acid_bg2(:,:)             !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_bg2(:,:)            !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_bg2(:,:)          !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_bg2(:,:)       !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_humus_2(:,:)               !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    !   parameters 
    real(dp),optional,intent(in)    :: LeafLitcoef(:,:,:)              !! Factor to spread non woody litterfall into yasso pools [ ]
    real(dp),optional,intent(in)    :: WoodLitcoef(:,:,:)              !! Factor to spread woody litterfall into yasso pools [ ]

    ! variables only needed with the original landuse transitions scheme (lcc_scheme==1)
    real(dp),optional,intent(out)   :: C_2_litterGreenPools(:)         !! Amount of carbon relocated by land cover change from green
                                                                       !!    and reserve pool to below and above ground green litter
    real(dp),optional,intent(out)   :: C_2_litterWoodPool_ag(:)        !! Amount of carbon relocated by land cover change from wood
                                                                       !!    pool to above ground woody litter pool
                                                                       !!    [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(out)   :: C_2_litterWoodPool_bg(:)        !! Amount of carbon relocated by land cover change from wood
                                                                       !!    pool to below ground woody litter pool
                                                                       !!    [mol(C)/m^2(vegetated area)]

    ! variables only needed with the landuse transitions scheme with antropogenic pools (lcc_scheme 2 or 3)
    real(dp),optional,intent(inout) ::  Cpool_onSite(:)                !! Amount of carbon remains in anthro annual pool
                                                                        !!   [mol(C)/m^2(vegetated area)] 
    real(dp),optional,intent(inout) ::  Cpool_paper(:)                  !! Amount of carbon remains in anthro decadal pool
    real(dp),optional,intent(inout) ::  Cpool_construction(:)           !! Amount of carbon remains in anthro centinnial pool
    real(dp),optional,intent(out)   ::  C_onSite_2_atmos(:)    !! C flux from anthro annual pool to atmos [mol(C)/m^2(veg)/day] 
    real(dp),optional,intent(out)   ::  C_paper_2_atmos(:)              !! Carbon flux from green anthro annual pool to atmos
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_construction_2_atmos(:)       !! Carbon flux from woody anthro annual pool to atmos
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_2_onSite(:)                   !! C flux to anthro annual pool [mol(C)/m^2(veg area)/day]
    real(dp),optional,intent(out)   ::  C_2_construction(:)             !! Carbon flux from woody anthro annual pool to paper
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_2_paper(:)                    !! Carbon flux from woody anthro annual pool to
                                                                        !!   construction [mol(C)/m^2(vegetated area)/day]
    
    !! locals

    real(dp) :: Cpool_green_loss(1:nidx,1:ntiles)    !! Carbon loss of green pool for each tile (per corrected vegetated area)
    real(dp) :: Cpool_reserve_loss(1:nidx,1:ntiles)  !! Carbon loss of reserve pool for each tile (per corrected vegetated area)
    real(dp) :: Cpool_woods_loss(1:nidx,1:ntiles)    !! Carbon loss of wood pool for each tile (per corrected vegetated area)
    real(dp) :: Cpool_litter_green_ag_loss(1:nidx,1:ntiles)
    real(dp) :: Cpool_litter_wood_ag_loss(1:nidx,1:ntiles)      
    real(dp) :: YCpool_acid_ag1_loss(1:nidx,1:ntiles)      
    real(dp) :: YCpool_water_ag1_loss(1:nidx,1:ntiles)      
    real(dp) :: YCpool_ethanol_ag1_loss(1:nidx,1:ntiles)      
    real(dp) :: YCpool_nonsoluble_ag1_loss(1:nidx,1:ntiles)      
    real(dp) :: YCpool_acid_ag2_loss(1:nidx,1:ntiles)      
    real(dp) :: YCpool_water_ag2_loss(1:nidx,1:ntiles)      
    real(dp) :: YCpool_ethanol_ag2_loss(1:nidx,1:ntiles)      
    real(dp) :: YCpool_nonsoluble_ag2_loss(1:nidx,1:ntiles)      
    real(dp) :: scale_fac(1:nidx,1:ntiles)

    integer  :: n, i
    real(dp) :: LitterGreen(1:nidx,1:ntiles)    !! Carbon gain of green litter pools (below and above ground)    
    real(dp) :: LitterWood(1:nidx,1:ntiles)     !! Carbon gain of wood litter pools (below and above ground)
    real(dp) :: Npool_green_loss(1:nidx,1:ntiles) 
    real(dp) :: Npool_wood_loss(1:nidx,1:ntiles) 
    real(dp) :: Npool_mobile_loss(1:nidx,1:ntiles)
    real(dp) :: LitterGreenN(1:nidx,1:ntiles)
    real(dp) :: LitterWoodN(1:nidx,1:ntiles)
    real(dp) :: SMINNgain(1:nidx,1:ntiles)
   
    logical  :: with_yasso
    logical :: nitrogenMode

    !! --- GO! ---------------------------------------

    !! Check for nitrogen mode
    nitrogenMode = .false.
    if(present(Npool_green)           .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(Npool_slow)            .or. &
       present(SMINN_pool) ) then
       nitrogenMode = .true.
       if (.not. (present(Npool_green)           .and. present(Npool_woods)           .and. &
                  present(Npool_mobile)          .and. present(Npool_litter_green_ag) .and. &
                  present(Npool_litter_wood_ag)  .and. present(Npool_litter_wood_bg)  .and. &
                  present(Npool_litter_green_bg) .and. present(Npool_slow)            .and. &
                  present(SMINN_pool) ) ) &
              CALL finish('relocate_CarbonAndNitrogen()','at least one variable missing to handle nitrogen dynamics') 
    end if

    !! Check if yasso or cbalance litter and soil pools should be handled
    with_yasso  = .FALSE.  
    IF (PRESENT(YCpool_acid_ag1)        .OR. PRESENT(YCpool_acid_bg1)            .OR. &
        PRESENT(YCpool_water_ag1)       .OR. PRESENT(YCpool_water_bg1)           .OR. &
        PRESENT(YCpool_ethanol_ag1)     .OR. PRESENT(YCpool_ethanol_bg1)         .OR. &
        PRESENT(YCpool_nonsoluble_ag1)  .OR. PRESENT(YCpool_nonsoluble_bg1)      .OR. &
        PRESENT(YCpool_humus_1)         .OR.                                          &
        PRESENT(YCpool_acid_ag2)        .OR. PRESENT(YCpool_acid_bg2)            .OR. &
        PRESENT(YCpool_water_ag2)       .OR. PRESENT(YCpool_water_bg2)           .OR. &
        PRESENT(YCpool_ethanol_ag2)     .OR. PRESENT(YCpool_ethanol_bg2)         .OR. &
        PRESENT(YCpool_nonsoluble_ag2)  .OR. PRESENT(YCpool_nonsoluble_bg2)      .OR. &
        PRESENT(YCpool_humus_2)         .OR. PRESENT(LeafLitcoef)                .OR. &
         PRESENT(WoodLitcoef)) THEN
        with_yasso=.TRUE.
        IF (.NOT. (PRESENT(YCpool_acid_ag1)        .AND. PRESENT(YCpool_acid_bg1)        .AND. &
                   PRESENT(YCpool_water_ag1)       .AND. PRESENT(YCpool_water_bg1)       .AND. &
                   PRESENT(YCpool_ethanol_ag1)     .AND. PRESENT(YCpool_ethanol_bg1)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag1)  .AND. PRESENT(YCpool_nonsoluble_bg1)  .AND. &
                   PRESENT(YCpool_humus_1)         .AND.                                       &
                   PRESENT(YCpool_acid_ag2)        .AND. PRESENT(YCpool_acid_bg2)        .AND. &
                   PRESENT(YCpool_water_ag2)       .AND. PRESENT(YCpool_water_bg2)       .AND. &
                   PRESENT(YCpool_ethanol_ag2)     .AND. PRESENT(YCpool_ethanol_bg2)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag2)  .AND. PRESENT(YCpool_nonsoluble_bg2)  .AND. &
                   PRESENT(YCpool_humus_2)          .AND. PRESENT(LeafLitcoef)           .AND. &
                   PRESENT(WoodLitcoef) ) &
            ) THEN
            CALL finish('relocate_carbon_desert()','at least one variable missing to handle yasso pools')
        END IF
    END IF
    
    ! Scaling factor for new cover fractions where tiles are shrinking
    FORALL(i=1:ntiles)
       scale_fac(:,i) = veg_fract_correction(:,i) * cf_old(:,i) * (1._dp-Tile_TransMtrx(:,i,i))
    END FORALL

    ! New sizes of plant pools
    ! Calculate carbon loss of all pools and update anthro pools if requested
    IF (.NOT. with_yasso) THEN
      if (nitrogenMode) then
       CALL C_loss_and_update_anthro_pools(nidx, ntiles, lcc_scheme, lctlib, surface,      &
                                              Cpool_green, Cpool_woods, Cpool_reserve,     &
              Cpool_litter_green_ag        =  Cpool_litter_green_ag,                       &
              Cpool_litter_wood_ag         =  Cpool_litter_wood_ag,                        &
              Cpool_onSite                 =  Cpool_onSite,                                &
              Cpool_paper                  =  Cpool_paper,                                 &
              Cpool_construction           =  Cpool_construction,                          &
              scale_fac                    =  scale_fac,                                   &
              C_2_onSite                   =  C_2_onSite,                                  &
              C_2_paper                    =  C_2_paper,                                   &
              C_2_construction             =  C_2_construction,                            &
              C_onSite_2_atmos             =  C_onSite_2_atmos,                            &
              C_paper_2_atmos              =  C_paper_2_atmos,                             &
              C_construction_2_atmos       =  C_construction_2_atmos,                      &
              C_2_atmos                    =  C_2_atmos,                                   &
              C_green_loss                 =  Cpool_green_loss,                            &
              C_wood_loss                  =  Cpool_woods_loss,                            &
              C_reserve_loss               =  Cpool_reserve_loss,                          &
              C_litter_green_ag_loss       =  Cpool_litter_green_ag_loss,                  &
              C_litter_wood_ag_loss        =  Cpool_litter_wood_ag_loss,                   &
              Npool_green                  =  Npool_green,                                 &
              Npool_woods                  =  Npool_woods,                                 &
              Npool_mobile                 =  Npool_mobile,                                &
              Npool_green_loss             =  Npool_green_loss,                            &
              Npool_wood_loss              =  Npool_wood_loss,                             &
              Npool_mobile_loss            =  Npool_mobile_loss)
    
      else

       CALL C_loss_and_update_anthro_pools(nidx, ntiles, lcc_scheme, lctlib, surface,      &
                                              Cpool_green, Cpool_woods, Cpool_reserve,     &
              Cpool_litter_green_ag        =  Cpool_litter_green_ag,                       &
              Cpool_litter_wood_ag         =  Cpool_litter_wood_ag,                        &
              Cpool_onSite                 =  Cpool_onSite,                                &
              Cpool_paper                  =  Cpool_paper,                                 &
              Cpool_construction           =  Cpool_construction,                          &
              scale_fac                    =  scale_fac,                                   &
              C_2_onSite                   =  C_2_onSite,                                  &
              C_2_paper                    =  C_2_paper,                                   &
              C_2_construction             =  C_2_construction,                            &
              C_onSite_2_atmos             =  C_onSite_2_atmos,                            &
              C_paper_2_atmos              =  C_paper_2_atmos,                             &
              C_construction_2_atmos       =  C_construction_2_atmos,                      &
              C_2_atmos                    =  C_2_atmos,                                   &
              C_green_loss                 =  Cpool_green_loss,                            & 
              C_wood_loss                  =  Cpool_woods_loss,                            &
              C_reserve_loss               =  Cpool_reserve_loss,                          &
              C_litter_green_ag_loss       =  Cpool_litter_green_ag_loss,                  &
              C_litter_wood_ag_loss        =  Cpool_litter_wood_ag_loss)
      end if ! nitrogenMode
    ELSE  !! with yasso
       CALL C_loss_and_update_anthro_pools(nidx, ntiles, lcc_scheme, lctlib, surface,      &
                                              Cpool_green, Cpool_woods, Cpool_reserve,     &
              YCpool_acid_ag1              =  YCpool_acid_ag1,                             &
              YCpool_water_ag1             =  YCpool_water_ag1,                            &
              YCpool_ethanol_ag1           =  YCpool_ethanol_ag1,                          &
              YCpool_nonsoluble_ag1        =  YCpool_nonsoluble_ag1,                       &
              YCpool_acid_ag2              =  YCpool_acid_ag2,                             &
              YCpool_water_ag2             =  YCpool_water_ag2,                            &
              YCpool_ethanol_ag2           =  YCpool_ethanol_ag2,                          &
              YCpool_nonsoluble_ag2        =  YCpool_nonsoluble_ag2,                       &
              Cpool_onSite                 =  Cpool_onSite,                                &
              Cpool_paper                  =  Cpool_paper,                                 &
              Cpool_construction           =  Cpool_construction,                          &
              scale_fac                    =  scale_fac,                                   &
              C_2_onSite                   =  C_2_onSite,                                  &
              C_2_paper                    =  C_2_paper,                                   &
              C_2_construction             =  C_2_construction,                            &
              C_onSite_2_atmos             =  C_onSite_2_atmos,                            &
              C_paper_2_atmos              =  C_paper_2_atmos,                             &
              C_construction_2_atmos       =  C_construction_2_atmos,                      &
              C_2_atmos                    =  C_2_atmos,                                   &
              C_green_loss                 =  Cpool_green_loss,                            &
              C_wood_loss                  =  Cpool_woods_loss,                            &
              C_reserve_loss               =  Cpool_reserve_loss,                          &
              C_acid_ag1_loss               =  YCpool_acid_ag1_loss,                       &
              C_water_ag1_loss              =  YCpool_water_ag1_loss,                      &
              C_ethanol_ag1_loss            =  YCpool_ethanol_ag1_loss,                    &
              C_nonsoluble_ag1_loss         =  YCpool_nonsoluble_ag1_loss,                 &
              C_acid_ag2_loss               =  YCpool_acid_ag2_loss,                       &
              C_water_ag2_loss              =  YCpool_water_ag2_loss,                      &
              C_ethanol_ag2_loss            =  YCpool_ethanol_ag2_loss,                    &
              C_nonsoluble_ag2_loss         =  YCpool_nonsoluble_ag2_loss)
    END IF

    !! Update main C-pools 
    WHERE(is_vegetation(:,:))
       Cpool_green            =   (cf_old * veg_fract_correction * Cpool_green           - Cpool_green_loss           ) &
                                / (cf_new * veg_fract_correction)
       Cpool_woods            =   (cf_old * veg_fract_correction * Cpool_woods           - Cpool_woods_loss           ) &
                                / (cf_new * veg_fract_correction)
       Cpool_reserve          =   (cf_old * veg_fract_correction * Cpool_reserve         - Cpool_reserve_loss         ) &
                                / (cf_new * veg_fract_correction)
    END WHERE

    if (nitrogenMode) then
       where(is_vegetation)
           Npool_green = ( cf_old * veg_fract_correction*Npool_green   - Npool_green_loss )   &
                                          / ( cf_new*veg_fract_correction )
           Npool_woods = ( cf_old * veg_fract_correction*Npool_woods   - Npool_wood_loss )   &
                                          / ( cf_new*veg_fract_correction )
           Npool_mobile = ( cf_old * veg_fract_correction*Npool_mobile   - Npool_mobile_loss )   &
                                          / ( cf_new*veg_fract_correction )
       end where
    end if

    IF (lcc_scheme == 2) THEN
       IF (.NOT. with_yasso) THEN
          WHERE(is_vegetation(:,:))
             Cpool_litter_green_ag  =   (cf_old * veg_fract_correction * Cpool_litter_green_ag - Cpool_litter_green_ag_loss ) &
                                      / (cf_new * veg_fract_correction)
             Cpool_litter_wood_ag   =   (cf_old * veg_fract_correction * Cpool_litter_wood_ag  - Cpool_litter_wood_ag_loss  ) &
                                      / (cf_new * veg_fract_correction)
          END WHERE
          FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
             Cpool_litter_green_bg(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                &
                                               * Cpool_litter_green_bg(n,:))                                                  &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))
             Cpool_litter_wood_bg (n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                &
                                               * Cpool_litter_wood_bg (n,:))                                                  &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))
          END FORALL
       ELSE  !! with yasso 
          WHERE(is_vegetation(:,:))
             YCpool_acid_ag1       =   (cf_old * veg_fract_correction * YCpool_acid_ag1       - YCpool_acid_ag1_loss )         &
                                    / (cf_new * veg_fract_correction) 
             YCpool_water_ag1      =   (cf_old * veg_fract_correction * YCpool_water_ag1      - YCpool_water_ag1_loss )        &
                                    / (cf_new * veg_fract_correction) 
             YCpool_ethanol_ag1    =   (cf_old * veg_fract_correction * YCpool_ethanol_ag1    - YCpool_ethanol_ag1_loss )      &
                                    / (cf_new * veg_fract_correction) 
             YCpool_nonsoluble_ag1 =   (cf_old * veg_fract_correction * YCpool_nonsoluble_ag1 - YCpool_nonsoluble_ag1_loss )   &
                                    / (cf_new * veg_fract_correction) 
             YCpool_acid_ag2       =   (cf_old * veg_fract_correction * YCpool_acid_ag2       - YCpool_acid_ag2_loss )         &
                                    / (cf_new * veg_fract_correction) 
             YCpool_water_ag2      =   (cf_old * veg_fract_correction * YCpool_water_ag2      - YCpool_water_ag2_loss )        &
                                    / (cf_new * veg_fract_correction) 
             YCpool_ethanol_ag2      =   (cf_old * veg_fract_correction * YCpool_ethanol_ag2  - YCpool_ethanol_ag2_loss )      &
                                    / (cf_new * veg_fract_correction) 
             YCpool_nonsoluble_ag2 =   (cf_old * veg_fract_correction * YCpool_nonsoluble_ag2 - YCpool_nonsoluble_ag2_loss )   &
                                    / (cf_new * veg_fract_correction) 
          END WHERE
          FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
           ! new pool content =  inputs due to cover fraction change 
           !                   + inputs from vegetation pools (only belowground & acid soluble part)
           !                     scaled with new fraction
            YCpool_acid_bg1(n,i)    = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                    &
                                               * YCpool_acid_bg1(n,:))                                                        &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,1))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_water_bg1(n,i)   = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                    &
                                               * YCpool_water_bg1(n,:))                                                       &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,2))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_ethanol_bg1(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                    &
                                               * YCpool_ethanol_bg1(n,:))                                                     &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,3))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_nonsoluble_bg1(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                 &
                                               * YCpool_nonsoluble_bg1(n,:))                                                  &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,4))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_humus_1(n,i)   =  (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                     &
                                               * YCpool_humus_1(n,:))                                                         &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,5))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_acid_bg2(n,i)  = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                      &
                                               * YCpool_acid_bg2(n,:))                                                        &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,1))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_water_bg2(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                      &
                                               * YCpool_water_bg2(n,:))                                                       &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,2))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_ethanol_bg2(n,i)=(SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                      &
                                               * YCpool_ethanol_bg2(n,:))                                                     &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,3))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_nonsoluble_bg2(n,i)= (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                  &
                                               * YCpool_nonsoluble_bg2(n,:))                                                  &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,4))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_humus_2(n,i)   =  (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                     &
                                               * YCpool_humus_2(n,:))                                                         &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,5))/ (cf_new(n,i) * veg_fract_correction(n,i))
          END FORALL
       END IF

    ELSE ! lcc_scheme=1

       !! Carbon loss to atmosphere
       C_2_atmos(:) =   frac_green_2_atmos   * SUM(Cpool_green_loss  (:,:),DIM=2) &
                      + frac_wood_2_atmos    * SUM(Cpool_woods_loss  (:,:),DIM=2) &
                      + frac_reserve_2_atmos * SUM(Cpool_reserve_loss(:,:),DIM=2)

       !!nitrogen loss to atmosphere
       IF (nitrogenMode) THEN
           Nitrogen_2_atmos(:) =  frac_green_2_atmos  * SUM(Npool_green_loss(:,:), DIM=2) &
                                + frac_wood_2_atmos   * SUM(Npool_wood_loss(:,:),  DIM=2) &
                                + frac_mobile_2_atmos * SUM(Npool_mobile_loss(:,:),DIM=2)
       END IF

       !! Carbon gains of litter pools
       LitterGreen(:,:) =  (1._dp-frac_green_2_atmos  ) * Cpool_green_loss  (:,:) &
                         + (1._dp-frac_reserve_2_atmos) * Cpool_reserve_loss(:,:)
       LitterWood(:,:)  =  (1._dp-frac_wood_2_atmos  )  * Cpool_woods_loss  (:,:)

       IF (nitrogenMode) THEN
          LitterGreenN(:,:) = (1.0_dp - frac_green_2_atmos) * Npool_green_loss(:,:)
          SMINNgain(:,:) = (1.0_dp - frac_mobile_2_atmos) * Npool_mobile_loss(:,:) + &
                           (1.0_dp - frac_wood_2_atmos)*Cpool_woods_loss(:,:)        &
                            * (1.0_dp / cn_woods - 1.0_dp / cn_litter_wood)
          LitterWoodN(:,:) = (1.0_dp - frac_wood_2_atmos) * Cpool_woods_loss(:,:) / cn_litter_wood
       END IF
    
       !! Compute output arrays for litter
       C_2_litterGreenPools (:) = SUM(LitterGreen(:,:),DIM=2)
       C_2_litterWoodPool_ag(:) = SUM(LitterWood (:,:),DIM=2) *        frac_wood_aboveGround
       C_2_litterWoodPool_bg(:) = SUM(LitterWood (:,:),DIM=2) * (1._dp-frac_wood_aboveGround)

       !! New sizes of litter pools
       IF (.NOT. with_yasso) THEN
          FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
             Cpool_litter_green_ag(n,i) = &
                  ( frac_green_aboveGround * LitterGreen(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_litter_green_ag(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
             Cpool_litter_green_bg(n,i) = &
                  ( (1._dp-frac_green_aboveGround) * LitterGreen(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_litter_green_bg(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
             Cpool_litter_wood_ag(n,i) = &
                  ( frac_wood_aboveGround * LitterWood(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_litter_wood_ag(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
             Cpool_litter_wood_bg(n,i) = &
                  ( (1._dp-frac_wood_aboveGround) * LitterWood(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_litter_wood_bg(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
          END FORALL
          if (nitrogenMode) then
              FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
                 Npool_litter_green_ag(n,i) = &
                  ( frac_green_aboveGround * LitterGreenN(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_litter_green_ag(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
                 Npool_litter_green_bg(n,i) = &
                  ( (1._dp-frac_green_aboveGround) * LitterGreenN(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_litter_green_bg(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
                 Npool_litter_wood_ag(n,i) = &
                  ( frac_wood_aboveGround * LitterWoodN(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_litter_wood_ag(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))	    
                 Npool_litter_wood_bg(n,i) = &
                  ( (1._dp-frac_wood_aboveGround) * LitterWoodN(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_litter_wood_bg(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))	          
             end forall
            
              FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
                 SMINN_pool(n,i) = SMINNgain(n,i)                   &
                              + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*SMINN_pool(n,:) ) &
                                / (cf_new(n,i)*veg_fract_correction(n,i))
              END FORALL
          end if
       ELSE
          FORALL (n=1:nidx, i=1:ntiles, is_vegetation(n,i))
             YCpool_acid_ag1(n,i)       =                                                                     &  ! Aboveground
                  (  frac_green_aboveGround * LitterGreen(n,i) * LeafLitcoef(n,i,1)                          &  ! Input leaf litter
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:) *veg_fract_correction(n,:)*YCpool_acid_ag1(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_water_ag1(n,i)      =                                                                                 &
                  (  frac_green_aboveGround * LitterGreen(n,i) * LeafLitcoef(n,i,2)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_water_ag1(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_ethanol_ag1(n,i)    =                                                                                 &
                  (  frac_green_aboveGround * LitterGreen(n,i) * LeafLitcoef(n,i,3)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_ethanol_ag1(n,:) )    &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_nonsoluble_ag1(n,i) =                                                                                 &
                  (  frac_green_aboveGround * LitterGreen(n,i) * LeafLitcoef(n,i,4)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_nonsoluble_ag1(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_acid_bg1(n,i)       =                                                                                 &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,1)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_acid_bg1(n,:) )       &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_water_bg1(n,i)      =                                                                                 &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,2)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_water_bg1(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_ethanol_bg1(n,i)    =                                                                                 &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,3)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_ethanol_bg1(n,:) )    &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_nonsoluble_bg1(n,i) =                                                                                 &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,4)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_nonsoluble_bg1(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_humus_1(n,i) =                                                                                        &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,5)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_humus_1(n,:) )        &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))

            YCpool_acid_ag2(n,i)       =                                                                     &  ! Aboveground
                   ( frac_wood_aboveGround  * LitterWood(n,i)  * WoodLitcoef(n,i,1)                          &  ! Input wood litter
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:) *veg_fract_correction(n,:)*YCpool_acid_ag2(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_water_ag2(n,i)      =                                                                                 &
                   ( frac_wood_aboveGround  * LitterWood(n,i)  * WoodLitcoef(n,i,2)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_water_ag2(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_ethanol_ag2(n,i)    =                                                                                 &
                   ( frac_wood_aboveGround  * LitterWood(n,i)  * WoodLitcoef(n,i,3)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_ethanol_ag2(n,:) )    &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_nonsoluble_ag2(n,i) =                                                                                 &
                   ( frac_wood_aboveGround  * LitterWood(n,i)  * WoodLitcoef(n,i,4)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_nonsoluble_ag2(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_acid_bg2(n,i)       =                                                                                 &
                   ( (1._dp - frac_wood_aboveGround) * LitterWood(n,i)  * WoodLitcoef(n,i,1)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_acid_bg2(n,:) )       &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_water_bg2(n,i)      =                                                                                 &
                   ( (1._dp - frac_wood_aboveGround) * LitterWood(n,i)  * WoodLitcoef(n,i,2)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_water_bg2(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_ethanol_bg2(n,i)    =                                                                                 &
                   ( (1._dp - frac_wood_aboveGround) * LitterWood(n,i)  * WoodLitcoef(n,i,3)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_ethanol_bg2(n,:) )    &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_nonsoluble_bg2(n,i) =                                                                                 &
                   ( (1._dp - frac_wood_aboveGround) * LitterWood(n,i)  * WoodLitcoef(n,i,4)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_nonsoluble_bg2(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_humus_2(n,i) =                                                                                        &
                  (  (1._dp - frac_green_aboveGround) * LitterWood(n,i) * WoodLitcoef(n,i,5)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_humus_2(n,:) )        &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
          END FORALL
       END IF
    ENDIF !lcc_scheme == 1

    !! New size of slow soil pool - same treatment for all lcc_schemes
    IF (.NOT. with_yasso) THEN
       FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
          Cpool_slow(n,i) = SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_slow(n,:) ) &
                             / (cf_new(n,i)*veg_fract_correction(n,i))
       END FORALL
   ! DSG: YASSO - treatment of humus is not the same for llc_schemes    
    END IF

  END SUBROUTINE C_relocation_from_LUtransitions

  ! --- C_relocation_from_harvest() -----------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of harvest for the carbon pools.
  !
  SUBROUTINE C_relocation_from_harvest(nidx,ntiles,lcc_scheme,lctlib,surface, &
                         cf,veg_fract_correction, is_naturalVeg,              &
                         Cpool_green,Cpool_woods,Cpool_reserve,               &
                         harvest,                                             &
                         C2litter, C2atmos,                                   &
                         Cpool_litter_green_bg, Cpool_litter_wood_bg,         &
                         YCpool_acid_bg1, YCpool_water_bg1, YCpool_ethanol_bg1,  &
                         YCpool_nonsoluble_bg1, YCpool_humus_1,                  &
                         YCpool_acid_bg2, YCpool_water_bg2, YCpool_ethanol_bg2,  &
                         YCpool_nonsoluble_bg2, YCpool_humus_2,                  &
                         LeafLit_coef, WoodLit_coef,                          &
                         Cpool_paper, Cpool_construction,                     &
                         Carbon_2_paper, Carbon_2_construction,               &
                         Carbon_paper_2_atmos, Carbon_construction_2_atmos,   &
                         frac_harvest_2_atmos,                                &
                         Npool_green, Npool_woods,Npool_mobile,               &
                         Npool_litter_green_bg, Npool_litter_wood_bg,         &
                         SMINN_pool,                                          &
                         N2atmos, frac_mobile_2_atmos)

    USE mo_jsbach_lctlib,   ONLY: lctlib_type
    USE mo_land_surface,    ONLY: land_surface_type
    USE mo_cbal_parameters, ONLY: tau_paper, tau_construction
    USE mo_exception,       ONLY: finish

    integer,  intent(in)                :: nidx                       !! Number of gridpoints in vectors
    integer,  intent(in)                :: ntiles                     !! Number of tiles
    integer,  intent(in)                :: lcc_scheme
    type(lctlib_type), intent(in)       :: lctlib
    type(land_surface_type), intent(in) :: surface
    real(dp), intent(in)                :: cf(:,:)                    !! Cover fractions before landcover change or dynveg [] 
    real(dp), intent(in)                :: veg_fract_correction(:,:)  !! Correction factor for cover fractions 1-exp(-LAI_max/2)
                                                                      !!   (accounts for sparseness of vegetation)
    logical,  intent(in)                :: is_naturalVeg(:,:)         !! Indicates tiles with natural (non-agricultural) vegetation
    real(dp), intent(inout)             :: Cpool_green(:,:)           !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp), intent(inout)             :: Cpool_woods(:,:)           !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp), intent(inout)             :: Cpool_reserve(:,:)         !! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    real(dp), intent(in)                :: harvest(:)                 !! Carbon harvested [mol(C)/m^2(corr. veget. area)/timestep]
    real(dp), intent(inout), optional   :: Cpool_litter_green_bg(:,:) !! Value of below ground green litter C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: Cpool_litter_wood_bg(:,:)  !! Value of below ground wood litter C pool [mol/m^2(canopy)]
                                                                      !! (e.g. by  using it as fuel wood within days)
                                                                                       !! <---DSG: what?!
    real(dp), intent(inout), optional   :: YCpool_acid_bg1(:,:)        !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_water_bg1(:,:)       !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_ethanol_bg1(:,:)     !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_nonsoluble_bg1(:,:)  !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_humus_1(:,:)         !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_acid_bg2(:,:)        !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_water_bg2(:,:)       !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_ethanol_bg2(:,:)     !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_nonsoluble_bg2(:,:)  !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_humus_2(:,:)         !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(in),    optional   :: LeafLit_coef(:,:,:)        !! Factor to seperate non woody litte into EWAN pools [ ]
    real(dp), intent(in),    optional   :: WoodLit_coef(:,:,:)        !! Factor to seperate woody litte into EWAN pools     [ ]
    real(dp), intent(out)               :: C2litter(:)                !! Amount of harvested carbon relocated to the below ground 
                                                                      !!     green litter pools in this timestep 
                                                                      !!     [mol(C)/m^2(vegetated area)/timestep]
    real(dp), intent(out)              :: C2atmos(:)                  !! Amount of harvested C directly emitted to atmosphere in 
                                                                      !!     this timestep [mol(C)/m^2(vegetated area)/timestep]
    real(dp), intent(in),    optional  :: frac_harvest_2_atmos        !! Fraction of harvested carbon immediately emitted to atmos. 
    real(dp), intent(inout), optional  :: Cpool_paper(:)              !! Anthropogenic pools & fluxes : ...
    real(dp), intent(inout), optional  :: Cpool_construction(:)
    real(dp), intent(out),   optional  :: Carbon_2_paper(:)
    real(dp), intent(out),   optional  :: Carbon_2_construction(:)
    real(dp), intent(out),   optional  :: Carbon_paper_2_atmos(:)
    real(dp), intent(out),   optional  :: Carbon_construction_2_atmos(:)

    real(dp),intent(inout),optional :: Npool_green(:,:)
    real(dp),intent(inout),optional :: Npool_woods(:,:)
    real(dp),intent(inout),optional :: Npool_mobile(:,:)
    real(dp),intent(inout),optional :: Npool_litter_green_bg(:,:)
    real(dp),intent(inout),optional :: Npool_litter_wood_bg(:,:)
    real(dp),intent(inout),optional :: SMINN_pool(:,:)
    real(dp),intent(out),optional   :: N2atmos(:)  
    real(dp),intent(in),optional    :: frac_mobile_2_atmos

    !! Local variables
    real(dp) :: Cpool_G_ag_corr(1:nidx,1:ntiles)     !! Above ground carbon in green pool (per corrected vegetated area)
    real(dp) :: Cpool_R_ag_corr(1:nidx,1:ntiles)     !! Above ground carbon in reserve pool (per corrected vegetated area)
    real(dp) :: Cpool_W_ag_corr(1:nidx,1:ntiles)     !! Above ground carbon in wood pool (per corrected vegetated area)

    real(dp) :: Cpool_nat_ag_corr(1:nidx)     !! Total above ground carbon from natural vegetation (per corrected vegetated area)
    real(dp) :: Harvest_G(1:nidx,1:ntiles)  !! Harvest from green pool
    real(dp) :: Harvest_R(1:nidx,1:ntiles)  !! Harvest from reserve pool
    real(dp) :: Harvest_W(1:nidx,1:ntiles)  !! Harvest from wood pools

    real(dp) :: Harvest_G_rel(1:nidx,1:ntiles)  !! fraction harvested from green pool (Harvest_G/Cpool_green)
    real(dp) :: Npool_mobile_loss(1:nidx,1:ntiles) !!amount of mobile N pools harvested

    real(dp) :: scale_factor(1:nidx,1:ntiles)

    integer  :: n, k, ct
    real(dp) :: Tmp
    
    logical  :: with_yasso 
    logical :: nitrogenMode

    !! --- GO! ---------------------------------------

    !! Check for nitrogen mode
    nitrogenMode=.false.
    if(present(Npool_green)           .or. present(Npool_woods) .or. &
       present(Npool_mobile)          .or. &
       present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. &
       present(SMINN_pool) ) then
       nitrogenMode=.true.
    end if

    !! Check if yasso or cbalance litter and soil pools should be handled
    with_yasso  = .FALSE.  
    IF (PRESENT(YCpool_acid_bg1)            .OR. &
        PRESENT(YCpool_water_bg1)           .OR. &
        PRESENT(YCpool_ethanol_bg1)         .OR. &
        PRESENT(YCpool_nonsoluble_bg1)      .OR. &
        PRESENT(YCpool_humus_2)             .OR. &
         PRESENT(YCpool_acid_bg2)            .OR. &
         PRESENT(YCpool_water_bg2)           .OR. &
         PRESENT(YCpool_ethanol_bg2)         .OR. &
         PRESENT(YCpool_nonsoluble_bg2)      .OR. &
        PRESENT(YCpool_humus_2)         .OR. PRESENT(LeafLit_coef)                .OR. &
         PRESENT(WoodLit_coef)) THEN
        with_yasso=.TRUE.
        IF (.NOT. (PRESENT(YCpool_acid_bg1)        .AND. &
                   PRESENT(YCpool_water_bg1)       .AND. &
                   PRESENT(YCpool_ethanol_bg1)     .AND. &
                   PRESENT(YCpool_nonsoluble_bg1)  .AND. &
                   PRESENT(YCpool_humus_1)         .AND. &
                   PRESENT(YCpool_acid_bg2)        .AND. &
                   PRESENT(YCpool_water_bg2)       .AND. &
                   PRESENT(YCpool_ethanol_bg2)     .AND. &
                   PRESENT(YCpool_nonsoluble_bg2)  .AND. &
                   PRESENT(YCpool_humus_2)          .AND. PRESENT(LeafLit_coef)           .AND. &
                   PRESENT(WoodLit_coef) ) &
            ) THEN
            CALL finish('relocate_carbon_desert()','at least one variable missing to handle yasso pools')
        END IF
    END IF

    !! Go
    !!dk: new function in cosmos-landveg with anthro pools depending on lcc_scheme
    !!dk: but N runs should be comparable with cmip5 runs, therefore no anthro pools considered yet in harvest

    !! Above ground carbon in vegetation pools
    Cpool_G_ag_corr(1:nidx,1:ntiles) = frac_green_aboveGround * veg_fract_correction(:,:) * cf(:,:) * Cpool_green  (:,:)
    Cpool_R_ag_corr(1:nidx,1:ntiles) = frac_green_aboveGround * veg_fract_correction(:,:) * cf(:,:) * Cpool_reserve(:,:)
    Cpool_W_ag_corr(1:nidx,1:ntiles) = frac_wood_aboveGround  * veg_fract_correction(:,:) * cf(:,:) * Cpool_woods  (:,:)

    !! Total above carbon in pools of natural vegetation of a grid cell
    Cpool_nat_ag_corr(1:nidx) = 0.0_dp
    Cpool_nat_ag_corr(1:nidx) = SUM(   Cpool_G_ag_corr(1:nidx,1:ntiles) &
                                   + Cpool_R_ag_corr(1:nidx,1:ntiles) &
                                   + Cpool_W_ag_corr(1:nidx,1:ntiles) &
                                            , mask=is_naturalVeg(1:nidx,1:ntiles), DIM=2)

    !! Harvest from vegetation pools assuming harvest proportional to above ground carbon
    Harvest_G(1:nidx,1:ntiles) = 0.0_dp
    Harvest_R(1:nidx,1:ntiles) = 0.0_dp
    Harvest_W(1:nidx,1:ntiles) = 0.0_dp
    Harvest_G_rel(1:nidx,1:ntiles) = 0.0_dp

    FORALL(n=1:nidx,k=1:ntiles, is_naturalVeg(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
       Harvest_G(n,k) = MIN(harvest(n),Cpool_nat_ag_corr(n))*Cpool_G_ag_corr(n,k)/Cpool_nat_ag_corr(n)
       Harvest_R(n,k) = MIN(harvest(n),Cpool_nat_ag_corr(n))*Cpool_R_ag_corr(n,k)/Cpool_nat_ag_corr(n)
       Harvest_W(n,k) = MIN(harvest(n),Cpool_nat_ag_corr(n))*Cpool_W_ag_corr(n,k)/Cpool_nat_ag_corr(n) 
    END FORALL

    WHERE (veg_fract_correction(:,:) > 0._dp)
       scale_factor(:,:) = 1._dp / (veg_fract_correction(:,:) * MAX(cf(:,:),fract_small))
    ELSEWHERE
       scale_factor(:,:) = 0._dp
    ENDWHERE

    FORALL(n=1:nidx,k=1:ntiles, is_naturalVeg(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp) &
                                                   .AND. Cpool_green(n,k) > EPSILON(1.0_dp))
       Harvest_G_rel(n,k) = Harvest_G(n,k) * scale_factor(n,k) / Cpool_green(n,k)
    END FORALL

    !! Depletion of living plant pools (note: veg_fract_correction(:,:)=1 in glacier boxes) 
    FORALL(n=1:nidx,k=1:ntiles, is_naturalVeg(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
       Cpool_green  (n,k) = Cpool_green  (n,k) - Harvest_G(n,k) * scale_factor(n,k)
       Cpool_reserve(n,k) = Cpool_reserve(n,k) - Harvest_R(n,k) * scale_factor(n,k)
       Cpool_woods  (n,k) = Cpool_woods  (n,k) - Harvest_W(n,k) * scale_factor(n,k)
    END FORALL
    
   if (nitrogenMode) then
     FORALL(n=1:nidx,k=1:ntiles, is_naturalVeg(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
       Npool_green(n,k)   = Npool_green(n,k)     - (Harvest_G(n,k)/cn_green) * scale_factor(n,k) 
       Npool_woods(n,k)   = Npool_woods(n,k)     - (Harvest_W(n,k)/cn_woods) * scale_factor(n,k)
       Npool_mobile_loss(n,k)  = Npool_mobile(n,k) * Harvest_G_rel(n,k)
                  !dk: relative decreased in mobileN is equal to relativ decrease in Cpool_green 
     END FORALL	    
   end if

    IF (lcc_scheme==1) THEN

       !! Increase in below ground litter pools
       IF (.NOT. with_yasso) THEN
          FORALL(n=1:nidx,k=1:ntiles, is_naturalVeg(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
             Cpool_litter_green_bg(n,k) = &
             Cpool_litter_green_bg(n,k) + (1._dp-frac_harvest_2_atmos) * scale_factor(n,k) * (Harvest_G(n,k)+Harvest_R(n,k)) 
             Cpool_litter_wood_bg (n,k) = &
             Cpool_litter_wood_bg (n,k) + (1._dp-frac_harvest_2_atmos) * scale_factor(n,k) *  Harvest_W(n,k)
          END FORALL

          if (nitrogenMode) then
            FORALL(n=1:nidx,k=1:ntiles, is_naturalVeg(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
                Npool_litter_wood_bg(n,k) = Npool_litter_wood_bg(n,k)  + &
                    (1._dp - frac_harvest_2_atmos) * scale_factor(n,k) * (Harvest_W(n,k)/cn_litter_wood)
                Npool_litter_green_bg(n,k) = Npool_litter_green_bg(n,k) + &
                    (1._dp - frac_harvest_2_atmos) * scale_factor(n,k) * (Harvest_G(n,k)/cn_green)
                SMINN_pool(n,k) = SMINN_pool(n,k) + (1._dp - frac_harvest_2_atmos) * scale_factor(n,k) * &
                    (Harvest_W(n,k) * (1._dp / cn_woods - 1._dp / cn_litter_wood))
                SMINN_pool(n,k) = SMINN_pool(n,k) + &
                    (1._dp - frac_mobile_2_atmos) * Npool_mobile_loss(n,k)
                Npool_mobile(n,k) = Npool_mobile(n,k) - Npool_mobile_loss(n,k)
            END FORALL
          end if

       ELSE
          FORALL(n=1:nidx,k=1:ntiles, is_naturalVeg(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
             YCpool_acid_bg1      (n,k) =                                                                  &
             YCpool_acid_bg1      (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                       * ( LeafLit_coef(n,k,1)    * (Harvest_G(n,k)+Harvest_R(n,k))  )       ! Non woody litter input
             YCpool_water_bg1     (n,k) =                                                                  &
             YCpool_water_bg1     (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                       * ( LeafLit_coef(n,k,2)     * (Harvest_G(n,k)+Harvest_R(n,k)) )       ! Non woody litter input
             YCpool_ethanol_bg1   (n,k) =                                                                  &
             YCpool_ethanol_bg1   (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                       * ( LeafLit_coef(n,k,3)     * (Harvest_G(n,k)+Harvest_R(n,k)) )       ! Non woody litter input
             YCpool_nonsoluble_bg1(n,k) =                                                                  &
             YCpool_nonsoluble_bg1(n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                       * ( LeafLit_coef(n,k,4)     * (Harvest_G(n,k)+Harvest_R(n,k)) )       ! Non woody litter input
             YCpool_humus_1        (n,k) =                                                                 &
             YCpool_humus_1        (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)               &
                                       * ( LeafLit_coef(n,k,5)     * (Harvest_G(n,k)+Harvest_R(n,k)) )       ! Non woody litter input


             YCpool_acid_bg2      (n,k) =                                                                  &
             YCpool_acid_bg2      (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                           * WoodLit_coef(n,k,1)  * Harvest_W(n,k)                          ! Woody litter input
             YCpool_water_bg2     (n,k) =                                                                  &
             YCpool_water_bg2     (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                           * WoodLit_coef(n,k,2)   * Harvest_W(n,k)                         ! Woody litter input
             YCpool_ethanol_bg2   (n,k) =                                                                  &
             YCpool_ethanol_bg2   (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                           * WoodLit_coef(n,k,3)   * Harvest_W(n,k)                         ! Woody litter input
             YCpool_nonsoluble_bg2(n,k) =                                                                  &
             YCpool_nonsoluble_bg2(n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                           * WoodLit_coef(n,k,4)   * Harvest_W(n,k)                         ! Woody litter input
             YCpool_humus_2        (n,k) =                                                                  &
             YCpool_humus_2        (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                           * WoodLit_coef(n,k,5)   * Harvest_W(n,k)                         ! Woody litter input
          END FORALL
       END IF
       !! Direct emissions to atmosphere
       C2atmos(1:nidx) = frac_harvest_2_atmos * &
            SUM(Harvest_G(1:nidx,1:ntiles) + Harvest_R(1:nidx,1:ntiles) + Harvest_W(1:nidx,1:ntiles), DIM=2)

       !! Harvest into litter pools
       C2litter(1:nidx) = (1._dp-frac_harvest_2_atmos) * &
                             SUM(Harvest_G(1:nidx,1:ntiles) + Harvest_R(1:nidx,1:ntiles) + Harvest_W(1:nidx,1:ntiles), DIM=2)

       IF (nitrogenMode) THEN
          N2atmos(1:nidx) = frac_harvest_2_atmos * &
             SUM(Harvest_G(1:nidx,1:ntiles) / cn_green + Harvest_W(1:nidx,1:ntiles) / cn_woods, DIM=2)
       END IF
!dk: harvested biomass only from above ground green pool, while green pool is distributed to _ag and _bg with constant ratio
!dk: might have implications for other calculations, as _bg pool is diminished and _ag pool is larger by applying constant ratio
  
    ELSE ! lcc_scheme==2

       !! Relocate Carbon into anthropogenic pools (paper and construction)
       Carbon_2_paper       (:) = 0._dp
       Carbon_2_construction(:) = 0._dp
       C2litter             (:) = 0._dp
       DO k = 1, ntiles
          DO n=1,nidx
             IF ( is_naturalVeg(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp)) THEN
                ct  = surface%cover_type(n,k)
                IF (lctlib%frac_wood_2_paper(ct) + lctlib%frac_wood_2_construction(ct) >  EPSILON(1.0_dp) ) THEN

                   Tmp =   (Harvest_G(n,k) + Harvest_R(n,k) + Harvest_W(n,k))                   &
                         / (lctlib%frac_wood_2_paper(ct) + lctlib%frac_wood_2_construction(ct))
                
                   Carbon_2_paper       (n) = Carbon_2_paper       (n) + Tmp * lctlib%frac_wood_2_paper       (ct)
                   Carbon_2_construction(n) = Carbon_2_construction(n) + Tmp * lctlib%frac_wood_2_construction(ct)
                ELSE
                   IF (.NOT. with_yasso) THEN
                      Cpool_litter_green_bg(n,k) = Cpool_litter_green_bg(n,k) + scale_factor(n,k) * (Harvest_G(n,k)+Harvest_R(n,k))
                      Cpool_litter_wood_bg (n,k) = Cpool_litter_wood_bg (n,k) + scale_factor(n,k) *  Harvest_W(n,k)
                   ELSE
                      YCpool_acid_bg1      (n,k) = YCpool_acid_bg1      (n,k) +  LeafLit_coef(n,k,1)          &     ! DSG: not tested
                                                  * scale_factor(n,k) * (Harvest_G(n,k)+Harvest_R(n,k))
                      YCpool_water_bg1     (n,k) = YCpool_water_bg1     (n,k) +  LeafLit_coef(n,k,2)          &
                                                  * scale_factor(n,k) * (Harvest_G(n,k)+Harvest_R(n,k))
                      YCpool_ethanol_bg1   (n,k) = YCpool_ethanol_bg1   (n,k) +  LeafLit_coef(n,k,3)          &
                                                  * scale_factor(n,k) * (Harvest_G(n,k)+Harvest_R(n,k))
                      YCpool_nonsoluble_bg1(n,k) = YCpool_nonsoluble_bg1(n,k) +  LeafLit_coef(n,k,4)          &
                                                  * scale_factor(n,k) * (Harvest_G(n,k)+Harvest_R(n,k))
                      YCpool_humus_1        (n,k) = YCpool_humus_1      (n,k) +  LeafLit_coef(n,k,5)          &
                                                  * scale_factor(n,k) * (Harvest_G(n,k)+Harvest_R(n,k))
                      YCpool_acid_bg2      (n,k) = YCpool_acid_bg2      (n,k) +  WoodLit_coef(n,k,1)          &
                                                  * scale_factor(n,k) * Harvest_W(n,k)
                      YCpool_water_bg2     (n,k) = YCpool_water_bg2     (n,k) +  WoodLit_coef(n,k,2)          &
                                                  * scale_factor(n,k) * Harvest_W(n,k)
                      YCpool_ethanol_bg2   (n,k) = YCpool_ethanol_bg2   (n,k) +  WoodLit_coef(n,k,3)          &
                                                  * scale_factor(n,k) * Harvest_W(n,k)
                      YCpool_nonsoluble_bg2(n,k) = YCpool_nonsoluble_bg2(n,k) +  WoodLit_coef(n,k,4)          &
                                                  * scale_factor(n,k) * Harvest_W(n,k)
                      YCpool_humus_2        (n,k) = YCpool_humus_2      (n,k) +  WoodLit_coef(n,k,5)          &
                                                  * scale_factor(n,k) * Harvest_W(n,k)
                   END IF
                   C2litter             (n  ) = C2litter             (n  )                          &
                                              * (Harvest_G(n,k) + Harvest_R(n,k) + Harvest_W(n,k))
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       !! Carbon from anthropogenic pools to atmosphere
       Carbon_paper_2_atmos       (:) = Cpool_paper       (:) / tau_paper
       Carbon_construction_2_atmos(:) = Cpool_construction(:) / tau_construction

       !! New sizes of anthropogenic pools
       Cpool_paper       (:) = Cpool_paper       (:) + Carbon_2_paper       (:) - Carbon_paper_2_atmos       (:)
       Cpool_construction(:) = Cpool_construction(:) + Carbon_2_construction(:) - Carbon_construction_2_atmos(:)

       !! Total CO2 emission from harvest
       C2atmos(1:nidx) = Carbon_paper_2_atmos(1:nidx) + Carbon_construction_2_atmos(1:nidx)

    ENDIF 

  END SUBROUTINE C_relocation_from_harvest

! Relocation of carbon to anthropogenic pools from landcover maps or landuse transitions
! and update of anthropogenic pools
 
  SUBROUTINE C_loss_and_update_anthro_pools(nidx,ntiles,lcc_scheme,lctlib,surface,        &
                                            Cpool_green,Cpool_wood,Cpool_reserve,         &
                                            Cpool_litter_green_ag,Cpool_litter_wood_ag,   &
                                            YCpool_acid_ag1, YCpool_water_ag1,            &
                                            YCpool_ethanol_ag1, YCpool_nonsoluble_ag1,    &
                                            YCpool_acid_ag2, YCpool_water_ag2,            &
                                            YCpool_ethanol_ag2, YCpool_nonsoluble_ag2,    &
                                            Cpool_onSite,Cpool_paper,Cpool_construction,  &
                                            scale_fac,                                    &
                                            C_2_onSite,C_2_paper,C_2_construction,        &
                                            C_onSite_2_atmos,                             &
                                            C_paper_2_atmos,C_construction_2_atmos,       &
                                            C_2_atmos,                                    &
                                            C_green_loss,                                 & 
                                            C_wood_loss,                                  &
                                            C_reserve_loss,                               &
                                            C_litter_green_ag_loss,                       &
                                            C_litter_wood_ag_loss,                        &
                                            !for N option
                                            Npool_green,                                  &
                                            Npool_woods,                                  &
                                            Npool_mobile,                                 &
                                            Npool_green_loss,                             &
                                            Npool_wood_loss,                              &
                                            Npool_mobile_loss,                            &
                                            !yasso
                                            C_acid_ag1_loss,                              &
                                            C_water_ag1_loss,                             &
                                            C_ethanol_ag1_loss,                           &
                                            C_nonsoluble_ag1_loss,                        &
                                            C_acid_ag2_loss,                              &
                                            C_water_ag2_loss,                             &
                                            C_ethanol_ag2_loss,                           &
                                            C_nonsoluble_ag2_loss)

  USE mo_cbal_parameters, ONLY : tau_onSite, tau_paper, tau_construction
  USE mo_jsbach_lctlib,   ONLY : lctlib_type
  USE mo_land_surface,    ONLY : land_surface_type
  USE mo_exception,       ONLY : finish

  INTEGER,                  INTENT(in)    :: nidx, ntiles     ! Grid size
  INTEGER,                  INTENT(in)    :: lcc_scheme       ! landcover change scheme number
  TYPE (lctlib_type),       INTENT(in)    :: lctlib           ! PFT specific constants
  TYPE (land_surface_type), INTENT(in)    :: surface          ! Surface cover types
  REAL(dp),                 INTENT(inout) ::               &
    Cpool_green(:,:), Cpool_reserve(:,:), Cpool_wood(:,:)     ! Carbon pools of living material
  REAL(dp), OPTIONAL,       INTENT(inout) ::               &  ! CBALANCE
    Cpool_litter_green_ag(:,:), Cpool_litter_wood_ag(:,:)     !   Above ground carbon litter pools
  REAL(dp), OPTIONAL,       INTENT(inout) ::               &  ! YASSO
    YCpool_acid_ag1(:,:), YCpool_water_ag1(:,:),             &  !   Above ground carbon litter pools
    YCpool_ethanol_ag1(:,:), YCpool_nonsoluble_ag1(:,:),     &  !   Above ground carbon litter pools
    YCpool_acid_ag2(:,:), YCpool_water_ag2(:,:),             &  !   Above ground carbon litter pools
    YCpool_ethanol_ag2(:,:), YCpool_nonsoluble_ag2(:,:)         !   Above ground carbon litter pools
  REAL(dp), OPTIONAL,       INTENT(inout) ::               &  ! Optional since optional to C_relocation_from_LUtransitions
    Cpool_onSite(:), Cpool_paper(:), Cpool_construction(:)    ! Anthropogenic carbon pools
  REAL(dp),                 INTENT(in)    ::               &
    scale_fac(:,:)                                            ! Scaling factor for C fluxes relating to lcc
  REAL(dp), OPTIONAL,       INTENT(out)   ::               &  ! Optional since optional to C_relocation_from_LUtransitions
    C_2_onSite(:), C_2_paper(:), C_2_construction(:),      &  ! Carbon floating into anthropogenic pools
    C_onSite_2_atmos(:),C_paper_2_atmos(:),C_Construction_2_atmos(:), & ! Carbon from anthropogenic pools to atmosphere
    C_2_atmos(:)                                              ! Total carbon floating from anthropogenic pools to atmosphere
! These variables are needed only for lcc specified by transition matrices
  REAL(dp), OPTIONAL,TARGET,INTENT(out)   ::               &
    C_green_loss(:,:),                                     &  ! Carbon lost from green pool
    C_wood_loss(:,:),                                      &  ! Carbon lost from wood pool
    C_reserve_loss(:,:),                                   &  ! Carbon lost from reserve pool
    C_litter_green_ag_loss(:,:),                           &  ! Carbon lost from green litter pool above ground
    C_litter_wood_ag_loss(:,:),                            &  ! Carbon lost from green litter pool above ground
    ! Yasso   
    C_acid_ag1_loss(:,:),                                   &  ! Carbon lost from acid soluble pool above ground
    C_water_ag1_loss(:,:),                                  &  ! Carbon lost from water soluble pool above ground
    C_ethanol_ag1_loss(:,:),                                &  ! Carbon lost from ethanol soluble pool above ground
    C_nonsoluble_ag1_loss(:,:),                             &  ! Carbon lost from non soluble pool above ground
    C_acid_ag2_loss(:,:),                                   &  ! Carbon lost from acid soluble pool above ground
    C_water_ag2_loss(:,:),                                  &  ! Carbon lost from water soluble pool above ground
    C_ethanol_ag2_loss(:,:),                                &  ! Carbon lost from ethanol soluble pool above ground
    C_nonsoluble_ag2_loss(:,:)                                 ! Carbon lost from non soluble pool above ground
! These variables are needed for land use transitions with nitrogen
  REAL(dp), OPTIONAL,       INTENT(in)    :: &
    Npool_green(:,:), Npool_woods(:,:), Npool_mobile(:,:)      ! Nitrogen pools of living material
  REAL(dp), OPTIONAL,TARGET,INTENT(out)   ::  &
    Npool_green_loss(:,:),                                 &   ! Nitrogen lost from green pool
    Npool_wood_loss(:,:),                                  &   ! Nitrogen lost from wood pool
    Npool_mobile_loss(:,:)                                     ! Nitrogen lost from mobile plant pool
! locals

  INTEGER  :: i, itile, ct
  LOGICAL  :: with_yasso  
  REAL(dp), POINTER :: CGreenLoss(:,:), CWoodLoss(:,:), CReserveLoss(:,:), &
                       NGreenLoss(:,:), NWoodLoss(:,:), NMobileLoss(:,:),  &
                       CLitterGreenAGLoss(:,:), CLitterWoodAGLoss(:,:)
  REAL(dp), POINTER :: CAcidAG1Loss(:,:),                                  &
                       CWaterAG1Loss(:,:),                                 &
                       CEthanolAG1Loss(:,:),                               &
                       CNonsolubleAG1Loss(:,:),                            &
                       CAcidAG2Loss(:,:),                                  &
                       CWaterAG2Loss(:,:),                                 &
                       CEthanolAG2Loss(:,:),                               &
                       CNonsolubleAG2Loss(:,:)

    !! Check if yasso or cbalance litter and soil pools should be handled
    with_yasso  = .FALSE.  
    IF (PRESENT(YCpool_acid_ag1)        .OR.  &
        PRESENT(YCpool_water_ag1)       .OR.  &
        PRESENT(YCpool_ethanol_ag1)     .OR.  &
        PRESENT(YCpool_nonsoluble_ag1)  .OR.  &
        PRESENT(YCpool_acid_ag2)        .OR.  &
        PRESENT(YCpool_water_ag2)       .OR.  &
        PRESENT(YCpool_ethanol_ag2)     .OR.  &
        PRESENT(YCpool_nonsoluble_ag2)        &
        ) THEN
        with_yasso=.TRUE.
        IF (.NOT. (PRESENT(YCpool_acid_ag1)        .AND. &
                   PRESENT(YCpool_water_ag1)       .AND. &
                   PRESENT(YCpool_ethanol_ag1)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag1)  .AND. &
                   PRESENT(YCpool_acid_ag2)        .AND. &
                   PRESENT(YCpool_water_ag2)       .AND. &
                   PRESENT(YCpool_ethanol_ag2)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag2)        &
                   ) &
            ) CALL finish('C_loss_and_update_anthro_pools()','at least one variable missing to handle yasso pools')
    END IF

    IF (PRESENT(C_green_loss)) THEN
       CGreenLoss         => C_green_loss
       CWoodLoss          => C_wood_loss
       CReserveLoss       => C_reserve_loss
    ELSE
       ALLOCATE(CGreenLoss(nidx,ntiles),CWoodLoss(nidx,ntiles),CReserveLoss(nidx,ntiles))
    ENDIF

    IF (PRESENT(Npool_green_loss)) THEN
      NGreenLoss         => Npool_green_loss
      NWoodLoss          => Npool_wood_loss
      NMobileLoss        => Npool_mobile_loss
    ELSE
      ALLOCATE(NGreenLoss(nidx,ntiles),NWoodLoss(nidx,ntiles),NMobileLoss(nidx,ntiles))
    ENDIF

    IF (.NOT. with_yasso) THEN
       IF (PRESENT(C_litter_green_ag_loss)) THEN 
          CLitterGreenAGLoss => C_litter_green_ag_loss
          CLitterWoodAGLoss  => C_litter_wood_ag_loss
       ELSE
          ALLOCATE(CLitterGreenAGLoss(nidx,ntiles), CLitterWoodAGLoss(nidx,ntiles))
       END IF
    ELSE
       IF (PRESENT(C_acid_ag1_loss)) THEN
          CAcidAG1Loss        => C_acid_ag1_loss
          CWaterAG1Loss       => C_water_ag1_loss
          CEthanolAG1Loss     => C_ethanol_ag1_loss
          CNonsolubleAG1Loss  => C_nonsoluble_ag1_loss
          CAcidAG2Loss        => C_acid_ag2_loss
          CWaterAG2Loss       => C_water_ag2_loss
          CEthanolAG2Loss     => C_ethanol_ag2_loss
          CNonsolubleAG2Loss  => C_nonsoluble_ag2_loss
       ELSE
          ALLOCATE(CAcidAG1Loss(nidx,ntiles),         &
                   CWaterAG1Loss(nidx,ntiles),        &
                   CEthanolAG1Loss(nidx,ntiles),      &
                   CNonsolubleAG1Loss(nidx,ntiles)  )
          ALLOCATE(CAcidAG2Loss(nidx,ntiles),         &
                   CWaterAG2Loss(nidx,ntiles),        &
                   CEthanolAG2Loss(nidx,ntiles),      &
                   CNonsolubleAG2Loss(nidx,ntiles)  )
       END IF
    ENDIF

    IF (PRESENT(C_green_loss) .OR. (lcc_scheme==2)) THEN
       CGreenLoss        (:,:) = scale_fac(:,:) * Cpool_green          (:,:) 
       CWoodLoss         (:,:) = scale_fac(:,:) * Cpool_wood           (:,:)
       CReserveLoss      (:,:) = scale_fac(:,:) * Cpool_reserve        (:,:)
       IF (.NOT. with_yasso) THEN
          CLitterGreenAGLoss(:,:) = scale_fac(:,:) * Cpool_litter_green_ag(:,:)  
          CLitterWoodAGLoss (:,:) = scale_fac(:,:) * Cpool_litter_wood_ag (:,:)  
       ELSE
          CAcidAG1Loss       (:,:) = scale_fac(:,:) * YCpool_acid_ag1       (:,:) 
          CWaterAG1Loss      (:,:) = scale_fac(:,:) * YCpool_water_ag1      (:,:)
          CEthanolAG1Loss    (:,:) = scale_fac(:,:) * YCpool_ethanol_ag1    (:,:)
          CNonsolubleAG1Loss (:,:) = scale_fac(:,:) * YCpool_nonsoluble_ag1 (:,:)
          CAcidAG2Loss       (:,:) = scale_fac(:,:) * YCpool_acid_ag2       (:,:) 
          CWaterAG2Loss      (:,:) = scale_fac(:,:) * YCpool_water_ag2      (:,:)
          CEthanolAG2Loss    (:,:) = scale_fac(:,:) * YCpool_ethanol_ag2    (:,:)
          CNonsolubleAG2Loss (:,:) = scale_fac(:,:) * YCpool_nonsoluble_ag2 (:,:)
       END IF
    ENDIF


    IF (PRESENT(Npool_green_loss)) THEN
      NGreenLoss  (:,:) = scale_fac(:,:) * Npool_green(:,:)
      NWoodLoss   (:,:) = scale_fac(:,:) * Npool_woods(:,:)
      NMobileLoss (:,:) = scale_fac(:,:) * Npool_mobile(:,:)
    ENDIF

    IF (lcc_scheme==2) THEN

      C_2_onSite      (:) = 0._dp
      C_2_paper       (:) = 0._dp
      C_2_construction(:) = 0._dp

      ! Add aboveground wood loss to anthropogenic pools
      DO itile=1,ntiles
        DO i=1,nidx
          ct=surface%cover_type(i,itile)
          C_2_construction(i) = C_2_construction(i) + CWoodLoss(i,itile)*frac_wood_aboveGround*lctlib%frac_wood_2_construction(ct)
          C_2_paper       (i) = C_2_paper       (i) + CWoodLoss(i,itile)*frac_wood_aboveGround*lctlib%frac_wood_2_paper(ct)
          C_2_onSite(i)       = C_2_onSite      (i) + CWoodLoss(i,itile)*frac_wood_aboveGround*lctlib%frac_wood_2_onSite(ct)  &
                                                    + (CGreenLoss(i,itile) + CReserveLoss(i,itile))*frac_green_aboveGround
        ENDDO
      ENDDO
      ! Add litter inputs
      IF (.NOT. with_yasso) THEN
        DO itile=1,ntiles
          DO i=1,nidx
            ct=surface%cover_type(i,itile)
            C_2_onSite(i) = C_2_onSite(i)                                                                       &  
                               +  CLitterGreenAGLoss(i,itile)                                                   &
                               +  CLitterWoodAGLoss (i,itile)                                                   
          ENDDO
        ENDDO
      ELSE  ! with_yasso
        DO itile=1,ntiles
          DO i=1,nidx
            ct=surface%cover_type(i,itile)
            C_2_onSite(i) = C_2_onSite(i)                                                                       &  
                               +   CAcidAG1Loss(i,itile)                                                        & 
                               +   CAcidAG2Loss(i,itile)                                                        & 
                               +   CWaterAG1Loss(i,itile)                                                       & 
                               +   CWaterAG2Loss(i,itile)                                                       & 
                               +   CEthanolAG1Loss(i,itile)                                                     & 
                               +   CEthanolAG2Loss(i,itile)                                                     & 
                               +   CNonsolubleAG1Loss(i,itile)                                                  & 
                               +   CNonsolubleAG2Loss(i,itile)                                                 
          ENDDO
        ENDDO
      END IF

      ! Fluxes out of anthropogenic pools
      C_construction_2_atmos(:) = Cpool_construction(:) / tau_construction ! ... construction
      C_paper_2_atmos       (:) = Cpool_paper       (:) / tau_paper        ! ... paper
      C_onSite_2_atmos      (:) = Cpool_onSite      (:) / tau_onSite

      ! Carbon loss to atmosphere
      C_2_atmos(:) = C_onSite_2_atmos      (:) + & 
                     C_paper_2_atmos       (:) + & 
                     C_construction_2_atmos(:)
          
      ! New sizes of anthropogenic pools
      Cpool_onSite      (:) = Cpool_onSite      (:) + C_2_onSite      (:) - C_onSite_2_atmos      (:)
      Cpool_paper       (:) = Cpool_paper       (:) + C_2_paper       (:) - C_paper_2_atmos       (:)
      Cpool_construction(:) = Cpool_construction(:) + C_2_construction(:) - C_construction_2_atmos(:)
    ENDIF

    IF (.NOT. PRESENT(C_green_loss)) &
       DEALLOCATE(CGreenLoss, CWoodLoss, CReserveLoss)
    IF (.NOT. with_yasso) THEN
       IF (.NOT. PRESENT(C_litter_green_ag_loss)) &
          DEALLOCATE(CLitterGreenAGLoss, CLitterWoodAGLoss) 
    ELSE
       IF (.NOT. PRESENT(C_acid_ag1_loss))  &
           DEALLOCATE(CAcidAG1Loss,       &
                      CWaterAG1Loss,      &
                      CEthanolAG1Loss,    &
                      CNonsolubleAG1Loss, & 
                      CAcidAG2Loss,       &
                      CWaterAG2Loss,      &
                      CEthanolAG2Loss,    &
                      CNonsolubleAG2Loss)
    END IF
    IF (.NOT. PRESENT(Npool_green_loss)) &
      DEALLOCATE(NGreenLoss,NWoodLoss,NMobileLoss)

  END SUBROUTINE C_loss_and_update_anthro_pools


  PURE SUBROUTINE yasso (Yasso_io_pools, Weather, litter, Lit_coefV,WoodLitterSize, Yasso_out, frac_aboveground, root_exudates)
    !! DSG: 11.03.2013
    !! modified Yasso soil C model using an Euler scheme and readable coding style
    !! JSBACH specific modification of model structure: AWEN pools are separated into aboveground and belowground pools
    !! 
    !! Herbivory flux is treated like normal litter, because the herbivory fluxes in JSBACH are not evaluated, yet.
    !! If the represenation of herbivory in JSBACH is evaluated I recommend to treat herbivory flux in YASSO different from general
    !! litter.
    !!
    !! The routine needs as input 1. daily litter flux
    !!                            2. Yasso pools
    !!                            3. PFT specific parameters
    !!
    !! the subroutine runs on a annual time step only
    !! output is on a daily timestep: respiration
    !! output pools: AWEN + humus 
    !! The routine must be called twice, first for non-woody litter
    !! then for woody litter
    !!-------------------------------------------------------------------------------------------------------------------------
    !! Yasso - Soil carbon model (Jari Liski)
    !!
    !! The code for the yasso model has been developed and is provided by the Finnish Environment Institute SYKE. It is 
    !! distributed as part of the MPI Earth System Model under the MPI-M license conditions as stated in the header of this module. 
    !!
    !! model to calculate the amount of soil organic carbon, changes in the amount of soil organic carbon and
    !! heterotrophic soil respiration.
    !! For documention of yasso see Liski et al 2005, Tuomi et al 2008,2009,2011
    !!
    !! implementation: Tea Thum (FMI), Petri Räisänen (FMI), Daniel Goll (MPI-M)
    !!-------------------------------------------------------------------------------------------------------------------------

    IMPLICIT NONE
    
    REAL(dp), DIMENSION(9),  INTENT(IN)  :: Yasso_io_pools    ! Yasso pools IN                                        [mol(c)/m2]
    REAL(dp),                INTENT(IN)  :: WoodLitterSize    ! Size of woody litter; 
                                                              !     is zero when the subroutine is called for non-woody litter
    REAL(dp), DIMENSION(2),  INTENT(IN)  :: Weather           ! climatic drivers: air temperature and precipitation, 15-d means
    REAL(dp),                INTENT(IN)  :: litter            ! fresh litter inputs  (above and belowground)         [mol(C)/m2/day]
    REAL(dp), DIMENSION(5),  INTENT(IN)  :: Lit_coefV         ! fractions to seperate incoming litter to the yasso pool    [ ]
    REAL(dp),                INTENT(IN)  :: frac_aboveground  ! parameter to seperate litter influx into above- and belowground 
                                                              !     part                                                   [ ]
    REAL(dp),                INTENT(IN)  :: root_exudates     ! root exudates                                        [mol(C)/m2/day]

    REAL(dp), DIMENSION(10), INTENT(OUT) :: Yasso_out         ! updated Yasso pools(1:9) [mol(c)/m2] & respiration(10) [mol(c)/m2/d]

    !local
    ! Yasso pools   [mol(c)/m2]
    REAL(dp)     :: YCpool_acid_ag    ! Aboveground  
    REAL(dp)     :: YCpool_water_ag  
    REAL(dp)     :: YCpool_ethanol_ag 
    REAL(dp)     :: YCpool_nonsoluble_ag
    REAL(dp)     :: YCpool_acid_bg    ! Belowground  
    REAL(dp)     :: YCpool_water_bg  
    REAL(dp)     :: YCpool_ethanol_bg 
    REAL(dp)     :: YCpool_nonsoluble_bg
    
    REAL(dp)     :: YCpool_humus   

    ! Yasso fluxes (internally used, therefore annual fluxes)  
    REAL(dp)     :: Cflx_litter_2_acid             ! Litter influx to acid soluble pools    [mol(c)/m2/a]   
    REAL(dp)     :: Cflx_litter_2_water            ! Litter influx to water soluble pools   [mol(c)/m2/a]   
    REAL(dp)     :: Cflx_litter_2_ethanol          ! Litter influx to ethanol soluble pools [mol(c)/m2/a]   
    REAL(dp)     :: Cflx_litter_2_nonsoluble       ! Litter influx to non soluble pools     [mol(c)/m2/a]   
    REAL(dp)     :: Cflx_litter_2_humus            ! Litter influx to humus soluble pool    [mol(c)/m2/a]   
  
    ! Aboveground
    REAL(dp)     :: Cflx_from_acid_ag              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_water_ag             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_ethanol_ag           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_nonsoluble_ag        ! Loss flux of non soluble pool         [mol(c)/m2/a]
    
    REAL(dp)     :: Cflx_2_acid_ag                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_water_ag                ! Gain flux of water soluble pool       [mol(c)/m2/a]              
    REAL(dp)     :: Cflx_2_ethanol_ag              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_nonsoluble_ag           ! Gain flux of non soluble pool         [mol(c)/m2/a]
     
    ! Belowground
    REAL(dp)     :: Cflx_from_acid_bg              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_water_bg             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_ethanol_bg           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_nonsoluble_bg        ! Loss flux of non soluble pool         [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_humus                ! Loss flux of humus pool               [mol(c)/m2/a]
    
    REAL(dp)     :: Cflx_2_acid_bg                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_water_bg                ! Gain flux of water soluble pool       [mol(c)/m2/a]              
    REAL(dp)     :: Cflx_2_ethanol_bg              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_nonsoluble_bg           ! Gain flux of non soluble pool         [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_humus                   ! Gain flux of humus pool               [mol(c)/m2/a]
    
    ! Respiration is output of Yasso therefore a daily flux 
    REAL(dp)     :: soilResp_rateYasso             ! Flux to the atmosphere                [mol(c)/m2/d]  

    ! Decomposition rates [1/a]
    REAL(dp)     :: d_acid
    REAL(dp)     :: d_water
    REAL(dp)     :: d_ethanol
    REAL(dp)     :: d_nonsoluble
    REAL(dp)     :: d_humus

    ! PFT-specific parameters
    REAL(dp)     :: frac_litter_2_acid        ! Fraction of incoming material which enters the acid-soluble fraction of litter
    REAL(dp)     :: frac_litter_2_water       ! Fraction of incoming material which enters the water-soluble fraction of litter
    REAL(dp)     :: frac_litter_2_ethanol     ! Fraction of incoming material which enters the ethanol-soluble fraction of litter
    REAL(dp)     :: frac_litter_2_nonsoluble  ! Fraction of incoming material which enters the nonsoluble fraction of litter
    REAL(dp)     :: frac_litter_2_humus       ! Fraction of incoming material which enters the humus pool

    ! Yasso parameters (Tuomi '11)
    ! reference decompositation rates for T = 0C and unlimited water availability  in [1/yr]... 
    REAL(dp), PARAMETER   :: ref_decomp_rate_acid       = -0.72    !  ... acid-soluble fraction 
    REAL(dp), PARAMETER   :: ref_decomp_rate_water      = -5.9     !  ... water-soluble fraction 
    REAL(dp), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28    !  ... ethanol-soluble fraction 
    REAL(dp), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031   !  ... nonsoluble-soluble fraction 
    REAL(dp), PARAMETER   :: ref_decomp_rate_humus      = -0.0016  !  ... humus 
    ! parameters for temperature dependence
    REAL(dp), PARAMETER   :: temp_p1                    = 0.095    ! first  order temperature dependence [1/C]
    REAL(dp), PARAMETER   :: temp_p2                    = -0.0014  ! second order temperature dependence [1/C]
    ! parameter for precipitation dependence
    REAL(dp), PARAMETER   :: precip_p1                  = -1.21    ! first  order precipitation dependence [a/m]
    ! parameters for size dependence 
    REAL(dp), PARAMETER   :: size_p1                    = -1.71    ! first  order size dependence [1/cm]
    REAL(dp), PARAMETER   :: size_p2                    = 0.86     ! second order size dependence [1/cm]
    REAL(dp), PARAMETER   :: size_p3                    = -0.306   ! size dependence power [ ]
    
    ! Fractions to split the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere 
    ! REMARK: there is no a priori (to the calibration procedure of YASSO) assumption about possible fluxes between pools, 
    ! therefore all fluxes are theoretically possible, however some of the fraction are zero which means the corresponding flux 
    ! does not exists. However, after a recalibration of the model this fraction may change.
    !
    ! fractions of decomposition flux of acid soluble pool going to the other AWEN pools [ ]
    ! (8,11,14)
    REAL(dp), PARAMETER   :: A_2_W                      =  0.99
    REAL(dp), PARAMETER   :: A_2_E                      =  0.00
    REAL(dp), PARAMETER   :: A_2_N                      =  0.00
    ! fractions of decomposition flux of water soluble pool going to the other AWEN pools [ ]
    ! (5,12,15)
    REAL(dp), PARAMETER   :: W_2_A                      = 0.48 
    REAL(dp), PARAMETER   :: W_2_E                      = 0.00 
    REAL(dp), PARAMETER   :: W_2_N                      = 0.015
    ! fractions of decomposition flux of ethanol soluble pool going to the other AWEN pools [ ]
    ! (6,9,16)
    REAL(dp), PARAMETER   :: E_2_A                      = 0.01 
    REAL(dp), PARAMETER   :: E_2_W                      = 0.00
    REAL(dp), PARAMETER   :: E_2_N                      = 0.95
    ! fractions of decomposition flux of non soluble pool going to the other AWEN pools [ ]
    ! (7,10,13)
    REAL(dp), PARAMETER   :: N_2_A                      = 0.83
    REAL(dp), PARAMETER   :: N_2_W                      = 0.01
    REAL(dp), PARAMETER   :: N_2_E                      = 0.02
    ! fraction of decomposition fluxes of AWEN pools which enters the humus pool [ ]
    REAL(dp), PARAMETER   :: AWEN_2_H                   = 0.0045

    ! climatic drivers
    REAL(dp)     :: precip                  ! 15 runnning mean of precipitation      [m/a]
    REAL(dp)     :: temp                    ! 15 runnning mean of 2m air temperature [C]

    REAL(dp)     :: d_temp                  ! term which accounts for the temperature influence on litter decomposition
    REAL(dp)     :: d_precip                ! term which accounts for the precipitation influence on litter decomposition
    REAL(dp)     :: d_size                  ! term which accounts for the litter size influence on litter decomposition
     

    ! time stepping; the parameterization of yasso is done on a annual time
    ! step, thus the fluxes have to be scaled up to annual rates. 
    REAL(dp)     :: dt                      ! time step of mo_update_CPools [year]             
    dt = 1.0/days_per_year   

    ! 0. Preparations
     ! initialize parameters
     frac_litter_2_acid       = Lit_coefV(1)
     frac_litter_2_water      = Lit_coefV(2)
     frac_litter_2_ethanol    = Lit_coefV(3)
     frac_litter_2_nonsoluble = Lit_coefV(4)
     frac_litter_2_humus      = Lit_coefV(5)

     ! initialize yasso pools (AWEN + humus)
     YCpool_acid_ag        = Yasso_io_pools(1)
     YCpool_water_ag       = Yasso_io_pools(2)
     YCpool_ethanol_ag     = Yasso_io_pools(3)
     YCpool_nonsoluble_ag  = Yasso_io_pools(4)
     
     YCpool_acid_bg        = Yasso_io_pools(5)
     YCpool_water_bg       = Yasso_io_pools(6)
     YCpool_ethanol_bg     = Yasso_io_pools(7)
     YCpool_nonsoluble_bg  = Yasso_io_pools(8)
     
     YCpool_humus          = Yasso_io_pools(9)
     
     ! Change units of climatic forcing variables
     precip = Weather(2)*sec_per_year/1000._dp  ! mm/s -> m/a 
     temp   = Weather(1) - Tmelt                ! K -> C
    
     ! Calculate annual litter influxes [mol(c)/m2/a]
     Cflx_litter_2_acid       = frac_litter_2_acid        *litter * days_per_year
     Cflx_litter_2_water      = frac_litter_2_water       *litter * days_per_year
     Cflx_litter_2_ethanol    = frac_litter_2_ethanol     *litter * days_per_year
     Cflx_litter_2_nonsoluble = frac_litter_2_nonsoluble  *litter * days_per_year
     Cflx_litter_2_humus      = frac_litter_2_humus       *litter * days_per_year


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! 1. Calculate decomposition rates 
     
     ! Temperature dependence of decomposition
     d_temp   = EXP(temp_p1*temp+temp_p2*temp**2.0)  
     ! Precipitation dependence of decomposition
     d_precip = 1.0-EXP(precip_p1*precip)
     ! Litter size dependence of decomposition -- no effect if WoodlitterSize = 0.0 
     d_size   = MIN(1.0_dp,(1.0 + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0)**size_p3)

     ! decomposition rates accounting for temperature, precipitation and the litter size 
     d_acid       =  ref_decomp_rate_acid       *d_temp *d_precip *d_size
     d_water      =  ref_decomp_rate_water      *d_temp *d_precip *d_size
     d_ethanol    =  ref_decomp_rate_ethanol    *d_temp *d_precip *d_size
     d_nonsoluble =  ref_decomp_rate_nonsoluble *d_temp *d_precip *d_size
     d_humus      =  ref_decomp_rate_humus      *d_temp *d_precip           ! no size effect on humus 
     
    ! 2. Calculate fluxes 
     
     ! loss fluxes (negative values): 
     Cflx_from_acid_ag         = YCpool_acid_ag       * d_acid        
     Cflx_from_water_ag        = YCpool_water_ag      * d_water                  
     Cflx_from_ethanol_ag      = YCpool_ethanol_ag    * d_ethanol             
     Cflx_from_nonsoluble_ag   = YCpool_nonsoluble_ag * d_nonsoluble       
     
     Cflx_from_acid_bg         = YCpool_acid_bg       * d_acid        
     Cflx_from_water_bg        = YCpool_water_bg      * d_water                  
     Cflx_from_ethanol_bg      = YCpool_ethanol_bg    * d_ethanol             
     Cflx_from_nonsoluble_bg   = YCpool_nonsoluble_bg * d_nonsoluble       
     
     Cflx_from_humus           = YCpool_humus         * d_humus         

     ! gain fluxes (positive values): 
     ! fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why) 
     Cflx_2_acid_ag           = ABS(Cflx_from_water_ag      * W_2_A  &    ! returns positive fluxes
                                  + Cflx_from_ethanol_ag    * E_2_A  &
                                  + Cflx_from_nonsoluble_ag * N_2_A) 

     Cflx_2_water_ag          = ABS(Cflx_from_acid_ag       * A_2_W  &
                                  + Cflx_from_ethanol_ag    * E_2_W  &
                                  + Cflx_from_nonsoluble_ag * N_2_W) 

     Cflx_2_ethanol_ag        = ABS(Cflx_from_acid_ag       * A_2_E  &
                                  + Cflx_from_water_ag      * W_2_E  &
                                  + Cflx_from_nonsoluble_ag * N_2_E) 

     Cflx_2_nonsoluble_ag     = ABS(Cflx_from_acid_ag       * A_2_N  &
                                  + Cflx_from_water_ag      * W_2_N  &
                                  + Cflx_from_ethanol_ag    * E_2_N) 

     Cflx_2_acid_bg           = ABS(Cflx_from_water_bg      * W_2_A  &    
                                  + Cflx_from_ethanol_bg    * E_2_A  &
                                  + Cflx_from_nonsoluble_bg * N_2_A) 

     Cflx_2_water_bg          = ABS(Cflx_from_acid_bg       * A_2_W  &
                                  + Cflx_from_ethanol_bg    * E_2_W  &
                                  + Cflx_from_nonsoluble_bg * N_2_W) 

     Cflx_2_ethanol_bg        = ABS(Cflx_from_acid_bg       * A_2_E  &
                                  + Cflx_from_water_bg      * W_2_E  &
                                  + Cflx_from_nonsoluble_bg * N_2_E) 

     Cflx_2_nonsoluble_bg     = ABS(Cflx_from_acid_bg       * A_2_N  &
                                  + Cflx_from_water_bg      * W_2_N  &
                                  + Cflx_from_ethanol_bg    * E_2_N) 

     Cflx_2_humus          = ABS(Cflx_from_acid_ag                &
                               + Cflx_from_water_ag               &
                               + Cflx_from_ethanol_ag             &
                               + Cflx_from_nonsoluble_ag          &
                               + Cflx_from_acid_bg                &
                               + Cflx_from_water_bg               &
                               + Cflx_from_ethanol_bg             &
                               + Cflx_from_nonsoluble_bg          &
                               ) * AWEN_2_H  

     ! the remaining fractions of the loss fluxes enter the atmosphere as respiration
     soilResp_rateYasso    =  (Cflx_from_acid_ag + Cflx_from_acid_bg)               &
                              * (1.0_dp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                              + (Cflx_from_water_ag + Cflx_from_water_bg)           &      
                              * (1.0_dp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                              + (Cflx_from_ethanol_ag + Cflx_from_ethanol_bg)       &     
                              * (1.0_dp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                              + (Cflx_from_nonsoluble_ag + Cflx_from_nonsoluble_bg) &  
                              * (1.0_dp - N_2_A - N_2_W - N_2_E - AWEN_2_H)         &
                              + Cflx_from_humus 

   ! 3. update Yasso pools
     YCpool_acid_ag           = MAX(0.0_dp,                                         &
                              YCpool_acid_ag                                        & ! old pool   
                              + (Cflx_from_acid_ag                                  & ! incoming flux from AWEN pools
                                + Cflx_2_acid_ag                                    & ! outgoing flux
                                + frac_aboveground * Cflx_litter_2_acid) *dt)         ! incoming flux from litter
     YCpool_water_ag          = MAX(0.0_dp,                                         & ! and so on .....
                              YCpool_water_ag                                       &
                              + (Cflx_from_water_ag                                 &
                                + Cflx_2_water_ag                                   &
                                + frac_aboveground * Cflx_litter_2_water) *dt)
     YCpool_ethanol_ag        = MAX(0.0_dp,                                         &
                              YCpool_ethanol_ag                                     &
                              + (Cflx_from_ethanol_ag                               &
                                + Cflx_2_ethanol_ag                                 &
                                + frac_aboveground * Cflx_litter_2_ethanol) *dt)
     YCpool_nonsoluble_ag     = MAX(0.0_dp,                                         &
                              YCpool_nonsoluble_ag                                  &
                              + (Cflx_from_nonsoluble_ag                            &
                                + Cflx_2_nonsoluble_ag                              &
                                + frac_aboveground * Cflx_litter_2_nonsoluble) *dt)
    

     YCpool_acid_bg           = MAX(0.0_dp,                                         &
                              YCpool_acid_bg                                        & ! old pool   
                              + (Cflx_from_acid_bg                                  & ! incoming flux from AWEN pools
                                + Cflx_2_acid_bg                                    & ! outgoing flux
                                + (1._dp - frac_aboveground) *                      & ! 
                                  Cflx_litter_2_acid) *dt)                            ! incoming flux from litter
     YCpool_water_bg          = MAX(0.0_dp,                                         & ! ...
                              YCpool_water_bg                                       &
                              + (Cflx_from_water_bg                                 &
                                + Cflx_2_water_bg                                   &
                                + (1._dp - frac_aboveground) *                      &
                                  Cflx_litter_2_water) *dt                          &
                              + root_exudates)                                    ! exudates are carbohydrates only and belowground 
     YCpool_ethanol_bg        = MAX(0.0_dp,                                         &
                              YCpool_ethanol_bg                                     &
                              + (Cflx_from_ethanol_bg                               &
                                + Cflx_2_ethanol_bg                                 &
                                + (1._dp - frac_aboveground) *                      &
                                     Cflx_litter_2_ethanol) *dt)
     YCpool_nonsoluble_bg     = MAX(0.0_dp,                                         &
                              YCpool_nonsoluble_bg                                  &
                              + (Cflx_from_nonsoluble_bg                            &
                                + Cflx_2_nonsoluble_bg                              &
                                + (1._dp - frac_aboveground) *                      &
                                  Cflx_litter_2_nonsoluble) *dt)

     YCpool_humus          = MAX(0.0_dp,                            &
                              YCpool_humus                          &
                              + (Cflx_from_humus                    &  
                                + Cflx_2_humus                      &
                                + Cflx_litter_2_humus) *dt)

   
   ! 4.1 Write pools into the output array
    Yasso_out(1)           = YCpool_acid_ag 
    Yasso_out(2)           = YCpool_water_ag
    Yasso_out(3)           = YCpool_ethanol_ag
    Yasso_out(4)           = YCpool_nonsoluble_ag
    Yasso_out(5)           = YCpool_acid_bg 
    Yasso_out(6)           = YCpool_water_bg
    Yasso_out(7)           = YCpool_ethanol_bg
    Yasso_out(8)           = YCpool_nonsoluble_bg
    Yasso_out(9)           = YCPool_humus
   ! 4.2 Write respiration into the output array
    Yasso_out(10)           = soilResp_rateYasso * dt ! convert back to daily

  END SUBROUTINE yasso

end module mo_cbal_cpools
