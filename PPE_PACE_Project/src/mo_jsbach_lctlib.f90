!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_jsbach_lctlib
  !
  ! This module provides the data from the landcover library file, which contains additional information on the 
  ! landcover data used in a particular run of JSBACH.
  ! The name of the landcover library file is specified in the JSBACH configuration file (run.def) under the 
  ! keyword "LCT_FILE".
  !
  ! Authors: Reiner Schnur, 7/2003
  !          Christian H. Reick, 8/2003
  !
  USE mo_mpi,       ONLY: p_parallel, p_parallel_io, p_bcast, p_io
  USE mo_kind,      ONLY: dp
  USE mo_exception, ONLY: message_text, finish
  USE mo_netcdf,    ONLY: nf_max_name

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART =======================================================================================================

  PUBLIC :: init_lctlib ! This subroutine initializes lctlib, i.e. it allocates memory and reads in the data


  ! --- lctlib ---------------------------------------------------------------------------------------------------------------------

  TYPE lctlib_type
     ! --- general parameters ------------------------------------------------------------------------------------------------------
     INTEGER           :: nlct                    !! Number of landcover types
     INTEGER           :: npft                    !! Number of vegetation types, including crops and pastures
     INTEGER, POINTER  :: LctNumber(:)            !! Unique number for each landcover type (1, 2, 3, ..)
     CHARACTER(LEN=16), POINTER  :: LctName(:)    !! Names of landcover types
     INTEGER, POINTER  :: LandcoverClass(:)       !! Landcover classes (not to be confused with landcover types!):
                                                  !! .. 0: Bare soil; 1:Glacier; 2: Lake; 
                                                  !! .. 3: Natural forest; 4: Natural Grassland; 
                                                  !! .. 5: Other natural vegetation; 6: crops; 7: pastures
     LOGICAL, POINTER  :: NaturalVegFlag(:)       !! "true" for natural vegetation landcover types, i.e. non-agricultural vegetation
     LOGICAL, POINTER  :: ForestFlag(:)           !! "true" for forest landcover types
     LOGICAL, POINTER  :: GrassFlag(:)            !! "True" for natural graslands     
     LOGICAL, POINTER  :: CropFlag(:)             !! "true" for croplands
     LOGICAL, POINTER  :: PastureFlag(:)          !! "true" for landcover of type "pasture"
     LOGICAL, POINTER  :: LakeFlag(:)             !! "true" for a land surface class that signifies lake or ocean
     LOGICAL, POINTER  :: GlacierFlag(:)          !! "true" for a land surface class that signifies glaciers
     LOGICAL, POINTER  :: BareSoilFlag(:)         !! "true" if none of the foregoing landcover classes

     ! --- Albedo ------------------------------------------------------------------------------------------------------------------
     REAL(dp), POINTER :: AlbedoSnowVisMin(:)      !! Minimum snow albedo in the visible range
     REAL(dp), POINTER :: AlbedoSnowVisMax(:)      !! Maximum snow albedo in the visible range
     REAL(dp), POINTER :: AlbedoSnowNirMin(:)      !! Minimum snow albedo in the NIR range
     REAL(dp), POINTER :: AlbedoSnowNirMax(:)      !! Maximum snow albedo in the NIR range
     REAL(dp), POINTER :: AlbedoSnowMin(:)         !! Minimum snow albedo
     REAL(dp), POINTER :: AlbedoSnowMax(:)         !! Maximum snow albedo
     REAL(dp), POINTER :: AlbedoCanopyVIS(:)       !! Albedo of the canopy (vegetation) in the visible range
     REAL(dp), POINTER :: AlbedoCanopyNIR(:)       !! Albedo of the canopy (vegetation) in the NIR range
     REAL(dp), POINTER :: AlbedoLitterVis(:)       !! Albedo of the leaf litter in the visible range
     REAL(dp), POINTER :: AlbedoLitterNir(:)       !! Albedo of the leaf litter in the NIR range

     ! --- parameters used by bethy ------------------------------------------------------------------------------------------------
     LOGICAL, POINTER  :: NitrogenScalingFlag(:)   !! Indicates, whether nitrogen scaling shall be applied to that vegetation type
     LOGICAL, POINTER  :: C4flag(:)                !! Photosynthetic pathway: C4=.true or C3=.false.
     REAL(dp), POINTER :: CarboxRate(:)            !! Maximum carboxilation rate at 25 degrees Celsius [1.E-6 * Mol(CO2)/m^2/s] ...
                                                   !! ... (Table 2.6 in Knorr)
     REAL(dp), POINTER :: ETransport(:)            !! Maximum electron transport rate at 25 degrees Celsius [1.E-6 * Mol/m^2/s] ...
                                                   !! ... (Table 2.6 in Knorr)
     REAL(dp), POINTER :: VegHeight(:)             !! Typical height of the vegetation classes [m]
     REAL(dp), POINTER :: VegRoughness(:)          !! Roughness length of the vegetation classes [m]
     REAL(dp), POINTER :: MinVegRoughness(:)       !! Minimal roughness length of the vegetation classes (LAI = 0) [m]
     REAL(dp), POINTER :: MaxVegRoughness(:)       !! Maximal roughness length of the vegetation classes (LAI = LAI_max) [m]
     ! --- Parameters used in LoGroP Phenology scheme --------------------------------------------------------------------------
     INTEGER, POINTER  :: PhenologyType(:)         !! Phenology type (only for natural vegetation):
                                                   !! ... none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4; crops: 5
     REAL(dp), POINTER :: MaxLAI(:)                !! Upper LAI boundary for LoGoP-Scheme (phenology) in [m2/m2]
     REAL(dp), POINTER :: specificLeafArea_C(:)    !! Carbon content per leaf area in [mol(Carbon)/m^2(leaf)]
 
     ! --- Parameters used in Knorr Phenology scheme ---------------------------------------------------------------------------
     REAL(dp), POINTER :: knorr_Tau_w(:)           !! Time before leaf shedding  [days]
     REAL(dp), POINTER :: knorr_T_phi(:)           !! Temperature trigger for leaf growth [deg C]
     REAL(dp), POINTER :: knorr_T_r(:)             !! "Spread" (sigma) of T_phi [deg C]
     REAL(dp), POINTER :: knorr_Day_c(:)           !! Day-length at leaf shedding [hours]
     REAL(dp), POINTER :: knorr_Day_r(:)           !! "Spread" (sigma) of Day_c [hours]
     REAL(dp), POINTER :: knorr_k_l(:)             !! Inverse of leaf longevity [days-1]
     REAL(dp), POINTER :: knorr_leaf_growth_rate(:) !! Initial leaf growth rate [days-1]
     REAL(dp), POINTER :: knorr_max_lai(:)         !! Maximum LAI

     ! --- Parameters used in the carbon balance model -----------------------------------------------------------------------------
     REAL(dp), POINTER :: reserveC2leafC(:)        !! Ratio of opt. C content of the reserve pool to the opt. C content of leaves 
     REAL(dp), POINTER :: frac_npp_2_woodPool(:)   !! Maximum fraction of NPP put into the wood pool of the carbon balance model
     REAL(dp), POINTER :: frac_npp_2_reservePool(:)!! Optimal fraction of NPP put into the reserve pool of the carbon balance model
     REAL(dp), POINTER :: frac_npp_2_exudates(:)   !! Optimal fraction of NPP put into the root exudates of the carbon balance model
     REAL(dp), POINTER :: frac_green_2_herbivory(:)!! Fraction of NPP that is grazed by herbivores
     REAL(dp), POINTER :: tau_Cpool_litter_leaf(:) !! Time constant by which leaf litter pool is depreciated  [days]
     REAL(dp), POINTER :: tau_Cpool_litter_wood(:) !! Time constant by which woody litter pool is depreciated  [days]
     REAL(dp), POINTER :: LAI_shed_constant(:)     !! 
     REAL(dp), POINTER :: frac_C_litter_green2atmos(:) !! Fraction of heterotrophic loss of green litter pools emitted to the 
                                                       !! .. atmosphere (rest enters slow pool)
     REAL(dp), POINTER :: Max_C_content_woods(:)   !! Maximum carbon content in the woody parts of plants [mol(C)/m^2(canopy)] (for
                                                   !! .. forests this is closely related to the yield, usually measured in 
                                                   !! .. m^3(wood)/hectar)
     REAL(dp), POINTER :: ClumpinessFactor(:)      !! Factor to calculate vegetation clumpiness: 
                                                   !!    veg_ratio=veg_ratio_max*(1-exp(-LAI_max/ClumpinessFactor))
                                                   !! converting natural PFTs into agricultural PFTs
     REAL(dp), POINTER :: frac_wood_2_onSite(:)    !! Fraction of wood pool to anthropogenically controlled onSite pool wenn
                                                   !! converting natural PFTs into agricultural PFTs
     REAL(dp), POINTER :: frac_wood_2_paper(:)     !! Fraction of wood pool to paper (intermediately longlived) pool wenn...
     REAL(dp), POINTER :: frac_wood_2_construction(:) !! Fract of wood pool to constructions (e.g. houses, furniture = longlived)
                                                   !!   pool when converting natural PFTs into agricultural PFTs

     ! --- Parameters used in the carbon balance model and the dynamic vegetation --------------------------------------------------
     REAL(dp), POINTER :: tau_Cpool_woods(:)           !! PFT-specific time scale for woody (lignified) plant tissue 

     ! --- Parameters used by the dynamic vegetation -------------------------------------------------------------------------------
     LOGICAL,  POINTER :: dynamic_PFT(:)               !! indicates those PFTs which shall take part in the vegetation dynamics
     LOGICAL,  POINTER :: woody_PFT(:)                 !! indicates those PFTs which are of woody type (in contrast to grasses)
     LOGICAL,  POINTER :: pasture_PFT(:)               !! indicates those PFTs which are pasture
     REAL(dp), POINTER :: bclimit_min_cold_mmtemp(:)   !! PFT-specific minimum coldest monthly mean temperature
     REAL(dp), POINTER :: bclimit_max_cold_mmtemp(:)   !! PFT-specific maximum coldest monthly mean temperature
     REAL(dp), POINTER :: bclimit_max_warm_mmtemp(:)   !! PFT-specific upper limit of warmest-month temperature
     REAL(dp), POINTER :: bclimit_min_temprange(:)     !! PFT-specific 20-year average min warmest - coldest month temperature range
     REAL(dp), POINTER :: bclimit_min_gdd(:)           !! PFT-specific minimum growing degree days (at or above 5 deg C)

     ! --- Parameters used in SPITFIRE (mo_disturbance_thonicke.f90) ---------------------------------------------------------------
     REAL(dp), POINTER :: moist_extinction(:)          !! PFT specific moisture of extinction
     REAL(dp), POINTER :: fuel_dens(:)                 !! fuel bulk density
     REAL(dp), POINTER :: flame_length_f(:)            !! f parameter for flame length(scorch height) parameter 45 in lpj
     REAL(dp), POINTER :: crown_length(:)              !! crown length parameter see table 1 in thonicke et al. 2010
     REAL(dp), POINTER :: bark_par1(:)                 !! bark thickness parameter 1 see table 1 and eq. 21 in thonicke et al. 2010
     REAL(dp), POINTER :: bark_par2(:)                 !! bark thickness parameter 2 see table 1 and eq. 21 in thonicke et al. 2010
     REAL(dp), POINTER :: RCK(:)                       !! resistance factor to crown damage, tab1 and eq22 in thonicke et al. 2010
     REAL(dp), POINTER :: mort_prob(:)                 !! mortality probability  see table 1 and eq. 22 in thonicke et al. 2010

     ! --- Parameters for the climate buffer ---------------------------------------------------------------------------------------
     REAL(dp), POINTER :: gdd_base(:)                  !! PFT-specific GDD base
     REAL(dp), POINTER :: upper_tlim(:)                !! PFT-specific base to calculate GDD_upper_tlim 

     ! --- Parameters used by the yasso soil carbon model -------------------------------------------------------------------------
     REAL(dp), POINTER :: LitVeg_coef(:,:)        !! Yasso coefficient for separating litter into chemical pools
     REAL(dp), POINTER :: LeafLit_coef(:,:)       !! Yasso coefficient for separating leaf litter into chemical pools
     REAL(dp), POINTER :: WoodLit_coef(:,:)       !! Yasso coefficient for separating woody litter into chemical pools
     REAL(dp), POINTER :: WoodLitterSize(:)       !! Yasso: size of the litter

     ! --- Other parameters --------------------------------------------------------------------------------------------------------
     REAL(dp), POINTER :: CanopyResistanceMin(:)  !! Minimum canopy resistance (optional, only used in VIC scheme
                                                  !! and if BETHY is not used)
     REAL(dp), POINTER :: StemArea(:)             !! Area of stems and branches of woody plants

  END TYPE lctlib_type
  PUBLIC :: lctlib_type                            

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE 


  ! Save nlct and npft for convenience
  INTEGER, SAVE             :: nlct ! The number of land cover types
  INTEGER, SAVE             :: npft ! Number of plant functional types, including crops and pastures

  ! === Private declarations =======================================================================================================

CONTAINS

  SUBROUTINE init_lctlib(lctlib_file_name, lctlib)
    !
    ! Get landcover library, i.e. lookup table for landcover types
    !
    ! The structure of the lookup table in the input file is as follows. It can contain only three types of lines in arbitrary
    ! order:
    !
    !     Comment lines: These contain before a '#' only blanks
    !       Blank lines: These contain only blanks
    !        Data lines: The first nonblank characters form a keyword which is follwed by the data (separated by blanks)
    !
    ! Keywords are the same as the names of the lctlib components, with one exception:
    !
    !        NLCT : This is the keyword after which the number of landcover types, i.e. the number of data columns
    !               in the file is expected. This number has to be identical with the number of landcover types
    !               from the landcover data file (keyword: LCTLIB_FILE in the JSBACH configuration file).
    !               
    ! After all other keywords NLCT columns of data are expected.
    !
    ! Example:
    !            # ---- LANDCOVER LIBRARY FILE -------
    !            NLCT 3
    !            LctNumber               3 5 9
    !            # Phenology types: 0= none 1=summergreen, 2=evergreen, 3=raingreen, 4=grasses
    !            PhenologyType          2 2 4
    !            # C4flag: 0=C3, 1=C4
    !            C4flag                 0 1 1
    !
    ! The first string on each line (except first line and comment lines) must correspond to the name of the lctlib component
    ! (case doesn't matter)
    ! For lctlib components of type LOGICAL use 0/1 in the file to indicate .FALSE./.TRUE.
    ! For lctlib components of type INTEGER the numbers on the line must be integer values
    ! For lctlib components of type REAL the numbers on the line can be either integer or floating point
    ! Comments can appear on any line, everything to the right of and including "#" is disregarded
    !
    ! The file can contain more keywords than needed --- therefore the same file can be used by several model components.
    !
    USE mo_util_string, ONLY: tolower

    CHARACTER(nf_max_name), INTENT(in) :: lctlib_file_name
    TYPE(lctlib_type),    INTENT(out) :: lctlib ! LCT library to initialize

    INTEGER, PARAMETER            :: lctlib_file_unit = 66

    CHARACTER(len=30)  :: key
    CHARACTER(len=256) :: line
    INTEGER            :: pos_comment, read_status
    INTEGER            :: pos,length
    CHARACTER(len=2)   :: blank_set = " "//achar(9) !! the blank characters: BLANK and TAB
    INTEGER,ALLOCATABLE:: itmp(:)                   !! temporary array used for input of logicals
    INTEGER            :: apu  ! THT 30.10.2012

    LOGICAL            :: exist_LctNumber              = .FALSE.
    LOGICAL            :: exist_LctName                = .FALSE.
    LOGICAL            :: exist_LandcoverClass         = .FALSE.
    LOGICAL            :: exist_AlbedoSnowVisMin       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowVisMax       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowNirMin       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowNirMax       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowMin          = .FALSE.
    LOGICAL            :: exist_AlbedoSnowMax          = .FALSE.
    LOGICAL            :: exist_AlbedoCanopyVIS        = .FALSE.
    LOGICAL            :: exist_AlbedoCanopyNIR        = .FALSE.
    LOGICAL            :: exist_AlbedoLitterVis        = .FALSE.
    LOGICAL            :: exist_AlbedoLitterNir        = .FALSE.
    LOGICAL            :: exist_nitrogenScalingFlag    = .FALSE.
    LOGICAL            :: exist_C4flag                 = .FALSE.
    LOGICAL            :: exist_CarboxRate             = .FALSE.
    LOGICAL            :: exist_ETransport             = .FALSE.
    LOGICAL            :: exist_VegHeight              = .FALSE.
    LOGICAL            :: exist_VegRoughness           = .FALSE.
    LOGICAL            :: exist_MinVegRoughness        = .FALSE.
    LOGICAL            :: exist_MaxVegRoughness        = .FALSE.
    LOGICAL            :: exist_PhenologyType          = .FALSE.
    LOGICAL            :: exist_CanopyResistanceMin    = .FALSE.
    LOGICAL            :: exist_MaxLAI                 = .FALSE.
    LOGICAL            :: exist_StemArea               = .FALSE.
    LOGICAL            :: exist_specificLeafArea_C     = .FALSE.
    LOGICAL            :: exist_knorr_Tau_w            = .FALSE.
    LOGICAL            :: exist_knorr_T_phi            = .FALSE.
    LOGICAL            :: exist_knorr_T_r              = .FALSE.
    LOGICAL            :: exist_knorr_Day_c            = .FALSE.
    LOGICAL            :: exist_knorr_Day_r            = .FALSE.
    LOGICAL            :: exist_knorr_k_l              = .FALSE.
    LOGICAL            :: exist_knorr_leaf_growth_rate = .FALSE.
    LOGICAL            :: exist_knorr_max_lai          = .FALSE.
    LOGICAL            :: exist_reserveC2leafC         = .FALSE.
    LOGICAL            :: exist_frac_npp_2_woodPool    = .FALSE.
    LOGICAL            :: exist_frac_npp_2_reservePool = .FALSE.
    LOGICAL            :: exist_frac_npp_2_exudates    = .FALSE.
    LOGICAL            :: exist_frac_green_2_herbivory = .FALSE.
    LOGICAL            :: exist_tau_Cpool_litter_leaf  = .FALSE.
    LOGICAL            :: exist_tau_Cpool_litter_wood  = .FALSE.
    LOGICAL            :: exist_tau_Cpool_woods        = .FALSE.
    LOGICAL            :: exist_LAI_shed_constant      = .FALSE.
    LOGICAL            :: exist_frac_C_litter_green2atmos = .FALSE.
    LOGICAL            :: exist_Max_C_content_woods       = .FALSE.
    LOGICAL            :: exist_ClumpinessFactor          = .FALSE.
    LOGICAL            :: exists_dynamic_PFT              = .FALSE.
    LOGICAL            :: exists_woody_PFT                = .FALSE.
    LOGICAL            :: exists_pasture_PFT              = .FALSE.
    LOGICAL            :: exists_bclimit_min_cold_mmtemp  = .FALSE.
    LOGICAL            :: exists_bclimit_max_cold_mmtemp  = .FALSE.
    LOGICAL            :: exists_bclimit_max_warm_mmtemp  = .FALSE.
    LOGICAL            :: exists_bclimit_min_temprange    = .FALSE.
    LOGICAL            :: exists_bclimit_min_gdd          = .FALSE.
    LOGICAL            :: exists_gdd_base                 = .FALSE.
    LOGICAL            :: exists_upper_tlim               = .FALSE.
    LOGICAL            :: exists_frac_wood_2_onSite       = .FALSE.
    LOGICAL            :: exists_frac_wood_2_paper        = .FALSE.
    LOGICAL            :: exists_frac_wood_2_construction = .FALSE.
    LOGICAL            :: exists_moist_extinction      = .FALSE.
    LOGICAL            :: exists_fuel_dens             = .FALSE.
    LOGICAL            :: exists_flame_length_f        = .FALSE.
    LOGICAL            :: exists_crown_length          = .FALSE.
    LOGICAL            :: exists_bark_par1             = .FALSE.
    LOGICAL            :: exists_bark_par2             = .FALSE.
    LOGICAL            :: exists_RCK                   = .FALSE.
    LOGICAL            :: exists_mort_prob             = .FALSE.
    LOGICAL            :: exists_woodlittersize_coef      = .FALSE.
    LOGICAL            :: exists_woodlit_coef             = .FALSE.
    LOGICAL            :: exists_leaflit_coef             = .FALSE.
    LOGICAL            :: exists_litveg_coef              = .FALSE.
    INTEGER            :: i

    ! Determine name of landcover library file

    ! Read lctlib data and eventually allocate also memory for lctlib data on io-processor

    IF (p_parallel_io) THEN

       OPEN(unit=lctlib_file_unit, file=lctlib_file_name, form='FORMATTED', status='OLD', iostat=read_status)
       IF (read_status /= 0) CALL finish('init_lctlib','Error opening landcover library file')

       ! --- find keyword NLCT in landcover library file

       DO
          READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) THEN
             CALL finish('init_lctlib','No keyword NLCT found in land cover library file '//TRIM(lctlib_file_name))
          END IF
          ! Look for comment
          pos_comment = SCAN(line,"#")                 ! Start position of comment, 0 if none present
          IF (pos_comment > 1) THEN
             line = TRIM(ADJUSTL(line(1:pos_comment))) ! Disregard everything to the right of "#", ..
             ! .. adjust to the left and trim to the right
          ELSE IF (pos_comment == 1) THEN
             CYCLE                                     ! Disregard whole line
          ELSE
             line = TRIM(ADJUSTL(line))                ! Only disregard blanks to the left and right
          ENDIF
          length=LEN_TRIM(line)
          IF(length== 0) CYCLE                 ! Line is empty
          
          pos = SCAN(line,blank_set)                   ! Position of first blank character
          IF(pos == 0)   CALL finish('init_lctlib',"Wrong syntax in lctlib definitions")
          
          READ(line(1:pos-1),'(A)') key
          
          IF(tolower(TRIM(key)) == 'nlct') THEN
             READ(line(pos:length),*,IOSTAT=read_status) nlct
             IF (read_status /= 0) THEN
                CALL finish('init_lctlib','Could not read number of landcover types (keyword: NLCT) from '//TRIM(lctlib_file_name))
             END IF
             EXIT ! found number of landcover types in landcover library file --- continue after loop
          END IF
       END DO

    END IF

    IF (p_parallel) CALL p_bcast(nlct, p_io)
    lctlib%nlct = nlct

    ALLOCATE(lctlib%LctNumber                (nlct), lctlib%LctName                (nlct), &
             lctlib%LandcoverClass           (nlct), lctlib%NaturalVegFlag         (nlct), &
             lctlib%ForestFlag               (nlct), lctlib%GrassFlag              (nlct), &
             lctlib%CropFlag                 (nlct), lctlib%PastureFlag            (nlct), &
             lctlib%LakeFlag                 (nlct), lctlib%GlacierFlag            (nlct), &
             lctlib%BareSoilFlag             (nlct), lctlib%AlbedoSnowVisMin       (nlct), &
             lctlib%AlbedoSnowVisMax         (nlct), lctlib%AlbedoSnowNirMin       (nlct), &
             lctlib%AlbedoSnowNirMax         (nlct), lctlib%AlbedoSnowMin          (nlct), &
             lctlib%AlbedoSnowMax            (nlct), lctlib%AlbedoCanopyVIS        (nlct), &
             lctlib%AlbedoCanopyNIR          (nlct), lctlib%AlbedoLitterVis        (nlct), &
             lctlib%AlbedoLitterNir          (nlct), lctlib%nitrogenScalingFlag    (nlct), &
             lctlib%C4flag                   (nlct), lctlib%CarboxRate             (nlct), &
             lctlib%ETransport               (nlct), lctlib%VegHeight              (nlct), &
             lctlib%VegRoughness             (nlct), lctlib%MinVegRoughness        (nlct), &
             lctlib%MaxVegRoughness          (nlct), lctlib%PhenologyType          (nlct), &
             lctlib%CanopyResistanceMin      (nlct), lctlib%MaxLAI                 (nlct), &
             lctlib%StemArea                 (nlct), lctlib%specificLeafArea_C     (nlct), &
             lctlib%knorr_Tau_w              (nlct), lctlib%knorr_T_phi            (nlct), &
             lctlib%knorr_T_r                (nlct), lctlib%knorr_Day_c            (nlct), &
             lctlib%knorr_Day_r              (nlct), lctlib%knorr_k_l              (nlct), &
             lctlib%knorr_leaf_growth_rate   (nlct), lctlib%knorr_max_lai          (nlct), &
             lctlib%reserveC2leafC           (nlct), lctlib%frac_npp_2_woodPool    (nlct), &
             lctlib%frac_npp_2_reservePool   (nlct), lctlib%frac_npp_2_exudates    (nlct), &
             lctlib%frac_green_2_herbivory   (nlct), lctlib%tau_Cpool_litter_leaf  (nlct), &
             lctlib%tau_Cpool_litter_wood    (nlct), lctlib%LAI_shed_constant      (nlct), &
             lctlib%frac_C_litter_green2atmos(nlct), lctlib%Max_C_content_woods    (nlct), &
             lctlib%ClumpinessFactor         (nlct), lctlib%dynamic_PFT            (nlct), &
             lctlib%woody_PFT                (nlct), lctlib%pasture_PFT            (nlct), &
             lctlib%bclimit_min_cold_mmtemp  (nlct), lctlib%bclimit_max_cold_mmtemp(nlct), &
             lctlib%bclimit_max_warm_mmtemp  (nlct), lctlib%bclimit_min_temprange  (nlct), &
             lctlib%bclimit_min_gdd          (nlct), lctlib%gdd_base               (nlct), &
             lctlib%upper_tlim               (nlct), lctlib%tau_Cpool_woods        (nlct), &
             lctlib%frac_wood_2_onSite       (nlct), lctlib%frac_wood_2_paper      (nlct), &
             lctlib%frac_wood_2_construction (nlct), lctlib%moist_extinction       (nlct), &
             lctlib%flame_length_f           (nlct), lctlib%crown_length           (nlct), &
             lctlib%bark_par1                (nlct), lctlib%bark_par2              (nlct), &
             lctlib%RCK                      (nlct), lctlib%mort_prob              (nlct), &
             lctlib%fuel_dens                (nlct),                                       &
             lctlib%LitVeg_coef            (nlct,5), lctlib%LeafLit_coef         (nlct,5), &
             lctlib%WoodLit_coef           (nlct,5), lctlib%WoodLitterSize       (nlct))

    IF (p_parallel_io) THEN

       REWIND(unit=lctlib_file_unit) ! Go back to beginning of landcover library file

       ALLOCATE(itmp(nlct)) 

       DO
          READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) EXIT                   ! Finished reading

          ! Look for comment
          pos_comment = SCAN(line,"#")                 ! Start position of comment, 0 if none present
          IF (pos_comment > 1) THEN
             line = TRIM(ADJUSTL(line(1:pos_comment))) ! Disregard everything to the right of "#", ..
             ! .. adjust to the left and trim to the right
          ELSE IF (pos_comment == 1) THEN
             CYCLE                                     ! Disregard whole line
          ELSE
             line = TRIM(ADJUSTL(line))                ! Only disregard blanks to the left and right
          ENDIF
          length=LEN_TRIM(line)
          IF(length== 0) CYCLE                 ! Line is empty
          
          pos = SCAN(line,blank_set)                   ! Position of first blank character
          IF(pos == 0)   CALL finish('init_lctlib',"Wrong syntax in lctlib definitions")
          
          READ(line(1:pos-1),'(A)') key
          key = tolower(key)

          IF(TRIM(key) == "nlct") CYCLE ! nlct already read above

          SELECT CASE (TRIM(key))
          CASE ('lctnumber')
             READ(line(pos:length),*) lctlib%LctNumber(1:nlct)
             exist_LctNumber = .TRUE.
          CASE ('lctname')                                           ! Name of landcover type 
             READ(line(pos:length),*) lctlib%LctName(1:nlct)
             exist_LctName = .TRUE.
          CASE ('landcoverclass')                                    ! Landcover class
             READ(line(pos:length),*) lctlib%LandcoverClass(1:nlct)
             exist_LandcoverClass = .TRUE.
          CASE ('albedosnowvismin')
             READ(line(pos:length),*) lctlib%AlbedoSnowVisMin(1:nlct)
             exist_AlbedoSnowVisMin = .TRUE.
          CASE ('albedosnowvismax')
             READ(line(pos:length),*) lctlib%AlbedoSnowVisMax(1:nlct)
             exist_AlbedoSnowVisMax = .TRUE.
          CASE ('albedosnownirmin')
             READ(line(pos:length),*) lctlib%AlbedoSnowNirMin(1:nlct)
             exist_AlbedoSnowNirMin = .TRUE.
          CASE ('albedosnownirmax')
             READ(line(pos:length),*) lctlib%AlbedoSnowNirMax(1:nlct)
             exist_AlbedoSnowNirMax = .TRUE.
          CASE ('albedosnowmin')
             READ(line(pos:length),*) lctlib%AlbedoSnowMin(1:nlct)
             exist_AlbedoSnowMin = .TRUE.
          CASE ('albedosnowmax')
             READ(line(pos:length),*) lctlib%AlbedoSnowMax(1:nlct)
             exist_AlbedoSnowMax = .TRUE.
          CASE ('albedocanopyvis')
             READ(line(pos:length),*) lctlib%AlbedoCanopyVIS(1:nlct)
             exist_AlbedoCanopyVIS = .TRUE.
          CASE ('albedocanopynir')
             READ(line(pos:length),*) lctlib%AlbedoCanopyNIR(1:nlct)
             exist_AlbedoCanopyNIR = .TRUE.
          CASE ('albedolittervis')
             READ(line(pos:length),*) lctlib%AlbedoLitterVis(1:nlct)
             exist_AlbedoLitterVis = .TRUE.
          CASE ('albedolitternir')
             READ(line(pos:length),*) lctlib%AlbedoLitterNir(1:nlct)
             exist_AlbedoLitterNir = .TRUE.
          CASE ('nitrogenscalingflag')        !Whether nitrogen scaling should be accounted for (.false. input as 0 and .true. as 1)
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%nitrogenScalingFlag(:) = .FALSE.
             ELSEWHERE
                lctlib%nitrogenScalingFlag(:) = .TRUE.
             END WHERE
             exist_nitrogenScalingFlag  = .TRUE.
          CASE ('c4flag')                              ! Photosynthetic pathway (C4: .true.; C3: .false.)
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 1) 
                lctlib%c4flag(:) = .TRUE.
             ELSEWHERE
                lctlib%c4flag(:) = .FALSE.
             END WHERE
             exist_C4flag = .TRUE.
          CASE ('carboxrate')                                         ! Carbox rate
             READ(line(pos:length),*) lctlib%CarboxRate(1:nlct)
             exist_CarboxRate = .TRUE.
          CASE ('etransport')                                         ! E-transport
             READ(line(pos:length),*) lctLib%ETransport(1:nlct)
           exist_ETransport = .TRUE.
          CASE ('vegheight')                                          ! typical vegetation height
             READ(line(pos:length),*) lctLib%VegHeight(1:nlct)
             exist_VegHeight = .TRUE.
          CASE ('vegroughness')                                       ! typical vegetation roughness length
             READ(line(pos:length),*) lctLib%VegRoughness(1:nlct)
             exist_VegRoughness = .TRUE.
          CASE ('minvegroughness')                                    ! typical vegetation roughness length at LAI=0
             READ(line(pos:length),*) lctLib%MinVegRoughness(1:nlct)
             exist_MinVegRoughness = .TRUE.
          CASE ('maxvegroughness')                                    ! typical vegetation roughness length at LAI-->inf
             READ(line(pos:length),*) lctLib%MaxVegRoughness(1:nlct)
             exist_MaxVegRoughness = .TRUE.
          CASE ('phenologytype')                                      ! Phenology type
             READ(line(pos:length),*) lctlib%PhenologyType(1:nlct)
             exist_PhenologyType = .TRUE.
          CASE ('canopyresistancemin')                                ! Minimum canopy resistance
             READ(line(pos:length),*) lctlib%CanopyResistanceMin(1:nlct)
             exist_CanopyResistanceMin = .TRUE.
          CASE ('maxlai')                                              ! Maximum LAI for LoGro-P
             READ(line(pos:length),*) lctlib%MaxLAI(1:nlct)
             exist_MaxLAI = .TRUE.
          CASE ('stemarea')                                            ! Area of stems and branches
             READ(line(pos:length),*) lctlib%StemArea(1:nlct)
             exist_StemArea = .TRUE.
          CASE ('specificleafarea_c')                                  ! Carbon content per leaf area for LoGro-P
             READ(line(pos:length),*) lctlib%specificLeafArea_C(1:nlct)
             exist_specificLeafArea_C = .TRUE.
          CASE ('knorr_tau_w')                                         ! time before leaf shedding
             READ(line(pos:length),*) lctlib%knorr_Tau_w(1:nlct)
             exist_knorr_tau_w = .TRUE.
          CASE ('knorr_t_phi')                                         ! Temperature trigger for leaf growth
             READ(line(pos:length),*) lctlib%knorr_T_phi(1:nlct)
             exist_knorr_t_phi = .TRUE.
          CASE ('knorr_t_r')                                           ! "Spread" (sigma) of T_phi
             READ(line(pos:length),*) lctlib%knorr_T_r(1:nlct)
             exist_knorr_t_r = .TRUE.
          CASE ('knorr_day_c')                                         ! Day-length at leaf shedding
             READ(line(pos:length),*) lctlib%knorr_Day_c(1:nlct)
             exist_knorr_day_c = .TRUE.
          CASE ('knorr_day_r')                                         ! "Spread" (sigma) of Day_c
             READ(line(pos:length),*) lctlib%knorr_Day_r(1:nlct)
             exist_knorr_day_r = .TRUE.
          CASE ('knorr_k_l')                                           ! Inverse of leaf longevity 
             READ(line(pos:length),*) lctlib%knorr_k_l(1:nlct)
             exist_knorr_k_l = .TRUE.
          CASE ('knorr_max_lai')                                       ! maximum LAI
             READ(line(pos:length),*) lctlib%knorr_max_lai(1:nlct)
             exist_knorr_max_lai = .TRUE.
          CASE ('knorr_leaf_growth_rate')                              ! Initial leaf growth rate
             READ(line(pos:length),*) lctlib%knorr_leaf_growth_rate(1:nlct)
             exist_knorr_leaf_growth_rate = .TRUE.
          CASE ('reservec2leafc')                                      ! ratio of C in reserve pool to C in leaves
             READ(line(pos:length),*) lctlib%reserveC2leafC(1:nlct)
             exist_reserveC2leafC = .TRUE.
          CASE ('frac_npp_2_woodpool')                                 ! fraction of NPP directed to wood pool (see cbalance)
             READ(line(pos:length),*) lctlib%frac_npp_2_woodPool(1:nlct)
             exist_frac_npp_2_woodPool = .TRUE.
          CASE ('frac_npp_2_reservepool')                              ! fraction of NPP directed to reserve pool (see cbalance)
             READ(line(pos:length),*) lctlib%frac_npp_2_reservePool(1:nlct)
             exist_frac_npp_2_reservePool = .TRUE.
          CASE ('frac_npp_2_exudates')                                 ! fraction of NPP directed to root exudates (see cbalance)
             READ(line(pos:length),*) lctlib%frac_npp_2_exudates(1:nlct)
             exist_frac_npp_2_exudates = .TRUE.
          CASE ('frac_green_2_herbivory')                              ! fraction of Green directed to herbivory (see cbalance)
             READ(line(pos:length),*) lctlib%frac_green_2_herbivory(1:nlct)
             exist_frac_green_2_herbivory = .TRUE.
          CASE ('tau_cpool_litter_leaf')                           ! Time constant by which leaf litter pool is depreciated  [days]
             READ(line(pos:length),*) lctlib%tau_Cpool_litter_leaf(1:nlct)
             exist_tau_Cpool_litter_leaf = .TRUE.
          CASE ('tau_cpool_litter_wood')                           ! Time constant by which woody litter pool is depreciated  [days]
             READ(line(pos:length),*) lctlib%tau_Cpool_litter_wood(1:nlct)
             exist_tau_Cpool_litter_wood = .TRUE.
          CASE ('tau_cpool_woods')                              ! PFT-specific time_scale for Cpool_woods (and vegetation dynamics)
             READ(line(pos:length),*) lctlib%tau_Cpool_woods(1:nlct)
             ! Conversion from years to days
             lctlib%tau_Cpool_woods(1:nlct) = lctlib%tau_Cpool_woods(1:nlct)*365._dp
             exist_tau_Cpool_woods = .TRUE.
          CASE ('lai_shed_constant')                               ! Time constant by which leaves are constantly shedded [days-1]
             READ(line(pos:length),*) lctlib%LAI_shed_constant(1:nlct)
             exist_LAI_shed_constant = .TRUE.
          CASE ('frac_c_litter_green2atmos')                    ! fraction of loss of fast pool emitted to atmosphere (see cbalance)
             READ(line(pos:length),*) lctlib%frac_C_litter_green2atmos(1:nlct)
             exist_frac_C_litter_green2atmos = .TRUE.
          CASE ('max_c_content_woods')                          ! Maximum carbon content of wood pool (see cbalance)
             READ(line(pos:length),*) lctlib%Max_C_content_woods(1:nlct)
             exist_Max_C_content_woods = .TRUE.
          CASE ('clumpinessfactor')                             ! Factor to calculate clumpiness of vegetation
             READ(line(pos:length),*) lctlib%ClumpinessFactor(1:nlct)
             exist_ClumpinessFactor = .TRUE.
          CASE ('dynamic_pft')                                  ! Indicates those PFTs that shall take part in vegetation dynamics
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%dynamic_PFT(:) = .FALSE.
             ELSEWHERE
                lctlib%dynamic_PFT(:) = .TRUE.
             END WHERE
             exists_dynamic_PFT = .TRUE.
          CASE ('woody_pft')                                    ! Indicates woody PFTs
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%woody_PFT(:) = .FALSE.
             ELSEWHERE
                lctlib%woody_PFT (:) = .TRUE.
            END WHERE
             exists_woody_PFT = .TRUE.
          CASE ('pasture_pft')                                  ! Indicates pasture PFTs
             READ(line(pos:length),*) itmp(1:nlct)
             WHERE(itmp(:) == 0) 
                lctlib%pasture_PFT(:) = .FALSE.
             ELSEWHERE
                lctlib%pasture_PFT (:) = .TRUE.
             END WHERE
             exists_pasture_PFT = .TRUE.
          CASE ('bclimit_min_cold_mmtemp')                      ! PFT-specific minimum coldest monthly mean temperature
             READ(line(pos:length),*) lctlib%bclimit_min_cold_mmtemp(1:nlct)
             exists_bclimit_min_cold_mmtemp = .TRUE.
          CASE ('bclimit_max_cold_mmtemp')                      ! PFT-specific maximum coldest monthly mean temperature
             READ(line(pos:length),*) lctlib%bclimit_max_cold_mmtemp(1:nlct)
             exists_bclimit_max_cold_mmtemp = .TRUE.
          CASE ('bclimit_max_warm_mmtemp')                      ! PFT-specific maximum warmest monthly mean temperature
             READ(line(pos:length),*) lctlib%bclimit_max_warm_mmtemp(1:nlct)
             exists_bclimit_max_warm_mmtemp = .TRUE.
          CASE ('bclimit_min_temprange')                        ! PFT-specific minimum temperature range between warmest and coldest
                                                                !     month (20-year average)
             READ(line(pos:length),*) lctlib%bclimit_min_temprange(1:nlct)
             exists_bclimit_min_temprange = .TRUE.
          CASE ('bclimit_min_gdd')                              ! PFT-specific minimum growing degree days (at or above 5 deg C)
             READ(line(pos:length),*) lctlib%bclimit_min_gdd(1:nlct)
             exists_bclimit_min_gdd = .TRUE.
          CASE ('gdd_base')                                     ! PFT-specific GDD base
             READ(line(pos:length),*) lctlib%gdd_base(1:nlct)
             exists_gdd_base = .TRUE.
          CASE ('upper_tlim')                                   ! PFT-specific upper temperature limit (used for gdd_upper_tlim)
             READ(line(pos:length),*) lctlib%upper_tlim(1:nlct)
             exists_upper_tlim = .TRUE.
          CASE ('frac_wood_2_onsite')                           ! PFT-specific fraction for wood to onSite conversion
             READ(line(pos:length),*) lctlib%frac_wood_2_onSite
             exists_frac_wood_2_onSite = .TRUE.
          CASE ('frac_wood_2_paper')                            ! PFT-specific fraction for wood to paper conversion
             READ(line(pos:length),*) lctlib%frac_wood_2_paper
             exists_frac_wood_2_paper = .TRUE.
          CASE ('frac_wood_2_construction')                     ! PFT-specific fraction for wood to construction conversion
             READ(line(pos:length),*) lctlib%frac_wood_2_construction
             exists_frac_wood_2_construction = .TRUE.
          CASE ('moist_extinction')
            READ(line(pos:length),*) lctlib%moist_extinction    ! PFT-specific moisture of extinction
             exists_moist_extinction = .TRUE.
          CASE ('fuel_dens')
            READ(line(pos:length),*) lctlib%fuel_dens           ! PFT-specific fuel density
             exists_fuel_dens = .TRUE.
          CASE ('flame_length_f')
            READ(line(pos:length),*) lctlib%flame_length_f      ! PFT-specific flame length parameter
             exists_flame_length_f = .TRUE.
          CASE ('crown_length')
            READ(line(pos:length),*) lctlib%crown_length        ! PFT-specific crown length parameter
             exists_crown_length = .TRUE.
          CASE ('bark_par1')
            READ(line(pos:length),*) lctlib%bark_par1           ! PFT-specific bark thickness parameter 1
             exists_bark_par1 = .TRUE.
          CASE ('bark_par2')
            READ(line(pos:length),*) lctlib%bark_par2           ! PFT-specific bark thickness parameter 2
             exists_bark_par2 = .TRUE.
          CASE ('rck')
            READ(line(pos:length),*) lctlib%RCK                 ! PFT-specific resistance to crown damage
             exists_RCK = .TRUE.
          CASE ('mort_prob')
            READ(line(pos:length),*) lctlib%mort_prob           ! PFT-specific parameter for mortality
             exists_mort_prob = .TRUE.
          CASE ('litveg_coef')                                  ! PFT-specific fractions to seperate total litter into yasso pools 
                                                                        !! (NOT NEEDED!)
             READ(line(pos:length),*) lctlib%LitVeg_coef(1:nlct,1)
             DO apu = 2, 5
                READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line  
                READ(line(1:length),*) lctlib%LitVeg_coef(1:nlct,apu)           
             END DO
             DO i = 1, nlct
                ! DSG: this check is important, if values don't add up to 1 carbon is not conserved
                IF (SUM(lctlib%LitVeg_coef(i,:)) .NE. 1.) &
                     CALL finish('init_lctlib','invalid value for LitVeg_coef in '//TRIM(lctlib_file_name))
             END DO
             exists_litveg_coef = .TRUE.
          CASE ('leaflit_coef')                                 ! PFT-specific fractions to seperate leaf litter into yasso pools
             READ(line(pos:length),*) lctlib%LeafLit_coef(1:nlct,1)
             DO apu = 2, 5
                READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line  
                READ(line(1:length),*) lctlib%LeafLit_coef(1:nlct,apu)               
             END DO
             DO i = 1, nlct
                ! DSG: this check is important, if values dont add up to 1 carbon is not conserved
                IF (SUM(lctlib%LeafLit_coef(i,:)) .NE. 1.) &
                     CALL finish('init_lctlib','invalid value for LeafLit_coef in '//TRIM(lctlib_file_name))
             END DO   
             exists_leaflit_coef = .TRUE.
          CASE ('woodlit_coef')                                 ! PFT-specific  fractions to seperate wood litter into yasso pools
             READ(line(pos:length),*) lctlib%WoodLit_coef(1:nlct,1)
             DO apu = 2, 5
                READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line  
                READ(line(1:length),*) lctlib%WoodLit_coef(1:nlct,apu)               
             END DO
             DO i = 1, nlct
                ! DSG: this check is important, if values dont add up to 1 carbon is not conserved
                IF (SUM(lctlib%WoodLit_coef(i,:)) .NE. 1.) &
                     CALL finish('init_lctlib','invalid value for WoodLit_coef_coef in '//TRIM(lctlib_file_name))
             END DO
             exists_woodlit_coef = .TRUE.
          CASE ('woodlittersize')                                ! PFT-specific wood liter size (needed for wood decomposition)
             READ(line(pos:length),*) lctlib%WoodLitterSize(1:nlct)
             exists_woodlittersize_coef = .TRUE.
          CASE default
             ! nothing to do
          END SELECT
       END DO
       DEALLOCATE(itmp)

       IF(.NOT. exist_LctNumber) &
            CALL finish('init_lctlib','No data for LctNumber found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_LctName) &
!            CALL finish('init_lctlib','No data for LctName found in '//TRIM(lctlib_file_name))
            ! Set default
            lctlib%LctName = "undefined"
       IF(.NOT. exist_LandcoverClass) &
            CALL finish('init_lctlib','No data for LandcoverClass found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowVisMin) &
            CALL finish('init_lctlib','No data for AlbedoSnowVisMin found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowVisMax) &
            CALL finish('init_lctlib','No data for AlbedoSnowVisMax found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowNirMin) &
            CALL finish('init_lctlib','No data for AlbedoSnowNirMin found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowNirMax) &
            CALL finish('init_lctlib','No data for AlbedoSnowNirMax found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowMin) &
            CALL finish('init_lctlib','No data for AlbedoSnowMin found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowMax) &
            CALL finish('init_lctlib','No data for AlbedoSnowMax found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoCanopyVIS) &
            CALL finish('init_lctlib','No data for AlbedoCanopyVIS found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoCanopyNIR) &
            CALL finish('init_lctlib','No data for albedoCanopyNIR found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoLitterVis) &
            CALL finish('init_lctlib','No data for AlbedoLitterVis found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoLitterNir) &
            CALL finish('init_lctlib','No data for AlbedoLitterNir found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_nitrogenScalingFlag) &
            CALL finish('init_lctLib','No data for NitrogenScalingFlag found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_C4flag) &
            CALL finish('init_lctLib','No data for C4flag found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_CarboxRate) &
            CALL finish('init_lctLib','No data for CarboxRate found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_ETransport) &
            CALL finish('init_lctLib','No data for ETransport found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_VegHeight) &
            CALL finish('init_lctLib','No data for VegHeight found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_VegRoughness) &
            CALL finish('init_lctLib','No data for VegRoughness found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_MinVegRoughness) &
            CALL finish('init_lctLib','No data for MinVegRoughness found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_MaxVegRoughness) &
            CALL finish('init_lctLib','No data for MaxVegRoughness found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_PhenologyType) &
            CALL finish('init_lctlib','No data for PhenologyType found in '//TRIM(lctlib_file_name))
       IF (.NOT. exist_CanopyResistanceMin) &
            lctlib%CanopyResistanceMin = -1.0_dp
       IF(.NOT. exist_MaxLAI) &
            CALL finish('init_lctlib','No data for MaxLAI found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_StemArea) &
            CALL finish('init_lctlib','No data for StemArea found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_specificLeafArea_C) &
            CALL finish('init_lctlib','No data for specificLeafArea_C found in '//TRIM(lctlib_file_name))

          IF(.NOT. exist_knorr_tau_w) &
               CALL finish('init_lctlib','No data for knorr_Tau_w found in '//TRIM(lctlib_file_name))
          IF(.NOT. exist_knorr_t_phi) &
               CALL finish('init_lctlib','No data for knorr_T_phi found in '//TRIM(lctlib_file_name))
          IF(.NOT. exist_knorr_t_r) &
               CALL finish('init_lctlib','No data for knorr_T_rfound in '//TRIM(lctlib_file_name))
          IF(.NOT. exist_knorr_day_c) &
               CALL finish('init_lctlib','No data for knorr_Day_c found in '//TRIM(lctlib_file_name))
          IF(.NOT. exist_knorr_day_r) &
               CALL finish('init_lctlib','No data for knorr_Day_r found in '//TRIM(lctlib_file_name))
          IF(.NOT. exist_knorr_k_l) &
               CALL finish('init_lctlib','No data for knorr_k_l found in '//TRIM(lctlib_file_name))
          IF(.NOT. exist_knorr_max_lai) &
               CALL finish('init_lctlib','No data for knorr_max_lai found in '//TRIM(lctlib_file_name))
          IF(.NOT. exist_knorr_leaf_growth_rate) &
               CALL finish('init_lctlib','No data for leaf_growth_rate found in '//TRIM(lctlib_file_name))

       IF(.NOT. exist_reserveC2leafC) &
            CALL finish('init_lctlib','No data for reserveC2leafC found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_frac_npp_2_woodPool) &
            CALL finish('init_lctlib','No data for frac_npp_2_woodPool found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_frac_npp_2_reservePool) &
            CALL finish('init_lctlib','No data for frac_npp_2_reservePool found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_frac_npp_2_exudates) &
            CALL finish('init_lctlib','No data for frac_npp_2_exudates found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_frac_green_2_herbivory) &
            CALL finish('init_lctlib','No data for frac_green_2_herbivory found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_tau_Cpool_litter_leaf) &
            CALL finish('init_lctlib','No data for tau_Cpool_litter_leaf found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_tau_Cpool_litter_wood) &
            CALL finish('init_lctlib','No data for tau_Cpool_litter_wood found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_tau_Cpool_woods) &
            CALL finish('init_lctlib','No data for tau_Cpool_woods found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_LAI_shed_constant) &
            CALL finish('init_lctlib','No data for LAI_shed_constant found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_frac_C_litter_green2atmos) &
            CALL finish('init_lctlib','No data for frac_C_litter_green2atmos found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_Max_C_content_woods) &
            CALL finish('init_lctlib','No data for Max_C_content_woods found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_ClumpinessFactor) &
            CALL finish('init_lctlib','No data for ClupinessFactor found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_dynamic_PFT) &
            CALL finish('init_lctlib','No data for DYNAMIC_PFT found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_woody_PFT) &
            CALL finish('init_lctlib','No data for WOODY_PFT found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_pasture_PFT) &
            CALL finish('init_lctlib','No data for PASTURE_PFT found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bclimit_min_cold_mmtemp) &
            CALL finish('init_lctlib','No data for BCLIMIT_MIN_COLD_MMTEMP found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bclimit_max_cold_mmtemp) &
            CALL finish('init_lctlib','No data for BCLIMIT_MAX_COLD_MMTEMP found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bclimit_max_warm_mmtemp) &
            CALL finish('init_lctlib','No data for BCLIMIT_MAX_WARM_MMTEMP found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bclimit_min_temprange ) &
            CALL finish('init_lctlib','No data for BCLIMIT_MIN_TEMPRANGE found in '//TRIM(lctlib_file_name))
       IF(.NOT.exists_bclimit_min_gdd ) &
            CALL finish('init_lctlib','No data for BCLIMIT_MIN_GDD found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_gdd_base) &
            CALL finish('init_lctlib','No data for GDD_BASE found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_upper_tlim) &
            CALL finish('init_lctlib','No data for UPPER_TLIM found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_frac_wood_2_onSite) &
            CALL finish('init_lctlib','No data for FRAC_WOOD_2_ONSITE found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_frac_wood_2_paper) &
            CALL finish('init_lctlib','No data for FRACT_WOOD_2_PAPER found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_frac_wood_2_construction) &
            CALL finish('init_lctlib','No data for FRACT_WOOD_2_CONSTRUCTION found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_moist_extinction) &
            CALL finish('init_lctlib','No data for moisture of extinction found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_fuel_dens) &
            CALL finish('init_lctlib','No data for bulk fuel density found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_flame_length_f) &
            CALL finish('init_lctlib','No data for flame length f found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_crown_length) &
            CALL finish('init_lctlib','No data for crown length found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bark_par1) &
            CALL finish('init_lctlib','No data for bark_par1 found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bark_par2) &
            CALL finish('init_lctlib','No data for bark_par2 found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_RCK) &
            CALL finish('init_lctlib','No data for RCK (resistance to crown damage) found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_mort_prob) &
            CALL finish('init_lctlib','No data for mortality probability found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_woodlittersize_coef) &
            CALL finish('init_lctlib','No data for WOODLITTERSIZE found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_woodlit_coef) &
            CALL finish('init_lctlib','No data for WOODLIT_COEF found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_leaflit_coef) &
            CALL finish('init_lctlib','No data for LEAFLIT_COEF found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_litveg_coef) &
            CALL finish('init_lctlib','No data for LITVEG_COEF found in '//TRIM(lctlib_file_name))

       CLOSE(unit=lctlib_file_unit) !! Closing access to lctlib file

       npft = 0
       !------------------------------------------------------------------------------------------
       ! kalle, 151004
       ! Set carbox rate to  [Mol(CO2)/m^2/s] 
       ! SET e-transport rate to [Mol(CO2)/m^2/s]       
       DO i=1,nlct
         lctLib%ETransport(i)=lctLib%ETransport(i) * 1.e-06_dp
         lctlib%CarboxRate(i)=lctlib%CarboxRate(i) * 1.e-06_dp
       END DO
       !------------------------------------------------------------------------------------------

       !! --- translate landcover class information into landcover class flags
       
       lctlib%BareSoilFlag(1:nlct)   = .FALSE.
       lctlib%GlacierFlag(1:nlct)    = .FALSE.
       lctlib%LakeFlag(1:nlct)       = .FALSE.
       lctlib%ForestFlag(1:nlct)     = .FALSE.
       lctlib%GrassFlag(1:nlct)      = .FALSE.
       lctlib%CropFlag(1:nlct)       = .FALSE.
       lctlib%PastureFlag(1:nlct)    = .FALSE.
       lctlib%NaturalVegFlag(1:nlct) = .FALSE.

       DO i=1,nlct
          SELECT CASE (lctlib%LandcoverClass(i))
          CASE (0)
             lctlib%BareSoilFlag(i) = .TRUE.
          CASE (1)
             lctlib%GlacierFlag(i) = .TRUE.
          CASE (2)
             lctlib%LakeFlag(i) = .TRUE.
          CASE (3)
             lctlib%ForestFlag(i) = .TRUE.
             lctlib%NaturalVegFlag(i) = .TRUE.
          CASE (4)
             lctlib%GrassFlag(i) = .TRUE.
             lctlib%NaturalVegFlag(i) = .TRUE.
          CASE (5)
             lctlib%NaturalVegFlag(i) = .TRUE.
          CASE (6)
             lctlib%CropFlag(i) = .TRUE.
          CASE (7)
             lctlib%PastureFlag(i) = .TRUE.
          CASE default
             WRITE(message_text,*) 'LandcoverClass ', lctlib%LandcoverClass(i),' not allowed in ', TRIM(lctlib_file_name)
             CALL finish('init_lctlib', message_text)
          END SELECT
       END DO

       DO i=1,nlct
          IF(lctlib%naturalVegFlag(i) .OR. lctlib%cropFlag(i) .OR. lctlib%PastureFlag(i)) THEN

             npft = npft+1 ! count PFTs

             ! Check consistency with PhenologyTypes
             IF(lctlib%naturalVegFlag(i) .AND. & 
                  (lctlib%PhenologyType(i) .LE.0 .OR. lctlib%PhenologyType(i) > 4) ) THEN
                CALL finish('init_lctlib',&
                     'Inconsistent entries for one LandcoverClass and PhenologyType in '//TRIM(lctlib_file_name))
             END IF
             IF(lctlib%CropFlag(i) .AND. .NOT. lctlib%PhenologyType(i)==5 ) THEN
                CALL finish('init_lctlib',&
                     'Inconsistent entries for LandcoverClass "crops" and PhenologyType in '//TRIM(lctlib_file_name))
             END IF

          END IF

          !! --- check consistency between landcover classes and specifications for the dynamic vegetation

          IF(lctlib%dynamic_pft(i)  .AND. .NOT. lctlib%NaturalVegFlag(i)) THEN
             CALL finish('init_lctlib',&
                     'Inconsistent entries: Every DYNAMIC_PFT must belong to LandcoverClass "natural" in '//TRIM(lctlib_file_name))
          END IF

          IF(lctlib%woody_pft(i)  .AND. .NOT. lctlib%NaturalVegFlag(i)) THEN
             CALL finish('init_lctlib',&
                     'Inconsistent entries: Every WOODY_PFT must belong to LandcoverClass "natural" in '//TRIM(lctlib_file_name))
          END IF

          IF( lctlib%GrassFlag(i) .NEQV. (lctlib%dynamic_pft(i) .AND. .NOT. lctlib%woody_PFT(i)) ) THEN
             CALL finish('init_lctlib',&
                  'Combination of DYNAMIC_PFT and WOODY_PFT is not consistent with LandcoverClass "grasses" in '//&
                                                                                                     TRIM(lctlib_file_name))
          END IF
       END DO

    END IF !p_parallel_io


    IF(p_parallel) THEN
       CALL p_bcast(lctlib%LctNumber, p_io)
       DO i=1,nlct
          CALL p_bcast(lctlib%LctName(i), p_io)
       ENDDO
       CALL p_bcast(lctlib%LandcoverClass, p_io)
       CALL p_bcast(lctlib%NaturalVegFlag,p_io)
       CALL p_bcast(lctlib%ForestFlag,p_io)
       CALL p_bcast(lctlib%GrassFlag,p_io)
       CALL p_bcast(lctlib%CropFlag,p_io)
       CALL p_bcast(lctlib%PastureFlag,p_io)
       CALL p_bcast(lctlib%LakeFlag, p_io)
       CALL p_bcast(lctlib%GlacierFlag, p_io)
       CALL p_bcast(lctlib%BareSoilFlag, p_io)
       CALL p_bcast(lctlib%AlbedoSnowVisMin, p_io)
       CALL p_bcast(lctlib%AlbedoSnowVisMax, p_io)
       CALL p_bcast(lctlib%AlbedoSnowNirMin, p_io)
       CALL p_bcast(lctlib%AlbedoSnowNirMax, p_io)
       CALL p_bcast(lctlib%AlbedoSnowMin, p_io)
       CALL p_bcast(lctlib%AlbedoSnowMax, p_io)
       CALL p_bcast(lctlib%AlbedoCanopyVIS, p_io)
       CALL p_bcast(lctlib%AlbedoCanopyNIR, p_io)
       CALL p_bcast(lctlib%AlbedoLitterVis, p_io)
       CALL p_bcast(lctlib%AlbedoLitterNir, p_io)
       CALL p_bcast(lctlib%NitrogenScalingFlag,p_io)
       CALL p_bcast(lctlib%C4flag, p_io)
       CALL p_bcast(lctlib%CarboxRate, p_io)
       CALL p_bcast(lctlib%ETransport, p_io)
       CALL p_bcast(lctlib%VegHeight, p_io)
       CALL p_bcast(lctlib%VegRoughness, p_io)
       CALL p_bcast(lctlib%MinVegRoughness, p_io)
       CALL p_bcast(lctlib%MaxVegRoughness, p_io)
       CALL p_bcast(lctlib%PhenologyType, p_io)
       CALL p_bcast(lctlib%CanopyResistanceMin, p_io)
       CALL p_bcast(lctlib%MaxLAI, p_io)
       CALL p_bcast(lctlib%StemArea, p_io)
       CALL p_bcast(lctlib%specificLeafArea_C, p_io)
       CALL p_bcast(lctlib%knorr_Tau_w, p_io)
       CALL p_bcast(lctlib%knorr_T_phi, p_io)
       CALL p_bcast(lctlib%knorr_T_r, p_io)
       CALL p_bcast(lctlib%knorr_Day_c, p_io)
       CALL p_bcast(lctlib%knorr_Day_r, p_io)
       CALL p_bcast(lctlib%knorr_k_l, p_io)
       CALL p_bcast(lctlib%knorr_leaf_growth_rate, p_io)
       CALL p_bcast(lctlib%knorr_max_lai, p_io)
       CALL p_bcast(lctlib%reserveC2leafC, p_io)
       CALL p_bcast(lctlib%frac_npp_2_woodPool, p_io)
       CALL p_bcast(lctlib%frac_npp_2_reservePool, p_io)
       CALL p_bcast(lctlib%frac_npp_2_exudates, p_io)
       CALL p_bcast(lctlib%frac_green_2_herbivory, p_io)
       CALL p_bcast(lctlib%tau_Cpool_litter_leaf, p_io)
       CALL p_bcast(lctlib%tau_Cpool_litter_wood, p_io)
       CALL p_bcast(lctlib%tau_Cpool_woods, p_io)
       CALL p_bcast(lctlib%LAI_shed_constant, p_io)
       CALL p_bcast(lctlib%frac_C_litter_green2atmos, p_io)
       CALL p_bcast(lctlib%Max_C_content_woods, p_io)
       CALL p_bcast(lctlib%ClumpinessFactor, p_io)
       CALL p_bcast(lctlib%dynamic_pft(:), p_io)
       CALL p_bcast(lctlib%woody_pft(:), p_io)
       CALL p_bcast(lctlib%pasture_pft(:), p_io)
       CALL p_bcast(lctlib%bclimit_min_cold_mmtemp(:), p_io)
       CALL p_bcast(lctlib%bclimit_max_cold_mmtemp(:), p_io)
       CALL p_bcast(lctlib%bclimit_max_warm_mmtemp(:), p_io)
       CALL p_bcast(lctlib%bclimit_min_temprange(:), p_io)
       CALL p_bcast(lctlib%bclimit_min_gdd(:), p_io)
       CALL p_bcast(lctlib%gdd_base(:), p_io)
       CALL p_bcast(lctlib%upper_tlim(:), p_io)
       CALL p_bcast(lctlib%frac_wood_2_onSite(:), p_io)
       CALL p_bcast(lctlib%frac_wood_2_paper(:), p_io)
       CALL p_bcast(lctlib%frac_wood_2_construction(:), p_io)
       CALL p_bcast(lctlib%moist_extinction(:), p_io)
       CALL p_bcast(lctlib%fuel_dens(:), p_io)
       CALL p_bcast(lctlib%flame_length_f(:), p_io)
       CALL p_bcast(lctlib%crown_length(:), p_io)
       CALL p_bcast(lctlib%bark_par1(:), p_io)
       CALL p_bcast(lctlib%bark_par2(:), p_io)
       CALL p_bcast(lctlib%RCK(:), p_io)
       CALL p_bcast(lctlib%mort_prob(:), p_io)
       CALL p_bcast(lctlib%LitVeg_coef(:,:), p_io)
       CALL p_bcast(lctlib%LeafLit_coef(:,:), p_io)
       CALL p_bcast(lctlib%WoodLit_coef(:,:), p_io)
       CALL p_bcast(lctlib%WoodLitterSize(:), p_io)
    END IF

    IF(npft>nlct) CALL finish('init_lctlib','Christian did bad programming!')

    lctlib%npft = npft

  END SUBROUTINE init_lctlib
  
END MODULE mo_jsbach_lctlib
