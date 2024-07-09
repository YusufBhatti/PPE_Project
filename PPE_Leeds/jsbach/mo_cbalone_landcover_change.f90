!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cbalone_landcover_change

  USE mo_kind,          ONLY: dp

  IMPLICIT none

  PUBLIC :: initLandCoverChange

  PRIVATE

CONTAINS

  SUBROUTINE initLandCoverChange(grid, ntiles, &
                                 with_landuse_transitions, with_nitrogen, &
                                 lcc, landuse_transitions, lcc_scheme)

    USE mo_jsbach_grid,               ONLY: grid_type
    USE mo_cbal_landcover_change,     ONLY: landuse_transitions_type, landcover_change_type
    USE mo_jsbach,                    ONLY: debug_Cconservation

    TYPE(grid_type),                  INTENT(in)    :: grid
    INTEGER,                          INTENT(in)    :: ntiles
    LOGICAL,                          INTENT(in)    :: with_landuse_transitions
    LOGICAL,                          INTENT(in)    :: with_nitrogen
    TYPE(landcover_change_type),      INTENT(inout) :: lcc
    TYPE(landuse_transitions_type),   INTENT(inout) :: landuse_transitions
    INTEGER,                          INTENT(in)    :: lcc_scheme

    ALLOCATE (lcc%LCC_flux_box_C2atmos(grid%nland))
    lcc%LCC_flux_box_C2atmos(:) = 0.0_dp
    IF(with_nitrogen) THEN
       ALLOCATE (lcc%LCC_flux_box_N2atmos(grid%nland))
       lcc%LCC_flux_box_N2atmos(:) = 0.0_dp
    ENDIF

    IF (lcc_scheme==1) THEN  ! standard JSBACH landcover scheme 
       allocate (lcc%LCC_flux_box_C2litterGreenPools(grid%nland), &
                 lcc%LCC_flux_box_C2litterWoodPool(grid%nland))
       lcc%LCC_flux_box_C2litterWoodPool  (:) = 0.0_dp
       lcc%LCC_flux_box_C2litterGreenPools(:) = 0.0_dp
       
       IF(with_nitrogen) THEN
          allocate (lcc%LCC_flux_box_N2litterGreenPools(grid%nland), &
                    lcc%LCC_flux_box_N2litterWoodPool(grid%nland),   &
                    lcc%LCC_flux_box_N2SMINNpool(grid%nland))
          lcc%LCC_flux_box_N2litterGreenPools(:) = 0.0_dp
          lcc%LCC_flux_box_N2litterWoodPool  (:) = 0.0_dp
          lcc%LCC_flux_box_N2SMINNpool       (:) = 0.0_dp
          
          ALLOCATE (lcc%LCC_testNconserv(grid%nland))
          
       END IF

    ENDIF
 
    IF (lcc_scheme==2) THEN 

      ALLOCATE(lcc%C_2_onSite(grid%nLand),                &
               lcc%C_2_paper(grid%nLand),                 &
               lcc%C_2_construction(grid%nLand),          &
               lcc%C_onSite_2_atmos(grid%nLand),          &
               lcc%C_paper_2_atmos(grid%nLand),           &
               lcc%C_construction_2_atmos(grid%nLand),    &
               lcc%boxC_2_onSite(grid%nLand),             &
               lcc%boxC_2_paper(grid%nLand),              &
               lcc%boxC_2_construction(grid%nLand),       &
               lcc%boxC_onSite_2_atmos(grid%nLand),       &
               lcc%boxC_paper_2_atmos(grid%nLand),        &
               lcc%boxC_construction_2_atmos(grid%nLand))

      lcc%boxC_2_onSite                    (:) = 0._dp
      lcc%boxC_2_paper                     (:) = 0._dp
      lcc%boxC_2_construction              (:) = 0._dp
      lcc%boxC_onSite_2_atmos              (:) = 0._dp
      lcc%boxC_paper_2_atmos               (:) = 0._dp
      lcc%boxC_construction_2_atmos        (:) = 0._dp

      IF (with_landuse_transitions) THEN
         ALLOCATE(lcc%C_2_paper_harv(grid%nLand),                 &
                  lcc%C_2_construction_harv(grid%nLand),          &
                  lcc%C_paper_harv_2_atmos(grid%nLand),           &
                  lcc%C_construction_harv_2_atmos(grid%nLand),    &
                  lcc%boxC_2_paper_harv(grid%nLand),              &
                  lcc%boxC_2_construction_harv(grid%nLand),       &
                  lcc%boxC_paper_harv_2_atmos(grid%nLand),        &
                  lcc%boxC_construction_harv_2_atmos(grid%nLand))

         lcc%boxC_2_paper_harv             (:) = 0._dp
         lcc%boxC_2_construction_harv      (:) = 0._dp
         lcc%boxC_paper_harv_2_atmos       (:) = 0._dp
         lcc%boxC_construction_harv_2_atmos(:) = 0._dp
      ENDIF

    ENDIF

    IF (debug_Cconservation) THEN
       ALLOCATE (lcc%LCC_testCconserv(grid%nland))
    END IF
    ALLOCATE (lcc%LCC_coverFract_target(grid%nland, ntiles))
    
    IF (with_landuse_transitions) THEN
       ALLOCATE (landuse_transitions%TransMtrx_CROP_2_PAST(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_PAST_2_CROP(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_NATL_2_PAST(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_PAST_2_NATL(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_CROP_2_NATL(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_NATL_2_CROP(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_FRST_2_PAST(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_PAST_2_FRST(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_GRAS_2_PAST(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_PAST_2_GRAS(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_FRST_2_CROP(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_CROP_2_FRST(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_GRAS_2_CROP(grid%nland))
       ALLOCATE (landuse_transitions%TransMtrx_CROP_2_GRAS(grid%nland))
       ALLOCATE (landuse_transitions%Grass_coverFract_lastYear(grid%nland))
       ALLOCATE (landuse_transitions%NatWood_coverFract_lastYear(grid%nland))
       ALLOCATE (landuse_transitions%Pasture_coverFract_lastYear(grid%nland))
       ALLOCATE (landuse_transitions%Crop_coverFract_lastYear(grid%nland))

       IF (debug_Cconservation) THEN
          ALLOCATE (landuse_transitions%Test_NATL_2_PAST(grid%nland))
          ALLOCATE (landuse_transitions%Test_PAST_2_NATL(grid%nland))
          ALLOCATE (landuse_transitions%Test_CROP_2_NATL(grid%nland))
          ALLOCATE (landuse_transitions%Test_NATL_2_CROP(grid%nland))
          ALLOCATE (landuse_transitions%Test_CROP_2_PAST(grid%nland))
          ALLOCATE (landuse_transitions%Test_PAST_2_CROP(grid%nland))
          ALLOCATE (landuse_transitions%Test_FRST_2_PAST(grid%nland))
          ALLOCATE (landuse_transitions%Test_PAST_2_FRST(grid%nland))
          ALLOCATE (landuse_transitions%Test_GRAS_2_PAST(grid%nland))
          ALLOCATE (landuse_transitions%Test_PAST_2_GRAS(grid%nland))
          ALLOCATE (landuse_transitions%Test_GRAS_2_CROP(grid%nland))
          ALLOCATE (landuse_transitions%Test_CROP_2_GRAS(grid%nland))
          ALLOCATE (landuse_transitions%Test_FRST_2_CROP(grid%nland))
          ALLOCATE (landuse_transitions%Test_CROP_2_FRST(grid%nland))
       END IF

       ALLOCATE (landuse_transitions%CROP_2_NATL_ignored(grid%nland))
       ALLOCATE (landuse_transitions%PAST_2_NATL_ignored(grid%nland))
       ALLOCATE (landuse_transitions%Box_harvest(grid%nland))

       ALLOCATE (landuse_transitions%Box_flux_harvest(grid%nland))        
       ALLOCATE (landuse_transitions%Box_flux_harvest_2atmos(grid%nland)) 
       ALLOCATE (landuse_transitions%C2atmos_LUtrans(grid%nland))
       ALLOCATE (landuse_transitions%C2atmos_harvest(grid%nland))
       ALLOCATE (landuse_transitions%N2atmos_LUtrans(grid%nland))
       ALLOCATE (landuse_transitions%N2atmos_harvest(grid%nland))
    ELSE
       ALLOCATE (lcc%C2atmos(grid%nland))
       lcc%C2atmos = 0._dp
    END IF

  END SUBROUTINE initLandCoverChange

END MODULE mo_cbalone_landcover_change
