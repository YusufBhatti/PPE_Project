!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_land_boundary

  USE mo_jsbach_grid, ONLY: kstart, kend, nidx
  USE mo_jsbach_lctlib, ONLY: lctlib_type
  USE mo_land_surface, ONLY: land_surface_type
  USE mo_soil, ONLY: soil_type
  USE mo_jsbach_veg, ONLY: vegetation_type
  USE mo_utils, ONLY: average_tiles
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  REAL(dp), PARAMETER :: blending_height = 100._dp     ! blending height [m]
  REAL(dp), PARAMETER :: roughness_snow = 0.001_dp     ! roughness length of surfaces covered by snow [m]
  REAL(dp), PARAMETER :: roughness_bare = 0.005_dp     ! roughness length of bare land [m]
  REAL(dp), PARAMETER :: roughness_lai_saturation = 0.4_dp ! factor in roughness length dependence on LAI
                                                           ! (should between 0.1 and 1)
  !tr REAL(dp), PARAMETER :: roughness_heat_oro = 0.5_dp   ! maximal value of z0,oro contributing to z0,heat
  !tr REAL(dp), PARAMETER :: roughness_momentum_to_heat = 3._dp ! factor, z0,heat is smaller than z0,momentum

CONTAINS

  SUBROUTINE update_land_boundary_up(lctlib, land_surface, soil, vegetation, UseRoughnessLAI, UseRoughnessOro, &
      lai, roughness_oro, roughness_heat, roughness_momentum, &
      surface_temperature, sat_surface_specific_humidity, evapotranspiration, evaporation_pot, &
      sensible_heat_flux_avg, latent_heat_flux_avg, radiative_temp_avg, &
      dry_static_energy_new_avg, soil_moisture_avg, snow_depth_avg, &
      skin_res_avg, veg_height_avg)

    USE mo_soil, ONLY: get_soil_diag

    TYPE(lctlib_type),       INTENT(in) :: lctlib
    TYPE(land_surface_type), INTENT(in) :: land_surface
    TYPE(soil_type),         INTENT(in) :: soil
    TYPE(vegetation_type),   INTENT(in) :: vegetation
    LOGICAL,                 INTENT(in) :: UseRoughnessLAI
    LOGICAL,                 INTENT(in) :: UseRoughnessOro
    REAL(dp),                INTENT(in), DIMENSION(nidx,land_surface%ntiles) :: &
         lai                                   !! LAI
    REAL(dp),                INTENT(in), DIMENSION(nidx) :: &
         roughness_oro                         !! Surface roughness due to orographie (z0,oro)
    REAL(dp),                INTENT(out), DIMENSION(nidx) :: &
         roughness_heat,                    &  !! Surface roughness of heat (z0,heat)
         roughness_momentum,                &  !! Surface roughness of momentum (z0,momentum)
         surface_temperature,               &  !! Surface temperature for grid box
         sat_surface_specific_humidity,     &  !! Saturated surface specific humidity for grid box
         evapotranspiration,                &  !! Evapotranspiration for grid box
         evaporation_pot,                   &  !! Potential evaporation for grid box
         sensible_heat_flux_avg,            &  !! sensible heat flux from surface (W/m2)
         latent_heat_flux_avg,              &  !! latent heat flux from surface (W/m2)
         radiative_temp_avg,                &
         dry_static_energy_new_avg,         &
         soil_moisture_avg, snow_depth_avg, &
         skin_res_avg, veg_height_avg

    REAL(dp) ::  roughness_sum(nidx)
    REAL(dp) ::  roughness_momentum_sum(nidx), roughness_heat_sum(nidx)
    REAL(dp) ::  bare_fract, vegetated_fract
    REAL(dp) ::  roughness_vegetation
    INTEGER  ::  kidx0, kidx1, ntiles, i, j, k

    kidx0 = kstart
    kidx1 = kend

    ntiles = land_surface%ntiles

    ! Surface roughness
    DO i = 1,nidx
      j = kidx0 - 1 + i
      IF (UseRoughnessLAI) THEN
        ! rouhgness length depending on LAI
        roughness_momentum_sum(i) = 0._dp
        roughness_heat_sum(i) = 0._dp
        DO k = 1,ntiles
          bare_fract = (1._dp - land_surface%veg_ratio_max(j)) * land_surface%cover_fract(j,k)
          vegetated_fract = land_surface%veg_ratio_max(j) * land_surface%cover_fract(j,k)
          ! contributions from bare land
          roughness_momentum_sum(i) = roughness_momentum_sum(i) + bare_fract / & ! bare land
                                      (LOG(blending_height / roughness_bare) ** 2)
! TR lower z0 for snow covered surfaces should be discussed again, if a multi-layer snow model is implemented.
          roughness_heat_sum(i) = roughness_heat_sum(i) + bare_fract * soil%snow_fract(j,k) / & ! snow covered bare land
                                  (LOG(blending_height / roughness_snow) ** 2)
          roughness_heat_sum(i) = roughness_heat_sum(i) + bare_fract * (1._dp - soil%snow_fract(j,k)) / & ! bare land without snow
                                  (LOG(blending_height / roughness_bare) ** 2)
          ! contributions from vegetated land
          roughness_vegetation = lctlib%MinVegRoughness(land_surface%cover_type(j,k)) + &
                                 (lctlib%MaxVegRoughness(land_surface%cover_type(j,k)) - &
                                 lctlib%MinVegRoughness(land_surface%cover_type(j,k))) * &
                                 TANH(lai(i,k) * roughness_lai_saturation)
          roughness_momentum_sum(i) = roughness_momentum_sum(i) + vegetated_fract / & ! vegetated land
                                      (LOG(blending_height / roughness_vegetation) ** 2)
          IF (land_surface%is_forest(j,k)) THEN
            roughness_heat_sum(i) = roughness_heat_sum(i) + vegetated_fract / & ! vegetated land
                                    (LOG(blending_height / roughness_vegetation) ** 2)
          ELSE
            roughness_heat_sum(i) = roughness_heat_sum(i) + vegetated_fract * soil%snow_fract(j,k) / & ! snow covered vegetated land
                                    (LOG(blending_height / roughness_snow) ** 2)
            roughness_heat_sum(i) = roughness_heat_sum(i) + vegetated_fract * (1._dp - soil%snow_fract(j,k)) / & ! veg. without snow
                                    (LOG(blending_height / roughness_vegetation) ** 2)
          END IF
        END DO
        ! TR  z0 of heat should be much smaller than z0 of momentum (bluff body effect). However for echam6.3 z0 of heat and z0 of
        !     momentum are kept equal.
        !     Change the comment status of all lines with "tr" to reduce z0 of heat.
        IF (UseRoughnessOro) THEN
          roughness_momentum(i) = SQRT((blending_height * EXP(-1._dp / SQRT(roughness_momentum_sum(i)))) ** 2 + &
                                  roughness_oro(i) ** 2)
          !tr roughness_heat(i) = SQRT((blending_height * EXP(-1._dp / SQRT(roughness_heat_sum(i)))) ** 2 + &
          !tr                         MIN(roughness_oro(i),roughness_heat_oro) ** 2) / roughness_momentum_to_heat
        ELSE
          roughness_momentum(i) = blending_height * EXP(-1._dp / SQRT(roughness_momentum_sum(i)))
          !tr roughness_heat(i) = blending_height * EXP(-1._dp / SQRT(roughness_heat_sum(i))) / roughness_momentum_to_heat
        END IF
        roughness_heat(i) = blending_height * EXP(-1._dp / SQRT(roughness_heat_sum(i)))  !tr
      ELSE ! UseRoughnessLAI = .false.
        ! roughness length without dependence on LAI
        roughness_sum(i) = 0._dp
        DO k = 1,ntiles
          roughness_sum(i) = roughness_sum(i) + ((1._dp - land_surface%veg_ratio_max(j)) *      & ! bare land
                             land_surface%cover_fract(j,k) /                                    &
                             (LOG(blending_height / roughness_bare) ** 2))
          roughness_sum(i) = roughness_sum(i) + (land_surface%veg_ratio_max(j) *                & ! vegetated land
                             land_surface%cover_fract(j,k) /                                    &
                             (LOG(blending_height / lctlib%VegRoughness(land_surface%cover_type(j,k))) ** 2))
        END DO
        IF (UseRoughnessOro) THEN
          roughness_momentum(i) = SQRT((blending_height * EXP(-1._dp / SQRT(roughness_sum(i)))) ** 2 + roughness_oro(i) ** 2)
        ELSE
          roughness_momentum(i) = blending_height * EXP(-1._dp / SQRT(roughness_sum(i)))
        END IF
        roughness_sum(i) = 0._dp
        DO k = 1,ntiles
          roughness_sum(i) = roughness_sum(i) + (soil%snow_fract(j,k) *                         & ! snow covered land
                             land_surface%cover_fract(j,k) /                                    &
                             (LOG(blending_height / roughness_snow) ** 2))
          roughness_sum(i) = roughness_sum(i) + ((1._dp - soil%snow_fract(j,k)) *               & ! not snow covered land
                             land_surface%cover_fract(j,k) /                                    &
                            (LOG(blending_height / MIN(1._dp,roughness_momentum(i))) ** 2))
        END DO
        roughness_heat(i) = blending_height * EXP(-1._dp / SQRT(roughness_sum(i)))
      END IF ! UseRoughnessLAI

    END DO

    surface_temperature(1:nidx) = get_soil_diag(nidx, 'surface_temperature')
    sat_surface_specific_humidity(1:nidx) = get_soil_diag(nidx, 'sat_surface_specific_humidity')
    radiative_temp_avg(1:nidx) = get_soil_diag(nidx, 'surface_radiative_temperature')

    CALL average_tiles(soil%evapotranspiration(kidx0:kidx1,1:ntiles), land_surface%is_present(kidx0:kidx1,1:ntiles) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,1:ntiles), evapotranspiration(:))
    CALL average_tiles(soil%evaporation_pot(kidx0:kidx1,1:ntiles), land_surface%is_present(kidx0:kidx1,1:ntiles) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,1:ntiles), evaporation_pot(:))
    CALL average_tiles(soil%latent_heat_flux(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),latent_heat_flux_avg(:))
    CALL average_tiles(soil%sensible_heat_flux(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),sensible_heat_flux_avg(:))
    CALL average_tiles(soil%dry_static_energy_new(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),dry_static_energy_new_avg(:))
    CALL average_tiles(soil%moisture(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),soil_moisture_avg(:))
    CALL average_tiles(soil%snow(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),snow_depth_avg(:))
    CALL average_tiles(soil%skin_reservoir(kidx0:kidx1,:), land_surface%is_present(kidx0:kidx1,:) &
         .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:),skin_res_avg(:))
    CALL average_tiles(vegetation%veg_height(kidx0:kidx1,:), &
         land_surface%is_present(kidx0:kidx1,:) .AND. .NOT. land_surface%is_lake(kidx0:kidx1,:), &
         land_surface%cover_fract(kidx0:kidx1,:), veg_height_avg(:))

  END SUBROUTINE update_land_boundary_up

END MODULE mo_land_boundary
