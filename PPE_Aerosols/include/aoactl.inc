! aoactl.inc
!
!-------------------------------------------------------------------------
NAMELIST /aoactl/ &
  dt_start_emi_winter_1,     &! start date of emission (age+winter spectrum)
  dt_start_emi_winter_2,     &! start date of emission (winter spectrum)
  dt_start_emi_winter_3,     &! start date of emission (winter spectrum)
  dt_start_emi_winter_4,     &! start date of emission (winter spectrum)
  dt_start_emi_summer_1,     &! start date of emission (summer spectrum)
  dt_start_emi_summer_2,     &! start date of emission (summer spectrum)
  dt_start_emi_summer_3,     &! start date of emission (summer spectrum)
  dt_start_emi_summer_4,     &! start date of emission (summer spectrum)
  conc_increase,             &! increase of tracer concentration in mmr
  emission_plev,             &! pressure height of emissions (hPa)
  emission_lat_s,            &! lower threshold latitude for emissions
  emission_lat_n              ! upper threshold latitude for emissions
!-------------------------------------------------------------------------
