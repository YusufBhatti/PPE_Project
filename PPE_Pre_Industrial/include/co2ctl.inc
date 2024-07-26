!------------------------------------------------------------------------------
!
! co2ctl.inc
! 
!------------------------------------------------------------------------------
!
! variable name   & description & default \\
!
!------------------------------------------------------------------------------
NAMELIST /co2ctl/ &
  lco2_flxcor,    &! switch on/off flux correction for exact mass balance, .T.
  lco2_mixpbl,    &! switch on/off mix CO2 in PBL, .T.
  lco2_2perc,     &! switch on/off 2% tendency limitation, .F.
  lco2_emis,      &! switch on/off CO2 emission from file, .T.
  lco2_clim,      &! switch on/off climatological CO2 (not transported)
  lco2_scenario    ! switch on/off for CO2 concentrations from external ghg scenario
