! submdiagctl
!
! controls diagnostic output for generic submodel variables and streams
! - vphysc : physical ECHAM variables not contained in standard ECHAM output
! - wetdep : variables used in the calculation of wet deposition and extra diagnostics
! - drydep : variables used in the calculation of dry deposition
! 
!---------------------------------------------------------------------------
NAMELIST /submdiagctl/ &
vphysc_lpost,       & ! true for output of vphysc stream, false otherwise
vphysc_tinterval,   & ! output interval
vphyscnam,          & ! names of variables you like to have in output
                      !      special names: ALL and DEFAULT
                      !!!
wetdep_lpost,       & ! true for output of wetdep stream, false otherwise
wetdep_tinterval,   & ! output interval
wetdepnam,          & ! names of variables (diagnostic quantities) you like to have in output
                      !      special names: ALL, DETAIL and DEFAULT
wetdep_gastrac,     & ! names of gas-phase tracers to be included in wetdep diagnostic output
                      ! aerosol tracers will always be output -- potentially aggregated (see below)
wetdep_keytype,     & ! aggregation level of output: BYTRACER, BYSPECIES, BYMODE
                      !!!
drydep_lpost,       & ! true for output of drydep stream, false otherwise
drydep_tinterval,   & ! output interval
drydepnam,          & ! names of variables (diagnostic quantities) you like to have in output
                      !      special names: ALL (= DETAIL) and DEFAULT
drydep_gastrac,     & ! names of gas-phase tracers to be included in drydep diagnostic output
                      ! aerosol tracers will always be output -- potentially aggregated (see below)
drydep_keytype,     & ! aggregation level of output: BYTRACER, BYSPECIES, BYMODE
drydep_ldetail,     & ! detailed output of dry deposition diagnostics   !++mgs 20140519
drydep_trac_detail, & ! tracer (only one at a run!) for detailed drydep output (drydep_ldetail=.true.)
                      !!!
sedi_lpost,         & ! true for output of sedi stream, false otherwise
sedi_tinterval,     & ! output interval
sedinam,            & ! names of variables (diagnostic quantities) you like to have in output
                      !      special names: ALL, DETAIL and DEFAULT     ### needed ?? ###
sedi_keytype,       & ! aggregation level of output: BYTRACER, BYSPECIES, BYMODE
emi_lpost,          & ! true for output of emi stream, false otherwise
emi_lpost_sector,   & ! true for output of detailed (per sector) emission diagnostics
emi_tinterval,      & ! output interval
eminam,             & ! names of variables (diagnostic quantities) you like to have in output
                      !      special names: ALL (= DETAIL) and DEFAULT
emi_gastrac,        & ! names of gas-phase tracers to be included in emi diagnostic output
                      ! aerosol tracers will always be output -- potentially aggregated (see below)
emi_keytype           ! aggregation level of output: BYTRACER, BYSPECIES, BYMODE
                      !!!

