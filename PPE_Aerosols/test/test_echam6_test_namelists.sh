SCR_DIR=$1
ODIR=$2
EXP=$3
YEAR=$4
NPROCA=$5
NPROMA=$6
LRERUN=$7
SUBM=$8
RERUNYEARS=$9
LFORCERERUN=${10}
if [ $LFORCERERUN == '.false.' ]; then
NO_CYCLES=2
else 
NO_CYCLES=1
fi
#read out_filetype from file
OUTFILETYPE=`cat ${SCR_DIR}/outfiletype.dat`
#evaluate submodel switch (l_methox,l_co2):
if [ ${#SUBM} -lt 2 ]; then
 SUBM=0${SUBM}
fi
if [ $SUBM == '00' ]; then
 CO2=.false.      # logical: switch on/off CO2 model
 METHOX=.false.
fi
if [ $SUBM == '10' ]; then
 CO2=.true.
 METHOX=.false.
fi
if [ $SUBM == '01' ]; then
 CO2=.false.
 METHOX=.true.
fi
if [ $SUBM == '11' ]; then
 CO2=.true.
 METHOX=.true.
fi
cd ${ODIR}/${EXP}
#---------------------------------------------------
# ECHAM6 namelist
#---------------------------------------------------
cat > namelist.echam << EOF
&parctl
  nproca       = ${NPROCA}
  nprocb       = ${NPROCB}
/
&runctl
  ldebugev     = .true.
  ldebugio     = .false.
!  default_output = .false.
  out_datapath = "${ODIR}/${EXP}/"
  out_expname  = "${EXP}"
  out_filetype = ${OUTFILETYPE}                 ! 1 - GRIB1, 2 - netCDF
  lresume      = $LRERUN
  lamip        = .true.
!  labort       = .false.
  dt_start     = ${YEAR},12,31,22,0,0
!  dt_start     = ${YEAR},01,01,12,0,0
  putdata      = 1, 'steps', 'first', 0
  no_steps     = 13
  no_cycles    = ${NO_CYCLES}
  putcheckpoint= 32,'years','last',0
  putrerun     = ${RERUNYEARS},'years','last',0
! putrerun     = 1,'steps','last',0
  lforcererun  = ${LFORCERERUN}
  nproma       = ${NPROMA} 
  lmidatm      = .true. 
  ltdiag       = .true.
!  ldebugs      = .true.
  lmeltpond    = .true. 
/
&tdiagctl
puttdiag = 1, 'steps', 'first', 0
tdiagnam = 'all'
/
&debugsctl
  putdebug_stream = 1, 'steps', 'first', 0
/
&gwsctl
lrmscon_lat    = .true.
!    lat_rmscon_lo = 2.5
!    lat_rmscon_hi = 20.
!    rmscon_lo     = 1.15
!    rmscon_hi     = 0.9
/
&radctl
  lradforcing = T,T
  isolrad      = 1
  ighg         = 1
  iaero        = 5
  io3          = 4
  ico2         = 4
  in2o         = 4
  ich4         = 4
  icfc         = 4
/
&mvstreamctl
!target = 'echam'
interval = 3, 'steps', 'first', 0
filetag = 'echamm'
source = 'g3b','sp','gl'
meannam = ''
sqrmeannam = 'svo','q','tsurf'
/
&submodelctl
  lmethox      = ${METHOX}
  lco2         = ${CO2}
!  lco2         = .true.
/
!&co2ctl
!lco2_emis       = .false.
!lco2_flxcor     = .false.
!lco2_2perc      = .false.
!/
&submdiagctl
  vphyscnam='all'
/
&cospctl
 locosp = .TRUE.
/
&cfdiagctl
 locfdiag = .true.
/
&stationctl
lostation=.TRUE.
/
&mvstreamctl
m_stream_name='tdiag'
/
EOF
cat > tdiag.nml <<EOF
&mvctl
putmean = 2,'steps','first',0
meannam = 'all'
/
EOF
#---------------------------------------------------
# JS-BACH namelists 
#---------------------------------------------------
cat > namelist.jsbach <<EOF
&jsbach_ctl
  standalone    = .false.
  ntiles        = 11                   ! --- number of tiles ---
  ! --- options to activate the different jsbach modules ---
  use_bethy     = .true. 
  use_phenology = .true. 
  use_albedo    = .true. 
  use_dynveg    = .false.
  use_disturbance = .true.
  with_nitrogen = .false.
  use_roughness_lai = .true.
  use_roughness_oro = .false.
  ! --- output options ---
  file_type     = $OUTFILETYPE
  file_ztype    = 0
  lpost_echam   = .false.
  veg_at_1200 = .true.
  debug         = .false.
/
&albedo_ctl
  use_albedocanopy = .false.
/
&cbalance_ctl
  read_cpools = .false.
/
&dynveg_ctl
  read_fpc = .false.
  dynveg_feedback = .false. 
/
&climbuf_ctl
  init_running_means = .false.
    read_climbuf = .false.
/
&disturbance_ctl
    fire_name = ''
/
&windbreak_jsbach_ctl
    wind_damage_scale = 0.01
/
&soil_ctl
    nsoil = 5
/
EOF
