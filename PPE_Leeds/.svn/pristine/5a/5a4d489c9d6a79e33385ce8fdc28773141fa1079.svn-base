#! /usr/bin/env ksh
#set -ex
if [[ "$MPI_MODULE" ]]; then
    . $MODULESHOME/init/ksh
    module add $MPI_MODULE
fi
#----------------------------------------------------
# Job file to run echam model on PC in test suite
#----------------------------------------------------
#
#
#----------------------------------------------------
# Variables set by calling script
#----------------------------------------------------
SCR_DIR=$1
KEY=$2
ODIR=$3           # string: base output and run directory
EXP=$4            # string: experiment identifier
MODEL=$5          # string: binary including absolute path
NPROMA=$6         # integer: block length nproma
NPROCA=$7         # integer: nproca processors (nprocb from environment)
LRERUN=$8         # logical: switch on/off restart
SUBM=$9           # binary encoded switches for CO2 submodel, methox submodel
RERUNYEARS=${10}  # years after which rerun is written/performed
LFORCERERUN=${11} # force internal rerun
#----------------------------------------------------
# Model configuration
#----------------------------------------------------
YEAR=1999        # integer: year of simulation
#YEAR=1978        # integer: year of simulation
HRES=31          # integer: horizontal res.
VRES=47          # integer: vertical   res.
#VRES=19          # integer: vertical   res.
ORES=GR30        # ocean resolution
# NPROCB from environment # integer: num. of cpus for lon-dimension
NPROC=$((NPROCA * NPROCB))
#----------------------------------------------------
# create links for input data
#----------------------------------------------------
if [ ! -d $ODIR ]; then
 echo "output directory ${ODIR} for experiment ${EXP} does not exist"
 exit 1
fi
 cd $ODIR
if [ -d $EXP ]; then
 cd $EXP
 if [ -f ${EXP}.err ]; then
  echo "experiment ${EXP} already exists, try to use this one"
  exit
 fi
else
 mkdir $EXP
 cd $EXP
fi
if [ -f ${SCR_DIR}/test_echam6_${KEY}_links.sh ]; then
 ${SCR_DIR}/test_echam6_${KEY}_links.sh ${SCR_DIR} ${HRES} ${VRES} ${ORES} ${YEAR} 1>test.log 2>&1
else
 echo "no file ${SCR_DIR}/test_echam6_${KEY}_links.sh for creating links"
 exit 1
fi
if [ -f ${SCR_DIR}/test_echam6_${KEY}_namelists.sh ]; then
 ${SCR_DIR}/test_echam6_${KEY}_namelists.sh ${SCR_DIR} ${ODIR} ${EXP} ${YEAR} ${NPROCA} ${NPROMA} ${LRERUN} ${SUBM} ${RERUNYEARS} ${LFORCERERUN} 1>>test.log 2>&1
else
 echo "no file ${SCR_DIR}/test_echam6_${KEY}_namelists.sh for creating namelists"
 exit 1
fi
 # In MPI command template, fill in %... place holders, then run command.
 MPIRUN=`sed "s|%n|$NPROC|g;s|%x|$MODEL|g" ${SCR_DIR}/mpirun.dat`
 eval ${MPIRUN} 1> ${EXP}.log 2> ${EXP}.err
