#!/bin/bash --login
#PBS -N PPE_test_run
#PBS -q genoa
#PBS -l walltime=120:00:00
#PBS -l select=11
#PBS -j oe
#PBS -A srsei9480 
#PBS -M y.bhatti@sron.nl
#PBS -m ae

#--- comments -------------------------------------------------------
# *******
# For real runs -l should be 8*2 and walltime should be number of experiments * 5hrs / 2
# E.g. 15 hours and 16 nodes for 6 experiments
#       For 13 experiments (with plenty of /work space): -l select 8*4 and walltime (5*14) / 4 = 18
# *******

# This script executes a PPE (Perturbed Parameter Ensemble) experiment.
# It will execute ${nexp_maxrunning} experiments at the same time, until
# ${nexp_maxran} experiments have been executed. These experiments are
# defined in a file called ${PPEvalues}.
#
# The names of the perturbed parameters, the experiments and the values of
# the pertubed parameters are stored in a file referred to as ${PPEvalues}.
#
# Basic structure of ${PPEvalues}:
# 0 <variable name> <variable name> <variable name> ... INSTALL START END
# <experiment name> <float> <float> <float> ... 1 0 0
# <experiment name> <float> <float> <float> ... 1 0 0
# <experiment name> <float> <float> <float> ... 1 0 0
# ...
#
# This scripts will only execute an experiment if INSTALL=1 
# (due to another script PPE_install.sh) and if the experiment has not been executed before
# (START=0). 
#
# As a result of its execution, the START values in ${PPEvalues} will 
# be set to 1. A file ${PPElog} will log actions by this script.
#
# Note that this script is written for parallel queue, using
# multiple aprun commands to execute the PPE as efficiently as possible.
# PBS -l select=N should be chosen so that N is equal to the maximum number of 
# experiments running concurrently times the number of nodes required for
# an individual experiment
#
#
# To execute: qsub -q "queue name" PPE_batch.sh
#
# Before execution, adapt:
# PBS -l walltime=00:10:00            # estimated time required by this job
# PBS -l select=4                     # nexp_maxrunning * nr of nodes
# PBS -A
# PBS -M 
# nexp_maxrunning=2                   # Maximum number of experiments run concurrently
# nexp_maxran=8                       # Maximum number of experiments run by this batch job
# rundir=".."                         # directory where you store your experimental setup
#
# N.A.J. Schutgens (schutgens@physics.ox.ac.uk)
# 2015/07/07
#--------------------------------------------------------------------

#--- User definitions -----------------------------------------------
nexp_maxrunning=10  # Maximum number of experiments run concurrently
nexp_maxran=40      # Maximum number of experiments run by this batch job

#--- run directory for ECHAM-HAM ------------------------------------
#rundir="/home/ybhatti/yusufb/Branches/PPE_Leeds/my_experiments/snellius/"
cwd=${PWD}/

#--- directories and files used by PPE scripts ----------------------
PPEdir='/home/ybhatti/yusufb/Branches/PPE_Aerosols/my_experiments/PI_PPE_Experiments/'
PPElog=$cwd'PPE_log.txt'
PPEtmp=$cwd'PPE_tmp.txt'
PPEdefaults=$cwd'PPE_Default'
PPEvalues=$cwd'PPE_values.txt'

#--- set maximum number of qsub runs ---
#max_qsub_runs=$((nexp_maxran / nexp_maxrunning - 1))
max_qsub_runs = 3 # I will want this script to run 4 times (0-3) = 40 * 4 simulations)
#--- maintain qsub run count ---
if [ ! -f "${cwd}/current_qsub_run.txt" ]; then
  echo 0 > "${cwd}/current_qsub_run.txt"
fi
current_qsub_run=$(cat "${cwd}/current_qsub_run.txt")
echo "Current Resub is $current_qsub_run out of $max_qsub_runs"

#echo 'Activating python environment'
#source activate master

cd $PPEdir

rm -f $PPElog
echo 'Starting PPE_batch script' >$PPElog

#--- start outer loop
nexp_ran=0 # Number of experiments that have already ran
batch_finished=0 # Flag: have we done everything we wanted to do?
while [ "$nexp_ran" -lt "$nexp_maxran" ] && [ "$batch_finished" == 0 ]; do
 
  #--- start inner loop: find nexp_maxrunning experiments and start them
  nexp_running=0 # Number of experiments currently running
  nexp_left=0 # Number of experiments that we still need to do
  while read line
  do

    # Read PPEvalues and determine parameter names and values
    elements=( $line )
    nparam=`expr ${#elements[@]} - 4` 
    expid=${elements[0]}
    if [ "$expid" == 0 ]; then
      for ((iparam=0; iparam<${nparam}; iparam++)); do
        param_names[iparam]=${elements[iparam+1]}
      done # iparam
      echo '  Found '$nparam' parameters: '${param_names[*]} >>$PPElog
      rm -f $PPEtmp
      echo "$expid ${param_names[*]} INSTALL START END" >$PPEtmp
    else
      for ((iparam=0; iparam<${nparam}; iparam++)); do
        param_values[iparam]=${elements[iparam+1]}
      done # iparam
    fi # expid == 0
    installed=${elements[nparam+1]}
    running=${elements[nparam+2]}
    ended=${elements[nparam+3]}

    # Check if we found an experiment has not run
    if [ "$installed" == 1 ] && [ "$running" == 0 ]; then

      nexp_left=`expr $nexp_left + 1`

     # Start experiment if not too many experiments already running
      if [ "$nexp_running" -lt "$nexp_maxrunning" ]; then
        # Launch job
#	srun --nodes=2 $expid/echam_jobscript_$expid.sh &
        $expid/echam_jobscript_$expid.sh &
        echo '    Started  experiment '$expid >>$PPElog
	echo $expid'/echam_jobscript_'$expid'.sh' >>$PPElog
        nexp_running=`expr $nexp_running + 1`
        running=1
	echo $nexp_running >>$PPElog

      fi # nexp_running < nexp_maxrunning
    #  echo 'THIS IS FOR FI #1' >>$PPElog
    fi # running == 0
   # echo 'THIS IS FOR FI #2' >>$PPElog

    if [ "$expid" != 0 ]; then
      echo "$expid ${param_values[*]} $installed $running $ended" >>$PPEtmp
      #echo 'THIS IS CHANGING RUNNING INSTALLED NUMBERS' >>$PPElog

    fi

  done < "$PPEvalues" # inner loop
  echo 'inner loop finished' >>$PPElog

  # Some admin at the end of each sequence of job submissions 
  if [ "$nexp_running" == "$nexp_left" ]; then
    batch_finished=0
  fi # nexp_running == nexp_left
  nexp_left=`expr $nexp_left - $nexp_running` 
  mv $PPEtmp $PPEvalues
  echo '  Started '$nexp_running' experiments, '$nexp_left' experiments left' >>$PPElog

  wait # wait until all recently started experiments are finished 
  nexp_ran=`expr $nexp_ran + $nexp_running`

done # nexp_ran > nexp_maxran (outer loop)

echo 'Finishing PPE_batch script' >>$PPElog
#. $cwd/../Scripts_for_config/mv_tracking_files.sh
cd $cwd


#--- run the batch job again if conditions are met ---
if [[ "$current_qsub_run" -lt "$max_qsub_runs" ]]; then
  current_qsub_run=$((current_qsub_run + 1))
  echo $current_qsub_run > "${cwd}/current_qsub_run.txt"
  qsub PPE_batch.sh
  echo "Submitted PPE_batch.sh for the $current_qsub_run time out of $max_qsub_runs times" >>$PPElog
else
  echo "Max qsub runs of $max_qsub_runs reached. Exiting." >>$PPElog
  rm ${cwd}/current_qsub_run.txt
fi
