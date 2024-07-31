#!/bin/bash --login
#PBS -N PPE_test_install
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=128
#PBS -j oe
#PBS -A srsei8447 
#PBS -M y.bhatti@sron.nl

#--- comments -------------------------------------------------------
# This script prepares for a PPE (Perturbed Parameter Ensemble)
# experiment. It
# 1) executes 'prepare -r' for each individual ensemble member
# 2) perturbs relevant model parameters for each individual member
# 3) executes 'jobsubm_echam -d' for each individual ensemble member
#
# The names of the perturbed parameters, the experiments and the values of
# the pertubed parameters are stored in a file referred to as $PPEvalues.
# The general settings_*, emi_* and symlinks_* files, which define the
# baseline experiment, should be stored in the directory $PPEdefaults
#
# Basic structure of $PPEvalues:
# 0 <variable name> <variable name> <variable name> ... INSTALL START END
# <experiment name> <float> <float> <float> ... 0 0 0
# <experiment name> <float> <float> <float> ... 0 0 0
# <experiment name> <float> <float> <float> ... 0 0 0
# ...
#
# As a result of its execution, the INSTALL values in $PPEvalues will 
# be set to 1. A file $PPElog will log actions by this script.
#
#
# To execute: qsub PPE_install.sh
#
# Before execution, adapt:
# PBS -l walltime=00:10:00          # estimated time required by this job
# PBS -A n02-NEL013746              # account on ARCHER
# PBS -M schutgens@physics.ox.ac.uk # your email addres
# rundir=".."                       # directory where you store your experimental setup
#
# N.A.J. Schutgens (n.a.j.schutgens@vu.nl)
# 2015/07/07
#--------------------------------------------------------------------
#--- run directory for ECHAM-HAM ------------------------------------
rundir="/home/ybhatti/yusufb/Branches/PPE_Aerosol/my_experiments/snellius/TESTER/"

#--- directories and files used by PPE scripts ----------------------
PPEdir=$rundir'PPE_Experiments/'
PPElog=$rundir'PPE_exp4_log.txt'
PPEtmp=$rundir'PPE_exp4_tmp.txt'
PPEdefaults=$rundir'PPE_Default'
PPEvalues=$rundir'PPE_values.txt'

cd $rundir
echo 'Starting PPE_install script' >$PPElog

while read line
do

  # Read PPEvalues and determine parameter names and values
  elements=( $line )
  nparam=`expr ${#elements[@]} - 4` 
  expid=${elements[0]}
  echo $expid $nparam
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
  echo $installed

  # Check if we found an experiment not yet installed
  if [ "$installed" == 0 ] && [ "$running" == 0 ]; then
      cd $PPEdir
      # Prepare job using JST and modify settings-file for parameter values
      prepare_run.sh -r ${PPEdefaults} $expid
      commands="sed \"-e s/=\s*${param_names[0]}\s*\!/= ${param_values[0]} \!/\""
      for ((iparam=1; iparam<${nparam}; iparam++)); do
        commands="$commands \"-e s/=\s*${param_names[iparam]}\s*\!/= ${param_values[iparam]} \!/\""
      done # iparam
      commands="$commands $expid/settings_$expid >$expid/settings_tmp"
      eval $commands
      mv $expid/settings_tmp $expid/settings_$expid
      echo '    Prepared experiment '$expid >>$PPElog
      echo '    Preparing to edit setting '
      ###########
      ### setting
#      sed -i "\
#      s|walltime=\"*.*|walltime=\"${walltime}\"|g; \
#      s|account=\"*.*|account=\"${account_name}\"|g; \
#      s|email=\"*.*|email=\"${user_email}\"|g; \
#      s|^ncpus=*.*|ncpus=${ncpus}|g; \
#      s|queue_name=\"*.*|queue_name=\"${queue_name}\"|g; \
#
#      " \
#      $expid/settings_${expid}


#      cd ${DIR_JOBKIT}machine_configs/${case_config}
#      sed -i "\
#      s|MEM_MAX_QSUB=\"*.*|MEM_MAX_QSUB=\"${MEM_MAX_QSUB}\"|g; \
#      s|logfile_prefix=\"*.*|logfile_prefix=\"${logfile_prefix}\"|g; \
#      s|ncpus_per_node=\"*.*|ncpus_per_node=\"${ncpus_per_node}\"|g; \
#      s|dir_ncdump=\"*.*|dir_ncdump=\"${DIR_NCDUMP}\"|g; \
#      " \
#      mh-${case_config}.sh
#
#      cd $rundir

      # Dry run job submission
      jobsubm_echam.sh -d $expid/settings_$expid
      echo '    Installed experiment '$expid >>$PPElog
 #     cd $rundir
      installed=1
      #mv $expid $PPEdir

  fi # installed ==0 && running == 0

  if [ "$expid" != 0 ]; then
    echo "$expid ${param_values[*]} $installed $running $ended" >>$PPEtmp
  fi

done < "$PPEvalues"

mv $PPEtmp $PPEvalues 

echo 'Finishing PPE_install script' >>$PPElog
