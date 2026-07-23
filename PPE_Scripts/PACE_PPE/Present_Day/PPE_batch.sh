#!/bin/bash --login
#SBATCH --job-name=PPE_run
#SBATCH -t 14:00:00
#SBATCH -p genoa
#SBATCH --nodes=65
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --account=srs25001

#--- user settings ---
nexp_maxrunning=50    # number of concurrent jobs
nexp_maxran=100      # max total jobs (big number = all)

cwd=${PWD}/
PPEdir='/home/ybhatti2/prjs1474/PPE_Project/PPE_PACE_Project/my_experiments/PPE_Experiments/'
PPElog=$cwd'PPE_log.txt'
PPEtmp=$cwd'PPE_tmp.txt'
PPEvalues=$cwd'PPE_values.txt'

cd $PPEdir
rm -f $PPElog
echo 'Starting PPE_batch script' >$PPElog

nexp_ran=0
batch_finished=0

while [ "$nexp_ran" -lt "$nexp_maxran" ] && [ "$batch_finished" == 0 ]; do
  nexp_running=0
  nexp_left=0
  rm -f $PPEtmp

  # take a snapshot copy to avoid truncating
  cp "$PPEvalues" "${PPEvalues}.work"

  while read line; do
    elements=( $line )
    nparam=$(( ${#elements[@]} - 4 ))
    expid=${elements[0]}

    if [ "$expid" == 0 ]; then
      # header line
      echo "$line" > $PPEtmp
      continue
    fi

    installed=${elements[nparam+1]}
    running=${elements[nparam+2]}
    ended=${elements[nparam+3]}

    if [ "$installed" == 1 ] && [ "$running" == 0 ]; then
      nexp_left=$((nexp_left + 1))

      if [ "$nexp_running" -lt "$nexp_maxrunning" ]; then
        $expid/echam_jobscript_$expid.sh &
        echo "    Started experiment $expid" >>$PPElog
        nexp_running=$((nexp_running + 1))
        running=1
      fi
    fi

    echo "$expid ${elements[@]:1:$nparam} $installed $running $ended" >> $PPEtmp
  done < "${PPEvalues}.work"

  rm "${PPEvalues}.work"
  mv $PPEtmp $PPEvalues

  wait   # wait for all launched jobs
  nexp_ran=$((nexp_ran + nexp_running))

  if [ "$nexp_running" -eq 0 ]; then
    batch_finished=1
  fi
done

echo 'Finishing PPE_batch script' >>$PPElog
cd $cwd

