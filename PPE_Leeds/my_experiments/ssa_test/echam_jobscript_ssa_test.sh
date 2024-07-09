#!/bin/bash
#SBATCH --job-name="ssa_test"
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=12:00:00
#SBATCH --partition rome
###SBATCH --output=slurm_ssa_test_%j.txt
###SBATCH --error=slurm_ssa_test_%j.txt
#SBATCH --output=/projects/0/prjs0937/yusufb/ECHAM/OUTPUT//snellius/DEBUG/ssa_test/ssa_test.out
#SBATCH --error=/projects/0/prjs0937/yusufb/ECHAM/OUTPUT//snellius/DEBUG/ssa_test/ssa_test.err
#SBATCH --account="srsei8447"
#!/bin/bash -l

set -eu 

#-----------------------------------------------------------------------------
# Sylvaine Ferrachat 2011 (after Doris Folini)
#
# Runscript for echam on CSCS Cray XE6 machine (rosa), as of 01.2012 (slurm)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Env Variables:

unset mc
ulimit -s unlimited

export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=536870912
export OMP_NUM_THREADS=1
export MPICH_MAX_SHORT_MSG_SIZE=
export MPICH_PTL_UNEX_EVENTS=
export MPICH_UNEX_BUFFER_SIZE=

#----------------------------------------------------------------------------------------
#--- Start the run:

cd /projects/0/prjs0937/yusufb/ECHAM/OUTPUT//snellius/DEBUG/ssa_test/

set +e
mpiexec -n 40 ./echam6

status_echam="$?"
#
if [[ "$status_echam" -ne "0" ]] && [[ "$status_echam" -ne "127" ]] && [[ "$status_echam" -ne "1" ]] ; then
    echo "ERROR: model run stopped with return value ${status_echam}."
    exit
fi

set -e

#----------------------------------------------------------------------------------------
#--- Launch some p-proc job here:

flag_p_proc=false                # flag to launch the p-proc

if $flag_p_proc ; then # this submits a job to be executed on the p-proc machine julier.cscs.ch. 
                       # This job drives the execution of :
                       # /path/to/some/post-processing script.
                       # If this script needs to know more variables than 
                       # 'exp', 'exp_dir' and 'p_proc_dir',
                       # you need to declare them as environment variables before
                       # launching the script

   declare -x exp=ssa_test
   declare -x exp_dir=/projects/0/prjs0937/yusufb/ECHAM/OUTPUT//snellius/DEBUG/ssa_test/
   declare -x p_proc_dir=/path/to/post_proc_disk/ybhatti/ssa_test

#   sbatch -M julier \
#                       --job-name="p-proc" \
#                       --nodes=1 \
#                       --ntasks=1 \
#                       --time=XX:XX:XX \
#                       --output=/gpfs/work3/0/prjs0937/yusufb/Branches/PPE_Leeds/my_experiments/ssa_test/slurm_p-proc_%j.txt \
#                       --error=/gpfs/work3/0/prjs0937/yusufb/Branches/PPE_Leeds/my_experiments/ssa_test/slurm_p-proc_%j.txt \
#                       --account="srsei8447" \
#                       --export=ALL \
#                             /path/to/some/post-processing script

   sbatch -M julier \
                       --job-name="p-proc" \
                       --nodes=1               \
                       --ntasks=1 \
                       --time=XX:XX:XX \
                       --output=/projects/0/prjs0937/yusufb/ECHAM/OUTPUT//snellius/DEBUG/ssa_test/ssa_test.out \
                       --error=/projects/0/prjs0937/yusufb/ECHAM/OUTPUT//snellius/DEBUG/ssa_test/ssa_test.err \
                       --account="srsei8447" \
                       --export=ALL \
                             /path/to/some/post-processing script


fi

#----------------------------------------------------------------------------------------
#--- Tar the rerun files:

#>>SF uncomment the following line if the path to ncdump is not none none by default
#     and if netcdf is a valid module of your system (else, you may hardcode the path to ncdump)
#module load netcdf   # to get the path to ncdump. Note: still safe if netcdf is already loaded 
#<<SF

rerun_echam="restart_ssa_test_echam.nc"

if [ -e ${rerun_echam} ] ; then
   rerun_date=`ncdump -h ${rerun_echam} | grep ":vdate =" | cut -c12-19`
   rerun_tar=restart_ssa_test_${rerun_date}.tar

   tar cvf $rerun_tar restart_ssa_test_[a-z]*
else
   echo "No rerun file was produced: Don't do any job chaining and stop execution here"
   exit # no need to prepare a potential namelist file for a next submission if there's no rerun file
fi

#----------------------------------------------------------------------------------------
#--- Adjust the rerun flag in the namelist for any subsequent job:

\mv namelist.echam namelist.echam.bak
cat namelist.echam.bak | sed -e 's/^ *[lL][rR][eE][sS][uU][mM][eE].*$/lresume=.true./' > namelist.echam

#----------------------------------------------------------------------------------------
#--- Submit the next job:

   if [[ "$status_echam" -eq "0" ]] ; then
      cd /gpfs/work3/0/prjs0937/yusufb/Branches/PPE_Leeds/my_experiments/ssa_test
      sbatch /gpfs/work3/0/prjs0937/yusufb/Branches/PPE_Leeds/my_experiments/ssa_test/echam_jobscript_ssa_test.sh
   fi

