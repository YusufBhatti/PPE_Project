#!/bin/bash
#SBATCH --job-name="PPE_Leeds"
#SBATCH --nodes=1
#SBATCH --ntasks=60
#SBATCH --time=12:00:00
#SBATCH --partition normal
###SBATCH --output=slurm_PPE_Leeds_%j.txt
###SBATCH --error=slurm_PPE_Leeds_%j.txt
#SBATCH --output=/home/ybhatti/yusufb/ECHAM/OUTPUT/PPE_Leeds/PPE_Leeds.out
#SBATCH --error=/home/ybhatti/yusufb/ECHAM/OUTPUT/PPE_Leeds/PPE_Leeds.err
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

cd /home/ybhatti/yusufb/ECHAM/OUTPUT/PPE_Leeds/

set +e
mpiexec -n 60 ./echam6

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

   declare -x exp=PPE_Leeds
   declare -x exp_dir=/home/ybhatti/yusufb/ECHAM/OUTPUT/PPE_Leeds/
   declare -x p_proc_dir=/path/to/post_proc_disk/ybhatti/PPE_Leeds

#   sbatch -M julier \
#                       --job-name="p-proc" \
#                       --nodes=1 \
#                       --ntasks=1 \
#                       --time=XX:XX:XX \
#                       --output=/gpfs/work3/0/prjs0937/yusufb/Branches/PPE_Leeds/my_experiments/PPE_Leeds/slurm_p-proc_%j.txt \
#                       --error=/gpfs/work3/0/prjs0937/yusufb/Branches/PPE_Leeds/my_experiments/PPE_Leeds/slurm_p-proc_%j.txt \
#                       --account="srsei8447" \
#                       --export=ALL \
#                             /path/to/some/post-processing script

   sbatch -M julier \
                       --job-name="p-proc" \
                       --nodes=1               \
                       --ntasks=1 \
                       --time=XX:XX:XX \
                       --output=/home/ybhatti/yusufb/ECHAM/OUTPUT/PPE_Leeds/PPE_Leeds.out \
                       --error=/home/ybhatti/yusufb/ECHAM/OUTPUT/PPE_Leeds/PPE_Leeds.err \
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

rerun_echam="restart_PPE_Leeds_echam.nc"

if [ -e ${rerun_echam} ] ; then
   rerun_date=`ncdump -h ${rerun_echam} | grep ":vdate =" | cut -c12-19`
   rerun_tar=restart_PPE_Leeds_${rerun_date}.tar

   tar cvf $rerun_tar restart_PPE_Leeds_[a-z]*
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
      cd /gpfs/work3/0/prjs0937/yusufb/Branches/PPE_Leeds/my_experiments/PPE_Leeds
      sbatch /gpfs/work3/0/prjs0937/yusufb/Branches/PPE_Leeds/my_experiments/PPE_Leeds/echam_jobscript_PPE_Leeds.sh
   fi

