#!/bin/bash
#SBATCH --job-name="obscrip"
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=10:00:00
#SBATCH --partition genoa
###SBATCH --output=slurm_obscrip_%j.txt
###SBATCH --error=slurm_obscrip_%j.txt
#SBATCH --output=/home/ybhatti2/prjs1474/Pace_PPE_Output/obscrip/obscrip.out
#SBATCH --error=/home/ybhatti2/prjs1474/Pace_PPE_Output/obscrip/obscrip.err
#SBATCH --account="srs25001"
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

cd /home/ybhatti2/prjs1474/Pace_PPE_Output/obscrip/

set +e
srun --exclusive --ntasks=128 ./echam6

status_echam="$?"
#
if [[ "$status_echam" -ne "0" ]] && [[ "$status_echam" -ne "127" ]] && [[ "$status_echam" -ne "1" ]] ; then
    echo "ERROR: model run stopped with return value ${status_echam}."
    exit
fi

set -e

#----------------------------------------------------------------------------------------
#--- Launch some p-proc job here:

## Add my CDO command here for P-processing the TAU 3D into 2D (vert int) modes
cwd=/home/ybhatti2/prjs1474/Pace_PPE_Output/obscrip/
DIR=/home/ybhatti1/prjs1076/ybhatti/yusufb
# Copy the CDO.sh template to a new script based on jobname
#cp -f ~/Scripts/CDO.sh ~/Scripts/CDO_PP_PI/obscrip.sh
cp ~/prjs1076/yusufb/Scripts/CDO.sh ~/prjs1076/yusufb/Scripts/CDO_PP/obscrip.sh
# Replace line 6 with the exp_dir variable value
#sed -i '6c\ name='"~/prjs1076/PPE_Output/PPE_Experiments/obscrip/" ~/prjs1076/yusufb/Scripts/CDO_PP_PI/obscrip.sh
sed -i "6c\ name=\"~/prjs1076/PPE_Output/PPE_Experiments/obscrip/\"" ~/prjs1076/yusufb/Scripts/CDO_PP/obscrip.sh
# Replace line 14 with the jobname variable value

sed -i "5c\ jobname='obscrip'" ~/Scripts/CDO_PI/obscrip.sh

sed -i "63c\mv -f ~/prjs1076/PPE_Output/Pace_PPE_Experiments/obscrip /archive/ybhatti1/Pace_PPE/PPE_Experiments/" ~/Scripts/CDO_PP/obscrip.sh

# Source the new script (if you want to run it in the current shell)
#. ~/Scripts/CDO_PP_PI/obscrip.sh

#mv -f ~/prjs1076/PPE_Output/Corrected_PI_PPE_Experiments/obscrip /archive/ybhatti1/Corrected_PI_PPE_Experiments/

flag_p_proc=false                # flag to launch the p-proc

if $flag_p_proc ; then # this submits a job to be executed on the p-proc machine julier.cscs.ch. 
                       # This job drives the execution of :
                       # /path/to/some/post-processing script.
                       # If this script needs to know more variables than 
                       # 'exp', 'exp_dir' and 'p_proc_dir',
                       # you need to declare them as environment variables before
                       # launching the script

   declare -x exp=obscrip
   declare -x exp_dir=/home/ybhatti2/prjs1474/Pace_PPE_Output/obscrip/
   declare -x p_proc_dir=/path/to/post_proc_disk/ybhatti2/obscrip

#   sbatch -M julier \
#                       --job-name="p-proc" \
#                       --nodes=1 \
#                       --ntasks=1 \
#                       --time=24:00:00 \
#                       --output=/gpfs/work3/0/prjs1474/ybhatti/PPE_Project/PPE_PACE_Project/my_experiments/obscrip/slurm_p-proc_%j.txt \
#                       --error=/gpfs/work3/0/prjs1474/ybhatti/PPE_Project/PPE_PACE_Project/my_experiments/obscrip/slurm_p-proc_%j.txt \
#                       --account="srs25001" \
#                       --export=ALL \
#                             /path/to/some/post-processing script

   sbatch -M julier \
                       --job-name="p-proc" \
                       --nodes=1               \
                       --ntasks=1 \
                       --time=24:00:00 \
                       --output=/home/ybhatti2/prjs1474/Pace_PPE_Output/obscrip/obscrip.out \
                       --error=/home/ybhatti2/prjs1474/Pace_PPE_Output/obscrip/obscrip.err \
                       --account="srs25001" \
                       --export=ALL \
                             /path/to/some/post-processing script


fi

#----------------------------------------------------------------------------------------
#--- Tar the rerun files:

#>>SF uncomment the following line if the path to ncdump is not none none by default
#     and if netcdf is a valid module of your system (else, you may hardcode the path to ncdump)
#module load netcdf   # to get the path to ncdump. Note: still safe if netcdf is already loaded 
#<<SF

rerun_echam="restart_obscrip_echam.nc"

if [ -e ${rerun_echam} ] ; then
   rerun_date=`/home/ybhatti/mpich-4/netcdf4-4.9.0/bin/ncdump -h ${rerun_echam} | grep ":vdate =" | cut -c12-19`
   rerun_tar=restart_obscrip_${rerun_date}.tar

   tar cvf $rerun_tar restart_obscrip_[a-z]*
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
      cd /gpfs/work3/0/prjs1474/ybhatti/PPE_Project/PPE_PACE_Project/my_experiments/obscrip
      sbatch /gpfs/work3/0/prjs1474/ybhatti/PPE_Project/PPE_PACE_Project/my_experiments/obscrip/echam_jobscript_obscrip.sh
   fi

