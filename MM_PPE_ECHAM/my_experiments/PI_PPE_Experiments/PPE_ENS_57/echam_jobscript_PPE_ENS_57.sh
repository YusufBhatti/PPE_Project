#!/bin/bash
#SBATCH --job-name="PPE_ENS_57"
#SBATCH --nodes=1
#SBATCH --ntasks=192
#SBATCH --time=14:00:00
#SBATCH --partition genoa
###SBATCH --output=slurm_PPE_ENS_57_%j.txt
###SBATCH --error=slurm_PPE_ENS_57_%j.txt
#SBATCH --output=/home/ybhatti2/prjs1474/Pace_PPE_Output/PI_PPE_Experiments/PPE_ENS_57/PPE_ENS_57.out
#SBATCH --error=/home/ybhatti2/prjs1474/Pace_PPE_Output/PI_PPE_Experiments/PPE_ENS_57/PPE_ENS_57.err
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

cd /home/ybhatti2/prjs1474/Pace_PPE_Output/PI_PPE_Experiments/PPE_ENS_57/

set +e
srun --exclusive --ntasks=192 ./echam6

status_echam="$?"
#
if [[ "$status_echam" -ne "0" ]] && [[ "$status_echam" -ne "127" ]] && [[ "$status_echam" -ne "1" ]] ; then
    echo "ERROR: model run stopped with return value ${status_echam}."
    exit
fi

set -e

#----------------------------------------------------------------------------------------
#--- Launch some p-proc job here:

### Add my CDO command here for P-processing the TAU 3D into 2D (vert int) modes
#cwd=/home/ybhatti2/prjs1474/Pace_PPE_Output/PI_PPE_Experiments/PPE_ENS_57/
## Copy the CDO.sh template to a new script based on jobname
#cp -f /home/ybhatti2/prjs1474/Utils/CDO_Scripts/PPE_ENS_Example.sh /home/ybhatti2/prjs1474/Utils/CDO_Scripts/CDO_PP/PPE_ENS_57.sh
#
## Replace line 6 with the exp_dir variable value
##sed -i '6c\ name='"~/prjs1076/PPE_Output/PPE_Experiments/PPE_ENS_57/" ~/prjs1076/yusufb/Scripts/CDO_PP_PI/PPE_ENS_57.sh
#sed -i "7c\ name=\"/home/ybhatti2/prjs1474/Pace_PPE_Output/PPE_Experiments/PPE_ENS_57/\"" /home/ybhatti2/prjs1474/Utils/CDO_Scripts/CDO_PP/PPE_ENS_57.sh
#
## Replace line 14 with the jobname variable value
#sed -i "6c\ jobname='PPE_ENS_57'" /home/ybhatti2/prjs1474/Utils/CDO_Scripts/CDO_PP/PPE_ENS_57.sh 
#sed -i "6c\mv -f ~/prjs1076/PPE_Output/Pace_PPE_Experiments/PPE_ENS_57 /archive/ybhatti1/Pace_PPE/PPE_Experiments/" ~/Scripts/CDO_PP/PPE_ENS_57.sh

# -----------------------------------------------------

# Define target script once
target="/home/ybhatti2/prjs1474/Utils/CDO_Scripts/CDO_PP_PI/PPE_ENS_57.sh"

# Define paths for reuse
output_dir="/home/ybhatti2/prjs1474/Pace_PPE_Output/PPE_Experiments/PPE_ENS_57/"
archive_dir="~/prjs1076/PPE_Output/Pace_PPE_Experiments/PPE_ENS_57"

# Copy the CDO.sh template to a new script based on jobname
#cp -f /home/ybhatti2/prjs1474/Utils/CDO_Scripts/PPE_ENS_Example.sh /home/ybhatti2/prjs1474/Utils/CDO_Scripts/CDO_PP/PPE_ENS_57.sh
cp -f /home/ybhatti2/prjs1474/Utils/CDO_Scripts/PPE_ENS_Defaults_PI.sh /home/ybhatti2/prjs1474/Utils/CDO_Scripts/CDO_PP_PI/PPE_ENS_57.sh

# Update placeholders inside the copied script in a single sed call
sed -i     -e "6c jobname='PPE_ENS_57'"     -e "7c name=\"/home/ybhatti2/prjs1474/Pace_PPE_Output/PI_PPE_Experiments/PPE_ENS_57\""     "/home/ybhatti2/prjs1474/Utils/CDO_Scripts/CDO_PP_PI/PPE_ENS_57.sh"

set +e   # allow errors
module load 2024                                                                                                                                                                                                                    
module load CDO/2.4.4-gompi-2024a
module load NCO/5.2.9-foss-2024a

. /home/ybhatti2/prjs1474/Utils/CDO_Scripts/CDO_PP_PI/PPE_ENS_57.sh

set -e   # restore strict mode
# -----------------------------------------------------

# Command line prompt tester
#jobname="PPE_ENS_1"; #target="/home/ybhatti2/prjs1474/Utils/CDO_Scripts/CDO_PP/PPE_ENS_57.sh"; #output_dir="/home/ybhatti2/prjs1474/Pace_PPE_Output/PPE_Experiments/PPE_ENS_57/"; #archive_dir="~/prjs1076/PPE_Output/Pace_PPE_Experiments/PPE_ENS_57"; #cp -f /home/ybhatti2/prjs1474/Utils/CDO_Scripts/PPE_ENS_Example.sh ""; #sed -i -e "6c jobname='PPE_ENS_57'" -e "7c name=\"\"" ""


flag_p_proc=false                # flag to launch the p-proc

if $flag_p_proc ; then # this submits a job to be executed on the p-proc machine julier.cscs.ch. 
                       # This job drives the execution of :
                       # /path/to/some/post-processing script.
                       # If this script needs to know more variables than 
                       # 'exp', 'exp_dir' and 'p_proc_dir',
                       # you need to declare them as environment variables before
                       # launching the script

   declare -x exp=PPE_ENS_57
   declare -x exp_dir=/home/ybhatti2/prjs1474/Pace_PPE_Output/PI_PPE_Experiments/PPE_ENS_57/
   declare -x p_proc_dir=/path/to/post_proc_disk/ybhatti2/PPE_ENS_57

#   sbatch -M julier \
#                       --job-name="p-proc" \
#                       --nodes=1 \
#                       --ntasks=1 \
#                       --time=24:00:00 \
#                       --output=/gpfs/work3/0/prjs1474/ybhatti/PPE_Project/PPE_PACE_Project/my_experiments/PI_PPE_Experiments/PPE_ENS_57/slurm_p-proc_%j.txt \
#                       --error=/gpfs/work3/0/prjs1474/ybhatti/PPE_Project/PPE_PACE_Project/my_experiments/PI_PPE_Experiments/PPE_ENS_57/slurm_p-proc_%j.txt \
#                       --account="srs25001" \
#                       --export=ALL \
#                             /path/to/some/post-processing script

   sbatch -M julier \
                       --job-name="p-proc" \
                       --nodes=1               \
                       --ntasks=1 \
                       --time=24:00:00 \
                       --output=/home/ybhatti2/prjs1474/Pace_PPE_Output/PI_PPE_Experiments/PPE_ENS_57/PPE_ENS_57.out \
                       --error=/home/ybhatti2/prjs1474/Pace_PPE_Output/PI_PPE_Experiments/PPE_ENS_57/PPE_ENS_57.err \
                       --account="srs25001" \
                       --export=ALL \
                             /path/to/some/post-processing script


fi

#----------------------------------------------------------------------------------------
#--- Tar the rerun files:

#>>SF uncomment the following line if the path to ncdump is not none none by default
#     and if netcdf is a valid module of your system (else, you may hardcode the path to ncdump)
#module load netcdf   # to get the path to ncdump. Note: still safe if netcdf is already loaded
module load 2024
module load netCDF/4.9.2-gompi-2024a
#<<SF

rerun_echam="restart_PPE_ENS_57_echam.nc"

if [ -e ${rerun_echam} ] ; then
   rerun_date=`ncdump -h ${rerun_echam} | grep ":vdate =" | cut -c12-19`
   rerun_tar=restart_PPE_ENS_57_${rerun_date}.tar

   tar cvf $rerun_tar restart_PPE_ENS_57_[a-z]*
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
      cd /gpfs/work3/0/prjs1474/ybhatti/PPE_Project/PPE_PACE_Project/my_experiments/PI_PPE_Experiments/PPE_ENS_57
      sbatch /gpfs/work3/0/prjs1474/ybhatti/PPE_Project/PPE_PACE_Project/my_experiments/PI_PPE_Experiments/PPE_ENS_57/echam_jobscript_PPE_ENS_57.sh
   fi

