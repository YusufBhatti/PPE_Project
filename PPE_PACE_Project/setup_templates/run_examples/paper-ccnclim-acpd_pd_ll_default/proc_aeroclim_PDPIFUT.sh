#!/bin/bash

#------------------------------------------------------------------------------------------------------
# Ulrike Proske 2023-08
#
# Main script to process echam raw output to produce data as input for AEROclim
#
# Usage, e.g.: ~/code/echam_aerosol_scripts/AEROclim/proc_aeroclim_PDPIFUT.sh aeroclim_ll_default_05 200301 201212 aeroclim_ll_default_10years
# while in: /scratch/snx3000/uproske/Pproc/aeroclim_ll_default_10years/Sandbox
#------------------------------------------------------------------------------------------------------

set -e

#-- Arguments
exp1="$1"         # exp name as in the echam raw output
yymm_start="$2"  # YYYYMM for start of the exp
yymm_stop="$3"   # YYYYMM for end of exp
exp2="$4"         # exp name for the processed output

raw_dir="/scratch/snx3000/uproske/${exp1}"
pproc_dir="/scratch/snx3000/uproske/Pproc/${exp2}"
output_dir="/project/s1144/uproske/aeroclim_climatologies/${exp2}"
ext=".nc"

cyyyymm=`date --utc "+%Y%m" --date "${yymm_start}01"`
cyyyymm_stop=`date --utc "+%Y%m" --date "${yymm_stop}01"`
cmm=`date --utc "+%Y%m" --date "${cyyyymm}01"`

while [[ $cyyyymm -le ${cyyyymm_stop} ]] ; do

	echo -e ${cyyyymm}
	echo -e ${cyyyymm_stop}
	echo -e ${cmm}

	eval cdo "-settunits,days -settaxis,${cyyyymm:0:4}-${cyyyymm:4:2}-01,00:00:00,1month -setcalendar,standard -setreftime,2002-10-01,00:00:00" ${raw_dir}/${exp1}_${cyyyymm:0:4}${cyyyymm:4:2}.01_tracerm.nc ${pproc_dir}/${exp2}_${cyyyymm:0:4}${cyyyymm:4:2}.01_tracerm.nc
	#eval cdo "-settunits,days -settaxis,${cyyyymm:0:4}-${cyyyymm:5:2}-01,00:00:00,1month -setcalendar,standard -setreftime,2002-10-01,00:00:00" ${raw_dir}/${exp_1}_${cyyyymm:0:4}${cyyyymm:5:2}.01_tracerm.nc ${exp2}/${exp2}_${cyyyymm:0:4}${cyyyymm:5:2}.01_tracerm.nc

        cyyyymm=`date --utc "+%Y%m" --date "${cyyyymm}01 + 1 month"`
        cmm=`date --utc "+%Y%m" --date "${cyyyymm}01"`
done

# Take the mean over all years
for i in {1..9}; do
	cdo ensmean ${pproc_dir}/${exp2}_20??0${i}.01_tracerm.nc ${output_dir}/${exp2}_2003-20120${i}_tracerm.nc
done
for i in {0..2}; do
	cdo ensmean ${pproc_dir}/${exp2}_20??1${i}.01_tracerm.nc ${output_dir}/${exp2}_2003-20121${i}_tracerm.nc
done

