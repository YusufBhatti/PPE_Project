#!/bin/bash

#---------------------------------------------------
#
# Sylvaine Ferrachat 2014-03 and Ulrike Proske 2023
#
# Script to produce 3D CCN input for the CCN climatology
# submodel, going over 10 years
#---------------------------------------------------


infile=$1 # In the case of this script, the infile is created in between with cdo. How the average over multiple years is coded below it makes sense to call this infile 2003-2012_ham.nc
outfile_rootname=$2 # 2003-2012 (makes sense like this because of hardcoding, as above)
refdate=$3  # e.g. 2000-01-01
timestep=$4 # 1month or 1day
exp=$5 # Name of the experiment
exp_out=$6 # Folder under which the experiment should be saved, e.g. ${exp}_10years

raw_dir="/scratch/snx3000/uproske/${exp}"
pproc_dir="/scratch/snx3000/uproske/Pproc/${exp}/CCNclim_correcttimestamp"
output_dir="/project/s1144/uproske/aeroclim_climatologies/${exp_out}"

#-- modified CCN:
infile_mod="${infile%??????}"activ.nc

# ------- Create 10 year files ---------

# First need to get the yearly netcdf file (mean over 10 years)
# Take the mean over all years
for i in {1..9}; do
	cdo ensmean ${raw_dir}/${exp}_200[23456789]0${i}.01_activ.nc ${raw_dir}/${exp}_201[012]0${i}.01_activ.nc ${pproc_dir}/${exp}_2003-20120${i}_activ.nc
	cdo ensmean ${raw_dir}/${exp}_200[23456789]0${i}.01_ham.nc ${raw_dir}/${exp}_201[012]0${i}.01_ham.nc ${pproc_dir}/${exp}_2003-20120${i}_ham.nc
done
for i in {0..2}; do
	cdo ensmean ${raw_dir}/${exp}_200[23456789]1${i}.01_activ.nc ${raw_dir}/${exp}_201[012]1${i}.01_activ.nc ${pproc_dir}/${exp}_2003-20121${i}_activ.nc
	cdo ensmean ${raw_dir}/${exp}_200[23456789]1${i}.01_ham.nc ${raw_dir}/${exp}_201[012]1${i}.01_ham.nc ${pproc_dir}/${exp}_2003-20121${i}_ham.nc
done
cdo copy ${pproc_dir}/${exp}_2003-2012??_activ.nc ${infile_mod}
cdo copy ${pproc_dir}/${exp}_2003-2012??_ham.nc ${infile}

echo "done!"

# Above way is smarter since it avoids creating very large files to then merge over
#for i in {3..9}; do
#	cdo copy ${raw_dir}/${exp}_200${i}??.01_ham.nc ${pproc_dir}/${exp}_200${i}_ham.nc
#	cdo copy ${raw_dir}/${exp}_200${i}??.01_activ.nc ${pproc_dir}/${exp}_200${i}_activ.nc
#done
#for i in {0..2}; do
#	cdo copy ${raw_dir}/${exp}_201${i}??.01_ham.nc ${pproc_dir}/${exp}_201${i}_ham.nc
#	cdo copy ${raw_dir}/${exp}_201${i}??.01_activ.nc ${pproc_dir}/${exp}_201${i}_activ.nc
#done
#cdo ensmean ${pproc_dir}/${exp}_20??_ham.nc ${infile}
#cdo ensmean ${pproc_dir}/${exp}_20??_activ.nc ${infile_mod}

# ----- Below here is original CCNclim creation -----------

\rm -r tmp_*.nc

#-- Set list of operators to apply for time and lev fixes:
operator_cdo_time_fix="-settunits,days -settaxis,$refdate,00:00:00,$timestep -setcalendar,standard -setreftime,$refdate,00:00:00"

operator_nco_lev_fix="-O -d .lev_2,lev -v .lev_2,lev"

#-- CCN:

outfile_ccn=${outfile_rootname}_CCN.nc

cdo chname,CCN_STRAT,CCN_an_strat_mo,CCN_CONV,CCN_an_conv_mo,FRACDUSOL,CCN_dust_mo -selvar,CCN_STRAT,CCN_CONV,FRACDUSOL $infile tmp_1.nc

cdo -chname,CCN_an_strat_mo,CCN_pi_strat_mo -setrtoc,-inf,inf,0 -selvar,CCN_an_strat_mo tmp_1.nc tmp_pi_strat.nc
cdo -chname,CCN_an_conv_mo,CCN_pi_conv_mo -setrtoc,-inf,inf,0 -selvar,CCN_an_conv_mo tmp_1.nc tmp_pi_conv.nc

cdo -chname,CCN_an_strat_mo,CCN_coarse_mo -setrtoc,-inf,inf,0 -selvar,CCN_an_strat_mo tmp_1.nc tmp_coarse.nc

cdo merge tmp_1.nc tmp_pi_strat.nc tmp_pi_conv.nc tmp_coarse.nc tmp_3d_full.nc

eval cdo $operator_cdo_time_fix tmp_3d_full.nc tmp_2.nc 
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_ccn

cp $outfile_ccn ${output_dir}/${outfile_ccn}

#-- IN:

outfile_in=${outfile_rootname}_IN.nc

#--- artificially create dummy rwetki (=1. It's needed for technical reasons),
#    artificially create dummy fracbcinsol (=0. It's needed for technical reasons),
#    and fetches the other vars

cat > tmp_exprf.txt << EOF
FRACDUSOL=FRACDUSOL;
FRACDUAI=FRACDUAI;
FRACDUCI=FRACDUCI;
FRACBCSOL=FRACBCSOL;
FRACBCINSOL=FRACBCSOL;
RWETKI=RWETAI;
RWETAI=RWETAI;
RWETCI=RWETCI;
RWETAS=RWETAS;
EOF
# Are these typos above? RWETKI, FRACBCINSOL? No, these are renamed because we don't have them for the new ccnclim

cdo -exprf,tmp_exprf.txt $infile tmp_1.nc

eval cdo $operator_cdo_time_fix tmp_1.nc tmp_2.nc
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_in

cp $outfile_in ${output_dir}/${outfile_in}

#-- modified CCN:

echo $infile_mod

if [ -e $infile_mod ] ; then 

# CCN only where ACT
echo 'ACTCCN'
outfile_ccn_mod=${outfile_rootname}_CCN_ACTCCN.nc
cdo chname,CCN_CONV,CCN_an_conv_mo,FRACDUSOL,CCN_dust_mo -selvar,CCN_CONV,FRACDUSOL $infile tmp_1b.nc
cdo selvar,CCN_STRAT_ACTCCN,CCN_STRAT_ACTCCN_CTR $infile_mod tmp.nc 
cdo -setmisstoc,0 -expr,'CCN_an_strat_mo= CCN_STRAT_ACTCCN / CCN_STRAT_ACTCCN_CTR' tmp.nc tmp_1a.nc
# below needed for model to read in data
cdo -setunit,'m-3' tmp_1a.nc tmp_2a.nc
ncatted -a long_name,CCN_an_strat_mo,a,c,'Number of activated aerosols in stratiform clouds (acc.) - modified ACTCCN' tmp_2a.nc tmp_3a.nc
cdo merge tmp_3a.nc tmp_1b.nc tmp_pi_strat.nc tmp_pi_conv.nc tmp_coarse.nc tmp_3d_full_mod.nc
eval cdo $operator_cdo_time_fix tmp_3d_full_mod.nc tmp_2.nc 
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_ccn_mod

cp $outfile_ccn_mod ${output_dir}/${outfile_ccn_mod}

cdo ensmax $outfile_ccn_mod $outfile_ccn ${output_dir}/${outfile_rootname}_CCN_maxofact.nc

fi
