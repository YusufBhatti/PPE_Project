#!/bin/bash

#---------------------------------------------------
#
# Sylvaine Ferrachat 2014-03 and Ulrike Proske 2023
#
# Script to produce 3D CCN input for the CCN climatology
# submodel
# Different variants of the climatology are produced here, 
# only some of which were used in sensitivity simulations in the paper
#---------------------------------------------------


infile=$1
outfile_rootname=$2
refdate=$3  # e.g. 2000-01-01
timestep=$4 # 1month or 1day

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

#-- modified CCN:

infile_mod="${infile%??????}"activ.nc

echo $infile_mod

if [ -e $infile_mod ] ; then 

# CDNCMIN
echo 'CDNCMIN'
outfile_ccn_mod=${outfile_rootname}_CCN_CDNCMIN.nc
cdo chname,CCN_CONV,CCN_an_conv_mo,FRACDUSOL,CCN_dust_mo -selvar,CCN_CONV,FRACDUSOL $infile tmp_1b.nc
cdo selvar,CCN_STRAT_CDNCMIN,CCN_STRAT_CDNCMIN_CTR $infile_mod tmp.nc
cdo -setmisstoc,0 -expr,'CCN_an_strat_mo=CCN_STRAT_CDNCMIN/CCN_STRAT_CDNCMIN_CTR' tmp.nc tmp_1a.nc
# below needed for model to read in data
cdo -setunit,'m-3' tmp_1a.nc tmp_2a.nc
ncatted -a long_name,CCN_an_strat_mo,a,c,'Number of activated aerosols in stratiform clouds (acc.) - modified CDNCMIN' tmp_2a.nc tmp_3a.nc
cdo merge tmp_3a.nc tmp_1b.nc tmp_pi_strat.nc tmp_pi_conv.nc tmp_coarse.nc tmp_3d_full_mod.nc
eval cdo $operator_cdo_time_fix tmp_3d_full_mod.nc tmp_2.nc 
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_ccn_mod

# NOSMALLER40
echo 'NOSMALLER40'
outfile_ccn_mod=${outfile_rootname}_CCN_NOSMALL40.nc
cdo selvar,CCN_STRAT_NOSMALL40,CCN_STRAT_NOSMALL40_CTR $infile_mod tmp.nc
cdo -setmisstoc,0 -expr,'CCN_an_strat_mo=CCN_STRAT_NOSMALL40/CCN_STRAT_NOSMALL40_CTR' tmp.nc tmp_1a.nc
# below needed for model to read in data
cdo -setunit,'m-3' tmp_1a.nc tmp_2a.nc
ncatted -a long_name,CCN_an_strat_mo,a,c,'Number of activated aerosols in stratiform clouds (acc.) - modified NOSMALL40' tmp_2a.nc tmp_3a.nc
cdo merge tmp_3a.nc tmp_1b.nc tmp_pi_strat.nc tmp_pi_conv.nc tmp_coarse.nc tmp_3d_full_mod.nc
eval cdo $operator_cdo_time_fix tmp_3d_full_mod.nc tmp_2.nc 
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_ccn_mod

# NOSMALLER1
echo 'NOSMALLER1'
outfile_ccn_mod=${outfile_rootname}_CCN_NOSMALL1.nc
cdo selvar,CCN_STRAT_NOSMALL1,CCN_STRAT_NOSMALL1_CTR $infile_mod tmp.nc
cdo -setmisstoc,0 -expr,'CCN_an_strat_mo=CCN_STRAT_NOSMALL1/CCN_STRAT_NOSMALL1_CTR' tmp.nc tmp_1a.nc
# below needed for model to read in data
cdo -setunit,'m-3' tmp_1a.nc tmp_2a.nc
ncatted -a long_name,CCN_an_strat_mo,a,c,'Number of activated aerosols in stratiform clouds (acc.) - modified NOSMALL1' tmp_2a.nc tmp_3a.nc
cdo merge tmp_3a.nc tmp_1b.nc tmp_pi_strat.nc tmp_pi_conv.nc tmp_coarse.nc tmp_3d_full_mod.nc
eval cdo $operator_cdo_time_fix tmp_3d_full_mod.nc tmp_2.nc 
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_ccn_mod

# CCN+MAX
echo 'CCNMAX'
outfile_ccn_mod=${outfile_rootname}_CCN_MAX.nc
cdo selvar,CCN_STRAT_MAX,CCN_STRAT_MAX_CTR $infile_mod tmp.nc 
cdo -setmisstoc,0 -expr,'CCN_an_strat_mo= CCN_STRAT_MAX / CCN_STRAT_MAX_CTR' tmp.nc tmp_1a.nc
# below needed for model to read in data
cdo -setunit,'m-3' tmp_1a.nc tmp_2a.nc
ncatted -a long_name,CCN_an_strat_mo,a,c,'Number of activated aerosols in stratiform clouds (acc.) - modified MAX' tmp_2a.nc tmp_3a.nc
cdo merge tmp_3a.nc tmp_1b.nc tmp_pi_strat.nc tmp_pi_conv.nc tmp_coarse.nc tmp_3d_full_mod.nc
eval cdo $operator_cdo_time_fix tmp_3d_full_mod.nc tmp_2.nc 
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_ccn_mod

# CCN+SO4KS
echo 'PLUSSO4KS'
outfile_ccn_mod=${outfile_rootname}_CCN_PLUSSO4KS.nc
cdo selvar,CCN_STRAT $infile tmp_ccn.nc
cdo selvar,CCN_STRAT_SO4KS $infile_mod tmp_so4ks.nc
cdo merge tmp_ccn.nc tmp_so4ks.nc tmp.nc
cdo -setmisstoc,0 -expr,'CCN_an_strat_mo= 0.8 * CCN_STRAT_SO4KS + CCN_STRAT' tmp.nc tmp_1a.nc
# below needed for model to read in data
cdo -setunit,'m-3' tmp_1a.nc tmp_2a.nc
ncatted -a long_name,CCN_an_strat_mo,a,c,'Number of activated aerosols in stratiform clouds (acc.) - modified plus SO4KS' tmp_2a.nc tmp_3a.nc
cdo merge tmp_3a.nc tmp_1b.nc tmp_pi_strat.nc tmp_pi_conv.nc tmp_coarse.nc tmp_3d_full_mod.nc
eval cdo $operator_cdo_time_fix tmp_3d_full_mod.nc tmp_2.nc 
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_ccn_mod

# CCN+SO4KS, only in SO
echo 'PLUSSO4KSSO'
outfile_ccn_mod=${outfile_rootname}_CCN_PLUSSO4KSSO.nc
cdo selvar,CCN_STRAT $infile tmp_ccn.nc
cdo -masklonlatbox,-180,180,-70,-30 -selvar,CCN_STRAT_SO4KS $infile_mod tmp_so4ks_cut.nc
cdo -setmisstoc,0 tmp_so4ks_cut.nc tmp_so4ks.nc
cdo merge tmp_ccn.nc tmp_so4ks.nc tmp.nc
cdo -setmisstoc,0 -expr,'CCN_an_strat_mo= 1.2 * CCN_STRAT_SO4KS + CCN_STRAT' tmp.nc tmp_1a.nc
# below needed for model to read in data
cdo -setunit,'m-3' tmp_1a.nc tmp_2a.nc
ncatted -a long_name,CCN_an_strat_mo,a,c,'Number of activated aerosols in stratiform clouds (acc.) - modified plus SO4KS in SO' tmp_2a.nc tmp_3a.nc
cdo merge tmp_3a.nc tmp_1b.nc tmp_pi_strat.nc tmp_pi_conv.nc tmp_coarse.nc tmp_3d_full_mod.nc
eval cdo $operator_cdo_time_fix tmp_3d_full_mod.nc tmp_2.nc 
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_ccn_mod

# CCN only where ACT
echo 'ACTCCN'
outfile_ccn_mod=${outfile_rootname}_CCN_ACTCCN.nc
cdo selvar,CCN_STRAT_ACTCCN,CCN_STRAT_ACTCCN_CTR $infile_mod tmp.nc 
cdo -setmisstoc,0 -expr,'CCN_an_strat_mo= CCN_STRAT_ACTCCN / CCN_STRAT_ACTCCN_CTR' tmp.nc tmp_1a.nc
# below needed for model to read in data
cdo -setunit,'m-3' tmp_1a.nc tmp_2a.nc
ncatted -a long_name,CCN_an_strat_mo,a,c,'Number of activated aerosols in stratiform clouds (acc.) - modified ACTCCN' tmp_2a.nc tmp_3a.nc
cdo merge tmp_3a.nc tmp_1b.nc tmp_pi_strat.nc tmp_pi_conv.nc tmp_coarse.nc tmp_3d_full_mod.nc
eval cdo $operator_cdo_time_fix tmp_3d_full_mod.nc tmp_2.nc 
eval ncrename $operator_nco_lev_fix tmp_2.nc $outfile_ccn_mod

fi
