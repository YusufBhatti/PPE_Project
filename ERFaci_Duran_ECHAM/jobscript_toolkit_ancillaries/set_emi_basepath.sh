#!/bin/bash

#--------------------------------------------------------------------------------------------------
#
# Sylvaine Ferrachat 2016-03
#
# This script is called by jobsubm_echam.sh
#
# It sets the correct path to the emission files which are specified in the emi_spec file.
#
#--------------------------------------------------------------------------------------------------

namelist_echam=${script_dir}/namelist_${exp}.echam

#-- Check if emi_basepath has already a valid value (existing directory) and fixes it
#   in case it does not 
emi_basepath="$(get_from_namelist $namelist_echam submodelctl emi_basepath)"

# Remove prohibited characters for checking if it's a directory
emi_basepath=$(sed -e 's|\[||g;s|\]||g;s| ||g;s|\"||g' <<< "$emi_basepath")

if [ ! -d "$emi_basepath" ] ; then # emi_basepath is not a directory

   emi_basepath=${input_basepath}/emissions_inventories/

   echo "----" | tee -a $log_file
   echo "Setting emi_basepath to:" | tee -a $log_file
   echo "$emi_basepath" | tee -a $log_file
   echo "in" | tee -a $log_file
   echo "./$(rel_path ${namelist_echam} $current_wd_at_beg)" | tee -a $log_file
   echo | tee -a $log_file

fi          

# Note that the command below will be harmless in case emi_basepath already existed as a valid 
# directory (and this allows to add security with quotes)
replace_in_namelist $namelist_echam submodelctl emi_basepath "\"$emi_basepath\"" 

