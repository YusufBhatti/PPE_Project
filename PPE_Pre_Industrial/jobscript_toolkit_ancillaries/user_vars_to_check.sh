#!/bin/bash

#--------------------------------------------------------------------------------------------------
#
# Sylvaine Ferrachat 2016-03
#
# This script defines the list of user-defined (and model-specific) variables that will be checked
#
#--------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
#-- Perform basic tests (mostly existence checks) on model-specific user-defined vars

list_var_basic_tests=( \
    "hres   | -z | Please set 'hres' in your settings file." \
    "vres   | -z | Please set 'vres' in your settings file." \
    "oceres | -z | Please set 'oceres' in your settings file." \
    "ntiles | -z | Please set 'ntiles' in your settings file." \
    "date_start | -z | Please set 'date_start' in your settings file." \
    "date_stop  | -z | Please set 'date_stop' in your settings file." \
                     )

if $this_model_flag_emi ; then
   list_var_basic_tests+=( \
       "emissions_spec | -z   | Please set 'emissions_spec' in your settings file." \
       "emissions_spec | ! -f | Please reset 'emissions_spec' in your settings file." \
                         )
fi

basic_var_tests "${list_var_basic_tests[@]}"

#----------------------------------------------------------------------
#-- Perform pattern-matching tests on model-specific user-defined vars

nproma_pattern="^[0-9]+$"
nproma_message="Please reset 'nproma' in your settings file (should be a number)."

hres_pattern="^T[0-9]+$"
hres_message="Please reset 'hres' in your settings file (format: TX, with X being a number)"

vres_pattern="^L[0-9]+$"
vres_message="Please reset 'vres' in your settings file (format: TX, with X being a number)"

ntiles_pattern="^[0-9]+$"
ntiles_message="Please reset 'ntiles' in your settings file (should be a number)."

date_start_pattern="^[0-9]{1,4},[0-9]{2},[0-9]{2},[0-9]{1,2},[0-9]{1,2},[0-9]{1,2}$"
date_start_message="Please reset 'date_start' in your settings file (format: YYYY,MM,DD,HH,mm,SS)"

date_stop_pattern="^[0-9]{1,4},[0-9]{2},[0-9]{2},[0-9]{1,2},[0-9]{1,2},[0-9]{1,2}$"
date_stop_message="Please reset 'date_stop' in your settings file (format: YYYY,MM,DD,HH,mm,SS)"

list_var_patt_tests=( \
    "nproma | $nproma_pattern | $nproma_message" \
    "hres | $hres_pattern | $hres_message" \
    "vres | $vres_pattern | $vres_message" \
    "ntiles | $ntiles_pattern | $ntiles_message" \
    "date_start | $date_start_pattern | $date_start_message" \
    "date_stop | $date_stop_pattern | $date_stop_message" \
                    )

pattern_var_tests "${list_var_patt_tests[@]}"

#----------------------------------------------------------
#-- check that exp_dir finishes with '/' 
#   (and correct it in the namelist if necessary)
((index_last_char=${#dir}-1))
last_char=${exp_dir:$index_last_char:1}
if [[ "$last_char" != "/" ]] ; then # last character of $exp_dir does not finish by '/'
                                    # need to correct that in the echam namelist

   replace_in_namelist ${script_dir}/namelist_${exp}.echam \
                       "runctl" \
                       "out_datapath" \
                       "\"${exp_dir}/\"" 
fi

