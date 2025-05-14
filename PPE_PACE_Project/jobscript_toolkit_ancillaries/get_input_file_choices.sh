#!/bin/bash

#--------------------------------------------------------------------------------------------------
#
# Sylvaine Ferrachat 2012-02
#
# This script is called by jobsubm_echam.sh
#
# It parses relevant information from the echam namelists to 
# find out which input files have to be used (e.g. time-dep ozone, nudging files, etc..)
#
#--------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#-- Sea surface temperature and sea ice coverage:

#-- parse from namelist (lamip)
flag_time_dep_sst_sic=$(parse_fortran_logical ${script_dir}/namelist_${exp}.echam runctl lamip)

#-- security in case 'lamip' was not found in namelist
if [ -z $flag_time_dep_sst_sic ] ; then # flag_time_dep_sst_sic is non-assigned (ie empty??)
   flag_time_dep_sst_sic=false          # this reproduces the default in echam6 and echam5
fi

#-- set sst_sic_dataset in case it is not defined
if [ -z $sst_sic_dataset ] ; then
   sst_sic_dataset="amip"
fi

#-------------------------------------------------------------------------------------
#-- Ozone:

#-- parse from namelist (io3)
io3=$(get_from_namelist ${script_dir}/namelist_${exp}.echam radctl io3)
 
#-- security in case io3 was not found in namelist
if [ -z $io3 ] ; then
   io3=3
fi

#-- set flag_CMIP5_ozon 
if [[ "$io3" -eq "4" ]] ; then # CMIP5 ozone
   flag_CMIP5_ozon=true
else # this also applies if io3 is not defined,
   flag_CMIP5_ozon=false
fi

#-------------------------------------------------------------------------------------
#-- Nudging:

#-- parse from namelist (lnudge)
flag_nudg=$(parse_fortran_logical ${script_dir}/namelist_${exp}.echam runctl lnudge)

#-- security in case 'lnudge' was not found in namelist
if [ -z $flag_nudg ] ; then
   flag_nudg=false
fi

#-- Additional check: removes nproma from namelist if lnudge = .true.
if $flag_nudg ; then
   echo
   echo "You're running in nudged mode, removing any remaining nproma definition in "\
"${script_dir}/namelist_${exp}.echam..."

   delete_from_namelist ${script_dir}/namelist_${exp}.echam runctl nproma
fi

#-------------------------------------------------------------------------------------
#-- Hydrology:

#-- parse from namelist (with_hd)
flag_hd=$(parse_fortran_logical ${script_dir}/namelist_${exp}.jsbach jsbach_ctl with_hd)

#-- security in case 'with_hd' was not found in namelist
if [ -z $flag_hd ] ; then
   flag_hd=false
fi

#-------------------------------------------------------------------------------------
#-- Time-dep solar irradiance 

isolrad=$(get_from_namelist ${script_dir}/namelist_${exp}.echam radctl isolrad)

#-- security in case 'isolrad' was not found in namelist
if [ -z $isolrad ] ; then
   isolrad=3  # default in echam6
fi  

#-- set flag_time_dep_sol_irr
if [[ "$isolrad" -eq "1" ]] ; then # time-dep read from file
   flag_time_dep_sol_irr=true
else
   flag_time_dep_sol_irr=false
fi

#-------------------------------------------------------------------------------------
#-- aerosols for radiation

iaero=$(get_from_namelist ${script_dir}/namelist_${exp}.echam radctl iaero)

#-- security in case iaero was not found in namelist
#SF note: the script will stop if iaero is not found. In echam6, default iaero is 2
#         but I don't want to set it here to 2, since iaero=2 is obsolete

if [ -z $iaero ] ; then
   echo "Error! iaero does not seem to be defined in ${script_dir}/namelist_${exp}.echam." \
        | tee -a ${log_file}
   echo "Exiting!"
   exit 1
fi

#-- set flag_kinne_aerosols, flag_stenchikov_aerosols, flag_crowley_aerosols, flag_hamclimatology_aerosols
flag_kinne_aerosols=false
flag_stenchikov_aerosols=false
flag_crowley_aerosols=false
flag_submclim_aerosols=false

case $iaero in
   3)
     flag_kinne_aerosols=true
   ;;
   5)
     flag_kinne_aerosols=true
     flag_stenchikov_aerosols=true
   ;;
   6)
     flag_kinne_aerosols=true
     flag_stenchikov_aerosols=true
     flag_submclim_aerosols=true
   ;;
   7)
     flag_kinne_aerosols=true
     flag_crowley_aerosols=true
   ;;               
esac

#-------------------------------------------------------------------------------------
#-- Nudging format

inudgformat=$(get_from_namelist ${script_dir}/namelist_${exp}.echam ndgctl inudgformat)

#-- set flag_nudg_netcdf
if [[ "inudgformat" -eq "2" ]] ; then
   flag_nudg_netcdf=true
else
   flag_nudg_netcdf=false
fi

#-------------------------------------------------------------------------------------
#-- HAMMOZ-specific

if $this_model_flag_with_hammoz ; then

   #-- parse from namelist (lham)
   flag_ham=$(parse_fortran_logical ${script_dir}/namelist_${exp}.echam submodelctl lham)

   #-- security in case 'lham' was not found in namelist
   if [ -z $flag_ham ] ; then
      flag_ham=false
   fi

   #-- parse from namelist (nham_subm)
   nham_subm=$(get_from_namelist ${script_dir}/namelist_${exp}.echam hamctl nham_subm)

   #-- security in case 'nham_subm' was not found in namelist
   if [ -z $nham_subm ] ; then
      nham_subm=2
   fi

   #-- set aero_micro_scheme
   case "$nham_subm" in
        2) # M7
           aero_micro_scheme="M7"
        ;;
        3) # salsa
           aero_micro_scheme="salsa"
        ;;
   esac

   #-- parse from namelist (lmoz)
   flag_moz=$(parse_fortran_logical ${script_dir}/namelist_${exp}.echam submodelctl lmoz)

   #-- security in case 'lmoz' was not found in namelist
   if [ -z $flag_moz ] ; then
      flag_moz=false
   fi

   #-- Add support for HAMMOZ
   flag_hammoz=$(parse_fortran_logical ${script_dir}/namelist_${exp}.echam submodelctl lhammoz)

   #-- security in case 'lhammoz' was not found in namelist
   if [ -z $flag_hammoz ] ; then
      flag_hammoz=false
   fi

   if $flag_hammoz ; then
      flag_ham=true
      flag_moz=true
   fi

   #-- set moz-dependent flags
   if $flag_moz ; then

      #-- parse from namelist (lphotolysis)
      flag_photolysis=$(parse_fortran_logical ${script_dir}/namelist_${exp}.echam mozctl lphotolysis)

      #-- security in case 'lphotolysis' was not found in namelist
      if [ -z $flag_photolysis ] ; then
         flag_photolysis=true
      fi

      #-- parse from namelist (lfastj)
      flag_fastj=$(parse_fortran_logical ${script_dir}/namelist_${exp}.echam mozctl lfastj)

      #-- security in case 'lfastj' was not found in namelist
      if [ -z $flag_fastj ] ; then
         flag_fastj=false
      fi

   fi

fi # end $flag_hammoz

#-------------------------------------------------------------------------------------
#-- end (normal)

echo | tee -a $log_file
