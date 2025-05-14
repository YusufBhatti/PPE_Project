#!/bin/bash

#----------------------------------------------------------------------------
# Sylvaine Ferrachat / Grazia Frontoso 2011
#
# Script to perform the linking of input files into the experiment directory
# for echam6-hammoz (and echam6) 
#
# This script is called jobsum_echam.sh
#
# Variables that are used in this script (and defined prior to executing it):
#
#    exp_dir="[path to your experiment directory here]"
#    log_file="[file to log the current process]"
#    hres="[horizontal resolution, e.g. T63]"
#    vres="[vertical resolution, e.g. L47]"
#    oceres="[ocean resolution, e.g. GR30]"
#    start_yyyymm="[starting date as YYYYMM]"
#    stop_yyyymm="[stopping date as YYYYMM]"
#    start_year_m1="[year before starting year]"
#    stop_year_p1="[year after stopping year]"
#    start_yyyymm_m1="[1 month before starting date as YYYYMM]"
#    stop_yyyymm_p1="[1 month after stopping date  as YYYYMM]"
#    flag_time_dep_sst_sic="[time-dependence of SST and SIC: 'false' (climatologic) or 'true']"
#    sst_sic_dataset="[SST/SIC dataset: amip (default) or XXXX _no other choice for the moment_]"
#    flag_CMIP5_ozon="[usage of CMIP5 ozon: 'false' (obsolete climatologic input) or 'true' (CMIP5)]"
#    flag_time_dep_sol_irr="[time-dependence of solar irradiance: 'false' (climatologic) or 'true' (time-dep)]"
#    flag_kinne_aerosols="[flag to use Kinne aerosols (radiative prop + CCN (time-dep)): 'false' (no Kinne aerosols) or 'true']"
#    flag_stenchikov_aerosols="[flag to use Stenchikov volcanic aerosols (radiative prop): 'false' (no volcano aerosols) or 'true']"
#    flag_crowley_aerosols="[flag to use crowley volcanic aerosols (radiative prop): 'false' (no volcano aerosols) or 'true']"
#    flag_submclim_aerosols="[flag to use submodel climatology volcanic aerosols (radiative prop): 'false' (no volcano aerosols) or 'true']"
#    scenario="[future climate scenarios for GHG, ozone and aerosol climatology: RCP45, XXX (case insensitive)]"
#    flag_nudg="[flag to switch on nudging: 'false' (no nudging) or 'true']"
#    flag_nudg_netcdf="[flag to switch on netcdf format for nudging: 'false' (binary) or 'true' (netcdf)]"
#    flag_hd="[flag to switch on hydrology in jsbach: 'false' (no hydrology) or 'true']"
#    aero_micro_scheme="[aerosol microphysics scheme: 'M7' or 'salsa']" --> only if echam6-hammoz
#    flag_ham="[flag to use HAM-specific files: 'false' or 'true']"
#    flag_moz="[flag to use MOZ-specific files: 'false' or 'true']"
#    flag_photolysis="[flag to use photolysis: 'false' or 'true']" --> only relevant for MOZ
#    flag_fastj="[flag to use fastj: 'false' or 'true']" --> only relevant for MOZ
#    trac_init_date="[date for tracer intialization: 'jan2003' or other date with similar pattern]" --> only relevant for MOZ
#
#----------------------------------------------------------------------------

#----------------------------------------------
# Machine-specific (customize your path here):
#----------------------------------------------

input_basepath=XX/${input_files_version}      # where all input files except nudging data is to be found
nudg_basepath=XX                              # where all nudging data is to be found

#--------------------------------------------------
# You shouldn't need to modify the following
# unless you use a non-standard file organization
# or additional options for specific datasets 
#--------------------------------------------------

echo "----" | tee -a $log_file
echo "Starting input files linking process into $exp_dir..." | tee -a $log_file

#------------------------------
# Initialization / Flag checks
#------------------------------

  # Note: the flag_check function is stored in [toolkit root path]/lib/parse_text_utils.sh
  #       it takes 3 or more arguments:
  #        * a logical flag for the check to be case insensitive (true) or not (false)
  #        * the string name of the variable to check (!!not the variable itself!!)
  #        * a list of all possible values

  #-- flag_time_dep_sst_sic

      # flag_time_dep_sst_sic is automatically set prior to executing this script, based on lamip:
      #   lamip=.false. --> flag_time_dep_sst_sic=false        
      #   lamip=.true.  --> flag_time_dep_sst_sic=true
      # uncomment the following line to override the automatic setting:
      # flag_time_dep_sst_sic=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_time_dep_sst_sic ] ; then
         flag_check $case_insensitive flag_time_dep_sst_sic "${possible_values[@]}"
      else
         echo "Error! flag_time_dep_sst_sic is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- sst_sic_dataset

      # sst_sic_dataset is set by default to 'amip', you may set an alternate value below:
      # uncomment the following line to override the automatic setting:
      sst_sic_dataset="amip" #SF tmp

      possible_values=( amip ) # add additional values in case you want to use an alternate sst/sic dataset
      case_insensitive=true
      if [ ! -z $sst_sic_dataset ] ; then
         flag_check $case_insensitive sst_sic_dataset "${possible_values[@]}"
      else
         echo "Error! sst_sic_dataset is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- flag_CMIP5_ozon

      # flag_CMIP5_ozon is automatically set prior to executing this script, based on io3:
      #  if io3=4   --> flag_CMIP5_ozon=true
      #  flag_CMIP5_ozon=false otherwise
      # uncomment the following line to override the automatic setting:
      # flag_CMIP5_ozon=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_CMIP5_ozon ] ; then
         flag_check $case_insensitive flag_CMIP5_ozon "${possible_values[@]}"
      else
         echo "Error! flag_CMIP5_ozon is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- flag_time_dep_sol_irr

      # flag_time_dep_sol_irr is automatically set prior to executing this script, based on isolrad:
      # if isolrad=1   --> flag_time_dep_sol_irr=true
      # flag_time_dep_sol_irr=false otherwise
      # uncomment the following line to override the automatic setting:
      # flag_time_dep_sol_irr=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_time_dep_sol_irr ] ; then
         flag_check $case_insensitive flag_time_dep_sol_irr "${possible_values[@]}"
      else
         echo "Error! flag_time_dep_sol_irr is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- flag_kinne_aerosols

      # flag_kinne_aerosols is automatically set prior to executing this script, based on iaero:
      # if iaero=3,5,6 or 7 --> flag_kinne_aerosols=true
      # flag_kinne_aerosols=false otherwise  
      # uncomment the following line to override the automatic setting:
      # flag_kinne_aerosols=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_kinne_aerosols ] ; then
         flag_check $case_insensitive flag_kinne_aerosols "${possible_values[@]}"
      else
         echo "Error! flag_kinne_aerosols is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- flag_stenchikov_aerosols

      # flag_stenchikov_aerosols is automatically set prior to executing this script, based on iaero:
      # if iaero=5 or 6 --> flag_stenchikov_aerosols=true
      # uncomment the following line to override the automatic setting:
      # flag_stenchikov_aerosols=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_stenchikov_aerosols ] ; then
         flag_check $case_insensitive flag_stenchikov_aerosols "${possible_values[@]}"
      else
         echo "Error! flag_stenchikov_aerosols is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi
   
  #-- flag_crowley_aerosols

      # flag_crowley_aerosols is automatically set prior to executing this script, based on iaero:
      # if iaero=7 --> flag_crowley_aerosols=true
      # uncomment the following line to override the automatic setting:
      # flag_crowley_aerosols=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_crowley_aerosols ] ; then
         flag_check $case_insensitive flag_crowley_aerosols "${possible_values[@]}"
      else
         echo "Error! flag_crowley_aerosols is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- flag_submclim_aerosols

      # flag_submclim_aerosols is automatically set prior to executing this script, based on iaero:
      # if iaero=6 --> flag_submclim_aerosols=true
      # uncomment the following line to override the automatic setting:
      # flag_submclim_aerosols=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_submclim_aerosols ] ; then
         flag_check $case_insensitive flag_submclim_aerosols "${possible_values[@]}"
      else
         echo "Error! flag_submclim_aerosols is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- scenario

      # scenario is defined in the settings file:

      possible_values=( "rcp[0-9][0-9]" "historic" ) # add additional patterns in case you want to use an 
                                          # alternate scenario
      case_insensitive=true
      if [ ! -z $scenario ] ; then
         flag_check $case_insensitive scenario "${possible_values[@]}"
      else
         echo "Error! scenario is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- flag_nudg

      # flag_nudg is automatically set prior to executing this script, based on lnudge:
      # if lnudge=.true.  --> flag_nudg=true
      # if lnudge=.false. --> flag_nudg=false
      # uncomment the following line to override the automatic setting:
      # flag_nudg=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_nudg ] ; then
         flag_check $case_insensitive flag_nudg "${possible_values[@]}"
      else
         echo "Error! flag_nudg is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- flag_nudg_netcdf

      # flag_nudg_netcdf is automatically set prior to executing this script, based on inudgformat:
      # if inudgformat=2 --> flag_nudg_netcdf=true
      # flag_nudg_netcdf=false otherwise
      # uncomment the following line to override the automatic setting:
      # flag_nudg_netcdf=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_nudg_netcdf ] ; then
         flag_check $case_insensitive flag_nudg_netcdf "${possible_values[@]}"
      else
         echo "Error! flag_nudg_netcdf is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

 #-- flag_hd

      # flag_hd is automatically set prior to executing this script, based on lnudge:
      # if with_hd=.true.  --> flag_hd=true
      # if with_hd=.false. --> flag_hd=false
      # uncomment the following line to override the automatic setting:
      # flag_hd=XXX

      possible_values=( true false )
      case_insensitive=true
      if [ ! -z $flag_hd ] ; then
         flag_check $case_insensitive flag_hd "${possible_values[@]}"
      else
         echo "Error! flag_hd is not defined! Please review your symlinks script!" | tee -a $log_file
         echo "Exiting"
         exit 1
      fi

  #-- echam6-hammoz specific:

      if $this_model_flag_with_hammoz ; then

         #-- aero_micro_scheme 

         # aero_micro_scheme is automatically set set prior to executing this script, based on nham_subm in the
         # settings file. Default value is 'M7'.
         possible_values=( M7 salsa ) 
         case_insensitive=false
         if [ ! -z $aero_micro_scheme ] ; then
            flag_check $case_insensitive aero_micro_scheme "${possible_values[@]}"
         else
            echo "Error! aero_micro_scheme is not defined! Please review your symlinks script!" | tee -a $log_file
            echo "Exiting"
            exit 1
         fi

         #-- flag_ham
   
         # flag_ham is automatically set prior to executing this script, based on lham:
         # if lham=.true.  --> flag_ham=true
         # if lham=.false. --> flag_ham=false
         # uncomment the following line to override the automatic setting:
         # flag_ham=XXX
   
            possible_values=( true false )
            case_insensitive=true
            if [ ! -z $flag_ham ] ; then
               flag_check $case_insensitive flag_ham "${possible_values[@]}"
            else
               echo "Error! flag_ham is not defined! Please review your symlinks script!" | tee -a $log_file
               echo "Exiting"
               exit 1
            fi
   
         #-- flag_moz
   
         # flag_moz is automatically set prior to executing this script, based on lmoz:
         # if lmoz=.true.  --> flag_moz=true
         # if lmoz=.false. --> flag_moz=false
         # uncomment the following line to override the automatic setting:
         # flag_moz=XXX
   
            possible_values=( true false )
            case_insensitive=true
            if [ ! -z $flag_moz ] ; then
               flag_check $case_insensitive flag_moz "${possible_values[@]}"
            else
               echo "Error! flag_moz is not defined! Please review your symlinks script!" | tee -a $log_file
               echo "Exiting"
               exit 1
            fi
   
            if $flag_moz ; then

               #-- flag_photolysis
         
               # flag_photolysis is automatically set prior to executing this script, based on lphotolysis:
               # if lphotolysis=.true.  --> flag_photolysis=true
               # if lphotolysis=.false. --> flag_photolysis=false
               # uncomment the following line to override the automatic setting:
               # flag_photolysis=XXX
         
               possible_values=( true false )
               case_insensitive=true
               if [ ! -z $flag_photolysis ] ; then
                  flag_check $case_insensitive flag_photolysis "${possible_values[@]}"
               else
                  echo "Error! flag_photolysis is not defined! Please review your symlinks script!" | tee -a $log_file
                  echo "Exiting"
                  exit 1
               fi

               #-- flag_fastj
         
               # flag_fastj is automatically set prior to executing this script, based on lfastj:
               # if lfastj=.true.  --> flag_fastj=true
               # if lfastj=.false. --> flag_fastj=false
               # uncomment the following line to override the automatic setting:
               # flag_fastj=XXX
         
               possible_values=( true false )
               case_insensitive=true
               if [ ! -z $flag_fastj ] ; then
                  flag_check $case_insensitive flag_fastj "${possible_values[@]}"
               else
                  echo "Error! flag_fastj is not defined! Please review your symlinks script!" | tee -a $log_file
                  echo "Exiting"
                  exit 1
               fi

#>>csld irrelevant for pi Control run
#               #-- trac_init_date 
#
#               # trac_init_date should be defined in the settings file.
#
#               possible_values=( "[a-z][a-z][a-z][0-9][0-9][0-9][0-9]" ) 
#               case_insensitive=true
#               if [ ! -z $trac_init_date ] ; then
#                  flag_check $case_insensitive trac_init_date "${possible_values[@]}"
#               else
#                  echo "Error! trac_init_date is not defined! Please review your symlinks script!" | tee -a $log_file
#                  echo "Exiting"
#                  exit 1
#               fi
#<<csld
            fi

      fi

#---------------------------------
# Summary of necessary variables:
#---------------------------------

echo                                                  | tee -a $log_file
echo "with:"                                          | tee -a $log_file
echo                                                  | tee -a $log_file
echo "input_basepath: $input_basepath"                | tee -a $log_file
echo "nudg_basepath: $nudg_basepath"                  | tee -a $log_file
echo "exp_dir: $exp_dir"                              | tee -a $log_file
echo "hres: $hres"                                    | tee -a $log_file
echo "vres: $vres"                                    | tee -a $log_file
echo "oceres: $oceres"                                | tee -a $log_file
echo "start_yyyymm: $start_yyyymm"                    | tee -a $log_file
echo "stop_yyyymm: $stop_yyyymm"                      | tee -a $log_file
echo "start_year_m1: $start_year_m1"                  | tee -a $log_file
echo "stop_year_p1: $stop_year_p1"                    | tee -a $log_file
echo "start_yyyymm_m1: $start_yyyymm_m1"              | tee -a $log_file
echo "stop_yyyymm_p1: $stop_yyyymm_p1"                | tee -a $log_file
echo "flag_time_dep_sst_sic: $flag_time_dep_sst_sic"  | tee -a $log_file
echo "sst_sic_dataset: $sst_sic_dataset"              | tee -a $log_file
echo "flag_CMIP5_ozon: $flag_CMIP5_ozon"              | tee -a $log_file
echo "flag_time_dep_sol_irr: $flag_time_dep_sol_irr"  | tee -a $log_file
echo "flag_kinne_aerosols: $flag_kinne_aerosols"      | tee -a $log_file
echo "flag_stenchikov_aerosols: $flag_stenchikov_aerosols" | tee -a $log_file
echo "flag_crowley_aerosols: $flag_crowley_aerosols"       | tee -a $log_file
echo "flag_submclim_aerosols: $flag_submclim_aerosols"     | tee -a $log_file
echo "scenario: $scenario"                            | tee -a $log_file
echo "flag_nudg: $flag_nudg"                          | tee -a $log_file
echo "flag_nudg_netcdf: $flag_nudg_netcdf"            | tee -a $log_file
echo "flag_hd: $flag_hd"                          | tee -a $log_file

if $this_model_flag_with_hammoz ; then # echam6-hammoz specific
   echo "aero_micro_scheme: $aero_micro_scheme"  | tee -a $log_file
   echo "flag_ham: $flag_ham"  | tee -a $log_file
   echo "flag_moz: $flag_moz"  | tee -a $log_file
   if $flag_moz ; then
      echo "flag_photolysis: $flag_photolysis"  | tee -a $log_file
      echo "flag_fastj: $flag_fastj"  | tee -a $log_file
#csld       echo "trac_init_date: $trac_init_date" | tee -a $log_file
   fi
fi
echo | tee -a $log_file

#---------
# Process
#---------

cd $exp_dir

#--------------------------------------
# Removing pre-existing symbolic links 
#--------------------------------------

#SF this should not be necessary, but is conserved for security (in case 'ln' is not GNU ln)
#SF exclude the links to restart files

find . -type l -and -not -name "${prefix_rerun_file}_${exp}_*${suffix_rerun_file}" -exec \rm -f {} \;

#--------
# ECHAM6  
#--------

ln -sf ${input_basepath}/echam6/${hres}/${hres}${vres}_jan_spec.nc           unit.23
ln -sf ${input_basepath}/echam6/${hres}/${hres}${oceres}_jan_surf.nc         unit.24
ln -sf ${input_basepath}/echam6/${hres}/${hres}${oceres}_VLTCLIM.nc          unit.90
ln -sf ${input_basepath}/echam6/${hres}/${hres}${oceres}_VGRATCLIM.nc        unit.91
ln -sf ${input_basepath}/echam6/${hres}/${hres}_TSLCLIM.nc                   unit.92

ln -sf  ${input_basepath}/echam6/rrtmg_sw.nc               rrtmg_sw.nc
ln -sf  ${input_basepath}/echam6/rrtmg_lw.nc               rrtmg_lw.nc
ln -sf  ${input_basepath}/echam6/ECHAM6_CldOptProps.nc     ECHAM6_CldOptProps.nc

ln -sf ${input_basepath}/echam6/greenhouse_${scenario}.nc  greenhouse_gases.nc

#---------
# JS-BACH  
#---------

ln -sf  ${input_basepath}/echam6/jsbach/lctlib_nlct21.def_rev7624 lctlib.def
#CSW
#ln -sf  ${input_basepath}/jsbach/${hres}/jsbach_${hres}${oceres}_${ntiles}tiles_5layers_2005.nc jsbach.nc
#ln -sf  ${input_basepath}/jsbach/${hres}/jsbach_${hres}${oceres}_${ntiles}tiles_5layers_1850.nc jsbach.nc
ln -sf  ${input_basepath}/jsbach/${hres}/jsbach_${hres}${oceres}_${ntiles}tiles_5layers_1976.nc jsbach.nc

if $flag_hd ; then

    ln -sf  ${input_basepath}/jsbach/HD/hdpara.nc          hdpara.nc
    ln -sf  ${input_basepath}/jsbach/HD/hdstart.nc         hdstart.nc

fi

#--------------------------
# Climatologic SST and SIC
#--------------------------
#CSW
ln -sf ${input_basepath}/echam6/${hres}/amip/${hres}_amipsst_1979-2008_mean.nc unit.20
ln -sf ${input_basepath}/echam6/${hres}/amip/${hres}_amipsic_1979-2008_mean.nc unit.96
#ln -sf ${input_basepath}/echam6/${hres}/${hres}${oceres}_piControl-LR_sst_1880-2379.nc unit.20
#ln -sf ${input_basepath}/echam6/${hres}/${hres}${oceres}_piControl-LR_sic_1880-2379.nc unit.96


#----------------------------
# Time-dep SST & SIC 
#----------------------------

if $flag_time_dep_sst_sic ; then

   if [[ "$sst_sic_dataset" == "amip" ]] ; then # use amip2 data

       for ((year=$start_year_m1; year<=$stop_year_p1; year++)) ; do
           ln -sf ${input_basepath}/echam6/${hres}/amip/${hres}_amipsst_${year}.nc sst${year}
           ln -sf ${input_basepath}/echam6/${hres}/amip/${hres}_amipsic_${year}.nc ice${year}
       done

   # uncomment the following in case of an alternate SST / SIC dataset:

   # elif [[ "$sst_sic_dataset" == "XXXX" ]] ; then # use XXXX data

       #for ((year=$start_year_m1; year<=$stop_year_p1; year++)) ; do
       #    ln -sf ${input_basepath}/echam6/${hres}/XXXX/XXXX.nc sst${year}
       #    ln -sf ${input_basepath}/echam6/${hres}/XXXX/XXXX.nc ice${year}
       #done
   fi

fi

#----------------
# Time-dep ozone
#----------------
#
# !!WARNING!! (SF, 2012.08) (SF update 2016.01)
# Due to a limitation in the current echam6 code, the only way to use climatologic CMIP5 ozone
# is to fake a time-dependent input, because the former reading of 'unit.21' which was normally meant for
# climatologic ozone, is coded in a way that is only compliant to a non-supported nput file 
# (former ${hres}_O3clim2.nc, which has been removed from the echam input distrib)! 
#
# In order to decipher whether the user wants a *climatology* or year-dependent input, I chose to bound it
# to temp-dependent SSTs: 
#     - if $flag_time_dep_sst_sic is true (ie lamip=.true.), then ozone files will be truely year-dependent.
#     - if $flag_time_dep_sst_sic is false (ie lamip=.false.), then ozone files will all point to the ozone 
#        CMIP5 climatology
#
# This should make sense, but please adapt it if this is bad strategy for purposes!!
#

if $flag_CMIP5_ozon ; then  # ie io3=4

   for ((year=$start_year_m1; year<=$stop_year_p1; year++)) ; do

       if $flag_time_dep_sst_sic ; then  # true year-dep input, see comment above

          if [[ "$year" -lt "2009" ]] ; then # no scenario
             scenario_str=""
          else                               # take user-defined scenario
             scenario_str=`echo $scenario | tr '[:lower:]' '[:upper:]'`"_"
          fi
          ln -sf ${input_basepath}/echam6/${hres}/ozone/${hres}_ozone_CMIP5_${scenario_str}${year}.nc ozon${year}

       else # fake year-dep input for using the CMIP5 climatology
          ln -sf ${input_basepath}/echam6/${hres}/${hres}_ozone_CMIP5_1979-1988.nc ozon${year}
       fi
   done

fi

#---------------------------
# Time-dep solar irradiance
#---------------------------

if $flag_time_dep_sol_irr ; then

   for ((year=$start_year_m1; year<=$stop_year_p1; year++)) ; do
       ln -sf ${input_basepath}/echam6/solar_irradiance/swflux_14band_${year}.nc swflux_${year}.nc
       ln -sf ${input_basepath}/hammoz/photolysis/etfphot_${year}.nc etfphot_${year}.nc
   done

fi

#-----------------------------------------------------
# Time-dep Kinne aerosols (including CCN climatology) 
#-----------------------------------------------------

if $flag_kinne_aerosols ; then

   for ((year=$start_year_m1; year<=$stop_year_p1; year++)) ; do
       if [[ "$year" -lt "2001" ]] ; then # no scenario
          scenario_str=""
       else                               # take user-defined scenario
          scenario_str=`echo $scenario | tr '[:upper:]' '[:lower:]' `"_"
       fi

       ln -sf ${input_basepath}/echam6/${hres}/aero/${hres}_aeropt_kinne_sw_b14_fin_${scenario_str}${year}.nc aero_fine_${year}.nc
       ln -sf ${input_basepath}/echam6/${hres}/aero/${hres}_aeropt_kinne_sw_b14_coa.nc aero_coarse_${year}.nc
       ln -sf ${input_basepath}/echam6/${hres}/aero/${hres}_aeropt_kinne_lw_b16_coa.nc aero_farir_${year}.nc

       #>> SFtemp: prepare link for CCN data when available
       #ln -sf ${input_basepath}/echam6/${hres}/aero/CCN_${hres}_1km-8km_${year}.nc ccn${year}.nc
       #<< SFtemp
   done

fi

#---------------------------------------
# Time-dep Stenchikov volcanic aerosols  
#---------------------------------------

if $flag_stenchikov_aerosols ; then

   for ((year=$start_year_m1; year<=$stop_year_p1; year++)) ; do

       ln -sf ${input_basepath}/echam6/${hres}/volcano_aerosols/strat_aerosol_sw_${hres}_${year}.nc strat_aerosol_sw_${year}.nc
       ln -sf ${input_basepath}/echam6/${hres}/volcano_aerosols/strat_aerosol_ir_${hres}_${year}.nc strat_aerosol_ir_${year}.nc

   done

fi

#---------------------------
# Crowley volcanic aerosols  
#---------------------------

if $flag_crowley_aerosols ; then

   ln -sf ${input_basepath}/echam6/volc_data aodreff_crow.dat  # SF note: volc_data is referred to as
                                                               # "ici5d-ad800-1999.asc" in the echam6
                             							       # users guide	
   ln -sf ${input_basepath}/echam6/b30w120 aero_volc_tables.dat

fi

#------------------------------------------------------------
# Submodel climatology for volcanic aerosols (e.g. from HAM)  
#------------------------------------------------------------

#SF: the following does not exist yet!
#if $flag_submclim_aerosols ; then

   #ln -sf ${input_basepath}/echam6/XXXX/XXXX.nc aoddz_ham_yyyy.nc

#fi

#-------
# ISCCP  
#-------
   
#>>SF: this is not re-implemented in echam6 yet
#ln -sf ${input_basepath}/echam6/invtau.formatted invtau.formatted
#ln -sf ${input_basepath}/echam6/tautab.formatted tautab.formatted
#<<SF

#-------
# CFMIP
#-------

ln -sf ${input_basepath}/echam6/CFMIP/pointlocations.txt	pointlocations.txt
   
#-----------------
# HAMMOZ-specific
#-----------------

if $this_model_flag_with_hammoz ; then # echam6-hammoz specific
                                       # Note that this does not necessarily implies that HAM and/or MOZ
                                       # are effectively used.

   #---------------------------------
   # soil and surface properties
   #---------------------------------

   ln -sf ${input_basepath}/hammoz/${hres}/soilpHfrac_${hres}.nc xtsoil.nc
   ln -sf ${input_basepath}/hammoz/${hres}/xtsurf_v2_${hres}.nc surface_properties.nc

   #-------------------------
   # Emissions specification 
   #-------------------------

   # --> All specs taken care of by emi_spec.txt.
   # --> The setting of the basepath, emi_basepath, is handled by jobsubm_echam.sh
   # --> Note that a few emission-dependent files are nevertheless necessary here, see below
   
   #--------------------
   # Volcanic emissions
   #--------------------
  
   ln -sf ${input_basepath}/emissions_inventories/aerocom_II/explosive_volcanos.dat explosive_volcanos.dat
   ln -sf ${input_basepath}/emissions_inventories/aerocom_II/continuous_volcanos.dat continuous_volcanos.dat
   
   #--------------------------
   # Oceanic emissions of DMS
   #--------------------------
   
   ln -sf ${input_basepath}/emissions_inventories/others/${hres}/emiss_fields_dms_sea_monthly_${hres}.nc conc_aerocom_DMS_sea.nc 
   
   #-----------------------------------------------------
   # MEGANV2 emission factors for biogenic VOC emissions
   #----------------------------------------------------
  
   #SF:
   #ToDo: this part might later be conditionned to whether it is really needed (ie having online biogenic emissions) 

   ln -sf ${input_basepath}/hammoz/${hres}/megan_ef_specific_${hres}.nc            megan_emission_factors.nc   
   ln -sf ${input_basepath}/hammoz/${hres}/megan_mksrf_pft_cesm4_2005_${hres}.nc   megan_clm_pft_map.nc
   ln -sf ${input_basepath}/hammoz/${hres}/megan_fraction_broadl_de_${hres}.nc     megan_fraction_broadl_de.nc
   ln -sf ${input_basepath}/hammoz/${hres}/megan_fraction_broadl_ev_${hres}.nc     megan_fraction_broadl_ev.nc
   ln -sf ${input_basepath}/hammoz/${hres}/megan_boreal_trees_limit_${hres}.nc     megan_boreal_trees_limit.nc
   ln -sf ${input_basepath}/hammoz/${hres}/megan_boreal_grasses_limit_${hres}.nc   megan_boreal_grasses_limit.nc

   #--------------
   # HAM-specific
   #--------------

   if $flag_ham ; then

      #---------------------------------------------------
      # Boundary conditions for the dust emissions scheme
      #---------------------------------------------------

      ln -sf ${input_basepath}/hammoz/${hres}/dust_preferential_sources_${hres}.nc  dust_preferential_sources.nc
      ln -sf ${input_basepath}/hammoz/${hres}/dust_potential_sources_${hres}.nc     dust_potential_sources.nc
      ln -sf ${input_basepath}/hammoz/${hres}/soil_type_all_${hres}.nc              soil_type_all.nc
      ln -sf ${input_basepath}/hammoz/${hres}/surface_rough_12m_${hres}.nc          surface_rough_12m.nc
      ln -sf ${input_basepath}/hammoz/${hres}/msg_pot_sources_${hres}.nc            dust_msg_pot_sources.nc
      ln -sf ${input_basepath}/hammoz/${hres}/dust_regions_${hres}.nc               dust_regions.nc

      #----------------------------------------------------------------
      # Prescribed oxidant concentrations for the HAM chemistry scheme 
      #----------------------------------------------------------------
   
      if ! $flag_moz ; then
         ln -sf ${input_basepath}/hammoz/${hres}/ham_oxidants_monthly_${hres}${vres}_macc.nc ham_oxidants_monthly.nc
      fi

      #--------------------------------------
      # Ion production rate from cosmic rays 
      #--------------------------------------
   
      ln -sf ${input_basepath}/hammoz/solmin.txt gcr_ipr_solmin.txt
      ln -sf ${input_basepath}/hammoz/solmax.txt gcr_ipr_solmax.txt
      
      #--------------------------------------------
      # PARNUC neutral and charge nucleation scheme
      #--------------------------------------------
      
      ln -sf ${input_basepath}/hammoz/parnuc.15H2SO4.A0.total.nc parnuc.15H2SO4.nc
      
      #--------------------------------------------------
      # Aerosol water uptake (using Kappa-Koehler theory)
      #--------------------------------------------------
      
      ln -sf ${input_basepath}/hammoz/lut_kappa.nc lut_kappa.nc
      
      #---------------------------
      # Aerosol optical properties
      #---------------------------
      
      ln -sf ${input_basepath}/hammoz/lut_optical_properties_${aero_micro_scheme}.nc    lut_optical_properties.nc
      ln -sf ${input_basepath}/hammoz/lut_optical_properties_lw_${aero_micro_scheme}.nc lut_optical_properties_lw.nc

   fi

   #--------------
   # MOZ-specific
   #--------------

   if $flag_moz ; then

      #-----------------------
      # Tracer initialization
      #-----------------------
#CSW
      #ln -sf  ${input_basepath}/hammoz/${hres}/moz_initial_jan2003_${hres}${vres}.nc tracer_ic.nc
      ln -sf  ${input_basepath}/hammoz/${hres}/moz_B1955WCN_transient_1980-01_tracer_ic_${hres}${vres}.nc tracer_ic.nc

      #---------------------------------
      # UV albedo data (green and white)
      #---------------------------------

      ln -sf ${input_basepath}/hammoz/${hres}/moz_uvalbedo_${hres}.nc moz_uvalbedo.${hres}.nc

      #--------------------
      # aerosol climatology
      #--------------------

      if ! $flag_ham ; then
         ln -sf ${input_basepath}/hammoz/${hres}/het_chem_2000-2010_${hres}${vres}.nc ham_aerosol_climatology.${hres}${vres}.nc
      fi
    
      #------------------------------------------
      # stratospheric aerosol density climatology
      #------------------------------------------

      ln -sf ${input_basepath}/hammoz/${hres}/moz_sad_sulf_1850-2100_c080220_${hres}${vres}.nc moz_sad_sulf.${hres}.nc

      #---------------------------------
      # generic lower boundary condition
      #---------------------------------

      scenario_str=`echo $scenario | tr '[:lower:]' '[:upper:]'`

      ln -sf ${input_basepath}/hammoz/${hres}/hammoz_lbc_CCMI_${scenario_str}_za_c130313_${hres}.nc moz_lbc.${hres}.nc

      #------------
      # Photolysis
      #------------

      if $flag_photolysis ; then
         ln -sf ${input_basepath}/hammoz/moz_effxs.nc     moz_effxs.nc
         ln -sf ${input_basepath}/hammoz/moz_xs_short.nc  moz_xs_short.nc

         if ! $flag_fastj ; then
            ln -sf ${input_basepath}/hammoz/moz_rsf_gt200nm.nc                     moz_rsf_gt200nm.nc
            ln -sf ${input_basepath}/hammoz/moz_temp_prs_GT200nm_JPL10_c130206.nc  moz_temp_prs_gt200nm.nc
         else
            ln -sf ${input_basepath}/hammoz/moz_FJX_scat.dat   FJX_scat.dat
            ln -sf ${input_basepath}/hammoz/moz_FJX_spec.dat   FJX_spec.dat
            ln -sf ${input_basepath}/hammoz/moz_chem_Js.dat    chem_Js.dat
         fi
      fi

   fi

fi # end echam6-hammoz specific

#---------
# Nudging
#---------

if $flag_nudg ; then

   imo=10#$exp_start_month
   iyear=$exp_start_year
   idate=$exp_start_yyyymm
   previous_prefix="" # cosmetics. Just for nice output!
   if $flag_nudg_netcdf ; then
      declare -a exts=( "nc" )
      dir_ext="_netcdf"
      format="netcdf"
      date_in_subdir_name=false
      subdir_tail="_${hres}${vres}${dir_ext}" # <-- underscore at the beg like latest pattern by Sebastian
   else
      declare -a exts=( "div" "sst" "stp" "vor" )
      dir_ext=""
      format="binary"
      date_in_subdir_name=true
   fi

   #-- Initialization of the default preferred order for the nudging prefixes, if not already set by user
   if [ -z $nudg_prefix_order ] ; then 
      declare -a nudg_prefix_order=( "eraia" "era40"  "ana" )
   fi

   #-- Loop over YYYYMM:
   while ((idate<=$stop_yyyymm)) ; do

       if $date_in_subdir_name ; then
          subdir_tail="${hres}${vres}_${idate}${dir_ext}" # <-- no underscore at the beg (historical way)
       fi

       #-- Loop through all potential nudging data prefixes to find the first valid one
       flag_not_found=true
       for prefix in ${nudg_prefix_order[*]} ; do

           nudg_subdir_attempt=${nudg_basepath}/${prefix}/${prefix}${subdir_tail}

           if [ -d ${nudg_subdir_attempt} ] ; then
              flag_not_found=false
              nudg_subdir=$nudg_subdir_attempt
              break
           fi
       done

       if $flag_not_found ; then
          #-- Try any pattern (*) and pick up the first found by alphabetical order:
          nudg_subdir=`find ${nudg_basepath} -maxdepth 2 -type d \
                       -name "*${subdir_tail}" 2> /dev/null | sort | head -1`

          #-- Raise an error if nothing found:
          if [[ "$nudg_subdir" == "" ]] ; then
             echo "ERROR! Nudging data for ${hres}${vres} at $idate in $format format does not seem to be available in ${nudg_basepath}." | tee -a $log_file
             echo "You may check if alternate format exists for this date and change the NDGCTL namelist parameters accordingly." | tee -a $log_file
             echo "Exiting" | tee -a $log_file
             exit 1
          fi

          #-- Get prefix
          prefix=`basename $nudg_subdir`
          prefix=${prefix%%${subdir_tail}}
       fi

       #-- Prints out the prefix, but only when it changes to avoid cluttering the log:
       if [[ "$prefix" != "$previous_prefix" ]]  ; then
          echo "Nudging dataset from $idate on: $prefix" | tee -a $log_file
          previous_prefix=$prefix
       fi

       nudg_name="${prefix}${hres}${vres}_${idate}"

       #-- make the links:
       for ext in ${exts[*]} ; do
           ln -sf ${nudg_subdir}/${nudg_name}.${ext} ndg${idate}.${ext}
       done

       #-- additional linking for the edges 
       #   (fake nudging files for exp_start_yyyymm_m1 and stop_yyyymm_p1)
       #   Note: this is possible due to the fact that any nudging file contains
       #         the necessary timsteps BEFORE and AFTER the current month
       #         for echam to perform the correct interpolation

       if [[ "$idate" -eq "$exp_start_yyyymm" ]] ; then

          for ext in ${exts[*]} ; do
              ln -sf ${nudg_subdir}/${nudg_name}.${ext} ndg${exp_start_yyyymm_m1}.${ext}
          done

       elif [[ "$idate" -eq "$stop_yyyymm" ]] ; then

          for ext in ${exts[*]} ; do
              ln -sf ${nudg_subdir}/${nudg_name}.${ext} ndg${stop_yyyymm_p1}.${ext}
          done

       fi

       #-- increment idate (YYYYMM):
       ((imo=$imo+1))
       if [[ "$imo" -eq "13" ]] ; then
          imo=1
          ((iyear=$iyear+1))
       fi
       idate=${iyear}`printf '%02d' $imo`

   done

fi # end nudging

#--------------------------------
# Get back to previous directory
#--------------------------------

cd $OLDPWD

echo  | tee -a $log_file
