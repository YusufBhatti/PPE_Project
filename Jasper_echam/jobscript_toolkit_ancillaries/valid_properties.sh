#!/bin/bash

#--------------------------------------------------------------------------------------------------
#
# Sylvaine Ferrachat 2016-01
#
# This script defines the valid model types and their properties, as well as the valid emissions
#
#--------------------------------------------------------------------------------------------------

#------------------------------------------------
#-- Valid model types 
#   (must not contain any ':' nor '=')

    e6="echam6"
    e6_comment='ECHAM6'

    e6_ham="echam6-ham"
    e6_ham_comment='ECHAM6-HAMMOZ, with only HAM switched on. This is the default model type.'

    e6_moz="echam6-moz"
    e6_moz_comment='ECHAM6-HAMMOZ, with only MOZ switched on.'

    e6_hammoz="echam6-hammoz"
    e6_hammoz_comment='echam6-hammoz, with both HAM and MOZ switched on.'

#------------------------------------------------
#-- Set model properties
#
#   properties construct (each property but the last one MUST be terminated with ':'):  
#    - model_name
#    - model_comment: any string that describes briefly the model
#    - model_suffix (currently mostly useful for fetching the proper emi_spec file)
#    - model_echam_type (0, 5, or 6) 0 means this is not applicable, (e.g. M7 box model)
#    - model_flag_with_hammoz (false/true) : false for echam, true for echam5.5-ham, echam6-ham, etc...
#    - model_flag_emi (false/true) : true if the model needs an emi_spec file

# WARNING: key names must start with 'model_' !!!

    e6_props="model_name=${e6}:\
              model_comment=${e6_comment}:\
              model_suffix=:\
              model_echam_type=6:\
              model_flag_with_hammoz=false:\
              model_flag_emi=false:"
    e6_props=`echo $e6_props | sed -e 's/: */:/g'`

    e6_ham_props="model_name=${e6_ham}:\
                  model_comment=${e6_ham_comment}:\
                  model_suffix=_ham:\
                  model_echam_type=6:\
                  model_flag_with_hammoz=true:\
                  model_flag_emi=true:"
    e6_ham_props=`echo $e6_ham_props | sed -e 's/: */:/g'`

    e6_moz_props="model_name=${e6_moz}:\
                  model_comment=${e6_moz_comment}:\
                  model_suffix=_moz:\
                  model_echam_type=6:\
                  model_flag_with_hammoz=true:\
                  model_flag_emi=true:"
    e6_moz_props=`echo $e6_moz_props | sed -e 's/: */:/g'`

    e6_hammoz_props="model_name=${e6_hammoz}:\
                     model_comment=${e6_hammoz_comment}:\
                     model_suffix=:\
                     model_echam_type=6:\
                     model_flag_with_hammoz=true:\
                     model_flag_emi=true:"
    e6_hammoz_props=`echo $e6_hammoz_props | sed -e 's/: */:/g'`

    model_array=( \
                "$e6_props"\
                "$e6_ham_props"
                "$e6_moz_props"\
                "$e6_hammoz_props"\
                )

    n_model=${#model_array[*]}
