#!/bin/bash

#--------------------------------------------------------------------------------------------------
#
# Sylvaine Ferrachat 2016-03
#
# This script defines the list of files that should be copied from script dir to exp dir
#
#--------------------------------------------------------------------------------------------------

list_files_to_copy+=( \
    "${script_dir}/namelist_${exp}.echam | namelist.echam | true"
    "${script_dir}/namelist_${exp}.jsbach | namelist.jsbach | true"
                    )

for namelist_file in ${script_dir}/*.nml ; do
    list_files_to_copy+=( "${namelist_file} | . | true" )
done

if $this_model_flag_emi ; then # emi_spec is relevant
   list_files_to_copy+=( "${emissions_spec} | emi_spec.txt | true" )
fi
