#!/bin/bash
#-------------------------------------------------------------------------------
# Sylvaine Ferrachat 2013-03
#
# This piece of script is read by the jobscript_toolkit scripts, in order to get
# important informations about the echam-hammoz version history so that
# backward compatibility is ensured.
#
# Please do not trash nor modify this file, otherwise it may yield 
# unexpected behaviours in the preparation of your run.
#
#-------------------------------------------------------------------------------

# Keep track of the new vs old setup regarding to emi_basepath definition in echam6-hammoz:
flag_old_emi_basepath=false

# Keep track of the new vs old biogenic emissions flag (MEGAN-related):
flag_old_bioemi=false

# Check for obsolete jobscript toolkit version:
minimum_version="1.2"

# Set flag for new naming convention (HAM):
new_naming_ham=true

# Set prefix and suffix for rerun files (new in echam6.3)
prefix_rerun_file="restart"
suffix_rerun_file=".nc"
