#!/bin/bash
# this script allows you to reset the PPE_values file start, install, and end to 0.
# this is very useful for debugging and testing your PPE script. 
# to run, simply . fix_values.sh
#


# Define the input file and the output file
input_file="PPE_values.txt"
output_file="PPE_valuei_new.txt"

# Use awk to process the file
awk 'NR==1 { print; next } { $(NF-2) = 0; $(NF-1) = 0; $NF = 0; print }' $input_file > $output_file
#'NR==1 { print; next }: This ensures that the first line (record number 1) is printed as is and then skips to the next line without processing.
#{ $(NF-2) = 0; $(NF-1) = 0; $NF = 0; print }: This processes each subsequent line, setting the last three fields to 0 and then printing the modified line.

mv $output_file $input_file

