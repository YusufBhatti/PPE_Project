#!/bin/bash

input_file="PPE_values.txt"
output_file="PPE_values_new.txt"

# Read the first line (header)
IFS= read -r header < "$input_file"

# Write the header to the new output file
echo "$header" > "$output_file"

# Convert the header into an array
header_array=($header)

# Read the rest of the lines, starting from the second line
tail -n +2 "$input_file" | while IFS= read -r line; do
    # Read the parameter values into an array
    param_values=($line)

    # Iterate over the header variables (starting from index 1 to skip the first column)
    for ((i = 1; i < ${#header_array[@]} - 3; i++)); do
        # Create a new experiment name based on the header variable
        new_exp_name="TEST_${header_array[i]}"

        # Copy the current param_values array
        new_param_values=("${param_values[@]}")

        # Set the corresponding value to 4 for the current variable
        new_param_values[i]=4

        # Set all other values after the first three variables to 1 if they are not already 4
        for ((j = 4; j < ${#new_param_values[@]} - 3; j++)); do
            if [ "${new_param_values[j]}" -ne 4 ]; then
                new_param_values[j]=1
            fi
        done

        # Construct the new line with the updated values
        new_line="$new_exp_name"
        for ((j = 1; j < ${#new_param_values[@]}; j++)); do
            new_line+=" ${new_param_values[j]}"
        done

        # Append the new experiment line to the output file
        echo "$new_line" >> "$output_file"

        # Stop the iteration if the new experiment name is TEST_V_SCALE_WETDEP_IC_BC_ONLY
        if [ "$new_exp_name" == "TEST_V_SCALE_WETDEP_IC_BC_ONLY" ]; then
            exit 0
        fi
    done
done

