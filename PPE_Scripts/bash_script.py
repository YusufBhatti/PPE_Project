
input_file = "PPE_values.txt"
output_file = "PPE_values_new.txt"
from subprocess import run
from os import rename
from re import search
import numpy as np
parameters_file = "LHC_Parameters.txt"
with open(parameters_file, 'r') as file:
    column_count = len(file.readline().strip().split())
    line_count = sum(1 for line in file)
Array = np.zeros((line_count+1, column_count)); Array[:]=np.nan
print(f'The script will use {line_count} simulations.')
print(f'The PPE contains {column_count} parameters.')

parameters = {}
with open(parameters_file, 'r') as f_params:
    for i, line in enumerate(f_params):
        #print(count)

        parts = line.strip().split()
      #  if header_array == parts[0]:
        parameter_name = 'TEST'
        scaling_factors = list(map(float, parts[:]))
       # for i in range(0, line_count):
        Array[i,:] = scaling_factors


# Read the first line (header)
with open(input_file, 'r') as f_in:
    header = f_in.readline().strip()

# Write the header to the new output file
with open(output_file, 'w') as f_out:
    with open(parameters_file, 'r') as f_params:
        f_out.write(header + "\n")
        # Convert the header into an array
        header_array = header.split()
        # Read the rest of the lines, starting from the second line
        with open(input_file, 'r') as f_in:
            for i, line in enumerate(f_params):
                parts = line.strip().split()
                param_values = line.strip().split()
                # Create a new experiment name based on the header variable
                new_exp_name = f"PPE_TEST_{i+1}"
                # Copy the current param_values array
               # new_param_values = param_values[:]
                # Set the corresponding value to 4 for the current variable
                new_param_values = [f"{float(val):.4e}" if float(val) < 0.1 else f"{float(val):.4f}" for val in param_values]

                # Set last 3 values to 1 
                new_param_values.extend(['0'] * 3)
                # Include the name of the simulation here
                new_param_values.insert(0, new_exp_name)

                # Ensure the new_param_values array has the same length as the header_array
                if len(new_param_values) < len(header_array):
                    print("STOP Your number of parameters does not match up TOO LITTLE with" 
                          "  the ppe_values.txt file header!. I am doing a workaround" 
                          "  for you and implementing a 'NAN' to make it match"
                          "  This is not a workable copy, only a runable one!")

                    while len(new_param_values) < len(header_array):
                        new_param_values.append("NAN")
                elif len(new_param_values) > len(header_array):
                    print("STOP Your number of parameters is TOO much than the"
                          "  the ppe_values.txt file header!. I am doing a workaround"
                          "  for you and REDUCING your parameter values to FORCE it match. "
                          "  This is not a workable copy, only a runable one!")
                    new_param_values = new_param_values[:len(header_array)]

                # Construct the new line with the updated values
                new_line = " ".join([new_exp_name] + new_param_values[1:])
                print(new_line)
                # Append the new experiment line to the output file
                f_out.write(new_line + "\n")
                #end
                # Stop the iteration if the new experiment name is TEST_V_SCALE_WATER
                #if new_exp_name == f"*{header_array[-4]}":
        #         if search(f".*{header_array[-4]}$", new_exp_name):

        #             break
            
        # break

rename(output_file, input_file)
run(["bash", "fix_values.sh"], check=True)
print(f"New PPE values generated and saved to {input_file}")

def check_column_duplicates(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Read the data into a 2D list
    data = [line.strip().split() for line in lines]
    
    # Transpose the data to work with columns
    columns = list(zip(*data))
    colum = len(line.strip().split())

    for col_index, column in enumerate(columns):
        duplicates = set([val for val in column if column.count(val) > 1])
        if duplicates:
            if col_index > colum:
                pass
            else:
                print(f"Duplicate values found in {header_array[col_index]}: {duplicates}")
                print(f"STOP AND CHANGE LINE 55 in the bash_script.py to increased decimal points.")

        else:
            print(f"No duplicate values found in {header_array[col_index]} columns")
            pass

# Replace 'parameters_file' with the path to your .txt file

print(f"Now we will check if any rows contain duplicates in the PPE values generated")

parameters_file = "/home/ybhatti/yusufb/Branches/PPE_Scripts/PPE_values.txt"
check_column_duplicates(parameters_file)


#############################
""" This following script is just when you want all values to be '1' except for the perturbed parameter"""
#############################

# #%time

# # Python script to modify PPE_values.txt based on scaling factors from parameters_for_script.txt
# def Parameter_values(header_array):
#     parameters_file = "parameters_for_script.txt"
#     # Read parameters from parameters_for_script.txt into a dictionary
#     scale=[]
#     parameters = {}
#     with open(parameters_file, 'r') as f_params:
#         for line in f_params:
#             parts = line.strip().split()
#             if header_array == parts[0]:
#                 parameter_name = parts[0]
#                 scaling_factors = list(map(float, parts[1:]))
#                 length = len(scaling_factors)
#                 for i in range(length):
#                     scale.append(scaling_factors[i])
#                 #scale = scaling_factors[0]

#     return scale
    
# input_file = "PPE_values.txt"
# output_file = "PPE_values_new.txt"
# from subprocess import run
# from os import rename
# from re import search

# # Read the first line (header)
# with open(input_file, 'r') as f_in:
#     header = f_in.readline().strip()

# # Write the header to the new output file
# with open(output_file, 'w') as f_out:
#     f_out.write(header + "\n")

#     # Convert the header into an array
#     header_array = header.split()

#     # Read the rest of the lines, starting from the second line
#     with open(input_file, 'r') as f_in:
#        # next(f_in)  # Skip the header line
#         for line in f_in:
#             print(line)
#             param_values = line.strip().split()

#             # Set the first three values 
#             """ Commented out as this is no longer needed
#             #param_values[1] = "1"
#             # param_values[2] = "0.17"
#             # param_values[3] = "0.08"
#             """
#             # Iterate over the header variables (starting from index 1)
#             for i in range(1, len(header_array)):
#                 scale_factor = Parameter_values(header_array[i])
#                 for r in range(0,len(scale_factor)):
#                     # Create a new experiment name based on the header variable
#                     new_exp_name = f"PPE_{header_array[i]}_{r+1}"

#                     # Copy the current param_values array
#                     new_param_values = param_values[:]
#                     # Set the corresponding value to 4 for the current variable
#                     new_param_values[i] = f"{scale_factor[r]}"
                    
#                     # Set all other values after the first three variables to 1 if they are not already 4
#                     for j in range(1, len(new_param_values)):
#                         if j != i:
#                             if j !=i:
#                                 new_param_values[j] = "1"
    
#                     # Ensure the new_param_values array has the same length as the header_array
#                     if len(new_param_values) < len(header_array):
#                         while len(new_param_values) < len(header_array):
#                             new_param_values.append("1")
#                     elif len(new_param_values) > len(header_array):
#                         new_param_values = new_param_values[:len(header_array)]
    
#                     # Construct the new line with the updated values
#                     new_line = " ".join([new_exp_name] + new_param_values[1:])
#                     print(new_line)
#                     # Append the new experiment line to the output file
#                     f_out.write(new_line + "\n")
#                     #end
#                     # Stop the iteration if the new experiment name is TEST_V_SCALE_WATER
#                     #if new_exp_name == f"*{header_array[-4]}":
#                     if search(f".*{header_array[-4]}$", new_exp_name):

#                         break
                
#             break

# rename(output_file, input_file)
# run(["bash", "fix_values.sh"], check=True)
# print(f"New PPE values generated and saved to {input_file}")
