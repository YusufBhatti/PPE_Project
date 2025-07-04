import numpy as np
import pandas as pd
from pyDOE import *
from scipy.stats import qmc
import sys
import os
cwd = os.getcwd()
utils_path = os.path.join(cwd, 'Scripts_for_config')
sys.path.append(utils_path)

from utils import create_param_ranges
from utils import PPE_values_modify
from utils import txt_to_csv
from scipy.stats import lognorm

os.chdir(cwd)


'''

Code Owner: Yusuf Bhatti (y.bhatti@sron.nl)

This file first edits your PPE_values_template.txt to ensure the headers match your parameters list.

We then read your parameter list, establish what each parameter is called, and its minimum value (left) and maximum value (right) as a list.

you then state the number of simulations you plan of doing (n_simulations). 

We then do a minmax LHC sampling across your n_simulations, which are distributed randomly across the column for each parameter, respectively. This allows you to run only the first 80% simulations to ensure you have the computational resources available for this. The LHC is an array (number of simulations x number of parameters)

This is then saved as a .txt. file. We also save a csv file for ease of analysis.

'''


input_file = f"{cwd}/PPE_values_template.txt"
parameters_file = f"{cwd}/parameters_for_script.txt"





PPE_values_modify(input_file, parameters_file)



parameters_file = f"{cwd}/parameters_for_script.txt"
Parameters_and_ranges = create_param_ranges(parameters_file)
param_ranges = Parameters_and_ranges.values()
# Number of samples
# --- Step 1: Print all available parameter names ---
param_names = list(Parameters_and_ranges.keys())
print("Available parameters:\n" + "\n".join(param_names))

# --- Ask user which parameters to treat with log-normal distribution ---
custom_dist_params = input(
    "\nEnter parameter(s) to use log-normal distribution (space-separated, e.g., 'NUC_FT EMI_DMS'): "
).split()


if __name__ == "__main__":
    try:
        # Prompt user for the number of simulations
        n_simulations = int(input("Enter number of simulations you wish to perform: "))
    except ValueError:
        print("Invalid input. Please enter an integer value.")
        exit(1)
    
#n_simulations = 200

lower = []
upper = []
names = []

for i, (low, high) in enumerate(param_ranges):
    lower.append(low)
    upper.append(high)
for i, name in enumerate(Parameters_and_ranges.items()):
    names.append(name[0])

l_bound = np.array(lower)
u_bound = np.array(upper)
names = np.array(names)

# sampler = qmc.LatinHypercube(d=23)
# sample = sampler.random(n=150)

lhs_sample = lhs(len(param_ranges), samples=n_simulations, criterion='maximin')

scale = qmc.scale(lhs_sample, l_bound, u_bound)


"""
I will make the variables below input values, so you can manually identify which parameters should be log-normal, vs uniform distributions, for the emulator 
"""

# # Identify the index for NUC_FT
# nuc_index = names.tolist().index('NUC_FT')

# # Parameters for log-normal (these are tunable):
# # This sets the mode near 1
# s = 2  # shape (σ)
# scale_param = 1  # scale = exp(μ) (log-normal distribution directly related to the median)
# dist = lognorm(s=s, scale=scale_param)

# # Rescale uniform LHS samples for NUC_FT via inverse CDF (ppf)
# # This transforms the uniform LHC sample to log-normal
# uniform_samples = lhs_sample[:, nuc_index]
# nuc_samples = dist.ppf(uniform_samples)

# # Clip values to [min, max] for physical validity of the parameters which need this.
# nuc_min, nuc_max = Parameters_and_ranges['NUC_FT']
# nuc_samples = np.clip(nuc_samples, nuc_min, nuc_max)

# # Replace in the final scaled LHC array
# scale[:, nuc_index] = nuc_samples
# --- Step 4: Apply log-normal transformation to selected parameters ---
for param in custom_dist_params:
    if param not in Parameters_and_ranges:
        print(f"Warning: '{param}' not found in parameter list. Skipping.")
        continue

    idx = names.tolist().index(param)
    low, high = Parameters_and_ranges[param]

    # Choose shape and scale for log-normal (tweakable)
    s = 1
    scale_param = 1  # Median will be at 1
    dist = lognorm(s=s, scale=scale_param)

    # Transform uniform LHC samples to log-normal
    transformed = dist.ppf(lhs_sample[:, idx])
    transformed = np.clip(transformed, low, high)

    # Replace column in LHC matrix
    scale[:, idx] = transformed




print('Your LHC has been saved as a .txt file and as a .csv file')
print(f'the LHC has {np.shape(scale)[0]} simulations and {np.shape(scale)[1]} parameters')
np.savetxt(f"{cwd}/parameter_values_data/LHC_Parameters.txt", scale, delimiter=" ")
txt_to_csv(f"{cwd}/parameter_values_data/LHC_Parameters.txt", f"{cwd}/parameter_values_data/LHC_Parameters.csv")
