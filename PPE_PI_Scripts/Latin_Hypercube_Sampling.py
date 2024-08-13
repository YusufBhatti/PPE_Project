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

os.chdir(cwd)


'''
This file first edits your PPE_values.txt to ensure the headers match your parameters list.

We then read your parameter list, establish what each parameter is called, and its minimum value (left) and maximum value (right) as a list.

you then state the number of simulations you plan of doing (n_simulations). 

We then do a minmax LHC sampling across your n_simulations, which are distributed randomly across the column for each parameter, respectively. This allows you to run only the first 80% simulations to ensure you have the computational resources available for this. The LHC is an array (number of simulations x number of parameters)

This is then saved as a .txt. file. We also save a csv file for ease of analysis.

'''


input_file = f"{cwd}/PPE_values.txt"
parameters_file = f"{cwd}/parameters_for_script.txt"





PPE_values_modify(input_file, parameters_file)



parameters_file = f"{cwd}/parameters_for_script.txt"
Parameters_and_ranges = create_param_ranges(parameters_file)
param_ranges = Parameters_and_ranges.values()
# Number of samples

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
print('Your LHC has been saved as a .txt file and as a .csv file')
print(f'the LHC has {np.shape(scale)[0]} simulations and {np.shape(scale)[1]} parameters')
np.savetxt(f"{cwd}/parameter_values_data/LHC_Parameters.txt", scale, delimiter=" ")
txt_to_csv(f"{cwd}/parameter_values_data/LHC_Parameters.txt", f"{cwd}/parameter_values_data/LHC_Parameters.csv")
