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

import numpy as np
import matplotlib.pyplot as plt

def trapezoidal_pdf(x, min_val, left, right, max_val, aL=2, aR=2, slope=0):
    """
    Compute trapezoidal PDF for given input `x` values.
    
    Parameters:
        x : np.ndarray
            Points at which to evaluate the PDF
        min_val : float
            Minimum value of the support (a)
        left : float
            Start of the flat top (b)
        right : float
            End of the flat top (c)
        max_val : float
            Maximum value of the support (d)
        aL : float
            Shape of the left slope (convex/concave)
        aR : float
            Shape of the right slope (convex/concave)
        slope : float
            Optional slope between the flat top
        
    Returns:
        pdf : np.ndarray
            Probability density function evaluated at x
    """
    pdf = np.zeros_like(x)

    # Normalize coefficients to ensure integral = 1
    def area():
        left_area = (left - min_val) / (aL + 1)
        middle_area = (right - left)
        right_area = (max_val - right) / (aR + 1)
        return left_area + middle_area + right_area

    norm = 1 / area()

    for i, xi in enumerate(x):
        if min_val <= xi < left:
            # Left slope
            u = (xi - min_val) / (left - min_val)
            pdf[i] = norm * u**aL
        elif left <= xi <= right:
            # Flat top or sloped top
            if slope == 0:
                pdf[i] = norm
            else:
                top_len = right - left
                pdf[i] = norm * (1 + slope * (xi - (left + right) / 2) / top_len)
        elif right < xi <= max_val:
            # Right slope
            u = (max_val - xi) / (max_val - right)
            pdf[i] = norm * u**aR
        else:
            pdf[i] = 0
    return pdf

def sample_from_trapezoidal(n_samples, min_val, left, right, max_val, aL=2, aR=2, slope=0):
    from scipy.interpolate import interp1d
    from scipy.integrate import cumulative_trapezoid

    # Step 1: Create fine x grid
    x = np.linspace(min_val, max_val, 10000)
    
    # Step 2: Evaluate PDF
    pdf = trapezoidal_pdf(x, min_val, left, right, max_val, aL, aR, slope)
    
    # Step 3: Compute normalized CDF
    cdf = cumulative_trapezoid(pdf, x, initial=0)
    cdf /= cdf[-1]  # Normalize to [0,1]

    # Step 4: Invert the CDF
    inverse_cdf = interp1d(cdf, x, bounds_error=False, fill_value=(min_val, max_val))

    # # Step 5: Sample uniform values and transform through inverse CDF
    # uniform_samples = np.random.uniform(0, 1, n_samples)
    # samples = inverse_cdf(uniform_samples)
    # LHS sampling in 1D on [0,1]
    lhs_samples = lhs(1, samples=n_samples, criterion='maximin').flatten()
    
    # Transform LHS samples through inverse CDF
    samples = inverse_cdf(lhs_samples)

    return samples


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
print('--')

# --- Ask user which parameters to treat with log-normal distribution ---
# custom_dist_params = input(
#     "\nEnter parameter(s) to use log-normal distribution (space-separated, e.g., 'NUC_FT EMI_DMS'): "
# ).split()
# print('--')


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

# # --- Step 4: Apply log-normal transformation to selected parameters ---
# for param in custom_dist_params:
#     if param not in Parameters_and_ranges:
#         print(f"Warning: '{param}' not found in parameter list. Skipping.")
#         continue

#     idx = names.tolist().index(param)
#     low, high = Parameters_and_ranges[param]

#     # Choose shape and scale for log-normal (tweakable)
#     s = 10
#     scale_param = 2  # Median will be at 1
#     dist = lognorm(s=s, scale=scale_param)

#     # Transform uniform LHC samples to log-normal
#     transformed = dist.ppf(lhs_sample[:, idx])
#     transformed = np.clip(transformed, low, high)

#     # Replace column in LHC matrix
#     scale[:, idx] = transformed

# Predefined trapezoidal shape settings for specific parameters:
# Format: param_name: (vertex1, vertex2, aL, aR)
trapezoidal_param_config = {
    'SO4_COATING':      (0.5, 2.5, 2, 2),  # monolayers, absolute
    'NUC_FT':           (0.1, 5, 2, 2),
    'DRYDEP_ACC':       (0.7, 4, 2, 2),
    'SCALE_EMI_BB_SO2': (0.5, 2, 2, 2),
    'SCALE_EMI_BB_BC':  (0.5, 2, 2, 2),
    'SCALE_EMI_BB_OC':  (0.5, 2, 2, 2)
}

# --- Step 4: Apply Trapezoidal transformation using inverse CDF sampling ---
# --- Apply trapezoidal transformation ---
#for param in custom_dist_params:
for param, (left, right, aL, aR) in trapezoidal_param_config.items():
    if param not in Parameters_and_ranges:
        print(f"Warning: '{param}' not found in parameter list. Skipping.")
        continue

    if param not in trapezoidal_param_config:
        print(f"Skipping trapezoidal distribution for '{param}'. Using uniform.")
        continue
    print('--')
        
    idx = names.tolist().index(param)
    a, d = Parameters_and_ranges[param]  # a = min_val, d = max_val

    left, right, aL, aR = trapezoidal_param_config[param]

    print(f"Applying trapezoidal distribution to '{param}' with:")
    print(f"  Vertex 1: {left}, Vertex 2: {right}, aL: {aL}, aR: {aR} | Range: ({a}, {d})")

    # Define inner plateau (can be customized)
    range_width = d - a
    slope = 0

    # Generate Trapezoidal samples using inverse CDF
    samples = sample_from_trapezoidal(
        n_samples=n_simulations,
        min_val=a,
        left=left,
        right=right,
        max_val=d,
        aL=aL,
        aR=aR,
        slope=slope
    )

    # Replace column in LHC
    scale[:, idx] = samples




print('Your LHC has been saved as a .txt file and as a .csv file')
print(f'the LHC has {np.shape(scale)[0]} simulations and {np.shape(scale)[1]} parameters')
# np.savetxt(f"{cwd}/parameter_values_data/LHC_Parameters.txt", scale, delimiter=" ")
# txt_to_csv(f"{cwd}/parameter_values_data/LHC_Parameters.txt", f"{cwd}/parameter_values_data/LHC_Parameters.csv")
