Author: Yusuf Bhatti, (based on Duncan Watson-Paris PPE Scripts), y.bhatti@sron.nl

Information on how to install and run an ECHAM-HAM PPE

Directory structure:
Branches
    Perturbed_Climate_model
	my_experiments
    PPE_Scripts	

load necessary python packages such as pandas.

Scripts (submitted to job queue):
PPE_install.sh  : installs directories for individual experiments that are part of PPE
PPE_batch.sh    : runs PPE in small batches of simultaneously run individual experiments

NOTE: further information can be found in the scripts themselves
NOTE: PPE_batch.sh can be run multiple times. Suppose HPC goes down after only part of the PPE has been run. Resubmission of PPE_batch.sh will lead to a continuation of the PPE. The PPE_values.txt file is used to keep track of which experiments have already been submitted and which not. The file does not contain information on succesful conclusion of an experiment (see also below)!!!

Initialisation files:
Latin_Hypercube_Sampling.py : used to create LHC_Parameters.txt, which is the sampling used based on your min and max parameters and the parameters being perturbed. Need to modify simulations number. when running this file, you will be prompted to input the number of simulations you wish.

bash_script.py : python script to execute. This writes your LHC into the PPE_values file for final prep of P  PE. It also detects for if a parameter contains duplicate values. If there is a duplicate value detected, you should modify the number of decimal places it uses.

Steps:
1.
$ python Latin_Hypercube_Sampling.py
you will then be prompt to enter how many ensemble members you want to use. Enter the number
2.
$ python bash_script.py
Will write into your PPE_values.txt and fill in your parameter values with the values from your LHC.
3.
Sort out directories within the PPE_install.sh and PPE_batch.sh
4.
qsub PPE_install.sh
5.
qsub PPE_batch.sh


NOTES: The Initialisation files rely on the parameters_for_script.txt file. The hypercube sampling script uses the parameters_for_script.txt to take in how many columns there will be (how many parameters will be perturbed in your PPE). You specify in the script (n_samples) to the amount of simulations you want. Typically Nu. parameters * 5-7. This will output an array (n_samples, number of parameters) called LHC_Parameters.txt. 

Additional files:
emi_spec_PPE_defaults.txt  : emi_spec template for PPE, stored in PPE_defaults 
settings_PPE_defaults      : settings template for PPE, stored in PPE_defaults
symlinks_PPE_defaults.sh   : symlinks template for PPE, stored in PPE_defaults
parameters_for_script.txt  : Parameter name AND the values which are scaling factors. This is core file which is needing to be modified based on PPE specified. The left value is minimum range, right value is maximum value.

NOTE: these additional files will need to be adapted to the sort of PPE a user wants to run

Files:
PPE_values.txt  : contains names of perturbed parameters and their perturbed values. This will be changed relative to the parameters within parameters_for_script.txt.


in parameters_for_script.txt, the name needs to be identical to the top row in PPE_values.txt. The order of them DO NOT MATTER. They will automatically fill in the necessary areas relative to the scaling factor provided
To fill in the values based on the values you are perturbing
run $python bash_values.py 
This will read the parameters_for_script.txt and use them values. The script is set up so that the amount of scale factors you use for one parameter will duplicate this variable replacing the scale factor. So you can have 3 scales for 1 variable but 5 scales for 5 variables. This will duplicate the 5 variables 5 times (5*5) and the 1 variable * 3. 

ISSUES:

1) Identical queue and account need to be specified by the user in both PPE_batch.sh and settings_PPE_defaults, ncpu in settings_PPE_defaults and PBS -l select in PPE_bathc.sh need to be consistent.

2) All run-time logging output from individual experiments is sent to a single file: PPE_log.txt

3) The value of END in PPE_values.txt currently serves no pruposes and will always be 0 (zero)




