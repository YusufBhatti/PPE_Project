Author: Nick Schutgens (based on Mirek Andrejczuk's earlier work), n.a.j.schutgens@vu.nl


Information on how to install and run an ECHAM-HAM PPE

Directory structure:
my_experiments | run
  PPE
    PPE_defaults

Scripts (submitted to job queue):
PPE_install.sh  : installs directories for individual experiments that are part of PPE
PPE_batch.sh    : runs PPE in small batches of simultaneously run individual experiments

NOTE: further information can be found in the scripts themselves
NOTE: PPE_batch.sh can be run multiple times. Suppose ARCHER goes down after only part of the PPE has been run. Resubmission of PPE_batch.sh will lead to a continuation of the PPE. The PPE_values.txt file is used to keep track of which experiments have already been submitted and which not. The file does not contain information on succesful conclusion of an experiment (see also below)!!!

Additional files:
PPE_values.txt  : contains names of perturbed parameters and their perturbed values
emi_spec_PPE_defaults.txt  : emi_spec template for PPE, stored in PPE_defaults 
settings_PPE_defaults      : settings template for PPE, stored in PPE_defaults
symlinks_PPE_defaults.sh   : symlinks template for PPE, stored in PPE_defaults

NOTE: these additional files will need to be adapted to the sort of PPE a user wants to run


ISSUES:
1) Identical queue and account need to be specified by the user in both PPE_batch.sh and settings_PPE_defaults, ncpu in settings_PPE_defaults and PBS -l select in PPE_bathc.sh need to be consistent.

2) All run-time logging output from individual experiments is sent to a single file: PPE_log.txt

3) The value of END in PPE_values.txt currently serves no pruposes and will always be 0 (zero)




