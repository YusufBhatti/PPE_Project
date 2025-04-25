Template created by Ulrike Proske (ulrike.proske@env.ethz.ch) on 2023-11-22_07.59.

Specificities of this template:
-------------------------------
Template for the publication
_Developing a climatological simplification of aerosols to enter the cloud microphysics of a global climate model_
Proske et al., 2023
ACPD
Specifically, this is the template for the present-day default simulations with the Lin&Leaitch scheme.

This template also contains the post-processing scripts to produce CCN/INclim and AEROclim
----------------------------------------------------
 Scripts to create the ccnclim climatology
 Author: Sylvaine Ferrachat and Ulrike Proske
 Date: 28 02 2023
----------------------------------------------------

**USAGE**

`proc_echam-ham_as_ccnclim_3d_new.sh 2003_ham.nc 2003 2003-01-01 1month` produces 2003_ccn.nc and 2003_IN.nc plus more modified ccnclims if 2003_activ.nc is present

For *proc_echam-ham_as_ccnclim_3d_new_PDPIFUT.sh* on daint:
`conda activate conda_process_for-ice`
`module load NCO`
`~/code/echam_aerosol_scripts/CCN_clim/proc_echam-ham_as_ccnclim_3d_new_PDPIFUT.sh 2003-2012_ham.nc 2003-2012 2003-01-01 1month aeroclim_ll_default_05 aeroclim_ll_default_10years`
while in `/project/s1144/uproske/aeroclim_climatologies/aeroclim_ll_default_05`

*ATTENTION:* The reference date needs to be set to a first of January for the model to use it correctly (otherwise there's a mistake in the time variable and the model uses the last December file until the October of the next year)!
