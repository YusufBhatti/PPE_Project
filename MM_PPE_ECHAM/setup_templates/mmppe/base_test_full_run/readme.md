**EMI**

CEDS: AGR and SLV are deactivated, since it seems that they don't include any aerosol emissions relevant for HAM-M7.
todo> here, a check is needed

CMIP: AWB is deactivated since we want to avoid double counting with GFED! 

**SYMLINKS/SETTINGS**

ERA5 for sst/sic is included. (However only 2025 data exists at the end of this year they might be some problem).

Scenario is fixed for ozon (vl) and ghg (m). Based on the setting file
, we used rcp45. However, the scenario is applied for the kinne climatology, %R0 in emi-matrix and MOZ lower boundary conditions. However based on our configuration all of these seems to be not activated at all. 
todo> check again if there aren't activated.

For the symlinks it . 