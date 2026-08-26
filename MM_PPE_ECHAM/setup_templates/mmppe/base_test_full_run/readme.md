**IMPORTANT**
This works only for L47 model since in the nudging file we hardcoded the highest level

**EMI**

CEDS: AGR and SLV are deactivated, since it seems that they don't include any aerosol emissions relevant for HAM-M7.

> a check with claude revealed no emissions (BC, OC, SO2) at all from AGR and SLV

CMIP: AWB is deactivated since we want to avoid double counting with GFED! 

**SYMLINKS/SETTINGS**

ERA5 for sst/sic is included. (However only 2025 data exists at the end of this year they might be some problem).

NDG: First 40 layers (until 2km) are nudged, in 6h steps. Surface pressure is also nudged implicitly since nothing was set, it's like (24 hourly):
 nudgp           = 1.1574           ! lnps

Scenario is fixed for ozon (vl) and ghg (m). Based on the setting file, we used rcp45. However, the scenario is applied for the kinne climatology, %R0 in emi-matrix and MOZ lower boundary conditions. However based on our configuration all of these seems to be not activated at all. 
todo (not of great importance) > check again if there aren't activated.