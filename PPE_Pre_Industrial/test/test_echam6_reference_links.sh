INI_ECHAM=/pool/data/ECHAM6/
#-----------------------------------------------------------------------
SCR_DIR=$1
HRES=$2
VRES=$3
ORES=$4
YEAR=$5
YEARM1=$(( YEAR - 1 ))
YEARP1=$(( YEAR + 1 ))
YEARP2=$(( YEAR + 2 ))
#-----------------------------------------------------------------------
#----------
# ECHAM6  
#----------
T_DIR=${INI_ECHAM}/T${HRES}
#
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_OZONE_cmip5_clim.nc          unit.21
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}L${VRES}_jan_spec.nc          unit.23
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}${ORES}_jan_surf.nc           unit.24
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}${ORES}_VLTCLIM.nc            unit.90
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}${ORES}_VGRATCLIM.nc          unit.91
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_TSLCLIM2.nc                  unit.92
T_DIR=${INI_ECHAM}
#
${SCR_DIR}/link.sh ${T_DIR}/ECHAM6_CldOptProps.nc  ECHAM6_CldOptProps.nc
${SCR_DIR}/link.sh ${T_DIR}/rrtmg_lw.nc            rrtmg_lw.nc
${SCR_DIR}/link.sh ${T_DIR}/rrtmg_sw.nc            rrtmg_sw.nc
#----------
# JS-BACH  
#----------
T_DIR=/pool/data/JSBACH/input/r0004/T${HRES}
#
${SCR_DIR}/link.sh  ${T_DIR}/jsbach_T${HRES}${ORES}_11tiles_5layers_1976.nc jsbach.nc
#
T_DIR=${INI_ECHAM}/jsbach
#
${SCR_DIR}/link.sh  ${T_DIR}/lctlib_nlct21.def_rev6931          lctlib.def
#${SCR_DIR}/link.sh  /scratch/local1/m218036/data/echam6_data/input/input_echam6/jsbach/lctlib_nlct21.def lctlib.def
#----------
# AMIP2
#----------
T_DIR=${INI_ECHAM}/T${HRES}/amip2
#
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_amip2sst_${YEARM1}.nc sst${YEARM1}
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_amip2sst_${YEAR}.nc   sst${YEAR}
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_amip2sst_${YEARP1}.nc sst${YEARP1}
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_amip2sst_${YEARP2}.nc sst${YEARP2}
#
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_amip2sic_${YEARM1}.nc ice${YEARM1}
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_amip2sic_${YEAR}.nc   ice${YEAR}
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_amip2sic_${YEARP1}.nc ice${YEARP1}
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_amip2sic_${YEARP2}.nc ice${YEARP2}
#----------------
# Kinne aerosols
#----------------
T_DIR=${INI_ECHAM}/T${HRES}/aero
#
${SCR_DIR}/link.sh ${T_DIR}/g30_coa_T${HRES}.nc           aero_coarse_${YEARM1}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_coa_T${HRES}.nc           aero_coarse_${YEAR}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_coa_T${HRES}.nc           aero_coarse_${YEARP1}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_coa_T${HRES}.nc           aero_coarse_${YEARP2}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_fin_T${HRES}_${YEARM1}.nc aero_fine_${YEARM1}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_fin_T${HRES}_${YEAR}.nc   aero_fine_${YEAR}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_fin_T${HRES}_${YEARP1}.nc aero_fine_${YEARP1}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_fin_T${HRES}_${YEARP2}.nc aero_fine_${YEARP2}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_fir_T${HRES}.nc           aero_farir_${YEARM1}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_fir_T${HRES}.nc           aero_farir_${YEAR}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_fir_T${HRES}.nc           aero_farir_${YEARP1}.nc
${SCR_DIR}/link.sh ${T_DIR}/g30_fir_T${HRES}.nc           aero_farir_${YEARP2}.nc
#-------------------------------
# Volcanic aerosols (Stenchikov)
#-------------------------------
T_DIR=${INI_ECHAM}/T${HRES}/volcano_aerosols
#
${SCR_DIR}/link.sh ${T_DIR}/strat_aerosol_sw_T${HRES}_1991.nc strat_aerosol_sw_${YEARM1}.nc
${SCR_DIR}/link.sh ${T_DIR}/strat_aerosol_sw_T${HRES}_1991.nc strat_aerosol_sw_${YEAR}.nc
${SCR_DIR}/link.sh ${T_DIR}/strat_aerosol_sw_T${HRES}_1991.nc strat_aerosol_sw_${YEARP1}.nc
${SCR_DIR}/link.sh ${T_DIR}/strat_aerosol_sw_T${HRES}_1991.nc strat_aerosol_sw_${YEARP2}.nc
${SCR_DIR}/link.sh ${T_DIR}/strat_aerosol_ir_T${HRES}_1991.nc strat_aerosol_ir_${YEARM1}.nc
${SCR_DIR}/link.sh ${T_DIR}/strat_aerosol_ir_T${HRES}_1991.nc strat_aerosol_ir_${YEAR}.nc
${SCR_DIR}/link.sh ${T_DIR}/strat_aerosol_ir_T${HRES}_1991.nc strat_aerosol_ir_${YEARP1}.nc
${SCR_DIR}/link.sh ${T_DIR}/strat_aerosol_ir_T${HRES}_1991.nc strat_aerosol_ir_${YEARP2}.nc
#-------------------------------
# Volcanic aerosols (predefined basic tables for HAM and Crowley aerosols)
#-------------------------------
#${SCR_DIR}/link.sh ${INI_ECHAM}/b30w120 aero_volc_tables.dat
#-------------------------------
# Volcanic aerosols based on HAM simulation
#-------------------------------
#${SCR_DIR}/link.sh /scratch/local1/m218036/data/echam6_data/input/input_echam6/volcano_tables/8mt_t31l39_zm_mm_aod_ham.nc aoddz_ham_${YEARM1}.nc
#${SCR_DIR}/link.sh /scratch/local1/m218036/data/echam6_data/input/input_echam6/volcano_tables/8mt_t31l39_zm_mm_aod_ham.nc aoddz_ham_${YEAR}.nc
#${SCR_DIR}/link.sh /scratch/local1/m218036/data/echam6_data/input/input_echam6/volcano_tables/8mt_t31l39_zm_mm_aod_ham.nc aoddz_ham_${YEARP1}.nc
#${SCR_DIR}/link.sh /scratch/local1/m218036/data/echam6_data/input/input_echam6/volcano_tables/8mt_t31l39_zm_mm_aod_ham.nc aoddz_ham_${YEARP2}.nc
#-------------------------------
# Volcanic aerosols Crowley data set
#-------------------------------
#${SCR_DIR}/link.sh /scratch/local1/m218036/data/echam6_data/input/input_echam6/volcano_tables/ici5d-ad800-1999.asc_modified aodreff_crow.dat
#${SCR_DIR}/link.sh ${INI_ECHAM}/ici5d-ad800-1999.asc_modified aodreff_crow.dat
#--------------------------
# Variable solar irradiance
#--------------------------
T_DIR=${INI_ECHAM}/solar_irradiance/
#
${SCR_DIR}/link.sh ${T_DIR}/swflux_14band_${YEARM1}.nc swflux_${YEARM1}.nc
${SCR_DIR}/link.sh ${T_DIR}/swflux_14band_${YEAR}.nc   swflux_${YEAR}.nc
${SCR_DIR}/link.sh ${T_DIR}/swflux_14band_${YEARP1}.nc swflux_${YEARP1}.nc
${SCR_DIR}/link.sh ${T_DIR}/swflux_14band_${YEARP2}.nc swflux_${YEARP2}.nc
#-----------------
# Greenhouse gases
#-----------------
${SCR_DIR}/link.sh ${INI_ECHAM}/greenhouse_rcp45.nc    greenhouse_gases.nc
#---------
# 3d Ozone
#---------
T_DIR=${INI_ECHAM}/T${HRES}/ozone
#
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_ozone_CMIP5_${YEARM1}.nc ozon${YEARM1}
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_ozone_CMIP5_${YEAR}.nc   ozon${YEAR}
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_ozone_CMIP5_${YEARP1}.nc ozon${YEARP1}
${SCR_DIR}/link.sh ${T_DIR}/T${HRES}_ozone_CMIP5_${YEARP2}.nc ozon${YEARP2}
#--------------------------
# Satellite diagnostic
#--------------------------
${SCR_DIR}/link.sh ${INI_ECHAM}/CFMIP/pointlocations.txt pointlocations.txt
