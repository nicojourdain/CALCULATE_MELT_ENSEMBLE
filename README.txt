##### GSW TOOLBOX #####
# First you may need the TEOS10 toolbox to convert in-situ to potential temperature (for WOA).
# To avoid issues with updates on the GSW-Fortran tools, I've cloned the 2016 GSW-Fortran
# in this repository. In case you want to check for updates (not recommended), you can still do:
# git clone https://github.com/TEOS-10/GSW-Fortran.git
cd GSW-Fortran/test
# edit makefile (in particular fill the FC and NETCDF_INCDIR variables)
make
./gsw_check  ## to check (some have to work, maybe not all of them)
cd ../..
# The gsw_data_v3_0.nc file provided in this repository should be fine, but in case you start again
# from the original file (i.e. from github.com/TEOS-10/GSW-Fortran.git), you will need to modify it
# to avoid NaN in some places (otherwise, skip the next 4 lines):
./compile.sh remove_NaN_from_gsw_data_v3_0.f90
./remove_NaN_from_gsw_data_v3_0
ncks -x -v ocean_ref,ndepth_ref,deltaSA_ref,SA_ref,SAAR_ref GSW-Fortran/test/gsw_data_v3_0.nc gsw_data_v3_0.nc
ncks -A gsw_data_v3_0_to_be_ncks-A.nc gsw_data_v3_0.nc

###############################################################################################

# Prepare ice shelf mask and attribute an ID and melt rate to each individual ice shelf :
ifort -c $NC_INC identify_ice_shelf_mask.f90 
ifort -o identify_ice_shelf_mask  identify_ice_shelf_mask.o $NC_LIB
./submit.sh identify_ice_shelf_mask 06 64

# Prepare observational T,S forcing for individual ice shelves :

ifort -c $NC_INC -I ./GSW-Fortran/modules calculate_TS_profiles_WOA2018_EN4_MEOP.f90
ifort -o calculate_TS_profiles_WOA2018_EN4_MEOP ./GSW-Fortran/modules/*.o ./GSW-Fortran/toolbox/*.o  calculate_TS_profiles_WOA2018_EN4_MEOP.o $NC_LIB
./submit.sh calculate_TS_profiles_WOA2018_EN4_MEOP 01 16

ifort -c $NC_INC calculate_TS_profiles_SCHMIDTKO2014.f90 
ifort -o calculate_TS_profiles_SCHMIDTKO2014 calculate_TS_profiles_SCHMIDTKO2014.o $NC_LIB
./calculate_TS_profiles_SCHMIDTKO2014

# Prepare CMIP5 T,S anomalies for individual ice shelves :

ifort -c $NC_INC calculate_TS_profiles_CMIP5_anom.f90
ifort -o calculate_TS_profiles_CMIP5_anom calculate_TS_profiles_CMIP5_anom.o $NC_LIB
./submit.sh calculate_TS_profiles_CMIP5_anom 03 32

ifort -c $NC_INC calculate_TS_profiles_CMIP5_anom_bottom.f90
ifort -o calculate_TS_profiles_CMIP5_anom_bottom calculate_TS_profiles_CMIP5_anom_bottom.o $NC_LIB
./submit.sh calculate_TS_profiles_CMIP5_anom_bottom 01 16

##################################################################################################
# Calculate ensemble melt rates :
# Everything based on a single file that can also be run manually: calculate_melt_rate_ensemble.f90

# Calculate box and plume geometries and a few first melt netcdf files for a quick check :
./submit_all_calculate_plume_box_geometries.sh

# once checked, remove all the mlet netcdf files:
rm melt_ensemble_[0-9][0-9][0-9].nc
rm melt_ensemble_[0-9][0-9][0-9]_?.nc
rm melt_ensemble_[0-9][0-9][0-9]_??.nc

# Calculate statistical ensembles :
./submit_all_calculate_melt_rate_ensemble.sh 

# Rebuild Ross and Ronne melt files :
# Adapt ID_isf manually to 10, 11, 21, then :
ifort -c $NC_INC rebuild_melt_large_ice_shelves.f90
ifort -o rebuild_melt_large_ice_shelves rebuild_melt_large_ice_shelves.o $NC_LIB
./submit_short.sh rebuild_melt_large_ice_shelves 20 64

# Rebuild patterns and statistics at the scale of Antarctica :

ifort -c $NC_INC rebuild_melt_ensemble_patterns.f90
ifort -o rebuild_melt_ensemble_patterns rebuild_melt_ensemble_patterns.o $NC_LIB
./submit_short.sh rebuild_melt_ensemble_patterns 15 48

ifort -c $NC_INC rebuild_melt_ensemble_statistics.f90
ifort -o rebuild_melt_ensemble_statistics rebuild_melt_ensemble_statistics.o $NC_LIB
./submit_short.sh rebuild_melt_ensemble_statistics 29 64   ## or execute directly

##################################################################################################
#  NOW WE REDO THE JOB WITH AN ADDITIONAL CORRECTION ON THERMAL FORCING

# Remove all the "DeltaT" melt netcdf files if they already exist :
rm melt_ensemble_DeltaT_[0-9][0-9][0-9].nc
rm melt_ensemble_DeltaT_[0-9][0-9][0-9]_?.nc
rm melt_ensemble_DeltaT_[0-9][0-9][0-9]_??.nc

# Calculate statistical ensembles :
./submit_all_calculate_melt_rate_ensemble_DeltaT.sh

# Rebuild Ross and Ronne melt files :
# Adapt ID_isf manually to 10, 11, 21, then :
ifort -c $NC_INC rebuild_melt_large_ice_shelves_DeltaT.f90
ifort -o rebuild_melt_large_ice_shelves_DeltaT rebuild_melt_large_ice_shelves_DeltaT.o $NC_LIB
./submit_short.sh rebuild_melt_large_ice_shelves_DeltaT 20 64

# Rebuild patterns and statistics at the scale of Antarctica :

ifort -c $NC_INC rebuild_melt_ensemble_patterns_DeltaT.f90
ifort -o rebuild_melt_ensemble_patterns_DeltaT rebuild_melt_ensemble_patterns_DeltaT.o $NC_LIB
./submit_short.sh rebuild_melt_ensemble_patterns_DeltaT 15 48

ifort -c $NC_INC rebuild_melt_ensemble_statistics_DeltaT.f90
ifort -o rebuild_melt_ensemble_statistics_DeltaT rebuild_melt_ensemble_statistics_DeltaT.o $NC_LIB
./submit_short.sh rebuild_melt_ensemble_statistics_DeltaT 29 64

##################################################################################################
# NOW WE INTERPOLATE ONTO A STEREOGRAPHIC GRID :

# rebuild isfmask for larger domain (just to plot the full Peninsula) :
cp -p identify_ice_shelf_mask.f90 tmp.f90
sed -e "s/ANT3/ANT/g ; s/Created using identify_ice_shelf_mask.f90/Created using identify_ice_shelf_mask_EXTENDED.f90/g" tmp.f90 > identify_ice_shelf_mask_EXTENDED.f90
ifort -c $NC_INC identify_ice_shelf_mask_EXTENDED.f90
ifort -o identify_ice_shelf_mask_EXTENDED  identify_ice_shelf_mask_EXTENDED.o $NC_LIB
./submit.sh identify_ice_shelf_mask_EXTENDED 06 64

ifort -c $NC_INC put_isfmask_on_stereographic_grid.f90
ifort -o put_isfmask_on_stereographic_grid  put_isfmask_on_stereographic_grid.o $NC_LIB
./submit.sh put_isfmask_on_stereographic_grid 03 16

ifort -c $NC_INC put_melt_pattern_on_stereographic_grid.f90
ifort -o put_melt_pattern_on_stereographic_grid  put_melt_pattern_on_stereographic_grid.o $NC_LIB
./submit.sh put_melt_pattern_on_stereographic_grid 20 64

# to build the mean melt pattern (averaged across all parameterizations):
ifort -c $NC_INC calculate_mean_melt_all_param_DeltaT.f90
ifort -o calculate_mean_melt_all_param_DeltaT calculate_mean_melt_all_param_DeltaT.o $NC_LIB
./submit_short.sh calculate_mean_melt_all_param_DeltaT 15 48

##################################################################################################
#  NOW WE REDO THE JOB WITH IDEALIZED WARMING (FROM 0.2 to 6.0 degC) :

# Remove all the "DeltaT" melt netcdf files if they already exist :
rm melt_ensemble_DeltaT_IDEALIZED_[0-9][0-9][0-9].nc
rm melt_ensemble_DeltaT_IDEALIZED_[0-9][0-9][0-9]_?.nc
rm melt_ensemble_DeltaT_IDEALIZED_[0-9][0-9][0-9]_??.nc

# Calculate statistical ensembles :
./submit_all_calculate_melt_rate_ensemble_DeltaT_IDEALIZED.sh

# Rebuild Ross and Ronne melt files :
# Adapt ID_isf manually to 10, 11, 21, then :
ifort -c $NC_INC rebuild_melt_large_ice_shelves_DeltaT_IDEALIZED.f90
ifort -o rebuild_melt_large_ice_shelves_DeltaT_IDEALIZED rebuild_melt_large_ice_shelves_DeltaT_IDEALIZED.o $NC_LIB
./submit_short.sh rebuild_melt_large_ice_shelves_DeltaT_IDEALIZED 20 64

# Rebuild patterns and statistics at the scale of Antarctica :

ifort -c $NC_INC rebuild_melt_ensemble_patterns_DeltaT_IDEALIZED.f90
ifort -o rebuild_melt_ensemble_patterns_DeltaT_IDEALIZED rebuild_melt_ensemble_patterns_DeltaT_IDEALIZED.o $NC_LIB
./submit_short.sh rebuild_melt_ensemble_patterns_DeltaT_IDEALIZED 15 48

ifort -c $NC_INC rebuild_melt_ensemble_statistics_DeltaT_IDEALIZED.f90
ifort -o rebuild_melt_ensemble_statistics_DeltaT_IDEALIZED rebuild_melt_ensemble_statistics_DeltaT_IDEALIZED.o $NC_LIB
./submit_short.sh rebuild_melt_ensemble_statistics_DeltaT_IDEALIZED 29 64
