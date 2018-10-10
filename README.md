# MODIS_processor
This repository contains Matlab codes I made to process different MODIS products. Currently it contains the following:
 1)Caculate the mean emissivity using MOD11A1 and MYD11A1 band 31 & 32 (so, four) images for a study area of a day (Emis_process.m and MODISimg.m);
 
 2)Caculate the fraction of a land cover type (from Land Cover Type 3 of the MCD12Q1 product) over regular grid cells with a user-specific resolution (LC_process.m, LC_Frac.m and MODISimg.m); and
 
 3)Read, mosaic, project, crop and resample MOD13Q1 and MYD13Q1 products for a study area of a day (VI_process.m and MODISimg.m);

Additional Notes:
 - For 2), if the grid cell resolution is finer than MCD12Q1, nearest interpolation is used;
 - The MODISimg.m function is evoked in 1), 2) and 3) to read, mosaic, project, crop and resample multiple MODIS tiles for a study area of a day using GDAL;
 - The MODISimg.m function can also work alone for the users' specific problems on reading, mosaicing, projecting, cropping and resampling multiple tiles of other products.
