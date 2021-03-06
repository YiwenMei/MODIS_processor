# MODIS_processor
This repository contains Matlab codes I made to process different MODIS products. Currently it contains the following:

 1)Caculate the mean emissivity using MOD11A1 and MYD11A1 band 31 & 32 (so, four) images for a study area of a day (Emis_process.m and MODISimg.m);
 
 2)Read, mosaic, project, crop and resample MCD12Q1 product for a study area of a year and then caculate the fraction of a land cover type (from the processed MCD12Q1) over regular grid cells with a user-specific resolution (LC_process.m, MODISimg.m, and LC_Frac.m); 
 
 3)Read, mosaic, project, crop and resample MOD13Q1 and MYD13Q1 products for a study area of a day (VI_process.m and MODISimg.m);
 
 4)Calculate the mean actural/potential evapotranspiration rate of MOD16A2 and MYD16A2 images for a study area of a time step (ET_process.m, and MODISimg.m); and
 
 5)Calculate the mean NDSI of MOD10A1 and MYD10A1 images for a study area of a time step (SC_process.m, and MODISimg.m);
 
 6)Read, mosaic, project, crop and resample MCD43A3 products for a study area of a day (ALB_process.m and MODISimg.m);
 
 7)Read, mosaic, project, crop and resample SRTM product for a study area (SRTMimg.m);
 
 8)Calculate the mean values of a variable for a time step using its two nearest neighbors in time (Tinterp2D.m);
 
Additional Notes:
 - For 2), if the grid cell resolution is finer than MCD12Q1, nearest interpolation is used;
 - The MODISimg.m function is evoked in all processes to read, mosaic, project, crop and resample multiple MODIS tiles for a study area of a day using GDAL;
 - The MODISimg.m function can also work alone for the users' specific problems on reading, mosaicing, projecting, cropping and resampling multiple tiles of other products.
