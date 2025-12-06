# r script for correcting aspect values for polar regions and creating northness and eastness rasters
# modified by: Alexandra Strang
# date: 2025

# Aspect correction code provided by Annie Schmidt 
# from: https://github.com/pointblue/adpe_subcol_success/blob/main/code/data_prep/aspect_correction.R

# set working directory
setwd("C:/Users/astra/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data")

# load packages
library(dplyr)
library(circular)
library(terra)
library(raster)
library(maptools)

#######################################################################################################################################################
# aspect correction function modified from matlab function by 
# Antarctic Science 19 (1), 129-130 (2007) ? Antarctic Science Ltd Printed in the UK DOI: 10.1017/S0954102007000181
# 129
# Short Note
# Correcting GIS-based slope aspect calculations for the Polar Regions
# GEOFF J.M. MORET and AUDREY D. HUERTA
#######################################################################################################################################################
trueaspect= function(data,dx,xll,yll){
  #TRUEASPECT corrects slope aspect data so that angles are relative to
  #geographic north instead of grid north. 
  # The inputs are: 
  # data: a matrix of slope aspect data, with NoData cells set to -9999 and cells without an aspect set to -1 
  # dx: the cell size 
  # xll: The x-coordinate of the lower left hand grid cell. This should be in polar stereographic projection, with false easting removed.
  # yll: The y-coordinate of the lower left hand grid cell. This should be in a polar stereographic projection, with false northing removed.
  # Original function included a variable hemi (Use hemi = 1 for the south pole and hemi = 2 for the north pole.) which I removed 
  # Also removed range (If range = 1, the output values will be between 0 and 360 degrees) because only want in 0-360
  # If range has any other value, they will be between -180 and +180 degrees.
  nrows=nrow(data)
  ncols=ncol(data)
  for(m in 1:nrows){
    for(n in 1:ncols){
      ifelse(data[m,n]==-1,9999,data[m,n])
    }
  }
  x=xll+c(0:(ncols-1))*dx
  y=yll+c(0:(nrows-1))*dx
  newdata=matrix(0,nrows,ncols)
  for (i in 1:nrows){
    for(j in 1:ncols){
      lon=atan2(x[j],y[i])*180/pi
      newdata[i,j]=data[i,j]-lon
    }
  }
  
  newdata[newdata< -1000]<- -9999
  newdata[newdata<0&newdata> -9999]<- newdata[newdata<0&newdata> -9999]+360
  newdata[newdata>1000] <- -1
  newdata[newdata>360] <- newdata[newdata>360]-360
  newdata
}

# correct crozier aspect

# read in aspect raster
croz_aspect_raw <- raster("Crozier_terrain_mesh/Cape_Crozier_aspect.tif") # 2m aspect raster
crs(croz_aspect_raw)
hist(croz_aspect_raw)
# get lower left coordinates from raster
c_xll <- xmin(croz_aspect_raw)
c_yll <- ymin(croz_aspect_raw)
# Convert to matrix
c_data <- raster::as.matrix(croz_aspect_raw, mode="numeric")
# replace NA with -9999
c_data[is.na(c_data)]<- -9999

# calculate corrected aspect
c_aspect_correct <- trueaspect(c_data,dx=2,xll=c_xll,yll=c_yll)
# check that distribution looks right
hist(c_aspect_correct[!c_aspect_correct==-9999])
# replace -9999 with NA (not sure this step is meaningful)
c_aspect_correct[c_aspect_correct==-9999]<- NA

# save as raster
c_aspect_corr_rast <- raster(c_aspect_correct,xmn=croz_aspect_raw@extent@xmin,
                             xmx=croz_aspect_raw@extent@xmax,
                             ymn=croz_aspect_raw@extent@ymin,
                             ymx=croz_aspect_raw@extent@ymax, crs=croz_aspect_raw@crs)
c_aspect_corr_rast@file@nodatavalue<- -9999
c_aspect_corr_rast@data@min<- 0

writeRaster(c_aspect_corr_rast,"Crozier_terrain_mesh/Cape_Crozier_aspect_corrected.tif", overwrite=TRUE)

# correct aspect
corr_aspect_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_aspect_corrected.tif")
plot(corr_aspect_raster)

# convert aspect raster to northness and eastness rasters

# convert aspect from degrees to radians
aspect_radians <- corr_aspect_raster * pi / 180

# calculate northness and eastness
northness <- cos(aspect_radians)
eastness  <- sin(aspect_radians)

# assign names
names(northness) <- "northness"
names(eastness)  <- "eastness"

# change crs to match other covariates
crs(northness) <- "EPSG:3031"
crs(northness)

crs(eastness) <- "EPSG:3031"
crs(eastness)

# inspect results
northness
eastness

plot(northness)
plot(eastness)

# save as raster (run once)
writeRaster(northness, "Crozier_terrain_mesh/Cape_Crozier_northness.tif", overwrite=TRUE)
writeRaster(eastness,  "Crozier_terrain_mesh/Cape_Crozier_eastness.tif",  overwrite=TRUE)
