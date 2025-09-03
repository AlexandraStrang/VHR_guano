# r script for testing LGCP models
# creator: Alexandra Strang
# created: 2025

# set working directory
setwd("/home/stranga/00_nesi_projects/landcare04225/Alexandra_Data/point_process_data/LGCP_nesi_data")

# load packages 
library(sf)
library(fmesher)
library(ggplot2)
library(INLA)
library(ggpubr)
library(inlabru) # for lgcp()
library(dplyr)
library(purrr)
library(tidyr)
library(terra) # for rasters

# sessionInfo() - important to note r version and versions of inlabru, INLA, fmesher
# R version 4.5.0 (2025-04-11 ucrt)

# EPSG 3031 is in metres but when plotting geom_sf will convert to degrees

# read in points from xy csv
# Cape Crozier 2020
Crozier_xy <- read.csv("Crozier_Points_2020_3031.csv")

# convert to sf object and check CRS
sf_Crozier <- st_as_sf(Crozier_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier)

# import Crozier coastline boundary
sf_Crozier_boundary <- st_read("Crozier_boundary.shp")
# transform shapefile so it has the right crs code
sf_Crozier_boundary <- st_transform(sf_Crozier_boundary, crs = st_crs(sf_Crozier))
st_crs(sf_Crozier_boundary)

# import Crozier UAV bounds shapefile
sf_Crozier_UAV_area <- st_read("Crozier_20201129_3031_UAV_area.shp")
# ensure shapefile has right crs code
sf_Crozier_UAV_area <- st_transform(sf_Crozier_UAV_area, crs = st_crs(sf_Crozier))
st_crs(sf_Crozier_UAV_area)

# mesh parameters
# this is 1/15 study size using x range
# x range is longer than y here
Crozier_max.edge <- diff(range(st_coordinates(sf_Crozier)[,1]))/(3*5)
print(Crozier_max.edge)
# ~150 metres

# expand outer layer out by 1/5
Crozier_bound.outer = diff(range(st_coordinates(sf_Crozier)[,1]))/5
print(Crozier_bound.outer)
# ~500 metres

# create Crozier mesh
Crozier_mesh <- fm_mesh_2d(boundary = sf_Crozier_boundary, # use coastline as boundary
                           max.edge = c(1,5)*Crozier_max.edge, # inner and outer max edge where outer layer has triangle density lower than inner
                           offset = c(Crozier_max.edge, Crozier_bound.outer),
                           cutoff = Crozier_max.edge/10,
                           crs = st_crs(sf_Crozier))
# saved mesh as polygon to crop terrain rasters

# Build LGCP model with covarites

# use continuous guano raster
percent_guano_raster <- rast("CrozierGuano_2m.tif")

# change guano crs to match 
crs(percent_guano_raster) <- "EPSG:3031"
crs(percent_guano_raster)

res(percent_guano_raster) # 2m x 2m

# terrain rasters clipped by mesh
slope_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_slope.tif") # 2m slope raster
aspect_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_aspect_corrected.tif") # 2m aspect raster
roughness_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_Roughness.tif") # 2m
TRI_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_TRI.tif") # 2m

# change corrected aspect crs to match 
crs(aspect_raster) <- "EPSG:3031"
crs(aspect_raster)

# match guano area extent to terrain rasters
percent_guano_raster <- extend(percent_guano_raster, slope_raster)

percent_guano_raster[values(percent_guano_raster) > 1] <- 0
percent_guano_raster[is.na(percent_guano_raster)] <- 0

# Extract values from raster
vals <- values(percent_guano_raster)

# Remove NAs and filter out 0s and 1s
filtered_vals <- vals[!is.na(vals) & vals > 0 & vals < 1]
mean(filtered_vals)

# Plot histogram of intermediate values only
hist(filtered_vals,
     main = "Histogram of Guano Percent Values (Excluding 0 & 1)",
     xlab = "Value",
     ylab = "Frequency",
     col = "skyblue",
     breaks = 100)

plot(percent_guano_raster)

# scale
# (mean of 0 sd of 1)
standardize <- function(r) {
  m <- global(r, fun = "mean", na.rm = TRUE)[1, 1]
  s <- global(r, fun = "sd", na.rm = TRUE)[1, 1]
  (r - m) / s
}

# apply to continuous variables
percent_guano_raster   <- standardize(percent_guano_raster)
slope_raster     <- standardize(slope_raster)
aspect_raster    <- standardize(aspect_raster)
roughness_raster <- standardize(roughness_raster)
TRI_raster       <- standardize(TRI_raster)

# check for misalignment
covariate_plot <- c(percent_guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)
plot(covariate_plot)

# check for collinearity between covariates
cov_stack <- c(percent_guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)

cov_values <- as.data.frame(cov_stack, na.rm = TRUE)

cor(cov_values)

# subdivide mesh
# splits triangles into subtriangles
mesh_sub <- fm_subdivide(Crozier_mesh,3) # try 1/9 of max.edge

# update model formula to include covariates
matern <- inla.spde2.pcmatern(mesh = mesh_sub,
                              prior.range = c(150, 0.5), # distance decay in metres
                              prior.sigma = c(5, 0.5)) # amount of spatial variation

# trial with just GA
Full_cmp <- geometry ~
  Intercept(1) + 
  percentguano(percent_guano_raster, model = "linear") +
  slope(slope_raster, model = "linear") +
  mySmooth(geometry, model = matern) # random effect

Full_model <- lgcp(Full_cmp, # formula
                   data = sf_Crozier, # locations
                   samplers = sf_Crozier_UAV_area, # sample area of UAV bounds
                   domain = list(geometry = mesh_sub), # mesh
                   options = list(
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
                   )
)

summary(Full_model)

print(summary(lambda$mean))
hist(lambda$mean, breaks = 100, main = "Predicted Intensity", xlab = "Penguins/m²")

# need to sum over prediction grid
# make prediction grid as sf points
grid_pts <- fm_pixels(
  mesh_sub,
  format = "sf",
  mask = sf_Crozier_UAV_area
)

# predict intensity per m2
lambda <- predict(
  Full_model,
  grid_pts, # use this instead of fm_pixels
  ~ exp(mySmooth + Intercept + percentguano + slope)
)

# cell spacing from the point coordinates
coords <- st_coordinates(grid_pts)
dx <- median(diff(sort(unique(coords[,1]))))
dy <- median(diff(sort(unique(coords[,2]))))
print(paste("Prediction grid resolution:", dx, "x", dy, "meters"))

cell_area <- dx * dy  # m2 per pixel

# multiply and sum
total_pred <- sum(lambda$mean * cell_area)
print(total_pred)

# integration
lambda_integrated <- predict(
  Full_model,
  grid_pts,
  ~ exp(mySmooth + Intercept + percentguano + slope),
  integrate = TRUE
)
print(sum(lambda_integrated$mean))

