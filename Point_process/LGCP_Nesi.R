# r script for testing LGCP models
# creator: Alexandra Strang
# created: 2025

# set working directory
setwd("/home/stranga/00_nesi_projects/landcare04225/Alexandra_Data/point_process_data/LGCP_nesi_data")

# load packages 
library(sf)
library(fmesher)
library(ggplot2)
library(INLA) # version 25.09.04
library(inlabru) # for lgcp() version 2.13.0
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

# add 100 m buffer around coastline boundary
buff_boundary <- st_buffer(sf_Crozier_boundary, dist = 100)

# mesh parameters
# this is 1/15 study size using x range
# x range is longer than y here
Crozier_max.edge <- diff(range(st_coordinates(sf_Crozier)[,1]))/(3*5)
print(Crozier_max.edge)
# ~150 metres

# expand outer layer
Crozier_bound.outer = diff(range(st_coordinates(sf_Crozier)[,1]))/5
print(Crozier_bound.outer)
# ~450 metres

# create Crozier mesh
Crozier_mesh <- fm_mesh_2d(boundary = buff_boundary, # use buffered coastline as boundary
                           max.edge = c(1,5)*Crozier_max.edge, # inner and outer max edge where outer layer has triangle density lower than inner
                           offset = c(Crozier_max.edge, Crozier_bound.outer),
                           cutoff = Crozier_max.edge/5,
                           crs = st_crs(sf_Crozier))
print(Crozier_mesh$n)

# plot boundary mesh
mesh_plot_Crozier <- ggplot() + 
  geom_fm(data = Crozier_mesh) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
  ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh_plot_Crozier

# null LGCP

# import guano area shapefile
sf_Crozier_guano <- st_read("Crozier_2020_1_3031_guano.shp")
# ensure shapefile has right crs code
sf_Crozier_guano <- st_transform(sf_Crozier_guano, crs = st_crs(sf_Crozier))

# create 100m buffer around guano area
sf_Crozier_guano_buffered <- st_buffer(sf_Crozier_guano, dist = 100)
sf_Crozier_guano_buffered <- st_sf(geometry = st_union(sf_Crozier_guano_buffered))
crs(sf_Crozier_guano_buffered)

GA_plot_Crozier <- ggplot() + 
  geom_fm(data = Crozier_mesh) +
  geom_sf(data = sf_Crozier_guano_buffered, fill = NA, color = "red") +
  geom_sf(data = sf_Crozier, color = "purple", size = 1.7, alpha = 0.5) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
  ) + 
  theme_minimal()
# geom_sf converts to degrees

GA_plot_Crozier

# define the SPDE priors (matern)
matern <- inla.spde2.pcmatern(mesh = Crozier_mesh,
                              prior.range = c(1000, 0.5), # distance decay in metres
                              prior.sigma = c(1, 0.5)) # amount of spatial variation

# formula specification of model components
# specify a model where for 2D models geometry is on the left of ~
# and an SPDE + Intercept(1) on the right
# define the domain of the LGCP and model components
# (spatial SPDE effect and Intercept)
null_cmp <- geometry ~
  Intercept(1) + # fixed effect (eventually add here the covariates of slope etc. using rasters)
  mySmooth(geometry, model = matern) # random effect

# formula (geometry is response = point locations)
# mysmooth for spatial autocorrelation

print("running null model with 1000, 0.5 range prior and 1, 0.5 sigma prior")

# use lgcp() with 2D model components, the sf points and the sf boundary
null_model <- lgcp(null_cmp, # formula
                   data = sf_Crozier, # locations
                   samplers = sf_Crozier_guano_buffered, # sample area
                   domain = list(geometry = Crozier_mesh), # mesh
                   options = list(
                     control.inla = list(verbose = TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
                   )
)
# takes less than ~2 mins to run

summary(null_model)

# predicting intensity 
# predict the spatial intensity surface
# need to sum over prediction grid
# make prediction grid as sf points
grid_pts <- fm_pixels(
  Crozier_mesh,
  format = "sf",
  mask = sf_Crozier_guano_buffered
)

# predict intensity per m2
null_lambda <- predict(
  null_model,
  grid_pts, # use this instead of fm_pixels
  ~ exp(mySmooth + Intercept)
)

# cell spacing from the point coordinates
coords <- st_coordinates(grid_pts)
dx <- median(diff(sort(unique(coords[,1]))))
dy <- median(diff(sort(unique(coords[,2]))))

cell_area <- dx * dy  # m2 per pixel

# multiply and sum
null_total_pred <- sum(null_lambda$mean * cell_area)
print(null_total_pred)

# plot log intensity
null_Intensity_plot <- ggplot() +
  geom_fm(data = Crozier_mesh) +
  gg(null_lambda, geom = "tile")

null_Intensity_plot


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
filtered_vals <- vals[!is.na(vals)]
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

hist(percent_guano_raster)
summary(percent_guano_raster)

# check for misalignment
covariate_plot <- c(percent_guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)
plot(covariate_plot)

# check for collinearity between covariates
cov_stack <- c(percent_guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)

cov_values <- as.data.frame(cov_stack, na.rm = TRUE)

cor(cov_values)

# subdivide mesh
# splits triangles into subtriangles
mesh_sub <- fm_subdivide(Crozier_mesh,3)
print(mesh_sub$n)

# need to adjust prior range as spatial autocorrelation wont be explaining as much

# update model formula to include covariates
matern <- inla.spde2.pcmatern(mesh = mesh_sub,
                              prior.range = c(25, 0.5), # distance decay in metres
                              prior.sigma = c(0.01, 0.01)) # amount of spatial variation

print("running full model with 25, 0.5 range prior and 0.01, 0.01 sigma prior and mesh sub 3")

# full
Full_cmp <- geometry ~
  Intercept(1) + 
  percentguano(percent_guano_raster, model = "linear") +
  slope(slope_raster, model = "linear") +
  mySmooth(geometry, model = matern) # random effect

Full_model <- lgcp(Full_cmp, # formula
                   data = sf_Crozier, # locations
                   samplers = sf_Crozier_guano_buffered, # sample area
                   domain = list(geometry = mesh_sub), # mesh
                   options = list(
                     control.inla = list(verbose = TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
)

summary(Full_model)

# need to sum over prediction grid
# make prediction grid as sf points
grid_pts <- fm_pixels(
  mesh_sub,
  format = "sf",
  mask = sf_Crozier_guano_buffered
)

# predict intensity per m2
full_lambda <- predict(
  Full_model,
  grid_pts, # use this instead of fm_pixels
  ~ exp(mySmooth + Intercept + percentguano + slope)
)

# cell spacing from the point coordinates
coords <- st_coordinates(grid_pts)
dx <- median(diff(sort(unique(coords[,1]))))
dy <- median(diff(sort(unique(coords[,2]))))
print(dx)
print(dy)

cell_area <- dx * dy  # m2 per pixel
print(cell_area)

# multiply and sum
full_total_pred <- sum(full_lambda$mean * cell_area)
print(full_total_pred)

# plot log intensity
full_Intensity_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(full_lambda, geom = "tile")

full_Intensity_plot

print("running guano model with 25, 0.5 range prior and 0.01, 0.01 sigma prior and mesh sub 3")

guano_cmp <- geometry ~
  Intercept(1) + 
  percentguano(percent_guano_raster, model = "linear") +
  mySmooth(geometry, model = matern) # random effect

guano_model <- lgcp(guano_cmp, # formula
                    data = sf_Crozier, # locations
                    samplers = sf_Crozier_guano_buffered, # sample area
                    domain = list(geometry = mesh_sub), # mesh
                    options = list(
                      control.inla = list(verbose = TRUE),
                      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
)

summary(guano_model)

# need to sum over prediction grid
# make prediction grid as sf points
grid_pts <- fm_pixels(
  mesh_sub,
  format = "sf",
  mask = sf_Crozier_guano_buffered
)

# predict intensity per m2
guano_lambda <- predict(
  guano_model,
  grid_pts, # use this instead of fm_pixels
  ~ exp(mySmooth + Intercept + percentguano)
)

# cell spacing from the point coordinates
coords <- st_coordinates(grid_pts)
dx <- median(diff(sort(unique(coords[,1]))))
dy <- median(diff(sort(unique(coords[,2]))))
print(dx)
print(dy)

cell_area <- dx * dy  # m2 per pixel
print(cell_area)

# multiply and sum
guano_total_pred <- sum(guano_lambda$mean * cell_area)
print(guano_total_pred)

# plot log intensity
guano_Intensity_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(guano_lambda, geom = "tile")

guano_Intensity_plot

# save model outputs
saveRDS(null_model, file = "null_model.rds")
saveRDS(Full_model, file = "full_model.rds")
saveRDS(guano_model, file = "guano_model.rds")

