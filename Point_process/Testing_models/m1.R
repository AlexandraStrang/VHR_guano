# m1: 250,0.9 and 1,0.5

# r script for testing LGCP models
# creator: Alexandra Strang
# created: 2025

# set working directory
# path to the data
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

# results nobackup path
results_path <- "/home/stranga/00_nesi_projects/landcare04225_nobackup/Alexandra_results"

# model path
m1_path <- file.path(results_path, "m1")

# ensure model folder exists
if (!dir.exists(m1_path)) dir.create(m1_path, recursive = TRUE)

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

# import guano area shapefile
sf_Crozier_guano <- st_read("Crozier_2020_1_3031_guano.shp")
# ensure shapefile has right crs code
sf_Crozier_guano <- st_transform(sf_Crozier_guano, crs = st_crs(sf_Crozier))

# create 100m buffer around guano area
sf_Crozier_guano_buffered <- st_buffer(sf_Crozier_guano, dist = 100)
sf_Crozier_guano_buffered <- st_sf(geometry = st_union(sf_Crozier_guano_buffered))
crs(sf_Crozier_guano_buffered)

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

# covariates

# need to adjust prior range as spatial autocorrelation wont be explaining as much

# SPDE priors
matern <- inla.spde2.pcmatern(mesh = mesh_sub,
                              prior.range = c(250, 0.9), 
                              prior.sigma = c(1, 0.5))

print("running guano model with 250, 0.9 range prior and 1, 0.5 sigma prior and mesh sub 3")

G_cmp <- geometry ~
  Intercept(1) + 
  percentguano(percent_guano_raster, model = "linear") +
  mySmooth(geometry, model = matern) # random effect

G_model <- lgcp(G_cmp, # formula
                data = sf_Crozier, # locations
                samplers = sf_Crozier_guano_buffered, # sample area
                domain = list(geometry = mesh_sub), # mesh
                options = list(
                  control.inla = list(verbose = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
)

summary(G_model)

# need to sum over prediction grid
# make prediction grid as sf points
grid_pts <- fm_pixels(
  mesh_sub,
  format = "sf",
  mask = sf_Crozier_guano_buffered
)

# predict intensity per m2
G_lambda <- predict(
  G_model,
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
G_total_pred <- sum(G_lambda$mean * cell_area)
print(G_total_pred)

# plot log intensity
G_Intensity_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(G_lambda, geom = "tile")

G_Intensity_plot

ggsave(file.path(m1_path, "G_Intensity_plot.png"), G_Intensity_plot, 
       width = 8, height = 5, units = "in", 
       dpi = 600)

print("running GS model with 250, 0.9 range prior and 1, 0.5 sigma prior and mesh sub 3")

# Guano and slope
GS_cmp <- geometry ~
  Intercept(1) + 
  percentguano(percent_guano_raster, model = "linear") +
  slope(slope_raster, model = "linear") +
  mySmooth(geometry, model = matern) # random effect

GS_model <- lgcp(GS_cmp, # formula
                 data = sf_Crozier, # locations
                 samplers = sf_Crozier_guano_buffered, # sample area
                 domain = list(geometry = mesh_sub), # mesh
                 options = list(
                   control.inla = list(verbose = TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
)

summary(GS_model)

# need to sum over prediction grid
# make prediction grid as sf points
grid_pts <- fm_pixels(
  mesh_sub,
  format = "sf",
  mask = sf_Crozier_guano_buffered
)

# predict intensity per m2
GS_lambda <- predict(
  GS_model,
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
GS_total_pred <- sum(GS_lambda$mean * cell_area)
print(GS_total_pred)

# plot log intensity
GS_Intensity_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GS_lambda, geom = "tile")

GS_Intensity_plot

ggsave(file.path(m1_path, "GS_Intensity_plot.png"), GS_Intensity_plot, 
       width = 8, height = 5, units = "in", 
       dpi = 600)

# save model outputs
saveRDS(null_model, file = file.path(m1_path, "null_model.rds"))
saveRDS(GS_model, file = file.path(m1_path, "GS_model.rds"))
saveRDS(G_model, file = file.path(m1_path, "G_model.rds"))
