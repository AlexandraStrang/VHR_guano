# r script for testing LGCP null models
# creator: Alexandra Strang
# created: 2025

# set working directory
setwd("D:/Points_test")

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

# import guano area shapefile (not needed until modelling spatial covariates)
sf_Crozier_guano <- st_read("Point_process_GA_boundaries/Crozier_2020_1_3031_guano.shp")
# ensure shapefile has right crs code
sf_Crozier_guano <- st_transform(sf_Crozier_guano, crs = st_crs(sf_Crozier))

# import Crozier coastline boundary
sf_Crozier_boundary <- st_read("Crozier_boundary.shp")
# transform shapefile so it has the right crs code
sf_Crozier_boundary <- st_transform(sf_Crozier_boundary, crs = st_crs(sf_Crozier))
st_crs(sf_Crozier_boundary)

# import Crozier UAV bounds shapefile
sf_Crozier_UAV_area <- st_read("Crozier_20201129_UAV/Crozier_20201129_3031_UAV_area.shp")
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

# plot boundary mesh
mesh_plot_Crozier <- ggplot() + 
  geom_fm(data = Crozier_mesh) + 
  geom_sf(data = sf_Crozier_boundary, fill = NA, color = "blue", linetype = "dashed") + 
  labs( 
    x = "Easting", 
    y = "Northing", 
    ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh_plot_Crozier

# plot boundary and samplers area
mesh_plot_Crozier2 <- ggplot() + 
  geom_fm(data = Crozier_mesh) + 
  geom_sf(data = sf_Crozier_boundary, fill = NA, color = "blue", linetype = "dashed") +
  geom_sf(data = sf_Crozier_UAV_area, fill = NA, color = "orange") + 
  #geom_sf(data = sf_Crozier_guano, fill = NA, color = "red", linetype = "dashed") + 
  #geom_sf(data = sf_Crozier, color = "purple", size = 1.7, alpha = 0.5) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
    ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh_plot_Crozier2

# plot boundary, GA, and points
mesh_plot_Crozier3 <- ggplot() + 
  geom_fm(data = Crozier_mesh) + 
  geom_sf(data = sf_Crozier_boundary, fill = NA, color = "blue", linetype = "dashed") + 
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "red", linetype = "dashed") + 
  geom_sf(data = sf_Crozier, color = "purple", size = 1.7, alpha = 0.5) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
    ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh_plot_Crozier3

# plot together
Together <- plot(ggarrange(mesh_plot_Crozier, 
                           mesh_plot_Crozier2,
                           mesh_plot_Crozier3,
                        ncol = 3, nrow = 1, labels=c("a","b","c")))
#annotate_figure(Together, left = "Northing", bottom = "Easting")


# use Cape Crozier 2019 GA and points
# for model training as it has better precison and
# foercasts rather than hindcasts

# create the spatial random field
# fit an LGCP model to the locations of the penguin nests

# Penalized Complexity (PC) priors
# the probability of the range exceeding that to be 0.5
# the prior for the variance explained by the spatial effect
# is set that the probability that the SD exceeds 1 is 0.5

# define the SPDE priors (matern)
matern <- inla.spde2.pcmatern(mesh = Crozier_mesh,
                              prior.range = c(250, 0.5), # distance decay in metres
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

# use lgcp() with 2D model components, the sf points and the sf boundary
null_model <- lgcp(null_cmp, # formula
                    data = sf_Crozier, # locations
                    samplers = sf_Crozier_UAV_area, # sample area of UAV bounds
                    domain = list(geometry = Crozier_mesh), # mesh
                    options = list(
                      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
                    )
)
# takes less than ~2 mins to run

summary(null_model)
# mean should be density per m 2 (once unlogged)

# predicting intensity 
# predict the spatial intensity surface
lambda <- predict(
  null_model, 
  fm_pixels(Crozier_mesh, format = "sf", mask = sf_Crozier_UAV_area), # use samplers 
  ~ exp(mySmooth + Intercept) # exp to unlog, mySmooth is field
  )
# sum of mean just gives density

# need to sum over prediction grid
# make prediction grid as sf points
grid_pts <- fm_pixels(
  Crozier_mesh,
  format = "sf",
  mask = sf_Crozier_UAV_area
)

# predict intensity per m2
lambda <- predict(
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
total_pred <- sum(lambda$mean * cell_area)
total_pred

# plot log intensity
Intensity_plot <- ggplot() +
  geom_fm(data = Crozier_mesh) +
  gg(lambda, geom = "tile")

Intensity_plot
# intensity is highest at 0.5 ish (seems right for density, maybe a little low)

# Next steps:
# sensitivity to priors
# covariates
# residual map

# analyse sensitivity of model to prior specifications

# start with the prior range first

# prior ranges to test (100, 250, 500, 750, 1000 m)
range_grid <- tibble(
  r0 = c(100, 250, 500, 750, 1000),  # prior ranges
  pr = 0.5,                       # probability that range < r0
  s0 = 1,                         # baseline sigma anchor
  ps = 0.5                        # prob(sigma > s0)
)

# function for fitting LGCP models for range priors
fit_LGCP <- function(r0, pr, s0, ps) {
  matern <- inla.spde2.pcmatern(
    mesh = Crozier_mesh,
    prior.range = c(r0, pr),
    prior.sigma = c(s0, ps)
  )
  
  cmp <- geometry ~ Intercept(1) + mySmooth(geometry, model = matern)
  
  model <- lgcp(
    cmp,
    data = sf_Crozier,
    samplers = sf_Crozier_UAV_area,
    domain = list(geometry = Crozier_mesh),
    options = list(control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
  )
  
  lambda <- predict(model, grid_pts, ~ exp(mySmooth + Intercept))
  total_pred <- sum(lambda$mean * cell_area)
  
  # hyperparameter summaries
  hp <- model$summary.hyperpar
  get_row <- function(pattern) {
    hp[grepl(pattern, rownames(hp), ignore.case = TRUE), , drop = FALSE]
  }
  range_row <- get_row("range")
  stdev_row <- get_row("stdev|std")
  
  tibble(
    r0 = r0, pr = pr, s0 = s0, ps = ps,
    WAIC = model$waic$waic,
    DIC = model$dic$dic,
    post_range_mean = range_row$mean[1],
    post_range_2.5 = range_row$`0.025quant`[1],
    post_range_97.5 = range_row$`0.975quant`[1],
    post_sd_mean = stdev_row$mean[1],
    post_sd_2.5 = stdev_row$`0.025quant`[1],
    post_sd_97.5 = stdev_row$`0.975quant`[1],
    total_pred = total_pred
  )
}

# Combine results
results_range <- pmap_dfr(range_grid, fit_LGCP)
print(results_range)

# plots

# Posterior Range vs prior r0
ggplot(results_range, aes(x = r0, y = post_range_mean)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = post_range_2.5, ymax = post_range_97.5), width = 30) +
  labs(x = "Prior range anchor r0 (m)",
       y = "Posterior mean range (95% CI)",
       title = "Posterior spatial range vs prior range") +
  theme_minimal()

# Posterior SD vs prior r0
ggplot(results_range, aes(x = r0, y = post_sd_mean)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = post_sd_2.5, ymax = post_sd_97.5), width = 30) +
  labs(x = "Prior range anchor r0 (m)",
       y = "Posterior SD (95% CI)",
       title = "Posterior SD vs prior range") +
  theme_minimal()

# WAIC vs prior r0
ggplot(results_range, aes(x = r0, y = WAIC)) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = "Prior range anchor r0 (m)",
       y = "WAIC",
       title = "Model WAIC vs prior range") +
  theme_minimal()

# DIC vs prior r0
ggplot(results_range, aes(x = r0, y = DIC)) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = "Prior range anchor r0 (m)",
       y = "DIC",
       title = "Model DIC vs prior range") +
  theme_minimal()

# Total predicted count vs prior r0
ggplot(results_range, aes(x = r0, y = total_pred)) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = "Prior range anchor r0 (m)",
       y = "Total predicted count",
       title = "Predicted total vs prior range") +
  theme_minimal()


# Build LGCP model with covarites

# Incorporate guano area and terrain variables as rasters
guano_raster <- rast("Rasters/cleaned_Crozier_2020_1_3031_guano.tif") # could use 2m guano raster
slope_raster <- rast("Rasters/Cape_Crozier_slope.tif") # 2m slope raster
aspect_raster <- rast("Rasters/Cape_Crozier_aspect.tif") # 2m aspect raster
roughness_raster <- rast("Rasters/Cape_Crozier_roughness.tif") # 2m
TRI_raster <- rast("Rasters/Cape_Crozier_TRI.tif") # 2m

# use plot() to view
# use summary() for stats
# use res() for resolution
# use crs() for crs

# zonal stats within UAV bounds
zones <- vect(sf_Crozier_boundary)
zonal(guano_raster, zones, fun = "mean")  # Mean raster value per zone
zonal(slope_raster, zones, fun = "mean")

# change guano crs to match 
crs(guano_raster) <- "EPSG:3031"
crs(guano_raster)

# remove old gunao from guano raster
rcl <- matrix(c(
  0, 0, 0,  # background stays 0
  1, 1, 1,  # fresh guano stays 1
  2, 2, 0   # old guano becomes 0
), ncol = 3, byrow = TRUE)
guano_raster <- classify(guano_raster, rcl, include.lowest = TRUE, right = NA)
plot(guano_raster)
# now should just be fresh guano 1 and everything else 0

# match res/ extent to GA
ref <- guano_raster
slope_raster     <- resample(slope_raster, ref)
aspect_raster    <- resample(aspect_raster, ref)
roughness_raster <- resample(roughness_raster, ref)
TRI_raster       <- resample(TRI_raster, ref)
# other rasters are res of 0.4 now

# scale
# (mean of 0 sd of 1)
standardize <- function(r) {
  m <- global(r, fun = "mean", na.rm = TRUE)[1, 1]
  s <- global(r, fun = "sd", na.rm = TRUE)[1, 1]
  (r - m) / s
}

# apply to continuous variables only
slope_raster     <- standardize(slope_raster)
aspect_raster    <- standardize(aspect_raster)
roughness_raster <- standardize(roughness_raster)
TRI_raster       <- standardize(TRI_raster)

# check for collinearity between covariates?
cov_stack <- c(guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)

cov_values <- as.data.frame(cov_stack, na.rm = TRUE)

cor(cov_values)

############################################################
# testing

# look into NA error
# Count NAs
n_na <- sum(is.na(values(guano_raster)))
cat("Number of NA cells:", n_na)

# Plot to see where they are
plot(is.na(guano_raster), main = "NA Locations in Guano Raster")

# check other covariates too
n_na <- sum(is.na(values(slope_raster)))
cat("Number of NA cells:", n_na)

# Plot to see where they are
plot(is.na(guano_raster), main = "NA Locations in Slope Raster")

# mask covariates to mesh size as a polygon

mesh_coords <- Crozier_mesh$loc[, 1:2]
mesh_points <- st_as_sf(data.frame(x = mesh_coords[,1], y = mesh_coords[,2]),
                        coords = c("x", "y"), crs = 3031)
mesh_hull <- st_convex_hull(st_union(mesh_points))
mesh_vect <- vect(mesh_hull)

guano_raster     <- mask(guano_raster, mesh_vect)
slope_raster     <- mask(slope_raster, mesh_vect)
aspect_raster    <- mask(aspect_raster, mesh_vect)
roughness_raster <- mask(roughness_raster, mesh_vect)
TRI_raster       <- mask(TRI_raster, mesh_vect)

valid_mask <- !is.na(guano_raster) & 
  !is.na(slope_raster) & 
  !is.na(aspect_raster) & 
  !is.na(roughness_raster) & 
  !is.na(TRI_raster)
guano_raster     <- mask(guano_raster, valid_mask, maskvalues = FALSE)
slope_raster     <- mask(slope_raster, valid_mask, maskvalues = FALSE)
aspect_raster    <- mask(aspect_raster, valid_mask, maskvalues = FALSE)
roughness_raster <- mask(roughness_raster, valid_mask, maskvalues = FALSE)
TRI_raster       <- mask(TRI_raster, valid_mask, maskvalues = FALSE)

# neither the polygon or the mask work

# see if any extra NAs
sum(is.na(values(slope_raster)))
sum(is.na(values(guano_raster)))
# check for misalignment
covariate_stack <- c(guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)
plot(covariate_stack)

# should aspect be northness instead?
northness_raster <- cos(aspect_raster * pi / 180)
northness_raster <- standardize(northness_raster)
# correct aspect
corr_aspect_raster <- rast("Rasters/Cape_Crozier_aspect_corrected.tif")
plot(corr_aspect_raster)

# crop rasters/ fill in NAs
slope_raster[is.na(slope_raster)] <- 0

guano_raster[is.na(guano_raster)] <- 0

guano_raster   <- standardize(guano_raster)
plot(guano_raster)

#####################################################
percent_guano_raster <- rast("CrozierGuano_2m.tif")
plot(percent_guano_raster) # doesn't include below 0s

res(percent_guano_raster) # 2m x 2m
unique(percent_guano_raster)

percent_guano_raster[values(percent_guano_raster) > 1] <- 0

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

# tesing done
############################################################

# subdivide mesh
# splits triangles into subtriangles
mesh_sub <- fm_subdivide(Crozier_mesh,3)

# update model formula to include covariates
matern <- inla.spde2.pcmatern(mesh = mesh_sub,
                              prior.range = c(250, 0.5), # distance decay in metres
                              prior.sigma = c(1, 0.5)) # amount of spatial variation

# trial
Full_cmp <- geometry ~
  Intercept(1) + 
  percentguano(percent_guano_raster, model = "linear") +
  #slope(slope_raster, model = "linear") +
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

# need to sum over prediction grid
# make prediction grid as sf points
grid_pts <- fm_pixels(
  Crozier_mesh,
  format = "sf",
  mask = sf_Crozier_UAV_area
)

# predict intensity per m2
lambda <- predict(
  Full_model,
  grid_pts, # use this instead of fm_pixels
  ~ exp(mySmooth + Intercept + percentguano)
)

# cell spacing from the point coordinates
coords <- st_coordinates(grid_pts)
dx <- median(diff(sort(unique(coords[,1]))))
dy <- median(diff(sort(unique(coords[,2]))))

cell_area <- dx * dy  # m2 per pixel

# multiply and sum
total_pred <- sum(lambda$mean * cell_area)
total_pred
