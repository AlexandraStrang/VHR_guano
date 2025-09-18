# r script for testing LGCP model sensitivity 
# creator: Alexandra Strang
# created: 2025

# needs updating

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

#################################################################################################################

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

#################################################################################################################

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

#################################################################################################################

# import guano area shapefile for samplers
sf_Crozier_guano <- st_read("Crozier_2020_1_3031_guano.shp")
# ensure shapefile has right crs code
sf_Crozier_guano <- st_transform(sf_Crozier_guano, crs = st_crs(sf_Crozier))

# create 100m buffer around guano area
sf_Crozier_guano_buffered <- st_buffer(sf_Crozier_guano, dist = 100)
sf_Crozier_guano_buffered <- st_sf(geometry = st_union(sf_Crozier_guano_buffered))
crs(sf_Crozier_guano_buffered)

#################################################################################################################

# add null model

#################################################################################################################

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
