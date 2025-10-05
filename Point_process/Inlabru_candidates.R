# r script for testing inlabru candidate models
# creator: Alexandra Strang
# created: 2025

# set working directory
setwd("C:/Users/astra/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data")

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

# load in mesh parameters

# import Crozier coastline boundary
sf_Crozier_boundary <- st_read("Crozier_boundary.shp")
# transform shapefile so it has the right crs code
sf_Crozier_boundary <- st_transform(sf_Crozier_boundary, crs = st_crs(sf_Crozier))
st_crs(sf_Crozier_boundary)

# add 100 m buffer around coastline boundary
buff_boundary <- st_buffer(sf_Crozier_boundary, dist = 100)

# mesh parameters
Crozier_max.edge <- 90
print(Crozier_max.edge)
# 90 metres

# expand outer layer
Crozier_bound.outer <- 100
print(Crozier_bound.outer)
# 100 metres

# create Crozier mesh
Crozier_mesh <- fm_mesh_2d(boundary = buff_boundary,
                           max.edge = c(1,3)*Crozier_max.edge,
                           offset = c(Crozier_max.edge, Crozier_bound.outer),
                           cutoff = 0.3,
                           crs = st_crs(sf_Crozier))
print(Crozier_mesh$n)

# subdivide mesh
# splits triangles into subtriangles
mesh_sub <- fm_subdivide(Crozier_mesh,3)
print(mesh_sub$n)

# load in covariates

# add covariates at 2m res (aggregated below)

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

summary(percent_guano_raster)

# have to run this again
percent_guano_raster[is.na(percent_guano_raster)] <- 0

# check for misalignment
covariate_plot <- c(percent_guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)

# stack covariates
cov_stack <- c(percent_guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)

# prepare response variable
count_raster <- 
  terra::rasterize(vect(sf_Crozier), cov_stack, fun = sum, background = 0) %>%
  terra::aggregate(fact = 5, fun = sum) %>%
  mask(vect(sf::st_geometry(buff_boundary)))

# counts of nests
plot(count_raster)

count_raster <- count_raster %>%
  cellSize(unit = "m") %>%
  c(count_raster)

res(count_raster) # 10 x 10 m
summary(count_raster)
# 149 penguins in 100 m2

# extract the coordinates for these pixels
counts_df <- crds(count_raster, df = TRUE, na.rm = TRUE) %>%
  bind_cols(values(count_raster, mat = TRUE, na.rm = TRUE)) %>%
  rename(count = sum) %>%
  mutate(present = (count > 0) *1L) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(sf_Crozier))

# aggregate covariates to match count raster
cov_stack <- terra::aggregate(cov_stack, fact = 5, fun = mean)

# add matern function 

# SPDE priors
matern <- inla.spde2.pcmatern(mesh = mesh_sub,
                              prior.range = c(100, 0.9), 
                              prior.sigma = c(0.01, 0.1))

# Candidate model codes (8):
# G - percent guano only
# GS - percent guano + slope
# GSA - percent guano + slope + aspect
# GR - percent guano + roughness
# GRA - percent guano + roughness + aspect
# GT - percent guano + TRI
# GTA - percent guano + TRI + aspect
# GA - percent guano + aspect
# N - Null model (spatial field only)

# Null model
N_cmp
N_model

# Guano model
G_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  field(geometry, model = matern)

G_model <- bru(
  G_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + 
      field,
    E = area
  )
)

summary(G_model)

# Guano + slope model
GS_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  slope(cov_stack$Cape_Crozier_slope, model = "linear") +
  field(geometry, model = matern)

GS_model <- bru(
  GS_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + slope +
      field,
    E = area
  )
)

summary(GS_model)

# Guano + slope + aspect model
GSA_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  slope(cov_stack$Cape_Crozier_slope, model = "linear") +
  aspect(cov_stack$Cape_Crozier_aspect_corrected, model = "linear") +
  field(geometry, model = matern)

GSA_model <- bru(
  GSA_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + slope + aspect +
      field,
    E = area
  )
)

summary(GSA_model)

# Guano + roughness model
GR_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  roughness(cov_stack$Cape_Crozier_Roughness, model = "linear") +
  field(geometry, model = matern)

GR_model <- bru(
  GR_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + roughness +
      field,
    E = area
  )
)

summary(GR_model)

# Guano + roughness + aspect model
GRA_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  roughness(cov_stack$Cape_Crozier_Roughness, model = "linear") +
  aspect(cov_stack$Cape_Crozier_aspect_corrected, model = "linear") +
  field(geometry, model = matern)

GRA_model <- bru(
  GRA_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + roughness + aspect +
      field,
    E = area
  )
)

summary(GRA_model)

# Guano + TRI model
GT_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  TRI(cov_stack$Cape_Crozier_TRI, model = "linear") +
  field(geometry, model = matern)

GT_model <- bru(
  GT_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + TRI +
      field,
    E = area
  )
)

summary(GT_model)

# Guano + TRI + aspect model
GTA_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  TRI(cov_stack$Cape_Crozier_TRI, model = "linear") +
  aspect(cov_stack$Cape_Crozier_aspect_corrected, model = "linear") +
  field(geometry, model = matern)

GTA_model <- bru(
  GTA_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + TRI + aspect +
      field,
    E = area
  )
)

summary(GTA_model)

# Guano + aspect model
GA_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  aspect(cov_stack$Cape_Crozier_aspect_corrected, model = "linear") +
  field(geometry, model = matern)

GA_model <- bru(
  GA_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + aspect +
      field,
    E = area
  )
)
# aspect has some NAs

summary(GA_model)

