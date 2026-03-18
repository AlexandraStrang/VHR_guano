# r script for testing inlabru candidate models
# creator: Alexandra Strang
# created: 2025

# set working directory
setwd("C:/Users/ajs424/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data")

# load packages 
library(sf)
library(fmesher)
library(ggplot2)
library(INLA) # version 25.09.19
library(inlabru) # for bru() version 2.13.0
library(dplyr)
library(terra) # for rasters

bru_options_set(control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))

# sessionInfo() - important to note r version and versions of inlabru, INLA, fmesher
# R version 4.5.1

# EPSG 3031 is in metres but when plotting geom_sf will convert to degrees

##############################################################################################
# Load data
##############################################################################################

# read in points from xy csv
# Cape Crozier 2020
Crozier_xy <- read.csv("Crozier_UAV_points/Reprojected_3031/Crozier_2020_Points_3031.csv")

# convert to sf object and check CRS
sf_Crozier <- st_as_sf(Crozier_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier)

# mesh and coastline boundary
mesh_sub <- readRDS("Inlabru_outputs/mesh_sub.rds")
buff_boundary <- readRDS("Inlabru_outputs/buff_boundary.rds")

##############################################################################################
# covariates
##############################################################################################

# add covariates at 2m res (aggregated below)

# use continuous guano raster
percent_guano_raster <- rast("CrozierGuano_2m.tif")

# change guano crs to match 
crs(percent_guano_raster) <- "EPSG:3031"
crs(percent_guano_raster)

res(percent_guano_raster) # 2m x 2m

# 2m terrain rasters clipped by mesh
slope_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_slope.tif")
northness_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_northness.tif")
eastness_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_eastness.tif")
roughness_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_Roughness.tif")
TRI_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_TRI.tif")

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
northness_raster    <- standardize(northness_raster)
eastness_raster    <- standardize(eastness_raster)
roughness_raster <- standardize(roughness_raster)
TRI_raster       <- standardize(TRI_raster)

summary(percent_guano_raster)

# have to run this again
percent_guano_raster[is.na(percent_guano_raster)] <- 0

# check for misalignment
covariate_plot <- c(percent_guano_raster, slope_raster, northness_raster, eastness_raster, roughness_raster, TRI_raster)
plot(covariate_plot)

# stack covariates
cov_stack <- c(percent_guano_raster, slope_raster, northness_raster, eastness_raster, roughness_raster, TRI_raster)

# check for colinearity between covariates
cov_values <- as.data.frame(cov_stack, na.rm = TRUE)
cor(cov_values)

##############################################################################################
# Count raster
##############################################################################################

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

# sum of area that has zero penguins
cell_areas <- count_raster[[1]]
zero_cells <- count_raster[[2]] == 0
zero_area <- mask(cell_areas, zero_cells, maskvalue = FALSE)
total_zero_area <- global(zero_area, fun = "sum", na.rm = TRUE)
print(total_zero_area)

total_raster_area <- global(cell_areas, fun = "sum", na.rm = TRUE)
print(total_raster_area)
# total - zero area = 585,594 m2 (area with penguins)

# extract the coordinates for these pixels
counts_df <- crds(count_raster, df = TRUE, na.rm = TRUE) %>%
  bind_cols(values(count_raster, mat = TRUE, na.rm = TRUE)) %>%
  rename(count = sum) %>%
  mutate(present = (count > 0) *1L) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(sf_Crozier))

# extract coordinates from sf geometry
counts_df <- counts_df %>%
  mutate(
    x_coord = st_coordinates(.)[, 1],
    y_coord = st_coordinates(.)[, 2]
  )

# save
saveRDS(counts_df, file = "Inlabru_outputs/counts_df.rds")

# aggregate covariates to match count raster
cov_stack <- terra::aggregate(cov_stack, fact = 5, fun = mean)

##############################################################################################
# run models
##############################################################################################

set.seed(28)

# add matern function 

# SPDE priors
matern <- inla.spde2.pcmatern(mesh = mesh_sub,
                              prior.range = c(100, 0.9), 
                              prior.sigma = c(0.01, 0.1))

# Candidate model codes (8):
# G - percent guano only
# GS - percent guano + slope
# GSNE - percent guano + slope + northness + eastness
# GR - percent guano + roughness
# GRNE - percent guano + roughness + northness + eastness
# GT - percent guano + TRI
# GTNE - percent guano + TRI + northness + eastness
# GNE - percent guano + northness + eastness
# N - Null model (spatial field only)

# northness and eastness have some NAs

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

saveRDS(G_model, file = "Inlabru_outputs/G_model.rds")

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

saveRDS(GS_model, file = "Inlabru_outputs/GS_model.rds")

# Guano + slope + northness + eastness model
GSNE_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  slope(cov_stack$Cape_Crozier_slope, model = "linear") +
  northness(cov_stack$northness, model = "linear") +
  eastness(cov_stack$eastness, model = "linear") +
  field(geometry, model = matern)

GSNE_model <- bru(
  GSNE_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + slope + northness + eastness +
      field,
    E = area
  )
)

summary(GSNE_model)

saveRDS(GSNE_model, file = "Inlabru_outputs/GSNE_model.rds")

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

saveRDS(GR_model, file = "Inlabru_outputs/GR_model.rds")

# Guano + roughness + northness + eastness model
GRNE_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  roughness(cov_stack$Cape_Crozier_Roughness, model = "linear") +
  northness(cov_stack$northness, model = "linear") +
  eastness(cov_stack$eastness, model = "linear") +
  field(geometry, model = matern)

GRNE_model <- bru(
  GRNE_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + roughness + northness + eastness +
      field,
    E = area
  )
)

summary(GRNE_model)

saveRDS(GRNE_model, file = "Inlabru_outputs/GRNE_model.rds")

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

saveRDS(GT_model, file = "Inlabru_outputs/GT_model.rds")

# Guano + TRI + northness + eastness model
GTNE_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  TRI(cov_stack$Cape_Crozier_TRI, model = "linear") +
  northness(cov_stack$northness, model = "linear") +
  eastness(cov_stack$eastness, model = "linear") +
  field(geometry, model = matern)

GTNE_model <- bru(
  GTNE_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + TRI + northness + eastness +
      field,
    E = area
  )
)

summary(GTNE_model)

saveRDS(GTNE_model, file = "Inlabru_outputs/GTNE_model.rds")

# Guano + northness + eastness model
GNE_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  northness(cov_stack$northness, model = "linear") +
  eastness(cov_stack$eastness, model = "linear") +
  field(geometry, model = matern)

GNE_model <- bru(
  GNE_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + northness + eastness +
      field,
    E = area
  )
)

summary(GNE_model)

saveRDS(GNE_model, file = "Inlabru_outputs/GNE_model.rds")

# Null model
N_cmp <- ~ Intercept(1) +
  field(geometry, model = matern)

N_model <- bru(
  N_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      field,
    E = area
  )
)

summary(N_model)

saveRDS(N_model, file = "Inlabru_outputs/N_model.rds")

#############################################################################################
# Predictions
##############################################################################################

# G model predictions
G_pred <- predict(
  G_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + 
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

saveRDS(G_pred, file = "Inlabru_outputs/G_pred.rds")

# GS model predictions
GS_pred <- predict(
  GS_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + slope +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

saveRDS(GS_pred, file = "Inlabru_outputs/GS_pred.rds")

# GSNE model predictions
GSNE_pred <- predict(
  GSNE_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + slope + northness + eastness +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

saveRDS(GSNE_pred, file = "Inlabru_outputs/GSNE_pred.rds")

# GR model predictions
GR_pred <- predict(
  GR_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + roughness +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

saveRDS(GR_pred, file = "Inlabru_outputs/GR_pred.rds")

# GRNE model predictions
GRNE_pred <- predict(
  GRNE_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + roughness + northness + eastness +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

saveRDS(GRNE_pred, file = "Inlabru_outputs/GRNE_pred.rds")

# GT model predictions
GT_pred <- predict(
  GT_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + TRI +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

saveRDS(GT_pred, file = "Inlabru_outputs/GT_pred.rds")

# GTNE model predictions
GTNE_pred <- predict(
  GTNE_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + TRI + northness + eastness +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

saveRDS(GTNE_pred, file = "Inlabru_outputs/GTNE_pred.rds")

# GNE model predictions
GNE_pred <- predict(
  GNE_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + northness + eastness +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

saveRDS(GNE_pred, file = "Inlabru_outputs/GNE_pred.rds")

# N model predictions
N_pred <- predict(
  N_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

saveRDS(N_pred, file = "Inlabru_outputs/N_pred.rds")

##############################################################################################
# Save
##############################################################################################

model_list <- list(G_model, 
                GS_model,
                GSNE_model,
                GR_model,
                GRNE_model,
                GT_model,
                GTNE_model,
                GNE_model,
                N_model
                )

saveRDS(model_list, file = "Inlabru_outputs/model_list.rds")

pred_list <- list(G_pred,
                  GS_pred,
                  GSNE_pred,
                  GR_pred,
                  GRNE_pred,
                  GT_pred,
                  GTNE_pred,
                  GNE_pred,
                  N_pred
)

saveRDS(pred_list, file = "Inlabru_outputs/pred_list.rds")
