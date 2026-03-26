# r script for testing inlabru candidate models
# creator: Alexandra Strang
# created: 2026

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

# get 2019 points from xy csv
# Cape Crozier 2019
Crozier_xy_2019 <- read.csv("Crozier_UAV_points/Reprojected_3031/Crozier_2019_Points_3031.csv")

observed_n_2019 <- nrow(Crozier_xy_2019)
# 245,918 breeding pairs

# convert to sf object and check CRS
sf_Crozier <- st_as_sf(Crozier_xy_2019, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier)

# mesh and coastline boundary
mesh_sub <- readRDS("Inlabru_outputs/mesh_sub.rds")
buff_boundary <- readRDS("Inlabru_outputs/buff_boundary.rds")

# plot to check overlap
# plot mesh with coastline boundary
plot_Crozier <- ggplot() + 
  geom_fm(data = mesh_sub) + 
  geom_sf(data = buff_boundary, fill = NA, color = "blue") +
  geom_sf(data = sf_Crozier, color = "purple", size = 1.7, alpha = 0.5) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
  ) + 
  theme_minimal()
# geom_sf converts to degrees

plot_Crozier

##############################################################################################
# covariates
##############################################################################################

# add covariates at 2m res (aggregated below)

# use continuous 2019 guano raster
percent_guano_2019_raster <- rast("2019_CrozierGuano_2m.tif") # select 2019 data

# change guano crs to match 
crs(percent_guano_2019_raster) <- "EPSG:3031"
crs(percent_guano_2019_raster)

res(percent_guano_2019_raster) # 2m x 2m

# 2m terrain rasters clipped by mesh
slope_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_slope.tif")
northness_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_northness.tif")
eastness_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_eastness.tif")
roughness_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_Roughness.tif")
TRI_raster <- rast("Crozier_terrain_mesh/Cape_Crozier_TRI.tif")

# match 2019 guano area extent to terrain rasters
percent_guano_2019_raster <- extend(percent_guano_2019_raster, slope_raster)

percent_guano_2019_raster[values(percent_guano_2019_raster) > 1] <- 0
percent_guano_2019_raster[is.na(percent_guano_2019_raster)] <- 0

# Extract values from raster
vals_2019 <- values(percent_guano_2019_raster)

# Remove NAs and filter out 0s and 1s
filtered_vals_2019 <- vals_2019[!is.na(vals_2019)]
mean(filtered_vals_2019)

# scale
# (mean of 0 sd of 1)
standardize <- function(r) {
  m <- global(r, fun = "mean", na.rm = TRUE)[1, 1]
  s <- global(r, fun = "sd", na.rm = TRUE)[1, 1]
  (r - m) / s
}

# apply to continuous variables
percent_guano_2019_raster   <- standardize(percent_guano_2019_raster)
slope_raster     <- standardize(slope_raster)
northness_raster    <- standardize(northness_raster)
eastness_raster    <- standardize(eastness_raster)
roughness_raster <- standardize(roughness_raster)
TRI_raster       <- standardize(TRI_raster)

summary(percent_guano_2019_raster)

# check for misalignment
covariate_plot <- c(percent_guano_2019_raster, slope_raster, northness_raster, eastness_raster, roughness_raster, TRI_raster)
plot(covariate_plot)

# stack covariates
cov_stack <- c(percent_guano_2019_raster, slope_raster, northness_raster, eastness_raster, roughness_raster, TRI_raster)

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
# 148 penguins in 100 m2

# sum of area that has zero penguins
cell_areas <- count_raster[[1]]
zero_cells <- count_raster[[2]] == 0
zero_area <- mask(cell_areas, zero_cells, maskvalue = FALSE)
total_zero_area <- global(zero_area, fun = "sum", na.rm = TRUE)
print(total_zero_area)

total_raster_area <- global(cell_areas, fun = "sum", na.rm = TRUE)
print(total_raster_area)

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

# aggregate covariates to match count raster
cov_stack <- terra::aggregate(cov_stack, fact = 5, fun = mean)
res(cov_stack) # 10 x 10 m

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

# N - Null model (spatial field only)

# northness and eastness have some NAs

# Guano model
G_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$`2019_CrozierGuano_2m`, model = "linear") +
  field(geometry, model = matern)

G_2019_model <- bru(
  G_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + 
      field,
    E = area
  )
)

summary(G_2019_model)

# Guano + slope model
GS_cmp <- ~ Intercept(1) +
  percentguano(cov_stack$`2019_CrozierGuano_2m`, model = "linear") +
  slope(cov_stack$Cape_Crozier_slope, model = "linear") +
  field(geometry, model = matern)

GS_2019_model <- bru(
  GS_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      percentguano + slope +
      field,
    E = area
  )
)

summary(GS_2019_model)


# Null model
N_cmp <- ~ Intercept(1) +
  field(geometry, model = matern)

N_2019_model <- bru(
  N_cmp,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~ Intercept +
      field,
    E = area
  )
)

summary(N_2019_model)

#############################################################################################
# Predictions
##############################################################################################

# G model predictions
G_2019_pred <- predict(
  G_2019_model, counts_df,
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

# abundance
G_expected_2019 <- G_2019_pred$expect
G_abundance_2019 <- sum(G_expected_2019$mean)

G_abundance_2019 - observed_n_2019
# 16,794.4 breeding pairs off

# GS model predictions
GS_2019_pred <- predict(
  GS_2019_model, counts_df,
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

# abundance
GS_expected_2019 <- GS_2019_pred$expect
GS_abundance_2019 <- sum(GS_expected_2019$mean)

GS_abundance_2019 - observed_n_2019
# 17,847.1 breeding pairs off

# N model predictions
N_2019_pred <- predict(
  N_2019_model, counts_df,
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

# abundance
N_expected_2019 <- N_2019_pred$expect
N_abundance_2019 <- sum(N_expected_2019$mean)

N_abundance_2019 - observed_n_2019
# 91,481.06 breeding pairs off

##############################################################################################
# Get 2020 guano area
##############################################################################################

# Need to use 2019 to predict 2020 now instead

# use continuous 2019 guano raster
percent_guano_2020_raster <- rast("CrozierGuano_2m.tif") # select 2020 data

# change guano crs to match 
crs(percent_guano_2020_raster) <- "EPSG:3031"
crs(percent_guano_2020_raster)

res(percent_guano_2020_raster) # 2m x 2m

# match 2020 guano area extent to terrain rasters
percent_guano_2020_raster <- extend(percent_guano_2020_raster, slope_raster)

percent_guano_2020_raster[values(percent_guano_2020_raster) > 1] <- 0
percent_guano_2020_raster[is.na(percent_guano_2020_raster)] <- 0

# Extract values from raster
vals_2020 <- values(percent_guano_2020_raster)

# Remove NAs and filter out 0s and 1s
filtered_vals_2020 <- vals_2020[!is.na(vals_2020)]
mean(filtered_vals_2020)

# use raw continuous 2019 guano raster for standardisation
raw_percent_guano_2019_raster <- rast("2019_CrozierGuano_2m.tif")
# change guano crs to match 
crs(raw_percent_guano_2019_raster) <- "EPSG:3031"
crs(raw_percent_guano_2019_raster)
res(raw_percent_guano_2019_raster) # 2m x 2m
# match 2019 guano area extent to terrain rasters
raw_percent_guano_2019_raster <- extend(raw_percent_guano_2019_raster, slope_raster)
raw_percent_guano_2019_raster[values(raw_percent_guano_2019_raster) > 1] <- 0
raw_percent_guano_2019_raster[is.na(raw_percent_guano_2019_raster)] <- 0

# scale 2020 guano to 2019 guano values (mean and sd)
standardize_2020 <- function(guano_2020, guano_2019) {
  m <- global(guano_2019, fun = "mean", na.rm = TRUE)[1, 1]
  s <- global(guano_2019, fun = "sd", na.rm = TRUE)[1, 1]
  (guano_2020 - m) / s
}

# apply to 2020 guano
percent_guano_2020_raster  <- standardize_2020(percent_guano_2020_raster, raw_percent_guano_2019_raster)
summary(percent_guano_2020_raster)

# check for misalignment
covariate_plot_2020 <- c(percent_guano_2020_raster, slope_raster, northness_raster, eastness_raster, roughness_raster, TRI_raster)
plot(covariate_plot_2020)

# stack covariates
cov_stack2 <- c(percent_guano_2020_raster, slope_raster, northness_raster, eastness_raster, roughness_raster, TRI_raster)

# aggregate covariates to match count raster
cov_stack2 <- terra::aggregate(cov_stack2, fact = 5, fun = mean)
# 2020 guano and 2019 guano extents need to match (with terrain covariates also)

##############################################################################################
# Predict for 2020
##############################################################################################

# read in points from xy csv
# Cape Crozier 2020
Crozier_xy_2020 <- read.csv("Crozier_UAV_points/Reprojected_3031/Crozier_2020_Points_3031.csv")
observed_n_2020 <- nrow(Crozier_xy_2020)

# predict 2020 counts with 2019 fitted models

# 2020 guano scaled with 2019 mean and sd

# swap for 2020 guano (standardised with 2019)
cov_stack$`2019_CrozierGuano_2m` <- cov_stack2$CrozierGuano_2m

# from 2019 counts
geom_data <- st_as_sf(counts_df) %>% select(geometry, area)

# 2020 G model predictions
G_pred_2020 <- predict(
  G_2019_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

# abundance
G_expected_2020 <- G_pred_2020$expect
G_abundance_2020 <- sum(G_expected_2020$mean)

G_abundance_2020 - observed_n_2020


# 2020 GS model predictions
GS_pred_2020 <- predict(
  GS_2019_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + slope + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

# abundance
GS_expected_2020 <- GS_pred_2020$expect
GS_abundance_2020 <- sum(GS_expected_2020$mean)

GS_abundance_2020 - observed_n_2020

