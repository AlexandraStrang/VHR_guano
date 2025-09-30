# r script for testing inlabru count models
# creator: Alexandra Strang
# created: 2025

# set working directory
# path to the data
#setwd("/home/stranga/00_nesi_projects/landcare04225/Alexandra_Data/point_process_data/LGCP_nesi_data")
#setwd("E:/Points_test/LGCP_nesi_data")
setwd("D:/PhD_Chap1_part2/LGCP_nesi_data")

# load packages 
library(sf)
library(fmesher)
library(ggplot2)
library(INLA) # version 25.09.04
library(inlabru) # version 2.13.0 (update)
library(dplyr)
library(purrr)
library(tidyr)
library(terra) # for rasters

bru_options_set(control.compute = list(cpo = TRUE))

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

# plot boundary mesh
mesh_plot_Crozier <- ggplot() + 
  geom_fm(data = mesh_sub) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
  ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh_plot_Crozier

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

# have to run this again
percent_guano_raster[is.na(percent_guano_raster)] <- 0

# check for misalignment
covariate_plot <- c(percent_guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)
plot(covariate_plot)

# stack covariates
cov_stack <- c(percent_guano_raster, slope_raster, aspect_raster, roughness_raster, TRI_raster)

# use inlabru to handle abundance model (positive integer response)

# species counts are recorded at each observed location
# the count model is a coarse aggregation of point process model

# rasterise the species counts to match the spatial covariates
# mask the regions outside of the study area (buff coastline boundary)

# prepare response variable
count_raster <- 
  terra::rasterize(vect(sf_Crozier), cov_stack, fun = sum, background = 0) %>%
  terra::aggregate(fact = 5, fun = sum) %>%
  mask(vect(sf::st_geometry(buff_boundary)))

plot(count_raster)

count_raster <- count_raster %>%
  cellSize(unit = "m") %>%
  c(count_raster)

res(count_raster)
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
plot(cov_stack)

# SPDE priors
matern <- inla.spde2.pcmatern(mesh = mesh_sub,
                              prior.range = c(100, 0.9), 
                              prior.sigma = c(0.01, 0.1))

# plot of non-zero counts
ggplot() +
  geom_fm(data = mesh_sub) +
  geom_sf(
    data = counts_df[counts_df$count > 0, ],
    aes(color = count),
    size = 1,
    pch = 4
  ) +
  theme_minimal()

# Poisson GLM
# Poisson model that links species counts per raster cell to spatial covariates

cmps <- ~ percentguano(cov_stack$CrozierGuano_2m, model = "linear") +
  slope(cov_stack$Cape_Crozier_slope, model = "linear") +
  field(geometry, model = matern) + Intercept(1)

fit_poi <- bru(
  cmps,
  bru_obs(
    family = "poisson", data = counts_df,
    formula = count ~
      percentguano + slope +
      field + Intercept,
    E = area
  )
)

summary(fit_poi)

pred_poi <- predict(
  fit_poi, counts_df,
  ~{
    expect <- exp(percentguano +
                    slope +
                    field + Intercept
                  ) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

expect_poi <- pred_poi$expect
expect_poi$pred_var <- expect_poi$mean + expect_poi$sd^2
expect_poi$log_score <- -log(pred_poi$obs_prob$mean)

predicted_abundance <- sum(expect_poi$mean)

ggplot() +
  geom_fm(data = mesh_sub) +
  gg(expect_poi, aes(fill = mean / area), geom = "tile") +
  #geom_sf(data = sf_Crozier, color = "red", size = 1, pch = 4, alpha = 0.2) +
  ggtitle("Nest intensity per ~ m")

# add covariates to cell raster

# Extract covariate values from aggregated cov_stack into counts_df
cov_values <- terra::extract(cov_stack, vect(counts_df))[, -1]
counts_df <- bind_cols(counts_df, cov_values)

# compute cell area from covariate raster
area_raster <- terra::cellSize(cov_stack[[1]], unit = "m")

# extract area values to prediction grid
area_vals <- terra::extract(area_raster, vect(pred_grid))[,2]
pred_grid$area <- area_vals

pred_map <- predict(
  fit_poi, pred_grid,
  ~{
    expect <- exp(percentguano +
                    slope +
                    field + Intercept
    ) * area
    list(
      expect = expect
    )
  },
  n.samples = 1000
)

pred_grid$expected_abundance <- pred_map$expect$mean
pred_raster <- terra::rasterize(vect(pred_grid), cov_stack, field = "expected_abundance")
plot(pred_raster)

total_abundance <- sum(pred_map$expect$mean, na.rm = TRUE)


pred_grid$linear_predictor <- pred_map$percentguano +
  pred_map$slope +
  pred_map$field +
  pred_map$Intercept

hist(pred_grid$linear_predictor)



# ZIP model?
# Zero-Inflated Poisson model

