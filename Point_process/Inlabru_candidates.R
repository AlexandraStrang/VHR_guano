# r script for testing inlabru candidate models
# creator: Alexandra Strang
# created: 2025

# set working directory
setwd("C:/Users/astra/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data")
setwd(setwd("C:/Users/ajs424/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data"))

# load packages 
library(sf)
library(fmesher)
library(ggplot2)
library(INLA) # version 25.09.19
library(inlabru) # for lgcp() version 2.13.0
library(dplyr)
library(purrr)
library(tidyr)
library(terra) # for rasters

bru_options_set(control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))

# sessionInfo() - important to note r version and versions of inlabru, INLA, fmesher
# R version 4.5.0 (2025-04-11 ucrt)

# EPSG 3031 is in metres but when plotting geom_sf will convert to degrees

# read in points from xy csv
# Cape Crozier 2020
Crozier_xy <- read.csv("Crozier_Points_2020_3031.csv")

# convert to sf object and check CRS
sf_Crozier <- st_as_sf(Crozier_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier)

##############################################################################################

# mesh

##############################################################################################

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
# GSA - percent guano + slope + aspect
# GR - percent guano + roughness
# GRA - percent guano + roughness + aspect
# GT - percent guano + TRI
# GTA - percent guano + TRI + aspect
# GA - percent guano + aspect
# N - Null model (spatial field only)

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

saveRDS(GSA_model, file = "Inlabru_outputs/GSA_model.rds")

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

saveRDS(GRA_model, file = "Inlabru_outputs/GRA_model.rds")

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

saveRDS(GTA_model, file = "Inlabru_outputs/GTA_model.rds")

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

saveRDS(GA_model, file = "Inlabru_outputs/GA_model.rds")

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

##############################################################################################

# predictions

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

G_expected <- G_pred$expect
G_expected$pred_var <- G_expected$mean + G_expected$sd^2
G_expected$log_score <- -log(G_pred$obs_prob$mean)

G_abundance <- sum(G_expected$mean)

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

GS_expected <- GS_pred$expect
GS_expected$pred_var <- GS_expected$mean + GS_expected$sd^2
GS_expected$log_score <- -log(GS_pred$obs_prob$mean)

GS_abundance <- sum(GS_expected$mean)

# GSA model predictions
GSA_pred <- predict(
  GSA_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + slope + aspect +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

GSA_expected <- GSA_pred$expect
GSA_expected$pred_var <- GSA_expected$mean + GSA_expected$sd^2
GSA_expected$log_score <- -log(GSA_pred$obs_prob$mean)

GSA_abundance <- sum(GSA_expected$mean)

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

GR_expected <- GR_pred$expect
GR_expected$pred_var <- GR_expected$mean + GR_expected$sd^2
GR_expected$log_score <- -log(GR_pred$obs_prob$mean)

GR_abundance <- sum(GR_expected$mean)

# GRA model predictions
GRA_pred <- predict(
  GRA_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + roughness + aspect +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

GRA_expected <- GRA_pred$expect
GRA_expected$pred_var <- GRA_expected$mean + GRA_expected$sd^2
GRA_expected$log_score <- -log(GRA_pred$obs_prob$mean)

GRA_abundance <- sum(GRA_expected$mean)

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

GT_expected <- GT_pred$expect
GT_expected$pred_var <- GT_expected$mean + GT_expected$sd^2
GT_expected$log_score <- -log(GT_pred$obs_prob$mean)

GT_abundance <- sum(GT_expected$mean)

# GTA model predictions
GTA_pred <- predict(
  GTA_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + TRI + aspect +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

GTA_expected <- GTA_pred$expect
GTA_expected$pred_var <- GTA_expected$mean + GTA_expected$sd^2
GTA_expected$log_score <- -log(GTA_pred$obs_prob$mean)

GTA_abundance <- sum(GTA_expected$mean)

# GA model predictions
GA_pred <- predict(
  GA_model, counts_df,
  ~{
    expect <- exp(Intercept + 
                    percentguano + aspect +
                    field) * area
    list(
      expect = expect,
      obs_prob = dpois(count, expect)
    )
  },
  n.samples = 1000
)

GA_expected <- GA_pred$expect
GA_expected$pred_var <- GA_expected$mean + GA_expected$sd^2
GA_expected$log_score <- -log(GA_pred$obs_prob$mean)

GA_abundance <- sum(GA_expected$mean)

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

N_expected <- N_pred$expect
N_expected$pred_var <- N_expected$mean + N_expected$sd^2
N_expected$log_score <- -log(N_pred$obs_prob$mean)

N_abundance <- sum(N_expected$mean)


##############################################################################################

# summary statistics (WAIC, DIC, Marginal log-Likelihood)

##############################################################################################

model_list <- c(G_model, 
                GS_model,
                GSA_model,
                GR_model,
                GRA_model,
                GT_model,
                GTA_model,
                GA_model,
                N_model
                )

pred_list <- list(G_pred,
                  GS_pred,
                  GSA_pred,
                  GR_pred,
                  GRA_pred,
                  GT_pred,
                  GTA_pred,
                  GA_pred,
                  N_pred
)

saveRDS(pred_list, file = "Inlabru_outputs/pred_list.rds")

# loop through models
for (i in model_list) {
  model <- i
  
  # NOT WORKING
  
  fixed_effects <- model$names.fixed
  
  print(fixed_effects)
  
  # grab summary statistics
  Intercept_mean <- model$summary.fixed$mean[1]
  Intercept_sd <- model$summary.fixed$sd[1]
  Intercept_0.025 <- model$summary.fixed$`0.025quant`[1]
  Intercept_0.975 <- model$summary.fixed$`0.975quant`[1]
  
  # grab diagnostic statistics
  WAIC <- model$waic$waic
  DIC <- model$dic$dic
  MLik <- model$mlik[2]
}


# table of predictions


# plot predictions
G_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(G_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per ~ m")
ggsave("Inlabru_outputs/G_predictions.png", G_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
       )

GS_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GS_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per ~ m")
ggsave("Inlabru_outputs/GS_predictions.png", GS_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GSA_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GSA_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per ~ m")
ggsave("Inlabru_outputs/GSA_predictions.png", GSA_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GR_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GR_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per ~ m")
ggsave("Inlabru_outputs/GR_predictions.png", GR_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GRA_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GRA_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per ~ m")
ggsave("Inlabru_outputs/GRA_predictions.png", GRA_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GT_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GT_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per ~ m")
ggsave("Inlabru_outputs/GT_predictions.png", GT_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GTA_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GTA_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per ~ m")
ggsave("Inlabru_outputs/GTA_predictions.png", GTA_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GA_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GA_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per ~ m")
ggsave("Inlabru_outputs/GA_predictions.png", GA_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

N_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(N_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per ~ m")
ggsave("Inlabru_outputs/N_predictions.png", N_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
) # shows up with nothing

# Plot residuals

