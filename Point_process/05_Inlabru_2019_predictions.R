# r script for using inlabru models to predict 2019 count
# creator: Alexandra Strang
# created: 2026

# set working directory
setwd("C:/Users/ajs424/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data")

# load packages 
library(sf)
library(ggplot2)
library(INLA) # version 25.09.19
library(inlabru) # for bru() version 2.13.0
library(dplyr)
library(terra) # for rasters
library(ggpubr)
library(scoringRules) # for CRPS

bru_options_set(control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))

# sessionInfo() - important to note r version and versions of inlabru, INLA, fmesher
# R version 4.5.1

##############################################################################################
# Load data
##############################################################################################

# mesh and coastline boundary
mesh_sub <- readRDS("Inlabru_outputs/mesh_sub.rds")
buff_boundary <- readRDS("Inlabru_outputs/buff_boundary.rds")

# 2020 count dataframe created in inlabru candidates
counts_df     <- readRDS("Inlabru_outputs/counts_df.rds")

##############################################################################################
# Load models
##############################################################################################

G_model    <- readRDS("Inlabru_outputs/G_model.rds")
GS_model   <- readRDS("Inlabru_outputs/GS_model.rds")
GSNE_model <- readRDS("Inlabru_outputs/GSNE_model.rds")
GR_model   <- readRDS("Inlabru_outputs/GR_model.rds")
GRNE_model <- readRDS("Inlabru_outputs/GRNE_model.rds")
GT_model   <- readRDS("Inlabru_outputs/GT_model.rds")
GTNE_model <- readRDS("Inlabru_outputs/GTNE_model.rds")
GNE_model  <- readRDS("Inlabru_outputs/GNE_model.rds")
N_model    <- readRDS("Inlabru_outputs/N_model.rds")

##############################################################################################
# Original covariates
##############################################################################################

# add covariates at 2m res (aggregated below)

# continuous guano raster (2020)
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
# Get 2019 guano area
##############################################################################################

# make 2019 newdata (2019 guano area but same terrain)

# use continuous 2019 guano raster
percent_guano_2019_raster <- rast("2019_CrozierGuano_2m.tif") # select 2019 data

# change guano crs to match 
crs(percent_guano_2019_raster) <- "EPSG:3031"
crs(percent_guano_2019_raster)

res(percent_guano_2019_raster) # 2m x 2m

# match 2019 guano area extent to terrain rasters
percent_guano_2019_raster <- extend(percent_guano_2019_raster, slope_raster)

percent_guano_2019_raster[values(percent_guano_2019_raster) > 1] <- 0
percent_guano_2019_raster[is.na(percent_guano_2019_raster)] <- 0

# Extract values from raster
vals_2019 <- values(percent_guano_2019_raster)

# Remove NAs and filter out 0s and 1s
filtered_vals_2019 <- vals_2019[!is.na(vals_2019)]
mean(filtered_vals_2019)

# scale 2019 guano to 2020 guano values (mean and sd)
standardize_2020 <- function(guano_2020, guano_2019) {
  m <- global(guano_2020, fun = "mean", na.rm = TRUE)[1, 1]
  s <- global(guano_2020, fun = "sd", na.rm = TRUE)[1, 1]
  (guano_2019 - m) / s
}

# apply to 2019 guano
percent_guano_2019_raster  <- standardize_2020(percent_guano_raster, percent_guano_2019_raster)

summary(percent_guano_2019_raster)

# have to run this again
percent_guano_2019_raster[is.na(percent_guano_2019_raster)] <- 0

# check for misalignment
covariate_plot <- c(percent_guano_2019_raster, slope_raster, northness_raster, eastness_raster, roughness_raster, TRI_raster)
plot(covariate_plot)

# stack covariates
cov_stack2 <- c(percent_guano_2019_raster, slope_raster, northness_raster, eastness_raster, roughness_raster, TRI_raster)

# check for colinearity between covariates again
cov_values2 <- as.data.frame(cov_stack2, na.rm = TRUE)
cor(cov_values2)

##############################################################################################
# Predict for 2019
##############################################################################################

# predict 2019 counts with 2020 fitted models

# 2019 guano scaled with 2020 mean and sd

# swap for 2019 guano (2020 standardised)
cov_stack$CrozierGuano_2m <- percent_guano_2019_raster

geom_data <- st_as_sf(counts_df) %>% select(geometry, area)

# 2019 G model predictions
G_pred_2019 <- predict(
  G_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

saveRDS(G_pred_2019, file = "Inlabru_outputs/G_pred_2019.rds")


# 2019 GS model predictions
GS_pred_2019 <- predict(
  GS_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + slope + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

saveRDS(GS_pred_2019, file = "Inlabru_outputs/GS_pred_2019.rds")


# 2019 GSNE model predictions
GSNE_pred_2019 <- predict(
  GSNE_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + slope + northness + eastness + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

saveRDS(GSNE_pred_2019, file = "Inlabru_outputs/GSNE_pred_2019.rds")


# 2019 GR model predictions
GR_pred_2019 <- predict(
  GR_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + roughness + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

saveRDS(GR_pred_2019, file = "Inlabru_outputs/GR_pred_2019.rds")


# 2019 GRNE model predictions
GRNE_pred_2019 <- predict(
  GRNE_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + roughness + northness + eastness + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

saveRDS(GRNE_pred_2019, file = "Inlabru_outputs/GRNE_pred_2019.rds")


# 2019 GT model predictions
GT_pred_2019 <- predict(
  GT_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + TRI + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

saveRDS(GT_pred_2019, file = "Inlabru_outputs/GT_pred_2019.rds")


# 2019 GTNE model predictions
GTNE_pred_2019 <- predict(
  GTNE_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + TRI + northness + eastness + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

saveRDS(GTNE_pred_2019, file = "Inlabru_outputs/GTNE_pred_2019.rds")


# 2019 GNE model predictions
GNE_pred_2019 <- predict(
  GNE_model,
  newdata = geom_data,
  ~{
    eta <- Intercept + percentguano + northness + eastness + field
    mu  <- exp(eta) * area
    list(expect = mu)
  },
  n.samples = 1000
)

saveRDS(GNE_pred_2019, file = "Inlabru_outputs/GNE_pred_2019.rds")


# calculate expected counts and overall abundance

G_expected_2019 <- G_pred_2019$expect
G_abundance_2019 <- sum(G_expected_2019$mean)

GS_expected_2019 <- GS_pred_2019$expect
GS_abundance_2019 <- sum(GS_expected_2019$mean)

GSNE_expected_2019 <- GSNE_pred_2019$expect
GSNE_abundance_2019 <- sum(GSNE_expected_2019$mean)

GR_expected_2019 <- GR_pred_2019$expect
GR_abundance_2019 <- sum(GR_expected_2019$mean)

GRNE_expected_2019 <- GRNE_pred_2019$expect
GRNE_abundance_2019 <- sum(GRNE_expected_2019$mean)

GT_expected_2019 <- GT_pred_2019$expect
GT_abundance_2019 <- sum(GT_expected_2019$mean)

GTNE_expected_2019 <- GTNE_pred_2019$expect
GTNE_abundance_2019 <- sum(GTNE_expected_2019$mean)

GNE_expected_2019 <- GNE_pred_2019$expect
GNE_abundance_2019 <- sum(GNE_expected_2019$mean)

# add 2019 abundance predictions
abundance_df_2019 <- data.frame(
  Model = c("G", "GS", "GSNE", "GR", "GRNE", "GT", "GTNE", "GNE"),
  predicted_abundance = c(G_abundance_2019, GS_abundance_2019, GSNE_abundance_2019, GR_abundance_2019,
                          GRNE_abundance_2019, GT_abundance_2019, GTNE_abundance_2019, GNE_abundance_2019)
)

# abundance 2019 predictions plot
a_pred_df_2019 <- abundance_df_2019 %>%
  select(Model, predicted_abundance) %>%
  distinct()

# get 2019 points from xy csv
# Cape Crozier 2019
Crozier_xy_2019 <- read.csv("Crozier_UAV_points/Reprojected_3031/Crozier_2019_Points_3031.csv")
observed_n_2019 <- nrow(Crozier_xy_2019)
# 245,918 breeding pairs

abundance_plot_2019 <- ggplot(a_pred_df_2019, aes(x = Model, y = predicted_abundance)) +
  geom_point(size = 3, color = "black") +
  geom_hline(aes(yintercept = observed_n_2019, linetype = "Observed"), color = "red") +
  scale_linetype_manual(name = "", values = c("Observed" = "dashed")) +
  labs(
    x = "Model",
    y = "Predicted count (BP)"
  ) +
  scale_y_continuous(limits = c(245000, 350000)) +
  theme_minimal() +
  theme(legend.position = "right") 

abundance_plot_2019

ggsave("Inlabru_outputs/Predicted_count_2019.png", abundance_plot_2019,
       width = 8, height = 5, units = "in",
       dpi = 600
)

# difference from observed plot 2019
diff_df_2019 <- abundance_df_2019 %>%
  select(Model, predicted_abundance) %>%
  distinct() %>%
  mutate(difference = predicted_abundance - observed_n_2019)

# Plot
difference_plot_2019 <- ggplot(diff_df_2019, aes(x = Model, y = difference)) +
  geom_point(size = 3, color = "black") +
  labs(
    x = "Model",
    y = "Predicted - observed count (BP)") +
  scale_y_continuous(limits = c(80000, 85000)) +
  theme_minimal()

difference_plot_2019

ggsave("Inlabru_outputs/Difference_plot_2019.png", difference_plot_2019,
       width = 8, height = 5, units = "in",
       dpi = 600
)

# save 2019 difference df
write.csv(diff_df_2019, file = "Inlabru_outputs/2019_predictions.csv", row.names = FALSE)

# compute CRPS for 2019 predictions

# get 2019 observed counts df

# convert 2019 to sf object and check CRS
sf_Crozier_2019 <- st_as_sf(Crozier_xy_2019, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier_2019)

# note cov_stack contains 2019 guano
count_raster_2019 <- 
  terra::rasterize(vect(sf_Crozier_2019), cov_stack, fun = sum, background = 0) %>%
  terra::aggregate(fact = 5, fun = sum) %>%
  mask(vect(sf::st_geometry(buff_boundary)))

# counts of nests
plot(count_raster_2019)

count_raster_2019 <- count_raster_2019 %>%
  cellSize(unit = "m") %>%
  c(count_raster_2019)

res(count_raster_2019) # 10 x 10 m
summary(count_raster_2019)
# 148 penguins in 100 m2

# extract the coordinates for these pixels
counts_df_2019 <- crds(count_raster_2019, df = TRUE, na.rm = TRUE) %>%
  bind_cols(values(count_raster_2019, mat = TRUE, na.rm = TRUE)) %>%
  rename(count = sum) %>%
  mutate(present = (count > 0) *1L) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(sf_Crozier_2019))

# extract coordinates from sf geometry
counts_df_2019 <- counts_df_2019 %>%
  mutate(
    x_coord = st_coordinates(.)[, 1],
    y_coord = st_coordinates(.)[, 2]
  )

CRPS_list_2019 <- list(G_pred_2019,
                       GS_pred_2019,
                       GSNE_pred_2019,
                       GR_pred_2019,
                       GRNE_pred_2019,
                       GT_pred_2019,
                       GTNE_pred_2019,
                       GNE_pred_2019
)

# observed counts
observed_counts_2019 <- counts_df_2019$count #y

# crps_pois(y = vector of observations, lambda = vector of non-negative means)
compute_crps_2019 <- function(pred_obj) {
  predicted_means <- pred_obj$expect$mean
  crps_values <- crps_pois(y = observed_counts_2019, lambda = predicted_means)
  as.list(summary(crps_values))
}

model_names <- c("G", "GS", "GSNE", "GR", "GRNE", "GT", "GTNE", "GNE")
crps_summary_list_2019 <- lapply(CRPS_list_2019, compute_crps_2019)

# table of crps values
crps_table_2019 <- bind_rows(crps_summary_list_2019) %>%
  mutate(Model = model_names) %>%
  select(Model, Min., `1st Qu.`, Median, Mean, `3rd Qu.`, Max.)

crps_plot_2019 <- ggtexttable(crps_table_2019, rows = NULL)
plot(crps_plot_2019)

# could use R2 predicted for predicting at Cape Crozier 2019
