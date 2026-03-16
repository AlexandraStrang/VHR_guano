# r script for testing inlabru candidate models
# creator: Alexandra Strang
# created: 2025

# set working directory
#setwd("C:/Users/astra/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data")
setwd("C:/Users/ajs424/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data")

# load packages 
library(sf)
library(fmesher)
library(ggplot2)
library(INLA) # version 25.09.19
library(inlabru) # for bru() version 2.13.0
library(dplyr)
library(purrr)
library(tidyr)
library(terra) # for rasters
library(ggpubr)
library(scoringRules) # for CRPS

bru_options_set(control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))

# sessionInfo() - important to note r version and versions of inlabru, INLA, fmesher
# R version 4.5.1

# EPSG 3031 is in metres but when plotting geom_sf will convert to degrees

# read in points from xy csv
# Cape Crozier 2020
Crozier_xy <- read.csv("Crozier_UAV_points/Reprojected_3031/Crozier_2020_Points_3031.csv")

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
# 2816

# subdivide mesh
# splits triangles into subtriangles
mesh_sub <- fm_subdivide(Crozier_mesh,3)
print(mesh_sub$n)
# 44519

# plot mesh with coastline boundary
mesh_plot_Crozier <- ggplot() + 
  geom_fm(data = mesh_sub) + 
  geom_sf(data = buff_boundary, fill = NA, color = "blue") +
  geom_sf(data = sf_Crozier, color = "purple", size = 1.7, alpha = 0.5) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
  ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh_plot_Crozier
ggsave("Inlabru_outputs/Mesh_plot.png", mesh_plot_Crozier,
       width = 8, height = 5, units = "in",
       dpi = 600
)

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
# summary statistics (Abundance, residuals, CRPS)
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

# calculate expected counts, log-scores and overall abundance

G_expected <- G_pred$expect
G_expected$pred_var <- G_expected$mean + G_expected$sd^2
G_expected$log_score <- -log(G_pred$obs_prob$mean)

G_abundance <- sum(G_expected$mean)


GS_expected <- GS_pred$expect
GS_expected$pred_var <- GS_expected$mean + GS_expected$sd^2
GS_expected$log_score <- -log(GS_pred$obs_prob$mean)

GS_abundance <- sum(GS_expected$mean)


GSNE_expected <- GSNE_pred$expect
GSNE_expected$pred_var <- GSNE_expected$mean + GSNE_expected$sd^2
GSNE_expected$log_score <- -log(GSNE_pred$obs_prob$mean)

GSNE_abundance <- sum(GSNE_expected$mean)

GR_expected <- GR_pred$expect
GR_expected$pred_var <- GR_expected$mean + GR_expected$sd^2
GR_expected$log_score <- -log(GR_pred$obs_prob$mean)

GR_abundance <- sum(GR_expected$mean)


GRNE_expected <- GRNE_pred$expect
GRNE_expected$pred_var <- GRNE_expected$mean + GRNE_expected$sd^2
GRNE_expected$log_score <- -log(GRNE_pred$obs_prob$mean)

GRNE_abundance <- sum(GRNE_expected$mean)


GT_expected <- GT_pred$expect
GT_expected$pred_var <- GT_expected$mean + GT_expected$sd^2
GT_expected$log_score <- -log(GT_pred$obs_prob$mean)

GT_abundance <- sum(GT_expected$mean)


GTNE_expected <- GTNE_pred$expect
GTNE_expected$pred_var <- GTNE_expected$mean + GTNE_expected$sd^2
GTNE_expected$log_score <- -log(GTNE_pred$obs_prob$mean)

GTNE_abundance <- sum(GTNE_expected$mean)


GNE_expected <- GNE_pred$expect
GNE_expected$pred_var <- GNE_expected$mean + GNE_expected$sd^2
GNE_expected$log_score <- -log(GNE_pred$obs_prob$mean)

GNE_abundance <- sum(GNE_expected$mean)


N_expected <- N_pred$expect
N_expected$pred_var <- N_expected$mean + N_expected$sd^2
N_expected$log_score <- -log(N_pred$obs_prob$mean)

N_abundance <- sum(N_expected$mean)


# extract summaries and results
results_list <- list()

model_names <- c("G", "GS", "GSNE", "GR", "GRNE", "GT", "GTNE", "GNE", "N")

# loop through models
for (i in seq_along(model_list)) {
  model <- model_list[[i]]
  model_name <- model_names[i]

  # print model name
  print(model_name)

  # hyperparameters and diagnostics
    Range_Mean <- model$summary.hyperpar$mean[1]
    Range_sd <- model$summary.hyperpar$sd[1]
    Range_0.025 <- model$summary.hyperpar$`0.025quant`[1]
    Range_0.975 <- model$summary.hyperpar$`0.975quant`[1]
    
    Stdev_Mean <- model$summary.hyperpar$mean[2]
    Stdev_sd <- model$summary.hyperpar$sd[2]
    Stdev_0.025 <- model$summary.hyperpar$`0.025quant`[2]
    Stdev_0.975 <- model$summary.hyperpar$`0.975quant`[2]
    
    WAIC <- model$waic$waic
    DIC <- model$dic$dic
    MLik <- model$mlik[2]
    
    # fixed effects
    for (j in seq_along(model$names.fixed)) {
      effect_name <- model$names.fixed[j]
      
      effect_row <- data.frame(
        Model = model_name,
        Effect = effect_name,
        Mean = model$summary.fixed$mean[j],
        SD = model$summary.fixed$sd[j],
        CI_0.025 = model$summary.fixed$`0.025quant`[j],
        CI_0.975 = model$summary.fixed$`0.975quant`[j],
        Range_Mean = Range_Mean,
        Range_SD = Range_sd,
        Range_0.025 = Range_0.025,
        Range_0.975 = Range_0.975,
        Stdev_Mean = Stdev_Mean,
        Stdev_SD = Stdev_sd,
        Stdev_0.025 = Stdev_0.025,
        Stdev_0.975 = Stdev_0.975,
        WAIC = WAIC,
        DIC = DIC,
        MLik = MLik,
        stringsAsFactors = FALSE
      )
      
      results_list[[length(results_list) + 1]] <- effect_row
    }
}
  
results_df <- do.call(rbind, results_list)

# add abundance predictions
abundance_df <- data.frame(
  Model = c("G", "GS", "GSNE", "GR", "GRNE", "GT", "GTNE", "GNE", "N"),
  predicted_abundance = c(G_abundance, GS_abundance, GSNE_abundance, GR_abundance,
                  GRNE_abundance, GT_abundance, GTNE_abundance, GNE_abundance, N_abundance)
)

# save abundance df
write.csv(abundance_df, file = "Inlabru_outputs/Abundance_predictions.csv", row.names = FALSE)

results_df <- merge(results_df, abundance_df, by = "Model", all.x = TRUE)

# save model outputs
write.csv(results_df, file = "Inlabru_outputs/Candidate_results.csv", row.names = FALSE)

# abundance predictions plot
a_pred_df <- abundance_df %>%
  filter(Model != "N") %>%
  select(Model, predicted_abundance) %>%
  distinct()

observed_n <- nrow(Crozier_xy)
# 249,007 UAV xy points

abundance_plot <- ggplot(a_pred_df, aes(x = Model, y = predicted_abundance)) +
  geom_point(size = 3, color = "black") +
  geom_hline(aes(yintercept = observed_n, linetype = "Observed"), color = "red") +
  scale_linetype_manual(name = "", values = c("Observed" = "dashed")) +
  labs(
    x = "Model",
    y = "Predicted count (BP)"
    ) +
  scale_y_continuous(limits = c(247500, 257500)) +
  theme_minimal() +
  theme(legend.position = "right") 

abundance_plot

ggsave("Inlabru_outputs/Predicted_count.png", abundance_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

# difference from observed plot
diff_df <- abundance_df %>%
  filter(Model != "N") %>%
  select(Model, predicted_abundance) %>%
  distinct() %>%
  mutate(difference = predicted_abundance - observed_n)

# Plot
difference_plot <- ggplot(diff_df, aes(x = Model, y = difference)) +
  geom_point(size = 3, color = "black") +
  labs(
    x = "Model",
    y = "Predicted - observed count (BP)") +
  scale_y_continuous(limits = c(5500, 6500)) +
  theme_minimal()

difference_plot

ggsave("Inlabru_outputs/Difference_plot.png", difference_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

# fixed effect plot
models_to_plot <- subset(results_df, Model != "N")

plot_list <- list()

for (m in unique(models_to_plot$Model)) {
  df_model <- models_to_plot %>% 
    filter(Model == m) %>%
    mutate(
      crosses_zero = factor(CI_0.025 <= 0 & CI_0.975 >= 0, levels = c(TRUE, FALSE)))
  
  p <- ggplot(df_model, aes(x = Mean, y = Effect)) +
    geom_point(aes(color = crosses_zero), size = 3) +
    geom_errorbar(aes(xmin = CI_0.025, xmax = CI_0.975, color = crosses_zero), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      labels = c("TRUE" = "TRUE", "FALSE" = "FALSE")
    ) + 
    labs(color = "Includes 0") +
    labs(title = paste(m, "Model"), x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      legend.position = "right" 
    )
  
  plot_list[[m]] <- p
}

effects_plots <- ggpubr::ggarrange(plotlist = plot_list, ncol = 4, nrow = 2,
                                   common.legend = TRUE,
                                   legend = "right")
effects_plots <- annotate_figure(effects_plots,
    left = text_grob("Effect", rot = 90, vjust = 1, face = "bold"),
    bottom = text_grob("Mean", face = "bold"))

effects_plots

ggexport(effects_plots, filename = "Inlabru_outputs/Effects_plots.png",
       width = 9600, height = 4800, units = "in",
       res = 600
)

# plot predictions
G_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(G_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per m2")

ggsave("Inlabru_outputs/G_predictions.png", G_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
       )

GS_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GS_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per m2")

ggsave("Inlabru_outputs/GS_predictions.png", GS_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GSNE_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GSNE_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per m2")

ggsave("Inlabru_outputs/GSNE_predictions.png", GSNE_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GR_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GR_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per m2")

ggsave("Inlabru_outputs/GR_predictions.png", GR_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GRNE_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GRNE_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per m2")

ggsave("Inlabru_outputs/GRNE_predictions.png", GRNE_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GT_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GT_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per m2")

ggsave("Inlabru_outputs/GT_predictions.png", GT_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GTNE_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GTNE_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per m2")

ggsave("Inlabru_outputs/GTNE_predictions.png", GTNE_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

GNE_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(GNE_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per m2")

ggsave("Inlabru_outputs/GNE_predictions.png", GNE_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

N_plot <- ggplot() +
  geom_fm(data = mesh_sub) +
  gg(N_expected, aes(fill = mean / area), geom = "tile") +
  ggtitle("Nest intensity per m2")

ggsave("Inlabru_outputs/N_predictions.png", N_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
) # shows up with nothing

Intensity_plots <- ggpubr::ggarrange(G_plot,
                                     GS_plot,
                                     GSNE_plot,
                                     GR_plot,
                                     GRNE_plot,
                                     GT_plot,
                                     GTNE_plot,
                                     GNE_plot,
                                     ncol = 2, nrow = 4, 
                                     labels=c("G","GS","GSNE","GR","GRNE","GT",
                                              "GTNE","GNE"))

Intensity_plots

ggsave("Inlabru_outputs/Intensity_plots.png", Intensity_plots,
       width = 7, height = 7.8, units = "in",
       dpi = 600
)

# Compute and plot cell-wise Pearson residuals for each model

residuals_list <- list()
resid_plots <- list()

model_resids_list <- list(G_model, 
                   GS_model,
                   GSNE_model,
                   GR_model,
                   GRNE_model,
                   GT_model,
                   GTNE_model,
                   GNE_model
)

# loop through models
for (i in seq_along(model_resids_list)) {
  model <- model_resids_list[[i]]
  model_name <- model_names[i]
  pred <- pred_list[[i]]
  
  # print model name
  print(model_name)

  # compute residuals
  resids <- counts_df %>%
    mutate(
      fitted = pred$expect$mean,
      residual = (count - fitted)/ sqrt(fitted),
      model = model_name
    )
    
  residuals_list[[model_name]] <- resids
  
  # plot residuals
  p <- ggplot() +
    gg(data = mesh_sub) +
    geom_point(data = resids, aes(x = x_coord, y = y_coord, color = residual), size = 2) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    coord_fixed() +
    labs(title = paste(model_name),
          color = "Residuals") +
    theme_minimal()
  
  assign(paste0(model_name, "_resids"), p)
  resid_plots[[model_name]] <- p
  
  # save each plot individually
  ggsave(
    filename = paste0("Inlabru_outputs/", model_name, "_resids.png"),
    plot = p,
    width = 8, height = 5, units = "in", dpi = 600
  )
  
}

resid_plots <- ggpubr::ggarrange(plotlist = resid_plots, ncol = 2, nrow = 4)

resid_plots

ggsave("Inlabru_outputs/Resid_plots.png", resid_plots,
       width = 7, height = 8, units = "in",
       dpi = 600
)


# model evaluation: CRPS (continuous ranked probability scores)
# use scoringRules package

CRPS_list <- list(G_pred,
                  GS_pred,
                  GSNE_pred,
                  GR_pred,
                  GRNE_pred,
                  GT_pred,
                  GTNE_pred,
                  GNE_pred
)

# observed counts
observed_counts <- counts_df$count #y

# crps_pois(y = vector of observations, lambda = vector of non-negative means)
compute_crps <- function(pred_obj) {
  predicted_means <- pred_obj$expect$mean
  crps_values <- crps_pois(y = observed_counts, lambda = predicted_means)
  as.list(summary(crps_values))
}

model_names <- c("G", "GS", "GSNE", "GR", "GRNE", "GT", "GTNE", "GNE")
crps_summary_list <- lapply(CRPS_list, compute_crps)

# table of crps values
crps_table <- bind_rows(crps_summary_list) %>%
  mutate(Model = model_names) %>%
  select(Model, Min., `1st Qu.`, Median, Mean, `3rd Qu.`, Max.)

crps_plot <- ggtexttable(crps_table, rows = NULL)
plot(crps_plot)

# add log-score values for checking predictive performance

logscore_list <- list(G_expected,
                      GS_expected,
                      GSNE_expected,
                      GR_expected,
                      GRNE_expected,
                      GT_expected,
                      GTNE_expected,
                      GNE_expected
)

compute_logscore <- function(expected_obj) {
  mean_logscore <- mean(expected_obj$log_score)
  as.list(mean_logscore)
}

mean_logscores <- lapply(logscore_list, compute_logscore)

logscore_vector <- unlist(mean_logscores)

crps_table <- bind_rows(crps_summary_list) %>%
  mutate(Model = model_names,
         LogScore = logscore_vector) %>%
  select(Model, Min., `1st Qu.`, Median, Mean, `3rd Qu.`, Max., LogScore)

crps_plot <- ggtexttable(crps_table, rows = NULL)
plot(crps_plot)

# make nicer table
model_eval <- crps_table
model_eval$CRPS_Mean <- model_eval$Mean
model_eval$CRPS_Min <- model_eval$Min.
model_eval$CRPS_Median <- model_eval$Median
model_eval$CRPS_Max <- model_eval$Max.
model_eval <- model_eval[,c("Model","CRPS_Min","CRPS_Median","CRPS_Mean","CRPS_Max","LogScore")]

model_eval_table <- ggtexttable(model_eval, rows = NULL)
plot(model_eval_table)

##############################################################################################
# summary statistics (WAIC, DIC, Marginal log-Likelihood)
##############################################################################################

# in-model diagnostics: not needed

# table of predictions and summaries
summary_table <- results_df %>%
  select(Model, predicted_abundance, WAIC, DIC, MLik) %>%
  distinct() %>%
  arrange(WAIC)

table_plot <- ggtexttable(summary_table, rows = NULL)
plot(table_plot)

# WAIC plot
waic_df <- results_df %>%
  filter(Model != "N") %>%
  select(Model, WAIC) %>%
  distinct()

waic_plot <- ggplot(waic_df, aes(x = Model, y = WAIC)) +
  geom_point(size = 3, color = "black") +
  labs(
    x = "Model",
    y = "WAIC"
  ) +
  theme_minimal()

waic_plot

# DIC plot
dic_df <- results_df %>%
  filter(Model != "N") %>%
  select(Model, DIC) %>%
  distinct()

dic_plot <- ggplot(dic_df, aes(x = Model, y = DIC)) +
  geom_point(size = 3, color = "black") +
  labs(
    x = "Model",
    y = "DIC"
  ) +
  theme_minimal()

dic_plot

# MLik plot
mlik_df <- results_df %>%
  filter(Model != "N") %>%
  select(Model, MLik) %>%
  distinct()

mlik_plot <- ggplot(mlik_df, aes(x = Model, y = MLik)) +
  geom_point(size = 3, color = "black") +
  labs(
    x = "Model",
    y = "MLik"
  ) +
  theme_minimal()

mlik_plot

# model selection ..........................................

# could use R2 predicted for predicting at Cape Crozier 2019

# select final model and then run partial predictions
# partial predictions

# vectors of full range of covariates
percent_guano <- data.frame(percent_guano = seq(-0.18382, 5.76802, length.out = 1000))
# other covariates depend on final model

# compute response curves

##############################################################################################
# Predict for 2019
##############################################################################################

# make 2019 newdata (2019 guano area but same terrain)
# add covariates at 2m res (aggregated below)

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
filtered_vals_2019 <- vals[!is.na(vals_2019)]
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

# match model component names
name_map <- c(
  percentguano = "2019_CrozierGuano_2m",
  slope        = "Cape_Crozier_slope",
  northness    = "northness",
  eastness     = "eastness",
  roughness    = "Cape_Crozier_Roughness",
  TRI          = "Cape_Crozier_TRI"
)

# extract 2019 guano values in each observation grid

# ensure count crs is the same as stack CRS
counts_sf <- sf::st_as_sf(counts_df)
counts_sf <- sf::st_transform(counts_sf, sf::st_crs(cov_stack2))

# convert to SpatVector and extract all covariates
coords <- terra::vect(counts_sf)     # SpatVector of points
ext <- terra::extract(cov_stack2, coords, ID = FALSE)

covariate_vals <- ext[, unname(name_map), drop = FALSE]
names(covariate_vals) <- names(name_map)

# don't need newdata section?
# build 2019 newdata for prediction
newdata2019 <- tibble(
  geometry = counts_df$geometry,
  area     = counts_df$area
) %>%
  bind_cols(covariate_vals) %>%
  select(geometry, area, all_of(names(name_map)))

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
