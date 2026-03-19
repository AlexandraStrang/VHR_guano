# r script for computing partial predictions from inlabru models
# creator: Alexandra Strang
# created: 2026

# set working directory
setwd("C:/Users/ajs424/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data")

# load packages 
library(ggplot2)
library(INLA) # version 25.09.19
library(inlabru) # for bru() version 2.13.0
library(terra) # for rasters
library(patchwork) # for partial prediction plots

bru_options_set(control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))

# sessionInfo() - important to note r version and versions of inlabru, INLA, fmesher
# R version 4.5.1

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

##############################################################################################
# Covariates
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

# save raw stack to get original raster values later
raw_stack <- c(percent_guano_raster, slope_raster, northness_raster,
               eastness_raster, roughness_raster, TRI_raster)
names(raw_stack) <- c("percentguano", "slope", "northness", "eastness", "roughness", "TRI")

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

##############################################################################################
# Extract raster values
##############################################################################################

# extract ranges from standardised rasters for partial prediction sequences
# use the 1st to 99th percentiles to avoid extremes (cliffs/edges)
raster_ranges <- lapply(as.list(cov_stack), function(r) {
  vals <- global(r, fun = quantile, probs = c(0.01, 0.99), na.rm = TRUE)
  list(min = vals[1, 1], max = vals[1, 2])
})
names(raster_ranges) <- names(cov_stack)
raster_ranges

# save raw stats for back-transforming axes in partial prediction plots
raw_stats <- lapply(as.list(raw_stack), function(r) {
  list(
    mean = global(r, fun = "mean", na.rm = TRUE)[1, 1],
    sd   = global(r, fun = "sd",   na.rm = TRUE)[1, 1]
  )
})
names(raw_stats) <- names(raw_stack)

# prediction sequence length
n_seq <- 1000

# make prediction df for focal covariate
make_pred_df <- function(focal_var, raster_ranges, n = n_seq) {
  
  seq_vals <- seq(raster_ranges[[focal_var]]$min,
                  raster_ranges[[focal_var]]$max,
                  length.out = n)
  
  df <- data.frame(
    percentguano = rep(0, n),
    slope        = rep(0, n),
    northness    = rep(0, n),
    eastness     = rep(0, n),
    roughness    = rep(0, n),
    TRI          = rep(0, n)
  )
  
  df[[focal_var]] <- seq_vals
  df
}

# match names to covariate names
names(raster_ranges) <- c("percentguano", "slope", "northness", "eastness", "roughness", "TRI")

##############################################################################################
# Partial predictions
##############################################################################################

set.seed(28)

# compute partial predictions for each model
# (dropping field from formula)

# G model
G_pred_guano <- predict(
  G_model,
  newdata = make_pred_df("percentguano", raster_ranges),
  formula  = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano)),
  include  = character(0),
  n.samples = 1000
)
saveRDS(G_pred_guano, file = "Inlabru_outputs/G_pred_guano.rds")

# GS model
GS_pred_guano <- predict(
  GS_model,
  newdata = make_pred_df("percentguano", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + slope_eval(slope)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GS_pred_guano, file = "Inlabru_outputs/GS_pred_guano.rds")

GS_pred_slope <- predict(
  GS_model,
  newdata = make_pred_df("slope", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + slope_eval(slope)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GS_pred_slope, file = "Inlabru_outputs/GS_pred_slope.rds")

# GSNE model
GSNE_pred_guano <- predict(
  GSNE_model,
  newdata = make_pred_df("percentguano", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + slope_eval(slope) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GSNE_pred_guano, file = "Inlabru_outputs/GSNE_pred_guano.rds")

GSNE_pred_slope <- predict(
  GSNE_model,
  newdata = make_pred_df("slope", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + slope_eval(slope) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GSNE_pred_slope, file = "Inlabru_outputs/GSNE_pred_slope.rds")

GSNE_pred_northness <- predict(
  GSNE_model,
  newdata = make_pred_df("northness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + slope_eval(slope) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GSNE_pred_northness, file = "Inlabru_outputs/GSNE_pred_northness.rds")

GSNE_pred_eastness <- predict(
  GSNE_model,
  newdata = make_pred_df("eastness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + slope_eval(slope) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GSNE_pred_eastness, file = "Inlabru_outputs/GSNE_pred_eastness.rds")

# GR model
GR_pred_guano <- predict(
  GR_model,
  newdata = make_pred_df("percentguano", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + roughness_eval(roughness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GR_pred_guano, file = "Inlabru_outputs/GR_pred_guano.rds")

GR_pred_roughness <- predict(
  GR_model,
  newdata = make_pred_df("roughness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + roughness_eval(roughness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GR_pred_roughness, file = "Inlabru_outputs/GR_pred_roughness.rds")

# GRNE model
GRNE_pred_guano <- predict(
  GRNE_model,
  newdata = make_pred_df("percentguano", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + roughness_eval(roughness) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GRNE_pred_guano, file = "Inlabru_outputs/GRNE_pred_guano.rds")

GRNE_pred_roughness <- predict(
  GRNE_model,
  newdata = make_pred_df("roughness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + roughness_eval(roughness) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GRNE_pred_roughness, file = "Inlabru_outputs/GRNE_pred_roughness.rds")

GRNE_pred_northness <- predict(
  GRNE_model,
  newdata = make_pred_df("northness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + roughness_eval(roughness) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GRNE_pred_northness, file = "Inlabru_outputs/GRNE_pred_northness.rds")

GRNE_pred_eastness <- predict(
  GRNE_model,
  newdata = make_pred_df("eastness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + roughness_eval(roughness) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GRNE_pred_eastness, file = "Inlabru_outputs/GRNE_pred_eastness.rds")

# GT model
GT_pred_guano <- predict(
  GT_model,
  newdata = make_pred_df("percentguano", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + TRI_eval(TRI)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GT_pred_guano, file = "Inlabru_outputs/GT_pred_guano.rds")

GT_pred_TRI <- predict(
  GT_model,
  newdata = make_pred_df("TRI", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + TRI_eval(TRI)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GT_pred_TRI, file = "Inlabru_outputs/GT_pred_TRI.rds")

# GTNE model
GTNE_pred_guano <- predict(
  GTNE_model,
  newdata = make_pred_df("percentguano", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + TRI_eval(TRI) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GTNE_pred_guano, file = "Inlabru_outputs/GTNE_pred_guano.rds")

GTNE_pred_TRI <- predict(
  GTNE_model,
  newdata = make_pred_df("TRI", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + TRI_eval(TRI) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GTNE_pred_TRI, file = "Inlabru_outputs/GTNE_pred_TRI.rds")

GTNE_pred_northness <- predict(
  GTNE_model,
  newdata = make_pred_df("northness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + TRI_eval(TRI) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GTNE_pred_northness, file = "Inlabru_outputs/GTNE_pred_northness.rds")

GTNE_pred_eastness <- predict(
  GTNE_model,
  newdata = make_pred_df("eastness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + TRI_eval(TRI) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GTNE_pred_eastness, file = "Inlabru_outputs/GTNE_pred_eastness.rds")

# GNE model
GNE_pred_guano <- predict(
  GNE_model,
  newdata = make_pred_df("percentguano", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GNE_pred_guano, file = "Inlabru_outputs/GNE_pred_guano.rds")

GNE_pred_northness <- predict(
  GNE_model,
  newdata = make_pred_df("northness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GNE_pred_northness, file = "Inlabru_outputs/GNE_pred_northness.rds")

GNE_pred_eastness <- predict(
  GNE_model,
  newdata = make_pred_df("eastness", raster_ranges),
  formula   = ~ exp(Intercept_eval(rep(1, nrow(.data.))) + percentguano_eval(percentguano) + northness_eval(northness) + eastness_eval(eastness)),
  include   = character(0),
  n.samples = 1000
)
saveRDS(GNE_pred_eastness, file = "Inlabru_outputs/GNE_pred_eastness.rds")

##############################################################################################
# Plot response curves
##############################################################################################

# back transform values to original
back_transform <- function(x, var) {
  x * raw_stats[[var]]$sd + raw_stats[[var]]$mean
}

# partial prediction plot function
plot_partial <- function(pred_obj, var, xlabel, model_name) {
  
  df <- as.data.frame(pred_obj)
  df$x_raw <- back_transform(df[[var]], var)
  
  df$mean   <- df$mean   * 100
  df$q0.025 <- df$q0.025 * 100
  df$q0.975 <- df$q0.975 * 100
  
  ggplot(df, aes(x = x_raw)) +
    geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "#648FFF", alpha = 0.2) +
    geom_line(aes(y = mean), color = "#648FFF", linewidth = 1) +
    labs(x = xlabel, y = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "gray80"),
      axis.text = element_text(color = "#4a4a4a", size = 9),
      axis.title.x = element_text(color = "#4a4a4a", size = 10),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
}

# x axis labels (get original units?)
xlabs <- list(
  percentguano = "Guano cover (proportion)",
  slope        = "Slope (°)",
  northness    = "Northness",
  eastness     = "Eastness",
  roughness    = "Roughness (m)",
  TRI          = "TRI (m)"
)

# G model
G_plot <- (
  plot_partial(G_pred_guano, "percentguano", xlabs$percentguano, "G")
) +
  plot_annotation(
    title = "G Model",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  ) &
  labs(y = "Expected count per 100 m²")

G_plot

ggsave("Inlabru_outputs/G_partial_predictions.png", G_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

# GS model
GS_plot <- (
  plot_partial(GS_pred_guano, "percentguano", xlabs$percentguano, "GS") |
    plot_partial(GS_pred_slope, "slope",        xlabs$slope,        "GS")
) +
  plot_annotation(
    title = "GS",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  ) &
  labs(y = "Expected count per 100 m²")

GS_plot

ggsave("Inlabru_outputs/GS_partial_predictions.png", GS_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

# GSNE model
GSNE_plot <- (
  plot_partial(GSNE_pred_guano,     "percentguano", xlabs$percentguano, "GSNE") |
    plot_partial(GSNE_pred_slope,     "slope",        xlabs$slope,        "GSNE") |
    plot_partial(GSNE_pred_northness, "northness",    xlabs$northness,    "GSNE") |
    plot_partial(GSNE_pred_eastness,  "eastness",     xlabs$eastness,     "GSNE")
) +
  plot_annotation(
    title = "GSNE",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  ) &
  labs(y = "Expected count per 100 m²")

GSNE_plot

ggsave("Inlabru_outputs/GSNE_partial_predictions.png", GSNE_plot,
       width = 10, height = 5, units = "in",
       dpi = 600
)

# GR model
GR_plot <- (
  plot_partial(GR_pred_guano,     "percentguano", xlabs$percentguano, "GR") |
    plot_partial(GR_pred_roughness, "roughness",    xlabs$roughness,    "GR")
) +
  plot_annotation(
    title = "GR",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  ) &
  labs(y = "Expected count per 100 m²")

GR_plot

ggsave("Inlabru_outputs/GR_partial_predictions.png", GR_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

# GRNE model
GRNE_plot <- (
  plot_partial(GRNE_pred_guano,     "percentguano", xlabs$percentguano, "GRNE") |
    plot_partial(GRNE_pred_roughness, "roughness",    xlabs$roughness,    "GRNE") |
    plot_partial(GRNE_pred_northness, "northness",    xlabs$northness,    "GRNE") |
    plot_partial(GRNE_pred_eastness,  "eastness",     xlabs$eastness,     "GRNE")
) +
  plot_annotation(
    title = "GRNE",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  ) &
  labs(y = "Expected count per 100 m²")

GRNE_plot

ggsave("Inlabru_outputs/GRNE_partial_predictions.png", GRNE_plot,
       width = 10, height = 5, units = "in",
       dpi = 600
)

# GT model
GT_plot <- (
  plot_partial(GT_pred_guano, "percentguano", xlabs$percentguano, "GT") |
    plot_partial(GT_pred_TRI,   "TRI",          xlabs$TRI,          "GT")
) +
  plot_annotation(
    title = "GT",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  ) &
  labs(y = "Expected count per 100 m²")

GT_plot

ggsave("Inlabru_outputs/GT_partial_predictions.png", GT_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)

# GTNE model
GTNE_plot <- (
  plot_partial(GTNE_pred_guano,     "percentguano", xlabs$percentguano, "GTNE") |
    plot_partial(GTNE_pred_TRI,       "TRI",          xlabs$TRI,          "GTNE") |
    plot_partial(GTNE_pred_northness, "northness",    xlabs$northness,    "GTNE") |
    plot_partial(GTNE_pred_eastness,  "eastness",     xlabs$eastness,     "GTNE")
) +
  plot_annotation(
    title = "GTNE",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  ) &
  labs(y = "Expected count per 100 m²")

GTNE_plot

ggsave("Inlabru_outputs/GTNE_partial_predictions.png", GTNE_plot,
       width = 10, height = 5, units = "in",
       dpi = 600
)

# GNE model
GNE_plot <- (
  plot_partial(GNE_pred_guano,     "percentguano", xlabs$percentguano, "GNE") |
    plot_partial(GNE_pred_northness, "northness",    xlabs$northness,    "GNE") |
    plot_partial(GNE_pred_eastness,  "eastness",     xlabs$eastness,     "GNE")
) +
  plot_annotation(
    title = "GNE",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  ) &
  labs(y = "Expected count per 100 m²")

GNE_plot

ggsave("Inlabru_outputs/GNE_partial_predictions.png", GNE_plot,
       width = 8, height = 5, units = "in",
       dpi = 600
)
