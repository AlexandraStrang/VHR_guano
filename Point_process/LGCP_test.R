# r script for testing LGCP null models
# creator: Alexandra Strang
# created: 2025

library(sf)
library(fmesher)
library(ggplot2)
library(INLA)
library(ggpubr)
library(inlabru) # for lgcp()

# sessionInfo() - important to note r version and versions of inlabru, INLA, fmesher
# R version 4.5.0 (2025-04-11 ucrt)

# EPSG 3031 is in metres but when plotting geom_sf will convert to degrees

# Construct file path
path <- file.path("nesi", "project", "landcare04225", "Alexandra_Data", "point_process_data")
print(path)

# set working dir
setwd(path)

# read in points from xy csv
# Cape Crozier 2020
Crozier_xy <- read.csv("Crozier_Points_2020_3031.csv")

# convert to sf object and check CRS
sf_Crozier <- st_as_sf(Crozier_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier)

# import guano area shapefile (not needed until modelling spatial covariates)
sf_Crozier_guano <- st_read("Crozier/Crozier_2020_GA/Crozier_2020_1_3031_guano.shp")
# ensure shapefile has right crs code
sf_Crozier_guano <- st_transform(sf_Crozier_guano, crs = st_crs(sf_Crozier))

# import Crozier coastline boundary
sf_Crozier_boundary <- st_read("Crozier/Crozier_boundary_coastline/Crozier_boundary.shp")
# transform shapefile so it has the right crs code
sf_Crozier_boundary <- st_transform(sf_Crozier_boundary, crs = st_crs(sf_Crozier))
st_crs(sf_Crozier_boundary)

# import Crozier UAV bounds shapefile
sf_Crozier_UAV_area <- st_read("Crozier/Crozier_2020_UAV_area/Crozier_20201129_3031_UAV_area.shp")
# ensure shapefile has right crs code
sf_Crozier_UAV_area <- st_transform(sf_Crozier_UAV_area, crs = st_crs(sf_Crozier))
st_crs(sf_Crozier_UAV_area)

# mesh parameters
# this is 1/15 study size using x range
# x range is longer than y here
Crozier_max.edge <- diff(range(st_coordinates(sf_Crozier)[,1]))/(3*5)
print("Printing Crozier max.edge:" Crozier_max.edge)
# ~150 metres

# expand outer layer out by 1/5
Crozier_bound.outer = diff(range(st_coordinates(sf_Crozier)[,1]))/5
print("Printing Crozier bound.outer:" Crozier_bound.outer)
# ~500 metres

# create finer Crozier mesh
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

# plot boundary and GA
mesh_plot_Crozier2 <- ggplot() + 
  geom_fm(data = Crozier_mesh) + 
  geom_sf(data = sf_Crozier_boundary, fill = NA, color = "blue", linetype = "dashed") + 
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "red", linetype = "dashed") + 
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

# start with prior for spatial range using half the study area
# [,1] for x range
0.5*(diff(range(st_coordinates(sf_Crozier)[,1])))
# 1255.526 metres

# the probability of the range exceeding that to be 0.5
# the prior for the variance explained by the spatial effect
# is set that the probability that the SD exceeds 1 is 0.5

# define the SPDE prior (matern)
matern <- inla.spde2.pcmatern(mesh = Crozier_mesh,
                            prior.range = c(1000, 0.5), # try 100 or 10m (might be like 5m?)
                            prior.sigma = c(1, 0.5))

# specify a model where for 2D models geometry is on the left of ~
# and an SPDE + Intercept(1) on the right
# define the domain of the LGCP and model components
# (spatial SPDE effect and Intercept)
null_cmp <- geometry ~
  Intercept(1) + # fixed effect (eventually add here the covariates of slope etc. using rasters)
  mySmooth(geometry, model = matern) # random effect

# formula (geometry is response = point locations)
# mysmooth spatial autocorrelation

# use lgcp() with 2D model components, the sf points and the sf boundary
null_model <- lgcp(null_cmp, # formula
                   data = sf_Crozier, # locations
                   samplers = sf_Crozier_UAV_area, # sample area of UAV bounds
                   domain = list(geometry = Crozier_mesh), # mesh
                   options = ) # no control fix
# took several hours to run (could indicate bad model)

summary(null_model)
# mean should be desnity per m 2

# predicting intensity
# predict the spatial intensity surface
lambda <- predict(
  null_model,
  fm_pixels(Crozier_mesh, format = "sf", mask = sf_Crozier_guano),
  ~ exp(mySmooth + Intercept) # to unlog, mySmooth is field
)

sum(lambda$mean) # check how big

# can plot log and and log lambda intensity

