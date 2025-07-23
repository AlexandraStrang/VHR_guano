# r script for testing LGCP null models
# creator: Alexandra Strang
# created: 2025

library(sf)
library(fmesher)
library(ggplot2)
library(INLA)
library(ggpubr)

# best mesh so far Crozier 2020
# use Crozier 2019 as testing 
# EPSG 3031 is in metres but geom_sf will convert to degrees
# mesh previews = meshbuilder()

# read in points from xy csv
# Cape Crozier 2020
Crozier_xy <- read.csv("Crozier_Points_3031.csv")

# convert to sf object and check CRS
sf_Crozier <- st_as_sf(Crozier_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier)

# import guano area shapefiles
sf_Crozier_guano <- st_read("Point_process_GA_boundaries/Crozier_2020_1_3031_guano.shp")

# ensure shapefile has right crs code
sf_Crozier_guano <- st_transform(sf_Crozier_guano, crs = st_crs(sf_Crozier))

# mesh parameters
# this is 1/15 study size using x range
# x range is longer than y here
Crozier_max.edge <- diff(range(st_coordinates(sf_Crozier)[,1]))/(3*5)
# 160 metres

# expand inner layer same amount (1/5)
Crozier_bound.outer = diff(range(st_coordinates(sf_Crozier)[,1]))/5
# 500 metres

# create Crozier mesh
Crozier_mesh5 <- fm_mesh_2d(boundary = sf_Crozier_guano,
                            max.edge = c(1,5)*Crozier_max.edge, # inner and outer max edge where outer layer has triangle density lower than inner
                            offset = c(Crozier_max.edge, Crozier_bound.outer/2),
                            cutoff = Crozier_max.edge/10, # reduced cut off?
                            crs = st_crs(sf_Crozier))

# plot GA boundary mesh with points
mesh5.5_plot_Crozier <- ggplot() + 
  geom_fm(data = Crozier_mesh5) + 
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "blue", linetype = "dashed") + 
  geom_sf(data = sf_Crozier, color = "purple", size = 1.7, alpha = 0.5) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
  ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh5.5_plot_Crozier

# create the spatial random field
