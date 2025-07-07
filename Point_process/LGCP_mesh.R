# r script for testing mesh method for LGCP model
# creator: Alexandra Strang
# created: 2025

library(sf)
library(fmesher)
library(ggplot2)

# read in points from xy csv
Royds_xy <- read.csv("Royds_Points_3031.csv")

# convert to sf object
sf_Royds <- st_as_sf(Royds_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Royds)

# plot points
ggplot(sf_Royds) +
  geom_sf(color = "black", fill = "purple", shape = 21, size = 1.7) +
  labs(
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()

# first try max.edge of 1/3 of the size of the study area
max.edge <- diff(range(st_coordinates(sf_Royds)[,1]))/(3*5)

# testing coarse mesh that uses location of points
mesh1 <- fm_mesh_2d(loc=st_coordinates(sf_Royds),
                   max.edge = max.edge,
                   crs = st_crs(sf_Royds))

# plot mesh with points
ggplot() +
  geom_fm(data = mesh1) +
  geom_sf(data=sf_Royds, color = 'purple', size=1.7, alpha=0.5) +
  labs(
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal()
# EPSG 3031 causes  plot distortion

# import guano area shapefile
sf_Royds_guano <- st_read("Point_process_GA_boundaries/Royds_2020_1_3031_guano.shp")

# use boundary of guano area polygon instead
mesh2 <- fm_mesh_2d(boundary = sf_Royds_guano,
                    max.edge = max.edge,
                    crs = st_crs(sf_Royds))

# just plot GA boundary mesh
ggplot() +
  geom_fm(data = mesh2) +
  geom_sf(data = sf_Royds_guano, fill = NA, color = "blue", linetype = "dashed") +
  labs(
    x = "Easting (m)",
    y = "Northing (m)",
  ) +
  theme_minimal()

# plot GA boundary mesh with points
ggplot() +
  geom_fm(data = mesh2) +
  geom_sf(data = sf_Royds_guano, fill = NA, color = "blue", linetype = "dashed") +
  geom_sf(data = sf_Royds, color = "purple", size = 1.7, alpha = 0.5) +
  labs(
    x = "Easting (m)",
    y = "Northing (m)",
  ) +
  theme_minimal()
# wayyy more points than mesh area (eek)
# have a look at Crozier
# create buffer too
