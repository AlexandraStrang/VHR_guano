# r script for testing mesh method for LGCP model
# creator: Alexandra Strang
# created: 2025

# set working directory
setwd("/home/stranga/00_nesi_projects/landcare04225/Alexandra_Data/point_process_data/LGCP_nesi_data")

# load packages 
library(sf)
library(fmesher)
library(ggplot2)
library(INLA)

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
# this is 1/15 study size using x range
# x range is longer than y here
Crozier_max.edge <- diff(range(st_coordinates(sf_Crozier)[,1]))/(3*5)
print(Crozier_max.edge)
# ~150 metres

# expand outer layer
Crozier_bound.outer = diff(range(st_coordinates(sf_Crozier)[,1]))/5
print(Crozier_bound.outer)
# ~450 metres

# create Crozier mesh
Crozier_mesh <- fm_mesh_2d(boundary = buff_boundary, # use buffered coastline as boundary
                           max.edge = c(1,5)*Crozier_max.edge, # inner and outer max edge where outer layer has triangle density lower than inner
                           offset = c(Crozier_max.edge, Crozier_bound.outer),
                           cutoff = Crozier_max.edge/5,
                           crs = st_crs(sf_Crozier))
print(Crozier_mesh$n)

# plot boundary mesh
mesh_plot_Crozier <- ggplot() + 
  geom_fm(data = Crozier_mesh) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
  ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh_plot_Crozier

# save mesh as shapefile - Koerich et al.
# to crop terrain variables to mesh boundary
# save Cozier mesh (do once)
vertices <- Crozier_mesh$loc  # Vertex coordinates
triangles <- Crozier_mesh$graph$tv  # Triangle indices

# create polygons from triangles
polygon_list <- lapply(1:nrow(triangles), function(i) {
  # Get the vertex indices for the current triangle
  tri <- triangles[i, ]
  
  # create a matrix of coordinates for the triangle
  coords <- vertices[tri, c(1, 2)]  
  
  # close the polygon by repeating the first point
  coords <- rbind(coords, coords[1, ])
  
  # create an sf polygon
  st_polygon(list(coords))
})

# combine all polygons into an sf object
mesh_polygons <- st_sf(
  geometry = st_sfc(polygon_list),
  crs = st_crs(sf_Crozier) # Set CRS
)

# save polygon
st_write(mesh_polygons, "Crozier_mesh/Crozier_mesh.shp")
