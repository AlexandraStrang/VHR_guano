# r script for testing mesh method for inlabru model
# creator: Alexandra Strang
# created: 2025

# set working directory
setwd("C:/Users/ajs424/OneDrive - University of Canterbury/ANTA - PhD/Data/Inlabru/Inlabru_data")

# load packages 
library(sf)
library(fmesher)
library(ggplot2)
library(INLA)

# sessionInfo() - important to note r version and versions of inlabru, INLA, fmesher
# R version 4.5.0 (2025-04-11 ucrt)

# EPSG 3031 is in metres but when plotting geom_sf will convert to degrees

##############################################################################################
# Load data
##############################################################################################

# read in points from xy csv
# Cape Crozier 2020
Crozier_xy <- read.csv("Crozier_UAV_points/Reprojected_3031/Crozier_2020_Points_3031.csv")

# convert to sf object and check CRS
sf_Crozier <- st_as_sf(Crozier_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier)

##############################################################################################
# Mesh
##############################################################################################

# import Crozier coastline boundary
sf_Crozier_boundary <- st_read("Crozier_boundary.shp")
# transform shapefile so it has the right crs code
sf_Crozier_boundary <- st_transform(sf_Crozier_boundary, crs = st_crs(sf_Crozier))
st_crs(sf_Crozier_boundary)

# add 100 m buffer around coastline boundary
buff_boundary <- st_buffer(sf_Crozier_boundary, dist = 100)
# get domain area
boundary_area <- st_area(buff_boundary)
print(boundary_area) # 4.17 km2

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

# subdivide mesh
# splits triangles into subtriangles
# mesh gets subdivided by 3 (max edge 30m)
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

# save outputs
saveRDS(mesh_sub,      file = "Inlabru_outputs/mesh_sub.rds")
saveRDS(buff_boundary, file = "Inlabru_outputs/buff_boundary.rds")

# save unsubdivded mesh as shapefile - Koerich et al.
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
# Note: now a smaller mesh to what was originally saved
