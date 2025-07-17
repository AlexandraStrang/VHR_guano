# r script for testing LGCP null models
# creator: Alexandra Strang
# created: 2025

library(sf)
library(fmesher)
library(ggplot2)
library(INLA)
library(ggpubr)

# include best meshes so far for Royds and Crozier


# EPSG 3031 is in metres but geom_sf will convert to degrees
# mesh previews = meshbuilder()

# read in points from xy csv

# Cape Royds 2020
Royds_xy <- read.csv("Royds_Points_3031.csv")
# Cape Crozier 2020
Crozier_xy <- read.csv("Crozier_Points_3031.csv")

# convert to sf object
sf_Royds <- st_as_sf(Royds_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Royds)

sf_Crozier <- st_as_sf(Crozier_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier)

# import guano area shapefile
sf_Royds_guano <- st_read("Point_process_GA_boundaries/Royds_2020_1_3031_guano.shp")
sf_Crozier_guano <- st_read("Point_process_GA_boundaries/Crozier_2020_1_3031_guano.shp")

# ensure shapefile has right crs code
sf_Royds_guano <- st_transform(sf_Royds_guano, crs = st_crs(sf_Royds))
sf_Crozier_guano <- st_transform(sf_Crozier_guano, crs = st_crs(sf_Crozier))

# use finer mesh
#1/15 max.edge = fine
Royds_max.edge <- diff(range(st_coordinates(sf_Royds)[,2]))/(3*5)
# 18 metres

# this is 1/15 study size using x range
# x range is longer than y here
Crozier_max.edge <- diff(range(st_coordinates(sf_Crozier)[,1]))/(3*5)
# ~160 m

# expand inner layer same amount (1/3)
Royds_bound.outer = diff(range(st_coordinates(sf_Royds)[,2]))/3
# 90 metres

# expand inner layer same amount (1/5)
Crozier_bound.outer = diff(range(st_coordinates(sf_Crozier)[,1]))/5
# 500 metres

# try finer mesh
Royds_mesh3 <- fm_mesh_2d(boundary = sf_Royds_guano,
                          max.edge = c(1,2)*Royds_max.edge, # inner and outer max edge where outer layer has triangle density two times lower than inner
                          offset = c(Royds_max.edge, Royds_bound.outer),
                          cutoff = Royds_max.edge/5, # cutoff 1/5 of max.edge
                          crs = st_crs(sf_Royds))
# coarser on the edge

# try finer mesh
Crozier_mesh5 <- fm_mesh_2d(boundary = sf_Crozier_guano,
                            #loc = st_coordinates(sf_Crozier), # more inner triangles?
                            max.edge = c(1,2)*Crozier_max.edge, # inner and outer max edge where outer layer has triangle density 2 times lower than inner
                            offset = c(Crozier_max.edge, Crozier_bound.outer),
                            cutoff = 0.5, # reduce cut off?
                            crs = st_crs(sf_Crozier))
# takes ~2 minutes 

