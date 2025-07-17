# r script for testing mesh method for LGCP model
# creator: Alexandra Strang
# created: 2025

library(sf)
library(fmesher)
library(ggplot2)
library(INLA)
library(ggpubr)

# EPSG 3031 is in metres but geom_sf will convert to degrees
# mesh previews = meshbuilder()

# Cape Royds 2020
# read in points from xy csv
Royds_xy <- read.csv("Royds_Points_3031.csv")

# convert to sf object
sf_Royds <- st_as_sf(Royds_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Royds)

# plot points
ggplot(sf_Royds) +
  geom_sf(color = "black", fill = "purple", shape = 21, size = 1.7) +
  labs(
    x = "Easting ()",
    y = "Northing ()"
  ) +
  theme_minimal()
# geom_sf converts to degrees

# first try max.edge of 1/3 of the size of the study area
max.edge <- diff(range(st_coordinates(sf_Royds)[,2]))/3
# 90 metres

# testing coarse mesh that uses location of points
Royds_mesh1 <- fm_mesh_2d(loc=st_coordinates(sf_Royds),
                   max.edge = max.edge,
                   crs = st_crs(sf_Royds))

# plot mesh with points
ggplot() +
  geom_fm(data = Royds_mesh1) +
  geom_sf(data=sf_Royds, color = 'purple', size=1.7, alpha=0.5) +
  labs(
    x = "Easting ()",
    y = "Northing ()"
  ) +
  theme_minimal()
# geom_sf converts to degrees
# use boundary instead

# import guano area shapefile
sf_Royds_guano <- st_read("Point_process_GA_boundaries/Royds_2020_1_3031_guano.shp")
# ensure shapefile has right crs code
sf_Royds_guano <- st_transform(sf_Royds_guano, crs = st_crs(sf_Royds))

# try finer mesh
#1/15 max.edge = fine
max.edge <- diff(range(st_coordinates(sf_Royds)[,2]))/(3*5)
# 18 metres

# use boundary of guano area polygon instead
Royds_mesh2 <- mesh2 <- fm_mesh_2d(boundary = sf_Royds_guano,
                                   max.edge = max.edge, 
                                   crs = st_crs(sf_Royds))

# just plot GA boundary mesh
mesh2_plot_Royds <- ggplot() +
  geom_fm(data = Royds_mesh2) +
  geom_sf(data = sf_Royds_guano, fill = NA, color = "blue", linetype = "dashed") +
  labs(
    x = "Easting ()",
    y = "Northing ()",
  ) +
  theme_minimal()
# geom_sf converts to degrees

mesh2_plot_Royds

# expand inner layer same amount as 1/3
bound.outer = diff(range(st_coordinates(sf_Royds)[,2]))/3
# 90 metres

Royds_mesh3 <- fm_mesh_2d(boundary = sf_Royds_guano,
                    max.edge = c(1,2)*max.edge, # inner and outer max edge where outer layer has triangle density two times lower than inner
                    offset = c(max.edge, bound.outer),
                    cutoff = max.edge/5, # cutoff 1/5 of max.edge
                    crs = st_crs(sf_Royds))

# just plot GA boundary mesh
mesh3_plot_Royds <- ggplot() +
  geom_fm(data = Royds_mesh3) +
  geom_sf(data = sf_Royds_guano, fill = NA, color = "blue", linetype = "dashed") +
  labs(
    x = "Easting ()",
    y = "Northing ()",
  ) +
  theme_minimal()
# geom_sf converts to degrees

mesh3_plot_Royds
# inner does not look finer than outer

# plot GA boundary mesh with points
mesh3.5_plot_Royds <- ggplot() + 
  geom_fm(data = Royds_mesh3) + 
  geom_sf(data = sf_Royds_guano, fill = NA, color = "blue", linetype = "dashed") + 
  geom_sf(data = sf_Royds, color = "purple", size = 1.7, alpha = 0.5) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
    ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh3.5_plot_Royds
# very coarse mesh beyond GA boundary
# way more points than GA boundary

# plot together
Together <- plot(ggarrange(mesh2_plot_Royds, 
                           mesh3_plot_Royds, 
                           mesh3.5_plot_Royds,
                           ncol = 3, nrow = 1, labels=c("a","b","c")))
#annotate_figure(Together, left = "Northing", bottom = "Easting")


# Cape Crozier 2020
# read in points from xy csv
Crozier_xy <- read.csv("Crozier_Points_3031.csv")

# convert to sf object
sf_Crozier <- st_as_sf(Crozier_xy, coords = c("x", "y"), crs = 3031)
st_crs(sf_Crozier)

# plot points
ggplot(sf_Crozier) +
  geom_sf(color = "black", fill = "purple", shape = 21, size = 1.7) +
  labs(
    x = "Easting ()",
    y = "Northing ()"
  ) +
  theme_minimal()
# takes ~60 seconds to show up
# geom_sf converts to degrees

# first try max.edge of 1/3 of the size of the study area
max.edge <- diff(range(st_coordinates(sf_Crozier)[,1]))/3
# ~800 m

# testing coarse mesh that uses location of points
Crozier_mesh1 <- fm_mesh_2d(loc=st_coordinates(sf_Crozier),
                    max.edge = max.edge,
                    crs = st_crs(sf_Crozier))
# takes ~3 minutes

# plot mesh with points
ggplot() +
  geom_fm(data = Crozier_mesh1) +
  geom_sf(data=sf_Crozier, color = 'purple', size=1.7, alpha=0.5) +
  labs(
    x = "Easting ()",
    y = "Northing ()"
  ) +
  theme_minimal()
# takes ~9 minutes to show up
# geom_sf converts to degrees

# import guano area shapefile
sf_Crozier_guano <- st_read("Point_process_GA_boundaries/Crozier_2020_1_3031_guano.shp")
# ensure shapefile has right crs code
sf_Crozier_guano <- st_transform(sf_Crozier_guano, crs = st_crs(sf_Crozier))

# use boundary of guano area polygon instead
Crozier_mesh2 <- fm_mesh_2d(boundary = sf_Crozier_guano,
                    max.edge = max.edge,
                    crs = st_crs(sf_Crozier))
# takes ~2 minutes

# just plot GA boundary mesh
mesh2_plot_Crozier <- ggplot() +
  geom_fm(data = Crozier_mesh2) +
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "blue", linetype = "dashed") +
  labs(
    x = "Easting ()",
    y = "Northing ()",
  ) +
  theme_minimal()
# takes ~40 seconds
# geom_sf converts to degrees

mesh2_plot_Crozier
# takes ~30 seconds

# plot GA boundary mesh with points
ggplot() +
  geom_fm(data = Crozier_mesh2) +
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "blue", linetype = "dashed") +
  geom_sf(data = sf_Crozier, color = "purple", size = 1.7, alpha = 0.5) +
  labs(
    x = "Easting ()",
    y = "Northing ()",
  ) +
  theme_minimal()
# takes ~3 minutes to show up
# geom_sf converts to degrees

# expand inner layer same amount as max.edge
bound.outer = diff(range(st_coordinates(sf_Crozier)[,1]))/3

Crozier_mesh3 <- fm_mesh_2d(boundary = sf_Crozier_guano,
                          max.edge = c(1,2)*max.edge, # inner and outer max edge where outer layer has triangle density two times lower than inner
                          offset = c(max.edge, bound.outer), # offset wrong?
                          cutoff = max.edge/5, # cutoff 1/5 of max.edge
                          crs = st_crs(sf_Crozier))
# takes ~20 seconds

# just plot GA boundary mesh
mesh3_plot_Crozier <- ggplot() +
  geom_fm(data = Crozier_mesh3) +
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "blue", linetype = "dashed") +
  labs(
    x = "Easting ()",
    y = "Northing ()",
  ) +
  theme_minimal()
# geom_sf converts to degrees

mesh3_plot_Crozier
# triangle density goes down
# too coarse

# plot GA boundary mesh with points
mesh3.5_plot_Crozier <- ggplot() + 
  geom_fm(data = Crozier_mesh3) + 
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "blue", linetype = "dashed") + 
  geom_sf(data = sf_Crozier, color = "purple", size = 1.7, alpha = 0.5) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
  ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh3.5_plot_Crozier
# takes ~80 seconds to show up
# very coarse mesh, would encompass points

# plot together
Together <- plot(ggarrange(mesh2_plot_Crozier, 
                           mesh3_plot_Crozier, 
                           mesh3.5_plot_Crozier,
                           ncol = 3, nrow = 1, labels=c("a","b","c")))
#annotate_figure(Together, left = "Northing", bottom = "Easting")

# try a finer mesh

# change max.edge
# this is 1/15 study size using x range
# x range is longer than y here
max.edge <- diff(range(st_coordinates(sf_Crozier)[,1]))/(3*5)
# ~160 m

# expand inner layer same amount (1/5)
bound.outer = diff(range(st_coordinates(sf_Crozier)[,1]))/5
# 500 metres

Crozier_mesh4 <- fm_mesh_2d(boundary = sf_Crozier_guano,
                            max.edge = c(1,2)*max.edge, # inner and outer max edge where outer layer has triangle density much lower than inner
                            offset = c(max.edge, bound.outer),
                            cutoff = max.edge/5, # cutoff 1/5 of max.edge
                            crs = st_crs(sf_Crozier))
# takes ~20 seconds

# just plot GA boundary mesh
mesh4_plot_Crozier <- ggplot() +
  geom_fm(data = Crozier_mesh4) +
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "blue", linetype = "dashed") +
  labs(
    x = "Easting ()",
    y = "Northing ()",
  ) +
  theme_minimal()
# geom_sf converts to degrees

mesh4_plot_Crozier
# inner mesh too coarse?

# plot GA boundary mesh with points
mesh4.5_plot_Crozier <- ggplot() + 
  geom_fm(data = Crozier_mesh4) + 
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "blue", linetype = "dashed") + 
  geom_sf(data = sf_Crozier, color = "purple", size = 1.7, alpha = 0.5) + 
  labs( 
    x = "Easting", 
    y = "Northing", 
  ) + 
  theme_minimal()
# geom_sf converts to degrees

mesh4.5_plot_Crozier

# try finer mesh
Crozier_mesh5 <- fm_mesh_2d(boundary = sf_Crozier_guano,
                            #loc = st_coordinates(sf_Crozier), # more inner triangles?
                            max.edge = c(1,2)*max.edge, # inner and outer max edge where outer layer has triangle density 2 times lower than inner
                            offset = c(max.edge, bound.outer),
                            cutoff = 0.5, # reduce cut off?
                            crs = st_crs(sf_Crozier))
# takes ~2 minutes 

# just plot GA boundary mesh
mesh5_plot_Crozier <- ggplot() +
  geom_fm(data = Crozier_mesh5) +
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "blue", linetype = "dashed") +
  labs(
    x = "Easting ()",
    y = "Northing ()",
  ) +
  theme_minimal()
# takes over 2 minutes
# geom_sf converts to degrees

mesh5_plot_Crozier
# takes over 10 minutes to show up
# inner mesh too coarse?

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

# look at finer max.edge and boundary of guano area polygon 
Crozier_mesh6 <- fm_mesh_2d(boundary = sf_Crozier_guano,
                            max.edge = max.edge,
                            crs = st_crs(sf_Crozier))

# just plot finer GA boundary mesh
mesh6_plot_Crozier <- ggplot() +
  geom_fm(data = Crozier_mesh6) +
  geom_sf(data = sf_Crozier_guano, fill = NA, color = "blue", linetype = "dashed") +
  labs(
    x = "Easting ()",
    y = "Northing ()",
  ) +
  theme_minimal()
# geom_sf converts to degrees

mesh6_plot_Crozier

# plot together
Together <- plot(ggarrange(mesh6_plot_Crozier,
                           mesh5_plot_Crozier, 
                           mesh5.5_plot_Crozier,
                           ncol = 3, nrow = 1, labels=c("a","b","c")))
#annotate_figure(Together, left = "Northing", bottom = "Easting")

