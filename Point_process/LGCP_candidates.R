# r script for testing LGCP candidate models
# creator: Alexandra Strang
# created: 2025

# set working directory
setwd("/home/stranga/00_nesi_projects/landcare04225/Alexandra_Data/point_process_data/LGCP_nesi_data")

# load packages 
library(sf)
library(fmesher)
library(ggplot2)
library(INLA) # version 25.09.04
library(inlabru) # for lgcp() version 2.13.0
library(dplyr)
library(purrr)
library(tidyr)
library(terra) # for rasters

# sessionInfo() - important to note r version and versions of inlabru, INLA, fmesher
# R version 4.5.0 (2025-04-11 ucrt)

# EPSG 3031 is in metres but when plotting geom_sf will convert to degrees

# read in points from xy csv
# Cape Crozier 2020
Crozier_xy <- read.csv("Crozier_Points_2020_3031.csv")

# load in mesh parameters

# load in covariates

# add matern function 

# Candidate model codes (8):
# G - percent guano only
# GS - percent guano + slope
# GSA - percent guano + slope + aspect
# GR - percent guano + roughness
# GRA - percent guano + roughness + aspect
# GT - percent guano + TRI
# GTA - percent guano + TRI + aspect
# GA - percent guano + aspect

# N - Null model (spatial field only)

# Guano model
G_cmp
G_model

# Guano + slope model
GS_cmp
GS_model

# Guano + slope + aspect model
GSA_cmp
GSA_model

# Guano + roughness model

# Guano + roughness + aspect model

# Guano + TRI model

# Guano + TRI + aspect model

# Guano + aspect model
GA_cmp
GA_model

