# python script for batch reprojecting imagery
# creator: Alexandra Strang
# created: 2025

# to run this script the python interpreter must be OSGeo4w to ensure that gdal can be imported from osgeo

import os
from osgeo import gdal

# reproject images in epsg=3031 to match PGC DEM (REMA) projection
input_dir = r'D:\Ortho\vhr_images'
output_dir = r'D:\Ortho\vhr_images_3031'

# create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# loop through all tif files in the input directory
for file in os.listdir(input_dir):
    if file.lower().endswith('.tif'):
        input_path = os.path.join(input_dir, file)
        output_path = os.path.join(output_dir, file)

        # reproject imagery with handling error
        try:
            gdal.Warp(output_path, input_path, dstSRS='EPSG:3031', resampleAlg='bilinear')
            print(f"Reprojected: {file}")
        except RuntimeError as e:
             print(f"Error reprojecting {file}: {e}")
