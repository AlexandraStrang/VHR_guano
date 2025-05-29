# Python script for estabilishing finescale relationships between nesting Adelie penguins, terrain and guano area
# creator: Alexandra Strang
# created: 2025

import os
import rasterio
import statistics
import subprocess
import numpy as np

# *pre-process rasters before to set minimum subcolony size?*

# directories and paths
cleaned_rasters = r"D:\Guano_rasters\Clean_guano_rasters"
terrain_files = r'D:\Terrain\Terrain_files'  # in seperate folders by colony name
INPUT_DIRS = [cleaned_rasters, terrain_files]
OUTPUT_DIRS = ["path/to/warped_colonies", "path/to/warped_dems"]

# check pixel resolution for each raster
cell_widths = []
cell_heights = []

for filename in os.listdir(cleaned_rasters):
    if filename.endswith(".tif"):
        filepath = os.path.join(cleaned_rasters, filename)
        with rasterio.open(filepath) as src:
            transform = src.transform
            cell_width = transform.a  # Pixel size in X direction (width)
            cell_height = -transform.e  # Pixel size in Y direction (height)

            cell_widths.append(cell_width)
            cell_heights.append(cell_height)

    # print(f"\n{filename}: Cell size (width, height): ({cell_width}, {cell_height})")

if cell_widths and cell_heights:
    min_width = min(cell_widths)
    max_width = max(cell_widths)
    min_height = min(cell_heights)
    max_height = max(cell_heights)
    median_width = statistics.median(cell_widths)
    median_height = statistics.median(cell_heights)

    print("\nOverall cell size range:")
    print(f"  Median: {median_width} and {median_height}")
    print(f"  Min: {min_width} and {min_height}")
    print(f"  Max: {max_width} and {max_height}")

# median cell size around 0.49 (gdalwarp to 0.50 for consistency?)
# DEM cell size 2 m (gdalwarp to 0.50 to match image rasters)
# need to ensure layers are aligned with same -te (target extent) and -tap (pixel alignment) when warping


def warp_rasters(input_path, output_path, target_res=(0.5, 0.5), target_extent=None, resampling='bilinear'):
    """
    function for warping rasters
    """
    cmd = [
        "gdalwarp",
        "-tr", str(target_res[0]), str(target_res[1]),
        "-r", resampling,
        "-overwrite"  # of the output directory
    ]

    # target extent part?

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


for INPUT_DIR, OUTPUT_DIR in zip(INPUT_DIRS, OUTPUT_DIRS):
    os.makedirs(OUTPUT_DIR, exist_ok=True)  # check if output dir exists
    for filename in os.listdir(INPUT_DIR):
        if filename.endswith(".tif"):
            input_path = os.path.join(INPUT_DIR, filename)
            output_path = os.path.join(OUTPUT_DIR, filename)
            warp_rasters(input_path, output_path, target_res=(0.5, 0.5))

# *add random sampler for selecting pixels below?*


def extract_values():
    """
    function for extracting number of nests, guano area, and terrain from pixels and storing in an array
    """
    
    pixel_values_array = np.array([])
