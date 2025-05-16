# python script for reprojecting colony AOI shapefiles to EPSG 3031
# creator: Alexandra Strang
# created: 2025

import fiona
from fiona.transform import transform_geom
import os

# find directories and paths
INPUT_DIR = r'D:\Terrain\Edited_RossSea_AOI'
OUTPUT_DIR = r'D:\Terrain\Reprojected_3031_AOI'
target_crs = "EPSG:3031"  # Antarctic Polar Stereographic

os.makedirs(OUTPUT_DIR, exist_ok=True)  # create output folder if it doesn't exist

# loop through all shapefiles
for filename in os.listdir(INPUT_DIR):
    if filename.endswith(".shp"):
        input_path = os.path.join(INPUT_DIR, filename)
        output_path = os.path.join(OUTPUT_DIR, filename)

        print(f"Reprojecting {filename}...")

        try:
            with fiona.open(input_path, "r") as src:
                src_crs = src.crs
                schema = src.schema

                with fiona.open(
                    output_path, "w",
                    driver=src.driver,
                    crs=target_crs,
                    schema=schema
                ) as dst:
                    for feature in src:
                        try:
                            reprojected_geom = transform_geom(
                                src_crs, target_crs, feature["geometry"]
                            )
                            dst.write({
                                "geometry": reprojected_geom,
                                "properties": feature["properties"]
                            })
                        except Exception as e:
                            print(f"Skipping feature in {filename} due to reprojection error: {e}")
        except Exception as e:
            print(f"Failed to process {filename}: {e}")