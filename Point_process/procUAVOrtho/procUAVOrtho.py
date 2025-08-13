# Python script for cleaning UAV orthomoasaic from Point Blue into polygon area
# Copyright (C) 2025 Alexandra Strang

import rioxarray
import geopandas
from geocube.api.core import make_geocube
from geocube.vector import vectorize
import os

# use conda env proc_env

# OLD 
# now split into three parts: resample, reclassify, vectorise


class Params:
    def __init__(self):
        # DATA FOLDER
        self.inputDataPath = r"D:\PhD_Chap1_part2\VHR_UAV_Tiles"

        # CROZIER 2020 UAV FILENAME
        self.Crozier2020UAVFName = os.path.join(self.inputDataPath,
            'Crozier_20201129_3031.tif')  # EPSG 3031 transformed version of UAV

        # OUTPUT DATA PATH
        self.outputDataPath = os.path.join(self.inputDataPath, 
            'Cleaned_crozier_UAV_raster')

        self.Crozier2020classifiedFName = os.path.join(
            self.outputDataPath, 'Classified_Temp_Crozier_2020.tif'
        )
        self.shapefile_path = os.path.join(
            self.outputDataPath, 'Crozier_20201129_UAV_Area.shp'
        )


params = Params()

# read in raster file
Crozier_UAV_raster = rioxarray.open_rasterio(params.Crozier2020UAVFName, chunks=True)

# raster info
print(Crozier_UAV_raster)
print(f"CRS of raster: {Crozier_UAV_raster.rio.crs}\n")

# # reclassify the raster
# # set values between 0 and 255 (exclusive of 0 and 255) to 1
# # set everything else to 0
# nodata_val = Crozier_UAV_raster.rio.nodata
# reclassified_raster = Crozier_UAV_raster.where(
#     (Crozier_UAV_raster > 0) & (Crozier_UAV_raster < 255), 1, 0
# )
# # nodata value is preserved if needed
# reclassified_raster = reclassified_raster.rio.write_nodata(0)

# print(reclassified_raster)

# # run reclassification and save as temp file
# reclassified_raster.rio.to_raster(
#     params.Crozier2020classifiedFName,
#     compress='LZW',
#     tiled=True
#     )

# open temp 4-band raster
classified_raster = rioxarray.open_rasterio(params.Crozier2020classifiedFName, chunks=True)

# select single band to vectorise (array has no name)
classified_vector_1_band = classified_raster.sel(band=1)

vectorized_gdf = vectorize(classified_vector_1_band)  # geodataframe

UAV_polygon = vectorized_gdf[vectorized_gdf.iloc[:, 0] == 1]
print(f"CRS of polygon: {UAV_polygon.crs}")

# save polygon
UAV_polygon.to_file(params.shapefile_path)
