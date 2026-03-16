#!/usr/bin/env python

########################################
########################################
# This file is part of Ph.D. work by Alexandra Strang on Adelie population dynamics in the Ross Sea
# Copyright (C) 2025 Alexandra Strang and Dean Anderson 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ADDITIONAL NOTES ABOUT THIS SCRIPT.
########################################
########################################

import os
import json
import numpy as np
from numba import njit
from osgeo import gdal, ogr, osr
from pathlib import Path

osr.UseExceptions()
gdal.UseExceptions()


class Params:
    def __init__(self):
        # DATA FOLDER
        self.inputDataPath = r"C:/Users/astra/OneDrive - University of Canterbury/ANTA - PhD\Data/Inlabru/Inlabru_data/"

        self.terrain_folder = os.path.join(self.inputDataPath, 'Crozier_terrain_mesh')

        # GUANO, DEM DATA FILES AND POINT SHAPEFILES
        self.CrozierGuanoFName = os.path.join(self.inputDataPath, 'cleaned_Crozier_2019_1_3031_guano.tif')  # changed to 2019

        self.CrozierDEMFName = os.path.join(self.terrain_folder, 'Cape_Crozier_clipped.tif')

        # OUTPUT DATA PATHS AND FILENAMES
        self.outputDataPath = self.inputDataPath
        self.CrozierGuano2m = os.path.join(self.outputDataPath, '2019_CrozierGuano_2m.tif')

        self.colonies = ['Crozier']  # using on 2019 guano area


class DataProcessor:
    def __init__(self, params):
        """
        CLASS OBJ TO PROCESS DATA
        """
        self.params = params

        #############################
        # RUN FUNCTIONS
        self.WarpRasters()

        # END RUNNING 
        #############################

    def WarpRasters(self):
        """
        ## READ IN AND REPROJECT RASTERS TO 2-M
        """
        for col in self.params.colonies:
            # GET ATTRIBUTE TO RIGHT COLONY
            guanoFName = getattr(self.params, '{}GuanoFName'.format(col))
            guanoOut2m = getattr(self.params, '{}Guano2m'.format(col))
            print('guanoFName', guanoFName, 'Out', guanoOut2m)
            # FIND IN COLONY DEM
            demFName = getattr(self.params, '{}DEMFName'.format(col))
            print('demFName', demFName)

            # READ IN GUANO RASTER AND RECLASS
            src_ds = gdal.Open(guanoFName)
            band = src_ds.GetRasterBand(1)
            data = band.ReadAsArray()
            data[data == 2] = 0

            # CREATE IN MEMORY RASTER WITH GEOREFERENCE
            driver = gdal.GetDriverByName('MEM')
            mem_ds = driver.Create('', src_ds.RasterXSize, src_ds.RasterYSize, 1, gdal.GDT_Byte)
            mem_ds.SetGeoTransform(src_ds.GetGeoTransform())
            mem_ds.SetProjection(src_ds.GetProjection())
            mem_ds.GetRasterBand(1).WriteArray(data)

            # OPEN DEM TO GET EXTENT AND BOUNDS
            dem_ds = gdal.Open(demFName)
            dem_gt = dem_ds.GetGeoTransform()
            x_min = dem_gt[0]
            y_max = dem_gt[3]
            pixel_width = dem_gt[1]
            pixel_height = dem_gt[5]
            x_size = dem_ds.RasterXSize
            y_size = dem_ds.RasterYSize
            # BOUNDING BOX
            x_max = x_min + (x_size * pixel_width)
            y_min = y_max + (y_size * pixel_height)

            # REPROJECT RASTERS TO 2 M
            gdal.Warp(destNameOrDestDS=guanoOut2m,
                    srcDSOrSrcDSTab=mem_ds,
                    format='GTiff',
                    xRes=2.0,
                    yRes=-2.0,
                    resampleAlg='average',
                    dstSRS='EPSG:3031',
                    outputBounds=(x_min, y_min, x_max, y_max),
                    dstNodata=-np.nan,
                    outputType=gdal.GDT_Float32,
                    creationOptions=["TILED=YES", "COMPRESS=LZW"],
                    warpOptions=["INIT_DEST=NO_DATA"],
                    multithread=True,
                    targetAlignedPixels=True,
                    outputBoundsSRS='EPSG:3031',
                    srcNodata=-9999)


def main():
    params = Params()
    processor = DataProcessor(params)


if __name__ == '__main__':
    main()
