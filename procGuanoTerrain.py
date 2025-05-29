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
import pandas as pd
from osgeo import gdal, ogr, osr
from osgeo import gdal_array
from datetime import datetime
from glob import glob
import resource

osr.UseExceptions()
gdal.UseExceptions()


class Params:
    def __init__(self):
        ## DATA FOLDERS
        self.inputDataPath = os.path.join(os.getenv('ADELIEPROJDIR', default = '.'), 
            'Alexandra_Data')
        self.guano_folder = os.path.join(self.inputDataPath, 'clean_guano_rasters')
        self.terrain_folder = os.path.join(self.inputDataPath, 'colony_terrain_rasters')

        ## GUANO, DEM DATA FILES AND POINT SHAPEFILES
        self.CrozierGuanoFName = os.path.join(self.guano_folder, 'cleaned_Crozier_2020_1_3031.tif')
        self.RoydsGuanoFName = os.path.join(self.guano_folder, 'cleaned_Royds_2020_1_3031.tif')
        self.CrozierDEMFName = os.path.join(self.terrainPath, 'Cape_Crozier',
            'Cape_Crozier_clipped.tif')
        self.RoydsDEMFName = os.path.join(self.terrainPath, 'Cape_Royds',
            'Cape_Royds_clipped.tif')
        self.CrozierPointsFName = os.path.join(self.inputDataPath, '2020_UAV_points',
            'croz_masked_labels_cleaned_coords_added_2020-11-29',
            'croz_masked_labels_cleaned_coords_added_2020-11-29.shp')
        self.RoydsPointsFName = os.path.join(self.inputDataPath, '2020_UAV_points',
            'royd_masked_labels_cleaned_coords_added_2020-12-01',
            'royd_masked_labels_cleaned_coords_added_2020-12-01.shp')

        ## OUTPUT DATA PATHS AND FILENAMES
        self.outputDataPath = os.path.join(os.getenv('ADELIEPROJDIR', default = '.'), 
            'Alexandra_Data', 'Results_GuanoTerrain', 'Rasters_2m')
        self.CrozierGuano2m = os.path.join(self.outputDataPath, 'CrozierGuano_2m.tif')
        self.RoydsGuano2m = os.path.join(self.outputDataPath, 'RoydsGuano_2m.tif')
        self.CrozierPenguin2m = os.path.join(self.outputDataPath, 'CrozierPenguinCounts_2m.tif')
        self.RoydsPenguin2m = os.path.join(self.outputDataPath, 'RoydsPenguinCounts_2m.tif')

        # ## MAKE NEW RESULTS DIRECTORY IF DOESN'T EXIST
        if not os.path.isdir(self.outputDataPath):
            os.makedirs(self.outputDataPath)

        self.colonies = ['Crozier', 'Royds']


class DataProcessor:
    def __init__(self, params):
        """
        CLASS OBJ TO PROCESS DATA
        """
        self.params = params

        #############################
        ## RUN FUNCTIONS
        self.WarpRasters()
        self.makePenguinCountRaster()

        ## END RUNNING 
        #############################

    def WarpRasters(self):
        """
        ## READ IN AND REPROJECT RASTERS TO 2-M
        """
        for col in self.params.colonies:
            ## GET ATTRIBUTE TO RIGHT COLONY
            guanoFName = getattr(self.params, '{}GuanoFName'.format(col))
            guanoOut2m = getattr(self.params, {}2mData.format(col)
            print('guanoFName', guanoFName, guanoOut2m)
            ## READ IN GUANO RASTER AND RECLASS
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

            ## READ IN COLONY DEM
            demFName = getattr(self.params, '{}DEMFName'.format(col))
            print('demFName', demFName)

            ## REPROJECT RASTERS TO 2 M
            gdal.Warp(destNameOrDestDS=guanoOut2m,
                    srcDSOrSrcDSTab=mem_ds,
                    format='GTiff',
                    xRes=2.0,
                    yRes=2.0,
                    resampleAlg='average',
                    dstSRS='EPSG:3031',
                    outputBounds=dem_ds.GetGeoTransform(),  # use bounds of DEM (we refine this below)
                    dstNodata=-np.nan,
                    outputType=gdal.GDT_Float32,
                    creationOptions=["TILED=YES", "COMPRESS=LZW"],
                    warpOptions=["INIT_DEST=NO_DATA"],
                    multithread=True,
                    options=["-tap"],  # tap: align pixels to target grid
                    targetAlignedPixels=True,
                    outputBoundsSRS='EPSG:3031',
                    srcNodata=-9999)

    def makePenguinCountRaster(self):
        """
        ## READ IN SHAPEFILE, AND GET COUNTS AT 2-M RESOLUTION
        """
        for col in self.params.colonies:    
            ## GET DATA NAMES AND PATHS
            demFName = getattr(self.params, '{}DEMFName'.format(col))
            print('demFName', demFName)
            ptShpFName = getattr(self.params, '{}PointsFName'.format(col))
            counts2mFName = getattr(self.params, '{}Penguin2m'.format(col))
            print('demFName:', demFName, 'shp name:', ptShpFName, 'count name:', counts2mFName)

            # Open DEM to get georeferencing and size
            dem_ds = gdal.Open(demFName)
            geotransform = dem_ds.GetGeoTransform()
            projection = dem_ds.GetProjection()
#            x_min = geotransform[0]
#            y_max = geotransform[3]
            x_res = dem_ds.RasterXSize
            y_res = dem_ds.RasterYSize

            # Create in-memory target raster
            nodata_value = 65535
            mem_driver = gdal.GetDriverByName("MEM")
            target_ds = mem_driver.Create("", x_res, y_res, 1, gdal.GDT_UInt16)
            target_ds.SetGeoTransform(geotransform)
            target_ds.SetProjection(projection)
            band = target_ds.GetRasterBand(1)
            band.Fill(nodata_value)
            band.SetNoDataValue(nodata_value)

            # Rasterize using point count per pixel
            shp_ds = ogr.Open(ptShpFName)
            layer = shp_ds.GetLayer()
            gdal.RasterizeLayer(target_ds, [1],
                                layer,
                                burn_values=[1],
                                options=["ALL_TOUCHED=TRUE"],  
                                mergeAlg=gdal.GRA_Add)  # Count all overlapping points
            # Save to disk
            gtiff_driver = gdal.GetDriverByName("GTiff")
            gtiff_driver.CreateCopy(counts2mFName, target_ds)


def main():
    params = Params()
    processor = DataProcessor(params)

    maxMem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print('Max Mem Usage in KB', maxMem)

if __name__ == '__main__':
    main()
