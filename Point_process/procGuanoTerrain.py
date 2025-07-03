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
import resource
from pathlib import Path

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
        self.CrozierGuanoFName = os.path.join(self.guano_folder, 
            'cleaned_Crozier_2020_1_3031_guano.tif')
        self.RoydsGuanoFName = os.path.join(self.guano_folder, 
            'cleaned_Royds_2020_1_3031_guano.tif')
        self.CrozierDEMFName = os.path.join(self.terrain_folder, 'Cape_Crozier',
            'Cape_Crozier_clipped.tif')
        self.RoydsDEMFName = os.path.join(self.terrain_folder, 'Cape_Royds',
            'Cape_Royds_clipped.tif')
        self.CrozierPointsFName = os.path.join(self.inputDataPath, '2020_UAV_points',
            'croz_masked_labels_cleaned_coords_added_2020-11-29',
            'croz_masked_labels_cleaned_coords_added_2020-11-29.shp')
        self.RoydsPointsFName = os.path.join(self.inputDataPath, '2020_UAV_points',
            'royd_masked_labels_cleaned_coords_added_2020-12-01',
            'royd_masked_labels_cleaned_coords_added_2020-12-01.shp')

        ## OUTPUT DATA PATHS AND FILENAMES
        self.outputDataPath = os.path.join(os.getenv('ADELIEPROJDIR', default = '.'), 
            'Dean_Data', 'Results_GuanoTerrain', 'Rasters_2m')
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
            guanoOut2m = getattr(self.params, '{}Guano2m'.format(col))
            print('guanoFName', guanoFName, 'Out', guanoOut2m)
            ## FIND IN COLONY DEM
            demFName = getattr(self.params, '{}DEMFName'.format(col))
            print('demFName', demFName)

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

            ## OPEN DEM TO GET EXTENT AND BOUNDS
            dem_ds = gdal.Open(demFName)
            dem_gt = dem_ds.GetGeoTransform()
            x_min = dem_gt[0]
            y_max = dem_gt[3]
            pixel_width = dem_gt[1]
            pixel_height = dem_gt[5]
            x_size = dem_ds.RasterXSize
            y_size = dem_ds.RasterYSize
            ## BOUNDING BOX
            x_max = x_min + (x_size * pixel_width)
            y_min = y_max + (y_size * pixel_height)

            ## REPROJECT RASTERS TO 2 M
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

    def makePenguinCountRaster(self):
        """
        ## READ IN SHAPEFILE, AND GET COUNTS AT 2-M RESOLUTION
        """
        for col in self.params.colonies:    
            ## GET DATA NAMES AND PATHS
            demFName = getattr(self.params, '{}DEMFName'.format(col))
            ptShpFName = getattr(self.params, '{}PointsFName'.format(col))
            counts2mFName = getattr(self.params, '{}Penguin2m'.format(col))
            print('shp name:', ptShpFName, 'count name:', counts2mFName)

            # OPEN DEM TO GET GEOREFERENCING
            dem_ds = gdal.Open(demFName)
            dem_band = dem_ds.GetRasterBand(1)
            dem_array = dem_band.ReadAsArray()
            dem_nodata = dem_band.GetNoDataValue()
            gt = dem_ds.GetGeoTransform()
            proj = dem_ds.GetProjection()
            rows, cols = dem_array.shape

            ## MAKE EMPTY ARRAY TO POPULATE
            count_array = np.zeros((rows, cols), dtype=np.uint16)

            ## OPEN POINT SHAPEFILE
            shp_ds = ogr.Open(ptShpFName)
            layer = shp_ds.GetLayer()
            ## TRANSFORM CRS TO EPSG:3031
            source_srs = layer.GetSpatialRef()
            target_srs = osr.SpatialReference()
            target_srs.ImportFromEPSG(3031)  # Your DEM is in EPSG:3031
            coord_transform = osr.CoordinateTransformation(source_srs, target_srs)

            ## GET X Y DATA OF PENGUIN LOCATIONS.
            x_list, y_list = [], []
            for feat in layer:
                geom = feat.GetGeometryRef()
                geom.Transform(coord_transform)
                x, y = geom.GetX(), geom.GetY()
                x_list.append(x)
                y_list.append(y)
            x_coords = np.array(x_list)
            y_coords = np.array(y_list)
            ## POPULATE EMPTY 2D ARRAY
            addCountsToCells(x_coords, y_coords, gt, count_array)

            ## SHOULD MATCH NUMBER OF FEATURES IN SHAPEFILE
            print(col, "Penguin count total:", np.sum(count_array))
            print(col, "Nonzero pixels:", np.count_nonzero(count_array))

            ## SAVE GeoTIFF
            driver = gdal.GetDriverByName('GTiff')
            out_ds = driver.Create(counts2mFName, cols, rows, 1, gdal.GDT_UInt16)
            out_ds.SetGeoTransform(gt)
            out_ds.SetProjection(proj)
            out_band = out_ds.GetRasterBand(1)
            nodata_val = 65535
            if dem_nodata is not None:
                count_array[dem_array == dem_nodata] = nodata_val
            out_band.WriteArray(count_array)
            out_band.SetNoDataValue(nodata_val)

            out_ds.FlushCache()
            out_ds = None

@njit
def addCountsToCells(x_coords, y_coords, gt, out_array):
    """
    ## NUMBA FUNCTION TO POPULATE COUNTS INTO EMPTY ARRAY
    """
    for i in range(len(x_coords)):
        x = x_coords[i]
        y = y_coords[i]
        col = int((x - gt[0]) / gt[1])
#        row = int((y - gt[3]) / gt[5])  # gt[5] is negative for north-up
        row = int((gt[3] - y) / abs(gt[5]))  

        if 0 <= row < out_array.shape[0] and 0 <= col < out_array.shape[1]:
            out_array[row, col] += 1




def main():
    params = Params()
    processor = DataProcessor(params)

    maxMem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print('Max Mem Usage in KB', maxMem)

if __name__ == '__main__':
    main()
