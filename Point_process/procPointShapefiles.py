# Python script for turning penguin points into XY data (EPSG:3031)
# Copyright (C) 2025 Alexandra Strang and Dean Anderson

import os
import pandas as pd
from osgeo import ogr, osr

# use PPP_conda_env in NeSI

class Params:
    def __init__(self):
        # DATA FOLDER
        self.inputDataPath = r"/nesi/project/landcare04225/Alexandra_Data"
        
        self.terrain_folder = os.path.join(self.inputDataPath, 'colony_terrain_rasters')

        # DEM FILE AND POINT SHAPEFILE
        self.RoydsDEMFName = os.path.join(self.terrain_folder, 'Cape_Royds',
            'Cape_Royds_clipped.tif')
        self.RoydsPointsFName = os.path.join(self.inputDataPath, '2020_UAV_points',
            'royd_masked_labels_cleaned_coords_added_2020-12-01',
            'royd_masked_labels_cleaned_coords_added_2020-12-01.shp')
        
        self.CrozierDEMFName = os.path.join(self.terrain_folder, 'Cape_Crozier',
            'Cape_Crozier_clipped.tif')
        self.CrozierPointsFName = os.path.join(self.inputDataPath, '2020_UAV_points',
            'Edit_croz_masked_labels_cleaned_coords_added_2020-11-29',
            'Edit_croz_masked_labels_cleaned_coords_added_2020-11-29.shp')
 
        ## OUTPUT DATA PATHS AND FILENAMES
        self.outputDataPath = os.path.join(os.getenv('ADELIEPROJDIR', default = '.'), 
            'Alexandra_Data', 'point_process_data')
 
        self.colonies = ['Royds', 'Crozier']


class DataProcessor:
    def __init__(self, params):
        """
        CLASS OBJ TO PROCESS DATA
        """
        self.params = params

        #############################
        # RUN FUNCTIONS
        self.makePenguinPointsXY()

        # END RUNNING
        #############################

    def makePenguinPointsXY(self):
        """
        ## READ IN SHAPEFILE, TRANSFORM TO 3031, EXTRACT XY AND EXPORT AS CSV
        """
        for col in self.params.colonies:
            # GET DATA NAMES AND PATHS
            demFName = getattr(self.params, f'{col}DEMFName')
            ptShpFName = getattr(self.params, f'{col}PointsFName')
            print(f'demFName: {demFName}, shp name: {ptShpFName}')

            # OPEN POINT SHAPEFILE
            shp_ds = ogr.Open(ptShpFName)
            layer = shp_ds.GetLayer()
            source_srs = layer.GetSpatialRef()

            # TRANSFORM CRS TO EPSG:3031
            target_srs = osr.SpatialReference()
            target_srs.ImportFromEPSG(3031)  # DEM is in EPSG:3031
            coord_transform = osr.CoordinateTransformation(source_srs, target_srs)

            # GET X Y DATA OF PENGUIN LOCATIONS
            x_list, y_list = [], []

            for feat in layer:
                geom = feat.GetGeometryRef()
                geom.Transform(coord_transform)
                x, y = geom.GetX(), geom.GetY()
                x_list.append(x)
                y_list.append(y)

            # CONVERT TO DATAFRAME
            df = pd.DataFrame({"x": x_list, "y": y_list})

            # WRITE CSV
            csv_path = os.path.join(self.params.outputDataPath, f'{col}_Points_2020_3031.csv')
            df.to_csv(csv_path, index=False)
            print(f'Saved: {csv_path}')


def main():
    params = Params()
    processor = DataProcessor(params)


if __name__ == '__main__':
    main()
