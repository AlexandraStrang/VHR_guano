# Python script for turning penguin points into XY data (EPSG:3031)
# Copyright (C) 2025 Alexandra Strang and Dean Anderson
# NOT DONE JUST TESTING

import os
import pandas as pd
from osgeo import ogr, osr


class Params:
    def __init__(self):
        # LOOK AT ROYDS FIRST

        # DATA FOLDER
        self.DataPath = r"D:\Royds_points_test"

        # DEM FILE AND POINT SHAPEFILE
        # COULD HAVE IMPORTED FROM ALEXANDRA DATA
        self.RoydsDEMFName = os.path.join(self.DataPath, 'Cape_Royds_clipped.tif')
        self.RoydsPointsFName = os.path.join(self.DataPath, 'royd_masked_labels_cleaned_coords_added_2020-12-01.shp')

        # OUTPUT DATA FILENAME
        self.RoydsPointsCSV = os.path.join(self.DataPath, 'Royds_Points_3031.csv')

        self.colonies = ['Royds']


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
            target_srs.ImportFromEPSG(3031)  # Your DEM is in EPSG:3031
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
            csv_path = os.path.join(self.params.DataPath, f'{col}_Points_3031.csv')
            df.to_csv(csv_path, index=False)
            print(f'Saved: {csv_path}')


def main():
    params = Params()
    processor = DataProcessor(params)


if __name__ == '__main__':
    main()
