# Python script for turning penguin points into XY data (EPSG:3031)
# Copyright (C) 2025 Alexandra Strang and Dean Anderson 
# NOT DONE JUST TESTING

class Params:
    def __init__(self):
        ## LOOK AT ROYDS FIRST
        
        ## DATA FOLDERS
        self.tempinputDataPath =
        
        ## DEM FILE AND POINT SHAPEFILES
        self.RoydsDEMFName = os.path.join(self.tempinputDataPath, 'Cape_Royds_clipped.tif')
        self.RoydsPointsFName = os.path.join(self.tempinputDataPath, '2020_UAV_points',
            'royd_masked_labels_cleaned_coords_added_2020-12-01',
            'royd_masked_labels_cleaned_coords_added_2020-12-01.shp')

        ## OUTPUT DATA PATHS AND FILENAMES
        self.tempoutputDataPath = 
        self.RoydsPointsCSV = os.path.join(self.tempoutputDataPath, 'Royds_Points_3031.csv')

        ## MAKE NEW RESULTS DIRECTORY IF DOESN'T EXIST
        if not os.path.isdir(self.tempoutputDataPath):
            os.makedirs(self.tempoutputDataPath)
        

def makePenguinPointsXY(self):
        """
        ## READ IN SHAPEFILE, ADD XY AND EXPORT AS CSV
        """
          for col in self.params.colonies:    
            ## GET DATA NAMES AND PATHS
            demFName = getattr(self.params, '{}DEMFName'.format(col))
            ptShpFName = getattr(self.params, '{}PointsFName'.format(col))
            counts2mFName = getattr(self.params, '{}Penguin2m'.format(col))
            print('shp name:', ptShpFName)

            # OPEN DEM TO GET GEOREFERENCING
            dem_ds = gdal.Open(demFName)
            dem_band = dem_ds.GetRasterBand(1)
            dem_array = dem_band.ReadAsArray()
            dem_nodata = dem_band.GetNoDataValue()
            gt = dem_ds.GetGeoTransform()
            proj = dem_ds.GetProjection()
            rows, cols = dem_array.shape
  
            ## MAKE EMPTY ARRAY TO POPULATE
            XY_array = np.zeros((rows, cols), dtype=np.uint16)

            ## OPEN POINT SHAPEFILE
            shp_ds = ogr.Open(ptShpFName)
            layer = shp_ds.GetLayer()
            ## TRANSFORM CRS TO EPSG:3031
            source_srs = layer.GetSpatialRef()
            target_srs = osr.SpatialReference()
            target_srs.ImportFromEPSG(3031)  # Your DEM is in EPSG:3031
            coord_transform = osr.CoordinateTransformation(source_srs, target_srs)

            ## GET X Y DATA OF PENGUIN LOCATIONS
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
            addXYPoints(x_coords, y_coords, gt, XY_array)

# NOT DONE
