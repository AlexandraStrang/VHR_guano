# python script for calculating terrain variables from REMA DEMs
# creator: Alexandra Strang
# created: 2025

# Use conda environmenet pgc_test
# REMA DEM projection in WGS 1984 Antartic Polar Stereographic (EPSG 3031)
# https://gdal.org/en/stable/programs/gdaldem.html
# Could adapt script to add values into CSV
# Mask by colony AOI for smaller tif files

import os
import typing
import numpy as np
from osgeo.gdal import DEMProcessing, UseExceptions, Dataset, Band
from dataclasses import dataclass
import fiona
import rasterio
from rasterio.mask import mask

# find directories and paths
DEM_FILES_DIR = r'D:\Terrain\REMA_DEM'
OUTPUT_DIR = r'D:\Terrain\Terrain_files'
AOI_DIR = r'D:\Terrain\Reprojected_3031_AOI'


@dataclass
class DemProcessDesc:
    kwargs: dict
    dataset: Dataset = None
    metrics: dict = None


class DemMetrics:
    """
    define dem proccess and metrics
    """
    def __init__(self, dem_path: os.PathLike, colony_name: str, output_base_dir: os.PathLike):
        self.original_dem_path = dem_path
        self.colony_name = colony_name
        self.out_dir = os.path.join(output_base_dir, colony_name)
        self.dem_processes = {'slope': DemProcessDesc({'alg': 'Horn'}),  # Use Horn's method (3x3 neighbourhood of elevation values), -9999 is no data (or flat?)
                              'aspect': DemProcessDesc({}),  # aspect might need correcting (off 180?) Azimuth aspect where 0 degress is north
                              'Roughness': DemProcessDesc({}),  # roughness takes the difference between the maximum and minimum elevation in a 3×3 window around each pixel (elevation range of local varaibility and more sensitive to extreme changes)
                              'TRI': DemProcessDesc({'alg': 'Riley'})}  # Terrain Ruggedness Index (TRI) using the method developed by Riley et al. (1999), quantifies how “rough” the terrain, measuring elevation differences between a central pixel and its 8 neighbors

    def _get_ds_metrics(self, ds: Dataset):
        band: Band = ds.GetRasterBand(1)
        array: np.ndarray = band.ReadAsArray()
        return {
            'min': float(array.min()),
            'max': float(array.max()),
            'mean': float(array.mean()),
            'std': float(array.std())
        }

    def _clip_dem_to_aoi(self):
        """
        clip the DEM to the matching AOI shapefile
        """
        aoi_path = os.path.join(AOI_DIR, f'{self.colony_name}.shp')
        if not os.path.exists(aoi_path):
            raise FileNotFoundError(f"AOI shapefile not found for {self.colony_name}: {aoi_path}")

        with fiona.open(aoi_path, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]

        with rasterio.open(self.original_dem_path) as src:
            try:
                out_image, out_transform = mask(src, shapes, crop=True)
            except ValueError as e:
                raise RuntimeError(f"Failed to mask {self.colony_name} — likely no overlap: {e}")

            out_meta = src.meta.copy()
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })

            os.makedirs(self.out_dir, exist_ok=True)
            masked_dem_path = os.path.join(self.out_dir, f'{self.colony_name}_clipped.tif')
            with rasterio.open(masked_dem_path, "w", **out_meta) as dest:
                dest.write(out_image)

            self.dem_path = masked_dem_path  # use clipped DEM

    def process(self, processes: list = ['slope', 'aspect', 'Roughness', 'TRI']) -> typing.List[DemProcessDesc]:
        print(f'Processing: {self.colony_name}')
        os.makedirs(self.out_dir, exist_ok=True)
        UseExceptions()

        # clip DEM before processing
        self._clip_dem_to_aoi()

        for p in processes:
            try:
                desc = self.dem_processes[p]
            except KeyError:
                raise RuntimeError(f"DEM process {p} not valid")

            print(f"\tDEM processing: {p}")
            desc.dataset = DEMProcessing(os.path.join(self.out_dir, f'{self.colony_name}_{p}.tif'),
                                        self.dem_path, p, **desc.kwargs)
            desc.metrics = self._get_ds_metrics(desc.dataset)
            desc.dataset = None
        return self.dem_processes

    @property
    def slope(self) -> DemProcessDesc:
        return self.dem_processes['slope']

    @property
    def aspect(self) -> DemProcessDesc:
        return self.dem_processes['aspect']

    @property
    def Roughness(self) -> DemProcessDesc:
        return self.dem_processes['Roughness']

    @property
    def TRI(self) -> DemProcessDesc:
        return self.dem_processes['TRI']


# colony to DEM dictionary
# match DEM to colony name
# might not work for Cape Crozier IKNONOS image to DEMs (cape_crozier_dem, ortho dem combines 17_33_2_2_2m_v2.0, 17_33_2_1_2m_v2.0, 17_33_1_2_2m_v2.0 and 17_33_1_1_2m_v2.0 and )

colony_to_dem = {
    'Beaufort_Island': '17_34_1_1_2m_v2.0',
    'Beaufort_North': '17_34_1_1_2m_v2.0',
    'Cape_Adare': 'cape_adare',  # ortho dem combines 10_34_2_1_2m_v2.0 and 10_34_2_2_2m_v2.0
    'Cape_Crozier': '17_33_2_2_2m_v2.0',  # not used for Crozier IKNONOS images
    'Cape_Bird': '17_34_1_1_2m_v2.0',
    'Cape_Hallet': '11_34_2_1_2m_v2.0',
    'Cape_Royds': '17_34_2_1_2m_v2.0',
    'Coulman_Middle': '13_34_1_1_2m_v2.0',
    'Coulman_North': '13_34_1_1_2m_v2.0',
    'Coulman_South': '13_34_1_1_2m_v2.0',
    'Downshire_Cliffs': '11_34_1_1_2m_v2.0',
    'Duke_Island': '11_34_1_1_2m_v2.0',
    'Inexpressible_Island': '15_35_1_2_2m_v2.0',
    'Possession_Island': '11_34_1_1_2m_v2.0',
    'Sven_Island': '11_34_1_1_2m_v2.0',
    'Terra_Nova': 'terra_nova'  # ortho dem combines 14_35_2_2_2m_v2.0 and 15_35_1_2_2m_v2.0
}

# main
if __name__ == "__main__":
    all_metrics = []

    for colony_name, dem_base_name in colony_to_dem.items():
        dem_path = os.path.join(DEM_FILES_DIR, f"{dem_base_name}_dem.tif")
        if os.path.isfile(dem_path):
            dem_metrics = DemMetrics(dem_path, colony_name, OUTPUT_DIR)
            dem_metrics.process(['slope', 'aspect', 'Roughness', 'TRI'])  # add dem proccesses here
            print(dem_metrics)
