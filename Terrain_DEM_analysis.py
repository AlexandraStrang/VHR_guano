# python script for calculating terrain variables from REMA DEMs
# creator: Alexandra Strang
# created: 2025

# Use conda environmenet pgc_test
# REMA DEM projection in WGS 1984 Antartic Polar Stereographic (EPSG 3031)
# https://gdal.org/en/stable/programs/gdaldem.html
# Could adapt script to add values into CSV
# Could mask by colony AOI for smaller tif files

import os
import typing
import numpy as np
from osgeo.gdal import DEMProcessing, UseExceptions, Dataset, Band
from dataclasses import dataclass

# find directories and paths
DEM_FILES_DIR = r'C:\Test_dem'
OUTPUT_DIR = r'C:\Test_dem'


@dataclass
class DemProcessDesc:
    kwargs: dict
    dataset: Dataset = None
    metrics: dict = None


class DemMetrics:
    """
    define dem proccess and algorithm
    """
    def __init__(self, dem_path: os.PathLike, colony_name: str, output_base_dir: os.PathLike):
        self.dem_path = dem_path
        self.colony_name = colony_name
        self.out_dir = os.path.join(output_base_dir, colony_name)
        self.dem_processes = {'slope': DemProcessDesc({'alg': 'Horn'}),  # Use Horn's method (3x3 neighbourhood of elevation values), -9999 is no data (or flat?)
                              'aspect': DemProcessDesc({}),  # aspect might need correcting (off 180?) Azimuth aspect where 0 degress is north
                              'Roughness': DemProcessDesc({}),  # roughness takes the difference between the maximum and minimum elevation in a 3×3 window around each pixel (elevation range of local varaibility and more sensitive to extreme changes)
                              'TRI': DemProcessDesc({'alg': 'Riley'})}  # Terrain Ruggedness Index (TRI) using the method developed by Riley et al. (1999), quantifies how “rough” the terrain, measuring elevation differences between a central pixel and its 8 neighbors

    def _get_ds_metrics(self, ds: Dataset):
        """
        dataset information and metrics
        """
        band: Band = ds.GetRasterBand(1)
        array: np.ndarray = band.ReadAsArray()
        return array.min(), array.max(), array.mean(), array.std()

    def process(self, processes: list = ['slope', 'aspect', 'Roughness', 'TRI']) -> typing.List[DemProcessDesc]:
        print(f'Processing: {self.colony_name}')
        os.makedirs(self.out_dir)
        UseExceptions()
        for p in processes:
            try:
                desc = self.dem_processes[p]
            except KeyError:
                raise RuntimeError(f"DEM process {p} not valid")

            print(f"\tDEM processing: {p}")
            desc.dataset = DEMProcessing(os.path.join(self.out_dir, f'{self.colony_name}_{p}.tif'),
                                         self.dem_path, p, **desc.kwargs)
            desc.metrics = self._get_ds_metrics(desc.dataset)
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
# does not work for Cape Crozier IKNONOS image to DEMs (cape_crozier_dem, ortho dem combines 17_33_2_2_2m_v2.0, 17_33_2_1_2m_v2.0, 17_33_1_2_2m_v2.0 and 17_33_1_1_2m_v2.0 and )

colony_to_dem = {
    'Beaufort_Island': '17_34_1_1_2m_v2.0',
    'Cape_Adare': 'cape_adare',  # ortho dem combines 10_34_2_1_2m_v2.0 and 10_34_2_2_2m_v2.0
    'Cape_Crozier': '17_33_2_2_2m_v2.0',  # not for Crozier IKNONOS images
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
