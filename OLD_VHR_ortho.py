# python script for batch orthorectification of VHR imagery
# creator: Alexandra Strang
# created: 2025

# use PGC GitHub: orthorectification script https://www.pgc.umn.edu/guides/pgc-coding-and-utilities/using-pgc-github-orthorectification/
# https://github.com/PolarGeospatialCenter/imagery_utils
# to run this script the python interpreter must ensure that gdal can be imported from osgeo
# use images projected in WGS 1984 APS to match REMA DEM projection

import os
from osgeo import gdal
import subprocess
import numpy as np

# find input and output directories
maxar_images = r'D:\Ortho\vhr_images_3031'
ortho_images = r'D:\Ortho\ortho_images'
dem_files = r'D:\Ortho\rema_dems'

# map colony to corresponding DEM files
colony_to_dem = {
    'Beaufort_Island': '17_34_1_1_2m_v2.0',
    'Cape_Adare': '10_34_2_1_2m_v2.0',
    'Cape_Crozier': '17_33_2_2_2m_v2.0',
    'Cape_Bird': '17_34_1_1_2m_v2.0',
    'Cape_Hallett': '11_34_2_1_2m_v2.0', # DEM wasn't in folder?
    'Cape_Royds': '17_34_2_1_2m_v2.0',
    'Coulman_Middle': '13_34_1_1_2m_v2.0',
    'Coulman_North': '13_34_1_1_2m_v2.0',
    'Coulman_South': '13_34_1_1_2m_v2.0',
    'Downshire_Cliffs': '11_34_1_1_2m_v2.0',
    'Duke_Island': '11_34_1_1_2m_v2.0',
    'Inexpressible_Island': '15_35_1_2_2m_v2.0',
    'Possession_Island': '11_34_1_1_2m_v2.0',
    'Sven_Island': '11_34_1_1_2m_v2.0',
    'Terra_Nova': '14_35_2_2_2m_v2.0',
}

# function to match images to DEMs
def lookup_dem(image):
    name = '_'.join([i for i in image.split('_') if i.isalpha()])
    matches = [v for k,v in colony_to_dem.items() if name in k]
    if len(matches) > 1:
        raise AssertionError(f'{image} matches more than one dem: {matches}')
    if not matches:
        raise AssertionError(f'{image} does not match any dem')
    return matches[0]

# check that colonies and DEMs match up 
#print(lookup_dem('Coulman_2011_1.tif')) # should be incorrect
print(lookup_dem('Crozier_2011_1.tif')) # should print Crozier DEM name

# list all VHR images in the input directory
vhr_files = [i for i in os.listdir(maxar_images) if i.lower().endswith('.tif')]

# find pgc_ortho.py script
PGC_ORTHO_IMAGE_UTILS_PATH = r"D:\Ortho\imagery_utils"
# pgc_ortho_path = os.path.join('D:', 'Ortho', "pgc_ortho.py") # path to PGC script

# find OSGeo4w path for python interpeter 
# env_dir = r'C:\OSGeo4W\bin\python.exe'

# run script for each VHR image
for image in vhr_files:
    image_path = os.path.join(maxar_images, image)
    try:
        python_path = r'C:\Users\astra\miniforge3\envs\pgc_test\python.exe'
        dem_path = os.path.join(dem_files, f'{lookup_dem(image)}_dem.tif')
        output_path = os.path.join(ortho_images)
        command = [python_path, os.path.join(PGC_ORTHO_IMAGE_UTILS_PATH, 'pgc_ortho.py')]
        # Command to run PGC orthorectification script
        command = [ python_path, os.path.join(PGC_ORTHO_IMAGE_UTILS_PATH, 'pgc_ortho.py'),
                   #image_path, output_path,
                r"D:\Maxar_Images_20220420\014666468010_01\014666468010_01_P001_PSH\10JAN19200728-S2AS-014666468010_01_P001.TIF", 
                output_path,
            "--epsg", "3031", # images need to match DEM projection of epsg=3031
            "--dem", dem_path,
            
        ]

        # run the script with handling error
        subprocess.run(command, check=True)
    except AssertionError as e:
        print(f"Skipping {image}: {e}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {image}: {e}")
    break
