# python script for batch orthorectification of VHR imagery
# use PGC GitHub: orthorectification script

import os
import numpy as np
#from osgeo import gdal
import subprocess
import pandas as pd # PS C:\Users\ajs424> py -m pip install pandas

# input and output directories
maxar_images = r'C:\Local\Ajs424\PhD ADPE project\Ortho\vhr_images'
ortho_images = r'C:\Local\Ajs424\PhD ADPE project\Ortho\ortho_images'
dem_files = r'C:\Local\Ajs424\PhD ADPE project\Ortho\rema_dems'

# map colony to corresponding DEM files
colony_to_dem = {
    'Beaufort_Island': '17_34_1_1_2m_v2.0',
    'Cape_Adare': '10_34_2_1_2m_v2.0',
    'Cape_Crozier': '17_33_2_2_2m_v2.0',
    'Cape_Bird': '17_34_1_1_2m_v2.0',
    'Cape_Hallett': '11_34_2_1_2m_v2.0',
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

def lookup_dem(image):
    name = "_".join([i for i in image.split('_') if i.isalpha()])
    matches = [v for k,v in colony_to_dem.items() if name in k]
    if len(matches) > 1:
        raise AssertionError(f"{image} matches more than one dem: {matches}")
    if not matches:
        raise AssertionError(f"{image} does not match any dem")
    return matches[0]

print(lookup_dem('Coulman_2011_1.tif')) # test if work (should be incorrect)

# list all VHR images in the input directory
vhr_files = [i for i in os.listdir(maxar_images) if i.endswith('.tif')]

# run orthorectification script for each VHR image
for image in vhr_files:
    image_path = os.path.join(maxar_images, image)
    dem_path = os.path.join(dem_files, lookup_dem(image))
    output_path = os.path.join(ortho_images, f'ortho_{image}')
# Command to run PGC orthorectification script
    command = f'C:\Local\Ajs424\PhD ADPE project\Ortho\pgc_ortho.py -input {image_path} -output {output_path} -dem {dem_path}'
    subprocess.run(command, shell=True)
        
