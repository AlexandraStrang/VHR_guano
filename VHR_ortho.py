# python script for batch orthorectification of VHR imagery
# creator: Alexandra Strang
# created: 2025

# use PGC GitHub: orthorectification script https://www.pgc.umn.edu/guides/pgc-coding-and-utilities/using-pgc-github-orthorectification/
# open folder https://github.com/PolarGeospatialCenter/imagery_utils
# to run this script the python interpreter come from imagery utils
# conda environmenet pgc_test
# reprojects image in WGS 1984 Antartic Polar Stereographic to match REMA DEM projection

import os
import xml.etree
import xml.etree.ElementTree
import subprocess
from dataclasses import dataclass
import xml

# find directories and paths
PGC_ORTHO_IMAGE_UTILS_PATH = r'D:\Ortho\imagery_utils'
DEM_FILES_DIR = r'D:\Ortho\rema_dems'
PYTHON_PATH = r'C:\Users\astra\miniforge3\envs\pgc_test\python.exe'
OUTPUT_DIR = r'D:\Ortho\ortho_16_rf'

# run for DG images
MAXAR_IMAGES_2022_DIR = r'D:\VHR_Images\Maxar_Images_20220420'
MAXAR_IMAGES_2024_DIR = r'D:\VHR_Images\Maxar_Images_20241219'

# find metadata paths
@dataclass
class TifXmlFilePaths:
    tif_path: str
    xml_metadata_path: str
    xml_readme_path: str

    def __post_init__(self): 
        assert all(os.path.isfile(f) for f in [self.tif_path, self.xml_metadata_path, self.xml_readme_path])

# command for running pgc ortho script
def run_pgc_ortho(dem_file_path, image_file_path, output_dir, dry_run: bool = False):
    output_dir = os.path.join(output_dir, os.path.splitext(os.path.basename(image_file_path))[0])
    os.makedirs(output_dir)
    command = [PYTHON_PATH,
               os.path.join(PGC_ORTHO_IMAGE_UTILS_PATH, 'pgc_ortho.py'),
               image_file_path,
               output_dir, # default resample is nearest neighbour
               "--epsg", "3031", # needs to match DEM projection of epsg=3031
               "-t", "UInt16", # 16-bit output
               "-c", "rf", # reflectance correction
               "--dem", dem_file_path]
    if dry_run:
        command.append("--dryrun")
    return subprocess.run(command).returncode == 0

# match tif and xml readme
def get_tif_xml_files(path) -> list[TifXmlFilePaths]:
    tif_files = []
    for root, _, files in os.walk(path):
        for file in files:
            if file.lower().endswith('.tif'):
                xml_file = os.path.splitext(file)[0] + '.xml'
                read_me_file = next(f for f in os.listdir(os.path.dirname(root)) if f.lower().endswith('readme.xml'))
                read_me_path = os.path.join(os.path.dirname(root), read_me_file)
                tif_xml_paths = TifXmlFilePaths(os.path.join(root, file), os.path.join(root, xml_file), read_me_path)
                tif_files.append(tif_xml_paths)
    return tif_files

# area name for DEM matching from xml readme
def get_area_name(xml_readme_path: str):
    xml_tree = xml.etree.ElementTree.parse(xml_readme_path)
    root = xml_tree.getroot()
    area_element = root.find('AREADESCRIPTION')
    return area_element.text.replace(' ', '_')

# match DEM to images using area name from xml readme
# image to DEM dictionary
def get_dem_file(tif_xml_paths: TifXmlFilePaths):
    colony_to_dem = {
        'Beaufort_Island': '17_34_1_1_2m_v2.0',
        'Cape_Adare': '10_34_2_1_2m_v2.0', # not needed
        'Cape_Crozier': '17_33_2_2_2m_v2.0',
        'Cape_Bird': '17_34_1_1_2m_v2.0',
        'Cape_Hallet': '11_34_2_1_2m_v2.0', # typo in order Hallett
        'Cape_Royds': '17_34_2_1_2m_v2.0',
        'Coulman_Middle': '13_34_1_1_2m_v2.0',
        'Coulman_North': '13_34_1_1_2m_v2.0',
        'Coulman_South': '13_34_1_1_2m_v2.0',
        'Downshire_Cliffs': '11_34_1_1_2m_v2.0',
        'Duke_Island': '11_34_1_1_2m_v2.0',
        'Inexpressible_Island': '15_35_1_2_2m_v2.0',
        'Possession_Island': '11_34_1_1_2m_v2.0',
        'Sven_Island': '11_34_1_1_2m_v2.0',
        'Terra_Nova': '14_35_2_2_2m_v2.0', # not needed
    }
    area_name = get_area_name(tif_xml_paths.xml_readme_path)
    try:
        dem_file = next(v for k,v in colony_to_dem.items() if k in area_name) + "_dem.tif"
    except StopIteration:
        raise RuntimeError(f"Area name not matched: {area_name}")
    dem_path = os.path.join(DEM_FILES_DIR, dem_file)
    assert os.path.isfile(dem_path), f'{dem_file} does not exist'
    return dem_path

# run for DG images
for f in get_tif_xml_files(MAXAR_IMAGES_2022_DIR) + get_tif_xml_files(MAXAR_IMAGES_2024_DIR):
    status = run_pgc_ortho(get_dem_file(f), f.tif_path, OUTPUT_DIR)
    print(f"{'OK' if status else 'FAILED':<10}" f"{f.tif_path}")

"""
Run Terra Nova images with mosaic DEM
TERRA_NOVA_DEM_FILE_DIR = "D:\Ortho\rema_dems\Terra_nova\terra_nova.tif"
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Terra_nova\terra_nova.tif D:\VHR_Images\Maxar_Images_20241219\050272141450_01\050272141450_01_P001_PSH\15DEC29211407-S2AS-050272141450_01_P001.tif D:\Ortho\ortho_16_rf\15DEC29211407-S2AS-050272141450_01_P001
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Terra_nova\terra_nova.tif D:\VHR_Images\Maxar_Images_20241219\050272141460_01\050272141460_01_P001_PSH\21JAN05215921-S2AS-050272141460_01_P001.tif D:\Ortho\ortho_16_rf\21JAN05215921-S2AS-050272141460_01_P001
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Terra_nova\terra_nova.tif D:\VHR_Images\Maxar_Images_20241219\050272141470_01\050272141470_01_P001_PSH\24FEB01211801-S2AS-050272141470_01_P001.tif D:\Ortho\ortho_16_rf\24FEB01211801-S2AS-050272141470_01_P001

Warning 1: The definition of projected CRS EPSG:3031 got from GeoTIFF keys is not the same as the one from the EPSG registry, which may cause issues during reprojection operations. Set GTIFF_SRS_SOURCE configuration option to EPSG to use official parameters (overriding the ones from GeoTIFF keys), or to GEOKEYS to use custom values from GeoTIFF keys and drop the EPSG code.

"""
"""
Run Cape Adare images with mosaic DEM
CAPE_ADARE_DEM_FILE_DIR = "D:\Ortho\rema_dems\Cape_adare\cape_adare.tif"
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_adare\cape_adare.tif D:\VHR_Images\Maxar_Images_20241219\050272141050_01\050272141050_01_P001_PSH\12JAN08210317-S2AS-050272141050_01_P001.tif D:\Ortho\ortho_16_rf\12JAN08210317-S2AS-050272141050_01_P001
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_adare\cape_adare.tif D:\VHR_Images\Maxar_Images_20241219\050272141060_01\050272141060_01_P001_PSH\15DEC11203632-S2AS-050272141060_01_P001.tif D:\Ortho\ortho_16_rf\15DEC11203632-S2AS-050272141060_01_P001
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_adare\cape_adare.tif D:\VHR_Images\Maxar_Images_20241219\050272141070_01\050272141070_01_P001_PSH\20DEC15215632-S2AS-050272141070_01_P001.tif D:\Ortho\ortho_16_rf\20DEC15215632-S2AS-050272141070_01_P001
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_adare\cape_adare.tif D:\VHR_Images\Maxar_Images_20241219\050272141080_01\050272141080_01_P001_PSH\21JAN04205608-S2AS-050272141080_01_P001.tif D:\Ortho\ortho_16_rf\21JAN04205608-S2AS-050272141080_01_P001
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_adare\cape_adare.tif D:\VHR_Images\Maxar_Images_20241219\050272141090_01\050272141090_01_P001_PSH\21JAN16220016-S2AS-050272141090_01_P001.tif D:\Ortho\ortho_16_rf\21JAN16220016-S2AS-050272141090_01_P001
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_adare\cape_adare.tif D:\VHR_Images\Maxar_Images_20241219\050272141100_01\050272141100_01_P001_PSH\21FEB09213924-S2AS-050272141100_01_P001.tif D:\Ortho\ortho_16_rf\21FEB09213924-S2AS-050272141100_01_P001
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_adare\cape_adare.tif D:\VHR_Images\Maxar_Images_20241219\050272141110_01\050272141110_01_P001_PSH\21FEB13210309-S2AS-050272141110_01_P001.tif D:\Ortho\ortho_16_rf\21FEB13210309-S2AS-050272141110_01_P001
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_adare\cape_adare.tif D:\VHR_Images\Maxar_Images_20241219\050272141120_01\050272141120_01_P001_PSH\23DEC03220040-S2AS-050272141120_01_P001.tif D:\Ortho\ortho_16_rf\23DEC03220040-S2AS-050272141120_01_P001

Warning 1: The definition of projected CRS EPSG:3031 got from GeoTIFF keys is not the same as the one from the EPSG registry, which may cause issues during reprojection operations. Set GTIFF_SRS_SOURCE configuration option to EPSG to use official parameters (overriding the ones from GeoTIFF keys), or to GEOKEYS to use custom values from GeoTIFF keys and drop the EPSG code.

"""
"""
Run IKO1 images (Cape Crozier)
IK01_IMAGE_1_DIR = r'D:\VHR_Images\Maxar_Ikonos_20241219\Cape_Crozier_1\2011120221030920000011607365\po_1495866_0000000\po_1495866_bgrn_0000000.tif'
IK01_IMAGE_2_DIR = r'D:\VHR_Images\Maxar_Ikonos_20241219\Cape_Crozier_4\2012022921450540000011624223\po_1505906_0000000\po_1505906_bgrn_0000000.tif'
IK01_DEM_FILE_DIR = r'D:\Ortho\rema_dems\Cape_crozier\cape_crozier.tif'

python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_crozier\cape_crozier.tif D:\VHR_Images\Maxar_Ikonos_20241219\Cape_Crozier_1\2011120221030920000011607365\po_1495866_0000000\po_1495866_bgrn_0000000.tif D:\Ortho\ortho_16_rf\2011120221030920000011607365
python pgc_ortho.py -p 3031 -t UInt16 -c rf -d D:\Ortho\rema_dems\Cape_crozier\cape_crozier.tif D:\VHR_Images\Maxar_Ikonos_20241219\Cape_Crozier_4\2012022921450540000011624223\po_1505906_0000000\po_1505906_bgrn_0000000.tif D:\Ortho\ortho_16_rf\2012022921450540000011624223

Warning 1: The definition of projected CRS EPSG:3031 got from GeoTIFF keys is not the same as the one from the EPSG registry, which may cause issues during reprojection operations. Set GTIFF_SRS_SOURCE configuration option to EPSG to use official parameters (overriding the ones from GeoTIFF keys), or to GEOKEYS to use custom values from GeoTIFF keys and drop the EPSG code.

Manually find IK metadata
 match = re.search(regex, os.path.basename(metafile.lower()))
    if match:
        siid = "INSERT Source Image ID"

"""
