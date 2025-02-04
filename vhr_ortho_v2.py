
import os
import xml.etree
import xml.etree.ElementTree
import subprocess
from dataclasses import dataclass
import xml

PGC_ORTHO_IMAGE_UTILS_PATH = r"D:\Ortho\imagery_utils"
IMAGE_FILES_DIR = r'D:\Ortho\dem_groups'
DEM_FILES_DIR = r'D:\Ortho\rema_dems'
PYTHON_PATH = r'C:\Users\astra\miniforge3\envs\pgc_test\python.exe'
OUTPUT_DIR = r'D:\Ortho\ortho_images'
MAXAR_IMAGES_2022_DIR = r'D:\VHR_Images\Maxar_Images_20220420'
MAXAR_IMAGES_2024_DIR = r'D:\VHR_Images\Maxar_Images_20241219'

@dataclass
class TifXmlFilePaths:
    tif_path: str
    xml_metadata_path: str
    xml_readme_path: str

    def __post_init__(self): 
        assert all(os.path.isfile(f) for f in [self.tif_path, self.xml_metadata_path, self.xml_readme_path])

def run_pgc_ortho(dem_file_path, image_file_path, output_dir, dry_run: bool = False):
    output_dir = os.path.join(output_dir, os.path.splitext(os.path.basename(image_file_path))[0])
    os.makedirs(output_dir)
    command = [PYTHON_PATH,
               os.path.join(PGC_ORTHO_IMAGE_UTILS_PATH, 'pgc_ortho.py'),
               image_file_path,
               output_dir,
               "--epsg", "3031", # images need to match DEM projection of epsg=3031
               "--dem", dem_file_path]
    if dry_run:
        command.append("--dryrun")
    return subprocess.run(command).returncode == 0

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

def get_area_name(xml_readme_path: str):
    xml_tree = xml.etree.ElementTree.parse(xml_readme_path)
    root = xml_tree.getroot()
    area_element = root.find('AREADESCRIPTION')
    return area_element.text.replace(' ', '_')

def get_dem_file(tif_xml_paths: TifXmlFilePaths):
    colony_to_dem = {
        'Beaufort_Island': '17_34_1_1_2m_v2.0',
        'Cape_Adare': '10_34_2_1_2m_v2.0',
        'Cape_Crozier': '17_33_2_2_2m_v2.0',
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
        'Terra_Nova': '14_35_2_2_2m_v2.0',
    }
    area_name = get_area_name(tif_xml_paths.xml_readme_path)
    try:
        dem_file = next(v for k,v in colony_to_dem.items() if k in area_name) + "_dem.tif"
    except StopIteration:
        raise RuntimeError(f"Area name not matched: {area_name}")
    dem_path = os.path.join(DEM_FILES_DIR, dem_file)
    assert os.path.isfile(dem_path), f'{dem_file} does not exist'
    return dem_path

for f in get_tif_xml_files(MAXAR_IMAGES_2022_DIR) + get_tif_xml_files(MAXAR_IMAGES_2024_DIR):
    status = run_pgc_ortho(get_dem_file(f), f.tif_path, OUTPUT_DIR)
    print(f"{'OK' if status else 'FAILED':<10}" f"{f.tif_path}")
