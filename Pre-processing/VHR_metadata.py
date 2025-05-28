# python script for extracting metadata from VHR imagery
# creator: Alexandra Strang
# created: 2025

# Extract satellite image metadata from image xml metadata and readmes (and txt for IKONOS images)
# Merge metadata for images with masterdata spreadsheet containing guano analysis data

import os
import xml.etree
import xml.etree.ElementTree
from dataclasses import dataclass
import xml
import csv
from collections import OrderedDict
import pandas as pd

# find directories
VHR_METADATA_OUTPUT_DIR = r"C:\Users\astra\OneDrive - University of Canterbury\ANTA - PhD\VHR_metadata"

# run for DG images
MAXAR_IMAGES_2022_DIR = r"C:\Users\astra\OneDrive - University of Canterbury\ANTA - PhD\Maxar images\Maxar_Images_20220420"
MAXAR_IMAGES_2024_DIR = r"C:\Users\astra\OneDrive - University of Canterbury\ANTA - PhD\Maxar images\Maxar_Images_20241219"
MAXAR_IMAGES_2024_IK02_DIR = r"C:\Users\astra\OneDrive - University of Canterbury\ANTA - PhD\Maxar images\Maxar_Ikonos_20241219"


# metadata and readme paths for DG images
@dataclass
class TifXmlFilePaths:
    tif_path: str
    xml_metadata_path: str
    xml_readme_path: str

    def __post_init__(self):
        """
        check if all paths exist
        """
        assert all(os.path.isfile(f) for f in [self.tif_path, self.xml_metadata_path, self.xml_readme_path])


# metadata for IKONOS-2 images
@dataclass
class TifTxtFilePaths:
    IK02_tif_path: str
    metadata_txt_path: str

    def __post_init__(self):
        """
        check if all paths exist
        """
        assert all(os.path.isfile(f) for f in [self.IK02_tif_path, self.metadata_txt_path]), \
            f"Missing file(s): {self.IK02_tif_path}, {self.metadata_txt_path}"


def get_tif_xml_files(path) -> list[TifXmlFilePaths]:
    """
    match tif with xml readme and xml metadata for DG images
    """
    tif_files = []
    for root, _, files in os.walk(path):
        for file in files:
            if file.lower().endswith('.tif'):
                xml_file = os.path.splitext(file)[0] + '.xml'
                read_me_file = next(
                    f for f in os.listdir(
                        os.path.dirname(root)) if f.lower().endswith('readme.xml'))
                read_me_path = os.path.join(
                    os.path.dirname(root), read_me_file)
                tif_xml_paths = TifXmlFilePaths(
                    os.path.join(root, file), os.path.join(root, xml_file), read_me_path)
                tif_files.append(tif_xml_paths)
    return tif_files


def get_IKO2_txt_files(path) -> list[TifTxtFilePaths]:
    """
    match tif and txt metadata for IKONOS-2 images
    """
    tif_files = []
    for root, _, files in os.walk(path):
        for file in files:
            if file.lower().endswith('bgrn_0000000.tif'):
                tif_path = os.path.join(root, file)
                # remove the 'bgrn_0000000' part to match the metadata file name pattern
                base_name = file.lower().replace('bgrn_0000000.tif', '').rstrip('_')
                txt_file = base_name + '_metadata.txt'
                txt_path = os.path.join(root, txt_file)

                # print statement to confirm matched file
                print(f"Matched image: {tif_path}")
                print(f"Matched metadata: {txt_path}")

                tif_txt_paths = TifTxtFilePaths(tif_path, txt_path)
                tif_files.append(tif_txt_paths)
    return tif_files


matched_files = get_IKO2_txt_files(MAXAR_IMAGES_2024_IK02_DIR)
print(f"\nTotal matched files: {len(matched_files)}")


def get_area_name(xml_readme_path: str):
    """
    get area name from xml readme for colonies
    """
    xml_tree = xml.etree.ElementTree.parse(xml_readme_path)
    root = xml_tree.getroot()
    area_element = root.find('AREADESCRIPTION')
    return area_element.text.replace(' ', '_')


def parse_metadata(xml_metadata_path: str) -> dict:
    """
    get metadata variables from xml files for DG images
    """
    xml_tree = xml.etree.ElementTree.parse(xml_metadata_path)
    root = xml_tree.getroot()
    return_value = OrderedDict()
    for metadata in ['CATID', 'SATID']:
        value = root.find('.//' + metadata)
        if value is None:
            raise RuntimeError(f"could not find {metadata}")
        return_value[metadata] = value.text
    for metadata in ['MEANCOLLECTEDGSD', 'MEANSUNAZ', 'MEANSUNEL',
                     'MEANOFFNADIRVIEWANGLE', 'MEANINTRACKVIEWANGLE',
                     'EXPOSUREDURATION', 'ROWUNCERTAINTY', 'COLUNCERTAINTY']:
        value = root.find('.//' + metadata)
        if value is None:
            raise RuntimeError(f"could not find {metadata}")
        return_value[metadata] = float(value.text)
    return return_value


def parse_ikonos_metadata(txt_path: str) -> dict:
    """
    get metadata variables from txt files for IKONOS-2 images
    """
    metadata = {}
    ms_cross = None
    ms_along = None

    with open(txt_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("Source Image ID:"):
            metadata["IK_ID"] = line.split(":")[1].strip()  # temporarily store
        elif line.startswith("Sensor:"):
            metadata["SATID"] = line.split(":")[1].strip()
        elif "MS Cross Scan" in line:
            ms_cross = float(line.split(":")[1].strip().split()[0])
        elif "MS Along Scan" in line:
            ms_along = float(line.split(":")[1].strip().split()[0])
        elif "Sun Angle Azimuth" in line:
            metadata["MEANSUNAZ"] = float(line.split(":")[1].strip().split()[0])
        elif "Sun Angle Elevation" in line:
            metadata["MEANSUNEL"] = float(line.split(":")[1].strip().split()[0])

    if "IK_ID" not in metadata:
        raise RuntimeError(f"Could not find 'Source Image ID' in {txt_path}")
    if "SATID" not in metadata:
        raise RuntimeError(f"Could not find 'Sensor' in {txt_path}")
    if ms_cross is None or ms_along is None:
        raise RuntimeError(f"Could not find MS Cross and/or Along Scan GSD in {txt_path}")
    if "MEANSUNAZ" not in metadata or "MEANSUNEL" not in metadata:
        raise RuntimeError(f"Could not find Sun Angles in {txt_path}")

    # have to calculate mean GSD
    meancollectedgsd = round((ms_cross + ms_along) / 2, 3)

    # CATID mapping based on known IK_IDs
    ik_to_dg_id = {
        "2011120221030920000011607365": "10600100072D0200",
        "2012022921450540000011624223": "106001000753CF00"
    }

    ik_id = metadata["IK_ID"]
    if ik_id not in ik_to_dg_id:
        raise RuntimeError(f"Could not find CATID match for IK_ID: {ik_id} in {txt_path}")
    metadata["CATID"] = ik_to_dg_id[ik_id]

    # manually add MEANOFFNADIRVIEWANGLE based on known values from MAXAR archive
    off_nadir_lookup = {
            "10600100072D0200": 28.6,
            "106001000753CF00": 27.1
        }
    metadata["MEANOFFNADIRVIEWANGLE"] = off_nadir_lookup.get(metadata["CATID"], None)

    # build OrderedDict to match DG format
    ordered = OrderedDict()
    ordered["CATID"] = metadata["CATID"]
    ordered["SATID"] = metadata["SATID"]
    ordered["MEANCOLLECTEDGSD"] = meancollectedgsd
    ordered["MEANSUNAZ"] = metadata["MEANSUNAZ"]
    ordered["MEANSUNEL"] = metadata["MEANSUNEL"]
    ordered["MEANOFFNADIRVIEWANGLE"] = metadata["MEANOFFNADIRVIEWANGLE"]
    ordered["MEANINTRACKVIEWANGLE"] = None
    ordered["EXPOSUREDURATION"] = None
    ordered["ROWUNCERTAINTY"] = None
    ordered["COLUNCERTAINTY"] = None

    return ordered


# assume this returns a list of TifTxtFilePaths with valid paths
ikonos_files = get_IKO2_txt_files(MAXAR_IMAGES_2024_IK02_DIR)


# export DG image metadata into CSV
with open(os.path.join(VHR_METADATA_OUTPUT_DIR, "VHR_metadata.csv"), 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    files = get_tif_xml_files(MAXAR_IMAGES_2022_DIR) + get_tif_xml_files(MAXAR_IMAGES_2024_DIR)
    csvwriter.writerow(['Area_desc'] + [j for j in parse_metadata(files[0].xml_metadata_path).keys()])
    for i in files:
        csvwriter.writerow([get_area_name(i.xml_readme_path)] + [j for j in parse_metadata(i.xml_metadata_path).values()])
# puts metadata values in csv matched to area name (colony) for every tif

# export ikonos metadata into csv
with open(os.path.join(VHR_METADATA_OUTPUT_DIR, "IKONOS_metadata.csv"), 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')

    # use the first parsed metadata dict to define CSV headers
    example_metadata = parse_ikonos_metadata(ikonos_files[0].metadata_txt_path)
    csvwriter.writerow(['Area_desc'] + list(example_metadata.keys()))

    # write one row per file
    for file in ikonos_files:
        Area_desc = 'Cape_Crozier'  # for mergering
        metadata_values = parse_ikonos_metadata(file.metadata_txt_path)
        csvwriter.writerow([Area_desc] + list(metadata_values.values()))

# merge with masterdata csv
# path to masterdata and merged csv output location
MASTERDATA_DIR = r"C:\Users\astra\OneDrive - University of Canterbury\ANTA - PhD\Data\Data sheets"

# ensure output directory exists
os.makedirs(MASTERDATA_DIR, exist_ok=True)

# full file paths
METADATA_PATH = os.path.join(VHR_METADATA_OUTPUT_DIR, "VHR_metadata.csv")
IKONOS_METADATA_PATH = os.path.join(VHR_METADATA_OUTPUT_DIR, "IKONOS_metadata.csv")
MASTERDATA_PATH = os.path.join(MASTERDATA_DIR, "Strang_PhD_masterdata.csv")
MERGED_CSV_PATH = os.path.join(MASTERDATA_DIR, "Merged_masterdata.csv")

# read in both CSVs
VHR_metadata = pd.read_csv(METADATA_PATH, encoding='ISO-8859-1')
Ikonos_metadata = pd.read_csv(IKONOS_METADATA_PATH, encoding='ISO-8859-1')
masterdata = pd.read_csv(MASTERDATA_PATH, encoding='ISO-8859-1')

# concatenate DG and IKonos metadata into one dataframe
all_metadata = pd.concat([VHR_metadata, Ikonos_metadata], ignore_index=True)

# merge with suffixes using masterdata["ID", "Area_desc"] and VHR_metadata["CATID", "Area_desc"]
merged = pd.merge(
    masterdata,
    all_metadata,
    left_on=["ID", "Area_desc"],
    right_on=["CATID", "Area_desc"],
    how='left',
    suffixes=('_masterdata', '_metadata'))

# check columns to see duplicates
print(merged.columns)

# drop the duplicated column from the masterdata csv
merged = merged.drop(columns=[
    'Sensor', 'ID'])

# rename SATID column to sensor
merged = merged.rename(columns={
    'SATID': 'Sensor',
    'CATID': 'ID'
    })

merged.to_csv(MERGED_CSV_PATH, index=False)
print(f"Merged data saved to: {MERGED_CSV_PATH}")
