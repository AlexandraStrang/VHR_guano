# python script for modifying guano rasters
# creator: Alexandra Strang
# created: 2025

# Use rasterio to modify exported gunao rasters from ArcGIS Pro
# Cell size of rasters differ and depend on image
# Uses class name to modify as raster values and class values not consistent

import os
import rasterio
import numpy as np
import csv
from dbfread import DBF

# directories and paths (on SEE hard-drive)
raster_folder = r"D:\Guano_rasters"
vat_path = r"D:\Guano_rasters\VAT_guano_rasters"
output_csv = r"D:\Guano_rasters\raster_classes.csv"
output_folder = r"D:\Guano_rasters\Clean_guano_rasters"  # need to create
os.makedirs(output_folder, exist_ok=True)

# map class name from VAT (Value Attrivute Table) to new raster value
class_remap = {
    "Fresh guano": 1,
    "Shadowed guano": 1,
    "Old guano": 2,
    # all others default to 0
}


def load_vat_mapping(raster_name):
    """
    find and read VAT files
    """
    base_name = raster_name.replace(".tif", "")
    vat_file = os.path.join(vat_path, base_name + ".tif.vat.dbf")
    if not os.path.exists(vat_file):
        print(f"VAT not found for {raster_name}")
        return {}

    try:
        table = DBF(vat_file, encoding="latin1")
        return {
            rec["Value"]: {
                "Classvalue": rec.get("Classvalue", None),
                "Class_name": rec.get("Class_name", "Unknown").strip(),
                "RemapValue": class_remap.get(
                    rec.get("Class_name", "").strip(), 0)
                # everything else is 0
            }
            for rec in table
        }
    except Exception as e:
        print(f"Error reading {vat_file}: {e}")
        return {}


csv_rows = []  # to add raster values into a csv file

# loop through rasters
for filename in os.listdir(raster_folder):
    if filename.endswith(".tif"):
        filepath = os.path.join(raster_folder, filename)
        with rasterio.open(filepath) as src:
            data = src.read(1)  # only one band
            profile = src.profile
            nodata = profile.get("nodata", 255)  # set no data value to 255

        unique_vals = np.unique(data)
        class_map = load_vat_mapping(filename)

        print(f"\n{filename}: unique values: {unique_vals}")
        print("Class values:")

        remapped = np.full_like(data, 0)

        for val in unique_vals:
            if class_map and val in class_map:
                entry = class_map[val]
                remap_val = entry["RemapValue"]
                classval = entry["Classvalue"]
                classname = entry["Class_name"]
                print(f"  {val}: {classval} - {classname}")
            else:
                remap_val = 0
                classval = "Nodata"
                classname = "Nodata"
                print(f"  {val}: Not found in VAT (Nodata)")

            csv_rows.append({  # create csv file
                "Raster": filename,
                "PixelValue": val,
                "Classvalue": classval,
                "Class_name": classname
            })

            if val == nodata:
                remapped[data == val] = nodata
            else:
                remapped[data == val] = remap_val

        # save modified raster
        out_path = os.path.join(output_folder, f"cleaned_{filename}")
        profile.update(dtype=rasterio.uint8)

        with rasterio.open(out_path, "w", **profile) as dst:
            dst.write(remapped.astype(rasterio.uint8), 1)

        print(f"Processed and saved: {out_path}")

# export raster values
with open(output_csv, mode='w', newline='', encoding='utf-8') as csvfile:
    fieldnames = ["Raster", "PixelValue", "Classvalue", "Class_name"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(csv_rows)

print(f"CSV export complete: {output_csv}")
