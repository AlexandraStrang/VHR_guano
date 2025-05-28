# python script for checking unique values in raster files and plotting
# creator: Alexandra Strang
# created: 2025

# Checks for unique raster values (after cleaning should be 0,1,2)
# set 255 = no data or background

import os
import rasterio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# directories and paths
raster_folder = r"D:\Guano_rasters"
vat_path = r"D:\Guano_rasters\VAT_guano_rasters"
output_csv = r"D:\Guano_rasters\raster_classes.csv"
cleaned_rasters = r"D:\Guano_rasters\Clean_guano_rasters"

# colormap for plotting values 0, 1, 2
cmap = ListedColormap(['lightgray', 'pink', 'blue'])
bounds = [-0.5, 0.5, 1.5, 2.5]
norm = plt.Normalize(vmin=0, vmax=2)

# loop through rasters
for filename in os.listdir(cleaned_rasters):
    if filename.endswith(".tif"):
        filepath = os.path.join(cleaned_rasters, filename)
        with rasterio.open(filepath) as src:
            data = src.read(1)

            # mask out value 255
            masked_data = np.ma.masked_where(data == 255, data)
            unique_vals = np.unique(masked_data.compressed())

            print(f"\n{filename}: unique values (excluding 255): {unique_vals}")

            # plot rasters
            plt.figure(figsize=(6, 6))
            plt.imshow(masked_data, cmap=cmap, norm=norm)
            plt.colorbar(ticks=[0, 1, 2], label='Class')
            plt.title(f"{filename}")
            plt.axis('off')
            plt.show()

# loop through rasters (option without plotting)
# for filename in os.listdir(cleaned_rasters):
#     if filename.endswith(".tif"):
#         filepath = os.path.join(cleaned_rasters, filename)
#         with rasterio.open(filepath) as src:
#             data = src.read(1)
#             unique_vals = np.unique(data)

#         print(f"\n{filename}: unique values: {unique_vals}")
