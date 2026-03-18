## Contains:

The scripts for completing analysis to address question 2 in Strang et al. (in prep).

Inlabru model analyses:
- `01_Inlabru_mesh.R`: This script contains the R code used to construst the 2D triangulated mesh.
- `02_Inlabru_candidates.R`: This script contains the R code to run 8 candidate models for 2020 data.
- `03_Inlabru_diagnostics.R`: This script contains the R code to calculated the expected abundances, residuals, and diagnostic statistics for the 8 candidate models.
- `04_Inlabru_partial_predictions.R`: This script contains the R code to compute partial predictions and plots for each of the 8 candidate models.
- `05_Inlabru_2019_predictions.R`: This script contains the R code to use the 2020 trained candidate models to predict abundance for 2019.

Pre-processing analyses:
- `procGuanoTerrain.py`: This script contains the python code to turn guano raster into percent cover at 2X2 m resolution.
- `procPointShapefiles.py`: This script contains the R code to turn UAV penguin points into xy data in EPSG:3031.
