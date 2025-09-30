## Contains:

The scripts for completing analysis in Strang et al. (in prep).

Inlabru model analyses:
- `Inlabru_test.R`: This script contains the R code to run the Poisson GLM. 
- `Inlabru_mesh.R`: This script contains the R code used to construst the 2D triangulated mesh.

Old:
- `LGCP_Nesi.R`: This script contains the R code currently being developed and run on NeSI - 20250910.
- `LGCP_test.R`: NEEDS UPDATING This script contains the R code for testing model sensitivity to the SPDE priors.
- `LGCP_candidates.R`: This script contains the R code for testing the LGCP candidate models.
- `LGCP_analysis.RMD`: NEEDS UPDATING

Pre-processing analyses:
- `getAnalysisData.py`: NOT COMPLETE/ NOT USED
- `procGuanoTerrain.py`: This script contains the python code to turn guano raster into percent cover at 2X2 m resolution.
- `procPointShapefiles.py`: This script contains the R code to turn UAV penguin points into xy data in EPSG:3031.
- `procUAVOrtho`: This folder contains the workflow to turn the UAV orthomosaic image into a shapefile/ polygon.
