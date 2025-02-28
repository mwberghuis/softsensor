# COMSOL + MATLAB Analysis

## Overview
This repository contains a COMSOL (Multiphysics 6.2) model and MATLAB (R2021a) scripts for postprocessing simulation results for the Goshtasbi et al. 2025 paper. Please consider citing this if you find the implementation useful. 

## Folder Structure
- **axisymsensor_results/**: Contains the (manually exported) .xlsx files with the result data of all model sweeps
- **comsol_model/**: Contains `.mph` file from COMSOL for one loadcase (load_mag=40kPa -> = F_max 4.8N) and fluid viscosity (mu=0.1 Pas = 100 cSt). Changing these parameters (e.g. in batch sweep mode) generates the result data.
- **axisymsensor_postprocessing_sweeps.m**: MATLAB script to load data from _results and generate figures.

## How to Use Comsol (Multiphysics 6.2) model and matlab script
0. Copy the unedited files provided before changing parameters. Then, change load/material/geometry etc. 
1. Run simulation in COMSOL.
2. Derived values (pressure at endpoints, gap, flux etc.) can be generated automatically by running the Sequence under Job configurations
3. Export results manually to `exported_results/` as .xlsx. Keep the same formatting as the provided result files. (exporting as .csv gives delimiter trouble while loading into matlab)
5. Copy the result file name into the matlab script. 
6. Run the MATLAB script in `scripts/` to generate figures.

## License
CC BY 4.0
