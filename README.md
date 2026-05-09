# UKESM1.1_MH_NorthAmerica

This repository contains code and data from Wilson et al., North American forest dieback simulated in response to warm and dry mid-Holocene summers.

## Notebooks

`Notebooks/` contains the code to produce the figures in this study. Each notebook is named for the figures it creates. \
Also included are notebooks that generate a compilation of climate reconstruction anomalies from the LegacyClimate dataset and calculate RMSE scores between model outputs and proxies.

**N.B.** Some notebooks reference files that are part of publicly available datasets used in this study but not inlcuded in this repository

## bin

Contains python modules with frequently used functions.

The `get_cpt` package (Bouziotas, 2020) is used for colormaps and can be found here: https://doi.org/10.5281/zenodo.3703160. \
`geometry_fix.py` applies a patch to fix issues with cartopy and shapely 2.x. See https://github.com/SciTools/cartopy/issues/2176 for more details.

## Data

Only new datasets created as part of this study are included here.

**Simulation Output:**
NetCDF outputs from UKESM1.1 mid-Holocene and pre-Industrial simulations are archived on Zenodo at \
`NETCDF_FILES.md` describes the files available on Zenodo \
**N.B.** Additional model outputs can be made available upon request.

`Data/Models/`:

- `GP_pr_ssp_anoms.xslx` -> Spreadsheet of monthly average Great Plains precipitation anomalies from SSP projection simulations.

NetCDF outputs from UKESM1.1 mid-Holocene and pre-Industrial simulations are archived on Zenodo at \
**N.B.** Additional model outputs can be made available upon request.

`Data/Proxies/`:

- `LC_climate_reconstruction_anomalies.xlsx` -> Our compilation of climate reconstruction anomalies from the LegacyClimate dataset.
- `pr_rmse_vals.xlsx` -> RMSE scores for model vs proxy precipitation anomalies.
- `tas_rmse_vals.xslx` -> RMSE scores for model vs proxy surface air temperature anomalies.
- `veg_rmse_vals.xslx` -> RMSE scores for model vs proxy vegetation surface fraction anomalies.

Publicly available files referenced in the code but not found in this repository include:

- Bartlein et al., (2011) proxy datasets -> Available at https://doi.org/10.1007/s00382-010-0904-1
  - `map_delta_06ka_ALL_grid_2x2_ex.nc`
  - `mat_delta_06ka_ALL_grid_2x2_ex.nc`

- LegacyClimate datasets -> Available at https://doi.pangaea.de/10.1594/PANGAEA.930512 (see 'Reconstruction files in .csv format')
  - `climate_reconstruction_<continent>.csv`

- Dawson et al., (2015) Vegetation Reconstructions -> Available at https://doi.org/10.5061/dryad.c2fqz61m5
  - `REVEALS_LCT_gridded.csv`

- CMIP6/PMIP4 Simulation Outputs downloaded from the Earth System Grid Federation (ESGF) nodes (https://esgf.github.io/nodes.html)
