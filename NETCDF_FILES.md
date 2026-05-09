# UKESM1.1_MH_NorthAmerica

This file describes netCDF files available in this repository. Each file contains detailed metadata on the variables contained within

**N.B.** Other simulation outputs are available upon request

`Data/Models`:

- `lsm_um13.2.nc` -> Land-Sea Mask
- `sat_vsmc_um13.2.nc` -> Volumetric soil moisture content at saturation [m3 m-3]

`Data/Models/Monthly_Means`:

- `MH_monthly_means_197001-204912.nc` -> Monthly mean outputs from mid-Holocene simulation representing an average over 80 model years
- `PI_monthly_means_187001-194912.nc` -> Monthly mean outputs from pre-Industrial simulation representing an average over 80 model years
- `PMIP4_Ensemble.nc` -> PMIP4 Ensemble average monthly mean surface air temperature (tas) and precipitation (pr) for mid-Holocene and pre-Indsutrial simulations

`Data/Models/Timeseries`:

The following complete timeseries are provided for both mid-Holocene (MH) and pre-Industrial simulations:

- `pr` -> Total precipitaiton [kg m-2 s-1]
- `q_plev` -> Specific humidity on atmospheric pressure levels [kg kg-1]
- `u_wind` -> U (eastward) component of wind speed on atmopsheric pressure levels [m s-1]
- `unfrozen_smc` -> Unfrozen soil moisture as a fraction of saturation on soil levels [\1]
- `v_wind` -> V (northward) component of wind speed on atmopsheric pressure levels [m s-1]
