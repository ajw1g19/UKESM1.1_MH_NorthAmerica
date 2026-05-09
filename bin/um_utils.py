#!/usr/bin/env python3

import numpy as np
import pandas as pd
import xarray as xr
import warnings
from xarray.coding.common import SerializationWarning
from dask.array.core import PerformanceWarning

#---------------------------------------------------
# Suppress frequenctly occurring non-consequential 
# warnings
#---------------------------------------------------

def WarningSuppress():
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in buffer")
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in intersects")
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in intersection")
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in covers")
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in disjoint")
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in union")
    warnings.filterwarnings("ignore", category=UserWarning, message="This figure includes Axes that are not compatible with tight_layout, so results might be incorrect.")
    warnings.filterwarnings("ignore", category=SerializationWarning)
    warnings.filterwarnings("ignore", category=PerformanceWarning)

#---------------------------------------------------
# Import UM data from netCDF file
#---------------------------------------------------

def importUMData(ncfile, coords=False):
    Data = xr.open_dataset(ncfile, engine="netcdf4", decode_times=False)

    c = {}
    for coord in list(Data.coords):
        c[coord] = Data[coord].data

    if len(Data.data_vars) == 1:
        Data = xr.load_dataarray(ncfile, engine="netcdf4", decode_times=False)
    else:
        Data = xr.load_dataset(ncfile, engine="netcdf4", decode_times=False)
    
    if coords:
        return Data, c
    else:
        return Data
    
#------------------------------------------------------------
# PFT-parameters as defined in UM namelists
#------------------------------------------------------------

def pft_params(parameter=None):
    a_ws = np.array([12, 13, 12, 10, 10, 1, 1, 1, 1, 1, 1, 13, 13]) # Ratio
    a_wl = np.array([0.78, 0.845, 0.78, 0.8, 0.65, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.13, 0.13]) # kgC m-2
    b_wl = np.array([1.667, 1.667, 1.667, 1.667, 1.667, 1.667, 1.667, 1.667, 1.667, 1.667, 1.667, 1.667, 1.667]) # Allometric Exponent
    eta_sl = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]) # kgC m-1 m-2 leaf

    dgl_dt = np.array([9, 9, 9, 9, 9, 0, 0, 0, 0, 0, 0, 9, 9]) # K-1
    tleaf = np.array([280, 278.15, 233.15, 278.15, 233.15, 278.15, 278.15, 278.15, 278.15, 278.15, 278.15, 280, 233.15]) # K
    g_leaf_0 = np.array([0.25, 0.25, 0.5, 0.25, 0.25, 3, 3, 3, 3, 3, 3, 0.25, 0.66]) # yr-1
    g_grow = np.array([20, 15, 15, 20, 15, 20, 20, 20, 20, 20, 20, 30, 15]) # yr-1
    g_root = np.array([0.15, 0.25, 0.25, 0.15, 0.15, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.15, 0.15]) # yr-1
    g_wood = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.05, 0.05]) # yr-1

    params = {
        "a_ws": a_ws, "a_wl": a_wl, "b_wl": b_wl, "eta_sl": eta_sl, 
        "dgl_dt": dgl_dt, "tleaf": tleaf,
        "g_leaf_0": g_leaf_0, "g_grow": g_grow, "g_root": g_root, "g_wood": g_wood
    }

    if parameter is None:
        return {k: v for k, v in params.items()}

    if isinstance(parameter, str):
        if parameter not in params:
            raise KeyError(f"Unknown parameter '{parameter}'. Available: {list(params.keys())}")
        return params[parameter]

    if isinstance(parameter, (list, tuple, set)):
        missing = [p for p in parameter if p not in params]
        if missing:
            raise KeyError(f"Unknown parameters requested: {missing}")
        return {p: params[p] for p in parameter}

    raise TypeError("parameter must be None, str, or list/tuple/set of str")

#------------------------------------------------------------
# Import reconstructions from Dawson REVEALS dataset
#------------------------------------------------------------

def importRevealsYrs(file, years):
    """
    Read REVEALS CSV and return a tuple of DataFrames for the requested years.

    - years: int or iterable of ints
    - Only (x,y) locations that contain all requested years are retained.
    - Returns: tuple of DataFrames in the same order as `years`.
    """
    df = pd.read_csv(file)
    if isinstance(years, (int, np.integer)):
        years = [int(years)]
    else:
        years = [int(y) for y in years]

    df["ages"] = df["ages"].astype(int)
    groups = df.groupby(["x", "y"])

    req_set = set(years)
    mask = groups["ages"].transform(lambda s: set(s) >= req_set)
    df_both = df[mask].copy()

    out = tuple(df_both[df_both["ages"] == y].reset_index(drop=True) for y in years)

    return out

#------------------------------------------------------------
# Bootstrapping functions for timeseries outputs
#------------------------------------------------------------

def bootstrap_years_anom(years_old, years_new, n_boot=2000, rng=42):
    rng = np.random.default_rng(rng)
    n_years = years_old.shape[0]
    boot_mon_mean = np.empty((n_boot, 12))
    for i in range(n_boot):
        idx = rng.integers(0, n_years, n_years)
        boot_mon_mean[i] = years_new[idx].mean(dim="year") - years_old[idx].mean(dim="year")

    mon_mean = boot_mon_mean.mean(axis=0)
    ci_upper = np.percentile(boot_mon_mean, 97.5, axis=0)
    ci_lower = np.percentile(boot_mon_mean, 2.5, axis=0)

    return mon_mean, ci_upper, ci_lower


def bootstrap_years(years, n_boot=2000, rng=42):
    rng = np.random.default_rng(rng)
    n_years = years.shape[0]
    boot_mon_mean = np.empty((n_boot, 12))
    for i in range(n_boot):
        idx = rng.integers(0, n_years, n_years)
        boot_mon_mean[i] = years[idx].mean(dim="year")

    mon_mean = boot_mon_mean.mean(axis=0)
    ci_upper = np.percentile(boot_mon_mean, 97.5, axis=0)
    ci_lower = np.percentile(boot_mon_mean, 2.5, axis=0)

    return mon_mean, ci_upper, ci_lower