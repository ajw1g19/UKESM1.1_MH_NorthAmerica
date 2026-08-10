# ============================================================================
#  npp_decomposition_mh_climstruct.py  —  MH vs PI change in
#  POTENTIAL NPP, split in log space into three terms that sum exactly to
#  ln(pot_M/pot_P):
#
#      moisture   = ln(beta_M / beta_P)
#      co2        = ln(1 + fco2)
#      climstruct = ln(pot_M/pot_P) - moisture - co2
#
#  Direct climate and vegetation structure are reported as ONE term: 
#  separating them needs a fixed-vegetation run
#
#  CO2 uses the mechanistic framework from co2_fertilisation_fields.py: 
#  the JULES Collatz/Jacobs leaf rate is evaluated at PI and MH CO2
#  on BOTH climates, and the two log-ratios averaged (geometric mean of 1+frac)
#  so the climate x CO2 interaction is split symmetrically.
#
#  The three terms decompose POTENTIAL NPP. Nitrogen limitation is recovered at
#  region level as the gap between the aggregate changes in npp_nlimit_pft and
#  pot_npp_pft — pure model output, independent of any component.
#
#  OUTPUT: NPP_decomposition.nc   fields on (time, PFTs, lat, lon)
#          NPP_decomposition.csv  Contributions averaged over the Great Plains
# ============================================================================

import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

import sys
sys.path.insert(0, "../bin")
import um_utils as um
import um_utils_plotting as umplt

from co2_fertilisation_fields import frac_on_climate, R_GROW, ZERODEGC

# ----------------------------------------------------------------------------
# 1. CONFIG
# ----------------------------------------------------------------------------

MEANS_DIR = Path("../Data/Models/Monthly_Means")
FILES = {
    "M": "MH_monthly_means.nc",    # MH 
    "P": "PI_monthly_means.nc",    # PI 
}

CO2_PPM = {"MH": 264.4, "PI": 284.4}

GP_PFTS = ["BrDe", "BrEvTe", "NeDe", "NeEv", "C3Gr", "C4Gr", "ShEv"]

OUT_NC  = "Data/Models/NPP_decomposition.nc"
OUT_CSV = "Data/Models/NPP_decomposition.csv"

# ----------------------------------------------------------------------------
# 2. LOAD
# ----------------------------------------------------------------------------

M = um.importUMData(MEANS_DIR / FILES["M"])
P = um.importUMData(MEANS_DIR / FILES["P"])

# ----------------------------------------------------------------------------
# 3. CO2 FERTILISATION
# ----------------------------------------------------------------------------

def climate(ds):
    """Gridbox climate driving the leaf biochemistry. ts is archived in Celsius."""
    return dict(tstar=ds["ts"] + ZERODEGC, pstar=ds["surfpress"],
                q1=ds["q"], ipar=0.45 * ds["down_ssw"])

def co2_fertilisation():
    """Per-PFT fractional NPP change from the PI -> MH CO2 drop, as the geometric
    mean of the effect evaluated on the PI climate and on the MH climate."""
    clim_pi, clim_mh = climate(P), climate(M)

    cue = (P["pot_npp_pft"] / P["gpp_pft"]).where(P["gpp_pft"] > 0.0)

    fracs = []
    for ft in range(P.sizes["PFTs"]):
        f_R = 1.0 - cue.isel(PFTs=ft, drop=True) / (1.0 - R_GROW[ft])
        a = frac_on_climate(ft, CO2_PPM["MH"], CO2_PPM["PI"], f_R, **clim_pi)
        b = frac_on_climate(ft, CO2_PPM["MH"], CO2_PPM["PI"], f_R, **clim_mh)
        fracs.append(np.sqrt((1.0 + a) * (1.0 + b)) - 1.0)

    return xr.concat(fracs, dim="PFTs").assign_coords(PFTs=P["PFTs"])

fco2 = co2_fertilisation()


# ----------------------------------------------------------------------------
# 4. COMPONENTS
# ----------------------------------------------------------------------------

# Every log in the identity must be finite; marginal cells are dropped rather
# than zero-filled, which would bias the regional aggregation.

valid = ((M["pot_npp_pft"] > 0) & (P["pot_npp_pft"] > 0)
         & (M["beta_factor"] > 0) & (P["beta_factor"] > 0)
         & np.isfinite(fco2))

TOTAL      = np.log(M["pot_npp_pft"] / P["pot_npp_pft"])
MOISTURE   = np.log(M["beta_factor"] / P["beta_factor"])
CO2        = np.log(1.0 + fco2)
CLIMSTRUCT = TOTAL - MOISTURE - CO2                     # residual, by definition
SUM        = MOISTURE + CO2 + CLIMSTRUCT                # == ln(pot_M / pot_P)

FIELDS = [
    ("moisture_cont",   MOISTURE,   "Log-fractional contribution of moisture (beta_factor) to potential NPP"),
    ("co2_cont",        CO2,        "Log-fractional contribution of CO2 fertilisation to potential NPP"),
    ("climstruct_cont", CLIMSTRUCT, "Log-fractional contribution of direct climate and vegetation structure (residual)"),
    ("combined",        SUM,        "Sum of all three components = ln(pot_npp_MH / pot_npp_PI)"),
]

out = xr.Dataset({name: da.where(valid).rename(name) for name, da, _ in FIELDS})
for name, _, long_name in FIELDS:
    out[name].attrs.update(long_name=long_name, units="1")
out.attrs["description"] = ("MH - PI potential-NPP decomposition in log space")
out.to_netcdf(OUT_NC)


# ----------------------------------------------------------------------------
# 5. AVERAGE OVER GREAT PLAINS
# Per cell the three terms sum to ln(pot_M/pot_P). Weighting each cell by
# area * logmean(pot_P, pot_M) makes the weighted terms sum exactly to the log
# change of the regional total potential NPP, so the carbon
# budget closes at region level rather than only per cell
# ----------------------------------------------------------------------------

AREA    = np.cos(np.deg2rad(M.latitude))
gp_mask = umplt.umRegionMask("GP", M.longitude, M.latitude)

def log_mean(a, b):
    """L(a,b) = (b-a)/(ln b - ln a); -> a as b -> a. The unique cell weight for
    which flux-weighted log-changes sum to the change in the regional total."""
    a = a.where(a > 0.0); b = b.where(b > 0.0)
    dln = np.log(b) - np.log(a)
    return xr.where(dln == 0.0, a, (b - a) / xr.where(dln == 0.0, 1.0, dln))

potM, potP = M["pot_npp_pft"], P["pot_npp_pft"]
W = (AREA * log_mean(potP, potM)).where(valid)

def wsum(field, sel):
    """Area-weighted regional sum over time/lat/lon for one PFT."""
    return float((AREA * field.isel(sel)).where(valid.isel(sel) & gp_mask).sum())

def wmean(field, sel):
    """Area-weighted regional mean over time/lat/lon for one PFT."""
    return float(field.isel(sel).where(valid.isel(sel) & gp_mask).weighted(AREA).mean())

def tile_frac(ds):
    """Surface fraction on the 13 PFT tiles, sharing the PFTs coordinate."""
    return (ds["surfacefrac"].isel(surface_tiles=slice(0, 13))
            .rename({"surface_tiles": "PFTs"}).assign_coords(PFTs=ds["PFTs"]))

COMPONENTS = [("moisture", MOISTURE), ("co2", CO2),
              ("climstruct", CLIMSTRUCT)]

fracM, fracP = tile_frac(M), tile_frac(P)
rows = {}

for i, pft in enumerate(M.PFTs.values):
    if pft not in GP_PFTS:
        continue
    sel = dict(PFTs=i)
    w = W.isel(sel).where(gp_mask)

    # Regional totals give the exact log-change and the LMDI normaliser.
    A_agg, B_agg = wsum(potP, sel), wsum(potM, sel)
    total_log = np.log(B_agg / A_agg)
    L_AB      = (B_agg - A_agg) / (np.log(B_agg) - np.log(A_agg))
    total_pct = (np.exp(total_log) - 1) * 100

    print(f"==={pft}===")
    S, contribs = 0.0, {}
    for name, comp in COMPONENTS:
        c_log = float((w * comp.isel(sel)).sum()) / L_AB
        S += c_log
        contribs[name] = (c_log / total_log) * total_pct
        print(f"{name}: {contribs[name]:.2f} % ({c_log / total_log * 100:.1f} % share)")

    dpct = (np.exp(S) - 1) * 100
    print(f"  potential dNPP [decomposition Sum comps] = {dpct:.2f} %")
    print(f"  potential dNPP [model pot_npp_pft]       = {total_pct:.2f} %   "
          f"(closure gap = {dpct - total_pct:+.2f} pp)")

    # N-limited change straight from the model, aggregated over the SAME valid
    # and GP cells, so the gap against the potential change is the N-limitation.
    A_nl, B_nl = wsum(P["npp_nlimit_pft"], sel), wsum(M["npp_nlimit_pft"], sel)
    actual_pct = (B_nl / A_nl - 1) * 100
    nitrogen   = actual_pct - total_pct
    print(f"  actual dNPP [model npp_nlimit_pft]       = {actual_pct:.2f} %")
    print(f"  N-limitation (model actual - model pot)  = {nitrogen:+.2f} pp")
    print(" ")

    rows[pft] = {
        "moisture":      contribs["moisture"],
        "co2":           contribs["co2"],
        "climstruct":    contribs["climstruct"],
        "nitrogen":      nitrogen,
        "potential_npp": total_pct,
        "actual_npp":    actual_pct,
        "surfacefrac":   (wmean(fracM, sel) - wmean(fracP, sel)) * 100,
    }

df = pd.DataFrame.from_dict(rows, orient="index",
                            columns=["moisture", "co2", "climstruct",
                                     "nitrogen", "potential_npp",
                                     "actual_npp", "surfacefrac"])
df.index.name = "PFT"
df.to_csv(OUT_CSV)
print(f"saved -> {OUT_NC} + {OUT_CSV}")
