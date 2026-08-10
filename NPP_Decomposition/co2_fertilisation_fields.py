# ============================================================================
#  Mechanistic CO2 fertilisation on NPP: the JULES leaf biochemistry evaluated
#  at two CO2 concentrations on a FIXED climate, giving the per-PFT fractional
#  NPP change attributable to the CO2 change alone.
#
#  This is a library module — npp_decomposition.py imports
#  `frac_on_climate`, `R_GROW` and `ZERODEGC`, supplies the PI and MH climates
#  and CO2 values, and averages the two endpoints itself.
#
#  METHOD
#      NPP = (1 - r_grow) * (Ag - R_c),   R_c = f_R * Ag(ref),
#      f_R = 1 - CUE/(1 - r_grow),        CUE = pot_npp / gpp.
#
#  Ag is the JULES Collatz-C3/C4 + Jacobs smoothed co-limited leaf gross rate.
#  R_c is CO2-independent: it is fixed at the reference CO2 but recomputed on
#  whichever climate is passed in, so only climate and CO2 move between calls.
#  The caller holds f_R at its reference (PI) value on both climates.
#
#  Because the N-limitation factor f_N cancels in a log-change, calibrating f_R
#  on POTENTIAL NPP is what correctly attributes the CO2 effect to the actual
#  (npp_nlimit) NPP change.
# ============================================================================

import numpy as np
import xarray as xr

from jules_qsat import qsat 

# ----------------------------------------------------------------------------
# 1. PHYSICAL CONSTANTS 
# ----------------------------------------------------------------------------

EPO2   = 1.106           # ccarbon_mod.F90:28   M_o2 /M_air
O2     = 0.23            # jules_surface_mod.F90:156  O2 mass mixing ratio
FWE_C3 = 0.5             # jules_surface_mod.F90:270
FWE_C4 = 20000.0         # jules_surface_mod.F90:271
BETA1  = 0.83            # jules_surface_mod.F90:266  co-limit carb & lite
BETA2  = 0.93            # jules_surface_mod.F90:267  co-limit (carb,lite) & export
CONPAR = 2.19e5          # sf_stom:250  W -> mol photons s-1  (J/mol)
ZERODEGC = 273.15

# ----------------------------------------------------------------------------
# 2. PER-PFT PARAMETERS
# ----------------------------------------------------------------------------

F0     = np.array([0.875, 0.875, 0.892, 0.875, 0.875, 0.931, 0.931, 0.931, 0.8, 0.8, 0.8, 0.875, 0.875])
DQCRIT = np.array([0.09, 0.09, 0.09, 0.041, 0.06, 0.051, 0.051, 0.051, 0.075, 0.075, 0.075, 0.03, 0.04])
ALPHA  = np.array([0.064, 0.064, 0.048, 0.08, 0.064, 0.048, 0.048, 0.048, 0.04, 0.04, 0.04, 0.064, 0.048])
C3     = np.array([1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1])
VSL    = np.array([31.22, 19.22, 28.4, 23.79, 18.15, 40.96, 40.96, 49.96, 20.48, 20.48, 20.48, 23.15, 23.15])
VINT   = np.array([7.56, 7.21, 3.9, 6.32, 6.32, 6.42, 6.42, 6.42, 0, 0, 0, 14.71, 14.71])
OMEGA  = np.array([0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.17, 0.17, 0.17, 0.15, 0.15])
Q10_LEAF = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
TUPP   = np.array([43, 43, 43, 26, 32, 32, 32, 32, 45, 45, 45, 40, 36])
TLOW   = np.array([10, 13, 13, -10, 0, 10, 10, 10, 13, 13, 13, 10, 0])
NMASS  = np.array([0.021, 0.017, 0.0144, 0.0186, 0.0115, 0.024, 0.024, 0.024, 0.0113, 0.0113, 0.0113, 0.0218, 0.0136])
LMA    = np.array([0.0823, 0.1039, 0.1403, 0.1006, 0.2263, 0.0495, 0.0495, 0.0495, 0.137, 0.137, 0.137, 0.0709, 0.1515])
R_GROW = np.array([0.25] * 13) 

# ----------------------------------------------------------------------------
# 3. LEAF BIOCHEMISTRY 
# ----------------------------------------------------------------------------

def temperature_terms(tstar, pstar):
    """Collatz C3 temperature terms (sf_stom:826-839)."""
    power = 0.1 * (tstar - ZERODEGC - 25.0)
    oa = O2 / EPO2 * pstar
    ccp = 0.5 * oa / (2600.0 * 0.57**power)                    # Pa (C4 -> set 0)
    K = (30.0 * 2.1**power) * (1.0 + oa / (30000.0 * 1.2**power))
    return ccp, K, power

def smoothmin(a, b, beta):
    """Smaller root of beta*s^2 - (a+b)s + a*b = 0  (leaf_jls_mod:170-182)."""
    return (a + b) / (2.0 * beta) - np.sqrt((a + b)**2 / (4.0 * beta**2) - a * b / beta)

def leaf_gross_field(ppm, ft, tstar, pstar, q1, ipar):
    """Gross leaf photosynthesis Ag (mol CO2 m-2 leaf s-1) at CO2 = `ppm` on the
    given (time, lat, lon) climate for PFT `ft`.  NaN where stomata are closed
    (dq >= dqcrit) or dark (ipar <= 0)."""
    ca = ppm * 1e-6 * pstar                                     # Pa
    ccp, K, power = temperature_terms(tstar, pstar)
    if C3[ft] == 0:
        ccp = ccp * 0.0                                         # C4: no photorespiration

    nleaf = NMASS[ft] * LMA[ft] * 1000.0                        # gN/m2 leaf
    tc = tstar - ZERODEGC
    denom = (1 + np.exp(0.3 * (tc - TUPP[ft]))) * (1 + np.exp(0.3 * (TLOW[ft] - tc)))
    vcmax = (VSL[ft] * nleaf + VINT[ft]) * 1e-6 * Q10_LEAF[ft] ** power / denom
    acr = (1.0 - OMEGA[ft]) * ipar / CONPAR

    dq = np.maximum(0.0, qsat(tstar, pstar) - q1)               # humidity deficit (sf_stom:707)
    eta = F0[ft] * (1.0 - dq / DQCRIT[ft])
    ci = eta * ca + (1.0 - eta) * ccp

    if C3[ft] == 1:
        wcarb = vcmax * (ci - ccp) / (ci + K)
        wlite = ALPHA[ft] * acr * (ci - ccp) / (ci + 2.0 * ccp)
        wexpt = FWE_C3 * vcmax + 0.0 * ca
    else:
        wcarb = vcmax + 0.0 * ca
        wlite = ALPHA[ft] * acr + 0.0 * ca
        wexpt = FWE_C4 * vcmax * ci / pstar
    wlite = np.maximum(wlite, np.finfo(float).tiny)

    ag = smoothmin(smoothmin(wcarb, wlite, BETA1), wexpt, BETA2)
    return xr.where((dq < DQCRIT[ft]) & (ipar > 0.0), ag, np.nan)

# ----------------------------------------------------------------------------
# 4. FRACTIONAL CHANGE CALCULATION
# ----------------------------------------------------------------------------

def frac_on_climate(ft, co2, ref_ppm, f_R, tstar, pstar, q1, ipar):
    """Fractional NPP change from the ref -> co2 CO2 change on a FIXED climate.
    R_c = f_R * Ag(ref) is CO2-independent (recomputed here at this climate)."""
    ag_ref = leaf_gross_field(ref_ppm, ft, tstar, pstar, q1, ipar)
    rc     = f_R * ag_ref
    npp_ref = (1.0 - R_GROW[ft]) * (ag_ref - rc)
    ag_c   = leaf_gross_field(co2, ft, tstar, pstar, q1, ipar)
    npp_c  = (1.0 - R_GROW[ft]) * (ag_c - rc)
    return (npp_c / npp_ref - 1.0).where(npp_ref > 0.0)
