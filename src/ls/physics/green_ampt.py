"""Slope-corrected 1D Green-Ampt infiltration, vectorised over cells, daily time step.

Infiltration capacity (proposal eq.):  f = Ksat * (z_w cos(beta) + psi_f) / z_w
Wetting-front depth:                   z_w = F / dtheta,  dtheta = n - theta_i
Water accumulated in the profile F recedes (drainage/redistribution) with time constant tau.
"""
import numpy as np


def rawls_psi_f(sand_pct, clay_pct, porosity):
    """Wetting-front suction head [m] from Rawls & Brakensiek (1983) regression (S, C in %)."""
    S = np.clip(sand_pct, 5, 70)
    C = np.clip(clay_pct, 5, 60)
    n = np.clip(porosity, 0.25, 0.65)
    ln_psi_cm = (6.5309 - 7.32561 * n + 0.001583 * C**2 + 3.809479 * n**2
                 + 0.000344 * S * C - 0.049837 * S * n + 0.001608 * S**2 * n**2
                 + 0.001602 * C**2 * n**2 - 0.0000136 * S**2 * C
                 - 0.003479 * C**2 * n - 0.000799 * S**2 * n)
    return np.clip(np.exp(ln_psi_cm), 1.0, 100.0) / 100.0


def green_ampt_series(rain_mm, theta_i, ksat_ms, porosity, psi_f_m, slope_rad,
                      soil_depth_m=2.0, tau_days=5.0, min_dtheta=0.02):
    """Run Green-Ampt over a daily rain series.

    rain_mm  : [N, T] daily rainfall (mm/day)
    theta_i  : [N, T] antecedent volumetric water content (m3/m3), e.g. ERA5-Land layer 2
    ksat_ms, porosity, psi_f_m, slope_rad : [N]
    Returns dict of [N, T] arrays: infil_mm, F_m (stored infiltrated water), z_w (wetting-front
    depth, m), h_w (saturated thickness above failure plane, m), dtheta.
    """
    rain = np.nan_to_num(np.asarray(rain_mm, dtype=np.float64), nan=0.0) / 1000.0  # m/day
    N, T = rain.shape
    ksat = np.asarray(ksat_ms, dtype=np.float64) * 86400.0                          # m/day
    n = np.asarray(porosity, dtype=np.float64)
    psi = np.asarray(psi_f_m, dtype=np.float64)
    cosb = np.cos(np.asarray(slope_rad, dtype=np.float64))
    H = float(soil_depth_m)
    decay = np.exp(-1.0 / float(tau_days))

    th = np.nan_to_num(np.asarray(theta_i, dtype=np.float64), nan=np.nanmean(theta_i))
    dtheta = np.maximum(n[:, None] - np.minimum(th, n[:, None]), min_dtheta)

    F = np.zeros(N)
    out = {k: np.zeros((N, T), dtype=np.float32) for k in ("infil_mm", "F_m", "z_w", "h_w")}
    for t in range(T):
        dth = dtheta[:, t]
        z_w = np.maximum(F / dth, 1e-3)
        # infiltration capacity over one day (slope-normal rain = rain*cos(beta))
        f_cap = ksat * (z_w * cosb + psi) / z_w
        rain_n = rain[:, t] * cosb
        infil = np.minimum(rain_n, f_cap)
        F = F * decay + infil
        z_w = F / dth
        # once the front passes H the excess is stored as a perched saturated layer above the
        # (less permeable) failure plane
        h_w = np.clip(z_w - H, 0.0, H)
        out["infil_mm"][:, t] = infil * 1000.0
        out["F_m"][:, t] = F
        out["z_w"][:, t] = np.minimum(z_w, H)
        out["h_w"][:, t] = h_w
    out["dtheta"] = dtheta.astype(np.float32)
    return out
