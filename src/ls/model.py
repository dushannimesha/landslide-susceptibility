"""Glue: run Green-Ampt + infinite slope for every cell over the 2020 daily series."""
import numpy as np

from .physics import factor_of_safety, failure_probability, green_ampt_series, rawls_psi_f


def run_physics(static, dyn, phys, cell_idx=None, return_series=False):
    """static: DataFrame with Ksat_ms, Cohesion_kPa, Porosity_n, Friction_Angle_deg, Slope_deg,
    Soil_Sand_pct, Soil_Clay_pct, Bulk_Density_gcm3.  dyn: dict with 'rain' and soil moisture.
    phys: physics config dict (soil_depth_m, drainage_tau_days, psi_scale, k_sigmoid, ...).
    """
    s = static if cell_idx is None else static.iloc[cell_idx]
    sel = (lambda a: a) if cell_idx is None else (lambda a: a[cell_idx])
    rain = sel(dyn["rain"])
    theta = sel(dyn[phys.get("soil_moisture_col", "SoilMoist_7_28cm")])

    slope = np.radians(s["Slope_deg"].to_numpy())
    n = s["Porosity_n"].to_numpy()
    psi = rawls_psi_f(s["Soil_Sand_pct"].to_numpy(), s["Soil_Clay_pct"].to_numpy(), n)
    psi = psi * phys.get("psi_scale", 1.0)
    H = phys["soil_depth_m"]

    ga = green_ampt_series(rain, theta, s["Ksat_ms"].to_numpy(), n, psi, slope,
                           soil_depth_m=H, tau_days=phys["drainage_tau_days"],
                           min_dtheta=phys.get("min_delta_theta", 0.02))
    fs = factor_of_safety(s["Cohesion_kPa"].to_numpy(), s["Friction_Angle_deg"].to_numpy(), slope,
                          s["Bulk_Density_gcm3"].to_numpy(), n, theta, ga["z_w"], ga["h_w"], psi,
                          soil_depth_m=H, gamma_w=phys.get("gamma_w", 9.81)).astype(np.float32)
    prob = failure_probability(fs, phys["k_sigmoid"]).astype(np.float32)

    out = {
        "fs_min": fs.min(axis=1),
        "fs_dry": fs[:, 0],
        "p_max": prob.max(axis=1),
        "t_pmax": prob.argmax(axis=1),
        "zw_max": ga["z_w"].max(axis=1),
        "hw_max": ga["h_w"].max(axis=1),
        "psi_f_m": psi,
    }
    if return_series:
        out.update(fs=fs, p=prob, **ga)
    return out
