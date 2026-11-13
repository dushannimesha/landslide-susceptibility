"""Step 5: run the calibrated Green-Ampt + infinite-slope model for every valid cell over 2020.

Outputs:
  outputs/susceptibility_points.csv  one row per cell with coordinates, the 4 parameters, FS, p, class
  outputs/daily_p.npy, daily_fs.npy  [N, T] float16 series (for date maps)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ls import load_config, p  # noqa: E402
from ls.io import load_processed, save_json  # noqa: E402
from ls.model import run_physics  # noqa: E402


def classify(pv, edges, names):
    return pd.Categorical(np.array(names)[np.digitize(pv, edges)], categories=names, ordered=True)


def main():
    cfg = load_config()
    phys = cfg["physics"]
    out = p(cfg["paths"]["outputs"])
    out.mkdir(parents=True, exist_ok=True)
    static, dyn, dates = load_processed(cfg)
    print("physics constants:", {k: phys[k] for k in ("soil_depth_m", "drainage_tau_days", "psi_scale", "k_sigmoid")})

    r = run_physics(static, dyn, phys, return_series=True)
    lab = pd.read_parquet(p(cfg["paths"]["processed"]) / "labels.parquet").set_index("cell_id")

    pts = static[["cell_id", "Latitude", "Longitude", "row", "col", "Elevation_m", "Slope_deg", "TWI",
                  "LULC_Class", "Lithology_Class", "Lithology_Code", "Soil_Sand_pct", "Soil_Clay_pct",
                  "Ksat_ms", "Cohesion_kPa", "Porosity_n", "Friction_Angle_deg"]].copy()
    pts["psi_f_m"] = r["psi_f_m"]
    pts["FS_dry"] = r["fs_dry"]
    pts["FS_min_2020"] = r["fs_min"]
    pts["p_max_2020"] = r["p_max"]
    pts["date_p_max"] = dates[r["t_pmax"]].strftime("%Y-%m-%d")
    pts["wetting_front_max_m"] = r["zw_max"]
    pts["perched_water_max_m"] = r["hw_max"]
    pts["rain_2020_mm"] = np.nansum(dyn["rain"], axis=1)
    pts["max_rain_24h_mm"] = np.nanmax(dyn["rain"], axis=1)
    pts["absolute_class"] = classify(pts.p_max_2020, cfg["classes"]["edges"], cfg["classes"]["names"])
    # relative susceptibility index (0-1): percentile rank of -FS_min (lowest FS = most susceptible)
    pts["susceptibility_index"] = (-pts.FS_min_2020).rank(pct=True)
    pts["susceptibility_class"] = classify(pts.susceptibility_index * 100, cfg["classes"]["relative_percentiles"],
                                           cfg["classes"]["names"])
    pts["glc_label"] = lab.loc[pts.cell_id, "label"].to_numpy()
    pts["dist_to_glc_event_m"] = lab.loc[pts.cell_id, "dist_to_event_m"].to_numpy()
    pts.to_csv(out / "susceptibility_points.csv", index=False, float_format="%.6g")

    np.save(out / "daily_p.npy", r["p"].astype(np.float16))
    np.save(out / "daily_fs.npy", np.clip(r["fs"], 0, 60).astype(np.float16))
    np.save(out / "dates.npy", dates.strftime("%Y%m%d").astype(int).to_numpy())

    # physical-consistency diagnostics
    fs = r["fs"]
    rain = np.nan_to_num(dyn["rain"])
    dry = rain < 1.0
    wet_rank = np.argsort(-rain.mean(axis=0))[:10]
    diag = {
        "n_cells": int(len(pts)),
        "relative_class_counts": pts.susceptibility_class.value_counts().reindex(cfg["classes"]["names"]).to_dict(),
        "absolute_class_counts": pts.absolute_class.value_counts().reindex(cfg["classes"]["names"]).to_dict(),
        "FS_min_quantiles": pts.FS_min_2020.quantile([0, .05, .5, .95, 1]).round(3).to_dict(),
        "frac_cells_FS_lt_1_any_day": float((r["fs_min"] < 1).mean()),
        "frac_cells_FS_lt_1_dry_start": float((r["fs_dry"] < 1).mean()),
        "median_FS_drop_dry_to_min": float(np.median(r["fs_dry"] - r["fs_min"])),
        "corr_FS_vs_AntRain30_daily_mean": float(np.corrcoef(fs.mean(axis=0), dyn["AntRain_30d_mm"].mean(axis=0))[0, 1]),
        "mean_FS_on_dry_cell_days": float(fs[dry].mean()),
        "mean_FS_on_10_wettest_days": float(fs[:, wet_rank].mean()),
        "params_in_physical_range": float((pts.Ksat_ms.between(1e-8, 1e-3) & pts.Cohesion_kPa.between(0, 40)
                                           & pts.Porosity_n.between(0.2, 0.7)
                                           & pts.Friction_Angle_deg.between(15, 45)).mean()),
    }
    save_json(diag, p(cfg["paths"]["reports"]) / "physics_diagnostics.json")
    print(diag)
    print(pts.sort_values("p_max_2020", ascending=False).head(10)[
        ["Latitude", "Longitude", "Slope_deg", "Lithology_Code", "FS_min_2020", "p_max_2020", "date_p_max"]].to_string())


if __name__ == "__main__":
    main()
