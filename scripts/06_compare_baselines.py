"""Step 6: benchmark the physics model against data-driven models with spatial-block CV.

Models (all scored on the same held-out spatial blocks against the GLC inventory):
  physics      : min FS(2020) from calibrated Green-Ampt + infinite slope (no learning)
  slope_only   : terrain slope (trivial baseline)
  random_forest, xgboost : static + rainfall summary features
Writes reports/metrics.json and reports/metrics.md
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ls import load_config, p  # noqa: E402
from ls.io import load_processed, save_json  # noqa: E402
from ls.metrics import block_folds, scores, spatial_blocks  # noqa: E402

FEATS = ["Elevation_m", "Slope_deg", "Eastness", "Northness", "Horiz_Curvature", "Vert_Curvature", "TPI",
         "TRI", "TWI", "SPI", "Soil_Sand_pct", "Soil_Silt_pct", "Soil_Clay_pct", "Bulk_Density_gcm3",
         "LULC_Class", "Lithology_Class", "Ksat_ms", "Cohesion_kPa", "Porosity_n", "Friction_Angle_deg"]


def main():
    cfg = load_config()
    static, dyn, _ = load_processed(cfg)
    pts = pd.read_csv(p(cfg["paths"]["outputs"]) / "susceptibility_points.csv").set_index("cell_id")
    pts = pts.loc[static.cell_id]
    y = pts.glc_label.to_numpy()
    use = y >= 0

    X = static[FEATS].copy()
    X["rain_2020"] = np.nansum(dyn["rain"], axis=1)
    X["max_rain_24h"] = np.nanmax(dyn["rain"], axis=1)
    X["max_ant30"] = np.nanmax(dyn["AntRain_30d_mm"], axis=1)
    X["sm_mean"] = np.nanmean(dyn[cfg["physics"]["soil_moisture_col"]], axis=1)
    X = X.to_numpy(dtype=np.float32)

    blocks = spatial_blocks(static.Latitude.to_numpy(), static.Longitude.to_numpy(), cfg["calibration"]["n_blocks"])
    folds = block_folds(blocks, y)
    preds = {k: np.full(len(y), np.nan) for k in ("physics_p", "physics_minFS", "slope_only", "random_forest", "xgboost")}
    preds["physics_p"] = pts.p_max_2020.to_numpy()
    preds["physics_minFS"] = -pts.FS_min_2020.to_numpy()   # ranking score
    preds["slope_only"] = static.Slope_deg.to_numpy()

    for f in range(folds.max() + 1):
        tr = use & (folds != f)
        te = use & (folds == f)
        if (y[tr] == 1).sum() == 0 or te.sum() == 0:
            continue
        w = (y[tr] == 0).sum() / max((y[tr] == 1).sum(), 1)
        rf = RandomForestClassifier(n_estimators=400, min_samples_leaf=5, class_weight="balanced_subsample",
                                    n_jobs=-1, random_state=0).fit(X[tr], y[tr])
        preds["random_forest"][te] = rf.predict_proba(X[te])[:, 1]
        xg = XGBClassifier(n_estimators=400, max_depth=4, learning_rate=0.05, subsample=0.8,
                           colsample_bytree=0.8, scale_pos_weight=w, n_jobs=16, eval_metric="auc").fit(X[tr], y[tr])
        preds["xgboost"][te] = xg.predict_proba(X[te])[:, 1]

    res = {}
    for k, s in preds.items():
        m = use & ~np.isnan(s)
        res[k] = {"pooled": scores(y[m], s[m]),
                  "per_fold_auc": [scores(y[m & (folds == f)], s[m & (folds == f)])["auc"]
                                   for f in range(folds.max() + 1)]}
        res[k]["auc_fold_mean"] = float(np.nanmean(res[k]["per_fold_auc"]))
    res["note"] = ("Labels = GLC 2007-2017 event locations (spatial inventory); physics forced by 2020 daily data. "
                   "F1/hit/false-alarm thresholds are only meaningful for probability outputs (physics_p, RF, XGB).")
    save_json(res, p(cfg["paths"]["reports"]) / "metrics.json")

    L = ["# Model comparison (spatial-block CV)", "", "| model | ROC-AUC | PR-AUC | mean fold AUC | hit@0.5 | false alarm@0.5 |",
         "|---|---|---|---|---|---|"]
    for k, v in res.items():
        if k == "note":
            continue
        s = v["pooled"]
        L.append(f"| {k} | {s['auc']:.3f} | {s['pr_auc']:.3f} | {v['auc_fold_mean']:.3f} | "
                 f"{s.get('hit_rate@0.5', float('nan')):.3f} | {s.get('false_alarm@0.5', float('nan')):.3f} |")
    L += ["", f"Positives: {int((y == 1).sum())}, negatives: {int((y == 0).sum())}", "", res["note"]]
    (p(cfg["paths"]["reports"]) / "metrics.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
