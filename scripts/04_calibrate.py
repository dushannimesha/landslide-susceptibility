"""Step 4: calibrate the few global physics constants (soil depth H, drainage tau, psi scale)
by spatial-block cross-validated ROC-AUC of max-p against the GLC inventory.
The 4 geotechnical parameters (Ksat, c', n, phi') are NOT fitted - they come from the
parameter CSV. Writes reports/calibration.csv and configs/calibrated.yaml.
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ls import ROOT, load_config, p  # noqa: E402
from ls.io import load_processed  # noqa: E402
from ls.metrics import block_folds, scores, spatial_blocks  # noqa: E402
from ls.model import run_physics  # noqa: E402


def main():
    cfg = load_config(use_calibrated=False)
    cal = cfg["calibration"]
    static, dyn, _ = load_processed(cfg)
    lab = pd.read_parquet(p(cfg["paths"]["processed"]) / "labels.parquet").set_index("cell_id")
    y = lab.loc[static.cell_id, "label"].to_numpy()
    use = y >= 0
    blocks = spatial_blocks(static.Latitude.to_numpy(), static.Longitude.to_numpy(), cal["n_blocks"])
    folds = block_folds(blocks, y)
    print(f"cells={len(y)} positives={(y == 1).sum()} folds={np.bincount(folds[use])}", flush=True)

    rows = []
    grid = list(itertools.product(cal["soil_depth_m"], cal["drainage_tau_days"], cal["psi_scale"]))
    for H, tau, ps in grid:
        phys = {**cfg["physics"], "soil_depth_m": H, "drainage_tau_days": tau, "psi_scale": ps}
        r = run_physics(static, dyn, phys)
        # rank by min FS (k-independent); p_max is monotone in fs_min
        s = -r["fs_min"]
        row = {"soil_depth_m": H, "drainage_tau_days": tau, "psi_scale": ps,
               "auc_all": scores(y[use], s[use])["auc"],
               "frac_fs_lt1": float((r["fs_min"] < 1).mean()),
               "fs_min_median": float(np.median(r["fs_min"]))}
        for f in range(folds.max() + 1):
            m = use & (folds == f)
            row[f"auc_fold{f}"] = scores(y[m], s[m])["auc"]
        rows.append(row)
        print({k: round(v, 3) if isinstance(v, float) else v for k, v in row.items()}, flush=True)
    res = pd.DataFrame(rows)
    fold_cols = [c for c in res if c.startswith("auc_fold")]
    res["auc_cv_mean"] = res[fold_cols].mean(axis=1)
    # cross-validated selection: for each fold choose params on the other folds, score on this one
    cv_scores = []
    for c in fold_cols:
        others = [o for o in fold_cols if o != c]
        best = res[others].mean(axis=1).idxmax()
        cv_scores.append(res.loc[best, c])
    res = res.sort_values("auc_cv_mean", ascending=False)
    rep = p(cfg["paths"]["reports"])
    rep.mkdir(parents=True, exist_ok=True)
    res.to_csv(rep / "calibration.csv", index=False)

    best = res.iloc[0]
    informative = bool(best.auc_cv_mean >= cal.get("min_informative_auc", 0.55))
    result = {"auc_cv_mean": float(best.auc_cv_mean), "auc_all": float(best.auc_all),
              "nested_cv_auc": float(np.nanmean(cv_scores)), "labels_informative": informative,
              "best_grid_point": {k: float(best[k]) for k in ("soil_depth_m", "drainage_tau_days", "psi_scale")}}
    if informative:
        phys = {**cfg["physics"], **result["best_grid_point"]}
    else:
        # labels do not discriminate (AUC ~<= 0.5) -> calibrating on them would be fitting noise/bias;
        # keep literature defaults from default.yaml
        phys = dict(cfg["physics"])
        result["note"] = ("Best CV AUC below threshold: GLC inventory not informative for the physics "
                          "constants (reporting bias toward urban cut slopes). Literature defaults kept.")
    out = {"physics": phys, "calibration_result": result}
    with open(ROOT / "configs" / "calibrated.yaml", "w") as f:
        yaml.safe_dump(out, f, sort_keys=False)
    print("\nBEST:", out)


if __name__ == "__main__":
    main()
