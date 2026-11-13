"""Step 1: read the 2.5 GB daily CSV + the 4-parameter grid CSV, join on coordinates and write
compact processed files:
  data/processed/daily.parquet   full table (float32) for ad-hoc analysis
  data/processed/static.parquet  one row per cell: coords, grid row/col, terrain, soil, LULC,
                                 lithology + Ksat_ms, Cohesion_kPa, Porosity_n, Friction_Angle_deg
  data/processed/dynamic.npz     [N_cells, T] arrays: rain, AntRain_*, SoilMoist_* + dates
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.csv as pacsv

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ls import load_config, p  # noqa: E402

STATIC_COLS = ["Elevation_m", "Slope_deg", "Aspect_deg", "Eastness", "Northness", "Horiz_Curvature",
               "Vert_Curvature", "TPI", "TRI", "TWI", "SPI", "Soil_Sand_pct", "Soil_Silt_pct",
               "Soil_Clay_pct", "Bulk_Density_gcm3", "LULC_Class", "Lithology_Class", "Lithology_Code"]
DYN_COLS = ["Rain_24h_mm", "AntRain_3d_mm", "AntRain_7d_mm", "AntRain_30d_mm",
            "SoilMoist_0_7cm", "SoilMoist_7_28cm"]
PARAM_COLS = ["Ksat_ms", "Cohesion_kPa", "Porosity_n", "Friction_Angle_deg"]


def main():
    cfg = load_config()
    t0 = time.time()
    out = p(cfg["paths"]["processed"])
    out.mkdir(parents=True, exist_ok=True)
    r = cfg["coord_round"]

    print("[1/4] reading main CSV with pyarrow ...", flush=True)
    tbl = pacsv.read_csv(p(cfg["paths"]["main_csv"]))
    df = tbl.to_pandas()
    del tbl
    for c in df.columns:
        if df[c].dtype == "float64" and c not in ("Latitude", "Longitude"):
            df[c] = df[c].astype("float32")
    df["Lithology_Code"] = df["Lithology_Code"].astype("category")
    print(f"      {len(df):,} rows x {df.shape[1]} cols  ({time.time()-t0:.0f}s)", flush=True)
    df.to_parquet(out / "daily.parquet", index=False)

    print("[2/4] joining 4 geotechnical parameters ...", flush=True)
    prm = pd.read_csv(p(cfg["paths"]["params_csv"]))
    for d in (df, prm):
        d["lat_r"] = d["Latitude"].round(r)
        d["lon_r"] = d["Longitude"].round(r)
    prm = prm.drop_duplicates(["lat_r", "lon_r"])

    static = (df.drop_duplicates(["lat_r", "lon_r"])[["lat_r", "lon_r", "Latitude", "Longitude"] + STATIC_COLS]
              .merge(prm[["lat_r", "lon_r"] + PARAM_COLS], on=["lat_r", "lon_r"], how="left"))
    unmatched = static[PARAM_COLS].isna().any(axis=1).sum()
    print(f"      cells: {len(static):,}   unmatched to params CSV: {unmatched}", flush=True)
    if unmatched:
        static[static[PARAM_COLS].isna().any(axis=1)][["Latitude", "Longitude"]].to_csv(
            out / "unmatched_cells.csv", index=False)

    # grid row/col from the regular 250 m spacing
    lats = np.sort(static["lat_r"].unique())
    lons = np.sort(static["lon_r"].unique())
    dlat = np.median(np.diff(lats))
    dlon = np.median(np.diff(lons))
    static["row"] = np.rint((static["lat_r"].max() - static["lat_r"]) / dlat).astype(int)
    static["col"] = np.rint((static["lon_r"] - static["lon_r"].min()) / dlon).astype(int)
    static = static.sort_values(["row", "col"]).reset_index(drop=True)
    static.insert(0, "cell_id", np.arange(len(static)))
    static.attrs = {}
    static.to_parquet(out / "static.parquet", index=False)
    print(f"      grid {static.row.max()+1} rows x {static.col.max()+1} cols, "
          f"dlat={dlat:.6f} dlon={dlon:.6f}", flush=True)

    print("[3/4] building [N, T] dynamic cube ...", flush=True)
    dates = np.sort(df["Date"].unique())
    cid = df[["lat_r", "lon_r"]].merge(static[["lat_r", "lon_r", "cell_id"]], on=["lat_r", "lon_r"], how="left")["cell_id"].to_numpy()
    tid = np.searchsorted(dates, df["Date"].to_numpy())
    N, T = len(static), len(dates)
    cube = {}
    for c in DYN_COLS:
        a = np.full((N, T), np.nan, dtype=np.float32)
        a[cid, tid] = df[c].to_numpy(dtype=np.float32)
        cube[c] = a
    dup = len(df) - len(np.unique(cid.astype(np.int64) * T + tid))
    print(f"      N={N} T={T}  duplicate (cell,date) rows: {dup}  "
          f"missing cell-days: {N*T - (len(df) - dup)}", flush=True)
    np.savez_compressed(out / "dynamic.npz", dates=dates.astype(np.int64), rain=cube["Rain_24h_mm"],
                        **{k: v for k, v in cube.items()})
    print(f"[4/4] done in {time.time()-t0:.0f}s -> {out}", flush=True)


if __name__ == "__main__":
    main()
