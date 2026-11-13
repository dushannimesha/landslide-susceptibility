"""Step 2: exploratory data analysis + quality filtering.

Writes reports/eda/eda_report.md (+ PNGs) and data/processed/valid_mask.npy.
Also recomputes antecedent rainfall from Rain_24h (the GEE export leaves AntRain NaN at the start
of the year) and saves it back into dynamic.npz.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ls import load_config, p  # noqa: E402
from ls.io import load_processed  # noqa: E402

LULC_NAMES = {10: "Tree cover", 20: "Shrubland", 30: "Grassland", 40: "Cropland", 50: "Built-up",
              60: "Bare/sparse", 70: "Snow/ice", 80: "Water", 90: "Herb. wetland", 95: "Mangroves", 100: "Moss"}


def grid_image(static, values):
    img = np.full((static.row.max() + 1, static.col.max() + 1), np.nan)
    img[static.row, static.col] = values
    return img


def extent(static):
    return [static.Longitude.min(), static.Longitude.max(), static.Latitude.min(), static.Latitude.max()]


def main():
    cfg = load_config()
    f = cfg["filters"]
    rep = p(cfg["paths"]["reports"]) / "eda"
    rep.mkdir(parents=True, exist_ok=True)
    static, dyn, dates = load_processed(cfg, filtered=False)
    N, T = dyn["rain"].shape
    L = ["# EDA & filtering report", "", f"Cells: **{N:,}**, days: **{T}** "
         f"({dates.min().date()} – {dates.max().date()}), cell-days: {N*T:,}", ""]

    # ---- missing values
    L += ["## Missing values", "", "| variable | missing | % |", "|---|---|---|"]
    for c in static.columns:
        m = static[c].isna().sum()
        if m:
            L.append(f"| {c} (static) | {m:,} | {100*m/N:.2f} |")
    for k, a in dyn.items():
        m = int(np.isnan(a).sum())
        L.append(f"| {k} (daily) | {m:,} | {100*m/a.size:.2f} |")

    # ---- recompute antecedent rainfall from Rain_24h (fills the January gap consistently)
    rain = np.nan_to_num(dyn["rain"], nan=0.0)
    cs = np.cumsum(np.pad(rain, ((0, 0), (1, 0))), axis=1)
    full = dict(np.load(p(cfg["paths"]["processed"]) / "dynamic.npz"))
    for w in (3, 7, 30):
        idx = np.arange(1, T + 1)
        ant = cs[:, idx] - cs[:, np.maximum(idx - w, 0)]
        orig = dyn[f"AntRain_{w}d_mm"]
        ok = ~np.isnan(orig)
        diff = np.nanmedian(np.abs(ant[ok] - orig[ok])) if ok.any() else np.nan
        L.append(f"| AntRain_{w}d recomputed from Rain_24h | median |Δ| vs GEE = {diff:.2f} mm | |")
        full[f"AntRain_{w}d_mm"] = ant.astype(np.float32)
    np.savez_compressed(p(cfg["paths"]["processed"]) / "dynamic.npz", **full)

    # ---- descriptive stats
    L += ["", "## Static variables", "", static.drop(columns=["cell_id", "row", "col", "lat_r", "lon_r"],
          errors="ignore").describe().T.round(4).to_markdown(), ""]
    rain_cell_total = rain.sum(axis=1)
    L += ["## Rainfall 2020", "",
          f"- Annual total per cell: mean {rain_cell_total.mean():.0f} mm, min {rain_cell_total.min():.0f}, max {rain_cell_total.max():.0f}",
          f"- Max daily rain anywhere: {np.nanmax(dyn['rain']):.1f} mm; negative values: {(dyn['rain'] < 0).sum()}",
          f"- Unique rain values per day (spatial resolution check): median {np.median([len(np.unique(rain[:, t])) for t in range(T)]):.0f}", ""]

    # ---- range / consistency checks
    tex = static.Soil_Sand_pct + static.Soil_Silt_pct + static.Soil_Clay_pct
    sm = dyn[cfg["physics"]["soil_moisture_col"]]
    checks = {
        "Slope outside 0–90°": int((~static.Slope_deg.between(0, 90)).sum()),
        "Sand+Silt+Clay not within 95–105 %": int((~tex.between(95, 105)).sum()),
        "Bulk density outside 0.8–2.1 g/cm³": int((~static.Bulk_Density_gcm3.between(0.8, 2.1)).sum()),
        "Soil moisture > porosity (cell-days)": int((sm > static.Porosity_n.to_numpy()[:, None]).sum()),
        "Soil moisture outside 0–1 (cell-days)": int(((sm < 0) | (sm > 1)).sum()),
        "Ksat outside 1e-8–1e-3 m/s": int((~static.Ksat_ms.between(1e-8, 1e-3)).sum()),
        "Duplicate coordinates": int(static.duplicated(["lat_r", "lon_r"]).sum()),
    }
    const = [c for c in static.columns if static[c].nunique() <= 1]
    L += ["## Consistency checks", "", "| check | count |", "|---|---|"]
    L += [f"| {k} | {v:,} |" for k, v in checks.items()]
    L += [f"| Constant columns | {', '.join(const) or 'none'} |", ""]

    L += ["## Categorical layers", "", "| LULC | cells |", "|---|---|"]
    L += [f"| {int(k)} {LULC_NAMES.get(int(k), '')} | {v:,} |" for k, v in static.LULC_Class.value_counts().items()]
    L += ["", "| Lithology | cells |", "|---|---|"]
    L += [f"| {k} | {v:,} |" for k, v in static.Lithology_Code.value_counts().items()]

    # ---- correlation
    num = static[["Elevation_m", "Slope_deg", "TPI", "TRI", "TWI", "SPI", "Soil_Sand_pct", "Soil_Clay_pct",
                  "Bulk_Density_gcm3", "Ksat_ms", "Cohesion_kPa", "Porosity_n", "Friction_Angle_deg"]].copy()
    num["Ksat_ms"] = np.log10(num["Ksat_ms"])
    num["AnnualRain"] = rain_cell_total
    corr = num.corr()
    hi = [(a, b, corr.loc[a, b]) for i, a in enumerate(corr) for b in corr.columns[i + 1:] if abs(corr.loc[a, b]) > 0.8]
    L += ["", "## Highly correlated pairs (|r| > 0.8)", ""] + [f"- {a} ~ {b}: r = {r:.2f}" for a, b, r in hi] + [""]
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr)), corr.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(len(corr)), corr.columns, fontsize=8)
    fig.colorbar(im)
    fig.tight_layout()
    fig.savefig(rep / "correlation.png", dpi=130)
    plt.close(fig)

    # ---- maps of key static layers
    layers = ["Elevation_m", "Slope_deg", "TWI", "Lithology_Class", "LULC_Class", "Ksat_ms",
              "Cohesion_kPa", "Friction_Angle_deg", "Porosity_n"]
    fig, axs = plt.subplots(3, 3, figsize=(14, 13))
    for ax, c in zip(axs.flat, layers):
        v = np.log10(static[c]) if c == "Ksat_ms" else static[c]
        im = ax.imshow(grid_image(static, v), extent=extent(static), cmap="viridis")
        ax.set_title(("log10 " if c == "Ksat_ms" else "") + c)
        ax.set_xlabel("Longitude (°E)")
        ax.set_ylabel("Latitude (°N)")
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(rep / "static_layers.png", dpi=110)
    plt.close(fig)

    # ---- rainfall time series
    fig, ax = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    ax[0].bar(dates, rain.mean(axis=0), color="tab:blue")
    ax[0].set_ylabel("Mean daily rain (mm)")
    ax[1].plot(dates, np.nanmean(dyn["SoilMoist_0_7cm"], axis=0), label="0–7 cm")
    ax[1].plot(dates, np.nanmean(dyn["SoilMoist_7_28cm"], axis=0), label="7–28 cm")
    ax[1].set_ylabel("Soil moisture (m³/m³)")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(rep / "rain_soilmoisture_timeseries.png", dpi=120)
    plt.close(fig)

    # ---- filters
    masks = {
        f"slope < {f['min_slope_deg']}°": static.Slope_deg < f["min_slope_deg"],
        "water (LULC)": static.LULC_Class.isin(f["water_lulc"]),
        "water/no-data lithology": static.Lithology_Class.isin(f["drop_lithology"]),
        "missing geotech params": static[["Ksat_ms", "Cohesion_kPa", "Porosity_n", "Friction_Angle_deg"]].isna().any(axis=1),
        "invalid soil texture": ~tex.between(90, 110),
        "no rainfall data": pd.Series(np.isnan(dyn["rain"]).all(axis=1)),
    }
    drop = np.zeros(N, bool)
    L += ["## Filters applied", "", "| filter | cells removed (alone) |", "|---|---|"]
    for k, m in masks.items():
        m = np.asarray(m)
        L.append(f"| {k} | {m.sum():,} |")
        drop |= m
    valid = ~drop
    np.save(p(cfg["paths"]["processed"]) / "valid_mask.npy", valid)
    urban = static.LULC_Class.isin(f["urban_lulc"]).to_numpy() & valid
    L += ["", f"**Valid cells kept: {valid.sum():,} / {N:,}** ({100*valid.mean():.1f} %). "
          f"Urban cells kept but flagged: {urban.sum():,}. "
          f"First {f['burn_in_days']} days are an antecedent-rain burn-in period.", ""]

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(grid_image(static, valid.astype(float)), extent=extent(static), cmap="Greys_r")
    ax.set_title("Valid cells (white) after filtering")
    fig.savefig(rep / "valid_mask.png", dpi=110)
    plt.close(fig)

    L += ["## Figures", "", "![](static_layers.png)", "![](rain_soilmoisture_timeseries.png)",
          "![](correlation.png)", "![](valid_mask.png)"]
    (rep / "eda_report.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:60]))
    print(f"\n... full report: {rep/'eda_report.md'}")


if __name__ == "__main__":
    main()
