# Physics-based rainfall-induced landslide susceptibility — Hulu Kelang / Selangor, Malaysia

This repo maps how susceptible each part of a 50 × 50 km box in Hulu Kelang / Selangor is to rainfall-induced landslides. The box covers 101.50–101.95 °E, 2.95–3.40 °N, on a 250 m grid with daily data for 2020.

The model combines two parts:
- **Green–Ampt infiltration** tracks how far rain soaks into the soil.
- **Infinite-slope Mohr–Coulomb stability** turns that into a factor of safety, FS(t).

FS is converted to a failure probability with `p = 1 / (1 + exp(k (FS − 1)))`.

## How the inputs are built

**Soil parameters.** Ksat, c′, n and φ′ come from pedotransfer and lithology relationships in `calculat.py`. They are not learned. The output file is `data/Malaysia_Unique_Grid_4_Geotech_Params.csv`.

**Predictors.** These come from the Google Earth Engine export, merged with GLiM lithology:
- terrain
- soil
- land use / land cover (LULC)
- daily rain
- antecedent rain
- soil moisture

**Labels.** Event locations come from the NASA Global Landslide Catalog. They cover 2007–2017, so they serve as a spatial inventory only.

## Pipeline

```bash
conda env create -f environment.yml   # or: pip install the packages listed there
bash run_all.sh                       # runs the tests, then steps 01 → 07
```

| step | script | output |
|---|---|---|
| 01 | `scripts/01_to_parquet.py` | joins the main CSV with the parameter CSV and writes `data/processed/` (static.parquet, dynamic.npz) |
| 02 | `scripts/02_eda_filter.py` | data analysis report and filtering mask in `reports/eda/` |
| 03 | `scripts/03_labels_glc.py` | NASA Global Landslide Catalog labels per cell (`reports/labels.md`) |
| 04 | `scripts/04_calibrate.py` | calibrates soil depth H, drainage time τ and ψ scale by spatial-block cross-validation (`configs/calibrated.yaml`) |
| 05 | `scripts/05_run_physics.py` | `outputs/susceptibility_points.csv`: lat/lon, parameters, minimum FS, p and class for each cell |
| 06 | `scripts/06_compare_baselines.py` | physics model compared with slope-only, random forest and XGBoost (`reports/metrics.md`) |
| 07 | `scripts/07_maps.py` | `outputs/maps/`: interactive HTML map, PNG, GeoTIFFs, GeoPackage |

The data files are not in the repo because they are too large. Put them in `data/`:
- `Malaysia_50x50km_250m_Daily_2020_WithLithology.csv`
- `Malaysia_Unique_Grid_4_Geotech_Params.csv`

## Current findings and limitations
- **No cell reaches FS < 1.** The minimum FS across all cells is 1.73. At 250 m the slopes are smoothed (max 34°), and c′ runs 8–20 kPa. So the map shows a **relative** susceptibility index: cells are ranked by their minimum FS, with class breaks at the 40/70/85/95th percentiles. The absolute p is kept in the CSV.
- **The catalog has no events in the 2020 forcing year.** Its events are biased toward urban cut slopes in KL. That makes the labels uninformative for calibration (physics AUC ≈ 0.39), and the literature defaults are kept.
- **Improving validation would need:**
  - a local inventory, for example from JMG or JKR (Malaysian geoscience and public works departments)
  - an Earth Engine re-export for the catalog years
