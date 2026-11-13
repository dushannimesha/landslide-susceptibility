import time
import numpy as np
import pandas as pd

t0 = time.time()

# -------------------------------------------------------------------------
# 1. INPUT & OUTPUT FILE PATHS
# -------------------------------------------------------------------------
csv_in = "Malaysia_50x50km_250m_Daily_2020_WithLithology.csv"

# Output 1: Full daily rows with ONLY Coordinates, Date, and the 4 Parameters
csv_out_daily = "Malaysia_2020_Only_4_Geotech_Params_Daily.csv"

# Output 2: Unique spatial grid points with the 4 Parameters (tiny, fast lookup file)
csv_out_unique = "Malaysia_Unique_Grid_4_Geotech_Params.csv"

print(f"[1/3] Reading required columns from {csv_in} ...")
# Load only the necessary columns to make it run super fast and save RAM
needed_cols = [
    "Year", "Date", "Latitude", "Longitude",
    "Soil_Sand_pct", "Soil_Clay_pct", "Bulk_Density_gcm3",
    "Lithology_Class", "LULC_Class"
]
df = pd.read_csv(csv_in, usecols=needed_cols)
print(f"      -> Loaded {len(df):,} rows in {time.time() - t0:.1f}s")

# -------------------------------------------------------------------------
# 2. CALCULATE THE 4 GEOTECHNICAL PARAMETERS
# -------------------------------------------------------------------------
print("[2/3] Calculating Porosity (n), Ksat, Cohesion (c'), and Friction Angle (phi')...")

sand = df["Soil_Sand_pct"].clip(1.0, 98.0)
clay = df["Soil_Clay_pct"].clip(1.0, 98.0)
bdod = df["Bulk_Density_gcm3"].clip(0.8, 2.1)

# 1. Porosity (n) [dimensionless fraction, 0 to 1]
df["Porosity_n"] = (1.0 - (bdod / 2.65)).clip(0.25, 0.65).round(4)

# 2. Saturated Hydraulic Conductivity (Ksat) [m/s]
log10_ksat_in_hr = -0.60 + (0.0126 * sand) - (0.0064 * clay)
ksat_mm_hr = 25.4 * (10.0 ** log10_ksat_in_hr)
df["Ksat_ms"] = (ksat_mm_hr / 3.6e6).clip(1e-8, 1e-3)

# 3. Effective Cohesion (c') [kPa]
base_cohesion_map = {
    7: 8.5,   # pa: Granite residual soil (Hulu Kelang / Titiwangsa)
    6: 10.5,  # mt: Metamorphic Schist / Phyllite
    3: 12.0,  # sm: Mixed Sedimentary (Kenny Hill Formation)
    1: 5.0,   # su: Unconsolidated Alluvium
    0: 8.5    # Fallback
}
c_base = df["Lithology_Class"].map(base_cohesion_map).fillna(8.5)
root_cohesion = np.where(df["LULC_Class"] == 10, 2.5, 0.5)  # +2.5 kPa for Forest/Trees
df["Cohesion_kPa"] = (c_base + 0.12 * (clay - 25.0) + root_cohesion).clip(2.0, 28.0).round(2)

# 4. Effective Internal Friction Angle (phi') [degrees]
base_friction_map = {
    7: 31.5,  # pa: Granitic sandy residual soil
    6: 26.5,  # mt: Metamorphic Schist
    3: 28.5,  # sm: Mixed Sedimentary
    1: 24.0,  # su: Unconsolidated Alluvium
    0: 28.0    # Fallback
}
phi_base = df["Lithology_Class"].map(base_friction_map).fillna(28.0)
df["Friction_Angle_deg"] = (phi_base + 0.10 * (sand - 45.0) - 0.08 * (clay - 25.0)).clip(18.0, 42.0).round(2)

# -------------------------------------------------------------------------
# 3. EXPORT SEPARATE CSV FILES WITH ONLY THE 4 OUTPUTS (+ ID KEYS)
# -------------------------------------------------------------------------
output_columns = [
    "Year", "Date", "Latitude", "Longitude",
    "Ksat_ms", "Cohesion_kPa", "Porosity_n", "Friction_Angle_deg"
]

print(f"[3/3] Saving separate daily 4-parameter CSV to: {csv_out_daily} ...")
df[output_columns].to_csv(csv_out_daily, index=False)

# Also save a small lookup CSV with just the unique (Latitude, Longitude) grid points
unique_cols = ["Latitude", "Longitude", "Ksat_ms", "Cohesion_kPa", "Porosity_n", "Friction_Angle_deg"]
df_unique = df[unique_cols].drop_duplicates(subset=["Latitude", "Longitude"])
df_unique.to_csv(csv_out_unique, index=False)

print(f"\nDone in {time.time() - t0:.1f}s!")
print(f" -> Saved full daily file: {csv_out_daily} ({len(df):,} rows)")
print(f" -> Saved unique grid lookup file: {csv_out_unique} ({len(df_unique):,} unique locations)\n")
print("Preview of the 4 calculated outputs:")
print(df_unique.head(10).to_string(index=False))