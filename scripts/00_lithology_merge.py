import pandas as pd
import geopandas as gpd
from shapely.geometry import box

# -------------------------------------------------------------------------
# 1. FILE PATHS
# -------------------------------------------------------------------------
gdb_path = "LiMW_GIS 2015.gdb"  # Or "LiMW_GIS 2015.gdb.zip"
csv_input_path = "Malaysia_50x50km_250m_Daily_2020.csv"
csv_output_path = "Malaysia_50x50km_250m_Daily_2020_WithLithology.csv"

# -------------------------------------------------------------------------
# 2. CREATE A PROPER EPSG:4326 MASK SO GEOPANDAS REPROJECTS AUTOMATICALLY
# -------------------------------------------------------------------------
print("Loading GLiM polygons for the 50x50 km Malaysia box...")
roi_poly = gpd.GeoDataFrame(
    geometry=[box(101.45, 2.90, 102.00, 3.45)], 
    crs="EPSG:4326"
)

# Passing mask=roi_poly automatically translates degrees to the GDB's native CRS!
gdf_glim = gpd.read_file(gdb_path, mask=roi_poly)
print(f"Loaded {len(gdf_glim)} geology polygons inside your Malaysia region!")

# Now convert those loaded polygons to EPSG:4326 (Lat/Lon degrees) to match your CSV
gdf_glim = gdf_glim.to_crs("EPSG:4326")

# -------------------------------------------------------------------------
# 3. EXTRACT THE FIRST 2 LETTERS OF THE 'xx' LITHOLOGY CODE
# -------------------------------------------------------------------------
glim_code_to_id = {
    "su": 1,   # Unconsolidated sediments (Alluvium)
    "ss": 2,   # Siliciclastic sedimentary rocks (Sandstone, Mudstone)
    "sm": 3,   # Mixed sedimentary rocks
    "sc": 4,   # Carbonate sedimentary rocks (Limestone)
    "ev": 5,   # Evaporites
    "mt": 6,   # Metamorphic rocks (Schist, Phyllite, Quartzite)
    "pa": 7,   # Acid plutonic rocks (Granite - Hulu Kelang / Titiwangsa hills)
    "pi": 8,   # Intermediate plutonic rocks
    "pb": 9,   # Basic plutonic rocks (Gabbro)
    "va": 10,  # Acid volcanic rocks
    "vi": 11,  # Intermediate volcanic rocks
    "vb": 12,  # Basic volcanic rocks (Basalt)
    "py": 13,  # Pyroclastics
    "wb": 0,   # Water bodies
    "nd": 0    # No data
}

# Find the 'xx' column (or 'Identity') and take the first 2 characters (.str[:2])
rock_col = [c for c in gdf_glim.columns if c.lower() in ["xx", "identity"]][0]
gdf_glim["Lithology_Code"] = gdf_glim[rock_col].astype(str).str.strip().str.lower().str[:2]
gdf_glim["Lithology_Class"] = gdf_glim["Lithology_Code"].map(glim_code_to_id).fillna(0).astype(int)

print("\nRock Classes Found in Your 50x50 km Selangor Study Area:")
print(gdf_glim[["Lithology_Code", "Lithology_Class"]].drop_duplicates().to_string(index=False))

# -------------------------------------------------------------------------
# 4. SPATIAL JOIN WITH YOUR GEE CSV POINTS (WITH NEAREST FALLBACK)
# -------------------------------------------------------------------------
print(f"\nLoading GEE CSV: {csv_input_path} ...")
df_csv = pd.read_csv(csv_input_path)

unique_coords = df_csv[["Latitude", "Longitude"]].drop_duplicates().copy()
gdf_points = gpd.GeoDataFrame(
    unique_coords,
    geometry=gpd.points_from_xy(unique_coords["Longitude"], unique_coords["Latitude"]),
    crs="EPSG:4326"
)

# Use sjoin_nearest so even points sitting right on polygon border lines get matched!
joined_points = gpd.sjoin_nearest(
    gdf_points,
    gdf_glim[["Lithology_Class", "Lithology_Code", "geometry"]],
    how="left"
)

# Drop any duplicate matches on borders
joined_points = joined_points[~joined_points.index.duplicated(keep="first")]

# -------------------------------------------------------------------------
# 5. MERGE BACK INTO FULL DAILY TIME-SERIES CSV & SAVE
# -------------------------------------------------------------------------
df_final = df_csv.merge(
    joined_points[["Latitude", "Longitude", "Lithology_Class", "Lithology_Code"]],
    on=["Latitude", "Longitude"],
    how="left"
)

df_final.to_csv(csv_output_path, index=False)
print(f"\nSuccess! Saved updated dataset to: {csv_output_path}")
print("\nValue Counts of Lithology_Code across your dataset:")
print(df_final["Lithology_Code"].value_counts())
print("\nFirst 10 rows preview:")
print(df_final[["Latitude", "Longitude", "Slope_deg", "Lithology_Class", "Lithology_Code"]].head(10))