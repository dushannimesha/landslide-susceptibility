"""NASA Global Landslide Catalog -> per-cell labels."""
import numpy as np
import pandas as pd

ACC_KM = {"exact": 0.0, "1km": 1.0, "5km": 5.0, "10km": 10.0, "25km": 25.0, "50km": 50.0}


def load_glc(path, bbox, pad_deg=0.02):
    df = pd.read_csv(path, low_memory=False)
    df = df.dropna(subset=["latitude", "longitude"])
    x0, y0, x1, y1 = bbox
    m = df.longitude.between(x0 - pad_deg, x1 + pad_deg) & df.latitude.between(y0 - pad_deg, y1 + pad_deg)
    df = df[m].copy()
    df["date"] = pd.to_datetime(df["event_date"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce")
    df["acc_km"] = df["location_accuracy"].map(ACC_KM).fillna(np.inf)
    return df.sort_values("date").reset_index(drop=True)


def haversine_m(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * 6371000.0 * np.arcsin(np.sqrt(a))


def cell_labels(static, events, pos_radius_m=500, neg_buffer_m=1500):
    """Return label array per cell: 1 positive, 0 negative, -1 excluded (ambiguous buffer),
    and distance to nearest event (m)."""
    lat = static["Latitude"].to_numpy()
    lon = static["Longitude"].to_numpy()
    dmin = np.full(len(static), np.inf)
    for ev in events.itertuples():
        dmin = np.minimum(dmin, haversine_m(lat, lon, ev.latitude, ev.longitude))
    y = np.where(dmin <= pos_radius_m, 1, np.where(dmin <= neg_buffer_m, -1, 0))
    return y.astype(np.int8), dmin
