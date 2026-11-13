"""Step 7: susceptibility maps with coordinates.

outputs/maps/
  susceptibility_map.html     interactive Leaflet map: coloured points with lat/lon popups, layer
                              switcher (classes, min FS, parameters, wettest-day p), GLC events,
                              mouse lat/lon readout, 0.05 deg graticule
  susceptibility_map.png      static figure with lat/lon axes, grid, legend, north arrow, scale bar
  susceptibility_maxp.tif     GeoTIFF (EPSG:4326) of max p ; susceptibility_class.tif
  susceptibility_points.gpkg  points layer for QGIS / ArcGIS
"""
import base64
import io
import sys
from pathlib import Path

import branca
import folium
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
from folium.plugins import Fullscreen, MeasureControl, MousePosition  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from rasterio.transform import from_origin  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ls import load_config, p  # noqa: E402

CLASS_COLORS = ["#1a9850", "#a6d96a", "#fee08b", "#f46d43", "#a50026"]
LITHO = {1: "su Unconsolidated sediments", 2: "ss Siliciclastic sed.", 3: "sm Mixed sedimentary",
         4: "sc Carbonate", 6: "mt Metamorphic", 7: "pa Acid plutonic (granite)", 0: "nd / water"}


def grid_geom(pts):
    lats = np.sort(pts.Latitude.unique())
    lons = np.sort(pts.Longitude.unique())
    dlat, dlon = np.median(np.diff(lats)), np.median(np.diff(lons))
    nr, nc = int(pts.row.max()) + 1, int(pts.col.max()) + 1
    lat_top = pts.Latitude.max() + dlat / 2
    lon_left = pts.Longitude.min() - dlon / 2
    bounds = [[lat_top - nr * dlat, lon_left], [lat_top, lon_left + nc * dlon]]
    return nr, nc, dlat, dlon, lat_top, lon_left, bounds


def to_img(pts, values, nr, nc):
    img = np.full((nr, nc), np.nan, dtype=np.float32)
    img[pts.row.to_numpy(), pts.col.to_numpy()] = values
    return img


def png_data_url(img, cmap, vmin, vmax, norm=None):
    rgba = (plt.get_cmap(cmap) if isinstance(cmap, str) else cmap)(
        norm(img) if norm is not None else (img - vmin) / max(vmax - vmin, 1e-12))
    rgba[np.isnan(img)] = 0
    buf = io.BytesIO()
    plt.imsave(buf, rgba, format="png")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def add_overlay(m, name, img, cmap, vmin, vmax, bounds, caption, show=False, norm=None):
    fg = folium.FeatureGroup(name=name, show=show)
    folium.raster_layers.ImageOverlay(png_data_url(img, cmap, vmin, vmax, norm), bounds=bounds, opacity=0.75,
                                      interactive=False, zindex=2).add_to(fg)
    fg.add_to(m)
    if norm is None:
        colors = [matplotlib.colors.to_hex(plt.get_cmap(cmap)(x)) for x in np.linspace(0, 1, 9)]
        cm = branca.colormap.LinearColormap(colors, vmin=vmin, vmax=vmax, caption=caption)
        return cm
    return None


def interactive_map(cfg, pts, events, nr, nc, bounds, wet_days, daily_p, dates, out):
    centre = [pts.Latitude.mean(), pts.Longitude.mean()]
    m = folium.Map(location=centre, zoom_start=11, tiles=None, prefer_canvas=True, control_scale=True)
    folium.TileLayer("OpenStreetMap", name="OpenStreetMap").add_to(m)
    folium.TileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
                     attr="Esri World Imagery", name="Satellite (Esri)").add_to(m)
    folium.TileLayer("https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png", attr="OpenTopoMap", name="Topographic").add_to(m)

    # ---- point layer: every cell, coloured by susceptibility class, popup with coordinates
    names = cfg["classes"]["names"]
    col = dict(zip(names, CLASS_COLORS))
    fields = ["Latitude", "Longitude", "susceptibility_class", "susceptibility_index", "FS_min_2020", "FS_dry",
              "date_p_max", "p_max_2020", "Slope_deg", "Lithology", "Ksat_ms", "Cohesion_kPa", "Friction_Angle_deg",
              "Porosity_n"]
    g = pts.assign(Lithology=pts.Lithology_Class.map(LITHO)).copy()
    for c, d in [("Latitude", 5), ("Longitude", 5), ("susceptibility_index", 3), ("FS_min_2020", 2),
                 ("FS_dry", 2), ("Slope_deg", 1), ("Cohesion_kPa", 1), ("Friction_Angle_deg", 1), ("Porosity_n", 3)]:
        g[c] = g[c].round(d)
    g["Ksat_ms"] = g["Ksat_ms"].map(lambda v: f"{v:.2e}")
    g["p_max_2020"] = g["p_max_2020"].map(lambda v: f"{v:.2e}")
    g["susceptibility_class"] = g["susceptibility_class"].astype(str)
    gdf = gpd.GeoDataFrame(g[fields], geometry=gpd.points_from_xy(pts.Longitude, pts.Latitude), crs="EPSG:4326")
    aliases = ["Lat (°N)", "Lon (°E)", "Susceptibility class", "Susceptibility index (0–1)", "Min FS 2020",
               "FS at start (dry)", "Date of min FS", "Abs. failure prob. p", "Slope (°)", "Lithology",
               "Ksat (m/s)", "c' (kPa)", "φ' (°)", "n"]

    for cls in names[::-1]:
        sub = gdf[gdf.susceptibility_class == cls]
        if sub.empty:
            continue
        fg = folium.FeatureGroup(name=f"Points: {cls} ({len(sub):,})", show=True)
        folium.GeoJson(
            sub.to_json(),
            marker=folium.CircleMarker(radius=3, weight=0, fill=True, fill_opacity=0.85),
            style_function=lambda f, c=col[cls]: {"fillColor": c, "color": c},
            popup=folium.GeoJsonPopup(fields=fields, aliases=aliases, localize=True),
            tooltip=folium.GeoJsonTooltip(fields=["Latitude", "Longitude", "susceptibility_class", "FS_min_2020"],
                                          aliases=["Lat", "Lon", "Class", "Min FS"]),
        ).add_to(fg)
        fg.add_to(m)

    # ---- raster overlays
    legends = []
    legends.append(add_overlay(m, "Raster: min FS 2020", to_img(pts, np.clip(pts.FS_min_2020, 0.5, 3), nr, nc),
                               "RdYlGn", 0.5, 3.0, bounds, "Min factor of safety (2020)"))
    legends.append(add_overlay(m, "Raster: log10 Ksat", to_img(pts, np.log10(pts.Ksat_ms), nr, nc), "Blues",
                               np.log10(pts.Ksat_ms).min(), np.log10(pts.Ksat_ms).max(), bounds, "log10 Ksat (m/s)"))
    for c, cm, cap in [("Cohesion_kPa", "YlOrBr", "Effective cohesion c' (kPa)"),
                       ("Friction_Angle_deg", "PuRd", "Friction angle φ' (°)"),
                       ("Porosity_n", "GnBu", "Porosity n"), ("Slope_deg", "magma_r", "Slope (°)")]:
        legends.append(add_overlay(m, f"Raster: {c}", to_img(pts, pts[c], nr, nc), cm, pts[c].min(), pts[c].max(),
                                   bounds, cap))
    for t in wet_days:
        d = pd.Timestamp(str(dates[t])).strftime("%Y-%m-%d")
        add_overlay(m, f"Raster: p on {d} (wet day)", to_img(pts, daily_p[:, t].astype(np.float32), nr, nc),
                    ListedColormap(CLASS_COLORS), 0, 1, bounds, "", norm=BoundaryNorm([0, .2, .4, .6, .8, 1.0], 5))

    # ---- GLC events
    fg = folium.FeatureGroup(name="NASA GLC landslides 2007–2017", show=True)
    for e in events.itertuples():
        folium.Marker([e.latitude, e.longitude],
                      icon=folium.Icon(color="black" if e.used_for_labels else "gray", icon="warning-sign"),
                      popup=folium.Popup(f"<b>{e.event_title}</b><br>{str(e.date)[:10]}<br>"
                                         f"{e.latitude:.5f} °N, {e.longitude:.5f} °E<br>accuracy: {e.location_accuracy}"
                                         f"<br>trigger: {e.landslide_trigger}", max_width=280)).add_to(fg)
    fg.add_to(m)

    # ---- graticule (0.05 deg) with labels
    fg = folium.FeatureGroup(name="Lat/Lon grid 0.05°", show=True)
    (s, w), (n, e) = bounds
    for lat in np.arange(np.ceil(s / 0.05) * 0.05, n, 0.05):
        folium.PolyLine([[lat, w], [lat, e]], color="#555", weight=0.6, opacity=0.6, dash_array="4").add_to(fg)
        folium.Marker([lat, w], icon=folium.DivIcon(html=f"<div style='font-size:10px;color:#222;background:#fff8;"
                                                         f"width:48px'>{lat:.2f}°N</div>")).add_to(fg)
    for lon in np.arange(np.ceil(w / 0.05) * 0.05, e, 0.05):
        folium.PolyLine([[s, lon], [n, lon]], color="#555", weight=0.6, opacity=0.6, dash_array="4").add_to(fg)
        folium.Marker([s, lon], icon=folium.DivIcon(html=f"<div style='font-size:10px;color:#222;background:#fff8;"
                                                         f"width:56px'>{lon:.2f}°E</div>")).add_to(fg)
    fg.add_to(m)
    folium.Rectangle(bounds, color="#000", weight=1.5, fill=False, tooltip="Study area").add_to(m)

    MousePosition(position="bottomleft", separator=" | ", prefix="Lat | Lon:", num_digits=5).add_to(m)
    Fullscreen().add_to(m)
    MeasureControl(primary_length_unit="kilometers").add_to(m)
    for cm in legends:
        if cm is not None:
            cm.add_to(m)
    counts = pts.susceptibility_class.value_counts()
    legend = ("<div style='position:fixed;top:12px;left:56px;z-index:9999;background:white;padding:8px 10px;"
              "border:1px solid #888;border-radius:6px;font:12px sans-serif'><b>Landslide susceptibility "
              "(2020)</b><br>" +
              "".join(f"<span style='display:inline-block;width:12px;height:12px;background:{c};margin-right:6px'>"
                      f"</span>{nm} ({counts.get(nm, 0):,})<br>" for nm, c in zip(names, CLASS_COLORS)) +
              "<span style='font-size:10px'>Classes = percentile of min FS(t) from Green–Ampt +<br>"
              "infinite slope (breaks 40/70/85/95 %). Click a point for<br>coordinates, FS, p and parameters."
              "</span></div>")
    m.get_root().html.add_child(folium.Element(legend))
    folium.LayerControl(collapsed=False).add_to(m)
    m.fit_bounds(bounds)
    m.save(out / "susceptibility_map.html")


def static_map(cfg, pts, events, nr, nc, bounds, out):
    names = cfg["classes"]["names"]
    cls = to_img(pts, pts.susceptibility_class.cat.codes.astype(float), nr, nc)
    (s, w), (n, e) = bounds
    fig, axs = plt.subplots(1, 2, figsize=(17, 8.2))
    ax = axs[0]
    ax.imshow(cls, extent=[w, e, s, n], cmap=ListedColormap(CLASS_COLORS), vmin=-0.5, vmax=4.5,
              interpolation="nearest")
    ax.set_title("Relative susceptibility class (percentile of min FS, 2020 forcing)", fontsize=12)
    ev_used = events[events.used_for_labels]
    ax.scatter(events.longitude, events.latitude, marker="^", s=40, facecolor="none", edgecolor="#444", lw=1,
               label="GLC event (other)")
    ax.scatter(ev_used.longitude, ev_used.latitude, marker="^", s=55, c="k", label="GLC event (≤1 km, used)")
    handles = [Patch(color=c, label=f"{nm} ({(pts.susceptibility_class == nm).sum():,})")
               for nm, c in zip(names, CLASS_COLORS)]
    ax.legend(handles=handles + ax.get_legend_handles_labels()[0], loc="lower right", fontsize=8, framealpha=0.9)

    ax = axs[1]
    im = ax.imshow(to_img(pts, np.clip(pts.FS_min_2020, 0.5, 3), nr, nc), extent=[w, e, s, n], cmap="RdYlGn",
                   vmin=0.5, vmax=3, interpolation="nearest")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Minimum factor of safety (2020)")
    ax.scatter(ev_used.longitude, ev_used.latitude, marker="^", s=55, c="k")
    ax.set_title("Minimum dynamic FS(t), 2020", fontsize=12)

    for ax in axs:
        ax.set_xlim(w, e)
        ax.set_ylim(s, n)
        ax.set_xlabel("Longitude (°E)")
        ax.set_ylabel("Latitude (°N)")
        ax.grid(True, ls="--", lw=0.5, color="k", alpha=0.4)
        ax.set_xticks(np.arange(np.ceil(w * 20) / 20, e, 0.05))
        ax.set_yticks(np.arange(np.ceil(s * 20) / 20, n, 0.05))
        ax.ticklabel_format(useOffset=False)
        # north arrow
        ax.annotate("N", xy=(0.06, 0.93), xytext=(0.06, 0.82), xycoords="axes fraction", ha="center",
                    fontsize=13, fontweight="bold", arrowprops=dict(facecolor="k", width=4, headwidth=12))
        # 10 km scale bar
        km_deg = 10 / (111.32 * np.cos(np.radians((s + n) / 2)))
        x0, y0 = w + 0.03 * (e - w), s + 0.04 * (n - s)
        ax.plot([x0, x0 + km_deg], [y0, y0], color="k", lw=4)
        ax.text(x0 + km_deg / 2, y0 + 0.01 * (n - s), "10 km", ha="center", fontsize=9, fontweight="bold")
    fig.suptitle("Hulu Kelang / Selangor, Malaysia — physics-based (Green–Ampt + infinite slope) landslide "
                 "susceptibility", fontsize=13)
    fig.tight_layout()
    fig.savefig(out / "susceptibility_map.png", dpi=220)
    plt.close(fig)


def write_rasters(pts, nr, nc, dlat, dlon, lat_top, lon_left, out):
    tr = from_origin(lon_left, lat_top, dlon, dlat)
    for name, vals, dt in [("susceptibility_maxp.tif", pts.p_max_2020, "float32"),
                           ("susceptibility_index.tif", pts.susceptibility_index, "float32"),
                           ("fs_min.tif", pts.FS_min_2020, "float32"),
                           ("susceptibility_class.tif", pts.susceptibility_class.cat.codes + 1, "float32")]:
        img = to_img(pts, vals.astype(np.float32), nr, nc)
        with rasterio.open(out / name, "w", driver="GTiff", height=nr, width=nc, count=1, dtype=dt,
                           crs="EPSG:4326", transform=tr, nodata=np.nan, compress="deflate") as dst:
            dst.write(img, 1)


def main():
    cfg = load_config()
    o = p(cfg["paths"]["outputs"])
    out = o / "maps"
    out.mkdir(parents=True, exist_ok=True)
    pts = pd.read_csv(o / "susceptibility_points.csv")
    pts["susceptibility_class"] = pd.Categorical(pts.susceptibility_class, categories=cfg["classes"]["names"],
                                                 ordered=True)
    events = pd.read_csv(p(cfg["paths"]["processed"]) / "glc_events.csv")
    events = events[events.inside_grid]
    daily_p = np.load(o / "daily_p.npy")
    dates = np.load(o / "dates.npy")
    # 3 days with the highest mean p
    wet_days = list(np.argsort(-daily_p.astype(np.float32).mean(axis=0))[:3])

    nr, nc, dlat, dlon, lat_top, lon_left, bounds = grid_geom(pts)
    print("rasters ...", flush=True)
    write_rasters(pts, nr, nc, dlat, dlon, lat_top, lon_left, out)
    print("geopackage ...", flush=True)
    gpd.GeoDataFrame(pts.assign(susceptibility_class=pts.susceptibility_class.astype(str)),
                     geometry=gpd.points_from_xy(pts.Longitude, pts.Latitude), crs="EPSG:4326"
                     ).to_file(out / "susceptibility_points.gpkg", layer="susceptibility", driver="GPKG")
    print("static map ...", flush=True)
    static_map(cfg, pts, events, nr, nc, bounds, out)
    print("interactive map ...", flush=True)
    interactive_map(cfg, pts, events, nr, nc, bounds, wet_days, daily_p, dates, out)
    for f in sorted(out.iterdir()):
        print(f"  {f.name:32s} {f.stat().st_size/1e6:7.2f} MB")


if __name__ == "__main__":
    main()
