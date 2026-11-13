"""Step 3: NASA Global Landslide Catalog (GLC) -> landslide inventory + per-cell labels.

The GLC export ends in 2017, so no event overlaps the 2020 forcing. Event *locations* are used as a
spatial inventory (standard for susceptibility mapping); the 2020 daily series is the physical
forcing applied to every cell.
Outputs: data/processed/glc_events.csv, data/processed/labels.parquet, reports/labels.md
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ls import load_config, p  # noqa: E402
from ls.io import load_processed  # noqa: E402
from ls.labels import cell_labels, load_glc  # noqa: E402


def main():
    cfg = load_config()
    lc = cfg["labels"]
    path = p(cfg["paths"]["glc_csv"])
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["curl", "-sL", "-o", str(path), cfg["paths"]["glc_url"]], check=True)

    static, _, dates = load_processed(cfg, filtered=True)
    ev = load_glc(path, cfg["bbox"])
    # keep only events that fall inside the actual grid footprint
    inside = (ev.latitude.between(static.Latitude.min(), static.Latitude.max())
              & ev.longitude.between(static.Longitude.min(), static.Longitude.max()))
    ev["inside_grid"] = inside
    ev["rain_triggered"] = ev.landslide_trigger.isin(lc["rain_triggers"])
    ev["used_for_labels"] = inside & ev.rain_triggered & (ev.acc_km <= lc["max_location_accuracy_km"])
    keep = ["event_id", "date", "latitude", "longitude", "location_accuracy", "acc_km", "landslide_trigger",
            "landslide_size", "fatality_count", "event_title", "inside_grid", "rain_triggered", "used_for_labels"]
    ev[keep].to_csv(p(cfg["paths"]["processed"]) / "glc_events.csv", index=False)

    use = ev[ev.used_for_labels]
    y, dmin = cell_labels(static, use, lc["positive_radius_m"], lc["negative_buffer_m"])
    lab = pd.DataFrame({"cell_id": static.cell_id, "label": y, "dist_to_event_m": dmin})
    lab.to_parquet(p(cfg["paths"]["processed"]) / "labels.parquet", index=False)

    yr = ev.date.dt.year
    L = ["# Landslide labels (NASA GLC)", "",
         f"- Events in bbox: {len(ev)}; inside grid: {inside.sum()}; rain-triggered: {ev.rain_triggered.sum()}",
         f"- Used for labels (inside grid, rain, accuracy <= {lc['max_location_accuracy_km']} km): **{len(use)}**",
         f"- Event years: {yr.min()}–{yr.max()}  |  events in forcing period "
         f"({dates.min().year}): {(yr == dates.min().year).sum()}",
         f"- Positive cells (<= {lc['positive_radius_m']} m): **{(y == 1).sum()}**, excluded buffer: {(y == -1).sum()}, "
         f"negatives: {(y == 0).sum()}", "",
         "Events per year:", "", ev.groupby(yr).size().to_frame("events").to_markdown(), ""]
    if (yr == dates.min().year).sum() == 0:
        L += ["> **Note:** no catalogued events in the forcing year. Labels are spatial (historical inventory "
              "2007–2017). For event-date (temporal) validation, re-export the GEE dataset for 2007–2017 "
              "with the same columns; the pipeline accepts any year range unchanged."]
    (p(cfg["paths"]["reports"]) / "labels.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
