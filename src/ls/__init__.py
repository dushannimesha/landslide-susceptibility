"""Physics-based rainfall-induced landslide susceptibility (Green-Ampt + infinite slope)."""
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


def load_config(path=None, use_calibrated=True):
    path = Path(path) if path else ROOT / "configs" / "default.yaml"
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cal = ROOT / "configs" / "calibrated.yaml"
    if use_calibrated and cal.exists():
        with open(cal) as f:
            c = yaml.safe_load(f)
        cfg["physics"].update(c.get("physics", {}))
        cfg["calibration_result"] = c.get("calibration_result", {})
    return cfg


def p(rel):
    """Resolve a project-relative path."""
    return ROOT / rel
