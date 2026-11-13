"""Loading the processed arrays produced by scripts/01_to_parquet.py."""
import json

import numpy as np
import pandas as pd

from . import p


def load_processed(cfg, filtered=True):
    """Return (static DataFrame, dynamic dict of [N, T] arrays, dates).

    static rows and dynamic first axis share the same order (cell_id).
    """
    d = p(cfg["paths"]["processed"])
    static = pd.read_parquet(d / "static.parquet")
    dyn = dict(np.load(d / "dynamic.npz"))
    dates = pd.to_datetime(dyn.pop("dates").astype(str), format="%Y%m%d")
    if filtered:
        mask_path = d / "valid_mask.npy"
        if mask_path.exists():
            m = np.load(mask_path)
            static = static[m].reset_index(drop=True)
            dyn = {k: v[m] for k, v in dyn.items()}
    return static, dyn, dates


def save_json(obj, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=float)
