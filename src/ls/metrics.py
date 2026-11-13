"""Metrics and spatial-block cross-validation splits."""
import numpy as np
from sklearn.metrics import average_precision_score, f1_score, roc_auc_score


def spatial_blocks(lat, lon, n_blocks=4):
    """Assign each cell to one of n_blocks x n_blocks spatial blocks (block id 0..n^2-1)."""
    qa = np.clip(((lat - lat.min()) / (np.ptp(lat) + 1e-12) * n_blocks).astype(int), 0, n_blocks - 1)
    qo = np.clip(((lon - lon.min()) / (np.ptp(lon) + 1e-12) * n_blocks).astype(int), 0, n_blocks - 1)
    return qa * n_blocks + qo


def block_folds(blocks, y, n_folds=4, seed=0):
    """Group blocks into n_folds folds, trying to spread positives evenly."""
    rng = np.random.default_rng(seed)
    ub = np.unique(blocks)
    pos = np.array([(y[blocks == b] == 1).sum() for b in ub])
    order = ub[np.lexsort((rng.random(len(ub)), -pos))]
    fold_of_block, load = {}, np.zeros(n_folds)
    for b in order:
        f = int(np.argmin(load))
        fold_of_block[b] = f
        load[f] += (y[blocks == b] == 1).sum() + 1e-3
    return np.array([fold_of_block[b] for b in blocks])


def scores(y, s, thr=0.5):
    y = np.asarray(y)
    s = np.asarray(s)
    if len(np.unique(y)) < 2:
        return {"auc": float("nan"), "pr_auc": float("nan"), "f1@0.5": float("nan"), "n_pos": int(y.sum())}
    yhat = (s >= thr).astype(int)
    return {
        "auc": float(roc_auc_score(y, s)),
        "pr_auc": float(average_precision_score(y, s)),
        "f1@0.5": float(f1_score(y, yhat, zero_division=0)),
        "hit_rate@0.5": float(yhat[y == 1].mean()),
        "false_alarm@0.5": float(yhat[y == 0].mean()),
        "n_pos": int(y.sum()),
        "n": int(len(y)),
    }
