"""Shared helpers for the experiment scripts."""
from __future__ import annotations

import json
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

try:
    import matplotlib  # noqa: E402
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402,F401
    HAVE_MPL = True
except Exception:  # pragma: no cover
    HAVE_MPL = False


def outdir(path: str) -> str:
    p = path if os.path.isabs(path) else os.path.join(ROOT, path)
    os.makedirs(p, exist_ok=True)
    return p


def dump_json(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=conv)


class Timer:
    def __init__(self):
        self.t0 = time.perf_counter()

    def __call__(self):
        return time.perf_counter() - self.t0


def vg_grid(vmin=-1.5, vmax=2.5, step=0.05):
    return np.round(np.arange(vmin, vmax + 1e-9, step), 6)
