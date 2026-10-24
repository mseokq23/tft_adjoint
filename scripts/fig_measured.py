"""Figure: open-source model calibrated to the digitised measurements of arXiv:2609.12184.

(a) 2 %-W calibration curves: measured vs calibrated DEVSIM model (and the paper's Step-5
    parameters in DEVSIM, to show that descriptor values do not transfer between simulators)
(b) hold-out conditions predicted by the calibrated model vs measured
(c) measured-data profile of ND,eff: Δχ² (σ_eff) of every row of results/e12_nd
(d) hold-out Vth predicted by every row of the profile (full simulated curves), relative to the median
    over the rows; the digitised hold-out curves have gaps around 100 nA, so a ΔVth against them
    would be an interpolation across the gap

Reads results/e10, results/e11_2W_v2 and results/e12_nd.
Usage:  python scripts/fig_measured.py [--out results/fig_measured.png]
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from _common import ROOT  # noqa: F401  (puts the repository root on sys.path)
from tftsim.metrics import extract

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"


def style(ax):
    ax.grid(True, color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)


def load_meas(path="data/meas.csv"):
    d = collections.defaultdict(list)
    for r in csv.DictReader(open(path)):
        d[(r["set"], float(r["L_um"]), float(r["vd"]))].append((float(r["vg"]), float(r["id"])))
    return {k: np.array(sorted(v)) for k, v in d.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit2w", default="results/e11_2W_v2/summary.json")
    ap.add_argument("--profile", default="results/e12_nd/summary.json")
    ap.add_argument("--out", default="results/fig_measured.png")
    args = ap.parse_args()
    meas = load_meas()
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "legend.fontsize": 7})
    fig, ax = plt.subplots(2, 2, figsize=(7.6, 6.0))
    a, b, c, d = ax.ravel()

    fit = json.load(open(args.fit2w))
    vg = np.array(fit["vg"])
    e10 = json.load(open("results/e10/summary.json")) if os.path.exists("results/e10/summary.json") else None
    for vd, col in (("1.0", BLUE), ("0.05", ORANGE)):
        m = meas[("2W", 1.0, float(vd))]
        a.semilogy(m[:, 0], m[:, 1], "s", mfc="none", color=col, ms=3.5, label=f"measured VD = {vd} V")
        a.semilogy(vg, fit["curves"][vd], "-", color=col, lw=1.6, label=f"calibrated (DEVSIM) VD = {vd} V")
        if e10 and "step5" in e10["curves"]:
            cv = e10["curves"]["step5"][vd]
            a.semilogy(cv[0], cv[1], ":", color=col, lw=1.2, label=f"paper Step-5 values in DEVSIM" if vd == "1.0" else None)
    a.set_ylim(1e-12, 2e-3)
    a.set_xlabel("VG (V)")
    a.set_ylabel("ID (A)")
    a.set_title("(a) 2%-W calibration curves (L = 1 µm)", loc="left")
    a.legend(loc="lower right", fontsize=6)
    style(a)

    cols = {"L1_VD0.1": AQUA, "L1_VD0.2": "#eda100", "L1_VD0.5": "#6250d6", "L2_VD0.05": ORANGE, "L2_VD1": BLUE}
    for k, h in fit.get("holdout", {}).items():
        if "curve" not in h:
            continue
        tag, vds = k.split("_VD")
        L = 2.0 if tag == "L2" else 1.0
        key = ("2W", L, float(vds))
        if key not in meas:
            continue
        m = meas[key]
        b.semilogy(m[:, 0], m[:, 1], "o", mfc="none", color=cols.get(k, INK2), ms=3)
        b.semilogy(vg, h["curve"], "-", color=cols.get(k, INK2), lw=1.4, label=k.replace("_", " "))
    b.set_ylim(1e-12, 2e-3)
    b.set_xlabel("VG (V)")
    b.set_title("(b) hold-out conditions: prediction vs measured", loc="left")
    b.legend(loc="lower right", fontsize=6)
    style(b)

    if os.path.exists(args.profile):
        pr = json.load(open(args.profile))
        rows = sorted(pr["rows"], key=lambda r: r["value"])
        x = [r["value"] for r in rows]
        c.plot(x, [r["dchi2_eff"] for r in rows], "o-", color=BLUE, lw=1.8, ms=5, label="measured calibration curves")
        c.axhline(3.84, color=INK2, ls=":", lw=1)
        c.text(x[-1], 3.95, "95 % threshold", fontsize=6.5, color=INK2, ha="right", va="bottom")
        w4 = np.log10(6.39e16)
        c.axvline(w4, color=INK2, ls="--", lw=0.8)
        c.text(w4 + 0.04, 2.3, "paper's 4%-W\nvalue", fontsize=6.5, color=INK2)
        c.set_ylim(-0.2, 5)
        c.set_xlim(min(x[0], w4) - 0.2, x[-1] + 0.15)
        c.set_xlabel("log10 ND,eff (cm⁻³) fixed; other 8 re-fitted")
        c.set_ylabel("Δχ² (model error included)")
        c.set_title("(c) measured-data profile of ND,eff", loc="left")
        c.legend(loc="upper left", fontsize=6.5)
        style(c)
        keys = ["L1_VD0.1", "L1_VD0.2", "L1_VD0.5", "L2_VD0.05", "L2_VD1"]
        sets = [json.load(open(pr["best"]))] + [json.load(open(f)) for f in sorted(
            glob.glob(os.path.join(os.path.dirname(args.profile), "log10_nd_eff_*", "summary.json")))]
        sets.sort(key=lambda s_: s_["theta_hat"]["log10_nd_eff"])
        xs = [s_["theta_hat"]["log10_nd_eff"] for s_ in sets]
        for k in keys:
            v = np.array([1e3 * extract(np.array(s_["vg"]), np.array(s_["holdout"][k]["curve"]))["vth"] for s_ in sets])
            d.plot(xs, v - np.median(v), "o-", color=cols.get(k, INK2), lw=1.4, ms=4, label=k.replace("_", " "))
        d.axhspan(-46.2, 46.2, color=GRID, alpha=0.5, lw=0)
        d.text(xs[0], 41, "for scale: paper's own hold-out max |ΔVth| 46.2 mV", fontsize=6.5, color=INK2, va="top")
        d.set_ylim(-52, 52)
        d.set_xlabel("log10 ND,eff (cm⁻³) of the re-fitted set")
        d.set_ylabel("predicted hold-out Vth − median (mV)")
        d.set_title("(d) hold-out tests cannot tell the sets apart", loc="left")
        d.legend(loc="center", bbox_to_anchor=(0.5, 0.74), fontsize=6, ncol=3)
        style(d)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    print("saved", args.out)


if __name__ == "__main__":
    main()
