"""Figure for the core claim: passing the transfer tests does not pin down the parameters.

(a) two parameter sets 160x apart in ND,eff, both fitted to the measured calibration curves: same
    calibration curves and same transfer-test curve (L = 2 µm)
(b) profile of ND,eff: Δχ² of the calibration fit (measured data, σ_eff; synthetic twin, σ = 0.02 dec)
(c) the same measured-fit sets: spread of their predicted Vth (each test centred on the median set) for
    the paper's five transfer tests and for the two recommended tests (thicker channel / gate oxide);
    Vth from the full simulated curves (the digitised transfer curves have gaps around 100 nA)
(d) linear check at the twin truth: predictive 1σ of each test given the calibration data, in units
    of the test's own noise (< 1: in the linear approximation the calibration-consistent sets agree
    within the test's noise)

Reads results/e5_nd, e11_2W_v3, e11_2W_v2, e12_nd, e13, e14 (no simulation).
Usage:  python scripts/fig_core_claim.py [--best results/e11_2W_v3/summary.json,results/e11_2W_v2/summary.json]
            [--out results/fig_core_claim.png]
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
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#d9d8d4"
ND_PAPER_2W, ND_PAPER_4W = np.log10(2.56e18), np.log10(6.39e16)
TRANSFER = ["L1_VD0.1", "L1_VD0.2", "L1_VD0.5", "L2_VD0.05", "L2_VD1"]


def style(ax):
    ax.grid(True, color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=7.5)


def load_meas(path="data/meas.csv"):
    d = collections.defaultdict(list)
    for r in csv.DictReader(open(path)):
        d[(r["set"], float(r["L_um"]), float(r["vd"]))].append((float(r["vg"]), float(r["id"])))
    return {k: np.array(sorted(v)) for k, v in d.items()}


def measured_sets(best_paths, profile="results/e12_nd"):
    sets = [json.load(open(f)) for f in best_paths.split(",")]
    sets += [json.load(open(f)) for f in sorted(glob.glob(f"{profile}/log10_nd_eff_*/summary.json"))]
    return sorted(sets, key=lambda s: s["theta_hat"]["log10_nd_eff"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--best", default="results/e11_2W_v3/summary.json,results/e11_2W_v2/summary.json")
    ap.add_argument("--out", default="results/fig_core_claim.png")
    args = ap.parse_args()
    meas = load_meas()
    sets = measured_sets(args.best)
    best = min((json.load(open(f)) for f in args.best.split(",")), key=lambda s: s["chi2"])
    nd = np.array([s["theta_hat"]["log10_nd_eff"] for s in sets])
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "legend.fontsize": 6.8, "axes.labelcolor": INK2})
    fig, ax = plt.subplots(2, 2, figsize=(7.6, 6.3))
    a, b, c, d = ax.ravel()

    # (a) curves of the two extreme sets
    lo, hi = sets[0], sets[-1]
    vg = np.array(lo["vg"])
    for key, curve_of, mk in (((1.0, 1.0), lambda s: s["curves"]["1.0"], "s"),
                              ((2.0, 1.0), lambda s: s["holdout"]["L2_VD1"]["curve"], "^"),
                              ((1.0, 0.05), lambda s: s["curves"]["0.05"], "s")):
        m = meas[("2W", key[0], key[1])]
        a.semilogy(m[:, 0], m[:, 1], mk, mfc="none", mec=ORANGE, ms=3.2, mew=0.8)
        a.semilogy(vg, curve_of(lo), "-", color=INK, lw=1.4)
        a.semilogy(vg, curve_of(hi), "--", color=AQUA, lw=1.4)
    a.text(1.55, 9e-4, "VD = 1 V", fontsize=6.5, color=INK2, va="center")
    a.text(1.55, 7e-6, "VD = 0.05 V", fontsize=6.5, color=INK2, va="center")
    a.plot([], [], "s", mfc="none", mec=ORANGE, ms=3.2, label="measured, L = 1 µm (calibration)")
    a.plot([], [], "^", mfc="none", mec=ORANGE, ms=3.2, label="measured, L = 2 µm (transfer test)")
    a.plot([], [], "-", color=INK, lw=1.4, label=f"fit, ND,eff = $10^{{{lo['theta_hat']['log10_nd_eff']:.2f}}}$ cm$^{{-3}}$")
    a.plot([], [], "--", color=AQUA, lw=1.4, label=f"fit, ND,eff = $10^{{{hi['theta_hat']['log10_nd_eff']:.2f}}}$ cm$^{{-3}}$")
    a.set_ylim(1e-11, 3e-3)
    a.set_xlim(-1.0, 2.0)
    a.set_xlabel("VG (V)")
    a.set_ylabel("ID (A)")
    a.set_title("(a) 160× different ND,eff, same measured curves", loc="left", color=INK)
    a.legend(loc="lower right", frameon=False)
    style(a)

    # (b) profiles: fixed-ND,eff fits as a line, the free fit (best) as a star
    scale = (best["sigma_eff"] / best["sigma_dec"]) ** 2
    chi_min = min(s["chi2"] for s in sets)
    prof = sorted((s for s in sets if "log10_nd_eff" in s.get("fixed", {})), key=lambda s: s["fixed"]["log10_nd_eff"])
    b.plot([s["fixed"]["log10_nd_eff"] for s in prof], [(s["chi2"] - chi_min) / scale for s in prof], "o-",
           color=ORANGE, lw=1.8, ms=4.5)
    b.plot([best["theta_hat"]["log10_nd_eff"]], [(best["chi2"] - chi_min) / scale], "*", color=ORANGE, ms=10, mec=INK, mew=0.5)
    b.text(best["theta_hat"]["log10_nd_eff"], 0.3, "free fit", fontsize=6.2, color=INK2, ha="center", va="bottom")
    b.text(16.9, 1.6, "measured data\n(σ_eff, model error incl.)", fontsize=6.5, color=INK2, ha="left")
    e5 = json.load(open("results/e5_nd/summary.json"))
    tw = sorted((r["value"], r["chi2"]) for r in e5["records"] if "chi2" in r)
    tx, ty = np.array([t[0] for t in tw]), np.array([t[1] for t in tw])
    b.plot(tx, np.minimum(ty, 6.0), "s-", color=AQUA, lw=1.4, ms=3.8)
    b.text(15.85, 0.3, "synthetic twin\n(σ = 0.02 dec)", fontsize=6.5, color=INK2, ha="left")
    if ty.max() > 6:
        b.annotate(f"χ² = {ty.max():.0f}", (tx[np.argmax(ty)], 6.0), (tx[np.argmax(ty)] - 0.6, 5.2), fontsize=6.3, color=INK2,
                   arrowprops={"arrowstyle": "->", "color": MUTED, "lw": 0.8})
    b.axhline(3.84, color=INK2, ls=":", lw=1)
    b.text(17.0, 3.95, "95 % threshold (3.84)", fontsize=6.3, color=INK2, va="bottom")
    for x, lab in ((ND_PAPER_2W, "paper\n2%-W"), (ND_PAPER_4W, "paper\n4%-W")):
        b.axvline(x, color=MUTED, ls="--", lw=0.7)
        b.text(x - 0.05, 5.95, lab, fontsize=6.2, color=INK2, ha="right", va="top")
    b.set_ylim(-0.25, 6.2)
    b.set_xlim(15.75, 19.55)
    b.set_xlabel("log10 ND,eff (cm⁻³), fixed; other 8 parameters re-fitted")
    b.set_ylabel("Δχ² of the calibration fit")
    b.set_title("(b) calibration data cannot fix ND,eff", loc="left", color=INK)
    style(b)

    # (c) spread of the predicted Vth across all measured-fit sets, per test (centred on the median set)
    cats, cols_c, vals_c = [], [], []
    for k, lab in zip(TRANSFER, ["VD\n0.1 V", "VD\n0.2 V", "VD\n0.5 V", "L 2 µm\n0.05 V", "L 2 µm\n1 V"]):
        cats.append(lab); cols_c.append(INK2)
        vals_c.append(np.array([1e3 * extract(np.array(s["vg"]), np.array(s["holdout"][k]["curve"]))["vth"] for s in sets]))
    e13p = "results/e13/summary.json"
    if os.path.exists(e13p):
        e13 = json.load(open(e13p))
        rows = e13["rows"]
        for cond, lab in (("t_ox 14 nm", "t_ox\n14 nm"), ("t_IWO 7 nm", "t_IWO\n7 nm")):
            v = np.array([rows[k][cond]["1"]["dVth_mV"] for k in rows], float)
            if np.all(np.isfinite(v)):
                cats.append(lab); cols_c.append(BLUE); vals_c.append(v)
    rng = np.random.default_rng(3)
    for j, (v, col) in enumerate(zip(vals_c, cols_c)):
        y = v - np.median(v)
        c.plot(j + rng.uniform(-0.12, 0.12, y.size), y, "o", color=col, ms=3.6, alpha=0.9, mec="none")
        sp = y.max() - y.min()
        c.text(j, y.max() + 12, f"{sp:.0f} mV", ha="center", va="bottom", fontsize=6.3, color=INK if col == BLUE else INK2)
    c.axhspan(-46.2, 46.2, color=GRID, alpha=0.35, lw=0)
    c.text(-0.5, -42, "paper's transfer-test error (±46 mV)", fontsize=6.2, color=INK2, ha="left", va="bottom")
    c.axhline(0, color=MUTED, lw=0.6)
    c.set_xticks(range(len(cats)))
    c.set_xticklabels(cats, fontsize=6.4)
    for lab, col in zip(c.get_xticklabels(), cols_c):
        lab.set_color(INK if col == BLUE else INK2)
    c.set_xlim(-0.6, len(cats) - 0.4)
    c.set_ylabel("predicted Vth − median over the sets (mV)")
    if os.path.exists(e13p):
        c.text(0.02, 0.03, f"{len(sets)} sets fitting the measured calibration\ncurves (Δχ²_eff ≤ "
               f"{max((s['chi2'] - chi_min) / scale for s in sets):.1f}); new tests at VD = 1 V;\n"
               f"C–V differs by up to {e13['tests']['C-V']['max_dev_pct']:.1f} % of C$_{{ox}}$",
               transform=c.transAxes, fontsize=5.8, color=INK2, va="bottom")
    lo_, hi_ = min(np.min(v - np.median(v)) for v in vals_c), max(np.max(v - np.median(v)) for v in vals_c)
    c.set_ylim(lo_ - 30, hi_ + 45)
    c.set_title("(c) same sets: transfer tests agree, recommended tests do not", loc="left", color=INK)
    style(c)
    c.grid(False, axis="x")

    # (d) linear predictive spread
    e14 = json.load(open("results/e14/summary.json"))["linear"]
    names = list(e14)
    vals = [e14[n]["rms_ratio"] for n in names]
    cols = [INK2 if "(paper)" in n else BLUE for n in names]
    yy = np.arange(len(names))[::-1]
    d.barh(yy, vals, color=cols, height=0.62)
    for y_, v_, n in zip(yy, vals, names):
        d.text(v_ * 1.15, y_, f"{v_:.2f}" if v_ < 1 else f"{v_:.0f}", va="center", fontsize=6.3, color=INK2)
    d.set_yticks(yy)
    d.set_yticklabels([n.replace(" (paper)", "") for n in names], fontsize=6.6)
    for lab, col in zip(d.get_yticklabels(), cols):
        lab.set_color(INK if col == BLUE else INK2)
    d.set_xscale("log")
    d.set_xlim(0.1, 2000)
    d.axvline(1.0, color=INK, lw=0.9, ls=":")
    d.text(1.1, yy[0] + 0.6, "= noise", fontsize=6.3, color=INK2)
    d.text(0.98, 0.98, "grey: paper's transfer tests\nblue: recommended tests", transform=d.transAxes,
           fontsize=6.3, color=INK2, ha="right", va="top")
    d.set_xlabel("predictive σ given calibration data / test noise")
    d.set_title("(d) test informativeness (twin, linear)", loc="left", color=INK)
    style(d)
    d.grid(False, axis="y")

    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    print("saved", args.out)


if __name__ == "__main__":
    main()
