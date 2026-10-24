"""One-row figure for the one-page conference abstract (KCS 2027). Reads stored results only.

(a) profile of ND,eff: Δχ² of the calibration fit with ND,eff fixed and the other eight parameters re-fitted,
    for the synthetic twin (E5: two calibration curves; E18: five drain-bias curves; σ = 0.02 dec) and the
    digitised measured curves (E12: two curves; E18: five curves; σ_eff of each free fit, model error included)
(b) the seven sets that fit the measured calibration curves equally well (E12 profile points and the free fits
    v2, v3): predicted Vth for the paper's transfer tests (full simulated curves, the digitised curves have gaps
    around 100 nA) and for the thicker-oxide / thicker-channel devices (E13), each centred on its median set
(c) linear check at the twin truth (E19): predictive σ of each test given the calibration data, divided by the
    test noise (curve RMS), for the 32 lateral/bias tests (VD 0.02–3 V at L = 1 µm; L = 0.5–10 µm) and
    vertical-structure tests (t_ox 5–20 nm, t_IWO 2–10 nm), each at VD = 0.05 and/or 1 V

Usage:  python scripts/fig_kcs_abstract.py [--out results/fig_kcs_abstract.png]
"""
from __future__ import annotations

import argparse
import glob
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from _common import ROOT  # noqa: F401  (puts the repository root on sys.path)
from tftsim.metrics import extract

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#d9d8d4"
ND_PAPER_2W, ND_PAPER_4W = np.log10(2.56e18), np.log10(6.39e16)
TRANSFER = [("L1_VD0.1", "0.1"), ("L1_VD0.2", "0.2"), ("L1_VD0.5", "0.5"), ("L2_VD0.05", "0.05"), ("L2_VD1", "1")]
YMAX = 7.4


def style(ax):
    ax.grid(True, color=GRID, lw=0.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
        ax.spines[s].set_linewidth(0.6)
    ax.tick_params(colors=INK2, labelsize=6.2, width=0.6, length=2.5)


def twin_profile(path):
    s = json.load(open(path))
    pts = sorted((r["value"], r["chi2"]) for r in s["records"] if "chi2" in r)
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])


def measured_profile(path):
    s = json.load(open(path))
    pts = sorted((r["value"], r["dchi2_eff"]) for r in s["rows"])
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])


def measured_sets(best="results/e11_2W_v3/summary.json,results/e11_2W_v2/summary.json", profile="results/e12_nd"):
    sets = [json.load(open(f)) for f in best.split(",")]
    sets += [json.load(open(f)) for f in sorted(glob.glob(f"{profile}/log10_nd_eff_*/summary.json"))]
    return sorted(sets, key=lambda s: s["theta_hat"]["log10_nd_eff"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/fig_kcs_abstract.png")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": 6.5, "axes.titlesize": 6.8, "legend.fontsize": 5.6, "axes.labelcolor": INK2,
                         "axes.labelsize": 6.4})
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(6.9, 2.5), gridspec_kw={"width_ratios": [0.97, 1.13, 0.98]})

    # (a) profiles
    series = [(twin_profile("results/e5_nd/summary.json"), AQUA, "-", "o", "twin, 2 curves (σ = 0.02 dec)"),
              (twin_profile("results/e18_joint_twin/summary.json"), AQUA, "--", "o", "twin, 5 VD curves"),
              (measured_profile("results/e12_nd/summary.json"), ORANGE, "-", "s", "measured, 2 curves (σ_eff)"),
              (measured_profile("results/e18_joint_meas/profile/summary.json"), ORANGE, "--", "s", "measured, 5 VD curves")]
    for (x, y), col, ls, mk, lab in series:
        filled = ls == "-"
        a.plot(x, np.minimum(y, YMAX), ls=ls, marker=mk, color=col, lw=1.1, ms=3.0, label=lab,
               mfc=col if filled else "white", mec=col, mew=0.8)
        if y.max() > YMAX:
            k = int(np.argmax(y))
            a.annotate(f"{y[k]:.0f}", (x[k], YMAX), (18.97, 5.25), fontsize=5.6, color=INK2, ha="right",
                       arrowprops={"arrowstyle": "->", "color": MUTED, "lw": 0.6})
    a.axhline(3.84, color=INK2, ls=":", lw=0.8)
    a.text(18.33, 3.95, "95 % threshold", fontsize=5.6, color=INK2, ha="right", va="bottom")
    for xv, lab, ha, dx in ((ND_PAPER_4W, "paper\n4%-W", "left", 0.05), (ND_PAPER_2W, "paper\n2%-W", "right", -0.05)):
        a.axvline(xv, color=MUTED, ls=(0, (3, 2)), lw=0.5)
        a.text(xv + dx, 2.75, lab, fontsize=5.3, color=INK2, ha=ha, va="bottom")
    a.set_xlim(15.75, 19.55)
    a.set_ylim(-0.25, YMAX + 0.1)
    a.set_xlabel("log10 ND,eff (cm⁻³), fixed")
    a.set_ylabel("Δχ² of the calibration fit")
    a.set_title("(a) calibration data do not fix ND,eff", loc="left", color=INK)
    leg = a.legend(loc="upper left", frameon=True, handlelength=2.2, borderaxespad=0.3, borderpad=0.3, labelspacing=0.3)
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_edgecolor("none")
    leg.get_frame().set_alpha(1.0)
    style(a)

    # (b) predicted Vth of the measured-fit sets
    sets = measured_sets()
    e13 = json.load(open("results/e13/summary.json"))["rows"]
    cats, cols, vals = [], [], []
    for k, lab in TRANSFER:
        cats.append(lab); cols.append(INK2)
        vals.append(np.array([1e3 * extract(np.array(s["vg"]), np.array(s["holdout"][k]["curve"]))["vth"] for s in sets]))
    for cond, lab in (("t_ox 14 nm", "t_ox\n14 nm"), ("t_IWO 7 nm", "t_IWO\n7 nm")):
        cats.append(lab); cols.append(BLUE)
        vals.append(np.array([e13[k][cond]["1"]["dVth_mV"] for k in e13], float))
    rng = np.random.default_rng(3)
    for j, (v, col) in enumerate(zip(vals, cols)):
        y = v - np.median(v)
        b.plot(j + rng.uniform(-0.13, 0.13, y.size), y, "o", color=col, ms=2.6, mec="white", mew=0.3)
        b.text(j, y.max() + 10, f"{y.max() - y.min():.0f} mV", ha="center", va="bottom", fontsize=5.6,
               color=INK if col == BLUE else INK2)
    b.axhspan(-46.2, 46.2, color=GRID, alpha=0.45, lw=0)
    b.text(-0.55, -52, "band: paper's transfer-test\nerror (±46 mV)", fontsize=5.4, color=INK2, ha="left", va="top")
    b.axhline(0, color=MUTED, lw=0.5)
    b.set_xticks(range(len(cats)))
    b.set_xticklabels(cats, fontsize=5.2)
    for t, col in zip(b.get_xticklabels(), cols):
        t.set_color(INK if col == BLUE else INK2)
    tr = matplotlib.transforms.blended_transform_factory(b.transData, b.transAxes)
    for x0, x1, lab in ((0, 2, "L 1 µm, VD (V)"), (3, 4, "L 2 µm")):
        b.plot([x0 - 0.3, x1 + 0.3], [-0.10, -0.10], color=MUTED, lw=0.5, transform=tr, clip_on=False)
        b.text((x0 + x1) / 2, -0.115, lab, transform=tr, ha="center", va="top", fontsize=5.2, color=INK2)
    b.set_xlim(-0.6, len(cats) - 0.4)
    lo = min(np.min(v - np.median(v)) for v in vals)
    hi = max(np.max(v - np.median(v)) for v in vals)
    b.set_ylim(lo - 25, hi + 45)
    b.set_ylabel("predicted Vth − median (mV)")
    b.set_title(f"(b) {len(sets)} equivalent measured-fit sets", loc="left", color=INK)
    b.text(0.03, 0.97, "grey: paper's transfer tests\nblue: thicker oxide / channel (VD 1 V)", transform=b.transAxes,
           fontsize=5.4, color=INK2, va="top")
    style(b)
    b.grid(False, axis="x")

    # (c) predictive sigma / noise for the E19 tests
    e19 = json.load(open("results/e19/summary.json"))
    tests = {k: v for k, v in e19["tests"].items() if k not in e19["calibration"]}
    groups = [("VD", lambda k: k.startswith("L1_"), INK2),
              ("L", lambda k: k.startswith("L") and not k.startswith("L1_"), INK2),
              ("t_ox", lambda k: k.startswith("tox"), BLUE),
              ("t_IWO", lambda k: k.startswith("tiwo"), BLUE)]
    n_used = 0
    for j, (lab, sel, col) in enumerate(groups):
        v = np.array([t["rms_ratio"] for k, t in tests.items() if sel(k)])
        n_used += v.size
        c.plot(j + rng.uniform(-0.15, 0.15, v.size), v, "o", color=col, ms=2.8, mec="white", mew=0.3)
        if v.max() < 1:
            c.text(j, v.min() / 1.6, f"{v.min():.2f}–{v.max():.2f}", ha="center", va="top", fontsize=5.6, color=INK2)
        else:
            c.text(j, v.max() * 1.6, f"{v.min():.0f}–{v.max():.0f}", ha="center", va="bottom", fontsize=5.6, color=INK)
    assert n_used == len(tests)
    c.axhline(1.0, color=INK, ls=":", lw=0.8)
    c.text(3.45, 1.12, "= test noise", fontsize=5.4, color=INK2, ha="right", va="bottom")
    c.set_yscale("log")
    c.set_ylim(0.05, 3000)
    c.set_xticks(range(len(groups)))
    c.set_xticklabels([g[0] for g in groups], fontsize=5.6)
    for t, g in zip(c.get_xticklabels(), groups):
        t.set_color(INK if g[2] == BLUE else INK2)
    c.set_xlim(-0.55, len(groups) - 0.45)
    c.set_ylabel("predictive σ / test noise")
    c.set_title(f"(c) {len(tests)} tests: lateral vs vertical (twin)", loc="left", color=INK)
    style(c)
    c.grid(False, axis="x")

    fig.tight_layout(pad=0.3, w_pad=0.8)
    fig.savefig(args.out, dpi=300)
    print("saved", args.out)


if __name__ == "__main__":
    main()
