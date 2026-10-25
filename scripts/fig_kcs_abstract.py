"""One-row figure for the one-page conference abstract (KCS 2027). Reads stored results only.

The abstract presents an identifiability check added to the agentic TCAD calibration workflow of the reference
paper (Jeong et al., arXiv:2609.12184). The earlier three-panel version (ND,eff profile, Vth spread of the
measured-fit sets, design rule) is in the git history (commit 645f871).

(a) schematic: the workflow of the reference paper (calibration loop, transfer tests) followed by the added
    identifiability check (diagnose, reparameterize, design, resolve) before the descriptors are interpreted
(b) design step, linear check at the twin truth (E19): predictive σ of each candidate test given the calibration
    data, divided by the test noise (curve RMS), for the 32 lateral/bias tests (VD 0.02–3 V at L = 1 µm;
    L = 0.5–10 µm) and vertical-structure tests (t_ox 5–20 nm, channel t_ch = t_IWO 2–10 nm)
(c) resolve step on the authors' open demo (E22): end points of the demo's own controller (trust radius 0.3–0.7,
    three starting gate work functions; 15 runs each) with the demo calibration set and with one t_ox = 14 nm
    curve added: largest held-out VD / L error against the Vth error predicted for a t_ox = 20 nm device

Usage:  python scripts/fig_kcs_abstract.py [--out results/fig_kcs_abstract.png]
"""
from __future__ import annotations

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

BLUE = "#2a78d6"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#d9d8d4"
WASH = "#eaf2fc"  # blue at about 10 % on white


def style(ax):
    ax.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
        ax.spines[s].set_linewidth(0.6)
    ax.tick_params(colors=INK2, labelsize=6.2, width=0.6, length=2.5)


def box(ax, x, y, w, h, text, edge, face="white", lw=0.8, size=5.6, weight="normal", ha="center", pad=0.012):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad={pad},rounding_size=0.025", fc=face, ec=edge,
                                lw=lw, transform=ax.transAxes, clip_on=False))
    tx = x + w / 2 if ha == "center" else x + 0.025
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=size, color=INK, weight=weight, transform=ax.transAxes,
            linespacing=1.25)


def arrow(ax, p0, p1, color=MUTED, style="-|>", rad=0.0, lw=0.8):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=6, color=color, lw=lw,
                                 connectionstyle=f"arc3,rad={rad}", transform=ax.transAxes, clip_on=False))


def panel_workflow(a):
    a.set_axis_off()
    a.set_title("(a) workflow with the added step", loc="left", color=INK, fontsize=6.8)
    x0, w = 0.02, 0.74
    xr = x0 + w + 0.045                      # x of the recalibration line
    box(a, x0, 0.855, w, 0.115, "Agentic calibration loop [1]\nresiduals → sensitivities → update", INK2)
    arrow(a, (x0 + w / 2, 0.84), (x0 + w / 2, 0.79))
    box(a, x0, 0.695, w, 0.085, r"Transfer tests [1]: $V_{\rm D}$, $L$ splits", INK2)
    arrow(a, (x0 + w / 2, 0.68), (x0 + w / 2, 0.635))
    # added step
    a.add_patch(FancyBboxPatch((x0, 0.165), w, 0.455, boxstyle="round,pad=0.012,rounding_size=0.03", fc=WASH,
                               ec=BLUE, lw=1.1, transform=a.transAxes, clip_on=False))
    a.text(x0 + 0.03, 0.583, "Identifiability check (this work)", fontsize=5.6, color=INK, weight="bold",
           transform=a.transAxes, va="center")
    steps = ["1  Diagnose: Fisher, profile",
             "2  Reparameterize: combinations",
             r"3  Design: $\sigma_{\rm pred}$/noise of tests",
             r"4  Resolve: $t_{\rm ox}$, $t_{\rm ch}$ split, C–V"]
    ys = [0.475, 0.38, 0.285, 0.19]
    for y, txt in zip(ys, steps):
        box(a, x0 + 0.03, y, w - 0.06, 0.07, txt, BLUE, lw=0.6, size=5.4, ha="left", pad=0.006)
    # recalibration loop: step 4 -> calibration loop
    y4 = ys[-1] + 0.035
    a.plot([x0 + w - 0.03 + 0.006, xr], [y4, y4], color=BLUE, lw=0.7, transform=a.transAxes, clip_on=False)
    a.plot([xr, xr], [y4, 0.912], color=BLUE, lw=0.7, transform=a.transAxes, clip_on=False)
    arrow(a, (xr, 0.912), (x0 + w + 0.014, 0.912), color=BLUE, lw=0.7)
    a.text(xr + 0.035, (y4 + 0.912) / 2, "recalibrate with split data", rotation=90, fontsize=5.3, color=INK2,
           ha="center", va="center", transform=a.transAxes)
    arrow(a, (x0 + w / 2, 0.15), (x0 + w / 2, 0.105))
    box(a, x0, 0.015, w, 0.085, "Interpret descriptors\n(process insight)", INK2, size=5.4)


def panel_design(b, rng):
    e19 = json.load(open("results/e19/summary.json"))
    tests = {k: v for k, v in e19["tests"].items() if k not in e19["calibration"]}
    groups = [(r"$V_{\rm D}$", lambda k: k.startswith("L1_"), INK2, "o"),
              (r"$L$", lambda k: k.startswith("L") and not k.startswith("L1_"), INK2, "o"),
              (r"$t_{\rm ox}$", lambda k: k.startswith("tox"), BLUE, "s"),
              (r"$t_{\rm ch}$", lambda k: k.startswith("tiwo"), BLUE, "s")]
    n_used = 0
    for j, (lab, sel, col, mk) in enumerate(groups):
        v = np.array([t["rms_ratio"] for k, t in tests.items() if sel(k)])
        n_used += v.size
        b.plot(j + rng.uniform(-0.15, 0.15, v.size), v, mk, color=col, ms=2.8, mec="white", mew=0.3)
        if v.max() < 1:
            b.text(j, v.min() / 1.6, f"{v.min():.2f}–{v.max():.2f}", ha="center", va="top", fontsize=5.6, color=INK2)
        else:
            # one decimal on the lower bound: the smallest ratio (40.6) is the number the text quotes
            b.text(j, v.max() * 1.6, f"{v.min():.1f}–{v.max():.0f}", ha="center", va="bottom", fontsize=5.6, color=INK)
    assert n_used == len(tests) == 32
    b.axhline(1.0, color=INK, ls=":", lw=0.8)
    b.text(3.45, 1.12, "= test noise", fontsize=5.4, color=INK2, ha="right", va="bottom")
    b.set_yscale("log")
    b.set_ylim(0.05, 3000)
    b.set_xticks(range(len(groups)))
    b.set_xticklabels([g[0] for g in groups], fontsize=5.8)
    for t, g in zip(b.get_xticklabels(), groups):
        t.set_color(INK if g[2] == BLUE else INK2)
    tr = matplotlib.transforms.blended_transform_factory(b.transData, b.transAxes)
    for x0, x1, lab in ((0, 1, "lateral, as in [1]"), (2, 3, "vertical (added)")):
        b.plot([x0 - 0.3, x1 + 0.3], [-0.10, -0.10], color=MUTED, lw=0.5, transform=tr, clip_on=False)
        b.text((x0 + x1) / 2, -0.115, lab, transform=tr, ha="center", va="top", fontsize=5.4,
               color=INK if x0 == 2 else INK2)
    b.set_xlim(-0.55, len(groups) - 0.45)
    b.set_ylabel(r"predictive $\sigma$ / test noise")
    b.set_title("(b) design: 32 tests (DEVSIM twin)", loc="left", color=INK, fontsize=6.8)
    style(b)
    b.grid(False, axis="x")


def panel_resolve(c):
    runs = json.load(open("results/e22/summary.json"))["D_controller_runs"]
    spec = (("demo set", INK2, "o", "demo calibration data"),
            ("+ t_ox 14 nm curve", BLUE, "s", r"+ one $t_{\rm ox}$ = 14 nm curve"))
    for tag, col, mk, lab in spec:
        rs = [x for x in runs if x["calibration"] == tag]
        assert len(rs) == 15
        hx = [max(x["heldout_drain_max_dVth_mV"], x["heldout_length_max_dVth_mV"]) for x in rs]
        hy = [x["t_ox 20 nm dVth vs truth (mV)"] for x in rs]
        c.scatter(hx, hy, s=15, marker=mk, color=col, edgecolors="white", linewidths=0.6, label=lab, zorder=3)
        qf = [x["Q_f"] / 1e12 for x in rs]
        if tag == "demo set":
            # the largest end point sits on the demo's upper bound of Q_f (1e13 cm^-2)
            assert abs(max(qf) - 10.0) < 1e-9
            c.text(15.6, -400, f"end-point $Q_{{\\rm f}}$\n{min(qf):.1f}–{max(qf):.1f}\n(10 = upper bound)", fontsize=5.4,
                   color=INK2, ha="right", va="center")
        else:
            c.text(12.35, -60, f"$Q_{{\\rm f}}$ {min(qf):.2f}–{max(qf):.2f}", fontsize=5.4, color=INK, ha="center",
                   va="top")
    c.text(10.7, -230, "$Q_{\\rm f}$ in 10$^{12}$ cm$^{-2}$\n(truth 6)", fontsize=5.3, color=INK2, ha="left", va="center")
    c.axhline(0.0, color=MUTED, lw=0.6, zorder=1)
    c.set_xlim(10.5, 18.5)
    c.set_ylim(-560, 80)
    c.set_xlabel(r"largest held-out $|\Delta V_{\rm th}|$, $V_{\rm D}$ and $L$ (mV)")
    c.set_ylabel(r"$V_{\rm th}$ error, $t_{\rm ox}$ = 20 nm device (mV)")
    c.set_title("(c) resolve: open demo [6], 15 runs", loc="left", color=INK, fontsize=6.8)
    c.legend(loc="lower left", frameon=False, handletextpad=0.3, borderaxespad=0.2)
    style(c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/fig_kcs_abstract.png")
    args = ap.parse_args()
    plt.rcParams.update({"font.size": 6.5, "axes.titlesize": 6.8, "legend.fontsize": 5.6, "axes.labelcolor": INK2,
                         "axes.labelsize": 6.4})
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(6.9, 2.5), gridspec_kw={"width_ratios": [1.12, 0.92, 0.98]})
    rng = np.random.default_rng(3)
    panel_workflow(a)
    panel_design(b, rng)
    panel_resolve(c)
    fig.tight_layout(pad=0.3, w_pad=0.9)
    fig.savefig(args.out, dpi=300, metadata={"Software": None})
    print("saved", args.out)


if __name__ == "__main__":
    main()
