#!/usr/bin/env python
"""E22 figure: ambiguity analysis on the open demo of the reference paper. Reads results/e22/summary.json only.

(a) profile likelihood of Q_f (other four parameters re-fitted) for the demo's calibration set and
    with one t_ox = 14 nm curve added; sigma_eff of each best fit (model-form error included);
    symmetric-log y axis, linear below 1
(b) end points of the demo's own controller (trust radius 0.3-0.7, three starting gate work
    functions each): largest held-out V_D / L error against the V_th error predicted for a
    t_ox = 20 nm device (held out in both cases)
(c) profile likelihood of N_tail when the interface-trap term C_it is a free parameter

Usage:  python scripts/fig_e22_demo.py [results/e22/summary.json] [results/e22/fig_e22_demo.png]
"""

from __future__ import annotations

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#ffffff"
YLIM = (-0.15, 2e5)


def style(ax):
    ax.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
        ax.spines[s].set_linewidth(0.6)
    ax.tick_params(colors=INK2, labelsize=6.2, width=0.6, length=2.5)


def main() -> None:
    src = sys.argv[1] if len(sys.argv) > 1 else "results/e22/summary.json"
    out = sys.argv[2] if len(sys.argv) > 2 else "results/e22/fig_e22_demo.png"
    r = json.load(open(src))
    plt.rcParams.update({"font.size": 6.5, "axes.titlesize": 6.8, "legend.fontsize": 5.6,
                         "axes.labelcolor": INK2, "axes.labelsize": 6.4})
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.0, 2.45), gridspec_kw={"width_ratios": [1.0, 1.0, 1.0]})
    truth = r["default_run"]["truth_params"]

    # (a) Q_f profile
    prof = r["C_profile_Q_f"]
    for (tag, p), col, mk, lab in zip(prof.items(), (BLUE, ORANGE), ("o", "s"),
                                      ("demo set", "+ t_ox 14 nm curve")):
        x = np.array([row["Q_f"] for row in p["rows"]]) / 1e12
        y = np.array([row["dchi2_eff"] for row in p["rows"]])
        keep = x >= 1.9
        a.plot(x[keep], y[keep], color=col, lw=1.2, marker=mk, ms=2.4, mfc=col, mec=SURF, mew=0.5,
               label=lab, zorder=3)
    a.set_yscale("symlog", linthresh=1.0, linscale=0.6)
    a.axhline(3.84, color=MUTED, lw=0.6, zorder=1)
    a.text(2.1, 2.1, "Δχ² = 3.84 (95 %)", color=INK2, fontsize=5.4, va="center")
    a.axvline(truth["Q_f"] / 1e12, color=INK2, lw=0.6, zorder=1)
    a.text(truth["Q_f"] / 1e12 + 0.12, 1.6, "truth", color=INK2, fontsize=5.4, ha="left")
    qd = r["default_run"]["final_params"]["Q_f"] / 1e12
    a.plot([qd], [0.0], marker="v", ms=4.5, color=INK, mec=SURF, mew=0.5, zorder=4, clip_on=False)
    a.text(qd, 0.45, "demo end point", fontsize=5.4, color=INK2, ha="center")
    a.set_xlim(1.9, 10.3)
    a.set_ylim(*YLIM)
    a.set_xlabel("Q$_f$ (10$^{12}$ cm$^{-2}$), fixed")
    a.set_ylabel("Δχ²$_{eff}$ (symlog)")
    a.set_title("(a) profile of Q$_f$", loc="left", color=INK)
    a.legend(loc="lower right", bbox_to_anchor=(1.08, 0.22), frameon=False, handlelength=1.3, handletextpad=0.4)
    style(a)

    # (b) controller end points
    runs = r["D_controller_runs"]
    for tag, col, mk, lab in (("demo set", BLUE, "o", "demo set"), ("+ t_ox 14 nm curve", ORANGE, "s",
                                                                    "+ t_ox 14 nm curve")):
        rs = [x for x in runs if x["calibration"] == tag]
        hx = [max(x["heldout_drain_max_dVth_mV"], x["heldout_length_max_dVth_mV"]) for x in rs]
        hy = [x["t_ox 20 nm dVth vs truth (mV)"] for x in rs]
        b.scatter(hx, hy, s=16, marker=mk, color=col, edgecolors=SURF, linewidths=0.8, label=lab, zorder=3)
    b.axhline(0.0, color=MUTED, lw=0.6, zorder=1)
    seen = set()
    for x in (x for x in runs if x["calibration"] == "demo set"):
        hxv = max(x["heldout_drain_max_dVth_mV"], x["heldout_length_max_dVth_mV"])
        hyv = x["t_ox 20 nm dVth vs truth (mV)"]
        key = (round(hxv, 2), round(hyv, 1))
        if key in seen:
            continue
        seen.add(key)
        right = x["trust_radius"] == 0.7
        b.text(hxv + (0.18 if right else -0.18), hyv, f"{x['trust_radius']:g}", fontsize=5.4, color=INK2,
               ha="left" if right else "right", va="center")
    b.text(10.2, -470, "numbers: trust radius\n(3 starting φ$_G$ each)", fontsize=5.4, color=INK2, va="top")
    b.set_xlim(10.0, 18.5)
    b.set_ylim(-560, 90)
    b.set_xlabel("largest held-out |ΔV$_{th}$|, V$_D$ and L tests (mV)")
    b.set_ylabel("V$_{th}$ error, t$_{ox}$ 20 nm device (mV)")
    b.set_title("(b) demo controller, 15 runs each", loc="left", color=INK)
    b.legend(loc="center left", bbox_to_anchor=(0.0, 0.47), frameon=False, handletextpad=0.3)
    style(b)

    # (c) N_tail profile with a free interface term
    pn = r["E_cit_vs_tail"]["profile_N_tail"]
    spec_c = (("demo set", BLUE, "-", "o", "demo set"), ("+ t_ox 14 nm curve", ORANGE, "--", "s", "+ t_ox 14 nm curve"),
              ("+ t_ch 7 nm curve", AQUA, "-", "^", "+ t_ch 7 nm curve"))
    for tag, col, ls, mk, lab in spec_c:
        p = pn[tag]
        x = np.array([row["log10_N_tail"] for row in p["rows"]])
        y = np.array([row["dchi2_eff"] for row in p["rows"]])
        c.plot(x, y, color=col, lw=1.2, ls=ls, marker=mk, ms=2.4, mfc=col, mec=SURF, mew=0.5,
               label=lab, zorder=3 if tag != "+ t_ox 14 nm curve" else 4)
    c.set_yscale("symlog", linthresh=1.0, linscale=0.6)
    c.axhline(3.84, color=MUTED, lw=0.6, zorder=1)
    c.axvline(np.log10(truth["N_tail"]), color=INK2, lw=0.6, zorder=1)
    c.text(np.log10(truth["N_tail"]) - 0.06, 6e4, "truth", color=INK2, fontsize=5.4, ha="right")
    c.set_xlim(15.9, 20.05)
    c.set_ylim(*YLIM)
    c.set_xlabel("log$_{10}$ N$_{tail}$ (cm$^{-3}$ eV$^{-1}$), fixed")
    c.set_ylabel("Δχ²$_{eff}$ (symlog)")
    c.set_title("(c) N$_{tail}$ when C$_{it}$ is also fitted", loc="left", color=INK)
    c.legend(loc="lower left", bbox_to_anchor=(0.0, 0.2), frameon=False, handlelength=1.8)
    style(c)

    fig.tight_layout(pad=0.4, w_pad=1.0)
    fig.savefig(out, dpi=300, metadata={"Software": None})
    print("wrote", out)


if __name__ == "__main__":
    main()
