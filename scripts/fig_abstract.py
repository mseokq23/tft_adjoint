"""Draft figure for the KCS 2027 abstract (reads finished results; no simulation).

(a) Profile likelihood of ND,eff: RMS deviation of the re-fitted calibration curves
    and of the predicted paper hold-out curves, against the digitisation noise.
(b) Calibration curves of parameter sets that differ ~30x in ND,eff.
(c) 1-sigma of ND,eff and the tail density N_TA after adding one measurement type.
(d) DOS detectability vs energy (acceptor / donor) and the Fermi-level window.

Usage:  python scripts/fig_abstract.py [--e5 results/e5_nd] [--out results/fig_abstract.png]
"""
from __future__ import annotations

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"


def style(ax):
    ax.grid(True, color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--e5", default="results/e5_nd")
    ap.add_argument("--out", default="results/fig_abstract.png")
    args = ap.parse_args()
    e5 = json.load(open(f"{args.e5}/summary.json"))
    e6 = json.load(open("results/e6/summary.json"))
    e7 = json.load(open("results/e7/summary.json"))
    sigma = e5["sigma_dec"]
    recs = e5["records"]
    v_true = e5["v_true"]

    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "legend.fontsize": 7})
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.8))
    a, b, c, d = ax.ravel()

    # (a) profile: RMS deviation (one unit: decades) of the fit and of the hold-out prediction
    off = [r for r in recs if r["chi2"] > 0.0]          # the truth itself has zero deviation by definition
    x = np.array([r["value"] for r in off])
    fit = np.array([max(r["rms_dec"], 1e-5) for r in off])
    ho = [(r["value"], max(max(q["rms_dec"] for q in r["holdout"].values()), 1e-5))
          for r in off if isinstance(r.get("holdout"), dict) and "error" not in r["holdout"]]
    a.semilogy(x, fit, "o-", color=BLUE, lw=2, ms=5, label="calibration curves (re-fitted)")
    if ho:
        a.semilogy([h[0] for h in ho], [h[1] for h in ho], "s--", color=ORANGE, lw=2, ms=5,
                   label="paper hold-out curves (predicted, worst of 5)")
    a.axhline(sigma, color=INK2, lw=1, ls=":")
    a.text(x.min(), sigma * 1.15, f"digitisation noise σ = {sigma} dec", color=INK2, fontsize=7,
           bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
    a.axvline(v_true, color=INK2, lw=0.8)
    a.text(v_true + 0.04, 2.2e-5, "2%-W fit\n(paper)", color=INK2, fontsize=6.5, va="center")
    w4 = np.log10(6.39e16)
    a.axvline(w4, color=INK2, lw=0.8, ls="--")
    a.text(w4 + 0.04, 2.2e-5, "4%-W fit\n(paper)", color=INK2, fontsize=6.5, va="center")
    a.set_xlim(min(x.min(), w4) - 0.25, x.max() + 0.15)
    a.set_xlabel("log10 ND,eff (cm⁻³), held fixed; other 8 parameters re-fitted")
    a.set_ylabel("RMS deviation from the truth (dec)")
    a.set_ylim(1e-5, 0.2)
    a.set_title("(a) ND,eff profile: fit and hold-out prediction", loc="left")
    a.legend(loc="upper left")
    style(a)

    # (b) calibration curves of the two most distant indistinguishable sets
    vg = np.array(e5["vg"])
    good = sorted([r for r in recs if r["chi2"] < 3.84], key=lambda r: r["value"])
    pick = [good[0], next(r for r in good if r["chi2"] == 0.0), good[-1]]
    styles = [(ORANGE, "--"), (INK, "-"), (BLUE, ":")]
    for (col, ls), r in zip(styles, pick):
        key = min(e5["curves"], key=lambda k: abs(float(k) - r["value"]))
        th = r["theta"]
        lab = f"ND {10 ** r['value']:.1e}, WF {th['wf_gate']:.3f} eV"
        for vd, i in e5["curves"][key].items():
            b.semilogy(vg, np.maximum(np.array(i), 1e-16), ls=ls, color=col, lw=1.6 if col == INK else 2.0,
                       label=lab if float(vd) == 1.0 else None)
    b.axhline(1e-11, color=INK2, lw=0.6, ls=":")
    b.text(2.45, 1.0e-3, "VD = 1 V", ha="right", fontsize=7, color=INK2)
    b.text(2.45, 5.0e-6, "VD = 0.05 V", ha="right", fontsize=7, color=INK2)
    b.set_ylim(1e-12, 3e-3)
    b.set_xlabel("VG (V)")
    b.set_ylabel("ID (A)")
    b.set_title(f"(b) {10 ** (pick[-1]['value'] - pick[0]['value']):.0f}× different ND,eff, same curves", loc="left")
    b.legend(loc="lower right", fontsize=6.5, title="re-fitted sets (L = 1 µm)", title_fontsize=6.5)
    style(b)

    # (c) 1-sigma of ND,eff and N_TA (decades) after adding one measurement type
    base = e7["base"]["std_phys"]
    cand = e7["candidates"]
    rows = [("calibration only", base), ("+ paper hold-outs (5)", cand["all paper hold-outs"]["std_phys"]),
            ("+ t_ox 14 nm", cand["tox14"]["std_phys"]), ("+ t_IWO 7 nm", cand["tch7"]["std_phys"]),
            ("+ C–V", cand["CV"]["std_phys"]), ("+ t_ox + t_IWO + C–V", cand["tox14+tch7+CV"]["std_phys"])]
    y = np.arange(len(rows))[::-1]
    h = 0.36
    c.barh(y + h / 2, [r[1]["log10_nd_eff"] for r in rows], height=h, color=BLUE, label="ND,eff")
    c.barh(y - h / 2, [r[1]["log10_n_ta"] for r in rows], height=h, color=ORANGE, label="N_TA (tail)")
    for yy, r in zip(y, rows):
        c.text(r[1]["log10_nd_eff"] * 1.08, yy + h / 2, f"{r[1]['log10_nd_eff']:.2f}", va="center", fontsize=6, color=INK2)
        c.text(r[1]["log10_n_ta"] * 1.08, yy - h / 2, f"{r[1]['log10_n_ta']:.2f}", va="center", fontsize=6, color=INK2)
    c.set_yticks(y)
    c.set_yticklabels([r[0] for r in rows], fontsize=7)
    c.set_xscale("log")
    c.set_xlim(0.08, 6)
    c.set_xlabel("1σ uncertainty (decades), σ_IV = 0.02 dec, σ_CV = 0.5 %")
    c.set_title("(c) Which extra measurement helps?", loc="left")
    c.legend(loc="lower right", fontsize=6.5, handlelength=1.2)
    style(c)

    # (d) DOS detectability
    E = np.array(e6["window_centres_eV"])
    for key, col in (("acceptor", BLUE), ("donor", ORANGE)):
        dos = np.array(e6["dos_model"][key])
        d.semilogy(E, np.where(dos > 0, dos, np.nan), "-", color=col, lw=2, label=f"{key}-like DOS (Step 5)")
    sig = np.array(e6["result"]["acceptor"]["sigma_profiled"])
    d.semilogy(E, sig, "--", color=INK, lw=1.5, label="1σ detection limit (either type)")
    d.text(0.58, 2.5e14, "donor ≡ acceptor + fixed charge\nN(1−f) = N − N·f", fontsize=6.3, color=INK2)
    fe = e6["fermi"]["0.05"]                         # Ec - EF (front interface) where ID crosses 10 pA
    fi, ff = np.log10(np.array(fe["id"])), np.array(fe["front"])
    k = int(np.argmax(fi >= -11))
    ef_floor = float(ff[k - 1] + (-11 - fi[k - 1]) / (fi[k] - fi[k - 1]) * (ff[k] - ff[k - 1]))
    d.axvspan(0, ef_floor, color=GRID, alpha=0.6, lw=0)
    d.text(0.02, 2e14, "EF swept while\nID ≥ 10 pA", fontsize=7, color=INK2)
    d.set_xlim(0, 1.2)
    d.set_ylim(1e14, 1e22)
    d.set_xlabel("Ec − E (eV)")
    d.set_ylabel("DOS / detection limit (cm⁻³ eV⁻¹)")
    d.set_title("(d) DOS is visible only where EF moves", loc="left")
    d.legend(loc="upper right", fontsize=6)
    style(d)

    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
