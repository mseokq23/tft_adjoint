"""E10 - the paper's reported parameter sets in our open-source model, against the digitised measurements.

Measured curves: data/meas.csv (scripts/digitize_paper.py; panel-to-panel σ ≈ 0.02 dec).
For every comparison the simulation is evaluated on the measured VG grid (−1 … 2 V, 0.05 V)
and the paper's metrics are extracted from both: Vth at 100 nA, SS between 0.1 and 100 nA,
Ion = max over the measured sweep, RMS log-current error above 10 pA.

Runs
  * Step 0 … 5 (2 %-W, L = 1 µm, VD = 0.05 / 1 V): objective J with the paper's weights and
    our own normalisation (the paper's a_i are unpublished), per-metric errors
  * Step 5 at the paper's hold-out conditions (VD = 0.1/0.2/0.5 V; L = 2 µm, VD = 0.05/1 V)
  * the 4 %-W recalibration against the 4 %-W data (VG ≥ −0.7 V, as in the paper)

Usage:  python scripts/e10_measured_comparison.py [--meas data/meas.csv] [--out results/e10]
"""
from __future__ import annotations

import argparse
import collections
import csv

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from e1_forward_reproduce import sweep_device
from tftsim import PAPER_SETS, PAPER_J, ConvergenceError, Geometry
from tftsim.metrics import I_FLOOR, extract, objective_J

PAPER_STEP5_VD1 = {"dVth_mV": 31.5, "dlogIon_dec": 0.022}      # Fig. 5(b) of the paper


def load_meas(path):
    d = collections.defaultdict(list)
    for r in csv.DictReader(open(path)):
        d[(r["set"], float(r["L_um"]), float(r["vd"]))].append((float(r["vg"]), float(r["id"])))
    return {k: (np.array([p[0] for p in sorted(v)]), np.array([p[1] for p in sorted(v)])) for k, v in d.items()}


def compare(sim_vg, sim_i, m_vg, m_i):
    s = np.interp(m_vg, sim_vg, sim_i)
    ok = m_i >= I_FLOOR
    e_log = float(np.sqrt(np.mean((np.log10(np.maximum(s[ok], 1e-30)) - np.log10(m_i[ok])) ** 2)))
    ms, mm = extract(m_vg, s), extract(m_vg, m_i)
    return {"E_log_dec": e_log, "dVth_mV": 1e3 * (ms["vth"] - mm["vth"]), "dSS_mVdec": ms["ss_mVdec"] - mm["ss_mVdec"],
            "dlogIon_dec": ms["log10_ion"] - mm["log10_ion"], "sim": ms, "meas": mm}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meas", default="data/meas.csv")
    ap.add_argument("--out", default="results/e10")
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    meas = load_meas(args.meas)
    vgs = vg_grid(-1.0, 2.0, 0.05)
    res, curves = {}, {}

    # ---- measured metrics (facts about the digitised data)
    mfacts = {f"{k[0]} L={k[1]:g} VD={k[2]:g}": extract(*meas[k]) for k in sorted(meas)}

    # ---- calibration path
    for name in ("step0", "step1", "step2", "step3", "step4", "step5"):
        try:
            _, cv = sweep_device(PAPER_SETS[name], Geometry(), vds=(0.05, 1.0), vgs=vgs)
        except ConvergenceError as exc:
            res[name] = {"error": str(exc)}
            continue
        curves[name] = cv
        mc = {vd: meas[("2W", 1.0, vd)] for vd in (0.05, 1.0)}
        J = objective_J(cv, mc)
        res[name] = {"J_ours": J["J"], "J_paper": PAPER_J[name], "terms": J["terms"],
                     "cmp": {str(vd): compare(*cv[vd], *mc[vd]) for vd in (0.05, 1.0)}}
        c1 = res[name]["cmp"]["1.0"]
        print(f"[{name}] J(ours) {J['J']:.3f} (paper {PAPER_J[name]})  VD=1: dVth {c1['dVth_mV']:+.0f} mV, "
              f"dSS {c1['dSS_mVdec']:+.0f} mV/dec, dlogIon {c1['dlogIon_dec']:+.3f}, E_log {c1['E_log_dec']:.3f} ({clock():.0f} s)")

    # ---- hold-out conditions with the final 2 %-W model
    ho = {}
    for geom, vds, tag in ((Geometry(), (0.1, 0.2, 0.5), "L1"), (Geometry(L_um=2.0), (0.05, 1.0), "L2")):
        _, cv = sweep_device(PAPER_SETS["step5"], geom, vds=vds, vgs=vgs)
        for vd in vds:
            key = ("2W", geom.L_um, vd)
            if key in meas:
                ho[f"{tag}_VD{vd:g}"] = compare(*cv[vd], *meas[key])
                curves[f"step5_{tag}_VD{vd:g}"] = {vd: cv[vd]}
    print(f"hold-out done ({clock():.0f} s)")

    # ---- 4 %-W recalibration
    w4 = {}
    vgs4 = vg_grid(-0.7, 2.0, 0.05)
    try:
        _, cv = sweep_device(PAPER_SETS["w4"], Geometry(), vds=(0.05, 1.0), vgs=vgs4)
        curves["w4"] = cv
        for vd in (0.05, 1.0):
            k = ("4W", 1.0, vd)
            if k in meas:
                m_vg, m_i = meas[k]
                sel = m_vg >= -0.7 - 1e-9
                w4[str(vd)] = compare(*cv[vd], m_vg[sel], m_i[sel])
    except ConvergenceError as exc:
        w4 = {"error": str(exc)}
    print(f"4%-W done ({clock():.0f} s)")

    dump_json({"measured_metrics": mfacts, "path": res, "holdout_step5": ho, "w4": w4,
               "paper_step5_vd1": PAPER_STEP5_VD1, "curves": {k: {str(vd): v for vd, v in cv.items()} for k, cv in curves.items()},
               "wall_s": clock()}, f"{out}/summary.json")

    lines = ["# E10: paper parameter sets in DEVSIM vs digitised measurements", "",
             "Our J uses the paper's weights with our own normalisation (paper a_i unpublished): compare trends, not absolute J.", "",
             "## Measured metrics (digitised)", "", "| curve | Vth (V) | SS (mV/dec) | Ion (A) |", "|---|---|---|---|"]
    for k, m in mfacts.items():
        lines.append(f"| {k} | {m['vth']:+.3f} | {m['ss_mVdec']:.1f} | {m['ion']:.3e} |")
    lines += ["", "## Calibration path (2 %-W, L = 1 µm)", "",
              "| step | J ours | J paper | VD | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |", "|---|---|---|---|---|---|---|---|"]
    for name, r in res.items():
        if "error" in r:
            lines.append(f"| {name} | failed | {PAPER_J[name]} | | | | | |")
            continue
        for vd in ("0.05", "1.0"):
            c = r["cmp"][vd]
            lines.append(f"| {name} | {r['J_ours']:.3f} | {r['J_paper']} | {vd} | {c['E_log_dec']:.3f} | {c['dVth_mV']:+.0f} | "
                         f"{c['dSS_mVdec']:+.1f} | {c['dlogIon_dec']:+.3f} |")
    lines += ["", f"Paper (Sentaurus) Step 5 at VD = 1 V: |ΔVth| = {PAPER_STEP5_VD1['dVth_mV']} mV, |ΔIon| = {PAPER_STEP5_VD1['dlogIon_dec']} dec", "",
              "## Step 5 at the paper's hold-out conditions", "",
              "| condition | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |", "|---|---|---|---|---|"]
    for k, c in ho.items():
        lines.append(f"| {k} | {c['E_log_dec']:.3f} | {c['dVth_mV']:+.0f} | {c['dSS_mVdec']:+.1f} | {c['dlogIon_dec']:+.3f} |")
    lines += ["", "Paper: drain-bias test avg/max |ΔVth| 24.9/36.1 mV, |ΔIon| 0.013/0.022 dec; L = 2 µm 27.4/46.2 mV, 0.050/0.062 dec", "",
              "## 4 %-W recalibration vs 4 %-W data (VG ≥ −0.7 V)", ""]
    if "error" in w4:
        lines.append(f"failed: {w4['error']}")
    else:
        for vd, c in w4.items():
            lines.append(f"- VD = {vd} V: E_log {c['E_log_dec']:.3f} dec, ΔVth {c['dVth_mV']:+.0f} mV, ΔSS {c['dSS_mVdec']:+.1f} mV/dec, "
                         f"Δlog Ion {c['dlogIon_dec']:+.3f} dec")
    lines += ["", "Paper 4 %-W: avg/max |ΔIon| 0.031/0.057 dec, |ΔVth| 28.2/47.3 mV"]
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
        for j, vd in enumerate((1.0, 0.05)):
            m_vg, m_i = meas[("2W", 1.0, vd)]
            ax[j].semilogy(m_vg, m_i, "s", mfc="none", color="k", ms=4, label="measured (digitised)")
            for name in ("step0", "step1", "step2", "step3", "step4", "step5"):
                if name in curves:
                    ax[j].semilogy(*curves[name][vd], lw=1.4, label=name)
            ax[j].set_ylim(1e-13, 1e-2)
            ax[j].set_xlabel("VG (V)")
            ax[j].set_ylabel("ID (A)")
            ax[j].set_title(f"2 %-W, L = 1 µm, VD = {vd} V: paper parameters in DEVSIM")
            ax[j].legend(fontsize=7)
        for k, c in ho.items():
            pass
        for key, lab in ((("2W", 2.0, 1.0), "L2 VD1"), (("2W", 2.0, 0.05), "L2 VD0.05"), (("4W", 1.0, 1.0), "4%-W VD1")):
            if key in meas:
                ax[2].semilogy(*meas[key], "o", mfc="none", ms=4, label=f"meas {lab}")
        for nm, lab in (("step5_L2_VD1", "Step 5 L2 VD1"), ("step5_L2_VD0.05", "Step 5 L2 VD0.05")):
            if nm in curves:
                vd = list(curves[nm])[0]
                ax[2].semilogy(*curves[nm][vd], lw=1.4, label=lab)
        if "w4" in curves:
            ax[2].semilogy(*curves["w4"][1.0], lw=1.4, label="4%-W set VD1")
        ax[2].set_ylim(1e-13, 1e-2)
        ax[2].set_xlabel("VG (V)")
        ax[2].legend(fontsize=7)
        ax[2].set_title("Hold-out and 4 %-W")
        fig.tight_layout()
        fig.savefig(f"{out}/e10_measured.png", dpi=150)
    print(f"done in {clock():.0f} s -> {out}")


if __name__ == "__main__":
    main()
