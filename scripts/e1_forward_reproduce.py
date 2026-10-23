"""E1 - forward reproduction of the calibration path of arXiv:2609.12184 in DEVSIM.

Runs the reported parameter sets (Step 0-5 for 2%-W, and the 4%-W
recalibration) at VD = 0.05 V and 1.0 V (L = 1 um), the final model at
L = 2 um (length-transfer test), and numerical-convergence checks
(DOS level count, mesh refinement). Extracts Vth (100 nA), SS (0.1-100 nA) and
Ion with the paper's definitions.

If digitized measurements are available, pass --meas data/meas.csv with columns
set,L_um,vd,vg,id (A) to compute metric errors and the objective J (own a_i).

Usage:  python scripts/e1_forward_reproduce.py [--quick] [--out results/e1]
        python scripts/e1_forward_reproduce.py --only w4 --out results/e1_w4     (re-run one set)
        python scripts/e1_forward_reproduce.py --report --merge results/e1_w4/summary.json
            (rebuild results/e1/summary.* with the separately run 4%-W set; no simulation)
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import replace

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from tftsim import PAPER_SETS, TFT, ConvergenceError, Geometry
from tftsim.metrics import extract, objective_J
from tftsim.params import default_dos_edges

VDS = (0.05, 1.0)


def sweep_device(params, geom=Geometry(), vds=VDS, vgs=None, **kw):
    vgs = vg_grid() if vgs is None else vgs
    t = TFT(params, geom, **kw).build()
    curves = {}
    for j, vd in enumerate(vds):
        seq = vgs if j % 2 == 0 else vgs[::-1]
        recs = t.transfer(vd, seq)
        recs.sort(key=lambda r: r["vg"])
        curves[vd] = (np.array([r["vg"] for r in recs]), np.array([r["id"] for r in recs]))
    return t, curves


def load_meas(path):
    data = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (r["set"], float(r["L_um"]))
            data.setdefault(key, {}).setdefault(float(r["vd"]), []).append((float(r["vg"]), float(r["id"])))
    out = {}
    for key, d in data.items():
        out[key] = {vd: (np.array([p[0] for p in sorted(v)]), np.array([p[1] for p in sorted(v)]))
                    for vd, v in d.items()}
    return out


def direction_checks(results):
    """Expected-direction checks against the paper's narrative, and numerical-convergence deltas.
    A check is skipped (not failed) when one of its runs is missing."""
    def m(name, vd=1.0):
        return results[name]["metrics"][str(vd)]
    check_defs = [
        ("Step1: lower S/D WF raises Ion", lambda: m("step1")["ion"] > m("step0")["ion"]),
        ("Step2: gate WF 4.4->4.0 eV shifts Vth negative (~0.4 V)", lambda: m("step2")["vth"] < m("step1")["vth"]),
        ("Step3: higher ND,eff shifts Vth negative", lambda: m("step3")["vth"] < m("step2")["vth"]),
        ("Step4: DOS changes the subthreshold slope", lambda: abs(m("step4")["ss_mVdec"] - m("step3")["ss_mVdec"]) > 1.0),
        ("Step5: mu0 10->24 raises Ion", lambda: m("step5")["ion"] > m("step4")["ion"]),
        ("4%-W: positive Vth shift vs 2%-W (paper Fig. 7a)", lambda: m("w4")["vth"] > m("step5")["vth"]),
        ("4%-W: lower current vs 2%-W (paper Fig. 7a)", lambda: m("w4")["ion"] < m("step5")["ion"]),
        ("L=2um: lower Ion than L=1um", lambda: m("step5_L2")["ion"] < m("step5")["ion"]),
    ]
    checks = []
    for txt, fn in check_defs:
        try:
            checks.append((txt, bool(fn())))
        except KeyError:
            pass
    convd = {}
    for name in ("step5_dosCoarse", "step5_dosFine", "step5_meshFine"):
        if name in results and "metrics" in results[name]:
            a, b = m(name), m("step5")
            convd[name] = {"dVth_mV": 1e3 * (a["vth"] - b["vth"]), "dSS_mVdec": a["ss_mVdec"] - b["ss_mVdec"],
                           "dlogIon_dec": a["log10_ion"] - b["log10_ion"]}
    return checks, convd


def write_summary(out, results, checks, convd, meas_report, wall_s, merged=None):
    dump_json({"results": results, "checks": checks, "convergence": convd, "meas": meas_report,
               "wall_s": wall_s, "merged": merged or {}}, f"{out}/summary.json")
    lines = ["# E1 forward reproduction", "", "| run | VD (V) | Vth (V) | SS (mV/dec) | Ion (A) |", "|---|---|---|---|---|"]
    for name, r in results.items():
        for vd, mm in r.get("metrics", {}).items():
            lines.append(f"| {name} | {vd} | {mm['vth']:+.3f} | {mm['ss_mVdec']:.1f} | {mm['ion']:.3e} |")
    if merged:
        lines += [""] + [f"- {name}: taken from {m['from']}" + (f" (the original run failed: {m['original_run_error']})"
                                                               if m.get("original_run_error") else "")
                         for name, m in merged.items()]
    lines += ["", "## Expected-direction checks", ""] + [f"- [{'x' if ok else ' '}] {txt}" for txt, ok in checks]
    lines += ["", "## Numerical convergence (vs default step5, VD = 1 V)", ""]
    lines += [f"- {k}: dVth {v['dVth_mV']:+.2f} mV, dSS {v['dSS_mVdec']:+.2f} mV/dec, dlogIon {v['dlogIon_dec']:+.4f} dec"
              for k, v in convd.items()]
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")


def report_only(out, merge):
    """Rebuild summary.md/json (checks, convergence) from stored results, no simulation.
    --merge adds runs from other summary.json files (e.g. a separate --only w4 run)."""
    s = json.load(open(f"{out}/summary.json"))
    results, merged = s["results"], dict(s.get("merged", {}))
    for path in merge:
        for name, r in json.load(open(path))["results"].items():
            if r.get("metrics"):
                note = {"from": path}
                if "error" in results.get(name, {}):              # keep why the original run has no result
                    note["original_run_error"] = results[name]["error"]
                elif isinstance(merged.get(name), dict):
                    note = dict(merged[name], **note)
                results[name] = r
                merged[name] = note
    checks, convd = direction_checks(results)
    write_summary(out, results, checks, convd, s.get("meas", {}), s.get("wall_s"), merged)
    print(f"{sum(ok for _, ok in checks)}/{len(checks)} direction checks pass")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e1")
    ap.add_argument("--quick", action="store_true", help="coarser VG grid, no convergence study")
    ap.add_argument("--meas", default=None)
    ap.add_argument("--only", default=None, help="comma-separated run names to (re)run")
    ap.add_argument("--report", action="store_true", help="rebuild the summary from stored results (no simulation)")
    ap.add_argument("--merge", nargs="*", default=[], help="with --report: summary.json files of separate runs to merge")
    args = ap.parse_args()
    out = outdir(args.out)
    if args.report:
        report_only(out, args.merge)
        return
    vgs = vg_grid(step=0.1) if args.quick else vg_grid()
    clock = Timer()

    runs = [(name, PAPER_SETS[name], Geometry(), {}) for name in
            ("step0", "step1", "step2", "step3", "step4", "step5", "w4")]
    runs.append(("step5_L2", PAPER_SETS["step5"], Geometry(L_um=2.0), {}))
    conv = []
    if not args.quick:
        conv = [
            ("step5_dosCoarse", PAPER_SETS["step5"], Geometry(), {"dos_edges": default_dos_edges(fine_step=0.05)}),
            ("step5_dosFine", PAPER_SETS["step5"], Geometry(), {"dos_edges": default_dos_edges(fine_step=0.0125)}),
            ("step5_meshFine", PAPER_SETS["step5"], Geometry(), {"mesh_scale": 0.5}),
        ]

    if args.only:
        keep = set(args.only.split(","))
        runs = [r for r in runs if r[0] in keep]
        conv = [r for r in conv if r[0] in keep]
    results, curves_all = {}, {}
    for name, p, g, kw in runs + conv:
        tm = Timer()
        vds = VDS if name in dict((r[0], 1) for r in runs) else (1.0,)
        try:
            t, curves = sweep_device(p, g, vds=vds, vgs=vgs, **kw)
        except ConvergenceError as exc:
            print(f"[{name}] FAILED: {exc}")
            results[name] = {"error": str(exc)}
            continue
        mets = {vd: extract(*curves[vd]) for vd in curves}
        results[name] = {"params": p.as_dict(), "L_um": g.L_um, "device": t.info(), "time_s": tm(),
                         "metrics": {str(vd): m for vd, m in mets.items()}}
        curves_all[name] = curves
        msg = "  ".join(f"VD={vd}: Vth={m['vth']:+.3f} V SS={m['ss_mVdec']:.1f} mV/dec Ion={m['ion']:.3e} A"
                        for vd, m in mets.items())
        print(f"[{name:16s}] {msg}  ({tm():.0f} s)")

    checks, convd = direction_checks(results)

    # ----- optional comparison with digitized measurements
    meas_report = {}
    if args.meas:
        meas = load_meas(args.meas)
        for (setname, L), mc in meas.items():
            key = setname if L == 1.0 else f"{setname}_L{int(L)}"
            if key in curves_all:
                meas_report[key] = objective_J(curves_all[key], mc)

    with open(f"{out}/curves.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run", "vd", "vg", "id_A"])
        for name, cv in curves_all.items():
            for vd, (vg, i) in cv.items():
                for a, b in zip(vg, i):
                    w.writerow([name, vd, a, b])
    write_summary(out, results, checks, convd, meas_report, clock())

    if HAVE_MPL and curves_all:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
        for j, vd in enumerate(VDS):
            for name in ("step0", "step1", "step2", "step3", "step4", "step5"):
                if name in curves_all:
                    vg, i = curves_all[name][vd]
                    ax[j].semilogy(vg, np.maximum(i, 1e-16), label=name)
            ax[j].axhline(100e-9, color="k", lw=0.6, ls=":")
            ax[j].axhline(10e-12, color="grey", lw=0.6, ls="--")
            ax[j].set_ylim(1e-14, 1e-2)
            ax[j].set_xlabel("VG (V)")
            ax[j].set_ylabel("ID (A)")
            ax[j].set_title(f"Calibration path, VD = {vd} V (L = 1 µm)")
            ax[j].legend(fontsize=8)
        for name, sty in (("step5", "-"), ("w4", "--"), ("step5_L2", ":")):
            if name in curves_all:
                vg, i = curves_all[name][1.0]
                ax[2].semilogy(vg, np.maximum(i, 1e-16), sty, label=name)
        ax[2].set_ylim(1e-14, 1e-2)
        ax[2].set_xlabel("VG (V)")
        ax[2].set_title("Transfer tests, VD = 1 V")
        ax[2].legend()
        fig.tight_layout()
        fig.savefig(f"{out}/e1_transfer.png", dpi=150)
    print(f"done in {clock():.0f} s -> {out}")


if __name__ == "__main__":
    main()
