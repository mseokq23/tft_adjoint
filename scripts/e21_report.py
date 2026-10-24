"""E21 report - robustness of the measured-data conclusions (no simulation; reads stored fits).

1. Channel permittivity: the two ends of the measured ND,eff profile (10^16.85 and 10^19.05) re-fitted with
   ε_ch = 15 instead of 10 (results/e21_eps15), as E16 did for Nc. Equivalent if the two ends stay within
   Δχ²_eff < 3.84 of each other (σ_eff of the better end).
2. Multi-start free fits (results/e21_multistart): free calibrations to the same measured 2 %-W pair
   (L = 1 µm, VD = 0.05/1 V, 40 evaluations each) from different starting points, next to the two earlier free
   fits (v2: path from Step 5 via the first-stage fit; v3: from the E12 10^19.05 point). Δχ²_eff is taken
   against the lowest χ² of all fits with the default material (E11/E12/E15/E17 and these), σ_eff of v3.

Usage:  python scripts/e21_report.py [--out results/e21]
"""
from __future__ import annotations

import argparse
import glob
import json
import os

from _common import dump_json, np, outdir
from e2_twin_calibration import BOUNDS
from tftsim.metrics import extract

SHOW = ["wf_sd", "wf_gate", "log10_nd_eff", "log10_mu0", "log10_n_gd", "e_gd", "s_gd", "log10_n_ta", "w_ta"]
START = {"step5": "paper Step 5 (2%-W)", "w4": "paper 4%-W values", "mid": "midpoints of the fit bounds",
         "low_nd": "E12 set at ND,eff = 10^16.85 (ND,eff free)"}
# scipy least_squares status codes
STATUS = {0: "evaluation limit", 1: "gtol", 2: "ftol", 3: "xtol", 4: "ftol+xtol"}
TESTS = [("0.05", "cal VD 0.05"), ("1.0", "cal VD 1"), ("L1_VD0.1", "VD 0.1"), ("L1_VD0.2", "VD 0.2"), ("L1_VD0.5", "VD 0.5"),
         ("L2_VD0.05", "L 2 µm VD 0.05"), ("L2_VD1", "L 2 µm VD 1")]


def load(p):
    return json.load(open(p))


def vth_full_mV(s):
    """Vth (100 nA) of every calibration and transfer-test curve, from the full simulated curve of the set
    (the digitised transfer-test curves have gaps around 100 nA, so set-to-set comparisons use the model curves)."""
    vg = np.array(s["vg"])
    out = {}
    for key, label in TESTS:
        if key in s["curves"]:
            cur = s["curves"][key]
        else:
            h = s["holdout"].get(key, {})
            cur = h.get("curve")
        if cur is not None:
            out[label] = 1e3 * extract(vg, np.array(cur))["vth"]
    return out


def at_bounds(theta, tol=1e-3, near=1e-2):
    """Descriptors at a fit bound (within tol of the bound range) or near one (within `near`)."""
    out = []
    for n in SHOW:
        lo, hi = BOUNDS[n]
        v, w = theta[n], hi - lo
        for side, d in (("lower", v - lo), ("upper", hi - v)):
            if d <= tol * w:
                out.append(f"{n} ({side})")
            elif d <= near * w:
                out.append(f"{n} (near {side}, {d:.3g} from it)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e21")
    args = ap.parse_args()
    out = outdir(args.out)
    v3 = load("results/e11_2W_v3/summary.json")
    scale = (v3["sigma_eff"] / v3["sigma_dec"]) ** 2
    default_fits = (["results/e11_2W_v2/summary.json", "results/e11_2W_v3/summary.json"]
                    + sorted(glob.glob("results/e12_nd/log10_nd_eff_*/summary.json"))
                    + sorted(glob.glob("results/e15_wta/*/summary.json")) + sorted(glob.glob("results/e17_cross/*/summary.json")))
    found = {os.path.basename(os.path.dirname(f)): f for f in glob.glob("results/e21_multistart/*/summary.json")}
    order = list(START) + sorted(k for k in found if k not in START)
    ms = {k: found[k] for k in order if k in found}
    chi_min = min(load(f)["chi2"] for f in default_fits + list(ms.values()))
    res = {"chi2_min": chi_min, "sigma_eff_v3": v3["sigma_eff"]}
    lines = ["# E21: robustness of the measured-data conclusions", ""]

    # 1. permittivity
    eps = sorted(glob.glob("results/e21_eps15/*/summary.json"))
    if eps:
        ss = [load(f) for f in eps]
        best = min(ss, key=lambda s: s["chi2"])
        sc = (best["sigma_eff"] / best["sigma_dec"]) ** 2
        rows = []
        lines += ["## 1. Channel permittivity ε_ch = 15 (default 10)", "",
                  "| fixed log10 ND,eff | χ² (ε_ch = 15) | χ² (ε_ch = 10, E12) | Δχ²_eff vs the better end | Δχ²_eff vs default-material χ²_min | "
                  "RMS (dec) | nfev | status | parameters at a bound |", "|---|---|---|---|---|---|---|---|---|"]
        for f, s in zip(eps, ss):
            nd = s["fixed"]["log10_nd_eff"]
            ref_f = f"results/e12_nd/log10_nd_eff_{nd:.2f}/summary.json"
            ref = load(ref_f) if os.path.exists(ref_f) else None
            r = {"file": f, "log10_nd_eff": nd, "chi2": s["chi2"], "chi2_eps10": ref["chi2"] if ref else None,
                 "dchi2_eff_vs_better_end": (s["chi2"] - best["chi2"]) / sc,
                 "dchi2_eff_vs_default_min": (s["chi2"] - chi_min) / scale, "rms_dec": s["rms_dec"], "nfev": s["nfev"],
                 "status": s["status"], "material": s["material"], "at_bounds": at_bounds(s["theta_hat"]), "theta": s["theta_hat"],
                 "theta_eps10": {n: ref["theta_hat"][n] for n in SHOW} if ref else None}
            rows.append(r)
            c10 = f"{r['chi2_eps10']:.1f}" if ref else "–"
            lines.append(f"| {nd:.2f} | {r['chi2']:.1f} | {c10} | {r['dchi2_eff_vs_better_end']:.2f} | "
                         f"{r['dchi2_eff_vs_default_min']:.2f} | {r['rms_dec']:.4f} | {r['nfev']} | {STATUS.get(r['status'], r['status'])} | "
                         f"{', '.join(r['at_bounds']) or '–'} |")
        res["eps15"] = rows
        lines += ["", "σ_eff of the better end is used for its column, σ_eff of v3 for the default-material column. "
                  "Both ends are equivalent (95 %) if the first column stays below 3.84. Each ε_ch = 15 fit starts from the "
                  "ε_ch = 10 point with the same ND,eff (E12).", "",
                  "Fitted descriptors with ε_ch = 10 (E12) and 15 (these fits):", "",
                  "| descriptor | " + " | ".join(f"{r['log10_nd_eff']:.2f}: ε_ch 10 | {r['log10_nd_eff']:.2f}: ε_ch 15" for r in rows) + " |",
                  "|---|" + "---|---|" * len(rows)]
        for n in SHOW:
            if n == "log10_nd_eff":
                continue
            lines.append(f"| {n} | " + " | ".join((f"{r['theta_eps10'][n]:.3f}" if r["theta_eps10"] else "–") + f" | {r['theta'][n]:.3f}"
                                                for r in rows) + " |")
        lines += [""]

    # 2. multi-start
    fits = {"v2 (earlier, Step 5 path)": "results/e11_2W_v2/summary.json",
            "v3 (earlier, from the 10^19.05 point)": "results/e11_2W_v3/summary.json"}
    fits.update({f"{k}: {START.get(k, k)}": f for k, f in ms.items()})
    if ms:
        lines += ["## 2. Free fits from different starting points (2 %-W pair, Qf at the Step-5 value)", "",
                  f"Δχ²_eff against χ²_min = {chi_min:.1f} (lowest of all default-material fits), σ_eff = {v3['sigma_eff']:.4f} dec (v3).", "",
                  "| fit | χ² | Δχ²_eff | RMS (dec) | nfev | status | " + " | ".join(SHOW) + " | parameters at a bound |",
                  "|---|---|---|---|---|---|" + "---|" * len(SHOW) + "---|"]
        rows = []
        for name, f in fits.items():
            s = load(f)
            th = s["theta_hat"]
            r = {"fit": name, "file": f, "chi2": s["chi2"], "dchi2_eff": (s["chi2"] - chi_min) / scale, "rms_dec": s["rms_dec"],
                 "nfev": s.get("nfev"), "status": s.get("status"), "theta": {n: th[n] for n in SHOW}, "at_bounds": at_bounds(th),
                 "vth_full_mV": vth_full_mV(s)}
            rows.append(r)
            lines.append(f"| {name} | {r['chi2']:.1f} | {r['dchi2_eff']:.2f} | {r['rms_dec']:.4f} | {r['nfev']} | "
                         f"{STATUS.get(r['status'], r['status'])} | " + " | ".join(f"{th[n]:.3f}" for n in SHOW)
                         + f" | {', '.join(r['at_bounds']) or '–'} |")
        res["multistart"] = rows
        sd_names = ["log10_nd_eff", "log10_n_ta", "log10_n_gd", "wf_gate", "w_ta", "s_gd"]
        lines += ["", "Calibration-only Fisher 1σ stored with each fit (σ_eff of that fit, linearised at its own end point):", "",
                  "| fit | " + " | ".join(sd_names) + " |", "|---|" + "---|" * len(sd_names)]
        for name, f in fits.items():
            sd = load(f).get("std_eff", {})
            lines.append(f"| {name} | " + " | ".join(f"{sd[n]:.3g}" if n in sd else "–" for n in sd_names) + " |")
        res["calibration_only_std_eff"] = {name: {n: load(f).get("std_eff", {}).get(n) for n in sd_names} for name, f in fits.items()}
        labels = [lab for _, lab in TESTS]
        lines += ["", "Predicted Vth (mV, 100 nA) of the calibration and transfer-test curves, from each fit's full simulated curve:", "",
                  "| fit | " + " | ".join(labels) + " |", "|---|" + "---|" * len(labels)]
        for r in rows:
            v = r["vth_full_mV"]
            lines.append(f"| {r['fit']} | " + " | ".join(f"{v[lab]:.1f}" if lab in v else "–" for lab in labels) + " |")
        eq = [r for r in rows if r["dchi2_eff"] < 3.84]
        if len(eq) > 1:
            spread = {n: max(r["theta"][n] for r in eq) - min(r["theta"][n] for r in eq) for n in SHOW}
            vspread = {lab: max(r["vth_full_mV"][lab] for r in eq) - min(r["vth_full_mV"][lab] for r in eq)
                       for lab in labels if all(lab in r["vth_full_mV"] for r in eq)}
            res["multistart_equivalent"] = {"fits": [r["fit"] for r in eq], "descriptor_spread": spread, "vth_spread_mV": vspread}
            lines += ["", f"{len(eq)} of {len(rows)} fits are within Δχ²_eff < 3.84 of the best: " + ", ".join(r["fit"].split(":")[0].split(" (")[0] for r in eq) + ".",
                      "", "- spread of their descriptors (max − min): " + ", ".join(f"{n} {v:.3g}" for n, v in spread.items()) + ".",
                      "- spread of their predicted Vth (mV): " + ", ".join(f"{lab} {v:.1f}" for lab, v in vspread.items()) + "."]
    dump_json(res, os.path.join(out, "summary.json"))
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
