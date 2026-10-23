"""E12 - profile likelihood of a parameter on the MEASURED 2 %-W curves.

Runs e11_calibrate_measured.py repeatedly with the parameter fixed, walking outward from the
best fit (continuation: each fit starts from the previous one). Then compares
  * Δχ² = χ²(value) − χ²(best), scaled with σ_eff of the best fit (model error included),
  * the fitted sets' errors at the paper's hold-out conditions against the MEASURED hold-out curves.

Usage:
  python scripts/e12_profile_measured.py --best results/e11_2W_v2/summary.json --param log10_nd_eff \
      --values 18.35,17.85,17.35,16.85 --out results/e12_nd
  python scripts/e12_profile_measured.py --best results/e11_2W_v3/summary.json --report --out results/e12_nd
(the stored table uses the longer free fit e11_2W_v3 as the best fit; see README)
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys

from _common import HAVE_MPL, dump_json, np, outdir


def run_side(best, param, values, out, extra):
    prev = best
    for v in values:
        tag = f"{param}_{v:g}".replace("-", "m")
        od = os.path.join(out, tag)
        cmd = [sys.executable, os.path.join(os.path.dirname(__file__), "e11_calibrate_measured.py"),
               "--set", "2W", "--start", prev, "--fix", f"{param}={v}", "--out", od] + extra
        print("running", " ".join(cmd[1:]), flush=True)
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(os.path.join(od, "summary.json")):
            print(f"[{param} = {v}] failed:\n{r.stdout[-800:]}\n{r.stderr[-800:]}", flush=True)
            break
        s = json.load(open(os.path.join(od, "summary.json")))
        print(f"[{param} = {v}] χ² {s['chi2']:.1f}  RMS {s['rms_dec']:.4f} dec", flush=True)
        prev = os.path.join(od, "summary.json")


def report(best_path, param, out):
    best = json.load(open(best_path))
    rows = [(best["theta_hat"][param], best, "best")]
    for f in sorted(glob.glob(os.path.join(out, f"{param}_*", "summary.json"))):
        s = json.load(open(f))
        rows.append((s["fixed"][param], s, os.path.basename(os.path.dirname(f))))
    rows.sort(key=lambda r: r[0])
    s2 = best["sigma_eff"] ** 2 / best["sigma_dec"] ** 2          # χ² scale factor for model error
    chi_min = min(s["chi2"] for _, s, _ in rows)                    # the free fit may stop short of the minimum
    best = dict(best, chi2_ref=chi_min)
    keys = [k for k in best.get("holdout", {}) if "error" not in best["holdout"][k]]
    lines = [f"# E12: profile likelihood of {param} on the measured 2 %-W calibration curves", "",
             f"- best fit: {param} = {best['theta_hat'][param]:.3f}, χ² = {best['chi2']:.1f} on {best['n_res']} residuals "
             f"(σ = {best['sigma_dec']} dec), σ_eff = {best['sigma_eff']:.4f} dec (model error included)",
             f"- Δχ²_eff = (χ² − χ²_min)·σ²/σ_eff² with χ²_min = {chi_min:.1f} (lowest over all fits); 95 % threshold 3.84", "",
             f"| {param} | Δχ²_eff | RMS (dec) | ΔVth VD=1 (mV) | " + " | ".join(f"hold-out {k}: E_log / ΔVth" for k in keys) + " |",
             "|---|---|---|---|" + "---|" * len(keys)]
    table = []
    for v, s, tag in rows:
        dchi = (s["chi2"] - chi_min) / s2
        ho = s.get("holdout", {})
        cells = [f"{ho[k]['E_log_dec']:.3f} / {ho[k]['dVth_mV']:+.0f}" if k in ho and "error" not in ho[k] else "–" for k in keys]
        lines.append(f"| {v:.3f} | {dchi:.2f} | {s['rms_dec']:.4f} | {s['cal']['1.0']['dVth_mV']:+.1f} | " + " | ".join(cells) + " |")
        table.append({"value": v, "dchi2_eff": dchi, "rms_dec": s["rms_dec"], "theta": s["theta_hat"],
                      "holdout": {k: {kk: vv for kk, vv in ho[k].items() if kk != "curve"} for k in keys if k in ho}})
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    dump_json({"param": param, "best": best_path, "rows": table}, os.path.join(out, "summary.json"))
    print("\n".join(lines))
    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
        ax[0].plot([t["value"] for t in table], [t["dchi2_eff"] for t in table], "o-")
        ax[0].axhline(3.84, color="r", ls="--", lw=0.8)
        ax[0].set_xlabel(param)
        ax[0].set_ylabel("Δχ² (σ_eff)")
        for k in keys:
            ax[1].plot([t["value"] for t in table], [t["holdout"].get(k, {}).get("E_log_dec", np.nan) for t in table], "o-", label=k)
        ax[1].set_xlabel(param)
        ax[1].set_ylabel("hold-out RMS error vs measured (dec)")
        ax[1].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "e12_profile.png"), dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--best", required=True)
    ap.add_argument("--param", default="log10_nd_eff")
    ap.add_argument("--values", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--max_nfev", default="20")
    args = ap.parse_args()
    out = outdir(args.out)
    if not args.report:
        vals = [float(v) for v in args.values.split(",") if v]
        run_side(args.best, args.param, vals, out, ["--max_nfev", args.max_nfev])
    report(args.best, args.param, out)


if __name__ == "__main__":
    main()
