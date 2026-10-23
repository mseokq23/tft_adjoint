"""E13 - are the parameter sets that fit the MEASURED data equally well physically different?

E12 re-fitted our model to the measured 2 %-W calibration curves with ND,eff fixed between
10^16.85 and 10^19.05 cm^-3: all sets fit equally well (Δχ²_eff ≤ 1.26 against the best free fit
e11_2W_v3; ≤ 0.77 against the reference used before e11_2W_v3 existed) and give the same errors
on the paper's transfer tests (drain bias VD = 0.1/0.2/0.5 V, channel length L = 2 µm). This script
asks whether those sets nevertheless predict different results for measurements outside the paper's
split - the measurements the OED analysis (E4/E7) recommends:
  * I-V of a device with a thicker channel (t_IWO = 7 nm), VD = 0.05 / 1.0 V
  * I-V of a device with a thicker gate oxide (t_ox = 14 nm), VD = 0.05 / 1.0 V
  * quasi-static gate C-V of the calibration device (VD = 0)
For every condition: spread of the predictions across the equivalent sets (Vth, RMS log-current
deviation from the reference set = lowest χ² on the measured data), next to the spread the same sets
show on the paper's own transfer tests (from E12). Thresholds: digitisation noise 0.02 dec and the
model error σ_eff of the measured fit (0.12 dec); for C-V 0.5 % and 2 % of C_ox·A.

The thick-channel device turns on near the start of the default sweep (Vth = -0.77 ... -1.15 V at
VD = 1 V across the sets), so its curves are re-run from VG = -3 V:
        python scripts/e13_measured_discrimination.py --recompute "t_IWO 7 nm" --vgmin_iv -3.0

Usage:  python scripts/e13_measured_discrimination.py [--profile results/e12_nd]
            [--best results/e11_2W_v3/summary.json,results/e11_2W_v2/summary.json] [--out results/e13] [--report]
"""
from __future__ import annotations

import argparse
import glob
import json
import os

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from e9_discrimination_test import cv, iv, params_of
from tftsim import Geometry
from tftsim.metrics import I_FLOOR, extract

CONDS = {"t_IWO 7 nm": Geometry(t_iwo_nm=7.0), "t_ox 14 nm": Geometry(t_gox_nm=14.0)}
TRANSFER = ["L1_VD0.1", "L1_VD0.2", "L1_VD0.5", "L2_VD0.05", "L2_VD1"]


def load_sets(profile, best):
    """best: one or more comma-separated free-fit summary.json files (all are kept as sets)."""
    sets = {}
    for f in best.split(","):
        b = json.load(open(f))
        sets[f"{b['theta_hat']['log10_nd_eff']:.3f}"] = b
    for f in sorted(glob.glob(os.path.join(profile, "log10_nd_eff_*", "summary.json"))):
        s = json.load(open(f))
        sets[f"{s['fixed']['log10_nd_eff']:.3f}"] = s
    return dict(sorted(sets.items(), key=lambda kv: float(kv[0])))


def simulate(sets, out, clock, recompute=(), vgmin_iv=-1.0):
    path = os.path.join(out, "predictions.json")
    pred = json.load(open(path)) if os.path.exists(path) else {}
    vgs_iv, vgs_cv = vg_grid(-1.0, 2.5, 0.1), vg_grid(-1.0, 2.5, 0.05)
    for key, s in sets.items():
        if key in pred:
            todo = [c for c in recompute if pred[key].get("vgmin", {}).get(c) != vgmin_iv]
            if todo:                                # re-run selected conditions on a wider VG grid
                p = params_of(s["theta_hat"])
                grid = vg_grid(vgmin_iv, 2.5, 0.1)
                for c in todo:
                    cur = iv(p, CONDS[c], grid)
                    pred[key]["iv"][c] = {f"{vd:g}": cur[vd].tolist() for vd in cur}
                    pred[key].setdefault("vg", {})[c] = grid.tolist()
                    pred[key].setdefault("vgmin", {})[c] = vgmin_iv
                dump_json(pred, path)
                print(f"[ND,eff = 10^{key}] re-ran {', '.join(todo)} from VG = {vgmin_iv} V ({clock():.0f} s)", flush=True)
            continue
        p = params_of(s["theta_hat"])
        row = {"iv": {}, "vg_iv": vgs_iv.tolist(), "vg_cv": (0.5 * (vgs_cv[1:] + vgs_cv[:-1])).tolist()}
        for c, g in CONDS.items():
            cur = iv(p, g, vgs_iv)
            row["iv"][c] = {f"{vd:g}": cur[vd].tolist() for vd in cur}
        c_v, coxa = cv(p, vgs_cv)
        row["cv"], row["coxa"] = c_v.tolist(), float(coxa)
        pred[key] = row
        dump_json(pred, path)                       # partial results survive an interruption
        print(f"[ND,eff = 10^{key}] done ({clock():.0f} s)", flush=True)
    return pred


def compare(sets, pred, best):
    sig, sig_eff = best["sigma_dec"], best["sigma_eff"]
    chi_min = min(s["chi2"] for s in sets.values())
    scale = (sig_eff / sig) ** 2
    ref = min(sets, key=lambda k: sets[k]["chi2"])
    rows = {}
    for key, s in sets.items():
        r = {"chi2": s["chi2"], "dchi2_eff": (s["chi2"] - chi_min) / scale, "rms_dec": s["rms_dec"],
             "transfer": {k: {m: s["holdout"][k][m] for m in ("dVth_mV", "dlogIon_dec", "E_log_dec")}
                          for k in TRANSFER if k in s.get("holdout", {}) and "error" not in s["holdout"][k]}}
        for c in CONDS:
            r[c] = {}
            vg = np.array(pred[ref].get("vg", {}).get(c, pred[ref]["vg_iv"]))
            for vd in ("0.05", "1"):
                a, b = np.array(pred[ref]["iv"][c][vd]), np.array(pred[key]["iv"][c][vd])
                m = (a >= I_FLOOR) & (b >= I_FLOOR)
                dev = np.log10(b[m]) - np.log10(a[m])
                ea, eb = extract(vg, a), extract(vg, b)
                r[c][vd] = {"vth_V": eb["vth"], "dVth_mV": 1e3 * (eb["vth"] - ea["vth"]),
                            "dlogIon_dec": eb["log10_ion"] - ea["log10_ion"],
                            "rms_dec": float(np.sqrt(np.mean(dev ** 2))), "max_dec": float(np.max(np.abs(dev))),
                            "chi2_sigma": float(np.sum((dev / sig) ** 2)), "chi2_sigma_eff": float(np.sum((dev / sig_eff) ** 2)),
                            "n": int(m.sum())}
        ca, cb = np.array(pred[ref]["cv"]), np.array(pred[key]["cv"])
        dc = (cb - ca) / pred[ref]["coxa"]
        r["C-V"] = {"max_dev_pct": float(100 * np.max(np.abs(dc))), "chi2_0.5pct": float(np.sum((dc / 0.005) ** 2)),
                    "chi2_2pct": float(np.sum((dc / 0.02) ** 2))}
        rows[key] = r
    return ref, rows


def spread(vals):
    vals = [v for v in vals if v is not None and np.isfinite(v)]
    return float(max(vals) - min(vals)) if vals else float("nan")


def report(sets, pred, best, out):
    ref, rows = compare(sets, pred, best)
    keys = list(rows)
    lines = ["# E13: do the sets that fit the measured data equally well differ physically?", "",
             f"Sets: measured 2 %-W profile of ND,eff (E12) and the free best fit; reference = lowest χ² (log10 ND,eff = {ref}). "
             f"σ = {best['sigma_dec']} dec (digitisation), σ_eff = {best['sigma_eff']:.3f} dec (model error of the measured fit).", "",
             "## Per set", "",
             "| log10 ND,eff | Δχ²_eff (calibration) | RMS calib (dec) | t_IWO 7 nm VD=1: ΔVth (mV) / RMS (dec) | t_ox 14 nm VD=1: ΔVth (mV) / RMS (dec) | C–V max ΔC (% C_ox) |",
             "|---|---|---|---|---|---|"]
    for k in keys:
        r = rows[k]
        lines.append(f"| {k} | {r['dchi2_eff']:.2f} | {r['rms_dec']:.4f} | "
                     f"{r['t_IWO 7 nm']['1']['dVth_mV']:+.1f} / {r['t_IWO 7 nm']['1']['rms_dec']:.3f} | "
                     f"{r['t_ox 14 nm']['1']['dVth_mV']:+.1f} / {r['t_ox 14 nm']['1']['rms_dec']:.3f} | {r['C-V']['max_dev_pct']:.2f} |")
    lines += ["", "## Spread across the equivalent sets (max − min)", "",
              "| test | in the paper's split? | Vth spread (mV) | Ion spread (dec) | largest RMS deviation between sets (dec) |",
              "|---|---|---|---|---|"]
    summ = {"reference": ref, "tests": {}}
    for k in TRANSFER:
        dv = [rows[s]["transfer"].get(k, {}).get("dVth_mV") for s in keys]
        di = [rows[s]["transfer"].get(k, {}).get("dlogIon_dec") for s in keys]
        summ["tests"][k] = {"paper": True, "vth_spread_mV": spread(dv), "ion_spread_dec": spread(di)}
        lines.append(f"| {k.replace('_', ' ')} | yes (transfer test) | {spread(dv):.1f} | {spread(di):.3f} | – |")
    for c in CONDS:
        for vd in ("0.05", "1"):
            dv = [rows[s][c][vd]["dVth_mV"] for s in keys]
            di = [rows[s][c][vd]["dlogIon_dec"] for s in keys]
            rm = max(rows[s][c][vd]["rms_dec"] for s in keys)
            tag = f"{c} VD={vd}"
            summ["tests"][tag] = {"paper": False, "vth_spread_mV": spread(dv), "ion_spread_dec": spread(di), "max_rms_dec": rm}
            lines.append(f"| {tag} | no | {spread(dv):.1f} | {spread(di):.3f} | {rm:.3f} |")
    cvm = max(rows[s]["C-V"]["max_dev_pct"] for s in keys)
    summ["tests"]["C-V"] = {"paper": False, "max_dev_pct": cvm}
    lines.append(f"| quasi-static C–V | no | – | – | max ΔC {cvm:.2f} % of C_ox·A |")
    lines += ["", "Transfer-test spreads are differences between the sets' errors against the MEASURED hold-out curves (E12); "
              "new-test spreads are differences between the sets' predictions (no data exist)."]
    dump_json({"sets": list(sets), "rows": rows, **summ}, os.path.join(out, "summary.json"))
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
        cmap = plt.get_cmap("viridis")
        for j, k in enumerate(keys):
            col = cmap(j / max(1, len(keys) - 1))
            for a, c in zip(ax[:2], CONDS):
                vg = np.array(pred[k].get("vg", {}).get(c, pred[k]["vg_iv"]))
                a.semilogy(vg, pred[k]["iv"][c]["1"], color=col, lw=1.3, label=f"10^{k}")
            ax[2].plot(pred[k]["vg_cv"], np.array(pred[k]["cv"]) / pred[k]["coxa"], color=col, lw=1.3)
        for a, c in zip(ax[:2], CONDS):
            a.set_title(f"{c}, VD = 1 V (predicted)")
            a.set_xlabel("VG (V)")
            a.set_ylim(1e-12, 1e-2)
        ax[0].set_ylabel("ID (A)")
        ax[0].legend(fontsize=7, title="ND,eff (cm⁻³)")
        ax[2].set_title("quasi-static C–V (predicted)")
        ax[2].set_xlabel("VG (V)")
        ax[2].set_ylabel("C / (C_ox A)")
        fig.tight_layout()
        fig.savefig(os.path.join(out, "e13_predictions.png"), dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="results/e12_nd")
    ap.add_argument("--best", default="results/e11_2W_v3/summary.json,results/e11_2W_v2/summary.json",
                    help="free-fit summary.json file(s), comma-separated; the lowest χ² sets σ_eff")
    ap.add_argument("--out", default="results/e13")
    ap.add_argument("--report", action="store_true", help="no new simulations")
    ap.add_argument("--recompute", default="", help="comma-separated conditions to re-run on the VG grid from --vgmin_iv")
    ap.add_argument("--vgmin_iv", type=float, default=-1.0)
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    sets = load_sets(args.profile, args.best)
    best = min((json.load(open(f)) for f in args.best.split(",")), key=lambda s: s["chi2"])
    print("sets:", ", ".join(sets), flush=True)
    recompute = [c.strip() for c in args.recompute.split(",") if c.strip()]
    pred = (json.load(open(os.path.join(out, "predictions.json"))) if args.report
            else simulate(sets, out, clock, recompute, args.vgmin_iv))
    report(sets, pred, best, out)
    print(f"done in {clock():.0f} s")


if __name__ == "__main__":
    main()
