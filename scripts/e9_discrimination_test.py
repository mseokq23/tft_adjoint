"""E9 - nonlinear check of the OED recommendation: do the extra measurements separate the
parameter sets that the paper's data (calibration + hold-out) cannot?

Takes the re-fitted sets of the ND,eff profile (E5) that fit the calibration curves
within the 95 % threshold, and predicts for each of them
  * I-V of a device with a thicker channel (t_IWO = 7 nm), VD = 0.05 / 1.0 V
  * I-V of a device with a thicker gate oxide (t_ox = 14 nm), VD = 0.05 / 1.0 V
  * quasi-static gate C-V of the calibration device (VD = 0)
and compares with the prediction of the Step-5 truth: RMS deviation and the chi2
those new data would add (sigma = 0.02 dec for I-V, 0.5 % of C_ox A for C-V).
chi2 >> 3.84 means the new measurement rejects that set.

Usage:  python scripts/e9_discrimination_test.py results/e5_nd_down results/e5_nd_down2 results/e5_nd_up
            [--out results/e9]
"""
from __future__ import annotations

import argparse
import json

from _common import Timer, dump_json, np, outdir, vg_grid
from tftsim import PAPER_SETS, TFT, Geometry
from tftsim.cv import GateCharge
from tftsim.metrics import I_FLOOR, extract
from tftsim.params import with_theta


def params_of(theta: dict):
    p = PAPER_SETS["step5"]
    for n, v in theta.items():
        p = with_theta(p, n, float(v))
    return p


def iv(params, geom, vgs, vds=(0.05, 1.0)):
    t = TFT(params, geom).build()
    out = {}
    for j, vd in enumerate(vds):
        seq = vgs if j % 2 == 0 else vgs[::-1]
        recs = sorted(t.transfer(vd, seq), key=lambda r: r["vg"])
        out[vd] = np.array([r["id"] for r in recs])
    return out


def cv(params, vgs):
    t = TFT(params, Geometry()).build()
    gq = GateCharge(t)
    t.goto(vd=0.0)
    t.goto(vg=vgs[0])
    q = []
    for vg in vgs:
        t.goto(vg=vg)
        q.append(gq.value())
    return np.diff(np.array(q)) / np.diff(vgs), gq.c_ox_area()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--sigma", type=float, default=0.02)
    ap.add_argument("--frac", type=float, default=0.005)
    ap.add_argument("--out", default="results/e9")
    ap.add_argument("--values", default=None, help="comma-separated log10 ND,eff values to run (default: all)")
    ap.add_argument("--combine", nargs="*", default=None, help="only merge the summary.json of earlier runs")
    args = ap.parse_args()
    out = outdir(args.out)
    if args.combine:
        res = {}
        for d in args.combine:
            res.update(json.load(open(f"{d}/summary.json"))["results"])
        report(res, args, out)
        return
    clock = Timer()
    recs = {}
    for d in args.dirs:
        for r in json.load(open(f"{d}/summary.json"))["records"]:
            if "theta" in r and r["chi2"] < 3.84:
                recs[round(r["value"], 4)] = r
    vgs_iv = vg_grid(-1.0, 2.5, 0.1)
    vgs_cv = vg_grid(-1.0, 2.5, 0.05)
    conds = {"t_IWO 7 nm": Geometry(t_iwo_nm=7.0), "t_ox 14 nm": Geometry(t_gox_nm=14.0)}
    truth_key = min(recs, key=lambda k: recs[k]["chi2"])
    ref = {c: iv(PAPER_SETS["step5"], g, vgs_iv) for c, g in conds.items()}
    c_ref, coxa = cv(PAPER_SETS["step5"], vgs_cv)
    print(f"reference predictions ready ({clock():.0f} s)")
    res = {}
    keep = None if args.values is None else {round(float(v), 4) for v in args.values.split(",")}
    for key in sorted(recs):
        r = recs[key]
        if key == truth_key or (keep is not None and key not in keep):
            continue
        p = params_of(r["theta"])
        row = {"chi2_calibration": r["chi2"], "holdout_worst_rms": None}
        if isinstance(r.get("holdout"), dict) and "error" not in r["holdout"]:
            row["holdout_worst_rms"] = max(q["rms_dec"] for q in r["holdout"].values())
        for c, g in conds.items():
            cur = iv(p, g, vgs_iv)
            dev = np.concatenate([np.log10(cur[vd][ref[c][vd] >= I_FLOOR]) - np.log10(ref[c][vd][ref[c][vd] >= I_FLOOR])
                                  for vd in ref[c]])
            a, b = extract(vgs_iv, ref[c][1.0]), extract(vgs_iv, cur[1.0])
            row[c] = {"rms_dec": float(np.sqrt(np.mean(dev ** 2))), "chi2": float(np.sum((dev / args.sigma) ** 2)),
                      "dVth_mV": 1e3 * (b["vth"] - a["vth"])}
        c_cur, _ = cv(p, vgs_cv)
        dc = (c_cur - c_ref) / coxa
        row["C-V"] = {"max_dev_pct": float(100 * np.max(np.abs(dc))), "chi2": float(np.sum((dc / args.frac) ** 2))}
        res[key] = row
        print(f"[ND = 10^{key}] " + "  ".join(f"{c}: χ² {row[c]['chi2']:.3g}" for c in list(conds) + ["C-V"])
              + f"  ({clock():.0f} s)")
    dump_json({"sigma_dec": args.sigma, "cv_frac": args.frac, "results": res}, f"{out}/summary.json")
    report(res, args, out)


def report(res, args, out):
    cols = ["t_IWO 7 nm", "t_ox 14 nm", "C-V"]
    lines = ["# E9: do the OED-recommended measurements reject the sets the paper's data cannot?", "",
             f"Re-fitted sets of the ND,eff profile (χ² < 3.84 on the calibration curves). χ² of the NEW data vs the truth "
             f"(I–V σ = {args.sigma} dec over VD = 0.05/1.0 V; C–V σ = {100 * args.frac:g} % of C_ox·A). > 3.84 → rejected.", "",
             "| log10 ND,eff | χ² calibration | worst hold-out RMS (dec) | " + " | ".join(f"χ² {c}" for c in cols)
             + " | ΔVth t_IWO 7 nm (mV) | ΔVth t_ox 14 nm (mV) |", "|---|---|---|" + "---|" * len(cols) + "---|---|"]
    for key, row in sorted(res.items(), key=lambda kv: float(kv[0])):
        ho = "–" if row["holdout_worst_rms"] is None else f"{row['holdout_worst_rms']:.4f}"
        lines.append(f"| {key} | {row['chi2_calibration']:.3g} | {ho} | " + " | ".join(f"{row[c]['chi2']:.3g}" for c in cols)
                     + f" | {row['t_IWO 7 nm']['dVth_mV']:+.1f} | {row['t_ox 14 nm']['dVth_mV']:+.1f} |")
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
