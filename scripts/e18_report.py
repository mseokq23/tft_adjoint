"""E18 report - does adding the drain-bias transfer curves to the calibration fix ND,eff?

Compares the ND,eff profiles obtained with the paper's calibration pair (L = 1 µm, VD = 0.05 / 1 V)
and with all five L = 1 µm drain biases fitted jointly (VD = 0.05 / 0.1 / 0.2 / 0.5 / 1 V):
  measured data   E12 (pair)   vs  results/e18_joint_meas (joint; free fit + profile points)
  synthetic twin  E5 (pair)    vs  results/e18_joint_twin (joint; truth = Step 5)
Criterion: the joint profile stays below the 95 % threshold (Δχ² < 3.84; σ_eff for measured data)
over at least a factor of 10 in ND,eff. No simulation; reads stored results only.

Note on the measured joint fit: the digitised VD = 0.1 / 0.2 / 0.5 V curves have almost no points
around threshold (overlapping markers), so they add mainly on-state points (see n_res_vd).
The L = 2 µm hold-out curves have a 0.4 V digitisation gap around 100 nA: the ΔVth against them
interpolates across the gap, so the spread between the sets is taken from the full simulated curves.

Usage:  python scripts/e18_report.py [--out results/e18]
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import json
import os

from _common import HAVE_MPL, ROOT, dump_json, np, outdir
from tftsim.metrics import extract

L2_KEYS = {"L2_VD0.05": 0.05, "L2_VD1": 1.0}


def load(p):
    with open(p if os.path.isabs(p) else os.path.join(ROOT, p)) as f:
        return json.load(f)


def load_meas_l2(path="data/meas.csv", set_="2W"):
    d = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(ROOT, path))):
        if r["set"] == set_ and float(r["L_um"]) == 2.0:
            d[float(r["vd"])].append((float(r["vg"]), float(r["id"])))
    return {vd: np.array(sorted(v)) for vd, v in d.items()}


def full_vth(s, meas_l2):
    """Vth (mV) of the stored full simulated curves, and the shift caused by sampling the model only at the
    measured VG points (what the ΔVth against the gapped measured curve sees)."""
    vg = np.array(s["vg"])
    full, shift = {}, {}
    for k, h in s.get("holdout", {}).items():
        if "error" in h or "curve" not in h:
            continue
        cur = np.array(h["curve"])
        full[k] = 1e3 * extract(vg, cur)["vth"]
        if k in L2_KEYS and L2_KEYS[k] in meas_l2:
            m = meas_l2[L2_KEYS[k]]
            shift[k] = 1e3 * extract(m[:, 0], np.interp(m[:, 0], vg, cur))["vth"] - full[k]
    for vd, c in s.get("curves", {}).items():
        full["cal_" + vd] = 1e3 * extract(vg, np.array(c))["vth"]
    return full, shift


def gap_around(m, ith=1e-7):
    """Measured points that bracket ith: (VG_lo, VG_hi, I_lo, I_hi)."""
    j = int(np.where((m[:-1, 1] < ith) & (m[1:, 1] >= ith))[0][0])
    return m[j, 0], m[j + 1, 0], m[j, 1], m[j + 1, 1]


def spread(vals):
    vals = [v for v in vals if v is not None and np.isfinite(v)]
    return float(max(vals) - min(vals)) if vals else float("nan")


def measured_profile(folder, best_files, meas_l2):
    """Rows (log10 ND,eff, χ², summary) of a measured profile folder plus its free fit(s)."""
    rows = []
    for f in best_files:
        s = load(f)
        rows.append((s["theta_hat"]["log10_nd_eff"], s, "free fit"))
    for f in sorted(glob.glob(os.path.join(ROOT, folder, "log10_nd_eff_*", "summary.json"))):
        s = load(f)
        rows.append((s["fixed"]["log10_nd_eff"], s, "fixed"))
    rows.sort(key=lambda r: r[0])
    ref = min(rows, key=lambda r: r[1]["chi2"])[1]                       # lowest χ² over all fits
    best_free = min((r[1] for r in rows if r[2] == "free fit"), key=lambda s: s["chi2"])
    scale = (best_free["sigma_eff"] / best_free["sigma_dec"]) ** 2
    out = []
    for v, s, kind in rows:
        vth_full, shift = full_vth(s, meas_l2)
        out.append({"log10_nd_eff": v, "kind": kind, "chi2": s["chi2"], "n_res": s["n_res"],
                    "dchi2_eff": (s["chi2"] - ref["chi2"]) / scale, "rms_dec": s["rms_dec"],
                    "n_res_vd": s.get("n_res_vd"), "nfev": s.get("nfev"), "status": s.get("status"),
                    "theta": s["theta_hat"],
                    "holdout_dVth_mV": {k: h["dVth_mV"] for k, h in s.get("holdout", {}).items() if "error" not in h},
                    "cal_dVth_mV": {k: c["dVth_mV"] for k, c in s.get("cal", {}).items()},
                    "vth_full_mV": vth_full, "sampling_shift_mV": shift})
    return {"rows": out, "chi2_min": ref["chi2"], "sigma_eff": best_free["sigma_eff"], "sigma_dec": best_free["sigma_dec"],
            "vds": best_free.get("vds", [0.05, 1.0])}


def twin_profile(path):
    s = load(path)
    rows = [{"log10_nd_eff": r["value"], "chi2": r["chi2"], "rms_dec": r["rms_dec"], "nfev": r.get("nfev"),
             "status": r.get("status"), "at_bound": r.get("at_bound"), "theta": r["theta"],
             "holdout_dVth_mV": {k: h["dVth_mV"] for k, h in r.get("holdout", {}).items()
                                 if isinstance(h, dict) and "dVth_mV" in h}}
            for r in s["records"]]
    return {"rows": sorted(rows, key=lambda r: r["log10_nd_eff"]), "vds": s.get("vds", [0.05, 1.0]),
            "n_res": s.get("n_res")}


def equivalent_range(rows, key, thr=3.84):
    ok = [r["log10_nd_eff"] for r in rows if r[key] < thr]
    return (min(ok), max(ok), 10 ** (max(ok) - min(ok))) if ok else (None, None, 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meas_joint", default="results/e18_joint_meas")
    ap.add_argument("--twin_joint", default="results/e18_joint_twin/summary.json")
    ap.add_argument("--out", default="results/e18")
    args = ap.parse_args()
    out = outdir(args.out)
    meas_l2 = load_meas_l2()
    res = {}
    lines = ["# E18: does fitting the drain-bias transfer curves as well fix ND,eff?", "",
             "Calibration data: the paper's pair (L = 1 µm, VD = 0.05 / 1 V) vs all five L = 1 µm drain biases fitted jointly "
             "(VD = 0.05 / 0.1 / 0.2 / 0.5 / 1 V). Criterion: below the 95 % threshold (3.84) over at least 10× in ND,eff.", ""]

    # measured
    mj_free = os.path.join(args.meas_joint, "free", "summary.json")
    if os.path.exists(os.path.join(ROOT, mj_free)):
        mj = measured_profile(os.path.join(args.meas_joint, "profile"), [mj_free], meas_l2)
        mp = measured_profile("results/e12_nd", ["results/e11_2W_v3/summary.json", "results/e11_2W_v2/summary.json"], meas_l2)
        res["measured_joint"], res["measured_pair"] = mj, mp
        lo, hi, fac = equivalent_range(mj["rows"], "dchi2_eff")
        lo2, hi2, fac2 = equivalent_range(mp["rows"], "dchi2_eff")
        res["measured_joint_range"] = {"lo": lo, "hi": hi, "factor": fac}
        sp = {name: {k: spread([r["vth_full_mV"].get(k) for r in d["rows"]]) for k in keys}
              for name, d, keys in (("joint", mj, ("L2_VD0.05", "L2_VD1", "cal_0.1", "cal_0.2", "cal_0.5")),
                                    ("pair", mp, ("L2_VD0.05", "L2_VD1", "L1_VD0.1", "L1_VD0.2", "L1_VD0.5")))}
        shifts = [v for d in (mj, mp) for r in d["rows"] for v in r["sampling_shift_mV"].values()]
        res["vth_spread_full_mV"] = sp
        res["l2_sampling_shift_mV"] = [min(shifts), max(shifts)]
        nvd = mj["rows"][0]["n_res_vd"]
        lines += ["## Measured (digitised) data", "",
                  f"Joint fit: {mj['rows'][0]['n_res']} residuals {nvd}; χ²_min = {mj['chi2_min']:.1f}, "
                  f"σ_eff = {mj['sigma_eff']:.4f} dec. Pair fit (E12): χ²_min = {mp['chi2_min']:.1f}, σ_eff = {mp['sigma_eff']:.4f} dec.", "",
                  "| log10 ND,eff | kind | joint: Δχ²_eff | joint: RMS (dec) | joint: nfev | L 2 µm predicted Vth VD 0.05 / 1 V "
                  "(mV, full curve) | L 2 µm ΔVth VD 0.05 / 1 V vs measured (mV, across the gap) |",
                  "|---|---|---|---|---|---|---|"]
        for r in mj["rows"]:
            h, f_ = r["holdout_dVth_mV"], r["vth_full_mV"]
            lines.append(f"| {r['log10_nd_eff']:.3f} | {r['kind']} | {r['dchi2_eff']:.2f} | {r['rms_dec']:.4f} | {r['nfev']} | "
                         f"{f_.get('L2_VD0.05', float('nan')):.1f} / {f_.get('L2_VD1', float('nan')):.1f} | "
                         f"{h.get('L2_VD0.05', float('nan')):+.1f} / {h.get('L2_VD1', float('nan')):+.1f} |")
        lines += ["", f"Equivalent at 95 % (joint): log10 ND,eff {lo:.2f}–{hi:.2f} ({fac:.0f}× range). "
                  f"Pair fit (E12): {lo2:.2f}–{hi2:.2f} ({fac2:.0f}×).", "",
                  "Spread of the predicted Vth between the sets (full simulated curves): joint sets "
                  + ", ".join(f"{k} {v:.1f}" for k, v in sp["joint"].items())
                  + " mV; pair sets (E12) " + ", ".join(f"{k} {v:.1f}" for k, v in sp["pair"].items()) + " mV "
                  "(cal_*: fitted in the joint calibration; L1_*: hold-out tests of the pair fit).", "",
                  "The digitised L = 2 µm curves have no points around 100 nA: "
                  + "; ".join("VD = {:g} V: {:.1e} A at {:.2f} V to {:.1e} A at {:.2f} V".format(vd, g[2], g[0], g[3], g[1])
                              for vd, g in ((vd, gap_around(meas_l2[vd])) for vd in sorted(meas_l2)))
                  + ". The ΔVth against them samples the model at the measured VG and interpolates both curves across "
                  "the gap; sampling alone moves the model Vth by "
                  f"{min(shifts):+.0f} to {max(shifts):+.0f} mV, so the absolute ΔVth is uncertain by about this much.", ""]

    # twin
    if os.path.exists(os.path.join(ROOT, args.twin_joint)):
        tj = twin_profile(args.twin_joint)
        tp = twin_profile("results/e5_nd/summary.json")
        res["twin_joint"], res["twin_pair"] = tj, tp
        lo, hi, fac = equivalent_range(tj["rows"], "chi2")
        lo2, hi2, fac2 = equivalent_range(tp["rows"], "chi2")
        res["twin_joint_range"] = {"lo": lo, "hi": hi, "factor": fac}
        pair = {round(r["log10_nd_eff"], 3): r for r in tp["rows"]}
        lines += ["## Synthetic twin (truth = Step 5, σ = 0.02 dec, perfect model)", "",
                  f"Joint data: {tj['n_res']} residuals (VD = {'/'.join(f'{v:g}' for v in tj['vds'])} V).", "",
                  "| log10 ND,eff | joint: χ² | joint: RMS (dec) | joint: nfev | pair (E5): χ² | joint: L 2 µm ΔVth VD 0.05 / 1 V vs truth (mV) |",
                  "|---|---|---|---|---|---|"]
        for r in tj["rows"]:
            p = pair.get(round(r["log10_nd_eff"], 3))
            h = r["holdout_dVth_mV"]
            pc = f"{p['chi2']:.3g}" if p else "–"
            l2 = (f"{h['L2_VD0.05']:+.2f} / {h['L2_VD1']:+.2f}" if "L2_VD0.05" in h and "L2_VD1" in h
                  else "– (truth)" if r["chi2"] == 0 else "–")
            lines.append(f"| {r['log10_nd_eff']:.3f} | {r['chi2']:.3g} | {r['rms_dec']:.4f} | {r['nfev']} | {pc} | {l2} |")
        lines += ["", f"Equivalent at 95 % (joint): log10 ND,eff {lo:.2f}–{hi:.2f} ({fac:.0f}× range). "
                  f"Pair (E5): {lo2:.2f}–{hi2:.2f} ({fac2:.0f}×).", ""]

    dump_json(res, os.path.join(out, "summary.json"))
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    if HAVE_MPL and res:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
        if "measured_joint" in res:
            for key, lab, mk in (("measured_pair", "pair (VD 0.05/1 V)", "o-"), ("measured_joint", "joint (5 VD)", "s-")):
                rr = res[key]["rows"]
                ax[0].plot([r["log10_nd_eff"] for r in rr], [r["dchi2_eff"] for r in rr], mk, label=lab)
            ax[0].axhline(3.84, color="k", ls=":", lw=1)
            ax[0].text(ax[0].get_xlim()[0], 3.9, " 95 % threshold (3.84)", fontsize=7, va="bottom")
            ax[0].set_ylim(-0.2, 4.6)
            ax[0].set_xlabel("log10 ND,eff (cm⁻³), fixed")
            ax[0].set_ylabel("Δχ²_eff (measured, σ_eff)")
            ax[0].set_title("(a) measured data", loc="left")
            ax[0].legend(fontsize=8)
        if "twin_joint" in res:
            for key, lab, mk in (("twin_pair", "pair (VD 0.05/1 V)", "o-"), ("twin_joint", "joint (5 VD)", "s-")):
                rr = res[key]["rows"]
                ax[1].plot([r["log10_nd_eff"] for r in rr], [max(r["chi2"], 1e-4) for r in rr], mk, label=lab)
            ax[1].axhline(3.84, color="k", ls=":", lw=1)
            ax[1].text(15.85, 4.5, "95 % threshold (3.84)", fontsize=7, va="bottom")
            ax[1].text(18.45, 1.2e-4, "truth (χ² = 0)", fontsize=7, va="bottom")
            ax[1].set_yscale("log")
            ax[1].set_xlabel("log10 ND,eff (cm⁻³), fixed")
            ax[1].set_ylabel("χ² (twin, σ = 0.02 dec)")
            ax[1].set_title("(b) synthetic twin", loc="left")
            ax[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "e18_profile.png"), dpi=150)


if __name__ == "__main__":
    main()
