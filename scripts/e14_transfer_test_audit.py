"""E14 - does passing the paper's transfer tests say anything about the parameters?

A test can discriminate between parameter sets only if the sets predict different outcomes for it.
For every test this script compares the spread of the predicted outcome over the calibration-
equivalent parameter sets with the size of the errors the test is judged by. No new simulations.

1. The paper's transfer tests in the paper's own metrics - avg/max |ΔVth| and |Δlog Ion| over the
   drain-bias split (L = 1 µm, VD = 0.05/0.1/0.2/0.5/1.0 V) and the length split (L = 2 µm,
   VD = 0.05/1.0 V):
   * synthetic twin (E5): every re-fitted set against the truth (Step 5), i.e. a perfect model;
   * measured data (E12 + the free fit): every re-fitted set against the digitised measurements;
   next to the values the paper reports for its Sentaurus model.
2. Spread of ΔVth across the equivalent sets (max - min) for every test: the paper's transfer tests
   against the OED-recommended tests (E9 synthetic, E13 measured, if present).
3. Linear version at the Step-5 truth (sensitivity rows of E3/E4/E7, Qf fixed): predictive 1σ of each
   test under the calibration-data posterior, σ_pred = sqrt(diag(S_test F_cal^-1 S_test^T)), relative
   to the noise of that test. σ_pred/σ << 1 means every calibration-consistent set passes the test, so
   the test carries no information about the parameters. (The log-det gain used in E3/E7 is
   Σ log(1 + λ_i(S_test F_cal^-1 S_test^T)/σ²), the same quantity aggregated.)

Usage:  python scripts/e14_transfer_test_audit.py [--out results/e14]
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import json
import os

from _common import HAVE_MPL, dump_json, np, outdir
from tftsim.metrics import extract

PAPER = {"drain": {"vth_avg": 24.9, "vth_max": 36.1, "ion_avg": 0.013, "ion_max": 0.022},
         "length": {"vth_avg": 27.4, "vth_max": 46.2, "ion_avg": 0.050, "ion_max": 0.062}}
DRAIN = ["cal_0.05", "L1_VD0.1", "L1_VD0.2", "L1_VD0.5", "cal_1"]
LENGTH = ["L2_VD0.05", "L2_VD1"]
HOLD = ["L1_VD0.1", "L1_VD0.2", "L1_VD0.5", "L2_VD0.05", "L2_VD1"]


def norm_key(k):
    tag, vd = k.split("_VD")
    return f"{tag}_VD{float(vd):g}"


def aggregate(err):
    """err: {condition: (dVth_mV, dlogIon_dec)} -> paper-style aggregates per split."""
    out = {}
    for split, keys in (("drain", DRAIN), ("length", LENGTH)):
        v = [abs(err[k][0]) for k in keys if k in err and np.isfinite(err[k][0])]
        i = [abs(err[k][1]) for k in keys if k in err and np.isfinite(err[k][1])]
        out[split] = {"vth_avg": float(np.mean(v)), "vth_max": float(np.max(v)),
                      "ion_avg": float(np.mean(i)), "ion_max": float(np.max(i)), "n": len(v)}
    return out


def synthetic_sets(path="results/e5_nd/summary.json"):
    s = json.load(open(path))
    vg = np.array(s["vg"])
    vt = float(s["v_true"])
    curves = {float(k): v for k, v in s["curves"].items()}
    truth = curves[min(curves, key=lambda k: abs(k - vt))]
    sets = {}
    for r in s["records"]:
        if "theta" not in r or not isinstance(r.get("holdout"), dict) or "error" in r["holdout"]:
            continue
        err = {norm_key(k): (h["dVth_mV"], h["dlogIon_dec"]) for k, h in r["holdout"].items()}
        cur = curves.get(min(curves, key=lambda k: abs(k - r["value"])))
        for vd, tag in (("0.05", "cal_0.05"), ("1.0", "cal_1")):
            a, b = extract(vg, np.array(truth[vd])), extract(vg, np.array(cur[vd]))
            err[tag] = (1e3 * (b["vth"] - a["vth"]), b["log10_ion"] - a["log10_ion"])
        sets[r["value"]] = {"chi2": r["chi2"], "err": err}
    return sets, vt


def load_meas_curves(path="data/meas.csv"):
    d = collections.defaultdict(list)
    for r in csv.DictReader(open(path)):
        d[(r["set"], float(r["L_um"]), float(r["vd"]))].append((float(r["vg"]), float(r["id"])))
    return {k: np.array(sorted(v)) for k, v in d.items()}


def elog_above(s, meas, floor=1e-9):
    """RMS log10 error of each hold-out prediction against the measured points above `floor`
    (the digitised curves reach the instrument floor around 1e-10 A; isolated floor points dominate an RMS)."""
    vg = np.array(s["vg"])
    out = {}
    for k, h in s["holdout"].items():
        if "error" in h or "curve" not in h:
            continue
        tag, vd = norm_key(k).split("_VD")
        m = meas[("2W", 2.0 if tag == "L2" else 1.0, float(vd))]
        ok = m[:, 1] >= floor
        sim = np.interp(m[ok, 0], vg, np.array(h["curve"]))
        out[norm_key(k)] = float(np.sqrt(np.mean((np.log10(sim) - np.log10(m[ok, 1])) ** 2)))
    return out


def measured_sets(profile="results/e12_nd", best="results/e11_2W_v3/summary.json,results/e11_2W_v2/summary.json"):
    files = best.split(",") + sorted(glob.glob(os.path.join(profile, "log10_nd_eff_*", "summary.json")))
    raw = [json.load(open(f)) for f in files]
    b = min(raw[:len(best.split(","))], key=lambda s: s["chi2"])
    scale = (b["sigma_eff"] / b["sigma_dec"]) ** 2
    chi_min = min(s["chi2"] for s in raw)
    meas = load_meas_curves()
    sets = {}
    for s in raw:
        err = {norm_key(k): (h["dVth_mV"], h["dlogIon_dec"]) for k, h in s["holdout"].items() if "error" not in h}
        for vd, tag in (("0.05", "cal_0.05"), ("1.0", "cal_1")):
            err[tag] = (s["cal"][vd]["dVth_mV"], s["cal"][vd]["dlogIon_dec"])
        sets[s["theta_hat"]["log10_nd_eff"]] = {"chi2": s["chi2"], "dchi2_eff": (s["chi2"] - chi_min) / scale,
                                                "rms_dec": s["rms_dec"], "J": s.get("J"), "err": err,
                                                "elog_1nA": elog_above(s, meas)}
    return sets


def spread(vals):
    vals = [v for v in vals if v is not None and np.isfinite(v)]
    return float(max(vals) - min(vals)) if len(vals) > 1 else float("nan")


def linear_predictive(sigma=0.02, cv_frac=0.005):
    """Predictive 1σ of each test under the calibration posterior (Qf fixed), at the Step-5 truth."""
    z3, z4, z7 = (np.load(f"results/{e}/rows.npz") for e in ("e3", "e4", "e7"))
    names = [str(n) for n in z3["names"]]
    keep = [k for k, n in enumerate(names) if n != "qf"]
    S0 = z3["S_base"][:, keep]
    Finv = np.linalg.inv(S0.T @ S0 / sigma ** 2)
    tests = {"VD 0.1 V (paper)": z3["S_L1_VD0.1"], "VD 0.2 V (paper)": z3["S_L1_VD0.2"], "VD 0.5 V (paper)": z3["S_L1_VD0.5"],
             "L 2 µm, VD 0.05 V (paper)": z3["S_L2_VD0.05"], "L 2 µm, VD 1 V (paper)": z3["S_L2_VD1.0"],
             "t_ox 10 nm": z4["S_tox10"], "t_ox 14 nm": z4["S_tox14"], "t_IWO 5 nm": z7["S_tch5"], "t_IWO 7 nm": z7["S_tch7"]}
    out = {}
    for name, S in tests.items():
        S = S[:, keep]
        v = np.einsum("ij,jk,ik->i", S, Finv, S)
        out[name] = {"rms_ratio": float(np.sqrt(np.mean(v)) / sigma), "max_ratio": float(np.sqrt(np.max(v)) / sigma),
                     "noise": f"{sigma} dec", "dlogdet": float(np.linalg.slogdet(np.eye(len(keep)) + Finv @ S.T @ S / sigma ** 2)[1])}
    Scv = z7["S_cv_unit"][:, keep]
    coxa = float(z7["coxa"])
    v = np.einsum("ij,jk,ik->i", Scv, Finv, Scv)
    sc = cv_frac * coxa
    out["quasi-static C–V"] = {"rms_ratio": float(np.sqrt(np.mean(v)) / sc), "max_ratio": float(np.sqrt(np.max(v)) / sc),
                               "noise": f"{100 * cv_frac:g} % C_ox·A",
                               "dlogdet": float(np.linalg.slogdet(np.eye(len(keep)) + Finv @ Scv.T @ Scv / sc ** 2)[1])}
    return out


def measured_ss(path="data/meas.csv"):
    d = collections.defaultdict(list)
    for r in csv.DictReader(open(path)):
        if r["set"] == "2W" and float(r["L_um"]) == 1.0:
            d[float(r["vd"])].append((float(r["vg"]), float(r["id"])))
    out = {}
    for vd in (0.05, 1.0):
        a = np.array(sorted(d[vd]))
        out[vd] = extract(a[:, 0], a[:, 1])["ss_mVdec"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e14")
    args = ap.parse_args()
    out = outdir(args.out)
    syn, vt = synthetic_sets()
    syn_eq = {v: s for v, s in syn.items() if s["chi2"] < 3.84}
    mea = measured_sets()
    ss = measured_ss()
    sig_vth = {vd: 0.02 * s for vd, s in ss.items()}            # digitisation noise expressed in Vth (mV)
    lines = ["# E14: the paper's transfer tests cannot tell calibration-equivalent parameter sets apart", ""]
    res = {"paper": PAPER, "synthetic": {}, "measured": {}, "spread": {}, "linear": {}}

    # 1. paper-style aggregates
    lines += ["## 1. Transfer tests in the paper's metrics", "",
              "Drain-bias split: L = 1 µm, VD = 0.05/0.1/0.2/0.5/1.0 V. Length split: L = 2 µm, VD = 0.05/1.0 V. "
              "avg/max of |ΔVth| (mV) and |Δlog10 Ion| (dec).", "",
              "| model / set | log10 ND,eff | fit to calibration curves | drain: Vth avg/max | drain: Ion avg/max | length: Vth avg/max | length: Ion avg/max |",
              "|---|---|---|---|---|---|---|"]
    p = PAPER
    lines.append(f"| paper (Sentaurus, Step 5) | 18.41 | J = 0.181 | {p['drain']['vth_avg']}/{p['drain']['vth_max']} | "
                 f"{p['drain']['ion_avg']}/{p['drain']['ion_max']} | {p['length']['vth_avg']}/{p['length']['vth_max']} | "
                 f"{p['length']['ion_avg']}/{p['length']['ion_max']} |")
    for v, s in sorted(syn.items()):
        a = aggregate(s["err"])
        res["synthetic"][f"{v:.2f}"] = {"chi2": s["chi2"], **a}
        lines.append(f"| twin: re-fit vs truth | {v:.2f} | χ² = {s['chi2']:.3g} | {a['drain']['vth_avg']:.2f}/{a['drain']['vth_max']:.2f} | "
                     f"{a['drain']['ion_avg']:.4f}/{a['drain']['ion_max']:.4f} | {a['length']['vth_avg']:.2f}/{a['length']['vth_max']:.2f} | "
                     f"{a['length']['ion_avg']:.4f}/{a['length']['ion_max']:.4f} |")
    for v, s in sorted(mea.items()):
        a = aggregate(s["err"])
        res["measured"][f"{v:.3f}"] = {"dchi2_eff": s["dchi2_eff"], "rms_dec": s["rms_dec"], **a}
        lines.append(f"| ours vs measured | {v:.2f} | Δχ²_eff = {s['dchi2_eff']:.2f} (RMS {s['rms_dec']:.4f} dec) | "
                     f"{a['drain']['vth_avg']:.1f}/{a['drain']['vth_max']:.1f} | {a['drain']['ion_avg']:.3f}/{a['drain']['ion_max']:.3f} | "
                     f"{a['length']['vth_avg']:.1f}/{a['length']['vth_max']:.1f} | {a['length']['ion_avg']:.3f}/{a['length']['ion_max']:.3f} |")
    crit = {v: 0.02 * np.sqrt(s["chi2"] / 3.84) for v, s in syn_eq.items() if s["chi2"] > 0}
    res["synthetic_noise_needed_dec"] = {f"{v:.2f}": c_ for v, c_ in crit.items()}
    lines += ["", "Twin: noise level at which the calibration curves alone would reject each set (χ² = 3.84, perfect model): "
              + ", ".join(f"10^{v:.2f}: {c_:.1e} dec" for v, c_ in sorted(crit.items())) + "."]
    agg_m = [aggregate(s["err"]) for s in mea.values()]
    rng = {f"{sp}_{m}": spread([a[sp][m] for a in agg_m]) for sp in ("drain", "length") for m in ("vth_avg", "vth_max", "ion_avg", "ion_max")}
    res["measured_aggregate_spread"] = rng
    lines += ["", f"Measured sets: spread of the aggregates across ND,eff = 10^{min(mea):.2f}–10^{max(mea):.2f}: "
              + ", ".join(f"{k} {v:.2g}" for k, v in rng.items()) + " (mV or dec).", ""]

    # 2. spread across the equivalent sets, per test
    js = [s["J"] for s in mea.values() if s.get("J") is not None]
    lines += [f"The paper's objective J (paper weights, our normalisation a_i) over the same measured sets: {min(js):.3f}–{max(js):.3f}.", ""]
    res["measured_J_range"] = [min(js), max(js)]
    lines += ["## 2. Spread across the equivalent sets (max − min)", "",
              f"Digitisation noise expressed in Vth: σ_logI × SS = 0.02 dec × {ss[0.05]:.0f}/{ss[1.0]:.0f} mV/dec = "
              f"{sig_vth[0.05]:.1f}/{sig_vth[1.0]:.1f} mV (VD = 0.05/1 V). The paper's own transfer errors: 25–46 mV. "
              "Measured E_log: RMS log10 error against the measured points above 1 nA (range over the sets, and its spread).", "",
              "| test | in the paper's split | twin: ΔVth spread over ND,eff 10^15.9–10^18.9 (mV) | measured: ΔVth spread over ND,eff "
              f"10^{min(mea):.2f}–10^{max(mea):.2f} (mV) | measured: E_log above 1 nA (dec) |",
              "|---|---|---|---|---|"]
    for k in HOLD + ["cal_0.05", "cal_1"]:
        a = spread([s["err"][k][0] for s in syn_eq.values() if k in s["err"]] + [0.0])
        b = spread([s["err"][k][0] for s in mea.values() if k in s["err"]])
        el = [s["elog_1nA"][k] for s in mea.values() if k in s["elog_1nA"]]
        res["spread"][k] = {"paper_split": True, "synthetic_mV": a, "measured_mV": b,
                            "measured_elog_range": [min(el), max(el)] if el else None}
        name = {"cal_0.05": "calibration VD 0.05 V", "cal_1": "calibration VD 1 V"}.get(k, k.replace("_", " "))
        elt = f"{min(el):.3f}–{max(el):.3f} (spread {max(el) - min(el):.3f})" if el else "(fitted)"
        lines.append(f"| {name} | yes | {a:.2f} | {b:.1f} | {elt} |")
    e9 = json.load(open("results/e9/summary.json"))["results"] if os.path.exists("results/e9/summary.json") else None
    if e9 is None:
        parts = [json.load(open(f"results/{d}/summary.json"))["results"] for d in ("e9a", "e9b") if os.path.exists(f"results/{d}/summary.json")]
        e9 = {k: v for p_ in parts for k, v in p_.items()}
    e13 = json.load(open("results/e13/summary.json")) if os.path.exists("results/e13/summary.json") else None
    for c in ("t_IWO 7 nm", "t_ox 14 nm"):
        a = spread([r[c]["dVth_mV"] for r in e9.values()] + [0.0])
        b = e13["tests"].get(f"{c} VD=1", {}).get("vth_spread_mV", float("nan")) if e13 else float("nan")
        rm = e13["tests"].get(f"{c} VD=1", {}).get("max_rms_dec", float("nan")) if e13 else float("nan")
        res["spread"][c] = {"paper_split": False, "synthetic_mV": a, "measured_mV": b, "measured_max_rms_between_sets_dec": rm}
        bt = f"{b:.0f}" if np.isfinite(b) else "(E13 pending)"
        lines.append(f"| {c}, VD 1 V | no (recommended) | {a:.0f} | {bt} | predictions differ by up to {rm:.2f} dec RMS |")
    a = max(r["C-V"]["max_dev_pct"] for r in e9.values())
    b = e13["tests"]["C-V"]["max_dev_pct"] if e13 else float("nan")
    res["spread"]["C-V"] = {"paper_split": False, "synthetic_max_dev_pct": a, "measured_max_dev_pct": b}
    lines.append(f"| quasi-static C–V (max ΔC, % of C_ox·A) | no (recommended) | {a:.1f} % | " + (f"{b:.1f} % | – |" if np.isfinite(b) else "(E13 pending) | – |"))
    lines += ["", "Twin: ΔVth of each re-fitted set relative to the truth (truth included as 0). Measured: ΔVth of each set against the "
              "measured curve; the spread is what can differ between the sets. For tests outside the paper's split no data exist: the "
              "spread is between the sets' predictions (E9 vs the truth, E13 between the measured-fit sets).", ""]

    # 3. linear predictive spread
    lin = linear_predictive()
    res["linear"] = lin
    lines += ["## 3. Linear check at the Step-5 truth: predictive 1σ of each test given the calibration data", "",
              "σ_pred/σ ≪ 1: every calibration-consistent parameter set predicts the same outcome within noise, so passing the "
              "test is automatic. Qf fixed (exactly degenerate with the gate WF). Calibration noise 0.02 dec.", "",
              "| test | noise σ | σ_pred/σ (RMS over the curve) | σ_pred/σ (max) | Δlog det F |", "|---|---|---|---|---|"]
    for name, r in lin.items():
        lines.append(f"| {name} | {r['noise']} | {r['rms_ratio']:.2f} | {r['max_ratio']:.2f} | {r['dlogdet']:.1f} |")
    # 4. robustness of the measured profile (other parameters, material constants, the paper's 4%-W values)
    rob = robustness(mea)
    res["robustness"] = rob
    if rob:
        lines += ["", "## 4. Robustness of the measured-data equivalence", "",
                  f"Reference χ²_min = {rob['chi2_min']:.1f} (lowest default-material fit), σ_eff = {rob['sigma_eff']:.4f} dec. "
                  "Fits are limited to 25–40 evaluations, so each χ² is an upper bound of the minimum under its constraints.", "",
                  "| check | fixed | free | χ² | Δχ²_eff | equivalent at 95 %? |", "|---|---|---|---|---|---|"]
        for r in rob["rows"]:
            lines.append(f"| {r['check']} | {r['fixed']} | {r['free']} | {r['chi2']:.1f} | {r['dchi2_eff']:.2f} | "
                         f"{'yes' if r['dchi2_eff'] < 3.84 else 'no'} |")
    dump_json(res, os.path.join(out, "summary.json"))
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def robustness(mea):
    b = json.load(open("results/e11_2W_v3/summary.json")) if os.path.exists("results/e11_2W_v3/summary.json") else None
    if b is None:
        return {}
    scale = (b["sigma_eff"] / b["sigma_dec"]) ** 2
    default_runs = glob.glob("results/e15_wta/*/summary.json") + glob.glob("results/e17_cross/*/summary.json")
    chi_min = min([s["chi2"] for s in mea.values()] + [json.load(open(f))["chi2"] for f in default_runs])
    rows = []

    def row(check, f, ref=chi_min, sc=scale):
        s = json.load(open(f))
        fixed = ", ".join(f"{k} = {v:g}" for k, v in s["fixed"].items())
        rows.append({"check": check, "file": f, "fixed": fixed or "–", "free": ", ".join(s["free"]), "chi2": s["chi2"],
                     "dchi2_eff": (s["chi2"] - ref) / sc, "status": s["status"]})

    for f in sorted(glob.glob("results/e15_wta/*/summary.json")):
        row("tail width at the paper's 4%-W value (2%-W data)", f)
    nc = sorted(glob.glob("results/e16_nc2e18/*/summary.json"))
    if nc:
        ss = [json.load(open(f)) for f in nc]
        ref = min(s["chi2"] for s in ss)
        sc = (min(ss, key=lambda s: s["chi2"])["sigma_eff"] / ss[0]["sigma_dec"]) ** 2
        for f in nc:
            row("Nc = 2e18 cm⁻³ instead of 5e18 (vs the better of the two ends)", f, ref, sc)
    for f, name in (("results/e17_cross/dos_w4/summary.json", "paper's 4%-W DOS and ND,eff, only WF_sd/WF_gate/µ0 fitted"),
                    ("results/e17_cross/dos_step5/summary.json", "paper's 2%-W (Step 5) DOS and ND,eff, only WF_sd/WF_gate/µ0 fitted")):
        if os.path.exists(f):
            row(name, f)
    return {"chi2_min": chi_min, "sigma_eff": b["sigma_eff"], "rows": rows}


if __name__ == "__main__":
    main()
