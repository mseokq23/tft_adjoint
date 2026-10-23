"""E3c - is the identifiability structure of E3 specific to one parameter set / material assumption?

Repeats the calibration-condition Fisher analysis (L = 1 um, VD = 0.05/1.0 V,
sigma = 0.02 dec, ID >= 10 pA) for
  * the paper's 4 %-W recalibration (w4),
  * the final 2 %-W model with a different assumed conduction-band DOS (Nc = 2e18 cm^-3),
  * the final 2 %-W model with a different channel permittivity (eps = 15),
and compares with the reference (Step 5, rows from E3). Reported: the exact
Qf-WF null eigenvalue, the condition number without Qf, and the 1-sigma of
every parameter. sigma only rescales all standard deviations linearly, so it is
not varied here.

Usage:  python scripts/e3c_robustness.py [--cases w4,nc2e18,eps15] [--out results/e3c]
"""
from __future__ import annotations

import argparse
from dataclasses import replace

from _common import Timer, dump_json, np, outdir, vg_grid
from e3_identifiability_oed_preview import NAMES, log_rows
from tftsim import PAPER_SETS, TFT, Geometry
from tftsim.metrics import I_FLOOR
from tftsim.params import Material, theta_scale
from tftsim.sensitivity import SensitivityEngine

CASES = {"w4": (PAPER_SETS["w4"], Material()),
         "nc2e18": (PAPER_SETS["step5"], Material(Nc_cm3=2.0e18)),
         "eps15": (PAPER_SETS["step5"], Material(eps_ch=15.0))}


def sweep(params, mat, vgs, vds=(0.05, 1.0)):
    t = TFT(params, Geometry(), mat).build()
    eng = SensitivityEngine(t, central=False)
    data = {}
    for j, vd in enumerate(vds):
        seq = vgs if j % 2 == 0 else vgs[::-1]

        def cb(tft, rec):
            if rec["id"] >= 0.1 * I_FLOOR:
                rec["g"] = eng.current_gradient(NAMES, "adjoint")

        recs = sorted(t.transfer(vd, seq, on_point=cb), key=lambda r: r["vg"])
        data[vd] = (np.array([r["vg"] for r in recs]), np.array([r["id"] for r in recs]),
                    np.array([r.get("g", np.zeros(len(NAMES))) for r in recs]))
    return data


def analyse(S, sigma, scale):
    F = S.T @ S / sigma ** 2
    w = np.linalg.eigvalsh(F)
    idx = [k for k, n in enumerate(NAMES) if n != "qf"]
    Fr = F[np.ix_(idx, idx)]
    wr, Vr = np.linalg.eigh(Fr)
    sd = np.sqrt(np.diag(np.linalg.inv(Fr))) * scale[idx]
    kn = [NAMES[k] for k in idx]
    top = np.argsort(-np.abs(Vr[:, 0]))[:4]
    return {"null_rel": float(w[0] / w[-1]), "cond_noqf": float(wr[-1] / wr[0]), "n_rows": int(S.shape[0]),
            "std_phys": dict(zip(kn, sd.tolist())),
            "sloppiest": {kn[k]: float(Vr[k, 0]) for k in top}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default="w4,nc2e18,eps15")
    ap.add_argument("--sigma", type=float, default=0.02)
    ap.add_argument("--out", default="results/e3c")
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    scale = theta_scale(NAMES)
    vgs = vg_grid(vmin=-1.0, step=0.05)
    res = {"step5 (reference)": analyse(np.load("results/e3/rows.npz")["S_base"], args.sigma, scale)}
    for name in args.cases.split(","):
        p, mat = CASES[name]
        d = sweep(p, mat, vgs)
        S = np.vstack([log_rows(*d[vd]) for vd in (0.05, 1.0)]) * scale
        res[name] = analyse(S, args.sigma, scale)
        np.save(f"{out}/S_{name}.npy", S)
        print(f"[{name}] rows {S.shape[0]}, null {res[name]['null_rel']:.1e}, cond(no qf) {res[name]['cond_noqf']:.2e} ({clock():.0f} s)")
    dump_json({"sigma_dec": args.sigma, "results": res, "wall_s": clock()}, f"{out}/summary.json")
    keys = list(res)
    kn = list(res[keys[0]]["std_phys"])
    lines = ["# E3c: robustness of the identifiability structure (calibration conditions, σ = %g dec)" % args.sigma, "",
             "| quantity | " + " | ".join(keys) + " |", "|---|" + "---|" * len(keys),
             "| rows (ID ≥ 10 pA) | " + " | ".join(str(res[k]["n_rows"]) for k in keys) + " |",
             "| Qf–WF null λmin/λmax | " + " | ".join(f"{res[k]['null_rel']:.1e}" for k in keys) + " |",
             "| condition number (Qf removed) | " + " | ".join(f"{res[k]['cond_noqf']:.1e}" for k in keys) + " |"]
    for n in kn:
        lines.append(f"| 1σ {n} | " + " | ".join(f"{res[k]['std_phys'][n]:.3g}" for k in keys) + " |")
    lines += ["", "Sloppiest direction (Qf removed, scaled units):", ""]
    for k in keys:
        lines.append(f"- {k}: " + ", ".join(f"{n} {v:+.2f}" for n, v in res[k]["sloppiest"].items()))
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
