"""E5 - profile likelihood: is a parameter identifiable from the paper's calibration data?

The Fisher analysis of E3 is local and linear. Here the question is asked of the
full nonlinear model: fix one parameter (default log10 ND,eff) at a value away
from the truth, re-fit all other parameters to the calibration curves, and see
how much the best achievable misfit grows.

Data: noise-free synthetic curves of the final 2%-W model (Step 5) at
L = 1 um, VD = 0.05 and 1.0 V (the conditions the paper calibrates on), log10 ID
above 10 pA, sigma = 0.02 dec (digitisation-level noise). Qf is held at its true
value because it is exactly degenerate with the gate work function (E3).

Profile chi2 below 3.84 (95 %, 1 dof) over a range of values means the data
cannot distinguish them: the parameter is practically non-identifiable.

Fits use the adjoint Jacobian (trust-region reflective LM, bounds as in E2) and
continuation outward from the truth on each side.

Usage:  python scripts/e5_profile_likelihood.py [--param log10_nd_eff]
            [--values 16.5,17,17.5,18,18.8,19.2] [--step 0.1] [--max_nfev 20]
Joint calibration (E18): --vds 0.05,0.1,0.2,0.5,1.0 uses all five L = 1 µm drain biases as the
calibration data; the hold-out predictions then cover only the remaining conditions (L = 2 µm).
"""
from __future__ import annotations

import argparse
import json

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from scipy.optimize import least_squares

from e1_forward_reproduce import sweep_device
from e2_twin_calibration import BOUNDS
from tftsim import PAPER_SETS, THETA_ALL, ConvergenceError, Geometry
from tftsim.calib import CurveModel, log_residual_problem
from tftsim.metrics import I_FLOOR, extract
from tftsim.params import theta_scale, theta_vector, with_theta

HOLDOUT = [(1.0, (0.1, 0.2, 0.5)), (2.0, (0.05, 1.0))]     # (L_um, VDs): the paper's validation conditions


def holdout_curves(params, vgs, holdout=HOLDOUT):
    """Transfer curves of a parameter set at the paper's hold-out conditions."""
    out = {}
    for L, vds in holdout:
        _, cv = sweep_device(params, Geometry(L_um=L), vds=vds, vgs=vgs)
        for vd, (vg, i) in cv.items():
            out[f"L{L:g}_VD{vd:g}"] = i
    return out


def holdout_compare(ref, cur, vgs):
    """RMS / max log10 deviation (points above 10 pA in the reference) and metric differences."""
    rep = {}
    for k in ref:
        m = ref[k] >= I_FLOOR
        dev = np.log10(np.maximum(cur[k][m], 1e-30)) - np.log10(ref[k][m])
        a, b = extract(vgs, ref[k]), extract(vgs, cur[k])
        rep[k] = {"rms_dec": float(np.sqrt(np.mean(dev ** 2))), "max_dec": float(np.max(np.abs(dev))),
                  "dVth_mV": 1e3 * (b["vth"] - a["vth"]), "dSS_mVdec": b["ss_mVdec"] - a["ss_mVdec"],
                  "dlogIon_dec": b["log10_ion"] - a["log10_ion"]}
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--param", default="log10_nd_eff")
    ap.add_argument("--values", default="16.5,17.0,17.5,18.0,18.8,19.2")
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--sigma", type=float, default=0.02)
    ap.add_argument("--max_nfev", type=int, default=20)
    ap.add_argument("--vds", default="0.05,1.0", help="drain biases (V) of the L = 1 µm calibration curves")
    ap.add_argument("--out", default=None)
    ap.add_argument("--no_holdout", action="store_true", help="skip predictions at the paper's hold-out conditions")
    ap.add_argument("--init", default=None, help="summary.json of an earlier profile run: each side starts from its nearest record")
    args = ap.parse_args()
    out = outdir(args.out or f"results/e5_{args.param}")
    clock = Timer()
    prof = args.param
    free = [n for n in THETA_ALL if n not in ("qf", prof)]
    names = free + [prof]                       # profiled parameter last
    scale = theta_scale(names)
    truth = PAPER_SETS["step5"]
    th_true = theta_vector(truth, names)
    vgs = vg_grid(-1.0, 2.5, args.step)
    lb = np.array([(BOUNDS[n][0] - th_true[k]) / scale[k] for k, n in enumerate(free)])
    ub = np.array([(BOUNDS[n][1] - th_true[k]) / scale[k] for k, n in enumerate(free)])

    fit_vds = tuple(float(v) for v in args.vds.split(","))
    holdout = [(L, tuple(vd for vd in vds if not (L == 1.0 and vd in fit_vds))) for L, vds in HOLDOUT]
    holdout = [(L, vds) for L, vds in holdout if vds]
    model = CurveModel(truth, names, vds=fit_vds, vgs=vgs)
    meas, _ = model.evaluate(th_true)
    fun_all, jac_all, masks = log_residual_problem(model, meas, args.sigma, th_true, scale)
    n_res = sum(int(m.sum()) for m in masks.values())
    print(f"synthetic data: {n_res} points above 10 pA ({clock():.0f} s)")
    ho_ref = None        # hold-out predictions are made after all fits: building another device resets DEVSIM

    def params_of(th):
        p = truth
        for n, v in zip(names, th):
            p = with_theta(p, n, float(v))
        return p

    values = sorted(float(v) for v in args.values.split(","))
    v_true = th_true[-1]
    sides = [[v for v in values if v > v_true], [v for v in values if v < v_true][::-1]]
    records = [{"value": float(v_true), "cost": 0.0, "chi2": 0.0, "rms_dec": 0.0, "max_dec": 0.0,
                "theta": dict(zip(names, th_true.tolist())), "nfev": 0, "njev": 0, "wall_s": 0.0}]
    curves = {float(v_true): {vd: meas[vd].tolist() for vd in meas}}
    init_recs = [r for r in json.load(open(args.init))["records"] if "theta" in r] if args.init else []
    for side in sides:
        z0 = np.zeros(len(free))
        if side and init_recs:                     # continue from an earlier run instead of the truth
            near = min(init_recs, key=lambda r: abs(r["value"] - side[0]))
            z0 = np.clip((np.array([near["theta"][n] for n in free]) - th_true[:-1]) / scale[:-1], lb + 1e-6, ub - 1e-6)
            print(f"side starting at {side[0]}: initialised from the fit at {prof} = {near['value']:.3f}")
        for v in side:
            zp = (v - v_true) / scale[-1]
            fun = lambda z, zp=zp: fun_all(np.append(z, zp))
            jac = lambda z, zp=zp: jac_all(np.append(z, zp))[:, :-1]
            tm = Timer()
            n0 = model.t.n_solves
            try:
                res = least_squares(fun, z0, jac=jac, bounds=(lb, ub), method="trf", x_scale=1.0,
                                    max_nfev=args.max_nfev, ftol=1e-8, xtol=1e-8, gtol=1e-8)
            except Exception as exc:          # convergence failure far from the truth
                print(f"[{prof} = {v:.3f}] fit failed: {exc}")
                records.append({"value": v, "error": str(exc)})
                break
            th = np.append(th_true[:-1] + res.x * scale[:-1], v)
            I, _ = model.evaluate(th)
            dev = np.concatenate([np.log10(I[vd][masks[vd]]) - np.log10(meas[vd][masks[vd]]) for vd in meas])
            chi2 = float(np.sum((dev / args.sigma) ** 2))
            rec = {"value": v, "cost": float(res.cost), "chi2": chi2, "rms_dec": float(np.sqrt(np.mean(dev ** 2))),
                   "max_dec": float(np.max(np.abs(dev))), "theta": dict(zip(names, th.tolist())),
                   "nfev": int(res.nfev), "njev": int(res.njev or 0), "status": int(res.status),
                   "pde_solves": model.t.n_solves - n0, "wall_s": tm(),
                   "at_bound": [free[k] for k in range(len(free))
                                if res.active_mask[k] != 0 or min(res.x[k] - lb[k], ub[k] - res.x[k]) < 1e-3]}
            records.append(rec)
            curves[v] = {vd: I[vd].tolist() for vd in I}
            print(f"[{prof} = {v:.3f}] chi2 = {rec['chi2']:.3f}  rms = {rec['rms_dec']:.4f} dec  max = {rec['max_dec']:.4f} dec  "
                  f"nfev {rec['nfev']} njev {rec['njev']}  at-bound {rec['at_bound']}  ({tm():.0f} s, total {clock():.0f} s)")
            z0 = res.x.copy()                  # continuation
            dump_json({"param": prof, "names": names, "sigma_dec": args.sigma, "records": records}, f"{out}/partial.json")
    if not args.no_holdout:
        del model                              # DEVSIM holds one device set; the hold-out builds reset it
        ho_ref = holdout_curves(truth, vgs, holdout)
        print(f"hold-out reference curves ready ({clock():.0f} s)")
        for r in records:
            if "error" in r or r["value"] == v_true:
                continue
            try:
                th = np.array([r["theta"][n] for n in names])
                r["holdout"] = holdout_compare(ho_ref, holdout_curves(params_of(th), vgs, holdout), vgs)
                worst = max(r["holdout"].values(), key=lambda q: q["rms_dec"])
                print(f"[{prof} = {r['value']:.3f}] hold-out: worst RMS {worst['rms_dec']:.4f} dec, max |dVth| "
                      f"{max(abs(q['dVth_mV']) for q in r['holdout'].values()):.1f} mV ({clock():.0f} s)")
            except ConvergenceError as exc:
                r["holdout"] = {"error": str(exc)}
    records = [r for r in records if "error" not in r] + [r for r in records if "error" in r]
    ok = sorted([r for r in records if "error" not in r], key=lambda r: r["value"])
    dump_json({"param": prof, "names": names, "sigma_dec": args.sigma, "vg": vgs, "n_res": n_res, "vds": list(fit_vds),
               "records": ok, "failed": [r for r in records if "error" in r],
               "curves": {str(k): v for k, v in curves.items()}, "wall_s": clock()}, f"{out}/summary.json")
    vd_txt = "0.05/1.0" if fit_vds == (0.05, 1.0) else "/".join(f"{vd:g}" for vd in fit_vds)
    lines = [f"# E5 profile likelihood of {prof} (Step 5 twin, L = 1 um, VD = {vd_txt} V, σ = {args.sigma} dec)", "",
             f"{n_res} residuals; profile χ² < 3.84 means statistically indistinguishable from the truth (95 %).", "",
             f"| {prof} | χ² | RMS misfit (dec) | max misfit (dec) | " + " | ".join(free) + " |",
             "|---|---|---|---|" + "---|" * len(free)]
    for r in ok:
        lines.append(f"| {r['value']:.3f} | {r['chi2']:.3g} | {r['rms_dec']:.4f} | {r['max_dec']:.4f} | "
                     + " | ".join(f"{r['theta'][n]:.4g}" for n in free) + " |")
    if ho_ref is not None:
        keys = list(ho_ref)
        lines += ["", "Prediction of the re-fitted sets at the paper's hold-out conditions (vs the truth):", "",
                  f"| {prof} | " + " | ".join(f"{k}: RMS dec / ΔVth mV" for k in keys) + " |", "|---|" + "---|" * len(keys)]
        for r in ok:
            h = r.get("holdout")
            if not h or "error" in h:
                continue
            lines.append(f"| {r['value']:.3f} | " + " | ".join(f"{h[k]['rms_dec']:.4f} / {h[k]['dVth_mV']:+.1f}" for k in keys) + " |")
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
        ax[0].plot([r["value"] for r in ok], [r["chi2"] for r in ok], "o-")
        ax[0].axhline(3.84, color="r", ls="--", lw=0.8, label="95 % threshold")
        ax[0].axvline(v_true, color="k", lw=0.6)
        ax[0].set_yscale("symlog", linthresh=0.1)
        ax[0].set_xlabel(prof)
        ax[0].set_ylabel("profile χ²")
        ax[0].legend()
        for k, (v, cv) in enumerate(sorted(curves.items())):
            for vd, i in cv.items():
                ax[1].semilogy(vgs, np.maximum(np.array(i), 1e-16), "-" if vd == 1.0 else "--", lw=1.0,
                               label=f"{prof} = {v:.2f}" if vd == 1.0 else None)
        ax[1].axhline(1e-11, color="grey", lw=0.6, ls=":")
        ax[1].set_ylim(1e-13, 1e-2)
        ax[1].set_xlabel("VG (V)")
        ax[1].set_ylabel("ID (A)")
        ax[1].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(f"{out}/e5_profile.png", dpi=150)
    print(f"done in {clock():.0f} s -> {out}")


if __name__ == "__main__":
    main()
