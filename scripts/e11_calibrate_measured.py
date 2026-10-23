"""E11 - calibrate the open-source model to the digitised measurements (adjoint Levenberg-Marquardt).

Data: data/meas.csv (2 %-W or 4 %-W, L = 1 µm, VD = 0.05 and 1 V, ID ≥ 10 pA), as in the
paper's nominal calibration split. Residuals r = (log10 ID_sim − log10 ID_meas)/σ with σ = 0.02 dec
(the digitisation scatter measured in digitize_paper.py). Qf is held fixed (exactly degenerate with
the gate work function); any other parameter can be fixed with --fix (used for profile likelihoods).

After the fit: χ², RMS misfit, the paper's metrics and J, a Gauss-Newton 1σ with σ inflated so that
χ²/dof = 1 (model error included), and - for 2 %-W - predictions at the paper's hold-out conditions
compared with the measured hold-out curves.

Usage:
  python scripts/e11_calibrate_measured.py --set 2W --out results/e11_2W
  python scripts/e11_calibrate_measured.py --set 2W --fix log10_nd_eff=17.4 --start results/e11_2W/summary.json --out ...
  python scripts/e11_calibrate_measured.py --set 4W --vgmin -0.7 --start w4 --out results/e11_4W
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os

from _common import ROOT, Timer, dump_json, np, outdir, vg_grid
from scipy.optimize import least_squares

from e1_forward_reproduce import sweep_device
from e2_twin_calibration import BOUNDS
from tftsim import PAPER_SETS, THETA_ALL, ConvergenceError, Geometry
from tftsim.calib import CurveModel, log_residual_problem
from tftsim.metrics import I_FLOOR, extract, objective_J
from tftsim.params import Material, theta_scale, theta_vector, with_theta


def load_meas(path, set_, L=1.0):
    d = collections.defaultdict(list)
    for r in csv.DictReader(open(path)):
        if r["set"] == set_ and float(r["L_um"]) == L:
            d[float(r["vd"])].append((float(r["vg"]), float(r["id"])))
    return {vd: (np.array([p[0] for p in sorted(v)]), np.array([p[1] for p in sorted(v)])) for vd, v in d.items()}


def params_from(base, names, th):
    p = base
    for n, v in zip(names, th):
        p = with_theta(p, n, float(v))
    return p


def compare(sim_vg, sim_i, m_vg, m_i):
    s = np.interp(m_vg, sim_vg, sim_i)
    ok = m_i >= I_FLOOR
    e = np.log10(np.maximum(s[ok], 1e-30)) - np.log10(m_i[ok])
    ms, mm = extract(m_vg, s), extract(m_vg, m_i)
    return {"E_log_dec": float(np.sqrt(np.mean(e ** 2))), "dVth_mV": 1e3 * (ms["vth"] - mm["vth"]),
            "dSS_mVdec": ms["ss_mVdec"] - mm["ss_mVdec"], "dlogIon_dec": ms["log10_ion"] - mm["log10_ion"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meas", default="data/meas.csv")
    ap.add_argument("--set", default="2W", choices=["2W", "4W"])
    ap.add_argument("--start", default=None, help="step5 | w4 | path to an earlier summary.json")
    ap.add_argument("--fix", action="append", default=[], help="name=value (theta units), repeatable")
    ap.add_argument("--vgmin", type=float, default=-1.0)
    ap.add_argument("--step", type=float, default=0.1, help="VG step of the fitted points")
    ap.add_argument("--sigma", type=float, default=0.02)
    ap.add_argument("--floor", type=float, default=1e-10, help="measured points below this are noise-floor dominated and not fitted")
    ap.add_argument("--max_nfev", type=int, default=25)
    ap.add_argument("--no_holdout", action="store_true")
    ap.add_argument("--nc", type=float, default=None, help="material: effective conduction-band DOS (cm^-3), default 5e18")
    ap.add_argument("--eps_ch", type=float, default=None, help="material: channel relative permittivity, default 10")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    mat = Material(**{k: v for k, v in (("Nc_cm3", args.nc), ("eps_ch", args.eps_ch)) if v is not None})

    base = PAPER_SETS["w4" if args.set == "4W" else "step5"]
    start_theta = {}
    if args.start in ("step5", "w4"):
        base = PAPER_SETS[args.start]
    elif args.start:
        prev = json.load(open(args.start))
        start_theta = prev["theta_hat"]
    fixed = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in args.fix}
    free = [n for n in THETA_ALL if n != "qf" and n not in fixed]
    names = free + list(fixed)
    scale = theta_scale(names)
    th0 = theta_vector(base, names)
    for k, n in enumerate(names):
        if n in start_theta:
            th0[k] = start_theta[n]
        if n in fixed:
            th0[k] = fixed[n]
    base = params_from(base, names, th0)

    meas_all = load_meas(args.meas, args.set)
    vgs = vg_grid(args.vgmin, 2.0, args.step)
    meas = {}
    for vd in (0.05, 1.0):
        m_vg, m_i = meas_all[vd]
        # drop the noise-floor region: everything up to the last point below the floor (isolated
        # floor spikes above 1e-10 A at the start of the sweep are measurement noise, not the device)
        below = np.where(m_i < args.floor)[0]
        if below.size:
            keep = m_vg > m_vg[below[-1]] if m_vg[below[-1]] < 0.5 else np.ones_like(m_vg, bool)
            m_vg, m_i = m_vg[keep], m_i[keep]
        meas[vd] = np.interp(vgs, m_vg, m_i, left=np.nan, right=np.nan)          # measured on the fit grid
        meas[vd] = np.where(np.isfinite(meas[vd]) & (meas[vd] >= args.floor), meas[vd], 1e-30)
    model = CurveModel(base, names, vds=(0.05, 1.0), vgs=vgs, mat=mat)
    fun_all, jac_all, masks = log_residual_problem(model, meas, args.sigma, th0, scale)
    nfix = len(fixed)
    lb = np.array([(BOUNDS[n][0] - th0[k]) / scale[k] for k, n in enumerate(free)])
    ub = np.array([(BOUNDS[n][1] - th0[k]) / scale[k] for k, n in enumerate(free)])
    zf = np.zeros(nfix)
    fun = lambda z: fun_all(np.concatenate([z, zf]))
    jac = lambda z: jac_all(np.concatenate([z, zf]))[:, :len(free)]
    r0 = fun(np.zeros(len(free)))
    print(f"start χ² = {float(r0 @ r0):.1f} on {r0.size} residuals ({clock():.0f} s)")
    res = least_squares(fun, np.zeros(len(free)), jac=jac, bounds=(lb, ub), method="trf", x_scale=1.0,
                        max_nfev=args.max_nfev, ftol=1e-8, xtol=1e-8, gtol=1e-8)
    th = th0.copy()
    th[:len(free)] = th0[:len(free)] + res.x * scale[:len(free)]
    r = fun(res.x)
    chi2 = float(r @ r)
    dof = max(1, r.size - len(free))
    Jz = res.jac
    sig_eff = args.sigma * np.sqrt(max(1.0, chi2 / dof))
    cov = np.linalg.pinv(Jz.T @ Jz, rcond=1e-12) * (sig_eff / args.sigma) ** 2
    std = np.sqrt(np.maximum(np.diag(cov), 0)) * scale[:len(free)]
    print(f"fit: χ² = {chi2:.1f} (dof {dof}), RMS = {args.sigma * np.sqrt(chi2 / r.size):.4f} dec, nfev {res.nfev} "
          f"({clock():.0f} s)")
    pfit = params_from(base, names, th)

    # ---- full-resolution evaluation, metrics, J
    del model
    vfull = vg_grid(args.vgmin, 2.0, 0.05)
    _, cv = sweep_device(pfit, Geometry(), vds=(0.05, 1.0), vgs=vfull, mat=mat)
    mc = {}
    for vd in (0.05, 1.0):
        m_vg, m_i = meas_all[vd]
        sel = m_vg >= args.vgmin - 1e-9
        mc[vd] = (m_vg[sel], m_i[sel])
    J = objective_J(cv, mc)
    cmp_cal = {str(vd): compare(*cv[vd], *mc[vd]) for vd in (0.05, 1.0)}
    ho = {}
    if args.set == "2W" and not args.no_holdout:
        for geom, vds, tag in ((Geometry(), (0.1, 0.2, 0.5), "L1"), (Geometry(L_um=2.0), (0.05, 1.0), "L2")):
            try:
                _, hc = sweep_device(pfit, geom, vds=vds, vgs=vfull, mat=mat)
            except ConvergenceError as exc:
                ho[tag] = {"error": str(exc)}
                continue
            mh = load_meas(args.meas, "2W", geom.L_um)
            for vd in vds:
                if vd in mh:
                    ho[f"{tag}_VD{vd:g}"] = compare(*hc[vd], *mh[vd])
                    ho[f"{tag}_VD{vd:g}"]["curve"] = hc[vd][1]
    summary = {"set": args.set, "material": {"Nc_cm3": mat.Nc_cm3, "eps_ch": mat.eps_ch},
               "names": names, "free": free, "fixed": fixed, "sigma_dec": args.sigma,
               "theta_hat": dict(zip(names, th.tolist())), "std_eff": dict(zip(free, std.tolist())), "sigma_eff": float(sig_eff),
               "chi2": chi2, "n_res": int(r.size), "dof": dof, "rms_dec": float(args.sigma * np.sqrt(chi2 / r.size)),
               "J": J["J"], "terms": J["terms"], "cal": cmp_cal, "holdout": ho, "nfev": int(res.nfev), "status": int(res.status),
               "vg": vfull, "curves": {str(vd): v[1] for vd, v in cv.items()}, "wall_s": clock()}
    dump_json(summary, f"{out}/summary.json")
    lines = [f"# E11 calibration to measured {args.set} curves (L = 1 µm, VD = 0.05/1 V, VG ≥ {args.vgmin} V)", "",
             f"- fixed: {fixed if fixed else 'Qf only (degenerate with WF)'}; start: "
             f"{(os.path.relpath(args.start, ROOT) if args.start and os.path.exists(args.start) else args.start) or 'paper set'}; "
             f"fitted points: ID ≥ {args.floor:g} A, VG step {args.step} V",
             f"- material: Nc = {mat.Nc_cm3:.3g} cm⁻³, ε_ch = {mat.eps_ch:g}" + (" (defaults)" if mat == Material() else " (changed)"),
             f"- χ² = {chi2:.1f} on {r.size} residuals (σ = {args.sigma} dec) → RMS misfit {summary['rms_dec']:.4f} dec; σ_eff = {sig_eff:.4f} dec",
             f"- J (paper weights, our a_i) = {J['J']:.3f}", "", "| parameter | value | 1σ (σ_eff) |", "|---|---|---|"]
    for k, n in enumerate(names):
        lines.append(f"| {n} | {th[k]:.4g} | {std[k]:.3g} |" if n in free else f"| {n} | {th[k]:.4g} (fixed) | – |")
    lines += ["", "| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |", "|---|---|---|---|---|"]
    for k, c in list(cmp_cal.items()) + [(k, c) for k, c in ho.items() if "error" not in c]:
        lines.append(f"| {k} | {c['E_log_dec']:.3f} | {c['dVth_mV']:+.1f} | {c['dSS_mVdec']:+.1f} | {c['dlogIon_dec']:+.3f} |")
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"done in {clock():.0f} s -> {out}")


if __name__ == "__main__":
    main()
