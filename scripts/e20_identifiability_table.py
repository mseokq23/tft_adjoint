"""E20 - what do the data identify? Fisher 1σ of every descriptor, with and without the recommended tests.

Linear analysis at two points:
  * measured: the joint free fit to the digitised 2 %-W curves (E18: L = 1 µm, VD = 0.05/0.1/0.2/0.5/1 V,
    grid points at measured VG only). The calibration rows are exactly the residual rows of that fit
    (same grid, mask and σ); their χ² is checked against the stored value. Noise for all I–V data: σ_eff of
    that fit (model error included); for the recommended I–V tests this assumes the same model error.
  * twin: the paper's Step 5 with σ = 0.02 dec (rows of E3/E4/E7, no new simulation).
Recommended tests (as in E4/E7/E13): t_ox = 14 nm and t_IWO = 7 nm devices at VD = 0.05/1 V (I–V, VG step
0.1 V like the calibration fit grid, ID >= 10 pA as in E4/E7, also with the 100 pA floor of the measured
calibration; the thick channel from VG = -3 V) and quasi-static gate C–V of the calibration device (VD = 0,
noise 0.5 % of C_ox·A, also 2 %). These are local linear values at one point of a flat valley.
Qf is fixed in every column (it is exactly degenerate with the gate WF unless t_ox changes); the "all three"
combination is also given with Qf free.

Usage:
  python scripts/e20_identifiability_table.py [--out results/e20]      # about 15 min (sensitivities at the E18 optimum)
  python scripts/e20_identifiability_table.py --report                 # tables from the saved rows, no simulation
"""
from __future__ import annotations

import argparse
import json
import os

from _common import Timer, dump_json, np, outdir, vg_grid
from e3_identifiability_oed_preview import log_rows, sweep_with_grad
from e7_oed_cv_tch import cv_sweep
from e11_calibrate_measured import load_meas, params_from
from tftsim import PAPER_SETS, THETA_ALL, Geometry
from tftsim.calib import CurveModel, log_residual_problem
from tftsim.metrics import I_FLOOR
from tftsim.params import THETA_DEFS, theta_scale, theta_vector

NAMES = THETA_ALL
FIT = "results/e18_joint_meas/free/summary.json"
UNIT = {n: ("dec" if THETA_DEFS[n][0] == "log10" else "cm⁻²" if n == "qf" else "eV") for n in NAMES}
ORDER = ["wf_gate", "wf_sd", "log10_mu0", "log10_nd_eff", "log10_n_ta", "w_ta", "log10_n_gd", "e_gd", "s_gd"]
MEAS_FLOOR = 1e-10      # current floor of the measured calibration rows (e11 --floor)


def measured_mask(fit, vgs, meas_path="data/meas.csv"):
    """Measured curves on the fit grid with exactly the rules of e11_calibrate_measured.py (floor, max_gap)."""
    meas_all = load_meas(meas_path, fit["set"])
    floor, max_gap = MEAS_FLOOR, fit.get("max_gap")
    meas = {}
    for vd in fit["vds"]:
        m_vg, m_i = meas_all[vd]
        below = np.where(m_i < floor)[0]
        if below.size:
            keep = m_vg > m_vg[below[-1]] if m_vg[below[-1]] < 0.5 else np.ones_like(m_vg, bool)
            m_vg, m_i = m_vg[keep], m_i[keep]
        meas[vd] = np.interp(vgs, m_vg, m_i, left=np.nan, right=np.nan)
        if max_gap is not None:
            j = np.searchsorted(m_vg, vgs - 1e-6)
            at_point = (j < m_vg.size) & (np.abs(m_vg[np.minimum(j, m_vg.size - 1)] - vgs) < 1e-6)
            inside = (j > 0) & (j < m_vg.size) & (m_vg[np.minimum(j, m_vg.size - 1)] - m_vg[np.maximum(j - 1, 0)]
                                                  <= max_gap + 1e-9)
            meas[vd] = np.where(at_point | inside, meas[vd], np.nan)
        meas[vd] = np.where(np.isfinite(meas[vd]) & (meas[vd] >= floor), meas[vd], 1e-30)
    return meas


def simulate(out, clock):
    fit = json.load(open(FIT))
    base = PAPER_SETS["step5"]
    th = theta_vector(base, NAMES)
    for k, n in enumerate(NAMES):
        if n in fit["theta_hat"]:
            th[k] = fit["theta_hat"][n]
    p = params_from(base, NAMES, th)
    scale = theta_scale(NAMES)
    sigma = fit["sigma_dec"]

    # calibration rows = residual rows of the E18 fit (same grid and mask), in scaled units, σ removed
    vgs_fit = vg_grid(-1.0, 2.0, 0.1)
    meas = measured_mask(fit, vgs_fit)
    model = CurveModel(p, NAMES, vds=tuple(fit["vds"]), vgs=vgs_fit)
    fun, jac, masks = log_residual_problem(model, meas, sigma, th, scale)
    z0 = np.zeros(len(NAMES))
    r = fun(z0)
    J = jac(z0) * sigma
    n_vd = [int(masks[vd].sum()) for vd in fit["vds"]]
    print(f"calibration rows {J.shape}, χ² {float(r @ r):.1f} (stored {fit['chi2']:.1f}) ({clock():.0f} s)", flush=True)
    del model

    rows, cur = {}, {}
    for tag, geom, vmin in (("tox14", Geometry(t_gox_nm=14.0), -1.5), ("tiwo7", Geometry(t_iwo_nm=7.0), -3.0)):
        d, _, _ = sweep_with_grad(p, geom, (0.05, 1.0), vg_grid(vmin, 2.0, 0.1))
        for vd in (0.05, 1.0):
            rows[f"{tag}_VD{vd:g}"] = log_rows(*d[vd]) * scale
            keep = d[vd][1] >= I_FLOOR                # the rows kept by log_rows, in the same order
            cur[f"I_{tag}_VD{vd:g}"], cur[f"VG_{tag}_VD{vd:g}"] = d[vd][1][keep], d[vd][0][keep]
            print(f"{tag} VD {vd:g}: {rows[f'{tag}_VD{vd:g}'].shape[0]} rows, ID at VG_min {d[vd][1][0]:.1e} A ({clock():.0f} s)",
                  flush=True)
    vgs_cv = vg_grid(vmin=-1.0, step=0.05)
    Q, Gq, Cqs, coxa, chk = cv_sweep(p, Geometry(), vgs_cv)
    S_cv_unit = (np.diff(Gq, axis=0) / np.diff(vgs_cv)[:, None]) * scale
    print(f"C–V done ({clock():.0f} s); C_ox·A = {coxa:.3e} F", flush=True)
    np.savez(os.path.join(out, "rows.npz"), names=np.array(NAMES), scale=scale, S_cal=J, n_vd=np.array(n_vd),
             vds=np.array(fit["vds"]), S_cv_unit=S_cv_unit, coxa=coxa, cv_vg=vgs_cv, cv_Q=Q, cv_C=Cqs,
             **{f"S_{k}": v for k, v in rows.items()}, **cur)
    dump_json({"fit": FIT, "theta": dict(zip(NAMES, th.tolist())), "chi2_check": [float(r @ r), fit["chi2"]],
               "cv_checks": chk, "wall_s": clock()}, os.path.join(out, "sim_meta.json"))


def sd_phys(F, names, scale, free_qf=False):
    idx = [k for k, n in enumerate(names) if free_qf or n != "qf"]
    C = np.linalg.inv(F[np.ix_(idx, idx)])
    sd = np.full(len(names), np.nan)
    sd[idx] = np.sqrt(np.diag(C)) * scale[idx]
    return dict(zip(names, sd.tolist()))


def report(out):
    z = np.load(os.path.join(out, "rows.npz"))
    meta = json.load(open(os.path.join(out, "sim_meta.json")))
    fit = json.load(open(FIT))
    names = [str(n) for n in z["names"]]
    scale = z["scale"]
    seff = fit["sigma_eff"]
    S_cal, n_vd, vds = z["S_cal"], list(z["n_vd"]), list(z["vds"])
    cuts = np.cumsum([0] + n_vd)
    pair = np.vstack([S_cal[cuts[i]:cuts[i + 1]] for i, vd in enumerate(vds) if vd in (0.05, 1.0)])
    F = lambda S, s: S.T @ S / s ** 2
    tests = {"t_ox 14 nm": F(np.vstack([z["S_tox14_VD0.05"], z["S_tox14_VD1"]]), seff),
             "t_IWO 7 nm": F(np.vstack([z["S_tiwo7_VD0.05"], z["S_tiwo7_VD1"]]), seff),
             "C–V 0.5 %": F(z["S_cv_unit"], 0.005 * float(z["coxa"])),
             "C–V 2 %": F(z["S_cv_unit"], 0.02 * float(z["coxa"]))}
    Fj = F(S_cal, seff)
    cols = {"pair (VD 0.05/1 V)": F(pair, seff), "joint (5 VD, E18)": Fj,
            "joint + t_ox 14 nm": Fj + tests["t_ox 14 nm"], "joint + t_IWO 7 nm": Fj + tests["t_IWO 7 nm"],
            "joint + C–V (0.5 %)": Fj + tests["C–V 0.5 %"], "joint + C–V (2 %)": Fj + tests["C–V 2 %"],
            "joint + all three": Fj + tests["t_ox 14 nm"] + tests["t_IWO 7 nm"] + tests["C–V 0.5 %"],
            "joint + all three (C–V 2 %)": Fj + tests["t_ox 14 nm"] + tests["t_IWO 7 nm"] + tests["C–V 2 %"]}
    meas_sd = {c: sd_phys(Fc, names, scale) for c, Fc in cols.items()}
    meas_sd["joint + all three, Qf free"] = sd_phys(cols["joint + all three"], names, scale, free_qf=True)
    # the same with the current floor of the measured calibration (100 pA) for the recommended I–V tests
    floor_sd = {}
    if all(f"I_{k}" in z.files for k in ("tox14_VD0.05", "tox14_VD1", "tiwo7_VD0.05", "tiwo7_VD1")):
        def rows_above(tag, floor):
            return np.vstack([z[f"S_{tag}_VD{vd}"][z[f"I_{tag}_VD{vd}"] >= floor] for vd in ("0.05", "1")])
        n_rows = {f"{tag} ({lab})": int(rows_above(tag, fl).shape[0]) for tag in ("tox14", "tiwo7")
                  for lab, fl in (("10 pA", 0.0), ("100 pA", MEAS_FLOOR))}
        t100 = {"t_ox 14 nm": F(rows_above("tox14", MEAS_FLOOR), seff), "t_IWO 7 nm": F(rows_above("tiwo7", MEAS_FLOOR), seff)}
        fcols = {"+ t_ox 14 nm, 10 pA": cols["joint + t_ox 14 nm"], "+ t_ox 14 nm, 100 pA": Fj + t100["t_ox 14 nm"],
                 "+ t_IWO 7 nm, 10 pA": cols["joint + t_IWO 7 nm"], "+ t_IWO 7 nm, 100 pA": Fj + t100["t_IWO 7 nm"],
                 "+ all three, 10 pA": cols["joint + all three"],
                 "+ all three, 100 pA": Fj + t100["t_ox 14 nm"] + t100["t_IWO 7 nm"] + tests["C–V 0.5 %"]}
        floor_sd = {c: sd_phys(Fc, names, scale) for c, Fc in fcols.items()}
    # twin at Step 5 from the stored rows of E3/E4/E7 (σ = 0.02 dec)
    z3, z4, z7 = (np.load(f"results/{e}/rows.npz") for e in ("e3", "e4", "e7"))
    assert [str(n) for n in z3["names"]] == names
    Ft0 = F(z3["S_base"], 0.02)
    Ft_all = Ft0 + F(z4["S_tox14"], 0.02) + F(z7["S_tch7"], 0.02) + F(z7["S_cv_unit"], 0.005 * float(z7["coxa"]))
    twin_sd = {"twin: pair": sd_phys(Ft0, names, scale), "twin: + all three": sd_phys(Ft_all, names, scale),
               "twin: + all three, Qf free": sd_phys(Ft_all, names, scale, free_qf=True)}
    stored = fit.get("std_eff", {})

    def fmt(v, n):
        if not np.isfinite(v):
            return "–"
        if n == "qf":
            return f"{v:.2g}"
        return f"{v:.3g}" if v < 100 else f"{v:.0f}"

    lines = ["# E20: what do the data identify? Fisher 1σ of every descriptor", "",
             f"Measured columns: linear analysis at the joint free fit to the digitised 2 %-W curves (E18, log10 ND,eff = "
             f"{fit['theta_hat']['log10_nd_eff']:.2f}), σ_eff = {seff:.4f} dec for all I–V data (model error included; for the "
             f"recommended I–V tests this assumes the same model error), C–V noise as a fraction of C_ox·A. Calibration rows "
             f"= the {S_cal.shape[0]} residual rows of that fit (χ² {meta['chi2_check'][0]:.1f}, stored {meta['chi2_check'][1]:.1f}). "
             "Twin columns: Step 5, σ = 0.02 dec, rows of E3/E4/E7. Qf fixed unless stated. "
             "Recommended I–V tests: VD = 0.05/1 V, ID ≥ 10 pA as in E4/E7 (t_ox 14 nm from VG = -1.5 V, t_IWO 7 nm from -3 V, "
             "0.1 V step; the measured calibration rows use a 100 pA floor, see the floor table below); "
             "C–V: VD = 0, VG = -1…2.5 V in 0.05 V steps. The columns with recommended tests add them to the joint (5-curve) "
             "calibration; the pair column is the VD 0.05/1 V subset of the joint rows at the same point.", "",
             "| descriptor | unit | " + " | ".join(list(meas_sd) + list(twin_sd)) + " |",
             "|---|---|" + "---|" * (len(meas_sd) + len(twin_sd))]
    for n in ORDER + ["qf"]:
        lines.append(f"| {n} | {UNIT[n]} | " + " | ".join(fmt(d[n], n) for d in list(meas_sd.values()) + list(twin_sd.values())) + " |")
    chk = {n: [meas_sd["joint (5 VD, E18)"][n], stored[n]] for n in stored}
    worst = max(abs(a / b - 1) for a, b in chk.values())
    lines += ["", f"Check: the joint-calibration column reproduces the 1σ stored with the E18 fit (largest relative "
              f"difference {worst:.1e}; the stored values use the Jacobian of the last accepted step)."]
    res = {"sigma_eff": seff, "measured": meas_sd, "twin": twin_sd, "stored_std_eff_check": chk,
           "n_rows_calibration": int(S_cal.shape[0]), "n_rows_vd": dict(zip([f"{v:g}" for v in vds], n_vd)),
           "chi2_check": meta["chi2_check"]}
    if floor_sd:
        lines += ["", "## Current floor of the recommended I–V tests", "",
                  "Measured point, joint calibration plus the test(s); C–V 0.5 % in the last two columns. Rows kept (both VD): "
                  + ", ".join(f"{k} {v}" for k, v in n_rows.items()) + ".", "",
                  "| descriptor | unit | " + " | ".join(floor_sd) + " |", "|---|---|" + "---|" * len(floor_sd)]
        for n in ORDER:
            lines.append(f"| {n} | {UNIT[n]} | " + " | ".join(fmt(d[n], n) for d in floor_sd.values()) + " |")
        res["floor_check"] = {"floor_A": MEAS_FLOOR, "n_rows": n_rows, "sd": floor_sd}
    dump_json(res, os.path.join(out, "summary.json"))
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e20")
    ap.add_argument("--report", action="store_true", help="only rebuild the tables from the saved rows")
    args = ap.parse_args()
    out = outdir(args.out)
    if not args.report:
        simulate(out, Timer())
    report(out)


if __name__ == "__main__":
    main()
