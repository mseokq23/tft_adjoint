"""E2 - synthetic-twin calibration: gradient (adjoint) LM vs finite-difference LM (vs DFO-LS).

Ground truth = final 2%-W model of arXiv:2609.12184 (Step 5). "Measured" curves
are simulated at L = 1 um, VD = 0.05/1.0 V and perturbed by Gaussian noise of
sigma decades in log10(ID). Each method starts from the same generic oxide-TFT
guess and minimises the log-current residuals (ID >= 10 pA) inside bounds.

Reported per method: final cost, PDE solves, wall time, parameter errors
relative to the truth (in typical-change units) and the Gauss-Newton 1-sigma at
the solution (is the error explained by the predicted uncertainty?).

Usage:
  python scripts/e2_twin_calibration.py --methods adjoint,fd [--names all|noqf]
         [--max_nfev 25] [--step 0.1] [--sigma 0.02] [--out results/e2]
"""
from __future__ import annotations

import argparse
from dataclasses import replace

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from scipy.optimize import least_squares

from tftsim import PAPER_SETS, THETA_ALL
from tftsim.calib import CurveModel, log_residual_problem
from tftsim.params import theta_scale, theta_vector, with_theta

START = dict(wf_sd=4.30, wf_gate=4.20, qf=0.0, nd_eff=1.0e18, mu0=10.0, n_gd=1.0e19,
             e_gd=0.20, s_gd=0.10, n_ta=3.0e19, w_ta=0.10)
BOUNDS = {"wf_sd": (3.90, 4.80), "wf_gate": (3.70, 4.80), "qf": (-2e12, 2e12),
          "log10_nd_eff": (15.0, 19.5), "log10_mu0": (0.0, 2.0), "log10_n_gd": (17.0, 21.0),
          "e_gd": (0.0, 0.5), "s_gd": (0.02, 0.30), "log10_n_ta": (18.0, 21.0), "w_ta": (0.03, 0.30)}


def run_method(method, start_params, names, theta0, scale, meas, sigma, vgs, max_nfev, lb, ub):
    model = CurveModel(start_params, names, vds=tuple(meas), vgs=vgs)
    n0 = model.t.n_solves
    hist = []
    fun, jac, masks = log_residual_problem(model, meas, sigma, theta0, scale, history=hist)
    clock = Timer()
    if method == "adjoint":
        res = least_squares(fun, np.zeros(len(names)), jac=jac, bounds=(lb, ub), method="trf",
                            x_scale=1.0, max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        x, Jz = res.x, res.jac
    elif method == "fd":
        res = least_squares(fun, np.zeros(len(names)), jac="2-point", diff_step=1e-3, bounds=(lb, ub),
                            method="trf", x_scale=1.0, max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        x, Jz = res.x, res.jac
    elif method == "dfols":
        import dfols  # pip install DFO-LS  (GPL-3.0-or-later; optional)
        res = dfols.solve(fun, np.zeros(len(names)), bounds=(lb, ub), maxfun=max_nfev * (len(names) + 1),
                          rhobeg=0.5, rhoend=1e-4)
        x, Jz = res.x, None
    else:
        raise ValueError(method)
    wall = clock()
    r = fun(x)
    if Jz is None:
        Jz = jac(x)
    JTJ = Jz.T @ Jz
    cov = np.linalg.pinv(JTJ, rcond=1e-12)
    std = np.sqrt(np.maximum(np.diag(cov), 0.0))
    return {"method": method, "z": x, "cost": 0.5 * float(r @ r), "n_res": int(r.size),
            "chi2_per_dof": float(r @ r) / max(1, r.size - len(names)),
            "pde_solves": model.t.n_solves - n0, "wall_s": wall, "history": hist,
            "std_z": std, "model_stats": model.stats, "engine": model.eng.report(),
            "nfev": getattr(res, "nfev", None), "njev": getattr(res, "njev", None)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--methods", default="adjoint,fd")
    ap.add_argument("--names", default="all", choices=["all", "noqf"])
    ap.add_argument("--max_nfev", type=int, default=25)
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--sigma", type=float, default=0.02)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="results/e2")
    args = ap.parse_args()
    out = outdir(args.out)
    names = THETA_ALL if args.names == "all" else [n for n in THETA_ALL if n != "qf"]
    scale = theta_scale(names)
    vgs = vg_grid(-1.0, 2.5, args.step)
    truth = PAPER_SETS["step5"]
    start = replace(truth, **START)
    if args.names == "noqf":
        start = replace(start, qf=truth.qf)
    theta_true = theta_vector(truth, names)
    theta0 = theta_vector(start, names)
    lb = np.array([(BOUNDS[n][0] - theta0[k]) / scale[k] for k, n in enumerate(names)])
    ub = np.array([(BOUNDS[n][1] - theta0[k]) / scale[k] for k, n in enumerate(names)])

    # ----- synthetic measurement
    clock = Timer()
    twin = CurveModel(truth, names, vds=(0.05, 1.0), vgs=vgs)
    I_true, _ = twin.evaluate(theta_true)
    rng = np.random.default_rng(args.seed)
    meas = {vd: I_true[vd] * 10 ** (args.sigma * rng.standard_normal(I_true[vd].size)) for vd in I_true}
    print(f"twin data ready ({clock():.0f} s)")

    results = []
    for m in args.methods.split(","):
        clock = Timer()
        try:
            r = run_method(m, start, names, theta0, scale, meas, args.sigma, vgs, args.max_nfev, lb, ub)
        except ImportError as exc:
            print(f"[{m}] skipped: {exc}")
            continue
        th = theta0 + r["z"] * scale
        err_z = (th - theta_true) / scale
        r.update({"theta_hat": dict(zip(names, th)), "err_typical_units": dict(zip(names, err_z)),
                  "z_score": dict(zip(names, err_z / np.where(r["std_z"] > 0, r["std_z"], np.inf)))})
        results.append(r)
        print(f"[{m:8s}] cost={r['cost']:.1f} chi2/dof={r['chi2_per_dof']:.2f} PDE solves={r['pde_solves']} "
              f"wall={r['wall_s']:.0f}s  max|err|={np.max(np.abs(err_z)):.2f} typical units ({clock():.0f} s)")

    summary = {"names": names, "truth": dict(zip(names, theta_true)), "start": dict(zip(names, theta0)),
               "scale": dict(zip(names, scale)), "sigma_dec": args.sigma, "vg_step": args.step,
               "results": results}
    dump_json(summary, f"{out}/summary.json")
    lines = ["# E2 synthetic-twin calibration", "",
             f"- truth = Step 5 (2%-W), start = generic guess {START}", f"- noise sigma = {args.sigma} dec, VG step {args.step} V",
             "", "| method | cost | χ²/dof | PDE solves | wall (s) | max abs error (typical units) |", "|---|---|---|---|---|---|"]
    for r in results:
        lines.append(f"| {r['method']} | {r['cost']:.1f} | {r['chi2_per_dof']:.2f} | {r['pde_solves']} | {r['wall_s']:.0f} | "
                     f"{max(abs(v) for v in r['err_typical_units'].values()):.2f} |")
    lines += ["", "| parameter | truth | start | " + " | ".join(f"{r['method']} estimate ± 1σ" for r in results) + " |",
              "|---|---|---|" + "---|" * len(results)]
    for k, n in enumerate(names):
        cells = []
        for r in results:
            cells.append(f"{r['theta_hat'][n]:.4g} ± {r['std_z'][k] * scale[k]:.2g}")
        lines.append(f"| {n} | {theta_true[k]:.4g} | {theta0[k]:.4g} | " + " | ".join(cells) + " |")
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    if HAVE_MPL and results:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
        for r in results:
            h = r["history"]
            best, xs, ys = np.inf, [], []
            for e in h:
                best = min(best, e["cost"])
                xs.append(e["n_solves"])
                ys.append(best)
            ax[0].semilogy(xs, ys, "o-", ms=3, label=r["method"])
            ax[1].bar(np.arange(len(names)) + 0.35 * results.index(r), [abs(r["err_typical_units"][n]) for n in names],
                      width=0.35, label=r["method"])
        ax[0].set_xlabel("cumulative nonlinear PDE solves")
        ax[0].set_ylabel("best cost so far")
        ax[0].legend()
        ax[1].set_xticks(np.arange(len(names)))
        ax[1].set_xticklabels(names, rotation=60, fontsize=8)
        ax[1].set_ylabel("|estimate - truth| (typical units)")
        ax[1].set_yscale("log")
        ax[1].legend()
        fig.tight_layout()
        fig.savefig(f"{out}/e2_twin.png", dpi=150)
    print(open(f"{out}/summary.md").read())


if __name__ == "__main__":
    main()
