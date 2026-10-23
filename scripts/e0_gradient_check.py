"""E0 - gradient verification on the BG-IWO TFT (final 2%-W model of arXiv:2609.12184).

For several bias points (subthreshold, near Vth, on-state, linear) it compares
dID/dtheta from
  (a) adjoint sensitivity   (one transposed solve for all parameters)
  (b) forward sensitivity   (one solve per parameter)
  (c) central finite differences with two full nonlinear re-solves (two steps)
and runs
  * a Taylor-remainder test along a random direction (should give slope 2),
  * an FD step-size scan (truncation vs. solver-noise trade-off of FD).

Usage:  python scripts/e0_gradient_check.py [--set step5] [--out results/e0]
"""
from __future__ import annotations

import argparse
import csv
import math

from _common import HAVE_MPL, Timer, dump_json, np, outdir
from tftsim import PAPER_SETS, THETA_ALL, TFT
from tftsim.params import theta_scale
from tftsim.sensitivity import SensitivityEngine

POINTS = [  # (VD, VG, label)
    (1.0, -0.60, "subthreshold"),
    (1.0, -0.20, "near-Vth"),
    (1.0, 2.00, "on-state"),
    (0.05, 0.50, "linear"),
]
FD_SMALL = {"lin": 1e-4, "qf": 1e8, "log10": 1e-4}
FD_LARGE = {"lin": 1e-3, "qf": 1e9, "log10": 1e-3}


def fd_step(name, table):
    if name == "qf":
        return table["qf"]
    if name.startswith("log10"):
        return table["log10"]
    return table["lin"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="step5")
    ap.add_argument("--out", default="results/e0")
    args = ap.parse_args()
    out = outdir(args.out)
    names = THETA_ALL + ["gate_bias"]

    clock = Timer()
    t = TFT(PAPER_SETS[args.set]).build()
    eng = SensitivityEngine(t)
    print(f"built {args.set}: {t.info()}  ({clock():.1f} s)")

    rows, per_point = [], []
    for vd, vg, label in POINTS:
        t.goto(vd=vd)
        t.goto(vg=vg)
        idv = t.drain_current()
        tm = Timer()
        ga = eng.current_gradient(names, "adjoint")
        t_adj = tm()
        tm = Timer()
        gf = eng.current_gradient(names, "forward")
        t_fwd = tm()
        fd_s, fd_l, t_fd = [], [], 0.0
        for nm in names:
            tm = Timer()
            fd_s.append(eng.resolve_fd(nm, h=fd_step(nm, FD_SMALL)))
            t_fd += tm()
            fd_l.append(eng.resolve_fd(nm, h=fd_step(nm, FD_LARGE)))
        fd_s, fd_l = np.array(fd_s), np.array(fd_l)
        for k, nm in enumerate(names):
            ref = fd_s[k]
            rows.append({
                "point": label, "vd": vd, "vg": vg, "id_A": idv, "param": nm,
                "adjoint": ga[k], "forward": gf[k], "fd_small": fd_s[k], "fd_large": fd_l[k],
                "rel_adj_fwd": abs(ga[k] - gf[k]) / max(abs(gf[k]), 1e-300),
                "rel_adj_fdsmall": abs(ga[k] - ref) / max(abs(ref), 1e-300),
                "rel_adj_fdlarge": abs(ga[k] - fd_l[k]) / max(abs(fd_l[k]), 1e-300),
                "rel_fdsmall_fdlarge": abs(fd_s[k] - fd_l[k]) / max(abs(fd_s[k]), 1e-300),
            })
        # normwise error: what matters for an optimiser / Fisher matrix at this point
        norm_err = float(np.max(np.abs(ga - fd_s)) / np.max(np.abs(fd_s)))
        per_point.append({"point": label, "vd": vd, "vg": vg, "id_A": idv,
                          "t_adjoint_s": t_adj, "t_forward_s": t_fwd, "t_fd_resolve_s": t_fd,
                          "n_params": len(names), "normwise_err_adj_vs_fd": norm_err})
        print(f"[{label:12s}] VD={vd} VG={vg:+.2f} ID={idv:.3e} A | adjoint {t_adj:.2f}s  "
              f"forward {t_fwd:.2f}s  FD re-solve {t_fd:.2f}s")

    # ---------------------------------------------------- Taylor remainder test
    vd, vg, _ = POINTS[1]
    t.goto(vd=vd)
    t.goto(vg=vg)
    snap = t.snapshot()
    th0 = np.array([t.get_theta(n) for n in THETA_ALL])
    sc = theta_scale(THETA_ALL)
    rng = np.random.default_rng(0)
    d = rng.standard_normal(len(THETA_ALL))
    d /= np.linalg.norm(d)
    g = eng.current_gradient(THETA_ALL, "adjoint")
    i0 = t.drain_current()
    taylor = []
    for eps in (0.3, 0.1, 0.03, 0.01, 0.003, 0.001):
        for k, n in enumerate(THETA_ALL):
            t.set_theta(n, th0[k] + eps * d[k] * sc[k])
        t.restore(snap)
        t.solve()
        i1 = t.drain_current()
        lin = eps * float(np.dot(g, d * sc))
        taylor.append({"eps": eps, "r0": abs(i1 - i0), "r1": abs(i1 - i0 - lin)})
        for k, n in enumerate(THETA_ALL):
            t.set_theta(n, th0[k])
        t.restore(snap)
    e = np.log10([r["eps"] for r in taylor])
    slope0 = float(np.polyfit(e[1:], np.log10([r["r0"] for r in taylor][1:]), 1)[0])
    slope1 = float(np.polyfit(e[1:4], np.log10([r["r1"] for r in taylor][1:4]), 1)[0])
    print(f"Taylor test: slope(r0)={slope0:.2f} (expect 1), slope(r1)={slope1:.2f} (expect 2)")

    # -------------------------------------------------------- FD step-size scan
    scan = []
    ref_g = eng.current_gradient(["w_ta", "log10_n_ta", "e_gd"], "adjoint")
    for k, nm in enumerate(["w_ta", "log10_n_ta", "e_gd"]):
        for h in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
            try:
                v = eng.resolve_fd(nm, h=h)
                scan.append({"param": nm, "h": h, "fd": v, "adjoint": ref_g[k],
                             "rel_err": abs(v - ref_g[k]) / abs(ref_g[k])})
            except Exception as exc:  # convergence failure at large steps
                scan.append({"param": nm, "h": h, "fd": None, "adjoint": ref_g[k], "rel_err": None,
                             "error": str(exc)[:80]})

    # ------------------------------------------------------------------ save
    with open(f"{out}/gradients.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    summary = {
        "set": args.set, "device": t.info(), "points": per_point, "taylor": taylor,
        "taylor_slopes": {"r0": slope0, "r1": slope1}, "fd_scan": scan,
        "max_rel_adj_fwd": max(r["rel_adj_fwd"] for r in rows),
        "median_rel_adj_fdsmall": float(np.median([r["rel_adj_fdsmall"] for r in rows])),
        "max_rel_adj_fdsmall": max(r["rel_adj_fdsmall"] for r in rows),
        "engine": eng.report(), "wall_s": clock(),
    }
    dump_json(summary, f"{out}/summary.json")
    lines = ["# E0 gradient check", "",
             f"- device: {t.info()}",
             f"- max |adjoint-forward|/|forward| = {summary['max_rel_adj_fwd']:.2e}",
             f"- adjoint vs FD(re-solve, small step): median {summary['median_rel_adj_fdsmall']:.2e}, "
             f"max {summary['max_rel_adj_fdsmall']:.2e}",
             f"- Taylor remainder slopes: r0 {slope0:.2f} (1 expected), r1 {slope1:.2f} (2 expected)", "",
             "| point | ID (A) | adjoint (s) | forward (s) | FD re-solve (s) | normwise err |",
             "|---|---|---|---|---|---|"]
    lines += [f"| {p['point']} | {p['id_A']:.3e} | {p['t_adjoint_s']:.2f} | {p['t_forward_s']:.2f} | "
              f"{p['t_fd_resolve_s']:.2f} | {p['normwise_err_adj_vs_fd']:.1e} |" for p in per_point]
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")

    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
        a = np.array([abs(r["adjoint"]) for r in rows])
        b = np.array([abs(r["fd_small"]) for r in rows])
        m = (a > 0) & (b > 0)
        ax[0].loglog(b[m], a[m], "o", ms=4)
        lim = [min(a[m].min(), b[m].min()), max(a[m].max(), b[m].max())]
        ax[0].plot(lim, lim, "k--", lw=0.8)
        ax[0].set_xlabel("|dID/dθ| FD re-solve")
        ax[0].set_ylabel("|dID/dθ| adjoint")
        ax[0].set_title("Adjoint vs re-solve FD")
        ax[1].loglog([r["eps"] for r in taylor], [r["r0"] for r in taylor], "o-", label="no gradient (slope 1)")
        ax[1].loglog([r["eps"] for r in taylor], [r["r1"] for r in taylor], "s-", label="with adjoint gradient (slope 2)")
        ax[1].set_xlabel("step ε (scaled)")
        ax[1].set_ylabel("Taylor remainder |ΔID| (A)")
        ax[1].legend()
        ax[1].set_title("Taylor remainder test")
        for nm in ("w_ta", "log10_n_ta", "e_gd"):
            s = [r for r in scan if r["param"] == nm and r["rel_err"] is not None]
            ax[2].loglog([r["h"] for r in s], [max(r["rel_err"], 1e-16) for r in s], "o-", label=nm)
        ax[2].set_xlabel("FD step h")
        ax[2].set_ylabel("|FD - adjoint| / |adjoint|")
        ax[2].set_title("FD needs step tuning; adjoint does not")
        ax[2].legend()
        fig.tight_layout()
        fig.savefig(f"{out}/e0_gradient_check.png", dpi=150)
    print(f"done in {clock():.1f} s -> {out}")


if __name__ == "__main__":
    main()
