"""E4b - energy-resolved trap sensitivity: which part of the DOS does a transfer curve probe?

For every discretised trap level k the residual derivative is local and analytic:
    acceptor level:  dR_i/dNtA_k = -b_i * f_k(n_i)
    donor level:     dR_i/dNtD_k = +b_i * (1 - f_k(n_i)),   f_k = n/(n + n1_k)
where b_i = dR_i/dND (= q*V_i on channel Poisson rows, 0 on S/D contact rows),
so with one adjoint solve per bias point
    dID/dN_k = sum_i lam[row_i] * dR_i/dN_k   for all 2K levels at once.
Self-check: ln(10) * sum_k NtA_k dID/dNtA_k must equal dID/dlog10(n_ta) (and
the same for the Gaussian donors).

Output: heat maps of |N_k dID/dN_k| / ID (relative current change per relative
change of the level density) versus VG and trap energy, at VD = 0.05 and 1.0 V.

Usage:  python scripts/e4_energy_sensitivity.py [--out results/e4energy]
"""
from __future__ import annotations

import argparse
import math

import devsim as ds

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from tftsim import PAPER_SETS, TFT
from tftsim.params import trap_levels
from tftsim.sensitivity import SensitivityEngine


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e4energy")
    ap.add_argument("--step", type=float, default=0.1)
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    t = TFT(PAPER_SETS["step5"]).build()
    eng = SensitivityEngine(t)
    d = t.dev
    rows = np.array(ds.get_equation_numbers(device=d, region="ch", equation="PotentialEquation"))
    dE, wA, wD, n1 = trap_levels(t.p, t.mat, t.edges)
    b_full = eng.dR("log10_nd_eff") / (t.p.nd_eff * math.log(10))
    b = b_full[rows]                               # q*V_i (0 on contact rows)
    scale_i = t.isign * t.width_cm
    vgs = vg_grid(-1.0, 2.5, args.step)
    res = {}
    checks = []
    for j, vd in enumerate((0.05, 1.0)):
        t.goto(vd=vd)
        seq = vgs if j % 2 == 0 else vgs[::-1]
        SA, SD, I = {}, {}, {}
        t.goto(vg=seq[0])
        for vg in seq:
            t.goto(vg=vg)
            lu, R0, n = eng.factorize()
            e = np.zeros(n)
            e[eng.iI] = 1.0
            lam = lu.solve(e, trans="T")[rows]
            ne = np.array(ds.get_node_model_values(device=d, region="ch", name="Electrons"))
            f = ne[:, None] / (ne[:, None] + n1[None, :])            # nodes x levels
            gA = scale_i * ((lam * b)[:, None] * (-f)).sum(axis=0)     # dID/dNtA_k
            gD = scale_i * ((lam * b)[:, None] * (1.0 - f)).sum(axis=0)
            idv = t.drain_current()
            SA[vg], SD[vg], I[vg] = gA, gD, idv
            if abs(vg - 0.0) < 1e-9 or abs(vg - 1.0) < 1e-9:
                g = eng.current_gradient(["log10_n_ta", "log10_n_gd"], "adjoint")
                checks.append({"vd": vd, "vg": vg,
                               "tail_ratio": float(math.log(10) * np.dot(wA, gA) / g[0]),
                               "gauss_ratio": float(math.log(10) * np.dot(wD, gD) / g[1])})
        vg_s = np.array(sorted(I))
        width = np.diff(t.edges)                    # eV; per-eV density removes the bin-width artefact
        res[vd] = {"vg": vg_s, "id": np.array([I[v] for v in vg_s]),
                   "relA": np.array([np.abs(wA * SA[v]) / I[v] / width for v in vg_s]),
                   "relD": np.array([np.abs(wD * SD[v]) / I[v] / width for v in vg_s])}
        print(f"VD={vd}: done ({clock():.0f} s)")
    for c in checks:
        print(f"self-check VD={c['vd']} VG={c['vg']}: tail {c['tail_ratio']:.6f}, gauss {c['gauss_ratio']:.6f} (expect 1)")
    dump_json({"dE_eV": dE, "checks": checks,
               "res": {str(k): {kk: vv for kk, vv in v.items()} for k, v in res.items()}, "wall_s": clock()},
              f"{out}/summary.json")
    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 2, figsize=(13, 7.5), sharey=True)
        for c, vd in enumerate((0.05, 1.0)):
            r = res[vd]
            for rr, (key, ttl) in enumerate((("relA", "acceptor-like tail"), ("relD", "donor-like Gaussian"))):
                z = np.log10(np.maximum(r[key].T, 1e-8))
                pc = ax[rr, c].pcolormesh(r["vg"], dE, z, shading="nearest", cmap="viridis", vmin=-5, vmax=1)
                ax[rr, c].set_title(f"{ttl}, VD = {vd} V")
                ax[rr, c].set_xlabel("VG (V)")
                fig.colorbar(pc, ax=ax[rr, c], label="log10 |N(E) δID/δN(E)| / ID  (per eV)")
            ax[0, c].set_ylim(1.0, 0.0)
            ax[1, c].set_ylim(1.0, 0.0)
        ax[0, 0].set_ylabel("Ec − E (eV)")
        ax[1, 0].set_ylabel("Ec − E (eV)")
        fig.tight_layout()
        fig.savefig(f"{out}/e4_energy_sensitivity.png", dpi=150)
    print(f"done in {clock():.0f} s -> {out}")


if __name__ == "__main__":
    main()
