"""E4 - spatial adjoint sensitivity map: where in the channel does a fixed charge change ID?

For a local donor-like (fixed positive) charge perturbation dN at channel node i,
the only residual entry that changes is the Poisson row of node i, so
    dID/dN_i = lam[row_i] * dR_i/dN_i ,  dR_i/dN_i = c * q * V_i ,
with lam the adjoint solution for the drain-current unknown and V_i the node
volume. The constant c (the assembly sign convention) is fixed by requiring
sum_i dID/dN_i = dID/dND_eff (global residual-difference gradient) - this is also
the self-consistency check printed by the script.

One transposed solve per bias point gives the full map (thousands of local
parameters), which is where the adjoint mode is much cheaper than forward/FD.

Usage:  python scripts/e4_sensitivity_map.py [--out results/e4map]
"""
from __future__ import annotations

import argparse
import math

import devsim as ds

from _common import HAVE_MPL, Timer, dump_json, np, outdir
from tftsim import PAPER_SETS, TFT
from tftsim.sensitivity import SensitivityEngine

POINTS = [(1.0, -0.60, "subthreshold"), (1.0, -0.20, "near-Vth"), (1.0, 2.00, "on-state")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e4map")
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    t = TFT(PAPER_SETS["step5"]).build()
    eng = SensitivityEngine(t)
    d = t.dev
    x = np.array(ds.get_node_model_values(device=d, region="ch", name="x")) * 1e4      # um
    y = np.array(ds.get_node_model_values(device=d, region="ch", name="y")) * 1e7      # nm
    vol = np.array(ds.get_node_model_values(device=d, region="ch", name="NodeVolume"))
    maps, checks = {}, []
    for vd, vg, label in POINTS:
        t.goto(vd=vd)
        t.goto(vg=vg)
        lu, R0, n = eng.factorize()
        e = np.zeros(n)
        e[eng.iI] = 1.0
        tm = Timer()
        lam = lu.solve(e, trans="T")
        t_adj = tm()
        rows = np.array(ds.get_equation_numbers(device=d, region="ch", equation="PotentialEquation"))
        # exact per-node residual derivative for a uniform fixed-charge (ND) change; it is q*V_i on
        # every channel Poisson row except the S/D contact nodes, whose rows are boundary conditions
        dR_nd = eng.dR("log10_nd_eff") / (t.p.nd_eff * math.log(10))
        m = lam[rows] * dR_nd[rows] * t.isign * t.width_cm     # dID per unit dN_i
        g_nd = eng.current_gradient(["log10_nd_eff"], "adjoint")[0] / (t.p.nd_eff * math.log(10))
        checks.append({"point": label, "sum_local_over_global": float(m.sum() / g_nd), "t_adjoint_s": t_adj})
        idv = t.drain_current()
        # relative change of ID per added elementary charge at the node: (dID/dN_i) / V_i / I
        s = m / vol / idv
        maps[label] = {"vd": vd, "vg": vg, "id": idv, "s_rel_per_charge": s}
        print(f"[{label}] ID={idv:.3e} A  sum(local)/global = {m.sum() / g_nd:.6f}  adjoint solve {t_adj * 1e3:.0f} ms")
    dump_json({"checks": checks, "x_um": x, "y_nm": y,
               "maps": {k: {kk: vv for kk, vv in v.items()} for k, v in maps.items()}, "wall_s": clock()},
              f"{out}/summary.json")
    if HAVE_MPL:
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
        tri = mtri.Triangulation(x, y)
        fig, ax = plt.subplots(len(maps), 1, figsize=(12, 2.6 * len(maps)), sharex=True)
        lov = t.geom.Lov_um
        for k, (label, mp) in enumerate(maps.items()):
            s = np.abs(mp["s_rel_per_charge"])
            s = np.log10(np.maximum(s / s.max(), 1e-6))
            pc = ax[k].tripcolor(tri, s, shading="gouraud", cmap="magma", vmin=-6, vmax=0)
            ax[k].axvline(lov, color="c", lw=0.8)
            ax[k].axvline(lov + t.geom.L_um, color="c", lw=0.8)
            ax[k].set_ylabel("y (nm)")
            ax[k].set_title(f"{label}: VG = {mp['vg']} V, VD = {mp['vd']} V, ID = {mp['id']:.2e} A  "
                            f"(log10 |∂ID/∂q(x,y)| normalised)")
            fig.colorbar(pc, ax=ax[k], fraction=0.02)
        ax[-1].set_xlabel("x (µm)   [cyan lines: source/drain contact edges]")
        fig.tight_layout()
        fig.savefig(f"{out}/e4_sensitivity_map.png", dpi=150)
        # 1-D profile: integrate over channel thickness
        fig, axp = plt.subplots(figsize=(8, 3.6))
        xb = np.unique(np.round(x, 5))
        for label, mp in maps.items():
            prof = np.array([np.sum(np.abs(mp["s_rel_per_charge"][np.isclose(x, xv)] * vol[np.isclose(x, xv)])) for xv in xb])
            axp.semilogy(xb, prof / prof.max(), label=label)
        axp.axvline(lov, color="grey", lw=0.6)
        axp.axvline(lov + t.geom.L_um, color="grey", lw=0.6)
        axp.set_xlabel("x (µm)")
        axp.set_ylabel("|∂ID/∂Q(x)| / max (thickness-integrated)")
        axp.legend()
        fig.tight_layout()
        fig.savefig(f"{out}/e4_sensitivity_profile.png", dpi=150)
    print(f"done in {clock():.0f} s -> {out}")


if __name__ == "__main__":
    main()
