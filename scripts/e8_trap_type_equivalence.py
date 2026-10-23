"""E8 - donor- vs acceptor-like character of a trap level is not identifiable when ND,eff is calibrated.

In quasi-equilibrium a donor-like level of density N carries +N(1 - f) and an
acceptor-like level -N f. Since  N(1 - f) = N - N f,  a donor-like DOS is
EXACTLY equivalent to an acceptor-like DOS of the same shape plus a uniform
fixed positive charge N_tot = integral of the DOS, i.e. ND,eff -> ND,eff + N_tot
(both uniform in the channel in this descriptor class).

Numerical proof: the Step-5 model (donor-like Gaussian) against a model whose
Gaussian is acceptor-like and whose ND,eff is raised by the Gaussian's total
density. Transfer curves at the calibration and hold-out conditions, and the
quasi-static C-V, must agree to solver tolerance.

With --cases thickness the same comparison is made for devices outside the paper's split (thicker
channel t_IWO = 5 / 7 nm, thicker gate oxide t_ox = 14 nm, and their C-V): the equivalence is a
property of the model equations, so no electrical measurement in this descriptor class can break it.

Usage:  python scripts/e8_trap_type_equivalence.py [--out results/e8]
        python scripts/e8_trap_type_equivalence.py --cases thickness --out results/e8_thickness
"""
from __future__ import annotations

import argparse
from dataclasses import replace

from _common import Timer, dump_json, np, outdir, vg_grid
import tftsim.device as tdev
from tftsim import PAPER_SETS, TFT, Geometry
from tftsim.cv import GateCharge
from tftsim.params import Material, default_dos_edges
from tftsim.params import trap_levels as trap_levels_orig

_FLIP = {"on": False}


def trap_levels_flipped(p, mat, edges):
    """Same DOS, but the Gaussian levels are acceptor-like (added to the tail)."""
    c, wA, wD, n1 = trap_levels_orig(p, mat, edges)
    if _FLIP["on"]:
        return c, wA + wD, np.zeros_like(wD), n1
    return c, wA, wD, n1


def curves(params, flip, geom, vds, vgs):
    _FLIP["on"] = flip
    t = TFT(params, geom).build()
    out = {}
    for j, vd in enumerate(vds):
        seq = vgs if j % 2 == 0 else vgs[::-1]
        recs = sorted(t.transfer(vd, seq), key=lambda r: r["vg"])
        out[vd] = np.array([r["id"] for r in recs])
    return out


CASES = {
    "paper": (("L1", Geometry(), (0.05, 0.1, 0.2, 0.5, 1.0)), ("L2", Geometry(L_um=2.0), (0.05, 1.0))),
    "thickness": (("tIWO5", Geometry(t_iwo_nm=5.0), (0.05, 1.0)), ("tIWO7", Geometry(t_iwo_nm=7.0), (0.05, 1.0)),
                  ("tox14", Geometry(t_gox_nm=14.0), (0.05, 1.0))),
}
CV_GEOMS = {"paper": {"L1": Geometry()},
            "thickness": {"tIWO7": Geometry(t_iwo_nm=7.0), "tox14": Geometry(t_gox_nm=14.0)}}


def cv_curve(params, flip, vgs, geom=None):
    _FLIP["on"] = flip
    t = TFT(params, geom or Geometry()).build()
    gq = GateCharge(t)
    t.goto(vd=0.0)
    q = []
    for vg in vgs:
        t.goto(vg=vg)
        q.append(gq.value())
    return np.array(q)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e8")
    ap.add_argument("--cases", default="paper", choices=list(CASES))
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    tdev.trap_levels = trap_levels_flipped              # the device model calls this name
    p = PAPER_SETS["step5"]
    _, _, wD, _ = trap_levels_orig(p, Material(), default_dos_edges())
    n_tot = float(np.sum(wD))
    p_flip = replace(p, nd_eff=p.nd_eff + n_tot)
    vgs = vg_grid(-1.0, 2.5, 0.05)
    rep = {"gauss_total_cm3": n_tot, "nd_eff_original": p.nd_eff, "nd_eff_flipped": p_flip.nd_eff, "cases": {}}
    rep["case_set"] = args.cases
    for name, geom, vds in CASES[args.cases]:
        a = curves(p, False, geom, vds, vgs)
        b = curves(p_flip, True, geom, vds, vgs)
        for vd in vds:
            m = a[vd] >= 1e-11
            dev = np.abs(np.log10(b[vd][m]) - np.log10(a[vd][m]))
            rep["cases"][f"{name}_VD{vd:g}"] = float(dev.max())
        print(f"[{name}] max |Δlog10 ID| = {max(rep['cases'][f'{name}_VD{vd:g}'] for vd in vds):.2e} dec ({clock():.0f} s)")
    rep["cv"] = {}
    for name, geom in CV_GEOMS[args.cases].items():
        qa, qb = cv_curve(p, False, vgs, geom), cv_curve(p_flip, True, vgs, geom)
        ca, cb = np.diff(qa), np.diff(qb)
        rep["cv"][name] = float(np.max(np.abs(cb - ca)) / np.max(np.abs(ca)))
        print(f"[C-V {name}] max |ΔC|/max C = {rep['cv'][name]:.2e} ({clock():.0f} s)")
    rep["cv_max_rel"] = max(rep["cv"].values())
    dump_json(rep, f"{out}/summary.json")
    lines = [f"# E8: donor-like Gaussian ≡ acceptor-like Gaussian + uniform fixed charge ({args.cases} conditions)", "",
             f"- Gaussian total density ∫N_GD dE = {n_tot:.3e} cm⁻³ → ND,eff {p.nd_eff:.3e} → {p_flip.nd_eff:.3e} cm⁻³",
             "- max |Δlog10 ID| (ID ≥ 10 pA): " + ", ".join(f"{k} {v:.1e}" for k, v in rep["cases"].items()),
             "- quasi-static C–V: max |ΔC| / max C = " + ", ".join(f"{k} {v:.1e}" for k, v in rep["cv"].items())]
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
