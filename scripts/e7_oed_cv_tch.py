"""E7 - extended OED: quasi-static C-V and channel-thickness variants as extra measurements.

Question: which ONE additional measurement would best constrain the descriptors
that the paper's calibration data leave sloppy (E3/E4/E5/E6)?

Candidates (all on top of the calibration data L = 1 um, VD = 0.05/1.0 V):
  * paper hold-out I-V: VD = 0.1/0.2/0.5 V (L = 1 um), L = 2 um          (rows from E3)
  * second gate-oxide thickness: t_ox = 10 / 14 nm                        (rows from E4)
  * second channel thickness: t_IWO = 5 / 7 nm (I-V, VD = 0.05/1.0 V)     (new)
  * quasi-static gate C-V of the calibration device (VD = 0, S/D shorted) (new)

C-V data model: C_{j+1/2} = (Q_G(V_{j+1}) - Q_G(V_j)) / dV on the VG grid, noise
sigma_C = frac * C_ox A (absolute), frac = 0.5 % by default (also 0.2/1/2 %).
dQ_G/dtheta comes from ONE adjoint solve per bias (tftsim.cv.GateCharge), so the
C-V sensitivities are exactly consistent with the discretised model.
Self-checks: the gate-charge functional against DEVSIM's contact charge, and its
adjoint gradients against re-solve finite differences.

Usage:  python scripts/e7_oed_cv_tch.py [--out results/e7] [--frac 0.005]
"""
from __future__ import annotations

import argparse

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from e3_identifiability_oed_preview import NAMES, log_rows, sweep_with_grad
from tftsim import PAPER_SETS, TFT, Geometry
from tftsim.cv import GateCharge
from tftsim.params import THETA_DEFS, theta_scale
from tftsim.sensitivity import SensitivityEngine

KEY = ["wf_gate", "qf", "log10_nd_eff", "log10_n_ta", "w_ta", "log10_n_gd", "e_gd", "s_gd", "wf_sd", "log10_mu0"]


def cv_sweep(params, geom, vgs, checks: bool = True):
    t = TFT(params, geom).build()
    gq = GateCharge(t)
    eng = SensitivityEngine(t, central=True)
    t.goto(vd=0.0)
    t.goto(vg=vgs[0])
    Q, G, C = [], [], []
    chk = []
    for vg in vgs:
        t.goto(vg=vg)
        lu, R0, n = eng.factorize()
        g, _ = eng.functional_gradient(gq.c_vector(n), list(NAMES) + ["gate_bias"], lu=lu, R0=R0)
        Q.append(gq.value())
        G.append(g[:-1])
        C.append(g[-1])
        if checks and any(abs(vg - v) < 1e-9 for v in (-0.5, 0.5)):
            q_dev = gq.devsim_value()
            row = {"vg": float(vg), "Q_functional": Q[-1], "Q_devsim": q_dev,
                   "rel_diff": abs(Q[-1] - q_dev) / abs(q_dev), "grad": {}}
            for nm in ("gate_bias", "wf_gate", "log10_nd_eff", "log10_n_ta"):
                fd = resolve_q(t, gq, nm, 1e-3)
                a = C[-1] if nm == "gate_bias" else g[NAMES.index(nm)]
                row["grad"][nm] = {"adjoint": float(a), "resolve_fd": float(fd), "rel": float(abs(a - fd) / abs(fd))}
            chk.append(row)
    return np.array(Q), np.array(G), np.array(C), gq.c_ox_area(), chk


def resolve_q(t, gq, name, h):
    v0 = t.get_theta(name)
    snap = t.snapshot()
    try:
        t.set_theta(name, v0 + h)
        t.solve()
        qp = gq.value()
        t.restore(snap)
        t.set_theta(name, v0 - h)
        t.solve()
        qm = gq.value()
    finally:
        t.set_theta(name, v0)
        t.restore(snap)
    return (qp - qm) / (2.0 * h)


def stats(F, scale):
    """Full-Fisher eigen-analysis + reduced (no-qf) log det + physical 1-sigma (inf if singular)."""
    w = np.linalg.eigvalsh(F)
    idx = [k for k, n in enumerate(NAMES) if n != "qf"]
    ld = float(np.linalg.slogdet(F[np.ix_(idx, idx)])[1])
    rel = float(w[0] / w[-1])
    if rel > 1e-12:
        sd = np.sqrt(np.diag(np.linalg.inv(F))) * scale
    else:                                      # Qf-WF degenerate: report the reduced model, qf = inf
        sdr = np.sqrt(np.diag(np.linalg.inv(F[np.ix_(idx, idx)]))) * scale[idx]
        sd = np.full(len(NAMES), np.inf)
        sd[idx] = sdr
    return {"rel_lam_min": rel, "logdet_noqf": ld, "std_phys": dict(zip(NAMES, sd.tolist()))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e7")
    ap.add_argument("--sigma", type=float, default=0.02, help="I-V noise, decades")
    ap.add_argument("--frac", type=float, default=0.005, help="C-V noise as a fraction of C_ox*A")
    ap.add_argument("--step", type=float, default=0.05)
    ap.add_argument("--reuse", action="store_true", help="reuse the saved C-V / t_ch rows (no new simulations)")
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    scale = theta_scale(NAMES)
    p = PAPER_SETS["step5"]
    vgs = vg_grid(vmin=-1.0, step=args.step)

    if args.reuse:                                   # recompute the tables/figure from saved sensitivities
        zs = np.load(f"{out}/rows.npz", allow_pickle=True)
        prev = __import__("json").load(open(f"{out}/summary.json"))
        S_cv_unit, coxa, chk = zs["S_cv_unit"], float(zs["coxa"]), prev["checks"]
        Q, Cqs = np.array(prev["cv_curve"]["Q"]), np.array(prev["cv_curve"]["C_qs"])
        rows = {k[2:]: zs[k] for k in zs.files if k.startswith("S_tch")}
    else:
        Q, Gq, Cqs, coxa, chk, S_cv_unit, rows = simulate(p, vgs, scale, clock)
    run_oed(args, out, clock, scale, vgs, Q, Cqs, coxa, chk, S_cv_unit, rows)


def simulate(p, vgs, scale, clock):
    # ---------------- C-V of the calibration device
    Q, Gq, Cqs, coxa, chk = cv_sweep(p, Geometry(), vgs)
    print(f"C-V sweep done ({clock():.0f} s); C_ox*A = {coxa:.3e} F")
    for c in chk:
        print(f"  check VG={c['vg']:+.2f}: Q functional vs DEVSIM rel {c['rel_diff']:.1e}; "
              + "; ".join(f"{k} rel {v['rel']:.1e}" for k, v in c["grad"].items()))
    dv = np.diff(vgs)
    S_cv_unit = (np.diff(Gq, axis=0) / dv[:, None]) * scale          # dC/dtheta per typical change (F)

    # ---------------- channel-thickness variants (I-V)
    rows = {}
    for tch in (5.0, 7.0):
        d, _, _ = sweep_with_grad(p, Geometry(t_iwo_nm=tch), (0.05, 1.0), vgs)
        rows[f"tch{tch:g}"] = np.vstack([log_rows(*d[vd]) for vd in (0.05, 1.0)]) * scale
        print(f"t_IWO = {tch} nm I-V sweeps done ({clock():.0f} s)")
    return Q, Gq, Cqs, coxa, chk, S_cv_unit, rows


def run_oed(args, out, clock, scale, vgs, Q, Cqs, coxa, chk, S_cv_unit, rows):
    rows = dict(rows)
    # ---------------- previously computed rows (same NAMES order, same scaling, sigma applied below)
    z3 = np.load("results/e3/rows.npz", allow_pickle=True)
    z4 = np.load("results/e4/rows.npz", allow_pickle=True)
    assert [str(n) for n in z3["names"]] == NAMES and [str(n) for n in z4["names"]] == NAMES
    S_base = z3["S_base"]
    for k in z3.files:
        if k.startswith("S_") and k != "S_base":
            rows[k[2:]] = z3[k]
    for k in ("tox10", "tox14"):
        rows[k] = z4[f"S_{k}"]
    rows["L2 (both VD)"] = z4["S_L2"]                # same 2-curve format as the t_ox / t_ch candidates
    F0 = S_base.T @ S_base / args.sigma ** 2
    st0 = stats(F0, scale)

    def fisher_cv(frac):
        s = frac * coxa
        return S_cv_unit.T @ S_cv_unit / s ** 2

    cands = {k: v.T @ v / args.sigma ** 2 for k, v in rows.items()}
    cands["CV"] = fisher_cv(args.frac)
    cats = {"L1_VD0.1": "paper hold-out", "L1_VD0.2": "paper hold-out", "L1_VD0.5": "paper hold-out",
            "L2_VD0.05": "paper hold-out", "L2_VD1.0": "paper hold-out", "tox10": "new: t_ox", "tox14": "new: t_ox",
            "tch5": "new: t_ch", "tch7": "new: t_ch", "CV": "new: C-V", "L2 (both VD)": "paper hold-out"}
    combos = {"tox14+CV": ["tox14", "CV"], "tox14+tch7": ["tox14", "tch7"], "tox14+tch7+CV": ["tox14", "tch7", "CV"],
              "all paper hold-outs": ["L1_VD0.1", "L1_VD0.2", "L1_VD0.5", "L2_VD0.05", "L2_VD1.0"]}
    res = {}
    for name, F in cands.items():
        res[name] = dict(stats(F0 + F, scale), category=cats.get(name, ""))
    for name, members in combos.items():
        res[name] = dict(stats(F0 + sum(cands[m] for m in members), scale), category="combination")
    cv_noise = {f: stats(F0 + fisher_cv(f), scale)["logdet_noqf"] - st0["logdet_noqf"] for f in (0.002, 0.005, 0.01, 0.02)}

    summary = {"sigma_dec": args.sigma, "cv_frac": args.frac, "cox_area_F": coxa, "checks": chk, "base": st0,
               "candidates": res, "cv_noise_dlogdet": cv_noise,
               "cv_curve": {"vg": vgs, "Q": Q, "C_qs": Cqs}, "wall_s": clock()}
    dump_json(summary, f"{out}/summary.json")
    np.savez(f"{out}/rows.npz", S_cv_unit=S_cv_unit, coxa=coxa, names=np.array(NAMES), scale=scale,
             **{f"S_{k}": v for k, v in rows.items() if k.startswith("tch")})

    fmt = lambda v: "∞" if not np.isfinite(v) else f"{v:.3g}"
    unit = {n: ("eV" if THETA_DEFS[n][0] == "lin" and n != "qf" else ("cm⁻²" if n == "qf" else "dec")) for n in NAMES}
    lines = ["# E7: which extra measurement constrains the sloppy descriptors best?", "",
             f"- base = calibration I-V (L = 1 µm, VD = 0.05/1.0 V), σ = {args.sigma} dec; C-V noise = {100 * args.frac:g} % of C_ox·A "
             f"({coxa:.3g} F)",
             f"- Δlog det F is computed without Qf (exact Qf–WF degeneracy unless t_ox varies); 'breaks Qf–WF' = λmin/λmax > 1e-12",
             "- where Qf–WF is not broken, σ wf_gate is that of the effective work function WF − qQf/C_ox (Qf merged)",
             "- tox/tch/'L2 (both VD)' add two I-V curves (VD = 0.05 and 1.0 V); L1_VD*/L2_VD* add one curve; CV adds one C-V curve", "",
             "| added measurement | type | Δlog det F | breaks Qf–WF | " + " | ".join(f"σ {n} ({unit[n]})" for n in KEY) + " |",
             "|---|---|---|---|" + "---|" * len(KEY),
             "| (calibration only) | – | 0 | no | " + " | ".join(fmt(st0["std_phys"][n]) for n in KEY) + " |"]
    for name, r in sorted(res.items(), key=lambda kv: -kv[1]["logdet_noqf"]):
        lines.append(f"| {name} | {r['category']} | {r['logdet_noqf'] - st0['logdet_noqf']:.2f} | "
                     f"{'yes' if r['rel_lam_min'] > 1e-12 else 'no'} | " + " | ".join(fmt(r["std_phys"][n]) for n in KEY) + " |")
    lines += ["", "C-V information gain vs assumed C-V noise: " + ", ".join(f"{100 * f:g} % → Δlogdet {v:.2f}" for f, v in cv_noise.items())]
    lines += ["", "Self-checks (gate-charge functional and adjoint gradients):"]
    for c in chk:
        lines.append(f"- VG = {c['vg']:+.2f} V: |Q_functional − Q_DEVSIM|/|Q| = {c['rel_diff']:.1e}; "
                     + ", ".join(f"{k}: adjoint vs re-solve FD rel {v['rel']:.1e}" for k, v in c["grad"].items()))
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(14, 4.6))
        ax[0].plot(vgs, Cqs / coxa, "k-")
        ax[0].set_xlabel("VG (V)")
        ax[0].set_ylabel("C_G / (C_ox A)")
        ax[0].set_title("Quasi-static gate C-V (VD = 0, Step-5 model)")
        singles = ["L2 (both VD)", "tox10", "tox14", "tch5", "tch7", "CV", "all paper hold-outs"]
        singles.sort(key=lambda k: res[k]["logdet_noqf"])
        colors = {"paper hold-out": "C7", "new: t_ox": "C0", "new: t_ch": "C2", "new: C-V": "C3", "combination": "0.4"}
        ax[1].barh(singles, [res[k]["logdet_noqf"] - st0["logdet_noqf"] for k in singles],
                   color=[colors[res[k]["category"]] for k in singles])
        ax[1].set_xlabel("Δ log det F (information gain, Qf excluded)")
        ax[1].set_title("Extra data (2 I-V curves, or 1 C-V curve); grey = paper's hold-out conditions")
        fig.tight_layout()
        fig.savefig(f"{out}/e7_oed.png", dpi=150)
    print(f"done in {clock():.0f} s -> {out}")


if __name__ == "__main__":
    main()
