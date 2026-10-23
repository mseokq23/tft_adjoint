"""E4 add-on - OED ranking of measurement conditions that are NOT in the paper's hold-out set.

Hypothesis tested: a fixed interface charge Qf and the gate work function shift
the transfer curve identically for one gate stack (both act as a VG offset), so
no VD or L variation can separate them; a device with a different gate-oxide
thickness can, because the Qf-induced shift scales with t_ox while the WF shift
does not.

Candidates (each evaluated at VD = 0.05 and 1.0 V on top of the calibration data):
  * t_gox = 10 nm (L = 1 um)
  * t_gox = 14 nm (L = 1 um)
  * L = 2 um (for reference; also in e3)

Usage:  python scripts/e4_oed_extra_conditions.py [--sigma 0.02] [--out results/e4]
"""
from __future__ import annotations

import argparse
from dataclasses import replace

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from e3_identifiability_oed_preview import NAMES, fim_stats, log_rows, sweep_with_grad
from tftsim import PAPER_SETS, Geometry
from tftsim.params import theta_scale


def full_stats(F):
    """Eigen-analysis of the full Fisher matrix without hiding structural null directions."""
    w, V = np.linalg.eigh(F)
    rel = w[0] / w[-1]
    out = {"lam_min": float(w[0]), "lam_max": float(w[-1]), "rel_lam_min": float(rel),
           "null_vec": dict(zip(NAMES, V[:, 0].round(4).tolist()))}
    if rel > 1e-12:                                   # identifiable: ordinary inverse is meaningful
        C = np.linalg.inv(F)
        sd = np.sqrt(np.diag(C))
        out["std_scaled"] = dict(zip(NAMES, sd.tolist()))
        iw, iq = NAMES.index("wf_gate"), NAMES.index("qf")
        out["corr_wfgate_qf"] = float(C[iw, iq] / (sd[iw] * sd[iq]))
    else:                                             # structurally singular: std is infinite
        out["std_scaled"] = None
        out["corr_wfgate_qf"] = None
    return out


def logdet_noqf(F):
    idx = [k for k, n in enumerate(NAMES) if n != "qf"]
    return float(np.linalg.slogdet(F[np.ix_(idx, idx)])[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, default=0.02)
    ap.add_argument("--step", type=float, default=0.05)
    ap.add_argument("--out", default="results/e4")
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    vgs = vg_grid(vmin=-1.0, step=args.step)
    scale = theta_scale(NAMES)
    p = PAPER_SETS["step5"]
    base, _, _ = sweep_with_grad(p, Geometry(), (0.05, 1.0), vgs)
    S_base = np.vstack([log_rows(*base[vd]) for vd in (0.05, 1.0)]) * scale
    F_base = S_base.T @ S_base / args.sigma ** 2
    st = full_stats(F_base)
    ld0 = logdet_noqf(F_base)
    print(f"[base] rel lam_min = {st['rel_lam_min']:.2e}  null vector = {st['null_vec']}  ({clock():.0f} s)")
    cands = {"tox10": Geometry(t_gox_nm=10.0), "tox14": Geometry(t_gox_nm=14.0), "L2": Geometry(L_um=2.0)}
    rows = {}
    for name, g in cands.items():
        d, _, _ = sweep_with_grad(p, g, (0.05, 1.0), vgs)
        rows[name] = np.vstack([log_rows(*d[vd]) for vd in (0.05, 1.0)]) * scale
        print(f"[{name}] sweep done ({clock():.0f} s)")
    np.savez(f"{out}/rows.npz", S_base=S_base, scale=scale, names=np.array(NAMES), **{f"S_{k}": v for k, v in rows.items()})
    combos = {k: [k] for k in rows}
    combos.update({"tox10+tox14": ["tox10", "tox14"], "tox10+L2": ["tox10", "L2"], "all": list(rows)})
    res = {}
    for name, members in combos.items():
        F = F_base + sum(rows[m].T @ rows[m] for m in members) / args.sigma ** 2
        r = full_stats(F)
        r["d_logdet_noqf"] = logdet_noqf(F) - ld0
        res[name] = r
    summary = {"sigma_dec": args.sigma, "base": st, "candidates": res, "wall_s": clock()}
    dump_json(summary, f"{out}/summary.json")
    lines = ["# E4 add-on: can a second gate-oxide thickness separate Qf from the gate work function?", "",
             f"- calibration data only (L = 1 um, t_ox = 7 nm): smallest/largest Fisher eigenvalue = {st['rel_lam_min']:.1e} "
             f"-> structurally non-identifiable; null direction (scaled units): "
             + ", ".join(f"{k} {v:+.2f}" for k, v in st["null_vec"].items() if abs(v) > 0.05), "",
             "| added data (VD = 0.05 + 1.0 V each) | λmin/λmax | σ(qf) [1e11 cm⁻²] | σ(wf_gate) [eV] | corr(wf_gate, qf) | Δlog det F (qf excluded) |",
             "|---|---|---|---|---|---|"]
    for name, r in res.items():
        if r["std_scaled"] is None:
            lines.append(f"| {name} | {r['rel_lam_min']:.1e} | ∞ | ∞ | – | {r['d_logdet_noqf']:.2f} |")
        else:
            sq = r["std_scaled"]["qf"] * scale[NAMES.index("qf")] / 1e11
            sw = r["std_scaled"]["wf_gate"] * scale[NAMES.index("wf_gate")]
            lines.append(f"| {name} | {r['rel_lam_min']:.1e} | {sq:.3g} | {sw:.3g} | {r['corr_wfgate_qf']:+.3f} | {r['d_logdet_noqf']:.2f} |")
    lines += ["", "1-sigma of every parameter (scaled units, typical change = 1) when the data identify all 10:", "",
              "| parameter | " + " | ".join(n for n in res if res[n]["std_scaled"] is not None) + " |",
              "|---|" + "---|" * sum(res[n]["std_scaled"] is not None for n in res)]
    for k in NAMES:
        lines.append(f"| {k} | " + " | ".join(f"{res[n]['std_scaled'][k]:.3g}" for n in res if res[n]["std_scaled"] is not None) + " |")
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"done in {clock():.0f} s -> {out}")


if __name__ == "__main__":
    main()
