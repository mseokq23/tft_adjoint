"""E3b - post-processing of E3 on the identifiable subspace (no new simulations).

E3 shows that Qf and the gate work function are exactly degenerate for one gate
stack (Fisher eigenvalue ~1e-17 of the largest). Any E-optimal criterion or
determinant that keeps that direction is dominated by round-off, so here the
OED ranking is recomputed with Qf removed (equivalently: Qf merged into the gate
work function, WF_eff = WF_gate - q Qf / C_ox).

Also reports the eigen-directions of the reduced Fisher matrix and the 1-sigma
of every parameter after adding each candidate condition.

Usage:  python scripts/e3b_identifiable_subspace.py [--e3 results/e3] [--sigma 0.02]
"""
from __future__ import annotations

import argparse
import itertools

from _common import dump_json, np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--e3", default="results/e3")
    ap.add_argument("--sigma", type=float, default=0.02)
    args = ap.parse_args()
    z = np.load(f"{args.e3}/rows.npz", allow_pickle=True)
    names = [str(n) for n in z["names"]]
    scale = z["scale"]
    keep = [k for k, n in enumerate(names) if n != "qf"]
    kn = [names[k] for k in keep]
    F0 = (z["S_base"].T @ z["S_base"])[np.ix_(keep, keep)] / args.sigma ** 2
    cands = {k[2:]: z[k][:, keep] for k in z.files if k.startswith("S_") and k != "S_base"}

    def stats(F):
        w, V = np.linalg.eigh(F)
        C = np.linalg.inv(F)
        return w, V, np.sqrt(np.diag(C)), float(np.linalg.slogdet(F)[1]), C

    w0, V0, sd0, ld0, C0 = stats(F0)
    res = {}
    for name, S in cands.items():
        w, V, sd, ld, _ = stats(F0 + S.T @ S / args.sigma ** 2)
        res[name] = {"d_logdet": ld - ld0, "lam_min_gain": float(w[0] / w0[0]), "std_ratio": dict(zip(kn, (sd / sd0).tolist()))}
    pairs = []
    for a, b in itertools.combinations(cands, 2):
        S = np.vstack([cands[a], cands[b]])
        w, V, sd, ld, _ = stats(F0 + S.T @ S / args.sigma ** 2)
        pairs.append({"pair": [a, b], "d_logdet": ld - ld0, "lam_min_gain": float(w[0] / w0[0])})
    pairs.sort(key=lambda r: -r["d_logdet"])
    S_all = np.vstack(list(cands.values()))
    w_all, _, sd_all, ld_all, _ = stats(F0 + S_all.T @ S_all / args.sigma ** 2)

    lines = ["# E3b: identifiability on the subspace without Qf (Qf ≡ gate-WF shift)", "",
             f"- reduced Fisher (9 parameters): condition number {w0[-1] / w0[0]:.2e}; eigenvalues "
             + ", ".join(f"{x:.3g}" for x in w0), "",
             "## Sloppy directions (scaled units; typical change = 1)", ""]
    for j in range(3):
        v = V0[:, j]
        top = np.argsort(-np.abs(v))[:5]
        lines.append(f"- λ = {w0[j]:.3g} (1σ along it = {1 / np.sqrt(w0[j]):.0f} typical units): "
                     + ", ".join(f"{kn[k]} {v[k]:+.2f}" for k in top))
    lines += ["", "## 1σ (physical units)", "", "| parameter | calibration data | + all paper hold-out conditions |", "|---|---|---|"]
    for k, n in enumerate(kn):
        s = scale[keep[k]]
        lines.append(f"| {n} | {sd0[k] * s:.3g} | {sd_all[k] * s:.3g} |")
    lines += ["", "## OED ranking of the paper's hold-out conditions (Qf removed)", "",
              "| condition | Δ log det F | λmin gain | " + " | ".join(f"σ ratio {n}" for n in kn) + " |",
              "|---|---|---|" + "---|" * len(kn)]
    for name, r in sorted(res.items(), key=lambda kv: -kv[1]["d_logdet"]):
        lines.append(f"| {name} | {r['d_logdet']:.2f} | ×{r['lam_min_gain']:.2f} | "
                     + " | ".join(f"{r['std_ratio'][n]:.2f}" for n in kn) + " |")
    lines += ["", "Best pairs: " + "; ".join(f"{'+'.join(p['pair'])} (Δlogdet {p['d_logdet']:.2f}, λmin ×{p['lam_min_gain']:.2f})"
                                             for p in pairs[:3]),
              f"All five conditions: Δlogdet {ld_all - ld0:.2f}, λmin ×{w_all[0] / w0[0]:.2f}"]
    with open(f"{args.e3}/subspace.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    dump_json({"names": kn, "eigvals": w0, "std_phys_base": dict(zip(kn, (sd0 * scale[keep]).tolist())),
               "std_phys_all": dict(zip(kn, (sd_all * scale[keep]).tolist())), "candidates": res, "pairs": pairs,
               "all": {"d_logdet": ld_all - ld0, "lam_min_gain": float(w_all[0] / w0[0])}}, f"{args.e3}/subspace.json")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
