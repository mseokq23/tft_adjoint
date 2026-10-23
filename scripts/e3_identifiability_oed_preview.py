"""E3/E4 preview - local identifiability and gradient-based OED at the calibrated model.

1. Sweeps the final 2%-W model (arXiv:2609.12184) at the calibration conditions
   (L = 1 um, VD = 0.05 and 1.0 V) and computes d log10(ID)/d theta at every
   bias point above the current floor with the adjoint engine.
2. Builds the Fisher information F = S^T S / sigma^2 (theta scaled by the
   typical changes in params.THETA_DEFS, sigma = assumed noise of digitized
   log10 ID) and reports eigen-spectrum, condition number, parameter standard
   deviations and correlations (local, linearised; a profile-likelihood check
   comes later).
3. OED preview: adds each hold-out condition published in the paper
   (VD = 0.1/0.2/0.5 V at L = 1 um; L = 2 um at VD = 0.05/1.0 V) and ranks them by
   D-optimal (log det) and E-optimal (smallest eigenvalue) information gain.

Usage:  python scripts/e3_identifiability_oed_preview.py [--sigma 0.02] [--out results/e3]
"""
from __future__ import annotations

import argparse
import itertools
import math

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from tftsim import PAPER_SETS, THETA_ALL, TFT, Geometry
from tftsim.metrics import I_FLOOR, LN10, extract
from tftsim.params import THETA_DEFS, theta_scale
from tftsim.sensitivity import SensitivityEngine

NAMES = THETA_ALL


def sweep_with_grad(params, geom, vds, vgs):
    t = TFT(params, geom).build()
    eng = SensitivityEngine(t, central=False)
    data = {}
    for j, vd in enumerate(vds):
        seq = vgs if j % 2 == 0 else vgs[::-1]

        def cb(tft, rec):
            if rec["id"] >= 0.1 * I_FLOOR:
                rec["g"] = eng.current_gradient(NAMES, "adjoint")

        recs = t.transfer(vd, seq, on_point=cb)
        recs.sort(key=lambda r: r["vg"])
        vg = np.array([r["vg"] for r in recs])
        i = np.array([r["id"] for r in recs])
        G = np.array([r.get("g", np.zeros(len(NAMES))) for r in recs])
        data[vd] = (vg, i, G)
    return data, eng.report(), t.info()


def log_rows(vg, i, G):
    m = i >= I_FLOOR
    return G[m] / (i[m, None] * LN10)


def fim_stats(F, scale):
    w, V = np.linalg.eigh(F)
    order = np.argsort(w)
    w, V = w[order], V[:, order]
    Finv = np.linalg.pinv(F, rcond=1e-14)
    std_scaled = np.sqrt(np.maximum(np.diag(Finv), 0))
    d = np.sqrt(np.outer(np.diag(Finv), np.diag(Finv)))
    corr = np.where(d > 0, Finv / np.where(d > 0, d, 1), 0)
    sign, logdet = np.linalg.slogdet(F + 1e-12 * np.eye(len(F)))
    return {"eigvals": w, "eigvecs": V, "cond": float(w[-1] / max(w[0], 1e-300)),
            "logdet": float(logdet), "std_scaled": std_scaled, "std_phys": std_scaled * scale, "corr": corr}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, default=0.02, help="noise of log10(ID) in decades")
    ap.add_argument("--out", default="results/e3")
    ap.add_argument("--step", type=float, default=0.05)
    ap.add_argument("--vgmin", type=float, default=-1.0, help="sweep start; deeper off-state is below the 10 pA floor")
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    vgs = vg_grid(vmin=args.vgmin, step=args.step)
    scale = theta_scale(NAMES)
    p = PAPER_SETS["step5"]

    # 1) calibration conditions
    base, rep_base, info = sweep_with_grad(p, Geometry(), (0.05, 1.0), vgs)
    print(f"base sweeps done ({clock():.0f} s) engine {rep_base}")
    S_base = np.vstack([log_rows(*base[vd]) for vd in (0.05, 1.0)]) * scale
    F_base = S_base.T @ S_base / args.sigma ** 2
    st = fim_stats(F_base, scale)

    # metric sensitivities at VD = 1 V (physical insight table)
    vg, i, G = base[1.0]
    mt = extract(vg, i, G)
    msens = {nm: {"dVth_mV": 1e3 * mt["d_vth"][k] * scale[k],
                  "dSS_mVdec": mt["d_ss_mVdec"][k] * scale[k],
                  "dlogIon_dec": mt["d_log10_ion"][k] * scale[k]} for k, nm in enumerate(NAMES)}

    # 2) hold-out candidates
    cand_data = {}
    extra, rep1, _ = sweep_with_grad(p, Geometry(), (0.1, 0.2, 0.5), vgs)
    for vd in (0.1, 0.2, 0.5):
        cand_data[f"L1_VD{vd}"] = extra[vd]
    l2, rep2, _ = sweep_with_grad(p, Geometry(L_um=2.0), (0.05, 1.0), vgs)
    for vd in (0.05, 1.0):
        cand_data[f"L2_VD{vd}"] = l2[vd]
    print(f"candidate sweeps done ({clock():.0f} s)")

    cands = {}
    for name, d in cand_data.items():
        S_c = log_rows(*d) * scale
        F_c = F_base + S_c.T @ S_c / args.sigma ** 2
        sc = fim_stats(F_c, scale)
        cands[name] = {"d_logdet": sc["logdet"] - st["logdet"],
                       "lam_min_ratio": float(sc["eigvals"][0] / st["eigvals"][0]),
                       "cond": sc["cond"],
                       "std_ratio": (sc["std_scaled"] / st["std_scaled"]).tolist(),
                       "n_rows": int(S_c.shape[0])}
    # best pair (greedy check)
    pairs = []
    keys = list(cand_data)
    for a, b in itertools.combinations(keys, 2):
        F_ab = F_base.copy()
        for k in (a, b):
            S_c = log_rows(*cand_data[k]) * scale
            F_ab += S_c.T @ S_c / args.sigma ** 2
        pairs.append({"pair": [a, b], "d_logdet": fim_stats(F_ab, scale)["logdet"] - st["logdet"]})
    pairs.sort(key=lambda r: -r["d_logdet"])

    sloppy = []
    for j in range(3):
        v = st["eigvecs"][:, j]
        top = np.argsort(-np.abs(v))[:4]
        sloppy.append({"eigval": float(st["eigvals"][j]),
                       "components": {NAMES[k]: float(v[k]) for k in top}})
    strong_corr = [(NAMES[a], NAMES[b], float(st["corr"][a, b]))
                   for a in range(len(NAMES)) for b in range(a + 1, len(NAMES)) if abs(st["corr"][a, b]) > 0.9]

    # sub-model identifiability: what if the perfectly translating qf is fixed (or merged)?
    def sub_std(F, keep):
        idx = [NAMES.index(n) for n in keep]
        Fs = F[np.ix_(idx, idx)]
        Finv = np.linalg.pinv(Fs, rcond=1e-14)
        return dict(zip(keep, np.sqrt(np.maximum(np.diag(Finv), 0))))
    no_qf = [n for n in NAMES if n != "qf"]
    no_dos_qf = [n for n in no_qf if n not in ("log10_n_gd", "e_gd", "s_gd", "log10_n_ta", "w_ta")]
    sub = {"no_qf": sub_std(F_base, no_qf), "no_qf_no_dos": sub_std(F_base, no_dos_qf),
           "dos_only": sub_std(F_base, ["log10_n_gd", "e_gd", "s_gd", "log10_n_ta", "w_ta"])}
    np.savez(f"{out}/rows.npz", S_base=S_base, scale=scale, names=np.array(NAMES),
             **{f"S_{k}": log_rows(*v) * scale for k, v in cand_data.items()})

    summary = {
        "sub_model_std_scaled": sub,
        "sigma_dec": args.sigma, "device": info, "names": NAMES,
        "scales": {n: THETA_DEFS[n][2] for n in NAMES},
        "base": {"n_rows": int(S_base.shape[0]), "eigvals": st["eigvals"], "cond": st["cond"],
                 "std_scaled": dict(zip(NAMES, st["std_scaled"])), "std_phys": dict(zip(NAMES, st["std_phys"])),
                 "sloppy_directions": sloppy, "strong_correlations": strong_corr,
                 "corr": st["corr"]},
        "metric_sensitivity_per_typical_change": msens,
        "candidates": cands, "best_pairs": pairs[:5],
        "engine": {"base": rep_base, "extra": rep1, "L2": rep2}, "wall_s": clock(),
    }
    dump_json(summary, f"{out}/summary.json")

    lines = ["# E3/E4 preview: identifiability and OED (final 2%-W model)", "",
             f"- assumed noise sigma = {args.sigma} dec in log10(ID); {S_base.shape[0]} calibration rows "
             f"(L = 1 um, VD = 0.05/1.0 V, ID >= 10 pA)",
             f"- Fisher condition number = {st['cond']:.2e}", "",
             "| parameter | typical change | 1-sigma (scaled) | 1-sigma (physical) |", "|---|---|---|---|"]
    for k, n in enumerate(NAMES):
        lines.append(f"| {n} | {scale[k]:g} | {st['std_scaled'][k]:.3g} | {st['std_phys'][k]:.3g} |")
    lines += ["", "## Sloppiest directions (smallest eigenvalues)", ""]
    for s in sloppy:
        comp = ", ".join(f"{k} {v:+.2f}" for k, v in s["components"].items())
        lines.append(f"- λ = {s['eigval']:.3g}: {comp}")
    lines += ["", "## |correlation| > 0.9", ""] + [f"- {a} vs {b}: {c:+.3f}" for a, b, c in strong_corr]
    lines += ["", "## 1-sigma (scaled) of sub-models", "", "| parameter | all | qf fixed | qf+DOS fixed | DOS only |", "|---|---|---|---|---|"]
    for k, n in enumerate(NAMES):
        cells = [f"{st['std_scaled'][k]:.3g}"] + [f"{sub[m][n]:.3g}" if n in sub[m] else "-" for m in ("no_qf", "no_qf_no_dos", "dos_only")]
        lines.append(f"| {n} | " + " | ".join(cells) + " |")
    lines += ["", "## OED preview: adding one hold-out condition", "",
              "| condition | Δ log det F | λmin gain | cond |", "|---|---|---|---|"]
    for name, c in sorted(cands.items(), key=lambda kv: -kv[1]["d_logdet"]):
        lines.append(f"| {name} | {c['d_logdet']:.2f} | ×{c['lam_min_ratio']:.2f} | {c['cond']:.2e} |")
    lines += ["", "Best pairs: " + "; ".join(f"{'+'.join(r['pair'])} ({r['d_logdet']:.2f})" for r in pairs[:3])]
    lines += ["", "## Metric sensitivity per typical change (VD = 1 V)", "",
              "| parameter | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |", "|---|---|---|---|"]
    for n, v in msens.items():
        lines.append(f"| {n} | {v['dVth_mV']:+.1f} | {v['dSS_mVdec']:+.2f} | {v['dlogIon_dec']:+.4f} |")
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")

    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
        ax[0].semilogy(range(1, len(NAMES) + 1), st["eigvals"], "o-", label="calibration curves only")
        F_all = F_base.copy()
        for d in cand_data.values():
            S_c = log_rows(*d) * scale
            F_all += S_c.T @ S_c / args.sigma ** 2
        ax[0].semilogy(range(1, len(NAMES) + 1), np.sort(np.linalg.eigvalsh(F_all)), "s--", label="+ all hold-out conditions")
        ax[0].set_xlabel("eigen-direction (sorted)")
        ax[0].set_ylabel("Fisher eigenvalue")
        ax[0].set_title("Information spectrum")
        ax[0].legend()
        im = ax[1].imshow(st["corr"], vmin=-1, vmax=1, cmap="RdBu_r")
        ax[1].set_xticks(range(len(NAMES)))
        ax[1].set_xticklabels(NAMES, rotation=70, fontsize=8)
        ax[1].set_yticks(range(len(NAMES)))
        ax[1].set_yticklabels(NAMES, fontsize=8)
        ax[1].set_title("Parameter correlation (calibration data)")
        fig.colorbar(im, ax=ax[1], fraction=0.046)
        names_c = sorted(cands, key=lambda k: -cands[k]["d_logdet"])
        ax[2].barh(names_c[::-1], [cands[k]["d_logdet"] for k in names_c[::-1]])
        ax[2].set_xlabel("Δ log det F (D-optimal gain)")
        ax[2].set_title("Which extra measurement helps most?")
        fig.tight_layout()
        fig.savefig(f"{out}/e3_identifiability_oed.png", dpi=150)
    print(f"done in {clock():.0f} s -> {out}")
    print(open(f"{out}/summary.md").read())


if __name__ == "__main__":
    main()
