"""E5b - merge profile-likelihood runs that were split into parallel halves (no new simulations).

e5_profile_likelihood.py walks outward from the truth on each side; running the
two sides as separate processes halves the wall time. This script merges their
summary.json files into one table and figure.

Usage:  python scripts/e5b_merge_profiles.py results/e5_nd_up results/e5_nd_down [--out results/e5_nd]
"""
from __future__ import annotations

import argparse
import json

from _common import HAVE_MPL, dump_json, np, outdir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--out", default="results/e5_nd")
    args = ap.parse_args()
    out = outdir(args.out)
    parts = [json.load(open(f"{d}/summary.json")) for d in args.dirs]
    prof, names, sigma, vg = parts[0]["param"], parts[0]["names"], parts[0]["sigma_dec"], np.array(parts[0]["vg"])
    recs, curves, failed = {}, {}, []
    for p in parts:
        assert p["param"] == prof and p["names"] == names
        for r in p["records"]:
            recs[round(r["value"], 6)] = r
        curves.update(p["curves"])
        failed += p.get("failed", [])
    ok = [recs[k] for k in sorted(recs)]
    free = [n for n in names if n != prof]
    v_true = next(r["value"] for r in ok if r["chi2"] == 0.0)
    lines = [f"# E5 profile likelihood of {prof} (merged; Step 5 twin, L = 1 µm, VD = 0.05/1.0 V, σ = {sigma} dec)", "",
             f"{parts[0]['n_res']} residuals. χ² < 3.84 → indistinguishable from the truth at 95 % (1 dof).", "",
             f"| {prof} | χ² | RMS misfit (dec) | max misfit (dec) | " + " | ".join(free) + " |",
             "|---|---|---|---|" + "---|" * len(free)]
    for r in ok:
        lines.append(f"| {r['value']:.3f} | {r['chi2']:.3g} | {r['rms_dec']:.4f} | {r['max_dec']:.4f} | "
                     + " | ".join(f"{r['theta'][n]:.4g}" for n in free) + " |")
    ho = [r for r in ok if isinstance(r.get("holdout"), dict) and "error" not in r["holdout"]]
    if ho:
        keys = list(ho[0]["holdout"])
        lines += ["", "Re-fitted sets predicting the paper's hold-out conditions (vs the truth): RMS (dec) / ΔVth (mV) / ΔSS (mV/dec)", "",
                  f"| {prof} | " + " | ".join(keys) + " |", "|---|" + "---|" * len(keys)]
        for r in ho:
            h = r["holdout"]
            lines.append(f"| {r['value']:.3f} | " + " | ".join(
                f"{h[k]['rms_dec']:.4f} / {h[k]['dVth_mV']:+.1f} / {h[k]['dSS_mVdec']:+.1f}" for k in keys) + " |")
    if failed:
        lines += ["", "Failed fits: " + ", ".join(f"{f['value']}: {f['error'][:80]}" for f in failed)]
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    dump_json({"param": prof, "names": names, "sigma_dec": sigma, "vg": vg, "records": ok, "failed": failed,
               "curves": curves, "v_true": v_true}, f"{out}/summary.json")
    print("\n".join(lines))
    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
        xv = [r["value"] for r in ok]
        ax[0].plot(xv, [r["chi2"] for r in ok], "o-")
        ax[0].axhline(3.84, color="r", ls="--", lw=0.8, label="95 % threshold (Δχ² = 3.84)")
        ax[0].axvline(v_true, color="k", lw=0.6)
        ax[0].set_xlabel(f"{prof} (fixed)")
        ax[0].set_ylabel("profile χ² (calibration curves)")
        ax[0].legend()
        for k, v in enumerate(sorted(curves, key=float)):
            for vd, i in curves[v].items():
                ax[1].semilogy(vg, np.maximum(np.array(i), 1e-16), "-" if float(vd) == 1.0 else "--", lw=1.0,
                               color=f"C{k}", label=f"{prof} = {float(v):.2f}" if float(vd) == 1.0 else None)
        ax[1].axhline(1e-11, color="grey", lw=0.6, ls=":")
        ax[1].set_ylim(1e-13, 1e-2)
        ax[1].set_xlabel("VG (V)")
        ax[1].set_ylabel("ID (A)  (solid: VD = 1 V, dashed: 0.05 V)")
        ax[1].legend(fontsize=7)
        if ho:
            for key in keys:
                ax[2].plot([r["value"] for r in ho], [r["holdout"][key]["rms_dec"] for r in ho], "o-", label=key)
            ax[2].axhline(sigma, color="r", ls="--", lw=0.8, label=f"σ = {sigma} dec")
            ax[2].set_xlabel(f"{prof} (fixed)")
            ax[2].set_ylabel("hold-out RMS deviation vs truth (dec)")
            ax[2].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(f"{out}/e5_profile.png", dpi=150)


if __name__ == "__main__":
    main()
