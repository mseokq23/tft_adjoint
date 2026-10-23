"""E6 - which part of the subgap DOS can the calibration transfer curves actually see?

Non-parametric version of the DOS identifiability question. The DOS is
represented by its 40 energy bins (acceptor-like and donor-like separately);
for every bin the adjoint gives dID/dN_k at every bias point with ONE transposed
solve (dR_i/dNtA_k = -b_i f_k, dR_i/dNtD_k = +b_i (1 - f_k), b_i = dR_i/dND).

1. Detection limit. For a DOS change of amplitude a (cm^-3 eV^-1, uniform in a
   0.1 eV window centred at Ec - E) the 1-sigma uncertainty from the
   log-current data (sigma = 0.02 dec) is
        sigma_a = 1 / || (I - P_n) s_a ||,
   where s_a is the noise-normalised sensitivity vector and P_n projects onto
   the nuisance directions that are calibrated anyway (S/D and gate WF, ND,eff,
   mu0). Without the projection the nuisance parameters are assumed known.
   A trap level whose occupancy never changes in the swept bias range acts as a
   fixed charge, i.e. exactly like ND,eff, and is projected out.
2. Fermi-level window. Ec - EF at the front (gate) and back interface,
   mid-channel, versus VG: the energies the sweep moves EF through.
3. Resolvable DOS components. Singular values of the projected, noise-normalised
   sensitivity matrix in the absolute-DOS basis (unit 1e18 cm^-3 eV^-1 per
   0.025 eV bin): components with singular value > 1 are constrained better
   than that unit.

Usage:  python scripts/e6_dos_detectability.py [--out results/e6] [--step 0.05]
"""
from __future__ import annotations

import argparse
import math

import devsim as ds

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from tftsim import PAPER_SETS, TFT
from tftsim.metrics import I_FLOOR, LN10
from tftsim.params import trap_levels
from tftsim.sensitivity import SensitivityEngine

NUIS = ["wf_sd", "wf_gate", "log10_nd_eff", "log10_mu0"]
G_REF = 1e18          # cm^-3 eV^-1, unit of the absolute DOS basis


def window_vectors(edges, width_eV=0.1):
    """Uniform-amplitude windows; returns (centres, M) with M[k, w] = overlap (eV) of bin k with window w."""
    lo = np.arange(0.0, edges[-1] - width_eV + 1e-9, width_eV / 2)          # 50 % overlapping windows
    M = np.zeros((len(edges) - 1, len(lo)))
    for w, a in enumerate(lo):
        b = a + width_eV
        M[:, w] = np.clip(np.minimum(edges[1:], b) - np.maximum(edges[:-1], a), 0, None)
    return lo + width_eV / 2, M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/e6")
    ap.add_argument("--step", type=float, default=0.05)
    ap.add_argument("--sigma", type=float, default=0.02)
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    t = TFT(PAPER_SETS["step5"]).build()
    eng = SensitivityEngine(t, central=True)
    d = t.dev
    rows = np.array(ds.get_equation_numbers(device=d, region="ch", equation="PotentialEquation"))
    dE, wA, wD, n1 = trap_levels(t.p, t.mat, t.edges)
    width = np.diff(t.edges)
    b = (eng.dR("log10_nd_eff") / (t.p.nd_eff * math.log(10)))[rows]
    scale_i = t.isign * t.width_cm
    x = np.array(ds.get_node_model_values(device=d, region="ch", name="x"))
    y = np.array(ds.get_node_model_values(device=d, region="ch", name="y"))
    xmid = t.geom.Lov_um * 1e-4 + 0.5 * t.geom.L_um * 1e-4
    col = np.where(np.abs(x - x[np.argmin(np.abs(x - xmid))]) < 1e-12)[0]
    i_front, i_back = col[np.argmin(y[col])], col[np.argmax(y[col])]
    vgs = vg_grid(-1.0, 2.5, args.step)
    data = {}
    for j, vd in enumerate((0.05, 1.0)):
        t.goto(vd=vd)
        seq = vgs if j % 2 == 0 else vgs[::-1]
        t.goto(vg=seq[0])
        rec = []
        for vg in seq:
            t.goto(vg=vg)
            lu, R0, n = eng.factorize()
            e = np.zeros(n)
            e[eng.iI] = 1.0
            lam_full = lu.solve(e, trans="T")
            lam = lam_full[rows]
            ne = np.array(ds.get_node_model_values(device=d, region="ch", name="Electrons"))
            f = ne[:, None] / (ne[:, None] + n1[None, :])
            gA = scale_i * ((lam * b)[:, None] * (-f)).sum(axis=0)
            gD = scale_i * ((lam * b)[:, None] * (1.0 - f)).sum(axis=0)
            gN = scale_i * np.array([lam_full @ eng.dR(nm, R0=R0) for nm in NUIS])
            rec.append((vg, t.drain_current(), gA, gD, gN, ne[i_front], ne[i_back]))
        rec.sort(key=lambda r: r[0])
        data[vd] = rec
        print(f"VD = {vd}: {len(rec)} bias points ({clock():.0f} s)")

    # ---------------------------------------------------------------- matrices
    blocks = {k: [] for k in ("A", "D", "N")}
    for vd, rec in data.items():
        for vg, i, gA, gD, gN, _, _ in rec:
            if i < I_FLOOR:
                continue
            c = 1.0 / (i * LN10 * args.sigma)          # d(log10 I)/sigma per unit
            blocks["A"].append(gA * c)
            blocks["D"].append(gD * c)
            blocks["N"].append(gN * c)
    SA, SD, SN = (np.array(blocks[k]) for k in ("A", "D", "N"))   # per cm^-3 of bin concentration (A, D); per unit theta (N)
    nrow = SA.shape[0]
    Q, _ = np.linalg.qr(SN)
    proj = lambda S: S - Q @ (Q.T @ S)

    centres, M = window_vectors(t.edges, 0.1)
    res = {}
    for key, S in (("acceptor", SA), ("donor", SD)):
        Sw = S @ M                                       # per (cm^-3 eV^-1) of window amplitude
        sig0 = 1.0 / np.maximum(np.linalg.norm(Sw, axis=0), 1e-300)
        sigp = 1.0 / np.maximum(np.linalg.norm(proj(Sw), axis=0), 1e-300)
        res[key] = {"sigma_known_nuisance": sig0, "sigma_profiled": sigp}
    # profile with ND,eff known but the other nuisance parameters free (shows the fixed-charge degeneracy)
    Q2, _ = np.linalg.qr(SN[:, [0, 1, 3]])
    for key, S in (("acceptor", SA), ("donor", SD)):
        Sw = S @ M
        res[key]["sigma_profiled_nd_known"] = 1.0 / np.maximum(np.linalg.norm(Sw - Q2 @ (Q2.T @ Sw), axis=0), 1e-300)

    # model DOS of Step 5 averaged over each window (cm^-3 eV^-1)
    gA_model = ((wA / width) @ M) / M.sum(axis=0)
    gD_model = ((wD / width) @ M) / M.sum(axis=0)
    # resolvable components in the absolute basis (unit G_REF per bin)
    fine = np.where(width < 0.03)[0]
    B = np.hstack([proj(SA[:, fine]) * width[fine] * G_REF, proj(SD[:, fine]) * width[fine] * G_REF])
    U, sv, Vt = np.linalg.svd(B, full_matrices=False)
    n_resolv = int(np.sum(sv > 1.0))

    # Fermi level window (Boltzmann: Ec - EF = Vt ln(Nc / n))
    ef = {}
    for vd, rec in data.items():
        vg = np.array([r[0] for r in rec])
        ef[vd] = {"vg": vg, "front": t.mat.Vt * np.log(t.mat.Nc_cm3 / np.array([r[5] for r in rec])),
                  "back": t.mat.Vt * np.log(t.mat.Nc_cm3 / np.array([r[6] for r in rec])),
                  "id": np.array([r[1] for r in rec])}

    np.savez(f"{out}/matrices.npz", SA=SA, SD=SD, SN=SN, dE=dE, width=width, wA=wA, wD=wD,
             centres=centres, M=M, sv=sv, Vt=Vt, nuis=np.array(NUIS))
    dump_json({"sigma_dec": args.sigma, "n_rows": nrow, "nuisance": NUIS, "window_centres_eV": centres,
               "result": res, "dos_model": {"acceptor": gA_model, "donor": gD_model},
               "singular_values": sv, "n_resolvable_unit_1e18_per_bin": n_resolv,
               "fermi": {str(k): v for k, v in ef.items()}, "wall_s": clock()}, f"{out}/summary.json")

    lines = ["# E6 DOS detectability from the calibration transfer curves (Step 5 model)", "",
             f"- {nrow} log-current residuals (L = 1 um, VD = 0.05/1.0 V, ID >= 10 pA), σ = {args.sigma} dec",
             f"- nuisance parameters profiled: {', '.join(NUIS)} (Qf is degenerate with the gate WF)",
             f"- singular values of the projected sensitivity (unit {G_REF:.0e} cm⁻³eV⁻¹ per 0.025 eV bin): "
             + ", ".join(f"{s:.3g}" for s in sv[:8]) + f" … → {n_resolv} components above 1", "",
             "1σ detection limit for a 0.1 eV-wide DOS change (cm⁻³eV⁻¹) and the Step-5 DOS in the same window:", "",
             "| Ec−E (eV) | acceptor: model DOS | σ (nuisance known) | σ (profiled) | σ (ND known) | donor: model DOS | σ (nuisance known) | σ (profiled) | σ (ND known) |",
             "|---|---|---|---|---|---|---|---|---|"]
    for w in range(0, len(centres), 2):
        a, dd = res["acceptor"], res["donor"]
        lines.append(f"| {centres[w]:.2f} | {gA_model[w]:.2g} | {a['sigma_known_nuisance'][w]:.2g} | {a['sigma_profiled'][w]:.2g} | "
                     f"{a['sigma_profiled_nd_known'][w]:.2g} | {gD_model[w]:.2g} | {dd['sigma_known_nuisance'][w]:.2g} | "
                     f"{dd['sigma_profiled'][w]:.2g} | {dd['sigma_profiled_nd_known'][w]:.2g} |")
    lines += ["", "Fermi level at mid-channel (VD = 0.05 V): Ec−EF front/back interface", ""]
    e = ef[0.05]
    for v in (-1.0, -0.5, 0.0, 0.5, 1.0, 2.5):
        k = int(np.argmin(np.abs(e["vg"] - v)))
        lines.append(f"- VG = {e['vg'][k]:+.2f} V: {e['front'][k]:+.3f} / {e['back'][k]:+.3f} eV (ID = {e['id'][k]:.2e} A)")
    with open(f"{out}/summary.md", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    if HAVE_MPL:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
        for k, (key, gm, c) in enumerate((("acceptor", gA_model, "C0"), ("donor", gD_model, "C3"))):
            r = res[key]
            ax[k].semilogy(centres, gm, "k-", lw=2, label="Step-5 model DOS")
            ax[k].semilogy(centres, r["sigma_known_nuisance"], c + ":", label="1σ, nuisance known")
            ax[k].semilogy(centres, r["sigma_profiled_nd_known"], c + "--", label="1σ, WF/µ0 profiled (ND known)")
            ax[k].semilogy(centres, r["sigma_profiled"], c + "-", label="1σ, WF/ND/µ0 profiled")
            ax[k].set_ylim(1e14, 1e23)
            ax[k].set_xlim(0, 1.6)
            ax[k].set_xlabel("Ec − E (eV)")
            ax[k].set_ylabel("DOS / detection limit (cm⁻³ eV⁻¹)")
            ax[k].set_title(f"{key}-like traps: detectable where σ < DOS")
            ax[k].legend(fontsize=8)
        for vd, sty in ((0.05, "-"), (1.0, "--")):
            e = ef[vd]
            ax[2].plot(e["vg"], e["front"], "C0" + sty, label=f"front interface, VD = {vd} V")
            ax[2].plot(e["vg"], e["back"], "C1" + sty, label=f"back interface, VD = {vd} V")
        ax[2].invert_yaxis()
        ax[2].set_xlabel("VG (V)")
        ax[2].set_ylabel("Ec − EF at mid-channel (eV)")
        ax[2].set_title("Fermi-level window swept by the transfer curve")
        ax[2].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(f"{out}/e6_dos_detectability.png", dpi=150)
    print(f"done in {clock():.0f} s -> {out}")


if __name__ == "__main__":
    main()
