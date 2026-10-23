"""Auxiliary analysis - a textbook 1-D gate-stack model with a Pao-Sah drain current.

Why: it explains analytically why drain-bias and channel-length transfer tests carry no
information about the descriptors (E14), and it is the planned analytic basis of E19.

Model (NOT the 2-D DEVSIM model used in the experiments; numpy only)
  stack      gate | HfO2 7 nm (eps 20) | IWO 3.5 nm (eps 10) | HfO2 passivation 5 nm, zero field on top
  charges    electrons n = Nc exp((psi - V)/Vt) (Boltzmann), fully ionised ND,eff, sheet charge Qf at the
             gate-dielectric/IWO interface, acceptor-like tail + donor-like Gaussian DOS with
             f = n/(n + n1) on a 0.005 eV grid; chi = 4.16 eV, Nc = 5e18 cm^-3 (as in tftsim)
  current    I = q mu (W/L) int_0^VD N_s(VG, V) dV (gradual channel, constant mobility, no contacts,
             no overlap regions, no 2-D effects)

Exact property of this model: every charge in the film is a function of (psi - V), and the fixed
charges (ND,eff, Qf) do not depend on psi, so N_s(VG, V) = N_s(VG - V, 0). Hence
    I(VG, VD) = q mu (W/L) int_{VG-VD}^{VG} N_s(u, 0) du.
All drain-bias and channel-length transfer curves follow from ONE function N_s(u): once the
calibration curves fix it over the swept range, any parameter set that reproduces it passes the
VD and L tests automatically (L only scales the current by 1/L). A thickness change (t_IWO, t_ox)
changes N_s(u) itself, and C-V sees a different functional (trapped charge included).

Runs
  1. Step 5 in the 1-D model vs the 2-D DEVSIM model (results/e1/summary.json) at L = 1 / 2 um,
     VD = 0.05 / 1 V, and the subthreshold slope without traps / with the tail only.
  2. The re-fitted twin sets of E5 (results/e5_nd/summary.json; theta fitted with DEVSIM) in the
     1-D model: RMS difference from the Step-5 truth on the calibration curves, dVth at L = 2 um,
     and dVth for t_IWO = 5 / 7 nm and t_ox = 14 nm (VD = 1 V).

Usage:  python scripts/aux_paosah_1d.py [--out results/aux_paosah_1d]      (about 10-20 min, 1 core)
"""
from __future__ import annotations

import argparse
import json
import math

from _common import Timer, dump_json, np, outdir

Q = 1.602176634e-19
EPS0 = 8.8541878128e-14
VT = 1.380649e-23 * 300.0 / Q
NC = 5e18
CHI = 4.16

STEP5 = dict(wf_g=4.00, qf=1.59e10, nd=2.56e18, n_gd=2.19e19, e_gd=0.12, s_gd=0.08,
             n_ta=7.8e19, w_ta=0.133, mu=24.0)


class Stack:
    """1-D Poisson solver for the gate stack (same sign conventions as tftsim/device.py:
    Ec = -q psi in the channel, gate: psi = VG - (WF_gate - chi))."""

    def __init__(self, t_ox=7.0, t_ch=3.5, t_pas=5.0, eps_ox=20.0, eps_ch=10.0, **p):
        self.t_ox, self.t_ch, self.t_pas = t_ox, t_ch, t_pas
        par = dict(STEP5)
        par.update(p)
        self.p = par
        nm = 1e-7
        y_ox = np.linspace(0, t_ox, int(round(t_ox / 0.25)) + 1)
        y_ch = np.linspace(t_ox, t_ox + t_ch, int(round(t_ch / 0.05)) + 1)
        y_pa = np.linspace(t_ox + t_ch, t_ox + t_ch + t_pas, int(round(t_pas / 0.25)) + 1)
        self.y = np.concatenate([y_ox, y_ch[1:], y_pa[1:]]) * nm
        y = self.y
        h = np.diff(y)
        mid = 0.5 * (y[:-1] + y[1:])
        in_ch_edge = (mid > t_ox * nm) & (mid < (t_ox + t_ch) * nm)
        self.eps_e = np.where(in_ch_edge, eps_ch, eps_ox) * EPS0
        self.h = h
        lo = np.concatenate([[y[0]], mid])                      # control volume of each node
        hi = np.concatenate([mid, [y[-1]]])
        a, b = t_ox * nm, (t_ox + t_ch) * nm
        self.vch = np.clip(np.minimum(hi, b) - np.maximum(lo, a), 0, None)   # channel part of it
        self.i_if = int(np.argmin(abs(y - a)))                  # gate-dielectric/IWO interface node
        self.ch = self.vch > 0
        dE = np.arange(0.0025, 1.6, 0.005)                      # trap grid, Ec - E
        g_a = par["n_ta"] * np.exp(-dE / par["w_ta"]) if par["n_ta"] > 0 else 0 * dE
        g_d = (par["n_gd"] * np.exp(-(dE - par["e_gd"]) ** 2 / (2 * par["s_gd"] ** 2))
               if par["n_gd"] > 0 else 0 * dE)
        self.wA = g_a * 0.005
        self.wD = g_d * 0.005
        self.n1 = NC * np.exp(-dE / VT)
        self.psi = None
        self.cox = eps_ox * EPS0 / (t_ox * nm)

    def charges(self, psi, V):
        n = np.zeros_like(psi)
        A = np.zeros_like(psi)
        D = np.zeros_like(psi)
        dndp = np.zeros_like(psi)
        dA = np.zeros_like(psi)
        dD = np.zeros_like(psi)
        c = self.ch
        nn = NC * np.exp(np.clip((psi[c] - V) / VT, -200, 60))
        n[c] = nn
        dndp[c] = nn / VT
        den = nn[:, None] + self.n1[None, :]
        f = nn[:, None] / den
        A[c] = (self.wA[None, :] * f).sum(1)                    # filled acceptors (negative)
        D[c] = (self.wD[None, :] * (1 - f)).sum(1)              # empty donors (positive)
        dfdn = self.n1[None, :] / den ** 2
        dA[c] = (self.wA[None, :] * dfdn).sum(1) * dndp[c]
        dD[c] = -(self.wD[None, :] * dfdn).sum(1) * dndp[c]
        return n, A, D, dndp, dA, dD

    def solve(self, VG, V=0.0, psi0=None, tol=1e-12, maxit=200):
        p = self.p
        N = len(self.y)
        psi = (psi0.copy() if psi0 is not None else
               (self.psi.copy() if self.psi is not None else np.full(N, VG - (p["wf_g"] - CHI))))
        psi[0] = VG - (p["wf_g"] - CHI)
        nd_vec = np.where(self.ch, p["nd"], 0.0)
        g = self.eps_e / self.h
        for _ in range(maxit):
            n, A, D, dndp, dA, dD = self.charges(psi, V)
            rho = Q * (nd_vec - n + D - A) * self.vch
            rho[self.i_if] += Q * p["qf"]
            drho = Q * (-dndp + dD - dA) * self.vch
            F = self.eps_e * np.diff(psi) / self.h              # flux between nodes
            R = np.zeros(N)
            R[1:-1] = F[1:] - F[:-1] + rho[1:-1]
            R[-1] = -F[-1] + rho[-1]                            # zero field on top of the passivation
            J = np.zeros((N, N))
            for i in range(1, N):
                gl = g[i - 1]
                gr = g[i] if i < N - 1 else 0.0
                J[i, i] = -gl - gr + drho[i]
                J[i, i - 1] = gl
                if i < N - 1:
                    J[i, i + 1] = gr
            J[0, 0] = 1.0
            R[0] = 0.0
            d = np.linalg.solve(J, -R)
            m = np.max(abs(d))
            if m > 0.1:
                d *= 0.1 / m
            psi += d
            if m < tol:
                break
        self.psi = psi
        return psi

    def sheet_n(self, VG, V=0.0):
        psi = self.solve(VG, V)
        n, A, D, *_ = self.charges(psi, V)
        return (n * self.vch).sum(), (A * self.vch).sum(), (D * self.vch).sum(), psi


def ns_curve(stack, us):
    """Equilibrium sheet densities (free electrons, filled acceptors, empty donors) vs u = VG - V."""
    stack.psi = None
    return np.array([stack.sheet_n(u, 0.0)[:3] for u in us])


def iv_from_ns(us, ns, vgs, VD, mu=24.0, W_um=5.0, L_um=1.0):
    """I(VG, VD) = q mu W/L * int_{VG-VD}^{VG} N_s(u) du (exact for this 1-D model)."""
    cum = np.concatenate([[0], np.cumsum(0.5 * (ns[1:] + ns[:-1]) * np.diff(us))])
    hi = np.interp(vgs, us, cum)
    lo = np.interp(vgs - VD, us, cum)
    return Q * mu * (W_um / L_um) * (hi - lo)


def crossing(vg, i, target):
    L = np.log10(np.maximum(i, 1e-300))
    lt = math.log10(target)
    idx = np.where((L[:-1] < lt) & (L[1:] >= lt))[0]
    if len(idx) == 0:
        return math.nan
    a = int(idx[0])
    t = (lt - L[a]) / (L[a + 1] - L[a])
    return vg[a] + t * (vg[a + 1] - vg[a])


def metrics(vg, i):
    """Vth at 100 nA, SS = average inverse slope between 0.1 and 100 nA (mV/dec), Ion = max."""
    vth = crossing(vg, i, 100e-9)
    vlo = crossing(vg, i, 0.1e-9)
    return {"vth": vth, "ss_mVdec": (vth - vlo) / 3 * 1e3, "ion": float(np.max(i))}


def pars(th, qf=1.59e10):
    return dict(wf_g=th["wf_gate"], qf=qf, nd=10 ** th["log10_nd_eff"], n_gd=10 ** th["log10_n_gd"],
                e_gd=th["e_gd"], s_gd=th["s_gd"], n_ta=10 ** th["log10_n_ta"], w_ta=th["w_ta"],
                mu=10 ** th["log10_mu0"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/aux_paosah_1d")
    ap.add_argument("--e1", default="results/e1/summary.json")
    ap.add_argument("--e5", default="results/e5_nd/summary.json")
    args = ap.parse_args()
    out = outdir(args.out)
    clock = Timer()
    us = np.arange(-3.2, 3.0001, 0.005)
    vgs = np.arange(-1.5, 2.5001, 0.01)
    res = {"model": "1-D gate stack + Pao-Sah current (illustration, not the 2-D DEVSIM model)"}

    # 1. Step 5: 1-D vs 2-D DEVSIM
    e1 = json.load(open(args.e1))["results"]
    ns5 = ns_curve(Stack(), us)[:, 0]
    rows = []
    for L, name in ((1.0, "step5"), (2.0, "step5_L2")):
        for vd in (0.05, 1.0):
            m1 = metrics(vgs, iv_from_ns(us, ns5, vgs, vd, L_um=L))
            m2 = e1[name]["metrics"][str(vd)]
            rows.append({"L_um": L, "vd": vd, "oned": m1, "devsim": {k: m2[k] for k in ("vth", "ss_mVdec", "ion")}})
    ss_variants = {}
    for tag, kw in (("no traps", dict(n_ta=0.0, n_gd=0.0)), ("tail only", dict(n_gd=0.0)), ("tail + Gaussian", {})):
        ns = ns_curve(Stack(**kw), us)[:, 0]
        ss_variants[tag] = metrics(vgs, iv_from_ns(us, ns, vgs, 1.0))["ss_mVdec"]
    res["step5"] = rows
    res["ss_vd1_by_dos"] = ss_variants
    print(f"step-5 comparison done ({clock():.0f} s)", flush=True)

    # 2. the E5 re-fitted sets in the 1-D model
    d = json.load(open(args.e5))
    curves = {}
    for r in d["records"]:
        p = pars(r["theta"])
        o = {}
        for tag, geo in (("base", {}), ("tch5", dict(t_ch=5.0)), ("tch7", dict(t_ch=7.0)), ("tox14", dict(t_ox=14.0))):
            ns = ns_curve(Stack(**geo, **p), us)[:, 0]
            for vd in (0.05, 1.0):
                o[(tag, vd)] = iv_from_ns(us, ns, vgs, vd, mu=p["mu"])
                if tag == "base":
                    o[("L2", vd)] = iv_from_ns(us, ns, vgs, vd, mu=p["mu"], L_um=2.0)
        curves[round(r["value"], 4)] = o
        print(f"[ND,eff = 10^{r['value']:.2f}] done ({clock():.0f} s)", flush=True)
    truth = min(curves, key=lambda v: abs(v - d["v_true"]))
    ref = curves[truth]
    sets = {}
    for v, o in sorted(curves.items()):
        row = {}
        cal = []
        for vd in (0.05, 1.0):
            a, b = o[("base", vd)], ref[("base", vd)]
            m = b >= 1e-11
            cal.append(np.log10(np.maximum(a[m], 1e-30)) - np.log10(b[m]))
        row["calib_rms_dec"] = float(np.sqrt(np.mean(np.concatenate(cal) ** 2)))
        for key in (("L2", 0.05), ("L2", 1.0), ("tch5", 1.0), ("tch7", 1.0), ("tox14", 1.0)):
            row[f"dVth_{key[0]}_VD{key[1]:g}_mV"] = 1e3 * (metrics(vgs, o[key])["vth"] - metrics(vgs, ref[key])["vth"])
        sets[f"{v:.3f}"] = row
    res["e5_sets"] = sets
    res["truth"] = truth
    dump_json(res, f"{out}/summary.json")

    lines = ["# Auxiliary: 1-D gate stack + Pao-Sah current", "",
             "Illustration model (numpy), not the 2-D DEVSIM model. Exact property: N_s(VG, V) = N_s(VG - V, 0), so "
             "I(VG, VD) = q µ W/L ∫_{VG-VD}^{VG} N_s(u) du for every VD and L.", "",
             "## Step 5: 1-D model vs 2-D DEVSIM (results/e1)", "",
             "| L (µm) | VD (V) | DEVSIM Vth / SS / Ion | 1-D Vth / SS / Ion |", "|---|---|---|---|"]
    for r in rows:
        a, b = r["devsim"], r["oned"]
        lines.append(f"| {r['L_um']:g} | {r['vd']:g} | {a['vth']:+.3f} V / {a['ss_mVdec']:.1f} / {a['ion']:.2e} A | "
                     f"{b['vth']:+.3f} V / {b['ss_mVdec']:.1f} / {b['ion']:.2e} A |")
    lines += ["", "SS at VD = 1 V (1-D): " + ", ".join(f"{k} {v:.0f} mV/dec" for k, v in ss_variants.items()), "",
              "## E5 re-fitted sets (θ from DEVSIM) in the 1-D model, relative to the Step-5 truth", "",
              "| log10 ND,eff | calibration RMS (dec) | ΔVth L 2 µm VD 0.05 / 1 V (mV) | ΔVth t_IWO 5 nm | ΔVth t_IWO 7 nm | ΔVth t_ox 14 nm |",
              "|---|---|---|---|---|---|"]
    for v, r in sets.items():
        lines.append(f"| {v} | {r['calib_rms_dec']:.4f} | {r['dVth_L2_VD0.05_mV']:+.1f} / {r['dVth_L2_VD1_mV']:+.1f} | "
                     f"{r['dVth_tch5_VD1_mV']:+.0f} | {r['dVth_tch7_VD1_mV']:+.0f} | {r['dVth_tox14_VD1_mV']:+.0f} |")
    lines += ["", "Thickness tests at VD = 1 V. The 2-D DEVSIM predictions for the same sets are in results/e9."]
    with open(f"{out}/summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"done in {clock():.0f} s")


if __name__ == "__main__":
    main()
