#!/usr/bin/env python
"""E22: parameter-ambiguity analysis on the open demo of the reference paper's workflow.

The demo (G. Jeong, "Agentic TCAD Calibration: Reproducible Open Demo",
https://github.com/gyulab/agentic-TCAD-calibration-public, MIT License, commit ccaa678) is a
closed-form compact-model version of the paper's calibration loop on synthetic data. It is used
as installed: nothing in it is modified, and its code is imported, not copied (the eligibility
checks repeated in JointReachController are noted there). This script adds, on the demo's own
data and model:

  A. the exact null direction of the calibration data at the demo's final parameters
     (gate work function phi_G vs. fixed charge Q_f) and its extent inside the approved bounds,
  B. the spread of held-out (V_D, L) and structure-split (t_ox, t_ch) predictions across it,
  C. profile likelihoods of Q_f with the demo's calibration set and with one extra
     t_ox = 14 nm curve,
  D. path dependence: the demo's own controller rerun from other starting points and trust
     radii, with and without the extra t_ox curve,
  E. an unknown interface-trap term (C_it) next to the tail-state density N_tail: their exact
     degeneracy at one channel thickness, the profile likelihood of N_tail with and without a
     t_ch = 7 nm curve, and the controller's model selection for two activation tests (the demo's,
     and the demo's OR a joint-reach test) and two update groupings.

Assumptions for the thickness splits (not part of the demo): C_it and Q_f are interface
properties, so they do not change with t_ox or t_ch; n_interface = C_it / C_ox therefore scales
with t_ox. The synthetic truth keeps the demo's gate-dependent mobility and DIBL term.

Environment (separate from this repository's DEVSIM environment; the demo pins Python 3.11):
    git clone https://github.com/gyulab/agentic-TCAD-calibration-public demo && git -C demo checkout ccaa678
    python3.11 -m venv demo_venv
    demo_venv/bin/pip install -r demo/requirements-lock.txt && demo_venv/bin/pip install -e demo --no-deps

Usage (from the repository root, about 2 minutes on 2 vCPU):
    demo_venv/bin/python -I scripts/e22_demo_ambiguity.py --demo demo --out results/e22
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, lsq_linear

from agentic_tcad_demo.config import load_config
from agentic_tcad_demo.controller import Controller, Reviewer, RulePlanner
from agentic_tcad_demo.device import Q, CompactTFT, Condition, synthesize_measurement
from agentic_tcad_demo.evaluate import evaluate_fixed, summarize
from agentic_tcad_demo.log import ActionLog
from agentic_tcad_demo.metrics import extract_metrics
from agentic_tcad_demo.objective import objective, score_curve
from agentic_tcad_demo.parameters import ParameterGroup, ParameterSpec
from agentic_tcad_demo.runner import all_conditions, run_demo
from agentic_tcad_demo.tools import Workspace

ACTIVE = frozenset({"fixed_charge", "tail_states"})
NAMES = ("phi_SD", "phi_G", "mu0", "Q_f", "N_tail")
SIGMA_LOG = 0.03 / np.log(10.0)  # 3 % lognormal noise of the demo, in decades
FLOOR_A = 1e-11  # the demo's 10 pA measurement floor for log-current residuals
VG_WIDE = np.round(np.arange(-3.0, 3.0 + 1e-9, 0.02), 10)  # prediction grid for split devices


@dataclasses.dataclass(frozen=True)
class GeoCondition(Condition):
    """Transfer-curve condition that also carries the gate-oxide and channel thickness."""

    t_ox_nm: float = 7.0
    t_ch_nm: float = 3.5

    @property
    def label(self) -> str:
        return f"L={self.length_um:g}um,VD={self.vd:g}V,tox={self.t_ox_nm:g}nm,tch={self.t_ch_nm:g}nm"


class GeoSimulator:
    """The demo's CompactTFT with per-condition t_ox and t_ch.

    n_interface = C_it / C_ox scales with t_ox (C_it fixed). With free_cit=True its value at the
    reference t_ox is the calibration parameter "n_it" instead of the demo's fixed constant.
    """

    def __init__(self, base, mobility_gamma=0.0, dibl=0.0, free_cit=False, bisection_iters=60):
        self.base, self.gamma, self.dibl = base, mobility_gamma, dibl
        self.free_cit, self.iters = free_cit, bisection_iters

    def model(self, cond, params) -> CompactTFT:
        t_ox = getattr(cond, "t_ox_nm", self.base.t_ox_nm)
        t_ch = getattr(cond, "t_ch_nm", self.base.t_ch_nm)
        n_it_ref = params["n_it"] if self.free_cit else self.base.n_interface
        scale = t_ox / self.base.t_ox_nm
        const = dataclasses.replace(self.base, t_ox_nm=t_ox, t_ch_nm=t_ch, n_interface=n_it_ref * scale)
        return CompactTFT(const, self.gamma, self.dibl, self.iters)

    def simulate(self, params, active, cond, vg):
        p = {k: v for k, v in params.items() if k != "n_it"}
        return self.model(cond, params).simulate(p, active, cond, vg)

    def recovered(self) -> "GeoSimulator":
        return GeoSimulator(self.base, self.gamma, self.dibl, self.free_cit, 2 * self.iters)


def G(vd, length=1.0, t_ox=7.0, t_ch=3.5) -> GeoCondition:
    return GeoCondition(float(vd), float(length), float(t_ox), float(t_ch))


def synth_extra(truth: GeoSimulator, params, conds, vg, seed, sigma, floor) -> dict:
    """Noisy truth curves for split conditions; one seed stream per condition geometry."""
    out = {}
    for c in conds:
        entropy = [seed, int(round(1e3 * c.vd)), int(round(1e3 * c.length_um)),
                   int(round(10 * c.t_ox_nm)), int(round(10 * c.t_ch_nm))]
        rng = np.random.default_rng(np.random.SeedSequence(entropy))
        res = truth.simulate(params, ACTIVE, c, vg)
        if not res.converged:
            raise RuntimeError(f"truth device failed for {c.label}")
        noisy = res.id * np.exp(sigma * rng.standard_normal(vg.size))
        out[c] = np.abs(noisy + floor * rng.standard_normal(vg.size))
    return out


def metrics_wide(sim, params, active, cond) -> dict:
    res = sim.simulate(params, active, cond, VG_WIDE)
    m = extract_metrics(VG_WIDE, res.id)
    return {"vth_V": m.vth, "ss_V_per_dec": m.ss, "log10_ion": float(np.log10(m.ion)), "converged": res.converged}


def on_null_line(params, qf, c_ox) -> dict:
    """Move Q_f to `qf` and phi_G so that V_T = phi_G - chi - q Q_f / C_ox is unchanged."""
    p = dict(params)
    p["phi_G"] = params["phi_G"] + Q * (qf - params["Q_f"]) / c_ox
    p["Q_f"] = qf
    return p


def on_cit_line(params, n_tail, t_ch_nm, c_ox) -> dict:
    """Move N_tail to `n_tail` and n_it so that n = 1 + n_it + q N_tail t_ch / C_ox is unchanged."""
    p = dict(params)
    p["n_it"] = params["n_it"] - Q * (n_tail - params["N_tail"]) * t_ch_nm * 1e-7 / c_ox
    p["N_tail"] = n_tail
    return p


class LogLik:
    """Gaussian log-current likelihood of the demo noise model (points above the 10 pA floor)."""

    def __init__(self, sim, measured, conds, vg, active=ACTIVE):
        self.sim, self.meas, self.conds, self.vg, self.active = sim, measured, list(conds), vg, active
        self.masks = {c: measured[c] >= FLOOR_A for c in self.conds}

    def resid(self, params) -> np.ndarray | None:
        out = []
        for c in self.conds:
            res = self.sim.simulate(params, self.active, c, self.vg)
            if not res.converged:
                return None
            m = self.masks[c]
            out.append(np.log10(np.maximum(res.id[m], 1e-15)) - np.log10(self.meas[c][m]))
        return np.concatenate(out) / SIGMA_LOG

    def chi2(self, params) -> float:
        r = self.resid(params)
        return float("inf") if r is None else float(r @ r)

    def fit(self, start, free, specs, fixed=None, max_nfev=400):
        """Minimize chi2 over `free` (normalized coordinates) with the other parameters fixed."""
        base = {**start, **(fixed or {})}
        x0 = np.clip([specs[n].to_x(base[n]) for n in free], 1e-9, 1 - 1e-9)

        def fun(x):
            p = {**base, **{n: specs[n].from_x(xi) for n, xi in zip(free, x)}}
            r = self.resid(p)
            return np.full(self.n_rows, 1e3) if r is None else r

        self.n_rows = sum(int(m.sum()) for m in self.masks.values())
        sol = least_squares(fun, x0, bounds=(0.0, 1.0), diff_step=1e-6, xtol=1e-12, ftol=1e-12,
                            gtol=1e-12, max_nfev=max_nfev, method="trf")
        p = {**base, **{n: specs[n].from_x(xi) for n, xi in zip(free, sol.x)}}
        return p, float(sol.fun @ sol.fun)

    def fisher(self, params, names, h=1e-6) -> np.ndarray:
        """Fisher matrix in physical units (phi in eV, Q_f in 1e12 cm^-2, densities in log10)."""
        cols = []
        for n in names:
            p1, p2 = dict(params), dict(params)
            if n in ("mu0", "N_tail"):
                p1[n] = params[n] * 10 ** (-h)
                p2[n] = params[n] * 10 ** h
                step = 2 * h
            elif n == "Q_f":
                p1[n], p2[n] = params[n] - 1e12 * h, params[n] + 1e12 * h
                step = 2 * h
            else:
                p1[n], p2[n] = params[n] - h, params[n] + h
                step = 2 * h
            cols.append((self.resid(p2) - self.resid(p1)) / step)
        s = np.column_stack(cols)
        return s.T @ s


def demo_j(sim, measured, conds, vg, spec, params, active=ACTIVE) -> float:
    scores = []
    for c in conds:
        res = sim.simulate(params, active, c, vg)
        scores.append(score_curve(vg, res.id, measured[c], spec))
    return objective(scores)


class JointReachController(Controller):
    """The demo's controller with one extra, OR-ed eligibility test for "insensitive" components.

    The demo's test asks whether any single active parameter, moved across its full range, could
    remove at least half of the residual to first order. The extra test asks the same of all
    active parameters moved jointly within their bounds (bounded linear least squares on the
    OFAT columns of that residual's rows). It differs from the original only when one parameter
    shifts the residual of every calibration condition alike but cannot fit their pattern, for
    example an interface term facing two channel thicknesses.

    The share and persistence checks below repeat those of the demo's Controller._eligible
    (agentic-TCAD-calibration-public, MIT License, Copyright (c) 2026 Gyujun Jeong).
    """

    def _eligible(self, ev, params, active, rankings, exhausted):
        out = super()._eligible(ev, params, active, rankings, exhausted)
        have = {c["name"] for c in out}
        terms = ev.weighted_terms
        for comp in self.components.values():
            if (comp.trigger != "insensitive" or comp.name in active or comp.name in have
                    or comp.name in self._rejected_components):
                continue
            share = terms[comp.addresses] / max(ev.j, 1e-30)
            if share < self.s.activation_min_share:
                continue
            hist = self._share_history[comp.addresses][-self.s.persistence_cycles:]
            persisted = len(hist) >= self.s.persistence_cycles and min(hist) >= self.s.activation_min_share
            if not (persisted or exhausted):
                continue
            complete = [rk for rk in rankings.values() if all(p["converged"] for p in rk.probes)]
            if not complete or len(complete) < len(rankings):
                continue
            rows = self.ws.metric_rows(comp.addresses)
            r_rows = complete[0].r[rows]
            r_norm = float(np.linalg.norm(r_rows))
            if r_norm == 0.0:
                continue
            cols, xs, seen = [], [], set()
            for rk in complete:
                for j, n in enumerate(rk.names):
                    if n not in seen:
                        seen.add(n)
                        cols.append(rk.jacobian[rows, j])
                        xs.append(rk.x[j])
            s_mat, x = np.column_stack(cols), np.array(xs)
            sol = lsq_linear(s_mat, -r_rows, bounds=(-x, 1.0 - x), method="bvls")
            remaining = float(np.linalg.norm(r_rows + s_mat @ sol.x)) / r_norm
            if remaining > 1.0 - self.s.insensitivity_ratio:
                out.append({"name": comp.name, "addresses": comp.addresses, "share_of_J": share,
                            "reason": f"active parameters moved jointly within bounds leave {remaining:.2f} "
                                      f"of the {comp.addresses} residual (first order)"})
        return out


def run_controller(cfg, sim, measured, calib, initial, specs=None, groups=None, settings=None, cls=Controller):
    log = ActionLog(None)
    ws = Workspace(sim, measured, cfg.vg, specs or cfg.specs, groups or cfg.groups, cfg.objective,
                   calib, log, cfg.ofat_step)
    ctl = cls(ws, cfg.components, RulePlanner(), Reviewer(cfg.review_mode, cfg.approved_components),
              settings or cfg.controller, log)
    out = ctl.run(initial)
    out.log = log
    return out, ws


def heldout_max_dvth(ws, params, active, conds) -> float:
    rows = evaluate_fixed(ws, params, active, conds)["rows"]
    return float(summarize(rows)["max_abs_dVth_mV"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--profile-points", type=int, default=29)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    res: dict = {"demo_commit": "ccaa678", "assumptions": [
        "C_it and Q_f do not change with t_ox or t_ch; n_interface = C_it/C_ox scales with t_ox",
        "split-device data come from the demo's synthetic truth (gate-dependent mobility, DIBL) "
        "with the demo's noise model; one seed stream per condition geometry",
        "thickness-test predictions use a V_G grid of -3..3 V (20 mV) so that all metrics are defined"]}

    # ---- 0. reproduce the default demo run --------------------------------------------------
    cfg = load_config(args.demo / "configs" / "synthetic_demo.yaml")
    run = run_demo(cfg, None)
    final = dict(run["result"].params)
    assert run["result"].active == ACTIVE
    base, vg, spec = cfg.constants, cfg.vg, cfg.objective
    c_ox = base.c_ox
    truth_p = dict(cfg.truth_params)
    res["default_run"] = {"J_initial": run["summary"]["J_initial"], "J_final": run["summary"]["J_final"],
                          "accepted": len(run["result"].accepted), "rejected": run["summary"]["rejected_trials"],
                          "final_params": final, "truth_params": truth_p,
                          "heldout_max_abs_dVth_mV": {k: v["final"]["max_abs_dVth_mV"]
                                                      for k, v in run["summary"]["evaluation"].items()}}

    model = GeoSimulator(base)
    truth = GeoSimulator(base, cfg.truth.mobility_gamma, cfg.truth.dibl_v_per_v_um)
    conds0 = all_conditions(cfg)
    meas = {G(c.vd, c.length_um): run["measured"][c] for c in conds0}
    # consistency checks: same model and same data as the demo
    for c in conds0:
        a = cfg.simulator.simulate(final, ACTIVE, c, vg).id
        b = model.simulate(final, ACTIVE, G(c.vd, c.length_um), vg).id
        assert np.array_equal(a, b), "GeoSimulator differs from the demo model"
    redo = synthesize_measurement(cfg.truth, cfg.truth_params, cfg.truth_active, conds0, vg, cfg.seed,
                                  cfg.noise["lognormal_sigma"], cfg.noise["floor_noise_A"])
    assert all(np.array_equal(redo[c], run["measured"][c]) for c in conds0)
    calib2 = [G(0.05), G(1.0)]
    heldout = {"drain": [G(0.1), G(0.2), G(0.5)], "length": [G(0.05, 2.0), G(1.0, 2.0)]}

    # extra split data (noisy truth) on the demo grid, used for calibration
    tox_curve = G(0.05, 1.0, 14.0)
    tch_curve = G(0.05, 1.0, 7.0, 7.0)
    meas.update(synth_extra(truth, truth_p, [tox_curve, tch_curve], vg, cfg.seed,
                            cfg.noise["lognormal_sigma"], cfg.noise["floor_noise_A"]))
    calib3 = calib2 + [tox_curve]

    tests = {"VD 0.1 V": G(0.1), "VD 0.2 V": G(0.2), "VD 0.5 V": G(0.5),
             "L 2 um, VD 0.05 V": G(0.05, 2.0), "L 2 um, VD 1 V": G(1.0, 2.0),
             "t_ox 5 nm": G(0.05, 1.0, 5.0), "t_ox 10 nm": G(0.05, 1.0, 10.0),
             "t_ox 14 nm": G(0.05, 1.0, 14.0), "t_ox 20 nm": G(0.05, 1.0, 20.0),
             "t_ch 2 nm": G(0.05, 1.0, 7.0, 2.0), "t_ch 7 nm": G(0.05, 1.0, 7.0, 7.0)}
    truth_m = {k: metrics_wide(truth, truth_p, ACTIVE, c) for k, c in tests.items()}

    # ---- A. exact null direction phi_G - Q_f --------------------------------------------------
    spg, sqf = cfg.specs["phi_G"], cfg.specs["Q_f"]
    q_lo = max(sqf.lower, final["Q_f"] + (spg.lower - final["phi_G"]) * c_ox / Q)
    q_hi = min(sqf.upper, final["Q_f"] + (spg.upper - final["phi_G"]) * c_ox / Q)
    valley = [on_null_line(final, q, c_ox) for q in np.linspace(q_lo, q_hi, 11)]
    valley_truthq = on_null_line(final, truth_p["Q_f"], c_ox)
    j_valley = [demo_j(model, meas, calib2, vg, spec, p) for p in valley]
    dlog = 0.0
    for p in valley:
        for c in calib2 + heldout["drain"] + heldout["length"]:
            a = model.simulate(final, ACTIVE, c, vg).id
            b = model.simulate(p, ACTIVE, c, vg).id
            dlog = max(dlog, float(np.max(np.abs(np.log10(b) - np.log10(a)))))
    # SVD of the demo's own weighted residual Jacobian (Eq. 2) at the final point, all 5 parameters
    ws_a = Workspace(model, meas, vg, cfg.specs, cfg.groups, spec, calib2, ActionLog(None), cfg.ofat_step)
    r0, w0 = ws_a.residual(final, ACTIVE)
    h = 1e-4
    cols = []
    for n in NAMES:
        x0 = cfg.specs[n].to_x(final[n])
        pp, pm = dict(final), dict(final)
        pp[n], pm[n] = cfg.specs[n].from_x(x0 + h), cfg.specs[n].from_x(x0 - h)
        cols.append((ws_a.residual(pp, ACTIVE)[0] - ws_a.residual(pm, ACTIVE)[0]) / (2 * h))
    a_mat = np.column_stack(cols) * np.sqrt(w0)[:, None]
    _, sv, vt = np.linalg.svd(a_mat, full_matrices=False)
    v_null = vt[-1]
    i_g, i_q = NAMES.index("phi_G"), NAMES.index("Q_f")
    dphig = v_null[i_g] * (spg.upper - spg.lower)
    dqf = v_null[i_q] * (sqf.upper - sqf.lower)
    res["A_null_phiG_Qf"] = {
        "C_ox_F_per_cm2": c_ox, "dphiG_per_1e12_Qf_mV": 1e3 * Q * 1e12 / c_ox,
        "valley_Q_f_cm2": [q_lo, q_hi], "valley_phi_G_eV": [valley[0]["phi_G"], valley[-1]["phi_G"]],
        "phi_G_at_truth_Q_f": valley_truthq["phi_G"],
        "J_along_valley_minmax": [min(j_valley), max(j_valley)],
        "max_abs_dlog10I_calib_and_heldout_dec": dlog,
        "eq2_jacobian_singular_values_rel": (sv / sv[0]).tolist(),
        "eq2_null_vector_normalized_x": dict(zip(NAMES, v_null.tolist())),
        "eq2_null_vector_dphiG_per_1e12_Qf_mV": 1e3 * dphig / (dqf / 1e12)}

    # ---- B. prediction spread across the valley ------------------------------------------------
    spread = {}
    for k, c in tests.items():
        vals = [metrics_wide(model, p, ACTIVE, c) for p in valley]
        vth = [v["vth_V"] for v in vals]
        ss = [v["ss_V_per_dec"] for v in vals]
        fin = metrics_wide(model, final, ACTIVE, c)
        tq = metrics_wide(model, valley_truthq, ACTIVE, c)
        spread[k] = {"vth_spread_mV": 1e3 * (max(vth) - min(vth)), "ss_spread_mV_dec": 1e3 * (max(ss) - min(ss)),
                     "final_minus_truth_vth_mV": 1e3 * (fin["vth_V"] - truth_m[k]["vth_V"]),
                     "valley_at_truth_Qf_minus_truth_vth_mV": 1e3 * (tq["vth_V"] - truth_m[k]["vth_V"]),
                     "final_minus_truth_ss_mV_dec": 1e3 * (fin["ss_V_per_dec"] - truth_m[k]["ss_V_per_dec"])}
    res["B_valley_prediction_spread"] = spread

    # ---- C. profile likelihood of Q_f ----------------------------------------------------------
    free4 = ("phi_SD", "phi_G", "mu0", "N_tail")
    prof = {}
    for tag, conds in (("demo set (L 1 um, VD 0.05/1 V)", calib2), ("+ t_ox 14 nm, VD 0.05 V", calib3)):
        ll = LogLik(model, meas, conds, vg)
        best, chi_best = ll.fit(final, NAMES, cfg.specs, max_nfev=800)
        if tag.startswith("+"):
            # also start from the truth-Q_f point of the null line, keep the better optimum
            b2, c2 = ll.fit(valley_truthq, NAMES, cfg.specs, max_nfev=800)
            if c2 < chi_best:
                best, chi_best = b2, c2
        n_rows = ll.n_rows
        grid = np.linspace(sqf.lower, sqf.upper, args.profile_points)
        extra = np.concatenate([np.linspace(best["Q_f"] - 2e12, best["Q_f"] + 2e12, 17),
                                np.linspace(best["Q_f"] - 3e11, best["Q_f"] + 3e11, 25),
                                np.linspace(q_lo - 4e11, q_lo + 4e11, 9)])
        grid = np.unique(np.concatenate([grid, extra[(extra >= sqf.lower) & (extra <= sqf.upper)]]))
        rows = []
        for q in grid:
            start = on_null_line(best, q, c_ox)
            start["phi_G"] = float(np.clip(start["phi_G"], spg.lower, spg.upper))
            p, chi = ll.fit(start, free4, cfg.specs, fixed={"Q_f": q})
            rows.append({"Q_f": q, "chi2": chi, "phi_G": p["phi_G"], "N_tail": p["N_tail"],
                         "J": demo_j(model, meas, conds, vg, spec, p)})
        chi_best = min(chi_best, min(r["chi2"] for r in rows))
        sigma_eff_dec = SIGMA_LOG * np.sqrt(chi_best / n_rows)
        for r in rows:
            r["dchi2_eff"] = (r["chi2"] - chi_best) * (SIGMA_LOG / sigma_eff_dec) ** 2
        fis = ll.fisher(best, NAMES)
        eig = np.linalg.eigvalsh(fis)
        cov = None
        try:
            cov_eff = np.linalg.inv(fis) * (sigma_eff_dec / SIGMA_LOG) ** 2
            cov = {n: float(np.sqrt(cov_eff[i, i])) for i, n in enumerate(NAMES)}
        except np.linalg.LinAlgError:
            pass
        ok = [r["Q_f"] for r in rows if r["dchi2_eff"] <= 3.84]
        prof[tag] = {"best": best, "chi2_best": chi_best, "n_rows": n_rows, "sigma_eff_dec": sigma_eff_dec,
                     "fisher_eig_ratio_min_max": float(eig[0] / eig[-1]),
                     "linear_1sigma_eff (Q_f in 1e12 cm^-2, densities in dec)": cov,
                     "Q_f_range_dchi2eff_le_3.84": [min(ok), max(ok)] if ok else None, "rows": rows}
    res["C_profile_Q_f"] = prof

    # ---- D. path dependence of the demo controller ---------------------------------------------
    runs = []
    for calib_tag, calib in (("demo set", calib2), ("+ t_ox 14 nm curve", calib3)):
        for delta in (0.3, 0.4, 0.5, 0.6, 0.7):
            for phig0 in (4.45, 4.55, 4.65):
                init = dict(cfg.initial)
                init["phi_G"] = phig0
                settings = dataclasses.replace(cfg.controller, trust_radius=delta)
                out, ws = run_controller(cfg, model, meas, calib, init, settings=settings)
                p = out.params
                act = out.active
                row = {"calibration": calib_tag, "trust_radius": delta, "phi_G_initial": phig0,
                       "J_final": out.j_history[-1], "accepted": len(out.accepted), "activations": out.activations,
                       "phi_G": p["phi_G"], "Q_f": p["Q_f"], "N_tail": p["N_tail"], "mu0": p["mu0"],
                       "phi_SD": p["phi_SD"],
                       "heldout_drain_max_dVth_mV": heldout_max_dvth(ws, p, act, heldout["drain"]),
                       "heldout_length_max_dVth_mV": heldout_max_dvth(ws, p, act, heldout["length"])}
                for k in ("t_ox 10 nm", "t_ox 14 nm", "t_ox 20 nm"):
                    m = metrics_wide(model, p, act, tests[k])
                    row[f"{k} dVth vs truth (mV)"] = 1e3 * (m["vth_V"] - truth_m[k]["vth_V"])
                runs.append(row)
    res["D_controller_runs"] = runs

    # ---- E. interface traps (C_it) vs tail states ----------------------------------------------
    model_cit = GeoSimulator(base, free_cit=True)
    specs_cit = dict(cfg.specs)
    specs_cit["n_it"] = ParameterSpec("n_it", 0.0, 1.0, "linear", "C_it/C_ox at t_ox = 7 nm")
    groups_cit = dict(cfg.groups)
    groups_cit["interface"] = ParameterGroup("interface", ("n_it",))
    fin_cit = {**final, "n_it": base.n_interface}
    t_ch0 = base.t_ch_nm
    m_tot = base.n_interface + Q * final["N_tail"] * t_ch0 * 1e-7 / c_ox
    n_max = (m_tot - 0.0) * c_ox / (Q * t_ch0 * 1e-7)
    n_lo = specs_cit["N_tail"].lower
    cit_valley = [on_cit_line(fin_cit, n, t_ch0, c_ox) for n in np.logspace(np.log10(n_lo), np.log10(n_max), 9)]
    cit_spread = {}
    for k, c in tests.items():
        vals = [metrics_wide(model_cit, p, ACTIVE, c) for p in cit_valley]
        ss = [v["ss_V_per_dec"] for v in vals]
        vth = [v["vth_V"] for v in vals]
        cit_spread[k] = {"ss_spread_mV_dec": 1e3 * (max(ss) - min(ss)), "vth_spread_mV": 1e3 * (max(vth) - min(vth))}
    jv = [demo_j(model_cit, meas, calib2, vg, spec, p) for p in cit_valley]
    # controller model selection: activation test x update grouping x calibration set
    groups_joint = dict(groups_cit)
    groups_joint["tail_states"] = ParameterGroup("tail_states", ("N_tail", "n_it"), "tail_states")
    sel = []
    for trig, cls in (("demo test", Controller), ("demo test OR joint reach", JointReachController)):
        for grouping, grp in (("separate", groups_cit), ("tail_states group = [N_tail, n_it]", groups_joint)):
            for tag, calib in (("demo set", calib2), ("+ t_ox 14 nm curve", calib3),
                               ("+ t_ch 7 nm curve", calib2 + [tch_curve]),
                               ("+ t_ch 7 nm + t_ox 14 nm curves", calib3 + [tch_curve])):
                init = {**cfg.initial, "n_it": base.n_interface}
                out, ws = run_controller(cfg, model_cit, meas, calib, init, specs=specs_cit, groups=grp, cls=cls)
                p = out.params
                k_ch = metrics_wide(model_cit, p, out.active, tests["t_ch 7 nm"])
                k_ox = metrics_wide(model_cit, p, out.active, tests["t_ox 20 nm"])
                acts = [e for e in out.log.events if e["event"] == "activation"]
                sel.append({"activation_test": trig, "grouping": grouping, "calibration": tag,
                            "J_final": out.j_history[-1], "accepted": len(out.accepted),
                            "activations": out.activations,
                            "activation_attempts": [{"component": e["component"], "kept": e["kept"]} for e in acts],
                            "n_it": p["n_it"], "N_tail": p["N_tail"] if "tail_states" in out.active else None,
                            "Q_f": p["Q_f"] if "fixed_charge" in out.active else None, "phi_G": p["phi_G"],
                            "heldout_drain_max_dVth_mV": heldout_max_dvth(ws, p, out.active, heldout["drain"]),
                            "heldout_length_max_dVth_mV": heldout_max_dvth(ws, p, out.active, heldout["length"]),
                            "t_ch 7 nm dSS vs truth (mV/dec)":
                                1e3 * (k_ch["ss_V_per_dec"] - truth_m["t_ch 7 nm"]["ss_V_per_dec"]),
                            "t_ox 20 nm dVth vs truth (mV)": 1e3 * (k_ox["vth_V"] - truth_m["t_ox 20 nm"]["vth_V"]),
                            "stop_reason": out.stop_reason})
    # profile likelihood of N_tail with the interface term free (both components active)
    names6 = NAMES + ("n_it",)
    free5 = ("phi_SD", "phi_G", "mu0", "Q_f", "n_it")
    prof_n = {}
    truth_cit = {**truth_p, "n_it": base.n_interface}
    for tag, conds in (("demo set", calib2), ("+ t_ox 14 nm curve", calib3),
                       ("+ t_ch 7 nm curve", calib2 + [tch_curve])):
        ll = LogLik(model_cit, meas, conds, vg)
        best, chi_best = ll.fit(fin_cit, names6, specs_cit, max_nfev=1200)
        b2, c2 = ll.fit({**best, "N_tail": truth_cit["N_tail"], "n_it": truth_cit["n_it"]}, names6, specs_cit,
                        max_nfev=1200)
        if c2 < chi_best:
            best, chi_best = b2, c2
        lg = np.log10(best["N_tail"])
        grid = np.unique(np.concatenate([np.linspace(16.0, 20.0, 17), np.linspace(lg - 0.15, lg + 0.15, 13)]))
        grid = grid[(grid >= 16.0) & (grid <= 20.0)]
        rows = []
        for g in grid:
            nt = 10 ** g
            start = on_cit_line(best, nt, t_ch0, c_ox)
            start["n_it"] = float(np.clip(start["n_it"], 0.0, 1.0))
            p, chi = ll.fit(start, free5, specs_cit, fixed={"N_tail": nt})
            rows.append({"log10_N_tail": g, "chi2": chi, "n_it": p["n_it"]})
        chi_best = min(chi_best, min(r["chi2"] for r in rows))
        sig_eff = SIGMA_LOG * np.sqrt(chi_best / ll.n_rows)
        for r in rows:
            r["dchi2_eff"] = (r["chi2"] - chi_best) * (SIGMA_LOG / sig_eff) ** 2
        # Q_f is left out: with no t_ox split it is exactly degenerate with phi_G
        names_noq = tuple(n for n in names6 if n != "Q_f")
        eig = np.linalg.eigvalsh(ll.fisher(best, names_noq))
        ok = [r["log10_N_tail"] for r in rows if r["dchi2_eff"] <= 3.84]
        prof_n[tag] = {"best_N_tail": best["N_tail"], "best_n_it": best["n_it"], "best_Q_f": best["Q_f"],
                       "best_phi_G": best["phi_G"], "sigma_eff_dec": sig_eff,
                       "fisher_eig_ratio_min_max_without_Q_f": float(eig[0] / eig[-1]),
                       "log10_N_tail_range_dchi2eff_le_3.84": [min(ok), max(ok)] if ok else None, "rows": rows}

    # the OR-ed joint-reach test leaves the default demo run unchanged
    out_d, _ = run_controller(cfg, model, meas, calib2, cfg.initial, cls=JointReachController)
    joint = {"case": "default demo (n_interface fixed)", "J_final": out_d.j_history[-1],
             "activations": out_d.activations, "params": out_d.params,
             "identical_to_demo": out_d.params == final}
    res["E_cit_vs_tail"] = {"identified_sum_n_it_plus_tail": m_tot, "valley_N_tail": [n_lo, n_max],
                            "valley_n_it": [cit_valley[0]["n_it"], cit_valley[-1]["n_it"]],
                            "J_along_valley_minmax": [min(jv), max(jv)], "prediction_spread": cit_spread,
                            "controller_model_selection": sel, "profile_N_tail": prof_n,
                            "joint_reach_controller": joint}

    res["truth_metrics_tests"] = truth_m
    (args.out / "summary.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(f"done in {time.time() - t0:.0f} s -> {args.out / 'summary.json'}")


if __name__ == "__main__":
    main()
