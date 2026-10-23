"""Transfer-curve model with warm-started re-evaluation, for calibration (E2).

CurveModel keeps one converged solution per (VD, VG) point. A new parameter
vector is reached by a short homotopy at the first bias point; every other
point is re-solved starting from its cached solution (2-4 Newton steps), with a
fall-back ramp from the neighbouring point if that fails. Gradients come from
the adjoint engine at each point.
"""
from __future__ import annotations

import time

import numpy as np

from .device import TFT, ConvergenceError
from .metrics import I_FLOOR, LN10
from .params import Geometry
from .sensitivity import SensitivityEngine


class CurveModel:
    def __init__(self, base_params, names, geom: Geometry = Geometry(), vds=(0.05, 1.0), vgs=None, **tft_kw):
        self.names = list(names)
        self.vds = tuple(vds)
        self.vgs = np.asarray(vgs, float)
        self.t = TFT(base_params, geom, **tft_kw).build()
        self.eng = SensitivityEngine(self.t, central=False)
        self.cache = {}
        self.stats = {"n_eval": 0, "n_grad_eval": 0, "t_eval": 0.0, "t_grad": 0.0}
        for vd in self.vds:
            self.t.goto(vd=vd)
            self.t.goto(vg=self.vgs[0])
            for vg in self.vgs:
                self.t.goto(vg=vg)
                self.cache[(vd, vg)] = self.t.snapshot()
        self.x_cache = self.theta()

    def theta(self):
        return np.array([self.t.get_theta(n) for n in self.names])

    def _set(self, x):
        for n, v in zip(self.names, x):
            self.t.set_theta(n, float(v))

    def _move_theta(self, x_new):
        key0 = (self.vds[0], self.vgs[0])
        x_old = self.x_cache.copy()
        self.t.restore(self.cache[key0])
        done, step = 0.0, 1.0
        while done < 1.0:
            f = min(1.0, done + step)
            snap = self.t.snapshot()
            self._set(x_old + f * (x_new - x_old))
            try:
                self.t.solve()
                done = f
                step = min(1.0, 2.0 * step)
            except ConvergenceError:
                self.t.restore(snap)
                step *= 0.5
                if step < 1.0 / 1024:
                    self._set(x_old)
                    raise

    def evaluate(self, x, grads: bool = False):
        """Return ({vd: ID array}, {vd: dID/dtheta array or None})."""
        t0 = time.perf_counter()
        x = np.asarray(x, float)
        self._move_theta(x)
        I, G, new = {}, {}, {}
        first = True
        for vd in self.vds:
            ids, gs = [], []
            for vg in self.vgs:
                key = (vd, vg)
                if not first:
                    prev = self.t.snapshot()
                    self.t.restore(self.cache[key])
                    try:
                        self.t.solve()
                    except ConvergenceError:
                        self.t.restore(prev)
                        self.t.goto(vd=vd)
                        self.t.goto(vg=vg)
                first = False
                new[key] = self.t.snapshot()
                i = self.t.drain_current()
                ids.append(i)
                if grads:
                    gs.append(self.eng.current_gradient(self.names, "adjoint") if i >= 0.1 * I_FLOOR
                              else np.zeros(len(self.names)))
            I[vd] = np.array(ids)
            G[vd] = np.array(gs) if grads else None
        self.cache.update(new)
        self.x_cache = x.copy()
        dt = time.perf_counter() - t0
        self.stats["n_eval"] += 1
        self.stats["t_eval"] += dt
        if grads:
            self.stats["n_grad_eval"] += 1
            self.stats["t_grad"] += dt
        return I, G


def log_residual_problem(model: CurveModel, meas: dict, sigma: float, theta_ref, scale, history=None):
    """Residuals r = (log10 I_sim - log10 I_meas)/sigma on points with I_meas >= 10 pA,
    in scaled coordinates z = (theta - theta_ref)/scale."""
    masks = {vd: meas[vd] >= I_FLOOR for vd in model.vds}
    theta_ref = np.asarray(theta_ref, float)
    scale = np.asarray(scale, float)

    def _res(I):
        return np.concatenate([(np.log10(np.maximum(I[vd][masks[vd]], 1e-30)) - np.log10(meas[vd][masks[vd]])) / sigma
                               for vd in model.vds])

    def fun(z):
        I, _ = model.evaluate(theta_ref + np.asarray(z) * scale, grads=False)
        r = _res(I)
        if history is not None:
            history.append({"kind": "f", "z": np.asarray(z).tolist(), "cost": 0.5 * float(r @ r),
                            "n_solves": model.t.n_solves, "t": model.stats["t_eval"]})
        return r

    def jac(z):
        I, G = model.evaluate(theta_ref + np.asarray(z) * scale, grads=True)
        rows = [G[vd][masks[vd]] / (I[vd][masks[vd], None] * LN10) / sigma for vd in model.vds]
        Jz = np.vstack(rows) * scale
        if history is not None:
            r = _res(I)
            history.append({"kind": "j", "z": np.asarray(z).tolist(), "cost": 0.5 * float(r @ r),
                            "n_solves": model.t.n_solves, "t": model.stats["t_eval"]})
        return Jz

    return fun, jac, masks
