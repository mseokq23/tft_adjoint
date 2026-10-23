"""Solver-consistent parameter sensitivities of the drain current.

At a converged bias point the Newton residual satisfies F(u, theta) = 0.
DEVSIM's get_matrix_and_rhs() returns the Newton Jacobian J = dF/du and
R = -F (the Newton right-hand side). Differentiating F(u(theta), theta) = 0:

    J du/dtheta = dR/dtheta                     (forward / direct sensitivity)
    dI/dtheta   = e_I^T du/dtheta = lam^T dR/dtheta,   J^T lam = e_I   (adjoint)

where I = u[e_I] is the circuit unknown "V1.I" (drain current per unit depth).

dR/dtheta is obtained by central differencing ONE residual assembly at the
fixed converged solution (no nonlinear re-solve), so the gradient inherits the
solver's discretisation exactly. It is therefore "solver-consistent" but not
symbolic; its error is O(h^2) plus round-off and is reported by e0_gradient_check.py.
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import devsim as ds

from .params import RES_STEP, RESOLVE_STEP, step_for


class ScaledLU:
    """LU of the row/column-equilibrated Jacobian with iterative refinement.

    The Newton matrix mixes potentials (V) and electron densities (up to ~1e20
    cm^-3); a plain sparse LU loses ~5-6 digits in the output component. With
    As = Dr J Dc factorised and two refinement sweeps, adjoint and forward
    results agree normwise to 1e-11 - 1e-9 at and above threshold and to ~2e-6 in
    deep subthreshold, where small components (e.g. dID/dWF_SD) differ by up to
    ~4e-3 relative (results/e0, gradients.csv).
    """

    def __init__(self, J, refine: int = 2):
        self.J = J.tocsc()
        self.JT = self.J.T.tocsc()
        A = J.tocsr()
        r = np.asarray(abs(A).max(axis=1).todense()).ravel()
        r[r == 0] = 1.0
        self.dr = 1.0 / r
        A1 = sp.diags(self.dr) @ A
        c = np.asarray(abs(A1).max(axis=0).todense()).ravel()
        c[c == 0] = 1.0
        self.dc = 1.0 / c
        self.lu = spla.splu((A1 @ sp.diags(self.dc)).tocsc())
        self.refine = refine

    def solve(self, b, trans: str = "N"):
        b = np.asarray(b, float)
        if trans == "N":
            x = self.dc * self.lu.solve(self.dr * b)
            for _ in range(self.refine):
                x = x + self.dc * self.lu.solve(self.dr * (b - self.J @ x))
        else:
            x = self.dr * self.lu.solve(self.dc * b, trans="T")
            for _ in range(self.refine):
                x = x + self.dr * self.lu.solve(self.dc * (b - self.JT @ x), trans="T")
        return x


class SensitivityEngine:
    def __init__(self, tft, central: bool = True):
        self.tft = tft
        self.central = central
        self.iI = ds.get_circuit_equation_number(node="V1.I")
        self.timing = {"assemble": 0.0, "factor": 0.0, "solve": 0.0}
        self.counts = {"assemble": 0, "factor": 0, "solve": 0}

    # ------------------------------------------------------------- assembly
    def assemble(self):
        t0 = time.perf_counter()
        m = ds.get_matrix_and_rhs(format="csc")["static"]
        n = len(m["rhs"])
        J = sp.csc_matrix((np.asarray(m["av"], float), np.asarray(m["ai"]), np.asarray(m["ap"])),
                          shape=(n, n))
        R = np.asarray(m["rhs"], float)
        self.timing["assemble"] += time.perf_counter() - t0
        self.counts["assemble"] += 1
        return J, R

    def residual(self):
        return self.assemble()[1]

    def dR(self, name: str, h: float | None = None, R0=None):
        """dR/dtheta by differencing one residual assembly (no re-solve)."""
        tft = self.tft
        v0 = tft.get_theta(name)
        h = step_for(name, RES_STEP) if h is None else h
        try:
            tft.set_theta(name, v0 + h)
            rp = self.residual()
            if self.central:
                tft.set_theta(name, v0 - h)
                rm = self.residual()
                d = (rp - rm) / (2.0 * h)
            else:
                r0 = self.residual() if R0 is None else R0
                d = (rp - r0) / h
        finally:
            tft.set_theta(name, v0)
        return d

    def factorize(self):
        J, R0 = self.assemble()
        t0 = time.perf_counter()
        lu = ScaledLU(J)
        self.timing["factor"] += time.perf_counter() - t0
        self.counts["factor"] += 1
        return lu, R0, J.shape[0]

    # ------------------------------------------------------------ gradients
    def current_gradient(self, names, mode: str = "adjoint", return_fields: bool = False):
        """dID/dtheta (A per unit of each theta) at the present bias point.

        mode = "adjoint": one transposed solve serves all parameters.
        mode = "forward": one solve per parameter (also returns du/dtheta fields).
        """
        lu, R0, n = self.factorize()
        dRs = [self.dR(nm, R0=R0) for nm in names]
        scale = self.tft.isign * self.tft.width_cm
        t0 = time.perf_counter()
        fields = None
        if mode == "adjoint":
            e = np.zeros(n)
            e[self.iI] = 1.0
            lam = lu.solve(e, trans="T")
            g = np.array([lam @ d for d in dRs])
            self.counts["solve"] += 1
            fields = lam
        elif mode == "forward":
            sols = [lu.solve(d) for d in dRs]
            g = np.array([s[self.iI] for s in sols])
            self.counts["solve"] += len(dRs)
            fields = sols
        else:
            raise ValueError(mode)
        self.timing["solve"] += time.perf_counter() - t0
        g = scale * g
        return (g, fields) if return_fields else g

    def functional_gradient(self, c, names, lu=None, R0=None):
        """d(c^T u)/dtheta for a linear output functional c^T u of the solution (one adjoint solve).

        Used for quantities other than the drain current, e.g. the gate charge
        (tftsim.cv.GateCharge). Returns (gradient, lam)."""
        if lu is None:
            lu, R0, n = self.factorize()
        t0 = time.perf_counter()
        lam = lu.solve(np.asarray(c, float), trans="T")
        self.counts["solve"] += 1
        self.timing["solve"] += time.perf_counter() - t0
        g = np.array([lam @ self.dR(nm, R0=R0) for nm in names])
        return g, lam

    def resolve_fd(self, name: str, h: float | None = None):
        """Reference: central difference with two full nonlinear re-solves."""
        tft = self.tft
        v0 = tft.get_theta(name)
        h = step_for(name, RESOLVE_STEP) if h is None else h
        snap = tft.snapshot()
        try:
            tft.set_theta(name, v0 + h)
            tft.solve()
            ip = tft.drain_current()
            tft.restore(snap)
            tft.set_theta(name, v0 - h)
            tft.solve()
            im = tft.drain_current()
        finally:
            tft.set_theta(name, v0)
            tft.restore(snap)
        return (ip - im) / (2.0 * h)

    def report(self) -> dict:
        return {"timing_s": {k: round(v, 4) for k, v in self.timing.items()}, "counts": dict(self.counts)}
