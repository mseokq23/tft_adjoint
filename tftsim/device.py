"""2-D bottom-gate oxide-semiconductor TFT in DEVSIM.

Model class (kept identical to the calibration descriptors of arXiv:2609.12184):
* electron-only drift-diffusion (Scharfetter-Gummel), Boltzmann statistics
* constant mobility mu0
* effective donor-like background ND,eff (fully ionised)
* donor-like Gaussian + acceptor-like exponential-tail DOS, discretised into
  trap levels on a fixed energy grid; trap occupancy f = n / (n + n1)
  (quasi-equilibrium with the conduction band, valid for DC sweeps)
* fixed sheet charge Qf at the gate-dielectric/channel interface
* work-function contacts: gate (Dirichlet on HfO2, psi = VG - (WFg - chi)),
  S/D (Dirichlet psi = V - (WFsd - chi), n = Nc exp(-(WFsd - chi)/Vt))

Potential reference: psi is shifted so that Ec = -q psi in the channel, i.e.
n = Nc exp((EFn + psi)/Vt); EF = 0 at zero bias.

The drain is driven by a circuit voltage source V1, so the drain current is a
solution unknown ("V1.I"). This makes terminal-current sensitivities a single
row of the Newton system (see sensitivity.py).
"""
from __future__ import annotations

import contextlib
import io
import math
import time

import numpy as np
import devsim as ds
from devsim.python_packages.model_create import (
    CreateContactNodeModel,
    CreateContinuousInterfaceModel,
    CreateEdgeModel,
    CreateEdgeModelDerivatives,
    CreateNodeModel,
    CreateNodeModelDerivative,
)
from devsim.python_packages.simple_dd import CreateBernoulli

from .params import (
    EPS0, Q, Geometry, Material, DeviceParams, THETA_DEFS,
    default_dos_edges, theta_value, trap_levels, with_theta,
)

NM = 1e-7   # cm
UM = 1e-4   # cm


class ConvergenceError(RuntimeError):
    pass


class TFT:
    """Bottom-gate TFT: build, bias, read currents, change parameters."""

    def __init__(self, params: DeviceParams, geom: Geometry = Geometry(), mat: Material = Material(),
                 *, dos_edges=None, mesh_scale: float = 1.0, solver: str = "umfpack",
                 name: str = "tft", quiet: bool = True, rel_tol: float = 1e-10):
        self.p = params
        self.geom = geom
        self.mat = mat
        self.edges = default_dos_edges() if dos_edges is None else np.asarray(dos_edges, float)
        self.K = len(self.edges) - 1 if params.dos_on else 0
        self.mesh_scale = mesh_scale
        self.solver = solver
        self.dev = name
        self.quiet = quiet
        self.vg = 0.0
        self.vd = 0.0
        self.isign = 1.0                      # set in build() so that ID > 0 for VD > 0
        self.tol = dict(absolute_error=1e30, relative_error=rel_tol, maximum_iterations=100)
        self.tol_p = dict(absolute_error=1e-12, relative_error=1e-12, maximum_iterations=80)
        self.n_solves = 0
        self.t_solve = 0.0
        self.n_nodes = 0

    # ------------------------------------------------------------------ build
    def build(self):
        ds.reset_devsim()
        self._set_solver()
        self._mesh()
        self._set_parameters()
        self._potential_only()
        self._solve_raw(self.tol_p)
        self._drift_diffusion()
        try:
            self._solve_raw(self.tol)
        except ds.error:
            # same relaxed retry as solve(): the relative error can stall just above 1e-10
            self._solve_raw(dict(self.tol, relative_error=max(1e-8, self.tol["relative_error"]), maximum_iterations=200))
            self.n_relaxed = getattr(self, "n_relaxed", 0) + 1
        self._calibrate_sign()
        return self

    @property
    def width_cm(self) -> float:
        return self.geom.W_um * UM

    def _set_solver(self):
        if self.solver == "umfpack":
            from devsim.umfpack import umfshim
            ds.set_parameter(name="direct_solver", value="custom")
            ds.set_parameter(name="solver_callback", value=umfshim.local_solver_callback)
        elif self.solver == "superlu":
            ds.set_parameter(name="direct_solver", value="superlu")
        else:
            ds.set_parameter(name="direct_solver", value="mkl_pardiso")

    def _mesh(self):
        g, s = self.geom, self.mesh_scale
        lov, L = g.Lov_um * UM, g.L_um * UM
        ltot = 2 * lov + L
        y1 = g.t_gox_nm * NM
        y2 = y1 + g.t_iwo_nm * NM
        y3 = y2 + g.t_pas_nm * NM
        y0 = -2.0 * NM                      # bottom of the (equation-free) gate-metal region
        self.y_if, self.y_top, self.ltot = y1, y2, ltot
        m = self.dev + "_mesh"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ds.create_2d_mesh(mesh=m)
            for pos, ps in ((0.0, 60 * NM), (lov, 4 * NM), (lov + L / 2, 40 * NM), (lov + L, 4 * NM), (ltot, 60 * NM)):
                ds.add_2d_mesh_line(mesh=m, dir="x", pos=pos, ps=ps * s)
            for pos, ps in ((y0, 2.0 * NM), (0.0, 2.0 * NM), (y1, 0.2 * NM), (y2, 0.3 * NM), (y3, 2.0 * NM)):
                ds.add_2d_mesh_line(mesh=m, dir="y", pos=pos, ps=ps * s)
            ds.add_2d_region(mesh=m, material="HfO2", region="gox", xl=0, xh=ltot, yl=0, yh=y1)
            ds.add_2d_region(mesh=m, material="IWO", region="ch", xl=0, xh=ltot, yl=y1, yh=y2)
            ds.add_2d_region(mesh=m, material="HfO2", region="pas", xl=lov, xh=lov + L, yl=y2, yh=y3)
            # metal regions: no equations; they only let the 2D mesher locate the contacts
            ds.add_2d_region(mesh=m, material="metal", region="gmetal", xl=0, xh=ltot, yl=y0, yh=0)
            ds.add_2d_region(mesh=m, material="metal", region="smetal", xl=0, xh=lov, yl=y2, yh=y3)
            ds.add_2d_region(mesh=m, material="metal", region="dmetal", xl=lov + L, xh=ltot, yl=y2, yh=y3)
            b = 1e-10
            ds.add_2d_interface(mesh=m, name="ch_gox", region0="ch", region1="gox", xl=0, xh=ltot, yl=y1, yh=y1, bloat=b)
            ds.add_2d_interface(mesh=m, name="ch_pas", region0="ch", region1="pas", xl=lov, xh=lov + L, yl=y2, yh=y2, bloat=b)
            ds.add_2d_contact(mesh=m, name="gate", material="metal", region="gox", xl=0, xh=ltot, yl=y0, yh=0, bloat=b)
            ds.add_2d_contact(mesh=m, name="source", material="metal", region="ch", xl=0, xh=lov, yl=y2, yh=y3, bloat=b)
            ds.add_2d_contact(mesh=m, name="drain", material="metal", region="ch", xl=lov + L, xh=ltot, yl=y2, yh=y3, bloat=b)
            ds.finalize_mesh(mesh=m)
            ds.create_device(mesh=m, device=self.dev)
        self.mesh_log = buf.getvalue()
        self.contact_nodes = {}
        for line in self.mesh_log.splitlines():
            if line.startswith("Contact ") and " with " in line:
                parts = line.split()
                self.contact_nodes[parts[1]] = int(parts[-2])
        self.n_nodes = sum(len(ds.get_node_model_values(device=self.dev, region=r, name="x"))
                           for r in ("gox", "ch", "pas"))

    def _set_parameters(self):
        d, m = self.dev, self.mat
        ds.set_parameter(device=d, name="V_t", value=m.Vt)
        ds.set_parameter(device=d, name="ElectronCharge", value=Q)
        ds.set_parameter(device=d, name="Nc", value=m.Nc_cm3)
        ds.set_parameter(device=d, name="gate_bias", value=0.0)
        ds.set_parameter(device=d, name="source_bias", value=0.0)
        for r, eps in (("gox", m.eps_hfo2), ("pas", m.eps_hfo2), ("ch", m.eps_ch)):
            ds.set_parameter(device=d, region=r, name="Permittivity", value=eps * EPS0)
        self.push_params(self.p)

    def push_params(self, p: DeviceParams):
        if p.dos_on and self.K == 0:
            raise ValueError("device was built without DOS levels; rebuild with dos_on=True")
        d, m = self.dev, self.mat
        ds.set_parameter(device=d, name="phi_gate", value=p.wf_gate - m.chi_eV)
        ds.set_parameter(device=d, name="phi_sd", value=p.wf_sd - m.chi_eV)
        ds.set_parameter(device=d, name="Qf", value=p.qf)
        ds.set_parameter(device=d, name="ND_eff", value=p.nd_eff)
        ds.set_parameter(device=d, name="mu_n", value=p.mu0)
        if self.K:
            _, wA, wD, n1 = trap_levels(p, m, self.edges)
            for k in range(self.K):
                ds.set_parameter(device=d, name=f"NtA{k}", value=float(wA[k]))
                ds.set_parameter(device=d, name=f"NtD{k}", value=float(wD[k]))
                ds.set_parameter(device=d, name=f"n1_{k}", value=float(n1[k]))
        self.p = p

    # ---------------------------------------------------------------- physics
    def _flux_models(self, r):
        d = self.dev
        ds.node_solution(device=d, region=r, name="Potential")
        ds.edge_from_node_model(device=d, region=r, node_model="Potential")
        ef = "(Potential@n0-Potential@n1)*EdgeInverseLength"
        CreateEdgeModel(d, r, "ElectricField", ef)
        CreateEdgeModelDerivatives(d, r, "ElectricField", ef, "Potential")
        CreateEdgeModel(d, r, "PotentialEdgeFlux", "Permittivity*ElectricField")
        CreateEdgeModelDerivatives(d, r, "PotentialEdgeFlux", "Permittivity*ElectricField", "Potential")

    def _trap_models(self, nname: str, tag: str, deriv_vars):
        d = self.dev
        if self.K == 0:
            ea, ed = "0", "0"
        else:
            ea = " + ".join(f"NtA{k}*{nname}/({nname} + n1_{k})" for k in range(self.K))
            ed = " + ".join(f"NtD{k}*n1_{k}/({nname} + n1_{k})" for k in range(self.K))
        CreateNodeModel(d, "ch", f"TrapA_{tag}", ea)
        CreateNodeModel(d, "ch", f"TrapD_{tag}", ed)
        CreateNodeModelDerivative(d, "ch", f"TrapA_{tag}", ea, *deriv_vars)
        CreateNodeModelDerivative(d, "ch", f"TrapD_{tag}", ed, *deriv_vars)

    def _make_qf_density(self):
        """Interface area per node volume at the gate-dielectric/channel interface."""
        d = self.dev
        x = np.array(ds.get_node_model_values(device=d, region="ch", name="x"))
        y = np.array(ds.get_node_model_values(device=d, region="ch", name="y"))
        vol = np.array(ds.get_node_model_values(device=d, region="ch", name="NodeVolume"))
        idx = np.where(np.abs(y - self.y_if) < 1e-3 * NM)[0]
        order = idx[np.argsort(x[idx])]
        xs = x[order]
        area = np.empty_like(xs)
        area[1:-1] = 0.5 * (xs[2:] - xs[:-2])
        area[0] = 0.5 * (xs[1] - xs[0])
        area[-1] = 0.5 * (xs[-1] - xs[-2])
        apv = np.zeros_like(x)
        apv[order] = area / vol[order]
        ds.node_solution(device=d, region="ch", name="QfAreaPerVolume")
        ds.set_node_values(device=d, region="ch", name="QfAreaPerVolume", values=apv.tolist())
        self.qf_check = float(np.sum(apv * vol) / self.ltot)   # should be 1.0

    def _potential_only(self):
        d = self.dev
        for r in ("gox", "pas", "ch"):
            self._flux_models(r)
        psi0 = self.mat.Vt * math.log(max(self.p.nd_eff, 1e10) / self.mat.Nc_cm3)
        for r in ("gox", "pas", "ch"):
            n = len(ds.get_node_model_values(device=d, region=r, name="x"))
            ds.set_node_values(device=d, region=r, name="Potential", values=[psi0] * n)
        for r in ("gox", "pas"):
            ds.equation(device=d, region=r, name="PotentialEquation", variable_name="Potential",
                        edge_model="PotentialEdgeFlux", variable_update="log_damp")
        self._make_qf_density()
        CreateNodeModel(d, "ch", "QfCharge", "Qf*QfAreaPerVolume")
        neq = "Nc*exp(Potential/V_t)"
        CreateNodeModel(d, "ch", "n_eq", neq)
        CreateNodeModelDerivative(d, "ch", "n_eq", neq, "Potential")
        self._trap_models("n_eq", "eq", ("Potential",))
        pic = "-ElectronCharge*(ND_eff - n_eq + TrapD_eq - TrapA_eq + QfCharge)"
        CreateNodeModel(d, "ch", "PotentialIntrinsicCharge", pic)
        CreateNodeModelDerivative(d, "ch", "PotentialIntrinsicCharge", pic, "Potential")
        ds.equation(device=d, region="ch", name="PotentialEquation", variable_name="Potential",
                    node_model="PotentialIntrinsicCharge", edge_model="PotentialEdgeFlux",
                    variable_update="log_damp")
        for itf in ("ch_gox", "ch_pas"):
            mname = CreateContinuousInterfaceModel(d, itf, "Potential")
            ds.interface_equation(device=d, interface=itf, name="PotentialEquation",
                                  interface_model=mname, type="continuous")
        # gate on HfO2
        CreateContactNodeModel(d, "gate", "gate_bc", "Potential - (gate_bias - phi_gate)")
        CreateContactNodeModel(d, "gate", "gate_bc:Potential", "1")
        ds.contact_equation(device=d, contact="gate", name="PotentialEquation",
                            node_model="gate_bc", edge_charge_model="PotentialEdgeFlux")
        # source (parameter) and drain (circuit node)
        CreateContactNodeModel(d, "source", "source_bc", "Potential - (source_bias - phi_sd)")
        CreateContactNodeModel(d, "source", "source_bc:Potential", "1")
        ds.contact_equation(device=d, contact="source", name="PotentialEquation",
                            node_model="source_bc", edge_charge_model="PotentialEdgeFlux")
        ds.circuit_element(name="V1", n1="drain_bias", n2=0, value=0.0)
        CreateContactNodeModel(d, "drain", "drain_bc", "Potential - (drain_bias - phi_sd)")
        CreateContactNodeModel(d, "drain", "drain_bc:Potential", "1")
        CreateContactNodeModel(d, "drain", "drain_bc:drain_bias", "-1")
        ds.contact_equation(device=d, contact="drain", name="PotentialEquation",
                            node_model="drain_bc", edge_charge_model="PotentialEdgeFlux",
                            circuit_node="drain_bias")

    def _drift_diffusion(self):
        d = self.dev
        ds.node_solution(device=d, region="ch", name="Electrons")
        ds.edge_from_node_model(device=d, region="ch", node_model="Electrons")
        ds.set_node_values(device=d, region="ch", name="Electrons", init_from="n_eq")
        self._trap_models("Electrons", "dd", ("Electrons",))
        pnc = "-ElectronCharge*(ND_eff - Electrons + TrapD_dd - TrapA_dd + QfCharge)"
        CreateNodeModel(d, "ch", "PotentialNodeCharge", pnc)
        CreateNodeModelDerivative(d, "ch", "PotentialNodeCharge", pnc, "Potential", "Electrons")
        ds.equation(device=d, region="ch", name="PotentialEquation", variable_name="Potential",
                    node_model="PotentialNodeCharge", edge_model="PotentialEdgeFlux",
                    variable_update="log_damp")
        CreateBernoulli(d, "ch")
        jn = ("ElectronCharge*mu_n*EdgeInverseLength*V_t*"
              "kahan3(Electrons@n1*Bern01, Electrons@n1*vdiff, -Electrons@n0*Bern01)")
        CreateEdgeModel(d, "ch", "ElectronCurrent", jn)
        for v in ("Electrons", "Potential"):
            CreateEdgeModelDerivatives(d, "ch", "ElectronCurrent", jn, v)
        ds.equation(device=d, region="ch", name="ElectronContinuityEquation", variable_name="Electrons",
                    edge_model="ElectronCurrent", variable_update="positive")
        for c in ("source", "drain"):
            CreateContactNodeModel(d, c, f"{c}_nbc", "Electrons - Nc*exp(-phi_sd/V_t)")
            CreateContactNodeModel(d, c, f"{c}_nbc:Electrons", "1")
            kw = dict(circuit_node="drain_bias") if c == "drain" else {}
            ds.contact_equation(device=d, contact=c, name="ElectronContinuityEquation",
                                node_model=f"{c}_nbc", edge_current_model="ElectronCurrent", **kw)

    # ----------------------------------------------------------------- solve
    def _solve_raw(self, tol):
        t0 = time.perf_counter()
        if self.quiet:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                ds.solve(type="dc", solver_type="direct", **tol)
        else:
            ds.solve(type="dc", solver_type="direct", **tol)
        self.n_solves += 1
        self.t_solve += time.perf_counter() - t0

    def snapshot(self):
        d = self.dev
        snap = {r: list(ds.get_node_model_values(device=d, region=r, name="Potential"))
                for r in ("gox", "pas", "ch")}
        if "Electrons" in ds.get_node_model_list(device=d, region="ch"):
            snap["Electrons"] = list(ds.get_node_model_values(device=d, region="ch", name="Electrons"))
        snap["circuit"] = {nd: ds.get_circuit_node_value(solution="dcop", node=nd)
                           for nd in ("drain_bias", "V1.I")}
        snap["bias"] = (self.vg, self.vd)
        return snap

    def restore(self, snap):
        d = self.dev
        for r in ("gox", "pas", "ch"):
            ds.set_node_values(device=d, region=r, name="Potential", values=snap[r])
        if "Electrons" in snap:
            ds.set_node_values(device=d, region="ch", name="Electrons", values=snap["Electrons"])
        for nd, v in snap["circuit"].items():
            ds.set_circuit_node_value(solution="dcop", node=nd, value=v)
        self._apply_bias(*snap["bias"])

    def _apply_bias(self, vg, vd):
        ds.set_parameter(device=self.dev, name="gate_bias", value=float(vg))
        ds.circuit_alter(name="V1", value=float(vd))
        self.vg, self.vd = float(vg), float(vd)

    def solve(self):
        """Newton solve at rel_tol; on failure one retry at a relaxed tolerance.

        In deep depletion (ID far below the 10 pA floor) the relative update can
        stall around 1e-9 because of nodes with vanishing electron density. The
        relaxed retry (rel 1e-8) is counted in self.n_relaxed; such points are
        not used for gradients (they are below the current floor).
        """
        pre = self.snapshot()
        try:
            self._solve_raw(self.tol)
        except ds.error as exc:
            try:
                self._restore_solution(pre)
                self._solve_raw(dict(self.tol, relative_error=max(1e-8, self.tol["relative_error"]),
                                     maximum_iterations=150))
                self.n_relaxed = getattr(self, "n_relaxed", 0) + 1
            except ds.error:
                self._restore_solution(pre)
                raise ConvergenceError(str(exc)) from exc

    def _restore_solution(self, snap):
        d = self.dev
        for r in ("gox", "pas", "ch"):
            ds.set_node_values(device=d, region=r, name="Potential", values=snap[r])
        if "Electrons" in snap:
            ds.set_node_values(device=d, region="ch", name="Electrons", values=snap["Electrons"])
        for nd, v in snap["circuit"].items():
            ds.set_circuit_node_value(solution="dcop", node=nd, value=v)

    def goto(self, vg=None, vd=None, max_step=0.1, min_frac=1.0 / 512):
        """Move bias to (vg, vd) with adaptive sub-stepping; restores state on failure."""
        tvg = self.vg if vg is None else float(vg)
        tvd = self.vd if vd is None else float(vd)
        frac = 1.0
        while abs(tvg - self.vg) > 1e-12 or abs(tvd - self.vd) > 1e-12:
            dvg, dvd = tvg - self.vg, tvd - self.vd
            nstep = max(1, math.ceil(max(abs(dvg), abs(dvd)) / max_step))
            step = min(1.0, frac / nstep)
            nvg, nvd = self.vg + step * dvg, self.vd + step * dvd
            snap = self.snapshot()
            try:
                self._apply_bias(nvg, nvd)
                self.solve()
                frac = min(1.0, frac * 2.0)
            except ConvergenceError:
                self.restore(snap)
                frac *= 0.5
                if frac < min_frac:
                    raise ConvergenceError(f"bias step failed near VG={nvg:.4f} VD={nvd:.4f}")

    def _calibrate_sign(self):
        snap = self.snapshot()
        self.goto(vd=0.05)
        i = ds.get_circuit_node_value(solution="dcop", node="V1.I")
        self.isign = -1.0 if i < 0 else 1.0
        self.goto(vd=0.0)
        self.restore(snap)

    # --------------------------------------------------------------- outputs
    def drain_current(self) -> float:
        """Total drain current in A (positive for VD > 0)."""
        return self.isign * ds.get_circuit_node_value(solution="dcop", node="V1.I") * self.width_cm

    def contact_current_check(self) -> float:
        """Drain contact electron current (A) from DEVSIM's contact integration, for cross-checking."""
        i2d = ds.get_contact_current(device=self.dev, contact="drain", equation="ElectronContinuityEquation")
        return i2d * self.width_cm

    # -------------------------------------------------------------- parameters
    def get_theta(self, name: str) -> float:
        if name == "gate_bias":
            return self.vg
        if name == "drain_bias":
            return self.vd
        return theta_value(self.p, name)

    def set_theta(self, name: str, value: float):
        """Change one parameter WITHOUT solving (used for residual differencing)."""
        if name == "gate_bias":
            ds.set_parameter(device=self.dev, name="gate_bias", value=float(value))
            self.vg = float(value)
        elif name == "drain_bias":
            ds.circuit_alter(name="V1", value=float(value))
            self.vd = float(value)
        elif name in THETA_DEFS:
            self.push_params(with_theta(self.p, name, value))
        else:
            raise KeyError(name)

    def transfer(self, vd, vgs, on_point=None, max_step=0.1):
        """Sweep VG at fixed VD. on_point(tft, rec) is called after each converged point."""
        self.goto(vd=vd, max_step=max_step)
        self.goto(vg=vgs[0], max_step=max_step)
        recs = []
        for vg in vgs:
            self.goto(vg=vg, max_step=max_step)
            rec = {"vg": float(vg), "vd": float(vd), "id": self.drain_current()}
            if on_point is not None:
                on_point(self, rec)
            recs.append(rec)
        return recs

    def info(self) -> dict:
        return {"nodes": self.n_nodes, "dos_levels": self.K, "mesh_scale": self.mesh_scale,
                "solver": self.solver, "qf_area_check": getattr(self, "qf_check", None),
                "n_solves": self.n_solves, "t_solve_s": round(self.t_solve, 3)}
