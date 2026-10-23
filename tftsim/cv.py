"""Quasi-static gate charge and capacitance with adjoint parameter sensitivities.

The gate electrode sits on the HfO2 gate dielectric, which carries no charge,
so the electrode charge is the dielectric flux leaving the gate-contact nodes:

    Q_G = W * sum_{edges (g, j), g on gate} eps_ox * couple_e * (psi_g - psi_j) / l_e ,

a LINEAR functional c^T u of the potential unknowns. Its parameter gradient is
therefore one adjoint solve, J^T lam = c, dQ_G/dtheta = lam^T dR/dtheta (same
residual derivatives as for the drain current). With theta = gate_bias this
gives the quasi-static capacitance C = dQ_G/dV_G exactly (solver-consistent).

The functional is checked against DEVSIM's own contact charge (get_contact_charge,
which uses the same edge flux) in scripts/e7_oed_cv_tch.py.

Quasi-static = all traps follow the gate voltage (low-frequency limit). Real
C-V of oxide TFTs is frequency dependent (trap response, channel RC); a
measurement must be compared with the matching limit.
"""
from __future__ import annotations

import devsim as ds
import numpy as np

NM = 1e-7


class GateCharge:
    def __init__(self, tft, region: str = "gox"):
        self.t = tft
        d = tft.dev
        ds.edge_from_node_model(device=d, region=region, node_model="node_index")
        n0 = np.array(ds.get_edge_model_values(device=d, region=region, name="node_index@n0")).astype(int)
        n1 = np.array(ds.get_edge_model_values(device=d, region=region, name="node_index@n1")).astype(int)
        couple = np.array(ds.get_edge_model_values(device=d, region=region, name="EdgeCouple"))
        inv_l = np.array(ds.get_edge_model_values(device=d, region=region, name="EdgeInverseLength"))
        eps = ds.get_parameter(device=d, region=region, name="Permittivity")
        y = np.array(ds.get_node_model_values(device=d, region=region, name="y"))
        on_gate = np.abs(y - y.min()) < 1e-3 * NM                 # gate contact = bottom face of the dielectric
        g0, g1 = on_gate[n0], on_gate[n1]
        k0 = g0 & ~g1                                            # gate node is n0
        k1 = g1 & ~g0                                            # gate node is n1
        w = tft.width_cm
        coef = eps * couple * inv_l * w
        # Q = sum coef*(psi_gate - psi_other)
        self.gate_nodes = np.concatenate([n0[k0], n1[k1]])
        self.other_nodes = np.concatenate([n1[k0], n0[k1]])
        self.coef = np.concatenate([coef[k0], coef[k1]])
        self.region = region
        self.n_gate_nodes = int(on_gate.sum())

    def value(self) -> float:
        psi = np.array(ds.get_node_model_values(device=self.t.dev, region=self.region, name="Potential"))
        return float(np.sum(self.coef * (psi[self.gate_nodes] - psi[self.other_nodes])))

    def devsim_value(self) -> float:
        """DEVSIM's contact charge (per unit depth) times the width, for checking value()."""
        return ds.get_contact_charge(device=self.t.dev, contact="gate", equation="PotentialEquation") * self.t.width_cm

    def c_vector(self, n: int) -> np.ndarray:
        eq = np.array(ds.get_equation_numbers(device=self.t.dev, region=self.region, equation="PotentialEquation"))
        c = np.zeros(n)
        np.add.at(c, eq[self.gate_nodes], self.coef)
        np.add.at(c, eq[self.other_nodes], -self.coef)
        return c

    def c_ox_area(self) -> float:
        """Geometric oxide capacitance of the whole gate (F): eps_ox * A / t_ox."""
        eps = ds.get_parameter(device=self.t.dev, region=self.region, name="Permittivity")
        return eps * self.t.ltot * self.t.width_cm / (self.t.geom.t_gox_nm * NM)
