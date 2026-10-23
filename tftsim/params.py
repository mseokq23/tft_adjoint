"""Geometry, material constants and calibration parameters for the BG-IWO TFT.

Anchor paper: G. Jeong et al., "Agentic TCAD Calibration Workflow for Oxide
Semiconductor Transistors", arXiv:2609.12184 (v2, 2026).

* Geometry (Table III of the paper): IWO 3.5 nm, HfO2 gate dielectric 7.0 nm,
  conformal HfO2 passivation 5.0 nm, L = 1.0 um (calibration) / 2.0 um (length
  test), W = 5.0 um, S/D overlap 1.5 um per side.
* Calibration descriptors (Sec. III-A and Sec. IV of the paper): S/D and gate
  effective work functions, gate-interface fixed charge Qf, effective donor-like
  background ND,eff, Gaussian (donor-like) and exponential tail (acceptor-like)
  effective DOS, and the mobility prefactor mu0.
* Material constants are NOT given in the paper. The defaults below are the
  widely used a-IGZO TCAD values of Fung et al. (J. Appl. Phys. 106, 084511,
  2009): chi = 4.16 eV, Eg = 3.05 eV, Nc = 5e18 cm^-3, eps = 10. HfO2 eps = 20 is
  an assumption. They are kept fixed (nuisance parameters) and are listed in
  every result file so that their influence can be studied separately.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, replace
import math

import numpy as np
from scipy.special import erf

Q = 1.602176634e-19       # C
KB = 1.380649e-23         # J/K
EPS0 = 8.8541878128e-14   # F/cm


@dataclass(frozen=True)
class Geometry:
    L_um: float = 1.0
    Lov_um: float = 1.5
    W_um: float = 5.0
    t_iwo_nm: float = 3.5
    t_gox_nm: float = 7.0
    t_pas_nm: float = 5.0


@dataclass(frozen=True)
class Material:
    chi_eV: float = 4.16      # electron affinity of the channel (assumed)
    Eg_eV: float = 3.05       # band gap (assumed); only bounds the DOS grid
    Nc_cm3: float = 5.0e18    # effective conduction-band DOS (assumed)
    eps_ch: float = 10.0      # channel relative permittivity (assumed)
    eps_hfo2: float = 20.0    # HfO2 relative permittivity (assumed)
    T_K: float = 300.0

    @property
    def Vt(self) -> float:
        return KB * self.T_K / Q


@dataclass(frozen=True)
class DeviceParams:
    """Calibration descriptors (same model class as arXiv:2609.12184)."""
    wf_sd: float = 4.70       # eV, S/D effective work function
    wf_gate: float = 4.40     # eV, gate work function
    qf: float = 0.0           # cm^-2, fixed charge at gate-dielectric/channel interface
    nd_eff: float = 1.0e18    # cm^-3, effective donor-like background density
    mu0: float = 10.0         # cm^2/Vs, constant mobility prefactor
    dos_on: bool = False      # Gaussian + tail DOS active?
    n_gd: float = 0.0         # cm^-3 eV^-1, donor-like Gaussian peak density
    e_gd: float = 0.12        # eV, Gaussian centre below Ec
    s_gd: float = 0.08        # eV, Gaussian width (sigma)
    n_ta: float = 0.0         # cm^-3 eV^-1, acceptor-like tail density at Ec
    w_ta: float = 0.133       # eV, tail characteristic energy

    def as_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# Parameter sets reported in arXiv:2609.12184
# ---------------------------------------------------------------------------
_STEP0 = DeviceParams()
_STEP1 = replace(_STEP0, wf_sd=4.16)
_STEP2 = replace(_STEP1, wf_gate=4.00, qf=1.59e10)
_STEP3 = replace(_STEP2, nd_eff=2.56e18)
_STEP4 = replace(_STEP3, dos_on=True, n_gd=2.19e19, e_gd=0.12, s_gd=0.08,
                 n_ta=7.8e19, w_ta=0.133)
_STEP5 = replace(_STEP4, mu0=24.0)              # final 2%-W calibrated model (J = 0.181)
# 4%-W recalibration (Sec. IV): ND,eff, gate WF, Qf, tail width and Gaussian DOS changed
_W4 = replace(_STEP5, nd_eff=6.39e16, wf_gate=4.45, qf=-1.0e12, w_ta=0.238,
              n_gd=4.40e19, e_gd=0.080, s_gd=0.224)

PAPER_SETS = {
    "step0": _STEP0,
    "step1": _STEP1,
    "step2": _STEP2,
    "step3": _STEP3,
    "step4": _STEP4,
    "step5": _STEP5,
    "w4": _W4,
}

# Objective J values reported in the paper for each step (2%-W, L = 1 um)
PAPER_J = {"step0": 2.589, "step1": 1.537, "step2": 0.574, "step3": 0.443,
           "step4": 0.282, "step5": 0.181}

# ---------------------------------------------------------------------------
# Calibration vector theta (what the optimiser / sensitivity engine sees)
# ---------------------------------------------------------------------------
# name -> (kind, attribute, typical scale used for scaled sensitivities)
THETA_DEFS = {
    "wf_sd":        ("lin",   "wf_sd",   0.05),    # eV
    "wf_gate":      ("lin",   "wf_gate", 0.05),    # eV
    "qf":           ("lin",   "qf",      1.0e11),  # cm^-2
    "log10_nd_eff": ("log10", "nd_eff",  0.05),    # decades
    "log10_mu0":    ("log10", "mu0",     0.05),
    "log10_n_gd":   ("log10", "n_gd",    0.05),
    "e_gd":         ("lin",   "e_gd",    0.01),    # eV
    "s_gd":         ("lin",   "s_gd",    0.01),    # eV
    "log10_n_ta":   ("log10", "n_ta",    0.05),
    "w_ta":         ("lin",   "w_ta",    0.01),    # eV
}
THETA_DOS = ["log10_n_gd", "e_gd", "s_gd", "log10_n_ta", "w_ta"]
THETA_ALL = list(THETA_DEFS)
THETA_NO_DOS = [t for t in THETA_ALL if t not in THETA_DOS]

# Residual finite-difference steps (no re-solve) and re-solve FD steps (reference)
RES_STEP = {"lin_eV": 1e-5, "qf": 1e6, "log10": 1e-6, "bias": 1e-5}
RESOLVE_STEP = {"lin_eV": 1e-3, "qf": 1e9, "log10": 1e-3, "bias": 1e-3}


def theta_value(p: DeviceParams, name: str) -> float:
    kind, attr, _ = THETA_DEFS[name]
    v = getattr(p, attr)
    return math.log10(v) if kind == "log10" else v


def with_theta(p: DeviceParams, name: str, value: float) -> DeviceParams:
    kind, attr, _ = THETA_DEFS[name]
    return replace(p, **{attr: (10.0 ** value if kind == "log10" else value)})


def theta_vector(p: DeviceParams, names) -> np.ndarray:
    return np.array([theta_value(p, n) for n in names])


def theta_scale(names) -> np.ndarray:
    return np.array([THETA_DEFS[n][2] for n in names])


def step_for(name: str, table: dict) -> float:
    if name in ("gate_bias", "drain_bias"):
        return table["bias"]
    kind, attr, _ = THETA_DEFS[name]
    if kind == "log10":
        return table["log10"]
    if attr == "qf":
        return table["qf"]
    return table["lin_eV"]


# ---------------------------------------------------------------------------
# Discretised DOS (fixed energy grid; only the weights depend on theta)
# ---------------------------------------------------------------------------
def default_dos_edges(fine_step=0.025, fine_until=0.8, coarse_step=0.1, emax=1.6):
    """Energy-bin edges in eV below Ec (Delta E = Ec - E >= 0)."""
    fine = np.arange(0.0, fine_until + 1e-12, fine_step)
    coarse = np.arange(fine_until + coarse_step, emax + 1e-12, coarse_step)
    return np.concatenate([fine, coarse])


def trap_levels(p: DeviceParams, mat: Material, edges: np.ndarray):
    """Return (dE_centres, N_acceptor_k, N_donor_k, n1_k) for the discretised DOS.

    * acceptor-like tail   g_A(dE) = n_ta * exp(-dE / w_ta)
    * donor-like Gaussian  g_D(dE) = n_gd * exp(-(dE - e_gd)^2 / (2 s_gd^2)),
      truncated to the gap (dE >= 0)
    Bin weights are exact integrals over each bin, so they are smooth functions
    of the DOS parameters. Occupancy uses n1_k = Nc exp(-dE_k / Vt), i.e.
    f_k = n / (n + n1_k) (traps in quasi-equilibrium with the conduction band).
    """
    a, b = edges[:-1], edges[1:]
    centres = 0.5 * (a + b)
    if not p.dos_on:
        z = np.zeros_like(centres)
        return centres, z, z, mat.Nc_cm3 * np.exp(-centres / mat.Vt)
    wA = p.n_ta * p.w_ta * (np.exp(-a / p.w_ta) - np.exp(-b / p.w_ta))
    s2 = p.s_gd * math.sqrt(2.0)
    wD = p.n_gd * p.s_gd * math.sqrt(math.pi / 2.0) * (erf((b - p.e_gd) / s2) - erf((a - p.e_gd) / s2))
    n1 = mat.Nc_cm3 * np.exp(-centres / mat.Vt)
    return centres, wA, wD, n1
