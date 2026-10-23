"""Open-source (DEVSIM) TFT simulation with solver-consistent parameter sensitivities."""
from .params import Geometry, Material, DeviceParams, PAPER_SETS, PAPER_J, THETA_ALL, THETA_NO_DOS
from .device import TFT, ConvergenceError
