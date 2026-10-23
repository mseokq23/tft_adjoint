"""Transfer-curve metrics (as defined in arXiv:2609.12184) and their gradients.

* Vth : VG at |ID| = 100 nA (log-linear interpolation between sweep points)
* SS  : average inverse slope between 0.1 nA and 100 nA  = (VG@100nA - VG@0.1nA)/3
* Ion : maximum |ID| of the sweep
* currents below the 10 pA floor are excluded from log-current errors

Gradients of interpolated crossings are exact derivatives of the reported
interpolated metric: for a crossing between points a and b with
t = (Lt - La)/(Lb - La), L = log10 I,
    dV/dtheta = -(Vb - Va)/(Lb - La) * ((1 - t) dLa/dtheta + t dLb/dtheta).
"""
from __future__ import annotations

import math

import numpy as np

I_TH = 100e-9
I_SS_LO = 0.1e-9
I_FLOOR = 10e-12
LN10 = math.log(10.0)


def crossing(vg, i, itarget, di=None):
    vg = np.asarray(vg, float)
    L = np.log10(np.maximum(np.asarray(i, float), 1e-300))
    lt = math.log10(itarget)
    idx = np.where((L[:-1] < lt) & (L[1:] >= lt))[0]
    if len(idx) == 0:
        return math.nan, None
    a = int(idx[0])
    b = a + 1
    t = (lt - L[a]) / (L[b] - L[a])
    v = vg[a] + t * (vg[b] - vg[a])
    if di is None:
        return float(v), None
    di = np.asarray(di, float)
    dla = di[a] / (i[a] * LN10)
    dlb = di[b] / (i[b] * LN10)
    dv = -(vg[b] - vg[a]) / (L[b] - L[a]) * ((1 - t) * dla + t * dlb)
    return float(v), dv


def extract(vg, i, di=None) -> dict:
    """Vth (V), SS (mV/dec), Ion (A) and, if di (npts x P) is given, their gradients."""
    vg = np.asarray(vg, float)
    i = np.asarray(i, float)
    vth, dvth = crossing(vg, i, I_TH, di)
    vlo, dvlo = crossing(vg, i, I_SS_LO, di)
    ss = (vth - vlo) / 3.0 * 1e3
    k = int(np.argmax(i))
    out = {"vth": vth, "ss_mVdec": ss, "ion": float(i[k]), "log10_ion": math.log10(i[k]),
           "vg_ion": float(vg[k])}
    if di is not None:
        di = np.asarray(di, float)
        out["d_vth"] = dvth
        out["d_ss_mVdec"] = None if (dvth is None or dvlo is None) else (dvth - dvlo) / 3.0 * 1e3
        out["d_log10_ion"] = di[k] / (i[k] * LN10)
    return out


def log_current_rows(i, di, floor=I_FLOOR, mask_from=None):
    """Rows of d log10(ID)/d theta for points above the current floor."""
    i = np.asarray(i, float)
    ref = i if mask_from is None else np.asarray(mask_from, float)
    m = ref >= floor
    return np.asarray(di, float)[m] / (i[m, None] * LN10), m


# Normalisation constants a_i: the paper does not publish them, so these are our
# own choices, fixed before any calibration and reported with every result.
A_NORM = {"log": 1.0, "lin": None, "th": 0.1, "on": 1.0, "ss": 0.1}   # dec, rel-to-Ion, V, dec, V/dec
W_PAPER = {"log": 0.25, "lin": 0.10, "th": 1.00, "on": 0.20, "ss": 0.25}


def objective_J(sim: dict, meas: dict, weights=W_PAPER, a=A_NORM) -> dict:
    """J of arXiv:2609.12184, Eq. (1), with our own normalisation a_i.

    sim/meas: {vd: (vg array, id array)} on the same VG grid (interpolate first).
    """
    terms = []
    for vd, (vg_m, i_m) in meas.items():
        vg_s, i_s = sim[vd]
        i_s = np.interp(vg_m, vg_s, i_s)
        m_log = i_m >= I_FLOOR
        e_log = float(np.sqrt(np.mean((np.log10(i_s[m_log]) - np.log10(i_m[m_log])) ** 2)))
        ion_m = float(np.max(i_m))
        m_lin = i_m >= ion_m / 10.0
        e_lin = float(np.sqrt(np.mean((i_s[m_lin] - i_m[m_lin]) ** 2))) / ion_m
        ms, mm = extract(vg_m, i_s), extract(vg_m, i_m)
        t = {
            "vd": vd,
            "E_log_dec": e_log,
            "E_lin_rel": e_lin,
            "dVth_V": ms["vth"] - mm["vth"],
            "dlogIon_dec": ms["log10_ion"] - mm["log10_ion"],
            "dSS_Vdec": (ms["ss_mVdec"] - mm["ss_mVdec"]) / 1e3,
        }
        t["J_vd"] = (weights["log"] * e_log / a["log"] + weights["lin"] * e_lin
                     + weights["th"] * abs(t["dVth_V"]) / a["th"] + weights["on"] * abs(t["dlogIon_dec"]) / a["on"]
                     + weights["ss"] * abs(t["dSS_Vdec"]) / a["ss"])
        terms.append(t)
    return {"J": float(np.mean([t["J_vd"] for t in terms])), "terms": terms}
