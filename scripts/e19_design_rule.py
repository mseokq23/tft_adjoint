"""E19 - design rule: which kinds of test can tell calibration-equivalent parameter sets apart?

Linear check at the twin truth (paper Step 5). For every candidate test: the predictive 1σ of the
test given the calibration data (L = 1 µm, VD = 0.05/1 V), in units of the test's own noise
(σ = 0.02 dec on log10 ID), as in E14 section 3 (Qf fixed: it is exactly degenerate with the gate
WF for one gate stack). Below 1, the sets that fit the calibration curves predict the same test
result within its noise (in the linear approximation), so passing the test carries almost no
information about the parameters.

Candidates (one curve each, VG from -1.5 to 2.5 V in 0.05 V steps, points with ID >= 10 pA):
  lateral / bias (vertical stack unchanged): L = 1 µm at VD 0.02/0.1/0.2/0.5/2/3 V,
      L = 0.5/2/5/10 µm at VD 0.05/1 V (L = 0.5 and 10 µm also at VD = 3 V)
  vertical (stack changed): t_ox = 5/10/14/20 nm and t_IWO = 2/5/7/10 nm at VD 0.05/1 V
The calibration curves themselves (L = 1 µm, VD 0.05/1 V) are recomputed in the same run and give
the reference rows; they are checked against results/e3/rows.npz.

For each test the report also gives
  * the largest single-point value inside the calibration window (VG where the calibration curves have
    points above 10 pA; a test with more current reaches further into subthreshold),
  * whole-curve factors: λ_max = the largest χ² the test curve can add per unit of calibration Δχ²
    (a set on the calibration's Δχ² = c boundary adds at most c·λ_max), and q = the same along the
    linearised ND,eff profile (other descriptors re-fitted; the two calibration curves sum to 1),
  * Δlog det F and the 1σ of every descriptor after adding the test, relative to calibration only
    - i.e. WHAT the test identifies.

The 1-D Pao–Sah model (scripts/aux_paosah_1d.py) is the analytic reference: there, every VD and L
curve is fixed by one function µ·N_s(u), so within the range of u where the calibration curves fix
that function, lateral and bias tests add nothing. E19 shows how far the 2D model departs from that
(deeper subthreshold or higher VD than the calibration, contacts and overlap regions, short channel).

Nonlinear check (--nonlinear): the re-fitted sets of the twin ND,eff profile (E5, χ² < 3.84 on the
calibration curves) are simulated at every lateral/bias test and compared with the Step-5 truth:
χ² of the test (σ = 0.02 dec, VG step 0.1 V and ID >= 10 pA as in E5/E9). χ² > 3.84 means the test
rejects that set. The vertical tests were checked the same way in E9 (t_IWO 7 nm, t_ox 14 nm).

Usage:
  python scripts/e19_design_rule.py --group lateral --out results/e19 > results/e19_lateral_log.txt
  python scripts/e19_design_rule.py --group vertical --out results/e19 > results/e19_vertical_log.txt
  python scripts/e19_design_rule.py --nonlinear --values 15.9,16.9,17.4 --tag a --out results/e19 > results/e19_nl_a_log.txt
  python scripts/e19_design_rule.py --nonlinear --values 17.9,18.9 --tag b --out results/e19 > results/e19_nl_b_log.txt
  python scripts/e19_design_rule.py --report --out results/e19            # tables and figure, no simulation
(runs can go in parallel; each writes its own rows_<group>.npz / nonlinear_<tag>.json as it goes)
"""
from __future__ import annotations

import argparse
import dataclasses
import glob
import json
import os

from _common import HAVE_MPL, Timer, dump_json, np, outdir, vg_grid
from tftsim import PAPER_SETS, THETA_ALL, TFT, ConvergenceError, Geometry
from tftsim.metrics import I_FLOOR, LN10, extract
from tftsim.params import theta_scale, with_theta
from tftsim.sensitivity import SensitivityEngine

NAMES = THETA_ALL
SIGMA = 0.02
CAL = ("L1_VD0.05", "L1_VD1")
PAPER = {"L1_VD0.1", "L1_VD0.2", "L1_VD0.5", "L2_VD0.05", "L2_VD1"}

# device label -> (Geometry, drain biases, group)
DEVICES = {
    "L1": (Geometry(), (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0), "lateral"),
    "L0.5": (Geometry(L_um=0.5), (0.05, 1.0, 3.0), "lateral"),
    "L2": (Geometry(L_um=2.0), (0.05, 1.0), "lateral"),
    "L5": (Geometry(L_um=5.0), (0.05, 1.0), "lateral"),
    "L10": (Geometry(L_um=10.0), (0.05, 1.0, 3.0), "lateral"),
    "tox5": (Geometry(t_gox_nm=5.0), (0.05, 1.0), "vertical"),
    "tox10": (Geometry(t_gox_nm=10.0), (0.05, 1.0), "vertical"),
    "tox14": (Geometry(t_gox_nm=14.0), (0.05, 1.0), "vertical"),
    "tox20": (Geometry(t_gox_nm=20.0), (0.05, 1.0), "vertical"),
    "tiwo2": (Geometry(t_iwo_nm=2.0), (0.05, 1.0), "vertical"),
    "tiwo5": (Geometry(t_iwo_nm=5.0), (0.05, 1.0), "vertical"),
    "tiwo7": (Geometry(t_iwo_nm=7.0), (0.05, 1.0), "vertical"),
    "tiwo10": (Geometry(t_iwo_nm=10.0), (0.05, 1.0), "vertical"),
}


def key(dev, vd):
    return f"{dev}_VD{vd:g}"


def sweep_device(params, geom, vds, vgs, log):
    """Transfer curves with adjoint d ID/d theta at every point above 0.1 x the current floor.
    A drain bias that does not converge is recorded and skipped; the others are kept."""
    t = TFT(params, geom).build()
    eng = SensitivityEngine(t, central=False)
    out, failed = {}, {}
    for j, vd in enumerate(vds):
        seq = vgs if j % 2 == 0 else vgs[::-1]

        def cb(tft, rec):
            if rec["id"] >= 0.1 * I_FLOOR:
                rec["g"] = eng.current_gradient(NAMES, "adjoint")

        try:
            recs = t.transfer(vd, seq, on_point=cb)
        except ConvergenceError as exc:
            failed[vd] = str(exc)
            log(f"  VD = {vd:g} V failed: {exc}")
            t = TFT(params, geom).build()                    # fresh device for the next drain bias
            eng = SensitivityEngine(t, central=False)
            continue
        recs.sort(key=lambda r: r["vg"])
        vg = np.array([r["vg"] for r in recs])
        i = np.array([r["id"] for r in recs])
        G = np.array([r.get("g", np.zeros(len(NAMES))) for r in recs])
        out[vd] = (vg, i, G)
    return out, failed, t.info()


def rows_of(vg, i, G, scale):
    m = i >= I_FLOOR
    return G[m] / (i[m, None] * LN10) * scale


def run_group(group, out, vgs, only=None):
    clock = Timer()
    scale = theta_scale(NAMES)
    p = PAPER_SETS["step5"]
    path = os.path.join(out, f"rows_{group}.npz")
    meta_path = os.path.join(out, f"meta_{group}.json")
    rows, meta = {}, {}
    if os.path.exists(path) and only:                     # re-run selected devices, keep the others
        z = np.load(path)
        rows = {k[2:]: z[k] for k in z.files if k.startswith("S_")}
        meta = json.load(open(meta_path))
    def log(msg):
        print(msg, flush=True)

    for dev, (geom, vds, grp) in DEVICES.items():
        if grp != group or (only and dev not in only):
            continue
        t0 = clock()
        try:
            data, failed, info = sweep_device(p, geom, vds, vgs, log)
        except ConvergenceError as exc:
            log(f"[{dev}] build failed: {exc}")
            meta[dev] = {"error": str(exc)}
            continue
        for vd, (vg, i, G) in data.items():
            k = key(dev, vd)
            rows[k] = rows_of(vg, i, G, scale)
            m = extract(vg, i)
            meta[k] = {"device": dev, "vd": vd, "geom": dataclasses.asdict(geom), "n_rows": int(rows[k].shape[0]),
                       "vth_V": m["vth"], "ss_mVdec": m["ss_mVdec"], "ion_A": m["ion"],
                       "i_at_vgmin_A": float(i[0]), "off_reached": bool(i[0] < I_FLOOR),
                       "vg_range_used": [float(vg[i >= I_FLOOR].min()), float(vg[i >= I_FLOOR].max())]}
        for vd, msg in failed.items():
            meta[key(dev, vd)] = {"device": dev, "vd": vd, "error": msg}
        meta[dev] = {"nodes": info["nodes"], "wall_s": round(clock() - t0, 1), "failed_vd": list(failed)}
        np.savez(path, names=np.array(NAMES), scale=scale, vgs=vgs, **{f"S_{k}": v for k, v in rows.items()})
        dump_json(meta, meta_path)
        log(f"[{dev}] {len(data)}/{len(vds)} curves, {info['nodes']} nodes, {clock() - t0:.0f} s (total {clock():.0f} s)")


def params_of(theta):
    p = PAPER_SETS["step5"]
    for n, v in theta.items():
        p = with_theta(p, n, float(v))
    return p


def forward(params, geom, vds, vgs):
    """Transfer curves without gradients; a drain bias that does not converge is skipped."""
    t = TFT(params, geom).build()
    out, failed = {}, {}
    for j, vd in enumerate(vds):
        seq = vgs if j % 2 == 0 else vgs[::-1]
        try:
            recs = sorted(t.transfer(vd, seq), key=lambda r: r["vg"])
        except ConvergenceError as exc:
            failed[vd] = str(exc)
            t = TFT(params, geom).build()
            continue
        out[vd] = np.array([r["id"] for r in recs])
    return out, failed


def run_nonlinear(out, sets_path, values, tag):
    """E5 sets (calibration χ² < 3.84) vs the Step-5 truth at every lateral/bias test (not the calibration curves)."""
    clock = Timer()
    recs = {round(r["value"], 4): r for r in json.load(open(sets_path))["records"] if "theta" in r and r["chi2"] < 3.84}
    truth = min(recs, key=lambda k: recs[k]["chi2"])
    keep = None if values is None else {round(float(v), 4) for v in values.split(",")}
    vgs = vg_grid(-1.5, 2.5, 0.1)
    conds = {dev: (g, tuple(vd for vd in vds if key(dev, vd) not in CAL)) for dev, (g, vds, grp) in DEVICES.items()
             if grp == "lateral"}
    ref = {}
    for dev, (g, vds) in conds.items():
        ref[dev], failed = forward(PAPER_SETS["step5"], g, vds, vgs)
        if failed:
            print(f"[truth {dev}] failed at VD {list(failed)}", flush=True)
    print(f"truth curves done ({clock():.0f} s)", flush=True)
    path = os.path.join(out, f"nonlinear_{tag}.json")
    res = {"sigma_dec": SIGMA, "vg_step": 0.1, "sets": sets_path, "results": {}}
    for k in sorted(recs):
        if k == truth or (keep is not None and k not in keep):
            continue
        p = params_of(recs[k]["theta"])
        row = {"chi2_calibration": recs[k]["chi2"]}
        for dev, (g, vds) in conds.items():
            cur, failed = forward(p, g, vds, vgs)
            for vd in vds:
                if vd not in cur or vd not in ref[dev]:
                    row[key(dev, vd)] = {"error": failed.get(vd, "reference failed")}
                    continue
                a, b = ref[dev][vd], cur[vd]
                m = a >= I_FLOOR
                dev_log = np.log10(np.maximum(b[m], 1e-30)) - np.log10(a[m])
                ea, eb = extract(vgs, a), extract(vgs, b)
                row[key(dev, vd)] = {"chi2": float(np.sum((dev_log / SIGMA) ** 2)), "rms_dec": float(np.sqrt(np.mean(dev_log ** 2))),
                                     "max_dec": float(np.max(np.abs(dev_log))), "n": int(m.sum()),
                                     "dVth_mV": 1e3 * (eb["vth"] - ea["vth"])}
        res["results"][f"{k:.4f}"] = row
        dump_json(res, path)
        worst = max((v["chi2"], kk) for kk, v in row.items() if isinstance(v, dict) and "chi2" in v)
        print(f"[ND = 10^{k}] largest χ² {worst[0]:.3g} ({worst[1]})  ({clock():.0f} s)", flush=True)


def load_rows(out):
    rows, meta = {}, {}
    for f in sorted(glob.glob(os.path.join(out, "rows_*.npz"))):
        z = np.load(f)
        rows.update({k[2:]: z[k] for k in z.files if k.startswith("S_")})
        names = [str(n) for n in z["names"]]
    for f in sorted(glob.glob(os.path.join(out, "meta_*.json"))):
        meta.update(json.load(open(f)))
    return rows, meta, names


def vg_of(meta_k, n):
    """VG of the stored rows (contiguous points above the floor on the 0.05 V grid), or None."""
    if not meta_k or "vg_range_used" not in meta_k:
        return None
    a, b = meta_k["vg_range_used"]
    vg = np.round(np.arange(a, b + 1e-9, 0.05), 6)
    return vg if vg.size == n else None


def analyse(rows, names, meta=None):
    meta = meta or {}
    keep = [k for k, n in enumerate(names) if n != "qf"]
    kn = [names[k] for k in keep]
    S0 = np.vstack([rows[c] for c in CAL])[:, keep]
    F0 = S0.T @ S0 / SIGMA ** 2
    Finv = np.linalg.inv(F0)
    sd0 = np.sqrt(np.diag(Finv))
    ld0 = np.linalg.slogdet(F0)[1]
    w0, V0 = np.linalg.eigh(F0)
    Fmh = V0 @ np.diag(w0 ** -0.5) @ V0.T                   # F^-1/2
    ind = kn.index("log10_nd_eff")
    dprof = Finv[:, ind] / Finv[ind, ind]                     # linearised ND,eff profile: other descriptors re-fitted

    def factors(S):
        """lam_max: largest test χ² per unit calibration Δχ² over all directions (a set on the calibration's Δχ² = c
        boundary adds at most c·lam_max to the test); q_prof: the same along the linearised ND,eff profile."""
        G = Fmh @ (S.T @ S / SIGMA ** 2) @ Fmh
        return float(np.linalg.eigvalsh(G)[-1]), float(np.sum((S @ dprof) ** 2) / SIGMA ** 2 * Finv[ind, ind])

    cal_vg = [vg_of(meta.get(c), rows[c].shape[0]) for c in CAL]
    vg_lo = max(v[0] for v in cal_vg) if all(v is not None for v in cal_vg) else None   # calibration window starts here
    res = {}
    for k, S in rows.items():
        S = S[:, keep]
        v = np.einsum("ij,jk,ik->i", S, Finv, S)
        F1 = F0 + S.T @ S / SIGMA ** 2
        sd1 = np.sqrt(np.diag(np.linalg.inv(F1)))
        r = np.sqrt(v) / SIGMA
        lam, qp = factors(S)
        res[k] = {"n_rows": int(S.shape[0]), "rms_ratio": float(np.sqrt(np.mean(v)) / SIGMA),
                  "max_ratio": float(r.max()), "dlogdet": float(np.linalg.slogdet(F1)[1] - ld0),
                  "lam_max": lam, "q_profile": qp, "sd_ratio": dict(zip(kn, (sd1 / sd0).tolist()))}
        vg = vg_of(meta.get(k), S.shape[0])
        if vg is not None:
            res[k]["vg_of_max"] = float(vg[np.argmax(r)])
            res[k]["vg_first"] = float(vg[0])
            if vg_lo is not None:
                w = vg >= vg_lo - 1e-9
                res[k]["max_ratio_in_cal_window"] = float(r[w].max())
                res[k]["n_below_cal_window"] = int((~w).sum())
    pairs = {}
    for dev in DEVICES:
        a, b = key(dev, 0.05), key(dev, 1.0)
        if a in rows and b in rows and dev != "L1":
            S = np.vstack([rows[a], rows[b]])[:, keep]
            v = np.einsum("ij,jk,ik->i", S, Finv, S)
            F1 = F0 + S.T @ S / SIGMA ** 2
            sd1 = np.sqrt(np.diag(np.linalg.inv(F1)))
            lam, qp = factors(S)
            pairs[dev] = {"rms_ratio": float(np.sqrt(np.mean(v)) / SIGMA), "max_ratio": float(np.sqrt(np.max(v)) / SIGMA),
                          "dlogdet": float(np.linalg.slogdet(F1)[1] - ld0), "lam_max": lam, "q_profile": qp,
                          "sd_ratio": dict(zip(kn, (sd1 / sd0).tolist()))}
    return res, pairs, dict(zip(kn, sd0.tolist())), S0, vg_lo


def report(out):
    rows, meta, names = load_rows(out)
    missing = [c for c in CAL if c not in rows]
    if missing:
        raise SystemExit(f"calibration rows missing: {missing} (run --group lateral first)")
    res, pairs, sd0, S0, vg_lo = analyse(rows, names, meta)
    vg_lo = float("nan") if vg_lo is None else vg_lo
    # consistency with the stored E3/E4/E7/E14 analysis
    checks = {}
    if os.path.exists("results/e3/rows.npz"):
        z3 = np.load("results/e3/rows.npz")
        keep = [k for k, n in enumerate(names) if n != "qf"]
        Sb = z3["S_base"][:, keep]
        checks["S_base_shape_e3"] = list(Sb.shape)
        checks["S_base_shape_e19"] = list(S0.shape)
        if Sb.shape == S0.shape:
            checks["S_base_max_rel_diff"] = float(np.max(np.abs(Sb - S0)) / np.max(np.abs(Sb)))
    e14 = json.load(open("results/e14/summary.json"))["linear"] if os.path.exists("results/e14/summary.json") else {}
    ref = {"L1_VD0.1": "VD 0.1 V (paper)", "L1_VD0.2": "VD 0.2 V (paper)", "L1_VD0.5": "VD 0.5 V (paper)",
           "L2_VD0.05": "L 2 µm, VD 0.05 V (paper)", "L2_VD1": "L 2 µm, VD 1 V (paper)"}
    ref_pairs = {"tox10": "t_ox 10 nm", "tox14": "t_ox 14 nm", "tiwo5": "t_IWO 5 nm", "tiwo7": "t_IWO 7 nm"}
    checks["vs_e14_rms_ratio"] = {k: [res[k]["rms_ratio"], e14[r]["rms_ratio"]] for k, r in ref.items() if k in res and r in e14}
    checks["vs_e14_rms_ratio"].update({k: [pairs[k]["rms_ratio"], e14[r]["rms_ratio"]]
                                       for k, r in ref_pairs.items() if k in pairs and r in e14})
    if os.path.exists("results/e7/rows.npz") and "quasi-static C–V" in e14:
        checks["cv_reference_e14"] = e14["quasi-static C–V"]

    groups = {k: DEVICES[k.split("_VD")[0]][2] for k in res}
    lat = {k: v for k, v in res.items() if groups[k] == "lateral" and k not in CAL}
    ver = {k: v for k, v in res.items() if groups[k] == "vertical"}
    show = ["log10_nd_eff", "log10_n_ta", "log10_n_gd", "wf_gate", "wf_sd", "log10_mu0"]

    def order(k):
        dev, vd = k.split("_VD")
        return (list(DEVICES).index(dev), float(vd))

    lines = ["# E19: which kinds of test can tell calibration-equivalent parameter sets apart?", "",
             "Linear check at the twin truth (Step 5). Calibration data: L = 1 µm, VD = 0.05/1 V. Noise σ = 0.02 dec on log10 ID "
             "for the calibration and for every test; Qf fixed (exactly degenerate with the gate WF for one gate stack). "
             "σ_pred/σ = predictive 1σ of the test given the calibration data, in units of the test's noise (RMS and max over "
             "the test curve). Below 1, the calibration-equivalent sets predict the same test result within its noise. "
             "sd ratio = 1σ of a descriptor after adding the test / 1σ from calibration only (1 = the test adds nothing about it). "
             "λ_max = the largest χ² the whole test curve can add per unit of calibration Δχ² (a set on the calibration's "
             "Δχ² = c boundary adds at most c·λ_max; below 1 the test cannot reject any set the calibration accepts). "
             "q = the same along the linearised ND,eff profile (other descriptors re-fitted); the two calibration curves "
             "sum to 1.", "",
             f"\"max in window\" = the largest single-point value at VG ≥ {vg_lo:+.2f} V, where the calibration curves have "
             "points above the 10 pA floor (a test whose current is larger reaches further into subthreshold).", "",
             "| test | group | rows (first VG) | σ_pred/σ RMS | σ_pred/σ max (at VG) | max in window | λ_max | q (ND,eff profile) | "
             "Δlog det F | sd ratio: " + " / ".join(show) + " |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for k in sorted(res, key=order):
        r = res[k]
        tag = " (calibration)" if k in CAL else " (paper)" if k in PAPER else ""
        lines.append(f"| {k}{tag} | {groups[k]} | {r['n_rows']} ({r.get('vg_first', float('nan')):+.2f}) | {r['rms_ratio']:.3g} | "
                     f"{r['max_ratio']:.3g} ({r.get('vg_of_max', float('nan')):+.2f}) | {r.get('max_ratio_in_cal_window', float('nan')):.3g} | "
                     f"{r['lam_max']:.3g} | {r['q_profile']:.3g} | {r['dlogdet']:.2f} | "
                     + " / ".join(f"{r['sd_ratio'][n]:.2f}" for n in show) + " |")
    lines += ["", "Two-curve tests (VD = 0.05 and 1 V together), as in E14:", "",
              "| device | σ_pred/σ RMS | σ_pred/σ max | λ_max | q (ND,eff profile) | Δlog det F | sd ratio: " + " / ".join(show) + " |",
              "|---|---|---|---|---|---|---|"]
    for dev, r in pairs.items():
        lines.append(f"| {dev} | {r['rms_ratio']:.3g} | {r['max_ratio']:.3g} | {r['lam_max']:.3g} | {r['q_profile']:.3g} | "
                     f"{r['dlogdet']:.2f} | " + " / ".join(f"{r['sd_ratio'][n]:.2f}" for n in show) + " |")
    lat_max = max(lat.values(), key=lambda r: r["rms_ratio"]) if lat else None
    verdict = {}
    if lat and ver:
        lk = max(lat, key=lambda k: lat[k]["rms_ratio"])
        lkm = max(lat, key=lambda k: lat[k]["max_ratio"])
        lkw = max(lat, key=lambda k: lat[k].get("max_ratio_in_cal_window", 0.0))
        vk = min(ver, key=lambda k: ver[k]["rms_ratio"])
        verdict = {"lateral_max_rms": [lk, lat[lk]["rms_ratio"]], "lateral_max_max": [lkm, lat[lkm]["max_ratio"]],
                   "lateral_max_in_window": [lkw, lat[lkw].get("max_ratio_in_cal_window")],
                   "vertical_min_rms": [vk, ver[vk]["rms_ratio"]],
                   "lateral_min_sd_ratio": {n: min(r["sd_ratio"][n] for r in lat.values()) for n in show},
                   "vertical_min_sd_ratio": {n: min(r["sd_ratio"][n] for r in ver.values()) for n in show}}
        lines += ["", "## Verdict", "",
                  f"- Lateral/bias tests (calibration curves excluded): σ_pred/σ RMS at most {lat[lk]['rms_ratio']:.3g} ({lk}); "
                  f"the largest single-point value is {lat[lkm]['max_ratio']:.3g} ({lkm}, VG = {lat[lkm].get('vg_of_max', float('nan')):+.2f} V); "
                  f"inside the calibration window (VG ≥ {vg_lo:+.2f} V) it is {lat[lkw].get('max_ratio_in_cal_window', float('nan')):.3g} ({lkw}).",
                  f"- Vertical tests: σ_pred/σ RMS at least {ver[vk]['rms_ratio']:.3g} ({vk}).",
                  "- Whole-curve factors: lateral/bias λ_max {:.3g}–{:.3g}, q {:.3g}–{:.3g}; vertical λ_max {:.3g}–{:.3g}, q {:.3g}–{:.3g}.".format(
                      min(r["lam_max"] for r in lat.values()), max(r["lam_max"] for r in lat.values()),
                      min(r["q_profile"] for r in lat.values()), max(r["q_profile"] for r in lat.values()),
                      min(r["lam_max"] for r in ver.values()), max(r["lam_max"] for r in ver.values()),
                      min(r["q_profile"] for r in ver.values()), max(r["q_profile"] for r in ver.values())),
                  "- Smallest sd ratio reached by any single lateral/bias test: "
                  + ", ".join(f"{n} {v:.2f}" for n, v in verdict["lateral_min_sd_ratio"].items()) + ".",
                  "- Smallest sd ratio reached by any single vertical test: "
                  + ", ".join(f"{n} {v:.2f}" for n, v in verdict["vertical_min_sd_ratio"].items()) + "."]
    bad = {k: v for k, v in meta.items() if isinstance(v, dict) and ("error" in v or v.get("off_reached") is False)}
    if bad:
        lines += ["", "Notes on the sweeps: " + "; ".join(
            f"{k}: " + (v["error"] if "error" in v else f"ID at VG = -1.5 V is {v['i_at_vgmin_A']:.1e} A (off-state not reached)")
            for k, v in bad.items()) + "."]
    nl = {}
    for f in sorted(glob.glob(os.path.join(out, "nonlinear_*.json"))):
        nl.update(json.load(open(f))["results"])
    if nl:
        ck = sorted({c for r in nl.values() for c, v in r.items() if isinstance(v, dict)}, key=order)
        lines += ["", "## Nonlinear check of the lateral/bias tests (twin E5 sets vs the Step-5 truth)", "",
                  "χ² that each test would add (σ = 0.02 dec, VG step 0.1 V, ID ≥ 10 pA, as in E5/E9); > 3.84 rejects the set. "
                  "ΔVth at 100 nA in brackets (mV).", "",
                  "| log10 ND,eff | calibration χ² | " + " | ".join(ck) + " | sum over the tests | max abs(ΔVth) (mV) |",
                  "|---|---|" + "---|" * (len(ck) + 2)]
        sums = {}
        for v in sorted(nl, key=float):
            r = nl[v]
            cells = [(f"{r[c]['chi2']:.3g} ({r[c]['dVth_mV']:+.2f})" if "chi2" in r.get(c, {}) else "–") for c in ck]
            got = [r[c] for c in ck if "chi2" in r.get(c, {})]
            sums[v] = sum(g["chi2"] for g in got)
            lines.append(f"| {float(v):.2f} | {r['chi2_calibration']:.3g} | " + " | ".join(cells)
                         + f" | {sums[v]:.3g} | {max(abs(g['dVth_mV']) for g in got):.2f} |")
        verdict["nonlinear_sum_chi2"] = sums
        allc = [(r[c]["chi2"], v, c) for v, r in nl.items() for c in ck if "chi2" in r.get(c, {})]
        top = max(allc)
        verdict["nonlinear_max_chi2"] = {"chi2": top[0], "set": top[1], "test": top[2]}
        cmp_txt = ""
        e9 = {}
        for f in sorted(glob.glob("results/e9[a-z]/summary.json")):          # E9 was run in two parts (e9a, e9b)
            e9.update(json.load(open(f))["results"])
        if e9:
            rg = {c: [min(r[c]["chi2"] for r in e9.values()), max(r[c]["chi2"] for r in e9.values())] for c in ("t_IWO 7 nm", "t_ox 14 nm")}
            verdict["e9_chi2_range"] = rg
            cmp_txt = " For comparison, E9 (vertical tests, same sets, VD 0.05/1 V together): " + ", ".join(
                f"{c} χ² {a:.3g}–{b:.3g}" for c, (a, b) in rg.items()) + "."
        lines += ["", f"Largest χ² over all sets and lateral/bias tests: {top[0]:.3g} (log10 ND,eff = {float(top[1]):.2f}, {top[2]})."
                  + cmp_txt]
        summary_nl = nl
    else:
        summary_nl = {}
    lines += ["", "Consistency with earlier runs: calibration rows vs results/e3/rows.npz max relative difference "
              + (f"{checks['S_base_max_rel_diff']:.1e}" if "S_base_max_rel_diff" in checks else "n/a")
              + "; σ_pred/σ RMS here vs E14: "
              + ", ".join(f"{k} {a:.3g} vs {b:.3g}" for k, (a, b) in checks["vs_e14_rms_ratio"].items()) + "."]
    summary = {"sigma_dec": SIGMA, "calibration": list(CAL), "sd_calibration_only_scaled": sd0, "tests": res,
               "two_curve_tests": pairs, "groups": groups, "verdict": verdict, "checks": checks, "meta": meta,
               "nonlinear": summary_nl}
    dump_json(summary, os.path.join(out, "summary.json"))
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if HAVE_MPL:
        figure(res, meta, out)


def figure(res, meta, out):
    """Row 1: pointwise σ_pred/σ (RMS and max over the curve). Row 2: whole-curve factors λ_max and q."""
    import matplotlib.pyplot as plt
    INK2, GRID, BLUE, ORANGE = "#52514e", "#d9d8d4", "#2a78d6", "#eb6834"
    fig, ax = plt.subplots(2, 3, figsize=(11.5, 6.6), sharey="row")

    def style(a, ylab=None):
        a.grid(True, color=GRID, lw=0.5)
        for s_ in ("top", "right"):
            a.spines[s_].set_visible(False)
        a.axhline(1.0, color="k", ls=":", lw=1)
        a.set_xscale("log")
        if ylab:
            a.set_ylabel(ylab)
        if a.get_legend_handles_labels()[0]:
            a.legend(fontsize=6.5, loc="upper left")

    def series(prefix, vd=None):
        """(x, keys) of the tests of one family, sorted by x."""
        if prefix == "L1":
            ks = sorted((k for k in res if k.startswith("L1_VD")), key=lambda k: float(k.split("_VD")[1]))
            return [float(k.split("_VD")[1]) for k in ks], ks
        ks = [k for k in res if k.startswith(prefix) and k.endswith(f"_VD{vd:g}") and not k.startswith("L1_")] if prefix == "L" \
            else [k for k in res if k.startswith(prefix) and k.endswith(f"_VD{vd:g}")]
        if prefix == "L":
            ks = ks + [key("L1", vd)] if key("L1", vd) in res else ks
            xs = [float(k.split("_VD")[0][1:]) for k in ks]
        else:
            t0 = 7.0 if prefix == "tox" else 3.5
            ks = ks + [key("L1", vd)]
            xs = [float(k.split("_VD")[0][len(prefix):]) / t0 if not k.startswith("L1_") else 1.0 for k in ks]
        o = np.argsort(xs)
        return list(np.array(xs)[o]), [ks[i] for i in o]

    for row, (m1, m2, lab1, lab2, ylab) in enumerate((("rms_ratio", "max_ratio", "RMS over the curve", "max over the curve",
                                                        "σ_pred / σ (pointwise)"),
                                                       ("lam_max", "q_profile", "λ_max (any direction)", "q (ND,eff profile)",
                                                        "test χ² per unit calibration Δχ²"))):
        a = ax[row, 0]
        x, ks = series("L1")
        a.plot(x, [res[k][m1] for k in ks], "o-", color=INK2, label=lab1)
        a.plot(x, [res[k][m2] for k in ks], "s--", color=INK2, mfc="none", label=lab2)
        for v in (0.05, 1.0):
            a.axvline(v, color=GRID, lw=3, zorder=0)
        a.set_xlabel("VD (V), L = 1 µm")
        style(a, ylab)
        a = ax[row, 1]
        for vd, mk in ((0.05, "o"), (1.0, "^")):
            x, ks = series("L", vd)
            if len(x) > 1:
                a.plot(x, [res[k][m1] for k in ks], mk + "-", color=INK2, label=f"{lab1.split(' (')[0]}, VD {vd:g} V")
                a.plot(x, [res[k][m2] for k in ks], mk + "--", color=INK2, mfc="none", label=f"{lab2.split(' (')[0]}, VD {vd:g} V")
        a.axvline(1.0, color=GRID, lw=3, zorder=0)
        a.set_xlabel("L (µm)")
        style(a)
        a = ax[row, 2]
        for prefix, lab, col in (("tox", "t_ox", BLUE), ("tiwo", "t_IWO", ORANGE)):
            for vd, mk in ((0.05, "o"), (1.0, "^")):
                x, ks = series(prefix, vd)
                if len(x) > 1:
                    a.plot(x, [res[k][m1] for k in ks], mk + "-", color=col, label=f"{lab}, VD {vd:g} V")
        a.axvline(1.0, color=GRID, lw=3, zorder=0)
        a.set_xlabel("thickness / calibration thickness")
        style(a)
    ax[0, 0].set_title("(a) drain bias", loc="left")
    ax[0, 1].set_title("(b) channel length", loc="left")
    ax[0, 2].set_title("(c) vertical stack (RMS)", loc="left")
    ax[1, 0].set_title("(d) drain bias: whole curve", loc="left")
    ax[1, 1].set_title("(e) channel length: whole curve", loc="left")
    ax[1, 2].set_title("(f) vertical stack: λ_max", loc="left")
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_ylim(0.05, 5000)
    ax[1, 0].set_yscale("log")
    ax[1, 0].set_ylim(0.1, 1e8)
    k05 = key("L0.5", 1.0)
    if k05 in res and res[k05]["max_ratio"] > 1 and "vg_of_max" in res[k05]:
        ax[0, 1].annotate(f"first point above 10 pA (VG {res[k05]['vg_of_max']:+.2f} V)", (0.5, res[k05]["max_ratio"]),
                          (0.6, 40), fontsize=6.3, color=INK2, arrowprops={"arrowstyle": "->", "color": INK2, "lw": 0.7})
    fig.text(0.01, 0.005, "grey bands: calibration conditions (L = 1 µm, VD = 0.05/1 V, t_ox = 7 nm, t_IWO = 3.5 nm; thickness ratio 1). "
             "Linear analysis at the twin truth (Step 5), σ = 0.02 dec, Qf fixed.", fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(os.path.join(out, "e19_design_rule.png"), dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", choices=["lateral", "vertical"])
    ap.add_argument("--only", default="", help="comma-separated device labels to (re)run, e.g. L10,tiwo10")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--nonlinear", action="store_true", help="E5 sets vs truth at every lateral/bias test")
    ap.add_argument("--sets", default="results/e5_nd/summary.json")
    ap.add_argument("--values", default=None, help="log10 ND,eff of the E5 sets to run (default: all with χ² < 3.84)")
    ap.add_argument("--tag", default="all")
    ap.add_argument("--vgmin", type=float, default=-1.5)
    ap.add_argument("--out", default="results/e19")
    args = ap.parse_args()
    out = outdir(args.out)
    if args.group:
        run_group(args.group, out, vg_grid(args.vgmin, 2.5, 0.05), [d for d in args.only.split(",") if d] or None)
    if args.nonlinear:
        run_nonlinear(out, args.sets, args.values, args.tag)
    if args.report:
        report(out)


if __name__ == "__main__":
    main()
