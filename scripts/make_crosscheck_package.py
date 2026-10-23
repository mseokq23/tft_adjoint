"""Cross-check package: the alternative 2 %-W descriptor sets as a table that another
simulator (e.g. the Sentaurus deck of arXiv:2609.12184) can run directly.

No simulation is run here; the table is assembled from stored results:
  results/e5_nd/summary.json       ND,eff profile on the Step-5 twin (re-fitted sets,
                                   calibration chi2/RMS/max, hold-out errors)
  results/e9a, results/e9b         nonlinear check of the discriminating measurements
                                   (t_IWO = 7 nm, t_ox = 14 nm, quasi-static C-V)
  results/e8/summary.json          donor -> acceptor + ND,eff identity (set F)

Rows
  Step5_reference          the paper's Step-5 values (twin truth)
  A-E                      fixed ND,eff = 10^15.9, 16.9, 17.4, 17.9, 18.9 cm^-3, other
                           descriptors re-fitted (Qf fixed, exactly degenerate with gate WF)
  G_rejected_by_calibration  ND,eff = 10^19.4: control rejected by the calibration curves
  F_acceptor_type_Gaussian   exact identity (analytic), I-V/C-V unchanged

Columns starting with devsim_ are results of the DEVSIM model relative to Step 5:
chi2 over the 64 calibration points (VD = 0.05/1 V, ID >= 10 pA, sigma = 0.02 dec),
RMS / max |dlog10 ID| on those points, worst RMS and max |dVth| over the five transfer
curves, dVth (VD = 1 V, 100 nA) for t_IWO = 7 nm and t_ox = 14 nm, max |dC| of the
quasi-static gate C-V in % of C_ox*A.

Usage:  python scripts/make_crosscheck_package.py [--out results/crosscheck]
        python scripts/make_crosscheck_package.py --stack_vth     (also simulates, about 3 min)
Output: <out>/alternative_parameter_sets.csv  (the README.txt next to it is hand-written)
        <out>/step5_stack_vth.json  with --stack_vth: Step-5 Vth (100 nA, VD = 1 V) of the nominal
        stack and of the two thicker stacks, quoted in README.txt for choosing the sweep range
"""
from __future__ import annotations

import argparse
import csv
import json

from _common import ROOT, dump_json, np, outdir, vg_grid
from tftsim import PAPER_SETS

LABELS = {15.9: "A", 16.9: "B", 17.4: "C", 17.9: "D", 18.9: "E", 19.4: "G_rejected_by_calibration"}
HEADER = ["set", "ND_eff_cm-3", "WF_SD_eV", "WF_gate_eV", "Qf_cm-2", "mu0_cm2_Vs", "N_GD_peak_cm-3_eV-1",
          "E_GD_below_Ec_eV", "sigma_GD_eV", "Gaussian_type", "N_TA_edge_cm-3_eV-1", "W_TA_eV",
          "devsim_chi2_calibration_sigma0.02dec", "devsim_RMS_dev_calibration_dec",
          "devsim_max_dev_calibration_dec", "devsim_worst_RMS_dev_transfer_dec",
          "devsim_max_abs_dVth_transfer_mV", "devsim_dVth_tIWO7nm_VD1_mV", "devsim_dVth_tox14nm_VD1_mV",
          "devsim_max_CV_dev_pct_Cox"]


def g4(x):
    return f"{x:.4g}"


def f4(x):
    return f"{x:.4f}"


def descriptor_cols(nd, wf_sd, wf_gate, qf, mu0, n_gd, e_gd, s_gd, kind, n_ta, w_ta):
    return [g4(nd), f4(wf_sd), f4(wf_gate), g4(qf), g4(mu0), g4(n_gd), f4(e_gd), f4(s_gd), kind, g4(n_ta), f4(w_ta)]


def load(path):
    with open(path if path.startswith("/") else f"{ROOT}/{path}") as f:
        return json.load(f)


def stack_vth(path):
    """Step-5 Vth (100 nA, VD = 1 V) for the nominal stack, t_IWO = 7 nm and t_ox = 14 nm (DEVSIM)."""
    from tftsim import TFT, Geometry
    from tftsim.metrics import extract
    vgs = vg_grid(-1.0, 2.5, 0.05)
    res = {}
    for tag, g in (("nominal", Geometry()), ("t_IWO 7 nm", Geometry(t_iwo_nm=7.0)), ("t_ox 14 nm", Geometry(t_gox_nm=14.0))):
        t = TFT(PAPER_SETS["step5"], g).build()
        recs = sorted(t.transfer(1.0, vgs), key=lambda r: r["vg"])
        i = np.array([r["id"] for r in recs])
        res[tag] = {"vth_V": extract(vgs, i)["vth"], "vg_min_V": float(vgs[0]), "id_at_vg_min_A": float(i[0])}
        print(f"{tag}: Vth = {res[tag]['vth_V']:+.3f} V, ID(VG = {vgs[0]:g} V) = {i[0]:.2e} A", flush=True)
    dump_json(res, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="results/e5_nd/summary.json")
    ap.add_argument("--e9", default="results/e9a/summary.json,results/e9b/summary.json")
    ap.add_argument("--e8", default="results/e8/summary.json")
    ap.add_argument("--out", default="results/crosscheck")
    ap.add_argument("--stack_vth", action="store_true", help="also simulate the Step-5 Vth of the thicker stacks")
    args = ap.parse_args()
    out = outdir(args.out)
    if args.stack_vth:
        stack_vth(f"{out}/step5_stack_vth.json")

    p5 = PAPER_SETS["step5"]
    rows = [["Step5_reference"] + descriptor_cols(p5.nd_eff, p5.wf_sd, p5.wf_gate, p5.qf, p5.mu0, p5.n_gd, p5.e_gd,
                                                  p5.s_gd, "donor", p5.n_ta, p5.w_ta)
            + ["0", f"{0:.2e}", f"{0:.2e}", f"{0:.2e}", f"{0:.2f}", "0", "0", "0"]]

    e9 = {}
    for p in args.e9.split(","):
        e9.update({round(float(k), 4): v for k, v in load(p)["results"].items()})

    for r in load(args.profile)["records"]:
        key = round(r["value"], 4)
        if key not in LABELS:
            continue
        th = r["theta"]
        ho = r["holdout"].values()
        row = [LABELS[key]] + descriptor_cols(10 ** th["log10_nd_eff"], th["wf_sd"], th["wf_gate"], p5.qf,
                                              10 ** th["log10_mu0"], 10 ** th["log10_n_gd"], th["e_gd"], th["s_gd"],
                                              "donor", 10 ** th["log10_n_ta"], th["w_ta"])
        row += [f"{r['chi2']:.3g}", f"{r['rms_dec']:.2e}", f"{r['max_dec']:.2e}",
                f"{max(q['rms_dec'] for q in ho):.2e}", f"{max(abs(q['dVth_mV']) for q in ho):.2f}"]
        if key in e9:
            x = e9[key]
            row += [f"{x['t_IWO 7 nm']['dVth_mV']:+.1f}", f"{x['t_ox 14 nm']['dVth_mV']:+.1f}",
                    f"{x['C-V']['max_dev_pct']:.2f}"]
        else:  # not evaluated: rejected by the calibration curves already
            row += ["n/a", "n/a", "n/a"]
        rows.append(row)

    e8 = load(args.e8)
    rows.append(["F_acceptor_type_Gaussian"]
                + descriptor_cols(e8["nd_eff_flipped"], p5.wf_sd, p5.wf_gate, p5.qf, p5.mu0, p5.n_gd, p5.e_gd,
                                  p5.s_gd, "acceptor", p5.n_ta, p5.w_ta)
                + ["identical (analytic)", "<1e-14", f"{max(e8['cases'].values()):.1e}", "<1e-14", f"{0:.2f}",
                   "0 (identical)", "0 (identical)", "0 (identical)"])

    order = ["Step5_reference", "A", "B", "C", "D", "E", "G_rejected_by_calibration", "F_acceptor_type_Gaussian"]
    rows.sort(key=lambda r: order.index(r[0]))
    with open(f"{out}/alternative_parameter_sets.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(HEADER)
        w.writerows(rows)
    print(f"wrote {out}/alternative_parameter_sets.csv ({len(rows)} sets)")


if __name__ == "__main__":
    main()
