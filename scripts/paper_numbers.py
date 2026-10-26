"""Key numbers for the abstract and the paper, each recomputed from the stored results with its source file.

No simulation. Writes docs/paper/key_numbers.md and docs/paper/key_numbers.json. Every number quoted in the
abstract or the paper should come from this table (or from the result file it names), so that the text can be
re-checked whenever a result changes.

Usage:  python scripts/paper_numbers.py [--out docs/paper]
"""
from __future__ import annotations

import argparse
import json
import math
import os

from _common import ROOT, dump_json, np  # noqa: F401
from e14_transfer_test_audit import measured_ss
from tftsim import PAPER_SETS, Geometry, Material

Q, EPS0 = 1.602176634e-19, 8.8541878128e-12


def load(p):
    return json.load(open(p))


def ub(x, nd):
    """Upper bound shown with nd decimals: rounded up, so that "≤ ub" stays true."""
    return f"{math.ceil(x * 10 ** nd - 1e-9) / 10 ** nd:.{nd}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/paper")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rows = []

    def add(group, key, text, value, source, note=""):
        rows.append({"group": group, "key": key, "text": text, "value": value, "source": source, "note": note})

    # A. structural degeneracies
    ev = load("results/e3/summary.json")["base"]["eigvals"]
    add("A structural", "qf_wf_eig_ratio", f"Qf ≡ gate WF: Fisher eigenvalue ratio λmin/λmax = {min(ev) / max(ev):.1e}",
        min(ev) / max(ev), "results/e3/summary.json (base.eigvals)", "single gate-oxide thickness, Step-5 point")
    geo, mat = Geometry(), Material()
    cox = mat.eps_hfo2 * EPS0 / (geo.t_gox_nm * 1e-9)
    add("A structural", "qf_per_1e11_mV", f"Qf 1e11 cm⁻² ≙ gate WF {Q * 1e15 / cox * 1e3:.1f} mV", Q * 1e15 / cox * 1e3,
        "tftsim/params.py (t_ox, ε_HfO2)")
    e8, e8t = load("results/e8/summary.json"), load("results/e8_thickness/summary.json")
    add("A structural", "donor_acceptor_iv", f"donor ≡ acceptor + ND,eff: max |Δlog10 ID| {max(e8['cases'].values()):.1e} dec "
        f"(7 curves), thickness variants {max(e8t['cases'].values()):.1e} dec", max(max(e8["cases"].values()), max(e8t["cases"].values())),
        "results/e8/summary.json, results/e8_thickness/summary.json (cases)")
    add("A structural", "donor_acceptor_cv", f"donor ≡ acceptor + ND,eff: max C–V difference {e8['cv_max_rel']:.1e} (rel.), "
        f"thickness variants {e8t['cv_max_rel']:.1e}", max(e8["cv_max_rel"], e8t["cv_max_rel"]), "results/e8*/summary.json (cv_max_rel)")
    s5, w4 = PAPER_SETS["step5"], PAPER_SETS["w4"]
    add("A structural", "paper_w4_step2", f"paper 4%-W recalibration moved gate WF {s5.wf_gate:.2f} → {w4.wf_gate:.2f} eV and Qf "
        f"{s5.qf:+.2e} → {w4.qf:+.2e} cm⁻² together (exactly degenerate with one t_ox)", [s5.wf_gate, w4.wf_gate, s5.qf, w4.qf],
        "tftsim/params.py PAPER_SETS (from the reference paper, Sec. IV)")

    # B. synthetic twin
    e5 = load("results/e5_nd/summary.json")
    rec = {round(r["value"], 2): r for r in e5["records"] if "chi2" in r}
    inside = [v for v in rec if 15.85 <= v <= 18.95]
    add("B twin", "twin_profile_max", f"twin profile: ND,eff {min(inside)}–{max(inside)} (×{10 ** (max(inside) - min(inside)):.0f}), "
        f"calibration χ² ≤ {ub(max(rec[v]['chi2'] for v in inside), 3)} (95 % threshold 3.84)",
        max(rec[v]["chi2"] for v in inside), "results/e5_nd/summary.json (records)")
    add("B twin", "twin_profile_rejected", f"twin profile: χ² {rec[19.15]['chi2']:.2f} at 10^19.15, {rec[19.4]['chi2']:.1f} at 10^19.4",
        [rec[19.15]["chi2"], rec[19.4]["chi2"]], "results/e5_nd/summary.json")
    e14 = load("results/e14/summary.json")
    syn = {k: v for k, v in e14["synthetic"].items() if float(k) <= 18.95}
    tmax = max(max(v["drain"]["vth_max"], v["length"]["vth_max"]) for v in syn.values())
    add("B twin", "twin_transfer_dvth", f"twin sets reproduce the paper's 5 transfer-test curves within |ΔVth| ≤ {ub(tmax, 2)} mV",
        tmax, "results/e14/summary.json (synthetic)")
    r194 = e14["synthetic"]["19.40"]
    add("B twin", "twin_rejected_set_transfer", f"the calibration-rejected 10^19.4 set still passes: max |ΔVth| "
        f"{max(r194['drain']['vth_max'], r194['length']['vth_max']):.1f} mV", max(r194["drain"]["vth_max"], r194["length"]["vth_max"]),
        "results/e14/summary.json (synthetic 19.40)")
    pap = e14["paper"]
    add("B twin", "paper_errors", f"paper's transfer errors: drain avg/max {pap['drain']['vth_avg']}/{pap['drain']['vth_max']} mV, "
        f"length {pap['length']['vth_avg']}/{pap['length']['vth_max']} mV", pap, "results/e14/summary.json (paper; reference paper Sec. III)")
    e18 = load("results/e18/summary.json")
    tj = [r["chi2"] for r in e18["twin_joint"]["rows"] if r.get("kind", "fixed") != "truth" and r["chi2"] > 0]
    add("B twin", "twin_joint_max", f"twin, five VD curves calibrated together: χ² ≤ {ub(max(tj), 3)} over ×"
        f"{e18['twin_joint_range']['factor']:.0f}", max(tj), "results/e18/summary.json (twin_joint)")
    th = {round(r["value"], 2): r["theta"] for r in e5["records"] if "theta" in r}
    dwf = (th[18.9]["wf_gate"] - th[15.9]["wf_gate"]) * 1e3
    dn = (10 ** 18.9 - 10 ** 15.9) * 1e6
    est = Q * dn * geo.t_iwo_nm * 1e-9 / cox * 1e3
    add("B twin", "mechanism", f"ND,eff 10^15.9 → 10^18.9: gate WF moves {dwf:+.0f} mV vs q·ΔN·t_IWO/C_ox = {est:+.0f} mV",
        [dwf, est], "results/e5_nd/summary.json (records.theta); tftsim/params.py")

    # C. measured (digitised)
    e12 = load("results/e12_nd/summary.json")
    v3 = load("results/e11_2W_v3/summary.json")
    d12 = [r["dchi2_eff"] for r in e12["rows"]]
    vals = [r["value"] for r in e12["rows"]]
    add("C measured", "meas_profile_max", f"measured profile: ND,eff 10^{min(vals):.2f}–10^{max(vals):.2f} (×{10 ** (max(vals) - min(vals)):.0f}), "
        f"Δχ²_eff ≤ {ub(max(d12), 2)} (max {max(d12):.3f}; χ²_min {v3['chi2']:.1f}, σ_eff {v3['sigma_eff']:.4f} dec)", max(d12), "results/e12_nd/summary.json")
    sp = e14["spread"]
    tr = {k: sp[k]["measured_mV"] for k in ("L1_VD0.1", "L1_VD0.2", "L1_VD0.5", "L2_VD0.05", "L2_VD1")}
    add("C measured", "meas_transfer_spread", "predicted transfer-test Vth spread over the 7 sets (full curves): "
        + ", ".join(f"{k} {v:.1f}" for k, v in tr.items()) + " mV", tr, "results/e14/summary.json (spread.*.measured_mV)")
    ss = measured_ss()
    me = [v3["sigma_eff"] * ss[vd] for vd in (0.05, 1.0)]
    add("C measured", "model_error_vth", f"model error σ_eff in Vth: {min(me):.1f}–{max(me):.1f} mV (σ_eff × measured SS "
        f"{ss[0.05]:.0f}/{ss[1.0]:.0f} mV/dec)", me, "results/e11_2W_v3/summary.json; data/meas.csv")
    e13 = load("results/e13/summary.json")
    tt = {k: e13["tests"][k]["vth_spread_mV"] for k in ("t_ox 14 nm VD=1", "t_IWO 7 nm VD=1")}
    ends = {c: abs(e13["rows"]["19.050"][c]["1"]["dVth_mV"] - e13["rows"]["16.850"][c]["1"]["dVth_mV"]) for c in ("t_ox 14 nm", "t_IWO 7 nm")}
    add("C measured", "meas_thickness_spread", f"same 7 sets, VD 1 V: t_ox 14 nm {tt['t_ox 14 nm VD=1']:.0f} mV, t_IWO 7 nm "
        f"{tt['t_IWO 7 nm VD=1']:.0f} mV; ND,eff ends only {ends['t_ox 14 nm']:.0f}/{ends['t_IWO 7 nm']:.0f} mV",
        {"all": tt, "ends": ends}, "results/e13/summary.json", "the largest differences are between v2 and v3 (different DOS)")
    mj = [r["dchi2_eff"] for r in e18["measured_joint"]["rows"]]
    add("C measured", "meas_joint_max", f"measured, five VD curves together: Δχ²_eff ≤ {ub(max(mj), 2)}; L 2 µm Vth spread "
        f"{e18['vth_spread_full_mV']['joint']['L2_VD0.05']:.1f}/{e18['vth_spread_full_mV']['joint']['L2_VD1']:.1f} mV",
        max(mj), "results/e18/summary.json (measured_joint, vth_spread_full_mV)")
    rob = {r["check"]: [] for r in e14["robustness"]["rows"]}
    for r in e14["robustness"]["rows"]:
        rob[r["check"]].append(round(r["dchi2_eff"], 2))
    add("C measured", "robustness", "robustness Δχ²_eff: " + "; ".join(f"{k}: {v}" for k, v in rob.items()), rob,
        "results/e14/summary.json (robustness)")
    e21 = load("results/e21/summary.json")
    eq = e21["multistart_equivalent"]
    ms_max = max(r["dchi2_eff"] for r in e21["multistart"])
    dsp, vsp = eq["descriptor_spread"], eq["vth_spread_mV"]
    tests_only = {k: v for k, v in vsp.items() if not k.startswith("cal")}
    add("C measured", "multistart", f"{len(e21['multistart'])} free fits from different starts: Δχ²_eff ≤ {ub(ms_max, 2)}; N_TA spread "
        f"{dsp['log10_n_ta']:.2f} dec (×{10 ** dsp['log10_n_ta']:.0f}), gate WF {dsp['wf_gate'] * 1e3:.0f} mV, σ_GD "
        f"{dsp['s_gd']:.3f} eV; transfer-test Vth spread {min(tests_only.values()):.1f}–{max(tests_only.values()):.1f} mV",
        {"dchi2_eff_max": ms_max, "spread": dsp, "vth_spread_mV": vsp}, "results/e21/summary.json", "all stopped at the evaluation limit")
    eps = e21["eps15"]
    add("C measured", "eps15", "ε_ch = 15: " + ", ".join(f"10^{r['log10_nd_eff']:.2f} χ² {r['chi2']:.1f} (ε 10: {r['chi2_eps10']:.1f})"
        for r in eps) + f"; ends Δχ²_eff {max(r['dchi2_eff_vs_better_end'] for r in eps):.2f}", eps and
        [r["chi2"] for r in eps], "results/e21/summary.json (eps15)")

    # D. design rule and identifiability
    e19 = load("results/e19/summary.json")
    tests = {k: v["rms_ratio"] for k, v in e19["tests"].items() if k not in e19["calibration"]}
    grp = {"VD (L 1 µm)": [v for k, v in tests.items() if k.startswith("L1_")],
           "L": [v for k, v in tests.items() if k.startswith("L") and not k.startswith("L1_")],
           "t_ox": [v for k, v in tests.items() if k.startswith("tox")], "t_IWO": [v for k, v in tests.items() if k.startswith("tiwo")]}
    add("D design rule", "e19_ratio", f"{len(tests)} tests, predictive σ / noise (curve RMS): "
        + ", ".join(f"{k} {min(v):.2f}–{max(v):.2f}" if max(v) < 1 else f"{k} {min(v):.1f}–{max(v):.0f}" for k, v in grp.items()),
        {k: [min(v), max(v)] for k, v in grp.items()}, "results/e19/summary.json (tests.*.rms_ratio)")
    lin = e14["linear"]
    pap_t = [v["rms_ratio"] for k, v in lin.items() if "(paper)" in k]
    rec_t = [v["rms_ratio"] for k, v in lin.items() if "(paper)" not in k]
    add("D design rule", "e14_linear", f"paper's transfer tests {min(pap_t):.2f}–{max(pap_t):.2f}, recommended tests "
        f"{min(rec_t):.0f}–{max(rec_t):.0f}", [min(pap_t), max(pap_t), min(rec_t), max(rec_t)], "results/e14/summary.json (linear)")
    e20 = load("results/e20/summary.json")["measured"]
    add("D design rule", "e20_nd", f"measured optimum (local, linear): ND,eff 1σ calibration only {e20['joint (5 VD, E18)']['log10_nd_eff']:.2f} dec "
        f"(2 curves {e20['pair (VD 0.05/1 V)']['log10_nd_eff']:.2f}); + C–V 0.5 % of the same device only: "
        f"{e20['joint + C–V (0.5 %)']['log10_nd_eff']:.3f} dec (2 %: {e20['joint + C–V (2 %)']['log10_nd_eff']:.3f}); "
        "+ t_ox 14 nm, t_IWO 7 nm, C–V 0.5 %: "
        f"{e20['joint + all three']['log10_nd_eff']:.3f} dec (C–V 2 %: {e20['joint + all three (C–V 2 %)']['log10_nd_eff']:.3f})",
        {k: e20[k]["log10_nd_eff"] for k in e20}, "results/e20/summary.json (measured)", "local values; strongly point dependent")

    # E. the authors' open demo (E22, analytic model, synthetic data)
    e22 = load("results/e22/summary.json")
    dr = e22["default_run"]
    add("E open demo", "demo_default_run", f"demo default run reproduced: J {dr['J_initial']:.2f} → {dr['J_final']:.3f}, held-out max "
        f"|ΔVth| {dr['heldout_max_abs_dVth_mV']['heldout_drain']:.1f}/{dr['heldout_max_abs_dVth_mV']['heldout_length']:.1f} mV "
        "(drain/length)", dr["J_final"], "results/e22/summary.json (default_run)", "demo commit ccaa678")
    a = e22["A_null_phiG_Qf"]
    sp = e22["B_valley_prediction_spread"]
    add("E open demo", "demo_null_valley", f"φG–Qf null ({a['dphiG_per_1e12_Qf_mV']:.1f} mV per 1e12 cm⁻²): Qf "
        f"{a['valley_Q_f_cm2'][0] / 1e12:.2f}–{a['valley_Q_f_cm2'][1] / 1e12:.2f} ×1e12, φG {a['valley_phi_G_eV'][0]:.2f}–"
        f"{a['valley_phi_G_eV'][1]:.2f} eV give the same J; Vth spread VD/L tests "
        f"{max(sp[k]['vth_spread_mV'] for k in sp if k.startswith(('VD', 'L '))):.3f} mV, t_ox 14 nm {sp['t_ox 14 nm']['vth_spread_mV']:.0f} mV",
        a["valley_Q_f_cm2"], "results/e22/summary.json (A_null_phiG_Qf, B_valley_prediction_spread)")
    for tag, key in (("demo set", "demo_endpoints"), ("+ t_ox 14 nm curve", "demo_endpoints_tox")):
        rs = [x for x in e22["D_controller_runs"] if x["calibration"] == tag]
        ho = [max(x["heldout_drain_max_dVth_mV"], x["heldout_length_max_dVth_mV"]) for x in rs]
        t20 = [x["t_ox 20 nm dVth vs truth (mV)"] for x in rs]
        qf = [x["Q_f"] / 1e12 for x in rs]
        add("E open demo", key, f"{tag}, {len(rs)} controller runs: held-out max |ΔVth| {min(ho):.2f}–{max(ho):.2f} mV, end-point Qf "
            f"{min(qf):.4f}–{max(qf):.4f} ×1e12 (truth 6), t_ox 20 nm Vth error {min(t20):.2f} to {max(t20):.2f} mV",
            {"heldout": [min(ho), max(ho)], "Qf_1e12": [min(qf), max(qf)], "tox20": [min(t20), max(t20)]},
            "results/e22/summary.json (D_controller_runs)", "trust radius 0.3–0.7 × three starting gate WFs")
    ec = e22["E_cit_vs_tail"]
    add("E open demo", "demo_cit_tail", f"interface term C_it unknown: N_tail {ec['valley_N_tail'][0]:.0e}–{ec['valley_N_tail'][1]:.2e} "
        f"equivalent at one t_ch; SS spread t_ch 7 nm {ec['prediction_spread']['t_ch 7 nm']['ss_spread_mV_dec']:.1f} mV/dec, "
        f"t_ox splits {max(ec['prediction_spread'][k]['ss_spread_mV_dec'] for k in ec['prediction_spread'] if k.startswith('t_ox')):.3f}",
        ec["valley_N_tail"], "results/e22/summary.json (E_cit_vs_tail)", "extension: the demo fixes C_it/C_ox = 0.10")

    lines =["# Key numbers (recomputed from the stored results)", "",
             "Generated by `scripts/paper_numbers.py` (no simulation). Quote numbers from this table or from the file it names.",
             "Bounds written as ≤ are rounded up at the shown precision.", ""]
    group = None
    for r in rows:
        if r["group"] != group:
            group = r["group"]
            lines += ["", f"## {group}", "", "| key | number | source | note |", "|---|---|---|---|"]
        cell = r["text"].replace("|", "\\|")
        lines.append(f"| {r['key']} | {cell} | {r['source']} | {r['note']} |")
    with open(os.path.join(args.out, "key_numbers.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    dump_json(rows, os.path.join(args.out, "key_numbers.json"))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
