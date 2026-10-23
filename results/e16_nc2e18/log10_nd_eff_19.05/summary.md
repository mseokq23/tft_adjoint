# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 19.05}; start: results/e12_nd/log10_nd_eff_19.05/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 2e+18 cm⁻³, ε_ch = 10 (changed)
- χ² = 1598.6 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1088 dec; σ_eff = 0.1179 dec
- J (paper weights, our a_i) = 0.395

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.325 | 0.0491 |
| wf_gate | 4.143 | 0.239 |
| log10_mu0 | 1.417 | 0.09 |
| log10_n_gd | 20.03 | 5 |
| e_gd | 0.1626 | 0.0389 |
| s_gd | 0.02 | 0.259 |
| log10_n_ta | 20.12 | 0.522 |
| w_ta | 0.05005 | 0.067 |
| log10_nd_eff | 19.05 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.640 | -4.3 | +22.0 | +0.018 |
| 1.0 | 0.307 | +39.4 | +19.1 | -0.016 |
| L1_VD0.1 | 0.624 | -34.4 | +3.7 | +0.004 |
| L1_VD0.2 | 0.084 | -40.6 | -5.1 | +0.017 |
| L1_VD0.5 | 0.801 | +78.3 | -34.0 | -0.002 |
| L2_VD0.05 | 0.779 | +4.7 | -36.1 | +0.113 |
| L2_VD1 | 0.808 | +90.0 | +nan | +0.084 |
