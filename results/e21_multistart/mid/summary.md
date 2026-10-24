# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: Qf only (degenerate with WF); start: results/e21_multistart/start_mid.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1624.8 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1097 dec; σ_eff = 0.1202 dec
- J (paper weights, our a_i) = 0.356

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.348 | 0.0646 |
| wf_gate | 4.137 | 0.756 |
| log10_nd_eff | 18.53 | 124 |
| log10_mu0 | 1.414 | 0.138 |
| log10_n_gd | 19.73 | 51.6 |
| e_gd | 0.05941 | 28.3 |
| s_gd | 0.1631 | 11.1 |
| log10_n_ta | 18.68 | 644 |
| w_ta | 0.04549 | 129 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.519 | -2.0 | +23.2 | +0.017 |
| 1.0 | 0.210 | +38.5 | +21.6 | -0.016 |
| L1_VD0.1 | 0.449 | -32.3 | +7.0 | +0.002 |
| L1_VD0.2 | 0.087 | -37.8 | -1.0 | +0.016 |
| L1_VD0.5 | 0.561 | +47.1 | -21.2 | -0.003 |
| L2_VD0.05 | 0.622 | +6.0 | -29.2 | +0.111 |
| L2_VD1 | 0.656 | +91.7 | +nan | +0.083 |
