# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: Qf only (degenerate with WF); start: results/e21_multistart/start_w4.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1625.0 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1097 dec; σ_eff = 0.1202 dec
- J (paper weights, our a_i) = 0.356

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.347 | 0.0542 |
| wf_gate | 4.14 | 0.506 |
| log10_nd_eff | 18.53 | 2.36 |
| log10_mu0 | 1.414 | 0.101 |
| log10_n_gd | 19.74 | 3.06 |
| e_gd | 0.05208 | 4.14 |
| s_gd | 0.1675 | 1.98 |
| log10_n_ta | 18.02 | 599 |
| w_ta | 0.09189 | 119 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.517 | -1.9 | +23.2 | +0.017 |
| 1.0 | 0.209 | +38.7 | +21.7 | -0.016 |
| L1_VD0.1 | 0.447 | -32.2 | +7.0 | +0.002 |
| L1_VD0.2 | 0.087 | -37.5 | -0.9 | +0.016 |
| L1_VD0.5 | 0.558 | +47.0 | -21.1 | -0.003 |
| L2_VD0.05 | 0.620 | +6.1 | -29.1 | +0.111 |
| L2_VD1 | 0.654 | +91.9 | +nan | +0.083 |
