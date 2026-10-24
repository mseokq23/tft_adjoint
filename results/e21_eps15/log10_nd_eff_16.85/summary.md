# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 16.85}; start: results/e12_nd/log10_nd_eff_16.85/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 15 (changed)
- χ² = 1631.4 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1099 dec; σ_eff = 0.1191 dec
- J (paper weights, our a_i) = 0.353

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.334 | 0.0747 |
| wf_gate | 4.033 | 0.806 |
| log10_mu0 | 1.396 | 0.124 |
| log10_n_gd | 19.81 | 1.99 |
| e_gd | 0.04053 | 1.9 |
| s_gd | 0.165 | 1.09 |
| log10_n_ta | 19.68 | 12.8 |
| w_ta | 0.03 | 0.34 |
| log10_nd_eff | 16.85 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.523 | -1.2 | +23.4 | +0.018 |
| 1.0 | 0.212 | +38.6 | +21.5 | -0.011 |
| L1_VD0.1 | 0.454 | -32.8 | +6.7 | +0.004 |
| L1_VD0.2 | 0.086 | -37.5 | -1.0 | +0.018 |
| L1_VD0.5 | 0.568 | +48.3 | -21.4 | -0.000 |
| L2_VD0.05 | 0.626 | +7.5 | -28.9 | +0.105 |
| L2_VD1 | 0.662 | +93.2 | +nan | +0.078 |
