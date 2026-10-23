# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 17.35}; start: results/e12_nd/log10_nd_eff_17.85/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- χ² = 1633.4 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1100 dec; σ_eff = 0.1192 dec
- J (paper weights, our a_i) = 0.351

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.349 | 0.0595 |
| wf_gate | 4.054 | 0.712 |
| log10_mu0 | 1.412 | 0.122 |
| log10_n_gd | 19.8 | 2.82 |
| e_gd | 3.203e-06 | 2.62 |
| s_gd | 0.184 | 1.33 |
| log10_n_ta | 19.27 | 36.1 |
| w_ta | 0.03 | 0.935 |
| log10_nd_eff | 17.35 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.525 | +0.6 | +23.7 | +0.013 |
| 1.0 | 0.214 | +38.7 | +21.1 | -0.020 |
| L1_VD0.1 | 0.458 | -31.3 | +6.9 | -0.002 |
| L1_VD0.2 | 0.084 | -38.3 | -1.7 | +0.012 |
| L1_VD0.5 | 0.574 | +48.8 | -22.3 | -0.008 |
| L2_VD0.05 | 0.629 | +8.8 | -28.5 | +0.109 |
| L2_VD1 | 0.664 | +93.6 | +nan | +0.081 |
