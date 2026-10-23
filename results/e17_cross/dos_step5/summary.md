# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 18.4082, 'log10_n_gd': 19.3404, 'e_gd': 0.12, 's_gd': 0.08, 'log10_n_ta': 19.8921, 'w_ta': 0.133}; start: results/e11_2W_v2/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 2065.6 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1237 dec; σ_eff = 0.1273 dec
- J (paper weights, our a_i) = 0.405

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 3.952 | 47.2 |
| wf_gate | 3.808 | 0.00871 |
| log10_mu0 | 1.246 | 0.108 |
| log10_nd_eff | 18.41 (fixed) | – |
| log10_n_gd | 19.34 (fixed) | – |
| e_gd | 0.12 (fixed) | – |
| s_gd | 0.08 (fixed) | – |
| log10_n_ta | 19.89 (fixed) | – |
| w_ta | 0.133 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.588 | -14.0 | +13.9 | -0.012 |
| 1.0 | 0.285 | +36.5 | +14.2 | -0.029 |
| L1_VD0.1 | 0.555 | -41.9 | -1.6 | -0.025 |
| L1_VD0.2 | 0.055 | -24.3 | -2.8 | -0.012 |
| L1_VD0.5 | 0.719 | +68.5 | -32.0 | -0.026 |
| L2_VD0.05 | 0.708 | -4.8 | -38.0 | +0.028 |
| L2_VD1 | 0.743 | +88.5 | +nan | +0.004 |
