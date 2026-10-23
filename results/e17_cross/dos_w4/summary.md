# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 16.8055, 'log10_n_gd': 19.6435, 'e_gd': 0.08, 's_gd': 0.224, 'log10_n_ta': 19.8921, 'w_ta': 0.238}; start: results/e11_2W_v2/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 2660.9 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1404 dec; σ_eff = 0.1445 dec
- J (paper weights, our a_i) = 0.701

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.342 | 0.0799 |
| wf_gate | 3.7 | 0.0259 |
| log10_mu0 | 1.409 | 0.141 |
| log10_nd_eff | 16.81 (fixed) | – |
| log10_n_gd | 19.64 (fixed) | – |
| e_gd | 0.08 (fixed) | – |
| s_gd | 0.224 (fixed) | – |
| log10_n_ta | 19.89 (fixed) | – |
| w_ta | 0.238 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.301 | +28.3 | +59.7 | +0.008 |
| 1.0 | 0.245 | +65.2 | +67.8 | -0.019 |
| L1_VD0.1 | 0.222 | -42.5 | +36.8 | -0.007 |
| L1_VD0.2 | 0.261 | -82.1 | +nan | +0.007 |
| L1_VD0.5 | 0.123 | -53.1 | +17.6 | -0.010 |
| L2_VD0.05 | 0.276 | +20.8 | +3.2 | +0.105 |
| L2_VD1 | 0.307 | +103.5 | +nan | +0.080 |
