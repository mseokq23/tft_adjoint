# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'w_ta': 0.238}; start: results/e11_2W_v2/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1621.9 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1096 dec; σ_eff = 0.1188 dec
- J (paper weights, our a_i) = 0.369

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.348 | 0.0559 |
| wf_gate | 4.188 | 0.533 |
| log10_nd_eff | 18.95 | 60.8 |
| log10_mu0 | 1.414 | 0.114 |
| log10_n_gd | 19.66 | 28.8 |
| e_gd | 0.1277 | 3.93 |
| s_gd | 0.09514 | 5.53 |
| log10_n_ta | 19.16 | 161 |
| w_ta | 0.238 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.556 | -1.8 | +22.9 | +0.017 |
| 1.0 | 0.235 | +39.9 | +21.0 | -0.016 |
| L1_VD0.1 | 0.501 | -34.0 | +5.4 | +0.003 |
| L1_VD0.2 | 0.082 | -38.3 | -2.3 | +0.016 |
| L1_VD0.5 | 0.631 | +56.1 | -24.9 | -0.003 |
| L2_VD0.05 | 0.670 | +3.1 | -32.2 | +0.111 |
| L2_VD1 | 0.702 | +89.9 | +nan | +0.083 |
