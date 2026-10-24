# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/0.1/0.2/0.5/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 19.05}; start: results/e12_nd/log10_nd_eff_19.05/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1857.5 on 93 residuals (σ = 0.02 dec) → RMS misfit 0.0894 dec; σ_eff = 0.0935 dec
- J (paper weights, our a_i) = 0.381
- fitted points per VD: {'0.05': 24, '0.1': 17, '0.2': 15, '0.5': 9, '1': 28} (grid points at measured VG only, gaps ≤ 0.05 V)

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.347 | 0.0407 |
| wf_gate | 4.148 | 0.308 |
| log10_mu0 | 1.415 | 0.0747 |
| log10_n_gd | 19.92 | 2.54 |
| e_gd | 0.1637 | 0.0425 |
| s_gd | 0.02 | 0.15 |
| log10_n_ta | 20.02 | 0.529 |
| w_ta | 0.0806 | 0.067 |
| log10_nd_eff | 19.05 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.676 | -2.0 | +17.5 | +0.015 |
| 0.1 | 0.688 | -28.6 | +1.5 | +0.001 |
| 0.2 | 0.087 | -27.7 | -6.0 | +0.014 |
| 0.5 | 0.904 | +97.6 | -41.0 | -0.005 |
| 1.0 | 0.367 | +39.2 | +13.0 | -0.018 |
| L2_VD0.05 | 0.828 | +6.7 | -36.8 | +0.109 |
| L2_VD1 | 0.856 | +92.0 | +nan | +0.080 |
