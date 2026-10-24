# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/0.1/0.2/0.5/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 16.85}; start: results/e12_nd/log10_nd_eff_16.85/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1885.9 on 93 residuals (σ = 0.02 dec) → RMS misfit 0.0901 dec; σ_eff = 0.0942 dec
- J (paper weights, our a_i) = 0.373
- fitted points per VD: {'0.05': 24, '0.1': 17, '0.2': 15, '0.5': 9, '1': 28} (grid points at measured VG only, gaps ≤ 0.05 V)

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.343 | 0.0494 |
| wf_gate | 4.054 | 0.718 |
| log10_mu0 | 1.401 | 0.0883 |
| log10_n_gd | 19.82 | 1.16 |
| e_gd | 0.05934 | 0.633 |
| s_gd | 0.1096 | 0.311 |
| log10_n_ta | 19.22 | 35.5 |
| w_ta | 0.03525 | 1.12 |
| log10_nd_eff | 16.85 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.655 | +3.0 | +19.1 | +0.010 |
| 0.1 | 0.650 | -25.7 | +2.7 | -0.004 |
| 0.2 | 0.070 | -26.5 | -5.1 | +0.010 |
| 0.5 | 0.849 | +88.1 | -37.8 | -0.009 |
| 1.0 | 0.337 | +37.6 | +13.0 | -0.020 |
| L2_VD0.05 | 0.800 | +11.9 | -33.9 | +0.102 |
| L2_VD1 | 0.830 | +96.2 | +nan | +0.075 |
