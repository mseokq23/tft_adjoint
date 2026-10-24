# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/0.1/0.2/0.5/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 17.85}; start: results/e12_nd/log10_nd_eff_17.85/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1884.6 on 93 residuals (σ = 0.02 dec) → RMS misfit 0.0900 dec; σ_eff = 0.0942 dec
- J (paper weights, our a_i) = 0.370
- fitted points per VD: {'0.05': 24, '0.1': 17, '0.2': 15, '0.5': 9, '1': 28} (grid points at measured VG only, gaps ≤ 0.05 V)

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.342 | 0.05 |
| wf_gate | 4.071 | 0.657 |
| log10_mu0 | 1.399 | 0.0872 |
| log10_n_gd | 19.81 | 1.39 |
| e_gd | 0.06028 | 0.935 |
| s_gd | 0.1098 | 0.394 |
| log10_n_ta | 18.87 | 73.8 |
| w_ta | 0.05815 | 6.09 |
| log10_nd_eff | 17.85 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.654 | +2.6 | +18.9 | +0.011 |
| 0.1 | 0.648 | -25.8 | +2.7 | -0.004 |
| 0.2 | 0.070 | -26.3 | -5.0 | +0.010 |
| 0.5 | 0.847 | +87.9 | -37.7 | -0.008 |
| 1.0 | 0.336 | +37.6 | +13.0 | -0.019 |
| L2_VD0.05 | 0.798 | +11.6 | -34.0 | +0.102 |
| L2_VD1 | 0.829 | +95.9 | +nan | +0.074 |
