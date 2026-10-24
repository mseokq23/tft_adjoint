# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/0.1/0.2/0.5/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 18.35}; start: results/e12_nd/log10_nd_eff_18.35/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1880.7 on 93 residuals (σ = 0.02 dec) → RMS misfit 0.0899 dec; σ_eff = 0.0941 dec
- J (paper weights, our a_i) = 0.374
- fitted points per VD: {'0.05': 24, '0.1': 17, '0.2': 15, '0.5': 9, '1': 28} (grid points at measured VG only, gaps ≤ 0.05 V)

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.344 | 0.0459 |
| wf_gate | 4.104 | 0.595 |
| log10_mu0 | 1.403 | 0.0819 |
| log10_n_gd | 19.79 | 0.702 |
| e_gd | 0.08488 | 0.489 |
| s_gd | 0.09588 | 0.258 |
| log10_n_ta | 19.04 | 44.4 |
| w_ta | 0.04828 | 2.65 |
| log10_nd_eff | 18.35 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.661 | +2.5 | +18.9 | +0.011 |
| 0.1 | 0.660 | -26.3 | +2.4 | -0.003 |
| 0.2 | 0.073 | -26.4 | -5.2 | +0.011 |
| 0.5 | 0.863 | +90.1 | -38.5 | -0.008 |
| 1.0 | 0.343 | +38.0 | +13.0 | -0.020 |
| L2_VD0.05 | 0.808 | +10.3 | -34.7 | +0.104 |
| L2_VD1 | 0.838 | +95.1 | +nan | +0.076 |
