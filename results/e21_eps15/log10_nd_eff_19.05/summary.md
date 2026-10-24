# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 19.05}; start: results/e12_nd/log10_nd_eff_19.05/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 15 (changed)
- χ² = 1608.4 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1092 dec; σ_eff = 0.1183 dec
- J (paper weights, our a_i) = 0.389

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.334 | 0.0714 |
| wf_gate | 4.068 | 0.36 |
| log10_mu0 | 1.401 | 0.112 |
| log10_n_gd | 19.78 | 1.93 |
| e_gd | 0.1822 | 0.0399 |
| s_gd | 0.02835 | 0.188 |
| log10_n_ta | 20.11 | 0.353 |
| w_ta | 0.08872 | 0.113 |
| log10_nd_eff | 19.05 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.613 | -4.3 | +22.1 | +0.021 |
| 1.0 | 0.283 | +39.6 | +19.6 | -0.009 |
| L1_VD0.1 | 0.584 | -35.1 | +4.0 | +0.007 |
| L1_VD0.2 | 0.082 | -39.8 | -4.2 | +0.020 |
| L1_VD0.5 | 0.746 | +72.1 | -30.9 | +0.003 |
| L2_VD0.05 | 0.743 | +5.9 | -34.5 | +0.108 |
| L2_VD1 | 0.775 | +91.2 | +nan | +0.081 |
