# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 16.85}; start: results/e12_nd/log10_nd_eff_16.85/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 2e+18 cm⁻³, ε_ch = 10 (changed)
- χ² = 1665.1 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1111 dec; σ_eff = 0.1203 dec
- J (paper weights, our a_i) = 0.383

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.322 | 0.0639 |
| wf_gate | 3.949 | 0.707 |
| log10_mu0 | 1.398 | 0.122 |
| log10_n_gd | 19.75 | 3.77 |
| e_gd | 4.984e-08 | 3.53 |
| s_gd | 0.1721 | 1.86 |
| log10_n_ta | 20 | 6.77 |
| w_ta | 0.03 | 0.194 |
| log10_nd_eff | 16.85 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.536 | +5.3 | +25.0 | +0.005 |
| 1.0 | 0.221 | +39.6 | +20.9 | -0.025 |
| L1_VD0.1 | 0.472 | -30.8 | +6.7 | -0.009 |
| L1_VD0.2 | 0.082 | -39.2 | -2.5 | +0.005 |
| L1_VD0.5 | 0.594 | +50.2 | -23.7 | -0.014 |
| L2_VD0.05 | 0.643 | +11.9 | -27.8 | +0.100 |
| L2_VD1 | 0.678 | +96.3 | +nan | +0.074 |
