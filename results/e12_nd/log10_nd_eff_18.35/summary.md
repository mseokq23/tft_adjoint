# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 18.35}; start: results/e11_2W_v2/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- χ² = 1627.3 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1098 dec; σ_eff = 0.1190 dec
- J (paper weights, our a_i) = 0.353

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.35 | 0.0539 |
| wf_gate | 4.101 | 0.563 |
| log10_mu0 | 1.415 | 0.108 |
| log10_n_gd | 19.74 | 1.56 |
| e_gd | 0.05945 | 2.07 |
| s_gd | 0.1531 | 1.21 |
| log10_n_ta | 18.86 | 75.1 |
| w_ta | 0.05911 | 6.79 |
| log10_nd_eff | 18.35 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.538 | -0.2 | +23.4 | +0.014 |
| 1.0 | 0.221 | +38.9 | +21.1 | -0.020 |
| L1_VD0.1 | 0.473 | -32.0 | +6.5 | -0.001 |
| L1_VD0.2 | 0.084 | -38.3 | -1.8 | +0.013 |
| L1_VD0.5 | 0.594 | +50.8 | -23.1 | -0.007 |
| L2_VD0.05 | 0.646 | +6.9 | -29.8 | +0.111 |
| L2_VD1 | 0.681 | +92.3 | +nan | +0.083 |
