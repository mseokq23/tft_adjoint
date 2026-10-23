# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 16.85}; start: results/e12_nd/log10_nd_eff_17.35/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- χ² = 1634.0 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1100 dec; σ_eff = 0.1192 dec
- J (paper weights, our a_i) = 0.348

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.349 | 0.0596 |
| wf_gate | 4.044 | 0.729 |
| log10_mu0 | 1.413 | 0.124 |
| log10_n_gd | 19.79 | 2.91 |
| e_gd | 1.75e-06 | 2.78 |
| s_gd | 0.1872 | 1.42 |
| log10_n_ta | 19.49 | 22 |
| w_ta | 0.03 | 0.553 |
| log10_nd_eff | 16.85 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.520 | +0.4 | +23.7 | +0.013 |
| 1.0 | 0.212 | +38.5 | +21.2 | -0.021 |
| L1_VD0.1 | 0.452 | -31.2 | +7.1 | -0.002 |
| L1_VD0.2 | 0.085 | -38.3 | -1.5 | +0.012 |
| L1_VD0.5 | 0.566 | +47.8 | -21.9 | -0.008 |
| L2_VD0.05 | 0.624 | +9.0 | -28.1 | +0.109 |
| L2_VD1 | 0.658 | +93.7 | +nan | +0.082 |
