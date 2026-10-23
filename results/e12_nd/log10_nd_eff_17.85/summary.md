# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 17.85}; start: results/e12_nd/log10_nd_eff_18.35/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- χ² = 1631.4 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1099 dec; σ_eff = 0.1191 dec
- J (paper weights, our a_i) = 0.349

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.35 | 0.057 |
| wf_gate | 4.072 | 0.666 |
| log10_mu0 | 1.415 | 0.118 |
| log10_n_gd | 19.81 | 2.67 |
| e_gd | 4.515e-06 | 2.53 |
| s_gd | 0.1844 | 1.29 |
| log10_n_ta | 18.81 | 104 |
| w_ta | 0.03 | 2.84 |
| log10_nd_eff | 17.85 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.524 | +0.3 | +23.6 | +0.013 |
| 1.0 | 0.214 | +38.7 | +21.1 | -0.021 |
| L1_VD0.1 | 0.457 | -31.5 | +6.9 | -0.001 |
| L1_VD0.2 | 0.085 | -38.3 | -1.6 | +0.012 |
| L1_VD0.5 | 0.573 | +48.6 | -22.2 | -0.007 |
| L2_VD0.05 | 0.629 | +8.1 | -28.7 | +0.110 |
| L2_VD1 | 0.663 | +93.1 | +nan | +0.082 |
