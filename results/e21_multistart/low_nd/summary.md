# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: Qf only (degenerate with WF); start: results/e12_nd/log10_nd_eff_16.85/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1628.1 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1098 dec; σ_eff = 0.1203 dec
- J (paper weights, our a_i) = 0.349

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.35 | 0.0648 |
| wf_gate | 4.096 | 0.833 |
| log10_nd_eff | 18.29 | 65.3 |
| log10_mu0 | 1.418 | 0.147 |
| log10_n_gd | 19.76 | 26.7 |
| e_gd | 0.00873 | 21.6 |
| s_gd | 0.1971 | 9.31 |
| log10_n_ta | 19.33 | 32.1 |
| w_ta | 0.03025 | 11.5 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.501 | -1.8 | +23.3 | +0.016 |
| 1.0 | 0.201 | +37.9 | +21.7 | -0.019 |
| L1_VD0.1 | 0.429 | -31.4 | +7.6 | +0.001 |
| L1_VD0.2 | 0.089 | -37.8 | -0.6 | +0.015 |
| L1_VD0.5 | 0.533 | +44.0 | -20.0 | -0.005 |
| L2_VD0.05 | 0.599 | +7.4 | -27.8 | +0.112 |
| L2_VD1 | 0.633 | +92.5 | +nan | +0.084 |
