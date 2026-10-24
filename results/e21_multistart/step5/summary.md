# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: Qf only (degenerate with WF); start: step5; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1620.1 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1095 dec; σ_eff = 0.1200 dec
- J (paper weights, our a_i) = 0.372

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.348 | 0.0701 |
| wf_gate | 4.175 | 0.704 |
| log10_nd_eff | 18.86 | 349 |
| log10_mu0 | 1.414 | 0.156 |
| log10_n_gd | 19.67 | 183 |
| e_gd | 0.1399 | 26.3 |
| s_gd | 0.09506 | 19 |
| log10_n_ta | 19.46 | 245 |
| w_ta | 0.08195 | 157 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.577 | -1.9 | +22.9 | +0.017 |
| 1.0 | 0.244 | +39.8 | +20.8 | -0.016 |
| L1_VD0.1 | 0.523 | -33.9 | +5.2 | +0.003 |
| L1_VD0.2 | 0.082 | -38.5 | -2.6 | +0.016 |
| L1_VD0.5 | 0.661 | +58.6 | -25.9 | -0.003 |
| L2_VD0.05 | 0.697 | +3.6 | -32.9 | +0.111 |
| L2_VD1 | 0.729 | +90.2 | +nan | +0.083 |
