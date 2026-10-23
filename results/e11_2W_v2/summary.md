# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: Qf only (degenerate with WF); start: results/e11_2W/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- χ² = 1620.5 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1096 dec; σ_eff = 0.1200 dec
- J (paper weights, our a_i) = 0.371

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.349 | 0.0678 |
| wf_gate | 4.175 | 0.591 |
| log10_nd_eff | 18.85 | 860 |
| log10_mu0 | 1.416 | 0.153 |
| log10_n_gd | 19.68 | 429 |
| e_gd | 0.1362 | 57.9 |
| s_gd | 0.09745 | 45.4 |
| log10_n_ta | 19.38 | 787 |
| w_ta | 0.08771 | 422 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.576 | -1.8 | +22.9 | +0.017 |
| 1.0 | 0.244 | +39.8 | +20.8 | -0.017 |
| L1_VD0.1 | 0.522 | -33.8 | +5.2 | +0.003 |
| L1_VD0.2 | 0.082 | -38.4 | -2.5 | +0.016 |
| L1_VD0.5 | 0.659 | +58.4 | -25.9 | -0.004 |
| L2_VD0.05 | 0.695 | +3.5 | -32.9 | +0.112 |
| L2_VD1 | 0.728 | +90.1 | +nan | +0.084 |
