# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: Qf only (degenerate with WF); start: step5; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- χ² = 33526.5 on 55 residuals (σ = 0.02 dec) → RMS misfit 0.4938 dec; σ_eff = 0.5399 dec
- J (paper weights, our a_i) = 0.368

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.351 | 0.25 |
| wf_gate | 4.068 | 2.63 |
| log10_nd_eff | 18.64 | 13.6 |
| log10_mu0 | 1.418 | 0.51 |
| log10_n_gd | 19.61 | 4.94 |
| e_gd | 0.09329 | 13.2 |
| s_gd | 0.1459 | 6.68 |
| log10_n_ta | 19.6 | 61 |
| w_ta | 0.1045 | 15.1 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.525 | +2.0 | +24.8 | +0.014 |
| 1.0 | 0.212 | +40.0 | +22.3 | -0.020 |
| L1_VD0.1 | 0.456 | -31.4 | +7.4 | -0.000 |
| L1_VD0.2 | 0.089 | -39.8 | -1.5 | +0.014 |
| L1_VD0.5 | 0.569 | +47.2 | -21.6 | -0.006 |
| L2_VD0.05 | 0.630 | +9.1 | -28.3 | +0.112 |
| L2_VD1 | 0.665 | +93.9 | +nan | +0.084 |
