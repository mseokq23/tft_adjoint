# E11 calibration to measured 4W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -0.7 V)

- fixed: Qf only (degenerate with WF); start: results/start_4W.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- χ² = 12076.3 on 44 residuals (σ = 0.02 dec) → RMS misfit 0.3313 dec; σ_eff = 0.3715 dec
- J (paper weights, our a_i) = 1.434

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.29 | 1.05 |
| wf_gate | 4.517 | 2.53 |
| log10_nd_eff | 16.67 | 0.877 |
| log10_mu0 | 1.173 | 0.669 |
| log10_n_gd | 20.56 | 1.93 |
| e_gd | 0.08577 | 0.146 |
| s_gd | 0.02001 | 0.123 |
| log10_n_ta | 19.42 | 30.1 |
| w_ta | 0.1379 | 5.22 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.634 | -34.8 | +14.2 | +0.147 |
| 1.0 | 1.308 | -180.0 | -57.7 | +0.026 |
