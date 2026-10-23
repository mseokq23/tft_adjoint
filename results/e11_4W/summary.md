# E11 calibration to measured 4W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -0.7 V)

- fixed: Qf only (degenerate with WF); start: w4; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- χ² = 10897.5 on 45 residuals (σ = 0.02 dec) → RMS misfit 0.3112 dec; σ_eff = 0.3480 dec
- J (paper weights, our a_i) = nan

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.445 | 0.0252 |
| wf_gate | 4.683 | 0.196 |
| log10_nd_eff | 19.3 | 0.279 |
| log10_mu0 | 1.468 | 0.322 |
| log10_n_gd | 21 | 2.92 |
| e_gd | 0.5 | 1.82 |
| s_gd | 0.06356 | 1.18 |
| log10_n_ta | 18 | 1.23e+03 |
| w_ta | 0.03 | 65.6 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.413 | -60.2 | +29.0 | +0.144 |
| 1.0 | 0.463 | -176.9 | +nan | -0.293 |
