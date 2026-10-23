# E3c: robustness of the identifiability structure (calibration conditions, σ = 0.02 dec)

| quantity | step5 (reference) | w4 | nc2e18 | eps15 |
|---|---|---|---|---|
| rows (ID ≥ 10 pA) | 126 | 108 | 124 | 126 |
| Qf–WF null λmin/λmax | 2.2e-17 | 4.5e-18 | -1.6e-17 | -2.8e-17 |
| condition number (Qf removed) | 8.0e+07 | 2.9e+14 | 1.6e+08 | 7.7e+07 |
| 1σ wf_sd | 0.282 | 0.279 | 0.106 | 0.304 |
| 1σ wf_gate | 0.979 | 0.464 | 0.502 | 0.876 |
| 1σ log10_nd_eff | 2.95 | 7.31e+03 | 5.31 | 3.25 |
| 1σ log10_mu0 | 0.0243 | 0.0409 | 0.0145 | 0.0322 |
| 1σ log10_n_gd | 1.92 | 26.3 | 1.28 | 1.45 |
| 1σ e_gd | 0.245 | 8.17 | 0.191 | 0.233 |
| 1σ s_gd | 0.229 | 5.75 | 0.212 | 0.2 |
| 1σ log10_n_ta | 2.37 | 17.8 | 0.695 | 1.49 |
| 1σ w_ta | 0.192 | 4.01 | 0.137 | 0.105 |

Sloppiest direction (Qf removed, scaled units):

- step5 (reference): log10_n_ta -0.61, log10_n_gd +0.45, log10_nd_eff +0.40, wf_gate +0.26
- w4: log10_nd_eff -1.00, e_gd -0.01, s_gd +0.00, log10_n_gd +0.00
- nc2e18: log10_nd_eff -0.94, log10_n_gd +0.20, s_gd +0.17, e_gd -0.16
- eps15: log10_nd_eff +0.85, log10_n_ta -0.32, log10_n_gd +0.24, wf_gate +0.23
