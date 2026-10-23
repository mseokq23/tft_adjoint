# E3/E4 preview: identifiability and OED (final 2%-W model)

- assumed noise sigma = 0.02 dec in log10(ID); 126 calibration rows (L = 1 um, VD = 0.05/1.0 V, ID >= 10 pA)
- Fisher condition number = 4.47e+16

| parameter | typical change | 1-sigma (scaled) | 1-sigma (physical) |
|---|---|---|---|
| wf_sd | 0.05 | 5.64 | 0.282 |
| wf_gate | 0.05 | 19.3 | 0.963 |
| qf | 1e+11 | 2.44 | 2.44e+11 |
| log10_nd_eff | 0.05 | 58.9 | 2.95 |
| log10_mu0 | 0.05 | 0.486 | 0.0243 |
| log10_n_gd | 0.05 | 38.3 | 1.92 |
| e_gd | 0.01 | 24.5 | 0.245 |
| s_gd | 0.01 | 22.9 | 0.229 |
| log10_n_ta | 0.05 | 47.4 | 2.37 |
| w_ta | 0.01 | 19.2 | 0.192 |

## Sloppiest directions (smallest eigenvalues)

- λ = 3.17e-13: qf +0.99, wf_gate +0.13, log10_n_ta -0.00, log10_n_gd +0.00
- λ = 0.000176: log10_n_ta -0.61, log10_n_gd +0.45, log10_nd_eff +0.40, e_gd -0.26
- λ = 0.000296: log10_nd_eff +0.87, log10_n_gd -0.29, e_gd +0.25, s_gd -0.25

## |correlation| > 0.9

- wf_sd vs log10_mu0: +0.988
- wf_gate vs qf: -1.000
- wf_gate vs log10_n_ta: -0.970
- wf_gate vs w_ta: +0.970
- qf vs log10_n_ta: +0.970
- qf vs w_ta: -0.970
- log10_n_gd vs e_gd: -0.971
- log10_n_gd vs s_gd: +0.974
- log10_n_gd vs log10_n_ta: -0.969
- e_gd vs s_gd: -0.986
- log10_n_ta vs w_ta: -0.926

## 1-sigma (scaled) of sub-models

| parameter | all | qf fixed | qf+DOS fixed | DOS only |
|---|---|---|---|---|
| wf_sd | 5.64 | 5.64 | 2.8 | - |
| wf_gate | 19.3 | 19.6 | 0.175 | - |
| qf | 2.44 | - | - | - |
| log10_nd_eff | 58.9 | 58.9 | 0.882 | - |
| log10_mu0 | 0.486 | 0.486 | 0.245 | - |
| log10_n_gd | 38.3 | 38.3 | - | 1.3 |
| e_gd | 24.5 | 24.5 | - | 3.95 |
| s_gd | 22.9 | 22.9 | - | 1.11 |
| log10_n_ta | 47.4 | 47.4 | - | 2.25 |
| w_ta | 19.2 | 19.2 | - | 3.26 |

## OED preview: adding one hold-out condition

| condition | Δ log det F | λmin gain | cond |
|---|---|---|---|
| L2_VD1.0 | 4.92 | ×0.36 | 1.93e+17 |
| L2_VD0.05 | 4.79 | ×-2.42 | 1.97e+304 |
| L1_VD0.1 | 3.62 | ×4.61 | 1.45e+16 |
| L1_VD0.5 | 3.31 | ×-5.21 | 2.17e+304 |
| L1_VD0.2 | 3.18 | ×-0.66 | 2.15e+304 |

Best pairs: L1_VD0.1+L2_VD1.0 (7.53); L2_VD0.05+L2_VD1.0 (7.50); L1_VD0.1+L2_VD0.05 (7.36)

## Metric sensitivity per typical change (VD = 1 V)

| parameter | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|
| wf_sd | +0.1 | +0.03 | -0.0049 |
| wf_gate | +50.1 | -0.71 | -0.0113 |
| qf | -6.3 | +0.09 | +0.0014 |
| log10_nd_eff | -8.3 | +0.45 | +0.0014 |
| log10_mu0 | -7.8 | -1.05 | +0.0500 |
| log10_n_gd | -5.5 | +3.01 | +0.0001 |
| e_gd | +3.1 | +1.26 | -0.0001 |
| s_gd | -3.6 | +2.18 | +0.0001 |
| log10_n_ta | +16.0 | +3.51 | -0.0053 |
| w_ta | +17.5 | +2.39 | -0.0037 |
