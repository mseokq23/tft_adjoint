# E3b: identifiability on the subspace without Qf (Qf ≡ gate-WF shift)

- reduced Fisher (9 parameters): condition number 7.96e+07; eigenvalues 0.000176, 0.000296, 0.0416, 0.127, 0.699, 12.7, 239, 621, 1.4e+04

## Sloppy directions (scaled units; typical change = 1)

- λ = 0.000176 (1σ along it = 75 typical units): log10_n_ta -0.61, log10_n_gd +0.45, log10_nd_eff +0.40, wf_gate +0.26, e_gd -0.26
- λ = 0.000296 (1σ along it = 58 typical units): log10_nd_eff -0.87, log10_n_gd +0.29, e_gd -0.25, s_gd +0.25, log10_n_ta -0.18
- λ = 0.0416 (1σ along it = 5 typical units): w_ta +0.58, e_gd -0.53, wf_sd -0.48, s_gd -0.24, wf_gate -0.23

## 1σ (physical units)

| parameter | calibration data | + all paper hold-out conditions |
|---|---|---|
| wf_sd | 0.282 | 0.0845 |
| wf_gate | 0.979 | 0.499 |
| log10_nd_eff | 2.95 | 1.37 |
| log10_mu0 | 0.0243 | 0.00627 |
| log10_n_gd | 1.92 | 1.11 |
| e_gd | 0.245 | 0.147 |
| s_gd | 0.229 | 0.135 |
| log10_n_ta | 2.37 | 1.31 |
| w_ta | 0.192 | 0.101 |

## OED ranking of the paper's hold-out conditions (Qf removed)

| condition | Δ log det F | λmin gain | σ ratio wf_sd | σ ratio wf_gate | σ ratio log10_nd_eff | σ ratio log10_mu0 | σ ratio log10_n_gd | σ ratio e_gd | σ ratio s_gd | σ ratio log10_n_ta | σ ratio w_ta |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L2_VD1.0 | 4.85 | ×1.60 | 0.45 | 0.76 | 0.73 | 0.40 | 0.83 | 0.86 | 0.85 | 0.80 | 0.78 |
| L2_VD0.05 | 4.70 | ×1.72 | 0.44 | 0.72 | 0.68 | 0.37 | 0.81 | 0.83 | 0.83 | 0.77 | 0.73 |
| L1_VD0.1 | 3.62 | ×1.37 | 0.83 | 0.85 | 0.81 | 0.82 | 0.86 | 0.86 | 0.86 | 0.86 | 0.85 |
| L1_VD0.2 | 3.49 | ×1.31 | 0.85 | 0.87 | 0.84 | 0.84 | 0.87 | 0.88 | 0.87 | 0.88 | 0.87 |
| L1_VD0.5 | 3.25 | ×1.20 | 0.93 | 0.91 | 0.86 | 0.92 | 0.90 | 0.90 | 0.89 | 0.91 | 0.91 |

Best pairs: L1_VD0.1+L2_VD1.0 (Δlogdet 7.47, λmin ×2.12); L2_VD0.05+L2_VD1.0 (Δlogdet 7.37, λmin ×2.12); L1_VD0.2+L2_VD1.0 (Δlogdet 7.35, λmin ×2.03)
All five conditions: Δlogdet 12.36, λmin ×3.34
