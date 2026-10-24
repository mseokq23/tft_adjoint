# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/0.1/0.2/0.5/1 V, VG ≥ -1.0 V)

- fixed: Qf only (degenerate with WF); start: results/e11_2W_v3/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1836.8 on 93 residuals (σ = 0.02 dec) → RMS misfit 0.0889 dec; σ_eff = 0.0935 dec
- J (paper weights, our a_i) = 0.399
- fitted points per VD: {'0.05': 24, '0.1': 17, '0.2': 15, '0.5': 9, '1': 28} (grid points at measured VG only, gaps ≤ 0.05 V)

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.347 | 0.0525 |
| wf_gate | 3.939 | 1.21 |
| log10_nd_eff | 18.49 | 7.3 |
| log10_mu0 | 1.417 | 0.133 |
| log10_n_gd | 20.04 | 2.95 |
| e_gd | 0.1627 | 0.113 |
| s_gd | 0.02 | 0.19 |
| log10_n_ta | 20.58 | 3.85 |
| w_ta | 0.03 | 0.197 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.700 | -4.9 | +16.5 | +0.012 |
| 0.1 | 0.727 | -28.7 | +1.3 | -0.002 |
| 0.2 | 0.104 | -30.7 | -7.5 | +0.012 |
| 0.5 | 0.958 | +104.7 | -44.0 | -0.007 |
| 1.0 | 0.396 | +39.1 | +12.4 | -0.020 |
| L2_VD0.05 | 0.858 | +11.1 | -36.4 | +0.108 |
| L2_VD1 | 0.886 | +93.8 | +nan | +0.080 |
