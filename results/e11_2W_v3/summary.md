# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: Qf only (degenerate with WF); start: results/e12_nd/log10_nd_eff_19.05/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- material: Nc = 5e+18 cm⁻³, ε_ch = 10 (defaults)
- χ² = 1589.4 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1085 dec; σ_eff = 0.1189 dec
- J (paper weights, our a_i) = 0.412

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.353 | 0.063 |
| wf_gate | 3.994 | 1.01 |
| log10_nd_eff | 18.71 | 3.96 |
| log10_mu0 | 1.432 | 0.157 |
| log10_n_gd | 20.06 | 5.08 |
| e_gd | 0.1774 | 0.107 |
| s_gd | 0.02 | 0.294 |
| log10_n_ta | 20.45 | 3.65 |
| w_ta | 0.0357 | 0.254 |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.655 | -7.6 | +20.5 | +0.017 |
| 1.0 | 0.326 | +39.2 | +18.2 | -0.019 |
| L1_VD0.1 | 0.650 | -34.4 | +3.4 | +0.003 |
| L1_VD0.2 | 0.089 | -42.9 | -6.6 | +0.016 |
| L1_VD0.5 | 0.837 | +84.5 | -36.3 | -0.004 |
| L2_VD0.05 | 0.798 | +8.0 | -35.8 | +0.118 |
| L2_VD1 | 0.827 | +90.9 | +nan | +0.089 |
