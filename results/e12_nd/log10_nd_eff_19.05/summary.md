# E11 calibration to measured 2W curves (L = 1 µm, VD = 0.05/1 V, VG ≥ -1.0 V)

- fixed: {'log10_nd_eff': 19.05}; start: results/e11_2W_v2/summary.json; fitted points: ID ≥ 1e-10 A, VG step 0.1 V
- χ² = 1606.3 on 54 residuals (σ = 0.02 dec) → RMS misfit 0.1091 dec; σ_eff = 0.1182 dec
- J (paper weights, our a_i) = 0.387

| parameter | value | 1σ (σ_eff) |
|---|---|---|
| wf_sd | 4.351 | 0.0508 |
| wf_gate | 4.117 | 0.34 |
| log10_mu0 | 1.422 | 0.0998 |
| log10_n_gd | 19.72 | 2.06 |
| e_gd | 0.1834 | 0.0509 |
| s_gd | 0.0319 | 0.219 |
| log10_n_ta | 19.99 | 0.67 |
| w_ta | 0.09474 | 0.0733 |
| log10_nd_eff | 19.05 (fixed) | – |

| curve | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| 0.05 | 0.610 | -3.7 | +22.0 | +0.017 |
| 1.0 | 0.280 | +39.7 | +19.5 | -0.018 |
| L1_VD0.1 | 0.579 | -33.6 | +4.3 | +0.002 |
| L1_VD0.2 | 0.080 | -39.4 | -4.2 | +0.016 |
| L1_VD0.5 | 0.740 | +70.4 | -30.7 | -0.004 |
| L2_VD0.05 | 0.740 | +5.1 | -34.5 | +0.114 |
| L2_VD1 | 0.771 | +90.5 | +nan | +0.085 |
