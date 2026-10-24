# E12: profile likelihood of log10_nd_eff on the measured 2 %-W calibration curves

- best fit: log10_nd_eff = 18.709, χ² = 1589.4 on 54 residuals (σ = 0.02 dec), σ_eff = 0.1189 dec (model error included)
- Δχ²_eff = (χ² − χ²_min)·σ²/σ_eff² with χ²_min = 1589.4 (lowest over all fits); 95 % threshold 3.84

| log10_nd_eff | Δχ²_eff | RMS (dec) | ΔVth VD=1 (mV) | hold-out L1_VD0.1: E_log / ΔVth | hold-out L1_VD0.2: E_log / ΔVth | hold-out L1_VD0.5: E_log / ΔVth | hold-out L2_VD0.05: E_log / ΔVth | hold-out L2_VD1: E_log / ΔVth |
|---|---|---|---|---|---|---|---|---|
| 16.850 | 1.26 | 0.1100 | +38.5 | 0.452 / -31 | 0.085 / -38 | 0.566 / +48 | 0.624 / +9 | 0.658 / +94 |
| 17.350 | 1.25 | 0.1100 | +38.7 | 0.458 / -31 | 0.084 / -38 | 0.574 / +49 | 0.629 / +9 | 0.664 / +94 |
| 17.850 | 1.19 | 0.1099 | +38.7 | 0.457 / -31 | 0.085 / -38 | 0.573 / +49 | 0.629 / +8 | 0.663 / +93 |
| 18.350 | 1.07 | 0.1098 | +38.9 | 0.473 / -32 | 0.084 / -38 | 0.594 / +51 | 0.646 / +7 | 0.681 / +92 |
| 18.709 | 0.00 | 0.1085 | +39.2 | 0.650 / -34 | 0.089 / -43 | 0.837 / +85 | 0.798 / +8 | 0.827 / +91 |
| 19.050 | 0.48 | 0.1091 | +39.7 | 0.579 / -34 | 0.080 / -39 | 0.740 / +70 | 0.740 / +5 | 0.771 / +90 |

Hold-out ΔVth: the model is sampled at the measured VG points and compared with the measured curve. The digitised hold-out curves have no points around 100 nA (VG gaps of 0.4–1.2 V), so these ΔVth are approximate (for L = 1 µm at VD 0.1/0.2/0.5 V the measured Vth is essentially undetermined). Comparisons between the sets use the full simulated curves (e13, e14, e18_report).
