# E12: profile likelihood of log10_nd_eff on the measured 2 %-W calibration curves (joint fit of VD = 0.05/0.1/0.2/0.5/1 V, L = 1 µm)

- best fit: log10_nd_eff = 18.488, χ² = 1836.8 on 93 residuals (σ = 0.02 dec), σ_eff = 0.0935 dec (model error included)
- Δχ²_eff = (χ² − χ²_min)·σ²/σ_eff² with χ²_min = 1836.8 (lowest over all fits); 95 % threshold 3.84

| log10_nd_eff | Δχ²_eff | RMS (dec) | ΔVth VD=1 (mV) | hold-out L2_VD0.05: E_log / ΔVth | hold-out L2_VD1: E_log / ΔVth |
|---|---|---|---|---|---|
| 16.850 | 2.25 | 0.0901 | +37.6 | 0.800 / +12 | 0.830 / +96 |
| 17.850 | 2.18 | 0.0900 | +37.6 | 0.798 / +12 | 0.829 / +96 |
| 18.350 | 2.00 | 0.0899 | +38.0 | 0.808 / +10 | 0.838 / +95 |
| 18.488 | 0.00 | 0.0889 | +39.1 | 0.858 / +11 | 0.886 / +94 |
| 19.050 | 0.95 | 0.0894 | +39.2 | 0.828 / +7 | 0.856 / +92 |

Hold-out ΔVth: the model is sampled at the measured VG points and compared with the measured curve. The digitised hold-out curves have no points around 100 nA (VG gaps of 0.4 V), so these ΔVth are approximate. Comparisons between the sets use the full simulated curves (e13, e14, e18_report).
