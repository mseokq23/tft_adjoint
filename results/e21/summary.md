# E21: robustness of the measured-data conclusions

## 1. Channel permittivity ε_ch = 15 (default 10)

| fixed log10 ND,eff | χ² (ε_ch = 15) | χ² (ε_ch = 10, E12) | Δχ²_eff vs the better end | Δχ²_eff vs default-material χ²_min | RMS (dec) | nfev | status | parameters at a bound |
|---|---|---|---|---|---|---|---|---|
| 16.85 | 1631.4 | 1634.0 | 0.66 | 1.19 | 0.1099 | 25 | evaluation limit | w_ta (lower) |
| 19.05 | 1608.4 | 1606.3 | 0.00 | 0.54 | 0.1092 | 25 | evaluation limit | – |

σ_eff of the better end is used for its column, σ_eff of v3 for the default-material column. Both ends are equivalent (95 %) if the first column stays below 3.84. Each ε_ch = 15 fit starts from the ε_ch = 10 point with the same ND,eff (E12).

Fitted descriptors with ε_ch = 10 (E12) and 15 (these fits):

| descriptor | 16.85: ε_ch 10 | 16.85: ε_ch 15 | 19.05: ε_ch 10 | 19.05: ε_ch 15 |
|---|---|---|---|---|
| wf_sd | 4.349 | 4.334 | 4.351 | 4.334 |
| wf_gate | 4.044 | 4.033 | 4.117 | 4.068 |
| log10_mu0 | 1.413 | 1.396 | 1.422 | 1.401 |
| log10_n_gd | 19.795 | 19.808 | 19.723 | 19.783 |
| e_gd | 0.000 | 0.041 | 0.183 | 0.182 |
| s_gd | 0.187 | 0.165 | 0.032 | 0.028 |
| log10_n_ta | 19.492 | 19.682 | 19.990 | 20.112 |
| w_ta | 0.030 | 0.030 | 0.095 | 0.089 |

## 2. Free fits from different starting points (2 %-W pair, Qf at the Step-5 value)

Δχ²_eff against χ²_min = 1589.4 (lowest of all default-material fits), σ_eff = 0.1189 dec (v3).

| fit | χ² | Δχ²_eff | RMS (dec) | nfev | status | wf_sd | wf_gate | log10_nd_eff | log10_mu0 | log10_n_gd | e_gd | s_gd | log10_n_ta | w_ta | parameters at a bound |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| v2 (earlier, Step 5 path) | 1620.5 | 0.88 | 0.1096 | 30 | evaluation limit | 4.349 | 4.175 | 18.847 | 1.416 | 19.677 | 0.136 | 0.097 | 19.379 | 0.088 | – |
| v3 (earlier, from the 10^19.05 point) | 1589.4 | 0.00 | 0.1085 | 40 | evaluation limit | 4.353 | 3.994 | 18.709 | 1.432 | 20.065 | 0.177 | 0.020 | 20.452 | 0.036 | s_gd (lower) |
| step5: paper Step 5 (2%-W) | 1620.1 | 0.87 | 0.1095 | 40 | evaluation limit | 4.348 | 4.175 | 18.860 | 1.414 | 19.674 | 0.140 | 0.095 | 19.456 | 0.082 | – |
| w4: paper 4%-W values | 1625.0 | 1.01 | 0.1097 | 40 | evaluation limit | 4.347 | 4.140 | 18.526 | 1.414 | 19.738 | 0.052 | 0.167 | 18.016 | 0.092 | log10_n_ta (near lower, 0.0158 from it) |
| mid: midpoints of the fit bounds | 1624.8 | 1.00 | 0.1097 | 40 | evaluation limit | 4.348 | 4.137 | 18.530 | 1.414 | 19.732 | 0.059 | 0.163 | 18.681 | 0.045 | – |
| low_nd: E12 set at ND,eff = 10^16.85 (ND,eff free) | 1628.1 | 1.10 | 0.1098 | 40 | evaluation limit | 4.350 | 4.096 | 18.290 | 1.418 | 19.761 | 0.009 | 0.197 | 19.329 | 0.030 | w_ta (lower) |

Calibration-only Fisher 1σ stored with each fit (σ_eff of that fit, linearised at its own end point):

| fit | log10_nd_eff | log10_n_ta | log10_n_gd | wf_gate | w_ta | s_gd |
|---|---|---|---|---|---|---|
| v2 (earlier, Step 5 path) | 860 | 787 | 429 | 0.591 | 422 | 45.4 |
| v3 (earlier, from the 10^19.05 point) | 3.96 | 3.65 | 5.08 | 1.01 | 0.254 | 0.294 |
| step5: paper Step 5 (2%-W) | 349 | 245 | 183 | 0.704 | 157 | 19 |
| w4: paper 4%-W values | 2.36 | 599 | 3.06 | 0.506 | 119 | 1.98 |
| mid: midpoints of the fit bounds | 124 | 644 | 51.6 | 0.756 | 129 | 11.1 |
| low_nd: E12 set at ND,eff = 10^16.85 (ND,eff free) | 65.3 | 32.1 | 26.7 | 0.833 | 11.5 | 9.31 |

Predicted Vth (mV, 100 nA) of the calibration and transfer-test curves, from each fit's full simulated curve:

| fit | cal VD 0.05 | cal VD 1 | VD 0.1 | VD 0.2 | VD 0.5 | L 2 µm VD 0.05 | L 2 µm VD 1 |
|---|---|---|---|---|---|---|---|
| v2 (earlier, Step 5 path) | -287.5 | -345.6 | -321.6 | -339.5 | -344.7 | -234.0 | -295.5 |
| v3 (earlier, from the 10^19.05 point) | -293.3 | -346.2 | -325.2 | -340.6 | -345.3 | -235.7 | -298.2 |
| step5: paper Step 5 (2%-W) | -287.5 | -345.6 | -321.6 | -339.5 | -344.7 | -234.0 | -295.5 |
| w4: paper 4%-W values | -287.5 | -346.7 | -322.6 | -340.7 | -345.8 | -232.0 | -296.0 |
| mid: midpoints of the fit bounds | -287.7 | -346.9 | -322.8 | -340.9 | -346.0 | -232.1 | -296.2 |
| low_nd: E12 set at ND,eff = 10^16.85 (ND,eff free) | -287.4 | -347.6 | -323.2 | -341.5 | -346.7 | -230.8 | -296.5 |

6 of 6 fits are within Δχ²_eff < 3.84 of the best: v2, v3, step5, w4, mid, low_nd.

- spread of their descriptors (max − min): wf_sd 0.00595, wf_gate 0.182, log10_nd_eff 0.569, log10_mu0 0.018, log10_n_gd 0.39, e_gd 0.169, s_gd 0.177, log10_n_ta 2.44, w_ta 0.0616.
- spread of their predicted Vth (mV): cal VD 0.05 5.8, cal VD 1 2.0, VD 0.1 3.6, VD 0.2 2.0, VD 0.5 1.9, L 2 µm VD 0.05 4.9, L 2 µm VD 1 2.7.
