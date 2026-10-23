# Auxiliary: 1-D gate stack + Pao-Sah current

Illustration model (numpy), not the 2-D DEVSIM model. Exact property: N_s(VG, V) = N_s(VG - V, 0), so I(VG, VD) = q µ W/L ∫_{VG-VD}^{VG} N_s(u) du for every VD and L.

## Step 5: 1-D model vs 2-D DEVSIM (results/e1)

| L (µm) | VD (V) | DEVSIM Vth / SS / Ion | 1-D Vth / SS / Ion |
|---|---|---|---|
| 1 | 0.05 | -0.130 V / 134.1 / 3.46e-05 A | -0.130 V / 134.0 / 3.51e-05 A |
| 1 | 1 | -0.177 V / 122.7 / 5.55e-04 A | -0.175 V / 122.8 / 5.61e-04 A |
| 2 | 0.05 | -0.076 V / 141.8 / 1.75e-05 A | -0.077 V / 141.6 / 1.75e-05 A |
| 2 | 1 | -0.128 V / 129.4 / 2.80e-04 A | -0.127 V / 129.5 / 2.80e-04 A |

SS at VD = 1 V (1-D): no traps 62 mV/dec, tail only 99 mV/dec, tail + Gaussian 123 mV/dec

## E5 re-fitted sets (θ from DEVSIM) in the 1-D model, relative to the Step-5 truth

| log10 ND,eff | calibration RMS (dec) | ΔVth L 2 µm VD 0.05 / 1 V (mV) | ΔVth t_IWO 5 nm | ΔVth t_IWO 7 nm | ΔVth t_ox 14 nm |
|---|---|---|---|---|---|
| 15.900 | 0.0005 | +0.2 / +0.1 | +31 | +82 | +58 |
| 16.900 | 0.0005 | +0.2 / +0.1 | +31 | +81 | +56 |
| 17.400 | 0.0006 | +0.2 / +0.1 | +23 | +63 | +39 |
| 17.900 | 0.0006 | +0.2 / +0.1 | +11 | +33 | +18 |
| 18.408 | 0.0000 | +0.0 / +0.0 | +0 | +0 | +0 |
| 18.900 | 0.0065 | -1.3 / -0.4 | -61 | -191 | -72 |
| 19.150 | 0.0226 | -3.9 / -0.6 | -188 | -540 | -190 |
| 19.400 | 0.0438 | -15.5 / -12.2 | -337 | -938 | -333 |

Thickness tests at VD = 1 V. The 2-D DEVSIM predictions for the same sets are in results/e9.
