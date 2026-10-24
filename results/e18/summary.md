# E18: does fitting the drain-bias transfer curves as well fix ND,eff?

Calibration data: the paper's pair (L = 1 µm, VD = 0.05 / 1 V) vs all five L = 1 µm drain biases fitted jointly (VD = 0.05 / 0.1 / 0.2 / 0.5 / 1 V). Criterion: below the 95 % threshold (3.84) over at least 10× in ND,eff.

## Measured (digitised) data

Joint fit: 93 residuals {'0.05': 24, '0.1': 17, '0.2': 15, '0.5': 9, '1': 28}; χ²_min = 1836.8, σ_eff = 0.0935 dec. Pair fit (E12): χ²_min = 1589.4, σ_eff = 0.1189 dec.

| log10 ND,eff | kind | joint: Δχ²_eff | joint: RMS (dec) | joint: nfev | L 2 µm predicted Vth VD 0.05 / 1 V (mV, full curve) | L 2 µm ΔVth VD 0.05 / 1 V vs measured (mV, across the gap) |
|---|---|---|---|---|---|---|
| 16.850 | fixed | 2.25 | 0.0901 | 20 | -222.9 / -293.8 | +11.9 / +96.2 |
| 17.850 | fixed | 2.18 | 0.0900 | 20 | -223.5 / -294.0 | +11.6 / +95.9 |
| 18.350 | fixed | 2.00 | 0.0899 | 20 | -224.8 / -294.0 | +10.3 / +95.1 |
| 18.488 | free fit | 0.00 | 0.0889 | 39 | -230.3 / -296.4 | +11.1 / +93.8 |
| 19.050 | fixed | 0.95 | 0.0894 | 20 | -232.3 / -295.8 | +6.7 / +92.0 |

Equivalent at 95 % (joint): log10 ND,eff 16.85–19.05 (158× range). Pair fit (E12): 16.85–19.05 (158×).

Spread of the predicted Vth between the sets (full simulated curves): joint sets L2_VD0.05 9.4, L2_VD1 2.5, cal_0.1 3.0, cal_0.2 1.3, cal_0.5 1.6 mV; pair sets (E12) L2_VD0.05 8.2, L2_VD1 3.6, L1_VD0.1 3.6, L1_VD0.2 1.2, L1_VD0.5 1.3 mV (cal_*: fitted in the joint calibration; L1_*: hold-out tests of the pair fit).

The digitised L = 2 µm curves have no points around 100 nA: VD = 0.05 V: 1.6e-09 A at -0.55 V to 2.4e-07 A at -0.15 V; VD = 1 V: 1.1e-08 A at -0.55 V to 1.3e-06 A at -0.15 V. The ΔVth against them samples the model at the measured VG and interpolates both curves across the gap; sampling alone moves the model Vth by +15 to +24 mV, so the absolute ΔVth is uncertain by about this much.

## Synthetic twin (truth = Step 5, σ = 0.02 dec, perfect model)

Joint data: 160 residuals (VD = 0.05/0.1/0.2/0.5/1 V).

| log10 ND,eff | joint: χ² | joint: RMS (dec) | joint: nfev | pair (E5): χ² | joint: L 2 µm ΔVth VD 0.05 / 1 V vs truth (mV) |
|---|---|---|---|---|---|
| 15.900 | 0.00582 | 0.0001 | 14 | 0.00319 | +0.02 / +0.03 |
| 17.900 | 0.00453 | 0.0001 | 20 | 0.00341 | +0.06 / +0.03 |
| 18.408 | 0 | 0.0000 | 0 | 0 | – (truth) |
| 18.900 | 0.644 | 0.0013 | 20 | 0.284 | -0.79 / -0.29 |

Equivalent at 95 % (joint): log10 ND,eff 15.90–18.90 (1000× range). Pair (E5): 15.90–18.90 (1000×).

