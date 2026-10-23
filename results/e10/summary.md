# E10: paper parameter sets in DEVSIM vs digitised measurements

Our J uses the paper's weights with our own normalisation (paper a_i unpublished): compare trends, not absolute J.

## Measured metrics (digitised)

| curve | Vth (V) | SS (mV/dec) | Ion (A) |
|---|---|---|---|
| 2W L=1 VD=0.05 | -0.286 | 123.4 | 2.245e-05 |
| 2W L=1 VD=0.1 | -0.242 | 146.7 | 4.693e-05 |
| 2W L=1 VD=0.2 | -0.084 | 208.0 | 8.651e-05 |
| 2W L=1 VD=0.5 | -0.138 | 204.1 | 2.121e-04 |
| 2W L=1 VD=1 | -0.385 | 111.3 | 3.700e-04 |
| 2W L=2 VD=0.05 | -0.220 | 178.7 | 1.028e-05 |
| 2W L=2 VD=1 | -0.367 | nan | 1.719e-04 |
| 4W L=1 VD=0.05 | +0.358 | 210.0 | 7.384e-06 |
| 4W L=1 VD=1 | +0.286 | 210.8 | 1.394e-04 |

## Calibration path (2 %-W, L = 1 µm)

| step | J ours | J paper | VD | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|---|---|---|
| step0 | nan | 2.589 | 0.05 | 6.408 | +nan | +nan | -4.469 |
| step0 | nan | 2.589 | 1.0 | 7.315 | +nan | +nan | -5.337 |
| step1 | 7.050 | 1.537 | 0.05 | 5.180 | +509 | -49.5 | -0.362 |
| step1 | 7.050 | 1.537 | 1.0 | 5.842 | +578 | -46.0 | -0.409 |
| step2 | 2.091 | 0.574 | 0.05 | 1.997 | +90 | -55.5 | -0.266 |
| step2 | 2.091 | 0.574 | 1.0 | 2.343 | +177 | -46.0 | -0.284 |
| step3 | 1.546 | 0.443 | 0.05 | 1.713 | +45 | -54.9 | -0.258 |
| step3 | 1.546 | 0.443 | 1.0 | 1.992 | +130 | -45.9 | -0.274 |
| step4 | 2.960 | 0.282 | 0.05 | 1.318 | +224 | +20.5 | -0.306 |
| step4 | 2.960 | 0.282 | 1.0 | 1.457 | +269 | +19.9 | -0.335 |
| step5 | 2.162 | 0.181 | 0.05 | 1.134 | +156 | +10.7 | +0.074 |
| step5 | 2.162 | 0.181 | 1.0 | 1.213 | +208 | +11.4 | +0.045 |

Paper (Sentaurus) Step 5 at VD = 1 V: |ΔVth| = 31.5 mV, |ΔIon| = 0.022 dec

## Step 5 at the paper's hold-out conditions

| condition | E_log (dec) | ΔVth (mV) | ΔSS (mV/dec) | Δlog Ion (dec) |
|---|---|---|---|---|
| L1_VD0.1 | 1.302 | +100 | -28.9 | +0.061 |
| L1_VD0.2 | 0.697 | +157 | -39.9 | +0.072 |
| L1_VD0.5 | 1.708 | +207 | -78.2 | +0.056 |
| L2_VD0.05 | 1.265 | +163 | -35.7 | +0.116 |
| L2_VD1 | 1.323 | +239 | +nan | +0.080 |

Paper: drain-bias test avg/max |ΔVth| 24.9/36.1 mV, |ΔIon| 0.013/0.022 dec; L = 2 µm 27.4/46.2 mV, 0.050/0.062 dec

## 4 %-W recalibration vs 4 %-W data (VG ≥ −0.7 V)

- VD = 0.05 V: E_log 1.214 dec, ΔVth +200 mV, ΔSS -36.5 mV/dec, Δlog Ion +0.308 dec
- VD = 1.0 V: E_log 1.893 dec, ΔVth +212 mV, ΔSS -32.0 mV/dec, Δlog Ion +0.150 dec

Paper 4 %-W: avg/max |ΔIon| 0.031/0.057 dec, |ΔVth| 28.2/47.3 mV
