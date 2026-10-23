# E1 forward reproduction

| run | VD (V) | Vth (V) | SS (mV/dec) | Ion (A) |
|---|---|---|---|---|
| step0 | 0.05 | +nan | nan | 7.708e-10 |
| step0 | 1.0 | +nan | nan | 2.066e-09 |
| step1 | 0.05 | +0.206 | 68.0 | 1.308e-05 |
| step1 | 1.0 | +0.193 | 65.2 | 2.044e-04 |
| step2 | 0.05 | -0.196 | 67.9 | 1.551e-05 |
| step2 | 1.0 | -0.208 | 65.3 | 2.529e-04 |
| step3 | 0.05 | -0.241 | 68.5 | 1.573e-05 |
| step3 | 1.0 | -0.255 | 65.3 | 2.573e-04 |
| step4 | 0.05 | -0.062 | 143.9 | 1.441e-05 |
| step4 | 1.0 | -0.116 | 131.2 | 2.312e-04 |
| step5 | 0.05 | -0.130 | 134.1 | 3.458e-05 |
| step5 | 1.0 | -0.177 | 122.7 | 5.549e-04 |
| step5_L2 | 0.05 | -0.076 | 141.8 | 1.747e-05 |
| step5_L2 | 1.0 | -0.128 | 129.4 | 2.799e-04 |
| step5_dosCoarse | 1.0 | -0.176 | 122.8 | 5.546e-04 |
| step5_dosFine | 1.0 | -0.177 | 122.7 | 5.550e-04 |
| step5_meshFine | 1.0 | -0.177 | 122.7 | 5.514e-04 |

## Expected-direction checks


## Numerical convergence (vs default step5, VD = 1 V)

- step5_dosCoarse: dVth +0.67 mV, dSS +0.08 mV/dec, dlogIon -0.0003 dec
- step5_dosFine: dVth -0.17 mV, dSS -0.02 mV/dec, dlogIon +0.0001 dec
- step5_meshFine: dVth +0.07 mV, dSS +0.03 mV/dec, dlogIon -0.0028 dec
