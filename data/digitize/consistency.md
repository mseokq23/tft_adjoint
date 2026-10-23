# Digitisation consistency: the same measured curve read from different panels

Differences above the 10 pA floor (log10 ID > -11); linear panels only above 2e-5 A.

| curve | panel A | panel B | common bins | RMS Δ (dec) | max |Δ| (dec) |
|---|---|---|---|---|---|
| 2W L=1 VD=1 | 5b_log | 5b_lin | 31 | 0.022 | 0.087 |
| 2W L=1 VD=1 | 5b_log | 6a_log | 11 | 0.017 | 0.030 |
| 2W L=1 VD=1 | 5b_log | 6b_log | 34 | 0.010 | 0.033 |
| 2W L=1 VD=1 | 5b_log | 6b_lin | 24 | 0.042 | 0.128 |
| 2W L=1 VD=1 | 5b_log | 7a_log | 44 | 0.017 | 0.075 |
| 2W L=1 VD=1 | 5b_log | 7a_lin | 23 | 0.040 | 0.094 |
| 2W L=1 VD=1 | 5b_log | 7c_log | 48 | 0.042 | 0.254 |
| 2W L=1 VD=1 | 5b_lin | 6b_log | 30 | 0.016 | 0.074 |
| 2W L=1 VD=1 | 5b_lin | 6b_lin | 29 | 0.032 | 0.133 |
| 2W L=1 VD=1 | 5b_lin | 7a_log | 34 | 0.013 | 0.032 |
| 2W L=1 VD=1 | 5b_lin | 7a_lin | 26 | 0.027 | 0.086 |
| 2W L=1 VD=1 | 5b_lin | 7c_log | 34 | 0.011 | 0.033 |
| 2W L=1 VD=1 | 6a_log | 6b_log | 7 | 0.014 | 0.024 |
| 2W L=1 VD=1 | 6a_log | 7a_log | 7 | 0.017 | 0.033 |
| 2W L=1 VD=1 | 6a_log | 7c_log | 9 | 0.013 | 0.030 |
| 2W L=1 VD=1 | 6b_log | 6b_lin | 23 | 0.042 | 0.134 |
| 2W L=1 VD=1 | 6b_log | 7a_log | 34 | 0.012 | 0.033 |
| 2W L=1 VD=1 | 6b_log | 7a_lin | 21 | 0.030 | 0.089 |
| 2W L=1 VD=1 | 6b_log | 7c_log | 35 | 0.018 | 0.075 |
| 2W L=1 VD=1 | 6b_lin | 7a_log | 27 | 0.037 | 0.138 |
| 2W L=1 VD=1 | 6b_lin | 7a_lin | 20 | 0.031 | 0.092 |
| 2W L=1 VD=1 | 6b_lin | 7c_log | 27 | 0.032 | 0.124 |
| 2W L=1 VD=1 | 7a_log | 7a_lin | 25 | 0.034 | 0.098 |
| 2W L=1 VD=1 | 7a_log | 7c_log | 45 | 0.011 | 0.029 |
| 2W L=1 VD=1 | 7a_lin | 7c_log | 25 | 0.028 | 0.072 |
| 2W L=1 VD=0.05 | 6a_log | 6b_log | 30 | 0.017 | 0.026 |
| 2W L=1 VD=0.05 | 6a_log | 7a_log | 16 | 0.044 | 0.167 |
| 2W L=1 VD=0.05 | 6a_log | 6a_lin | 8 | 0.052 | 0.056 |
| 2W L=1 VD=0.05 | 6a_log | 6b_lin | 6 | 0.033 | 0.037 |
| 2W L=1 VD=0.05 | 6b_log | 7a_log | 18 | 0.016 | 0.048 |
| 2W L=1 VD=0.05 | 6b_log | 6a_lin | 8 | 0.039 | 0.044 |
| 2W L=1 VD=0.05 | 6b_log | 6b_lin | 6 | 0.018 | 0.025 |
| 2W L=1 VD=0.05 | 6a_lin | 6b_lin | 6 | 0.021 | 0.028 |
| 2W L=1 VD=0.1 | 6a_log | 6a_lin | 18 | 0.039 | 0.055 |
| 2W L=1 VD=0.2 | 6a_log | 6a_lin | 24 | 0.027 | 0.055 |
| 2W L=2 VD=1 | 6b_log | 6b_lin | 18 | 0.007 | 0.014 |
| 4W L=1 VD=1 | 7c_log | 7a_log | 22 | 0.057 | 0.188 |
| 4W L=1 VD=1 | 7c_log | 7a_lin | 12 | 0.014 | 0.032 |
| 4W L=1 VD=1 | 7a_log | 7a_lin | 11 | 0.025 | 0.046 |

All pairs: RMS 0.028 dec → single-reading σ ≈ 0.020 dec

Points removed by the monotonicity check (transfer curves must increase with VG above the 10 pA floor):

- 2W L=1 VD=1: none
- 2W L=1 VD=0.05: VG +0.15, VG +0.10
- 2W L=1 VD=0.1: VG -1.00
- 2W L=1 VD=0.2: none
- 2W L=1 VD=0.5: none
- 2W L=2 VD=1: VG -0.95
- 2W L=2 VD=0.05: none
- 4W L=1 VD=1: VG -0.50
- 4W L=1 VD=0.05: VG -0.20, VG +0.75, VG +1.75, VG +1.70, VG +1.65, VG +1.60, VG +1.55
