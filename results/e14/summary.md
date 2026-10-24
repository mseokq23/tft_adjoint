# E14: the paper's transfer tests cannot tell calibration-equivalent parameter sets apart

## 1. Transfer tests in the paper's metrics

Drain-bias split: L = 1 µm, VD = 0.05/0.1/0.2/0.5/1.0 V. Length split: L = 2 µm, VD = 0.05/1.0 V. avg/max of |ΔVth| (mV) and |Δlog10 Ion| (dec).

| model / set | log10 ND,eff | fit to calibration curves | drain: Vth avg/max | drain: Ion avg/max | length: Vth avg/max | length: Ion avg/max |
|---|---|---|---|---|---|---|
| paper (Sentaurus, Step 5) | 18.41 | J = 0.181 | 24.9/36.1 | 0.013/0.022 | 27.4/46.2 | 0.05/0.062 |
| twin: re-fit vs truth | 15.90 | χ² = 0.00319 | 0.01/0.02 | 0.0001/0.0001 | 0.02/0.02 | 0.0001/0.0001 |
| twin: re-fit vs truth | 16.90 | χ² = 0.00307 | 0.01/0.02 | 0.0001/0.0001 | 0.02/0.02 | 0.0001/0.0001 |
| twin: re-fit vs truth | 17.40 | χ² = 0.00361 | 0.01/0.03 | 0.0000/0.0000 | 0.05/0.06 | 0.0002/0.0002 |
| twin: re-fit vs truth | 17.90 | χ² = 0.00341 | 0.02/0.03 | 0.0000/0.0000 | 0.05/0.07 | 0.0003/0.0003 |
| twin: re-fit vs truth | 18.90 | χ² = 0.284 | 0.18/0.25 | 0.0005/0.0007 | 0.52/0.76 | 0.0043/0.0043 |
| twin: re-fit vs truth | 19.15 | χ² = 4.18 | 0.64/0.94 | 0.0004/0.0006 | 1.50/2.44 | 0.0160/0.0160 |
| twin: re-fit vs truth | 19.40 | χ² = 188 | 8.85/10.94 | 0.0206/0.0215 | 11.98/13.69 | 0.0388/0.0404 |
| ours vs measured | 16.85 | Δχ²_eff = 1.26 (RMS 0.1100 dec) | 31.2/47.8 | 0.011/0.021 | 51.4/93.7 | 0.095/0.109 |
| ours vs measured | 17.35 | Δχ²_eff = 1.25 (RMS 0.1100 dec) | 31.5/48.8 | 0.011/0.020 | 51.2/93.6 | 0.095/0.109 |
| ours vs measured | 17.85 | Δχ²_eff = 1.19 (RMS 0.1099 dec) | 31.4/48.6 | 0.011/0.021 | 50.6/93.1 | 0.096/0.110 |
| ours vs measured | 18.35 | Δχ²_eff = 1.07 (RMS 0.1098 dec) | 32.1/50.8 | 0.011/0.020 | 49.6/92.3 | 0.097/0.111 |
| ours vs measured | 18.71 | Δχ²_eff = 0.00 (RMS 0.1085 dec) | 41.7/84.5 | 0.012/0.019 | 49.5/90.9 | 0.103/0.118 |
| ours vs measured | 18.85 | Δχ²_eff = 0.88 (RMS 0.1096 dec) | 34.4/58.4 | 0.011/0.017 | 46.8/90.1 | 0.098/0.112 |
| ours vs measured | 19.05 | Δχ²_eff = 0.48 (RMS 0.1091 dec) | 37.4/70.4 | 0.012/0.018 | 47.8/90.5 | 0.100/0.114 |

Twin: noise level at which the calibration curves alone would reject each set (χ² = 3.84, perfect model): 10^15.90: 5.8e-04 dec, 10^16.90: 5.7e-04 dec, 10^17.40: 6.1e-04 dec, 10^17.90: 6.0e-04 dec, 10^18.90: 5.4e-03 dec.

Measured sets: spread of the Ion aggregates across ND,eff = 10^16.85–10^19.05: drain_ion_avg 0.001, drain_ion_max 0.0037, length_ion_avg 0.0081, length_ion_max 0.0089 (dec). For Vth, the sets are compared on their full simulated curves (section 2); the Vth aggregates above depend on the digitisation gaps (their spread is kept in summary.json as measured_aggregate_spread_vs_points).

Caution: all five digitised transfer curves have no points around 100 nA (VG gaps: L1_VD0.1 0.55 V, L1_VD0.2 1.10 V, L1_VD0.5 1.20 V, L2_VD0.05 0.40 V, L2_VD1 0.40 V), so the errors against the measured curves (model sampled at the measured VG points, Vth interpolated across the gap) are approximate. For L = 1 µm at VD 0.1/0.2/0.5 V the measured Vth is essentially undetermined, so the drain-split Vth columns are not meaningful; for L = 2 µm, sampling alone moves the model Vth by +16 to +24 mV.

The paper's objective J (paper weights, our normalisation a_i) over the same measured sets: 0.348–0.412.

## 2. Spread across the equivalent sets (max − min)

Digitisation noise expressed in Vth: σ_logI × SS = 0.02 dec × 123/111 mV/dec = 2.5/2.2 mV (VD = 0.05/1 V). The paper's own transfer errors: 25–46 mV. Measured E_log: RMS log10 error against the measured points above 1 nA (range over the sets, and its spread).

| test | in the paper's split | twin: ΔVth spread over ND,eff 10^15.9–10^18.9 (mV) | measured-fit sets: spread of the predicted Vth over ND,eff 10^16.85–10^19.05 (mV) | measured: E_log above 1 nA (dec) |
|---|---|---|---|---|
| L1 VD0.1 | yes | 0.14 | 3.6 | 0.038–0.039 (spread 0.001) |
| L1 VD0.2 | yes | 0.27 | 1.2 | 0.021–0.023 (spread 0.002) |
| L1 VD0.5 | yes | 0.20 | 1.3 | 0.017–0.020 (spread 0.003) |
| L2 VD0.05 | yes | 0.83 | 8.2 | 0.112–0.118 (spread 0.006) |
| L2 VD1 | yes | 0.32 | 3.6 | 0.157–0.161 (spread 0.004) |
| calibration VD 0.05 V | yes | 0.24 | 8.2 | (fitted) |
| calibration VD 1 V | yes | 0.14 | 1.3 | (fitted) |
| t_IWO 7 nm, VD 1 V | no (recommended) | 273 | 376 | predictions differ by up to 0.81 dec RMS |
| t_ox 14 nm, VD 1 V | no (recommended) | 129 | 133 | predictions differ by up to 0.36 dec RMS |
| quasi-static C–V (max ΔC, % of C_ox·A) | no (recommended) | 6.9 % | 8.2 % | – |

Twin: ΔVth of each re-fitted set relative to the truth (truth included as 0). Measured-fit sets: spread of the Vth that the sets predict on their full simulated curves (the measured transfer curves have digitisation gaps around 100 nA; the spread of the errors against the measured points is kept in summary.json as measured_vs_points_mV). For tests outside the paper's split no data exist: the spread is between the sets' predictions (E9 vs the truth, E13 between the measured-fit sets).

## 3. Linear check at the Step-5 truth: predictive 1σ of each test given the calibration data

σ_pred/σ < 1: in the linear approximation the calibration-consistent parameter sets predict the same outcome within the test's noise, so passing the test says little about the parameters (a set at the edge of the 95 % region can deviate by about 2·σ_pred). Qf fixed (exactly degenerate with the gate WF). Calibration noise 0.02 dec.

| test | noise σ | σ_pred/σ (RMS over the curve) | σ_pred/σ (max) | Δlog det F |
|---|---|---|---|---|
| VD 0.1 V (paper) | 0.02 dec | 0.27 | 0.70 | 3.6 |
| VD 0.2 V (paper) | 0.02 dec | 0.27 | 0.88 | 3.5 |
| VD 0.5 V (paper) | 0.02 dec | 0.25 | 0.90 | 3.2 |
| L 2 µm, VD 0.05 V (paper) | 0.02 dec | 0.43 | 0.64 | 4.7 |
| L 2 µm, VD 1 V (paper) | 0.02 dec | 0.43 | 0.91 | 4.9 |
| t_ox 10 nm | 0.02 dec | 55.08 | 172.61 | 18.5 |
| t_ox 14 nm | 0.02 dec | 118.29 | 394.35 | 22.4 |
| t_IWO 5 nm | 0.02 dec | 101.97 | 391.42 | 22.7 |
| t_IWO 7 nm | 0.02 dec | 239.33 | 919.95 | 27.3 |
| quasi-static C–V | 0.5 % C_ox·A | 32.30 | 93.56 | 21.0 |

## 4. Robustness of the measured-data equivalence

Reference χ²_min = 1589.4 (lowest default-material fit), σ_eff = 0.1189 dec. Fits are limited to 25–40 evaluations, so each χ² is an upper bound of the minimum under its constraints.

| check | fixed | free | χ² | Δχ²_eff | equivalent at 95 %? |
|---|---|---|---|---|---|
| tail width at the paper's 4%-W value (2%-W data) | w_ta = 0.238 | wf_sd, wf_gate, log10_nd_eff, log10_mu0, log10_n_gd, e_gd, s_gd, log10_n_ta | 1621.9 | 0.92 | yes |
| Nc = 2e18 cm⁻³ instead of 5e18 (vs the better of the two ends) | log10_nd_eff = 16.85 | wf_sd, wf_gate, log10_mu0, log10_n_gd, e_gd, s_gd, log10_n_ta, w_ta | 1665.1 | 1.91 | yes |
| Nc = 2e18 cm⁻³ instead of 5e18 (vs the better of the two ends) | log10_nd_eff = 19.05 | wf_sd, wf_gate, log10_mu0, log10_n_gd, e_gd, s_gd, log10_n_ta, w_ta | 1598.6 | 0.00 | yes |
| ε_ch = 15 instead of 10 (vs the better of the two ends) | log10_nd_eff = 16.85 | wf_sd, wf_gate, log10_mu0, log10_n_gd, e_gd, s_gd, log10_n_ta, w_ta | 1631.4 | 0.66 | yes |
| ε_ch = 15 instead of 10 (vs the better of the two ends) | log10_nd_eff = 19.05 | wf_sd, wf_gate, log10_mu0, log10_n_gd, e_gd, s_gd, log10_n_ta, w_ta | 1608.4 | 0.00 | yes |
| paper's 4%-W DOS and ND,eff, only WF_sd/WF_gate/µ0 fitted | log10_nd_eff = 16.8055, log10_n_gd = 19.6435, e_gd = 0.08, s_gd = 0.224, log10_n_ta = 19.8921, w_ta = 0.238 | wf_sd, wf_gate, log10_mu0 | 2660.9 | 30.34 | no |
| paper's 2%-W (Step 5) DOS and ND,eff, only WF_sd/WF_gate/µ0 fitted | log10_nd_eff = 18.4082, log10_n_gd = 19.3404, e_gd = 0.12, s_gd = 0.08, log10_n_ta = 19.8921, w_ta = 0.133 | wf_sd, wf_gate, log10_mu0 | 2065.6 | 13.48 | no |
