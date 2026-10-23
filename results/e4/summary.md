# E4 add-on: can a second gate-oxide thickness separate Qf from the gate work function?

- calibration data only (L = 1 um, t_ox = 7 nm): smallest/largest Fisher eigenvalue = 2.2e-17 -> structurally non-identifiable; null direction (scaled units): wf_gate +0.13, qf +0.99

| added data (VD = 0.05 + 1.0 V each) | λmin/λmax | σ(qf) [1e11 cm⁻²] | σ(wf_gate) [eV] | corr(wf_gate, qf) | Δlog det F (qf excluded) |
|---|---|---|---|---|---|
| tox10 | 2.6e-08 | 4.53 | 0.00969 | -0.964 | 18.54 |
| tox14 | 5.2e-08 | 2.48 | 0.00537 | -0.960 | 22.37 |
| L2 | -6.2e-17 | ∞ | ∞ | – | 7.37 |
| tox10+tox14 | 4.8e-08 | 2.45 | 0.00531 | -0.963 | 25.14 |
| tox10+L2 | 2.8e-08 | 3.91 | 0.00846 | -0.966 | 22.74 |
| all | 4.9e-08 | 2.07 | 0.00455 | -0.966 | 28.12 |

1-sigma of every parameter (scaled units, typical change = 1) when the data identify all 10:

| parameter | tox10 | tox14 | tox10+tox14 | tox10+L2 | all |
|---|---|---|---|---|---|
| wf_sd | 2.72 | 2.66 | 2.33 | 1.33 | 1.28 |
| wf_gate | 0.194 | 0.107 | 0.106 | 0.169 | 0.091 |
| qf | 4.53 | 2.48 | 2.45 | 3.91 | 2.07 |
| log10_nd_eff | 31.6 | 22 | 18.8 | 25.2 | 16.4 |
| log10_mu0 | 0.245 | 0.233 | 0.201 | 0.104 | 0.0993 |
| log10_n_gd | 11.4 | 8.12 | 6.91 | 8.89 | 5.87 |
| e_gd | 10 | 7.9 | 6.62 | 7.7 | 5.54 |
| s_gd | 9.37 | 6.69 | 5.68 | 7.51 | 4.97 |
| log10_n_ta | 7.42 | 5.64 | 4.76 | 5.67 | 3.91 |
| w_ta | 2.98 | 1.88 | 1.66 | 2.49 | 1.47 |
