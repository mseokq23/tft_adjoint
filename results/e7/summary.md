# E7: which extra measurement constrains the sloppy descriptors best?

- base = calibration I-V (L = 1 µm, VD = 0.05/1.0 V), σ = 0.02 dec; C-V noise = 0.5 % of C_ox·A (5.06e-13 F)
- Δlog det F is computed without Qf (exact Qf–WF degeneracy unless t_ox varies); 'breaks Qf–WF' = λmin/λmax > 1e-12
- where Qf–WF is not broken, σ wf_gate is that of the effective work function WF − qQf/C_ox (Qf merged)
- tox/tch/'L2 (both VD)' add two I-V curves (VD = 0.05 and 1.0 V); L1_VD*/L2_VD* add one curve; CV adds one C-V curve

| added measurement | type | Δlog det F | breaks Qf–WF | σ wf_gate (eV) | σ qf (cm⁻²) | σ log10_nd_eff (dec) | σ log10_n_ta (dec) | σ w_ta (eV) | σ log10_n_gd (dec) | σ e_gd (eV) | σ s_gd (eV) | σ wf_sd (eV) | σ log10_mu0 (dec) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| (calibration only) | – | 0 | no | 0.979 | ∞ | 2.95 | 2.37 | 0.192 | 1.92 | 0.245 | 0.229 | 0.282 | 0.0243 |
| tox14+tch7+CV | combination | 37.34 | yes | 0.00113 | 1.06e+10 | 0.438 | 0.14 | 0.0101 | 0.178 | 0.0406 | 0.0285 | 0.00412 | 0.00193 |
| tox14+CV | combination | 32.48 | yes | 0.00363 | 1.69e+11 | 0.503 | 0.172 | 0.0126 | 0.215 | 0.0483 | 0.0338 | 0.0125 | 0.00273 |
| tox14+tch7 | combination | 31.40 | yes | 0.00118 | 1.09e+10 | 0.586 | 0.18 | 0.0121 | 0.234 | 0.0512 | 0.0377 | 0.0436 | 0.00532 |
| tch7 | new: t_ch | 27.34 | no | 0.00102 | ∞ | 0.863 | 0.266 | 0.017 | 0.346 | 0.0758 | 0.0558 | 0.0513 | 0.00686 |
| tch5 | new: t_ch | 22.66 | no | 0.00166 | ∞ | 1.05 | 0.299 | 0.0184 | 0.406 | 0.0841 | 0.0658 | 0.0873 | 0.00957 |
| tox14 | new: t_ox | 22.37 | yes | 0.00537 | 2.48e+11 | 1.1 | 0.282 | 0.0188 | 0.406 | 0.079 | 0.0669 | 0.133 | 0.0116 |
| CV | new: C-V | 20.97 | no | 0.022 | ∞ | 0.766 | 0.284 | 0.0205 | 0.343 | 0.0777 | 0.053 | 0.0191 | 0.00374 |
| tox10 | new: t_ox | 18.54 | yes | 0.00969 | 4.53e+11 | 1.58 | 0.371 | 0.0298 | 0.569 | 0.1 | 0.0937 | 0.136 | 0.0123 |
| all paper hold-outs | combination | 12.36 | no | 0.499 | ∞ | 1.37 | 1.31 | 0.101 | 1.11 | 0.147 | 0.135 | 0.0845 | 0.00627 |
| L2 (both VD) | paper hold-out | 7.37 | no | 0.633 | ∞ | 1.78 | 1.64 | 0.128 | 1.39 | 0.183 | 0.169 | 0.106 | 0.00736 |
| L2_VD1.0 | paper hold-out | 4.85 | no | 0.74 | ∞ | 2.14 | 1.9 | 0.15 | 1.59 | 0.211 | 0.194 | 0.127 | 0.00965 |
| L2_VD0.05 | paper hold-out | 4.70 | no | 0.7 | ∞ | 2.01 | 1.82 | 0.141 | 1.55 | 0.204 | 0.19 | 0.124 | 0.00906 |
| L1_VD0.1 | paper hold-out | 3.62 | no | 0.833 | ∞ | 2.4 | 2.04 | 0.163 | 1.65 | 0.21 | 0.197 | 0.234 | 0.02 |
| L1_VD0.2 | paper hold-out | 3.49 | no | 0.849 | ∞ | 2.46 | 2.07 | 0.168 | 1.67 | 0.215 | 0.199 | 0.24 | 0.0205 |
| L1_VD0.5 | paper hold-out | 3.25 | no | 0.891 | ∞ | 2.55 | 2.16 | 0.175 | 1.73 | 0.22 | 0.204 | 0.262 | 0.0224 |

C-V information gain vs assumed C-V noise: 0.2 % → Δlogdet 29.31, 0.5 % → Δlogdet 20.97, 1 % → Δlogdet 16.17, 2 % → Δlogdet 12.08

Self-checks (gate-charge functional and adjoint gradients):
- VG = -0.50 V: |Q_functional − Q_DEVSIM|/|Q| = 6.7e-16; gate_bias: adjoint vs re-solve FD rel 2.8e-07, wf_gate: adjoint vs re-solve FD rel 2.8e-07, log10_nd_eff: adjoint vs re-solve FD rel 6.7e-07, log10_n_ta: adjoint vs re-solve FD rel 2.7e-07
- VG = +0.50 V: |Q_functional − Q_DEVSIM|/|Q| = 3.5e-16; gate_bias: adjoint vs re-solve FD rel 8.8e-08, wf_gate: adjoint vs re-solve FD rel 8.8e-08, log10_nd_eff: adjoint vs re-solve FD rel 6.2e-07, log10_n_ta: adjoint vs re-solve FD rel 1.8e-06
