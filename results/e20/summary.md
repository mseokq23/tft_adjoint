# E20: what do the data identify? Fisher 1σ of every descriptor

Measured columns: linear analysis at the joint free fit to the digitised 2 %-W curves (E18, log10 ND,eff = 18.49), σ_eff = 0.0935 dec for all I–V data (model error included; for the recommended I–V tests this assumes the same model error), C–V noise as a fraction of C_ox·A. Calibration rows = the 93 residual rows of that fit (χ² 1836.8, stored 1836.8). Twin columns: Step 5, σ = 0.02 dec, rows of E3/E4/E7. Qf fixed unless stated. Recommended I–V tests: VD = 0.05/1 V, ID ≥ 10 pA as in E4/E7 (t_ox 14 nm from VG = -1.5 V, t_IWO 7 nm from -3 V, 0.1 V step; the measured calibration rows use a 100 pA floor, see the floor table below); C–V: VD = 0, VG = -1…2.5 V in 0.05 V steps. The columns with recommended tests add them to the joint (5-curve) calibration; the pair column is the VD 0.05/1 V subset of the joint rows at the same point.

| descriptor | unit | pair (VD 0.05/1 V) | joint (5 VD, E18) | joint + t_ox 14 nm | joint + t_IWO 7 nm | joint + C–V (0.5 %) | joint + C–V (2 %) | joint + all three | joint + all three (C–V 2 %) | joint + all three, Qf free | twin: pair | twin: + all three | twin: + all three, Qf free |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| wf_gate | eV | 1.32 | 1.21 | 0.0107 | 0.00698 | 0.0278 | 0.0491 | 0.00442 | 0.00513 | 0.00694 | 0.979 | 0.000858 | 0.00113 |
| wf_sd | eV | 0.0634 | 0.0525 | 0.0391 | 0.0324 | 0.0215 | 0.037 | 0.00579 | 0.0117 | 0.00594 | 0.282 | 0.00412 | 0.00412 |
| log10_mu0 | dec | 0.166 | 0.133 | 0.0747 | 0.0695 | 0.0405 | 0.0737 | 0.013 | 0.026 | 0.0132 | 0.0243 | 0.00193 | 0.00193 |
| log10_nd_eff | dec | 8 | 7.3 | 0.537 | 0.555 | 0.315 | 0.717 | 0.176 | 0.344 | 0.176 | 2.95 | 0.437 | 0.438 |
| log10_n_ta | dec | 4.22 | 3.85 | 0.503 | 0.569 | 0.234 | 0.622 | 0.174 | 0.341 | 0.174 | 2.37 | 0.14 | 0.14 |
| w_ta | eV | 0.216 | 0.197 | 0.0392 | 0.0424 | 0.0172 | 0.0459 | 0.0134 | 0.0261 | 0.0134 | 0.192 | 0.0101 | 0.0101 |
| log10_n_gd | dec | 3.2 | 2.95 | 0.817 | 0.779 | 0.229 | 0.715 | 0.195 | 0.454 | 0.196 | 1.92 | 0.178 | 0.178 |
| e_gd | eV | 0.122 | 0.113 | 0.0271 | 0.0283 | 0.0105 | 0.0302 | 0.00923 | 0.018 | 0.00923 | 0.245 | 0.0405 | 0.0406 |
| s_gd | eV | 0.206 | 0.19 | 0.0504 | 0.0488 | 0.0159 | 0.0469 | 0.013 | 0.0288 | 0.013 | 0.229 | 0.0285 | 0.0285 |
| qf | cm⁻² | – | – | – | – | – | – | – | – | 7e+10 | – | – | 1.1e+10 |

Check: the joint-calibration column reproduces the 1σ stored with the E18 fit (largest relative difference 1.3e-09; the stored values use the Jacobian of the last accepted step).

## Current floor of the recommended I–V tests

Measured point, joint calibration plus the test(s); C–V 0.5 % in the last two columns. Rows kept (both VD): tox14 (10 pA) 60, tox14 (100 pA) 59, tiwo7 (10 pA) 66, tiwo7 (100 pA) 64.

| descriptor | unit | + t_ox 14 nm, 10 pA | + t_ox 14 nm, 100 pA | + t_IWO 7 nm, 10 pA | + t_IWO 7 nm, 100 pA | + all three, 10 pA | + all three, 100 pA |
|---|---|---|---|---|---|---|---|
| wf_gate | eV | 0.0107 | 0.011 | 0.00698 | 0.00701 | 0.00442 | 0.0049 |
| wf_sd | eV | 0.0391 | 0.0391 | 0.0324 | 0.0325 | 0.00579 | 0.00604 |
| log10_mu0 | dec | 0.0747 | 0.0747 | 0.0695 | 0.07 | 0.013 | 0.0134 |
| log10_nd_eff | dec | 0.537 | 0.537 | 0.555 | 0.641 | 0.176 | 0.176 |
| log10_n_ta | dec | 0.503 | 0.503 | 0.569 | 0.618 | 0.174 | 0.174 |
| w_ta | eV | 0.0392 | 0.0392 | 0.0424 | 0.0472 | 0.0134 | 0.0134 |
| log10_n_gd | dec | 0.817 | 0.819 | 0.779 | 1.18 | 0.195 | 0.197 |
| e_gd | eV | 0.0271 | 0.0271 | 0.0283 | 0.0307 | 0.00923 | 0.00924 |
| s_gd | eV | 0.0504 | 0.0505 | 0.0488 | 0.0702 | 0.013 | 0.0131 |
