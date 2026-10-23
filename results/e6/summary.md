# E6 DOS detectability from the calibration transfer curves (Step 5 model)

- 126 log-current residuals (L = 1 um, VD = 0.05/1.0 V, ID >= 10 pA), σ = 0.02 dec
- nuisance parameters profiled: wf_sd, wf_gate, log10_nd_eff, log10_mu0 (Qf is degenerate with the gate WF)
- singular values of the projected sensitivity (unit 1e+18 cm⁻³eV⁻¹ per 0.025 eV bin): 1.38, 0.437, 0.174, 0.0787, 0.0382, 0.0209, 0.00986, 0.00429 … → 1 components above 1

1σ detection limit for a 0.1 eV-wide DOS change (cm⁻³eV⁻¹) and the Step-5 DOS in the same window:

| Ec−E (eV) | acceptor: model DOS | σ (nuisance known) | σ (profiled) | σ (ND known) | donor: model DOS | σ (nuisance known) | σ (profiled) | σ (ND known) |
|---|---|---|---|---|---|---|---|---|
| 0.05 | 5.5e+19 | 7.9e+17 | 9.1e+18 | 1.5e+18 | 1.5e+19 | 1.7e+17 | 9.1e+18 | 9.4e+17 |
| 0.15 | 2.6e+19 | 3.4e+17 | 9.7e+17 | 6.1e+17 | 1.9e+19 | 2e+17 | 9.7e+17 | 5.1e+17 |
| 0.25 | 1.2e+19 | 2.2e+17 | 6.3e+17 | 6.1e+17 | 6.4e+18 | 3.4e+17 | 6.3e+17 | 5.6e+17 |
| 0.35 | 5.7e+18 | 1.7e+17 | 2.1e+18 | 1.9e+18 | 5.3e+17 | 1.6e+18 | 2.1e+18 | 2e+18 |
| 0.45 | 2.7e+18 | 1.6e+17 | 3.3e+19 | 2.6e+18 | 1e+16 | 2.6e+19 | 3.3e+19 | 3.2e+19 |
| 0.55 | 1.3e+18 | 1.6e+17 | 8.5e+20 | 2.6e+18 | 4.5e+13 | 7e+20 | 8.5e+20 | 8.3e+20 |
| 0.65 | 6e+17 | 1.6e+17 | 1.6e+22 | 2.6e+18 | 4.3e+10 | 1.3e+22 | 1.6e+22 | 1.5e+22 |
| 0.75 | 2.8e+17 | 1.6e+17 | 2.1e+23 | 2.6e+18 | 9.2e+06 | 1.8e+23 | 2.1e+23 | 2e+23 |
| 0.85 | 1.3e+17 | 1.6e+17 | 1.1e+25 | 2.6e+18 | 0 | 9.4e+24 | 1.1e+25 | 1e+25 |
| 0.95 | 6.3e+16 | 1.6e+17 | 5e+26 | 2.6e+18 | 0 | 4.5e+26 | 5e+26 | 5e+26 |
| 1.05 | 3e+16 | 1.6e+17 | 1.4e+28 | 2.6e+18 | 0 | 2.1e+28 | 2.4e+28 | 2.4e+28 |
| 1.15 | 1.4e+16 | 1.6e+17 | 2.1e+28 | 2.6e+18 | 0 | 1e+30 | 1.2e+30 | 1.1e+30 |
| 1.25 | 6.6e+15 | 1.6e+17 | 2.2e+28 | 2.6e+18 | 0 | 4.9e+31 | 5.5e+31 | 5.5e+31 |
| 1.35 | 3.1e+15 | 1.6e+17 | 2.2e+28 | 2.6e+18 | 0 | 2.4e+33 | 2.7e+33 | 2.7e+33 |
| 1.45 | 1.5e+15 | 1.6e+17 | 2.2e+28 | 2.6e+18 | 0 | 1.2e+35 | 1.3e+35 | 1.3e+35 |
| 1.55 | 6.9e+14 | 1.6e+17 | 2.2e+28 | 2.6e+18 | 0 | 5.9e+36 | 6.6e+36 | 6.5e+36 |

Fermi level at mid-channel (VD = 0.05 V): Ec−EF front/back interface

- VG = -1.00 V: +0.708 / +0.635 eV (ID = 1.55e-17 A)
- VG = -0.50 V: +0.263 / +0.218 eV (ID = 2.11e-10 A)
- VG = +0.00 V: -0.002 / +0.055 eV (ID = 4.93e-07 A)
- VG = +0.50 V: -0.095 / +0.026 eV (ID = 6.01e-06 A)
- VG = +1.00 V: -0.125 / +0.022 eV (ID = 1.30e-05 A)
- VG = +2.50 V: -0.160 / +0.019 eV (ID = 3.46e-05 A)
