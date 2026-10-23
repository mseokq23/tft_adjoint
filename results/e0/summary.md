# E0 gradient check

- device: {'nodes': 5889, 'dos_levels': 40, 'mesh_scale': 1.0, 'solver': 'umfpack', 'qf_area_check': 1.0000000000000002, 'n_solves': 334, 't_solve_s': 67.272}
- max |adjoint-forward|/|forward| = 3.93e-03
- adjoint vs FD(re-solve, small step): median 1.68e-08, max 4.33e-03
- Taylor remainder slopes: r0 1.00 (1 expected), r1 2.00 (2 expected)

| point | ID (A) | adjoint (s) | forward (s) | FD re-solve (s) | normwise err |
|---|---|---|---|---|---|
| subthreshold | 2.270e-11 | 0.75 | 0.72 | 5.79 | 4.4e-06 |
| near-Vth | 7.147e-08 | 0.72 | 0.71 | 3.08 | 1.0e-06 |
| on-state | 4.107e-04 | 0.76 | 0.73 | 3.12 | 3.0e-08 |
| linear | 6.008e-06 | 0.73 | 0.80 | 2.88 | 1.5e-08 |
