# control_safety_slip

## Description

`control_safety_slip index switch` (manual 6.369): slip-surface safety
factor calculation in the time steps of control index `index`. The slip
surfaces are given with `safety_slip_multi_linear` (6.1034-6.1037).

```
safety_slip_multi_linear index x_first,0 y_first,0 x_first,1 y_first,1 ... [x_last,0 y_last,0 ...]
safety_slip_multi_linear_n index n          (number of lines between first and last)
safety_slip_multi_linear_segment_n index n  (points per line piece, default 10)
safety_slip_multi_linear_result index x_0 y_0 x_1 y_1 ... safety_factor   (output)
```

## Effect on the calculation

- At the end of every time step of the control index (switch `-yes`),
  for every slip surface, sampled in points k with weight A_k:
  Fs = sum A_k tau_mc_k / sum A_k |tau_k|, with
  tau_mc = c - sigma_n tan(phi) from `group_materi_plasti_mohr_coul`
  (or `_mohr_coul_direct`) of the element group, sigma_n = n.sigma.n and
  tau the shear stress along the surface.
- The stress in a point is the average over the elements containing the
  point of the element nodal stresses (`element_dof`) interpolated in the
  element.
- Points outside the mesh, in elements without a Mohr-Coulomb record, or
  with sigma_n > 0 (tension) do not count. A surface without counting
  points gets no result.
- Line pieces: `segment_n` equidistant points per piece (the pieces of a
  polyline do not share points), trapezoid weights.
- With `safety_slip_multi_linear_n` the record holds the first and the
  last polyline (same number of points); n polylines are interpolated
  equidistantly between them. Without it the record is one polyline.
- `safety_slip_multi_linear_result` gets the critical polyline of the
  index and its safety factor.
- The evaluated surfaces are written as `safety_slip_surface_areas`,
  `_elements`, `_nodes`, `_critical` (the overall critical surface),
  `_factor`, `_stresses` (sigma_n |tau| per point), `_stress_limit`
  (tau_mc per point), `_data_index`, `_data_name` (index = surface number).
- Not evaluated in control indices without the record; switch `-no` turns
  it off. The other surface types, sets, the refinement (switch = number)
  and the GiD plot follow in later F5 T6 units.

Same results as the Professional 02-08-2026 in elastic states (to the
solver accuracy).

## Example

```
control_safety_slip 10 -yes
safety_slip_multi_linear 0 0.5 2. 3.5 2.
```

Uniform sigxx = sigyy = -1, sigxy = 0.5, c 1, phi 30 deg:
`safety_slip_multi_linear_result 0 0.5 2 3.5 2 3.1547005384`.

Self-checking inputs in `validation-suite/test-2026`:
`safety_slip_multi_linear1` (analytic Fs, weights),
`safety_slip_multi_linear2` (points outside the mesh, polyline),
`safety_slip_multi_linear3` (tension points),
`safety_slip_multi_linear4` (_n, _segment_n, critical line),
`safety_slip_time1` (evaluated every step),
`safety_slip_error1` (first/last lines of different size).
