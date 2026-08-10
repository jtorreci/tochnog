# group_groundflow_permeability_vertical_stress

## Files and functions

- `groundda.cc` — `groundflow_data()`:
  - Signature extended to `(element, gr, old_unknowns, new_unknowns,
    coord_ip, pe, C, h[], nnol)` to allow interpolation of the nodal stress.
  - After reading the base `GROUP_GROUNDFLOW_PERMEABILITY` and the capacity,
    if `GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS` is active, it reads the
    5 parameters (`a b sig0 minimum maximum`), determines the vertical
    direction from `force_gravity_calculate()` (largest component), and
    interpolates the effective vertical stress to the integration point:
    `sigv = sum h[inol] * ( sig_vertical[inol] - pres[inol] )`.
  - `kp = a / pow(|sigv|/sig0, b)`, clamped to `[minimum, maximum]`, applied
    to all components of `pe`.
- `groundfl.cc:61` and `general.cc:89` — calls updated to pass `h` and
  `nnol`.
- `tochnog.h` / `tochnog-mod.h` — enum `GROUP_GROUNDFLOW_PERMEABILITY_VERTICAL_STRESS`
  and the extended `groundflow_data` declaration (kept in sync).
- `database.cc` — `group_groundflow_permeability_vertical_stress`: DOUBLE,
  length 5, GROUNDFLOW class.
- `check.cc` — requires `groundflow_pressure` and `materi_stress`.

## Implementation details

- **Nodal dof layout** (tochnog): `nuknwn = npuknwn*nder`; the dof `k` of
  the `inol`-th element node is `new_unknowns[inol*nuknwn + k]`, where
  `k = stres_indx + stress_indx(idim,jdim)*nder`. `stress_indx(1,1)` is 3 in
  2D (not 4).
- **Vertical direction**: from the gravity vector. In tochnog 2D the gravity
  axis is `x` (index 0) by `force_gravity_calculate()`; the code picks the
  component with the largest magnitude, so it is geometry-independent.
- **Effective stress**: `sigv = sig_vertical - pres` at each node,
  interpolated with the shape functions `h[inol]`.
- **Division by zero guard**: `|sigv| <= 1e-10` keeps the base
  `group_groundflow_permeability` value.
- `stress_indx(idim,jdim)` is the symmetric tensor index (miscel.cc): for
  (1,1) in 2D → 3, for (0,0) → 0.

## Validation status (pending)

The record is implemented and the mechanism runs (the vertical effective
stress is read correctly, e.g. `sigv=-101` for `sig=-100, pres=1`). The
end-to-end flow test (`/tmp/gpv.dat`, a two-element column with prescribed
pressure 1 at the bottom and 0 at the top) is NOT yet a passing regression:

- Expected middle-node pressure `0.5` (uniform kp) is NOT reached
  (observed `0.714`).
- Instrumentation showed some integration points read `sigv=0`, which keeps
  the base permeability there — the interpolated nodal stress is not
  uniform across all integration points in the two-element mesh.

Remaining work (P4-E2b): resolve why some integration points interpolate
`sigv=0` (the `h` interpolation over the two-element quad4 mesh), then
calibrate the regression test. The constitutive formula and parameter
reading are verified working.

## External dependencies

- `force_gravity_calculate()` (miscel.cc) for the vertical direction.
- `get_group_data()` with `GET_AND_CHECK` (needs `ldum=5` for the 5 params).
- Globals `stres_indx`, `pres_indx`, `npuknwn`, `nder`, `ndim`, `MDIM`.

## Hardcoded parameters / pending refactorings

- The tiny-stress threshold `1.e-10` is hardcoded (no user record).
- The `h` interpolation is performed inside `groundflow_data`; the caller
  must pass the shape functions and node count (both call sites updated).
- Validation of the flow result is pending (see above).
