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

## Validation (P4-E2b, completed)

- The record is implemented and VALIDATED: with `sigv=-100` the computed
  `kp = 0.1`; with `sigv=-200` the computed `kp = 0.05` — exactly the law
  `kp = a/(sigv/sig0)^b` with `a=0.1, b=1, sig0=100`.
- Regression test `validation-suite/test-2014/groundperm_vs.dat` (two-element
  column, base permeability 0.1, `a=0.1 b=1 sig0=100 min=0.01 max=1`,
  `sigv=-100`): middle-node pressure `0.498756` (expected `0.5` for uniform
  kp, within 0.25%). Passes.

### Key fix (P4-E2b)

The nodal stress is read via `db_dbl(NODE_DOF, inod, VERSION_NEW)` — NOT by
manually indexing `new_unknowns`. Two earlier attempts were wrong:

1. Indexing `new_unknowns[inol*nuknwn + dof]` with the ELEMENT-local `inol`
   reads the wrong global node (only node-local index 0 coincidentally
   matched).
2. The `h[]` passed to `groundflow_data` is a node-basis indicator
   (`[1,0,0,0]`, `[0,1,0,0]`, ...), NOT the integration-point shape
   functions, so interpolating with it is invalid.

The effective vertical stress is averaged over the element nodes (adequate
for a smoothly varying `sigv`). The vertical direction is the gravity axis
(largest `force_gravity` component).

## External dependencies

- `force_gravity_calculate()` (miscel.cc) for the vertical direction.
- `get_group_data()` with `GET_AND_CHECK` (needs `ldum=5` for the 5 params).
- Globals `stres_indx`, `pres_indx`, `npuknwn`, `nder`, `ndim`, `MDIM`.

## Hardcoded parameters / pending refactorings

- The tiny-stress threshold `1.e-10` is hardcoded (no user record).
- The `h` interpolation is performed inside `groundflow_data`; the caller
  must pass the shape functions and node count (both call sites updated).
- Validation of the flow result is pending (see above).
