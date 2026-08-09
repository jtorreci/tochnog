# bounda_dof_radial

## Files and functions

- `bounda.cc` — within `bounda()`.
  - Local `bounda_dof_radial[3]` declared at line 48.
  - Array zeroed and read per `iboun`: `db( BOUNDA_DOF_RADIAL, iboun, idum,
    bounda_dof_radial, ldum, VERSION_NORMAL, GET_IF_EXISTS )` (lines 204-206).
  - Projection applied in the velocity branch of the dof loop (lines 553-571),
    right after the prescribed dof value is computed, before derivatives.
- `database.cc` — keyword registration (lines 206-210): name
  `bounda_dof_radial`, `type = DOUBLE_PRECISION`, `data_length = 3`,
  `data_class = BOUNDA`, `data_required = BOUNDA_DOF`.
- `tochnog.h` — enum `BOUNDA_DOF_RADIAL` (line 145).
- `tochnog-mod.h` — mirror enum `BOUNDA_DOF_RADIAL` (line 138), must stay in sync.

## Implementation details

- The center point is read with `GET_IF_EXISTS`; the projection is inactive
  when no `bounda_dof_radial` record exists for `iboun`.
- Active only when at least one component is non-zero AND the dof is a velocity
  dof (`iuknwn>=vel_indx && iuknwn<vel_indx+ndim*nder`, line 556). It applies to
  the dofs of `materi_velocity`.
- For each velocity dof, the in-plane dimension is derived as
  `idim_r = (iuknwn - vel_indx) / nder` (line 557), and the node coordinates
  are fetched via `db( NODE, inod, ... )` (lines 559-560).
- The distance is computed as `r = sqrt(sum_k (coord_k - point_k)^2)` over the
  first `ndim` coordinates (lines 561-565).
- The prescribed value is scaled by `(coord[idim_r] - point[idim_r]) / r`
  (lines 567-569). This decomposes the prescribed radial magnitude onto the
  `idim_r` axis. A `r>0.` guard avoids division by zero when the node coincides
  with the center point.
- Because the scale factor is applied to each velocity dof individually, the
  three dofs of a node form a vector parallel to `(node - point)`, i.e. purely
  radial.
- The factor is applied after the `bounda_factor`/`bounda_factor_parabolic_x`/
  `bounda_water` load computation and before the derivative update, so it
  multiplies the final prescribed velocity.

## External dependencies

None. Internal database API only; reuses `materi_velocity`, the velocity dof
offset `vel_indx`, the number of unknowns per dof `nder` and the global `ndim`.

## Hardcoded parameters / pending refactorings

- The center point is hardcoded to 3 components; only the first `ndim` are used
  in the distance computation.
- The projection is only implemented for `bounda_dof_radial`.
- PENDING: `bounda_dof_cylindrical` is registered in `database.cc` (lines
  200-204, `data_length = 6`) and its array is read in `bounda.cc` (lines
  207-209), but no logic uses it yet — the keyword is non-functional.
