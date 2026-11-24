# bounda_dof_cylindrical

## Files and functions

- `bounda.cc` — within `bounda()`.
  - Local `bounda_dof_cylindrical[6]` declared at line 49.
  - Array zeroed and read per `iboun`: `db( BOUNDA_DOF_CYLINDRICAL, iboun,
    idum, bounda_dof_cylindrical, ldum, VERSION_NORMAL, GET_IF_EXISTS )`
    (lines 207-209).
  - Projection applied in the velocity branch of the dof loop (lines 572-612),
    right after the prescribed dof value is computed, before derivatives.
- `database.cc` — keyword registration (lines 200-204): name
  `bounda_dof_cylindrical`, `type = DOUBLE_PRECISION`, `data_length = 6`,
  `data_class = BOUNDA`, `data_required = BOUNDA_DOF`.
- `tochnog.h` — enum `BOUNDA_DOF_CYLINDRICAL` (line 144).
- `tochnog-mod.h` — mirror enum `BOUNDA_DOF_CYLINDRICAL` (line 137), must stay in sync.

## Implementation details

- The two line points are read with `GET_IF_EXISTS`; the projection is inactive
  when no `bounda_dof_cylindrical` record exists for `iboun`.
- Active only when at least one of the first three components (the first point)
  is non-zero AND the dof is a velocity dof
  (`iuknwn>=vel_indx && iuknwn<vel_indx+ndim*nder`, line 577). It applies to
  the dofs of `materi_velocity`.
- For each velocity dof, the in-plane dimension is derived as
  `idim_c = (iuknwn - vel_indx) / nder` (line 578), and the node coordinates
  are fetched via `db( NODE, inod, ... )` (lines 586-587).
- The line axis is built as `axis = p2 - p1` over `MDIM` components
  (lines 581-585); the first `ndim` components of the axis are used below.
- The node is projected onto the line with the parameter
  `t = (rc . axis) / (axis . axis)` where `rc = coord - p1` (lines 588-598).
  A `aa>0.` guard skips the projection when the two points coincide.
- The projection point is `proj = p1 + t*axis`; the radial distance is
  `r = sqrt(sum_k (coord_k - proj_k)^2)` over the first `ndim` coordinates
  (lines 599-605).
- The prescribed value is scaled by `(coord[idim_c] - proj[idim_c]) / r`
  (lines 607-609). This decomposes the prescribed radial magnitude onto the
  `idim_c` axis. A `r>0.` guard avoids division by zero when the node lies on
  the axis.
- Because the scale factor is applied to each velocity dof individually, the
  velocity vector of a node is parallel to `(node - proj)`, i.e. purely radial
  to the line.
- In 2D only the `x, y` dofs exist, so only the in-plane projection is used and
  the line behaves as perpendicular to the plane; in 3D all three components of
  both points participate.
- The factor is applied after the `bounda_factor`/`bounda_factor_parabolic_x`/
  `bounda_water` load computation and before the derivative update, so it
  multiplies the final prescribed velocity.

## External dependencies

None. Internal database API only; reuses `materi_velocity`, the velocity dof
offset `vel_indx`, the number of unknowns per dof `nder` and the global `ndim`.

## Hardcoded parameters / pending refactorings

- The line is hardcoded to 6 components (two points xyz); only the first `ndim`
  are used in the projection, and the activity check only inspects the first
  point.
- The projection is only implemented for `bounda_dof_cylindrical`; the
  analogous `bounda_dof_radial` (point axis) shares the same structure and
  could be refactored into a common "project onto an axis" helper.
- In 2D the user must supply consistent points; nothing validates that the two
  points differ (`aa>0.` silently disables the projection).
