# group_beam_direction_z + group_beam_direction_z_reference_point + element_beam_direction_z (developer)

## Where

- `beam.cc`, `beam_direction_z_record()`: returns the record that applies
  to a beam element and its 3 values: `ELEMENT_BEAM_DIRECTION_Z`, else
  `GROUP_BEAM_DIRECTION_Z`, else `GROUP_BEAM_DIRECTION_Z_REFERENCE_POINT`,
  else 0 (priority measured on the Professional 02-08-2026, probe `prio`).
- `beam.cc`, `beam_3d()`: first statement; when a record applies it calls
  `beam_local()` and returns (2D and 3D). Otherwise the classic path is
  unchanged (`beam_2d()` in 2D, the `group_beam_plane` mapping in 3D).
- `beam.cc`, `beam_local()`: the 12-dof local-axes beam.
  - Geometry as `beam_2d()`: initial, old and new coordinates (3
    components, z = 0 in 2D); `-total_linear` / `-updated_without_rotation`
    keep the initial coordinates. The materi_velocity_integrated path uses
    `initial = coord - veli` (FOLLOW_MATERIAL), `old = initial + veli`
    (`beam_2d()` adds `coord` twice there; not copied).
  - Axes from the new coordinates: x = node 2 - node 1; z = record vector
    (direction) or `point - middle point` (reference point), minus its x
    component, normalized; y = z x x (`array_outproduct_3D`).
  - Checks: a direction record with `|cos(v, x0)| > 1e-3` on the initial
    geometry stops with `Error, ELEMENT_BEAM_DIRECTION_Z should be
    perpendicular to the beam length direction.` (the Professional text;
    `GROUP_...` for the group record); a reference point on the beam axis
    stops too.
  - Local stiffness (dofs per node ux uy uz rx ry rz): Euler-Bernoulli
    bending E*Izz on (uy, rz), E*Iyy on (uz, ry) with the sign of
    `ry = -dw/dx`, torsion G*J/L on rx; no axial term. L = initial length.
  - `Kg = T^t K T` (T = 4 blocks of the 3x3 axes matrix). State:
    `ELEMENT_BEAM_FORCE_MOMENT` (local, VERSION_NORMAL) +
    `K * T * du`, with du = vel*dtime for translations and new - old for
    rotations; global forces `T^t f`. Assembly as `beam_2d()` (matrix
    columns of the translations times dtime, lhside = diagonal).
  - 2D: only velx, vely and the single rotation (rotz) are assembled.
  - Writes `ELEMENT_BEAM_FORCE_MOMENT` (VERSION_NEW); does NOT write
    `ELEMENT_BEAM_DIRECTION` (2-value record of the planar beam) nor
    `ELEMENT_BEAM_MOMENT`.
- `database.cc`: `GROUP_BEAM_DIRECTION_Z`,
  `GROUP_BEAM_DIRECTION_Z_REFERENCE_POINT` (DOUBLE, length 3, class BEAM,
  required GROUP_TYPE) next to `GROUP_BEAM_PLANE`;
  `ELEMENT_BEAM_DIRECTION_Z` (DOUBLE, 3, class ELEMENT) next to
  `ELEMENT_BEAM_MOMENT`.
- `check.cc`: the three records need materi_displacement or
  materi_velocity_integrated, beam_rotation and ndim 2..3.
- `print_beam_force_moment.cc`: a beam element without
  `ELEMENT_BEAM_MOMENT` takes the 12 local components from
  `ELEMENT_BEAM_FORCE_MOMENT` (axial columns = truss force).
- Enums appended before `LAST_DUMMY` in `tochnog.h` / `tochnog-mod.h`;
  prototypes of `beam_direction_z_record()` and `beam_local()` in both.

## Hardcoded parameters

- `BEAM_LOCAL_PERPENDICULAR_TOL 1.e-3`: |cos| allowed between a direction
  record and the initial beam axis.
- `BEAM_LOCAL_ZERO_TOL 1.e-10`: relative length below which the local z
  axis of a reference point is zero.

## Verification

Professional 02-08-2026 on omen (`~/tochnog-runs/v4c`), 3D cantilevers:
`group_beam_direction_z1`, `_reference_point1`, `_priority1`,
`element_beam_direction_z1` (14 directions, prescribed tip vely with
rotz free: bending about both axes, torsion G*J with E 2, G 0.5),
`element_beam_force_moment2`; GNU = Professional within 1e-5. Default
axes of the Professional measured (probes `dflty`, `obl`): z = global z
projected perpendicular to x, global y for a beam along z.

## Pending refactorings / open items

- C37: the 3D beam without a record could use `beam_local()` with the
  Professional default axes (owner decision; changes existing results).
- BEAM-MEMORY-DEFAULT: GNU default `-updated`, Professional
  `-total_linear` (owner decision).
- BEAM-DIRECTION-LENGTH: `ELEMENT_BEAM_DIRECTION` is 2 values (Pro 9).
- BEAM-AXIAL-FREE: zero axial stiffness makes the direct solver fail.
- `print_el.cc` (`control_print_element -element_beam_moment`) still
  reads `ELEMENT_BEAM_MOMENT`, so it skips the local-axes beams.
- `group_beam_plane` and a local z record together: the record wins
  silently.

## Differences with the Professional version

See the user page: C37, BEAM-MEMORY-DEFAULT, BEAM-DIRECTION-LENGTH,
BEAM-AXIAL-FREE, and the 1e-3 perpendicularity tolerance.
