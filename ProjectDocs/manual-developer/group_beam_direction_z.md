# group_beam_direction_z + group_beam_direction_z_reference_point + element_beam_direction_z (developer)

## Where

- `beam.cc`, `beam_direction_z_record()`: returns the record that applies
  to a beam element and its 3 values: `ELEMENT_BEAM_DIRECTION_Z`, else
  `GROUP_BEAM_DIRECTION_Z`, else `GROUP_BEAM_DIRECTION_Z_REFERENCE_POINT`,
  else 0 (priority measured on the Professional 02-08-2026, probe `prio`).
- `beam.cc`, `beam_3d()`: calls `beam_local()` when a record applies
  (2D and 3D) and for every 3D beam (record 0, F5 T2 V4d, owner decision
  C37); otherwise (2D without a record) `beam_2d()`. The former 3D path
  that mapped the beam onto `beam_2d()` in the `group_beam_plane` plane
  was removed.
- `beam.cc`, `beam_local()`: the 12-dof local-axes beam.
  - Geometry as `beam_2d()`: initial, old and new coordinates (3
    components, z = 0 in 2D); `-total_linear` / `-updated_without_rotation`
    keep the initial coordinates. The materi_velocity_integrated path uses
    `initial = coord - veli` (FOLLOW_MATERIAL), `old = initial + veli`
    (`beam_2d()` adds `coord` twice there; not copied).
  - Axes from the new coordinates: x = node 2 - node 1; z = record vector
    (direction) or `point - middle point` (reference point), minus its x
    component, normalized; y = z x x (`array_outproduct_3D`).
  - Record 0 (3D, no record): v = the `group_beam_plane` normal a x b
    when the group has that GNU-only record (a beam along the normal,
    `|cos| > 1 - 1e-3`, stops with `Error, the beam is perpendicular to
    its GROUP_BEAM_PLANE.`), else global z, or global y when the
    horizontal part of the initial beam direction is <= 1e-6
    (`BEAM_LOCAL_VERTICAL_TOL`, hardcoded; Professional default measured
    on probe `obl`, threshold not known).
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
  - Writes `ELEMENT_BEAM_FORCE_MOMENT` and `ELEMENT_BEAM_DIRECTION` (the
    3x3 axes matrix, rows x, y, z; F5 T2 V4d) (VERSION_NEW); does NOT
    write `ELEMENT_BEAM_MOMENT`.
  - `group_beam_memory` default `-total_linear` (F5 T2 V4d), mapped onto
    `-updated_without_rotation`; see [group_beam_memory](group_beam_memory.md).
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

- BEAM-AXIAL-FREE: zero axial stiffness makes the direct solver fail.
- `print_el.cc` (`control_print_element -element_beam_moment`) still
  reads `ELEMENT_BEAM_MOMENT`, so it skips the local-axes beams, which
  since F5 T2 V4d are every 3D beam.
- `group_beam_plane` and a local z record together: the record wins
  silently.

## Differences with the Professional version

See the user page: `group_beam_plane` (GNU-only), BEAM-AXIAL-FREE, and
the 1e-3 perpendicularity tolerance. C37, BEAM-MEMORY-DEFAULT and
BEAM-DIRECTION-LENGTH were resolved in F5 T2 V4d (2026-11-02).
