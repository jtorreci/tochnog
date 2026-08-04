# control_mesh_rotate_angle

## Files and functions

- `mesh.cc` — `void mesh_rotate_2d( double angle_deg )` (line 339). Rotates the 2D mesh around the z-axis.
- `top.cc` — called in `step_start` (lines 525–531) when `db_active_index( CONTROL_MESH_ROTATE_ANGLE, icontrol, VERSION_NORMAL )` is true; the angle is read and passed to `mesh_rotate_2d`.
- `tochnog.h` — enum `CONTROL_MESH_ROTATE_ANGLE` (line 260) and prototype (line 1366).
- `database.cc` — keyword registration (lines 771–774).
- `tochnog.h` — `PIRAD` constant (line 98), used for the degrees-to-radians conversion.

## Implementation details

- Converts the angle to radians: `a = angle_deg * PIRAD / 180.`, then `ca = cos(a)`, `sa = sin(a)`.
- Iterates over all active nodes `0..max_node` (`db_max_index`); for each node reads its coordinates and applies the planar rotation:
  - `new_coords[0] = ca*coords[0] - sa*coords[1]`
  - `new_coords[1] = sa*coords[0] + ca*coords[1]`
  - if `ndim == 3`, `new_coords[2] = coords[2]` (z unchanged).
- The same rotation is applied to `NODE` and, when present, to `NODE_START_REFINED`.
- Element connectivity is not touched.
- Closes with `mesh_has_changed( VERSION_NORMAL )`.
- `database.cc`: `type = DOUBLE_PRECISION`, `data_length = 1`, `data_class = CONTROL`.
- The `data_required` dependency to `CONTROL_MESH_ROTATE` was removed so the record works on its own.

## External dependencies

`cos`/`sin` from `<math.h>` and the `PIRAD` constant. Uses the internal database API (`db`, `db_active_index`, `db_max_index`).

## Hardcoded parameters / pending refactorings

- The rotation is hardcoded as a planar rotation around the z-axis on the first two coordinates; the 2D to 3D `control_mesh_rotate` remains unimplemented.
- When `ndim == 3` the rotation only affects `x`/`y` and leaves `z` untouched.
- `rot_length` and `idum_r` are read but not validated.
