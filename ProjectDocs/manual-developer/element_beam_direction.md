# element_beam_direction (developer)

## Implementation

- `database.cc`: `ELEMENT_BEAM_DIRECTION` DOUBLE, `data_length` 3*3
  (2 until F5 T2 V4d), `version_all`, `print_only`, class ELEMENT.
  Allocated VERSION_NEW in `top.cc` (`db_allocate(..., max_element,
  VERSION_NEW, MINIMAL)`).
- `beam.cc`, `beam_local()`: PUT the axes matrix `R` (rows local x, y,
  z) at the end of the routine.
- `beam.cc`, `beam_2d()`: PUT `c s 0 -s c 0 0 0 1` with (c, s) =
  `new_beam_direction` (output only: `beam_2d()` recomputes its
  direction from the coordinates every step).
- Layout = Professional 02-08-2026 (probes `dfltx`, `obl`).

## Tests

- `test-2026/beam_3d_default_axes1` (14 directions, = Professional),
  `group_beam_plane_3d1`, `element_beam_force_moment_2d1` (2D).

## Pending refactorings

- None.
