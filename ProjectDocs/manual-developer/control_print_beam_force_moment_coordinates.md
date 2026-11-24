# control_print_beam_force_moment_coordinates

## Implementación

- Registered in `database.cc`:
  `CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES` (DOUBLE_PRECISION,
  `data_length = DATA_ITEM_SIZE`, `fixed_length = 0` so the same record
  holds 4 values in 2D and 6 in 3D; `data_required =
  CONTROL_PRINT_BEAM_FORCE_MOMENT`).
- Read inside `print_beam_force_moment()` (print_beam_force_moment.cc):
  ```
  db( CONTROL_PRINT_BEAM_FORCE_MOMENT_COORDINATES, icontrol, idum,
    coordinates, ncoord, VERSION_NORMAL, GET );
  if ( ncoord!=2*ndim ) db_error( ... );
  ```
- The cut segment is built in 3D arrays (z = 0 in 2D) and its length
  is checked: a zero-length cut is an error.

## Diseño / decisiones

- The coordinates record is MANDATORY: `db_error` when absent. Unlike
  `control_print_interface_stress` (which has a default cut), there is
  no sensible default cut for beams.
- The distance of the first output column is the projection of the
  closest point between the element axis and the cut onto the cut
  direction, i.e. `s * cut_length` with `s` the cut parameter returned
  by `segment_segment_distance()`.

## Pendiente

- None.
