# control_print_dof_line_n

## Implementación

- Stored as a CONTROL INTEGER record. Read in `print_dof_line_point()`
  of `print_dl.cc`:
  ```
  if ( db_active_index( CONTROL_PRINT_DOF_LINE_N, icontrol, ... ) ) {
    db( CONTROL_PRINT_DOF_LINE_N, icontrol, &n, ddum, ldum, ... );
    if ( n<1 ) db_error( CONTROL_PRINT_DOF_LINE_N, icontrol );
  }
  else
    n = DOF_LINE_DEFAULT_N;   // 5
  ```

## Diseño / decisiones

- Default `n = 5`: the legacy GNU `POST_LINE_N` default (post.cc: `else
  n = 5;`). The Professional manual does not state a default, so the
  existing GNU convention was reused (documented in the user manual).
- `n < 1` is an error (`db_error`): the record was given but is invalid;
  the legacy `POST_LINE` silently skipped `n <= 0`, but a print with zero
  points is meaningless.
- `n == 1` prints only the start point of the polyline (the first
  vertex).
- The points are equidistant over the TOTAL polyline length (see
  `control_print_dof_line_coordinates`).

## Detalles

- The line points are computed in `dof_line_point_position()`
  (`point ipoint` at fraction `ipoint/(n-1)` of the total length; the
  last point is exactly the last vertex).

## Pendiente

- None.
