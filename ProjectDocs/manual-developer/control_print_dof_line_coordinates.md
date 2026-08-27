# control_print_dof_line_coordinates

## Implementación

- Stored as one CONTROL record (DOUBLE_PRECISION, `data_length`
  DATA_ITEM_SIZE, `fixed_length 0`) in `database.cc`. Read in
  `print_dof_line_point()` of `print_dl.cc`:
  ```
  db( CONTROL_PRINT_DOF_LINE_COORDINATES, icontrol, idum, coordinates,
    ncoord, VERSION_NORMAL, GET );
  if ( ncoord%ndim!=0 || ncoord<2*ndim )
    db_error( CONTROL_PRINT_DOF_LINE_COORDINATES, icontrol );
  ```
- `dof_line_point_position()` walks the polyline: total length =
  sum of segment lengths; point `ipoint` sits at fraction
  `ipoint/(n-1)` of the total length (`n == 1` -> start point; zero
  total length -> first vertex).

## Diseño / decisiones

- Variable-length record (the manual gives `x_0 y_0 z_0 x_1 y_1 z_1 ...`
  without a fixed count). The parser stops at the next keyword, so the
  record length is exactly `nvertices*ndim`.
- Equidistant over the TOTAL length was chosen over per-segment
  distribution (the manual does not specify; equal spacing over the whole
  polyline is the natural reading of "how many points will be printed
  along the line").
- Validation: at least 2 vertices (`ncoord >= 2*ndim`) and
  `ncoord % ndim == 0`.

## Detalles

- The same `coordinates` buffer is updated in place by `_move` and PUT
  back, so the record holds the CURRENT (possibly moved) vertices.

## Pendiente

- None.
