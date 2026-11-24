# control_print_dof_point_coordinates

## Implementación

- Stored as a CONTROL DOUBLE record (`data_length` DATA_ITEM_SIZE,
  `fixed_length 0`). Read in `print_dof_line_point()` of `print_dl.cc`
  with the `is_line=0` branch:
  ```
  db( CONTROL_PRINT_DOF_POINT_COORDINATES, icontrol, idum, coordinates,
    ncoord, VERSION_NORMAL, GET );
  if ( ncoord<ndim ) db_error( CONTROL_PRINT_DOF_POINT_COORDINATES, icontrol );
  ```
- The point is used directly as the single interpolation point
  (`npoints = 1`; the driver does not walk any polyline).

## Diseño / decisiones

- The record holds `ndim` values (x in 1D, x y in 2D, x y z in 3D). The
  parser reads until the next keyword, so a longer record is accepted but
  only the first `ndim` values are used; `ncoord < ndim` is an error.
- No move/method/group/eps records exist for the point family (manual
  6.281-6.283); the line defaults apply.

## Detalles

- The buffer is shared with the line machinery (`coordinates`); for the
  point it is never modified (no move).

## Pendiente

- None.
