# control_print_node_angular

## Implementación

- CONTROL INTEGER record (up to 3 values: switch_x switch_y switch_z),
  read in `print_node()` (`print_node.cc`). `data_required =
  CONTROL_PRINT_NODE` in `database.cc`.
- Validation (manual 6.331):
  - 1D -> `db_error` ("In 1D you cannot use this record").
  - 2D: exactly 2 values and both `-yes` ("you should not specify
    switch_z and you can only use -yes -yes").
  - 3D: exactly 3 values, one of the three combinations
    `-yes -yes -no` / `-no -yes -yes` / `-yes -no -yes`; anything else
    -> `db_error`.

## Diseño / decisiones

- **Angle in DEGREES** (decision; the manual says "the number of
  degrees"): `180/pi * atan2(dy, dx)` where the axis pair is selected
  by the switches (helper `node_angle_degrees()` in `print_node.cc`):
  - axis 0 (`-yes -yes -no`): dx = x−xm, dy = y−ym (from +x to +y);
  - axis 1 (`-no -yes -yes`): dx = y−ym, dy = z−zm (from +y to +z);
  - axis 2 (`-yes -no -yes`): dx = x−xm, dy = z−zm (from +x to +z).
- The middle point is `control_print_node_angular_middle` (2 values in
  2D, 3 in 3D; missing -> `(0,0,0)`, decision: the manual says it
  "should be specified" but does not define an error).
- The angle REPLACES the coordinates in the line: `angle <value>`.
- With angular, `control_print_node_sort` accepts ONLY `-angle`
  (manual: "In case you use -angular ... you can set the sort_method to
  -angle"): any other sort method -> `db_error`.

## Detalles

- The sort key of the angular files is the same `node_angle_degrees`
  value written to the line, so `-angle` sorts by the printed column.

## Pendiente

- None.
