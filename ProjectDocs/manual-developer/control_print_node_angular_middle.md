# control_print_node_angular_middle

## Implementación

- CONTROL DOUBLE_PRECISION record (2 or 3 values), read in
  `print_node()` (`print_node.cc`) into the `middle[MDIM]` array
  (default `{0,0,0}` when the record is not given).
  `data_required = CONTROL_PRINT_NODE_ANGULAR` in `database.cc`.
- Validation: 2D requires exactly 2 values, 3D exactly 3
  (manual 6.332: "In 2D you should not specify z_middle"); 1D is
  rejected by `control_print_node_angular` before this point.

## Diseño / decisiones

- Default middle point `(0,0,0)` when the record is absent (decision:
  the manual says the middle "should be specified" with the angular
  record but does not define the behaviour without it).
- Only the used number of coordinates is copied (`array_move(md,
  middle, md_len)`).

## Detalles

- Consumed exclusively by `node_angle_degrees()` through the `middle`
  array; not used by any other part of `print_node`.

## Pendiente

- None.
